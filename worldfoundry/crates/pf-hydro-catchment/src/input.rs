//! Reader adapted from retained Q20 pf-hydro-mr (see lane BASE.json).
//! Direct, read-only C7 snapshot and native trace ingestion. No pf-surface dependency.
use crate::{Grid, Result, NONE};
use serde::{Deserialize, Deserializer};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File},
    io::{BufReader, Read},
    path::{Path, PathBuf},
};

fn bits<'de, D: Deserializer<'de>>(d: D) -> std::result::Result<f64, D::Error> {
    let s = String::deserialize(d)?;
    let b = u64::from_str_radix(&s, 16).map_err(serde::de::Error::custom)?;
    let x = f64::from_bits(b);
    if s.len() != 16 || !x.is_finite() {
        return Err(serde::de::Error::custom("nonfinite/invalid binary64"));
    }
    Ok(x)
}
fn bitvec<'de, D: Deserializer<'de>>(d: D) -> std::result::Result<Vec<f64>, D::Error> {
    #[derive(Deserialize)]
    struct V(#[serde(deserialize_with = "bits")] f64);
    Ok(Vec::<V>::deserialize(d)?.into_iter().map(|v| v.0).collect())
}
#[derive(Deserialize)]
pub struct Clock {
    pub step: u64,
    #[serde(deserialize_with = "bits")]
    pub start_my: f64,
    #[serde(deserialize_with = "bits")]
    pub end_my: f64,
}
#[derive(Deserialize, Clone, Debug)]
pub struct Row {
    pub cell: usize,
    pub receiver: Option<usize>,
    #[serde(deserialize_with = "bits")]
    pub incoming_m3: f64,
    #[serde(deserialize_with = "bits")]
    pub detached_m3: f64,
    #[serde(deserialize_with = "bits")]
    pub retained_m3: f64,
    #[serde(deserialize_with = "bits")]
    pub outgoing_m3: f64,
}
#[derive(Deserialize)]
pub struct Budget {
    #[serde(deserialize_with = "bitvec")]
    pub final_dz_m: Vec<f64>,
}
#[derive(Deserialize)]
pub struct Frame {
    pub clock: Clock,
    pub wet: Vec<bool>,
    pub ocean: Vec<bool>,
    pub receivers: Vec<Option<usize>>,
    pub rows: Vec<Row>,
    pub surface_budget: Budget,
}
#[derive(Deserialize)]
struct Trace {
    schema: String,
    grid: [usize; 2],
    frames: Vec<Frame>,
    sources: Value,
}
#[derive(Deserialize)]
struct Archive {
    trace: Trace,
    state_digest: String,
    recipe_digest: String,
}
#[derive(Deserialize)]
struct Meta {
    nlat: usize,
    nlon: usize,
    seed: u64,
    t_my: f64,
    steps: u64,
}

pub struct Input {
    pub grid: Grid,
    pub seed: u64,
    pub elevation: Vec<f64>,
    pub initial_z: Vec<f64>,
    pub relief: Vec<f64>,
    pub sediment: Vec<f64>,
    pub thickness: Vec<f64>,
    pub frame: Frame,
    pub runoff_m_yr: f64,
    pub published_flow: Vec<usize>,
    pub published_q: Vec<f64>,
    pub published_area: Vec<f64>,
    pub sea: f64,
    pub provenance: Value,
}
pub fn floats(path: &Path, n: usize) -> Result<Vec<f64>> {
    let b = fs::read(path)?;
    if b.len() != 4 * n {
        return Err(format!("{}: expected {n} f32 cells", path.display()).into());
    }
    let v: Vec<_> = b
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()) as f64)
        .collect();
    if v.iter().any(|v| !v.is_finite()) {
        return Err(format!("{}: nonfinite field", path.display()).into());
    }
    Ok(v)
}
pub fn uints(path: &Path, n: usize) -> Result<Vec<usize>> {
    let b = fs::read(path)?;
    if b.len() != 4 * n {
        return Err("u32 shape mismatch".into());
    }
    Ok(b.chunks_exact(4)
        .map(|b| {
            let v = u32::from_le_bytes(b.try_into().unwrap());
            if v == u32::MAX {
                NONE
            } else {
                v as usize
            }
        })
        .collect())
}
pub fn digest(path: &Path) -> Result<String> {
    let mut f = File::open(path)?;
    let mut h = Sha256::new();
    let mut b = [0; 65536];
    loop {
        let n = f.read(&mut b)?;
        if n == 0 {
            break;
        }
        h.update(&b[..n]);
    }
    Ok(format!("{:x}", h.finalize()))
}
impl Input {
    pub fn read(root: &Path) -> Result<Self> {
        let mut snapshots = Vec::new();
        for entry in fs::read_dir(root.join("snaps"))? {
            let p = entry?.path().join("meta.json");
            if p.is_file() {
                let m: Meta = serde_json::from_reader(File::open(&p)?)?;
                snapshots.push((m, p));
            }
        }
        snapshots.sort_by(|a, b| a.0.t_my.total_cmp(&b.0.t_my));
        let (meta, meta_path) = snapshots.pop().ok_or("no snapshot metadata")?;
        let snap = meta_path.parent().unwrap();
        let pin: Value = serde_json::from_reader(File::open(root.join("terrain-inputs.json"))?)?;
        let recipe: Value = serde_json::from_reader(File::open(root.join("recipe.json"))?)?;
        let terrain: Value =
            serde_json::from_reader(File::open(root.join("terrain/artifact.json"))?)?;
        let radius = pin["radius_m"].as_f64().ok_or("missing planet radius")?;
        let grid = Grid::new(meta.nlat, meta.nlon, radius)?;
        let n = grid.len();
        if pin["model_time_my"].as_f64() != Some(meta.t_my)
            || pin["seed_path"]["seed"].as_u64() != Some(meta.seed)
        {
            return Err("snapshot/terrain pin epoch or seed mismatch".into());
        }
        let archive: Archive =
            serde_json::from_reader(BufReader::new(File::open(root.join("fluvial-trace.json"))?))?;
        if archive.trace.schema != "vz.native-fluvial-trace/v2"
            || archive.trace.grid != [meta.nlat, meta.nlon]
            || archive.state_digest != pin["state_digest"].as_str().unwrap_or("")
            || archive.recipe_digest != pin["recipe_digest"].as_str().unwrap_or("")
        {
            return Err("native trace schema/grid/world pin mismatch".into());
        }
        let trace_sources = archive.trace.sources;
        let frame = archive
            .trace
            .frames
            .into_iter()
            .find(|f| f.clock.step == meta.steps)
            .ok_or("no native trace at latest snapshot step")?;
        if frame.clock.end_my != meta.t_my
            || frame.wet.len() != n
            || frame.ocean.len() != n
            || frame.receivers.len() != n
            || frame.surface_budget.final_dz_m.len() != n
        {
            return Err("trace shape/clock mismatch".into());
        }
        let elevation = floats(&snap.join("elevation_m.f32"), n)?;
        let sediment = floats(&snap.join("sediment_m.f32"), n)?;
        let thickness = floats(&snap.join("thickness_m.f32"), n)?;
        let relief = floats(&root.join("subgrid_relief_m.f32"), n)?;
        if sediment
            .iter()
            .chain(&thickness)
            .chain(&relief)
            .any(|x| *x < 0.)
        {
            return Err("negative sediment/crust/relief".into());
        }
        let initial_z = elevation
            .iter()
            .zip(&frame.surface_budget.final_dz_m)
            .map(|(z, dz)| z - dz)
            .collect();
        let surface = &recipe["inputs"]["profile"]["surface"];
        if surface["rainfall_mode"] != "Uniform" {
            return Err("unsupported C7 runoff model: need native per-cell runoff archive".into());
        }
        let runoff_m_yr = surface["uniform_runoff_m_per_yr"]
            .as_f64()
            .ok_or("missing runoff")?;
        if !runoff_m_yr.is_finite() || runoff_m_yr < 0. {
            return Err("invalid runoff".into());
        }
        let sea = terrain["header"]["domain"]["sea_level_m"]
            .as_f64()
            .ok_or("missing sea datum")?;
        let published_flow = uints(&root.join("hydrology/fields/flow_to_element.u32"), n)?;
        let published_q = floats(&root.join("hydrology/fields/discharge_proxy_m3_s.f32"), n)?;
        let published_area = floats(&root.join("hydrology/fields/contributing_area_m2.f32"), n)?;
        let paths: Vec<PathBuf> = vec![
            meta_path.clone(),
            snap.join("elevation_m.f32"),
            snap.join("sediment_m.f32"),
            snap.join("thickness_m.f32"),
            root.join("subgrid_relief_m.f32"),
            root.join("fluvial-trace.json"),
            root.join("terrain-inputs.json"),
            root.join("recipe.json"),
            root.join("terrain/artifact.json"),
            root.join("hydrology/fields/flow_to_element.u32"),
            root.join("hydrology/fields/discharge_proxy_m3_s.f32"),
            root.join("hydrology/fields/contributing_area_m2.f32"),
        ];
        let mut hashes = Vec::new();
        for p in paths {
            hashes.push(json!({"path":p,"sha256":digest(&p)?}));
        }
        let changes = (0..n)
            .filter(|&i| frame.receivers[i].unwrap_or(NONE) != published_flow[i])
            .count();
        let provenance = json!({"schema":"pf.hydro-mr-input/v1","root":root,"snapshot":snap,"seed":meta.seed,"shape":[meta.nlat,meta.nlon],"radius_m":radius,"epoch_start_my":frame.clock.start_my,"epoch_end_my":frame.clock.end_my,"state_digest":archive.state_digest,"recipe_digest":archive.recipe_digest,"native_source_digests":trace_sources,"published_vs_native_receiver_changes":changes,"files":hashes,"pre_surface_reconstruction":"snapshot minus native final surface dz; conditional, not a saved start-of-pass terrain"});
        Ok(Self {
            grid,
            seed: meta.seed,
            elevation,
            initial_z,
            relief,
            sediment,
            thickness,
            frame,
            runoff_m_yr,
            published_flow,
            published_q,
            published_area,
            sea,
            provenance,
        })
    }
}
