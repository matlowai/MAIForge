//! native pipeline: one physical catchment solve, independent of request charts.
//! native pipeline supplies the priority-flood/metric-D8 kernel unchanged. Native coastal
//! outlet cells partition work, not material authority or permitted geography.
//! Internal T0 edges are ordinary mesh edges. Routing heads never become bed.
use crate::{planet_relief::Inputs, R};
use pf_hydro_catchment::{
    domain::Hydro,
    input::{EndpointInput, Input},
    terrain::{distance, terrain},
    Mesh, Terminal, NONE, YEAR_S,
};
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
    sync::{Arc, Mutex},
    time::Instant,
};

pub const MODEL: &str = "planet-catchment-flow-hillslope/6";
pub const LIMITS: &str = "Conditional initial field: monotone C1 reconstruction of native endpoint centre observations plus unchanged native pipeline spherical geological roughness. Native cell means are not imposed on final fine terrain. Then 2 My native implicit n=1 stream-power with conservative multiple-downhill flow on native pipeline heads, and conservative finite-volume hillslope diffusion on a catchment-joint 128-subdivision lattice (about 543 m north-south for the 576 by 288 reference grid); 32m/2m outputs sample that field. Two-native-cell halo with explicit zero external inflow and fixed endpoint boundary heads. Native water-body beds/stages and coarse shorelines retained. Diagnostic fills are not lakes or eroded bed; incision in dry filled depressions is unresolved and reported. Hillslope bulk flux and boundary exchange are conserved; fluvial sediment destinations, evolved meandering and new tectonic history are unresolved. Finite material uses the separately versioned radial-volume lift contract.";
pub const RESOLUTION: usize = 128;
const STEPS: usize = 4;
const DT_MY: f64 = 0.5;

/// Only the actual water fields required by refinement. A historical routing
/// frame and its pre-pass reconstruction are deliberately not represented here.
pub struct EndpointHydrology {
    pub body: Vec<usize>,
    pub fill: Vec<f64>,
    pub records: BTreeMap<usize, pf_hydro_catchment::domain::Body>,
    pub provenance: serde_json::Value,
}
impl EndpointHydrology {
    fn from_legacy(h: Hydro) -> Self {
        Self {
            body: h.body,
            fill: h.fill,
            records: h.records,
            provenance: h.provenance,
        }
    }
}

/// Load a native endpoint without requiring snapshot or fluvial-trace archives.
/// Every artifact is verified; the auxiliary subgrid field is separately hashed
/// and its original manifest pin is checked when one is present.
pub fn load_published(root: &Path) -> R<(Inputs, EndpointInput, EndpointHydrology)> {
    use pf_contracts::artifact::{
        read_artifact_header, read_descriptor_bytes, verify_artifact, Dtype,
    };
    use pf_contracts::domain::{HydrologyDomainV0, WaterBodyKindV0};
    use pf_contracts::terrain::StratigraphyDomainV0;
    use serde_json::{json, Value};
    let direct = root.join("terrain/artifact.json").is_file();
    let wrapped = root.join("run/terrain/artifact.json").is_file();
    let published = match (direct, wrapped) {
        (true, false) => root.to_path_buf(),
        (false, true) => root.join("run"),
        _ => return Err("endpoint: expected one native root or run/ wrapper".into()),
    };
    let td = published.join("terrain");
    let hd = published.join("hydrology");
    let t = read_artifact_header::<Value>(&td).map_err(|e| e.to_string())?;
    let h = read_artifact_header::<Value>(&hd).map_err(|e| e.to_string())?;
    if verify_artifact::<Value>(&td).map_err(|e| e.to_string())? != t.digest
        || verify_artifact::<Value>(&hd).map_err(|e| e.to_string())? != h.digest
    {
        return Err("endpoint: source changed during verification".into());
    }
    if t.header.schema != "pf.terrain/v0"
        || h.header.schema != "pf.hydrology/v0"
        || !h.header.semantic.inputs.contains(&t.digest)
        || t.header.semantic.seed_path != h.header.semantic.seed_path
    {
        return Err("endpoint: native source identity mismatch".into());
    }
    let tf = t
        .header
        .semantic
        .fields
        .iter()
        .find(|f| f.name == "elevation_m" && f.dtype == Dtype::F32Le && f.units == "m")
        .ok_or("endpoint elevation contract")?;
    if tf.shape.len() != 2 || tf.order != "C" {
        return Err("endpoint grid layout".into());
    }
    let rows = usize::try_from(tf.shape[0]).map_err(|_| "endpoint rows")?;
    let cols = usize::try_from(tf.shape[1]).map_err(|_| "endpoint columns")?;
    let domain = t
        .header
        .semantic
        .domains
        .iter()
        .find(|d| d.domain_id == tf.domain)
        .ok_or("endpoint domain")?;
    if domain.discretization.backend != "latlon-v0" || domain.datum.center != [0.0; 3] {
        return Err("endpoint requires planet-centred latitude/longitude grid".into());
    }
    let grid = pf_hydro_catchment::Grid::new(rows, cols, domain.datum.radius_m)
        .map_err(|e| e.to_string())?;
    let n = grid.len();
    let bytes = |dir: &Path,
                 a: &pf_contracts::artifact::ArtifactOnDisk<Value>,
                 name: &str,
                 dtype: Dtype|
     -> R<Vec<u8>> {
        let f = a
            .header
            .semantic
            .fields
            .iter()
            .find(|f| f.name == name && f.dtype == dtype)
            .ok_or_else(|| format!("endpoint missing {name}"))?;
        let d = a
            .header
            .semantic
            .domains
            .iter()
            .find(|d| d.domain_id == f.domain)
            .ok_or("endpoint field domain")?;
        let expected_unit = match name {
            "elevation_m" | "fill_level_m" => "m",
            "flow_to_element" | "water_body_id" => "id",
            "contributing_area_m2" => "m2",
            "discharge_proxy_m3_s" => "m3/s",
            _ => return Err("unsupported endpoint field".into()),
        };
        if f.units != expected_unit {
            return Err(format!("endpoint field unit mismatch:{name}"));
        }
        if f.shape != tf.shape
            || f.order != "C"
            || d.datum != domain.datum
            || d.discretization.topology_digest != domain.discretization.topology_digest
        {
            return Err("endpoint field grid mismatch".into());
        }
        read_descriptor_bytes(dir, &a.header.semantic, f).map_err(|e| e.to_string())
    };
    let f32s = |b: Vec<u8>| -> R<Vec<f32>> {
        let v: Vec<f32> = b
            .chunks_exact(4)
            .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
            .collect();
        if v.len() != n || v.iter().any(|x| !x.is_finite()) {
            return Err("endpoint nonfinite/incorrect field".into());
        }
        Ok(v)
    };
    let indices = |b: Vec<u8>| -> R<Vec<usize>> {
        let v: Vec<usize> = b
            .chunks_exact(4)
            .map(|b| {
                let i = u32::from_le_bytes(b.try_into().unwrap());
                if i == u32::MAX {
                    NONE
                } else {
                    i as usize
                }
            })
            .collect();
        if v.len() != n {
            return Err("endpoint index count".into());
        }
        Ok(v)
    };
    let elevation = f32s(bytes(&td, &t, "elevation_m", Dtype::F32Le)?)?;
    let flow = indices(bytes(&hd, &h, "flow_to_element", Dtype::U32Le)?)?;
    let body = indices(bytes(&hd, &h, "water_body_id", Dtype::U32Le)?)?;
    let fill = f32s(bytes(&hd, &h, "fill_level_m", Dtype::F32Le)?)?;
    let q = f32s(bytes(&hd, &h, "discharge_proxy_m3_s", Dtype::F32Le)?)?;
    let area = f32s(bytes(&hd, &h, "contributing_area_m2", Dtype::F32Le)?)?;
    let epoch = t.header.domain["t_my"]
        .as_f64()
        .filter(|v| v.is_finite())
        .ok_or("endpoint epoch")?;
    let sea = t.header.domain["sea_level_m"]
        .as_f64()
        .filter(|v| v.is_finite())
        .ok_or("endpoint sea level")?;
    let native_hydro: HydrologyDomainV0 =
        serde_json::from_value(h.header.domain.clone()).map_err(|e| e.to_string())?;
    native_hydro.validate().map_err(|e| e.to_string())?;
    if native_hydro.t_my != epoch {
        return Err("endpoint epoch mismatch".into());
    }
    let records: BTreeMap<usize, pf_hydro_catchment::domain::Body> =
        serde_json::from_value::<Vec<pf_hydro_catchment::domain::Body>>(
            h.header.domain["water_bodies"].clone(),
        )
        .map_err(|e| e.to_string())?
        .into_iter()
        .map(|b| (b.body_id, b))
        .collect();
    if body.iter().any(|&i| i != NONE && !records.contains_key(&i))
        || records.values().any(|b| b.seed_element >= n)
    {
        return Err("endpoint water identity".into());
    }
    for (i, &b) in body.iter().enumerate() {
        if b != NONE {
            let stage = records[&b].water_surface_m;
            let tol =
                2.0 * f32::EPSILON as f64 * (elevation[i] as f64).abs().max(stage.abs()).max(1.0);
            if elevation[i] as f64 > stage + tol {
                return Err("endpoint water bed above stage".into());
            }
        }
    }
    if native_hydro
        .water_bodies
        .iter()
        .any(|b| b.kind == WaterBodyKindV0::Marine && (b.water_surface_m - sea).abs() > 1e-6)
    {
        return Err("endpoint marine datum".into());
    }
    let runoff = t.header.domain["assumptions"]
        .as_array()
        .and_then(|a| a.iter().find(|a| a["id"] == "A-RAIN-UNITS"))
        .and_then(|a| a["value"].as_f64())
        .ok_or("endpoint native runoff constant")?;
    if !runoff.is_finite() || runoff < 0.0 {
        return Err("endpoint runoff".into());
    }
    // This sampler's accepted construction uses uniform runoff. Refuse source
    // data whose published discharge/area contradict that assumption.
    for (&q, &a) in q.iter().zip(&area) {
        if q < 0.0 || a < 0.0 {
            return Err("endpoint negative discharge/area".into());
        }
        let expected = a as f64 * runoff / YEAR_S;
        if (q as f64 - expected).abs() > 1e-4 * expected.abs().max(q as f64) + 1e-12 {
            return Err("endpoint requires uniform runoff; published discharge/area differ".into());
        }
    }
    let sg = root.join("subgrid_relief_m.f32");
    let sg_bytes = std::fs::read(&sg).map_err(|e| format!("endpoint subgrid relief:{e}"))?;
    let relief = f32s(sg_bytes.clone())?;
    // Pin the consumed bytes, not a second read that could observe a different file.
    let sg_sha = pf_hydro_catchment::input::sha256_bytes(&sg_bytes);
    let mut auxiliary_binding="supplied native auxiliary field; bytes pinned to this endpoint, no historical trace claimed";
    let pin_path = root.join("terrain-inputs.json");
    if pin_path.is_file() {
        let pin: Value =
            serde_json::from_slice(&std::fs::read(pin_path).map_err(|e| e.to_string())?)
                .map_err(|e| e.to_string())?;
        if pin["subgrid_sha256"] != sg_sha
            || pin["artifact_digests"]["terrain"] != t.digest
            || pin["artifact_digests"]["hydrology"] != h.digest
            || pin["model_time_my"].as_f64() != Some(epoch)
            || pin["seed_path"]["seed"].as_u64() != Some(t.header.semantic.seed_path.seed)
        {
            return Err("endpoint auxiliary source manifest mismatch".into());
        }
        auxiliary_binding = "original native auxiliary manifest verified";
    }
    let sd = t
        .header
        .semantic
        .sidecars
        .iter()
        .find(|s| s.logical_name == "stratigraphy")
        .ok_or("endpoint requires native stratigraphy for erosion lithology")?;
    let sb = std::fs::read(td.join(&sd.path)).map_err(|e| e.to_string())?;
    if sb.len() as u64 != sd.byte_length || pf_contracts::sim_module::digest_bytes(&sb) != sd.digest
    {
        return Err("endpoint strata changed".into());
    }
    let strata: StratigraphyDomainV0 = serde_json::from_slice(&sb).map_err(|e| e.to_string())?;
    strata.validate().map_err(|e| e.to_string())?;
    if strata.nlat as usize != rows || strata.nlon as usize != cols || strata.t_my != epoch {
        return Err("endpoint strata grid/epoch".into());
    }
    let mut lithology = Vec::with_capacity(n);
    for c in 0..n {
        let a = strata.layer_offset[c] as usize;
        if a == strata.layer_offset[c + 1] as usize {
            return Err("endpoint empty native column".into());
        }
        lithology.push(strata.lith[a]);
    }
    let erosion = &t.header.domain["erosion"];
    let law = |key: &str| -> R<f64> {
        erosion[key]
            .as_f64()
            .filter(|x| x.is_finite() && *x > 0.)
            .ok_or_else(|| format!("endpoint invalid law:{key}"))
    };
    if law("stream_power_n")? != 1.0 {
        return Err("endpoint requires native n=1 stream power".into());
    }
    let head = t.header.domain["assumptions"]
        .as_array()
        .and_then(|a| a.iter().find(|a| a["id"] == "A-SUBGRID-CHANNEL-HEAD"))
        .and_then(|a| a["value"].as_f64())
        .filter(|v| v.is_finite() && *v > 0.)
        .ok_or("endpoint channel head length")?;
    let seed = t.header.semantic.seed_path.seed;
    let provenance = json!({"schema":"maiforge.native-endpoint/v1","seed":seed,"grid":[rows,cols],"radiusM":grid.radius,"epochMy":epoch,"terrainArtifact":t.digest,"hydrologyArtifact":h.digest,"stratigraphy":sd.digest,"subgridSha256":sg_sha,"auxiliaryBinding":auxiliary_binding,"history":"No pre-pass elevation or native trace was fabricated or required."});
    let input = EndpointInput {
        grid,
        seed,
        elevation: elevation.iter().map(|&x| x as f64).collect(),
        relief: relief.iter().map(|&x| x as f64).collect(),
        runoff_m_yr: runoff,
        published_flow: flow.clone(),
        sea,
        provenance,
    };
    input.validate().map_err(|e| e.to_string())?;
    let source = Inputs {
        seed,
        rows,
        cols,
        radius_m: grid.radius,
        elevation_m: elevation.clone(),
        relief_m: relief,
        receiver: flow
            .iter()
            .map(|&x| if x == NONE { u32::MAX } else { x as u32 })
            .collect(),
        land: body
            .iter()
            .enumerate()
            .map(|(i, &b)| b == NONE && elevation[i] as f64 > sea)
            .collect(),
        lithology,
        laws: crate::planet_relief::Laws {
            stream_power_k_per_my: law("stream_power_k_per_my")?,
            diffusivity_m2_per_my: law("hillslope_diffusivity_m2_per_my")?,
            channel_head_length_m: head,
        },
    };
    let hydro = EndpointHydrology {
        body,
        fill: fill.into_iter().map(f64::from).collect(),
        records,
        provenance: json!({"artifact":h.digest,"terminalRule":"native water bodies and published receivers; no invented historical frame"}),
    };
    Ok((source, input, hydro))
}

#[derive(Clone, Serialize)]
pub struct BlockWork {
    pub outlet_cell: usize,
    pub core_cells: usize,
    pub halo_cells: usize,
    pub physical_nodes: usize,
}
/// Canonical process-node values. The transfer exporter may quantize to f32,
/// but must report that error and retain the native pole-cap convention.
#[derive(Serialize)]
pub struct ProcessBlock {
    pub schema: &'static str,
    pub model: &'static str,
    pub cell: u32,
    pub native_grid: [usize; 2],
    pub process_grid: [usize; 2],
    pub resolution: usize,
    pub apron: usize,
    pub n: usize,
    pub origin_post: [isize; 2],
    pub height_m: Vec<f64>,
    pub pole_height_m: [f64; 2],
    pub work: Vec<BlockWork>,
    pub process_receipts: serde_json::Value,
    pub source: serde_json::Value,
    pub limits: &'static str,
}

#[derive(Serialize)]
pub struct Receipt {
    pub model: &'static str,
    pub outlet_cell: usize,
    pub terminal: Terminal,
    pub core_cells: usize,
    pub halo_cells: usize,
    pub physical_nodes: usize,
    pub resolution: usize,
    pub north_south_spacing_m: f64,
    pub evolution_my: f64,
    pub eroded_nodes: usize,
    pub detached_bulk_m3: f64,
    pub maximum_incision_m: f64,
    pub maximum_deposition_m: f64,
    pub hillslope_boundary_import_m3: f64,
    pub hillslope_redistributed_m3: f64,
    pub maximum_diffusion_balance_residual: f64,
    pub maximum_diffusion_iterations: usize,
    pub geometric_bulk_relative_residual: f64,
    pub cross_native_cell_receivers: usize,
    pub divergent_coarse_receivers: usize,
    pub initially_filled: usize,
    pub initial_maximum_fill_m: f64,
    pub final_dry_filled_nodes: usize,
    pub final_maximum_dry_fill_m: f64,
    pub final_routing: pf_hydro_catchment::Audit,
    pub final_transport: crate::planet_transport::Budget,
    pub boundary_assumption: &'static str,
    pub sediment_account: &'static str,
    pub seconds: f64,
}
struct Surface {
    blocks: Vec<usize>,
    bed: Vec<f64>,
    receipt: Receipt,
}
pub struct World {
    input: EndpointInput,
    hydro: EndpointHydrology,
    outlet: Vec<usize>,
    cache: Mutex<BTreeMap<usize, Arc<Surface>>>,
    k: f64,
    diffusivity: f64,
    lithology: Vec<u8>,
}
impl World {
    pub fn open(root: &Path, source: &Inputs) -> R<Self> {
        let legacy = Input::read(root).map_err(|e| format!("catchment_input:{e}"))?;
        let hydro = Hydro::read(root, &legacy).map_err(|e| format!("catchment_hydrology:{e}"))?;
        let input = EndpointInput::from_legacy(&legacy).map_err(|e| e.to_string())?;
        Self::from_endpoint(input, EndpointHydrology::from_legacy(hydro), source)
    }
    pub fn from_endpoint(
        input: EndpointInput,
        hydro: EndpointHydrology,
        source: &Inputs,
    ) -> R<Self> {
        input.validate().map_err(|e| e.to_string())?;
        let n = input.grid.len();
        if [
            source.laws.stream_power_k_per_my,
            source.laws.diffusivity_m2_per_my,
            source.laws.channel_head_length_m,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err("invalid endpoint process laws".into());
        }
        if hydro.body.len() != n
            || hydro.fill.len() != n
            || hydro.fill.iter().any(|v| !v.is_finite())
            || hydro
                .body
                .iter()
                .any(|&b| b != NONE && !hydro.records.contains_key(&b))
            || hydro.records.iter().any(|(&id, b)| {
                id != b.body_id || b.seed_element >= n || !b.water_surface_m.is_finite()
            })
            || source.elevation_m.len() != n
            || source.relief_m.len() != n
            || source.receiver.len() != n
            || source.land.len() != n
            || source.lithology.len() != n
            || source.lithology.iter().any(|&l| l >= 12)
            || input
                .published_flow
                .iter()
                .zip(&source.receiver)
                .any(|(&a, &b)| a != if b == u32::MAX { NONE } else { b as usize })
        {
            return Err("endpoint/source shape, topology or lithology mismatch".into());
        }
        if input.grid.ny != source.rows
            || input.grid.nx != source.cols
            || input.grid.radius != source.radius_m
            || input.seed != source.seed
            || input
                .elevation
                .iter()
                .zip(&source.elevation_m)
                .any(|(a, b)| *a != f64::from(*b))
            || input
                .relief
                .iter()
                .zip(&source.relief_m)
                .any(|(a, b)| *a != f64::from(*b))
        {
            return Err("catchment_endpoint_differs_from_streamed_source".into());
        }
        // A receiving lake/ocean keeps its true body identity. Distinct coastal
        // entries are independent workload keys, not fictitious water bodies.
        let n = input.grid.len();
        let mut outlet = vec![NONE; n];
        for c in 0..n {
            let mut path = Vec::new();
            let mut i = c;
            while outlet[i] == NONE && hydro.body[i] == NONE && input.published_flow[i] != NONE {
                path.push(i);
                i = input.published_flow[i];
                if i >= n || path.len() > n {
                    return Err("catchment_native_receiver_cycle".into());
                }
            }
            let key = if outlet[i] == NONE { i } else { outlet[i] };
            outlet[i] = key;
            for c in path {
                outlet[c] = key;
            }
        }
        Ok(Self {
            input,
            hydro,
            outlet,
            cache: Mutex::new(BTreeMap::new()),
            k: source.laws.stream_power_k_per_my,
            diffusivity: source.laws.diffusivity_m2_per_my,
            lithology: source.lithology.clone(),
        })
    }
    fn bed(&self, c: usize, u: f64, v: f64) -> f64 {
        // native pipeline's declared coarse wet control volumes, not invented fine lakes.
        if self.hydro.body[c] != NONE {
            self.input.elevation[c]
        } else {
            let original_datum = self.input.grid.sample(&self.input.elevation, c, u, v);
            let reconstructed =
                crate::planet_initial::elevation(self.input.grid, &self.input.elevation, c, u, v);
            terrain(&self.input, c, u, v) - original_datum + reconstructed
        }
    }
    fn surface(&self, key: usize) -> R<Arc<Surface>> {
        // Serialize solves as well as publication: bounded resident memory and
        // request order cannot alter a domain's physical initial conditions.
        let mut cache = self.cache.lock().map_err(|_| "catchment_cache_poisoned")?;
        if let Some(s) = cache.get(&key) {
            return Ok(s.clone());
        }
        let s = Arc::new(self.build(key)?);
        if cache.len() >= 8 {
            if let Some(k) = cache.keys().next().copied() {
                cache.remove(&k);
            }
        }
        cache.insert(key, s.clone());
        Ok(s)
    }
    fn build(&self, key: usize) -> R<Surface> {
        let start = Instant::now();
        let g = self.input.grid;
        let r = RESOLUTION;
        let core: Vec<bool> = self.outlet.iter().map(|&v| v == key).collect();
        let mut included = core.clone();
        // Derived topology halo; no per-cell Dirichlet boundaries.
        for _ in 0..2 {
            let before = included.clone();
            for c in 0..g.len() {
                if before[c] {
                    for dy in -1..=1 {
                        for dx in -1..=1 {
                            if let Some(j) = g.neighbor(c, dx, dy) {
                                included[j] = true;
                            }
                        }
                    }
                }
            }
        }
        let cells: Vec<_> = (0..g.len()).filter(|&c| included[c]).collect();
        let n = cells.len() * r * r;
        if n > 12_000_000 {
            return Err(format!(
                "catchment_data_exceeds_memory_budget:outlet={key}:nodes={n}"
            ));
        }
        let mut blocks = vec![NONE; g.len()];
        let mut mesh = Mesh {
            bed: Vec::with_capacity(n),
            source: Vec::with_capacity(n),
            point: Vec::with_capacity(n),
            cell: Vec::with_capacity(n),
            edges: vec![vec![NONE; 8]; n],
            fixed: Vec::with_capacity(n),
            physical_len: n,
            radius: g.radius,
        };
        let mut areas = Vec::with_capacity(n);
        for (block, &c) in cells.iter().enumerate() {
            blocks[c] = block * r * r;
            for y in 0..r {
                for x in 0..r {
                    let (u, v) = ((x as f64 + 0.5) / r as f64, (y as f64 + 0.5) / r as f64);
                    let z = self.bed(c, u, v);
                    let b = self.hydro.body[c];
                    // A coarse body tag alone is not a submerged fine bed.
                    // native pipeline applies this same stage check; never raise a water
                    // stage to make an invalid Dirichlet condition pass.
                    let fixed = if b != NONE && z <= self.hydro.records[&b].water_surface_m {
                        Some((
                            self.hydro.records[&b].water_surface_m,
                            Terminal::WaterBody(b),
                        ))
                    } else {
                        None
                    };
                    let area = g.area(
                        c,
                        y as f64 / r as f64,
                        (y + 1) as f64 / r as f64,
                        1.0 / r as f64,
                    );
                    mesh.bed.push(z);
                    mesh.source.push(area * self.input.runoff_m_yr / YEAR_S);
                    mesh.point.push(g.point(c, u, v));
                    mesh.cell.push(c);
                    mesh.fixed.push(fixed);
                    areas.push(area);
                }
            }
        }
        // Declared native closed terminals remain closed; a minimum is selected
        // from actual bed, not a placed anchor or forced path to the ocean.
        for &c in &cells {
            if self.hydro.body[c] == NONE && self.input.published_flow[c] == NONE {
                let lo = blocks[c];
                let i = (lo..lo + r * r)
                    .min_by(|&a, &b| mesh.bed[a].total_cmp(&mesh.bed[b]))
                    .unwrap();
                mesh.fixed[i] = Some((mesh.bed[i], Terminal::ClosedFloor(c)));
            }
        }
        let width = g.nx * r;
        let height = g.ny * r;
        let mut ghosts = BTreeMap::new();
        for i in 0..n {
            let c = mesh.cell[i];
            let local = i - blocks[c];
            let gx = c % g.nx * r + local % r;
            let gy = c / g.nx * r + local / r;
            let mut slot = 0;
            for dy in -1isize..=1 {
                for dx in -1isize..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let s = slot;
                    slot += 1;
                    let yy = gy as isize + dy;
                    if yy < 0 || yy >= height as isize {
                        continue;
                    }
                    let xx = (gx as isize + dx).rem_euclid(width as isize) as usize;
                    let yy = yy as usize;
                    let cc = yy / r * g.nx + xx / r;
                    let id = yy * width + xx;
                    let j = if blocks[cc] != NONE {
                        blocks[cc] + yy % r * r + xx % r
                    } else if let Some(&j) = ghosts.get(&id) {
                        j
                    } else {
                        let (u, v) = (
                            (xx % r) as f64 / r as f64 + 0.5 / r as f64,
                            (yy % r) as f64 / r as f64 + 0.5 / r as f64,
                        );
                        let z = self.bed(cc, u, v);
                        let head = z.max(g.sample(&self.hydro.fill, cc, u, v));
                        let j = mesh.bed.len();
                        mesh.bed.push(z);
                        mesh.source.push(0.0);
                        mesh.point.push(g.point(cc, u, v));
                        mesh.cell.push(cc);
                        mesh.fixed.push(Some((head, Terminal::Boundary(id))));
                        mesh.edges.push(Vec::new());
                        ghosts.insert(id, j);
                        j
                    };
                    mesh.edges[i][s] = j;
                    if j >= n {
                        mesh.edges[j].push(i);
                    }
                }
            }
        }
        let initial = mesh.bed.clone();
        let mut routing =
            pf_hydro_catchment::route(&mesh).map_err(|e| format!("catchment_route:{e}"))?;
        let mut flow = crate::planet_transport::Flow::build(&mesh, &routing)?;
        let initially_filled = routing.audit.filled_nodes;
        let initial_maximum_fill_m = routing.audit.max_fill_m;
        let rocks = crate::t1b_litho::LithologyV1::standin_v1();
        let hillslope = crate::planet_hillslope::Operator::new(
            &mesh,
            &areas,
            g.ny * r,
            g.nx * r,
            self.diffusivity,
            DT_MY,
        )?;
        let mut detached = 0.0;
        let mut boundary_import = 0.0;
        let mut redistributed = 0.0;
        let mut diffusion_residual = 0.0_f64;
        let mut diffusion_iterations = 0;
        for _ in 0..STEPS {
            let mut next = mesh.bed.clone();
            for &i in &routing.order {
                if i >= n || mesh.fixed[i].is_some() || routing.head[i] - mesh.bed[i] > 0.002 {
                    continue;
                }
                let area = if self.input.runoff_m_yr > 0.0 {
                    flow.discharge[i] * YEAR_S / self.input.runoff_m_yr
                } else {
                    0.0
                };
                let mut gradient_weight = 0.0;
                let mut receiver_term = 0.0;
                for (slot, &j) in mesh.edges[i].iter().enumerate() {
                    let w = flow.weights[i][slot];
                    if w == 0.0 {
                        continue;
                    }
                    // Only actual downhill bed or known terminal-stage gradients
                    // erode rock. No dry filled head is substituted as a bed.
                    let base = mesh.fixed[j].as_ref().map_or(next[j], |(h, _)| *h);
                    if base >= mesh.bed[i] {
                        continue;
                    }
                    let d = distance(mesh.point[i], mesh.point[j], g.radius);
                    gradient_weight += w / d;
                    receiver_term += w * base / d;
                }
                let factor = DT_MY * self.k * rocks.k(self.lithology[mesh.cell[i]]) * area.sqrt();
                next[i] = (mesh.bed[i] + factor * receiver_term) / (1.0 + factor * gradient_weight);
            }
            for i in 0..n {
                detached += (mesh.bed[i] - next[i]) * areas[i];
            }
            mesh.bed = next;
            let budget = hillslope.evolve(&mut mesh)?;
            boundary_import += budget.boundary_import_m3;
            redistributed += budget.redistributed_m3;
            diffusion_residual = diffusion_residual.max(budget.relative_residual);
            diffusion_iterations = diffusion_iterations.max(budget.iterations);
            routing = pf_hydro_catchment::route(&mesh)
                .map_err(|e| format!("catchment_evolved_route:{e}"))?;
            flow = crate::planet_transport::Flow::build(&mesh, &routing)?;
        }
        let mut bed_change = 0.0;
        let mut max_deposition = 0.0_f64;
        let mut max_cut: f64 = 0.0;
        let mut eroded = 0;
        let mut cross = 0;
        let mut divergent = 0;
        for i in 0..n {
            let cut = initial[i] - mesh.bed[i];
            bed_change -= cut * areas[i];
            max_deposition = max_deposition.max(-cut);
            if cut > 0.0 {
                eroded += 1;
                max_cut = max_cut.max(cut);
            }
            let j = routing.receiver[i];
            if j != NONE && mesh.cell[i] != mesh.cell[j] {
                cross += 1;
                if self.input.published_flow[mesh.cell[i]] != mesh.cell[j] {
                    divergent += 1;
                }
            }
        }
        let bulk_residual = (bed_change + detached - boundary_import).abs()
            / detached.abs().max(boundary_import.abs()).max(1.0);
        if bulk_residual > 1e-4 {
            return Err(format!("catchment_bulk_conservation:{bulk_residual}"));
        }
        let terminal = if self.hydro.body[key] != NONE {
            Terminal::WaterBody(self.hydro.body[key])
        } else {
            Terminal::ClosedFloor(key)
        };
        let core_cells = core.iter().filter(|v| **v).count();
        let dry_fills: Vec<_> = (0..n)
            .filter(|&i| mesh.fixed[i].is_none())
            .map(|i| routing.head[i] - mesh.bed[i])
            .filter(|&v| v > 0.002)
            .collect();
        let receipt=Receipt {model:MODEL,outlet_cell:key,terminal,core_cells,halo_cells:cells.len()-core_cells,physical_nodes:n,
            resolution:r,north_south_spacing_m:g.radius*std::f64::consts::PI/(g.ny*r) as f64,evolution_my:STEPS as f64*DT_MY,
            eroded_nodes:eroded,detached_bulk_m3:detached,maximum_incision_m:max_cut,cross_native_cell_receivers:cross,
            maximum_deposition_m:max_deposition,hillslope_boundary_import_m3:boundary_import,hillslope_redistributed_m3:redistributed,
            maximum_diffusion_balance_residual:diffusion_residual,maximum_diffusion_iterations:diffusion_iterations,
            geometric_bulk_relative_residual:bulk_residual,
            divergent_coarse_receivers:divergent,initially_filled,initial_maximum_fill_m,
            final_dry_filled_nodes:dry_fills.len(),final_maximum_dry_fill_m:dry_fills.iter().copied().fold(0.0,f64::max),final_routing:routing.audit,final_transport:flow.budget,
            boundary_assumption:"two native-cell halo; zero incoming external flux; endpoint outer routing heads fixed; core not constrained to old receivers",
            sediment_account:"stream-power detachment and conservative hillslope bulk flux in m3; fluvial sediment destinations and layer chemistry transport unresolved; no canonical stock is created",
            seconds:start.elapsed().as_secs_f64()};
        mesh.bed.truncate(n);
        Ok(Surface {
            blocks,
            bed: mesh.bed,
            receipt,
        })
    }
    fn node(&self, x: isize, y: usize) -> R<f64> {
        let g = self.input.grid;
        let r = RESOLUTION;
        let x = x.rem_euclid((g.nx * r) as isize) as usize;
        let y = y.min(g.ny * r - 1);
        let c = y / r * g.nx + x / r;
        if self.hydro.body[c] != NONE
            && self.input.elevation[c] <= self.hydro.records[&self.hydro.body[c]].water_surface_m
        {
            return Ok(self.input.elevation[c]);
        }
        let s = self.surface(self.outlet[c])?;
        Ok(s.bed[s.blocks[c] + y % r * r + x % r])
    }
    pub fn sample(&self, lat: f64, lon: f64) -> R<f64> {
        if !lat.is_finite() || !lon.is_finite() || lat.abs() > 90.0 || lon.abs() > 180.0 {
            return Err("planet_position_invalid".into());
        }
        let g = self.input.grid;
        let r = RESOLUTION;
        let x = (lon + 180.0).rem_euclid(360.0) / 360.0 * (g.nx * r) as f64 - 0.5;
        let y = (lat + 90.0) / 180.0 * (g.ny * r) as f64 - 0.5;
        let i = x.floor() as isize;
        let tx = x - x.floor();
        let row =
            |y: usize| -> R<f64> { Ok((1.0 - tx) * self.node(i, y)? + tx * self.node(i + 1, y)?) };
        // Collapse polar longitude to one native datum, as in the atlas.
        if y < 0.0 || y > (g.ny * r - 1) as f64 {
            let north = y >= 0.0;
            let coarse = if north { g.ny - 1 } else { 0 };
            let pole = self.input.elevation[coarse * g.nx..(coarse + 1) * g.nx]
                .iter()
                .sum::<f64>()
                / g.nx as f64;
            let mix = if north {
                2.0 * (y - (g.ny * r - 1) as f64)
            } else {
                -2.0 * y
            };
            return Ok(pole * mix + row(if north { g.ny * r - 1 } else { 0 })? * (1.0 - mix));
        }
        let j = y.floor() as usize;
        let ty = y - y.floor();
        Ok(row(j)? * (1.0 - ty) + row((j + 1).min(g.ny * r - 1))? * ty)
    }
    pub fn receipts(&self, cells: &BTreeSet<u32>) -> R<serde_json::Value> {
        let keys: BTreeSet<_> = cells.iter().map(|&c| self.outlet[c as usize]).collect();
        let cache = self.cache.lock().map_err(|_| "catchment_cache_poisoned")?;
        serde_json::to_value(
            keys.iter()
                .filter_map(|k| cache.get(k).map(|s| &s.receipt))
                .collect::<Vec<_>>(),
        )
        .map_err(|e| e.to_string())
    }
    /// Work estimation is a resource guard, never a geographic admission list.
    pub fn block_work(&self, cell: u32) -> R<Vec<BlockWork>> {
        let g = self.input.grid;
        if cell as usize >= g.len() {
            return Err("process block cell outside source grid".into());
        }
        let r = RESOLUTION;
        let n = r + 2;
        let x0 = (cell as usize % g.nx * r) as isize - 1;
        let y0 = (cell as usize / g.nx * r) as isize - 1;
        let mut keys = BTreeSet::new();
        for j in 0..n {
            for i in 0..n {
                let x = (x0 + i as isize).rem_euclid((g.nx * r) as isize) as usize;
                let y = (y0 + j as isize).clamp(0, (g.ny * r - 1) as isize) as usize;
                let c = y / r * g.nx + x / r;
                if self.hydro.body[c] == NONE
                    || self.input.elevation[c]
                        > self.hydro.records[&self.hydro.body[c]].water_surface_m
                {
                    keys.insert(self.outlet[c]);
                }
            }
        }
        let mut result = Vec::new();
        for key in keys {
            let core: Vec<bool> = self.outlet.iter().map(|&v| v == key).collect();
            let mut included = core.clone();
            for _ in 0..2 {
                let before = included.clone();
                for c in 0..g.len() {
                    if before[c] {
                        for dy in -1..=1 {
                            for dx in -1..=1 {
                                if let Some(j) = g.neighbor(c, dx, dy) {
                                    included[j] = true;
                                }
                            }
                        }
                    }
                }
            }
            let core_cells = core.iter().filter(|&&x| x).count();
            let count = included.iter().filter(|&&x| x).count();
            result.push(BlockWork {
                outlet_cell: key,
                core_cells,
                halo_cells: count - core_cells,
                physical_nodes: count * r * r,
            });
        }
        Ok(result)
    }
    pub fn process_block(&self, cell: u32, max_work_nodes: usize) -> R<ProcessBlock> {
        let work = self.block_work(cell)?;
        let total = work
            .iter()
            .try_fold(0usize, |a, w| a.checked_add(w.physical_nodes))
            .ok_or("process block workload overflow")?;
        if max_work_nodes == 0 || total > max_work_nodes {
            return Err(format!(
                "process block requires {total} native solve nodes; caller budget {max_work_nodes}"
            ));
        }
        // Pin every owning solve needed by the apron during this operation.
        // Otherwise a cache eviction could repeatedly regenerate alternating
        // boundary owners. Never sample a neighbour's halo as its own solution.
        let mut surfaces = BTreeMap::new();
        for w in &work {
            surfaces.insert(w.outlet_cell, self.surface(w.outlet_cell)?);
        }
        let g = self.input.grid;
        let r = RESOLUTION;
        let n = r + 2;
        let origin = [
            (cell as usize % g.nx * r) as isize - 1,
            (cell as usize / g.nx * r) as isize - 1,
        ];
        let mut height = Vec::with_capacity(n * n);
        for j in 0..n {
            for i in 0..n {
                let x = (origin[0] + i as isize).rem_euclid((g.nx * r) as isize) as usize;
                let y = (origin[1] + j as isize).clamp(0, (g.ny * r - 1) as isize) as usize;
                let c = y / r * g.nx + x / r;
                let z = if self.hydro.body[c] != NONE
                    && self.input.elevation[c]
                        <= self.hydro.records[&self.hydro.body[c]].water_surface_m
                {
                    self.input.elevation[c]
                } else {
                    let s = surfaces
                        .get(&self.outlet[c])
                        .ok_or("canonical block owner missing")?;
                    let b = *s.blocks.get(c).ok_or("canonical block address missing")?;
                    if b == NONE {
                        return Err("canonical block tried to use an unrelated halo".into());
                    }
                    s.bed[b + y % r * r + x % r]
                };
                if !z.is_finite() {
                    return Err("nonfinite canonical block height".into());
                }
                height.push(z);
            }
        }
        let pole = |row: usize| {
            self.input.elevation[row * g.nx..(row + 1) * g.nx]
                .iter()
                .sum::<f64>()
                / g.nx as f64
        };
        let receipts =
            serde_json::to_value(surfaces.values().map(|s| &s.receipt).collect::<Vec<_>>())
                .map_err(|e| e.to_string())?;
        Ok(ProcessBlock {
            schema: "maiforge.native-process-block/v1",
            model: MODEL,
            cell,
            native_grid: [g.ny, g.nx],
            process_grid: [g.ny * r, g.nx * r],
            resolution: r,
            apron: 1,
            n,
            origin_post: origin,
            height_m: height,
            pole_height_m: [pole(0), pole(g.ny - 1)],
            work,
            process_receipts: receipts,
            source: self.provenance(),
            limits: LIMITS,
        })
    }

    pub fn provenance(&self) -> serde_json::Value {
        serde_json::json!({"input":self.input.provenance,"hydrology":self.hydro.provenance,
            "initial_elevation_law":crate::planet_initial::LAW,
            "initial_elevation_meaning":"native finite-volume endpoint heights treated as centre observations; no new observations, no post-solve tile blending; native wet floors unchanged"})
    }
}
