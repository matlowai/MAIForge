//! T1b lane VE (2026-10-15; `unpublished-design-record`): eruption styles and a keyed eruption history
//! for the volcano edifices of lanes CO / CF. LABEL: pre-signoff, fidelity tier Preview, derived, never canon; a construction
//! (no SimTime). Everything is opt-in under `params.ls.tb.tc.ve` (inside params.ls.tb, so lane TB's paired keying keeps every
//! other keyed draw); absent = lane TH's behaviour byte-identical.
//!
//! Reduced-order model (design 2026-10-15: "consider real volcanos and their types of erruptions"; physics-first rule):
//! - the record kind selects an eruptive family (composite stratovolcano, basaltic shield, rift fissure field, flood-basalt
//!   plateau); the record gives no per-eruption detail, so the youngest mappable units are drawn as a keyed statistical
//!   closure over the recorded activity span (counts, VEI and size distributions declared in the rule strings);
//! - each unit is a process-shaped deposit or cut: lava flows that descend the edifice and are captured by the network
//!   valleys downstream (flat valley fill confined by the valley walls, lobate fronts at angle of repose), cinder cones with
//!   their measured shape ratios, pyroclastic aprons (gentle fan surface filling the lower flank and ring plain, exponential
//!   fall thinning), lahar valley fills, lava domes, craters, calderas, sector-collapse amphitheatres, maars with tuff rings,
//!   fissure ramparts, pit craters, stacked plateau sheets;
//! - every unit is lowered by the record-age erosion law (E at the vent x unit age); units are applied oldest first, the
//!   exposed unit is the youngest of mappable thickness;
//! - evaluation is a pure function of chart position and the post's own pre-unit height (windows agree exactly).

use crate::t0::T0;
use crate::t1b::{CHANNEL, FREE, GUIDE, LAKE, MARINE};
use crate::t1b_coast::{ang_diff, bil_at, cone_cf_opt, seg_nearest, smoothstep01, CfParamsV1, CoastCtx, VolcanoParamsV1, VolcanoSite, OVER_EDIFICE};
use crate::t1b_rivernet::{Curve, EvalLat};
use crate::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use vz_keyed_rng::Scope;
use vz_tilekey::chart::Chart;

pub const PURPOSE_VE: &str = "t1b.volcano.ve";
/// edifice_style codes (u8).
pub const ST_NONE: u8 = 0;
pub const ST_LAVA: u8 = 1;
pub const ST_FRONT: u8 = 2;
pub const ST_CINDER: u8 = 3;
pub const ST_APRON: u8 = 4;
pub const ST_LAHAR: u8 = 5;
pub const ST_DOME: u8 = 6;
pub const ST_WALL: u8 = 7;
pub const ST_SHIELD: u8 = 8;
pub const ST_FISSURE: u8 = 9;
pub const ST_FLOOR: u8 = 10;
pub const ST_TUFF: u8 = 11;
pub const ST_BODY: u8 = 12;
pub const ST_PLATEAU: u8 = 13;
pub const ST_NECK: u8 = 14;
/// edifice_unit sentinels (u16).
pub const UNIT_BODY: u16 = 0;
pub const UNIT_NONE: u16 = u16::MAX;
/// Path sample levels: drape on the edifice body / drape on open ground (otherwise the captured valley bed).
pub const LV_BODY: f64 = -1.0e9;
pub const LV_OPEN: f64 = -2.0e9;

pub const STYLE_NOTE: &str = "lane VE (T1b, opt-in ve): the volcanic unit exposed at the post, by eruption style: 0 none (not an edifice post, or water), 1 lava flow (andesite/dacite flow body), 2 lobe front (snout of a lava flow past its end point plus its last min(1.5 w, 400 m): the steep lobate front zone), 3 cinder (scoria) cone, 4 pyroclastic apron (fall and pyroclastic-density-current fan), 5 lahar (volcaniclastic valley fill / fan), 6 lava dome, 7 crater / caldera / amphitheatre / maar wall, 8 shield flow (basaltic pahoehoe / aa; also the shield body), 9 fissure / spatter rampart, 10 crater / caldera / maar / pit floor (closed depressions are filled to their spill level by T1b's final fill: a declared lake surface), 11 maar tuff ring, 12 composite body (undifferentiated older flows and breccias), 13 plateau flow (flood basalt sheet), 14 neck / dike (residual of an eroded edifice). The header's prep.net.volcanoes[k].ve lists every unit (kind, style, age, VEI, size). A unit is exposed where it adds (or cuts) at least ve.t_map_mm and no younger unit does; wave-cut edifice posts (over 4) keep their unit";
pub const UNIT_NOTE: &str = "lane VE (T1b, opt-in ve): index of the exposed volcanic unit in the volcano's unit list (header prep.net.volcanoes[volcano_id].ve.units[k - 1] for k >= 1; the volcano is the post's volcano_id); 0 = the edifice body (older undifferentiated edifice or a residual neck); 65535 = not an edifice post, or water";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct VeParamsV1 {
    /// 1 = eruption styles and keyed history on (0 = lane TH behaviour; the struct is then normally absent).
    pub on: u32,
    /// Mappable thickness: the exposed unit at a post is the youngest that adds or cuts at least this.
    pub t_map_mm: i64,
    /// Lava / lahar path step, valley capture radius, stop distance from T1a curves.
    pub path_step_mm: i64,
    pub capture_mm: i64,
    pub t1a_stop_mm: i64,
    /// CF gullies (composite) cut the body and the units older than this; younger units are not gullied.
    pub gully_age_yr: i64,
    /// Units per volcano cap (all families).
    pub max_units: u32,
    pub style_rule: String,
    pub history_rule: String,
    pub size_rule: String,
    pub shape_rule: String,
    pub path_rule: String,
    pub deposit_rule: String,
    /// Lane VE2 options (absent = ve=1 byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ve2: Option<Ve2ParamsV1>,
    /// Lane VE3 options (absent = ve2 byte-identical; needs ve2 steer).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ve3: Option<Ve3ParamsV1>,
    /// Lane VE4 options (absent = ve3 byte-identical; needs ve3 rough).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ve4: Option<Ve4ParamsV1>,
    /// Lane VE5 options (absent = ve4 byte-identical; needs ve4 lobes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ve5: Option<Ve5ParamsV1>,
    /// Lane VH-I (opt-in `ve_wind=1`, window-only; absent = byte-identical): tephra dispersal on the T0 wind.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wind: Option<VeWindParamsV1>,
    /// CX-EDIFICE-EREC: record-rate degradation, layer AND windows. Absent = old bytes and reads.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub erec: Option<VeErecParamsV1>,
    /// CX-CONE-STARBURST: native drainage-driven flank, layer AND windows.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flank: Option<FlankParamsV1>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flank2: Option<Flank2Params>,
}

/// Continuous constructional envelope plus finite-duration detachment-limited erosion.
/// Sources support the functional laws, not the project calibration constants.
pub const FLANK_RULE: &str = "ve_flank_v1: retain CF exponential concave envelope z=top-s_top*lambda*(1-exp(-(r-rc)/lambda))*scale and crater (exponential stratovolcano profiles: Cosburn and Roy 2020, doi:10.1016/j.jvolgeores.2020.107123). Replace CO straight V-sector cuts and CF young template gullies with edifice-local priority flood on a periodic polar mesh, centred metric gradient and angular two-receiver partition (partition follows Tarboton 1997 doi:10.1029/96WR03137; polar discretization is this implementation) on the envelope plus inherited CF lava emplacement lobes. Uniform unit runoff gives conserved accumulated area. For erosion A is the contributing discharge-area integrated across a fixed128m-wide radial gate (fractional end cells, only flux crossing radial faces), a declared physical support width [H] independent of mesh size; drainage joins by actual receivers, never keyed gully angles. Sixty receiver-first implicit steps over 30000 yr, E=K*A^0.5*S, K=20/My [H] (Whipple and Tucker 1999 doi:10.1029/1999JB900120; Braun and Willett 2013 doi:10.1016/j.geomorph.2012.10.008). Incision begins smoothly at 0.05-0.10 km2 [H]. Bank retreat and deposition use conservative finite-volume hillslope transport dz/dt=div(kappa grad z), implicit radial and periodic angular solves, zero normal flux at the protected crater/foot. kappa=b(A)^2/(2T)*(K/20), b=32m+80m*(A/1km2)^0.4; b is the diffusion length over T=30ky [H, hydraulic-geometry motivated bank mobility, not a measured soil diffusivity; Leopold and Maddock 1953 doi:10.3133/pp252]. Internal edge flux is shared by adjacent polar cells; crept material is deposited on the flank and recorded separately from stream-power export. A fixed192m crater/foot buffer and128m transition preserve the physical domain under refinement. Stream-detached volume is exported to local routing terminals; hillslope-transported volume stays as deposition on the edifice, no downstream sediment or groundwater dynamics invented; uniform runoff and rock export are audited. Global record-age scale remains the inherited bulk-degradation initial condition; this declared 30ky interval is sub-record flank incision, not 2My-resolved eruption history. 64m radial spacing, angular spacing <=64m at the foot, chart-independent vent-centred polar grid; C2 cubic B-spline reconstruction of the signed erosion correction; actual VE lava units and roughness retained. No new random draws. Numerical grid and calibration are explicit model limits.";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct FlankParamsV1 { pub rule: String }
impl FlankParamsV1 { pub fn v1() -> Self { Self { rule: FLANK_RULE.into() } } }

pub const FLANK2_SURFACE_RULE: &str = "ve_flank2_surface_v1: tangent exponential native apron; thickness-dependent terrain contact; existing native VE units and retained inferred buried-flow geometry; no Arc1 polar erosion";
/// Surface mode emits a typed, hashed native reference extract. Every retained
/// byte comes from WindowRun; checks still run on the full native surface.
/// This avoids storing unrelated CT/material fields in each intermediate.
pub fn publish_flank2_surface(
    p: &crate::t1b_rivernet::RiverNetParamsV1,
    run: &crate::t1b_window::WindowRun,
    out: &std::path::Path,
) -> R<bool> {
    if !ve(p)
        .and_then(|v| v.flank2.as_ref())
        .is_some_and(|f| f.phase == 0)
    {
        return Ok(false);
    }
    let keep = [
        "height_mm",
        "fixed",
        "runoff_weight",
        "volcano_id",
        "edifice_style",
        "edifice_unit",
        "drainage_area_m2",
        "process_geometry_fixed", "process_native_link", "process_edifice_style", "process_edifice_unit",
    ];
    let mut header = serde_json::to_value(&run.header).map_err(|e| e.to_string())?;
    header["schema"] = "vz.flank2-native-surface/v1".into();
    if run.header.recipe.upstream.iter().any(|(name,_)|name=="ls.native_registered_channels_candidate") {header["schema"]="vz.surface-owner-construction/v1".into();}
    header["native_source_window_digest"] = run.raw.digest.clone().into();
    let fields = header["fields"]
        .as_array_mut()
        .ok_or("reference field descriptors")?;
    fields.retain(|f| keep.contains(&f["name"].as_str().unwrap_or("")));
    let files: std::collections::BTreeSet<String> = fields
        .iter()
        .map(|f| f["file"].as_str().unwrap().to_string())
        .collect();
    std::fs::create_dir_all(out.join("fields")).map_err(|e| e.to_string())?;
    for (file, bytes) in &run.blobs {
        if files.contains(file) {
            std::fs::write(out.join(file), bytes).map_err(|e| e.to_string())?;
        }
    }
    let digest = pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())?;
    std::fs::write(
        out.join("window.json"),
        serde_json::to_vec_pretty(&serde_json::json!({"header":header,"digest":digest}))
            .map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    Ok(true)
}

/// Arc2 is a reference-bound native terrain process, not an angular template.
pub const FLANK2_RULE: &str = "ve_flank2_v1: surface mode publishes the native constructed cone, actual VE units and terrain contacts. Evolve mode requires those hashed native32m windows from the same T0/site/realization; it routes on their heights and uses their runoff/fixed fields. No polar mesh, radial head gate or fixed-width discharge averaging. Freeman-type MFD S^1.1 on the native priority-flood DAG (Freeman1991; routing sensitivity Shelef2013 doi:10.1002/jgrf.20127). Channel initiation sqrt(A)*S>12m (Montgomery/Dietrich1989 doi:10.1029/WR025i008p01907 motivates slope-area dependence, not this coefficient). Detachment E=32/My*max(sqrt(A)*S-12m,0) times a lava-resistance factor; finite volcanic substrate supply. Flood width=3.5m*sqrt(Q/(1m3/s)), Q=A*(5m3/s/km2) [H effective flood], valley banks relax toward material repose after incision. Annual sediment capacity=0.02m/yr*A*S^2 [H], with supply-limited erosion and deposition when load exceeds capacity; fans therefore respond to actual slope drops. Dry crater floors have zero erosible stock but remain flow-connected. All water and sediment reach native fixed outlets or the reference-domain boundary; no discharge clipping. Historical lava paths are inferred native sub-record construction contacts, not observed buried strata or fresh deposits.30ky in120steps [H], detachment of exposed native units scaled by min(unit_age/30ky,1), an exposure-time approximation rather than chronological emplacement. The exponential core joins a tangent exponential apron at20percent remaining height, blended C1 over0.1Rf; actual terrain burial/drape weight depends on deposited thickness, not a radius band. Parent ve_flank alone retains Arc1 exactly. Limits:32m resolution, inferred buried contacts, effective flood and sediment coefficients; native VE scale split: regional E/F construction fit is retained as initial condition, because that solve sees body only. The process increment is applied with actual eruption units in window evaluation; literal whole-cell F6 still checks the final surface. Later local T1b constraints can modify the increment.";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Flank2Params {
    pub rule: String,
    /// 0 native construction surface; 1 process its supplied reference.
    pub phase: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference: Option<Flank2Reference>,
    #[serde(skip)]
    pub runtime: Flank2Runtime,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Flank2Reference {
    pub digest: String,
    pub windows: u32,
    pub input_recipes: Vec<String>,
    pub ledgers: serde_json::Value,
}
#[derive(Clone, Default)]
pub struct Flank2Runtime(std::sync::Arc<Vec<Flank2Patch>>);
impl std::fmt::Debug for Flank2Runtime {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Flank2Runtime({} patches)", self.0.len())
    }
}
impl PartialEq for Flank2Runtime {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}
impl Eq for Flank2Runtime {}
#[derive(Clone, Debug)]
struct Flank2Patch {
    gx: i64,
    gy: i64,
    n: usize,
    delta: Vec<f32>,
}
#[derive(Clone, Debug, Serialize)]
pub struct Flank2Site {
    pub history: Vec<VeUnit>,
    pub history_rule: String,
}

/// A tangent exponential apron gives continuous radial slope and no finite-radius foot.
pub fn flank2_height(v: &VolcanoSite, r: f64) -> f64 {
    flank2_profile(
        v.h_m,
        v.denudation_m,
        v.s_top,
        v.lambda_m,
        v.crater_radius_m,
        v.r_main_m,
        r,
    )
}
fn flank2_profile(h0: f64, d: f64, s: f64, lambda: f64, rc: f64, rf: f64, r: f64) -> f64 {
    let scale = ((h0 - d) / h0).max(0.0);
    let h = h0 * scale;
    let rho = (r - rc).max(0.0);
    let c = s * lambda * scale;
    if h <= 0.0 {
        return 0.0;
    }
    let hj = 0.2 * h;
    let e = ((hj + c - h) / c).clamp(1e-12, 1.0);
    let join = -lambda * e.ln();
    let sj = s * scale * e;
    let len = hj / sj.max(1e-12);
    let core = h - c * (1.0 - (-rho / lambda).exp());
    let tail = hj * (-(rho - join) / len).exp();
    let blend = smoothstep01((rho - join + 0.05 * rf) / (0.1 * rf));
    ((1.0 - blend) * core + blend * tail).max(0.0)
}

fn flank2_bytes(p: &std::path::Path) -> R<Vec<u8>> {
    if p.exists() {
        return std::fs::read(p).map_err(|e| e.to_string());
    }
    let gz = std::path::PathBuf::from(format!("{}.gz", p.display()));
    let o = std::process::Command::new("gzip")
        .arg("-dc")
        .arg(&gz)
        .output()
        .map_err(|e| e.to_string())?;
    if !o.status.success() {
        return Err(format!(
            "cannot read {}: {}",
            gz.display(),
            String::from_utf8_lossy(&o.stderr)
        ));
    }
    Ok(o.stdout)
}

impl Flank2Params {
    pub fn bind(
        mode: &str,
        root: Option<&str>,
        t0: &T0,
        site_path: &std::path::Path,
        realization: u32,
        configuration: &serde_json::Value,
        cache: Option<&str>,
    ) -> R<Self> {
        if mode == "surface" {
            if root.is_some() {
                return Err("surface mode does not consume a reference".into());
            }
            return Ok(Self {
                rule: FLANK2_SURFACE_RULE.into(),
                phase: 0,
                reference: None,
                runtime: Default::default(),
            });
        }
        if mode != "1" {
            return Err("expected ve_flank2=0|surface|1".into());
        }
        let root = std::path::Path::new(
            root.ok_or("ve_flank2=1 needs ve_flank2_ref=<native surface windows>")?,
        );
        let site = crate::t0::SiteRecord::load(site_path)?;
        struct Tile {
            x: i64,
            y: i64,
            gx: i64,
            gy: i64,
            n: usize,
            h: Vec<f64>,
            fix: Vec<u8>,
            rain: Vec<f32>,
            vid: Vec<u8>,
            style: Vec<u8>,
            unit: Vec<u16>,
            volcs: Vec<serde_json::Value>,
        }
        let mut tiles = Vec::new();
        let mut hasher = blake3::Hasher::new();
        let mut recipes = std::collections::BTreeSet::new();
        let mut paths: Vec<_> = std::fs::read_dir(root)
            .map_err(|e| e.to_string())?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.is_dir() && p.file_name().unwrap().to_string_lossy().starts_with("w_"))
            .collect();
        paths.sort();
        for p in paths {
            let bytes = flank2_bytes(&p.join("window.json"))?;
            hasher.update(&bytes);
            let j: serde_json::Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
            let h = &j["header"];
            if j["digest"]
                != pf_contracts::artifact::digest_canonical(h).map_err(|e| e.to_string())?
            {
                return Err("reference header hash mismatch".into());
            }
            let params = &h["recipe"]["params"]["rivernet"];
            for key in [
                "confinement",
                "volcano",
                "lithology",
                "channel_head_area_m2",
                "solve_spacing_mm",
                "margin_mm",
                "realization",
            ] {
                if params[key] != configuration[key] {
                    return Err(format!("reference construction parameter {key} mismatch"));
                }
            }

            let mut reference_ve = params["ls"]["tb"]["tc"]["ve"].clone();
            let mut current_ve = configuration["ls"]["tb"]["tc"]["ve"].clone();
            reference_ve
                .as_object_mut()
                .ok_or("reference VE parameters")?
                .remove("flank2");
            current_ve
                .as_object_mut()
                .ok_or("current VE parameters")?
                .remove("flank2");
            if reference_ve != current_ve {
                return Err("reference volcanic construction parameters mismatch".into());
            }
            let fp = &params["ls"]["tb"]["tc"]["ve"]["flank2"];
            if fp["phase"] != 0 || fp["rule"] != FLANK2_SURFACE_RULE {
                return Err("reference is not this native construction-stage law".into());
            }
            if h["recipe"]["seed_scope"]["seed"].as_u64() != Some(site.seed)
                || params["realization"].as_u64() != Some(realization as u64)
            {
                return Err("reference seed/realization mismatch".into());
            }
            let upstream = h["recipe"]["upstream"]
                .as_array()
                .ok_or("reference upstream missing")?;
            for (role, want) in [
                ("site", &site.bytes_digest),
                ("t0.terrain", &t0.terrain_digest),
                ("t0.hydrology", &t0.hydrology_digest),
                ("t0.subgrid_relief_sidecar", &t0.subgrid_digest),
            ] {
                if !upstream.iter().any(|a| a[0] == role && a[1] == *want) {
                    return Err(format!("reference {role} mismatch"));
                }
            }
            recipes.insert(
                h["recipe_digest"]
                    .as_str()
                    .ok_or("reference recipe digest")?
                    .to_string(),
            );
            let lat = &h["lattice"];
            let n = lat["n"].as_u64().ok_or("reference n")? as usize;
            if n != 641 || lat["halo_posts"] != 64 {
                return Err("reference must use native32m641-post lattice".into());
            }
            let field = |name: &str| -> R<Vec<u8>> {
                let f = h["fields"]
                    .as_array()
                    .ok_or("reference fields")?
                    .iter()
                    .find(|f| f["name"] == name)
                    .ok_or_else(|| format!("reference field {name} missing"))?;
                let b = flank2_bytes(&p.join(f["file"].as_str().ok_or("field file")?))?;
                if f["digest"] != format!("blake3:{}", blake3::hash(&b).to_hex()) {
                    return Err(format!("reference {name} hash mismatch"));
                }
                Ok(b)
            };
            let z = field("height_mm")?;
            let fix = field("fixed")?;
            let rw = field("runoff_weight")?;
            let vid = field("volcano_id")?;
            let style = field("edifice_style")?;
            let unit = field("edifice_unit")?;
            if z.len() != n * n * 4
                || fix.len() != n * n
                || rw.len() != n * n * 4
                || vid.len() != n * n
                || style.len() != n * n
                || unit.len() != n * n * 2
            {
                return Err("reference field sizes".into());
            }
            for b in [&z, &fix, &rw, &vid, &style, &unit] {
                hasher.update(b);
            }
            tiles.push(Tile {
                x: lat["rx"].as_i64().ok_or("rx")?,
                y: lat["ry"].as_i64().ok_or("ry")?,
                gx: lat["gx0"].as_i64().ok_or("gx0")?,
                gy: lat["gy0"].as_i64().ok_or("gy0")?,
                n,
                h: z.chunks_exact(4)
                    .map(|b| i32::from_le_bytes(b.try_into().unwrap()) as f64 / 1000.0)
                    .collect(),
                fix,
                rain: rw
                    .chunks_exact(4)
                    .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
                    .collect(),
                vid,
                style,
                unit: unit
                    .chunks_exact(2)
                    .map(|b| u16::from_le_bytes(b.try_into().unwrap()))
                    .collect(),
                volcs: h["prep"]["net"]["volcanoes"]
                    .as_array()
                    .ok_or("reference volcanoes")?
                    .clone(),
            });
        }
        if tiles.is_empty() {
            return Err("empty native surface reference".into());
        }
        let input_digest = format!("blake3:{}", hasher.finalize().to_hex());
        let mut kh = blake3::Hasher::new();
        kh.update(include_bytes!("t1b_ve.rs"));
        kh.update(include_bytes!("t1b_lem.rs"));
        let kernel_digest = format!("blake3:{}", kh.finalize().to_hex());
        if let Some(cache) = cache {
            let base = std::path::Path::new(cache);
            let meta = base.join("model.json");
            if meta.exists() {
                let j: serde_json::Value =
                    serde_json::from_slice(&std::fs::read(&meta).map_err(|e| e.to_string())?)
                        .map_err(|e| e.to_string())?;
                if j["input_digest"] != input_digest || j["kernel_digest"] != kernel_digest {
                    return Err("flank2 cache input/kernel mismatch; use a new cache path".into());
                }
                let mut pp = Vec::new();
                for p in j["patches"].as_array().ok_or("cache patches")? {
                    let data = std::fs::read(base.join(p["file"].as_str().ok_or("cache file")?))
                        .map_err(|e| e.to_string())?;
                    if p["digest"] != format!("blake3:{}", blake3::hash(&data).to_hex()) {
                        return Err("flank2 cache checksum".into());
                    }
                    let n = p["n"].as_u64().ok_or("cache n")? as usize;
                    if data.len() != n * n * 4 {
                        return Err("flank2 cache shape".into());
                    }
                    pp.push(Flank2Patch {
                        gx: p["gx"].as_i64().ok_or("cache gx")?,
                        gy: p["gy"].as_i64().ok_or("cache gy")?,
                        n,
                        delta: data
                            .chunks_exact(4)
                            .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
                            .collect(),
                    });
                }
                return Ok(Self {
                    rule: FLANK2_RULE.into(),
                    phase: 1,
                    reference: Some(Flank2Reference {
                        digest: input_digest,
                        windows: tiles.len() as u32,
                        input_recipes: recipes.into_iter().collect(),
                        ledgers: j["ledgers"].clone(),
                    }),
                    runtime: Flank2Runtime(std::sync::Arc::new(pp)),
                });
            }
        }
        let mut seen = vec![false; tiles.len()];
        let mut patches = Vec::new();
        let mut ledgers = Vec::new();
        for start in 0..tiles.len() {
            if seen[start] {
                continue;
            }
            let mut group = vec![start];
            seen[start] = true;
            let mut k = 0;
            while k < group.len() {
                let i = group[k];
                for j in 0..tiles.len() {
                    if !seen[j]
                        && (tiles[i].x - tiles[j].x).abs() + (tiles[i].y - tiles[j].y).abs() == 1
                    {
                        seen[j] = true;
                        group.push(j);
                    }
                }
                k += 1;
            }
            let gx = group.iter().map(|&i| tiles[i].gx).min().unwrap();
            let gy = group.iter().map(|&i| tiles[i].gy).min().unwrap();
            let nx = (group
                .iter()
                .map(|&i| tiles[i].gx + tiles[i].n as i64)
                .max()
                .unwrap()
                - gx) as usize;
            let ny = (group
                .iter()
                .map(|&i| tiles[i].gy + tiles[i].n as i64)
                .max()
                .unwrap()
                - gy) as usize;
            if nx != ny {
                return Err("reference components must be complete square window sets".into());
            }
            let n = nx;
            let nn = n * n;
            let mut z = vec![f64::NAN; nn];
            let mut fixed = vec![0u8; nn];
            let mut rain = vec![0f32; nn];
            let mut stock = vec![0f32; nn];
            let mut material = vec![0f32; nn];
            let mut quality = vec![0u8; nn];
            let mut exposure = vec![1f32; nn];
            let mut volcanoes = BTreeMap::<String, serde_json::Value>::new();
            let mut overlap_differences = 0usize;
            for &ti in &group {
                let t = &tiles[ti];
                for v in &t.volcs {
                    volcanoes
                        .entry(v["id"].as_str().ok_or("volcano id")?.into())
                        .or_insert_with(|| v.clone());
                }
                for y in 0..t.n {
                    for x in 0..t.n {
                        let s = y * t.n + x;
                        let j =
                            (t.gy - gy + y as i64) as usize * n + (t.gx - gx + x as i64) as usize;
                        if z[j].is_finite() && z[j] != t.h[s] {
                            overlap_differences += 1;
                        }
                        let q = if x >= 64 && x <= 576 && y >= 64 && y <= 576 {
                            2
                        } else {
                            1
                        };
                        if q < quality[j] {
                            continue;
                        }
                        quality[j] = q;
                        z[j] = t.h[s];
                        fixed[j] = t.fix[s];
                        rain[j] = t.rain[s];
                        if let Some(v) = t.volcs.get(t.vid[s] as usize) {
                            if let Some(u) = v["ve"]["units"].as_array().and_then(|uu| {
                                uu.iter()
                                    .find(|u| u["id"].as_u64() == Some(t.unit[s] as u64))
                            }) {
                                exposure[j] =
                                    (u["age_yr"].as_f64().ok_or("native unit age")? / 30_000.0)
                                        .clamp(0.0, 1.0) as f32;
                            }

                            if family_of(v["kind"].as_str().unwrap_or(""))
                                == Some(Family::Composite)
                            {
                                // Reduced-order mobile stock: the remaining constructed volcanic thickness.
                                // The buried basal contact is not available in the planet record.
                                let h = v["h_m"].as_f64().ok_or("volcano height")?;
                                let d = v["denudation_m"].as_f64().ok_or("denudation")?;
                                let r = ((t.gx + x as i64) as f64 * 32.0
                                    - v["vent_xy_m"][0].as_f64().ok_or("vent x")?)
                                .hypot(
                                    (t.gy + y as i64) as f64 * 32.0
                                        - v["vent_xy_m"][1].as_f64().ok_or("vent y")?,
                                );
                                stock[j] = flank2_profile(
                                    h,
                                    d,
                                    v["s_top"].as_f64().ok_or("s_top")?,
                                    v["lambda_m"].as_f64().ok_or("lambda")?,
                                    v["crater_radius_m"].as_f64().ok_or("crater")?,
                                    v["r_main_m"].as_f64().ok_or("radius")?,
                                    r,
                                ) as f32;
                            }
                        }
                        if t.style[s] == ST_FLOOR {
                            // A dry crater floor is not a lake or a drainage outlet.
                            // Protect its substrate while retaining native connected flow.
                            stock[j] = 0.0;
                        }
                        if matches!(t.style[s], ST_LAVA | ST_FRONT | ST_DOME | ST_CINDER) {
                            material[j] = 30.0;
                        }
                    }
                }
            }
            if z.iter().any(|v| !v.is_finite()) {
                return Err("native reference has holes".into());
            }
            if overlap_differences > 0 {
                return Err(format!(
                    "reference height seams disagree at {overlap_differences} posts"
                ));
            }
            // Buried contacts are an explicit stationary-emplacement inference from native VE history.
            for v in volcanoes.values() {
                if let Some(hist) = v["ve"]["flank2"]["history"].as_array() {
                    for u in hist {
                        let g = &u["geom"];
                        let Some(pts) = g["pts"].as_array() else {
                            continue;
                        };
                        let w = g["w0_m"].as_f64().unwrap_or(0.0);
                        let t = g["t_m"].as_f64().unwrap_or(0.0);
                        if w <= 0.0 || t <= 0.0 {
                            continue;
                        }
                        let remain = (1.0
                            - v["denudation_m"].as_f64().unwrap_or(0.0)
                                / v["h_m"].as_f64().unwrap_or(1.0))
                        .clamp(0.0, 1.0);
                        for seg in pts.windows(2) {
                            let (ax, ay, bx, by) = (
                                seg[0][0].as_f64().unwrap(),
                                seg[0][1].as_f64().unwrap(),
                                seg[1][0].as_f64().unwrap(),
                                seg[1][1].as_f64().unwrap(),
                            );
                            let x0 =
                                ((ax.min(bx) - w) / 32.0 - gx as f64).floor().max(0.0) as usize;
                            let x1 = ((ax.max(bx) + w) / 32.0 - gx as f64)
                                .ceil()
                                .min((n - 1) as f64) as usize;
                            let y0 =
                                ((ay.min(by) - w) / 32.0 - gy as f64).floor().max(0.0) as usize;
                            let y1 = ((ay.max(by) + w) / 32.0 - gy as f64)
                                .ceil()
                                .min((n - 1) as f64) as usize;
                            if x0 > x1 || y0 > y1 {
                                continue;
                            }
                            for y in y0..=y1 {
                                for x in x0..=x1 {
                                    let i = y * n + x;
                                    if stock[i] <= 0.0 {
                                        continue;
                                    }
                                    let p = [
                                        (gx + x as i64) as f64 * 32.0,
                                        (gy + y as i64) as f64 * 32.0,
                                    ];
                                    let q = seg_nearest(p, [ax, ay], [bx, by]);
                                    let d = (p[0] - q[0]).hypot(p[1] - q[1]);
                                    if d < w {
                                        material[i] = material[i].max(
                                            (t * remain * (1.0 - (d / w).powi(2)).powi(2)) as f32,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
            let (delta, ledger) = flank2_evolve(
                n, &z, &fixed, &rain, &stock, &material, &exposure, 32.0, 30_000.0, 120,
            )?;
            ledgers.push(ledger);
            patches.push(Flank2Patch { gx, gy, n, delta });
        }
        if let Some(cache) = cache {
            let base = std::path::Path::new(cache);
            std::fs::create_dir_all(base).map_err(|e| e.to_string())?;
            let mut desc = Vec::new();
            for (i, p) in patches.iter().enumerate() {
                let file = format!("patch_{i}.f32");
                let bytes: Vec<u8> = p.delta.iter().flat_map(|v| v.to_le_bytes()).collect();
                std::fs::write(base.join(&file), &bytes).map_err(|e| e.to_string())?;
                desc.push(serde_json::json!({"gx":p.gx,"gy":p.gy,"n":p.n,"file":file,"digest":format!("blake3:{}",blake3::hash(&bytes).to_hex())}));
            }
            std::fs::write(base.join("model.json"),serde_json::to_vec_pretty(&serde_json::json!({"input_digest":input_digest,"kernel_digest":kernel_digest,"patches":desc,"ledgers":ledgers})).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
        }
        let digest = input_digest;
        Ok(Self {
            rule: FLANK2_RULE.into(),
            phase: 1,
            reference: Some(Flank2Reference {
                digest,
                windows: tiles.len() as u32,
                input_recipes: recipes.into_iter().collect(),
                ledgers: serde_json::Value::Array(ledgers),
            }),
            runtime: Flank2Runtime(std::sync::Arc::new(patches)),
        })
    }
    pub fn at(&self, gx: i64, gy: i64, dx_m: f64) -> Option<f64> {
        let x = gx as f64 * dx_m / 32.0;
        let y = gy as f64 * dx_m / 32.0;
        for p in self.runtime.0.iter() {
            let (u, v) = (x - p.gx as f64, y - p.gy as f64);
            if u < 0.0 || v < 0.0 || u >= (p.n - 1) as f64 || v >= (p.n - 1) as f64 {
                continue;
            }
            let (ix, iy) = (u.floor() as usize, v.floor() as usize);
            let (a, b) = (u - ix as f64, v - iy as f64);
            let i = iy * p.n + ix;
            return Some(
                (1.0 - b) * ((1.0 - a) * p.delta[i] as f64 + a * p.delta[i + 1] as f64)
                    + b * ((1.0 - a) * p.delta[i + p.n] as f64 + a * p.delta[i + p.n + 1] as f64),
            );
        }
        None
    }
}

/// Finite-volume bedrock detachment, sediment transport and slope failure on an
/// actual native reference surface. No radius or azimuth enters this solve.
fn flank2_evolve(
    n: usize,
    z0: &[f64],
    fixed0: &[u8],
    rain: &[f32],
    stock0: &[f32],
    material: &[f32],
    exposure: &[f32],
    dx: f64,
    years: f64,
    steps: usize,
) -> R<(Vec<f32>, serde_json::Value)> {
    let nn = n * n;
    if n < 3
        || steps == 0
        || !dx.is_finite()
        || dx <= 0.0
        || !years.is_finite()
        || years < 0.0
        || z0.len() != nn
        || fixed0.len() != nn
        || rain.len() != nn
        || stock0.len() != nn
        || material.len() != nn
        || exposure.len() != nn
        || !exposure
            .iter()
            .all(|x| x.is_finite() && *x >= 0.0 && *x <= 1.0)
        || !z0.iter().all(|x| x.is_finite())
        || !rain
            .iter()
            .chain(stock0)
            .chain(material)
            .all(|x| x.is_finite() && *x >= 0.0)
    {
        return Err("invalid native flank process inputs".into());
    }
    let area_cell = dx * dx;
    let dt = years / steps as f64;
    let mut fixed = fixed0.to_vec();
    for i in 0..nn {
        if i % n == 0 || i % n == n - 1 || i / n == 0 || i / n == n - 1 {
            fixed[i] = crate::t1b::RING;
        }
    }
    let mut z = z0.to_vec();
    let mut stock: Vec<f64> = stock0.iter().map(|&x| x as f64).collect();
    let mut eroded = 0.;
    let mut deposited = 0.;
    let mut exported = 0.;
    let mut slide_volume = 0.;
    let mut water_error = 0f64;
    let mut max_area = 0f64;
    let mut max_width = 0f64;
    let mut confluences = 0usize;
    let source_area = (0..nn)
        .filter(|&i| fixed[i] == FREE)
        .map(|i| area_cell * rain[i] as f64)
        .sum::<f64>();
    for _step in 0..steps {
        if _step % 10 == 0 {
            eprintln!("ve_flank2 native solve {n}x{n} step {_step}/{steps}");
        }
        let route = crate::t1b_lem::route(n, &z, &fixed, dx);
        let mut weights = vec![[0f64; 8]; nn];
        let mut pos = vec![0u32; nn];
        for (k, &i) in route.order.iter().enumerate() {
            pos[i as usize] = k as u32;
        }
        for i in 0..nn {
            if fixed[i] != FREE {
                continue;
            }
            let x = (i % n) as i64;
            let y = (i / n) as i64;
            let mut sum = 0.;
            for (k, (ox, oy)) in crate::t1b_lem::D8N.iter().copied().enumerate() {
                let j = (y + oy) as usize * n + (x + ox) as usize;
                if pos[j] >= pos[i] || route.fill[j] >= route.fill[i] {
                    continue;
                }
                let len = dx * (ox as f64).hypot(oy as f64);
                let s = (route.fill[i] - route.fill[j]) / len;
                let w = s.powf(1.1) / len;
                weights[i][k] = w;
                sum += w;
            }
            if sum > 0. {
                let sum = weights[i].iter().map(|&w| w as f64).sum::<f64>();
                if !sum.is_finite() || sum <= 0.0 {
                    return Err(format!("invalid routing weights at step {_step} post {i}, z={} fill={} weights={:?}",z[i],route.fill[i],weights[i]));
                }
                for w in &mut weights[i] {
                    *w /= sum;
                }
            } else {
                let j = route.receiver[i] as usize;
                let ox = (j % n) as i64 - x;
                let oy = (j / n) as i64 - y;
                let k = crate::t1b_lem::D8N
                    .iter()
                    .position(|&a| a == (ox, oy))
                    .ok_or("invalid native flood receiver")?;
                weights[i][k] = 1.;
            }
        }
        let mut a: Vec<f64> = (0..nn)
            .map(|i| {
                if fixed[i] == FREE {
                    area_cell * rain[i] as f64
                } else {
                    0.
                }
            })
            .collect();
        let mut donors = vec![0u8; nn];
        for &i in route.order.iter().rev() {
            let i = i as usize;
            if fixed[i] != FREE {
                continue;
            }
            let sum = weights[i].iter().map(|&w| w as f64).sum::<f64>();
            if !sum.is_finite() || sum <= 0.0 {
                return Err(format!(
                    "invalid routing weights at step {_step} post {i}, z={} fill={} weights={:?}",
                    z[i], route.fill[i], weights[i]
                ));
            }
            for (k, (ox, oy)) in crate::t1b_lem::D8N.iter().copied().enumerate() {
                let w = weights[i][k] as f64 / sum;
                if w <= 0. {
                    continue;
                }
                let j = ((i / n) as i64 + oy) as usize * n + ((i % n) as i64 + ox) as usize;
                a[j] += w * a[i];
                if a[i] > 4096. && w > 0.1 {
                    donors[j] = donors[j].saturating_add(1);
                }
            }
        }
        let sinks = (0..nn)
            .filter(|&i| fixed[i] != FREE)
            .map(|i| a[i])
            .sum::<f64>();
        if !sinks.is_finite() || !a.iter().all(|x| x.is_finite()) {
            return Err(format!("nonfinite water at step {_step}"));
        }
        water_error = water_error.max((sinks - source_area).abs() / source_area.max(1.));
        confluences = donors.iter().filter(|&&x| x > 1).count();
        max_area = a.iter().copied().fold(0., f64::max);
        let mut sediment = vec![0f64; nn];
        // Donor-first sediment routing permits deposition when carrying capacity falls.
        for &i in route.order.iter().rev() {
            let i = i as usize;
            if !sediment[i].is_finite() {
                return Err(format!(
                    "nonfinite sediment at step {_step} post {i}, area={} weights={:?}",
                    a[i], weights[i]
                ));
            }
            if fixed[i] != FREE {
                exported += sediment[i];
                continue;
            }
            let sum = weights[i].iter().map(|&w| w as f64).sum::<f64>();
            if !sum.is_finite() || sum <= 0.0 {
                return Err(format!(
                    "invalid routing weights at step {_step} post {i}, z={} fill={} weights={:?}",
                    z[i], route.fill[i], weights[i]
                ));
            }
            let mut slope = 0.;
            for (k, (ox, oy)) in crate::t1b_lem::D8N.iter().copied().enumerate() {
                let w = weights[i][k] as f64 / sum;
                if w <= 0. {
                    continue;
                }
                let j = ((i / n) as i64 + oy) as usize * n + ((i % n) as i64 + ox) as usize;
                slope += w * (route.fill[i] - route.fill[j]).max(0.)
                    / (dx * (ox as f64).hypot(oy as f64));
            }
            let q = a[i] * 5e-6;
            let width = 3.5 * q.sqrt();
            max_width = max_width.max(width);
            let capacity = 0.02 * a[i] * slope * slope * dt;
            let available = (capacity - sediment[i]).max(0.);
            if z[i] >= route.fill[i] - 1e-9 && stock[i] > 0. {
                // A/W*S is unit stream power; W proportional sqrt(Q) gives the familiar sqrt(A)*S law.
                let power = if width > 0. { a[i] * slope / width } else { 0. };
                let width_coefficient = 3.5 * (5e-6_f64).sqrt();
                let resistance = 1.0 / (1.0 + material[i] as f64 / 20.0);
                let potential = 32e-6
                    * width_coefficient
                    * (power - 12.0 / width_coefficient).max(0.)
                    * dt
                    * resistance
                    * exposure[i] as f64;
                let depth = potential.min(stock[i]).min(available / area_cell);
                z[i] -= depth;
                stock[i] -= depth;
                let dv = depth * area_cell;
                sediment[i] += dv;
                eroded += dv;
            }
            if sediment[i] > capacity {
                // Finite adaptation to transport capacity over a200m settling length [H].
                // Unsettled load remains in transport, never disappears.
                let dh = (sediment[i] - capacity) / area_cell * (1.0 - (-dx / 200.0).exp());
                z[i] += dh;
                stock[i] += dh;
                let dv = dh * area_cell;
                sediment[i] -= dv;
                deposited += dv;
            }
            for (k, (ox, oy)) in crate::t1b_lem::D8N.iter().copied().enumerate() {
                let w = weights[i][k] as f64 / sum;
                if w <= 0. {
                    continue;
                }
                let j = ((i / n) as i64 + oy) as usize * n + ((i % n) as i64 + ox) as usize;
                sediment[j] += w * sediment[i];
            }
        }
        // Conservative granular bank/fan relaxation. Repose varies with native lava contacts.
        // Only mobile volcanic material moves; no rock is created or clipped.
        for _ in 0..8 {
            let mut moved = 0f64;
            for parity in 0..2 {
                for i in 0..nn {
                    if fixed[i] != FREE || (i % n + i / n) % 2 != parity {
                        continue;
                    }
                    for (ox, oy) in [(1i64, 0i64), (0, 1), (1, 1), (-1, 1)] {
                        let j = ((i / n) as i64 + oy) as usize * n + ((i % n) as i64 + ox) as usize;
                        if fixed[j] != FREE {
                            continue;
                        }
                        let (hi, lo) = if z[i] > z[j] { (i, j) } else { (j, i) };
                        let len = dx * (ox as f64).hypot(oy as f64);
                        let sc = 0.7 + 0.5 * (material[hi] as f64 / (material[hi] as f64 + 20.));
                        let dh = (0.5 * (z[hi] - z[lo] - sc * len)).max(0.).min(stock[hi]);
                        if dh > 0. {
                            z[hi] -= dh;
                            z[lo] += dh;
                            stock[hi] -= dh;
                            stock[lo] += dh;
                            moved += dh * area_cell;
                        }
                    }
                }
            }
            slide_volume += moved;
            if moved < 1e-6 {
                break;
            }
        }
    }
    let change = z
        .iter()
        .zip(z0)
        .map(|(a, b)| (a - b) * area_cell)
        .sum::<f64>();
    let residual = (change + exported).abs() / exported.max(1.);
    if !z.iter().all(|x| x.is_finite())
        || !exported.is_finite()
        || !residual.is_finite()
        || !water_error.is_finite()
        || water_error > 1e-9
        || residual > 1e-8
    {
        return Err(format!(
            "flank2 conservation/nonfinite failure water={water_error} sediment={residual}"
        ));
    }
    let mut out = Vec::with_capacity(nn);
    let (mut maxcut, mut maxdep) = (0f64, 0f64);
    for i in 0..nn {
        let d = z[i] - z0[i];
        maxcut = maxcut.max(-d);
        maxdep = maxdep.max(d);
        out.push(d as f32);
    }
    Ok((
        out,
        serde_json::json!({"source":"hashed native construction windows","n":n,"spacing_m":dx,"duration_yr":years,"steps":steps,"rain_source_area_m2":source_area,"water_relative_residual":water_error,"eroded_m3":eroded,"deposited_m3":deposited,"bank_redistributed_m3":slide_volume,"exported_m3":exported,"net_volume_change_m3":change,"volume_relative_residual":residual,"max_cut_m":maxcut,"max_deposit_m":maxdep,"max_contributing_area_m2":max_area,"max_bankfull_width_m":max_width,"confluence_posts":confluences}),
    ))
}

/// Constant-rate denudation is an approximation to integrated erosion, not an edifice erosion history.
pub const EREC_AGE_RULE: &str = "edifice_record_erosion_v1 (CX-EDIFICE-EREC, ve_erec=1, layer AND windows): the existing CO body D = erosion_mult x E x age and VE unit age laws read E_rec instead of fitted regional E. E_rec is the bound T0 frame-200 recent bedrock erosion rate (398-400 My for C6/C7), with the verbatim record loader's lifetime fallback and 5 m/My floor. Hat-blend record rates in degrees at the record point for pre-placement residual clearance, then at the actual cleared vent for body and unit degradation. No solved-table mask, coast taper or fitted E enters this read. The recent 2 My rate is extrapolated over the record dormancy/unit age [H: constant-rate approximation, not measured edifice-specific erosion history]. Geometry, ages, thresholds, floor, erosion multiplier and VE2 young-unit weathering law are unchanged. No solve variable is changed directly; changed edifice relief still participates in the existing regional solve. Draws retain the existing ls.tb paired keying.";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct VeErecParamsV1 {
    pub rule: String,
    pub record: crate::t0_erate::ErecRecordV1,
}

/// Record erosion at an edifice point, independent of the fitted solve, including its zero-bound cells.
/// With the option absent, call the original interpolation exactly.
pub fn erosion_at(ve: Option<&VeParamsV1>, t0: &T0, chart: &Chart, table: &BTreeMap<u32, f64>, x: f64, y: f64) -> R<f64> {
    let Some(p) = ve.and_then(|v| v.erec.as_ref()) else { return crate::t1b_coast::table_at(t0, chart, table, x, y); };
    let (la, lo) = chart.inverse_f64(x * 1000.0, y * 1000.0).map_err(|e| e.to_string())?;
    let (cells, weights) = crate::t1b_net::hat4(t0, la.to_degrees(), lo.to_degrees());
    let mut e = 0.0;
    for k in 0..4 {
        if weights[k] > 0.0 {
            let rate = *p.record.rec.e.get(cells[k] as usize).ok_or("ve_erec: record table does not cover vent")?;
            e += weights[k] * rate;
        }
    }
    Ok(e)
}

/// Lane VH-I: the climate wind field (in memory only, never serialized; the recipe carries the artifact digest).
#[derive(Clone, Default)]
pub struct VeWindField { pub u: std::sync::Arc<Vec<f32>>, pub v: std::sync::Arc<Vec<f32>> }
impl std::fmt::Debug for VeWindField { fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "VeWindField({} cells)", self.u.len()) } }
impl PartialEq for VeWindField { fn eq(&self, _: &Self) -> bool { true } }
impl Eq for VeWindField {}

pub const WIND_RULE: &str = "ve_wind_t0_v1 (lane VH-I, interim of the eruption-history design I2; [H] where marked): the dispersal axis of every composite explosive event (VE's apron unit dir_rad: the fall isopach ellipse and the density-current sector) is the T0 surface wind azimuth plus a keyed deviation, theta = atan2(d_y, d_x) + 2 atan(((1 - rho) / (1 + rho)) tan(pi (u - 1/2))) (inverse CDF of a wrapped Cauchy, a von Mises-like closed form), u the SAME keyed uniform that drew the uniform axis before (no draw count changes). Wind: (wind_u_m_s, wind_v_m_s) of the climate= artifact (seed-climate-v1, digest below) at the volcano record's T0 cell (a planet property: no chart or window enters it), as d = u e + v n with e, n the chart-frame unit east and north vectors at the record point (finite differences of the chart projection); the axis points downwind (where the wind blows to). rho = 0.698 is the mean resultant length of a von Mises with kappa 2 [H]; below calm (0.5 m/s) [H] the uniform axis is kept. The surface wind stands in for the plume-height wind (fair in mid-latitude westerlies, recalled). The record's wind and axis are published in prep.net.volcanoes[k].ve.wind";

/// Lane VH-I ve_wind_t0_v1 parameters.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct VeWindParamsV1 {
    pub version: u32,
    pub rho_permille: i64,
    pub calm_mm_per_s: i64,
    /// Digest of the climate artifact the wind was read from (must be the climate= artifact).
    pub climate_digest: String,
    pub rule: String,
    #[serde(skip)]
    pub field: VeWindField,
}

impl VeWindParamsV1 {
    /// Read wind_u_m_s / wind_v_m_s of a verified pf.climate artifact directory.
    pub fn load(dir: &std::path::Path) -> R<VeWindParamsV1> {
        let digest = pf_contracts::artifact::verify_artifact::<serde_json::Value>(dir).map_err(|e| format!("verify climate {}: {e}", dir.display()))?;
        let rd = |n: &str| -> R<Vec<f32>> { let v = pf_contracts::artifact::read_field_f32(dir, n).map_err(|e| format!("climate {n}: {e}"))?; if v.iter().any(|x| !x.is_finite()) { return Err(format!("climate {n}: nonfinite value")); } Ok(v) };
        let (u, v) = (rd("wind_u_m_s")?, rd("wind_v_m_s")?);
        if u.len() != v.len() { return Err("climate wind fields differ in length".into()); }
        Ok(VeWindParamsV1 { version: 1, rho_permille: 698, calm_mm_per_s: 500, climate_digest: digest, rule: WIND_RULE.into(), field: VeWindField { u: std::sync::Arc::new(u), v: std::sync::Arc::new(v) } })
    }
}

/// Lane VH-I: the T0 wind at a volcano record (published in the header).
#[derive(Clone, Debug, Serialize)]
pub struct VeWindAt {
    pub cell: u32,
    pub u_m_s: f64,
    pub v_m_s: f64,
    /// Chart-frame downwind azimuth (rad, atan2(y, x)); None when calm (uniform axis kept).
    pub chart_azimuth_rad: Option<f64>,
    pub rho: f64,
}

/// Lane VH-I: wrapped-Cauchy deviation (rad) from a uniform u in [0, 1) for mean resultant length rho.
#[inline]
pub fn wind_deviation(rho: f64, u: f64) -> f64 { 2.0 * (((1.0 - rho) / (1.0 + rho)) * (std::f64::consts::PI * (u - 0.5)).tan()).atan() }

/// Lane VH-I: the T0 wind at the record of `v` (its record cell) in the chart frame.
pub fn wind_at(w: &VeWindParamsV1, t0: &T0, chart: &Chart, v: &VolcanoSite) -> R<VeWindAt> {
    let n = t0.nlat * t0.nlon;
    if w.field.u.len() != n || w.field.v.len() != n { return Err("ve_wind: the wind field is not attached or not on the T0 grid (give ve_wind=1 with climate=<dir>)".into()); }
    match t0.climate.as_ref() { Some(c) if c.digest == w.climate_digest => {} _ => return Err("ve_wind: the wind must come from the climate= artifact attached to T0".into()) }
    let c = v.record_cell as usize;
    if c >= n { return Err("ve_wind: record cell outside the T0 grid".into()); }
    let (u, vv) = (w.field.u[c] as f64, w.field.v[c] as f64);
    let (la, lo) = chart.inverse_f64(v.record_xy_m[0] * 1000.0, v.record_xy_m[1] * 1000.0).map_err(|e| e.to_string())?;
    let d = 1e-5;
    let p0 = chart.forward_f64(la, lo).map_err(|e| e.to_string())?;
    let (lo1, se) = if lo + d <= std::f64::consts::PI { (lo + d, 1.0) } else { (lo - d, -1.0) };
    let (la1, sn) = if la + d <= std::f64::consts::FRAC_PI_2 { (la + d, 1.0) } else { (la - d, -1.0) };
    let pe = chart.forward_f64(la, lo1).map_err(|e| e.to_string())?;
    let pn = chart.forward_f64(la1, lo).map_err(|e| e.to_string())?;
    let unit = |a: [f64; 2], s: f64| { let l = a[0].hypot(a[1]).max(1e-300); [s * a[0] / l, s * a[1] / l] };
    let e = unit([pe[0] - p0[0], pe[1] - p0[1]], se);
    let no = unit([pn[0] - p0[0], pn[1] - p0[1]], sn);
    let (dx, dy) = (u * e[0] + vv * no[0], u * e[1] + vv * no[1]);
    let calm = u.hypot(vv) * 1000.0 < w.calm_mm_per_s as f64;
    Ok(VeWindAt { cell: v.record_cell, u_m_s: u, v_m_s: vv, chart_azimuth_rad: if calm { None } else { Some(dy.atan2(dx)) }, rho: w.rho_permille as f64 / 1000.0 })
}

/// Lane VE5 (2026-10-15; `unpublished-design-record`): exposure-age lobes, summit taper of flows.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Ve5ParamsV1 {
    pub exposure: u32,
    pub summit: u32,
    pub exposure_rule: String,
    pub summit_rule: String,
}

impl Ve5ParamsV1 {
    pub fn v1() -> Self {
        Ve5ParamsV1 {
            exposure: 1, summit: 1,
            exposure_rule: "ve5_exposure_age_v1 (replaces ve4's log-uniform lobe age): the visible surface lobes are the youngest unburied ones; older lobes are buried by younger units, not eroded in place, so a lobe's age is its surface exposure age. Resurfacing from the edifice's own construction: the record height H was built within the edifice lifetime L_e with a mean thickness of H / 3 over the body (a cone; also declared for the shield), so the mean accretion rate is R = H / (3 L_e) and a surface post is resurfaced by a lobe of thickness t_med every tau_e = t_med / R on average (Poisson); a visible lobe's exposure age = record age + T_h + tau_e x (-ln(1 - u)) (just older than the mapped window; u = ve4's keyed age uniform), then the fresh-lava erosion law (weathering until the lag, then the cell E)".into(),
            summit_rule: "ve5_summit_flows_v1: lava-path deposits are multiplied by smoothstep((r - r_t) / (0.5 r_t)), r the distance from the vent, r_t = max(the body's crater radius, the caldera unit radius): flows that start at the crater or caldera rim do not stand above it (full thickness from 1.5 r_t); steering keeps the untapered drape (declared)".into(),
        }
    }
    pub fn off() -> Self { let mut p = Self::v1(); p.exposure = 0; p.summit = 0; p }
}

/// Lane VE4 (2026-10-15; `unpublished-design-record`): lobes by slope and age with a summit taper,
/// lahar dammed ponds, the front_mm snout field. Each item has its own switch (0 off).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Ve4ParamsV1 {
    pub lobes: u32,
    pub lahar_pond: u32,
    pub front: u32,
    pub front_cap_mm: i64,
    pub lobes_rule: String,
    pub lahar_pond_rule: String,
    pub front_rule: String,
}

impl Ve4ParamsV1 {
    pub fn v1() -> Self {
        Ve4ParamsV1 {
            lobes: 1, lahar_pond: 1, front: 1, front_cap_mm: 64_000,
            lobes_rule: "ve4_lobes_v1 (replaces ve3's uniform lobes): each roughness cell is occupied with keyed probability 0.8 (category 21, 8 draws per cell); a lobe's age is keyed log-uniform between the mappable window T_h and the edifice lifetime L_e (older unmapped flows) and it is eroded by the fresh-lava law (weathering until the lag, then the cell E); size and amplitude scale with the local body slope s by Jeffreys' viscous sheet thickness h ~ (q mu / (g sin theta))^(1/3) (recalled: thicker, wider lobes on gentler lower slopes): f = clamp((s_med / s)^(1/3), 0.6, 1.6), s_med = the median body slope over the occupied cells; amplitude max(0, t_med U(0.3, 1.0) f - erosion), half-width w_med U(0.6, 1.2) f; summit taper x smoothstep((r - r_c) / (2 r_c)), r_c = the body's crater (caldera) radius, so lobes do not stack on the crater rim; search 7 x 7 cells".into(),
            lahar_pond_rule: "ve4_lahar_pond_v1: a steered lahar whose next step rises fills the closed low with its remaining volume (remaining length x 2 w0 x t) by the ve3 pond rule (flood from the bottom, to the spill level or the volume level, overflow continues); published flat, style 5 (lahar), the lahar's unit id; above sea level only".into(),
            front_rule: "ve4_front_field_v1: window field front_mm (i32 mm): on posts within front_cap of the snout outline of a lava flow (the region past its path end point), the signed horizontal distance to that outline, > 0 inside the snout; the outline is where the drape thickness reaches 0, radius w_end (1 + 0.25 lobes(angle)); the youngest flow wins; i32::MIN elsewhere. Consumers: the zero contour of the bilinear field is the front at any spacing; thickness behind it t smoothstep(front_mm / m), m = max(t / 0.7, 1.5 x 32 m) with ve3 aa (t from the unit record, after erosion)".into(),
        }
    }
    pub fn off() -> Self { let mut p = Self::v1(); p.lobes = 0; p.lahar_pond = 0; p.front = 0; p }
}

/// Lane VE4 front_mm sentinel.
pub const FRONT_NONE: i32 = i32::MIN;
pub const FRONT_NOTE: &str = "lane VE4 (T1b, opt-in ve4 front): signed horizontal distance in mm from the post to the snout outline of the youngest lava flow whose snout (the region past its path end point) lies within the recipe's front_cap of the post; > 0 inside the snout, < 0 outside; -2147483648 (i32::MIN) elsewhere. The outline is where the flow's drape thickness reaches 0 (radius w_end (1 + 0.25 lobes), unit record in prep.net.volcanoes[..].ve.units). Consumer use: interpolate bilinearly; the zero isoline is the front at any spacing; the thickness behind it is t smoothstep(front_mm / m), m = max(t / 0.7, 48 m) with ve3 aa";

/// Lane VE3 (2026-10-15; `unpublished-design-record`): sub-resolution lobe roughness, ponded lava,
/// anti-aliased flow margins. Each item has its own switch (0 off).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Ve3ParamsV1 {
    pub rough: u32,
    pub pond: u32,
    pub aa: u32,
    pub pond_radius_mm: i64,
    pub ponds_max: u32,
    pub rough_rule: String,
    pub pond_rule: String,
    pub aa_rule: String,
}

impl Ve3ParamsV1 {
    pub fn v1() -> Self {
        Ve3ParamsV1 {
            rough: 1, pond: 1, aa: 1, pond_radius_mm: 1_500_000, ponds_max: 3,
            rough_rule: "ve3_lobe_roughness_v1 (declared sub-resolution closure: the older, unmapped flows that make up the edifice body): lobes on a chart grid of cell L = 2 w_med (w_med = median half-width of the volcano's own lava-flow units), one per cell inside the body at a keyed position (purpose t1b.volcano.ve, category 20, counter = cell index), half-width w_med U(0.6, 1.2), half-length 3 x that, long axis along the body fall line at the lobe centre, amplitude U(0.3, 1.0) x t_med (median thickness of its own flow units after erosion), profile a sqrt(1 - q^2) with q the normalised elliptical distance, overlaps take the max; added to the steering surface and to the published body (style body, unit 0)".into(),
            pond_rule: "ve3_ponded_lava_v1: lava flows only (not lahars), above sea level: where a steered flow's next step would rise (checked at four quarter sub-steps so a 128 m step cannot hop an older levee), the closed low of the steering surface (body, roughness, older units, this flow's earlier ponds) is filled: priority flood on the global 32 m post grid from the bottom of the low (8-neighbour descent from the stop post), within pond_radius; the spill level is the highest post taken before the first neighbour lower than the current level (the spill post); the fill level is the spill level or, if the low holds more than the flow's remaining volume (remaining length x 2 w0 x t), the level holding exactly that volume (bisection); the flow continues from the spill post with the length reduced by the pond volume / (2 w0 t), at most ponds_max ponds per flow, ending when the volume is used. Published: z = max(z, min(level, z + 3 t)) on the pond posts, style 10 (flat floor: ponded lava), the flow's unit id. Younger flows do not see older ponds (declared)".into(),
            aa_rule: "ve3_margin_aa_v1: flow and steering drape margins m = max(t / 0.7, 1.5 x post spacing): the cell average of a ~35 deg front at 32 m posts (a point-sampled 1.3-post margin aliases into a serrated edge, measured); thin flows get a 48 m edge at 32 m posts (resolution limit)".into(),
        }
    }
    pub fn off() -> Self { let mut p = Self::v1(); p.rough = 0; p.pond = 0; p.aa = 0; p }
}

/// Lane VE3 pond tracing settings (half-width, thickness, radius, ponds per flow).
#[derive(Clone, Copy, Debug)]
pub struct PondCfg { pub w: f64, pub t: f64, pub radius: f64, pub max: u32, pub edge_spill: bool }

/// Lane VE3: one ponded-lava fill (global 32 m post rows: [gy, gx0, gx1]).
#[derive(Clone, Debug, Serialize)]
pub struct Pond { pub level_m: f64, pub volume_m3: f64, pub cells: u32, pub rows: Vec<[i64; 3]> }

/// Lane VE3: fill the closed low around `p` on the surface `g` (global 32 m grid, within `radius`), at most `vol` m3.
/// Returns the pond and the spill point (None: no outlet inside the radius).
pub fn fill_pond(g: &dyn Fn([f64; 2]) -> f64, p: [f64; 2], radius: f64, vol: f64) -> (Pond, Option<[f64; 2]>) { fill_pond_e(g, p, radius, vol, false) }

/// As `fill_pond`; `edge_spill` (lane VE5 bug fix): a low that is not closed inside the radius is not a closed low, so reaching
/// the search box edge (in the descent walk or the flood) is a spill there (ve3 / ve4 filled against the box edge: rectangles).
pub fn fill_pond_e(g: &dyn Fn([f64; 2]) -> f64, p: [f64; 2], radius: f64, vol: f64, edge_spill: bool) -> (Pond, Option<[f64; 2]>) {
    use std::collections::BinaryHeap;
    use std::cmp::Reverse;
    let d = 32.0;
    let cs = ((p[0] / d).round() as i64, (p[1] / d).round() as i64);
    let rc = (radius / d).ceil() as i64;
    let key = |z: f64| -> i64 { (z * 1000.0).round() as i64 };
    let mut zc: BTreeMap<(i64, i64), f64> = BTreeMap::new();
    let val = |c: (i64, i64), zc: &mut BTreeMap<(i64, i64), f64>| -> f64 { *zc.entry(c).or_insert_with(|| g([c.0 as f64 * d, c.1 as f64 * d])) };
    // Walk down to the bottom of the low first (8-neighbour steepest descent on the grid, within the radius).
    let mut c0 = cs;
    for _ in 0..(4 * rc) {
        let z = val(c0, &mut zc);
        let mut best = (z, c0);
        for (ox, oy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)] {
            let n = (c0.0 + ox, c0.1 + oy);
            if (n.0 - cs.0).abs() > rc || (n.1 - cs.1).abs() > rc { continue; }
            let zn = val(n, &mut zc);
            if zn < best.0 { best = (zn, n); }
        }
        if best.1 == c0 { break; }
        c0 = best.1;
        if edge_spill && ((c0.0 - cs.0).abs() == rc || (c0.1 - cs.1).abs() == rc) {
            return (Pond { level_m: val(c0, &mut zc), volume_m3: 0.0, cells: 0, rows: Vec::new() }, Some([c0.0 as f64 * d, c0.1 as f64 * d]));
        }
    }
    let mut taken: Vec<(i64, i64)> = Vec::new();
    let mut seen: std::collections::BTreeSet<(i64, i64)> = std::collections::BTreeSet::new();
    let mut heap = BinaryHeap::new();
    let z0 = val(c0, &mut zc);
    heap.push(Reverse((key(z0), c0.1, c0.0)));
    seen.insert(c0);
    let mut level = z0;
    let mut spill: Option<[f64; 2]> = None;
    while let Some(Reverse((_, cy, cx))) = heap.pop() {
        let c = (cx, cy);
        let z = val(c, &mut zc);
        if !z.is_finite() { break; }
        if z < level - 1e-6 { spill = Some([cx as f64 * d, cy as f64 * d]); break; }
        if edge_spill && ((cx - cs.0).abs() == rc || (cy - cs.1).abs() == rc) { level = level.max(z); spill = Some([cx as f64 * d, cy as f64 * d]); break; }
        level = level.max(z);
        taken.push(c);
        for (ox, oy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)] {
            let n = (cx + ox, cy + oy);
            if (n.0 - cs.0).abs() > rc || (n.1 - cs.1).abs() > rc || seen.contains(&n) { continue; }
            seen.insert(n);
            let zn = val(n, &mut zc);
            if zn.is_finite() { heap.push(Reverse((key(zn), n.1, n.0))); }
        }
        if taken.len() > 40_000 { break; }
    }
    let vol_at = |l: f64| -> f64 { taken.iter().map(|c| (l - zc[c]).max(0.0) * d * d).sum() };
    let mut lv = level;
    if vol_at(level) > vol {
        let (mut a, mut b) = (taken.iter().map(|c| zc[c]).fold(f64::INFINITY, f64::min), level);
        for _ in 0..60 { let m = 0.5 * (a + b); if vol_at(m) > vol { b = m; } else { a = m; } }
        lv = a;
    }
    let mut cells: Vec<(i64, i64)> = taken.iter().copied().filter(|c| zc[c] < lv).collect();
    cells.sort_by_key(|c| (c.1, c.0));
    let mut rows: Vec<[i64; 3]> = Vec::new();
    for c in &cells { match rows.last_mut() { Some(r) if r[0] == c.1 && r[2] + 1 == c.0 => r[2] = c.0, _ => rows.push([c.1, c.0, c.0]) } }
    let v = vol_at(lv);
    (Pond { level_m: lv, volume_m3: v, cells: cells.len() as u32, rows }, if lv >= level - 1e-6 { spill } else { None })
}

/// Lane VE3 roughness (see rough_rule): lobes and their cell index (never serialized: rebuilt from the keyed draws).
#[derive(Clone, Debug, Serialize)]
pub struct RoughSpec {
    pub w_med_m: f64,
    pub t_med_m: f64,
    pub cell_m: f64,
    pub lobes_n: u32,
    /// (x, y, amplitude, half-width, half-length, axis angle)
    #[serde(skip)]
    pub lobes: Vec<[f64; 6]>,
    #[serde(skip)]
    pub grid: BTreeMap<(i64, i64), u32>,
    /// Cells searched around a point (2: ve3; 3: ve4 lobes, larger on gentle slopes).
    #[serde(skip)]
    pub reach: i64,
}

/// Lane VE3: lobe roughness height at chart (x, y).
pub fn rough_at(r: &RoughSpec, x: f64, y: f64) -> f64 {
    let (cx, cy) = ((x / r.cell_m).floor() as i64, (y / r.cell_m).floor() as i64);
    let mut h: f64 = 0.0;
    for gy in cy - r.reach..=cy + r.reach { for gx in cx - r.reach..=cx + r.reach {
        let Some(&k) = r.grid.get(&(gx, gy)) else { continue };
        let l = r.lobes[k as usize];
        let (dx, dy) = (x - l[0], y - l[1]);
        let (c, s) = (l[5].cos(), l[5].sin());
        let (al, la) = (dx * c + dy * s, -dx * s + dy * c);
        let q2 = (al / l[4]).powi(2) + (la / l[3]).powi(2);
        if q2 < 1.0 { h = h.max(l[2] * (1.0 - q2).sqrt()); }
    } }
    h
}

/// Lane VE2 (2026-10-15; `unpublished-design-record`): opt-in additions on ve=1. Every item has its own
/// switch (0 off).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Ve2ParamsV1 {
    pub steer: u32,
    pub erode: u32,
    pub dome: u32,
    pub taper: u32,
    /// Fresh-lava infiltration lag (yr): basaltic units, andesite / dacite composite flows and domes; weathering rate (mm/My).
    pub lag_basalt_yr: i64,
    pub lag_andesite_yr: i64,
    pub weather_mm_per_my: i64,
    pub steer_rule: String,
    pub erode_rule: String,
    pub dome_rule: String,
    pub taper_rule: String,
}

impl Ve2ParamsV1 {
    pub fn v1() -> Self {
        Ve2ParamsV1 {
            steer: 1, erode: 1, dome: 1, taper: 1, lag_basalt_yr: 50_000, lag_andesite_yr: 10_000, weather_mm_per_my: 5_000,
            steer_rule: "REVISION 1: the composite gully term applies at every unit age (the CF radial gully set stands for the cone's drainage at any time; without it radial flows on the axisymmetric cone ran straight). ve2_steer_v1: paths are traced in age order (oldest first) on G = body + (composite) - CF gully cut + the drape heights of the older units already built (lava / lahar paths: their deposit thickness with the deposit rule's profile, cinder cones, domes, fissures); steepest descent of G by 64 m differences (off the body: T0 bilinear gradient + the older units' gradient); the keyed sinusoidal wander is off; a step that would rise on G (a closed low between older levees, a cone or a dome) ends the path (ponded front). Older units with a lateral gradient above the local G slope (a 15-50 m flow with ~35 deg margins on a 5-30 deg flank) deflect a younger flow along them".into(),
            erode_rule: "ve2_fresh_lava_v1: coherent lava units (lava, cone, shield, rift and plateau flows, domes, cinder cones, fissure ramparts) lose E_w min(a, tau) + E_cell max(0, a - tau), E_w = weather (5 m/My), tau = lag_basalt (50 ky) for basaltic units (shield, rift field, plateau, cinder cones and their flows) and lag_andesite (10 ky) for composite lava flows and domes; recalled, unverified: young Hawaiian flows keep primary surface morphology and show no channel network for 10^4-10^5 yr because rain infiltrates (Jefferson, Grant and Lewis 2010; Lohse and Dietrich 2005; Ferrier et al. 2013 Nature), cinder cones degrade over 10^5-10^6 yr (Wood 1980; Hooper and Sheridan 1998); unconsolidated units (aprons, lahars) and depression depths (craters, calderas, pits, maars: infill) keep the cell E".into(),
            dome_rule: "ve2_dome_cap_v1: a composite dome lands in the summit crater with declared share 0.5 (centre r < 0.5 r_c, one extra keyed draw), else on the upper flank (0.1-0.35 R_f); a summit dome's radius is at most the crater radius; every dome's top is at most the crater rim (the body surface at r_c, record-age lowered) when its centre is within 2 r_c of the vent, so stacked summit domes fill the crater instead of standing above the record summit".into(),
            taper_rule: "ve2_taper_v1: lahar ribbons take a margin of 0.5 w (thin sides) instead of t / 0.7; apron density-current onset ramps over 3 km beyond r_a (was 1.5 km), outer taper over 0.5 of the reach (was 0.3)".into(),
        }
    }
    pub fn off() -> Self { let mut p = Self::v1(); p.steer = 0; p.erode = 0; p.dome = 0; p.taper = 0; p }
}

impl VeParamsV1 {
    pub fn v1() -> Self {
        VeParamsV1 {
            on: 1, t_map_mm: 2_000, path_step_mm: 128_000, capture_mm: 600_000, t1a_stop_mm: 300_000, gully_age_yr: 3_000, max_units: 96, ve2: None, ve3: None, ve4: None, ve5: None, wind: None, erec: None, flank: None, flank2: None,
            style_rule: "ve_style_v1 (declared closure; the record has no composition): island_arc, continental_arc -> composite stratovolcano (alternating effusive andesite/dacite flows and explosive Vulcanian/Plinian VEI 3-6 eruptions, pyroclastic aprons, lahars, domes, summit crater, sector-collapse amphitheatre, caldera for a large old system, parasitic cinder cones, maars where wet; e.g. Fuji, Mayon, Rainier, St Helens 1980, Mazama/Crater Lake, Taranaki ring plain); hotspot -> basaltic shield (long thin flows, 2-3 rift zones with fissures, spatter/cinder cones and pit craters, summit caldera; Mauna Loa, Kilauea); rift -> rift fissure field (low elongate fissure-fed shield on a keyed rift axis: the record carries no extension direction; en-echelon fissure rows with ramparts, cinder cones and flow fields; maars where wet; Laki, Afar, harrats); large_igneous_province -> flood plateau (stacked sheets with step fronts; Deccan, Columbia River). All literature values recalled, unverified".into(),
            history_rule: "ve_history_v1 (keyed statistical closure, purpose t1b.volcano.ve, counter = blake3(id)[0..4] << 24 | category << 16 | k): edifice lifetime L_e composite 0.5 My, shield 1.0, rift field 2.0, plateau 1.0; mappable-unit window T_vis 40 / 10 / 100 / 1000 ky; history window T_h = min(T_vis, last - onset) (onset absent: T_vis); unit age a = record age + u T_h; counts Poisson(rate T_h), rates per ky: composite lava flows 0.5, explosive VEI >= 3 0.25 with P(VEI >= v) = 10^(-(v - 3) / 2) (v 3..6), domes 0.05 continental arc / 0.025 island arc, flank vents 0.05; shield flows 4.0, rift cones 0.5, fissures 0.6, pit craters 2-4; rift field rows 3-6 with cones 3-10 and flows 2-4 per row; plateau sheets clamp(H / 40 m, 3, 12); lahars 1 (VEI 3), 1-2 (VEI 4), 2-3 (VEI 5-6) per explosive event; sector collapse P 0.15 per composite; caldera at a VEI 6 event of a large old composite system (span >= 1 My, H >= 1500 m); a shield always has a summit caldera; each family capped at max_units (flows first cut). Erosion: every unit loses E(vent) x a (thickness, height or depth; gone when nothing mappable is left), the record-age law of lane CO".into(),
            size_rule: "ve_sizes_v1 (recalled, unverified; REVISION 1 aprons by volume): bulk tephra per explosive event V = 10^(VEI + 4.5) m3 (geometric mean of the VEI class: VEI 4 1e8-1e9 m3); apron density-current share 30% over the apron annulus (r_a to reach, sector), thickness cap = min(cap, 2 x mean) (valley ponding), proximal fall share 20% as T0 exp(-r_e / b) with elliptical isopachs (r_e = hypot(along / 2.5 downwind or along / 0.6 upwind, cross); T0 = 0.2 V / (1.55 2 pi b^2)); density currents weighted by a smooth cos lobe of the declared half-angle, onset smoothstep over 1.5 km beyond r_a, outer taper 0.3 reach;  andesite flows thickness 15-50 m, half-width 100-300 m (x1.5 downstream), length lognormal median 3.5 km sigma 0.5 in [1.5, 8] km, 60% from the summit crater, 40% flank vents at 0.15-0.5 R_f; shield flows 3-12 m, 200-800 m, median 12 km sigma 0.6 in [3, 50] km, 40% summit caldera rim, 60% rift vents; rift-field flows 3-10 m, 200-600 m, median 6 km in [2, 20] km; breached cinder-cone flows 5-15 m, 80-200 m, 1-5 km; cinder cones basal width W_co lognormal median 0.9 km sigma 0.4 in [0.3, 2.0] km (Wood 1980), shield/rift spatter-cinder cones median 0.6 km in [0.2, 1.5]; domes D median 0.8 km in [0.4, 1.5], H/D 0.2-0.35 (St Helens 1980-86 ~0.27); maars D median 0.7 km in [0.3, 1.5], depth 0.15 D, rim 0.04 D; pit craters D 0.3-1.0 km, depth 50-150 m; composite caldera R 2.5-4 km, depth 0.2-0.3 H below the rim, wall 0.8; shield caldera R 1.5-2.5 km, depth 100-200 m, wall 1.0; amphitheatre half-width 0.8-1.2 km, depth 0.15-0.25 H, floor slope 0.15, length 0.5 R_f; aprons VEI 4/5/6: fan slope 0.06/0.04/0.03, thickness cap at most 15/40/100 m (volume rule above), reach 1.2/1.5/2.0 R_f, half-angle 45-90 / 90-180 / 180 deg, fall decay length b 2/4/8 km (Pyle 1989 exponential thinning); lahars thickness 5-20 m, half-width 250-600 m, length median 10 km in [5, 25]; fissures length 0.5-3 km, height 5-15 m, half-width 60 m; plateau sheet radius 25 km x (1 - 0.5 k / N) with 15% keyed lobes".into(),
            shape_rule: "ve_shapes_v1: composite body = lane CF cone (H / W_B 0.15) without CF's lava lobes (VE flows replace them), CF gullies kept as a cut at gully_age; shield body h = H exp(-(r_e / R0)^2), W_B = H / 0.05 clamped [10, 100] km at h = 0.05 H, R0 = (W_B / 2) / sqrt(ln 20), r_e = r (1 - 0.35 sum_k exp(-(lat_k / 2.5 km)^2) [along_k > 0]) along 2 or 3 keyed rift azimuths (rift-stretched shield, Kilauea-type ridges); rift field the same Gaussian with H / W_B 0.08 stretched 2.5x along the keyed axis; broad bodies are at most 4 H x shape thick over the post's own pre-edifice ground (only the distal flank, shape < ~0.4 near deep valleys, is affected: thin distal lavas drape valleys far below the T0-bilinear base instead of standing as a foot wall; build8's 1 x cap made the whole shield follow the 1.5 km T1b relief, measured); plateau no body (sheets); burial and foot fade as lane CF (inner 0.55, outer 0.95, taper 0.15 of the foot radius). Cinder cone: H_co = 0.18 W_co, crater W_cr = 0.40 W_co (Porter 1972 Mauna Kea; Wood 1980), crater depth 0.3 H_co, linear flanks (~31 deg), 50% breached toward the downslope side with a lava flow from the breach; dome h = H_d sqrt(1 - (r / R_d)^4); maar bowl depth d below the ground inside R_m rising to a rim 0.04 D, tuff ring decaying exp(-(r - R_m) / (0.25 D)); pit crater -d smoothstep((R - r) / (d / 0.9)); fissure rampart h exp(-(lat / w)^2) with 200 m end tapers; caldera min(z, floor + wall max(0, r - r_floor)); amphitheatre min(z, floor0 - s along + wall max(0, |lat| - half_w) + wall max(0, -back - along)); plateau sheet t smoothstep((R_k(theta) - r) / m)".into(),
            path_rule: "REVISION 1: units are built on the published window surface only; the layer E / F solve evaluations see the family body without units (the edifice is outside the F6 mean). ve_paths_v1 (chart-wide at site build; pure geometry): from the source a step of path_step along the steepest descent of the analytic body surface (finite differences of 64 m; flat: radially outward) turned by a keyed sinusoidal wander (0.35 x cos(2 pi s / lambda + phi), lambda 1.5-3 km); once the body surface at a network curve point within capture radius is no higher than its bed + 5 m (the old valley floor is exposed) and the bed is not above the path's current level + 5 m, the path is captured (nearest such point; ties by curve then point index) and follows the curve downstream, chain to chain (a chain continues at the curve starting at its last point), its level = the curve bed; uncaptured off the body it descends the T0 bilinear gradient (500 m differences) and stops on a flat; stops at its length, within t1a_stop of a T1a curve (lane CO's corridor keeps T1a valleys, F1), at the sea (bilinear or bed <= sea), at a chain end without continuation".into(),
            deposit_rule: "ve_deposits_v1: per post the nearest path segment gives along-path s, lateral distance and, past the end, the distance to the end point (snout, radius x (1 + 0.25 keyed lobes)); half-width w(s) = w0 (1 + 0.5 s / L), margin m = clamp(t / 0.7, spacing, w) (35 deg fronts), thickness t smoothstep((w - lat) / m) x min(1, s / 200 m); top = level(s) + thickness in a captured valley (fills the floor, the walls confine it), z + thickness on the body or open ground (lahars: only beyond 0.8 of the body radius or in a captured valley: lahar fans on the lower flank and ring plain); z = max(z, min(top, z + 3 t)). Cones, domes, maars, fissures, sheets add to the post's height; calderas, amphitheatres, pits cut with min(); aprons: z += fall, then z = max(z, min(F, z + t_cap g taper)) with F = body(r_a, theta) - s_f (r - r_a), r_a where the body slope falls to 0.15 (at most 0.85 R_f), g the sector weight; units oldest first; posts raised by units obey lane CO's T1a corridor cap with the body (relief over the pre-edifice surface <= cap_slope (d - mask_inner)); T1a channel, lake and non-finite posts untouched; network channel posts raised > 1 m are buried (free); marine posts raised above sea become land; over = 1, volcano_id, style and unit on posts a unit changes by >= t_map (thinner veneers change the height only)".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.flank2.is_some() && self.flank.is_none() { return Err("ve_flank2 needs ve_flank=1".into()); }
        if self.flank.as_ref().is_some_and(|f| f.rule != FLANK_RULE) { return Err("ve_flank: unknown rule".into()); }
        if let Some(e) = &self.erec {
            if e.rule != EREC_AGE_RULE || e.record.frame_step != 200 || e.record.floor_m_per_my != 5 { return Err("ve_erec: unknown age rule, frame or floor".into()); }
            e.record.validate()?;
        }
        if self.ve3.is_some() && self.ve2.as_ref().map(|v| v.steer != 1).unwrap_or(true) { return Err("ve3 needs ve2 steer".into()); }
        if self.ve4.is_some() && self.ve3.as_ref().map(|v| v.rough != 1).unwrap_or(true) { return Err("ve4 needs ve3 rough".into()); }
        if self.ve5.is_some() && self.ve4.as_ref().map(|v| v.lobes != 1).unwrap_or(true) { return Err("ve5 needs ve4 lobes".into()); }
        if self.on > 1 || self.t_map_mm < 0 || self.path_step_mm < 16_000 || self.capture_mm < 0 || self.t1a_stop_mm < 0 || self.max_units == 0 { return Err("ve: bounds".into()); }
        Ok(())
    }
}

pub fn ve(p: &crate::t1b_rivernet::RiverNetParamsV1) -> Option<&VeParamsV1> { crate::t1b_ls::tc(p).and_then(|t| t.ve.as_ref()).filter(|v| v.on == 1) }

// ------------------------------------------------------------------ per-volcano history

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Family { Composite, Shield, Rift, Plateau }

pub fn family_of(kind: &str) -> Option<Family> {
    match kind {
        "island_arc" | "continental_arc" => Some(Family::Composite),
        "hotspot" => Some(Family::Shield),
        "rift" => Some(Family::Rift),
        "large_igneous_province" => Some(Family::Plateau),
        _ => None,
    }
}
impl Family {
    pub fn name(self) -> &'static str { match self { Family::Composite => "composite", Family::Shield => "shield", Family::Rift => "rift_field", Family::Plateau => "flood_plateau" } }
    /// (edifice lifetime My, mappable window ky).
    pub fn windows(self) -> (f64, f64) { match self { Family::Composite => (0.5, 40.0), Family::Shield => (1.0, 10.0), Family::Rift => (2.0, 100.0), Family::Plateau => (1.0, 1000.0) } }
    pub fn body_style(self) -> u8 { match self { Family::Composite => ST_BODY, Family::Shield | Family::Rift => ST_SHIELD, Family::Plateau => ST_PLATEAU } }
    /// Broad families: CO clearance plus the shield clearance (not CF's foot clearance).
    pub fn broad(self) -> bool { !matches!(self, Family::Composite) }
}

/// Unit geometry (chart metres; absolute heights in metres).
#[derive(Clone, Debug, Serialize)]
#[serde(tag = "g")]
pub enum VeGeom {
    /// Path samples (x, y, s, level): level is the captured valley bed, LV_BODY / LV_OPEN for drape.
    Path { pts: Vec<[f64; 4]>, len_m: f64, w0_m: f64, t_m: f64, lahar: bool, lobes: [f64; 4], captured_at_m: f64, end: String,
        /// Lane VE3 ponded lava along the path (absent without ve3 pond).
        #[serde(skip_serializing_if = "Vec::is_empty")]
        ponds: Vec<Pond> },
    Cone { x: f64, y: f64, rb_m: f64, h_m: f64, rc_m: f64, dc_m: f64, breach_rad: Option<f64> },
    Dome { x: f64, y: f64, r_m: f64, h_m: f64,
        /// Lane VE2 dome cap: the crater rim the dome top may not exceed.
        #[serde(skip_serializing_if = "Option::is_none")]
        top_m: Option<f64> },
    Maar { x: f64, y: f64, r_m: f64, d_m: f64, rim_m: f64, ring_b_m: f64 },
    Pit { x: f64, y: f64, r_m: f64, d_m: f64 },
    Fissure { x0: f64, y0: f64, x1: f64, y1: f64, h_m: f64, w_m: f64 },
    Caldera { x: f64, y: f64, r_m: f64, floor_m: f64, r_floor_m: f64, wall: f64 },
    Breach { x: f64, y: f64, dir_rad: f64, half_w_m: f64, floor0_m: f64, slope: f64, len_m: f64, back_m: f64, wall: f64 },
    Apron { x: f64, y: f64, dir_rad: f64, half_rad: f64, r_a_m: f64, r_max_m: f64, s_f: f64, t_cap_m: f64, t0_m: f64, b_m: f64, za: Vec<f64> },
    Sheet { x: f64, y: f64, r_m: f64, t_m: f64, lobes: [f64; 4] },
    /// CF radial gullies (composite), cut at gully_age.
    Gullies { n: u32 },
}

#[derive(Clone, Debug, Serialize)]
pub struct VeUnit {
    pub id: u16,
    pub kind: String,
    pub style: u8,
    pub age_yr: f64,
    /// 0 = effusive or not an eruption product.
    pub vei: u8,
    /// Declared size before erosion (thickness, height or depth, m) and the erosion depth applied.
    pub size_m: f64,
    pub erosion_m: f64,
    pub geom: VeGeom,
    /// Lane VE2 steering: a path traced after the history is sorted (start, length, wander); never serialized.
    #[serde(skip)]
    pub pending: Option<([f64; 2], f64, [f64; 3])>,
}

#[derive(Clone, Debug, Serialize)]
pub struct VeVolcano {
    pub family: String,
    /// Arc2 construction/erosion mode. Absent preserves all Arc1 bytes.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flank2: Option<Flank2Site>,
    /// Native, chart-independent composite drainage field. Absent in legacy output.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flank: Option<crate::t1b_coast::FlankSurface>,
    pub span_my: Option<f64>,
    pub t_h_ky: f64,
    pub lifetime_my: f64,
    /// Body: foot radius (CF R_f or shield W_B / 2), Gaussian scale, stretch axes (rad), stretch factor, axis half-width.
    pub r_body_m: f64,
    pub r0_m: f64,
    pub axes_rad: Vec<f64>,
    pub stretch: f64,
    pub axis_w_m: f64,
    pub elong: f64,
    pub q_built: f64,
    /// Units reach (from the vent) for window selection and the T1a corridor segments.
    pub r_units_m: f64,
    /// CF gullies moved here (composite): (angle, r0, depth, wander amplitude, wander length, phase).
    pub gullies: Vec<[f64; 6]>,
    /// Lane VE3 lobe roughness (absent without ve3 rough).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rough: Option<RoughSpec>,
    pub counts: BTreeMap<String, u32>,
    pub dropped_eroded: u32,
    pub dropped_cap: u32,
    pub units: Vec<VeUnit>,
    /// Lane VH-I (absent without ve_wind): the T0 wind at the record and the downwind axis.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind: Option<VeWindAt>,
}

// ------------------------------------------------------------------ keyed draws

struct Kr<'a> { sc: &'a Scope, base: u64, k: u64, real: u32 }
impl<'a> Kr<'a> {
    fn new(sc: &'a Scope, id: &str, cat: u64, real: u32) -> Self {
        let hb = blake3::hash(id.as_bytes());
        let h = u32::from_le_bytes([hb.as_bytes()[0], hb.as_bytes()[1], hb.as_bytes()[2], hb.as_bytes()[3]]) as u64;
        Kr { sc, base: (h << 24) | (cat << 16), k: 0, real }
    }
    fn u(&mut self) -> R<f64> { let d = self.sc.draw(PURPOSE_VE, self.base + self.k, self.real).map_err(|e| e.to_string())?; self.k += 1; Ok(d.uniform_f64()) }
    fn n(&mut self) -> R<f64> { let d = self.sc.draw(PURPOSE_VE, self.base + self.k, self.real).map_err(|e| e.to_string())?; self.k += 1; Ok(d.normal_f64()) }
    fn range(&mut self, a: f64, b: f64) -> R<f64> { Ok(a + (b - a) * self.u()?) }
    fn logn(&mut self, med: f64, sig: f64, lo: f64, hi: f64) -> R<f64> { Ok((med * (sig * self.n()?).exp()).clamp(lo, hi)) }
    /// Poisson by inversion of one uniform.
    fn poisson(&mut self, mean: f64) -> R<u32> {
        let u = self.u()?;
        let (mut k, mut p) = (0u32, (-mean).exp());
        let mut c = p;
        while u > c && k < 10_000 { k += 1; p *= mean / k as f64; c += p; }
        Ok(k)
    }
    fn lobes(&mut self) -> R<[f64; 4]> { Ok([self.u()? * std::f64::consts::TAU, self.u()? * std::f64::consts::TAU, self.u()? * std::f64::consts::TAU, self.u()? * std::f64::consts::TAU]) }
}

#[inline]
fn lobe(l: &[f64; 4], th: f64) -> f64 { ((2.0 * th + l[0]).cos() / 2.0 + (3.0 * th + l[1]).cos() / 3.0 + (4.0 * th + l[2]).cos() / 4.0 + (5.0 * th + l[3]).cos() / 5.0) / 1.2833 }

// ------------------------------------------------------------------ body surfaces

/// Sloped base of lane CF (R1a; C1-smooth min, k = 100 m).
#[inline]
fn sloped(bxy: f64, base: f64) -> f64 { let t = bxy - base; let k = 100.0; if t >= 0.0 { 0.0 } else if t > -2.0 * k { -t * t / (4.0 * k) } else { t + k } }

/// Normalised body height (0..1) of a broad family at offset (dx, dy) from the vent.
fn broad_shape(g: &VeVolcano, dx: f64, dy: f64) -> f64 {
    let r = dx.hypot(dy);
    if g.elong > 1.0 {
        let a = g.axes_rad.first().copied().unwrap_or(0.0);
        let (al, la) = (dx * a.cos() + dy * a.sin(), -dx * a.sin() + dy * a.cos());
        let re = (al / g.elong).hypot(la);
        return (-(re / g.r0_m).powi(2)).exp();
    }
    let mut f = 0.0;
    for &a in &g.axes_rad {
        let (al, la) = (dx * a.cos() + dy * a.sin(), -dx * a.sin() + dy * a.cos());
        if al > 0.0 { f += (-(la / g.axis_w_m).powi(2)).exp(); }
    }
    let re = r * (1.0 - g.stretch * f.min(1.0));
    (-(re / g.r0_m).powi(2)).exp()
}

/// Lane VE body at chart (x, y): (E surface, burial weight w, foot fade u), as lane CF's `cone_cf`. None = no body here.
pub fn body(vp: &VolcanoParamsV1, v: &VolcanoSite, cf: &CfParamsV1, x: f64, y: f64, bxy: f64) -> Option<(f64, f64, f64)> {
    let g = v.ve.as_ref()?;
    match family_of(&v.kind)? {
        Family::Composite => cone_cf_opt(vp, v, cf, x, y, bxy, false).map(|(e, w, u)| (e + g.rough.as_ref().map(|r| rough_at(r, x, y)).unwrap_or(0.0), w, u)),
        Family::Plateau => None,
        Family::Shield | Family::Rift => {
            let (dx, dy) = (x - v.vent_xy_m[0], y - v.vent_xy_m[1]);
            let r = dx.hypot(dy);
            let rf = g.r_body_m;
            let scale = if v.h_m > 0.0 { ((v.h_m - v.denudation_m) / v.h_m).max(0.0) } else { 0.0 };
            let e = v.base_m + v.h_m * scale * broad_shape(g, dx, dy) + if cf.sloped_base == 1 { sloped(bxy, v.base_m) } else { 0.0 };
            let e = match g.rough.as_ref() { Some(r) => e + rough_at(r, x, y), None => e };
            // The body's own extent along a stretched axis is longer: burial and fade by the equivalent radius of the
            // normalised height (r for an unstretched Gaussian).
            let rn = if g.elong > 1.0 || g.stretch > 0.0 { let s = broad_shape(g, dx, dy).max(1e-300); g.r0_m * (-s.ln()).sqrt() } else { r };
            let (bi, bo) = (cf.bury_inner_permille as f64 / 1000.0 * rf, cf.bury_outer_permille as f64 / 1000.0 * rf);
            let w = if rn <= bi { 1.0 } else if rn >= bo { 0.0 } else { smoothstep01((bo - rn) / (bo - bi)) };
            let th = dy.atan2(dx);
            let lob = if v.lobe_phase_rad.is_empty() { 0.0 } else { v.lobe_phase_rad.iter().enumerate().map(|(j, p)| ((j + 2) as f64 * th + p).cos() / (j + 2) as f64).sum::<f64>() / v.lobe_norm };
            let r_end = rf * (1.0 + cf.foot_taper_permille as f64 / 1000.0 * (1.0 + 0.8 * lob));
            let u = if rn <= rf { 1.0 } else if r_end <= rf + 1.0 { 0.0 } else { smoothstep01((r_end - rn) / (r_end - rf)) };
            Some((e, w, u))
        }
    }
}

/// Broad families: the body's thickness over the post's own pre-edifice ground is at most 4 H x shape (the Gaussian height
/// above its base): distal lavas drape the ground they reach instead of standing at the base level over deep valleys
/// (no foot wall where the T1b ground lies far below the T0 bilinear). Infinite for composites and plateaus.
pub fn thickness_cap(v: &VolcanoSite, x: f64, y: f64) -> f64 {
    let Some(g) = v.ve.as_ref() else { return f64::INFINITY };
    match family_of(&v.kind) {
        Some(Family::Shield) | Some(Family::Rift) => {
            let scale = if v.h_m > 0.0 { ((v.h_m - v.denudation_m) / v.h_m).max(0.0) } else { 0.0 };
            4.0 * v.h_m * scale * broad_shape(g, x - v.vent_xy_m[0], y - v.vent_xy_m[1])
        }
        _ => f64::INFINITY,
    }
}

/// The CF gully cut (depth, m) of a composite at chart (x, y) (lane CF's radial_gullies_young_v1 formula).
pub fn gully_cut(g: &VeVolcano, v: &VolcanoSite, cf: &CfParamsV1, x: f64, y: f64) -> f64 {
    if g.flank.is_some() || g.flank2.is_some() { return 0.0; }
    let (dx, dy) = (x - v.vent_xy_m[0], y - v.vent_xy_m[1]);
    let r = dx.hypot(dy);
    let th = dy.atan2(dx);
    let rf = g.r_body_m;
    let mut cut: f64 = 0.0;
    for q in &g.gullies {
        if r <= q[1] || r >= rf { continue; }
        let d = ang_diff(th, q[0]);
        if d.abs() >= std::f64::consts::FRAC_PI_2 { continue; }
        let off = q[3] * (std::f64::consts::TAU * r / q[4].max(1.0) + q[5]).sin();
        let lat = (r * d.sin() - off).abs();
        let hw = cf.gully_halfwidth_mm as f64 / 1000.0 + cf.gully_halfwidth_ppm_of_r as f64 / 1e6 * r;
        if lat >= hw { continue; }
        let depth = q[2] * (std::f64::consts::PI * (r - q[1]) / (rf - q[1])).sin();
        cut = cut.max(depth * (1.0 - lat / hw));
    }
    cut
}

// ------------------------------------------------------------------ build (chart-wide, at site construction)

/// Everything the history build reads from the chart.
pub struct BuildIn<'a> {
    pub vp: &'a VolcanoParamsV1,
    pub cf: &'a CfParamsV1,
    pub ve: &'a VeParamsV1,
    pub t0: &'a T0,
    pub chart: &'a Chart,
    pub curves: &'a [Curve],
    pub table: &'a BTreeMap<u32, f64>,
    pub seed: u64,
    pub recipe_digest: &'a str,
    pub realization: u32,
    /// Record onset (My) per volcano id.
    pub onset: &'a BTreeMap<String, f64>,
    pub t_my: f64,
}

/// Network curve index for path capture.
struct CurveIx { cell: f64, grid: BTreeMap<(i64, i64), Vec<(u32, u32)>>, starts: BTreeMap<(i64, i64), Vec<u32>>, t1a: Vec<[f64; 4]> }
impl CurveIx {
    fn new(curves: &[Curve]) -> Self {
        let cell = 1000.0;
        let mut grid: BTreeMap<(i64, i64), Vec<(u32, u32)>> = BTreeMap::new();
        let mut starts: BTreeMap<(i64, i64), Vec<u32>> = BTreeMap::new();
        let mut t1a = Vec::new();
        for (ci, c) in curves.iter().enumerate() {
            if c.code == CHANNEL { for w in c.pts.windows(2) { t1a.push([w[0][0], w[0][1], w[1][0], w[1][1]]); } }
            if c.code != GUIDE && c.code != CHANNEL { continue; }
            for (pi, p) in c.pts.iter().enumerate() { grid.entry(((p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64)).or_default().push((ci as u32, pi as u32)); }
            if c.code == GUIDE { if let Some(p) = c.pts.first() { starts.entry((p[0].round() as i64, p[1].round() as i64)).or_default().push(ci as u32); } }
        }
        CurveIx { cell, grid, starts, t1a }
    }
    fn t1a_dist(&self, p: [f64; 2]) -> f64 {
        let mut d = f64::INFINITY;
        for s in &self.t1a { let q = seg_nearest(p, [s[0], s[1]], [s[2], s[3]]); d = d.min((q[0] - p[0]).hypot(q[1] - p[1])); }
        d
    }
}

struct BodyEval<'a> { b: &'a BuildIn<'a>, v: &'a VolcanoSite }
impl<'a> BodyEval<'a> {
    fn bil(&self, x: f64, y: f64) -> f64 { bil_at(self.b.t0, self.b.chart, x, y).unwrap_or(f64::NAN) }
    /// Body surface (with the sloped base) or None outside its fade.
    fn e(&self, x: f64, y: f64) -> Option<f64> {
        let bx = self.bil(x, y);
        if !bx.is_finite() { return None; }
        // Lane T1B-PEDESTAL: the sloped base keeps t = B(post) - B(vent) (ped_shift 0 without the key).
        body(self.b.vp, self.v, self.b.cf, x, y, bx - self.v.ped_shift()).and_then(|(e, _, u)| if u > 0.0 { Some(e) } else { None })
    }
    fn e_any(&self, x: f64, y: f64) -> Option<f64> { let bx = self.bil(x, y); body(self.b.vp, self.v, self.b.cf, x, y, bx - self.v.ped_shift()).map(|t| t.0) }
    fn grad(&self, x: f64, y: f64) -> Option<[f64; 2]> {
        let h = 64.0;
        let (a, b, c, d) = (self.e_any(x + h, y)?, self.e_any(x - h, y)?, self.e_any(x, y + h)?, self.e_any(x, y - h)?);
        Some([(a - b) / (2.0 * h), (c - d) / (2.0 * h)])
    }
}

/// Trace a path (see path_rule). Returns (pts, captured_at_m, end reason).
#[allow(clippy::too_many_arguments)]
fn trace(be: &BodyEval, ix: &CurveIx, curves: &[Curve], start: [f64; 2], len: f64, wander: [f64; 3], r_body: f64, sea: f64) -> (Vec<[f64; 4]>, f64, String) {
    trace_g(be, ix, curves, start, len, wander, r_body, sea, None)
}

/// As `trace`; with `extra` (lane VE2 steering) the descent runs on body + extra (older units, gullies) and a rising step ponds.
#[allow(clippy::too_many_arguments)]
fn trace_g(be: &BodyEval, ix: &CurveIx, curves: &[Curve], start: [f64; 2], len: f64, wander: [f64; 3], r_body: f64, sea: f64, extra: Option<&dyn Fn(f64, f64) -> f64>) -> (Vec<[f64; 4]>, f64, String) {
    let ve = be.b.ve;
    let step = ve.path_step_mm as f64 / 1000.0;
    let cap = ve.capture_mm as f64 / 1000.0;
    let t1s = ve.t1a_stop_mm as f64 / 1000.0;
    let vent = be.v.vent_xy_m;
    let mut p = start;
    let mut s = 0.0;
    let on_body = |p: [f64; 2]| be.e(p[0], p[1]).is_some();
    let lv0 = if on_body(p) { LV_BODY } else { LV_OPEN };
    let mut pts = vec![[p[0], p[1], 0.0, lv0]];
    let mut level = be.e(p[0], p[1]).unwrap_or_else(|| be.bil(p[0], p[1]));
    let mut dir_prev = { let d = [p[0] - vent[0], p[1] - vent[1]]; let l = d[0].hypot(d[1]); if l > 1.0 { [d[0] / l, d[1] / l] } else { [1.0, 0.0] } };
    let mut guard = 0;
    while s < len && guard < 4000 {
        guard += 1;
        // Capture (off the upper body).
        let rr = (p[0] - vent[0]).hypot(p[1] - vent[1]);
        if rr > 0.5 * r_body {
            let (cx0, cy0) = ((p[0] / ix.cell).floor() as i64, (p[1] / ix.cell).floor() as i64);
            let mut best: Option<(f64, u32, u32)> = None;
            let mut t1a_hit = false;
            for gy in cy0 - 1..=cy0 + 1 { for gx in cx0 - 1..=cx0 + 1 {
                let Some(l) = ix.grid.get(&(gx, gy)) else { continue };
                for &(ci, pi) in l {
                    let c = &curves[ci as usize];
                    let q = c.pts[pi as usize];
                    let d = (q[0] - p[0]).hypot(q[1] - p[1]);
                    if c.code == CHANNEL { if d < t1s { t1a_hit = true; } continue; }
                    if d > cap { continue; }
                    if q[2] > level + 5.0 { continue; }
                    if let Some(eb) = be.e(q[0], q[1]) { if eb > q[2] + 5.0 { continue; } }
                    if best.map(|b| (d, ci, pi) < b).unwrap_or(true) { best = Some((d, ci, pi)); }
                }
            } }
            if t1a_hit || ix.t1a_dist(p) < t1s { return (pts, f64::NAN, "t1a_corridor".into()); }
            if let Some((_, ci, pi)) = best {
                let captured_at = s;
                let (mut ci, mut pi) = (ci as usize, pi as usize);
                let mut hops = 0;
                loop {
                    let c = &curves[ci];
                    for q in &c.pts[pi..] {
                        let d = (q[0] - p[0]).hypot(q[1] - p[1]);
                        if d < 1e-6 { continue; }
                        if s + d >= len {
                            let f = (len - s) / d;
                            pts.push([p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), len, q[2]]);
                            return (pts, captured_at, "length".into());
                        }
                        s += d;
                        p = [q[0], q[1]];
                        pts.push([q[0], q[1], s, q[2]]);
                        if q[2] <= sea { return (pts, captured_at, "sea".into()); }
                        if ix.t1a_dist(p) < t1s { return (pts, captured_at, "t1a_corridor".into()); }
                    }
                    let last = c.pts[c.pts.len() - 1];
                    let key = (last[0].round() as i64, last[1].round() as i64);
                    let next = ix.starts.get(&key).and_then(|v| v.iter().copied().find(|&k| k as usize != ci));
                    hops += 1;
                    match next { Some(k) if hops < 400 => { ci = k as usize; pi = 0; } _ => return (pts, captured_at, "chain_end".into()) }
                }
            }
        }
        // Steepest descent (body, else the T0 bilinear), with wander.
        let lower = (p[0] - vent[0]).hypot(p[1] - vent[1]) >= 0.8 * r_body;
        let gx = |x: f64, y: f64| -> [f64; 2] { match extra { Some(f) => { let h = 64.0; [(f(x + h, y) - f(x - h, y)) / (2.0 * h), (f(x, y + h) - f(x, y - h)) / (2.0 * h)] } None => [0.0, 0.0] } };
        let (g, lv) = match (on_body(p), be.grad(p[0], p[1])) {
            (true, Some(g)) => { let e = gx(p[0], p[1]); ([g[0] + e[0], g[1] + e[1]], if lower { LV_OPEN } else { LV_BODY }) }
            _ => {
                let h = 500.0;
                let mut g = [(be.bil(p[0] + h, p[1]) - be.bil(p[0] - h, p[1])) / (2.0 * h), (be.bil(p[0], p[1] + h) - be.bil(p[0], p[1] - h)) / (2.0 * h)];
                let e = gx(p[0], p[1]);
                g = [g[0] + e[0], g[1] + e[1]];
                if !g[0].is_finite() || !g[1].is_finite() || g[0].hypot(g[1]) < 1e-4 { return (pts, f64::NAN, "flat".into()); }
                (g, LV_OPEN)
            }
        };
        let gl = g[0].hypot(g[1]);
        let mut d = if gl < 1e-6 { dir_prev } else { [-g[0] / gl, -g[1] / gl] };
        let wv = wander[0] * (std::f64::consts::TAU * s / wander[1].max(1.0) + wander[2]).cos();
        d = [d[0] - wv * d[1], d[1] + wv * d[0]];
        let dl = d[0].hypot(d[1]);
        d = [d[0] / dl, d[1] / dl];
        dir_prev = d;
        let st = step.min(len - s);
        // Lane VE2: a step that rises on body + older units is a closed low: the flow ponds there.
        if let (Some(f), true) = (extra, on_body(p)) {
            let gz = |q: [f64; 2]| be.e_any(q[0], q[1]).unwrap_or(f64::NAN) + f(q[0], q[1]);
            let (z0, z1) = (gz(p), gz([p[0] + st * d[0], p[1] + st * d[1]]));
            if z0.is_finite() && z1.is_finite() && z1 > z0 && s > 0.0 { return (pts, f64::NAN, "ponded".into()); }
        }
        p = [p[0] + st * d[0], p[1] + st * d[1]];
        s += st;
        pts.push([p[0], p[1], s, lv]);
        level = match be.e(p[0], p[1]) { Some(e) => e, None => be.bil(p[0], p[1]) };
        if be.bil(p[0], p[1]) <= sea - 5.0 && !on_body(p) { return (pts, f64::NAN, "sea".into()); }
    }
    (pts, f64::NAN, "length".into())
}

/// Lane VE3: as `trace_g` on the lifted steering surface, with ponds (flat fill to the spill level or the volume, overflow).
#[allow(clippy::too_many_arguments)]
fn trace_g3(be: &BodyEval, ix: &CurveIx, curves: &[Curve], start: [f64; 2], len: f64, wander: [f64; 3], r_body: f64, sea: f64, extra: Option<&dyn Fn(f64, f64) -> f64>, pcfg: Option<PondCfg>) -> (Vec<[f64; 4]>, f64, String, Vec<Pond>) {
    let mut ponds: Vec<Pond> = Vec::new();
    let mut npond = 0u32;
    let mut lift: BTreeMap<(i64, i64), f64> = BTreeMap::new();
    // The steering surface with the ponds already filled by this flow.
    let gsurf = |q: [f64; 2], lift: &BTreeMap<(i64, i64), f64>| -> f64 {
        let z = be.e_any(q[0], q[1]).unwrap_or(f64::NAN) + extra.map(|f| f(q[0], q[1])).unwrap_or(0.0);
        match lift.get(&((q[0] / 32.0).round() as i64, (q[1] / 32.0).round() as i64)) { Some(l) if *l > z => *l, _ => z }
    };
    let ve = be.b.ve;
    let step = ve.path_step_mm as f64 / 1000.0;
    let cap = ve.capture_mm as f64 / 1000.0;
    let t1s = ve.t1a_stop_mm as f64 / 1000.0;
    let vent = be.v.vent_xy_m;
    let mut p = start;
    let mut s = 0.0;
    let on_body = |p: [f64; 2]| be.e(p[0], p[1]).is_some();
    let lv0 = if on_body(p) { LV_BODY } else { LV_OPEN };
    let mut pts = vec![[p[0], p[1], 0.0, lv0]];
    let mut level = be.e(p[0], p[1]).unwrap_or_else(|| be.bil(p[0], p[1]));
    let mut dir_prev = { let d = [p[0] - vent[0], p[1] - vent[1]]; let l = d[0].hypot(d[1]); if l > 1.0 { [d[0] / l, d[1] / l] } else { [1.0, 0.0] } };
    let mut guard = 0;
    while s < len && guard < 4000 {
        guard += 1;
        // Capture (off the upper body).
        let rr = (p[0] - vent[0]).hypot(p[1] - vent[1]);
        if rr > 0.5 * r_body {
            let (cx0, cy0) = ((p[0] / ix.cell).floor() as i64, (p[1] / ix.cell).floor() as i64);
            let mut best: Option<(f64, u32, u32)> = None;
            let mut t1a_hit = false;
            for gy in cy0 - 1..=cy0 + 1 { for gx in cx0 - 1..=cx0 + 1 {
                let Some(l) = ix.grid.get(&(gx, gy)) else { continue };
                for &(ci, pi) in l {
                    let c = &curves[ci as usize];
                    let q = c.pts[pi as usize];
                    let d = (q[0] - p[0]).hypot(q[1] - p[1]);
                    if c.code == CHANNEL { if d < t1s { t1a_hit = true; } continue; }
                    if d > cap { continue; }
                    if q[2] > level + 5.0 { continue; }
                    if let Some(eb) = be.e(q[0], q[1]) { if eb > q[2] + 5.0 { continue; } }
                    if best.map(|b| (d, ci, pi) < b).unwrap_or(true) { best = Some((d, ci, pi)); }
                }
            } }
            if t1a_hit || ix.t1a_dist(p) < t1s { return (pts, f64::NAN, "t1a_corridor".into(), ponds); }
            if let Some((_, ci, pi)) = best {
                let captured_at = s;
                let (mut ci, mut pi) = (ci as usize, pi as usize);
                let mut hops = 0;
                loop {
                    let c = &curves[ci];
                    for q in &c.pts[pi..] {
                        let d = (q[0] - p[0]).hypot(q[1] - p[1]);
                        if d < 1e-6 { continue; }
                        if s + d >= len {
                            let f = (len - s) / d;
                            pts.push([p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), len, q[2]]);
                            return (pts, captured_at, "length".into(), ponds);
                        }
                        s += d;
                        p = [q[0], q[1]];
                        pts.push([q[0], q[1], s, q[2]]);
                        if q[2] <= sea { return (pts, captured_at, "sea".into(), ponds); }
                        if ix.t1a_dist(p) < t1s { return (pts, captured_at, "t1a_corridor".into(), ponds); }
                    }
                    let last = c.pts[c.pts.len() - 1];
                    let key = (last[0].round() as i64, last[1].round() as i64);
                    let next = ix.starts.get(&key).and_then(|v| v.iter().copied().find(|&k| k as usize != ci));
                    hops += 1;
                    match next { Some(k) if hops < 400 => { ci = k as usize; pi = 0; } _ => return (pts, captured_at, "chain_end".into(), ponds) }
                }
            }
        }
        // Steepest descent (body, else the T0 bilinear), with wander.
        let lower = (p[0] - vent[0]).hypot(p[1] - vent[1]) >= 0.8 * r_body;
        let gx = |x: f64, y: f64| -> [f64; 2] { match extra { Some(f) => { let h = 64.0; [(f(x + h, y) - f(x - h, y)) / (2.0 * h), (f(x, y + h) - f(x, y - h)) / (2.0 * h)] } None => [0.0, 0.0] } };
        let (g, lv) = match (on_body(p), be.grad(p[0], p[1])) {
            (true, Some(_)) => { let h = 64.0; let f = |q: [f64; 2]| gsurf(q, &lift);
                ([(f([p[0] + h, p[1]]) - f([p[0] - h, p[1]])) / (2.0 * h), (f([p[0], p[1] + h]) - f([p[0], p[1] - h])) / (2.0 * h)], if lower { LV_OPEN } else { LV_BODY }) }
            _ => {
                let h = 500.0;
                let mut g = [(be.bil(p[0] + h, p[1]) - be.bil(p[0] - h, p[1])) / (2.0 * h), (be.bil(p[0], p[1] + h) - be.bil(p[0], p[1] - h)) / (2.0 * h)];
                let e = gx(p[0], p[1]);
                g = [g[0] + e[0], g[1] + e[1]];
                if !g[0].is_finite() || !g[1].is_finite() || g[0].hypot(g[1]) < 1e-4 { return (pts, f64::NAN, "flat".into(), ponds); }
                (g, LV_OPEN)
            }
        };
        let gl = g[0].hypot(g[1]);
        let mut d = if gl < 1e-6 { dir_prev } else { [-g[0] / gl, -g[1] / gl] };
        let wv = wander[0] * (std::f64::consts::TAU * s / wander[1].max(1.0) + wander[2]).cos();
        d = [d[0] - wv * d[1], d[1] + wv * d[0]];
        let dl = d[0].hypot(d[1]);
        d = [d[0] / dl, d[1] / dl];
        dir_prev = d;
        let st = step.min(len - s);
        // Lane VE2: a step that rises on body + older units is a closed low: the flow ponds there.
        if on_body(p) {
            // Four 1/4 sub-steps: a 128 m step on a steep flank drops more than an older flow's margin is high, so it could hop
            // over a levee; any rising sub-step is a closed low (stop at the last falling sub-step).
            let mut z0 = gsurf(p, &lift);
            let mut z1 = z0;
            for k in 1..=4 {
                let q = [p[0] + st * d[0] * k as f64 / 4.0, p[1] + st * d[1] * k as f64 / 4.0];
                z1 = gsurf(q, &lift);
                if z1.is_finite() && z0.is_finite() && z1 > z0 { if k > 1 { let qb = [p[0] + st * d[0] * (k - 1) as f64 / 4.0, p[1] + st * d[1] * (k - 1) as f64 / 4.0]; s += st * (k - 1) as f64 / 4.0; p = qb; pts.push([p[0], p[1], s, lv]); } z0 = gsurf(p, &lift); break; }
                if k < 4 { z0 = z1; }
            }
            if z0.is_finite() && z1.is_finite() && z1 > z0 && s > 0.0 {
                // Lane VE3: fill the closed low (flat, to the spill level or the volume), then overflow from the spill point.
                let Some(pc) = pcfg else { return (pts, f64::NAN, "ponded".into(), ponds) };
                if z0 <= sea + 1.0 { return (pts, f64::NAN, "sea".into(), ponds); }
                if npond >= pc.max { return (pts, f64::NAN, "ponded".into(), ponds); }
                npond += 1;
                let vol_av = (len - s).max(0.0) * 2.0 * pc.w * pc.t;
                let (pond, spill) = fill_pond_e(&|q: [f64; 2]| gsurf(q, &lift), p, pc.radius, vol_av, pc.edge_spill);
                for r in &pond.rows { for ixx in r[1]..=r[2] { lift.insert((ixx, r[0]), pond.level_m); } }
                let used = pond.volume_m3 / (2.0 * pc.w * pc.t);
                let full = pond.volume_m3 >= vol_av * (1.0 - 1e-9);
                if pond.cells > 0 { ponds.push(pond); }
                match spill { Some(q) if !full => {
                    s += used;
                    if s >= len { return (pts, f64::NAN, "ponded".into(), ponds); }
                    p = q; s += 1e-6;
                    pts.push([p[0], p[1], s, lv]);
                    continue;
                }
                _ => return (pts, f64::NAN, "ponded".into(), ponds) }
            }
        }
        p = [p[0] + st * d[0], p[1] + st * d[1]];
        s += st;
        pts.push([p[0], p[1], s, lv]);
        level = match be.e(p[0], p[1]) { Some(e) => e, None => be.bil(p[0], p[1]) };
        if be.bil(p[0], p[1]) <= sea - 5.0 && !on_body(p) { return (pts, f64::NAN, "sea".into(), ponds); }
    }
    (pts, f64::NAN, "length".into(), ponds)
}

/// Lane VE2: the drape height a built unit adds at (x, y) (paths: the deposit profile; cones, domes, fissures), for steering.
pub fn unit_drape(u: &VeUnit, x: f64, y: f64) -> f64 { unit_drape_m(u, x, y, 32.0) }

/// As `unit_drape`, with the path margin floor `mfloor` (lane VE3 aa: 1.5 x 32 m).
pub fn unit_drape_m(u: &VeUnit, x: f64, y: f64, mfloor: f64) -> f64 {
    let a = u.erosion_m;
    match &u.geom {
        VeGeom::Path { pts, len_m, w0_m, t_m, lobes, .. } => {
            let t = t_m - a;
            if t <= 0.0 || pts.len() < 2 { return 0.0; }
            let wmax = w0_m * 1.5 * 1.3 + 64.0;
            let n = pts.len() - 1;
            let (mut best, mut bs, mut bb, mut ba) = (f64::INFINITY, 0.0, 0.0, 0.0);
            for k in 0..n {
                let (p, q) = (pts[k], pts[k + 1]);
                if x < p[0].min(q[0]) - wmax || x > p[0].max(q[0]) + wmax || y < p[1].min(q[1]) - wmax || y > p[1].max(q[1]) + wmax { continue; }
                let (dx, dy) = (q[0] - p[0], q[1] - p[1]);
                let l2 = dx * dx + dy * dy;
                if l2 <= 0.0 { continue; }
                let ln = l2.sqrt();
                let tt = ((x - p[0]) * dx + (y - p[1]) * dy) / l2;
                let (d, s, be, an) = if k == n - 1 && tt > 1.0 { let (ex, ey) = (x - q[0], y - q[1]); (ex.hypot(ey), q[2], (ex * dx + ey * dy) / ln, ey.atan2(ex) - dy.atan2(dx)) }
                    else { let tc = tt.clamp(0.0, 1.0); ((x - p[0] - tc * dx).hypot(y - p[1] - tc * dy), p[2] + tc * (q[2] - p[2]), 0.0, 0.0) };
                if d < best { best = d; bs = s; bb = be; ba = an; }
            }
            if !best.is_finite() { return 0.0; }
            let w = w0_m * (1.0 + 0.5 * (bs / len_m.max(1.0)).min(1.0));
            let weff = if bb > 0.0 { w * (1.0 + 0.25 * lobe(lobes, ba)) } else { w };
            let m = (t / 0.7).max(mfloor).min(weff.max(mfloor));
            t * smoothstep01((weff - best) / m) * (bs / 200.0).min(1.0)
        }
        VeGeom::Cone { x: cx, y: cy, rb_m, h_m, rc_m, dc_m, .. } => {
            let hh = h_m - a; if hh <= 0.0 { return 0.0; }
            let r = (x - cx).hypot(y - cy); if r >= *rb_m { return 0.0; }
            if r >= *rc_m { hh * ((rb_m - r) / (rb_m - rc_m)).clamp(0.0, 1.0) } else { hh - (dc_m - a).max(0.0) * (1.0 - (r / rc_m).powi(2)) }
        }
        VeGeom::Dome { x: cx, y: cy, r_m, h_m, .. } => { let hh = h_m - a; let r = (x - cx).hypot(y - cy); if hh <= 0.0 || r >= *r_m { 0.0 } else { hh * (1.0 - (r / r_m).powi(4)).max(0.0).sqrt() } }
        VeGeom::Fissure { x0, y0, x1, y1, h_m, w_m } => {
            let hh = h_m - a; if hh <= 0.0 { return 0.0; }
            let (dx, dy) = (x1 - x0, y1 - y0); let l = dx.hypot(dy).max(1.0);
            let al = ((x - x0) * dx + (y - y0) * dy) / l; if al < 0.0 || al > l { return 0.0; }
            let la = ((x - x0) * dy - (y - y0) * dx).abs() / l;
            hh * (-(la / w_m).powi(2)).exp() * smoothstep01(al.min(l - al) / 200.0)
        }
        _ => 0.0,
    }
}

/// Build the history of one volcano (sets `v.ve`); `v` is a constructional / dissected edifice (not a residual).
pub fn build(b: &BuildIn, v: &mut VolcanoSite) -> R<()> {
    let Some(fam) = family_of(&v.kind) else { return Ok(()) };
    if v.denudation_m >= v.h_m || v.h_m <= 0.0 { return Ok(()); }
    let layer = Scope::root(b.seed).layer(b.recipe_digest).map_err(|e| e.to_string())?.chart(0).map_err(|e| e.to_string())?;
    let tile = ((v.record_xy_m[0] * 1000.0 / SPACING_MM as f64).floor().div_euclid(128.0) as i64, (v.record_xy_m[1] * 1000.0 / SPACING_MM as f64).floor().div_euclid(128.0) as i64);
    let sc = layer.tile(8, tile.0, tile.1).map_err(|e| e.to_string())?;
    let sea = b.t0.sea_level_m;
    let (lifetime, t_vis) = fam.windows();
    let span = b.onset.get(&v.id).map(|o| (b.t_my - v.age_my - o).max(0.0));
    let span_my = b.onset.get(&v.id).map(|o| (b.t_my - v.age_my) - o);
    let t_h_ky = match span { Some(sp) => (sp * 1000.0).min(t_vis), None => t_vis }.max(1.0);
    let e_vent = v.e_vent_m_per_my;
    let age0_yr = v.age_my * 1e6;
    let tau = std::f64::consts::TAU;
    // Body geometry.
    let mut g = VeVolcano { flank2: None, flank: None, family: fam.name().into(), span_my, t_h_ky, lifetime_my: lifetime, r_body_m: 0.0, r0_m: 0.0, axes_rad: vec![], stretch: 0.0, axis_w_m: 2500.0, elong: 1.0,
        q_built: 0.0, r_units_m: 0.0, gullies: vec![], rough: None, counts: BTreeMap::new(), dropped_eroded: 0, dropped_cap: 0, units: vec![], wind: None };
    // Lane VH-I ve_wind: the T0 wind at the record (no draw).
    if let Some(w) = b.ve.wind.as_ref() { g.wind = Some(wind_at(w, b.t0, b.chart, v)?); }
    let wind_axis = g.wind.as_ref().and_then(|w| w.chart_azimuth_rad.map(|a| (a, w.rho)));
    let mut kb = Kr::new(&sc, &v.id, 0, b.realization);
    match fam {
        Family::Composite => {
            let cfv = v.cf.as_ref().ok_or("ve composite needs the CF cone")?;
            g.r_body_m = cfv.r_foot_m;
            g.q_built = cfv.q_built;
            // CF's lobes are replaced by the VE flows and its gullies become a cut at gully_age (the CF cone itself stays as
            // published for the solve).
            if b.ve.flank2.is_some() { g.flank2 = Some(Flank2Site { history: vec![], history_rule: "Representative buried lava contacts from native seeded VE history; positions are inferred, not measured. Surviving exposed units retain their actual native geometry and ages.".into() }); }
            else if b.ve.flank.is_some() { g.flank = Some(crate::t1b_coast::FlankSurface::build(v)); }
            else { g.gullies = cfv.gullies.clone(); }
        }
        Family::Shield | Family::Rift => {
            let q = if fam == Family::Shield { 0.05 } else { 0.08 };
            let wb = (v.h_m / q).clamp(10_000.0, 100_000.0);
            g.r_body_m = wb / 2.0;
            g.r0_m = g.r_body_m / 20f64.ln().sqrt();
            g.q_built = v.h_m / wb;
            if fam == Family::Shield {
                let n = if kb.u()? < 0.5 { 2 } else { 3 };
                let a0 = kb.u()? * tau;
                g.axes_rad.push(a0);
                if n == 2 { g.axes_rad.push(a0 + std::f64::consts::PI + kb.range(-0.7, 0.7)?); } else { g.axes_rad.push(a0 + tau / 3.0 + kb.range(-0.35, 0.35)?); g.axes_rad.push(a0 + 2.0 * tau / 3.0 + kb.range(-0.35, 0.35)?); }
                g.stretch = 0.35;
            } else {
                g.axes_rad.push(kb.u()? * std::f64::consts::PI);
                g.elong = 2.5;
            }
        }
        Family::Plateau => { g.r_body_m = 25_000.0; g.q_built = v.h_m / 50_000.0; }
    }
    v.ve = Some(g);
    let rbody = v.ve.as_ref().unwrap().r_body_m;
    // Influence of the body (broad families).
    let elong = v.ve.as_ref().unwrap().elong;
    let r_body_reach = elong * rbody * (1.0 + b.cf.foot_taper_permille as f64 / 1000.0 * 1.8);
    let ix = CurveIx::new(b.curves);
    let vent = v.vent_xy_m;
    let snapshot = v.clone();
    let be = BodyEval { b, v: &snapshot };
    let mut units: Vec<VeUnit> = Vec::new();
    let mut counts: BTreeMap<String, u32> = BTreeMap::new();
    let bump = |k: &str, counts: &mut BTreeMap<String, u32>| { *counts.entry(k.to_string()).or_insert(0) += 1; };
    let age = |u: f64| age0_yr + u * t_h_ky * 1000.0;
    let ero = |a_yr: f64| e_vent * a_yr / 1e6;
    let pt_at = |r: f64, th: f64| [vent[0] + r * th.cos(), vent[1] + r * th.sin()];
    let mut kw = Kr::new(&sc, &v.id, 9, b.realization);
    let v2 = b.ve.ve2.clone().unwrap_or_else(Ve2ParamsV1::off);
    let path_unit = |kind: &str, style: u8, a: f64, vei: u8, start: [f64; 2], len: f64, w0: f64, t: f64, lahar: bool, kw: &mut Kr<'_>| -> R<Option<VeUnit>> {
        let wander = [0.35, kw.range(1500.0, 3000.0)?, kw.u()? * tau];
        let lobes = kw.lobes()?;
        if v2.steer == 1 {
            // Lane VE2: traced after the history is sorted, on the older units, without the wander.
            return Ok(Some(VeUnit { pending: Some((start, len, [0.0, wander[1], wander[2]])), id: 0, kind: kind.into(), style, age_yr: a, vei, size_m: t, erosion_m: ero(a),
                geom: VeGeom::Path { pts: vec![[start[0], start[1], 0.0, LV_BODY]], len_m: len, w0_m: w0, t_m: t, lahar, lobes, captured_at_m: -1.0, end: "pending".into(), ponds: Vec::new() } }));
        }
        let (pts, cap_at, end) = trace(&be, &ix, b.curves, start, len, wander, rbody, sea);
        if pts.len() < 2 { return Ok(None); }
        let l = pts[pts.len() - 1][2];
        Ok(Some(VeUnit { pending: None, id: 0, kind: kind.into(), style, age_yr: a, vei, size_m: t, erosion_m: ero(a), geom: VeGeom::Path { pts, len_m: l, w0_m: w0, t_m: t, lahar, lobes, captured_at_m: if cap_at.is_finite() { cap_at } else { -1.0 }, end, ponds: Vec::new() } }))
    };
    // A wet vent (maar): a network / T1a curve point within 1.5 km whose bed is at most 30 m below the vent surface, or the vent
    // surface within 30 m of sea level.
    let wet = |p: [f64; 2], zs: f64| -> bool {
        if zs <= sea + 30.0 { return true; }
        let (cx0, cy0) = ((p[0] / ix.cell).floor() as i64, (p[1] / ix.cell).floor() as i64);
        for gy in cy0 - 2..=cy0 + 2 { for gx in cx0 - 2..=cx0 + 2 {
            if let Some(l) = ix.grid.get(&(gx, gy)) { for &(ci, pi) in l { let q = b.curves[ci as usize].pts[pi as usize]; if (q[0] - p[0]).hypot(q[1] - p[1]) <= 1500.0 && q[2] >= zs - 30.0 { return true; } } }
        } }
        false
    };
    let surf = |p: [f64; 2]| be.e(p[0], p[1]).unwrap_or_else(|| be.bil(p[0], p[1]));
    // Monogenetic vent: a cinder cone (maybe breached with a flow) or a maar where wet.
    let mut kc = Kr::new(&sc, &v.id, 3, b.realization);
    let vent_unit = |p: [f64; 2], wmed: f64, wlo: f64, whi: f64, a: f64, flow: (f64, f64, f64, f64, f64, f64), units: &mut Vec<VeUnit>, counts: &mut BTreeMap<String, u32>, kc: &mut Kr<'_>, kw: &mut Kr<'_>| -> R<()> {
        let zs = surf(p);
        let um = kc.u()?;
        if wet(p, zs) && um < 0.7 {
            let dm = kc.logn(700.0, 0.4, 300.0, 1500.0)?;
            units.push(VeUnit { pending: None, id: 0, kind: "maar".into(), style: ST_FLOOR, age_yr: a, vei: 2, size_m: 0.15 * dm, erosion_m: ero(a), geom: VeGeom::Maar { x: p[0], y: p[1], r_m: dm / 2.0, d_m: 0.15 * dm, rim_m: 0.04 * dm, ring_b_m: 0.25 * dm } });
            bump("maar", counts);
            return Ok(());
        }
        let wco = kc.logn(wmed, 0.4, wlo, whi)?;
        let breached = kc.u()? < 0.5;
        // Breach toward the downslope side (body gradient, else outward).
        let gdir = be.grad(p[0], p[1]).map(|g| (-g[1]).atan2(-g[0])).unwrap_or((p[1] - vent[1]).atan2(p[0] - vent[0]));
        units.push(VeUnit { pending: None, id: 0, kind: "cinder_cone".into(), style: ST_CINDER, age_yr: a, vei: 2, size_m: 0.18 * wco, erosion_m: ero(a),
            geom: VeGeom::Cone { x: p[0], y: p[1], rb_m: wco / 2.0, h_m: 0.18 * wco, rc_m: 0.2 * wco, dc_m: 0.3 * 0.18 * wco, breach_rad: if breached { Some(gdir) } else { None } } });
        bump("cinder_cone", counts);
        if breached {
            let (lmin, lmax, wl, wh, tl, th) = flow;
            let len = kc.range(lmin, lmax)?;
            let st = [p[0] + 0.5 * wco * gdir.cos(), p[1] + 0.5 * wco * gdir.sin()];
            let (w0, t) = (kc.range(wl, wh)?, kc.range(tl, th)?);
            if let Some(u) = path_unit("cone_flow", ST_LAVA, a - 1.0, 0, st, len, w0, t, false, kw)? { units.push(u); bump("cone_flow", counts); }
        }
        Ok(())
    };
    let (rc, h) = (v.crater_radius_m, v.h_m);
    match fam {
        Family::Composite => {
            let rf = rbody;
            let mut kf = Kr::new(&sc, &v.id, 1, b.realization);
            let mut ke = Kr::new(&sc, &v.id, 2, b.realization);
            let mut kd = Kr::new(&sc, &v.id, 4, b.realization);
            // Lava flows.
            let nf = kf.poisson(0.5 * t_h_ky)?;
            for _ in 0..nf {
                let a = age(kf.u()?);
                let th = kf.u()? * tau;
                let start = if kf.u()? < 0.6 { pt_at(rc + 50.0, th) } else { pt_at(kf.range(0.15, 0.5)? * rf, th) };
                let len = kf.logn(3500.0, 0.5, 1500.0, 8000.0)?;
                let (w0, t) = (kf.range(100.0, 300.0)?, kf.range(15.0, 50.0)?);
                if let Some(u) = path_unit("lava_flow", ST_LAVA, a, 0, start, len, w0, t, false, &mut kw)? { units.push(u); bump("lava_flow", &mut counts); }
            }
            // Explosive eruptions: VEI, aprons, lahars, caldera.
            let large_old = span.map(|s| s >= 1.0).unwrap_or(false) && h >= 1500.0;
            let ne = ke.poisson(0.25 * t_h_ky)?;
            let mut caldera_done = false;
            for _ in 0..ne {
                let a = age(ke.u()?);
                let uv = ke.u()?.max(1e-12);
                let vei = (3.0 + (2.0 * -uv.log10()).floor()).min(6.0) as u8;
                bump(&format!("explosive_vei{vei}"), &mut counts);
                if vei >= 4 {
                    let (s_f, cap_max, reach, h0, h1, bb) = match vei { 4 => (0.06, 15.0, 1.2, 45.0, 90.0, 2000.0), 5 => (0.04, 40.0, 1.5, 90.0, 180.0, 4000.0), _ => (0.03, 100.0, 2.0, 180.0, 180.0, 8000.0) };
                    // Lane VH-I ve_wind: the T0 downwind azimuth plus a wrapped-Cauchy deviation from the same draw.
                    let dir = match wind_axis { Some((az, rho)) => az + wind_deviation(rho, ke.u()?), None => ke.u()? * tau };
                    let half = ke.range(h0, h1)?.to_radians();
                    // Volume budget: bulk tephra V = 10^(VEI + 4.5) m3 (class geometric mean); 30% proximal density currents
                    // over the apron annulus (thickness cap 6 x the mean: valleys take most of it), 20% proximal fall.
                    let vbulk = 10f64.powf(vei as f64 + 4.5);
                    let frac = (2.0 * half).min(tau) / tau;
                    let area = frac * std::f64::consts::PI * ((reach * rf).powi(2) - (0.85 * rf).powi(2)).max(1.0);
                    let t_cap = (2.0 * 0.3 * vbulk / area).min(cap_max);
                    // Fall: elliptical isopachs elongated downwind (x2.5 downwind, x0.6 upwind of the vent): area factor 1.55.
                    let t0 = 0.2 * vbulk / (1.55 * tau * bb * bb);
                    // Apex: where the body slope falls to 0.15 (at most 0.85 R_f).
                    let mut ra = 0.85 * rf;
                    let mut rr = rc + 100.0;
                    while rr < 0.85 * rf {
                        let p0 = pt_at(rr, dir);
                        if let Some(gr) = be.grad(p0[0], p0[1]) { if gr[0].hypot(gr[1]) <= 0.15 { ra = rr; break; } }
                        rr += 100.0;
                    }
                    let za: Vec<f64> = (0..72).map(|k| { let p0 = pt_at(ra, k as f64 * tau / 72.0); be.e_any(p0[0], p0[1]).unwrap_or(v.base_m) }).collect();
                    units.push(VeUnit { pending: None, id: 0, kind: "pyroclastic_apron".into(), style: ST_APRON, age_yr: a, vei, size_m: t_cap, erosion_m: ero(a),
                        geom: VeGeom::Apron { x: vent[0], y: vent[1], dir_rad: dir, half_rad: half, r_a_m: ra, r_max_m: reach * rf, s_f, t_cap_m: t_cap, t0_m: t0, b_m: bb, za } });
                    bump("pyroclastic_apron", &mut counts);
                }
                if vei >= 6 && large_old && !caldera_done {
                    caldera_done = true;
                    let rcal = ke.range(2500.0, 4000.0)?.min(0.45 * rf);
                    let zr: f64 = (0..36).map(|k| { let p0 = pt_at(rcal, k as f64 * tau / 36.0); be.e_any(p0[0], p0[1]).unwrap_or(v.base_m) }).fold(f64::INFINITY, f64::min);
                    let depth = ke.range(0.2, 0.3)? * h;
                    let wall = 0.8;
                    let floor = zr - depth;
                    units.push(VeUnit { pending: None, id: 0, kind: "caldera".into(), style: ST_WALL, age_yr: a - 1.0, vei, size_m: depth, erosion_m: ero(a), geom: VeGeom::Caldera { x: vent[0], y: vent[1], r_m: rcal, floor_m: floor, r_floor_m: (rcal - depth / wall).max(0.0), wall } });
                    bump("caldera", &mut counts);
                }
                let nl = match vei { 3 => 1, 4 => 1 + (ke.u()? < 0.5) as u32, _ => 2 + (ke.u()? < 0.5) as u32 };
                for _ in 0..nl {
                    let th = ke.u()? * tau;
                    let start = pt_at(ke.range(0.55, 0.75)? * rf, th);
                    let len = ke.logn(10_000.0, 0.4, 5000.0, 25_000.0)?;
                    let (w0, t) = (ke.range(250.0, 600.0)?, ke.range(5.0, 20.0)?);
                    if let Some(u) = path_unit("lahar", ST_LAHAR, a - 10.0, vei, start, len, w0, t, true, &mut kw)? { units.push(u); bump("lahar", &mut counts); }
                }
            }
            // Sector collapse / lateral blast amphitheatre.
            if kd.u()? < 0.15 {
                let a = age(kd.u()?);
                let dir = v.breach_dir_rad + kd.range(-0.5, 0.5)?;
                let depth = kd.range(0.15, 0.25)? * h;
                let ztop = be.e_any(vent[0], vent[1]).unwrap_or(v.summit_built());
                units.push(VeUnit { pending: None, id: 0, kind: "amphitheatre".into(), style: ST_WALL, age_yr: a, vei: 5, size_m: depth, erosion_m: ero(a),
                    geom: VeGeom::Breach { x: vent[0], y: vent[1], dir_rad: dir, half_w_m: kd.range(800.0, 1200.0)?, floor0_m: ztop - depth, slope: 0.15, len_m: 0.5 * rf, back_m: 300.0, wall: 0.8 } });
                bump("amphitheatre", &mut counts);
            }
            // Domes.
            let rate = if v.kind == "continental_arc" { 0.05 } else { 0.025 };
            let nd = kd.poisson(rate * t_h_ky)?;
            for _ in 0..nd {
                let a = age(kd.u()?);
                let p = pt_at(kd.range(0.0, 0.35)? * rf, kd.u()? * tau);
                let mut dd = kd.logn(800.0, 0.35, 400.0, 1500.0)?;
                let hd = kd.range(0.2, 0.35)? * dd;
                let (mut p, mut top) = (p, None);
                if v2.dome == 1 {
                    // Lane VE2 dome cap: summit share 0.5 inside the crater; tops at most the rim near the summit.
                    let us = kd.u()?;
                    if us < 0.5 { p = pt_at(0.5 * rc * us / 0.5, kd.u()? * tau); dd = dd.min(2.0 * rc); }
                    else { let r = (p[0] - vent[0]).hypot(p[1] - vent[1]); if r < 0.1 * rf { p = pt_at(0.1 * rf, (p[1] - vent[1]).atan2(p[0] - vent[0])); } }
                    if (p[0] - vent[0]).hypot(p[1] - vent[1]) < 2.0 * rc {
                        let rim = (0..36).map(|k| { let q = pt_at(rc, k as f64 * tau / 36.0); be.e_any(q[0], q[1]).unwrap_or(v.summit_built()) }).fold(f64::INFINITY, f64::min);
                        top = Some(rim);
                    }
                }
                units.push(VeUnit { pending: None, id: 0, kind: "dome".into(), style: ST_DOME, age_yr: a, vei: 0, size_m: hd, erosion_m: ero(a), geom: VeGeom::Dome { x: p[0], y: p[1], r_m: dd / 2.0, h_m: hd, top_m: top } });
                bump("dome", &mut counts);
            }
            // Flank vents.
            let nv = kd.poisson(0.05 * t_h_ky)?;
            for _ in 0..nv {
                let a = age(kd.u()?);
                let p = pt_at(kd.range(0.3, 0.9)? * rf, kd.u()? * tau);
                vent_unit(p, 900.0, 300.0, 2000.0, a, (1000.0, 5000.0, 80.0, 200.0, 5.0, 15.0), &mut units, &mut counts, &mut kc, &mut kw)?;
            }
            // CF gullies as a cut at gully_age.
            if !v.ve.as_ref().unwrap().gullies.is_empty() {
                let n = v.ve.as_ref().unwrap().gullies.len() as u32;
                units.push(VeUnit { pending: None, id: 0, kind: "gully_incision".into(), style: ST_NONE, age_yr: age0_yr + b.ve.gully_age_yr as f64, vei: 0, size_m: 0.0, erosion_m: 0.0, geom: VeGeom::Gullies { n } });
            }
        }
        Family::Shield => {
            let rb = rbody;
            let axes = v.ve.as_ref().unwrap().axes_rad.clone();
            let mut kf = Kr::new(&sc, &v.id, 1, b.realization);
            let mut kr = Kr::new(&sc, &v.id, 2, b.realization);
            // Summit caldera.
            let rcal = kr.range(1500.0, 2500.0)?;
            let depth = kr.range(100.0, 200.0)?;
            let zr: f64 = (0..36).map(|k| { let p0 = pt_at(rcal, k as f64 * tau / 36.0); be.e_any(p0[0], p0[1]).unwrap_or(v.base_m) }).fold(f64::INFINITY, f64::min);
            let ac = age(kr.u()?);
            units.push(VeUnit { pending: None, id: 0, kind: "caldera".into(), style: ST_WALL, age_yr: ac, vei: 0, size_m: depth, erosion_m: ero(ac), geom: VeGeom::Caldera { x: vent[0], y: vent[1], r_m: rcal, floor_m: zr - depth, r_floor_m: (rcal - depth).max(0.0), wall: 1.0 } });
            bump("caldera", &mut counts);
            let rift_pt = |k: &mut Kr, a0: f64, a1: f64| -> R<[f64; 2]> {
                let ax = axes[(k.u()? * axes.len() as f64).floor().min(axes.len() as f64 - 1.0) as usize];
                let al = k.range(a0, a1)? * rb;
                let la = k.range(-300.0, 300.0)?;
                Ok([vent[0] + al * ax.cos() - la * ax.sin(), vent[1] + al * ax.sin() + la * ax.cos()])
            };
            // Flows.
            let nf = kf.poisson(4.0 * t_h_ky)?.min(40);
            for _ in 0..nf {
                let a = age(kf.u()?);
                let start = if kf.u()? < 0.4 { pt_at(rcal + 50.0, kf.u()? * tau) } else { rift_pt(&mut kf, 0.1, 0.7)? };
                let len = kf.logn(12_000.0, 0.6, 3000.0, 50_000.0)?;
                let (w0, t) = (kf.range(200.0, 800.0)?, kf.range(3.0, 12.0)?);
                if let Some(u) = path_unit("shield_flow", ST_SHIELD, a, 0, start, len, w0, t, false, &mut kw)? { units.push(u); bump("shield_flow", &mut counts); }
            }
            // Rift fissures, cones, pit craters.
            let nfi = kr.poisson(0.6 * t_h_ky)?.min(8);
            for _ in 0..nfi {
                let a = age(kr.u()?);
                let ax = axes[(kr.u()? * axes.len() as f64).floor().min(axes.len() as f64 - 1.0) as usize];
                let al = kr.range(0.05, 0.6)? * rb;
                let l = kr.range(500.0, 3000.0)?;
                let p0 = [vent[0] + al * ax.cos(), vent[1] + al * ax.sin()];
                let p1 = [p0[0] + l * ax.cos(), p0[1] + l * ax.sin()];
                let hf = kr.range(5.0, 15.0)?;
                units.push(VeUnit { pending: None, id: 0, kind: "fissure".into(), style: ST_FISSURE, age_yr: a, vei: 1, size_m: hf, erosion_m: ero(a), geom: VeGeom::Fissure { x0: p0[0], y0: p0[1], x1: p1[0], y1: p1[1], h_m: hf, w_m: 60.0 } });
                bump("fissure", &mut counts);
            }
            let nc = kr.poisson(0.5 * t_h_ky)?.min(8);
            for _ in 0..nc {
                let a = age(kr.u()?);
                let p = rift_pt(&mut kr, 0.1, 0.8)?;
                vent_unit(p, 600.0, 200.0, 1500.0, a, (2000.0, 8000.0, 150.0, 400.0, 3.0, 10.0), &mut units, &mut counts, &mut kc, &mut kw)?;
            }
            let np = 2 + (kr.u()? * 3.0).floor() as u32;
            for _ in 0..np {
                let a = age(kr.u()?);
                let p = rift_pt(&mut kr, 0.05, 0.35)?;
                let dp = kr.range(300.0, 1000.0)?;
                let d = kr.range(50.0, 150.0)?;
                units.push(VeUnit { pending: None, id: 0, kind: "pit_crater".into(), style: ST_FLOOR, age_yr: a, vei: 0, size_m: d, erosion_m: ero(a), geom: VeGeom::Pit { x: p[0], y: p[1], r_m: dp / 2.0, d_m: d } });
                bump("pit_crater", &mut counts);
            }
        }
        Family::Rift => {
            let ax = v.ve.as_ref().unwrap().axes_rad[0];
            let r0 = v.ve.as_ref().unwrap().r0_m;
            let mut kr = Kr::new(&sc, &v.id, 1, b.realization);
            let rows = 3 + (kr.u()? * 4.0).floor() as u32;
            for _ in 0..rows {
                let a = age(kr.u()?);
                let off = kr.range(-1.0, 1.0)? * r0;
                let al = kr.range(-1.0, 1.0)? * 2.5 * r0;
                let l = kr.range(2000.0, 8000.0)?;
                let c = [vent[0] + al * ax.cos() - off * ax.sin(), vent[1] + al * ax.sin() + off * ax.cos()];
                let (p0, p1) = ([c[0] - 0.5 * l * ax.cos(), c[1] - 0.5 * l * ax.sin()], [c[0] + 0.5 * l * ax.cos(), c[1] + 0.5 * l * ax.sin()]);
                let hf = kr.range(5.0, 15.0)?;
                units.push(VeUnit { pending: None, id: 0, kind: "fissure".into(), style: ST_FISSURE, age_yr: a, vei: 1, size_m: hf, erosion_m: ero(a), geom: VeGeom::Fissure { x0: p0[0], y0: p0[1], x1: p1[0], y1: p1[1], h_m: hf, w_m: 60.0 } });
                bump("fissure", &mut counts);
                let ncone = 3 + (kr.u()? * 8.0).floor() as u32;
                for _ in 0..ncone {
                    let f = kr.u()?;
                    let p = [p0[0] + f * (p1[0] - p0[0]), p0[1] + f * (p1[1] - p0[1])];
                    vent_unit(p, 600.0, 200.0, 1500.0, a - kr.u()? * 500.0, (2000.0, 8000.0, 150.0, 400.0, 3.0, 10.0), &mut units, &mut counts, &mut kc, &mut kw)?;
                }
                let nfl = 2 + (kr.u()? * 3.0).floor() as u32;
                for _ in 0..nfl {
                    let f = kr.u()?;
                    let p = [p0[0] + f * (p1[0] - p0[0]), p0[1] + f * (p1[1] - p0[1])];
                    let len = kr.logn(6000.0, 0.5, 2000.0, 20_000.0)?;
                    let (w0, t) = (kr.range(200.0, 600.0)?, kr.range(3.0, 10.0)?);
                    if let Some(u) = path_unit("rift_flow", ST_SHIELD, a - 1.0, 0, p, len, w0, t, false, &mut kw)? { units.push(u); bump("rift_flow", &mut counts); }
                }
            }
        }
        Family::Plateau => {
            let mut kp = Kr::new(&sc, &v.id, 1, b.realization);
            let n = ((h / 40.0).round() as u32).clamp(3, 12);
            for k in 0..n {
                let a = age0_yr + (n - k) as f64 / n as f64 * t_h_ky * 1000.0;
                let t = h / n as f64 * kp.range(0.7, 1.3)?;
                let r = 25_000.0 * (1.0 - 0.5 * k as f64 / n as f64);
                units.push(VeUnit { pending: None, id: 0, kind: "plateau_flow".into(), style: ST_PLATEAU, age_yr: a, vei: 0, size_m: t, erosion_m: ero(a), geom: VeGeom::Sheet { x: vent[0], y: vent[1], r_m: r, t_m: t, lobes: kp.lobes()? } });
                bump("plateau_flow", &mut counts);
            }
        }
    }
    // Lane VE2 fresh-lava erosion: coherent lava units erode by weathering only until the infiltration lag has passed.
    if v2.erode == 1 {
        let ew = v2.weather_mm_per_my as f64 / 1000.0; // m/My
        for u in units.iter_mut() {
            let basaltic = fam != Family::Composite || matches!(u.kind.as_str(), "cinder_cone" | "cone_flow");
            let lava = matches!(u.kind.as_str(), "lava_flow" | "cone_flow" | "shield_flow" | "rift_flow" | "plateau_flow" | "dome" | "cinder_cone" | "fissure");
            if !lava { continue; }
            let tau = if basaltic { v2.lag_basalt_yr } else { v2.lag_andesite_yr } as f64;
            u.erosion_m = ew * u.age_yr.min(tau) / 1e6 + e_vent * (u.age_yr - tau).max(0.0) / 1e6;
        }
    }
    let mut buried_history = Vec::new();
    if b.ve.flank2.is_some() && fam == Family::Composite {
        for u in &units {
            if let VeGeom::Path { lahar: false, .. } = u.geom {
                let mut old = u.clone();
                if let Some((start,len,wander)) = old.pending.take() {
                    let (pts,cap_at,end) = trace_g(&be,&ix,b.curves,start,len,wander,rbody,sea,None);
                    if let VeGeom::Path { pts: pp, len_m, captured_at_m, end: en, .. } = &mut old.geom {
                        *len_m=pts.last().map(|p|p[2]).unwrap_or(0.0);*pp=pts;*captured_at_m=if cap_at.is_finite(){cap_at}else{-1.0};*en=end;
                    }
                }
                buried_history.push(old);
            }
        }
    }
    // Erosion and the cap; oldest first; ids.
    let tmap = b.ve.t_map_mm as f64 / 1000.0;
    let mut kept: Vec<VeUnit> = Vec::new();
    let mut dropped_eroded = 0;
    for u in units {
        let gone = match u.geom { VeGeom::Gullies { .. } => false, _ => u.size_m - u.erosion_m < tmap };
        if gone { dropped_eroded += 1; } else { kept.push(u); }
    }
    kept.sort_by(|a, b| b.age_yr.total_cmp(&a.age_yr).then(a.kind.cmp(&b.kind)));
    let mut dropped_cap = 0;
    if kept.len() > b.ve.max_units as usize {
        // Keep the youngest (the oldest are the most buried); the gully cut always stays.
        let excess = kept.len() - b.ve.max_units as usize;
        let mut removed = 0;
        kept.retain(|u| { if removed < excess && !matches!(u.geom, VeGeom::Gullies { .. }) { removed += 1; false } else { true } });
        dropped_cap = removed as u32;
    }
    // Lane VE2 steering: trace the pending paths oldest first on the body plus the units already built.
    // Lane VE3: the lobe roughness from this volcano's own flows, and the body with it for the steered tracing.
    let v3 = b.ve.ve3.clone().unwrap_or_else(Ve3ParamsV1::off);
    let ve3_on = b.ve.ve3.is_some() && v2.steer == 1;
    let v4 = b.ve.ve4.clone().unwrap_or_else(Ve4ParamsV1::off);
    let v5 = b.ve.ve5.clone().unwrap_or_else(Ve5ParamsV1::off);
    let mut rough: Option<RoughSpec> = None;
    if ve3_on && v3.rough == 1 {
        let mut ws: Vec<f64> = Vec::new();
        let mut ts: Vec<f64> = Vec::new();
        for u in &kept { if let VeGeom::Path { w0_m, t_m, lahar: false, .. } = &u.geom { ws.push(*w0_m); ts.push((t_m - u.erosion_m).max(0.0)); } }
        if !ws.is_empty() {
            ws.sort_by(|a, b| a.total_cmp(b)); ts.sort_by(|a, b| a.total_cmp(b));
            let (wm, tm) = (ws[ws.len() / 2], ts[ts.len() / 2]);
            let cell = 2.0 * wm;
            let elong = v.ve.as_ref().map(|g| g.elong).unwrap_or(1.0);
            let ext = rbody * elong * 1.3;
            let (c0x, c1x) = (((vent[0] - ext) / cell).floor() as i64, ((vent[0] + ext) / cell).ceil() as i64);
            let (c0y, c1y) = (((vent[1] - ext) / cell).floor() as i64, ((vent[1] + ext) / cell).ceil() as i64);
            let mut lobes: Vec<[f64; 6]> = Vec::new();
            let mut grid: BTreeMap<(i64, i64), u32> = BTreeMap::new();
            let mut kr = Kr::new(&sc, &v.id, 20, b.realization);
            for cy in c0y..=c1y { for cx in c0x..=c1x {
                // Keyed per cell (counter = cell index in the volcano's grid, 6 draws).
                kr.k = (((cy - c0y) * (c1x - c0x + 1) + (cx - c0x)) as u64) * 6;
                let (ux, uy) = (kr.u()?, kr.u()?);
                let (x, y) = ((cx as f64 + ux) * cell, (cy as f64 + uy) * cell);
                if be.e(x, y).is_none() { continue; }
                let hw = wm * kr.range(0.6, 1.2)?;
                let amp = tm * kr.range(0.3, 1.0)?;
                let axis = be.grad(x, y).map(|g| g[1].atan2(g[0])).unwrap_or(0.0);
                grid.insert((cx, cy), lobes.len() as u32);
                lobes.push([x, y, amp, hw, 3.0 * hw, axis]);
            } }
            rough = Some(RoughSpec { w_med_m: wm, t_med_m: tm, cell_m: cell, lobes_n: lobes.len() as u32, lobes, grid, reach: 2 });
            // Lane VE4: lobes by slope (Jeffreys) and age (eroded), occupancy 0.8, summit taper.
            if v4.lobes == 1 {
                let ncx = c1x - c0x + 1;
                let mut kr4 = Kr::new(&sc, &v.id, 21, b.realization);
                let mut cand: Vec<(i64, i64, [f64; 2], f64, [f64; 4])> = Vec::new();
                for cy in c0y..=c1y { for cx in c0x..=c1x {
                    kr4.k = (((cy - c0y) * ncx + (cx - c0x)) as u64) * 8;
                    let occ = kr4.u()?;
                    let (ux, uy) = (kr4.u()?, kr4.u()?);
                    let (x, y) = ((cx as f64 + ux) * cell, (cy as f64 + uy) * cell);
                    let rest = [kr4.u()?, kr4.u()?, kr4.u()?, kr4.u()?];
                    if occ >= 0.8 || be.e(x, y).is_none() { continue; }
                    let slope = be.grad(x, y).map(|g| g[0].hypot(g[1])).unwrap_or(0.0).max(1e-3);
                    cand.push((cx, cy, [x, y], slope, rest));
                } }
                let mut sl: Vec<f64> = cand.iter().map(|c| c.3).collect();
                sl.sort_by(|a, b| a.total_cmp(b));
                let smed = if sl.is_empty() { 0.1 } else { sl[sl.len() / 2] };
                let (t_lo, t_hi) = (t_h_ky * 1000.0, (lifetime * 1e6).max(t_h_ky * 1000.0 * 1.01));
                let ew = v2.weather_mm_per_my as f64 / 1000.0;
                let tau = if fam == Family::Composite { v2.lag_andesite_yr } else { v2.lag_basalt_yr } as f64;
                let rc = v.crater_radius_m.max(1.0);
                let mut lobes4: Vec<[f64; 6]> = Vec::new();
                let mut grid4: BTreeMap<(i64, i64), u32> = BTreeMap::new();
                for (cx, cy, p, slope, r4) in cand {
                    let f = (smed / slope).powf(1.0 / 3.0).clamp(0.6, 1.6);
                    let age_l = if v5.exposure == 1 {
                        // Lane VE5: exposure age (burial by younger units hides older lobes).
                        let rate = (v.h_m / 3.0) / (lifetime * 1e6);
                        let tau_e = tm / rate.max(1e-9);
                        age0_yr + t_lo + tau_e * (-(1.0 - r4[0]).max(1e-12).ln())
                    } else { age0_yr + t_lo * (t_hi / t_lo).powf(r4[0]) };
                    let er = if v2.erode == 1 { ew * age_l.min(tau) / 1e6 + e_vent * (age_l - tau).max(0.0) / 1e6 } else { e_vent * age_l / 1e6 };
                    let rr = (p[0] - vent[0]).hypot(p[1] - vent[1]);
                    let taper = smoothstep01((rr - rc) / (2.0 * rc));
                    let amp = ((tm * (0.3 + 0.7 * r4[1]) * f) - er).max(0.0) * taper;
                    if amp <= 0.05 { continue; }
                    let hw = wm * (0.6 + 0.6 * r4[2]) * f;
                    let axis = be.grad(p[0], p[1]).map(|g| g[1].atan2(g[0])).unwrap_or(0.0);
                    grid4.insert((cx, cy), lobes4.len() as u32);
                    lobes4.push([p[0], p[1], amp, hw, 3.0 * hw, axis]);
                }
                rough = Some(RoughSpec { w_med_m: wm, t_med_m: tm, cell_m: cell, lobes_n: lobes4.len() as u32, lobes: lobes4, grid: grid4, reach: 3 });
            }
        }
    }
    let snap3 = { let mut s3 = snapshot.clone(); if let Some(g) = s3.ve.as_mut() { g.rough = rough.clone(); } s3 };
    let be3: Option<BodyEval> = if ve3_on { Some(BodyEval { b, v: &snap3 }) } else { None };
    if v2.steer == 1 {
        let gully_age = b.ve.gully_age_yr as f64 + age0_yr;
        let mut built: Vec<VeUnit> = Vec::new();
        for mut u in std::mem::take(&mut kept) {
            if let Some((start, len, wander)) = u.pending.take() {
                let older = &built;
                // REVISION 1 (VE2): the CF gully set is the cone's drainage at every age (a sub-record closure), so every
                // composite flow sees it; the gully_age cut still marks when the published gullies incise the units.
                let _ = gully_age;
                let young = true;
                let gref = &snapshot;
                let extra = |x: f64, y: f64| -> f64 {
                    let mut h = 0.0;
                    for o in older.iter() { h += unit_drape(o, x, y); }
                    if young && fam == Family::Composite { if let Some(g) = gref.ve.as_ref() { h -= gully_cut(g, gref, b.cf, x, y); } }
                    h
                };
                let (pts, cap_at, end, pds) = match &be3 {
                    None => { let (a1, a2, a3) = trace_g(&be, &ix, b.curves, start, len, wander, rbody, sea, Some(&extra)); (a1, a2, a3, Vec::new()) }
                    Some(bb) => {
                        let mf = if v3.aa == 1 { 48.0 } else { 32.0 };
                        let extra3 = |x: f64, y: f64| -> f64 {
                            let mut h = 0.0;
                            for o in older.iter() { h += unit_drape_m(o, x, y, mf); }
                            if young && fam == Family::Composite { if let Some(g) = gref.ve.as_ref() { h -= gully_cut(g, gref, b.cf, x, y); } }
                            h
                        };
                        let (w0, t) = match &u.geom { VeGeom::Path { w0_m, t_m, .. } => (*w0_m, (t_m - u.erosion_m).max(0.1)), _ => (100.0, 10.0) };
                        let is_lahar = matches!(u.geom, VeGeom::Path { lahar: true, .. });
                        let pc = if (v3.pond == 1 && !is_lahar) || (v4.lahar_pond == 1 && is_lahar) { Some(PondCfg { w: w0, t, radius: v3.pond_radius_mm as f64 / 1000.0, max: v3.ponds_max, edge_spill: b.ve.ve5.is_some() }) } else { None };
                        trace_g3(bb, &ix, b.curves, start, len, wander, rbody, sea, Some(&extra3), pc)
                    }
                };
                if pts.len() < 2 { continue; }
                let l = pts[pts.len() - 1][2];
                if let VeGeom::Path { pts: pp, len_m, captured_at_m, end: en, ponds, .. } = &mut u.geom { *pp = pts; *len_m = l; *captured_at_m = if cap_at.is_finite() { cap_at } else { -1.0 }; *en = end; *ponds = pds; }
            }
            built.push(u);
        }
        kept = built;
    }
    let mut reach: f64 = rbody.max(r_body_reach);
    for (k, u) in kept.iter_mut().enumerate() {
        u.id = (k + 1) as u16;
        let rr = match &u.geom {
            VeGeom::Path { pts, w0_m, .. } => pts.iter().map(|p| (p[0] - vent[0]).hypot(p[1] - vent[1])).fold(0.0, f64::max) + 2.0 * w0_m + 200.0,
            VeGeom::Cone { x, y, rb_m, .. } => (x - vent[0]).hypot(y - vent[1]) + rb_m,
            VeGeom::Dome { x, y, r_m, .. } | VeGeom::Pit { x, y, r_m, .. } => (x - vent[0]).hypot(y - vent[1]) + r_m,
            VeGeom::Maar { x, y, r_m, ring_b_m, .. } => (x - vent[0]).hypot(y - vent[1]) + r_m + 4.0 * ring_b_m,
            VeGeom::Fissure { x0, y0, x1, y1, w_m, .. } => (x0 - vent[0]).hypot(y0 - vent[1]).max((x1 - vent[0]).hypot(y1 - vent[1])) + 3.0 * w_m,
            VeGeom::Caldera { r_m, .. } => *r_m,
            VeGeom::Breach { len_m, .. } => len_m + 2000.0,
            VeGeom::Apron { r_max_m, .. } => *r_max_m,
            VeGeom::Sheet { r_m, .. } => r_m * 1.2,
            VeGeom::Gullies { .. } => rbody,
        };
        reach = reach.max(rr);
    }
    let g = v.ve.as_mut().unwrap();
    if let Some(f) = g.flank2.as_mut() { f.history = buried_history; }
    g.rough = rough;
    g.units = kept;
    g.counts = counts;
    g.dropped_eroded = dropped_eroded;
    g.dropped_cap = dropped_cap;
    g.r_units_m = reach;
    Ok(())
}

/// Shield / rift / plateau clearance from T1a channels (vent placement, lane VE): the vent moves until the body height at the
/// nearest T1a curve is at most 20% of H (Gaussian: r = R0 sqrt(ln 5)), at least lane CO's clearance.
pub fn broad_clearance(kind: &str, h: f64) -> Option<f64> {
    let fam = family_of(kind)?;
    if !fam.broad() { return None; }
    let q = match fam { Family::Shield => 0.05, Family::Rift => 0.08, _ => return Some(0.0) };
    let wb = (h / q).clamp(10_000.0, 100_000.0);
    let r0 = (wb / 2.0) / 20f64.ln().sqrt();
    Some(r0 * 5f64.ln().sqrt())
}

/// Build every volcano's history (only those whose family is known and whose edifice is not a residual).
pub fn build_all(b: &BuildIn, vs: &mut [VolcanoSite]) -> R<()> {
    for v in vs.iter_mut() { build(b, v)?; }
    Ok(())
}

// ------------------------------------------------------------------ application on a lattice

#[derive(Clone, Debug, Default, Serialize)]
pub struct VeCoastStats {
    pub volcanoes_with_units: u32,
    pub units_in_lattice: u32,
    pub posts_by_style: BTreeMap<String, u32>,
    pub unit_raised_posts: u32,
    pub unit_cut_posts: u32,
    pub corridor_capped_posts: u32,
    pub buried_channel_posts: u32,
    /// Declared lakes: closed depressions of craters, calderas, maars and pits the final fill raised (posts, max raise m).
    pub declared_lake_posts: u32,
    pub declared_lake_raise_max_m: f64,
    pub other_fill_raised_posts: u32,
}

fn style_name(s: u8) -> &'static str {
    match s { ST_NONE => "0_none", ST_LAVA => "1_lava_flow", ST_FRONT => "2_lobe_front", ST_CINDER => "3_cinder_cone", ST_APRON => "4_pyroclastic_apron", ST_LAHAR => "5_lahar", ST_DOME => "6_dome",
        ST_WALL => "7_crater_wall", ST_SHIELD => "8_shield_flow", ST_FISSURE => "9_fissure", ST_FLOOR => "10_crater_floor", ST_TUFF => "11_tuff_ring", ST_BODY => "12_composite_body", ST_PLATEAU => "13_plateau_flow", _ => "14_neck" }
}

/// Posts of the lattice within `r` of (x, y) (index ranges clipped to the lattice).
fn post_range(lat: &EvalLat, x0: f64, y0: f64, x1: f64, y1: f64) -> Option<(i64, i64, i64, i64)> {
    let dx = lat.dx();
    let (gx0, gy0) = (lat.gx0, lat.gy0);
    let (gx1, gy1) = (gx0 + lat.n as i64 - 1, gy0 + lat.n as i64 - 1);
    let (a0, a1) = (((x0 / dx).floor() as i64).max(gx0), ((x1 / dx).ceil() as i64).min(gx1));
    let (b0, b1) = (((y0 / dx).floor() as i64).max(gy0), ((y1 / dx).ceil() as i64).min(gy1));
    if a0 > a1 || b0 > b1 { return None; }
    Some((a0, a1, b0, b1))
}

/// Apply every unit (oldest first) of the volcanoes `vsi` after the body loop. `z_pre`: the surface before any edifice;
/// `z`: after the body (in place). Writes style / unit per post and updates over / volcano index / fixed.
#[allow(clippy::too_many_arguments)]
pub fn apply_units(lat: &EvalLat, vp: &VolcanoParamsV1, cx: &CoastCtx, vsi: &[usize], z_pre: &[f64], z: &mut [f64], fixed: &mut [u8], t1a_post: &[bool],
    over: &mut [u8], vid: &mut [u8], style: &mut [u8], unit: &mut [u16], sea: f64, st: &mut VeCoastStats) {
    apply_units_f(lat, vp, cx, vsi, z_pre, z, fixed, t1a_post, over, vid, style, unit, sea, st, &mut [])
}

/// As `apply_units`, with the lane VE4 front field (`front` empty = off; NaN = none).
#[allow(clippy::too_many_arguments)]
pub fn apply_units_f(lat: &EvalLat, vp: &VolcanoParamsV1, cx: &CoastCtx, vsi: &[usize], z_pre: &[f64], z: &mut [f64], fixed: &mut [u8], t1a_post: &[bool],
    over: &mut [u8], vid: &mut [u8], style: &mut [u8], unit: &mut [u16], sea: f64, st: &mut VeCoastStats, front: &mut [f64]) {
    let (nn, dx) = (lat.nodes(), lat.dx());
    let Some(cf) = cx.cf.as_ref() else { return };
    let Some(ve) = cx.ve.as_ref() else { return };
    let tmap = ve.t_map_mm as f64 / 1000.0;
    let v2 = ve.ve2.clone().unwrap_or_else(Ve2ParamsV1::off);
    let v3 = ve.ve3.clone().unwrap_or_else(Ve3ParamsV1::off);
    let v4 = ve.ve4.clone().unwrap_or_else(Ve4ParamsV1::off);
    let v5 = ve.ve5.clone().unwrap_or_else(Ve5ParamsV1::off);
    let pos = |i: usize| { let (gx, gy) = lat.g(i); (gx as f64 * dx, gy as f64 * dx) };
    let skip = |i: usize, fixed: &[u8], z: &[f64]| fixed[i] == CHANNEL || fixed[i] == LAKE || t1a_post[i] || !z[i].is_finite();
    let z_body: Vec<f64> = z.to_vec();
    // Body styles.
    for i in 0..nn {
        if over[i] != OVER_EDIFICE || vid[i] == 255 { continue; }
        let Some(v) = cx.volcanoes.get(vid[i] as usize) else { continue };
        let (x, y) = pos(i);
        style[i] = match (&v.ve, family_of(&v.kind)) {
            (Some(_), Some(f)) => {
                let r = (x - v.vent_xy_m[0]).hypot(y - v.vent_xy_m[1]);
                if f == Family::Composite && r < 0.65 * v.crater_radius_m { ST_FLOOR } else if f == Family::Composite && r < 1.1 * v.crater_radius_m { ST_WALL } else { f.body_style() }
            }
            _ => if v.denudation_m >= v.h_m { ST_NECK } else { family_of(&v.kind).map(|f| f.body_style()).unwrap_or(ST_BODY) },
        };
        unit[i] = UNIT_BODY;
    }
    let mut best_d = vec![f64::INFINITY; 0];
    let mut best_s = vec![0f64; 0];
    let mut best_b = vec![0f64; 0];
    let mut best_lv = vec![0f64; 0];
    let mut best_ang = vec![0f64; 0];
    let mut touched: Vec<usize> = Vec::new();
    for &k in vsi {
        let v = &cx.volcanoes[k];
        let Some(g) = v.ve.as_ref() else { continue };
        if g.units.is_empty() { continue; }
        st.volcanoes_with_units += 1;
        let kid = k.min(254) as u8;
        let vent = v.vent_xy_m;
        // Lane VE5: summit taper radius (crater or caldera unit).
        let r_t = g.units.iter().filter_map(|u| match u.geom { VeGeom::Caldera { r_m, .. } => Some(r_m), _ => None }).fold(v.crater_radius_m, f64::max).max(1.0);
        for u in &g.units {
            let a = u.erosion_m;
            // Each unit: compute a new height for the posts it touches.
            let set = |i: usize, znew: f64, sty: u8, z: &mut [f64], fixed: &mut [u8], over: &mut [u8], vid: &mut [u8], style: &mut [u8], unit: &mut [u16], st: &mut VeCoastStats| {
                if skip(i, fixed, z) { return; }
                let d = znew - z[i];
                if d.abs() < 1e-9 { return; }
                if d > 0.0 { st.unit_raised_posts += 1; } else { st.unit_cut_posts += 1; }
                // Mappable change: the unit is exposed (edifice post); thinner veneers change the height only.
                if d.abs() >= tmap { over[i] = OVER_EDIFICE; vid[i] = kid; style[i] = sty; unit[i] = u.id; }
                if fixed[i] == GUIDE && znew > z[i] + 1.0 { fixed[i] = FREE; st.buried_channel_posts += 1; }
                if fixed[i] == MARINE && znew > sea { fixed[i] = FREE; }
                z[i] = znew;
            };
            match &u.geom {
                VeGeom::Path { pts, len_m, w0_m, t_m, lahar, lobes, ponds, .. } => {
                    let t = (t_m - a).max(0.0);
                    if t <= 0.0 { continue; }
                    let fcap = if v4.front == 1 && !front.is_empty() { v4.front_cap_mm as f64 / 1000.0 } else { 0.0 };
                    let wmax = w0_m * 1.5 * 1.3 + 2.0 * dx + fcap;
                    let nseg = pts.len() - 1;
                    if best_d.len() != nn { best_d = vec![f64::INFINITY; nn]; best_s = vec![0.0; nn]; best_b = vec![0.0; nn]; best_lv = vec![0.0; nn]; best_ang = vec![0.0; nn]; }
                    touched.clear();
                    for kseg in 0..nseg {
                        let (p, q) = (pts[kseg], pts[kseg + 1]);
                        let (ddx, ddy) = (q[0] - p[0], q[1] - p[1]);
                        let l2 = ddx * ddx + ddy * ddy;
                        if l2 <= 0.0 { continue; }
                        let ln = l2.sqrt();
                        let Some((a0, a1, b0, b1)) = post_range(lat, p[0].min(q[0]) - wmax, p[1].min(q[1]) - wmax, p[0].max(q[0]) + wmax, p[1].max(q[1]) + wmax) else { continue };
                        for gy in b0..=b1 { for gx in a0..=a1 {
                            let Some(i) = lat.idx(gx, gy) else { continue };
                            let (x, y) = (gx as f64 * dx, gy as f64 * dx);
                            let tt = ((x - p[0]) * ddx + (y - p[1]) * ddy) / l2;
                            let last = kseg == nseg - 1;
                            let (dist, s, beyond, ang) = if last && tt > 1.0 {
                                let (ex, ey) = (x - q[0], y - q[1]);
                                (ex.hypot(ey), q[2], (ex * ddx + ey * ddy) / ln, ey.atan2(ex) - ddy.atan2(ddx))
                            } else {
                                let tc = tt.clamp(0.0, 1.0);
                                let (cx_, cy_) = (p[0] + tc * ddx, p[1] + tc * ddy);
                                ((x - cx_).hypot(y - cy_), p[2] + tc * (q[2] - p[2]), 0.0, 0.0)
                            };
                            if dist > wmax { continue; }
                            if !(dist < best_d[i]) { continue; }
                            if best_d[i].is_infinite() { touched.push(i); }
                            best_d[i] = dist; best_s[i] = s; best_b[i] = beyond; best_ang[i] = ang;
                            // Level: both ends captured -> the bed interpolated; else the nearer end's level (bed or drape kind).
                            let tc = tt.clamp(0.0, 1.0);
                            best_lv[i] = if p[3] > -1e8 && q[3] > -1e8 { p[3] + tc * (q[3] - p[3]) } else if tc <= 0.5 { p[3] } else { q[3] };
                        } }
                    }
                    let l = *len_m;
                    for &i in &touched {
                        let (dist, s, beyond) = (best_d[i], best_s[i], best_b[i]);
                        best_d[i] = f64::INFINITY;
                        let lv = best_lv[i];
                        if *lahar && (lv - LV_BODY).abs() < 1.0 { continue; }
                        let w = w0_m * (1.0 + 0.5 * (s / l.max(1.0)).min(1.0));
                        let weff = if beyond > 0.0 { w * (1.0 + 0.25 * lobe(lobes, best_ang[i])) } else { w };
                        let m = if *lahar && v2.taper == 1 { (0.5 * weff).max(dx) } else if v3.aa == 1 { (t / 0.7).max(1.5 * dx).min(weff.max(1.5 * dx)) } else { (t / 0.7).max(dx).min(weff.max(dx)) };
                        // Lane VE4 front field: signed distance to the snout outline (youngest flow wins: units run oldest first).
                        if fcap > 0.0 && !*lahar && beyond > 0.0 && (weff - dist).abs() <= fcap && !skip(i, fixed, z) { front[i] = weff - dist; }
                        let prof = smoothstep01((weff - dist) / m);
                        if prof <= 0.0 { continue; }
                        let ramp = (s / 200.0).min(1.0);
                        let th = t * prof * ramp * if v5.summit == 1 { let (px, py) = pos(i); smoothstep01(((px - vent[0]).hypot(py - vent[1]) - r_t) / (0.5 * r_t)) } else { 1.0 };
                        if th <= 1e-6 { continue; }
                        let top = if lv > -1e8 { lv + th } else { z[i] + th };
                        let znew = z[i].max(top.min(z[i] + 3.0 * t));
                        let front = !*lahar && (beyond > 0.0 || s >= l - (1.5 * w).min(400.0));
                        let sty = if front { ST_FRONT } else { u.style };
                        if znew > z[i] { set(i, znew, sty, z, fixed, over, vid, style, unit, st); }
                    }
                    // Lane VE3 ponded lava: a flat top at the pond level (global 32 m posts; window lattices only).
                    if lat.dx_mm == 32_000 {
                        for pd in ponds {
                            for r in &pd.rows { for gx in r[1]..=r[2] {
                                let Some(i) = lat.idx(gx, r[0]) else { continue };
                                let zn = z[i].max(pd.level_m.min(z[i] + 3.0 * t));
                                if zn > z[i] { set(i, zn, if *lahar { ST_LAHAR } else { ST_FLOOR }, z, fixed, over, vid, style, unit, st); }
                            } }
                        }
                    }
                }
                VeGeom::Cone { x: cx0, y: cy0, rb_m, h_m, rc_m, dc_m, breach_rad } => {
                    let hh = h_m - a;
                    if hh <= 0.0 { continue; }
                    let dc = (dc_m - a).max(0.0);
                    let Some((a0, a1, b0, b1)) = post_range(lat, cx0 - rb_m, cy0 - rb_m, cx0 + rb_m, cy0 + rb_m) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let (x, y) = (gx as f64 * dx, gy as f64 * dx);
                        let r = (x - cx0).hypot(y - cy0);
                        if r >= *rb_m { continue; }
                        let mut hc = if r >= *rc_m { hh * ((rb_m - r) / (rb_m - rc_m)).clamp(0.0, 1.0) } else { hh - dc * (1.0 - (r / rc_m).powi(2)) };
                        if let Some(bd) = breach_rad {
                            let dd = ang_diff((y - cy0).atan2(x - cx0), *bd).abs();
                            if dd < 0.45 { let notch = (hh - dc) * ((rb_m - r) / rb_m).clamp(0.0, 1.0); hc = hc.min(notch.max(0.0)); }
                        }
                        if hc <= 0.0 { continue; }
                        let sty = if r < 0.5 * rc_m { ST_FLOOR } else if r < *rc_m { ST_WALL } else { ST_CINDER };
                        set(i, z[i] + hc, sty, z, fixed, over, vid, style, unit, st);
                    } }
                }
                VeGeom::Dome { x: cx0, y: cy0, r_m, h_m, top_m } => {
                    let hh = h_m - a;
                    if hh <= 0.0 { continue; }
                    let Some((a0, a1, b0, b1)) = post_range(lat, cx0 - r_m, cy0 - r_m, cx0 + r_m, cy0 + r_m) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let r = (gx as f64 * dx - cx0).hypot(gy as f64 * dx - cy0);
                        if r >= *r_m { continue; }
                        let hd = hh * (1.0 - (r / r_m).powi(4)).max(0.0).sqrt();
                        // Lane VE2 dome cap: the top stays at or below the crater rim.
                        let zn = match top_m { Some(t) => z[i].max((z[i] + hd).min(*t)), None => z[i] + hd };
                        if zn > z[i] { set(i, zn, ST_DOME, z, fixed, over, vid, style, unit, st); }
                    } }
                }
                VeGeom::Maar { x: cx0, y: cy0, r_m, d_m, rim_m, ring_b_m } => {
                    let d = d_m - a;
                    if d <= 0.0 { continue; }
                    let rim = (rim_m - a).max(0.0);
                    let ro = r_m + 4.0 * ring_b_m;
                    let Some((a0, a1, b0, b1)) = post_range(lat, cx0 - ro, cy0 - ro, cx0 + ro, cy0 + ro) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let r = (gx as f64 * dx - cx0).hypot(gy as f64 * dx - cy0);
                        if r >= ro { continue; }
                        let (hm, sty) = if r < *r_m { (rim - (rim + d) * (1.0 - (r / r_m).powi(2)), if r < 0.7 * r_m { ST_FLOOR } else { ST_WALL }) } else { (rim * (-(r - r_m) / ring_b_m).exp(), ST_TUFF) };
                        if hm.abs() < 1e-6 { continue; }
                        set(i, z[i] + hm, sty, z, fixed, over, vid, style, unit, st);
                    } }
                }
                VeGeom::Pit { x: cx0, y: cy0, r_m, d_m } => {
                    let d = d_m - a;
                    if d <= 0.0 { continue; }
                    let Some((a0, a1, b0, b1)) = post_range(lat, cx0 - r_m, cy0 - r_m, cx0 + r_m, cy0 + r_m) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let r = (gx as f64 * dx - cx0).hypot(gy as f64 * dx - cy0);
                        if r >= *r_m { continue; }
                        let hp = -d * smoothstep01((r_m - r) / (d / 0.9));
                        let sty = if hp < -0.9 * d { ST_FLOOR } else { ST_WALL };
                        set(i, z[i] + hp, sty, z, fixed, over, vid, style, unit, st);
                    } }
                }
                VeGeom::Fissure { x0, y0, x1, y1, h_m, w_m } => {
                    let hh = h_m - a;
                    if hh <= 0.0 { continue; }
                    let rr = 3.0 * w_m;
                    let (ddx, ddy) = (x1 - x0, y1 - y0);
                    let l = ddx.hypot(ddy).max(1.0);
                    let Some((a0, a1, b0, b1)) = post_range(lat, x0.min(*x1) - rr, y0.min(*y1) - rr, x0.max(*x1) + rr, y0.max(*y1) + rr) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let (x, y) = (gx as f64 * dx, gy as f64 * dx);
                        let al = ((x - x0) * ddx + (y - y0) * ddy) / l;
                        if al < 0.0 || al > l { continue; }
                        let la = ((x - x0) * ddy - (y - y0) * ddx).abs() / l;
                        let tap = smoothstep01(al.min(l - al) / 200.0);
                        let hf = hh * (-(la / w_m).powi(2)).exp() * tap;
                        if hf < 0.05 { continue; }
                        set(i, z[i] + hf, ST_FISSURE, z, fixed, over, vid, style, unit, st);
                    } }
                }
                VeGeom::Caldera { x: cx0, y: cy0, r_m, floor_m, r_floor_m, wall } => {
                    let fl = floor_m + a;
                    let Some((a0, a1, b0, b1)) = post_range(lat, cx0 - r_m, cy0 - r_m, cx0 + r_m, cy0 + r_m) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let r = (gx as f64 * dx - cx0).hypot(gy as f64 * dx - cy0);
                        if r > *r_m { continue; }
                        let zc = fl + wall * (r - r_floor_m).max(0.0);
                        if zc < z[i] { let sty = if r <= *r_floor_m { ST_FLOOR } else { ST_WALL }; set(i, zc, sty, z, fixed, over, vid, style, unit, st); }
                    } }
                }
                VeGeom::Breach { x: cx0, y: cy0, dir_rad, half_w_m, floor0_m, slope, len_m, back_m, wall } => {
                    let (c, s) = (dir_rad.cos(), dir_rad.sin());
                    let rr = len_m + half_w_m + 2000.0;
                    let Some((a0, a1, b0, b1)) = post_range(lat, cx0 - rr, cy0 - rr, cx0 + rr, cy0 + rr) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let (x, y) = (gx as f64 * dx - cx0, gy as f64 * dx - cy0);
                        let (al, la) = (x * c + y * s, (-x * s + y * c).abs());
                        if al > *len_m { continue; }
                        let fl = floor0_m + a - slope * al.max(0.0);
                        let zc = fl + wall * (la - half_w_m).max(0.0) + wall * (-back_m - al).max(0.0);
                        if zc < z[i] { let sty = if zc - fl < 10.0 { ST_FLOOR } else { ST_WALL }; set(i, zc, sty, z, fixed, over, vid, style, unit, st); }
                    } }
                }
                VeGeom::Apron { x: cx0, y: cy0, dir_rad, half_rad, r_a_m, r_max_m, s_f, t_cap_m, t0_m, b_m, za } => {
                    let tc = (t_cap_m - a).max(0.0);
                    let t0 = (t0_m - a).max(0.0);
                    if tc <= 0.0 && t0 <= 0.0 { continue; }
                    let Some((a0, a1, b0, b1)) = post_range(lat, cx0 - r_max_m, cy0 - r_max_m, cx0 + r_max_m, cy0 + r_max_m) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let (x, y) = (gx as f64 * dx - cx0, gy as f64 * dx - cy0);
                        let r = x.hypot(y);
                        if r >= *r_max_m { continue; }
                        let th = y.atan2(x);
                        let dth = ang_diff(th, *dir_rad);
                        // Density currents: smooth directional weight (cos lobe of the declared half-angle; full ring at 180 deg).
                        let gw = if *half_rad >= std::f64::consts::PI - 1e-6 { 1.0 } else { (0.5 + 0.5 * (dth * std::f64::consts::PI / (2.0 * half_rad)).cos()).clamp(0.0, 1.0).powi(2) * (dth.abs() < 2.0 * half_rad) as i32 as f64 };
                        let tap = smoothstep01((r_max_m - r) / (if v2.taper == 1 { 0.5 } else { 0.3 } * r_max_m));
                        // Fall: elliptical isopachs downwind (dir_rad = dispersal axis).
                        let (al, cr) = (r * dth.cos(), r * dth.sin());
                        let re = (if al > 0.0 { al / 2.5 } else { al / 0.6 }).hypot(cr);
                        let mut zn = z[i] + t0 * (-re / b_m).exp() * tap;
                        if r >= *r_a_m && tc > 0.0 && gw > 0.0 {
                            let onset = smoothstep01((r - r_a_m) / if v2.taper == 1 { 3000.0 } else { 1500.0 });
                            let fk = th.rem_euclid(std::f64::consts::TAU) / std::f64::consts::TAU * 72.0;
                            let (k0, f) = (fk.floor() as usize % 72, fk - fk.floor());
                            let zr = za[k0] * (1.0 - f) + za[(k0 + 1) % 72] * f;
                            let fsurf = zr - s_f * (r - r_a_m);
                            zn = zn.max(fsurf.min(zn + tc * gw * tap * onset));
                        }
                        if zn > z[i] { set(i, zn, ST_APRON, z, fixed, over, vid, style, unit, st); }
                    } }
                }
                VeGeom::Sheet { x: cx0, y: cy0, r_m, t_m, lobes } => {
                    let t = t_m - a;
                    if t <= 0.0 { continue; }
                    let ro = r_m * 1.3;
                    let Some((a0, a1, b0, b1)) = post_range(lat, cx0 - ro, cy0 - ro, cx0 + ro, cy0 + ro) else { continue };
                    let m = (t / 0.7).max(dx);
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        let (x, y) = (gx as f64 * dx - cx0, gy as f64 * dx - cy0);
                        let r = x.hypot(y);
                        let rk = r_m * (1.0 + 0.15 * lobe(lobes, y.atan2(x)));
                        let th = t * smoothstep01((rk - r) / m);
                        if th <= 1e-6 { continue; }
                        set(i, z[i] + th, ST_PLATEAU, z, fixed, over, vid, style, unit, st);
                    } }
                }
                VeGeom::Gullies { .. } => {
                    let Some((a0, a1, b0, b1)) = post_range(lat, vent[0] - g.r_body_m, vent[1] - g.r_body_m, vent[0] + g.r_body_m, vent[1] + g.r_body_m) else { continue };
                    for gy in b0..=b1 { for gx in a0..=a1 {
                        let Some(i) = lat.idx(gx, gy) else { continue };
                        if over[i] != OVER_EDIFICE || vid[i] != kid { continue; }
                        let c = gully_cut(g, v, cf, gx as f64 * dx, gy as f64 * dx);
                        if c <= 0.0 { continue; }
                        // The gully cuts the cone, never below the pre-edifice surface.
                        let zn = (z[i] - c).max(z_pre[i].min(z[i]));
                        if zn < z[i] { let keep = (style[i], unit[i]); set(i, zn, keep.0, z, fixed, over, vid, style, unit, st); style[i] = keep.0; unit[i] = keep.1; }
                    } }
                }
            }
            st.units_in_lattice += 1;
        }
    }
    let esc_trace = crate::t1b_esc_deps::tracing();
    // T1a corridor cap on unit-raised posts (with the body): relief over the pre-edifice surface <= cap_slope (d - mask_inner).
    let (mi, cs) = (vp.mask_inner_mm as f64 / 1000.0, vp.cap_slope_ppm as f64 / 1e6);
    for i in 0..nn {
        if z[i] <= z_body[i] + 1e-9 || !z_pre[i].is_finite() { continue; }
        let rel = z[i] - z_pre[i];
        if rel <= 0.0 { continue; }
        let (x, y) = pos(i);
        let dneed = mi + rel / cs;
        let mut d2 = dneed * dneed;
        if esc_trace { crate::t1b_esc_deps::cap_begin(1,x,y,dneed); }
        for s in &cx.t1a_segs {
            if (s[0].min(s[2]) - x) > dneed || (x - s[0].max(s[2])) > dneed || (s[1].min(s[3]) - y) > dneed || (y - s[1].max(s[3])) > dneed { continue; }
            if esc_trace { crate::t1b_esc_deps::cap_segment(s); }
            let q = seg_nearest([x, y], [s[0], s[1]], [s[2], s[3]]);
            d2 = d2.min((q[0] - x).powi(2) + (q[1] - y).powi(2));
        }
        if esc_trace { crate::t1b_esc_deps::geometry_end(); }
        let allowed = z_pre[i] + cs * (d2.sqrt() - mi).max(0.0);
        let lim = allowed.max(z_body[i]);
        if z[i] > lim { z[i] = lim; st.corridor_capped_posts += 1; if z[i] - z_body[i] < tmap && unit[i] != UNIT_BODY && unit[i] != UNIT_NONE && z_body[i] > z_pre[i] { unit[i] = UNIT_BODY; } }
    }
    for i in 0..nn { if style[i] != ST_NONE { *st.posts_by_style.entry(style_name(style[i]).into()).or_insert(0) += 1; } }
}

/// After the final fill: declared lakes (posts the fill raised that lie on a crater / caldera / maar / pit floor or wall).
pub fn lake_stats(z_before: &[f64], z: &[f64], style: &[u8], st: &mut VeCoastStats) {
    for i in 0..z.len() {
        let r = z[i] - z_before[i];
        if r <= 1e-9 || !r.is_finite() { continue; }
        if style[i] == ST_FLOOR || style[i] == ST_WALL { st.declared_lake_posts += 1; st.declared_lake_raise_max_m = st.declared_lake_raise_max_m.max(r); } else { st.other_fill_raised_posts += 1; }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wind_deviation_is_a_wrapped_cauchy_with_mean_resultant_rho() {
        // Lane VH-I: median 0, monotone in u, E[cos] = rho (wrapped Cauchy), within (-pi, pi].
        assert!(wind_deviation(0.698, 0.5).abs() < 1e-15);
        let n = 200_000;
        let (mut c, mut prev) = (0.0, -10.0);
        for k in 0..n { let d = wind_deviation(0.698, (k as f64 + 0.5) / n as f64); assert!(d > prev && d.abs() <= std::f64::consts::PI); prev = d; c += d.cos(); }
        assert!((c / n as f64 - 0.698).abs() < 2e-3, "mean resultant {}", c / n as f64);
        let off = VeParamsV1::v1();
        assert!(!serde_json::to_string(&off).unwrap().contains("\"wind\""), "off: params JSON unchanged");
    }
    #[test]
    fn family_mapping_and_clearance() {
        assert_eq!(family_of("continental_arc"), Some(Family::Composite));
        assert_eq!(family_of("island_arc"), Some(Family::Composite));
        assert_eq!(family_of("hotspot"), Some(Family::Shield));
        assert_eq!(family_of("rift"), Some(Family::Rift));
        assert_eq!(family_of("large_igneous_province"), Some(Family::Plateau));
        assert!(broad_clearance("continental_arc", 3000.0).is_none());
        // Shield H 3000 m: W_B 60 km, R0 = 30 / sqrt(ln 20) km; clearance R0 sqrt(ln 5).
        let c = broad_clearance("hotspot", 3000.0).unwrap();
        let r0 = 30_000.0 / 20f64.ln().sqrt();
        assert!((c - r0 * 5f64.ln().sqrt()).abs() < 1e-6);
        // At that radius the Gaussian is 20% of H.
        assert!(((-(c / r0).powi(2)).exp() - 0.2).abs() < 1e-9);
    }

    #[test]
    fn lobe_is_bounded() {
        let l = [0.3, 1.1, 2.0, 4.0];
        let m = (0..3600).map(|k| lobe(&l, k as f64 * std::f64::consts::TAU / 3600.0).abs()).fold(0.0, f64::max);
        assert!(m <= 1.0 + 1e-9, "max {m}");
    }

    /// Lane VE3: a synthetic older flow (a 30 m ridge across a planar 10% slope) deflects a younger steepest-descent path.
    #[test]
    fn older_ridge_deflects_a_younger_path() {
        // Plane falling toward +y at 0.1; an older flow laid along a line at 45 deg across the fall line.
        let older = VeUnit { pending: None, id: 1, kind: "lava_flow".into(), style: ST_LAVA, age_yr: 2000.0, vei: 0, size_m: 30.0, erosion_m: 0.0,
            geom: VeGeom::Path { pts: (0..=40).map(|k| { let s = k as f64 * 100.0; [-2000.0 + s * 0.7071, 1000.0 + s * 0.7071, s, LV_BODY] }).collect(),
                len_m: 4000.0, w0_m: 150.0, t_m: 30.0, lahar: false, lobes: [0.0; 4], captured_at_m: -1.0, end: "length".into(), ponds: Vec::new() } };
        let g = |q: [f64; 2]| -0.1 * q[1] + unit_drape_m(&older, q[0], q[1], 48.0);
        let mut p = [0.0f64, 0.0];
        let mut xs = Vec::new();
        for _ in 0..80 {
            let h = 64.0;
            let gr = [(g([p[0] + h, p[1]]) - g([p[0] - h, p[1]])) / (2.0 * h), (g([p[0], p[1] + h]) - g([p[0], p[1] - h])) / (2.0 * h)];
            let l = gr[0].hypot(gr[1]);
            p = [p[0] - 64.0 * gr[0] / l, p[1] - 64.0 * gr[1] / l];
            xs.push(p);
        }
        // Without the ridge the path would stay at x = 0 (straight downslope). The ridge (from (-2000, 1000) rising to the
        // upper right) lies across the fall line ahead; the path must be pushed sideways (|x| > 300 m) before passing it.
        let maxdx = xs.iter().map(|q| q[0].abs()).fold(0.0, f64::max);
        assert!(maxdx > 300.0, "deflection {maxdx}");
        // And it never climbs over the ridge's crest line while beside it: every point stays on one side of the ridge axis.
        let side = |q: [f64; 2]| (q[0] + 2000.0) * 0.7071 - (q[1] - 1000.0) * 0.7071; // >0 below the ridge (toward +x/-y)
        let s0 = side(xs[0]);
        assert!(xs.iter().all(|q| side(*q) * s0 > 0.0), "path crossed the older ridge");
    }

    /// Lane VE5 bug fix: on an open slope (no closed low inside the radius) the old flood filled against the search box
    /// (a rectangle); with edge_spill it spills at the box edge and builds no pond.
    #[test]
    fn open_slope_is_not_a_closed_low() {
        let g = |q: [f64; 2]| -0.05 * q[1] - 0.01 * q[0];
        let (old, _) = fill_pond_e(&g, [0.0, 0.0], 500.0, 1e12, false);
        assert!(old.cells > 0, "the ve3 / ve4 behaviour fills against the box edge");
        let (new, sp) = fill_pond_e(&g, [0.0, 0.0], 500.0, 1e12, true);
        assert!(new.cells == 0 && sp.is_some(), "cells {} spill {:?}", new.cells, sp);
    }

    #[test]
    fn pond_fills_a_closed_low_to_its_spill_level() {
        // A bowl (depth 20 m, radius 400 m) on a 5% plane: the fill reaches the spill level, the spill is downslope.
        let g = |q: [f64; 2]| -0.05 * q[1] - 20.0 * (1.0 - (q[0].hypot(q[1]) / 400.0).powi(2)).max(0.0);
        let (pd, spill) = fill_pond(&g, [0.0, 0.0], 1500.0, 1e12);
        assert!(pd.cells > 50 && spill.is_some(), "cells {} spill {:?}", pd.cells, spill);
        assert!(spill.unwrap()[1] > 0.0, "spills downslope (+y)");
        // Volume-limited: a small volume gives a lower level and no spill.
        let (pd2, sp2) = fill_pond(&g, [0.0, 0.0], 1500.0, 1e5);
        assert!(pd2.level_m < pd.level_m && sp2.is_none() && (pd2.volume_m3 - 1e5).abs() < 1e3, "{} {} {:?}", pd2.level_m, pd2.volume_m3, sp2);
    }

    fn lat(n: usize) -> EvalLat { EvalLat { gx0: -(n as i64) / 2, gy0: -(n as i64) / 2, n, dx_mm: 32_000 } }

    /// Lane VE5: front_mm from apply_units_f on a synthetic straight flow: > 0 exactly where the snout deposit is > 0, < 0 just
    /// outside, so its zero contour is the drape outline.
    #[test]
    fn front_field_zero_contour_is_the_snout_outline() {
        use crate::t1b_coast::{CfParamsV1, CoastCtx, VolcanoParamsV1, VolcanoSite};
        let l = lat(161);
        let nn = l.nodes();
        let pts: Vec<[f64; 4]> = (0..=25).map(|k| [-2000.0 + 100.0 * k as f64, 0.0, 100.0 * k as f64, LV_BODY]).collect();
        let unit = VeUnit { pending: None, id: 1, kind: "lava_flow".into(), style: ST_LAVA, age_yr: 1000.0, vei: 0, size_m: 20.0, erosion_m: 0.0,
            geom: VeGeom::Path { pts, len_m: 2500.0, w0_m: 150.0, t_m: 20.0, lahar: false, lobes: [0.3, 1.0, 2.0, 4.0], captured_at_m: -1.0, end: "length".into(), ponds: Vec::new() } };
        let g = VeVolcano { flank2: None, flank: None, family: "composite".into(), span_my: None, t_h_ky: 40.0, lifetime_my: 0.5, r_body_m: 5000.0, r0_m: 0.0, axes_rad: vec![], stretch: 0.0,
            axis_w_m: 2500.0, elong: 1.0, q_built: 0.15, r_units_m: 6000.0, gullies: vec![], rough: None, counts: BTreeMap::new(), dropped_eroded: 0, dropped_cap: 0, units: vec![unit], wind: None };
        let v = VolcanoSite { id: "t".into(), kind: "continental_arc".into(), shape: "stratovolcano".into(), record_xy_m: [-2000.0, 0.0], vent_xy_m: [-2000.0, 0.0], moved_m: 0.0,
            record_cell: 0, vent_cell: 0, summit_record_m: 1000.0, base_m: 0.0, h_m: 1000.0, lambda_m: 1000.0, r_main_m: 5000.0, r_influence_m: 1.0, s_top: 0.65, s_base: 0.09,
            crater_radius_m: 350.0, crater_depth_m: 150.0, age_my: 0.0, e_vent_m_per_my: 0.0, denudation_m: 0.0, stage: String::new(), breached: false, breach_dir_rad: 0.0,
            gully_dirs_rad: vec![], dike_dirs_rad: vec![], neck_relief_m: 0.0, lobe_phase_rad: vec![], lobe_norm: 1.0, cf: None, ve: Some(g), ped_shift_m: None, built_summit_m: None };
        let mut vep = VeParamsV1::v1();
        vep.ve4 = Some(Ve4ParamsV1::v1());
        let cx = CoastCtx { cf: Some(CfParamsV1::v1()), ve: Some(vep), ve_units: true, volcanoes: vec![v], ..Default::default() };
        let vp = VolcanoParamsV1::v1();
        let z_pre = vec![100.0; nn];
        let mut z = z_pre.clone();
        let (mut fixed, t1a, mut over, mut vid, mut sty, mut un) = (vec![FREE; nn], vec![false; nn], vec![0u8; nn], vec![255u8; nn], vec![ST_NONE; nn], vec![UNIT_NONE; nn]);
        let mut front = vec![f64::NAN; nn];
        let mut st = VeCoastStats::default();
        apply_units_f(&l, &vp, &cx, &[0], &z_pre, &mut z, &mut fixed, &t1a, &mut over, &mut vid, &mut sty, &mut un, 0.0, &mut st, &mut front);
        let (mut pos, mut neg) = (0, 0);
        for i in 0..nn {
            let f = front[i];
            if !f.is_finite() { continue; }
            let raised = z[i] > 100.0 + 1e-9;
            assert_eq!(f > 1e-9, raised, "post {i}: front {f} raised {raised} dz {}", z[i] - 100.0);
            if f > 0.0 { pos += 1; } else { neg += 1; }
            let (gx, _) = l.g(i);
            assert!(gx as f64 * 32.0 > 500.0 - 1e-6, "front only past the path end");
        }
        assert!(pos > 10 && neg > 10, "both sides sampled: {pos} inside, {neg} outside");
    }

    #[test]
    fn cinder_cone_shape_ratios_on_the_lattice() {
        // A 900 m cone on a flat surface: built H / W and crater / base widths (32 m sampling).
        let l = lat(81);
        let nn = l.nodes();
        let g = VeGeom::Cone { x: 0.0, y: 0.0, rb_m: 450.0, h_m: 162.0, rc_m: 180.0, dc_m: 48.6, breach_rad: None };
        let mut z = vec![100.0; nn];
        let (dx, n) = (l.dx(), l.n);
        let mut hmax: f64 = 0.0;
        let (mut rim_r, mut base_r): (f64, f64) = (0.0, 0.0);
        if let VeGeom::Cone { rb_m, h_m, rc_m, dc_m, .. } = g {
            for i in 0..nn {
                let (gx, gy) = l.g(i);
                let r = (gx as f64 * dx).hypot(gy as f64 * dx);
                if r >= rb_m { continue; }
                let hc = if r >= rc_m { h_m * ((rb_m - r) / (rb_m - rc_m)).clamp(0.0, 1.0) } else { h_m - dc_m * (1.0 - (r / rc_m).powi(2)) };
                z[i] += hc;
                if hc > hmax { hmax = hc; rim_r = r; }
                if hc > 1.0 { base_r = base_r.max(r); }
            }
        }
        let _ = n;
        let ratio = hmax / (2.0 * base_r);
        assert!(ratio > 0.16 && ratio < 0.20, "H/W {ratio}");
        let cr = rim_r / base_r;
        assert!(cr > 0.3 && cr < 0.5, "crater/base {cr}");
    }
}

#[cfg(test)]
mod flank2_tests {
    use super::*;
    #[test]
    fn native_flank2_conserves_sediment_and_fixed_posts() {
        let n = 65;
        let dx = 32.;
        let nn = n * n;
        let z: Vec<_> = (0..nn)
            .map(|i| {
                let x = i % n;
                let y = i / n;
                let r = ((x as f64 - 32.).powi(2) + (y as f64 - 32.).powi(2)).sqrt() * dx;
                600. * (-r / 430.).exp() + 0.01 * x as f64
            })
            .collect();
        let mut fix = vec![FREE; nn];
        fix[32 * n + 32] = LAKE;
        let rain = vec![1.; nn];
        let stock = vec![600.; nn];
        let mat = vec![0.; nn];
        let (d, l) =
            flank2_evolve(n, &z, &fix, &rain, &stock, &mat, &rain, dx, 30_000., 120).unwrap();
        assert_eq!(d[32 * n + 32], 0.);
        for i in 0..nn {
            if i % n == 0 || i / n == 0 || i % n == n - 1 || i / n == n - 1 {
                assert_eq!(d[i], 0.);
            }
        }
        assert!(l["eroded_m3"].as_f64().unwrap() > 0.);
        assert!(l["deposited_m3"].as_f64().unwrap() > 0.);
        assert!(l["water_relative_residual"].as_f64().unwrap() < 1e-9);
        assert!(l["volume_relative_residual"].as_f64().unwrap() < 1e-8);
    }
    #[test]
    fn native_flank2_requires_material_and_rejects_nonfinite() {
        let n = 9;
        let nn = n * n;
        let mut z: Vec<_> = (0..nn).map(|i| (nn - i) as f64).collect();
        let zero = vec![0.; nn];
        let one = vec![1.; nn];
        let fixed = vec![FREE; nn];
        let (d, _) = flank2_evolve(n, &z, &fixed, &one, &zero, &zero, &one, 32., 1000., 4).unwrap();
        assert!(d.iter().all(|&x| x == 0.));
        z[20] = f64::NAN;
        assert!(flank2_evolve(n, &z, &fixed, &one, &one, &zero, &one, 32., 1000., 4).is_err());
    }
    #[test]
    fn native_flank2_subnormal_transport_has_finite_audit() {
        // Depositional tails on flat native terrain can create tiny nonzero slopes.
        // Casting routing weights to f32 before normalization previously made0/0.
        let n = 9;
        let nn = n * n;
        let z: Vec<_> = (0..nn).map(|i| (nn - i) as f64 * 1e-90).collect();
        let fixed = vec![FREE; nn];
        let one = vec![1.; nn];
        let zero = vec![0.; nn];
        let (_, l) = flank2_evolve(n, &z, &fixed, &one, &zero, &zero, &one, 32., 1000., 4).unwrap();
        assert!(l["exported_m3"].as_f64().unwrap().is_finite());
        assert!(l["volume_relative_residual"].as_f64().unwrap() < 1e-8);
    }
    #[test]
    fn native_flank2_smooth_apron_and_parent_serialization() {
        let f = |r| flank2_profile(2500., 720., 0.65, 4711.8, 350., 8333.3, r);
        for r in [350. + 4990., 8333.3, 4. * 8333.3] {
            let a = (f(r) - f(r - 0.1)) / 0.1;
            let b = (f(r + 0.1) - f(r)) / 0.1;
            assert!((a - b).abs() < 0.0001);
            assert!(f(r) > 0.);
        }
        let p = VeParamsV1::v1();
        assert!(serde_json::to_value(p).unwrap().get("flank2").is_none());
    }
}

/// Arc3 owner-driven process on an already closed, supplied native domain.
/// This entry point never creates a perimeter outlet or fills physical ground.
/// The old inferred-stock/material proxy is explicit; this is not qualified
/// geological stock, eruption contacts, strength or exposure history.
pub const FLANK3_OWNER_RULE:&str="ve_flank3_owner_v1: supplied closed native domain, LS MFD slope^1.1/link-length and registered graph at each accepted revision; 30ky/120 steps baseline, 240-step refinement; retained 32/My max(sqrt(A)S-12m,0)/(1+M/20m) exposure proxy; retained 0.02m/yr A S^2 transport capacity and 200m settling length; finite supplied proxy stock, mixed by source identity; conservative bank transfer at retained 0.7+0.5M/(M+20m) repose, at most eight sweeps; immutable native boundaries and supplied bank floors; no physical post-process fill; final millimetre projection separately accounted and rerouted by owner; no qualified geology or transient lake claim";
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Flank3SourceLedger{pub source:u32,pub name:String,pub initial_stock_m3:f64,pub stream_detached_m3:f64,pub deposited_m3:f64,pub bank_transferred_m3:f64,pub exported_m3:f64,pub remaining_stock_m3:f64,pub absolute_residual_m3:f64,pub relative_residual:f64}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Flank3StepReport{pub step:usize,pub input_route_digest:String,pub accepted_domain_digest:String,pub route_digest:String,pub registry_digest:String,pub water_source_m2:f64,pub water_export_m2:f64,pub water_absolute_residual_m2:f64,pub water_relative_residual:f64,pub sediment_absolute_residual_m3:f64,pub sediment_relative_residual:f64,pub withheld_detachment_m3:f64}
pub struct Flank3Proxy<'a>{pub stock_depth_m:&'a[f64],pub material_m:&'a[f64],pub exposure_fraction:&'a[f64],pub floor_m:&'a[f64],pub source_id:&'a[u32],pub source_names:&'a BTreeMap<u32,String>}
pub struct Flank3Output{
    pub physical_ground_m:Vec<f64>,pub published_height_mm:Vec<i32>,pub final_domain:crate::t1b_surface_domain::Domain,
    pub route:crate::t1b_ls::DetailRoute,pub registry:crate::t1b_surface_channels::Registry,
    pub native_graph:crate::t1b_surface_channels::NativeGraphProduct,
    pub remaining_source_stock_m3:Vec<BTreeMap<u32,f64>>,
    pub sources:Vec<Flank3SourceLedger>,pub steps:Vec<Flank3StepReport>,pub report:serde_json::Value,
}
/// Transfer actual available source volumes. Geometry uses the actual sum taken,
/// not an assumed target. No provenance pruning or hidden deficit replacement.
fn flank3_take(stock:&mut BTreeMap<u32,f64>,target:f64)->BTreeMap<u32,f64>{
    let total: f64=stock.values().sum();let mut out=BTreeMap::new();if target<=0.0||total<=0.0{return out;}
    for (&s,v) in stock.iter_mut(){let take=(target*(*v/total)).min(*v);if take>0.0{*v-=take;out.insert(s,take);}}
    stock.retain(|_,v|*v!=0.0);out
}
fn flank3_revision(mut domain:crate::t1b_surface_domain::Domain,z:&[f64])->R<crate::t1b_surface_domain::Domain>{
    use crate::t1b_surface_domain::{component_digests,content_digest,VERSION};
    if z.len()!=domain.posts.len(){return Err("flank3 geometry dimensions".into());}
    let previous=domain.digest.clone();
    for (p,&v) in domain.posts.iter_mut().zip(z){if !v.is_finite()||(p.boundary.is_some()&&v!=p.ground_m){return Err("flank3 nonfinite or changed physical boundary".into());}p.ground_m=v;}
    domain.pin.predecessor=previous;domain.pin.ground=component_digests(&domain.posts).map_err(|e|e.to_string())?.0;
    domain.digest=content_digest(&(VERSION,&domain.pin,domain.policy,&domain.posts)).map_err(|e|e.to_string())?;Ok(domain)
}
/// Native reporting runoff/threshold are supplied independently of the erosion
/// flood-width hypothesis. A reporting threshold is never tuned by this law.
pub fn flank3_evolve_owner(initial:&crate::t1b_surface_domain::Domain,proxy:Flank3Proxy<'_>,steps:usize,runoff_reference_m_s:f64,threshold_m3_s:f64,on_step:&mut dyn FnMut(&Flank3StepReport)->R<()>)->R<Flank3Output>{
    use crate::t1b_surface_domain::{Mode,Policy};
    if ![120,240].contains(&steps)||!matches!(initial.policy,Policy::NativeDetailV1|Policy::NativePhysicalV1){return Err("flank3 requires native owner policy and 120 or 240 steps".into());}
    let n=initial.posts.len();let dx=initial.pin.spacing_mm as f64/1000.0;let dt=30000.0/steps as f64;
    if dx!=32.0||[proxy.stock_depth_m.len(),proxy.material_m.len(),proxy.exposure_fraction.len(),proxy.floor_m.len(),proxy.source_id.len()].iter().any(|&k|k!=n){return Err("flank3 requires native 32m geometry and complete proxy fields".into());}
    let mut stock=vec![BTreeMap::<u32,f64>::new();n];let mut ledger=BTreeMap::<u32,Flank3SourceLedger>::new();let z0:Vec<f64>=initial.posts.iter().map(|p|p.ground_m).collect();let mut z=z0.clone();
    for i in 0..n{
        if !proxy.stock_depth_m[i].is_finite()||proxy.stock_depth_m[i]<0.0||!proxy.material_m[i].is_finite()||proxy.material_m[i]<0.0||!proxy.exposure_fraction[i].is_finite()||!(0.0..=1.0).contains(&proxy.exposure_fraction[i])||!proxy.floor_m[i].is_finite()||proxy.floor_m[i]>z[i]{return Err("invalid proxy material, exposure, stock or bank floor".into());}
        if initial.posts[i].boundary.is_some()&&proxy.stock_depth_m[i]!=0.0{return Err("physical boundary cannot carry mobile flank stock".into());}
        if proxy.stock_depth_m[i]>0.0{
            let source=proxy.source_id[i];if source==u32::MAX{return Err("positive stock has unresolved source identity".into());}let name=proxy.source_names.get(&source).filter(|s|!s.is_empty()).ok_or("positive stock has no declared source identity")?;
            let volume=proxy.stock_depth_m[i]*initial.posts[i].support_m2;if !volume.is_finite(){return Err("nonfinite proxy stock volume".into());}
            stock[i].insert(source,volume);ledger.entry(source).or_insert(Flank3SourceLedger{source,name:name.clone(),initial_stock_m3:0.0,stream_detached_m3:0.0,deposited_m3:0.0,bank_transferred_m3:0.0,exported_m3:0.0,remaining_stock_m3:0.0,absolute_residual_m3:0.0,relative_residual:0.0}).initial_stock_m3+=volume;
        }
    }
    let kernel=blake3::hash(include_str!("t1b_ve.rs").split_once("pub const FLANK3_OWNER_RULE").ok_or("flank3 source anchor")?.1.as_bytes()).to_hex().to_string();
    let process_input=crate::t1b_surface_domain::content_digest(&(FLANK3_OWNER_RULE,initial.digest(),steps,runoff_reference_m_s,threshold_m3_s,proxy.stock_depth_m,proxy.material_m,proxy.exposure_fraction,proxy.floor_m,proxy.source_id,proxy.source_names,&kernel)).map_err(|e|e.to_string())?;
    let mut domain=initial.clone();domain.pin.recipe=process_input.clone();domain.pin.stage="flank3-physical".into();domain=flank3_revision(domain,&z)?;
    let mut reports=Vec::new();let mut max_water=0.0f64;let mut max_sediment=0.0f64;
    let mut route=crate::t1b_ls::route_detail_surface(&domain,Mode::MultipleFlowSlope11V1).map_err(|e|e.to_string())?;
    {
        let initial_registry=crate::t1b_ls::register_detail_channels(&domain,&route,runoff_reference_m_s,threshold_m3_s).map_err(|e|e.to_string())?;
        crate::t1b_surface_channels::audit(&domain,&route,&initial_registry).map_err(|e|e.to_string())?;
    }
    for step in 0..steps{
        let mut sediment=vec![BTreeMap::<u32,f64>::new();n];let mut withheld=0.0;
        for &iu in route.downstream_first.iter().rev(){let i=iu as usize;let a=route.effective_area_m2[i];let support=domain.posts[i].support_m2;let mut slope=0.0;
            for edge in &route.edges[i]{let j=edge.receiver as usize;let p=domain.posts[i].id;let q=domain.posts[j].id;let len=dx*((p.gx-q.gx)as f64).hypot((p.gy-q.gy)as f64);slope+=edge.weight*(route.routing_m[i]-route.routing_m[j]).max(0.0)/len;}
            let capacity=0.02*a*slope*slope*dt;if !capacity.is_finite(){return Err("flank3 nonfinite sediment capacity".into());}let incoming: f64=sediment[i].values().sum();
            if z[i]>=route.routing_m[i]-1e-9{
                let potential=32e-6*(a.sqrt()*slope-12.0).max(0.0)*dt/(1.0+proxy.material_m[i]/20.0)*proxy.exposure_fraction[i]*support;
                if !potential.is_finite(){return Err("flank3 nonfinite potential detachment".into());}
                let requested=potential.min((capacity-incoming).max(0.0)).min((z[i]-proxy.floor_m[i]).max(0.0)*support);
                let taken=flank3_take(&mut stock[i],requested);let actual:f64=taken.values().sum();withheld+=(potential-actual).max(0.0);z[i]-=actual/support;
                for (source,v) in taken{ledger.get_mut(&source).unwrap().stream_detached_m3+=v;*sediment[i].entry(source).or_default()+=v;}
            }
            let load:f64=sediment[i].values().sum();if load>capacity{
                let settled=flank3_take(&mut sediment[i],(load-capacity)*(1.0-(-dx/200.0).exp()));let volume:f64=settled.values().sum();z[i]+=volume/support;
                for (source,v) in settled{ledger.get_mut(&source).unwrap().deposited_m3+=v;*stock[i].entry(source).or_default()+=v;}
            }
            let outgoing=std::mem::take(&mut sediment[i]);for edge in &route.edges[i]{let j=edge.receiver as usize;for (&source,&v) in &outgoing{*sediment[j].entry(source).or_default()+=edge.weight*v;}}
        }
        for i in 0..n{if domain.posts[i].boundary.is_some(){for (&source,&v) in &sediment[i]{ledger.get_mut(&source).unwrap().exported_m3+=v;}}else if !sediment[i].is_empty(){return Err("flank3 sediment remained outside a supported boundary".into());}}
        for _ in 0..8{let mut moved=0.0;
            for parity in 0..2{for i in 0..n{let p=domain.posts[i].id;if domain.posts[i].boundary.is_some()||(p.gx+p.gy).rem_euclid(2)!=parity{continue;}
                for (ox,oy) in [(1,0),(0,1),(1,1),(-1,1)]{
                    let Some(j)=domain.index(crate::t1b_surface_domain::PostId{gy:p.gy+oy,gx:p.gx+ox})else{continue};if domain.posts[j].boundary.is_some(){continue;}
                    // Bank exchange is allowed only across a mutually open
                    // physical edge. It cannot bypass a wall or locked conduit.
                    let k=crate::t1b_surface_domain::D8.iter().position(|&(x,y)|x==ox&&y==oy).ok_or("bank edge missing from D8")?;
                    if domain.edge(i,k)!=Some(j)||domain.edge(j,(k+4)%8)!=Some(i){continue;}
                    let(hi,lo)=if z[i]>z[j]{(i,j)}else{(j,i)};let support_hi=domain.posts[hi].support_m2;let support_lo=domain.posts[lo].support_m2;
                    let sc=0.7+0.5*proxy.material_m[hi]/(proxy.material_m[hi]+20.0);let len=dx*(ox as f64).hypot(oy as f64);
                    let target=((z[hi]-z[lo]-sc*len)/(1.0/support_hi+1.0/support_lo)).max(0.0).min((z[hi]-proxy.floor_m[hi]).max(0.0)*support_hi);
                    let transfer=flank3_take(&mut stock[hi],target);let volume:f64=transfer.values().sum();if volume==0.0{continue;}z[hi]-=volume/support_hi;z[lo]+=volume/support_lo;moved+=volume;
                    for(source,v)in transfer{ledger.get_mut(&source).unwrap().bank_transferred_m3+=v;*stock[lo].entry(source).or_default()+=v;}
                }
            }}if moved<1e-6{break;}
        }
        for l in ledger.values_mut(){l.remaining_stock_m3=0.0;}
        for map in &stock{for(&s,&v)in map{if !v.is_finite()||v<0.0{return Err("flank3 invalid remaining source stock".into());}ledger.get_mut(&s).unwrap().remaining_stock_m3+=v;}}
        let mut sed_abs=0.0f64;let mut sed_rel=0.0f64;
        for l in ledger.values_mut(){if [l.initial_stock_m3,l.remaining_stock_m3,l.exported_m3,l.stream_detached_m3,l.deposited_m3,l.bank_transferred_m3].iter().any(|x|!x.is_finite()){return Err("flank3 nonfinite per-source ledger".into());}l.absolute_residual_m3=(l.initial_stock_m3-l.remaining_stock_m3-l.exported_m3).abs();l.relative_residual=l.absolute_residual_m3/l.initial_stock_m3.max(1.0);sed_abs=sed_abs.max(l.absolute_residual_m3);sed_rel=sed_rel.max(l.relative_residual);}
        // Qualify the changed geometry before accepting this revision. The
        // accepted route is reused as the next step's input, including at the
        // final physical step before the separate millimetre projection.
        let input_route_digest=route.digest.clone();
        domain=flank3_revision(domain,&z)?;
        route=crate::t1b_ls::route_detail_surface(&domain,Mode::MultipleFlowSlope11V1).map_err(|e|e.to_string())?;
        let registered=crate::t1b_ls::register_detail_channels(&domain,&route,runoff_reference_m_s,threshold_m3_s).map_err(|e|e.to_string())?;
        let audit=crate::t1b_surface_channels::audit(&domain,&route,&registered).map_err(|e|e.to_string())?;
        let water_abs=(audit.effective_source_m2-audit.effective_export_m2).abs();let water_rel=water_abs/audit.effective_source_m2.max(1.0);
        if !withheld.is_finite()||!sed_rel.is_finite()||!water_rel.is_finite()||sed_rel>1e-9||water_rel>1e-9||z.iter().any(|z|!z.is_finite()){return Err(format!("flank3 rejected step {step}: water={water_rel}, per-source sediment={sed_rel}"));}
        let row=Flank3StepReport{step:step+1,input_route_digest,accepted_domain_digest:domain.digest().into(),route_digest:route.digest.clone(),registry_digest:registered.digest,water_source_m2:audit.effective_source_m2,water_export_m2:audit.effective_export_m2,water_absolute_residual_m2:water_abs,water_relative_residual:water_rel,sediment_absolute_residual_m3:sed_abs,sediment_relative_residual:sed_rel,withheld_detachment_m3:withheld};
        on_step(&row)?;reports.push(row);max_water=max_water.max(water_rel);max_sediment=max_sediment.max(sed_rel);
    }
    let exported:f64=ledger.values().map(|l|l.exported_m3).sum();let detached:f64=ledger.values().map(|l|l.stream_detached_m3).sum();
    let physical_change:f64=(0..n).map(|i|(z[i]-z0[i])*domain.posts[i].support_m2).sum();let geom_residual=(physical_change+exported).abs();let processed:f64=ledger.values().map(|l|l.stream_detached_m3+l.deposited_m3+l.bank_transferred_m3).sum();
    if !processed.is_finite()||!geom_residual.is_finite()||geom_residual>1e-9*processed.max(exported).max(1.0){return Err(format!("flank3 geometric sediment ledger fails: {geom_residual} m3"));}
    let published_height_mm:Vec<i32>=z.iter().map(|v|i32::try_from(vz_tilekey::chart::quantize_mm(*v*1000.0)).map_err(|_|"flank3 projected height overflow".to_string())).collect::<R<_>>()?;
    let projected:Vec<f64>=published_height_mm.iter().map(|v|*v as f64/1000.0).collect();
    // Quantization is representation error, explicitly retained rather than
    // claimed as physical material creation. Fixed boundaries must already be
    // native millimetre inputs for final publication.
    for i in 0..n{if domain.posts[i].boundary.is_some()&&projected[i]!=z[i]{return Err("flank3 needs initially mm-quantized immutable boundaries".into());}}
    let physical_revision=domain.digest.clone();
    let projection:f64=(0..n).map(|i|(projected[i]-z[i])*domain.posts[i].support_m2).sum();
    domain.pin.stage="flank3-mm-projection".into();domain=flank3_revision(domain,&projected)?;
    let route=crate::t1b_ls::route_detail_surface(&domain,Mode::MultipleFlowSlope11V1).map_err(|e|e.to_string())?;
    let registry=crate::t1b_ls::register_detail_channels(&domain,&route,runoff_reference_m_s,threshold_m3_s).map_err(|e|e.to_string())?;
    crate::t1b_surface_channels::audit(&domain,&route,&registry).map_err(|e|e.to_string())?;
    let report=serde_json::json!({"rule":FLANK3_OWNER_RULE,"process_input_digest":process_input,"kernel_digest":kernel,"duration_yr":30000,"steps":steps,"initial_domain":initial.digest(),"physical_pre_projection_revision":physical_revision,"final_projected_domain":domain.digest(),"final_route":route.digest,"final_registry":registry.digest,"max_water_relative_residual":max_water,"max_per_source_sediment_relative_residual":max_sediment,"physical_net_volume_change_m3":physical_change,"exported_m3":exported,"stream_detached_m3":detached,"physical_geometry_absolute_residual_m3":geom_residual,"processed_volume_m3":processed,"physical_geometry_relative_residual":geom_residual/processed.max(exported).max(1.0),"projection_representation_delta_m3":projection,"projection_note":"explicit representation discrepancy, not water or sediment supplied","geology_stock_qualified":false});
    let native_graph=crate::t1b_surface_channels::native_graph_product(&domain,&route,&registry).map_err(|e|e.to_string())?;
    Ok(Flank3Output{physical_ground_m:z,published_height_mm,final_domain:domain,route,registry,native_graph,remaining_source_stock_m3:stock,sources:ledger.into_values().collect(),steps:reports,report})
}

/// Actual native initial-condition proxy fields for the owner process. Values
/// retain arc2's declared inference; missing native metadata is never defaulted.
/// This adapter is not a qualification of geology contacts/strength/exposure.
pub struct Flank3NativeProxy {
    pub stock_depth_m:Vec<f64>,pub material_m:Vec<f64>,pub exposure_fraction:Vec<f64>,pub floor_m:Vec<f64>,
    pub source_id:Vec<u32>,pub source_names:BTreeMap<u32,String>,pub input_digest:String,pub report:serde_json::Value,
}
impl Flank3NativeProxy {
    pub fn view(&self)->Flank3Proxy<'_>{Flank3Proxy{stock_depth_m:&self.stock_depth_m,material_m:&self.material_m,exposure_fraction:&self.exposure_fraction,floor_m:&self.floor_m,source_id:&self.source_id,source_names:&self.source_names}}
}
pub fn flank3_native_proxy(lat:&EvalLat,z:&[f64],geometry_fixed:&[u8],vid:&[u8],style:&[u8],unit:&[u16],volcanoes:&[serde_json::Value],source_surface_digest:&str)->R<Flank3NativeProxy>{
    let n=lat.nodes();if lat.dx_mm!=crate::t1b::T1B_SPACING_MM||[z.len(),geometry_fixed.len(),vid.len(),style.len(),unit.len()].iter().any(|&m|m!=n)||source_surface_digest.is_empty(){return Err("native flank proxy requires complete pinned32m fields".into());}
    let field=|v:&serde_json::Value,key:&str|->R<f64>{v[key].as_f64().filter(|x|x.is_finite()).ok_or_else(||format!("native proxy missing finite {key}"))};
    let mut stock=vec![0.0;n];let mut material=vec![0.0f64;n];let mut exposure=vec![1.0;n];let mut floor=z.to_vec();let mut ids=vec![u32::MAX;n];let mut names=BTreeMap::new();let mut dry=0;let mut positive=0;let mut edge_stock=0;let mut total=0.0;
    for i in 0..n{
        if !z[i].is_finite(){return Err("native proxy nonfinite ground".into());}
        if geometry_fixed[i]!=FREE||vid[i]==255{continue;}
        let v=volcanoes.get(vid[i]as usize).ok_or("native proxy volcano index")?;let kind=v["kind"].as_str().ok_or("native proxy volcano kind")?;
        if family_of(kind)!=Some(Family::Composite){continue;}
        if style[i]==ST_FLOOR{dry+=1;continue;}
        let id=v["id"].as_str().ok_or("native proxy volcano identity")?;let(gx,gy)=lat.g(i);let vx=v["vent_xy_m"][0].as_f64().ok_or("native proxy vent x")?;let vy=v["vent_xy_m"][1].as_f64().ok_or("native proxy vent y")?;
        let r=(gx as f64*32.0-vx).hypot(gy as f64*32.0-vy);
        let depth=flank2_profile(field(v,"h_m")?,field(v,"denudation_m")?,field(v,"s_top")?,field(v,"lambda_m")?,field(v,"crater_radius_m")?,field(v,"r_main_m")?,r);
        if !depth.is_finite()||depth<0.0{return Err("native proxy invalid inferred stock".into());}if depth==0.0{continue;}
        if unit[i]==UNIT_NONE{return Err("positive native volcanic proxy has unresolved unit identity".into());}
        let name=if unit[i]==UNIT_BODY{format!("legacy-native-proxy:{id}:body")}else{
            let u=v["ve"]["units"].as_array().and_then(|a|a.iter().find(|u|u["id"].as_u64()==Some(unit[i]as u64))).ok_or("positive native proxy unit is missing from actual history")?;
            let age=field(u,"age_yr")?;if age<0.0{return Err("native unit age is negative".into());}exposure[i]=(age/30000.0).clamp(0.0,1.0);format!("legacy-native-proxy:{id}:unit-{}",unit[i])
        };
        let hash=blake3::hash(name.as_bytes());let source=u32::from_le_bytes(hash.as_bytes()[..4].try_into().unwrap());if source==u32::MAX{return Err("native proxy source identifier collides with unresolved sentinel".into());}
        if let Some(old)=names.insert(source,name.clone()){if old!=name{return Err("native proxy source identifier collision".into());}}
        ids[i]=source;stock[i]=depth;floor[i]=z[i]-depth;positive+=1;total+=depth*1024.0;
        if i%lat.n==0||i/lat.n==0||i%lat.n+1==lat.n||i/lat.n+1==lat.n{edge_stock+=1;}
        if matches!(style[i],ST_LAVA|ST_FRONT|ST_DOME|ST_CINDER){material[i]=30.0;}
    }
    // Retained stationary buried-path proxy, now bound to actual native history
    // bytes and source lineage. No synthetic alternating layers or rock priors.
    for v in volcanoes{let Some(history)=v.pointer("/ve/flank2/history").and_then(serde_json::Value::as_array)else{continue};
        let h=field(v,"h_m")?;if h<=0.0{return Err("native proxy invalid edifice height".into());}let remain=(1.0-field(v,"denudation_m")?/h).clamp(0.0,1.0);
        for u in history{let g=&u["geom"];let Some(points)=g["pts"].as_array()else{continue};let w=field(g,"w0_m")?;let t=field(g,"t_m")?;if w<=0.0||t<=0.0{return Err("native buried path has invalid thickness/width".into());}
            for pair in points.windows(2){let a=[pair[0][0].as_f64().ok_or("path x")?,pair[0][1].as_f64().ok_or("path y")?];let b=[pair[1][0].as_f64().ok_or("path x")?,pair[1][1].as_f64().ok_or("path y")?];
                let x0=((a[0].min(b[0])-w)/32.0-lat.gx0 as f64).floor().max(0.0)as usize;let x1=((a[0].max(b[0])+w)/32.0-lat.gx0 as f64).ceil().min((lat.n-1)as f64).max(0.0)as usize;
                let y0=((a[1].min(b[1])-w)/32.0-lat.gy0 as f64).floor().max(0.0)as usize;let y1=((a[1].max(b[1])+w)/32.0-lat.gy0 as f64).ceil().min((lat.n-1)as f64).max(0.0)as usize;
                if x0>x1||y0>y1||x0>=lat.n||y0>=lat.n{continue;}
                for y in y0..=y1{for x in x0..=x1{let i=y*lat.n+x;if stock[i]<=0.0{continue;}let(gx,gy)=lat.g(i);let point=[gx as f64*32.0,gy as f64*32.0];let q=seg_nearest(point,a,b);let d=(point[0]-q[0]).hypot(point[1]-q[1]);if d<w{material[i]=material[i].max(t*remain*(1.0-(d/w).powi(2)).powi(2));}}}
            }
        }
    }
    // Native T1a bank constraint is a lower bound on removal, never an
    // unaccounted post-process raise. A missing neighbour is reported in scope.
    for i in 0..n{if stock[i]<=0.0{continue;}let x=i%lat.n;let y=i/lat.n;
        for (ox,oy)in crate::t1b_lem::D8N{let xx=x as i64+ox;let yy=y as i64+oy;if xx<0||yy<0||xx>=lat.n as i64||yy>=lat.n as i64{continue;}let j=yy as usize*lat.n+xx as usize;if geometry_fixed[j]==CHANNEL{floor[i]=floor[i].max(z[j]+0.01);}}
        if floor[i]>z[i]+1e-9{return Err("native proxy input already violates the immutable T1a bank floor".into());}
    }
    if !total.is_finite(){return Err("native proxy nonfinite total stock".into());}
    let input_digest=crate::t1b_surface_domain::content_digest(&("native_flank_proxy_binding_v1",source_surface_digest,(lat.gx0,lat.gy0,lat.n),(volcanoes,z,geometry_fixed,vid,style,unit),(&stock,&material,&exposure,&floor,&ids,&names))).map_err(|e|e.to_string())?;
    let report=serde_json::json!({"rule":"native_flank_proxy_binding_v1","source_surface":source_surface_digest,"positive_mobile_posts":positive,"dry_floor_posts":dry,"positive_stock_on_reference_edge":edge_stock,"mobile_proxy_volume_m3":total,"sources":names.len(),"geology_stock_qualified":false,"stock_support_contained":edge_stock==0,"routing_closure_not_established_by_this_adapter":true,"note":"Actual native mapped units and retained inferred-body/path law. Positive edge stock requires more native coverage, not a perimeter outlet or zeroing the stock. Exposure is the legacy age proxy, not geological exposure truth."});
    Ok(Flank3NativeProxy{stock_depth_m:stock,material_m:material,exposure_fraction:exposure,floor_m:floor,source_id:ids,source_names:names,input_digest,report})
}


/// Content-bound native input for the physical flank process. A publication
/// role (including a dam-lake display label) never creates a physical terminal.
pub struct Flank3NativeWindow {
    pub native_context_pin:String,
    pub lat:EvalLat,pub ground_m:Vec<f64>,pub geometry_fixed:Vec<u8>,
    pub own_source_m2:Vec<f64>,pub native_link:Vec<u32>,pub runoff_reference_m_s:f64,pub threshold_m3_s:f64,pub proxy:Flank3NativeProxy,pub report:serde_json::Value,
}
pub fn flank3_load_native_window(path:&std::path::Path)->R<Flank3NativeWindow>{
    use crate::t1b_surface_domain::content_digest;
    let bytes=flank2_bytes(&path.join("window.json"))?;
    let header:serde_json::Value=serde_json::from_slice(&bytes).map_err(|e|e.to_string())?;
    let h=&header["header"];
    if header["digest"]!=pf_contracts::artifact::digest_canonical(h).map_err(|e|e.to_string())?{return Err("native process header hash mismatch".into());}
    if h["schema"]!="vz.surface-owner-construction/v1"{return Err("native process needs role-aware owner construction".into());}
    let l=&h["lattice"];let lat=EvalLat{gx0:l["gx0"].as_i64().ok_or("native process gx0")?,gy0:l["gy0"].as_i64().ok_or("native process gy0")?,n:l["n"].as_u64().ok_or("native process n")?as usize,dx_mm:32000};
    if lat.n!=641||l["halo_posts"]!=64{return Err("native process requires native641-post reference windows".into());}
    let field=|name:&str,size:usize|->R<Vec<u8>>{
        let f=h["fields"].as_array().and_then(|a|a.iter().find(|f|f["name"]==name)).ok_or_else(||format!("native process missing {name}"))?;
        let file=f["file"].as_str().ok_or("native process field file")?;
        if std::path::Path::new(file).components().any(|c|!matches!(c,std::path::Component::Normal(_))){return Err("native process invalid field locator".into());}
        let b=flank2_bytes(&path.join(file))?;
        if b.len()!=lat.nodes()*size||f["digest"]!=format!("blake3:{}",blake3::hash(&b).to_hex()){return Err(format!("native process {name} size/hash mismatch"));}Ok(b)
    };
    let ground_m:Vec<f64>=field("height_mm",4)?.chunks_exact(4).map(|b|i32::from_le_bytes(b.try_into().unwrap())as f64/1000.0).collect();
    let publication_fixed=field("fixed",1)?;let geometry_fixed=field("process_geometry_fixed",1)?;let vid=field("volcano_id",1)?;let style=field("process_edifice_style",1)?;
    let native_link:Vec<u32>=field("process_native_link",4)?.chunks_exact(4).map(|b|u32::from_le_bytes(b.try_into().unwrap())).collect();
    let unit:Vec<u16>=field("process_edifice_unit",2)?.chunks_exact(2).map(|b|u16::from_le_bytes(b.try_into().unwrap())).collect();
    let weight:Vec<f64>=field("runoff_weight",4)?.chunks_exact(4).map(|b|f32::from_le_bytes(b.try_into().unwrap())as f64).collect();
    if weight.iter().any(|x|!x.is_finite()||*x<0.0){return Err("native process invalid own runoff weight".into());}
    let pin=h.pointer("/prep/net/surface_area_owner/digest").and_then(|v|v.as_str()).ok_or("native process missing area owner")?;
    if pin.len()!=64||!pin.bytes().all(|b|b.is_ascii_hexdigit()){return Err("native process invalid area owner digest".into());}
    let root=path.parent().ok_or("native process reference parent")?;
    let mut out:crate::t1b_surface_area::Output=serde_json::from_slice(&std::fs::read(root.join("surface-area").join(format!("{pin}.json"))).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
    if out.digest!=pin{return Err("native process wrong owner artifact".into());}
    out.digest.clear();if content_digest(&out).map_err(|e|e.to_string())?!=pin{return Err("native process changed owner artifact".into());}
    let cert=&out.certificate;
    if !cert.pass||cert.unresolved_requested_posts!=0||!cert.missing.is_empty()||serde_json::to_value(cert).map_err(|e|e.to_string())?!=h["prep"]["net"]["surface_area_owner"]["certificate"]{return Err("native process unqualified/mismatched area proof".into());}
    let reg=out.registry.as_ref().ok_or("native process needs registered native input")?;
    let mut regcopy=reg.clone();regcopy.digest.clear();
    if content_digest(&regcopy).map_err(|e|e.to_string())?!=reg.digest||reg.source_pin!=cert.source_pin||reg.route_digest!=cert.route_digest||h["prep"]["net"]["surface_channel_registry"]["digest"]!=reg.digest{return Err("native process stale registration binding".into());}
    let upstream=h["recipe"]["upstream"].as_array().ok_or("native process upstream pins")?;
    if !upstream.iter().any(|v|v[0]=="ls.canonical_area_candidate"&&v[1]==pin)||!upstream.iter().any(|v|v[0]=="ls.native_registered_channels_candidate"&&v[1]==reg.digest){return Err("native process missing recipe pins".into());}
    let registered:std::collections::BTreeSet<_>=reg.channels.iter().map(|c|c.id).collect();
    let mut present=vec![false;lat.nodes()];let mut guides=0;let mut lakes=0;
    for p in &out.posts{
        let i=lat.idx(p.id.gx,p.id.gy).ok_or("native process owner post outside window")?;
        if present[i]||!p.area_m2.is_finite()||p.area_m2<0.0||p.receiver.is_none(){return Err("native process duplicate or invalid owner source".into());}present[i]=true;
        if geometry_fixed[i]!=FREE{return Err("native process stored geometry conflicts with owner free post".into());}
        match publication_fixed[i]{
            FREE=>{},crate::t1b::GUIDE if registered.contains(&p.id)=>{guides+=1;},
            crate::t1b::LAKE=>{lakes+=1;},
            _=>return Err("native process role conflicts with authoritative free geometry".into()),
        }
    }
    if geometry_fixed.iter().enumerate().any(|(i,&f)|f==FREE&&!present[i]){return Err("native process missing free-post owner coverage".into());}
    let own_source_m2:Vec<_>=geometry_fixed.iter().enumerate().map(|(i,&f)|if f==FREE{1024.0*weight[i]}else{0.0}).collect();
    let volcanoes=h.pointer("/prep/net/volcanoes").and_then(|v|v.as_array()).ok_or("native process volcano history")?;
    let input_pin=header["digest"].as_str().ok_or("native process header digest")?;
    let proxy=flank3_native_proxy(&lat,&ground_m,&geometry_fixed,&vid,&style,&unit,volcanoes,input_pin)?;
    let report=serde_json::json!({"header_digest":input_pin,"owner_digest":pin,"registry_digest":reg.digest,"restored_free_guide_roles":guides,"restored_free_lake_roles":lakes,"own_source_m2":own_source_m2.iter().sum::<f64>(),"proxy":proxy.report,"scope":"content-bound native process input; construction routing proof does not qualify a changed geometry or a closed material domain"});
    Ok(Flank3NativeWindow{native_context_pin:cert.source_pin.clone(),lat,ground_m,geometry_fixed,own_source_m2,native_link,runoff_reference_m_s:reg.runoff_m_per_yr/crate::SECONDS_PER_YEAR,threshold_m3_s:reg.threshold_m3_s,proxy,report})
}

/// Closure census for mutable native terrain. Full D8 connectivity through
/// physical free posts is a sufficient support condition for future donor,
/// capture, flat, bank-transfer and sediment changes. It is conservative, not
/// a claim every fetched post will influence the requested erosion. No current
/// receiver mask or coarse catchment label freezes the future watershed.
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Flank3SupportReport {
    pub rule:String,pub native_context_pin:String,pub input_windows:usize,pub unique_posts:usize,pub duplicate_reads:usize,
    pub positive_stock_posts:usize,pub visited_free_posts:usize,pub touched_fixed_posts:usize,
    pub positive_window_edge_posts:usize,pub positive_stock_frontier_posts:usize,
    pub closed_positive_components:Vec<(crate::t1b_surface_domain::PostId,usize,usize)>,pub source_names:BTreeMap<u32,String>,
    pub missing:Vec<crate::t1b_surface_domain::PostId>,pub required_windows:Vec<(i64,i64)>,
    pub input_digests:Vec<String>,pub complete:bool,pub digest:String,
}
#[derive(Clone,Copy,PartialEq)]
struct Flank3SupportPost {height:f64,fixed:u8,link:u32,own:f64,stock:f64,material:f64,exposure:f64,source:u32}
pub fn flank3_native_support(paths:&[std::path::PathBuf],max_posts:usize)->R<Flank3SupportReport>{
    use crate::t1b_surface_domain::{PostId,D8,content_digest};
    use std::collections::{BTreeMap,BTreeSet};
    let mut context:Option<String>=None;
    let mut posts=BTreeMap::<PostId,Flank3SupportPost>::new();let mut names=BTreeMap::new();let mut digests=Vec::new();let mut duplicates=0;let mut edge_positive=BTreeSet::new();
    for path in paths {
        let w=flank3_load_native_window(path)?;digests.push(w.proxy.input_digest.clone());
        if let Some(pin)=&context{if *pin!=w.native_context_pin{return Err("mutable support mixes different native contexts".into());}}else{context=Some(w.native_context_pin.clone());}
        if w.geometry_fixed.iter().any(|v|!matches!(*v,crate::t1b::FREE|crate::t1b::CHANNEL|crate::t1b::GUIDE|crate::t1b::MARINE|crate::t1b::LAKE)){return Err("unsupported physical support boundary (crop rings are not terminals)".into());}
        for(&id,name)in &w.proxy.source_names{if let Some(old)=names.insert(id,name.clone()){if old!=*name{return Err("support source identity collision".into());}}}
        for i in 0..w.lat.nodes(){let(gx,gy)=w.lat.g(i);let id=PostId{gx,gy};let p=Flank3SupportPost{height:w.ground_m[i],fixed:w.geometry_fixed[i],link:w.native_link[i],own:w.own_source_m2[i],stock:w.proxy.stock_depth_m[i],material:w.proxy.material_m[i],exposure:w.proxy.exposure_fraction[i],source:w.proxy.source_id[i]};
            if let Some(old)=posts.get(&id){if *old!=p{return Err(format!("native support overlap differs at {id:?}"));}duplicates+=1;}else{
                if posts.len()>=max_posts{return Err(format!("NeedsDomain resource: native support exceeds {max_posts} unique posts while reading {} at {id:?}; no process run",path.display()));}
                posts.insert(id,p);
            }
            if p.stock>0.0&&(i%w.lat.n==0||i/w.lat.n==0||i%w.lat.n+1==w.lat.n||i/w.lat.n+1==w.lat.n){edge_positive.insert(id);}
        }
    }
    let positive=posts.values().filter(|p|p.stock>0.0).count();let mut seen=BTreeSet::new();let mut touched=BTreeSet::new();let mut missing=BTreeSet::new();let mut stack=Vec::new();let mut stock_frontier=BTreeSet::new();let mut closed=Vec::new();
    // Walk each physical component once, including its real donor and bank
    // support. The seed iteration does not enqueue millions of duplicates.
    for(&seed,p)in &posts{if p.stock<=0.0||seen.contains(&seed){continue;}stack.push(seed);let mut component_posts=0;let mut component_positive=0;let mut component_open=false;
        while let Some(id)=stack.pop(){if !seen.insert(id){continue;}
            let p=posts.get(&id).ok_or("support traversal lost an actual source")?;
            if p.fixed!=FREE{return Err("mobile stock at a fixed native boundary".into());}
            component_posts+=1;if p.stock>0.0{component_positive+=1;}
            for k in 0..D8.len(){let q=id.neighbour(k);
                let half=vz_tilekey::CHART_HALF_EXTENT_MM/crate::t1b::T1B_SPACING_MM;
                // The native finite-chart contract has no exterior land edge;
                // skipping that nonexistent edge does not create an outlet.
                if q.gx.abs()>half||q.gy.abs()>half{continue;}
                match posts.get(&q){
                None=>{missing.insert(q);component_open=true;if p.stock>0.0{stock_frontier.insert(id);}},Some(v) if v.fixed!=FREE=>{touched.insert(q);},Some(_)=>{if !seen.contains(&q){stack.push(q);}},
            }}
        }
        if !component_open {closed.push((seed,component_posts,component_positive));}
    }
    let required:BTreeSet<_>=missing.iter().map(|p|(p.gx.div_euclid(512),p.gy.div_euclid(512))).collect();digests.sort();
    let mut out=Flank3SupportReport{rule:"native_mutable_full_d8_support_v1: physical free connectivity from all supplied positive stock; real native fixed terminals stop the component; no T0 basin wall; missing posts are requests, not exports; sufficient conservative closure for mutable MFD/banks/deposition, not minimum influence".into(),native_context_pin:context.ok_or("native support needs actual windows")?,input_windows:paths.len(),unique_posts:posts.len(),duplicate_reads:duplicates,positive_stock_posts:positive,visited_free_posts:seen.len(),touched_fixed_posts:touched.len(),positive_window_edge_posts:edge_positive.len(),positive_stock_frontier_posts:stock_frontier.len(),closed_positive_components:closed,source_names:names,missing:missing.into_iter().collect(),required_windows:required.iter().copied().collect(),input_digests:digests,complete:required.is_empty(),digest:String::new()};
    out.digest=content_digest(&out).map_err(|e|e.to_string())?;Ok(out)
}

/// Compact process source from the canonical native evaluator. It does not
/// require a construction-SFD registry to sample geometry, material and own
/// runoff. It carries the actual construction proofs and no mutable-flow claim.
pub struct Flank3MaterialTile {pub header:serde_json::Value,pub fields:BTreeMap<String,Vec<u8>>}
/// The source-solid adapter lives beside this process, inside pf-region. It
/// consumes E1/E3 records; it never constructs another geological column.
#[path = "t1b_surface_material.rs"]
pub mod material;
pub fn flank3_material_tile(ctx:&crate::t1b_rivernet::RnContext<'_>,owner:&crate::t1b_surface_stages::NativeOwner<'_,'_>,tx:i64,ty:i64)->R<Flank3MaterialTile>{
    use crate::t1b_surface_domain::content_digest;
    use crate::t1b_surface_completion::Rect;
    let half=vz_tilekey::CHART_HALF_EXTENT_MM/crate::t1b::T1B_SPACING_MM;
    let x0=tx.checked_mul(512).ok_or("material tile x overflow")?;let y0=ty.checked_mul(512).ok_or("material tile y overflow")?;
    if x0 < -half||y0 < -half||x0>half||y0>half{return Err("material tile outside native finite chart".into());}
    let rect=Rect{gx0:x0,gy0:y0,width:(half+1-x0).min(512)as usize,height:(half+1-y0).min(512)as usize};
    let n=1153;let lat=EvalLat{gx0:(x0-320).clamp(-half,half-n as i64+1),gy0:(y0-320).clamp(-half,half-n as i64+1),n,dx_mm:32000};
    let verified=ctx.skel.verify_tiles_for_rect([lat.gx0*32000,lat.gy0*32000],[(lat.gx0+n as i64-1)*32000,(lat.gy0+n as i64-1)*32000])?;
    let mut lp=crate::t1b_rivernet::lat_prep(lat,&ctx.chart,ctx.t0,ctx.extras,ctx.params,ctx.strata,ctx.anchor)?;lp.publish_margin_posts=256;
    crate::t1b_esc_deps::begin_window(tx,ty,&ctx.rn_recipe_digest,&ctx.layer_digest);
    let s=crate::t1b_rivernet::evaluate_owned(&lp,&ctx.net,&ctx.curves,ctx.t0,&ctx.extras.laws,&ctx.table,ctx.params,&ctx.coast,Some(owner),None)?.surface()?;
    if [s.fixed.len(),s.vid.len(),s.ve_style.len(),s.ve_unit.len(),s.chan_link.len(),lp.rw.len()].iter().any(|&m|m!=lat.nodes()){return Err("native material source needs complete physical fields and explicit own runoff".into());}
    let proof=s.seam.as_ref().ok_or("native material source lacks construction proof")?;
    if proof.owner.is_empty(){return Err("native material source has no canonical stage ownership".into());}
    let z:Vec<f64>=s.z.iter().map(|v|{
        if !v.is_finite(){return Err("nonfinite native material ground".to_string());}
        i32::try_from(vz_tilekey::chart::quantize_mm(*v*1000.0)).map(|v|v as f64/1000.0).map_err(|_|"native material height overflows mm".into())
    }).collect::<R<_>>()?;
    let volcanoes=serde_json::to_value(&ctx.coast.volcanoes).map_err(|e|e.to_string())?;
    let source=content_digest(&("native_material_source_v1",owner.source_pin(),&ctx.rn_recipe_digest,&ctx.layer_digest,&ctx.chart.decl,&proof.owner)).map_err(|e|e.to_string())?;
    // The proxy needs one physical bank-neighbour beyond the owned core,
    // not material identities in the outer numerical evaluation margin.
    // At the actual chart edge the existing finite-graph contract is one-sided.
    let pn=rect.width.max(rect.height)+2;
    let pl=EvalLat{gx0:(x0-1).clamp(-half,half-pn as i64+1),gy0:(y0-1).clamp(-half,half-pn as i64+1),n:pn,dx_mm:32000};
    let pi:Vec<usize>=(0..pl.nodes()).map(|i|{let(gx,gy)=pl.g(i);lat.idx(gx,gy).ok_or("material bank support outside actual evaluation")}).collect::<Result<_,_>>()?;
    let pz:Vec<f64>=pi.iter().map(|&i|z[i]).collect();let pf:Vec<u8>=pi.iter().map(|&i|s.fixed[i]).collect();let pv:Vec<u8>=pi.iter().map(|&i|s.vid[i]).collect();let ps:Vec<u8>=pi.iter().map(|&i|s.ve_style[i]).collect();let pu:Vec<u16>=pi.iter().map(|&i|s.ve_unit[i]).collect();
    let proxy=flank3_native_proxy(&pl,&pz,&pf,&pv,&ps,&pu,volcanoes.as_array().ok_or("native volcano metadata")?,&source)?;
    let count=rect.count().map_err(|e|e.to_string())?;let indices:Vec<usize>=(0..count).map(|i|{let p=rect.id(i);lat.idx(p.gx,p.gy).ok_or("material source core outside evaluation")}).collect::<Result<_,_>>()?;
    let proxy_indices:Vec<usize>=(0..count).map(|i|{let p=rect.id(i);pl.idx(p.gx,p.gy).ok_or("material core outside bank support")}).collect::<Result<_,_>>()?;
    let mut fields:BTreeMap<String,Vec<u8>>=BTreeMap::new();
    fields.insert("height_mm.i32".into(),indices.iter().flat_map(|&i|((z[i]*1000.0).round()as i32).to_le_bytes()).collect());
    fields.insert("process_geometry_fixed.u8".into(),indices.iter().map(|&i|s.fixed[i]).collect());
    fields.insert("process_native_link.u32".into(),indices.iter().flat_map(|&i|s.chan_link[i].to_le_bytes()).collect());
    fields.insert("volcano_id.u8".into(),indices.iter().map(|&i|s.vid[i]).collect());
    fields.insert("process_edifice_style.u8".into(),indices.iter().map(|&i|s.ve_style[i]).collect());
    fields.insert("process_edifice_unit.u16".into(),indices.iter().flat_map(|&i|s.ve_unit[i].to_le_bytes()).collect());
    fields.insert("runoff_weight.f32".into(),indices.iter().flat_map(|&i|lp.rw[i].to_le_bytes()).collect());
    fields.insert("proxy_source_id.u32".into(),proxy_indices.iter().flat_map(|&i|proxy.source_id[i].to_le_bytes()).collect());
    for(name,v)in [("proxy_stock_depth_m.f64",&proxy.stock_depth_m),("proxy_material_m.f64",&proxy.material_m),("proxy_exposure_fraction.f64",&proxy.exposure_fraction),("proxy_floor_m.f64",&proxy.floor_m)]{
        fields.insert(name.into(),proxy_indices.iter().flat_map(|&i|v[i].to_le_bytes()).collect());
    }
    let mut basins=BTreeMap::new();for &i in &indices {let c=lp.cell[i];if !basins.contains_key(&c){basins.insert(c,ctx.t0.final_basin(c)?);}}
    fields.insert("native_t0_cell.u32".into(),indices.iter().flat_map(|&i|lp.cell[i].to_le_bytes()).collect());
    fields.insert("native_t0_basin.u32".into(),indices.iter().flat_map(|&i|basins[&lp.cell[i]].to_le_bytes()).collect());
    let descriptors:Vec<_>=fields.iter().map(|(name,b)|serde_json::json!({"file":name,"bytes":b.len(),"digest":blake3::hash(b).to_hex().to_string()})).collect();
    let mut header=serde_json::json!({"schema":"vz.native-physical-material/v1","sampling_kernel_digest":blake3::hash(include_str!("t1b_ve.rs").as_bytes()).to_hex().to_string(),"rect":rect,"spacing_mm":32000,"support_m2":1024.0,"runoff_reference_m_s":ctx.t0.runoff_m_per_yr/crate::SECONDS_PER_YEAR,"threshold_m3_s":crate::t1b_rivernet::window_params(ctx.params,0).threshold_m3_s(),"own_source_rule":"FREE: support_m2 times stored native f32 runoff_weight; physical fixed terminals: zero local own source; accumulated area is never a source","native_context_pin":owner.source_pin(),"source_digest":source,"rn_recipe_digest":ctx.rn_recipe_digest,"layer_digest":ctx.layer_digest,"chart":ctx.chart.decl,"construction_proofs":proof.owner,"verified_t1a_tiles":verified,"fields":descriptors,"source_names":proxy.source_names,"proxy_rule":proxy.report["rule"],"geology_stock_qualified":false,"mutable_support_closed":false,"note":"Canonical half-open512-post source owner, clipped only at the declared finite chart extent. Physical fixed/style/unit survive display-only lake/channel roles. Proxy quantities retain declared inferred stock/material/exposure, not geological contacts. Native T0 basin labels are observations, not routing barriers. No own source inferred from accumulated area."});
    if let Some(env)=&ctx.params.lithology.envelope {
        header["source_envelope"]=serde_json::json!({"digest":env.digest,"rule":env.rule,"ownership_rule":material::RULE,"stock_admitted":false});
    }
    let digest=content_digest(&header).map_err(|e|e.to_string())?;header["digest"]=digest.into();Ok(Flank3MaterialTile{header,fields})
}
pub fn flank3_publish_material_tile(tile:&Flank3MaterialTile,path:&std::path::Path)->R<()>{
    use std::io::Write;
    std::fs::create_dir(path).map_err(|e|e.to_string())?;
    for(name,bytes)in &tile.fields {std::fs::OpenOptions::new().write(true).create_new(true).open(path.join(name)).and_then(|mut f|f.write_all(bytes)).map_err(|e|e.to_string())?;}
    std::fs::OpenOptions::new().write(true).create_new(true).open(path.join("material.json")).and_then(|mut f|f.write_all(&serde_json::to_vec_pretty(&tile.header).unwrap())).map_err(|e|e.to_string())
}
