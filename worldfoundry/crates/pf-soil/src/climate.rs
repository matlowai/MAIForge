//! The canonical T1b climate (`seed-climate-v1`, the `climate=` key artifact) on the T0 lat/lon grid, reduced to sea
//! level the way pf-flight's climate_look_proxy_v1 does: T_sl(c) = temperature_mean_c(c) + lapse x max(max(elevation,
//! fill_level)(c) - sea, 0) (undoing the tier's own lapse cooling at the cell's elevation). At a post: bilinear between
//! T0 cell centres (all cells; the pf-region `bilinear` with land_only = false), then lapsed to the post's height.
//! The climate-coarse-v0 run climate (T0 run/climate) is refused: CLIMATE-SOURCE (2026-10-21) showed it is a different
//! model with a first-windward-cell rain spike; the T1b chain reads seed-climate-v1.

use pf_contracts::artifact::{read_artifact_header, read_field_f32, verify_artifact};
use serde_json::Value;
use std::path::Path;

pub const CLIMATE_TIER_REQUIRED: &str = "seed-climate-v1";

pub struct Climate {
    pub nlat: usize,
    pub nlon: usize,
    pub sea: f64,
    pub lapse: f64,
    pub t_sl: Vec<f64>,
    pub half: Vec<f64>,
    pub p: Vec<f64>,
    pub ai: Vec<f64>,
    pub climate_digest: String,
    pub terrain_digest: String,
    pub hydrology_digest: String,
    pub tier: String,
    /// The artifact's axial tilt when it declares one (recorded, not used: the params carry tilt_deg).
    pub axial_tilt_deg: Option<f64>,
}

fn find_key<'v>(v: &'v Value, key: &str) -> Option<&'v Value> {
    match v {
        Value::Object(m) => { if let Some(x) = m.get(key) { return Some(x); } m.values().find_map(|x| find_key(x, key)) }
        Value::Array(a) => a.iter().find_map(|x| find_key(x, key)),
        _ => None,
    }
}

impl Climate {
    /// `climate_dir`: the seed-climate-v1 artifact; `t0_run`: the T0 run dir holding `terrain/` and `hydrology/`.
    pub fn load(climate_dir: &Path, t0_run: &Path) -> Result<Climate, String> {
        let climate_digest = verify_artifact::<Value>(climate_dir).map_err(|e| format!("verify climate {}: {e}", climate_dir.display()))?;
        let tdir = t0_run.join("terrain");
        let hdir = t0_run.join("hydrology");
        let terrain_digest = verify_artifact::<Value>(&tdir).map_err(|e| format!("verify terrain: {e}"))?;
        let hydrology_digest = verify_artifact::<Value>(&hdir).map_err(|e| format!("verify hydrology: {e}"))?;
        let ch = read_artifact_header::<Value>(climate_dir).map_err(|e| format!("climate header: {e}"))?;
        let tier = ch.header.semantic.backend.name.clone();
        if tier != CLIMATE_TIER_REQUIRED { return Err(format!("climate tier {tier}; pf-soil reads {CLIMATE_TIER_REQUIRED} only (the T1b climate= key artifact); refusing")); }
        for (role, d) in [("terrain", &terrain_digest), ("hydrology", &hydrology_digest)] {
            if !ch.header.semantic.inputs.iter().any(|x| x == d) { return Err(format!("climate does not cite the T0 {role} digest {d}; refusing")); }
        }
        let shape = ch.header.semantic.fields.first().ok_or("climate has no fields")?.shape.clone();
        if shape.len() != 2 { return Err("climate field not 2-D".into()); }
        let (nlat, nlon) = (shape[0] as usize, shape[1] as usize);
        let n = nlat * nlon;
        let head: Value = serde_json::from_slice(&std::fs::read(climate_dir.join("artifact.json")).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let lapse = find_key(&head, "lapse_rate_k_per_m").and_then(|v| v.get("value")).and_then(Value::as_f64).ok_or("climate: no lapse_rate_k_per_m.value")?;
        if !(lapse > 0.0) { return Err(format!("climate lapse {lapse} not positive")); }
        let axial_tilt_deg = find_key(&head, "axial_tilt_deg").and_then(|v| v.get("value").or(Some(v))).and_then(Value::as_f64);
        let th = read_artifact_header::<Value>(&tdir).map_err(|e| e.to_string())?;
        let sea = th.header.domain.get("sea_level_m").and_then(Value::as_f64).ok_or("terrain sea_level_m")?;
        let f = |dir: &Path, name: &str| -> Result<Vec<f64>, String> {
            let v = read_field_f32(dir, name).map_err(|e| format!("{name}: {e}"))?;
            if v.len() != n { return Err(format!("{name}: {} values, expected {n}", v.len())); }
            if v.iter().any(|x| !x.is_finite()) { return Err(format!("{name}: nonfinite value")); }
            Ok(v.into_iter().map(|x| x as f64).collect())
        };
        let (tm, half, p, ai) = (f(climate_dir, "temperature_mean_c")?, f(climate_dir, "temperature_half_range_k")?, f(climate_dir, "precipitation_mm_per_yr")?, f(climate_dir, "aridity_index")?);
        let el = f(&tdir, "elevation_m")?;
        let fl = f(&hdir, "fill_level_m")?;
        let t_sl = (0..n).map(|c| tm[c] + lapse * (el[c].max(fl[c]) - sea).max(0.0)).collect();
        Ok(Climate { nlat, nlon, sea, lapse, t_sl, half, p, ai, climate_digest, terrain_digest, hydrology_digest, tier, axial_tilt_deg })
    }

    /// Bilinear between cell centres (all cells), identical in form to pf_region::t1b::bilinear(.., land_only = false).
    pub fn bilinear(&self, field: &[f64], lat_deg: f64, lon_deg: f64) -> f64 {
        let dlat = 180.0 / self.nlat as f64;
        let dlon = 360.0 / self.nlon as f64;
        let fi = (lat_deg + 90.0) / dlat - 0.5;
        let fj = (lon_deg + 180.0) / dlon - 0.5;
        let (i0, j0) = (fi.floor(), fj.floor());
        let (ti, tj) = (fi - i0, fj - j0);
        let mut acc = 0.0;
        for (di, wi) in [(0i64, 1.0 - ti), (1, ti)] {
            let ii = (i0 as i64 + di).clamp(0, self.nlat as i64 - 1) as usize;
            for (dj, wj) in [(0i64, 1.0 - tj), (1, tj)] {
                let jj = (j0 as i64 + dj).rem_euclid(self.nlon as i64) as usize;
                acc += wi * wj * field[ii * self.nlon + jj];
            }
        }
        acc
    }
}
