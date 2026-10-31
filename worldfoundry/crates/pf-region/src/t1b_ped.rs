//! Lane T1B-PEDESTAL (2026-10-21): measurement dump for the T1b-below-T0 pedestal under composite edifices.
//!
//! Env var `VZ_PEDESTAL_DUMP=<dir>` (read-only, never digested, no output bytes): every `evaluate` call (layer E-solve
//! evaluations and window lattices) writes one folder `<dir>/<gx0>_<gy0>_<n>_c<k>` (k = call counter in this process)
//! with the heights before the coast / volcano stage and after the whole evaluation, the over / fixed codes, the T0
//! bilinear, the T0 cell and the per-post E. An existing folder is left untouched and the dump skipped.

use crate::t1b_rivernet::EvalLat;
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::sync::atomic::{AtomicU32, Ordering};

static CALLS: AtomicU32 = AtomicU32::new(0);

/// The dump directory when the env var is set.
pub fn dump_dir() -> Option<std::ffi::OsString> { std::env::var_os("VZ_PEDESTAL_DUMP") }

#[allow(clippy::too_many_arguments)]
pub fn dump(dir: &std::ffi::OsStr, lat: &EvalLat, z_pre: &[f64], z_after: &[f64], over: &[u8], fixed: &[u8], bil: &[f32], cell: &[u32], e: &[f32]) {
    let k = CALLS.fetch_add(1, Ordering::SeqCst);
    let d = std::path::Path::new(dir).join(format!("{}_{}_{}_c{}", lat.gx0, lat.gy0, lat.n, k));
    if d.exists() { println!("  VZ_PEDESTAL_DUMP: {} exists, dump skipped (never overwritten)", d.display()); return; }
    let go = || -> std::io::Result<()> {
        std::fs::create_dir_all(&d)?;
        let wb = |name: &str, b: &[u8]| -> std::io::Result<()> { std::fs::OpenOptions::new().write(true).create_new(true).open(d.join(name))?.write_all(b) };
        let f64s = |v: &[f64]| -> Vec<u8> { v.iter().flat_map(|x| (*x as f32).to_le_bytes()).collect() };
        wb("z_pre.f32", &f64s(z_pre))?;
        wb("z_after.f32", &f64s(z_after))?;
        wb("over.u8", over)?;
        wb("fixed.u8", fixed)?;
        wb("bil.f32", &bil.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>())?;
        wb("cell.u32", &cell.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>())?;
        wb("e.f32", &e.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>())?;
        let j = format!("{{\"gx0\": {}, \"gy0\": {}, \"n\": {}, \"dx_mm\": {}, \"call\": {}}}", lat.gx0, lat.gy0, lat.n, lat.dx_mm, k);
        wb("lat.json", j.as_bytes())
    };
    match go() { Ok(()) => println!("  VZ_PEDESTAL_DUMP: wrote {} (measurement only)", d.display()), Err(e) => println!("  VZ_PEDESTAL_DUMP: {} failed: {e}", d.display()) }
}

pub const LABEL: &str = "EXPLORATORY lane T1B-PEDESTAL 2026-10-21; layer AND windows";
pub const SOLVE_RULE: &str = "edifice_ground_in_solve_v1: the E solve's cell means (hat-weighted and plain) count edifice posts (over = edifice) at their pre-edifice ground (the height before the coast / volcano stage) against the T0 bilinear there, instead of leaving them out of both M and T: the T0 cell mean is ground (the volcano record's summit is a projection above it, not crust), so the ground the cone stands on belongs to the cell mean; posts whose pre-edifice ground is at or below sea level are left out as before";

pub const BASE_RULE: &str = "edifice_local_base_v1: a composite cone (lane CF scaling; lane VE composite body) stands on its own pre-edifice T1b ground instead of the T0 bilinear at the vent: base G = mean over the foot ring 0.8 R_f <= r <= 1.2 R_f of the T1b ground before the coast / volcano stage (marine posts at sea level), R_f = clamp((summit - G) / q, W_min, W_max) / 2 (lane CF's scaling law, fixed point by 4 iterations from the bilinear base); the record summit is kept (H = summit - G), so the cone's flanks run down to the local ground instead of standing on a pedestal at the bilinear; the vent placement (clearance from the bilinear-based H) is unchanged; the CF sloped base keeps the T0 surface tilt (t = B(post) - B(vent)); G is measured by the layer on its final E-solve evaluation (128 m chart lattice; the solve evaluations keep the bilinear bases, as lane VE keeps CF's cone, so the E / F tables equal the key-off ones) and published (header lane_ped), windows read it (every window the same number: seams unchanged); applies to land vents (bilinear base above sea) only; floor H >= H_bilinear / 2";

pub const RELIEF_RULE: &str = "edifice_local_relief_v1 (lane VOLCANIC-SPLIT, design 2026-10-21 21:04: keep the record relief above local ground, summits ~900 m lower; option E of lane T1B-PEDESTAL; needs ped_base): as edifice_local_base_v1, except the cone keeps the record RELIEF instead of the record summit: H = H_bilinear = record subgrid summit - T0 bilinear at the vent (the relief T0 published above its smooth surface), built summit = G + H (lower than the record summit by the pedestal, bilinear - G); R_f = clamp(H / q, W_min, W_max) / 2 is the bilinear-based foot (so the foot ring that measures G is the ring the cone actually stands on: G is the ring mean in one pass, no fixed point, no floor); the vent placement and the CF sloped base (t = B(post) - B(vent)) are unchanged; the record summit stays in the record (VolcanoSite.summit_record_m) and the built summit is published with G (lane_ped.ground[].built_summit_m)";

/// Lane T1B-PEDESTAL params (opt-in; None = unchanged).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct PedestalParamsV1 {
    pub on: u32,
    /// 1 = edifice_ground_in_solve_v1 (measurement experiment).
    pub solve_ground: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solve_rule: Option<String>,
    /// 1 = edifice_local_base_v1.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub base: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub base_rule: Option<String>,
    /// Lane VOLCANIC-SPLIT: 1 = edifice_local_relief_v1 (with base 1).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub relief: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relief_rule: Option<String>,
    pub label: String,
}

fn is_zero(v: &u32) -> bool { *v == 0 }

impl PedestalParamsV1 {
    pub fn v1() -> Self { PedestalParamsV1 { on: 1, solve_ground: 0, solve_rule: None, base: 0, base_rule: None, relief: 0, relief_rule: None, label: LABEL.into() } }
    pub fn with_relief(mut self, on: bool) -> Self { self.relief = on as u32; self.relief_rule = if on { Some(RELIEF_RULE.into()) } else { None }; self }
    pub fn with_solve(mut self, on: bool) -> Self { self.solve_ground = on as u32; self.solve_rule = if on { Some(SOLVE_RULE.into()) } else { None }; self }
    pub fn with_base(mut self, on: bool) -> Self { self.base = on as u32; self.base_rule = if on { Some(BASE_RULE.into()) } else { None }; self }
    pub fn any(&self) -> bool { self.solve_ground == 1 || self.base == 1 }
    pub fn validate(&self) -> Result<(), String> {
        if self.on != 1 { return Err("ped: on must be 1 (absent = off)".into()); }
        if self.solve_ground > 1 || self.base > 1 || self.relief > 1 { return Err("ped_solve / ped_base / ped_relief must be 0 or 1".into()); }
        if self.relief == 1 && self.base != 1 { return Err("ped_relief=1 needs ped_base=1".into()); }
        if !self.any() { return Err("ped: no rule selected (ped_solve=1 or ped_base=1)".into()); }
        Ok(())
    }
}

/// True when the E solve counts edifice posts at their pre-edifice ground.
pub fn solve_ground(p: &crate::t1b_rivernet::RiverNetParamsV1) -> bool { p.pedestal.as_ref().map(|q| q.solve_ground == 1).unwrap_or(false) }
/// True when composite cones stand on the published local ground (edifice_local_base_v1).
pub fn local_base(p: &crate::t1b_rivernet::RiverNetParamsV1) -> bool { p.pedestal.as_ref().map(|q| q.base == 1).unwrap_or(false) }
/// True when the eligible cones keep the record relief over G (edifice_local_relief_v1).
pub fn local_relief(p: &crate::t1b_rivernet::RiverNetParamsV1) -> bool { p.pedestal.as_ref().map(|q| q.base == 1 && q.relief == 1).unwrap_or(false) }
/// True when evaluate must keep the heights before the coast / volcano stage.
pub fn needs_ground(p: &crate::t1b_rivernet::RiverNetParamsV1) -> bool { solve_ground(p) || local_base(p) }

/// Whether edifice_local_base_v1 applies to a record of this kind with this bilinear base (the cone must be a lane CF
/// scaled cone; with lane VE only the composite family).
pub fn base_applies(kind: &str, base_bil: f64, sea: f64, cf: bool, ve: bool) -> bool {
    cf && base_bil > sea && (!ve || crate::t1b_ve::family_of(kind) == Some(crate::t1b_ve::Family::Composite))
}

/// Per-volcano local ground record (published in the layer header lane_ped.ground).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalGround { pub id: String, pub ground_m: f64, pub base_bilinear_m: f64, pub summit_m: f64, pub h_bilinear_m: f64, pub h_m: f64, pub r_foot_m: f64, pub ring_posts: u64, pub ring_marine_posts: u64, pub floored: bool,
    /// Lane VOLCANIC-SPLIT (ped_relief only): built summit G + H.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub built_summit_m: Option<f64> }

/// edifice_local_base_v1: G per eligible volcano from `ground` (heights before the coast / volcano stage) on `lat`.
/// `sites` are the volcano sites as built (their base_m is the bilinear base when built without a ground map, else
/// base_m + ped_shift_m is).
pub fn local_ground(lat: &EvalLat, ground: &[f64], sea: f64, sites: &[crate::t1b_coast::VolcanoSite], cf: &crate::t1b_coast::CfParamsV1, ve: bool) -> Vec<LocalGround> {
    local_ground2(lat, ground, sea, sites, cf, ve, false)
}

/// As `local_ground`; `relief` = edifice_local_relief_v1 (H = the bilinear-based relief, R_f from it, G one ring mean).
pub fn local_ground2(lat: &EvalLat, ground: &[f64], sea: f64, sites: &[crate::t1b_coast::VolcanoSite], cf: &crate::t1b_coast::CfParamsV1, ve: bool, relief: bool) -> Vec<LocalGround> {
    let (n, dx) = (lat.n, lat.dx());
    let q = cf.cone_hw_permille as f64 / 1000.0;
    let (wmin, wmax) = (cf.cone_wb_min_mm as f64 / 1000.0, cf.cone_wb_max_mm as f64 / 1000.0);
    let mut out = Vec::new();
    for v in sites {
        let base_bil = v.base_m + v.ped_shift_m.unwrap_or(0.0);
        if v.cf.is_none() || q <= 0.0 || !base_applies(&v.kind, base_bil, sea, true, ve) { continue; }
        let summit = v.summit_record_m;
        let h0 = summit - base_bil;
        let ring = |g: f64| -> (f64, u64, u64, f64) {
            let rf = (if relief { h0 } else { summit - g } / q).clamp(wmin, wmax) / 2.0;
            let (r0, r1) = (0.8 * rf, 1.2 * rf);
            let (gx0, gy0) = (((v.vent_xy_m[0] - r1) / dx).floor() as i64, ((v.vent_xy_m[1] - r1) / dx).floor() as i64);
            let (gx1, gy1) = (((v.vent_xy_m[0] + r1) / dx).ceil() as i64, ((v.vent_xy_m[1] + r1) / dx).ceil() as i64);
            let (mut sum, mut cnt, mut mar) = (0.0, 0u64, 0u64);
            for gy in gy0.max(lat.gy0)..=gy1.min(lat.gy0 + n as i64 - 1) {
                for gx in gx0.max(lat.gx0)..=gx1.min(lat.gx0 + n as i64 - 1) {
                    let r = (gx as f64 * dx - v.vent_xy_m[0]).hypot(gy as f64 * dx - v.vent_xy_m[1]);
                    if r < r0 || r > r1 { continue; }
                    let i = (gy - lat.gy0) as usize * n + (gx - lat.gx0) as usize;
                    let z = ground[i];
                    if !z.is_finite() { continue; }
                    if z <= sea { mar += 1; }
                    sum += z.max(sea); cnt += 1;
                }
            }
            (if cnt > 0 { sum / cnt as f64 } else { f64::NAN }, cnt, mar, rf)
        };
        let mut g = base_bil;
        let mut last = (f64::NAN, 0u64, 0u64, 0.0);
        for _ in 0..if relief { 1 } else { 4 } { last = ring(g); if !last.0.is_finite() { break; } g = last.0; }
        if !last.0.is_finite() { continue; }
        if relief {
            let rf = (h0 / q).clamp(wmin, wmax) / 2.0;
            out.push(LocalGround { id: v.id.clone(), ground_m: g, base_bilinear_m: base_bil, summit_m: summit, h_bilinear_m: h0, h_m: h0, r_foot_m: rf, ring_posts: last.1, ring_marine_posts: last.2, floored: false, built_summit_m: Some(g + h0) });
            continue;
        }
        let floored = summit - g < 0.5 * h0;
        if floored { g = summit - 0.5 * h0; }
        let rf = ((summit - g) / q).clamp(wmin, wmax) / 2.0;
        out.push(LocalGround { id: v.id.clone(), ground_m: g, base_bilinear_m: base_bil, summit_m: summit, h_bilinear_m: h0, h_m: summit - g, r_foot_m: rf, ring_posts: last.1, ring_marine_posts: last.2, floored, built_summit_m: None });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::t1b_coast::{CfParamsV1, CfVolcano, VolcanoSite};
    fn site(base_bil: f64, summit: f64) -> VolcanoSite {
        VolcanoSite { id: "v".into(), kind: "continental_arc".into(), shape: "stratovolcano".into(), record_xy_m: [0.0, 0.0], vent_xy_m: [0.0, 0.0], moved_m: 0.0, record_cell: 0, vent_cell: 0,
            summit_record_m: summit, base_m: base_bil, h_m: summit - base_bil, lambda_m: 3000.0, r_main_m: 8000.0, r_influence_m: 10000.0, s_top: 0.65, s_base: 0.09,
            crater_radius_m: 350.0, crater_depth_m: 150.0, age_my: 0.0, e_vent_m_per_my: 0.0, denudation_m: 0.0, stage: String::new(), breached: false, breach_dir_rad: 0.0,
            gully_dirs_rad: vec![], dike_dirs_rad: vec![], neck_relief_m: 0.0, lobe_phase_rad: vec![], lobe_norm: 1.0,
            cf: Some(CfVolcano { wb_m: 16000.0, r_foot_m: 8000.0, q_built: 0.15, lambda_m: 3000.0, s_base: 0.1, r_end_max_m: 10000.0, clearance_m: 8000.0, flows: vec![], gullies: vec![], vent_to_t1a_m: 8000.0 }),
            ve: None, ped_shift_m: None, built_summit_m: None }
    }
    #[test]
    fn local_base_rule_ring_mean_floor_and_scope() {
        let lat = EvalLat { gx0: -40, gy0: -40, n: 81, dx_mm: 500_000 };
        let cf = CfParamsV1::v1();
        // Flat ground at 500 m under a vent whose bilinear base is 1500 m: G = 500, H = 3500, R_f = clamp(H / 0.15) / 2.
        let g = local_ground(&lat, &vec![500.0; 81 * 81], 0.0, &[site(1500.0, 4000.0)], &cf, false);
        assert_eq!(g.len(), 1);
        assert!((g[0].ground_m - 500.0).abs() < 1e-9 && (g[0].h_m - 3500.0).abs() < 1e-9 && !g[0].floored);
        assert!((g[0].r_foot_m - 3500.0 / 0.15 / 2.0).abs() < 1e-6 && g[0].ring_posts > 0 && g[0].ring_marine_posts == 0);
        // Ground near the summit: H floored at half the bilinear H.
        let g = local_ground(&lat, &vec![3900.0; 81 * 81], 0.0, &[site(1500.0, 4000.0)], &cf, false);
        assert!(g[0].floored && (g[0].ground_m - (4000.0 - 0.5 * 2500.0)).abs() < 1e-9);
        // Marine ring posts count at sea level; a submarine vent (bilinear base below sea) is out of scope.
        let g = local_ground(&lat, &vec![-300.0; 81 * 81], 0.0, &[site(1500.0, 4000.0)], &cf, false);
        assert!(g[0].ground_m.abs() < 1e-9 && g[0].ring_marine_posts == g[0].ring_posts);
        assert!(local_ground(&lat, &vec![500.0; 81 * 81], 0.0, &[site(-100.0, 2400.0)], &cf, false).is_empty());
        // With lane VE only composites; key-off serialization has no new fields.
        assert!(!base_applies("hotspot", 1500.0, 0.0, true, true) && base_applies("hotspot", 1500.0, 0.0, true, false) && base_applies("continental_arc", 1500.0, 0.0, true, true));
        let j = serde_json::to_value(site(1500.0, 4000.0)).unwrap();
        assert!(j.get("ped_shift_m").is_none());
        let p = PedestalParamsV1::v1();
        assert!(p.validate().is_err() && p.clone().with_base(true).validate().is_ok());
        let j = serde_json::to_value(p.clone().with_solve(true)).unwrap();
        assert!(j.get("base").is_none() && j.get("base_rule").is_none() && j.get("relief").is_none());
        // Lane VOLCANIC-SPLIT edifice_local_relief_v1: H = bilinear relief, built summit G + H, R_f from H; needs ped_base.
        assert!(p.clone().with_relief(true).validate().is_err() && p.clone().with_base(true).with_relief(true).validate().is_ok());
        let g = local_ground2(&lat, &vec![500.0; 81 * 81], 0.0, &[site(1500.0, 4000.0)], &cf, false, true);
        assert!((g[0].ground_m - 500.0).abs() < 1e-9 && (g[0].h_m - 2500.0).abs() < 1e-9 && (g[0].built_summit_m.unwrap() - 3000.0).abs() < 1e-9);
        assert!((g[0].r_foot_m - (2500.0f64 / 0.15).clamp(cf.cone_wb_min_mm as f64 / 1000.0, cf.cone_wb_max_mm as f64 / 1000.0) / 2.0).abs() < 1e-6);
        assert!(serde_json::to_value(&local_ground(&lat, &vec![500.0; 81 * 81], 0.0, &[site(1500.0, 4000.0)], &cf, false)[0]).unwrap().get("built_summit_m").is_none());
    }
}
