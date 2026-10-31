//! Lane EDIFICE-FOOT (2026-10-21): opt-in window-only key `ef_cap=1` (sub-keys `ef_k_permille=`, `ef_fill_mm=`);
//! absent = byte-identical.
//!
//! PROBLEM. A composite edifice (lane CF cone, lane VE composite body) is built as an absolute surface E referenced to the
//! T0 bilinear at the vent (base_m, sloped base); the deposit over a post's own pre-edifice T1b ground is E - g, unbounded.
//! Where the T1b ground lies far below that base (C7 site 5: the ground under volcano 2's flanks is ~100-700 m against a
//! 1726 m base) the cone stands on a pedestal: measured deposit over pre-edifice ground mean 1083 m at r / R_f 0.8-1.0 and
//! up to 1530 m in the foot fade ring (records/diag_ef0_on_s5.txt), a ~1 km foot wall across every valley it meets.
//! Lane VE already bounds the broad families (shield, rift field) by `thickness_cap` = 4 H x shape (ve_shapes_v1: thin
//! distal lavas drape valleys instead of standing as a foot wall); composites were left unbounded (INFINITY).
//!
//! RULE (`edifice_foot_drape_v1`). With the key, the composite body (and the CF cone when ve is off) is at most
//! max(t_fill, k x h_flat(r)) thick over the post's own pre-edifice ground, h_flat(r) = the cone's own height above its
//! base at r on flat ground (lane CO's concave profile, record-age lowered, crater bowl; 0 at and beyond R_f). k = 4 is
//! lane VE's broad-family factor (one law for every family). t_fill (default 100 m): where the cone is thin the deposit is
//! what flows and lahars leave in the ground they reach: they run down valleys and spread (andesite flows 15-50 m thick,
//! lahars 5-20 m, VE's deposit rule stacks at most 3 t; a few stacked units ~100 m) instead of stacking into a wall.
//! Applied as E <- min(E, g + cap) before the burial / softmax blend (exactly where VE applies its broad cap), so ridges
//! the cone buries are unchanged and the capped surface is continuous (min of two continuous surfaces). VE units are built
//! on the capped body as before (their own thickness bounds). Broad families keep VE's own cap.
//!
//! Measurement dump (env var `VZ_EDIFICE_DUMP=<dir>`, read-only, never digested, no output bytes): the evaluation-lattice
//! heights before the coast / volcano stage (pre-edifice ground), before its final fill and after it, with the over,
//! edifice-style, fixed and volcano-index codes, one folder per evaluation lattice (`<dir>/<gx0>_<gy0>_<n>/`; an existing
//! folder is left untouched and the dump skipped).

use crate::t1b_coast::{edifice_at, EdificeAt, VolcanoParamsV1, VolcanoSite};
use crate::t1b_rivernet::EvalLat;
use serde::{Deserialize, Serialize};
use std::io::Write;

pub const LABEL: &str = "EXPLORATORY lane EDIFICE-FOOT 2026-10-21; window-only";
pub const RULE: &str = "edifice_foot_drape_v1: with ef_cap the composite body (lane VE) / CF cone is at most max(t_fill, k x h_flat(r)) thick over the post's own pre-edifice ground, h_flat(r) = the cone's own height above its base at r on flat ground (lane CO profile, record-age lowered, crater bowl; <= 0 at and beyond R_f), k = lane VE's broad-family factor 4 (ve_shapes_v1 thickness_cap: one law for every family), t_fill = what flows and lahars leave where the cone is thin (a few stacked 15-50 m flows / 5-20 m lahars, declared 100 m): they run down the valleys they reach and spread instead of stacking into a foot wall; E <- min(E, g + cap) before the burial / softmax blend (where VE caps broad bodies); units built on the capped body as before; with dam_lake, posts where the cap bound are pre-edifice ground under a veneer for its wholly-on-edifice test (a valley the cone dams on its draped foot is a dammed valley)";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct EdificeFootParamsV1 {
    pub on: u32,
    /// k in permille (default 4000 = lane VE's broad-family 4 H x shape).
    pub k_permille: i64,
    /// t_fill in mm (default 100 000).
    pub fill_mm: i64,
    pub rule: String,
    pub label: String,
}

impl EdificeFootParamsV1 {
    pub fn v1() -> Self { EdificeFootParamsV1 { on: 1, k_permille: 4000, fill_mm: 100_000, rule: RULE.into(), label: LABEL.into() } }
    pub fn validate(&self) -> Result<(), String> {
        if self.on != 1 { return Err("ef_cap: on must be 1 (absent = off)".into()); }
        if self.k_permille < 0 || self.fill_mm < 0 { return Err("ef_cap: k_permille / fill_mm < 0".into()); }
        Ok(())
    }
    /// The cap (m) of volcano `v` at chart (x, y).
    pub fn cap(&self, vp: &VolcanoParamsV1, v: &VolcanoSite, x: f64, y: f64) -> f64 {
        let hf = match edifice_at(vp, v, x, y) { EdificeAt::Abs(z, _) => z - v.base_m, EdificeAt::Add(dz) => dz, EdificeAt::None => 0.0 };
        (self.fill_mm as f64 / 1000.0).max(self.k_permille as f64 / 1000.0 * hf)
    }
}

/// Per-lattice stats (serialized in prep.net.coast.ef only with the key).
#[derive(Clone, Debug, Default, Serialize)]
pub struct EfCoastStats {
    /// Posts whose body surface the cap lowered (by more than 1 mm), the largest and the summed lowering (m, m3).
    pub capped_posts: u32,
    pub cut_max_m: f64,
    pub cut_volume_m3: f64,
}

/// Write the measurement dump for one evaluation lattice (see the module doc). Errors are printed, never returned.
pub fn dump(dir: &std::ffi::OsStr, lat: &EvalLat, z_pre_edifice: &[f64], z_pre_fill: &[f64], z_after: &[f64], over: &[u8], style: &[u8], fixed: &[u8], vid: &[u8]) {
    let d = std::path::Path::new(dir).join(format!("{}_{}_{}", lat.gx0, lat.gy0, lat.n));
    if d.exists() { println!("  VZ_EDIFICE_DUMP: {} exists, dump skipped (never overwritten)", d.display()); return; }
    let go = || -> std::io::Result<()> {
        std::fs::create_dir_all(&d)?;
        let wf = |name: &str, v: &[f64]| -> std::io::Result<()> {
            let mut f = std::fs::OpenOptions::new().write(true).create_new(true).open(d.join(name))?;
            let b: Vec<u8> = v.iter().flat_map(|x| (*x as f32).to_le_bytes()).collect();
            f.write_all(&b)
        };
        let wb = |name: &str, v: &[u8]| -> std::io::Result<()> { std::fs::OpenOptions::new().write(true).create_new(true).open(d.join(name))?.write_all(v) };
        wf("z_pre_edifice.f32", z_pre_edifice)?;
        wf("z_pre_fill.f32", z_pre_fill)?;
        wf("z_after.f32", z_after)?;
        wb("over.u8", over)?;
        wb("style.u8", style)?;
        wb("fixed.u8", fixed)?;
        wb("vid.u8", vid)?;
        let j = format!("{{\"gx0\": {}, \"gy0\": {}, \"n\": {}, \"dx_mm\": {}, \"style_len\": {}}}", lat.gx0, lat.gy0, lat.n, lat.dx_mm, style.len());
        std::fs::OpenOptions::new().write(true).create_new(true).open(d.join("lat.json"))?.write_all(j.as_bytes())
    };
    match go() { Ok(()) => println!("  VZ_EDIFICE_DUMP: wrote {} (measurement only)", d.display()), Err(e) => println!("  VZ_EDIFICE_DUMP: {} failed: {e}", d.display()) }
}
