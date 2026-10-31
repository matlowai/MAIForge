//! Exposed-layer provenance (lane PV, node ROCK-LOOK-PROVENANCE, 2026-10-15). LABEL: pre-signoff, fidelity tier
//! Preview, derived, never canon.
//!
//! design 2026-10-15: "i feel like the different look should be attached to the providence atleast if it matters" and
//! "geology-derived generation". This module carries the RECORDED history of the rock a post exposes (the T0 strata
//! column: lithology, depositional environment, age span, deepest burial, metamorphic grade, flags) to the T1b window
//! and the T2 tile, as a stable per-post id plus a table of records. It reads only upstream records and invents
//! nothing: unknown values are JSON null.
//!
//! Used by pf-region (T1b `rnwin prov=1`) and pf-detail (T2 `block prov=1`; a crate dependency on pf-region since lane BS,
//! 2026-10-15, formerly included by `#[path]`), so both tiers encode the same id for the same layer. Depends on std,
//! serde_json and blake3 only.
//!
//! Id encoding (u32, little endian per post):
//! - strata layer: `(t0_cell << 8) | k`, k = index of the layer in the cell's published FULL column, top = 0 (the
//!   T0 `stratafull_*` CSR order, basement last). Stable for a given T0 record. Needs t0_cell < 2^20 (the L1 grid has
//!   165 888 cells) and k < 256.
//! - `0xF000_0000 | v`: volcanic edifice of the volcano at index v of the T1b window header's `prep.net.volcanoes`
//!   list (`v` = 255: edifice, volcano not identified).
//! - `0xE000_0000 | d`: unconsolidated present-day deposit built by the T1b coast stage: d = 2 fan, 3 beach, 0 kind
//!   not identified.
//! - `0xD000_0000 | t0_cell`: no recorded column for the cell (oceanic or empty column): no history.
//! - `0xFFFF_FFFF`: none (marine or lake water).

use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::Path;

pub const PROV_NONE: u32 = u32::MAX;
pub const PROV_KIND_MASK: u32 = 0xF000_0000;
pub const PROV_EDIFICE: u32 = 0xF000_0000;
pub const PROV_DEPOSIT: u32 = 0xE000_0000;
pub const PROV_NOCOL: u32 = 0xD000_0000;
pub const DEPOSIT_UNKNOWN: u32 = 0;
pub const DEPOSIT_FAN: u32 = 2;
pub const DEPOSIT_BEACH: u32 = 3;
/// T0 flag bit (pf-contracts STRATA_FLAG_CUT_AT_DEPTH): the layer continues below the record's publish depth.
const FLAG_CUT_AT_DEPTH: u8 = 16;

pub const PROV_RULE: &str = "exposed_layer_provenance_v1 (lane PV): per post the id of the exposed unit (the layer of the T0 column that contains the published surface, by the same cut that sets the exposed lithology; the edifice or deposit over-code where it applies). The record of each id (the header / tile table) copies the T0 stratigraphy record of that layer unchanged: lith, env, age_base_my, age_top_my (model time in My since the run start; present = present_t_my), max_burial_m, grade (0 unmetamorphosed, 1 low, 2 high), flags (1 truncated_top, 2 merged, 4 organic_rich, 8 unconformity_at_base, 16 continues_below_publish_depth, 32 basement_residual), thickness_m; edifice records copy the volcano record (kind, shape, onset_my, last_eruption_my, active, stage, denudation_m); deposit records are T1b's present-day fans and beaches (undated). JSON null = not recorded. Nothing is drawn or inferred.";

pub const ID_NOTE: &str = "lane PV (opt-in prov=1): stable id of the exposed unit per post: (t0_cell << 8) | k for layer k (top = 0) of the T0 cell's full column; 0xF0000000 | v volcanic edifice of volcano index v (header prep.net.volcanoes; 255 unknown); 0xE0000000 | d present-day deposit (2 fan, 3 beach, 0 unknown); 0xD0000000 | t0_cell no recorded column; 0xFFFFFFFF none (water). The record of every id present is in the provenance table (header.provenance.rows / tile.json provenance.rows); JSON null = not recorded";

pub const UNKNOWN_NOTE: &str = "not recorded upstream, published as null: the column's last_exposure_my (pf-geodynamics StrataStack keeps it, publish() and the stratafull dumps drop it; for a column that is land now it equals the present, so it would not date an outcrop's weathering either); the protolith of metamorphic layers (lith becomes 2, env kept); sub-layer facies, grain size, cement, colour, fossils; burial, grade and flags of edifice rock; ages of deposits";

#[inline]
pub fn strata_id(cell: u32, k: u16) -> u32 { (cell << 8) | (k as u32 & 0xFF) }

fn rd(p: &Path) -> Result<Vec<u8>, String> { std::fs::read(p).map_err(|e| format!("{}: {e}", p.display())) }
fn blake(b: &[u8]) -> String { format!("blake3:{}", blake3::hash(b).to_hex()) }
fn f32s(b: &[u8]) -> Vec<f32> { b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
fn u32s(b: &[u8]) -> Vec<u32> { b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
#[inline]
fn r3(v: f64) -> Value { if v.is_finite() { json!((v * 1000.0).round() / 1000.0) } else { Value::Null } }

/// The per-layer record (CSR, top down, as `stratigraphy::publish` writes it).
pub struct StrataProv {
    /// "stratafull_dumps" (full column) or "stratigraphy_json_clipped" (StratigraphyDomainV0 to its publish depth).
    pub source: String,
    pub present_t_my: f64,
    pub publish_depth_m: f64,
    pub unconformity_gap_my: f64,
    off: Vec<u32>,
    env: Vec<u8>,
    grade: Vec<u8>,
    flags: Vec<u8>,
    age_base: Vec<f32>,
    age_top: Vec<f32>,
    burial: Vec<f32>,
    /// (role, digest) of every byte source read.
    pub digests: Vec<(String, String)>,
}

impl StrataProv {
    /// Full dumps in `final_dir` (stratafull_{off,lith,th,env,grade,flags,agebase,agetop,burial} + stratafull_meta.json)
    /// when `stratafull_env.u8` exists; else the StratigraphyDomainV0 JSON `json`, verified layer by layer against the
    /// full column the cut uses (`full_*`: lith equal; thickness equal unless the layer is cut at the publish depth).
    pub fn load(final_dir: &Path, json: Option<&Path>, full_off: &[u32], full_lith: &[u8], full_th: &[f32]) -> Result<StrataProv, String> {
        if final_dir.join("stratafull_env.u8").exists() {
            let mut digests = Vec::new();
            let mut get = |role: &str, f: &str| -> Result<Vec<u8>, String> { let b = rd(&final_dir.join(f))?; digests.push((role.to_string(), blake(&b))); Ok(b) };
            let off = u32s(&get("t0.stratigraphy.full.offset", "stratafull_off.u32")?);
            let lith = get("t0.stratigraphy.full.lith", "stratafull_lith.u8")?;
            let th = f32s(&get("t0.stratigraphy.full.thickness", "stratafull_th.f32")?);
            let env = get("t0.stratigraphy.full.env", "stratafull_env.u8")?;
            let grade = get("t0.stratigraphy.full.grade", "stratafull_grade.u8")?;
            let flags = get("t0.stratigraphy.full.flags", "stratafull_flags.u8")?;
            let age_base = f32s(&get("t0.stratigraphy.full.age_base", "stratafull_agebase.f32")?);
            let age_top = f32s(&get("t0.stratigraphy.full.age_top", "stratafull_agetop.f32")?);
            let burial = f32s(&get("t0.stratigraphy.full.burial", "stratafull_burial.f32")?);
            let meta: Value = serde_json::from_slice(&get("t0.stratigraphy.full.meta", "stratafull_meta.json")?).map_err(|e| format!("stratafull_meta.json: {e}"))?;
            let present_t_my = meta.get("t_my").and_then(Value::as_f64).ok_or("stratafull_meta.json: t_my")?;
            let unconformity_gap_my = meta.get("unconformity_gap_my").and_then(Value::as_f64).unwrap_or(f64::NAN);
            let m = lith.len();
            if [th.len(), env.len(), grade.len(), flags.len(), age_base.len(), age_top.len(), burial.len()].iter().any(|&l| l != m) || off.last().copied() != Some(m as u32) {
                return Err("stratafull provenance dump sizes differ".into());
            }
            if off != full_off || lith != full_lith || th != full_th { return Err("stratafull provenance dumps differ from the column the cut reads".into()); }
            digests.sort();
            return Ok(StrataProv { source: "stratafull_dumps".into(), present_t_my, publish_depth_m: f64::INFINITY, unconformity_gap_my, off, env, grade, flags, age_base, age_top, burial, digests });
        }
        let jp = json.ok_or_else(|| format!("{}: no stratafull_env.u8 (env/age/burial/grade/flags not dumped); give prov_strata=<StratigraphyDomainV0 JSON> of the same T0 run", final_dir.display()))?;
        let bytes = rd(jp)?;
        Self::from_json_bytes(&bytes, full_off, full_lith, full_th).map_err(|e| format!("{}: {e}", jp.display()))
    }

    /// A StratigraphyDomainV0 JSON (depth-clipped), verified against the full column (see `load`).
    pub fn from_json_bytes(bytes: &[u8], full_off: &[u32], full_lith: &[u8], full_th: &[f32]) -> Result<StrataProv, String> {
        #[derive(serde::Deserialize)]
        struct D { schema: String, t_my: f64, publish_depth_m: f64, unconformity_gap_my: f64, layer_offset: Vec<u32>, lith: Vec<u8>, env: Vec<u8>, grade: Vec<u8>, flags: Vec<u8>,
                   thickness_m: Vec<f32>, age_base_my: Vec<f32>, age_top_my: Vec<f32>, max_burial_m: Vec<f32> }
        let d: D = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
        if d.schema != "pf.stratigraphy/v0" { return Err(format!("schema {}", d.schema)); }
        let m = d.lith.len();
        if [d.env.len(), d.grade.len(), d.flags.len(), d.thickness_m.len(), d.age_base_my.len(), d.age_top_my.len(), d.max_burial_m.len()].iter().any(|&l| l != m)
            || d.layer_offset.last().copied() != Some(m as u32) || d.layer_offset.len() != full_off.len() {
            return Err("stratigraphy JSON sizes differ from the full column grid".into());
        }
        // Verify: the clipped record is the top of the full column, cell by cell.
        let mut bad = 0u64;
        for c in 0..full_off.len() - 1 {
            let (a, b) = (d.layer_offset[c] as usize, d.layer_offset[c + 1] as usize);
            let (fa, fb) = (full_off[c] as usize, full_off[c + 1] as usize);
            if b - a > fb - fa { bad += 1; continue; }
            for j in 0..b - a {
                let cut = d.flags[a + j] & FLAG_CUT_AT_DEPTH != 0;
                if d.lith[a + j] != full_lith[fa + j] || (!cut && d.thickness_m[a + j] != full_th[fa + j]) || (cut && d.thickness_m[a + j] > full_th[fa + j] + 0.01) { bad += 1; }
            }
        }
        if bad > 0 { return Err(format!("stratigraphy JSON disagrees with the full column on {bad} layers or cells; refusing")); }
        Ok(StrataProv { source: "stratigraphy_json_clipped".into(), present_t_my: d.t_my, publish_depth_m: d.publish_depth_m, unconformity_gap_my: d.unconformity_gap_my,
            off: d.layer_offset, env: d.env, grade: d.grade, flags: d.flags, age_base: d.age_base_my, age_top: d.age_top_my, burial: d.max_burial_m,
            digests: vec![("t0.stratigraphy.record_json".into(), blake(bytes))] })
    }

    /// Header block describing the source.
    pub fn header(&self) -> Value {
        json!({"rule": PROV_RULE, "id_note": ID_NOTE, "unknown": UNKNOWN_NOTE, "source": self.source, "present_t_my": self.present_t_my,
               "source_publish_depth_m": if self.publish_depth_m.is_finite() && self.publish_depth_m < 1.0e8 { json!(self.publish_depth_m) } else { Value::Null },
               "unconformity_gap_my": r3(self.unconformity_gap_my), "upstream": self.digests})
    }

    /// The record of layer k (top = 0) of `cell`'s full column (`full_lith`, `full_th`: that layer in the full column the
    /// cut reads). Layers below a clipped source's publish depth carry lith and thickness only.
    pub fn strata_row(&self, cell: u32, k: u16, full_lith: u8, full_th: f32, full_layers: usize) -> Value {
        let (a, b) = (self.off[cell as usize] as usize, self.off[cell as usize + 1] as usize);
        let i = a + k as usize;
        let mut r = json!({"id": strata_id(cell, k), "src": "strata", "cell": cell, "layer": k, "column_layers": full_layers, "lith": full_lith, "thickness_m": r3(full_th as f64),
            "env": null, "age_base_my": null, "age_top_my": null, "max_burial_m": null, "grade": null, "flags": null, "last_exposure_my": null});
        if i < b {
            r["env"] = json!(self.env[i]);
            r["age_base_my"] = r3(self.age_base[i] as f64);
            r["age_top_my"] = r3(self.age_top[i] as f64);
            r["max_burial_m"] = r3(self.burial[i] as f64);
            r["grade"] = json!(self.grade[i]);
            // the cut-at-depth bit belongs to the clipped record, not the layer: the full column the cut reads has no cut
            r["flags"] = json!(self.flags[i] & !FLAG_CUT_AT_DEPTH);
        } else {
            r["below_source_depth"] = json!(true);
        }
        r
    }
}

/// Volcano facts for an edifice record: from the T1b window list (id, kind, shape, age_my = present - last eruption,
/// stage, denudation_m) and, when the T0 volcano records are given, onset_my / last_eruption_my / active by id.
#[derive(Clone, Debug, Default)]
pub struct VolcInfo {
    pub id: String,
    pub kind: String,
    pub shape: String,
    pub age_my: f64,
    pub stage: String,
    pub denudation_m: f64,
}

/// id -> (onset_my, last_eruption_my, active) from a T0 volcanoes.json, and its digest.
pub fn volcano_records(path: &Path) -> Result<(BTreeMap<String, (Option<f64>, Option<f64>, Option<bool>)>, String), String> {
    let b = rd(path)?;
    let v: Value = serde_json::from_slice(&b).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut m = BTreeMap::new();
    for f in v.get("features").and_then(Value::as_array).ok_or("volcanoes.json: features")? {
        let Some(id) = f.get("id").and_then(Value::as_str) else { continue };
        m.insert(id.to_string(), (f.get("onset_my").and_then(Value::as_f64), f.get("last_eruption_my").and_then(Value::as_f64), f.get("active").and_then(Value::as_bool)));
    }
    Ok((m, blake(&b)))
}

pub fn edifice_row(v: u32, info: Option<&VolcInfo>, rec: Option<&(Option<f64>, Option<f64>, Option<bool>)>, present_t_my: f64) -> Value {
    let mut r = json!({"id": PROV_EDIFICE | (v & 0xFF), "src": "edifice", "volcano_index": v, "lith": 3, "env": 2,
        "volcano": null, "kind": null, "shape": null, "stage": null, "denudation_m": null, "active": null,
        "age_base_my": null, "age_top_my": null, "max_burial_m": null, "grade": null, "flags": null, "last_exposure_my": null});
    if let Some(i) = info {
        r["volcano"] = json!(i.id); r["kind"] = json!(i.kind); r["shape"] = json!(i.shape); r["stage"] = json!(i.stage); r["denudation_m"] = r3(i.denudation_m);
        r["age_top_my"] = r3(present_t_my - i.age_my);
    }
    if let Some((on, last, act)) = rec {
        if let Some(o) = on { r["age_base_my"] = r3(*o); }
        if let Some(l) = last { r["age_top_my"] = r3(*l); }
        if let Some(a) = act { r["active"] = json!(a); }
    }
    r
}

pub fn deposit_row(d: u32) -> Value {
    let kind = match d { DEPOSIT_FAN => "fan", DEPOSIT_BEACH => "beach", _ => "fan_or_beach" };
    json!({"id": PROV_DEPOSIT | d, "src": "deposit", "kind": kind, "lith": 10, "env": 5, "present_day_construct": true,
           "age_base_my": null, "age_top_my": null, "max_burial_m": null, "grade": null, "flags": null, "last_exposure_my": null})
}

pub fn nocol_row(cell: u32, fallback_lith: Option<u8>) -> Value {
    json!({"id": PROV_NOCOL | cell, "src": "no_column", "cell": cell, "lith": fallback_lith, "env": null, "age_base_my": null, "age_top_my": null,
           "max_burial_m": null, "grade": null, "flags": null, "last_exposure_my": null})
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn id_ranges_do_not_collide() {
        let max_strata = strata_id((1 << 20) - 1, 255);
        assert!(max_strata < PROV_NOCOL && PROV_NOCOL < PROV_DEPOSIT && PROV_DEPOSIT < PROV_EDIFICE);
        assert_eq!(strata_id(42569, 3) >> 8, 42569);
        assert_eq!(strata_id(42569, 3) & 0xFF, 3);
        assert_eq!(PROV_EDIFICE | 255, 0xF000_00FF);
        assert_ne!(PROV_EDIFICE | 255, PROV_NONE);
    }
    #[test]
    fn clipped_json_is_verified_against_the_full_column() {
        // full column of 2 cells: cell 0 = [shale 10, basalt 5000], cell 1 = [oceanic 7000]
        let full_off = vec![0u32, 2, 3];
        let full_lith = vec![5u8, 3, 11];
        let full_th = vec![10.0f32, 5000.0, 7000.0];
        let doc = |th1: f32, lith0: u8| serde_json::to_vec(&json!({"schema": "pf.stratigraphy/v0", "t_my": 400.0, "publish_depth_m": 3000.0, "unconformity_gap_my": 5.0,
            "layer_offset": [0, 2, 3], "lith": [lith0, 3, 11], "env": [6, 2, 10], "grade": [0, 0, 0], "flags": [1, 16, 48], "thickness_m": [10.0, th1, 3000.0],
            "age_base_my": [380.0, 370.0, 300.0], "age_top_my": [390.0, 375.0, 300.0], "max_burial_m": [10.0, 5010.0, 7000.0]})).unwrap();
        let s = StrataProv::from_json_bytes(&doc(2990.0, 5), &full_off, &full_lith, &full_th).unwrap();
        let r = s.strata_row(0, 1, 3, 5000.0, 2);
        assert_eq!(r["env"], json!(2));
        assert_eq!(r["flags"], json!(0), "the clipped cut bit is not a layer property");
        assert_eq!(r["max_burial_m"], json!(5010.0));
        assert_eq!(s.strata_row(0, 0, 5, 10.0, 2)["flags"], json!(1));
        // a lith disagreement is refused; a cut layer thicker than the full layer is refused
        assert!(StrataProv::from_json_bytes(&doc(2990.0, 4), &full_off, &full_lith, &full_th).is_err());
        assert!(StrataProv::from_json_bytes(&doc(6000.0, 5), &full_off, &full_lith, &full_th).is_err());
    }
}
