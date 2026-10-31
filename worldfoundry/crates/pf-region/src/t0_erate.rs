//! Lane HANG-1 (2026-10-22): the read-only loader of the T0 record's per-cell recent erosion rate E_rec, ported from lane
//! MS-S0S1's preserved `t0_erate.rs` (`unpublished-design-record`; loader only: none of S1's hillslope
//! switches, kappa reading, crest cap or L3 are ported). Read by exactly one law: the hanging-gully depth of the ceil sub-key
//! `ceil_hang=erec` (CEIL_HANG_RULE in t1b_rivernet.rs). No hillslope law reads it. Absent key = every byte unchanged.

use crate::t0::T0;
use crate::R;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Arc;

pub const EREC_RECORD_RULE: &str = "record_erosion_rate_v1 loader (lane MS-S0S1, ported read-only by lane HANG-1): E_rec per T0 cell = the T0 fluvial-trace frame `frame_step` surface_budget.final_counters.bedrock_removed_m / the frame's dt (the per-cell recent bedrock rate AGE-RECORD-TEST validated); where that is 0, the lifetime mean cumulative_erosion_m / min(crust_age_my, the frame's end time); at least floor_m_per_my [H: project assumption]; bound to T0 by the terrain digest and to the run by status.json's state / recipe digests; hat-blended at a node with T0's hat weights, a cell whose solved table value is 0 contributing 0 (as e_at)";

/// The record table (not serialized; bound by `source_digest`).
#[derive(Clone, Default)]
pub struct ErecTable {
    /// E_rec (m/My) per T0 element.
    pub e: Arc<Vec<f64>>,
    /// 0 recent (frame), 1 lifetime mean, 2 floor.
    pub src: Arc<Vec<u8>>,
}
impl std::fmt::Debug for ErecTable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "ErecTable({} cells)", self.e.len()) }
}

/// The E_rec record as a params block (digests serialized; the table itself is not).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ErecRecordV1 {
    pub rule: String,
    /// blake3 of fluvial-trace.json bytes (the upstream record; the terrain fields are covered by `terrain_digest`).
    pub source_digest: String,
    pub terrain_digest: String,
    pub frame_step: u32,
    pub floor_m_per_my: i64,
    #[serde(skip)]
    pub rec: ErecTable,
}
impl PartialEq for ErecRecordV1 {
    fn eq(&self, o: &Self) -> bool {
        self.rule == o.rule && self.source_digest == o.source_digest && self.terrain_digest == o.terrain_digest && self.frame_step == o.frame_step && self.floor_m_per_my == o.floor_m_per_my
    }
}
impl Eq for ErecRecordV1 {}

#[derive(Deserialize)]
struct Tr { recipe_digest: String, state_digest: String, trace: TrIn }
#[derive(Deserialize)]
struct TrIn { frames: Vec<TrFrame> }
#[derive(Deserialize)]
struct TrFrame { clock: TrClock, surface_budget: TrSb }
#[derive(Deserialize)]
struct TrClock { start_my: String, end_my: String, step: u32 }
#[derive(Deserialize)]
struct TrSb { final_counters: TrFc }
#[derive(Deserialize)]
struct TrFc { bedrock_removed_m: Vec<String> }

fn hex_f64(s: &str) -> R<f64> {
    if s.len() != 16 { return Err(format!("erec: bad float hex {s}")); }
    u64::from_str_radix(s, 16).map(f64::from_bits).map_err(|e| format!("erec: bad float hex {s}: {e}"))
}

impl ErecRecordV1 {
    /// Load and bind to T0: the run's status.json must cite T0's terrain digest and the trace's state / recipe digests.
    pub fn load(dir: &Path, t0: &T0, frame_step: u32, floor_m_per_my: i64) -> R<ErecRecordV1> {
        let tdir = dir.join("terrain");
        let terrain_digest = pf_contracts::artifact::verify_artifact::<Value>(&tdir).map_err(|e| format!("erec: verify terrain: {e}"))?;
        if terrain_digest != t0.terrain_digest { return Err(format!("erec: terrain {terrain_digest} is not T0's {}; refusing", t0.terrain_digest)); }
        let st: Value = serde_json::from_slice(&std::fs::read(dir.join("status.json")).map_err(|e| format!("erec: status.json: {e}"))?).map_err(|e| format!("erec: status.json: {e}"))?;
        if st.pointer("/artifact_digests/terrain").and_then(Value::as_str) != Some(t0.terrain_digest.as_str()) { return Err("erec: status.json cites another terrain; refusing".into()); }
        let bytes = std::fs::read(dir.join("fluvial-trace.json")).map_err(|e| format!("erec: fluvial-trace.json: {e}"))?;
        let source_digest = format!("blake3:{}", blake3::hash(&bytes).to_hex());
        let tr: Tr = serde_json::from_slice(&bytes).map_err(|e| format!("erec: fluvial-trace.json parse: {e}"))?;
        drop(bytes);
        if Some(tr.state_digest.as_str()) != st.get("state_digest").and_then(Value::as_str) || Some(tr.recipe_digest.as_str()) != st.get("recipe_digest").and_then(Value::as_str) {
            return Err("erec: fluvial-trace state / recipe digest differs from status.json; refusing".into());
        }
        let f = tr.trace.frames.iter().find(|f| f.clock.step == frame_step).ok_or_else(|| format!("erec: the trace has no frame at step {frame_step}"))?;
        let (t_start, t_end) = (hex_f64(&f.clock.start_my)?, hex_f64(&f.clock.end_my)?);
        let dt = t_end - t_start;
        if !(dt > 0.0) { return Err("erec: frame dt <= 0".into()); }
        let n = t0.nlat * t0.nlon;
        let br = &f.surface_budget.final_counters.bedrock_removed_m;
        if br.len() != n { return Err(format!("erec: bedrock_removed_m has {} cells, T0 {n}", br.len())); }
        let rd = |name: &str| -> R<Vec<f32>> {
            let v = pf_contracts::artifact::read_field_f32(&tdir, name).map_err(|e| format!("erec: {name}: {e}"))?;
            if v.len() != n { return Err(format!("erec: {name} size")); }
            Ok(v)
        };
        let cum = rd("cumulative_erosion_m")?;
        let age = rd("crust_age_my")?;
        let floor = floor_m_per_my as f64;
        let (mut e, mut src) = (vec![0.0f64; n], vec![0u8; n]);
        for c in 0..n {
            let recent = hex_f64(&br[c])? / dt;
            let life = (age[c] as f64).min(t_end);
            let lifetime = if life > 0.0 { (cum[c] as f64).max(0.0) / life } else { 0.0 };
            let (v, s) = if recent > 0.0 { (recent, 0u8) } else { (lifetime, 1u8) };
            if v >= floor { e[c] = v; src[c] = s; } else { e[c] = floor; src[c] = 2; }
        }
        Ok(ErecRecordV1 { rule: EREC_RECORD_RULE.into(), source_digest, terrain_digest, frame_step, floor_m_per_my, rec: ErecTable { e: Arc::new(e), src: Arc::new(src) } })
    }
    pub fn validate(&self) -> R<()> {
        if self.rule != EREC_RECORD_RULE { return Err("erec: unknown record rule".into()); }
        if self.rec.e.is_empty() { return Err("erec: the record table is not loaded (ErecRecordV1::load)".into()); }
        if self.floor_m_per_my <= 0 { return Err("erec: floor must be > 0".into()); }
        Ok(())
    }
    /// E_rec at (la, lo) over T0's hat weights (a cell whose table value is 0 contributes 0, as in `e_at`).
    pub fn at(&self, t0: &T0, table: &BTreeMap<u32, f64>, la: f64, lo: f64) -> f64 {
        let (c, w) = crate::t1b_net::hat4(t0, la, lo);
        let mut v = 0.0f64;
        for k in 0..4 { if w[k] > 0.0 && table.get(&c[k]).copied().unwrap_or(0.0) > 0.0 { v += w[k] * self.rec.e[c[k] as usize]; } }
        v
    }
}
