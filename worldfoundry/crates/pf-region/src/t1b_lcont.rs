//! Lane LITHO-CONTINUITY (2026-10-21): rock layers that do not jump at T0 cell edges. LABEL: pre-signoff, fidelity tier
//! Preview, derived, never canon. Opt-in (`lcont=1`, layer AND windows, strata= rock source only); absent = every byte
//! unchanged.
//!
//! Design principle (operator-side review 2026-10-21, binding for the lane): "Geological formations should have
//! planet-wide identities and histories, while tiles merely sample their local geometry." No column cross-fade: every
//! unit is one 3D body whose thickness (and so its top and base contacts) is sampled at T0 cell centres and interpolated
//! between them, with the same identity on both sides of a cell edge; a unit ends only where its thickness goes to zero
//! (pinch-out toward a cell whose column has no matching unit), where the surface cuts it (erosion), or at a declared
//! conflict pinch-out (below). No fault or intrusion boundary is produced (the record locates none inside a 70 km cell).
//!
//! What the record carries (measured on C7, see the lane REPORT): the T0 strata dumps the cut reads (stratafull_off / lith
//! / th, strata_grad_e / n) carry per cell an ordered list of layers with a rock code and a thickness only; the clipped
//! stratigraphy JSON adds env, age span, burial, grade and flags down to the 3000 m publish depth. There is NO formation,
//! unit or event id shared between cells: each column is the stack of the cell's dominant parcel. Identity is therefore
//! DEFINED here by a deterministic matching rule (`IDENTITY_RULE`) over the dumps the cut already reads; its limits are
//! stated in the rule text.

use crate::t0::{CellWater, T0};
use crate::t1b_rivernet::StrataCols;
use crate::t1b_volc::{SplitCol, VolcanicSplitV1, FAC_LAVA, FAC_NONE};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const LABEL: &str = "EXPLORATORY lane LITHO-CONTINUITY 2026-10-21 (operator-side principle 2026-10-21: planet-wide formation identities, tiles sample local geometry; no column cross-fade); layer AND windows";

pub const IDENTITY_RULE: &str = "lcont_identity_v2: kept bodies = the strata cut's layers (>= 1 m, top down) of every T0 cell, each run of consecutive layers of one rock code taken as one body (recording layer: the run's top layer; the record's internal contacts of such a run, e.g. an unconformity between two basalt layers, are not kept); a cell is a land column iff T0 water = Land and its kept column is not empty; for every pair of 4-neighbour land columns a monotone (order-preserving, one-to-one) alignment pairs bodies of the SAME rock code maximizing the sum of min(thickness) (dynamic programming; traceback prefers a pair, then skipping the lower-id cell's body); the pairs are candidate unions, taken planet-wide strongest first (min thickness descending; ties edge key, index) and accepted iff the merged body still holds at most one layer per cell and every 2x2 quad of cell centres holding a cell of each side keeps an acyclic stacking order (ORDER_RULE); the planet-wide formation id is the smallest global kept index of the merged body. Limits: rock code and depth order only (the dumps the cut reads carry no age or environment; the clipped JSON's ages are not used), so two different events of one rock code in one order position read as one body; a rejected pair leaves two bodies that pinch out against each other (cause: identity conflict)";

pub const ORDER_RULE: &str = "lcont_order_v1: an edge's stacking order is its matched layers in order with, between two matches, the unmatched layers of both columns merged by u = the layer's mid-depth below its own column top (m) (ties: lower cell id, then layer index); a quad's order is the topological order of its columns' and edges' orders (ties: smallest u over the body's members, then node index), so on every quad edge it restricts to that edge's order";

pub const GEOMETRY_RULE: &str = "lcont_geometry_v1: at a post (lat, lon) the quad of the 4 surrounding T0 cell centres (the bilinear stencil) gives weights w = smoothstep (3t^2 - 2t^3) of the bilinear fractions, kept only on land columns and renormalized (no land column, or a polar half row: the post's own cell alone); datum D = sum w_c (E_c + g_c . (east, north) of the post from cell c's centre on the sphere), E_c the T0 cell elevation, g_c the T0 uplift gradient clamped to dip_clamp (each column's tilted plane anchored at ITS cell centre: planet-keyed, not the site); each body's thickness at the post = sum w_c thickness_c (0 in a column without it); bodies stack down from D in the quad order; a body thinner than min_unit at the post is absorbed by the next body below (its band takes the lower body's identity); the exposed code at height z is the body containing depth D - z (above D: the first body; below the last base: the last). C1 across quad edges in weights; contacts continuous everywhere except at the border of a fallback stencil (counted)";

pub const VOLCANIC_RULE: &str = "lcont_volcanic_v1 (with vs=1): a lithology-3 body is split as ONE body keyed by its planet-wide formation id: the volcanic_facies_split_v2 package stack of the representative layer (cell r, layer k of the formation root; the body is in scope iff ANY of its member layers' T0 cells is, one decision per body) of thickness T_r, mapped proportionally into the body at the post (depth e = (d - top) T_r / t, t the body's thickness at the post; riser minimum scaled the same way); the post's first body keeps the stack's upward extension and the last body the downward one; the riser ladder walks the contiguous run of split bodies; lens (lateral pinch-out) positions are keyed in the body's own planet frame (east, north m on the local tangent plane at the representative layer's T0 cell centre), not chart metres, so every chart, window and lattice sees the same lenses; package draws are keyed by the formation's representative layer and depth within the body (material coordinates in depth)";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct LithoContinuityV1 {
    pub on: u32,
    /// Bodies thinner than this at a post are absorbed by the body below (mm) [1000, the cut's kept-layer minimum].
    pub min_unit_mm: i64,
    pub identity_rule: String,
    pub order_rule: String,
    pub geometry_rule: String,
    pub volcanic_rule: String,
    pub label: String,
    /// Lane LITHO-AGE-GATE (opt-in `lcont_record=<json>`): age and environment gate on the identity from the T0
    /// stratigraphy record (RECORD_RULE); absent = unchanged (not serialized).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub record: Option<LcontRecordV1>,
}

impl LithoContinuityV1 {
    pub fn v1() -> Self {
        LithoContinuityV1 { on: 1, min_unit_mm: 1_000, identity_rule: IDENTITY_RULE.into(), order_rule: ORDER_RULE.into(), geometry_rule: GEOMETRY_RULE.into(),
            volcanic_rule: VOLCANIC_RULE.into(), label: LABEL.into(), record: None }
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.on != 1 { return Err("lcont: on must be 1 (absent = off)".into()); }
        if self.min_unit_mm < 0 { return Err("lcont: min_unit_mm >= 0".into()); }
        if let Some(r) = self.record.as_ref() {
            if r.age_tol_kyr.is_some_and(|k| k < 0) { return Err("lcont_record: age_tol_kyr >= 0".into()); }
            if r.m2o_rule.as_ref().is_some_and(|v| v != M2O_RULE) { return Err("lcont_m2o: unknown rule".into()); }
            if r.env_tol > 3 { return Err("lcont_record: env_tol in 0..=3".into()); }
        }
        Ok(())
    }
}

// ------------------------------------------------------------------ lane LITHO-AGE-GATE (2026-10-22)

pub const RECORD_RULE: &str = "lcont_record_v1 (lane LITHO-AGE-GATE, opt-in lcont_record=<StratigraphyDomainV0 JSON of the same T0 run>, its blake3 digest in params): the clipped T0 stratigraphy record, verified layer by layer against the full column the cut reads (lith equal; thickness equal unless the layer is cut at the publish depth) and mapped by the record's layer_offset CSR (not layer_count, which counts the full column), gives each layer above the publish depth its age span [age_base_my, age_top_my], environment and flags; deeper layers carry no record and are never gated. (1) Bodies: a run of consecutive kept layers of one rock code (lcont_identity_v2) is split where the record puts an unconformity (flag 8 unconformity_at_base on any full layer from the upper kept layer down to, not including, the lower one) or where adding the layer would make the body's environments incompatible (2); a body's record is the union of its recorded layers' age spans and environment tracts. (2) Gate: two bodies may be one body only if their age spans overlap within the age tolerance tol (base_a <= top_b + tol and base_b <= top_a + tol; tol = age_tol_kyr, default the record's unconformity_gap_my, the gap below which the record itself calls a contact conformable) and their environments are compatible: the same class among basement 0, pluton 1, subaerial volcanic 2, oceanic crust 10 and depositional; depositional environments sit on a facies tract (fluvial 3, lacustrine 4, evaporitic 8, swamp 9 = 0; deltaic / coastal 5 = 1; shelf 6 = 2; deep marine 7 = 3) and are compatible iff their tract indices differ by at most env_tol (default 1: laterally adjacent belts, Walther's law); a body with no recorded layer passes. The gate applies to the pairwise alignment (a gated pair is not a candidate, so the alignment may pick another same-code pair) and to every union (a merged body keeps the intersection of its members' widened spans [base, top + tol] and the range of their tracts; a union is refused when the intersection is empty or the range exceeds env_tol, which by Helly's theorem in one dimension is exactly 'every member pair passes the gate'). Everything else is lcont_identity_v2 unchanged.";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct LcontRecordV1 {
    /// blake3 of the record JSON bytes ("blake3:<hex>").
    pub digest: String,
    /// Age tolerance (kyr); None = the record's unconformity_gap_my.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub age_tol_kyr: Option<i64>,
    /// Largest facies-tract distance between joined depositional environments [1].
    pub env_tol: u32,
    pub rule: String,
    /// Opt-in conservative subdivision; absent keeps the prior recipe bytes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub m2o_rule: Option<String>,
    /// Where the record is read (not digested; the digest identifies the bytes and is re-checked at use).
    #[serde(skip)]
    pub path: Option<String>,
}

impl LcontRecordV1 {
    /// From the record file (reads it once for the digest).
    pub fn from_path(p: &str) -> Result<LcontRecordV1, String> {
        let b = std::fs::read(p).map_err(|e| format!("lcont_record {p}: {e}"))?;
        Ok(LcontRecordV1 { digest: format!("blake3:{}", blake3::hash(&b).to_hex()), age_tol_kyr: None, env_tol: 1, rule: RECORD_RULE.into(), m2o_rule: None, path: Some(p.to_string()) })
    }
}

/// The clipped record per FULL layer (global index `cols.off[c] + k`), verified against the full column.
pub struct Record {
    pub known: Vec<bool>,
    pub age_base: Vec<f32>,
    pub age_top: Vec<f32>,
    pub env: Vec<u8>,
    pub flags: Vec<u8>,
    pub unconformity_gap_my: f64,
}

const FLAG_UNCONFORMITY_BASE: u8 = 8;
const FLAG_CUT_AT_DEPTH: u8 = 16;

impl Record {
    pub fn load(r: &LcontRecordV1, cols: &StrataCols) -> Result<Record, String> {
        let p = r.path.as_ref().ok_or("lcont_record: the record path is not known (params read from a recipe?); give lcont_record=<json>")?;
        let bytes = std::fs::read(p).map_err(|e| format!("lcont_record {p}: {e}"))?;
        let dg = format!("blake3:{}", blake3::hash(&bytes).to_hex());
        if dg != r.digest { return Err(format!("lcont_record {p}: digest {dg} differs from the params' {}", r.digest)); }
        Self::from_json_bytes(&bytes, &cols.off, &cols.lith, &cols.th).map_err(|e| format!("lcont_record {p}: {e}"))
    }
    pub fn from_json_bytes(bytes: &[u8], off: &[u32], lith: &[u8], th: &[f32]) -> Result<Record, String> {
        #[derive(Deserialize)]
        struct D { schema: String, unconformity_gap_my: f64, layer_offset: Vec<u32>, lith: Vec<u8>, env: Vec<u8>, flags: Vec<u8>, thickness_m: Vec<f32>, age_base_my: Vec<f32>, age_top_my: Vec<f32> }
        let d: D = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
        if d.schema != "pf.stratigraphy/v0" { return Err(format!("schema {}", d.schema)); }
        let m = d.lith.len();
        if [d.env.len(), d.flags.len(), d.thickness_m.len(), d.age_base_my.len(), d.age_top_my.len()].iter().any(|&l| l != m)
            || d.layer_offset.last().copied() != Some(m as u32) || d.layer_offset.len() != off.len() {
            return Err("record sizes differ from the full column grid".into());
        }
        let n = th.len();
        let mut rec = Record { known: vec![false; n], age_base: vec![0.0; n], age_top: vec![0.0; n], env: vec![0; n], flags: vec![0; n], unconformity_gap_my: d.unconformity_gap_my };
        let mut bad = 0u64;
        for c in 0..off.len() - 1 {
            let (a, b) = (d.layer_offset[c] as usize, d.layer_offset[c + 1] as usize);
            let (fa, fb) = (off[c] as usize, off[c + 1] as usize);
            if b - a > fb - fa { bad += 1; continue; }
            for j in 0..b - a {
                let cut = d.flags[a + j] & FLAG_CUT_AT_DEPTH != 0;
                if d.lith[a + j] != lith[fa + j] || (!cut && d.thickness_m[a + j] != th[fa + j]) || (cut && d.thickness_m[a + j] > th[fa + j] + 0.01) { bad += 1; continue; }
                let g = fa + j;
                rec.known[g] = true; rec.age_base[g] = d.age_base_my[a + j]; rec.age_top[g] = d.age_top_my[a + j]; rec.env[g] = d.env[a + j];
                rec.flags[g] = d.flags[a + j] & !FLAG_CUT_AT_DEPTH;
            }
        }
        if bad > 0 { return Err(format!("record disagrees with the full column on {bad} layers or cells; refusing")); }
        if !(rec.unconformity_gap_my >= 0.0) { return Err("record unconformity_gap_my".into()); }
        Ok(rec)
    }
}

const ANY: u8 = u8::MAX;
const CONFLICT: u8 = u8::MAX - 1;

/// (environment class, facies tract) of a recorded environment (RECORD_RULE (2)); tract ANY = none.
fn env_class(env: u8) -> (u8, u8) {
    match env {
        3 | 4 | 8 | 9 => (3, 0),
        5 => (3, 1),
        6 => (3, 2),
        7 => (3, 3),
        e => (e, ANY),
    }
}

/// A body's (or a merged body's) record: widened age interval [lo, hi], environment class and tract range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sig { pub lo: f64, pub hi: f64, pub class: u8, pub tmin: u8, pub tmax: u8 }

impl Sig {
    /// No record: passes every gate.
    pub const NONE: Sig = Sig { lo: f64::NEG_INFINITY, hi: f64::INFINITY, class: ANY, tmin: ANY, tmax: 0 };
    /// One recorded layer (age span widened by tol at its top).
    pub fn layer(age_base: f64, age_top: f64, env: u8, tol: f64) -> Sig {
        let (class, t) = env_class(env);
        Sig { lo: age_base, hi: age_top + tol, class, tmin: t, tmax: if t == ANY { 0 } else { t } }
    }
    fn class_meet(a: u8, b: u8) -> u8 { if a == ANY { b } else if b == ANY || a == b { a } else { CONFLICT } }
    fn tracts(&self, o: &Sig) -> (u8, u8) { (self.tmin.min(o.tmin), self.tmax.max(o.tmax)) }
    /// One cell's conformable layers taken together (union of spans).
    pub fn hull(&self, o: &Sig) -> Sig {
        if self.lo == f64::NEG_INFINITY && self.hi == f64::INFINITY && self.class == ANY && self.tmin == ANY { return *o; }
        if o.lo == f64::NEG_INFINITY && o.hi == f64::INFINITY && o.class == ANY && o.tmin == ANY { return *self; }
        let (tmin, tmax) = self.tracts(o);
        Sig { lo: self.lo.min(o.lo), hi: self.hi.max(o.hi), class: Self::class_meet(self.class, o.class), tmin, tmax }
    }
    /// Bodies of different cells joined (intersection of widened spans).
    pub fn meet(&self, o: &Sig) -> Sig {
        let (tmin, tmax) = self.tracts(o);
        Sig { lo: self.lo.max(o.lo), hi: self.hi.min(o.hi), class: Self::class_meet(self.class, o.class), tmin, tmax }
    }
    pub fn age_ok(&self) -> bool { self.lo <= self.hi }
    pub fn env_ok(&self, env_tol: u32) -> bool { self.class != CONFLICT && (self.tmin == ANY || (self.tmax - self.tmin) as u32 <= env_tol) }
}

/// Gate statistics (stdout only; never digested).
#[derive(Clone, Debug, Default, Serialize)]
pub struct GateStats {
    pub age_tol_my: f64,
    pub env_tol: u32,
    /// Kept bodies before / after splitting same-code runs at recorded boundaries.
    pub bodies_ungated_runs: u64,
    pub bodies: u64,
    pub runs_split_unconformity: u64,
    pub runs_split_env: u64,
    pub bodies_with_record: u64,
    pub bodies_without_record: u64,
    /// Pairs of the ungated alignment (same bodies) that the gate refuses.
    pub alignment_pairs_gated: u64,
    pub rejected_union_age: u64,
    pub rejected_union_env: u64,
}
/// Identity statistics (stdout only; never digested).
#[derive(Clone, Debug, Default, Serialize)]
pub struct IdentityStats {
    pub land_columns: u64,
    pub land_edges: u64,
    /// Pairwise alignment matches (candidate unions).
    pub candidates: u64,
    pub accepted: u64,
    /// Candidates whose two layers were already one body through other edges.
    pub already_joined: u64,
    pub rejected_two_layers_one_cell: u64,
    pub rejected_order: u64,
    /// Surviving same-body pairs across land edges.
    pub edge_pairs_final: u64,
    pub formations: u64,
    pub formations_multi_cell: u64,
    pub largest_formation_cells: u64,
    pub seconds: f64,
}

/// The planet-wide identity (built once per process from the strata dumps and the T0 water mask).
pub struct Identity {
    pub nlat: usize,
    pub nlon: usize,
    /// Per cell: CSR offsets into the kept bodies (runs of one rock code).
    pub koff: Vec<u32>,
    pub kth: Vec<f64>,
    pub kcode: Vec<u8>,
    /// Index in the cell's FULL column (top = 0) of the run's top layer.
    pub kidx: Vec<u16>,
    /// Mid-depth below the column top (m).
    pub ku: Vec<f64>,
    /// Full column length per cell.
    pub full: Vec<u32>,
    pub land: Vec<bool>,
    /// (lower cell id, higher cell id) -> same-body pairs (kept index in the lower-id cell, in the other) (dump / stats).
    pub edges: BTreeMap<(u32, u32), Vec<(u16, u16)>>,
    /// Planet-wide formation id (smallest global kept index of the body) per kept layer.
    pub root: Vec<u32>,
    pub stats: IdentityStats,
    /// Lane LITHO-AGE-GATE: the gate's statistics (None: ungated).
    pub gate: Option<GateStats>,
}

const QE: [(usize, usize); 4] = [(0, 1), (0, 2), (1, 3), (2, 3)];

fn align(a: &[(f64, u8)], b: &[(f64, u8)]) -> Vec<(u16, u16)> { align_g(a, b, &|_, _| true) }

/// `align` with a gate: pair (i, j) allowed iff same code and `gate(i, j)` (lane LITHO-AGE-GATE).
fn align_g(a: &[(f64, u8)], b: &[(f64, u8)], gate: &dyn Fn(usize, usize) -> bool) -> Vec<(u16, u16)> {
    let (n, m) = (a.len(), b.len());
    let w = m + 1;
    let mut s = vec![0.0f64; (n + 1) * w];
    let ok = |i: usize, j: usize| a[i].1 == b[j].1 && gate(i, j);
    for i in (0..n).rev() {
        for j in (0..m).rev() {
            let mut best = s[(i + 1) * w + j].max(s[i * w + j + 1]);
            if ok(i, j) { best = best.max(a[i].0.min(b[j].0) + s[(i + 1) * w + j + 1]); }
            s[i * w + j] = best;
        }
    }
    let (mut i, mut j, mut out) = (0usize, 0usize, Vec::new());
    while i < n && j < m {
        if ok(i, j) && s[i * w + j] == a[i].0.min(b[j].0) + s[(i + 1) * w + j + 1] { out.push((i as u16, j as u16)); i += 1; j += 1; }
        else if s[i * w + j] == s[(i + 1) * w + j] { i += 1; }
        else { j += 1; }
    }
    out
}

/// One cell's kept bodies ((thickness, code, full index of the top layer), top down) from its kept layers (depth bases,
/// codes, full indices; `fo` = the cell's offset in the full CSR): runs of one rock code; with the gate (lane
/// LITHO-AGE-GATE, RECORD_RULE (1)) runs are split at recorded unconformities and incompatible environments, and each body
/// gets its record. Counts: [same-code runs, splits at an unconformity, splits at an environment].
fn cell_bodies(base: &[f64], code: &[u8], idx: &[u16], fo: usize, gate: Option<(&Record, f64, u32)>) -> (Vec<(f64, u8, u16)>, Vec<Sig>, [u64; 3]) {
    let mut prev = 0.0;
    let mut out: Vec<(f64, u8, u16)> = Vec::new();
    let mut osig: Vec<Sig> = Vec::new();
    let mut cnt = [0u64; 3];
    let mut last_ix = 0u16;
    for ((&b, &k), &ix) in base.iter().zip(code).zip(idx) {
        let th = b - prev; prev = b;
        let Some((rec, tol, et)) = gate else {
            // a run of consecutive kept layers of one rock code is one body (the cut cannot tell them apart)
            match out.last_mut() { Some(l) if l.1 == k => l.0 += th, _ => out.push((th, k, ix)) }
            continue;
        };
        // lane LITHO-AGE-GATE: ... unless the record puts an unconformity or an incompatible environment between
        let g = fo + ix as usize;
        let ls = if rec.known[g] { Sig::layer(rec.age_base[g] as f64, rec.age_top[g] as f64, rec.env[g], tol) } else { Sig::NONE };
        let same = matches!(out.last(), Some(l) if l.1 == k);
        let mut join = same;
        if same {
            let unc = (last_ix as usize..ix as usize).any(|q| rec.known[fo + q] && rec.flags[fo + q] & FLAG_UNCONFORMITY_BASE != 0);
            let h = osig.last().unwrap().hull(&ls);
            if unc { join = false; cnt[1] += 1; } else if !h.env_ok(et) { join = false; cnt[2] += 1; }
            if join { *osig.last_mut().unwrap() = h; }
        } else { cnt[0] += 1; }
        if join { out.last_mut().unwrap().0 += th; } else { out.push((th, k, ix)); osig.push(ls); }
        last_ix = ix;
    }
    (out, osig, cnt)
}

// ------------------------------------------------------------------ lane CX-LITHO-M2O (2026-10-22)

pub const M2O_RULE: &str = "lcont_m2o_v2: opt-in with lcont_record. Before identity alignment, a recorded body may be subdivided to correspond to two or more ordered same-code bodies of a neighbouring land column. Every matched member must have a record; their widened age intervals must have a common intersection with the thick body's interval and their environment ranges must pass the unchanged RECORD_RULE (Helly gate). Target age endpoints must be nonincreasing down-section. No recorded unconformity may occur anywhere between the first and last target, including intervening or discarded thin layers. The complete target span including intervening facies must be no thicker than the source body, so this models thinning and lateral pinch-out, not amplification. Source thickness is partitioned by recorded time under a constant mean accumulation-rate approximation: internal time cuts are the midpoints between each upper target age-base and the next lower target age-top, bounded by the source interval; cuts must be strictly decreasing and every resulting slice at least 1 m. Slice thickness = source thickness times slice duration / source duration; the last slice takes the residual. Zero-duration or reversed chronologies are refused, never replaced by thickness fractions. Each slice inherits the source environment and a narrowed age interval between its cuts (plus the original tolerance), and must still pass the corresponding target record gate. Total source thickness and its full age envelope are conserved. Intervening facies have zero thickness in the source cell. All slices retain their real source layer. The time cuts are explicitly inferred internal chronology, not new T0 event observations; original recorded ages are never overwritten. Choose proposals by matched thickness minus thickest target (descending), then cell/index keys; a source is subdivided only once and is never also a target; an unsplit target may anchor subdivisions in multiple neighbouring source cells. Re-run the existing one-to-one alignment on refined bodies and retain every original age, environment, one-member-per-cell and acyclic-quad union gate. Subdivision does not guarantee the proposed joins survive global conflicts. Geological basis: conformable lateral facies succession (Walther; SEPM Strata Terminology Walther), sedimentary interbeds pinching out within basalt piles (USGS WRI 84-4071). This is a reduced-order correlation hypothesis, not proof of common eruption or a deposition simulation. Walther applies to conformable sedimentary facies, not a license to join igneous events. Unknown deep bodies never initiate or join a subdivision. Units: thickness m, age My.";

#[derive(Clone, Debug, Default, Serialize)]
struct M2oStats {
    proposals: u64,
    selected: u64,
    added_slices: u64,
    rejected_overlap: u64,
    max_column_thickness_error_m: f64,
}

/// Conservative one-to-many correspondence represented as several distinct, ordered identities in the thick cell.
/// Keeping them distinct preserves intervening shale rather than making a cyclic basalt-shale-basalt identity.
fn refine_m2o(nlat: usize, nlon: usize, cols: &mut Vec<Vec<(f64, u8, u16)>>, sig: &mut Vec<Sig>,
              boundaries: &[Vec<bool>], land: &[bool], et: u32, tol: f64) -> M2oStats {
    #[derive(Clone)]
    struct Proposal { c: usize, a: usize, d: usize, targets: Vec<usize>, score: f64, parts: Vec<(f64, Sig)> }
    let mut stats = M2oStats::default();
    let mut offsets = vec![0usize];
    for col in cols.iter() { offsets.push(offsets.last().unwrap() + col.len()); }
    let mut proposals = Vec::new();
    for c in 0..cols.len() {
        if !land[c] { continue; }
        let (i, j) = (c / nlon, c % nlon);
        let mut neighbours = vec![i * nlon + (j + 1) % nlon, i * nlon + (j + nlon - 1) % nlon];
        if i > 0 { neighbours.push(c - nlon); }
        if i + 1 < nlat { neighbours.push(c + nlon); }
        neighbours.sort(); neighbours.dedup();
        for d in neighbours {
            if d == c || !land[d] { continue; }
            for (a, &(thick, code, _)) in cols[c].iter().enumerate() {
                let source = sig[offsets[c] + a];
                if source == Sig::NONE { continue; }
                for start in 0..cols[d].len() {
                    if cols[d][start].1 != code { continue; }
                    let mut meet = source;
                    let mut targets = Vec::new();
                    let (mut span, mut total, mut largest) = (0.0f64, 0.0f64, 0.0f64);
                    let mut previous: Option<Sig> = None;
                    for b in start..cols[d].len() {
                        if b > start && boundaries[d][b - 1] { break; }
                        span += cols[d][b].0;
                        if span > thick { break; }
                        if cols[d][b].1 != code { continue; }
                        let next = sig[offsets[d] + b];
                        if next == Sig::NONE { break; }
                        if previous.is_some_and(|p| next.lo > p.lo || next.hi > p.hi) { break; }
                        meet = meet.meet(&next);
                        if !meet.age_ok() || !meet.env_ok(et) { break; }
                        previous = Some(next);
                        targets.push(b); total += cols[d][b].0; largest = largest.max(cols[d][b].0);
                        if targets.len() > 1 {
                            let duration = source.hi - tol - source.lo;
                            if !(duration > 0.0) { continue; }
                            let mut cuts = vec![source.hi - tol];
                            for pair in targets.windows(2) {
                                let upper = sig[offsets[d] + pair[0]];
                                let lower = sig[offsets[d] + pair[1]];
                                cuts.push(0.5 * (upper.lo + lower.hi - tol));
                            }
                            cuts.push(source.lo);
                            if cuts.windows(2).any(|w| w[0] <= w[1]) { continue; }
                            let mut parts = Vec::new();
                            let mut remaining = thick;
                            for (i, &target) in targets.iter().enumerate() {
                                let t = if i + 1 == targets.len() { remaining } else { thick * ((cuts[i] - cuts[i + 1]) / duration) };
                                let part = Sig { lo: cuts[i + 1], hi: cuts[i] + tol, ..source };
                                let m = part.meet(&sig[offsets[d] + target]);
                                if t < 1.0 || !t.is_finite() || !m.age_ok() || !m.env_ok(et) { break; }
                                parts.push((t, part)); remaining -= t;
                            }
                            if parts.len() == targets.len() {
                                proposals.push(Proposal { c, a, d, targets: targets.clone(), score: total - largest, parts });
                            }
                        }
                    }
                }
            }
        }
    }
    stats.proposals = proposals.len() as u64;
    proposals.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.c.cmp(&b.c)).then(a.a.cmp(&b.a))
        .then(a.d.cmp(&b.d)).then(a.targets.cmp(&b.targets)));
    let mut sources = BTreeSet::new();
    let mut targets = BTreeSet::new();
    let mut splits: BTreeMap<(usize, usize), Vec<(f64, Sig)>> = BTreeMap::new();
    for p in proposals {
        if sources.contains(&(p.c, p.a)) || targets.contains(&(p.c, p.a)) || p.targets.iter().any(|&b| sources.contains(&(p.d, b))) {
            stats.rejected_overlap += 1; continue;
        }
        let parts = p.parts;
        sources.insert((p.c, p.a));
        for &b in &p.targets { targets.insert((p.d, b)); }
        stats.selected += 1; stats.added_slices += parts.len() as u64 - 1;
        splits.insert((p.c, p.a), parts);
    }
    let mut refined = Vec::with_capacity(cols.len());
    let mut refined_sig = Vec::new();
    for (c, col) in cols.iter().enumerate() {
        let mut out = Vec::new();
        for (a, &(th, code, ix)) in col.iter().enumerate() {
            if let Some(parts) = splits.get(&(c, a)) {
                for &(t, age) in parts { out.push((t, code, ix)); refined_sig.push(age); }
            } else { out.push((th, code, ix)); refined_sig.push(sig[offsets[c] + a]); }
        }
        let err = (out.iter().map(|x| x.0).sum::<f64>() - col.iter().map(|x| x.0).sum::<f64>()).abs();
        stats.max_column_thickness_error_m = stats.max_column_thickness_error_m.max(err);
        assert!(err < 1e-8, "lcont_m2o: thickness conservation failed");
        assert!(out.len() <= u16::MAX as usize, "lcont_m2o: too many bodies");
        refined.push(out);
    }
    *cols = refined; *sig = refined_sig;
    stats
}

/// A quad's bodies in stacking order: per body the member kept index (local to the cell) in each slot.
pub type QuadOrder = Vec<[Option<u16>; 4]>;

#[inline]
fn find_ro(par: &[u32], mut x: u32) -> u32 { while par[x as usize] != x { x = par[x as usize]; } x }

impl Identity {
    /// Build from per-cell kept columns ((thickness, code, full index) top down), full lengths and the land mask.
    pub fn from_columns(nlat: usize, nlon: usize, cols: &[Vec<(f64, u8, u16)>], full: Vec<u32>, land_water: &[bool]) -> Identity {
        Self::from_columns_g(nlat, nlon, cols, full, land_water, None)
    }

    /// `from_columns` with the lane LITHO-AGE-GATE gate: per kept body (CSR order of `cols`) its record and the env tolerance;
    /// None = ungated (identical to `from_columns`).
    pub fn from_columns_g(nlat: usize, nlon: usize, cols: &[Vec<(f64, u8, u16)>], full: Vec<u32>, land_water: &[bool], gate: Option<(&[Sig], u32)>) -> Identity {
        let t = std::time::Instant::now();
        let n = nlat * nlon;
        let mut koff = Vec::with_capacity(n + 1);
        let (mut kth, mut kcode, mut kidx, mut ku) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        koff.push(0u32);
        for col in cols.iter().take(n) {
            let mut d = 0.0;
            for &(th, code, ix) in col { kth.push(th); kcode.push(code); kidx.push(ix); ku.push(d + 0.5 * th); d += th; }
            koff.push(kth.len() as u32);
        }
        let land: Vec<bool> = (0..n).map(|c| land_water[c] && koff[c + 1] > koff[c]).collect();
        let nk = kth.len();
        let mut id = Identity { nlat, nlon, koff, kth, kcode, kidx, ku, full, land, edges: BTreeMap::new(), root: Vec::new(), stats: IdentityStats::default(), gate: None };
        let mut gs = GateStats::default();
        if let Some((sig, et)) = gate {
            assert_eq!(sig.len(), nk, "lcont_record: one record per kept body");
            gs.env_tol = et;
            gs.bodies = nk as u64;
            gs.bodies_with_record = sig.iter().filter(|s| **s != Sig::NONE).count() as u64;
            gs.bodies_without_record = gs.bodies - gs.bodies_with_record;
        }
        let pass = |x: usize, y: usize| -> bool { match gate { None => true, Some((sig, et)) => { let m = sig[x].meet(&sig[y]); m.age_ok() && m.env_ok(et) } } };
        id.stats.land_columns = id.land.iter().filter(|&&l| l).count() as u64;
        // Candidate unions: pairwise alignments of 4-neighbour land columns, strongest first.
        let mut cand: Vec<(f64, (u32, u32), u16, u16)> = Vec::new();
        let mut keys: Vec<(u32, u32)> = Vec::new();
        for i in 0..nlat {
            for j in 0..nlon {
                let c = (i * nlon + j) as u32;
                if !id.land[c as usize] { continue; }
                let mut nb = vec![(i * nlon + (j + 1) % nlon) as u32];
                if i + 1 < nlat { nb.push(((i + 1) * nlon + j) as u32); }
                for d in nb {
                    if d == c || !id.land[d as usize] { continue; }
                    let key = (c.min(d), c.max(d));
                    keys.push(key);
                    let (ca, cb) = (id.col(key.0), id.col(key.1));
                    let (ga, gb) = (id.koff[key.0 as usize] as usize, id.koff[key.1 as usize] as usize);
                    if gate.is_some() { gs.alignment_pairs_gated += align(&ca, &cb).iter().filter(|&&(a, b)| !pass(ga + a as usize, gb + b as usize)).count() as u64; }
                    for (a, b) in align_g(&ca, &cb, &|a, b| pass(ga + a, gb + b)) {
                        let s = id.kth[id.g(key.0, a)].min(id.kth[id.g(key.1, b)]);
                        cand.push((s, key, a, b));
                    }
                }
            }
        }
        keys.sort(); keys.dedup();
        id.stats.land_edges = keys.len() as u64;
        cand.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)).then(x.3.cmp(&y.3)));
        cand.dedup_by(|x, y| x.1 == y.1 && x.2 == y.2 && x.3 == y.3);
        id.stats.candidates = cand.len() as u64;
        // Greedy union: keep each body one layer per cell and every quad's stacking order acyclic.
        let mut par: Vec<u32> = (0..nk as u32).collect();
        let mut minrep: Vec<u32> = (0..nk as u32).collect();
        let mut cells_of: Vec<Vec<u32>> = (0..n as u32).flat_map(|c| (id.koff[c as usize]..id.koff[c as usize + 1]).map(move |_| vec![c])).collect();
        // lane LITHO-AGE-GATE: per root the merged body's record (meet of its members)
        let mut rsig: Vec<Sig> = match gate { Some((sig, _)) => sig.to_vec(), None => Vec::new() };
        for &(_, key, a, b) in &cand {
            let (x, y) = (id.g(key.0, a) as u32, id.g(key.1, b) as u32);
            let (rx, ry) = (find_ro(&par, x), find_ro(&par, y));
            if rx == ry { id.stats.already_joined += 1; continue; }
            let (sm, lg) = if cells_of[rx as usize].len() < cells_of[ry as usize].len() || (cells_of[rx as usize].len() == cells_of[ry as usize].len() && rx > ry) { (rx, ry) } else { (ry, rx) };
            let mut merged = None;
            if let Some((_, et)) = gate {
                let m = rsig[sm as usize].meet(&rsig[lg as usize]);
                if !m.age_ok() { gs.rejected_union_age += 1; continue; }
                if !m.env_ok(et) { gs.rejected_union_env += 1; continue; }
                merged = Some(m);
            }
            let lgset = &cells_of[lg as usize];
            if cells_of[sm as usize].iter().any(|c| lgset.binary_search(c).is_ok()) { id.stats.rejected_two_layers_one_cell += 1; continue; }
            // quads holding a cell of each body
            let mut quads: BTreeSet<u32> = BTreeSet::new();
            for &c in &cells_of[sm as usize] {
                let (i, j) = ((c as usize / nlon) as i64, (c as usize % nlon) as i64);
                for (di, dj) in [(-1i64, -1i64), (-1, 0), (0, -1), (0, 0)] {
                    let (qi, qj) = (i + di, (j + dj).rem_euclid(nlon as i64));
                    if qi < 0 || qi + 1 >= nlat as i64 { continue; }
                    let q = (qi as usize * nlon + qj as usize) as u32;
                    if id.quad_cells(q).iter().any(|cc| lgset.binary_search(cc).is_ok()) { quads.insert(q); }
                }
            }
            let fid = |k: usize| -> u32 { let r = find_ro(&par, k as u32); if r == sm { lg } else { r } };
            if quads.iter().any(|&q| id.quad_order(q, &fid).is_none()) { id.stats.rejected_order += 1; continue; }
            par[sm as usize] = lg;
            if let Some(m) = merged { rsig[lg as usize] = m; }
            minrep[lg as usize] = minrep[lg as usize].min(minrep[sm as usize]);
            let s = std::mem::take(&mut cells_of[sm as usize]);
            let l = std::mem::take(&mut cells_of[lg as usize]);
            let mut m = Vec::with_capacity(s.len() + l.len());
            let (mut p, mut r) = (0, 0);
            while p < s.len() || r < l.len() { if r >= l.len() || (p < s.len() && s[p] < l[r]) { m.push(s[p]); p += 1; } else { m.push(l[r]); r += 1; } }
            cells_of[lg as usize] = m;
            id.stats.accepted += 1;
        }
        id.root = (0..nk as u32).map(|x| minrep[find_ro(&par, x) as usize]).collect();
        for &key in &keys {
            let (a0, a1) = (id.koff[key.0 as usize], id.koff[key.0 as usize + 1]);
            let (b0, b1) = (id.koff[key.1 as usize], id.koff[key.1 as usize + 1]);
            let mut v = Vec::new();
            for x in a0..a1 { for y in b0..b1 { if id.root[x as usize] == id.root[y as usize] { v.push(((x - a0) as u16, (y - b0) as u16)); } } }
            id.stats.edge_pairs_final += v.len() as u64;
            id.edges.insert(key, v);
        }
        let mut size: BTreeMap<u32, u64> = BTreeMap::new();
        for &r in &id.root { *size.entry(r).or_insert(0) += 1; }
        id.stats.formations = size.len() as u64;
        id.stats.formations_multi_cell = size.values().filter(|&&v| v > 1).count() as u64;
        id.stats.largest_formation_cells = size.values().copied().max().unwrap_or(0);
        id.stats.seconds = t.elapsed().as_secs_f64();
        if gate.is_some() { id.gate = Some(gs); }
        id
    }

    /// From the strata dumps and T0.
    pub fn build(cols: &StrataCols, t0: &T0) -> Identity { Self::build_g(cols, t0, None) }

    /// `build` with the lane LITHO-AGE-GATE gate (record, age tolerance My, env tolerance); None = `build`.
    pub fn build_g(cols: &StrataCols, t0: &T0, gate: Option<(&Record, f64, u32)>) -> Identity {
        Self::build_g_m2o(cols, t0, gate, false)
    }

    fn build_g_m2o(cols: &StrataCols, t0: &T0, gate: Option<(&Record, f64, u32)>, m2o: bool) -> Identity {
        let n = t0.nlat * t0.nlon;
        let mut kc = Vec::with_capacity(n);
        let mut full = Vec::with_capacity(n);
        let mut sigs: Vec<Sig> = Vec::new();
        let (mut runs, mut split_unc, mut split_env) = (0u64, 0u64, 0u64);
        for c in 0..n as u32 {
            full.push(cols.off[c as usize + 1] - cols.off[c as usize]);
            let fo = cols.off[c as usize] as usize;
            kc.push(cols.kept(c, |(base, code, idx)| {
                let (out, osig, cnt) = cell_bodies(base, code, idx, fo, gate);
                sigs.extend(osig);
                runs += cnt[0]; split_unc += cnt[1]; split_env += cnt[2];
                out
            }));
        }
        let lw: Vec<bool> = (0..n as u32).map(|c| t0.water(c) == CellWater::Land).collect();
        let Some((_, tol, et)) = gate else {
            let id = Identity::from_columns(t0.nlat, t0.nlon, &kc, full, &lw);
            println!("  lane LITHO-CONTINUITY identity: {}", serde_json::to_string(&id.stats).unwrap_or_default());
            if let Some(p) = std::env::var_os("VZ_LCONT_DUMP") { id.dump(std::path::Path::new(&p)); }
            return id;
        };
        if m2o {
            let rec = gate.unwrap().0;
            // A boundary after body b includes all full layers up to the next kept body, including discarded layers.
            let boundaries: Vec<Vec<bool>> = kc.iter().enumerate().map(|(c, col)| col.iter().enumerate().map(|(b, &(_, _, ix))| {
                let end = col.get(b + 1).map(|x| x.2 as usize).unwrap_or(full[c] as usize);
                (ix as usize..end).any(|k| { let g = cols.off[c] as usize + k; rec.known[g] && rec.flags[g] & FLAG_UNCONFORMITY_BASE != 0 })
            }).collect()).collect();
            let ms = refine_m2o(t0.nlat, t0.nlon, &mut kc, &mut sigs, &boundaries, &lw, et, tol);
            println!("  lane CX-LITHO-M2O subdivision: {}", serde_json::to_string(&ms).unwrap());
        }
        let mut id = Identity::from_columns_g(t0.nlat, t0.nlon, &kc, full, &lw, Some((&sigs, et)));
        if let Some(g) = id.gate.as_mut() { g.age_tol_my = tol; g.bodies_ungated_runs = runs; g.runs_split_unconformity = split_unc; g.runs_split_env = split_env; }
        println!("  lane LITHO-AGE-GATE identity: {} gate: {}", serde_json::to_string(&id.stats).unwrap_or_default(), serde_json::to_string(&id.gate).unwrap_or_default());
        if let Some(p) = std::env::var_os("VZ_LCONT_GATE_DUMP") {
            let p = std::path::PathBuf::from(p);
            id.dump(&p.with_extension("edges"));
            id.dump_bodies(&p.with_extension("bodies"));
            if m2o {
                // Explicit thicknesses are essential: several slices can share the same honest source-layer index.
                let file = std::fs::OpenOptions::new().write(true).create_new(true).open(p.with_extension("virtual.json"))
                    .expect("lcont_m2o: create measurement dump");
                serde_json::to_writer(file, &kc).expect("lcont_m2o: write measurement dump");
                let ages: Vec<_> = sigs.iter().map(|s| if *s == Sig::NONE { None } else { Some((s.lo, s.hi, s.class, s.tmin, s.tmax)) }).collect();
                let file = std::fs::OpenOptions::new().write(true).create_new(true).open(p.with_extension("signatures.json"))
                    .expect("lcont_m2o: create inferred-age dump");
                serde_json::to_writer(file, &ages).expect("lcont_m2o: write inferred-age dump");
            }
        }
        id
    }

    /// Measurement only (lane LITHO-AGE-GATE): per cell with kept bodies, the full index of each body's top layer (never
    /// digested; an existing file is kept).
    fn dump_bodies(&self, p: &std::path::Path) {
        use std::io::Write;
        let mut s = String::new();
        for c in 0..self.nlat * self.nlon {
            let (a, b) = (self.koff[c] as usize, self.koff[c + 1] as usize);
            if a == b { continue; }
            s.push_str(&c.to_string());
            for k in a..b { s.push_str(&format!(" {}", self.kidx[k])); }
            s.push('\n');
        }
        let r = std::fs::OpenOptions::new().write(true).create_new(true).open(p).and_then(|mut f| f.write_all(s.as_bytes()));
        match r { Ok(()) => println!("  VZ_LCONT_GATE_DUMP: wrote {} (measurement only)", p.display()), Err(e) => println!("  VZ_LCONT_GATE_DUMP: {} not written: {e}", p.display()) }
    }

    /// Lane LITHO-AGE-GATE: the gated identity for `r`, built once per process per (record, tolerances, strata, land mask).
    pub fn gated(r: &LcontRecordV1, cols: &StrataCols, t0: &T0) -> std::sync::Arc<Identity> {
        static GATED: std::sync::Mutex<Vec<(String, std::sync::Arc<Identity>)>> = std::sync::Mutex::new(Vec::new());
        let land: Vec<u8> = (0..(t0.nlat * t0.nlon) as u32).map(|c| (t0.water(c) == CellWater::Land) as u8).collect();
        let key = format!("{}|{:?}|{}|{:?}|{}|{:?}", r.digest, r.age_tol_kyr, r.env_tol, cols.digests, blake3::hash(&land).to_hex(), r.m2o_rule);
        let mut g = GATED.lock().unwrap_or_else(|e| e.into_inner());
        if let Some((_, id)) = g.iter().find(|(k, _)| *k == key) { return id.clone(); }
        let rec = Record::load(r, cols).unwrap_or_else(|e| panic!("{e}"));
        let tol = r.age_tol_kyr.map(|k| k as f64 / 1000.0).unwrap_or(rec.unconformity_gap_my);
        let id = std::sync::Arc::new(Identity::build_g_m2o(cols, t0, Some((&rec, tol, r.env_tol)), r.m2o_rule.is_some()));
        g.push((key, id.clone()));
        id
    }

    /// Measurement only: env `VZ_LCONT_DUMP=<file>` writes the same-body pairs per land edge (never digested; an existing file is kept).
    fn dump(&self, p: &std::path::Path) {
        use std::io::Write;
        let mut s = String::new();
        for (&(a, b), ms) in &self.edges {
            s.push_str(&format!("{a} {b}"));
            for &(i, j) in ms { s.push_str(&format!(" {i}:{j}")); }
            s.push('\n');
        }
        let r = std::fs::OpenOptions::new().write(true).create_new(true).open(p).and_then(|mut f| f.write_all(s.as_bytes()));
        match r { Ok(()) => println!("  VZ_LCONT_DUMP: wrote {} (measurement only)", p.display()), Err(e) => println!("  VZ_LCONT_DUMP: {} not written: {e}", p.display()) }
    }

    fn col(&self, c: u32) -> Vec<(f64, u8)> {
        let (a, b) = (self.koff[c as usize] as usize, self.koff[c as usize + 1] as usize);
        (a..b).map(|k| (self.kth[k], self.kcode[k])).collect()
    }
    #[inline] fn nk(&self, c: u32) -> usize { (self.koff[c as usize + 1] - self.koff[c as usize]) as usize }
    #[inline] fn g(&self, c: u32, i: u16) -> usize { self.koff[c as usize] as usize + i as usize }

    pub fn quad_cells(&self, q: u32) -> [u32; 4] {
        let (i0, j0) = (q as usize / self.nlon, q as usize % self.nlon);
        let j1 = (j0 + 1) % self.nlon;
        [(i0 * self.nlon + j0) as u32, (i0 * self.nlon + j1) as u32, ((i0 + 1) * self.nlon + j0) as u32, ((i0 + 1) * self.nlon + j1) as u32]
    }

    /// The quad's body order under the final identity (None: inconsistent; never after `from_columns`).
    pub fn check(&self, q: u32) -> Option<QuadOrder> { self.quad_order(q, &|k| self.root[k]) }

    /// The quad's body order with body id `fid(global kept index)`: None if a body holds two layers of one cell, an edge's
    /// common bodies are in different orders, or the constraints (columns, edge orders) are cyclic (ORDER_RULE).
    pub fn quad_order(&self, q: u32, fid: &dyn Fn(usize) -> u32) -> Option<QuadOrder> {
        let cells = self.quad_cells(q);
        let land: [bool; 4] = std::array::from_fn(|s| self.land[cells[s] as usize] && (s == 0 || !cells[..s].contains(&cells[s])));
        let mut members: BTreeMap<u32, [Option<u16>; 4]> = BTreeMap::new();
        let mut chains: [Vec<u32>; 4] = Default::default();
        for s in 0..4 {
            if !land[s] { continue; }
            for a in 0..self.nk(cells[s]) {
                let f = fid(self.g(cells[s], a as u16));
                let e = members.entry(f).or_insert([None; 4]);
                if e[s].is_some() { return None; }
                e[s] = Some(a as u16);
                chains[s].push(f);
            }
        }
        let mut succ: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
        let mut add = |u: u32, v: u32| { if u != v { succ.entry(u).or_default().insert(v); } };
        for ch in &chains { for w in ch.windows(2) { add(w[0], w[1]); } }
        for &(p, r) in QE.iter() {
            if !(land[p] && land[r]) { continue; }
            let (cp, cr) = (&chains[p], &chains[r]);
            let common_p: Vec<usize> = (0..cp.len()).filter(|&a| members[&cp[a]][r].is_some()).collect();
            let common_r: Vec<usize> = (0..cr.len()).filter(|&b| members[&cr[b]][p].is_some()).collect();
            if common_p.iter().map(|&a| cp[a]).ne(common_r.iter().map(|&b| cr[b])) { return None; }
            let mut seq: Vec<u32> = Vec::new();
            let (mut pi, mut pj) = (0usize, 0usize);
            let anchors: Vec<(usize, usize)> = common_p.iter().zip(&common_r).map(|(&a, &b)| (a, b)).chain(std::iter::once((cp.len(), cr.len()))).collect();
            for (mi, mj) in anchors {
                let mut gap: Vec<(f64, u32, usize, u32)> = Vec::new();
                for a in pi..mi { gap.push((self.ku[self.g(cells[p], a as u16)], cells[p], a, cp[a])); }
                for b in pj..mj { gap.push((self.ku[self.g(cells[r], b as u16)], cells[r], b, cr[b])); }
                gap.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
                seq.extend(gap.iter().map(|g| g.3));
                if mi < cp.len() { seq.push(cp[mi]); }
                pi = mi + 1; pj = mj + 1;
            }
            for w in seq.windows(2) { add(w[0], w[1]); }
        }
        let mut indeg: BTreeMap<u32, usize> = members.keys().map(|&k| (k, 0)).collect();
        for vs in succ.values() { for v in vs { *indeg.get_mut(v)? += 1; } }
        let key = |f: u32| -> f64 {
            let m = &members[&f];
            let mut u = f64::INFINITY;
            for s in 0..4 { if let Some(a) = m[s] { u = u.min(self.ku[self.g(cells[s], a)]); } }
            u
        };
        let mut avail: Vec<u32> = indeg.iter().filter(|(_, &d)| d == 0).map(|(&k, _)| k).collect();
        let mut order: QuadOrder = Vec::with_capacity(members.len());
        while !avail.is_empty() {
            let (bi, _) = avail.iter().enumerate().min_by(|a, b| key(*a.1).total_cmp(&key(*b.1)).then(a.1.cmp(b.1)))?;
            let u = avail.swap_remove(bi);
            order.push(members[&u]);
            if let Some(vs) = succ.get(&u) { for &v in vs { let d = indeg.get_mut(&v)?; *d -= 1; if *d == 0 { avail.push(v); } } }
        }
        if order.len() != members.len() { return None; }
        Some(order)
    }
}

/// One body of a stencil (top down).
#[derive(Clone, Debug)]
pub struct Body {
    pub th: [f64; 4],
    pub code: u8,
    /// Member per slot: (cell, full index) or None.
    pub src: [Option<(u32, u16)>; 4],
    /// Representative layer (formation root): (cell, full index, kept thickness).
    pub rep: (u32, u16, f64),
    /// Planet-wide formation id (global kept index of the root).
    pub root: u32,
    /// Volcanic split column index (u32::MAX: not split).
    pub vc: u32,
}

pub struct Stencil { pub cells: [u32; 4], pub bodies: Vec<Body> }

#[derive(Clone, Debug, Default, Serialize)]
pub struct ContStats { pub posts: u64, pub quad_posts: u64, pub fallback_no_land_posts: u64, pub fallback_polar_posts: u64, pub fallback_inconsistent_stencils: u64, pub stencils: u64 }

/// Per-post access for one lattice.
pub struct ContPosts {
    pub min_unit: f64,
    pub stencils: Vec<Stencil>,
    pub st: Vec<u32>,
    pub w: Vec<[f64; 4]>,
    pub datum: Vec<f64>,
    pub grad: Vec<(f64, f64)>,
    /// Per post (lat, lon) degrees: lateral detail of a split body is keyed in the body's own planet frame (bxy).
    pub ll: Vec<(f64, f64)>,
    pub radius_m: f64,
    pub nlat: usize,
    pub nlon: usize,
    pub own: Vec<u32>,
    pub vrule: Option<VolcanicSplitV1>,
    /// Per split body: (column clipped at the body base with the upward extension, column extended downward).
    pub vcols: Vec<(SplitCol, SplitCol)>,
    pub stats: ContStats,
}

/// One unit of a post's column: body index, top and base depth below the datum, first / last flags.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Unit { pub k: usize, pub top: f64, pub base: f64, pub first: bool, pub last: bool }

#[inline] fn smooth(t: f64) -> f64 { t * t * (3.0 - 2.0 * t) }

impl ContPosts {
    /// Weights and stencil key at (lat, lon): Ok((quad id, weights)) or Err(own cell) (single-cell fallback; the flag says why:
    /// 1 polar, 2 no land column).
    fn stencil_at(id: &Identity, la: f64, lo: f64) -> Result<(u32, [f64; 4]), u8> {
        let dlat = 180.0 / id.nlat as f64;
        let dlon = 360.0 / id.nlon as f64;
        let fi = (la + 90.0) / dlat - 0.5;
        let fj = (lo + 180.0) / dlon - 0.5;
        let i0 = fi.floor();
        let j0 = fj.floor();
        if i0 < 0.0 || i0 as usize + 1 >= id.nlat { return Err(1); }
        let (si, sj) = (smooth(fi - i0), smooth(fj - j0));
        let q = (i0 as usize * id.nlon + (j0 as i64).rem_euclid(id.nlon as i64) as usize) as u32;
        let cells = id.quad_cells(q);
        let mut w = [(1.0 - si) * (1.0 - sj), (1.0 - si) * sj, si * (1.0 - sj), si * sj];
        let mut sum = 0.0;
        for s in 0..4 { if !id.land[cells[s] as usize] { w[s] = 0.0; } sum += w[s]; }
        if !(sum > 0.0) { return Err(2); }
        for v in w.iter_mut() { *v /= sum; }
        Ok((q, w))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn build(rule: &LithoContinuityV1, id: &Identity, cols: &StrataCols, t0: &T0, rock_family: &[u8], dip_clamp: f64, vrule: Option<&VolcanicSplitV1>,
                 pos: &[(f64, f64, f64, f64)]) -> ContPosts {
        // Lane LITHO-AGE-GATE: with lcont_record the gated identity replaces the caller's (built once per process).
        let gated = rule.record.as_ref().map(|r| Identity::gated(r, cols, t0));
        let id: &Identity = gated.as_deref().unwrap_or(id);
        let mut cp = ContPosts { min_unit: rule.min_unit_mm as f64 / 1000.0, stencils: Vec::new(), st: Vec::with_capacity(pos.len()), w: Vec::with_capacity(pos.len()),
            datum: Vec::with_capacity(pos.len()), grad: Vec::with_capacity(pos.len()), ll: Vec::with_capacity(pos.len()), radius_m: t0.radius_m, nlat: id.nlat, nlon: id.nlon, own: Vec::with_capacity(pos.len()),
            vrule: vrule.cloned(), vcols: Vec::new(), stats: ContStats::default() };
        let mut skey: BTreeMap<u64, u32> = BTreeMap::new();
        // Split scope per formation: in scope iff any member layer's T0 cell is in scope (one decision for the whole body).
        let fscope: BTreeSet<u32> = match vrule {
            Some(v) => (0..id.nlat * id.nlon).filter(|&c| v.cell_in_scope(rock_family[c])).flat_map(|c| (id.koff[c]..id.koff[c + 1]).map(|k| id.root[k as usize])).collect(),
            None => BTreeSet::new(),
        };
        let mut vkey: BTreeMap<u32, u32> = BTreeMap::new();
        let r = t0.radius_m;
        let dlat = 180.0 / id.nlat as f64;
        let dlon = 360.0 / id.nlon as f64;
        let plane = |c: u32, la: f64, lo: f64| -> f64 {
            let (i, j) = (c as usize / id.nlon, c as usize % id.nlon);
            let (lac, loc) = (-90.0 + (i as f64 + 0.5) * dlat, -180.0 + (j as f64 + 0.5) * dlon);
            let dlo = (lo - loc + 540.0).rem_euclid(360.0) - 180.0;
            let east = r * lac.to_radians().cos() * dlo.to_radians();
            let north = r * (la - lac).to_radians();
            let (ge, gn) = (cols.grad_e[c as usize] as f64, cols.grad_n[c as usize] as f64);
            let m = ge.hypot(gn);
            let s = if m > dip_clamp { dip_clamp / m } else { 1.0 };
            t0.elevation_m[c as usize] as f64 + ge * s * east + gn * s * north
        };
        let datum_at = |la: f64, lo: f64| -> f64 {
            match Self::stencil_at(id, la, lo) {
                Ok((q, w)) => { let cells = id.quad_cells(q); (0..4).filter(|&s| w[s] > 0.0).map(|s| w[s] * plane(cells[s], la, lo)).sum() }
                _ => plane(t0.cell_of_deg(la, lo), la, lo),
            }
        };
        let rep_of = |gk: usize| -> ((u32, u16, f64), u32) {
            let rt = id.root[gk] as usize;
            let c = id.koff.partition_point(|&o| o as usize <= rt) - 1;
            ((c as u32, id.kidx[rt], id.kth[rt]), rt as u32)
        };
        let single = |c: u32| -> Stencil {
            let (a, b) = (id.koff[c as usize] as usize, id.koff[c as usize + 1] as usize);
            Stencil { cells: [c; 4], bodies: (a..b).map(|g| { let (rep, root) = rep_of(g); Body { th: [id.kth[g], 0.0, 0.0, 0.0], code: id.kcode[g], src: [Some((c, id.kidx[g])), None, None, None], rep, root, vc: u32::MAX } }).collect() }
        };
        for &(la, lo, x, y) in pos {
            let own = t0.cell_of_deg(la, lo);
            cp.stats.posts += 1;
            let (key, w) = match Self::stencil_at(id, la, lo) {
                Ok((q, w)) => { cp.stats.quad_posts += 1; (q as u64, w) }
                Err(f) => { if f == 1 { cp.stats.fallback_polar_posts += 1 } else { cp.stats.fallback_no_land_posts += 1 }; ((1u64 << 32) | own as u64, [1.0, 0.0, 0.0, 0.0]) }
            };
            let si = match skey.get(&key) {
                Some(&s) => s,
                None => {
                    let mut stc = if key >> 32 == 0 {
                        let q = key as u32;
                        let cells = id.quad_cells(q);
                        match id.check(q) {
                            Some(order) => Stencil { cells, bodies: order.iter().map(|m| {
                                let mut b = Body { th: [0.0; 4], code: 0, src: [None; 4], rep: (0, 0, 0.0), root: 0, vc: u32::MAX };
                                for s in 0..4 { if let Some(a) = m[s] { let g = id.g(cells[s], a); b.th[s] = id.kth[g]; b.code = id.kcode[g]; b.src[s] = Some((cells[s], id.kidx[g])); let (r, rt) = rep_of(g); b.rep = r; b.root = rt; } }
                                b
                            }).collect() },
                            None => { cp.stats.fallback_inconsistent_stencils += 1; single(own) }
                        }
                    } else { single(key as u32) };
                    if let Some(v) = vrule {
                        for b in stc.bodies.iter_mut() {
                            if b.code != crate::t1b_litho::LITHO_BASALT || !fscope.contains(&b.root) { continue; }
                            let ix = *vkey.entry(b.root).or_insert_with(|| {
                                let t = b.rep.2;
                                let mid = crate::t1b_volc::split_column(v, b.rep.0, true, &[t, t + 1.0], &[crate::t1b_litho::LITHO_BASALT, 255], &[b.rep.1, b.rep.1]).unwrap_or_default();
                                let bot = crate::t1b_volc::split_column(v, b.rep.0, true, &[t], &[crate::t1b_litho::LITHO_BASALT], &[b.rep.1]).unwrap_or_default();
                                cp.vcols.push((mid, bot));
                                (cp.vcols.len() - 1) as u32
                            });
                            b.vc = ix;
                        }
                    }
                    cp.stencils.push(stc);
                    let s = (cp.stencils.len() - 1) as u32;
                    skey.insert(key, s);
                    s
                }
            };
            let cells = cp.stencils[si as usize].cells;
            let d = if key >> 32 == 0 { (0..4).filter(|&s| w[s] > 0.0).map(|s| w[s] * plane(cells[s], la, lo)).sum() } else { plane(own, la, lo) };
            let dd = 1e-4;
            let gn = (datum_at(la + dd, lo) - datum_at(la - dd, lo)) / (2.0 * r * dd.to_radians());
            let ge = (datum_at(la, lo + dd) - datum_at(la, lo - dd)) / (2.0 * r * la.to_radians().cos().max(1e-6) * dd.to_radians());
            cp.st.push(si); cp.w.push(w); cp.datum.push(d); cp.grad.push((ge, gn)); cp.ll.push((la, lo)); let _ = (x, y); cp.own.push(own);
        }
        cp.stats.stencils = cp.stencils.len() as u64;
        cp
    }

    #[inline]
    fn bth(&self, i: usize, b: &Body) -> f64 { let w = &self.w[i]; w[0] * b.th[0] + w[1] * b.th[1] + w[2] * b.th[2] + w[3] * b.th[3] }

    /// The unit containing depth d (below the datum) at post i (see GEOMETRY_RULE); None for an empty stencil.
    pub fn unit_at(&self, i: usize, d: f64) -> Option<Unit> {
        let st = &self.stencils[self.st[i] as usize];
        let (mut cum, mut prev, mut first) = (0.0, 0.0, true);
        let mut found: Option<Unit> = None;
        let mut last: Option<Unit> = None;
        for (k, b) in st.bodies.iter().enumerate() {
            let t = self.bth(i, b);
            cum += t;
            if t < self.min_unit { continue; }
            if let Some(f) = found.as_mut() { f.last = false; return Some(*f); }
            let u = Unit { k, top: prev, base: cum, first, last: true };
            if cum > d { found = Some(u); } else { last = Some(u); }
            prev = cum; first = false;
        }
        if found.is_some() { return found; }
        if last.is_some() { return last; }
        // every body thinner than min_unit here: the thickest one
        let k = (0..st.bodies.len()).max_by(|&a, &b| self.bth(i, &st.bodies[a]).total_cmp(&self.bth(i, &st.bodies[b])).then(b.cmp(&a)))?;
        Some(Unit { k, top: 0.0, base: cum, first: true, last: true })
    }

    /// All units of post i, top down.
    pub fn units(&self, i: usize) -> Vec<Unit> {
        let st = &self.stencils[self.st[i] as usize];
        let (mut cum, mut prev) = (0.0, 0.0);
        let mut out: Vec<Unit> = Vec::new();
        for (k, b) in st.bodies.iter().enumerate() {
            let t = self.bth(i, b);
            cum += t;
            if t < self.min_unit { continue; }
            out.push(Unit { k, top: prev, base: cum, first: out.is_empty(), last: false });
            prev = cum;
        }
        if out.is_empty() {
            if let Some(u) = self.unit_at(i, 0.0) { out.push(u); }
        }
        if let Some(l) = out.last_mut() { l.last = true; }
        out
    }

    /// Post i's position (east, north m) in body b's own frame: the local tangent plane at its representative layer's T0 cell
    /// centre (planet-keyed: the same for every chart, window and lattice).
    pub fn bxy(&self, i: usize, b: &Body) -> (f64, f64) {
        let (la, lo) = self.ll[i];
        let c = b.rep.0 as usize;
        let (dlat, dlon) = (180.0 / self.nlat as f64, 360.0 / self.nlon as f64);
        let (lac, loc) = (-90.0 + ((c / self.nlon) as f64 + 0.5) * dlat, -180.0 + ((c % self.nlon) as f64 + 0.5) * dlon);
        let dlo = (lo - loc + 540.0).rem_euclid(360.0) - 180.0;
        (self.radius_m * lac.to_radians().cos() * dlo.to_radians(), self.radius_m * (la - lac).to_radians())
    }
    #[inline] pub fn body(&self, i: usize, k: usize) -> &Body { &self.stencils[self.st[i] as usize].bodies[k] }

    /// (cell, full index) of the layer that records unit k at post i: the post's own cell's member, else the member with the
    /// largest weight.
    pub fn src(&self, i: usize, k: usize) -> Option<(u32, u16)> {
        let b = self.body(i, k);
        let w = &self.w[i];
        let st = &self.stencils[self.st[i] as usize];
        if let Some(s) = (0..4).find(|&s| st.cells[s] == self.own[i] && b.src[s].is_some()) { return b.src[s]; }
        (0..4).filter(|&s| b.src[s].is_some()).max_by(|&a, &c| w[a].total_cmp(&w[c]).then(c.cmp(&a))).and_then(|s| b.src[s])
    }

    /// Volcanic facies of post i at depth d: Some((slot j, facies, riser)) when the unit is a split body.
    pub fn volc_at(&self, i: usize, d: f64) -> Option<(u8, bool)> {
        let v = self.vrule.as_ref()?;
        let u = self.unit_at(i, d)?;
        let b = self.body(i, u.k);
        if b.vc == u32::MAX { return None; }
        let t = u.base - u.top;
        if !(t > 0.0) { return None; }
        let scale = b.rep.2 / t;
        let (cm, cb) = &self.vcols[b.vc as usize];
        let col = if u.last { cb } else { cm };
        if col.code.is_empty() { return None; }
        let (x, y) = self.bxy(i, b);
        let (_, f, riser) = col.facies_at((d - u.top) * scale, x, y, v.riser_min() * scale);
        Some((f, riser))
    }

    /// The post's column with split bodies expanded into their packages (post depths): (bases, codes, (cell, full index) of
    /// the recording layer, is-split per unit).
    pub fn expanded(&self, i: usize) -> (Vec<f64>, Vec<u8>, Vec<(u32, u16)>) {
        let us = self.units(i);
        let (mut base, mut code, mut src) = (Vec::new(), Vec::new(), Vec::new());
        for u in &us {
            let b = self.body(i, u.k);
            let s = self.src(i, u.k).unwrap_or((self.own[i], 0));
            if self.vrule.is_some() && b.vc != u32::MAX && u.base > u.top {
                let scale = b.rep.2 / (u.base - u.top);
                let (cm, cb) = &self.vcols[b.vc as usize];
                let col = if u.last { cb } else { cm };
                let mut any = false;
                for j in 0..col.code.len() {
                    if col.fac[j] == FAC_NONE { continue; }
                    if !u.first && col.base[j] <= 0.0 { continue; }
                    let bb = u.top + col.base[j] / scale;
                    base.push(if u.last { bb } else { bb.min(u.base) });
                    code.push(b.code); src.push(s); any = true;
                }
                if any { if !u.last { *base.last_mut().unwrap() = u.base; } continue; }
            }
            base.push(u.base); code.push(b.code); src.push(s);
        }
        (base, code, src)
    }

    /// volcanic_ct_v2 ladder at post i, depth d (see VOLCANIC_RULE): None = not a split unit or the split's ct off.
    pub fn ct_ladder(&self, i: usize, z: f64) -> Option<Option<([f64; 4], f64)>> {
        let v = self.vrule.as_ref()?;
        let datum = self.datum[i];
        let d = datum - z;
        let Some(sc) = self.ladder_col(i, z)? else { return Some(None) };
        Some(sc.ladder(d, 0.0, 0.0, v.riser_min()).and_then(|(l, lava)| {
            let (s1, s0) = (datum - l[1], datum - l[2]);
            if !(s1 > s0) { return None; }
            Some(([datum - l[3], s0, s1, datum - l[0]], (lava / (s1 - s0)).clamp(0.0, 1.0)))
        }))
    }

    /// Lane SPLIT-T1B volcanic_ct_blend_v1 at post i, height z (see t1b_volc::VolcPosts::ct_ladder_alts).
    pub fn ct_ladder_alts(&self, i: usize, z: f64) -> Option<Vec<crate::t1b_volc::LadderAlt>> {
        let v = self.vrule.as_ref()?;
        let datum = self.datum[i];
        let Some(sc) = self.ladder_col(i, z)? else { return Some(vec![(None, 1.0)]) };
        Some(match sc.ladder_alts(datum - z, 0.0, 0.0, v.riser_min(), v.band()) {
            None => vec![(None, 1.0)],
            Some(a) => crate::t1b_volc::heights_of(&a, datum),
        })
    }

    /// The contiguous run of split bodies around height z at post i as one split column in post depths (lens phases shifted
    /// to the post's position in each body's frame): None = not a split unit or the split's ct off; Some(None) = empty.
    fn ladder_col(&self, i: usize, z: f64) -> Option<Option<SplitCol>> {
        let v = self.vrule.as_ref()?;
        if v.ct != 1 { return None; }
        let datum = self.datum[i];
        let d = datum - z;
        let u0 = self.unit_at(i, d)?;
        if self.body(i, u0.k).vc == u32::MAX { return None; }
        let us = self.units(i);
        let p = us.iter().position(|u| u.k == u0.k)?;
        let (mut r0, mut r1) = (p, p);
        while r0 > 0 && self.body(i, us[r0 - 1].k).vc != u32::MAX { r0 -= 1; }
        while r1 + 1 < us.len() && self.body(i, us[r1 + 1].k).vc != u32::MAX { r1 += 1; }
        let mut sc = SplitCol { mid: v.lens_mid_permille as f64 / 1000.0, amp: v.lens_amp_permille as f64 / 1000.0, ..Default::default() };
        for u in &us[r0..=r1] {
            let b = self.body(i, u.k);
            if !(u.base > u.top) { continue; }
            let scale = b.rep.2 / (u.base - u.top);
            let (cm, cb) = &self.vcols[b.vc as usize];
            let col = if u.last { cb } else { cm };
            for j in 0..col.code.len() {
                if col.fac[j] == FAC_NONE { continue; }
                if !u.first && col.base[j] <= 0.0 { continue; }
                let top = u.top + col.top[j] / scale;
                let base = if u.last { u.top + col.base[j] / scale } else { (u.top + col.base[j] / scale).min(u.base) };
                if !(base > top) { continue; }
                // lens phases shifted so the ladder's evaluation at (0, 0) reads each body at the post's position in its own frame
                let (bx, by) = self.bxy(i, b);
                let mut l = col.lens[j];
                for qd in 0..3 { l[3 * qd + 2] += l[3 * qd] * bx + l[3 * qd + 1] * by; }
                sc.top.push(top); sc.base.push(base); sc.code.push(col.code[j]); sc.idx.push(col.idx[j]); sc.fac.push(col.fac[j]); sc.lens.push(l);
            }
        }
        if sc.code.is_empty() { return Some(None); }
        Some(Some(sc))
    }

    /// Lane SPLIT-T1B volcanic_aa_v1: (slope, top depth) of the piece just above depth d at post i: inside a split body its
    /// lava / volcaniclastic piece (body depth mapped as `volc_at`), else the unit's class slope up to the unit top.
    fn piece_up(&self, i: usize, d: f64, v: &VolcanicSplitV1, sc: &[f64; 12]) -> (f64, f64) {
        let dq = d - 1e-6;
        let Some(u) = self.unit_at(i, dq) else { return (sc[11], f64::NEG_INFINITY) };
        let b = self.body(i, u.k);
        let top_u = if u.first { f64::NEG_INFINITY } else { u.top };
        let t = u.base - u.top;
        if b.vc != u32::MAX && t > 0.0 {
            let scale = b.rep.2 / t;
            let (cm, cb) = &self.vcols[b.vc as usize];
            let col = if u.last { cb } else { cm };
            if !col.code.is_empty() {
                let (x, y) = self.bxy(i, b);
                let (s, pe) = col.piece_up((d - u.top) * scale, x, y, v.riser_min() * scale, v.band() * scale, v.mass(), v.face(), &|_| v.mass());
                return (s, (u.top + pe / scale).max(top_u));
            }
        }
        (sc[(b.code as usize).min(11)], top_u)
    }

    /// (critical slope, K multiplier, code) at post i, height z; `table` the code table; `fallback` the stand-in basement code.
    #[inline]
    pub fn at(&self, i: usize, z: f64, sc: &[f64; 12], k: &[f64; 12], fallback: u8) -> (f64, f64, u8) {
        let d = self.datum[i] - z;
        let Some(u) = self.unit_at(i, d) else { let c = (fallback as usize).min(11); return (sc[c], k[c], fallback); };
        let code = self.body(i, u.k).code;
        if let (Some((f, riser)), Some(v)) = (self.volc_at(i, d), self.vrule.as_ref()) {
            // Lane SPLIT-T1B volcanic_aa_v1: the exact climb rate over aa_dx through the post's units (K read at d).
            if v.aa.is_some() {
                let s = crate::t1b_volc::climb_rate(d, v.aa_dx(), &|dd| self.piece_up(i, dd, v, sc));
                return (s, if f == FAC_LAVA { v.lava_k_permille as f64 / 1000.0 } else { v.vc_k_permille as f64 / 1000.0 }, code);
            }
            return match f {
                FAC_LAVA => (if riser { v.face() } else { v.mass() }, v.lava_k_permille as f64 / 1000.0, code),
                _ => (v.mass(), v.vc_k_permille as f64 / 1000.0, code),
            };
        }
        let c = (code as usize).min(11);
        (sc[c], k[c], code)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 3 x 3 cells; columns differ so that edges need alignment and pinch-outs.
    fn ident() -> Identity {
        let (nlat, nlon) = (3, 4);
        let c = |v: &[(f64, u8)]| -> Vec<(f64, u8, u16)> { v.iter().enumerate().map(|(k, &(t, x))| (t, x, k as u16)).collect() };
        let mut cols = vec![c(&[(100.0, 5), (500.0, 3), (2000.0, 1), (8000.0, 11)]); nlat * nlon];
        cols[1] = c(&[(300.0, 3), (50.0, 6), (2500.0, 1), (9000.0, 11)]);
        cols[5] = c(&[(40.0, 6), (700.0, 3), (100.0, 4), (1500.0, 1), (7000.0, 11)]);
        cols[6] = c(&[(900.0, 3), (3000.0, 2), (7000.0, 11)]);
        let full = cols.iter().map(|v| v.len() as u32).collect();
        Identity::from_columns(nlat, nlon, &cols, full, &vec![true; nlat * nlon])
    }

    #[test]
    fn alignment_is_monotone_and_same_code() {
        let a = [(100.0, 5), (500.0, 3), (2000.0, 1)];
        let b = [(300.0, 3), (50.0, 6), (2500.0, 1)];
        assert_eq!(align(&a, &b), vec![(1, 0), (2, 2)]);
    }

    #[test]
    fn fixpoint_leaves_every_quad_consistent_and_edges_agree() {
        let id = ident();
        for q in 0..(2 * 4) as u32 { assert!(id.check(q).is_some(), "quad {q}"); }
        // The two quads sharing the edge cells (1, 5) order the bodies of 1 and 5 identically.
        let restrict = |q: u32, a: u32, b: u32| -> Vec<(Option<u16>, Option<u16>)> {
            let cells = id.quad_cells(q);
            let sa = cells.iter().position(|&c| c == a).unwrap();
            let sb = cells.iter().position(|&c| c == b).unwrap();
            id.check(q).unwrap().iter().map(|m| (m[sa], m[sb])).filter(|p| p.0.is_some() || p.1.is_some()).collect()
        };
        // cell 1 = (0,1), cell 5 = (1,1): quads (0,0) [slots b,d] and (0,1) [slots a,c]
        assert_eq!(restrict(0, 1, 5), restrict(1, 1, 5));
    }

    #[test]
    fn columns_agree_on_a_quad_edge_line() {
        let id = ident();
        // Emulate two posts on the line between cells 1 and 5 seen from quads 0 and 1 (weights only on 1 and 5).
        let mk = |q: u32, wa: f64| -> Vec<(f64, u8)> {
            let cells = id.quad_cells(q);
            let order = id.check(q).unwrap();
            let mut w = [0.0; 4];
            for s in 0..4 { if cells[s] == 1 { w[s] = wa; } if cells[s] == 5 { w[s] = 1.0 - wa; } }
            let mut cum = 0.0;
            let mut out = Vec::new();
            for m in &order {
                let mut t = 0.0; let mut code = 0;
                for s in 0..4 { if let Some(a) = m[s] { let g = id.g(cells[s], a); t += w[s] * id.kth[g]; code = id.kcode[g]; } }
                cum += t;
                if t >= 1.0 { out.push((cum, code)); }
            }
            out
        };
        for wa in [0.0, 0.3, 0.7, 1.0] { assert_eq!(mk(0, wa), mk(1, wa), "wa {wa}"); }
    }

    // ---------------------------------------------------------- lane LITHO-AGE-GATE

    /// 2 x 4 cells, only the listed cells land; per body (thickness, code, record).
    fn gated(cols: &[(usize, Vec<(f64, u8, Sig)>)], et: Option<u32>) -> Identity {
        let (nlat, nlon) = (2, 4);
        let mut cc = vec![Vec::new(); nlat * nlon];
        let mut land = vec![false; nlat * nlon];
        let mut sg = vec![Vec::new(); nlat * nlon];
        for (c, v) in cols {
            cc[*c] = v.iter().enumerate().map(|(k, x)| (x.0, x.1, k as u16)).collect();
            sg[*c] = v.iter().map(|x| x.2).collect();
            land[*c] = true;
        }
        let sig: Vec<Sig> = sg.into_iter().flatten().collect();
        let full = cc.iter().map(|v| v.len() as u32).collect();
        match et { Some(e) => Identity::from_columns_g(nlat, nlon, &cc, full, &land, Some((&sig, e))), None => Identity::from_columns(nlat, nlon, &cc, full, &land) }
    }
    fn same(id: &Identity, a: (usize, usize), b: (usize, usize)) -> bool { id.root[id.koff[a.0] as usize + a.1] == id.root[id.koff[b.0] as usize + b.1] }

    #[test]
    fn age_gate_splits_a_false_join_of_two_basalt_events() {
        // site 1 pattern: 398-400 My basalt beside 286-324 My basalt, nothing between them; same code, same order position.
        let bas = |a, b| Sig::layer(a, b, 2, 5.0);
        let gra = Sig::NONE;
        let cols = vec![(0, vec![(36.0, 3, bas(398.0, 400.0)), (20000.0, 1, gra)]), (1, vec![(3686.0, 3, bas(286.0, 324.0)), (23000.0, 1, gra)])];
        let off = gated(&cols, None);
        assert!(same(&off, (0, 0), (1, 0)), "lcont joins them by code and order");
        let on = gated(&cols, Some(1));
        assert!(!same(&on, (0, 0), (1, 0)), "the gate splits them");
        assert!(same(&on, (0, 1), (1, 1)), "the unrecorded granite is still one body");
        assert_eq!(on.gate.as_ref().unwrap().alignment_pairs_gated, 1);
    }

    #[test]
    fn age_gate_lets_the_alignment_pick_the_right_event() {
        // left: young basalt over old basalt (split at a recorded unconformity); right: old basalt only.
        let bas = |a, b| Sig::layer(a, b, 2, 5.0);
        let cols = vec![(0, vec![(500.0, 3, bas(398.0, 400.0)), (100.0, 3, bas(286.0, 324.0)), (9000.0, 1, Sig::NONE)]),
                        (1, vec![(800.0, 3, bas(290.0, 320.0)), (9000.0, 1, Sig::NONE)])];
        let on = gated(&cols, Some(1));
        assert!(same(&on, (0, 1), (1, 0)) && !same(&on, (0, 0), (1, 0)));
        for q in 0..4 { assert!(on.check(q).is_some()); }
    }

    #[test]
    fn env_gate_uses_facies_tract_distance() {
        let sh = |e| Sig::layer(390.0, 400.0, e, 5.0);
        let pair = |e1, e2| { let id = gated(&vec![(0, vec![(300.0, 5, sh(e1))]), (1, vec![(300.0, 5, sh(e2))])], Some(1)); same(&id, (0, 0), (1, 0)) };
        assert!(!pair(4, 6), "lacustrine and shelf shale are two bodies");
        assert!(pair(6, 7), "shelf and deep-marine shale are one body");
        assert!(!pair(2, 6), "a volcanic and a depositional environment never join");
    }

    #[test]
    fn union_gate_is_pairwise_for_every_member() {
        // A [0, 10], B [8, 20], C [18, 30], tolerance 0: A-B and B-C overlap, A-C does not, so A, B, C are never one body.
        let s = |a, b| Sig::layer(a, b, 6, 0.0);
        let cols = vec![(0, vec![(100.0, 4, s(0.0, 10.0))]), (1, vec![(100.0, 4, s(8.0, 20.0))]), (2, vec![(100.0, 4, s(18.0, 30.0))])];
        let off = gated(&cols, None);
        assert!(same(&off, (0, 0), (2, 0)));
        let on = gated(&cols, Some(1));
        assert!(!same(&on, (0, 0), (2, 0)));
        assert!(same(&on, (0, 0), (1, 0)) ^ same(&on, (1, 0), (2, 0)));
        assert_eq!(on.gate.as_ref().unwrap().rejected_union_age, 1);
    }

    #[test]
    fn runs_split_at_a_recorded_unconformity_only_with_the_gate() {
        // one cell, full column: basalt 36 m (398-400, unconformity at base), basalt 3686 m (286-324), granite (no record)
        let off = [0u32, 3]; let lith = [3u8, 3, 1]; let th = [36.0f32, 3686.0, 20000.0];
        let js = r#"{"schema":"pf.stratigraphy/v0","unconformity_gap_my":5.0,"layer_offset":[0,2],"lith":[3,3],"env":[2,2],"flags":[8,16],"thickness_m":[36.0,2964.0],"age_base_my":[398.0,286.0],"age_top_my":[400.0,324.0]}"#;
        let rec = Record::from_json_bytes(js.as_bytes(), &off, &lith, &th).unwrap();
        assert!(rec.known[1] && !rec.known[2] && rec.flags[1] == 0);
        let base = [36.0, 3722.0, 23722.0]; let idx = [0u16, 1, 2];
        let (b0, s0, _) = cell_bodies(&base, &lith, &idx, 0, None);
        assert_eq!(b0.len(), 2);
        assert!(s0.is_empty());
        let (b1, s1, n1) = cell_bodies(&base, &lith, &idx, 0, Some((&rec, 5.0, 1)));
        assert_eq!(b1.iter().map(|b| (b.1, b.2)).collect::<Vec<_>>(), vec![(3, 0), (3, 1), (1, 2)]);
        assert_eq!(n1, [2, 1, 0]);
        assert_eq!(s1[2], Sig::NONE);
        // a record that disagrees with the column is refused
        let bad = js.replace("\"lith\":[3,3]", "\"lith\":[3,5]");
        assert!(Record::from_json_bytes(bad.as_bytes(), &off, &lith, &th).is_err());
    }

    #[test]
    fn record_key_absent_is_not_serialized() {
        let v = serde_json::to_value(LithoContinuityV1::v1()).unwrap();
        assert!(v.get("record").is_none());
    }
    fn m2o_fixture(change: &str) -> (Vec<Vec<(f64, u8, u16)>>, Vec<Sig>, M2oStats) {
        let mut c = vec![vec![]; 8];
        c[0] = vec![(600.0, 3, 0), (1000.0, 1, 1)];
        c[1] = vec![(100.0, 3, 0), (20.0, 5, 1), (200.0, 3, 2), (1000.0, 1, 3)];
        let b = |lo, hi| Sig::layer(lo, hi, 2, 5.0);
        let mut sg = vec![b(370.0, 400.0), Sig::NONE, b(394.0, 400.0), Sig::layer(392.0, 394.0, 6, 5.0), b(390.0, 394.0), Sig::NONE];
        let mut bd = vec![vec![]; 8]; bd[0] = vec![false; 2]; bd[1] = vec![false; 4];
        match change {
            "unconformity" => bd[1][1] = true,
            "age" => sg[4] = b(360.0, 370.0),
            "environment" => sg[4] = Sig::layer(390.0, 394.0, 6, 5.0),
            "unknown" => sg[4] = Sig::NONE,
            "order" => sg[4] = b(396.0, 402.0),
            "thickness" => c[0][0].0 = 310.0,
            _ => (),
        }
        let st = refine_m2o(2, 4, &mut c, &mut sg, &bd, &[true, true, false, false, false, false, false, false], 1, 5.0);
        (c, sg, st)
    }

    #[test]
    fn m2o_conserves_column_and_keeps_interbed_with_distinct_ordered_identities() {
        let (c, sg, st) = m2o_fixture("");
        assert_eq!(st.selected, 1); assert_eq!(st.added_slices, 1);
        assert_eq!(c[0], vec![(120.0, 3, 0), (480.0, 3, 0), (1000.0, 1, 1)]);
        assert_eq!(st.max_column_thickness_error_m, 0.0);
        let id = Identity::from_columns_g(2, 4, &c, vec![2, 4, 0, 0, 0, 0, 0, 0],
            &[true, true, false, false, false, false, false, false], Some((&sg, 1)));
        assert!(same(&id, (0, 0), (1, 0))); assert!(same(&id, (0, 1), (1, 2)));
        assert!(!same(&id, (0, 0), (0, 1)));
        assert_eq!(id.kth[id.g(1, 1)], 20.0); // shale is preserved in its actual source cell
        for q in 0..4 { assert!(id.check(q).is_some()); }
        // Age gate still applies to every union after refinement.
        for a in 0..id.root.len() { for b in 0..id.root.len() {
            if id.root[a] == id.root[b] { let m = sg[a].meet(&sg[b]); assert!(m.age_ok() && m.env_ok(1)); }
        } }
    }

    #[test]
    fn m2o_refuses_unconformity_age_environment_unknown_reversal_and_excess_span() {
        for reason in ["unconformity", "age", "environment", "unknown", "order", "thickness"] {
            let (c, _, st) = m2o_fixture(reason);
            assert_eq!(st.selected, 0, "{reason}"); assert_eq!(c[0].len(), 2, "{reason}");
        }
    }

    #[test]
    fn m2o_key_off_keeps_record_recipe_bytes() {
        let js = r#"{"digest":"blake3:test","env_tol":1,"rule":"old"}"#;
        let r: LcontRecordV1 = serde_json::from_str(js).unwrap();
        assert!(r.m2o_rule.is_none());
        assert_eq!(serde_json::to_string(&r).unwrap(), js);
    }

    #[test]
    fn m2o_target_stack_can_anchor_both_neighbouring_thick_columns() {
        let mut c = vec![vec![]; 8];
        c[0] = vec![(600.0, 3, 0)]; c[2] = c[0].clone();
        c[1] = vec![(100.0, 3, 0), (20.0, 5, 1), (200.0, 3, 2)];
        let b = |lo, hi| Sig::layer(lo, hi, 2, 5.0);
        let mut sg = vec![b(370.0, 400.0), b(394.0, 400.0), Sig::layer(392.0, 394.0, 6, 5.0), b(390.0, 394.0), b(370.0, 400.0)];
        let bd: Vec<Vec<bool>> = c.iter().map(|v| vec![false; v.len()]).collect();
        let land = [true, true, true, false, false, false, false, false];
        let st = refine_m2o(2, 4, &mut c, &mut sg, &bd, &land, 1, 5.0);
        assert_eq!(st.selected, 2);
        let id = Identity::from_columns_g(2, 4, &c, vec![1, 3, 1, 0, 0, 0, 0, 0], &land, Some((&sg, 1)));
        for a in [0, 1] { assert!(same(&id, (0, a), (2, a))); }
        for q in 0..4 { assert!(id.check(q).is_some()); }
    }

    #[test]
    fn m2o_time_partition_is_conservative_and_zero_duration_is_not_guessed() {
        let (c, sg, _) = m2o_fixture("");
        assert_eq!((sg[0].lo, sg[0].hi), (394.0, 405.0));
        assert_eq!((sg[1].lo, sg[1].hi), (370.0, 399.0));
        assert_eq!(c[0][0].0 + c[0][1].0, 600.0);
        let mut cc = vec![vec![]; 8];
        cc[0] = vec![(600.0, 3, 0)]; cc[1] = vec![(100.0, 3, 0), (20.0, 5, 1), (200.0, 3, 2)];
        let zero = Sig::layer(400.0, 400.0, 2, 5.0);
        let mut ss = vec![zero, zero, Sig::layer(400.0, 400.0, 6, 5.0), zero];
        let bd = cc.iter().map(|v| vec![false; v.len()]).collect::<Vec<_>>();
        let st = refine_m2o(2, 4, &mut cc, &mut ss, &bd, &[true, true, false, false, false, false, false, false], 1, 5.0);
        assert_eq!(st.selected, 0); assert_eq!(cc[0].len(), 1);
    }

}

/// E3 adapts source envelopes to the existing rock-column interface without invoking
/// lcont's private matching, absorption or interpolation. A stencil has ONE owner.
pub fn envelope_rock<'a>(cfg: &crate::t1b_litho::EnvelopeConfig, table: &'a crate::t1b_litho::LithologyV1,
    t0: &T0, cols: &'a StrataCols, pos: &[(f64,f64,f64,f64)]) -> Result<crate::t1b_rivernet::RockPosts<'a>, String> {
    let env = cfg.load()?;
    env.validate_inputs(t0, cols)?;
    let mut cp = ContPosts { min_unit: 0., stencils: Vec::new(), st: Vec::with_capacity(pos.len()), w: Vec::with_capacity(pos.len()),
        datum: Vec::with_capacity(pos.len()), grad: Vec::with_capacity(pos.len()), ll: Vec::with_capacity(pos.len()), radius_m: t0.radius_m,
        nlat: t0.nlat, nlon: t0.nlon, own: Vec::with_capacity(pos.len()), vrule: None, vcols: Vec::new(), stats: ContStats::default() };
    let mut st = BTreeMap::new();
    for &(la,lo,_,_) in pos {
        let c = t0.cell_of_deg(la,lo);
        let e = env.cells.get(&c);
        if t0.water(c) == CellWater::Land && e.is_none() { return Err(format!("env missing {c}")); }
        let si = *st.entry(c).or_insert_with(|| {
            // Nonland is outside E1 support. Read its actual full strata, never invent a fallback unit.
            let units: Vec<(u16,f64,u8)> = match e {
                Some(e) => e.layers.iter().filter(|l| l.h>0.).map(|l| (l.index,l.h,l.code)).collect(),
                None => (cols.off[c as usize]..cols.off[c as usize+1]).filter(|&k| cols.th[k as usize]>0.).map(|k| ((k-cols.off[c as usize]) as u16,cols.th[k as usize] as f64,cols.lith[k as usize])).collect(),
            };
            let bodies = units.into_iter().map(|(k,h,code)| Body { th:[h,0.,0.,0.],code,src:[Some((c,k)),None,None,None],rep:(c,k,h),root:cols.off[c as usize]+k as u32,vc:u32::MAX }).collect();
            cp.stencils.push(Stencil { cells:[c;4],bodies }); (cp.stencils.len()-1) as u32
        });
        if cp.stencils[si as usize].bodies.is_empty() { return Err(format!("env adapter: no recorded rock in nonland cell {c}")); }
        let (datum,grad) = if let Some(e)=e { (e.mean_m+e.tilt(c,la,lo),(e.grad[0],e.grad[1])) } else { (t0.elevation_m[c as usize] as f64,(0.,0.)) };
        cp.st.push(si);cp.w.push([1.,0.,0.,0.]);cp.datum.push(datum);cp.grad.push(grad);cp.ll.push((la,lo));cp.own.push(c);
    }
    cp.stats.posts=pos.len() as u64;cp.stats.stencils=cp.stencils.len() as u64;
    let mut sc=[0.;12];let mut k=[0.;12];for j in 0..12 {sc[j]=table.sc(j as u8);k[j]=table.k(j as u8);}
    Ok(crate::t1b_rivernet::RockPosts::Strata {table,cols,cell:cp.own.clone(),datum:cp.datum.clone(),fallback:vec![255;pos.len()],sc,k,volc:None,cont:Some(Box::new(cp))})
}
