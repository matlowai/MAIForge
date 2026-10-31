//! Lane CHANNEL-CENTRELINE (2026-10-21): opt-in centreline of the T1a channels in T1b RIVERNET windows.
//! LABEL: EXPLORATORY, window-only, every constant a project assumption awaiting review. Absent (`params.cl` None) =
//! byte-identical.
//!
//! Why (lanes VALLEY-FLOOR and DELTA-FLOOR-LOOK): a T1a link is drawn as a Catmull-Rom curve through its two 512 m
//! lattice nodes, and the T1a node chain is the 512 m rasterization (`skeleton::rasterize`, round to the nearest node,
//! 8-connected, taut) of the T1a keyed planform curve. Collinear D8 nodes give straight reaches at exactly 0 or 45 deg
//! and the rasterization adds 0/45 deg jogs. Neither is a property of the T1a valley axis: it is lattice quantization.
//!
//! Two switches (keys of pf-region-relief rnwin; both window-only, the layer refuses them):
//! * `cl=1`, CURVE_RULE: every T1a link follows the T1a planform curve itself, re-evaluated from the segment records
//!   published in the digested T1a header (no draws). This is the T1a path F1 refers to; the lattice nodes are its
//!   rasterization. It is NOT a physical meander model (the T1a planform is T1a's declared keyed stand-in).
//! * `cl_mn=1` (needs cl=1), MEANDER_RULE: a DECLARED DECORATIVE CLOSURE, keyed and labelled, explicitly NOT meander
//!   evolution (AGENTS.md PHYSICS-FIRST; integrator ruling 2026-10-15 06:02; real meanders are RIVER-PILOT /
//!   RIVER-INTEGRATE): a sine-generated meander belt on the T1a curve, confined to the valley floor.
//! * `cl_mq=1` (lane MEANDER-Q, needs cl=1, exclusive with cl_mn): the sine belt is replaced by a kinematic meander-migration
//!   model with cutoffs and a swept-envelope floodplain (t1b_mq.rs).
//!
//! Why not the steepest-descent path on the 32 m relief (the lane's first option): T1b RIVERNET is a river-network-first
//! construction (module doc of t1b_rivernet.rs): the valley floor is the union of discs around the drawn channel posts and
//! the hillslopes are an eikonal from them, so the 32 m surface has no valley axis independent of the curve; a descent
//! path on it reproduces the curve (circular). The only independent sub-node planform information is the T1a curve.

use serde::{Deserialize, Serialize};
use serde_json::Value;

pub const LABEL: &str = "EXPLORATORY lane CHANNEL-CENTRELINE 2026-10-21; window-only; constants project assumptions awaiting review";
pub const CURVE_RULE: &str = "t1a_curve_centreline_v1: each T1a segment's planform curve is re-evaluated from its published record (skeleton.rs curve(): P(t) = start + t (end - start) + amplitude v(t) n, v = smoothstep between the control offsets (published in permille, so offsets are within 0.0005 amplitude of the T1a draw), n the left unit normal of the chord), sampled at chord steps <= 4 m; every T1a node of a segment is projected to its nearest curve sample (coarse 128 m scan, then +-256 m fine scan; downstream monotone, a non-monotone link is drawn straight and counted); a T1a link n -> r of one segment is the curve between the two projections; a link into another segment (junction) or a node without a segment follows its own curve to the sample nearest the receiver's projected point and ends there with a straight connector (reported); beds are the T1b node levels at both ends, linear in path length; network channels and conduits that end or start at a T1a node use its projected point. F1 with the key reads 'every T1a channel node has a channel post of its own link' instead of 'its lattice post is a channel post' (the node is the 512 m rasterization of the curve, up to ~360 m from it)";
pub const MEANDER_RULE: &str = "sine_generated_meander_belt_v1 (DECLARED DECORATIVE CLOSURE, keyed, NOT meander evolution: no flow, bank erosion, deposition or cutoff is modeled; stands in until RIVER-PILOT / RIVER-INTEGRATE): along each T1a segment curve (arc length s) the channel is offset by y(s) = T(s) R(lambda) lambda yhat(omega, phi) along the curve normal; yhat is the sine-generated curve of unit valley wavelength (direction theta = omega cos(2 pi s' / M); Langbein and Leopold 1966 USGS PP 422-H, recalled: sinuosity = 1 / J0(omega)); wavelength lambda = c_lambda w (Leopold and Wolman 1960: lambda ~ 10.9 w^1.01, web-verified abstract, earth-numbers-20261020 C; c_lambda = 11), bankfull width w = the valley-floor-fix hydraulic geometry 3.5 Q^0.5 (pf-detail B6, Q = A r_ref / yr) of the link's T1a area; phase phi = integral ds / lambda; sinuosity k keyed per T1a node in [k_min, k_max] (the existing t1b.rivernet.meander draw of the node, linear between nodes; 1.3-2.0 from the lane brief, meandering rivers about 1.5 and above), capped so the belt half-width lambda max(yhat) <= f_belt (h - w / 2) with h the T1b floor half-width, or near a volcano edifice (curve point within r_influence + W/2 max of a vent) the corridor cap's mask_inner (the edifice flank buries the floor beyond it), wavelength smoothed by a moving average and room by a running minimum + moving average over half a wavelength each side (continuous offsets); T = smoothstep over one wavelength at both segment ends (junctions, mouths); R = smoothstep((lambda - 256 m) / 256 m) (bends shorter than 16 T1b posts are below the lattice resolution, left to finer tiers); the floor extent stays the band around the unmeandered curve (the valley walls do not follow the bends), floor values come from the nearest meandered channel post";

pub const WAVELENGTH_WIDTHS_PERMILLE: i64 = 11_000;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ClParamsV1 {
    pub on: u32,
    pub meander: u32,
    pub wavelength_widths_permille: i64,
    pub sinuosity_min_permille: i64,
    pub sinuosity_max_permille: i64,
    pub belt_fill_permille: i64,
    pub curve_rule: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub meander_rule: String,
    pub label: String,
    /// Lane MEANDER-Q (opt-in cl_mq=1, t1b_mq.rs): kinematic meander migration replacing the declared sine belt; None = absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mq: Option<crate::t1b_mq::MqParamsV1>,
}

impl ClParamsV1 {
    pub fn v1() -> Self {
        ClParamsV1 { on: 1, meander: 0, wavelength_widths_permille: WAVELENGTH_WIDTHS_PERMILLE, sinuosity_min_permille: 1300, sinuosity_max_permille: 2000,
            belt_fill_permille: 700, curve_rule: CURVE_RULE.into(), meander_rule: String::new(), label: LABEL.into(), mq: None }
    }
    pub fn with_meander(mut self) -> Self { self.meander = 1; self.meander_rule = MEANDER_RULE.into(); self }
    pub fn validate(&self) -> Result<(), String> {
        if self.on != 1 { return Err("cl: on must be 1 (absent = off)".into()); }
        if self.meander > 1 { return Err("cl_mn is 0 or 1".into()); }
        if !(1000..=2000).contains(&self.sinuosity_min_permille) || !(self.sinuosity_min_permille..=2050).contains(&self.sinuosity_max_permille) {
            return Err("cl sinuosity range must lie in [1.0, 2.05] (the sine-generated belt needs omega < pi/2)".into());
        }
        if !(1000..=30_000).contains(&self.wavelength_widths_permille) || !(100..=1000).contains(&self.belt_fill_permille) { return Err("cl constants".into()); }
        if let Some(m) = &self.mq { m.validate()?; if self.meander == 1 { return Err("cl_mq (physical planform) and cl_mn (decorative sine belt) are exclusive".into()); } }
        Ok(())
    }
    /// A planform belt (the declared sine belt cl_mn or lane MEANDER-Q's migration cl_mq) is on: the volcano corridor confines it.
    pub fn meander_on(&self) -> bool { self.meander == 1 || self.mq.is_some() }
    /// The declared decorative sine belt (cl_mn=1).
    pub fn sine_on(&self) -> bool { self.meander == 1 }
    /// Lane MEANDER-Q kinematic migration (cl_mq=1).
    pub fn mq_on(&self) -> bool { self.mq.is_some() }
}

pub fn cl(p: &crate::t1b_rivernet::RiverNetParamsV1) -> Option<&ClParamsV1> { p.cl.as_ref().filter(|c| c.on == 1) }

/// One T1a segment's planform curve as published (chart metres).
#[derive(Clone, Debug)]
pub struct T1aSegCurve { pub id: u32, pub start: [f64; 2], pub end: [f64; 2], pub amp: f64, pub ctrl: Vec<f64>,
    /// Lane MEANDER-Q: the segment's world address (seg/<from>-<to> | head/<cell>, T0 L-grid element ids) and the world seed of
    /// the T1a recipe (its seed_scope), for world-keyed draws; empty / None when the header lacks them.
    pub address: String, pub seed: Option<u64> }

impl T1aSegCurve {
    pub fn from_header(h: &Value) -> Result<Vec<T1aSegCurve>, String> {
        let Some(a) = h.get("segments").and_then(Value::as_array) else { return Ok(Vec::new()) };
        let seed = h.pointer("/recipe/seed_scope/seed").and_then(Value::as_u64);
        let mut out = Vec::with_capacity(a.len());
        for (i, s) in a.iter().enumerate() {
            let id = s.get("id").and_then(Value::as_u64).ok_or("T1a segment id")? as u32;
            if id as usize != i { return Err(format!("T1a segment {i} has id {id} (ids must equal positions)")); }
            let pt = |k: &str| -> Result<[f64; 2], String> {
                let v = s.get(k).and_then(Value::as_array).ok_or_else(|| format!("T1a segment {k}"))?;
                Ok([v[0].as_i64().ok_or("mm")? as f64 / 1000.0, v[1].as_i64().ok_or("mm")? as f64 / 1000.0])
            };
            let ctrl = s.get("control_offsets_permille").and_then(Value::as_array).ok_or("T1a segment controls")?
                .iter().map(|v| v.as_i64().map(|x| x as f64 / 1000.0).ok_or("control")).collect::<Result<Vec<f64>, _>>()?;
            out.push(T1aSegCurve { id, start: pt("start_mm")?, end: pt("end_mm")?, amp: s.get("amplitude_mm").and_then(Value::as_i64).ok_or("amplitude")? as f64 / 1000.0, ctrl,
                address: s.get("address").and_then(Value::as_str).unwrap_or("").to_string(), seed });
        }
        Ok(out)
    }
    pub fn chord(&self) -> f64 { (self.end[0] - self.start[0]).hypot(self.end[1] - self.start[1]) }
    /// skeleton.rs curve() at parameter t (metres).
    pub fn at(&self, t: f64) -> [f64; 2] {
        let d = [self.end[0] - self.start[0], self.end[1] - self.start[1]];
        let len = self.chord();
        if len == 0.0 || self.ctrl.len() < 2 { return [self.start[0] + t * d[0], self.start[1] + t * d[1]]; }
        let n_int = self.ctrl.len() - 1;
        let pos = t * n_int as f64;
        let k = (pos.floor().max(0.0) as usize).min(n_int - 1);
        let f = pos - k as f64;
        let sm = f * f * (3.0 - 2.0 * f);
        let v = self.ctrl[k] + (self.ctrl[k + 1] - self.ctrl[k]) * sm;
        let off = self.amp * v;
        [self.start[0] + t * d[0] - off * d[1] / len, self.start[1] + t * d[1] + off * d[0] / len]
    }
}

/// Sampled segment: positions, arc length, unit left normals; `y` the meander offset (empty without the meander).
pub struct SegTab { pub x: Vec<[f64; 2]>, pub s: Vec<f64>, pub nrm: Vec<[f64; 2]>, pub y: Vec<f64> }

pub const SAMPLE_M: f64 = 4.0;

impl SegTab {
    pub fn new(c: &T1aSegCurve) -> SegTab {
        let ns = ((c.chord() / SAMPLE_M).ceil() as usize).max(1);
        let x: Vec<[f64; 2]> = (0..=ns).map(|q| c.at(q as f64 / ns as f64)).collect();
        let mut s = vec![0.0; x.len()];
        for k in 1..x.len() { s[k] = s[k - 1] + (x[k][0] - x[k - 1][0]).hypot(x[k][1] - x[k - 1][1]); }
        let nrm = (0..x.len()).map(|k| {
            let (a, b) = (x[k.saturating_sub(1)], x[(k + 1).min(x.len() - 1)]);
            let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
            let l = dx.hypot(dy).max(1e-12);
            [-dy / l, dx / l]
        }).collect();
        SegTab { x, s, nrm, y: Vec::new() }
    }
    /// Nearest sample to p: coarse scan every 32 samples (128 m), then +-64 samples.
    pub fn project(&self, p: [f64; 2]) -> usize {
        let d2 = |k: usize| { let (dx, dy) = (self.x[k][0] - p[0], self.x[k][1] - p[1]); dx * dx + dy * dy };
        let mut best = 0usize;
        let mut k = 0usize;
        while k < self.x.len() { if d2(k) < d2(best) { best = k; } k += 32; }
        if d2(self.x.len() - 1) < d2(best) { best = self.x.len() - 1; }
        let (lo, hi) = (best.saturating_sub(64), (best + 64).min(self.x.len() - 1));
        let mut b2 = best;
        for k in lo..=hi { if d2(k) < d2(b2) { b2 = k; } }
        b2
    }
    /// Nearest sample to p among k >= k0.
    pub fn project_from(&self, p: [f64; 2], k0: usize) -> usize {
        let d2 = |k: usize| { let (dx, dy) = (self.x[k][0] - p[0], self.x[k][1] - p[1]); dx * dx + dy * dy };
        let mut best = k0;
        for k in k0..self.x.len() { if d2(k) < d2(best) { best = k; } }
        best
    }
    /// The channel point at sample k (meandered when the offsets exist).
    pub fn pt(&self, k: usize) -> [f64; 2] {
        match self.y.get(k) { Some(&y) => [self.x[k][0] + y * self.nrm[k][0], self.x[k][1] + y * self.nrm[k][1]], None => self.x[k] }
    }
}

// ------------------------------------------------------------------ sine-generated curve tables

const OMEGA_STEPS: usize = 160; // omega grid 0 .. 1.59 rad by 0.01
const X_STEPS: usize = 512;

/// yhat(x) for x in [0, 1) on X_STEPS points, valley wavelength 1, direction theta = omega cos(2 pi s / M); and its
/// sinuosity M / 1 = 1 / J0(omega).
fn sg_shape(omega: f64) -> (Vec<f64>, f64) {
    let m = 8192usize;
    let (mut xs, mut ys) = (vec![0.0; m + 1], vec![0.0; m + 1]);
    for q in 0..m {
        let th = |u: f64| omega * (2.0 * std::f64::consts::PI * u).cos();
        let (u0, u1) = (q as f64 / m as f64, (q + 1) as f64 / m as f64);
        let tm = th(0.5 * (u0 + u1));
        xs[q + 1] = xs[q] + tm.cos() / m as f64;
        ys[q + 1] = ys[q] + tm.sin() / m as f64;
    }
    let adv = xs[m];
    let mut out = vec![0.0; X_STEPS];
    let mut q = 0usize;
    for (j, o) in out.iter_mut().enumerate() {
        let xt = j as f64 / X_STEPS as f64 * adv;
        while q + 1 < m && xs[q + 1] < xt { q += 1; }
        let f = if xs[q + 1] > xs[q] { (xt - xs[q]) / (xs[q + 1] - xs[q]) } else { 0.0 };
        *o = (ys[q] + f * (ys[q + 1] - ys[q])) / adv;
    }
    (out, 1.0 / adv)
}

pub struct SgTables { shape: Vec<Vec<f64>>, sinuosity: Vec<f64>, ymax: Vec<f64> }

impl SgTables {
    pub fn new() -> SgTables {
        let (mut shape, mut sinuosity, mut ymax) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..OMEGA_STEPS {
            let (t, k) = sg_shape(i as f64 * 0.01);
            ymax.push(t.iter().cloned().fold(0.0, f64::max));
            shape.push(t);
            sinuosity.push(k);
        }
        SgTables { shape, sinuosity, ymax }
    }
    fn lerp_idx(v: &[f64], x: f64) -> f64 {
        // v increasing; returns fractional omega index with v(idx) = x (clamped).
        if x <= v[0] { return 0.0; }
        for i in 1..v.len() { if v[i] >= x { return (i - 1) as f64 + (x - v[i - 1]) / (v[i] - v[i - 1]).max(1e-15); } }
        (v.len() - 1) as f64
    }
    /// omega index for a sinuosity.
    pub fn omega_for_sinuosity(&self, k: f64) -> f64 { Self::lerp_idx(&self.sinuosity, k) }
    /// omega index for a maximum offset (unit wavelength).
    pub fn omega_for_ymax(&self, y: f64) -> f64 { Self::lerp_idx(&self.ymax, y) }
    pub fn ymax(&self, oi: f64) -> f64 { let i = (oi.floor() as usize).min(OMEGA_STEPS - 2); let f = oi - i as f64; self.ymax[i] + f * (self.ymax[i + 1] - self.ymax[i]) }
    pub fn sinuosity(&self, oi: f64) -> f64 { let i = (oi.floor() as usize).min(OMEGA_STEPS - 2); let f = oi - i as f64; self.sinuosity[i] + f * (self.sinuosity[i + 1] - self.sinuosity[i]) }
    /// yhat at omega index oi and phase (cycles).
    pub fn y(&self, oi: f64, phase: f64) -> f64 {
        let i = (oi.floor() as usize).min(OMEGA_STEPS - 2);
        let f = oi - i as f64;
        let ph = phase - phase.floor();
        let xj = ph * X_STEPS as f64;
        let j = (xj.floor() as usize).min(X_STEPS - 1);
        let g = xj - j as f64;
        let at = |t: &Vec<f64>| t[j] + g * (t[(j + 1) % X_STEPS] - t[j]);
        at(&self.shape[i]) * (1.0 - f) + at(&self.shape[i + 1]) * f
    }
}

impl Default for SgTables { fn default() -> Self { Self::new() } }

fn smooth01(x: f64) -> f64 { let t = x.clamp(0.0, 1.0); t * t * (3.0 - 2.0 * t) }

/// A T1a node on a segment: sample index, area of its link (m2, T1a area), keyed meander uniform.
pub struct SegNode { pub k: usize, pub area_m2: f64, pub u: f64 }

/// Meander belt statistics of one segment.
#[derive(Clone, Debug, Default, Serialize)]
pub struct SegMeander { pub segment: u32, pub nodes: u32, pub belt_max_m: f64, pub capped_samples: u64, pub confined_samples: u64, pub resolved_samples: u64, pub samples: u64 }

/// Fill `tab.y` (the meander offsets) from the segment's nodes (sorted by k).
/// `hyd(A) -> (w, h)`: bankfull width and floor half-width (m) for a drainage area; `confine(p)`: the largest floor
/// half-width available at a curve point from other constructions (a volcano edifice corridor cap; INFINITY = none).
/// Wavelength and room are smoothed along the curve (moving average / running minimum over one half wavelength each
/// side) so the offset stays continuous where the area steps at a junction or the room changes.
pub fn meander(tab: &mut SegTab, nodes: &[SegNode], cp: &ClParamsV1, sg: &SgTables, hyd: &dyn Fn(f64) -> (f64, f64), confine: &dyn Fn([f64; 2]) -> f64, seg_id: u32) -> SegMeander {
    let nk = tab.x.len();
    let mut st = SegMeander { segment: seg_id, nodes: nodes.len() as u32, samples: nk as u64, ..Default::default() };
    tab.y = vec![0.0; nk];
    if nodes.is_empty() || nk < 2 { return st; }
    let (kmin, kmax) = (cp.sinuosity_min_permille as f64 / 1000.0, cp.sinuosity_max_permille as f64 / 1000.0);
    let cl = cp.wavelength_widths_permille as f64 / 1000.0;
    let fill = cp.belt_fill_permille as f64 / 1000.0;
    let s_end = tab.s[nk - 1];
    // pass 1: raw wavelength, room and keyed sinuosity per sample
    let (mut lam, mut room, mut ks) = (vec![0.0f64; nk], vec![0.0f64; nk], vec![0.0f64; nk]);
    let mut ni = 0usize;
    for k in 0..nk {
        while ni + 1 < nodes.len() && nodes[ni + 1].k <= k { ni += 1; }
        let (w, h) = hyd(nodes[ni].area_m2);
        lam[k] = (cl * w).max(1e-6);
        let hc = h.min(confine(tab.x[k]));
        if hc < h { st.confined_samples += 1; }
        room[k] = (fill * (hc - 0.5 * w)).max(0.0);
        ks[k] = if k <= nodes[0].k { nodes[0].u } else if ni + 1 < nodes.len() {
            let (k0, k1) = (nodes[ni].k, nodes[ni + 1].k);
            let f = if k1 > k0 { (k - k0) as f64 / (k1 - k0) as f64 } else { 0.0 };
            nodes[ni].u + f * (nodes[ni + 1].u - nodes[ni].u)
        } else { nodes[ni].u };
    }
    // pass 2: smoothing over +-W samples (half the mean wavelength)
    let lam_mean = lam.iter().sum::<f64>() / nk as f64;
    let wsm = ((0.5 * lam_mean / SAMPLE_M).round() as usize).max(1);
    let avg = |v: &[f64]| -> Vec<f64> {
        let mut c = vec![0.0f64; v.len() + 1];
        for k in 0..v.len() { c[k + 1] = c[k] + v[k]; }
        (0..v.len()).map(|k| { let (a, b) = (k.saturating_sub(wsm), (k + wsm).min(v.len() - 1)); (c[b + 1] - c[a]) / (b + 1 - a) as f64 }).collect()
    };
    let lam = avg(&lam);
    let mut rmin = vec![0.0f64; nk];
    {
        // running minimum over +-wsm (monotone deque)
        let mut dq: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
        let mut next = 0usize;
        for k in 0..nk {
            let hi = (k + wsm).min(nk - 1);
            while next <= hi { while dq.back().map(|&b| room[b] >= room[next]).unwrap_or(false) { dq.pop_back(); } dq.push_back(next); next += 1; }
            while dq.front().map(|&f| f + wsm < k).unwrap_or(false) { dq.pop_front(); }
            rmin[k] = room[*dq.front().unwrap()];
        }
    }
    let room = avg(&rmin);
    // pass 3: phase, capped omega, offset
    let mut phase = 0.0f64;
    for k in 0..nk {
        let lk = lam[k];
        if k > 0 { phase += (tab.s[k] - tab.s[k - 1]) / lk; }
        let mut oi = sg.omega_for_sinuosity(kmin + (kmax - kmin) * ks[k]);
        if lk * sg.ymax(oi) > room[k] { oi = sg.omega_for_ymax(room[k] / lk).min(oi); st.capped_samples += 1; }
        let res = smooth01((lk - 256.0) / 256.0);
        if res > 0.0 { st.resolved_samples += 1; }
        let taper = smooth01(tab.s[k] / lk) * smooth01((s_end - tab.s[k]) / lk);
        let y = taper * res * lk * sg.y(oi, phase);
        st.belt_max_m = st.belt_max_m.max(y.abs());
        tab.y[k] = y;
    }
    st
}

/// Window statistics of the centreline (published under prep.net.cl).
#[derive(Clone, Debug, Default, Serialize)]
pub struct ClStats {
    pub t1a_links: u32,
    pub links_on_curve: u32,
    pub junction_links: u32,
    pub nonmonotone_links_straight: u32,
    pub links_without_segment: u32,
    pub connector_max_m: f64,
    pub node_shift_max_m: f64,
    pub node_shift_mean_m: f64,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub meander: Vec<SegMeander>,
    /// Lane MEANDER-Q (cl_mq=1): per-reach planform statistics and the oxbow features.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mq: Option<crate::t1b_mq::MqStats>,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sine_generated_sinuosity_is_one_over_j0() {
        let sg = SgTables::new();
        // J0(1.0) = 0.7652, J0(1.5) = 0.5118 (tables)
        assert!((sg.sinuosity(100.0) - 1.0 / 0.765198).abs() < 1e-3, "{}", sg.sinuosity(100.0));
        assert!((sg.sinuosity(150.0) - 1.0 / 0.511828).abs() < 2e-3, "{}", sg.sinuosity(150.0));
        // zero offset at half wavelengths, symmetric lobes
        assert!(sg.y(150.0, 0.0).abs() < 1e-3 && sg.y(150.0, 0.5).abs() < 1e-3);
        assert!((sg.y(150.0, 0.25) + sg.y(150.0, 0.75)).abs() < 1e-3 && sg.y(150.0, 0.25) > 0.1);
        let oi = sg.omega_for_sinuosity(1.5);
        assert!((sg.sinuosity(oi) - 1.5).abs() < 1e-3);
    }
    #[test]
    fn curve_matches_skeleton_formula_and_projection_is_nearest() {
        let c = T1aSegCurve { id: 0, start: [0.0, 0.0], end: [16384.0, 0.0], amp: 1536.0, ctrl: vec![0.0, 0.5, -0.5, 0.25, 0.0], address: String::new(), seed: None };
        assert_eq!(c.at(0.0), [0.0, 0.0]);
        let p = c.at(0.25);
        assert!((p[0] - 4096.0).abs() < 1e-9 && (p[1] - 768.0).abs() < 1e-9);
        let t = SegTab::new(&c);
        let k = t.project([4100.0, 900.0]);
        assert!((t.x[k][0] - 4100.0).abs() < 60.0);
    }
}
