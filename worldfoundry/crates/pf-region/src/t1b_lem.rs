//! The T1b relaxation solver on a square lattice (lane TR). Same law and
//! schemes as `pf_surface::erosion` (stream-power-n1-linear-diffusion-v0):
//!
//! - STREAM POWER, n = 1, IMPLICIT receiver-first (Braun and Willett 2013,
//!   Geomorphology 180:170): z_i' = (z_i + f b) / (1 + f), f = K A_i^m dt / L_ir,
//!   b = the receiver's new elevation, or its fill level when the receiver is
//!   under a fill. Unconditionally stable; never lowers a post below its
//!   receiver. Posts under a fill (z < fill) get no fluvial change.
//!   Detachment-limited: detached volume leaves along the flow path to the
//!   fixed post the path ends on (exported through an outlet).
//! - HILLSLOPE DIFFUSION dz/dt = D lap z, finite-volume 4-neighbour flux form,
//!   explicit, substeps = ceil(4 D dt / dx^2 / limit) with the declared limit
//!   0.5 (pf-surface's), so each substep is a convex combination. Fixed posts
//!   are Dirichlet; flux into them is exported.
//! - Uplift (the relaxation source term) is added to free posts first.
//!
//! pf-surface's functions are bound to the spherical `GridSpec`, so this is a
//! re-implementation on the planar lod-4 lattice; known-answer tests pin it
//! (1-D steady stream-power profile, steady diffusion parabola).
//!
//! Determinism: one thread; a fixed neighbour order; the priority queue key
//! is (fill, post index) with `f64::total_cmp`.

use crate::t1b::{CHANNEL, FREE, GUIDE, LAKE, MARINE, PLAIN, RING};
use crate::*;
use serde::Serialize;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Clone, Copy, Debug)]
pub struct LemConfig {
    pub k_per_my: f64,
    pub m: f64,
    pub d_m2_per_my: f64,
    pub stability_limit: f64,
    pub dt_my: f64,
    pub steps: usize,
    pub stats_every_steps: usize,
    pub dx_m: f64,
}

#[derive(Clone, Copy, PartialEq)]
struct Key(f64, u32);
impl Eq for Key {}
impl PartialOrd for Key { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for Key {
    // Reversed: BinaryHeap is a max-heap; we pop the lowest (fill, index).
    fn cmp(&self, o: &Self) -> Ordering { o.0.total_cmp(&self.0).then(o.1.cmp(&self.1)) }
}

/// Square lattice n x n, D8 neighbour offsets in a fixed order.
pub const D8N: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];

pub struct Routing {
    /// Posts in priority-flood pop order (receivers before donors).
    pub order: Vec<u32>,
    /// Receiver of each free post (NONE for fixed posts).
    pub receiver: Vec<u32>,
    /// Distance to the receiver (m).
    pub dist: Vec<f64>,
    pub fill: Vec<f64>,
    /// Posts flooded by the unconstrained second pass (0 without zones).
    pub unconfined: u32,
}

/// Priority flood from all fixed posts; steepest descent on the fill surface.
pub fn route(n: usize, z: &[f64], fixed: &[u8], dx: f64) -> Routing { route_zoned(n, z, fixed, dx, None) }

/// Zone confinement for `route_zoned` (lane TW, t0_cell_confined_routing_v1): `zone[i]` is the
/// T0 cell of post i and `drain[i]` the T0 receiver cell of that cell (NONE if none). Flow from j
/// to i is allowed iff zone[i] == zone[j], or i is fixed and zone[i] == drain[j].
#[derive(Clone, Copy)]
pub struct Zones<'a> { pub zone: &'a [u32], pub drain: &'a [u32] }

/// `route` with optional zone confinement. Posts the confined flood cannot reach are flooded in a
/// second, unconstrained pass (their count is `Routing::unconfined`, 0 without zones). With
/// `zones = None` this is exactly the 0.4.0 routing.
pub fn route_zoned(n: usize, z: &[f64], fixed: &[u8], dx: f64, zones: Option<Zones>) -> Routing {
    let nn = n * n;
    let allowed = |j: usize, i: usize| -> bool {
        match zones { None => true, Some(zz) => zz.zone[i] == zz.zone[j] || (fixed[i] != FREE && zz.zone[i] == zz.drain[j]) }
    };
    let mut fill = vec![0f64; nn];
    let mut parent = vec![NONE; nn];
    let mut seen = vec![false; nn];
    let mut heap = BinaryHeap::with_capacity(nn / 4);
    for i in 0..nn {
        if fixed[i] != FREE { seen[i] = true; fill[i] = z[i]; heap.push(Key(z[i], i as u32)); }
    }
    let mut order = Vec::with_capacity(nn);
    let mut pos = vec![u32::MAX; nn];
    let mut phase2 = vec![false; nn];
    let mut unconfined = 0u32;
    for phase in 0..2 {
        if phase == 1 {
            if order.len() == nn || zones.is_none() { break; }
            // Re-seed with everything flooded so far; the rest is flooded unconstrained.
            for &i in &order { heap.push(Key(fill[i as usize], i)); }
        }
        while let Some(Key(f, i)) = heap.pop() {
            let i = i as usize;
            if pos[i] == u32::MAX { pos[i] = order.len() as u32; order.push(i as u32); }
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            for (ox, oy) in D8N {
                let (xx, yy) = (x + ox, y + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let j = yy as usize * n + xx as usize;
                if seen[j] { continue; }
                if phase == 0 && !allowed(j, i) { continue; }
                seen[j] = true;
                if phase == 1 { phase2[j] = true; unconfined += 1; }
                fill[j] = if z[j] > f { z[j] } else { f };
                parent[j] = i as u32;
                heap.push(Key(fill[j], j as u32));
            }
        }
    }
    let diag = dx * std::f64::consts::SQRT_2;
    let mut receiver = vec![NONE; nn];
    let mut dist = vec![0f64; nn];
    for i in 0..nn {
        if fixed[i] != FREE { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let mut best: Option<(f64, usize, f64)> = None;
        for (ox, oy) in D8N {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if fill[j] < fill[i] && pos[j] < pos[i] && (phase2[i] || allowed(i, j)) {
                let d = if ox != 0 && oy != 0 { diag } else { dx };
                let s = (fill[i] - fill[j]) / d;
                let better = match best { None => true, Some((bs, bj, _)) => s > bs || (s == bs && j < bj) };
                if better { best = Some((s, j, d)); }
            }
        }
        match best {
            Some((_, j, d)) => { receiver[i] = j as u32; dist[i] = d; }
            None => {
                let p = parent[i] as usize;
                let (px, py) = ((p % n) as i64, (p / n) as i64);
                receiver[i] = p as u32;
                dist[i] = if px != x && py != y { diag } else { dx };
            }
        }
    }
    Routing { order, receiver, dist, fill, unconfined }
}

/// Drainage area (m2) by accumulation in reverse pop order, plus any
/// boundary inflow area per post.
pub fn accumulate(r: &Routing, cell_area: f64, inflow: Option<&[f64]>) -> Vec<f64> {
    let mut acc = vec![cell_area; r.receiver.len()];
    if let Some(q) = inflow { for (a, b) in acc.iter_mut().zip(q) { *a += *b; } }
    for &i in r.order.iter().rev() {
        let i = i as usize;
        let rc = r.receiver[i];
        if rc != NONE { acc[rc as usize] += acc[i]; }
    }
    acc
}

/// The fixed post each post's flow path ends on.
pub fn terminals(r: &Routing) -> Vec<u32> {
    let mut t = vec![NONE; r.receiver.len()];
    for &i in &r.order {
        let i = i as usize;
        t[i] = if r.receiver[i] == NONE { i as u32 } else { t[r.receiver[i] as usize] };
    }
    t
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Ledger {
    pub uplift_m3: f64,
    pub stream_power_eroded_m3: f64,
    /// Stream-power export by the class of the fixed post the path ends on
    /// [free(unused), channel, marine, ring, lake].
    pub stream_power_export_by_class_m3: [f64; 5],
    pub diffusion_eroded_m3: f64,
    pub diffusion_deposited_m3: f64,
    pub diffusion_export_m3: f64,
    pub diffusion_export_by_class_m3: [f64; 5],
    /// Stream-power and diffusion export into chart-wide guide posts (class 5; absent without guides,
    /// so 0.4.0 ledgers serialize exactly as before).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stream_power_export_to_guides_m3: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub diffusion_export_to_guides_m3: Option<f64>,
    /// Volume added by a declared floor (t1a_trunk_floor_v1; guide layer only; absent otherwise).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub floor_fill_m3: Option<f64>,
    /// Stream-power and creep export into the fixed reclaimed coastal plain (class 6, coast_mode 2; absent otherwise).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub export_to_plain_m3: Option<f64>,
    /// Volume removed by threshold-slope failure and exported (lane T3; absent without lithology).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub failure_eroded_m3: Option<f64>,
    /// Nonlinear creep (0.7.0): edge-steps whose |grad z|/Sc hit the clamp, and edge-steps at or above Sc (absent otherwise).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub creep_clamped_edge_steps: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub creep_over_sc_edge_steps: Option<u64>,
    /// Nonlinear creep solver: most Newton iterations on one line, and lines stopped at the iteration cap.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub creep_newton_iters_max: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub creep_lines_not_converged: Option<u64>,
    pub eroded_m3: f64,
    pub deposited_m3: f64,
    pub exported_m3: f64,
    pub free_volume_start_m3: f64,
    pub free_volume_end_m3: f64,
    /// |eroded - deposited - exported| / eroded.
    pub f7_relative_residual: f64,
    /// |dV - (uplift - exported)| / (uplift + exported).
    pub volume_relative_residual: f64,
    /// |diffusion export by edge flux - diffusion export by net change| / diffusion export.
    pub diffusion_flux_relative_residual: f64,
    pub diffusion_substeps_per_step: u32,
    pub diffusion_number_per_substep: f64,
    pub submerged_post_steps: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct SeriesPoint {
    /// t1b-relaxation-not-history SimTime (My).
    pub t_relax_my: f64,
    pub core_free_mean_m: f64,
    pub core_free_max_m: f64,
    /// Mean |dz/dt| over core free posts since the previous point (m/My).
    pub core_mean_abs_rate_m_per_my: f64,
    pub core_max_abs_rate_m_per_my: f64,
    pub core_pits: u32,
}

pub struct LemResult {
    pub z: Vec<f64>,
    pub ledger: Ledger,
    pub series: Vec<SeriesPoint>,
}

/// Run the relaxation. `core` marks posts for the convergence series.
pub fn run(n: usize, fixed: &[u8], z0: &[f64], uplift_m_per_my: &[f64], inflow: &[f64], core: &[bool], cfg: &LemConfig) -> LemResult {
    run_zoned(n, fixed, z0, uplift_m_per_my, inflow, core, cfg, None, None)
}

/// `run` with zone-confined routing (`route_zoned`); `None` is exactly `run`.
#[allow(clippy::too_many_arguments)]
pub fn run_zoned(n: usize, fixed: &[u8], z0: &[f64], uplift_m_per_my: &[f64], inflow: &[f64], core: &[bool], cfg: &LemConfig, zones: Option<Zones>, floor: Option<&[f64]>) -> LemResult {
    run_litho(n, fixed, z0, uplift_m_per_my, inflow, core, cfg, zones, floor, None)
}

/// Threshold-slope failure (lane T3, threshold-slope-lithology-v0): lower every free post so that
/// z_i <= z_j + sc_i d_ij for all 8 neighbours j. Exact min-plus closure: a Dijkstra pass seeded with
/// the lower neighbours of initially violating posts, key (z, index), so the result is the unique
/// closure and deterministic. Returns the removed height summed over posts (m).
pub fn threshold_cap(n: usize, z: &mut [f64], fixed: &[u8], sc: &[f64], dx: f64) -> f64 {
    const EPS: f64 = 1e-9;
    let nn = n * n;
    let diag = dx * std::f64::consts::SQRT_2;
    let mut heap = BinaryHeap::new();
    for i in 0..nn {
        if fixed[i] != FREE { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8N {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            let d = if ox != 0 && oy != 0 { diag } else { dx };
            if z[i] > z[j] + sc[i] * d + EPS { heap.push(Key(z[j], j as u32)); }
        }
    }
    let mut removed = 0.0;
    while let Some(Key(v, j)) = heap.pop() {
        let j = j as usize;
        if v != z[j] { continue; } // stale entry: j was lowered after this push
        let (x, y) = ((j % n) as i64, (j / n) as i64);
        for (ox, oy) in D8N {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let i = yy as usize * n + xx as usize;
            if fixed[i] != FREE { continue; }
            let d = if ox != 0 && oy != 0 { diag } else { dx };
            let cap = z[j] + sc[i] * d;
            if z[i] > cap + EPS { removed += z[i] - cap; z[i] = cap; heap.push(Key(cap, i as u32)); }
        }
    }
    removed
}

/// `run_zoned` with exposed lithology (lane T3, T1b 0.6.0 / guides 0.4.0): per-post K multiplier and
/// threshold-slope failure from `litho` (evaluated on the surface after the source term, once per step).
/// `None` is exactly `run_zoned`.
#[allow(clippy::too_many_arguments)]
pub fn run_litho(n: usize, fixed: &[u8], z0: &[f64], uplift_m_per_my: &[f64], inflow: &[f64], core: &[bool], cfg: &LemConfig, zones: Option<Zones>, floor: Option<&[f64]>,
    litho: Option<&dyn crate::t1b_litho::ExposedLithology>) -> LemResult {
    run_full(n, fixed, z0, uplift_m_per_my, inflow, core, cfg, zones, floor, litho, None)
}

/// Nonlinear (Roering-type) hillslope creep for T1b 0.7.0 (lane T3b, TER-REGION v4).
#[derive(Clone, Copy, Debug)]
pub struct RoeringConfig {
    /// Above S = ratio_max x Sc the transport factor 1 / (1 - (S/Sc)^2) continues along its tangent (finite, still
    /// increasing, so the flux stays monotone in S).
    pub ratio_max: f64,
}

/// Scratch for `roering_diffuse` (reused across steps).
pub struct RoeringScratch { crx: Vec<f64>, cry: Vec<f64>, scx: Vec<f64>, scy: Vec<f64>, a: Vec<f64>, b: Vec<f64>, c: Vec<f64>, d: Vec<f64>, idx: Vec<usize>, zl: Vec<f64>, z0: Vec<f64>, zb: Vec<f64> }
impl RoeringScratch {
    pub fn new(n: usize) -> Self {
        let v = |k: usize| vec![0.0; k];
        RoeringScratch { crx: v(n * n), cry: v(n * n), scx: v(n * n), scy: v(n * n), a: v(n), b: v(n), c: v(n), d: v(n), idx: Vec::with_capacity(n), zl: v(n + 2), z0: v(n + 2), zb: v(n + 2) }
    }
}

/// Per-step volume bookkeeping of `roering_diffuse` (m3) and solver diagnostics.
#[derive(Default, Clone, Copy, Debug)]
pub struct RoeringStepLedger { pub eroded: f64, pub deposited: f64, pub net_export: f64, pub flux_export: f64, pub export_by_class: [f64; 8], pub clamped_edges: u64, pub over_sc_edges: u64, pub newton_iters_max: u32, pub lines_not_converged: u32 }

/// Transport factor f(S) and its derivative (see `RoeringConfig`).
#[inline]
fn creep_factor(s: f64, sc: f64, rm: f64) -> (f64, f64) {
    let r = s / sc;
    if r < rm {
        let q = 1.0 - r * r;
        (1.0 / q, 2.0 * r / sc / (q * q))
    } else {
        let q = 1.0 - rm * rm;
        let (fm, dfm) = (1.0 / q, 2.0 * rm / sc / (q * q));
        (fm + dfm * (s - rm * sc), dfm)
    }
}

/// One step of nonlinear creep, q = D S / (1 - (S/Sc)^2) (Roering et al. 1999 WRR 35:853 form), on free posts.
/// Finite volume on the 4-neighbour edges; S = |grad z| at the edge (axis difference, plus the cross-axis central
/// difference averaged over the edge's two posts, lagged from the start of the sweep), so the limit is isotropic; Sc of
/// an edge is the critical slope of its upper free post at the start of the sweep. Time: locally one-dimensional
/// backward Euler (x sweep, then y sweep), each sweep FULLY implicit in the along-line gradient: Newton iterations with
/// a tridiagonal Jacobian per run of free posts and a backtracking line search on max |R|, to |dz| <= 1e-9 m or
/// 200 iterations (counted). Fixed posts are Dirichlet; flux into them is exported. The lattice edge is closed.
pub fn roering_diffuse(n: usize, z: &mut [f64], fixed: &[u8], sc: &[f64], d: f64, dt: f64, dx: f64, rc: &RoeringConfig, s: &mut RoeringScratch, area: f64) -> RoeringStepLedger {
    let mut led = RoeringStepLedger::default();
    let rm = rc.ratio_max;
    let k_dt = dt / (dx * dx);
    for sweep in 0..2usize {
        // Lagged cross gradients and edge Sc for this sweep's edges (edge i = (i, i + stride)).
        let stride = if sweep == 0 { 1 } else { n };
        for y in 0..n { for x in 0..n {
            let i = y * n + x;
            let has = if sweep == 0 { x + 1 < n } else { y + 1 < n };
            if !has { continue; }
            let j = i + stride;
            let (fi, fj) = (fixed[i] == FREE, fixed[j] == FREE);
            if !fi && !fj { continue; }
            let up = if fi && fj { if z[i] >= z[j] { i } else { j } } else if fi { i } else { j };
            let cross_at = |k: usize| -> f64 {
                let (kx, ky) = (k % n, k / n);
                if sweep == 0 {
                    let (u, dd, h) = if ky + 1 < n && ky > 0 { (k + n, k - n, 2.0) } else if ky + 1 < n { (k + n, k, 1.0) } else { (k, k - n, 1.0) };
                    (z[u] - z[dd]) / (h * dx)
                } else {
                    let (u, dd, h) = if kx + 1 < n && kx > 0 { (k + 1, k - 1, 2.0) } else if kx + 1 < n { (k + 1, k, 1.0) } else { (k, k - 1, 1.0) };
                    (z[u] - z[dd]) / (h * dx)
                }
            };
            let cr = 0.5 * (cross_at(i) + cross_at(j));
            if sweep == 0 { s.crx[i] = cr; s.scx[i] = sc[up]; } else { s.cry[i] = cr; s.scy[i] = sc[up]; }
        } }
        for line in 0..n {
            let start = if sweep == 0 { line * n } else { line };
            let mut k = 0usize;
            while k < n {
                if fixed[start + k * stride] != FREE { k += 1; continue; }
                let mut k2 = k;
                s.idx.clear();
                while k2 < n && fixed[start + k2 * stride] == FREE { s.idx.push(start + k2 * stride); k2 += 1; }
                let m = s.idx.len();
                let left = (k > 0).then(|| start + (k - 1) * stride);
                let right = (k2 < n).then(|| start + k2 * stride);
                // Local line: zl[0] = left fixed (if any), zl[1..=m] the run, zl[m+1] = right fixed (if any).
                for t in 0..m { s.zl[t + 1] = z[s.idx[t]]; s.z0[t + 1] = z[s.idx[t]]; }
                if let Some(j) = left { s.zl[0] = z[j]; }
                if let Some(j) = right { s.zl[m + 1] = z[j]; }
                // Edge e (0..=m) joins zl[e] and zl[e+1]; edge 0 exists iff left, edge m iff right. Its lagged data
                // lives at the lower-index post of the pair in lattice order.
                let edge_at = |s: &RoeringScratch, e: usize| -> Option<(f64, f64)> {
                    if (e == 0 && left.is_none()) || (e == m && right.is_none()) { return None; }
                    let lo = if e == 0 { left.unwrap() } else { s.idx[e - 1] };
                    Some(if sweep == 0 { (s.crx[lo], s.scx[lo]) } else { (s.cry[lo], s.scy[lo]) })
                };
                // flux on edge e toward zl[e] (mass from e+1 to e when positive) and dF/dG.
                let flux = |s: &RoeringScratch, e: usize| -> Option<(f64, f64, f64)> {
                    let (cr, scv) = edge_at(s, e)?;
                    let g = (s.zl[e + 1] - s.zl[e]) / dx;
                    let sm = (g * g + cr * cr).sqrt();
                    let (f, df) = creep_factor(sm, scv, rm);
                    let dfdg = if sm > 0.0 { f + df * g * g / sm } else { f };
                    Some((d * f * g, d * dfdg, sm / scv))
                };
                // max |R| on the current line state.
                let resid_max = |s: &RoeringScratch| -> f64 {
                    let mut fl = flux(s, 0).map_or(0.0, |v| v.0);
                    let mut mx = 0.0f64;
                    for t in 0..m {
                        let fr = flux(s, t + 1).map_or(0.0, |v| v.0);
                        mx = mx.max((s.zl[t + 1] - s.z0[t + 1] - dt / dx * (fr - fl)).abs());
                        fl = fr;
                    }
                    mx
                };
                let mut iters = 0u32;
                loop {
                    iters += 1;
                    let mut r0 = 0.0f64;
                    let mut fl = flux(s, 0);
                    for t in 0..m {
                        let fr = flux(s, t + 1);
                        let (fl_v, gl) = fl.map_or((0.0, 0.0), |v| (v.0, v.1));
                        let (fr_v, gr) = fr.map_or((0.0, 0.0), |v| (v.0, v.1));
                        // R = z - z0 - dt/dx (F_r - F_l)
                        s.d[t] = -(s.zl[t + 1] - s.z0[t + 1] - dt / dx * (fr_v - fl_v));
                        r0 = r0.max(s.d[t].abs());
                        s.b[t] = 1.0 + k_dt * (gl + gr);
                        s.a[t] = if t > 0 { -k_dt * gl } else { 0.0 };
                        s.c[t] = if t + 1 < m { -k_dt * gr } else { 0.0 };
                        fl = fr;
                    }
                    for t in 1..m { let w = s.a[t] / s.b[t - 1]; s.b[t] -= w * s.c[t - 1]; s.d[t] -= w * s.d[t - 1]; }
                    let mut nxt = 0.0;
                    let mut big = 0.0f64;
                    for t in (0..m).rev() {
                        let v = if t + 1 == m { s.d[t] / s.b[t] } else { (s.d[t] - s.c[t] * nxt) / s.b[t] };
                        nxt = v; s.d[t] = v; big = big.max(v.abs());
                    }
                    // Backtracking line search on max |R| (halve the step until the residual falls; at most 30 halvings).
                    for t in 0..m { s.zb[t + 1] = s.zl[t + 1]; }
                    let mut lam = 1.0f64;
                    loop {
                        for t in 0..m { s.zl[t + 1] = s.zb[t + 1] + lam * s.d[t]; }
                        if resid_max(s) < r0 || lam < 1e-9 { break; }
                        lam *= 0.5;
                    }
                    if big * lam <= 1e-9 || r0 <= 1e-10 { break; }
                    if iters >= 200 { led.lines_not_converged += 1; break; }
                }
                led.newton_iters_max = led.newton_iters_max.max(iters);
                for t in 0..m {
                    let i = s.idx[t];
                    let dz = s.zl[t + 1] - z[i];
                    if dz < 0.0 { led.eroded -= dz * area; } else { led.deposited += dz * area; }
                    led.net_export -= dz * area;
                    z[i] = s.zl[t + 1];
                }
                // Diagnostics and export into the fixed ends (nonlinear flux on the final line state).
                for e in 0..=m {
                    if let Some((f, _, r)) = flux(s, e) {
                        if r >= 1.0 { led.over_sc_edges += 1; }
                        if r >= rm { led.clamped_edges += 1; }
                        // f > 0 moves mass from zl[e+1] to zl[e]: out of the run at the left end, into it at the right end.
                        if e == 0 { let j = left.unwrap(); let q = f * dt / dx * area; led.flux_export += q; led.export_by_class[fixed[j] as usize] += q; }
                        if e == m { let j = right.unwrap(); let q = -f * dt / dx * area; led.flux_export += q; led.export_by_class[fixed[j] as usize] += q; }
                    }
                }
                k = k2;
            }
        }
    }
    led
}

/// `run_litho` with an optional nonlinear-creep law (T1b 0.7.0): `roering = Some(..)` replaces the explicit linear
/// diffusion and the threshold cap with `roering_diffuse` (needs `litho` for Sc). `None` is exactly `run_litho`.
#[allow(clippy::too_many_arguments)]
pub fn run_full(n: usize, fixed: &[u8], z0: &[f64], uplift_m_per_my: &[f64], inflow: &[f64], core: &[bool], cfg: &LemConfig, zones: Option<Zones>, floor: Option<&[f64]>,
    litho: Option<&dyn crate::t1b_litho::ExposedLithology>, roering: Option<RoeringConfig>) -> LemResult {
    assert!(roering.is_none() || litho.is_some(), "nonlinear creep needs exposed lithology (Sc)");
    let nn = n * n;
    let mut sc = if litho.is_some() { vec![0f64; nn] } else { Vec::new() };
    let mut kmul = if litho.is_some() { vec![1f64; nn] } else { Vec::new() };
    let a = cfg.dx_m * cfg.dx_m;
    let mut z = z0.to_vec();
    let mut led = Ledger::default();
    let free_vol = |z: &[f64]| -> f64 { (0..nn).filter(|&i| fixed[i] == FREE).map(|i| z[i] * a).sum() };
    led.free_volume_start_m3 = free_vol(&z);
    let dt = cfg.dt_my;
    let nu_step = 4.0 * cfg.d_m2_per_my * dt / (cfg.dx_m * cfg.dx_m);
    let nsub = ((nu_step / cfg.stability_limit).ceil() as u32).max(1);
    let c = cfg.d_m2_per_my * (dt / nsub as f64) / (cfg.dx_m * cfg.dx_m);
    led.diffusion_substeps_per_step = nsub;
    led.diffusion_number_per_substep = 4.0 * c;
    if roering.is_some() {
        // Implicit (no substeps); report the base diffusion number of one step.
        led.diffusion_substeps_per_step = 0;
        led.diffusion_number_per_substep = 4.0 * cfg.d_m2_per_my * dt / (cfg.dx_m * cfg.dx_m);
    }
    let mut series = Vec::new();
    let mut z_prev = z.clone();
    let mut dz = vec![0f64; nn];
    let mut wet = vec![false; nn];
    let mut diff_flux_export = 0.0;
    let mut last_done = 0usize;
    let mut rscratch = roering.map(|_| RoeringScratch::new(n));
    let core_stats = |z: &[f64], z_prev: &[f64], dt_span: f64, pits: u32, t: f64| -> SeriesPoint {
        let (mut s, mut mx, mut cnt, mut rs, mut rmx) = (0.0, f64::MIN, 0usize, 0.0, 0.0f64);
        for i in 0..nn {
            if !core[i] || fixed[i] != FREE { continue; }
            s += z[i]; mx = mx.max(z[i]); cnt += 1;
            if dt_span > 0.0 { let r = ((z[i] - z_prev[i]) / dt_span).abs(); rs += r; rmx = rmx.max(r); }
        }
        SeriesPoint { t_relax_my: t, core_free_mean_m: s / cnt.max(1) as f64, core_free_max_m: mx,
            core_mean_abs_rate_m_per_my: rs / cnt.max(1) as f64, core_max_abs_rate_m_per_my: rmx, core_pits: pits }
    };
    series.push(core_stats(&z, &z_prev, 0.0, 0, 0.0));
    for step in 0..cfg.steps {
        // 1. Relaxation source term.
        for i in 0..nn {
            if fixed[i] == FREE { let u = uplift_m_per_my[i] * dt; z[i] += u; led.uplift_m3 += u * a; }
        }
        // 1b. Exposed lithology on the surface after the source term.
        if let Some(l) = litho {
            for i in 0..nn { if fixed[i] == FREE { let (s, k, _) = l.at(i, z[i]); sc[i] = s; kmul[i] = k; } }
        }
        // 2. Routing.
        let r = route_zoned(n, &z, fixed, cfg.dx_m, zones);
        let acc = accumulate(&r, a, Some(inflow));
        let term = terminals(&r);
        // 3. Stream power, implicit, receivers first. "Wet" (under a fill) is
        //    decided on the surface BEFORE this pass: a receiver lowered earlier
        //    in the pass is not wet (lane TR bug fix, 2026-10-12, about 20:27 EDT).
        for i in 0..nn { wet[i] = z[i] < r.fill[i]; }
        for &i in &r.order {
            let i = i as usize;
            if fixed[i] != FREE { continue; }
            if wet[i] { led.submerged_post_steps += 1; continue; }
            let rc = r.receiver[i] as usize;
            let b = if wet[rc] { r.fill[rc] } else { z[rc] };
            if z[i] <= b { continue; }
            let kk = if litho.is_some() { cfg.k_per_my * kmul[i] } else { cfg.k_per_my };
            let f = kk * acc[i].powf(cfg.m) * dt / r.dist[i];
            let zn = (z[i] + f * b) / (1.0 + f);
            let v = (z[i] - zn) * a;
            led.stream_power_eroded_m3 += v;
            let cls = fixed[term[i] as usize];
            if cls == GUIDE { *led.stream_power_export_to_guides_m3.get_or_insert(0.0) += v; }
            else if cls == PLAIN { *led.export_to_plain_m3.get_or_insert(0.0) += v; }
            else { led.stream_power_export_by_class_m3[cls as usize] += v; }
            z[i] = zn;
        }
        // 4. Diffusion: nonlinear creep (0.7.0) or explicit linear substeps.
        if let Some(rc) = roering.as_ref() {
            let st = roering_diffuse(n, &mut z, fixed, &sc, cfg.d_m2_per_my, dt, cfg.dx_m, rc, rscratch.as_mut().unwrap(), a);
            led.diffusion_eroded_m3 += st.eroded;
            led.diffusion_deposited_m3 += st.deposited;
            led.diffusion_export_m3 += st.net_export;
            diff_flux_export += st.flux_export;
            for (k, v) in st.export_by_class.iter().enumerate() {
                if *v == 0.0 { continue; }
                if k == GUIDE as usize { *led.diffusion_export_to_guides_m3.get_or_insert(0.0) += v; }
                else if k == PLAIN as usize { *led.export_to_plain_m3.get_or_insert(0.0) += v; }
                else { led.diffusion_export_by_class_m3[k] += v; }
            }
            *led.creep_clamped_edge_steps.get_or_insert(0) += st.clamped_edges;
            *led.creep_over_sc_edge_steps.get_or_insert(0) += st.over_sc_edges;
            let it = led.creep_newton_iters_max.get_or_insert(0); *it = (*it).max(st.newton_iters_max);
            *led.creep_lines_not_converged.get_or_insert(0) += st.lines_not_converged as u64;
        }
        for _ in 0..(if roering.is_some() { 0 } else { nsub }) {
            for y in 0..n {
                for x in 0..n {
                    let i = y * n + x;
                    if fixed[i] != FREE { dz[i] = 0.0; continue; }
                    let zi = z[i];
                    // No flux across the lattice edge (closed ring).
                    let nbr = [(x + 1 < n).then(|| i + 1), (x > 0).then(|| i - 1), (y + 1 < n).then(|| i + n), (y > 0).then(|| i - n)];
                    let mut sum = 0.0;
                    for j in nbr.into_iter().flatten() { sum += z[j] - zi; }
                    dz[i] = c * sum;
                    for j in nbr.into_iter().flatten() {
                        let fj = fixed[j];
                        if fj != FREE {
                            let q = c * (zi - z[j]) * a; // volume from i into fixed j
                            if fj == GUIDE { *led.diffusion_export_to_guides_m3.get_or_insert(0.0) += q; }
                            else { led.diffusion_export_by_class_m3[fj as usize] += q; }
                            diff_flux_export += q;
                        }
                    }
                }
            }
            for i in 0..nn {
                if fixed[i] != FREE { continue; }
                let d = dz[i];
                if d < 0.0 { led.diffusion_eroded_m3 -= d * a; } else { led.diffusion_deposited_m3 += d * a; }
                led.diffusion_export_m3 -= d * a;
                z[i] += d;
            }
        }
        // 4b. Threshold-slope failure (exported, like detachment-limited stream power).
        if litho.is_some() && roering.is_none() {
            let v = threshold_cap(n, &mut z, fixed, &sc, cfg.dx_m) * a;
            *led.failure_eroded_m3.get_or_insert(0.0) += v;
        }
        // 5. Declared floor (guide layer only): free posts never below floor[i].
        if let Some(fl) = floor {
            for i in 0..nn {
                if fixed[i] == FREE && z[i] < fl[i] { *led.floor_fill_m3.get_or_insert(0.0) += (fl[i] - z[i]) * a; z[i] = fl[i]; }
            }
        }
        let done = step + 1;
        if done % cfg.stats_every_steps == 0 || done == cfg.steps {
            let span = (done - last_done) as f64 * dt;
            last_done = done;
            let pits = (0..nn).filter(|&i| core[i] && fixed[i] == FREE && z[i] < r.fill[i]).count() as u32;
            series.push(core_stats(&z, &z_prev, span, pits, done as f64 * dt));
            z_prev.copy_from_slice(&z);
        }
    }
    led.free_volume_end_m3 = free_vol(&z);
    let fail = led.failure_eroded_m3.unwrap_or(0.0);
    led.eroded_m3 = led.stream_power_eroded_m3 + led.diffusion_eroded_m3 + fail;
    led.deposited_m3 = led.diffusion_deposited_m3;
    led.exported_m3 = led.stream_power_eroded_m3 + led.diffusion_export_m3 + fail;
    led.f7_relative_residual = if led.eroded_m3 > 0.0 { (led.eroded_m3 - led.deposited_m3 - led.exported_m3).abs() / led.eroded_m3 } else { 0.0 };
    let dv = led.free_volume_end_m3 - led.free_volume_start_m3;
    let ff = led.floor_fill_m3.unwrap_or(0.0);
    let denom = led.uplift_m3 + ff + led.exported_m3.abs();
    led.volume_relative_residual = if denom > 0.0 { (dv - (led.uplift_m3 + ff - led.exported_m3)).abs() / denom } else { 0.0 };
    led.diffusion_flux_relative_residual = if led.diffusion_export_m3.abs() > 0.0 {
        (diff_flux_export - led.diffusion_export_m3).abs() / led.diffusion_export_m3.abs()
    } else { (diff_flux_export - led.diffusion_export_m3).abs() };
    let _ = (CHANNEL, MARINE, RING, LAKE);
    LemResult { z, ledger: led, series }
}
