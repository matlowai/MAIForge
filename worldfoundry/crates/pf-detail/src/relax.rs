//! B5 tile-local gully relaxation (SimTime `t2-relaxation-not-history`).
//!
//! Per step: MFD drainage area (Freeman 1991 exponent with Quinn-type contour
//! weights), implicit detachment-limited stream power on the D8 steepest-descent
//! receiver (Braun and Willett 2013 ordering) where A >= A_c, then explicit
//! 4-neighbour nonlinear creep q = D S / (1 - (S/Sc)^2) (Roering et al. 1999 form)
//! with the factor capped at S/Sc = cap. Lithology (Sc, K multiplier) is re-evaluated
//! at the current surface every step. Fixed posts never change.

use crate::{DetailParamsV1, POST_M};

const D8: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];

pub struct RelaxStats {
    pub steps: usize,
    pub max_change_m: f64,
    pub mean_change_m: f64,
    pub max_area_m2: f64,
}

/// `n` x `n` domain, row-major (index = j * n + i, j north). `extra_area` is added to each post's own area.
/// Lane DE DE3 (`bed_k` Some only with de=1 and de_sp): the stream-power K blends from the litho K (v3 W4 window) to the bed K
/// returned by `bed_k` with drainage area, w = smoothstep(ln(A / a_lo) / ln(a_hi / a_lo)) (geometric interpolation of K).
/// Lane RIB (`runoff` Some only with de_runoff=1): runoff(post, z) = (f, t): each post contributes cell area x f to the drainage
/// area (runoff-weighted area), and passes on t x the area it receives (transmission loss into a permeable substrate along the
/// flow path, t = exp(-2 m / L) per post); the extra area entering at the domain ring is unchanged where it enters.
#[allow(clippy::too_many_arguments)]
pub fn relax(n: usize, z: &mut [f64], fixed: &[bool], extra_area: &[f64], litho: &dyn Fn(usize, f64) -> (f64, f64), bed_k: Option<&dyn Fn(usize, f64) -> f64>,
             runoff: Option<&dyn Fn(usize, f64) -> (f64, f64)>, p: &DetailParamsV1) -> RelaxStats {
    relax_dx(n, POST_M, z, fixed, extra_area, litho, bed_k, runoff, p)
}

/// `relax` on a lattice of spacing `dx` m (lane RILL de_b5_chain: the coarse block-free pass; `relax` is `relax_dx` at POST_M,
/// the same arithmetic).
#[allow(clippy::too_many_arguments)]
pub fn relax_dx(n: usize, dx: f64, z: &mut [f64], fixed: &[bool], extra_area: &[f64], litho: &dyn Fn(usize, f64) -> (f64, f64), bed_k: Option<&dyn Fn(usize, f64) -> f64>,
             runoff: Option<&dyn Fn(usize, f64) -> (f64, f64)>, p: &DetailParamsV1) -> RelaxStats {
    relax_dx_nf(n, dx, z, fixed, None, extra_area, litho, bed_k, runoff, p)
}

/// `relax_dx` with lane CARRY2's no-flux posts (`noflux` Some only with ct_sharp_riser_noflux=1): a 4-neighbour creep edge with a
/// no-flux post at either end carries no flux (held riser posts are bare rock: a zero-flux boundary of the creeping surface instead
/// of a fixed-height source / sink); stream power and drainage area unchanged (no-flux posts are also `fixed`).
#[allow(clippy::too_many_arguments)]
pub fn relax_dx_nf(n: usize, dx: f64, z: &mut [f64], fixed: &[bool], noflux: Option<&[bool]>, extra_area: &[f64], litho: &dyn Fn(usize, f64) -> (f64, f64),
             bed_k: Option<&dyn Fn(usize, f64) -> f64>, runoff: Option<&dyn Fn(usize, f64) -> (f64, f64)>, p: &DetailParamsV1) -> RelaxStats {
    let nn = n * n;
    let z0: Vec<f64> = z.to_vec();
    let cell_a = dx * dx;
    let dist = [dx, dx * std::f64::consts::SQRT_2];
    // Quinn et al. 1991 contour-length weights (0.5 dx orthogonal, 0.354 dx diagonal).
    let clen = [0.5 * dx, 0.354 * dx];
    let mut order: Vec<usize> = (0..nn).collect();
    let mut area = vec![0.0f64; nn];
    let mut rcv = vec![usize::MAX; nn];
    let mut rlen = vec![0.0f64; nn];
    let mut sc = vec![0.0f64; nn];
    let mut kk = vec![0.0f64; nn];
    let mut dz = vec![0.0f64; nn];
    let mut tr = if runoff.is_some() { vec![1.0f64; nn] } else { Vec::new() };
    let mut max_area: f64 = 0.0;
    let head_law = p.v3p().filter(|w| w.head_law_on).map(|_| p.channel_head_area_m2);
    let nb = |i: usize, d: usize| -> Option<usize> {
        let (x, y) = ((i % n) as i64 + D8[d].0, (i / n) as i64 + D8[d].1);
        if x < 0 || y < 0 || x >= n as i64 || y >= n as i64 { None } else { Some(y as usize * n + x as usize) }
    };
    for _step in 0..p.steps {
        for i in 0..nn { let (s, k) = litho(i, z[i]); sc[i] = s; kk[i] = k; }
        // Highest first; ties by index (deterministic).
        order.sort_unstable_by(|&a, &b| z[b].partial_cmp(&z[a]).unwrap().then(a.cmp(&b)));
        match runoff {
            None => for i in 0..nn { area[i] = cell_a + extra_area[i]; },
            Some(r) => for i in 0..nn { let (f, t) = r(i, z[i]); area[i] = cell_a * f + extra_area[i]; tr[i] = t; },
        }
        for &i in &order {
            // steepest-descent receiver
            let mut best = (0.0f64, usize::MAX, 0.0f64);
            let mut wsum = 0.0;
            let mut w = [0.0f64; 8];
            for d in 0..8 {
                if let Some(j) = nb(i, d) {
                    let s = (z[i] - z[j]) / dist[d & 1];
                    if s > 0.0 {
                        if s > best.0 { best = (s, j, dist[d & 1]); }
                        w[d] = clen[d & 1] * s.powf(p.mfd_exponent);
                        wsum += w[d];
                    }
                }
            }
            rcv[i] = best.1;
            rlen[i] = best.2;
            if fixed[i] || wsum <= 0.0 { continue; }
            let a = if runoff.is_some() { area[i] * tr[i] } else { area[i] };
            for d in 0..8 {
                if w[d] > 0.0 { let j = nb(i, d).unwrap(); area[j] += a * w[d] / wsum; }
            }
        }
        // Implicit stream power, lowest first (receivers are lower at step start).
        for &i in order.iter().rev() {
            max_area = max_area.max(area[i]);
            if fixed[i] || rcv[i] == usize::MAX || area[i] < p.channel_head_area_m2 { continue; }
            // v3 W8 (lane T2C revision 1): channel-head area-slope law A_c = A_1 / S^2 (Montgomery and Dietrich 1988),
            // A_1 = the v2 constant head area at S = 1; gentler ground needs a larger area before it incises
            if let Some(w) = head_law {
                let s = (z[i] - z[rcv[i]]) / rlen[i];
                if s <= 0.0 || area[i] * s * s < w { continue; }
            }
            let kf = match (bed_k, p.de.as_ref().and_then(|d| d.sp_area_m2)) {
                (Some(bk), Some([alo, ahi])) if area[i] > alo => {
                    let w = crate::smoothstep((area[i] / alo).ln() / (ahi / alo).ln());
                    (kk[i].ln() * (1.0 - w) + bk(i, z[i]).ln() * w).exp()
                }
                _ => kk[i],
            };
            let f = p.k_sp_per_step * kf * area[i].sqrt() / rlen[i];
            let zr = z[rcv[i]];
            let zn = (z[i] + f * zr) / (1.0 + f);
            if zn < z[i] { z[i] = zn; }
        }
        // Explicit nonlinear creep on 4-neighbour edges.
        for v in dz.iter_mut() { *v = 0.0; }
        let cap = p.creep_ratio_cap;
        let fmax = 1.0 / (1.0 - cap * cap);
        let coef = p.d_creep_m2_per_step / (dx * dx);
        for j in 0..n { for i in 0..n {
            let a = j * n + i;
            for b in [if i + 1 < n { Some(a + 1) } else { None }, if j + 1 < n { Some(a + n) } else { None }].into_iter().flatten() {
                let diff = z[a] - z[b];
                if diff == 0.0 { continue; }
                if let Some(nf) = noflux { if nf[a] || nf[b] { continue; } }
                let up = if diff > 0.0 { a } else { b };
                let r = (diff.abs() / dx) / sc[up];
                let fac = if r >= cap { fmax } else { 1.0 / (1.0 - r * r) };
                let flux = coef * fac * diff;
                dz[a] -= flux;
                dz[b] += flux;
            }
        } }
        for i in 0..nn { if !fixed[i] { z[i] += dz[i]; } }
    }
    let mut mx: f64 = 0.0;
    let mut sum = 0.0;
    for i in 0..nn { let d = (z[i] - z0[i]).abs(); mx = mx.max(d); sum += d; }
    RelaxStats { steps: p.steps, max_change_m: mx, mean_change_m: sum / nn as f64, max_area_m2: max_area }
}

/// Lane RILL (de_b5_chain): MFD drainage area (m2) on a fixed surface, the same routing as one `relax_dx` step (Freeman exponent,
/// Quinn contour weights, fixed posts pass nothing on, runoff (f, t) per post when given). Returns (area, receiver weights are
/// recomputed by the caller where needed).
pub fn mfd_area(n: usize, dx: f64, z: &[f64], fixed: &[bool], runoff: Option<(&[f64], &[f64])>, p: &DetailParamsV1) -> Vec<f64> {
    let nn = n * n;
    let cell_a = dx * dx;
    let dist = [dx, dx * std::f64::consts::SQRT_2];
    let clen = [0.5 * dx, 0.354 * dx];
    let mut order: Vec<usize> = (0..nn).collect();
    order.sort_unstable_by(|&a, &b| z[b].partial_cmp(&z[a]).unwrap().then(a.cmp(&b)));
    let mut area: Vec<f64> = match runoff { None => vec![cell_a; nn], Some((f, _)) => (0..nn).map(|i| cell_a * f[i]).collect() };
    let nb = |i: usize, d: usize| -> Option<usize> {
        let (x, y) = ((i % n) as i64 + D8[d].0, (i / n) as i64 + D8[d].1);
        if x < 0 || y < 0 || x >= n as i64 || y >= n as i64 { None } else { Some(y as usize * n + x as usize) }
    };
    for &i in &order {
        if fixed[i] { continue; }
        let mut w = [0.0f64; 8];
        let mut wsum = 0.0;
        for d in 0..8 {
            if let Some(j) = nb(i, d) {
                let s = (z[i] - z[j]) / dist[d & 1];
                if s > 0.0 { w[d] = clen[d & 1] * s.powf(p.mfd_exponent); wsum += w[d]; }
            }
        }
        if wsum <= 0.0 { continue; }
        let a = match runoff { Some((_, t)) => area[i] * t[i], None => area[i] };
        for d in 0..8 { if w[d] > 0.0 { let j = nb(i, d).unwrap(); area[j] += a * w[d] / wsum; } }
    }
    area
}
