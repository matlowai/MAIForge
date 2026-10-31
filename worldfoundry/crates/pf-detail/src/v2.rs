//! pf-detail v2 (lane T2B, 2026-10-14; generator 0.2.0, switch `version=2`). PREREG:
//! `unpublished-design-record` (stages V1-V7).
//!
//! V1 channel curves traced on the T1b channel posts and smoothed (no 32 m raster trenches);
//! V2 valley reprojection (the valley axis follows the curve); V3 carve with banks graded at the
//! material's limiting slope; V4 T1b over-codes (edifice, deposits) honoured; V5 hard/soft bed
//! pairs conserving the class rock-mass critical slope; V6 slope clamp of the remap; V7 amplitude
//! bound tied to the T1b local relief (ADR 043 condition 7).
//!
//! Everything here is a pure function of the T1b window, the T1a skeleton, the T0 column and the
//! chart position, so a post gets the same value whichever tiles are requested (F9).

use crate::channels::ChannelHit;
use crate::input::{T1a, T1bWindow};
use crate::strata::{CellStrata, Column, StrataModel};
use crate::{k_of, sc_of, smoothstep, DetailParamsV1, R, SECONDS_PER_YEAR, T1B_M};
use serde::Serialize;
use std::collections::HashMap;
use vz_keyed_rng::Scope;

pub const GENERATOR_VERSION_V2: &str = "0.2.0";
pub const PURPOSE_BEDDING_V2: &str = "t2.closure.bedding.v2";
/// Over-class codes (per T1b post and per 2 m post).
pub const OVER_NONE: u8 = 0;
pub const OVER_EDIFICE: u8 = 1;
pub const OVER_DEPOSIT: u8 = 2;
/// Lithology code of deposits and of the edifice (LithologyCodeV1).
pub const LITHO_DEPOSIT: u8 = 10;
pub const LITHO_VOLCANIC: u8 = 3;
/// Pseudo cell id of the edifice column in the strata model.
pub const EDIFICE_CELL: u32 = u32::MAX - 1;

/// V5 hard/soft bed pair medians (m) per LithologyCodeV1 class; (0, 0) = unbedded. Project assumptions (PREREG V5).
pub const BED_PAIR_M: [(f64, f64); 12] = [
    (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (15.0, 3.0), (4.0, 1.0), (0.6, 2.4),
    (3.0, 1.0), (4.0, 1.5), (3.0, 1.0), (2.0, 2.0), (0.0, 0.0), (0.0, 0.0),
];
/// Edifice lava / pyroclastic bed medians (m).
pub const EDIFICE_PAIR_M: (f64, f64) = (10.0, 4.0);

#[derive(Clone, Debug, Serialize)]
pub struct V2Params {
    pub generator_version: String,
    // V1 curves
    pub smooth_passes: u32,
    pub chaikin_passes: u32,
    pub spur_min_posts: usize,
    // V2 valley
    pub valley_on: bool,
    pub valley_rv_m: f64,
    pub valley_rmin_m: f64,
    pub valley_reach_m: f64,
    // V3 banks
    pub bank_sc_alluvium: f64,
    pub bank_reach_m: f64,
    pub bed_floor_radius_m: f64,
    // V4 over-codes
    pub overcodes_on: bool,
    pub edifice_raise_m: f64,
    /// Opt-in local feature support. None retains certified legacy outputs.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub edifice_support: Option<bool>,
    pub edifice_box_half_posts: i64,
    pub edifice_max_depth_m: f64,
    pub deposit_sc: f64,
    pub deposit_k: f64,
    // V5 bedding pairs
    pub bed_hard_mult: f64,
    pub bed_sigma_ln: f64,
    pub bed_clamp_lo: f64,
    pub bed_clamp_hi: f64,
    pub bed_soft_min_frac: f64,
    pub bed_soft_mult_min: f64,
    // V6, V7
    pub slope_clamp_on: bool,
    pub amp_frac_of_relief: f64,
    pub amp_cap_m: f64,
    pub relief_half_posts: i64,
    pub remap_after_relax: bool,
    /// OPT-IN, default off (changes the F6/F8 reference; needs a ruling): T1b channel posts above every free neighbour
    /// (perched beds, request 02) are replaced by their free neighbours' mean in the BASE itself.
    pub deperch_on: bool,
    /// v3 (lane T2C, generator 0.3.0): None = v2 (its params JSON, and so its recipe digest, are unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v3: Option<crate::v3::V3Params>,
}

impl Default for V2Params {
    fn default() -> Self {
        V2Params {
            generator_version: GENERATOR_VERSION_V2.into(),
            smooth_passes: 2,
            chaikin_passes: 2,
            spur_min_posts: 4,
            valley_on: true,
            valley_rv_m: 96.0,
            valley_rmin_m: 32.0,
            valley_reach_m: 96.0,
            bank_sc_alluvium: 0.625,
            bank_reach_m: 48.0,
            bed_floor_radius_m: 48.0,
            overcodes_on: true,
            edifice_raise_m: 0.5,
            edifice_support: None,
            edifice_box_half_posts: 4,
            edifice_max_depth_m: 4000.0,
            deposit_sc: 0.625,
            deposit_k: 3.0,
            bed_hard_mult: 2.0,
            bed_sigma_ln: 0.5,
            bed_clamp_lo: 0.3,
            bed_clamp_hi: 3.0,
            bed_soft_min_frac: 0.1,
            bed_soft_mult_min: 0.15,
            slope_clamp_on: true,
            amp_frac_of_relief: 0.04,
            amp_cap_m: 8.0,
            relief_half_posts: 4,
            remap_after_relax: true,
            deperch_on: false,
            v3: None,
        }
    }
}

// ------------------------------------------------------------------ V5 bed pairs

/// Hard/soft bed pairs from the top of a column of layers (lith, base depth). Returns (bed base depths, Sc, K).
/// Each pair's thickness-weighted mean of 1/Sc equals the class value (soft multiplier solved, clamped).
pub fn beds_v2(lith: &[u8], base: &[f64], pair_of: &dyn Fn(u8) -> (f64, f64), scope: &Scope, tag: &str, max_depth: f64, v: &V2Params) -> R<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let (mut bb, mut bs, mut bk) = (Vec::new(), Vec::new(), Vec::new());
    let mut top = 0.0;
    for (k, (&l, &b)) in lith.iter().zip(base.iter()).enumerate() {
        let (hh, hs) = pair_of(l);
        if hh <= 0.0 || top >= max_depth {
            bb.push(b); bs.push(sc_of(l)); bk.push(k_of(l));
            top = b;
            continue;
        }
        let sk = scope.address(&format!("{tag}/layer/{k}")).map_err(|e| e.0)?;
        let mut d = top;
        let mut n = 0u64;
        while d < b {
            let dr = |c: u64| -> R<f64> { Ok(sk.draw(PURPOSE_BEDDING_V2, c, 0).map_err(|e| e.0)?.normal_f64()) };
            let th = (hh * (v.bed_sigma_ln * dr(2 * n)?).exp()).clamp(v.bed_clamp_lo * hh, v.bed_clamp_hi * hh);
            let ts = (hs * (v.bed_sigma_ln * dr(2 * n + 1)?).exp()).clamp(v.bed_clamp_lo * hs, v.bed_clamp_hi * hs).max(v.bed_soft_min_frac * th);
            // soft multiplier: (th / rh + ts / rs) = th + ts
            let rh = v.bed_hard_mult;
            let rs = (ts / (th + ts - th / rh)).max(v.bed_soft_mult_min);
            // hard bed on top of each pair (cap over its recessive base)
            for (t, m) in [(th, rh), (ts, rs)] {
                if d >= b { break; }
                d = (d + t).min(b);
                if d > max_depth { d = b; }
                bb.push(d); bs.push(sc_of(l) * m); bk.push(k_of(l) / m);
            }
            n += 1;
        }
        top = b;
    }
    Ok((bb, bs, bk))
}

/// Bed pairs of the top class continued ABOVE the column datum (terrain higher than the T0 cell mean is still the top
/// layer; v1 gave it one unbedded bed). Returns (bed top heights above the datum, Sc, K), lowest bed first.
pub fn up_beds(top_class: u8, pair: (f64, f64), scope: &Scope, tag: &str, height: f64, v: &V2Params) -> R<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let pf = |_l: u8| pair;
    beds_v2(&[top_class], &[height], &pf, scope, &format!("{tag}/up"), height, v)
}

// ------------------------------------------------------------------ V1 curves

#[derive(Clone, Debug)]
pub struct Seg {
    pub a: (f64, f64),
    pub b: (f64, f64),
    pub za: f64,
    pub zb: f64,
    pub hwa: f64,
    pub hwb: f64,
    pub da: f64,
    pub db: f64,
    pub t1a: bool,
    /// Polyline (stream) id: the valley, carve and bank use each line's NEAREST point only.
    pub line: u32,
    /// Lane INTS2n4 arm_cl: both ends on the marine arm (always false without the key).
    pub arm: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct Near {
    pub r: f64,
    pub bed: f64,
    pub hw: f64,
    pub depth: f64,
    pub t1a: bool,
    pub line: u32,
    /// Lane INTS2n4 arm_cl: the segment lies on the marine arm (always false without the key).
    pub arm: bool,
}

#[derive(Default, Debug, Serialize, Clone)]
pub struct CurveStats {
    pub channel_posts: usize,
    pub pruned_spur_posts: usize,
    pub streams: usize,
    pub segments: usize,
    pub snapped_tributaries: usize,
    pub max_offset_from_raster_m: f64,
    pub mean_offset_from_raster_m: f64,
    pub bed_points_raised_to_local_min: usize,
    pub bed_raise_max_m: f64,
    /// Lane INTS2n4 arm_cl (None without the key, so the stats JSON is unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub arm: Option<crate::chn::ArmStats>,
}

pub struct Curves {
    /// Window post index of each traced channel post -> distance from the post to its own smoothed curve point (m).
    pub own_off: HashMap<usize, f64>,
    pub segs: Vec<Seg>,
    pub polylines: Vec<Vec<[f64; 3]>>,
    /// Raw traced post chains (chart m, T1b height), for diagnostics.
    pub raw_lines: Vec<Vec<[f64; 3]>>,
    grid: HashMap<(i64, i64), Vec<u32>>,
    cell: f64,
    pub stats: CurveStats,
    /// Lane VALLEY-FLOOR-FIX b8_levee: a second segment grid with the levee reach (empty without the key, so every other
    /// query is unchanged).
    grid_lev: HashMap<(i64, i64), Vec<u32>>,
    /// Lane INTS2n4 arm_cl: window-indexed marine arm mask (None without the key).
    pub arm: Option<Vec<bool>>,
}

fn binomial_pass(p: &mut [[f64; 2]]) {
    let k = p.len();
    if k < 3 { return; }
    let q = p.to_vec();
    for i in 1..k - 1 {
        let m = 2.min(i).min(k - 1 - i);
        p[i] = if m >= 2 {
            let f = |c: usize| (q[i - 2][c] + 4.0 * q[i - 1][c] + 6.0 * q[i][c] + 4.0 * q[i + 1][c] + q[i + 2][c]) / 16.0;
            [f(0), f(1)]
        } else {
            let f = |c: usize| (q[i - 1][c] + 2.0 * q[i][c] + q[i + 1][c]) / 4.0;
            [f(0), f(1)]
        };
    }
}

/// Chaikin corner cutting on (x, y, bed, half width, depth, t1a flag) with ends pinned.
fn chaikin6(p: &[[f64; 6]]) -> Vec<[f64; 6]> {
    if p.len() < 3 { return p.to_vec(); }
    let m = p.len() - 1;
    let mut q = Vec::with_capacity(2 * p.len());
    q.push(p[0]);
    for k in 0..m {
        let (a, b) = (p[k], p[k + 1]);
        let mix = |u: f64| { let mut o = [0.0; 6]; for c in 0..6 { o[c] = a[c] + (b[c] - a[c]) * u; } o };
        if k > 0 { q.push(mix(0.25)); }
        if k + 1 < m { q.push(mix(0.75)); }
    }
    q.push(p[m]);
    q
}

impl Curves {
    /// Trace and smooth every channel chain of the whole T1b window (independent of the requested tiles), then
    /// index the segments that come within `reach` of the bbox.
    pub fn build(w: &T1bWindow, t1a: &T1a, p: &DetailParamsV1, v: &V2Params, bbox: (f64, f64, f64, f64)) -> Curves {
        let n = w.n;
        let nn = n * n;
        // lane INTS2n4 arm_cl: marine arm posts are channel posts (chn.rs ARM_CL_RULE); without the key every closure is as before
        let armp = p.chn.as_ref().and_then(|c| c.arm.as_ref());
        let mut arm_st = crate::chn::ArmStats::default();
        let (arm, arm_d) = match armp { Some(_) => { let (a, d) = crate::chn::arm_mask(w, &mut arm_st); (Some(a), Some(d)) } None => (None, None) };
        let is_arm = |i: usize| arm.as_ref().map(|a| a[i]).unwrap_or(false);
        let is_ch = |i: usize| w.fixed[i] == 1 || w.fixed[i] == 5 || is_arm(i);
        let is_water = |i: usize| (w.fixed[i] == 2 || w.fixed[i] == 4) && !is_arm(i);
        let zmm = |i: usize| (w.height_m[i] * 1000.0).round() as i64;
        // key (z, -area, gy, gx) as a sortable tuple; an arm post's second element is 1 + its steps to open water
        let key = |i: usize| (zmm(i), match &arm_d { Some(d) if is_arm(i) => 1 + d[i].min(u32::MAX - 1) as i64, _ => -(w.drainage_area_m2[i] as f64 * 16.0).round() as i64 }, (i / n) as i64, (i % n) as i64);
        let nb = |i: usize, dx: i64, dy: i64| -> Option<usize> {
            let (x, y) = ((i % n) as i64 + dx, (i / n) as i64 + dy);
            if x < 0 || y < 0 || x >= n as i64 || y >= n as i64 { None } else { Some(y as usize * n + x as usize) }
        };
        const NONE: usize = usize::MAX;
        let mut rcv = vec![NONE; nn];
        let mut st = CurveStats::default();
        for i in 0..nn {
            if !is_ch(i) { continue; }
            st.channel_posts += 1;
            let ki = key(i);
            let mut best: Option<((i64, i64, i64, i64), usize)> = None;
            for dy in -1..=1 { for dx in -1..=1 {
                if dx == 0 && dy == 0 { continue; }
                if let Some(j) = nb(i, dx, dy) {
                    if !(is_ch(j) || is_water(j)) { continue; }
                    let kj = key(j);
                    if kj < ki && best.map(|b| kj < b.0).unwrap_or(true) { best = Some((kj, j)); }
                }
            } }
            if let Some((_, j)) = best { rcv[i] = j; }
        }
        let mut alive: Vec<bool> = (0..nn).map(is_ch).collect();
        let donors_of = |alive: &Vec<bool>| -> Vec<u32> {
            let mut d = vec![0u32; nn];
            for i in 0..nn { if alive[i] && rcv[i] != NONE && alive.get(rcv[i]).copied().unwrap_or(false) { d[rcv[i]] += 1; } }
            d
        };
        // spur pruning: heads whose chain reaches a junction within spur_min_posts posts
        let don = donors_of(&alive);
        for h in 0..nn {
            if !alive[h] || don[h] != 0 { continue; }
            let mut path = vec![h];
            let mut cur = h;
            let mut spur = false;
            while path.len() < v.spur_min_posts {
                let r = rcv[cur];
                if r == NONE || !is_ch(r) { break; }
                if don[r] >= 2 { spur = true; break; }
                path.push(r);
                cur = r;
            }
            if spur { for &q in &path { alive[q] = false; st.pruned_spur_posts += 1; } }
        }
        let don = donors_of(&alive);
        let mut main = vec![NONE; nn];
        for i in 0..nn {
            if !alive[i] || rcv[i] == NONE { continue; }
            let r = rcv[i];
            if !alive[r] { continue; }
            let better = match main[r] {
                NONE => true,
                m => { let (am, ai) = (w.drainage_area_m2[m], w.drainage_area_m2[i]); ai > am || (ai == am && i < m) }
            };
            if better { main[r] = i; }
        }
        let pos = |i: usize| -> [f64; 2] { [((i % n) as i64 + w.gx0) as f64 * T1B_M, ((i / n) as i64 + w.gy0) as f64 * T1B_M] };
        // streams (raw post lists); the last element may be a junction / water post owned by another stream
        let mut streams: Vec<(Vec<usize>, bool)> = Vec::new(); // (posts, last_is_foreign)
        for h in 0..nn {
            if !alive[h] || don[h] != 0 { continue; }
            let mut posts = vec![h];
            let mut cur = h;
            let mut foreign = false;
            loop {
                let r = rcv[cur];
                if r == NONE { break; }
                if is_water(r) || !alive[r] { posts.push(r); foreign = true; break; }
                if main[r] != cur { posts.push(r); foreign = true; break; }
                posts.push(r);
                cur = r;
                if posts.len() > nn { break; }
            }
            streams.push((posts, foreign));
        }
        st.streams = streams.len();
        // pass 1: smooth with raw ends pinned
        // lane INTS2n4 arm_cl: a stream that contains arm posts takes the arm's smoothing passes
        let passes = |ps: &Vec<usize>, foreign: bool| -> u32 {
            let m = if foreign { ps.len() - 1 } else { ps.len() };
            match armp { Some(a) if ps[..m].iter().any(|&i| is_arm(i)) => a.smooth_passes, _ => v.smooth_passes } };
        let smooth = |pts: &mut Vec<[f64; 2]>, np: u32| { for _ in 0..np { binomial_pass(pts); } };
        let mut sm: Vec<Vec<[f64; 2]>> = streams.iter().map(|(ps, foreign)| { let mut q: Vec<[f64; 2]> = ps.iter().map(|&i| pos(i)).collect(); smooth(&mut q, passes(ps, *foreign)); q }).collect();
        // smoothed position of each post at its own (non-foreign) occurrence
        let mut own: HashMap<usize, [f64; 2]> = HashMap::new();
        for (s, (ps, foreign)) in streams.iter().enumerate() {
            let m = if *foreign { ps.len() - 1 } else { ps.len() };
            for k in 0..m { own.insert(ps[k], sm[s][k]); }
        }
        // pass 2: tributaries re-smoothed with the end snapped to the main stream's smoothed position
        for (s, (ps, foreign)) in streams.iter().enumerate() {
            if !*foreign { continue; }
            if let Some(&e) = own.get(ps.last().unwrap()) {
                let mut q: Vec<[f64; 2]> = ps.iter().map(|&i| pos(i)).collect();
                *q.last_mut().unwrap() = e;
                smooth(&mut q, passes(ps, *foreign));
                sm[s] = q;
                st.snapped_tributaries += 1;
            }
        }
        let q_of = |i: usize| -> f64 {
            let (gx, gy) = ((i % n) as i64 + w.gx0, (i / n) as i64 + w.gy0);
            if w.fixed[i] == 1 || is_arm(i) { if let Some(q) = t1a.channel_q_near(gx, gy) { return q; } }
            w.drainage_area_m2[i] as f64 * p.runoff_m_per_yr / SECONDS_PER_YEAR
        };
        let hg = |q: f64| -> (f64, f64) {
            (0.5 * (p.hg_width_coeff * q.powf(p.hg_width_exp)).max(p.min_width_m), (p.hg_depth_coeff * q.powf(p.hg_depth_exp)).max(p.min_depth_m))
        };
        let mut segs = Vec::new();
        let mut polylines = Vec::new();
        let mut own_off: HashMap<usize, f64> = HashMap::new();
        for (s, (ps, foreign)) in streams.iter().enumerate() {
            let m = if *foreign { ps.len() - 1 } else { ps.len() };
            for k in 0..m { own_off.insert(ps[k], (sm[s][k][0] - pos(ps[k])[0]).hypot(sm[s][k][1] - pos(ps[k])[1])); }
        }
        let (mut off_sum, mut off_n, mut off_max) = (0.0f64, 0usize, 0.0f64);
        let (mut raised_n, mut raised_max) = (0usize, 0.0f64);
        for (s, (ps, foreign)) in streams.iter().enumerate() {
            let mut bed = f64::INFINITY;
            let mut pts: Vec<[f64; 6]> = Vec::with_capacity(ps.len());
            for (k, &i) in ps.iter().enumerate() {
                // a water end point is not a bed (T1b holds construction water at sea - F there): keep the last land bed
                if !is_water(i) || k == 0 { bed = bed.min(w.height_m[i]); }
                let qi = if is_ch(i) { q_of(i) } else { q_of(ps[k.saturating_sub(1)]) };
                let (hw, d) = hg(qi);
                let o = (sm[s][k][0] - pos(i)[0]).hypot(sm[s][k][1] - pos(i)[1]);
                off_sum += o; off_n += 1; off_max = off_max.max(o);
                // flag column: 1 = T1a channel post, 2 = marine arm post (lane INTS2n4 arm_cl; only with the key)
                pts.push([sm[s][k][0], sm[s][k][1], bed, hw, d, if is_arm(i) { 2.0 } else if w.fixed[i] == 1 { 1.0 } else { 0.0 }]);
            }
            // lane INTS2n4 arm_cl: funnel-shaped tidal width seaward of the arm head (ARM_CL_RULE)
            if let Some(a) = armp {
                let m_own = if *foreign { ps.len() - 1 } else { ps.len() };
                if let Some(h) = ps[..m_own].iter().position(|&i| is_arm(i)) {
                    arm_st.arm_streams += 1;
                    let mut u = 0.0f64;
                    let hr = pts[h][3];
                    let lo = if arm_st.arm_streams == 1 { hr } else { arm_st.hw_river_m_min_max[0].min(hr) };
                    arm_st.hw_river_m_min_max = [lo, arm_st.hw_river_m_min_max[1].max(hr)];
                    for k in h..pts.len() {
                        if k > h { u += (pts[k][0] - pts[k - 1][0]).hypot(pts[k][1] - pts[k - 1][1]); }
                        let hw = (pts[k][3] * (u / a.b_m).exp()).min(a.hw_max_m.max(pts[k][3]));
                        pts[k][3] = hw;
                        arm_st.hw_tidal_m_max = arm_st.hw_tidal_m_max.max(hw);
                    }
                    arm_st.arm_length_m_max = arm_st.arm_length_m_max.max(u);
                    arm_st.streams.push([pts[h][0], pts[h][1], q_of(ps[h]), hr, u, pts[pts.len() - 1][3]]);
                }
            }
            for _ in 0..v.chaikin_passes { pts = chaikin6(&pts); }
            // the valley bottom T1b actually has there: a curve point's bed is never below the lowest T1b post within
            // bed_floor_radius_m (a perched or up-and-down T1b stretch is not cut down to the running minimum)
            for q in pts.iter_mut() {
                let rr = v.bed_floor_radius_m;
                let (gxa, gxb) = (((q[0] - rr) / T1B_M).ceil() as i64, ((q[0] + rr) / T1B_M).floor() as i64);
                let (gya, gyb) = (((q[1] - rr) / T1B_M).ceil() as i64, ((q[1] + rr) / T1B_M).floor() as i64);
                let mut lo = f64::INFINITY;
                for gy in gya..=gyb { for gx in gxa..=gxb {
                    if (gx as f64 * T1B_M - q[0]).hypot(gy as f64 * T1B_M - q[1]) > rr { continue; }
                    if let Some(i) = w.idx(gx, gy) { if !is_water(i) { lo = lo.min(w.height_m[i]); } }
                } }
                if lo.is_finite() && lo > q[2] { raised_max = raised_max.max(lo - q[2]); q[2] = lo; raised_n += 1; }
            }
            polylines.push(pts.iter().map(|q| [q[0], q[1], q[2]]).collect());
            for wv in pts.windows(2) {
                let (a, b) = (wv[0], wv[1]);
                segs.push(Seg { a: (a[0], a[1]), b: (b[0], b[1]), za: a[2], zb: b[2], hwa: a[3], hwb: b[3], da: a[4], db: b[4], t1a: a[5] + b[5] > 1.0, arm: a[5] + b[5] > 3.0, line: s as u32 });
            }
        }
        st.segments = segs.len();
        st.bed_points_raised_to_local_min = raised_n;
        st.bed_raise_max_m = raised_max;
        if armp.is_some() { st.arm = Some(arm_st); }
        st.mean_offset_from_raster_m = if off_n > 0 { off_sum / off_n as f64 } else { 0.0 };
        st.max_offset_from_raster_m = off_max;
        let cell = 32.0;
        let reach = v.valley_reach_m.max(v.bank_reach_m + 100.0);
        let (x0, y0, x1, y1) = bbox;
        let mut grid: HashMap<(i64, i64), Vec<u32>> = HashMap::new();
        let mut grid_lev: HashMap<(i64, i64), Vec<u32>> = HashMap::new();
        if let Some(q) = p.vff.as_ref().filter(|q| q.levee) {
            let reach = q.levee_width_max_m + 120.0;
            for (k, sgm) in segs.iter().enumerate() {
                let (mnx, mxx) = (sgm.a.0.min(sgm.b.0) - reach, sgm.a.0.max(sgm.b.0) + reach);
                let (mny, mxy) = (sgm.a.1.min(sgm.b.1) - reach, sgm.a.1.max(sgm.b.1) + reach);
                if mxx < x0 || mnx > x1 || mxy < y0 || mny > y1 { continue; }
                for cy in (mny.max(y0 - cell) / cell).floor() as i64..=(mxy.min(y1 + cell) / cell).floor() as i64 {
                    for cx in (mnx.max(x0 - cell) / cell).floor() as i64..=(mxx.min(x1 + cell) / cell).floor() as i64 {
                        grid_lev.entry((cx, cy)).or_default().push(k as u32);
                    }
                }
            }
        }
        for (k, sgm) in segs.iter().enumerate() {
            let (mnx, mxx) = (sgm.a.0.min(sgm.b.0) - reach, sgm.a.0.max(sgm.b.0) + reach);
            let (mny, mxy) = (sgm.a.1.min(sgm.b.1) - reach, sgm.a.1.max(sgm.b.1) + reach);
            if mxx < x0 || mnx > x1 || mxy < y0 || mny > y1 { continue; }
            for cy in (mny.max(y0 - cell) / cell).floor() as i64..=(mxy.min(y1 + cell) / cell).floor() as i64 {
                for cx in (mnx.max(x0 - cell) / cell).floor() as i64..=(mxx.min(x1 + cell) / cell).floor() as i64 {
                    grid.entry((cx, cy)).or_default().push(k as u32);
                }
            }
        }
        let raw_lines = streams.iter().map(|(ps, _)| ps.iter().map(|&i| { let q = pos(i); [q[0], q[1], w.height_m[i]] }).collect()).collect();
        Curves { own_off, segs, polylines, raw_lines, grid, cell, stats: st, grid_lev, arm }
    }

    /// Lane INTS2n4 arm_cl: is window post `i` a marine arm post (false without the key).
    #[inline]
    pub fn is_arm(&self, i: usize) -> bool { self.arm.as_ref().map(|a| a[i]).unwrap_or(false) }

    /// Visit every segment within `reach` m of (x, y) (the grid holds segments within the build reach).
    #[inline]
    pub fn near(&self, x: f64, y: f64, reach: f64, f: &mut dyn FnMut(Near)) {
        let Some(list) = self.grid.get(&((x / self.cell).floor() as i64, (y / self.cell).floor() as i64)) else { return };
        for &k in list {
            let s = &self.segs[k as usize];
            let (vx, vy) = (s.b.0 - s.a.0, s.b.1 - s.a.1);
            let l2 = vx * vx + vy * vy;
            let t = if l2 > 0.0 { (((x - s.a.0) * vx + (y - s.a.1) * vy) / l2).clamp(0.0, 1.0) } else { 0.0 };
            let r = (x - (s.a.0 + t * vx)).hypot(y - (s.a.1 + t * vy));
            if r > reach { continue; }
            f(Near { r, bed: s.za + t * (s.zb - s.za), hw: s.hwa + t * (s.hwb - s.hwa), depth: s.da + t * (s.db - s.da), t1a: s.t1a, line: s.line, arm: s.arm });
        }
    }

    /// Lane VALLEY-FLOOR-FIX b8_levee_v1: the largest levee thickness at (x, y) over the lines near it (0 without the key).
    /// `surf(x, y)` = the T1b surface (bicubic) for the bank height above the bed.
    pub fn levee_at(&self, x: f64, y: f64, q: &crate::VffDetailParamsV1, surf: &dyn Fn(f64, f64) -> f64) -> f64 {
        let Some(list) = self.grid_lev.get(&((x / self.cell).floor() as i64, (y / self.cell).floor() as i64)) else { return 0.0 };
        let mut best: Vec<(u32, f64, f64, f64, f64, (f64, f64), (f64, f64))> = Vec::new(); // line, r, hw, depth, bed, point, unit normal
        for &k in list {
            let s = &self.segs[k as usize];
            // lane INTS2n4 arm_cl: the tidal arm builds no levee (ARM_CL_RULE)
            if s.arm { continue; }
            let (vx, vy) = (s.b.0 - s.a.0, s.b.1 - s.a.1);
            let l2 = vx * vx + vy * vy;
            let t = if l2 > 0.0 { (((x - s.a.0) * vx + (y - s.a.1) * vy) / l2).clamp(0.0, 1.0) } else { 0.0 };
            let r = (x - (s.a.0 + t * vx)).hypot(y - (s.a.1 + t * vy));
            let (hw, dep) = (s.hwa + t * (s.hwb - s.hwa), s.da + t * (s.db - s.da));
            let l = l2.sqrt().max(1e-9);
            let rec = (s.line, r, hw, dep, s.za + t * (s.zb - s.za), (s.a.0 + t * vx, s.a.1 + t * vy), (-vy / l, vx / l));
            match best.iter_mut().find(|o| o.0 == s.line) {
                Some(o) => { if r < o.1 { *o = rec; } }
                None => best.push(rec),
            }
        }
        let mut lev = 0.0f64;
        for (_, r, hw, dep, bed, pt, nrm) in best {
            if r <= hw || hw <= 0.0 || dep <= 0.0 { continue; }
            // bank height above the bed: the carve also removes the T1b surface standing above the bed at the banks
            let hb = 0.5 * (surf(pt.0 + hw * nrm.0, pt.1 + hw * nrm.1) + surf(pt.0 - hw * nrm.0, pt.1 - hw * nrm.1)) - bed;
            let (e, h) = q.levee_shape(hw, dep, hb.max(0.0));
            let u = (r - hw) / e;
            if u < 1.0 { lev = lev.max(h * (1.0 - u) * (1.0 - u)); }
        }
        lev
    }

    /// The nearest point of each polyline within `reach` (a line's other segments are not separate channels).
    pub fn near_lines(&self, x: f64, y: f64, reach: f64) -> Vec<Near> {
        let mut out: Vec<Near> = Vec::new();
        self.near(x, y, reach, &mut |q: Near| {
            match out.iter_mut().find(|o| o.line == q.line) {
                Some(o) => { if q.r < o.r { *o = q; } }
                None => out.push(q),
            }
        });
        out
    }

    /// v1-compatible channel hit (interior parabola, corridor, excess) from the curves.
    pub fn hit(&self, x: f64, y: f64, p: &DetailParamsV1) -> ChannelHit {
        let mut h = ChannelHit { excess_m: 1e9, ..Default::default() };
        // T2D thalweg switch: a channel narrower than the 2 m post grid is carved at least `hmin` half-width, so its interior
        // posts form an 8-connected line instead of a dotted one (corridor width unchanged)
        let hmin = p.v3p().and_then(|w| w.thalweg_hw_m).unwrap_or(0.0);
        for q in self.near_lines(x, y, 250.0) {
            let rc = q.hw + p.bank_width_m;
            h.excess_m = h.excess_m.min(q.r - rc);
            if q.r <= rc { h.corridor = true; }
            let hwi = q.hw.max(hmin);
            if q.arm && q.r <= q.hw + p.bank_width_m.max(48.0) { h.arm = true; }
            if q.r <= hwi {
                let v = q.bed - q.depth * (1.0 - (q.r / hwi).powi(2));
                if h.interior.map(|u| v < u).unwrap_or(true) { h.interior = Some(v); h.t1a = q.t1a; }
            }
        }
        h
    }

    /// Nearest-curve distance (m) and that curve's (bed, hw) within `reach`, if any.
    pub fn nearest(&self, x: f64, y: f64, reach: f64) -> Option<Near> {
        let mut best: Option<Near> = None;
        self.near(x, y, reach, &mut |q: Near| { if best.map(|b| q.r < b.r).unwrap_or(true) { best = Some(q); } });
        best
    }
}

// ------------------------------------------------------------------ T1b-lattice auxiliary fields

/// Per-T1b-post auxiliary fields over a box of the lattice.
pub struct Aux {
    pub gx0: i64,
    pub gy0: i64,
    pub nx: usize,
    pub ny: usize,
    pub over: Vec<u8>,
    /// De-slotted heights (V2): slot channel posts filled to the mean of their banks.
    pub zmod: Vec<f64>,
    /// The base lattice: T1b heights, or with deperch_on the perched channel posts replaced by their free neighbours' mean.
    pub zbase: Vec<f64>,
    /// Valley side gradient g (NaN where undefined).
    pub g: Vec<f64>,
    /// Edifice bedding datum S (box mean of the T1b height).
    pub s_ed: Vec<f64>,
    /// T1b local relief (max - min over the relief box).
    pub relief: Vec<f64>,
    pub counts: AuxCounts,
    /// Lane DATUM (ed_datum=1): the pre-wave-cut T1b surface at the aux posts (None = off).
    pub zpre: Option<Vec<f64>>,
}

#[derive(Default, Debug, Serialize, Clone)]
pub struct AuxCounts {
    pub posts: usize,
    pub land_posts: usize,
    pub edifice: usize,
    pub deposit: usize,
    pub strata_agree: usize,
    pub strata_compared: usize,
    pub deslotted_channel_posts: usize,
    /// T1b channel posts higher than every free 8-neighbour (perched beds: a T1b defect, not carved by v2).
    pub perched_channel_posts: usize,
    pub perched_max_m: f64,
}

/// The producer clears volcano_id outside its constructed volcanic support.
/// Do not infer a body from height, lithology, or distance to a vent.
fn supported_feature(w: &T1bWindow, wi: usize) -> Option<u8> {
    let id = *w.volcano_id.as_ref()?.get(wi)?;
    if id == 255 { return None; }
    let feature = w.volcanoes.get(id as usize)?;
    (!feature.id.is_empty() && !feature.kind.is_empty()).then_some(id)
}

impl Aux {
    #[inline]
    pub fn idx(&self, gx: i64, gy: i64) -> Option<usize> {
        let (x, y) = (gx - self.gx0, gy - self.gy0);
        if x < 0 || y < 0 || x >= self.nx as i64 || y >= self.ny as i64 { None } else { Some(y as usize * self.nx + x as usize) }
    }
    #[inline]
    fn clamp_get(&self, f: &[f64], gx: i64, gy: i64) -> f64 {
        let x = (gx - self.gx0).clamp(0, self.nx as i64 - 1) as usize;
        let y = (gy - self.gy0).clamp(0, self.ny as i64 - 1) as usize;
        f[y * self.nx + x]
    }

    /// Build over the T1b post box [gx0, gx1] x [gy0, gy1] (inclusive).
    /// `pre` (lane DATUM ed_datum=1): a full-window lattice surface (the pre-wave-cut T1b) for the edifice datum S; None = height.
    #[allow(clippy::too_many_arguments)]
    pub fn build(w: &T1bWindow, strata: &StrataModel, curves: &Curves, v: &V2Params, gx0: i64, gy0: i64, gx1: i64, gy1: i64, pre: Option<&[f64]>) -> Aux {
        let (nx, ny) = ((gx1 - gx0 + 1) as usize, (gy1 - gy0 + 1) as usize);
        let np = nx * ny;
        let mut c = AuxCounts { posts: np, ..Default::default() };
        let mut over = vec![OVER_NONE; np];
        let mut g = vec![f64::NAN; np];
        let mut zmod = vec![0.0f64; np];
        let mut zbase = vec![0.0f64; np];
        let mut s_ed = vec![0.0f64; np];
        let mut relief = vec![0.0f64; np];
        let base_m = |i: usize| -> f64 { w.base_m.as_ref().map(|b| b[i]).unwrap_or(w.height_m[i]) };
        // the datum surface: the published height, or (ed_datum) the pre-cut surface, clamped to the lattice as w.h
        let hd = |gx: i64, gy: i64| -> f64 { match pre {
            None => w.h(gx, gy),
            Some(f) => { let x = (gx - w.gx0).clamp(0, w.n as i64 - 1) as usize; let y = (gy - w.gy0).clamp(0, w.n as i64 - 1) as usize; f[y * w.n + x] }
        } };
        let mut zpre = pre.map(|_| vec![0.0f64; np]);
        let mut nearc: Vec<Option<Near>> = vec![None; np];
        for j in 0..ny { for i in 0..nx {
            let k = j * nx + i;
            let (gx, gy) = (gx0 + i as i64, gy0 + j as i64);
            let (x, y) = (gx as f64 * T1B_M, gy as f64 * T1B_M);
            let wi = w.idx(gx, gy);
            let z = w.h(gx, gy);
            zmod[k] = z;
            zbase[k] = z;
            // over-codes
            if let Some(wi) = wi {
                let code = w.fixed[wi];
                if code != 2 && code != 4 {
                    c.land_posts += 1;
                    let t1b_l = w.exposed_litho.as_ref().map(|e| e[wi]).unwrap_or(255);
                    let own = strata.cells.get(&w.t0_cell[wi]).map(|cs| strata.exposed(cs, x, y, z).1).unwrap_or(255);
                    // Height above the construction baseline is not evidence of a volcanic body.
                    // Require the producer's local feature support and a resolvable record.
                    // Native volcanic columns can be raised by the hillslope/fill stages.
                    let supported_edifice = supported_feature(w, wi).is_some();
                    if v.overcodes_on && (v.edifice_support != Some(true) || supported_edifice) && t1b_l == LITHO_VOLCANIC && (own != LITHO_VOLCANIC || z - base_m(wi) > v.edifice_raise_m) { over[k] = OVER_EDIFICE; c.edifice += 1; }
                    else if v.overcodes_on && t1b_l == LITHO_DEPOSIT && own != LITHO_DEPOSIT { over[k] = OVER_DEPOSIT; c.deposit += 1; }
                    else { c.strata_compared += 1; if own == t1b_l { c.strata_agree += 1; } }
                }
            }
            // relief and edifice datum
            let (mut lo, mut hi, mut sum, mut cnt) = (f64::MAX, f64::MIN, 0.0, 0.0);
            let (rh, eh) = (v.relief_half_posts, v.edifice_box_half_posts);
            let hmax = rh.max(eh);
            for dy in -hmax..=hmax { for dx in -hmax..=hmax {
                let zz = w.h(gx + dx, gy + dy);
                if dx.abs() <= rh && dy.abs() <= rh { lo = lo.min(zz); hi = hi.max(zz); }
                if dx.abs() <= eh && dy.abs() <= eh { sum += if pre.is_some() { hd(gx + dx, gy + dy) } else { zz }; cnt += 1.0; }
            } }
            if let Some(zp) = zpre.as_mut() { zp[k] = hd(gx, gy); }
            relief[k] = hi - lo;
            s_ed[k] = sum / cnt;
            nearc[k] = if v.valley_on { curves.nearest(x, y, v.valley_reach_m) } else { None };
            if let (Some(q), Some(wi)) = (nearc[k], wi) {
                if w.fixed[wi] == 0 { g[k] = (z - q.bed).max(0.0) / (q.r - q.hw).max(v.valley_rmin_m); }
            }
        } }
        if v.valley_on {
            let g0 = g.clone();
            for j in 0..ny { for i in 0..nx {
                let k = j * nx + i;
                let (gx, gy) = (gx0 + i as i64, gy0 + j as i64);
                let Some(wi) = w.idx(gx, gy) else { continue };
                if w.fixed[wi] != 1 && w.fixed[wi] != 5 && !curves.is_arm(wi) { continue; }
                let Some(_q) = nearc[k] else { continue };
                let (mut s, mut m) = (0.0, 0.0);
                let (mut zs, mut zn, mut zmax) = (0.0, 0.0, f64::MIN);
                for dy in -1i64..=1 { for dx in -1i64..=1 {
                    let (a, b) = (i as i64 + dx, j as i64 + dy);
                    if a < 0 || b < 0 || a >= nx as i64 || b >= ny as i64 { continue; }
                    let v0 = g0[b as usize * nx + a as usize];
                    if v0.is_finite() { s += v0; m += 1.0; }
                    if let Some(wj) = w.idx(gx + dx, gy + dy) {
                        if w.fixed[wj] == 0 { let zz = w.height_m[wj]; zs += zz; zn += 1.0; zmax = zmax.max(zz); }
                    }
                } }
                let gg = if m > 0.0 { s / m } else { 0.0 };
                g[k] = gg;
                let zc = w.height_m[wi];
                if zn > 0.0 && zc > zmax {
                    c.perched_channel_posts += 1; c.perched_max_m = c.perched_max_m.max(zc - zmax);
                    if v.deperch_on { zbase[k] = zs / zn; }
                }
                // the raster slot is filled to its banks (mean of the free 8-neighbours; this also caps a perched T1b bed);
                // the valley is then re-cut along the smoothed curve in the 2 m valley stage (Ctx::valley)
                let _ = curves.own_off.get(&wi);
                // only a slot (below the mean of its banks) is filled; a channel post at or above its banks (perched, a T1b
                // defect) keeps its T1b height: T2 does not remove resolved T1b relief (F6, condition 7)
                zmod[k] = if zn > 0.0 && zc < zs / zn { zs / zn } else { zbase[k] };
                c.deslotted_channel_posts += 1;
            } }
        }
        Aux { gx0, gy0, nx, ny, over, zmod, zbase, g, s_ed, relief, counts: c, zpre }
    }

    /// Catmull-Rom bicubic of a lattice field at chart (x, y) m.
    pub fn bicubic(&self, f: &[f64], x: f64, y: f64) -> f64 {
        let fx = x / T1B_M;
        let fy = y / T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let cr = |t: f64| { let (t2, t3) = (t * t, t * t * t); [(-t3 + 2.0 * t2 - t) * 0.5, (3.0 * t3 - 5.0 * t2 + 2.0) * 0.5, (-3.0 * t3 + 4.0 * t2 + t) * 0.5, (t3 - t2) * 0.5] };
        let (wx, wy) = (cr(tx), cr(ty));
        let mut acc = 0.0;
        for (j, wyj) in wy.iter().enumerate() {
            let mut row = 0.0;
            for (i, wxi) in wx.iter().enumerate() { row += wxi * self.clamp_get(f, ix - 1 + i as i64, iy - 1 + j as i64); }
            acc += wyj * row;
        }
        acc
    }

    /// Bilinear of a lattice field; NaN entries are skipped (weights renormalized); None if all four are NaN.
    pub fn bilinear(&self, f: &[f64], x: f64, y: f64) -> Option<f64> {
        let fx = x / T1B_M;
        let fy = y / T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let (mut s, mut ws) = (0.0, 0.0);
        for (dx, dy, wgt) in [(0, 0, (1.0 - tx) * (1.0 - ty)), (1, 0, tx * (1.0 - ty)), (0, 1, (1.0 - tx) * ty), (1, 1, tx * ty)] {
            let v0 = self.clamp_get(f, ix + dx, iy + dy);
            if v0.is_finite() && wgt > 0.0 { s += wgt * v0; ws += wgt; }
        }
        if ws > 0.0 { Some(s / ws) } else if tx == 0.0 && ty == 0.0 { None } else { None }
    }

    /// Dominant actually contributing feature, using the same bilinear support
    /// as over_w. This never borrows an unrelated nearby feature for provenance.
    pub fn supported_feature_at(&self, w: &T1bWindow, x: f64, y: f64) -> Option<u8> {
        let (fx, fy) = (x / T1B_M, y / T1B_M);
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let mut weights = std::collections::BTreeMap::<u8, f64>::new();
        for (dx, dy, weight) in [(0, 0, (1.0-tx)*(1.0-ty)), (1, 0, tx*(1.0-ty)), (0, 1, (1.0-tx)*ty), (1, 1, tx*ty)] {
            if weight <= 0.0 { continue; }
            let gx = (ix + dx).clamp(self.gx0, self.gx0 + self.nx as i64 - 1);
            let gy = (iy + dy).clamp(self.gy0, self.gy0 + self.ny as i64 - 1);
            if self.over[self.idx(gx, gy)?] != OVER_EDIFICE { continue; }
            if let Some(id) = w.idx(gx, gy).and_then(|wi| supported_feature(w, wi)) {
                *weights.entry(id).or_default() += weight;
            }
        }
        let mut best: Option<(u8, f64)> = None;
        for (id, weight) in weights {
            if best.map(|(_, b)| weight > b + 1e-12).unwrap_or(true) { best = Some((id, weight)); }
        }
        best.map(|(id, _)| id)
    }

    /// Bilinear indicator weights (edifice, deposit) at (x, y).
    pub fn over_w(&self, x: f64, y: f64) -> (f64, f64) {
        let fx = x / T1B_M;
        let fy = y / T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let (mut e, mut d) = (0.0, 0.0);
        for (dx, dy, wgt) in [(0, 0, (1.0 - tx) * (1.0 - ty)), (1, 0, tx * (1.0 - ty)), (0, 1, (1.0 - tx) * ty), (1, 1, tx * ty)] {
            let gx = (ix + dx - self.gx0).clamp(0, self.nx as i64 - 1) as usize;
            let gy = (iy + dy - self.gy0).clamp(0, self.ny as i64 - 1) as usize;
            match self.over[gy * self.nx + gx] { OVER_EDIFICE => e += wgt, OVER_DEPOSIT => d += wgt, _ => {} }
        }
        (e, d)
    }
}

// ------------------------------------------------------------------ the v2 lithology model at a 2 m post

/// What a 2 m post is made of.
#[derive(Clone, Copy, Debug)]
pub struct OverAt { pub w_ed: f64, pub w_dep: f64, /** lane BS: dominant volcano index of the edifice part (beds=2; 255 otherwise) */ pub vd: u8 }
impl OverAt {
    #[inline] pub fn class(&self) -> u8 { if self.w_dep >= 0.5 { OVER_DEPOSIT } else if self.w_ed >= 0.5 { OVER_EDIFICE } else { OVER_NONE } }
}

/// Edifice column (all class 3) with bed pairs by depth below the datum S.
pub fn edifice_strata(scope: &Scope, p: &DetailParamsV1, v: &V2Params) -> R<CellStrata> {
    let depth = v.edifice_max_depth_m;
    let lith = vec![LITHO_VOLCANIC];
    let base = vec![depth];
    let pair = |_l: u8| EDIFICE_PAIR_M;
    let (bb, mut bs, mut bk) = match &v.v3 {
        Some(v3) if v3.pkg_on => crate::v3::beds_v3(&lith, &base, &pair, &|_l: u8| crate::v3::PKG_MEDIAN_EDIFICE_M, scope, "edifice", depth, v, v3)?,
        _ => beds_v2(&lith, &base, &pair, scope, "edifice", depth, v)?,
    };
    let (ut, us, uk) = match &v.v3 {
        Some(v3) if v3.pkg_on => crate::v3::beds_v3(&[LITHO_VOLCANIC], &[3000.0], &pair, &|_l: u8| crate::v3::PKG_MEDIAN_EDIFICE_M, scope, "edifice/up", 3000.0, v, v3)?,
        _ => up_beds(LITHO_VOLCANIC, EDIFICE_PAIR_M, scope, "edifice", 3000.0, v)?,
    };
    let n_down = bb.len();
    bs.extend(us); bk.extend(uk);
    let col = Column { lith, base_depth: base, bed_base: bb, bed_sc: bs, bed_k: bk, up_top: ut, n_down, bed_fac: vec![], bed_unit: vec![], unit_top: vec![] };
    let mut cs = CellStrata { cell: EDIFICE_CELL, top_m: 0.0, g: (0.0, 0.0), g_raw: (0.0, 0.0), col, u0: 0.0, du: p.remap_table_step_m, psi: Vec::new(), max_drift_m: 0.0, sc_win: Vec::new(), k_win: Vec::new(), cum_is: Vec::new() };
    crate::strata::build_remap_pub(&mut cs, p)?;
    Ok(cs)
}

/// Bounded displacement (V7): B tanh(d / B).
#[inline]
pub fn bound(d: f64, b: f64) -> f64 { if b <= 0.0 { 0.0 } else { b * (d / b).tanh() } }

/// V6 remap weight from the base slope G0 and the layer's class Sc.
#[inline]
pub fn slope_alpha(g0: f64, sc_layer: f64, rh: f64) -> f64 {
    if g0 <= sc_layer { return 1.0; }
    ((sc_layer * rh - g0) / (g0 * (rh - 1.0))).clamp(0.0, 1.0)
}

pub fn smooth_w(w: f64) -> f64 { smoothstep((w - 0.25) / 0.5) }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bed_pairs_conserve_class_rock_mass_sc() {
        let v = V2Params::default();
        let scope = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let pair = |l: u8| BED_PAIR_M[l as usize];
        let (bb, bs, _) = beds_v2(&[3, 5], &[460.0, 932.0], &pair, &scope, "cell/test", 7000.0, &v).unwrap();
        // thickness-weighted mean of 1/Sc over the basalt layer equals 1/Sc(basalt) where no soft clamp bound
        let (mut top, mut acc, mut h) = (0.0, 0.0, 0.0);
        for (b, s) in bb.iter().zip(bs.iter()) {
            if *b > 460.0 { break; }
            acc += (b - top) / s; h += b - top; top = *b;
        }
        let mean_inv = acc / h;
        assert!((mean_inv * sc_of(3) - 1.0).abs() < 0.05, "basalt rock-mass 1/Sc ratio {}", mean_inv * sc_of(3));
        assert!(bs.iter().take(10).any(|&s| s > 3.0) && bs.iter().take(10).any(|&s| s < 1.0), "hard and soft beds alternate");
    }

    #[test]
    fn slope_alpha_never_pushes_hard_beds_past_their_limit() {
        for &g0 in &[0.1, 0.5, 1.0, 1.5, 1.732, 2.5, 3.464, 5.0] {
            let a = slope_alpha(g0, 1.732, 2.0);
            assert!((0.0..=1.0).contains(&a));
            // hard-bed slope after the remap <= its own critical slope (2 x 1.732) whenever the base is below it
            if g0 < 2.0 * 1.732 { assert!(g0 * (1.0 - a + 2.0 * a) <= 2.0 * 1.732 + 1e-9, "g0 {g0} a {a}"); }
        }
        assert_eq!(slope_alpha(1.0, 1.732, 2.0), 1.0);
        assert_eq!(slope_alpha(4.0, 1.732, 2.0), 0.0);
    }

    #[test]
    fn bound_is_bounded_and_odd() {
        for &d in &[-100.0, -3.0, 0.0, 2.0, 50.0] { let b = bound(d, 4.0); assert!(b.abs() < 4.0 + 1e-12); assert!((bound(-d, 4.0) + b).abs() < 1e-12); }
        assert_eq!(bound(5.0, 0.0), 0.0);
    }

    #[test]
    fn chaikin_keeps_ends_and_monotone_beds() {
        let p = vec![[0.0, 0.0, 10.0, 1.0, 0.2, 0.0], [32.0, 0.0, 9.0, 1.0, 0.2, 0.0], [32.0, 32.0, 7.0, 1.0, 0.2, 0.0], [64.0, 32.0, 6.0, 1.0, 0.2, 0.0]];
        let mut q = p.clone();
        for _ in 0..2 { q = chaikin6(&q); }
        assert_eq!(q[0], p[0]);
        assert_eq!(*q.last().unwrap(), *p.last().unwrap());
        for w in q.windows(2) { assert!(w[1][2] <= w[0][2] + 1e-12); }
    }
}

#[cfg(test)]
mod support_tests {
    use super::*;
    #[test]
    fn actual_window_feature_support_and_falsifiers() {
        let Ok(dir) = std::env::var("VZ_ED_SUPPORT_WINDOW") else {
            eprintln!("SKIP actual-window support regression: set VZ_ED_SUPPORT_WINDOW"); return;
        };
        let mut w = T1bWindow::load(std::path::Path::new(&dir)).unwrap();
        let i = (0..w.height_m.len()).find(|&i| supported_feature(&w, i).is_some()).unwrap();
        let n = w.volcano_id.as_ref().unwrap().iter().position(|&v| v == 255).unwrap();
        assert!(supported_feature(&w, n).is_none());
        let id = w.volcano_id.as_ref().unwrap()[i];
        w.volcano_id.as_mut().unwrap()[i] = 255;
        assert!(supported_feature(&w, i).is_none());
        w.volcano_id.as_mut().unwrap()[i] = 254;
        assert!(supported_feature(&w, i).is_none());
        w.volcano_id.as_mut().unwrap()[i] = id;
        w.volcanoes[id as usize].id.clear();
        assert!(supported_feature(&w, i).is_none());
        w.volcano_id = None;
        assert!(supported_feature(&w, i).is_none());
    }
    #[test]
    fn support_key_off_retains_exact_params() {
        let mut p = DetailParamsV1::v3();
        let before = serde_json::to_vec(&p).unwrap();
        p.set("ed_support", "0").unwrap();
        assert_eq!(before, serde_json::to_vec(&p).unwrap());
        p.set("ed_support", "1").unwrap();
        assert_ne!(before, serde_json::to_vec(&p).unwrap());
    }
}
