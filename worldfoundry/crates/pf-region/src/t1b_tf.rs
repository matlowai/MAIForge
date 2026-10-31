//! T1b 0.9.4 lane TF (2026-10-15; `unpublished-design-record`): walk-clean fixes on walk config C.
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon; construction (no SimTime). Everything is opt-in under
//! `params.ls.tb.tc.tf` (inside params.ls.tb, so lane TB's paired keying keeps every keyed draw); absent = lane TE's
//! behaviour byte-identical.
//!
//! Diagnosis (measured on TE's kCr0 / kCr1, PREREG section 1):
//! - D1: TD5's soft-rock bluff classifies a post by the rock at its CURRENT height and lowers it to the cap, which can lie
//!   inside a hard unit below the soft cover; hard posts are skipped and stay as a rim, so the lowered interior becomes a
//!   closed basin that the final fill fills flat (1-2.8 km2 basalt plateaus in kC).
//! - D3: the T0 bilinear has slope breaks along T0 cell-centre lines; CF's R1a sloped base carries the break into the cone.
//! - D4: coast distances measured from sea POST centres alternate by up to a post along an oblique shore, so a cliff a post
//!   wide becomes a row of lattice teeth; and every shore is cut alike, sheltered or open.
//! Items: T1 fill guard and basin report (report only), T2 soft_bluff_v2, T3 sub-post shore distance, T4 fetch exposure,
//! T5 C1 cone base.

use crate::t0::T0;
use crate::t1b_rivernet::EvalLat;
use crate::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct TfParamsV1 {
    /// T1: fill guard and closed-basin report (window statistics only; no height changes).
    pub report: u32,
    /// T2: soft_bluff_v2 (flag; replaces TD5's bluff when TD5 is on).
    pub bluff_v2: u32,
    /// T3: coast distances from the sub-post sea-level contour (flag).
    pub subpost: u32,
    /// T4: wave exposure from fetch scales cliff retreat and beach width (flag); rays and fetch cap (mm).
    pub exposure: u32,
    pub fetch_dirs: u32,
    pub fetch_cap_mm: i64,
    /// T5: the cone's R1a sloped base on the C1 (Catmull-Rom) T0 surface (flag).
    pub cone_c1: u32,
    /// REVISION 1: half-width (mm) of the box average of the T0 bilinear used as the cone base.
    pub cone_smooth_mm: i64,
    pub report_rule: String,
    pub bluff_rule: String,
    pub subpost_rule: String,
    pub exposure_rule: String,
    pub cone_rule: String,
}

impl TfParamsV1 {
    pub fn v1() -> Self {
        TfParamsV1 {
            report: 1, bluff_v2: 1, subpost: 1, exposure: 1, fetch_dirs: 16, fetch_cap_mm: 4_000_000, cone_c1: 1, cone_smooth_mm: 2_000_000,
            report_rule: "fill_guard_v1 (report only): both priority floods (pre-coast and final) carry each post's outlet origin (the eval lattice ring or a fixed outlet, inherited from the post it was reached from); per window and fill: ring-origin posts RAISED (a depression filled to a spill on or through the ring) or free BELOW SEA (sub-sea ground drained to the ring), the largest raise, and their smallest distance to the published core + halo box (guard pass = none within 8 km, the TD5 bluff reach); closed basins = 8-connected components of raised posts: count, area, volume, deepest raise, spill height, and those with a pre-fill floor below sea level (no lake model: they stay filled, declared)".into(),
            bluff_rule: "soft_bluff_v2 (replaces soft_bluff_v1's distance and limits): d = Euclidean distance to the shore carried by vector propagation through SOFT land posts only (k_mult > cliff_k_max at the post's height; hard posts, channels, lakes, T1a posts and no-cut fan zones do not pass it), posts settled in increasing d; cap z_b = sea + beach_slope w + Sc (d - w) as v1 (w the beach width at the shore point d comes from); new z = max(z_b, top of the first hard layer under the post between z_b and its height (bench; 2 m scan then 12 bisections), the new height of the post d was carried from); z = min(z, new z): soft ground behind a hard rim is not reached, the cut never enters the hard unit, and no closed basin is cut".into(),
            subpost_rule: "shore_subpost_v1: the shore is the sea-level contour of the surface at the wave-cut step: on every 4-neighbour link between a sea post and a non-sea post the crossing at t = (sea - z_sea) / (z_land - z_sea) (clamped, linear), joined per lattice cell into segments (marching squares; saddles by the 4-corner mean); posts within 2 posts of a segment take the exact distance to it, others by vector propagation (each post keeps its nearest shore point; heap key (d, global gy, gx)); cliff platform/talus/wall, beach and bluff all use this distance".into(),
            exposure_rule: "fetch_exposure_v1 (recalled, unverified: fetch-limited wave height H ~ F^0.5, Hasselmann et al. 1973 JONSWAP; wave power ~ H^2 ~ F; cliff retreat ~ wave power, Sunamura 1992): per shore sea post (a sea post with a non-sea 4-neighbour) fetch_dirs rays at equal angles from 0 (east) over contiguous sea posts (nearest post every spacing), each fetch capped at fetch_cap (a ray leaving the lattice stops there); exposure e = sum of fetches / (fetch_dirs / 2 x fetch_cap), clamped to [0, 1] (a straight open coast = 1, a narrow inlet < 1); a shore segment takes the mean e of its two crossings' sea posts, a land post the e of the shore point its distance comes from; cliff platform width (cliff_retreat x k_mult) and beach width (min(beach_width x k_mult, beach_width_max)) are multiplied by e".into(),
            cone_rule: "REVISION 1 cone_c1_base_v2: CF's R1a sloped base uses the T0 bilinear box-averaged over half-width cone_smooth (t1b_rivernet::box_mean weights clamp(r/dx + 1/2 - |k|, 0, 1) on the lattice posts; at the vent the same separable weights over the global 32 m posts around the vent's nearest post): t = B(post) - B(vent). The bilinear has slope breaks along T0 cell-centre lines that ran down the cone; the box average is C1 across them and equals the bilinear wherever it is linear over the box (unlike a bicubic, which changed the cone's whole footprint: lane TF tA sets)".into(),
        }
    }
    pub fn off() -> Self { let mut p = Self::v1(); p.report = 0; p.bluff_v2 = 0; p.subpost = 0; p.exposure = 0; p.cone_c1 = 0; p }
    pub fn validate(&self) -> R<()> {
        if self.fetch_dirs < 4 || self.fetch_cap_mm <= 0 { return Err("tf: fetch_dirs >= 4 and fetch_cap > 0".into()); }
        Ok(())
    }
}

/// REVISION 1: the box-averaged bilinear on a lattice (B per post).
pub fn smooth_bil(lat: &EvalLat, bil: &[f32], half_m: f64) -> Vec<f32> {
    let v: Vec<f64> = bil.iter().map(|b| *b as f64).collect();
    crate::t1b_rivernet::box_mean(lat, &v, &|_| true, &|_| half_m).iter().map(|b| *b as f32).collect()
}

/// REVISION 1: B at a vent: separable box weights over the global 32 m posts around the vent's nearest post; `bil(gx, gy)` is
/// the T0 bilinear at global 32 m post (gx, gy).
pub fn smooth_bil_at(bil: &dyn Fn(i64, i64) -> R<f64>, x_m: f64, y_m: f64, half_m: f64) -> R<f64> {
    let dx = 32.0;
    let (cx, cy) = ((x_m / dx).round() as i64, (y_m / dx).round() as i64);
    let h = half_m / dx + 0.5;
    let k = h.ceil() as i64;
    let (mut s, mut ws) = (0.0, 0.0);
    for oy in -k..=k { let wy = (h - oy.abs() as f64).clamp(0.0, 1.0); if wy <= 0.0 { continue; }
        for ox in -k..=k { let wx = (h - ox.abs() as f64).clamp(0.0, 1.0); if wx <= 0.0 { continue; }
            s += wx * wy * bil(cx + ox, cy + oy)?; ws += wx * wy; } }
    Ok(s / ws)
}

pub fn tf(p: &crate::t1b_rivernet::RiverNetParamsV1) -> Option<&TfParamsV1> { crate::t1b_ls::tc(p).and_then(|t| t.tf.as_ref()) }

// ------------------------------------------------------------------ T5: C1 T0 surface

#[inline]
fn cr(t: f64) -> [f64; 4] {
    let (t2, t3) = (t * t, t * t * t);
    [0.5 * (-t3 + 2.0 * t2 - t), 0.5 * (3.0 * t3 - 5.0 * t2 + 2.0), 0.5 * (-3.0 * t3 + 4.0 * t2 + t), 0.5 * (t3 - t2)]
}

/// Catmull-Rom bicubic of a T0 cell field at (lat, lon) degrees on the same cell-centre grid as `t1b::bilinear` (rows
/// clamped at the poles, columns wrapped). C1, interpolating at cell centres.
pub fn bicubic_t0(t0: &T0, field: &dyn Fn(u32) -> f64, lat: f64, lon: f64) -> f64 {
    let dlat = 180.0 / t0.nlat as f64;
    let dlon = 360.0 / t0.nlon as f64;
    let fi = (lat + 90.0) / dlat - 0.5;
    let fj = (lon + 180.0) / dlon - 0.5;
    let (i0, j0) = (fi.floor(), fj.floor());
    let (wi, wj) = (cr(fi - i0), cr(fj - j0));
    let mut acc = 0.0;
    for (a, wa) in wi.iter().enumerate() {
        let ii = (i0 as i64 + a as i64 - 1).clamp(0, t0.nlat as i64 - 1) as usize;
        for (b, wb) in wj.iter().enumerate() {
            let jj = (j0 as i64 + b as i64 - 1).rem_euclid(t0.nlon as i64) as usize;
            acc += wa * wb * field((ii * t0.nlon + jj) as u32);
        }
    }
    acc
}

// ------------------------------------------------------------------ T1: fill with origins

#[derive(Clone, Copy)]
struct HK { z: f64, gy: i64, gx: i64, i: u32 }
impl PartialEq for HK { fn eq(&self, o: &Self) -> bool { self.cmp(o) == Ordering::Equal } }
impl Eq for HK {}
impl PartialOrd for HK { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for HK { fn cmp(&self, o: &Self) -> Ordering { o.z.total_cmp(&self.z).then(o.gy.cmp(&self.gy)).then(o.gx.cmp(&self.gx)) } }

/// One closed basin of a fill (an 8-connected component of raised posts).
#[derive(Clone, Debug, Default, Serialize)]
pub struct Basin { pub posts: u32, pub area_m2: f64, pub volume_m3: f64, pub max_raise_m: f64, pub floor_m: f64, pub spill_m: f64, pub x_m: f64, pub y_m: f64, pub ring_origin: bool }

/// Fill guard and basin report of one priority flood (published-box numbers; the box is the lattice minus the margin).
#[derive(Clone, Debug, Default, Serialize)]
pub struct FillReport {
    pub raised_posts: u32,
    pub ring_origin_raised_posts: u32,
    pub ring_origin_subsea_posts: u32,
    pub ring_origin_max_raise_m: f64,
    /// Smallest distance (m) from a flagged post to the published box (-1: none flagged). Flagged (fix B2): ring-origin posts
    /// whose flood path's highest point is the ring seed (the spill is ON the ring), raised or free below sea.
    pub ring_flag_min_dist_to_box_m: f64,
    /// Chart position (m) of that flagged post, and the ring-origin raised posts within the guard distance of the box.
    pub ring_flag_nearest_xy_m: [f64; 2],
    pub ring_raised_within_guard_posts: u32,
    pub ring_subsea_within_guard_posts: u32,
    pub guard_pass: bool,
    /// Basins with a post in the published box: all, and those whose pre-fill floor is below sea level.
    pub basins_in_box: u32,
    pub basins_in_box_posts: u32,
    pub basins_in_box_volume_m3: f64,
    pub subsea_basins_in_box: u32,
    pub subsea_basins_in_box_area_m2: f64,
    pub subsea_basins_in_box_volume_m3: f64,
    /// The 6 largest basins (by area) touching the box.
    pub largest: Vec<Basin>,
}

/// `t1b_rivernet::fill` with origins (identical heights: same seeds, same heap key, same order) and the report.
/// `margin_posts`: the published box is [margin, n - 1 - margin] in both axes (0: the whole lattice).
pub fn fill_report(lat: &EvalLat, z: &mut [f64], outlet: &[bool], sea: f64, fixed_free: &dyn Fn(usize) -> bool, margin_posts: i64, guard_m: f64) -> (u32, f64, FillReport) {
    let n = lat.n;
    let nn = lat.nodes();
    let dx = lat.dx();
    let z0 = z.to_vec();
    let mut done = vec![false; nn];
    let mut ring = vec![false; nn];
    // Fix B2: the spill of a post toward its seed is the highest point on its flood path (before it); a ring-origin post's
    // spill comes THROUGH the ring only when that highest point is the ring seed itself.
    let mut pmax = vec![f64::NEG_INFINITY; nn];
    let mut seedz = vec![f64::NEG_INFINITY; nn];
    let mut heap = BinaryHeap::new();
    for i in 0..nn {
        let (x, y) = (i % n, i / n);
        let on_ring = x == 0 || y == 0 || x == n - 1 || y == n - 1;
        if outlet[i] || on_ring {
            done[i] = true;
            ring[i] = !outlet[i];
            seedz[i] = z[i];
            let (gx, gy) = lat.g(i);
            heap.push(HK { z: z[i], gy, gx, i: i as u32 });
        }
    }
    let (mut raised, mut maxr) = (0u32, 0f64);
    while let Some(HK { z: v, i, .. }) = heap.pop() {
        let i = i as usize;
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if done[j] { continue; }
            done[j] = true;
            ring[j] = ring[i];
            seedz[j] = seedz[i];
            pmax[j] = pmax[i].max(z[i]);
            if z[j] < v { raised += 1; maxr = maxr.max(v - z[j]); z[j] = v; }
            let (gx, gy) = lat.g(j);
            heap.push(HK { z: z[j], gy, gx, i: j as u32 });
        }
    }
    // Report.
    let m = margin_posts.max(0);
    let box_dist = |i: usize| -> f64 {
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let (lo, hi) = (m, n as i64 - 1 - m);
        let ddx = if x < lo { lo - x } else if x > hi { x - hi } else { 0 };
        let ddy = if y < lo { lo - y } else if y > hi { y - hi } else { 0 };
        ((ddx * ddx + ddy * ddy) as f64).sqrt() * dx
    };
    let mut r = FillReport { ring_flag_min_dist_to_box_m: -1.0, ..Default::default() };
    let mut up = vec![false; nn];
    for i in 0..nn {
        if z[i] > z0[i] { up[i] = true; r.raised_posts += 1; }
        if !ring[i] || pmax[i] > seedz[i] + 1e-9 { continue; }
        let flag_raise = up[i];
        let flag_sub = !up[i] && fixed_free(i) && z[i] <= sea;
        if flag_raise { r.ring_origin_raised_posts += 1; r.ring_origin_max_raise_m = r.ring_origin_max_raise_m.max(z[i] - z0[i]); }
        if flag_sub { r.ring_origin_subsea_posts += 1; }
        if flag_raise || flag_sub {
            let d = box_dist(i);
            if d <= guard_m { if flag_raise { r.ring_raised_within_guard_posts += 1; } else { r.ring_subsea_within_guard_posts += 1; } }
            if r.ring_flag_min_dist_to_box_m < 0.0 || d < r.ring_flag_min_dist_to_box_m {
                r.ring_flag_min_dist_to_box_m = d;
                let (gx, gy) = lat.g(i);
                r.ring_flag_nearest_xy_m = [gx as f64 * dx, gy as f64 * dx];
            }
        }
    }
    r.guard_pass = r.ring_flag_min_dist_to_box_m < 0.0 || r.ring_flag_min_dist_to_box_m > guard_m;
    // Basins: 8-connected components of raised posts (deterministic scan order).
    let mut seen = vec![false; nn];
    let mut basins: Vec<Basin> = Vec::new();
    let mut stack = Vec::new();
    for s in 0..nn {
        if !up[s] || seen[s] { continue; }
        seen[s] = true;
        stack.push(s);
        let mut b = Basin { floor_m: f64::INFINITY, ..Default::default() };
        let (mut sx, mut sy, mut in_box) = (0f64, 0f64, false);
        while let Some(i) = stack.pop() {
            b.posts += 1;
            let dz = z[i] - z0[i];
            b.volume_m3 += dz * dx * dx;
            b.max_raise_m = b.max_raise_m.max(dz);
            b.floor_m = b.floor_m.min(z0[i]);
            b.spill_m = b.spill_m.max(z[i]);
            b.ring_origin |= ring[i];
            let (gx, gy) = lat.g(i);
            sx += gx as f64 * dx; sy += gy as f64 * dx;
            if box_dist(i) == 0.0 { in_box = true; }
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            for (ox, oy) in D8 {
                let (xx, yy) = (x + ox, y + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let j = yy as usize * n + xx as usize;
                if up[j] && !seen[j] { seen[j] = true; stack.push(j); }
            }
        }
        if !in_box { continue; }
        b.area_m2 = b.posts as f64 * dx * dx;
        b.x_m = sx / b.posts as f64; b.y_m = sy / b.posts as f64;
        r.basins_in_box += 1; r.basins_in_box_posts += b.posts; r.basins_in_box_volume_m3 += b.volume_m3;
        if b.floor_m < sea { r.subsea_basins_in_box += 1; r.subsea_basins_in_box_area_m2 += b.area_m2; r.subsea_basins_in_box_volume_m3 += b.volume_m3; }
        basins.push(b);
    }
    basins.sort_by(|a, b| b.posts.cmp(&a.posts).then(a.x_m.total_cmp(&b.x_m)).then(a.y_m.total_cmp(&b.y_m)));
    basins.truncate(6);
    r.largest = basins;
    (raised, maxr, r)
}

// ------------------------------------------------------------------ T3 / T4: shore points, exposure, vector propagation

/// Exposure per post (non-zero only on shore sea posts).
pub fn exposure(lat: &EvalLat, sea: &dyn Fn(usize) -> bool, dirs: u32, cap_m: f64) -> Vec<f32> {
    let (n, nn, dx) = (lat.n, lat.nodes(), lat.dx());
    let mut e = vec![0f32; nn];
    let steps = (cap_m / dx).floor() as i64;
    let dv: Vec<(f64, f64)> = (0..dirs).map(|k| { let a = std::f64::consts::TAU * k as f64 / dirs as f64; (a.cos(), a.sin()) }).collect();
    for i in 0..nn {
        if !sea(i) { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let shore = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)].iter().any(|&(ox, oy)| {
            let (xx, yy) = (x + ox, y + oy);
            xx >= 0 && yy >= 0 && xx < n as i64 && yy < n as i64 && !sea(yy as usize * n + xx as usize)
        });
        if !shore { continue; }
        let mut sum = 0.0;
        for &(c, s) in &dv {
            let mut k = 1;
            while k <= steps {
                let (xx, yy) = ((x as f64 + c * k as f64).round() as i64, (y as f64 + s * k as f64).round() as i64);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { break; }
                if !sea(yy as usize * n + xx as usize) { break; }
                k += 1;
            }
            sum += ((k - 1) as f64 * dx).min(cap_m);
        }
        e[i] = (sum / (dirs as f64 / 2.0 * cap_m)).clamp(0.0, 1.0) as f32;
    }
    e
}

/// Shore segments (global post coordinates, the segment's exposure) per lattice cell, from the sea mask and the surface.
pub fn shore_segments(lat: &EvalLat, z: &[f64], sea_lvl: f64, sea: &dyn Fn(usize) -> bool, ex: Option<&[f32]>) -> Vec<(usize, [f64; 4], f32)> {
    let n = lat.n;
    let mut out = Vec::new();
    if n < 2 { return out; }
    for y in 0..n - 1 { for x in 0..n - 1 {
        let c = [y * n + x, y * n + x + 1, (y + 1) * n + x + 1, (y + 1) * n + x];
        let s = [sea(c[0]), sea(c[1]), sea(c[2]), sea(c[3])];
        let ns = s.iter().filter(|&&b| b).count();
        if ns == 0 || ns == 4 { continue; }
        let (gx, gy) = lat.g(c[0]);
        let pc = [(gx as f64, gy as f64), (gx as f64 + 1.0, gy as f64), (gx as f64 + 1.0, gy as f64 + 1.0), (gx as f64, gy as f64 + 1.0)];
        // Crossing on edge k (corner k -> corner k+1).
        let cross = |k: usize| -> ([f64; 2], f32) {
            let (a, b) = (k, (k + 1) % 4);
            let (za, zb) = (z[c[a]], z[c[b]]);
            let t = if (zb - za).abs() > 1e-9 { ((sea_lvl - za) / (zb - za)).clamp(0.0, 1.0) } else { 0.5 };
            let e = ex.map(|v| if s[a] { v[c[a]] } else { v[c[b]] }).unwrap_or(1.0);
            ([pc[a].0 + t * (pc[b].0 - pc[a].0), pc[a].1 + t * (pc[b].1 - pc[a].1)], e)
        };
        let ed: Vec<usize> = (0..4).filter(|&k| s[k] != s[(k + 1) % 4]).collect();
        let cell = y * n + x;
        let mut seg = |p: ([f64; 2], f32), q: ([f64; 2], f32)| out.push((cell, [p.0[0], p.0[1], q.0[0], q.0[1]], 0.5 * (p.1 + q.1)));
        if ed.len() == 2 { seg(cross(ed[0]), cross(ed[1])); }
        else if ed.len() == 4 {
            // Saddle: corners 0 and 2 share a class; the centre's class (4-corner mean against sea level) decides.
            let zc = 0.25 * (z[c[0]] + z[c[1]] + z[c[2]] + z[c[3]]);
            let centre_sea = zc <= sea_lvl;
            if centre_sea == s[0] { seg(cross(0), cross(1)); seg(cross(2), cross(3)); } else { seg(cross(3), cross(0)); seg(cross(1), cross(2)); }
        }
    } }
    out
}

#[inline]
fn seg_near(p: [f64; 2], s: &[f64; 4]) -> [f64; 2] {
    let (dx, dy) = (s[2] - s[0], s[3] - s[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 { (((p[0] - s[0]) * dx + (p[1] - s[1]) * dy) / l2).clamp(0.0, 1.0) } else { 0.0 };
    [s[0] + t * dx, s[1] + t * dy]
}

/// Result of a vector propagation: distance (m, INFINITY beyond the cap or unreached), the carried exposure, the post the
/// distance was carried from (u32::MAX for initial posts) and the settle order.
pub struct VProp { pub d: Vec<f64>, pub e: Vec<f32>, pub parent: Vec<u32>, pub order: Vec<u32> }

#[derive(Clone, Copy)]
struct DK { d: f64, gy: i64, gx: i64, i: u32 }
impl PartialEq for DK { fn eq(&self, o: &Self) -> bool { self.cmp(o) == Ordering::Equal } }
impl Eq for DK {}
impl PartialOrd for DK { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for DK { fn cmp(&self, o: &Self) -> Ordering { o.d.total_cmp(&self.d).then(o.gy.cmp(&self.gy)).then(o.gx.cmp(&self.gx)) } }

/// Vector propagation of the nearest shore point through posts with `pass` (heap key (d, global gy, gx)). Initial posts:
/// passable posts within 2 posts of a segment get the exact nearest point of their nearby segments.
pub fn vprop_segments(lat: &EvalLat, segs: &[(usize, [f64; 4], f32)], pass: &dyn Fn(usize) -> bool, cap_m: f64) -> VProp {
    let (n, nn, dx) = (lat.n, lat.nodes(), lat.dx());
    let mut d = vec![f64::INFINITY; nn];
    let mut e = vec![0f32; nn];
    let mut q = vec![[0f64; 2]; nn];
    let mut sid = vec![u32::MAX; nn];
    let mut parent = vec![u32::MAX; nn];
    let gpos = |i: usize| { let (gx, gy) = lat.g(i); [gx as f64, gy as f64] };
    // Segments by cell (segs are in cell order): start offsets.
    let mut start = vec![0u32; nn + 1];
    for (cell, _, _) in segs { start[*cell + 1] += 1; }
    for c in 0..nn { start[c + 1] += start[c]; }
    for (k, (cell, s, se)) in segs.iter().enumerate() {
        let (cx, cy) = ((cell % n) as i64, (cell / n) as i64);
        for oy in -2..=3 { for ox in -2..=3 {
            let (xx, yy) = (cx + ox, cy + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if !pass(j) { continue; }
            let p = gpos(j);
            let a = seg_near(p, s);
            let dd = (p[0] - a[0]).hypot(p[1] - a[1]) * dx;
            if dd < d[j] || (dd == d[j] && (a[1], a[0]) < (q[j][1], q[j][0])) { d[j] = dd; q[j] = a; e[j] = *se; sid[j] = k as u32; }
        } }
    }
    let mut heap = BinaryHeap::new();
    for i in 0..nn { if d[i] <= cap_m { let (gx, gy) = lat.g(i); heap.push(DK { d: d[i], gy, gx, i: i as u32 }); } else { d[i] = f64::INFINITY; } }
    let mut done = vec![false; nn];
    let mut order = Vec::new();
    while let Some(DK { d: v, i, .. }) = heap.pop() {
        let i = i as usize;
        if done[i] || v > d[i] { continue; }
        done[i] = true;
        order.push(i as u32);
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if done[j] || !pass(j) { continue; }
            let p = gpos(j);
            // The neighbour's shore segment and the segments of the 3 x 3 cells around it: exact nearest point
            // (vector propagation by segment).
            let s0 = segs[sid[i] as usize].0;
            let (sx, sy) = ((s0 % n) as i64, (s0 / n) as i64);
            let (mut c, mut a, mut kb) = (f64::INFINITY, [0f64; 2], u32::MAX);
            for oy in -1..=1 { for ox in -1..=1 {
                let (xx, yy) = (sx + ox, sy + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let cc = yy as usize * n + xx as usize;
                for k in start[cc]..start[cc + 1] {
                    let b = seg_near(p, &segs[k as usize].1);
                    let dd = (p[0] - b[0]).hypot(p[1] - b[1]) * dx;
                    if dd < c { c = dd; a = b; kb = k; }
                }
            } }
            if c > cap_m { continue; }
            if c < d[j] - 1e-9 {
                d[j] = c; q[j] = a; e[j] = segs[kb as usize].2; sid[j] = kb; parent[j] = i as u32;
                let (gx, gy) = lat.g(j);
                heap.push(DK { d: c, gy, gx, i: j as u32 });
            }
        }
    }
    for i in 0..nn { if !done[i] { d[i] = f64::INFINITY; } }
    VProp { d, e, parent, order }
}

/// Segments from sea POST centres (each shore sea post a zero-length segment; exposure of that post): the v1 distance
/// definition carried by the same propagation (used by soft_bluff_v2 without the sub-post shore).
pub fn post_segments(lat: &EvalLat, sea: &dyn Fn(usize) -> bool, ex: Option<&[f32]>) -> Vec<(usize, [f64; 4], f32)> {
    let n = lat.n;
    let mut out = Vec::new();
    for i in 0..lat.nodes() {
        if !sea(i) { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let shore = D8.iter().any(|&(ox, oy)| { let (xx, yy) = (x + ox, y + oy); xx >= 0 && yy >= 0 && xx < n as i64 && yy < n as i64 && !sea(yy as usize * n + xx as usize) });
        if !shore { continue; }
        let (gx, gy) = lat.g(i);
        out.push((i, [gx as f64, gy as f64, gx as f64, gy as f64], ex.map(|v| v[i]).unwrap_or(1.0)));
    }
    out
}

/// Top of the first hard layer (k_mult <= kmax) in the post's column between lo and hi (hi soft), if any.
pub fn hard_top(k_at: &dyn Fn(f64) -> f64, lo: f64, hi: f64, kmax: f64) -> Option<f64> {
    if !(hi > lo) { return None; }
    let mut h = hi;
    let step = 2.0;
    loop {
        let nh = (h - step).max(lo);
        if k_at(nh) <= kmax {
            // Contact between nh (hard) and h (soft): bisect.
            let (mut a, mut b) = (nh, h);
            for _ in 0..12 { let m = 0.5 * (a + b); if k_at(m) <= kmax { a = m; } else { b = m; } }
            return Some(a);
        }
        if nh <= lo { return None; }
        h = nh;
    }
}

/// Lane TF coast statistics (window header, only with tf).
#[derive(Clone, Debug, Default, Serialize)]
pub struct TfCoastStats {
    pub shore_segments: u32,
    pub exposure_shore_posts: u32,
    pub exposure_mean: f64,
    pub bluff_v2_cut_posts: u32,
    pub bluff_v2_bench_posts: u32,
    pub bluff_v2_chain_posts: u32,
    pub bluff_v2_cut_max_m: f64,
    pub cut_posts: u32,
    pub beach_posts: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub final_fill: Option<FillReport>,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fill_report_matches_fill_and_flags_ring() {
        let lat = EvalLat { gx0: 0, gy0: 0, n: 7, dx_mm: 32_000 };
        let mut z = vec![10.0; 49];
        z[24] = 1.0; // interior pit
        z[3 * 7] = -5.0; z[3 * 7 + 1] = -5.0; // a low reaching the ring (x 0)
        let outlet = vec![false; 49];
        let mut a = z.clone();
        let (r0, m0) = crate::t1b_rivernet::fill(&lat, &mut a, &outlet);
        let mut b = z.clone();
        let (r1, m1, rep) = fill_report(&lat, &mut b, &outlet, 0.0, &|_| true, 1, 8000.0);
        assert_eq!((r0, m0), (r1, m1));
        assert_eq!(a, b);
        assert_eq!(rep.basins_in_box, 1);
        assert!(rep.ring_origin_subsea_posts >= 2, "sub-sea ring-drained posts flagged");
        assert!(!rep.guard_pass);
    }
    #[test]
    fn shore_distance_is_smooth_along_an_oblique_coast() {
        // A plane rising at slope 1 per post along (2, 1): sea where z <= 0. Post-centre distance alternates; the
        // contour distance is the exact plane distance.
        let n = 40;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 };
        let nr = (5f64).sqrt();
        let z: Vec<f64> = (0..n * n).map(|i| { let (x, y) = ((i % n) as f64, (i / n) as f64); (2.0 * x + y - 40.0) / nr }).collect();
        let sea = |i: usize| z[i] <= 0.0;
        let segs = shore_segments(&lat, &z, 0.0, &sea, None);
        let v = vprop_segments(&lat, &segs, &|i| !sea(i), 1e9);
        for i in 0..n * n {
            let (x, y) = (i % n, i / n);
            if sea(i) || x < 3 || y < 3 || x > n - 4 || y > n - 4 { continue; }
            // Only posts whose perpendicular foot on the line lies inside the lattice.
            let t = z[i] / nr;
            let (fx, fy) = (x as f64 - 2.0 * t, y as f64 - t);
            if fx < 1.0 || fy < 1.0 || fx > (n - 2) as f64 || fy > (n - 2) as f64 { continue; }
            let exact = z[i] * 32.0;
            assert!((v.d[i] - exact).abs() < 0.05, "post {x},{y}: {} vs {exact}", v.d[i]);
        }
    }
    #[test]
    fn exposure_open_vs_inlet() {
        let n = 60;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 };
        // Sea for x < 30, plus a 3-post-wide inlet along y 28..30 into the land up to x 50.
        let sea = |i: usize| { let (x, y) = (i % n, i / n); x < 30 || ((28..=30).contains(&y) && x < 50) };
        let e = exposure(&lat, &sea, 16, 640.0);
        let open = e[10 * n + 29] as f64;
        let inlet = e[28 * n + 45] as f64;
        assert!(open > 0.8, "open {open}");
        assert!(inlet < 0.5, "inlet {inlet}");
    }
    #[test]
    fn hard_top_finds_contact() {
        let k = |h: f64| if h < 123.4 { 0.6 } else { 2.0 };
        let t = hard_top(&k, 50.0, 300.0, 0.8).unwrap();
        assert!((t - 123.4).abs() < 0.01, "{t}");
        assert!(hard_top(&k, 130.0, 300.0, 0.8).is_none());
    }
}

/// Per-lattice inputs of the lane TF coast steps: the C1 T0 surface per post (T5; empty when off) and the published-box
/// margin in posts (T1 report; 0 = the whole lattice, as in the chart-wide E solve).
pub struct TfIn<'a> { pub bil_c1: &'a [f32], pub margin_posts: i64 }
