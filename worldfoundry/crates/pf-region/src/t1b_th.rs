//! T1b 0.9.4 lane TH (2026-10-15; `unpublished-design-record`): lane TG ideas 1 and 2 on walk config tG2.
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon; construction (no SimTime). Everything is opt-in under
//! `params.ls.tb.tc.th` (inside params.ls.tb, so lane TB's paired keying keeps every keyed draw); absent = lane TG's
//! behaviour byte-identical.
//!
//! - H1 alongshore exposure smoothing: lane TF's fetch exposure is computed per shore sea post from 16 rays on the post
//!   lattice, so it changes from post to post along a coast; the cliff platform width (cliff_retreat x k x e) and the beach
//!   width inherit that noise (lane TG: a dashed cliff wall at 60-80 m). H1 averages e along the coast (geodesic through shore
//!   sea posts, triangular weights, declared length) before it is attached to the shore segments.
//! - H2 `uncut_mm`: an additive window field with the height each cliff-class post had immediately before the wave cut, so a
//!   consumer can evaluate T1b's own cliff rule min(U, P(s)) at its own spacing instead of estimating U.

use crate::t1b_rivernet::EvalLat;
use crate::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ThParamsV1 {
    /// H1: alongshore smoothing of the fetch exposure (flag; needs tf exposure) and its declared length L (mm, the full
    /// support of the triangular kernel: half-width L / 2).
    pub expo_smooth: u32,
    pub expo_len_mm: i64,
    /// H2: publish uncut_mm (flag; needs the CF cliff profile).
    pub uncut: u32,
    pub expo_rule: String,
    pub uncut_rule: String,
}

impl ThParamsV1 {
    pub fn v1() -> Self {
        ThParamsV1 {
            expo_smooth: 1, expo_len_mm: 384_000, uncut: 1,
            expo_rule: "expo_smooth_v1: before the exposure e is attached to the shore segments, every shore sea post (a sea post with a non-sea 4-neighbour, as fetch_exposure_v1) takes e_s = sum w_j e_j / sum w_j over the shore sea posts j reachable from it THROUGH shore sea posts (8-neighbour links of length dx and sqrt 2 dx; Dijkstra geodesic g_j, heap key (g, global gy, gx); sums in global (gy, gx) order) with g_j < h, w_j = 1 - g_j / h, h = expo_len / 2. Declared expo_len 384 m: the lattice-scale noise (2-3 posts) is removed (triangle response sinc^2(h / lambda): 0 at 96 and 64 m, 0.016 at 80 m) while bay / headland contrast at >= 1 km keeps 88%; recalled, unverified: refraction and diffraction spread swell energy alongshore over about 1-3 wavelengths (156 m for a 10 s swell, deep-water L = g T^2 / 2 pi). Along the coast, not Euclidean: the two sides of a headland do not mix. e_s replaces e in the cliff platform width, the beach width and soft_bluff_v2's beach width; lane TF's exposure statistics stay those of the raw e".into(),
            uncut_rule: "uncut_v1 (additive window field uncut_mm, heights unchanged): on cliff-class posts (the coast stage evaluated the CF C8 hard-rock cliff profile there: k <= cliff_k_max, cliff_profile on, finite shore distance > 0, not marine, lake, channel, T1a or fan no-cut; includes the posts the cut sends to sea) the post's height immediately before the wave cut (after edifices and fans); i32::MIN elsewhere. With shore_mm (s) and the CF constants a consumer evaluates the rule at any spacing: P(s) = talus_slope s for 0 <= s < talus_width, talus_slope talus_width + wall_slope (s - talus_width) beyond; z = min(U, P(s)) with U and s bilinear (published surface outside the class, blended by the bilinear class flag). Exact at a post only where the platform width varies slowly (s = d - pw); measured, not assumed".into(),
        }
    }
    pub fn off() -> Self { let mut p = Self::v1(); p.expo_smooth = 0; p.uncut = 0; p }
    pub fn validate(&self) -> R<()> {
        if self.expo_len_mm <= 0 { return Err("th: expo_len_mm > 0".into()); }
        Ok(())
    }
}

pub fn th(p: &crate::t1b_rivernet::RiverNetParamsV1) -> Option<&ThParamsV1> { crate::t1b_ls::tc(p).and_then(|t| t.th.as_ref()) }
pub fn th_on(p: &crate::t1b_rivernet::RiverNetParamsV1, f: fn(&ThParamsV1) -> u32) -> bool { th(p).map(|t| f(t) == 1).unwrap_or(false) }

/// Sentinel of uncut_mm on posts outside the cliff class.
pub const UNCUT_NONE: i32 = i32::MIN;

#[derive(Clone, Copy)]
struct GK { g: f64, gy: i64, gx: i64, i: u32 }
impl PartialEq for GK { fn eq(&self, o: &Self) -> bool { self.cmp(o) == Ordering::Equal } }
impl Eq for GK {}
impl PartialOrd for GK { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for GK { fn cmp(&self, o: &Self) -> Ordering { o.g.total_cmp(&self.g).then(o.gy.cmp(&self.gy)).then(o.gx.cmp(&self.gx)) } }

/// Shore sea posts of a mask (a sea post with a non-sea 4-neighbour inside the lattice), as `t1b_tf::exposure` defines them.
pub fn shore_posts(lat: &EvalLat, sea: &dyn Fn(usize) -> bool) -> Vec<bool> {
    let n = lat.n;
    (0..lat.nodes()).map(|i| {
        if !sea(i) { return false; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)].iter().any(|&(ox, oy)| {
            let (xx, yy) = (x + ox, y + oy);
            xx >= 0 && yy >= 0 && xx < n as i64 && yy < n as i64 && !sea(yy as usize * n + xx as usize)
        })
    }).collect()
}

/// H1 statistics (window header, only with th expo_smooth).
#[derive(Clone, Debug, Default, Serialize)]
pub struct ThExpoStats {
    pub shore_posts: u32,
    /// Mean |e_i - e_j| over 8-adjacent shore post pairs, raw and smoothed.
    pub adj_abs_diff_raw: f64,
    pub adj_abs_diff_smoothed: f64,
    pub mean_abs_change: f64,
    pub max_abs_change: f64,
}

/// H1: alongshore triangular smoothing of a per-post exposure over the shore sea posts (see `expo_rule`). Posts that are not
/// shore posts keep their value.
pub fn smooth_exposure(lat: &EvalLat, e: &[f32], shore: &[bool], len_m: f64) -> (Vec<f32>, ThExpoStats) {
    let (n, nn, dx) = (lat.n, lat.nodes(), lat.dx());
    let h = 0.5 * len_m;
    let mut out = e.to_vec();
    let mut g = vec![f64::INFINITY; nn];
    let mut done = vec![false; nn];
    let mut touched: Vec<usize> = Vec::new();
    let mut ball: Vec<(i64, i64, usize, f64)> = Vec::new();
    let mut heap = BinaryHeap::new();
    let mut st = ThExpoStats::default();
    let (mut ch_sum, mut ch_max) = (0f64, 0f64);
    for s in 0..nn {
        if !shore[s] { continue; }
        st.shore_posts += 1;
        ball.clear();
        let (sgx, sgy) = lat.g(s);
        g[s] = 0.0; touched.push(s);
        heap.push(GK { g: 0.0, gy: sgy, gx: sgx, i: s as u32 });
        while let Some(GK { g: v, gy, gx, i }) = heap.pop() {
            let i = i as usize;
            if done[i] || v > g[i] { continue; }
            done[i] = true;
            ball.push((gy, gx, i, v));
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            for (ox, oy) in D8 {
                let (xx, yy) = (x + ox, y + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let j = yy as usize * n + xx as usize;
                if !shore[j] || done[j] { continue; }
                let w = if ox != 0 && oy != 0 { dx * std::f64::consts::SQRT_2 } else { dx };
                let c = v + w;
                if c >= h || c >= g[j] { continue; }
                if g[j].is_infinite() { touched.push(j); }
                g[j] = c;
                let (jgx, jgy) = lat.g(j);
                heap.push(GK { g: c, gy: jgy, gx: jgx, i: j as u32 });
            }
        }
        ball.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        let (mut sw, mut se) = (0f64, 0f64);
        for &(_, _, j, gj) in &ball { let w = 1.0 - gj / h; sw += w; se += w * e[j] as f64; }
        out[s] = (se / sw) as f32;
        let d = (out[s] as f64 - e[s] as f64).abs();
        ch_sum += d; ch_max = ch_max.max(d);
        for &k in &touched { g[k] = f64::INFINITY; done[k] = false; }
        touched.clear();
    }
    // Statistics: adjacent shore pairs (each unordered pair once: forward half of D8).
    let (mut a_raw, mut a_sm, mut np) = (0f64, 0f64, 0u64);
    for i in 0..nn {
        if !shore[i] { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in [(1i64, 0i64), (1, 1), (0, 1), (-1, 1)] {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if !shore[j] { continue; }
            a_raw += (e[i] as f64 - e[j] as f64).abs(); a_sm += (out[i] as f64 - out[j] as f64).abs(); np += 1;
        }
    }
    if np > 0 { st.adj_abs_diff_raw = a_raw / np as f64; st.adj_abs_diff_smoothed = a_sm / np as f64; }
    if st.shore_posts > 0 { st.mean_abs_change = ch_sum / st.shore_posts as f64; }
    st.max_abs_change = ch_max;
    (out, st)
}

/// Lane TH coast statistics (window header, only with th).
#[derive(Clone, Debug, Default, Serialize)]
pub struct ThCoastStats {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expo: Option<ThExpoStats>,
    /// H2: cliff-class posts, of which the cut lowered by more than 0.5 m (kept land) and sent to sea.
    pub uncut_class_posts: u32,
    pub uncut_lowered_posts: u32,
    pub uncut_to_sea_posts: u32,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn smoothing_keeps_constants_and_removes_alternation() {
        // A straight coast along y (sea x < 20): shore posts at x 19. Alternating e on the shore column.
        let n = 60;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 };
        let sea = |i: usize| (i % n) < 20;
        let shore = shore_posts(&lat, &sea);
        assert_eq!(shore.iter().filter(|&&b| b).count(), n);
        let flat: Vec<f32> = (0..n * n).map(|i| if shore[i] { 0.7 } else { 0.0 }).collect();
        let (s, _) = smooth_exposure(&lat, &flat, &shore, 384.0);
        for i in 0..n * n { if shore[i] { assert!((s[i] - 0.7).abs() < 1e-6); } }
        let alt: Vec<f32> = (0..n * n).map(|i| if shore[i] { if (i / n) % 2 == 0 { 0.4 } else { 1.0 } } else { 0.0 }).collect();
        let (s, st) = smooth_exposure(&lat, &alt, &shore, 384.0);
        assert!((st.adj_abs_diff_raw - 0.6).abs() < 1e-6);
        assert!(st.adj_abs_diff_smoothed < 0.05 * st.adj_abs_diff_raw, "{}", st.adj_abs_diff_smoothed);
        for y in 10..50 { assert!((s[y * n + 19] - 0.7).abs() < 0.03, "{}", s[y * n + 19]); }
    }
    #[test]
    fn smoothing_is_lattice_independent_and_follows_the_coast() {
        // Two lattices offset by (7, 5) posts over the same global mask: identical smoothed values at shared shore posts. A thin
        // peninsula: its two sides are 3 posts apart in Euclidean distance but far along the coast, so they do not mix.
        let mask = |gx: i64, gy: i64| -> bool { !(gx >= 20 && (gy < 30 || (gx <= 23 && gy < 70))) };
        let ev = |gx: i64, gy: i64| -> f32 { let b: f32 = if gx < 22 { 0.2 } else { 0.9 }; b + 0.01 * ((gx * 7 + gy * 13) % 5) as f32 };
        let mk = |gx0: i64, gy0: i64, n: usize| {
            let lat = EvalLat { gx0, gy0, n, dx_mm: 32_000 };
            let sea = |i: usize| { let (gx, gy) = lat.g(i); mask(gx, gy) };
            let shore = shore_posts(&lat, &sea);
            let e: Vec<f32> = (0..n * n).map(|i| { let (gx, gy) = lat.g(i); if shore[i] { ev(gx, gy) } else { 0.0 } }).collect();
            let (s, _) = smooth_exposure(&lat, &e, &shore, 384.0);
            (lat, shore, s)
        };
        let (la, sa, va) = mk(0, 0, 100);
        let (lb, sb, vb) = mk(7, 5, 100);
        let mut shared = 0;
        for i in 0..la.nodes() {
            let (gx, gy) = la.g(i);
            let Some(j) = lb.idx(gx, gy) else { continue };
            if !sa[i] || !sb[j] { continue; }
            // Away from both lattices' edges (the 6-post ball fits).
            if gx < 7 + 7 || gy < 5 + 7 || gx > 93 || gy > 93 { continue; }
            assert_eq!(va[i].to_bits(), vb[j].to_bits(), "{gx} {gy}");
            shared += 1;
        }
        assert!(shared > 50);
        // West side of the peninsula (x 19, low e) vs east side (x 24, high e) at y 50: no mixing across.
        let w = va[la.idx(19, 50).unwrap()]; let e = va[la.idx(24, 50).unwrap()];
        assert!(w < 0.3 && e > 0.8, "{w} {e}");
    }
}
