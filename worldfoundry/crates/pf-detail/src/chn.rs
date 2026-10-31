//! Lane INTS2n4 (2026-10-21, EXPLORATORY; constants are project assumptions awaiting review): channel-network fixes at 2 m.
//!
//! `arm_cl=1` (needs version=2 or 3): T1b's marine arm posts (its estuary rule fixes a T1a channel post whose bed is at or
//! below sea level as marine, so a drowned trunk is a one-post-wide line of sea posts) are traced and carved as a channel on a
//! smooth centreline (the v2 curves) instead of being the nearest-T1b-post staircase of fixed sea. `ringflow_ch=1` (needs
//! de_ringflow=1): the T1b drainage area that enters a B5 domain across its ring excludes T1b channel / marine posts.

use crate::input::T1bWindow;
use serde::Serialize;

pub const ARM_CL_RULE: &str = "arm_cl_v1 (lane INTS2n4): marine arm posts = T1b fixed-2 posts outside the 3x3 morphological opening of the marine mask (fixed 2 or 4, window edge clamped), in 8-connected components of >= 4 posts that touch a T1b channel post (fixed 1 or 5): the drowned T1a trunk of T1b's estuary rule. Under the key they are channel posts of the v2 curves: traced toward the open sea (receiver = the arm / water neighbour with the lowest (height, steps to open water) key), smoothed with smooth_passes binomial passes (a stream that contains arm posts) then the usual Chaikin passes; width from the B6 hydraulic-geometry law on the T1a discharge near the post (channel_q_near, else the T1b area x runoff), widened seaward as a funnel-shaped tidal channel, half-width hw(u) = min(hw_r exp(u / b), hw_max) with u = distance along the smoothed line from the arm head (first arm post of the stream), depth from the law unchanged (estuary depth roughly constant along the channel); bed = the T1b height (sea level) so the carved interior lies below sea level. The arm posts are de-slotted in the valley lattice (zmod = mean of the free 8-neighbours, as T1b channel posts), a 2 m post is fixed sea only if it is outside every arm line's half-width and the de-slotted surface is within fixed_above_sea_m of sea level, or its nearest T1b post is marine and not an arm post; arm lines do not reproject the valley (V2: the arm's valley is T1b's floor), so the B6 carve holds the whole arm cross-section; with b8_levee the arm carve and banks (posts within hw + max(bank_width_m, 48 m) of an arm line) are booked from the B1 base (the T1b surface, which already holds the 32 m sea trench), so B8 restores every other stage, including the de-slot fill on the free posts beside the arm, and leaves only the arm carve beyond T1b's trench in place, and the bed floor (bed_floor_radius_m) counts arm posts. The arm builds no b8_levee levee (a segment with both ends on arm posts is skipped): T1b holds the arm post as sea water with no bed (T1b publishes no bathymetry), so the carve below the T1b floor there has no T1b mass to return; it stays in the tile stats' carve_mean_m (B8 restores every other stage as before), so the deficit is reportable. B1 base (T1b bicubic) and so the B8 tile-mean target are unchanged: the arm's 32 m trench is in the base as a T1b channel slot is";
pub const RINGFLOW_CH_RULE: &str = "ringflow_ch_v1 (lane INTS2n4): de_ringflow=1 places T1b's upslope drainage area on the B5 domain ring as the inward flux of the bilinear T1b area; with this key the bilinear runs over the T1b posts that are not channel or marine (fixed 0; weights renormalised over those corners, 0 if none), because a channel post's area is the whole catchment carried by the B6 carved channel (a fixed corridor), not hillslope flow crossing the ring. Measured cause (preview-client west fan, tile 38_29): a ring-1 post 26 m from a T1b guide post took ~0.2 of the guide's area and B5 cut a 1.1 m deep, 160 m straight gully parallel to the guide; de_ringflow=0 removes it";

#[derive(Clone, Debug, Serialize)]
pub struct ArmParams {
    pub rule: String,
    /// Binomial smoothing passes of a stream containing arm posts (each pass: variance 1 post^2).
    pub smooth_passes: u32,
    /// Funnel convergence length b (m) of the tidal width.
    pub b_m: f64,
    /// Cap of the tidal half-width (m).
    pub hw_max_m: f64,
}

impl ArmParams {
    pub fn new() -> Self { ArmParams { rule: ARM_CL_RULE.into(), smooth_passes: 8, b_m: 4000.0, hw_max_m: 150.0 } }
}

#[derive(Clone, Debug, Serialize)]
pub struct ChnParams {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub arm: Option<ArmParams>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ringflow_ch: Option<String>,
}

#[derive(Default, Debug, Serialize, Clone)]
pub struct ArmStats {
    pub arm_posts: usize,
    pub components: usize,
    pub components_kept: usize,
    pub posts_without_open_water_path: usize,
    pub arm_streams: usize,
    pub hw_river_m_min_max: [f64; 2],
    pub hw_tidal_m_max: f64,
    pub arm_length_m_max: f64,
    /// Per arm stream: head x, y (chart m), discharge at the head (m3/s), river half-width hw_r (m), arm length (m), mouth hw (m).
    pub streams: Vec<[f64; 6]>,
}

/// The marine arm mask of a T1b window (ARM_CL_RULE), window-indexed, and its steps to open water (u32::MAX if none).
pub fn arm_mask(w: &T1bWindow, st: &mut ArmStats) -> (Vec<bool>, Vec<u32>) { arm_mask_codes(w.n, &w.fixed, st) }

/// arm_mask on a bare n x n fixed-code lattice (row-major).
pub fn arm_mask_codes(nu: usize, fixed: &[u8], st: &mut ArmStats) -> (Vec<bool>, Vec<u32>) {
    struct W<'a> { fixed: &'a [u8] }
    let w = W { fixed };
    let n = nu as i64;
    let at = |x: i64, y: i64| -> usize { (y.clamp(0, n - 1) * n + x.clamp(0, n - 1)) as usize };
    let marine: Vec<bool> = w.fixed.iter().map(|&c| c == 2 || c == 4).collect();
    let nn = (n * n) as usize;
    let mut ero = vec![false; nn];
    for y in 0..n { for x in 0..n {
        let mut all = true;
        for dy in -1..=1 { for dx in -1..=1 { if !marine[at(x + dx, y + dy)] { all = false; } } }
        ero[(y * n + x) as usize] = all;
    } }
    let mut cand = vec![false; nn];
    for y in 0..n { for x in 0..n {
        let i = (y * n + x) as usize;
        if w.fixed[i] != 2 { continue; }
        let mut any = false;
        for dy in -1..=1 { for dx in -1..=1 { if ero[at(x + dx, y + dy)] { any = true; } } }
        cand[i] = !any;
    } }
    // 8-connected components, kept if >= 4 posts and touching a T1b channel post
    let mut arm = vec![false; nn];
    let mut seen = vec![false; nn];
    for s in 0..nn {
        if !cand[s] || seen[s] { continue; }
        st.components += 1;
        let mut comp = vec![s];
        seen[s] = true;
        let mut k = 0;
        let mut touch = false;
        while k < comp.len() {
            let i = comp[k]; k += 1;
            let (x, y) = ((i as i64) % n, (i as i64) / n);
            for dy in -1..=1 { for dx in -1..=1 {
                let (a, b) = (x + dx, y + dy);
                if a < 0 || b < 0 || a >= n || b >= n { continue; }
                let j = (b * n + a) as usize;
                if w.fixed[j] == 1 || w.fixed[j] == 5 { touch = true; }
                if cand[j] && !seen[j] { seen[j] = true; comp.push(j); }
            } }
        }
        if comp.len() >= 4 && touch {
            st.components_kept += 1;
            for &i in &comp { arm[i] = true; }
        }
    }
    st.arm_posts = arm.iter().filter(|&&a| a).count();
    // steps to open water (marine, not arm) over arm posts, 8-connected BFS
    let mut dist = vec![u32::MAX; nn];
    let mut q = std::collections::VecDeque::new();
    for i in 0..nn {
        if !arm[i] { continue; }
        let (x, y) = ((i as i64) % n, (i as i64) / n);
        let mut open = false;
        for dy in -1..=1 { for dx in -1..=1 {
            let (a, b) = (x + dx, y + dy);
            if a < 0 || b < 0 || a >= n || b >= n { continue; }
            let j = (b * n + a) as usize;
            if marine[j] && !arm[j] { open = true; }
        } }
        if open { dist[i] = 0; q.push_back(i); }
    }
    while let Some(i) = q.pop_front() {
        let (x, y) = ((i as i64) % n, (i as i64) / n);
        for dy in -1..=1 { for dx in -1..=1 {
            let (a, b) = (x + dx, y + dy);
            if a < 0 || b < 0 || a >= n || b >= n { continue; }
            let j = (b * n + a) as usize;
            if arm[j] && dist[j] == u32::MAX { dist[j] = dist[i] + 1; q.push_back(j); }
        } }
    }
    st.posts_without_open_water_path = (0..nn).filter(|&i| arm[i] && dist[i] == u32::MAX).count();
    (arm, dist)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn arm_keys_are_opt_in() {
        let mut p = crate::DetailParamsV1::v3();
        let before = serde_json::to_string(&p).unwrap();
        assert!(!before.contains("chn"));
        assert!(crate::DetailParamsV1::v1().set("arm_cl", "1").is_err(), "needs version 2/3");
        assert!(p.set("arm_cl_b_m", "100").is_err(), "needs arm_cl first");
        assert!(p.set("ringflow_ch", "1").is_err(), "needs de_ringflow=1");
        p.set("arm_cl", "1").unwrap();
        p.set("arm_cl_b_m", "3000").unwrap();
        assert_eq!(p.chn.as_ref().unwrap().arm.as_ref().unwrap().b_m, 3000.0);
        p.set("arm_cl", "0").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), before);
    }
    #[test]
    fn arm_mask_keeps_the_drowned_trunk_only() {
        // 16 x 16: a 6-wide bay on the west (x 0..5), an arm from the bay east along a staircase to a T1a post (code 1), and a
        // one-post coastal strip (y = 14) off the bay that touches no channel
        let n = 16usize;
        let mut f = vec![0u8; n * n];
        for y in 0..n { for x in 0..6 { f[y * n + x] = 2; } }
        let arm = [(6, 5), (7, 5), (8, 6), (9, 6), (10, 7), (11, 7), (12, 8)];
        for &(x, y) in &arm { f[y * n + x] = 2; }
        f[8 * n + 13] = 1;
        for x in 6..9 { f[14 * n + x] = 2; }
        let mut st = ArmStats::default();
        let (m, d) = arm_mask_codes(n, &f, &mut st);
        for &(x, y) in &arm { assert!(m[y * n + x], "arm post {x},{y}"); }
        for x in 6..9 { assert!(!m[14 * n + x], "strip post {x} is not an arm"); }
        for y in 0..n { for x in 0..6 { assert!(!m[y * n + x], "bay post {x},{y}"); } }
        assert_eq!(st.components_kept, 1);
        assert_eq!(d[5 * n + 6], 0);
        assert_eq!(d[8 * n + 12], 6);
    }
}
