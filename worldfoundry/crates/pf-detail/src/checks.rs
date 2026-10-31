//! Per-tile checks: F6 tile mean, slope histogram (G2), F10 walkability share,
//! F5 pits, F1 channel descent, F8c invention magnitude.

use crate::tile::TileOut;
use crate::{POST_M, TILE_POSTS};
use serde::Serialize;

#[derive(Serialize, Debug, Clone)]
pub struct TileChecks {
    pub key: String,
    pub f6_t2_mean_m: f64,
    pub f6_base_mean_m: f64,
    pub f6_diff_m: f64,
    pub f6_pass: bool,
    /// T2 mean minus the plain T1b bicubic mean (differs from f6_diff_m only with deperch_on).
    pub f6_raw_t1b_diff_m: f64,
    pub slope_hist_1deg: Vec<f64>,
    pub slope_max_bin_share: f64,
    pub slope_max_bin_deg: usize,
    pub slope_peaks_deg: Vec<usize>,
    pub slope_bimodal: bool,
    pub slope_mean_deg: f64,
    pub base_slope_mean_deg: f64,
    pub base_slope_max_bin_share: f64,
    pub f10_share_gt_600permille: f64,
    pub f10_base_share_gt_600permille: f64,
    pub f10_points_over_base: f64,
    pub f5_interior_pits: u64,
    pub f1_t1a_interior_posts: u64,
    pub f1_paths_leaving_corridor: u64,
    pub f1_paths_ending_in_pit: u64,
    pub f1_pass: bool,
    pub f8c_abs_diff_p99_m: f64,
    pub f8c_abs_diff_max_m: f64,
    pub f8c_bound_m: f64,
    pub f8c_pass: bool,
    pub litho_share: Vec<(u8, f64)>,
    pub layers_exposed: usize,
}

fn slopes_deg(h: &[i32]) -> Vec<f64> {
    let n = TILE_POSTS;
    let mut out = Vec::with_capacity((n - 2) * (n - 2));
    for j in 1..n - 1 { for i in 1..n - 1 {
        let z = |a: usize, b: usize| h[b * n + a] as f64 / 1000.0;
        let gx = (z(i + 1, j) - z(i - 1, j)) / (2.0 * POST_M);
        let gy = (z(i, j + 1) - z(i, j - 1)) / (2.0 * POST_M);
        out.push(gx.hypot(gy).atan().to_degrees());
    } }
    out
}

fn hist(s: &[f64]) -> Vec<f64> {
    let mut h = vec![0.0f64; 90];
    for &v in s { h[(v.floor() as usize).min(89)] += 1.0; }
    let n = s.len() as f64;
    h.iter().map(|c| c / n).collect()
}

/// G2 bimodality: 3-bin running mean has two local maxima >= 10 deg apart with a minimum between <= 0.8 x the lower peak.
pub fn bimodal(h: &[f64]) -> (bool, Vec<usize>) {
    let m: Vec<f64> = (0..h.len()).map(|i| {
        let a = i.saturating_sub(1); let b = (i + 1).min(h.len() - 1);
        h[a..=b].iter().sum::<f64>() / (b - a + 1) as f64
    }).collect();
    let mut peaks = Vec::new();
    for i in 0..m.len() {
        let l = if i == 0 { 0.0 } else { m[i - 1] };
        let r = if i + 1 == m.len() { 0.0 } else { m[i + 1] };
        if m[i] > 0.002 && m[i] >= l && m[i] > r { peaks.push(i); }
    }
    for a in 0..peaks.len() { for b in a + 1..peaks.len() {
        let (p, q) = (peaks[a], peaks[b]);
        if q - p < 10 { continue; }
        let lo = m[p..=q].iter().cloned().fold(f64::MAX, f64::min);
        if lo <= 0.8 * m[p].min(m[q]) { return (true, vec![p, q]); }
    } }
    (false, peaks)
}

pub fn check_tile(t: &TileOut, f8_bound_m: f64) -> TileChecks {
    let n = TILE_POSTS;
    let np = n * n;
    let mean = t.height_mm.iter().map(|&v| v as f64).sum::<f64>() / np as f64 / 1000.0;
    let s = slopes_deg(&t.height_mm);
    let sb = slopes_deg(&t.base_mm);
    let h = hist(&s);
    let hb = hist(&sb);
    let (mi, mx) = h.iter().enumerate().fold((0, 0.0), |a, (i, &v)| if v > a.1 { (i, v) } else { a });
    let (bim, peaks) = bimodal(&h);
    let gt = |v: &[f64]| v.iter().filter(|&&d| d.to_radians().tan() > 0.6).count() as f64 / v.len() as f64;
    // F5 interior pits (strictly lower than all 8 neighbours, free posts)
    let mut pits = 0u64;
    for j in 1..n - 1 { for i in 1..n - 1 {
        let k = j * n + i;
        if t.fixed[k] != 0 { continue; }
        let z = t.height_mm[k];
        let mut low = true;
        for (dx, dy) in [(1i64, 0i64), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)] {
            let q = (j as i64 + dy) as usize * n + (i as i64 + dx) as usize;
            if t.height_mm[q] <= z { low = false; break; }
        }
        if low { pits += 1; }
    } }
    // F1: steepest descent from every T1a-interior post (fixed == 1 marks any channel interior; we test all)
    let (mut f1n, mut leave, mut pit) = (0u64, 0u64, 0u64);
    for k in 0..np {
        if t.fixed[k] != 1 { continue; }
        f1n += 1;
        let mut cur = k;
        for _ in 0..4 * n {
            let (i, j) = ((cur % n) as i64, (cur / n) as i64);
            if i == 0 || j == 0 || i == n as i64 - 1 || j == n as i64 - 1 { break; }
            let z = t.height_mm[cur] as f64;
            let mut best = (0.0f64, usize::MAX);
            for (dx, dy) in [(1i64, 0i64), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)] {
                let q = (j + dy) as usize * n + (i + dx) as usize;
                let d = if dx != 0 && dy != 0 { std::f64::consts::SQRT_2 } else { 1.0 };
                let sl = (z - t.height_mm[q] as f64) / d;
                if sl > best.0 { best = (sl, q); }
            }
            if best.1 == usize::MAX {
                // flat or pit: allowed if a channel neighbour at equal height continues (flat bed); else pit
                pit += 1; break;
            }
            if t.fixed[best.1] == 0 { leave += 1; break; }
            if t.fixed[best.1] == 3 { break; }
            cur = best.1;
        }
    }
    // F8c: |T2 - base| outside corridors (fixed 0)
    let mut d: Vec<f64> = (0..np).filter(|&k| t.fixed[k] == 0).map(|k| ((t.height_mm[k] - t.base_mm[k]) as f64 / 1000.0).abs()).collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p99 = if d.is_empty() { 0.0 } else { d[((d.len() - 1) as f64 * 0.99).round() as usize] };
    let dmax = d.last().copied().unwrap_or(0.0);
    let mut ls = std::collections::BTreeMap::new();
    for &c in &t.litho { *ls.entry(c).or_insert(0.0) += 1.0 / np as f64; }
    let layers: std::collections::BTreeSet<u8> = t.layer.iter().cloned().collect();
    let base_mean = t.base_mm.iter().map(|&v| v as f64).sum::<f64>() / np as f64 / 1000.0;
    TileChecks {
        key: format!("0.{}.{}", t.ix, t.iy),
        f6_t2_mean_m: mean, f6_base_mean_m: base_mean, f6_diff_m: mean - t.base_mean_m, f6_pass: (mean - t.base_mean_m).abs() <= 0.25,
        f6_raw_t1b_diff_m: mean - t.raw_t1b_mean_m,
        slope_max_bin_share: mx, slope_max_bin_deg: mi, slope_peaks_deg: peaks, slope_bimodal: bim,
        slope_mean_deg: s.iter().sum::<f64>() / s.len() as f64, base_slope_mean_deg: sb.iter().sum::<f64>() / sb.len() as f64,
        base_slope_max_bin_share: hb.iter().cloned().fold(0.0, f64::max),
        slope_hist_1deg: h,
        f10_share_gt_600permille: gt(&s), f10_base_share_gt_600permille: gt(&sb), f10_points_over_base: 100.0 * (gt(&s) - gt(&sb)),
        f5_interior_pits: pits,
        f1_t1a_interior_posts: f1n, f1_paths_leaving_corridor: leave, f1_paths_ending_in_pit: pit, f1_pass: leave == 0 && pit == 0,
        f8c_abs_diff_p99_m: p99, f8c_abs_diff_max_m: dmax, f8c_bound_m: f8_bound_m, f8c_pass: dmax <= f8_bound_m,
        litho_share: ls.into_iter().collect(), layers_exposed: layers.len(),
    }
}
