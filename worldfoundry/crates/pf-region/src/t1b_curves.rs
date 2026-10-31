//! Curved guide valleys for T1b 0.7.0 (lane T3b, TER-REGION v4, 2026-10-13).
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon.
//!
//! The chart-wide guide network (`t1b_guides`) is a D8 graph on 512 m nodes;
//! rasterizing each link as a straight 16-post line gives straight 45 degree
//! slots. Here the network is cut into SEGMENTS between break nodes (heads,
//! junctions with two or more guide donors, and the non-guide node a guide
//! drains into), and each segment's polyline (x, y, bed) is smoothed by
//! Chaikin corner cutting with its two ends pinned (Chaikin 1974, Computer
//! Graphics and Image Processing 3:346). Junctions and the T1a channel nodes
//! the guides end on stay exactly where they are; bends between them become
//! curves. Pure geometry, no draws: identical for every window.
//! Chaikin of a sequence that never rises keeps it non-rising (each new point
//! is a convex combination of two neighbours, in order), so beds still fall
//! downstream along each segment.

use crate::t1b_guides::GuideView;
use crate::*;

/// One smoothed segment: points in lod-4 post units (global gx, gy) with the bed (m).
pub struct GuideCurve { pub pts: Vec<[f64; 3]>, pub terminal: u32, pub end_node: usize }

/// The chart-wide smoothed guide segments, in a fixed order (by start node index).
pub fn guide_curves(gv: &GuideView, iterations: u32) -> Vec<GuideCurve> {
    let mut donors = vec![0u32; NODES];
    for n in 0..NODES { if gv.is_guide(n) { let r = gv.receiver[n]; if r != NONE { donors[r as usize] += 1; } } }
    let mut out = Vec::new();
    for s in 0..NODES {
        if !gv.is_guide(s) || donors[s] == 1 || gv.receiver[s] == NONE { continue; }
        let mut chain = vec![s];
        let mut cur = s;
        loop {
            let r = gv.receiver[cur] as usize;
            chain.push(r);
            if !gv.is_guide(r) || donors[r] != 1 || gv.receiver[r] == NONE { break; }
            cur = r;
        }
        let mut pts: Vec<[f64; 3]> = chain.iter().map(|&n| {
            let (gx8, gy8) = node_g(n);
            [(gx8 * crate::t1b::L8_TO_L4) as f64, (gy8 * crate::t1b::L8_TO_L4) as f64, gv.bed_m(n)]
        }).collect();
        for _ in 0..iterations { pts = chaikin_pinned(&pts); }
        out.push(GuideCurve { pts, terminal: gv.terminal[s], end_node: *chain.last().unwrap() });
    }
    out
}

/// One Chaikin pass on an open polyline with both ends kept.
pub fn chaikin_pinned(p: &[[f64; 3]]) -> Vec<[f64; 3]> {
    if p.len() < 3 { return p.to_vec(); }
    let mut q = Vec::with_capacity(2 * p.len());
    q.push(p[0]);
    let m = p.len() - 1;
    for k in 0..m {
        let (a, b) = (p[k], p[k + 1]);
        let mix = |u: f64| [a[0] * (1.0 - u) + b[0] * u, a[1] * (1.0 - u) + b[1] * u, a[2] * (1.0 - u) + b[2] * u];
        // The first and last edges keep their pinned end (no cut next to it).
        if k > 0 { q.push(mix(0.25)); }
        if k + 1 < m { q.push(mix(0.75)); }
    }
    q.push(p[m]);
    q
}

/// A river mouth carried from the old (T0 cell-edge) coast down to the zero-contour shore (coast_mode 2).
pub struct MouthExtension { pub curve: GuideCurve, pub from_channel: bool, pub length_m: f64, pub reached_sea: bool }

/// Mouth extensions, chart-wide and in a fixed order (see `ReliefParamsV1::with_coast`, rule t0_zero_contour_coast_v2).
pub fn mouth_extensions(t0: &crate::t0::T0, sk: &crate::t1b_view::SkeletonView, gv: Option<&GuideView>, chart: &vz_tilekey::chart::Chart, passes: u32, plain_ratio: f64) -> R<Vec<MouthExtension>> {
    use crate::t0::CellWater;
    use vz_tilekey::chart::RAD2DEG;
    let sea = t0.sea_level_m;
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let lim = (G_MAX * crate::t1b::L8_TO_L4) as f64 - 1.0;
    let bil = |x: f64, y: f64| -> R<f64> {
        let (la, lo) = chart.inverse_mm((x * crate::t1b::T1B_SPACING_MM as f64).round() as i64, (y * crate::t1b::T1B_SPACING_MM as f64).round() as i64).map_err(|e| e.to_string())?;
        Ok(crate::t1b::bilinear(t0, &elev, la * RAD2DEG, lo * RAD2DEG, false))
    };
    let marine_cell = |n: usize| t0.water(sk.t0_cell[n]) == CellWater::Marine;
    // (end node, from_channel) -> (start level, terminal node)
    let mut starts: std::collections::BTreeMap<(usize, bool), (f64, u32)> = std::collections::BTreeMap::new();
    for n8 in 0..NODES {
        if !(sk.kind[n8] == 1 && sk.channel[n8] == 1) { continue; }
        let r = sk.receiver[n8];
        if r == NONE || r == OFF_CHART { continue; }
        let r = r as usize;
        if !(marine_cell(r) || sk.kind[r] == 3) { continue; }
        let lev = sk.trunk_level_m[r] as f64;
        let e = starts.entry((r, true)).or_insert((lev, n8 as u32));
        if lev < e.0 { *e = (lev, n8 as u32); }
    }
    if let Some(gv) = gv {
        for c in guide_curves(gv, passes) {
            let r = c.end_node;
            if gv.is_guide(r) || !marine_cell(r) { continue; }
            let lev = gv.bed_m(r);
            let e = starts.entry((r, false)).or_insert((lev, r as u32));
            if lev < e.0 { *e = (lev, r as u32); }
        }
    }
    let mut out = Vec::new();
    for (&(node, from_channel), &(l0, terminal)) in &starts {
        let (gx8, gy8) = node_g(node);
        let (mut x, mut y) = ((gx8 * crate::t1b::L8_TO_L4) as f64, (gy8 * crate::t1b::L8_TO_L4) as f64);
        if bil(x, y)? <= sea { continue; }
        let mut pts = vec![[x, y, bil(x, y)?]];
        let mut reached = false;
        for _ in 0..12_500 {
            if x.abs() > lim || y.abs() > lim { break; }
            let gxv = (bil(x + 1.0, y)? - bil(x - 1.0, y)?) / 2.0;
            let gyv = (bil(x, y + 1.0)? - bil(x, y - 1.0)?) / 2.0;
            let g = (gxv * gxv + gyv * gyv).sqrt();
            if g < 1e-9 { break; }
            x -= 0.5 * gxv / g; y -= 0.5 * gyv / g;
            let b = bil(x, y)?;
            pts.push([x, y, b]);
            if b <= sea { reached = true; break; }
        }
        let m = (pts.len() - 1).max(1) as f64;
        let l1 = l0.min(sea);
        // Bed: linear from the end node's level to the shore, and never above HALF the reclaimed-plain height there (the
        // river is incised into the plain; both fall along the path, so the bed never rises).
        let curve = GuideCurve { pts: pts.iter().enumerate().map(|(k, q)| [q[0], q[1], (l0 + (l1 - l0) * k as f64 / m).min(sea + 0.5 * plain_ratio * (q[2] - sea).max(0.0))]).collect(), terminal, end_node: node };
        out.push(MouthExtension { curve, from_channel, length_m: (pts.len() - 1) as f64 * 0.5 * crate::t1b::T1B_SPACING_M, reached_sea: reached });
    }
    Ok(out)
}

/// Rasterize a curve at the 32 m posts: samples every <= 1/4 post along each piece, nearest post (round half away
/// from zero), bed linear along the piece. Consecutive samples land on the same or an 8-adjacent post, so the
/// rasterized valley is 8-connected. Calls `put(gx, gy, bed)` for every sample.
pub fn rasterize(c: &GuideCurve, put: &mut dyn FnMut(i64, i64, f64)) {
    for w in c.pts.windows(2) {
        let (a, b) = (w[0], w[1]);
        let len = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
        let ns = ((len * 4.0).ceil() as usize).max(1);
        for s in 0..=ns {
            let t = s as f64 / ns as f64;
            put((a[0] + t * (b[0] - a[0])).round() as i64, (a[1] + t * (b[1] - a[1])).round() as i64, a[2] + t * (b[2] - a[2]));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn chaikin_keeps_ends_and_monotone_beds_and_rounds_a_stair() {
        // A 512 m stair (E, N, E, N) in lod-4 units with falling beds.
        let p = vec![[0.0, 0.0, 100.0], [16.0, 0.0, 90.0], [16.0, 16.0, 80.0], [32.0, 16.0, 70.0], [32.0, 32.0, 60.0]];
        let mut q = p.clone();
        for _ in 0..3 { q = chaikin_pinned(&q); }
        assert_eq!(q[0], p[0]);
        assert_eq!(*q.last().unwrap(), *p.last().unwrap());
        for w in q.windows(2) { assert!(w[1][2] <= w[0][2] + 1e-12, "bed rises"); }
        // The interior corners are cut: no output point sits on the corner (16, 0).
        assert!(q.iter().all(|v| !(v[0] == 16.0 && v[1] == 0.0)));
        // 8-connected raster.
        let c = GuideCurve { pts: q, terminal: 0, end_node: 0 };
        let mut last: Option<(i64, i64)> = None;
        rasterize(&c, &mut |x, y, _| { if let Some((lx, ly)) = last { assert!((x - lx).abs() <= 1 && (y - ly).abs() <= 1); } last = Some((x, y)); });
    }
}
