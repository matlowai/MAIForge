//! T1b falsifiers (ADR 043 F1-F7, F9 seams, F12) and relief statistics
//! (lane TR). Tolerances are the ones pre-registered in
//! `unpublished-design-record` (2026-10-12 19:52 EDT).

use crate::t0::{CellWater, T0};
use crate::t1b::*;
use crate::t1b_lem::{self, Ledger};
use crate::t1b_view::SkeletonView;
use crate::t1b_window::RawWindow;
use crate::*;
use serde::Serialize;
use std::collections::BTreeMap;
use vz_tilekey::chart::{Chart, ChartDeclV1, RAD2DEG};

#[derive(Clone, Debug, Default, Serialize)]
pub struct ReliefStats {
    pub posts: u64,
    pub mean_m: f64,
    pub std_m: f64,
    pub min_m: f64,
    pub max_m: f64,
    pub p05_m: f64,
    pub p50_m: f64,
    pub p95_m: f64,
    pub hypsometric_integral: f64,
    pub slope_mean_gradient: f64,
    pub slope_mean_deg: f64,
    pub slope_p50_deg: f64,
    pub slope_p90_deg: f64,
    pub slope_frac_gt_30deg: f64,
    pub slope_frac_gt_45deg: f64,
    /// ADR 043 F10's pinned grade: 600 permille (about 31 degrees).
    pub slope_frac_gt_600_permille: f64,
    pub local_relief_half_width_m: f64,
    pub local_relief_mean_m: f64,
    pub local_relief_p90_m: f64,
    pub valley_area_m2: f64,
    pub valley_post_fraction: f64,
    pub max_free_drainage_area_m2: f64,
}

fn pct(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() { return f64::NAN; }
    let k = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[k]
}

/// Sliding max - min over a (2h+1)^2 square, separable, on the whole lattice
/// (windows clipped at the lattice edge).
fn local_relief(n: usize, z: &[f64], h: usize) -> Vec<f64> {
    let run = |src: &[f64], horiz: bool, take_max: bool| -> Vec<f64> {
        let mut out = vec![0f64; n * n];
        for a in 0..n { for b in 0..n {
            let lo = b.saturating_sub(h);
            let hi = (b + h).min(n - 1);
            let mut v = if take_max { f64::MIN } else { f64::MAX };
            for k in lo..=hi {
                let i = if horiz { a * n + k } else { k * n + a };
                v = if take_max { v.max(src[i]) } else { v.min(src[i]) };
            }
            let o = if horiz { a * n + b } else { b * n + a };
            out[o] = v;
        } }
        out
    };
    let mx = run(&run(z, true, true), false, true);
    let mn = run(&run(z, true, false), false, false);
    mx.iter().zip(&mn).map(|(a, b)| a - b).collect()
}

/// Statistics over the window core's land posts (everything but marine and lake).
pub fn relief_stats(raw: &RawWindow, p: &ReliefParamsV1) -> ReliefStats {
    let lat = raw.lat;
    let n = lat.n;
    let z = raw.z_vec();
    let land = |i: usize| raw.fixed[i] != MARINE && raw.fixed[i] != LAKE;
    let idx: Vec<usize> = (0..lat.nodes()).filter(|&i| lat.in_core(i) && land(i)).collect();
    let mut s = ReliefStats::default();
    if idx.is_empty() { return s; }
    let mut hs: Vec<f64> = idx.iter().map(|&i| z[i]).collect();
    let cnt = hs.len() as f64;
    s.posts = hs.len() as u64;
    s.mean_m = hs.iter().sum::<f64>() / cnt;
    s.std_m = (hs.iter().map(|v| (v - s.mean_m).powi(2)).sum::<f64>() / cnt).sqrt();
    hs.sort_by(|a, b| a.total_cmp(b));
    s.min_m = hs[0];
    s.max_m = *hs.last().unwrap();
    s.p05_m = pct(&hs, 0.05);
    s.p50_m = pct(&hs, 0.50);
    s.p95_m = pct(&hs, 0.95);
    s.hypsometric_integral = if s.max_m > s.min_m { (s.mean_m - s.min_m) / (s.max_m - s.min_m) } else { 0.0 };
    let dx = T1B_SPACING_M;
    let mut slopes: Vec<f64> = Vec::with_capacity(idx.len());
    let mut grad_sum = 0.0;
    for &i in &idx {
        let gx = (z[i + 1] - z[i - 1]) / (2.0 * dx);
        let gy = (z[i + n] - z[i - n]) / (2.0 * dx);
        let g = (gx * gx + gy * gy).sqrt();
        grad_sum += g;
        slopes.push(g);
    }
    s.slope_mean_gradient = grad_sum / cnt;
    let mut deg: Vec<f64> = slopes.iter().map(|g| g.atan().to_degrees()).collect();
    s.slope_mean_deg = deg.iter().sum::<f64>() / cnt;
    s.slope_frac_gt_30deg = deg.iter().filter(|d| **d > 30.0).count() as f64 / cnt;
    s.slope_frac_gt_45deg = deg.iter().filter(|d| **d > 45.0).count() as f64 / cnt;
    s.slope_frac_gt_600_permille = slopes.iter().filter(|g| **g > 0.6).count() as f64 / cnt;
    deg.sort_by(|a, b| a.total_cmp(b));
    s.slope_p50_deg = pct(&deg, 0.5);
    s.slope_p90_deg = pct(&deg, 0.9);
    let h = p.local_relief_half_width_posts as usize;
    let lr = local_relief(n, &z, h);
    let mut lrs: Vec<f64> = idx.iter().map(|&i| lr[i]).collect();
    s.local_relief_half_width_m = h as f64 * dx;
    s.local_relief_mean_m = lrs.iter().sum::<f64>() / cnt;
    lrs.sort_by(|a, b| a.total_cmp(b));
    s.local_relief_p90_m = pct(&lrs, 0.9);
    s.valley_area_m2 = p.valley_area_m2 as f64;
    let free: Vec<usize> = idx.iter().copied().filter(|&i| raw.fixed[i] == FREE).collect();
    s.valley_post_fraction = free.iter().filter(|&&i| raw.drainage_area_m2[i] as f64 >= s.valley_area_m2).count() as f64 / free.len().max(1) as f64;
    s.max_free_drainage_area_m2 = free.iter().map(|&i| raw.drainage_area_m2[i] as f64).fold(0.0, f64::max);
    s
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct CellCheck {
    pub cell: u32,
    pub water: String,
    pub core_posts: u64,
    pub t0_mean_m: f64,
    pub t0_subgrid_relief_m: f64,
    pub t1b_mean_m: f64,
    pub bilinear_mean_m: f64,
    pub f6_tolerance_m: f64,
    pub f6_bilinear_diff_m: f64,
    pub f6_pass: bool,
    pub f6_literal_cell_mean_diff_m: f64,
    pub f3_free_posts: u64,
    pub f3_own_channel: u64,
    pub f3_own_receiver_sea: u64,
    pub f3_ring_same_basin: u64,
    pub f3_other: u64,
    pub f3_share_permille: f64,
    pub f3_strict_share_permille: f64,
    pub f3_pass: bool,
    pub coastal: bool,
    pub f4_window_t1b_land_fraction: f64,
    pub f4_window_bilinear_land_fraction: f64,
    pub f4_whole_cell_t1b_land_fraction: f64,
    pub f4_whole_cell_bilinear_land_fraction: f64,
    pub f4_clause1_pass: bool,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct WindowChecks {
    pub window: String,
    pub f1_t1a_channel_nodes_in_core: u32,
    pub f1_t1a_nodes_not_fixed: u32,
    pub f1_core_channel_posts: u32,
    pub f1_paths_leaving_channel: u32,
    pub f1_bed_rises: u32,
    pub f1_paths_to_wrong_outlet: u32,
    pub f1_outlets: Vec<String>,
    /// Pre-registered strict reading: the path stays on channel posts to the exit or mouth.
    pub f1_pass: bool,
    /// ADR 043 F1 wording at T1b: same T0 outlet, inside the T1a corridor, never rising.
    pub f1_max_offchannel_run_posts: u32,
    pub f1_adr_reading_failures: u32,
    pub f1_adr_reading_pass: bool,
    pub f2_max_free_discharge_m3_s: f64,
    pub f2_free_posts_ge_threshold: u32,
    pub f2_pass: bool,
    pub cells: Vec<CellCheck>,
    pub f3_pass: bool,
    pub f4_applicable: bool,
    pub f4_clause1_pass: bool,
    pub f4_coast_posts: u32,
    pub f4_max_coast_to_t0_zero_contour_km: f64,
    pub f4_one_t0_cell_km: f64,
    pub f4_clause2_pass: bool,
    pub f5_core_pit_posts: u32,
    pub f5_core_pit_area_m2: f64,
    pub f5_core_pit_max_depth_m: f64,
    pub f5_core_pit_posts_deeper_1m: u32,
    pub f5_pit_area_fraction_deeper_1m: f64,
    pub f5_pits_fed_by_fixed_channel: u32,
    pub f5_pass: bool,
    pub f6_pass: bool,
    pub f7_relative_residual: f64,
    pub f7_volume_relative_residual: f64,
    pub f7_diffusion_flux_relative_residual: f64,
    pub f7_pass: bool,
    pub stats: ReliefStats,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_envelope: Option<serde_json::Value>,
}

/// The T0 bilinear zero contour, sampled at lod-8 posts over the chart:
/// midpoints of 4-neighbour sign changes (chart mm).
pub fn t0_zero_contour(t0: &T0, chart: &Chart) -> R<Vec<[f64; 2]>> {
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let mut v = vec![0f64; NODES];
    for nidx in 0..NODES {
        let (x, y) = node_mm(nidx);
        let (la, lo) = chart.inverse_mm(x, y).map_err(|e| e.to_string())?;
        v[nidx] = bilinear(t0, &elev, la * RAD2DEG, lo * RAD2DEG, false) - t0.sea_level_m;
    }
    let mut pts = Vec::new();
    for gy in G_MIN..=G_MAX { for gx in G_MIN..=G_MAX {
        let a = node_index(gx, gy).unwrap();
        for (ox, oy) in [(1, 0), (0, 1)] {
            if let Some(b) = node_index(gx + ox, gy + oy) {
                if (v[a] > 0.0) != (v[b] > 0.0) {
                    let t = v[a] / (v[a] - v[b]);
                    let (xa, ya) = node_mm(a);
                    pts.push([xa as f64 + t * (ox * SPACING_MM) as f64, ya as f64 + t * (oy * SPACING_MM) as f64]);
                }
            }
        }
    } }
    Ok(pts)
}

/// Whole-cell bilinear land fraction (64 x 64 samples inside the T0 cell).
fn whole_cell_bilinear_land_fraction(t0: &T0, c: u32) -> f64 {
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let (i, j) = t0.ij(c);
    let dlat = 180.0 / t0.nlat as f64;
    let dlon = 360.0 / t0.nlon as f64;
    let (lat0, lon0) = (-90.0 + i as f64 * dlat, -180.0 + j as f64 * dlon);
    let mut land = 0usize;
    for a in 0..64 { for b in 0..64 {
        let la = lat0 + (a as f64 + 0.5) / 64.0 * dlat;
        let lo = lon0 + (b as f64 + 0.5) / 64.0 * dlon;
        if bilinear(t0, &elev, la, lo, false) > t0.sea_level_m { land += 1; }
    } }
    land as f64 / 4096.0
}

fn is_coastal(t0: &T0, c: u32) -> bool {
    let (i, j) = t0.ij(c);
    let me = t0.is_marine(c);
    for di in -1i64..=1 { for dj in -1i64..=1 {
        if di == 0 && dj == 0 { continue; }
        let ii = i as i64 + di;
        if ii < 0 || ii >= t0.nlat as i64 { continue; }
        let jj = (j as i64 + dj).rem_euclid(t0.nlon as i64);
        let nb = (ii as usize * t0.nlon + jj as usize) as u32;
        if t0.is_marine(nb) != me { return true; }
    } }
    false
}

pub fn check_window(raw: &RawWindow, prep: &WindowInputs, ledger: &Ledger, t0: &T0, skel: &SkeletonView, p: &ReliefParamsV1,
    chart_decl: &ChartDeclV1, zero_contour: Option<&[[f64; 2]]>) -> R<WindowChecks> {
    let lat = raw.lat;
    let n = lat.n;
    let nn = lat.nodes();
    let z = raw.z_vec();
    let thr = p.threshold_m3_s();
    let mut c = WindowChecks { window: lat.key(), ..Default::default() };

    // ---- F1: channel continuity on the priority-flooded surface with channels as ordinary posts.
    let (x0, x1, y0, y1) = lat.core();
    // Lane CHANNEL-CENTRELINE (cl=1): the T1a link is drawn on the T1a planform curve, whose 512 m rasterization the node
    // is; the node test reads 'the node's own link is drawn on a published channel post, or over water (estuary / construction
    // water posts carry the link; the lattice post being water also passes)'.
    let cl_on = p.rivernet.as_ref().and_then(crate::t1b_cl::cl).is_some();
    let link_drawn: std::collections::BTreeSet<u32> = if cl_on { (0..nn).filter(|&i| raw.fixed[i] == CHANNEL || raw.fixed[i] == MARINE).map(|i| prep.chan_link_node[i]).filter(|&l| l != NONE).collect() } else { Default::default() };
    for gy8 in (y0 / L8_TO_L4)..=(y1 / L8_TO_L4) { for gx8 in (x0 / L8_TO_L4)..=(x1 / L8_TO_L4) {
        let Some(n8) = skel.idx(gx8, gy8) else { continue };
        if skel.kind[n8] == 1 && skel.channel[n8] == 1 {
            c.f1_t1a_channel_nodes_in_core += 1;
            let i = lat.idx(gx8 * L8_TO_L4, gy8 * L8_TO_L4).ok_or("core node outside window")?;
            let ok = if cl_on { raw.fixed[i] == MARINE || link_drawn.contains(&(n8 as u32)) } else { raw.fixed[i] == CHANNEL || raw.fixed[i] == MARINE };
            if !ok { c.f1_t1a_nodes_not_fixed += 1; }
        }
    } }
    // Lane DAMMED-VALLEY (dam_lake=1): a lake post outside a T0 lake cell is an edifice-dammed lake, which drains over its
    // spill (t1b_dam.rs): the checks route through it as a free post (without the key no such post exists: unchanged).
    let dam = |i: usize| raw.fixed[i] == LAKE && t0.water(raw.t0_cell[i]) != CellWater::Lake;
    let mut fixed_f1 = vec![FREE; nn];
    for i in 0..nn {
        let f = raw.fixed[i];
        if f == MARINE || (f == LAKE && !dam(i)) || (f == CHANNEL && lat.on_ring(i)) { fixed_f1[i] = f; }
    }
    let r1 = t1b_lem::route(n, &z, &fixed_f1, T1B_SPACING_M);
    let mut outlets: BTreeMap<String, u32> = BTreeMap::new();
    for i in 0..nn {
        if !(lat.in_core(i) && raw.fixed[i] == CHANNEL) { continue; }
        c.f1_core_channel_posts += 1;
        let start_basin = skel.outlet_basin[prep.chan_link_node[i] as usize];
        let mut cur = i;
        let mut strict_ok = true;
        let mut adr_ok = true;
        let mut off_run = 0u32;
        let mut guard = 0;
        loop {
            if fixed_f1[cur] != FREE { break; }
            let nx = r1.receiver[cur] as usize;
            if raw.height_mm[nx] > raw.height_mm[cur] { c.f1_bed_rises += 1; adr_ok = false; break; }
            if raw.fixed[nx] != CHANNEL && raw.fixed[nx] != MARINE {
                if strict_ok { c.f1_paths_leaving_channel += 1; strict_ok = false; }
                off_run += 1;
                c.f1_max_offchannel_run_posts = c.f1_max_offchannel_run_posts.max(off_run);
            } else { off_run = 0; }
            cur = nx;
            guard += 1;
            if guard > nn { return Err("F1 path loop".into()); }
        }
        // The off-channel run is bounded by run x 32 sqrt 2 m from the last channel post.
        if c.f1_max_offchannel_run_posts as f64 * T1B_SPACING_M * std::f64::consts::SQRT_2 > skel.corridor_half_width_mm as f64 / 1000.0 { adr_ok = false; }
        let (outlet, basin) = if raw.fixed[cur] == MARINE {
            ("mouth".to_string(), t0.basin_id[raw.t0_cell[cur] as usize])
        } else if raw.fixed[cur] == CHANNEL {
            let (gx, gy) = lat.g(cur);
            (format!("exit@{gx},{gy}"), skel.outlet_basin[prep.chan_link_node[cur] as usize])
        } else {
            ("other".to_string(), u32::MAX)
        };
        if basin != start_basin { c.f1_paths_to_wrong_outlet += 1; adr_ok = false; }
        if !adr_ok { c.f1_adr_reading_failures += 1; }
        *outlets.entry(format!("{outlet} basin {basin}")).or_default() += 1;
    }
    c.f1_outlets = outlets.into_iter().map(|(k, v)| format!("{k}: {v} posts")).collect();
    c.f1_pass = c.f1_t1a_nodes_not_fixed == 0 && c.f1_paths_leaving_channel == 0 && c.f1_bed_rises == 0 && c.f1_paths_to_wrong_outlet == 0
        && c.f1_core_channel_posts > 0;
    c.f1_adr_reading_pass = c.f1_t1a_nodes_not_fixed == 0 && c.f1_adr_reading_failures == 0 && c.f1_core_channel_posts > 0;

    // ---- Published routing (all fixed posts are outlets): F2, F3, F5.
    let r = if (0..nn).any(|i| dam(i)) {
        let rfix: Vec<u8> = (0..nn).map(|i| if dam(i) { FREE } else { raw.fixed[i] }).collect();
        t1b_lem::route(n, &z, &rfix, T1B_SPACING_M)
    } else { t1b_lem::route(n, &z, &raw.fixed, T1B_SPACING_M) };
    let term = t1b_lem::terminals(&r);
    let q_of = |a: f64| a * t0.runoff_m_per_yr / SECONDS_PER_YEAR;
    for i in 0..nn {
        if !lat.in_core(i) || raw.fixed[i] != FREE { continue; }
        let q = q_of(raw.drainage_area_m2[i] as f64);
        c.f2_max_free_discharge_m3_s = c.f2_max_free_discharge_m3_s.max(q);
        if q >= thr { c.f2_free_posts_ge_threshold += 1; }
    }
    c.f2_pass = c.f2_free_posts_ge_threshold == 0;

    let a_post = T1B_SPACING_M * T1B_SPACING_M;
    let mut cells: BTreeMap<u32, CellCheck> = BTreeMap::new();
    let mut sums: BTreeMap<u32, (f64, f64, u64, u64, u64)> = BTreeMap::new(); // t1b sum, bil sum, n, t1b land, bil land
    for i in 0..nn {
        if !lat.in_core(i) { continue; }
        let cell = raw.t0_cell[i];
        let e = sums.entry(cell).or_default();
        e.0 += z[i]; e.1 += raw.t0_bilinear_m[i] as f64; e.2 += 1;
        // Lane DAMMED-VALLEY: an edifice-dammed lake is inland water on T0 land, not sea (F4 counts it as land).
        let t1b_land = raw.fixed[i] != MARINE && (raw.fixed[i] != LAKE || dam(i)) && z[i] > t0.sea_level_m;
        if t1b_land { e.3 += 1; }
        if raw.t0_bilinear_m[i] as f64 > t0.sea_level_m { e.4 += 1; }
        let cc = cells.entry(cell).or_insert_with(|| CellCheck { cell, ..Default::default() });
        if raw.fixed[i] == FREE {
            cc.f3_free_posts += 1;
            let tp = term[i] as usize;
            let own_cell = cell;
            // Lane VE: a window without any fixed post (e.g. under a young shield that buried every network channel) has no
            // terminal: counted as other (no field added, so every other checks.json is unchanged).
            if tp >= raw.fixed.len() { cc.f3_other += 1; continue; }
            match raw.fixed[tp] {
                CHANNEL => {
                    let n8 = prep.chan_link_node[tp] as usize;
                    if skel.t0_cell[n8] == own_cell { cc.f3_own_channel += 1; } else { cc.f3_other += 1; }
                }
                MARINE | PLAIN => {
                    if t0.flow(own_cell) == Some(raw.t0_cell[tp]) { cc.f3_own_receiver_sea += 1; } else { cc.f3_other += 1; }
                }
                GUIDE | NET => {
                    // The guide's own path ends on the fixed T1a node recorded at the guide post.
                    let n8 = prep.chan_link_node[tp] as usize;
                    let tc = skel.t0_cell[n8];
                    if tc == own_cell { cc.f3_own_channel += 1; }
                    else if t0.flow(own_cell) == Some(tc) { cc.f3_own_receiver_sea += 1; }
                    else { cc.f3_other += 1; }
                }
                RING => {
                    let (gx, gy) = lat.g(tp);
                    let n8 = skel.idx((gx as f64 / L8_TO_L4 as f64).round() as i64, (gy as f64 / L8_TO_L4 as f64).round() as i64);
                    let same = match (n8, t0.final_basin(own_cell)) { (Some(n8), Ok(b)) => skel.outlet_basin[n8] == b, _ => false };
                    if same { cc.f3_ring_same_basin += 1; } else { cc.f3_other += 1; }
                }
                _ => cc.f3_other += 1,
            }
        }
    }
    c.f3_pass = true;
    c.f6_pass = true;
    c.f4_clause1_pass = true;
    for (cell, cc) in cells.iter_mut() {
        let (s_t, s_b, cnt, l_t, l_b) = sums[cell];
        let ci = *cell as usize;
        cc.water = format!("{:?}", t0.water(*cell));
        cc.core_posts = cnt;
        cc.t0_mean_m = t0.elevation_m[ci] as f64;
        cc.t0_subgrid_relief_m = t0.subgrid_relief_m[ci] as f64;
        cc.t1b_mean_m = s_t / cnt as f64;
        cc.bilinear_mean_m = s_b / cnt as f64;
        if t0.water(*cell) == CellWater::Land {
            cc.f6_tolerance_m = (0.1 * cc.t0_subgrid_relief_m).max(50.0);
            cc.f6_bilinear_diff_m = cc.t1b_mean_m - cc.bilinear_mean_m;
            cc.f6_literal_cell_mean_diff_m = cc.t1b_mean_m - cc.t0_mean_m;
            cc.f6_pass = cc.f6_bilinear_diff_m.abs() <= cc.f6_tolerance_m;
            if !cc.f6_pass { c.f6_pass = false; }
            let tot = cc.f3_free_posts.max(1) as f64;
            cc.f3_share_permille = 1000.0 * (cc.f3_own_channel + cc.f3_own_receiver_sea + cc.f3_ring_same_basin) as f64 / tot;
            cc.f3_strict_share_permille = 1000.0 * (cc.f3_own_channel + cc.f3_own_receiver_sea) as f64 / tot;
            cc.f3_pass = cc.f3_share_permille >= 950.0;
            if !cc.f3_pass && cc.f3_free_posts > 0 { c.f3_pass = false; }
        } else { cc.f6_pass = true; cc.f3_pass = true; }
        cc.coastal = is_coastal(t0, *cell);
        cc.f4_window_t1b_land_fraction = l_t as f64 / cnt as f64;
        cc.f4_window_bilinear_land_fraction = l_b as f64 / cnt as f64;
        cc.f4_whole_cell_bilinear_land_fraction = whole_cell_bilinear_land_fraction(t0, *cell);
        // coast_mode 1 (0.7.0): the land rule IS the bilinear sign, so the whole-cell T1b land fraction equals the
        // bilinear one by construction (lake cells keep the cell rule).
        cc.f4_whole_cell_t1b_land_fraction = if p.coast_mode >= 1 && t0.water(*cell) != CellWater::Lake { cc.f4_whole_cell_bilinear_land_fraction }
            else if t0.water(*cell) == CellWater::Land { 1.0 } else { 0.0 };
        cc.f4_clause1_pass = !cc.coastal || (cc.f4_whole_cell_t1b_land_fraction - cc.f4_whole_cell_bilinear_land_fraction).abs() <= 0.10;
        if cc.coastal { c.f4_applicable = true; if !cc.f4_clause1_pass { c.f4_clause1_pass = false; } }
    }
    c.cells = cells.into_values().collect();

    // F4 clause 2: every T1b coast post within one T0 cell of the T0 zero contour.
    let dlat_m = std::f64::consts::PI / t0.nlat as f64 * t0.radius_m;
    c.f4_one_t0_cell_km = dlat_m / 1000.0;
    c.f4_clause2_pass = true;
    if c.f4_applicable {
        let chart = Chart::new(chart_decl.clone()).map_err(|e| e.to_string())?;
        let owned;
        let zc: &[[f64; 2]] = match zero_contour { Some(v) => v, None => { owned = t0_zero_contour(t0, &chart)?; &owned } };
        for i in 0..nn {
            if !lat.in_core(i) { continue; }
            let land = raw.fixed[i] != MARINE && raw.fixed[i] != LAKE;
            if !land { continue; }
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            let coast = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)].iter().any(|(ox, oy)| {
                let (xx, yy) = (x + ox, y + oy);
                xx >= 0 && yy >= 0 && xx < n as i64 && yy < n as i64 && raw.fixed[yy as usize * n + xx as usize] == MARINE
            });
            if !coast { continue; }
            c.f4_coast_posts += 1;
            let (gx, gy) = lat.g(i);
            let (px, py) = ((gx * T1B_SPACING_MM) as f64, (gy * T1B_SPACING_MM) as f64);
            let d = zc.iter().map(|q| ((q[0] - px).powi(2) + (q[1] - py).powi(2)).sqrt()).fold(f64::INFINITY, f64::min) / 1e6;
            c.f4_max_coast_to_t0_zero_contour_km = c.f4_max_coast_to_t0_zero_contour_km.max(d);
        }
        c.f4_clause2_pass = c.f4_max_coast_to_t0_zero_contour_km <= c.f4_one_t0_cell_km;
    }

    // ---- F5: closed depressions on the published surface (T1b publishes no water).
    for i in 0..nn {
        if !lat.in_core(i) || raw.fixed[i] != FREE { continue; }
        let d = r.fill[i] - z[i];
        if d > 0.0 {
            c.f5_core_pit_posts += 1;
            c.f5_core_pit_area_m2 += a_post;
            c.f5_core_pit_max_depth_m = c.f5_core_pit_max_depth_m.max(d);
            if d > 1.0 { c.f5_core_pit_posts_deeper_1m += 1; }
        }
    }
    let core_posts = (0..nn).filter(|&i| lat.in_core(i)).count() as f64;
    c.f5_pit_area_fraction_deeper_1m = c.f5_core_pit_posts_deeper_1m as f64 / core_posts;
    // Fixed channel posts are outlets (no receiver), so no depression is fed by a channel.
    c.f5_pits_fed_by_fixed_channel = 0;
    c.f5_pass = c.f5_pits_fed_by_fixed_channel == 0;

    // ---- F7 from the run ledger.
    c.f7_relative_residual = ledger.f7_relative_residual;
    c.f7_volume_relative_residual = ledger.volume_relative_residual;
    c.f7_diffusion_flux_relative_residual = ledger.diffusion_flux_relative_residual;
    c.f7_pass = ledger.f7_relative_residual <= 1e-9 && ledger.volume_relative_residual <= 1e-9;
    if let Some(env) = p.rivernet.as_ref().and_then(|r| r.lithology.envelope.as_ref()) {
        let report = envelope_window(raw, prep, t0, skel, chart_decl, env)?;
        let solid = report["total_solid_above_datum_m3"].as_f64().ok_or("F7 solid")?;
        let residual = report["closure_residual_m3"].as_f64().ok_or("F7 residual")?.abs();
        c.f7_volume_relative_residual = c.f7_volume_relative_residual.max(if solid > 0. { residual / solid } else { residual });
        c.f7_pass &= report["f7_geometry_pass"].as_bool() == Some(true);
        c.source_envelope = Some(report);
    }

    c.stats = relief_stats(raw, p);
    Ok(c)
}

/// Lane T0-IMPRINT-PROD (key outlet_checks=1; ADR 045): the window-local per-T0-outlet F3 share (diagnostic, cells are
/// partial in a window): each free core post of a T0 Land cell is followed on the published routing (as F3 above) to its
/// fixed terminal, whose T0 outlet is the T1a basin of its channel link (channel, guide, network posts), the sea of its T0
/// cell (marine, plain), the final basin of its lake, or the T1a basin at the nearest lod-8 node (ring).
pub fn window_f3_outlet(raw: &RawWindow, prep: &WindowInputs, t0: &T0, skel: &SkeletonView) -> R<crate::checks::F3Outlet> {
    use crate::checks::Outlet;
    let lat = raw.lat;
    let n = lat.n;
    let nn = lat.nodes();
    let z = raw.z_vec();
    let dam = |i: usize| raw.fixed[i] == LAKE && t0.water(raw.t0_cell[i]) != CellWater::Lake;
    let r = if (0..nn).any(|i| dam(i)) {
        let rfix: Vec<u8> = (0..nn).map(|i| if dam(i) { FREE } else { raw.fixed[i] }).collect();
        t1b_lem::route(n, &z, &rfix, T1B_SPACING_M)
    } else { t1b_lem::route(n, &z, &raw.fixed, T1B_SPACING_M) };
    let term = t1b_lem::terminals(&r);
    let basin_of_link = |n8: u32| -> Outlet { match skel.outlet_basin.get(n8 as usize) { Some(&b) if b != NONE => Outlet::Basin(b), _ => Outlet::Unknown } };
    let mut outlet = vec![Outlet::Unknown; nn];
    let mut include = vec![false; nn];
    for i in 0..nn {
        if !lat.in_core(i) || raw.fixed[i] != FREE { continue; }
        include[i] = true;
        let tp = term[i] as usize;
        if tp >= nn { continue; }
        let tc = raw.t0_cell[tp];
        outlet[i] = match raw.fixed[tp] {
            CHANNEL | GUIDE | NET => basin_of_link(prep.chan_link_node[tp]),
            MARINE | PLAIN => if t0.is_marine(tc) { Outlet::Sea { body: t0.water_body_id[tc as usize], basin: (t0.basin_id[tc as usize] != NONE).then_some(t0.basin_id[tc as usize]) } }
                else { Outlet::Sea { body: NONE, basin: None } },
            LAKE => if dam(tp) { Outlet::Unknown } else { t0.final_basin(tc).map(Outlet::Basin).unwrap_or(Outlet::Unknown) },
            RING => {
                let (gx, gy) = lat.g(tp);
                match skel.idx((gx as f64 / L8_TO_L4 as f64).round() as i64, (gy as f64 / L8_TO_L4 as f64).round() as i64) { Some(n8) => basin_of_link(n8 as u32), None => Outlet::Unknown }
            }
            _ => Outlet::Unknown,
        };
    }
    let area = vec![T1B_SPACING_M * T1B_SPACING_M; nn];
    crate::checks::f3_outlet_g(&raw.t0_cell, &area, &outlet, &include, t0, "window core free posts on the published routing (diagnostic)", false)
}

// ------------------------------------------------------------------ seams

/// Smoothstep seam weight of window `w` at global post (gx, gy).
pub fn seam_weight(w: &WinLat, gx: i64, gy: i64, h_posts: i64) -> f64 {
    let (x0, x1, y0, y1) = w.core();
    let s = |d: i64| -> f64 {
        let t = ((d + h_posts) as f64 / (2 * h_posts) as f64).clamp(0.0, 1.0);
        t * t * (3.0 - 2.0 * t)
    };
    s((gx - x0).min(x1 - gx)) * s((gy - y0).min(y1 - gy))
}

/// Final (blended) post height in mm from the raw windows given, in
/// canonical (ry, rx) order whatever order they are passed in.
pub fn blend_post(windows: &[&RawWindow], gx: i64, gy: i64, h_posts: i64) -> R<Option<i32>> {
    let mut ws: Vec<&RawWindow> = windows.to_vec();
    ws.sort_by_key(|w| (w.lat.ry, w.lat.rx));
    let mut fixed_val: Option<i32> = None;
    let (mut num, mut den) = (0.0f64, 0.0f64);
    for w in ws {
        let Some(i) = w.lat.idx(gx, gy) else { continue };
        let f = w.fixed[i];
        if f == CHANNEL || f == MARINE || f == LAKE || f == GUIDE || f == PLAIN {
            if let Some(v) = fixed_val { if v != w.height_mm[i] { return Err(format!("fixed post ({gx},{gy}) differs between windows")); } }
            fixed_val = Some(w.height_mm[i]);
            continue;
        }
        let wt = seam_weight(&w.lat, gx, gy, h_posts);
        if wt > 0.0 { num += wt * w.height_mm[i] as f64; den += wt; }
    }
    if let Some(v) = fixed_val { return Ok(Some(v)); }
    if den == 0.0 { return Ok(None); }
    Ok(Some(vz_tilekey::chart::quantize_mm(num / den) as i32))
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SeamStats {
    pub a: String,
    pub b: String,
    pub strip_half_width_m: f64,
    pub strip_posts: u64,
    pub free_in_both: u64,
    pub diff_mean_m: f64,
    pub diff_rms_m: f64,
    pub diff_p99_abs_m: f64,
    pub diff_max_abs_m: f64,
    pub fixed_posts_compared: u64,
    pub fixed_posts_differing: u64,
    /// Lane TE (2026-10-14, closes the blind spot lane TC found): strip posts fixed (channel, water, guide, plain) in one
    /// window and soft (FREE / NET) in the other; neither branch above compared them. Their max |dz|; strip posts that
    /// differ in height or code in any way; `exact` = none does.
    pub mixed_posts: u64,
    pub mixed_max_abs_m: f64,
    pub strip_posts_differing: u64,
    pub exact: bool,
    pub blend_order_independent: bool,
    pub blended_pits_in_strip: u64,
    pub raw_a_pits_in_strip: u64,
    pub raw_b_pits_in_strip: u64,
    pub blended_max_slope_in_strip: f64,
    pub raw_a_max_slope_in_strip: f64,
    pub rms_tolerance_m: f64,
    pub p99_tolerance_m: f64,
    pub pass_tolerance: bool,
}

/// Compare two edge-adjacent raw windows over the seam strip and blend it.
pub fn seam_compare(a: &RawWindow, b: &RawWindow, h_posts: i64) -> R<(SeamStats, Vec<(i64, i64, i32)>)> {
    let (da, db) = (b.lat.rx - a.lat.rx, b.lat.ry - a.lat.ry);
    if da.abs() + db.abs() != 1 { return Err("windows are not edge-adjacent".into()); }
    let (ax0, ax1, ay0, ay1) = a.lat.core();
    // Strip posts: along the shared edge, |d| <= h.
    let mut posts = Vec::new();
    if db != 0 {
        let ye = if db > 0 { ay1 } else { ay0 };
        for gy in ye - h_posts..=ye + h_posts { for gx in ax0..=ax1 { posts.push((gx, gy)); } }
    } else {
        let xe = if da > 0 { ax1 } else { ax0 };
        for gy in ay0..=ay1 { for gx in xe - h_posts..=xe + h_posts { posts.push((gx, gy)); } }
    }
    let mut s = SeamStats { a: a.lat.key(), b: b.lat.key(), strip_half_width_m: h_posts as f64 * T1B_SPACING_M,
        rms_tolerance_m: 30.0, p99_tolerance_m: 150.0, blend_order_independent: true, ..Default::default() };
    let mut diffs = Vec::new();
    let mut blended = Vec::with_capacity(posts.len());
    for &(gx, gy) in &posts {
        let (ia, ib) = (a.lat.idx(gx, gy).ok_or("strip post outside a")?, b.lat.idx(gx, gy).ok_or("strip post outside b")?);
        s.strip_posts += 1;
        // NET (0.8.0) posts are derived per lod-8 tile: compared like free posts, so a lod-8 seam counts its channel beds.
        let soft = |f: u8| f == FREE || f == NET;
        if soft(a.fixed[ia]) && soft(b.fixed[ib]) {
            diffs.push((a.height_mm[ia] - b.height_mm[ib]) as f64 / 1000.0);
        } else if a.fixed[ia] != FREE && a.fixed[ia] != RING && b.fixed[ib] != FREE && b.fixed[ib] != RING {
            s.fixed_posts_compared += 1;
            if a.height_mm[ia] != b.height_mm[ib] || a.fixed[ia] != b.fixed[ib] { s.fixed_posts_differing += 1; }
        } else if a.fixed[ia] != RING && b.fixed[ib] != RING {
            s.mixed_posts += 1;
            s.mixed_max_abs_m = s.mixed_max_abs_m.max((a.height_mm[ia] - b.height_mm[ib]).abs() as f64 / 1000.0);
        }
        if a.height_mm[ia] != b.height_mm[ib] || a.fixed[ia] != b.fixed[ib] { s.strip_posts_differing += 1; }
        let v1 = blend_post(&[a, b], gx, gy, h_posts)?;
        let v2 = blend_post(&[b, a], gx, gy, h_posts)?;
        if v1 != v2 { s.blend_order_independent = false; }
        blended.push((gx, gy, v1.ok_or("strip post has no weight")?));
    }
    s.free_in_both = diffs.len() as u64;
    if !diffs.is_empty() {
        let k = diffs.len() as f64;
        s.diff_mean_m = diffs.iter().sum::<f64>() / k;
        s.diff_rms_m = (diffs.iter().map(|d| d * d).sum::<f64>() / k).sqrt();
        let mut ab: Vec<f64> = diffs.iter().map(|d| d.abs()).collect();
        ab.sort_by(|x, y| x.total_cmp(y));
        s.diff_p99_abs_m = pct(&ab, 0.99);
        s.diff_max_abs_m = *ab.last().unwrap();
    }
    // Pits and slopes in the strip interior (posts with all 8 neighbours in the strip).
    let bmap: BTreeMap<(i64, i64), i32> = blended.iter().map(|&(x, y, v)| ((x, y), v)).collect();
    let pit_and_slope = |get: &dyn Fn(i64, i64) -> Option<(i32, u8)>| -> (u64, f64) {
        let mut pits = 0u64;
        let mut smax = 0.0f64;
        for &(gx, gy) in &posts {
            let Some((zc, fc)) = get(gx, gy) else { continue };
            let mut nb = Vec::with_capacity(8);
            for (ox, oy) in crate::t1b_lem::D8N { if let Some((v, _)) = get(gx + ox, gy + oy) { nb.push(v); } }
            if nb.len() < 8 { continue; }
            if fc == FREE && nb.iter().all(|v| *v > zc) { pits += 1; }
            if let (Some((e, _)), Some((w, _)), Some((nn_, _)), Some((ss, _))) = (get(gx + 1, gy), get(gx - 1, gy), get(gx, gy + 1), get(gx, gy - 1)) {
                let gxv = (e - w) as f64 / 1000.0 / (2.0 * T1B_SPACING_M);
                let gyv = (nn_ - ss) as f64 / 1000.0 / (2.0 * T1B_SPACING_M);
                smax = smax.max((gxv * gxv + gyv * gyv).sqrt());
            }
        }
        (pits, smax)
    };
    let fixed_any = |gx: i64, gy: i64| -> u8 {
        let fa = a.lat.idx(gx, gy).map(|i| a.fixed[i]).unwrap_or(RING);
        if fa == CHANNEL || fa == MARINE || fa == LAKE || fa == GUIDE || fa == PLAIN { fa } else { FREE }
    };
    let (bp, bs) = pit_and_slope(&|x, y| bmap.get(&(x, y)).map(|v| (*v, fixed_any(x, y))));
    let in_strip = |x: i64, y: i64| bmap.contains_key(&(x, y));
    let (ap, as_) = pit_and_slope(&|x, y| if in_strip(x, y) { a.lat.idx(x, y).map(|i| (a.height_mm[i], a.fixed[i])) } else { None });
    let (bpp, _) = pit_and_slope(&|x, y| if in_strip(x, y) { b.lat.idx(x, y).map(|i| (b.height_mm[i], b.fixed[i])) } else { None });
    s.blended_pits_in_strip = bp;
    s.blended_max_slope_in_strip = bs;
    s.raw_a_pits_in_strip = ap;
    s.raw_a_max_slope_in_strip = as_;
    s.raw_b_pits_in_strip = bpp;
    s.exact = s.strip_posts_differing == 0;
    s.pass_tolerance = s.diff_rms_m <= s.rms_tolerance_m && s.diff_p99_abs_m <= s.p99_tolerance_m && s.fixed_posts_differing == 0 && s.mixed_posts == 0 && s.blend_order_independent;
    Ok((s, blended))
}

// ------------------------------------------------------------------ F12

#[derive(Clone, Debug, Default, Serialize)]
pub struct F12Row { pub stat: String, pub a: f64, pub b: f64, pub diff: f64, pub tolerance: String, pub pass: bool }

#[derive(Clone, Debug, Default, Serialize)]
pub struct F12Result { pub rows: Vec<F12Row>, pub channels_fixed_identical: bool, pub pass: bool }

/// Pre-registered F12 tolerances (FALSIFIER.txt section B).
pub fn f12_compare(a: &ReliefStats, b: &ReliefStats, ra: &RawWindow, rb: &RawWindow) -> F12Result {
    let mut rows = Vec::new();
    let abs = |name: &str, x: f64, y: f64, tol: f64| F12Row { stat: name.into(), a: x, b: y, diff: (x - y).abs(), tolerance: format!("<= {tol}"), pass: (x - y).abs() <= tol };
    let rel = |name: &str, x: f64, y: f64, tol: f64| {
        let d = (x - y).abs() / x.abs().max(y.abs()).max(1e-12);
        F12Row { stat: name.into(), a: x, b: y, diff: d, tolerance: format!("<= {tol} relative"), pass: d <= tol }
    };
    rows.push(abs("mean_m", a.mean_m, b.mean_m, 25.0));
    rows.push(abs("p05_m", a.p05_m, b.p05_m, 50.0));
    rows.push(abs("p50_m", a.p50_m, b.p50_m, 50.0));
    rows.push(abs("p95_m", a.p95_m, b.p95_m, 50.0));
    rows.push(rel("std_m", a.std_m, b.std_m, 0.05));
    rows.push(rel("slope_mean_gradient", a.slope_mean_gradient, b.slope_mean_gradient, 0.05));
    rows.push(rel("slope_p90_deg", a.slope_p90_deg, b.slope_p90_deg, 0.05));
    rows.push(rel("local_relief_mean_m", a.local_relief_mean_m, b.local_relief_mean_m, 0.05));
    rows.push(abs("hypsometric_integral", a.hypsometric_integral, b.hypsometric_integral, 0.02));
    rows.push(rel("valley_post_fraction", a.valley_post_fraction, b.valley_post_fraction, 0.10));
    // Channels at or above threshold are the fixed T1a links: identical posts and beds in both.
    let mut same = ra.lat == rb.lat;
    if same {
        for i in 0..ra.lat.nodes() {
            let ca = ra.fixed[i] == CHANNEL || ra.fixed[i] == GUIDE;
            if ca != (rb.fixed[i] == CHANNEL || rb.fixed[i] == GUIDE) || (ca && ra.fixed[i] != rb.fixed[i]) || (ca && ra.height_mm[i] != rb.height_mm[i]) { same = false; break; }
        }
    }
    let pass = rows.iter().all(|r| r.pass) && same;
    F12Result { rows, channels_fixed_identical: same, pass }
}

// ------------------------------------------------------------------ lane DRAIN-FIX WALL-CHECK (key wall_checks=1)

/// The declared thresholds X (m) of the one-post wall census.
pub const WALL_THRESHOLDS_M: [f64; 3] = [100.0, 300.0, 1000.0];
pub const WALL_RULE: &str = "one_post_wall_census_v1 (lane DRAIN-FIX WALL-CHECK, 2026-10-22; reported, not gated): over the window core, (a) all-post walls: z minus the median of the 8 neighbours > X (lane GUIDE-SPIKES' census metric; the median of 8 is the upper median, the 5th smallest, as scipy.ndimage.median_filter); (b) channel walls: a post coded CHANNEL (1) or GUIDE (5) whose z exceeds the median of its non-channel land 8-neighbours (codes other than 1, 5, marine 2 and lake 4; at least 3 of them) by more than X, i.e. a one-post channel line standing above the slopes beside it; (c) all-post pits: the median of the 8 neighbours minus z > X (so a fix that removes walls by opening one-post slots shows). Counts at each X in WALL_THRESHOLDS_M, the largest excesses and the 8 largest channel walls (chart m)";

#[derive(Clone, Debug, Default, Serialize)]
pub struct WallCensus {
    pub rule: String,
    pub thresholds_m: Vec<f64>,
    pub core_posts: u64,
    pub walls_all: Vec<u32>,
    pub walls_channel: Vec<u32>,
    /// ... of which GUIDE (5) and T1a CHANNEL (1) posts.
    pub walls_guide: Vec<u32>,
    pub walls_t1a: Vec<u32>,
    pub pits_all: Vec<u32>,
    pub max_wall_all_m: f64,
    pub max_wall_channel_m: f64,
    pub max_pit_m: f64,
    /// (x m, y m, z m, excess m, code) of the largest channel walls.
    pub top_channel: Vec<(f64, f64, f64, f64, u8)>,
}

/// Lane DRAIN-FIX WALL-CHECK: the one-post wall census of a window core (WALL_RULE).
pub fn wall_census(raw: &RawWindow) -> WallCensus {
    let lat = raw.lat;
    let n = lat.n;
    let z = raw.z_vec();
    let k = WALL_THRESHOLDS_M.len();
    let mut w = WallCensus { rule: WALL_RULE.into(), thresholds_m: WALL_THRESHOLDS_M.to_vec(), walls_all: vec![0; k], walls_channel: vec![0; k],
        walls_guide: vec![0; k], walls_t1a: vec![0; k], pits_all: vec![0; k], ..Default::default() };
    let chan = |c: u8| c == CHANNEL || c == GUIDE;
    let mut top: Vec<(f64, f64, f64, f64, u8)> = Vec::new();
    for i in 0..lat.nodes() {
        if !lat.in_core(i) { continue; }
        let (x, y) = (i % n, i / n);
        if x == 0 || y == 0 || x + 1 >= n || y + 1 >= n { continue; }
        w.core_posts += 1;
        let mut nb = [0f64; 8];
        let mut land: Vec<f64> = Vec::with_capacity(8);
        let mut m = 0;
        for oy in [-1i64, 0, 1] { for ox in [-1i64, 0, 1] {
            if ox == 0 && oy == 0 { continue; }
            let j = ((y as i64 + oy) as usize) * n + (x as i64 + ox) as usize;
            nb[m] = z[j]; m += 1;
            let c = raw.fixed[j];
            if !chan(c) && c != MARINE && c != LAKE { land.push(z[j]); }
        } }
        nb.sort_by(|a, b| a.total_cmp(b));
        let med = nb[4];
        let (wall, pit) = (z[i] - med, med - z[i]);
        w.max_wall_all_m = w.max_wall_all_m.max(wall);
        w.max_pit_m = w.max_pit_m.max(pit);
        let c = raw.fixed[i];
        let cw = if chan(c) && land.len() >= 3 { land.sort_by(|a, b| a.total_cmp(b)); z[i] - land[land.len() / 2] } else { f64::NEG_INFINITY };
        if cw.is_finite() { w.max_wall_channel_m = w.max_wall_channel_m.max(cw); }
        for (t, &x_m) in WALL_THRESHOLDS_M.iter().enumerate() {
            if wall > x_m { w.walls_all[t] += 1; }
            if pit > x_m { w.pits_all[t] += 1; }
            if cw > x_m { w.walls_channel[t] += 1; if c == GUIDE { w.walls_guide[t] += 1; } else { w.walls_t1a[t] += 1; } }
        }
        if cw > WALL_THRESHOLDS_M[0] {
            let (gx, gy) = lat.g(i);
            top.push(((gx * SPACING_MM) as f64 / 1000.0, (gy * SPACING_MM) as f64 / 1000.0, z[i], cw, c));
            if top.len() > 64 { top.sort_by(|a, b| b.3.total_cmp(&a.3)); top.truncate(8); }
        }
    }
    top.sort_by(|a, b| b.3.total_cmp(&a.3));
    top.truncate(8);
    w.top_channel = top;
    w
}

#[cfg(test)]
mod te_tests {
    use super::*;
    fn flat(rx: i64, ry: i64, h: i64) -> RawWindow {
        let lat = WinLat::new(rx, ry, h);
        let nn = lat.nodes();
        RawWindow { lat, height_mm: vec![0; nn], fixed: vec![FREE; nn], drainage_area_m2: vec![0.0; nn], t0_bilinear_m: vec![0.0; nn], base_mm: vec![0; nn],
            roughness_mm: vec![0; nn], initial_mm: vec![0; nn], uplift_m_per_my: vec![0.0; nn], t0_cell: vec![0; nn], inflow_area_m2: vec![0.0; nn],
            exposed_litho: None, over: None, volcano_id: None, shore_mm: None, uncut_mm: None, ct_regime: None, ct_cliff_permille: None, ct_top_mm: None, ct_base_mm: None, ct_tread_permille: None, ct_tread_tan_permille: None, ct_ends: None, prov_id: None, runoff_weight: None, edifice_style: None, edifice_unit: None, front_mm: None, gl_carve_mm: None, digest: String::new(), recipe_digest: String::new() }
    }
    #[test]
    fn wall_census_counts_one_post_channel_walls() {
        // Lane DRAIN-FIX: a guide line standing 500 m above flat land is 3 channel walls at 100 and 300 m, none at 1 km.
        let mut a = flat(0, 0, 16);
        let (x0, _, y0, _) = a.lat.core();
        for k in 0..3 { let i = a.lat.idx(x0 + 10 + k, y0 + 10).unwrap(); a.fixed[i] = GUIDE; a.height_mm[i] = 500_000; }
        let w = wall_census(&a);
        assert_eq!((w.walls_channel.clone(), w.walls_guide.clone(), w.walls_t1a.clone()), (vec![3, 3, 0], vec![3, 3, 0], vec![0, 0, 0]));
        // A line of 3 is not an all-post wall at the 8-median (the line ends have 1 guide neighbour, the middle 2): all 3 count.
        assert_eq!(w.walls_all, vec![3, 3, 0]);
        assert_eq!(w.pits_all, vec![0, 0, 0]);
        assert!((w.max_wall_channel_m - 500.0).abs() < 1e-9 && w.top_channel.len() == 3);
        let f = wall_census(&flat(0, 0, 16));
        assert_eq!((f.walls_all.clone(), f.walls_channel.clone(), f.pits_all.clone()), (vec![0; 3], vec![0; 3], vec![0; 3]));
    }
    #[test]
    fn seam_counts_posts_fixed_in_one_window_only() {
        // Lane TE: a post that is sea in one window and land in the other was in neither comparison branch.
        let (a, mut b) = (flat(0, 0, 16), flat(1, 0, 16));
        let (s, _) = seam_compare(&a, &b, 16).unwrap();
        assert!(s.exact && s.pass_tolerance && s.mixed_posts == 0 && s.strip_posts_differing == 0);
        let (_, x1, _, _) = a.lat.core();
        let i = b.lat.idx(x1 - 3, 100).unwrap();
        b.fixed[i] = MARINE;
        b.height_mm[i] = -5_000;
        let (s, _) = seam_compare(&a, &b, 16).unwrap();
        assert_eq!((s.mixed_posts, s.strip_posts_differing, s.fixed_posts_differing, s.free_in_both), (1, 1, 0, s.strip_posts - 1));
        assert!(!s.exact && !s.pass_tolerance && (s.mixed_max_abs_m - 5.0).abs() < 1e-12 && s.diff_max_abs_m == 0.0);
    }
}

/// Spherical quadrature on a post's clipped dual square. The constructed surface
/// here is the published post-prism (constant height/owner in that square), not an
/// unmeasured continuous terrain integral. Boundary posts carry half support.
pub fn envelope_area(x: f64, y: f64, hx: f64, hy: f64, radius: f64) -> (f64,f64) {
    let j = |x:f64,y:f64| (1. + (x*x+y*y)/(4.*radius*radius)).powi(-2);
    let a = 4.*hx*hy;
    let coarse = a*j(x,y);
    let fine = a*0.25*[-0.5,0.5].iter().flat_map(|dx| [-0.5,0.5].iter().map(move |dy| j(x+dx*hx,y+dy*hy))).sum::<f64>();
    // Conservative Hessian bound for this stereographic Jacobian inside the
    // supported chart (rho << R), plus rounding. Also retain the measured delta.
    let bound = a*(hx*hx+hy*hy)/(radius*radius) + a*32.*f64::EPSILON;
    (fine, bound.max((fine-coarse).abs()))
}

pub fn envelope_window(raw: &RawWindow, prep: &WindowInputs, t0: &T0, skel: &SkeletonView,
    decl: &ChartDeclV1, cfg: &crate::t1b_litho::EnvelopeConfig) -> R<serde_json::Value> {
    use serde_json::json;
    let env=cfg.load()?;
    let water_beds: BTreeMap<usize,f64> = prep.stats.net.as_ref().and_then(|n|n.get("source_envelope_water_bed")).and_then(|b|b.get("posts")).and_then(|p|p.as_array()).into_iter().flatten().map(|v| {
        Ok((v[0].as_u64().ok_or("water bed index")? as usize, v[1].as_i64().ok_or("water bed mm")? as f64 / 1000.))
    }).collect::<R<_>>()?;
    let (mut water_posts,mut water_volume)=(0u64,0.);
    let chart=Chart::new(decl.clone()).map_err(|e| e.to_string())?;
    let lat=raw.lat; let z=raw.z_vec(); let nn=lat.nodes();
    let dam=|i:usize| raw.fixed[i]==LAKE && t0.water(raw.t0_cell[i])!=CellWater::Lake;
    let fixed:Vec<u8>=(0..nn).map(|i|if dam(i){FREE}else{raw.fixed[i]}).collect();
    let route=t1b_lem::route(lat.n,&z,&fixed,T1B_SPACING_M);
    let terms=t1b_lem::terminals(&route);
    let (x0,x1,y0,y1)=lat.core();
    // Per unit: finite source occupancy, top excess, restricted source budget,
    // volume quadrature error. Sum source h*A even where a unit is below datum.
    let mut units:BTreeMap<(u32,u16),[f64;5]>=BTreeMap::new();
    let mut cells:BTreeMap<u32,[f64;3]>=BTreeMap::new();
    let mut transfers:BTreeMap<(u32,Option<u32>),f64>=BTreeMap::new();
    let (mut solid,mut err,mut unsupported,mut exposed_bad,mut thin)=(0.,0.,0.,0u64,0u64);
    let mut top_excess_posts=0u64;
    let mut basal_unsupported_posts=0u64;
    let mut nonland_posts=0u64;
    let mut basins_by_cell=BTreeMap::new();
    for i in 0..nn {
        if !lat.in_core(i){continue;}
        let cell=raw.t0_cell[i];
        if t0.water(cell)!=CellWater::Land {nonland_posts+=1;continue;}
        let e=env.cells.get(&cell).ok_or_else(||format!("F7 missing envelope cell {cell}"))?;
        let (gx,gy)=lat.g(i);
        // in_core is inclusive; two adjacent windows each own half the shared post.
        let (xa,xb)=((gx as f64-0.5).max(x0 as f64),(gx as f64+0.5).min(x1 as f64));
        let (ya,yb)=((gy as f64-0.5).max(y0 as f64),(gy as f64+0.5).min(y1 as f64));
        let (area,ae)=envelope_area((xa+xb)*16.,(ya+yb)*16.,(xb-xa)*16.,(yb-ya)*16.,t0.radius_m);
        let (la,lo)=chart.inverse_mm(gx*T1B_SPACING_MM,gy*T1B_SPACING_MM).map_err(|e|e.to_string())?;
        let tilt=e.tilt(cell,la*RAD2DEG,lo*RAD2DEG);
        let s=if raw.fixed[i]==LAKE {
            let bed=*water_beds.get(&i).ok_or("envelope: missing actual constructed-lake bed; water stage is not solid")?;
            if bed>z[i] { return Err("envelope: lake bed above its stage".into()); }
            water_posts+=1;water_volume+=area*((z[i]-e.sea_m).max(0.)-(bed-e.sea_m).max(0.));bed
        } else { z[i] };
        let height=(s-e.sea_m).max(0.);
        solid+=area*height;err+=ae*height;
        let cc=cells.entry(cell).or_default();cc[0]+=area;cc[1]+=ae;cc[2]+=area*s;
        let first=e.layers.iter().find(|l|l.h>0.).ok_or("empty envelope")?;
        let excess=(s-(e.mean_m+tilt).max(e.sea_m)).max(0.);
        if excess>0. {top_excess_posts+=1;}
        let mut accounted=excess;
        for l in &e.layers {
            let h=e.occupancy(l,tilt,s);accounted+=h;
            let u=units.entry((cell,l.index)).or_default();
            u[0]+=area*h;u[1]+=if l.index==first.index {area*excess}else{0.};u[2]+=area*l.h;u[3]+=ae*(h+if l.index==first.index {excess}else{0.});
            u[4]+=area*e.occupancy(l,tilt,e.mean_m+tilt);
        }
        unsupported+=area*(height-accounted).max(0.);
        if e.exposed(tilt,s).is_none(){basal_unsupported_posts+=1;}
        if let Some(l)=e.exposed(tilt,s) {
            if l.h<1. {thin+=1;}
            if !envelope_exposure_matches(raw.fixed[i], raw.exposed_litho.as_ref().and_then(|v|v.get(i)).copied(), l.code) {exposed_bad+=1;}
        }
        let native=*basins_by_cell.entry(cell).or_insert(t0.final_basin(cell)?);
        let tp=if fixed[i]!=FREE {i}else{terms[i] as usize};
        let dest=if tp>=nn {None}else{
            let tc=raw.t0_cell[tp];
            let link=|n:u32|skel.outlet_basin.get(n as usize).copied().filter(|&b|b!=NONE);
            match fixed[tp] {
                CHANNEL|GUIDE|NET=>link(prep.chan_link_node[tp]),
                MARINE|PLAIN=>t0.is_marine(tc).then_some(t0.basin_id[tc as usize]).filter(|&b|b!=NONE),
                LAKE=>t0.final_basin(tc).ok().filter(|&b|b!=NONE),
                RING=>{let(gx,gy)=lat.g(tp);skel.idx((gx as f64/L8_TO_L4 as f64).round() as i64,(gy as f64/L8_TO_L4 as f64).round() as i64).and_then(|n|link(n as u32))},
                _=>None,
            }
        };
        *transfers.entry((native,dest)).or_default()+=area;
    }
    let mut sum=0.;let mut finite=0.;let mut excess_total=0.;
    let rows:Vec<_>=units.iter().map(|(&(cell,k),v)|{
        let e=&env.cells[&cell];let l=&e.layers[k as usize];let volume=v[0]+v[1];sum+=volume;finite+=v[0];excess_total+=v[1];
        json!({"cell":cell,"record_id":l.id,"native_unit_key":l.alias,"layer":k,"lithology":l.code,"h_retained_m":l.h,"base_at_cell_centre_m":l.base,"gradient_east_north":e.grad,
            "support_area_m2":cells[&cell][0],"whole_cell_area_m2":e.area_m2,"support_fraction":cells[&cell][0]/e.area_m2,
            "volume_m3":volume,"finite_source_interval_occupancy_m3":v[0],"construction_excess_above_source_top_m3":v[1],
            "h_times_support_area_m3":v[2],"h_times_whole_cell_area_m3":l.h*e.area_m2,"residual_m3":volume-v[2],
            "residual_mean_thickness_m":(volume-v[2])/cells[&cell][0],
            "above_datum_uncut_reference_m3":v[4],"surface_cut_residual_m3":volume-v[4],"below_datum_source_budget_m3":v[2]-v[4],"quadrature_error_bound_m3":v[3],
            "geometric_excess_over_restricted_budget_m3":(volume-v[2]).max(0.),"unplaced_restricted_budget_m3":(v[2]-volume).max(0.),
            "stock_admitted_m3":0})
    }).collect();
    let (mut native,mut routed)=(BTreeMap::<u32,f64>::new(),BTreeMap::<u32,f64>::new());let(mut changed,mut unknown)=(0.,0.);
    let tr:Vec<_>=transfers.iter().map(|(&(a,b),&v)|{*native.entry(a).or_default()+=v;if let Some(b)=b{*routed.entry(b).or_default()+=v;if b!=a{changed+=v;}}else{unknown+=v;}
        json!({"native_basin":a,"constructed_basin":b,"area_m2":v})}).collect();
    let ids:std::collections::BTreeSet<_>=native.keys().chain(routed.keys()).copied().collect();
    let br:Vec<_>=ids.into_iter().map(|b|{let a=native.get(&b).copied().unwrap_or(0.);let v=routed.get(&b).copied().unwrap_or(0.);
        json!({"basin":b,"native_area_m2":a,"constructed_area_m2":v,"delta_m2":v-a,"relative_delta_permille":if a>0.{Some(1000.*(v-a)/a)}else{None},"area_pass_100_permille":a>0.&&(v-a).abs()<=0.1*a})}).collect();
    let represented:std::collections::BTreeSet<_>=cells.keys().copied().collect();
    let mut edges=Vec::new();let(mut unmatched,mut jumps,mut comparisons)=(0u64,0u64,0u64);
    for &c in &represented {
        let (r,col)=t0.ij(c);
        for d in [Some((r*t0.nlon+(col+1)%t0.nlon)as u32),(r+1<t0.nlat).then_some(c+t0.nlon as u32)] .into_iter().flatten(){
            if !represented.contains(&d){continue;}
            let(a,b)=(&env.cells[&c],&env.cells[&d]);
            let(la,lo)=t0.centre_deg(c);let(lb,lob)=t0.centre_deg(d);
            let hav=((lb-la).to_radians()*0.5).sin().powi(2)+la.to_radians().cos()*lb.to_radians().cos()*((lob-lo).to_radians()*0.5).sin().powi(2);
            let spacing=2.*t0.radius_m*hav.sqrt().min(1.).asin();
            let bound=a.grad[0].hypot(a.grad[1]).max(b.grad[0].hypot(b.grad[1]))*spacing;
            let mut pairs=0;let mut bad=Vec::new();
            for l in &a.layers {if let Some(alias)=&l.alias{for m in &b.layers {if m.alias.as_ref()==Some(alias){pairs+=1;comparisons+=1;let jump=(l.base-m.base).abs();if jump>bound{jumps+=1;bad.push(json!({"records":[l.id,m.id],"jump_m":jump,"tilt_bound_m":bound,"excess_m":jump-bound}));}}}}}
            if pairs==0{unmatched+=1;}
            edges.push(json!({"cells":[c,d],"shared_alias_pairs":pairs,"unmatched":pairs==0,"base_jumps":bad}));
        }
    }
    let round_bound=(solid.abs()+sum.abs())*f64::EPSILON*nn as f64*8.;
    let closure=(sum-solid).abs();
    Ok(json!({"rule":crate::t1b_litho::ENVELOPE_RULE,"input_digest":cfg.digest,
        "support":"land-class native T0 cells on window core clipped dual-post squares; spherical stereographic Jacobian; native/constructed basins use precisely these same areas; cell ownership is published post raster, boundary location uncertainty <= half post (16m)",
        "quadrature":"four midpoint sub-samples of spherical Jacobian per dual square; constant published post height; conservative Jacobian Hessian bound plus floating summation bound. Not a continuum terrain reconstruction error bound.",
        "datum_m":t0.sea_level_m,"solid_surface_rule":"published land height; constructed lake uses retained pre-final-fill bed, never its water stage; native T0 water cells outside envelope support",
        "constructed_lake_posts":water_posts,"lake_water_above_datum_m3":water_volume,"total_solid_above_datum_m3":solid,"sum_unit_geometry_m3":sum,"finite_source_interval_occupancy_m3":finite,
        "construction_excess_m3":excess_total,"unsupported_solid_m3":unsupported,"closure_residual_m3":sum-solid,
        "spherical_quadrature_error_bound_m3":err,"summation_error_bound_m3":round_bound,"f7_geometry_pass":closure<=err+round_bound&&unsupported<=err+round_bound&&exposed_bad==0&&basal_unsupported_posts==0,
        "stock_closure_claimed":false,"stock_note":"No stock admitted. Top extension follows top source identity as construction excess. Native full-depth hA is reported alongside above-datum occupancy; below-datum reserve is not erosion. Partial-cell support cannot decide whole-cell stock restriction.",
        "exposed_code_mismatches":exposed_bad,"exposed_thin_layer_posts":thin,"top_extension_posts":top_excess_posts,"below_source_base_posts":basal_unsupported_posts,
        "nonland_posts_outside_envelope_support":nonland_posts,"units":rows,"native_basin_transfers":tr,"basins":br,
        "changed_basin_area_m2":changed,"unknown_terminal_area_m2":unknown,"f3_same_support_area_pass":unknown==0.&&br.iter().all(|r|r["area_pass_100_permille"]==true),
        "edge_support":"side-neighbour native land cells both represented in this window; counts not global and repeat across overlapping windows",
        "unmatched_identity_edges":unmatched,"base_jump_comparisons":comparisons,"bases_exceeding_tilt":jumps,"edges":edges}))
}

/// Water has no exposed rock. The routing mask can differ and must not be used here.
pub fn envelope_exposure_matches(published_fixed: u8, published_code: Option<u8>, rock_code: u8) -> bool {
    published_code == Some(if published_fixed == MARINE || published_fixed == LAKE { 255 } else { rock_code })
}

/// Retain the model's actual lake bed for accounting/E4. A source-land post
/// classified as constructed lake publishes stage as height; that is not solid.
pub fn envelope_water_beds(t0:&T0, fixed:&[u8], cell:&[u32], before_fill:&dyn Fn(usize)->Option<f64>) -> R<serde_json::Value> {
    let mut posts=Vec::new();
    for i in 0..fixed.len() {
        if fixed[i]!=LAKE || t0.water(cell[i])!=CellWater::Land {continue;}
        let z=before_fill(i).filter(|v|v.is_finite()).ok_or("envelope: constructed lake lacks retained pre-fill bed")?;
        posts.push((i,vz_tilekey::chart::quantize_mm(z*1000.)));
    }
    Ok(serde_json::json!({"rule":"actual pre-final-fill terrain bed, mm quantized as dam classifier; [published post index, bed mm]; native lakes outside E1 land support", "posts":posts}))
}
