//! T1b 0.9.4 lane TG (2026-10-15; `unpublished-design-record`): lane TF ideas 1, 2, 4, 5 on walk config C.
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon; construction (no SimTime). Everything is opt-in under
//! `params.ls.tb.tc.tg` (inside params.ls.tb, so lane TB's paired keying keeps every keyed draw); absent = lane TF's
//! behaviour byte-identical.
//!
//! - G1 `shore_mm`: a signed sub-post distance to the sea shoreline, published as an additive window field (presentation can
//!   draw the coast and cliff at sub-post precision instead of the 32 m marine-mask staircase).
//! - G2 sill grade: TF's soft_bluff_v2 rule (c) holds a post at least at the height of the post its shore distance came from;
//!   with G2 it is held at that height plus the declared valley-floor slope times the link length (no flat terrace behind a
//!   hard sill).
//! - G3 smoothed bilinear: the T0 bilinear box-averaged over 2 km on a global 128 m chart grid replaces the plain bilinear where
//!   it shapes the surface (E taper, construction water, T1a shelf grade, fan subaqueous test): C1 across T0 cell-centre lines.
//! - G4 F deadband: the flooding offset moves only by the land-fraction error outside a declared inner band.

use crate::t1b_rivernet::EvalLat;
use crate::*;
use serde::{Deserialize, Serialize};
use std::sync::{Arc, Mutex, OnceLock};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct TgParamsV1 {
    /// G1: publish shore_mm (flag) and its clamp (mm).
    pub shore: u32,
    pub shore_cap_mm: i64,
    /// G2: grade TF's T2(c) chain at the declared valley-floor slope (params.floor_slope_ppm) (flag; needs tf bluff_v2).
    pub sill_grade: u32,
    /// G3: the box-averaged bilinear (flag), its half-width and the global grid spacing it is sampled on (mm).
    pub smooth_bil: u32,
    pub smooth_half_mm: i64,
    pub smooth_grid_mm: i64,
    /// G4: F4 deadband (flag) and its half-width (permille of land fraction).
    pub f_deadband: u32,
    pub f_deadband_permille: i64,
    pub shore_rule: String,
    pub sill_rule: String,
    pub smooth_rule: String,
    pub deadband_rule: String,
}

impl TgParamsV1 {
    pub fn v1() -> Self {
        TgParamsV1 {
            shore: 1, shore_cap_mm: 2_048_000, sill_grade: 1, smooth_bil: 1, smooth_half_mm: 2_000_000, smooth_grid_mm: 128_000, f_deadband: 0, f_deadband_permille: 15,
            shore_rule: "shore_field_v1: phi per eval-lattice post after the coast stage, before marine posts are published at sea level: phi = z - sea, except posts cut to sea by the wave cut take the cut profile extended seaward, sea + talus_slope (d - platform width) (d = the coast stage's shore distance; the kept posts follow the same talus line, so a cut / kept link crosses zero at d = platform width); sign forced to the published mask (marine: phi <= -1 mm, every other post >= +1 mm; lakes are not sea). Shoreline = the zero contour of phi by marching squares (linear crossings on 4-neighbour links, saddles by the 4-corner mean: lane TF shore_segments). shore_mm = the Euclidean distance (mm) from the post to that polyline (lane TF segment vector propagation through every post), negative on marine posts, clamped to +-shore_cap; a consumer's bilinear interpolation of shore_mm has its zero isoline on the sub-post shoreline".into(),
            sill_rule: "sill_grade_v1: in soft_bluff_v2 a post ends no lower than z(parent) + s x L (parent = the post its shore distance came from, L = the link length 32 or 45.25 m, s = params.floor_slope_ppm, the declared valley-floor slope of nearest_channel_floor_v1) instead of z(parent): soft ground behind a hard sill drains toward the sill at the valley-floor slope (rising inland strictly, so no closed basin is cut)".into(),
            smooth_rule: "smooth_bilinear_v1: B = the T0 bilinear (as lat_prep: chart inverse, degrees, plain bilinear between L1 cell centres) sampled on the GLOBAL chart grid at smooth_grid (coordinates clamped to the chart), box-averaged separably with box_mean weights clamp(h - |k|, 0, 1), h = half / grid + 1/2, summed in fixed order k = -K..K, then bilinearly interpolated: one function for nodes, the solve lattice and every window (window-independent). B replaces the bilinear in the E taper tau (posts and nodes), the construction water (post water class and node K_MARINE at <= sea - shelf), the CF C6b T1a shelf grade (final reach and cap) and the fan subaqueous test; unchanged: the F4 target land fraction, the F6 reference, coastal flags, vent placement, the TC clip, the published t0_bilinear_m and lane TF's cone base".into(),
            deadband_rule: "f4_deadband_v1: the flooding-offset step uses ex - clamp(ex, -db, db) (ex = land fraction - bilinear land fraction, db = f_deadband_permille / 1000) instead of ex; the convergence test (F4 within flood_tolerance) is unchanged".into(),
        }
    }
    pub fn off() -> Self { let mut p = Self::v1(); p.shore = 0; p.sill_grade = 0; p.smooth_bil = 0; p.f_deadband = 0; p }
    pub fn validate(&self) -> R<()> {
        if self.shore_cap_mm <= 0 || self.smooth_half_mm <= 0 || self.smooth_grid_mm <= 0 || self.f_deadband_permille < 0 { return Err("tg: positive cap, half-width, grid; deadband >= 0".into()); }
        Ok(())
    }
}

pub fn tg(p: &crate::t1b_rivernet::RiverNetParamsV1) -> Option<&TgParamsV1> { crate::t1b_ls::tc(p).and_then(|t| t.tg.as_ref()) }
pub fn tg_on(p: &crate::t1b_rivernet::RiverNetParamsV1, f: fn(&TgParamsV1) -> u32) -> bool { tg(p).map(|t| f(t) == 1).unwrap_or(false) }

// ------------------------------------------------------------------ G3: smoothed bilinear on a global grid

/// B on a rectangle of the global grid (indices ix0..ix0+nx, iy0..iy0+ny at `step_mm`).
pub struct SmoothBil { pub step_mm: i64, pub ix0: i64, pub iy0: i64, pub nx: usize, pub ny: usize, pub v: Vec<f64> }

impl SmoothBil {
    /// B at chart position (mm, may be fractional), bilinear on the grid (clamped to the rectangle).
    pub fn at_mm(&self, x_mm: f64, y_mm: f64) -> f64 {
        let fx = x_mm / self.step_mm as f64 - self.ix0 as f64;
        let fy = y_mm / self.step_mm as f64 - self.iy0 as f64;
        let i = (fx.floor() as i64).clamp(0, self.nx as i64 - 2) as usize;
        let j = (fy.floor() as i64).clamp(0, self.ny as i64 - 2) as usize;
        let (tx, ty) = ((fx - i as f64).clamp(0.0, 1.0), (fy - j as f64).clamp(0.0, 1.0));
        let g = |a: usize, b: usize| self.v[b * self.nx + a];
        g(i, j) * (1.0 - tx) * (1.0 - ty) + g(i + 1, j) * tx * (1.0 - ty) + g(i, j + 1) * (1.0 - tx) * ty + g(i + 1, j + 1) * tx * ty
    }
}

/// Box weights (k = -K..K) of a half-width over a spacing (box_mean's kernel).
pub fn box_weights(half_mm: i64, step_mm: i64) -> Vec<f64> {
    let h = half_mm as f64 / step_mm as f64 + 0.5;
    let k = h.ceil() as i64;
    (-k..=k).map(|o| (h - o.abs() as f64).clamp(0.0, 1.0)).collect()
}

/// Build B over the grid rectangle covering [x0, x1] x [y0, y1] (mm). `raw(x_mm, y_mm)` = the T0 bilinear at a chart point.
pub fn build_smooth(raw: &dyn Fn(i64, i64) -> R<f64>, x0: f64, x1: f64, y0: f64, y1: f64, step_mm: i64, half_mm: i64) -> R<SmoothBil> {
    let w = box_weights(half_mm, step_mm);
    let k = (w.len() / 2) as i64;
    let ws: f64 = w.iter().sum();
    let (ix0, ix1) = ((x0 / step_mm as f64).floor() as i64 - 1, (x1 / step_mm as f64).ceil() as i64 + 1);
    let (iy0, iy1) = ((y0 / step_mm as f64).floor() as i64 - 1, (y1 / step_mm as f64).ceil() as i64 + 1);
    let (nx, ny) = ((ix1 - ix0 + 1) as usize, (iy1 - iy0 + 1) as usize);
    let (rnx, rny) = (nx + 2 * k as usize, ny + 2 * k as usize);
    let h = vz_tilekey::CHART_HALF_EXTENT_MM;
    let mut rawv = vec![0f64; rnx * rny];
    for b in 0..rny { for a in 0..rnx {
        let (gx, gy) = (ix0 - k + a as i64, iy0 - k + b as i64);
        rawv[b * rnx + a] = raw((gx * step_mm).clamp(-h, h), (gy * step_mm).clamp(-h, h))?;
    } }
    // Along x (every raw row, output columns), then along y; fixed summation order k = -K..K.
    let mut tx = vec![0f64; nx * rny];
    for b in 0..rny { for a in 0..nx {
        let mut s = 0.0;
        for (o, wk) in w.iter().enumerate() { s += wk * rawv[b * rnx + a + o]; }
        tx[b * nx + a] = s / ws;
    } }
    let mut v = vec![0f64; nx * ny];
    for b in 0..ny { for a in 0..nx {
        let mut s = 0.0;
        for (o, wk) in w.iter().enumerate() { s += wk * tx[(b + o) * nx + a]; }
        v[b * nx + a] = s / ws;
    } }
    Ok(SmoothBil { step_mm, ix0, iy0, nx, ny, v })
}

/// The chart-wide B (whole chart, cached per process: build_network and node_ctx both need it).
pub fn chart_smooth(t0: &crate::t0::T0, chart: &vz_tilekey::chart::Chart, half_mm: i64, step_mm: i64) -> R<Arc<SmoothBil>> {
    static CACHE: OnceLock<Mutex<Vec<(String, Arc<SmoothBil>)>>> = OnceLock::new();
    let o = chart.inverse_mm(0, 0).map_err(|e| e.to_string())?;
    let key = format!("{:p} {} {} {} {:e} {:e} {}", t0 as *const _, t0.nlat, t0.nlon, t0.elevation_m.len(), o.0, o.1, half_mm * 1_000_000_007 + step_mm);
    let c = CACHE.get_or_init(|| Mutex::new(Vec::new()));
    if let Some((_, b)) = c.lock().unwrap().iter().find(|(k, _)| *k == key) { return Ok(b.clone()); }
    let h = vz_tilekey::CHART_HALF_EXTENT_MM as f64;
    let b = Arc::new(build_smooth(&raw_bilinear(t0, chart), -h, h, -h, h, step_mm, half_mm)?);
    c.lock().unwrap().push((key, b.clone()));
    Ok(b)
}

/// The T0 bilinear at a chart point (mm), exactly as lat_prep computes it (before its f32 store).
pub fn raw_bilinear<'a>(t0: &'a crate::t0::T0, chart: &'a vz_tilekey::chart::Chart) -> impl Fn(i64, i64) -> R<f64> + 'a {
    move |x, y| {
        let (la, lo) = chart.inverse_mm(x, y).map_err(|e| e.to_string())?;
        let elev = |c: u32| t0.elevation_m[c as usize] as f64;
        Ok(crate::t1b::bilinear(t0, &elev, la * vz_tilekey::chart::RAD2DEG, lo * vz_tilekey::chart::RAD2DEG, false))
    }
}

// ------------------------------------------------------------------ G1: shore field

/// Signed distance (m) to the zero contour of phi = zc - sea, sign forced to the marine mask, clamped to +-cap.
pub fn shore_field(lat: &EvalLat, zc: &[f64], sea: f64, marine: &dyn Fn(usize) -> bool, cap_m: f64) -> Vec<f64> {
    let nn = lat.nodes();
    let zp: Vec<f64> = (0..nn).map(|i| {
        let ph = if zc[i].is_finite() { zc[i] - sea } else { 0.0 };
        sea + if marine(i) { ph.min(-1e-3) } else { ph.max(1e-3) }
    }).collect();
    let segs = crate::t1b_tf::shore_segments(lat, &zp, sea, marine, None);
    let v = crate::t1b_tf::vprop_segments(lat, &segs, &|_| true, cap_m);
    (0..nn).map(|i| { let d = v.d[i].min(cap_m); if marine(i) { -d } else { d } }).collect()
}

/// G4: the part of the land-fraction error that moves F.
#[inline]
pub fn deadband(ex: f64, db: f64) -> f64 { ex - ex.clamp(-db, db) }

/// Lane TG coast statistics (window header, only with tg).
#[derive(Clone, Debug, Default, Serialize)]
pub struct TgCoastStats {
    pub sill_graded_posts: u32,
    pub sill_graded_raise_max_m: f64,
    pub cut_ext_posts: u32,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn smooth_is_window_independent_and_exact_on_planes() {
        // A tent (crease at x = 0) and a plane: B equals the plane away from the crease, is C1 across it, and two rectangles
        // give bit-identical values at shared posts.
        let raw = |x: i64, y: i64| -> R<f64> { let (x, y) = (x as f64 / 1000.0, y as f64 / 1000.0); Ok(-0.1 * x.abs() + 0.02 * y) };
        let a = build_smooth(&raw, -20_000_000.0, 5_000_000.0, -3_000_000.0, 3_000_000.0, 128_000, 2_000_000).unwrap();
        let b = build_smooth(&raw, -6_000_000.0, 20_000_000.0, -1_000_000.0, 9_000_000.0, 128_000, 2_000_000).unwrap();
        for gx in -150..150 { for gy in -30..90 {
            let (x, y) = (gx as f64 * 32_000.0, gy as f64 * 32_000.0);
            if x < -5_900_000.0 || x > 4_900_000.0 || y < -900_000.0 || y > 2_900_000.0 { continue; }
            assert_eq!(a.at_mm(x, y).to_bits(), b.at_mm(x, y).to_bits(), "{gx} {gy}");
        } }
        // Far from the crease: the plane.
        let p = a.at_mm(-4_000_000.0, 1_000_000.0);
        assert!((p - (-0.1 * 4000.0 + 0.02 * 1000.0)).abs() < 1e-6, "{p}");
        // At the crease: lowered by the box mean of |x| (about h / 4 x 0.1 x 2), slope continuous (C1).
        let s = |x: f64| (a.at_mm(x + 64_000.0, 0.0) - a.at_mm(x - 64_000.0, 0.0)) / 128.0;
        assert!(s(0.0).abs() < 1e-3, "slope at crease {}", s(0.0));
        let jumps: f64 = (-20..20).map(|k| (s(k as f64 * 128_000.0 + 64_000.0) - s(k as f64 * 128_000.0 - 64_000.0)).abs()).fold(0.0, f64::max);
        assert!(jumps < 0.02, "max slope change per 128 m {jumps} (raw tent 0.2)");
    }
    #[test]
    fn shore_field_signs_and_subpost_zero() {
        // An oblique plane coast: the zero isoline of shore_mm is the plane's zero line (sub-post), signs follow the mask.
        let n = 40;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 };
        let nr = (5f64).sqrt();
        let z: Vec<f64> = (0..n * n).map(|i| { let (x, y) = ((i % n) as f64, (i / n) as f64); (2.0 * x + y - 40.3) / nr }).collect();
        let sf = shore_field(&lat, &z, 0.0, &|i| z[i] <= 0.0, 2048.0);
        for i in 0..n * n {
            let (x, y) = (i % n, i / n);
            assert_eq!(sf[i] < 0.0, z[i] <= 0.0);
            if x < 4 || y < 4 || x > n - 5 || y > n - 5 { continue; }
            let t = z[i] / nr; let (fx, fy) = (x as f64 - 2.0 * t, y as f64 - t);
            if fx < 1.0 || fy < 1.0 || fx > (n - 2) as f64 || fy > (n - 2) as f64 { continue; }
            assert!((sf[i] - z[i] * 32.0).abs() < 0.05, "post {x},{y}: {} vs {}", sf[i], z[i] * 32.0);
        }
    }
    #[test]
    fn deadband_moves_only_outside() {
        assert_eq!(deadband(0.01, 0.015), 0.0);
        assert!((deadband(0.02, 0.015) - 0.005).abs() < 1e-12);
        assert!((deadband(-0.04, 0.015) + 0.025).abs() < 1e-12);
    }
}
