//! The 32 m tier: zonation on a T1b window lattice. Every quantity is a function of the post's own inputs and of a
//! stencil of at most `STENCIL` posts, evaluated on the whole published lattice (core + halo); a post is valid when its
//! stencil lies inside the lattice. The published product is the window core (always valid: halo >> stencil), so
//! adjacent windows give identical values on shared posts whenever their published inputs agree there.

use crate::climate::Climate;
use crate::params::SoilParams;
use crate::window::{Volc, Window, POST_M};
use crate::zonation::{age_factor, insolation_ratio, post_t_gs, zone, PostIn, PostOut};
use vz_tilekey::chart::{Chart, ChartDeclV1, RAD2DEG};

pub const FIXED_CHANNEL_T1A: u8 = 1;
pub const FIXED_MARINE: u8 = 2;
pub const FIXED_LAKE: u8 = 4;
pub const FIXED_CHANNEL_NET: u8 = 5;
pub const OVER_EDIFICE: u8 = 1;
pub const OVER_BEACH: u8 = 3;

/// Base colonisation time scale by lane VE edifice style.
pub fn tau_of_style(p: &SoilParams, style: u8) -> f64 {
    match style {
        3 | 11 => p.tau_cinder_yr,
        4 | 10 => p.tau_tephra_yr,
        5 => p.tau_lahar_yr,
        6 | 7 | 14 => p.tau_dome_yr,
        _ => p.tau_lava_yr, // 1 lava, 2 lobe front, 8 shield, 9 fissure, 12 composite body, 13 plateau
    }
}

/// Vent and tephra kill factor at chart (x, y) m (1 = no kill).
pub fn kill_factor(p: &SoilParams, vs: &[Volc], x: f64, y: f64, t_gs: f64, ai: f64) -> f64 {
    let mut f = 1.0;
    for v in vs {
        if !v.vent[0].is_finite() { continue; }
        let r = (x - v.vent[0]).hypot(y - v.vent[1]);
        let active = v.youngest_yr() < p.active_age_yr;
        if active { f *= crate::zonation::smoothstep(v.crater_r_m, v.crater_r_m + p.vent_ring_m, r); }
        for u in &v.units {
            if u.vei < 3 { continue; }
            let rk = match u.vei { 3 => p.kill_r_vei3_m, 4 => p.kill_r_vei4_m, 5 => p.kill_r_vei5_m, _ => p.kill_r_vei6_m };
            if r >= rk { continue; }
            let inside = 1.0 - crate::zonation::smoothstep(0.5 * rk, rk, r);
            let rec = age_factor(p, u.age_yr, p.tau_tephra_yr, t_gs, ai);
            f *= 1.0 - (1.0 - rec) * inside;
        }
    }
    f
}

/// Run `f(k)` for k in 0..n on `threads` scoped threads; results in index order (each post independent: deterministic).
pub fn par_map<T: Send + Default + Clone, F: Fn(usize) -> T + Sync>(n: usize, threads: usize, f: F) -> Vec<T> {
    let mut out = vec![T::default(); n];
    let th = threads.max(1);
    let chunk = n.div_ceil(th).max(1);
    std::thread::scope(|s| {
        for (c, sl) in out.chunks_mut(chunk).enumerate() {
            let f = &f;
            s.spawn(move || { for (k, o) in sl.iter_mut().enumerate() { *o = f(c * chunk + k); } });
        }
    });
    out
}

pub const STENCIL: usize = 2;

pub struct Tier32 {
    pub n: usize,
    pub gx0: i64,
    pub gy0: i64,
    pub valid: Vec<bool>,
    pub twi: Vec<f64>,
    pub r_sun: Vec<f64>,
    pub slope_deg: Vec<f64>,
    pub f_age: Vec<f64>,
    pub soil_age_f: Vec<f64>,
    pub lake_stage: Vec<f64>,
    pub beach: Vec<f64>,
    pub shore_m: Option<Vec<f64>>,
    pub out: Vec<PostOut>,
}

/// Climate inputs at a chart point (lat/lon from the chart inverse).
pub fn climate_at(cl: &Climate, chart: &Chart, x_mm: i64, y_mm: i64) -> Result<(f64, f64, f64, f64, f64), String> {
    let (la, lo) = chart.inverse_mm(x_mm, y_mm).map_err(|e| format!("chart inverse: {e:?}"))?;
    let (lad, lod) = (la * RAD2DEG, lo * RAD2DEG);
    Ok((la, cl.bilinear(&cl.t_sl, lad, lod), cl.bilinear(&cl.half, lad, lod), cl.bilinear(&cl.p, lad, lod), cl.bilinear(&cl.ai, lad, lod)))
}

pub fn chart_of(w: &Window) -> Result<Chart, String> {
    Chart::new(ChartDeclV1::new(w.origin_lat_e7, w.origin_lon_e7, w.radius_mm)).map_err(|e| format!("chart: {e:?}"))
}

pub fn run(p: &SoilParams, w: &Window, cl: &Climate, threads: usize) -> Result<Tier32, String> {
    let n = w.n;
    let nn = n * n;
    let chart = chart_of(w)?;
    let pm = (POST_M * 1000.0) as i64;
    let valid: Vec<bool> = (0..nn).map(|k| { let (i, j) = (k % n, k / n); i >= STENCIL && j >= STENCIL && i + STENCIL < n && j + STENCIL < n }).collect();
    // climate per post
    let clim: Vec<Option<(f64, f64, f64, f64, f64)>> = par_map(nn, threads, |k| { let (gx, gy) = w.g(k); climate_at(cl, &chart, gx * pm, gy * pm).ok() });
    if clim.iter().any(|c| c.is_none()) { return Err("chart inverse failed on a window post".into()); }
    let clim: Vec<(f64, f64, f64, f64, f64)> = clim.into_iter().map(|c| c.unwrap()).collect();
    let z = &w.height_m;
    let grad = |k: usize| -> (f64, f64) {
        if !valid[k] { return (0.0, 0.0); }
        ((z[k + 1] - z[k - 1]) / (2.0 * POST_M), (z[k + n] - z[k - n]) / (2.0 * POST_M))
    };
    let slope_deg: Vec<f64> = (0..nn).map(|k| { let (a, b) = grad(k); a.hypot(b).atan().to_degrees() }).collect();
    let r_sun: Vec<f64> = par_map(nn, threads, |k| { let (a, b) = grad(k); insolation_ratio(p, clim[k].0, a, b) });
    // TWI with the geometric mean of the effective area over the 3 x 3 posts: a stand-in for multiple-flow-direction
    // dispersion (Quinn et al. 1991, recalled); single-direction (D8-like) areas stripe on flats (seen in the first run)
    let twi: Vec<f64> = (0..nn).map(|k| {
        if !valid[k] { return f64::NAN; }
        let (a, b) = grad(k);
        let tb = a.hypot(b).max(p.twi_tan_min);
        let mut la = 0.0;
        for dj in [-(n as i64), 0, n as i64] { for di in [-1i64, 0, 1] { la += w.area[(k as i64 + dj + di) as usize].max(1.0).ln(); } }
        la / 9.0 - (POST_M * tb).ln()
    }).collect();
    let r = p.lake_marsh_posts.max(0.0) as i64;
    if r as usize > STENCIL { return Err(format!("lake_marsh_posts {r} exceeds the stencil {STENCIL}")); }
    let lake_stage: Vec<f64> = (0..nn).map(|k| {
        if !valid[k] { return f64::NAN; }
        let (i, j) = ((k % n) as i64, (k / n) as i64);
        let mut st = f64::NAN;
        for dj in -r..=r { for di in -r..=r {
            let q = ((j + dj) as usize) * n + (i + di) as usize;
            if w.fixed[q] == FIXED_LAKE && !(st >= z[q]) { st = z[q]; }
        } }
        st
    }).collect();
    let beach: Vec<f64> = w.over.iter().map(|&o| if o == OVER_BEACH { 1.0 } else { 0.0 }).collect();
    // per-post zonation
    let base_in = |k: usize| -> PostIn {
        let c = clim[k];
        let fx = w.fixed[k];
        let h = z[k];
        let lake_fringe = lake_stage[k].is_finite() && fx != FIXED_LAKE && h - lake_stage[k] > 0.0 && h - lake_stage[k] <= p.lake_marsh_m;
        PostIn { h, sea: cl.sea, water: fx == FIXED_MARINE || fx == FIXED_LAKE, t_sl: c.1, half: c.2, p: c.3, ai: c.4, lapse: cl.lapse,
            r_sun: r_sun[k], slope_deg: slope_deg[k], regime: w.ct_regime[k], reg_factor: None, twi: twi[k], f_age: 1.0, f_kill: 1.0, soil_age_f: 1.0, debris: 0.0,
            beach: w.over[k] == OVER_BEACH, riparian: fx == FIXED_CHANNEL_T1A || fx == FIXED_CHANNEL_NET, lake_fringe,
            flat_slope_deg: slope_deg[k], shore_m: w.shore_m.as_ref().map(|s| s[k]) }
    };
    let res: Vec<(f64, f64, PostOut)> = par_map(nn, threads, |k| {
        let mut x = base_in(k);
        let t_gs = post_t_gs(p, &x);
        let (mut fa, mut sa) = (1.0, 1.0);
        if w.over[k] == OVER_EDIFICE {
            let vid = w.volcano_id[k] as usize;
            if let Some(v) = w.volcanoes.get(vid) {
                let u = w.edifice_unit[k];
                let (age, style) = if u == 0 || u == 65535 || (u as usize) > v.units.len() { (v.t_h_yr, 12u8) } else {
                    let un = &v.units[u as usize - 1];
                    (un.age_yr, if w.edifice_style[k] != 0 { w.edifice_style[k] } else { un.style })
                };
                fa = age_factor(p, age, tau_of_style(p, style), t_gs, x.ai);
                sa = 1.0 - (-age.max(0.0) / p.tau_soil_yr).exp();
            }
        }
        let (gx, gy) = w.g(k);
        x.f_age = fa;
        x.soil_age_f = sa;
        x.f_kill = kill_factor(p, &w.volcanoes, gx as f64 * POST_M, gy as f64 * POST_M, t_gs, x.ai);
        (fa, sa, zone(p, &x))
    });
    let f_age = res.iter().map(|r| r.0).collect();
    let soil_age_f = res.iter().map(|r| r.1).collect();
    let out = res.into_iter().map(|r| r.2).collect();
    Ok(Tier32 { n, gx0: w.gx0, gy0: w.gy0, valid, twi, r_sun, slope_deg, f_age, soil_age_f, lake_stage, beach, shore_m: w.shore_m.clone(), out })
}

impl Tier32 {
    /// Bilinear sample of a 32 m field at chart (x, y) m; None when a corner is not valid.
    pub fn sample(&self, f: &[f64], x: f64, y: f64) -> Option<f64> {
        let (u, v) = (x / POST_M - self.gx0 as f64, y / POST_M - self.gy0 as f64);
        let (i0, j0) = (u.floor(), v.floor());
        let (tu, tv) = (u - i0, v - j0);
        let (i0, j0) = (i0 as i64, j0 as i64);
        if i0 < 0 || j0 < 0 || i0 + 1 >= self.n as i64 || j0 + 1 >= self.n as i64 { return None; }
        let k = j0 as usize * self.n + i0 as usize;
        let ks = [k, k + 1, k + self.n, k + self.n + 1];
        if ks.iter().any(|&q| !self.valid[q]) { return None; }
        Some((1.0 - tv) * ((1.0 - tu) * f[ks[0]] + tu * f[ks[1]]) + tv * ((1.0 - tu) * f[ks[2]] + tu * f[ks[3]]))
    }
    /// Max over the 4 bilinear corners of the finite values (NaN when none).
    pub fn corner_max(&self, f: &[f64], x: f64, y: f64) -> f64 {
        let (u, v) = (x / POST_M - self.gx0 as f64, y / POST_M - self.gy0 as f64);
        let (i0, j0) = (u.floor() as i64, v.floor() as i64);
        let mut m = f64::NAN;
        for (di, dj) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
            let (i, j) = (i0 + di, j0 + dj);
            if i < 0 || j < 0 || i >= self.n as i64 || j >= self.n as i64 { continue; }
            let q = f[j as usize * self.n + i as usize];
            if q.is_finite() && !(m >= q) { m = q; }
        }
        m
    }
}
