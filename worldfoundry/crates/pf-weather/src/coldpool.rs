//! COLDPOOL-1: dawn cold-air pool and radiation-fog top over a T1b window.
//!
//! A derived presentation/perception layer (never canon by itself): a pure
//! function of the window's canonical heights, the planet location of every
//! post and the WEATHER-V0 record of the night. Model (design review REPORT 1.4):
//!
//! 1. Valley floor: square min filter of the 32 m heights (half-width
//!    `pool_floor_radius_m`), then a box mean (half-width
//!    `pool_smooth_radius_m`) -> a smooth floor sheet P0. Cold air drains to
//!    the lowest ground within ~1 km and pools; in steep valleys the floor
//!    a little downstream is lower, so pools there are shallower.
//! 2. Night cooling potential Q (0..1) = mean over the night (sunset to
//!    sunrise, hourly) of the WEATHER-V0 `radiative_cooling` (clear^2 * calm *
//!    dark, zero while it precipitates) at the T0 cell centres, blended
//!    bilinearly to the post. Pool depth D = pool_depth_max * Q (published only
//!    when >= pool_depth_min); inversion strength dT = pool_inversion_max * Q.
//!    Before sunrise both grow as sqrt(elapsed night fraction): the pool is
//!    deepest at sunrise.
//! 3. Fog: inside the pool the temperature rises linearly from (T_free(P) - dT)
//!    at the floor to T_free(P) at the pool top P = P0 + D. With the
//!    near-surface dewpoint depression DD at elevation P (WEATHER-V0 evaluated
//!    at that surface elevation: the T0 humidity is a near-surface climatology),
//!    the air is saturated below F = P0 + D (1 - DD / dT): that is the fog top.
//!    DD <= 0 fills the whole pool; DD >= dT gives a clear cold pool (frost
//!    pocket, inversion) with no fog.
//! 4. Burn-off: after sunrise the fog top (and the pool lid) falls at b = burn_rate_max *
//!    sin(noon sun elevation) m/h, so thin fog at the edges and in small
//!    valleys clears first and the sheet shrinks inward; the formed state is
//!    frozen at sunrise (no double count with the warming free air).
//!
//! Window locality: floor + smoothing stencil is 63 posts, inside the 64-post
//! halo, and the weather blend is a function of the post's planet location
//! only, so two windows that agree on heights agree on the product at shared
//! posts (measured, see the lane report).

use std::path::Path;

use serde::Serialize;
use vz_tilekey::chart::{Chart, ChartDeclV1};

use crate::state::WeatherModel;
use crate::sun::{solar_clock, sunrise_world_hour, sunset_world_hour, DEG};

pub const SCHEMA: &str = "pf.coldpool-window/v0";
pub const T1B_SPACING_M: f64 = 32.0;

pub struct Window {
    pub n: usize,
    pub halo: usize,
    pub gx0: i64,
    pub gy0: i64,
    pub rx: i64,
    pub ry: i64,
    pub height_m: Vec<f64>,
    pub fixed: Vec<u8>,
    pub chart: Chart,
    pub digest: String,
    pub window_key: String,
}

fn json(p: &Path) -> Result<serde_json::Value, String> {
    let s = std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
    serde_json::from_str(&s).map_err(|e| format!("{}: {e}", p.display()))
}

impl Window {
    pub fn load(dir: &Path) -> Result<Window, String> {
        let v = json(&dir.join("window.json"))?;
        let h = &v["header"];
        let lat = &h["lattice"];
        let gi = |k: &str| lat[k].as_i64().ok_or(format!("lattice.{k}"));
        let n = gi("n")? as usize;
        let ch = &h["recipe"]["chart"];
        if ch["kind"].as_str() != Some(vz_tilekey::chart::CHART_KIND_V1) {
            return Err(format!("unsupported chart kind {:?}", ch["kind"]));
        }
        let decl = ChartDeclV1::new(
            ch["origin_lat_e7"].as_i64().ok_or("origin_lat_e7")?,
            ch["origin_lon_e7"].as_i64().ok_or("origin_lon_e7")?,
            ch["datum_radius_mm"].as_i64().ok_or("datum_radius_mm")?,
        );
        let chart = Chart::new(decl).map_err(|e| format!("chart: {e:?}"))?;
        let hb = std::fs::read(dir.join("fields/height_mm.i32")).map_err(|e| e.to_string())?;
        let fb = std::fs::read(dir.join("fields/fixed.u8")).map_err(|e| e.to_string())?;
        if hb.len() != n * n * 4 || fb.len() != n * n {
            return Err("window field sizes".into());
        }
        let height_m = hb.chunks_exact(4).map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as f64 / 1000.0).collect();
        Ok(Window {
            n,
            halo: gi("halo_posts")? as usize,
            gx0: gi("gx0")?,
            gy0: gi("gy0")?,
            rx: gi("rx")?,
            ry: gi("ry")?,
            height_m,
            fixed: fb,
            chart,
            digest: v["digest"].as_str().unwrap_or("unknown").to_string(),
            window_key: h["window"].as_str().unwrap_or("").to_string(),
        })
    }

    /// Planet (lat, lon) in degrees of post i.
    pub fn post_latlon(&self, i: usize) -> (f64, f64) {
        let (gx, gy) = ((i % self.n) as i64 + self.gx0, (i / self.n) as i64 + self.gy0);
        let (la, lo) = self.chart.inverse_mm(gx * 32_000, gy * 32_000).expect("post inside chart");
        (la / DEG, lo / DEG)
    }
}

/// Sliding min over a line with half-width r (clamped at the ends), O(n).
fn line_min(src: &[f64], r: usize, out: &mut [f64]) {
    let n = src.len();
    let mut dq: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    let mut next = 0usize;
    for i in 0..n {
        let hi = (i + r).min(n - 1);
        while next <= hi {
            while let Some(&b) = dq.back() {
                if src[b] >= src[next] { dq.pop_back(); } else { break; }
            }
            dq.push_back(next);
            next += 1;
        }
        let lo = i.saturating_sub(r);
        while let Some(&f) = dq.front() {
            if f < lo { dq.pop_front(); } else { break; }
        }
        out[i] = src[*dq.front().expect("non-empty")];
    }
}

/// Box mean with half-width r (clamped window), via prefix sums.
fn line_mean(src: &[f64], r: usize, out: &mut [f64]) {
    let n = src.len();
    let mut pre = vec![0.0f64; n + 1];
    for i in 0..n { pre[i + 1] = pre[i] + src[i]; }
    for i in 0..n {
        let lo = i.saturating_sub(r);
        let hi = (i + r).min(n - 1);
        out[i] = (pre[hi + 1] - pre[lo]) / (hi + 1 - lo) as f64;
    }
}

fn separable(z: &[f64], n: usize, r: usize, f: fn(&[f64], usize, &mut [f64])) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    let mut row = vec![0.0; n];
    for y in 0..n {
        f(&z[y * n..(y + 1) * n], r, &mut row);
        a[y * n..(y + 1) * n].copy_from_slice(&row);
    }
    let mut b = vec![0.0; n * n];
    let mut col = vec![0.0; n];
    let mut colo = vec![0.0; n];
    for x in 0..n {
        for y in 0..n { col[y] = a[y * n + x]; }
        f(&col, r, &mut colo);
        for y in 0..n { b[y * n + x] = colo[y]; }
    }
    b
}

/// The weather of one night at one T0 cell centre (cell T0 elevation).
#[derive(Debug, Clone, Serialize)]
pub struct NightSummary {
    pub cell: u32,
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub z_ref_m: f64,
    pub sunset_h: f64,
    pub sunrise_h: f64,
    pub cooling_q: f64,
    pub night_hours_raining: u32,
    pub mean_wind_m_s: f64,
    pub mean_cloud: f64,
    /// Near-surface air at the cell's T0 elevation at the formation time (min(t, sunrise)).
    pub t_ref_c: f64,
    pub td_ref_c: f64,
    pub rh_ref: f64,
    pub burn_rate_m_h: f64,
    /// Near-surface (T, Td) at the formation time for surface elevations
    /// 0, 100, .. 6000 m: the T0 humidity is a near-surface climatology, so a
    /// valley at 300 m has its own dewpoint depression, not the free air of a
    /// 1600 m cell mean carried down the lapse rate.
    pub profile_t_td: Vec<[f64; 2]>,
}

pub const PROFILE_STEP_M: f64 = 100.0;
pub const PROFILE_LEVELS: usize = 61;

impl NightSummary {
    /// Dewpoint depression of near-surface air at elevation z (linear in z).
    pub fn dd_at(&self, z: f64) -> f64 {
        let f = (z.max(0.0) / PROFILE_STEP_M).min((PROFILE_LEVELS - 1) as f64);
        let i = (f.floor() as usize).min(PROFILE_LEVELS - 2);
        let w = f - i as f64;
        let a = self.profile_t_td[i];
        let b = self.profile_t_td[i + 1];
        (1.0 - w) * (a[0] - a[1]) + w * (b[0] - b[1])
    }
}

pub fn night_summary(m: &WeatherModel, cell: u32, local_day: i64, t_h: f64) -> NightSummary {
    let (lat, lon) = m.climate.cell_centre(cell);
    let tilt = m.k.axial_tilt_deg;
    let sr = sunrise_world_hour(local_day, lat, lon, tilt);
    let ss = sunset_world_hour(local_day - 1, lat, lon, tilt);
    let mut q = 0.0;
    let mut nq = 0u32;
    let mut rain = 0u32;
    let (mut w, mut cl) = (0.0, 0.0);
    let mut h = ss.ceil();
    while h <= sr {
        let s = m.state_at(lat, lon, h, None);
        q += s.radiative_cooling;
        w += s.wind_speed_m_s;
        cl += s.cloud_cover;
        if s.precip_rate_mm_h > 0.0 { rain += 1; }
        nq += 1;
        h += 1.0;
    }
    let nqf = nq.max(1) as f64;
    let tf = t_h.min(sr);
    let s = m.state_at(lat, lon, tf, None);
    let sc = solar_clock(sr, lat, lon, tilt);
    let profile_t_td = (0..PROFILE_LEVELS)
        .map(|l| {
            let p = m.state_at(lat, lon, tf, Some(l as f64 * PROFILE_STEP_M));
            [p.air_temperature_c, p.dewpoint_c]
        })
        .collect();
    NightSummary {
        cell,
        lat_deg: lat,
        lon_deg: lon,
        z_ref_m: s.elevation_m,
        sunset_h: ss,
        sunrise_h: sr,
        cooling_q: q / nqf,
        night_hours_raining: rain,
        mean_wind_m_s: w / nqf,
        mean_cloud: cl / nqf,
        t_ref_c: s.air_temperature_c,
        td_ref_c: s.dewpoint_c,
        rh_ref: s.relative_humidity,
        burn_rate_m_h: m.k.burn_rate_max_m_h * (sc.noon_elevation_deg.max(0.0) * DEG).sin(),
        profile_t_td,
    }
}

pub struct ColdPoolProduct {
    pub n: usize,
    pub fog_top_m: Vec<f32>,
    pub pool_top_m: Vec<f32>,
    /// 0 none, 1 fog, 2 cold pool without fog, 3 marine (no radiation fog)
    pub mask: Vec<u8>,
    pub cells: Vec<NightSummary>,
    pub fog_core_fraction: f64,
    pub pool_core_fraction: f64,
    pub fog_top_core_min_m: f64,
    pub fog_top_core_max_m: f64,
}

pub fn coldpool_window(m: &WeatherModel, w: &Window, local_day: i64, t_h: f64) -> ColdPoolProduct {
    let k = &m.k;
    let n = w.n;
    let rf = (k.pool_floor_radius_m / T1B_SPACING_M).round() as usize;
    let rs = (k.pool_smooth_radius_m / T1B_SPACING_M).round() as usize;
    let floor = separable(&w.height_m, n, rf, line_min);
    let p0 = separable(&floor, n, rs, line_mean);
    // Per-post planet location and the T0 cells they blend.
    let ll: Vec<(f64, f64)> = (0..n * n).map(|i| w.post_latlon(i)).collect();
    let mut cell_ids: Vec<u32> = Vec::new();
    for &(la, lo) in &ll {
        for (c, wt) in m.climate.bilinear(la, lo) {
            if wt > 0.0 && !cell_ids.contains(&c) { cell_ids.push(c); }
        }
    }
    cell_ids.sort_unstable();
    let cells: Vec<NightSummary> = cell_ids.iter().map(|&c| night_summary(m, c, local_day, t_h)).collect();
    let idx = |c: u32| cell_ids.binary_search(&c).expect("cell");
    let mut fog_top = vec![f32::NAN; n * n];
    let mut pool_top = vec![f32::NAN; n * n];
    let mut mask = vec![0u8; n * n];
    let (mut nf, mut np, mut nc) = (0usize, 0usize, 0usize);
    let (mut fmin, mut fmax) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..n * n {
        let (la, lo) = ll[i];
        let wts = m.climate.bilinear(la, lo);
        let (mut q, mut burn, mut ss, mut sr) = (0.0, 0.0, 0.0, 0.0);
        for (c, wt) in wts {
            if wt == 0.0 { continue; }
            let s = &cells[idx(c)];
            q += wt * s.cooling_q;
            burn += wt * s.burn_rate_m_h;
            ss += wt * s.sunset_h;
            sr += wt * s.sunrise_h;
        }
        let grow = if t_h >= sr { 1.0 } else { ((t_h - ss) / (sr - ss)).clamp(0.0, 1.0).sqrt() };
        let dfull = k.pool_depth_max_m * q;
        let d = if dfull >= k.pool_depth_min_m { dfull * grow } else { 0.0 };
        let dt = k.pool_inversion_max_k * q * grow;
        let z = w.height_m[i];
        let top = p0[i] + d;
        let gx = (i % n) as i64 + w.gx0;
        let gy = (i / n) as i64 + w.gy0;
        let core = gx >= w.rx * 512 && gx <= w.rx * 512 + 512 && gy >= w.ry * 512 && gy <= w.ry * 512 + 512;
        if core { nc += 1; }
        if w.fixed[i] == 2 {
            mask[i] = 3;
            continue;
        }
        if d <= 0.0 { continue; }
        // After sunrise the inversion lid and the fog top fall together.
        let lid = if t_h > sr { top - burn * (t_h - sr) } else { top };
        if lid <= p0[i].min(z) { continue; }
        pool_top[i] = lid as f32;
        let mut dd = 0.0;
        for (c, wt) in wts {
            if wt > 0.0 { dd += wt * cells[idx(c)].dd_at(top); }
        }
        let mut f = if dd <= 0.0 { top } else if dt > dd { p0[i] + d * (1.0 - dd / dt) } else { f64::NEG_INFINITY };
        if t_h > sr && f.is_finite() { f -= burn * (t_h - sr); }
        if f.is_finite() && z < f {
            mask[i] = 1;
            fog_top[i] = f as f32;
            if core {
                nf += 1;
                fmin = fmin.min(f);
                fmax = fmax.max(f);
            }
        } else if z < lid {
            mask[i] = 2;
            if core { np += 1; }
        }
    }
    ColdPoolProduct {
        n,
        fog_top_m: fog_top,
        pool_top_m: pool_top,
        mask,
        cells,
        fog_core_fraction: nf as f64 / nc.max(1) as f64,
        pool_core_fraction: (nf + np) as f64 / nc.max(1) as f64,
        fog_top_core_min_m: fmin,
        fog_top_core_max_m: fmax,
    }
}
