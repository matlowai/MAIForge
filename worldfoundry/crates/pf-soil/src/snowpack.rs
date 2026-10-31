//! SNOWPACK-1 (lane SNOWPACK-1, 2026-10-22; model specification "seasonal snowpack"): a seasonal snowpack per post
//! from the same climate pf-soil reads, so the ground can be drawn on any day of the year.
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon; a presentation / ecology proxy (the kernel never reads
//! it). Opt-in: a separate product (`vz.snowpack/v1`, bin `pf-snow`); pf-soil's own outputs do not change.
//!
//! MODEL (per post, daily step, constants and source classes in `SnowParams`):
//! 1. Temperature: the pf-soil post temperature (seed-climate-v1 reduced to sea level, bilinear, lapsed to the post:
//!    `zonation::zone`'s tm without the aspect term, because aspect enters through radiation in 3) on the climate's one
//!    harmonic, phased as pf-weather (`season_peak_day_north`, southern hemisphere half a year later):
//!    T_d = T_mean + A cos(2 pi (d - d_peak) / 365), d the 0-based day of the 365-day calendar.
//! 2. Day-to-day spread: daily mean temperature ~ N(T_d, sigma_t^2). Expected positive degree-days per day
//!    PDD = sigma phi(T/sigma) + T Phi(T/sigma) (Calov and Greve 2005), expected snow share of precipitation by
//!    pf-weather's own precipitation-type thresholds (snow at T <= t_snow_c, mixed half and half below t_rain_c, rain above).
//! 3. Melt: Hock (1999) Eq. 2 (hourly M = (MF / 24 + a_snow I) T) summed over the day, M = ms (MF + 24 a_snow I_d) PDD, with I_d the daily-mean
//!    clear-sky potential direct solar irradiance on the post's slope (Hock 1999's formula: I0 psi_a^(p / (p0 cos Z))
//!    cos(theta), no horizon shading, circular orbit as the climate artifact declares).
//! 4. Accumulation: annual precipitation spread evenly over the year (DECLARED: the T0 / seed-climate-v1 tier has no
//!    precipitation seasonality; no pattern is invented) times the snow share.
//! 5. Snowpack: SWE += snowfall; SWE -= min(SWE, M); capped at swe_cap_mm (firn / ice conversion not modelled).
//!    Spin-up from SWE 0 on the warmest day for spin_years years; outputs are the last year's.
//! 6. Snow-cover fraction: fsc = tanh(SWE / swe_cover_mm) (a depletion-curve form after Yang et al. 1997 / Niu and Yang
//!    2007 with a SWE scale instead of depth; the scale is a declared presentation-level assumption).
//! 7. Annual climatic balance b = sum(snowfall) - sum(M) of a snow surface (mm w.e./yr): b = 0 is the model's
//!    equilibrium line; `perennial` = SWE stays > 0 all of the last spin-up year.
//! 8. Ohmura anchor (default on): the melt multiplier ms of each post's climate is the one that puts b = 0 on flat ground
//!    exactly at that climate's Ohmura et al. 1992 line (pf-soil nival_margin = 0, GLACIER-1's firn line), so the
//!    late-summer snow line agrees with the perennial line by construction and aspect / season come from the model.
//!    Lane SNOWPACK-1 gate 1: the free model (ms = 1, Hock's Storglaciaren values) puts the line +241 m above Ohmura at
//!    site 5 and +643 m at site 1; the anchored ms is 0.55 (site 5) and 0.24 (site 1, below the literature DDF range:
//!    flat precipitation makes half the precipitation rain at a warm maritime equilibrium line).
//!
//! Determinism: no draws; every post is a pure function of its height, its gradient stencil and its chart position;
//! f64 in a fixed order; transcendental functions are the platform's (as pf-weather); outputs quantized.

use serde_json::{json, Value};
use std::f64::consts::{PI, TAU};

pub const SCHEMA: &str = "vz.snowpack/v1";
pub const GENERATOR: &str = "pf-snow";
pub const GENERATOR_VERSION: &str = "0.1.0";
pub const LABEL: &str = "pre-signoff; fidelity tier Preview; derived, never canon; seasonal snowpack proxy (degree-day + potential direct radiation) from the pf-soil climate (lane SNOWPACK-1)";
pub const DAYS: usize = 365;

macro_rules! snow_params {
    ($( $name:ident : $val:expr, $unit:expr, $class:expr, $basis:expr; )*) => {
        /// Snowpack constants (all f64).
        #[derive(Clone, Debug, PartialEq)]
        pub struct SnowParams { $(pub $name: f64,)* }
        impl SnowParams {
            pub fn v1() -> Self { SnowParams { $($name: $val,)* } }
            pub fn set(&mut self, k: &str, v: f64) -> Result<(), String> {
                if !v.is_finite() { return Err(format!("param {k}: non-finite value")); }
                match k { $(stringify!($name) => { self.$name = v; Ok(()) })* _ => Err(format!("unknown param {k}")) }
            }
            pub fn table(&self) -> Value {
                let mut v: Vec<Value> = Vec::new();
                $( v.push(json!({"name": stringify!($name), "value": self.$name, "unit": $unit, "class": $class, "basis": $basis})); )*
                Value::Array(v)
            }
        }
    };
}

snow_params! {
    // ---- melt (Hock 1999 enhanced temperature index) ----
    mf_mm_d_k: 1.8, "mm w.e. d-1 K-1", "R", "Hock 1999 (J. Glaciol. 45(149) 101-111, Table 1 model 2) melt factor MF, Storglaciaren fit (CHECKED by lane SNOWPACK-1, records/literature.txt); Eq. 2 is M = (MF / 24 + a I) T per hour";
    a_snow: 0.0006, "mm w.e. h-1 (W m-2)-1 K-1", "R", "Hock 1999 Table 1 model 2 a_snow (CHECKED); the daily form 24 a_snow x the daily-mean irradiance is the day sum of Eq. 2 with the hour-scale T-I correlation ignored [A] (absorbed by the melt multiplier)";
    melt_scale: 1.0, "1", "A", "one multiplier on the melt term; v1 value calibrated so the model equilibrium line matches the Ohmura et al. 1992 line on flat ground over a climate grid (lane SNOWPACK-1 gate 1); 1.0 = Hock's Storglaciaren values";
    i0_w_m2: 1368.0, "W m-2", "R", "solar constant used by Hock 1999 (CHECKED)";
    psi_a: 0.75, "1", "R", "clear-sky atmospheric transmissivity, Hock 1999 value (CHECKED; literature range 0.6-0.9 quoted there)";
    p_scale_h_m: 8434.0, "m", "R", "barometric scale height for p / p0 = exp(-z / H) in the air-mass term (isothermal 288 K atmosphere, recalled)";
    tilt_deg: 23.44, "deg", "V", "axial tilt; the climate artifact declares 23.44 (axial_tilt_deg); declination by Cooper 1969 as pf-weather sun.rs";
    // ---- temperature ----
    season_peak_day_north: 200.0, "day", "V", "pf-weather constants.rs season_peak_day_north (0-based day of the warmest day, northern hemisphere; southern hemisphere minus 182.5): project consistency with the game's weather";
    sigma_t_k: 3.5, "K", "R", "standard deviation of the daily mean temperature about the seasonal cycle; Seguinot 2013 (J. Glaciol. 59(218)) ERA-Interim daily means about monthly means: 0.32-12.63 K, global mean 2.02 K, higher over continents and in winter (CHECKED); ice-sheet PDD convention 4.5-5.5 K includes the diurnal cycle (Fausto et al. 2011, CHECKED); 3.5 K assumed for mid / high-latitude land [A]; with the Ohmura anchor the line does not depend on it";
    t_snow_c: 0.5, "degC", "V", "pf-weather snow_max_t_c (precipitation is snow at or below this air temperature): project consistency";
    t_rain_c: 2.5, "degC", "V", "pf-weather rain_min_t_c (all rain above; mixed half and half between)";
    // ---- snowpack ----
    swe_cover_mm: 15.0, "mm w.e.", "A", "fsc = tanh(SWE / this): the tanh form of Niu and Yang 2007 (CHECKED via CTSM / Noah-MP code) with a SWE scale; ~5 cm of settled snow hides low ground cover; Noah SNUP 0.02-0.04 (grass / tundra) for 100 % cover (unit NOT CHECKED); value assumed";
    swe_cap_mm: 10000.0, "mm w.e.", "A", "SWE ceiling (perennial snow turns to firn / ice; not modelled here, GLACIER-1 owns ice)";
    spin_years: 3.0, "yr", "A", "spin-up years from SWE 0 on the warmest day; outputs from the last year";
    // ---- Ohmura anchor (consistency with pf-soil nival_margin and GLACIER-1's firn line) ----
    anchor_ohmura: 1.0, "flag", "A", "1: per post, melt_scale is replaced by the value that puts the flat-ground balance b = 0 exactly at the post's own Ohmura line (b is linear in the melt multiplier: ms = sum S / sum M of a flat surface at that height); aspect and season then come from the model. 0: the free model with melt_scale (gate 1 reports both)";
    ohm_a: 645.0, "mm/yr", "V", "copy of pf-soil params ohm_a (Ohmura et al. 1992; a unit test asserts equality)";
    ohm_b: 296.0, "mm/yr/K", "V", "copy of pf-soil params ohm_b";
    ohm_c: 9.0, "mm/yr/K2", "V", "copy of pf-soil params ohm_c";
    jja_factor: 0.900316316157106, "1", "V", "copy of pf-soil params jja_factor (warmest-quarter mean of the cosine cycle)";
    ms_min: 0.05, "1", "A", "anchor clamp: lowest per-post melt multiplier (counted when hit)";
    ms_max: 5.0, "1", "A", "anchor clamp: highest per-post melt multiplier (counted when hit)";
    // ---- quadrature (exact counts, not physics) ----
    decl_nodes: 25.0, "1", "V", "declination nodes over [-tilt, +tilt] of the per-post daily irradiance table (linear in between)";
    hours_n: 48.0, "1", "V", "hour-angle samples per day of the irradiance quadrature";
    // ---- 2 m tier ----
    insol_baseline_2m_posts: 8.0, "post", "V", "2 m tier: slope for the irradiance from central differences over +-8 posts (16 m half baseline), the same canopy-scale aspect baseline as pf-soil's insol_baseline_2m_posts";
}

pub const RULES: [(&str, &str); 7] = [
    ("temperature", "t_post_v1: T_mean = T_sl(bilinear seed-climate-v1, pf-soil climate.rs) - lapse max(h - sea, 0) (no aspect term); T_d = T_mean + A cos(2 pi (d - d_peak) / 365), d_peak = season_peak_day_north (north) or minus 182.5 (south); daily mean ~ N(T_d, sigma_t^2)"),
    ("pdd", "pdd_normal_v1: PDD_d = s phi(T_d / s) + T_d Phi(T_d / s) (Calov and Greve 2005); Phi via erfc (Numerical Recipes erfcc, fractional error < 1.2e-7)"),
    ("snowfall", "snowfall_v1: S_d = (P / 365) [Phi((t_snow - T_d) / s) + 0.5 (Phi((t_rain - T_d) / s) - Phi((t_snow - T_d) / s))]; P the bilinear annual precipitation, flat in season (declared)"),
    ("melt", "hock_eti_daily_v1: M_d = melt_scale (MF + 24 a_snow I_d) PDD_d; I_d = mean over hours_n hour angles of I0 psi_a^(exp(-z / H) / sin e) max(0, s . n) with the sun up (e > 0), z = max(h - sea, 0), declination delta_d = tilt sin(2 pi (285 + d) / 365) (pf-weather sun.rs); I tabulated per post at decl_nodes declinations, linear in delta"),
    ("pack", "swe_v1: start on the warmest day with SWE 0; each day SWE = min(SWE + S_d, cap); SWE -= min(SWE, M_d); spin_years years; outputs SWE on the requested days of the last year, fsc = tanh(SWE / swe_cover_mm), snow_days = days of the last year with fsc >= 0.5, perennial = min SWE of the last year > 0, balance = sum S_d - sum M_d"),
    ("anchor", "ohmura_anchor_v1 (anchor_ohmura = 1): z* = (T_sl + jja_factor A - T_ela(P)) / lapse with P = ohm_a + ohm_b T_ela + ohm_c T_ela^2 (pf-soil nival_margin = 0 on flat ground, no sea clamp); ms = clamp(sum_d S_d(z*) / sum_d M_d(z*; ms = 1, flat), ms_min, ms_max); the post's melt is ms (MF + 24 a_snow I_d) PDD_d"),
    ("slope", "gradient_v1: 32 m central differences over +-1 post; 2 m central differences over +-insol_baseline_2m_posts posts (pf-soil's canopy-scale aspect baseline)"),
];

pub fn rules_json() -> Value {
    Value::Object(RULES.iter().map(|(k, v)| (k.to_string(), Value::String(v.to_string()))).collect())
}

/// erfc, Numerical Recipes `erfcc` (fractional error < 1.2e-7 everywhere).
pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t * (-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806
        + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))).exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}

/// Standard normal CDF and PDF.
pub fn ncdf(z: f64) -> f64 { 0.5 * erfc(-z / std::f64::consts::SQRT_2) }
pub fn npdf(z: f64) -> f64 { (-0.5 * z * z).exp() / (TAU).sqrt() }

/// Expected positive part of N(t, s^2) (degree-days per day).
pub fn pdd(t: f64, s: f64) -> f64 {
    if s <= 1e-9 { return t.max(0.0); }
    let z = t / s;
    s * npdf(z) + t * ncdf(z)
}

/// Expected snow share of precipitation at daily mean t (pf-weather's type thresholds under N(t, s^2)).
pub fn snow_share(p: &SnowParams, t: f64) -> f64 {
    let s = p.sigma_t_k.max(1e-9);
    let a = ncdf((p.t_snow_c - t) / s);
    let b = ncdf((p.t_rain_c - t) / s);
    a + 0.5 * (b - a)
}

/// Declination (rad) on 0-based day d (pf-weather sun.rs: tilt sin(2 pi (284 + d + 1) / 365)).
pub fn declination(p: &SnowParams, d: usize) -> f64 {
    p.tilt_deg.to_radians() * (TAU * (284.0 + d as f64 + 1.0) / 365.0).sin()
}

/// 0-based warmest day for a latitude (deg).
pub fn peak_day(p: &SnowParams, lat_rad: f64) -> f64 {
    if lat_rad >= 0.0 { p.season_peak_day_north } else { p.season_peak_day_north - 182.5 }
}

/// Daily-mean clear-sky potential direct irradiance (W m-2) on a surface of gradient (zx east, zy north) at latitude
/// lat_rad, declination decl, pressure ratio pr = p / p0.
pub fn daily_direct(p: &SnowParams, lat_rad: f64, zx: f64, zy: f64, decl: f64, pr: f64) -> f64 {
    let nh = p.hours_n.max(1.0) as usize;
    let nn = (zx * zx + zy * zy + 1.0).sqrt();
    let (nx, ny, nz) = (-zx / nn, -zy / nn, 1.0 / nn);
    let (sp, cp) = lat_rad.sin_cos();
    let (sd, cd) = decl.sin_cos();
    let lnpsi = p.psi_a.ln();
    let mut acc = 0.0;
    for j in 0..nh {
        let hh = -PI + TAU * (j as f64 + 0.5) / nh as f64;
        let (sh, ch) = hh.sin_cos();
        let up = sp * sd + cp * cd * ch;
        if up <= 0.0 { continue; }
        let e = -cd * sh;
        let n = cp * sd - sp * cd * ch;
        let ci = e * nx + n * ny + up * nz;
        if ci <= 0.0 { continue; }
        acc += p.i0_w_m2 * (lnpsi * pr / up).exp() * ci;
    }
    acc / nh as f64
}

/// Per-day forcing shared by every post of one hemisphere: cos of the seasonal phase and (node, weight) of the
/// declination in the irradiance table.
pub struct Calendar { pub cos_n: Vec<f64>, pub cos_s: Vec<f64>, pub node: Vec<(usize, f64)>, pub decl_nodes: Vec<f64> }

impl Calendar {
    pub fn new(p: &SnowParams) -> Calendar {
        let nd = (p.decl_nodes.max(2.0)) as usize;
        let tilt = p.tilt_deg.to_radians();
        let decl_nodes: Vec<f64> = (0..nd).map(|k| -tilt + 2.0 * tilt * k as f64 / (nd - 1) as f64).collect();
        let mut node = Vec::with_capacity(DAYS);
        for d in 0..DAYS {
            let x = if tilt > 0.0 { ((declination(p, d) + tilt) / (2.0 * tilt) * (nd - 1) as f64).clamp(0.0, (nd - 1) as f64) } else { 0.0 };
            let k = (x.floor() as usize).min(nd - 2);
            node.push((k, x - k as f64));
        }
        let cos_at = |peak: f64| (0..DAYS).map(|d| (TAU * (d as f64 - peak) / 365.0).cos()).collect::<Vec<f64>>();
        Calendar { cos_n: cos_at(p.season_peak_day_north), cos_s: cos_at(p.season_peak_day_north - 182.5), node, decl_nodes }
    }
}

/// Inputs of one post.
#[derive(Clone, Copy, Debug, Default)]
pub struct SnowIn {
    pub lat_rad: f64,
    /// Sea-level-reduced annual mean (C), half range (K), annual precipitation (mm), lapse (K/m).
    pub t_sl: f64,
    pub half: f64,
    pub p_mm: f64,
    pub lapse: f64,
    pub h: f64,
    pub sea: f64,
    pub zx: f64,
    pub zy: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SnowOut {
    /// SWE (mm) on each requested day, last spin-up year.
    pub swe: Vec<f64>,
    pub snow_days: u16,
    pub perennial: bool,
    pub balance_mm: f64,
    pub min_swe: f64,
    pub max_swe: f64,
    /// Annual mean daily-mean irradiance on the slope / on flat ground (diagnostic).
    pub irr_ratio: f64,
    /// The melt multiplier used (melt_scale, or the Ohmura-anchored value) and whether the anchor clamp was hit.
    pub ms: f64,
    pub ms_clamped: bool,
}

pub fn fsc(p: &SnowParams, swe: f64) -> f64 { if swe <= 0.0 { 0.0 } else { (swe / p.swe_cover_mm.max(1e-9)).tanh() } }

/// The melt multiplier of a post's climate: melt_scale (anchor off), or the value that puts the flat-ground balance
/// b = 0 at the post's Ohmura line z* (b is linear in the multiplier). Returns (ms, clamp hit).
pub fn melt_multiplier(p: &SnowParams, cal: &Calendar, lat_rad: f64, t_sl: f64, half: f64, p_mm: f64, lapse: f64) -> (f64, bool) {
    if p.anchor_ohmura < 0.5 { return (p.melt_scale, false); }
    let a = half.max(0.0);
    let pmm = p_mm.max(0.0);
    let disc = p.ohm_b * p.ohm_b - 4.0 * p.ohm_c * (p.ohm_a - pmm);
    let t_ela = (-p.ohm_b + disc.max(0.0).sqrt()) / (2.0 * p.ohm_c);
    let zs = (t_sl + p.jja_factor * a - t_ela) / lapse;
    let tz = t_sl - lapse * zs;
    let pr = (-zs.max(0.0) / p.p_scale_h_m).exp();
    let tf: Vec<f64> = cal.decl_nodes.iter().map(|&dl| daily_direct(p, lat_rad, 0.0, 0.0, dl, pr)).collect();
    let cs = if lat_rad >= 0.0 { &cal.cos_n } else { &cal.cos_s };
    let pday = pmm / DAYS as f64;
    let (mut ss, mut mm) = (0.0, 0.0);
    for d in 0..DAYS {
        let t = tz + a * cs[d];
        let (k, w) = cal.node[d];
        ss += pday * snow_share(p, t);
        mm += (p.mf_mm_d_k + 24.0 * p.a_snow * (tf[k] * (1.0 - w) + tf[k + 1] * w)) * pdd(t, p.sigma_t_k);
    }
    let raw = if mm > 0.0 { ss / mm } else { p.ms_max };
    let c = raw.clamp(p.ms_min, p.ms_max);
    (c, c != raw)
}

/// The snowpack of one post; `days` are 0-based days of the calendar.
pub fn post(p: &SnowParams, cal: &Calendar, x: &SnowIn, days: &[usize]) -> SnowOut {
    let z = (x.h - x.sea).max(0.0);
    let tm = x.t_sl - x.lapse * z;
    let a = x.half.max(0.0);
    let pr_at = |zz: f64| (-zz.max(0.0) / p.p_scale_h_m).exp();
    let table = |zx: f64, zy: f64, pr: f64| -> Vec<f64> { cal.decl_nodes.iter().map(|&dl| daily_direct(p, x.lat_rad, zx, zy, dl, pr)).collect() };
    let irr_day = |tab: &[f64], d: usize| { let (k, w) = cal.node[d]; tab[k] * (1.0 - w) + tab[k + 1] * w };
    let tab = table(x.zx, x.zy, pr_at(z));
    let tab_flat = if x.zx == 0.0 && x.zy == 0.0 { tab.clone() } else { table(0.0, 0.0, pr_at(z)) };
    let cs = if x.lat_rad >= 0.0 { &cal.cos_n } else { &cal.cos_s };
    let pday = x.p_mm.max(0.0) / DAYS as f64;
    let s = p.sigma_t_k;
    // melt multiplier: free, or anchored to the post's own Ohmura line on flat ground
    let (ms, ms_clamped) = melt_multiplier(p, cal, x.lat_rad, x.t_sl, a, x.p_mm, x.lapse);
    let mut snow = [0.0f64; DAYS];
    let mut melt = [0.0f64; DAYS];
    let (mut irr_sum, mut flat_sum) = (0.0, 0.0);
    for d in 0..DAYS {
        let t = tm + a * cs[d];
        let irr = irr_day(&tab, d);
        irr_sum += irr;
        flat_sum += irr_day(&tab_flat, d);
        snow[d] = pday * snow_share(p, t);
        melt[d] = ms * (p.mf_mm_d_k + 24.0 * p.a_snow * irr) * pdd(t, s);
    }
    let lat_peak = peak_day(p, x.lat_rad);
    let d0 = (lat_peak.rem_euclid(365.0)).round() as usize % DAYS;
    let years = p.spin_years.max(1.0) as usize;
    let mut swe = 0.0f64;
    let mut last = [0.0f64; DAYS];
    for y in 0..years {
        for i in 0..DAYS {
            let d = (d0 + i) % DAYS;
            swe = (swe + snow[d]).min(p.swe_cap_mm);
            swe -= swe.min(melt[d]);
            if y + 1 == years { last[d] = swe; }
        }
    }
    let balance: f64 = snow.iter().sum::<f64>() - melt.iter().sum::<f64>();
    let min_swe = last.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_swe = last.iter().cloned().fold(0.0, f64::max);
    let snow_days = last.iter().filter(|&&w| fsc(p, w) >= 0.5).count() as u16;
    SnowOut { swe: days.iter().map(|&d| last[d % DAYS]).collect(), snow_days, perennial: min_swe > 0.0, balance_mm: balance, min_swe, max_swe,
        irr_ratio: if flat_sum > 0.0 { irr_sum / flat_sum } else { 1.0 }, ms, ms_clamped }
}

/// Annual climatic balance of a snow surface at a free height z (no clamp at sea level), flat ground (for the gates).
pub fn balance_at(p: &SnowParams, cal: &Calendar, lat_rad: f64, t_sl: f64, half: f64, p_mm: f64, lapse: f64, z: f64, zx: f64, zy: f64) -> f64 {
    // the multiplier belongs to the climate, not to the shifted column below sea level
    let (ms, _) = melt_multiplier(p, cal, lat_rad, t_sl, half, p_mm, lapse);
    let mut q = p.clone();
    q.anchor_ohmura = 0.0;
    q.melt_scale = ms;
    let x = SnowIn { lat_rad, t_sl: t_sl - lapse * z.min(0.0), half, p_mm, lapse, h: z.max(0.0), sea: 0.0, zx, zy };
    post(&q, cal, &x, &[]).balance_mm
}

/// Height (m) where the balance crosses 0 (bisection on [lo, hi]; balance rises with height). None if no sign change.
pub fn ela(p: &SnowParams, cal: &Calendar, lat_rad: f64, t_sl: f64, half: f64, p_mm: f64, lapse: f64, zx: f64, zy: f64, lo: f64, hi: f64) -> Option<f64> {
    let f = |z: f64| balance_at(p, cal, lat_rad, t_sl, half, p_mm, lapse, z, zx, zy);
    let (mut a, mut b) = (lo, hi);
    let (fa, fb) = (f(a), f(b));
    if fa > 0.0 || fb < 0.0 { return None; }
    for _ in 0..60 {
        let m = 0.5 * (a + b);
        if f(m) < 0.0 { a = m; } else { b = m; }
        if b - a < 0.01 { break; }
    }
    Some(0.5 * (a + b))
}

/// Gradient (zx east, zy north) of a planar slope of `slope_deg` facing `aspect_deg` (downslope direction, clockwise
/// from north: 0 = north-facing).
pub fn facing(slope_deg: f64, aspect_deg: f64) -> (f64, f64) {
    let g = slope_deg.to_radians().tan();
    let a = aspect_deg.to_radians();
    // downslope direction (sin a, cos a) east / north: the height falls that way
    (-g * a.sin(), -g * a.cos())
}
