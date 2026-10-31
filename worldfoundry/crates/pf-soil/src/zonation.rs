//! Per-post zonation rules: pure functions of the post's inputs (no neighbourhood, no draws), so every tier and
//! every window gives the same answer for the same inputs.

use crate::params::{SoilParams, CLOSED_FW, OPEN_FW, WOODLAND_FW};

/// veg_class codes (u8). 255 = no data (outside the stencil-valid area).
pub const C_WATER: u8 = 0;
pub const C_CLOSED: u8 = 1;
pub const C_WOODLAND: u8 = 2;
pub const C_OPEN: u8 = 3;
pub const C_ECOTONE: u8 = 4;
pub const C_ALPINE_MEADOW: u8 = 5;
pub const C_FELLFIELD: u8 = 6;
pub const C_NIVAL: u8 = 7;
pub const C_BARE_ROCK: u8 = 8;
pub const C_SCREE: u8 = 9;
pub const C_FRESH_VOLCANIC: u8 = 10;
pub const C_KILL_ZONE: u8 = 11;
pub const C_RIPARIAN: u8 = 12;
pub const C_MARSH: u8 = 13;
pub const C_SALTMARSH: u8 = 14;
pub const C_INTERTIDAL: u8 = 15;
pub const C_BEACH: u8 = 16;
pub const C_DRY_SPARSE: u8 = 17;
pub const C_NODATA: u8 = 255;

pub const CLASS_NAMES: [(u8, &str, &str); 19] = [
    (C_WATER, "water", "sea, lake or carved channel (2 m) / large river"),
    (C_CLOSED, "closed_lowland", "below the tree line, humid (f_W >= 0.8): closed woody vegetation potential (forest / tall scrub)"),
    (C_WOODLAND, "woodland_shrub", "below the tree line, 0.4 <= f_W < 0.8: open woodland / shrubland"),
    (C_OPEN, "grassland_steppe", "below the tree line, 0.2 <= f_W < 0.4: grassland / steppe"),
    (C_ECOTONE, "treeline_ecotone", "within ecotone_m above the tree line: krummholz, tree islands, tall alpine shrub"),
    (C_ALPINE_MEADOW, "alpine_meadow", "above the ecotone, summer warmth index >= 12: closed alpine grassland / dwarf shrub"),
    (C_FELLFIELD, "fellfield_sparse", "above the ecotone, SWI < 12: sparse fell-field / cushion plants, cover from the CAVM curve"),
    (C_NIVAL, "nival", "above the Ohmura snow line (or the soft ramp's upper half): perennial snow / ice or frost-shattered barren ground"),
    (C_BARE_ROCK, "bare_rock", "cliff regime or slope factor < 0.3: rock outcrop"),
    (C_SCREE, "scree", "talus regime with cover < 0.5, or fresh rockfall debris"),
    (C_FRESH_VOLCANIC, "fresh_volcanic", "volcanic unit too young to be colonised (f_age < 0.5)"),
    (C_KILL_ZONE, "vent_tephra_kill", "active vent ring or recent tephra burial (f_kill < 0.5)"),
    (C_RIPARIAN, "riparian", "channel post (32 m), bank / corridor or band beside a carved channel (2 m), or TWI >= riparian_twi"),
    (C_MARSH, "marsh_mire", "TWI >= marsh_twi on flats, or lake fringe up to lake_marsh_m above the stage"),
    (C_SALTMARSH, "salt_marsh", "tide_m < height above sea <= tide_m + saltmarsh_band_m, flat, near the shore"),
    (C_INTERTIDAL, "intertidal", "land up to tide_m above sea near the shore: foreshore, unvegetated"),
    (C_BEACH, "beach", "T1b beach over-code"),
    (C_DRY_SPARSE, "dry_sparse", "below the tree line, f_W < 0.2: semi-desert / desert"),
    (C_NODATA, "no_data", "outside the stencil-valid area (2 m ring, 32 m lattice edge)"),
];

/// limit_code: which factor limits cover most (when the smallest factor is < 0.95).
pub const LIMIT_NAMES: [&str; 6] = ["none", "temperature (treeline / alpine / snow line)", "moisture (aridity)", "substrate (cliff / talus / slope / debris)", "volcanic surface age", "vent / tephra kill"];

/// Inputs of one post.
#[derive(Clone, Debug, Default)]
pub struct PostIn {
    pub h: f64,
    pub sea: f64,
    pub water: bool,
    /// Climate at the post (bilinear between T0 cell centres): sea-level-reduced annual mean (C), seasonal half range (K),
    /// precipitation (mm/yr), aridity index P / PET.
    pub t_sl: f64,
    pub half: f64,
    pub p: f64,
    pub ai: f64,
    pub lapse: f64,
    /// Warm-half-year direct-beam ratio slope / flat.
    pub r_sun: f64,
    /// Stability slope (deg) and regime (0 soil, 1 talus, 2 cliff, 3 tread).
    pub slope_deg: f64,
    pub regime: u8,
    /// Regime cover factor override (2 m tier: the factor box-averaged over substrate_smooth_2m_posts); None = from `regime`.
    pub reg_factor: Option<f64>,
    pub twi: f64,
    /// Volcanic colonisation factor (1 off the edifice) and the vent / tephra kill factor.
    pub f_age: f64,
    pub f_kill: f64,
    /// Soil-age factor 1 - exp(-t / tau_soil) (1 on non-volcanic ground).
    pub soil_age_f: f64,
    /// pf-talus debris cover (0..1; 0 where absent).
    pub debris: f64,
    pub beach: bool,
    /// Riparian forced (channel post at 32 m; bank / near channel at 2 m).
    pub riparian: bool,
    /// Fresh marsh forced (lake fringe).
    pub lake_fringe: bool,
    /// Gentle-ground slope for the wet classes (deg).
    pub flat_slope_deg: f64,
    /// Distance to the sea shoreline (m); None when unknown.
    pub shore_m: Option<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PostOut {
    pub cover: f64,
    pub class: u8,
    pub limit: u8,
    pub tl_margin_m: f64,
    pub niv_margin_m: f64,
    pub soil_m: f64,
    pub t_gs: f64,
}

pub fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Growing season of the sinusoidal cycle T(t) = tm + a cos(2 pi t): (season fraction of the year, mean T over the
/// season). No season: (0, -inf).
pub fn season(tm: f64, a: f64, t0: f64) -> (f64, f64) {
    if a <= 1e-12 { return if tm > t0 { (1.0, tm) } else { (0.0, f64::NEG_INFINITY) }; }
    let c = (t0 - tm) / a;
    if c <= -1.0 { return (1.0, tm); }
    if c >= 1.0 { return (0.0, f64::NEG_INFINITY); }
    let th = c.acos();
    (th / std::f64::consts::PI, tm + a * th.sin() / th)
}

/// Paulsen-Koerner tree criterion at annual mean tm and half range a.
pub fn tree_ok(p: &SoilParams, tm: f64, a: f64) -> bool {
    let (f, tgs) = season(tm, a, p.t_season_c);
    f * 365.25 >= p.lgs_min_days && tgs >= p.t_gs_min_c
}

/// The lowest annual mean temperature at which trees pass, for half range a (bisection; monotone in tm).
pub fn tm_treeline(p: &SoilParams, a: f64) -> f64 {
    let (mut lo, mut hi) = (-80.0f64, 80.0f64);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        if tree_ok(p, mid, a) { hi = mid; } else { lo = mid; }
    }
    hi
}

/// Summer warmth index: sum of the positive monthly means of the twelve samples tm + a cos(2 pi m / 12).
pub fn swi(tm: f64, a: f64) -> f64 {
    let mut s = 0.0;
    for m in 0..12 {
        let v = tm + a * (std::f64::consts::TAU * m as f64 / 12.0).cos();
        if v > 0.0 { s += v; }
    }
    s
}

/// CAVM cover curve (piecewise linear).
pub fn cavm_cover(p: &SoilParams, s: f64) -> f64 {
    let pts = [(0.0, 0.0), (p.cav_swi_a, p.cav_cov_a), (p.cav_swi_b, p.cav_cov_b), (p.cav_swi_c, p.cav_cov_c), (p.cav_swi_d, p.cav_cov_d), (p.cav_swi_e, 1.0)];
    if s <= 0.0 { return 0.0; }
    for k in 1..pts.len() {
        let (x0, y0) = pts[k - 1];
        let (x1, y1) = pts[k];
        if s <= x1 { return y0 + (y1 - y0) * (s - x0) / (x1 - x0); }
    }
    1.0
}

/// Ohmura et al. 1992: JJA temperature at the equilibrium line for annual precipitation pmm.
pub fn ohmura_t_ela(p: &SoilParams, pmm: f64) -> f64 {
    let pmm = pmm.max(0.0);
    let disc = p.ohm_b * p.ohm_b - 4.0 * p.ohm_c * (p.ohm_a - pmm);
    (-p.ohm_b + disc.max(0.0).sqrt()) / (2.0 * p.ohm_c)
}

/// Warm-half-year direct-beam ratio (slope / flat) at latitude lat_rad for a surface of gradient (zx east, zy north).
pub fn insolation_ratio(p: &SoilParams, lat_rad: f64, zx: f64, zy: f64) -> f64 {
    let nd = p.insol_decl_n as usize;
    let nh = p.insol_hours_n as usize;
    let s = if lat_rad < 0.0 { -1.0 } else { 1.0 };
    let tilt = p.tilt_deg.to_radians();
    let nn = (zx * zx + zy * zy + 1.0).sqrt();
    let (nx, ny, nz) = (-zx / nn, -zy / nn, 1.0 / nn);
    let (sp, cp) = lat_rad.sin_cos();
    let (mut sl, mut fl) = (0.0, 0.0);
    for k in 0..nd {
        let u = -std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (k as f64 + 0.5) / nd as f64;
        let d = s * tilt * u.cos();
        let (sd, cd) = d.sin_cos();
        for j in 0..nh {
            let hh = -std::f64::consts::PI + std::f64::consts::TAU * (j as f64 + 0.5) / nh as f64;
            let (sh, ch) = hh.sin_cos();
            let up = sp * sd + cp * cd * ch;
            if up <= 0.0 { continue; }
            let e = -cd * sh;
            let n = cp * sd - sp * cd * ch;
            fl += up;
            let ci = e * nx + n * ny + up * nz;
            if ci > 0.0 { sl += ci; }
        }
    }
    if fl > 0.0 { sl / fl } else { 1.0 }
}

/// Cover factor of a cliff-and-talus regime code (0 soil, 1 talus, 2 cliff, 3 tread).
pub fn regime_factor(p: &SoilParams, regime: u8) -> f64 { match regime { 1 => p.f_talus, 2 => p.f_cliff, 3 => p.f_tread, _ => 1.0 } }

/// Volcanic colonisation factor for a surface of age t (yr) and base time scale tau (yr) at growing-season temperature
/// t_gs and aridity ai.
pub fn age_factor(p: &SoilParams, t: f64, tau: f64, t_gs: f64, ai: f64) -> f64 {
    let tg = if t_gs.is_finite() { t_gs.max(0.0) } else { 0.0 };
    let tau_c = tau * p.q10.powf((p.t_ref_c - tg) / 10.0) * (p.ai_hi / ai.max(1e-6)).max(1.0);
    1.0 - (-(t.max(0.0) / tau_c).powf(p.age_weibull_k)).exp()
}

/// Growing-season temperature at the post (for the age factor of a neighbourhood field): T_gs of the lapsed mean.
pub fn post_t_gs(p: &SoilParams, x: &PostIn) -> f64 {
    let tm = x.t_sl - x.lapse * (x.h - x.sea).max(0.0) + (p.k_sun_k * (x.r_sun - 1.0)).clamp(-p.dt_sun_max_k, p.dt_sun_max_k);
    season(tm, x.half, p.t_season_c).1
}

/// The zonation of one post.
pub fn zone(p: &SoilParams, x: &PostIn) -> PostOut {
    let dt_sun = (p.k_sun_k * (x.r_sun - 1.0)).clamp(-p.dt_sun_max_k, p.dt_sun_max_k);
    let tm = x.t_sl - x.lapse * (x.h - x.sea).max(0.0) + dt_sun;
    let a = x.half.max(0.0);
    let tl_margin = (tm - tm_treeline(p, a)) / x.lapse;
    let t_jja = tm + p.jja_factor * a;
    let niv_margin = (t_jja - ohmura_t_ela(p, x.p)) / x.lapse;
    let t_gs = season(tm, a, p.t_season_c).1;
    if x.water {
        return PostOut { cover: 0.0, class: C_WATER, limit: 0, tl_margin_m: tl_margin, niv_margin_m: niv_margin, soil_m: 0.0, t_gs };
    }
    let sw = swi(tm, a);
    let f_t = if tl_margin >= 0.0 { 1.0 } else { cavm_cover(p, sw) };
    let f_niv = smoothstep(-p.nival_ramp_m, p.nival_ramp_m, niv_margin);
    let ft = f_t * f_niv;
    let ai_loc = x.ai * (1.0 + p.twi_wet_gain * (x.twi - p.twi_wet_ref).max(0.0)) / (1.0 + p.pet_sun_gain * (x.r_sun - 1.0)).max(0.1);
    let f_w = ((ai_loc - p.ai_lo) / (p.ai_hi - p.ai_lo)).clamp(0.0, 1.0);
    let f_reg = x.reg_factor.unwrap_or_else(|| regime_factor(p, x.regime));
    let f_slope = 1.0 - (1.0 - p.slope_floor) * smoothstep(p.slope_soil_lo_deg, p.slope_soil_hi_deg, x.slope_deg);
    let f_s = f_reg.min(f_slope) * (1.0 - p.debris_gain * x.debris.clamp(0.0, 1.0));
    let f_age = x.f_age.clamp(0.0, 1.0);
    let f_kill = x.f_kill.clamp(0.0, 1.0);
    let base = ft * f_w * f_s * f_age * f_kill;
    // wet classes
    let dz = x.h - x.sea;
    let near = |lim: f64| x.shore_m.map(|s| s <= lim).unwrap_or(false);
    let intertidal = dz <= p.tide_m && near(p.intertidal_shore_m);
    let saltmarsh = !intertidal && dz > p.tide_m && dz <= p.tide_m + p.saltmarsh_band_m && x.flat_slope_deg <= p.saltmarsh_slope_deg && near(p.saltmarsh_shore_m);
    let marsh = x.lake_fringe || (x.twi >= p.marsh_twi && x.flat_slope_deg <= p.marsh_slope_deg);
    let riparian = x.riparian || x.twi >= p.riparian_twi;
    let wet = saltmarsh || marsh || riparian;
    let cover = if intertidal || x.beach { 0.0 } else if wet { base.max(p.wet_cover * ft * f_age * f_kill) } else { base };
    // limit code
    // (wet classes are not moisture- or slope-limited: those factors do not enter their cover)
    let facs = if wet { [ft, 1.0, 1.0, f_age, f_kill] } else { [ft, f_w, f_s, f_age, f_kill] };
    let mut limit = 0u8;
    let mut mn = 0.95;
    for (k, &f) in facs.iter().enumerate() { if f < mn { mn = f; limit = k as u8 + 1; } }
    // class (priority order)
    let class = if x.beach { C_BEACH }
        else if intertidal { C_INTERTIDAL }
        else if niv_margin < 0.0 { C_NIVAL }
        else if f_kill < 0.5 { C_KILL_ZONE }
        else if f_age < 0.5 { C_FRESH_VOLCANIC }
        else if !wet && (x.regime == 2 || f_slope < 0.3) { C_BARE_ROCK }
        else if !wet && ((x.regime == 1 && cover < 0.5) || x.debris >= 0.5) { C_SCREE }
        else if saltmarsh { C_SALTMARSH }
        else if marsh { C_MARSH }
        else if riparian { C_RIPARIAN }
        else if tl_margin >= 0.0 {
            if f_w >= CLOSED_FW { C_CLOSED } else if f_w >= WOODLAND_FW { C_WOODLAND } else if f_w >= OPEN_FW { C_OPEN } else { C_DRY_SPARSE }
        } else if tl_margin > -p.ecotone_m { C_ECOTONE }
        else if sw >= p.swi_meadow { C_ALPINE_MEADOW }
        else { C_FELLFIELD };
    // soil depth proxy
    let reg_soil = match x.regime { 1 => 0.2, 2 => 0.0, 3 => 0.6, _ => 1.0 };
    let soil = if intertidal || x.beach { 0.0 } else {
        p.soil_dmax_m * x.soil_age_f.clamp(0.0, 1.0) * (1.0 - smoothstep(p.soil_slope_lo_deg, p.soil_slope_hi_deg, x.slope_deg)) * reg_soil
            * (1.0 - x.debris.clamp(0.0, 1.0)) * (0.25 + 0.75 * ft)
    };
    PostOut { cover: cover.clamp(0.0, 1.0), class, limit, tl_margin_m: tl_margin, niv_margin_m: niv_margin, soil_m: soil, t_gs }
}
