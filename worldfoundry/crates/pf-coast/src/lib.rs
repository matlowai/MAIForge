//! Coastal potential classification. See MODEL.md for validity limits and sources.
//! No terrain modification, sediment transport, tide solver or vegetation creation.
pub mod publish;
use serde::{Deserialize, Serialize};

pub const G: f64 = 9.81;
pub const CLASS_NAMES: [&str; 8] = [
    "none",
    "rocky_cliff",
    "rocky_platform",
    "gravel_beach",
    "sand_beach",
    "tidal_flat",
    "salt_marsh",
    "delta_front",
];
// Native LithologyV1 priors, NOT laboratory measurements of coastal abrasion.
pub const K: [f64; 12] = [0.5, 0.4, 0.8, 0.6, 0.5, 2.0, 0.8, 0.9, 3.0, 1.5, 3.0, 0.7];

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Config {
    pub tide_range_m: f64,
    pub fetch_cap_m: f64,
    pub supply_reach_m: f64,
    pub reference_area_m2: f64,
}
impl Config {
    pub fn new(tide_range_m: f64) -> Result<Self, String> {
        let c = Self {
            tide_range_m,
            fetch_cap_m: 128_000.0,
            supply_reach_m: 512.0,
            reference_area_m2: 1_000_000.0,
        };
        c.validate()?;
        Ok(c)
    }
    pub fn validate(&self) -> Result<(), String> {
        if !self.tide_range_m.is_finite()
            || self.tide_range_m <= 0.0
            || self.tide_range_m > 20.0
            || !self.fetch_cap_m.is_finite()
            || !(32.0..=128_000.0).contains(&self.fetch_cap_m)
            || !self.supply_reach_m.is_finite()
            || !(32.0..=512.0).contains(&self.supply_reach_m)
            || !self.reference_area_m2.is_finite()
            || self.reference_area_m2 <= 0.0
        {
            return Err("invalid coast configuration".into());
        }
        Ok(())
    }
}
/// SPM 1984 eqs. 3-39/40, deep-water limit; UA = 0.71 U10^1.23 (SI).
/// Sustained wind, no swell, refraction, bottom friction or duration limitation.
pub fn wind_wave(u10: f64, fetch_m: f64) -> (f64, f64) {
    if u10 <= 0.0 || fetch_m <= 0.0 {
        return (0.0, 0.0);
    }
    let ua = 0.71 * u10.powf(1.23);
    let x = G * fetch_m / (ua * ua);
    (
        0.283 * ua * ua / G * (0.00565 * x.sqrt()).tanh(),
        7.54 * ua / G * (0.0379 * x.cbrt()).tanh(),
    )
}
/// Stockdon et al. 2006 general R2 expression. Natural sand beaches only.
pub fn runup_m(h: f64, period: f64, slope: f64) -> f64 {
    let hl = h * G * period * period / (2.0 * std::f64::consts::PI);
    1.1 * (0.35 * slope * hl.sqrt() + 0.5 * (hl * (0.563 * slope * slope + 0.004)).sqrt())
}
/// Index of potential supply, not a sediment flux. Maximum mouth, not summed
/// upstream/downstream samples of one channel. Area is the native effective area.
pub fn supply_index(area_m2: f64, erodibility: f64, c: &Config) -> f64 {
    let a = area_m2 * erodibility;
    a / (a + c.reference_area_m2)
}
#[derive(Clone, Copy, Debug, Default)]
pub struct Point {
    pub height_m: f64, // relative to sea datum
    pub slope: f64,    // tangent at this post, local 32 m gradient
    pub litho: u8,
    pub shore_m: f64,
    pub fixed: u8,
    /// Land-only path from a marine-touching shore anchor, with supratidal support.
    pub shore_path_m: f64,
    pub shore_connected: bool,
    pub tide_connected: bool,
    pub touches_marine: bool,
    pub mouth_area_m2: f64,
    pub wave_height_m: f64,
    pub wave_period_s: f64,
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Coast {
    pub class: u8,
    pub width_m: f32,
    pub supply: f32,
    pub wave_m: f32,
    pub slope: f32,
    pub flags: u16,
}
pub const INFERRED_SUBSTRATE: u16 = 1;
pub const MARSH_POTENTIAL: u16 = 2;
pub const DELTA_POTENTIAL: u16 = 4;
pub const WIDTH_CENSORED: u16 = 8;
pub const COARSE_FETCH: u16 = 16;
pub const FETCH_CAPPED: u16 = 32;
pub const NO_MARINE_BED: u16 = 64;
pub const ASSUMED_TIDE: u16 = 128;
pub const WIND_AS_U10: u16 = 256;
pub const UNSUPPORTED_SHORE: u16 = 512;
pub const CHANNEL_BED: u16 = 1024;
pub const TIDALLY_ISOLATED: u16 = 2048;

/// Balke et al. (2016), Fig. 5: pioneer elevation relative to MHW,
/// Pio_h(cm) = -108.23 log10(MTR(cm)) + 163.21. R²=0.39.
/// The configured effective tide is treated as MTR for this screening relation.
/// Do not clamp an out-of-frame prediction into an artificial marsh habitat.
pub fn marsh_lower_m(c: &Config) -> f64 {
    c.tide_range_m / 2.0 + (-108.23 * (100.0 * c.tide_range_m).log10() + 163.21) / 100.0
}

pub fn classify(p: Point, c: &Config) -> Result<Coast, String> {
    c.validate()?;
    if [
        p.height_m,
        p.slope,
        p.shore_m,
        p.shore_path_m,
        p.mouth_area_m2,
        p.wave_height_m,
        p.wave_period_s,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || p.shore_path_m < 0.0
        || p.slope < 0.0
        || p.mouth_area_m2 < 0.0
        || p.wave_height_m < 0.0
        || p.wave_period_s < 0.0
    {
        return Err("invalid/nonfinite coastal point".into());
    }
    let mut o = Coast {
        wave_m: p.wave_height_m as f32,
        slope: p.slope as f32,
        flags: ASSUMED_TIDE | WIND_AS_U10,
        ..Coast::default()
    };
    // Marine height is water surface, not bed. Do not invent intertidal seabed.
    if p.fixed == 2 {
        o.flags |= NO_MARINE_BED;
        return Ok(o);
    }
    if p.fixed == 1 || p.fixed == 5 {
        o.flags |= CHANNEL_BED;
        return Ok(o);
    }
    if p.fixed == 4 || p.shore_m < 0.0 || p.shore_m > 512.0 {
        return Ok(o);
    }
    if !p.shore_connected || p.shore_path_m > 512.0 {
        o.flags |= UNSUPPORTED_SHORE;
        return Ok(o);
    }
    if p.height_m < -c.tide_range_m / 2.0 {
        return Ok(o);
    }
    if p.height_m <= c.tide_range_m / 2.0 && !p.tide_connected {
        o.flags |= TIDALLY_ISOLATED;
        return Ok(o);
    }
    let k = *K.get(p.litho as usize).ok_or("unknown land lithology")?;
    let supply = supply_index(p.mouth_area_m2, k, c);
    o.supply = supply as f32;
    let cliff = 30_f64.to_radians().tan();
    let gravel = 10_f64.to_radians().tan();
    let flat = 3_f64.to_radians().tan();
    let rtr = if p.wave_height_m > 0.0 {
        c.tide_range_m / p.wave_height_m
    } else {
        f64::INFINITY
    };
    // No strength or sediment-stock measurements: these are eligibility rules.
    // NOAA ESI morphological descriptions, NOT an inverse grain-size measurement.
    let mobile = p.litho == 10
        || p.litho == 5
        || p.litho == 6
        || p.litho == 7
        || p.litho == 9
        || supply >= 0.5;
    let r2 = if p.slope < gravel {
        runup_m(p.wave_height_m, p.wave_period_s, p.slope)
    } else {
        p.wave_height_m
    };
    let raw_width = if p.slope > 0.0 {
        (c.tide_range_m + r2) / p.slope
    } else {
        f64::INFINITY
    };
    let width = raw_width.min(512.0);
    if raw_width > 512.0 {
        o.flags |= WIDTH_CENSORED;
    }
    // Rock-face annotation can represent the supratidal face at the actual
    // marine interface. Sediment and flat labels get no first-post exception.
    let face = p.slope >= cliff && p.touches_marine;
    if !face && (p.shore_path_m > width || p.height_m > c.tide_range_m / 2.0 + r2) {
        return Ok(o);
    }
    o.class = if face {
        1
    } else if !mobile {
        2
    } else if p.mouth_area_m2 >= c.reference_area_m2
        && p.touches_marine
        && p.slope <= flat
        && p.height_m <= c.tide_range_m / 2.0
    {
        o.flags |= DELTA_POTENTIAL;
        7
    } else if rtr > 15.0 && p.slope <= flat && p.height_m <= c.tide_range_m / 2.0 {
        if p.height_m >= marsh_lower_m(c) {
            o.flags |= MARSH_POTENTIAL;
            6
        } else {
            5
        }
    } else if p.slope >= gravel || p.litho == 7 {
        3
    } else {
        4
    };
    if mobile {
        o.flags |= INFERRED_SUBSTRATE;
    }
    o.width_m = if o.class == 1 { 0.0 } else { width as f32 };
    Ok(o)
}
