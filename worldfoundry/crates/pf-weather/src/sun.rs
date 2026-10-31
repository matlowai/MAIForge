//! Calendar and solar geometry over the world hour index.
//!
//! The world hour index `t` (f64 hours, integer = top of the hour) counts from
//! the world epoch: hour 0 is 00:00 at longitude 0 on day 0 of year 0 of a
//! 365-day calendar. The kernel's `WorldTime` tick has no calendar meaning yet
//! (`world/crates/vz-protocol/src/time.rs`), so the tick -> hour mapping is an
//! adapter decision that stays OPEN; this crate only takes hours.
//!
//! Declination: Cooper (1969) `delta = tilt * sin(2 pi (284 + n) / 365)`
//! (recalled), with the world's axial tilt.

use std::f64::consts::PI;

pub const DEG: f64 = PI / 180.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarClock {
    /// Local solar hour in [0, 24).
    pub solar_hour: f64,
    /// Local calendar day index since the epoch (integer days).
    pub local_day: i64,
    /// Day of year of the local day, 0..365.
    pub doy: f64,
    pub declination_rad: f64,
    /// Sunrise and sunset in local solar hours (polar day: 0 / 24, polar night: 12 / 12).
    pub sunrise_h: f64,
    pub sunset_h: f64,
    /// Sun elevation (deg) now and at local noon.
    pub elevation_deg: f64,
    pub noon_elevation_deg: f64,
}

pub fn solar_clock(t_h: f64, lat_deg: f64, lon_deg: f64, tilt_deg: f64) -> SolarClock {
    let local = t_h + lon_deg / 15.0;
    let local_day = (local / 24.0).floor() as i64;
    let solar_hour = local - local_day as f64 * 24.0;
    let doy = local_day.rem_euclid(365) as f64;
    let decl = tilt_deg * DEG * (2.0 * PI * (284.0 + doy + 1.0) / 365.0).sin();
    let phi = lat_deg * DEG;
    let x = -(phi.tan() * decl.tan());
    let h0 = if x <= -1.0 { PI } else if x >= 1.0 { 0.0 } else { x.acos() };
    let half_day_h = h0 / (15.0 * DEG);
    let ha = (solar_hour - 12.0) * 15.0 * DEG;
    let sin_e = phi.sin() * decl.sin() + phi.cos() * decl.cos() * ha.cos();
    let noon = 90.0 - (lat_deg - decl / DEG).abs();
    SolarClock {
        solar_hour,
        local_day,
        doy,
        declination_rad: decl,
        sunrise_h: 12.0 - half_day_h,
        sunset_h: 12.0 + half_day_h,
        elevation_deg: sin_e.clamp(-1.0, 1.0).asin() / DEG,
        noon_elevation_deg: noon,
    }
}

/// World hour of local solar sunrise on a local day at a place.
pub fn sunrise_world_hour(local_day: i64, lat_deg: f64, lon_deg: f64, tilt_deg: f64) -> f64 {
    // Evaluate the clock at local noon of that day, then step back.
    let noon_t = local_day as f64 * 24.0 + 12.0 - lon_deg / 15.0;
    let c = solar_clock(noon_t, lat_deg, lon_deg, tilt_deg);
    noon_t - (12.0 - c.sunrise_h)
}

/// World hour of local solar sunset on a local day at a place.
pub fn sunset_world_hour(local_day: i64, lat_deg: f64, lon_deg: f64, tilt_deg: f64) -> f64 {
    let noon_t = local_day as f64 * 24.0 + 12.0 - lon_deg / 15.0;
    let c = solar_clock(noon_t, lat_deg, lon_deg, tilt_deg);
    noon_t + (c.sunset_h - 12.0)
}
