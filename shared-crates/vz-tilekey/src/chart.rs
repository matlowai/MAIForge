//! `local_oblique_stereographic_v1`: the one chart of a generated world (ADR
//! 044 section 4, option B). Spherical oblique stereographic projection,
//! centred on the declared origin, scale factor 1 at the origin, x east and y
//! north in integer millimetres. f64 internally (with `vz.detmath/v1` trig, so
//! the JS twin gives the same bits), quantized once to integer mm with
//! `nearest_mm_half_up_v1` (floor(v + 0.5)).
//!
//! Formulas (sphere of radius R, origin (phi1, lam0)):
//! - forward: with dl = lam - lam0 normalized to (-pi, pi],
//!   e = cos(phi) sin(dl), n = cos(phi1) sin(phi) - sin(phi1) cos(phi) cos(dl),
//!   u = sin(phi1) sin(phi) + cos(phi1) cos(phi) cos(dl);
//!   x = 2R e / (1 + u), y = 2R n / (1 + u).
//! - inverse (closed form, no asin, well conditioned at the poles): with
//!   a = x / 2R, b = y / 2R, q = a^2 + b^2, d = 1 + q:
//!   local unit vector (east, north, up) = (2a/d, 2b/d, (1 - q)/d); rotate by
//!   phi1 into the lam0 meridian frame; phi = atan2(pz, sqrt(px^2 + py^2)),
//!   lam = lam0 + atan2(py, px), normalized to (-pi, pi].
//!
//! The chart's projection domain is the CLOSED square |x|, |y| <= 131.072 km
//! (edge-inclusive posts reach the boundary); key ownership is half-open
//! (see `TileKeyV2::of_point`).

// Negated comparisons deliberately send NaN to the refusal branch; literals
// are written to the digit so the JS twin parses identical doubles.
#![allow(clippy::neg_cmp_op_on_partial_ord, clippy::excessive_precision, clippy::approx_constant)]

use crate::detmath::{atan2, cos, sin, PI};
use crate::{TileKeyError, CHART_HALF_EXTENT_MM};
use serde::{Deserialize, Serialize};

pub const CHART_KIND_V1: &str = "local_oblique_stereographic_v1";
pub const QUANTIZATION_V1: &str = "nearest_mm_half_up_v1";
pub const TWO_PI: f64 = 6.283_185_307_179_586_232e0;
/// The f64 nearest pi/180 (0.017453292519943295).
pub const DEG2RAD: f64 = 1.745_329_251_994_329_547_4e-2;
/// The f64 nearest 180/pi (57.29577951308232).
pub const RAD2DEG: f64 = 5.729_577_951_308_232_286_5e1;

/// Same field names and integer types as vz-protocol's `GroundProjection`, so
/// the declaration can ride in a ground recipe's `projection` field (ADR 044
/// amendment C5).
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct ChartDeclV1 {
    pub kind: String,
    pub origin_lat_e7: i64,
    pub origin_lon_e7: i64,
    pub datum_radius_mm: i64,
}

impl ChartDeclV1 {
    pub fn new(origin_lat_e7: i64, origin_lon_e7: i64, datum_radius_mm: i64) -> Self {
        ChartDeclV1 { kind: CHART_KIND_V1.into(), origin_lat_e7, origin_lon_e7, datum_radius_mm }
    }

    /// Same bounds as the existing ground projection checks: origin within 89
    /// degrees latitude and 180 degrees longitude, radius 1..=10^13 mm. The
    /// radius must also be at least 10x the chart half extent.
    pub fn problem(&self) -> Option<&'static str> {
        if self.kind != CHART_KIND_V1 { return Some("unsupported_chart_kind"); }
        if !(-890_000_000..=890_000_000).contains(&self.origin_lat_e7)
            || !(-1_800_000_000..=1_800_000_000).contains(&self.origin_lon_e7) {
            return Some("chart_origin_out_of_range");
        }
        if !(10 * CHART_HALF_EXTENT_MM..=10_000_000_000_000).contains(&self.datum_radius_mm) {
            return Some("chart_radius_out_of_range");
        }
        None
    }
}

/// Integer degrees times 10^k to radians: (v / 10^k) * DEG2RAD.
pub fn deg_scaled_to_rad(v: i64, scale: f64) -> f64 { (v as f64 / scale) * DEG2RAD }
pub fn e7_to_rad(v: i64) -> f64 { deg_scaled_to_rad(v, 1e7) }
pub fn e9_to_rad(v: i64) -> f64 { deg_scaled_to_rad(v, 1e9) }
/// Radians to nano-degrees, nearest half up: floor(rad * RAD2DEG * 1e9 + 0.5).
pub fn rad_to_e9(r: f64) -> i64 { (((r * RAD2DEG) * 1e9) + 0.5).floor() as i64 }
/// `nearest_mm_half_up_v1`.
pub fn quantize_mm(v: f64) -> i64 { (v + 0.5).floor() as i64 }

#[derive(Clone, Debug)]
pub struct Chart {
    pub decl: ChartDeclV1,
    lat0: f64,
    lon0: f64,
    s1: f64,
    c1: f64,
    r: f64,
    two_r: f64,
}

impl Chart {
    pub fn new(decl: ChartDeclV1) -> Result<Chart, TileKeyError> {
        if let Some(p) = decl.problem() { return Err(TileKeyError(p)); }
        let lat0 = e7_to_rad(decl.origin_lat_e7);
        let lon0 = e7_to_rad(decl.origin_lon_e7);
        let r = decl.datum_radius_mm as f64;
        Ok(Chart { s1: sin(lat0), c1: cos(lat0), lat0, lon0, r, two_r: 2.0 * r, decl })
    }

    pub fn origin_rad(&self) -> (f64, f64) { (self.lat0, self.lon0) }
    pub fn radius_mm(&self) -> f64 { self.r }

    /// Unquantized chart mm. `lat` in [-pi/2, pi/2], `lon` in [-pi, pi].
    pub fn forward_f64(&self, lat: f64, lon: f64) -> Result<[f64; 2], TileKeyError> {
        if !(lat.abs() <= PI / 2.0 + 1e-15) || !(lon.abs() <= PI + 1e-15) {
            return Err(TileKeyError("geodetic_out_of_range"));
        }
        let mut dl = lon - self.lon0;
        if dl > PI { dl -= TWO_PI; }
        if dl <= -PI { dl += TWO_PI; }
        let (sp, cp) = (sin(lat), cos(lat));
        let (sd, cd) = (sin(dl), cos(dl));
        let e = cp * sd;
        let n = self.c1 * sp - (self.s1 * cp) * cd;
        let u = self.s1 * sp + (self.c1 * cp) * cd;
        let den = 1.0 + u;
        if !(den > 1e-6) { return Err(TileKeyError("outside_chart")); }
        let x = (self.two_r * e) / den;
        let y = (self.two_r * n) / den;
        Ok([x, y])
    }

    /// Chart mm, quantized once; refuses points outside the closed domain.
    pub fn forward_mm(&self, lat: f64, lon: f64) -> Result<[i64; 2], TileKeyError> {
        let [x, y] = self.forward_f64(lat, lon)?;
        let h = CHART_HALF_EXTENT_MM as f64;
        if !(x.abs() <= h + 0.5) || !(y.abs() <= h + 0.5) { return Err(TileKeyError("outside_chart")); }
        let q = [quantize_mm(x), quantize_mm(y)];
        if q.iter().any(|v| v.abs() > CHART_HALF_EXTENT_MM) { return Err(TileKeyError("outside_chart")); }
        Ok(q)
    }

    /// Geodetic (lat, lon) in radians of a chart point (f64 mm; integer mm
    /// callers pass `as f64`). Refuses points outside the closed domain.
    pub fn inverse_f64(&self, x: f64, y: f64) -> Result<(f64, f64), TileKeyError> {
        let h = CHART_HALF_EXTENT_MM as f64;
        if !(x.abs() <= h) || !(y.abs() <= h) { return Err(TileKeyError("outside_chart")); }
        if x == 0.0 && y == 0.0 { return Ok((self.lat0, self.lon0)); }
        let a = x / self.two_r;
        let b = y / self.two_r;
        let q = a * a + b * b;
        let d = 1.0 + q;
        let ve = (2.0 * a) / d;
        let vn = (2.0 * b) / d;
        let vu = (1.0 - q) / d;
        let px = self.c1 * vu - self.s1 * vn;
        let py = ve;
        let pz = self.s1 * vu + self.c1 * vn;
        let lat = atan2(pz, (px * px + py * py).sqrt());
        let mut lon = self.lon0 + atan2(py, px);
        if lon > PI { lon -= TWO_PI; }
        if lon <= -PI { lon += TWO_PI; }
        Ok((lat, lon))
    }

    pub fn inverse_mm(&self, x_mm: i64, y_mm: i64) -> Result<(f64, f64), TileKeyError> {
        self.inverse_f64(x_mm as f64, y_mm as f64)
    }

    /// Point scale factor of the projection at a chart point: 1 + rho^2/(4R^2).
    pub fn scale_factor(&self, x: f64, y: f64) -> f64 {
        let a = x / self.two_r;
        let b = y / self.two_r;
        1.0 + (a * a + b * b)
    }
}
