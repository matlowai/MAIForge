//! Rigid rotation on the sphere (brief K1/K2; acceptance: rigid rotation
//! preserves intra-plate distances; antimeridian/pole rotation sane).
//!
//! Conventions (declared, stable):
//! - Geocentric spherical coordinates: lat in [-90, 90] deg, lon in
//!   [-180, 180] deg. Unit vectors in a right-handed frame:
//!   `x = cos(lat)cos(lon), y = cos(lat)sin(lon), z = sin(lat)`.
//! - Rotation is the right-handed Rodrigues rotation about a UNIT axis by
//!   a signed angle in radians. An `EulerPoleV0` with positive
//!   `rate_deg_per_my` rotates positions counter-clockwise about its pole
//!   axis; a negative rate reverses (equivalently, rotates about the
//!   antipodal pole).
//! - `great_circle_rad` uses the atan2(|a x b|, a . b) form: accurate near
//!   0 and pi where acos(dot) loses precision.
//!
//! No silent clamping: lat outside [-90, 90] or lon outside [-180, 180] is
//! a loud contract error, not a wrapped coordinate.

use pf_contracts::domain::EulerPoleV0;

/// Loud error type for spherical helpers: every case names its field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RotationError {
    LatitudeOutOfRange { lat_deg: f64 },
    LongitudeOutOfRange { lon_deg: f64 },
    NonUnitAxis { norm: f64, expected: f64 },
    NonFinite { component: &'static str, value: f64 },
    DegeneratePoleAxis { lat_deg: f64, lon_deg: f64 },
}

impl std::fmt::Display for RotationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RotationError::LatitudeOutOfRange { lat_deg } => {
                write!(f, "lat_deg {lat_deg:?} outside [-90, 90]")
            }
            RotationError::LongitudeOutOfRange { lon_deg } => {
                write!(f, "lon_deg {lon_deg:?} outside [-180, 180]")
            }
            RotationError::NonUnitAxis { norm, expected } => {
                write!(
                    f,
                    "rotation axis norm {norm:?} != {expected:?} (unit required)"
                )
            }
            RotationError::NonFinite { component, value } => {
                write!(f, "nonfinite {component}: {value:?}")
            }
            RotationError::DegeneratePoleAxis { lat_deg, lon_deg } => {
                write!(
                    f,
                    "Euler pole at degenerate lat/lon ({lat_deg:?}, {lon_deg:?})"
                )
            }
        }
    }
}

impl std::error::Error for RotationError {}

const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

/// Geodetic-free geocentric lat/lon -> unit vector.
pub fn lat_lon_to_unit(lat_deg: f64, lon_deg: f64) -> Result<[f64; 3], RotationError> {
    if !lat_deg.is_finite() {
        return Err(RotationError::NonFinite {
            component: "lat_deg",
            value: lat_deg,
        });
    }
    if !lon_deg.is_finite() {
        return Err(RotationError::NonFinite {
            component: "lon_deg",
            value: lon_deg,
        });
    }
    if !(-90.0..=90.0).contains(&lat_deg) {
        return Err(RotationError::LatitudeOutOfRange { lat_deg });
    }
    if !(-180.0..=180.0).contains(&lon_deg) {
        return Err(RotationError::LongitudeOutOfRange { lon_deg });
    }
    let lat = lat_deg * DEG_TO_RAD;
    let lon = lon_deg * DEG_TO_RAD;
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_lon, cos_lon) = lon.sin_cos();
    Ok([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])
}

/// Unit vector -> (lat_deg, lon_deg). Exact at poles (lon = 0 by
/// convention there; the input vector's tiny x/y noise is not a contract
/// concern because pole inputs reproduce exactly).
pub fn unit_to_lat_lon(p: [f64; 3]) -> Result<(f64, f64), RotationError> {
    for v in p.iter() {
        if !v.is_finite() {
            return Err(RotationError::NonFinite {
                component: "unit vector component",
                value: *v,
            });
        }
    }
    let norm = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    if !(0.99..=1.01).contains(&norm) {
        return Err(RotationError::NonUnitAxis {
            norm,
            expected: 1.0,
        });
    }
    let z = (p[2] / norm).clamp(-1.0, 1.0);
    let lat = z.asin() / DEG_TO_RAD;
    let lon = p[1].atan2(p[0]) / DEG_TO_RAD;
    Ok((lat, lon))
}

/// Right-handed Rodrigues rotation of `p` about UNIT axis `k` by signed
/// `angle_rad`. Numerically standard form; precision documented as f64
/// roundoff (tests hold it to 1e-12 relative on unit norms).
pub fn rotate_about_axis(
    p: [f64; 3],
    k: [f64; 3],
    angle_rad: f64,
) -> Result<[f64; 3], RotationError> {
    Ok(AxisRotation::new(k, angle_rad)?.apply(p))
}

/// SEED-PERF (lane SP, 2026-10-13): a validated Rodrigues rotation with its
/// axis and `sin_cos` evaluated once. `apply` is the exact arithmetic of
/// `rotate_about_axis` (which is now `new` + `apply`), so a rotation built
/// once per plate and applied per cell is bit for bit the per-call form.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisRotation {
    k: [f64; 3],
    s: f64,
    c: f64,
}

impl AxisRotation {
    /// Validates exactly as `rotate_about_axis` did (finite angle, unit axis).
    pub fn new(k: [f64; 3], angle_rad: f64) -> Result<Self, RotationError> {
        if !angle_rad.is_finite() {
            return Err(RotationError::NonFinite {
                component: "angle_rad",
                value: angle_rad,
            });
        }
        let norm = (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]).sqrt();
        if (norm - 1.0).abs() > 1.0e-12 {
            return Err(RotationError::NonUnitAxis {
                norm,
                expected: 1.0,
            });
        }
        let (s, c) = angle_rad.sin_cos();
        Ok(Self { k, s, c })
    }

    /// Rotate `p` (no validation of `p`, as before).
    #[inline]
    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        let (k, s, c) = (self.k, self.s, self.c);
        let cross = [
            k[1] * p[2] - k[2] * p[1],
            k[2] * p[0] - k[0] * p[2],
            k[0] * p[1] - k[1] * p[0],
        ];
        let dot = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
        [
            p[0] * c + cross[0] * s + k[0] * dot * (1.0 - c),
            p[1] * c + cross[1] * s + k[1] * dot * (1.0 - c),
            p[2] * c + cross[2] * s + k[2] * dot * (1.0 - c),
        ]
    }

    /// The prepared form of `rotate_by_euler(_, pole, dt_my)`: same checks,
    /// same errors (they never depend on the rotated point).
    pub fn euler(pole: &EulerPoleV0, dt_my: f64) -> Result<Self, RotationError> {
        if !dt_my.is_finite() || dt_my < 0.0 {
            return Err(RotationError::NonFinite {
                component: "dt_my (must be finite >= 0)",
                value: dt_my,
            });
        }
        let axis = euler_pole_axis(pole)?;
        let angle_rad = pole.rate_deg_per_my * dt_my * DEG_TO_RAD;
        Self::new(axis, angle_rad)
    }
}

/// Great-circle angular distance in radians between two unit vectors,
/// atan2 form (accurate at 0 and pi).
pub fn great_circle_rad(a: [f64; 3], b: [f64; 3]) -> Result<f64, RotationError> {
    for p in [a, b] {
        for v in p.iter() {
            if !v.is_finite() {
                return Err(RotationError::NonFinite {
                    component: "vector component",
                    value: *v,
                });
            }
        }
    }
    let cross_norm = {
        let c = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
    };
    let dot = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]).clamp(-1.0, 1.0);
    Ok(cross_norm.atan2(dot))
}

/// Unit axis of an Euler pole (its lat/lon as a point on the sphere).
pub fn euler_pole_axis(pole: &EulerPoleV0) -> Result<[f64; 3], RotationError> {
    pole.validate()
        .map_err(|_| RotationError::DegeneratePoleAxis {
            lat_deg: pole.lat_deg,
            lon_deg: pole.lon_deg,
        })?;
    lat_lon_to_unit(pole.lat_deg, pole.lon_deg)
}

/// Advance one position by a plate's Euler rotation over `dt_my`:
/// angle = rate_deg_per_my * dt_my (signed), right-handed about the pole.
pub fn rotate_by_euler(
    p: [f64; 3],
    pole: &EulerPoleV0,
    dt_my: f64,
) -> Result<[f64; 3], RotationError> {
    Ok(AxisRotation::euler(pole, dt_my)?.apply(p))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1.0e-12;

    #[test]
    fn lat_lon_roundtrip_on_cardinal_points() {
        for (lat, lon) in [
            (0.0, 0.0),
            (0.0, 90.0),
            (0.0, -179.9),
            (45.0, 30.0),
            (90.0, 0.0),
            (-90.0, 0.0),
            (12.5, -74.3),
        ] {
            let u = lat_lon_to_unit(lat, lon).unwrap();
            let (lat2, lon2) = unit_to_lat_lon(u).unwrap();
            assert!((lat - lat2).abs() < TOL, "lat {lat} -> {lat2}");
            // lon compares modulo wrap for points exactly at the branch cut
            let dlon = (lon - lon2 + 540.0).rem_euclid(360.0) - 180.0;
            assert!(dlon.abs() < TOL, "lon {lon} -> {lon2}");
        }
    }

    #[test]
    fn out_of_range_coordinates_fail_loudly() {
        assert!(matches!(
            lat_lon_to_unit(91.0, 0.0),
            Err(RotationError::LatitudeOutOfRange { .. })
        ));
        assert!(matches!(
            lat_lon_to_unit(0.0, 181.0),
            Err(RotationError::LongitudeOutOfRange { .. })
        ));
        assert!(matches!(
            lat_lon_to_unit(f64::NAN, 0.0),
            Err(RotationError::NonFinite { .. })
        ));
    }

    #[test]
    fn rotation_by_zero_is_identity() {
        let p = lat_lon_to_unit(33.0, -117.0).unwrap();
        let q = rotate_about_axis(p, [0.0, 0.0, 1.0], 0.0).unwrap();
        for i in 0..3 {
            assert!((p[i] - q[i]).abs() < TOL);
        }
    }

    #[test]
    fn known_quarter_turns() {
        // +90 deg about +z: (1,0,0) -> (0,1,0); (0,1,0) -> (-1,0,0).
        let q = rotate_about_axis(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        assert!(
            (q[0] - 0.0).abs() < TOL && (q[1] - 1.0).abs() < TOL && (q[2] - 0.0).abs() < TOL,
            "{q:?}"
        );
        let q = rotate_about_axis(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        assert!(
            (q[0] + 1.0).abs() < TOL && (q[1] - 0.0).abs() < TOL,
            "{q:?}"
        );
        // +90 deg about +x: (0,1,0) -> (0,0,1) (right-handed).
        let q = rotate_about_axis(
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        assert!((q[2] - 1.0).abs() < TOL, "{q:?}");
    }

    #[test]
    fn rotation_preserves_norm_and_distances() {
        let k = lat_lon_to_unit(20.0, 65.0).unwrap();
        let mut pts = Vec::new();
        for i in 0..40 {
            pts.push(
                lat_lon_to_unit(-80.0 + 4.0 * i as f64, (17.0 * i as f64) % 360.0 - 180.0).unwrap(),
            );
        }
        let rotated: Vec<_> = pts
            .iter()
            .map(|p| rotate_about_axis(*p, k, 0.77).unwrap())
            .collect();
        for (p, q) in pts.iter().zip(&rotated) {
            let np = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let nq = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
            assert!((np - nq).abs() < TOL, "norm changed {np} -> {nq}");
        }
        for i in 0..pts.len() {
            for j in (i + 1)..pts.len() {
                let d0 = great_circle_rad(pts[i], pts[j]).unwrap();
                let d1 = great_circle_rad(rotated[i], rotated[j]).unwrap();
                assert!((d0 - d1).abs() < TOL, "pair ({i},{j}) distance drifted");
            }
        }
    }

    #[test]
    fn great_circle_extremes() {
        let north = lat_lon_to_unit(90.0, 0.0).unwrap();
        let south = lat_lon_to_unit(-90.0, 0.0).unwrap();
        let eq0 = lat_lon_to_unit(0.0, 0.0).unwrap();
        let eq90 = lat_lon_to_unit(0.0, 90.0).unwrap();
        assert!((great_circle_rad(north, south).unwrap() - std::f64::consts::PI).abs() < TOL);
        assert!(great_circle_rad(north, north).unwrap().abs() < TOL);
        assert!((great_circle_rad(eq0, eq90).unwrap() - std::f64::consts::FRAC_PI_2).abs() < TOL);
    }

    #[test]
    fn euler_rotation_moves_longitude_along_latitude_circle_about_z() {
        // A z-axis pole rotates points east/west along their latitude circle.
        let pole = EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: 1.0,
        };
        let p = lat_lon_to_unit(30.0, 10.0).unwrap();
        let q = rotate_by_euler(p, &pole, 25.0).unwrap();
        let (lat, lon) = unit_to_lat_lon(q).unwrap();
        assert!(
            (lat - 30.0).abs() < TOL,
            "latitude must be preserved, got {lat}"
        );
        let dlon = (lon - 35.0 + 540.0).rem_euclid(360.0) - 180.0;
        assert!(
            dlon.abs() < TOL,
            "longitude must advance by rate*dt = 25 deg, got {lon}"
        );
    }

    #[test]
    fn antimeridian_and_pole_positions_rotate_safely() {
        let pole = EulerPoleV0 {
            lat_deg: -15.0,
            lon_deg: 172.0,
            rate_deg_per_my: 3.3,
        };
        for (lat, lon) in [
            (0.0, 179.9),
            (0.0, -179.9),
            (5.0, 180.0 - 1.0e-9),
            (90.0, 0.0),
            (-90.0, 45.0),
            (89.9, -0.1),
        ] {
            let p = lat_lon_to_unit(lat, lon).unwrap();
            let q = rotate_by_euler(p, &pole, 12.5).unwrap();
            let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
            assert!((n - 1.0).abs() < TOL, "norm lost at ({lat},{lon})");
            // And rotating the full history back by -dt returns the start.
            let back_pole = EulerPoleV0 {
                rate_deg_per_my: -pole.rate_deg_per_my,
                ..pole
            };
            let r = rotate_by_euler(q, &back_pole, 12.5).unwrap();
            let d = great_circle_rad(p, r).unwrap();
            assert!(d < 1.0e-12, "round trip drifted {d} at ({lat},{lon})");
        }
    }

    #[test]
    fn non_unit_axis_is_refused() {
        let err = rotate_about_axis([1.0, 0.0, 0.0], [0.0, 0.0, 2.0], 0.5).unwrap_err();
        assert!(matches!(err, RotationError::NonUnitAxis { .. }), "{err}");
        let err = rotate_by_euler(
            [1.0, 0.0, 0.0],
            &EulerPoleV0 {
                lat_deg: f64::NAN,
                lon_deg: 0.0,
                rate_deg_per_my: 1.0,
            },
            1.0,
        )
        .unwrap_err();
        assert!(
            matches!(err, RotationError::DegeneratePoleAxis { .. }),
            "{err}"
        );
    }
}
