//! native pipeline: geological column tilt is measured from its T0 cell centre,
//! never from the requested chart or a site. Tangent chord offsets in metres
//! are nonsingular at poles and wrap continuously across the antimeridian.
pub fn offset(
    rows: usize,
    cols: usize,
    radius_m: f64,
    cell: u32,
    lat: f64,
    lon: f64,
) -> (f64, f64) {
    let phi = (-90.0 + (cell as usize / cols) as f64 * 180.0 / rows as f64 + 90.0 / rows as f64)
        .to_radians();
    let lam = (-180.0 + (cell as usize % cols) as f64 * 360.0 / cols as f64 + 180.0 / cols as f64)
        .to_radians();
    let p = lat.to_radians();
    let dl = lon.to_radians() - lam;
    (
        radius_m * p.cos() * dl.sin(),
        radius_m * (phi.cos() * p.sin() - phi.sin() * p.cos() * dl.cos()),
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn centre_and_dateline() {
        let (x, y) = offset(288, 576, 6371000.0, 0, -89.6875, -179.6875);
        assert!(x.abs() < 1e-8 && y.abs() < 1e-8);
        let a = offset(288, 576, 6371000.0, 288 * 576 - 1, 89.9, -180.0);
        let b = offset(288, 576, 6371000.0, 288 * 576 - 1, 89.9, 180.0);
        assert!((a.0 - b.0).abs() < 1e-7 && (a.1 - b.1).abs() < 1e-7);
        let pole = offset(288, 576, 6371000.0, 288 * 576 - 1, 90.0, 0.0);
        assert!(pole.0.is_finite() && pole.1.is_finite());
    }
}
