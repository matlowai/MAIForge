//! native pipeline initial-condition reconstruction, before routing or erosion.
//! Native finite-volume elevations are conditionally interpreted as centre
//! observations. Separable monotone cubic Hermite interpolation retains those
//! observations and stays inside the four bracketing values. Unlike bilinear
//! interpolation it does not put a slope discontinuity at every native centre.
//! This changes the declared initial-condition basis, not solved tile boundaries,
//! geological class boundaries, water stages, or post-erosion geometry.
use pf_hydro_catchment::terrain::Grid;
pub const LAW: &str = "native-centre-observations/monotone-c1/1";
fn tangent(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else {
        2.0 * a * b / (a + b)
    }
}
fn interpolate(p: [f64; 4], t: f64) -> f64 {
    let m0 = tangent(p[1] - p[0], p[2] - p[1]);
    let m1 = tangent(p[2] - p[1], p[3] - p[2]);
    let t2 = t * t;
    let t3 = t2 * t;
    (2.0 * t3 - 3.0 * t2 + 1.0) * p[1]
        + (t3 - 2.0 * t2 + t) * m0
        + (-2.0 * t3 + 3.0 * t2) * p[2]
        + (t3 - t2) * m1
}
pub fn elevation(g: Grid, z: &[f64], cell: usize, u: f64, v: f64) -> f64 {
    let x = (cell % g.nx) as f64 + u - 0.5;
    let y = ((cell / g.nx) as f64 + v - 0.5).clamp(0.0, g.ny as f64 - 1.0);
    let i = x.floor() as isize;
    let j = y.floor() as isize;
    let at = |x: isize, y: isize| {
        z[y.clamp(0, g.ny as isize - 1) as usize * g.nx + x.rem_euclid(g.nx as isize) as usize]
    };
    let row = |y| {
        interpolate(
            [at(i - 1, y), at(i, y), at(i + 1, y), at(i + 2, y)],
            x - x.floor(),
        )
    };
    interpolate([row(j - 1), row(j), row(j + 1), row(j + 2)], y - y.floor())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn initial_field_preserves_observations_and_slope_without_overshoot() {
        let g = Grid::new(4, 8, 6371000.0).unwrap();
        let z: Vec<_> = (0..g.len()).map(|i| ((i * 7) % 11) as f64).collect();
        for c in 0..g.len() {
            assert_eq!(elevation(g, &z, c, 0.5, 0.5), z[c]);
        }
        let p = [-100.0, 1.0, 2.0, 100.0];
        for k in 0..101 {
            assert!((1.0..=2.0).contains(&interpolate(p, k as f64 / 100.0)));
        }
        let p = [0.0, 1.0, 3.0, 4.0, 10.0];
        let e = 1e-6;
        let left = (p[2] - interpolate([p[0], p[1], p[2], p[3]], 1.0 - e)) / e;
        let right = (interpolate([p[1], p[2], p[3], p[4]], e) - p[2]) / e;
        assert!((left - right).abs() < 1e-5);
    }
}
