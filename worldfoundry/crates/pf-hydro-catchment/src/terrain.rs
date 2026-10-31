//! Grid and conditional geological relief adapted from retained Q20 pf-hydro-mr.
//! Source hashes and original revision are recorded in the lane BASE.json.
use crate::Result;
const PI: f64 = std::f64::consts::PI;
#[derive(Clone, Copy)]
pub struct Grid {
    pub ny: usize,
    pub nx: usize,
    pub radius: f64,
}
impl Grid {
    pub fn new(ny: usize, nx: usize, radius: f64) -> Result<Self> {
        if ny < 2
            || nx < 4
            || ny.checked_mul(nx).is_none_or(|n| n > 4_194_304)
            || !radius.is_finite()
            || radius <= 0.
        {
            return Err("unsupported grid dimensions/radius".into());
        }
        Ok(Self { ny, nx, radius })
    }
    pub fn len(&self) -> usize {
        self.ny * self.nx
    }
    pub fn is_empty(&self) -> bool {
        false
    }
    pub fn point(&self, cell: usize, u: f64, v: f64) -> [f64; 2] {
        [(cell % self.nx) as f64 + u, (cell / self.nx) as f64 + v].map_with_grid(self)
    }
    pub fn area(&self, cell: usize, v0: f64, v1: f64, du: f64) -> f64 {
        let y = (cell / self.nx) as f64;
        let lo = (y + v0) / self.ny as f64 * PI - PI / 2.;
        let hi = (y + v1) / self.ny as f64 * PI - PI / 2.;
        self.radius.powi(2) * 2. * PI / self.nx as f64 * du * (hi.sin() - lo.sin())
    }
    pub fn neighbor(&self, cell: usize, dx: isize, dy: isize) -> Option<usize> {
        let y = (cell / self.nx) as isize + dy;
        if y < 0 || y >= self.ny as isize {
            return None;
        }
        Some(
            y as usize * self.nx
                + ((cell % self.nx) as isize + dx).rem_euclid(self.nx as isize) as usize,
        )
    }
    pub fn sample(&self, a: &[f64], cell: usize, u: f64, v: f64) -> f64 {
        let x = (cell % self.nx) as f64 + u - 0.5;
        let y = ((cell / self.nx) as f64 + v - 0.5).clamp(0., self.ny as f64 - 1.);
        let ix = x.floor() as isize;
        let iy = y.floor() as usize;
        let tx = x - x.floor();
        let ty = y - iy as f64;
        let get = |x: isize, y: usize| {
            a[y.min(self.ny - 1) * self.nx + x.rem_euclid(self.nx as isize) as usize]
        };
        (1. - ty) * ((1. - tx) * get(ix, iy) + tx * get(ix + 1, iy))
            + ty * ((1. - tx) * get(ix, iy + 1) + tx * get(ix + 1, iy + 1))
    }
}
trait GridPoint {
    fn map_with_grid(self, g: &Grid) -> [f64; 2];
}
impl GridPoint for [f64; 2] {
    fn map_with_grid(self, g: &Grid) -> [f64; 2] {
        [
            self[0] / g.nx as f64 * 2. * PI - PI,
            self[1] / g.ny as f64 * PI - PI / 2.,
        ]
    }
}
pub fn distance(a: [f64; 2], b: [f64; 2], r: f64) -> f64 {
    let h = ((a[1] - b[1]) / 2.).sin().powi(2)
        + a[1].cos() * b[1].cos() * ((a[0] - b[0]) / 2.).sin().powi(2);
    2. * r * h.clamp(0., 1.).sqrt().asin()
}

fn hash(mut x: u64) -> f64 {
    x = x.wrapping_add(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    ((x ^ (x >> 31)) >> 11) as f64 / (1u64 << 53) as f64 * 2. - 1.
}
fn noise3(p: [f64; 3], seed: u64) -> f64 {
    let floor = p.map(|x| x.floor() as i64);
    let t = p.map(|x| {
        let t = x - x.floor();
        t * t * (3. - 2. * t)
    });
    let mut v = 0.;
    for k in 0..8 {
        let d = [k & 1, (k >> 1) & 1, (k >> 2) & 1];
        let mut key = seed;
        let mut w = 1.;
        for j in 0..3 {
            key ^= (floor[j] + d[j] as i64) as u64;
            key = key.wrapping_mul(0x9e3779b97f4a7c15).rotate_left(17);
            w *= if d[j] == 0 { 1. - t[j] } else { t[j] };
        }
        v += w * hash(key);
    }
    v
}
/// Q24 uses pinned endpoint elevation, not Q20 pre-pass reconstruction.
/// Conditional initial relief law at any planet address. No river geometry is sampled.
/// Native relief is mean-to-trunk, not RMS: 0.2 is an explicit global assumed RMS fraction.
pub fn terrain(input: &crate::input::Input, cell: usize, u: f64, v: f64) -> f64 {
    let g = input.grid;
    let p = g.point(cell, u, v);
    let xyz = [p[1].cos() * p[0].cos(), p[1].cos() * p[0].sin(), p[1].sin()];
    let mut r = 0.;
    for (scale, weight) in [(40000., 1.), (20000., 0.5), (10000., 0.25)] {
        r += weight * noise3(xyz.map(|x| x * g.radius / scale), input.seed);
    }
    g.sample(&input.elevation, cell, u, v) + g.sample(&input.relief, cell, u, v) * 0.2 * r / 0.4
}
