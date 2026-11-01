//! World-anchored first-touch relief, native pipeline / model design.
//!
//! This is a versioned reduced-order construction, not the chart-conditioned
//! RIVERNET realization and not simulated erosion history. The global T0
//! receiver graph supplies drainage segments. Published elevation is the
//! long-wavelength datum. Subgrid relief and the native lithology/creep laws
//! supply the coarse boundary field. World-addressed cell tributaries and
//! native hillslopes refine it (planet_cell). Draws use world seed and physical
//! cell/node identity. No chart, site or request order enters construction.
//!
//! Cross section: q = E*x, q = D*S/(1-(S/Sc)^2), where x is distance from
//! a divide, E is m/My, D is m^2/My, Sc is dimensionless. Integrate S from
//! divide to channel, choosing half-width L so relief H(L) = published R.
//! A segment lowers the datum by H(max(L-distance,0)). The union of cuts
//! is their maximum, not a blend of separately solved chart heights.
//! The coarse field supplies shared boundary conditions, not final fine
//! terrain. Cell solve receipts report mean residuals. No meander evolution
//! or sediment transport is claimed.
//! No material is canonically moved by constructing this candidate surface.
use crate::R;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
/// Shared sampling projection, also consumed by the detail tier.
pub use vz_tilekey::chart::{Chart, RAD2DEG};

pub use crate::planet_catchment::{LIMITS, MODEL};
const BIN_M: f64 = 100_000.0;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Laws {
    pub stream_power_k_per_my: f64,
    pub diffusivity_m2_per_my: f64,
    pub channel_head_length_m: f64,
}

pub struct Inputs {
    pub seed: u64,
    pub rows: usize,
    pub cols: usize,
    pub radius_m: f64,
    pub elevation_m: Vec<f32>,
    pub relief_m: Vec<f32>,
    pub receiver: Vec<u32>,
    pub land: Vec<bool>,
    pub lithology: Vec<u8>,
    pub laws: Laws,
}

#[derive(Clone, Copy, Debug)]
struct Segment {
    a: [f64; 3],
    tangent: [f64; 3],
    angle: f64,
    half_width_m: f64,
    erosion_m_per_my: f64,
    critical_slope: f64,
    from: u32,
    to: u32,
}

pub struct WorldRelief {
    pub input: Inputs,
    joint: Option<crate::planet_catchment::World>,
    segments: Vec<Segment>,
    axis: Vec<Option<usize>>,
    bins: HashMap<[i32; 3], Vec<u32>>,
    poles: [f64; 2],
    cells: std::sync::Mutex<HashMap<u32, std::sync::Arc<crate::planet_cell::CellSurface>>>,
}
#[derive(Clone, Copy, Debug, Serialize)]
pub struct Sample {
    pub height_m: f64,
    pub datum_m: f64,
    pub incision_m: f64,
    pub receiver_segment: Option<[u32; 2]>,
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn unit(lat: f64, lon: f64) -> [f64; 3] {
    let (p, l) = (lat.to_radians(), lon.to_radians());
    if lat.abs() == 90.0 {
        return [0.0, 0.0, lat.signum()];
    }
    [p.cos() * l.cos(), p.cos() * l.sin(), p.sin()]
}
fn bin(p: [f64; 3], radius: f64) -> [i32; 3] {
    p.map(|v| (v * radius / BIN_M).floor() as i32)
}

/// Integral of the native nonlinear-creep slope from the divide to x.
/// Stable at both zero erosion and arbitrarily small x; no numerical PDE solve.
pub fn creep_height(x: f64, e: f64, d: f64, sc: f64) -> f64 {
    if x <= 0.0 || e == 0.0 {
        return 0.0;
    }
    let b = 2.0 * e / (d * sc);
    let t = b * x;
    let u_minus_one = t * t / (1.0 + (1.0 + t * t).sqrt());
    sc / b * (u_minus_one - (u_minus_one * 0.5).ln_1p())
}
fn half_width(relief: f64, e: f64, d: f64, sc: f64) -> f64 {
    if relief <= 0.0 || e <= 0.0 {
        return 0.0;
    }
    let mut lo = 0.0;
    let mut hi = 2.0 * (relief / sc + (2.0 * d * relief / e).sqrt());
    for _ in 0..36 {
        let mid = (lo + hi) * 0.5;
        if creep_height(mid, e, d, sc) < relief {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) * 0.5
}
impl WorldRelief {
    pub fn new(input: Inputs) -> R<Self> {
        let n = input
            .rows
            .checked_mul(input.cols)
            .ok_or("planet_grid_overflow")?;
        if input.rows < 2
            || input.cols < 4
            || !(input.radius_m.is_finite() && input.radius_m > 0.0)
            || input.elevation_m.len() != n
            || input.relief_m.len() != n
            || input.receiver.len() != n
            || input.land.len() != n
            || input.lithology.len() != n
        {
            return Err("planet_relief_grid".into());
        }
        let l = &input.laws;
        if [
            l.stream_power_k_per_my,
            l.diffusivity_m2_per_my,
            l.channel_head_length_m,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
            || input.elevation_m.iter().any(|v| !v.is_finite())
            || input.relief_m.iter().any(|v| !v.is_finite() || *v < 0.0)
            || input
                .receiver
                .iter()
                .any(|&v| v != u32::MAX && v as usize >= n)
            || input.lithology.iter().any(|&v| v >= 12)
        {
            return Err("planet_relief_inputs".into());
        }
        let pole = |row: usize| {
            input.elevation_m[row * input.cols..(row + 1) * input.cols]
                .iter()
                .map(|&v| v as f64)
                .sum::<f64>()
                / input.cols as f64
        };
        let poles = [pole(0), pole(input.rows - 1)];
        let mut world = Self {
            input,
            joint: None,
            segments: Vec::new(),
            axis: vec![None; n],
            bins: HashMap::new(),
            poles,
            cells: std::sync::Mutex::new(HashMap::new()),
        };
        let rocks = crate::t1b_litho::LithologyV1::standin_v1();
        for c in 0..n {
            let r = world.input.receiver[c];
            if r == u32::MAX || r as usize == c || !world.input.land[c] {
                continue;
            }
            let (la, lo) = world.centre(c);
            let (lb, ln) = world.centre(r as usize);
            let a = unit(la, lo);
            let z = unit(lb, ln);
            let cos = dot(a, z).clamp(-1.0, 1.0);
            let angle = cos.acos();
            if angle < 1e-12 {
                continue;
            }
            if angle > std::f64::consts::FRAC_PI_2 {
                return Err(format!("planet_nonlocal_receiver:{c}:{r}"));
            }
            let tangent = std::array::from_fn(|i| (z[i] - cos * a[i]) / angle.sin());
            let relief = world.input.relief_m[c] as f64;
            // Same E0 closure as t1b::closure_rate_m_per_my, using the
            // physical spherical cell area, not a chart's measured hat mass.
            let dlat = std::f64::consts::PI / world.input.rows as f64;
            let area = world.input.radius_m.powi(2) * 2.0 * std::f64::consts::PI
                / world.input.cols as f64
                * ((la.to_radians() + dlat / 2.0).sin() - (la.to_radians() - dlat / 2.0).sin());
            let ratio = area.sqrt() / world.input.laws.channel_head_length_m;
            if ratio <= 1.0 {
                return Err(format!("planet_channel_head_exceeds_cell:{c}"));
            }
            let e = world.input.laws.stream_power_k_per_my * relief / ratio.ln();
            let sc = rocks.sc(world.input.lithology[c]);
            let width = half_width(relief, e, world.input.laws.diffusivity_m2_per_my, sc);
            if width == 0.0 {
                continue;
            }
            let seg = Segment {
                a,
                tangent,
                angle,
                half_width_m: width,
                erosion_m_per_my: e,
                critical_slope: sc,
                from: c as u32,
                to: r,
            };
            // Exact conservative bound for an arc: chord AABB plus sagitta,
            // expanded by its compact physical support. No neighbour-list
            // approximation or chart/window-specific selection.
            let pad = width + world.input.radius_m * (1.0 - (angle * 0.5).cos()) + 1.0;
            let low: [i32; 3] = std::array::from_fn(|i| {
                ((a[i].min(z[i]) * world.input.radius_m - pad) / BIN_M).floor() as i32
            });
            let high: [i32; 3] = std::array::from_fn(|i| {
                ((a[i].max(z[i]) * world.input.radius_m + pad) / BIN_M).floor() as i32
            });
            let id = world.segments.len() as u32;
            for x in low[0]..=high[0] {
                for y in low[1]..=high[1] {
                    for z in low[2]..=high[2] {
                        world.bins.entry([x, y, z]).or_default().push(id);
                    }
                }
            }
            world.axis[c] = Some(world.segments.len());
            world.segments.push(seg);
        }
        Ok(world)
    }
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }
    /// Active native pipeline path. The v2 constructor remains for historical comparison.
    pub fn new_joint(input: Inputs, root: &std::path::Path) -> R<Self> {
        let joint = crate::planet_catchment::World::open(root, &input)?;
        let mut world = Self::new(input)?;
        world.joint = Some(joint);
        Ok(world)
    }
    /// Load verified published endpoint data without private snapshot/trace layout.
    pub fn from_published(root: &std::path::Path) -> R<Self> {
        let (input, endpoint, hydro) = crate::planet_catchment::load_published(root)?;
        let joint = crate::planet_catchment::World::from_endpoint(endpoint, hydro, &input)?;
        let mut world = Self::new(input)?;
        world.joint = Some(joint);
        Ok(world)
    }
    pub fn process_block_work(&self, cell: u32) -> R<Vec<crate::planet_catchment::BlockWork>> {
        self.joint
            .as_ref()
            .ok_or("process blocks require the joint native model")?
            .block_work(cell)
    }
    pub fn process_block(
        &self,
        cell: u32,
        max_work_nodes: usize,
    ) -> R<crate::planet_catchment::ProcessBlock> {
        self.joint
            .as_ref()
            .ok_or("process blocks require the joint native model")?
            .process_block(cell, max_work_nodes)
    }
    pub fn process_provenance(&self) -> serde_json::Value {
        self.joint
            .as_ref()
            .map_or(serde_json::Value::Null, |w| w.provenance())
    }
    /// Distance only to this cell's published receiver axis. Its tributaries
    /// retain that actual downstream identity, not the nearest arbitrary basin.
    pub fn distance_to_cell_axis(&self, c: u32, lat: f64, lon: f64) -> f64 {
        let p = unit(lat, lon);
        self.axis[c as usize]
            .map(|id| {
                let s = &self.segments[id];
                let t = dot(p, s.tangent).atan2(dot(p, s.a)).clamp(0.0, s.angle);
                let chord = (0..3)
                    .map(|i| (p[i] - s.a[i] * t.cos() - s.tangent[i] * t.sin()).powi(2))
                    .sum::<f64>()
                    .sqrt();
                2.0 * self.input.radius_m * (chord * 0.5).min(1.0).asin()
            })
            .unwrap_or(f64::INFINITY)
    }
    fn cell_surface(&self, c: u32) -> R<std::sync::Arc<crate::planet_cell::CellSurface>> {
        if let Some(s) = self
            .cells
            .lock()
            .map_err(|_| "planet_cell_cache_poisoned")?
            .get(&c)
        {
            return Ok(s.clone());
        }
        let surface = std::sync::Arc::new(crate::planet_cell::CellSurface::build(self, c)?);
        let mut cells = self
            .cells
            .lock()
            .map_err(|_| "planet_cell_cache_poisoned")?;
        // Regenerable pure-function memo only. Durable published products and
        // receipts are never removed. Bounded for continuous planet travel.
        if cells.len() >= 128 {
            if let Some(k) = cells.keys().min().copied() {
                cells.remove(&k);
            }
        }
        Ok(cells.entry(c).or_insert(surface).clone())
    }
    pub fn cell_receipts(&self, cells: &std::collections::BTreeSet<u32>) -> R<serde_json::Value> {
        if let Some(w) = &self.joint {
            return w.receipts(cells);
        }
        let memo = self
            .cells
            .lock()
            .map_err(|_| "planet_cell_cache_poisoned")?;
        let rows = cells
            .iter()
            .filter_map(|c| memo.get(c).map(|s| &s.receipt))
            .collect::<Vec<_>>();
        serde_json::to_value(rows).map_err(|e| e.to_string())
    }
    pub fn sample(&self, lat: f64, lon: f64) -> R<Sample> {
        if let Some(w) = &self.joint {
            let height = w.sample(lat, lon)?;
            let datum = self.datum(lat, lon);
            return Ok(Sample {
                height_m: height,
                datum_m: datum,
                incision_m: datum - height,
                receiver_segment: None,
            });
        }
        let coarse = self.coarse_sample(lat, lon)?;
        let c = self.cell(lat, lon);
        let (clat, clon) = self.centre(c as usize);
        let da = 180.0 / self.input.rows as f64;
        let dl = 360.0 / self.input.cols as f64;
        let dlon = (lon - clon + 180.0).rem_euclid(360.0) - 180.0;
        let u = (dlon / dl + 0.5).clamp(0.0, 1.0);
        let v = ((lat - clat) / da + 0.5).clamp(0.0, 1.0);
        let height = if self.input.land[c as usize]
            && self.input.relief_m[c as usize] > 0.0
            && self.input.receiver[c as usize] != u32::MAX
        {
            self.cell_surface(c)?.sample(u, v)
        } else {
            // Unrefined ocean/terminal cells use the same global physical
            // lattice as refined cells, including shared edge vertices.
            let x = u * 256.0;
            let y = v * 256.0;
            let i = (x.floor() as usize).min(255);
            let j = (y.floor() as usize).min(255);
            let fx = x - i as f64;
            let fy = y - j as f64;
            let at = |i: usize, j: usize| -> R<f64> {
                let la = clat - da * 0.5 + j as f64 * da / 256.0;
                let lo =
                    (clon - dl * 0.5 + i as f64 * dl / 256.0 + 180.0).rem_euclid(360.0) - 180.0;
                Ok(self.coarse_sample(la, lo)?.height_m as f32 as f64)
            };
            (1.0 - fy) * ((1.0 - fx) * at(i, j)? + fx * at(i + 1, j)?)
                + fy * ((1.0 - fx) * at(i, j + 1)? + fx * at(i + 1, j + 1)?)
        };
        Ok(Sample {
            height_m: height,
            datum_m: coarse.datum_m,
            incision_m: coarse.datum_m - height,
            receiver_segment: coarse.receiver_segment,
        })
    }

    pub fn centre(&self, c: usize) -> (f64, f64) {
        (
            -90.0
                + (c / self.input.cols) as f64 * 180.0 / self.input.rows as f64
                + 90.0 / self.input.rows as f64,
            -180.0
                + (c % self.input.cols) as f64 * 360.0 / self.input.cols as f64
                + 180.0 / self.input.cols as f64,
        )
    }
    pub fn cell(&self, lat: f64, lon: f64) -> u32 {
        let row = (((lat + 90.0) / 180.0 * self.input.rows as f64).floor() as usize)
            .min(self.input.rows - 1);
        let col = (((lon + 180.0).rem_euclid(360.0) / 360.0 * self.input.cols as f64).floor()
            as usize)
            .min(self.input.cols - 1);
        (row * self.input.cols + col) as u32
    }
    /// Bilinear interpolation of source observations; caps taper the last row
    /// to one pole observation (its longitude mean), making longitude truly
    /// irrelevant at a pole. This is interpolation of ONE global datum, not
    /// blending of already solved chart surfaces.
    pub fn datum(&self, lat: f64, lon: f64) -> f64 {
        let x = (lon + 180.0).rem_euclid(360.0) / 360.0 * self.input.cols as f64 - 0.5;
        let j = x.floor() as i64;
        let fx = x - x.floor();
        let row = |i: usize| {
            let get = |j: i64| {
                self.input.elevation_m
                    [i * self.input.cols + j.rem_euclid(self.input.cols as i64) as usize]
                    as f64
            };
            get(j) * (1.0 - fx) + get(j + 1) * fx
        };
        let y = (lat + 90.0) / 180.0 * self.input.rows as f64 - 0.5;
        if y < 0.0 {
            return self.poles[0] * (-2.0 * y) + row(0) * (1.0 + 2.0 * y);
        }
        if y >= (self.input.rows - 1) as f64 {
            let f = 2.0 * (y - (self.input.rows - 1) as f64);
            return row(self.input.rows - 1) * (1.0 - f) + self.poles[1] * f;
        }
        let i = y.floor() as usize;
        let fy = y - y.floor();
        row(i) * (1.0 - fy) + row(i + 1) * fy
    }
    pub fn coarse_sample(&self, lat: f64, lon: f64) -> R<Sample> {
        if !lat.is_finite() || !lon.is_finite() || lat.abs() > 90.0 || lon.abs() > 180.0 {
            return Err("planet_position_invalid".into());
        }
        let p = unit(lat, lon);
        let datum = self.datum(lat, lon);
        let mut cut = 0.0;
        let mut source = None;
        if let Some(ids) = self.bins.get(&bin(p, self.input.radius_m)) {
            for &id in ids {
                let s = self.segments[id as usize];
                let t = dot(p, s.tangent).atan2(dot(p, s.a)).clamp(0.0, s.angle);
                let q: [f64; 3] =
                    std::array::from_fn(|i| s.a[i] * t.cos() + s.tangent[i] * t.sin());
                let chord = (0..3).map(|i| (p[i] - q[i]).powi(2)).sum::<f64>().sqrt();
                // Arc distance on the sphere, valid at poles and dateline.
                let distance = 2.0 * self.input.radius_m * (chord * 0.5).min(1.0).asin();
                let h = creep_height(
                    (s.half_width_m - distance).max(0.0),
                    s.erosion_m_per_my,
                    self.input.laws.diffusivity_m2_per_my,
                    s.critical_slope,
                );
                if h > cut {
                    cut = h;
                    source = Some([s.from, s.to]);
                }
            }
        }
        Ok(Sample {
            height_m: datum - cut,
            datum_m: datum,
            incision_m: cut,
            receiver_segment: source,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn creep_integral_obeys_flux_and_relief() {
        let (e, d, sc) = (100.0, 10000.0, 0.7);
        let width = half_width(250.0, e, d, sc);
        assert!((creep_height(width, e, d, sc) - 250.0).abs() < 1e-5);
        let x = width * 0.5;
        let dx = 0.01;
        let slope = (creep_height(x + dx, e, d, sc) - creep_height(x - dx, e, d, sc)) / (2.0 * dx);
        let flux = d * slope / (1.0 - (slope / sc).powi(2));
        assert!((flux - e * x).abs() / (e * x) < 1e-6);
    }
    #[test]
    fn field_wraps_and_pole_has_one_value() {
        let rows = 8;
        let cols = 16;
        let n = rows * cols;
        let w = WorldRelief::new(Inputs {
            seed: 42,
            rows,
            cols,
            radius_m: 6371000.0,
            elevation_m: (0..n).map(|i| i as f32).collect(),
            relief_m: vec![0.0; n],
            receiver: vec![u32::MAX; n],
            land: vec![true; n],
            lithology: vec![6; n],
            laws: Laws {
                stream_power_k_per_my: 1.0,
                diffusivity_m2_per_my: 10000.0,
                channel_head_length_m: 100.0,
            },
        })
        .unwrap();
        assert!(
            (w.sample(35.0, -180.0).unwrap().height_m - w.sample(35.0, 180.0).unwrap().height_m)
                .abs()
                < 1e-10
        );
        assert!(
            (w.sample(90.0, 0.0).unwrap().height_m - w.sample(90.0, 171.0).unwrap().height_m).abs()
                < 1e-10
        );
        assert!(w.sample(91.0, 0.0).is_err());
    }
}
