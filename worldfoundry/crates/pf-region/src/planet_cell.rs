//! Native, world-addressed tributary construction. model design/native pipeline.
//!
//! Each T0 finite-volume cell has a fixed 256x256 subgrid, a first-passage
//! drainage tree to its published coarse receiver axis, stream-power channel
//! beds, and nonlinear-creep hillslopes (fast marching). The scalar E solve
//! targets the cell's coarse datum mean. Shared cell edges are Dirichlet data
//! from the SAME global physical field, never stitched chart outputs.
//!
//! Approximation: subgrid tributaries are confined to the source T0 cell and
//! discharge to its published coarse edge. The cell boundaries are fixed
//! coarse constraints, not asserted to be discovered watershed divides. No
//! meander evolution, sediment transport or canonical material transfer.
use crate::{planet_relief::WorldRelief, R};
use serde::Serialize;
use std::{cmp::Ordering, collections::BinaryHeap, time::Instant};
use vz_keyed_rng::Scope;
const SUB: usize = 256;
const N: usize = SUB + 1;
const D8: [(isize, isize); 8] = [
    (-1, -1),
    (0, -1),
    (1, -1),
    (-1, 0),
    (1, 0),
    (-1, 1),
    (0, 1),
    (1, 1),
];
#[derive(Clone, Copy)]
struct Key(f64, usize);
impl PartialEq for Key {
    fn eq(&self, b: &Self) -> bool {
        self.cmp(b) == Ordering::Equal
    }
}
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
impl Ord for Key {
    fn cmp(&self, b: &Self) -> Ordering {
        b.0.total_cmp(&self.0).then(b.1.cmp(&self.1))
    }
}
#[derive(Serialize)]
pub struct Receipt {
    pub cell: u32,
    pub nodes: usize,
    pub channel_nodes: usize,
    pub roots: usize,
    pub assigned_coarse_receiver: Option<u32>,
    pub cell_area_m2: f64,
    pub routed_area_m2: f64,
    pub erosion_m_per_my: f64,
    pub mean_target_m: f64,
    pub mean_surface_m: f64,
    pub mean_error_m: f64,
    pub mean_tolerance_m: f64,
    pub mean_converged: bool,
    pub evaluations: u32,
    pub seconds: f64,
    pub receiver_digest: String,
}
pub struct CellSurface {
    pub height: Vec<f32>,
    pub receipt: Receipt,
}
fn neighbour(i: usize, dx: isize, dy: isize) -> Option<usize> {
    let x = i % N;
    let y = i / N;
    let xx = x as isize + dx;
    let yy = y as isize + dy;
    if xx < 0 || yy < 0 || xx >= N as isize || yy >= N as isize {
        None
    } else {
        Some(yy as usize * N + xx as usize)
    }
}
fn boundary(i: usize) -> bool {
    i % N == 0 || i % N == SUB || i / N == 0 || i / N == SUB
}
fn fmm(z: &mut [f64], fixed: &[bool], dx: f64, dy: f64, slope: f64) {
    let mut done = fixed.to_vec();
    let mut heap = BinaryHeap::new();
    // Native soft-channel cap: proposed tributary beds are upper bounds,
    // not immovable pillars. A lower connected creep front may erode them.
    // Primary axes and shared physical boundary data remain Dirichlet.
    for i in 0..z.len() {
        if !fixed[i] && z[i].is_finite() {
            heap.push(Key(z[i], i));
        }
    }
    let update = |i: usize, z: &[f64], done: &[bool], heap: &mut BinaryHeap<Key>| {
        if done[i] {
            return;
        }
        let v = |x: isize, y: isize| {
            neighbour(i, x, y)
                .filter(|&k| done[k])
                .map(|k| z[k])
                .unwrap_or(f64::INFINITY)
        };
        let a = v(-1, 0).min(v(1, 0));
        let b = v(0, -1).min(v(0, 1));
        let mut value = (a + dx * slope).min(b + dy * slope);
        if a.is_finite() && b.is_finite() {
            let (wx, wy) = (1.0 / (dx * dx), 1.0 / (dy * dy));
            let aa = wx + wy;
            let bb = wx * a + wy * b;
            let cc = wx * a * a + wy * b * b - slope * slope;
            let disc = bb * bb - aa * cc;
            if disc >= 0.0 {
                let q = (bb + disc.sqrt()) / aa;
                if q >= a.max(b) {
                    value = value.min(q);
                }
            }
        }
        if value.is_finite() {
            heap.push(Key(value, i));
        }
    };
    for i in 0..z.len() {
        if done[i] {
            for (x, y) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                if let Some(j) = neighbour(i, x, y) {
                    update(j, z, &done, &mut heap);
                }
            }
        }
    }
    while let Some(Key(v, i)) = heap.pop() {
        if done[i] {
            continue;
        }
        done[i] = true;
        z[i] = v;
        for (x, y) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
            if let Some(j) = neighbour(i, x, y) {
                update(j, z, &done, &mut heap);
            }
        }
    }
}
impl CellSurface {
    pub fn build(w: &WorldRelief, c: u32) -> R<Self> {
        let start = Instant::now();
        let inp = &w.input;
        let (lat, lon) = w.centre(c as usize);
        let da = 180.0 / inp.rows as f64;
        let dl = 360.0 / inp.cols as f64;
        let radius = inp.radius_m;
        let dy = radius * da.to_radians() / SUB as f64;
        let dx = radius * lat.to_radians().cos() * dl.to_radians() / SUB as f64;
        if dx <= 0.0 {
            return Err("planet_polar_cell_metric".into());
        }
        let nn = N * N;
        let mut base = Vec::with_capacity(nn);
        let mut datum = Vec::with_capacity(nn);
        let mut seeds = vec![false; nn];
        let half_step = 0.6 * dx.hypot(dy);
        for j in 0..N {
            for i in 0..N {
                let la = lat - da * 0.5 + j as f64 * da / SUB as f64;
                let lo =
                    (lon - dl * 0.5 + i as f64 * dl / SUB as f64 + 180.0).rem_euclid(360.0) - 180.0;
                let s = w.coarse_sample(la, lo)?;
                base.push(s.height_m);
                datum.push(s.datum_m);
                seeds[j * N + i] = w.distance_to_cell_axis(c, la, lo) < half_step;
            }
        }
        // A declared terminal without a receiver has no generated tributary
        // solve. Caller only asks for cells with an actual receiver axis.
        if !seeds.iter().any(|v| *v) {
            return Err(format!("planet_cell_no_receiver_axis:{c}"));
        }
        let scope = Scope::root(inp.seed)
            .layer("planet-cell-route/1")
            .and_then(|s| s.address(&format!("cell/{c}")))
            .map_err(|e| e.to_string())?;
        let cost = (0..nn)
            .map(|i| {
                scope
                    .draw("t1b.rivernet.fpp.cost", i as u64, 0)
                    .map(|d| 0.5 + d.uniform_f64())
                    .map_err(|e| e.to_string())
            })
            .collect::<R<Vec<_>>>()?;
        let mut dist = vec![f64::INFINITY; nn];
        let mut recv = vec![usize::MAX; nn];
        let mut order = Vec::with_capacity(nn);
        let mut done = vec![false; nn];
        let mut heap = BinaryHeap::new();
        for i in 0..nn {
            if seeds[i] {
                dist[i] = 0.0;
                heap.push(Key(0.0, i));
            }
        }
        while let Some(Key(v, i)) = heap.pop() {
            if done[i] || v > dist[i] {
                continue;
            }
            done[i] = true;
            order.push(i);
            for (x, y) in D8 {
                if let Some(j) = neighbour(i, x, y) {
                    if done[j] {
                        continue;
                    }
                    let step = (x as f64 * dx).hypot(y as f64 * dy);
                    let q = v + step * (cost[i] + cost[j]) * 0.5;
                    if q < dist[j] {
                        dist[j] = q;
                        recv[j] = i;
                        heap.push(Key(q, j));
                    }
                }
            }
        }
        if order.len() != nn {
            return Err("planet_cell_unrouted".into());
        }
        let cell_area = radius.powi(2)
            * dl.to_radians()
            * ((lat + da * 0.5).to_radians().sin() - (lat - da * 0.5).to_radians().sin());
        // Trapezoidal node areas sum to exactly one spherical source cell.
        let mut area = vec![0.0; nn];
        let mut sum = 0.0;
        for j in 0..N {
            for i in 0..N {
                let wx = if i == 0 || i == SUB { 0.5 } else { 1.0 };
                let wy = if j == 0 || j == SUB { 0.5 } else { 1.0 };
                area[j * N + i] = wx * wy;
                sum += wx * wy;
            }
        }
        for a in &mut area {
            *a *= cell_area / sum;
        }
        let own_area = area.clone();
        for &i in order.iter().rev() {
            if recv[i] != usize::MAX {
                area[recv[i]] += area[i];
            }
        }
        let routed: f64 = (0..nn).filter(|&i| seeds[i]).map(|i| area[i]).sum();
        if (routed - cell_area).abs() / cell_area > 1e-8 {
            return Err("planet_cell_flow_accounting".into());
        }
        let rocks = crate::t1b_litho::LithologyV1::standin_v1();
        let sc = rocks.sc(inp.lithology[c as usize]);
        let km = rocks.k(inp.lithology[c as usize]);
        let channel = (0..nn)
            .map(|i| seeds[i] || area[i] >= 500_000.0)
            .collect::<Vec<_>>();
        let target = datum.iter().zip(&own_area).map(|(z, a)| z * a).sum::<f64>() / cell_area;
        let evaluate = |e: f64| {
            let mut bed = base.clone();
            for &i in &order {
                let r = recv[i];
                if r == usize::MAX {
                    continue;
                }
                let x = i % N;
                let y = i / N;
                let xr = r % N;
                let yr = r / N;
                let step = ((x as f64 - xr as f64) * dx).hypot((y as f64 - yr as f64) * dy);
                let slope = (e / (inp.laws.stream_power_k_per_my * km * area[i].sqrt().max(1.0)))
                    .min(0.6 * sc);
                bed[i] = bed[r] + slope * step;
            }
            let mut z = vec![f64::INFINITY; nn];
            let mut fixed = vec![false; nn];
            for i in 0..nn {
                if boundary(i) {
                    z[i] = base[i];
                    fixed[i] = true;
                } else if channel[i] {
                    z[i] = bed[i];
                    fixed[i] = seeds[i];
                }
            }
            let beta = 2.0 * e * 400.0 / (inp.laws.diffusivity_m2_per_my * sc);
            let slope = sc * beta / (1.0 + (1.0 + beta * beta).sqrt());
            fmm(&mut z, &fixed, dx, dy, slope);
            let mean = z.iter().zip(&own_area).map(|(z, a)| z * a).sum::<f64>() / cell_area;
            (z, mean)
        };
        let mut low = 0.0;
        let mut high = 20_000.0;
        let e0 = inp.laws.stream_power_k_per_my * inp.relief_m[c as usize] as f64
            / (cell_area.sqrt() / inp.laws.channel_head_length_m).ln();
        let mut e = e0.clamp(0.0, high);
        let mut best_error = f64::INFINITY;
        let mut best = Vec::new();
        let mut best_mean = 0.0;
        let mut best_e = e;
        let tolerance = 20.0;
        let mut evals = 0;
        for _ in 0..10 {
            let (z, mean) = evaluate(e);
            evals += 1;
            if !mean.is_finite() {
                return Err("planet_cell_nonfinite_surface".into());
            }
            if (mean - target).abs() < best_error {
                best_error = (mean - target).abs();
                best_mean = mean;
                best = z;
                best_e = e;
            }
            if best_error <= tolerance {
                break;
            }
            if mean > target {
                high = e;
            } else {
                low = e;
            }
            e = if low == 0.0 {
                high * 0.5
            } else {
                (low * high).sqrt()
            };
        }
        let receiver_bytes: Vec<u8> = recv
            .iter()
            .flat_map(|&v| (if v == usize::MAX { u32::MAX } else { v as u32 }).to_le_bytes())
            .collect();
        Ok(Self {
            height: best.into_iter().map(|v| v as f32).collect(),
            receipt: Receipt {
                cell: c,
                nodes: nn,
                channel_nodes: channel.iter().filter(|v| **v).count(),
                roots: seeds.iter().filter(|v| **v).count(),
                assigned_coarse_receiver: Some(inp.receiver[c as usize]),
                cell_area_m2: cell_area,
                routed_area_m2: routed,
                erosion_m_per_my: best_e,
                mean_target_m: target,
                mean_surface_m: best_mean,
                mean_error_m: best_mean - target,
                mean_tolerance_m: tolerance,
                mean_converged: best_error <= tolerance,
                evaluations: evals,
                seconds: start.elapsed().as_secs_f64(),
                receiver_digest: format!("blake3:{}", blake3::hash(&receiver_bytes).to_hex()),
            },
        })
    }
    pub fn sample(&self, u: f64, v: f64) -> f64 {
        let x = (u * SUB as f64).clamp(0.0, SUB as f64);
        let y = (v * SUB as f64).clamp(0.0, SUB as f64);
        let i = (x.floor() as usize).min(SUB - 1);
        let j = (y.floor() as usize).min(SUB - 1);
        let fx = x - i as f64;
        let fy = y - j as f64;
        let get = |i, j| self.height[j * N + i] as f64;
        (1.0 - fy) * ((1.0 - fx) * get(i, j) + fx * get(i + 1, j))
            + fy * ((1.0 - fx) * get(i, j + 1) + fx * get(i + 1, j + 1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn anisotropic_march_matches_a_physical_plane() {
        let mut z = vec![f64::INFINITY; N * N];
        let mut fixed = vec![false; N * N];
        for j in 0..N {
            z[j * N] = 1000.0;
            fixed[j * N] = true;
        }
        fmm(&mut z, &fixed, 20.0, 70.0, 0.3);
        for (i, j) in [(1, 0), (100, 128), (SUB, SUB)] {
            assert!((z[j * N + i] - (1000.0 + 6.0 * i as f64)).abs() < 1e-5);
        }
    }
}
