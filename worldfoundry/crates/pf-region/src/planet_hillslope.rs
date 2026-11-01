//! Conservative implicit finite-volume hillslope diffusion. D is the pinned
//! native m2/My law. This changes physical bed by transported bulk, not display
//! filtering. Wet/outer fixed beds exchange signed bulk through explicit ports.
use crate::R;
use pf_hydro_catchment::{Mesh, NONE};

const SIDES: [usize; 4] = [1, 3, 4, 6]; // south, west, east, north in native pipeline's D8 mesh
pub struct Operator {
    coefficients: Vec<[f64; 4]>,
    diagonal: Vec<f64>,
    free: Vec<bool>,
    area: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hillslope_moves_a_bump_without_creating_bulk() {
        let g = pf_hydro_catchment::terrain::Grid::new(3, 4, 1000.0).unwrap();
        let n = g.len();
        let mut mesh = Mesh {
            bed: vec![0.0; n],
            source: vec![0.0; n],
            point: Vec::new(),
            cell: (0..n).collect(),
            edges: vec![vec![NONE; 8]; n],
            fixed: vec![None; n],
            physical_len: n,
            radius: g.radius,
        };
        let mut area = Vec::new();
        for i in 0..n {
            mesh.point.push(g.point(i, 0.5, 0.5));
            area.push(g.area(i, 0.0, 1.0, 1.0));
            let mut slot = 0;
            for dy in -1..=1 {
                for dx in -1..=1 {
                    if dx != 0 || dy != 0 {
                        mesh.edges[i][slot] = g.neighbor(i, dx, dy).unwrap_or(NONE);
                        slot += 1;
                    }
                }
            }
        }
        mesh.bed[5] = 100.0;
        let before = 100.0 * area[5];
        let op = Operator::new(&mesh, &area, 3, 4, 100000.0, 0.5).unwrap();
        let budget = op.evolve(&mut mesh).unwrap();
        let after = mesh.bed.iter().zip(&area).map(|(z, a)| z * a).sum::<f64>();
        assert!((after - before).abs() / before < 1e-8);
        assert!(mesh.bed[5] < 100.0 && mesh.bed[4] > 0.0);
        assert!(mesh.bed.iter().all(|&z| z >= 0.0 && z <= 100.0));
        assert_eq!(budget.boundary_import_m3, 0.0);
        assert!(budget.redistributed_m3 > 0.0 && budget.relative_residual < 1e-4);
    }
}
pub struct Budget {
    pub boundary_import_m3: f64,
    pub redistributed_m3: f64,
    pub relative_residual: f64,
    pub iterations: usize,
}
impl Operator {
    pub fn new(mesh: &Mesh, area: &[f64], nlat: usize, nlon: usize, d: f64, dt: f64) -> R<Self> {
        let n = mesh.physical_len;
        if area.len() != n
            || nlat == 0
            || nlon == 0
            || !d.is_finite()
            || d <= 0.0
            || !dt.is_finite()
            || dt <= 0.0
        {
            return Err("hillslope_operator_inputs".into());
        }
        let dy = mesh.radius * std::f64::consts::PI / nlat as f64;
        let dl = 2.0 * std::f64::consts::PI / nlon as f64;
        let free: Vec<_> = (0..n).map(|i| mesh.fixed[i].is_none()).collect();
        let mut coefficients = vec![[0.0; 4]; n];
        let mut diagonal = area.to_vec();
        for i in 0..n {
            for (side, &slot) in SIDES.iter().enumerate() {
                let j = mesh.edges[i][slot];
                if j == NONE {
                    continue;
                }
                let face = if side == 1 || side == 2 {
                    dy
                } else {
                    mesh.radius * dl * ((mesh.point[i][1] + mesh.point[j][1]) * 0.5).cos()
                };
                let length = pf_hydro_catchment::terrain::distance(
                    mesh.point[i],
                    mesh.point[j],
                    mesh.radius,
                );
                let c = dt * d * face / length;
                if !c.is_finite() || c < 0.0 {
                    return Err("hillslope_metric".into());
                }
                coefficients[i][side] = c;
                diagonal[i] += c;
            }
        }
        Ok(Self {
            coefficients,
            diagonal,
            free,
            area: area.to_vec(),
        })
    }
    fn apply(&self, mesh: &Mesh, x: &[f64], out: &mut [f64]) {
        for i in 0..self.area.len() {
            out[i] = 0.0;
            if !self.free[i] {
                continue;
            }
            out[i] = self.diagonal[i] * x[i];
            for (side, &slot) in SIDES.iter().enumerate() {
                let j = mesh.edges[i][slot];
                if j < self.area.len() && self.free[j] {
                    out[i] -= self.coefficients[i][side] * x[j];
                }
            }
        }
    }
    pub fn evolve(&self, mesh: &mut Mesh) -> R<Budget> {
        let n = self.area.len();
        let mut delta = vec![0.0; n];
        let mut residual = vec![0.0; n];
        for i in 0..n {
            if self.free[i] {
                for (side, &slot) in SIDES.iter().enumerate() {
                    let j = mesh.edges[i][slot];
                    if j != NONE {
                        residual[i] += self.coefficients[i][side] * (mesh.bed[j] - mesh.bed[i]);
                    }
                }
            }
        }
        // Solve for height CHANGE to avoid subtracting two large absolute
        // elevation terms. Jacobi-preconditioned CG of the symmetric FV matrix.
        let mut z: Vec<_> = residual
            .iter()
            .zip(&self.diagonal)
            .map(|(r, a)| r / a)
            .collect();
        let mut p = z.clone();
        let mut ap = vec![0.0; n];
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
        let mut rz = dot(&residual, &z);
        let mut iterations = 0;
        let norm = |r: &[f64]| {
            r.iter()
                .zip(&self.diagonal)
                .map(|(r, d)| r.abs() / d)
                .fold(0.0, f64::max)
        };
        while norm(&residual) > 1e-8 {
            if iterations >= 128 {
                return Err("hillslope_implicit_solve_not_converged".into());
            }
            self.apply(mesh, &p, &mut ap);
            let den = dot(&p, &ap);
            if !den.is_finite() || den <= 0.0 {
                return Err("hillslope_matrix_not_positive".into());
            }
            let alpha = rz / den;
            for i in 0..n {
                delta[i] += alpha * p[i];
                residual[i] -= alpha * ap[i];
                z[i] = residual[i] / self.diagonal[i];
            }
            let next = dot(&residual, &z);
            let beta = next / rz;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
            rz = next;
            iterations += 1;
        }
        let old = mesh.bed.clone();
        for i in 0..n {
            mesh.bed[i] += delta[i];
        }
        let mut boundary = 0.0;
        let mut redistributed = 0.0;
        let mut scale = 0.0;
        let mut change = 0.0;
        for i in 0..n {
            if self.free[i] {
                let volume = self.area[i] * (mesh.bed[i] - old[i]);
                change += volume;
                scale += volume.abs();
                for (side, &slot) in SIDES.iter().enumerate() {
                    let j = mesh.edges[i][slot];
                    if j == NONE {
                        continue;
                    }
                    let flux = self.coefficients[i][side] * (mesh.bed[j] - mesh.bed[i]);
                    if j >= n || !self.free[j] {
                        boundary += flux;
                        scale += flux.abs();
                    } else if i < j {
                        redistributed += flux.abs();
                    }
                }
            }
        }
        let relative = (change - boundary).abs() / scale.max(1.0);
        if relative > 1e-4 {
            return Err(format!("hillslope_bulk_balance:{relative}"));
        }
        Ok(Budget {
            boundary_import_m3: boundary,
            redistributed_m3: redistributed,
            relative_residual: relative,
            iterations,
        })
    }
}
