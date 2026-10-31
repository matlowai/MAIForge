//! MC-Q2 diagnosis: WHY the coupled Stokes solve floors on an `l = 1`
//! right-hand side, decided by a DENSE DIRECT solve of the same system.
//!
//! **This module is a diagnosis, not a fix.** Nothing in [`super::stokes`],
//! [`super::stokes_solve`], [`super::stokes_saddle`],
//! [`super::qualify_thieulot`], [`super::benchmark_thieulot`] or `pf-grid` is
//! modified; this module reads their public surface only, and `mod.rs` gained
//! one line. Every operator it measures is applied through
//! [`SaddleSystem`]'s own methods, so a dense column here is by construction
//! the same column the matrix-free solver applies - and that identity is
//! CHECKED against random vectors every run
//! ([`DenseSaddle::matrix_free_disagreement`]) rather than assumed.
//!
//! # The question this module was built to answer
//!
//! `MCQ1_THIEULOT_FINDINGS.md` section 6 established, by measurement, that
//! the FGMRES + block-preconditioner solve reaches `1e-8` on `l = 0` and
//! `l = 2` forcings and floors at `1.3e-4` on `l = 1`, that the floor lives
//! in the continuity row, and that it is not the stall detector, the deferred
//! correction, the inner CG tolerance or the radial profile. It did not say
//! why. Two hypotheses survive and they have opposite fixes:
//!
//! - **discretization**: the discrete system with an `l = 1` right-hand side
//!   has no solution at all - `b` is not in the range of the saddle matrix -
//!   in which case any solver floors and the fix is structural;
//! - **solver**: the discrete system is perfectly solvable and the FGMRES
//!   iteration, or the preconditioner, cannot reach the solution.
//!
//! A dense direct solve separates them and nothing else does. At `n = 2` the
//! system carries 72 velocity unknowns (120 faces, 48 of them on the two
//! walls and pinned) and 48 pressure unknowns; at `n = 4`, 576 and 288. Both
//! are small enough to assemble column by column and factor.
//!
//! # The third hypothesis, and why the least-squares form is the right tool
//!
//! There is a third possibility that a square LU solve alone cannot see, and
//! it is the one this module found. The block preconditioner re-projects the
//! `A`-block solve's output onto `null(K)^perp_{M_u}` on every call
//! (`PreconditionerOptions::reproject_inner_in_m_u`, default true, at
//! `worldfoundry/crates/pf-geodynamics/src/mantle/stokes_saddle.rs:1063`).
//! Right-preconditioned FGMRES builds its iterate as `x = sum_j y_j z_j` out
//! of preconditioner OUTPUTS, so the velocity iterate lives, exactly, in the
//! `M_u`-orthogonal complement of the six `A`-null modes - which is what the
//! `solved mode content: tra 3.9e-13 rot 1.7e-13` line of the MC-Q1 run
//! reports. The achievable residual is therefore not
//! `min over all x of |M x - b|` but
//!
//! ```text
//! min over { (u, p) : <n_j, u>_{M_u} = 0 for the six modes } of |M (u,p) - b|
//! ```
//!
//! and those two minima are different numbers. MC-Q1 section 6.3 tested the
//! six-mode hypothesis by evaluating the residual AT ONE POINT of that
//! subspace (the exact free-slip field minus its six-mode part, `4.184e-3`)
//! and refuted it on magnitude. That point is not the minimiser, and a
//! subspace of dimension `nu - 6` has a great deal of room to do better than
//! any single member of it, so the refutation does not follow. This module
//! computes the constrained MINIMUM exactly.
//!
//! # What is deliberately NOT done here
//!
//! No tolerance is moved, no operator is changed, and no run is reported as
//! converged that is not. Where a rank has to be decided from a spectrum, the
//! residual is printed as a FUNCTION of the cut so the cut is visible as a
//! choice.

use super::stokes_saddle::{GradientMode, ModeBasis, SaddleSystem};
// Only the test drivers name a deflation set; cfg-gated so the library build
// carries no unused import (the `qualify_thieulot` pattern).
#[cfg(test)]
use super::stokes_saddle::Deflation;
use super::MantleError;

// ---------------------------------------------------------------------
// Dense linear algebra, column major. Written out because this crate has no
// linear-algebra dependency and a diagnosis packet must not add one.

/// A column-major dense matrix. Column major because every kernel below
/// (Householder QR, one-sided Jacobi) works on whole columns.
#[derive(Clone, Debug)]
pub struct Dense {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Dense {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.set(i, i, 1.0);
        }
        m
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    #[inline]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[j * self.rows + i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[j * self.rows + i] = v;
    }

    pub fn col(&self, j: usize) -> &[f64] {
        &self.data[j * self.rows..(j + 1) * self.rows]
    }

    pub fn col_mut(&mut self, j: usize) -> &mut [f64] {
        &mut self.data[j * self.rows..(j + 1) * self.rows]
    }

    /// `A x`.
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.cols, "dimension mismatch in Dense::mul_vec");
        let mut out = vec![0.0f64; self.rows];
        for j in 0..self.cols {
            let xj = x[j];
            if xj == 0.0 {
                continue;
            }
            for (slot, a) in out.iter_mut().zip(self.col(j).iter()) {
                *slot += a * xj;
            }
        }
        out
    }

    /// `A^T y`.
    pub fn mul_transpose_vec(&self, y: &[f64]) -> Vec<f64> {
        assert_eq!(
            y.len(),
            self.rows,
            "dimension mismatch in Dense::mul_transpose_vec"
        );
        (0..self.cols)
            .map(|j| self.col(j).iter().zip(y.iter()).map(|(a, b)| a * b).sum())
            .collect()
    }

    /// Frobenius norm of the block `rows in [r0,r1) x cols in [c0,c1)`.
    pub fn block_norm(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> f64 {
        let mut acc = 0.0;
        for j in c0..c1 {
            for i in r0..r1 {
                let v = self.at(i, j);
                acc += v * v;
            }
        }
        acc.sqrt()
    }
}

fn nrm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Householder QR with column pivoting, with `b` carried through the
/// reflectors. The whole point is `residual_at_rank`: for a rank-deficient
/// least-squares problem the minimum residual is `|(Q^T b)[r..]|`, and
/// printing it against `r` shows what the rank cut is worth instead of hiding
/// the cut inside one number.
#[derive(Debug, Clone)]
pub struct PivotedQr {
    /// Upper triangular factor, in place (only the upper triangle is used).
    pub r: Dense,
    /// `Q^T b`.
    pub qtb: Vec<f64>,
    /// Column permutation: `perm[k]` is the original column in position `k`.
    pub perm: Vec<usize>,
    pub m: usize,
    pub n: usize,
}

impl PivotedQr {
    /// `|R_kk|` down the diagonal, in pivot order (non-increasing in exact
    /// arithmetic).
    pub fn rdiag(&self) -> Vec<f64> {
        (0..self.n.min(self.m)).map(|k| self.r.at(k, k).abs()).collect()
    }

    /// The least-squares residual `min |A x - b|` if the numerical rank is
    /// taken to be `rank`.
    pub fn residual_at_rank(&self, rank: usize) -> f64 {
        nrm2(&self.qtb[rank.min(self.m)..])
    }

    /// Numerical rank at a relative cut on `|R_kk|`.
    pub fn rank_at(&self, rel_cut: f64) -> usize {
        let d = self.rdiag();
        let top = d.first().copied().unwrap_or(0.0);
        d.iter().take_while(|v| **v > rel_cut * top).count()
    }

    /// The basic least-squares solution at a stated rank: the free columns
    /// are set to zero. Not the minimum-norm solution; where the null space
    /// matters the SVD path below is used instead and the two are compared.
    pub fn solve(&self, rank: usize) -> Vec<f64> {
        let r = rank.min(self.n).min(self.m);
        let mut y = vec![0.0f64; self.n];
        for k in (0..r).rev() {
            let mut acc = self.qtb[k];
            for j in (k + 1)..r {
                acc -= self.r.at(k, j) * y[j];
            }
            let d = self.r.at(k, k);
            y[k] = if d == 0.0 { 0.0 } else { acc / d };
        }
        let mut x = vec![0.0f64; self.n];
        for k in 0..self.n {
            x[self.perm[k]] = y[k];
        }
        x
    }
}

/// Householder QR with column pivoting of `a`, applying the reflectors to
/// `b` as it goes.
pub fn pivoted_qr(a: &Dense, b: &[f64]) -> PivotedQr {
    let (m, n) = (a.rows(), a.cols());
    assert_eq!(b.len(), m, "right-hand side length must match the row count");
    let mut r = a.clone();
    let mut qtb = b.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut colnorm: Vec<f64> = (0..n).map(|j| nrm2(r.col(j))).collect();
    let steps = m.min(n);
    for k in 0..steps {
        // Pivot on the largest remaining column norm, recomputed rather than
        // downdated: downdating loses accuracy exactly where the rank
        // decision is being made, and this matrix is small.
        for j in k..n {
            let mut acc = 0.0;
            for i in k..m {
                let v = r.at(i, j);
                acc += v * v;
            }
            colnorm[j] = acc.sqrt();
        }
        let mut best = k;
        for j in (k + 1)..n {
            if colnorm[j] > colnorm[best] {
                best = j;
            }
        }
        if best != k {
            perm.swap(k, best);
            colnorm.swap(k, best);
            for i in 0..m {
                let t = r.at(i, k);
                let u = r.at(i, best);
                r.set(i, k, u);
                r.set(i, best, t);
            }
        }
        // Householder reflector for column k, rows k..m.
        let mut x: Vec<f64> = (k..m).map(|i| r.at(i, k)).collect();
        let alpha = nrm2(&x);
        if alpha == 0.0 {
            continue;
        }
        let sign = if x[0] >= 0.0 { 1.0 } else { -1.0 };
        x[0] += sign * alpha;
        let vnorm = nrm2(&x);
        if vnorm == 0.0 {
            continue;
        }
        for v in x.iter_mut() {
            *v /= vnorm;
        }
        for j in k..n {
            let mut s = 0.0;
            for (t, i) in (k..m).enumerate() {
                s += x[t] * r.at(i, j);
            }
            s *= 2.0;
            for (t, i) in (k..m).enumerate() {
                let v = r.at(i, j) - s * x[t];
                r.set(i, j, v);
            }
        }
        let mut s = 0.0;
        for (t, i) in (k..m).enumerate() {
            s += x[t] * qtb[i];
        }
        s *= 2.0;
        for (t, i) in (k..m).enumerate() {
            qtb[i] -= s * x[t];
        }
    }
    PivotedQr {
        r,
        qtb,
        perm,
        m,
        n,
    }
}

/// A singular value decomposition by ONE-SIDED Jacobi.
///
/// One-sided Jacobi and not `A^T A`: squaring the matrix squares the
/// condition number, and the small singular values are the entire subject of
/// this module. One-sided Jacobi computes them to high RELATIVE accuracy,
/// which is what a claim of "rank deficient beyond the pressure constant"
/// needs.
#[derive(Debug, Clone)]
pub struct Svd {
    /// Descending.
    pub sigma: Vec<f64>,
    /// Left singular vectors as columns, `m x n`.
    pub u: Dense,
    /// Right singular vectors as columns, `n x n`.
    pub v: Dense,
    pub sweeps: usize,
    pub converged: bool,
}

/// One-sided Jacobi SVD. `max_sweeps` is a budget and the report says whether
/// it was reached; a run that exhausts it is reported as unconverged, never
/// silently used.
pub fn jacobi_svd(a: &Dense, max_sweeps: usize) -> Svd {
    let (m, n) = (a.rows(), a.cols());
    let mut w = a.clone();
    let mut v = Dense::identity(n);
    let eps = 1e-15;
    let mut sweeps = 0;
    let mut converged = false;
    while sweeps < max_sweeps {
        sweeps += 1;
        let mut rotated = 0usize;
        for p in 0..n.saturating_sub(1) {
            for q in (p + 1)..n {
                let alpha = dot(w.col(p), w.col(p));
                let beta = dot(w.col(q), w.col(q));
                let gamma = dot(w.col(p), w.col(q));
                if gamma == 0.0 || gamma.abs() <= eps * (alpha * beta).sqrt() {
                    continue;
                }
                rotated += 1;
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = if zeta >= 0.0 {
                    1.0 / (zeta + (1.0 + zeta * zeta).sqrt())
                } else {
                    -1.0 / (-zeta + (1.0 + zeta * zeta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for i in 0..m {
                    let wp = w.at(i, p);
                    let wq = w.at(i, q);
                    w.set(i, p, c * wp - s * wq);
                    w.set(i, q, s * wp + c * wq);
                }
                for i in 0..n {
                    let vp = v.at(i, p);
                    let vq = v.at(i, q);
                    v.set(i, p, c * vp - s * vq);
                    v.set(i, q, s * vp + c * vq);
                }
            }
        }
        if rotated == 0 {
            converged = true;
            break;
        }
    }
    // Column norms are the singular values; normalise to get U.
    let mut order: Vec<usize> = (0..n).collect();
    let norms: Vec<f64> = (0..n).map(|j| nrm2(w.col(j))).collect();
    order.sort_by(|a, b| norms[*b].partial_cmp(&norms[*a]).unwrap());
    let mut sigma = vec![0.0f64; n];
    let mut u = Dense::zeros(m, n);
    let mut vs = Dense::zeros(n, n);
    for (k, &j) in order.iter().enumerate() {
        sigma[k] = norms[j];
        let inv = if norms[j] > 0.0 { 1.0 / norms[j] } else { 0.0 };
        for i in 0..m {
            u.set(i, k, w.at(i, j) * inv);
        }
        for i in 0..n {
            vs.set(i, k, v.at(i, j));
        }
    }
    Svd {
        sigma,
        u,
        v: vs,
        sweeps,
        converged,
    }
}

/// Dense LU with partial pivoting. Returns the solution and the smallest
/// absolute pivot, which is the cheapest available rank signal.
pub fn lu_solve(a: &Dense, b: &[f64]) -> Option<(Vec<f64>, f64)> {
    let n = a.rows();
    assert_eq!(a.cols(), n, "LU needs a square matrix");
    assert_eq!(b.len(), n);
    let mut lu = a.clone();
    let mut x = b.to_vec();
    let mut smallest = f64::MAX;
    for k in 0..n {
        let mut piv = k;
        let mut best = lu.at(k, k).abs();
        for i in (k + 1)..n {
            let v = lu.at(i, k).abs();
            if v > best {
                best = v;
                piv = i;
            }
        }
        smallest = smallest.min(best);
        if best == 0.0 {
            return None;
        }
        if piv != k {
            for j in 0..n {
                let t = lu.at(k, j);
                let u = lu.at(piv, j);
                lu.set(k, j, u);
                lu.set(piv, j, t);
            }
            x.swap(k, piv);
        }
        let d = lu.at(k, k);
        for i in (k + 1)..n {
            let f = lu.at(i, k) / d;
            if f == 0.0 {
                continue;
            }
            lu.set(i, k, f);
            for j in (k + 1)..n {
                let v = lu.at(i, j) - f * lu.at(k, j);
                lu.set(i, j, v);
            }
            x[i] -= f * x[k];
        }
    }
    for k in (0..n).rev() {
        let mut acc = x[k];
        for j in (k + 1)..n {
            acc -= lu.at(k, j) * x[j];
        }
        x[k] = acc / lu.at(k, k);
    }
    Some((x, smallest))
}

// ---------------------------------------------------------------------
// The dense saddle matrix

/// The saddle-point system of [`SaddleSystem`], assembled densely and
/// SYMMETRICALLY WEIGHTED so that Euclidean norms in this module are the
/// norms the solver actually measures.
///
/// Row weights: `sqrt(M_u_f)` on a momentum row, `theta sqrt(V_c)` on a
/// continuity row, so `|R r|_2` is exactly
/// [`SaddleSystem::residual_norm`]. Column weights are the same numbers on
/// the matching unknowns, so `|C x|_2` is
/// `sqrt(|u|^2_{M_u} + theta^2 |p|^2_{M_p})`. The scaled matrix is
/// `Mhat = R M C^-1` and the scaled system is `Mhat y = R b`, `y = C x`.
///
/// Boundary faces carry no unknown and no equation: `M_u` is zero there,
/// [`SaddleSystem::viscous`] and both gradient halves return zero there, and
/// the residual norm cannot see them. They are simply absent from the dense
/// index map, which is why `nu` is the interior face count and not
/// `face_count`.
#[derive(Debug, Clone)]
pub struct DenseSaddle {
    /// Face index of each velocity unknown.
    pub iface: Vec<usize>,
    pub nu: usize,
    pub np: usize,
    /// `(nu + np) x (nu + np)`, weighted.
    pub m: Dense,
    pub row_scale: Vec<f64>,
    pub col_scale: Vec<f64>,
    pub grad: GradientMode,
}

impl DenseSaddle {
    pub fn dim(&self) -> usize {
        self.nu + self.np
    }

    /// Assemble column by column through [`SaddleSystem`]'s own operator
    /// methods. `sys.apply` is not used for the columns because it also
    /// removes the pressure mean from the continuity row, which is a
    /// projection of the OUTPUT and would make the dense matrix a different
    /// operator from the block form; the mean removal is exactly a projection
    /// onto the complement of the `M_p`-constant, and `D u` is already
    /// `M_p`-mean-free by telescoping, so the two agree to roundoff. That
    /// agreement is what [`Self::matrix_free_disagreement`] measures.
    pub fn assemble(sys: &SaddleSystem, grad: GradientMode) -> Result<Self, MantleError> {
        let topo = sys.topology();
        let iface: Vec<usize> = topo
            .faces()
            .iter()
            .enumerate()
            .filter(|(_, f)| f.neighbor.is_some())
            .map(|(i, _)| i)
            .collect();
        let nu = iface.len();
        let np = sys.cell_count();
        let n = nu + np;
        let theta = sys.theta();
        let mut row_scale = vec![0.0f64; n];
        let mut col_scale = vec![0.0f64; n];
        for (j, &f) in iface.iter().enumerate() {
            let s = sys.m_u()[f].sqrt();
            row_scale[j] = s;
            col_scale[j] = s;
        }
        for c in 0..np {
            let s = theta * sys.m_p()[c].sqrt();
            row_scale[nu + c] = s;
            col_scale[nu + c] = s;
        }
        let mut m = Dense::zeros(n, n);
        // Velocity columns: (B e_j, D e_j).
        let mut e = sys.zeros_u();
        for (j, &f) in iface.iter().enumerate() {
            e[f] = 1.0;
            let bu = sys.viscous(&e)?;
            let du = sys.divergence(&e)?;
            e[f] = 0.0;
            let inv = 1.0 / col_scale[j];
            for (i, &fi) in iface.iter().enumerate() {
                let v = row_scale[i] * bu[fi] * inv;
                if v != 0.0 {
                    m.set(i, j, v);
                }
            }
            for c in 0..np {
                let v = row_scale[nu + c] * du[c] * inv;
                if v != 0.0 {
                    m.set(nu + c, j, v);
                }
            }
        }
        // Pressure columns: (G e_c, 0).
        let mut q = sys.zeros_p();
        for c in 0..np {
            q[c] = 1.0;
            let g = sys.gradient(&q, grad)?;
            q[c] = 0.0;
            let inv = 1.0 / col_scale[nu + c];
            for (i, &fi) in iface.iter().enumerate() {
                let v = row_scale[i] * g[fi] * inv;
                if v != 0.0 {
                    m.set(i, nu + c, v);
                }
            }
        }
        Ok(Self {
            iface,
            nu,
            np,
            m,
            row_scale,
            col_scale,
            grad,
        })
    }

    /// Scale a full-face-length velocity and a cell pressure into the dense
    /// unknown vector `y = C x`.
    pub fn pack(&self, u: &[f64], p: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0f64; self.dim()];
        for (j, &f) in self.iface.iter().enumerate() {
            y[j] = self.col_scale[j] * u[f];
        }
        for c in 0..self.np {
            y[self.nu + c] = self.col_scale[self.nu + c] * p[c];
        }
        y
    }

    /// Inverse of [`Self::pack`]: back to a full-face-length velocity (zero
    /// on the walls) and a cell pressure.
    pub fn unpack(&self, y: &[f64], face_count: usize) -> (Vec<f64>, Vec<f64>) {
        let mut u = vec![0.0f64; face_count];
        for (j, &f) in self.iface.iter().enumerate() {
            u[f] = y[j] / self.col_scale[j];
        }
        let mut p = vec![0.0f64; self.np];
        for c in 0..self.np {
            p[c] = y[self.nu + c] / self.col_scale[self.nu + c];
        }
        (u, p)
    }

    /// Scale a right-hand side: `R b`.
    pub fn pack_rhs(&self, b_u: &[f64], b_p: &[f64]) -> Vec<f64> {
        let mut r = vec![0.0f64; self.dim()];
        for (i, &f) in self.iface.iter().enumerate() {
            r[i] = self.row_scale[i] * b_u[f];
        }
        for c in 0..self.np {
            r[self.nu + c] = self.row_scale[self.nu + c] * b_p[c];
        }
        r
    }

    /// GUARD AGAINST THE INSTRUMENT. Worst relative disagreement between
    /// `Mhat (C x)` and `R (matrix-free apply)(x)` over `trials`
    /// deterministic pseudo-random `(u, p)` pairs. If the dense assembly were
    /// wrong - a transposed index, a missed boundary face, a scale applied
    /// twice - every number this module prints would be about a different
    /// matrix, and this is the only check that says otherwise.
    pub fn matrix_free_disagreement(
        &self,
        sys: &SaddleSystem,
        modes: &ModeBasis,
        trials: usize,
        seed: u64,
    ) -> Result<f64, MantleError> {
        let mut state = seed | 1;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            ((state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
        };
        let mut worst = 0.0f64;
        for _ in 0..trials {
            let mut u = sys.zeros_u();
            for &f in self.iface.iter() {
                u[f] = next();
            }
            let mut p: Vec<f64> = (0..self.np).map(|_| next()).collect();
            sys.remove_mean_pressure(&mut p)?;
            let (ru, rp) = sys.apply(modes, self.grad, &u, &p)?;
            let want = self.pack_rhs(&ru, &rp);
            let got = self.m.mul_vec(&self.pack(&u, &p));
            let scale = nrm2(&want).max(1e-300);
            let diff: f64 = nrm2(
                &want
                    .iter()
                    .zip(got.iter())
                    .map(|(a, b)| a - b)
                    .collect::<Vec<f64>>(),
            );
            worst = worst.max(diff / scale);
        }
        Ok(worst)
    }

    /// The six (or three) `A`-null velocity modes packed into the dense
    /// SCALED coordinates, Euclidean-orthonormal there.
    ///
    /// They already are orthonormal: the modes are `M_u`-orthonormal and the
    /// velocity column scale is `sqrt(M_u)`, so `q_i . q_j = <n_i, n_j>_{M_u}`.
    /// The re-orthonormalisation below therefore costs nothing and makes the
    /// projector below exact to roundoff rather than to the accuracy of the
    /// `M_u` Gram-Schmidt that built the modes.
    pub fn packed_modes(&self, modes: &ModeBasis) -> Vec<Vec<f64>> {
        let n = self.dim();
        let mut qs: Vec<Vec<f64>> = Vec::new();
        for j in 0..modes.dim() {
            let mode = modes.mode(j);
            let mut q = vec![0.0f64; n];
            for (i, &f) in self.iface.iter().enumerate() {
                q[i] = self.col_scale[i] * mode[f];
            }
            for prev in qs.iter() {
                let c = dot(prev, &q);
                for (slot, x) in q.iter_mut().zip(prev.iter()) {
                    *slot -= c * x;
                }
            }
            let nq = nrm2(&q);
            if nq > 0.0 {
                for slot in q.iter_mut() {
                    *slot /= nq;
                }
                qs.push(q);
            }
        }
        qs
    }

    /// `Mhat` with the velocity block restricted to the `M_u`-orthogonal
    /// complement of `modes`: `Mhat P` with `P = I - sum_j q_j q_j^T`.
    ///
    /// The minimum of `|Mhat P z - bhat|` over all `z` IS the minimum of
    /// `|Mhat y - bhat|` over the constrained set, which is the residual
    /// floor a solver whose iterate lives in that set can reach.
    pub fn restricted(&self, modes: &ModeBasis) -> (Dense, Vec<Vec<f64>>) {
        let qs = self.packed_modes(modes);
        (restrict_columns(&self.m, &qs), qs)
    }
}

/// `A (I - sum_j q_j q_j^T)`: the operator seen by a solver whose iterate is
/// confined to the complement of `qs`.
pub fn restrict_columns(a: &Dense, qs: &[Vec<f64>]) -> Dense {
    let n = a.cols();
    let mut out = a.clone();
    for q in qs.iter() {
        let aq = out.mul_vec(q);
        for j in 0..n {
            let qj = q[j];
            if qj == 0.0 {
                continue;
            }
            for i in 0..a.rows() {
                let v = out.at(i, j) - aq[i] * qj;
                out.set(i, j, v);
            }
        }
    }
    out
}

/// `(I - sum_j q_j q_j^T) A`: the system a `Deflation` variant actually
/// solves, with the mode-tested momentum rows DELETED. `fgmres_saddle`
/// applies exactly this projection to both the operator output and the
/// right-hand side (`ModeBasis::project` on `b_u` and inside
/// `SaddleSystem::apply`), so deleting rows here and projecting `b` there is
/// the same discrete system - a DIFFERENT one from the undeflated system, as
/// `Deflation::changes_the_system` says.
pub fn project_rows(a: &Dense, qs: &[Vec<f64>]) -> Dense {
    let mut out = a.clone();
    for q in qs.iter() {
        // row_j <- row_j - q_j (q_j^T A)
        let qta = out.mul_transpose_vec(q);
        for j in 0..out.cols() {
            let c = qta[j];
            if c == 0.0 {
                continue;
            }
            for i in 0..out.rows() {
                let v = out.at(i, j) - q[i] * c;
                out.set(i, j, v);
            }
        }
    }
    out
}

/// `(I - sum_j q_j q_j^T) v`.
pub fn project_vec(v: &[f64], qs: &[Vec<f64>]) -> Vec<f64> {
    let mut out = v.to_vec();
    for q in qs.iter() {
        let c = dot(q, &out);
        for (slot, x) in out.iter_mut().zip(q.iter()) {
            *slot -= c * x;
        }
    }
    out
}


// ---------------------------------------------------------------------
// The augmented A block, `K_gamma = K + gamma D^T M_p D`
//
// Added by the follow-up packet (sections 9-10 of the findings). Nothing in
// `stokes_saddle.rs` is touched: the augmentation is expressed through
// `SaddleSystem`'s existing operator methods, which is exactly what makes it
// screenable without writing a solver.

/// `B_gamma u = B u + gamma M_u^-1 D^T M_p D u`, matrix free.
///
/// [`SaddleSystem::gradient_implicit`] is documented as
/// `G_imp p = -M_u^-1 D^T M_p p`, EXACTLY the negative weighted adjoint of
/// [`SaddleSystem::divergence`], so the augmentation term is
/// `-gamma G_imp(D u)` and needs no new assembler.
///
/// **The implicit half is the right one and the corrected one is not.** The
/// augmentation has to be the exact adjoint pair `D^T M_p D` or it is not
/// symmetric positive semidefinite, and the deferred non-orthogonality /
/// skewness correction is not part of that adjoint. The adjoint identity
/// `<G_imp p, u>_{M_u} = -<p, D u>_{M_p}` is MEASURED every run by
/// [`adjoint_disagreement`] rather than trusted to the doc comment.
///
/// `<B_gamma u, u>_{M_u} = <B u, u>_{M_u} + gamma |D u|^2_{M_p} >= 0`, so the
/// augmentation can only raise a Rayleigh quotient, never lower one. On a
/// mode with `K n = 0` it raises it by exactly `gamma |D n|^2_{M_p}`, which
/// is the whole question this screening asks.
pub fn augmented_a_apply(
    sys: &SaddleSystem,
    gamma: f64,
    u: &[f64],
) -> Result<Vec<f64>, MantleError> {
    let mut out = sys.viscous(u)?;
    if gamma != 0.0 {
        let d = sys.divergence(u)?;
        let g = sys.gradient_implicit(&d)?;
        for (slot, x) in out.iter_mut().zip(g.iter()) {
            *slot -= gamma * x;
        }
    }
    Ok(out)
}

/// GUARD AGAINST THE INSTRUMENT. Worst relative violation of
/// `<G_imp p, u>_{M_u} + <p, D u>_{M_p} = 0` over `trials` deterministic
/// pseudo-random pairs.
///
/// If this is not roundoff then the "augmentation" assembled above is not
/// `D^T M_p D`, it is not symmetric, and every eigenvalue below is about a
/// different operator.
pub fn adjoint_disagreement(
    sys: &SaddleSystem,
    iface: &[usize],
    trials: usize,
    seed: u64,
) -> Result<f64, MantleError> {
    let mut state = seed | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        ((state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
    };
    let mut worst = 0.0f64;
    for _ in 0..trials {
        let mut u = sys.zeros_u();
        for &f in iface.iter() {
            u[f] = next();
        }
        let p: Vec<f64> = (0..sys.cell_count()).map(|_| next()).collect();
        let g = sys.gradient_implicit(&p)?;
        let d = sys.divergence(&u)?;
        let lhs = sys.dot_u(&g, &u);
        let rhs = -sys.dot_p(&p, &d);
        let scale = lhs.abs().max(rhs.abs()).max(1e-300);
        worst = worst.max((lhs - rhs).abs() / scale);
    }
    Ok(worst)
}

/// The `A` block and its augmentation term, assembled densely in the SAME
/// symmetric `sqrt(M_u)` scaling [`DenseSaddle`] uses, so a singular value
/// here is comparable with the saddle spectrum of section 1.2 and a Rayleigh
/// quotient here is one in `M_u`.
///
/// `b` is `M_u^-1/2 K M_u^-1/2` and `c` is `M_u^-1/2 D^T M_p D M_u^-1/2`.
/// Both are symmetric by construction, and both symmetries are MEASURED
/// ([`DenseABlock::symmetry_defect`]) because "symmetric by construction" is
/// a claim about the code that assembles them.
#[derive(Debug, Clone)]
pub struct DenseABlock {
    pub iface: Vec<usize>,
    pub nu: usize,
    pub b: Dense,
    pub c: Dense,
    /// `sqrt(M_u_f)` per velocity unknown.
    pub scale: Vec<f64>,
}

impl DenseABlock {
    pub fn assemble(sys: &SaddleSystem) -> Result<Self, MantleError> {
        let topo = sys.topology();
        let iface: Vec<usize> = topo
            .faces()
            .iter()
            .enumerate()
            .filter(|(_, f)| f.neighbor.is_some())
            .map(|(i, _)| i)
            .collect();
        let nu = iface.len();
        let scale: Vec<f64> = iface.iter().map(|&f| sys.m_u()[f].sqrt()).collect();
        let mut b = Dense::zeros(nu, nu);
        let mut c = Dense::zeros(nu, nu);
        let mut e = sys.zeros_u();
        for (j, &f) in iface.iter().enumerate() {
            e[f] = 1.0;
            let bu = sys.viscous(&e)?;
            let d = sys.divergence(&e)?;
            let g = sys.gradient_implicit(&d)?;
            e[f] = 0.0;
            let inv = 1.0 / scale[j];
            for (i, &fi) in iface.iter().enumerate() {
                let vb = scale[i] * bu[fi] * inv;
                if vb != 0.0 {
                    b.set(i, j, vb);
                }
                // `+ D^T M_p D` is `-G_imp(D .)`.
                let vc = -scale[i] * g[fi] * inv;
                if vc != 0.0 {
                    c.set(i, j, vc);
                }
            }
        }
        Ok(Self {
            iface,
            nu,
            b,
            c,
            scale,
        })
    }

    /// `Bhat + gamma Chat`.
    pub fn combine(&self, gamma: f64) -> Dense {
        let mut out = self.b.clone();
        if gamma != 0.0 {
            for j in 0..self.nu {
                for i in 0..self.nu {
                    let v = out.at(i, j) + gamma * self.c.at(i, j);
                    out.set(i, j, v);
                }
            }
        }
        out
    }

    /// `|X - X^T|_F / |X|_F` for the two blocks.
    pub fn symmetry_defect(&self) -> (f64, f64) {
        let one = |m: &Dense| {
            let mut num = 0.0f64;
            let mut den = 0.0f64;
            for j in 0..m.cols() {
                for i in 0..m.rows() {
                    let a = m.at(i, j);
                    let t = m.at(j, i);
                    num += (a - t) * (a - t);
                    den += a * a;
                }
            }
            if den > 0.0 {
                (num / den).sqrt()
            } else {
                0.0
            }
        };
        (one(&self.b), one(&self.c))
    }

    /// Pack a full-face-length velocity into the scaled dense coordinates.
    pub fn pack(&self, u: &[f64]) -> Vec<f64> {
        self.iface
            .iter()
            .enumerate()
            .map(|(j, &f)| self.scale[j] * u[f])
            .collect()
    }
}

// ---------------------------------------------------------------------
// Exactly solenoidal poloidal fields of a stated angular degree

/// `P_l(x)` for `l <= 3`.
pub fn legendre_p(l: usize, x: f64) -> f64 {
    match l {
        0 => 1.0,
        1 => x,
        2 => 0.5 * (3.0 * x * x - 1.0),
        3 => 0.5 * (5.0 * x * x * x - 3.0 * x),
        4 => (35.0 * x.powi(4) - 30.0 * x * x + 3.0) / 8.0,
        5 => (63.0 * x.powi(5) - 70.0 * x.powi(3) + 15.0 * x) / 8.0,
        _ => panic!("legendre_p is written for l <= 5 only; l = {l} was asked for"),
    }
}

/// `dP_l/dx` for `l <= 3`.
pub fn legendre_dp(l: usize, x: f64) -> f64 {
    match l {
        0 => 0.0,
        1 => 1.0,
        2 => 3.0 * x,
        3 => 0.5 * (15.0 * x * x - 3.0),
        4 => (140.0 * x.powi(3) - 60.0 * x) / 8.0,
        5 => (315.0 * x.powi(4) - 210.0 * x * x + 15.0) / 8.0,
        _ => panic!("legendre_dp is written for l <= 5 only; l = {l} was asked for"),
    }
}

/// An axisymmetric poloidal velocity of angular degree `l` that is EXACTLY
/// divergence free in the continuum and impermeable on both walls.
///
/// Derivation, in one place. For `v_r = g(r) P_l(cos th)` and
/// `v_th = F(r) dP_l(cos th)/dth`, spherical incompressibility is
///
/// ```text
/// (1/r^2) d(r^2 g)/dr P_l + (F/(r sin th)) d(sin th dP_l/dth)/dth = 0
/// ```
///
/// and the angular factor is `-l(l+1) P_l`, so `div v = 0` for every `r`
/// exactly when `F = (r^2 g)' / (l(l+1) r)`. `g` is chosen to vanish at both
/// walls, so `v . n = 0` there. The construction is CHECKED, not trusted:
/// [`continuum_divergence`] differences the Cartesian field and every test
/// prints the result before using the field for anything.
///
/// `l = 0` is refused, and the refusal is the physics: `l = 0` forces
/// `(r^2 g)' = 0`, i.e. `g = C / r^2`, a net radial flux, which impermeable
/// walls forbid unless `C = 0`. **There is no nontrivial solenoidal `l = 0`
/// field on a closed shell**, so the `l = 0` row of any comparison like this
/// is necessarily the zero field.
pub fn solenoidal_poloidal(l: usize, r1: f64, r2: f64) -> impl Fn([f64; 3]) -> [f64; 3] {
    assert!(
        l >= 1,
        "l = 0 has no nontrivial solenoidal poloidal field on a closed shell"
    );
    let k = std::f64::consts::PI / (r2 - r1);
    let ll = (l * (l + 1)) as f64;
    move |x: [f64; 3]| {
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        let s = (x[0] * x[0] + x[1] * x[1]).sqrt();
        let (ct, st) = (x[2] / r, s / r);
        let (cp, sp) = (x[0] / s, x[1] / s);
        let g = (k * (r - r1)).sin();
        let gp = k * (k * (r - r1)).cos();
        let f = (2.0 * g + r * gp) / ll;
        let vr = g * legendre_p(l, ct);
        // d/dth P_l(cos th) = -sin th P_l'(cos th).
        let vt = f * (-st * legendre_dp(l, ct));
        [
            vr * st * cp + vt * ct * cp,
            vr * st * sp + vt * ct * sp,
            vr * ct - vt * st,
        ]
    }
}

/// A deliberately NON-solenoidal control of degree `l`: the same radial
/// profile put entirely in `v_r` with no tangential return flow. Its
/// continuum divergence is `O(1)`, and it exists so every table that reports
/// a small `|D u|` also reports what a large one looks like.
pub fn radial_only(l: usize, r1: f64, r2: f64) -> impl Fn([f64; 3]) -> [f64; 3] {
    let k = std::f64::consts::PI / (r2 - r1);
    move |x: [f64; 3]| {
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        let ct = x[2] / r;
        let g = (k * (r - r1)).sin() * legendre_p(l, ct);
        [g * x[0] / r, g * x[1] / r, g * x[2] / r]
    }
}

/// Central-difference `div v` at a point, step `h`. The instrument check for
/// [`solenoidal_poloidal`]: a hand-derived tangential profile is exactly the
/// kind of thing that is wrong by a factor of `l(l+1)` and still looks
/// plausible.
pub fn continuum_divergence(v: &impl Fn([f64; 3]) -> [f64; 3], x: [f64; 3], h: f64) -> f64 {
    let mut acc = 0.0;
    for d in 0..3 {
        let mut xp = x;
        let mut xm = x;
        xp[d] += h;
        xm[d] -= h;
        acc += (v(xp)[d] - v(xm)[d]) / (2.0 * h);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::super::qualify_thieulot::{
        area_mean_face_field, body_force_faces, face_labels, pair_slopes, reference_topology,
        slope, FreeSlipCase,
    };
    use super::super::benchmark_thieulot::ThieulotCase;
    use super::super::stokes::ViscousOperator;
    use super::super::stokes_saddle::{
        divergence_leak, fgmres_saddle, BlockPreconditioner, CoarseSpace, FacetClass,
        FgmresOptions, NullCoarse, PreconditionerOptions,
    };
    use super::super::stokes_solve::sampled_face_field;
    use super::*;
    use pf_grid::shell_ops::{FaceGradientRule, ShellTopology};

    /// The rungs this module runs. `(2,2)` always; `(4,3)` under
    /// `L1_DIAG_N4=1`, because the dense work there is `864^3` and the box is
    /// memory bound. Whichever rungs ran is printed with every table.
    fn rungs() -> Vec<(usize, usize)> {
        let mut v = vec![(2usize, 2usize)];
        if std::env::var("L1_DIAG_N4").is_ok() {
            v.push((4, 3));
        }
        v
    }

    /// Body force `rho(r) P_l(cos th) e_r` with Thieulot's radial profile,
    /// point-sampled on face normals exactly as
    /// [`body_force_faces`] does for `l = 1`.
    fn force_of_degree(topo: &ShellTopology, case: &ThieulotCase, l: usize) -> Vec<f64> {
        body_force_faces(topo, move |r, th| {
            // `case.density(r, th)` is `mu0 * shape(r) * cos(th)`, so dividing
            // by cos(th) and multiplying by P_l gives the same radial profile
            // at a different degree. Done through the public surface rather
            // than by copying `density_shape`, which is private.
            let c = th.cos();
            let shape = case.density(r, 0.0);
            shape * legendre_p(l, c)
        })
    }

    struct Rung<'a> {
        sys: SaddleSystem<'a>,
        none: ModeBasis,
        six: ModeBasis,
        tra: ModeBasis,
        rot: ModeBasis,
    }

    fn build_rung<'a>(topo: &'a ShellTopology, op: &'a ViscousOperator, mu0: f64) -> Rung<'a> {
        let sys = SaddleSystem::build(topo, op, mu0, FaceGradientRule::NonOrthogonalAndSkew)
            .expect("saddle system builds on the nondimensional shell");
        let none = ModeBasis::build(topo, &sys, Deflation::None).expect("empty basis");
        let six = ModeBasis::build(topo, &sys, Deflation::Six).expect("six modes");
        let tra = ModeBasis::build(topo, &sys, Deflation::Translations).expect("translations");
        let rot = ModeBasis::build(topo, &sys, Deflation::Rotations).expect("rotations");
        Rung {
            sys,
            none,
            six,
            tra,
            rot,
        }
    }

    /// Minimum residual of `|A y - b|` reported both ways: the pivoted-QR
    /// rank ladder and, at the small rung, the SVD. Printing the ladder is
    /// the honesty requirement - the "rank" of a saddle matrix with an
    /// `O(h^-2)` spread of singular values is a CHOICE, and the residual as a
    /// function of that choice shows how much the choice is worth.
    fn lsq_ladder(tag: &str, a: &Dense, b: &[f64]) -> (f64, usize) {
        let qr = pivoted_qr(a, b);
        let bn = nrm2(b);
        let d = qr.rdiag();
        let top = d.first().copied().unwrap_or(0.0);
        let full = a.cols().min(a.rows());
        let mut chosen = (0.0f64, 0usize);
        let mut line = String::new();
        for cut in [1e-8f64, 1e-10, 1e-12, 1e-14] {
            let rank = d.iter().take_while(|v| **v > cut * top).count();
            let res = qr.residual_at_rank(rank) / bn;
            line.push_str(&format!(" cut {cut:.0e}: rank {rank:>4} res {res:.3e} |"));
            if cut == 1e-12 {
                chosen = (res, rank);
            }
        }
        println!(
            "  {tag:<34} n {full:>4} |R00| {top:.3e} |Rnn| {:.3e} |{line}",
            d.last().copied().unwrap_or(0.0)
        );
        chosen
    }

    /// THE DECISIVE EXPERIMENT.
    ///
    /// Solve the same `l = 1` saddle-point system densely and directly, three
    /// ways: unconstrained least squares (is the system solvable at all?),
    /// least squares restricted to the subspace the preconditioner's output
    /// actually spans (is that restriction the floor?), and a gauge-fixed LU
    /// (an independent factorisation, so the answer does not rest on one
    /// kernel). The `l = 2` and `l = 3` forcings are run identically as the
    /// control, because a number that is small only for `l = 1` proves
    /// nothing without the degrees that converge.
    #[test]
    fn the_dense_direct_solve_separates_the_discretization_from_the_solver() {
        let case = ThieulotCase::isoviscous_reference();
        for (n, layers) in rungs() {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            let dense = DenseSaddle::assemble(&r.sys, GradientMode::Corrected)
                .expect("dense assembly of the saddle matrix");
            let worst = dense
                .matrix_free_disagreement(&r.sys, &r.none, 4, 0x5eed_1234)
                .expect("matrix-free cross-check");
            println!(
                "\n=== n {n} layers {layers}: faces {} interior {} cells {} dense {}x{} \
                 theta {:.6e} h_ref {:.6e}",
                topo.face_count(),
                dense.nu,
                dense.np,
                dense.dim(),
                dense.dim(),
                r.sys.theta(),
                r.sys.mesh_scale(),
            );
            println!(
                "  INSTRUMENT CHECK  worst |Mhat C x - R apply(x)| / |R apply(x)| over 4 random \
                 pairs = {worst:.3e}   (this must be roundoff or nothing below is about the \
                 solver's matrix)"
            );
            assert!(
                worst < 1e-10,
                "the dense assembly is not the matrix-free operator: {worst:.3e}"
            );
            println!(
                "  block Frobenius norms: B {:.4e}  G {:.4e}  D {:.4e}",
                dense.m.block_norm(0, dense.nu, 0, dense.nu),
                dense.m.block_norm(0, dense.nu, dense.nu, dense.dim()),
                dense.m.block_norm(dense.nu, dense.dim(), 0, dense.nu),
            );

            for l in [1usize, 2, 3, 4, 5] {
                let f = force_of_degree(&topo, &case, l);
                let zp = r.sys.zeros_p();
                let bhat = dense.pack_rhs(&f, &zp);
                let bn = nrm2(&bhat);
                println!("\n  --- forcing degree l = {l},  |b| = {bn:.6e}");
                let (res_free, rank_free) = lsq_ladder("UNCONSTRAINED least squares", &dense.m, &bhat);

                // The subspace the preconditioner's velocity output spans.
                for (label, modes) in [
                    ("restricted: u perp SIX      ", &r.six),
                    ("restricted: u perp trans(3) ", &r.tra),
                    ("restricted: u perp rot(3)   ", &r.rot),
                ] {
                    let (a, _) = dense.restricted(modes);
                    lsq_ladder(label, &a, &bhat);
                }

                // Independent factorisation: gauge-fix the square system by
                // replacing the LAST continuity row (the continuity rows have
                // an exact left null vector, the M_p-constant, so one of them
                // is redundant) with the mean-free pressure gauge.
                let mut sq = dense.m.clone();
                let mut rhs = bhat.clone();
                let last = dense.dim() - 1;
                for j in 0..dense.dim() {
                    sq.set(last, j, 0.0);
                }
                for c in 0..dense.np {
                    // sum_c V_c p_c = 0 in scaled coordinates: the coefficient
                    // on y_{nu+c} is V_c / col_scale = V_c / (theta sqrt(V_c)).
                    let v = r.sys.m_p()[c] / dense.col_scale[dense.nu + c];
                    sq.set(last, dense.nu + c, v);
                }
                rhs[last] = 0.0;
                match lu_solve(&sq, &rhs) {
                    None => println!("  LU gauge-fixed: EXACTLY SINGULAR (zero pivot)"),
                    Some((y, smallest)) => {
                        let res = dense.m.mul_vec(&y);
                        let diff: Vec<f64> =
                            res.iter().zip(bhat.iter()).map(|(a, b)| a - b).collect();
                        println!(
                            "  LU gauge-fixed: smallest |pivot| {smallest:.3e}  |M y - b|/|b| \
                             {:.3e}  |y| {:.3e}",
                            nrm2(&diff) / bn,
                            nrm2(&y)
                        );
                        // Where does the LU answer sit relative to the six modes?
                        let (u, _p) = dense.unpack(&y, topo.face_count());
                        let comps = r.six.components(&r.sys, &u);
                        let tra = comps[0..3].iter().fold(0.0f64, |a, c| a.max(c.abs()));
                        let rot = comps[3..6].iter().fold(0.0f64, |a, c| a.max(c.abs()));
                        println!(
                            "                  direct solution mode content: tra {tra:.3e} \
                             rot {rot:.3e}  (an FGMRES iterate is 1e-13 in both by construction)"
                        );
                    }
                }
                let _ = (res_free, rank_free);
            }
        }
    }

    /// The singular spectrum of the `l = 1` saddle matrix, and what the
    /// right-hand side has in its near-null LEFT singular directions.
    ///
    /// If the discrete system were inconsistent - `b` outside `range(M)` -
    /// the floor would be exactly `|U_null^T b| / |b|`, and this test would
    /// show it. `n = 2` always; `n = 4` costs `864^3` per sweep and is opt-in.
    #[test]
    fn the_singular_spectrum_of_the_l1_saddle_matrix() {
        let case = ThieulotCase::isoviscous_reference();
        let (n, layers) = (2usize, 2usize);
        let topo = reference_topology(&case, n, layers).expect("topology");
        let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
        let r = build_rung(&topo, &op, case.mu0());
        let dense = DenseSaddle::assemble(&r.sys, GradientMode::Corrected).expect("dense");
        let svd = jacobi_svd(&dense.m, 60);
        let nn = svd.sigma.len();
        println!(
            "\n=== SVD n {n} layers {layers}: {}x{} sweeps {} converged {}",
            dense.dim(),
            dense.dim(),
            svd.sweeps,
            svd.converged
        );
        assert!(svd.converged, "one-sided Jacobi did not converge in 60 sweeps");
        println!(
            "  sigma_max {:.6e}  sigma_min {:.6e}  ratio {:.3e}",
            svd.sigma[0],
            svd.sigma[nn - 1],
            svd.sigma[0] / svd.sigma[nn - 1].max(1e-300)
        );
        let tail: Vec<String> = svd.sigma[nn - 10..].iter().map(|s| format!("{s:.3e}")).collect();
        println!("  ten smallest singular values: {tail:?}");
        // WHAT the small directions are. For each of the ten smallest, the
        // fraction of the right singular vector's energy that sits in the
        // pressure block, and its overlap with the six `A`-null velocity
        // modes. This is the packet's "are they l = 1, and what do they look
        // like in velocity vs pressure" question, answered directly.
        let q_six = dense.packed_modes(&r.six);
        let q_tra = dense.packed_modes(&r.tra);
        let q_rot = dense.packed_modes(&r.rot);
        println!(
            "  character of the ten smallest right singular vectors \
             (pressure fraction, then overlap with the mode sets):"
        );
        for k in (nn - 10)..nn {
            let v = svd.v.col(k);
            let pe: f64 = v[dense.nu..].iter().map(|x| x * x).sum();
            let ov = |qs: &Vec<Vec<f64>>| -> f64 {
                qs.iter().map(|q| dot(q, v).powi(2)).sum::<f64>().sqrt()
            };
            println!(
                "    sigma {:.4e}  pressure energy {:.4}  |proj tra| {:.4}  |proj rot| {:.4} \
                 |proj six| {:.4}",
                svd.sigma[k],
                pe,
                ov(&q_tra),
                ov(&q_rot),
                ov(&q_six)
            );
        }
        // The Gram spectrum the designed remedy actually inverts. `NullCoarse`
        // fits the six-mode amplitude against the CONTINUITY row alone, so its
        // Gram matrix is `<D n_i, D n_j>_{M_p}`; the rotations carry
        // `|D r| = O(h^2)`, so those rows are `O(h^4)` and the 6x6 solve is
        // the ill-conditioned object. Printed rather than argued.
        for which in [CoarseSpace::Translations, CoarseSpace::Six] {
            let nc = NullCoarse::build(&topo, &r.sys, which).expect("coarse space");
            let ev: Vec<String> = nc.eigenvalues().iter().map(|e| format!("{e:.3e}")).collect();
            let lo = nc.eigenvalues().iter().cloned().fold(f64::MAX, f64::min);
            let hi = nc.eigenvalues().iter().cloned().fold(0.0f64, f64::max);
            println!(
                "  NullCoarse {:>15}: Gram <D n_i, D n_j>_Mp eigenvalues {ev:?}  spread {:.3e}",
                which.label(),
                hi / lo.max(1e-300)
            );
        }

        for l in [1usize, 2, 3] {
            let f = force_of_degree(&topo, &case, l);
            let zp = r.sys.zeros_p();
            let bhat = dense.pack_rhs(&f, &zp);
            let bn = nrm2(&bhat);
            // |U_j^T b| / |b| for the ten smallest, and the residual that
            // remains if every direction with sigma below a cut is
            // unreachable.
            let mut over: Vec<String> = Vec::new();
            for k in (nn - 10)..nn {
                let c = dot(svd.u.col(k), &bhat).abs() / bn;
                over.push(format!("{c:.3e}"));
            }
            println!("  l = {l}: |u_k^T b|/|b| for the ten smallest sigma: {over:?}");
            for cut in [1e-8f64, 1e-10, 1e-12, 1e-14] {
                let thresh = cut * svd.sigma[0];
                let mut acc = 0.0;
                for k in 0..nn {
                    if svd.sigma[k] <= thresh {
                        let c = dot(svd.u.col(k), &bhat);
                        acc += c * c;
                    }
                }
                // Everything outside span(U) as well: |b|^2 - sum_k (u_k^T b)^2.
                let mut inside = 0.0;
                for k in 0..nn {
                    let c = dot(svd.u.col(k), &bhat);
                    inside += c * c;
                }
                let outside = (bn * bn - inside).max(0.0);
                println!(
                    "        cut {cut:.0e}: unreachable-from-small-sigma {:.3e}  \
                     outside-span(U) {:.3e}",
                    acc.sqrt() / bn,
                    outside.sqrt() / bn
                );
            }
        }
    }

    /// The continuity row, examined directly. `|D u|` for fields that are
    /// EXACTLY divergence free in the continuum, by angular degree, against a
    /// deliberately non-solenoidal control at the same amplitude.
    ///
    /// The construction of each field is checked by finite differences before
    /// it is used, and the check is printed.
    #[test]
    fn the_discrete_divergence_of_continuum_solenoidal_fields_by_degree() {
        let case = ThieulotCase::isoviscous_reference();
        let (r1, r2) = (case.r_inner(), case.r_outer());
        for (n, layers) in rungs() {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            println!("\n=== |D u| by angular degree, n {n} layers {layers}");
            println!(
                "  l = 0 is ABSENT and that is the physics: (r^2 g)' = 0 forces g = C/r^2, a net \
                 radial flux, which impermeable walls forbid. The only solenoidal l = 0 poloidal \
                 field on a closed shell is zero."
            );
            for l in [1usize, 2, 3] {
                let v = solenoidal_poloidal(l, r1, r2);
                // Instrument check: is this field actually solenoidal?
                let mut worst_div = 0.0f64;
                let mut scale = 0.0f64;
                for c in topo.cell_centroids_m().expect("centroids").iter() {
                    let d = continuum_divergence(&v, *c, 1e-5);
                    let s = {
                        let f = v(*c);
                        (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt()
                    };
                    worst_div = worst_div.max(d.abs());
                    scale = scale.max(s);
                }
                let (dofs, geom) = area_mean_face_field(&topo, &v).expect("area-mean sampling");
                assert!(geom < 1e-12, "facet vertex cycle mismatch {geom:.3e}");
                assert!(
                    dofs.iter().all(|x| x.is_finite()),
                    "a facet sampled a non-finite field: a quadrature point is on the polar axis"
                );
                let du = r.sys.divergence(&dofs).expect("divergence");
                let un = r.sys.norm_u(&dofs);
                let dn = r.sys.norm_p(&du);
                // The non-solenoidal control at the same degree.
                let w = radial_only(l, r1, r2);
                let (wd, _) = area_mean_face_field(&topo, &w).expect("area-mean sampling");
                let wdu = r.sys.divergence(&wd).expect("divergence");
                let wun = r.sys.norm_u(&wd);
                let wdn = r.sys.norm_p(&wdu);
                println!(
                    "  l = {l}: continuum |div v| <= {worst_div:.3e} (|v| <= {scale:.3e}) | \
                     solenoidal h|D u|_Mp/|u|_Mu {:.4e} | radial-only control {:.4e}",
                    r.sys.mesh_scale() * dn / un,
                    r.sys.mesh_scale() * wdn / wun,
                );
            }
            // And the exact free-slip solution itself, which is the l = 1
            // field the driven solve is actually chasing. MC-Q1 section 6.3
            // measured 3.348e-3 for |D u|_Mp / |f| at n = 4; this is the same
            // quantity recomputed independently here.
            let free = FreeSlipCase::new(case).expect("free-slip reference");
            let (dofs, geom) =
                area_mean_face_field(&topo, |x| free.velocity_cartesian(x)).expect("sampling");
            assert!(geom < 1e-12);
            let du = r.sys.divergence(&dofs).expect("divergence");
            let f = body_force_faces(&topo, |rad, th| case.density(rad, th));
            let fnorm = r.sys.residual_norm(&f, &r.sys.zeros_p());
            println!(
                "  exact free-slip (l = 1): h|D u|_Mp/|u|_Mu {:.4e}   |D u|_Mp/|f| {:.4e}",
                r.sys.mesh_scale() * r.sys.norm_p(&du) / r.sys.norm_u(&dofs),
                r.sys.norm_p(&du) / fnorm,
            );
        }
    }

    /// CONTRADICTION HUNT, kept because it is the check that decides whether
    /// section 2's whole argument holds.
    ///
    /// The claim is that the FGMRES velocity iterate lies in
    /// `null(K)^perp_{M_u}` and therefore cannot beat the constrained dense
    /// minimum. A `restart = 200` run at `n = 4` reached `1.0700e-4`, BELOW
    /// the computed constrained minimum `1.313e-4`. Either the claim is
    /// wrong, or the iterate is not in the subspace. This test measures the
    /// iterate's six-mode content directly, and re-evaluates the residual
    /// after projecting that content out - if the projected residual jumps
    /// back to the constrained minimum, the content is real and load bearing.
    #[test]
    fn where_the_restart_200_iterate_actually_lives() {
        // n = 4 work: a 1056^3 dense QR plus three 2000-outer FGMRES runs.
        // Opt-in with the other n = 4 work so the default debug suite stays
        // cheap, exactly as `qualify_thieulot` gates its ladder.
        if std::env::var("L1_DIAG_N4").is_err() {
            println!("  skipped: set L1_DIAG_N4=1 to run the n = 4 restart contradiction hunt");
            return;
        }
        let case = ThieulotCase::isoviscous_reference();
        let (n, layers) = (4usize, 3usize);
        let topo = reference_topology(&case, n, layers).expect("topology");
        let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
        let r = build_rung(&topo, &op, case.mu0());
        let dense = DenseSaddle::assemble(&r.sys, GradientMode::Corrected).expect("dense");
        let f = force_of_degree(&topo, &case, 1);
        let zp = r.sys.zeros_p();
        let bhat = dense.pack_rhs(&f, &zp);
        let bn = nrm2(&bhat);
        for restart in [60usize, 200, 400] {
            let mut pre = BlockPreconditioner::build(
                &topo,
                &op,
                &r.sys,
                PreconditionerOptions {
                    coarse: CoarseSpace::None,
                    ..PreconditionerOptions::default()
                },
            )
            .expect("preconditioner");
            let (u, p, rep) = fgmres_saddle(
                &r.sys,
                &mut pre,
                &r.none,
                &r.six,
                GradientMode::Corrected,
                &f,
                &r.sys.zeros_p(),
                FgmresOptions {
                    max_iterations: 2000,
                    stall_window: 2000,
                    restart,
                    ..FgmresOptions::default()
                },
            )
            .expect("FGMRES runs");
            let comps = r.six.components(&r.sys, &u);
            let tra = comps[0..3].iter().fold(0.0f64, |a, c| a.max(c.abs()));
            let rot = comps[3..6].iter().fold(0.0f64, |a, c| a.max(c.abs()));
            // Residual through the DENSE matrix, as an independent evaluation
            // of the same quantity the report gives.
            let y = dense.pack(&u, &p);
            let res = dense.m.mul_vec(&y);
            let d: Vec<f64> = res.iter().zip(bhat.iter()).map(|(a, b)| a - b).collect();
            // ... and after removing the six-mode part of the velocity.
            let mut u2 = u.clone();
            r.six.project(&r.sys, &mut u2);
            let y2 = dense.pack(&u2, &p);
            let res2 = dense.m.mul_vec(&y2);
            let d2: Vec<f64> = res2.iter().zip(bhat.iter()).map(|(a, b)| a - b).collect();
            println!(
                "  restart {restart:>3}: outer {:>4} cycles {:>3} true_rel(report) {:.4e} \
                 true_rel(dense) {:.4e}  after removing six-mode part {:.4e}  \
                 mode content tra {tra:.3e} rot {rot:.3e}  |u|_Mu {:.4e}",
                rep.outer_iterations,
                rep.restart_cycles,
                rep.true_rel_residual,
                nrm2(&d) / bn,
                nrm2(&d2) / bn,
                r.sys.norm_u(&u),
            );
        }
        // And the constrained minimum this is being compared against,
        // recomputed in the same run so the two numbers cannot drift apart.
        let (a, _) = dense.restricted(&r.six);
        let qr = pivoted_qr(&a, &bhat);
        let d = qr.rdiag();
        let top = d.first().copied().unwrap_or(0.0);
        let rank = d.iter().take_while(|v| **v > 1e-12 * top).count();
        println!(
            "  constrained dense minimum over u perp six: rank {rank} res {:.4e}",
            qr.residual_at_rank(rank) / bn
        );
    }

    /// The same question at `n = 8`, where no dense matrix is affordable, so
    /// the residual is evaluated through `SaddleSystem::apply` instead. The
    /// point is one number: does the `restart = 200` iterate that reaches
    /// `1.46e-4` live INSIDE the preconditioner's subspace (in which case
    /// `restart = 60`'s `6.98e-2` was restart stagnation and nothing else) or
    /// outside it (in which case the comparison proves less)?
    #[test]
    fn the_n8_restart_gain_is_inside_or_outside_the_subspace() {
        if std::env::var("L1_DIAG_N8").is_err() {
            println!("  skipped: set L1_DIAG_N8=1 to run the n = 8 rung");
            return;
        }
        let case = ThieulotCase::isoviscous_reference();
        let (n, layers) = (8usize, 6usize);
        let topo = reference_topology(&case, n, layers).expect("topology");
        let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
        let r = build_rung(&topo, &op, case.mu0());
        let f = force_of_degree(&topo, &case, 1);
        let b_norm = r.sys.residual_norm(&f, &r.sys.zeros_p());
        for restart in [60usize, 200] {
            let mut pre = BlockPreconditioner::build(
                &topo,
                &op,
                &r.sys,
                PreconditionerOptions {
                    coarse: CoarseSpace::None,
                    ..PreconditionerOptions::default()
                },
            )
            .expect("preconditioner");
            let (u, p, rep) = fgmres_saddle(
                &r.sys,
                &mut pre,
                &r.none,
                &r.six,
                GradientMode::Corrected,
                &f,
                &r.sys.zeros_p(),
                FgmresOptions {
                    max_iterations: 2000,
                    stall_window: 2000,
                    restart,
                    ..FgmresOptions::default()
                },
            )
            .expect("FGMRES runs");
            let comps = r.six.components(&r.sys, &u);
            let tra = comps[0..3].iter().fold(0.0f64, |a, c| a.max(c.abs()));
            let rot = comps[3..6].iter().fold(0.0f64, |a, c| a.max(c.abs()));
            // Residual after removing the six-mode part of the velocity,
            // through the same operator the solver measures itself against.
            let mut u2 = u.clone();
            r.six.project(&r.sys, &mut u2);
            let (mu2, mp2) = r.sys.apply(&r.none, GradientMode::Corrected, &u2, &p).expect("apply");
            let ru: Vec<f64> = f.iter().zip(mu2.iter()).map(|(a, c)| a - c).collect();
            let rp: Vec<f64> = mp2.iter().map(|c| -c).collect();
            println!(
                "  n {n} restart {restart:>3}: outer {:>4} cycles {:>3} true_rel {:.4e}  \
                 after removing six-mode part {:.4e}  mode content tra {tra:.3e} rot {rot:.3e}",
                rep.outer_iterations,
                rep.restart_cycles,
                rep.true_rel_residual,
                r.sys.residual_norm(&ru, &rp) / b_norm,
            );
        }
    }

    /// Is the `n = 8` number a FLOOR or an unconverged iterate?
    ///
    /// `MCQ1_THIEULOT_FINDINGS.md` section 6.6 reports `6.981e-2` after 2000
    /// outer at `n = 8, layers = 6` and reads the sequence
    /// `1.638e-4, 1.338e-4, 6.981e-2` as "the floor gets worse with
    /// refinement". A floor and a slow descent are different things and the
    /// tail of the residual history distinguishes them, so the tail is
    /// printed here. Opt-in (`L1_DIAG_N8=1`) because it is the expensive run.
    #[test]
    fn the_n8_number_is_a_floor_or_it_is_not() {
        if std::env::var("L1_DIAG_N8").is_err() {
            println!("  skipped: set L1_DIAG_N8=1 to run the n = 8 rung");
            return;
        }
        let case = ThieulotCase::isoviscous_reference();
        for (n, layers) in [(2usize, 2usize), (4, 3), (8, 6)] {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            let f = force_of_degree(&topo, &case, 1);
            // RESTART is the discriminator. A flat residual proves stagnation,
            // not that the minimum over the reachable subspace was reached:
            // restarted GMRES(m) is well known to flatline above its own
            // minimum. Widening the Krylov subspace separates the two - a run
            // sitting on the true constrained minimum does not move, a run
            // stagnating on the restart does.
            for restart in [60usize, 200] {
                let mut pre = BlockPreconditioner::build(
                    &topo,
                    &op,
                    &r.sys,
                    PreconditionerOptions {
                        coarse: CoarseSpace::None,
                        ..PreconditionerOptions::default()
                    },
                )
                .expect("preconditioner");
                let (_u, _p, rep) = fgmres_saddle(
                    &r.sys,
                    &mut pre,
                    &r.none,
                    &r.six,
                    GradientMode::Corrected,
                    &f,
                    &r.sys.zeros_p(),
                    FgmresOptions {
                        max_iterations: 2000,
                        stall_window: 2000,
                        restart,
                        ..FgmresOptions::default()
                    },
                )
                .expect("FGMRES runs");
                let tail: Vec<String> = rep
                    .history
                    .iter()
                    .rev()
                    .take(4)
                    .rev()
                    .map(|v| format!("{v:.4e}"))
                    .collect();
                println!(
                    "  n {n} layers {layers} restart {restart:>3}: outer {:>4} cycles {:>3} \
                     inner_cg {:>7} true_rel {:.4e} mom {:.4e} cont {:.4e} mom/cont {:.3} \
                     last4 {tail:?}",
                    rep.outer_iterations,
                    rep.restart_cycles,
                    rep.inner_cg_iterations,
                    rep.true_rel_residual,
                    rep.momentum_rel,
                    rep.continuity_rel,
                    rep.momentum_rel / rep.continuity_rel,
                );
            }
        }
    }

    /// CANDIDATE REMEDIES, evaluated before anyone writes solver code.
    ///
    /// Two independent axes, crossed:
    ///
    /// - ROWS: which mode-tested momentum rows the discrete system deletes.
    ///   `Deflation::Rotations` and friends already do this in
    ///   `fgmres_saddle`, and each variant is a DIFFERENT discrete system,
    ///   which is why the answer is compared against the analytic free-slip
    ///   field for every row and not only the residual.
    /// - COLUMNS: which directions the solver's iterate can reach. `none` is
    ///   a hypothetical solver whose preconditioner can produce the six
    ///   modes; `six` is the solver we actually have, whose preconditioner
    ///   re-projects every `A`-block solve onto `null(K)^perp_{M_u}`.
    ///
    /// The `(rows none, cols six)` cell is today's solver. Any proposed fix
    /// has to show a small residual AND a velocity that has not moved away
    /// from the reference, and this table shows both for every cell.
    #[test]
    fn the_candidate_remedies_evaluated_densely() {
        let case = ThieulotCase::isoviscous_reference();
        for (n, layers) in rungs() {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            let dense = DenseSaddle::assemble(&r.sys, GradientMode::Corrected).expect("dense");
            let free = FreeSlipCase::new(case).expect("free-slip reference");
            let (u_free, geom) =
                area_mean_face_field(&topo, |x| free.velocity_cartesian(x)).expect("sampling");
            assert!(geom < 1e-12, "facet vertex cycle mismatch {geom:.3e}");
            let free_norm = r.sys.norm_u(&u_free);
            let q_six = dense.packed_modes(&r.six);
            let q_rot = dense.packed_modes(&r.rot);
            let q_tra = dense.packed_modes(&r.tra);
            println!(
                "\n=== candidate remedies, n {n} layers {layers}, l = 1 forcing. \
                 |u_free|_Mu {free_norm:.4e}"
            );
            println!(
                "  rows = which mode-tested momentum rows the system DELETES (a different \
                 discrete system each time)"
            );
            println!(
                "  cols = which directions the solver's velocity iterate can REACH \
                 (`six` is today's block preconditioner)"
            );
            let f = force_of_degree(&topo, &case, 1);
            let zp = r.sys.zeros_p();
            let bhat0 = dense.pack_rhs(&f, &zp);
            for (row_label, row_q) in [
                ("rows none    ", &Vec::new()),
                ("rows rot(3)  ", &q_rot),
                ("rows tra(3)  ", &q_tra),
                ("rows six     ", &q_six),
            ] {
                let a_row = project_rows(&dense.m, row_q);
                let b_row = project_vec(&bhat0, row_q);
                let bn = nrm2(&b_row);
                for (col_label, col_q) in
                    [("cols all ", &Vec::new()), ("cols six^perp", &q_six)]
                {
                    let a = restrict_columns(&a_row, col_q);
                    let qr = pivoted_qr(&a, &b_row);
                    let d = qr.rdiag();
                    let top = d.first().copied().unwrap_or(0.0);
                    let rank = d.iter().take_while(|v| **v > 1e-12 * top).count();
                    let res = qr.residual_at_rank(rank) / bn;
                    let y = project_vec(&qr.solve(rank), col_q);
                    let (u, mut p) = dense.unpack(&y, topo.face_count());
                    r.sys.remove_mean_pressure(&mut p).expect("gauge");
                    let comps = r.six.components(&r.sys, &u);
                    let tra = comps[0..3].iter().fold(0.0f64, |a, c| a.max(c.abs()));
                    let rot = comps[3..6].iter().fold(0.0f64, |a, c| a.max(c.abs()));
                    let diff: Vec<f64> =
                        u.iter().zip(u_free.iter()).map(|(a, b)| a - b).collect();
                    let err = r.sys.norm_u(&diff) / free_norm;
                    // The same error after removing the rigid-rotation part of
                    // the solved field: the continuum free-slip problem is
                    // only determined up to a rotation and the reference took
                    // the b = 0 gauge, so a solved field carrying a rotation
                    // is penalised by |u - u_free| for a gauge choice.
                    let mut u_norot = u.clone();
                    r.rot.project(&r.sys, &mut u_norot);
                    let diff2: Vec<f64> = u_norot
                        .iter()
                        .zip(u_free.iter())
                        .map(|(a, b)| a - b)
                        .collect();
                    let err2 = r.sys.norm_u(&diff2) / free_norm;
                    println!(
                        "  {row_label} {col_label}: rank {rank:>4} res {res:.3e}  \
                         |u-u_free|/|u_free| {err:.4e}  (rotation-gauged {err2:.4e})  \
                         mode content tra {tra:.2e} rot {rot:.2e}"
                    );
                }
            }
        }
    }

    /// The floor the actual solver reaches, at the same rungs, so every
    /// number above has the thing it is being compared against measured in
    /// the same run rather than quoted from another packet's table.
    #[test]
    fn the_fgmres_floor_at_the_same_rungs() {
        let case = ThieulotCase::isoviscous_reference();
        for (n, layers) in rungs() {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            for l in [1usize, 2, 3, 4, 5] {
                let f = force_of_degree(&topo, &case, l);
                let mut pre = BlockPreconditioner::build(
                    &topo,
                    &op,
                    &r.sys,
                    PreconditionerOptions {
                        coarse: CoarseSpace::None,
                        ..PreconditionerOptions::default()
                    },
                )
                .expect("preconditioner");
                let (u, _p, rep) = fgmres_saddle(
                    &r.sys,
                    &mut pre,
                    &r.none,
                    &r.six,
                    GradientMode::Corrected,
                    &f,
                    &r.sys.zeros_p(),
                    FgmresOptions {
                        max_iterations: 2000,
                        stall_window: 2000,
                        ..FgmresOptions::default()
                    },
                )
                .expect("FGMRES runs");
                let comps = r.six.components(&r.sys, &u);
                println!(
                    "  n {n} l {l}: outer {:>4} converged {:>5} true_rel {:.3e} mom {:.3e} \
                     cont {:.3e}  solved mode content tra {:.2e} rot {:.2e}",
                    rep.outer_iterations,
                    rep.converged,
                    rep.true_rel_residual,
                    rep.momentum_rel,
                    rep.continuity_rel,
                    comps[0..3].iter().fold(0.0f64, |a, c| a.max(c.abs())),
                    comps[3..6].iter().fold(0.0f64, |a, c| a.max(c.abs())),
                );
            }
        }
    }

    // -----------------------------------------------------------------
    // Section 9: is the 33 % velocity error a POINT-SAMPLED right-hand side?

    /// The rungs the dense Task-1 comparison runs on: `(2,2)` always, `(4,3)`
    /// under `L1_DIAG_N4`, `(8,6)` under `AL_RHS_N8_DENSE`. The `n = 8` rung
    /// is declined by default and the findings say why rather than omitting
    /// it silently.
    fn rhs_rungs() -> Vec<(usize, usize)> {
        let mut v = rungs();
        if std::env::var("AL_RHS_N8_DENSE").is_ok() {
            v.push((8, 6));
        }
        v
    }

    /// `M_u`-weighted `L2` and `L-infinity` of `u - u_ref` over a face
    /// subset.
    ///
    /// **This is deliberately the same convention as `qualify_thieulot`'s own
    /// private `subset_error`**, so the rows below are comparable with MC-Q1's
    /// tables rather than a second, differently normalised definition of the
    /// same words: both relatives use the REFERENCE field's size on the SAME
    /// subset, because a global denominator lets a subset where the reference
    /// is tiny look accurate. That choice does make the relative numbers
    /// incomparable BETWEEN subsets, so the ABSOLUTE `L-infinity` and the
    /// subset's share of the global weighted squared error are returned too -
    /// those are what answer "is the outlier the wall or the seam".
    struct SubsetError {
        l2_rel: f64,
        linf_rel: f64,
        linf_abs: f64,
        /// `sum_{f in S} M_u_f (u_f - u_ref_f)^2`, unnormalised.
        sq: f64,
        count: usize,
    }

    fn subset_error(
        sys: &SaddleSystem,
        keep: &dyn Fn(usize) -> bool,
        u: &[f64],
        u_ref: &[f64],
    ) -> SubsetError {
        let m = sys.m_u();
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        let mut linf = 0.0f64;
        let mut rmax = 0.0f64;
        let mut count = 0usize;
        for f in 0..u.len() {
            if !keep(f) {
                continue;
            }
            count += 1;
            let d = u[f] - u_ref[f];
            num += m[f] * d * d;
            den += m[f] * u_ref[f] * u_ref[f];
            linf = linf.max(d.abs());
            rmax = rmax.max(u_ref[f].abs());
        }
        SubsetError {
            l2_rel: if den > 0.0 { (num / den).sqrt() } else { f64::NAN },
            linf_rel: if rmax > 0.0 { linf / rmax } else { f64::NAN },
            linf_abs: linf,
            sq: num,
            count,
        }
    }

    /// Worst absolute deviation between two face fields, relative to the
    /// second one's largest entry. Used only for the CONTROLS, where the
    /// expected answer is exactly zero.
    fn face_deviation(a: &[f64], b: &[f64]) -> f64 {
        let scale = b.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
        a.iter()
            .zip(b.iter())
            .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()))
            / scale
    }

    /// TASK 1. Does an AREA-MEAN right-hand side move the 33 % velocity
    /// error, or is the 33 % somewhere else entirely?
    ///
    /// `L1_FLOOR_DIAGNOSIS.md` section 8 item 8 left this open: the solver's
    /// `b` is [`sampled_face_field`], a POINT value at the facet centroid,
    /// while the velocity unknown is defined as an area MEAN and
    /// [`area_mean_face_field`] computes that mean exactly by 4x4 Gauss on
    /// the bilinear patch. A point sample differs from the mean by `O(h^2)`,
    /// so it perturbs the ANSWER at `O(h^2)` - which is the same order as the
    /// discretization error the 33 % is being blamed on.
    ///
    /// **Both right-hand sides are solved DENSELY AND DIRECTLY**, not by
    /// FGMRES, so the `1.3e-4` preconditioner floor of sections 1-2 cannot
    /// contaminate the comparison and no iteration count enters it. The
    /// point-sampled column is therefore the `rows none cols all` row of
    /// section 3, and reproducing that row's numbers is this test's strongest
    /// instrument check.
    ///
    /// **The negative criterion.** Section 3 measures the exact discrete
    /// solution at `3.3180e-1` (rotation-gauged) at `n = 4`. Nothing here
    /// changes the discretization - only the right-hand side quadrature - so
    /// a LARGE improvement would be evidence of a broken instrument before it
    /// was evidence of a fix.
    #[test]
    fn the_area_mean_right_hand_side_against_the_point_sample() {
        let case = ThieulotCase::isoviscous_reference();
        let mut hs: Vec<f64> = Vec::new();
        let mut e_point: Vec<f64> = Vec::new();
        let mut e_area: Vec<f64> = Vec::new();
        let mut g_point: Vec<f64> = Vec::new();
        let mut g_area: Vec<f64> = Vec::new();
        let mut rhs_gap: Vec<f64> = Vec::new();
        let mut sol_gap: Vec<f64> = Vec::new();
        for (n, layers) in rhs_rungs() {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            let dense = DenseSaddle::assemble(&r.sys, GradientMode::Corrected).expect("dense");
            let mf = dense
                .matrix_free_disagreement(&r.sys, &r.none, 4, 0x5eed_1234)
                .expect("matrix-free cross-check");
            let free = FreeSlipCase::new(case).expect("free-slip reference");
            let (u_free, geom_v) =
                area_mean_face_field(&topo, |x| free.velocity_cartesian(x)).expect("sampling");
            assert!(geom_v < 1e-12, "facet vertex cycle mismatch {geom_v:.3e}");
            let free_norm = r.sys.norm_u(&u_free);
            let lab = face_labels(&topo);

            // ONE closure, two quadratures. `body_force_faces` builds exactly
            // this integrand internally (`rho e_r`, `rho = case.density`), so
            // the point-sampled column below IS the solver's own b.
            let force = |x: [f64; 3]| {
                let (rr, th) = ThieulotCase::spherical_of(x);
                let rho = case.density(rr, th);
                [rho * x[0] / rr, rho * x[1] / rr, rho * x[2] / rr]
            };
            let b_point = sampled_face_field(&topo, force);
            let (b_area, geom_f) = area_mean_face_field(&topo, force).expect("area-mean forcing");
            // THE CONTROL THE PACKET ASKED FOR, printed every run: the same
            // forcing through the two routes the rest of the tree already
            // trusts. `body_force_faces` is the qualification module's
            // assembler and `force_of_degree` is this module's l = 1 forcing;
            // both must be BIT-IDENTICAL to the local point sample, because
            // all three are the same point sample of the same integrand.
            let b_trusted = body_force_faces(&topo, |rr, th| case.density(rr, th));
            let b_module = force_of_degree(&topo, &case, 1);
            let ctl_trusted = face_deviation(&b_point, &b_trusted);
            let ctl_module = face_deviation(&b_point, &b_module);

            let npoint = r.sys.norm_u(&b_point);
            let narea = r.sys.norm_u(&b_area);
            let bdiff: Vec<f64> = b_area
                .iter()
                .zip(b_point.iter())
                .map(|(a, b)| a - b)
                .collect();
            let bgap = r.sys.norm_u(&bdiff) / npoint;

            println!(
                "\n=== area-mean vs point-sampled RHS, n {n} layers {layers}: \
                 interior {} cells {} dense {}x{} h_ref {:.6e}",
                dense.nu,
                dense.np,
                dense.dim(),
                dense.dim(),
                r.sys.mesh_scale()
            );
            println!(
                "  CONTROLS  matrix-free vs dense {mf:.3e} | facet vertex cycle vel {geom_v:.3e} \
                 force {geom_f:.3e} | point sample vs body_force_faces {ctl_trusted:.3e} \
                 vs force_of_degree {ctl_module:.3e}"
            );
            assert!(
                ctl_trusted == 0.0 && ctl_module == 0.0,
                "the local point sample is not the module's own forcing"
            );
            println!(
                "  |b|_Mu point {npoint:.6e} area {narea:.6e} ratio {:.6}  \
                 |b_area - b_point|_Mu / |b_point|_Mu {bgap:.4e}",
                narea / npoint
            );

            // DENSE DIRECT SOLVE, one per right-hand side. Unconstrained
            // least squares with the same rank ladder section 1.1 showed to
            // be invariant over six decades of cut.
            let zp = r.sys.zeros_p();
            let solve = |b_u: &[f64]| -> (Vec<f64>, Vec<f64>, f64, usize) {
                let bhat = dense.pack_rhs(b_u, &zp);
                let bn = nrm2(&bhat);
                let qr = pivoted_qr(&dense.m, &bhat);
                let d = qr.rdiag();
                let top = d.first().copied().unwrap_or(0.0);
                let rank = d.iter().take_while(|v| **v > 1e-12 * top).count();
                let res = qr.residual_at_rank(rank) / bn;
                let y = qr.solve(rank);
                let (u, mut p) = dense.unpack(&y, topo.face_count());
                r.sys.remove_mean_pressure(&mut p).expect("pressure gauge");
                (u, p, res, rank)
            };
            let (u_pt, p_pt, res_pt, rank_pt) = solve(&b_point);
            let (u_ar, p_ar, res_ar, rank_ar) = solve(&b_area);

            let report = |tag: &str, u: &[f64], p: &[f64], res: f64, rank: usize| -> (f64, f64) {
                let diff: Vec<f64> = u.iter().zip(u_free.iter()).map(|(a, b)| a - b).collect();
                let err = r.sys.norm_u(&diff) / free_norm;
                let mut u_norot = u.to_vec();
                r.rot.project(&r.sys, &mut u_norot);
                let diff2: Vec<f64> = u_norot
                    .iter()
                    .zip(u_free.iter())
                    .map(|(a, b)| a - b)
                    .collect();
                let err2 = r.sys.norm_u(&diff2) / free_norm;
                let comps = r.six.components(&r.sys, u);
                println!(
                    "  {tag}: rank {rank:>4} res {res:.3e}  |u-u_free|/|u_free| {err:.4e}  \
                     (rotation-gauged {err2:.4e})  |u|_Mu {:.4e} |p|_Mp {:.4e}  \
                     mode content tra {:.2e} rot {:.2e}",
                    r.sys.norm_u(u),
                    r.sys.norm_p(p),
                    comps[0..3].iter().fold(0.0f64, |a, c| a.max(c.abs())),
                    comps[3..6].iter().fold(0.0f64, |a, c| a.max(c.abs())),
                );
                (err, err2)
            };
            let (ep, gp) = report("point-sampled RHS", &u_pt, &p_pt, res_pt, rank_pt);
            let (ea, ga) = report("area-mean     RHS", &u_ar, &p_ar, res_ar, rank_ar);

            let udiff: Vec<f64> = u_ar.iter().zip(u_pt.iter()).map(|(a, b)| a - b).collect();
            let pdiff: Vec<f64> = p_ar.iter().zip(p_pt.iter()).map(|(a, b)| a - b).collect();
            let sgap = r.sys.norm_u(&udiff) / r.sys.norm_u(&u_pt);
            println!(
                "  the two SOLUTIONS differ by |u_area - u_point|_Mu / |u_point|_Mu {sgap:.4e}  \
                 (/|u_free|_Mu {:.4e})   pressure |dp|_Mp/|p_point|_Mp {:.4e}",
                r.sys.norm_u(&udiff) / free_norm,
                r.sys.norm_p(&pdiff) / r.sys.norm_p(&p_pt)
            );

            // Split by cell class. MC-Q1 found the L-infinity outlier is the
            // WALL, not the seam, so the wall row is the one to watch.
            let subsets: Vec<(&str, Box<dyn Fn(usize) -> bool + '_>)> = vec![
                ("global   ", Box::new(|i: usize| lab.interior[i])),
                (
                    "panel-int",
                    Box::new(|i: usize| {
                        lab.interior[i] && !lab.wall[i] && lab.class[i] == FacetClass::PanelInterior
                    }),
                ),
                (
                    "seam     ",
                    Box::new(|i: usize| {
                        lab.interior[i] && !lab.wall[i] && lab.class[i] != FacetClass::PanelInterior
                    }),
                ),
                (
                    "wall-adj ",
                    Box::new(|i: usize| lab.interior[i] && lab.wall[i]),
                ),
            ];
            let global_sq = subset_error(
                &r.sys,
                &|i: usize| lab.interior[i],
                &u_pt,
                &u_free,
            )
            .sq;
            let global_sq_a =
                subset_error(&r.sys, &|i: usize| lab.interior[i], &u_ar, &u_free).sq;
            for (label, keep) in subsets.iter() {
                let sp = subset_error(&r.sys, keep.as_ref(), &u_pt, &u_free);
                let sa = subset_error(&r.sys, keep.as_ref(), &u_ar, &u_free);
                let sd = subset_error(&r.sys, keep.as_ref(), &u_ar, &u_pt);
                println!(
                    "    {label} faces {:>5}  L2rel point {:.4e} area {:.4e} \
                     (change {:+.3e})  Linf_rel point {:.4e} area {:.4e}  \
                     Linf_ABS point {:.4e} area {:.4e}  share of global sq err \
                     point {:.4} area {:.4}  |u_a-u_p|/|u_p| here {:.3e}",
                    sp.count,
                    sp.l2_rel,
                    sa.l2_rel,
                    sa.l2_rel - sp.l2_rel,
                    sp.linf_rel,
                    sa.linf_rel,
                    sp.linf_abs,
                    sa.linf_abs,
                    sp.sq / global_sq,
                    sa.sq / global_sq_a,
                    sd.l2_rel
                );
            }

            hs.push(r.sys.mesh_scale());
            e_point.push(ep);
            e_area.push(ea);
            g_point.push(gp);
            g_area.push(ga);
            rhs_gap.push(bgap);
            sol_gap.push(sgap);
        }

        // Convergence order. With two rungs this is ONE consecutive-pair
        // slope and not a fitted exponent, and it is labelled as such; the
        // aggregate fit is printed only when three or more rungs ran, and it
        // names them.
        println!(
            "\n=== convergence, rungs h_ref {:?}",
            hs.iter().map(|v| format!("{v:.4e}")).collect::<Vec<_>>()
        );
        for (tag, ys) in [
            ("|u-u_free| point-sampled", &e_point),
            ("|u-u_free| area-mean    ", &e_area),
            ("gauged     point-sampled", &g_point),
            ("gauged     area-mean    ", &g_area),
            ("|b_a-b_p|/|b_p|         ", &rhs_gap),
            ("|u_a-u_p|/|u_p|         ", &sol_gap),
        ] {
            let ps = pair_slopes(&hs, ys);
            let fit = if hs.len() >= 3 {
                format!(" fitted exponent {:+.3} over {} rungs", slope(&hs, ys), hs.len())
            } else {
                format!(
                    " (2 rungs: pair slope only, no fit offered on {} points)",
                    hs.len()
                )
            };
            println!(
                "  {tag}: {:?}  pair slopes {:?}{fit}",
                ys.iter().map(|v| format!("{v:.4e}")).collect::<Vec<_>>(),
                ps.iter().map(|v| format!("{v:+.3}")).collect::<Vec<_>>()
            );
        }
        if !std::env::var("AL_RHS_N8_DENSE").is_ok() {
            println!(
                "  n = 8 dense DECLINED: 8832x8832, ~620 MB resident and ~1.3 GB live \
                 during a pivoted QR, and TWO factorisations are needed here (one per \
                 right-hand side) at ~6.9e11 flops each. Set AL_RHS_N8_DENSE=1 to run it. \
                 The cheap n = 8 number - the size of the RHS perturbation itself - is in \
                 the_rhs_quadrature_gap_scaling_in_h."
            );
        }
    }

    /// The cheap half of the `n = 8` question: how big is the point-sample /
    /// area-mean gap in the FORCING, on a ladder that reaches `n = 8` and
    /// `n = 16` without any dense factorisation at all.
    ///
    /// This does not answer what the gap does to the ANSWER - only a solve
    /// does that, and section 9 has it at two rungs. It answers whether the
    /// perturbation is `O(h^2)` as the quadrature theory says, which is what
    /// makes the two-rung solve result extrapolatable.
    #[test]
    fn the_rhs_quadrature_gap_scaling_in_h() {
        let case = ThieulotCase::isoviscous_reference();
        let mut ladder = vec![(2usize, 2usize), (4, 3)];
        if std::env::var("L1_DIAG_N8").is_ok() {
            ladder.push((8, 6));
        }
        if std::env::var("AL_DIAG_N16").is_ok() {
            ladder.push((16, 12));
        }
        let mut hs: Vec<f64> = Vec::new();
        let mut gaps: Vec<f64> = Vec::new();
        println!("\n=== |b_area - b_point| on the refinement ladder (no solve)");
        for (n, layers) in ladder {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let sys = SaddleSystem::build(&topo, &op, case.mu0(), FaceGradientRule::NonOrthogonalAndSkew)
                .expect("saddle system");
            let force = |x: [f64; 3]| {
                let (rr, th) = ThieulotCase::spherical_of(x);
                let rho = case.density(rr, th);
                [rho * x[0] / rr, rho * x[1] / rr, rho * x[2] / rr]
            };
            let b_point = sampled_face_field(&topo, force);
            let (b_area, geom) = area_mean_face_field(&topo, force).expect("area-mean forcing");
            let b_trusted = body_force_faces(&topo, |rr, th| case.density(rr, th));
            let ctl = face_deviation(&b_point, &b_trusted);
            assert!(ctl == 0.0, "the local point sample is not body_force_faces");
            let diff: Vec<f64> = b_area
                .iter()
                .zip(b_point.iter())
                .map(|(a, b)| a - b)
                .collect();
            let gap = sys.norm_u(&diff) / sys.norm_u(&b_point);
            println!(
                "  n {n:>2} layers {layers:>2} h_ref {:.6e}  |b_point|_Mu {:.6e} \
                 |b_area|_Mu {:.6e}  gap {gap:.4e}   CONTROLS vertex cycle {geom:.2e} \
                 vs body_force_faces {ctl:.2e}",
                sys.mesh_scale(),
                sys.norm_u(&b_point),
                sys.norm_u(&b_area),
            );
            hs.push(sys.mesh_scale());
            gaps.push(gap);
        }
        let ps = pair_slopes(&hs, &gaps);
        println!(
            "  pair slopes {:?}{}",
            ps.iter().map(|v| format!("{v:+.3}")).collect::<Vec<_>>(),
            if hs.len() >= 3 {
                format!(
                    "  fitted exponent {:+.3} over {} rungs",
                    slope(&hs, &gaps),
                    hs.len()
                )
            } else {
                "  (2 rungs: pair slope only)".to_string()
            }
        );
    }

    // -----------------------------------------------------------------
    // Section 10: the augmented-Lagrangian screening section 7.2 specified

    /// `gamma` grid, in the units section 7.2 asked for: multiples of
    /// `eta / h_ref^2`. The two decades below `1e-2` are added because the
    /// interesting question - at what `gamma` does a mode LEAVE the null
    /// space - has its answer at the bottom of the range, not the top, and a
    /// sweep that starts where the modes have already lifted cannot see it.
    /// The half-decade points exist because the first run showed
    /// `sigma_min(B_gamma)` SATURATES while `sigma_max` keeps growing
    /// linearly, so `kappa` has an interior minimum, and a decade grid can
    /// only bracket it to a decade.
    const GAMMA_DECADES: [f64; 10] =
        [1e-4, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1e0, 3e0, 1e1, 1e2];

    /// TASK 2. Does `K_gamma = K + gamma D^T M_p D` put the six `A`-null
    /// modes back into the preconditioner's reachable set - and in
    /// particular, does it do so for the ROTATIONS?
    ///
    /// **No solver code changes.** The augmented operator is applied through
    /// [`augmented_a_apply`], which is `SaddleSystem::viscous` minus
    /// `gamma SaddleSystem::gradient_implicit(SaddleSystem::divergence(.))`;
    /// `stokes_saddle.rs` is untouched. `g = 0` on this problem, so there is
    /// no right-hand-side change to model either.
    ///
    /// The concern section 7.2 stated: augmentation is nonsingular on the
    /// three truncated translations because `D t = O(1)`, but `D r = O(h^2)`
    /// for the rotations, so their Rayleigh quotient under the augmentation
    /// is only `gamma |D r|^2 = O(gamma h^4)` and they may stay effectively
    /// null. That is checked here mode by mode, named, and the `h` scaling of
    /// the rotation three is the subject of
    /// [`the_rotation_divergence_leak_scaling_in_h`].
    #[test]
    fn the_augmented_lagrangian_screening() {
        let case = ThieulotCase::isoviscous_reference();
        for (n, layers) in rungs() {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            let ab = DenseABlock::assemble(&r.sys).expect("dense A block");
            let h = r.sys.mesh_scale();
            let eta = r.sys.eta();
            let unit = eta / (h * h);

            // ---- CONTROLS, printed every run ----
            // 1. the adjoint identity the augmentation is built on
            let adj = adjoint_disagreement(&r.sys, &ab.iface, 4, 0x9e37_79b9).expect("adjoint");
            // 2. both dense blocks symmetric
            let (symb, symc) = ab.symmetry_defect();
            // 3. the dense blocks reproduce the matrix-free apply
            let mut worst_mf = 0.0f64;
            {
                let mut state: u64 = 0x1234_5678 | 1;
                let mut next = || {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    ((state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
                };
                for g in [0.0f64, 1.0, 37.0] {
                    let mut u = r.sys.zeros_u();
                    for &f in ab.iface.iter() {
                        u[f] = next();
                    }
                    let want_full = augmented_a_apply(&r.sys, g, &u).expect("augmented apply");
                    let want: Vec<f64> = ab
                        .iface
                        .iter()
                        .enumerate()
                        .map(|(j, &f)| ab.scale[j] * want_full[f])
                        .collect();
                    let got = ab.combine(g).mul_vec(&ab.pack(&u));
                    let d: Vec<f64> = want.iter().zip(got.iter()).map(|(a, b)| a - b).collect();
                    worst_mf = worst_mf.max(nrm2(&d) / nrm2(&want).max(1e-300));
                }
            }
            // 4. the local |D n_j| against the tree's own divergence_leak
            let leak = divergence_leak(&r.sys, &r.six).expect("divergence leak");
            let mut worst_leak = 0.0f64;
            let mut dnorm = Vec::with_capacity(6);
            for j in 0..r.six.dim() {
                let d = r.sys.divergence(r.six.mode(j)).expect("divergence");
                let v = r.sys.norm_p(&d) / r.sys.norm_u(r.six.mode(j));
                dnorm.push(v);
                worst_leak = worst_leak.max((v - leak.absolute[j]).abs() / leak.absolute[j].max(1e-300));
            }

            println!(
                "\n=== augmented-Lagrangian screening, n {n} layers {layers}: nu {} \
                 h_ref {h:.6e} eta {eta:.3e}  gamma unit eta/h^2 = {unit:.6e}",
                ab.nu
            );
            println!(
                "  CONTROLS  adjoint <G_imp p,u>_Mu + <p,Du>_Mp {adj:.3e} | \
                 |Bhat-Bhat^T|/|Bhat| {symb:.3e} |Chat-Chat^T|/|Chat| {symc:.3e} | \
                 dense vs matrix-free B_gamma {worst_mf:.3e} | local |D n| vs \
                 divergence_leak {worst_leak:.3e}"
            );
            assert!(adj < 1e-10, "the augmentation is not the adjoint pair D^T M_p D");
            assert!(symb < 1e-12 && symc < 1e-12, "an augmented block is not symmetric");

            // Which basis vector is which. The six-mode basis is built from
            // raw modes 0..3 = translations, 3..6 = rotations and then
            // Gram-Schmidted, so the labels below are CHECKED against the
            // translation-only and rotation-only bases rather than assumed.
            let names = ["tra_x", "tra_y", "tra_z", "rot_x", "rot_y", "rot_z"];
            for j in 0..r.six.dim() {
                let m = r.six.mode(j);
                let ot = (0..r.tra.dim())
                    .fold(0.0f64, |a, i| a.max(r.sys.dot_u(r.tra.mode(i), m).abs()));
                let orr = (0..r.rot.dim())
                    .fold(0.0f64, |a, i| a.max(r.sys.dot_u(r.rot.mode(i), m).abs()));
                println!(
                    "  mode {j} -> {:<5}  max|<.,translations>| {ot:.4}  \
                     max|<.,rotations>| {orr:.4}  |D n|_Mp/|n|_Mu {:.4e}  \
                     h|D n| {:.4e}  boundary share {:.4}",
                    names[j], dnorm[j], leak.dimensionless[j], leak.boundary_share[j]
                );
            }

            // The unaugmented baseline: what does the A block's spectrum look
            // like before anything is added?
            let svd0 = jacobi_svd(&ab.b, 30);
            let s0 = &svd0.sigma;
            let nu = ab.nu;
            let smax0 = s0[0];
            let s7 = s0[nu - 7];
            let svdc = jacobi_svd(&ab.c, 30);
            let cmax = svdc.sigma[0];
            // dim null(Chat) = the divergence-free velocity fields. On that
            // subspace the augmentation adds NOTHING, so `sigma_min(B_gamma)`
            // cannot exceed the smallest Rayleigh quotient of `Bhat` there no
            // matter how large `gamma` is - which is why `sigma_min`
            // saturates below and `kappa` has an interior minimum.
            let cnull = svdc.sigma.iter().filter(|v| **v <= 1e-12 * cmax).count();
            println!(
                "  gamma = 0 baseline: converged {} sweeps {}  sigma_max {smax0:.4e}  \
                 six smallest {:?}  smallest NONZERO (sigma_{}) {s7:.4e}  \
                 kappa on null^perp {:.4e}",
                svd0.converged,
                svd0.sweeps,
                s0[nu - 6..].iter().map(|v| format!("{v:.3e}")).collect::<Vec<_>>(),
                nu - 6,
                smax0 / s7
            );
            println!(
                "  sigma_max(Chat) {cmax:.4e}  dim null(Chat) {cnull} of {nu} \
                 (the discretely divergence-free fields; the augmentation is BLIND there)  \
                 =>  gamma that balances the two blocks \
                 sigma_max(Bhat)/sigma_max(Chat) = {:.4e}  = {:.4e} x eta/h^2 \
                 = {:.4e} x eta",
                smax0 / cmax,
                smax0 / cmax / unit,
                smax0 / cmax / eta
            );

            // ---- the sweep ----
            // Matrix-free mode norms at every gamma (cheap at any rung);
            // the dense spectrum at a subset unless AL_GAMMA_FULL is set,
            // because a 768x768 one-sided Jacobi per gamma is the whole cost
            // of this test.
            let dense_at: Vec<f64> = if std::env::var("AL_GAMMA_FULL").is_ok() || nu <= 256 {
                GAMMA_DECADES.to_vec()
            } else {
                vec![1e-2, 1e0, 1e2]
            };
            println!(
                "  dense spectrum computed at gamma/(eta/h^2) in {:?}{}",
                dense_at,
                if dense_at.len() < GAMMA_DECADES.len() {
                    "  (set AL_GAMMA_FULL=1 for the whole grid; the mode norms below are \
                     computed at EVERY gamma regardless)"
                } else {
                    ""
                }
            );
            println!(
                "  g/(eta/h2)  gamma       |B_g n|_Mu per mode: {}   then Rayleigh \
                 <B_g n,n>/<n,n> per mode",
                names.join("      ")
            );
            let mut sweep: Vec<(f64, f64, f64, f64)> = Vec::new();
            for &gd in GAMMA_DECADES.iter() {
                let gamma = gd * unit;
                let mut norms = Vec::with_capacity(6);
                let mut rayl = Vec::with_capacity(6);
                for j in 0..r.six.dim() {
                    let m = r.six.mode(j);
                    let y = augmented_a_apply(&r.sys, gamma, m).expect("augmented apply");
                    norms.push(r.sys.norm_u(&y) / r.sys.norm_u(m));
                    rayl.push(r.sys.dot_u(&y, m) / r.sys.dot_u(m, m));
                }
                println!(
                    "  {gd:>9.0e}  {gamma:.4e}  |B_g n| {:?}",
                    norms.iter().map(|v| format!("{v:.4e}")).collect::<Vec<_>>()
                );
                println!(
                    "                          Rayleigh {:?}",
                    rayl.iter().map(|v| format!("{v:.4e}")).collect::<Vec<_>>()
                );
                if dense_at.contains(&gd) {
                    let a = ab.combine(gamma);
                    let svd = jacobi_svd(&a, 30);
                    let s = &svd.sigma;
                    println!(
                        "                          dense: converged {} sweeps {}  \
                         sigma_max {:.4e} sigma_min {:.4e} kappa {:.4e}  \
                         eight smallest {:?}",
                        svd.converged,
                        svd.sweeps,
                        s[0],
                        s[nu - 1],
                        s[0] / s[nu - 1].max(1e-300),
                        s[nu - 8..].iter().map(|v| format!("{v:.3e}")).collect::<Vec<_>>()
                    );
                    // Are the three smallest singular values THE ROTATIONS?
                    // `<B_gamma r, r>_{M_u} = gamma |D r|^2` exactly (because
                    // `K r = 0`), so if the three smallest match that value
                    // the rotations ARE the bottom of the spectrum and the
                    // `kappa` a solver that gauges them away would see is
                    // `sigma_max / sigma_{nu-4}`. That identification is
                    // CHECKED, not assumed, and the deflated `kappa` is
                    // printed only when it holds.
                    let want_rot = gamma * dnorm[3] * dnorm[3];
                    let bottom_is_rot = (nu - 3..nu)
                        .all(|k| (s[k] - want_rot).abs() <= 0.25 * want_rot.max(1e-300));
                    if bottom_is_rot {
                        println!(
                            "                          the three smallest ARE the rotations \
                             (gamma|D r|^2 = {want_rot:.4e}); kappa with the rotations \
                             GAUGED AWAY = sigma_max / sigma_{} = {:.4e}  \
                             (baseline at gamma = 0 is {:.4e})",
                            nu - 4,
                            s[0] / s[nu - 4],
                            smax0 / s7
                        );
                    } else {
                        println!(
                            "                          the three smallest are NOT the \
                             rotations (gamma|D r|^2 = {want_rot:.4e} vs sigma_min \
                             {:.4e}); no rotation-gauged kappa reported at this gamma",
                            s[nu - 1]
                        );
                    }
                    sweep.push((gd, s[nu - 1], s[0], s[0] / s[nu - 1].max(1e-300)));
                }
            }
            // THE COST SIDE, in one line. The baseline to preserve is
            // `kappa` on `null(K)^perp`, which is what today's deflated inner
            // CG actually sees.
            if let Some(best) = sweep
                .iter()
                .min_by(|a, b| a.3.partial_cmp(&b.3).expect("finite kappa"))
            {
                let (gd, smin, smax, kap) = *best;
                println!(
                    "  COST: baseline kappa on null^perp {:.4e} (gamma = 0). Best kappa over \
                     the grid {kap:.4e} at gamma = {gd:.0e} x eta/h^2, where sigma_min \
                     {smin:.4e} sigma_max {smax:.4e}. Augmentation multiplies the \
                     conditioning the inner CG sees by {:.4e} AT BEST, and sigma_min \
                     saturates (largest sigma_min over the grid {:.4e}) while sigma_max \
                     grows linearly in gamma.",
                    smax0 / s7,
                    kap / (smax0 / s7),
                    sweep.iter().fold(0.0f64, |a, b| a.max(b.1))
                );
            }
        }
    }

    /// The rotation question, isolated and taken to `n = 8`: does `|D r|`
    /// fall like `h^2`?
    ///
    /// This is the number that decides whether the augmented Lagrangian can
    /// work on the rotations at all. `<B_gamma r, r>_{M_u} = gamma |D r|^2`
    /// exactly (because `K r = 0`), so a `|D r| ~ h^2` makes the augmented
    /// Rayleigh quotient `O(gamma h^4)` and the rotations stay effectively
    /// null under refinement for any fixed `gamma`. Matrix free, so `n = 8`
    /// costs a topology build and six divergences and is affordable here even
    /// though no dense factorisation at that rung is.
    #[test]
    fn the_rotation_divergence_leak_scaling_in_h() {
        let case = ThieulotCase::isoviscous_reference();
        let mut ladder = vec![(2usize, 2usize), (4, 3)];
        if std::env::var("L1_DIAG_N8").is_ok() {
            ladder.push((8, 6));
        }
        if std::env::var("AL_DIAG_N16").is_ok() {
            ladder.push((16, 12));
        }
        let names = ["tra_x", "tra_y", "tra_z", "rot_x", "rot_y", "rot_z"];
        let mut hs: Vec<f64> = Vec::new();
        let mut cols: Vec<Vec<f64>> = vec![Vec::new(); 6];
        println!("\n=== |D n_j|_Mp / |n_j|_Mu on the refinement ladder (matrix free)");
        for (n, layers) in ladder {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            let leak = divergence_leak(&r.sys, &r.six).expect("divergence leak");
            hs.push(r.sys.mesh_scale());
            for j in 0..6 {
                cols[j].push(leak.absolute[j]);
            }
            println!(
                "  n {n:>2} layers {layers:>2} h_ref {:.6e}  {:?}",
                r.sys.mesh_scale(),
                leak.absolute
                    .iter()
                    .map(|v| format!("{v:.4e}"))
                    .collect::<Vec<_>>()
            );
        }
        // The ratio is what actually sets the penalty. `<B_gamma n, n> =
        // gamma |D n|^2` exactly on an `A`-null mode, so the rotations are
        // weaker than the translations under the SAME augmentation by the
        // SQUARE of this ratio, and it is that square whose exponent section
        // 7.2 predicted to be `h^4`.
        let ratio: Vec<f64> = (0..hs.len()).map(|i| cols[3][i] / cols[0][i]).collect();
        let sq: Vec<f64> = ratio.iter().map(|v| v * v).collect();
        println!(
            "  |D rot| / |D tra| {:?}  pair slopes {:?}",
            ratio.iter().map(|v| format!("{v:.4e}")).collect::<Vec<_>>(),
            pair_slopes(&hs, &ratio)
                .iter()
                .map(|v| format!("{v:+.3}"))
                .collect::<Vec<_>>()
        );
        println!(
            "  its SQUARE (= Rayleigh_rot / Rayleigh_tra under augmentation, at any gamma) \
             {:?}  pair slopes {:?}{}",
            sq.iter().map(|v| format!("{v:.4e}")).collect::<Vec<_>>(),
            pair_slopes(&hs, &sq)
                .iter()
                .map(|v| format!("{v:+.3}"))
                .collect::<Vec<_>>(),
            if hs.len() >= 3 {
                format!(
                    "  fitted exponent {:+.3} over {} rungs (section 7.2 predicted h^4)",
                    slope(&hs, &sq),
                    hs.len()
                )
            } else {
                String::new()
            }
        );
        for j in 0..6 {
            let ps = pair_slopes(&hs, &cols[j]);
            println!(
                "  {:<5} pair slopes in h {:?}{}",
                names[j],
                ps.iter().map(|v| format!("{v:+.3}")).collect::<Vec<_>>(),
                if hs.len() >= 3 {
                    format!(
                        "  fitted exponent {:+.3} over {} rungs",
                        slope(&hs, &cols[j]),
                        hs.len()
                    )
                } else {
                    "  (2 rungs: pair slope only)".to_string()
                }
            );
        }
    }
}
