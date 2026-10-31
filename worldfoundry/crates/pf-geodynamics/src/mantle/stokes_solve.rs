//! MC0 solver prerequisite: the spectral conditioning of the assembled
//! viscous `A` block, and a deflated conjugate-gradient solve for it.
//!
//! [`super::stokes`] builds `A = -M_u^-1 E^T W E` and measures the LOCAL
//! coercivity of the strain fit. Section 3.5 of
//! `docs/reviews/planet-foundry-b1/MC0_CONDITIONING_GATES.md` names the one
//! thing that gate did not establish: the spectral condition number of the
//! ASSEMBLED block across refinement. A uniform local `c1` makes the standard
//! elliptic `O(h^-2)` the EXPECTED behaviour; expected is not measured, and
//! this module measures it. The write-up is
//! `docs/reviews/planet-foundry-b1/MC0_ABLOCK_SOLVE_FINDINGS.md`.
//!
//! Nothing in [`super::stokes`] is modified. This module only reads its
//! public surface (`stiffness_apply`, `cell_stencil`, `face_mass`,
//! `cell_weights`, `dense_stiffness`, `unknown_faces`).
//!
//! # One inner product, carried through every object
//!
//! `A` is self-adjoint in `M_u`, not in the Euclidean product, so every
//! object here carries the metric it lives in rather than assuming one:
//!
//! ```text
//! B     = P^-1 K,   K = E^T W E,   P a positive diagonal
//! <x,y>_P = sum_{f interior} P_f x_f y_f
//! <x, B y>_P = x^T K y = <B x, y>_P                 self-adjoint, always
//! range(B) = { v : <n_j, v>_P = 0 for the six null modes n_j }
//! ```
//!
//! With `P = M_u` ([`Metric::FaceMass`]) `B = -A` exactly, so the Lanczos
//! spectrum below is the spectrum of `A` and the CG iteration count below is
//! the count for UNPRECONDITIONED CG on `A`. With `P = diag(K)`
//! ([`Metric::StiffnessDiagonal`]) the same code is Jacobi-preconditioned CG,
//! and the same Lanczos code returns the preconditioned spectrum. That is the
//! whole reason the metric is a parameter: the two runs then differ in one
//! vector and in nothing else.
//!
//! The identity `range(B) = null(B)^perp_P` is what makes the deflation
//! correct, and it is an identity only in the `P` product. Projecting with a
//! DIFFERENT product leaves a persistent component that nothing fails on -
//! the same class of mistake `pf_grid::shell_ops::rigid_rotation_component`
//! documents for the rotation-fit weights. It is a negative control here
//! ([`tests::the_wrong_inner_product_leaves_a_persistent_null_component`])
//! and not a comment.
//!
//! # What is deflated, and why six
//!
//! Gate 1 measured `dim null(A) = 6` in the space that already excludes the
//! boundary faces: three rotations (rigid) and three TRUNCATED translations
//! (a translation on the interior faces, zero on the boundary faces - not a
//! rigid motion, and not divergence free). [`NullSpace::raw_modes`] builds
//! exactly those six and nothing else, and [`NullSpace::build`] makes them
//! `P`-orthonormal by twice-repeated modified Gram-Schmidt.
//!
//! **The near-null rotation modes of the SADDLE POINT are not touched here.**
//! The rotations are only `O(h^2)` divergence free, which is a property of
//! `[-A, G; D, 0]`, not of `A`; they are exactly null for `A` and are
//! deflated as such.
//!
//! # What was measured, in one block
//!
//! `M_u` inner product, i.e. `A` itself, over the self-similar ladder
//! `(n, layers) = (4,3), (8,6), (16,12), (32,24)`. Full tables and every
//! tolerance in `docs/reviews/planet-foundry-b1/MC0_ABLOCK_SOLVE_FINDINGS.md`.
//!
//! ```text
//!   n  lay  unknowns   h (m)     lambda_max   lambda_min      kappa   CG iters to 1e-8
//!   2    2       120  2.391e6    3.0743e-12   9.9260e-15   3.097e2         49
//!   4    3       768  1.406e6    1.1021e-11   2.6945e-14   4.090e2        131
//!   8    6      6528  7.150e5    3.2843e-11   1.4328e-14   2.292e3        286
//!  16   12     53760  3.590e5    1.2237e-10   3.1527e-14   3.882e3        553
//!  32   24    436224  1.797e5    4.7621e-10   3.8691e-14   1.231e4        990
//!
//! finest pair (16 -> 32):  lambda_max h^-1.963   lambda_min h^-0.296
//!                          kappa      h^-1.667   CG iters   h^-0.841
//!                          sqrt(kappa) predicts CG at h^-0.834
//! ```
//!
//! **`kappa(A)` is `O(h^-2)` and NOT worse**, so the `A` block does not need
//! multigrid on conditioning grounds; whether it needs multigrid on COST
//! grounds is a separate question and section 6 of the findings has the
//! extrapolated counts. The CG iteration count tracks `sqrt(kappa)` to `0.8%`
//! at the fine end, which is the check that the deflation is complete and the
//! inner product consistent.
//!
//! Every `lambda_min` above is an UPPER bound and every `kappa` a LOWER
//! bound - see [`SmallestEigenReport::residual_rel`] for the measured reason
//! no lower bound on `lambda_min` is claimed.
//!
//! # Honest scope
//!
//! - The manufactured-solution check ([`tests::the_manufactured_solution_is_recovered_modulo_the_null_space`])
//!   is a SELF-CONSISTENCY check of the solver against the operator. It
//!   verifies that CG inverts the operator this repository built. It says
//!   nothing about whether that operator is the right physics; the analytic
//!   comparison is MC-Q1 ([`super::benchmark_thieulot`]).
//! - Lanczos Ritz values are interior to the spectrum, so an unconverged run
//!   UNDERESTIMATES `lambda_max` and OVERESTIMATES `lambda_min`, i.e. it
//!   underestimates `kappa` in both directions at once. Every report
//!   therefore carries the Ritz residual bound `beta_m |e_m^T s_i|`, which is
//!   a rigorous bracket, alongside the half-length rerun, which is not.

use super::stokes::{strain_contraction, ViscousOperator};
use super::MantleError;
use pf_grid::shell_ops::ShellTopology;

/// The dimension of `null(A)` in this discrete space, measured in gate 1:
/// three rotations plus three truncated translations.
pub const NULL_DIM: usize = 6;

/// Which positive diagonal defines the inner product a solve runs in.
///
/// This is not a tuning knob. `A` is self-adjoint in exactly one of these
/// (`M_u`); the other is the Jacobi preconditioner, and running it means
/// solving the same `K u = b` in a different geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Metric {
    /// `M_u`, the over-relaxed face mass
    /// ([`pf_grid::shell_ops::ShellTopology::face_masses_over_relaxed`]).
    /// `B = P^-1 K = -A`, so this run measures `A` itself.
    FaceMass,
    /// `diag(K) = diag(E^T W E)`: the cheapest possible preconditioner, one
    /// number per face, assembled directly from the strain stencil.
    StiffnessDiagonal,
}

impl Metric {
    pub fn label(self) -> &'static str {
        match self {
            Metric::FaceMass => "M_u",
            Metric::StiffnessDiagonal => "diag(K)",
        }
    }
}

/// `diag(E^T W E)`, one entry per face, boundary entries zero.
///
/// `K_ff = sum_c W_c (dEps_c/du_f : dEps_c/du_f)`, read straight off the
/// stored strain stencil, so it costs one pass over the operator's own
/// coefficients rather than `n_faces` operator applications.
/// [`tests::the_stiffness_diagonal_matches_the_operator_it_claims_to_be`]
/// checks it against `stiffness_apply` on unit vectors.
pub fn stiffness_diagonal(op: &ViscousOperator) -> Vec<f64> {
    let e = op.strain();
    let mut d = vec![0.0f64; e.face_count()];
    for c in 0..e.cell_count() {
        let w = op.cell_weights()[c];
        for (f, coeff) in e.cell_stencil(c) {
            d[f as usize] += w * strain_contraction(coeff, coeff);
        }
    }
    d
}

/// The `A` block as a self-adjoint positive semidefinite operator in one
/// stated inner product, applied matrix-free.
///
/// Vectors are always FULL face-length with the boundary entries pinned to
/// zero, which is what the discrete space contains: the inner-product weight
/// is zero on a boundary face, so a boundary entry can neither be read nor
/// contribute, and there is no separate index map to get wrong.
#[derive(Debug)]
pub struct ABlock<'a> {
    op: &'a ViscousOperator,
    metric: Metric,
    /// `P_f`, zero on boundary faces.
    weight: Vec<f64>,
    /// `1 / P_f`, zero on boundary faces.
    inv_weight: Vec<f64>,
    unknowns: usize,
}

impl<'a> ABlock<'a> {
    /// Build the block in one metric. Fails loud if the metric's diagonal is
    /// not strictly positive on an interior face, which would mean the inner
    /// product is not one.
    pub fn new(
        topo: &ShellTopology,
        op: &'a ViscousOperator,
        metric: Metric,
    ) -> Result<Self, MantleError> {
        let raw = match metric {
            Metric::FaceMass => op.face_mass().to_vec(),
            Metric::StiffnessDiagonal => stiffness_diagonal(op),
        };
        if raw.len() != topo.face_count() {
            return Err(MantleError::Field {
                field: "inner-product diagonal",
                value: format!("length {}, want {}", raw.len(), topo.face_count()),
                why: "does not match the shell this operator was built for",
            });
        }
        let mut weight = vec![0.0f64; raw.len()];
        let mut inv_weight = vec![0.0f64; raw.len()];
        let mut unknowns = 0usize;
        for (index, face) in topo.faces().iter().enumerate() {
            if face.neighbor.is_none() {
                continue;
            }
            let w = raw[index];
            if !w.is_finite() || w <= 0.0 {
                return Err(MantleError::NonFinite {
                    field: "inner-product diagonal on an interior face",
                    index,
                    value: w,
                });
            }
            weight[index] = w;
            inv_weight[index] = 1.0 / w;
            unknowns += 1;
        }
        Ok(Self {
            op,
            metric,
            weight,
            inv_weight,
            unknowns,
        })
    }

    pub fn metric(&self) -> Metric {
        self.metric
    }

    pub fn face_count(&self) -> usize {
        self.weight.len()
    }

    /// How many face unknowns the solve actually carries.
    pub fn unknown_count(&self) -> usize {
        self.unknowns
    }

    /// A zero vector of the right length, boundary entries included.
    pub fn zeros(&self) -> Vec<f64> {
        vec![0.0f64; self.weight.len()]
    }

    /// `<a, b>_P`, summed over the interior faces only.
    pub fn dot(&self, a: &[f64], b: &[f64]) -> f64 {
        assert_eq!(a.len(), self.weight.len(), "left vector is not face-length");
        assert_eq!(
            b.len(),
            self.weight.len(),
            "right vector is not face-length"
        );
        a.iter()
            .zip(b.iter())
            .zip(self.weight.iter())
            .map(|((x, y), w)| w * x * y)
            .sum()
    }

    pub fn norm(&self, a: &[f64]) -> f64 {
        self.dot(a, a).sqrt()
    }

    /// `B u = P^-1 E^T W E u`, the positive semidefinite operator. With
    /// [`Metric::FaceMass`] this is exactly `-A`.
    ///
    /// Fails loud on a non-finite entry rather than returning a plausible
    /// vector: a solver that swallows a NaN reports a converged residual it
    /// computed from a NaN.
    pub fn apply(&self, u: &[f64]) -> Result<Vec<f64>, MantleError> {
        let mut k = self.op.stiffness_apply(u)?;
        for (index, (slot, iw)) in k.iter_mut().zip(self.inv_weight.iter()).enumerate() {
            *slot *= *iw;
            if !slot.is_finite() {
                return Err(MantleError::NonFinite {
                    field: "A-block operator application",
                    index,
                    value: *slot,
                });
            }
        }
        Ok(k)
    }

    /// `<u, B u>_P = u^T K u = Phi(u)`, the discrete dissipation. Equal to
    /// [`ViscousOperator::energy`] by construction; the equality is the
    /// cheapest available check that the metric and the operator agree.
    pub fn energy_form(&self, u: &[f64]) -> Result<f64, MantleError> {
        self.op.energy(u)
    }

    /// Fail loud on a non-finite entry of an iterate, naming the face.
    pub fn check_finite(&self, what: &'static str, v: &[f64]) -> Result<(), MantleError> {
        for (index, x) in v.iter().enumerate() {
            if !x.is_finite() {
                return Err(MantleError::NonFinite {
                    field: what,
                    index,
                    value: *x,
                });
            }
        }
        Ok(())
    }
}

/// The six exact null directions of `A`, orthonormal in one stated metric.
#[derive(Debug, Clone)]
pub struct NullSpace {
    modes: Vec<Vec<f64>>,
    metric: Metric,
    /// Smallest ratio `|mode after orthogonalisation| / |mode before|` over
    /// the six. A collapse here would mean the six raw modes are numerically
    /// dependent in this metric, which is a finding and not a warning.
    smallest_pivot: f64,
}

impl NullSpace {
    /// The six raw modes of gate 1, on full face-length vectors with the
    /// boundary entries zero:
    ///
    /// - three TRUNCATED translations, `u_f = (A_f . e_a) / |A_f|`, the exact
    ///   facet-mean normal component of a unit translation (the same
    ///   functional the operator's own design row is), truncated to zero on
    ///   the boundary faces;
    /// - three rotations, from
    ///   [`pf_grid::shell_ops::ShellTopology::rigid_rotation_field`], which
    ///   uses the exact facet moment and is already zero on the boundary.
    pub fn raw_modes(topo: &ShellTopology) -> Vec<Vec<f64>> {
        let mut out = Vec::with_capacity(NULL_DIM);
        for axis in 0..3 {
            out.push(
                topo.faces()
                    .iter()
                    .map(|f| {
                        if f.neighbor.is_none() {
                            0.0
                        } else {
                            f.area_vector_m2[axis] / f.area_m2
                        }
                    })
                    .collect(),
            );
        }
        for axis in 0..3 {
            let mut omega = [0.0f64; 3];
            omega[axis] = 1.0;
            out.push(topo.rigid_rotation_field(omega));
        }
        out
    }

    /// `P`-orthonormalise the six raw modes by modified Gram-Schmidt, run
    /// TWICE per vector: one pass is not enough when the raw modes are as
    /// differently scaled as a unit translation and a rotation whose entries
    /// carry the planet's radius, and "twice is enough" is the standard
    /// result for this.
    pub fn build(topo: &ShellTopology, block: &ABlock) -> Result<Self, MantleError> {
        let mut modes = Self::raw_modes(topo);
        let mut smallest_pivot = f64::MAX;
        for j in 0..NULL_DIM {
            let before = block.norm(&modes[j]);
            if !before.is_finite() || before <= 0.0 {
                return Err(MantleError::field(
                    "null mode norm",
                    before,
                    "a raw null mode has no content in this inner product",
                ));
            }
            let (head, tail) = modes.split_at_mut(j);
            let target = &mut tail[0];
            for _pass in 0..2 {
                for prev in head.iter() {
                    let c = block.dot(prev, target);
                    for (slot, x) in target.iter_mut().zip(prev.iter()) {
                        *slot -= c * x;
                    }
                }
            }
            let after = block.norm(target);
            let pivot = after / before;
            smallest_pivot = smallest_pivot.min(pivot);
            if !pivot.is_finite() || pivot <= 1e-8 {
                return Err(MantleError::field(
                    "null mode Gram-Schmidt pivot",
                    pivot,
                    "the six null modes are numerically dependent in this inner product",
                ));
            }
            let inv = 1.0 / after;
            for slot in target.iter_mut() {
                *slot *= inv;
            }
        }
        Ok(Self {
            modes,
            metric: block.metric(),
            smallest_pivot,
        })
    }

    pub fn dim(&self) -> usize {
        self.modes.len()
    }

    pub fn metric(&self) -> Metric {
        self.metric
    }

    pub fn smallest_pivot(&self) -> f64 {
        self.smallest_pivot
    }

    pub fn mode(&self, index: usize) -> &[f64] {
        &self.modes[index]
    }

    /// `v <- v - sum_j <n_j, v>_P n_j`, the `P`-orthogonal projection onto
    /// `range(B)`.
    pub fn project(&self, block: &ABlock, v: &mut [f64]) {
        for m in self.modes.iter() {
            let c = block.dot(m, v);
            for (slot, x) in v.iter_mut().zip(m.iter()) {
                *slot -= c * x;
            }
        }
    }

    /// `max_j |<n_j, v>_P| / |v|_P`: how much of `v` is still in the null
    /// space, relative to its own size. Zero vector reports zero.
    pub fn leak(&self, block: &ABlock, v: &[f64]) -> f64 {
        let nv = block.norm(v);
        if !nv.is_finite() || nv <= 0.0 {
            return 0.0;
        }
        self.modes
            .iter()
            .map(|m| (block.dot(m, v) / nv).abs())
            .fold(0.0f64, f64::max)
    }
}

// ---------------------------------------------------------------------
// Spectrum

/// One deflated Lanczos sweep, with the evidence needed to believe it.
#[derive(Debug, Clone)]
pub struct LanczosReport {
    pub metric: Metric,
    pub unknowns: usize,
    /// Lanczos steps actually taken.
    pub iterations: usize,
    pub lambda_max: f64,
    /// Smallest Ritz value, i.e. the smallest NON-null eigenvalue estimate.
    pub lambda_min: f64,
    pub kappa: f64,
    /// `beta_m |e_m^T s_i|` for the two extremal Ritz pairs, each divided by
    /// the Ritz value it brackets: a rigorous bracket, since some eigenvalue
    /// of `B` lies within `beta_m |e_m^T s_i|` of the Ritz value. So
    /// `bound_min_rel = 3e-4` reads "lambda_min is pinned to 3 parts in
    /// 10^4", and a bound of order 1 means the estimate is worthless.
    pub bound_max_rel: f64,
    pub bound_min_rel: f64,
    /// The same three numbers from the first half of the same run, as a
    /// convergence series. Not rigorous; reported next to the bound.
    pub half_lambda_max: f64,
    pub half_lambda_min: f64,
    pub half_kappa: f64,
    /// Worst null-space leak of any Lanczos vector over the sweep.
    pub worst_leak: f64,
    /// Step at which the recurrence broke down, if it did.
    pub breakdown_at: Option<usize>,
}

/// Deflated, fully reorthogonalised Lanczos for the extremal spectrum of
/// `B = P^-1 K` in the `P` inner product.
///
/// Full reorthogonalisation (twice per step) rather than a plain three-term
/// recurrence: without it the Ritz values duplicate the converged extremes
/// and the smallest Ritz value stops being an estimate of `lambda_min`. The
/// six null modes are projected out of every Lanczos vector, so the sweep
/// runs on `range(B)` and the smallest Ritz value is the smallest NON-null
/// eigenvalue rather than a roundoff-sized zero.
///
/// **The start vector is an INPUT, because it is part of the measurement.**
/// A purely random start overlaps the lowest eigenvector by about
/// `1/sqrt(N)`, so at the fine end of a ladder `lambda_min` converges last
/// and the run silently underestimates `kappa`. [`lanczos_start`] adds a
/// smooth component for exactly that reason, and the Ritz residual bound in
/// [`LanczosReport`] is what decides whether it worked.
pub fn lanczos_spectrum(
    block: &ABlock,
    null: &NullSpace,
    max_iters: usize,
    start: &[f64],
) -> Result<LanczosReport, MantleError> {
    if null.metric() != block.metric() {
        return Err(MantleError::Field {
            field: "Lanczos deflation metric",
            value: format!("{} vs {}", null.metric().label(), block.metric().label()),
            why: "the null space and the operator must live in the same inner product",
        });
    }
    let cap = block.unknown_count().saturating_sub(NULL_DIM);
    let m_max = max_iters.min(cap).max(1);

    let mut v = start.to_vec();
    block.check_finite("Lanczos start vector", &v)?;
    null.project(block, &mut v);
    let n0 = block.norm(&v);
    if !n0.is_finite() || n0 <= 0.0 {
        return Err(MantleError::field(
            "Lanczos start vector norm",
            n0,
            "the start vector is entirely inside the null space",
        ));
    }
    for slot in v.iter_mut() {
        *slot /= n0;
    }

    let mut basis: Vec<Vec<f64>> = vec![v];
    let mut alpha: Vec<f64> = Vec::with_capacity(m_max);
    let mut beta: Vec<f64> = Vec::with_capacity(m_max);
    let mut worst_leak = 0.0f64;
    let mut breakdown_at = None;

    for j in 0..m_max {
        let mut w = block.apply(&basis[j])?;
        null.project(block, &mut w);
        let a = block.dot(&basis[j], &w);
        if !a.is_finite() {
            return Err(MantleError::NonFinite {
                field: "Lanczos alpha",
                index: j,
                value: a,
            });
        }
        alpha.push(a);
        for (slot, x) in w.iter_mut().zip(basis[j].iter()) {
            *slot -= a * x;
        }
        if j > 0 {
            let b = beta[j - 1];
            for (slot, x) in w.iter_mut().zip(basis[j - 1].iter()) {
                *slot -= b * x;
            }
        }
        for _pass in 0..2 {
            for bk in basis.iter() {
                let c = block.dot(bk, &w);
                for (slot, x) in w.iter_mut().zip(bk.iter()) {
                    *slot -= c * x;
                }
            }
        }
        null.project(block, &mut w);
        block.check_finite("Lanczos vector", &w)?;
        worst_leak = worst_leak.max(null.leak(block, &w));

        let b = block.norm(&w);
        beta.push(b);
        let scale = alpha.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
        if b <= 1e-13 * scale {
            breakdown_at = Some(j + 1);
            break;
        }
        if alpha.len() == m_max {
            break;
        }
        for slot in w.iter_mut() {
            *slot /= b;
        }
        basis.push(w);
    }

    let m = alpha.len();
    let (full_vals, full_last) = tridiagonal_ritz(&alpha[..m], &beta[..m]);
    let half = (m / 2).max(1);
    let (half_vals, _) = tridiagonal_ritz(&alpha[..half], &beta[..half]);

    let lambda_max = full_vals[m - 1];
    let lambda_min = full_vals[0];
    let beta_m = beta[m - 1];
    let kappa = if lambda_min > 0.0 {
        lambda_max / lambda_min
    } else {
        f64::INFINITY
    };
    let half_kappa = if half_vals[0] > 0.0 {
        half_vals[half - 1] / half_vals[0]
    } else {
        f64::INFINITY
    };

    Ok(LanczosReport {
        metric: block.metric(),
        unknowns: block.unknown_count(),
        iterations: m,
        lambda_max,
        lambda_min,
        kappa,
        bound_max_rel: beta_m * full_last[m - 1] / lambda_max,
        bound_min_rel: beta_m * full_last[0] / lambda_min.abs().max(f64::MIN_POSITIVE),
        half_lambda_max: half_vals[half - 1],
        half_lambda_min: half_vals[0],
        half_kappa,
        worst_leak,
        breakdown_at,
    })
}

/// Eigenvalues of the Lanczos tridiagonal, ascending, plus `|e_m^T s_i|` for
/// each: the last component of every unit eigenvector, which is what the
/// residual bound needs. `beta` carries the trailing `beta_m` in its last
/// slot; only `beta[..m-1]` enters the matrix.
fn tridiagonal_ritz(alpha: &[f64], beta: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let m = alpha.len();
    let mut t = vec![0.0f64; m * m];
    for (i, a) in alpha.iter().enumerate() {
        t[i * m + i] = *a;
    }
    for i in 0..m.saturating_sub(1) {
        t[i * m + i + 1] = beta[i];
        t[(i + 1) * m + i] = beta[i];
    }
    symmetric_eigen_last_row(&mut t, m)
}

/// Cyclic Jacobi eigenvalues of a dense symmetric matrix (row-major, `n x n`,
/// overwritten), ascending, together with `|V[n-1][i]|` for each eigenvector.
///
/// Only the LAST ROW of the eigenvector matrix is accumulated, because that
/// is all the Lanczos residual bound uses; that makes the accumulation
/// `O(n^2)` instead of `O(n^3)`. Jacobi rather than a shifted QL: it resolves
/// small eigenvalues of a positive definite matrix to high RELATIVE accuracy,
/// which is the entire point when the smallest one is `1e-5` of the largest.
fn symmetric_eigen_last_row(a: &mut [f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut z = vec![0.0f64; n];
    z[n - 1] = 1.0;
    for _sweep in 0..60 {
        let mut off = 0.0f64;
        let mut diag = 0.0f64;
        for i in 0..n {
            diag += a[i * n + i] * a[i * n + i];
            for j in (i + 1)..n {
                off += a[i * n + j] * a[i * n + j];
            }
        }
        if off <= 1e-32 * diag.max(f64::MIN_POSITIVE) {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq.abs() <= 1e-17 * (a[p * n + p].abs() * a[q * n + q].abs()).sqrt() {
                    a[p * n + q] = 0.0;
                    a[q * n + p] = 0.0;
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (kp, kq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * kp - s * kq;
                    a[k * n + q] = s * kp + c * kq;
                }
                for k in 0..n {
                    let (pk, qk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * pk - s * qk;
                    a[q * n + k] = s * pk + c * qk;
                }
                let (zp, zq) = (z[p], z[q]);
                z[p] = c * zp - s * zq;
                z[q] = s * zp + c * zq;
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|x, y| {
        a[x * n + x]
            .partial_cmp(&a[y * n + y])
            .expect("finite Ritz values")
    });
    let vals = order.iter().map(|k| a[k * n + k]).collect();
    let last = order.iter().map(|k| z[*k].abs()).collect();
    (vals, last)
}

// ---------------------------------------------------------------------
// Solve

/// Tolerances and guards of one CG run. Every field is a stated choice; none
/// of them is tuned per resolution.
#[derive(Debug, Clone, Copy)]
pub struct CgOptions {
    /// Stop when `|r|_P <= tol |b|_P`.
    pub tol: f64,
    pub max_iters: usize,
    /// Re-project the iterate and the residual onto `range(B)` every this
    /// many steps. The projection is mathematically redundant - `B` maps
    /// `range(B)` to itself - so this measures roundoff drift rather than
    /// correcting an error.
    pub reproject_every: usize,
    /// Declare stagnation if the residual has not fallen by
    /// `stagnation_factor` over this many steps.
    pub stagnation_window: usize,
    pub stagnation_factor: f64,
}

impl Default for CgOptions {
    fn default() -> Self {
        Self {
            // 1e-8 sits well above the attainable floor: at the worst kappa
            // measured here (1.23e4 at n = 32) the CG residual cannot be
            // driven below about kappa * eps = 2.7e-12, so this reports the
            // asymptotic rate rather than the roundoff plateau.
            tol: 1e-8,
            max_iters: 20_000,
            reproject_every: 25,
            stagnation_window: 200,
            stagnation_factor: 1.05,
        }
    }
}

/// What one deflated CG run did, including the ways it can fail.
#[derive(Debug, Clone)]
pub struct CgReport {
    pub metric: Metric,
    pub iterations: usize,
    pub converged: bool,
    pub stagnated: bool,
    /// `|r|_P / |b|_P` from the recurrence at the last step.
    pub rel_residual: f64,
    /// The same ratio with `r = b - B x` recomputed from scratch. A gap
    /// between this and `rel_residual` is recurrence drift, and is reported
    /// rather than hidden by returning only one of them.
    pub true_rel_residual: f64,
    /// First iteration at which `|r|_P / |b|_P` fell below `1e-2, 1e-4,
    /// 1e-6, 1e-8`. `None` where it never did.
    pub milestones: [Option<usize>; 4],
    /// Worst null leak of the iterate observed just BEFORE a reprojection.
    pub worst_leak: f64,
}

/// Conjugate gradients for `B x = b` in the `P` inner product, with the six
/// null modes projected out of the right-hand side, of every residual and of
/// every iterate.
///
/// `B` is positive semidefinite and the system is consistent exactly when `b`
/// is `P`-orthogonal to `null(B)`, so the projection of `b` is part of the
/// problem statement and not a convenience. Non-finite iterates and
/// non-positive curvature are errors, not warnings.
pub fn deflated_cg(
    block: &ABlock,
    null: &NullSpace,
    rhs: &[f64],
    opts: CgOptions,
) -> Result<(Vec<f64>, CgReport), MantleError> {
    if null.metric() != block.metric() {
        return Err(MantleError::Field {
            field: "CG deflation metric",
            value: format!("{} vs {}", null.metric().label(), block.metric().label()),
            why: "the null space and the operator must live in the same inner product",
        });
    }
    if opts.reproject_every == 0 {
        return Err(MantleError::Field {
            field: "CG reprojection interval",
            value: "0".to_string(),
            why: "the interval is a modulus and zero would divide by zero",
        });
    }
    block.check_finite("CG right-hand side", rhs)?;

    let mut b = rhs.to_vec();
    null.project(block, &mut b);
    let b_norm = block.norm(&b);
    if !b_norm.is_finite() || b_norm <= 0.0 {
        return Err(MantleError::field(
            "CG right-hand side norm",
            b_norm,
            "the right-hand side is entirely inside the null space, so the system is trivial",
        ));
    }

    let mut x = block.zeros();
    let mut r = b.clone();
    let mut p = r.clone();
    let mut rho = block.dot(&r, &r);
    let mut milestones: [Option<usize>; 4] = [None; 4];
    let thresholds = [1e-2f64, 1e-4, 1e-6, 1e-8];
    let mut worst_leak = 0.0f64;
    let mut converged = false;
    let mut stagnated = false;
    // The relative residual starts at exactly 1, so the window reference is
    // 1 and not `|b|_P`: comparing an absolute norm against a relative
    // residual would make the first window's stagnation test unfireable.
    let mut window_ref = 1.0f64;
    let mut window_start = 0usize;
    let mut k = 0usize;
    let mut rel = 1.0f64;

    while k < opts.max_iters {
        k += 1;
        let mut q = block.apply(&p)?;
        null.project(block, &mut q);
        let pq = block.dot(&p, &q);
        if !pq.is_finite() || pq <= 0.0 {
            return Err(MantleError::field(
                "CG curvature <p, B p>_P",
                pq,
                "the operator is not positive definite on the deflated space",
            ));
        }
        let alpha = rho / pq;
        for ((xs, ps), (rs, qs)) in x.iter_mut().zip(p.iter()).zip(r.iter_mut().zip(q.iter())) {
            *xs += alpha * ps;
            *rs -= alpha * qs;
        }
        block.check_finite("CG iterate", &x)?;
        block.check_finite("CG residual", &r)?;
        if k.is_multiple_of(opts.reproject_every) {
            worst_leak = worst_leak.max(null.leak(block, &x));
            null.project(block, &mut x);
            null.project(block, &mut r);
        }
        let rho_new = block.dot(&r, &r);
        rel = rho_new.sqrt() / b_norm;
        for (slot, t) in milestones.iter_mut().zip(thresholds.iter()) {
            if slot.is_none() && rel <= *t {
                *slot = Some(k);
            }
        }
        if rel <= opts.tol {
            converged = true;
            break;
        }
        if k - window_start >= opts.stagnation_window {
            if window_ref / rel.max(f64::MIN_POSITIVE) < opts.stagnation_factor {
                stagnated = true;
                break;
            }
            window_ref = rel;
            window_start = k;
        }
        let beta = rho_new / rho;
        rho = rho_new;
        for (ps, rs) in p.iter_mut().zip(r.iter()) {
            *ps = rs + beta * *ps;
        }
    }

    worst_leak = worst_leak.max(null.leak(block, &x));
    null.project(block, &mut x);
    let mut true_r = block.apply(&x)?;
    for (slot, bs) in true_r.iter_mut().zip(b.iter()) {
        *slot = bs - *slot;
    }
    null.project(block, &mut true_r);
    let true_rel = block.norm(&true_r) / b_norm;

    Ok((
        x,
        CgReport {
            metric: block.metric(),
            iterations: k,
            converged,
            stagnated,
            rel_residual: rel,
            true_rel_residual: true_rel,
            milestones,
            worst_leak,
        },
    ))
}

/// Tolerances of the inverse-iteration estimate of `lambda_min`.
#[derive(Debug, Clone, Copy)]
pub struct InverseIterationOptions {
    /// Relative residual for the inner CG solve. Loose is fine: an inexact
    /// solve slows the OUTER convergence but cannot corrupt the answer,
    /// because the Rayleigh quotient and its residual are both recomputed
    /// with the exact operator on the vector actually returned.
    pub cg_tol: f64,
    pub cg_max_iters: usize,
    pub max_iterations: usize,
    /// Stop when `|B y - theta y|_P <= tol * theta`.
    pub tol: f64,
}

impl Default for InverseIterationOptions {
    fn default() -> Self {
        Self {
            cg_tol: 1e-8,
            cg_max_iters: 20_000,
            max_iterations: 30,
            tol: 1e-6,
        }
    }
}

/// The smallest non-null eigenvalue, with a bracket that does not depend on
/// believing the iteration converged.
#[derive(Debug, Clone)]
pub struct SmallestEigenReport {
    pub lambda_min: f64,
    /// `|B y - theta y|_P / theta`.
    ///
    /// **What this does and does not bound.** For a self-adjoint operator
    /// SOME eigenvalue lies within `theta * residual_rel` of `theta`, and a
    /// Rayleigh quotient never falls below `lambda_min`. Together those give
    /// `lambda_min <= theta` - rigorously - and they do NOT give
    /// `lambda_min >= theta (1 - residual_rel)`, because the eigenvalue the
    /// residual localises need not be the smallest one. Measured
    /// counterexample in this build: at `diag(K)`, `n = 16` the iteration
    /// settled on `2.9103e-3` with a residual of `5.3e-5` while a Lanczos
    /// Ritz value - which is also an upper bound on `lambda_min` - sat at
    /// `2.8840e-3` BELOW it. The bottom of this spectrum is clustered at
    /// about the percent level and the iteration locked onto the wrong
    /// member of the cluster.
    ///
    /// So `lambda_min` is bounded ABOVE by every estimate produced here and
    /// bounded below by none of them, and `kappa` built from any of them is a
    /// LOWER bound on the true condition number.
    pub residual_rel: f64,
    pub inverse_iterations: usize,
    pub cg_iterations_total: usize,
    pub converged: bool,
    pub leak: f64,
    /// The Rayleigh quotient after each inverse iteration.
    pub history: Vec<f64>,
}

/// `lambda_min` of the deflated block by inverse iteration, each solve done
/// by [`deflated_cg`].
///
/// Lanczos is the wrong tool for this end of the spectrum here and the
/// measurement says so rather than the theory: a 600-step sweep at `n = 16`
/// left `lambda_min` bracketed only to `+-20%`, because the lowest
/// eigenvector's share of any admissible start vector falls like
/// `1/sqrt(N)`. Inverse iteration contracts by `(lambda_1/lambda_2)^k`, a
/// ratio that is a property of the CONTINUUM problem and therefore does not
/// degrade with refinement, at the cost of one CG solve per step.
///
/// Both estimators are reported side by side in the findings. Both are upper
/// bounds on `lambda_min`, so the SMALLER of the two is the better bound
/// whichever has the tidier residual - and where they disagree, that
/// disagreement is itself the measurement that one of them locked onto the
/// wrong end of a cluster.
pub fn smallest_eigenvalue(
    block: &ABlock,
    null: &NullSpace,
    start: &[f64],
    opts: InverseIterationOptions,
) -> Result<SmallestEigenReport, MantleError> {
    if null.metric() != block.metric() {
        return Err(MantleError::Field {
            field: "inverse iteration metric",
            value: format!("{} vs {}", null.metric().label(), block.metric().label()),
            why: "the null space and the operator must live in the same inner product",
        });
    }
    let mut w = start.to_vec();
    block.check_finite("inverse iteration start vector", &w)?;
    null.project(block, &mut w);
    let n0 = block.norm(&w);
    if !n0.is_finite() || n0 <= 0.0 {
        return Err(MantleError::field(
            "inverse iteration start vector norm",
            n0,
            "the start vector is entirely inside the null space",
        ));
    }
    for slot in w.iter_mut() {
        *slot /= n0;
    }

    let mut history = Vec::with_capacity(opts.max_iterations);
    let mut theta = f64::NAN;
    let mut residual_rel = f64::INFINITY;
    let mut cg_total = 0usize;
    let mut used = 0usize;
    let mut converged = false;

    for it in 1..=opts.max_iterations {
        used = it;
        let cg_opts = CgOptions {
            tol: opts.cg_tol,
            max_iters: opts.cg_max_iters,
            ..CgOptions::default()
        };
        let (mut y, rep) = deflated_cg(block, null, &w, cg_opts)?;
        cg_total += rep.iterations;
        if !rep.converged {
            return Err(MantleError::field(
                "inverse iteration inner solve",
                rep.rel_residual,
                "the inner CG did not converge, so the eigenvalue estimate would be a guess",
            ));
        }
        null.project(block, &mut y);
        let ny = block.norm(&y);
        if !ny.is_finite() || ny <= 0.0 {
            return Err(MantleError::field(
                "inverse iteration vector norm",
                ny,
                "the inverse iterate collapsed into the null space",
            ));
        }
        for slot in y.iter_mut() {
            *slot /= ny;
        }
        let by = block.apply(&y)?;
        theta = block.dot(&y, &by);
        if !theta.is_finite() || theta <= 0.0 {
            return Err(MantleError::field(
                "inverse iteration Rayleigh quotient",
                theta,
                "the deflated operator is not positive definite on this vector",
            ));
        }
        let r: Vec<f64> = by
            .iter()
            .zip(y.iter())
            .map(|(b, v)| b - theta * v)
            .collect();
        residual_rel = block.norm(&r) / theta;
        history.push(theta);
        w = y;
        if residual_rel <= opts.tol {
            converged = true;
            break;
        }
    }

    Ok(SmallestEigenReport {
        lambda_min: theta,
        residual_rel,
        inverse_iterations: used,
        cg_iterations_total: cg_total,
        converged,
        leak: null.leak(block, &w),
        history,
    })
}

// ---------------------------------------------------------------------
// Fields

/// The EXACT face data of the affine field `v(x) = b + L x`, using the exact
/// facet tensor moments, zero on the boundary faces. Same construction the
/// operator's own design row is, so `E` reproduces `sym(L)` from it exactly.
pub fn affine_face_field(
    topo: &ShellTopology,
    b: [f64; 3],
    l: [[f64; 3]; 3],
) -> Result<Vec<f64>, MantleError> {
    let moments = topo
        .facet_tensor_moments()
        .map_err(|e| MantleError::Shell {
            what: "facet tensor moments",
            detail: e.to_string(),
        })?;
    Ok(topo
        .faces()
        .iter()
        .zip(moments.iter())
        .map(|(face, tf)| {
            if face.neighbor.is_none() {
                return 0.0;
            }
            let mut acc = face
                .area_vector_m2
                .iter()
                .zip(b.iter())
                .map(|(a, bv)| a * bv)
                .sum::<f64>();
            for (d, lrow) in l.iter().enumerate() {
                for (e, lv) in lrow.iter().enumerate() {
                    acc += lv * tf[d][e];
                }
            }
            acc / face.area_m2
        })
        .collect())
}

/// The face-normal component of an arbitrary smooth vector field, sampled at
/// the facet centroid, zero on the boundary faces.
///
/// Point sampling is deliberate and harmless here: this is used to
/// manufacture a member of the DISCRETE space, not to represent a continuum
/// field accurately. Any vector in the discrete space is an admissible exact
/// solution for a self-consistency check.
pub fn sampled_face_field<F>(topo: &ShellTopology, v: F) -> Vec<f64>
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    topo.faces()
        .iter()
        .map(|face| {
            if face.neighbor.is_none() {
                return 0.0;
            }
            let value = v(face.centroid_m);
            value
                .iter()
                .zip(face.normal.iter())
                .map(|(a, b)| a * b)
                .sum()
        })
        .collect()
}

/// The Lanczos start vector: a smooth field plus a quarter of a random one,
/// each normalised first.
///
/// Neither half alone is adequate and the reason is measured, not assumed. A
/// random vector has `O(1/sqrt(N))` overlap with the lowest eigenvector, so
/// at `n = 16` a 400-step random-start sweep bracketed `lambda_min` only to
/// `+-47%` (and therefore underestimated `kappa`); a purely smooth vector can
/// be orthogonal to a mode by symmetry and would have no generic
/// high-frequency content for `lambda_max`. The mixture has `O(1)` overlap
/// with the smooth end and generic overlap with everything else.
pub fn lanczos_start(topo: &ShellTopology, block: &ABlock, seed: u64) -> Vec<f64> {
    let k = std::f64::consts::PI / topo.spec().r_top_m();
    let mut smooth = sampled_face_field(topo, |x| {
        [
            (k * x[2]).sin() + 0.6 * (k * x[1]).cos(),
            (k * x[0]).sin() + 0.6 * (k * x[2]).cos(),
            (k * x[1]).sin() + 0.6 * (k * x[0]).cos(),
        ]
    });
    let mut noise = pseudo_random_field(block, seed);
    let ns = block.norm(&smooth);
    let nn = block.norm(&noise);
    for slot in smooth.iter_mut() {
        *slot /= ns;
    }
    for slot in noise.iter_mut() {
        *slot /= nn;
    }
    for (slot, x) in smooth.iter_mut().zip(noise.iter()) {
        *slot += 0.25 * x;
    }
    smooth
}

/// Deterministic pseudo-random face field in `[-1, 1)`, zero on the boundary
/// faces. A named generator, so "a random right-hand side" is reproducible.
pub fn pseudo_random_field(block: &ABlock, seed: u64) -> Vec<f64> {
    let mut state = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    let mut out = vec![0.0f64; block.face_count()];
    for (slot, w) in out.iter_mut().zip(block.weight.iter()) {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        if *w > 0.0 {
            *slot = ((state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_grid::shell::ShellSpec;

    const R_CMB: f64 = 3.480e6;
    const R_TOP: f64 = 6.271e6;
    /// `eta = 1`, not the mantle's `1e21`. `A` is linear in `eta`, so every
    /// ratio reported here is invariant under it
    /// ([`the_conditioning_does_not_depend_on_the_viscosity_scale`] measures
    /// that rather than asserting it), and unit viscosity keeps the printed
    /// eigenvalues readable.
    const ETA: f64 = 1.0;

    fn mini(n: usize, layers: usize) -> ShellTopology {
        let spec = ShellSpec::new(n, layers, R_CMB, R_TOP).expect("valid shell");
        ShellTopology::build(&spec).expect("topology builds")
    }

    fn op(n: usize, layers: usize) -> (ShellTopology, ViscousOperator) {
        let t = mini(n, layers);
        let a = ViscousOperator::isoviscous(&t, ETA).expect("operator builds");
        (t, a)
    }

    /// The mesh scale of one rung: the cube root of the MEAN polyhedral cell
    /// volume. One number, measured from the geometry the operator actually
    /// uses, rather than `1/n` - the ladder is not exactly self-similar
    /// (`layers` is an integer) and a fitted power against a nominal `1/n`
    /// would hide that.
    fn mesh_scale(t: &ShellTopology) -> f64 {
        let v = t
            .polyhedral_cell_volumes_m3()
            .expect("polyhedral cell volumes");
        (v.iter().sum::<f64>() / v.len() as f64).cbrt()
    }

    /// Least-squares slope of `log(y)` against `log(x)`.
    fn slope(xs: &[f64], ys: &[f64]) -> f64 {
        let lx: Vec<f64> = xs.iter().map(|v| v.ln()).collect();
        let ly: Vec<f64> = ys.iter().map(|v| v.ln()).collect();
        let mx = lx.iter().sum::<f64>() / lx.len() as f64;
        let my = ly.iter().sum::<f64>() / ly.len() as f64;
        let num: f64 = lx
            .iter()
            .zip(ly.iter())
            .map(|(a, b)| (a - mx) * (b - my))
            .sum();
        let den: f64 = lx.iter().map(|a| (a - mx) * (a - mx)).sum();
        num / den
    }

    /// `(n, layers, Lanczos steps)`. The default two rungs run inside the
    /// normal debug suite; `MC0_ABLOCK_LADDER=full` adds `n = 8` and
    /// `n = 16`, which are release-mode runs and are where the evidence
    /// file's refinement powers come from. The coverage is identical at every
    /// rung; only the top rungs are opt-in.
    ///
    /// `(4, 3)`, `(8, 6)`, `(16, 12)` is an exactly self-similar factor-two
    /// ladder. `(2, 2)` is not - `layers = 1.5` does not exist - so it is
    /// reported and is excluded from the fitted powers, which is stated in
    /// the printed table rather than left to be noticed.
    fn ablock_ladder() -> Vec<(usize, usize, usize)> {
        let mut v = vec![(2usize, 2usize, 110usize), (4, 3, 150)];
        if matches!(std::env::var("MC0_ABLOCK_LADDER").as_deref(), Ok("full")) {
            v.push((8, 6, 300));
            v.push((16, 12, 400));
            // n = 32 gets a SHORT sweep. A fully reorthogonalised Lanczos
            // basis is m vectors of 448,512 doubles, so m = 600 would be
            // 2 GB; 80 steps is 290 MB and is far more than lambda_max
            // needs (its Ritz bracket is printed, so this is checked and not
            // assumed), while lambda_min comes from inverse iteration, whose
            // cost is a handful of CG solves and whose storage is a handful
            // of vectors.
            v.push((32, 24, 80));
        }
        v
    }

    // -----------------------------------------------------------------
    // 1. The pieces, before anything is measured with them

    #[test]
    fn the_stiffness_diagonal_matches_the_operator_it_claims_to_be() {
        // diag(K) is read off the stored strain stencil rather than probed
        // with n_faces operator applications. If that shortcut is wrong, the
        // Jacobi runs below are measuring a different matrix and nothing
        // else in this module would notice.
        let (t, a) = op(2, 2);
        let d = stiffness_diagonal(&a);
        let unknowns = ViscousOperator::unknown_faces(&t);
        let mut probe = vec![0.0f64; a.strain().face_count()];
        let mut worst = 0.0f64;
        for &f in unknowns.iter() {
            probe[f as usize] = 1.0;
            let col = a.stiffness_apply(&probe).expect("stiffness apply");
            probe[f as usize] = 0.0;
            let exact = col[f as usize];
            assert!(exact > 0.0, "face {f}: K_ff = {exact:e} is not positive");
            worst = worst.max(((d[f as usize] - exact) / exact).abs());
        }
        println!(
            "diag(K) against probed K_ff over {} interior faces: worst relative difference {worst:.3e}",
            unknowns.len()
        );
        assert!(
            worst < 1e-12,
            "diag(K) disagrees with the operator by {worst:e}"
        );
    }

    #[test]
    fn the_six_modes_are_null_and_their_projection_is_idempotent() {
        for (n, layers) in [(2usize, 2usize), (4, 3)] {
            let (t, a) = op(n, layers);
            for metric in [Metric::FaceMass, Metric::StiffnessDiagonal] {
                let block = ABlock::new(&t, &a, metric).expect("block builds");
                let null = NullSpace::build(&t, &block).expect("null space builds");

                // (a) orthonormality of the deflation basis in ITS metric.
                let mut worst_gram = 0.0f64;
                for i in 0..null.dim() {
                    for j in 0..null.dim() {
                        let g = block.dot(null.mode(i), null.mode(j));
                        let want = if i == j { 1.0 } else { 0.0 };
                        worst_gram = worst_gram.max((g - want).abs());
                    }
                }

                // (b) the modes really are null: |B n_j|_P against the size
                // of the operator, which is what "null" has to mean.
                let scale = {
                    let mut probe = pseudo_random_field(&block, 7);
                    null.project(&block, &mut probe);
                    let nb = block.norm(&probe);
                    let bp = block.apply(&probe).expect("apply");
                    block.norm(&bp) / nb
                };
                let mut worst_null = 0.0f64;
                for j in 0..null.dim() {
                    let bn = block.apply(null.mode(j)).expect("apply");
                    worst_null = worst_null.max(block.norm(&bn) / scale);
                }

                // (c) idempotence, and the leak of a projected vector.
                let mut v = pseudo_random_field(&block, 11);
                null.project(&block, &mut v);
                let once = v.clone();
                null.project(&block, &mut v);
                let drift = {
                    let d: Vec<f64> = v.iter().zip(once.iter()).map(|(a, b)| a - b).collect();
                    block.norm(&d) / block.norm(&once)
                };
                let leak = null.leak(&block, &once);

                println!(
                    "n={n} metric={:<7} pivot={:.3e} gram-I={:.3e} |B n_j|/|B|={:.3e} \
                     idempotence drift={:.3e} leak={:.3e}",
                    metric.label(),
                    null.smallest_pivot(),
                    worst_gram,
                    worst_null,
                    drift,
                    leak
                );
                assert!(
                    worst_gram < 1e-12,
                    "deflation basis is not orthonormal in {}: worst |Gram - I| = {worst_gram:e}",
                    metric.label()
                );
                assert!(
                    worst_null < 1e-10,
                    "a deflated mode is not null in {}: |B n_j|/|B| = {worst_null:e}",
                    metric.label()
                );
                assert!(
                    drift < 1e-12,
                    "the projection is not idempotent in {}: drift {drift:e}",
                    metric.label()
                );
                assert!(
                    leak < 1e-12,
                    "a projected vector still has a null component in {}: {leak:e}",
                    metric.label()
                );
            }
        }
    }

    #[test]
    fn the_operator_is_self_adjoint_and_its_quadratic_form_is_the_dissipation() {
        let (t, a) = op(4, 3);
        for metric in [Metric::FaceMass, Metric::StiffnessDiagonal] {
            let block = ABlock::new(&t, &a, metric).expect("block builds");
            let u = pseudo_random_field(&block, 3);
            let v = pseudo_random_field(&block, 5);
            let bu = block.apply(&u).expect("apply");
            let bv = block.apply(&v).expect("apply");
            let left = block.dot(&u, &bv);
            let right = block.dot(&bu, &v);
            let sym = (left - right).abs() / left.abs().max(right.abs());
            let quad = block.dot(&u, &bu);
            let energy = block.energy_form(&u).expect("energy");
            let agree = (quad - energy).abs() / energy.abs();
            println!(
                "metric={:<7} <u,Bv>_P vs <Bu,v>_P: {sym:.3e} | <u,Bu>_P vs Phi(u): {agree:.3e}",
                metric.label()
            );
            assert!(
                sym < 1e-12,
                "B is not self-adjoint in {}: {sym:e}",
                metric.label()
            );
            assert!(
                agree < 1e-12,
                "the quadratic form is not the dissipation in {}: {agree:e}",
                metric.label()
            );
        }
    }

    // -----------------------------------------------------------------
    // 2. The Lanczos sweep, checked against a dense eigendecomposition

    /// Dense eigenvalues of `B = P^-1 K` at a resolution small enough to
    /// assemble it: the similarity transform `P^-1/2 K P^-1/2` is symmetric
    /// with the same spectrum, so a Jacobi eigendecomposition of it is the
    /// exact answer the matrix-free sweep must reproduce.
    fn dense_spectrum(t: &ShellTopology, a: &ViscousOperator, metric: Metric) -> Vec<f64> {
        let unknowns = ViscousOperator::unknown_faces(t);
        let nu = unknowns.len();
        let s = a.dense_stiffness(&unknowns).expect("dense stiffness");
        let diag = match metric {
            Metric::FaceMass => a.face_mass().to_vec(),
            Metric::StiffnessDiagonal => stiffness_diagonal(a),
        };
        let inv_sqrt: Vec<f64> = unknowns
            .iter()
            .map(|&f| 1.0 / diag[f as usize].sqrt())
            .collect();
        let mut b = vec![0.0f64; nu * nu];
        for (i, row) in s.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                b[i * nu + j] = 0.5 * (v + s[j][i]) * inv_sqrt[i] * inv_sqrt[j];
            }
        }
        let (vals, _) = symmetric_eigen_last_row(&mut b, nu);
        vals
    }

    #[test]
    fn the_deflated_lanczos_spectrum_matches_a_dense_eigendecomposition() {
        // The whole conditioning table below rests on this sweep. At n = 2
        // the matrix fits, so the sweep can be checked against the answer
        // rather than against a plausibility argument.
        let (t, a) = op(2, 2);
        for metric in [Metric::FaceMass, Metric::StiffnessDiagonal] {
            let block = ABlock::new(&t, &a, metric).expect("block builds");
            let null = NullSpace::build(&t, &block).expect("null space builds");
            let exact = dense_spectrum(&t, &a, metric);
            let nu = exact.len();
            let dense_max = exact[nu - 1];
            let dense_min = exact[NULL_DIM];
            let dense_kappa = dense_max / dense_min;
            let start = lanczos_start(&t, &block, 1);
            let r = lanczos_spectrum(&block, &null, 110, &start).expect("lanczos runs");

            let err_max = ((r.lambda_max - dense_max) / dense_max).abs();
            let err_min = ((r.lambda_min - dense_min) / dense_min).abs();
            println!(
                "n=2 metric={:<7} unknowns={nu} dense: lmax={dense_max:.6e} l7={dense_min:.6e} \
                 kappa={dense_kappa:.4e} | null eigenvalues 0..5 max |l|/lmax = {:.3e}",
                metric.label(),
                exact[..NULL_DIM]
                    .iter()
                    .fold(0.0f64, |acc, v| acc.max(v.abs()))
                    / dense_max
            );
            println!(
                "n=2 metric={:<7} lanczos m={} lmax={:.6e} lmin={:.6e} kappa={:.4e} \
                 rel err (max,min)=({err_max:.2e},{err_min:.2e}) bounds=({:.2e},{:.2e}) leak={:.2e}",
                metric.label(),
                r.iterations,
                r.lambda_max,
                r.lambda_min,
                r.kappa,
                r.bound_max_rel,
                r.bound_min_rel,
                r.worst_leak
            );
            assert!(
                err_max < 1e-9,
                "{}: Lanczos lambda_max is {err_max:e} off the dense value",
                metric.label()
            );
            assert!(
                err_min < 1e-6,
                "{}: Lanczos lambda_min is {err_min:e} off the dense value",
                metric.label()
            );
            assert!(
                r.worst_leak < 1e-10,
                "{}: a Lanczos vector leaked {:e} into the null space",
                metric.label(),
                r.worst_leak
            );

            // The second estimator of lambda_min, against the same known
            // answer. It is the one the ladder relies on at the fine end, so
            // it has to be checked where the answer is available.
            let inv =
                smallest_eigenvalue(&block, &null, &start, InverseIterationOptions::default())
                    .expect("inverse iteration runs");
            let err_inv = ((inv.lambda_min - dense_min) / dense_min).abs();
            println!(
                "n=2 metric={:<7} inverse iteration: lmin={:.6e} rel err={err_inv:.2e} \
                 bracket={:.2e} outer={} inner cg={} converged={} leak={:.1e}",
                metric.label(),
                inv.lambda_min,
                inv.residual_rel,
                inv.inverse_iterations,
                inv.cg_iterations_total,
                inv.converged,
                inv.leak
            );
            assert!(
                err_inv < 1e-8,
                "{}: inverse iteration lambda_min is {err_inv:e} off the dense value",
                metric.label()
            );
            assert!(
                inv.lambda_min >= dense_min * (1.0 - 1e-9),
                "{}: a Rayleigh quotient came out BELOW the true lambda_min ({:e} < {dense_min:e}), \
                 which is impossible and means the deflation is leaking",
                metric.label(),
                inv.lambda_min
            );
        }
    }

    // -----------------------------------------------------------------
    // 3. The solve

    #[test]
    fn deflated_cg_solves_the_system_and_the_iterates_stay_in_the_complement() {
        for (n, layers) in [(2usize, 2usize), (4, 3)] {
            let (t, a) = op(n, layers);
            for metric in [Metric::FaceMass, Metric::StiffnessDiagonal] {
                let block = ABlock::new(&t, &a, metric).expect("block builds");
                let null = NullSpace::build(&t, &block).expect("null space builds");
                let rhs = pseudo_random_field(&block, 17);
                let opts = CgOptions {
                    max_iters: 5000,
                    ..CgOptions::default()
                };
                let (x, rep) = deflated_cg(&block, &null, &rhs, opts).expect("cg runs");
                let final_leak = null.leak(&block, &x);
                println!(
                    "n={n} metric={:<7} cg iters={} converged={} stagnated={} \
                     recursive res={:.3e} true res={:.3e} worst leak={:.3e} final leak={:.3e} \
                     milestones={:?}",
                    metric.label(),
                    rep.iterations,
                    rep.converged,
                    rep.stagnated,
                    rep.rel_residual,
                    rep.true_rel_residual,
                    rep.worst_leak,
                    final_leak,
                    rep.milestones
                );
                assert!(
                    rep.converged,
                    "n={n} {}: CG did not reach 1e-8 in {} iterations (residual {:e}, stagnated={})",
                    metric.label(),
                    rep.iterations,
                    rep.rel_residual,
                    rep.stagnated
                );
                assert!(
                    rep.true_rel_residual < 10.0 * opts.tol,
                    "n={n} {}: the recomputed residual {:e} does not match the recursive one {:e}",
                    metric.label(),
                    rep.true_rel_residual,
                    rep.rel_residual
                );
                assert!(
                    final_leak < 1e-10,
                    "n={n} {}: the solution has a null component {final_leak:e}",
                    metric.label()
                );
            }
        }
    }

    #[test]
    fn the_solver_refuses_a_non_finite_or_degenerate_problem_loudly() {
        let (t, a) = op(2, 2);
        let block = ABlock::new(&t, &a, Metric::FaceMass).expect("block builds");
        let null = NullSpace::build(&t, &block).expect("null space builds");

        // (a) a NaN in the right-hand side is named, not solved around.
        let mut bad = pseudo_random_field(&block, 23);
        let victim = ViscousOperator::unknown_faces(&t)[5] as usize;
        bad[victim] = f64::NAN;
        let err = deflated_cg(&block, &null, &bad, CgOptions::default())
            .expect_err("a NaN right-hand side must be an error");
        println!("NaN right-hand side -> {err}");
        assert!(matches!(err, MantleError::NonFinite { index, .. } if index == victim));

        // (b) a right-hand side that is entirely inside the null space is a
        // degenerate problem, not a zero-iteration success.
        let mut trivial = block.zeros();
        for (slot, x) in trivial.iter_mut().zip(null.mode(0).iter()) {
            *slot = *x;
        }
        let err = deflated_cg(&block, &null, &trivial, CgOptions::default())
            .expect_err("a null right-hand side must be an error");
        println!("null right-hand side -> {err}");

        // (c) a null space built in the wrong inner product is refused
        // rather than silently used.
        let jac = ABlock::new(&t, &a, Metric::StiffnessDiagonal).expect("block builds");
        let err = deflated_cg(
            &jac,
            &null,
            &pseudo_random_field(&jac, 29),
            CgOptions::default(),
        )
        .expect_err("a metric mismatch must be an error");
        println!("metric mismatch -> {err}");
    }

    #[test]
    fn the_wrong_inner_product_leaves_a_persistent_null_component() {
        // [NEG-CONTROL] The failure this reproduces is not a crash. Deflate
        // with coefficients taken in M_u while the operator and the residual
        // live in diag(K), and the projection is simply not a projection:
        // what is left keeps a null component that no residual test can see,
        // because the residual is orthogonal to it in the wrong product.
        let (t, a) = op(4, 3);
        let right = ABlock::new(&t, &a, Metric::StiffnessDiagonal).expect("block builds");
        let wrong = ABlock::new(&t, &a, Metric::FaceMass).expect("block builds");
        let null_right = NullSpace::build(&t, &right).expect("null space builds");
        let null_wrong = NullSpace::build(&t, &wrong).expect("null space builds");

        let v = pseudo_random_field(&right, 31);

        let mut good = v.clone();
        null_right.project(&right, &mut good);
        let good_leak = null_right.leak(&right, &good);

        // The mismatch: subtract the M_u-orthonormal modes with M_u
        // coefficients, then ask what is left in the product the solve
        // actually runs in.
        let mut mixed = v.clone();
        null_wrong.project(&wrong, &mut mixed);
        let mixed_leak = null_right.leak(&right, &mixed);

        println!(
            "consistent projection leak = {good_leak:.3e}; \
             M_u coefficients used in a diag(K) solve leak = {mixed_leak:.3e} \
             (ratio {:.3e})",
            mixed_leak / good_leak.max(f64::MIN_POSITIVE)
        );
        assert!(
            mixed_leak > 1e-3,
            "the negative control did not reproduce the mismatch: leak {mixed_leak:e}"
        );
        assert!(
            good_leak < 1e-12,
            "the consistent projection is not clean: leak {good_leak:e}"
        );
    }

    #[test]
    fn the_conditioning_does_not_depend_on_the_viscosity_scale() {
        // Every ratio in the tables below is claimed to be independent of
        // eta. That is a claim about the code, not only about the algebra:
        // eta enters as a per-cell weight and could in principle interact
        // with a tolerance.
        let t = mini(2, 2);
        let mut kappas = Vec::new();
        for eta in [1.0f64, 1.0e21] {
            let a = ViscousOperator::isoviscous(&t, eta).expect("operator builds");
            let block = ABlock::new(&t, &a, Metric::FaceMass).expect("block builds");
            let null = NullSpace::build(&t, &block).expect("null space builds");
            let start = lanczos_start(&t, &block, 1);
            let r = lanczos_spectrum(&block, &null, 110, &start).expect("lanczos runs");
            println!(
                "eta={eta:.1e}: lmax={:.6e} lmin={:.6e} kappa={:.6e}",
                r.lambda_max, r.lambda_min, r.kappa
            );
            kappas.push(r.kappa);
        }
        let drift = (kappas[0] - kappas[1]).abs() / kappas[0];
        println!("kappa drift across 21 decades of eta: {drift:.3e}");
        assert!(drift < 1e-10, "kappa moved with eta by {drift:e}");
    }

    // -----------------------------------------------------------------
    // 4. Manufactured solution - a SELF-CONSISTENCY check, not physics

    #[test]
    fn the_manufactured_solution_is_recovered_modulo_the_null_space() {
        // This verifies that the solver inverts the operator this repository
        // built. It says NOTHING about whether that operator is the right
        // physics: u_exact is manufactured in the discrete space and the
        // right-hand side is A u_exact, so the continuum never enters. The
        // analytic comparison is MC-Q1 (benchmark_thieulot).
        for (n, layers) in [(2usize, 2usize), (4, 3)] {
            let (t, a) = op(n, layers);
            let block = ABlock::new(&t, &a, Metric::FaceMass).expect("block builds");
            let null = NullSpace::build(&t, &block).expect("null space builds");

            // (i) an affine symmetric traceless strain, exact facet data.
            let s = 1.0 / R_TOP;
            let l = [
                [1.0 * s, 0.5 * s, 0.25 * s],
                [0.5 * s, -0.6 * s, 0.3 * s],
                [0.25 * s, 0.3 * s, -0.4 * s],
            ];
            let affine = affine_face_field(&t, [0.0; 3], l).expect("affine face field");

            // (ii) a smooth non-affine member of the discrete space.
            let k = std::f64::consts::PI / R_TOP;
            let trig = sampled_face_field(&t, |x| {
                [(k * x[2]).sin(), (k * x[0]).sin(), (k * x[1]).sin()]
            });

            for (name, u_raw) in [("affine-strain", affine), ("trig-sampled", trig)] {
                let mut u_exact = u_raw.clone();
                null.project(&block, &mut u_exact);
                let rhs = block.apply(&u_exact).expect("apply");
                // A tighter tolerance than the ladder's 1e-8: the SOLUTION
                // error of a CG run stopped at a residual tol is bounded by
                // about kappa * tol, so a solution-error check has to buy
                // itself the headroom rather than report kappa as an error.
                let opts = CgOptions {
                    tol: 1e-12,
                    max_iters: 5000,
                    ..CgOptions::default()
                };
                let (x, rep) = deflated_cg(&block, &null, &rhs, opts).expect("cg runs");
                let diff: Vec<f64> = x.iter().zip(u_exact.iter()).map(|(a, b)| a - b).collect();
                let err = block.norm(&diff) / block.norm(&u_exact);
                println!(
                    "n={n} {name:<14} cg iters={:<5} residual={:.2e} \
                     |x - u|_Mu / |u|_Mu = {err:.3e} (null-projected u)",
                    rep.iterations, rep.true_rel_residual
                );
                assert!(
                    rep.converged,
                    "n={n} {name}: CG did not converge ({} iterations, residual {:e})",
                    rep.iterations, rep.rel_residual
                );
                assert!(
                    err < 1e-6,
                    "n={n} {name}: manufactured field not recovered, error {err:e}"
                );
            }
        }
    }

    // -----------------------------------------------------------------
    // 5. The measurement this packet exists for

    #[test]
    fn the_spectral_conditioning_and_the_cg_iteration_count_across_refinement() {
        let ladder = ablock_ladder();
        for metric in [Metric::FaceMass, Metric::StiffnessDiagonal] {
            let mut hs = Vec::new();
            let mut ns = Vec::new();
            let mut iters8 = Vec::new();
            let mut spec_h = Vec::new();
            let mut spec_n: Vec<usize> = Vec::new();
            let mut lmax = Vec::new();
            let mut lmin = Vec::new();
            let mut lmin_lanczos: Vec<f64> = Vec::new();
            // CG counts restricted to the rungs that also produced a
            // spectrum, so that `kappas[i]` and this share an index. `iters8`
            // is indexed by CG rung and the two lists differ the moment a
            // rung runs CG without a sweep.
            let mut spec_iters: Vec<f64> = Vec::new();
            let mut kappas = Vec::new();
            let mut bnd_min: Vec<f64> = Vec::new();

            println!(
                "\n=== metric {} : deflated spectrum and unpreconditioned CG ===",
                metric.label()
            );
            println!(
                "{:>3} {:>3} {:>8} {:>10} {:>6} {:>11} {:>11} {:>10} {:>9} {:>9} {:>10} {:>8}",
                "n",
                "lay",
                "unknown",
                "h (m)",
                "m",
                "lambda_max",
                "lambda_min",
                "kappa",
                "bnd_max",
                "bnd_min",
                "half_kappa",
                "leak"
            );
            for &(n, layers, steps) in ladder.iter() {
                let (t, a) = op(n, layers);
                let block = ABlock::new(&t, &a, metric).expect("block builds");
                let null = NullSpace::build(&t, &block).expect("null space builds");
                let h = mesh_scale(&t);
                let spectrum = if steps > 0 {
                    let start = lanczos_start(&t, &block, 1);
                    let r = lanczos_spectrum(&block, &null, steps, &start).expect("lanczos runs");
                    println!(
                        "{n:>3} {layers:>3} {:>8} {h:>10.4e} {:>6} {:>11.4e} {:>11.4e} {:>10.4e} \
                         {:>9.2e} {:>9.2e} {:>10.4e} {:>8.1e}",
                        r.unknowns,
                        r.iterations,
                        r.lambda_max,
                        r.lambda_min,
                        r.kappa,
                        r.bound_max_rel,
                        r.bound_min_rel,
                        r.half_kappa,
                        r.worst_leak
                    );
                    assert!(
                        r.breakdown_at.is_none(),
                        "n={n}: the Lanczos recurrence broke down at step {:?}",
                        r.breakdown_at
                    );
                    assert!(
                        r.lambda_min > 0.0,
                        "n={n}: the smallest deflated Ritz value is {:e}; the deflation is incomplete",
                        r.lambda_min
                    );
                    let inv = smallest_eigenvalue(
                        &block,
                        &null,
                        &start,
                        InverseIterationOptions {
                            // Loose inner solve on purpose: the outer
                            // residual is recomputed exactly, so this only
                            // costs outer steps.
                            cg_tol: 1e-6,
                            cg_max_iters: 20_000,
                            max_iterations: 30,
                            tol: 1e-6,
                        },
                    )
                    .expect("inverse iteration runs");
                    // BOTH estimators are upper bounds on lambda_min - a
                    // Ritz value and a Rayleigh quotient both lie inside the
                    // spectrum - so the SMALLER one is the better bound and
                    // is the one kappa is built from. Taking the inverse
                    // iteration's answer unconditionally would have thrown
                    // away the better number at diag(K), n = 4.
                    let best = inv.lambda_min.min(r.lambda_min);
                    println!(
                        "      lambda_min by inverse iteration: {:.4e} (bracket {:.1e}, \
                         outer {}, inner cg {}, converged {}) | lanczos said {:.4e} \
                         (bracket {:.1e}) | best = {:.4e} from {} -> kappa = {:.4e}",
                        inv.lambda_min,
                        inv.residual_rel,
                        inv.inverse_iterations,
                        inv.cg_iterations_total,
                        inv.converged,
                        r.lambda_min,
                        r.bound_min_rel,
                        best,
                        if inv.lambda_min <= r.lambda_min {
                            "inverse iteration"
                        } else {
                            "lanczos"
                        },
                        r.lambda_max / best
                    );
                    Some((r, inv))
                } else {
                    println!(
                        "{n:>3} {layers:>3} {:>8} {h:>10.4e} {:>6} {:>11} {:>11} {:>10} \
                         {:>9} {:>9} {:>10} {:>8}",
                        block.unknown_count(),
                        "-",
                        "(CG only)",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-"
                    );
                    None
                };

                let rhs = pseudo_random_field(&block, 101);
                let opts = CgOptions {
                    max_iters: 20_000,
                    ..CgOptions::default()
                };
                let (_x, rep) = deflated_cg(&block, &null, &rhs, opts).expect("cg runs");
                println!(
                    "      cg: iters={:<6} converged={} stagnated={} true res={:.2e} \
                     milestones(1e-2,1e-4,1e-6,1e-8)={:?} worst leak={:.1e}",
                    rep.iterations,
                    rep.converged,
                    rep.stagnated,
                    rep.true_rel_residual,
                    rep.milestones,
                    rep.worst_leak
                );

                hs.push(h);
                ns.push(n);
                iters8.push(rep.iterations as f64);
                if let Some((r, inv)) = spectrum {
                    let best = inv.lambda_min.min(r.lambda_min);
                    spec_h.push(h);
                    spec_n.push(n);
                    lmax.push(r.lambda_max);
                    lmin.push(best);
                    lmin_lanczos.push(r.lambda_min);
                    spec_iters.push(rep.iterations as f64);
                    kappas.push(r.lambda_max / best);
                    bnd_min.push(if inv.lambda_min <= r.lambda_min {
                        inv.residual_rel
                    } else {
                        r.bound_min_rel
                    });
                }
            }

            // The fitted powers use the SELF-SIMILAR rungs only (n >= 4).
            // n = 2 is (2, 2) because layers = 1.5 does not exist, so its
            // cell aspect ratio differs from the rest of the ladder and its
            // 48 cells are nowhere near an asymptotic regime; it is reported
            // and excluded, rather than quietly averaged in.
            let spec_keep: Vec<usize> = (0..spec_n.len()).filter(|i| spec_n[*i] >= 4).collect();
            let cg_keep: Vec<usize> = (0..ns.len()).filter(|i| ns[*i] >= 4).collect();
            if spec_keep.len() >= 2 {
                let sel = |v: &[f64]| -> Vec<f64> { spec_keep.iter().map(|i| v[*i]).collect() };
                let hh = sel(&spec_h);
                println!(
                    "observed powers of h over n = {:?}: lambda_max h^{:+.3} \
                     lambda_min h^{:+.3} kappa h^{:+.3}",
                    spec_keep.iter().map(|i| spec_n[*i]).collect::<Vec<_>>(),
                    slope(&hh, &sel(&lmax)),
                    slope(&hh, &sel(&lmin)),
                    slope(&hh, &sel(&kappas))
                );
                let p_kappa = slope(&hh, &sel(&kappas));
                let cgsel = |v: &[f64]| -> Vec<f64> { cg_keep.iter().map(|i| v[*i]).collect() };
                println!(
                    "cg iterations over n = {:?}: h^{:+.3}; sqrt(kappa) predicts h^{:+.3}",
                    cg_keep.iter().map(|i| ns[*i]).collect::<Vec<_>>(),
                    slope(&cgsel(&hs), &cgsel(&iters8)),
                    0.5 * p_kappa
                );
                // The pairwise powers matter more than the fit: the coarse
                // pairs are pre-asymptotic and a single least-squares number
                // averages that away.
                println!(
                    "lambda_min, best of two estimators {:?}; lanczos alone said {:?}",
                    spec_keep
                        .iter()
                        .map(|i| format!("{:.4e}", lmin[*i]))
                        .collect::<Vec<_>>(),
                    spec_keep
                        .iter()
                        .map(|i| format!("{:.4e}", lmin_lanczos[*i]))
                        .collect::<Vec<_>>()
                );
                for i in 1..spec_keep.len() {
                    let (a, b) = (spec_keep[i - 1], spec_keep[i]);
                    let r = (spec_h[b] / spec_h[a]).ln();
                    println!(
                        "   spectrum {} -> {}: lambda_max h^{:+.3}  lambda_min h^{:+.3}  \
                         kappa h^{:+.3}   (lambda_min bracket {:.1e} -> {:.1e})",
                        spec_n[a],
                        spec_n[b],
                        (lmax[b] / lmax[a]).ln() / r,
                        (lmin[b] / lmin[a]).ln() / r,
                        (kappas[b] / kappas[a]).ln() / r,
                        bnd_min[a],
                        bnd_min[b]
                    );
                }
            } else {
                println!(
                    "only {} self-similar spectrum rung(s) in the default ladder: run with \
                     MC0_ABLOCK_LADDER=full for the fitted powers",
                    spec_keep.len()
                );
            }

            // The CG bound is 0.5 sqrt(kappa) ln(2/tol) iterations, so the
            // constant to watch is iters / (0.5 sqrt(kappa) ln(2/tol)) - a
            // bare iters/sqrt(kappa) hides a factor of 19 that is pure
            // tolerance and would read as a disagreement that is not one.
            let logfac = (2.0f64 / 1e-8).ln();
            for (i, k) in kappas.iter().enumerate() {
                println!(
                    "   n={:<3} sqrt(kappa)={:>9.2e}  cg iters={:>7}  \
                     iters/sqrt(kappa)={:>7.3}  iters/(0.5 sqrt(kappa) ln(2/tol))={:.3}",
                    spec_n[i],
                    k.sqrt(),
                    spec_iters[i] as usize,
                    spec_iters[i] / k.sqrt(),
                    spec_iters[i] / (0.5 * k.sqrt() * logfac)
                );
            }
            for i in 1..ns.len() {
                println!(
                    "   cg rung {} -> {}: h x{:.4}, iters x{:.4} (power h^{:+.3})",
                    ns[i - 1],
                    ns[i],
                    hs[i] / hs[i - 1],
                    iters8[i] / iters8[i - 1],
                    (iters8[i] / iters8[i - 1]).ln() / (hs[i] / hs[i - 1]).ln()
                );
            }

            // The verdict this packet exists to reach is a NUMBER, not a
            // pass. Assert only what would invalidate the measurement - a
            // kappa that fails to grow is a finding to report, not a test
            // failure, so it is printed and not asserted.
            for (i, k) in kappas.iter().enumerate() {
                assert!(
                    k.is_finite() && *k > 1.0,
                    "n={}: kappa = {k:e} is not a condition number",
                    spec_n[i]
                );
            }
        }
    }
}
