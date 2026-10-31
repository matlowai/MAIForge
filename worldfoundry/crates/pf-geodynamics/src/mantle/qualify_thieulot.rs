//! MC-Q1: the coupled Stokes solve, driven by Thieulot's analytic body force
//! and compared against his closed-form solution at `m = -1`.
//!
//! **This is the first check in this project against an independent statement
//! of what mantle flow does.** Everything before it - the adjoint identities,
//! the manufactured saddle pair, the conditioning ladders - is
//! self-consistency: the operators inverting operators this repository wrote.
//! Here the right-hand side comes from a published paper and the answer is
//! not ours to choose.
//!
//! # THE HEADLINE, STATED BEFORE THE CODE
//!
//! **Thieulot's solution is not a solution of the boundary-value problem this
//! discretization solves, and no amount of mesh refinement will make it one.**
//! The two problems share a PDE and a body force and differ in the tangential
//! boundary condition:
//!
//! | | Thieulot (2017) | this discretization |
//! |---|---|---|
//! | normal velocity at `R1`, `R2` | `v_r = 0` (eq. 18) | `u_f = 0` on boundary facets - the same |
//! | tangential condition | the ANALYTIC VELOCITY is prescribed (Sect. 5.3, VERIFIED-PRIMARY in the research file) | none imposed, so the variational NATURAL condition applies: zero tangential traction, i.e. FREE SLIP |
//!
//! Thieulot's field does not have zero tangential traction. Its shear stress
//! at a boundary is `sigma_r_theta = (mu/r)(r f' - f - g) sin(theta)`, and at
//! `r = R1` or `R2` where `g = 0` this is `-3 mu alpha sin(theta) / r^3`,
//! which vanishes only if `alpha = 0` - i.e. only for the motionless case
//! this fixture refuses to build. Measured here at the reference parameters:
//! `f' - f/r = 30.2966` at `R1` and `3.7871` at `R2`, against a velocity of
//! order 1. So the traction the benchmark needs at the wall is of the same
//! order as the whole flow, not a small correction.
//!
//! Two further structural facts, both measured in this module rather than
//! assumed:
//!
//! 1. **The benchmark is not meridional.** Eq. (17) gives
//!    `v_phi = f(r) sin(theta) = v_theta`, and that toroidal component carries
//!    **47.77 % of the kinetic energy** (`v_theta` another 47.77 %, `v_r` only
//!    4.46 %). Its projection on a rigid rotation about `z` has amplitude
//!    `Omega_z = 0.556726`, which is **24.03 % of `|u|`** in the `L2` sense.
//!    The hazard note that "Thieulot's solution is meridional, so naively it
//!    should have little rigid-rotation content" is therefore FALSE at the
//!    reference parameters, which is exactly why it was to be checked.
//! 2. **Under free slip the toroidal part collapses to an undetermined rigid
//!    rotation.** `v_phi = f_t(r) sin(theta)` satisfies
//!    `r^2 f_t'' + 2 r f_t' - 2 f_t = 0`, so `f_t = a/r^2 + b r`; zero
//!    tangential traction forces `a = 0`, and `b r` is the rigid rotation,
//!    which is a null mode of the free-slip problem and is fixed at zero by
//!    the preconditioner. Half the benchmark's energy is therefore not merely
//!    hard to compute here - it is not part of the problem being solved.
//!
//! # SO WHAT IS THE PHYSICS CHECK
//!
//! Refusing to compare would leave MC-Q1 unanswered, and tuning a tolerance
//! until the mismatch passed would be worse. [`FreeSlipCase`] is the third
//! option: the **exact closed-form solution of the free-slip problem with
//! Thieulot's own body force**, derived here and verified here. The
//! qualification then asks the question the solver can actually be wrong
//! about:
//!
//! > given the published body force, does our discrete Stokes solve reproduce
//! > the analytic solution of the boundary-value problem it is actually
//! > discretizing?
//!
//! and reports the Thieulot comparison alongside it as the measurement of how
//! far the two boundary-value problems are apart.
//!
//! The construction, for the record. For `l = 1` and constant viscosity the
//! poloidal scalar `P(r)` (with `v_r = (2/r) P cos(theta)`,
//! `v_theta = -(P' + P/r) sin(theta)`) satisfies `D^2 P = source`, with
//! `D = d2/dr2 + (2/r) d/dr - 2/r^2`, so the homogeneous solutions are
//! `P in {r, r^3, 1, r^-2}` - four of them, and four boundary conditions
//! (`v_r = 0` and `sigma_r_theta = 0` at each of `R1`, `R2`). Thieulot's own
//! field is one particular solution (its poloidal scalar is
//! `P_T = r g_T(r) / 2`, and `-(P_T' + P_T/r)` reproduces `f_T` exactly,
//! which is checked in [`tests::the_free_slip_reference_is_a_solution_of_the_same_problem`]).
//! Adding the four homogeneous modes with coefficients fixed by the free-slip
//! conditions gives the free-slip answer to the same `rho`.
//!
//! # NEGATIVE CONTROL: the built-in one is DEGENERATE at this case
//!
//! [`super::benchmark_thieulot::ThieulotCase::density_as_printed_in_eq59`] is
//! the paper's printed, mislabelled density, and it is the right control for
//! `m != -1`. At `m = -1` with `mu0 = 1` - which is exactly
//! [`super::benchmark_thieulot::ThieulotCase::isoviscous_reference`] - it
//! returns **bit-identical values to the corrected density**: eq. (50) has no
//! `r^m` factor to lose (correction 1 does not apply) and `mu0 = 1` makes
//! correction 2 invisible. Asserted, not asserted-about, in
//! [`tests::the_built_in_negative_control_is_degenerate_here`]. A control that
//! cannot fail is not a control, so this module adds two that can: the same
//! comparison at `mu0 = 2`, where the printed form is half the correct one,
//! and a sign-flipped body force.
//!
//! [ADR 026]: ../../../../../docs/adr/026_physical_mantle_simulation_critical_path.md

use super::benchmark_thieulot::ThieulotCase;
use super::MantleError;

/// Number of Simpson intervals used for every radial quadrature here. Even,
/// and large enough that the quadrature is not the error being reported:
/// the integrands are smooth on `[R1, R2]` and the closed-form check in
/// [`tests::the_free_slip_reference_is_a_solution_of_the_same_problem`]
/// reproduces [`ThieulotCase::v_rms`] to 13 digits with it.
const QUAD_INTERVALS: usize = 8192;

/// The analytic solution of the FREE-SLIP problem carrying Thieulot's body
/// force: same PDE, same `rho`, same `v_r = 0` walls, but zero tangential
/// traction instead of a prescribed tangential velocity.
///
/// Only `m = -1` (constant viscosity). The homogeneous poloidal basis used
/// here is the constant-viscosity one; a radial power-law viscosity has a
/// different one and this constructor refuses rather than returning a
/// plausible wrong answer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeSlipCase {
    case: ThieulotCase,
    /// Coefficients of the homogeneous poloidal modes `P = r, r^3, 1, r^-2`.
    c: [f64; 4],
}

/// `f`, `g`, `h` of the four homogeneous poloidal modes, and the traction
/// functional `f' - f/r` each of them contributes. Derived in the module
/// docs; every entry is re-checked numerically by
/// [`tests::the_homogeneous_poloidal_modes_carry_no_body_force`].
fn hom_f(j: usize, r: f64) -> f64 {
    match j {
        0 => -2.0,
        1 => -4.0 * r * r,
        2 => -1.0 / r,
        _ => 1.0 / (r * r * r),
    }
}
fn hom_g(j: usize, r: f64) -> f64 {
    match j {
        0 => 2.0,
        1 => 2.0 * r * r,
        2 => 2.0 / r,
        _ => 2.0 / (r * r * r),
    }
}
/// `h` per unit `mu0`.
fn hom_h_over_mu0(j: usize, r: f64) -> f64 {
    match j {
        0 => 0.0,
        1 => 20.0 * r,
        2 => 2.0 / (r * r),
        _ => 0.0,
    }
}
/// `f' - f/r`, the traction functional at a wall where `g = 0`.
fn hom_traction(j: usize, r: f64) -> f64 {
    match j {
        0 => 2.0 / r,
        1 => -4.0 * r,
        2 => 2.0 / (r * r),
        _ => -4.0 / (r * r * r * r),
    }
}

impl FreeSlipCase {
    /// Solve the four boundary conditions for the four homogeneous
    /// coefficients.
    pub fn new(case: ThieulotCase) -> Result<Self, MantleError> {
        if case.m() != -1 {
            return Err(MantleError::Field {
                field: "freeslip.m",
                value: format!("{}", case.m()),
                why: "the homogeneous poloidal basis {r, r^3, 1, r^-2} is the CONSTANT-viscosity \
                      one; a radial power-law viscosity has a different one and reusing this \
                      basis would return a plausible wrong answer",
            });
        }
        let (r1, r2) = (case.r_inner(), case.r_outer());
        let alpha = case.alpha();
        // Rows: v_r = 0 at R1 and R2 (via g), then zero tangential traction
        // at R1 and R2. The Thieulot part already has g(R1) = g(R2) = 0, so
        // the first two right-hand sides are its residual there (roundoff),
        // carried rather than assumed zero.
        let mut m = [[0.0f64; 4]; 4];
        for (j, col) in m[0].iter_mut().enumerate() {
            *col = hom_g(j, r1);
        }
        for (j, col) in m[1].iter_mut().enumerate() {
            *col = hom_g(j, r2);
        }
        for (j, col) in m[2].iter_mut().enumerate() {
            *col = hom_traction(j, r1);
        }
        for (j, col) in m[3].iter_mut().enumerate() {
            *col = hom_traction(j, r2);
        }
        let rhs = [
            -case.g(r1),
            -case.g(r2),
            3.0 * alpha / (r1 * r1 * r1),
            3.0 * alpha / (r2 * r2 * r2),
        ];
        let c = solve4(m, rhs)?;
        Ok(Self { case, c })
    }

    pub fn thieulot(&self) -> ThieulotCase {
        self.case
    }

    /// The four homogeneous coefficients, in the order `r, r^3, 1, r^-2`.
    pub fn coefficients(&self) -> [f64; 4] {
        self.c
    }

    /// The tangential shape `f(r)`: `v_theta = f(r) sin(theta)`.
    pub fn f(&self, r: f64) -> f64 {
        self.case.f(r) + (0..4).map(|j| self.c[j] * hom_f(j, r)).sum::<f64>()
    }

    /// The radial shape `g(r)`: `v_r = g(r) cos(theta)`.
    pub fn g(&self, r: f64) -> f64 {
        self.case.g(r) + (0..4).map(|j| self.c[j] * hom_g(j, r)).sum::<f64>()
    }

    /// The pressure shape `h(r)`: `p = h(r) cos(theta)`.
    pub fn h(&self, r: f64) -> f64 {
        let mu0 = self.case.mu0();
        // The Thieulot particular part's h is eq. (60); the module reaches it
        // through the public pressure at theta = 0, where cos(theta) = 1.
        self.case.pressure(r, 0.0)
            + mu0
                * (0..4)
                    .map(|j| self.c[j] * hom_h_over_mu0(j, r))
                    .sum::<f64>()
    }

    /// `(v_r, v_theta, v_phi)`. The toroidal component is ZERO: free slip
    /// admits only `f_t = b r`, the rigid rotation, and the gauge chosen here
    /// is `b = 0` - the same gauge the block preconditioner imposes.
    pub fn velocity_spherical(&self, r: f64, theta: f64) -> (f64, f64, f64) {
        (self.g(r) * theta.cos(), self.f(r) * theta.sin(), 0.0)
    }

    pub fn pressure(&self, r: f64, theta: f64) -> f64 {
        self.h(r) * theta.cos()
    }

    /// Cartesian velocity. Undefined exactly on the polar axis, like
    /// [`ThieulotCase::velocity_cartesian`], and callers must stay off it.
    pub fn velocity_cartesian(&self, p: [f64; 3]) -> [f64; 3] {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        let s = (p[0] * p[0] + p[1] * p[1]).sqrt();
        let (ct, st) = (p[2] / r, s / r);
        let (cp, sp) = (p[0] / s, p[1] / s);
        let (vr, vt, _) = self.velocity_spherical(r, st.atan2(ct));
        [
            vr * st * cp + vt * ct * cp,
            vr * st * sp + vt * ct * sp,
            vr * ct - vt * st,
        ]
    }

    pub fn pressure_cartesian(&self, p: [f64; 3]) -> f64 {
        let (r, theta) = ThieulotCase::spherical_of(p);
        self.pressure(r, theta)
    }

    /// `sqrt( (4 pi / (3 V)) * (B + 2 A) )` - the same eq. (63) structure as
    /// [`ThieulotCase::v_rms`], with `2 A` and not `4 A` because there is no
    /// `v_phi` here. Integrated by Simpson rather than in closed form.
    pub fn v_rms(&self) -> f64 {
        let a = simpson(self.case.r_inner(), self.case.r_outer(), |r| {
            let f = self.f(r);
            f * f * r * r
        });
        let b = simpson(self.case.r_inner(), self.case.r_outer(), |r| {
            let g = self.g(r);
            g * g * r * r
        });
        let v = self.case.volume();
        (4.0 * std::f64::consts::PI / (3.0 * v) * (b + 2.0 * a)).sqrt()
    }
}

/// Composite Simpson on `[a, b]` with [`QUAD_INTERVALS`] intervals.
pub fn simpson<F: Fn(f64) -> f64>(a: f64, b: f64, f: F) -> f64 {
    let n = QUAD_INTERVALS;
    let hh = (b - a) / n as f64;
    let mut acc = f(a) + f(b);
    for i in 1..n {
        let x = a + hh * i as f64;
        acc += if i % 2 == 1 { 4.0 * f(x) } else { 2.0 * f(x) };
    }
    acc * hh / 3.0
}

/// 4x4 Gaussian elimination with partial pivoting. Fails loud on a singular
/// system rather than returning an infinity that would later read as physics.
fn solve4(mut m: [[f64; 4]; 4], mut rhs: [f64; 4]) -> Result<[f64; 4], MantleError> {
    for col in 0..4 {
        let mut piv = col;
        for row in (col + 1)..4 {
            if m[row][col].abs() > m[piv][col].abs() {
                piv = row;
            }
        }
        if m[piv][col].abs() <= 1e-14 {
            return Err(MantleError::field(
                "freeslip boundary-condition matrix pivot",
                m[piv][col],
                "the four free-slip conditions are numerically dependent",
            ));
        }
        m.swap(col, piv);
        rhs.swap(col, piv);
        for row in (col + 1)..4 {
            let factor = m[row][col] / m[col][col];
            let (head, tail) = m.split_at_mut(row);
            for (slot, x) in tail[0].iter_mut().zip(head[col].iter()).skip(col) {
                *slot -= factor * x;
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    let mut out = [0.0f64; 4];
    for row in (0..4).rev() {
        let mut acc = rhs[row];
        for k in (row + 1)..4 {
            acc -= m[row][k] * out[k];
        }
        out[row] = acc / m[row][row];
    }
    Ok(out)
}

// ---------------------------------------------------------------------
// The grid harness

use super::stokes_saddle::FacetClass;
use super::stokes_solve::sampled_face_field;
use pf_grid::shell::ShellSpec;
use pf_grid::shell_ops::ShellTopology;

// Used only by the qualification tests below. Kept at module scope so the
// tests can reach them through `use super::*`, and `cfg(test)`-gated so the
// library build carries no unused import.
#[cfg(test)]
use super::stokes::ViscousOperator;
#[cfg(test)]
use super::stokes_saddle::{
    fgmres_saddle, BlockPreconditioner, CoarseSpace, Deflation, FgmresOptions, FgmresReport,
    GradientMode, ModeBasis, PreconditionerOptions, SaddleSystem,
};
#[cfg(test)]
use pf_grid::shell_ops::FaceGradientRule;

/// Build the NONDIMENSIONAL shell the reference lives on: `R1` and `R2` are
/// the case's own radii, not Earth's.
///
/// [`ThieulotCase::sample_on_shell`] refuses any other radii, and it is right
/// to: the solution is dimensionless and a metre-scaled shell would give a
/// silently meaningless field.
pub fn reference_topology(
    case: &ThieulotCase,
    n: usize,
    layers: usize,
) -> Result<ShellTopology, MantleError> {
    let spec = ShellSpec::new(n, layers, case.r_inner(), case.r_outer())
        .map_err(|e| MantleError::shell("reference shell spec", e))?;
    ShellTopology::build(&spec).map_err(|e| MantleError::Shell {
        what: "reference shell topology",
        detail: e.to_string(),
    })
}

/// `(face, i, j, k)` from a cell index. The inverse of
/// [`ShellSpec::flat`], written out here so this module does not have to
/// carry `pf_grid`'s `SpatialError` through every call.
fn fijk(spec: &ShellSpec, cell: usize) -> (usize, usize, usize, usize) {
    let (n, layers) = (spec.n(), spec.layers());
    let k = cell % layers;
    let rest = cell / layers;
    let j = rest % n;
    let rest = rest / n;
    let i = rest % n;
    (rest / n, i, j, k)
}

/// Panel-interior / seam / corner, by whether `i`, `j` sit on a panel edge.
/// The same rule [`super::stokes_saddle::boundary_moment_report`] uses.
fn cell_class(spec: &ShellSpec, cell: usize) -> FacetClass {
    let (_, i, j, _) = fijk(spec, cell);
    let n = spec.n();
    let ei = i == 0 || i + 1 == n;
    let ej = j == 0 || j + 1 == n;
    match (ei, ej) {
        (true, true) => FacetClass::Corner,
        (false, false) => FacetClass::PanelInterior,
        _ => FacetClass::Seam,
    }
}

fn worse(a: FacetClass, b: FacetClass) -> FacetClass {
    let rank = |c: FacetClass| match c {
        FacetClass::PanelInterior => 0,
        FacetClass::Seam => 1,
        FacetClass::Corner => 2,
    };
    if rank(a) >= rank(b) {
        a
    } else {
        b
    }
}

/// Whether a cell touches the inner or outer spherical wall.
fn cell_is_wall(spec: &ShellSpec, cell: usize) -> bool {
    let (_, _, _, k) = fijk(spec, cell);
    k == 0 || k + 1 == spec.layers()
}

/// Per-face class and wall-adjacency, computed once per rung.
pub struct FaceLabels {
    pub class: Vec<FacetClass>,
    pub wall: Vec<bool>,
    pub interior: Vec<bool>,
}

pub fn face_labels(topo: &ShellTopology) -> FaceLabels {
    let spec = topo.spec();
    let mut class = Vec::with_capacity(topo.face_count());
    let mut wall = Vec::with_capacity(topo.face_count());
    let mut interior = Vec::with_capacity(topo.face_count());
    for f in topo.faces() {
        let o = f.owner as usize;
        let (c, w) = match f.neighbor {
            Some(nb) => (
                worse(cell_class(spec, o), cell_class(spec, nb as usize)),
                cell_is_wall(spec, o) || cell_is_wall(spec, nb as usize),
            ),
            None => (cell_class(spec, o), true),
        };
        class.push(c);
        wall.push(w);
        interior.push(f.neighbor.is_some());
    }
    FaceLabels {
        class,
        wall,
        interior,
    }
}

/// One error series: `L2` in the `M_u` (or `M_p`) metric and `L-infinity`,
/// each relative to the reference field's own norm in the same metric.
#[derive(Debug, Clone, Copy)]
pub struct ErrorPair {
    pub l2_rel: f64,
    pub linf_rel: f64,
    pub count: usize,
}

/// Weighted `L2` and `L-infinity` of `a - b` over a selected face subset,
/// each relative to the norm of `b` over the SAME subset. Reporting the
/// relative error against a global denominator would let a subset with a tiny
/// reference value look accurate.
#[cfg(test)]
fn subset_error(weights: &[f64], a: &[f64], b: &[f64], keep: impl Fn(usize) -> bool) -> ErrorPair {
    let (mut num, mut den, mut linf, mut bmax, mut count) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0usize);
    for i in 0..a.len() {
        if !keep(i) {
            continue;
        }
        let d = a[i] - b[i];
        num += weights[i] * d * d;
        den += weights[i] * b[i] * b[i];
        linf = linf.max(d.abs());
        bmax = bmax.max(b[i].abs());
        count += 1;
    }
    ErrorPair {
        l2_rel: if den > 0.0 {
            (num / den).sqrt()
        } else {
            f64::NAN
        },
        linf_rel: if bmax > 0.0 { linf / bmax } else { f64::NAN },
        count,
    }
}

/// Volume-weighted RMS of a reconstructed cell velocity. The reconstruction
/// ([`ShellTopology::reconstruct_cell_velocity`]) is exact only for a
/// CONSTANT field, so this number carries the reconstruction's own error;
/// that is why every table below reports it for the exact face DOFs too, so
/// the reconstruction error and the solve error are separable.
pub fn reconstructed_v_rms(topo: &ShellTopology, u: &[f64]) -> Result<f64, MantleError> {
    let v = topo
        .reconstruct_cell_velocity(u)
        .map_err(|e| MantleError::Shell {
            what: "cell velocity reconstruction",
            detail: e.to_string(),
        })?;
    let vol = topo.cell_volumes_m3().map_err(|e| MantleError::Shell {
        what: "cell volumes",
        detail: e.to_string(),
    })?;
    let mut num = 0.0;
    let mut den = 0.0;
    for (c, vc) in v.iter().enumerate() {
        num += vol[c] * (vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);
        den += vol[c];
    }
    Ok((num / den).sqrt())
}

/// The body force `b = rho e_r` sampled onto face normals, zero on the walls.
///
/// SIGN, derived and then measured. Thieulot's balance is
/// `-grad p + div(2 mu eps) + rho e_r = 0`, so
/// `-div(2 mu eps) + grad p = rho e_r`; this module's operator pair is
/// `B u + G p` with `B = -A` the positive viscous block and `G` the gradient,
/// so the right-hand side is `rho e_r` with a PLUS sign. If that were wrong
/// the consistency residual of
/// [`tests::the_discrete_operators_are_consistent_with_the_continuum`] would
/// be `O(1)` in the panel interiors and would not fall with `h`; it falls.
pub fn body_force_faces(topo: &ShellTopology, density: impl Fn(f64, f64) -> f64) -> Vec<f64> {
    sampled_face_field(topo, |x| {
        let (r, theta) = ThieulotCase::spherical_of(x);
        let rho = density(r, theta);
        [rho * x[0] / r, rho * x[1] / r, rho * x[2] / r]
    })
}

/// A cell-centred field sampled at the POLYHEDRAL cell centroids - the same
/// points the corrected gradient differentiates about.
pub fn cell_field(
    topo: &ShellTopology,
    f: impl Fn([f64; 3]) -> f64,
) -> Result<Vec<f64>, MantleError> {
    let cent = topo.cell_centroids_m().map_err(|e| MantleError::Shell {
        what: "cell centroids",
        detail: e.to_string(),
    })?;
    Ok(cent.into_iter().map(f).collect())
}

/// The EXACT area-mean face DOF of an analytic vector field:
/// `u_f = (1/|A_f|) integral_F v . n dA`, by tensor-product Gauss-Legendre on
/// the bilinear facet patch.
///
/// **Why this and not a facet-centroid point sample.** The unknown IS the
/// area mean ([`super::stokes`] module docs). A point sample differs from it
/// by `O(h^2)`, which is harmless when the field itself is being compared -
/// but `B` is a second-derivative operator, so it AMPLIFIES a DOF error by
/// `h^-2` and an `O(h^2)` sampling error becomes an `O(1)` momentum residual.
/// Measuring operator consistency with point samples therefore cannot
/// distinguish an inconsistent operator from a second-order interpolation,
/// and the first version of this module made exactly that mistake.
///
/// Returns the DOFs and the worst relative disagreement between the vector
/// area this routine integrates and the one [`pf_grid`] recorded - a live
/// check that the facet vertex cycle replicated here is the same one
/// `pf_grid` used, since `ShellSpec::facet_vertices` is private.
pub fn area_mean_face_field(
    topo: &ShellTopology,
    v: impl Fn([f64; 3]) -> [f64; 3],
) -> Result<(Vec<f64>, f64), MantleError> {
    // 4-point Gauss-Legendre on [0,1], enough that the quadrature is not the
    // error: the integrand is smooth on a facet whose diameter is h.
    const GX: [f64; 4] = [
        0.069_431_844_202_973_71,
        0.330_009_478_207_571_87,
        0.669_990_521_792_428_1,
        0.930_568_155_797_026_3,
    ];
    const GW: [f64; 4] = [
        0.173_927_422_568_727,
        0.326_072_577_431_273,
        0.326_072_577_431_273,
        0.173_927_422_568_727,
    ];
    let spec = topo.spec().clone();
    let mut out = vec![0.0f64; topo.face_count()];
    let mut worst = 0.0f64;
    for c in 0..topo.cell_count() {
        let (fc, i, j, k) = fijk(&spec, c);
        let verts = spec
            .cell_vertices(fc, i, j, k)
            .map_err(|e| MantleError::shell("cell vertices", e))?;
        for (t, (fidx, sgn)) in topo.cell_faces(c as u32).into_iter().enumerate() {
            // Geometry is recorded from the OWNER's side, so integrate once,
            // from the owner, and let the stored area vector arbitrate.
            if sgn != 1 {
                continue;
            }
            let face = &topo.faces()[fidx as usize];
            if face.neighbor.is_none() {
                continue; // boundary facets carry no unknown
            }
            let lat = |a: usize, b: usize| [verts[a], verts[b], verts[b + 4], verts[a + 4]];
            let q = match t {
                0 => lat(0, 3),                                // Lateral(AlphaMinus)
                1 => lat(1, 2),                                // Lateral(AlphaPlus)
                2 => lat(0, 1),                                // Lateral(BetaMinus)
                3 => lat(3, 2),                                // Lateral(BetaPlus)
                4 => [verts[0], verts[1], verts[2], verts[3]], // Down
                _ => [verts[4], verts[5], verts[6], verts[7]], // Up
            };
            let mut flux = 0.0f64;
            let mut vec_area = [0.0f64; 3];
            for (a, sa) in GX.iter().enumerate() {
                for (b, sb) in GX.iter().enumerate() {
                    let (su, tv) = (*sa, *sb);
                    let mut x = [0.0f64; 3];
                    let mut xs = [0.0f64; 3];
                    let mut xt = [0.0f64; 3];
                    for d in 0..3 {
                        x[d] = (1.0 - su) * (1.0 - tv) * q[0][d]
                            + su * (1.0 - tv) * q[1][d]
                            + su * tv * q[2][d]
                            + (1.0 - su) * tv * q[3][d];
                        xs[d] = (1.0 - tv) * (q[1][d] - q[0][d]) + tv * (q[2][d] - q[3][d]);
                        xt[d] = (1.0 - su) * (q[3][d] - q[0][d]) + su * (q[2][d] - q[1][d]);
                    }
                    let nda = [
                        xs[1] * xt[2] - xs[2] * xt[1],
                        xs[2] * xt[0] - xs[0] * xt[2],
                        xs[0] * xt[1] - xs[1] * xt[0],
                    ];
                    let w = GW[a] * GW[b];
                    let vv = v(x);
                    for d in 0..3 {
                        flux += w * vv[d] * nda[d];
                        vec_area[d] += w * nda[d];
                    }
                }
            }
            // Live check against pf_grid's own vector area, and the sign the
            // stored normal was oriented to.
            let dot: f64 = (0..3).map(|d| vec_area[d] * face.area_vector_m2[d]).sum();
            let orient = if dot < 0.0 { -1.0 } else { 1.0 };
            let mismatch: f64 = (0..3)
                .map(|d| (orient * vec_area[d] - face.area_vector_m2[d]).powi(2))
                .sum::<f64>()
                .sqrt()
                / face.area_m2;
            worst = worst.max(mismatch);
            out[fidx as usize] = orient * flux / face.area_m2;
        }
    }
    Ok((out, worst))
}

/// Least-squares slope of `log y` against `log x`.
pub fn slope(xs: &[f64], ys: &[f64]) -> f64 {
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

/// Consecutive-pair slopes, which an aggregate fit can hide. The MC0 packet's
/// evidence README records exactly that failure, so both are always reported.
pub fn pair_slopes(xs: &[f64], ys: &[f64]) -> Vec<f64> {
    (1..xs.len())
        .map(|i| (ys[i] / ys[i - 1]).ln() / (xs[i] / xs[i - 1]).ln())
        .collect()
}

/// The refinement ladder. Default two rungs run in the normal debug suite;
/// `MCQ1_LADDER=full` adds `(8, 6)`, `deep` adds `(16, 12)` on top. Same
/// opt-in-cost / identical-coverage pattern as `MC0_SADDLE_LADDER`.
/// Outer FGMRES budget. The reported runs need 2000 to reach the residual
/// FLOOR this module measures; the default debug suite must not spend two
/// minutes proving it again, so the budget is opt-in with the same env var
/// as the ladder. Every table in the findings file was produced with
/// `MCQ1_LADDER` set, and every run states its own achieved residual, so a
/// short default run is never mistaken for a converged one.
pub fn max_outer() -> usize {
    if std::env::var("MCQ1_LADDER").is_ok() {
        2000
    } else {
        120
    }
}

pub fn ladder() -> Vec<(usize, usize)> {
    let mut v = vec![(2usize, 2usize), (4, 3)];
    match std::env::var("MCQ1_LADDER").as_deref() {
        Ok("full") => v.push((8, 6)),
        Ok("deep") => {
            v.push((8, 6));
            v.push((16, 12));
        }
        _ => {}
    }
    v
}

// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// `-grad p + mu lap v + rho e_r`, by central differences in Cartesian
    /// coordinates. The vector Laplacian of a Cartesian component IS the
    /// scalar Laplacian of that component, and `div(2 mu eps) = mu lap v` for
    /// constant `mu` and `div v = 0`, so no hand-transcribed spherical
    /// tensor divergence enters - which is where a transcription error would
    /// otherwise hide.
    fn fd_residual(
        v: impl Fn([f64; 3]) -> [f64; 3],
        p: impl Fn([f64; 3]) -> f64,
        force: impl Fn([f64; 3]) -> [f64; 3],
        mu: f64,
        x: [f64; 3],
        eps: f64,
    ) -> [f64; 3] {
        let v0 = v(x);
        let mut lap = [0.0f64; 3];
        let mut gp = [0.0f64; 3];
        for d in 0..3 {
            let mut xp = x;
            let mut xm = x;
            xp[d] += eps;
            xm[d] -= eps;
            let (vp, vm) = (v(xp), v(xm));
            for c in 0..3 {
                lap[c] += (vp[c] - 2.0 * v0[c] + vm[c]) / (eps * eps);
            }
            gp[d] = (p(xp) - p(xm)) / (2.0 * eps);
        }
        let b = force(x);
        [
            -gp[0] + mu * lap[0] + b[0],
            -gp[1] + mu * lap[1] + b[1],
            -gp[2] + mu * lap[2] + b[2],
        ]
    }

    /// Sample points well inside the shell and well off the polar axis, where
    /// `e_theta` and `e_phi` are undefined.
    fn probe_points(case: &ThieulotCase) -> Vec<[f64; 3]> {
        let mut out = Vec::new();
        for it in 0..24u32 {
            let t = 0.35 + 2.2 * f64::from(it % 8) / 8.0; // polar angle, off the poles
            let ph = 0.7 * f64::from(it);
            let rad = case.r_inner()
                + (case.r_outer() - case.r_inner()) * (0.12 + 0.7 * f64::from(it % 5) / 5.0);
            out.push([
                rad * t.sin() * ph.cos(),
                rad * t.sin() * ph.sin(),
                rad * t.cos(),
            ]);
        }
        out
    }

    fn norm3(v: [f64; 3]) -> f64 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    // -----------------------------------------------------------------
    // 1. The reference itself, before any grid

    #[test]
    fn the_reference_numbers_are_reproduced_from_the_formulas() {
        let case = ThieulotCase::isoviscous_reference();
        println!(
            "m = {}  R1 = {}  R2 = {}  gamma = {}  mu0 = {}",
            case.m(),
            case.r_inner(),
            case.r_outer(),
            case.gamma(),
            case.mu0()
        );
        println!(
            "alpha = {:.16}   (FINDINGS 1.9: -1.262358160777843)",
            case.alpha()
        );
        println!(
            "beta  = {:.16}   (FINDINGS 1.9:  3.000000000000000)",
            case.beta()
        );
        println!(
            "v_rms = {:.16}   (FINDINGS 1.9:  1.542073684542304)",
            case.v_rms()
        );
        // Recomputed independently by quadrature rather than trusting the
        // DERIVED table, as the research file instructs.
        let a = simpson(case.r_inner(), case.r_outer(), |r| {
            case.f(r) * case.f(r) * r * r
        });
        let b = simpson(case.r_inner(), case.r_outer(), |r| {
            case.g(r) * case.g(r) * r * r
        });
        let quad = (4.0 * PI / (3.0 * case.volume()) * (b + 4.0 * a)).sqrt();
        println!("v_rms by Simpson quadrature of eqs. (64)/(65) = {quad:.16}");
        assert!(
            (quad - case.v_rms()).abs() <= 1e-12 * case.v_rms(),
            "closed form {} and quadrature {} disagree",
            case.v_rms(),
            quad
        );
        assert!(
            (case.v_rms() - 1.542073684542304).abs() <= 1e-12,
            "v_rms {} does not reproduce the DERIVED reference 1.542073684542304",
            case.v_rms()
        );

        // Energy split. eq. (63)'s B + 4A is (4pi/3)*(B) radial plus
        // (8pi/3)*A for EACH of v_theta and v_phi.
        let e_r = 4.0 * PI / 3.0 * b;
        let e_t = 8.0 * PI / 3.0 * a;
        let tot = e_r + 2.0 * e_t;
        println!(
            "kinetic energy split: v_r {:.4}%  v_theta {:.4}%  v_phi {:.4}%",
            100.0 * e_r / tot,
            100.0 * e_t / tot,
            100.0 * e_t / tot
        );

        // Rigid-rotation content about z. The mode is v_phi = r sin(theta).
        let num = simpson(case.r_inner(), case.r_outer(), |r| case.f(r) * r * r * r);
        let den = simpson(case.r_inner(), case.r_outer(), |r| r * r * r * r);
        let omega = num / den;
        let e_rot = omega * omega * 2.0 * PI * (4.0 / 3.0) * den;
        println!(
            "rigid rotation about z: Omega_z = {omega:.9}, {:.4}% of the kinetic energy, \
             {:.4}% of |u| in L2",
            100.0 * e_rot / tot,
            100.0 * (e_rot / tot).sqrt()
        );
        assert!(
            (e_rot / tot).sqrt() > 0.2,
            "the benchmark was expected to carry SUBSTANTIAL rigid-rotation content; \
             measured only {:.4} of |u|",
            (e_rot / tot).sqrt()
        );

        // The tangential traction at the walls: NOT zero, which is the whole
        // finding. sigma_r_theta / (mu sin theta) = f' - f/r - g/r, and g = 0
        // at both walls.
        let fd = |r: f64| {
            let e = 1e-6;
            (case.f(r + e) - case.f(r - e)) / (2.0 * e)
        };
        for r in [case.r_inner(), case.r_outer()] {
            let t = fd(r) - case.f(r) / r;
            println!(
                "wall r = {r}: g = {:.3e},  (f' - f/r) = {t:.6}  (closed form -3 alpha / r^3 = {:.6})",
                case.g(r),
                -3.0 * case.alpha() / (r * r * r)
            );
            assert!(
                t.abs() > 1.0,
                "Thieulot's solution was expected to violate free slip at r = {r}; \
                 measured traction functional {t:.3e}"
            );
        }
    }

    // -----------------------------------------------------------------
    // 2. The negative control, and the fact that the built-in one is
    //    degenerate at exactly this case

    #[test]
    fn the_built_in_negative_control_is_degenerate_here() {
        // At m = -1 there is no r^m factor to lose (correction 1 does not
        // apply) and mu0 = 1 hides correction 2, so the "as printed" density
        // is BIT-IDENTICAL to the corrected one. Asserted, so nobody can
        // report a passing run against it as evidence of anything.
        let case = ThieulotCase::isoviscous_reference();
        let mut worst_bits = 0u64;
        for x in probe_points(&case) {
            let (r, th) = ThieulotCase::spherical_of(x);
            let a = case.density(r, th);
            let b = case.density_as_printed_in_eq59(r, th);
            worst_bits = worst_bits.max((a.to_bits() ^ b.to_bits()).count_ones() as u64);
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "at r = {r}, theta = {th}: corrected {a} vs as-printed {b}"
            );
        }
        println!(
            "[NEG-CONTROL] m = -1, mu0 = 1: density == density_as_printed_in_eq59 bit for bit \
             over {} probes (worst differing bits {worst_bits}). THE BUILT-IN CONTROL CANNOT \
             FAIL HERE and proves nothing about this qualification.",
            probe_points(&case).len()
        );

        // A control that CAN fail: same m = -1 solution, mu0 = 2. The printed
        // equations drop mu0 (FINDING B), so the as-printed density is then
        // exactly half the correct one.
        let case2 = ThieulotCase::new(-1, 2.0, 0.5, 1.0, -1.0).expect("mu0 = 2 case is admissible");
        let mut ratio_min = f64::MAX;
        let mut ratio_max = f64::MIN;
        for x in probe_points(&case2) {
            let (r, th) = ThieulotCase::spherical_of(x);
            let ratio = case2.density_as_printed_in_eq59(r, th) / case2.density(r, th);
            ratio_min = ratio_min.min(ratio);
            ratio_max = ratio_max.max(ratio);
        }
        println!(
            "[NEG-CONTROL] m = -1, mu0 = 2: as-printed / corrected in [{ratio_min:.15}, \
             {ratio_max:.15}] - the control discriminates, by exactly the factor mu0."
        );
        assert!((ratio_min - 0.5).abs() < 1e-15 && (ratio_max - 0.5).abs() < 1e-15);
    }

    // -----------------------------------------------------------------
    // 3. The free-slip reference

    #[test]
    fn the_homogeneous_poloidal_modes_carry_no_body_force() {
        // Each of P = r, r^3, 1, r^-2 must satisfy -grad p + lap v = 0 with
        // NO body force, or the free-slip construction is adding physics.
        let case = ThieulotCase::isoviscous_reference();
        let labels = ["P = r", "P = r^3", "P = 1", "P = r^-2"];
        for (j, label) in labels.into_iter().enumerate() {
            let v = |x: [f64; 3]| {
                let r = norm3(x);
                let s = (x[0] * x[0] + x[1] * x[1]).sqrt();
                let (ct, st) = (x[2] / r, s / r);
                let (cp, sp) = (x[0] / s, x[1] / s);
                let (vr, vt) = (hom_g(j, r) * ct, hom_f(j, r) * st);
                [
                    vr * st * cp + vt * ct * cp,
                    vr * st * sp + vt * ct * sp,
                    vr * ct - vt * st,
                ]
            };
            let p = |x: [f64; 3]| {
                let (r, th) = ThieulotCase::spherical_of(x);
                hom_h_over_mu0(j, r) * th.cos()
            };
            let mut worst = 0.0f64;
            let mut scale = 0.0f64;
            for x in probe_points(&case) {
                let res = fd_residual(v, p, |_: [f64; 3]| [0.0; 3], 1.0, x, 3e-5);
                worst = worst.max(norm3(res));
                scale = scale.max(norm3(v(x)));
            }
            println!("{label:8} max |-grad p + lap v| = {worst:.3e}   (field scale {scale:.3e})");
            assert!(
                worst <= 1e-3 * scale.max(1.0),
                "{label} is not a homogeneous Stokes solution: residual {worst:.3e}"
            );
        }
    }

    #[test]
    fn the_free_slip_reference_is_a_solution_of_the_same_problem() {
        let case = ThieulotCase::isoviscous_reference();
        let fs = FreeSlipCase::new(case).expect("free-slip case builds");
        println!(
            "homogeneous coefficients (r, r^3, 1, r^-2) = {:?}",
            fs.coefficients()
        );

        // (a) the boundary conditions it was built to satisfy, exactly
        for r in [case.r_inner(), case.r_outer()] {
            let e = 1e-6;
            let fp = (fs.f(r + e) - fs.f(r - e)) / (2.0 * e);
            let traction = fp - fs.f(r) / r - fs.g(r) / r;
            println!(
                "wall r = {r}: g_free = {:.3e}   traction functional = {traction:.3e}   \
                 h_free = {:.6}",
                fs.g(r),
                fs.h(r)
            );
            assert!(fs.g(r).abs() < 1e-12, "v_r != 0 at r = {r}");
            assert!(traction.abs() < 1e-8, "tangential traction != 0 at r = {r}");
        }

        // (b) it solves the SAME momentum equation with the SAME rho
        let mut worst = 0.0f64;
        let mut rho_scale = 0.0f64;
        for x in probe_points(&case) {
            let res = fd_residual(
                |y| fs.velocity_cartesian(y),
                |y| fs.pressure_cartesian(y),
                |y| {
                    let (r, th) = ThieulotCase::spherical_of(y);
                    let rho = case.density(r, th);
                    [rho * y[0] / r, rho * y[1] / r, rho * y[2] / r]
                },
                case.mu0(),
                x,
                3e-5,
            );
            worst = worst.max(norm3(res));
            let (r, th) = ThieulotCase::spherical_of(x);
            rho_scale = rho_scale.max(case.density(r, th).abs());
        }
        println!(
            "free-slip field vs Thieulot's OWN body force: max residual {worst:.3e} against \
             |rho| up to {rho_scale:.3e}  (relative {:.3e})",
            worst / rho_scale
        );
        assert!(
            worst <= 1e-5 * rho_scale,
            "the free-slip field does not solve the same momentum equation: {worst:.3e}"
        );

        // (c) Thieulot's own field, same check, same machinery - so a sign
        //     error in (b) could not pass unnoticed.
        let mut worst_t = 0.0f64;
        for x in probe_points(&case) {
            let res = fd_residual(
                |y| case.velocity_cartesian(y),
                |y| case.pressure_cartesian(y),
                |y| {
                    let (r, th) = ThieulotCase::spherical_of(y);
                    let rho = case.density(r, th);
                    [rho * y[0] / r, rho * y[1] / r, rho * y[2] / r]
                },
                case.mu0(),
                x,
                3e-5,
            );
            worst_t = worst_t.max(norm3(res));
        }
        println!("Thieulot's own field, same check: max residual {worst_t:.3e}");
        assert!(worst_t <= 1e-5 * rho_scale);

        // (d) the two answers, and how far apart the two problems are
        let vt = case.v_rms();
        let vf = fs.v_rms();
        println!(
            "v_rms: Thieulot {vt:.9}   free slip {vf:.9}   ratio {:.6}",
            vf / vt
        );
        let a_d = simpson(case.r_inner(), case.r_outer(), |r| {
            let d = fs.f(r) - case.f(r);
            d * d * r * r
        });
        let b_d = simpson(case.r_inner(), case.r_outer(), |r| {
            let d = fs.g(r) - case.g(r);
            d * d * r * r
        });
        let c_d = simpson(case.r_inner(), case.r_outer(), |r| {
            case.f(r) * case.f(r) * r * r
        });
        let diff = (4.0 * PI / (3.0 * case.volume()) * (b_d + 2.0 * a_d + 2.0 * c_d)).sqrt();
        println!(
            "|v_freeslip - v_Thieulot|_L2 / |v_Thieulot|_L2 = {:.6}  \
             (poloidal part {:.6}, toroidal part {:.6})",
            diff / vt,
            (4.0 * PI / (3.0 * case.volume()) * (b_d + 2.0 * a_d)).sqrt() / vt,
            (4.0 * PI / (3.0 * case.volume()) * (2.0 * c_d)).sqrt() / vt
        );
    }

    // -----------------------------------------------------------------
    // 4. The discrete operators against the continuum, with NO solve

    /// The four face subsets every table below is split by. `wall` means the
    /// face touches a cell that owns an inner or outer boundary facet, which
    /// is where the two boundary-value problems differ.
    /// A named subset of indices: a label and a predicate.
    type Subsets<'a> = Vec<(&'static str, Box<dyn Fn(usize) -> bool + 'a>)>;

    fn face_subsets<'a>(lab: &'a FaceLabels) -> Subsets<'a> {
        vec![
            ("global   ", Box::new(move |i: usize| lab.interior[i])),
            (
                "panel-int",
                Box::new(move |i: usize| {
                    lab.interior[i] && !lab.wall[i] && lab.class[i] == FacetClass::PanelInterior
                }),
            ),
            (
                "seam     ",
                Box::new(move |i: usize| {
                    lab.interior[i] && !lab.wall[i] && lab.class[i] != FacetClass::PanelInterior
                }),
            ),
            (
                "wall-adj ",
                Box::new(move |i: usize| lab.interior[i] && lab.wall[i]),
            ),
        ]
    }

    /// Fractions of `|v|_{M_u}` in the three truncated translations, the three
    /// rotations, and everything else.
    fn mode_split(six: &ModeBasis, sys: &SaddleSystem, v: &[f64]) -> (f64, f64, f64) {
        let c = six.components(sys, v);
        let tra = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        let rot = (c[3] * c[3] + c[4] * c[4] + c[5] * c[5]).sqrt();
        let rest = (1.0 - tra * tra - rot * rot).max(0.0).sqrt();
        (tra, rot, rest)
    }

    struct Rung<'a> {
        sys: SaddleSystem<'a>,
        labels: FaceLabels,
        none: ModeBasis,
        six: ModeBasis,
    }

    fn build_rung<'a>(topo: &'a ShellTopology, op: &'a ViscousOperator, mu0: f64) -> Rung<'a> {
        let sys = SaddleSystem::build(topo, op, mu0, FaceGradientRule::NonOrthogonalAndSkew)
            .expect("saddle system builds on the nondimensional shell");
        let none = ModeBasis::build(topo, &sys, Deflation::None).expect("empty basis");
        let six = ModeBasis::build(topo, &sys, Deflation::Six).expect("six modes");
        let labels = face_labels(topo);
        Rung {
            sys,
            labels,
            none,
            six,
        }
    }

    /// Face DOFs of an analytic velocity: the EXACT area mean, which is what
    /// the unknown is. See [`area_mean_face_field`] for why a facet-centroid
    /// point sample is not adequate here.
    fn analytic_faces(topo: &ShellTopology, v: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<f64> {
        let (out, geom) = area_mean_face_field(topo, v).expect("area-mean sampling");
        assert!(
            geom < 1e-12,
            "the facet vertex cycle replicated in area_mean_face_field disagrees with pf_grid's \
             own vector area by {geom:.3e} relative - the mapping is wrong"
        );
        for (i, x) in out.iter().enumerate() {
            assert!(
                x.is_finite(),
                "face {i} sampled a non-finite analytic velocity - a facet touches the polar \
                 axis, where e_theta and e_phi do not exist"
            );
        }
        out
    }

    #[test]
    fn the_discrete_operators_are_consistent_with_the_continuum() {
        // NO SOLVE. Substitute each analytic solution into the discrete
        // equations and measure what is left over. This is the truncation
        // error of the operator pair on a genuine continuum Stokes solution -
        // the first time this project has measured that - and it localises
        // the boundary-value-problem mismatch to the cells that touch a wall.
        let case = ThieulotCase::isoviscous_reference();
        let fs = FreeSlipCase::new(case).expect("free-slip case");
        println!(
            "consistency residual of an EXACT analytic solution in the discrete equations.\n\
             momentum rows relative to |rho e_r|_Mu; continuity row as theta |D u|_Mp / |f|."
        );
        let mut hs: Vec<f64> = Vec::new();
        let mut fs_global: Vec<f64> = Vec::new();
        let mut fs_panel: Vec<f64> = Vec::new();
        let mut fs_seam: Vec<f64> = Vec::new();
        let mut fs_wall: Vec<f64> = Vec::new();
        for (n, layers) in ladder() {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            let sys = &r.sys;
            let m_u = sys.m_u().to_vec();
            let f = body_force_faces(&topo, |rad, th| case.density(rad, th));
            let f_norm = sys.norm_u(&f);
            println!(
                "\n--- n = {n}, layers = {layers}, h = {:.4e}, u_unknowns = {}, cells = {} \
                 |rho e_r|_Mu = {f_norm:.6e}",
                sys.mesh_scale(),
                sys.unknown_count(),
                sys.cell_count()
            );
            // The control for the sampling question itself: how far the
            // facet-centroid point sample is from the exact area mean, and
            // what that same error does after `B` (which is why the exact
            // mean is used everywhere below).
            {
                let exact = analytic_faces(&topo, |x| fs.velocity_cartesian(x));
                let point = sampled_face_field(&topo, |x| fs.velocity_cartesian(x));
                let d: Vec<f64> = point.iter().zip(exact.iter()).map(|(a, b)| a - b).collect();
                let bu = sys.viscous(&d).expect("B on the sampling error");
                println!(
                    "  sampling control: |point-sample - area-mean|_Mu / |u|_Mu = {:.3e}, and \
                     |B (that difference)|_Mu / |f|_Mu = {:.3e}",
                    sys.norm_u(&d) / sys.norm_u(&exact),
                    sys.norm_u(&bu) / f_norm
                );
            }
            for (label, u_fn, p_fn) in [
                (
                    "Thieulot ",
                    Box::new(|x: [f64; 3]| case.velocity_cartesian(x))
                        as Box<dyn Fn([f64; 3]) -> [f64; 3]>,
                    Box::new(|x: [f64; 3]| case.pressure_cartesian(x))
                        as Box<dyn Fn([f64; 3]) -> f64>,
                ),
                (
                    "free slip",
                    Box::new(|x: [f64; 3]| fs.velocity_cartesian(x)),
                    Box::new(|x: [f64; 3]| fs.pressure_cartesian(x)),
                ),
            ] {
                let u_ex = analytic_faces(&topo, &u_fn);
                let mut p_ex = cell_field(&topo, &p_fn).expect("pressure sample");
                sys.remove_mean_pressure(&mut p_ex).expect("gauge");
                let (mom, cont) = sys
                    .apply(&r.none, GradientMode::Corrected, &u_ex, &p_ex)
                    .expect("operator application");
                let res: Vec<f64> = mom.iter().zip(f.iter()).map(|(a, b)| a - b).collect();
                // Which term is big: the viscous one, the pressure gradient,
                // or the forcing. An O(1) residual made of three O(1) terms
                // that nearly cancel is a different statement from one made
                // of a single term that is simply wrong.
                let bu = sys.viscous(&u_ex).expect("B u_ex");
                let gp = sys
                    .gradient(&p_ex, GradientMode::Corrected)
                    .expect("G p_ex");
                println!(
                    "  {label}  term sizes / |f|_Mu:  |B u_ex| {:.3e}   |G p_ex| {:.3e}   \
                     |B u_ex + G p_ex| {:.3e}   |residual| {:.3e}",
                    sys.norm_u(&bu) / f_norm,
                    sys.norm_u(&gp) / f_norm,
                    sys.norm_u(&mom) / f_norm,
                    sys.norm_u(&res) / f_norm
                );
                let zero = vec![0.0f64; res.len()];
                print!("  {label}  momentum:");
                for (name, keep) in face_subsets(&r.labels) {
                    let e = subset_error(&m_u, &res, &zero, &keep);
                    // b = 0 there, so l2_rel is NaN; report the absolute
                    // weighted norm against |f| instead.
                    let mut num = 0.0;
                    let mut linf = 0.0f64;
                    for i in 0..res.len() {
                        if keep(i) {
                            num += m_u[i] * res[i] * res[i];
                            linf = linf.max(res[i].abs());
                        }
                    }
                    print!(
                        "  {name} L2 {:.3e} Linf {:.3e} ({})",
                        num.sqrt() / f_norm,
                        linf,
                        e.count
                    );
                }
                println!();
                let cont_rel = sys.theta() * sys.norm_p(&cont) / f_norm;
                println!("  {label}  continuity: theta |D u|_Mp / |f| = {cont_rel:.3e}");
                if label == "free slip" {
                    let mut vals = Vec::new();
                    for (_, keep) in face_subsets(&r.labels) {
                        let mut num = 0.0;
                        for i in 0..res.len() {
                            if keep(i) {
                                num += m_u[i] * res[i] * res[i];
                            }
                        }
                        vals.push(num.sqrt() / f_norm);
                    }
                    fs_global.push(vals[0]);
                    fs_panel.push(vals[1]);
                    fs_seam.push(vals[2]);
                    fs_wall.push(vals[3]);
                }
            }
            hs.push(sys.mesh_scale());
        }
        if hs.len() >= 2 {
            println!(
                "\nfree-slip momentum consistency, fitted over ALL {} rungs listed above:",
                hs.len()
            );
            for (name, series) in [
                ("global   ", &fs_global),
                ("panel-int", &fs_panel),
                ("seam     ", &fs_seam),
                ("wall-adj ", &fs_wall),
            ] {
                println!(
                    "  {name} ~ h^{:+.3}   consecutive-pair slopes {:?}",
                    slope(&hs, series),
                    pair_slopes(&hs, series)
                        .iter()
                        .map(|v| format!("{v:+.3}"))
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    // -----------------------------------------------------------------
    // 5. The driven solve - the actual MC-Q1 question

    /// Cell subsets: global, wall-adjacent, and the rest.
    fn cell_subsets(spec: &ShellSpec) -> Subsets<'_> {
        vec![
            ("global   ", Box::new(move |_: usize| true)),
            (
                "non-wall ",
                Box::new(move |c: usize| !cell_is_wall(spec, c)),
            ),
            ("wall-adj ", Box::new(move |c: usize| cell_is_wall(spec, c))),
        ]
    }

    struct Solved {
        u: Vec<f64>,
        p: Vec<f64>,
        rep: FgmresReport,
    }

    /// One driven solve. `coarse` is the only thing varied between the two
    /// configurations reported below; every tolerance is
    /// [`PreconditionerOptions::default`] / [`FgmresOptions::default`] with
    /// `max_iterations` raised to 2000, matching the MC0 ladder run so the
    /// counts are comparable with that packet's table.
    #[allow(clippy::too_many_arguments)]
    fn drive(
        topo: &ShellTopology,
        op: &ViscousOperator,
        sys: &SaddleSystem,
        six: &ModeBasis,
        none: &ModeBasis,
        force: &[f64],
        coarse: CoarseSpace,
        stall_window: usize,
        grad: GradientMode,
        inner_tol: f64,
    ) -> Solved {
        let mut pre = BlockPreconditioner::build(
            topo,
            op,
            sys,
            PreconditionerOptions {
                coarse,
                inner_tol,
                ..PreconditionerOptions::default()
            },
        )
        .expect("preconditioner builds");
        let (u, p, rep) = fgmres_saddle(
            sys,
            &mut pre,
            none,
            six,
            grad,
            force,
            &sys.zeros_p(),
            FgmresOptions {
                max_iterations: max_outer(),
                stall_window,
                ..FgmresOptions::default()
            },
        )
        .expect("FGMRES runs");
        Solved { u, p, rep }
    }

    #[test]
    fn the_stall_on_the_physical_right_hand_side_is_real() {
        // The solve stops short of 1e-8 on this right-hand side. Before that
        // is reported as a solver property it has to be separated from the
        // stall DETECTOR: a stall_window of 40 declares a stall after 40
        // Arnoldi steps without a 1.05x fall, and a slow-but-real descent
        // would trip it. Re-run with the detector effectively disabled and
        // report where the residual actually goes.
        let case = ThieulotCase::isoviscous_reference();
        let (n, layers) = (4usize, 3usize);
        let topo = reference_topology(&case, n, layers).expect("topology");
        let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
        let r = build_rung(&topo, &op, case.mu0());
        let f = body_force_faces(&topo, |rad, th| case.density(rad, th));
        for (coarse, window) in [
            (CoarseSpace::None, 40usize),
            (CoarseSpace::None, max_outer()),
            (CoarseSpace::Translations, max_outer()),
            (CoarseSpace::Six, max_outer()),
        ] {
            let sol = drive(
                &topo,
                &op,
                &r.sys,
                &r.six,
                &r.none,
                &f,
                coarse,
                window,
                GradientMode::Corrected,
                1e-2,
            );
            let tail: Vec<String> = sol
                .rep
                .history
                .iter()
                .rev()
                .take(4)
                .rev()
                .map(|v| format!("{v:.3e}"))
                .collect();
            println!(
                "coarse {:>15} window {window:>4}: outer {:>4} cycles {} converged {} stalled {} true_rel {:.3e} mom {:.3e} cont {:.3e} last4 {tail:?}",
                coarse.label(),
                sol.rep.outer_iterations,
                sol.rep.restart_cycles,
                sol.rep.converged,
                sol.rep.stalled,
                sol.rep.true_rel_residual,
                sol.rep.momentum_rel,
                sol.rep.continuity_rel
            );
            // WHERE the floor lives: the final residual resolved onto the six
            // A-null directions, each already relative to |b|.
            println!(
                "                 null residual: translations {:?} rotations {:?}",
                sol.rep.null_residual[0..3]
                    .iter()
                    .map(|v| format!("{v:.3e}"))
                    .collect::<Vec<_>>(),
                sol.rep.null_residual[3..6]
                    .iter()
                    .map(|v| format!("{v:.3e}"))
                    .collect::<Vec<_>>()
            );
            let fs = FreeSlipCase::new(case).expect("fs");
            let u_f = analytic_faces(&topo, |x| fs.velocity_cartesian(x));
            let m_u = r.sys.m_u().to_vec();
            let e = subset_error(&m_u, &sol.u, &u_f, |i| r.labels.interior[i]);
            let (tra, rot, _) = mode_split(&r.six, &r.sys, &sol.u);
            let rec = reconstructed_v_rms(&topo, &sol.u).expect("recon");
            println!(
                "                 |u_h - u_free|/|u_free| = {:.4e}  v_rms {rec:.6} (analytic free-slip {:.6})  solved mode content: tra {tra:.3e} rot {rot:.3e}",
                e.l2_rel,
                fs.v_rms()
            );
        }

        // Bisect the floor. Four one-variable changes against the row above
        // (coarse none, detector off, corrected gradient, inner_tol 1e-2,
        // true_rel 1.338e-4).
        let bench = super::super::stokes_saddle::buoyancy_shaped_force(&topo);
        for (label, force, grad, inner_tol) in [
            (
                "MC0's own l=2 forcing, same shell ",
                bench.clone(),
                GradientMode::Corrected,
                1e-2f64,
            ),
            (
                "rho e_r, IMPLICIT gradient only   ",
                f.clone(),
                GradientMode::Implicit,
                1e-2,
            ),
            (
                "rho e_r, inner CG tol 1e-10       ",
                f.clone(),
                GradientMode::Corrected,
                1e-10,
            ),
            (
                "rho e_r, theta-independent (l=0)  ",
                body_force_faces(&topo, |rad, _| case.density(rad, 0.0)),
                GradientMode::Corrected,
                1e-2,
            ),
            (
                "MC0 radial profile x cos(th) (l=1)",
                body_force_faces(&topo, |rad, th| {
                    let d = std::f64::consts::PI * (rad - case.r_inner())
                        / (case.r_outer() - case.r_inner());
                    d.sin() * th.cos()
                }),
                GradientMode::Corrected,
                1e-2,
            ),
            (
                "eq.50 radial profile x P2    (l=2)",
                body_force_faces(&topo, |rad, th| {
                    let c = th.cos();
                    case.density(rad, 0.0) * (1.5 * c * c - 0.5)
                }),
                GradientMode::Corrected,
                1e-2,
            ),
        ] {
            let sol = drive(
                &topo,
                &op,
                &r.sys,
                &r.six,
                &r.none,
                &force,
                CoarseSpace::None,
                max_outer(),
                grad,
                inner_tol,
            );
            println!(
                "{label}: outer {:>5} converged {:>5} true_rel {:.3e} mom {:.3e} cont {:.3e}",
                sol.rep.outer_iterations,
                sol.rep.converged,
                sol.rep.true_rel_residual,
                sol.rep.momentum_rel,
                sol.rep.continuity_rel
            );
        }

        // The floor is entirely in the CONTINUITY row and only for the l = 1
        // forcing. This is the candidate cause, measured directly: the
        // preconditioner's velocity output is M_u-orthogonal to all six
        // A-null modes, so the solved field has ZERO content there; the true
        // solution does not, and what is missing shows up as a continuity
        // residual because `D t` is O(1) on the truncated translations.
        let fs = FreeSlipCase::new(case).expect("fs");
        let u_f = analytic_faces(&topo, |x| fs.velocity_cartesian(x));
        let f_norm = r.sys.norm_u(&f);
        let comps = r.six.components(&r.sys, &u_f);
        let mut in_span = vec![0.0f64; u_f.len()];
        let nu = r.sys.norm_u(&u_f);
        for (j, c) in comps.iter().enumerate() {
            for (slot, x) in in_span.iter_mut().zip(r.six.mode(j).iter()) {
                *slot += c * nu * x;
            }
        }
        let perp: Vec<f64> = u_f.iter().zip(in_span.iter()).map(|(a, b)| a - b).collect();
        let cont = |v: &[f64]| {
            r.sys.theta() * r.sys.norm_p(&r.sys.divergence(v).expect("divergence")) / f_norm
        };
        println!(
            "exact free-slip DOFs: six-mode components {:?}",
            comps
                .iter()
                .map(|v| format!("{v:+.4e}"))
                .collect::<Vec<_>>()
        );
        println!(
            "continuity residual theta |D u|_Mp / |f|:  exact field {:.3e}   \
             exact MINUS its six-mode part {:.3e}   only its six-mode part {:.3e}   \
             (the solver floors at 1.34e-4)",
            cont(&u_f),
            cont(&perp),
            cont(&in_span)
        );
    }

    #[test]
    fn the_driven_solve_against_the_analytic_references() {
        let case = ThieulotCase::isoviscous_reference();
        let fs = FreeSlipCase::new(case).expect("free-slip case");
        println!(
            "MC-Q1. Driving the coupled saddle point with rho e_r from eq. (50).\n\
             GAUGES FIXED: pressure by volume-weighted zero mean over the shell (M_p), applied \
             to BOTH the solved and the analytic field; velocity rotation gauge is fixed at \
             zero by the block preconditioner (it projects all six A-null modes out of every \
             preconditioner output), which is reported, not chosen here.\n"
        );
        let mut hs: Vec<f64> = Vec::new();
        let mut e_glob: Vec<f64> = Vec::new();
        let mut e_panel: Vec<f64> = Vec::new();
        let mut e_seam: Vec<f64> = Vec::new();
        let mut e_wall: Vec<f64> = Vec::new();
        let mut e_pres: Vec<f64> = Vec::new();
        let mut e_vrms: Vec<f64> = Vec::new();
        for (n, layers) in ladder() {
            let topo = reference_topology(&case, n, layers).expect("topology");
            let op = ViscousOperator::isoviscous(&topo, case.mu0()).expect("operator");
            let r = build_rung(&topo, &op, case.mu0());
            let sys = &r.sys;
            let m_u = sys.m_u().to_vec();
            let m_p = sys.m_p().to_vec();
            let f = body_force_faces(&topo, |rad, th| case.density(rad, th));
            let sol = drive(
                &topo,
                &op,
                sys,
                &r.six,
                &r.none,
                &f,
                CoarseSpace::None,
                max_outer(),
                GradientMode::Corrected,
                1e-2,
            );

            println!(
                "=== n = {n} layers = {layers}  h = {:.4e}  u_unk = {}  p_unk = {}",
                sys.mesh_scale(),
                sys.unknown_count(),
                sys.cell_count()
            );
            println!(
                "  solver: outer {} cycles {} innerCG {} maxIn {} notconv {} converged {} \
                 stalled {} true_rel {:.3e} mom_rel {:.3e} cont_rel {:.3e}",
                sol.rep.outer_iterations,
                sol.rep.restart_cycles,
                sol.rep.inner_cg_iterations,
                sol.rep.inner_cg_max,
                sol.rep.inner_cg_not_converged,
                sol.rep.converged,
                sol.rep.stalled,
                sol.rep.true_rel_residual,
                sol.rep.momentum_rel,
                sol.rep.continuity_rel
            );

            let u_t = analytic_faces(&topo, |x| case.velocity_cartesian(x));
            let u_f = analytic_faces(&topo, |x| fs.velocity_cartesian(x));
            let mut p_t = cell_field(&topo, |x| case.pressure_cartesian(x)).expect("p_T");
            let mut p_f = cell_field(&topo, |x| fs.pressure_cartesian(x)).expect("p_F");
            let mut p_h = sol.p.clone();
            sys.remove_mean_pressure(&mut p_t).expect("gauge T");
            sys.remove_mean_pressure(&mut p_f).expect("gauge F");
            sys.remove_mean_pressure(&mut p_h).expect("gauge h");

            // --- the headline scalar
            let spec_owned = ShellSpec::new(n, layers, case.r_inner(), case.r_outer()).unwrap();
            println!(
                "  v_rms  analytic Thieulot {:.6}   analytic free-slip {:.6}   \
                 centroid-sampled Thieulot {:.6}",
                case.v_rms(),
                fs.v_rms(),
                case.v_rms_on_shell(&spec_owned).expect("v_rms on shell")
            );
            let rec_t = reconstructed_v_rms(&topo, &u_t).expect("recon T");
            let rec_f = reconstructed_v_rms(&topo, &u_f).expect("recon F");
            let rec_h = reconstructed_v_rms(&topo, &sol.u).expect("recon h");
            println!(
                "  v_rms  from face DOFs via cell reconstruction:  exact-Thieulot {rec_t:.6}   \
                 exact-free-slip {rec_f:.6}   SOLVED {rec_h:.6}"
            );
            println!(
                "  v_rms  SOLVED / analytic Thieulot = {:.6}      SOLVED / analytic free-slip \
                 = {:.6}      SOLVED / reconstructed-free-slip = {:.6}",
                rec_h / case.v_rms(),
                rec_h / fs.v_rms(),
                rec_h / rec_f
            );

            // --- velocity error, split by cell class
            for (name, reference) in [("vs Thieulot ", &u_t), ("vs free-slip", &u_f)] {
                print!("  velocity {name}:");
                for (sub, keep) in face_subsets(&r.labels) {
                    let e = subset_error(&m_u, &sol.u, reference, &keep);
                    print!("  {sub} L2 {:.4e} Linf {:.4e}", e.l2_rel, e.linf_rel);
                }
                println!();
            }
            // --- pressure error after the gauge fix
            for (name, reference) in [("vs Thieulot ", &p_t), ("vs free-slip", &p_f)] {
                print!("  pressure {name}:");
                for (sub, keep) in cell_subsets(topo.spec()) {
                    let e = subset_error(&m_p, &p_h, reference, &keep);
                    print!("  {sub} L2 {:.4e} Linf {:.4e}", e.l2_rel, e.linf_rel);
                }
                println!();
            }

            // --- the rotation-space decomposition (hazard 2)
            for (name, v) in [
                ("exact Thieulot ", u_t.clone()),
                ("exact free-slip", u_f.clone()),
                ("SOLVED         ", sol.u.clone()),
                (
                    "error vs Thieu ",
                    sol.u.iter().zip(u_t.iter()).map(|(a, b)| a - b).collect(),
                ),
                (
                    "error vs free  ",
                    sol.u.iter().zip(u_f.iter()).map(|(a, b)| a - b).collect(),
                ),
            ] {
                let (tra, rot, rest) = mode_split(&r.six, sys, &v);
                println!(
                    "  modes {name}: translations {tra:.4e}  rotations {rot:.4e}  \
                     everything else {rest:.4e}   (fractions of its own |.|_Mu)"
                );
            }
            let (_, rot_h, _) = mode_split(&r.six, sys, &sol.u);
            assert!(
                rot_h < 1e-8,
                "the solved velocity was expected to carry EXACTLY no rotation - the \
                 preconditioner projects all six A-null modes out of every output - but \
                 measured {rot_h:.3e}"
            );

            hs.push(sys.mesh_scale());
            let mut vals = Vec::new();
            for (_, keep) in face_subsets(&r.labels) {
                vals.push(subset_error(&m_u, &sol.u, &u_f, &keep).l2_rel);
            }
            e_glob.push(vals[0]);
            e_panel.push(vals[1]);
            e_seam.push(vals[2]);
            e_wall.push(vals[3]);
            e_pres.push(subset_error(&m_p, &p_h, &p_f, |_| true).l2_rel);
            e_vrms.push((rec_h - fs.v_rms()).abs() / fs.v_rms());
        }
        if hs.len() >= 2 {
            println!(
                "\nconvergence against the FREE-SLIP reference, fitted over ALL {} rungs run \
                 (h = {:?}):",
                hs.len(),
                hs.iter().map(|v| format!("{v:.4e}")).collect::<Vec<_>>()
            );
            for (name, series) in [
                ("velocity L2 global   ", &e_glob),
                ("velocity L2 panel-int", &e_panel),
                ("velocity L2 seam     ", &e_seam),
                ("velocity L2 wall-adj ", &e_wall),
                ("pressure L2 global   ", &e_pres),
                ("v_rms scalar         ", &e_vrms),
            ] {
                println!(
                    "  {name} ~ h^{:+.3}   pairs {:?}   values {:?}",
                    slope(&hs, series),
                    pair_slopes(&hs, series)
                        .iter()
                        .map(|v| format!("{v:+.3}"))
                        .collect::<Vec<_>>(),
                    series
                        .iter()
                        .map(|v| format!("{v:.4e}"))
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    // -----------------------------------------------------------------
    // 6. Negative controls on the SOLVE, not just on the formula

    #[test]
    fn the_negative_controls_move_the_solved_answer() {
        // The point of hazard 1: if a wrong right-hand side gives the same
        // verdict as the right one, the qualification is measuring nothing.
        // Three controls, all at one rung, all reported against the SAME
        // free-slip reference.
        let case = ThieulotCase::isoviscous_reference();
        let fs = FreeSlipCase::new(case).expect("free-slip case");
        let case2 = ThieulotCase::new(-1, 2.0, 0.5, 1.0, -1.0).expect("mu0 = 2 case");
        let (n, layers) = (4usize, 3usize);
        let topo = reference_topology(&case, n, layers).expect("topology");

        for (label, mu0, density) in [
            (
                "corrected rho, mu0 = 1 (the real run)",
                1.0f64,
                Box::new(|rad: f64, th: f64| case.density(rad, th)) as Box<dyn Fn(f64, f64) -> f64>,
            ),
            (
                "[NEG-CONTROL] corrected rho, mu0 = 2 ",
                2.0,
                Box::new(|rad: f64, th: f64| case2.density(rad, th)),
            ),
            (
                "[NEG-CONTROL] AS-PRINTED rho, mu0 = 2",
                2.0,
                Box::new(|rad: f64, th: f64| case2.density_as_printed_in_eq59(rad, th)),
            ),
            (
                "[NEG-CONTROL] sign-flipped rho, mu0=1",
                1.0,
                Box::new(|rad: f64, th: f64| -case.density(rad, th)),
            ),
        ] {
            let op = ViscousOperator::isoviscous(&topo, mu0).expect("operator");
            let r = build_rung(&topo, &op, mu0);
            let sys = &r.sys;
            let m_u = sys.m_u().to_vec();
            let f = body_force_faces(&topo, &density);
            let sol = drive(
                &topo,
                &op,
                sys,
                &r.six,
                &r.none,
                &f,
                CoarseSpace::None,
                max_outer(),
                GradientMode::Corrected,
                1e-2,
            );
            let u_f = analytic_faces(&topo, |x| fs.velocity_cartesian(x));
            let e = subset_error(&m_u, &sol.u, &u_f, |i| r.labels.interior[i]);
            let rec = reconstructed_v_rms(&topo, &sol.u).expect("recon");
            println!(
                "{label}: outer {:>4} true_rel {:.3e} converged {:>5}  |u_h - u_free|/|u_free| \
                 = {:.4e}   v_rms(solved) = {rec:.6}  (free-slip analytic {:.6})",
                sol.rep.outer_iterations,
                sol.rep.true_rel_residual,
                sol.rep.converged,
                e.l2_rel,
                fs.v_rms()
            );
        }
        println!(
            "Expected, from the algebra: mu0 = 2 with the corrected rho reproduces the SAME \
             velocity (both rho and mu scale), the as-printed rho at mu0 = 2 is exactly half of \
             it (FINDING B), and the sign flip is exactly minus it. A qualification that gave \
             the same verdict for all four would be measuring nothing."
        );
    }
}
