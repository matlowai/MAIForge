//! MC0: the coupled saddle-point Stokes solve. Velocity and pressure are
//! solved TOGETHER here for the first time in this project.
//!
//! ```text
//! [ B   G ] [u]   [f]          B = M_u^-1 E^T W E = -A     (viscous, PSD)
//! [ D   0 ] [p] = [g]          G = G_imp + G_exp           (corrected gradient)
//!                              D = conservative divergence
//! ```
//!
//! Nothing in [`super::stokes`] or [`super::stokes_solve`] is modified. This
//! module reads their public surface only, and its `A`-block solves ARE
//! [`super::stokes_solve::deflated_cg`] - the conditioning measured in
//! `docs/reviews/planet-foundry-b1/MC0_ABLOCK_SOLVE_FINDINGS.md`
//! (`kappa(A) = O(h^-2)`, not worse) is the licence for using a plain
//! Jacobi-preconditioned CG as the inner solve rather than multigrid.
//!
//! The write-up is `docs/reviews/planet-foundry-b1/MC0_SADDLE_FINDINGS.md`.
//!
//! # The sign convention, once
//!
//! `B` is the POSITIVE semidefinite viscous block ([`super::stokes_solve::ABlock`]
//! in [`Metric::FaceMass`], where `B = -A` exactly), so the momentum row here
//! is `B u + G p = f` and `f` is the body force divided by `M_u`. Writing the
//! system around `-A` rather than `A` is what makes the block preconditioner
//! below a positive-definite object; it is bookkeeping, not physics.
//!
//! # Which modes are EXACT and which are NEAR-null: the whole difficulty
//!
//! [`super::stokes_solve::NULL_DIM`] is six, and all six are exact null modes
//! of the `A` BLOCK - three truncated translations and three rotations,
//! annihilated to `4e-15` relative. **Not one of them is an exact null mode
//! of the COUPLED system**, and they fail to be for two completely different
//! reasons, which is the single most important fact in this module:
//!
//! | mode | `K n = 0`? | `D n = 0`? | status for the saddle point |
//! |---|---|---|---|
//! | 3 truncated translations | exact | **no, `O(1)`** | strongly non-null; three well-conditioned equations |
//! | 3 rotations | exact | **no, `O(h^2)`** | NEAR-null: the boundary flux free slip deletes |
//! | 1 constant pressure | - | `D^T M_p 1 = 0` exact | the ONE exact null mode of the coupled system |
//!
//! So the coupled system has exactly ONE exact null direction (the pressure
//! constant), and it is the only thing this module projects out of the
//! equations. **The six `A`-null modes are deflated inside the `A`-block
//! PRECONDITIONER, where they are exact, and nowhere else by default.**
//!
//! Projecting a near-null mode out of the discrete equations does not solve
//! the system faster, it solves a DIFFERENT system: the three rotation-tested
//! momentum rows `-(D r_j)^T M_p p = <r_j, f>_{M_u}` are real equations with
//! small coefficients, not redundancies. [`Deflation`] exists so that
//! removing them can be MEASURED against not removing them
//! ([`tests::the_near_null_rotations_and_what_deflating_them_costs`]), and
//! every variant that removes a velocity mode is labelled in the report as
//! changing the discrete system.
//!
//! # Where the rotation leakage comes from, exactly
//!
//! For a rigid rotation `u = omega x x`, the flux through a facet is
//! `q_F = integral_F (omega x x) . n dS = omega . m_F` with
//! `m_F = integral_F (x x n) dS` - which is exactly
//! [`pf_grid::shell_ops::Face::moment_m3`], already stored for every face
//! including the boundary ones. A closed polyhedron has `sum_f s(c,f) m_f = 0`
//! identically, so
//!
//! ```text
//! (D r)_c = -(1 / V_c) sum over BOUNDARY facets of c of s(c,f) omega . m_f
//! ```
//!
//! is an EXACT identity, not an estimate: the divergence a rotation carries is
//! precisely the flux that impermeability deletes, it is nonzero only in
//! boundary cells, and it is zero if and only if `m_F = 0` on the boundary
//! facets. On the true sphere `n = x / r` so `x x n` vanishes pointwise; on a
//! flat facet it does not, and `m_F = |F| (xbar_F x n_F)` is the residue.
//! [`tests::the_deleted_boundary_flux_predicts_the_rotation_divergence`]
//! verifies the identity to roundoff and reports `|m_F|` split by cell class.
//!
//! # Honest scope
//!
//! - [`tests::the_manufactured_saddle_pair_is_recovered`] is a
//!   SELF-CONSISTENCY check: the pair is manufactured in the DISCRETE space
//!   and the right-hand side is formed by applying these operators to it, so
//!   it verifies that the solver inverts the operators this repository built
//!   and says NOTHING about whether they are the right physics. The physics
//!   check is MC-Q1 against [`super::benchmark_thieulot`] and it is the next
//!   packet, not this one.
//! - No claim is made here about convergence order, analytic error, or any
//!   quantity that requires a continuum solution.

use super::stokes::ViscousOperator;
use super::stokes_solve::{
    deflated_cg, stiffness_diagonal, ABlock, CgOptions, Metric, NullSpace, NULL_DIM,
};
use super::MantleError;
use pf_grid::shell_ops::{FaceCorrection, FaceGradientRule, ShellTopology};

fn shell(what: &'static str, e: pf_grid::shell_ops::OpsError) -> MantleError {
    // The grid's own message, unedited, same as `MantleError::shell` does for
    // `ShellError`. `OpsError` is a different type, so it needs its own wrap.
    MantleError::Shell {
        what,
        detail: e.to_string(),
    }
}

// ---------------------------------------------------------------------
// The system

/// Which gradient the Krylov operator applies.
///
/// The deferred non-orthogonal + skewness correction is EXPLICIT by
/// construction in `pf_grid` - it never enters the matrix, so it cannot
/// perturb the adjoint identity `G_imp = -M_u^-1 D^T M_p`. That leaves a
/// choice this module measures rather than assumes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GradientMode {
    /// `G = G_imp` only. The correction is then carried by an outer
    /// fixed-point loop ([`deferred_correction_solve`]), lagged or
    /// under-relaxed, which is the classical finite-volume treatment.
    Implicit,
    /// `G = G_imp + G_exp`, the FULL corrected gradient, applied inside the
    /// Krylov operator. This costs the symmetry of the saddle point - which
    /// FGMRES does not need - and removes the fixed point entirely.
    Corrected,
}

/// The assembled saddle-point system, matrix free.
///
/// Velocity vectors are FULL face length with boundary entries pinned to
/// zero, exactly as [`ABlock`] carries them; pressure vectors are cell
/// length. `M_u` is the over-relaxed face mass
/// ([`pf_grid::shell_ops::ShellTopology::face_masses_over_relaxed`], the one
/// [`ViscousOperator`] already uses) and `M_p = diag(V_c)` is the spherical
/// cell volume - the same geometry the buoyancy forcing will use.
#[derive(Debug)]
pub struct SaddleSystem<'a> {
    topo: &'a ShellTopology,
    /// `B = M_u^-1 K = -A`, reused from `stokes_solve` rather than rebuilt.
    mu_block: ABlock<'a>,
    corr: Vec<FaceCorrection>,
    volumes: Vec<f64>,
    centroids: Vec<[f64; 3]>,
    m_u: Vec<f64>,
    eta: f64,
    rule: FaceGradientRule,
    h_ref: f64,
    theta: f64,
}

impl<'a> SaddleSystem<'a> {
    /// Build the system for one shell and one isoviscous operator.
    ///
    /// `eta` is passed separately because [`ViscousOperator`] stores only
    /// `W_c = 2 eta V_c` and the Schur scaling below needs `eta` itself. It
    /// is CHECKED against the operator's own weights rather than trusted
    /// ([`MantleError::FieldPair`] if they disagree), so a caller cannot
    /// silently hand in a viscosity the operator was not built with.
    pub fn build(
        topo: &'a ShellTopology,
        op: &'a ViscousOperator,
        eta: f64,
        rule: FaceGradientRule,
    ) -> Result<Self, MantleError> {
        if !eta.is_finite() || eta <= 0.0 {
            return Err(MantleError::field(
                "eta",
                eta,
                "viscosity must be finite and strictly positive",
            ));
        }
        let volumes = topo
            .cell_volumes_m3()
            .map_err(|e| shell("cell volumes", e))?;
        // W_c = 2 eta V_c is the ONLY place eta enters the operator, so this
        // is an exact identity and not a tolerance question.
        for (c, (w, v)) in op.cell_weights().iter().zip(volumes.iter()).enumerate() {
            let want = 2.0 * eta * v;
            if (w - want).abs() > 1e-10 * want.abs() {
                return Err(MantleError::field_pair(
                    "eta given to SaddleSystem",
                    eta,
                    "eta implied by the operator's cell weight",
                    w / (2.0 * v),
                    "the operator was built with a different viscosity; the Schur scaling would be wrong",
                ));
            }
            let _ = c;
        }
        let corr = topo
            .face_corrections()
            .map_err(|e| shell("face corrections", e))?;
        let centroids = topo
            .cell_centroids_m()
            .map_err(|e| shell("cell centroids", e))?;
        let v_poly = topo
            .polyhedral_cell_volumes_m3()
            .map_err(|e| shell("polyhedral cell volumes", e))?;
        let h_ref = (v_poly.iter().sum::<f64>() / v_poly.len() as f64).cbrt();
        let m_u = op.face_mass().to_vec();
        let mu_block = ABlock::new(topo, op, Metric::FaceMass)?;
        Ok(Self {
            topo,
            mu_block,
            corr,
            volumes,
            centroids,
            m_u,
            eta,
            rule,
            h_ref,
            // theta converts a continuity residual (1/s) into the units of a
            // momentum residual (Pa/m) so the two halves of the residual norm
            // can be added. It is UNIT CONSISTENCY, derived: for the true
            // solution u ~ f h^2 / eta and D u ~ f h / eta, so theta = eta / h
            // makes theta |D u| ~ |f| ~ |B u|. It is not a tuning knob, and
            // `tests::the_residual_scale_is_not_a_tuning_knob` measures what
            // changing it by two decades in each direction does.
            theta: eta / h_ref,
        })
    }

    /// Multiply the residual norm's continuity weight `theta` by `mult`.
    ///
    /// Exists so that the derived value can be MEASURED against neighbouring
    /// ones rather than asserted to be right
    /// ([`tests::the_residual_scale_is_not_a_tuning_knob`]). Nothing in the
    /// library calls it with anything but `1.0`.
    pub fn with_residual_scale(mut self, mult: f64) -> Result<Self, MantleError> {
        if !mult.is_finite() || mult <= 0.0 {
            return Err(MantleError::field(
                "residual scale multiplier",
                mult,
                "the continuity weight must stay strictly positive",
            ));
        }
        self.theta *= mult;
        Ok(self)
    }

    pub fn topology(&self) -> &ShellTopology {
        self.topo
    }

    pub fn face_count(&self) -> usize {
        self.m_u.len()
    }

    pub fn cell_count(&self) -> usize {
        self.volumes.len()
    }

    pub fn unknown_count(&self) -> usize {
        self.mu_block.unknown_count()
    }

    pub fn eta(&self) -> f64 {
        self.eta
    }

    /// The mesh scale: cube root of the MEAN polyhedral cell volume, the same
    /// definition the `A`-block findings fit their refinement powers against.
    pub fn mesh_scale(&self) -> f64 {
        self.h_ref
    }

    pub fn theta(&self) -> f64 {
        self.theta
    }

    pub fn m_u(&self) -> &[f64] {
        &self.m_u
    }

    pub fn m_p(&self) -> &[f64] {
        &self.volumes
    }

    pub fn a_block(&self) -> &ABlock<'a> {
        &self.mu_block
    }

    pub fn zeros_u(&self) -> Vec<f64> {
        vec![0.0f64; self.face_count()]
    }

    pub fn zeros_p(&self) -> Vec<f64> {
        vec![0.0f64; self.cell_count()]
    }

    /// `<a, b>_{M_u}`, interior faces only. Delegates to [`ABlock::dot`], so
    /// the velocity inner product here IS the one the `A`-block spectrum was
    /// measured in.
    pub fn dot_u(&self, a: &[f64], b: &[f64]) -> f64 {
        self.mu_block.dot(a, b)
    }

    pub fn norm_u(&self, a: &[f64]) -> f64 {
        self.mu_block.norm(a)
    }

    /// `<a, b>_{M_p} = sum_c V_c a_c b_c`.
    pub fn dot_p(&self, a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .zip(self.volumes.iter())
            .map(|((x, y), v)| v * x * y)
            .sum()
    }

    pub fn norm_p(&self, a: &[f64]) -> f64 {
        self.dot_p(a, a).sqrt()
    }

    /// `|(r_u, r_p)|^2 = <r_u,r_u>_{M_u} + theta^2 <r_p,r_p>_{M_p}`.
    pub fn residual_norm(&self, ru: &[f64], rp: &[f64]) -> f64 {
        (self.dot_u(ru, ru) + self.theta * self.theta * self.dot_p(rp, rp)).sqrt()
    }

    /// `B u = M_u^-1 K u`, positive semidefinite. Exactly `-A`.
    pub fn viscous(&self, u: &[f64]) -> Result<Vec<f64>, MantleError> {
        self.mu_block.apply(u)
    }

    /// `G_imp p = lambda_f (p_nb - p_own) / d_f = -M_u^-1 D^T M_p p`, the
    /// implicit half. Exactly the negative weighted adjoint of [`Self::divergence`].
    pub fn gradient_implicit(&self, p: &[f64]) -> Result<Vec<f64>, MantleError> {
        self.topo
            .gradient_implicit(p, &self.corr)
            .map_err(|e| shell("implicit gradient", e))
    }

    /// The deferred correction `G_exp p`: non-orthogonality plus skewness,
    /// built on the least-squares cell gradient (the Green-Gauss one is
    /// measured in `pf_grid` NOT to converge on this mesh).
    pub fn gradient_explicit(&self, p: &[f64]) -> Result<Vec<f64>, MantleError> {
        let g = self
            .topo
            .cell_gradient_least_squares(p, &self.centroids)
            .map_err(|e| shell("least-squares cell gradient", e))?;
        self.topo
            .gradient_nonorthogonal_correction(p, &g, &self.corr, self.rule)
            .map_err(|e| shell("deferred gradient correction", e))
    }

    /// `G p = G_imp p + G_exp p`.
    pub fn gradient_corrected(&self, p: &[f64]) -> Result<Vec<f64>, MantleError> {
        let mut out = self.gradient_implicit(p)?;
        let ex = self.gradient_explicit(p)?;
        for (slot, x) in out.iter_mut().zip(ex.iter()) {
            *slot += x;
        }
        Ok(out)
    }

    pub fn gradient(&self, p: &[f64], mode: GradientMode) -> Result<Vec<f64>, MantleError> {
        match mode {
            GradientMode::Implicit => self.gradient_implicit(p),
            GradientMode::Corrected => self.gradient_corrected(p),
        }
    }

    /// `(D u)_c = (1/V_c) sum_f s(c,f) u_f |A_f|`, boundary facets excluded.
    pub fn divergence(&self, u: &[f64]) -> Result<Vec<f64>, MantleError> {
        self.topo
            .divergence(u, &self.volumes)
            .map_err(|e| shell("divergence", e))
    }

    /// Volume-weighted mean removal: the `M_p`-orthogonal projection onto the
    /// mean-free pressure space. Returns the mean that was removed.
    pub fn remove_mean_pressure(&self, p: &mut [f64]) -> Result<f64, MantleError> {
        self.topo
            .remove_mean_pressure(p, &self.volumes)
            .map_err(|e| shell("mean pressure removal", e))
    }

    /// One application of the coupled operator:
    /// `(B u + G p, D u)`, with the momentum row projected by `modes` (empty
    /// by default) and the continuity row made mean free.
    ///
    /// `D u` is ALREADY exactly `M_p`-mean-free - every interior facet
    /// appears in two cells with opposite sign, so the volume-weighted sum
    /// telescopes to zero - so the projection here measures roundoff rather
    /// than correcting anything. `tests::the_saddle_pieces_are_the_adjoint_pair`
    /// reports the drift.
    pub fn apply(
        &self,
        modes: &ModeBasis,
        mode: GradientMode,
        u: &[f64],
        p: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), MantleError> {
        let mut ru = self.viscous(u)?;
        let g = self.gradient(p, mode)?;
        for (slot, x) in ru.iter_mut().zip(g.iter()) {
            *slot += x;
        }
        modes.project(self, &mut ru);
        let mut rp = self.divergence(u)?;
        self.remove_mean_pressure(&mut rp)?;
        Ok((ru, rp))
    }
}

// ---------------------------------------------------------------------
// Velocity modes: which of the six, and what removing them means

/// Which of the six `A`-null velocity modes are projected out of the COUPLED
/// equations.
///
/// **[`Deflation::None`] is the default and the one the reported solve uses.**
/// The other three exist to be measured against it, and each of them changes
/// the discrete system: the mode's momentum row is a real equation with a
/// small coefficient, not a redundancy. Reporting a run with any variant
/// other than `None` without saying so would be reporting the residual of a
/// different problem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Deflation {
    /// Nothing removed. The only exact null direction of the coupled system
    /// is the pressure constant, and [`SaddleSystem::apply`] removes that.
    None,
    /// The three truncated translations. `K t = 0` exactly, but `D t` is
    /// `O(1)`: these are the STRONGLY non-null modes, and removing them
    /// discards three well-conditioned equations.
    Translations,
    /// The three rotations. `K r = 0` exactly and `D r` is `O(h^2)`: these
    /// are the NEAR-null modes. Removing them is the step that needs an
    /// explicit decision, because a solve that then converges is converging
    /// to a system whose rotation-tested momentum rows were deleted.
    Rotations,
    /// All six.
    Six,
}

impl Deflation {
    /// Indices into [`NullSpace::raw_modes`]: `0..3` truncated translations,
    /// `3..6` rotations.
    pub fn indices(self) -> &'static [usize] {
        match self {
            Deflation::None => &[],
            Deflation::Translations => &[0, 1, 2],
            Deflation::Rotations => &[3, 4, 5],
            Deflation::Six => &[0, 1, 2, 3, 4, 5],
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Deflation::None => "none",
            Deflation::Translations => "translations(3)",
            Deflation::Rotations => "rotations(3)",
            Deflation::Six => "six",
        }
    }

    /// Whether this variant alters the discrete system being solved. Only
    /// [`Deflation::None`] does not.
    pub fn changes_the_system(self) -> bool {
        !matches!(self, Deflation::None)
    }
}

/// A subset of the six raw modes, `M_u`-orthonormal.
///
/// Built from [`NullSpace::raw_modes`] and orthonormalised here rather than
/// by taking a slice of [`NullSpace::build`]'s output, because Gram-Schmidt
/// over all six mixes the sets: the orthonormal "rotation" vectors of a
/// six-mode basis are rotations MINUS their translation components, and
/// projecting with those would not be the orthogonal projection onto the
/// rotation span. The equality of the two spans for the full six is measured
/// in [`tests::the_mode_basis_agrees_with_the_a_block_null_space`].
#[derive(Debug, Clone)]
pub struct ModeBasis {
    modes: Vec<Vec<f64>>,
    which: Deflation,
    smallest_pivot: f64,
}

impl ModeBasis {
    pub fn build(
        topo: &ShellTopology,
        sys: &SaddleSystem,
        which: Deflation,
    ) -> Result<Self, MantleError> {
        Self::from_indices(topo, sys, which, which.indices())
    }

    fn from_indices(
        topo: &ShellTopology,
        sys: &SaddleSystem,
        which: Deflation,
        indices: &[usize],
    ) -> Result<Self, MantleError> {
        let raw = NullSpace::raw_modes(topo);
        let mut modes: Vec<Vec<f64>> = indices.iter().map(|i| raw[*i].clone()).collect();
        let mut smallest_pivot = f64::MAX;
        for j in 0..modes.len() {
            let before = sys.norm_u(&modes[j]);
            if !before.is_finite() || before <= 0.0 {
                return Err(MantleError::field(
                    "velocity mode norm",
                    before,
                    "a raw velocity mode has no content in the M_u inner product",
                ));
            }
            let (head, tail) = modes.split_at_mut(j);
            let target = &mut tail[0];
            // Twice-repeated modified Gram-Schmidt, for the same reason
            // NullSpace::build gives: a unit translation and a rotation
            // carrying the planet's radius are scaled 6 decades apart.
            for _pass in 0..2 {
                for prev in head.iter() {
                    let c = sys.dot_u(prev, target);
                    for (slot, x) in target.iter_mut().zip(prev.iter()) {
                        *slot -= c * x;
                    }
                }
            }
            let after = sys.norm_u(target);
            let pivot = after / before;
            smallest_pivot = smallest_pivot.min(pivot);
            if !pivot.is_finite() || pivot <= 1e-8 {
                return Err(MantleError::field(
                    "velocity mode Gram-Schmidt pivot",
                    pivot,
                    "the chosen velocity modes are numerically dependent in M_u",
                ));
            }
            let inv = 1.0 / after;
            for slot in target.iter_mut() {
                *slot *= inv;
            }
        }
        if modes.is_empty() {
            smallest_pivot = 1.0;
        }
        Ok(Self {
            modes,
            which,
            smallest_pivot,
        })
    }

    /// The three rotations alone, `M_u`-orthonormal, regardless of what this
    /// basis is used to project. Needed by the near-null coarse space, which
    /// SUPPLIES the rotations to the preconditioner instead of deleting them.
    pub fn rotations(topo: &ShellTopology, sys: &SaddleSystem) -> Result<Self, MantleError> {
        Self::from_indices(topo, sys, Deflation::Rotations, &[3, 4, 5])
    }

    pub fn which(&self) -> Deflation {
        self.which
    }

    pub fn dim(&self) -> usize {
        self.modes.len()
    }

    pub fn smallest_pivot(&self) -> f64 {
        self.smallest_pivot
    }

    pub fn mode(&self, index: usize) -> &[f64] {
        &self.modes[index]
    }

    pub fn project(&self, sys: &SaddleSystem, v: &mut [f64]) {
        for m in self.modes.iter() {
            let c = sys.dot_u(m, v);
            for (slot, x) in v.iter_mut().zip(m.iter()) {
                *slot -= c * x;
            }
        }
    }

    /// `<n_j, v>_{M_u} / |v|_{M_u}` for each mode: how much of `v` lies in
    /// each mode direction, relative to its own size.
    pub fn components(&self, sys: &SaddleSystem, v: &[f64]) -> Vec<f64> {
        let nv = sys.norm_u(v);
        if !nv.is_finite() || nv <= 0.0 {
            return vec![0.0; self.modes.len()];
        }
        self.modes.iter().map(|m| sys.dot_u(m, v) / nv).collect()
    }
}

/// The measured divergence leak of a set of `M_u`-normalised velocity modes.
///
/// `delta_j = h_ref |D n_j|_{M_p} / |n_j|_{M_u}` is dimensionless and is the
/// number to compare across rungs: for a GENERIC velocity field `D u ~ u / h`,
/// so `delta` is `O(1)` for a field that is not close to divergence free and
/// `O(h^2)` for the rotations. The raw `|D n_j|_{M_p}` is carried too, because
/// the residual floor of section 4 is compared against an absolute quantity.
#[derive(Debug, Clone)]
pub struct DivergenceLeak {
    pub which: Deflation,
    /// `|D n_j|_{M_p}` with `|n_j|_{M_u} = 1`, units 1/m in the sense above.
    pub absolute: Vec<f64>,
    /// `h_ref` times the above: dimensionless.
    pub dimensionless: Vec<f64>,
    /// Fraction of `|D n_j|_{M_p}^2` contributed by cells that touch a
    /// boundary facet. For the rotations this is exactly 1 by the identity in
    /// the module docs, and it is measured rather than assumed.
    pub boundary_share: Vec<f64>,
}

/// `|D n_j|` for each mode of a basis, with the boundary-cell share.
pub fn divergence_leak(
    sys: &SaddleSystem,
    basis: &ModeBasis,
) -> Result<DivergenceLeak, MantleError> {
    let topo = sys.topology();
    let mut boundary_cell = vec![false; sys.cell_count()];
    for face in topo.faces().iter() {
        if face.neighbor.is_none() {
            boundary_cell[face.owner as usize] = true;
        }
    }
    let mut absolute = Vec::with_capacity(basis.dim());
    let mut dimensionless = Vec::with_capacity(basis.dim());
    let mut boundary_share = Vec::with_capacity(basis.dim());
    for j in 0..basis.dim() {
        let d = sys.divergence(basis.mode(j))?;
        let total = sys.dot_p(&d, &d);
        let mut inside = 0.0f64;
        for (c, value) in d.iter().enumerate() {
            if boundary_cell[c] {
                inside += sys.m_p()[c] * value * value;
            }
        }
        let n = total.sqrt();
        absolute.push(n);
        dimensionless.push(sys.mesh_scale() * n);
        boundary_share.push(if total > 0.0 { inside / total } else { 0.0 });
    }
    Ok(DivergenceLeak {
        which: basis.which(),
        absolute,
        dimensionless,
        boundary_share,
    })
}

// ---------------------------------------------------------------------
// The block preconditioner

/// Whether the three near-null rotations are SUPPLIED to the preconditioner
/// as a coarse space.
///
/// This is the alternative to deleting them. The `A` block cannot produce a
/// rotation component at all (`K r = 0`, so a deflated CG returns exactly
/// zero there), which means that without a coarse space the Krylov velocity
/// iterate lives permanently in `null(K)^perp` and the three rotation-tested
/// momentum rows are reachable only through the pressure. With a coarse
/// space, the rotation amplitudes are chosen by a 3x3 least-squares fit
/// against the CONTINUITY residual - which is the only equation that sees
/// them - and the preconditioner can produce them directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoarseSpace {
    /// Nothing. The Krylov velocity then lives permanently in
    /// `null(K)^perp_{M_u}`, i.e. the rotation gauge is implicitly fixed at
    /// zero and the truncated translations are unreachable.
    None,
    /// The three near-null rotations.
    Rotations,
    /// The three truncated translations. These are NOT near-null for the
    /// coupled system - `D t` is `O(1)` - so a coarse space for them is a
    /// straightforwardly well-conditioned 3x3 fit, and the preconditioner
    /// needs one because the `A` block cannot produce them either.
    Translations,
    /// All six.
    Six,
}

impl CoarseSpace {
    pub fn label(self) -> &'static str {
        match self {
            CoarseSpace::None => "none",
            CoarseSpace::Rotations => "rotations(3)",
            CoarseSpace::Translations => "translations(3)",
            CoarseSpace::Six => "six",
        }
    }

    /// Indices into [`NullSpace::raw_modes`].
    pub fn indices(self) -> &'static [usize] {
        match self {
            CoarseSpace::None => &[],
            CoarseSpace::Rotations => &[3, 4, 5],
            CoarseSpace::Translations => &[0, 1, 2],
            CoarseSpace::Six => &[0, 1, 2, 3, 4, 5],
        }
    }
}

/// The near-nullspace coarse space, and the evidence for how ill conditioned
/// it is.
///
/// The `A` block returns exactly zero on `null(K)`, so a block preconditioner
/// built on it cannot produce ANY of the six modes. The only equation that
/// sees them is continuity, so the amplitudes are chosen by a least-squares
/// fit of `D(sum c_j n_j)` against the continuity residual. For the three
/// truncated translations that fit is well conditioned (`D t` is `O(1)`); for
/// the three rotations it divides by `|D r|^2 = O(h^4)`, which is exactly why
/// [`Self::eigenvalues`] is reported rather than assumed.
///
/// This is the treatment the packet's amendment ALLOWS for a near-null mode:
/// the rotations are SUPPLIED to the preconditioner, not deleted from the
/// equations.
#[derive(Debug, Clone)]
pub struct NullCoarse {
    which: CoarseSpace,
    modes: Vec<Vec<f64>>,
    div: Vec<Vec<f64>>,
    gram: Vec<Vec<f64>>,
    eigenvalues: Vec<f64>,
}

impl NullCoarse {
    #[allow(clippy::needless_range_loop)]
    pub fn build(
        topo: &ShellTopology,
        sys: &SaddleSystem,
        which: CoarseSpace,
    ) -> Result<Self, MantleError> {
        let basis = ModeBasis::from_indices(topo, sys, Deflation::Six, which.indices())?;
        let modes: Vec<Vec<f64>> = (0..basis.dim()).map(|j| basis.mode(j).to_vec()).collect();
        let mut div = Vec::with_capacity(modes.len());
        for m in modes.iter() {
            div.push(sys.divergence(m)?);
        }
        let d = modes.len();
        let mut gram = vec![vec![0.0f64; d]; d];
        for j in 0..d {
            for k in 0..d {
                gram[j][k] = sys.dot_p(&div[j], &div[k]);
            }
        }
        let h2 = sys.mesh_scale() * sys.mesh_scale();
        let scaled: Vec<Vec<f64>> = gram
            .iter()
            .map(|row| row.iter().map(|v| v * h2).collect())
            .collect();
        let eigenvalues = if d == 0 {
            Vec::new()
        } else {
            jacobi_symmetric(scaled).0
        };
        Ok(Self {
            which,
            modes,
            div,
            gram,
            eigenvalues,
        })
    }

    pub fn which(&self) -> CoarseSpace {
        self.which
    }

    pub fn dim(&self) -> usize {
        self.modes.len()
    }

    /// Eigenvalues of `h_ref^2 <D n_j, D n_k>_{M_p}`, ascending and
    /// dimensionless. For the rotations these are `(h |D r|)^2`, i.e. the
    /// square of the leak measured in [`divergence_leak`].
    pub fn eigenvalues(&self) -> &[f64] {
        &self.eigenvalues
    }

    pub fn gram(&self) -> &[Vec<f64>] {
        &self.gram
    }

    /// Add `sum_j c_j n_j` to `du`, with `c` the least-squares minimiser of
    /// `|rp - D du - D(sum c_j n_j)|_{M_p}`. Returns the amplitudes.
    fn augment(
        &self,
        sys: &SaddleSystem,
        du: &mut [f64],
        rp: &[f64],
    ) -> Result<Vec<f64>, MantleError> {
        if self.modes.is_empty() {
            return Ok(Vec::new());
        }
        let d_du = sys.divergence(du)?;
        let mut rhs = vec![0.0f64; self.modes.len()];
        for (j, slot) in rhs.iter_mut().enumerate() {
            let mut acc = 0.0;
            for (c, (a, b)) in rp.iter().zip(d_du.iter()).enumerate() {
                acc += sys.m_p()[c] * self.div[j][c] * (a - b);
            }
            *slot = acc;
        }
        let coeff = solve_symmetric(&self.gram, &rhs)?;
        for (j, m) in self.modes.iter().enumerate() {
            if coeff[j] == 0.0 {
                continue;
            }
            for (slot, x) in du.iter_mut().zip(m.iter()) {
                *slot += coeff[j] * x;
            }
        }
        Ok(coeff)
    }
}

/// Tolerances and structure of the block preconditioner. Every field is a
/// stated choice.
#[derive(Debug, Clone, Copy)]
pub struct PreconditionerOptions {
    /// Relative tolerance of the inner `A`-block CG, in the inner metric.
    /// Deliberately LOOSE - this is a preconditioner, and FGMRES is the
    /// outer method precisely so that it may be approximate and variable.
    pub inner_tol: f64,
    pub inner_max_iters: usize,
    /// Which inner product the inner CG runs in.
    /// [`Metric::StiffnessDiagonal`] is Jacobi-preconditioned CG, measured in
    /// the `A`-block packet at a constant `1.6x` speedup over `M_u`.
    pub inner_metric: Metric,
    /// `S~^-1 = -alpha I` with `alpha = eta / c_schur`.
    pub c_schur: f64,
    pub coarse: CoarseSpace,
    /// Re-project the inner solve's output onto `null(K)^perp_{M_u}`.
    ///
    /// **Default true, and turning it off is a NEGATIVE CONTROL, not an
    /// option.** `deflated_cg` returns a vector orthogonal to the six modes in
    /// its OWN metric; with [`Metric::StiffnessDiagonal`] that is not
    /// `M_u`-orthogonality, and the residue lands in the near-null rotation
    /// directions of the coupled system, where the outer Krylov method
    /// amplifies it. Measured with it off, in
    /// [`tests::the_solve_does_not_depend_on_the_viscosity_scale`] and
    /// [`tests::the_manufactured_saddle_pair_is_recovered`].
    ///
    /// This changes the PRECONDITIONER only. No equation is altered and the
    /// residual is always measured against the full untouched system.
    pub reproject_inner_in_m_u: bool,
}

impl Default for PreconditionerOptions {
    fn default() -> Self {
        Self {
            inner_tol: 1e-2,
            inner_max_iters: 4_000,
            inner_metric: Metric::StiffnessDiagonal,
            // The TEXTBOOK value, taken literally and not fitted: the
            // standard result is `B_fe K^-1 B_fe^T ~ (1/eta) M_p`, which in
            // this operator scaling is `S ~ (1/eta) I`, i.e. c_schur = 1.
            // It is used at EVERY rung, so no per-resolution tuning enters,
            // and it was fixed before any residual curve was looked at.
            //
            // What it costs is measured rather than hoped for. The Schur
            // spectrum this approximates is NOT bounded on this operator pair
            // (`tests::the_schur_complement_is_a_scaled_identity_and_its_low_end_is_identified`
            // measures eta*lambda_max at 3.06, 46.4, 410, 657 over n = 2, 4,
            // 8, 16), so a scalar cannot follow it and the outer count grows.
            // That is the honest MC0 result, and the fix - a weighted BFBT
            // Schur - is explicitly out of scope for this packet.
            // `tests::the_schur_scaling_is_not_a_tuning_knob` sweeps c_schur
            // over two decades so the value used is visible as a choice.
            c_schur: 1.0,
            coarse: CoarseSpace::None,
            reproject_inner_in_m_u: true,
        }
    }
}

/// Upper block triangular preconditioner
///
/// ```text
/// P   = [ B~  G_imp ]        P^-1 (a, b) :  y = S~^-1 b
///       [ 0   S~    ]                       x = B~^-1 (a - G_imp y)
/// ```
///
/// **Why upper and not lower.** With EXACT blocks (`B~ = B` on `range(B)`,
/// `S~ = S = D B^-1 G`) the right-preconditioned operator is
///
/// ```text
/// M P^-1 = [ I        0  ]        and    (M P^-1)^2 = I
///          [ D B^-1  -I  ]
/// ```
///
/// so its minimal polynomial is `(t-1)(t+1)` and GMRES terminates in TWO
/// iterations. The lower triangular arrangement gives
/// `[I + G S^-1 D B^-1, -G S^-1; D B^-1, 0]`, which has no such two-term
/// annihilating polynomial for this system, and needs the same `B` solve. The
/// upper form additionally needs only ONE `A`-block solve and one gradient
/// application per preconditioner call.
/// [`tests::the_exact_block_preconditioner_terminates_in_two_iterations`]
/// measures the two-iteration property directly, which is the check that the
/// blocks are wired the way this comment claims.
///
/// **The Schur block.** For the finite-element pair the standard result is
/// `B_fe K^-1 B_fe^T ~ (1/eta) M_p` with `B_fe = M_p D` the divergence tested
/// against the pressure basis. In the `M_u^-1`-scaled operator form used
/// here, `S = D K^-1 D^T M_p = M_p^-1 (M_p D K^-1 D^T M_p) ~ (1/eta) I`:
/// **the pressure mass matrix cancels and the viscosity-weighted pressure-mass
/// Schur approximation degenerates to a SCALAR.** That is a property of this
/// scaling and of isoviscosity, not a simplification; MC2's variable
/// viscosity puts `eta_c` back inside and the scalar becomes a diagonal.
#[derive(Debug)]
pub struct BlockPreconditioner<'a> {
    block: ABlock<'a>,
    null: NullSpace,
    metric_diag: Vec<f64>,
    opts: PreconditionerOptions,
    /// All six `A`-null modes, `M_u`-orthonormal. Used to re-project the
    /// inner solve's output, NOT to alter any equation.
    six: ModeBasis,
    coarse: Option<NullCoarse>,
    /// The EXACT Schur block, assembled densely. Only ever used at the two
    /// smallest rungs, and only to verify that the block algebra this
    /// module's doc comment claims is the algebra it implements.
    dense_schur: Option<DenseSchur>,
    /// Cumulative inner CG iterations, so an outer count is never reported
    /// without the Krylov work behind it (MC-Q5).
    inner_iterations: usize,
    inner_calls: usize,
    inner_max: usize,
    inner_not_converged: usize,
    worst_coarse_amplitude: f64,
    worst_inner_leak: f64,
}

impl<'a> BlockPreconditioner<'a> {
    pub fn build(
        topo: &ShellTopology,
        op: &'a ViscousOperator,
        sys: &SaddleSystem,
        opts: PreconditionerOptions,
    ) -> Result<Self, MantleError> {
        if !opts.c_schur.is_finite() || opts.c_schur <= 0.0 {
            return Err(MantleError::field(
                "Schur scaling c_schur",
                opts.c_schur,
                "the Schur approximation must be a positive multiple of the identity",
            ));
        }
        let block = ABlock::new(topo, op, opts.inner_metric)?;
        let null = NullSpace::build(topo, &block)?;
        let metric_diag = match opts.inner_metric {
            Metric::FaceMass => op.face_mass().to_vec(),
            Metric::StiffnessDiagonal => stiffness_diagonal(op),
        };
        let six = ModeBasis::build(topo, sys, Deflation::Six)?;
        let coarse = match opts.coarse {
            CoarseSpace::None => None,
            other => Some(NullCoarse::build(topo, sys, other)?),
        };
        Ok(Self {
            block,
            null,
            metric_diag,
            opts,
            six,
            coarse,
            dense_schur: None,
            inner_iterations: 0,
            inner_calls: 0,
            inner_max: 0,
            inner_not_converged: 0,
            worst_coarse_amplitude: 0.0,
            worst_inner_leak: 0.0,
        })
    }

    pub fn options(&self) -> PreconditionerOptions {
        self.opts
    }

    pub fn inner_iterations(&self) -> usize {
        self.inner_iterations
    }

    pub fn inner_calls(&self) -> usize {
        self.inner_calls
    }

    pub fn inner_max(&self) -> usize {
        self.inner_max
    }

    pub fn inner_not_converged(&self) -> usize {
        self.inner_not_converged
    }

    pub fn coarse_space(&self) -> Option<&NullCoarse> {
        self.coarse.as_ref()
    }

    /// Worst `M_u` null leak observed in an inner solve's raw output, before
    /// the re-projection. This is the inner-product trap of the `A`-block
    /// findings, measured live: a `diag(K)`-deflated solve is NOT
    /// `M_u`-orthogonal to the six modes, and the leftover lands in the
    /// near-null rotation directions where the coupled system amplifies it.
    pub fn worst_inner_leak(&self) -> f64 {
        self.worst_inner_leak
    }

    pub fn worst_coarse_amplitude(&self) -> f64 {
        self.worst_coarse_amplitude
    }

    pub fn reset_counters(&mut self) {
        self.inner_iterations = 0;
        self.inner_calls = 0;
        self.inner_max = 0;
        self.inner_not_converged = 0;
        self.worst_coarse_amplitude = 0.0;
        self.worst_inner_leak = 0.0;
    }

    /// `B~^-1 a`: solve `K x = M_u a` with the deflated CG of
    /// [`super::stokes_solve`], in this preconditioner's metric.
    ///
    /// The rhs conversion is the one place a metric mistake could hide, so it
    /// is written out: `deflated_cg` solves `P^-1 K x = b`, therefore `b` must
    /// be `P^-1 M_u a`. The consistency of that rhs is metric independent -
    /// `<n_j, P^-1 M_u a>_P = n_j^T M_u a = <n_j, a>_{M_u}` - so a rhs that is
    /// `M_u`-orthogonal to the null space is in the range whichever metric the
    /// inner solve runs in. `tests::the_inner_solve_is_metric_consistent`
    /// measures that.
    pub fn solve_a_block(
        &mut self,
        sys: &SaddleSystem,
        a: &[f64],
    ) -> Result<Vec<f64>, MantleError> {
        let mut b = vec![0.0f64; a.len()];
        for (f, slot) in b.iter_mut().enumerate() {
            let d = self.metric_diag[f];
            if d > 0.0 {
                *slot = sys.m_u()[f] * a[f] / d;
            }
        }
        self.null.project(&self.block, &mut b);
        if self.block.norm(&b) <= 0.0 {
            // The whole rhs was inside null(K). There is nothing for the A
            // block to do and deflated_cg would (correctly) refuse; returning
            // zero is the right answer, not a fallback.
            self.inner_calls += 1;
            return Ok(vec![0.0f64; a.len()]);
        }
        let cg = CgOptions {
            tol: self.opts.inner_tol,
            max_iters: self.opts.inner_max_iters,
            ..CgOptions::default()
        };
        let (mut x, report) = deflated_cg(&self.block, &self.null, &b, cg)?;
        // Re-project in M_u. `deflated_cg` returns a vector orthogonal to the
        // six modes in the INNER metric; with Metric::StiffnessDiagonal that
        // is not M_u-orthogonality, and the leftover component sits in the
        // near-null rotation directions of the coupled system, where GMRES
        // will happily amplify it by orders of magnitude to shave a little
        // residual. Measured at 1.8e-2 relative in the A-block findings'
        // negative control; measured live here as `worst_inner_leak`.
        //
        // This changes the PRECONDITIONER, not any equation: the residual is
        // still measured against the full untouched system.
        let leak = self.six.components(sys, &x);
        self.worst_inner_leak = leak
            .iter()
            .fold(self.worst_inner_leak, |acc, c| acc.max(c.abs()));
        if self.opts.reproject_inner_in_m_u {
            self.six.project(sys, &mut x);
        }
        self.inner_iterations += report.iterations;
        self.inner_calls += 1;
        self.inner_max = self.inner_max.max(report.iterations);
        if !report.converged {
            self.inner_not_converged += 1;
        }
        Ok(x)
    }

    /// `P^-1 (a, b)`.
    pub fn apply(
        &mut self,
        sys: &SaddleSystem,
        a: &[f64],
        b: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), MantleError> {
        let mut dp = match self.dense_schur.as_ref() {
            Some(dense) => dense.apply(sys, b),
            None => {
                let alpha = sys.eta() / self.opts.c_schur;
                b.iter().map(|v| -alpha * v).collect()
            }
        };
        sys.remove_mean_pressure(&mut dp)?;
        let g = sys.gradient_implicit(&dp)?;
        let rhs: Vec<f64> = a.iter().zip(g.iter()).map(|(x, y)| x - y).collect();
        let mut du = self.solve_a_block(sys, &rhs)?;
        // `as_ref` and a deferred write, not a clone: RotationCoarse holds
        // six full-length vectors and this runs once per Krylov step.
        let mut seen: Option<f64> = None;
        if let Some(rot) = self.coarse.as_ref() {
            let coeff = rot.augment(sys, &mut du, b)?;
            let scale = sys.norm_u(&du);
            if scale > 0.0 {
                seen = Some(coeff.iter().fold(0.0f64, |acc, c| acc.max(c.abs())) / scale);
            }
        }
        if let Some(worst) = seen {
            self.worst_coarse_amplitude = self.worst_coarse_amplitude.max(worst);
        }
        Ok((du, dp))
    }
}

/// The exact Schur complement `S_0 = D B^+ G_imp`, assembled column by column.
///
/// One `A`-block solve per pressure cell, so this is affordable only at the
/// smallest rungs and is used for exactly one purpose: proving that the
/// upper block triangular preconditioner of [`BlockPreconditioner`] really
/// does give `(M P^-1)^2 = I` and a two-step GMRES. If the wiring were wrong
/// - a sign, a transpose, a metric - the two-step property would fail and no
///   iteration table elsewhere in this module would say so.
///
/// **It is self-adjoint only when the inner solve's deflation metric is
/// `M_u`.** `M_u G q = -D^T M_p q` has an `O(1)` component inside `null(K)`
/// (through the truncated translations, whose `D t` is `O(1)`), so the
/// `A`-block right-hand side is genuinely INCONSISTENT in three directions and
/// which pseudo-inverse is used changes the answer at `O(1)`, not at roundoff.
/// Measured: with the `diag(K)` deflation plus an `M_u` re-projection the three
/// translation-aligned Schur eigenvalues move from `+5.5e-4` to `-3.9e-3` and
/// the matrix stops being positive semidefinite. Both are honest descriptions
/// of a particular PRECONDITIONER; only the `M_u` one is the `V_0`-restricted
/// Schur complement, so that is the one used wherever a spectrum is reported.
///
/// Stored as the eigendecomposition of `M_p^{1/2} (-S_0) M_p^{-1/2}`, which is
/// symmetric because `-S_0` is then self-adjoint in `M_p`; the constant pressure
/// mode is dropped by a RELATIVE eigenvalue cut, and how many modes were
/// dropped is recorded so a second near-null pressure mode could not be
/// silently discarded along with it.
#[derive(Debug, Clone)]
pub struct DenseSchur {
    vectors: Vec<Vec<f64>>,
    inverse_values: Vec<f64>,
    /// Eigenvalues of `eta * M_p^{1/2} (-S_0) M_p^{-1/2}`, ascending.
    pub eigenvalues: Vec<f64>,
    pub dropped: usize,
    /// `max |S_ij - S_ji|` after the symmetrising similarity, relative to the
    /// largest entry: the measured self-adjointness of the Schur complement.
    pub asymmetry: f64,
}

impl DenseSchur {
    /// Assemble it. `cutoff` is RELATIVE to the largest eigenvalue.
    #[allow(clippy::needless_range_loop)]
    pub fn build(
        sys: &SaddleSystem,
        precon: &mut BlockPreconditioner,
        cutoff: f64,
    ) -> Result<Self, MantleError> {
        let n = sys.cell_count();
        let root: Vec<f64> = sys.m_p().iter().map(|v| v.sqrt()).collect();
        let mut a = vec![vec![0.0f64; n]; n];
        for j in 0..n {
            let mut e = vec![0.0f64; n];
            e[j] = 1.0;
            let col = schur_apply(sys, precon, &e)?;
            for (i, row) in a.iter_mut().enumerate() {
                row[j] = root[i] * col[i] / root[j];
            }
        }
        let scale = a
            .iter()
            .flat_map(|r| r.iter())
            .fold(0.0f64, |m, v| m.max(v.abs()));
        let mut asymmetry = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                asymmetry = asymmetry.max((a[i][j] - a[j][i]).abs());
            }
        }
        asymmetry = if scale > 0.0 { asymmetry / scale } else { 0.0 };
        // Symmetrise the roundoff before the eigensolver sees it; the size of
        // what was symmetrised away is reported above.
        for i in 0..n {
            for j in (i + 1)..n {
                let m = 0.5 * (a[i][j] + a[j][i]);
                a[i][j] = m;
                a[j][i] = m;
            }
        }
        let (values, vectors) = jacobi_symmetric(a);
        let top = values.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let mut inverse_values = vec![0.0f64; n];
        let mut dropped = 0usize;
        for (i, v) in values.iter().enumerate() {
            if v.abs() <= cutoff * top {
                dropped += 1;
            } else {
                inverse_values[i] = 1.0 / v;
            }
        }
        Ok(Self {
            vectors,
            inverse_values,
            eigenvalues: values.iter().map(|v| sys.eta() * v).collect(),
            dropped,
            asymmetry,
        })
    }

    /// `vectors[row][index]`: the `row`-th component of the `index`-th
    /// eigenvector, in the `M_p^{1/2}`-symmetrised basis.
    pub fn eigenvector(&self, row: usize, index: usize) -> f64 {
        self.vectors[row][index]
    }

    /// `S_0^+ b = -(-S_0)^+ b`, the exact Schur solve.
    fn apply(&self, sys: &SaddleSystem, b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let root: Vec<f64> = sys.m_p().iter().map(|v| v.sqrt()).collect();
        let scaled: Vec<f64> = (0..n).map(|i| root[i] * b[i]).collect();
        let mut out = vec![0.0f64; n];
        for (i, inv) in self.inverse_values.iter().enumerate() {
            if *inv == 0.0 {
                continue;
            }
            let dot: f64 = (0..n).map(|r| self.vectors[r][i] * scaled[r]).sum();
            let c = dot * inv;
            for (r, slot) in out.iter_mut().enumerate() {
                *slot += c * self.vectors[r][i];
            }
        }
        // The minus is the sign of S_0 itself, which is negative definite.
        (0..n).map(|i| -out[i] / root[i]).collect()
    }
}

impl BlockPreconditioner<'_> {
    /// Replace the scalar Schur approximation with the exact dense one.
    pub fn set_dense_schur(&mut self, dense: DenseSchur) {
        self.dense_schur = Some(dense);
    }

    pub fn dense_schur(&self) -> Option<&DenseSchur> {
        self.dense_schur.as_ref()
    }
}

// ---------------------------------------------------------------------
// FGMRES

/// Tolerances and guards of one FGMRES run.
#[derive(Debug, Clone, Copy)]
pub struct FgmresOptions {
    /// Stop when the TRUE relative residual in the weighted norm of
    /// [`SaddleSystem::residual_norm`] falls below this.
    pub tol: f64,
    /// Total Arnoldi steps across all restart cycles.
    pub max_iterations: usize,
    /// Krylov dimension per cycle. Storage is `2 * restart` full vectors
    /// (the Arnoldi basis and the flexible preconditioned basis), which is
    /// why it is a stated number and not `max_iterations`.
    pub restart: usize,
    /// Declare a stall if the relative residual has not fallen by
    /// `stall_factor` over this many Arnoldi steps. A stall is REPORTED, not
    /// converted into a smaller tolerance.
    pub stall_window: usize,
    pub stall_factor: f64,
}

impl Default for FgmresOptions {
    fn default() -> Self {
        Self {
            // 1e-8, the same target the A-block CG ladder used, chosen before
            // any residual curve was looked at and not moved afterwards. Where
            // a run does not reach it, the run is reported as not reaching it.
            tol: 1e-8,
            max_iterations: 400,
            restart: 60,
            stall_window: 40,
            stall_factor: 1.05,
        }
    }
}

/// What one FGMRES run did, with the inner Krylov work it cost.
///
/// MC-Q5 asks for outer counts reported ALONGSIDE Krylov counts, so
/// [`Self::outer_iterations`] and [`Self::inner_cg_iterations`] are separate
/// fields and no method collapses them.
#[derive(Debug, Clone)]
pub struct FgmresReport {
    pub outer_iterations: usize,
    pub restart_cycles: usize,
    pub inner_cg_iterations: usize,
    pub inner_cg_calls: usize,
    pub inner_cg_max: usize,
    pub inner_cg_not_converged: usize,
    pub converged: bool,
    pub stalled: bool,
    /// From the Arnoldi recurrence at the last step.
    pub rel_residual: f64,
    /// Recomputed from scratch as `|b - M x|` in the same norm. A gap between
    /// this and [`Self::rel_residual`] is Arnoldi drift and is reported
    /// rather than hidden by returning only one of them.
    pub true_rel_residual: f64,
    /// Separate halves of the true final residual, each relative to `|b|`:
    /// a converged momentum row with a stalled continuity row is a different
    /// failure from the reverse and one number cannot say which happened.
    pub momentum_rel: f64,
    pub continuity_rel: f64,
    /// First step at which the relative residual fell below
    /// `1e-2, 1e-4, 1e-6, 1e-8`.
    pub milestones: [Option<usize>; 4],
    /// `|<n_j, r_u>_{M_u}| / |b|` for the six `A`-null modes at the end:
    /// three truncated translations then three rotations. This is the
    /// task-4 measurement - which part of a residual floor lives in a mode
    /// the `A` block cannot see.
    pub null_residual: [f64; NULL_DIM],
    pub history: Vec<f64>,
}

impl FgmresReport {
    /// Largest final residual component among the three ROTATION directions.
    pub fn rotation_residual(&self) -> f64 {
        self.null_residual[3..6]
            .iter()
            .fold(0.0f64, |a, b| a.max(*b))
    }

    /// Largest final residual component among the three TRUNCATED TRANSLATION
    /// directions.
    pub fn translation_residual(&self) -> f64 {
        self.null_residual[0..3]
            .iter()
            .fold(0.0f64, |a, b| a.max(*b))
    }
}

/// Right-preconditioned FGMRES for the coupled saddle point.
///
/// FGMRES and not MINRES, and the reason is structural rather than a
/// preference: the preconditioner's `A` block is an inexact CG solve whose
/// iteration count varies from call to call, so the preconditioner is a
/// DIFFERENT linear operator at each Krylov step. A short-recurrence method
/// has no way to represent that; the flexible variant stores the
/// preconditioned basis `z_j` explicitly and is the standard answer.
/// [`GradientMode::Corrected`] additionally destroys the symmetry outright,
/// since the deferred correction is not the adjoint of anything.
///
/// Modified Gram-Schmidt run TWICE per vector. One pass is not enough once
/// the preconditioner is inexact, and the cost is one extra pass over `j`
/// vectors against an operator application that costs far more.
#[allow(clippy::too_many_arguments)]
pub fn fgmres_saddle(
    sys: &SaddleSystem,
    precon: &mut BlockPreconditioner,
    modes: &ModeBasis,
    six: &ModeBasis,
    grad: GradientMode,
    rhs_u: &[f64],
    rhs_p: &[f64],
    opts: FgmresOptions,
) -> Result<(Vec<f64>, Vec<f64>, FgmresReport), MantleError> {
    if opts.restart == 0 || opts.max_iterations == 0 {
        return Err(MantleError::field(
            "FGMRES restart / max_iterations",
            opts.restart as f64,
            "both must be strictly positive",
        ));
    }
    let mut b_u = rhs_u.to_vec();
    let mut b_p = rhs_p.to_vec();
    modes.project(sys, &mut b_u);
    sys.remove_mean_pressure(&mut b_p)?;
    let b_norm = sys.residual_norm(&b_u, &b_p);
    if !b_norm.is_finite() || b_norm <= 0.0 {
        return Err(MantleError::field(
            "FGMRES right-hand side norm",
            b_norm,
            "the right-hand side is empty after the pressure-mean projection",
        ));
    }

    let mut x_u = sys.zeros_u();
    let mut x_p = sys.zeros_p();
    let mut history: Vec<f64> = Vec::new();
    let mut milestones: [Option<usize>; 4] = [None; 4];
    let thresholds = [1e-2f64, 1e-4, 1e-6, 1e-8];
    let mut total = 0usize;
    let mut cycles = 0usize;
    let mut converged = false;
    let mut stalled = false;
    let mut rel = 1.0f64;
    let mut window_ref = 1.0f64;
    let mut window_start = 0usize;

    'outer: while total < opts.max_iterations {
        let (m_u, m_p) = sys.apply(modes, grad, &x_u, &x_p)?;
        let r_u: Vec<f64> = b_u.iter().zip(m_u.iter()).map(|(a, c)| a - c).collect();
        let r_p: Vec<f64> = b_p.iter().zip(m_p.iter()).map(|(a, c)| a - c).collect();
        let beta = sys.residual_norm(&r_u, &r_p);
        rel = beta / b_norm;
        if cycles == 0 {
            history.push(rel);
        }
        if rel <= opts.tol {
            converged = true;
            break;
        }
        cycles += 1;

        let inv = 1.0 / beta;
        let mut v_u: Vec<Vec<f64>> = vec![r_u.iter().map(|x| x * inv).collect()];
        let mut v_p: Vec<Vec<f64>> = vec![r_p.iter().map(|x| x * inv).collect()];
        let mut z_u: Vec<Vec<f64>> = Vec::new();
        let mut z_p: Vec<Vec<f64>> = Vec::new();
        let m = opts.restart.min(opts.max_iterations - total);
        let mut h = vec![vec![0.0f64; m]; m + 1];
        let mut cs = vec![0.0f64; m];
        let mut sn = vec![0.0f64; m];
        let mut g = vec![0.0f64; m + 1];
        g[0] = beta;
        let mut used = 0usize;

        for j in 0..m {
            let (mut w_u, mut w_p) = precon.apply(sys, &v_u[j], &v_p[j])?;
            z_u.push(std::mem::take(&mut w_u));
            z_p.push(std::mem::take(&mut w_p));
            let (mut a_u, mut a_p) = sys.apply(modes, grad, &z_u[j], &z_p[j])?;
            for _pass in 0..2 {
                for i in 0..=j {
                    let c = arnoldi_dot(sys, &a_u, &a_p, &v_u[i], &v_p[i]);
                    h[i][j] += c;
                    for (slot, x) in a_u.iter_mut().zip(v_u[i].iter()) {
                        *slot -= c * x;
                    }
                    for (slot, x) in a_p.iter_mut().zip(v_p[i].iter()) {
                        *slot -= c * x;
                    }
                }
            }
            let hn = sys.residual_norm(&a_u, &a_p);
            h[j + 1][j] = hn;
            // Grow the Arnoldi basis BEFORE the Givens sweep overwrites
            // h[j+1][j]: `hn` is captured, but the basis vector itself has to
            // exist for step j+1 to have a `v_{j+1}` to start from.
            if hn > 0.0 {
                let s = 1.0 / hn;
                v_u.push(a_u.iter().map(|x| x * s).collect());
                v_p.push(a_p.iter().map(|x| x * s).collect());
            }
            total += 1;
            used = j + 1;
            // Givens: apply the accumulated rotations to the new column, then
            // eliminate the subdiagonal entry.
            for i in 0..j {
                let t = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = t;
            }
            let denom = (h[j][j] * h[j][j] + h[j + 1][j] * h[j + 1][j]).sqrt();
            if !denom.is_finite() {
                return Err(MantleError::field(
                    "FGMRES Givens denominator",
                    denom,
                    "the Arnoldi recurrence produced a non-finite entry",
                ));
            }
            if denom == 0.0 {
                cs[j] = 1.0;
                sn[j] = 0.0;
            } else {
                cs[j] = h[j][j] / denom;
                sn[j] = h[j + 1][j] / denom;
            }
            h[j][j] = cs[j] * h[j][j] + sn[j] * h[j + 1][j];
            h[j + 1][j] = 0.0;
            g[j + 1] = -sn[j] * g[j];
            g[j] *= cs[j];
            rel = g[j + 1].abs() / b_norm;
            history.push(rel);
            for (slot, t) in milestones.iter_mut().zip(thresholds.iter()) {
                if slot.is_none() && rel <= *t {
                    *slot = Some(total);
                }
            }
            if rel <= opts.tol {
                converged = true;
            }
            if total - window_start >= opts.stall_window {
                if window_ref / rel.max(f64::MIN_POSITIVE) < opts.stall_factor {
                    stalled = true;
                }
                window_ref = rel;
                window_start = total;
            }
            // A happy breakdown (hn == 0) means the Krylov space is
            // invariant; continuing would divide by zero.
            if converged || stalled || hn <= 0.0 || total >= opts.max_iterations {
                break;
            }
        }

        // Back-substitute the small triangular system and update x.
        let mut y = vec![0.0f64; used];
        for i in (0..used).rev() {
            let mut acc = g[i];
            for k in (i + 1)..used {
                acc -= h[i][k] * y[k];
            }
            y[i] = if h[i][i] == 0.0 { 0.0 } else { acc / h[i][i] };
        }
        for (i, yi) in y.iter().enumerate() {
            for (slot, x) in x_u.iter_mut().zip(z_u[i].iter()) {
                *slot += yi * x;
            }
            for (slot, x) in x_p.iter_mut().zip(z_p[i].iter()) {
                *slot += yi * x;
            }
        }
        sys.remove_mean_pressure(&mut x_p)?;
        if converged || stalled {
            break 'outer;
        }
    }

    // The residual reported is recomputed from scratch, never the recurrence.
    let (m_u, m_p) = sys.apply(modes, grad, &x_u, &x_p)?;
    let r_u: Vec<f64> = b_u.iter().zip(m_u.iter()).map(|(a, c)| a - c).collect();
    let r_p: Vec<f64> = b_p.iter().zip(m_p.iter()).map(|(a, c)| a - c).collect();
    let true_rel = sys.residual_norm(&r_u, &r_p) / b_norm;
    let momentum_rel = sys.norm_u(&r_u) / b_norm;
    let continuity_rel = sys.theta() * sys.norm_p(&r_p) / b_norm;
    let mut null_residual = [0.0f64; NULL_DIM];
    for (j, slot) in null_residual.iter_mut().enumerate() {
        *slot = (sys.dot_u(six.mode(j), &r_u) / b_norm).abs();
    }

    let report = FgmresReport {
        outer_iterations: total,
        restart_cycles: cycles,
        inner_cg_iterations: precon.inner_iterations(),
        inner_cg_calls: precon.inner_calls(),
        inner_cg_max: precon.inner_max(),
        inner_cg_not_converged: precon.inner_not_converged(),
        converged: converged && true_rel <= opts.tol,
        stalled,
        rel_residual: rel,
        true_rel_residual: true_rel,
        momentum_rel,
        continuity_rel,
        milestones,
        null_residual,
        history,
    };
    Ok((x_u, x_p, report))
}

/// The FGMRES inner product: [`SaddleSystem::residual_norm`]'s bilinear form.
fn arnoldi_dot(sys: &SaddleSystem, a_u: &[f64], a_p: &[f64], b_u: &[f64], b_p: &[f64]) -> f64 {
    sys.dot_u(a_u, b_u) + sys.theta() * sys.theta() * sys.dot_p(a_p, b_p)
}

// ---------------------------------------------------------------------
// The deferred correction as an outer fixed point

/// The lagged / under-relaxed treatment of the deferred gradient correction.
#[derive(Debug, Clone, Copy)]
pub struct DeferredOptions {
    /// `p_lag <- (1 - w) p_lag + w p`. `w = 1` is the unrelaxed fixed point.
    pub relax: f64,
    pub max_sweeps: usize,
    /// Target on the FULL corrected residual, not on the implicit one the
    /// inner FGMRES actually minimises.
    pub tol: f64,
    /// Declare divergence when the full residual exceeds this multiple of the
    /// best value seen. Reported, never converted into a looser tolerance.
    pub divergence_factor: f64,
    pub inner: FgmresOptions,
}

impl Default for DeferredOptions {
    fn default() -> Self {
        Self {
            relax: 1.0,
            max_sweeps: 40,
            tol: 1e-8,
            divergence_factor: 1e3,
            inner: FgmresOptions {
                // The inner solve must be TIGHTER than the fixed point's own
                // target, or the measured contraction factor would be the
                // inner tolerance rather than the fixed point's spectral
                // radius, which is the quantity this exists to measure.
                tol: 1e-11,
                ..FgmresOptions::default()
            },
        }
    }
}

/// What the fixed point did. The contraction factor is the point of it.
#[derive(Debug, Clone)]
pub struct DeferredReport {
    pub relax: f64,
    pub sweeps: usize,
    pub converged: bool,
    pub diverged: bool,
    pub outer_iterations: usize,
    pub inner_cg_iterations: usize,
    /// Full-system relative residual after each sweep.
    pub history: Vec<f64>,
    /// Geometric mean of the last three residual ratios: the MEASURED
    /// spectral radius of the fixed point. Above 1 means it diverges, and no
    /// tolerance choice can rescue it.
    pub contraction: f64,
    pub final_rel_residual: f64,
}

/// Solve the corrected system by lagging `G_exp` and iterating.
///
/// This is the classical finite-volume treatment and it is here to be
/// MEASURED against [`GradientMode::Corrected`], where the same correction is
/// applied inside the Krylov operator and there is no fixed point at all.
/// Because everything in sight is linear, the sweep map
/// `p_lag -> p` is affine, so the ratio of successive residuals converges to
/// the spectral radius of `(1 - w) I + w T` and
/// [`DeferredReport::contraction`] measures it directly.
#[allow(clippy::too_many_arguments)]
pub fn deferred_correction_solve(
    sys: &SaddleSystem,
    precon: &mut BlockPreconditioner,
    modes: &ModeBasis,
    six: &ModeBasis,
    rhs_u: &[f64],
    rhs_p: &[f64],
    opts: DeferredOptions,
) -> Result<(Vec<f64>, Vec<f64>, DeferredReport), MantleError> {
    if !opts.relax.is_finite() || opts.relax <= 0.0 || opts.relax > 1.0 {
        return Err(MantleError::field(
            "deferred-correction relaxation",
            opts.relax,
            "under-relaxation must lie in (0, 1]",
        ));
    }
    let mut b_u = rhs_u.to_vec();
    let mut b_p = rhs_p.to_vec();
    modes.project(sys, &mut b_u);
    sys.remove_mean_pressure(&mut b_p)?;
    let b_norm = sys.residual_norm(&b_u, &b_p);

    let mut p_lag = sys.zeros_p();
    let mut x_u = sys.zeros_u();
    let mut x_p = sys.zeros_p();
    let mut history: Vec<f64> = Vec::new();
    let mut outer_iterations = 0usize;
    let mut converged = false;
    let mut diverged = false;
    let mut best = f64::MAX;
    let mut sweeps = 0usize;

    for _ in 0..opts.max_sweeps {
        let ex = sys.gradient_explicit(&p_lag)?;
        let source: Vec<f64> = b_u.iter().zip(ex.iter()).map(|(a, c)| a - c).collect();
        let (u, p, rep) = fgmres_saddle(
            sys,
            precon,
            modes,
            six,
            GradientMode::Implicit,
            &source,
            &b_p,
            opts.inner,
        )?;
        outer_iterations += rep.outer_iterations;
        sweeps += 1;
        // The residual of the FULL corrected system at the new iterate.
        let (m_u, m_p) = sys.apply(modes, GradientMode::Corrected, &u, &p)?;
        let r_u: Vec<f64> = b_u.iter().zip(m_u.iter()).map(|(a, c)| a - c).collect();
        let r_p: Vec<f64> = b_p.iter().zip(m_p.iter()).map(|(a, c)| a - c).collect();
        let rel = sys.residual_norm(&r_u, &r_p) / b_norm;
        history.push(rel);
        x_u = u;
        x_p = p;
        for (slot, x) in p_lag.iter_mut().zip(x_p.iter()) {
            *slot = (1.0 - opts.relax) * *slot + opts.relax * x;
        }
        best = best.min(rel);
        if rel <= opts.tol {
            converged = true;
            break;
        }
        if !rel.is_finite() || rel > best * opts.divergence_factor {
            diverged = true;
            break;
        }
    }

    let contraction = if history.len() >= 2 {
        let take = history.len().min(4);
        let tail = &history[history.len() - take..];
        let mut acc = 1.0f64;
        for w in tail.windows(2) {
            acc *= w[1] / w[0];
        }
        acc.powf(1.0 / (take - 1) as f64)
    } else {
        f64::NAN
    };

    let final_rel = *history.last().unwrap_or(&f64::NAN);
    Ok((
        x_u,
        x_p,
        DeferredReport {
            relax: opts.relax,
            sweeps,
            converged,
            diverged,
            outer_iterations,
            inner_cg_iterations: precon.inner_iterations(),
            history,
            contraction,
            final_rel_residual: final_rel,
        },
    ))
}

// ---------------------------------------------------------------------
// The Schur complement, measured rather than assumed

/// The measured spectrum of `eta * (-D B^+ G_imp)` in the `M_p` inner
/// product, on the mean-free pressure space.
///
/// `B^+` is the deflated pseudo-inverse on `null(K)^perp` - the same operator
/// the preconditioner's `A` block approximates - so this IS the Schur
/// complement the preconditioner is approximating, and not an idealisation of
/// it. Two things fall out of one sweep:
///
/// - `lambda_min` is the SQUARE of the discrete inf-sup constant of this
///   velocity/pressure pair. A checkerboard-type spurious pressure mode would
///   appear here as a `lambda_min` collapsing under refinement, and nothing
///   else in this repository would notice.
/// - `c_schur` for [`PreconditionerOptions`] is read off the middle of this
///   spectrum. Since it turns out to be `O(1)` and nearly `h`-independent,
///   ONE constant serves every rung.
#[derive(Debug, Clone)]
pub struct SchurReport {
    pub steps: usize,
    pub cells: usize,
    /// `eta * lambda`, dimensionless.
    pub lambda_min: f64,
    pub lambda_max: f64,
    /// Rigorous Ritz brackets `beta_m |e_m^T s_i| / lambda_i`.
    pub bracket_min: f64,
    pub bracket_max: f64,
    /// `|G_imp 1|_{M_u} / (|1|_{M_p} * scale)`: the constant pressure mode is
    /// an EXACT null mode of the gradient, and this measures it rather than
    /// assuming it.
    pub constant_mode_leak: f64,
    pub inner_cg_iterations: usize,
}

/// One deflated Lanczos sweep on the Schur complement.
pub fn schur_spectrum(
    sys: &SaddleSystem,
    precon: &mut BlockPreconditioner,
    steps: usize,
    seed: u64,
) -> Result<SchurReport, MantleError> {
    let n = sys.cell_count();
    let steps = steps.min(n.saturating_sub(1)).max(1);
    // |G_imp 1| against |G_imp q| for a random q of the SAME M_p norm: the
    // ratio is the constant mode's leak relative to the gradient's own size,
    // which is the only scale-free way to say "exactly null".
    let ones = vec![1.0f64; n];
    let g_one = sys.gradient_implicit(&ones)?;
    let mut probe = pseudo_random_cells(sys, seed ^ 0x5eed);
    let ratio = sys.norm_p(&ones) / sys.norm_p(&probe);
    for slot in probe.iter_mut() {
        *slot *= ratio;
    }
    let scale = sys.norm_u(&sys.gradient_implicit(&probe)?);
    let constant_mode_leak = if scale > 0.0 {
        sys.norm_u(&g_one) / scale
    } else {
        sys.norm_u(&g_one)
    };

    let mut v = pseudo_random_cells(sys, seed);
    sys.remove_mean_pressure(&mut v)?;
    let nv = sys.norm_p(&v);
    for slot in v.iter_mut() {
        *slot /= nv;
    }
    let mut basis: Vec<Vec<f64>> = vec![v];
    let mut alpha: Vec<f64> = Vec::new();
    let mut beta: Vec<f64> = Vec::new();

    for j in 0..steps {
        let mut w = schur_apply(sys, precon, &basis[j])?;
        sys.remove_mean_pressure(&mut w)?;
        let a = sys.dot_p(&basis[j], &w);
        alpha.push(a);
        // Full reorthogonalisation, two passes: same standing reason as the
        // A-block sweep, a three-term recurrence returns ghost copies of a
        // converged extreme instead of new information.
        for _pass in 0..2 {
            for q in basis.iter() {
                let c = sys.dot_p(q, &w);
                for (slot, x) in w.iter_mut().zip(q.iter()) {
                    *slot -= c * x;
                }
            }
        }
        // The mean is projected out AGAIN after reorthogonalisation. Without
        // it the constant mode leaks back in through roundoff once the Krylov
        // space approaches the full dimension, and the sweep then reports the
        // constant's eigenvalue (roundoff-sized) as lambda_min. Measured at
        // n = 2, where cells - 1 steps is the whole space.
        sys.remove_mean_pressure(&mut w)?;
        let b = sys.norm_p(&w);
        let cut = 1e-13 * alpha.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        if j + 1 == steps || b <= cut {
            beta.push(b);
            break;
        }
        beta.push(b);
        for slot in w.iter_mut() {
            *slot /= b;
        }
        basis.push(w);
    }

    let m = alpha.len();
    let mut t = vec![vec![0.0f64; m]; m];
    for i in 0..m {
        t[i][i] = alpha[i];
        if i + 1 < m {
            t[i][i + 1] = beta[i];
            t[i + 1][i] = beta[i];
        }
    }
    let (values, vectors) = jacobi_symmetric(t);
    let beta_m = *beta.last().unwrap_or(&0.0);
    let lo = 0usize;
    let hi = m - 1;
    let bracket = |i: usize, lam: f64| {
        if lam.abs() > 0.0 {
            beta_m * vectors[m - 1][i].abs() / lam.abs()
        } else {
            f64::INFINITY
        }
    };
    Ok(SchurReport {
        steps: m,
        cells: n,
        lambda_min: sys.eta() * values[lo],
        lambda_max: sys.eta() * values[hi],
        bracket_min: bracket(lo, values[lo]),
        bracket_max: bracket(hi, values[hi]),
        constant_mode_leak,
        inner_cg_iterations: precon.inner_iterations(),
    })
}

/// `q -> -D B^+ G_imp q`, positive semidefinite in `M_p`.
///
/// The sign is not a convention: `<q, D v>_{M_p} = -<G_imp q, v>_{M_u}` is the
/// adjoint identity, so with `v = B^+ G_imp q` the quadratic form is
/// `-<q, D B^+ G q>_{M_p} = <G q, B^+ G q>_{M_u} >= 0`.
pub fn schur_apply(
    sys: &SaddleSystem,
    precon: &mut BlockPreconditioner,
    q: &[f64],
) -> Result<Vec<f64>, MantleError> {
    let g = sys.gradient_implicit(q)?;
    let v = precon.solve_a_block(sys, &g)?;
    let d = sys.divergence(&v)?;
    Ok(d.iter().map(|x| -x).collect())
}

/// Deterministic pseudo-random mean-free cell field.
pub fn pseudo_random_cells(sys: &SaddleSystem, seed: u64) -> Vec<f64> {
    let mut state = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    let mut out = vec![0.0f64; sys.cell_count()];
    for slot in out.iter_mut() {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        *slot = ((state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0;
    }
    out
}

// ---------------------------------------------------------------------
// The deleted boundary flux: the exact diagnostic for the rotation leak

/// Which part of the cubed sphere a boundary facet sits on.
///
/// Split out because "the seams are where a cubed sphere goes wrong" is an
/// intuition this project has already measured to be BACKWARDS for facet
/// non-orthogonality (build brief, correction of 2026-09-08). Reporting the
/// deleted flux by class rather than as one average is how that gets checked
/// again rather than assumed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FacetClass {
    /// Neither `i` nor `j` on a panel edge.
    PanelInterior,
    /// Exactly one of `i`, `j` on a panel edge.
    Seam,
    /// Both on a panel edge.
    Corner,
}

impl FacetClass {
    pub fn label(self) -> &'static str {
        match self {
            FacetClass::PanelInterior => "panel-interior",
            FacetClass::Seam => "seam",
            FacetClass::Corner => "corner",
        }
    }
}

/// `m_F = integral_F (x x n) dS` on the boundary facets, and the divergence
/// it forces a rigid rotation to carry.
///
/// The identity, from the module docs and exact for any closed polyhedron:
///
/// ```text
/// sum over ALL facets of c of s(c,f) m_f = 0
/// (D r_omega)_c = -(1 / V_c) sum over BOUNDARY facets of s(c,f) omega . m_f
/// ```
///
/// so `m_F` is not a proxy for the rotation leak, it IS the rotation leak.
/// [`Self::identity_error`] is the measured gap between the two sides.
#[derive(Debug, Clone)]
pub struct BoundaryMomentReport {
    pub inner_faces: usize,
    pub outer_faces: usize,
    /// `(class, boundary kind, count, max |m_F|, rms |m_F|)`, with `|m_F|`
    /// divided by `|A_F| * r` so it is dimensionless and comparable across
    /// rungs: for an exactly spherical facet it would be zero, and for a flat
    /// facet it is the sine of the angle between the facet centroid and its
    /// own normal.
    pub by_class: Vec<(FacetClass, bool, usize, f64, f64)>,
    /// `max |m_F| / (|A_F| r)` over all boundary facets.
    pub max_normalised: f64,
    pub rms_normalised: f64,
    /// Induced `|D r|_{M_p}` for rotations about x, y, z, with each rotation
    /// mode normalised to unit `M_u` norm.
    pub induced_divergence: [f64; 3],
    /// The same, times `h_ref`: dimensionless, and the number to fit against
    /// `h`.
    pub induced_dimensionless: [f64; 3],
    /// `max_c |(D r)_c - predicted_c| / max_c |(D r)_c|` over the three
    /// rotations, where `predicted` is the boundary-moment sum above. This is
    /// the claim "m_F accounts for the whole leak", measured.
    pub identity_error: f64,
    /// Fraction of `|D r|^2` that sits in cells with NO boundary facet. The
    /// identity says exactly zero.
    pub interior_share: f64,
}

/// Measure `m_F` on the boundary and check it predicts `D r` exactly.
#[allow(clippy::needless_range_loop)]
pub fn boundary_moment_report(sys: &SaddleSystem) -> Result<BoundaryMomentReport, MantleError> {
    let topo = sys.topology();
    let spec = topo.spec();
    let n = spec.n();
    let r_top = spec.r_top_m();
    let mut by: Vec<(FacetClass, bool, usize, f64, f64)> = Vec::new();
    for class in [
        FacetClass::PanelInterior,
        FacetClass::Seam,
        FacetClass::Corner,
    ] {
        for outer in [false, true] {
            by.push((class, outer, 0, 0.0, 0.0));
        }
    }
    let mut inner_faces = 0usize;
    let mut outer_faces = 0usize;
    let mut max_norm = 0.0f64;
    let mut sum_sq = 0.0f64;

    // predicted[c] = -(1/V_c) sum_{f boundary} s(c,f) m_f, one 3-vector per
    // cell; dotted with omega it is the whole of (D r_omega)_c.
    let mut predicted = vec![[0.0f64; 3]; sys.cell_count()];
    for c in 0..sys.cell_count() {
        for (f, s) in topo.cell_faces(c as u32) {
            let face = &topo.faces()[f as usize];
            if face.neighbor.is_some() {
                continue;
            }
            let inv = -(s as f64) / sys.m_p()[c];
            for (slot, m) in predicted[c].iter_mut().zip(face.moment_m3.iter()) {
                *slot += inv * m;
            }
        }
    }

    for face in topo.faces().iter() {
        if face.neighbor.is_some() {
            continue;
        }
        let (_, i, j, _) = spec
            .fijk(face.owner as u64)
            .expect("owner index is within the cell count");
        let edge_i = i == 0 || i + 1 == n;
        let edge_j = j == 0 || j + 1 == n;
        let class = match (edge_i, edge_j) {
            (true, true) => FacetClass::Corner,
            (false, false) => FacetClass::PanelInterior,
            _ => FacetClass::Seam,
        };
        let outer = matches!(face.kind, pf_grid::shell_ops::FaceKind::BoundaryOuter);
        if outer {
            outer_faces += 1;
        } else {
            inner_faces += 1;
        }
        let m = (face.moment_m3.iter().map(|x| x * x).sum::<f64>()).sqrt();
        let normalised = m / (face.area_m2 * r_top);
        max_norm = max_norm.max(normalised);
        sum_sq += normalised * normalised;
        let slot = by
            .iter_mut()
            .find(|e| e.0 == class && e.1 == outer)
            .expect("every (class, side) pair was seeded above");
        slot.2 += 1;
        slot.3 = slot.3.max(normalised);
        slot.4 += normalised * normalised;
    }
    let total = (inner_faces + outer_faces) as f64;
    for slot in by.iter_mut() {
        if slot.2 > 0 {
            slot.4 = (slot.4 / slot.2 as f64).sqrt();
        }
    }

    let rotations = ModeBasis::rotations(topo, sys)?;
    let mut induced = [0.0f64; 3];
    let mut induced_nd = [0.0f64; 3];
    let mut identity_error = 0.0f64;
    let mut interior_share = 0.0f64;
    let mut boundary_cell = vec![false; sys.cell_count()];
    for face in topo.faces().iter() {
        if face.neighbor.is_none() {
            boundary_cell[face.owner as usize] = true;
        }
    }
    for axis in 0..3 {
        let mode = rotations.mode(axis);
        let d = sys.divergence(mode)?;
        induced[axis] = sys.norm_p(&d);
        induced_nd[axis] = sys.mesh_scale() * induced[axis];
        // The rotation modes were M_u-orthonormalised, so recover the omega
        // each one carries by projecting onto the raw rotation fields. The
        // orthonormal mode is a linear combination of the three raw
        // rotations, so build the prediction from the SAME combination.
        let omega = rotation_weights(topo, sys, mode)?;
        let mut worst = 0.0f64;
        let mut scale = 0.0f64;
        let mut inside = 0.0f64;
        let mut total_sq = 0.0f64;
        for (c, value) in d.iter().enumerate() {
            let pred = predicted[c]
                .iter()
                .zip(omega.iter())
                .map(|(a, b)| a * b)
                .sum::<f64>();
            worst = worst.max((value - pred).abs());
            scale = scale.max(value.abs());
            total_sq += sys.m_p()[c] * value * value;
            if !boundary_cell[c] {
                inside += sys.m_p()[c] * value * value;
            }
        }
        if scale > 0.0 {
            identity_error = identity_error.max(worst / scale);
        }
        if total_sq > 0.0 {
            interior_share = interior_share.max(inside / total_sq);
        }
    }

    Ok(BoundaryMomentReport {
        inner_faces,
        outer_faces,
        by_class: by,
        max_normalised: max_norm,
        rms_normalised: (sum_sq / total).sqrt(),
        induced_divergence: induced,
        induced_dimensionless: induced_nd,
        identity_error,
        interior_share,
    })
}

/// Recover the angular velocity a face field represents, by least squares
/// against the three unit rotation fields.
///
/// Used only to build the boundary-moment PREDICTION for an orthonormalised
/// rotation mode. It is a 3x3 solve on the Gram matrix of the raw rotation
/// fields in `M_u`, and it is exact for any field that IS a rigid rotation,
/// which every mode here is by construction.
#[allow(clippy::needless_range_loop)]
fn rotation_weights(
    topo: &ShellTopology,
    sys: &SaddleSystem,
    field: &[f64],
) -> Result<[f64; 3], MantleError> {
    let raw: Vec<Vec<f64>> = (0..3)
        .map(|a| {
            let mut omega = [0.0f64; 3];
            omega[a] = 1.0;
            topo.rigid_rotation_field(omega)
        })
        .collect();
    let mut gram = [[0.0f64; 3]; 3];
    let mut rhs = [0.0f64; 3];
    for a in 0..3 {
        rhs[a] = sys.dot_u(&raw[a], field);
        for b in 0..3 {
            gram[a][b] = sys.dot_u(&raw[a], &raw[b]);
        }
    }
    let g: Vec<Vec<f64>> = gram.iter().map(|r| r.to_vec()).collect();
    let out = solve_symmetric(&g, &rhs)?;
    Ok([out[0], out[1], out[2]])
}

// ---------------------------------------------------------------------
// Small dense linear algebra

/// Cyclic Jacobi eigendecomposition of a small symmetric matrix.
///
/// Returns `(eigenvalues ascending, eigenvectors)` with `vectors[r][i]` the
/// `r`-th component of the `i`-th eigenvector - so `vectors[m-1][i]` is the
/// last component the Ritz residual bound needs.
// Index loops throughout: a Jacobi rotation reads and writes `a[k][p]`,
// `a[k][q]`, `a[p][k]` and `a[q][k]` in the same pass, which iterator form
// cannot express without splitting borrows and which would obscure the one
// thing a reader needs to check here - that the rotation is applied on both
// sides and accumulated into `v`.
#[allow(clippy::needless_range_loop)]
fn jacobi_symmetric(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let m = a.len();
    let mut v = vec![vec![0.0f64; m]; m];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _sweep in 0..60 {
        let mut off = 0.0f64;
        for i in 0..m {
            for j in (i + 1)..m {
                off += a[i][j] * a[i][j];
            }
        }
        if off.sqrt() <= 1e-300 {
            break;
        }
        let mut moved = false;
        for p in 0..m {
            for q in (p + 1)..m {
                if a[p][q].abs() <= 1e-18 * (a[p][p].abs() + a[q][q].abs()).max(1e-300) {
                    continue;
                }
                moved = true;
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..m {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..m {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let vp = row[p];
                    let vq = row[q];
                    row[p] = c * vp - s * vq;
                    row[q] = s * vp + c * vq;
                }
            }
        }
        if !moved {
            break;
        }
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|x, y| {
        a[*x][*x]
            .partial_cmp(&a[*y][*y])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let values: Vec<f64> = order.iter().map(|i| a[*i][*i]).collect();
    let vectors: Vec<Vec<f64>> = (0..m)
        .map(|r| order.iter().map(|i| v[r][*i]).collect())
        .collect();
    (values, vectors)
}

/// Solve a symmetric positive semidefinite system through its
/// eigendecomposition, refusing loudly when it is singular RELATIVE TO ITSELF.
///
/// The relative test is the point: the rotation coarse space's Gram matrix has
/// entries of order `h^4` and is perfectly well conditioned internally. An
/// absolute pivot test would reject it at fine resolution and quietly hand
/// back a zero coarse correction, which is the exact failure this module is
/// supposed to be measuring.
#[allow(clippy::needless_range_loop)]
fn solve_symmetric(m: &[Vec<f64>], rhs: &[f64]) -> Result<Vec<f64>, MantleError> {
    let n = rhs.len();
    let a: Vec<Vec<f64>> = m.iter().map(|r| r.to_vec()).collect();
    let (values, vectors) = jacobi_symmetric(a);
    let scale = values.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
    if scale <= 0.0 {
        return Err(MantleError::field(
            "symmetric coarse system scale",
            scale,
            "the matrix is entirely zero, so the fit has no content",
        ));
    }
    let mut out = vec![0.0f64; n];
    for i in 0..n {
        if values[i].abs() <= 1e-12 * scale {
            return Err(MantleError::field(
                "symmetric coarse system relative eigenvalue",
                values[i] / scale,
                "the fit directions are numerically dependent",
            ));
        }
        let dot: f64 = (0..n).map(|r| vectors[r][i] * rhs[r]).sum();
        let coeff = dot / values[i];
        for (r, slot) in out.iter_mut().enumerate() {
            *slot += coeff * vectors[r][i];
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------
// Forcings and manufactured pairs

/// A buoyancy-SHAPED forcing: a radial body force with a degree-two lateral
/// pattern, projected onto the face normals.
///
/// **This is a forcing with the right shape, not a physical buoyancy field.**
/// It exists so the iteration tables below are driven by something smoother
/// and more representative than white noise; no physical claim is attached to
/// it, and MC-Q1 supplies the real one.
pub fn buoyancy_shaped_force(topo: &ShellTopology) -> Vec<f64> {
    let r_top = topo.spec().r_top_m();
    super::stokes_solve::sampled_face_field(topo, move |x| {
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt().max(1e-30);
        let c = x[2] / r;
        // P2(cos theta) laterally, one half-wavelength radially: a shape a
        // mantle solve will actually meet, and not an eigenfunction of
        // anything here.
        let depth = (std::f64::consts::PI * (r - topo.spec().r_bottom_m())
            / (r_top - topo.spec().r_bottom_m()))
        .sin();
        let amp = depth * (1.5 * c * c - 0.5);
        [amp * x[0] / r, amp * x[1] / r, amp * x[2] / r]
    })
}

/// How big the deferred correction actually is, against the implicit term it
/// corrects.
///
/// `|G_exp p|_{M_u} / |G_imp p|_{M_u}`, measured on a smooth pressure and on
/// a random one. The build brief records this reaching `1.00` at `n = 32`,
/// which is the edge of where a fixed point on the correction converges at
/// all; this function is how that inherited number is re-measured here rather
/// than carried forward on trust.
pub fn correction_parity(sys: &SaddleSystem, p: &[f64]) -> Result<f64, MantleError> {
    let imp = sys.gradient_implicit(p)?;
    let exp = sys.gradient_explicit(p)?;
    let d = sys.norm_u(&imp);
    if d <= 0.0 {
        return Ok(f64::NAN);
    }
    Ok(sys.norm_u(&exp) / d)
}

/// A smooth mean-free pressure field, for parity and manufactured checks.
pub fn smooth_cell_field(sys: &SaddleSystem) -> Result<Vec<f64>, MantleError> {
    let topo = sys.topology();
    let centroids = topo
        .cell_centroids_m()
        .map_err(|e| shell("cell centroids", e))?;
    let k = std::f64::consts::PI / topo.spec().r_top_m();
    let mut out: Vec<f64> = centroids
        .iter()
        .map(|x| (k * x[0]).sin() * (k * x[1]).cos() + 0.5 * (k * x[2]).sin())
        .collect();
    sys.remove_mean_pressure(&mut out)?;
    Ok(out)
}

/// A manufactured DISCRETE pair `(u, p)` and the right-hand side the
/// operators of this module give it.
///
/// Everything is built in the discrete space and the continuum never enters,
/// which is what makes this a self-consistency instrument and not a physics
/// one.
#[derive(Debug, Clone)]
pub struct ManufacturedPair {
    pub u: Vec<f64>,
    pub p: Vec<f64>,
    pub rhs_u: Vec<f64>,
    pub rhs_p: Vec<f64>,
}

/// Build the manufactured pair. `u` is a smooth sampled face field plus a
/// little noise, `p` is smooth and mean free, and the right-hand side is
/// `(B u + G p, D u)` - so the pair is the exact solution of the discrete
/// system BY CONSTRUCTION and any error is the solver's.
/// `restrict` says which velocity modes the manufactured SOLUTION is allowed
/// to contain; `equations` says which momentum rows the right-hand side is
/// projected onto. **They are different questions and conflating them hides
/// the only interesting failure here.** A solution built in the space the
/// preconditioner can represent tests that the solver inverts the operator; a
/// solution that also carries the six `A`-null components tests whether the
/// solver can recover them, which for the three rotations it structurally
/// cannot.
pub fn manufactured_pair(
    sys: &SaddleSystem,
    restrict: &ModeBasis,
    equations: &ModeBasis,
    mode: GradientMode,
    seed: u64,
) -> Result<ManufacturedPair, MantleError> {
    let topo = sys.topology();
    let k = std::f64::consts::PI / topo.spec().r_top_m();
    let mut u = super::stokes_solve::sampled_face_field(topo, |x| {
        [
            (k * x[2]).sin() + 0.4 * (k * x[1]).cos(),
            (k * x[0]).sin() + 0.4 * (k * x[2]).cos(),
            (k * x[1]).sin() + 0.4 * (k * x[0]).cos(),
        ]
    });
    let noise = super::stokes_solve::pseudo_random_field(sys.a_block(), seed);
    let nu = sys.norm_u(&u);
    let nn = sys.norm_u(&noise);
    for (slot, x) in u.iter_mut().zip(noise.iter()) {
        *slot = *slot / nu + 0.1 * x / nn;
    }
    // If the caller's solve deflates velocity modes, the manufactured
    // solution must live in the same space, or the comparison would charge
    // the solver for a component the problem statement removed.
    restrict.project(sys, &mut u);

    // BALANCE THE PAIR, and this is load bearing rather than cosmetic.
    //
    // `u` and `p` are manufactured independently, so their natural scales are
    // unrelated: with SI geometry, eta = 1 and both fields O(1) in their own
    // units, `|G p|_{M_u}` came out FIVE decades above `|B u|_{M_u}`. The
    // right-hand side is then essentially the pressure gradient alone, the
    // residual is blind to the velocity, and a solve can return a velocity
    // 3e4 times too large at a relative residual of 7e-13 - measured, before
    // this rescale was added. A real Stokes solution has the two terms in
    // balance because that is what the momentum equation says, so the
    // manufactured pair is put in balance too.
    let mut p = smooth_cell_field(sys)?;
    let bu = sys.norm_u(&sys.viscous(&u)?);
    let gp = sys.norm_u(&sys.gradient(&p, mode)?);
    if !gp.is_finite() || gp <= 0.0 || !bu.is_finite() || bu <= 0.0 {
        return Err(MantleError::field(
            "manufactured pair balance",
            gp,
            "one half of the manufactured momentum balance is empty",
        ));
    }
    let factor = bu / gp;
    for slot in p.iter_mut() {
        *slot *= factor;
    }
    let (rhs_u, rhs_p) = sys.apply(equations, mode, &u, &p)?;
    Ok(ManufacturedPair { u, p, rhs_u, rhs_p })
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_grid::shell::ShellSpec;

    const R_CMB: f64 = 3.480e6;
    const R_TOP: f64 = 6.271e6;
    /// `eta = 1`, not the mantle's `1e21`. Every ratio reported here is
    /// invariant under it, and
    /// [`the_solve_does_not_depend_on_the_viscosity_scale`] measures that
    /// rather than asserting it.
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

    /// The refinement ladder.
    ///
    /// The default two rungs run inside the normal debug suite.
    /// `MC0_SADDLE_LADDER=full` adds `(8, 6)` and `(16, 12)`;
    /// `MC0_SADDLE_LADDER=deep` adds `(32, 24)` on top. Following the
    /// `MC0_ABLOCK_LADDER` precedent exactly: the coverage is identical at
    /// every rung and only the cost is opt-in. `(4,3), (8,6), (16,12),
    /// (32,24)` is an exactly self-similar factor-two ladder; `(2,2)` is not
    /// (`layers = 1.5` does not exist) so it is reported and excluded from
    /// every fitted power.
    fn ladder() -> Vec<(usize, usize)> {
        let mut v = vec![(2usize, 2usize), (4, 3)];
        match std::env::var("MC0_SADDLE_LADDER").as_deref() {
            Ok("full") => {
                v.push((8, 6));
                v.push((16, 12));
            }
            Ok("deep") => {
                v.push((8, 6));
                v.push((16, 12));
                v.push((32, 24));
            }
            _ => {}
        }
        v
    }

    // -----------------------------------------------------------------
    // 1. The pieces, before anything is solved with them

    #[test]
    fn the_saddle_pieces_are_the_adjoint_pair_the_algebra_claims() {
        // Everything downstream rests on three exact identities. If any of
        // them is only approximately true, the block preconditioner is
        // preconditioning a different operator and the iteration tables mean
        // nothing.
        for (n, layers) in [(2usize, 2usize), (4, 3)] {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let u = super::super::stokes_solve::pseudo_random_field(sys.a_block(), 11);
            let p = pseudo_random_cells(&sys, 22);

            // (a) G_imp is the exact negative M_u-weighted adjoint of D.
            let g = sys.gradient_implicit(&p).expect("implicit gradient");
            let d = sys.divergence(&u).expect("divergence");
            let lhs = sys.dot_u(&u, &g);
            let rhs = -sys.dot_p(&d, &p);
            let adjoint = (lhs - rhs).abs() / lhs.abs().max(rhs.abs());
            assert!(
                adjoint < 1e-12,
                "n={n}: <u, G p>_Mu = {lhs:.6e} and -<D u, p>_Mp = {rhs:.6e} \
                 differ by {adjoint:.3e} relative; the implicit gradient is not \
                 -M_u^-1 D^T M_p and the whole block structure is void"
            );

            // (b) D u is EXACTLY M_p-mean-free, by pairwise facet
            // cancellation. Nothing enforces this; it either holds or the
            // continuity equations are inconsistent.
            let total: f64 = d
                .iter()
                .zip(sys.m_p().iter())
                .map(|(x, v)| x * v)
                .sum::<f64>();
            let d_scale = sys.norm_p(&d) * (sys.m_p().iter().sum::<f64>()).sqrt();
            let drift = total.abs() / d_scale;
            assert!(
                drift < 1e-12,
                "n={n}: sum_c V_c (D u)_c = {total:.6e}, i.e. {drift:.3e} of its own \
                 scale; the continuity row is not automatically mean free"
            );

            // (c) The constant pressure is an exact null mode of G, which is
            // why it is the ONE thing this module projects out of the
            // equations.
            let ones = vec![1.0f64; sys.cell_count()];
            let g_one = sys
                .gradient_implicit(&ones)
                .expect("gradient of a constant");
            let leak = sys.norm_u(&g_one) / sys.norm_u(&g) * sys.norm_p(&p) / sys.norm_p(&ones);
            assert!(
                leak < 1e-14,
                "n={n}: G_imp of a constant pressure has relative size {leak:.3e}, \
                 so the pressure constant is not an exact null mode"
            );
            println!(
                "n={n:<3} adjoint_err={adjoint:.3e} mean_free_drift={drift:.3e} \
                 constant_mode_leak={leak:.3e}"
            );
        }
    }

    #[test]
    fn the_mode_basis_agrees_with_the_a_block_null_space() {
        // ModeBasis re-orthonormalises the raw modes rather than slicing
        // NullSpace's output, so the two could in principle span different
        // subspaces. For the FULL six they must not, and the projections must
        // agree to roundoff.
        for (n, layers) in [(2usize, 2usize), (4, 3)] {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let six = ModeBasis::build(&t, &sys, Deflation::Six).expect("six modes");
            let null = NullSpace::build(&t, sys.a_block()).expect("A-block null space");
            let base = super::super::stokes_solve::pseudo_random_field(sys.a_block(), 33);
            let mut mine = base.clone();
            six.project(&sys, &mut mine);
            let mut theirs = base.clone();
            null.project(sys.a_block(), &mut theirs);
            let diff: f64 = {
                let d: Vec<f64> = mine.iter().zip(theirs.iter()).map(|(x, y)| x - y).collect();
                sys.norm_u(&d) / sys.norm_u(&base)
            };
            assert!(
                diff < 1e-12,
                "n={n}: the six-mode projector here and NullSpace's differ by {diff:.3e}; \
                 they do not span the same subspace"
            );
            println!(
                "n={n:<3} projector_gap={diff:.3e} pivot={:.6} null_pivot={:.6}",
                six.smallest_pivot(),
                null.smallest_pivot()
            );
        }
    }

    #[test]
    fn the_six_a_block_modes_are_exact_for_k_and_near_null_for_the_saddle_point() {
        // The table this module's doc comment asserts, measured. It is the
        // evidence for treating the six as exact INSIDE the A-block
        // preconditioner and as not-null in the coupled equations.
        println!("  n  lay        h        mode   |K n|/|K|   |D n| (abs)   h|D n|   bnd share");
        let mut hs: Vec<f64> = Vec::new();
        let mut rot: Vec<f64> = Vec::new();
        let mut tra: Vec<f64> = Vec::new();
        for (n, layers) in ladder() {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let translations =
                ModeBasis::build(&t, &sys, Deflation::Translations).expect("translations");
            let rotations = ModeBasis::build(&t, &sys, Deflation::Rotations).expect("rotations");
            let probe = super::super::stokes_solve::pseudo_random_field(sys.a_block(), 44);
            let k_scale = sys.norm_u(&sys.viscous(&probe).expect("B probe")) / sys.norm_u(&probe);

            let lt = divergence_leak(&sys, &translations).expect("translation leak");
            let lr = divergence_leak(&sys, &rotations).expect("rotation leak");
            for (label, basis, leak) in [
                ("translation", &translations, &lt),
                ("rotation", &rotations, &lr),
            ] {
                for j in 0..3 {
                    let kn = sys.norm_u(&sys.viscous(basis.mode(j)).expect("B mode")) / k_scale;
                    println!(
                        "{n:>3} {layers:>4} {:.4e} {label:>11}{j}  {kn:.3e}   {:.6e}  {:.4e}  {:.4}",
                        sys.mesh_scale(),
                        leak.absolute[j],
                        leak.dimensionless[j],
                        leak.boundary_share[j]
                    );
                    assert!(
                        kn < 1e-12,
                        "n={n} {label}{j}: |B n|/|B| = {kn:.3e}; this mode is NOT an exact \
                         null mode of the A block and deflating it inside the preconditioner \
                         would be wrong"
                    );
                }
            }
            let worst_rot = lr.dimensionless.iter().cloned().fold(0.0f64, f64::max);
            let worst_tra = lt.dimensionless.iter().cloned().fold(0.0f64, f64::max);
            // The GAP, not its size: the claim is that translations are O(1)
            // and rotations O(h^2), so the ratio must grow with refinement
            // and the fitted slopes below are where that is checked. Here the
            // only assertion is that the two classes are already separated at
            // the coarsest rung, where the ratio is smallest.
            println!(
                "{n:>3} {layers:>4} {:.4e}   leak ratio translation/rotation = {:.4}",
                sys.mesh_scale(),
                worst_tra / worst_rot
            );
            assert!(
                worst_tra > 10.0 * worst_rot,
                "n={n}: translations leak {worst_tra:.3e} and rotations {worst_rot:.3e}, \
                 a ratio of only {:.2}; the two classes are not separated at all and \
                 the near-null story of this module would be wrong",
                worst_tra / worst_rot
            );
            for j in 0..3 {
                assert!(
                    lr.boundary_share[j] > 1.0 - 1e-10,
                    "n={n} rotation{j}: only {:.6} of |D r|^2 sits in boundary cells, but \
                     the closed-polyhedron identity says all of it does",
                    lr.boundary_share[j]
                );
            }
            hs.push(sys.mesh_scale());
            rot.push(worst_rot);
            tra.push(worst_tra);
        }
        if hs.len() >= 3 {
            // Skip (2,2): it is not on the self-similar series.
            let s_rot = slope(&hs[1..], &rot[1..]);
            let s_tra = slope(&hs[1..], &tra[1..]);
            println!("fitted: h|D r| ~ h^{s_rot:+.3}   h|D t| ~ h^{s_tra:+.3}");
            assert!(
                s_rot > 1.5,
                "the rotation divergence falls at h^{s_rot:+.3}, not the O(h^2) the \
                 boundary-facet argument predicts"
            );
        }
    }

    // -----------------------------------------------------------------
    // 2. The deleted boundary flux (the exact diagnostic)

    #[test]
    fn the_deleted_boundary_flux_predicts_the_rotation_divergence() {
        // The claim under test, and it is an IDENTITY rather than a model:
        //
        //   q_F(omega) = int_F (omega x x) . n dS = omega . m_F,  m_F = int_F (x x n) dS
        //   sum over ALL facets of a closed cell of s(c,f) m_f = 0
        //   => (D r_omega)_c = -(1/V_c) sum over BOUNDARY facets of s(c,f) omega . m_f
        //
        // If m_F accounts for the whole of D r, the rotation leak is a
        // boundary-geometry approximation with an exact diagnostic. If it
        // does not, something else is wrong and that is the bigger finding.
        println!(
            "  n  lay        h    innerF  outerF   max|m|/(|A|r)  rms|m|/(|A|r)   \
             h|D r|      identity_err  interior_share"
        );
        let mut hs: Vec<f64> = Vec::new();
        let mut mx: Vec<f64> = Vec::new();
        let mut dr: Vec<f64> = Vec::new();
        let mut rows: Vec<(usize, BoundaryMomentReport)> = Vec::new();
        for (n, layers) in ladder() {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let r = boundary_moment_report(&sys).expect("boundary moments");
            println!(
                "{n:>3} {layers:>4} {:.4e} {:>7} {:>7}   {:.6e}   {:.6e}   {:.4e}  {:.3e}   {:.3e}",
                sys.mesh_scale(),
                r.inner_faces,
                r.outer_faces,
                r.max_normalised,
                r.rms_normalised,
                r.induced_dimensionless[0],
                r.identity_error,
                r.interior_share
            );
            // The identity is exact arithmetic on the same stored moments, so
            // the only admissible discrepancy is roundoff.
            assert!(
                r.identity_error < 1e-10,
                "n={n}: the boundary-moment prediction misses D r by {:.3e} relative. \
                 The identity is exact, so this means the leak has a second source and \
                 the boundary geometry does NOT account for it",
                r.identity_error
            );
            assert!(
                r.interior_share < 1e-20,
                "n={n}: {:.3e} of |D r|^2 sits in cells with no boundary facet, but the \
                 closed-polyhedron identity forbids any",
                r.interior_share
            );
            hs.push(sys.mesh_scale());
            mx.push(r.max_normalised);
            dr.push(r.induced_dimensionless[0]);
            rows.push((n, r));
        }
        println!();
        println!("  n   class            side   count    max|m|/(|A|r)   rms|m|/(|A|r)");
        for (n, r) in rows.iter() {
            for (class, outer, count, max, rms) in r.by_class.iter() {
                if *count == 0 {
                    continue;
                }
                println!(
                    "{n:>3}   {:<15} {:<6} {count:>6}   {max:.6e}    {rms:.6e}",
                    class.label(),
                    if *outer { "outer" } else { "inner" }
                );
            }
        }
        if hs.len() >= 3 {
            let s_m = slope(&hs[1..], &mx[1..]);
            let s_d = slope(&hs[1..], &dr[1..]);
            println!();
            println!("fitted (self-similar rungs only): max|m|/(|A|r) ~ h^{s_m:+.3}   h|D r| ~ h^{s_d:+.3}");
        }
    }

    // -----------------------------------------------------------------
    // 3. The preconditioner

    #[test]
    fn the_inner_solve_is_metric_consistent() {
        // The A-block solve converts its right-hand side into the inner
        // metric. That conversion is the one place a metric mistake could
        // hide, and the identity that makes it safe is
        //   <n_j, P^-1 M_u a>_P = n_j^T M_u a = <n_j, a>_{M_u}
        // i.e. rhs consistency is metric independent. Measured by solving the
        // SAME problem in both metrics and comparing the answers.
        let (t, a) = op(4, 3);
        let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
            .expect("system builds");
        let six = ModeBasis::build(&t, &sys, Deflation::Six).expect("six modes");
        let mut rhs = super::super::stokes_solve::pseudo_random_field(sys.a_block(), 55);
        six.project(&sys, &mut rhs);
        let mut answers: Vec<Vec<f64>> = Vec::new();
        for metric in [Metric::FaceMass, Metric::StiffnessDiagonal] {
            let mut pre = BlockPreconditioner::build(
                &t,
                &a,
                &sys,
                PreconditionerOptions {
                    inner_tol: 1e-12,
                    inner_metric: metric,
                    ..PreconditionerOptions::default()
                },
            )
            .expect("preconditioner builds");
            let mut x = pre.solve_a_block(&sys, &rhs).expect("A-block solve");
            six.project(&sys, &mut x);
            // Residual of the ORIGINAL equation B x = rhs, in M_u.
            let bx = sys.viscous(&x).expect("B x");
            let r: Vec<f64> = rhs.iter().zip(bx.iter()).map(|(p, q)| p - q).collect();
            let rel = sys.norm_u(&r) / sys.norm_u(&rhs);
            println!(
                "metric={:<9} inner_iters={:<5} residual of B x = rhs in M_u: {rel:.3e}",
                metric.label(),
                pre.inner_iterations()
            );
            assert!(
                rel < 1e-9,
                "metric={}: solving K x = M_u a in this metric left {rel:.3e} of the \
                 M_u residual; the right-hand side conversion is wrong",
                metric.label()
            );
            answers.push(x);
        }
        let d: Vec<f64> = answers[0]
            .iter()
            .zip(answers[1].iter())
            .map(|(x, y)| x - y)
            .collect();
        let gap = sys.norm_u(&d) / sys.norm_u(&answers[0]);
        println!("the two metrics return the same vector to {gap:.3e}");
        assert!(
            gap < 1e-7,
            "the two inner metrics returned different solutions ({gap:.3e} apart); \
             they are solving different systems"
        );
    }

    #[test]
    fn the_exact_block_preconditioner_terminates_in_two_iterations() {
        // The structural claim of BlockPreconditioner's doc comment: with
        // exact blocks, (M P^-1)^2 = I, so GMRES terminates in two steps.
        // This is the check that the sign, the transpose and the metric are
        // wired the way the comment says. It is affordable only because the
        // Schur complement is assembled densely, one A-solve per cell, which
        // is why it runs at n = 2 only.
        let (t, a) = op(2, 2);
        let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
            .expect("system builds");
        let six = ModeBasis::build(&t, &sys, Deflation::Six).expect("six modes");
        let none = ModeBasis::build(&t, &sys, Deflation::None).expect("empty basis");
        let mut pre = BlockPreconditioner::build(
            &t,
            &a,
            &sys,
            PreconditionerOptions {
                // "Exact" here means as exact as CG gets; 1e-14 is at the
                // roundoff floor for this kappa and is not a tuned number.
                inner_tol: 1e-14,
                // M_u and not the cheaper Jacobi metric: only with the
                // deflation metric equal to M_u is the A-block pseudo-inverse
                // self-adjoint in M_u, and only then is the assembled Schur
                // complement self-adjoint in M_p. See DenseSchur's doc.
                inner_metric: Metric::FaceMass,
                ..PreconditionerOptions::default()
            },
        )
        .expect("preconditioner builds");
        let dense = DenseSchur::build(&sys, &mut pre, 1e-10).expect("dense Schur");
        println!(
            "dense Schur: {} cells, {} mode(s) dropped, asymmetry {:.3e}",
            sys.cell_count(),
            dense.dropped,
            dense.asymmetry
        );
        println!(
            "eta*lambda(-S_0): min(kept)={:.6e} max={:.6e}",
            dense.eigenvalues[dense.dropped],
            dense.eigenvalues[dense.eigenvalues.len() - 1]
        );
        assert_eq!(
            dense.dropped, 1,
            "exactly one pressure mode (the constant) should be null; {} were dropped, \
             so there is a second near-null pressure mode and this discretization has a \
             spurious mode",
            dense.dropped
        );
        assert!(
            dense.asymmetry < 1e-8,
            "the assembled Schur complement is asymmetric by {:.3e} relative; it is \
             supposed to be self-adjoint in M_p",
            dense.asymmetry
        );
        pre.set_dense_schur(dense);
        pre.reset_counters();
        let f = buoyancy_shaped_force(&t);
        let (_, _, rep) = fgmres_saddle(
            &sys,
            &mut pre,
            &none,
            &six,
            GradientMode::Implicit,
            &f,
            &sys.zeros_p(),
            FgmresOptions {
                tol: 1e-10,
                max_iterations: 40,
                ..FgmresOptions::default()
            },
        )
        .expect("FGMRES runs");
        println!(
            "exact blocks: outer={} true_rel={:.3e} history={:?}",
            rep.outer_iterations,
            rep.true_rel_residual,
            rep.history
                .iter()
                .map(|v| format!("{v:.2e}"))
                .collect::<Vec<_>>()
        );
        assert!(
            rep.outer_iterations <= 3,
            "with exact blocks the preconditioned operator satisfies T^2 = I and GMRES \
             must terminate in two steps; it took {}. The block algebra is not what the \
             doc comment claims",
            rep.outer_iterations
        );
    }

    #[test]
    fn the_schur_complement_is_a_scaled_identity_and_its_low_end_is_identified() {
        // Two questions in one sweep.
        //
        // (a) Is `S_0 = D B^+ G` really `~ (c/eta) I` in this operator
        //     scaling, so that the viscosity-weighted pressure-mass Schur
        //     approximation degenerates to a scalar? The pressure mass matrix
        //     cancels ALGEBRAICALLY here, so what is left to measure is the
        //     spread of the spectrum, which is what a scalar cannot follow.
        //
        // (b) What is at the bottom of it? `lambda_min` is the SQUARE of the
        //     discrete inf-sup constant, so a collapse under refinement would
        //     be a spurious pressure mode. But `B^+` is the DEFLATED inverse,
        //     so it also deletes whatever part of `G q` lies in null(K) - and
        //     the truncated translations carry an O(1) share of exactly that.
        //     The eigenvectors are therefore projected onto span{D t_j} and
        //     span{D r_j} to say which explanation is the right one.
        for (n, layers) in [(2usize, 2usize), (4, 3)] {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let mut pre = BlockPreconditioner::build(
                &t,
                &a,
                &sys,
                PreconditionerOptions {
                    inner_tol: 1e-12,
                    // See DenseSchur's doc: self-adjointness of the assembled
                    // Schur needs the M_u deflation metric.
                    inner_metric: Metric::FaceMass,
                    ..PreconditionerOptions::default()
                },
            )
            .expect("preconditioner builds");
            let dense = DenseSchur::build(&sys, &mut pre, 1e-10).expect("dense Schur");
            let m = dense.eigenvalues.len();
            println!(
                "n={n} cells={m} dropped={} asymmetry={:.3e} inner_cg={}",
                dense.dropped,
                dense.asymmetry,
                pre.inner_iterations()
            );
            let low: Vec<String> = dense.eigenvalues[..8.min(m)]
                .iter()
                .map(|v| format!("{v:.4e}"))
                .collect();
            let high: Vec<String> = dense.eigenvalues[m - 4..]
                .iter()
                .map(|v| format!("{v:.4e}"))
                .collect();
            println!("  eta*lambda(-S_0) lowest : {}", low.join("  "));
            println!("  eta*lambda(-S_0) highest: {}", high.join("  "));

            // Identify the low modes. The dense Schur is stored in the
            // M_p^{1/2}-symmetrised basis, so a pressure vector q maps to
            // sqrt(V) q before the overlap is taken.
            let translations =
                ModeBasis::build(&t, &sys, Deflation::Translations).expect("translations");
            let rotations = ModeBasis::build(&t, &sys, Deflation::Rotations).expect("rotations");
            let mut probes: Vec<(&str, Vec<Vec<f64>>)> = Vec::new();
            for (label, basis) in [("D t", &translations), ("D r", &rotations)] {
                let mut set = Vec::new();
                for j in 0..3 {
                    let mut d = sys.divergence(basis.mode(j)).expect("divergence of a mode");
                    sys.remove_mean_pressure(&mut d).expect("mean free");
                    let nrm = sys.norm_p(&d);
                    for slot in d.iter_mut() {
                        *slot /= nrm;
                    }
                    set.push(d);
                }
                probes.push((label, set));
            }
            for i in dense.dropped..(dense.dropped + 5).min(m) {
                let vec_i: Vec<f64> = (0..m).map(|r| dense.eigenvector(r, i)).collect();
                let mut line = format!("  mode {i}: eta*lambda={:.4e}", dense.eigenvalues[i]);
                for (label, set) in probes.iter() {
                    let mut share = 0.0f64;
                    for q in set.iter() {
                        // <q, v> with v in the symmetrised basis: sqrt(V) q.
                        let dot: f64 = (0..m).map(|c| sys.m_p()[c].sqrt() * q[c] * vec_i[c]).sum();
                        share += dot * dot;
                    }
                    line.push_str(&format!("  |proj onto span{{{label}}}|^2={share:.4}"));
                }
                println!("{line}");
            }
            let kept_min = dense.eigenvalues[dense.dropped];
            let kept_max = dense.eigenvalues[m - 1];
            println!(
                "  spread (kept): {:.4e} .. {:.4e}, ratio {:.3e}",
                kept_min,
                kept_max,
                kept_max / kept_min
            );

            // The matrix-free Lanczos on the same operator, as an independent
            // estimator. It sees the same extremes or one of the two is wrong.
            pre.reset_counters();
            let free = schur_spectrum(&sys, &mut pre, 80.min(m - 1), 909).expect("Schur Lanczos");
            println!(
                "  matrix-free Lanczos m={}: lambda_min={:.4e} (bracket {:.1e}) \
                 lambda_max={:.4e} (bracket {:.1e}) constant_leak={:.1e}",
                free.steps,
                free.lambda_min,
                free.bracket_min,
                free.lambda_max,
                free.bracket_max,
                free.constant_mode_leak
            );
            let gap = (free.lambda_max - kept_max).abs() / kept_max;
            assert!(
                gap < 1e-6,
                "n={n}: the dense Schur says lambda_max={kept_max:.6e} and the matrix-free \
                 Lanczos says {:.6e}, a gap of {gap:.3e}. One of the two is not the \
                 operator it claims to be",
                free.lambda_max
            );
            assert!(
                free.constant_mode_leak < 1e-14,
                "n={n}: G_imp of a constant pressure has relative size {:.3e}",
                free.constant_mode_leak
            );
        }

        // The refinement behaviour, matrix free. The dense assembly costs one
        // A-solve per CELL and stops being affordable after n = 4; the Lanczos
        // sweep costs one per STEP and its lambda_max bracket is 1e-19, so it
        // carries the h-scaling on its own.
        //
        // What is being tested: the pressure-mass Schur approximation is
        // spectrally equivalent to `S_0` only if this spectrum is BOUNDED.
        // Fact 4 of this packet says `tr(E u)` is not `D u`, so the discrete
        // inequality `|D v|_Mp <= C sqrt(Phi(v)/eta)` that would bound it is
        // not available, and whether the bound holds anyway is a measurement.
        println!();
        println!("  n  lay        h   eta*lambda_min  bracket   eta*lambda_max  bracket   innerCG");
        let mut hs: Vec<f64> = Vec::new();
        let mut lmax: Vec<f64> = Vec::new();
        for (n, layers) in ladder() {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let mut pre = BlockPreconditioner::build(
                &t,
                &a,
                &sys,
                PreconditionerOptions {
                    inner_tol: 1e-10,
                    inner_metric: Metric::FaceMass,
                    ..PreconditionerOptions::default()
                },
            )
            .expect("preconditioner builds");
            let steps = 60.min(sys.cell_count() / 2).max(1);
            let r = schur_spectrum(&sys, &mut pre, steps, 909).expect("Schur Lanczos");
            println!(
                "{n:>3} {layers:>4} {:.3e}   {:.6e}  {:.1e}   {:.6e}  {:.1e}   {}",
                sys.mesh_scale(),
                r.lambda_min,
                r.bracket_min,
                r.lambda_max,
                r.bracket_max,
                r.inner_cg_iterations
            );
            hs.push(sys.mesh_scale());
            lmax.push(r.lambda_max);
        }
        if hs.len() >= 3 {
            println!(
                "fitted over the self-similar rungs: eta*lambda_max ~ h^{:+.3}",
                slope(&hs[1..], &lmax[1..])
            );
        }
    }

    // -----------------------------------------------------------------
    // 4. The coupled solve

    /// One configuration of the coupled solve, run end to end.
    #[allow(clippy::too_many_arguments)]
    fn run_coupled(
        n: usize,
        layers: usize,
        deflation: Deflation,
        coarse: CoarseSpace,
        grad: GradientMode,
        inner_tol: f64,
        c_schur: f64,
        opts: FgmresOptions,
    ) -> (f64, usize, FgmresReport, Vec<f64>, f64) {
        let (t, a) = op(n, layers);
        let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
            .expect("system builds");
        let six = ModeBasis::build(&t, &sys, Deflation::Six).expect("six modes");
        let modes = ModeBasis::build(&t, &sys, deflation).expect("deflation basis");
        let mut pre = BlockPreconditioner::build(
            &t,
            &a,
            &sys,
            PreconditionerOptions {
                inner_tol,
                c_schur,
                coarse,
                ..PreconditionerOptions::default()
            },
        )
        .expect("preconditioner builds");
        let f = buoyancy_shaped_force(&t);
        let (_, p, rep) =
            fgmres_saddle(&sys, &mut pre, &modes, &six, grad, &f, &sys.zeros_p(), opts)
                .expect("FGMRES runs");
        let rot = ModeBasis::rotations(&t, &sys).expect("rotations");
        let leak = divergence_leak(&sys, &rot).expect("rotation leak");
        // The predicted rotation-direction residual floor:
        //   <r_j, r_u>_Mu = -(D r_j)^T M_p p  <=  |D r_j|_Mp |p|_Mp
        // divided by |b|. This is an UPPER bound on that component, and the
        // measurement is whether the achieved value tracks it.
        let f_norm = sys.residual_norm(&f, &sys.zeros_p());
        let predicted = leak.absolute[0] * sys.norm_p(&p) / f_norm;
        (
            sys.mesh_scale(),
            sys.unknown_count(),
            rep,
            leak.dimensionless,
            predicted,
        )
    }

    #[test]
    fn the_coupled_solve_across_the_refinement_ladder() {
        // Task 1 and 2: the system solved, and the block preconditioner's
        // outer count reported ALONGSIDE the inner Krylov count, which is
        // what MC-Q5 asks for and what a single number would destroy.
        //
        // Deflation::None: no velocity mode is removed from the equations.
        // The six A-null modes are deflated inside the preconditioner, where
        // they are exact, and nowhere else.
        println!(
            "  n  lay   u_unk   p_unk        h   outer  cyc   innerCG  in/out  maxIn  \
             notconv   true_rel   mom_rel  cont_rel   parity  1e-2 1e-4 1e-6 1e-8"
        );
        let mut hs: Vec<f64> = Vec::new();
        let mut outers: Vec<f64> = Vec::new();
        let mut inners: Vec<f64> = Vec::new();
        for (n, layers) in ladder() {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let parity = correction_parity(&sys, &smooth_cell_field(&sys).expect("smooth p"))
                .expect("parity");
            let (h, uunk, rep, _, _) = run_coupled(
                n,
                layers,
                Deflation::None,
                CoarseSpace::None,
                GradientMode::Corrected,
                1e-2,
                PreconditionerOptions::default().c_schur,
                // 2000 and not the default 400: at n = 8 the solve needs more
                // than 400 steps, and a capped count is not a measurement. The
                // milestone columns are the "outer iterations to a FIXED
                // relative residual" the packet asks for, so they have to be
                // reachable at every rung or the growth cannot be fitted.
                FgmresOptions {
                    max_iterations: 2000,
                    ..FgmresOptions::default()
                },
            );
            let ms = |i: usize| {
                rep.milestones[i]
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            };
            println!(
                "{n:>3} {layers:>4} {uunk:>7} {:>7} {h:.3e} {:>6} {:>4} {:>9} {:>7.1} {:>6} \
                 {:>8} {:.3e} {:.3e} {:.3e} {parity:>8.4}  {:>4} {:>4} {:>4} {:>4}",
                sys.cell_count(),
                rep.outer_iterations,
                rep.restart_cycles,
                rep.inner_cg_iterations,
                rep.inner_cg_iterations as f64 / rep.outer_iterations as f64,
                rep.inner_cg_max,
                rep.inner_cg_not_converged,
                rep.true_rel_residual,
                rep.momentum_rel,
                rep.continuity_rel,
                ms(0),
                ms(1),
                ms(2),
                ms(3)
            );
            assert!(
                (rep.rel_residual - rep.true_rel_residual).abs()
                    <= 0.5 * rep.true_rel_residual.max(rep.rel_residual) + 1e-14,
                "n={n}: the Arnoldi recurrence says {:.3e} and the recomputed residual says \
                 {:.3e}; the recurrence has drifted and the reported count is meaningless",
                rep.rel_residual,
                rep.true_rel_residual
            );
            hs.push(h);
            outers.push(rep.outer_iterations as f64);
            inners.push(rep.inner_cg_iterations as f64);
        }
        if hs.len() >= 3 {
            println!(
                "fitted over the self-similar rungs: outer ~ h^{:+.3}   innerCG ~ h^{:+.3}",
                slope(&hs[1..], &outers[1..]),
                slope(&hs[1..], &inners[1..])
            );
        }
    }

    #[test]
    fn the_deferred_correction_unrelaxed_against_relaxed_against_fully_coupled() {
        // Task 3. Three treatments of the same explicit correction, on the
        // same problem, with everything else held fixed:
        //   (a) lagged fixed point, UNRELAXED (w = 1)
        //   (b) lagged fixed point, under-relaxed
        //   (c) applied inside the Krylov operator - no fixed point at all
        // The build brief records the correction's realised norm reaching
        // 1.00 of the implicit operator's at n = 32, which is the edge of
        // where (a) converges at all. `correction_parity` re-measures that
        // here rather than carrying it forward on trust.
        println!(
            "  n  lay        h   parity(smooth) parity(random)   scheme        sweeps  \
             outer   innerCG   contraction   final_rel  verdict"
        );
        for (n, layers) in ladder() {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let six = ModeBasis::build(&t, &sys, Deflation::Six).expect("six modes");
            let none = ModeBasis::build(&t, &sys, Deflation::None).expect("empty basis");
            let smooth = smooth_cell_field(&sys).expect("smooth pressure");
            let random = pseudo_random_cells(&sys, 4242);
            let p_smooth = correction_parity(&sys, &smooth).expect("parity");
            let p_random = correction_parity(&sys, &random).expect("parity");
            let f = buoyancy_shaped_force(&t);

            // The lagged schemes run only at n <= 4, and the reason is a
            // MEASUREMENT rather than cost. A fixed point's contraction factor
            // is only the fixed point's contraction factor if each sweep's
            // inner solve actually converges. At n = 8 the implicit FGMRES
            // hits its 400-step cap on EVERY sweep (measured: 6000 outer steps
            // over 15 sweeps, i.e. 400 every time), so the sweep map is not the
            // map whose spectral radius is wanted, and the 0.824 it reported
            // would be an artefact of the cap rather than a property of the
            // deferred correction. Rather than publish a polluted number the
            // lagged rows stop where the inner solve converges, and the parity
            // column - cheap, clean, and the quantity the build brief's
            // "1.00 at n = 32" refers to - is carried at every rung instead.
            let relaxations: &[f64] = if n <= 4 { &[1.0, 0.5, 0.25] } else { &[] };
            for relax in relaxations.iter().copied() {
                let mut pre =
                    BlockPreconditioner::build(&t, &a, &sys, PreconditionerOptions::default())
                        .expect("preconditioner builds");
                let (_, _, rep) = deferred_correction_solve(
                    &sys,
                    &mut pre,
                    &none,
                    &six,
                    &f,
                    &sys.zeros_p(),
                    DeferredOptions {
                        relax,
                        max_sweeps: 15,
                        ..DeferredOptions::default()
                    },
                )
                .expect("deferred correction runs");
                let verdict = if rep.diverged {
                    "DIVERGED"
                } else if rep.converged {
                    "converged"
                } else {
                    "stalled/capped"
                };
                println!(
                    "{n:>3} {layers:>4} {:.3e}  {p_smooth:>13.5} {p_random:>14.5}   \
                     lagged w={relax:<5} {:>5} {:>7} {:>9} {:>13.5} {:.3e}  {verdict}",
                    sys.mesh_scale(),
                    rep.sweeps,
                    rep.outer_iterations,
                    rep.inner_cg_iterations,
                    rep.contraction,
                    rep.final_rel_residual
                );
            }

            let mut pre =
                BlockPreconditioner::build(&t, &a, &sys, PreconditionerOptions::default())
                    .expect("preconditioner builds");
            let (_, _, rep) = fgmres_saddle(
                &sys,
                &mut pre,
                &none,
                &six,
                GradientMode::Corrected,
                &f,
                &sys.zeros_p(),
                FgmresOptions::default(),
            )
            .expect("FGMRES runs");
            println!(
                "{n:>3} {layers:>4} {:.3e}  {p_smooth:>13.5} {p_random:>14.5}   \
                 coupled(Krylov)  {:>5} {:>7} {:>9} {:>13}  {:.3e}  {}",
                sys.mesh_scale(),
                "-",
                rep.outer_iterations,
                rep.inner_cg_iterations,
                "-",
                rep.true_rel_residual,
                if rep.converged {
                    "converged"
                } else {
                    "NOT converged"
                }
            );
        }
    }

    #[test]
    fn the_near_null_rotations_and_what_deflating_them_costs() {
        // Task 4, and the reason for the packet's amendment. Four velocity
        // treatments of the SAME problem:
        //
        //   none          - nothing removed. The reported solve. Every
        //                   equation is kept, including the three
        //                   rotation-tested momentum rows whose coefficients
        //                   are O(h^2).
        //   rotations(3)  - the near-null rows DELETED. This solves a
        //                   DIFFERENT discrete system and the number below is
        //                   the residual of that other system.
        //   translations  - the strongly non-null rows deleted; also a
        //                   different system, and a control that separates
        //                   "deflation helps" from "deleting equations helps".
        //   six           - both.
        //   none + coarse - nothing deleted, and the three rotations SUPPLIED
        //                   to the preconditioner as a near-nullspace. This is
        //                   the treatment that is legitimate for a mode that
        //                   is not exactly null.
        println!(
            "  n  lay        h  h|D r|      deflation        coarse       outer  innerCG  \
             true_rel   mom_rel  cont_rel   rot_resid  tra_resid  predicted_floor"
        );
        for (n, layers) in ladder() {
            for (deflation, coarse) in [
                (Deflation::None, CoarseSpace::None),
                (Deflation::Rotations, CoarseSpace::None),
                (Deflation::Translations, CoarseSpace::None),
                (Deflation::Six, CoarseSpace::None),
                (Deflation::None, CoarseSpace::Rotations),
                (Deflation::None, CoarseSpace::Translations),
                (Deflation::None, CoarseSpace::Six),
            ] {
                let (h, _, rep, leak, predicted) = run_coupled(
                    n,
                    layers,
                    deflation,
                    coarse,
                    GradientMode::Corrected,
                    1e-2,
                    PreconditionerOptions::default().c_schur,
                    FgmresOptions::default(),
                );
                println!(
                    "{n:>3} {layers:>4} {h:.3e} {:.2e} {:>15} {:>13} {:>7} {:>8} {:.3e} \
                     {:.3e} {:.3e} {:.3e} {:.3e} {predicted:.3e}{}",
                    leak[0],
                    deflation.label(),
                    coarse.label(),
                    rep.outer_iterations,
                    rep.inner_cg_iterations,
                    rep.true_rel_residual,
                    rep.momentum_rel,
                    rep.continuity_rel,
                    rep.rotation_residual(),
                    rep.translation_residual(),
                    if deflation.changes_the_system() {
                        "   [DIFFERENT SYSTEM]"
                    } else {
                        ""
                    }
                );
            }
        }
        // The coarse spaces' own conditioning, which is the "measure their
        // eigenvalues" half of the amendment. These are the eigenvalues of
        // `h^2 <D n_j, D n_k>_{M_p}` with the modes M_u-normalised, so for the
        // rotations they are `(h |D r|)^2` and directly comparable with the
        // leak table above.
        println!();
        println!(
            "  n  lay  coarse space      eig(h^2 <D n_j, D n_k>_Mp) ascending                               condition"
        );
        for (n, layers) in ladder() {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            for which in [
                CoarseSpace::Rotations,
                CoarseSpace::Translations,
                CoarseSpace::Six,
            ] {
                let rc = NullCoarse::build(&t, &sys, which).expect("coarse space");
                let e = rc.eigenvalues();
                let shown: Vec<String> = e.iter().map(|v| format!("{v:.4e}")).collect();
                println!(
                    "{n:>3} {layers:>4}  {:<15}  {}   {:.4e}",
                    which.label(),
                    shown.join(" "),
                    e[e.len() - 1] / e[0]
                );
                assert!(
                    e[0] > 0.0,
                    "n={n} {}: the coarse Gram matrix has a non-positive eigenvalue \
                     {:.3e}; the modes are linearly dependent in the divergence and the \
                     coarse space is smaller than it claims",
                    which.label(),
                    e[0]
                );
            }
        }
    }

    // -----------------------------------------------------------------
    // 5. Self-consistency, NOT physics

    #[test]
    fn the_manufactured_saddle_pair_is_recovered() {
        // SELF-CONSISTENCY, NOT PHYSICS. The pair is manufactured in the
        // DISCRETE space and the right-hand side is formed by applying THESE
        // operators to it, so this measures whether the solver inverts the
        // operator this repository built. It says nothing about whether that
        // operator is the right physics; the analytic comparison is MC-Q1
        // against `mantle::benchmark_thieulot` and is the NEXT packet.
        //
        // Three rows per rung, and the difference between them is the whole
        // content:
        //
        //  in V_0   - the manufactured velocity carries no A-null component.
        //             The solver can represent it, and this is the row that
        //             says whether the solver inverts the operator.
        //  full     - the manufactured velocity carries all six. The
        //             preconditioner's A solve is deflated, so the Krylov
        //             velocity lives in V_0 and those six components are
        //             unreachable; the row measures what that costs, per
        //             direction, and where the residual then stops.
        //  full/six - the same solution with the six momentum rows DELETED.
        //             A DIFFERENT discrete system, reported as a control.
        //
        // Tolerance 1e-12 rather than the ladder's 1e-8, for the standing
        // reason: a run stopped at residual `tol` has SOLUTION error up to
        // about kappa * tol, so a solution-error check has to buy that
        // headroom. Capped at n <= 8 on COST.
        println!(
            "  n  lay  solution  equations  outer  innerCG   true_rel   verdict     \
             |du|/|u|   |dp|/|p|  |du_rot|/|u| |du_tra|/|u| |du_rest|/|u|"
        );
        for (n, layers) in ladder() {
            if n > 8 {
                println!("  (n = {n} skipped: cost, see the doc comment)");
                continue;
            }
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds");
            let six = ModeBasis::build(&t, &sys, Deflation::Six).expect("six modes");
            let none = ModeBasis::build(&t, &sys, Deflation::None).expect("empty basis");
            for (label, restrict, equations) in [
                ("in V_0", &six, &none),
                ("full", &none, &none),
                ("full", &none, &six),
            ] {
                let pair =
                    manufactured_pair(&sys, restrict, equations, GradientMode::Corrected, 77)
                        .expect("manufactured pair");
                let mut pre = BlockPreconditioner::build(
                    &t,
                    &a,
                    &sys,
                    PreconditionerOptions {
                        inner_tol: 1e-4,
                        ..PreconditionerOptions::default()
                    },
                )
                .expect("preconditioner builds");
                let (u, p, rep) = fgmres_saddle(
                    &sys,
                    &mut pre,
                    equations,
                    &six,
                    GradientMode::Corrected,
                    &pair.rhs_u,
                    &pair.rhs_p,
                    FgmresOptions {
                        // 1e-11 and not 1e-12: the recomputed residual bottoms
                        // out at 1.6e-12 on this problem (measured), so 1e-12
                        // sits ON the roundoff floor and "converged" would be
                        // reporting the floor. The solution errors below are
                        // the measurement; the residual target only has to be
                        // small enough that kappa * tol leaves them meaningful,
                        // and kappa(A) is 3.1e2 here.
                        tol: 1e-11,
                        max_iterations: 800,
                        restart: 100,
                        // A wider stall window than the default, so that a
                        // plateau is reported as a plateau rather than as an
                        // early stop that could be mistaken for one.
                        stall_window: 120,
                        ..FgmresOptions::default()
                    },
                )
                .expect("FGMRES runs");
                let du: Vec<f64> = u.iter().zip(pair.u.iter()).map(|(x, y)| x - y).collect();
                let mut dp: Vec<f64> = p.iter().zip(pair.p.iter()).map(|(x, y)| x - y).collect();
                sys.remove_mean_pressure(&mut dp).expect("mean free");
                let scale = sys.norm_u(&pair.u);
                let eu = sys.norm_u(&du) / scale;
                let ep = sys.norm_p(&dp) / sys.norm_p(&pair.p);
                let comps: Vec<f64> = (0..NULL_DIM).map(|j| sys.dot_u(six.mode(j), &du)).collect();
                let tra = (comps[0..3].iter().map(|c| c * c).sum::<f64>()).sqrt() / scale;
                let rot = (comps[3..6].iter().map(|c| c * c).sum::<f64>()).sqrt() / scale;
                let mut rest = du.clone();
                six.project(&sys, &mut rest);
                let rest = sys.norm_u(&rest) / scale;
                let verdict = if rep.converged {
                    "converged"
                } else if rep.stalled {
                    "STALLED"
                } else {
                    "capped"
                };
                println!(
                    "{n:>3} {layers:>4}  {label:<8}  {:<9} {:>6} {:>8} {:.3e}  {verdict:<10} \
                     {eu:.4e}  {ep:.4e}  {rot:.4e}   {tra:.4e}   {rest:.4e}{}",
                    equations.which().label(),
                    rep.outer_iterations,
                    rep.inner_cg_iterations,
                    rep.true_rel_residual,
                    if equations.which().changes_the_system() {
                        "  [DIFFERENT SYSTEM]"
                    } else {
                        ""
                    }
                );
                if label == "in V_0" {
                    // The one row that is a genuine pass/fail: the solution
                    // lies in the space the solver can represent, so nothing
                    // structural stands between it and the answer.
                    // The assertion is on the SOLUTION error and not on the
                    // residual: the residual has a measured floor on this
                    // problem and asserting below a floor would be asserting
                    // on roundoff. It is also CONDITIONAL on the solve having
                    // got somewhere, because a solution error is only a
                    // statement about the solver once the residual is small.
                    //
                    // Measured floors, printed in the `true_rel` and `verdict`
                    // columns rather than tuned around: 1.6e-12 at n = 2,
                    // 6.7e-10 at n = 4, and at n = 8 the solve STALLS at
                    // 3.4e-1 and never gets started. That last one is a real
                    // result and section 8 of the findings reports it as one -
                    // the manufactured right-hand side carries a nonzero
                    // continuity source `D u_ex`, which is a far harder
                    // problem for this preconditioner than the driven case
                    // with `g = 0`, and at the h^-2 outer growth of section 5
                    // the budget here does not reach it.
                    if rep.true_rel_residual < 1e-8 {
                        assert!(
                            eu < 1e-5 && ep < 1e-5,
                            "n={n}: |du|/|u| = {eu:.3e} and |dp|/|p| = {ep:.3e} at a \
                             residual of {:.3e}; the solver does not invert its own \
                             operator",
                            rep.true_rel_residual
                        );
                    } else {
                        println!(
                            "      ^ n={n} NOT CONVERGED ({:.3e}): the solution errors on \
                             this row are not a statement about the solver, only about \
                             where it stopped",
                            rep.true_rel_residual
                        );
                        // The rungs where it DOES converge are what give this
                        // test teeth; if those stopped converging the test
                        // would silently become vacuous.
                        assert!(
                            n > 4,
                            "n={n}: the V_0 manufactured solve reached only {:.3e}. At \
                             n <= 4 it converges to 1e-11, so a stall here is a \
                             regression, not the known fine-rung cost",
                            rep.true_rel_residual
                        );
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // 6. The two scales that are NOT tuning knobs

    #[test]
    fn the_solve_does_not_depend_on_the_viscosity_scale() {
        // eta enters the operator only as `W_c = 2 eta V_c` and enters the
        // solver only through `theta = eta / h` and `alpha = eta / c_schur`,
        // all three of which cancel out of the preconditioned operator. So in
        // exact arithmetic these runs are the SAME iteration and their residual
        // histories are the same sequence.
        //
        // WHY THIS TEST COMPARES HISTORIES AND NOT ITERATION COUNTS. An earlier
        // version compared final counts and got 84 / 110 / 85 at
        // eta = 1e-6 / 1 / 1e21, and its assertion message blamed "a hidden
        // absolute scale". THAT MESSAGE WAS WRONG, and its own numbers refute
        // it: a hidden absolute scale is monotone in eta, and this was
        // symmetric about eta = 1 with both extremes - 27 decades apart -
        // agreeing to one iteration. Two candidate causes were separated by
        // measurement rather than by argument:
        //
        //   (a) an absolute tolerance or floor binding near eta = 1. Ruled out
        //       twice: every threshold on this path is relative or a sign test
        //       (the only absolute constants in the module, the 1e-300 floors
        //       in `jacobi_symmetric`, are reached only through the coarse
        //       space and the dense Schur, neither of which is enabled here),
        //       and empirically the counts below are EQUAL, not merely close,
        //       once (b) is removed. A binding absolute threshold would still
        //       bind.
        //   (b) THE INNER SOLVE'S NULL LEAK. `deflated_cg` runs in the
        //       `diag(K)` metric and returns a vector orthogonal to the six
        //       modes IN THAT METRIC, which is not `M_u`-orthogonality. The
        //       residue sits in the near-null rotation directions, where the
        //       outer method amplifies it, and the resulting iteration is
        //       genuinely unstable rather than merely roundoff-sensitive. The
        //       84 / 110 / 85 run predates the `M_u` re-projection that now
        //       removes it.
        //
        // The last block is the NEGATIVE CONTROL that decides between them: it
        // turns the re-projection off and changes nothing else. If (b) is the
        // cause, the eta dependence comes back there and only there.
        //
        // Two things that were CHECKED AND ARE NOT the cause, recorded so the
        // next reader does not re-derive them:
        //
        //   - Restarting. It was the first hypothesis and it is wrong: with
        //     `restart = 20`, i.e. six restarts in a 134-step run, the outer
        //     counts are still identical across all three etas.
        //   - The inner CG's stopping ambiguity. Tightening `inner_tol` from
        //     1e-2 to 1e-10 brings the inner counts to within one iteration of
        //     each other (7351 / 7350 / 7350) and the histories STILL part at
        //     step 15. So the parting is ordinary floating-point divergence of
        //     the Arnoldi recurrence, and it does not move the outer count.
        //
        // The invariance claim itself is the outer iteration count, and it is
        // asserted as EQUALITY, not as a ratio.
        let (n, layers) = (4usize, 3usize);
        let etas = [1e-6f64, 1.0, 1e21];

        let run =
            |eta: f64, restart: usize, inner_tol: f64, reproject: bool| -> (FgmresReport, f64) {
                let t = mini(n, layers);
                let a = ViscousOperator::isoviscous(&t, eta).expect("operator builds");
                let sys = SaddleSystem::build(&t, &a, eta, FaceGradientRule::NonOrthogonalAndSkew)
                    .expect("system builds");
                let six = ModeBasis::build(&t, &sys, Deflation::Six).expect("six modes");
                let none = ModeBasis::build(&t, &sys, Deflation::None).expect("empty basis");
                let mut pre = BlockPreconditioner::build(
                    &t,
                    &a,
                    &sys,
                    PreconditionerOptions {
                        inner_tol,
                        reproject_inner_in_m_u: reproject,
                        ..PreconditionerOptions::default()
                    },
                )
                .expect("preconditioner builds");
                let f = buoyancy_shaped_force(&t);
                let (_, _, rep) = fgmres_saddle(
                    &sys,
                    &mut pre,
                    &none,
                    &six,
                    GradientMode::Corrected,
                    &f,
                    &sys.zeros_p(),
                    FgmresOptions {
                        restart,
                        ..FgmresOptions::default()
                    },
                )
                .expect("FGMRES runs");
                let leak = pre.worst_inner_leak();
                (rep, leak)
            };

        // First step at which two residual histories differ by more than
        // `1e-9` relative. `None` means they never do.
        let part = |a: &[f64], b: &[f64]| -> Option<usize> {
            a.iter()
                .zip(b.iter())
                .position(|(x, y)| (x - y).abs() > 1e-9 * x.abs().max(y.abs()))
        };

        for (restart, inner_tol, reproject, label) in [
            (400usize, 1e-2f64, true, "no restart, inner 1e-2"),
            (400, 1e-10, true, "no restart, inner 1e-10"),
            (20, 1e-2, true, "restart 20, inner 1e-2"),
            (60, 1e-2, false, "[NEG-CONTROL] no reproj"),
        ] {
            let pairs: Vec<(FgmresReport, f64)> = etas
                .iter()
                .map(|e| run(*e, restart, inner_tol, reproject))
                .collect();
            let leaks: Vec<f64> = pairs.iter().map(|p| p.1).collect();
            let reports: Vec<FgmresReport> = pairs.into_iter().map(|p| p.0).collect();
            for (eta, rep) in etas.iter().zip(reports.iter()) {
                println!(
                    "{label:<24} eta={eta:>8.1e}  outer={:<4} innerCG={:<6} true_rel={:.6e} \
                     converged={}",
                    rep.outer_iterations,
                    rep.inner_cg_iterations,
                    rep.true_rel_residual,
                    rep.converged
                );
            }
            println!("{label:<24} worst M_u null leak of an inner solve: {leaks:.3?}");
            let firsts: Vec<Option<usize>> = (1..etas.len())
                .map(|k| part(&reports[0].history, &reports[k].history))
                .collect();
            println!(
                "{label:<24} outer counts {:?}; inner CG counts {:?}; histories part at \
                 step {:?} of {}",
                reports
                    .iter()
                    .map(|r| r.outer_iterations)
                    .collect::<Vec<_>>(),
                reports
                    .iter()
                    .map(|r| r.inner_cg_iterations)
                    .collect::<Vec<_>>(),
                firsts,
                reports[0].history.len()
            );
            // Every run must reach the SAME stopping criterion, or "the
            // extremes look better" would only mean they stopped earlier
            // against a laxer test. Checked explicitly rather than inferred
            // from the counts.
            for (eta, rep) in etas.iter().zip(reports.iter()) {
                assert!(
                    rep.converged && rep.true_rel_residual <= FgmresOptions::default().tol,
                    "eta={eta:.1e} ({label}) stopped at {:.6e} against a requested {:.1e}: \
                     the runs are not being held to the same criterion",
                    rep.true_rel_residual,
                    FgmresOptions::default().tol
                );
            }
            if reproject {
                // THE invariance claim: the same number of outer iterations,
                // not merely a similar number, across 27 decades of eta.
                for (eta, rep) in etas.iter().zip(reports.iter()) {
                    assert_eq!(
                        rep.outer_iterations, reports[0].outer_iterations,
                        "eta={eta:.1e} ({label}): {} outer iterations against {} at \
                         eta=1e-6. eta cancels out of the preconditioned operator, so a \
                         difference here is a real scale leak",
                        rep.outer_iterations, reports[0].outer_iterations
                    );
                }
            } else {
                // The negative control. It is allowed to differ - that is the
                // point - and what is asserted is that the leak it reports is
                // real, so the row is evidence rather than an anecdote.
                let worst = leaks.iter().fold(0.0f64, |a, b| a.max(*b));
                assert!(
                    worst > 1e-6,
                    "the negative control reported a worst M_u null leak of {worst:.3e}, \
                     i.e. no leak. Then it is not a control for anything and the \
                     explanation above is unsupported"
                );
            }
        }
    }
    #[test]
    fn the_residual_scale_is_not_a_tuning_knob() {
        // theta = eta / h is derived from unit consistency, and it changes
        // the norm the residual is minimised in, so it changes the iterates.
        // The measurement is how much: a scale that had to be tuned would
        // show a sharp optimum.
        println!("  theta multiplier   outer   innerCG    true_rel    mom_rel   cont_rel");
        let (n, layers) = (4usize, 3usize);
        for mult in [1e-2f64, 1e-1, 1.0, 1e1, 1e2] {
            let (t, a) = op(n, layers);
            let sys = SaddleSystem::build(&t, &a, ETA, FaceGradientRule::NonOrthogonalAndSkew)
                .expect("system builds")
                .with_residual_scale(mult)
                .expect("scaled system");
            let six = ModeBasis::build(&t, &sys, Deflation::Six).expect("six modes");
            let none = ModeBasis::build(&t, &sys, Deflation::None).expect("empty basis");
            let mut pre =
                BlockPreconditioner::build(&t, &a, &sys, PreconditionerOptions::default())
                    .expect("preconditioner builds");
            let f = buoyancy_shaped_force(&t);
            let (_, _, rep) = fgmres_saddle(
                &sys,
                &mut pre,
                &none,
                &six,
                GradientMode::Corrected,
                &f,
                &sys.zeros_p(),
                FgmresOptions::default(),
            )
            .expect("FGMRES runs");
            println!(
                "{mult:>17.0e} {:>7} {:>9}  {:.3e}  {:.3e}  {:.3e}",
                rep.outer_iterations,
                rep.inner_cg_iterations,
                rep.true_rel_residual,
                rep.momentum_rel,
                rep.continuity_rel
            );
        }
    }

    #[test]
    fn the_schur_scaling_is_not_a_tuning_knob() {
        // c_schur is read off the measured Schur spectrum ONCE and used at
        // every rung. This sweeps it over two decades to show what a wrong
        // value costs, so that the value used is visible as a choice rather
        // than as a fit.
        println!("  c_schur   outer   innerCG    true_rel   converged");
        let (n, layers) = (4usize, 3usize);
        for c in [0.05f64, 0.15, 0.5, 1.5, 5.0] {
            let (h, _, rep, _, _) = run_coupled(
                n,
                layers,
                Deflation::None,
                CoarseSpace::None,
                GradientMode::Corrected,
                1e-2,
                c,
                FgmresOptions::default(),
            );
            let _ = h;
            println!(
                "{c:>9.3} {:>7} {:>9}  {:.3e}   {}",
                rep.outer_iterations, rep.inner_cg_iterations, rep.true_rel_residual, rep.converged
            );
        }
    }

    #[test]
    fn the_inner_tolerance_trades_outer_count_for_krylov_work() {
        // MC-Q5 wants both numbers. This is the trade, measured: a looser
        // inner solve costs outer iterations and buys inner ones, and the
        // product is what a runtime estimate has to be built from.
        println!("  inner_tol   outer   innerCG   inner/outer   maxInner   true_rel");
        let (n, layers) = (4usize, 3usize);
        for tol in [3e-1f64, 1e-1, 1e-2, 1e-4, 1e-6] {
            let (_, _, rep, _, _) = run_coupled(
                n,
                layers,
                Deflation::None,
                CoarseSpace::None,
                GradientMode::Corrected,
                tol,
                PreconditionerOptions::default().c_schur,
                FgmresOptions::default(),
            );
            println!(
                "{tol:>11.0e} {:>7} {:>9} {:>13.1} {:>10} {:.3e}",
                rep.outer_iterations,
                rep.inner_cg_iterations,
                rep.inner_cg_iterations as f64 / rep.outer_iterations as f64,
                rep.inner_cg_max,
                rep.true_rel_residual
            );
        }
    }
}
