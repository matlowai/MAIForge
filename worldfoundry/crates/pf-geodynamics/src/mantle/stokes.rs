//! MC0: the isoviscous viscous operator of the Stokes system, built as the
//! Gram form of a discrete strain-rate operator. **There is no solver here.**
//!
//! ```text
//! Phi(u) = sum_c  2 eta V_c  eps_c(u) : eps_c(u)        discrete strain energy
//! A      = -M_u^-1 E^T W E                              the viscous operator
//! ```
//!
//! # Why the energy form, and not a term-by-term stress difference
//!
//! Brief rev 5 "MC design and qualification" forbids a naive two-point
//! viscous Laplacian, and the cubed-sphere review
//! (`docs/research/2026-09-07-cubed-sphere-stokes-operators/FINDINGS.md`
//! section 3.1) says why the alternative is not a matter of taste. Written as
//! the Gram form of ANY linear `E`, with ANY positive-definite `W`, four
//! properties are unconditional - they hold for every interpolation choice
//! inside `E`, including whatever happens at a cube seam:
//!
//! 1. `M_u A = -E^T W E` is symmetric, so `A` is self-adjoint in the `M_u`
//!    inner product and the block system `[-A, G; D, 0]` is a genuine
//!    symmetric saddle point.
//! 2. `u . M_u A u = -Phi(u) <= 0` for every `u`: the discretization cannot
//!    manufacture energy, whatever the mesh does.
//! 3. `null(A) = null(E)`. Nothing else can enter the null space, so if `E`
//!    annihilates rigid motions the viscous term can never manufacture spin.
//! 4. `D` and `G` are untouched: the compatibility identity
//!    `G = -M_u^-1 D^T M_p` that [`pf_grid::shell_ops`] proves is unaffected
//!    by anything done here.
//!
//! Differencing the stress tensor term by term gives none of the four, and on
//! a non-orthogonal grid there is no cheap way to recover them afterwards.
//!
//! The four solver-safety gates raised against this operator in review -
//! local conditioning and coercivity, the volumetric mismatch, the null-space
//! statement, and the certification of the rotation residual - are measured in
//! the `gate1_`..`gate4_` tests below and written up in
//! `docs/reviews/planet-foundry-b1/MC0_CONDITIONING_GATES.md`.
//!
//! # The strain operator `E`, and why it is a per-cell affine fit
//!
//! The unknowns are face-normal scalars `u_f` (the facet's AREA-MEAN normal
//! component, matching [`pf_grid::shell_ops::ShellTopology::rigid_rotation_field`]).
//! A cell has six of them and a velocity gradient has nine components, so no
//! construction can build a cell strain rate from that cell alone: `E` must
//! reach into the neighbours. That is a property of the unknown layout, not
//! of this implementation.
//!
//! # THE WIDE STENCIL IS STRUCTURAL, NOT AN IMPLEMENTATION DETAIL
//!
//! Stated here because a reader would otherwise reasonably assume that the
//! per-cell construction is what makes this work, and it is not. Measured in
//! `tests::gate3_the_single_cell_moment_map_is_rank_three_and_the_defect_is_real`
//! and written up in `docs/reviews/planet-foundry-b1/MC0_CONDITIONING_GATES.md`:
//!
//! - On a box, the map from the six constant strain modes to the six
//!   face-normal INTEGRATED fluxes has exactly three zero singular values.
//!   All three shears - `S12`, `S13`, `S23` - are invisible. Exact facet
//!   integration does not help: the integral of `y` over the `x = +-a` face
//!   is exactly zero. One normal-flux scalar per face on a six-face cell
//!   does not contain the information.
//! - On our cubed-sphere cells that map has exactly ONE singular value at
//!   roundoff (relative `1.6e-16` at worst over `n = 4 .. 32`) and two more
//!   at `1.1e-2` to `1.9e-2` relative, decaying like `h^0.7` to `h^0.9`. So
//!   curvature buys back two of the three shears, weakly, and refinement
//!   takes them away again.
//!
//! **A per-cell face-normal-only viscous strain operator therefore cannot
//! exist on this degree-of-freedom set.** The Brezzi/Lipnikov/Shashkov
//! per-cell `M_E` framework remains valid for the pressure/flux pair, where
//! one scalar per face is the right amount of information; it does not carry
//! over to the constitutive strain operator, and that is the concrete reason
//! the curved-face MFD literature enriches cells with tangential components.
//!
//! What makes this operator work is the STENCIL: cell plus face neighbours,
//! 28 to 34 faces. A neighbour face carries a lever arm of order `2 h` in a
//! direction transverse to its own normal, which is exactly the information a
//! shear needs. Measured on the operator's own stencil, the smallest
//! rigid-projected strain singular value is `1.06e-1` to `2.31e-1` and does
//! not move with refinement (`h^-0.10` to `h^+0.15` over `n = 4, 8, 16, 32`).
//! **Narrowing the stencil to the cell would not degrade this operator, it
//! would destroy it.**
//!
//! What is chosen here is the *exactness class*. `E` fits the affine field
//!
//! ```text
//! v(x) = b + (G / h_c) (x - x_c)
//! ```
//!
//! to the face data of the cell and of its face neighbours by weighted least
//! squares, and returns `eps_c = sym(G) / h_c`. The row of the design matrix
//! for face `f` is the EXACT functional the datum `u_f` is:
//!
//! ```text
//! u_f = (1 / |A_f|) [ A_f . b + (1/h_c) sum_de G[d][e] (T_f[d][e] - A_f[d] x_c[e]) ]
//! ```
//!
//! with `T_f[a][d] = integral_f n_a x_d dA` the exact facet tensor moment
//! ([`pf_grid::shell_ops::ShellTopology::facet_tensor_moments`]). Because the
//! row is exact, the fit has zero residual on any affine field and therefore
//! returns its gradient exactly - on any mesh, for any positive weights, in
//! every cell including cube corners and seam-touching cells. Everything the
//! test suite of this module measures follows from that one property:
//!
//! - a uniform translation has `G = 0`, so `E u = 0` and `A u = 0`;
//! - a rigid rotation `omega x x` has `G` exactly antisymmetric, so
//!   `sym(G) = 0` and `A u = 0` **to roundoff, independent of `h`** - not to
//!   `O(h^2)`, which is what a point-sampled facet moment would give
//!   (measured in `shell_ops`: 5.6e-3 at n=8, still 4.6e-4 at n=32);
//! - a constant strain rate is reproduced exactly, again everywhere.
//!
//! **Deviation from the review's section 3.2, stated.** That section sketches
//! a two-stage `E`: reconstruct a cell velocity vector with
//! [`pf_grid::shell_ops::ShellTopology::reconstruct_cell_velocity`], build a
//! face vector from it, then Green-Gauss. That reconstruction is exact only
//! for a CONSTANT field, so on a rigid rotation - which is not constant
//! inside a cell - the face vector is wrong at `O(h)` and the resulting
//! strain does not vanish at machine zero. The exactness class of the whole
//! operator is set by its weakest stage. A single affine fit is the cheapest
//! stage that is exact for the fields whose annihilation the brief demands,
//! so it is what this module does. The structural rule of section 3.1 is
//! unaffected: this is still `A = -M_u^-1 E^T W E`.
//!
//! # Where the polyhedral volume went
//!
//! [`pf_grid::shell_ops::ShellTopology::polyhedral_cell_volumes_m3`]
//! (`V_poly = trace(sum_f s(c,f) T_f) / 3`) is the volume a Green-Gauss
//! gradient must divide by, because the isotropy identity
//! `sum_f s(c,f) T_f = V_poly I` is what makes such a gradient exact. This
//! `E` never forms a Green-Gauss sum, so no volume appears as a divisor; the
//! exactness comes directly from the exact `T_f` in the design rows. `V_poly`
//! still appears twice, both times honestly:
//!
//! 1. as the fit's length scale `h_c = V_poly^(1/3)`, which only conditions
//!    the normal equations and cannot affect the answer;
//! 2. as the reason the discrete strain trace and the discrete divergence are
//!    NOT the same operator. For an affine field the exact facet fluxes sum to
//!    `V_poly tr(G/h)`, while [`pf_grid::shell_ops::ShellTopology::divergence`]
//!    divides that sum by the spherical frustum volume, so
//!    `tr(eps_c) = (V_sphere / V_poly) (D u)_c`, a ratio of `1 + O(h^2)`
//!    (1.9e-2 at n=8, 1.3e-3 at n=32, measured in `shell_ops`). It is
//!    measured here too, in
//!    `tests::the_strain_trace_and_the_discrete_divergence_differ_by_the_volume_ratio`,
//!    because the solver has to know that its two "incompressibility" numbers
//!    are not the same number.
//!
//!    **That relation holds on AFFINE data only and is not an operator
//!    identity.** `E` reads 28 to 34 faces and `D` reads six; they agree
//!    because a zero-residual fit reproduces the Green-Gauss flux sum, not
//!    because one is a multiple of the other. On random data the two differ
//!    by `7.4e-1` of the trace scale
//!    (`tests::gate2_the_trace_identity_holds_on_affine_data_only_and_is_not_an_operator_identity`).
//!    So `D u = 0` does NOT force `tr(E u) = 0`, and the strain energy does
//!    carry a volumetric penalty the pressure cannot see. Measured on
//!    manufactured exactly-solenoidal fields, that penalty is a convergent
//!    truncation error and not a bulk viscosity: it is at most `6.7e-5` of
//!    the viscous energy at `n = 4` and falls at `h^2.0` to `h^3.3`
//!    (`tests::gate2_the_volumetric_mismatch_is_a_convergent_truncation_not_a_bulk_viscosity`).
//!
//! The energy weight uses the SPHERICAL volume `V_c`, the same geometry as
//! `M_p` and as the buoyancy forcing that will drive it. Nothing structural
//! depends on that choice - any positive `W` keeps all four properties above -
//! and it is recorded rather than left implicit.
//!
//! # Free slip, as the natural condition
//!
//! Impermeability (`u . n = 0`) is an ESSENTIAL constraint on the velocity
//! SPACE - not a property of `D`, which annihilates a rigid translation just
//! as happily as the continuum divergence does. Here it is imposed by
//! EXCLUDING the boundary faces from the unknown vector, and this module
//! additionally never reads a boundary face: they are absent from every fit
//! stencil.
//! Zero tangential traction is then the NATURAL condition, which means
//! concretely: the energy `Phi` is a sum over cells with no surface term, so
//! the stationarity condition `dPhi/du_f = 0` that defines `A` produces no
//! boundary contribution, and nothing constrains the tangential velocity at
//! the wall. Doing nothing IS free slip. Two consequences worth stating,
//! since "natural" is often said and rarely cashed out:
//!
//! - `A u` is unchanged by any value written into a boundary entry of `u`
//!   (`tests::boundary_face_values_never_enter_the_operator`), so there is no
//!   ghost value to choose and no ghost-value convention to get wrong;
//! - the tangential slip modes at the wall are NOT in the null space - they
//!   cost energy through the strain of the boundary cells.
//!
//! # The null space, as it actually is in THIS space
//!
//! Measured in
//! `tests::gate1_the_null_space_is_six_and_three_of_its_directions_are_not_rigid_motions`.
//! Our unknown vector has already had the boundary degrees of freedom
//! removed, and the null space of `A` is nevertheless SIX dimensional:
//!
//! - three ROTATIONS, which are rigid and are impermeable on the sphere to
//!   `O(h^2)` (their exact normal component on a boundary FACET is `1.06e-2`
//!   of their interior scale at `n = 4` and `2.80e-3` at `n = 8`, order
//!   `+1.92` - the flat facet's departure from the sphere);
//! - three TRUNCATED TRANSLATIONS: a uniform translation on the interior
//!   faces and zero on the boundary faces. That is not a rigid motion. The
//!   datum it was truncated at is `9.65e-1` of its own interior scale, so it
//!   is not even close to one, and it is not divergence free: `D u` is
//!   `1.41` (normalised) in every boundary cell. `A` cannot tell, only
//!   because `A` never reads a boundary face.
//!
//! `D` therefore has rank SIX on `null(A)`, not three, with a two-order gap:
//! the three translation directions are `O(1)` and the three rotation
//! directions are `1.04e-2` at `n = 4`, falling at order `+2.09`. So even the
//! rotations are only APPROXIMATELY in the null space of the saddle point -
//! for a Krylov solver a near-null mode is the more awkward of the two,
//! because it cannot be deflated away.
//!
//! **Consequence for the solver: a multigrid preconditioner on the `A` block
//! alone must deflate SIX modes.** Our code is already in the
//! after-restriction regime and still has six; the count is a property of the
//! free-slip stencil exclusion, not of the boundary bookkeeping.
//!
//! This is a property of the `eps` form and NOT of the `eta grad^2 u` form,
//! whose natural condition is zero normal derivative of each Cartesian
//! velocity component - which on a curved boundary is not zero shear
//! traction. That is a second, independent reason to build the operator from
//! a strain energy (review section 3.4).
//!
//! # Isoviscous, deliberately
//!
//! One constant `eta` (MC0). Variable viscosity is MC2 and later, and the
//! seam it will be inserted at is already there: `eta` enters ONLY as the
//! per-cell weight `W_c = 2 eta V_c` in [`ViscousOperator::from_cell_weights`],
//! never inside the strain assembly. A later packet replaces one vector and
//! touches no geometry.
//!
//! # Units
//!
//! The operator is homogeneous: with SI input (`eta` in Pa s, geometry in m,
//! `u` in m/s) `E u` is a strain rate in 1/s and `A u` is a viscous force
//! density `div(2 eta eps)` in Pa/m. Fed the nondimensional geometry of
//! [`super::benchmark_thieulot`] it returns that module's nondimensional
//! quantities instead; nothing here carries a hidden SI constant.

use super::MantleError;
use pf_grid::shell_ops::{OpsError, ShellTopology};

/// Number of independent components of a symmetric 3x3 tensor.
pub const VOIGT_LEN: usize = 6;

/// The `(row, col)` of each Voigt slot: `xx, yy, zz, xy, xz, yz`.
pub const VOIGT_INDEX: [(usize, usize); VOIGT_LEN] =
    [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

/// Multiplicities that make [`strain_contraction`] the true Frobenius
/// contraction `eps : eps` of the full 3x3 tensor.
const VOIGT_METRIC: [f64; VOIGT_LEN] = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];

/// Unknowns of the per-cell affine fit: 3 offsets plus 9 gradient entries.
const FIT_DOF: usize = 12;

/// A symmetric strain-rate tensor in Voigt order, 1/s.
pub type StrainVoigt = [f64; VOIGT_LEN];

/// `a : b` for two symmetric tensors in Voigt order.
pub fn strain_contraction(a: &StrainVoigt, b: &StrainVoigt) -> f64 {
    a.iter()
        .zip(b.iter())
        .zip(VOIGT_METRIC.iter())
        .map(|((x, y), m)| m * x * y)
        .sum()
}

/// Expand a Voigt strain to the full symmetric 3x3 tensor.
pub fn strain_tensor(e: &StrainVoigt) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for (slot, (a, d)) in e.iter().zip(VOIGT_INDEX.iter()) {
        t[*a][*d] = *slot;
        t[*d][*a] = *slot;
    }
    t
}

/// `trace(eps)`, the discrete dilation rate of a cell.
pub fn strain_trace(e: &StrainVoigt) -> f64 {
    e[0] + e[1] + e[2]
}

fn shell(what: &'static str, e: OpsError) -> MantleError {
    MantleError::Shell {
        what,
        detail: e.to_string(),
    }
}

fn wrong_length(field: &'static str, got: usize, want: usize) -> MantleError {
    MantleError::Field {
        field,
        value: format!("length {got}, want {want}"),
        why: "does not match the shell this operator was built for",
    }
}

/// Which faces one cell's affine fit is allowed to read.
///
/// The operator uses exactly one of these
/// ([`FitStencil::CellAndFaceNeighbours`]). The other two exist because the
/// conditioning question of
/// `docs/reviews/planet-foundry-b1/MC0_CONDITIONING_GATES.md` is a question
/// about WHICH stencil and not about the fit: on a box, a cell's own six
/// face-normal scalars cannot see a constant shear AT ALL, however exactly
/// they are integrated, and the gate is whether the wider stencil escapes
/// that degeneracy or merely postpones it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitStencil {
    /// All six facets of the cell, boundary facets included. **Not a legal
    /// operator stencil** - a boundary facet carries no unknown - but it is
    /// the local MFD observability map, and it is the geometry the collapse
    /// question is actually about.
    CellFacets,
    /// The cell's own INTERIOR faces: [`FitStencil::CellFacets`] minus the
    /// impermeable ones, i.e. what a purely local operator could legally
    /// read here.
    CellInteriorFaces,
    /// Every interior face of the cell plus every interior face of each of
    /// its face neighbours. **What the operator uses.**
    CellAndFaceNeighbours,
}

/// The three geometry arrays every per-cell fit needs, fetched once per
/// shell instead of once per cell.
#[derive(Debug, Clone)]
pub struct FitGeometry {
    centroids: Vec<[f64; 3]>,
    moments: Vec<[[f64; 3]; 3]>,
    v_poly: Vec<f64>,
}

impl FitGeometry {
    pub fn build(topo: &ShellTopology) -> Result<Self, MantleError> {
        Ok(Self {
            centroids: topo
                .cell_centroids_m()
                .map_err(|e| shell("cell centroids", e))?,
            moments: topo
                .facet_tensor_moments()
                .map_err(|e| shell("facet tensor moments", e))?,
            v_poly: topo
                .polyhedral_cell_volumes_m3()
                .map_err(|e| shell("polyhedral cell volumes", e))?,
        })
    }

    pub fn cell_centroids_m(&self) -> &[[f64; 3]] {
        &self.centroids
    }

    pub fn facet_tensor_moments(&self) -> &[[[f64; 3]; 3]] {
        &self.moments
    }

    pub fn polyhedral_cell_volumes_m3(&self) -> &[f64] {
        &self.v_poly
    }
}

/// One cell's local weighted least-squares system, in the dimensionless
/// coordinates the fit actually solves in.
///
/// The model is `v(x) = b + (G / h_c) (x - x_c)`, so all twelve parameters
/// carry the units of `u_f` and every design entry is `O(1)`: the `b` block
/// is a unit normal and the `G` block is `n_d (mean_f x_e - x_c,e) / h_c`.
/// Nothing here needs a further nondimensionalisation before its singular
/// values mean something.
#[derive(Debug, Clone, Default)]
pub struct CellFit {
    /// Stencil faces, in assembly order.
    pub faces: Vec<u32>,
    /// The exact functional each `u_f` is:
    /// `[n_x, n_y, n_z, (T_f[d][e] - A_f[d] x_c[e]) / (|A_f| h_c)]`.
    pub rows: Vec<[f64; FIT_DOF]>,
    /// `w_f = h_c^2 / |x_f - x_c|^2`, dimensionless.
    pub weights: Vec<f64>,
    /// `h_c = V_poly^(1/3)`.
    pub h: f64,
    pub centroid_m: [f64; 3],
}

impl CellFit {
    /// `sum_f w_f r_f r_f^T`, the normal matrix, accumulated in stencil
    /// order so that it is bit-identical to the one [`StrainOperator::build`]
    /// used to form inline.
    pub fn gram(&self) -> [[f64; FIT_DOF]; FIT_DOF] {
        let mut gram = [[0.0f64; FIT_DOF]; FIT_DOF];
        for (r, w) in self.rows.iter().zip(self.weights.iter()) {
            for (a, grow) in gram.iter_mut().enumerate() {
                for (b, slot) in grow.iter_mut().enumerate() {
                    *slot += w * r[a] * r[b];
                }
            }
        }
        gram
    }

    /// `sum_de M[d][e] r_f[3 + 3d + e]`: the face datum that the affine
    /// field `v(x) = M (x - x_c) / h_c` produces on stencil face `index`.
    /// This is the row of the strain-visibility map, and it is taken from
    /// the SAME design row the operator inverts.
    pub fn gradient_response(&self, index: usize, m: &[[f64; 3]; 3]) -> f64 {
        let r = &self.rows[index];
        let mut acc = 0.0;
        for (d, mrow) in m.iter().enumerate() {
            for (e, mv) in mrow.iter().enumerate() {
                acc += mv * r[3 + 3 * d + e];
            }
        }
        acc
    }
}

/// Assemble one cell's fit system into `out`, reusing its buffers.
///
/// This is the single assembly path: [`StrainOperator::build`] calls it, so a
/// diagnostic that reads a [`CellFit`] is reading the operator's own rows and
/// weights rather than a re-derivation that could drift from them.
pub fn cell_fit_system(
    topo: &ShellTopology,
    geo: &FitGeometry,
    cell: usize,
    stencil: FitStencil,
    out: &mut CellFit,
) {
    out.faces.clear();
    out.rows.clear();
    out.weights.clear();
    // The stencil: every interior face of this cell, plus (for the operator's
    // choice) every interior face of each of its face neighbours. Boundary
    // faces are excluded entirely - that omission IS the natural free-slip
    // condition, and it is also what keeps a rigid rotation exactly null in a
    // boundary cell, since the discrete rotation field is zeroed on boundary
    // faces and would otherwise be read as data that no affine field can
    // match.
    for (f, _) in topo.cell_faces(cell as u32) {
        let face = &topo.faces()[f as usize];
        let Some(nb) = face.neighbor else {
            if stencil == FitStencil::CellFacets {
                push_unique(&mut out.faces, f);
            }
            continue;
        };
        push_unique(&mut out.faces, f);
        if stencil != FitStencil::CellAndFaceNeighbours {
            continue;
        }
        let other = if face.owner as usize == cell {
            nb
        } else {
            face.owner
        };
        for (g, _) in topo.cell_faces(other) {
            if topo.faces()[g as usize].neighbor.is_some() {
                push_unique(&mut out.faces, g);
            }
        }
    }

    let x_c = geo.centroids[cell];
    let h = geo.v_poly[cell].cbrt();
    out.h = h;
    out.centroid_m = x_c;
    for &f in out.faces.iter() {
        let face = &topo.faces()[f as usize];
        let t = &geo.moments[f as usize];
        let inv_a = 1.0 / face.area_m2;
        let mut r = [0.0f64; FIT_DOF];
        for (d, slot) in r[..3].iter_mut().enumerate() {
            *slot = face.normal[d];
        }
        for (d, trow) in t.iter().enumerate() {
            for (e, tv) in trow.iter().enumerate() {
                r[3 + 3 * d + e] = (tv - face.area_vector_m2[d] * x_c[e]) * inv_a / h;
            }
        }
        // Any positive weight leaves the fit exact for affine fields; this
        // one is the dimensionless inverse-square distance the scalar
        // least-squares gradient in `shell_ops` already uses, so a near face
        // is trusted more than a far one.
        let dx = sub(face.centroid_m, x_c);
        let w = h * h / dot(dx, dx);
        out.rows.push(r);
        out.weights.push(w);
    }
}

/// `E`: face-normal velocity unknowns -> one symmetric strain-rate tensor per
/// cell, assembled once from geometry and then applied matrix-free.
///
/// Storage is the sparse stencil itself: for every cell, the faces its strain
/// reads and one symmetric 3x3 coefficient per face. That is
/// `6 * stencil_len` doubles per cell (about 1.7 kB at the measured stencil
/// width of 36), and it exists because the alternative - re-solving a 12x12
/// least-squares problem per cell per operator application - would dominate
/// every Krylov iteration.
#[derive(Debug, Clone)]
pub struct StrainOperator {
    cell_count: usize,
    face_count: usize,
    /// `offsets[c] .. offsets[c + 1]` indexes this cell's stencil entries.
    offsets: Vec<u32>,
    stencil_faces: Vec<u32>,
    /// `d eps_c / d u_f`, Voigt, 1/(s per m/s) = 1/m.
    coeffs: Vec<StrainVoigt>,
    max_stencil: usize,
    worst_gram_ratio: f64,
}

impl StrainOperator {
    /// Assemble `E` from the shell geometry.
    ///
    /// Fails loud if any cell's fit is rank deficient: a strain operator with
    /// a singular local system would silently return a plausible wrong
    /// gradient, which is exactly the failure mode this module exists to
    /// avoid.
    pub fn build(topo: &ShellTopology) -> Result<Self, MantleError> {
        let geo = FitGeometry::build(topo)?;

        let cell_count = topo.cell_count();
        let face_count = topo.face_count();
        let mut offsets = Vec::with_capacity(cell_count + 1);
        let mut stencil_faces: Vec<u32> = Vec::with_capacity(cell_count * 36);
        let mut coeffs: Vec<StrainVoigt> = Vec::with_capacity(cell_count * 36);
        let mut max_stencil = 0usize;
        let mut worst_gram_ratio = 0.0f64;

        let mut fit = CellFit::default();

        for c in 0..cell_count {
            offsets.push(stencil_faces.len() as u32);
            cell_fit_system(topo, &geo, c, FitStencil::CellAndFaceNeighbours, &mut fit);
            let h = fit.h;
            let gram = fit.gram();

            let l = cholesky(&gram).ok_or_else(|| MantleError::Field {
                field: "strain fit normal matrix",
                value: format!("cell {c}, {} stencil faces", fit.faces.len()),
                why: "is not positive definite: the face normals in this stencil do not determine an affine velocity field",
            })?;
            let (dmax, dmin) = l
                .iter()
                .enumerate()
                .fold((0.0f64, f64::MAX), |(x, n), (i, row)| {
                    (x.max(row[i]), n.min(row[i]))
                });
            worst_gram_ratio = worst_gram_ratio.max((dmax / dmin).powi(2));

            for (&f, (r, w)) in fit
                .faces
                .iter()
                .zip(fit.rows.iter().zip(fit.weights.iter()))
            {
                let y = cholesky_solve(&l, r);
                let mut coeff = [0.0f64; VOIGT_LEN];
                for (slot, (a, d)) in coeff.iter_mut().zip(VOIGT_INDEX.iter()) {
                    let g_ad = y[3 + 3 * a + d];
                    let g_da = y[3 + 3 * d + a];
                    *slot = w * 0.5 * (g_ad + g_da) / h;
                }
                stencil_faces.push(f);
                coeffs.push(coeff);
            }
            max_stencil = max_stencil.max(fit.faces.len());
        }
        offsets.push(stencil_faces.len() as u32);

        for (i, c) in coeffs.iter().enumerate() {
            if let Some(v) = c.iter().find(|v| !v.is_finite()) {
                return Err(MantleError::NonFinite {
                    field: "strain operator coefficient",
                    index: i,
                    value: *v,
                });
            }
        }

        Ok(Self {
            cell_count,
            face_count,
            offsets,
            stencil_faces,
            coeffs,
            max_stencil,
            worst_gram_ratio,
        })
    }

    pub fn cell_count(&self) -> usize {
        self.cell_count
    }

    pub fn face_count(&self) -> usize {
        self.face_count
    }

    /// The faces one cell's strain reads, and the coefficient of each.
    pub fn cell_stencil(&self, cell: usize) -> impl Iterator<Item = (u32, &StrainVoigt)> {
        let lo = self.offsets[cell] as usize;
        let hi = self.offsets[cell + 1] as usize;
        self.stencil_faces[lo..hi]
            .iter()
            .copied()
            .zip(self.coeffs[lo..hi].iter())
    }

    /// The widest cell stencil: how many face values one cell's strain rate
    /// touches. Reported because the solver's `A`-block bandwidth follows
    /// from it.
    pub fn max_stencil_faces(&self) -> usize {
        self.max_stencil
    }

    /// Total stored stencil entries, i.e. the number of nonzeros of `E`
    /// divided by six.
    pub fn stencil_entries(&self) -> usize {
        self.stencil_faces.len()
    }

    /// Worst per-cell ratio `(max diag L / min diag L)^2` of the fit's
    /// Cholesky factor: a cheap condition-number proxy for the local normal
    /// equations, reported so a badly-shaped cell cannot hide.
    pub fn worst_gram_ratio(&self) -> f64 {
        self.worst_gram_ratio
    }

    /// `E u`: one symmetric strain-rate tensor per cell, 1/s.
    pub fn apply(&self, u: &[f64]) -> Result<Vec<StrainVoigt>, MantleError> {
        if u.len() != self.face_count {
            return Err(wrong_length("face velocity", u.len(), self.face_count));
        }
        let mut out = vec![[0.0f64; VOIGT_LEN]; self.cell_count];
        for (c, eps) in out.iter_mut().enumerate() {
            for (f, coeff) in self.cell_stencil(c) {
                let uf = u[f as usize];
                for (slot, k) in eps.iter_mut().zip(coeff.iter()) {
                    *slot += k * uf;
                }
            }
        }
        Ok(out)
    }
}

/// `A = -M_u^-1 E^T W E`, the isoviscous viscous operator, applied
/// matrix-free.
///
/// `M_u` is the OVER-RELAXED face mass
/// ([`pf_grid::shell_ops::ShellTopology::face_masses_over_relaxed`]), not the
/// plain `|A_f| d_f`: that is the inner product in which the corrected
/// gradient is the exact negative adjoint of the divergence, so it is the
/// inner product the whole saddle point lives in. Using a different one here
/// would make `A` self-adjoint in a space nothing else uses.
#[derive(Debug, Clone)]
pub struct ViscousOperator {
    strain: StrainOperator,
    /// `W_c = 2 eta V_c`, Pa s m^3.
    cell_weights: Vec<f64>,
    /// The diagonal of `M_u`, m^3.
    face_mass: Vec<f64>,
}

impl ViscousOperator {
    /// Build the isoviscous operator for one constant viscosity, in the units
    /// of the shell it is given.
    pub fn isoviscous(topo: &ShellTopology, eta: f64) -> Result<Self, MantleError> {
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
        let weights: Vec<f64> = volumes.iter().map(|v| 2.0 * eta * v).collect();
        Self::from_cell_weights(topo, weights)
    }

    /// The variable-viscosity seam, present but unused at MC0: the ONLY place
    /// viscosity enters is this per-cell weight `W_c = 2 eta_c V_c`. MC2 fills
    /// this vector differently and changes nothing else.
    pub fn from_cell_weights(
        topo: &ShellTopology,
        cell_weights: Vec<f64>,
    ) -> Result<Self, MantleError> {
        if cell_weights.len() != topo.cell_count() {
            return Err(wrong_length(
                "cell viscous weight",
                cell_weights.len(),
                topo.cell_count(),
            ));
        }
        for (index, w) in cell_weights.iter().enumerate() {
            if !w.is_finite() || *w <= 0.0 {
                return Err(MantleError::NonFinite {
                    field: "cell viscous weight (2 eta V)",
                    index,
                    value: *w,
                });
            }
        }
        let corr = topo
            .face_corrections()
            .map_err(|e| shell("face corrections", e))?;
        let face_mass = topo
            .face_masses_over_relaxed(&corr)
            .map_err(|e| shell("over-relaxed face masses", e))?;
        Ok(Self {
            strain: StrainOperator::build(topo)?,
            cell_weights,
            face_mass,
        })
    }

    pub fn strain(&self) -> &StrainOperator {
        &self.strain
    }

    /// The diagonal of `M_u`, the inner product `A` is self-adjoint in.
    pub fn face_mass(&self) -> &[f64] {
        &self.face_mass
    }

    pub fn cell_weights(&self) -> &[f64] {
        &self.cell_weights
    }

    /// `E^T W E u`, the negative of `M_u A u`. This is the symmetric,
    /// positive semidefinite half of the operator; every structural property
    /// is a statement about this vector, and dividing by `M_u` afterwards
    /// cannot break any of them.
    pub fn stiffness_apply(&self, u: &[f64]) -> Result<Vec<f64>, MantleError> {
        let eps = self.strain.apply(u)?;
        let mut out = vec![0.0f64; self.strain.face_count()];
        for (c, e) in eps.iter().enumerate() {
            let w = self.cell_weights[c];
            let mut scaled = [0.0f64; VOIGT_LEN];
            for (slot, (ev, m)) in scaled.iter_mut().zip(e.iter().zip(VOIGT_METRIC.iter())) {
                *slot = w * m * ev;
            }
            for (f, coeff) in self.strain.cell_stencil(c) {
                out[f as usize] += coeff
                    .iter()
                    .zip(scaled.iter())
                    .map(|(k, s)| k * s)
                    .sum::<f64>();
            }
        }
        Ok(out)
    }

    /// `A u = -M_u^-1 E^T W E u`, the viscous force density `div(2 eta eps)`
    /// per unit face mass. Boundary faces carry no unknown and get exactly
    /// zero.
    pub fn apply(&self, u: &[f64]) -> Result<Vec<f64>, MantleError> {
        let mut out = self.stiffness_apply(u)?;
        for (index, (slot, m)) in out.iter_mut().zip(self.face_mass.iter()).enumerate() {
            *slot = -*slot / m;
            if !slot.is_finite() {
                return Err(MantleError::NonFinite {
                    field: "viscous operator",
                    index,
                    value: *slot,
                });
            }
        }
        Ok(out)
    }

    /// `Phi(u) = sum_c 2 eta V_c eps_c : eps_c >= 0`, the discrete viscous
    /// dissipation. Equal to `u . E^T W E u` by construction, which
    /// `tests::the_energy_and_the_gram_form_are_the_same_number` checks.
    pub fn energy(&self, u: &[f64]) -> Result<f64, MantleError> {
        let eps = self.strain.apply(u)?;
        Ok(eps
            .iter()
            .zip(self.cell_weights.iter())
            .map(|(e, w)| w * strain_contraction(e, e))
            .sum())
    }

    /// `u . M_u A u`, the rate of change of kinetic energy the viscous term
    /// produces. Always `-energy(u) <= 0`; returned as its own function
    /// because a solver should be able to assert the sign cheaply.
    pub fn dissipation(&self, u: &[f64]) -> Result<f64, MantleError> {
        Ok(-self.energy(u)?)
    }

    /// The interior faces, in ascending order: the actual unknown vector of
    /// the Stokes system, and the index map every dense assembly below uses.
    pub fn unknown_faces(topo: &ShellTopology) -> Vec<u32> {
        topo.faces()
            .iter()
            .enumerate()
            .filter(|(_, f)| f.neighbor.is_some())
            .map(|(i, _)| i as u32)
            .collect()
    }

    /// Assemble `E^T W E` densely on the interior faces, for tests and for a
    /// small-shell spectrum. `O(n_faces)` operator applications, so this is a
    /// diagnostic, not a solver path.
    pub fn dense_stiffness(&self, unknowns: &[u32]) -> Result<Vec<Vec<f64>>, MantleError> {
        let mut out = Vec::with_capacity(unknowns.len());
        let mut e = vec![0.0f64; self.strain.face_count()];
        for &f in unknowns {
            e[f as usize] = 1.0;
            let col = self.stiffness_apply(&e)?;
            e[f as usize] = 0.0;
            out.push(unknowns.iter().map(|&g| col[g as usize]).collect());
        }
        Ok(out)
    }
}

/// Append `f` if it is not already present. Linear in a list that never
/// exceeds a few dozen entries, which beats a hash set at this size.
fn push_unique(list: &mut Vec<u32>, f: u32) {
    if !list.contains(&f) {
        list.push(f);
    }
}

/// Lower Cholesky factor of a symmetric matrix, or `None` if it is not
/// positive definite. No pivoting: for an SPD normal matrix none is needed,
/// and a failure here is a rank deficiency the caller must report rather than
/// paper over.
fn cholesky(m: &[[f64; FIT_DOF]; FIT_DOF]) -> Option<[[f64; FIT_DOF]; FIT_DOF]> {
    let mut l = [[0.0f64; FIT_DOF]; FIT_DOF];
    for i in 0..FIT_DOF {
        for j in 0..=i {
            let s = m[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if s <= 0.0 || s.is_nan() {
                    return None;
                }
                l[i][i] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    Some(l)
}

/// Solve `L L^T y = rhs` for a lower Cholesky factor `L`.
fn cholesky_solve(l: &[[f64; FIT_DOF]; FIT_DOF], rhs: &[f64; FIT_DOF]) -> [f64; FIT_DOF] {
    let mut y = [0.0f64; FIT_DOF];
    for i in 0..FIT_DOF {
        y[i] = (rhs[i] - (0..i).map(|k| l[i][k] * y[k]).sum::<f64>()) / l[i][i];
    }
    for i in (0..FIT_DOF).rev() {
        y[i] = (y[i] - ((i + 1)..FIT_DOF).map(|k| l[k][i] * y[k]).sum::<f64>()) / l[i][i];
    }
    y
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_grid::shell::ShellSpec;

    const R_CMB: f64 = 3.480e6;
    const R_TOP: f64 = 6.271e6;
    /// Upper-mantle order of magnitude. The operator is linear in `eta`, so
    /// every ratio reported below is independent of this number; it is here
    /// so the printed magnitudes are recognisable.
    const ETA: f64 = 1.0e21;

    fn mini(n: usize, layers: usize) -> ShellTopology {
        let spec = ShellSpec::new(n, layers, R_CMB, R_TOP).expect("valid shell");
        ShellTopology::build(&spec).expect("topology builds")
    }

    fn op(n: usize, layers: usize) -> (ShellTopology, ViscousOperator) {
        let t = mini(n, layers);
        let a = ViscousOperator::isoviscous(&t, ETA).expect("operator builds");
        (t, a)
    }

    fn max_abs(v: &[f64]) -> f64 {
        v.iter().fold(0.0f64, |a, b| a.max(b.abs()))
    }

    /// The EXACT face data of the affine field `v(x) = b + L x`:
    /// `u_f = (A_f . b + sum_de L[d][e] T_f[d][e]) / |A_f|`, using the exact
    /// facet tensor moments. Boundary faces carry no unknown and are set to
    /// zero, which is what the discrete space actually contains.
    fn affine_face_field(t: &ShellTopology, b: [f64; 3], l: [[f64; 3]; 3]) -> Vec<f64> {
        let moments = t.facet_tensor_moments().expect("facet tensor moments");
        t.faces()
            .iter()
            .zip(moments.iter())
            .map(|(face, tf)| {
                if face.neighbor.is_none() {
                    return 0.0;
                }
                let mut acc = dot(face.area_vector_m2, b);
                for (d, lrow) in l.iter().enumerate() {
                    for (e, lv) in lrow.iter().enumerate() {
                        acc += lv * tf[d][e];
                    }
                }
                acc / face.area_m2
            })
            .collect()
    }

    /// The same field sampled at the facet centroid instead of integrated
    /// over the facet: the control that proves the machine-zero results below
    /// are not a tolerance artefact.
    fn point_sampled_face_field(t: &ShellTopology, b: [f64; 3], l: [[f64; 3]; 3]) -> Vec<f64> {
        t.faces()
            .iter()
            .map(|face| {
                if face.neighbor.is_none() {
                    return 0.0;
                }
                let x = face.centroid_m;
                let mut v = b;
                for (d, slot) in v.iter_mut().enumerate() {
                    *slot += dot(l[d], x);
                }
                dot(v, face.normal)
            })
            .collect()
    }

    /// `L` of the rigid rotation `v = omega x x`.
    fn rotation_l(omega: [f64; 3]) -> [[f64; 3]; 3] {
        [
            [0.0, -omega[2], omega[1]],
            [omega[2], 0.0, -omega[0]],
            [-omega[1], omega[0], 0.0],
        ]
    }

    /// The six rigid-body modes of the discrete space: three translations
    /// (exact face normals) and three rotations (exact facet moments).
    fn rigid_modes(t: &ShellTopology) -> Vec<Vec<f64>> {
        let mut out = Vec::with_capacity(6);
        for axis in 0..3 {
            let mut b = [0.0f64; 3];
            b[axis] = 1.0;
            out.push(affine_face_field(t, b, [[0.0; 3]; 3]));
        }
        for axis in 0..3 {
            let mut omega = [0.0f64; 3];
            omega[axis] = 1.0;
            out.push(t.rigid_rotation_field(omega));
        }
        out
    }

    /// Deterministic pseudo-random numbers in [-1, 1). A named generator, so
    /// "many random fields" is reproducible rather than a mood.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
        }

        fn field(&mut self, t: &ShellTopology) -> Vec<f64> {
            t.faces()
                .iter()
                .map(|f| {
                    if f.neighbor.is_none() {
                        0.0
                    } else {
                        self.next()
                    }
                })
                .collect()
        }
    }

    /// Cells touching a cube seam: any cell on the border of its cube face.
    /// The brief (rev 5, "MC design and qualification", 2026-09-08
    /// correction) requires three series - face interior, seam touching and
    /// global - because on this grid the two populations fail differently.
    fn seam_touching(t: &ShellTopology, c: usize) -> bool {
        let n = t.spec().n();
        let (_, i, j, _) = t.spec().fijk(c as u64).expect("cell index in range");
        i == 0 || j == 0 || i + 1 == n || j + 1 == n
    }

    /// A cube-corner cell: on the border in BOTH lateral directions.
    fn cube_corner(t: &ShellTopology, c: usize) -> bool {
        let n = t.spec().n();
        let (_, i, j, _) = t.spec().fijk(c as u64).expect("cell index in range");
        (i == 0 || i + 1 == n) && (j == 0 || j + 1 == n)
    }

    // ---------------------------------------------------------------
    // 1. Symmetry

    #[test]
    fn the_operator_is_symmetric_in_the_face_mass_inner_product() {
        // A is self-adjoint in the M_u inner product iff M_u A is symmetric.
        // Assembled densely and measured entrywise: this is the property the
        // Gram form buys unconditionally, so a measurable asymmetry would
        // mean the assembly is not the Gram form it claims to be.
        // `ShellSpec` refuses an odd `n`, so the two shells here are the
        // two smallest legal ones; both are entirely seam and corner cells,
        // which is where an asymmetry would be most likely to hide.
        for (n, layers) in [(2usize, 2usize), (4, 2)] {
            let (t, a) = op(n, layers);
            let unknowns = ViscousOperator::unknown_faces(&t);
            let s = a.dense_stiffness(&unknowns).expect("dense stiffness");
            let mass = a.face_mass();

            let scale = s
                .iter()
                .flat_map(|r| r.iter())
                .fold(0.0f64, |acc, v| acc.max(v.abs()));
            let mut worst_stiffness = 0.0f64;
            let mut worst_operator = 0.0f64;
            for (i, row) in s.iter().enumerate() {
                for (j, sij) in row.iter().enumerate() {
                    worst_stiffness = worst_stiffness.max((sij - s[j][i]).abs() / scale);
                    // A_ij = -S_ij / m_i, so (M_u A)_ij = -S_ij: the same
                    // number, formed the long way round so the reported
                    // quantity is literally the packet's "asymmetry of A in
                    // the M_u inner product".
                    let a_ij = -sij / mass[unknowns[i] as usize];
                    let a_ji = -s[j][i] / mass[unknowns[j] as usize];
                    let m_i = mass[unknowns[i] as usize];
                    let m_j = mass[unknowns[j] as usize];
                    worst_operator = worst_operator.max((m_i * a_ij - m_j * a_ji).abs() / scale);
                }
            }
            println!(
                "n={n} layers={layers} unknowns={} worst relative asymmetry: E^T W E {worst_stiffness:.3e}  M_u A {worst_operator:.3e}",
                unknowns.len()
            );
            // 1e-14 is the accumulated-roundoff bar for a sum of this length,
            // not a fitted bound: the measured values are ~1e-16, and the
            // failure this guards against (a non-variational assembly) is
            // O(1), not O(1e-13).
            assert!(
                worst_stiffness < 1e-14,
                "n={n}: E^T W E is not symmetric: {worst_stiffness:e}"
            );
            assert!(
                worst_operator < 1e-14,
                "n={n}: A is not M_u-self-adjoint: {worst_operator:e}"
            );
        }
    }

    // ---------------------------------------------------------------
    // 2. Dissipativity

    #[test]
    fn the_operator_is_dissipative_for_every_field_and_strictly_so_off_the_null_space() {
        let (t, a) = op(4, 3);
        let mut rng = Lcg(0x5eed_1234_9abc_def0);
        let mut worst_positive = 0.0f64;
        let mut weakest_ratio = f64::MAX;
        let modes = rigid_modes(&t);
        for trial in 0..300 {
            let u = rng.field(&t);
            let d = a.dissipation(&u).expect("dissipation");
            // Not a tolerance: Phi is a sum of positive weights times
            // squares, so a positive dissipation is impossible in floating
            // point, not merely unlikely.
            assert!(d <= 0.0, "trial {trial}: dissipation {d:e} is positive");
            worst_positive = worst_positive.max(d.max(0.0));
            // Strictness: a random field is not a rigid mode, and its
            // dissipation must be a real number, not a rounding residue.
            // Compare against the energy a unit-strain field of the same
            // face-norm would carry.
            let norm2: f64 = u.iter().map(|x| x * x).sum();
            let ratio = -d / norm2;
            weakest_ratio = weakest_ratio.min(ratio);
            assert!(
                ratio > 0.0,
                "trial {trial}: a random field is in the null space, ratio {ratio:e}"
            );
        }
        // The same statement for the rigid modes, which sit EXACTLY on the
        // boundary of the inequality. Measured against the weakest random
        // field's Rayleigh quotient, so "zero" is expressed in the operator's
        // own units rather than as an absolute number nobody can calibrate.
        let mut worst_rigid = 0.0f64;
        for (k, m) in modes.iter().enumerate() {
            let d = a.dissipation(m).expect("dissipation");
            let norm2: f64 = m.iter().map(|x| x * x).sum();
            let relative = (-d / norm2) / weakest_ratio;
            worst_rigid = worst_rigid.max(relative);
            println!(
                "rigid mode {k}: dissipation {d:.4e}, Rayleigh quotient {:.4e} = {relative:.3e} of the weakest random field's",
                -d / norm2
            );
        }
        assert!(
            worst_rigid < 1e-20,
            "a rigid mode dissipates {worst_rigid:e} of a random field's Rayleigh quotient; \
             that is not a null vector"
        );
        println!(
            "300 random fields: worst positive dissipation {worst_positive:.3e}, \
             weakest -Phi/|u|^2 {weakest_ratio:.4e}"
        );
        assert_eq!(
            worst_positive, 0.0,
            "dissipation was never allowed to be positive"
        );
    }

    // ---------------------------------------------------------------
    // 3 and 4. The rigid modes

    #[test]
    fn rigid_translation_is_annihilated_exactly() {
        // A uniform translation has zero strain everywhere. This is the test
        // that catches a broken strain assembly fastest: any error in the
        // design rows or in the sign convention shows up here immediately.
        let b = [0.7, -0.2, 0.45];
        for (n, layers) in [(4usize, 3usize), (8, 6)] {
            let (t, a) = op(n, layers);
            let u = affine_face_field(&t, b, [[0.0; 3]; 3]);
            let eps = a.strain().apply(&u).expect("strain");
            let h = t.spec().cell_angle_rad() * R_TOP;
            let strain_scale = max_abs(&b) / h;
            let worst_strain = eps.iter().fold(0.0f64, |acc, e| acc.max(max_abs(e))) / strain_scale;

            let au = a.apply(&u).expect("apply");
            // Reference magnitude: the same operator on a field with a real
            // strain rate of comparable size.
            let l = [[0.0, 0.6, 0.0], [0.6, 0.0, 0.0], [0.0, 0.0, 0.0]];
            let reference = a
                .apply(&affine_face_field(
                    &t,
                    [0.0; 3],
                    scale_matrix(l, max_abs(&b) / h),
                ))
                .expect("apply");
            let ratio = max_abs(&au) / max_abs(&reference);
            println!(
                "n={n:<2} translation: worst |eps| / (|b|/h) = {worst_strain:.4e}, \
                 |A u| / |A u_strain| = {ratio:.4e}"
            );
            assert!(
                worst_strain < 1e-12,
                "n={n}: translation is not strain free: {worst_strain:e}"
            );
            assert!(
                ratio < 1e-12,
                "n={n}: A does not annihilate translation: {ratio:e}"
            );
        }
    }

    #[test]
    fn rigid_rotation_is_annihilated_at_machine_zero_and_not_at_second_order() {
        // The whole point of the exact facet tensor moments. If any stage of
        // E point-samples a facet, this residual becomes O(h^2) - visibly
        // falling by 4x per refinement and starting near 1e-3 - instead of
        // roundoff. The point-sampled control in the same loop shows exactly
        // that failure, so the tolerance below is separated from the failure
        // mode by nine orders of magnitude rather than by a judgement call.
        let omega = [0.7, -1.3, 0.4];
        let l = rotation_l(omega);
        let scale = max_abs(&omega);
        let mut exact_series = Vec::new();
        let mut sampled_series = Vec::new();
        for (n, layers) in [(4usize, 3usize), (8, 6), (16, 10)] {
            let (t, a) = op(n, layers);
            let u = t.rigid_rotation_field(omega);
            // The exact-moment face field and shell_ops' own rotation field
            // must be the same vector; if they ever diverge, the two
            // constructions have drifted apart and every result below is
            // about a different field than the solver will see.
            let from_moments = affine_face_field(&t, [0.0; 3], l);
            let drift = u
                .iter()
                .zip(from_moments.iter())
                .fold(0.0f64, |acc, (x, y)| acc.max((x - y).abs()))
                / max_abs(&u);
            assert!(
                drift < 1e-14,
                "n={n}: rigid_rotation_field and the affine construction differ by {drift:e}"
            );

            let eps = a.strain().apply(&u).expect("strain");
            let exact = eps.iter().fold(0.0f64, |acc, e| acc.max(max_abs(e))) / scale;

            let sampled_u = point_sampled_face_field(&t, [0.0; 3], l);
            let sampled = a
                .strain()
                .apply(&sampled_u)
                .expect("strain")
                .iter()
                .fold(0.0f64, |acc, e| acc.max(max_abs(e)))
                / scale;

            let au = a.apply(&u).expect("apply");
            let reference = a
                .apply(&affine_face_field(
                    &t,
                    [0.0; 3],
                    [[0.0, scale, 0.0], [scale, 0.0, 0.0], [0.0, 0.0, 0.0]],
                ))
                .expect("apply");
            let ratio = max_abs(&au) / max_abs(&reference);
            println!(
                "n={n:<2} rotation: worst |eps|/|omega| exact moments {exact:.4e}  \
                 point sample {sampled:.4e}   |A u|/|A u_strain| {ratio:.4e}"
            );
            assert!(
                exact < 1e-12,
                "n={n}: rotation leaks strain {exact:e}; the strain assembly is point sampling somewhere"
            );
            assert!(
                ratio < 1e-11,
                "n={n}: A does not annihilate rigid rotation: {ratio:e}"
            );
            exact_series.push(exact);
            sampled_series.push(sampled);
        }
        // The control must be visibly worse AND second order, which is what
        // "not O(h^2)" means for the exact series next to it.
        for (k, s) in sampled_series.iter().enumerate() {
            assert!(
                *s > 1e-5,
                "the point-sampled control at index {k} did not leak: {s:e}"
            );
        }
        // The control converges at FIRST order, not second, and the
        // difference is informative rather than incidental: `shell_ops`
        // measures the same substitution at second order because it feeds the
        // point sample straight into a Green-Gauss sum, whereas `E`
        // DIFFERENTIATES its face data, so an O(h^2) error in `u_f` becomes an
        // O(h) error in the strain. Either way it is between nine and eleven
        // orders of magnitude above the exact series next to it.
        for w in sampled_series.windows(2) {
            let order = (w[0] / w[1]).log2();
            println!(
                "point-sampled control order {order:+.2}: {:e} -> {:e}",
                w[0], w[1]
            );
            assert!(
                order > 0.8,
                "the point-sampled control should still converge: {:e} -> {:e} is {order:+.2}",
                w[0],
                w[1]
            );
        }
        // The exact series is REPORTED here and certified elsewhere.
        //
        // This test used to argue that the series is roundoff because it
        // RISES with resolution, "which a discretization error cannot do".
        // That inference is false, and this build contains the
        // counterexample: the inconsistent two-point gradient's error also
        // grows under refinement
        // (`shell_ops::tests::raw_two_point_dg_laplacian_diverges_under_refinement`)
        // and it is pure discretization error. A rising series is evidence of
        // nothing on its own. The certification that actually implies
        // roundoff - the residual falls like the machine epsilon at FIXED h,
        // slope +1.003 over a 29-bit range - is
        // `gate4_the_rotation_residual_scales_with_the_precision_and_not_with_h`.
        let apparent = (exact_series[0] / exact_series[2]).log2() / 2.0;
        println!(
            "exact-moment rotation residual across the series: {:.3e} -> {:.3e} -> {:.3e} \
             (apparent order {apparent:+.2}; NOT a certification either way - see gate 4)",
            exact_series[0], exact_series[1], exact_series[2]
        );
        // What this test can still certify on its own: the separation from
        // the point-sampled control does not close under refinement. That is
        // a comparison between two measured series rather than an inference
        // from the shape of one.
        for (k, (exact, sampled)) in exact_series.iter().zip(sampled_series.iter()).enumerate() {
            let separation = sampled / exact;
            println!("  rung {k}: point-sampled / exact = {separation:.3e}");
            assert!(
                separation > 1e10,
                "rung {k}: the exact path is only {separation:e} times better than point \
                 sampling; the two have converged towards each other and the exact-moment \
                 construction has stopped buying anything"
            );
        }
    }

    fn scale_matrix(l: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
        let mut out = l;
        for row in out.iter_mut() {
            for v in row.iter_mut() {
                *v *= s;
            }
        }
        out
    }

    // ---------------------------------------------------------------
    // 5. The null space, exactly

    /// Cyclic Jacobi eigendecomposition of a dense symmetric matrix, stored
    /// row-major. Returns `(eigenvalues, eigenvectors as columns)`. Written
    /// here rather than pulled in as a dependency: it is thirty lines, it is
    /// used by one test, and the workspace has no linear-algebra crate.
    fn jacobi_eigen(mut a: Vec<f64>, n: usize) -> (Vec<f64>, Vec<f64>) {
        let mut v = vec![0.0f64; n * n];
        for i in 0..n {
            v[i * n + i] = 1.0;
        }
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
                    // Skip a rotation that cannot change the answer: once an
                    // off-diagonal entry is below the roundoff of its own
                    // diagonal pair, rotating on it is pure cost. Without
                    // this the later sweeps dominate the run time.
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
                    for k in 0..n {
                        let (kp, kq) = (v[k * n + p], v[k * n + q]);
                        v[k * n + p] = c * kp - s * kq;
                        v[k * n + q] = s * kp + c * kq;
                    }
                }
            }
        }
        ((0..n).map(|i| a[i * n + i]).collect(), v)
    }

    #[test]
    fn the_null_space_is_exactly_the_six_rigid_body_modes() {
        // null(A) = null(E) is structural. What is NOT structural, and is
        // therefore measured here, is whether E has a null vector beyond the
        // rigid modes - a least-squares fit filters its stencil's residual
        // space, and a spurious near-null mode there would become an
        // hourglass mode the solver's preconditioner has to fight.
        //
        // The similarity transform B = M_u^-1/2 (E^T W E) M_u^-1/2 is
        // symmetric with the same eigenvalues as -A, so the spectrum is
        // computable directly.
        for (n, layers) in [(2usize, 2usize)] {
            let (t, a) = op(n, layers);
            let unknowns = ViscousOperator::unknown_faces(&t);
            let nu = unknowns.len();
            let s = a.dense_stiffness(&unknowns).expect("dense stiffness");
            let inv_sqrt_m: Vec<f64> = unknowns
                .iter()
                .map(|&f| 1.0 / a.face_mass()[f as usize].sqrt())
                .collect();
            let mut b = vec![0.0f64; nu * nu];
            for (i, row) in s.iter().enumerate() {
                for (j, v) in row.iter().enumerate() {
                    b[i * nu + j] = 0.5 * (v + s[j][i]) * inv_sqrt_m[i] * inv_sqrt_m[j];
                }
            }
            let (eig, vecs) = jacobi_eigen(b, nu);
            let lambda_max = eig.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
            let mut order: Vec<usize> = (0..nu).collect();
            order.sort_by(|x, y| eig[*x].partial_cmp(&eig[*y]).expect("finite eigenvalues"));
            let sorted: Vec<f64> = order.iter().map(|k| eig[*k] / lambda_max).collect();
            println!(
                "n={n} layers={layers} unknowns={nu}: smallest normalised eigenvalues {:.3e} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e} | 7th {:.3e} | largest {:.3e}",
                sorted[0], sorted[1], sorted[2], sorted[3], sorted[4], sorted[5], sorted[6], sorted[nu - 1]
            );
            // Zero eigenvalues of a positive semidefinite matrix come out at
            // roundoff times the largest one; a genuine mode is separated
            // from that by orders of magnitude, which is why the count is
            // taken with a gap test rather than at a fixed cut.
            let null_count = sorted.iter().filter(|v| v.abs() < 1e-10).count();
            assert_eq!(
                null_count, 6,
                "n={n} layers={layers}: expected exactly the 6 rigid modes, found {null_count}; \
                 7th eigenvalue is {:.3e} of the largest",
                sorted[6]
            );
            assert!(
                sorted[6] > 1e-6,
                "n={n} layers={layers}: the 7th eigenvalue {:.3e} is too close to the null space",
                sorted[6]
            );

            // The six null directions are the rigid modes and nothing else:
            // every rigid mode, mapped into the same scaled space, must lie
            // inside the computed null eigenspace.
            let modes = rigid_modes(&t);
            let mut worst_leak = 0.0f64;
            for (k, m) in modes.iter().enumerate() {
                let y: Vec<f64> = unknowns
                    .iter()
                    .zip(inv_sqrt_m.iter())
                    .map(|(&f, im)| m[f as usize] / im)
                    .collect();
                let norm = y.iter().map(|x| x * x).sum::<f64>().sqrt();
                let mut residual = y.clone();
                for col in order.iter().take(6) {
                    let proj: f64 = (0..nu).map(|i| vecs[i * nu + col] * y[i]).sum();
                    for (i, r) in residual.iter_mut().enumerate() {
                        *r -= proj * vecs[i * nu + col];
                    }
                }
                let leak = residual.iter().map(|x| x * x).sum::<f64>().sqrt() / norm;
                worst_leak = worst_leak.max(leak);
                assert!(
                    leak < 1e-6,
                    "n={n}: rigid mode {k} is {leak:e} outside the computed null space"
                );
            }
            println!("n={n} layers={layers}: worst rigid-mode leak out of the null eigenspace {worst_leak:.3e}");
        }
    }

    /// In-place Cholesky of a dense symmetric matrix. Returns `false` at the
    /// first non-positive pivot, which by Sylvester's law of inertia is a
    /// certificate that the matrix has an eigenvalue at or below zero.
    fn dense_cholesky_succeeds(m: &mut [f64], n: usize) -> bool {
        for i in 0..n {
            for j in 0..=i {
                let s = m[i * n + j] - (0..j).map(|k| m[i * n + k] * m[j * n + k]).sum::<f64>();
                if i == j {
                    if s <= 0.0 || s.is_nan() {
                        return false;
                    }
                    m[i * n + i] = s.sqrt();
                } else {
                    m[i * n + j] = s / m[j * n + j];
                }
            }
        }
        true
    }

    /// Largest eigenvalue of a dense symmetric matrix by power iteration,
    /// from a deterministic start.
    fn dense_lambda_max(m: &[f64], n: usize) -> f64 {
        let mut v: Vec<f64> = (0..n).map(|i| 1.0 + (i % 7) as f64).collect();
        let mut lambda = 0.0f64;
        for _ in 0..300 {
            let w: Vec<f64> = (0..n)
                .map(|i| (0..n).map(|j| m[i * n + j] * v[j]).sum::<f64>())
                .collect();
            lambda = w.iter().map(|x| x * x).sum::<f64>().sqrt();
            if lambda == 0.0 {
                return 0.0;
            }
            v = w.iter().map(|x| x / lambda).collect();
        }
        lambda
    }

    #[test]
    fn no_spurious_null_mode_survives_deflation_of_the_rigid_body_modes() {
        // The same claim as the spectrum test above, certified on shells too
        // large to diagonalize in a unit test, and certified rather than
        // sampled. Deflate the six rigid modes by adding `lambda_max Q Q^T`,
        // subtract `tau I`, and attempt a Cholesky: it succeeds if and only
        // if every eigenvalue of `B` outside the rigid-body span exceeds
        // `tau` (Sylvester). A spurious near-null mode - the hourglass mode a
        // least-squares stencil could in principle admit - would make the
        // factorization fail, at any `tau` above it.
        for (n, layers) in [(2usize, 3usize), (4, 2), (4, 3)] {
            let (t, a) = op(n, layers);
            let unknowns = ViscousOperator::unknown_faces(&t);
            let nu = unknowns.len();
            let s = a.dense_stiffness(&unknowns).expect("dense stiffness");
            let inv_sqrt_m: Vec<f64> = unknowns
                .iter()
                .map(|&f| 1.0 / a.face_mass()[f as usize].sqrt())
                .collect();
            let mut b = vec![0.0f64; nu * nu];
            for (i, row) in s.iter().enumerate() {
                for (j, v) in row.iter().enumerate() {
                    b[i * nu + j] = 0.5 * (v + s[j][i]) * inv_sqrt_m[i] * inv_sqrt_m[j];
                }
            }
            let lambda_max = dense_lambda_max(&b, nu);

            // The rigid modes in the same scaled space, orthonormalized.
            let mut q: Vec<Vec<f64>> = Vec::new();
            for m in rigid_modes(&t) {
                let mut y: Vec<f64> = unknowns
                    .iter()
                    .zip(inv_sqrt_m.iter())
                    .map(|(&f, im)| m[f as usize] / im)
                    .collect();
                for prev in q.iter() {
                    let ip: f64 = y.iter().zip(prev.iter()).map(|(x, p)| x * p).sum();
                    for (yi, pi) in y.iter_mut().zip(prev.iter()) {
                        *yi -= ip * pi;
                    }
                }
                let norm = y.iter().map(|x| x * x).sum::<f64>().sqrt();
                assert!(
                    norm > 1e-6,
                    "n={n}: the six rigid modes are not independent on this shell"
                );
                for yi in y.iter_mut() {
                    *yi /= norm;
                }
                // Each one really is null, measured before it is deflated.
                let residual: f64 = (0..nu)
                    .map(|i| {
                        let bi: f64 = (0..nu).map(|j| b[i * nu + j] * y[j]).sum();
                        bi * bi
                    })
                    .sum::<f64>()
                    .sqrt()
                    / lambda_max;
                assert!(
                    residual < 1e-12,
                    "n={n}: a rigid mode is not in the null space, |B y| / lambda_max = {residual:e}"
                );
                q.push(y);
            }

            let tau = 1e-6 * lambda_max;
            for (i, row) in b.chunks_mut(nu).enumerate() {
                for (j, slot) in row.iter_mut().enumerate() {
                    let deflate: f64 = q.iter().map(|c| c[i] * c[j]).sum();
                    *slot += lambda_max * deflate;
                    if i == j {
                        *slot -= tau;
                    }
                }
            }
            let ok = dense_cholesky_succeeds(&mut b, nu);
            println!(
                "n={n} layers={layers} unknowns={nu}: lambda_max {lambda_max:.4e}; \
                 B - {tau:.3e} I is positive definite after deflating the 6 rigid modes: {ok}"
            );
            assert!(
                ok,
                "n={n} layers={layers}: a seventh null (or near-null) mode exists below \
                 1e-6 of lambda_max"
            );
        }
    }

    // ---------------------------------------------------------------
    // 6. Constant strain, in three series

    #[test]
    fn constant_strain_is_reproduced_in_every_cell_three_series() {
        // The consistency half of the operator: fed the exact face data of an
        // affine field, E must return that field's symmetric gradient in
        // EVERY cell. Reported as the brief's three series (face interior,
        // seam touching, global) plus the eight cube corners on their own,
        // because a cubed-sphere failure that averages away is not a pass.
        //
        // The affine field carries a rotation and a dilation as well as a
        // shear, so this one measurement also pins that E separates them.
        let b = [3.1e-9, -0.7e-9, 2.2e-9];
        let l = [
            [0.4e-15, -1.1e-15, 0.3e-15],
            [0.9e-15, 0.25e-15, -0.6e-15],
            [-0.2e-15, 0.8e-15, 1.4e-15],
        ];
        let mut expected = [0.0f64; VOIGT_LEN];
        for (slot, (a, d)) in expected.iter_mut().zip(VOIGT_INDEX.iter()) {
            *slot = 0.5 * (l[*a][*d] + l[*d][*a]);
        }
        let scale = expected.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));

        for (n, layers) in [(4usize, 3usize), (8, 6), (16, 10)] {
            let t = mini(n, layers);
            let e = StrainOperator::build(&t).expect("strain operator");
            let u = affine_face_field(&t, b, l);
            let eps = e.apply(&u).expect("strain");
            let (mut interior, mut seam, mut corner, mut global) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            let mut l2 = 0.0f64;
            for (c, got) in eps.iter().enumerate() {
                let err = got
                    .iter()
                    .zip(expected.iter())
                    .fold(0.0f64, |acc, (g, w)| acc.max((g - w).abs()))
                    / scale;
                global = global.max(err);
                l2 += err * err;
                if seam_touching(&t, c) {
                    seam = seam.max(err);
                } else {
                    interior = interior.max(err);
                }
                if cube_corner(&t, c) {
                    corner = corner.max(err);
                }
            }
            l2 = (l2 / eps.len() as f64).sqrt();
            println!(
                "n={n:<2} constant strain reproduction, relative Linf: interior {interior:.3e}  \
                 seam-touching {seam:.3e}  cube-corner {corner:.3e}  global {global:.3e}  L2 {l2:.3e}"
            );
            // Exactness is structural here (a zero-residual least squares
            // returns the model exactly), so the only thing this bound is
            // allowed to admit is roundoff, including the cancellation in
            // T_f - A_f (x) x_c which costs about a digit and a half at the
            // largest shells.
            assert!(
                global < 1e-10,
                "n={n}: constant strain is not reproduced: worst relative error {global:e}"
            );
        }
    }

    // ---------------------------------------------------------------
    // Free slip, the natural condition

    #[test]
    fn boundary_face_values_never_enter_the_operator() {
        // What "natural" means, concretely: no boundary term exists, so no
        // ghost value can be read, so there is no boundary convention to get
        // wrong. Writing garbage into the impermeable entries changes
        // nothing.
        let (t, a) = op(4, 3);
        let mut rng = Lcg(0x00c0_ffee_1234_5678);
        let u = rng.field(&t);
        let mut poisoned = u.clone();
        let mut touched = 0usize;
        for (idx, face) in t.faces().iter().enumerate() {
            if face.neighbor.is_none() {
                poisoned[idx] = 1.0e6 * rng.next();
                touched += 1;
            }
        }
        let a0 = a.apply(&u).expect("apply");
        let a1 = a.apply(&poisoned).expect("apply");
        let diff = a0
            .iter()
            .zip(a1.iter())
            .fold(0.0f64, |acc, (x, y)| acc.max((x - y).abs()));
        println!(
            "{touched} boundary faces poisoned with 1e6-scale values; |A u| changed by {diff:e}"
        );
        assert_eq!(diff, 0.0, "boundary values entered the operator");
        // And the tangential slip they permit is NOT free: only rigid modes
        // are null, so a field concentrated at the wall must dissipate.
        let mut wall = vec![0.0f64; t.face_count()];
        for (idx, face) in t.faces().iter().enumerate() {
            let touches = t
                .cell_faces(face.owner)
                .iter()
                .any(|(g, _)| t.faces()[*g as usize].neighbor.is_none());
            if face.neighbor.is_some() && touches {
                wall[idx] = 1.0;
            }
        }
        let d = a.dissipation(&wall).expect("dissipation");
        assert!(
            d < 0.0,
            "a boundary-layer field must cost energy under free slip, got {d:e}"
        );
        println!("boundary-layer slip field dissipation {d:.4e} (strictly negative: slip is not free of strain)");
    }

    // ---------------------------------------------------------------
    // The two "incompressibility" numbers are not the same number

    #[test]
    fn the_strain_trace_and_the_discrete_divergence_differ_by_the_volume_ratio() {
        // For an affine field the exact facet fluxes sum to V_poly tr(L),
        // while `divergence` divides that sum by the spherical frustum
        // volume. So tr(eps) = (V_sphere / V_poly) (D u) exactly, and the two
        // differ by 1 + O(h^2). The solver has to know this: the constraint
        // it enforces is D u = 0, which is NOT tr(eps) = 0 entrywise.
        let l = [
            [1.3e-15, -0.4e-15, 0.2e-15],
            [0.1e-15, 0.7e-15, -0.9e-15],
            [0.5e-15, 0.3e-15, -0.6e-15],
        ];
        let trace = l[0][0] + l[1][1] + l[2][2];
        for (n, layers) in [(4usize, 3usize), (8, 6), (16, 10)] {
            let t = mini(n, layers);
            let e = StrainOperator::build(&t).expect("strain operator");
            let u = affine_face_field(&t, [0.0; 3], l);
            let eps = e.apply(&u).expect("strain");
            let v_sphere = t.cell_volumes_m3().expect("volumes");
            let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
            let div = t.divergence(&u, &v_sphere).expect("divergence");
            let mut worst_identity = 0.0f64;
            let mut worst_gap = 0.0f64;
            let mut worst_boundary = 0.0f64;
            for (c, e_c) in eps.iter().enumerate() {
                let touches_boundary = t
                    .cell_faces(c as u32)
                    .iter()
                    .any(|(f, _)| t.faces()[*f as usize].neighbor.is_none());
                let gap = (strain_trace(e_c) - div[c]).abs() / trace.abs();
                if touches_boundary {
                    // The impermeable condition deletes a real flux there, so
                    // neither the identity nor the small gap is expected: the
                    // strain still sees the affine field, the divergence does
                    // not. Reported, never averaged in.
                    worst_boundary = worst_boundary.max(gap);
                    continue;
                }
                worst_gap = worst_gap.max(gap);
                let predicted = v_sphere[c] / v_poly[c] * div[c];
                worst_identity =
                    worst_identity.max((strain_trace(e_c) - predicted).abs() / trace.abs());
            }
            println!(
                "n={n:<2} tr(eps) = (V_sphere/V_poly) (D u): residual {worst_identity:.3e}; \
                 interior |tr(eps) - D u| / |tr L| = {worst_gap:.3e}; \
                 boundary-cell gap {worst_boundary:.3e}"
            );
            assert!(
                worst_identity < 1e-9,
                "n={n}: the trace identity does not hold: {worst_identity:e}"
            );
        }
    }

    // ---------------------------------------------------------------
    // Bookkeeping the solver packet will want

    #[test]
    fn the_energy_and_the_gram_form_are_the_same_number() {
        let (t, a) = op(4, 3);
        let mut rng = Lcg(0xfeed_face_0bad_beef);
        for trial in 0..20 {
            let u = rng.field(&t);
            let s = a.stiffness_apply(&u).expect("stiffness");
            let quadratic: f64 = u.iter().zip(s.iter()).map(|(x, y)| x * y).sum();
            let energy = a.energy(&u).expect("energy");
            let rel = (quadratic - energy).abs() / energy;
            assert!(
                rel < 1e-12,
                "trial {trial}: u.(E^T W E)u = {quadratic:e} but Phi(u) = {energy:e} (rel {rel:e})"
            );
            let d = a.dissipation(&u).expect("dissipation");
            assert_eq!(d, -energy, "dissipation is not -Phi");
        }
    }

    #[test]
    fn the_stencil_width_is_measured_not_assumed() {
        // MC0's solver needs the bandwidth. Reported for a real shell, and
        // separately for the A-block row, which is one cell wider than the
        // strain stencil because A = E^T (...) E composes two of them.
        let (_t, a) = op(4, 3);
        let e = a.strain();
        let mean = e.stencil_entries() as f64 / e.cell_count() as f64;
        println!(
            "strain stencil: max {} faces per cell, mean {mean:.2}, {} entries for {} cells / {} faces",
            e.max_stencil_faces(),
            e.stencil_entries(),
            e.cell_count(),
            e.face_count()
        );
        println!(
            "worst local normal-matrix condition proxy (max/min Cholesky pivot)^2: {:.3e}",
            e.worst_gram_ratio()
        );

        let (ts, small) = op(4, 3);
        let unknowns = ViscousOperator::unknown_faces(&ts);
        let s = small.dense_stiffness(&unknowns).expect("dense stiffness");
        let scale = s
            .iter()
            .flat_map(|r| r.iter())
            .fold(0.0f64, |acc, v| acc.max(v.abs()));
        let mut max_row = 0usize;
        let mut total = 0usize;
        for row in &s {
            let nz = row.iter().filter(|v| v.abs() > 1e-14 * scale).count();
            max_row = max_row.max(nz);
            total += nz;
        }
        println!(
            "A-block row (n=4, layers=3, {} unknowns): max {max_row} nonzeros, mean {:.1}",
            unknowns.len(),
            total as f64 / unknowns.len() as f64
        );
        assert!(
            e.max_stencil_faces() <= 36,
            "the strain stencil grew past cell + face neighbours: {}",
            e.max_stencil_faces()
        );
    }

    // ---------------------------------------------------------------
    // Task D: the nullspace projection and the operator must share an
    // inner product

    #[test]
    fn the_rigid_rotation_projection_is_orthogonal_in_the_operators_inner_product() {
        // `rigid_rotation_component` fits the spin in a WEIGHTED least
        // squares. If those weights are not the M_u this operator is
        // self-adjoint in, the projection is orthogonal in a different inner
        // product than the operator it protects, and the leftover component
        // is a small persistent spin with no visible cause. This test is the
        // reason `rigid_rotation_component` was moved onto
        // `face_masses_over_relaxed`.
        let (t, a) = op(4, 3);
        let mass = a.face_mass();
        let mut rng = Lcg(0x1234_5678_9abc_def0);
        let modes = rigid_modes(&t);
        let mut worst = 0.0f64;
        let mut worst_at = (0usize, 0usize);
        for trial in 0..20 {
            let mut u = rng.field(&t);
            t.remove_rigid_rotation(&mut u).expect("removes");
            let self_norm: f64 = u
                .iter()
                .zip(mass.iter())
                .map(|(x, m)| m * x * x)
                .sum::<f64>()
                .sqrt();
            for (k, m) in modes.iter().skip(3).enumerate() {
                let ip: f64 = u
                    .iter()
                    .zip(m.iter().zip(mass.iter()))
                    .map(|(x, (y, w))| w * x * y)
                    .sum();
                let mode_norm: f64 = m
                    .iter()
                    .zip(mass.iter())
                    .map(|(y, w)| w * y * y)
                    .sum::<f64>()
                    .sqrt();
                let rel = ip.abs() / (self_norm * mode_norm);
                if rel > worst {
                    worst = rel;
                    worst_at = (trial, k);
                }
            }
        }
        println!(
            "worst residual M_u-overlap with a rotation after remove_rigid_rotation: {worst:.4e} \
             (trial {}, axis {})",
            worst_at.0, worst_at.1
        );

        // The two measurements the existing shell_ops rotation tests make,
        // printed here so the effect of changing the fit weights is visible
        // as a number rather than as a pass/fail.
        let t4 = mini(4, 3);
        let omega = [1.3e-15, -4.1e-16, 2.7e-15];
        let field = t4.rigid_rotation_field(omega);
        let found = t4.rigid_rotation_component(&field).expect("fit");
        let recovery = (0..3)
            .map(|i| (found[i] - omega[i]).abs() / omega[i].abs())
            .fold(0.0f64, f64::max);
        let t2 = mini(4, 2);
        let radial: Vec<f64> = t2
            .faces()
            .iter()
            .map(|f| {
                if f.neighbor.is_none() {
                    0.0
                } else {
                    let r = dot(f.centroid_m, f.centroid_m).sqrt();
                    dot(
                        [
                            f.centroid_m[0] / r,
                            f.centroid_m[1] / r,
                            f.centroid_m[2] / r,
                        ],
                        f.normal,
                    )
                }
            })
            .collect();
        let spin = t2.rigid_rotation_component(&radial).expect("fit");
        let induced = t2.rigid_rotation_field(spin);
        let reference = t2.rigid_rotation_field([1.0, 1.0, 1.0]);
        println!(
            "shell_ops rotation metrics with the CURRENT weights: recovery relative error {recovery:.6e}, \
             radial-flow spurious spin {:.6e}",
            max_abs(&induced) / max_abs(&reference)
        );

        // Asserted last, so that a failure still prints every number above:
        // this is the pin for the `rigid_rotation_component` weight change,
        // and the before/after values are the evidence for it.
        assert!(
            worst < 1e-12,
            "the projected field still has M_u-overlap {worst:e} with a rigid rotation \
             (trial {}, axis {}); the rotation fit and this operator are using different weights",
            worst_at.0,
            worst_at.1
        );
    }

    // ---------------------------------------------------------------
    // Gate 3: is a real strain OBSERVABLE in face-normal data, and does it
    // stay observable as the cells become box-like?

    const S2: f64 = std::f64::consts::FRAC_1_SQRT_2;

    /// Frobenius-orthonormal basis of the symmetric 3x3 tensors: the six
    /// constant strain modes, each with `|S|_F = 1`.
    fn strain_basis() -> [[[f64; 3]; 3]; 6] {
        [
            [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, S2, 0.0], [S2, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, 0.0, S2], [0.0, 0.0, 0.0], [S2, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, S2], [0.0, S2, 0.0]],
        ]
    }

    /// The matching orthonormal basis of the antisymmetric tensors: the
    /// three rigid spin modes.
    fn spin_basis() -> [[[f64; 3]; 3]; 3] {
        [
            [[0.0, S2, 0.0], [-S2, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, 0.0, S2], [0.0, 0.0, 0.0], [-S2, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, S2], [0.0, -S2, 0.0]],
        ]
    }

    /// Frobenius inner product.
    fn frob(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> f64 {
        let mut s = 0.0;
        for (ra, rb) in a.iter().zip(b.iter()) {
            for (x, y) in ra.iter().zip(rb.iter()) {
                s += x * y;
            }
        }
        s
    }

    /// The normalised local observability map of one cell: the twelve affine
    /// modes, in the order `[6 strain | 3 spin | 3 translation]`, against the
    /// weighted face data each one produces. Row-major, `rows x 12`.
    ///
    /// Row `f` carries the factor `sqrt(w_f / sum w)`, so a column's
    /// Euclidean norm is the weighted RMS face datum of that mode per unit
    /// velocity amplitude - a pure number, independent of `h`, of the
    /// planet's radius, of the units, and of how many faces the stencil
    /// happens to have. A mode's amplitude is `h_c |eps|_F` for a strain and
    /// `|b|` for a translation, which are the same physical size, so the
    /// columns are directly comparable. THAT is the nondimensionalisation
    /// this gate needs: what is left is observability and nothing else.
    fn observability_map(fit: &CellFit) -> Vec<f64> {
        let total: f64 = fit.weights.iter().sum();
        let sym = strain_basis();
        let skew = spin_basis();
        let mut m = vec![0.0f64; fit.rows.len() * 12];
        for (f, w) in fit.weights.iter().enumerate() {
            let s = (w / total).sqrt();
            for (k, mode) in sym.iter().chain(skew.iter()).enumerate() {
                m[f * 12 + k] = s * fit.gradient_response(f, mode);
            }
            for d in 0..3 {
                m[f * 12 + 9 + d] = s * fit.rows[f][d];
            }
        }
        m
    }

    /// One-sided Jacobi SVD of a tall `rows x cols` matrix, row-major:
    /// rotate pairs of columns until they are mutually orthogonal, after
    /// which the column norms ARE the singular values. Returns them
    /// ascending and leaves `m` rotated.
    ///
    /// Chosen over an eigendecomposition of `M^T M` deliberately. Forming
    /// the Gram squares the condition number, and the first attempt at this
    /// gate did exactly that and returned `sigma_min` as an exact `0` for
    /// every cell class - the small singular value had been annihilated by
    /// the cancellation in `N_ss - N_sr N_rr^-1 N_rs`, not by the geometry.
    /// One-sided Jacobi computes small singular values to high RELATIVE
    /// accuracy, and the smallest one is the entire question here.
    fn one_sided_jacobi(m: &mut [f64], rows: usize, cols: usize) -> Vec<f64> {
        for _sweep in 0..60 {
            let mut rotated = false;
            for p in 0..cols {
                for q in (p + 1)..cols {
                    let (mut app, mut aqq, mut apq) = (0.0f64, 0.0f64, 0.0f64);
                    for i in 0..rows {
                        let (x, y) = (m[i * cols + p], m[i * cols + q]);
                        app += x * x;
                        aqq += y * y;
                        apq += x * y;
                    }
                    if apq.abs() <= 1e-17 * (app * aqq).sqrt() {
                        continue;
                    }
                    rotated = true;
                    let theta = (aqq - app) / (2.0 * apq);
                    let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                    let c = 1.0 / (t * t + 1.0).sqrt();
                    let s = t * c;
                    for i in 0..rows {
                        let (x, y) = (m[i * cols + p], m[i * cols + q]);
                        m[i * cols + p] = c * x - s * y;
                        m[i * cols + q] = s * x + c * y;
                    }
                }
            }
            if !rotated {
                break;
            }
        }
        let mut sv: Vec<f64> = (0..cols)
            .map(|p| {
                (0..rows)
                    .map(|i| m[i * cols + p] * m[i * cols + p])
                    .sum::<f64>()
                    .sqrt()
            })
            .collect();
        sv.sort_by(|a, b| a.partial_cmp(b).expect("finite singular values"));
        sv
    }

    /// Split the observability map into its strain columns and the same
    /// columns with the rigid ones projected out.
    ///
    /// The projection is done in the DATA space by modified Gram-Schmidt
    /// with a second pass, never through a normal matrix: `x^T S x` for the
    /// Schur complement `S` is the square of a quantity we need three or
    /// four digits of, so it may not be formed by subtracting two nearly
    /// equal squares. A rigid direction that produces no face data at all -
    /// which is what a spin is on a box - is dropped rather than
    /// regularised, because a mode that generates nothing cannot cancel
    /// anything either.
    fn strain_columns(m: &[f64], rows: usize) -> (Vec<f64>, Vec<f64>, usize) {
        let col_norm = |src: &[f64], c: usize, stride: usize| {
            (0..rows)
                .map(|i| src[i * stride + c] * src[i * stride + c])
                .sum::<f64>()
                .sqrt()
        };
        let scale = (6..12).fold(0.0f64, |a, c| a.max(col_norm(m, c, 12)));
        let mut q: Vec<Vec<f64>> = Vec::new();
        let mut dropped = 0usize;
        for c in 6..12 {
            let mut v: Vec<f64> = (0..rows).map(|i| m[i * 12 + c]).collect();
            for _pass in 0..2 {
                for u in q.iter() {
                    let ip: f64 = v.iter().zip(u.iter()).map(|(a, b)| a * b).sum();
                    for (a, b) in v.iter_mut().zip(u.iter()) {
                        *a -= ip * b;
                    }
                }
            }
            let nrm = v.iter().map(|a| a * a).sum::<f64>().sqrt();
            if nrm <= 1e-12 * scale {
                dropped += 1;
                continue;
            }
            for a in v.iter_mut() {
                *a /= nrm;
            }
            q.push(v);
        }
        let mut raw = vec![0.0f64; rows * 6];
        for i in 0..rows {
            for c in 0..6 {
                raw[i * 6 + c] = m[i * 12 + c];
            }
        }
        let mut proj = raw.clone();
        for c in 0..6 {
            for _pass in 0..2 {
                for u in q.iter() {
                    let ip: f64 = (0..rows).map(|i| proj[i * 6 + c] * u[i]).sum();
                    for i in 0..rows {
                        proj[i * 6 + c] -= ip * u[i];
                    }
                }
            }
        }
        (raw, proj, dropped)
    }

    /// `|B x|`: the weighted RMS face datum a unit strain mode produces,
    /// evaluated in the data space rather than through a quadratic form.
    fn mode_visibility(block: &[f64], rows: usize, m: &[[f64; 3]; 3]) -> f64 {
        let basis = strain_basis();
        let x: Vec<f64> = basis.iter().map(|s| frob(s, m)).collect();
        (0..rows)
            .map(|i| {
                let v: f64 = (0..6).map(|c| block[i * 6 + c] * x[c]).sum();
                v * v
            })
            .sum::<f64>()
            .sqrt()
    }

    /// The cell's own box frame, built from its facet centroids rather than
    /// assumed: `[e_1, e_2, e_r]` from (alpha+ - alpha-), (beta+ - beta-),
    /// (up - down), Gram-Schmidt'd radial-first. Also returns the raw
    /// non-orthogonality `|cos|` of the two lateral directions, so a claim
    /// about "the shear that is invisible on a box" can be checked against
    /// how box-like the cell actually is.
    fn cell_frame(t: &ShellTopology, c: usize) -> ([[f64; 3]; 3], f64) {
        let f = t.cell_faces(c as u32);
        let ctr = |idx: usize| t.faces()[f[idx].0 as usize].centroid_m;
        let a = sub(ctr(1), ctr(0));
        let b = sub(ctr(3), ctr(2));
        let r = sub(ctr(5), ctr(4));
        let cos = dot(a, b) / (dot(a, a) * dot(b, b)).sqrt();
        let unit = |v: [f64; 3]| {
            let s = 1.0 / dot(v, v).sqrt();
            [v[0] * s, v[1] * s, v[2] * s]
        };
        let ur = unit(r);
        let proj = |v: [f64; 3], u: [f64; 3]| {
            let d = dot(v, u);
            [v[0] - d * u[0], v[1] - d * u[1], v[2] - d * u[2]]
        };
        let u1 = unit(proj(a, ur));
        let u2 = unit(proj(proj(b, ur), u1));
        ([u1, u2, ur], cos.abs())
    }

    /// `(u_a (x) u_b + u_b (x) u_a) / sqrt(2)`: a unit-Frobenius shear in the
    /// cell's own frame. `mode(0, 1)` is the exact local analogue of
    /// `u = (y, x, 0)`, the strain a box cannot see.
    fn frame_mode(frame: &[[f64; 3]; 3], a: usize, b: usize) -> [[f64; 3]; 3] {
        let mut m = [[0.0f64; 3]; 3];
        let s = if a == b { 1.0 } else { S2 };
        for d in 0..3 {
            for e in 0..3 {
                m[d][e] = s * (frame[a][d] * frame[b][e] + frame[b][d] * frame[a][e]);
            }
        }
        if a == b {
            for row in m.iter_mut() {
                for v in row.iter_mut() {
                    *v *= 0.5;
                }
            }
        }
        m
    }

    /// Cell classes. The four lateral classes are reported only for cells
    /// that are NOT in a radial extreme layer, so that a lateral effect and
    /// a radial one cannot be confused; the two radial classes cover every
    /// lateral position. A cell can therefore be in at most two classes.
    #[derive(Clone, Copy, PartialEq, Eq, Debug)]
    enum Klass {
        PanelCentre,
        Interior,
        Seam,
        Corner,
        RadialInner,
        RadialOuter,
    }

    const KLASSES: [(Klass, &str); 6] = [
        (Klass::PanelCentre, "panel centre"),
        (Klass::Interior, "panel interior"),
        (Klass::Seam, "cube-face seam"),
        (Klass::Corner, "cube corner"),
        (Klass::RadialInner, "innermost layer"),
        (Klass::RadialOuter, "outermost layer"),
    ];

    fn classes_of(t: &ShellTopology, c: usize) -> Vec<Klass> {
        let n = t.spec().n();
        let layers = t.spec().layers();
        let (_, i, j, k) = t.spec().fijk(c as u64).expect("cell index in range");
        let mut out = Vec::with_capacity(2);
        if k == 0 {
            out.push(Klass::RadialInner);
        }
        if k + 1 == layers {
            out.push(Klass::RadialOuter);
        }
        if !out.is_empty() {
            return out;
        }
        let seam = i == 0 || j == 0 || i + 1 == n || j + 1 == n;
        if seam {
            out.push(Klass::Seam);
            if (i == 0 || i + 1 == n) && (j == 0 || j + 1 == n) {
                out.push(Klass::Corner);
            }
        } else {
            out.push(Klass::Interior);
            let mid = n / 2;
            if (i == mid || i + 1 == mid) && (j == mid || j + 1 == mid) {
                out.push(Klass::PanelCentre);
            }
        }
        out
    }

    /// One class's worst case at one resolution. Everything is a worst case
    /// over the cells of the class, never a mean: a cubed-sphere failure
    /// that averages away is not a pass.
    #[derive(Clone, Copy, Debug)]
    struct ClassStat {
        cells: usize,
        h: f64,
        /// `B_E` as the review defines it: the six Frobenius-orthonormal
        /// strain modes against the face data they produce, with NO rigid
        /// motion allowed to compensate.
        raw_min: f64,
        raw_max: f64,
        raw_kappa: f64,
        raw_rank: usize,
        /// The same after the rigid columns are projected out: the least
        /// data a strain can produce when a rigid motion of ANY amplitude is
        /// free to imitate it. Strictly smaller, and the honest coercivity
        /// constant wherever the fit is overdetermined.
        min: f64,
        max: f64,
        kappa: f64,
        rank: usize,
        /// Visibility of the local `(e1, e2)` shear - the exact analogue of
        /// `u = (y, x, 0)`, the mode a box cannot see - and of the radial
        /// normal strain, which a box CAN see, as the contrast.
        raw_shear: f64,
        shear: f64,
        raw_normal: f64,
        /// `|cos|` of the two lateral facet-centroid directions: how far
        /// from a box the cell actually is. The whole gate turns on whether
        /// this shrinks under refinement.
        skew: f64,
    }

    impl ClassStat {
        fn empty() -> Self {
            Self {
                cells: 0,
                h: 0.0,
                raw_min: f64::MAX,
                raw_max: 0.0,
                raw_kappa: 0.0,
                raw_rank: 6,
                min: f64::MAX,
                max: 0.0,
                kappa: 0.0,
                rank: 6,
                raw_shear: f64::MAX,
                shear: f64::MAX,
                raw_normal: f64::MAX,
                skew: 0.0,
            }
        }
    }

    /// Least-squares slope of `log y` against `log x`: the observed power of
    /// `h`. Reported instead of a two-point ratio because the ladder's
    /// radial refinement (3, 6, 10, 20 layers) is not an exact halving.
    fn log_slope(x: &[f64], y: &[f64]) -> f64 {
        let n = x.len() as f64;
        let lx: Vec<f64> = x.iter().map(|v| v.ln()).collect();
        let ly: Vec<f64> = y.iter().map(|v| v.ln()).collect();
        let mx = lx.iter().sum::<f64>() / n;
        let my = ly.iter().sum::<f64>() / n;
        let num: f64 = lx
            .iter()
            .zip(ly.iter())
            .map(|(a, b)| (a - mx) * (b - my))
            .sum();
        let den: f64 = lx.iter().map(|a| (a - mx) * (a - mx)).sum();
        num / den
    }

    /// The refinement ladder for the conditioning gates. `n = 4, 8, 16` is
    /// the ladder every other series in this module uses and it runs inside
    /// the normal debug suite; `MC0_GATE_LADDER=full` adds `n = 32`
    /// (122,880 cells), which is a release-mode run and is where the
    /// evidence file's `n = 32` column comes from.
    ///
    /// This is the runtime knob rather than a deleted case: the coverage is
    /// identical at every rung and only the top rung is opt-in.
    fn gate_ladder() -> Vec<(usize, usize)> {
        let mut v = vec![(4usize, 3usize), (8, 6), (16, 10)];
        if matches!(std::env::var("MC0_GATE_LADDER").as_deref(), Ok("full")) {
            v.push((32, 20));
        }
        v
    }

    /// Rank of an ascending singular-value list at a relative cut.
    fn numeric_rank(sig: &[f64]) -> usize {
        sig.iter()
            .filter(|v| **v > 1e-10 * sig[sig.len() - 1])
            .count()
    }

    /// Every per-class statistic of one stencil at one resolution.
    fn scan_shell(t: &ShellTopology, stencil: FitStencil) -> Vec<ClassStat> {
        let geo = FitGeometry::build(t).expect("fit geometry");
        let mut stats = vec![ClassStat::empty(); KLASSES.len()];
        let mut fit = CellFit::default();
        for c in 0..t.cell_count() {
            cell_fit_system(t, &geo, c, stencil, &mut fit);
            let rows = fit.rows.len();
            let m = observability_map(&fit);
            let (raw, proj, _dropped) = strain_columns(&m, rows);
            let sig = one_sided_jacobi(&mut proj.clone(), rows, 6);
            let raw_sig = one_sided_jacobi(&mut raw.clone(), rows, 6);
            let (frame, skew) = cell_frame(t, c);
            let shear_mode = frame_mode(&frame, 0, 1);
            let normal_mode = frame_mode(&frame, 2, 2);
            for k in classes_of(t, c) {
                let idx = KLASSES.iter().position(|(x, _)| *x == k).expect("class");
                let s = &mut stats[idx];
                s.cells += 1;
                s.h += fit.h;
                s.raw_min = s.raw_min.min(raw_sig[0]);
                s.raw_max = s.raw_max.max(raw_sig[5]);
                s.raw_kappa = s.raw_kappa.max(raw_sig[5] / raw_sig[0]);
                s.raw_rank = s.raw_rank.min(numeric_rank(&raw_sig));
                s.min = s.min.min(sig[0]);
                s.max = s.max.max(sig[5]);
                s.kappa = s.kappa.max(sig[5] / sig[0]);
                s.rank = s.rank.min(numeric_rank(&sig));
                s.raw_shear = s.raw_shear.min(mode_visibility(&raw, rows, &shear_mode));
                s.shear = s.shear.min(mode_visibility(&proj, rows, &shear_mode));
                s.raw_normal = s.raw_normal.min(mode_visibility(&raw, rows, &normal_mode));
                s.skew = s.skew.max(skew);
            }
        }
        for s in stats.iter_mut() {
            if s.cells > 0 {
                s.h /= s.cells as f64;
            }
        }
        stats
    }

    #[test]
    fn gate3_strain_observability_and_its_scaling_under_refinement() {
        // GATE 3. Exact reproduction of an affine field at finite `h` is NOT
        // evidence that this operator is safe under a Krylov solver, and this
        // test exists to stop that inference from being made.
        //
        // The counterexample the gate is built on: on a centred cube of
        // half-width `a`, the constant shear `u = (y, x, 0)` has `u . n = +-y`
        // on the `x = +-a` faces and `+-x` on the `y = +-a` faces, and the
        // INTEGRAL of `y` (resp. `x`) over those faces is exactly zero. Exact
        // facet integration does not rescue it: one normal-flux scalar per
        // face on a six-face box cannot see any of the three shear modes,
        // however exactly it is integrated. A cubed-sphere cell is not a box,
        // so a shear is weakly visible - but refinement makes the cell MORE
        // box-like, so the lever arm that visibility is borrowed from
        // shrinks. **The verdict is the scaling and nothing else.**
        //
        // Measured per cell: the observability map, i.e. the six
        // Frobenius-orthonormal constant strain modes against the weighted
        // face data each one produces, normalised so that a mode's amplitude
        // is `h_c |eps|_F` (a velocity) and a column's norm is a weighted RMS
        // face datum (the same velocity). Both `h` and `R` divide out, so
        // what is left is observability. Two versions:
        //
        // - `raw` is `B_E` exactly as the review defines it: strain modes to
        //   face-normal degrees of freedom, no rigid motion allowed to
        //   compensate.
        // - the projected version additionally lets a rigid motion of ANY
        //   amplitude imitate the strain's data, which is the sharper
        //   question wherever the fit is overdetermined - and is vacuous on a
        //   six-face cell, where twelve affine parameters cannot be
        //   determined from six numbers at all. That count, not the geometry,
        //   is why the local `rank` column reads 1: reported, and not leaned
        //   on.
        //
        // Three stencils, because the gate is a question about WHICH stencil
        // and not about the fit:
        //   1. the cell's own six facets - the local MFD map;
        //   2. the cell's own interior faces - a legal purely local operator;
        //   3. cell plus face neighbours - what this module's operator uses.
        let ladder = gate_ladder();
        /// Per stencil: the worst (least collapsing) class, and the generic
        /// panel-interior class on its own, because a cube corner and a
        /// panel interior answer this question differently and averaging
        /// them would hide both answers.
        struct Verdict {
            name: String,
            worst_raw: f64,
            interior_raw: f64,
            interior_shear: f64,
            worst_proj: f64,
            proj_floor: f64,
        }
        let mut summary: Vec<Verdict> = Vec::new();
        for (stencil, name) in [
            (
                FitStencil::CellFacets,
                "B_E: the cell's own six facets, the local MFD map",
            ),
            (
                FitStencil::CellInteriorFaces,
                "the cell's own interior faces: a legal purely-local operator",
            ),
            (
                FitStencil::CellAndFaceNeighbours,
                "cell + face neighbours: WHAT THIS OPERATOR USES",
            ),
        ] {
            println!("\n=== {name} ===");
            println!(
                "{:<16}{:>3}{:>7}{:>10}| {:>10}{:>10}{:>10}{:>5} | {:>10}{:>10}{:>10}{:>5} | {:>10}{:>10}{:>10}{:>9}",
                "class", "n", "cells", "mean h",
                "raw s_min", "raw s_max", "raw kappa", "rk",
                "s_min", "s_max", "kappa", "rk",
                "raw shear", "proj shear", "raw normal", "|cos|"
            );
            let mut series: Vec<Vec<ClassStat>> = Vec::new();
            for (n, layers) in ladder.iter().copied() {
                let t = mini(n, layers);
                let stats = scan_shell(&t, stencil);
                for (idx, (_, label)) in KLASSES.iter().enumerate() {
                    let s = stats[idx];
                    if s.cells == 0 {
                        continue;
                    }
                    println!(
                        "{label:<16}{n:>3}{:>7}{:>10.3e}| {:>10.3e}{:>10.3e}{:>10.3e}{:>5} | \
                         {:>10.3e}{:>10.3e}{:>10.3e}{:>5} | {:>10.3e}{:>10.3e}{:>10.3e}{:>9.2e}",
                        s.cells,
                        s.h,
                        s.raw_min,
                        s.raw_max,
                        s.raw_kappa,
                        s.raw_rank,
                        s.min,
                        s.max,
                        s.kappa,
                        s.rank,
                        s.raw_shear,
                        s.shear,
                        s.raw_normal,
                        s.skew
                    );
                }
                series.push(stats);
            }
            println!("-- observed powers of h (least squares over the ladder) --");
            let live = |v: &[f64]| v.iter().all(|x| *x > 1e-30);
            let fmt = |ok: bool, v: f64| {
                if ok {
                    format!("h^{v:+.2}")
                } else {
                    "NUMERICALLY ZERO".to_string()
                }
            };
            let (mut worst_raw, mut worst_proj) = (f64::MAX, f64::MAX);
            let mut min_proj_min = f64::MAX;
            let (mut interior_raw, mut interior_shear) = (f64::NAN, f64::NAN);
            for (idx, (_, label)) in KLASSES.iter().enumerate() {
                if series.iter().any(|s| s[idx].cells == 0) {
                    continue;
                }
                let col = |f: fn(&ClassStat) -> f64| -> Vec<f64> {
                    series.iter().map(|s| f(&s[idx]).max(1e-300)).collect()
                };
                let h = col(|s| s.h);
                let raw = col(|s| s.raw_min);
                let rawk = col(|s| s.raw_kappa);
                let rshear = col(|s| s.raw_shear);
                let rnormal = col(|s| s.raw_normal);
                let pmin = col(|s| s.min);
                let skew = col(|s| s.skew);
                let last = h.len() - 1;
                if live(&raw) {
                    worst_raw = worst_raw.min(log_slope(&h, &raw));
                }
                if live(&pmin) {
                    worst_proj = worst_proj.min(log_slope(&h, &pmin));
                    min_proj_min = min_proj_min.min(pmin.iter().cloned().fold(f64::MAX, f64::min));
                }
                if KLASSES[idx].0 == Klass::Interior {
                    interior_raw = log_slope(&h, &raw);
                    interior_shear = log_slope(&h, &rshear);
                }
                println!(
                    "{label:<16} raw: s_min {} ({:.3e} -> {:.3e})  kappa {}  shear {}  normal {}  \
                     | projected: c1 = s_min^2 {:.3e} -> {:.3e}, s_min {}  | cell |cos| {}",
                    fmt(live(&raw), log_slope(&h, &raw)),
                    raw[0],
                    raw[last],
                    fmt(live(&rawk), log_slope(&h, &rawk)),
                    fmt(live(&rshear), log_slope(&h, &rshear)),
                    fmt(live(&rnormal), log_slope(&h, &rnormal)),
                    pmin[0] * pmin[0],
                    pmin[last] * pmin[last],
                    fmt(live(&pmin), log_slope(&h, &pmin)),
                    fmt(live(&skew), log_slope(&h, &skew)),
                );
            }
            summary.push(Verdict {
                name: format!("{stencil:?}"),
                worst_raw,
                interior_raw,
                interior_shear,
                worst_proj,
                proj_floor: min_proj_min,
            });
        }

        println!("\n-- verdict --");
        for v in summary.iter() {
            let proj = if v.worst_proj == f64::MAX {
                "numerically zero in every class".to_string()
            } else {
                format!(
                    "worst class h^{:+.2}, never below {:.3e}",
                    v.worst_proj, v.proj_floor
                )
            };
            println!(
                "{:<24} raw B_E sigma_min: panel interior h^{:+.2}, least-collapsing class \
                 h^{:+.2}   raw local (e1,e2) shear, panel interior h^{:+.2}   | projected \
                 sigma_min: {proj}",
                v.name, v.interior_raw, v.worst_raw, v.interior_shear
            );
        }

        // The assertions. Both directions are pinned, because both are
        // findings: the local map DOES degenerate, and the stencil the
        // operator actually uses does not. Neither bound is fitted - a
        // collapse is a power near or above 1 with the constant falling
        // through the floor, and these are separated from that by the whole
        // interval between them.
        let local = &summary[0];
        assert!(
            local.interior_raw > 2.0,
            "the local six-facet map's smallest singular value is NOT collapsing in the panel \
             interior (h^{:+.2}, measured h^+3.26 when this bound was written); the box \
             counterexample no longer holds on this grid and every conclusion drawn from it \
             needs re-deriving",
            local.interior_raw
        );
        assert!(
            local.worst_raw > 0.8,
            "some cell class escaped the local degeneracy entirely (h^{:+.2}); the least \
             collapsing class measured h^+1.12, the cube corners, and it collapses because \
             even a corner cell's distortion cannot outrun refinement",
            local.worst_raw
        );
        assert!(
            local.interior_shear > 1.2,
            "the local visibility of the (e1,e2) shear - the exact analogue of u = (y, x, 0) - \
             is NOT collapsing in the panel interior (h^{:+.2}, measured h^+1.83)",
            local.interior_shear
        );
        let op = &summary[2];
        let op_power = op.worst_proj;
        let op_floor = op.proj_floor;
        assert!(
            op_power.abs() < 0.3,
            "the OPERATOR's strain-subspace coercivity constant is scaling with h \
             (sigma_min ~ h^{op_power:+.2}). If that power is positive and near 1, c1 is \
             collapsing under refinement and a Krylov solver must NOT be wrapped around this \
             operator until the local constitutive reconstruction is enriched."
        );
        assert!(
            op_floor > 0.05,
            "the OPERATOR's smallest rigid-projected strain singular value fell to {op_floor:e}; \
             it was 1.07e-1 at the worst class and resolution when this bound was written"
        );
    }

    /// The single-cell map exactly as the review's own independent check
    /// defines it: for the linear field `u(x) = G x` the exact integrated
    /// flux through facet `f` is `G : T_f`, so the map from the six
    /// Frobenius-orthonormal symmetric strain modes to the six facet fluxes
    /// has rows `sym(T_f)` contracted with each mode. No fit, no
    /// normalisation, and - when `centred` is false - no centroid
    /// subtraction, so the two implementations can be compared number for
    /// number.
    ///
    /// `centred` replaces `T_f` with `T_f - A_f (x) x_c`, i.e. it asks about
    /// the field `G (x - x_c)` instead of `G x`. That is exactly one change,
    /// and it is the difference between "is the local strain visible" and
    /// "is the strain visible given that the same field also carries the
    /// translation `G x_c`, which on a shell is `R/h` times larger than
    /// anything local".
    fn review_moment_map(
        t: &ShellTopology,
        geo: &FitGeometry,
        c: usize,
        centred: bool,
    ) -> Vec<f64> {
        let basis = strain_basis();
        let x_c = geo.cell_centroids_m()[c];
        let mut m = vec![0.0f64; 36];
        for (row, (f, _)) in t.cell_faces(c as u32).into_iter().enumerate() {
            let tf = geo.facet_tensor_moments()[f as usize];
            let af = t.faces()[f as usize].area_vector_m2;
            let mut g = [[0.0f64; 3]; 3];
            for (d, grow) in g.iter_mut().enumerate() {
                for (e, slot) in grow.iter_mut().enumerate() {
                    *slot = tf[d][e] - if centred { af[d] * x_c[e] } else { 0.0 };
                }
            }
            for (k, mode) in basis.iter().enumerate() {
                // `G : T = sym(G) : sym(T)` for a symmetric `G`, so
                // contracting the mode against `T` already restricts to the
                // symmetric part.
                m[row * 6 + k] = frob(mode, &g);
            }
        }
        m
    }

    /// The same map for a centred unit cube, built analytically: facet
    /// `+-x` has `T = integral n (x) x dA = |A| a (e_x (x) e_x)`, and
    /// likewise for `y` and `z`, so the six rows are three duplicated pairs.
    /// This is the review's validation case and it is here so that the two
    /// implementations are known to be computing the same object before
    /// their cubed-sphere numbers are compared.
    fn unit_cube_moment_map() -> Vec<f64> {
        let basis = strain_basis();
        let a = 0.5f64;
        let area = 1.0f64;
        let mut m = vec![0.0f64; 36];
        for axis in 0..3 {
            for side in 0..2 {
                let mut g = [[0.0f64; 3]; 3];
                g[axis][axis] = area * a;
                let row = 2 * axis + side;
                for (k, mode) in basis.iter().enumerate() {
                    m[row * 6 + k] = frob(mode, &g);
                }
            }
        }
        m
    }

    #[test]
    fn gate3_the_single_cell_moment_map_is_rank_three_and_the_defect_is_real() {
        // The review ran this check independently in Python, through the
        // normal equations, and reported "sigma_min either exactly 0 or of
        // order 1e-8 against a sigma_max of 10 to 80, at EVERY resolution",
        // with the honest caveat that a Gram-plus-Jacobi route cannot
        // distinguish 1e-8 relative from zero. This test answers the
        // question that caveat leaves open, with a one-sided Jacobi SVD that
        // resolves small singular values to high RELATIVE accuracy.
        //
        // The control comes first: on a centred unit cube the map must have
        // exactly three zeros and three singular values of `1/sqrt(2)`. If
        // that does not reproduce, nothing below is about the same object.
        let cube = one_sided_jacobi(&mut unit_cube_moment_map(), 6, 6);
        println!(
            "centred unit cube, singular values: {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}",
            cube[0], cube[1], cube[2], cube[3], cube[4], cube[5]
        );
        for (k, s) in cube.iter().take(3).enumerate() {
            assert!(
                *s < 1e-16,
                "unit cube singular value {k} is {s:e}, not zero: the three shear modes must be \
                 EXACTLY invisible on a box"
            );
        }
        for (k, s) in cube.iter().skip(3).enumerate() {
            assert!(
                (s - S2).abs() < 1e-15,
                "unit cube singular value {} is {s:e}, not 1/sqrt(2)",
                k + 3
            );
        }

        // Now the real cells. Reported per class is the LEAST degenerate
        // cell - the one whose third-smallest singular value is largest
        // relative to its largest - because that cell is the strongest
        // available counterexample to "the map is rank three", and its whole
        // relative spectrum is printed so that "genuinely zero or merely
        // tiny" is answered for all three of the small values rather than
        // for the smallest one.
        println!(
            "\n{:<16}{:>3}{:>4}{:>11}{:>11}{:>11}{:>11}{:>11}{:>11}{:>11}",
            "class", "n", "cvn", "s1/s6", "s2/s6", "s3/s6", "s4/s6", "s5/s6", "s6", "worst kappa"
        );
        let ladder = gate_ladder();
        let mut series: Vec<Vec<[f64; 6]>> = Vec::new();
        let mut hs: Vec<Vec<f64>> = Vec::new();
        for (n, layers) in ladder.iter().copied() {
            let t = mini(n, layers);
            let geo = FitGeometry::build(&t).expect("fit geometry");
            // Per class, per convention: the relative spectrum of the least
            // degenerate cell, and the cell's sigma_max.
            let mut best = vec![[[0.0f64; 6]; 2]; KLASSES.len()];
            let mut smax = vec![[0.0f64; 2]; KLASSES.len()];
            // The strongest counterexample in the other direction: the cell
            // whose SMALLEST relative singular value is largest. If that is
            // still at roundoff, the null direction is exact everywhere.
            let mut worst_null = vec![0.0f64; KLASSES.len()];
            let mut h_class = vec![(0.0f64, 0usize); KLASSES.len()];
            for c in 0..t.cell_count() {
                let spectra = [
                    one_sided_jacobi(&mut review_moment_map(&t, &geo, c, false), 6, 6),
                    one_sided_jacobi(&mut review_moment_map(&t, &geo, c, true), 6, 6),
                ];
                for k in classes_of(&t, c) {
                    let idx = KLASSES.iter().position(|(x, _)| *x == k).expect("class");
                    worst_null[idx] = worst_null[idx].max(spectra[0][0] / spectra[0][5]);
                    for (cv, sig) in spectra.iter().enumerate() {
                        if sig[2] / sig[5] > best[idx][cv][2] {
                            for j in 0..6 {
                                best[idx][cv][j] = sig[j] / sig[5];
                            }
                            smax[idx][cv] = sig[5];
                        }
                    }
                    h_class[idx].0 += geo.polyhedral_cell_volumes_m3()[c].cbrt();
                    h_class[idx].1 += 1;
                }
            }
            for (idx, (_, label)) in KLASSES.iter().enumerate() {
                if h_class[idx].1 == 0 {
                    continue;
                }
                for (cv, tag) in ["G x", "ctr"].iter().enumerate() {
                    let b = best[idx][cv];
                    println!(
                        "{label:<16}{n:>3}{tag:>4}{:>11.3e}{:>11.3e}{:>11.3e}{:>11.3e}{:>11.3e}{:>11.3e}{:>11.3e}",
                        b[0], b[1], b[2], b[3], b[4], smax[idx][cv], 1.0 / b[0]
                    );
                }
            }
            hs.push(
                h_class
                    .iter()
                    .map(|(sum, n)| if *n == 0 { 0.0 } else { sum / *n as f64 })
                    .collect(),
            );
            series.push(
                best.iter()
                    .zip(worst_null.iter())
                    .map(|(b, w)| [b[0][2], b[1][2], *w, 0.0, 0.0, 0.0])
                    .collect(),
            );
        }
        println!(
            "-- the THIRD smallest relative singular value (the last of the three a box cannot \
             see), and its power of h --"
        );
        let mut worst_null = 0.0f64;
        let mut weakest_pair = f64::MAX;
        let mut slowest_decay = f64::MAX;
        for (idx, (_, label)) in KLASSES.iter().enumerate() {
            if hs.iter().any(|h| h[idx] == 0.0) {
                continue;
            }
            let h: Vec<f64> = hs.iter().map(|x| x[idx]).collect();
            let raw: Vec<f64> = series.iter().map(|s| s[idx][0].max(1e-300)).collect();
            let ctr: Vec<f64> = series.iter().map(|s| s[idx][1].max(1e-300)).collect();
            worst_null = worst_null.max(series.iter().map(|s| s[idx][2]).fold(0.0, f64::max));
            weakest_pair = weakest_pair.min(raw.iter().cloned().fold(f64::MAX, f64::min));
            slowest_decay = slowest_decay.min(log_slope(&h, &raw));
            println!(
                "{label:<16} G x : {:.3e} -> {:.3e}  ~ h^{:+.2}   |   G (x - x_c) : {:.3e} -> \
                 {:.3e}  ~ h^{:+.2}",
                raw[0],
                raw[raw.len() - 1],
                log_slope(&h, &raw),
                ctr[0],
                ctr[ctr.len() - 1],
                log_slope(&h, &ctr)
            );
        }
        println!(
            "\nmeasured: exactly ONE relative singular value at roundoff (worst cell over every \
             class and rung: {worst_null:.3e}), and a weak PAIR at {weakest_pair:.3e} or above, \
             decaying no faster than h^{slowest_decay:+.2}"
        );

        // Two claims are pinned here, and they are different claims.
        //
        // (1) The null space of the single-cell map is EXACTLY ONE
        //     dimensional on our cells, not three. On a box it is three (the
        //     control above), so the curvature does buy back two of the
        //     three shear directions - but only two, and only weakly. The
        //     review's independent Python run reported "rank 3 to within
        //     1e-8"; that is refuted as stated, and the cause is the normal
        //     equations, which cannot separate one exact zero from three
        //     values at 1e-2. The practical conclusion is unaffected and is
        //     if anything stronger: an exact null direction is not a
        //     conditioning problem that refinement can improve.
        //
        // (2) The two directions that curvature does buy back are weak and
        //     they DO collapse under refinement. So the single-cell map is
        //     unusable at every resolution and becomes more so.
        assert!(
            worst_null < 1e-14,
            "the single-cell moment map has NO exactly null direction on some cell \
             (largest smallest-relative-singular-value {worst_null:e}); measured 5.1e-17 over \
             every class and rung when this bound was written, and the whole rank argument \
             needs re-deriving if it has moved"
        );
        assert!(
            weakest_pair > 1e-3,
            "the second and third smallest singular values have fallen to {weakest_pair:e}, i.e. \
             to the level of the exactly null direction: the map is now rank 3 rather than \
             rank 5 and the report of two weakly-visible shears no longer holds"
        );
        assert!(
            slowest_decay > 0.4,
            "the weak pair is NOT collapsing under refinement (h^{slowest_decay:+.2}); measured \
             h^+0.62 at the slowest class"
        );
    }

    // ---------------------------------------------------------------
    // Gate 1: what the null space of A actually is IN THIS discrete space

    #[test]
    fn gate1_the_null_space_is_six_and_three_of_its_directions_are_not_rigid_motions() {
        // The findings file said the null space is six dimensional, and that
        // the impermeable condition is carried by `D` rather than by `A`.
        // The second half is wrong and this test is the measurement that
        // replaces it.
        //
        // Impermeability is an ESSENTIAL boundary constraint on the velocity
        // SPACE, not a property of the divergence operator: a rigid
        // translation is divergence free, so `D u = 0` cannot be what
        // excludes it. In a conforming velocity space one would impose
        // `u . n = 0` on the boundary and translations would leave the space.
        // Here the constraint is imposed by EXCLUDING the boundary faces from
        // the unknown vector, and then the interesting thing happens:
        //
        //   raw symmetric-gradient A on a space with boundary unknowns
        //       6 rigid modes (3 translations + 3 rotations)
        //   OUR space, boundary faces already excluded
        //       STILL 6 - and this is the correction the review asked for.
        //       Three of them are rotations, which really are rigid and
        //       really are impermeable. The other three are TRUNCATED
        //       translations: a translation on the interior faces and zero on
        //       the boundary faces. That is not a rigid motion, it is not
        //       divergence free, and A cannot see the difference only because
        //       A never reads a boundary face.
        //   full free-slip spherical-shell Stokes
        //       3 rotations + 1 constant pressure, because `D u = 0` removes
        //       the truncated translations - measured below as the rank of D
        //       on the six-dimensional null space of A.
        //
        // So our code is ALREADY in the after-restriction regime and it has
        // six, not three. A multigrid preconditioner built on the A block
        // alone must deflate SIX.
        let mut weak: Vec<f64> = Vec::new();
        let mut boundary_defect: Vec<f64> = Vec::new();
        for (n, layers) in [(4usize, 3usize), (8, 6)] {
            let (t, a) = op(n, layers);
            let v_sphere = t.cell_volumes_m3().expect("volumes");
            let moments = t.facet_tensor_moments().expect("facet tensor moments");
            let modes = rigid_modes(&t);
            let names = [
                "translation x",
                "translation y",
                "translation z",
                "rotation x",
                "rotation y",
                "rotation z",
            ];
            println!(
                "\n{:<14}{:>3}{:>12}{:>14}{:>14}{:>14}{:>14}",
                "mode",
                "n",
                "|A u|/ref",
                "bnd u.n / int",
                "|D u| int",
                "|D u| bnd",
                "cells with Du"
            );
            let reference = a
                .apply(&affine_face_field(
                    &t,
                    [0.0; 3],
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                ))
                .expect("apply");
            // Each mode is normalised to unit peak face value BEFORE its
            // divergence is taken. Without that the comparison is not one:
            // a unit-omega rotation field has face values of order |omega| R
            // while a unit translation has order 1, and the Gram would be
            // reporting the planet's radius rather than the two modes'
            // relative incompressibility.
            let raw_scale: Vec<f64> = modes.iter().map(|m| max_abs(m)).collect();
            let modes: Vec<Vec<f64>> = modes
                .iter()
                .zip(raw_scale.iter())
                .map(|(m, s)| m.iter().map(|v| v / s).collect())
                .collect();
            let mut d_rows: Vec<Vec<f64>> = Vec::with_capacity(6);
            for (k, m) in modes.iter().enumerate() {
                let au = a.apply(m).expect("apply");
                let annihilation = max_abs(&au) / max_abs(&reference) * raw_scale[k];
                // The value this mode was TRUNCATED to zero at, on the
                // impermeable faces, relative to its own interior scale. For
                // a rotation the true normal component is already almost
                // zero, so the truncation is a small correction; for a
                // translation it is the whole datum, which is why the
                // resulting mode is not a rigid motion.
                // The mode is normalised, but the exact boundary datum below
                // is built in the mode's own physical units, so it is divided
                // by the amplitude that normalisation removed.
                let interior_scale = raw_scale[k];
                let mut boundary_truth = 0.0f64;
                for (idx, face) in t.faces().iter().enumerate() {
                    if face.neighbor.is_some() {
                        continue;
                    }
                    let truth = if k < 3 {
                        // A translation's exact facet-mean normal component.
                        let mut b = [0.0f64; 3];
                        b[k] = 1.0;
                        dot(b, face.normal)
                    } else {
                        // A rotation's, from the exact facet moment.
                        let mut omega = [0.0f64; 3];
                        omega[k - 3] = 1.0;
                        let mut acc = 0.0;
                        for (d, trow) in moments[idx].iter().enumerate() {
                            for (e, tv) in trow.iter().enumerate() {
                                acc += rotation_l(omega)[d][e] * tv;
                            }
                        }
                        acc / face.area_m2
                    };
                    boundary_truth = boundary_truth.max(truth.abs());
                }
                let div = t.divergence(m, &v_sphere).expect("divergence");
                let (mut interior, mut boundary, mut hit) = (0.0f64, 0.0f64, 0usize);
                let scale = 1.0 / (v_sphere[0].cbrt());
                for (c, d) in div.iter().enumerate() {
                    let touches = t
                        .cell_faces(c as u32)
                        .iter()
                        .any(|(f, _)| t.faces()[*f as usize].neighbor.is_none());
                    if touches {
                        boundary = boundary.max(d.abs() / scale);
                        if d.abs() / scale > 1e-10 {
                            hit += 1;
                        }
                    } else {
                        interior = interior.max(d.abs() / scale);
                    }
                }
                println!(
                    "{:<14}{n:>3}{annihilation:>12.3e}{:>14.3e}{interior:>14.3e}{boundary:>14.3e}{hit:>14}",
                    names[k],
                    boundary_truth / interior_scale
                );
                if k == 3 {
                    boundary_defect.push(boundary_truth / interior_scale);
                }
                assert!(
                    annihilation < 1e-11,
                    "{}: A does not annihilate it ({annihilation:e}), so the null space is not \
                     what this test claims",
                    names[k]
                );
                d_rows.push(div);
            }

            // The rank of D restricted to null(A): the dimension of the part
            // of the A-block null space that survives the incompressibility
            // constraint. Measured as the spectrum of the 6x6 Gram of the six
            // `D u_k` in the M_p inner product.
            let mut gram = vec![0.0f64; 36];
            for i in 0..6 {
                for j in 0..6 {
                    gram[i * 6 + j] = (0..t.cell_count())
                        .map(|c| v_sphere[c] * d_rows[i][c] * d_rows[j][c])
                        .sum();
                }
            }
            // The Gram's singular values are the SQUARES of D's on this
            // subspace, so report the square roots: those are relative sizes
            // of `D u` itself.
            let sig = one_sided_jacobi(&mut gram.clone(), 6, 6);
            let rel: Vec<f64> = sig.iter().map(|s| (s / sig[5]).sqrt()).collect();
            let gap = rel[2] / rel[3];
            println!(
                "n={n}: D restricted to null(A), relative singular values: {:.3e} {:.3e} {:.3e} \
                 | {:.3e} {:.3e} {:.3e}   (gap {gap:.3e})",
                rel[0], rel[1], rel[2], rel[3], rel[4], 1.0
            );
            // NOT rank 3. This is the second correction the gate produced and
            // it was not predicted: D has rank SIX on the null space of A,
            // with a five-order gap. Three strong directions are the
            // truncated translations, whose deleted boundary flux is the
            // whole datum. Three WEAK directions at ~1e-5 relative are the
            // rotations: the exact rotation is tangent to the SPHERE, but a
            // boundary facet is a flat panel that is not, so the rotation
            // carries a real O(h^2/R^2) flux through it - which the
            // impermeable condition then deletes. No field in this discrete
            // space is both a rotation and exactly divergence free in a
            // boundary cell.
            //
            // For a solver that matters more than the exact modes do: an
            // exact null mode can be deflated, a NEAR-null one has to be
            // resolved or it sets the iteration count.
            assert!(
                gap < 1e-1,
                "n={n}: D's action on null(A) no longer splits into three strong and three weak \
                 directions (gap {gap:e}); measured 8.7e-3 at n=4 when this was written"
            );
            for (k, r) in rel.iter().take(3).enumerate() {
                assert!(
                    *r > 1e-9,
                    "n={n}: rotation direction {k} is EXACTLY divergence free ({r:e}); that would \
                     be a better world than the measured one and the finding needs re-deriving"
                );
            }
            weak.push(rel[2]);

            // And the pressure side of the same statement.
            let constant = vec![7.0f64; t.cell_count()];
            let g = t.gradient(&constant).expect("gradient");
            assert_eq!(
                max_abs(&g),
                0.0,
                "n={n}: the constant pressure mode is not in the null space of G"
            );
        }

        // Both defects are second order, so they converge away - but they are
        // present at every resolution a solver will ever run at.
        let order = |v: &[f64]| (v[0] / v[1]).log2();
        println!(
            "\nrotation boundary-normal datum, relative: {:.3e} -> {:.3e} (order {:+.2})",
            boundary_defect[0],
            boundary_defect[1],
            order(&boundary_defect)
        );
        println!(
            "weak (rotation) singular value of D on null(A): {:.3e} -> {:.3e} (order {:+.2})",
            weak[0],
            weak[1],
            order(&weak)
        );
        assert!(
            order(&weak) > 1.5,
            "the rotations' residual discrete divergence is not second order ({:+.2}); it is the \
             deleted boundary flux and it should converge away with the facet's departure from \
             the sphere",
            order(&weak)
        );
        assert!(
            order(&boundary_defect) > 1.5,
            "the rotation's boundary-normal datum is not converging at second order \
             ({:+.2}); it is the polyhedral facet's departure from the sphere and it should",
            order(&boundary_defect)
        );
    }

    // ---------------------------------------------------------------
    // Gate 4: certifying the rotation residual as roundoff, by an argument
    // that implies it

    /// Round `x` to `bits` significant mantissa bits, round to nearest.
    /// `bits >= 53` is the identity.
    fn round_bits(x: f64, bits: u32) -> f64 {
        if bits >= 53 || x == 0.0 || !x.is_finite() {
            return x;
        }
        let shift = 53 - bits;
        let half = 1u64 << (shift - 1);
        // IEEE-754 is sign-magnitude, so adding half an ulp to the raw bits
        // rounds the MAGNITUDE up whatever the sign, and a mantissa carry
        // correctly increments the exponent.
        f64::from_bits(x.to_bits().wrapping_add(half) & (!0u64 << shift))
    }

    /// One cell's strain, recomputed from geometry with EVERY intermediate
    /// rounded to `bits` mantissa bits: rows, weights, normal matrix,
    /// Cholesky, solve, coefficients, face data and the final scatter.
    ///
    /// This is a deliberate second implementation of the assembly. A
    /// precision experiment has to be able to vary the precision, and the
    /// operator's own path is fixed at f64; the guard against the two
    /// drifting apart is that at `bits = 53` this function must reproduce
    /// `StrainOperator::apply` on the same field, which the test asserts
    /// before it reads anything from the ladder.
    fn rotation_strain_at_precision(
        t: &ShellTopology,
        geo: &FitGeometry,
        fit: &CellFit,
        cell: usize,
        l: [[f64; 3]; 3],
        bits: u32,
    ) -> StrainVoigt {
        let r = |v: f64| round_bits(v, bits);
        let x_c = geo.cell_centroids_m()[cell];
        let h = r(r(geo.polyhedral_cell_volumes_m3()[cell]).cbrt());
        let m = fit.faces.len();
        let mut rows = vec![[0.0f64; FIT_DOF]; m];
        let mut weights = vec![0.0f64; m];
        let mut data = vec![0.0f64; m];
        let mut gram = [[0.0f64; FIT_DOF]; FIT_DOF];
        for (i, &f) in fit.faces.iter().enumerate() {
            let face = &t.faces()[f as usize];
            let tf = &geo.facet_tensor_moments()[f as usize];
            let inv_a = r(1.0 / r(face.area_m2));
            for (d, slot) in rows[i][..3].iter_mut().enumerate() {
                *slot = r(face.normal[d]);
            }
            let mut datum = 0.0f64;
            for (d, trow) in tf.iter().enumerate() {
                for (e, tv) in trow.iter().enumerate() {
                    rows[i][3 + 3 * d + e] = r(r(r(r(
                        r(*tv) - r(r(face.area_vector_m2[d]) * r(x_c[e]))
                    ) * inv_a)
                        / h));
                    datum = r(datum + r(r(l[d][e]) * r(*tv)));
                }
            }
            data[i] = r(datum * inv_a);
            let dx = [
                r(r(face.centroid_m[0]) - r(x_c[0])),
                r(r(face.centroid_m[1]) - r(x_c[1])),
                r(r(face.centroid_m[2]) - r(x_c[2])),
            ];
            let d2 = r(r(r(dx[0] * dx[0]) + r(dx[1] * dx[1])) + r(dx[2] * dx[2]));
            weights[i] = r(r(h * h) / d2);
            for a in 0..FIT_DOF {
                for b in 0..FIT_DOF {
                    gram[a][b] = r(gram[a][b] + r(weights[i] * r(rows[i][a] * rows[i][b])));
                }
            }
        }
        // Cholesky and solve, rounded at every step.
        let mut lo = [[0.0f64; FIT_DOF]; FIT_DOF];
        for i in 0..FIT_DOF {
            for j in 0..=i {
                let mut s = gram[i][j];
                for (a, b) in lo[i].iter().take(j).zip(lo[j].iter()) {
                    s = r(s - r(a * b));
                }
                if i == j {
                    lo[i][i] = r(s.sqrt());
                } else {
                    lo[i][j] = r(s / lo[j][j]);
                }
            }
        }
        let mut eps = [0.0f64; VOIGT_LEN];
        for i in 0..m {
            let mut y = [0.0f64; FIT_DOF];
            for a in 0..FIT_DOF {
                let mut s = rows[i][a];
                for k in 0..a {
                    s = r(s - r(lo[a][k] * y[k]));
                }
                y[a] = r(s / lo[a][a]);
            }
            for a in (0..FIT_DOF).rev() {
                let mut s = y[a];
                for k in (a + 1)..FIT_DOF {
                    s = r(s - r(lo[k][a] * y[k]));
                }
                y[a] = r(s / lo[a][a]);
            }
            for (slot, (a, d)) in eps.iter_mut().zip(VOIGT_INDEX.iter()) {
                let g = r(r(y[3 + 3 * a + d]) + r(y[3 + 3 * d + a]));
                let coeff = r(r(r(weights[i] * 0.5) * g) / h);
                *slot = r(*slot + r(coeff * data[i]));
            }
        }
        eps
    }

    #[test]
    fn gate4_the_rotation_residual_scales_with_the_precision_and_not_with_h() {
        // REPLACES a false inference. The previous certification argued that
        // the rigid-rotation residual is roundoff because the series RISES
        // with resolution, "which a discretization error cannot do". That
        // inference is false and this build contains its counterexample: the
        // inconsistent two-point gradient's error also grows under
        // refinement (`shell_ops`:
        // `raw_two_point_dg_laplacian_diverges_under_refinement`) and it is
        // pure discretization error. A growing series is evidence of
        // NOTHING on its own.
        //
        // The replacement argument: at FIXED h, vary the arithmetic
        // precision and watch the residual. A discretization error is a
        // property of the mesh and the scheme and is precision-independent
        // by definition; a roundoff error is proportional to the machine
        // epsilon. So the residual is roundoff if and only if it falls like
        // `2^-bits` while `h` is held fixed. Measured over a 29-bit range,
        // which is a factor of 5.4e8 in machine epsilon.
        //
        // LIMITS, stated. (a) Rounding every intermediate to `bits` mantissa
        // bits is not identical to computing in a native `bits`-wide format:
        // each operation is rounded twice, once to f64 and once to `bits`.
        // Double rounding changes the constant in front of the error, not
        // its proportionality to the epsilon, which is what is claimed here.
        // (b) The geometry INPUTS are f64 values from `shell_ops`, rounded to
        // `bits` before use. So this certifies the fit's arithmetic, not the
        // construction of `T_f` - and that is the right scope, because the
        // rotation annihilation is exact for whatever `T_f` values are
        // supplied, the same numbers appearing in the design row and in the
        // face datum. An error in `T_f` itself cancels between them.
        // (c) The ladder is not run above 53 bits, so this shows the residual
        // falling to the f64 floor rather than continuing below it.
        let omega = [0.7, -1.3, 0.4];
        let l = rotation_l(omega);
        let scale = max_abs(&omega);
        let ladder = [24u32, 29, 34, 39, 44, 48, 53];
        // The guard field: a general affine gradient with an O(1) symmetric
        // part. The two assemblies must agree on THAT, because it is a real
        // number rather than a cancellation residue - comparing two
        // roundoff-level residuals would only ever compare two summation
        // orders.
        let guard_l = [[0.4, -1.1, 0.3], [0.9, 0.25, -0.6], [-0.2, 0.8, 1.4]];
        let mut guard_expected = [0.0f64; VOIGT_LEN];
        for (slot, (a, d)) in guard_expected.iter_mut().zip(VOIGT_INDEX.iter()) {
            *slot = 0.5 * (guard_l[*a][*d] + guard_l[*d][*a]);
        }
        let guard_scale = max_abs(&guard_expected);
        for (n, layers) in [(4usize, 3usize), (8, 6)] {
            let (t, a) = op(n, layers);
            let geo = FitGeometry::build(&t).expect("fit geometry");
            let u = t.rigid_rotation_field(omega);
            let operator = a
                .strain()
                .apply(&u)
                .expect("strain")
                .iter()
                .fold(0.0f64, |acc, e| acc.max(max_abs(e)))
                / scale;
            let mut fit = CellFit::default();
            let guard_ref = a
                .strain()
                .apply(&affine_face_field(&t, [0.0; 3], guard_l))
                .expect("strain");
            let mut guard_worst = 0.0f64;
            for (c, want) in guard_ref.iter().enumerate() {
                cell_fit_system(&t, &geo, c, FitStencil::CellAndFaceNeighbours, &mut fit);
                let got = rotation_strain_at_precision(&t, &geo, &fit, c, guard_l, 53);
                for (j, v) in got.iter().enumerate() {
                    guard_worst = guard_worst.max((v - want[j]).abs() / guard_scale);
                    guard_worst = guard_worst.max((v - guard_expected[j]).abs() / guard_scale);
                }
            }
            println!(
                "\nn={n}: precision-ladder assembly vs StrainOperator on an O(1) affine strain, \
                 worst relative difference {guard_worst:.3e}"
            );
            assert!(
                guard_worst < 1e-12,
                "n={n}: the precision-ladder assembly is not the operator's computation - it \
                 differs by {guard_worst:e} on a field where both should return sym(L) exactly"
            );
            let mut series = Vec::new();
            for bits in ladder {
                let mut worst = 0.0f64;
                for c in 0..t.cell_count() {
                    cell_fit_system(&t, &geo, c, FitStencil::CellAndFaceNeighbours, &mut fit);
                    let eps = rotation_strain_at_precision(&t, &geo, &fit, c, l, bits);
                    worst = worst.max(max_abs(&eps) / scale);
                }
                series.push(worst);
            }
            let eps_of = |bits: u32| 2.0f64.powi(-(bits as i32));
            println!(
                "\nn={n} layers={layers}: rigid-rotation strain residual against arithmetic \
                 precision, h FIXED"
            );
            for (bits, res) in ladder.iter().zip(series.iter()) {
                println!(
                    "  {bits:>2} mantissa bits (eps {:.2e}): residual {res:.4e}   residual/eps \
                     {:.3e}",
                    eps_of(*bits),
                    res / eps_of(*bits)
                );
            }
            let x: Vec<f64> = ladder.iter().map(|b| eps_of(*b)).collect();
            let slope = log_slope(&x, &series);
            println!(
                "  slope of log(residual) against log(machine epsilon): {slope:+.3}  \
                 (1.0 = pure roundoff, 0.0 = precision-independent, i.e. discretization)"
            );
            println!(
                "  operator's own f64 residual {operator:.4e} vs this path at 53 bits {:.4e}",
                series[series.len() - 1]
            );

            // At 53 bits the two paths differ only in which exact
            // construction supplies the rotation's face data - the stored
            // antisymmetric moment vector in `rigid_rotation_field` against
            // the full tensor moment here - so the two residuals are two
            // roundoff residues of the same zero and agree only in order of
            // magnitude. That is why the guard above is on `sym(L)` instead.
            let agreement = series[series.len() - 1] / operator;
            assert!(
                (0.1..10.0).contains(&agreement),
                "n={n}: the two f64 residuals differ by {agreement:e}x, which is more than two \
                 summation orders of the same zero can explain"
            );
            assert!(
                slope > 0.9,
                "n={n}: the rotation residual does NOT scale with the machine epsilon \
                 (slope {slope:+.3}); it is therefore not certified as roundoff, and since a \
                 rising series proves nothing either, nothing certifies it"
            );
            // And the other half of the same statement, which is what makes
            // the slope meaningful: a discretization error would have to be
            // there at every precision.
            assert!(
                series[0] / series[series.len() - 1] > 1e6,
                "n={n}: the residual barely moved between 24 and 53 bits, so a \
                 precision-independent floor - a discretization error - dominates it"
            );
        }
    }

    // ---------------------------------------------------------------
    // Gate 2: is `tr(E u) != D u` an artificial bulk viscosity?

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    /// The four vertices of each facet, in [`Facet::ALL`] order, as indices
    /// into `ShellSpec::cell_vertices`. Copied from `shell_ops::facet_quad`,
    /// which is private - and checked against the shipped geometry on every
    /// facet by [`quadrature_face_field`], which asserts that the two
    /// triangles' vector areas sum to the stored `area_vector_m2`. If the
    /// upstream table ever changes, that check fails rather than this
    /// silently integrating over a different surface.
    const FACET_QUAD_PICK: [[usize; 4]; 6] = [
        [0, 3, 7, 4],
        [1, 2, 6, 5],
        [0, 1, 5, 4],
        [3, 2, 6, 7],
        [0, 1, 2, 3],
        [4, 5, 6, 7],
    ];

    /// Exact face data for a field that is a polynomial of degree <= 2:
    /// `u_f = (1/|A_f|) integral_f v . n dA`, integrated over the SAME two
    /// triangles per facet that `facet_tensor_moments` uses, with the
    /// three-point edge-midpoint rule, which is exact for a quadratic.
    ///
    /// The point of integrating over that exact surface is the divergence
    /// theorem: for a closed cell the signed facet fluxes then sum to
    /// `integral_cell div v dV` identically, so a solenoidal `v` gives
    /// `D u = 0` to ROUNDOFF rather than to `O(h^2)`. That is what makes a
    /// statement about `tr(E u)` on the constraint set possible at all.
    fn quadrature_face_field(
        t: &ShellTopology,
        v: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> (Vec<f64>, f64) {
        let mut u = vec![0.0f64; t.face_count()];
        let mut worst_area = 0.0f64;
        for c in 0..t.cell_count() {
            let (f, i, j, k) = t.spec().fijk(c as u64).expect("cell index in range");
            let vert = t.spec().cell_vertices(f, i, j, k).expect("cell vertices");
            for (idx, (face_id, sign)) in t.cell_faces(c as u32).into_iter().enumerate() {
                if sign != 1 {
                    continue; // this cell is the neighbour, not the owner
                }
                let face = &t.faces()[face_id as usize];
                let q: [[f64; 3]; 4] = FACET_QUAD_PICK[idx].map(|p| vert[p]);
                let tri = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
                    let av = cross(sub(b, a), sub(c, a));
                    let av = [0.5 * av[0], 0.5 * av[1], 0.5 * av[2]];
                    let mid = |x: [f64; 3], y: [f64; 3]| {
                        [
                            0.5 * (x[0] + y[0]),
                            0.5 * (x[1] + y[1]),
                            0.5 * (x[2] + y[2]),
                        ]
                    };
                    let mut mean = [0.0f64; 3];
                    for p in [mid(a, b), mid(b, c), mid(c, a)] {
                        let val = v(p);
                        for (m, x) in mean.iter_mut().zip(val.iter()) {
                            *m += x / 3.0;
                        }
                    }
                    (av, dot(av, mean))
                };
                let (a0, f0) = tri(q[0], q[1], q[2]);
                let (a1, f1) = tri(q[0], q[2], q[3]);
                let winding = [a0[0] + a1[0], a0[1] + a1[1], a0[2] + a1[2]];
                let s = if dot(winding, face.area_vector_m2) < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                // The copied vertex table is only trustworthy if it produces
                // the shipped facet.
                let mismatch = (0..3)
                    .map(|d| (s * winding[d] - face.area_vector_m2[d]).abs())
                    .fold(0.0f64, f64::max)
                    / face.area_m2;
                worst_area = worst_area.max(mismatch);
                if face.neighbor.is_some() {
                    u[face_id as usize] = s * (f0 + f1) / face.area_m2;
                }
            }
        }
        (u, worst_area)
    }

    #[test]
    fn gate2_the_volumetric_mismatch_is_a_convergent_truncation_not_a_bulk_viscosity() {
        // The physics of the worry. Incompressible Newtonian Stokes relies on
        // `tr(eps) = div u = 0`, so the energy `2 eta eps:eps` carries no
        // independent volumetric penalty. If the discrete velocity satisfies
        // `D u = 0` while `tr(E u) != 0`, then `E^T W E` penalises a
        // GEOMETRIC fake volume change that the pressure constraint cannot
        // see: a small artificial bulk viscosity inside an incompressible
        // mantle.
        //
        // Two manufactured fields, both exactly solenoidal and both of degree
        // two so the facet quadrature is exact:
        //   toroidal   v = (x z, -y z, y^2 - x^2) / R,  which also satisfies
        //              v . x = 0, so it is tangent to every sphere and its
        //              flux through the impermeable boundary is only the
        //              polyhedral facet's departure from the sphere;
        //   generic    v = (y z, x z, -2 x y) / R,  which is not tangent, so
        //              its boundary cells carry a real deleted flux and are
        //              reported separately rather than averaged in.
        let r = R_TOP;
        let toroidal = move |x: [f64; 3]| {
            [
                x[0] * x[2] / r,
                -x[1] * x[2] / r,
                (x[1] * x[1] - x[0] * x[0]) / r,
            ]
        };
        let generic = move |x: [f64; 3]| [x[1] * x[2] / r, x[0] * x[2] / r, -2.0 * x[0] * x[1] / r];
        for (name, field) in [
            ("toroidal", &toroidal as &dyn Fn([f64; 3]) -> [f64; 3]),
            ("generic", &generic),
        ] {
            println!("\n=== {name} solenoidal field ===");
            println!(
                "{:>3}{:>12}{:>12}{:>12}{:>12}{:>12}{:>12}{:>12}",
                "n",
                "|D u| int",
                "|D u| bnd",
                "tr Linf",
                "tr seam",
                "tr rms",
                "vol/Phi int",
                "(Phi-Phi~)/Phi"
            );
            let mut hs = Vec::new();
            let mut traces = Vec::new();
            let mut traces_rms = Vec::new();
            let mut fractions = Vec::new();
            for (n, layers) in [(4usize, 3usize), (8, 6), (16, 10)] {
                let (t, a) = op(n, layers);
                let (u, area_mismatch) = quadrature_face_field(&t, field);
                assert!(
                    area_mismatch < 1e-12,
                    "n={n}: the facet vertex table copied from shell_ops no longer produces the \
                     shipped facet (relative area-vector mismatch {area_mismatch:e})"
                );
                let v_sphere = t.cell_volumes_m3().expect("volumes");
                let div = t.divergence(&u, &v_sphere).expect("divergence");
                let eps = a.strain().apply(&u).expect("strain");
                let w = a.cell_weights();

                // Scale: the RMS strain of this field, so every ratio below
                // is against the strain that is actually present.
                let eps_rms = (eps.iter().map(|e| strain_contraction(e, e)).sum::<f64>()
                    / eps.len() as f64)
                    .sqrt();
                let (mut div_int, mut div_bnd, mut tr_int) = (0.0f64, 0.0f64, 0.0f64);
                let (mut tr_sq, mut tr_cells) = (0.0f64, 0usize);
                let (mut tr_seam, mut tr_face) = (0.0f64, 0.0f64);
                let (mut vol_int, mut tot_int) = (0.0f64, 0.0f64);
                let (mut vol_all, mut tot_all, mut compat) = (0.0f64, 0.0f64, 0.0f64);
                for (c, e) in eps.iter().enumerate() {
                    let touches = t
                        .cell_faces(c as u32)
                        .iter()
                        .any(|(f, _)| t.faces()[*f as usize].neighbor.is_none());
                    let tr = strain_trace(e);
                    let total = w[c] * strain_contraction(e, e);
                    let volumetric = w[c] * tr * tr / 3.0;
                    // The trace made compatible with the SAME discrete
                    // divergence the pressure constraint uses.
                    let compatible = w[c] * div[c] * div[c] / 3.0;
                    tot_all += total;
                    vol_all += volumetric;
                    compat += total - volumetric + compatible;
                    if touches {
                        div_bnd = div_bnd.max(div[c].abs() / eps_rms);
                    } else {
                        div_int = div_int.max(div[c].abs() / eps_rms);
                        tr_int = tr_int.max(tr.abs() / eps_rms);
                        if seam_touching(&t, c) {
                            tr_seam = tr_seam.max(tr.abs() / eps_rms);
                        } else {
                            tr_face = tr_face.max(tr.abs() / eps_rms);
                        }
                        tr_sq += tr * tr;
                        tr_cells += 1;
                        tot_int += total;
                        vol_int += volumetric;
                    }
                }
                let tr_rms = (tr_sq / tr_cells as f64).sqrt() / eps_rms;
                println!(
                    "{n:>3}{div_int:>12.3e}{div_bnd:>12.3e}{tr_face:>12.3e}{tr_seam:>12.3e}{tr_rms:>12.3e}{:>12.3e}{:>12.3e}",
                    vol_int / tot_int,
                    (tot_all - compat) / tot_all
                );
                let _ = vol_all / tot_all;
                hs.push(t.spec().cell_angle_rad() * R_TOP);
                traces.push(tr_int);
                traces_rms.push(tr_rms);
                fractions.push(vol_int / tot_int);
                // The premise of the whole gate: the field really is
                // discretely divergence free, so `tr(E u) != 0` cannot be
                // blamed on the field.
                assert!(
                    div_int < 1e-12,
                    "n={n} {name}: the manufactured field is not discretely divergence free in \
                     the interior ({div_int:e} of the RMS strain); the quadrature or the \
                     triangulation is wrong and nothing below means anything"
                );
            }
            let order = |v: &[f64]| log_slope(&hs, v);
            println!(
                "  max |tr(E u)| / |eps|_rms over interior cells ~ h^{:+.2}   ({:.3e} -> {:.3e})",
                order(&traces),
                traces[0],
                traces[2]
            );
            println!(
                "  RMS |tr(E u)| / |eps|_rms over interior cells ~ h^{:+.2}   ({:.3e} -> {:.3e})",
                order(&traces_rms),
                traces_rms[0],
                traces_rms[2]
            );
            println!(
                "  spurious volumetric FRACTION of the viscous energy ~ h^{:+.2}   ({:.3e} -> \
                 {:.3e})",
                order(&fractions),
                fractions[0],
                fractions[2]
            );
            // The worst cell's spurious trace converges at FIRST order, not
            // second - measured h^+0.90 - which is the same Linf-versus-L2
            // split `shell_ops` already reports for the corrected Laplacian.
            // It converges, which is the question the gate asks; the rate is
            // recorded rather than assumed.
            // The worst CELL's spurious trace converges slowly and at a rate
            // that depends on the field - measured h^+0.90 for the toroidal
            // one and h^+0.35 for the non-tangent one - so what is asserted
            // is that it falls at every rung, not a rate this packet has not
            // established. Whether that Linf series plateaus, as the seam
            // Linf of the corrected Laplacian does in `shell_ops`, is stated
            // as unresolved in the evidence file rather than papered over
            // with a bound that happens to hold on three points.
            for w in traces.windows(2) {
                assert!(
                    w[1] < w[0],
                    "{name}: the worst-cell spurious trace stopped falling under refinement \
                     ({:.3e} -> {:.3e}); a non-convergent trace IS a persistent artificial bulk \
                     viscosity and the operator would have to change",
                    w[0],
                    w[1]
                );
            }
            assert!(
                order(&traces_rms) > 0.85,
                "{name}: the RMS spurious trace is not converging at first order (h^{:+.2}); it \
                 is the RMS that enters the energy, and measured h^+1.60 (toroidal) and h^+0.97 \
                 (generic)",
                order(&traces_rms)
            );
            assert!(
                order(&fractions) > 1.8,
                "{name}: the spurious volumetric ENERGY fraction is not converging (h^{:+.2}); \
                 measured h^+3.26 (toroidal) and h^+2.00 (generic), which is the square of the \
                 RMS trace in both cases",
                order(&fractions)
            );
            assert!(
                fractions.iter().all(|f| *f < 1e-3),
                "{name}: the spurious volumetric energy reached {:.3e} of the viscous energy; \
                 measured 6.7e-5 at worst when this bound was written, and a percent-level \
                 value would be a real bulk viscosity rather than a truncation error",
                fractions.iter().cloned().fold(0.0, f64::max)
            );
        }
    }

    #[test]
    fn gate2_the_trace_identity_holds_on_affine_data_only_and_is_not_an_operator_identity() {
        // The findings file states `tr(eps_c) = (V_sphere/V_poly) (D u)_c`
        // as though it were an identity between operators. It is not: it is
        // a statement about AFFINE data, where the least-squares fit has a
        // zero residual and therefore reproduces the Green-Gauss flux sum
        // exactly. `E` reads 28 to 34 faces and `D` reads six; on data that
        // is not affine there is no reason for the two to agree, and they do
        // not.
        //
        // This matters for gate 2 rather than being a pedantic correction:
        // if the identity DID hold as an operator identity, then `D u = 0`
        // would force `tr(E u) = 0` exactly and there would be no volumetric
        // mismatch to worry about at all.
        let (t, a) = op(4, 3);
        let v_sphere = t.cell_volumes_m3().expect("volumes");
        let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
        let mut rng = Lcg(0x0bad_1dea_0bad_1dea);
        let mut worst_affine = 0.0f64;
        let mut worst_random = 0.0f64;
        for trial in 0..20 {
            let u = rng.field(&t);
            let eps = a.strain().apply(&u).expect("strain");
            let div = t.divergence(&u, &v_sphere).expect("divergence");
            let scale = eps
                .iter()
                .fold(0.0f64, |acc, e| acc.max(strain_trace(e).abs()));
            for (c, e) in eps.iter().enumerate() {
                if t.cell_faces(c as u32)
                    .iter()
                    .any(|(f, _)| t.faces()[*f as usize].neighbor.is_none())
                {
                    continue;
                }
                let predicted = v_sphere[c] / v_poly[c] * div[c];
                worst_random = worst_random.max((strain_trace(e) - predicted).abs() / scale);
            }
            let _ = trial;
        }
        // The same measurement on affine data, which is where the identity
        // does hold, so the two numbers sit side by side.
        let l = [
            [1.3e-15, -0.4e-15, 0.2e-15],
            [0.1e-15, 0.7e-15, -0.9e-15],
            [0.5e-15, 0.3e-15, -0.6e-15],
        ];
        let u = affine_face_field(&t, [0.0; 3], l);
        let eps = a.strain().apply(&u).expect("strain");
        let div = t.divergence(&u, &v_sphere).expect("divergence");
        let scale = (l[0][0] + l[1][1] + l[2][2]).abs();
        for (c, e) in eps.iter().enumerate() {
            if t.cell_faces(c as u32)
                .iter()
                .any(|(f, _)| t.faces()[*f as usize].neighbor.is_none())
            {
                continue;
            }
            let predicted = v_sphere[c] / v_poly[c] * div[c];
            worst_affine = worst_affine.max((strain_trace(e) - predicted).abs() / scale);
        }
        println!(
            "tr(E u) - (V_sphere/V_poly) D u, worst interior cell, relative to the trace scale:\n  \
             affine data {worst_affine:.3e}   random data {worst_random:.3e}"
        );
        assert!(
            worst_random > 1e-2,
            "the trace identity now holds on random data too ({worst_random:e}); it would then \
             be an operator identity after all and gate 2 collapses to nothing"
        );
        assert!(
            worst_affine < 1e-9,
            "the trace identity no longer holds on affine data ({worst_affine:e})"
        );
    }

    // ---------------------------------------------------------------
    // Fail loud

    #[test]
    fn wrongly_sized_and_inadmissible_inputs_are_refused() {
        let t = mini(2, 2);
        let a = ViscousOperator::isoviscous(&t, ETA).expect("operator");
        let err = a
            .apply(&[0.0; 3])
            .expect_err("a short field must be refused");
        assert!(format!("{err}").contains("length 3"), "error was: {err}");
        let err = ViscousOperator::isoviscous(&t, 0.0).expect_err("zero viscosity is refused");
        assert!(format!("{err}").contains("eta"), "error was: {err}");
        let err = ViscousOperator::isoviscous(&t, f64::NAN).expect_err("NaN viscosity is refused");
        assert!(format!("{err}").contains("eta"), "error was: {err}");
        let err = ViscousOperator::from_cell_weights(&t, vec![1.0; 3])
            .expect_err("a short weight vector is refused");
        assert!(format!("{err}").contains("length 3"), "error was: {err}");
    }
}
