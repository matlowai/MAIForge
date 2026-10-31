//! Mimetic discrete operators on the cubed-sphere shell (MC prerequisite).
//!
//! Brief rev 5 "MC design and qualification" / [ADR 026]. This module owns
//! the face list and the compatible divergence / gradient pair that a
//! staggered Stokes discretization needs. It contains NO physics: no
//! viscosity, no buoyancy, no time. That separation is deliberate, because
//! every claim here is provable from geometry alone and is therefore
//! testable without believing anything about the mantle.
//!
//! # Two gradients live here. Use [`ShellTopology::gradient_corrected`].
//!
//! [`ShellTopology::gradient`], the plain two-point normal difference, is
//! compatible but NOT consistent, and it is kept only because pinning that
//! down is worth a test. `D` and `G` are exactly adjoint, and that turns out
//! to say nothing about whether `D G` approximates a Laplacian. It does not.
//! With the diagonal `M_u = diag(|A_f| d_f)`, the two-point gradient is only
//! ZEROTH order accurate on this grid, because the facet normal and the
//! centroid-to-centroid line point in different directions and a diagonal
//! Hodge star cannot correct for it (Thuburn and Cotter 2012, eq. 5.2).
//!
//! Measured on `p = (x/R)^2`, whose exact Laplacian is the constant `2/R^2`,
//! as relative error against it, at `n = 8, 16, 32` (layers 6, 10, 18). Three
//! series, because the brief requires it and because a two-series report
//! would miss half the story:
//!
//! ```text
//!                            n=8      n=16     n=32     observed order
//! raw two-point   interior   0.21924  0.35219  0.42584  -0.68  -0.27
//!                 seam       1.58261  3.78715  8.12281  -1.26  -1.10
//!                 L2         0.35257  0.50210  0.69292  -0.51  -0.46
//! corrected       interior   0.06595  0.02769  0.01372  +1.25  +1.01
//!                 seam       0.05879  0.04704  0.04849  +0.32  -0.04
//!                 L2         0.02184  0.00810  0.00327  +1.43  +1.31
//! control, exact face-normal derivative at the facet centroid:
//!                 interior   0.01973  0.00521  0.00135  +1.92  +1.95
//!                 L2         0.01684  0.00424  0.00106  +1.99  +2.00
//! ```
//!
//! The control isolates the original cause: same divergence, same volumes,
//! same area vectors, only the gradient replaced, and the error becomes clean
//! second order. So the defect was entirely in `G`, and therefore in the
//! diagonal `M_u` behind it.
//!
//! Counter-intuitively the SEAMS are the healthy part of this mesh for
//! non-orthogonality. Across a cube edge the two cells are exact mirror
//! images through the plane of the shared edge and the origin, so the
//! centroid segment is exactly perpendicular to the shared facet, measured at
//! 1e-6 degrees, which is roundoff. The non-orthogonality that breaks the
//! gradient lives in the FACE INTERIORS and grows toward 30 degrees (measured
//! 20.6, 25.2, 27.6 degrees at n = 8, 16, 32). What the seams have instead is
//! permanent SKEWNESS: the facet centroid sits about a quarter of the
//! centroid distance off the line joining the centroids and stays there under
//! refinement. Any error report that looks only at seams samples the cells
//! that are healthiest in one metric and sickest in the other.
//!
//! ## The fix, and exactly how far it goes
//!
//! `M_u` is FREE, which is what makes it affordable: `M_u,ff = |A_f| d_f /
//! lambda_f` yields `G_f = lambda_f (p_nb - p_own) / d_f` for any positive
//! `lambda`, and `G = -M_u^-1 D^T M_p` still holds EXACTLY. So the implicit
//! operator can be rescaled per face at no cost in compatibility, and the
//! remainder is carried explicitly, on the right-hand side, where it cannot
//! touch the adjoint identity.
//!
//! With the over-relaxed split `S = Delta + k`, `Delta = (S.S / (d.S)) d`
//! (Jasak), `|Delta| / |S| = lambda_f = 1 / (n_hat . d_hat)` and
//!
//! ```text
//! n_hat . grad p  =  lambda_f (p_nb - p_own) / d_f              [implicit]
//!                 +  (k_f / |A_f|) . (grad p)_face              [explicit]
//!                 +  lambda_f (g_nb - g_own) . r_f / d_f        [explicit]
//! ```
//!
//! where the third term re-centres the two-point difference from the midpoint
//! of the centroid segment - the point a difference quotient actually
//! differentiates about - to the facet centroid, `r_f` being the skewness
//! vector. The first two terms are the review's F1; the third is needed
//! because on this grid the seam skewness does not vanish under refinement.
//! `tests::the_residual_seam_error_is_skewness_not_non_orthogonality` is the
//! control that separates them: an EXACT face-normal derivative evaluated at
//! the segment midpoint instead of the facet centroid still saturates at 0.15
//! on the seam series while staying second order in the interior.
//!
//! **What this achieves, stated honestly.** At n=32 the corrected operator is
//! 31x better than the raw one in the interior Linf, 168x on the seam series
//! and 212x in L2. The L2 norm converges at order 1.3; the interior Linf at
//! order 1.0; the seam Linf falls and then SATURATES near 0.048. It does not
//! reach clean second order and this is not a tolerance that could be
//! loosened into one. Two facts locate the remainder:
//!
//! 1. Fed an exact cell gradient, the same operator IS clean second order in
//!    every series (interior 0.0197 -> 0.0052 -> 0.0014, L2 0.0169 -> 0.0043
//!    -> 0.0011), reproducing the exact-derivative control. The
//!    discretization is consistent.
//! 2. The compact cell gradient is first order in Linf and that is what is
//!    left: 1.79e-2 -> 7.88e-3 -> 3.85e-3 for the skew-seeded Green-Gauss
//!    form, 2.56e-2 -> 1.30e-2 -> 6.59e-3 for least squares. Re-seeding
//!    Green-Gauss with its own output converges to a fixed point 15 percent
//!    better and still first order, so iteration is not the missing
//!    ingredient. `|A_f| / V_c ~ 1/h` turns that first-order gradient error
//!    into the Linf plateau above.
//!
//! For calibration: Thuburn, Cotter and Dubos 2014 Table 4 report, for a
//! CONSISTENT non-diagonal Hodge star on a SMOOTHED cubed sphere, Linf rising
//! 0.21 -> 0.45 with resolution and "the L2 error is almost first order,
//! while the Linf error does not decrease". A saturating Linf with a
//! converging L2 is the published state of the art on this grid; the plateau
//! here is an order of magnitude below theirs. MC-Q2's pass criteria should
//! be set from these numbers, not from an assumption of second order.
//!
//! **Cost, for the solver.** The explicit term is large: `|k_f| / |A_f| =
//! tan(theta_f)`, measured max 0.3766, 0.4709, 0.5225 at n = 8, 16, 32,
//! rising toward `tan 30 deg = 0.577`. On the actual test field the realised
//! `|explicit| / |implicit|` at facets carrying a significant gradient
//! reaches 0.98-1.00. A deferred correction of that size is not something to
//! run unrelaxed; expect to need under-relaxation or to lag the correction
//! across outer iterations.
//!
//! # The unknowns
//!
//! Pressure lives at cell centres. Velocity lives on cell facets as the
//! single scalar normal component `u_f`, with the face's stored unit normal
//! fixing its sign. A cell-centred Cartesian velocity vector exists only as
//! a DERIVED reconstruction ([`ShellTopology::reconstruct_cell_velocity`]),
//! never as a primary unknown. This is what keeps the pressure and velocity
//! spaces compatible instead of collocated, so there is no checkerboard mode
//! to cure afterwards with Rhie-Chow machinery.
//!
//! # The compatibility identity
//!
//! With the mass matrices
//!
//! ```text
//! M_p = diag(cell volume)
//! M_u = diag(face area * centroid distance)
//! ```
//!
//! and the conservative divergence
//!
//! ```text
//! (D u)_c = (1 / V_c) * sum_f  s(c,f) * u_f * |A_f|
//! ```
//!
//! the weighted negative adjoint
//!
//! ```text
//! G = -M_u^-1 D^T M_p
//! ```
//!
//! works out to the ordinary two-point normal difference
//! `(p_neighbour - p_owner) / d_f`. That is not a coincidence to be relied
//! on silently: [`tests::gradient_is_the_weighted_negative_adjoint_of_divergence`]
//! assembles both operators as explicit matrices on a small shell and checks
//! the identity to machine precision, so if either operator is ever changed
//! the pair either stays compatible or the build fails.
//!
//! # Boundaries and nullspaces
//!
//! The inner and outer spherical surfaces are impermeable: their facets
//! carry no velocity unknown, and they are listed as boundary faces rather
//! than quietly omitted. Two nullspaces survive that and are removed
//! explicitly rather than left for a Krylov method to pick a representative
//! of:
//!
//! 1. pressure is determined only up to a constant
//!    ([`ShellTopology::remove_mean_pressure`]);
//! 2. with impermeable free-slip spherical boundaries, whole-mantle rigid
//!    rotation is physically null
//!    ([`ShellTopology::remove_rigid_rotation`]).
//!
//! The second one matters beyond tidiness: plate motion is later fitted to
//! the solved mantle velocity, so a spurious numerical spin would be
//! inherited as real plate motion and would show up in the world as a planet
//! whose plates all drift one way for no reason.
//!
//! [ADR 026]: ../../../../docs/adr/026_physical_mantle_simulation_critical_path.md

use crate::shell::{Facet, Lateral, ShellError, ShellSpec};

/// What kind of interface a face is. Boundary faces carry no unknown but are
/// listed, because an operator that silently drops them cannot be audited.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceKind {
    /// Between two cells in the same radial layer.
    Lateral,
    /// Between two radially adjacent cells.
    Radial,
    /// The core-mantle boundary surface.
    BoundaryInner,
    /// The lithosphere-base surface.
    BoundaryOuter,
}

impl FaceKind {
    pub fn is_boundary(self) -> bool {
        matches!(self, FaceKind::BoundaryInner | FaceKind::BoundaryOuter)
    }
}

/// The over-relaxed non-orthogonal split of one face: the geometry the
/// corrected gradient needs, and the diagnosis of why the raw one fails.
///
/// Built by [`ShellTopology::face_corrections`]. Boundary faces get the
/// identity entry (`lambda = 1`, zero tangential part) because they carry no
/// velocity unknown, so `face_masses_over_relaxed` agrees with `face_masses`
/// there and nothing downstream has to special-case them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceCorrection {
    /// `lambda_f = 1 / (n_hat_f . d_hat_f) >= 1`, the over-relaxation factor.
    /// Absorbed into `M_u,ff = |A_f| d_f / lambda_f`, which is why the
    /// rescaled implicit gradient is still EXACTLY `-M_u^-1 D^T M_p`.
    pub lambda: f64,
    /// The tangential remainder `k_f = S_f - Delta_f` of the area vector,
    /// m^2, where `Delta_f = (S.S / (d.S)) d` is the over-relaxed part
    /// parallel to the centroid line. `|k_f| = |S_f| tan(theta_f)`, so this
    /// vector's size relative to `|S_f|` is the size of the explicit
    /// correction relative to the implicit term.
    pub tangential_m2: [f64; 3],
    /// Unit owner-to-neighbour centroid direction `d_hat_f`.
    pub centroid_dir: [f64; 3],
    /// Linear interpolation weight for the OWNER value at the facet
    /// centroid: `phi_f = w phi_own + (1 - w) phi_nb`.
    pub owner_weight: f64,
    /// The INTERPOLATION offset `x_facet_centroid - (w x_own + (1-w) x_nb)`,
    /// m: how far the facet centroid sits off the point on the centroid
    /// segment that `owner_weight` actually interpolates to. A face value
    /// built by linear interpolation lands at that point, not at the facet
    /// centroid, and on this grid the gap does not close at the seams.
    pub interp_offset_m: [f64; 3],
    /// The SKEWNESS vector `x_facet_centroid - (x_own + x_nb)/2`, m: how far
    /// the facet centroid sits off the midpoint of the centroid segment.
    ///
    /// A two-point difference differentiates about the midpoint of its own
    /// segment, so this is the displacement between where the implicit term
    /// is accurate and where the flux is actually wanted. It matters here
    /// because on this grid it does NOT shrink relative to `d_f` at the
    /// seams: measured 0.236 -> 0.262 -> 0.275 of `d_f` at n = 8, 16, 32,
    /// rising to a limit near 0.28, while the in-face value does shrink.
    pub skew_m: [f64; 3],
    /// `acos(n_hat_f . d_hat_f)`, radians. Reported rather than hidden,
    /// because it is the entire reason this struct exists: it reaches 30
    /// degrees in the cube-face interiors and does not shrink with h.
    pub non_orthogonality_rad: f64,
}

/// How the face value of the cell-centred gradient is built for the explicit
/// non-orthogonal correction. Both rules are measured in the test suite; this
/// is a stated choice, not a tuning knob.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FaceGradientRule {
    /// Review F1 exactly as ranked: `(k_f / |A_f|) . grad_f` with `grad_f`
    /// the distance-weighted average of the two cell gradients. Corrects the
    /// NON-ORTHOGONALITY of the facet and nothing else.
    NonOrthogonalOnly,
    /// F1 with `grad_f` first shifted so its component along the centroid
    /// line equals the two-point difference. Measured, and it is not an
    /// improvement here: it injects the skew-biased two-point value into the
    /// tangential term. Kept because the measurement is worth keeping.
    NonOrthogonalCentroidLine,
    /// F1 plus the SKEWNESS term `lambda_f (g_nb - g_own) . r_f / d_f`, which
    /// moves the two-point difference from the midpoint of the centroid
    /// segment to the facet centroid where the flux is actually wanted.
    ///
    /// Beyond F1 as ranked, and it is what makes the operator converge: F1
    /// alone corrects the direction of the gradient but not the POINT it is
    /// evaluated at, and on this grid the seam skewness does not vanish under
    /// refinement, so the uncorrected part saturates. Still purely explicit,
    /// so the adjoint identity is untouched.
    #[default]
    NonOrthogonalAndSkew,
}

/// One face of the mesh, stored once. Geometry is recorded from the OWNER
/// cell's side: `normal` and `area_vector_m2` point out of `owner`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Face {
    pub owner: u32,
    /// `None` for a boundary face.
    pub neighbor: Option<u32>,
    pub kind: FaceKind,
    pub area_m2: f64,
    pub area_vector_m2: [f64; 3],
    /// Unit normal, out of `owner`.
    pub normal: [f64; 3],
    /// Owner-to-neighbour centroid distance; for a boundary face, the
    /// distance from the owner centroid to the face centroid.
    pub distance_m: f64,
    /// Face centroid, the mean of its four vertices. A declared
    /// approximation: for a curved facet this is not the exact centroid, and
    /// it is used only for interpolation weights and diagnostics, never for
    /// a conserved quantity.
    pub centroid_m: [f64; 3],
    /// The exact surface moment `integral of (x cross n) dA` over the facet,
    /// m^3, oriented out of `owner`.
    ///
    /// This exists so that a rigid rotation is EXACTLY divergence free on
    /// this mesh rather than nearly so. Point-sampling `(omega x x) . n` at
    /// a facet centroid is not the facet's mean value, because the facet is
    /// a bilinear patch and its vertex-average centroid is not its area
    /// centroid; the resulting flux error is small but it is a spurious
    /// source, and MC removes the rotation nullspace on every iteration, so
    /// that error would be injected into the continuity constraint over and
    /// over. Integrating the moment exactly over the two triangles of the
    /// facet instead makes the closed-cell sum vanish identically, because
    /// `surface integral of (x cross n) dA` over any closed polyhedron is
    /// zero by the divergence theorem componentwise. Since
    /// `(omega x x) . n == omega . (x cross n)`, the rotation flux is then
    /// just `omega . moment`.
    pub moment_m3: [f64; 3],
}

#[derive(Debug, thiserror::Error)]
pub enum OpsError {
    #[error("shell geometry error while building topology: {0}")]
    Geometry(#[from] ShellError),
    #[error("field length {got} does not match {want} {what}")]
    Length {
        what: &'static str,
        got: usize,
        want: usize,
    },
    #[error("cell {cell} facet {facet:?} has no recorded face")]
    MissingFace { cell: u32, facet: Facet },
    #[error("{what} produced a non-finite value at index {index}")]
    NonFinite { what: &'static str, index: usize },
    #[error("rigid-rotation Gram matrix is singular (determinant {det:e}); the face normals do not span rotations")]
    SingularRotationGram { det: f64 },
    #[error("face {face} has normal/centroid alignment {cos:e}; a two-point flux needs a strictly positive one")]
    DegenerateFace { face: usize, cos: f64 },
}

/// The face list and cell-to-face incidence for one shell.
///
/// Memory: six entries per cell in the incidence table plus one [`Face`] per
/// interface. At the ADR 026 preview shape that is roughly 650k faces; at the
/// standard shape roughly 4M. Both are recorded here rather than discovered
/// later on a machine that runs out of memory.
#[derive(Debug, Clone)]
pub struct ShellTopology {
    spec: ShellSpec,
    faces: Vec<Face>,
    /// `cell_faces[c * 6 + t]` is the face index for facet `t` of cell `c`,
    /// in [`Facet::ALL`] order.
    cell_faces: Vec<u32>,
    /// Parallel to `cell_faces`: +1 when the stored normal points out of
    /// this cell, -1 when it points in.
    cell_signs: Vec<i8>,
    interior_face_count: usize,
}

impl ShellTopology {
    /// Build the face list. O(cells) with a small constant; every facet is
    /// visited once from each side and stored once.
    pub fn build(spec: &ShellSpec) -> Result<Self, OpsError> {
        let cells = spec.cells();
        let mut faces: Vec<Face> = Vec::new();
        let mut cell_faces = vec![u32::MAX; cells * 6];
        let mut cell_signs = vec![0i8; cells * 6];

        for c in 0..cells {
            let (face, i, j, k) = spec
                .fijk(c as u64)
                .expect("loop index is within the cell count");
            for (t, facet) in Facet::ALL.into_iter().enumerate() {
                if cell_faces[c * 6 + t] != u32::MAX {
                    continue; // already recorded from the other side
                }
                let area_vector = spec.facet_area_vector_m2(face, i, j, k, facet)?;
                let area = norm(area_vector);
                let centroid = facet_centroid(spec, face, i, j, k, facet)?;
                let here = cell_centroid(spec, face, i, j, k)?;

                let (kind, neighbor) = match facet {
                    Facet::Down if k == 0 => (FaceKind::BoundaryInner, None),
                    Facet::Up if k + 1 == spec.layers() => (FaceKind::BoundaryOuter, None),
                    Facet::Down => (FaceKind::Radial, Some(spec.flat(face, i, j, k - 1) as u32)),
                    Facet::Up => (FaceKind::Radial, Some(spec.flat(face, i, j, k + 1) as u32)),
                    Facet::Lateral(dir) => {
                        let (nf, ni, nj) = spec.lateral_neighbor(face, i, j, dir)?;
                        (FaceKind::Lateral, Some(spec.flat(nf, ni, nj, k) as u32))
                    }
                };

                let distance = match neighbor {
                    Some(n) => {
                        let (nf, ni, nj, nk) = spec
                            .fijk(n as u64)
                            .expect("neighbour index is within the cell count");
                        let there = cell_centroid(spec, nf, ni, nj, nk)?;
                        norm(sub(there, here))
                    }
                    None => norm(sub(centroid, here)),
                };

                let moment = antisymmetric_moment(facet_tensor_moment(
                    spec,
                    face,
                    i,
                    j,
                    k,
                    facet,
                    area_vector,
                )?);
                let idx = faces.len() as u32;
                faces.push(Face {
                    owner: c as u32,
                    neighbor,
                    kind,
                    area_m2: area,
                    area_vector_m2: area_vector,
                    normal: scale(area_vector, 1.0 / area),
                    distance_m: distance,
                    centroid_m: centroid,
                    moment_m3: moment,
                });
                cell_faces[c * 6 + t] = idx;
                cell_signs[c * 6 + t] = 1;

                if let Some(n) = neighbor {
                    let n_usize = n as usize;
                    let back = Self::opposite_facet(spec, n, c as u32, facet)?;
                    let slot = n_usize * 6 + facet_index(back);
                    cell_faces[slot] = idx;
                    cell_signs[slot] = -1;
                }
            }
        }

        for (idx, f) in faces.iter().enumerate() {
            let finite = f.area_m2.is_finite()
                && f.distance_m.is_finite()
                && f.area_vector_m2.iter().all(|v| v.is_finite())
                && f.normal.iter().all(|v| v.is_finite())
                && f.centroid_m.iter().all(|v| v.is_finite())
                && f.moment_m3.iter().all(|v| v.is_finite());
            if !finite || f.area_m2 <= 0.0 || f.distance_m <= 0.0 {
                return Err(OpsError::NonFinite {
                    what: "face geometry",
                    index: idx,
                });
            }
        }

        let interior_face_count = faces.iter().filter(|f| f.neighbor.is_some()).count();
        Ok(Self {
            spec: spec.clone(),
            faces,
            cell_faces,
            cell_signs,
            interior_face_count,
        })
    }

    /// Which facet of `neighbor` looks back at `cell`.
    fn opposite_facet(
        spec: &ShellSpec,
        neighbor: u32,
        cell: u32,
        facet: Facet,
    ) -> Result<Facet, OpsError> {
        match facet {
            Facet::Up => Ok(Facet::Down),
            Facet::Down => Ok(Facet::Up),
            Facet::Lateral(_) => {
                let (nf, ni, nj, _) = spec
                    .fijk(neighbor as u64)
                    .expect("neighbour index is within the cell count");
                let (cf, ci, cj, _) = spec
                    .fijk(cell as u64)
                    .expect("cell index is within the cell count");
                for dir in Lateral::ALL {
                    let t = spec.lateral_neighbor(nf, ni, nj, dir)?;
                    if t == (cf, ci, cj) {
                        return Ok(Facet::Lateral(dir));
                    }
                }
                Err(OpsError::MissingFace {
                    cell: neighbor,
                    facet,
                })
            }
        }
    }

    pub fn spec(&self) -> &ShellSpec {
        &self.spec
    }

    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Faces with a velocity unknown: every face that is not on a boundary.
    pub fn interior_face_count(&self) -> usize {
        self.interior_face_count
    }

    pub fn cell_count(&self) -> usize {
        self.spec.cells()
    }

    /// The six (face index, sign) pairs of a cell, in [`Facet::ALL`] order.
    pub fn cell_faces(&self, cell: u32) -> [(u32, i8); 6] {
        let base = cell as usize * 6;
        let mut out = [(0u32, 0i8); 6];
        for (t, slot) in out.iter_mut().enumerate() {
            *slot = (self.cell_faces[base + t], self.cell_signs[base + t]);
        }
        out
    }

    /// Cell volumes, the diagonal of `M_p`.
    pub fn cell_volumes_m3(&self) -> Result<Vec<f64>, OpsError> {
        let mut v = Vec::with_capacity(self.cell_count());
        for c in 0..self.cell_count() {
            let (f, i, j, k) = self
                .spec
                .fijk(c as u64)
                .expect("loop index is within the cell count");
            v.push(self.spec.cell_volume_m3(f, i, j, k)?);
        }
        Ok(v)
    }

    /// Face mass weights `|A_f| * d_f`, the diagonal of `M_u`. Boundary
    /// faces get their weight too, so the vector is indexable by face id;
    /// operators simply never read those entries.
    pub fn face_masses(&self) -> Vec<f64> {
        self.faces
            .iter()
            .map(|f| f.area_m2 * f.distance_m)
            .collect()
    }

    // -----------------------------------------------------------------
    // The operators

    /// Conservative divergence: `(D u)_c = (1/V_c) sum_f s(c,f) u_f |A_f|`.
    ///
    /// Boundary faces are impermeable and contribute nothing, which is the
    /// free-slip normal condition applied at the operator rather than bolted
    /// on afterwards.
    pub fn divergence(&self, u: &[f64], volumes: &[f64]) -> Result<Vec<f64>, OpsError> {
        self.check_len("faces", u.len(), self.face_count())?;
        self.check_len("cells", volumes.len(), self.cell_count())?;
        let mut out = vec![0.0; self.cell_count()];
        for (c, slot) in out.iter_mut().enumerate() {
            let mut sum = 0.0;
            for (f, s) in self.cell_faces(c as u32) {
                let face = &self.faces[f as usize];
                if face.neighbor.is_none() {
                    continue;
                }
                sum += s as f64 * u[f as usize] * face.area_m2;
            }
            *slot = sum / volumes[c];
        }
        Ok(out)
    }

    /// The raw compatible gradient `G = -M_u^-1 D^T M_p` with
    /// `M_u = diag(|A_f| d_f)`, which reduces to the two-point normal
    /// difference `(p_neighbour - p_owner) / d_f`.
    ///
    /// **INCONSISTENT. Do not build a solver on this** - use
    /// [`Self::gradient_corrected`]. It is kept because the defect is worth
    /// pinning: `tests::raw_two_point_dg_laplacian_diverges_under_refinement`
    /// measures `D G p` growing under refinement in all three series, and the
    /// module docs carry the numbers.
    ///
    /// Derived rather than chosen: see the module docs and the identity
    /// test. Boundary faces get zero, matching the divergence.
    pub fn gradient(&self, p: &[f64]) -> Result<Vec<f64>, OpsError> {
        self.check_len("cells", p.len(), self.cell_count())?;
        let mut out = vec![0.0; self.face_count()];
        for (idx, face) in self.faces.iter().enumerate() {
            let Some(n) = face.neighbor else { continue };
            out[idx] = (p[n as usize] - p[face.owner as usize]) / face.distance_m;
        }
        Ok(out)
    }

    // -----------------------------------------------------------------
    // Exact facet tensor moments, and the polyhedral geometry built on them

    /// The exact facet tensor moments `T_f[a][d] = integral_f n_a x_d dA`,
    /// m^3, one per face, oriented OUT OF the face's owner.
    ///
    /// Recomputed from the mesh rather than stored on [`Face`]: at the ADR
    /// 026 standard shape there are roughly 4M faces, and a stored 3x3 would
    /// add 288 MB for a quantity every consumer wants summed per cell anyway.
    pub fn facet_tensor_moments(&self) -> Result<Vec<[[f64; 3]; 3]>, OpsError> {
        let mut out = vec![[[0.0f64; 3]; 3]; self.face_count()];
        for c in 0..self.cell_count() {
            let (f, i, j, k) = self
                .spec
                .fijk(c as u64)
                .expect("loop index is within the cell count");
            for (t, facet) in Facet::ALL.into_iter().enumerate() {
                let idx = self.cell_faces[c * 6 + t] as usize;
                if self.cell_signs[c * 6 + t] != 1 {
                    continue; // this cell is the neighbour, not the owner
                }
                let area_vector = self.spec.facet_area_vector_m2(f, i, j, k, facet)?;
                out[idx] = facet_tensor_moment(&self.spec, f, i, j, k, facet, area_vector)?;
            }
        }
        Ok(out)
    }

    /// The closed-cell tensor moment `sum_f s(c,f) T_f`, which for any closed
    /// polyhedron equals `V_poly * I` exactly.
    ///
    /// The isotropy is the geometric identity a Green-Gauss gradient needs in
    /// order to annihilate rigid-body rotation, and it is measured rather
    /// than assumed
    /// ([`tests::facet_tensor_moments_are_isotropic_over_a_closed_cell`]).
    pub fn cell_tensor_moments(&self) -> Result<Vec<[[f64; 3]; 3]>, OpsError> {
        let facet = self.facet_tensor_moments()?;
        let mut out = vec![[[0.0f64; 3]; 3]; self.cell_count()];
        for (c, cell) in out.iter_mut().enumerate() {
            for (f, s) in self.cell_faces(c as u32) {
                let t = facet[f as usize];
                let sgn = s as f64;
                for (a, row) in cell.iter_mut().enumerate() {
                    for (d, slot) in row.iter_mut().enumerate() {
                        *slot += sgn * t[a][d];
                    }
                }
            }
        }
        Ok(out)
    }

    /// The POLYHEDRAL cell volume `V_poly = trace(sum_f s(c,f) T_f) / 3`.
    ///
    /// This is the volume a Green-Gauss gradient must divide by, and it is
    /// NOT [`ShellSpec::cell_volume_m3`]: that one is the exact spherical
    /// frustum volume and remains the mass geometry `M_p`, while this one is
    /// the volume of the flat-facetted polyhedron whose area vectors the
    /// gradient actually sums over. They differ at `O(h^2)`
    /// ([`tests::polyhedral_volume_differs_from_the_spherical_frustum_at_second_order`]).
    /// Dividing a Green-Gauss sum by the spherical volume scales the gradient
    /// by `1 + O(h^2)` and would leave a spurious strain rate on a rigid
    /// rotation, which is the exact failure the moment construction exists to
    /// remove.
    pub fn polyhedral_cell_volumes_m3(&self) -> Result<Vec<f64>, OpsError> {
        Ok(self
            .cell_tensor_moments()?
            .iter()
            .map(|t| (t[0][0] + t[1][1] + t[2][2]) / 3.0)
            .collect())
    }

    /// Cell centroids in the global Cartesian frame, m.
    pub fn cell_centroids_m(&self) -> Result<Vec<[f64; 3]>, OpsError> {
        let mut out = Vec::with_capacity(self.cell_count());
        for c in 0..self.cell_count() {
            let (f, i, j, k) = self
                .spec
                .fijk(c as u64)
                .expect("loop index is within the cell count");
            out.push(cell_centroid(&self.spec, f, i, j, k)?);
        }
        Ok(out)
    }

    /// The over-relaxed non-orthogonal split of every face (Jasak): with area
    /// vector `S_f` and owner-to-neighbour centroid vector `d_f`, write
    /// `S = Delta + k` with `Delta = (S.S / (d.S)) d` parallel to `d`. Then
    ///
    /// ```text
    /// S . grad p = |Delta| (p_nb - p_own)/|d|  +  k . (grad p)_face
    /// ```
    ///
    /// and `|Delta| / |S| = 1 / (n_hat . d_hat) = lambda_f`, so the first term
    /// is exactly the two-point difference rescaled by `lambda_f` - which is
    /// absorbed into `M_u` and therefore costs no adjointness at all. The
    /// second term is explicit and never enters the matrix.
    pub fn face_corrections(&self) -> Result<Vec<FaceCorrection>, OpsError> {
        let centroids = self.cell_centroids_m()?;
        let mut out = Vec::with_capacity(self.face_count());
        for (idx, face) in self.faces.iter().enumerate() {
            let Some(n) = face.neighbor else {
                // A boundary face is the `w = 1` case of everything below: it
                // interpolates to the owner centroid, so its interpolation
                // offset is the owner-to-facet-centroid vector. That is not
                // cosmetic - it is what lets the Green-Gauss gradient close
                // its boundary facet to second order instead of assuming a
                // zero normal gradient there.
                out.push(FaceCorrection {
                    lambda: 1.0,
                    tangential_m2: [0.0; 3],
                    centroid_dir: [0.0; 3],
                    owner_weight: 1.0,
                    interp_offset_m: sub(face.centroid_m, centroids[face.owner as usize]),
                    skew_m: [0.0; 3],
                    non_orthogonality_rad: 0.0,
                });
                continue;
            };
            let x_own = centroids[face.owner as usize];
            let x_nb = centroids[n as usize];
            let d = sub(x_nb, x_own);
            let dhat = scale(d, 1.0 / norm(d));
            let cos = dot(dhat, face.normal);
            if !cos.is_finite() || cos <= 0.0 {
                return Err(OpsError::DegenerateFace { face: idx, cos });
            }
            let s = face.area_vector_m2;
            let delta = scale(d, dot(s, s) / dot(d, s));
            let own = norm(sub(face.centroid_m, x_own));
            let nb = norm(sub(x_nb, face.centroid_m));
            let midpoint = scale(
                [x_own[0] + x_nb[0], x_own[1] + x_nb[1], x_own[2] + x_nb[2]],
                0.5,
            );
            let w = nb / (own + nb);
            let interp_point = [
                w * x_own[0] + (1.0 - w) * x_nb[0],
                w * x_own[1] + (1.0 - w) * x_nb[1],
                w * x_own[2] + (1.0 - w) * x_nb[2],
            ];
            out.push(FaceCorrection {
                lambda: 1.0 / cos,
                tangential_m2: sub(s, delta),
                centroid_dir: dhat,
                owner_weight: w,
                interp_offset_m: sub(face.centroid_m, interp_point),
                skew_m: sub(face.centroid_m, midpoint),
                non_orthogonality_rad: cos.min(1.0).acos(),
            });
        }
        Ok(out)
    }

    /// The rescaled face mass weights `|A_f| d_f / lambda_f`, the diagonal of
    /// the `M_u` that makes [`Self::gradient_implicit`] the exact weighted
    /// negative adjoint of [`Self::divergence`].
    ///
    /// `M_u` is free: any positive diagonal gives an exactly adjoint pair, so
    /// the over-relaxation is bought for nothing.
    /// [`tests::gradient_is_the_weighted_negative_adjoint_of_divergence`]
    /// checks BOTH pairs as explicit matrices.
    pub fn face_masses_over_relaxed(&self, corr: &[FaceCorrection]) -> Result<Vec<f64>, OpsError> {
        self.check_len("faces", corr.len(), self.face_count())?;
        Ok(self
            .faces
            .iter()
            .zip(corr.iter())
            .map(|(f, c)| f.area_m2 * f.distance_m / c.lambda)
            .collect())
    }

    /// The IMPLICIT half of the corrected gradient:
    /// `G_f = lambda_f (p_nb - p_own) / d_f`.
    ///
    /// Exactly `-M_u^-1 D^T M_p` with
    /// `M_u = diag(|A_f| d_f / lambda_f)` ([`Self::face_masses_over_relaxed`]),
    /// so it is as compatible as the raw two-point gradient and no more
    /// expensive. It is NOT consistent on its own: `lambda` is a rescale, and
    /// consistency comes only from adding
    /// [`Self::gradient_nonorthogonal_correction`].
    pub fn gradient_implicit(
        &self,
        p: &[f64],
        corr: &[FaceCorrection],
    ) -> Result<Vec<f64>, OpsError> {
        self.check_len("cells", p.len(), self.cell_count())?;
        self.check_len("faces", corr.len(), self.face_count())?;
        let mut out = vec![0.0; self.face_count()];
        for (idx, face) in self.faces.iter().enumerate() {
            let Some(n) = face.neighbor else { continue };
            out[idx] =
                corr[idx].lambda * (p[n as usize] - p[face.owner as usize]) / face.distance_m;
        }
        Ok(out)
    }

    /// Cell-centred Green-Gauss gradient of a scalar:
    /// `(grad p)_c = (1 / V_poly_c) sum_f s(c,f) p_f A_f`.
    ///
    /// `p_f` is the linear interpolation of the two cell values to the facet
    /// centroid ([`FaceCorrection::owner_weight`]). A boundary facet is the
    /// `w = 1` case: without a seed it takes the owner value, which is a zero
    /// normal-gradient assumption and only zeroth order; with a seed it
    /// extrapolates `p_own + g_own . (x_facet - x_own)`, which is second
    /// order. The difference is not confined to the boundary layer, because
    /// the skewness term of
    /// [`FaceGradientRule::NonOrthogonalAndSkew`] differences neighbouring
    /// cell gradients and divides by `d_f`, so an `O(1)` gradient error in the
    /// boundary row lands as an `O(1)` Laplacian error one row inside it.
    ///
    /// `skew_seed` is what makes this operator usable. Without it the face
    /// value lands on the centroid SEGMENT rather than at the facet centroid,
    /// and on this grid the gap does not close under refinement at the seams
    /// (measured 0.236 -> 0.275 of `d_f`), so the resulting `O(h)` face bias
    /// is multiplied by `|A_f| / V_c ~ 1/h` and the gradient does not
    /// converge at all. With a seed gradient the face value is shifted to the
    /// facet centroid, `p_f += g_f . interp_offset_f`, and the gradient is
    /// second order. Pass [`Self::cell_gradient_least_squares`] as the seed.
    ///
    /// The volume is the POLYHEDRAL one. See
    /// [`Self::polyhedral_cell_volumes_m3`].
    pub fn cell_gradient_green_gauss(
        &self,
        p: &[f64],
        v_poly: &[f64],
        corr: &[FaceCorrection],
        skew_seed: Option<&[[f64; 3]]>,
    ) -> Result<Vec<[f64; 3]>, OpsError> {
        self.check_len("cells", p.len(), self.cell_count())?;
        self.check_len("cells", v_poly.len(), self.cell_count())?;
        self.check_len("faces", corr.len(), self.face_count())?;
        if let Some(seed) = skew_seed {
            self.check_len("cells", seed.len(), self.cell_count())?;
        }
        let mut out = vec![[0.0f64; 3]; self.cell_count()];
        for (c, slot) in out.iter_mut().enumerate() {
            let mut acc = [0.0f64; 3];
            for (f, s) in self.cell_faces(c as u32) {
                let face = &self.faces[f as usize];
                let fc = &corr[f as usize];
                let w = fc.owner_weight;
                let other = face.neighbor.unwrap_or(face.owner) as usize;
                let mut value = w * p[face.owner as usize] + (1.0 - w) * p[other];
                if let Some(seed) = skew_seed {
                    let go = seed[face.owner as usize];
                    let gn = seed[other];
                    let gf = [
                        w * go[0] + (1.0 - w) * gn[0],
                        w * go[1] + (1.0 - w) * gn[1],
                        w * go[2] + (1.0 - w) * gn[2],
                    ];
                    value += dot(gf, fc.interp_offset_m);
                }
                let contribution = scale(face.area_vector_m2, s as f64 * value);
                acc = [
                    acc[0] + contribution[0],
                    acc[1] + contribution[1],
                    acc[2] + contribution[2],
                ];
            }
            *slot = scale(acc, 1.0 / v_poly[c]);
        }
        Ok(out)
    }

    /// Cell-centred gradient of a scalar by weighted least squares over the
    /// face neighbours: minimise `sum_nb w_nb (g . (x_nb - x_c) - (p_nb - p_c))^2`
    /// with `w_nb = 1 / |x_nb - x_c|^2`.
    ///
    /// The scalar analogue of [`Self::reconstruct_cell_velocity`], and the
    /// gradient the deferred correction should actually use. It is EXACT for
    /// an affine field on any mesh and for any positive weights, which is the
    /// property that matters here and is not a function of the weight choice.
    ///
    /// [`Self::cell_gradient_green_gauss`] is the alternative and it does NOT
    /// converge on this grid: its face value is a linear interpolation
    /// between two cell centres, and the facet centroid sits permanently
    /// about a quarter of the centroid distance off that line at the seams
    /// (measured skewness rises to 0.28 and does not shrink with h), so the
    /// face value carries an `O(h)` bias that the `|A_f| / V_c ~ 1/h` factor
    /// turns into an `O(1)` gradient error. Measured in
    /// [`tests::cell_gradient_convergence_separates_least_squares_from_green_gauss`].
    ///
    /// Boundary cells need no special case: they simply have five neighbours
    /// instead of six, which still over-determines a 3-vector.
    pub fn cell_gradient_least_squares(
        &self,
        p: &[f64],
        centroids: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>, OpsError> {
        self.check_len("cells", p.len(), self.cell_count())?;
        self.check_len("cells", centroids.len(), self.cell_count())?;
        let mut out = Vec::with_capacity(self.cell_count());
        for c in 0..self.cell_count() {
            let mut m = [[0.0f64; 3]; 3];
            let mut rhs = [0.0f64; 3];
            for (f, _) in self.cell_faces(c as u32) {
                let face = &self.faces[f as usize];
                let Some(n) = face.neighbor else { continue };
                let other = if face.owner as usize == c {
                    n as usize
                } else {
                    face.owner as usize
                };
                let dx = sub(centroids[other], centroids[c]);
                let w = 1.0 / dot(dx, dx);
                let dp = p[other] - p[c];
                for (a, row) in m.iter_mut().enumerate() {
                    rhs[a] += w * dx[a] * dp;
                    for (b, slot) in row.iter_mut().enumerate() {
                        *slot += w * dx[a] * dx[b];
                    }
                }
            }
            out.push(solve3(m, rhs)?);
        }
        Ok(out)
    }

    /// Green-Gauss cell gradient TENSOR from per-facet integrals of
    /// `n (x) phi`: `(grad phi)_c[a][d] = (1/V_poly_c) sum_f s(c,f) I_f[a][d]`
    /// with `I_f[a][d] = integral_f n_a phi_d dA` oriented out of the owner.
    ///
    /// Taking the facet integrals as the input rather than face values is the
    /// whole point: supplied EXACTLY (which
    /// [`Self::facet_tensor_moments`] makes possible for any affine field)
    /// this operator is exact for a linear field and annihilates rigid-body
    /// rotation to machine zero; supplied as the point sample
    /// `A_f n_hat_f (x) phi(x_f)` it leaks a spurious strain rate at
    /// `O(h^2)`. Both are measured in
    /// [`tests::green_gauss_gradient_annihilates_rigid_rotation_only_with_exact_moments`].
    pub fn green_gauss_gradient(
        &self,
        facet_integrals: &[[[f64; 3]; 3]],
        v_poly: &[f64],
    ) -> Result<Vec<[[f64; 3]; 3]>, OpsError> {
        self.check_len("faces", facet_integrals.len(), self.face_count())?;
        self.check_len("cells", v_poly.len(), self.cell_count())?;
        let mut out = vec![[[0.0f64; 3]; 3]; self.cell_count()];
        for (c, cell) in out.iter_mut().enumerate() {
            for (f, s) in self.cell_faces(c as u32) {
                let i = facet_integrals[f as usize];
                let sgn = s as f64;
                for (a, row) in cell.iter_mut().enumerate() {
                    for (d, slot) in row.iter_mut().enumerate() {
                        *slot += sgn * i[a][d];
                    }
                }
            }
            let inv = 1.0 / v_poly[c];
            for row in cell.iter_mut() {
                for slot in row.iter_mut() {
                    *slot *= inv;
                }
            }
        }
        Ok(out)
    }

    /// The EXPLICIT half of the corrected gradient:
    /// `(k_f / |A_f|) . (grad p)_face`, the deferred non-orthogonal
    /// correction.
    ///
    /// It is a right-hand-side source, not a matrix entry: it never passes
    /// through `G`, so it cannot perturb the adjoint identity. That is the
    /// whole reason F1 is affordable.
    pub fn gradient_nonorthogonal_correction(
        &self,
        p: &[f64],
        cell_gradient: &[[f64; 3]],
        corr: &[FaceCorrection],
        rule: FaceGradientRule,
    ) -> Result<Vec<f64>, OpsError> {
        self.check_len("cells", p.len(), self.cell_count())?;
        self.check_len("cells", cell_gradient.len(), self.cell_count())?;
        self.check_len("faces", corr.len(), self.face_count())?;
        let mut out = vec![0.0; self.face_count()];
        for (idx, face) in self.faces.iter().enumerate() {
            let Some(n) = face.neighbor else { continue };
            let c = &corr[idx];
            let go = cell_gradient[face.owner as usize];
            let gn = cell_gradient[n as usize];
            let w = c.owner_weight;
            let mut gf = [
                w * go[0] + (1.0 - w) * gn[0],
                w * go[1] + (1.0 - w) * gn[1],
                w * go[2] + (1.0 - w) * gn[2],
            ];
            if rule == FaceGradientRule::NonOrthogonalCentroidLine {
                let two_point = (p[n as usize] - p[face.owner as usize]) / face.distance_m;
                let along = dot(gf, c.centroid_dir);
                let fix = scale(c.centroid_dir, two_point - along);
                gf = [gf[0] + fix[0], gf[1] + fix[1], gf[2] + fix[2]];
            }
            let skew = if rule == FaceGradientRule::NonOrthogonalAndSkew {
                // (p_nb + g_nb . r) - (p_own + g_own . r), over d_f, times
                // lambda: the two-point difference re-centred on the facet
                // centroid. Explicit, so it never touches the matrix.
                c.lambda * dot(sub(gn, go), c.skew_m) / face.distance_m
            } else {
                0.0
            };
            out[idx] = dot(c.tangential_m2, gf) / face.area_m2 + skew;
        }
        Ok(out)
    }

    /// The corrected face-normal pressure gradient: implicit rescaled
    /// two-point difference plus explicit non-orthogonal correction.
    ///
    /// This is the operator a solver should apply as `G`; in a real solve the
    /// two halves live on opposite sides of the equals sign (the first in the
    /// matrix, the second in the residual), and this method exists so the
    /// consistency of their SUM can be measured directly.
    ///
    /// The cell gradient is an ARGUMENT, like `volumes` is to
    /// [`Self::divergence`], because which one you pass changes the answer by
    /// an order of magnitude: use [`Self::cell_gradient_least_squares`].
    /// [`Self::cell_gradient_green_gauss`] is measured and does not converge
    /// on this mesh.
    pub fn gradient_corrected(
        &self,
        p: &[f64],
        cell_gradient: &[[f64; 3]],
        corr: &[FaceCorrection],
        rule: FaceGradientRule,
    ) -> Result<Vec<f64>, OpsError> {
        let implicit = self.gradient_implicit(p, corr)?;
        let explicit = self.gradient_nonorthogonal_correction(p, cell_gradient, corr, rule)?;
        Ok(implicit
            .iter()
            .zip(explicit.iter())
            .map(|(a, b)| a + b)
            .collect())
    }

    // -----------------------------------------------------------------
    // Nullspaces

    /// Remove the constant pressure mode, volume-weighted so the projection
    /// is orthogonal in the `M_p` inner product.
    pub fn remove_mean_pressure(&self, p: &mut [f64], volumes: &[f64]) -> Result<f64, OpsError> {
        self.check_len("cells", p.len(), self.cell_count())?;
        self.check_len("cells", volumes.len(), self.cell_count())?;
        let total: f64 = volumes.iter().sum();
        let mean: f64 = p
            .iter()
            .zip(volumes.iter())
            .map(|(pi, v)| pi * v)
            .sum::<f64>()
            / total;
        for pi in p.iter_mut() {
            *pi -= mean;
        }
        Ok(mean)
    }

    /// The face-normal velocity field of a rigid rotation with angular
    /// velocity `omega` (rad/s), as the facet's exact mean normal component
    /// `omega . moment / area` rather than a centroid point sample. See
    /// [`Face::moment_m3`] for why the distinction is load bearing.
    pub fn rigid_rotation_field(&self, omega: [f64; 3]) -> Vec<f64> {
        self.faces
            .iter()
            .map(|f| {
                if f.neighbor.is_none() {
                    0.0
                } else {
                    dot(omega, f.moment_m3) / f.area_m2
                }
            })
            .collect()
    }

    /// The three rigid-rotation basis fields evaluated on one face, matching
    /// [`Self::rigid_rotation_field`] exactly so the fit and the subtraction
    /// live in the same space.
    fn rotation_basis(face: &Face) -> [f64; 3] {
        [
            face.moment_m3[0] / face.area_m2,
            face.moment_m3[1] / face.area_m2,
            face.moment_m3[2] / face.area_m2,
        ]
    }

    /// The angular velocity of the rigid rotation that best fits `u` in the
    /// `M_u` inner product.
    ///
    /// This is a 3x3 weighted least-squares solve. It exists so that the
    /// rotation component can be REPORTED as well as removed: a solver whose
    /// spurious spin is growing is telling us something, and silently
    /// projecting it away every iteration would hide it.
    ///
    /// **`M_u` here is [`Self::face_masses_over_relaxed`], not
    /// [`Self::face_masses`].** Those two were the same vector until the
    /// corrected gradient moved onto the over-relaxed mass, and the fit has
    /// to move with it: the projection `u - omega . phi` is orthogonal in
    /// whichever inner product weights this solve, and if that is not the
    /// inner product the operators live in, the removed component is not the
    /// component the operator sees. What is left behind is a small persistent
    /// spin with no visible cause - measured at 7.6e-4 of the field norm on a
    /// random field before this change and 1.9e-16 after, in
    /// `pf_geodynamics::mantle::stokes`'s
    /// `the_rigid_rotation_projection_is_orthogonal_in_the_operators_inner_product`.
    /// The two rotation tests below keep their verdicts and move only in
    /// their last digits (recovery error 3.34e-15 -> 2.92e-15, radial-flow
    /// spurious spin 3.14e-25 -> 3.03e-25): they fit fields that a rotation
    /// reproduces exactly, and an exact fit is weight independent.
    pub fn rigid_rotation_component(&self, u: &[f64]) -> Result<[f64; 3], OpsError> {
        self.check_len("faces", u.len(), self.face_count())?;
        let masses = self.face_masses_over_relaxed(&self.face_corrections()?)?;
        let mut gram = [[0.0f64; 3]; 3];
        let mut rhs = [0.0f64; 3];
        for (idx, face) in self.faces.iter().enumerate() {
            if face.neighbor.is_none() {
                continue;
            }
            let w = masses[idx];
            let phi = Self::rotation_basis(face);
            for a in 0..3 {
                rhs[a] += w * phi[a] * u[idx];
                for b in 0..3 {
                    gram[a][b] += w * phi[a] * phi[b];
                }
            }
        }
        solve3(gram, rhs)
    }

    /// Project the whole-mantle rigid rotation out of a face-normal velocity
    /// field, putting the result in a non-rotating mantle frame. Returns the
    /// angular velocity that was removed.
    pub fn remove_rigid_rotation(&self, u: &mut [f64]) -> Result<[f64; 3], OpsError> {
        let omega = self.rigid_rotation_component(u)?;
        let field = self.rigid_rotation_field(omega);
        for (ui, fi) in u.iter_mut().zip(field.iter()) {
            *ui -= fi;
        }
        Ok(omega)
    }

    // -----------------------------------------------------------------
    // Derived quantities

    /// Reconstruct a cell-centred Cartesian velocity from face normals by
    /// weighted least squares: minimise `sum_f w_f (v . n_f - u_f)^2`.
    ///
    /// Exact for a constant field, which the tests check. This is a DERIVED
    /// diagnostic and interpolation quantity; it is never an unknown, and no
    /// conserved quantity is computed from it.
    pub fn reconstruct_cell_velocity(&self, u: &[f64]) -> Result<Vec<[f64; 3]>, OpsError> {
        self.check_len("faces", u.len(), self.face_count())?;
        let mut out = Vec::with_capacity(self.cell_count());
        for c in 0..self.cell_count() {
            let mut m = [[0.0f64; 3]; 3];
            let mut rhs = [0.0f64; 3];
            for (f, s) in self.cell_faces(c as u32) {
                let face = &self.faces[f as usize];
                if face.neighbor.is_none() {
                    continue;
                }
                // Orient both the normal and the value into this cell's
                // outward convention so the fit is sign-consistent.
                let sgn = s as f64;
                let n = scale(face.normal, sgn);
                let value = sgn * u[f as usize];
                let w = face.area_m2;
                for a in 0..3 {
                    rhs[a] += w * n[a] * value;
                    for b in 0..3 {
                        m[a][b] += w * n[a] * n[b];
                    }
                }
            }
            out.push(solve3(m, rhs)?);
        }
        Ok(out)
    }

    /// Fail loud on a non-finite entry, naming the index. Exposed so a
    /// solver can check its own iterates at the seams it chooses rather than
    /// paying for a scan inside every operator application.
    pub fn check_finite(what: &'static str, values: &[f64]) -> Result<(), OpsError> {
        for (index, v) in values.iter().enumerate() {
            if !v.is_finite() {
                return Err(OpsError::NonFinite { what, index });
            }
        }
        Ok(())
    }

    fn check_len(&self, what: &'static str, got: usize, want: usize) -> Result<(), OpsError> {
        if got != want {
            return Err(OpsError::Length { what, got, want });
        }
        Ok(())
    }
}

fn facet_index(facet: Facet) -> usize {
    Facet::ALL
        .iter()
        .position(|f| *f == facet)
        .expect("Facet::ALL is exhaustive")
}

fn cell_centroid(
    spec: &ShellSpec,
    face: usize,
    i: usize,
    j: usize,
    k: usize,
) -> Result<[f64; 3], ShellError> {
    let d = spec.cell_direction(face, i, j)?;
    let r = spec.center_radius_m(k)?;
    Ok(scale(d, r))
}

/// The exact facet tensor moment `T[a][d] = integral of n_a x_d dA`, m^3,
/// by splitting the facet into the same two triangles whose vector areas sum
/// to the facet's stored area vector.
///
/// On each planar triangle `n` is constant and `x` is linear, so
/// `integral n_a x_d dA = area_vector_a * centroid_d` EXACTLY; summing the
/// two triangles and applying the same winding-vs-stored-orientation sign
/// flip gives the facet's exact tensor. Two facts fall out of this and
/// nothing else does:
///
/// 1. `sum_f s(c,f) T_f = V_poly_c * I` exactly for a closed cell, so a
///    Green-Gauss gradient built on it annihilates rigid rotation to machine
///    zero instead of to `O(h^2)`
///    ([`tests::facet_tensor_moments_are_isotropic_over_a_closed_cell`]).
///    The point-sampled substitute `A_f n_hat_f (x) x_f` is isotropic only to
///    `O(h^2)` and that error is a spurious strain rate.
/// 2. Its antisymmetric part is the rotation moment
///    `integral of (x cross n) dA` the divergence-free rotation basis needs,
///    so there is ONE construction here, not two that could drift apart.
fn facet_tensor_moment(
    spec: &ShellSpec,
    face: usize,
    i: usize,
    j: usize,
    k: usize,
    facet: Facet,
    area_vector: [f64; 3],
) -> Result<[[f64; 3]; 3], ShellError> {
    let q = facet_quad(spec, face, i, j, k, facet)?;
    let tri = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
        let av = scale(cross(sub(b, a), sub(c, a)), 0.5);
        let centroid = scale(
            [a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]],
            1.0 / 3.0,
        );
        (av, centroid)
    };
    let (a0, c0) = tri(q[0], q[1], q[2]);
    let (a1, c1) = tri(q[0], q[2], q[3]);
    let mut t = [[0.0f64; 3]; 3];
    for (a, row) in t.iter_mut().enumerate() {
        for (d, slot) in row.iter_mut().enumerate() {
            *slot = a0[a] * c0[d] + a1[a] * c1[d];
        }
    }
    // The stored area vector was oriented outward, possibly against the
    // vertex winding; the moment must follow the same orientation or the
    // closed-cell sum will not cancel.
    let winding = [a0[0] + a1[0], a0[1] + a1[1], a0[2] + a1[2]];
    if dot(winding, area_vector) < 0.0 {
        for row in t.iter_mut() {
            for slot in row.iter_mut() {
                *slot = -*slot;
            }
        }
    }
    Ok(t)
}

/// The antisymmetric part of a facet tensor moment, as the vector
/// `integral of (x cross n) dA`: `m_i = eps_ijk T[k][j]`.
fn antisymmetric_moment(t: [[f64; 3]; 3]) -> [f64; 3] {
    [t[2][1] - t[1][2], t[0][2] - t[2][0], t[1][0] - t[0][1]]
}

fn facet_quad(
    spec: &ShellSpec,
    face: usize,
    i: usize,
    j: usize,
    k: usize,
    facet: Facet,
) -> Result<[[f64; 3]; 4], ShellError> {
    let v = spec.cell_vertices(face, i, j, k)?;
    let pick: [usize; 4] = match facet {
        Facet::Lateral(Lateral::AlphaMinus) => [0, 3, 7, 4],
        Facet::Lateral(Lateral::AlphaPlus) => [1, 2, 6, 5],
        Facet::Lateral(Lateral::BetaMinus) => [0, 1, 5, 4],
        Facet::Lateral(Lateral::BetaPlus) => [3, 2, 6, 7],
        Facet::Down => [0, 1, 2, 3],
        Facet::Up => [4, 5, 6, 7],
    };
    Ok([v[pick[0]], v[pick[1]], v[pick[2]], v[pick[3]]])
}

fn facet_centroid(
    spec: &ShellSpec,
    face: usize,
    i: usize,
    j: usize,
    k: usize,
    facet: Facet,
) -> Result<[f64; 3], ShellError> {
    let q = facet_quad(spec, face, i, j, k, facet)?;
    let mut acc = [0.0f64; 3];
    for v in q {
        acc = [acc[0] + v[0], acc[1] + v[1], acc[2] + v[2]];
    }
    Ok(scale(acc, 0.25))
}

/// Solve a 3x3 system by Cramer's rule. Fails loud on a singular matrix
/// rather than returning a plausible-looking wrong answer.
fn solve3(m: [[f64; 3]; 3], rhs: [f64; 3]) -> Result<[f64; 3], OpsError> {
    let det = det3(m);
    let scale_ref = m
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0f64, |a, b| a.max(b.abs()));
    if det.abs() <= 1e-12 * scale_ref.powi(3).max(f64::MIN_POSITIVE) {
        return Err(OpsError::SingularRotationGram { det });
    }
    let mut out = [0.0f64; 3];
    for (col, slot) in out.iter_mut().enumerate() {
        let mut mc = m;
        for row in 0..3 {
            mc[row][col] = rhs[row];
        }
        *slot = det3(mc) / det;
    }
    Ok(out)
}

fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A gradient operator under test, as a closure over one topology.
    type GradFn<'a> = &'a dyn Fn(&[f64]) -> Vec<f64>;

    const R_CMB: f64 = 3.480e6;
    const R_TOP: f64 = 6.271e6;

    fn mini(n: usize, layers: usize) -> ShellTopology {
        let spec = ShellSpec::new(n, layers, R_CMB, R_TOP).expect("valid shell");
        ShellTopology::build(&spec).expect("topology builds")
    }

    #[test]
    fn face_counts_match_the_analytic_expectation() {
        // Lateral faces: 4 per cell shared by 2 -> 12 n^2 layers.
        // Radial interior: 6 n^2 (layers - 1). Boundary: 2 * 6 n^2.
        for (n, layers) in [(2usize, 2usize), (4, 3), (6, 4)] {
            let t = mini(n, layers);
            let n2 = n * n;
            let lateral = 12 * n2 * layers;
            let radial = 6 * n2 * (layers - 1);
            let boundary = 12 * n2;
            assert_eq!(t.face_count(), lateral + radial + boundary, "n={n}");
            assert_eq!(t.interior_face_count(), lateral + radial, "n={n}");
            assert_eq!(
                t.faces()
                    .iter()
                    .filter(|f| f.kind == FaceKind::BoundaryInner)
                    .count(),
                6 * n2
            );
            assert_eq!(
                t.faces()
                    .iter()
                    .filter(|f| f.kind == FaceKind::BoundaryOuter)
                    .count(),
                6 * n2
            );
        }
    }

    #[test]
    fn every_cell_facet_is_recorded_exactly_once_from_each_side() {
        let t = mini(4, 3);
        let mut seen = vec![0usize; t.face_count()];
        for c in 0..t.cell_count() {
            for (f, s) in t.cell_faces(c as u32) {
                assert!(s == 1 || s == -1, "cell {c} facet has no sign");
                seen[f as usize] += 1;
                let face = &t.faces()[f as usize];
                if s == 1 {
                    assert_eq!(face.owner, c as u32);
                } else {
                    assert_eq!(face.neighbor, Some(c as u32));
                }
            }
        }
        for (idx, count) in seen.iter().enumerate() {
            let want = if t.faces()[idx].neighbor.is_some() {
                2
            } else {
                1
            };
            assert_eq!(*count, want, "face {idx} referenced {count} times");
        }
    }

    #[test]
    fn divergence_annihilates_a_constant_vector_field() {
        // The payoff of the closed-cell area-vector invariant: a uniform
        // translation has exactly zero discrete divergence everywhere, so no
        // spurious source or sink can be manufactured by the grid.
        let t = mini(6, 3);
        let volumes = t.cell_volumes_m3().expect("volumes");
        let field = [0.7, -0.2, 0.45];
        let u: Vec<f64> = t.faces().iter().map(|f| dot(field, f.normal)).collect();
        let div = t.divergence(&u, &volumes).expect("divergence");

        // Interior cells see only interior faces, so they must vanish to
        // roundoff. Cells touching a boundary legitimately do not, because
        // the impermeable condition removes real flux; check those against
        // the flux they are missing rather than pretending they are zero.
        let scale_ref = field.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        for (c, d) in div.iter().enumerate() {
            let touches_boundary = t
                .cell_faces(c as u32)
                .iter()
                .any(|(f, _)| t.faces()[*f as usize].neighbor.is_none());
            if touches_boundary {
                continue;
            }
            let cell_size = volumes[c].cbrt();
            assert!(
                d.abs() * cell_size / scale_ref < 1e-12,
                "cell {c} divergence {d:e} is not zero"
            );
        }
    }

    #[test]
    fn gradient_is_the_weighted_negative_adjoint_of_divergence() {
        // Assemble D and G as explicit dense matrices on a small shell and
        // check G == -M_u^-1 D^T M_p entrywise. This is the compatibility
        // that makes the pressure and velocity spaces a matched pair; if a
        // future change breaks it, the checkerboard mode comes back and this
        // test is what catches it.
        for (n, layers) in [(2usize, 2usize), (4, 3)] {
            let t = mini(n, layers);
            let nc = t.cell_count();
            let nf = t.face_count();
            let volumes = t.cell_volumes_m3().expect("volumes");
            let masses = t.face_masses();

            let mut d_mat = vec![vec![0.0f64; nf]; nc];
            for (c, row) in d_mat.iter_mut().enumerate() {
                let mut e = vec![0.0; nf];
                for (f_idx, slot) in row.iter_mut().enumerate() {
                    e[f_idx] = 1.0;
                    let div = t.divergence(&e, &volumes).expect("divergence");
                    *slot = div[c];
                    e[f_idx] = 0.0;
                }
            }

            // Both M_u choices are checked, because M_u is FREE: any positive
            // diagonal gives an exactly adjoint pair, and the over-relaxed
            // rescale that makes the gradient consistent is bought for nothing.
            // If a future change breaks either, the checkerboard mode comes back
            // and this test is what catches it.
            let corr = t.face_corrections().expect("corrections");
            let relaxed = t.face_masses_over_relaxed(&corr).expect("relaxed masses");
            let pairs: [(&str, GradFn, &[f64]); 2] = [
                (
                    "raw two-point, M_u = |A| d",
                    &|e: &[f64]| t.gradient(e).expect("gradient"),
                    &masses,
                ),
                (
                    "over-relaxed implicit, M_u = |A| d / lambda",
                    &|e: &[f64]| t.gradient_implicit(e, &corr).expect("implicit gradient"),
                    &relaxed,
                ),
            ];

            for (name, grad, mu) in pairs {
                let mut g_mat = vec![vec![0.0f64; nc]; nf];
                for c in 0..nc {
                    let mut e = vec![0.0; nc];
                    e[c] = 1.0;
                    let g = grad(&e);
                    for (f_idx, row) in g_mat.iter_mut().enumerate() {
                        row[c] = g[f_idx];
                    }
                }

                let mut worst = 0.0f64;
                for (f_idx, face) in t.faces().iter().enumerate() {
                    if face.neighbor.is_none() {
                        continue;
                    }
                    for c in 0..nc {
                        let want = -d_mat[c][f_idx] * volumes[c] / mu[f_idx];
                        let got = g_mat[f_idx][c];
                        let scale_ref = want.abs().max(got.abs()).max(f64::MIN_POSITIVE);
                        worst = worst.max((want - got).abs() / scale_ref);
                    }
                }
                println!("n={n} adjoint identity: {name}: worst rel {worst:e}");
                assert!(
                    worst < 1e-12,
                    "adjoint identity violated for {name}, worst rel {worst:e}"
                );
            }

            // And the rescale is not the identity: lambda actually differs from 1
            // on this mesh, so the check above is not vacuous.
            let worst_lambda = corr
                .iter()
                .zip(t.faces())
                .filter(|(_, f)| f.neighbor.is_some())
                .fold(0.0f64, |a, (c, _)| a.max(c.lambda));
            println!("n={n} max lambda {worst_lambda}");
            assert!(
                worst_lambda - 1.0 > 1e-3,
                "n={n}: max lambda {worst_lambda} - the two M_u choices are too close for this \
             check to have any power (a dropped lambda would have to show up above 1e-12)"
            );
        }
    }

    #[test]
    fn the_pressure_laplacian_is_symmetric_and_negative_semidefinite() {
        // L = D G in the M_p inner product must be symmetric negative
        // semidefinite with the constants as its nullspace. That is the
        // structural property a Schur-complement solve leans on.
        let t = mini(2, 2);
        let nc = t.cell_count();
        let volumes = t.cell_volumes_m3().expect("volumes");

        let mut l = vec![vec![0.0f64; nc]; nc];
        for c in 0..nc {
            let mut e = vec![0.0; nc];
            e[c] = 1.0;
            let g = t.gradient(&e).expect("gradient");
            let d = t.divergence(&g, &volumes).expect("divergence");
            for (r, row) in l.iter_mut().enumerate() {
                row[c] = d[r];
            }
        }

        let mut worst_sym = 0.0f64;
        for a in 0..nc {
            for b in 0..nc {
                let x = l[a][b] * volumes[a];
                let y = l[b][a] * volumes[b];
                let scale_ref = x.abs().max(y.abs()).max(f64::MIN_POSITIVE);
                worst_sym = worst_sym.max((x - y).abs() / scale_ref);
            }
        }
        assert!(worst_sym < 1e-10, "M_p L is not symmetric: {worst_sym:e}");

        // Constants are in the nullspace.
        let ones = vec![1.0; nc];
        let g = t.gradient(&ones).expect("gradient");
        assert!(g.iter().all(|v| v.abs() < 1e-30));

        // Negative semidefinite on a few deterministic non-constant vectors.
        for shift in 0..4usize {
            let x: Vec<f64> = (0..nc)
                .map(|i| ((i * 7 + shift * 13) % 11) as f64 - 5.0)
                .collect();
            let mut quad = 0.0;
            for a in 0..nc {
                for b in 0..nc {
                    quad += x[a] * volumes[a] * l[a][b] * x[b];
                }
            }
            assert!(quad <= 1e-6, "x^T M_p L x = {quad:e} is not <= 0");
            assert!(
                quad < 0.0,
                "a non-constant vector should give a strict sign"
            );
        }
    }

    #[test]
    fn mean_pressure_removal_is_an_idempotent_volume_weighted_projection() {
        let t = mini(4, 2);
        let volumes = t.cell_volumes_m3().expect("volumes");
        let mut p: Vec<f64> = (0..t.cell_count()).map(|i| (i % 17) as f64).collect();
        let removed = t.remove_mean_pressure(&mut p, &volumes).expect("removes");
        assert!(removed.abs() > 0.0);

        let total: f64 = volumes.iter().sum();
        let mean: f64 = p.iter().zip(&volumes).map(|(pi, v)| pi * v).sum::<f64>() / total;
        assert!(mean.abs() < 1e-9, "mean not removed: {mean:e}");

        let again = t.remove_mean_pressure(&mut p, &volumes).expect("removes");
        assert!(
            again.abs() < 1e-9,
            "projection is not idempotent: {again:e}"
        );
    }

    #[test]
    fn rigid_rotation_is_recovered_and_removed_exactly() {
        // A whole-mantle spin is physically null under free-slip impermeable
        // boundaries. It must be recoverable (so it can be reported) and
        // removable (so MF never fits plate motion to a numerical artefact).
        let t = mini(4, 3);
        let omega = [1.3e-15, -4.1e-16, 2.7e-15];
        let mut u = t.rigid_rotation_field(omega);

        let found = t.rigid_rotation_component(&u).expect("fit");
        for a in 0..3 {
            let rel = (found[a] - omega[a]).abs() / omega[a].abs();
            assert!(
                rel < 1e-9,
                "component {a} recovered as {} not {}",
                found[a],
                omega[a]
            );
        }

        let removed = t.remove_rigid_rotation(&mut u).expect("removes");
        let _ = removed;
        let residual = u.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        let original = t
            .rigid_rotation_field(omega)
            .iter()
            .fold(0.0f64, |a, b| a.max(b.abs()));
        assert!(
            residual / original < 1e-9,
            "rotation not removed: residual {residual:e} of {original:e}"
        );

        // Idempotent: a second pass finds nothing.
        let again = t.remove_rigid_rotation(&mut u).expect("removes");
        let mag = norm(again) / norm(omega);
        assert!(mag < 1e-8, "second pass found {mag:e} of the original spin");
    }

    #[test]
    fn facet_moments_close_over_a_cell() {
        // The exact statement behind an exactly divergence-free rotation:
        // the surface integral of (x cross n) dA over any closed polyhedron
        // vanishes, so a cell's six oriented moments must cancel.
        let t = mini(4, 3);
        for c in 0..t.cell_count() {
            let mut sum = [0.0f64; 3];
            let mut scale_ref = 0.0f64;
            for (f, s) in t.cell_faces(c as u32) {
                let face = &t.faces()[f as usize];
                let m = scale(face.moment_m3, s as f64);
                sum = [sum[0] + m[0], sum[1] + m[1], sum[2] + m[2]];
                scale_ref += norm(face.moment_m3);
            }
            let rel = norm(sum) / scale_ref;
            assert!(rel < 1e-12, "cell {c} moments do not close: rel {rel:e}");
        }
    }

    #[test]
    fn rigid_rotation_is_divergence_free_in_the_interior() {
        // A whole-mantle spin moves no mass, so removing it must not disturb
        // continuity. In the interior that now holds EXACTLY, because the
        // rotation flux is omega . moment and a cell's oriented moments
        // cancel identically.
        //
        // At the impermeable boundary it does not hold exactly, and the
        // reason is physical rather than sloppy: a flat facet approximating
        // the spherical surface is not quite perpendicular to the rotation,
        // so the true rotation has a small normal component there which the
        // impermeable condition sets to zero. That is a real modelling
        // residual of a polyhedral sphere, so it is MEASURED and shown to
        // converge rather than hidden inside a loose global tolerance. MC
        // will want this number.
        let omega = [2.0e-15, 1.0e-15, -3.0e-15];
        let mut boundary_leak: Vec<f64> = Vec::new();
        for n in [4usize, 8, 16] {
            let t = mini(n, 3);
            let volumes = t.cell_volumes_m3().expect("volumes");
            let u = t.rigid_rotation_field(omega);
            let div = t.divergence(&u, &volumes).expect("divergence");
            let scale_ref = u.iter().fold(0.0f64, |a, b| a.max(b.abs()));

            let mut worst_boundary = 0.0f64;
            for (c, d) in div.iter().enumerate() {
                let cell_size = volumes[c].cbrt();
                let dimensionless = d.abs() * cell_size / scale_ref;
                let touches_boundary = t
                    .cell_faces(c as u32)
                    .iter()
                    .any(|(f, _)| t.faces()[*f as usize].neighbor.is_none());
                if touches_boundary {
                    worst_boundary = worst_boundary.max(dimensionless);
                } else {
                    assert!(
                        dimensionless < 1e-12,
                        "n={n} interior cell {c} spin divergence {dimensionless:e} is not zero"
                    );
                }
            }
            boundary_leak.push(worst_boundary);
        }

        for w in boundary_leak.windows(2) {
            let ratio = w[0] / w[1];
            assert!(
                ratio > 3.0,
                "boundary spin leak must converge: {:e} -> {:e} is a ratio of {ratio}",
                w[0],
                w[1]
            );
        }
    }

    #[test]
    fn a_field_with_no_spin_is_left_alone() {
        let t = mini(4, 2);
        // Purely radial flow carries no angular momentum about any axis.
        let u: Vec<f64> = t
            .faces()
            .iter()
            .map(|f| {
                if f.neighbor.is_none() {
                    0.0
                } else {
                    let r = norm(f.centroid_m);
                    dot(scale(f.centroid_m, 1.0 / r), f.normal)
                }
            })
            .collect();
        let omega = t.rigid_rotation_component(&u).expect("fit");
        let reference = t.rigid_rotation_field([1.0, 1.0, 1.0]);
        let ref_scale = reference.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        let induced = t.rigid_rotation_field(omega);
        let mag = induced.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        assert!(
            mag / ref_scale < 1e-9,
            "spurious spin {mag:e} found in radial flow"
        );
    }

    #[test]
    fn cell_velocity_reconstruction_is_exact_for_a_constant_field() {
        let t = mini(6, 3);
        let field = [1.7, -0.6, 2.4];
        let u: Vec<f64> = t.faces().iter().map(|f| dot(field, f.normal)).collect();
        let v = t.reconstruct_cell_velocity(&u).expect("reconstructs");
        let scale_ref = norm(field);
        for (c, vi) in v.iter().enumerate() {
            let touches_boundary = t
                .cell_faces(c as u32)
                .iter()
                .any(|(f, _)| t.faces()[*f as usize].neighbor.is_none());
            if touches_boundary {
                continue;
            }
            let err = norm(sub(*vi, field)) / scale_ref;
            assert!(err < 1e-10, "cell {c} reconstructed {vi:?} not {field:?}");
        }
    }

    #[test]
    fn non_finite_fields_are_reported_with_their_index() {
        assert!(ShellTopology::check_finite("u", &[1.0, 2.0, 3.0]).is_ok());
        let err = ShellTopology::check_finite("u", &[1.0, f64::NAN, 3.0])
            .expect_err("NaN must be refused");
        assert!(format!("{err}").contains("index 1"), "error was: {err}");
    }

    /// Which gradient the Laplacian measurement is built from.
    #[derive(Debug, Clone, Copy, PartialEq)]
    enum GradKind {
        /// The shipped two-point difference. Compatible, inconsistent.
        Raw,
        /// Not an operator: the exact `grad p . n_hat` at the facet centroid.
        /// The control that localises the defect in `G`.
        ExactControl,
        /// Not an operator: the exact `grad p . n_hat` evaluated at the
        /// MIDPOINT of the owner-to-neighbour centroid segment rather than at
        /// the facet centroid. The isolation control for skewness - it is the
        /// point the two-point difference actually differentiates about.
        MidpointControl,
        /// Over-relaxed implicit part plus explicit non-orthogonal
        /// correction, built from the named cell gradient.
        Corrected(CellGradKind, FaceGradientRule),
    }

    #[derive(Debug, Clone, Copy, PartialEq)]
    enum CellGradKind {
        GreenGauss,
        GreenGaussSkewSeeded,
        LeastSquares,
        /// Not an operator: the exact analytic gradient of the test field at
        /// each cell centroid. The control that separates "how good is the
        /// cell gradient" from "what does the deferred correction ceiling
        /// look like with a perfect one".
        ExactAnalytic,
    }

    /// Relative errors of `D G p` against the exact Laplacian of
    /// `p = (x/R)^2`, which is the constant `2/R^2`.
    ///
    /// Radial boundary layers (`k == 0`, `k == layers - 1`) are excluded from
    /// every series: the impermeable condition removes real flux there and
    /// the Green-Gauss cell gradient falls to first order at those cells, so
    /// they are a separate story with a separate fix. "Seam-touching" means a
    /// cell on the border of a cube face; "interior" means the rest. Global
    /// is over the union of the two, so it is their max, and the L2 is
    /// volume-weighted over the same set.
    #[derive(Debug, Clone, Copy)]
    struct LapError {
        interior: f64,
        seam: f64,
        global: f64,
        l2: f64,
    }

    /// The three refinement levels every convergence series in this module
    /// uses, built once: a `mini(32, 18)` shell is 110592 cells and rebuilding
    /// it per variant dominated the suite runtime.
    fn refinements() -> [(ShellTopology, usize, usize); 3] {
        [
            (mini(8, 6), 8, 6),
            (mini(16, 10), 16, 10),
            (mini(32, 18), 32, 18),
        ]
    }

    fn dg_laplacian_error(t: &ShellTopology, n: usize, layers: usize, kind: GradKind) -> LapError {
        let volumes = t.cell_volumes_m3().expect("volumes");
        let r = t.spec().r_top_m();
        let mut p = Vec::with_capacity(t.cell_count());
        for c in 0..t.cell_count() {
            let (f, i, j, k) = t.spec().fijk(c as u64).expect("in range");
            let d = t.spec().cell_direction(f, i, j).expect("in range");
            let rr = t.spec().center_radius_m(k).expect("in range");
            p.push((d[0] * rr / r).powi(2));
        }
        let g = match kind {
            GradKind::Raw => t.gradient(&p).expect("gradient"),
            GradKind::ExactControl => {
                // Exact grad p . n at each face centroid: grad p = 2x/R^2 xhat.
                t.faces()
                    .iter()
                    .map(|f| {
                        if f.neighbor.is_none() {
                            0.0
                        } else {
                            2.0 * f.centroid_m[0] / (r * r) * f.normal[0]
                        }
                    })
                    .collect()
            }
            GradKind::MidpointControl => {
                let centroids = t.cell_centroids_m().expect("centroids");
                t.faces()
                    .iter()
                    .map(|f| match f.neighbor {
                        None => 0.0,
                        Some(nb) => {
                            let x = scale(
                                [
                                    centroids[f.owner as usize][0] + centroids[nb as usize][0],
                                    centroids[f.owner as usize][1] + centroids[nb as usize][1],
                                    centroids[f.owner as usize][2] + centroids[nb as usize][2],
                                ],
                                0.5,
                            );
                            2.0 * x[0] / (r * r) * f.normal[0]
                        }
                    })
                    .collect()
            }
            GradKind::Corrected(cell_kind, rule) => {
                let corr = t.face_corrections().expect("corrections");
                let cell_gradient = match cell_kind {
                    CellGradKind::GreenGauss => {
                        let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
                        t.cell_gradient_green_gauss(&p, &v_poly, &corr, None)
                            .expect("cell gradient")
                    }
                    CellGradKind::GreenGaussSkewSeeded => {
                        let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
                        let centroids = t.cell_centroids_m().expect("centroids");
                        let seed = t
                            .cell_gradient_least_squares(&p, &centroids)
                            .expect("cell gradient");
                        t.cell_gradient_green_gauss(&p, &v_poly, &corr, Some(&seed))
                            .expect("cell gradient")
                    }
                    CellGradKind::LeastSquares => {
                        let centroids = t.cell_centroids_m().expect("centroids");
                        t.cell_gradient_least_squares(&p, &centroids)
                            .expect("cell gradient")
                    }
                    CellGradKind::ExactAnalytic => t
                        .cell_centroids_m()
                        .expect("centroids")
                        .iter()
                        .map(|x| [2.0 * x[0] / (r * r), 0.0, 0.0])
                        .collect(),
                };
                t.gradient_corrected(&p, &cell_gradient, &corr, rule)
                    .expect("corrected gradient")
            }
        };
        let lap = t.divergence(&g, &volumes).expect("divergence");
        let exact = 2.0 / (r * r);
        let last = n - 1;
        let (mut interior, mut seam) = (0.0f64, 0.0f64);
        let (mut sq, mut vol) = (0.0f64, 0.0f64);
        for c in 0..t.cell_count() {
            let (_, i, j, k) = t.spec().fijk(c as u64).expect("in range");
            if k == 0 || k + 1 == layers {
                continue; // radial boundary cells are a separate story
            }
            let err = (lap[c] - exact).abs() / exact;
            sq += err * err * volumes[c];
            vol += volumes[c];
            if i == 0 || i == last || j == 0 || j == last {
                seam = seam.max(err);
            } else {
                interior = interior.max(err);
            }
        }
        LapError {
            interior,
            seam,
            global: interior.max(seam),
            l2: (sq / vol).sqrt(),
        }
    }

    /// `log2` of the error ratio between two consecutive 2x refinements: the
    /// OBSERVED convergence order. Negative means the operator is getting
    /// worse under refinement.
    fn order(coarse: f64, fine: f64) -> f64 {
        (coarse / fine).log2()
    }

    fn report(label: &str, e: [LapError; 3]) {
        println!(
            "{label:<28} interior {:.5} {:.5} {:.5}  orders {:+.2} {:+.2}",
            e[0].interior,
            e[1].interior,
            e[2].interior,
            order(e[0].interior, e[1].interior),
            order(e[1].interior, e[2].interior)
        );
        println!(
            "{label:<28} seam     {:.5} {:.5} {:.5}  orders {:+.2} {:+.2}",
            e[0].seam,
            e[1].seam,
            e[2].seam,
            order(e[0].seam, e[1].seam),
            order(e[1].seam, e[2].seam)
        );
        println!(
            "{label:<28} global   {:.5} {:.5} {:.5}  orders {:+.2} {:+.2}",
            e[0].global,
            e[1].global,
            e[2].global,
            order(e[0].global, e[1].global),
            order(e[1].global, e[2].global)
        );
        println!(
            "{label:<28} L2       {:.5} {:.5} {:.5}  orders {:+.2} {:+.2}",
            e[0].l2,
            e[1].l2,
            e[2].l2,
            order(e[0].l2, e[1].l2),
            order(e[1].l2, e[2].l2)
        );
    }

    fn series_on(shells: &[(ShellTopology, usize, usize); 3], kind: GradKind) -> [LapError; 3] {
        [
            dg_laplacian_error(&shells[0].0, shells[0].1, shells[0].2, kind),
            dg_laplacian_error(&shells[1].0, shells[1].1, shells[1].2, kind),
            dg_laplacian_error(&shells[2].0, shells[2].1, shells[2].2, kind),
        ]
    }

    fn series(kind: GradKind) -> [LapError; 3] {
        series_on(&refinements(), kind)
    }

    #[test]
    fn raw_two_point_dg_laplacian_diverges_under_refinement() {
        // CHARACTERIZATION OF A KNOWN DEFECT, still true and still worth
        // pinning: `ShellTopology::gradient` is kept as the compatible but
        // INCONSISTENT operator, and this is what keeps it labelled.
        //
        // Compatibility is not consistency. D and G are exactly adjoint, and
        // that is worth nothing here: on a non-orthogonal grid a DIAGONAL
        // M_u = diag(|A_f| d_f) is only zeroth-order accurate, because the
        // facet normal and the centroid line point in different directions
        // (Thuburn and Cotter 2012, eq. 5.2). The cubed sphere's in-face
        // non-orthogonality grows toward 30 degrees, so the two-point
        // gradient mis-states the pressure force by roughly 50 percent at
        // the worst facets AT EVERY RESOLUTION.
        let e = series(GradKind::Raw);
        report("raw two-point", e);

        // Interior error does not shrink; it grows toward a fixed fraction.
        assert!(
            e[0].interior > 0.15 && e[1].interior > e[0].interior && e[2].interior > e[1].interior,
            "interior: {} {} {}",
            e[0].interior,
            e[1].interior,
            e[2].interior
        );
        assert!(
            e[2].interior > 0.4,
            "interior error saturates near 46 percent: {}",
            e[2].interior
        );
        // Seam-touching cells are worse and diverge roughly as 1/h.
        assert!(
            e[1].seam > 1.8 * e[0].seam && e[2].seam > 1.8 * e[1].seam,
            "seam: {} {} {}",
            e[0].seam,
            e[1].seam,
            e[2].seam
        );
    }

    #[test]
    fn corrected_dg_laplacian_convergence_measured_in_three_series() {
        // THE MEASUREMENT MC-Q2 IS BUILT ON. Three series plus a
        // volume-weighted L2, reported separately because the build brief
        // requires it and because on this grid the seams are the HEALTHY
        // region for non-orthogonality: cross-seam facets are exactly
        // orthogonal by mirror symmetry, while the cube-face interiors carry
        // the 20-30 degree non-orthogonality. A seam-only report would look
        // at the wrong cells.
        //
        // Every combination is measured and printed, including the ones that
        // are worse. Nothing here is tuned to pass: the assertions below
        // state what was OBSERVED, and the observation is a PARTIAL success -
        // L2 converges at better than first order, both Linf series
        // eventually saturate. Thuburn, Cotter and Dubos 2014 Table 4 reports
        // exactly that shape for a consistent non-diagonal Hodge star on a
        // SMOOTHED cubed sphere (Linf 0.21 -> 0.45 rising, L2 almost first
        // order), so a saturating Linf is the published state of the art
        // here, not a tolerance to loosen.
        let shells = refinements();
        let mut shipped: Option<[LapError; 3]> = None;
        let mut ceiling: Option<[LapError; 3]> = None;
        for cell in [
            CellGradKind::LeastSquares,
            CellGradKind::GreenGauss,
            CellGradKind::GreenGaussSkewSeeded,
            CellGradKind::ExactAnalytic,
        ] {
            for face in [
                FaceGradientRule::NonOrthogonalOnly,
                FaceGradientRule::NonOrthogonalCentroidLine,
                FaceGradientRule::NonOrthogonalAndSkew,
            ] {
                let e = series_on(&shells, GradKind::Corrected(cell, face));
                report(&format!("corrected {cell:?}/{face:?}"), e);
                if face == FaceGradientRule::NonOrthogonalAndSkew {
                    if cell == CellGradKind::GreenGaussSkewSeeded {
                        shipped = Some(e);
                    }
                    if cell == CellGradKind::ExactAnalytic {
                        ceiling = Some(e);
                    }
                }
            }
        }
        let e = shipped.expect("the shipped combination was measured");
        let ceiling = ceiling.expect("the exact-gradient ceiling was measured");

        // 1. THE FORM IS CONSISTENT. Fed a perfect cell gradient, the same
        //    implicit-plus-explicit operator is clean second order in every
        //    series - it reproduces the exact-face-normal-derivative control
        //    to three digits. So what is left below is the accuracy of a
        //    compact cell-gradient reconstruction, not a defect in the
        //    discretization, and that is where any further work goes.
        for (label, a, b, c) in [
            (
                "ceiling interior",
                ceiling[0].interior,
                ceiling[1].interior,
                ceiling[2].interior,
            ),
            ("ceiling l2", ceiling[0].l2, ceiling[1].l2, ceiling[2].l2),
        ] {
            assert!(
                order(a, b) > 1.8 && order(b, c) > 1.8,
                "{label} with an exact cell gradient must be second order: \
                 {a:e} -> {b:e} -> {c:e} (orders {:+.2}, {:+.2})",
                order(a, b),
                order(b, c)
            );
        }

        // 2. The shipped operator converges in the interior series and in L2.
        //    Orders are asserted at a floor and PRINTED at their measured
        //    value; the floor is first order because that is what a compact
        //    Green-Gauss cell gradient delivers in Linf on an irregular mesh,
        //    which is measured separately in
        //    `cell_gradient_convergence_separates_least_squares_from_green_gauss`.
        for (label, a, b, c) in [
            ("interior", e[0].interior, e[1].interior, e[2].interior),
            ("l2", e[0].l2, e[1].l2, e[2].l2),
        ] {
            assert!(
                b < a && c < b,
                "{label} must decrease under refinement: {a:e} -> {b:e} -> {c:e} \
                 (observed orders {:+.2}, {:+.2})",
                order(a, b),
                order(b, c)
            );
            assert!(
                order(b, c) > 0.9,
                "{label} observed order {:+.2} between n=16 and n=32 ({b:e} -> {c:e})",
                order(b, c)
            );
        }

        // 3. The seam series SATURATES. It is pinned at its measured plateau
        //    rather than asserted to converge, because it does not: 0.0588 ->
        //    0.0470 -> 0.0485, i.e. it falls then flattens. The raw operator
        //    it replaces DIVERGES there (1.58 -> 3.79 -> 8.12), so this is a
        //    change of kind, and the number is what MC-Q2's seam pass
        //    criterion should be set from. If it ever converges, this
        //    assertion fires and the finding gets rewritten.
        assert!(
            e[2].seam < 0.06 && e[2].seam < e[0].seam,
            "seam plateau moved: {:e} -> {:e} -> {:e} (orders {:+.2}, {:+.2})",
            e[0].seam,
            e[1].seam,
            e[2].seam,
            order(e[0].seam, e[1].seam),
            order(e[1].seam, e[2].seam)
        );

        // 4. Against the operator it replaces, at the finest resolution
        //    measured, in every series.
        let raw = series_on(&shells, GradKind::Raw);
        for (label, corrected, raw) in [
            ("interior", e[2].interior, raw[2].interior),
            ("seam", e[2].seam, raw[2].seam),
            ("global", e[2].global, raw[2].global),
            ("l2", e[2].l2, raw[2].l2),
        ] {
            println!(
                "n=32 {label:<9} corrected {corrected:.5} vs raw {raw:.5}  ({:.0}x better)",
                raw / corrected
            );
            assert!(
                corrected < raw / 10.0,
                "{label}: corrected {corrected:e} is not decisively better than raw {raw:e}"
            );
        }
    }

    #[test]
    fn the_residual_seam_error_is_skewness_not_non_orthogonality() {
        // ISOLATION CONTROL for the part of the defect F1 does not fix.
        //
        // Replace the gradient with the EXACT face-normal derivative, but
        // evaluate it at the midpoint of the owner-to-neighbour centroid
        // segment instead of at the facet centroid. That is the only
        // difference from `the_gradient_is_the_sole_defect...`'s control,
        // which is clean second order - and it is exactly the approximation
        // the implicit two-point difference makes, because a difference
        // quotient differentiates about the midpoint of its own segment.
        //
        // If skewness is what is left, this control must reproduce the
        // corrected operator's seam behaviour while its interior stays clean:
        // in-face skewness shrinks with h (0.033 -> 0.016 -> ...), seam
        // skewness does not (0.236 -> 0.262 -> 0.275, rising to a limit near
        // 0.28), because a seam facet's centroid sits permanently about a
        // quarter of the centroid distance off the line joining the centroids.
        let e = series(GradKind::MidpointControl);
        report("midpoint control", e);
        assert!(
            e[1].interior < e[0].interior && e[2].interior < e[1].interior,
            "the midpoint control's INTERIOR must stay clean - in-face skewness \
             does shrink: {:e} -> {:e} -> {:e}",
            e[0].interior,
            e[1].interior,
            e[2].interior
        );
        assert!(
            e[2].seam > 5.0 * e[2].interior,
            "the midpoint control must reproduce a seam-localised error even with an \
             exact derivative: seam {:e} vs interior {:e}",
            e[2].seam,
            e[2].interior
        );
    }

    #[test]
    fn the_gradient_is_the_sole_defect_not_the_divergence_or_the_geometry() {
        // The control that makes the diagnosis unambiguous: same D, same cell
        // volumes, same facet area vectors, same everything, with ONLY the
        // gradient replaced by the exact face-normal derivative. If the
        // divergence or the mesh geometry were at fault this would be just as
        // bad. It is instead clean second order, which localises the defect
        // entirely in G and therefore in the diagonal M_u behind it. It is
        // also the accuracy ceiling the corrected operator is aiming at.
        let e = series(GradKind::ExactControl);
        report("exact-normal control", e);
        let c = [e[0].interior, e[1].interior, e[2].interior];
        assert!(c[0] < 0.05, "control should already be small: {}", c[0]);
        assert!(
            c[1] < c[0] * 0.6,
            "control must converge: {} -> {}",
            c[0],
            c[1]
        );
        assert!(
            c[2] < c[1] * 0.7,
            "control must keep converging: {} -> {}",
            c[1],
            c[2]
        );
    }

    #[test]
    fn facet_tensor_moments_are_isotropic_over_a_closed_cell() {
        // The geometric identity a Green-Gauss gradient needs in order to
        // annihilate rigid-body rotation: sum_f s(c,f) T_f must be V * I.
        // Exact facet moments give machine zero; the point sample
        // A_f n_hat_f (x) x_f - the obvious substitute - is isotropic only to
        // O(h^2), and that anisotropy IS a spurious strain rate.
        for (n, layers) in [(8usize, 6usize), (16, 10), (32, 18)] {
            let t = mini(n, layers);
            let exact = t.cell_tensor_moments().expect("tensor moments");

            let mut sampled = vec![[[0.0f64; 3]; 3]; t.cell_count()];
            for (c, cell) in sampled.iter_mut().enumerate() {
                for (f, s) in t.cell_faces(c as u32) {
                    let face = &t.faces()[f as usize];
                    let sgn = s as f64;
                    for (a, row) in cell.iter_mut().enumerate() {
                        for (d, slot) in row.iter_mut().enumerate() {
                            *slot += sgn * face.area_vector_m2[a] * face.centroid_m[d];
                        }
                    }
                }
            }

            let anisotropy = |ts: &[[[f64; 3]; 3]]| {
                ts.iter().fold(0.0f64, |acc, m| {
                    let v = (m[0][0] + m[1][1] + m[2][2]) / 3.0;
                    let mut worst = 0.0f64;
                    for (a, row) in m.iter().enumerate() {
                        for (d, x) in row.iter().enumerate() {
                            let want = if a == d { v } else { 0.0 };
                            worst = worst.max((x - want).abs() / v);
                        }
                    }
                    acc.max(worst)
                })
            };

            let (e, s) = (anisotropy(&exact), anisotropy(&sampled));
            println!("n={n:<3} tensor-moment anisotropy: exact {e:.4e}  point-sample {s:.4e}");
            assert!(
                e < 1e-12,
                "n={n}: exact tensor moments must be isotropic to machine zero, got {e:e}"
            );
            assert!(
                s > 100.0 * e,
                "n={n}: the point sample is supposed to be visibly worse; \
                 exact {e:e} vs sampled {s:e}"
            );
        }
    }

    #[test]
    fn polyhedral_volume_differs_from_the_spherical_frustum_at_second_order() {
        // V_poly = trace(T)/3 is the volume the Green-Gauss gradient divides
        // by; cell_volume_m3 is the exact spherical frustum and stays the
        // mass geometry M_p. They are NOT the same number, and the gap is
        // measured here rather than assumed negligible.
        let mut worst = Vec::new();
        for (n, layers) in [(8usize, 6usize), (16, 10), (32, 18)] {
            let t = mini(n, layers);
            let poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
            let sphere = t.cell_volumes_m3().expect("volumes");
            let w = poly
                .iter()
                .zip(sphere.iter())
                .fold(0.0f64, |a, (p, s)| a.max((p - s).abs() / s));
            println!("n={n:<3} max |V_poly - V_sphere| / V_sphere = {w:.4e}");
            worst.push(w);
        }
        for w in worst.windows(2) {
            let ord = order(w[0], w[1]);
            assert!(
                ord > 1.7,
                "the polyhedral/spherical gap must close at second order: \
                 {:e} -> {:e} is order {ord:+.2}",
                w[0],
                w[1]
            );
        }
        assert!(
            worst[0] > 1e-3,
            "the gap is supposed to be a real O(h^2) difference, not roundoff: {:e}",
            worst[0]
        );
    }

    /// Per-facet exact integrals of `n (x) phi` for the affine vector field
    /// `phi(x) = b + L x`, built from the exact facet tensor moments:
    /// `I_f[a][d] = A_f[a] b_d + sum_e L[d][e] T_f[a][e]`.
    fn affine_facet_integrals(
        t: &ShellTopology,
        b: [f64; 3],
        l: [[f64; 3]; 3],
    ) -> Vec<[[f64; 3]; 3]> {
        let moments = t.facet_tensor_moments().expect("facet tensor moments");
        t.faces()
            .iter()
            .zip(moments.iter())
            .map(|(face, tf)| {
                let mut out = [[0.0f64; 3]; 3];
                for (a, row) in out.iter_mut().enumerate() {
                    for (d, slot) in row.iter_mut().enumerate() {
                        let mut acc = face.area_vector_m2[a] * b[d];
                        for (e, le) in l[d].iter().enumerate() {
                            acc += le * tf[a][e];
                        }
                        *slot = acc;
                    }
                }
                out
            })
            .collect()
    }

    /// The same integrals approximated by a single point sample at the stored
    /// facet centroid: `I_f[a][d] = A_f[a] phi_d(x_f)`.
    fn point_sampled_facet_integrals(
        t: &ShellTopology,
        b: [f64; 3],
        l: [[f64; 3]; 3],
    ) -> Vec<[[f64; 3]; 3]> {
        t.faces()
            .iter()
            .map(|face| {
                let x = face.centroid_m;
                let mut phi = b;
                for (d, p) in phi.iter_mut().enumerate() {
                    *p += dot(l[d], x);
                }
                let mut out = [[0.0f64; 3]; 3];
                for (a, row) in out.iter_mut().enumerate() {
                    for (d, slot) in row.iter_mut().enumerate() {
                        *slot = face.area_vector_m2[a] * phi[d];
                    }
                }
                out
            })
            .collect()
    }

    #[test]
    fn green_gauss_gradient_is_exact_for_a_linear_field() {
        // Given exact facet integrals, the Green-Gauss gradient returns the
        // coefficient matrix of an affine field to machine precision, in
        // EVERY cell including seam and cube-corner cells. That is a direct
        // consequence of sum_f s(c,f) T_f = V_poly I and it is the property
        // that makes V_poly, not the spherical volume, the right divisor.
        let t = mini(8, 6);
        let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
        let b = [3.1, -0.7, 2.2];
        let l = [[0.4, -1.1, 0.3], [0.9, 0.25, -0.6], [-0.2, 0.8, 1.4]];
        let integrals = affine_facet_integrals(&t, b, l);
        let grad = t
            .green_gauss_gradient(&integrals, &v_poly)
            .expect("green-gauss");
        let scale_ref = l
            .iter()
            .flat_map(|r| r.iter())
            .fold(0.0f64, |a, x| a.max(x.abs()));
        let mut worst = 0.0f64;
        for g in &grad {
            for (a, row) in g.iter().enumerate() {
                for (d, x) in row.iter().enumerate() {
                    // green_gauss_gradient returns [a][d] = d(phi_d)/d(x_a),
                    // the transpose of the coefficient matrix L.
                    worst = worst.max((x - l[d][a]).abs() / scale_ref);
                }
            }
        }
        println!("green-gauss linear-field worst rel error {worst:.4e}");
        assert!(worst < 1e-12, "not exact for a linear field: {worst:e}");
    }

    #[test]
    fn green_gauss_gradient_annihilates_rigid_rotation_only_with_exact_moments() {
        // v = omega x x has zero strain rate everywhere. With exact facet
        // tensor moments the discrete strain rate is EXACTLY zero; with the
        // point sample it is not, and the residue is a spurious strain the
        // viscous operator would turn into a spurious stress - which MF would
        // then inherit as real plate motion.
        let omega = [0.7, -1.3, 0.4];
        // v_d = (omega x x)_d = eps_dpq omega_p x_q.
        let l = [
            [0.0, -omega[2], omega[1]],
            [omega[2], 0.0, -omega[0]],
            [-omega[1], omega[0], 0.0],
        ];
        let scale_ref = norm(omega);
        let mut sampled_series = Vec::new();
        for (n, layers) in [(8usize, 6usize), (16, 10), (32, 18)] {
            let t = mini(n, layers);
            let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
            let strain = |integrals: &[[[f64; 3]; 3]]| {
                let grad = t
                    .green_gauss_gradient(integrals, &v_poly)
                    .expect("green-gauss");
                grad.iter().fold(0.0f64, |acc, g| {
                    let mut worst = 0.0f64;
                    for (a, row) in g.iter().enumerate() {
                        for (d, x) in row.iter().enumerate() {
                            worst = worst.max(0.5 * (x + g[d][a]).abs() / scale_ref);
                        }
                    }
                    acc.max(worst)
                })
            };
            let exact = strain(&affine_facet_integrals(&t, [0.0; 3], l));
            let sampled = strain(&point_sampled_facet_integrals(&t, [0.0; 3], l));
            println!(
                "n={n:<3} rigid-rotation spurious strain: exact moments {exact:.4e}  \
                 point sample {sampled:.4e}"
            );
            // 1e-12 is the machine-zero bar here, not a fitted bound: the
            // measured values are 9.3e-15 / 3.1e-14 / 1.2e-13, growing with
            // the length of the summation, and the point-sample alternative
            // is 5.6e-3 / 1.6e-3 / 4.6e-4 - nine orders of magnitude away, so
            // no plausible tightening of this line would change a verdict.
            assert!(
                exact < 1e-12,
                "n={n}: exact moments must annihilate rigid rotation, got {exact:e}"
            );
            assert!(
                sampled > 1e-5,
                "n={n}: the point sample is supposed to leak visibly, got {sampled:e}"
            );
            sampled_series.push(sampled);
        }
        // The leak is O(h^2): real, converging, and still a source term that
        // the rotation projection would inject on every iteration.
        for w in sampled_series.windows(2) {
            let ord = order(w[0], w[1]);
            assert!(
                ord > 1.5,
                "point-sample leak should be second order: {:e} -> {:e} is {ord:+.2}",
                w[0],
                w[1]
            );
        }
    }

    #[test]
    fn cell_gradient_convergence_separates_least_squares_from_green_gauss() {
        // The gradient that feeds the explicit correction, measured against
        // the exact one on p = (x/R)^2, in units of |grad p| so the numbers
        // are relative errors.
        //
        // The least-squares gradient converges. The Green-Gauss gradient with
        // interpolated face values DOES NOT, and that is a characterisation,
        // not a bug to be fixed later: its face value is a linear
        // interpolation between two cell centres, while the facet centroid
        // sits about a quarter of the centroid distance off that line at the
        // seams and stays there under refinement. The resulting O(h) face bias
        // is multiplied by |A_f| / V_c ~ 1/h. This is why
        // `gradient_corrected` takes the cell gradient as an argument.
        let mut lsq = Vec::new();
        let mut gg = Vec::new();
        let mut seeded = Vec::new();
        for (n, layers) in [(8usize, 6usize), (16, 10), (32, 18)] {
            let t = mini(n, layers);
            let r = t.spec().r_top_m();
            let corr = t.face_corrections().expect("corrections");
            let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
            let centroids = t.cell_centroids_m().expect("centroids");
            let p: Vec<f64> = centroids.iter().map(|x| (x[0] / r).powi(2)).collect();
            let lsq_grad = t
                .cell_gradient_least_squares(&p, &centroids)
                .expect("cell gradient");
            let grads = [
                ("least-squares", lsq_grad.clone()),
                (
                    "green-gauss",
                    t.cell_gradient_green_gauss(&p, &v_poly, &corr, None)
                        .expect("cell gradient"),
                ),
                (
                    "gg-skew-seeded",
                    t.cell_gradient_green_gauss(&p, &v_poly, &corr, Some(&lsq_grad))
                        .expect("cell gradient"),
                ),
            ];
            for (label, grad) in grads {
                let mut w = 0.0f64;
                for (c, g) in grad.iter().enumerate() {
                    let (_, _, _, k) = t.spec().fijk(c as u64).expect("in range");
                    if k == 0 || k + 1 == layers {
                        continue; // first-order boundary closure, excluded and declared
                    }
                    let want = [2.0 * centroids[c][0] / (r * r), 0.0, 0.0];
                    // Normalised by max |grad p| = 2 r_top / R^2 = 2 / R.
                    w = w.max(norm(sub(*g, want)) * r / 2.0);
                }
                println!("n={n:<3} cell gradient {label:<14} Linf rel = {w:.4e}");
                match label {
                    "least-squares" => lsq.push(w),
                    "green-gauss" => gg.push(w),
                    _ => seeded.push(w),
                }
            }
        }
        for x in lsq.windows(2) {
            let ord = order(x[0], x[1]);
            assert!(
                ord > 0.85,
                "the least-squares cell gradient converges at FIRST order in Linf: \
                 {:e} -> {:e} is order {ord:+.2}",
                x[0],
                x[1]
            );
        }
        for x in seeded.windows(2) {
            let ord = order(x[0], x[1]);
            assert!(
                ord > 0.95,
                "the skew-seeded Green-Gauss cell gradient - the one the shipped \
                 operator uses - must converge: {:e} -> {:e} is order {ord:+.2}",
                x[0],
                x[1]
            );
        }
        // It converges at FIRST order in Linf, like the least-squares one and
        // for the same reason: a compact reconstruction from cell-centre
        // values on a mesh whose centres are not symmetrically placed leaves
        // an O(h^2) face residue that |A_f| / V_c ~ 1/h turns into O(h).
        // Re-seeding it with its own output converges to a fixed point that is
        // 15 percent better and still first order (measured 1.684e-2 ->
        // 6.789e-3 -> 3.227e-3 over five passes), so iteration is not the
        // missing ingredient; a wider stencil or a quadratic fit would be.
        // It is nonetheless about 1.7x better than the least-squares seed at
        // every level, which is why it is what `gradient_corrected` is fed.
        for (s, l) in seeded.iter().zip(lsq.iter()) {
            assert!(
                s < l,
                "the skew-seeded gradient should beat its own seed: {s:e} vs {l:e}"
            );
        }
        assert!(
            gg[1] > gg[0] && gg[2] > gg[1],
            "the Green-Gauss cell gradient is expected NOT to converge here; if it now \
             does, this characterisation is stale and the docs on \
             cell_gradient_green_gauss must be rewritten: {:e} {:e} {:e}",
            gg[0],
            gg[1],
            gg[2]
        );
    }

    #[test]
    fn the_explicit_correction_is_a_large_fraction_of_the_implicit_term() {
        // How big the deferred correction is decides whether a solver needs
        // under-relaxation. Two numbers: the purely geometric bound
        // |k_f| / |A_f| = tan(theta_f), and the ratio actually realised on
        // p = (x/R)^2 at facets where the implicit term is not itself small.
        for (n, layers) in [(8usize, 6usize), (16, 10), (32, 18)] {
            let t = mini(n, layers);
            let r = t.spec().r_top_m();
            let corr = t.face_corrections().expect("corrections");
            let centroids = t.cell_centroids_m().expect("centroids");
            let p: Vec<f64> = centroids.iter().map(|x| (x[0] / r).powi(2)).collect();
            let implicit = t.gradient_implicit(&p, &corr).expect("implicit");
            let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
            let seed = t
                .cell_gradient_least_squares(&p, &centroids)
                .expect("cell gradient");
            let cell_grad = t
                .cell_gradient_green_gauss(&p, &v_poly, &corr, Some(&seed))
                .expect("cell gradient");
            let explicit = t
                .gradient_nonorthogonal_correction(
                    &p,
                    &cell_grad,
                    &corr,
                    FaceGradientRule::NonOrthogonalAndSkew,
                )
                .expect("correction");

            let mut max_tan = 0.0f64;
            let mut max_theta = 0.0f64;
            let mut max_lambda = 0.0f64;
            // The two geometric claims the whole diagnosis rests on, measured
            // here rather than quoted: in-face lateral facets carry all the
            // non-orthogonality, cross-seam facets carry none and carry the
            // skewness instead.
            let mut in_face_theta = 0.0f64;
            let mut seam_theta = 0.0f64;
            let mut in_face_skew = 0.0f64;
            let mut seam_skew = 0.0f64;
            for (idx, face) in t.faces().iter().enumerate() {
                let Some(nb) = face.neighbor else { continue };
                max_tan = max_tan.max(norm(corr[idx].tangential_m2) / face.area_m2);
                max_theta = max_theta.max(corr[idx].non_orthogonality_rad);
                max_lambda = max_lambda.max(corr[idx].lambda);
                if face.kind != FaceKind::Lateral {
                    continue;
                }
                let (fo, ..) = t.spec().fijk(face.owner as u64).expect("in range");
                let (fn_, ..) = t.spec().fijk(nb as u64).expect("in range");
                let skew = norm(corr[idx].skew_m) / face.distance_m;
                if fo == fn_ {
                    in_face_theta = in_face_theta.max(corr[idx].non_orthogonality_rad);
                    in_face_skew = in_face_skew.max(skew);
                } else {
                    seam_theta = seam_theta.max(corr[idx].non_orthogonality_rad);
                    seam_skew = seam_skew.max(skew);
                }
            }
            println!(
                "n={n:<3} lateral facets: in-face max theta {:.4} deg skew {in_face_skew:.4} | \
                 cross-seam max theta {:.3e} deg skew {seam_skew:.4}",
                in_face_theta.to_degrees(),
                seam_theta.to_degrees()
            );
            assert!(
                seam_theta.to_degrees() < 1e-4,
                "n={n}: cross-seam facets are supposed to be EXACTLY orthogonal by mirror \
                 symmetry; measured {} deg",
                seam_theta.to_degrees()
            );
            assert!(
                in_face_theta > 0.3,
                "n={n}: in-face non-orthogonality collapsed to {} rad - the whole \
                 diagnosis assumed it does not shrink",
                in_face_theta
            );
            assert!(
                seam_skew > 0.2 && seam_skew > 4.0 * in_face_skew,
                "n={n}: seam skewness {seam_skew} vs in-face {in_face_skew} - the residual \
                 error mechanism assumed the seam value is large and does not shrink"
            );
            let ref_implicit = implicit.iter().fold(0.0f64, |a, b| a.max(b.abs()));
            let mut max_ratio = 0.0f64;
            for idx in 0..t.face_count() {
                if implicit[idx].abs() < 0.1 * ref_implicit {
                    continue;
                }
                max_ratio = max_ratio.max(explicit[idx].abs() / implicit[idx].abs());
            }
            println!(
                "n={n:<3} max non-orthogonality {:.2} deg  max lambda {max_lambda:.4}  \
                 geometric |k|/|A| = tan(theta) {max_tan:.4}  \
                 realised |explicit|/|implicit| {max_ratio:.4}",
                max_theta.to_degrees()
            );
            assert!(
                max_tan > 0.3,
                "n={n}: the correction is supposed to be large; tan(theta) max is {max_tan}"
            );
            assert!(
                (max_lambda - 1.0 / max_theta.cos()).abs() < 1e-9,
                "lambda and theta must be the same geometry: {max_lambda} vs {}",
                1.0 / max_theta.cos()
            );
        }
    }

    #[test]
    fn operators_reject_wrongly_sized_fields() {
        let t = mini(2, 2);
        let volumes = t.cell_volumes_m3().expect("volumes");
        assert!(t.divergence(&[0.0; 3], &volumes).is_err());
        assert!(t.gradient(&[0.0; 3]).is_err());
        assert!(t.reconstruct_cell_velocity(&[0.0; 3]).is_err());
        let mut short = [0.0; 3];
        assert!(t.remove_mean_pressure(&mut short, &volumes).is_err());

        let corr = t.face_corrections().expect("corrections");
        let v_poly = t.polyhedral_cell_volumes_m3().expect("polyhedral volumes");
        assert!(t.gradient_implicit(&[0.0; 3], &corr).is_err());
        assert!(t.face_masses_over_relaxed(&corr[..1]).is_err());
        assert!(t
            .cell_gradient_green_gauss(&[0.0; 3], &v_poly, &corr, None)
            .is_err());
        let centroids = t.cell_centroids_m().expect("centroids");
        assert!(t
            .cell_gradient_least_squares(&[0.0; 3], &centroids)
            .is_err());
        assert!(t
            .gradient_corrected(
                &[0.0; 3],
                &[[0.0; 3]; 3],
                &corr,
                FaceGradientRule::NonOrthogonalAndSkew
            )
            .is_err());
        assert!(t
            .green_gauss_gradient(&[[[0.0; 3]; 3]; 3], &v_poly)
            .is_err());
    }
}
