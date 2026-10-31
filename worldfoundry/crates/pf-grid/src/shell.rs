//! Cubed-sphere spherical shell: the mantle volume domain (MA).
//!
//! Brief rev 5 section M / [ADR 026]. This is the discretization the solved
//! mantle lives on: six gnomonic equiangular faces, `n x n` cells each,
//! `layers` radial shells between the core-mantle boundary and the base of
//! the lithosphere. It publishes through the existing `DomainDeclaration`
//! as `Volume3D` / `CubedSphere`, so it needs no artifact schema change.
//!
//! Why cubed sphere and not another lat/lon box: a lat/lon volume inherits
//! the pole problem, where cells shrink to nothing and the declared cap
//! policy has to defend the drainage graph from its own geometry. On the
//! cubed sphere every cell has exactly four lateral neighbors and no cell
//! degenerates, so there is no cap policy to declare and no polar timestep
//! collapse in the transport step. The cost is the eight cube corners, where
//! three cells meet instead of four; that is a mild metric distortion, not a
//! topological hole, and it is measured by the tests below rather than
//! assumed.
//!
//! Invariants validated at construction, fail-loud, never papered over:
//! the six faces' solid angles sum to 4*pi, and the cell volumes sum to the
//! shell volume. Those are this module's analogue of `GridSpec`'s
//! `sum area = 4*pi*R^2` gate.
//!
//! Conventions (also written into `transform_conventions`, so an artifact
//! reader never has to guess):
//! - element = ((face * n + i) * n + j) * layers + k, C-order, radial
//!   fastest, so one column is contiguous. Melting, geotherms and column
//!   integrals all walk k.
//! - **k = 0 is the DEEPEST layer** (at the core-mantle boundary);
//!   k = layers-1 is the shallowest (at the lithosphere base).
//! - face 0..5 = +X, -X, +Y, -Y, +Z, -Z by outward axis.
//! - i indexes the face's alpha angle, j its beta angle, both increasing.
//!
//! [ADR 026]: ../../../../docs/adr/026_physical_mantle_simulation_critical_path.md

use serde::{Deserialize, Serialize};

use pf_contracts::spatial::{
    CoordinateSystem, DatumDescriptor, DiscretizationDescriptor, DomainDeclaration, ElementId,
    NeighborLink, Position, SpatialDomainKind, SpatialError,
};

#[derive(Debug, thiserror::Error)]
pub enum ShellError {
    #[error("shell n={n} layers={layers} invalid: {why}")]
    Invalid {
        n: usize,
        layers: usize,
        why: &'static str,
    },
    #[error("shell radii invalid: bottom={bottom_m} top={top_m}: {why}")]
    Radii {
        bottom_m: f64,
        top_m: f64,
        why: &'static str,
    },
    #[error("shell geometry invariant broken: {what} relative error {rel_err:e} exceeds {tol:e}")]
    Invariant {
        what: &'static str,
        rel_err: f64,
        tol: f64,
    },
    #[error("cell index out of range: face={face} i={i} j={j} k={k} on n={n} layers={layers}")]
    Index {
        face: usize,
        i: usize,
        j: usize,
        k: usize,
        n: usize,
        layers: usize,
    },
}

/// Which lateral direction a neighbor lies in, in the source cell's own face
/// coordinates. Cross-face links keep the source-side label: "the neighbor
/// past my +alpha edge" is well defined even when that neighbor's own face
/// calls the shared edge something else.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Lateral {
    AlphaMinus,
    AlphaPlus,
    BetaMinus,
    BetaPlus,
}

impl Lateral {
    pub const ALL: [Lateral; 4] = [
        Lateral::AlphaMinus,
        Lateral::AlphaPlus,
        Lateral::BetaMinus,
        Lateral::BetaPlus,
    ];
}

/// One cell's link to a neighbor, with the geometry a finite-volume step
/// needs: the shared face area and the centroid separation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellLink {
    pub face: usize,
    pub i: usize,
    pub j: usize,
    pub k: usize,
    /// Which facet of the SOURCE cell this link crosses.
    pub facet: Facet,
    /// Shared interface area, m^2. Always `|area_vector_m2|`, so the scalar
    /// and the vector can never drift apart.
    pub area_m2: f64,
    /// Outward area vector of the shared facet, m^2, pointing OUT of the
    /// source cell. This is the primitive a conservative divergence uses;
    /// the scalar above is derived from it, not the other way round.
    pub area_vector_m2: [f64; 3],
    /// Centroid-to-centroid distance, m.
    pub distance_m: f64,
}

/// One of a cell's six facets, named from the source cell's point of view.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Facet {
    Lateral(Lateral),
    /// Toward the core (r decreasing).
    Down,
    /// Toward the surface (r increasing).
    Up,
}

impl Facet {
    pub const ALL: [Facet; 6] = [
        Facet::Lateral(Lateral::AlphaMinus),
        Facet::Lateral(Lateral::AlphaPlus),
        Facet::Lateral(Lateral::BetaMinus),
        Facet::Lateral(Lateral::BetaPlus),
        Facet::Down,
        Facet::Up,
    ];
}

/// Fields are PRIVATE: every construction path validates dimensions, radii
/// and both geometry invariants. There is no way to hold an invalid
/// `ShellSpec`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ShellSpec {
    n: usize,
    layers: usize,
    r_bottom_m: f64,
    r_top_m: f64,
    /// Ladder level when constructed via `preview()` / `standard()`; None
    /// for off-ladder shells (mini/test), recorded in the domain descriptor
    /// as `resolution_level = u32::MAX`, never silently defaulted.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    level: Option<u8>,
}

impl<'de> Deserialize<'de> for ShellSpec {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct Raw {
            n: usize,
            layers: usize,
            r_bottom_m: f64,
            r_top_m: f64,
            #[serde(default)]
            level: Option<u8>,
        }
        let raw = Raw::deserialize(d)?;
        let mut spec = ShellSpec::new(raw.n, raw.layers, raw.r_bottom_m, raw.r_top_m)
            .map_err(|e| serde::de::Error::custom(format!("invalid ShellSpec: {e}")))?;
        if let Some(level) = raw.level {
            // A claimed ladder level must match the physical shape exactly:
            // an L1 claim on an L0 shape is a contract violation, not a
            // silent mislabel (same rule as GridSpec).
            let (want_n, want_layers) = match level {
                0 => (PREVIEW_N, PREVIEW_LAYERS),
                1 => (STANDARD_N, STANDARD_LAYERS),
                _ => {
                    return Err(serde::de::Error::custom(format!(
                        "invalid ShellSpec: level={level} is not on the ADR 026 ladder"
                    )))
                }
            };
            if spec.n != want_n || spec.layers != want_layers {
                return Err(serde::de::Error::custom(format!(
                    "invalid ShellSpec: level={level} inconsistent with {}x{}x{}",
                    spec.n, spec.n, spec.layers
                )));
            }
        }
        spec.level = raw.level;
        Ok(spec)
    }
}

/// ADR 026's declared preview shape: 6 x 48 x 48 x 16 = 221,184 cells.
pub const PREVIEW_N: usize = 48;
pub const PREVIEW_LAYERS: usize = 16;
/// ADR 026's declared standard shape: 6 x 96 x 96 x 24 = 1,327,104 cells.
pub const STANDARD_N: usize = 96;
pub const STANDARD_LAYERS: usize = 24;

/// Relative tolerance on the construction-time geometry invariants. This is
/// a numerical-summation bound, not a physics knob: the analytic answers are
/// exactly 4*pi and the shell volume, and anything beyond this means the
/// geometry code is wrong.
const GEOMETRY_TOL: f64 = 1e-9;

impl ShellSpec {
    /// ADR 026 preview shape at the given radii.
    pub fn preview(r_bottom_m: f64, r_top_m: f64) -> Result<Self, ShellError> {
        let mut spec = Self::new(PREVIEW_N, PREVIEW_LAYERS, r_bottom_m, r_top_m)?;
        spec.level = Some(0);
        Ok(spec)
    }

    /// ADR 026 standard shape at the given radii.
    pub fn standard(r_bottom_m: f64, r_top_m: f64) -> Result<Self, ShellError> {
        let mut spec = Self::new(STANDARD_N, STANDARD_LAYERS, r_bottom_m, r_top_m)?;
        spec.level = Some(1);
        Ok(spec)
    }

    pub fn new(n: usize, layers: usize, r_bottom_m: f64, r_top_m: f64) -> Result<Self, ShellError> {
        if n < 2 {
            return Err(ShellError::Invalid {
                n,
                layers,
                why: "face resolution n must be >= 2",
            });
        }
        if !n.is_multiple_of(2) {
            // Even n keeps the face centre on a cell corner rather than
            // inside a cell, which keeps the six faces mutually symmetric.
            return Err(ShellError::Invalid {
                n,
                layers,
                why: "face resolution n must be even",
            });
        }
        // Upper bound mirrors pf-contracts' fidelity.radial_layers range.
        if !(1..=24).contains(&layers) {
            return Err(ShellError::Invalid {
                n,
                layers,
                why: "radial layers must be in 1..=24 (pf-contracts fidelity range)",
            });
        }
        if !r_bottom_m.is_finite() || !r_top_m.is_finite() {
            return Err(ShellError::Radii {
                bottom_m: r_bottom_m,
                top_m: r_top_m,
                why: "radii must be finite",
            });
        }
        if r_bottom_m <= 0.0 {
            return Err(ShellError::Radii {
                bottom_m: r_bottom_m,
                top_m: r_top_m,
                why: "bottom radius must be positive",
            });
        }
        if r_top_m <= r_bottom_m {
            return Err(ShellError::Radii {
                bottom_m: r_bottom_m,
                top_m: r_top_m,
                why: "top radius must exceed bottom radius",
            });
        }
        let spec = Self {
            n,
            layers,
            r_bottom_m,
            r_top_m,
            level: None,
        };
        spec.validate_geometry()?;
        Ok(spec)
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn layers(&self) -> usize {
        self.layers
    }

    pub fn r_bottom_m(&self) -> f64 {
        self.r_bottom_m
    }

    pub fn r_top_m(&self) -> f64 {
        self.r_top_m
    }

    pub fn ladder_level(&self) -> Option<u8> {
        self.level
    }

    /// Cells per face (lateral only, one radial layer).
    pub fn cells_per_face(&self) -> usize {
        self.n * self.n
    }

    /// Total element count over six faces and all radial layers.
    pub fn cells(&self) -> usize {
        6 * self.n * self.n * self.layers
    }

    // -----------------------------------------------------------------
    // Angular geometry

    /// Angular width of one cell in alpha or beta, radians.
    pub fn cell_angle_rad(&self) -> f64 {
        std::f64::consts::FRAC_PI_2 / self.n as f64
    }

    /// Angle of the i-th cell EDGE (i in 0..=n), radians in [-pi/4, pi/4].
    fn edge_angle(&self, i: usize) -> f64 {
        -std::f64::consts::FRAC_PI_4 + i as f64 * self.cell_angle_rad()
    }

    /// Angle of the i-th cell CENTER (i in 0..n), radians.
    fn center_angle(&self, i: usize) -> f64 {
        self.edge_angle(i) + 0.5 * self.cell_angle_rad()
    }

    /// Unit direction on `face` at face-local angles (alpha, beta).
    ///
    /// Gnomonic equiangular: the point is `out + tan(alpha) * a +
    /// tan(beta) * b` on the cube face, projected to the sphere. Angles
    /// beyond +-pi/4 are still meaningful directions (they simply lie past
    /// the face edge), which the neighbor search relies on.
    fn direction(&self, face: usize, alpha: f64, beta: f64) -> [f64; 3] {
        let (out, a, b) = face_basis(face);
        let (ta, tb) = (alpha.tan(), beta.tan());
        normalize([
            out[0] + ta * a[0] + tb * b[0],
            out[1] + ta * a[1] + tb * b[1],
            out[2] + ta * a[2] + tb * b[2],
        ])
    }

    /// Unit direction of a cell's lateral center.
    pub fn cell_direction(&self, face: usize, i: usize, j: usize) -> Result<[f64; 3], ShellError> {
        self.check_lateral(face, i, j)?;
        Ok(self.direction(face, self.center_angle(i), self.center_angle(j)))
    }

    /// The four unit corner directions of a cell, counter-clockwise in face
    /// coordinates: (i,j), (i+1,j), (i+1,j+1), (i,j+1).
    fn corners(&self, face: usize, i: usize, j: usize) -> [[f64; 3]; 4] {
        let (a0, a1) = (self.edge_angle(i), self.edge_angle(i + 1));
        let (b0, b1) = (self.edge_angle(j), self.edge_angle(j + 1));
        [
            self.direction(face, a0, b0),
            self.direction(face, a1, b0),
            self.direction(face, a1, b1),
            self.direction(face, a0, b1),
        ]
    }

    /// Solid angle subtended by one lateral cell, steradians. Exact: the
    /// spherical quadrilateral is split into two spherical triangles and
    /// each is evaluated with the Van Oosterom and Strackee formula, which
    /// stays well conditioned for the small triangles a cubed sphere makes.
    pub fn solid_angle_sr(&self, face: usize, i: usize, j: usize) -> Result<f64, ShellError> {
        self.check_lateral(face, i, j)?;
        let c = self.corners(face, i, j);
        Ok(triangle_solid_angle(c[0], c[1], c[2]) + triangle_solid_angle(c[0], c[2], c[3]))
    }

    /// Total solid angle over all six faces. Analytically 4*pi.
    pub fn total_solid_angle_sr(&self) -> f64 {
        let mut sum = 0.0;
        for face in 0..6 {
            for i in 0..self.n {
                for j in 0..self.n {
                    // Indices are in range by construction of the loops.
                    sum += self
                        .solid_angle_sr(face, i, j)
                        .expect("loop indices are in range");
                }
            }
        }
        sum
    }

    // -----------------------------------------------------------------
    // Radial geometry

    /// Radius of the k-th layer INTERFACE (k in 0..=layers), meters, with
    /// k = 0 at the core-mantle boundary. Uniform in radius: a declared
    /// approximation, replaceable by a depth-weighted layering without
    /// touching anything above this line.
    pub fn interface_radius_m(&self, k: usize) -> f64 {
        let t = k as f64 / self.layers as f64;
        self.r_bottom_m + t * (self.r_top_m - self.r_bottom_m)
    }

    /// Radius of the k-th layer's centroid shell, meters. This is the
    /// volume-bisecting radius, not the arithmetic midpoint, so that a
    /// cell's centroid actually sits at its centre of volume.
    pub fn center_radius_m(&self, k: usize) -> Result<f64, ShellError> {
        if k >= self.layers {
            return Err(ShellError::Index {
                face: 0,
                i: 0,
                j: 0,
                k,
                n: self.n,
                layers: self.layers,
            });
        }
        let r0 = self.interface_radius_m(k);
        let r1 = self.interface_radius_m(k + 1);
        Ok((0.5 * (r0.powi(3) + r1.powi(3))).cbrt())
    }

    /// Radial thickness of layer k, meters.
    pub fn layer_thickness_m(&self, k: usize) -> Result<f64, ShellError> {
        if k >= self.layers {
            return Err(ShellError::Index {
                face: 0,
                i: 0,
                j: 0,
                k,
                n: self.n,
                layers: self.layers,
            });
        }
        Ok(self.interface_radius_m(k + 1) - self.interface_radius_m(k))
    }

    /// Cell volume, m^3. Exact given the solid angle: the cell is a
    /// spherical frustum, volume = solid_angle * (r1^3 - r0^3) / 3.
    pub fn cell_volume_m3(
        &self,
        face: usize,
        i: usize,
        j: usize,
        k: usize,
    ) -> Result<f64, ShellError> {
        self.check(face, i, j, k)?;
        let omega = self.solid_angle_sr(face, i, j)?;
        let r0 = self.interface_radius_m(k);
        let r1 = self.interface_radius_m(k + 1);
        Ok(omega * (r1.powi(3) - r0.powi(3)) / 3.0)
    }

    /// Analytic shell volume, m^3.
    pub fn shell_volume_m3(&self) -> f64 {
        4.0 / 3.0 * std::f64::consts::PI * (self.r_top_m.powi(3) - self.r_bottom_m.powi(3))
    }

    /// Sum of all cell volumes. Equals `shell_volume_m3` analytically.
    pub fn total_volume_m3(&self) -> f64 {
        let mut sum = 0.0;
        for face in 0..6 {
            for i in 0..self.n {
                for j in 0..self.n {
                    let omega = self
                        .solid_angle_sr(face, i, j)
                        .expect("loop indices are in range");
                    for k in 0..self.layers {
                        let r0 = self.interface_radius_m(k);
                        let r1 = self.interface_radius_m(k + 1);
                        sum += omega * (r1.powi(3) - r0.powi(3)) / 3.0;
                    }
                }
            }
        }
        sum
    }

    /// Both construction-time geometry gates. Fail-loud with the measured
    /// relative error, never a silent clamp.
    pub fn validate_geometry(&self) -> Result<(), ShellError> {
        let solid = self.total_solid_angle_sr();
        let want_solid = 4.0 * std::f64::consts::PI;
        let rel = ((solid - want_solid) / want_solid).abs();
        if rel.is_nan() || rel > GEOMETRY_TOL {
            return Err(ShellError::Invariant {
                what: "sum of cell solid angles != 4*pi",
                rel_err: rel,
                tol: GEOMETRY_TOL,
            });
        }
        let vol = self.total_volume_m3();
        let want_vol = self.shell_volume_m3();
        let rel = ((vol - want_vol) / want_vol).abs();
        if rel.is_nan() || rel > GEOMETRY_TOL {
            return Err(ShellError::Invariant {
                what: "sum of cell volumes != shell volume",
                rel_err: rel,
                tol: GEOMETRY_TOL,
            });
        }
        Ok(())
    }

    // -----------------------------------------------------------------
    // Indexing

    fn check_lateral(&self, face: usize, i: usize, j: usize) -> Result<(), ShellError> {
        if face >= 6 || i >= self.n || j >= self.n {
            return Err(ShellError::Index {
                face,
                i,
                j,
                k: 0,
                n: self.n,
                layers: self.layers,
            });
        }
        Ok(())
    }

    fn check(&self, face: usize, i: usize, j: usize, k: usize) -> Result<(), ShellError> {
        self.check_lateral(face, i, j)?;
        if k >= self.layers {
            return Err(ShellError::Index {
                face,
                i,
                j,
                k,
                n: self.n,
                layers: self.layers,
            });
        }
        Ok(())
    }

    /// Flat element index: radial fastest, so a column is contiguous.
    pub fn flat(&self, face: usize, i: usize, j: usize, k: usize) -> u64 {
        (((face * self.n + i) * self.n + j) * self.layers + k) as u64
    }

    /// Decode a flat element index; fails loud out of range.
    pub fn fijk(&self, element: u64) -> Result<(usize, usize, usize, usize), SpatialError> {
        let count = self.cells() as u64;
        let err = SpatialError {
            element,
            element_count: count,
        };
        let e = usize::try_from(element).map_err(|_| err)?;
        if e >= self.cells() {
            return Err(SpatialError {
                element,
                element_count: count,
            });
        }
        let k = e % self.layers;
        let rest = e / self.layers;
        let j = rest % self.n;
        let rest = rest / self.n;
        let i = rest % self.n;
        let face = rest / self.n;
        Ok((face, i, j, k))
    }

    /// Which lateral cell contains a direction. Total on the unit sphere:
    /// the dominant axis picks the face, the face-local angles pick (i, j).
    pub fn locate(&self, v: [f64; 3]) -> Result<(usize, usize, usize), ShellError> {
        let v = normalize(v);
        let (ax, ay, az) = (v[0].abs(), v[1].abs(), v[2].abs());
        let face = if ax >= ay && ax >= az {
            if v[0] > 0.0 {
                0
            } else {
                1
            }
        } else if ay >= az {
            if v[1] > 0.0 {
                2
            } else {
                3
            }
        } else if v[2] > 0.0 {
            4
        } else {
            5
        };
        let (out, a, b) = face_basis(face);
        let p = dot(v, out);
        if p.is_nan() || p <= 0.0 {
            // Only reachable for a non-finite or zero input vector; the
            // dominant-axis choice guarantees p > 0 otherwise.
            return Err(ShellError::Invariant {
                what: "direction does not project onto its own face",
                rel_err: p,
                tol: 0.0,
            });
        }
        let alpha = (dot(v, a) / p).atan();
        let beta = (dot(v, b) / p).atan();
        let d = self.cell_angle_rad();
        let idx = |ang: f64| -> usize {
            let raw = (ang + std::f64::consts::FRAC_PI_4) / d;
            let f = raw.floor();
            if f < 0.0 {
                0
            } else if f >= self.n as f64 {
                self.n - 1
            } else {
                f as usize
            }
        };
        Ok((face, idx(alpha), idx(beta)))
    }

    // -----------------------------------------------------------------
    // Topology

    /// The lateral neighbor across one of a cell's four edges.
    ///
    /// In-face steps are index arithmetic. A step across a cube edge is
    /// resolved geometrically: reflect the cell centre through the midpoint
    /// of the shared edge and ask which cell contains the result. The six
    /// faces divide every shared edge into the same n equal angular
    /// intervals, so the reflected point lands inside the mirror cell rather
    /// than near a boundary. Extending the source face's own angle past
    /// pi/4 would NOT work: the neighbouring face's beta coordinate is
    /// distorted by the projection, and near the ends of an edge the naive
    /// step lands exactly on a cell boundary. The symmetry test below is
    /// what keeps this honest.
    pub fn lateral_neighbor(
        &self,
        face: usize,
        i: usize,
        j: usize,
        dir: Lateral,
    ) -> Result<(usize, usize, usize), ShellError> {
        self.check_lateral(face, i, j)?;
        let last = self.n - 1;
        let in_face = match dir {
            Lateral::AlphaMinus if i > 0 => Some((face, i - 1, j)),
            Lateral::AlphaPlus if i < last => Some((face, i + 1, j)),
            Lateral::BetaMinus if j > 0 => Some((face, i, j - 1)),
            Lateral::BetaPlus if j < last => Some((face, i, j + 1)),
            _ => None,
        };
        if let Some(t) = in_face {
            return Ok(t);
        }
        // Crossing a cube edge.
        let (a_c, b_c) = (self.center_angle(i), self.center_angle(j));
        let mid = match dir {
            Lateral::AlphaMinus => self.direction(face, self.edge_angle(i), b_c),
            Lateral::AlphaPlus => self.direction(face, self.edge_angle(i + 1), b_c),
            Lateral::BetaMinus => self.direction(face, a_c, self.edge_angle(j)),
            Lateral::BetaPlus => self.direction(face, a_c, self.edge_angle(j + 1)),
        };
        let c = self.direction(face, a_c, b_c);
        let refl = normalize([
            2.0 * dot(c, mid) * mid[0] - c[0],
            2.0 * dot(c, mid) * mid[1] - c[1],
            2.0 * dot(c, mid) * mid[2] - c[2],
        ]);
        let (nf, ni, nj) = self.locate(refl)?;
        Ok((nf, ni, nj))
    }

    /// Arc length of the edge shared with the neighbor in `dir`, on the unit
    /// sphere, radians. Seam diagnostics use it to check that two faces
    /// divide a shared cube edge identically; the flux geometry uses the
    /// facet area vectors, not this.
    pub fn shared_edge_angle(&self, face: usize, i: usize, j: usize, dir: Lateral) -> f64 {
        let (p, q) = match dir {
            Lateral::AlphaMinus => (
                self.direction(face, self.edge_angle(i), self.edge_angle(j)),
                self.direction(face, self.edge_angle(i), self.edge_angle(j + 1)),
            ),
            Lateral::AlphaPlus => (
                self.direction(face, self.edge_angle(i + 1), self.edge_angle(j)),
                self.direction(face, self.edge_angle(i + 1), self.edge_angle(j + 1)),
            ),
            Lateral::BetaMinus => (
                self.direction(face, self.edge_angle(i), self.edge_angle(j)),
                self.direction(face, self.edge_angle(i + 1), self.edge_angle(j)),
            ),
            Lateral::BetaPlus => (
                self.direction(face, self.edge_angle(i), self.edge_angle(j + 1)),
                self.direction(face, self.edge_angle(i + 1), self.edge_angle(j + 1)),
            ),
        };
        angle_between(p, q)
    }

    /// The eight vertices of a cell: the four lateral corner directions at
    /// the layer's bottom radius (indices 0..4), then the same four at its
    /// top radius (indices 4..8), in the corner order
    /// (i,j), (i+1,j), (i+1,j+1), (i,j+1).
    pub fn cell_vertices(
        &self,
        face: usize,
        i: usize,
        j: usize,
        k: usize,
    ) -> Result<[[f64; 3]; 8], ShellError> {
        self.check(face, i, j, k)?;
        let c = self.corners(face, i, j);
        let r0 = self.interface_radius_m(k);
        let r1 = self.interface_radius_m(k + 1);
        Ok([
            scale(c[0], r0),
            scale(c[1], r0),
            scale(c[2], r0),
            scale(c[3], r0),
            scale(c[0], r1),
            scale(c[1], r1),
            scale(c[2], r1),
            scale(c[3], r1),
        ])
    }

    /// The four vertices of one facet, in a cycle.
    fn facet_vertices(
        &self,
        face: usize,
        i: usize,
        j: usize,
        k: usize,
        facet: Facet,
    ) -> Result<[[f64; 3]; 4], ShellError> {
        let v = self.cell_vertices(face, i, j, k)?;
        // Bottom corners are v[0..4]; the same corner at the top radius is
        // v[c + 4].
        let lateral = |a: usize, b: usize| [v[a], v[b], v[b + 4], v[a + 4]];
        Ok(match facet {
            Facet::Lateral(Lateral::AlphaMinus) => lateral(0, 3),
            Facet::Lateral(Lateral::AlphaPlus) => lateral(1, 2),
            Facet::Lateral(Lateral::BetaMinus) => lateral(0, 1),
            Facet::Lateral(Lateral::BetaPlus) => lateral(3, 2),
            Facet::Down => [v[0], v[1], v[2], v[3]],
            Facet::Up => [v[4], v[5], v[6], v[7]],
        })
    }

    /// Outward area vector of one facet, m^2.
    ///
    /// The cell is treated as a hexahedron on its eight vertices and each
    /// facet's vector area is the triangulation-independent form
    /// `0.5 * (v2 - v0) x (v3 - v1)`. That choice is what makes the six
    /// outward vectors of a closed cell sum to zero to roundoff, which in
    /// turn makes the discrete divergence annihilate a constant field
    /// exactly. The orientation is fixed against a geometric reference
    /// direction rather than trusted to vertex winding, and
    /// `closed_cell_area_vectors_sum_to_zero` is the test that keeps that
    /// honest at seams and cube corners as well as in face interiors.
    ///
    /// NOTE the deliberate split of duties: these polyhedral areas are the
    /// FLUX geometry, while `solid_angle_sr` and `cell_volume_m3` remain the
    /// exact spherical MASS geometry. They differ at O(h^2), which
    /// `polyhedral_radial_area_converges_to_the_sphere` measures rather than
    /// hides. Conservation does not depend on the two agreeing: a discrete
    /// divergence is exact for a constant field whenever the outward area
    /// vectors close, whatever volume divides them.
    pub fn facet_area_vector_m2(
        &self,
        face: usize,
        i: usize,
        j: usize,
        k: usize,
        facet: Facet,
    ) -> Result<[f64; 3], ShellError> {
        let q = self.facet_vertices(face, i, j, k, facet)?;
        let d1 = sub(q[2], q[0]);
        let d2 = sub(q[3], q[1]);
        let mut a = scale_raw(cross(d1, d2), 0.5);
        let reference = self.facet_outward_reference(face, i, j, facet)?;
        if dot(a, reference) < 0.0 {
            a = scale_raw(a, -1.0);
        }
        Ok(a)
    }

    /// A direction guaranteed to lie on the outward side of a facet, used
    /// only to orient the vector area.
    fn facet_outward_reference(
        &self,
        face: usize,
        i: usize,
        j: usize,
        facet: Facet,
    ) -> Result<[f64; 3], ShellError> {
        let here = self.cell_direction(face, i, j)?;
        Ok(match facet {
            Facet::Up => here,
            Facet::Down => scale_raw(here, -1.0),
            Facet::Lateral(dir) => {
                let (a_c, b_c) = (self.center_angle(i), self.center_angle(j));
                let mid = match dir {
                    Lateral::AlphaMinus => self.direction(face, self.edge_angle(i), b_c),
                    Lateral::AlphaPlus => self.direction(face, self.edge_angle(i + 1), b_c),
                    Lateral::BetaMinus => self.direction(face, a_c, self.edge_angle(j)),
                    Lateral::BetaPlus => self.direction(face, a_c, self.edge_angle(j + 1)),
                };
                sub(mid, here)
            }
        })
    }

    /// Every neighbor of a cell: up to four lateral and up to two radial,
    /// each carrying the outward area vector of the facet crossed.
    pub fn neighbors(
        &self,
        face: usize,
        i: usize,
        j: usize,
        k: usize,
    ) -> Result<Vec<ShellLink>, ShellError> {
        self.check(face, i, j, k)?;
        let mut links = Vec::with_capacity(6);
        let r_c = self.center_radius_m(k)?;
        let here = self.cell_direction(face, i, j)?;

        for dir in Lateral::ALL {
            let (nf, ni, nj) = self.lateral_neighbor(face, i, j, dir)?;
            let there = self.cell_direction(nf, ni, nj)?;
            let facet = Facet::Lateral(dir);
            let av = self.facet_area_vector_m2(face, i, j, k, facet)?;
            links.push(ShellLink {
                face: nf,
                i: ni,
                j: nj,
                k,
                facet,
                area_m2: norm(av),
                area_vector_m2: av,
                distance_m: angle_between(here, there) * r_c,
            });
        }

        if k > 0 {
            let av = self.facet_area_vector_m2(face, i, j, k, Facet::Down)?;
            links.push(ShellLink {
                face,
                i,
                j,
                k: k - 1,
                facet: Facet::Down,
                area_m2: norm(av),
                area_vector_m2: av,
                distance_m: (r_c - self.center_radius_m(k - 1)?).abs(),
            });
        }
        if k + 1 < self.layers {
            let av = self.facet_area_vector_m2(face, i, j, k, Facet::Up)?;
            links.push(ShellLink {
                face,
                i,
                j,
                k: k + 1,
                facet: Facet::Up,
                area_m2: norm(av),
                area_vector_m2: av,
                distance_m: (self.center_radius_m(k + 1)? - r_c).abs(),
            });
        }
        Ok(links)
    }

    // -----------------------------------------------------------------
    // Contract surface

    /// Domain declaration for this shell as a `Volume3D` cubed-sphere
    /// domain, with computed topology and metric digests.
    ///
    /// The datum radius is the shell's OUTER radius, so the metric digest
    /// changes when the planet changes size even if the cell counts do not.
    pub fn domain_declaration(
        &self,
        domain_id: &str,
    ) -> Result<DomainDeclaration, pf_contracts::ContractError> {
        let resolution_level = self.level.map(u32::from).unwrap_or(u32::MAX);
        let mut decl = DomainDeclaration {
            domain_id: domain_id.into(),
            kind: SpatialDomainKind::Volume3D,
            discretization: DiscretizationDescriptor {
                topology_id: format!("mantle.cubedsphere.{}x{}x{}", self.n, self.n, self.layers),
                backend: "cubedsphere-v0".into(),
                coordinate_system: CoordinateSystem::CubedSphere,
                topology_digest: String::new(),
                metric_digest: String::new(),
                element_count: self.cells() as u64,
                resolution_level,
                transform_conventions: format!(
                    "element=((face*n+i)*n+j)*layers+k (C-order, radial fastest); \
                     n={n}; layers={layers}; k=0 is the DEEPEST layer at r_bottom, \
                     k=layers-1 the shallowest at r_top; faces 0..5 = +X,-X,+Y,-Y,+Z,-Z \
                     by outward axis; equiangular gnomonic, i indexes alpha and j beta, \
                     both increasing over [-pi/4, pi/4]; radial layering uniform in \
                     radius; every cell has exactly 4 lateral neighbours (no pole caps)",
                    n = self.n,
                    layers = self.layers
                ),
            },
            datum: DatumDescriptor {
                center: [0.0, 0.0, 0.0],
                radius_m: self.r_top_m,
            },
        };
        decl.compute_digests()?;
        Ok(decl)
    }

    fn spatial_err(&self, element: u64) -> SpatialError {
        SpatialError {
            element,
            element_count: self.cells() as u64,
        }
    }
}

impl pf_contracts::spatial::SpatialDiscretization for ShellSpec {
    fn domain(&self) -> DomainDeclaration {
        // domain_id is semantic; the backend supplies a conventional
        // default and artifacts may declare their own.
        self.domain_declaration("mantle.primary")
            .expect("ShellSpec invariants validated at construction")
    }

    fn element_count(&self) -> u64 {
        self.cells() as u64
    }

    fn neighbors(&self, element: ElementId) -> Result<Vec<NeighborLink>, SpatialError> {
        let (face, i, j, k) = self.fijk(element.0)?;
        let links = self
            .neighbors(face, i, j, k)
            .map_err(|_| self.spatial_err(element.0))?;
        Ok(links
            .into_iter()
            .map(|l| NeighborLink {
                element: ElementId(self.flat(l.face, l.i, l.j, l.k)),
                distance_m: l.distance_m,
            })
            .collect())
    }

    fn centroid(&self, element: ElementId) -> Result<Position, SpatialError> {
        let (face, i, j, k) = self.fijk(element.0)?;
        let d = self
            .cell_direction(face, i, j)
            .map_err(|_| self.spatial_err(element.0))?;
        let r = self
            .center_radius_m(k)
            .map_err(|_| self.spatial_err(element.0))?;
        Ok([r * d[0], r * d[1], r * d[2]])
    }
}

impl pf_contracts::spatial::VolumeMetrics for ShellSpec {
    fn volume_m3(&self, element: ElementId) -> Result<f64, SpatialError> {
        let (face, i, j, k) = self.fijk(element.0)?;
        self.cell_volume_m3(face, i, j, k)
            .map_err(|_| self.spatial_err(element.0))
    }
}

impl pf_contracts::spatial::FaceTopology for ShellSpec {
    fn faces(&self, element: ElementId) -> Result<Vec<(ElementId, f64)>, SpatialError> {
        let (face, i, j, k) = self.fijk(element.0)?;
        let links = self
            .neighbors(face, i, j, k)
            .map_err(|_| self.spatial_err(element.0))?;
        Ok(links
            .into_iter()
            .map(|l| (ElementId(self.flat(l.face, l.i, l.j, l.k)), l.area_m2))
            .collect())
    }
}

impl pf_contracts::spatial::MetricSpace for ShellSpec {
    fn distance_m(&self, a: ElementId, b: ElementId) -> Result<f64, SpatialError> {
        use pf_contracts::spatial::SpatialDiscretization;
        let pa = self.centroid(a)?;
        let pb = self.centroid(b)?;
        Ok(((pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2) + (pa[2] - pb[2]).powi(2)).sqrt())
    }
}

// ---------------------------------------------------------------------------
// Vector helpers (private; this module owns its own small kit rather than
// exporting a vector type that would then need a contract of its own)

/// Outward axis and the two face-local axes, in the geocentric frame. The
/// pairings are chosen so that adjacent faces divide their shared edge into
/// the same n equal angular intervals; `cross_face_edges_are_shared` proves
/// it rather than trusting the table.
fn face_basis(face: usize) -> ([f64; 3], [f64; 3], [f64; 3]) {
    match face {
        0 => ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        1 => ([-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]),
        2 => ([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        3 => ([0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        4 => ([0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]),
        _ => ([0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]),
    }
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

fn scale_raw(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Scale a UNIT direction to a radius.
fn scale(v: [f64; 3], r: f64) -> [f64; 3] {
    scale_raw(v, r)
}

fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let m = dot(v, v).sqrt();
    [v[0] / m, v[1] / m, v[2] / m]
}

/// Angle between two unit vectors, numerically stable near zero via atan2
/// rather than acos of a dot product.
fn angle_between(a: [f64; 3], b: [f64; 3]) -> f64 {
    let c = cross(a, b);
    dot(c, c).sqrt().atan2(dot(a, b))
}

/// Van Oosterom and Strackee: the solid angle of the spherical triangle on
/// unit vectors a, b, c is 2 * atan2(|a . (b x c)|, 1 + a.b + b.c + c.a).
fn triangle_solid_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let num = dot(a, cross(b, c)).abs();
    let den = 1.0 + dot(a, b) + dot(b, c) + dot(c, a);
    2.0 * num.atan2(den)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_contracts::spatial::{FaceTopology, MetricSpace, SpatialDiscretization, VolumeMetrics};

    // Earth-like mantle shell: core-mantle boundary to lithosphere base.
    // These are test radii, not the project's physical constants; the
    // parameter profile owns those.
    const R_CMB: f64 = 3.480e6;
    const R_TOP: f64 = 6.271e6;

    fn mini() -> ShellSpec {
        ShellSpec::new(8, 4, R_CMB, R_TOP).expect("mini shell is valid")
    }

    #[test]
    fn solid_angles_sum_to_four_pi() {
        for n in [2usize, 4, 8, 16, 48] {
            let s = ShellSpec::new(n, 2, R_CMB, R_TOP).expect("valid");
            let rel = (s.total_solid_angle_sr() - 4.0 * std::f64::consts::PI).abs()
                / (4.0 * std::f64::consts::PI);
            assert!(rel < 1e-12, "n={n} relative error {rel:e}");
        }
    }

    #[test]
    fn volumes_sum_to_shell_volume() {
        let s = ShellSpec::new(16, 6, R_CMB, R_TOP).expect("valid");
        let rel = (s.total_volume_m3() - s.shell_volume_m3()).abs() / s.shell_volume_m3();
        assert!(rel < 1e-12, "relative error {rel:e}");
    }

    #[test]
    fn every_cell_has_four_lateral_neighbors_and_no_caps() {
        // The property that motivates the cubed sphere: unlike a lat/lon
        // volume there is no degenerate row, so there is no cap policy.
        let s = mini();
        for face in 0..6 {
            for i in 0..s.n() {
                for j in 0..s.n() {
                    let mut seen = std::collections::HashSet::new();
                    for dir in Lateral::ALL {
                        let t = s.lateral_neighbor(face, i, j, dir).expect("in range");
                        assert_ne!(t, (face, i, j), "cell is its own neighbour");
                        assert!(seen.insert(t), "duplicate lateral neighbour {t:?}");
                    }
                    assert_eq!(seen.len(), 4);
                }
            }
        }
    }

    #[test]
    fn lateral_adjacency_is_symmetric_including_across_cube_edges() {
        // This is the test that keeps the geometric cross-face search
        // honest: if the reflection landed in the wrong cell anywhere, some
        // pair would disagree about being neighbours.
        for n in [2usize, 4, 6, 8, 12, 16] {
            let s = ShellSpec::new(n, 1, R_CMB, R_TOP).expect("valid");
            let mut cross_face = 0usize;
            for face in 0..6 {
                for i in 0..n {
                    for j in 0..n {
                        for dir in Lateral::ALL {
                            let (nf, ni, nj) = s.lateral_neighbor(face, i, j, dir).expect("range");
                            if nf != face {
                                cross_face += 1;
                            }
                            let back = Lateral::ALL
                                .iter()
                                .map(|d| s.lateral_neighbor(nf, ni, nj, *d).expect("range"))
                                .any(|t| t == (face, i, j));
                            assert!(
                                back,
                                "n={n}: ({face},{i},{j}) -> ({nf},{ni},{nj}) is not mutual"
                            );
                        }
                    }
                }
            }
            // 12 cube edges, n cells along each, counted from both sides.
            assert_eq!(cross_face, 24 * n, "n={n} cross-face link count");
        }
    }

    #[test]
    fn cross_face_edges_are_shared_exactly() {
        // Adjacent faces must divide a shared cube edge into the SAME n
        // intervals, otherwise no neighbour mapping can be exact. Check the
        // shared edge arc lengths agree from both sides.
        let s = ShellSpec::new(8, 1, R_CMB, R_TOP).expect("valid");
        let n = s.n();
        for face in 0..6 {
            for j in 0..n {
                let (nf, ni, nj) = s
                    .lateral_neighbor(face, n - 1, j, Lateral::AlphaPlus)
                    .expect("range");
                assert_ne!(nf, face);
                let mine = s.shared_edge_angle(face, n - 1, j, Lateral::AlphaPlus);
                // Find the neighbour's own label for the same edge.
                let theirs = Lateral::ALL
                    .iter()
                    .find(|d| {
                        s.lateral_neighbor(nf, ni, nj, **d).expect("range") == (face, n - 1, j)
                    })
                    .map(|d| s.shared_edge_angle(nf, ni, nj, *d))
                    .expect("mutual link exists");
                assert!(
                    (mine - theirs).abs() < 1e-12,
                    "edge arc mismatch {mine} vs {theirs}"
                );
            }
        }
    }

    #[test]
    fn locate_round_trips_every_cell_centre() {
        let s = ShellSpec::new(16, 1, R_CMB, R_TOP).expect("valid");
        for face in 0..6 {
            for i in 0..s.n() {
                for j in 0..s.n() {
                    let d = s.cell_direction(face, i, j).expect("range");
                    assert_eq!(s.locate(d).expect("locates"), (face, i, j));
                }
            }
        }
    }

    #[test]
    fn element_index_round_trips() {
        let s = mini();
        for e in 0..s.cells() as u64 {
            let (f, i, j, k) = s.fijk(e).expect("in range");
            assert_eq!(s.flat(f, i, j, k), e);
        }
        assert!(s.fijk(s.cells() as u64).is_err());
    }

    #[test]
    fn radial_convention_puts_layer_zero_deepest() {
        let s = mini();
        let deep = s.center_radius_m(0).expect("k=0");
        let shallow = s.center_radius_m(s.layers() - 1).expect("top k");
        assert!(deep < shallow);
        assert!(deep > s.r_bottom_m() && shallow < s.r_top_m());
        // Centroid radii bisect their layer's volume, so they sit above the
        // arithmetic midpoint.
        let mid = 0.5 * (s.interface_radius_m(0) + s.interface_radius_m(1));
        assert!(deep > mid);
    }

    #[test]
    fn all_volumes_positive_and_finite() {
        let s = mini();
        for e in 0..s.cells() as u64 {
            let v = s.volume_m3(ElementId(e)).expect("in range");
            assert!(v.is_finite() && v > 0.0, "element {e} volume {v}");
        }
    }

    #[test]
    fn cell_volumes_are_within_a_bounded_factor_of_uniform() {
        // Cubed-sphere cells are not equal-area; the corner-to-centre ratio
        // is the price of having no poles. Measure it rather than assume it.
        let s = ShellSpec::new(16, 1, R_CMB, R_TOP).expect("valid");
        let mut min = f64::INFINITY;
        let mut max: f64 = 0.0;
        for face in 0..6 {
            for i in 0..s.n() {
                for j in 0..s.n() {
                    let a = s.solid_angle_sr(face, i, j).expect("range");
                    min = min.min(a);
                    max = max.max(a);
                }
            }
        }
        let ratio = max / min;
        assert!(
            ratio > 1.0 && ratio < 2.0,
            "solid-angle max/min ratio {ratio} outside the expected cubed-sphere range"
        );
    }

    #[test]
    fn trait_face_areas_are_shared_consistently_between_neighbors() {
        // A conservative finite-volume step needs both cells to agree on the
        // interface area, in both the lateral and the radial direction.
        let s = mini();
        for e in 0..s.cells() as u64 {
            for (other, area) in s.faces(ElementId(e)).expect("in range") {
                let back = s
                    .faces(other)
                    .expect("in range")
                    .into_iter()
                    .find(|(id, _)| *id == ElementId(e))
                    .map(|(_, a)| a)
                    .expect("mutual face");
                let rel = (area - back).abs() / area;
                assert!(rel < 1e-12, "interface area mismatch rel {rel:e}");
            }
        }
    }

    #[test]
    fn polyhedral_radial_area_converges_to_the_sphere() {
        // The flux geometry (polyhedral facets) and the mass geometry
        // (exact spherical) differ at O(h^2). Measure that gap and show it
        // shrinks with resolution rather than asserting the two are equal,
        // which they are not and must not be assumed to be.
        let mut prev: Option<f64> = None;
        for n in [8usize, 16, 32] {
            let s = ShellSpec::new(n, 4, R_CMB, R_TOP).expect("valid");
            let k = 1usize;
            let r = s.interface_radius_m(k);
            let mut sum = 0.0;
            for face in 0..6 {
                for i in 0..s.n() {
                    for j in 0..s.n() {
                        let links = s.neighbors(face, i, j, k).expect("range");
                        let down = links
                            .iter()
                            .find(|l| l.facet == Facet::Down)
                            .expect("layer 1 has a layer below");
                        sum += down.area_m2;
                    }
                }
            }
            let want = 4.0 * std::f64::consts::PI * r * r;
            let rel = (want - sum) / want;
            assert!(rel > 0.0, "flat facets must underestimate the sphere");
            if let Some(p) = prev {
                let ratio = p / rel;
                assert!(
                    ratio > 3.0 && ratio < 5.0,
                    "n={n}: second-order convergence expected, saw ratio {ratio}"
                );
            }
            prev = Some(rel);
        }
    }

    #[test]
    fn closed_cell_area_vectors_sum_to_zero() {
        // Discrete Gauss: the six outward area vectors of a closed cell must
        // cancel, or the divergence operator will not annihilate a constant
        // field and every conservation claim above it is decoration. Cube
        // corners and cross-face seams are included, and asserted to be
        // included, because that is exactly where a plausible-looking grid
        // stops being closed.
        let s = ShellSpec::new(8, 4, R_CMB, R_TOP).expect("valid");
        let mut corner_cells = 0usize;
        let mut seam_cells = 0usize;
        let last = s.n() - 1;
        for face in 0..6 {
            for i in 0..s.n() {
                for j in 0..s.n() {
                    if (i == 0 || i == last) && (j == 0 || j == last) {
                        corner_cells += 1;
                    }
                    if i == 0 || i == last || j == 0 || j == last {
                        seam_cells += 1;
                    }
                    for k in 0..s.layers() {
                        let mut sum = [0.0f64; 3];
                        let mut scale = 0.0f64;
                        for facet in Facet::ALL {
                            let a = s
                                .facet_area_vector_m2(face, i, j, k, facet)
                                .expect("in range");
                            sum = [sum[0] + a[0], sum[1] + a[1], sum[2] + a[2]];
                            scale += norm(a);
                        }
                        let rel = norm(sum) / scale;
                        assert!(
                            rel < 1e-12,
                            "cell ({face},{i},{j},{k}) is not closed: relative residual {rel:e}"
                        );
                    }
                }
            }
        }
        assert_eq!(corner_cells, 24, "four corner cells on each of six faces");
        assert_eq!(seam_cells, 6 * (4 * s.n() - 4));
    }

    #[test]
    fn paired_internal_facets_have_equal_area_and_opposite_normals() {
        // Two cells sharing a facet must see the same interface with
        // opposite orientation, including across cube seams where the two
        // sides reach the same vertices through different face
        // parameterizations.
        let s = ShellSpec::new(8, 3, R_CMB, R_TOP).expect("valid");
        let mut checked_cross_face = 0usize;
        for e in 0..s.cells() as u64 {
            let (face, i, j, k) = s.fijk(e).expect("in range");
            for link in s.neighbors(face, i, j, k).expect("in range") {
                let back = s
                    .neighbors(link.face, link.i, link.j, link.k)
                    .expect("in range")
                    .into_iter()
                    .find(|l| (l.face, l.i, l.j, l.k) == (face, i, j, k))
                    .expect("mutual facet");
                let rel_area = (link.area_m2 - back.area_m2).abs() / link.area_m2;
                assert!(rel_area < 1e-12, "facet area mismatch rel {rel_area:e}");
                let sum = [
                    link.area_vector_m2[0] + back.area_vector_m2[0],
                    link.area_vector_m2[1] + back.area_vector_m2[1],
                    link.area_vector_m2[2] + back.area_vector_m2[2],
                ];
                let rel_vec = norm(sum) / link.area_m2;
                assert!(
                    rel_vec < 1e-12,
                    "facet normals are not opposite: relative residual {rel_vec:e}"
                );
                if link.face != face {
                    checked_cross_face += 1;
                }
            }
        }
        assert!(checked_cross_face > 0, "seam facets must be exercised");
    }

    #[test]
    fn internal_flux_contributions_cancel_pairwise() {
        // Summing a constant field's flux over every cell must leave only
        // the shell's outer and inner boundaries, which for a closed shell
        // cancel too. Any leak means an interface is counted with the wrong
        // sign or the wrong area somewhere.
        let s = ShellSpec::new(6, 3, R_CMB, R_TOP).expect("valid");
        let field = [0.37, -0.81, 0.44];
        let mut total = 0.0;
        let mut scale = 0.0;
        for e in 0..s.cells() as u64 {
            let (face, i, j, k) = s.fijk(e).expect("in range");
            for facet in Facet::ALL {
                let a = s
                    .facet_area_vector_m2(face, i, j, k, facet)
                    .expect("in range");
                total += dot(field, a);
                scale += norm(a);
            }
        }
        let rel = total.abs() / scale;
        assert!(
            rel < 1e-12,
            "constant-field flux does not cancel: rel {rel:e}"
        );
    }

    #[test]
    fn adr_026_declared_shapes_have_the_declared_element_counts() {
        // Guards the numbers written into ADR 026, the master spec section
        // 58 and the build brief: if the ladder changes, the docs are wrong
        // and this fails rather than the docs quietly drifting.
        let preview = ShellSpec::preview(R_CMB, R_TOP).expect("preview valid");
        assert_eq!(preview.cells(), 221_184);
        assert_eq!(preview.ladder_level(), Some(0));
        let standard = ShellSpec::standard(R_CMB, R_TOP).expect("standard valid");
        assert_eq!(standard.cells(), 1_327_104);
        assert_eq!(standard.ladder_level(), Some(1));
    }

    #[test]
    fn construction_rejects_invalid_shapes_and_radii() {
        assert!(ShellSpec::new(1, 4, R_CMB, R_TOP).is_err(), "n too small");
        assert!(ShellSpec::new(7, 4, R_CMB, R_TOP).is_err(), "n odd");
        assert!(ShellSpec::new(8, 0, R_CMB, R_TOP).is_err(), "zero layers");
        assert!(ShellSpec::new(8, 25, R_CMB, R_TOP).is_err(), "layers > 24");
        assert!(
            ShellSpec::new(8, 4, R_TOP, R_CMB).is_err(),
            "inverted radii"
        );
        assert!(ShellSpec::new(8, 4, 0.0, R_TOP).is_err(), "zero bottom");
        assert!(
            ShellSpec::new(8, 4, R_CMB, f64::NAN).is_err(),
            "nonfinite radius"
        );
    }

    #[test]
    fn domain_declaration_is_a_volume_domain_with_layered_digests() {
        let s = mini();
        let a = s.domain_declaration("mantle.primary").expect("declares");
        assert_eq!(a.kind, SpatialDomainKind::Volume3D);
        assert_eq!(
            a.discretization.coordinate_system,
            CoordinateSystem::CubedSphere
        );
        assert_eq!(a.discretization.element_count, s.cells() as u64);
        assert_eq!(a.discretization.resolution_level, u32::MAX, "off-ladder");
        assert!(!a.discretization.topology_digest.is_empty());
        assert_ne!(
            a.discretization.topology_digest,
            a.discretization.metric_digest
        );

        // Same topology, different planet radius: topology digest holds,
        // metric digest moves.
        let bigger = ShellSpec::new(8, 4, R_CMB, R_TOP * 1.1).expect("valid");
        let b = bigger
            .domain_declaration("mantle.primary")
            .expect("declares");
        assert_eq!(
            a.discretization.topology_digest,
            b.discretization.topology_digest
        );
        assert_ne!(
            a.discretization.metric_digest,
            b.discretization.metric_digest
        );
    }

    #[test]
    fn deserialization_validates_and_rejects_mislabelled_levels() {
        let s = ShellSpec::preview(R_CMB, R_TOP).expect("valid");
        let json = serde_json::to_string(&s).expect("serializes");
        let back: ShellSpec = serde_json::from_str(&json).expect("round trips");
        assert_eq!(back, s);

        // An L0 claim on a shape that is not the L0 shape must fail.
        let bad =
            format!(r#"{{"n":8,"layers":4,"r_bottom_m":{R_CMB},"r_top_m":{R_TOP},"level":0}}"#);
        assert!(serde_json::from_str::<ShellSpec>(&bad).is_err());

        // So must an invalid shape, even without a level claim.
        let bad = format!(r#"{{"n":7,"layers":4,"r_bottom_m":{R_CMB},"r_top_m":{R_TOP}}}"#);
        assert!(serde_json::from_str::<ShellSpec>(&bad).is_err());
    }

    #[test]
    fn spatial_discretization_surface_agrees_with_the_native_one() {
        let s = mini();
        assert_eq!(s.element_count(), s.cells() as u64);
        let e = ElementId(s.flat(2, 3, 4, 1));
        let (f, i, j, k) = s.fijk(e.0).expect("in range");
        assert_eq!((f, i, j, k), (2, 3, 4, 1));

        let c = SpatialDiscretization::centroid(&s, e).expect("centroid");
        let r = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        assert!((r - s.center_radius_m(1).expect("k=1")).abs() < 1e-6);

        let links = SpatialDiscretization::neighbors(&s, e).expect("neighbours");
        assert_eq!(links.len(), 6, "interior cell: 4 lateral + 2 radial");
        for l in &links {
            assert!(l.distance_m.is_finite() && l.distance_m > 0.0);
        }

        // Distance is symmetric and matches the straight-line separation.
        let other = links[0].element;
        let d1 = MetricSpace::distance_m(&s, e, other).expect("distance");
        let d2 = MetricSpace::distance_m(&s, other, e).expect("distance");
        assert!((d1 - d2).abs() < 1e-9);
    }

    #[test]
    fn boundary_layers_have_one_radial_neighbor() {
        let s = mini();
        let bottom = ElementId(s.flat(0, 0, 0, 0));
        let top = ElementId(s.flat(0, 0, 0, s.layers() - 1));
        assert_eq!(
            SpatialDiscretization::neighbors(&s, bottom)
                .expect("cmb")
                .len(),
            5
        );
        assert_eq!(
            SpatialDiscretization::neighbors(&s, top)
                .expect("top")
                .len(),
            5
        );
    }

    #[test]
    fn out_of_range_elements_fail_loud() {
        let s = mini();
        let bad = ElementId(s.cells() as u64);
        assert!(SpatialDiscretization::neighbors(&s, bad).is_err());
        assert!(SpatialDiscretization::centroid(&s, bad).is_err());
        assert!(VolumeMetrics::volume_m3(&s, bad).is_err());
        assert!(FaceTopology::faces(&s, bad).is_err());
    }
}
