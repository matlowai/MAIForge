//! Spatial topology as an explicit, described contract (B1 rev 3).
//!
//! THE SERIALIZED CONTRACT AND THE RUST TRAIT ARE DIFFERENT THINGS.
//! `DomainDeclaration` (with its digests) is the durable inter-module/data
//! contract persisted in every artifact manifest. `SpatialDiscretization`
//! and the capability traits are Rust-side conveniences we may rewrite any
//! number of times without changing artifacts on disk.
//!
//! Identity is split three ways (review round 3):
//! - `domain_id` — logical/semantic identity ("surface.primary"); stable
//!   across discretization changes; fields reference THIS.
//! - `topology_id` / `topology_digest` — discretization/connectivity
//!   identity ("surface.latlon.L0"; two grids with identical connectivity
//!   share these).
//! - `metric_digest` — physical metric identity (datum/radius-dependent
//!   areas and distances; differs when radii differ).
//!
//! All persisted counts are fixed-width (u64/u32/u8), never `usize`.

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Persisted contract

/// Logical domain kind. Volume discretizations never masquerade as surface
/// grids; capability traits below carry the metric each kind actually has.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SpatialDomainKind {
    Surface2D,
    Volume3D,
    Curve1D,
    PointSet,
}

/// Keep the coordinate-system enum open enough for the planned backend
/// family; a descriptor string can replace it later without artifact breaks
/// only via an additive v1 (A2.1 policy).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CoordinateSystem {
    SphericalLatLon,
    CubedSphere,
    GeocentricCartesian,
    CartesianEnu,
    UnstructuredSphericalShell,
}

/// How a domain is discretized: connectivity identity + conventions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DiscretizationDescriptor {
    /// Discretization identity, e.g. "surface.latlon.L0".
    pub topology_id: String,
    /// Implementation family/version, e.g. "latlon-v0". Records the BACKEND;
    /// consumers rely on this descriptor, never backend internals.
    pub backend: String,
    pub coordinate_system: CoordinateSystem,
    /// blake3 over {topology_id, coordinate_system, element_count,
    /// resolution_level, transform_conventions} — backend INDEPENDENT (the
    /// backend string is excluded, so identical descriptors across backend
    /// versions share identity). B1 semantics: backend-independent
    /// topology-DESCRIPTOR identity — not yet a literal connectivity hash;
    /// unstructured-mesh backends will fold connectivity blob identity into
    /// this when they arrive.
    pub topology_digest: String,
    /// blake3 over {topology_digest, datum}. Physical metric identity.
    pub metric_digest: String,
    pub element_count: u64,
    /// Backend ladder level; u32::MAX = explicitly off-ladder (mini/test
    /// grids) — recorded, not hidden.
    pub resolution_level: u32,
    /// Indexing/axes conventions, human+machine readable.
    pub transform_conventions: String,
}

/// Planetary datum: geocentric frame origin + radius.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DatumDescriptor {
    pub center: [f64; 3],
    pub radius_m: f64,
}

/// A declared spatial domain: the unit fields bind to by `domain_id`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DomainDeclaration {
    /// Logical/semantic identity, e.g. "surface.primary". Stable across
    /// discretization upgrades; `FieldDescriptor::domain` points here.
    pub domain_id: String,
    pub kind: SpatialDomainKind,
    pub discretization: DiscretizationDescriptor,
    pub datum: DatumDescriptor,
}

impl DomainDeclaration {
    /// Fill in empty digests (computed, never caller-invented). Digests
    /// exclude themselves; topology excludes metric/datum; metric includes
    /// topology + datum.
    pub fn compute_digests(&mut self) -> Result<(), crate::ContractError> {
        if self.discretization.topology_digest.is_empty() {
            let mut topo = self.discretization.clone();
            topo.topology_digest = String::new();
            topo.metric_digest = String::new();
            // Backend deliberately excluded: this is backend-independent
            // topology-descriptor identity (see field docs).
            topo.backend = String::new();
            self.discretization.topology_digest = crate::artifact::digest_canonical(&topo)?;
        }
        if self.discretization.metric_digest.is_empty() {
            let metric = MetricIdentity {
                topology_digest: self.discretization.topology_digest.clone(),
                datum: self.datum.clone(),
            };
            self.discretization.metric_digest = crate::artifact::digest_canonical(&metric)?;
        }
        Ok(())
    }
}

#[derive(Serialize)]
struct MetricIdentity {
    topology_digest: String,
    datum: DatumDescriptor,
}

// ---------------------------------------------------------------------------
// Rust-side convenience interfaces (rewritable; NOT the durable contract)

/// Element identifier inside one domain. Fixed-width in persisted form.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct ElementId(pub u64);

/// Position in the datum's geocentric frame, meters.
pub type Position = [f64; 3];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NeighborLink {
    pub element: ElementId,
    pub distance_m: f64,
}

#[derive(Debug, thiserror::Error)]
#[error("element {element} out of range (element_count={element_count})")]
pub struct SpatialError {
    pub element: u64,
    pub element_count: u64,
}

/// Neutral core interface. No surface assumptions: volume/curve/point
/// backends implement this without masquerading.
pub trait SpatialDiscretization {
    fn domain(&self) -> DomainDeclaration;
    fn element_count(&self) -> u64;
    fn neighbors(&self, element: ElementId) -> Result<Vec<NeighborLink>, SpatialError>;
    fn centroid(&self, element: ElementId) -> Result<Position, SpatialError>;
}

/// Optional capabilities — implement per kind:

pub trait SurfaceMetrics {
    fn area_m2(&self, element: ElementId) -> Result<f64, SpatialError>;
}

pub trait VolumeMetrics {
    fn volume_m3(&self, element: ElementId) -> Result<f64, SpatialError>;
}

pub trait MetricSpace {
    fn distance_m(&self, a: ElementId, b: ElementId) -> Result<f64, SpatialError>;
}

pub trait FaceTopology {
    /// Adjacent element + shared face area (m^2).
    fn faces(&self, element: ElementId) -> Result<Vec<(ElementId, f64)>, SpatialError>;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn declaration() -> DomainDeclaration {
        DomainDeclaration {
            domain_id: "surface.primary".into(),
            kind: SpatialDomainKind::Surface2D,
            discretization: DiscretizationDescriptor {
                topology_id: "surface.latlon.L0".into(),
                backend: "latlon-v0".into(),
                coordinate_system: CoordinateSystem::SphericalLatLon,
                topology_digest: String::new(),
                metric_digest: String::new(),
                element_count: 41472,
                resolution_level: 0,
                transform_conventions: "cell=i*nlon+j, C-order; lat row 0=-90deg; lon wraps".into(),
            },
            datum: DatumDescriptor {
                center: [0.0, 0.0, 0.0],
                radius_m: 6.371e6,
            },
        }
    }

    #[test]
    fn digests_are_deterministic_and_layered() {
        let mut a = declaration();
        a.compute_digests().unwrap();
        let mut b = declaration();
        b.compute_digests().unwrap();
        assert_eq!(
            a.discretization.topology_digest,
            b.discretization.topology_digest
        );
        assert_eq!(
            a.discretization.metric_digest,
            b.discretization.metric_digest
        );
        assert!(!a.discretization.topology_digest.is_empty());
        assert_ne!(
            a.discretization.topology_digest,
            a.discretization.metric_digest
        );
    }

    #[test]
    fn same_topology_different_radius_differs_only_in_metric_digest() {
        let mut small = declaration();
        small.compute_digests().unwrap();
        let mut big = declaration();
        big.datum.radius_m = 8.0e6;
        big.compute_digests().unwrap();
        assert_eq!(
            small.discretization.topology_digest,
            big.discretization.topology_digest
        );
        assert_ne!(
            small.discretization.metric_digest,
            big.discretization.metric_digest
        );
    }

    #[test]
    fn topology_digest_ignores_backend_string() {
        let mut a = declaration();
        a.compute_digests().unwrap();
        let mut b = declaration();
        b.discretization.backend = "latlon-v9-fork".into();
        b.compute_digests().unwrap();
        assert_eq!(
            a.discretization.topology_digest,
            b.discretization.topology_digest
        );
    }
}
