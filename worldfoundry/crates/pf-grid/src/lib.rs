//! pf-grid: spherical lat/lon grid geometry and topology (B1 rev 2 §G).
//!
//! Owns cell bounds/centers/areas, neighbor topology including antimeridian
//! wrap and the declared pole-cap policy, great-circle distances, and the
//! grid level ladder. Area invariant (sum = 4*pi*R^2) is validated at
//! construction — fail-loud, never papered over.
//!
//! Pole-cap policy (declared, tested): the polar rows have NO neighbors
//! beyond their own row — along-row east/west adjacency remains, plus the
//! three neighbors toward the equator. A cap row is therefore a normal part
//! of the drainage graph (water can move along-row and pole-to-equator) and
//! is NEVER an invented barrier or spurious outlet.
//!
//! The mantle volume domain lives in [`shell`]: a cubed-sphere spherical
//! shell (MA of brief rev 5 / ADR 026). It is a second backend of the same
//! contracts, never a redefinition of them.

pub mod shell;
pub mod shell_ops;

use serde::{Deserialize, Serialize};

#[derive(Debug, thiserror::Error)]
pub enum GridError {
    #[error("grid nlat={nlat} nlon={nlon} invalid: {why}")]
    Invalid {
        nlat: usize,
        nlon: usize,
        why: &'static str,
    },
    #[error("cell index out of range: i={i} j={j} on nlat={nlat} nlon={nlon}")]
    Index {
        i: usize,
        j: usize,
        nlat: usize,
        nlon: usize,
    },
}

/// Fields are PRIVATE: every construction path (constructors AND
/// deserialization) validates dimensions, radius, and the total-area
/// invariant. There is no way to hold an invalid GridSpec.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct GridSpec {
    nlat: usize,
    nlon: usize,
    radius_m: f64,
    /// Ladder level when constructed via `level()`; None for off-ladder
    /// grids (mini/test) — recorded explicitly in the domain descriptor as
    /// resolution_level u32::MAX, never silently defaulted.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    level: Option<u8>,
}

impl<'de> Deserialize<'de> for GridSpec {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct Raw {
            nlat: usize,
            nlon: usize,
            radius_m: f64,
            #[serde(default)]
            level: Option<u8>,
        }
        let raw = Raw::deserialize(d)?;
        let mut spec = GridSpec::new(raw.nlat, raw.nlon, raw.radius_m)
            .map_err(|e| serde::de::Error::custom(format!("invalid GridSpec: {e}")))?;
        if let Some(level) = raw.level {
            // A claimed ladder level must match the physical dimensions
            // exactly — L3 claims on an L0 shape are contract violations,
            // not silent mislabels (review hardening).
            if level > 3 || spec.nlat != BASE_NLAT << level || spec.nlon != BASE_NLON << level {
                return Err(serde::de::Error::custom(format!(
                    "invalid GridSpec: level={level} inconsistent with {}x{}",
                    spec.nlat, spec.nlon
                )));
            }
        }
        spec.level = raw.level;
        Ok(spec)
    }
}

/// B1 level ladder: L0 = 288x144, each level doubles both axes.
pub const BASE_NLON: usize = 288;
pub const BASE_NLAT: usize = 144;

impl GridSpec {
    /// Level ladder constructor: level 0..=3 => 288x144 .. 2304x1152.
    pub fn level(level: u8, radius_m: f64) -> Result<Self, GridError> {
        if level > 3 {
            return Err(GridError::Invalid {
                nlat: 0,
                nlon: 0,
                why: "level cap is 3",
            });
        }
        let f = 1usize << level;
        let mut spec = GridSpec::new(BASE_NLAT * f, BASE_NLON * f, radius_m)?;
        spec.level = Some(level);
        Ok(spec)
    }

    pub fn new(nlat: usize, nlon: usize, radius_m: f64) -> Result<Self, GridError> {
        if nlat < 4 || nlat % 2 != 0 {
            return Err(GridError::Invalid {
                nlat,
                nlon,
                why: "nlat must be even and >= 4",
            });
        }
        if nlon < 8 || nlon % 2 != 0 {
            return Err(GridError::Invalid {
                nlat,
                nlon,
                why: "nlon must be even and >= 8",
            });
        }
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(GridError::Invalid {
                nlat,
                nlon,
                why: "radius must be finite positive",
            });
        }
        let spec = Self {
            nlat,
            nlon,
            radius_m,
            level: None,
        };
        // Area invariant is enforced AT CONSTRUCTION, not just tested:
        // no code path can hold a grid whose cells do not tile the sphere.
        spec.validate_total_area()?;
        Ok(spec)
    }

    pub fn cells(&self) -> usize {
        self.nlat * self.nlon
    }

    // Accessors: fields are private so external crates cannot construct or
    // mutate a grid bypassing validation.
    pub fn nlat(&self) -> usize {
        self.nlat
    }
    pub fn nlon(&self) -> usize {
        self.nlon
    }
    pub fn radius_m(&self) -> f64 {
        self.radius_m
    }
    pub fn ladder_level(&self) -> Option<u8> {
        self.level
    }

    /// Cell center in degrees: latitude from -90 (row 0, south) to +90.
    pub fn center_lat_deg(&self, i: usize) -> Result<f64, GridError> {
        if i >= self.nlat {
            return Err(GridError::Index {
                i,
                j: 0,
                nlat: self.nlat,
                nlon: self.nlon,
            });
        }
        Ok(-90.0 + (i as f64 + 0.5) * 180.0 / self.nlat as f64)
    }

    /// Cell center longitude in degrees: -180 (col 0) to +180, exclusive.
    pub fn center_lon_deg(&self, j: usize) -> Result<f64, GridError> {
        if j >= self.nlon {
            return Err(GridError::Index {
                i: 0,
                j,
                nlat: self.nlat,
                nlon: self.nlon,
            });
        }
        Ok(-180.0 + (j as f64 + 0.5) * 360.0 / self.nlon as f64)
    }

    /// Spherical rectangle area in m^2: R^2 * dlon * (sin lat_top - sin lat_bot).
    pub fn cell_area_m2(&self, i: usize) -> Result<f64, GridError> {
        if i >= self.nlat {
            return Err(GridError::Index {
                i,
                j: 0,
                nlat: self.nlat,
                nlon: self.nlon,
            });
        }
        let dlat = std::f64::consts::PI / self.nlat as f64;
        let dlon = 2.0 * std::f64::consts::PI / self.nlon as f64;
        let lat_bot = -std::f64::consts::FRAC_PI_2 + i as f64 * dlat;
        let lat_top = lat_bot + dlat;
        Ok(self.radius_m * self.radius_m * dlon * (lat_top.sin() - lat_bot.sin()))
    }

    /// Total area by summation (all nlon cells per row band); compared
    /// against 4*pi*R^2 — this invariant caught a missing nlon factor once
    /// already, it stays.
    pub fn total_area_m2(&self) -> f64 {
        (0..self.nlat)
            .map(|i| self.cell_area_m2(i).unwrap() * self.nlon as f64)
            .sum()
    }

    /// Fail-loud area invariant (§G): relative error must be < 1e-9.
    pub fn validate_total_area(&self) -> Result<(), GridError> {
        let expected = 4.0 * std::f64::consts::PI * self.radius_m * self.radius_m;
        let got = self.total_area_m2();
        let rel = ((got - expected) / expected).abs();
        if !rel.is_finite() || rel > 1.0e-9 {
            return Err(GridError::Invalid {
                nlat: self.nlat,
                nlon: self.nlon,
                why: "total cell area != 4*pi*R^2 beyond 1e-9 relative",
            });
        }
        Ok(())
    }

    /// Great-circle distance between two cell centers (meters).
    pub fn center_distance_m(
        &self,
        i0: usize,
        j0: usize,
        i1: usize,
        j1: usize,
    ) -> Result<f64, GridError> {
        let lat0 = self.center_lat_deg(i0)?.to_radians();
        let lon0 = self.center_lon_deg(j0)?.to_radians();
        let lat1 = self.center_lat_deg(i1)?.to_radians();
        let lon1 = self.center_lon_deg(j1)?.to_radians();
        let dlat = lat1 - lat0;
        let dlon = lon1 - lon0;
        let a = (dlat / 2.0).sin().powi(2) + lat0.cos() * lat1.cos() * (dlon / 2.0).sin().powi(2);
        // Clamp guards fp excursions just above 1; result stays >= 0.
        let a = a.clamp(0.0, 1.0);
        let c = 2.0 * a.sqrt().asin();
        Ok(self.radius_m * c)
    }
}

/// One D8-style neighbor: grid offset (with antimeridian wrap on dj), plus
/// the center-to-center distance so downstream stages use real metrics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Neighbor {
    pub di: i32,
    pub dj: i32,
    pub distance_m: f64,
}

impl GridSpec {
    /// 8-neighborhood with antimeridian wrap; declared pole-cap policy:
    /// polar rows keep along-row E/W and the three equatorward neighbors.
    /// Returns (target i, target j, Neighbor) triples; j wraps mod nlon.
    pub fn neighbors(
        &self,
        i: usize,
        j: usize,
    ) -> Result<Vec<(usize, usize, Neighbor)>, GridError> {
        if i >= self.nlat || j >= self.nlon {
            return Err(GridError::Index {
                i,
                j,
                nlat: self.nlat,
                nlon: self.nlon,
            });
        }
        // SEED-PERF (lane SP, 2026-10-13): the distances come from a
        // per-row cache of exactly these `center_distance_m` values (same
        // calls, same bits), filled the first time a row is asked for.
        if let Some(row) = nb_cache::row(self, i) {
            let d = &row[j];
            let i0 = i as i32;
            let mut out = Vec::with_capacity(8);
            let mut slot = 0;
            for di in [-1, 0, 1] {
                let ti = i0 + di;
                if ti < 0 || ti >= self.nlat as i32 {
                    continue;
                }
                for dj in [-1, 0, 1] {
                    if di == 0 && dj == 0 {
                        continue;
                    }
                    let tj = ((j as i32 + dj).rem_euclid(self.nlon as i32)) as usize;
                    out.push((ti as usize, tj, Neighbor { di, dj, distance_m: d[slot] }));
                    slot += 1;
                }
            }
            return Ok(out);
        }
        self.neighbors_uncached(i, j)
    }

    /// The direct evaluation behind `neighbors` (and its cache rows).
    fn neighbors_uncached(
        &self,
        i: usize,
        j: usize,
    ) -> Result<Vec<(usize, usize, Neighbor)>, GridError> {
        let i0 = i as i32;
        let mut out = Vec::with_capacity(8);
        for di in [-1, 0, 1] {
            // Cap policy: no neighbor rows beyond the pole rows.
            let ti = i0 + di;
            if ti < 0 || ti >= self.nlat as i32 {
                continue;
            }
            for dj in [-1, 0, 1] {
                if di == 0 && dj == 0 {
                    continue;
                }
                let tj = ((j as i32 + dj).rem_euclid(self.nlon as i32)) as usize;
                let ti_u = ti as usize;
                let distance_m = self.center_distance_m(i, j, ti_u, tj)?;
                out.push((ti_u, tj, Neighbor { di, dj, distance_m }));
            }
        }
        Ok(out)
    }
}

impl GridSpec {
    /// Flat (element) index for cell (i, j): row-major C-order.
    pub fn flat(&self, i: usize, j: usize) -> u64 {
        (i * self.nlon + j) as u64
    }

    /// Decode a flat element index to (i, j); fails loud out of range.
    pub fn ij(&self, element: u64) -> Result<(usize, usize), pf_contracts::spatial::SpatialError> {
        let count = self.cells() as u64;
        let e = usize::try_from(element).map_err(|_| pf_contracts::spatial::SpatialError {
            element,
            element_count: count,
        })?;
        if e >= self.cells() {
            return Err(pf_contracts::spatial::SpatialError {
                element,
                element_count: count,
            });
        }
        Ok((e / self.nlon, e % self.nlon))
    }

    /// Domain declaration for this grid as a Surface2D lat/lon domain, with
    /// computed topology + metric digests (review round 3 identity split).
    pub fn domain_declaration(
        &self,
        domain_id: &str,
    ) -> Result<pf_contracts::spatial::DomainDeclaration, pf_contracts::ContractError> {
        use pf_contracts::spatial::{
            CoordinateSystem, DatumDescriptor, DiscretizationDescriptor, DomainDeclaration,
            SpatialDomainKind,
        };
        let resolution_level = self.level.map(u32::from).unwrap_or(u32::MAX);
        let mut decl = DomainDeclaration {
            domain_id: domain_id.into(),
            kind: SpatialDomainKind::Surface2D,
            discretization: DiscretizationDescriptor {
                topology_id: format!("surface.latlon.{}x{}", self.nlat, self.nlon),
                backend: "latlon-v0".into(),
                coordinate_system: CoordinateSystem::SphericalLatLon,
                topology_digest: String::new(),
                metric_digest: String::new(),
                element_count: self.cells() as u64,
                resolution_level,
                transform_conventions: format!(
                    "element=i*nlon+j (C-order); lat row 0 = -90deg; lon wraps at antimeridian; pole rows have 5 neighbors (declared cap policy)"
                ),
            },
            datum: DatumDescriptor { center: [0.0, 0.0, 0.0], radius_m: self.radius_m },
        };
        decl.compute_digests()?;
        Ok(decl)
    }
}

impl pf_contracts::spatial::SpatialDiscretization for GridSpec {
    fn domain(&self) -> pf_contracts::spatial::DomainDeclaration {
        // domain_id is semantic; the grid backend supplies a conventional
        // default — artifacts may declare their own via domain_declaration().
        self.domain_declaration("surface.primary")
            .expect("GridSpec invariants validated at construction")
    }

    fn element_count(&self) -> u64 {
        self.cells() as u64
    }

    fn neighbors(
        &self,
        element: pf_contracts::spatial::ElementId,
    ) -> Result<Vec<pf_contracts::spatial::NeighborLink>, pf_contracts::spatial::SpatialError> {
        let (i, j) = self.ij(element.0)?;
        let links = self
            .neighbors(i, j)
            .map_err(|_| pf_contracts::spatial::SpatialError {
                element: element.0,
                element_count: self.cells() as u64,
            })?
            .into_iter()
            .map(|(ti, tj, n)| pf_contracts::spatial::NeighborLink {
                element: pf_contracts::spatial::ElementId(self.flat(ti, tj)),
                distance_m: n.distance_m,
            })
            .collect::<Vec<_>>();
        Ok(links)
    }

    fn centroid(
        &self,
        element: pf_contracts::spatial::ElementId,
    ) -> Result<pf_contracts::spatial::Position, pf_contracts::spatial::SpatialError> {
        let err = || pf_contracts::spatial::SpatialError {
            element: element.0,
            element_count: self.cells() as u64,
        };
        let (i, j) = self.ij(element.0)?;
        let lat = self.center_lat_deg(i).map_err(|_| err())?.to_radians();
        let lon = self.center_lon_deg(j).map_err(|_| err())?.to_radians();
        let r = self.radius_m;
        Ok([
            r * lat.cos() * lon.cos(),
            r * lat.cos() * lon.sin(),
            r * lat.sin(),
        ])
    }
}

impl pf_contracts::spatial::SurfaceMetrics for GridSpec {
    fn area_m2(
        &self,
        element: pf_contracts::spatial::ElementId,
    ) -> Result<f64, pf_contracts::spatial::SpatialError> {
        let (i, _) = self.ij(element.0)?;
        self.cell_area_m2(i)
            .map_err(|_| pf_contracts::spatial::SpatialError {
                element: element.0,
                element_count: self.cells() as u64,
            })
    }
}

impl pf_contracts::spatial::MetricSpace for GridSpec {
    fn distance_m(
        &self,
        a: pf_contracts::spatial::ElementId,
        b: pf_contracts::spatial::ElementId,
    ) -> Result<f64, pf_contracts::spatial::SpatialError> {
        let _ = self.ij(a.0)?;
        let _ = self.ij(b.0)?;
        let (i0, j0) = self.ij(a.0).unwrap();
        let (i1, j1) = self.ij(b.0).unwrap();
        self.center_distance_m(i0, j0, i1, j1)
            .map_err(|_| pf_contracts::spatial::SpatialError {
                element: a.0,
                element_count: self.cells() as u64,
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn earth() -> GridSpec {
        GridSpec::level(0, 6.371e6).unwrap()
    }

    #[test]
    fn level_ladder() {
        assert_eq!(GridSpec::level(0, 1.0).unwrap().cells(), 288 * 144);
        assert_eq!(GridSpec::level(1, 1.0).unwrap().cells(), 576 * 288);
        assert_eq!(GridSpec::level(2, 1.0).unwrap().cells(), 1152 * 576);
        assert!(GridSpec::level(4, 1.0).is_err());
    }

    #[test]
    fn total_area_matches_sphere() {
        for spec in [earth(), GridSpec::level(1, 6.371e6).unwrap()] {
            spec.validate_total_area().unwrap();
        }
    }

    #[test]
    fn midrow_has_8_neighbors_polerows_have_5() {
        let g = earth();
        assert_eq!(g.neighbors(0, 5).unwrap().len(), 5); // south cap
        assert_eq!(g.neighbors(g.nlat - 1, 5).unwrap().len(), 5); // north cap
        assert_eq!(g.neighbors(10, 100).unwrap().len(), 8);
    }

    #[test]
    fn deserialization_rejects_inconsistent_level_and_roundtrips_valid() {
        let bad = format!(
            "{{\"nlat\":{},\"nlon\":{},\"radius_m\":6371000.0,\"level\":3}}",
            BASE_NLAT, BASE_NLON
        );
        let err = serde_json::from_str::<GridSpec>(&bad).unwrap_err();
        assert!(err.to_string().contains("inconsistent"), "{err}");
        let good = format!(
            "{{\"nlat\":{},\"nlon\":{},\"radius_m\":6371000.0,\"level\":1}}",
            BASE_NLAT * 2,
            BASE_NLON * 2
        );
        let spec = serde_json::from_str::<GridSpec>(&good).unwrap();
        assert_eq!(spec.ladder_level(), Some(1));
        // and it re-serializes
        let json = serde_json::to_string(&spec).unwrap();
        let back: GridSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(back, spec);
    }

    #[test]
    fn antimeridian_wraps() {
        let g = earth();
        let ns = g.neighbors(10, 0).unwrap();
        // west of col 0 is col nlon-1 (across the antimeridian)
        assert!(ns
            .iter()
            .any(|(i, j, n)| *i == 10 && *j == g.nlon - 1 && n.dj == -1));
        let ns2 = g.neighbors(10, g.nlon - 1).unwrap();
        assert!(ns2.iter().any(|(i, j, n)| *i == 10 && *j == 0 && n.dj == 1));
    }

    #[test]
    fn cap_is_not_a_barrier_along_row_flow_exists() {
        let g = earth();
        for (i, j) in [(0usize, 3usize), (g.nlat - 1, 3)] {
            let ns = g.neighbors(i, j).unwrap();
            // along-row = same latitude row (ti == i), adjacent column
            assert!(
                ns.iter().any(|n| n.0 == i && (n.2.dj == -1 || n.2.dj == 1)),
                "no along-row neighbor at {i},{j}"
            );
            assert!(
                ns.iter().any(|n| n.0 != i),
                "no off-row neighbor at {i},{j}"
            );
        }
    }

    #[test]
    fn great_circle_quarter_equator() {
        let g = earth();
        let i_eq = g.nlat / 2;
        // quarter of the equator = pi*R/2
        let j0 = 0;
        let jq = g.nlon / 4;
        let d = g.center_distance_m(i_eq, j0, i_eq, jq).unwrap();
        let expect = std::f64::consts::PI * g.radius_m / 2.0;
        let cell = 2.0 * std::f64::consts::PI * g.radius_m / g.nlon as f64;
        assert!(
            (d - expect).abs() < cell,
            "d={d} expect={expect} cell={cell}"
        );
    }

    #[test]
    fn areas_are_positive_and_poles_are_small() {
        let g = earth();
        assert!(g.cell_area_m2(0).unwrap() > 0.0);
        assert!(g.cell_area_m2(g.nlat / 2).unwrap() > g.cell_area_m2(0).unwrap());
    }

    #[test]
    fn domain_declaration_digests_are_layered_and_stable() {
        let g = earth();
        let d1 = g.domain_declaration("surface.primary").unwrap();
        let d2 = g.domain_declaration("surface.primary").unwrap();
        assert_eq!(d1, d2);
        assert!(!d1.discretization.topology_digest.is_empty());
        assert_ne!(
            d1.discretization.topology_digest,
            d1.discretization.metric_digest
        );
        // Same connectivity, different radius => same topology, different metric.
        let big = GridSpec::level(0, 8.0e6).unwrap();
        let d3 = big.domain_declaration("surface.primary").unwrap();
        assert_eq!(
            d1.discretization.topology_digest,
            d3.discretization.topology_digest
        );
        assert_ne!(
            d1.discretization.metric_digest,
            d3.discretization.metric_digest
        );
        assert_eq!(d1.discretization.resolution_level, 0);
    }

    #[test]
    fn off_ladder_grid_records_max_resolution_level() {
        let mini = GridSpec::new(72, 36, 6.371e6).unwrap();
        let d = mini.domain_declaration("surface.mini").unwrap();
        assert_eq!(d.discretization.resolution_level, u32::MAX);
    }

    #[test]
    fn trait_topology_matches_inherent_counts() {
        use pf_contracts::spatial::{ElementId, SpatialDiscretization};
        let g = earth();
        let mid = ElementId(g.flat(10, 100));
        let pole = ElementId(g.flat(0, 5));
        assert_eq!(SpatialDiscretization::neighbors(&g, mid).unwrap().len(), 8);
        assert_eq!(SpatialDiscretization::neighbors(&g, pole).unwrap().len(), 5);
        // centroid lies on the sphere
        let c = g.centroid(mid).unwrap();
        let norm = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        assert!(
            (norm - g.radius_m).abs() < 1.0,
            "centroid off sphere: {norm}"
        );
    }

    #[test]
    fn surface_metrics_total_area_via_trait() {
        use pf_contracts::spatial::{SpatialDiscretization, SurfaceMetrics};
        let g = earth();
        let count = g.element_count();
        let total: f64 = (0..count)
            .map(|e| g.area_m2(pf_contracts::spatial::ElementId(e)).unwrap())
            .sum();
        let expect = 4.0 * std::f64::consts::PI * g.radius_m * g.radius_m;
        assert!(((total - expect) / expect).abs() < 1.0e-9);
    }
}

/// SEED-PERF (lane SP, 2026-10-13): process-wide cache of the 8-neighbour
/// centre distances, per grid (nlat, nlon, radius bits) and per row, so
/// `GridSpec::neighbors` stops re-evaluating 8 haversines on every call.
/// Values are produced by the uncached `neighbors` itself (bit identical);
/// a row whose evaluation fails is not cached (the caller then takes the
/// uncached path and its exact error). Bounded: at most 16 grids of at most
/// 2^23 cells each; other grids are simply not cached.
mod nb_cache {
    use super::GridSpec;
    use std::sync::{Arc, OnceLock, RwLock};

    type Key = (usize, usize, u64);
    pub(super) struct Table {
        rows: Vec<OnceLock<Box<[[f64; 8]]>>>,
    }
    static REG: RwLock<Vec<(Key, Arc<Table>)>> = RwLock::new(Vec::new());
    const MAX_GRIDS: usize = 16;
    const MAX_CELLS: usize = 1 << 23;

    thread_local! {
        static LAST: std::cell::RefCell<Option<(Key, Arc<Table>)>> = const { std::cell::RefCell::new(None) };
    }

    fn table(g: &GridSpec) -> Option<Arc<Table>> {
        let key = (g.nlat, g.nlon, g.radius_m.to_bits());
        if let Some(t) = LAST.with(|l| l.borrow().as_ref().filter(|(k, _)| *k == key).map(|(_, t)| t.clone())) {
            return Some(t);
        }
        let found = REG.read().ok()?.iter().find(|(k, _)| *k == key).map(|(_, t)| t.clone());
        let t = match found {
            Some(t) => t,
            None => {
                if g.nlat.saturating_mul(g.nlon) > MAX_CELLS {
                    return None;
                }
                let mut reg = REG.write().ok()?;
                if let Some((_, t)) = reg.iter().find(|(k, _)| *k == key) {
                    t.clone()
                } else {
                    if reg.len() >= MAX_GRIDS {
                        return None;
                    }
                    let t = Arc::new(Table { rows: (0..g.nlat).map(|_| OnceLock::new()).collect() });
                    reg.push((key, t.clone()));
                    t
                }
            }
        };
        LAST.with(|l| *l.borrow_mut() = Some((key, t.clone())));
        Some(t)
    }

    /// Row `i`'s distances per column, in `neighbors` order (slot k = k-th
    /// emitted neighbour), or None (not cacheable / evaluation failed).
    pub(super) fn row(g: &GridSpec, i: usize) -> Option<RowRef> {
        let t = table(g)?;
        let cell = t.rows.get(i)?;
        if cell.get().is_none() {
            let mut v = Vec::with_capacity(g.nlon);
            for j in 0..g.nlon {
                let nb = g.neighbors_uncached(i, j).ok()?;
                let mut d = [0.0f64; 8];
                for (k, (_, _, n)) in nb.iter().enumerate() {
                    d[k] = n.distance_m;
                }
                v.push(d);
            }
            let _ = cell.set(v.into_boxed_slice());
        }
        Some(RowRef { t, i })
    }

    pub(super) struct RowRef {
        t: Arc<Table>,
        i: usize,
    }
    impl std::ops::Index<usize> for RowRef {
        type Output = [f64; 8];
        fn index(&self, j: usize) -> &[f64; 8] {
            &self.t.rows[self.i].get().expect("row filled")[j]
        }
    }
}

#[cfg(test)]
mod nb_cache_tests {
    use super::*;

    /// SEED-PERF (lane SP): the cached `neighbors` is the uncached one, bit
    /// for bit, on ladder and mini grids (radius part of the key), from many threads.
    #[test]
    fn cached_neighbors_are_the_uncached_ones() {
        let grids = [
            GridSpec::new(4, 8, 6.371e6).unwrap(),
            GridSpec::new(8, 16, 6.371e6).unwrap(),
            GridSpec::new(18, 36, 1.0e3).unwrap(),
            GridSpec::level(0, 6.371e6).unwrap(),
        ];
        std::thread::scope(|s| {
            for _ in 0..4 {
                s.spawn(|| {
                    for g in &grids {
                        for i in 0..g.nlat {
                            for j in 0..g.nlon {
                                let a = g.neighbors(i, j).unwrap();
                                let b = g.neighbors_uncached(i, j).unwrap();
                                assert_eq!(a.len(), b.len());
                                for (x, y) in a.iter().zip(&b) {
                                    assert_eq!((x.0, x.1, x.2.di, x.2.dj), (y.0, y.1, y.2.di, y.2.dj));
                                    assert_eq!(x.2.distance_m.to_bits(), y.2.distance_m.to_bits());
                                }
                            }
                        }
                        assert!(g.neighbors(g.nlat, 0).is_err());
                    }
                });
            }
        });
    }
}
