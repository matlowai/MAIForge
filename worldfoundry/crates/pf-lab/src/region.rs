//! Regional extraction for the settlement seam (ADR 032 deliverable 4):
//! a lat/lon centre and half-extent become a district-local INTEGER MILLIMETRE
//! grid with terrain, drainage and resource potential resampled from the
//! published global artifacts, tied to their digests, with explicit
//! `unavailable` states for what this tier does not supply.
//!
//! Frame: ADR 029 section 3 (district-local integer MILLIMETRES, axes
//! east/north/up, translation-only, recorded origin), shaped as a superset
//! of the client tile fixture `vz.terrain-tile/v1`: local tangent
//! equirectangular projection at the geodetic origin, cell centres at
//! half-spacing offsets, `order` east_fastest. Heights stay in the f32
//! `elevation_m` field (metres, on the planet datum).
//!
//! Authority: the source cell size at the origin latitude is recorded as
//! `source_cell_size_m`; anything finer is interpolation. Bilinear on the
//! source cell centres for continuous fields (exact for linear fields,
//! pinned by test), nearest source cell for categorical ones.

use pf_contracts::artifact::{
    digest_canonical, read_artifact_header, read_descriptor_bytes, digest_header,
    verify_artifact, BackendId, Tier,
};
use pf_contracts::domain::HydrologyDomainV0;
use pf_contracts::spatial::{
    CoordinateSystem, DatumDescriptor, DiscretizationDescriptor, DomainDeclaration,
    SpatialDomainKind,
};
use pf_contracts::terrain::{
    AvailabilityV0, RegionDomainV1, RegionFrameV1, RegionReachV1, RegionSourceV0, ResamplingV0,
    ResourceFamilyKindV0, ResourcePotentialDomainV0, TerrainDomainV0,
};
use pf_contracts::{units, ArtifactKind, ContractError};
use pf_grid::GridSpec;
use std::collections::BTreeMap;
use std::path::Path;
use std::time::Instant;

use crate::receipt::{f32_field, publish_typed, receipt, u32_field, u8_field, FieldBlob};

pub const REGION_BACKEND_NAME: &str = "region-extract-v0";
pub const REGION_BACKEND_VERSION: &str = "0.3.0";
pub const REGION_DOMAIN: &str = "region.local";
pub const PROJECTION: &str = "local-tangent-equirectangular";
pub const AXES: &str = "east_north_up";
pub const ORDER: &str = "east_fastest";
pub const HEIGHT_UNIT: &str = "m";
pub const SAMPLE_POSITION: &str = "cell_centre";
pub const MM_PER_M: f64 = 1000.0;

#[derive(Debug)]
pub enum RegionError {
    Contract(ContractError),
    Source(String),
    Args(String),
}

impl std::fmt::Display for RegionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RegionError::Contract(e) => write!(f, "contract: {e}"),
            RegionError::Source(w) => write!(f, "source artifact: {w}"),
            RegionError::Args(w) => write!(f, "arguments: {w}"),
        }
    }
}

impl std::error::Error for RegionError {}

impl From<ContractError> for RegionError {
    fn from(e: ContractError) -> Self {
        RegionError::Contract(e)
    }
}

#[derive(Debug, Clone)]
pub struct RegionRequest {
    pub lat_deg: f64,
    pub lon_deg: f64,
    /// Half-extent from the geodetic origin along each axis (mm).
    pub half_extent_mm: i64,
    /// Post spacing (mm); the client tile fixture uses 4000.
    pub spacing_mm: u32,
}

/// Source lat/lon grid geometry recovered from a published surface domain.
pub struct SourceGrid {
    pub grid: GridSpec,
    pub topology_id: String,
}

impl SourceGrid {
    /// Parse "surface.latlon.{nlat}x{nlon}" (pf-grid's declared format) and
    /// the datum radius. Anything else is refused: this extractor knows
    /// exactly one source discretization.
    pub fn from_domain(d: &DomainDeclaration) -> Result<Self, RegionError> {
        let id = &d.discretization.topology_id;
        let dims = id
            .strip_prefix("surface.latlon.")
            .ok_or_else(|| RegionError::Source(format!("unsupported topology {id:?}")))?;
        let (a, b) = dims
            .split_once('x')
            .ok_or_else(|| RegionError::Source(format!("cannot parse dims from {id:?}")))?;
        let nlat: usize = a
            .parse()
            .map_err(|_| RegionError::Source(format!("bad nlat in {id:?}")))?;
        let nlon: usize = b
            .parse()
            .map_err(|_| RegionError::Source(format!("bad nlon in {id:?}")))?;
        if d.discretization.coordinate_system != CoordinateSystem::SphericalLatLon {
            return Err(RegionError::Source("source is not SphericalLatLon".into()));
        }
        let grid = GridSpec::new(nlat, nlon, d.datum.radius_m)
            .map_err(|e| RegionError::Source(e.to_string()))?;
        if grid.cells() as u64 != d.discretization.element_count {
            return Err(RegionError::Source(format!(
                "element_count {} disagrees with {}x{}",
                d.discretization.element_count, nlat, nlon
            )));
        }
        Ok(Self {
            grid,
            topology_id: id.clone(),
        })
    }

    pub fn dlat_deg(&self) -> f64 {
        180.0 / self.grid.nlat() as f64
    }

    pub fn dlon_deg(&self) -> f64 {
        360.0 / self.grid.nlon() as f64
    }

    /// Fractional cell-centre coordinates of a lat/lon: (fi, fj) such that
    /// integer values sit exactly on cell centres; fj wraps.
    pub fn fractional_index(&self, lat_deg: f64, lon_deg: f64) -> (f64, f64) {
        let fi = (lat_deg + 90.0) / self.dlat_deg() - 0.5;
        let fj = (lon_deg + 180.0) / self.dlon_deg() - 0.5;
        (fi, fj)
    }

    /// Bilinear sample on cell centres (latitude clamped to the pole rows,
    /// longitude periodic). Exact for fields linear in (fi, fj).
    pub fn bilinear(&self, field: &[f32], lat_deg: f64, lon_deg: f64) -> f64 {
        let nlat = self.grid.nlat() as i64;
        let nlon = self.grid.nlon() as i64;
        let (fi, fj) = self.fractional_index(lat_deg, lon_deg);
        let i0 = fi.floor();
        let j0 = fj.floor();
        let ti = fi - i0;
        let tj = fj - j0;
        let at = |i: i64, j: i64| -> f64 {
            let ic = i.clamp(0, nlat - 1) as usize;
            let jc = j.rem_euclid(nlon) as usize;
            field[ic * nlon as usize + jc] as f64
        };
        let (i0, j0) = (i0 as i64, j0 as i64);
        let v00 = at(i0, j0);
        let v01 = at(i0, j0 + 1);
        let v10 = at(i0 + 1, j0);
        let v11 = at(i0 + 1, j0 + 1);
        (1.0 - ti) * ((1.0 - tj) * v00 + tj * v01) + ti * ((1.0 - tj) * v10 + tj * v11)
    }

    /// Nearest source cell (categorical fields).
    pub fn nearest(&self, lat_deg: f64, lon_deg: f64) -> usize {
        let nlat = self.grid.nlat() as i64;
        let nlon = self.grid.nlon() as i64;
        let (fi, fj) = self.fractional_index(lat_deg, lon_deg);
        let i = (fi.round() as i64).clamp(0, nlat - 1) as usize;
        let j = (fj.round() as i64).rem_euclid(nlon) as usize;
        i * nlon as usize + j
    }

    /// Source cell size (m) at `lat_deg`: the LARGER of the two cell sides
    /// (north-south is constant; east-west shrinks toward the poles), so
    /// the recorded authority scale never claims the finer dimension.
    pub fn cell_size_m(&self, lat_deg: f64) -> f64 {
        let r = self.grid.radius_m();
        let ns = r * self.dlat_deg().to_radians();
        let ew = r * lat_deg.to_radians().cos() * self.dlon_deg().to_radians();
        ns.max(ew)
    }
}

/// The local frame: integer-millimetre grid centred on the origin.
pub struct LocalFrame {
    pub frame: RegionFrameV1,
    pub radius_m: f64,
}

impl LocalFrame {
    pub fn new(req: &RegionRequest, radius_m: f64) -> Result<Self, RegionError> {
        if req.spacing_mm == 0 || req.half_extent_mm <= 0 {
            return Err(RegionError::Args(
                "spacing and half extent must be > 0".into(),
            ));
        }
        if !(-90.0..=90.0).contains(&req.lat_deg) || !(-180.0..=180.0).contains(&req.lon_deg) {
            return Err(RegionError::Args("centre outside lat/lon range".into()));
        }
        if req.lat_deg.abs() > 89.0 {
            return Err(RegionError::Args(
                "local tangent frame is not defined within 1 degree of a pole at this tier".into(),
            ));
        }
        if req.spacing_mm % 2 != 0 || !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(RegionError::Args(
                "even millimetre spacing and positive finite radius required".into(),
            ));
        }
        let width = req
            .half_extent_mm
            .checked_mul(2)
            .ok_or_else(|| RegionError::Args("extent overflow".into()))? as u64;
        if width % req.spacing_mm as u64 != 0 {
            return Err(RegionError::Args(
                "extent must be an exact multiple of spacing; no hidden footprint expansion".into(),
            ));
        }
        let n = width / req.spacing_mm as u64;
        // Includes all output arrays, projection scratch and publication copies.
        if n < 2
            || n > 4096
            || n.checked_mul(n)
                .and_then(|x| x.checked_mul(256))
                .map(|x| x > 512 * 1024 * 1024)
                .unwrap_or(true)
        {
            return Err(RegionError::Args(
                "region exceeds 512 MiB estimated working budget or has fewer than two samples"
                    .into(),
            ));
        }
        let extent_m = req.half_extent_mm as f64 / MM_PER_M;
        let dlat = (extent_m / radius_m).to_degrees();
        if req.lat_deg - dlat <= -89.0
            || req.lat_deg + dlat >= 89.0
            || extent_m >= std::f64::consts::FRAC_PI_2 * radius_m * req.lat_deg.to_radians().cos()
        {
            return Err(RegionError::Args(
                "region footprint leaves supported local projection".into(),
            ));
        }
        let n = n as u32;
        let x_min_mm = -req.half_extent_mm;
        let half_spacing = (req.spacing_mm / 2) as i64;
        Ok(Self {
            frame: RegionFrameV1 {
                height_datum: "planet_datum".into(),
                origin_lat_deg: req.lat_deg,
                origin_lon_deg: req.lon_deg,
                datum_radius_m: radius_m,
                projection: PROJECTION.into(),
                axes: AXES.into(),
                unit: "mm".into(),
                district_key: None,
                spacing_mm: req.spacing_mm,
                cells_east: n,
                cells_north: n,
                origin_mm: [x_min_mm + half_spacing, x_min_mm + half_spacing, 0],
                half_extent_mm: req.half_extent_mm,
                height_unit: HEIGHT_UNIT.into(),
                sample_position: SAMPLE_POSITION.into(),
                order: ORDER.into(),
            },
            radius_m,
        })
    }

    pub fn cells(&self) -> usize {
        (self.frame.cells_east as usize) * (self.frame.cells_north as usize)
    }

    /// South-west grid corner in local mm (origin_mm is half a spacing in).
    fn x_min_mm(&self) -> i64 {
        -self.frame.half_extent_mm
    }

    /// Centre of local cell (ix, iy) in local METRES (f64, for projection).
    pub fn cell_centre_m(&self, ix: u32, iy: u32) -> (f64, f64) {
        let spacing = self.frame.spacing_mm as i64;
        (
            (self.frame.origin_mm[0] + ix as i64 * spacing) as f64 / MM_PER_M,
            (self.frame.origin_mm[1] + iy as i64 * spacing) as f64 / MM_PER_M,
        )
    }

    /// Local metres to the frame's integer millimetres (east, north).
    pub fn to_mm(x_m: f64, y_m: f64) -> [i64; 2] {
        [
            (x_m * MM_PER_M).round() as i64,
            (y_m * MM_PER_M).round() as i64,
        ]
    }

    /// Local metres -> lat/lon (degrees), equirectangular about the origin.
    pub fn to_lat_lon(&self, x_m: f64, y_m: f64) -> (f64, f64) {
        let lat0 = self.frame.origin_lat_deg.to_radians();
        let lat = self.frame.origin_lat_deg + (y_m / self.radius_m).to_degrees();
        let lon = self.frame.origin_lon_deg + (x_m / (self.radius_m * lat0.cos())).to_degrees();
        let lon = ((lon + 180.0).rem_euclid(360.0)) - 180.0;
        (lat, lon)
    }

    /// Lat/lon -> local metres (inverse of `to_lat_lon`, longitude unwrapped
    /// toward the origin).
    pub fn from_lat_lon(&self, lat_deg: f64, lon_deg: f64) -> (f64, f64) {
        let lat0 = self.frame.origin_lat_deg.to_radians();
        let dlon = ((lon_deg - self.frame.origin_lon_deg + 540.0).rem_euclid(360.0)) - 180.0;
        let x = self.radius_m * lat0.cos() * dlon.to_radians();
        let y = self.radius_m * (lat_deg - self.frame.origin_lat_deg).to_radians();
        (x, y)
    }

    /// Local cell index for local metres, if inside the grid.
    pub fn cell_of(&self, x_m: f64, y_m: f64) -> Option<usize> {
        let r = self.frame.spacing_mm as f64;
        let x0 = self.x_min_mm() as f64;
        let ix = ((x_m * MM_PER_M - x0) / r).floor();
        let iy = ((y_m * MM_PER_M - x0) / r).floor();
        if ix < 0.0
            || iy < 0.0
            || ix >= self.frame.cells_east as f64
            || iy >= self.frame.cells_north as f64
        {
            return None;
        }
        Some(iy as usize * self.frame.cells_east as usize + ix as usize)
    }

    pub fn inside(&self, x_m: f64, y_m: f64) -> bool {
        self.cell_of(x_m, y_m).is_some()
    }
}

/// Loaded and VERIFIED source artifacts.
pub struct SourceRun {
    pub terrain_digest: String,
    pub hydrology_digest: String,
    pub resource_digest: String,
    pub terrain: TerrainDomainV0,
    pub hydrology: HydrologyDomainV0,
    pub resource: ResourcePotentialDomainV0,
    pub backend: BackendId,
    pub source: SourceGrid,
    pub fields_f32: BTreeMap<String, Vec<f32>>,
    pub fields_u8: BTreeMap<String, Vec<u8>>,
    pub fields_u32: BTreeMap<String, Vec<u32>>,
}

impl SourceRun {
    pub fn load(run_root: &Path) -> Result<Self, RegionError> {
        let t_dir = run_root.join("terrain");
        let h_dir = run_root.join("hydrology");
        let r_dir = run_root.join("resource");
        let terrain_digest = verify_artifact::<TerrainDomainV0>(&t_dir)?;
        let hydrology_digest = verify_artifact::<HydrologyDomainV0>(&h_dir)?;
        let resource_digest = verify_artifact::<ResourcePotentialDomainV0>(&r_dir)?;
        let t = read_artifact_header::<TerrainDomainV0>(&t_dir)?;
        let h = read_artifact_header::<HydrologyDomainV0>(&h_dir)?;
        let r = read_artifact_header::<ResourcePotentialDomainV0>(&r_dir)?;
        if t.digest != terrain_digest || digest_header(&t.header)? != terrain_digest
            || h.digest != hydrology_digest || digest_header(&h.header)? != hydrology_digest
            || r.digest != resource_digest || digest_header(&r.header)? != resource_digest {
            return Err(RegionError::Source("source header changed during verification".into()));
        }
        // Dependency chain must be intact: hydrology cites terrain, resource cites both.
        if !h.header.semantic.inputs.contains(&terrain_digest) {
            return Err(RegionError::Source(
                "hydrology artifact does not cite the terrain digest".into(),
            ));
        }
        if !r.header.semantic.inputs.contains(&terrain_digest)
            || !r.header.semantic.inputs.contains(&hydrology_digest)
        {
            return Err(RegionError::Source(
                "resource artifact does not cite terrain and hydrology digests".into(),
            ));
        }
        t.header
            .domain
            .validate()
            .map_err(|e| RegionError::Source(e.to_string()))?;
        h.header
            .domain
            .validate()
            .map_err(|e| RegionError::Source(e.to_string()))?;
        r.header
            .domain
            .validate()
            .map_err(|e| RegionError::Source(e.to_string()))?;
        let domain = t
            .header
            .semantic
            .domains
            .first()
            .ok_or_else(|| RegionError::Source("terrain declares no domain".into()))?;
        let source = SourceGrid::from_domain(domain)?;
        if t.header.semantic.domains.len() != 1
            || h.header.semantic.domains != t.header.semantic.domains
            || r.header.semantic.domains != t.header.semantic.domains
            || h.header.domain.t_my != t.header.domain.t_my
            || r.header.domain.t_my != t.header.domain.t_my
            || h.header.semantic.backend != t.header.semantic.backend
            || r.header.semantic.backend != t.header.semantic.backend
            || t.header.semantic.backend.tier != Tier::Preview
            || t.header.semantic.backend.name != "kinematic-preview-v0"
            || t.header.semantic.backend.version != "0.2.0"
        {
            return Err(RegionError::Source(
                "source geometry/epoch/backend mismatch or unsupported producer version".into(),
            ));
        }
        let shape = vec![source.grid.nlat() as u64, source.grid.nlon() as u64];
        for manifest in [&t.header.semantic, &h.header.semantic, &r.header.semantic] {
            for f in &manifest.fields {
                if f.shape != shape || f.domain != domain.domain_id || f.order != "C" {
                    return Err(RegionError::Source(format!(
                        "field {} has incompatible shape/domain/order",
                        f.name
                    )));
                }
            }
        }

        let mut fields_f32: BTreeMap<String, Vec<f32>> = BTreeMap::new();
        let mut fields_u8 = BTreeMap::new();
        let mut fields_u32: BTreeMap<String, Vec<u32>> = BTreeMap::new();
        let mut names: std::collections::BTreeSet<String> = Default::default();
        for (dir, manifest) in [
            (&t_dir, &t.header.semantic),
            (&h_dir, &h.header.semantic),
            (&r_dir, &r.header.semantic),
        ] {
            for f in &manifest.fields {
                insert_unique_name(&mut names, &f.name, &manifest.artifact_id)?;
                let bytes = read_descriptor_bytes(dir, manifest, f)?;
                match f.dtype {
                    pf_contracts::Dtype::F32Le => {
                        fields_f32.insert(f.name.clone(), bytes.chunks_exact(4).map(|b| f32::from_le_bytes(b.try_into().expect("f32 size"))).collect());
                    }
                    pf_contracts::Dtype::U8 => {
                        fields_u8.insert(f.name.clone(), bytes);
                    }
                    pf_contracts::Dtype::U32Le => {
                        fields_u32.insert(f.name.clone(), bytes.chunks_exact(4).map(|b| u32::from_le_bytes(b.try_into().expect("u32 size"))).collect());
                    }
                }
            }
        }
        let n = source.grid.cells();
        for (name, values) in &fields_f32 {
            if values.len() != n || values.iter().any(|v| !v.is_finite()) {
                return Err(RegionError::Source(format!(
                    "nonfinite or wrong-sized {name}"
                )));
            }
        }
        let required_f32 = [
            ("elevation_m", units::M),
            ("crust_thickness_m", units::M),
            ("crust_age_my", units::MY),
            ("sediment_thickness_m", units::M),
            ("cumulative_erosion_m", units::M),
            ("uplift_rate_m_per_my", units::M_PER_MY),
            ("discharge_proxy_m3_s", units::M3_S),
            ("fill_level_m", units::M),
        ];
        for (name, unit) in required_f32.into_iter().chain(
            ResourceFamilyKindV0::ALL
                .into_iter()
                .map(|k| (k.field_name(), units::FRACTION)),
        ) {
            let desc = [&t.header.semantic, &h.header.semantic, &r.header.semantic]
                .into_iter()
                .flat_map(|m| &m.fields)
                .find(|f| f.name == name)
                .ok_or_else(|| RegionError::Source(format!("missing {name}")))?;
            if desc.dtype != pf_contracts::Dtype::F32Le || desc.units != unit {
                return Err(RegionError::Source(format!("wrong units/dtype: {name}")));
            }
        }
        for name in [
            "plate_id",
            "basin_id",
            "reach_id",
            "flow_to_element",
            "water_body_id",
        ] {
            let desc = [&t.header.semantic, &h.header.semantic]
                .into_iter()
                .flat_map(|m| &m.fields)
                .find(|f| f.name == name)
                .ok_or_else(|| RegionError::Source(format!("missing {name}")))?;
            if desc.dtype != pf_contracts::Dtype::U32Le || desc.units != units::ID {
                return Err(RegionError::Source(format!("wrong id field: {name}")));
            }
        }
        for name in ["crust_kind", "rock_family", "boundary_class"] {
            let desc = t
                .header
                .semantic
                .fields
                .iter()
                .find(|f| f.name == name)
                .ok_or_else(|| RegionError::Source(format!("missing {name}")))?;
            if desc.dtype != pf_contracts::Dtype::U8 || desc.units != units::CODE {
                return Err(RegionError::Source(format!("wrong code field: {name}")));
            }
        }
        let cf = h
            .header
            .domain
            .cell_fields
            .as_ref()
            .ok_or_else(|| RegionError::Source("hydrology missing field conventions".into()))?;
        if cf != &pf_contracts::domain::HydrologyCellFieldsV0::v0() {
            return Err(RegionError::Source(
                "unsupported hydrology field conventions".into(),
            ));
        }
        let flow = &fields_u32["flow_to_element"];
        if flow
            .iter()
            .any(|&v| v != cf.terminal_code && v as usize >= n)
        {
            return Err(RegionError::Source("receiver out of source range".into()));
        }
        let mut colour = vec![0u8; n];
        for start in 0..n {
            let mut c = start;
            let mut path = Vec::new();
            while colour[c] == 0 {
                colour[c] = 1;
                path.push(c);
                let to = flow[c];
                if to == cf.terminal_code {
                    break;
                }
                c = to as usize;
                if colour[c] == 1 {
                    return Err(RegionError::Source("cycle in source drainage".into()));
                }
            }
            for c in path {
                colour[c] = 2;
            }
        }
        let bodies: std::collections::BTreeSet<_> = h
            .header
            .domain
            .water_bodies
            .iter()
            .map(|b| b.body_id)
            .collect();
        if fields_u32["water_body_id"]
            .iter()
            .any(|v| *v != u32::MAX && !bodies.contains(v))
        {
            return Err(RegionError::Source("unknown water body membership".into()));
        }
        for body in &h.header.domain.water_bodies {
            if body.seed_element.0 >= n as u64
                || body.outlet_edges.iter().flatten().any(|e| e.0 >= n as u64)
            {
                return Err(RegionError::Source("water body edge outside source".into()));
            }
        }
        Ok(Self {
            terrain_digest,
            hydrology_digest,
            resource_digest,
            terrain: t.header.domain,
            hydrology: h.header.domain,
            resource: r.header.domain,
            backend: t.header.semantic.backend.clone(),
            source,
            fields_f32,
            fields_u8,
            fields_u32,
        })
    }

    fn f32(&self, name: &str) -> Result<&Vec<f32>, RegionError> {
        self.fields_f32
            .get(name)
            .ok_or_else(|| RegionError::Source(format!("missing f32 field {name}")))
    }

    fn u8(&self, name: &str) -> Result<&Vec<u8>, RegionError> {
        self.fields_u8
            .get(name)
            .ok_or_else(|| RegionError::Source(format!("missing u8 field {name}")))
    }

    fn u32(&self, name: &str) -> Result<&Vec<u32>, RegionError> {
        self.fields_u32
            .get(name)
            .ok_or_else(|| RegionError::Source(format!("missing u32 field {name}")))
    }
}

/// Field names across the three source artifacts share one namespace in
/// the region; a repeat is refused loudly instead of overwriting (F7).
pub fn insert_unique_name(
    names: &mut std::collections::BTreeSet<String>,
    name: &str,
    artifact_id: &str,
) -> Result<(), RegionError> {
    if !names.insert(name.to_string()) {
        return Err(RegionError::Source(format!(
            "field name {name:?} in artifact {artifact_id:?} collides with a field of another source artifact"
        )));
    }
    Ok(())
}

/// The extracted region in memory.
pub struct RegionOutputs {
    pub domain: RegionDomainV1,
    pub fields: Vec<FieldBlob>,
    pub seed_path: pf_contracts::params::SeedPath,
    pub params_digest: String,
}

/// Extract a region from a verified run.
pub fn extract(run: &SourceRun, req: &RegionRequest) -> Result<RegionOutputs, RegionError> {
    let src = &run.source;
    let frame = LocalFrame::new(req, src.grid.radius_m())?;
    let nx = frame.frame.cells_east;
    let ny = frame.frame.cells_north;
    let n = frame.cells();
    let shape = [ny as u64, nx as u64];
    let cell_size = src.cell_size_m(req.lat_deg);
    let extreme_lat = req.lat_deg.abs()
        + (req.half_extent_mm as f64 / MM_PER_M / src.grid.radius_m()).to_degrees();
    let min_spacing = src.grid.radius_m()
        * src
            .dlat_deg()
            .to_radians()
            .min(extreme_lat.to_radians().cos() * src.dlon_deg().to_radians());
    if req.spacing_mm as f64 / MM_PER_M > min_spacing {
        return Err(RegionError::Args(
            "coarsening source cells is unsupported; region retains source drainage".into(),
        ));
    }

    // Per-cell lat/lon.
    let mut lat = Vec::with_capacity(n);
    let mut lon = Vec::with_capacity(n);
    for iy in 0..ny {
        for ix in 0..nx {
            let (x, y) = frame.cell_centre_m(ix, iy);
            let (la, lo) = frame.to_lat_lon(x, y);
            lat.push(la);
            lon.push(lo);
        }
    }

    let mut fields: Vec<FieldBlob> = Vec::new();
    // Continuous fields: bilinear.
    for (name, unit) in [
        ("elevation_m", units::M),
        ("crust_thickness_m", units::M),
        ("crust_age_my", units::MY),
        ("sediment_thickness_m", units::M),
        ("cumulative_erosion_m", units::M),
        ("uplift_rate_m_per_my", units::M_PER_MY),
        ("discharge_proxy_m3_s", units::M3_S),
    ] {
        let f = run.f32(name)?;
        let data: Vec<f32> = (0..n)
            .map(|c| src.bilinear(f, lat[c], lon[c]) as f32)
            .collect();
        fields.push(f32_field(name, REGION_DOMAIN, &shape, unit, &data)?);
    }
    for kind in ResourceFamilyKindV0::ALL {
        let f = run.f32(kind.field_name())?;
        let data: Vec<f32> = (0..n)
            .map(|c| src.bilinear(f, lat[c], lon[c]).clamp(0.0, 1.0) as f32)
            .collect();
        fields.push(f32_field(
            kind.field_name(),
            REGION_DOMAIN,
            &shape,
            units::FRACTION,
            &data,
        )?);
    }
    // Categorical fields: nearest.
    let nearest: Vec<usize> = (0..n).map(|c| src.nearest(lat[c], lon[c])).collect();
    for name in ["crust_kind", "rock_family", "boundary_class"] {
        let f = run.u8(name)?;
        let data: Vec<u8> = nearest.iter().map(|&s| f[s]).collect();
        fields.push(u8_field(name, REGION_DOMAIN, &shape, units::CODE, &data)?);
    }
    for name in ["plate_id", "basin_id", "water_body_id"] {
        let f = run.u32(name)?;
        let data: Vec<u32> = nearest.iter().map(|&s| f[s]).collect();
        fields.push(u32_field(name, REGION_DOMAIN, &shape, units::ID, &data)?);
    }
    fields.push(u32_field(
        "source_element",
        REGION_DOMAIN,
        &shape,
        units::ID,
        &nearest.iter().map(|&i| i as u32).collect::<Vec<_>>(),
    )?);
    fields.push(f32_field(
        "fill_level_m",
        REGION_DOMAIN,
        &shape,
        units::M,
        &nearest
            .iter()
            .map(|&i| run.fields_f32["fill_level_m"][i])
            .collect::<Vec<_>>(),
    )?);
    let flow = run.u32("flow_to_element")?;
    let reach_of = run.u32("reach_id")?;
    let cf = run
        .hydrology
        .cell_fields
        .as_ref()
        .ok_or_else(|| RegionError::Source("missing cell fields".into()))?;
    let src_xy = |e: usize| -> (f64, f64) {
        let (i, j) = src.grid.ij(e as u64).expect("validated source element");
        frame.from_lat_lon(
            src.grid.center_lat_deg(i).expect("row"),
            src.grid.center_lon_deg(j).expect("column"),
        )
    };
    let records: BTreeMap<_, _> = run
        .hydrology
        .rivers
        .iter()
        .map(|v| (v.reach_id.0, v))
        .collect();
    let mut fragments = BTreeMap::<u32, u32>::new();
    let mut reaches = Vec::new();
    let half = req.half_extent_mm as f64 / MM_PER_M;
    let period = 2.0 * std::f64::consts::PI * src.grid.radius_m() * req.lat_deg.to_radians().cos();
    for c in 0..flow.len() {
        let to = flow[c];
        let rid = reach_of[c];
        if to == cf.terminal_code || rid == cf.no_reach_code {
            continue;
        }
        let (x0, y0) = src_xy(c);
        let (mut x1, y1) = src_xy(to as usize);
        x1 = x0 + ((x1 - x0 + period / 2.).rem_euclid(period) - period / 2.);
        if let Some((a, b)) = clip_segment([x0, y0], [x1, y1], half) {
            let record = records
                .get(&(rid as u64))
                .ok_or_else(|| RegionError::Source(format!("unknown reach {rid}")))?;
            let a = LocalFrame::to_mm(a[0], a[1]);
            let b = LocalFrame::to_mm(b[0], b[1]);
            if a == b {
                continue;
            }
            let index = fragments.entry(rid).or_default();
            reaches.push(RegionReachV1 {
                reach_id: record.reach_id,
                fragment_index: *index,
                source_edges: vec![[
                    pf_contracts::spatial::ElementId(c as u64),
                    pf_contracts::spatial::ElementId(to as u64),
                ]],
                strahler_order: record.strahler_order,
                mean_discharge_proxy_m3_s: record.mean_discharge_proxy_m3_s,
                polyline_xy_mm: vec![a, b],
                downstream: record.downstream,
            });
            *index += 1;
        }
    }
    let members: std::collections::BTreeSet<_> = nearest
        .iter()
        .map(|&c| run.fields_u32["water_body_id"][c])
        .collect();
    let water_bodies = run
        .hydrology
        .water_bodies
        .iter()
        .filter(|b| members.contains(&b.body_id))
        .cloned()
        .collect();

    let mut availability = BTreeMap::new();
    let avail = |field: &str| AvailabilityV0::Available {
        field: field.into(),
        authority_scale_m: cell_size,
    };
    availability.insert("terrain".to_string(), avail("elevation_m"));
    availability.insert("drainage".to_string(), avail("source_element"));
    availability.insert(
        "resource_potential".to_string(),
        avail(ResourceFamilyKindV0::ConvergentMagmatic.field_name()),
    );
    let unavailable = |why: &str| AvailabilityV0::Unavailable { reason: why.into() };
    availability.insert(
        "local_flow_directions".into(),
        unavailable(
            "not resolved: source_element and clipped edges retain the coarse source graph",
        ),
    );
    availability.insert(
        "coast".into(),
        unavailable(
            "sub-cell coastline, shelf morphology and sea-level history are unresolved at Preview",
        ),
    );
    availability.insert("lakes".into(), avail("water_body_id"));
    availability.insert("lithology".to_string(), unavailable("not produced at tier Preview: rock_family is a coarse process class, not lithology (pf-geology material engine, spec 55, not built)"));
    availability.insert(
        "soil".to_string(),
        unavailable("not produced at any tier yet: no weathering/soil model (spec 54.8)"),
    );
    availability.insert("flood_risk".to_string(), unavailable("not produced: discharge_proxy_m3_s is a mean rainfall-area proxy with no seasonality or return periods"));
    availability.insert(
        "seasonal_flow".to_string(),
        unavailable("not produced: uniform steady rainfall only (RainfallMode::Uniform)"),
    );
    availability.insert("channel_cross_sections".to_string(), unavailable(&format!("source flow directions are {cell_size:.0} m coarse; local channels require an explicitly authored or refined candidate (PLANET_SETTLEMENT_ART_INTEGRATION.md)")));

    let domain = RegionDomainV1 {
        frame: frame.frame.clone(),
        source: RegionSourceV0 {
            terrain_digest: run.terrain_digest.clone(),
            hydrology_digest: run.hydrology_digest.clone(),
            resource_digest: run.resource_digest.clone(),
            source_topology_id: src.topology_id.clone(),
            source_cell_size_m: cell_size,
            source_backend_name: run.backend.name.clone(),
            source_backend_tier: run.backend.tier,
            epoch_t_my: run.terrain.t_my,
            sea_level_m: run.terrain.sea_level_m,
        },
        resampling: ResamplingV0 {
            continuous: "bilinear-on-source-cell-centres (latitude clamped at pole rows, longitude periodic)".into(),
            categorical: "nearest-source-cell-centre; fill_level_m also nearest, water_body_id distinguishes true water from dry routing surface".into(),
            source_resolution_deg_lat: src.dlat_deg(),
            source_resolution_deg_lon: src.dlon_deg(),
        },
        reaches,
        water_bodies,
        availability,
        state: "proposed".into(),
    };
    domain
        .validate()
        .map_err(|e| RegionError::Source(format!("region domain: {e}")))?;
    let params_digest = digest_canonical(&serde_json::json!({
        "lat_deg": req.lat_deg,
        "lon_deg": req.lon_deg,
        "half_extent_mm": req.half_extent_mm,
        "spacing_mm": req.spacing_mm,
    }))?;
    Ok(RegionOutputs {
        domain,
        fields,
        seed_path: pf_contracts::params::SeedPath {
            seed: 0,
            derivation: "region-extract".into(),
            salt: format!("{:.6},{:.6}", req.lat_deg, req.lon_deg),
        },
        params_digest,
    })
}

/// Liang-Barsky clipping of a directed source edge. Never join separated
/// fragments through the box or drop a crossing just because its vertices lie outside.
pub fn clip_segment(a: [f64; 2], b: [f64; 2], half: f64) -> Option<([f64; 2], [f64; 2])> {
    let d = [b[0] - a[0], b[1] - a[1]];
    let mut lo = 0.0f64;
    let mut hi = 1.0f64;
    for axis in 0..2 {
        if d[axis] == 0.0 {
            if a[axis] < -half || a[axis] > half {
                return None;
            }
        } else {
            let mut enter = (-half - a[axis]) / d[axis];
            let mut exit = (half - a[axis]) / d[axis];
            if enter > exit {
                std::mem::swap(&mut enter, &mut exit);
            }
            lo = lo.max(enter);
            hi = hi.min(exit);
            if lo >= hi {
                return None;
            }
        }
    }
    Some((
        [a[0] + lo * d[0], a[1] + lo * d[1]],
        [a[0] + hi * d[0], a[1] + hi * d[1]],
    ))
}

/// Publish the region as `pf.region/v1` under `out_root` with the given
/// artifact id; inputs cite the three source digests; tier inherits
/// Preview (the minimum of its inputs, and the only tier that exists).
pub fn publish(
    out: &RegionOutputs,
    run: &SourceRun,
    out_root: &Path,
    artifact_id: &str,
    started: Instant,
) -> Result<String, RegionError> {
    let f = &out.domain.frame;
    let mut decl = DomainDeclaration {
        domain_id: REGION_DOMAIN.into(),
        kind: SpatialDomainKind::Surface2D,
        discretization: DiscretizationDescriptor {
            topology_id: format!("region.enu.{}x{}@{}mm", f.cells_north, f.cells_east, f.spacing_mm),
            backend: "region-extract-v0".into(),
            coordinate_system: CoordinateSystem::CartesianEnu,
            topology_digest: String::new(),
            metric_digest: String::new(),
            element_count: (f.cells_east as u64) * (f.cells_north as u64),
            resolution_level: u32::MAX,
            transform_conventions: format!(
                "{}; unit mm; axes {}; origin lat {:.6} lon {:.6}; {}; sample (ie,in) centre = origin_mm + (ie, in) * spacing_mm (ADR 029 s3)",
                ORDER, AXES, f.origin_lat_deg, f.origin_lon_deg, PROJECTION
            ),
        },
        datum: DatumDescriptor {
            center: [0.0, 0.0, 0.0],
            radius_m: f.datum_radius_m,
        },
    };
    decl.compute_digests()?;
    let tier = match run.backend.tier {
        Tier::Preview => Tier::Preview,
        other => other,
    };
    let backend = BackendId {
        name: REGION_BACKEND_NAME.into(),
        tier,
        version: REGION_BACKEND_VERSION.into(),
    };
    let constraints_digest = digest_canonical(&serde_json::json!({"constraints": []}))?;
    let fields: Vec<FieldBlob> = out
        .fields
        .iter()
        .map(|b| FieldBlob {
            descriptor: b.descriptor.clone(),
            bytes: b.bytes.clone(),
        })
        .collect();
    let digest = publish_typed(
        out_root,
        ArtifactKind::Region,
        artifact_id,
        &out.seed_path,
        &out.params_digest,
        &constraints_digest,
        vec![
            run.terrain_digest.clone(),
            run.hydrology_digest.clone(),
            run.resource_digest.clone(),
        ],
        backend,
        vec![decl],
        fields,
        &out.domain,
        receipt(started),
    )?;
    Ok(digest)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source(nlat: usize, nlon: usize) -> SourceGrid {
        let grid = GridSpec::new(nlat, nlon, 6.371e6).unwrap();
        SourceGrid {
            topology_id: format!("surface.latlon.{nlat}x{nlon}"),
            grid,
        }
    }

    #[test]
    fn producer_height_conversion_checks_range_before_integer_cast() {
        assert_eq!(quantize_height_mm(12.5).unwrap(), 12500);
        assert_eq!(quantize_height_mm(-12.5).unwrap(), -12500);
        for v in [1e30, -1e30, f32::INFINITY, f32::NEG_INFINITY, f32::NAN] {
            assert!(quantize_height_mm(v).is_err());
        }
    }

    #[test]
    fn bilinear_is_exact_for_linear_fields_and_matches_centres() {
        let src = source(18, 36);
        let n = src.grid.cells();
        // f = 3 i + 2 j (linear in the fractional index, away from the wrap).
        let field: Vec<f32> = (0..n)
            .map(|c| (3 * (c / 36) + 2 * (c % 36)) as f32)
            .collect();
        for (i, j) in [(4usize, 5usize), (9, 17), (12, 30)] {
            let lat = src.grid.center_lat_deg(i).unwrap();
            let lon = src.grid.center_lon_deg(j).unwrap();
            let v = src.bilinear(&field, lat, lon);
            assert!(
                (v - field[i * 36 + j] as f64).abs() < 1e-9,
                "centre sample off at ({i},{j}): {v}"
            );
            // Quarter-cell offsets: exact linear interpolation.
            let v2 = src.bilinear(
                &field,
                lat + 0.25 * src.dlat_deg(),
                lon + 0.5 * src.dlon_deg(),
            );
            let expect = 3.0 * (i as f64 + 0.25) + 2.0 * (j as f64 + 0.5);
            assert!(
                (v2 - expect).abs() < 1e-9,
                "linear exactness broken: {v2} vs {expect}"
            );
        }
        assert_eq!(
            src.nearest(
                src.grid.center_lat_deg(9).unwrap(),
                src.grid.center_lon_deg(17).unwrap()
            ),
            9 * 36 + 17
        );
    }

    #[test]
    fn duplicate_field_names_across_sources_are_refused() {
        let mut names = std::collections::BTreeSet::new();
        insert_unique_name(&mut names, "elevation_m", "terrain").unwrap();
        insert_unique_name(&mut names, "flow_to_element", "hydrology").unwrap();
        let err = insert_unique_name(&mut names, "elevation_m", "resource")
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("elevation_m") && err.contains("collides"),
            "{err}"
        );
    }

    #[test]
    fn local_frame_roundtrips_and_records_origin() {
        let req = RegionRequest {
            lat_deg: 47.25,
            lon_deg: -12.5,
            half_extent_mm: 20_000_000,
            spacing_mm: 1_000_000,
        };
        let fr = LocalFrame::new(&req, 6.371e6).unwrap();
        assert_eq!(fr.frame.cells_east, 40);
        assert_eq!(
            fr.frame.origin_mm,
            [-19_500_000, -19_500_000, 0],
            "origin is the SW cell centre in mm"
        );
        assert_eq!(fr.frame.origin_lat_deg, 47.25);
        assert_eq!(
            (
                fr.frame.unit.as_str(),
                fr.frame.axes.as_str(),
                fr.frame.order.as_str()
            ),
            ("mm", "east_north_up", "east_fastest")
        );
        assert_eq!(LocalFrame::to_mm(12.3456, -0.0004), [12_346, 0]);
        let (x, y) = (12_345.0, -6_789.0);
        let (la, lo) = fr.to_lat_lon(x, y);
        let (x2, y2) = fr.from_lat_lon(la, lo);
        assert!(
            (x - x2).abs() < 1e-6 && (y - y2).abs() < 1e-6,
            "roundtrip {x2},{y2}"
        );
        assert_eq!(fr.cell_of(-20_000.0, -20_000.0), Some(0));
        assert_eq!(fr.cell_of(19_999.0, 19_999.0), Some(40 * 40 - 1));
        assert_eq!(fr.cell_of(20_000.0, 0.0), None);
        // Origin sits at the centre of the grid.
        let (cx, cy) = fr.cell_centre_m(20, 20);
        assert_eq!((cx, cy), (500.0, 500.0));
        assert!(LocalFrame::new(
            &RegionRequest {
                lat_deg: 89.5,
                ..req.clone()
            },
            6.371e6
        )
        .is_err());
    }
}

/// Bounds are checked before the float-to-integer cast, so even extreme
/// finite input cannot saturate to i64::MIN and bypass a later abs check.
fn quantize_height_mm(height_m: f32) -> Result<i64, RegionError> {
    let mm = height_m as f64 * 1000.;
    if !mm.is_finite() || mm.abs() > 100_000_000. {
        return Err(RegionError::Args(
            "height exceeds producer preview's 100 km bound".into(),
        ));
    }
    Ok(mm.round() as i64)
}

/// Projection of a verified/published region for the producer-preview reader.
/// Cell centres become posts without spatial resampling; mm heights quantize
/// the region's f32 metre values. This is never fixture provenance or canon.
pub fn producer_tile(
    verified: &pf_contracts::region::VerifiedRegion,
    run: &SourceRun,
    id: &str,
) -> Result<serde_json::Value, RegionError> {
    let artifact = verified.artifact();
    let domain = &artifact.header.domain;
    let source = &domain.source;
    if source.terrain_digest != run.terrain_digest || source.hydrology_digest != run.hydrology_digest
        || source.resource_digest != run.resource_digest || source.source_topology_id != run.source.topology_id
        || source.source_backend_name != run.backend.name || source.source_backend_tier != run.backend.tier
        || domain.frame.datum_radius_m != run.source.grid.radius_m()
        || pf_contracts::region::declared_source_count(domain)? != run.source.grid.cells() as u64 {
        return Err(RegionError::Source("verified region does not match verified source run".into()));
    }
    // Internal region consistency is not source authority. Recompute scale at
    // this region's origin and bind temporal/datum claims to the cited source.
    for (name, claimed, actual) in [
        ("source_cell_size_m", source.source_cell_size_m, run.source.cell_size_m(domain.frame.origin_lat_deg)),
        ("epoch_t_my", source.epoch_t_my, run.terrain.t_my),
        ("sea_level_m", source.sea_level_m, run.terrain.sea_level_m),
    ] {
        if claimed != actual {
            return Err(RegionError::Source(format!("verified region {name} does not match verified source run")));
        }
    }
    // This is the actual verified upstream bound, not merely the region's declaration.
    for bytes in verified.field("source_element").expect("closed region fields").chunks_exact(4) {
        if u32::from_le_bytes(bytes.try_into().expect("u32")) as usize >= run.source.grid.cells() {
            return Err(RegionError::Source("source_element outside verified source run".into()));
        }
    }
    let region_digest = &artifact.digest;
    let f = &domain.frame;
    let n = f.cells_east;
    if n != f.cells_north || n < 2 || n as u64 * n as u64 > 262144 || (n - 1) % 25 != 0 {
        return Err(RegionError::Args(
            "producer tile requires <=262144 samples and (samples-1) divisible by chunk edge 25"
                .into(),
        ));
    }
    let height = verified.field("elevation_m").expect("closed region fields");
    let heights: Vec<i64> = height
        .chunks_exact(4)
        .map(|b| quantize_height_mm(f32::from_le_bytes(b.try_into().expect("f32 chunk"))))
        .collect::<Result<Vec<_>, _>>()?;
    let min = *heights.iter().min().expect("nonempty");
    let max = *heights.iter().max().expect("nonempty");
    if min.abs() > 100000000
        || max.abs() > 100000000
        || f.spacing_mm < 250
        || f.spacing_mm > 64000000
        || f.origin_mm.iter().any(|v| v.abs() > 1000000000)
    {
        return Err(RegionError::Args(
            "producer preview tile outside supported coordinate/height/spacing bounds".into(),
        ));
    }
    let unavailable: BTreeMap<_, _> = domain
        .availability
        .iter()
        .filter_map(|(k, v)| match v {
            AvailabilityV0::Unavailable { reason } => Some((k.clone(), reason.clone())),
            _ => None,
        })
        .collect();
    Ok(serde_json::json!({
        "format":"vz.producer-terrain-tile/v2","id":id,"worldAdmission":false,
        "label":"PRODUCER PREVIEW - tier Preview - NOT CANON - state proposed",
        "provenance":{"kind":"producer","region_schema":"pf.region/v1","backend":{"name":REGION_BACKEND_NAME,"tier":"Preview","version":REGION_BACKEND_VERSION},"source_backend":run.backend,"source_digests":{"terrain":run.terrain_digest,"hydrology":run.hydrology_digest,"resource":run.resource_digest,"region":region_digest},"authority_scale_m":domain.source.source_cell_size_m,"epoch_t_my":domain.source.epoch_t_my,"sea_level_m":domain.source.sea_level_m,"state":"proposed","height_datum":"planet_datum","sample_mapping":"posts-equal-region-cell-centres","height_quantum_mm":1,"region_half_extent_mm":f.half_extent_mm},
        "frame":{"unit":"mm","axes":"east_north_up","district_key":f.district_key,"origin_mm":f.origin_mm,"note":"Region centres become mesh posts; displayed extent excludes half a source-region spacing on each side; no district admission"},
        "grid":{"spacing_mm":f.spacing_mm,"cells_east":n-1,"cells_north":n-1,"samples_east":n,"samples_north":n,"chunk_cells":25},
        "heights_mm":{"order":"east_fastest","count":heights.len(),"min":min,"max":max,"values":heights},
        "unavailable":unavailable,"drainage":{"reaches":domain.reaches,"water_bodies":domain.water_bodies}
    }))
}
