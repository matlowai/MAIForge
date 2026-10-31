//! E1 source join, not a surface generator or stock admission.
//! R1-R5, integrator 2026-10-25: retain published records including thin/zero
//! layers; declare the mean anchor; ADR 047 references remain unresolved.
use crate::planet::{PlanetKey, PlanetLattice};
use serde::{Deserialize, Serialize};

pub const SCHEMA: &str = "pf.source-column-envelope/v1";
pub const MATERIAL_UNRESOLVED: &str = "unresolved: ADR 047 not accepted";
pub const MEAN_ANCHOR: &str = "declared mean-anchored construction; not native carrier elev_ref_m";
pub const ROUTING: &str = "published C7 mean routing; E2 floor routing unresolved";
pub type Result<T> = std::result::Result<T, String>;

/// An exact sum of nonnegative IEEE binary32 metre values. Little-endian limbs
/// count units of 2^-149 m. Five limbs hold even u32::MAX maximum finite f32s.
/// This is numerical encoding, not a physical length scale or tolerance.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExactThickness {
    pub binary32_quantum_limbs_le: [u64; 5],
}
impl ExactThickness {
    pub fn add_f32(&mut self, value: f32) -> Result<()> {
        if !value.is_finite() || value < 0.0 {
            return Err("invalid thickness".into());
        }
        if value == 0.0 {
            return Ok(());
        }
        let bits = value.to_bits();
        let e = (bits >> 23) & 255;
        let mantissa = (bits & 0x7fffff) | if e == 0 { 0 } else { 1 << 23 };
        let shift = e.saturating_sub(1) as usize;
        let word = shift / 64;
        let mut carry = (mantissa as u128) << (shift % 64);
        for limb in &mut self.binary32_quantum_limbs_le[word..] {
            let sum = *limb as u128 + (carry & u64::MAX as u128);
            *limb = sum as u64;
            carry = (carry >> 64) + (sum >> 64);
        }
        if carry != 0 {
            return Err("exact thickness overflow".into());
        }
        Ok(())
    }
    pub fn sum(values: impl IntoIterator<Item = f32>) -> Result<Self> {
        let mut result = Self::default();
        for v in values {
            result.add_f32(v)?;
        }
        Ok(result)
    }
    /// Display/evaluation approximation only. Exact closure uses the limbs.
    pub fn metres(self) -> f64 {
        self.binary32_quantum_limbs_le
            .iter()
            .enumerate()
            .rev()
            .map(|(i, &v)| v as f64 * 2f64.powi(64 * i as i32 - 149))
            .sum()
    }
}

/// Target wire schema matches pf-inventory::rockspec_join::MaterialRef.
/// Defining the wire type here does not adopt the proposed ADR or mint a spec.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialRef {
    pub spec_id: String,
    pub spec_version: String,
    pub realization_id: String,
    pub body_id: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PendingMaterialRef {
    pub status: String,
    pub value: Option<MaterialRef>,
}
impl PendingMaterialRef {
    /// ADR 047 is accepted. Retain the wire schema and historical unresolved
    /// records, while requiring a reason for every newly unresolved reference.
    pub fn validate(&self) -> Result<()> {
        match (&*self.status, &self.value) {
            ("resolved", Some(v)) if [&v.spec_id, &v.spec_version, &v.realization_id, &v.body_id]
                .iter().all(|s| !s.trim().is_empty()) => Ok(()),
            (s, None) if s.strip_prefix("unresolved: ").is_some_and(|r| !r.trim().is_empty()) => Ok(()),
            _ => Err("inconsistent material status/value".into()),
        }
    }
    pub fn unresolved(reason: impl Into<String>) -> Self {
        Self { status: format!("unresolved: {}", reason.into()), value: None }
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConditionalComposite {
    pub body_id: String,
    pub episode: String,
    pub realization_id: String,
    /// Native f64 contribution, not a replacement for published f32 h_k.
    pub native_retained_m: f64,
    pub deposited_m: Option<f64>,
    pub declared_initial_m: Option<f64>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Layer {
    /// The full published record's identity, never a lcont kept-layer index.
    pub record_id: String,
    pub source_csr_index: u32,
    pub column_index: u32,
    pub lithology_code: u8,
    pub rock_key: String,
    pub environment_code: u8,
    pub grade: u8,
    pub flags: u8,
    pub formation_interval_my: [f32; 2],
    pub maximum_burial_m: f32,
    pub h_retained_m: f32,
    pub zero_thickness: bool,
    /// Published JSON has no event ID. This optional native alias is supplied
    /// only by the digest-bound retained provenance export, never by matching rock.
    pub native_unit_key: Option<String>,
    pub native_episode_refs: Vec<String>,
    pub conditional_composites: Vec<ConditionalComposite>,
    pub material_ref: PendingMaterialRef,
    /// Known only for a whole layer consisting of one captured deposited episode.
    /// An initial condition is carried separately in conditional_composites.
    pub t_deposited_m: Option<f64>,
    pub depth_base_exact: ExactThickness,
    /// b_k(x) = mean_anchor - exact cumulative h + gradient dot local east/north.
    /// This f64 value is its cell-centre approximation, not an erosion surface.
    pub base_at_cell_centre_m: f64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Support {
    pub shape: [u32; 2],
    pub radius_m: f64,
    pub latitude_bounds_deg: [f64; 2],
    pub longitude_bounds_deg: [f64; 2],
    pub whole_cell_area_m2: f64,
    pub land_area_m2: f64,
    pub land_support_rule: String,
    pub sea_level_m: f64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Moments {
    pub mean_z_m: f64,
    pub trunk_c_m: f64,
    pub relief_r_m: f64,
    pub uplift_rate_m_per_my: f64,
    pub rock_uplift_m: f64,
    pub uplift_gradient_east_north: [f64; 2],
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NativeDatum {
    pub carrier: u64,
    pub elev_ref_m: Option<f64>,
    pub source: String,
    pub native_minus_mean_m: Option<f64>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpillBound {
    pub neighbour: PlanetKey,
    pub neighbour_water_kind: String,
    pub neighbour_trunk_or_stage_m: f64,
    pub side_or_corner: String,
    pub saddle_lower_bound_m: f64,
    /// Construction must provide a saddle >= lower bound. None is not a pass.
    pub constructed_saddle_m: Option<f64>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Hydraulics {
    pub routing_basis: String,
    pub published_receiver: Option<PlanetKey>,
    pub published_basin: Option<u32>,
    pub published_terminal: PlanetKey,
    pub published_final_basin: Option<u32>,
    pub floor_receiver: Option<PlanetKey>,
    pub floor_final_basin: Option<u32>,
    pub spill_bounds: Vec<SpillBound>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExposedLayer {
    pub full_record_id: String,
    pub published_h_m: f32,
    pub flags: u8,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExposureFit {
    pub source: String,
    pub publish_depth_m: f64,
    /// Actual clipped export, including the published basement top when present.
    pub layers: Vec<ExposedLayer>,
    pub thickness_exact: ExactThickness,
    pub requested_interval_m: [f64; 2],
    pub fits_translated_interval: bool,
    pub fits_mean_anchored_interval: bool,
    /// Additional non-basement cover result, not substituted for the depth gate.
    pub nonbasement_thickness_exact: ExactThickness,
    pub nonbasement_fits_translated_interval: bool,
    pub verdict_scope: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LcontDifference {
    pub dropped_positive_record_indices: Vec<u32>,
    pub dropped_positive_thickness_exact: ExactThickness,
    pub zero_record_indices: Vec<u32>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceColumnEnvelope {
    pub schema: String,
    pub planet_key: PlanetKey,
    pub stratigraphy_sha256: String,
    pub full_export_publish_depth_m: f64,
    pub support: Support,
    pub moments: Moments,
    pub construction: String,
    pub native_datum: Option<NativeDatum>,
    pub layers: Vec<Layer>,
    pub total_thickness_exact: ExactThickness,
    pub exposure_fit: ExposureFit,
    pub hydraulics: Hydraulics,
    pub lcont_difference: LcontDifference,
}
impl SourceColumnEnvelope {
    pub fn validate(&self) -> Result<()> {
        if self.schema != SCHEMA || self.construction != MEAN_ANCHOR {
            return Err("envelope schema/construction".into());
        }
        let lattice = PlanetLattice {
            kind: crate::planet::PLANET_KIND.into(),
            rows: self.support.shape[0],
            cols: self.support.shape[1],
            subdivisions: 1,
        };
        lattice.validate_key(&self.planet_key)?;
        if !matches!(self.planet_key, PlanetKey::Cell { .. }) {
            return Err("envelope needs planet cell".into());
        }
        let cell = match self.planet_key {
            PlanetKey::Cell { cell } => cell,
            _ => unreachable!(),
        };
        if self.stratigraphy_sha256.len() != 64
            || !self
                .stratigraphy_sha256
                .bytes()
                .all(|b| b.is_ascii_hexdigit())
        {
            return Err("invalid source digest".into());
        }
        let m = &self.moments;
        if [
            m.mean_z_m,
            m.trunk_c_m,
            m.relief_r_m,
            m.uplift_rate_m_per_my,
            m.rock_uplift_m,
            m.uplift_gradient_east_north[0],
            m.uplift_gradient_east_north[1],
            self.support.radius_m,
            self.support.whole_cell_area_m2,
            self.support.sea_level_m,
            self.full_export_publish_depth_m,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || m.relief_r_m < 0.
            || self.support.radius_m <= 0.
            || self.support.whole_cell_area_m2 <= 0.
            || self.support.land_area_m2 != self.support.whole_cell_area_m2
            || m.trunk_c_m != m.mean_z_m - m.relief_r_m
        {
            return Err("invalid envelope moments/support".into());
        }
        let mut sum = ExactThickness::default();
        let mut dropped = ExactThickness::default();
        let (mut thin, mut zero) = (Vec::new(), Vec::new());
        for (k, l) in self.layers.iter().enumerate() {
            if l.column_index as usize != k
                || l.zero_thickness != (l.h_retained_m == 0.)
                || l.rock_key.is_empty()
                || l.record_id
                    != format!(
                        "stratigraphy:{}:cell:{cell}:layer:{k}",
                        self.stratigraphy_sha256
                    )
                || k > 0 && l.source_csr_index != self.layers[k - 1].source_csr_index + 1
                || l.lithology_code > 11
                || l.material_ref.validate().is_err()
                || l.material_ref.value.as_ref().is_some_and(|r| r.body_id != l.record_id
                    && !l.conditional_composites.iter().any(|c| c.body_id == r.body_id))
                || l.t_deposited_m.is_some_and(|t| !t.is_finite() || t < 0.)
            {
                return Err("invalid layer identity/material".into());
            }
            sum.add_f32(l.h_retained_m)?;
            if sum != l.depth_base_exact || l.base_at_cell_centre_m != m.mean_z_m - sum.metres() {
                return Err("layer structural base does not follow declared anchor".into());
            }
            if l.h_retained_m > 0. && l.h_retained_m < 1. {
                thin.push(k as u32);
                dropped.add_f32(l.h_retained_m)?;
            }
            if l.h_retained_m == 0. {
                zero.push(k as u32);
            }
        }
        if sum != self.total_thickness_exact
            || thin != self.lcont_difference.dropped_positive_record_indices
            || zero != self.lcont_difference.zero_record_indices
            || dropped != self.lcont_difference.dropped_positive_thickness_exact
        {
            return Err("exact thickness/lcont difference does not close".into());
        }
        let e = &self.exposure_fit;
        let es = ExactThickness::sum(e.layers.iter().map(|l| l.published_h_m))?;
        let ns = ExactThickness::sum(
            e.layers
                .iter()
                .filter(|l| l.flags & crate::terrain::STRATA_FLAG_BASEMENT == 0)
                .map(|l| l.published_h_m),
        )?;
        if !e.publish_depth_m.is_finite()
            || e.publish_depth_m <= 0.
            || es != e.thickness_exact
            || ns != e.nonbasement_thickness_exact
            || e.layers.len() > self.layers.len()
            || e.layers.iter().enumerate().any(|(k, l)| {
                let Some(f) = self.layers.get(k) else {
                    return true;
                };
                l.full_record_id != f.record_id
                    || l.published_h_m > f.h_retained_m
                    || (l.flags & !crate::terrain::STRATA_FLAG_CUT_AT_DEPTH)
                        != (f.flags & !crate::terrain::STRATA_FLAG_CUT_AT_DEPTH)
                    || l.flags & crate::terrain::STRATA_FLAG_CUT_AT_DEPTH == 0
                        && l.published_h_m != f.h_retained_m
            })
            || e.requested_interval_m != [m.trunk_c_m, m.mean_z_m + m.relief_r_m]
            || e.fits_translated_interval != (es.metres() <= 2. * m.relief_r_m)
            || e.fits_mean_anchored_interval != (es.metres() <= m.relief_r_m)
            || e.nonbasement_fits_translated_interval != (ns.metres() <= 2. * m.relief_r_m)
        {
            return Err("exposure fit inconsistent".into());
        }
        if self.hydraulics.routing_basis != ROUTING
            || self.hydraulics.floor_receiver.is_some()
            || self.hydraulics.floor_final_basin.is_some()
        {
            return Err("unsupported floor routing".into());
        }
        for key in [
            &self.hydraulics.published_receiver,
            &Some(self.hydraulics.published_terminal.clone()),
        ]
        .into_iter()
        .flatten()
        {
            lattice.validate_key(key)?;
            if !matches!(key, PlanetKey::Cell { .. }) {
                return Err("hydrology must name cells".into());
            }
        }
        let mut neighbours = std::collections::BTreeSet::new();
        for b in &self.hydraulics.spill_bounds {
            if b.neighbour == self.planet_key || !neighbours.insert(b.neighbour.clone()) {
                return Err("duplicate/self spill edge".into());
            }
            lattice.validate_key(&b.neighbour)?;
            if !b.neighbour_trunk_or_stage_m.is_finite()
                || !b.saddle_lower_bound_m.is_finite()
                || b.saddle_lower_bound_m != m.trunk_c_m.max(b.neighbour_trunk_or_stage_m)
                || b.constructed_saddle_m
                    .is_some_and(|s| !s.is_finite() || s < b.saddle_lower_bound_m)
            {
                return Err("spill bound violated".into());
            }
        }
        if let Some(n) = &self.native_datum {
            if n.elev_ref_m.is_some_and(|x| !x.is_finite())
                || n.native_minus_mean_m != n.elev_ref_m.map(|x| x - m.mean_z_m)
            {
                return Err("native datum silently changed".into());
            }
        }
        Ok(())
    }
}

/// VS4 is a replacement of retained genesis rock, never a fresh magma inlet.
/// These additive records do not change the E1 v1 wire format or its fingerprint.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RetainedVolcanicBatch {
    pub source_book_sha256: String,
    pub stratigraphy_sha256: String,
    pub batch_id: String,
    pub carrier: u64,
    pub parent_sources: Vec<String>,
    pub native_episode: String,
    pub interval_my: [f64; 2],
    pub admitted: crate::volcanic_supply::MassUg,
    pub remaining: crate::volcanic_supply::MassUg,
    pub dry_density_kg_m3: f64,
    pub extrusive_efficiency: f64,
    pub efficiency_evidence: String,
    pub retention_model: String,
}
impl RetainedVolcanicBatch {
    pub fn validate(&self) -> Result<()> {
        if self.source_book_sha256.len() != 64
            || !self.source_book_sha256.bytes().all(|b| b.is_ascii_hexdigit())
            || self.stratigraphy_sha256.len()!=64 || !self.stratigraphy_sha256.bytes().all(|b|b.is_ascii_hexdigit())
            || self.batch_id.is_empty() || self.native_episode.is_empty()
            || self.efficiency_evidence.is_empty() || self.retention_model.is_empty()
            || self.parent_sources.is_empty()
            || self.parent_sources.iter().any(String::is_empty)
            || self.parent_sources.iter().collect::<std::collections::BTreeSet<_>>().len() != self.parent_sources.len()
            || self.interval_my.iter().any(|v| !v.is_finite() || *v < 0.)
            || self.interval_my[0] >= self.interval_my[1]
            || self.remaining > self.admitted
            || !self.dry_density_kg_m3.is_finite() || self.dry_density_kg_m3 <= 0.
            || !self.extrusive_efficiency.is_finite()
            || !(0.0..=1.0).contains(&self.extrusive_efficiency)
        { return Err("invalid retained volcanic batch".into()); }
        Ok(())
    }
}

/// One occupied source unit on an entire native spherical cell. Below-sea rock
/// is included in occupied volume; water is not. No layer above its roof is stock.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnitRestriction {
    pub record_id: String,
    pub owner: String,
    pub occupied_m3: f64,
    pub expected_m3: f64,
    pub integration_error_m3: f64,
    pub dry_density_kg_m3: f64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CellRestriction {
    pub cell: u32,
    pub geometry_pin: String,
    pub stratigraphy_sha256: String,
    pub whole_cell_area_m2: f64,
    pub measured_area_m2: f64,
    pub original_mean_m: f64,
    pub measured_solid_mean_m: f64,
    pub mean_error_m: f64,
    pub gross_transfer_m3: f64,
    pub previous_resolution_m3: f64,
    pub refined_resolution_m3: f64,
    /// A modeled displacement receipt is required even for a declared zero.
    pub host_deformation_model: String,
    pub host_volume_residual_m3: f64,
    pub host_integration_error_m3: f64,
    pub units: Vec<UnitRestriction>,
    /// Missing marine or external support must be a reason, not a zero volume.
    pub unknown: Vec<String>,
}
impl CellRestriction {
    /// Fail independently of a large planetary denominator. Error is bounded
    /// against this transfer, not crust volume or a loose inherited F6 tolerance.
    pub fn failures(&self) -> Vec<String> {
        let mut f = self.unknown.clone();
        let finite = [self.whole_cell_area_m2,self.measured_area_m2,self.original_mean_m,
            self.measured_solid_mean_m,self.mean_error_m,self.gross_transfer_m3,
            self.previous_resolution_m3,self.refined_resolution_m3,self.host_volume_residual_m3,self.host_integration_error_m3]
            .iter().all(|v| v.is_finite());
        if !finite { f.push("nonfinite restriction".into()); return f; }
        if self.geometry_pin.is_empty() || self.host_deformation_model.is_empty() || self.units.is_empty() {
            f.push("missing geometry, host model or unit evidence".into());
        }
        if self.whole_cell_area_m2 <= 0. || self.measured_area_m2 < 0.
            || (self.measured_area_m2-self.whole_cell_area_m2).abs() > 1e-10*self.whole_cell_area_m2 {
            f.push("incomplete whole-cell spherical support".into());
        }
        if self.mean_error_m < 0. || self.gross_transfer_m3 <= 0.
            || self.previous_resolution_m3 < 0. || self.refined_resolution_m3 < 0. {
            f.push("invalid local error denominator".into());
        }
        let numerical = 0.01*self.gross_transfer_m3;
        if (self.refined_resolution_m3-self.previous_resolution_m3).abs() > numerical
            || self.mean_error_m*self.whole_cell_area_m2 > numerical {
            f.push("geometry integration not converged to one percent of gross transfer".into());
        }
        if (self.measured_solid_mean_m-self.original_mean_m).abs() > self.mean_error_m {
            f.push("original total solid mean not preserved; datum fitting forbidden".into());
        }
        if self.host_integration_error_m3 < 0. || self.host_integration_error_m3 > numerical
            || self.host_volume_residual_m3.abs() > self.host_integration_error_m3 {
            f.push("host deformation does not conserve occupied volume within its measured integration bound".into());
        }
        let mut ids = std::collections::BTreeSet::new();
        for u in &self.units {
            if u.owner.is_empty() || u.record_id.is_empty() || !ids.insert(&u.record_id)
                || [u.occupied_m3,u.expected_m3,u.integration_error_m3,u.dry_density_kg_m3].iter().any(|v| !v.is_finite() || *v < 0.)
                || u.dry_density_kg_m3 == 0. || u.integration_error_m3 > numerical
                || (u.occupied_m3-u.expected_m3).abs() > u.integration_error_m3 {
                f.push(format!("unit volume/dry-mass restriction failed: {}",u.record_id));
            }
        }
        f
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FundedDeposit {
    pub id: String,
    pub cell: u32,
    pub record_id: String,
    pub eruption_branch: crate::emplacement::EruptionBranch,
    pub mass: crate::volcanic_supply::MassUg,
    pub occupied_m3: f64,
    /// Measured geometric error, never a mass-account balancing term.
    pub integration_error_m3: f64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiffuseWithdrawal {
    pub cell: u32,
    pub record_id: String,
    pub mass: crate::volcanic_supply::MassUg,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RetainedReplacement {
    pub id: String,
    pub event_id: String,
    pub batch_id: String,
    pub source_book_sha256: String,
    pub observed_interval_my: [f64; 2],
    pub event_timing: crate::emplacement::TimingStatus,
    pub placement_efficiency: f64,
    pub placement_efficiency_evidence: String,
    pub gross: crate::volcanic_supply::MassUg,
    pub diffuse: Vec<DiffuseWithdrawal>,
    /// Includes off-cone recipients. Nothing disappears into an unrecorded sink.
    pub deposits: Vec<FundedDeposit>,
    pub restrictions: Vec<CellRestriction>,
}

impl RetainedReplacement {
    /// Conditional reconstruction identity is independent of chart, mesh and
    /// deposit height. Changing a tile or transfer ID cannot replay one event.
    pub fn expected_event_id(&self) -> Result<String> {
        let bytes=serde_json::to_vec(&("pf.retained-reconstruction-event/v1",&self.source_book_sha256,&self.batch_id,&self.event_timing))
            .map_err(|e|e.to_string())?;
        Ok(format!("blake3:{}",blake3::hash(&bytes).to_hex()))
    }
}

/// Local proposal ledger. The engine kernel remains the only persistent mining
/// authority. All validation precedes every mutation, including replay rejection.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReplacementLedger {
    batches: std::collections::BTreeMap<String, RetainedVolcanicBatch>,
    remaining: std::collections::BTreeMap<String, crate::volcanic_supply::MassUg>,
    diffuse: std::collections::BTreeMap<(String,u32,String), crate::volcanic_supply::MassUg>,
    transfers: std::collections::BTreeSet<String>,
    events: std::collections::BTreeSet<String>,
    deposits: std::collections::BTreeMap<String, (String,String,FundedDeposit)>,
}
impl ReplacementLedger {
    /// Opening diffuse ownership must partition a finite batch exactly once.
    /// A cell copy of a carrier is never an independent batch.
    pub fn register(&mut self, batch: RetainedVolcanicBatch, owners: Vec<DiffuseWithdrawal>) -> Result<()> {
        use crate::volcanic_supply::MassUg;
        batch.validate()?;
        if self.batches.contains_key(&batch.batch_id) { return Err("duplicate finite batch".into()); }
        let mut sum = MassUg::ZERO;
        let mut distinct = std::collections::BTreeSet::new();
        for o in &owners {
            if o.record_id.is_empty() || !distinct.insert((o.cell,&o.record_id)) {
                return Err("duplicate or unnamed diffuse owner".into());
            }
            sum = sum.checked_add(o.mass)?;
        }
        if sum != batch.remaining { return Err("diffuse ownership does not partition finite retained batch".into()); }
        for o in owners { self.diffuse.insert((batch.batch_id.clone(),o.cell,o.record_id),o.mass); }
        self.remaining.insert(batch.batch_id.clone(),batch.remaining);
        self.batches.insert(batch.batch_id.clone(),batch);
        Ok(())
    }
    pub fn apply(&mut self, r: &RetainedReplacement) -> Result<()> {
        use crate::volcanic_supply::MassUg;
        let b = self.batches.get(&r.batch_id).ok_or("unknown finite batch")?;
        if r.id.is_empty() || r.event_id.is_empty() || self.transfers.contains(&r.id) || self.events.contains(&r.event_id)
            || r.source_book_sha256 != b.source_book_sha256 || r.observed_interval_my != b.interval_my
            || r.event_id != r.expected_event_id()? || r.placement_efficiency_evidence.is_empty()
            || !r.placement_efficiency.is_finite() || !(0.0..=1.0).contains(&r.placement_efficiency)
            || r.gross == MassUg::ZERO || r.gross > self.remaining[&r.batch_id] {
            return Err("replayed, unbound or overdrawn retained replacement".into());
        }
        match &r.event_timing {
            crate::emplacement::TimingStatus::Conditional{observed_interval,realization_id,assumption_ref}
                if [observed_interval.start_my,observed_interval.end_my]==b.interval_my
                    && !realization_id.is_empty() && !assumption_ref.is_empty()=>(),
            _=>return Err("retained reconstruction cannot assert a recovered eruption schedule".into()),
        }
        let gross = r.gross.micrograms() as f64 * 1e-9;
        // Efficiency describes the finite retained portion selected for this
        // reconstruction, not an inferred magma-production law.
        let efficiency_mass = b.remaining.micrograms() as f64*1e-9*r.placement_efficiency;
        let tolerance = 0.001f64.max(1e-9*gross);
        if (gross-efficiency_mass).abs() > tolerance { return Err("declared placement efficiency does not match finite withdrawal".into()); }
        let mut cells = std::collections::BTreeSet::new();
        let mut debit = MassUg::ZERO;
        let mut keys = std::collections::BTreeSet::new();
        for o in &r.diffuse {
            let key=(r.batch_id.clone(),o.cell,o.record_id.clone());
            if !keys.insert(key.clone()) || self.diffuse.get(&key).is_none_or(|m| *m < o.mass) {
                return Err("unknown, duplicate or overdrawn diffuse owner".into());
            }
            debit=debit.checked_add(o.mass)?;cells.insert(o.cell);
        }
        let mut credit = MassUg::ZERO;
        let mut ids = std::collections::BTreeSet::new();
        for d in &r.deposits {
            if d.id.is_empty() || d.record_id.is_empty()
                || matches!(d.eruption_branch,crate::emplacement::EruptionBranch::Intrusion | crate::emplacement::EruptionBranch::VolatileRelease)
                || !ids.insert(&d.id) || self.deposits.contains_key(&d.id) || d.mass == MassUg::ZERO
                || !d.occupied_m3.is_finite() || d.occupied_m3 <= 0.
                || !d.integration_error_m3.is_finite() || d.integration_error_m3 < 0.
                || d.integration_error_m3 > 0.01*d.occupied_m3 {
                return Err("invalid or duplicate positive deposit".into());
            }
            let kg=d.mass.micrograms() as f64*1e-9;
            let geometry_kg=d.occupied_m3*b.dry_density_kg_m3;
            let geometry_error_kg=d.integration_error_m3*b.dry_density_kg_m3;
            if !geometry_kg.is_finite() || !geometry_error_kg.is_finite()
                || (kg-geometry_kg).abs() > tolerance+geometry_error_kg {
                return Err("deposit occupied volume and dry mass disagree".into());
            }
            credit=credit.checked_add(d.mass)?;cells.insert(d.cell);
        }
        if debit != r.gross || credit != r.gross { return Err("per-transfer exact dry mass does not close; missing sink or source".into()); }
        let mut measured=std::collections::BTreeSet::new();
        for c in &r.restrictions {
            if c.stratigraphy_sha256 != b.stratigraphy_sha256 {return Err("replacement geometry/source stratigraphy pin mismatch".into());}
            if !measured.insert(c.cell) || !c.failures().is_empty() { return Err(format!("whole-cell restriction failed for {}: {:?}",c.cell,c.failures())); }
            let withdrawals=r.diffuse.iter().filter(|o|o.cell==c.cell);
            let deposits=r.deposits.iter().filter(|d|d.cell==c.cell);
            let mut local_debit=MassUg::ZERO;
            let mut local_credit=MassUg::ZERO;
            for o in withdrawals {
                local_debit=local_debit.checked_add(o.mass)?;
                let key=(r.batch_id.clone(),o.cell,o.record_id.clone());
                let remaining=self.diffuse[&key].checked_sub(o.mass)?;
                let u=c.units.iter().find(|u|u.record_id==o.record_id && u.owner==r.batch_id)
                    .ok_or("withdrawn diffuse owner missing from restricted geometry")?;
                if u.dry_density_kg_m3 != b.dry_density_kg_m3
                    || (u.expected_m3*b.dry_density_kg_m3-remaining.micrograms() as f64*1e-9).abs()>tolerance {
                    return Err("diffuse restricted volume does not match closing ownership".into());
                }
            }
            for d in deposits {
                local_credit=local_credit.checked_add(d.mass)?;
                if r.diffuse.iter().any(|o|o.cell==d.cell && o.record_id==d.record_id){return Err("replacement needs distinct refined ownership".into());}
                let u=c.units.iter().find(|u|u.record_id==d.record_id && u.owner==r.batch_id)
                    .ok_or("positive deposit missing from restricted geometry")?;
                if u.dry_density_kg_m3 != b.dry_density_kg_m3
                    || (u.expected_m3*b.dry_density_kg_m3-d.mass.micrograms() as f64*1e-9).abs()>tolerance
                    || (u.occupied_m3-d.occupied_m3).abs()>u.integration_error_m3+d.integration_error_m3 {
                    return Err("deposit restriction and transfer disagree".into());
                }
            }
            // An in-cell transfer is counted once, not as both a debit and a
            // credit that doubles the permitted numerical error.
            let local_gross=local_debit.max(local_credit);
            if (c.gross_transfer_m3*b.dry_density_kg_m3-local_gross.micrograms() as f64*1e-9).abs()>tolerance {
                return Err("cell numerical denominator is not its actual gross transfer".into());
            }
        }
        if measured != cells { return Err("missing or unrelated affected cell restriction".into()); }
        // All fallible validation finished. Checked arithmetic above proves these.
        for o in &r.diffuse {
            let key=(r.batch_id.clone(),o.cell,o.record_id.clone());
            let old=self.diffuse[&key];self.diffuse.insert(key,old.checked_sub(o.mass)?);
        }
        self.remaining.insert(r.batch_id.clone(),self.remaining[&r.batch_id].checked_sub(r.gross)?);
        for d in &r.deposits {self.deposits.insert(d.id.clone(),(r.batch_id.clone(),r.event_id.clone(),d.clone()));}
        self.transfers.insert(r.id.clone());self.events.insert(r.event_id.clone());
        Ok(())
    }
    pub fn funded(&self, id: &str) -> Option<&(String,String,FundedDeposit)> { self.deposits.get(id) }
}
