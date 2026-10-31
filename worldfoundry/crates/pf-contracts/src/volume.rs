//! Volume ABI v1. Bricks are projections, never stock accounts or event authority.
//! Rest coordinates are fixed once in a pinned source frame. Chart/window/camera
//! IDs are not part of a brick's identity. W1 alone decides chart ownership.
use crate::{
    excavation::{check_elements, Accounts, PhaseElements},
    planet::*,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub type Result<T> = std::result::Result<T, String>;
pub const TICKS_PER_METRE: i64 = 32_000;
pub const EDGE: i64 = 16;
pub const MIXED: u16 = u16::MAX;
pub const SCHEMA: &str = "vz.rock-volume/1";

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TickBox {
    pub min: [i64; 3],
    pub max: [i64; 3],
}
impl TickBox {
    pub fn validate(&self) -> Result<()> {
        if (0..3).any(|i| {
            self.min[i] >= self.max[i]
                || self.min[i].unsigned_abs() > 1u64 << 48
                || self.max[i].unsigned_abs() > 1u64 << 48
        }) {
            return Err("invalid_volume_bounds".into());
        }
        Ok(())
    }
    pub fn contains(&self, p: [i64; 3]) -> bool {
        (0..3).all(|i| p[i] >= self.min[i] && p[i] < self.max[i])
    }
    /// Half-open ownership. A boundary point belongs to exactly one child.
    pub fn clip(&self, b: &Self) -> Option<Self> {
        let r = Self {
            min: std::array::from_fn(|i| self.min[i].max(b.min[i])),
            max: std::array::from_fn(|i| self.max[i].min(b.max[i])),
        };
        (0..3).all(|i| r.min[i] < r.max[i]).then_some(r)
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VolumeIdentity {
    pub schema: String,
    pub world: WorldScope,
    pub address: PlanetKey,
    pub frame_id: String,
    pub frame_digest: String,
    pub domain_id: String,
}
impl VolumeIdentity {
    pub fn validate(&self, lattice: &PlanetLattice) -> Result<()> {
        self.world.validate()?;
        lattice.validate_key(&self.address)?;
        validate_digest(&self.frame_digest)?;
        if self.schema != SCHEMA || self.frame_id.is_empty() || self.domain_id.is_empty() {
            return Err("invalid_volume_identity".into());
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReadGuard {
    pub identity: VolumeIdentity,
    pub revision: u64,
    pub state_digest: String,
    pub read_set_digest: String,
}
impl ReadGuard {
    pub fn validate(&self, lattice: &PlanetLattice) -> Result<()> {
        self.identity.validate(lattice)?;
        validate_digest(&self.state_digest)?;
        validate_digest(&self.read_set_digest)
    }
    pub fn require(&self, expected: &Self) -> Result<()> {
        if self != expected {
            return Err("stale_or_foreign_volume_snapshot".into());
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BrickCoord {
    pub xyz: [i64; 3],
    pub level: u8,
}
impl BrickCoord {
    /// h = 2^level / 32 metres. Version 1 allows 3.125 cm through 2 m.
    pub fn step(level: u8) -> Result<i64> {
        if level > 6 {
            return Err("unsupported_volume_level".into());
        }
        Ok(1000i64 << level)
    }
    pub fn at(p: [i64; 3], level: u8) -> Result<Self> {
        let w = Self::step(level)? * EDGE;
        Ok(Self {
            xyz: p.map(|v| v.div_euclid(w)),
            level,
        })
    }
    pub fn bounds(&self) -> Result<TickBox> {
        let w = Self::step(self.level)? * EDGE;
        let mut min = [0; 3];
        let mut max = [0; 3];
        for i in 0..3 {
            min[i] = self.xyz[i].checked_mul(w).ok_or("brick_overflow")?;
            max[i] = min[i].checked_add(w).ok_or("brick_overflow")?;
        }
        let b = TickBox { min, max };
        b.validate()?;
        Ok(b)
    }
    /// Local halo vertex indices -1..=17; shared vertices use identical integers.
    pub fn vertex(&self, local: [i64; 3]) -> Result<[i64; 3]> {
        if local.iter().any(|v| !(-1..=17).contains(v)) {
            return Err("invalid_halo_index".into());
        }
        let b = self.bounds()?;
        let h = Self::step(self.level)?;
        Ok(std::array::from_fn(|i| b.min[i] + local[i] * h))
    }
}
/// Source status survives geometry projection. Construction excess is not stock.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SourceStatus {
    Qualified,
    Conditional,
    ConstructionExcess,
    BelowRecord,
    Unresolved,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceRef {
    pub body_id: String,
    pub record_id: String,
    pub material_id: String,
    pub stock_owner: String,
    pub source_snapshot: String,
    pub status: SourceStatus,
    pub qualification: String,
}
impl SourceRef {
    pub fn validate(&self) -> Result<()> {
        if [
            &self.body_id,
            &self.record_id,
            &self.material_id,
            &self.stock_owner,
            &self.source_snapshot,
            &self.qualification,
        ]
        .iter()
        .any(|s| s.is_empty())
        {
            return Err("incomplete_source_reference".into());
        }
        Ok(())
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ThinStatus {
    Resolved,
    AnalyticAttachment,
    RefinementRequired,
    Unresolved,
}
/// No per-voxel mass. The exact transfer ledger remains outside this type.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BrickEnvelope {
    pub guard: ReadGuard,
    pub coord: BrickCoord,
    pub palette: Vec<SourceRef>,
    pub thin_status: ThinStatus,
    pub boundary_digest: String,
}
/// A projection chart is a view only. Never compare a legacy frame.chart to W1.
pub fn chart_owns(atlas: &Atlas, view_chart: u32, geo: GeoPoint) -> Result<bool> {
    atlas.centre(view_chart)?;
    Ok(atlas.owner(geo)? == view_chart)
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ChangeKind {
    Removal,
    Deposition,
    Weakening,
    Split,
    Motion,
    Exposure,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AccountTransfer {
    pub source: SourceRef,
    pub debit_account: String,
    pub credit_account: String,
    pub minerals: Accounts,
    pub phase_elements: PhaseElements,
}
/// JOINT/SEEPAGE producer output, not an admitted event. Geometry is a referenced
/// 3D payload, never a height subtraction. Only VOL-EVENTS may activate it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChangeProposal {
    pub proposal_id: String,
    pub guard: ReadGuard,
    pub kind: ChangeKind,
    pub domain: TickBox,
    pub geometry_digest: String,
    pub geometry_bytes: u64,
    pub interval_start_us: i64,
    pub interval_end_us: i64,
    pub process_receipt_digest: String,
    pub transfers: Vec<AccountTransfer>,
    pub dependency_ranges: Vec<TickBox>,
    pub unresolved: Vec<String>,
}
impl ChangeProposal {
    /// Structural checks only. This neither debits accounts nor grants admission.
    pub fn validate(&self, expected: &ReadGuard, lattice: &PlanetLattice) -> Result<()> {
        self.guard.validate(lattice)?;
        self.guard.require(expected)?;
        self.domain.validate()?;
        for s in [
            &self.proposal_id,
            &self.geometry_digest,
            &self.process_receipt_digest,
        ] {
            validate_digest(s)?;
        }
        if self.geometry_bytes == 0
            || self.interval_start_us > self.interval_end_us
            || self.dependency_ranges.is_empty()
        {
            return Err("incomplete_process_proposal".into());
        }
        for b in &self.dependency_ranges {
            b.validate()?;
        }
        if matches!(self.kind, ChangeKind::Removal | ChangeKind::Deposition)
            && self.transfers.is_empty()
        {
            return Err("geometry_without_transfers".into());
        }
        let mut seen = BTreeSet::new();
        for t in &self.transfers {
            t.source.validate()?;
            if !matches!(
                t.source.status,
                SourceStatus::Qualified | SourceStatus::Conditional
            ) || t.debit_account.is_empty()
                || t.credit_account.is_empty()
                || t.debit_account == t.credit_account
                || t.minerals.is_empty()
                || t.minerals.values().any(|v| v.0 < 0)
                || !t.minerals.values().any(|v| v.0 > 0)
                || !seen.insert((&t.source.stock_owner, &t.debit_account, &t.credit_account))
            {
                return Err("invalid_source_transfer".into());
            }
            check_elements(&t.minerals, &t.phase_elements)?;
        }
        Ok(())
    }
}

/// Ordered native command. The entire prior state is a conservative range guard:
/// even an insertion into a previously empty dependency range invalidates it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OrderedCommand {
    pub before: ReadGuard,
    pub domain: crate::excavation::BoxMm,
    pub request: crate::excavation::ExcavationContract,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OrderedEdit {
    pub before: ReadGuard,
    pub domain: crate::excavation::BoxMm,
    /// Disjoint finite pieces in domain-local metres, added to or subtracted
    /// from the void. Recomputed by native review, never trusted as authority.
    pub pieces: String,
    /// Round-trippable f64 chart-volume diagnostic in m3, not exact stock.
    pub chart_volume_m3: String,
}
