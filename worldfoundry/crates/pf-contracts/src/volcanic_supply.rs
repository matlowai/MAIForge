//! VS1 opt-in volcanic source contracts. No producer, geometry-to-stock
//! conversion, mantle solve, or kernel admission is implemented here.
//!
//! J0 stock is exact integer micrograms (1 ug = 1e-9 kg). Measurements and
//! conversion uncertainty belong to the producer's evidence; transferring an
//! already quantified mass introduces no rounding. Unknown is never zero.
use crate::{
    column_envelope::MaterialRef,
    planet::{PlanetKey, PlanetLattice, WorldScope},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const SCHEMA: &str = "pf.volcanic-supply/v1";
pub type Result<T> = std::result::Result<T, String>;

pub(crate) fn named(s: &str) -> Result<()> {
    if s.trim().is_empty() {
        Err("empty identity/evidence".into())
    } else {
        Ok(())
    }
}
pub(crate) fn nonnegative(v: f64) -> Result<()> {
    if !v.is_finite() || v < 0. {
        Err("invalid nonnegative quantity".into())
    } else {
        Ok(())
    }
}

/// Required tagged field: omission, null, or an empty unknown reason refuses.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum Knowledge<T> {
    Known { value: T },
    Unknown { reason: String },
}
impl<T> Knowledge<T> {
    pub fn validate_with(&self, check: impl FnOnce(&T) -> Result<()>) -> Result<()> {
        match self {
            Self::Known { value } => check(value),
            Self::Unknown { reason } => named(reason),
        }
    }
    pub fn require(&self) -> Result<&T> {
        match self {
            Self::Known { value } => Ok(value),
            Self::Unknown { .. } => Err("unknown quantity cannot fund transfer".into()),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
// Wire-compatible with pf-geology::quantities::MassUg (A1/E4). Kept here
// without a reverse dependency from contracts into geology.
#[serde(try_from = "MassWire", into = "MassWire")]
pub struct MassUg(i128);
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MassWire {
    unit: String,
    version: u8,
    decimal: String,
}
impl TryFrom<MassWire> for MassUg {
    type Error = String;
    fn try_from(w: MassWire) -> Result<Self> {
        if w.version != 1 || w.unit != "ug" {
            return Err("unsupported mass schema/unit".into());
        }
        let n: i128 = w.decimal.parse().map_err(|_| "invalid mass integer")?;
        if n.to_string() != w.decimal {
            return Err("noncanonical mass integer".into());
        }
        Self::new(n)
    }
}
impl From<MassUg> for MassWire {
    fn from(m: MassUg) -> Self {
        Self {
            unit: "ug".into(),
            version: 1,
            decimal: m.0.to_string(),
        }
    }
}
impl MassUg {
    pub const ZERO: Self = Self(0);
    pub fn new(value: i128) -> Result<Self> {
        if value < 0 {
            Err("negative stock".into())
        } else {
            Ok(Self(value))
        }
    }
    pub fn micrograms(self) -> i128 {
        self.0
    }
    pub fn checked_add(self, rhs: Self) -> Result<Self> {
        Self::new(self.0.checked_add(rhs.0).ok_or("mass overflow")?)
    }
    pub fn checked_sub(self, rhs: Self) -> Result<Self> {
        Self::new(self.0.checked_sub(rhs.0).ok_or("mass overflow")?)
    }
}

/// Forward simulation My from the pinned generation epoch, not age-before-now.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimeInterval {
    pub start_my: f64,
    pub end_my: f64,
}
impl TimeInterval {
    pub fn validate(&self) -> Result<()> {
        nonnegative(self.start_my)?;
        nonnegative(self.end_my)?;
        if self.start_my > self.end_my {
            return Err("reversed simulation interval".into());
        }
        Ok(())
    }
    pub fn contains(&self, other: &Self) -> bool {
        self.start_my <= other.start_my && other.end_my <= self.end_my
    }
}

/// Literature estimates may have no published interval. Do not invent one.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum EvidenceInterval {
    Bounds {
        lower: f64,
        upper: f64,
        interpretation: String,
    },
    NotReported {
        reason: String,
    },
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EvidenceValue {
    pub value: f64,
    pub unit: String,
    pub interval: EvidenceInterval,
    pub conditions: String,
    pub reference: String,
    pub evidence_kind: String,
}
impl EvidenceValue {
    pub fn validate(&self) -> Result<()> {
        if !self.value.is_finite() {
            return Err("nonfinite evidence value".into());
        }
        for s in [
            &self.unit,
            &self.conditions,
            &self.reference,
            &self.evidence_kind,
        ] {
            named(s)?;
        }
        match &self.interval {
            EvidenceInterval::Bounds {
                lower,
                upper,
                interpretation,
            } => {
                named(interpretation)?;
                if !lower.is_finite()
                    || !upper.is_finite()
                    || *lower > self.value
                    || self.value > *upper
                {
                    return Err("invalid evidence interval".into());
                }
            }
            EvidenceInterval::NotReported { reason } => named(reason)?,
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DensityKind {
    LegacyUniformDry,
    DenseRockEquivalentDry,
    HotMagmaDry,
    DepositBulkDry,
}
/// Descriptive basis only. VS1 offers no volume-to-mass conversion API.
/// Equality is deliberately strict: a changed basis requires reconciliation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DensityBasis {
    pub id: String,
    pub kind: DensityKind,
    pub density_kg_m3: EvidenceValue,
    pub intrinsic_pore_fraction: Knowledge<f64>,
    pub interparticle_pore_fraction: Knowledge<f64>,
}
impl DensityBasis {
    pub fn validate(&self) -> Result<()> {
        named(&self.id)?;
        self.density_kg_m3.validate()?;
        if self.density_kg_m3.unit != "kg/m3" || self.density_kg_m3.value <= 0. {
            return Err("invalid density basis".into());
        }
        if let EvidenceInterval::Bounds { lower, .. } = self.density_kg_m3.interval {
            if lower <= 0. {
                return Err("nonpositive density interval".into());
            }
        }
        for pore in [
            &self.intrinsic_pore_fraction,
            &self.interparticle_pore_fraction,
        ] {
            pore.validate_with(|p| {
                if p.is_finite() && (0.0..1.0).contains(p) {
                    Ok(())
                } else {
                    Err("invalid porosity".into())
                }
            })?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SupportCell {
    pub cell: PlanetKey,
    /// Full spherical-cell denominator, including ocean where present.
    pub whole_cell_area_m2: f64,
    pub represented_area_m2: f64,
    /// Exact share of the recorded source/support allocation, not area inferred
    /// from a rendered chart. Sum of weights equals denominator exactly.
    pub weight: u64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NativeSupport {
    pub radius_m: f64,
    pub denominator: u64,
    pub cells: Vec<SupportCell>,
    pub selection_rule: String,
    pub evidence_ref: String,
}
impl NativeSupport {
    pub fn validate(&self, lattice: &PlanetLattice) -> Result<()> {
        lattice.validate()?;
        named(&self.selection_rule)?;
        named(&self.evidence_ref)?;
        if !self.radius_m.is_finite()
            || self.radius_m <= 0.
            || self.denominator == 0
            || self.cells.is_empty()
        {
            return Err("invalid native support".into());
        }
        let mut seen = BTreeSet::new();
        let mut sum = 0u64;
        for c in &self.cells {
            lattice.validate_key(&c.cell)?;
            if !matches!(c.cell, PlanetKey::Cell { .. })
                || !seen.insert(&c.cell)
                || !c.whole_cell_area_m2.is_finite()
                || c.whole_cell_area_m2 <= 0.
                || !c.represented_area_m2.is_finite()
                || c.represented_area_m2 <= 0.
                || c.represented_area_m2 > c.whole_cell_area_m2
                || c.weight == 0
            {
                return Err("invalid/duplicate support cell".into());
            }
            sum = sum.checked_add(c.weight).ok_or("support weight overflow")?;
        }
        if sum != self.denominator {
            return Err("incomplete support weights".into());
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceIdentity {
    pub world: WorldScope,
    pub anchor: PlanetKey,
    /// Native front/slab or plume lineage, never the nearest displayed volcano.
    pub lineage: String,
    pub interval: TimeInterval,
}
impl SourceIdentity {
    pub fn id(&self) -> Result<String> {
        self.world.validate()?;
        named(&self.lineage)?;
        self.interval.validate()?;
        let scope = self.world.key("volcanic_supply", &self.anchor)?;
        let mut h = blake3::Hasher::new();
        h.update(b"pf.volcanic-source/v1\0");
        // Canonicalize signed zero for physical time identity.
        let times = [self.interval.start_my, self.interval.end_my].map(|v| {
            if v == 0. {
                0f64.to_bits()
            } else {
                v.to_bits()
            }
        });
        for b in [
            scope.as_bytes(),
            self.lineage.as_bytes(),
            &times[0].to_le_bytes(),
            &times[1].to_le_bytes(),
        ] {
            h.update(&(b.len() as u64).to_le_bytes());
            h.update(b);
        }
        Ok(format!("blake3:{}", h.finalize().to_hex()))
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum SupplyOrigin {
    /// Reduced-order boundary law; no claim of a finite mantle donor.
    BoundaryInflux {
        law_version: String,
        evidence_ref: String,
    },
    /// A later source producer must independently qualify this donor debit.
    FiniteExtraction {
        donor_account: String,
        debit_receipt_ref: String,
    },
    /// Retained/initial rock is not labeled a newly erupted supply.
    RetainedImport {
        body_account: String,
        debit_receipt_ref: String,
    },
}
impl SupplyOrigin {
    fn validate(&self) -> Result<()> {
        let (a, b) = match self {
            Self::BoundaryInflux {
                law_version,
                evidence_ref,
            } => (law_version, evidence_ref),
            Self::FiniteExtraction {
                donor_account,
                debit_receipt_ref,
            } => (donor_account, debit_receipt_ref),
            Self::RetainedImport {
                body_account,
                debit_receipt_ref,
            } => (body_account, debit_receipt_ref),
        };
        named(a)?;
        named(b)
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProductionPartition {
    pub gross_potential: Knowledge<MassUg>,
    pub admitted: Knowledge<MassUg>,
    pub capped_or_rejected: Knowledge<MassUg>,
    pub unextracted: Knowledge<MassUg>,
}
impl ProductionPartition {
    pub fn validate(&self) -> Result<()> {
        for m in [
            &self.gross_potential,
            &self.admitted,
            &self.capped_or_rejected,
            &self.unextracted,
        ] {
            m.validate_with(|_| Ok(()))?;
        }
        let mut known_sum = MassUg::ZERO;
        let mut complete = true;
        for v in [&self.admitted, &self.capped_or_rejected, &self.unextracted] {
            match v {
                Knowledge::Known { value } => known_sum = known_sum.checked_add(*value)?,
                Knowledge::Unknown { .. } => complete = false,
            }
        }
        if let Knowledge::Known { value } = self.gross_potential {
            if known_sum > value || (complete && known_sum != value) {
                return Err("production partition does not close".into());
            }
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MeltSupply {
    pub id: String,
    pub identity: SourceIdentity,
    pub model_version: String,
    pub origin: SupplyOrigin,
    pub support: NativeSupport,
    pub dry_rock: ProductionPartition,
    pub volatile: ProductionPartition,
    pub density_basis: Knowledge<DensityBasis>,
    pub material: Knowledge<MaterialRef>,
    pub provenance_refs: Vec<String>,
    /// Includes any quantization/integration receipt. No silent calibration.
    pub quantity_evidence_refs: Vec<String>,
    pub constants: Vec<EvidenceValue>,
    /// ORE-1 optional descriptive decomposition of admitted dry + volatile mass.
    /// Absence means unknown. This is not additional stock or transfer authority.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub composition: Option<crate::source_chemistry::ElementComposition>,
}
pub(crate) fn material_valid(m: &MaterialRef) -> Result<()> {
    for s in [&m.spec_id, &m.spec_version, &m.realization_id, &m.body_id] {
        named(s)?;
    }
    Ok(())
}
pub(crate) fn refs_valid(refs: &[String]) -> Result<()> {
    if refs.is_empty() {
        return Err("missing evidence references".into());
    }
    let mut seen = BTreeSet::new();
    for r in refs {
        named(r)?;
        if !seen.insert(r) {
            return Err("duplicate evidence reference".into());
        }
    }
    Ok(())
}
impl MeltSupply {
    pub fn validate(&self, world: &WorldScope, lattice: &PlanetLattice) -> Result<()> {
        world.validate()?;
        if &self.identity.world != world {
            return Err("wrong world pin".into());
        }
        lattice.validate_key(&self.identity.anchor)?;
        if self.id != self.identity.id()? {
            return Err("source identity mismatch".into());
        }
        named(&self.model_version)?;
        self.origin.validate()?;
        self.support.validate(lattice)?;
        self.dry_rock.validate()?;
        self.volatile.validate()?;
        if let Some(c) = &self.composition {
            c.validate_mass(
                self.dry_rock
                    .admitted
                    .require()?
                    .checked_add(*self.volatile.admitted.require()?)?,
            )?;
        }
        self.density_basis.validate_with(DensityBasis::validate)?;
        self.material.validate_with(material_valid)?;
        refs_valid(&self.provenance_refs)?;
        refs_valid(&self.quantity_evidence_refs)?;
        for c in &self.constants {
            c.validate()?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Component {
    DryRock,
    Volatile,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StockPhase {
    Melt,
    Crystal,
    SolidRock,
    AirborneSolid,
    DissolvedVolatile,
    ExsolvedVolatile,
}
impl StockPhase {
    fn component(self) -> Component {
        match self {
            Self::DissolvedVolatile | Self::ExsolvedVolatile => Component::Volatile,
            _ => Component::DryRock,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StockAccount {
    pub id: String,
    pub source_id: String,
    pub component: Component,
    pub phase: StockPhase,
    /// Reservoirs are explicitly unknown (not yet erupted). Emplaced lots keep
    /// one episode each, so later transport cannot erase per-eruption quantities.
    pub origin_episode_id: Knowledge<String>,
    pub material: Knowledge<MaterialRef>,
    /// Source dry-rock basis; volatile accounts carry this only as provenance,
    /// never as a gas-density or volatile volume conversion.
    pub density_basis: Knowledge<DensityBasis>,
    pub support: NativeSupport,
    pub initial: Knowledge<MassUg>,
    pub closing: Knowledge<MassUg>,
}
impl StockAccount {
    pub fn validate(&self, lattice: &PlanetLattice) -> Result<()> {
        named(&self.id)?;
        named(&self.source_id)?;
        if self.component != self.phase.component() {
            return Err("dry/volatile phase mismatch".into());
        }
        self.material.validate_with(material_valid)?;
        self.origin_episode_id.validate_with(|s| named(s))?;
        self.density_basis.validate_with(DensityBasis::validate)?;
        self.initial.validate_with(|_| Ok(()))?;
        self.closing.validate_with(|_| Ok(()))?;
        self.support.validate(lattice)
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MagmaReservoir {
    pub id: String,
    /// Source mixtures are the per-source exact account balances, not a second
    /// independently spendable aggregate. No inferred default chemistry.
    pub account_ids: Vec<String>,
    pub temperature_k: Knowledge<EvidenceValue>,
    pub volatile_model: Knowledge<String>,
    /// Closing snapshot of all referenced account masses, including volatiles.
    /// Only a projection: VS1 replay still checks bulk mass, not chemical transfers.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub composition: Option<crate::source_chemistry::ElementComposition>,
}
impl MagmaReservoir {
    fn validate(&self) -> Result<()> {
        named(&self.id)?;
        refs_valid(&self.account_ids)?;
        self.volatile_model.validate_with(|s| named(s))?;
        if let Some(c) = &self.composition {
            c.validate()?;
        }
        self.temperature_k.validate_with(|v| {
            v.validate()?;
            if v.unit != "K" || v.value <= 0. {
                return Err("invalid temperature".into());
            }
            Ok(())
        })
    }
}

/// Replayable contract receipt, not a mutable gameplay ledger. Validation is
/// pure: failure cannot partially debit an account. The kernel still owns
/// atomic persistence and admission. Source facts require producer evidence.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SupplyBook {
    pub schema: String,
    pub world: WorldScope,
    pub lattice: PlanetLattice,
    pub initial_revision: u64,
    pub committed_revision: u64,
    pub sources: Vec<MeltSupply>,
    pub reservoirs: Vec<MagmaReservoir>,
    pub accounts: Vec<StockAccount>,
    pub episodes: Vec<crate::emplacement::EruptionEpisode>,
    pub transfers: Vec<crate::emplacement::EmplacementTransfer>,
}
impl SupplyBook {
    pub fn validate(&self) -> Result<()> {
        use crate::emplacement::{TimingStatus, TransferPurpose};
        if self.schema != SCHEMA {
            return Err("unsupported supply schema".into());
        }
        self.world.validate()?;
        self.lattice.validate()?;
        let mut sources = BTreeMap::new();
        for s in &self.sources {
            s.validate(&self.world, &self.lattice)?;
            if self.sources[0].support.radius_m != s.support.radius_m {
                return Err("source/planet radius mismatch".into());
            }
            if sources.insert(s.id.as_str(), s).is_some() {
                return Err("duplicate source".into());
            }
        }
        if sources.is_empty() {
            return Err("missing sources".into());
        }
        let mut accounts = BTreeMap::new();
        let mut balances = BTreeMap::new();
        for a in &self.accounts {
            a.validate(&self.lattice)?;
            let s = sources
                .get(a.source_id.as_str())
                .ok_or("unknown account source")?;
            if a.support.radius_m != s.support.radius_m {
                return Err("account/source planet radius mismatch".into());
            }
            if a.material != s.material || a.density_basis != s.density_basis {
                return Err(
                    "incompatible density basis or material; explicit reconciliation required"
                        .into(),
                );
            }
            if accounts.insert(a.id.as_str(), a).is_some() {
                return Err("duplicate account".into());
            }
            if matches!(a.origin_episode_id, Knowledge::Known { .. })
                && a.initial
                    != (Knowledge::Known {
                        value: MassUg::ZERO,
                    })
            {
                return Err("episode stock must originate in a recorded transfer".into());
            }
            balances.insert(a.id.as_str(), a.initial.clone());
        }
        // All admitted source stock has exactly one opening owner, including
        // any reserved/unallocated quantity. Potential/capped stock cannot fund it.
        for s in &self.sources {
            for (component, admitted) in [
                (Component::DryRock, &s.dry_rock.admitted),
                (Component::Volatile, &s.volatile.admitted),
            ] {
                let aa: Vec<_> = self
                    .accounts
                    .iter()
                    .filter(|a| a.source_id == s.id && a.component == component)
                    .collect();
                let mut sum = MassUg::ZERO;
                let mut unknown = false;
                for a in aa {
                    match &a.initial {
                        Knowledge::Known { value } => sum = sum.checked_add(*value)?,
                        Knowledge::Unknown { .. } => unknown = true,
                    }
                }
                match admitted {
                    Knowledge::Known { value } if !unknown && sum == *value => (),
                    Knowledge::Unknown { .. } if unknown && sum == MassUg::ZERO => (),
                    _ => return Err("opening stock does not match admitted source branch".into()),
                }
            }
        }
        let mut reservoirs = BTreeMap::new();
        let mut reservoir_accounts = BTreeSet::new();
        for r in &self.reservoirs {
            r.validate()?;
            if reservoirs.insert(r.id.as_str(), r).is_some() {
                return Err("duplicate reservoir".into());
            }
            for id in &r.account_ids {
                let a = accounts
                    .get(id.as_str())
                    .ok_or("missing reservoir account")?;
                if !reservoir_accounts.insert(id.as_str()) {
                    return Err("account owned by two reservoirs".into());
                }
                if matches!(a.origin_episode_id, Knowledge::Known { .. }) {
                    return Err("reservoir account cannot own erupted episode stock".into());
                }
                if !matches!(
                    a.phase,
                    StockPhase::Melt
                        | StockPhase::Crystal
                        | StockPhase::DissolvedVolatile
                        | StockPhase::ExsolvedVolatile
                ) {
                    return Err("retained rock is not eruptible reservoir stock".into());
                }
                if matches!(
                    sources[a.source_id.as_str()].origin,
                    SupplyOrigin::RetainedImport { .. }
                ) {
                    return Err(
                        "retained import requires modeled remelting before reservoir admission"
                            .into(),
                    );
                }
            }
            if let Some(c) = &r.composition {
                let closing = r.account_ids.iter().try_fold(MassUg::ZERO, |sum, id| {
                    sum.checked_add(*accounts[id.as_str()].closing.require()?)
                })?;
                c.validate_mass(closing)?;
            }
        }
        let mut episodes = BTreeMap::new();
        let mut reserved = BTreeMap::<(&str, &str, Component), MassUg>::new();
        for e in &self.episodes {
            e.validate(&self.world, &self.lattice)?;
            if episodes.insert(e.id.as_str(), e).is_some() {
                return Err("duplicate episode".into());
            }
            let r = reservoirs
                .get(e.reservoir_id.as_str())
                .ok_or("missing episode reservoir")?;
            for b in &e.allocations {
                let s = sources
                    .get(b.source_id.as_str())
                    .ok_or("unknown episode source")?;
                // Stored magma can erupt after its production interval. Only
                // conditional reconstruction of unresolved subepisodes is
                // confined to the observed source interval below.
                if e.interval.start_my < s.identity.interval.start_my {
                    return Err("episode precedes source interval".into());
                }
                if let TimingStatus::Conditional {
                    observed_interval, ..
                } = &e.timing
                {
                    if !s.identity.interval.contains(observed_interval) {
                        return Err("timing evidence outside source interval".into());
                    }
                }
                if !r.account_ids.iter().any(|id| {
                    accounts[id.as_str()].source_id == b.source_id
                        && accounts[id.as_str()].component == b.component
                }) {
                    return Err("episode source is not in reservoir".into());
                }
                let available = r
                    .account_ids
                    .iter()
                    .map(|id| accounts[id.as_str()])
                    .filter(|a| a.source_id == b.source_id && a.component == b.component)
                    .try_fold(MassUg::ZERO, |sum, a| {
                        sum.checked_add(*a.initial.require()?)
                    })?;
                let reservation = reserved
                    .entry((&e.reservoir_id, &b.source_id, b.component))
                    .or_insert(MassUg::ZERO);
                *reservation = reservation.checked_add(b.amount)?;
                if *reservation > available {
                    return Err("reservoir allocations exceed admitted stock".into());
                }
            }
        }
        for a in &self.accounts {
            if let Knowledge::Known { value } = &a.origin_episode_id {
                if !episodes.contains_key(value.as_str()) {
                    return Err("unknown account episode".into());
                }
            }
        }
        let mut seen = BTreeSet::new();
        let mut revision = self.initial_revision;
        let mut allocations = BTreeMap::new();
        for t in &self.transfers {
            t.validate(&self.world, &self.lattice)?;
            if !seen.insert(t.id.as_str()) {
                return Err("duplicate transfer".into());
            }
            if t.input_revision != revision
                || t.committed_revision != revision.checked_add(1).ok_or("revision overflow")?
            {
                return Err("stale or nonsequential transfer revision".into());
            }
            let from = accounts
                .get(t.from_account.as_str())
                .ok_or("unknown debit account")?;
            let to = accounts
                .get(t.to_account.as_str())
                .ok_or("unknown credit account")?;
            if from.source_id != t.source_id
                || to.source_id != t.source_id
                || from.component != t.component
                || to.component != t.component
            {
                return Err("transfer changes source/component".into());
            }
            if from.material != to.material
                || from.density_basis != to.density_basis
                || t.density_basis != from.density_basis
            {
                return Err("incompatible density basis or material".into());
            }
            if t.interval.start_my < sources[t.source_id.as_str()].identity.interval.start_my {
                return Err("transfer precedes source".into());
            }
            if t.recipient_support != to.support {
                return Err("transfer support does not match recipient".into());
            }
            match &t.purpose {
                TransferPurpose::Eruption { episode_id, branch } => {
                    let e = episodes
                        .get(episode_id.as_str())
                        .ok_or("unknown eruption episode")?;
                    let r = reservoirs[e.reservoir_id.as_str()];
                    if !e.interval.contains(&t.interval)
                        || !r.account_ids.contains(&t.from_account)
                        || reservoir_accounts.contains(t.to_account.as_str())
                    {
                        return Err("invalid eruption debit/recipient/time".into());
                    }
                    if to.origin_episode_id.require()? != episode_id {
                        return Err("transfer changes eruption provenance".into());
                    }
                    if e.placement_support.radius_m != to.support.radius_m
                        || to.support.cells.iter().any(|c| {
                            !e.placement_support.cells.iter().any(|p| {
                                p.cell == c.cell
                                    && p.whole_cell_area_m2 == c.whole_cell_area_m2
                                    && c.represented_area_m2 <= p.represented_area_m2
                            })
                        })
                    {
                        return Err("eruption recipient outside declared placement support".into());
                    }
                    if !branch.accepts(t.component, to.phase) {
                        return Err("incompatible eruption branch".into());
                    }
                    let allocation = e
                        .allocations
                        .iter()
                        .find(|b| {
                            b.source_id == t.source_id
                                && b.component == t.component
                                && b.branch == *branch
                        })
                        .ok_or("unallocated eruption branch")?;
                    let key = (
                        episode_id.as_str(),
                        t.source_id.as_str(),
                        t.component,
                        *branch,
                    );
                    let used = allocations.entry(key).or_insert(MassUg::ZERO);
                    *used = used.checked_add(t.amount)?;
                    if *used > allocation.amount {
                        return Err("eruption allocation exceeded".into());
                    }
                }
                TransferPurpose::Redistribution { evidence_ref: _ } => {
                    if from.origin_episode_id != to.origin_episode_id {
                        return Err("redistribution changes eruption provenance".into());
                    }
                    // No disguised second eruption or return to a molten reservoir.
                    if reservoir_accounts.contains(t.from_account.as_str())
                        || reservoir_accounts.contains(t.to_account.as_str())
                    {
                        return Err("redistribution cannot tap or refill reservoir".into());
                    }
                    if from.phase != to.phase
                        && !(from.phase == StockPhase::AirborneSolid
                            && to.phase == StockPhase::SolidRock)
                    {
                        return Err("redistribution requires explicit phase-change model".into());
                    }
                }
                TransferPurpose::ReservoirProcess { law_version: _ } => {
                    if !self.reservoirs.iter().any(|r| {
                        r.account_ids.contains(&t.from_account)
                            && r.account_ids.contains(&t.to_account)
                    }) {
                        return Err("reservoir process must stay in one reservoir".into());
                    }
                }
            }
            // Check both sides before publishing even the local replay update.
            let debit = balances[t.from_account.as_str()]
                .require()?
                .checked_sub(t.amount)?;
            let credit = balances[t.to_account.as_str()]
                .require()?
                .checked_add(t.amount)?;
            balances.insert(&t.from_account, Knowledge::Known { value: debit });
            balances.insert(&t.to_account, Knowledge::Known { value: credit });
            revision = t.committed_revision;
        }
        if revision != self.committed_revision {
            return Err("closing revision mismatch".into());
        }
        for a in &self.accounts {
            if balances[a.id.as_str()] != a.closing {
                return Err("closing account does not match replay".into());
            }
        }
        Ok(())
    }
}
