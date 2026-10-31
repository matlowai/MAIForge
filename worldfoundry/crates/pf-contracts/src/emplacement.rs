//! Finite native body input. This contract confers no terrain or mining authority.
//! Future enrichment transfers debit a restricted parent reserve exactly; they
//! never create stock. Geometry and material must use one bound realization.
use serde::{Deserialize, Serialize};
pub const SCHEMA: &str = "vz.native-body-source/v1";
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BodyInput {
    pub schema: String,
    pub manifest_sha256: String,
}
impl BodyInput {
    pub fn validate(&self) -> Result<(), String> {
        if self.schema != SCHEMA || !digest(&self.manifest_sha256) {
            return Err("unsupported body input schema or invalid manifest SHA256".into());
        }
        Ok(())
    }
}
fn digest(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
}
/// Unknown future strings refuse decoding; never silently become a sheet.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum BodyKind { Sheet, Bed, Sill, Dyke, Pluton, Edifice, VeinSet, OreShell, Unresolved }
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum Support {
    Native { event_refs: Vec<String> },
    Conditional { assumption_refs: Vec<String> },
    Unresolved { reason: String },
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Assumption {
    pub id: String,
    pub value: f64,
    pub unit: String,
    pub interval: [f64; 2],
    pub interval_kind: String,
    pub conditions: String,
    pub reference: String,
    pub evidence_kind: String,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "vertices", rename_all = "snake_case")]
pub enum Element { Tetrahedron([u32; 4]), TriangularPrism([u32; 6]) }
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FiniteGeometry {
    pub frame_id: String,
    pub vertices_m: Vec<[f64; 3]>,
    pub elements: Vec<Element>,
    pub coverage_min_m: [f64; 3],
    pub coverage_max_m: [f64; 3],
    pub coordinate_error_m: f64,
}
impl FiniteGeometry {
    /// Structural validation only. Shared facets, nonoverlap and source closure
    /// require independent provider/producer certification before admission.
    pub fn validate(&self) -> Result<(), String> {
        if self.frame_id.is_empty() || self.vertices_m.is_empty() || self.elements.is_empty()
            || !self.coordinate_error_m.is_finite() || self.coordinate_error_m < 0.
            || (0..3).any(|i| !self.coverage_min_m[i].is_finite() || !self.coverage_max_m[i].is_finite()
                || self.coverage_min_m[i] >= self.coverage_max_m[i]) {
            return Err("invalid finite body coverage/frame".into());
        }
        for p in &self.vertices_m {
            if (0..3).any(|i| !p[i].is_finite() || p[i]<self.coverage_min_m[i] || p[i]>self.coverage_max_m[i]) {
                return Err("body vertex outside finite coverage".into());
            }
        }
        for e in &self.elements {
            let tets:Vec<[u32;4]>=match e {
                Element::Tetrahedron(t)=>vec![*t],
                Element::TriangularPrism(v)=>vec![[v[0],v[1],v[2],v[3]],[v[1],v[2],v[3],v[4]],[v[2],v[3],v[4],v[5]]],
            };
            for t in tets {
                if t.iter().any(|&i|i as usize>=self.vertices_m.len()) {return Err("invalid body vertex index".into())}
                let p=t.map(|i|self.vertices_m[i as usize]);
                let a=std::array::from_fn::<_,3,_>(|i|p[1][i]-p[0][i]);
                let b=std::array::from_fn::<_,3,_>(|i|p[2][i]-p[0][i]);
                let c=std::array::from_fn::<_,3,_>(|i|p[3][i]-p[0][i]);
                let det=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]);
                if !det.is_finite() || det<=0. {return Err("nonpositive body tetrahedron".into())}
            }
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Contact {
    pub id: String,
    pub oriented_vertices: Vec<u32>,
    pub adjacent_body_or_medium: Option<String>,
    pub lineage: Support,
    pub exposure_interval_my: Option<[f64; 2]>,
    pub exposure: Support,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BodyRecord {
    pub id: String,
    pub kind: BodyKind,
    pub parent_body: Option<String>,
    pub source_refs: Vec<String>,
    pub event_refs: Vec<String>,
    pub material_realization_id: String,
    pub assumptions: Vec<Assumption>,
    pub geometry: FiniteGeometry,
    pub contacts: Vec<Contact>,
    pub quantity_bases: Vec<QuantityBasis>,
    pub thermal: Support,
    pub consolidation: Support,
    pub terrain_join: TerrainJoin,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QuantityBasis { NativeThinShell, EuclideanEnvelope, Solid, DryMass, InterparticlePore, IntrinsicPore, Water }
/// A1 has no joined variant. A future version needs explicit terrain evidence.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TerrainJoin { Pending }
impl BodyRecord {
    /// Contract structure only; this is deliberately not a qualification gate.
    pub fn validate(&self) -> Result<(),String> {
        use std::collections::BTreeSet;
        self.geometry.validate()?;
        if self.id.is_empty() || self.parent_body.as_ref().is_some_and(|s|s.is_empty()||s==&self.id)
            || self.source_refs.is_empty() || self.source_refs.iter().any(String::is_empty)
            || self.material_realization_id.is_empty() || self.quantity_bases.is_empty()
            || self.event_refs.iter().any(String::is_empty) {
            return Err("body lacks source/material/quantity identity".into());
        }
        let mut assumptions=BTreeSet::new();
        for a in &self.assumptions {
            if a.id.is_empty() || !assumptions.insert(a.id.as_str()) || !a.value.is_finite()
                || a.interval.iter().any(|x|!x.is_finite()) || a.interval[0]>a.value || a.value>a.interval[1]
                || [&a.unit,&a.interval_kind,&a.conditions,&a.reference,&a.evidence_kind].iter().any(|s|s.is_empty()) {
                return Err("invalid conditional body assumption".into());
            }
        }
        let support=|s:&Support|->Result<(),String> {
            match s {
                Support::Native{event_refs} if !event_refs.is_empty() && event_refs.iter().all(|r|self.event_refs.contains(r))=>Ok(()),
                Support::Conditional{assumption_refs} if !assumption_refs.is_empty() && assumption_refs.iter().all(|r|assumptions.contains(r.as_str()))=>Ok(()),
                Support::Unresolved{reason} if !reason.is_empty()=>Ok(()),
                _=>Err("body support lacks bound evidence".into()),
            }
        };
        support(&self.thermal)?;support(&self.consolidation)?;
        let mut contacts=BTreeSet::new();
        for c in &self.contacts {
            if c.id.is_empty() || !contacts.insert(c.id.as_str()) || c.oriented_vertices.len()<3
                || c.oriented_vertices.iter().any(|&i|i as usize>=self.geometry.vertices_m.len())
                || c.oriented_vertices.iter().collect::<BTreeSet<_>>().len()!=c.oriented_vertices.len()
                || c.adjacent_body_or_medium.as_ref().is_some_and(String::is_empty)
                || c.exposure_interval_my.is_some_and(|v|v.iter().any(|x|!x.is_finite())||v[0]>v[1]) {
                return Err("invalid finite body contact".into());
            }
            support(&c.lineage)?;support(&c.exposure)?;
        }
        Ok(())
    }
}

// VS1 additions below. The existing native body schema above is unchanged.
use crate::planet::{GeoPoint, PlanetKey, PlanetLattice, WorldScope};
use crate::volcanic_supply::{
    self as supply, Component, DensityBasis, Knowledge, MassUg, NativeSupport, StockPhase,
    TimeInterval,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EruptionBranch {
    Intrusion,
    Lava,
    ProximalPyroclastic,
    AirborneAsh,
    VolatileRelease,
}
impl EruptionBranch {
    pub fn accepts(self, component: Component, phase: StockPhase) -> bool {
        match self {
            Self::VolatileRelease => {
                component == Component::Volatile && phase == StockPhase::ExsolvedVolatile
            }
            Self::AirborneAsh => {
                component == Component::DryRock && phase == StockPhase::AirborneSolid
            }
            _ => component == Component::DryRock && phase == StockPhase::SolidRock,
        }
    }
    fn component(self) -> Component {
        if self == Self::VolatileRelease {
            Component::Volatile
        } else {
            Component::DryRock
        }
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum TimingStatus {
    Captured {
        event_ref: String,
    },
    Conditional {
        observed_interval: TimeInterval,
        realization_id: String,
        assumption_ref: String,
    },
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EpisodeAllocation {
    pub source_id: String,
    pub component: Component,
    pub branch: EruptionBranch,
    pub amount: MassUg,
}
/// Allocation ceilings, not new stock. Unspent allocations remain in reservoir
/// accounts. Multiple source rows preserve mixtures without choosing chemistry.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EruptionEpisode {
    pub id: String,
    pub world: WorldScope,
    pub anchor: PlanetKey,
    pub reservoir_id: String,
    pub vent_id: String,
    pub vent: GeoPoint,
    pub parent_event_ref: String,
    pub interval: TimeInterval,
    pub timing: TimingStatus,
    pub allocations: Vec<EpisodeAllocation>,
    pub emplacement_law_version: String,
    pub placement_support: NativeSupport,
}
impl EruptionEpisode {
    /// Physical address and immutable parent/time, never chart/render keys.
    pub fn expected_id(&self) -> supply::Result<String> {
        let identity = supply::SourceIdentity {
            world: self.world.clone(),
            anchor: self.anchor.clone(),
            lineage: serde_json::to_string(&(
                &self.reservoir_id,
                &self.vent_id,
                &self.parent_event_ref,
            ))
            .map_err(|e| e.to_string())?,
            interval: self.interval.clone(),
        };
        Ok(format!("episode:{}", identity.id()?))
    }
    pub fn validate(&self, world: &WorldScope, lattice: &PlanetLattice) -> supply::Result<()> {
        world.validate()?;
        if &self.world != world {
            return Err("wrong world pin".into());
        }
        lattice.validate_key(&self.anchor)?;
        self.vent.validate()?;
        self.interval.validate()?;
        for s in [
            &self.reservoir_id,
            &self.vent_id,
            &self.parent_event_ref,
            &self.emplacement_law_version,
        ] {
            supply::named(s)?;
        }
        if self.id != self.expected_id()? {
            return Err("episode identity mismatch".into());
        }
        match &self.timing {
            TimingStatus::Captured { event_ref } => supply::named(event_ref)?,
            TimingStatus::Conditional {
                observed_interval,
                realization_id,
                assumption_ref,
            } => {
                observed_interval.validate()?;
                supply::named(realization_id)?;
                supply::named(assumption_ref)?;
                if !observed_interval.contains(&self.interval) {
                    return Err("conditional episode outside observed interval".into());
                }
            }
        }
        self.placement_support.validate(lattice)?;
        if self.allocations.is_empty() {
            return Err("missing episode allocations".into());
        }
        let mut seen = std::collections::BTreeSet::new();
        for a in &self.allocations {
            supply::named(&a.source_id)?;
            if a.amount == MassUg::ZERO
                || a.component != a.branch.component()
                || !seen.insert((&a.source_id, a.component, a.branch))
            {
                return Err("invalid/duplicate episode allocation".into());
            }
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum TransferPurpose {
    Eruption {
        episode_id: String,
        branch: EruptionBranch,
    },
    /// Includes ash settling/export. This spends previously emplaced stock,
    /// never the original reservoir a second time.
    Redistribution {
        evidence_ref: String,
    },
    ReservoirProcess {
        law_version: String,
    },
}
/// An edge transfers one source and component. A multi-source physical event
/// uses several exact edges, with revisions ordered in the containing receipt.
/// Geometry/refinement references confer no independent credit or stock.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EmplacementTransfer {
    pub id: String,
    pub world: WorldScope,
    pub source_id: String,
    pub from_account: String,
    pub to_account: String,
    pub component: Component,
    pub amount: MassUg,
    pub interval: TimeInterval,
    pub purpose: TransferPurpose,
    /// The source dry-rock reference basis, not a density law for volatile gas.
    pub density_basis: Knowledge<DensityBasis>,
    pub recipient_support: NativeSupport,
    pub input_revision: u64,
    pub committed_revision: u64,
    pub evidence_refs: Vec<String>,
}
impl EmplacementTransfer {
    pub fn validate(&self, world: &WorldScope, lattice: &PlanetLattice) -> supply::Result<()> {
        world.validate()?;
        if &self.world != world {
            return Err("wrong world pin".into());
        }
        for s in [
            &self.id,
            &self.source_id,
            &self.from_account,
            &self.to_account,
        ] {
            supply::named(s)?;
        }
        if self.from_account == self.to_account || self.amount == MassUg::ZERO {
            return Err("empty/self transfer".into());
        }
        self.interval.validate()?;
        self.recipient_support.validate(lattice)?;
        self.density_basis.validate_with(DensityBasis::validate)?;
        supply::refs_valid(&self.evidence_refs)?;
        match &self.purpose {
            TransferPurpose::Eruption { episode_id, .. } => supply::named(episode_id),
            TransferPurpose::Redistribution { evidence_ref } => supply::named(evidence_ref),
            TransferPurpose::ReservoirProcess { law_version } => supply::named(law_version),
        }
    }
}

/// A body is a view of existing accounts, not an additional stock owner.
/// Geometry contacts and restriction still need the later VS4 qualification.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DepositBody {
    pub id: String,
    pub world: WorldScope,
    pub account_ids: Vec<String>,
    pub episode_ids: Vec<String>,
    pub transfer_ids: Vec<String>,
    pub material: Knowledge<crate::column_envelope::MaterialRef>,
    pub formation_interval: TimeInterval,
    pub occupied_support: NativeSupport,
    pub density_basis: Knowledge<DensityBasis>,
    /// Reference only. No volume is interpreted as supply or recoverable mass.
    pub geometry_ref: Knowledge<String>,
    pub contact_refs: Knowledge<Vec<String>>,
    pub cooling_history_ref: Knowledge<String>,
    pub exposure_history_ref: Knowledge<String>,
}
impl DepositBody {
    /// Read-only reference check against a valid supply receipt. This does not
    /// certify occupied geometry, authorize mining, or resolve ADR 047.
    pub fn validate(&self, book: &supply::SupplyBook) -> supply::Result<()> {
        book.validate()?;
        if self.world != book.world {
            return Err("wrong world pin".into());
        }
        supply::named(&self.id)?;
        supply::refs_valid(&self.account_ids)?;
        supply::refs_valid(&self.episode_ids)?;
        supply::refs_valid(&self.transfer_ids)?;
        self.material.validate_with(supply::material_valid)?;
        self.density_basis.validate_with(DensityBasis::validate)?;
        self.formation_interval.validate()?;
        self.occupied_support.validate(&book.lattice)?;
        self.geometry_ref.validate_with(|s| supply::named(s))?;
        self.contact_refs.validate_with(|v| supply::refs_valid(v))?;
        self.cooling_history_ref
            .validate_with(|s| supply::named(s))?;
        self.exposure_history_ref
            .validate_with(|s| supply::named(s))?;
        let mut inherited_episodes = std::collections::BTreeSet::new();
        for id in &self.account_ids {
            let a = book
                .accounts
                .iter()
                .find(|a| &a.id == id)
                .ok_or("unknown deposit account")?;
            inherited_episodes.insert(a.origin_episode_id.require()?.clone());
            if a.component != Component::DryRock
                || a.phase != StockPhase::SolidRock
                || a.density_basis != self.density_basis
                || a.material != self.material
                || a.support != self.occupied_support
                || book.reservoirs.iter().any(|r| r.account_ids.contains(id))
            {
                return Err("deposit account basis/support mismatch".into());
            }
            if !self.transfer_ids.iter().any(|t| {
                book.transfers
                    .iter()
                    .any(|t0| &t0.id == t && &t0.to_account == id)
            }) {
                return Err("deposit account has no credited transfer".into());
            }
        }
        if inherited_episodes != self.episode_ids.iter().cloned().collect() {
            return Err("deposit episode provenance mismatch".into());
        }
        for id in &self.episode_ids {
            let e = book
                .episodes
                .iter()
                .find(|e| &e.id == id)
                .ok_or("unknown deposit episode")?;
            if !self.formation_interval.contains(&e.interval) {
                return Err("deposit formation interval mismatch".into());
            }
        }
        for id in &self.transfer_ids {
            let t = book
                .transfers
                .iter()
                .find(|t| &t.id == id)
                .ok_or("unknown deposit transfer")?;
            if !self.account_ids.contains(&t.to_account) {
                return Err("transfer does not credit deposit".into());
            }
            if let TransferPurpose::Eruption { episode_id, .. } = &t.purpose {
                if !self.episode_ids.contains(episode_id) {
                    return Err("deposit omits eruption provenance".into());
                }
            }
        }
        Ok(())
    }
}
