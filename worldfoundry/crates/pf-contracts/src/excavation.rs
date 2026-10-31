//! Opt-in excavation contracts. Mass is decimal-encoded i128 micrograms.
//! Geometry and accounts are separate: exact ledger closure is not a volume error bound.
use crate::planet::{PlanetKey, WorldScope};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
pub type Result<T> = std::result::Result<T, String>;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub struct Mass(pub i128);
impl TryFrom<String> for Mass {
    type Error = String;
    fn try_from(s: String) -> Result<Self> {
        let v = s.parse::<i128>().map_err(|_| "invalid_mass")?;
        if v < 0 || v.to_string() != s {
            return Err("negative_or_noncanonical_mass".into());
        }
        Ok(Self(v))
    }
}
impl From<Mass> for String {
    fn from(v: Mass) -> Self {
        v.0.to_string()
    }
}
pub type Accounts = BTreeMap<String, Mass>;
/// Checked floor(a*b/d), without an overflowing wide product. No float stock.
pub fn proportional(a: i128, b: i128, d: i128) -> Result<i128> {
    if a < 0 || b < 0 || d <= 0 || a > d {
        return Err("invalid_ratio".into());
    }
    let (a, b, d) = (a as u128, b as u128, d as u128);
    let (mut q, mut r) = (0u128, 0u128);
    for bit in (0..128).rev() {
        q = q.checked_mul(2).ok_or("mass_overflow")?;
        if r >= d - r {
            r -= d - r;
            q += 1;
        } else {
            r *= 2;
        }
        if (b >> bit) & 1 != 0 {
            if r >= d - a {
                r -= d - a;
                q += 1;
            } else {
                r += a;
            }
        }
    }
    i128::try_from(q).map_err(|_| "mass_overflow".into())
}
pub fn debit(before: &Accounts, removed: &Accounts) -> Result<Accounts> {
    if before.is_empty() || before.keys().ne(removed.keys()) {
        return Err("mineral_set_mismatch".into());
    }
    before
        .iter()
        .map(|(m, v)| {
            let d = removed[m].0;
            if v.0 < 0 || d < 0 {
                return Err("negative_stock".into());
            }
            let after =
                v.0.checked_sub(d)
                    .filter(|v| *v >= 0)
                    .ok_or("negative_stock")?;
            Ok((m.clone(), Mass(after)))
        })
        .collect()
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoxMm {
    pub min: [i64; 3],
    pub max: [i64; 3],
}
impl BoxMm {
    pub fn validate(&self) -> Result<()> {
        if (0..3).any(|i| {
            self.min[i] >= self.max[i]
                || self.min[i].unsigned_abs() > 1_000_000_000_000
                || self.max[i].unsigned_abs() > 1_000_000_000_000
        }) {
            return Err("invalid_box".into());
        }
        Ok(())
    }
    pub fn contains(&self, p: [i64; 3]) -> bool {
        (0..3).all(|i| p[i] >= self.min[i] && p[i] < self.max[i])
    }
    pub fn overlaps(&self, b: &Self) -> bool {
        (0..3).all(|i| self.min[i] < b.max[i] && b.min[i] < self.max[i])
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Operation {
    Dig,
    /// Versioned energy-free extraction to the initiating actor's inventory.
    /// The digest binds the complete wire input, including aiming metadata.
    MineV1 { request_digest: String },
    Blast {
        energy_j: u64,
    },
    /// Whole-asset emplacement into its original excavated void. No compaction model.
    Backfill {
        asset_id: String,
    },
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExcavationContract {
    pub request_id: String,
    /// Canonical actor UUID; the native kernel checks that it is alive.
    pub actor_id: String,
    pub world: WorldScope,
    pub address: PlanetKey,
    pub body_id: String,
    pub frame_id: String,
    pub bounds: BoxMm,
    pub operation: Operation,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    pub body_id: String,
    pub unit_id: String,
    pub formation_event: String,
    pub material_realization: String,
    pub snapshot: String,
    pub packet: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MiningAdmission {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub blast: Option<BlastAdmission>,
    pub world: WorldScope,
    pub address: PlanetKey,
    pub frame_id: String,
    pub source: Provenance,
    pub initial: Accounts,
    /// Floor of E4's allocated lower physical volume, mm3, decimal i128.
    pub allocated_volume_mm3: Mass,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialAsset {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elements: Option<PhaseElements>,
    pub asset_id: String,
    pub class: String,
    pub source: Provenance,
    pub minerals: Accounts,
    /// Original cut bounds. MineV1 assets are actor-owned inventory, not an
    /// in-world deposit at these bounds; old Dig assets retain their meaning.
    pub rest_bounds: BoxMm,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Excavation {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub volume: Option<crate::volume::OrderedEdit>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub blast: Option<serde_json::Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chemistry: Option<ElementTransfer>,
    pub request: ExcavationContract,
    pub physical_volume_mm3: [Mass; 2],
    pub removed: Accounts,
    pub assets: Vec<MaterialAsset>,
    pub fines: Accounts,
    pub halo_mm: u64,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MiningState {
    pub admission: MiningAdmission,
    pub remaining: Accounts,
    pub edits: Vec<Excavation>,
}
impl MiningState {
    pub fn validate(&self) -> Result<()> {
        self.admission.world.validate()?;
        // Reuse only a previously successful fold over this exact immutable
        // prefix. Every existing check below still runs for each new edit.
        let prefix = crate::mining_history::prefix(self);
        // A complete successful prefix already includes the exact owner fold.
        // Do not deep-copy its entire owner inventory merely to discard it.
        // A new suffix or any changed owner field still runs the same checks.
        let owners=if prefix.as_ref().is_some_and(|p|p.state.edits.len()==self.edits.len()) {
            None
        } else {Some(self.owner_inventory()?)};
        let mut remaining = prefix.as_ref().map_or_else(
            || self.admission.initial.clone(),
            |p| p.state.remaining.clone(),
        );
        let mut elements = if let Some(p) = &prefix {
            p.elements.clone()
        } else {
            self.admission.initial_elements()?
        };
        let ordered = self.edits.first().is_some_and(|e| e.volume.is_some());
        let mut consumed = prefix
            .as_ref()
            .map_or_else(std::collections::BTreeSet::new, |p| p.consumed.clone());
        let start = prefix.as_ref().map_or(0, |p| p.state.edits.len());
        for (i, e) in self.edits.iter().enumerate().skip(start) {
            e.request.bounds.validate()?;
            crate::player_inventory::validate_edit(e)?;
            if e.volume.is_some() != ordered {
                return Err("mixed_legacy_ordered_history".into());
            }
            if let Some(v) = &e.volume {
                v.domain.validate()?;
                if v.before.identity.world != self.admission.world
                    || v.before.identity.address != self.admission.address
                    || v.before.identity.frame_id != self.admission.frame_id
                    || v.before.identity.domain_id != self.admission.source.body_id
                    || (!crate::world_admission::dynamic(&self.admission) && self
                        .edits
                        .first()
                        .and_then(|x| x.volume.as_ref())
                        .is_some_and(|x| x.domain != v.domain))
                {
                    return Err("ordered_domain_mismatch".into());
                }
                if v.before.state_digest != crate::mining_history::state_digest_parts(
                    &self.admission, &remaining, &self.edits[..i],
                )? {
                    return Err("ordered_state_chain_mismatch".into());
                }
            }
            if e.request.world != self.admission.world
                || e.request.address != self.admission.address
                || e.request.body_id != self.admission.source.body_id
                || e.request.frame_id != self.admission.frame_id
                || (e.physical_volume_mm3[0].0 <= 0 && !(e.physical_volume_mm3 == [Mass(0);2]
                    && e.removed.values().all(|n|n.0==0)
                    && e.chemistry.as_ref().and_then(|t|t.owners.as_ref()).is_some_and(|t|!t.transfers.is_empty())))
                || e.physical_volume_mm3[0] > e.physical_volume_mm3[1]
                || e.halo_mm != 0
            {
                return Err("edit_admission_mismatch".into());
            }
            if self.edits[..i].iter().any(|old| {
                old.request.request_id == e.request.request_id
                    || (!ordered && old.request.bounds.overlaps(&e.request.bounds))
            }) {
                return Err("duplicate_or_overlapping_edit".into());
            }
            if let Operation::Backfill { asset_id } = &e.request.operation {
                if !ordered || !consumed.insert(asset_id.clone()) {
                    return Err("asset_already_consumed_or_legacy_fill".into());
                }
                let origin = self.edits[..i]
                    .iter()
                    .find(|x| x.assets.iter().any(|a| &a.asset_id == asset_id))
                    .ok_or("unknown_backfill_asset")?;
                let asset = origin
                    .assets
                    .iter()
                    .find(|a| &a.asset_id == asset_id)
                    .unwrap();
                let t = e.chemistry.as_ref().ok_or("backfill_elements_required")?;
                if !matches!(origin.request.operation, Operation::Dig)
                    || origin.assets.len() != 1
                    || origin.fines.values().any(|v| v.0 != 0)
                    || e.request.bounds != origin.request.bounds
                    || e.removed != asset.minerals
                    || Some(&t.removed) != asset.elements.as_ref()
                    || !e.assets.is_empty()
                    || e.blast.is_some()
                    || e.fines.keys().ne(e.removed.keys())
                    || e.fines.values().any(|v| v.0 != 0)
                    || t.fines.values().flat_map(|x| x.values()).any(|v| v.0 != 0)
                    || !t.class_allocations.is_empty()
                    || e.physical_volume_mm3 != origin.physical_volume_mm3
                {
                    return Err("invalid_backfill_transfer".into());
                }
                check_elements(&e.removed, &t.removed)?;
                check_elements(&e.fines, &t.fines)?;
                add(&mut remaining, &e.removed)?;
                let stock = elements.as_mut().ok_or("backfill_elements_required")?;
                for (m, v) in &t.removed {
                    add(stock.get_mut(m).ok_or("element_phase_set")?, v)?;
                }
                continue;
            }
            if e.fines.values().any(|v| v.0 < 0)
                || (e.assets.is_empty() && e.fines.values().all(|v| v.0 == 0)
                    && !e.chemistry.as_ref().and_then(|t|t.owners.as_ref()).is_some_and(|t|!t.transfers.is_empty()))
            {
                return Err("invalid_material_credits".into());
            }
            let mut credits = e.fines.clone();
            for (j, a) in e.assets.iter().enumerate() {
                if a.source != self.admission.source
                    || a.rest_bounds != e.request.bounds
                    || a.minerals.keys().ne(credits.keys())
                    || a.asset_id != format!("{}/asset/{j}", e.request.request_id)
                {
                    return Err("asset_provenance_mismatch".into());
                }
                for (m, v) in &a.minerals {
                    if v.0 < 0 {
                        return Err("negative_stock".into());
                    }
                    credits.get_mut(m).unwrap().0 =
                        credits[m].0.checked_add(v.0).ok_or("mass_overflow")?;
                }
            }
            if credits != e.removed {
                return Err("unbalanced_minerals".into());
            }
            if matches!(e.request.operation, Operation::Blast { .. }) != e.blast.is_some()
                || (e.blast.is_some() && self.admission.blast.is_none())
            {
                return Err("blast_receipt_required".into());
            }
            match (&mut elements, &e.chemistry) {
                (Some(stock), Some(t)) => {
                    check_elements(&e.removed, &t.removed)?;
                    check_elements(&e.fines, &t.fines)?;
                    let mut credit = t.fines.clone();
                    for a in &e.assets {
                        let ae = a.elements.as_ref().ok_or("asset_elements_required")?;
                        check_elements(&a.minerals, ae)?;
                        for (m, v) in ae {
                            add(credit.get_mut(m).ok_or("element_phase_set")?, v)?;
                        }
                    }
                    if credit != t.removed {
                        return Err("unbalanced_elements".into());
                    }
                    if t.class_allocations.keys().ne(t.removed.keys()) {
                        return Err("missing_class_receipts".into());
                    }
                    for (m, records) in &t.class_allocations {
                        if records.keys().ne(t.removed[m].keys()) {
                            return Err("missing_class_element_receipt".into());
                        }
                        for (element, receipt) in records {
                            let weights = receipt
                                .terms
                                .iter()
                                .map(|r| (r.key.clone(), r.weight))
                                .collect();
                            if receipt.total.micrograms() != t.removed[m][element].0
                                || crate::source_chemistry::allocate(receipt.total, weights)?
                                    != *receipt
                            {
                                return Err("invalid_class_allocation_receipt".into());
                            }
                            for term in &receipt.terms {
                                let amount = if term.key == "fines" {
                                    t.fines[m][element]
                                } else {
                                    let a = e
                                        .assets
                                        .iter()
                                        .find(|a| a.asset_id == term.key)
                                        .ok_or("unknown_credit_asset")?;
                                    a.elements.as_ref().ok_or("asset_elements_required")?[m]
                                        [element]
                                };
                                if amount.0 != term.amount.micrograms() {
                                    return Err("class_receipt_credit_mismatch".into());
                                }
                            }
                        }
                    }
                    for (m, v) in &t.removed {
                        stock.insert(m.clone(), debit(&stock[m], v)?);
                    }
                }
                (None, None) if e.assets.iter().all(|a| a.elements.is_none()) => (),
                _ => return Err("element_admission_mismatch".into()),
            }
            remaining = debit(&remaining, &e.removed)?;
        }
        if remaining != self.remaining || remaining.values().any(|x| x.0 < 0) {
            return Err("remaining_stock_mismatch".into());
        }
        // An exact full-prefix hit already owns this successful fold. The
        // remaining-stock comparison above still runs; avoid comparing all
        // retained bytes a second time merely to discover the same cache entry.
        if prefix.is_none() || start < self.edits.len() {
            crate::mining_history::remember(self, elements, consumed, owners.expect("new prefix owner fold"));
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_wide_ratios_and_negative_stock() {
        assert_eq!(
            proportional(i128::MAX - 1, i128::MAX, i128::MAX).unwrap(),
            i128::MAX - 1
        );
        assert_eq!(proportional(1, i128::MAX, 3).unwrap(), i128::MAX / 3);
        assert!(proportional(-1, 1, 1).is_err());
        let a = BTreeMap::from([("m".into(), Mass(4))]);
        assert!(debit(&a, &BTreeMap::from([("m".into(), Mass(5))])).is_err());
        assert!(debit(&a, &BTreeMap::from([("m".into(), Mass(-1))])).is_err());
        assert!(serde_json::from_str::<Mass>("\"-1\"").is_err());
    }
}

/// Each phase's elements describe that same mineral mass, not extra stock.
pub type PhaseElements = BTreeMap<String, Accounts>;
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BlastAdmission {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resource: Option<crate::resource_admission::ResourceAdmission>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub world_support: Option<crate::world_admission::WorldSupport>,
    /// Additive capability, absent on all legacy admissions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub owner_transfer_version: Option<String>,
    pub model: serde_json::Value,
    pub resolved_inputs: serde_json::Value,
    pub chemistry: BTreeMap<String, crate::source_chemistry::StoichiometryReceipt>,
    pub evidence_kind: String,
    pub source: Provenance,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElementTransfer {
    /// Separate owner debits. `removed` below remains the bedrock debit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub owners: Option<OwnerTransfers>,
    pub removed: PhaseElements,
    pub fines: PhaseElements,
    pub class_allocations:
        BTreeMap<String, BTreeMap<String, crate::source_chemistry::AllocationReceipt<String>>>,
}
pub fn add(to: &mut Accounts, value: &Accounts) -> Result<()> {
    if to.keys().ne(value.keys()) {
        return Err("element_set_mismatch".into());
    }
    for (k, v) in value {
        if v.0 < 0 {
            return Err("negative_stock".into());
        }
        to.get_mut(k).unwrap().0 = to[k].0.checked_add(v.0).ok_or("mass_overflow")?;
    }
    Ok(())
}
pub fn check_elements(phases: &Accounts, elements: &PhaseElements) -> Result<()> {
    if phases.keys().ne(elements.keys()) {
        return Err("element_phase_set".into());
    }
    for (m, stock) in elements {
        if stock.is_empty() {
            return Err("unknown_element_stock".into());
        }
        let sum = stock.values().try_fold(0i128, |s, v| {
            if v.0 < 0 {
                return Err("negative_stock");
            }
            s.checked_add(v.0).ok_or("mass_overflow")
        })?;
        if sum != phases[m].0 {
            return Err("element_phase_mass_mismatch".into());
        }
    }
    Ok(())
}
impl MiningAdmission {
    pub fn initial_elements(&self) -> Result<Option<PhaseElements>> {
        let Some(a) = &self.blast else {
            return Ok(None);
        };
        if a.source != self.source
            || a.evidence_kind.is_empty()
            || !a.model.is_object()
            || !a.resolved_inputs.is_object()
            || (a.resource.is_none() && a.chemistry.keys().ne(self.initial.keys()))
        {
            return Err("invalid_blast_admission".into());
        }
        if let Some(resource) = &a.resource {
            if !a.chemistry.is_empty() { return Err("resource_and_legacy_chemistry_conflict".into()); }
            return resource.validate_for(self).map(Some);
        }
        let mut result = PhaseElements::new();
        for (m, r) in &a.chemistry {
            r.composition.validate()?;
            if r.request.species != *m || r.composition.phase_mass.micrograms() != self.initial[m].0
            {
                return Err("chemistry_source_mass_mismatch".into());
            }
            result.insert(
                m.clone(),
                r.composition
                    .elements
                    .iter()
                    .map(|e| (format!("{:?}", e.element), Mass(e.mass.micrograms())))
                    .collect(),
            );
        }
        check_elements(&self.initial, &result)?;
        Ok(Some(result))
    }
}
impl MiningState {
    pub fn remaining_elements(&self) -> Result<Option<PhaseElements>> {
        self.validate()?;
        let Some(mut result) = self.admission.initial_elements()? else {
            return Ok(None);
        };
        for e in &self.edits {
            for (m, v) in &e
                .chemistry
                .as_ref()
                .ok_or("element_transfer_required")?
                .removed
            {
                if matches!(e.request.operation, Operation::Backfill { .. }) {
                    add(result.get_mut(m).ok_or("element_phase_set")?, v)?;
                } else {
                    result.insert(m.clone(), debit(&result[m], v)?);
                }
            }
        }
        Ok(Some(result))
    }
}
/// Deterministic prefix apportionment. Each bin gets a difference of cumulative
/// floors. Sum is EXACT, no float stock and no rounding microgram disappears.
pub fn partition(total: i128, weights: &Accounts) -> Result<Accounts> {
    let d = weights.values().try_fold(0i128, |s, v| {
        if v.0 < 0 {
            Err("negative_weight")
        } else {
            s.checked_add(v.0).ok_or("mass_overflow")
        }
    })?;
    if total < 0 || d <= 0 {
        return Err("invalid_partition".into());
    }
    let (mut cumulative, mut previous) = (0i128, 0i128);
    weights
        .iter()
        .map(|(k, w)| {
            cumulative += w.0;
            let next = proportional(cumulative, total, d)?;
            let amount = next - previous;
            previous = next;
            Ok((k.clone(), Mass(amount)))
        })
        .collect()
}

#[cfg(test)]
mod blast_accounting_tests {
    use super::*;
    #[test]
    fn prefix_partition_closes_at_i128_limit_and_cannot_create_stock() {
        let weights = Accounts::from([("a".into(), Mass(i128::MAX - 7)), ("b".into(), Mass(7))]);
        for n in [0, 1, 2, 23_039_943_219_000, i128::MAX - 1, i128::MAX] {
            let parts = partition(n, &weights).unwrap();
            assert_eq!(parts.values().map(|v| v.0).sum::<i128>(), n);
            assert!(parts.iter().all(|(k, v)| *v <= weights[k]));
            debit(&weights, &parts).unwrap();
        }
        assert!(partition(-1, &weights).is_err());
        assert!(partition(1, &Accounts::from([("a".into(), Mass(-1))])).is_err());
    }
    #[test]
    fn transported_element_vectors_close_across_three_fragment_classes() {
        let initial = Accounts::from([
            ("O".into(), Mass(1000000000000000001)),
            ("Si".into(), Mass(999999999999999999)),
        ]);
        let selected = partition(23039943219, &initial).unwrap();
        let w = Accounts::from([
            ("coarse".into(), Mass(41)),
            ("fine".into(), Mass(7)),
            ("large".into(), Mass(52)),
        ]);
        let mut recovered = initial
            .keys()
            .map(|k| (k.clone(), Mass(0)))
            .collect::<Accounts>();
        for (element, mass) in &selected {
            for (_, part) in partition(mass.0, &w).unwrap() {
                recovered.get_mut(element).unwrap().0 += part.0;
            }
        }
        assert_eq!(selected, recovered);
        let after = debit(&initial, &selected).unwrap();
        add(&mut recovered, &after).unwrap();
        assert_eq!(initial, recovered);
    }
}

/// Versioned re-excavation capability. Geometry is independently recomputed by
/// pf-volume at candidate review and again before canonization.
pub const OWNER_TRANSFERS_V1: &str = "vz.rubble-owner-transfers/1";
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OwnerTransfer {
    pub owner_id: String,
    pub source: Provenance,
    pub minerals: Accounts,
    pub elements: PhaseElements,
    pub assets: Vec<MaterialAsset>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OwnerTransfers {
    pub schema: String,
    pub transfers: Vec<OwnerTransfer>,
    pub removed_bounds: BoxMm,
    /// Exact serialized replacement projections keyed by original request ID.
    /// Immutable historical products remain available; current views use these.
    pub residual_deposits: BTreeMap<String, String>,
}
impl MiningState {
    /// Fold only actual asset inventories. Never resample the parent chemistry.
    /// Used by validation and proposal code; it does not itself authorize writes.
    pub fn owner_inventory(&self) -> Result<BTreeMap<String, MaterialAsset>> {
        let prefix=crate::mining_history::prefix(self);
        let mut owners=prefix.as_ref().map_or_else(BTreeMap::new,|p|p.owners.clone());
        let start=prefix.as_ref().map_or(0,|p|p.state.edits.len());
        if self.admission.blast.as_ref().and_then(|a| a.owner_transfer_version.as_deref())
            .is_some_and(|v| v != OWNER_TRANSFERS_V1) {
            return Err("unknown_owner_transfer_version".into());
        }
        for edit in &self.edits[start..] {
            if let Some(batch) = edit.chemistry.as_ref().and_then(|t| t.owners.as_ref()) {
                if batch.schema != OWNER_TRANSFERS_V1
                    || self.admission.blast.as_ref().and_then(|a| a.owner_transfer_version.as_deref()) != Some(OWNER_TRANSFERS_V1)
                    || !matches!(edit.request.operation, Operation::Blast { .. } | Operation::MineV1 { .. }) {
                    return Err("owner_transfer_not_admitted".into());
                }
                batch.removed_bounds.validate()?;
                let mut seen = std::collections::BTreeSet::new();
                // No transfer may debit an output born in the same candidate.
                let mut outputs = vec![];
                for t in &batch.transfers {
                    if !seen.insert(&t.owner_id) { return Err("duplicate_owner_debit".into()); }
                    let owner = owners.get_mut(&t.owner_id).ok_or("unknown_rubble_owner")?;
                    if owner.source != t.source || !owner.class.starts_with("rubble/v1/") {
                        return Err("rubble_owner_provenance_mismatch".into());
                    }
                    check_elements(&t.minerals, &t.elements)?;
                    let total = t.minerals.values().try_fold(0i128, |s,n| s.checked_add(n.0).ok_or("mass_overflow"))?;
                    if total <= 0 { return Err("empty_owner_debit".into()); }
                    let total_before = owner.minerals.values().try_fold(0i128, |s,n| s.checked_add(n.0).ok_or("mass_overflow"))?;
                    if total > total_before || partition(total, &owner.minerals)? != t.minerals {
                        return Err("owner_composition_mismatch".into());
                    }
                    let stock = owner.elements.as_mut().ok_or("rubble_elements_required")?;
                    for (m, n) in &t.minerals {
                        // Every element comes from this owner's current stock.
                        let own = stock.get(m).ok_or("element_phase_set")?;
                        let expected = if n.0 == 0 { own.keys().map(|k|(k.clone(),Mass(0))).collect() }
                            else { partition(n.0, own)? };
                        if expected != t.elements[m] {
                            return Err("owner_element_composition_mismatch".into());
                        }
                        stock.insert(m.clone(), debit(&stock[m], &t.elements[m])?);
                    }
                    owner.minerals = debit(&owner.minerals, &t.minerals)?;
                    let mut credit: Accounts = t.minerals.keys().map(|k|(k.clone(),Mass(0))).collect();
                    let mut ce: PhaseElements = t.elements.iter().map(|(m,es)|
                        (m.clone(), es.keys().map(|e|(e.clone(),Mass(0))).collect())).collect();
                    for asset in &t.assets {
                        let credit_id = if matches!(edit.request.operation, Operation::MineV1 { .. }) {
                            crate::player_inventory::rubble_lot_id(&edit.request.request_id, &t.owner_id)
                        } else { rubble_credit_id(&edit.request.request_id, &t.source, &asset.class)? };
                        let class_ok = if matches!(edit.request.operation, Operation::MineV1 { .. }) {
                            asset.class == crate::player_inventory::owner_class(&edit.request.actor_id)
                        } else { asset.class.starts_with("rubble/v1/") };
                        if asset.source != t.source || asset.rest_bounds != edit.request.bounds
                            || asset.asset_id != credit_id || !class_ok {
                            return Err("owner_credit_provenance_mismatch".into());
                        }
                        let es = asset.elements.as_ref().ok_or("rubble_elements_required")?;
                        check_elements(&asset.minerals, es)?;
                        add(&mut credit, &asset.minerals)?;
                        for (m, v) in es { add(ce.get_mut(m).ok_or("element_phase_set")?, v)?; }
                        outputs.push(asset.clone());
                    }
                    if credit != t.minerals || ce != t.elements { return Err("unbalanced_owner_transfer".into()); }
                }
                for asset in outputs {
                    if seen.contains(&asset.asset_id) { return Err("owner_cycle".into()); }
                    if let Some(existing) = owners.get_mut(&asset.asset_id) {
                        if existing.source != asset.source || existing.class != asset.class || existing.rest_bounds != asset.rest_bounds {
                            return Err("incompatible_owner_merge".into());
                        }
                        add(&mut existing.minerals, &asset.minerals)?;
                        let es = existing.elements.as_mut().ok_or("rubble_elements_required")?;
                        for (m, v) in asset.elements.as_ref().ok_or("rubble_elements_required")? {
                            add(es.get_mut(m).ok_or("element_phase_set")?, v)?;
                        }
                    } else { owners.insert(asset.asset_id.clone(), asset); }
                }
            }
            for asset in &edit.assets {
                if owners.insert(asset.asset_id.clone(), asset.clone()).is_some() { return Err("duplicate_asset_owner".into()); }
            }
            if let Operation::Backfill { asset_id } = &edit.request.operation {
                owners.remove(asset_id).ok_or("unknown_backfill_asset")?;
            }
        }
        Ok(owners)
    }
}

/// Merge only equal source provenance and equal material/size classes. The
/// per-input transfer contributions stay in the event, preserving the DAG.
pub fn rubble_credit_id(request: &str, source: &Provenance, class: &str) -> Result<String> {
    let bytes = serde_json::to_vec(&(source,class)).map_err(|e|e.to_string())?;
    Ok(format!("{request}/asset/group/{}",blake3::hash(&bytes)))
}
