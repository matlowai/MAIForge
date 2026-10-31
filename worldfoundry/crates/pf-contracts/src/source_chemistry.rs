//! ORE-1: descriptive chemistry of the SAME phase mass, never a second stock.
//! No source initialization, reaction law, assay, or kernel admission is implied.
//! Integer closure is exact; geological uncertainty and atomic-weight uncertainty
//! remain in evidence. All endmember fractions are MASS fractions, not molar An/Fo.
use crate::volcanic_supply::{named, EvidenceValue, MassUg, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const SCHEMA: &str = "pf.source-chemistry/v1";
pub const TABLE_SCHEMA: &str = "pf.ore-chemistry-table/v1";

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Element {
    H,
    C,
    N,
    O,
    Na,
    Mg,
    Al,
    Si,
    S,
    Cl,
    K,
    Ca,
    Ti,
    Mn,
    Fe,
    Cu,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElementMass {
    pub element: Element,
    pub mass: MassUg,
}

/// Complete decomposition within a declared realization. An absent vector means
/// unknown, not zero. Omitted elements in a PRESENT vector are zero only within
/// that declared model; this is not a claim about natural trace concentrations.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElementComposition {
    pub schema: String,
    pub realization_id: String,
    pub evidence_ref: String,
    pub phase_mass: MassUg,
    pub elements: Vec<ElementMass>,
}
impl ElementComposition {
    pub fn validate(&self) -> Result<()> {
        if self.schema != SCHEMA {
            return Err("unsupported composition schema".into());
        }
        named(&self.realization_id)?;
        named(&self.evidence_ref)?;
        if self.elements.is_empty() {
            return Err("empty composition is unknown".into());
        }
        let mut seen = BTreeSet::new();
        let mut total = MassUg::ZERO;
        for row in &self.elements {
            if !seen.insert(row.element) {
                return Err("duplicate element".into());
            }
            total = total.checked_add(row.mass)?;
        }
        if total != self.phase_mass {
            return Err("element mass differs from phase mass".into());
        }
        Ok(())
    }
    pub fn validate_mass(&self, expected: MassUg) -> Result<()> {
        self.validate()?;
        if self.phase_mass != expected {
            return Err("composition snapshot mass mismatch".into());
        }
        Ok(())
    }
}

/// Oxide reporting is a separate basis. Equivalent oxide mass may include
/// hypothetical oxygen and MUST NOT be added to an elemental/phase inventory.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OxideBasis {
    ActualOxide,
    AnalyticalEquivalent,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OxideMass {
    pub formula: String,
    pub mass: MassUg,
    pub basis: OxideBasis,
    pub evidence_ref: String,
}
impl OxideMass {
    pub fn validate(&self) -> Result<()> {
        if !matches!(
            self.formula.as_str(),
            "SiO2"
                | "Al2O3"
                | "FeO"
                | "Fe2O3"
                | "MgO"
                | "CaO"
                | "Na2O"
                | "K2O"
                | "TiO2"
                | "MnO"
                | "CuO"
                | "SO3"
                | "CO2"
                | "H2O"
        ) {
            return Err("unsupported oxide formula".into());
        }
        named(&self.evidence_ref)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtomicWeight {
    pub element: Element,
    /// Exact encoding of the chosen conventional atomic weight in 1e-6 u.
    pub micro_u: u64,
    pub evidence: EvidenceValue,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElementWeight {
    pub element: Element,
    pub weight: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Stoichiometry {
    /// Integer atom counts may be a common multiple of a fractional formula.
    Formula { atoms: Vec<ElementWeight> },
    /// Noncrystalline organic matter has no universal mineral formula.
    DeclaredMass { weights: Vec<ElementWeight> },
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Endmember {
    pub id: String,
    pub label: String,
    pub stoichiometry: Stoichiometry,
    /// Formula: value 1, unit formula unit; declared mass: value 1, unit kg/kg.
    /// Conditions describe substitutions excluded by the idealization.
    pub evidence: EvidenceValue,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpeciesGroup {
    pub id: String,
    pub allowed_endmembers: Vec<String>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChemistryTable {
    pub schema: String,
    pub version: String,
    pub atomic_weights: Vec<AtomicWeight>,
    pub endmembers: Vec<Endmember>,
    pub groups: Vec<SpeciesGroup>,
}
impl ChemistryTable {
    pub fn from_json(text: &str) -> Result<Self> {
        let table: Self = serde_json::from_str(text).map_err(|e| e.to_string())?;
        table.validate()?;
        Ok(table)
    }
    pub fn validate(&self) -> Result<()> {
        if self.schema != TABLE_SCHEMA {
            return Err("unsupported chemistry table".into());
        }
        named(&self.version)?;
        let mut elements = BTreeSet::new();
        for a in &self.atomic_weights {
            a.evidence.validate()?;
            if a.micro_u == 0
                || !elements.insert(a.element)
                || a.evidence.unit != "u"
                || (a.evidence.value * 1e6 - a.micro_u as f64).abs() > 0.01
            {
                return Err("invalid atomic weight encoding".into());
            }
        }
        let mut ids = BTreeSet::new();
        if self.endmembers.is_empty() {
            return Err("empty chemistry table".into());
        }
        for e in &self.endmembers {
            named(&e.id)?;
            named(&e.label)?;
            e.evidence.validate()?;
            if !ids.insert(e.id.as_str()) {
                return Err("duplicate endmember".into());
            }
            self.element_weights(e)?;
        }
        for g in &self.groups {
            named(&g.id)?;
            if !ids.insert(&g.id) || g.allowed_endmembers.is_empty() {
                return Err("duplicate/empty species group".into());
            }
            let mut seen = BTreeSet::new();
            for e in &g.allowed_endmembers {
                if !seen.insert(e) || !self.endmembers.iter().any(|m| &m.id == e) {
                    return Err("invalid group endmember".into());
                }
            }
        }
        Ok(())
    }
    /// Content pin, not just a caller-supplied version string. A changed formula,
    /// reference or atomic-weight encoding invalidates previous receipts.
    pub fn digest(&self) -> Result<String> {
        self.validate()?;
        let bytes = serde_json::to_vec(self).map_err(|e| e.to_string())?;
        Ok(format!("blake3:{}", blake3::hash(&bytes).to_hex()))
    }
    fn element_weights(&self, endmember: &Endmember) -> Result<BTreeMap<Element, u64>> {
        let (rows, atoms) = match &endmember.stoichiometry {
            Stoichiometry::Formula { atoms } => (atoms, true),
            Stoichiometry::DeclaredMass { weights } => (weights, false),
        };
        let mut out = BTreeMap::new();
        for r in rows {
            let a = self
                .atomic_weights
                .iter()
                .find(|a| a.element == r.element)
                .ok_or("unknown element atomic weight")?;
            let weight = r
                .weight
                .checked_mul(if atoms { a.micro_u } else { 1 })
                .ok_or("formula weight overflow")?;
            if weight == 0 || out.insert(r.element, weight).is_some() {
                return Err("zero/duplicate formula constituent".into());
            }
        }
        sum_weights(out.values().copied())?;
        Ok(out)
    }
    pub fn decompose(&self, request: &PhaseChemistryRequest) -> Result<StoichiometryReceipt> {
        let table_digest = self.digest()?;
        named(&request.species)?;
        named(&request.realization_id)?;
        named(&request.evidence_ref)?;
        let group = self.groups.iter().find(|g| g.id == request.species);
        let weights = match &request.split {
            None => {
                if group.is_some() {
                    return Err("species group requires a declared endmember split".into());
                }
                if !self.endmembers.iter().any(|e| e.id == request.species) {
                    return Err("unknown species".into());
                }
                BTreeMap::from([(request.species.clone(), 1)])
            }
            Some(split) => {
                named(&split.label)?;
                named(&split.evidence_ref)?;
                if split.basis != "mass_fraction" {
                    return Err("only mass fractions supported".into());
                }
                let mut weights = BTreeMap::new();
                for e in &split.endmembers {
                    let allowed = match group {
                        Some(g) => g.allowed_endmembers.contains(&e.endmember),
                        None => e.endmember == request.species,
                    };
                    if !allowed
                        || e.numerator == 0
                        || weights.insert(e.endmember.clone(), e.numerator).is_some()
                    {
                        return Err("invalid/duplicate split endmember".into());
                    }
                }
                if sum_weights(weights.values().copied())? != split.denominator {
                    return Err("endmember fractions do not close".into());
                }
                weights
            }
        };
        let endmember_allocation = allocate(request.phase_mass, weights)?;
        let mut components = Vec::new();
        let mut totals = BTreeMap::<Element, MassUg>::new();
        for part in &endmember_allocation.terms {
            let e = self
                .endmembers
                .iter()
                .find(|e| e.id == part.key)
                .ok_or("unknown endmember")?;
            let allocation = allocate(part.amount, self.element_weights(e)?)?;
            for term in &allocation.terms {
                let entry = totals.entry(term.key).or_insert(MassUg::ZERO);
                *entry = entry.checked_add(term.amount)?;
            }
            components.push(EndmemberReceipt {
                endmember: e.id.clone(),
                allocation,
            });
        }
        let composition = ElementComposition {
            schema: SCHEMA.into(),
            realization_id: request.realization_id.clone(),
            evidence_ref: request.evidence_ref.clone(),
            phase_mass: request.phase_mass,
            elements: totals
                .into_iter()
                .map(|(element, mass)| ElementMass { element, mass })
                .collect(),
        };
        composition.validate()?;
        Ok(StoichiometryReceipt {
            table_digest,
            request: request.clone(),
            endmember_allocation,
            components,
            composition,
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EndmemberFraction {
    pub endmember: String,
    pub numerator: u64,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EndmemberSplit {
    pub basis: String,
    pub label: String,
    pub evidence_ref: String,
    pub denominator: u64,
    pub endmembers: Vec<EndmemberFraction>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhaseChemistryRequest {
    pub species: String,
    pub phase_mass: MassUg,
    pub realization_id: String,
    pub evidence_ref: String,
    pub split: Option<EndmemberSplit>,
}

/// floor + remainder_numerator/denominator is the exact pre-rounding ug value.
/// Largest remainder gets +1 ug; ties use the key's Ord order. Each term's
/// quantization error is strictly less than 1 ug. No unnamed residue exists.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AllocationTerm<K> {
    pub key: K,
    pub weight: u64,
    pub floor: MassUg,
    pub remainder_numerator: u64,
    pub rounded_up: bool,
    pub amount: MassUg,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AllocationReceipt<K> {
    pub total: MassUg,
    pub denominator: u64,
    pub terms: Vec<AllocationTerm<K>>,
}
fn sum_weights(mut weights: impl Iterator<Item = u64>) -> Result<u64> {
    let sum = weights.try_fold(0u64, |sum, w| sum.checked_add(w).ok_or("weight overflow"))?;
    if sum == 0 {
        return Err("empty/zero weights".into());
    }
    Ok(sum)
}
pub fn allocate<K: Ord + Clone>(
    total: MassUg,
    weights: BTreeMap<K, u64>,
) -> Result<AllocationReceipt<K>> {
    let denominator = sum_weights(weights.values().copied())?;
    let n = total.micrograms();
    let d = i128::from(denominator);
    let mut terms = Vec::new();
    let mut floors = MassUg::ZERO;
    for (key, weight) in weights {
        if weight == 0 {
            return Err("zero weight".into());
        }
        // n*w can overflow for planet-scale i128 stocks. Divide n first; the
        // remainder product is at most (u64::MAX)^2 and fits u128 exactly.
        let small = (n % d) as u128 * u128::from(weight);
        let floor = MassUg::new(
            (n / d)
                .checked_mul(i128::from(weight))
                .and_then(|v| v.checked_add((small / u128::from(denominator)) as i128))
                .ok_or("allocation overflow")?,
        )?;
        floors = floors.checked_add(floor)?;
        terms.push(AllocationTerm {
            key,
            weight,
            floor,
            remainder_numerator: (small % u128::from(denominator)) as u64,
            rounded_up: false,
            amount: floor,
        });
    }
    let residual = usize::try_from(total.checked_sub(floors)?.micrograms())
        .map_err(|_| "residual overflow")?;
    if residual >= terms.len() {
        return Err("invalid rounding residual".into());
    }
    let mut order: Vec<_> = (0..terms.len()).collect();
    order.sort_by(|&a, &b| {
        terms[b]
            .remainder_numerator
            .cmp(&terms[a].remainder_numerator)
            .then_with(|| terms[a].key.cmp(&terms[b].key))
    });
    for &i in order.iter().take(residual) {
        terms[i].rounded_up = true;
        terms[i].amount = terms[i].floor.checked_add(MassUg::new(1)?)?;
    }
    Ok(AllocationReceipt {
        total,
        denominator,
        terms,
    })
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EndmemberReceipt {
    pub endmember: String,
    pub allocation: AllocationReceipt<Element>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StoichiometryReceipt {
    pub table_digest: String,
    pub request: PhaseChemistryRequest,
    pub endmember_allocation: AllocationReceipt<String>,
    pub components: Vec<EndmemberReceipt>,
    pub composition: ElementComposition,
}
impl StoichiometryReceipt {
    /// Recompute every floor, remainder, rounding bit, split and composition.
    /// This verifies arithmetic/model consistency, not geological provenance.
    pub fn validate(&self, table: &ChemistryTable) -> Result<()> {
        if table.decompose(&self.request)? != *self {
            return Err("stoichiometry receipt mismatch".into());
        }
        Ok(())
    }
}
