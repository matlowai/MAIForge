//! Exclusive integer quantity partitions. Optical recipes never supply stock.
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Amounts {
    /// Decimal i128 micrograms, separate from bulk volume.
    pub dry_mass_ug: String,
    pub minerals_ug: BTreeMap<String, String>,
    pub elements_ug: BTreeMap<String, String>,
}
fn integer(s: &str) -> Result<i128, String> {
    let n = s.parse::<i128>().map_err(|_| "invalid i128 quantity")?;
    if n < 0 || n.to_string() != s {
        return Err("negative or noncanonical quantity".into());
    }
    Ok(n)
}
/// Exact proportional partition of an admitted integer account. Weights are
/// independently recorded, quantized equivalent volumes, not deposition counts.
/// Stable ID order receives the division remainder, which is recorded separately.
/// Refuse arithmetic overflow instead of saturating or switching to floats.
pub fn partition(
    total: &str,
    weights: &BTreeMap<String, String>,
) -> Result<(BTreeMap<String, String>, String), String> {
    let total = integer(total)?;
    if weights.is_empty() {
        return Err("empty partition".into());
    }
    let w = weights
        .iter()
        .map(|(k, v)| Ok((k.clone(), integer(v)?)))
        .collect::<Result<BTreeMap<_, _>, String>>()?;
    let sum = w
        .values()
        .try_fold(0i128, |a, &b| a.checked_add(b).ok_or("weight overflow"))?;
    if sum <= 0 {
        return Err("zero weights".into());
    }
    let mut out = BTreeMap::new();
    let mut allocated = 0i128;
    for (id, n) in &w {
        let q = total.checked_mul(*n).ok_or("partition multiply overflow")? / sum;
        allocated = allocated.checked_add(q).ok_or("partition sum overflow")?;
        out.insert(id.clone(), q);
    }
    let mut remainder = total - allocated;
    let receipt = remainder.to_string();
    for (id, &n) in &w {
        if remainder == 0 {
            break;
        }
        if n > 0 {
            *out.get_mut(id).unwrap() += 1;
            remainder -= 1;
        }
    }
    if remainder != 0 {
        return Err("partition remainder invariant".into());
    }
    Ok((
        out.into_iter().map(|(k, v)| (k, v.to_string())).collect(),
        receipt,
    ))
}
pub fn validate_sum(
    parent: &str,
    mut children: impl Iterator<Item = String>,
) -> Result<(), String> {
    let total = children.try_fold(0i128, |a, b| {
        a.checked_add(integer(&b)?)
            .ok_or_else(|| String::from("sum overflow"))
    })?;
    if total != integer(parent)? {
        return Err("child accounts do not sum to parent".into());
    }
    Ok(())
}
impl Amounts {
    pub fn validate(&self) -> Result<(), String> {
        integer(&self.dry_mass_ug)?;
        validate_sum(&self.dry_mass_ug, self.minerals_ug.values().cloned())?;
        validate_sum(&self.dry_mass_ug, self.elements_ug.values().cloned())
    }
}

/// An explicitly selected interpretation, never the homogeneous native record.
pub const ORDERED_REVISION: &str = "pf-inventory/ordered-superposition/1";
pub const ORDERED_LABEL: &str = "conditional proxy: ordered superposition, erosion from top";
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Band {
    pub id: String,
    pub contribution_ids: Vec<String>,
    pub offsets_m: [f64; 2],
    pub bulk_mm3: String,
    pub material_ref: crate::column_envelope::PendingMaterialRef,
    pub redbed: bool,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OrderedStack {
    pub revision: String,
    pub qualification: String,
    pub parent_prov_id: u32,
    pub parent_body_id: String,
    pub parent_bulk_mm3: String,
    pub retained_thickness_m: f64,
    pub bands: Vec<Band>,
}
impl OrderedStack {
    pub fn validate(&self) -> Result<(), String> {
        if self.revision != ORDERED_REVISION
            || self.qualification != ORDERED_LABEL
            || self.parent_body_id.is_empty()
            || !self.retained_thickness_m.is_finite()
            || self.retained_thickness_m <= 0.
            || self.bands.is_empty()
        {
            return Err("invalid conditional stack header".into());
        }
        let mut last = 0.;
        let mut ids = std::collections::BTreeSet::new();
        let mut contributions = std::collections::BTreeSet::new();
        for b in &self.bands {
            b.material_ref.validate()?;
            if b.id.is_empty()
                || !ids.insert(&b.id)
                || b.contribution_ids.is_empty()
                || b.contribution_ids
                    .iter()
                    .any(|id| !contributions.insert(id))
                || b.offsets_m[0] != last
                || !b.offsets_m[1].is_finite()
                || b.offsets_m[1] <= last
                || integer(&b.bulk_mm3)? == 0
            {
                return Err("invalid conditional band partition".into());
            }
            last = b.offsets_m[1];
        }
        if last != self.retained_thickness_m {
            return Err("stack height does not close".into());
        }
        validate_sum(
            &self.parent_bulk_mm3,
            self.bands.iter().map(|b| b.bulk_mm3.clone()),
        )
    }
    /// Half-open contacts: the upper band owns an exact internal contact.
    /// Above the retained top is construction excess, never new stock.
    pub fn band_at(&self, height_above_base_m: f64) -> Option<&Band> {
        self.bands
            .iter()
            .find(|b| height_above_base_m >= b.offsets_m[0] && height_above_base_m < b.offsets_m[1])
    }
}
/// Exact conditional *bulk* ledger. No species, density or ore are invented.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConditionalBulkLedger {
    pub stack_pin: String,
    pub remaining_mm3: BTreeMap<String, String>,
    pub removed_mm3: BTreeMap<String, String>,
}
impl ConditionalBulkLedger {
    pub fn new(stack: &OrderedStack, pin: String) -> Result<Self, String> {
        stack.validate()?;
        crate::planet::validate_digest(&pin)?;
        Ok(Self {
            stack_pin: pin,
            remaining_mm3: stack
                .bands
                .iter()
                .map(|b| (b.id.clone(), b.bulk_mm3.clone()))
                .collect(),
            removed_mm3: stack
                .bands
                .iter()
                .map(|b| (b.id.clone(), "0".into()))
                .collect(),
        })
    }
    pub fn validate(&self, stack: &OrderedStack) -> Result<(), String> {
        stack.validate()?;
        if self.remaining_mm3.len() != stack.bands.len()
            || self.removed_mm3.len() != stack.bands.len()
        {
            return Err("ledger keys".into());
        }
        for b in &stack.bands {
            validate_sum(
                &b.bulk_mm3,
                [
                    self.remaining_mm3
                        .get(&b.id)
                        .ok_or("missing balance")?
                        .clone(),
                    self.removed_mm3
                        .get(&b.id)
                        .ok_or("missing removal")?
                        .clone(),
                ]
                .into_iter(),
            )?;
        }
        Ok(())
    }
    pub fn debit(&mut self, stack: &OrderedStack, id: &str, amount: &str) -> Result<(), String> {
        self.validate(stack)?;
        let n = integer(amount)?;
        if n == 0 {
            return Err("empty debit".into());
        }
        let r = integer(self.remaining_mm3.get(id).ok_or("unknown band")?)?;
        let old = integer(self.removed_mm3.get(id).ok_or("unknown band")?)?;
        if n > r {
            return Err("band exhausted".into());
        }
        let removed = old.checked_add(n).ok_or("debit overflow")?;
        self.remaining_mm3.insert(id.into(), (r - n).to_string());
        self.removed_mm3.insert(id.into(), removed.to_string());
        self.validate(stack)
    }
}

/// Explicit minor extension. The /1 stack, bulk accounts and material IDs stay
/// unchanged; only this selected book revision publishes deposit-scale contacts.
pub const BED_REVISION: &str = "pf-inventory/ordered-superposition/1.1";
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bed {
    pub contribution_id: String,
    /// [bottom, top], metres above this material band's base, not the world datum.
    pub offsets_above_band_base_m: [f64; 2],
    pub thickness_m: f64,
    pub environment_code: u8,
    pub deposition_step: u64,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BandBeds {
    pub band_id: String,
    pub beds: Vec<Bed>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BedHiatus {
    pub step: u64,
    /// Chronological boundary between surviving contributions bracketing the
    /// hiatus. None if erosion removed one side. Never a gap or extra stock;
    /// this projection is not evidence of a physical unconformity surface.
    pub offset_above_parent_base_m: Option<f64>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BedContacts {
    pub bands: Vec<BandBeds>,
    pub hiatus_markers: Vec<BedHiatus>,
}
impl BedContacts {
    pub fn hiatus_offset(&self, stack: &OrderedStack, step: u64) -> Option<f64> {
        let mut previous = None;
        for (band, list) in stack.bands.iter().zip(&self.bands) {
            for bed in &list.beds {
                if previous.is_some_and(|s| s < step) && step < bed.deposition_step {
                    return Some(band.offsets_m[0] + bed.offsets_above_band_base_m[0]);
                }
                previous = Some(bed.deposition_step);
            }
        }
        None
    }
    pub fn validate(&self, stack: &OrderedStack) -> Result<(), String> {
        stack.validate()?;
        if self.bands.len() != stack.bands.len() {
            return Err("bed band coverage".into());
        }
        let mut previous_step = None;
        for (band, list) in stack.bands.iter().zip(&self.bands) {
            if list.band_id != band.id || list.beds.len() != band.contribution_ids.len() {
                return Err("bed contribution coverage".into());
            }
            let mut top = 0.;
            let mut sum = 0.;
            for (bed, id) in list.beds.iter().zip(&band.contribution_ids) {
                let [lo, hi] = bed.offsets_above_band_base_m;
                if &bed.contribution_id != id
                    || lo != top
                    || !hi.is_finite()
                    || hi <= lo
                    || !bed.thickness_m.is_finite()
                    || bed.thickness_m != hi - lo
                    || previous_step.is_some_and(|s| s > bed.deposition_step)
                {
                    return Err("invalid bed order/contact/thickness".into());
                }
                sum += bed.thickness_m;
                top = hi;
                previous_step = Some(bed.deposition_step);
            }
            let height = band.offsets_m[1] - band.offsets_m[0];
            // Exact binary64 closure in the stated bottom-to-top summation order.
            if top != height || sum != height {
                return Err("bed thicknesses do not close exactly to band".into());
            }
        }
        let mut previous_hiatus = None;
        for h in &self.hiatus_markers {
            if previous_hiatus.is_some_and(|s| s >= h.step)
                || h.offset_above_parent_base_m != self.hiatus_offset(stack, h.step)
            {
                return Err("invalid chronological hiatus projection".into());
            }
            previous_hiatus = Some(h.step);
        }
        Ok(())
    }
    /// Half-open bed selection. Returns (band index, bed index, fraction [0,1)).
    /// The caller first evaluates the existing tilted base and finite support.
    pub fn bed_at(&self, stack: &OrderedStack, height_m: f64) -> Option<(usize, usize, f64)> {
        let band = stack
            .bands
            .iter()
            .position(|b| height_m >= b.offsets_m[0] && height_m < b.offsets_m[1])?;
        let base = stack.bands[band].offsets_m[0];
        self.bands[band].beds.iter().enumerate().find_map(|(i, b)| {
            // Compare in published parent coordinates. Subtracting the band
            // base first can round an exact contact into the lower bed.
            let lo = base + b.offsets_above_band_base_m[0];
            let hi = if i + 1 == self.bands[band].beds.len() {
                stack.bands[band].offsets_m[1]
            } else {
                base + b.offsets_above_band_base_m[1]
            };
            (height_m >= lo && height_m < hi).then_some((band, i, (height_m - lo) / (hi - lo)))
        })
    }
}

/// Opt-in realization-only extension. Stack accounts and /1.1 contacts are
/// unchanged. These are optical prior parameters, never chemical reserves.
pub const BED_PROPERTY_REVISION: &str = "pf-inventory/ordered-superposition/1.2";
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BedRealization {
    pub source_digest: String,
    pub paleolatitude_deg: Option<f64>,
    pub aridity_strength: Option<f64>,
    pub hematite_mass_fraction: Option<f64>,
    pub clay_mass_fraction: Option<f64>,
    pub grain_diameter_m: Option<f64>,
    pub original_deposit_thickness_m: f64,
    pub retained_parting_thickness_m: Option<f64>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BedProperties {
    pub revision: String,
    /// Same band/bed ordering as BedContacts. None fields mean no supported
    /// per-deposit delta; the labelled host catalogue prior remains in force.
    pub bands: Vec<Vec<BedRealization>>,
}
impl BedProperties {
    pub fn validate(&self, beds: &BedContacts) -> Result<(), String> {
        if self.revision != BED_PROPERTY_REVISION || self.bands.len() != beds.bands.len() {
            return Err("bed realization revision/coverage".into());
        }
        for (props, band) in self.bands.iter().zip(&beds.bands) {
            if props.len() != band.beds.len() {
                return Err("bed realization count".into());
            }
            for (p, b) in props.iter().zip(&band.beds) {
                crate::planet::validate_digest(&p.source_digest)?;
                if p.paleolatitude_deg
                    .is_some_and(|v| !v.is_finite() || v.abs() > 90.)
                    || p.aridity_strength
                        .is_some_and(|v| !v.is_finite() || !(0. ..=1.).contains(&v))
                    || p.hematite_mass_fraction
                        .is_some_and(|v| !v.is_finite() || !(0.002..=0.007).contains(&v))
                    || p.grain_diameter_m.is_some()
                    || p.clay_mass_fraction.is_some()
                    || !p.original_deposit_thickness_m.is_finite()
                    || p.original_deposit_thickness_m <= 0.
                    || p.original_deposit_thickness_m
                        + 1e-8 * p.original_deposit_thickness_m.max(1.)
                        < b.thickness_m
                    || p.retained_parting_thickness_m
                        .is_some_and(|t| !t.is_finite() || t < 0. || t > b.thickness_m || t > 0.002)
                {
                    return Err("unsupported/invalid bed realization".into());
                }
                // Revision1.2 has no recorded energy, distance, grain or clay
                // data. Future quantitative mapping needs a new revision.
                if p.aridity_strength.is_some() != p.paleolatitude_deg.is_some()
                    || (p.hematite_mass_fraction.is_some() && p.paleolatitude_deg.is_none())
                {
                    return Err("bed proxy lacks deposition latitude".into());
                }
            }
        }
        Ok(())
    }
}
