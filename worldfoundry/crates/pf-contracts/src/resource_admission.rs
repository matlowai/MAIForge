//! D-0141/D-0146/D-0148/D-0149. Additive finite native-resource admission.
//! Canon owns reservations through the existing source registry; no side ledger.
use crate::{excavation::*, player_inventory::{accumulate, add_elements}};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PoolStock {
    pub bulk_mm3: Mass,
    pub minerals_ug: Accounts,
    pub phase_elements_ug: PhaseElements,
}
impl PoolStock {
    pub fn validate(&self) -> Result<()> {
        if self.bulk_mm3.0 <= 0 { return Err("resource_pool_volume_invalid".into()); }
        check_elements(&self.minerals_ug, &self.phase_elements_ug)
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reservation {
    pub pool_id: String,
    pub reservoir_id: String,
    pub role: String,
    pub pool: PoolStock,
    pub reserved: PoolStock,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResourceAdmission {
    pub schema: String,
    pub product_pin: String,
    pub physical_recipe: String,
    pub native_parent_body_id: String,
    pub geometry_role: String,
    /// Parts per million upward from the served native layer base.
    pub interval_ppm: [u32; 2],
    /// Declared raster envelope bulk basis. Native pool volumes remain separate;
    /// raster aliases spread one finite carrier, never repeat its full mass.
    pub geometry_bulk_mm3: Mass,
    pub reservations: Vec<Reservation>,
    pub origin: Value,
    pub material: Value,
}
impl ResourceAdmission {
    pub fn stock(&self) -> Result<PoolStock> {
        if self.schema != "vz.history-resource-admission/1"
            || self.native_parent_body_id.is_empty()
            || self.interval_ppm[0] >= self.interval_ppm[1] || self.interval_ppm[1] > 1_000_000
            || self.geometry_bulk_mm3.0 <= 0 || self.reservations.is_empty()
            || self.origin["schema"] != "pf.resource-origin-facts/v1"
            || self.origin["product_pin"] != self.product_pin
            || self.origin["physical_recipe"] != self.physical_recipe {
            return Err("resource_admission_identity".into());
        }
        crate::planet::validate_digest(&self.product_pin)?;
        crate::planet::validate_digest(&self.physical_recipe)?;
        let mut out = PoolStock { bulk_mm3: Mass(0), minerals_ug: Accounts::new(), phase_elements_ug: PhaseElements::new() };
        let mut seen = BTreeSet::new();
        for r in &self.reservations {
            crate::planet::validate_digest(&r.pool_id)?;
            crate::planet::validate_digest(&r.reservoir_id)?;
            if !seen.insert(&r.pool_id) || !matches!(r.role.as_str(), "host" | "donor" | "background") {
                return Err("resource_pool_duplicate_or_role".into());
            }
            r.pool.validate()?; r.reserved.validate()?;
            if r.reserved.bulk_mm3.0 > r.pool.bulk_mm3.0 { return Err("resource_pool_volume_overdraw".into()); }
            debit(&r.pool.minerals_ug, &r.reserved.minerals_ug)?;
            if r.pool.phase_elements_ug.keys().ne(r.reserved.phase_elements_ug.keys()) { return Err("resource_pool_phase_set".into()); }
            for (m, es) in &r.reserved.phase_elements_ug {
                debit(&r.pool.phase_elements_ug[m], es)?;
                for (e, n) in es {
                    let expected = proportional(r.reserved.bulk_mm3.0, r.pool.phase_elements_ug[m][e].0, r.pool.bulk_mm3.0)?;
                    // D-0147: tolerance for once-quantized conditional pool partitions.
                    if (n.0 - expected).abs() > (expected / 10_000).max(1) {
                        return Err("resource_pool_partition_composition".into());
                    }
                }
            }
            out.bulk_mm3.0 = out.bulk_mm3.0.checked_add(r.reserved.bulk_mm3.0).ok_or("resource_volume_overflow")?;
            accumulate(&mut out.minerals_ug, &r.reserved.minerals_ug)?;
            add_elements(&mut out.phase_elements_ug, &r.reserved.phase_elements_ug)?;
        }
        out.validate()?;
        Ok(out)
    }
    pub fn validate_for(&self, admission: &MiningAdmission) -> Result<PhaseElements> {
        let stock = self.stock()?;
        if self.physical_recipe != admission.world.generation_pin
            || self.geometry_bulk_mm3 != admission.allocated_volume_mm3 || stock.minerals_ug != admission.initial
            || self.material["material_ref"]["realization_id"] != admission.source.material_realization
            || self.material["material_ref"]["body_id"] != admission.source.body_id {
            return Err("resource_admission_source_stock_mismatch".into());
        }
        let mode=self.material["mineral_mass_fractions"].as_object().ok_or("resource_material_mode_missing")?;
        if mode.keys().ne(stock.minerals_ug.keys()) {return Err("resource_material_phase_set".into());}
        let mass=stock.minerals_ug.values().try_fold(0i128,|a,v|a.checked_add(v.0).ok_or("resource_mass_overflow"))? as f64;
        if mass<=0. {return Err("resource_material_empty".into());}
        for (m,n) in &stock.minerals_ug {
            let f=mode[m].as_f64().ok_or("resource_material_mode_invalid")?;
            if !f.is_finite() || (f-n.0 as f64/mass).abs()>1e-10 {return Err("resource_material_stock_mismatch".into());}
        }
        let density=self.material["effective_dry_density_kg_m3"].as_f64().ok_or("resource_material_density_missing")?;
        if !density.is_finite() || (density-mass/self.geometry_bulk_mm3.0 as f64).abs()>density.abs()*1e-10 {return Err("resource_material_density_mismatch".into());}
        Ok(stock.phase_elements_ug)
    }
}
/// All raster aliases and native layer slices reserve from ONE pool. A second
/// product or mode cannot coexist with an old full background owner in a root.
pub fn validate_registry<'a>(states: impl Iterator<Item=&'a MiningState>) -> Result<()> {
    let states: Vec<_> = states.collect();
    let mut pools: BTreeMap<&str, (&PoolStock, PoolStock)> = BTreeMap::new();
    let mut versions = BTreeMap::new();
    let mut roles: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    let mut parents: BTreeMap<&str, Vec<[u32;2]>> = BTreeMap::new();
    for s in &states {
        let Some(r) = s.admission.blast.as_ref().and_then(|b|b.resource.as_ref()) else { continue };
        r.validate_for(&s.admission)?;
        if states.iter().any(|s|s.admission.source.body_id == r.native_parent_body_id) {
            return Err("resource_old_full_parent_double_stock".into());
        }
        let spans = parents.entry(&r.native_parent_body_id).or_default();
        if spans.iter().any(|p|p[0] < r.interval_ppm[1] && r.interval_ppm[0] < p[1]) {
            return Err("resource_parent_geometry_overlap".into());
        }
        spans.push(r.interval_ppm);
        for p in &r.reservations {
            if let Some(v) = versions.insert(&p.reservoir_id, &r.product_pin) {
                if v != &r.product_pin { return Err("resource_reservoir_version_changed".into()); }
            }
            let rs = roles.entry(&p.reservoir_id).or_default(); rs.insert(&p.role);
            if rs.contains("background") && rs.len() > 1 { return Err("resource_background_plus_enriched_double_stock".into()); }
            let (initial, used) = pools.entry(&p.pool_id).or_insert_with(||(&p.pool, PoolStock {bulk_mm3:Mass(0), minerals_ug:Accounts::new(), phase_elements_ug:PhaseElements::new()}));
            if *initial != &p.pool { return Err("resource_pool_identity_changed".into()); }
            used.bulk_mm3.0 = used.bulk_mm3.0.checked_add(p.reserved.bulk_mm3.0).ok_or("resource_volume_overflow")?;
            accumulate(&mut used.minerals_ug, &p.reserved.minerals_ug)?;
            add_elements(&mut used.phase_elements_ug, &p.reserved.phase_elements_ug)?;
            if used.bulk_mm3.0 > initial.bulk_mm3.0 { return Err("resource_shared_pool_exhausted".into()); }
            debit(&initial.minerals_ug, &used.minerals_ug).map_err(|_|"resource_shared_pool_exhausted")?;
            for (m, es) in &used.phase_elements_ug { debit(&initial.phase_elements_ug[m], es).map_err(|_|"resource_shared_pool_exhausted")?; }
        }
    }
    Ok(())
}
