//! Lane CLIMATE-RUNOFF (phase A, 2026-10-20): an OPT-IN climate runoff field for the region chain (T1a skeleton discharge,
//! T1b rivernet areas, T2 pf-detail rill areas). LABEL: pre-signoff, fidelity tier Preview, derived, never canon.
//!
//! PHYSICS. Stream power acts on discharge, not on area: E = K Q^m S^n (the classic Q^m form; Howard 1994, Whipple and
//! Tucker 1999, recalled). The region chain stores "drainage area" everywhere; with this key that area is the EFFECTIVE
//! area A_eff = Q / r_ref = sum over the catchment of a_i R_i / r_ref, where R_i is the local climatic runoff (m/yr) and
//! r_ref is the T0 genesis runoff (terrain assumption A-RAIN-UNITS, 1.0 m/yr on the current worlds). Every existing
//! "area = Q / r_ref" conversion stays as it is, so channel heads, bed slopes k_s A^-theta and channel widths become
//! discharge-driven: a wet cell grows channels from smaller catchments and grades them gentler; a dry cell needs a larger
//! catchment for the same discharge. Without the key every byte is unchanged (A_eff = A, R = r_ref everywhere).
//!
//! RUNOFF RULE (`land_p_minus_e_floor_v1`). Per T0 cell, from a `seed-climate-v1` climate artifact of the SAME world
//! (it must cite the T0 terrain and hydrology digests): R_cell = max(P - E, f P) / 1000 m/yr, f = floor_permille / 1000.
//! The floor is a GEOMORPHIC floor, not a water balance: in drylands the annual P - E understates the channel-forming
//! discharge, because floods there come from short intense storms that make Hortonian overland flow even when the annual
//! budget is nearly closed (ephemeral and flash-flood channels still cut). f = 0.1 is a project assumption awaiting
//! review (recalled order of magnitude of dryland event runoff coefficients, not cited). With seed-climate-v1 the floor is
//! inactive on land cells: that tier caps land E/P at 0.9 (its Budyko cap), so P - E >= 0.1 P already; the floor matters
//! for any climate tier that lets E reach P.
//!
//! INTERPOLATION (`land_only_bilinear_v1`). R(lat, lon) is the bilinear interpolation of R_cell between T0 cell CENTRES
//! with the weights renormalised over LAND cells (marine and lake cells carry their own water balance, not land runoff),
//! falling back to the all-cell bilinear where no corner is land. Continuous across cell edges: no cell-edge steps.
//!
//! LAKES (`lake_net_evaporation_v1`, ADR 043 section 5 amendment B2 with local climate). A T0 lake body keeps water when
//! its R-weighted catchment inflow Q_in (T0 flow graph, land cells R_cell x area) covers its net open-water evaporation
//! (PET - P) x lake area (PET as the open-water evaporation proxy, project assumption); otherwise the rule calls it dry
//! (a playa). The lake's outflow carries max(Q_in + (P - PET) x area, 0). Phase A REPORTS the rule per lake (and uses
//! the outflow in the climate T0 discharge that sets chart-entry inflow); it does not yet turn dry lakes into land
//! (that changes the T1a segment structure; T0 genesis coupling is a later lane).

use crate::t0::{CellWater, T0, NONE32};
use crate::{R, SECONDS_PER_YEAR};
use pf_contracts::artifact::{read_artifact_header, read_field_f32, verify_artifact};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::Path;

/// The climate tier this phase reads.
pub const CLIMATE_TIER_REQUIRED: &str = "seed-climate-v1";
pub const RUNOFF_RULE: &str = "land_p_minus_e_floor_v1: R_cell = max(P - E, f P) / 1000 m/yr per T0 cell from the seed-climate-v1 artifact of the same world (precipitation_mm_per_yr, evaporation_mm_per_yr), f = floor_permille / 1000 (geomorphic floor: dryland channel-forming floods come from intense storms even when the annual P - E is near zero; f a project assumption awaiting review; inactive on seed-climate-v1 land, whose Budyko cap keeps E/P <= 0.9)";
pub const INTERPOLATION_RULE: &str = "land_only_bilinear_v1: R(lat, lon) bilinear between T0 cell centres, weights renormalised over land cells; all-cell bilinear where no corner is land; continuous (no cell-edge steps)";
pub const REFERENCE_RULE: &str = "effective_area_q_over_r_ref_v1: stream power on discharge (E = K Q^m S^n); every chain area is A_eff = Q / r_ref = sum a_i R_i / r_ref with r_ref the T0 genesis runoff (terrain assumption A-RAIN-UNITS); existing area = Q / r_ref conversions unchanged; node / post own area weighted by R(x) / r_ref";
pub const LAKE_RULE: &str = "lake_net_evaporation_v1 (ADR 043 B2 with local climate): a T0 lake body keeps water when its R-weighted catchment inflow covers (PET - P) x lake area (PET = open-water evaporation proxy, project assumption); outflow = max(inflow + (P - PET) x area, 0); phase A reports the rule and uses the outflow in the climate T0 discharge, dry lakes are not yet turned into land";

/// Declared parameters (integers and rule strings; exact canonical JSON). Serialized into a recipe only when the key is on.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ClimateRunoffParamsV1 {
    /// Geomorphic runoff floor as permille of P [default 100: project assumption awaiting review].
    pub floor_permille: i64,
    pub runoff_rule: String,
    pub interpolation_rule: String,
    pub reference_rule: String,
    pub lake_rule: String,
}

impl ClimateRunoffParamsV1 {
    pub fn v1() -> Self {
        ClimateRunoffParamsV1 { floor_permille: 100, runoff_rule: RUNOFF_RULE.into(), interpolation_rule: INTERPOLATION_RULE.into(),
            reference_rule: REFERENCE_RULE.into(), lake_rule: LAKE_RULE.into() }
    }
    pub fn validate(&self) -> R<()> {
        if !(0..=1000).contains(&self.floor_permille) { return Err(format!("climate floor_permille {} outside 0..=1000", self.floor_permille)); }
        Ok(())
    }
    pub fn floor(&self) -> f64 { self.floor_permille as f64 / 1000.0 }
}

/// One T0 lake body under the climate lake rule.
#[derive(Clone, Debug, Serialize)]
pub struct LakeBudget {
    pub body_id: u32,
    pub cells: u32,
    pub area_km2: f64,
    pub inflow_m3_s: f64,
    pub net_evaporation_m3_s: f64,
    pub outflow_m3_s: f64,
    pub keeps_water: bool,
    pub uniform_runoff_inflow_m3_s: f64,
}

/// The loaded climate runoff field on the T0 grid.
pub struct ClimateRunoff {
    pub params: ClimateRunoffParamsV1,
    /// Digest of the verified climate artifact (an upstream of every recipe that reads it).
    pub digest: String,
    pub tier: String,
    /// R_cell (m/yr) of the land rule, evaluated on every cell (interpolation masks non-land cells).
    pub r_cell_m_per_yr: Vec<f64>,
    pub p_mm_per_yr: Vec<f32>,
    pub e_mm_per_yr: Vec<f32>,
    pub pet_mm_per_yr: Vec<f32>,
    /// The T0 genesis runoff (A-RAIN-UNITS), m/yr.
    pub r_ref_m_per_yr: f64,
}

impl ClimateRunoff {
    /// Load and verify `<dir>` (a published `pf.climate/v0` artifact directory, tier seed-climate-v1) for the T0 run `t0`.
    pub fn load(dir: &Path, t0: &T0, params: ClimateRunoffParamsV1) -> R<ClimateRunoff> {
        params.validate()?;
        let digest = verify_artifact::<Value>(dir).map_err(|e| format!("verify climate {}: {e}", dir.display()))?;
        let h = read_artifact_header::<Value>(dir).map_err(|e| format!("climate header: {e}"))?;
        let tier = h.header.semantic.backend.name.clone();
        if tier != CLIMATE_TIER_REQUIRED { return Err(format!("climate artifact tier {tier}; this rule reads {CLIMATE_TIER_REQUIRED}; refusing")); }
        for (role, d) in [("terrain", &t0.terrain_digest), ("hydrology", &t0.hydrology_digest)] {
            if !h.header.semantic.inputs.iter().any(|x| x == d) {
                return Err(format!("climate artifact does not cite the T0 {role} digest {d} (another world or run); refusing"));
            }
        }
        let n = t0.nlat * t0.nlon;
        let f = |name: &str| -> R<Vec<f32>> {
            let v = read_field_f32(dir, name).map_err(|e| format!("climate {name}: {e}"))?;
            if v.len() != n { return Err(format!("climate {name}: {} values, expected {n}", v.len())); }
            if v.iter().any(|x| !x.is_finite()) { return Err(format!("climate {name}: nonfinite value")); }
            Ok(v)
        };
        let (p, e, pet) = (f("precipitation_mm_per_yr")?, f("evaporation_mm_per_yr")?, f("pet_mm_per_yr")?);
        let fl = params.floor();
        let r_cell_m_per_yr = (0..n).map(|c| {
            let (pp, ee) = (p[c] as f64, e[c] as f64);
            (pp - ee).max(fl * pp).max(0.0) / 1000.0
        }).collect();
        Ok(ClimateRunoff { params, digest, tier, r_cell_m_per_yr, p_mm_per_yr: p, e_mm_per_yr: e, pet_mm_per_yr: pet, r_ref_m_per_yr: t0.runoff_m_per_yr })
    }

    /// R (m/yr) at (lat, lon) degrees: land-only bilinear between T0 cell centres (all-cell bilinear where no corner is land).
    pub fn r_at_deg(&self, t0: &T0, lat: f64, lon: f64) -> f64 {
        let dlat = 180.0 / t0.nlat as f64;
        let dlon = 360.0 / t0.nlon as f64;
        let fi = (lat + 90.0) / dlat - 0.5;
        let fj = (lon + 180.0) / dlon - 0.5;
        let (i0, j0) = (fi.floor(), fj.floor());
        let (ti, tj) = (fi - i0, fj - j0);
        let (mut acc, mut wsum, mut all) = (0.0, 0.0, 0.0);
        for (di, wi) in [(0i64, 1.0 - ti), (1, ti)] {
            let ii = (i0 as i64 + di).clamp(0, t0.nlat as i64 - 1) as usize;
            for (dj, wj) in [(0i64, 1.0 - tj), (1, tj)] {
                let jj = (j0 as i64 + dj).rem_euclid(t0.nlon as i64) as usize;
                let c = ii * t0.nlon + jj;
                let w = wi * wj;
                let r = self.r_cell_m_per_yr[c];
                all += w * r;
                if t0.water(c as u32) == CellWater::Land { acc += w * r; wsum += w; }
            }
        }
        if wsum > 0.0 { acc / wsum } else { all }
    }

    /// R / r_ref at (lat, lon) degrees: the factor on a node's or post's own area.
    pub fn weight_at_deg(&self, t0: &T0, lat: f64, lon: f64) -> f64 { self.r_at_deg(t0, lat, lon) / self.r_ref_m_per_yr }

    /// Climate discharge on the T0 flow graph (m3/s per cell): land cells add R_cell x area; lake cells add (P - PET) x area
    /// (may be negative: open-water loss) and pass on at least zero; marine cells pass nothing on. Returns the per-cell
    /// discharge and the lake budgets (`lake_net_evaporation_v1`).
    pub fn t0_discharge(&self, t0: &T0) -> R<(Vec<f64>, Vec<LakeBudget>)> {
        let n = t0.nlat * t0.nlon;
        let mut indeg = vec![0u32; n];
        for c in 0..n { let r = t0.flow_to[c]; if r != NONE32 { indeg[r as usize] += 1; } }
        let mut stack: Vec<usize> = (0..n).filter(|&c| indeg[c] == 0).rev().collect();
        let mut order = Vec::with_capacity(n);
        while let Some(c) = stack.pop() {
            order.push(c);
            let r = t0.flow_to[c];
            if r != NONE32 { let r = r as usize; indeg[r] -= 1; if indeg[r] == 0 { stack.push(r); } }
        }
        if order.len() != n { return Err("T0 flow graph has a cycle".into()); }
        let mut q = vec![0f64; n];
        let mut q_uniform = vec![0f64; n];
        // Inflow into each lake body from non-lake cells (and from lake cells of other bodies).
        let mut lake_in: std::collections::BTreeMap<u32, (f64, f64)> = std::collections::BTreeMap::new();
        for &c in &order {
            let a = t0.cell_area_m2(c as u32);
            match t0.water(c as u32) {
                CellWater::Land => { q[c] += self.r_cell_m_per_yr[c] * a / SECONDS_PER_YEAR; q_uniform[c] += self.r_ref_m_per_yr * a / SECONDS_PER_YEAR; }
                CellWater::Lake => {
                    q[c] += (self.p_mm_per_yr[c] as f64 - self.pet_mm_per_yr[c] as f64) / 1000.0 * a / SECONDS_PER_YEAR;
                    q_uniform[c] += self.r_ref_m_per_yr * a / SECONDS_PER_YEAR;
                }
                CellWater::Marine => { continue; }
            }
            if t0.water(c as u32) == CellWater::Lake { q[c] = q[c].max(0.0); }
            let r = t0.flow_to[c];
            if r == NONE32 { continue; }
            let r = r as usize;
            if t0.water(r as u32) == CellWater::Lake && t0.water_body_id[r] != t0.water_body_id[c] {
                let e = lake_in.entry(t0.water_body_id[r]).or_default();
                e.0 += q[c];
                e.1 += q_uniform[c];
            }
            q[r] += q[c];
            q_uniform[r] += q_uniform[c];
        }
        let mut lakes = Vec::new();
        for (&body, &(kind, _)) in &t0.bodies {
            if kind != pf_contracts::domain::WaterBodyKindV0::Lake { continue; }
            let cells: Vec<usize> = (0..n).filter(|&c| t0.water_body_id[c] == body).collect();
            if cells.is_empty() { continue; }
            let area: f64 = cells.iter().map(|&c| t0.cell_area_m2(c as u32)).sum();
            let net_evap: f64 = cells.iter().map(|&c| (self.pet_mm_per_yr[c] as f64 - self.p_mm_per_yr[c] as f64) / 1000.0 * t0.cell_area_m2(c as u32) / SECONDS_PER_YEAR).sum();
            let (qin, quni) = lake_in.get(&body).copied().unwrap_or((0.0, 0.0));
            lakes.push(LakeBudget { body_id: body, cells: cells.len() as u32, area_km2: area / 1e6, inflow_m3_s: qin, net_evaporation_m3_s: net_evap,
                outflow_m3_s: (qin - net_evap).max(0.0), keeps_water: qin >= net_evap, uniform_runoff_inflow_m3_s: quni });
        }
        Ok((q, lakes))
    }
}
