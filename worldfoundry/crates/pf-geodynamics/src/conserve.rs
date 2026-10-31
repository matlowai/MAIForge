//! SEED-CARRY-CONSERVE (lane CV, 2026-10-13): `continental-parcel-carry-v1`
//! (a conservative continental carry) and `contact-length-pole-share-v1`
//! (continent-continent suturing). Pre-registered in
//! `unpublished-design-record`.
//!
//! Why: the legacy carry (`kinematic.rs`, A-ADVECTION) resamples plate id,
//! kind and thickness every step from the nearest back-traced material. That
//! is not conservative: continental parcels fall through, foreign plate ids
//! speckle inside continents and punch oceanic holes into them, and the loss
//! scales with the number of steps (lane SC: seed 42 L0 continental area at
//! 300 My 0.109 / 0.219 / 0.249 at dt 1 / 2 / 4 My).
//!
//! `continental-parcel-carry-v1`:
//! - continental crust is a set of Lagrangian parcels (position, plate,
//!   AREA, thickness, sediment, lineage fields, birth time), seeded one per
//!   continental cell at its centroid; every step each parcel rotates
//!   rigidly with its plate (exact, never resampled), so the parcel area is
//!   conserved by construction;
//! - the grid is a PROJECTION: each parcel spreads its area over nearby
//!   cells with a tent kernel of radius one row spacing; a cell whose
//!   deposited area is at least half the cell is continental, takes the
//!   plate with the largest deposit and the deposit-weighted mean column;
//!   oceanic crust keeps the legacy Eulerian carry, which the projection
//!   overrides;
//! - at the end of the step `write_back` hands every change the step made to
//!   a projected continental cell to the parcels in proportion to their
//!   deposit (sum of parcel area x change = sum of cell area x change);
//! - the only area terms: rift births that end the step continental become
//!   new parcels (+), and collision consumption (-): where two plates each
//!   cover half a cell and the total exceeds `1 + overlap_cell_fraction`,
//!   the smaller depositors' parcels whose nearest cell it is are consumed
//!   and their volume is the belt's collision supply;
//! - plate-id pieces smaller than `fragment_min_area_fraction` of the
//!   sphere fold into the plate they touch most (parcels follow).
//!
//! `contact-length-pole-share-v1` (`suture_pass`): when the converging
//! continent-continent contact between two plates reaches the declared
//! length, the smaller plate group takes the larger's Euler pole.
//!
//! Determinism: parcels keep a fixed order (seed order, then births in
//! ascending cell order); every reduction runs in ascending parcel or cell
//! order; ties break to the lowest plate / parcel / cell id; no HashMap.

use crate::boundaries::normal_rate_at;
use crate::kinematic::{nearest_cell, AdvanceLedger, CrustState, KinematicError, NO_SOURCE};
use crate::rotation::{great_circle_rad, rotate_by_euler};
use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;
use std::collections::BTreeMap;

/// Constants of `continental-parcel-carry-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParcelParams {
    pub overlap_cell_fraction: f64,
    pub fragment_min_area_fraction: f64,
}

impl ParcelParams {
    pub fn validate(&self) -> Result<(), String> {
        if !(self.overlap_cell_fraction > 0.0 && self.overlap_cell_fraction.is_finite()) {
            return Err(format!("parcel overlap fraction {}", self.overlap_cell_fraction));
        }
        if !(self.fragment_min_area_fraction >= 0.0 && self.fragment_min_area_fraction < 1.0) {
            return Err(format!(
                "fragment min area fraction {}",
                self.fragment_min_area_fraction
            ));
        }
        Ok(())
    }
}

/// Proposed recording-only carrier handles. Not material or grain ancestry.
/// IDs are local to a recipe namespace supplied by the publisher.
pub const MAX_CARRIER_ID_RECORDS: u64 = 10_000_000;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarrierBirth {
    pub id: u64,
    pub model_time_my: f64,
    pub event_cell: usize,
    pub plate: u32,
    pub initial: bool,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarrierRetirement {
    pub id: u64,
    pub model_time_my: f64,
    pub event_cell: usize,
    pub plate: u32,
}
#[derive(Debug, Clone, PartialEq)]
pub struct CarrierIdentity {
    record_limit: u64,
    grid_cells: usize,
    started_at_my: f64,
    initial_count: u64,
    next_id: u64,
    ids: Vec<u64>,
    births: Vec<CarrierBirth>,
    retirements: Vec<CarrierRetirement>,
    retired_lookup: std::collections::BTreeSet<u64>,
}
impl CarrierIdentity {
    pub fn validate_limit(limit: u64) -> Result<(), String> {
        if !(1..=MAX_CARRIER_ID_RECORDS).contains(&limit) {
            return Err("carrier record limit must be 1..=10000000".into());
        }
        Ok(())
    }
    pub fn initial(
        grid_cells: usize,
        time: f64,
        cells: &[(usize, u32)],
        limit: u64,
    ) -> Result<Self, String> {
        Self::validate_limit(limit)?;
        if grid_cells == 0 || !time.is_finite() || cells.len() as u64 > limit {
            return Err("invalid initial carrier identity support/budget".into());
        }
        if cells.windows(2).any(|w| w[0].0 >= w[1].0) || cells.iter().any(|(c, _)| *c >= grid_cells)
        {
            return Err("initial carrier cells must be sorted, unique and in grid".into());
        }
        let births = cells
            .iter()
            .enumerate()
            .map(|(i, (cell, plate))| CarrierBirth {
                id: i as u64,
                model_time_my: time,
                event_cell: *cell,
                plate: *plate,
                initial: true,
            })
            .collect();
        Ok(Self {
            record_limit: limit,
            grid_cells,
            started_at_my: time,
            initial_count: cells.len() as u64,
            next_id: cells.len() as u64,
            ids: (0..cells.len() as u64).collect(),
            births,
            retirements: Vec::new(),
            retired_lookup: std::collections::BTreeSet::new(),
        })
    }
    pub fn ids(&self) -> &[u64] {
        &self.ids
    }
    pub fn births(&self) -> &[CarrierBirth] {
        &self.births
    }
    pub fn retirements(&self) -> &[CarrierRetirement] {
        &self.retirements
    }
    pub fn next_id(&self) -> u64 {
        self.next_id
    }
    pub fn record_limit(&self) -> u64 {
        self.record_limit
    }
    pub fn initial_count(&self) -> u64 {
        self.initial_count
    }
    pub fn check_context(&self, grid_cells: usize, active: usize) -> Result<(), String> {
        if self.grid_cells != grid_cells
            || self.ids.len() != active
            || self.births.len() != self.next_id as usize
            || self.births.len() != self.ids.len() + self.retirements.len()
        {
            return Err("carrier identity/grid/count mismatch".into());
        }
        Ok(())
    }
    /// Conservative logical-record reservation; not an allocator/RSS bound.
    /// Call before any physical step mutation with a bound on all births/retirements.
    pub fn preflight(&self, births: usize, retirements: usize) -> Result<(), String> {
        let required = self
            .births
            .len()
            .checked_add(self.retirements.len())
            .and_then(|n| n.checked_add(births))
            .and_then(|n| n.checked_add(retirements))
            .ok_or("carrier record count overflow")?;
        if required as u64 > self.record_limit || self.next_id.checked_add(births as u64).is_none()
        {
            return Err("carrier record budget exhausted before mutation".into());
        }
        Ok(())
    }
    pub fn reserve(&mut self, births: usize, retirements: usize) -> Result<(), String> {
        self.preflight(births, retirements)?;
        self.births
            .try_reserve(births)
            .map_err(|_| "carrier birth allocation")?;
        self.retirements
            .try_reserve(retirements)
            .map_err(|_| "carrier retirement allocation")?;
        self.ids
            .try_reserve(births)
            .map_err(|_| "carrier active-ID allocation")?;
        Ok(())
    }
    fn time_cell(&self, time: f64, cell: usize) -> Result<(), String> {
        if !time.is_finite() || time < self.started_at_my || cell >= self.grid_cells {
            return Err("invalid carrier event time/cell".into());
        }
        Ok(())
    }
    pub(crate) fn allocate(&mut self, time: f64, cell: usize, plate: u32) -> Result<u64, String> {
        self.time_cell(time, cell)?;
        if self.births.last().is_some_and(|b| time < b.model_time_my) {
            return Err("carrier creation time reversed".into());
        }
        self.reserve(1, 0)?;
        let id = self.next_id;
        self.births.push(CarrierBirth {
            id,
            model_time_my: time,
            event_cell: cell,
            plate,
            initial: false,
        });
        self.ids.push(id);
        self.next_id += 1;
        Ok(id)
    }
    pub(crate) fn retire(
        &mut self,
        index: usize,
        time: f64,
        cell: usize,
        plate: u32,
    ) -> Result<(), String> {
        self.time_cell(time, cell)?;
        let id = *self
            .ids
            .get(index)
            .ok_or("carrier retirement index out of range")?;
        if self.retired_lookup.contains(&id)
            || time < self.births[id as usize].model_time_my
            || self
                .retirements
                .last()
                .is_some_and(|r| time < r.model_time_my)
        {
            return Err("duplicate or time-reversed carrier retirement".into());
        }
        self.reserve(0, 1)?;
        self.retirements.push(CarrierRetirement {
            id,
            model_time_my: time,
            event_cell: cell,
            plate,
        });
        self.retired_lookup.insert(id);
        Ok(())
    }
    pub(crate) fn retain(&mut self, keep: &[bool]) -> Result<(), String> {
        if keep.len() != self.ids.len()
            || self
                .ids
                .iter()
                .zip(keep)
                .any(|(id, k)| *k == self.retired_lookup.contains(id))
        {
            return Err("carrier retain mask lacks matching retirement".into());
        }
        let mut i = 0;
        self.ids.retain(|_| {
            let k = keep[i];
            i += 1;
            k
        });
        Ok(())
    }
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({"schema":"vz.carrier-identities/v1","record_limit":self.record_limit,
            "grid_cells":self.grid_cells,"started_at_my":self.started_at_my,
            "initial_count":self.initial_count,"next_id":self.next_id,"ids":self.ids,
            "births":self.births.iter().map(|b|(b.id,b.model_time_my,b.event_cell,b.plate,b.initial)).collect::<Vec<_>>(),
            "retirements":self.retirements.iter().map(|r|(r.id,r.model_time_my,r.event_cell,r.plate)).collect::<Vec<_>>()})
    }
    pub fn from_json(value: serde_json::Value) -> Result<Self, String> {
        let keys = [
            "schema",
            "record_limit",
            "grid_cells",
            "started_at_my",
            "initial_count",
            "next_id",
            "ids",
            "births",
            "retirements",
        ];
        let object = value
            .as_object()
            .ok_or("carrier archive must be an object")?;
        if object.len() != keys.len()
            || keys.iter().any(|k| !object.contains_key(*k))
            || value["schema"] != "vz.carrier-identities/v1"
        {
            return Err("carrier archive schema/fields mismatch".into());
        }
        let limit = value["record_limit"]
            .as_u64()
            .ok_or("carrier limit is not an integer")?;
        Self::validate_limit(limit)?;
        for name in ["ids", "births", "retirements"] {
            if value[name]
                .as_array()
                .ok_or("carrier array required")?
                .len() as u64
                > limit
            {
                return Err("carrier archive exceeds declared record limit".into());
            }
        }
        let ids: Vec<u64> =
            serde_json::from_value(value["ids"].clone()).map_err(|e| e.to_string())?;
        let births: Vec<(u64, f64, usize, u32, bool)> =
            serde_json::from_value(value["births"].clone()).map_err(|e| e.to_string())?;
        let retired: Vec<(u64, f64, usize, u32)> =
            serde_json::from_value(value["retirements"].clone()).map_err(|e| e.to_string())?;
        let mut result = Self {
            record_limit: limit,
            grid_cells: serde_json::from_value(value["grid_cells"].clone())
                .map_err(|e| e.to_string())?,
            started_at_my: value["started_at_my"]
                .as_f64()
                .ok_or("carrier start time required")?,
            initial_count: value["initial_count"]
                .as_u64()
                .ok_or("carrier initial count required")?,
            next_id: value["next_id"]
                .as_u64()
                .ok_or("carrier allocator required")?,
            ids,
            births: births
                .into_iter()
                .map(|(id, t, c, p, i)| CarrierBirth {
                    id,
                    model_time_my: t,
                    event_cell: c,
                    plate: p,
                    initial: i,
                })
                .collect(),
            retirements: retired
                .into_iter()
                .map(|(id, t, c, p)| CarrierRetirement {
                    id,
                    model_time_my: t,
                    event_cell: c,
                    plate: p,
                })
                .collect(),
            retired_lookup: std::collections::BTreeSet::new(),
        };
        if result.grid_cells == 0
            || !result.started_at_my.is_finite()
            || result.initial_count > result.next_id
            || result.births.len() as u64 != result.next_id
            || result.births.len() + result.retirements.len() > limit as usize
            || result.ids.len() + result.retirements.len() != result.births.len()
        {
            return Err("carrier archive count identity failed".into());
        }
        let mut previous = result.started_at_my;
        for (i, b) in result.births.iter().enumerate() {
            result.time_cell(b.model_time_my, b.event_cell)?;
            if b.id != i as u64
                || b.initial != (b.id < result.initial_count)
                || b.model_time_my < previous
                || (b.initial && b.model_time_my != result.started_at_my)
            {
                return Err("carrier birth sequence invalid".into());
            }
            previous = b.model_time_my;
        }
        previous = result.started_at_my;
        for r in &result.retirements {
            result.time_cell(r.model_time_my, r.event_cell)?;
            if r.id >= result.next_id
                || r.model_time_my < result.births[r.id as usize].model_time_my
                || r.model_time_my < previous
                || !result.retired_lookup.insert(r.id)
            {
                return Err("carrier retirement sequence invalid".into());
            }
            previous = r.model_time_my;
        }
        let mut seen = std::collections::BTreeSet::new();
        if result.ids.iter().any(|id| {
            *id >= result.next_id || result.retired_lookup.contains(id) || !seen.insert(*id)
        }) {
            return Err("duplicate, retired or unknown active carrier ID".into());
        }
        Ok(result)
    }
}
#[cfg(test)]
mod carrier_identity_tests {
    use super::*;
    fn sample() -> CarrierIdentity {
        CarrierIdentity::initial(20, 0., &[(1, 2), (4, 2), (9, 3)], 20).unwrap()
    }
    #[test]
    fn compaction_does_not_rename_survivors_and_ids_are_not_reused() {
        let mut x = sample();
        x.reserve(2, 3).unwrap();
        x.retire(1, 2., 5, 2).unwrap();
        x.retain(&[true, false, true]).unwrap();
        assert_eq!(x.ids(), &[0, 2]);
        assert_eq!(x.allocate(2., 10, 4).unwrap(), 3);
        assert_eq!(x.ids(), &[0, 2, 3]);
        assert_eq!(x.next_id(), 4);
        x.check_context(20, 3).unwrap();
        assert_eq!(x.births().len() - x.retirements().len(), x.ids().len());
    }
    #[test]
    fn configured_refusal_preserves_logical_state() {
        let mut x = sample();
        let before = x.clone();
        assert!(x.reserve(18, 0).is_err());
        assert_eq!(x, before);
        assert!(x.allocate(f64::NAN, 2, 0).is_err());
        assert_eq!(x, before);
        assert!(x.retire(10, 2., 1, 0).is_err());
        assert_eq!(x, before);
        assert!(x.retain(&[true, false, true]).is_err());
        assert_eq!(x, before);
    }
    #[test]
    fn checkpoint_roundtrip_and_continued_allocation_match() {
        let mut a = sample();
        a.retire(0, 2., 2, 2).unwrap();
        a.retain(&[false, true, true]).unwrap();
        let encoded = serde_json::to_vec(&a.to_json()).unwrap();
        let mut b = CarrierIdentity::from_json(serde_json::from_slice(&encoded).unwrap()).unwrap();
        assert_eq!(a, b);
        for x in [&mut a, &mut b] {
            x.allocate(4., 12, 8).unwrap();
            x.retire(1, 6., 13, 3).unwrap();
            x.retain(&[true, false, true]).unwrap();
        }
        assert_eq!(a.to_json(), b.to_json());
    }
    #[test]
    fn malformed_checkpoint_is_refused() {
        let original = sample().to_json();
        for (field, value) in [
            ("ids", serde_json::json!([0, 0, 2])),
            ("next_id", serde_json::json!(4)),
            ("record_limit", serde_json::json!(2)),
        ] {
            let mut bad = original.clone();
            bad[field] = value;
            assert!(CarrierIdentity::from_json(bad).is_err());
        }
        let mut bad = original.clone();
        bad["births"][0][2] = serde_json::json!(20);
        assert!(CarrierIdentity::from_json(bad).is_err());
        let mut bad = original;
        bad["material_ancestry"] = serde_json::json!(true);
        assert!(CarrierIdentity::from_json(bad).is_err());
    }
    #[test]
    fn duplicate_retirement_and_bad_context_are_refused() {
        let mut x = sample();
        x.retire(0, 2., 2, 2).unwrap();
        let before = x.clone();
        assert!(x.retire(0, 2., 2, 2).is_err());
        assert_eq!(x, before);
        x.retain(&[false, true, true]).unwrap();
        assert!(x.check_context(21, 2).is_err());
        assert!(x.check_context(20, 3).is_err());
        assert!(CarrierIdentity::initial(5, 0., &[(2, 0), (2, 1)], 10).is_err());
    }
}


/// The continental parcels (struct of arrays, one entry per parcel).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ParcelLayer {
    pub pos: Vec<[f64; 3]>,
    pub plate: Vec<PlateId>,
    pub area_m2: Vec<f64>,
    pub thickness_m: Vec<f64>,
    pub sediment_m: Vec<f64>,
    pub eroded_total_m: Vec<f64>,
    pub arc_added_m: Vec<f64>,
    pub shortening_added_m: Vec<f64>,
    pub birth_t_my: Vec<f64>,
    /// SEED-STRATIGRAPHY (lane ST, 2026-10-13) `parcel-stacks-v1`: one
    /// stack per parcel (same order); `None` under every other law.
    pub strata: Option<Vec<crate::stratigraphy::StrataStack>>,
    /// Recording-only lifecycle. Not material or grain ancestry.
    pub identity: Option<CarrierIdentity>,
    pub history: Option<crate::carrier_history::StackHistory>,
    pub layer_provenance: Option<Vec<pf_contracts::layer_provenance::StackProvenance>>,
}

impl ParcelLayer {
    pub fn len(&self) -> usize {
        self.pos.len()
    }
    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    fn push_cell(&mut self, state: &CrustState, cents: &[[f64; 3]], areas: &[f64], c: usize) -> Result<(), KinematicError> {
        if let Some(ids) = self.identity.as_mut() {
            ids.allocate(state.t_my, c, state.plate[c].0).map_err(KinematicError::Grid)?;
        }
        if let Some(p)=self.layer_provenance.as_mut(){p.push(Default::default());}
        self.pos.push(cents[c]);
        self.plate.push(state.plate[c]);
        self.area_m2.push(areas[c]);
        self.thickness_m.push(state.thickness_m[c]);
        self.sediment_m.push(state.sediment_m[c]);
        self.eroded_total_m.push(state.eroded_total_m[c]);
        self.arc_added_m.push(state.arc_added_m[c]);
        self.shortening_added_m.push(state.shortening_added_m[c]);
        self.birth_t_my.push(state.birth_t_my[c]);
        // SEED-STRATIGRAPHY: a placeholder stack (the stratigraphy pass of
        // `write_back_strata` gives births their real stack).
        if let Some(v) = self.strata.as_mut() {
            v.push(crate::stratigraphy::StrataStack::basement_only(
                state.thickness_m[c] + state.sediment_m[c],
                crate::stratigraphy::LITH_BASEMENT,
                state.t_my,
            ));
        }
        Ok(())
    }

    /// Initialize recording at seeding, never retrospectively infer formation ages.
    pub fn seed_with_identity(grid: &GridSpec, state: &CrustState, limit: u64) -> Result<Self, KinematicError> {
        let mut layer = Self::seed(grid, state)?;
        let cells: Vec<_> = (0..state.cells())
            .filter(|&c| state.kind[c] == CrustKindCodeV0::Continental)
            .map(|c| (c, state.plate[c].0)).collect();
        layer.identity = Some(CarrierIdentity::initial(grid.cells(), state.t_my, &cells, limit)
            .map_err(KinematicError::Grid)?);
        Ok(layer)
    }
    pub fn reserve_identity_events(&mut self, grid_cells: usize, births: usize, retirements: usize) -> Result<(), KinematicError> {
        let active = self.len();
        if let Some(ids) = self.identity.as_mut() {
            ids.check_context(grid_cells, active).map_err(KinematicError::Grid)?;
            ids.reserve(births, retirements).map_err(KinematicError::Grid)?;
        }
        Ok(())
    }

    /// One parcel per continental cell of `state`, at the centroid, with
    /// the cell's area and column (ascending cell order).
    pub fn seed(grid: &GridSpec, state: &CrustState) -> Result<Self, KinematicError> {
        let cents = crate::kinematic::centroids(grid)?;
        let areas = crate::kinematic::cell_areas(grid)?;
        let mut l = ParcelLayer::default();
        for c in 0..state.cells() {
            if state.kind[c] == CrustKindCodeV0::Continental {
                l.push_cell(state, &cents, &areas, c)?;
            }
        }
        Ok(l)
    }

    fn retain(&mut self, keep: &[bool]) -> Result<(), KinematicError> {
        if let Some(ids) = self.identity.as_mut() {
            ids.retain(keep).map_err(KinematicError::Grid)?;
        }
        fn r<T: Copy>(v: &mut Vec<T>, keep: &[bool]) {
            let mut k = 0;
            v.retain(|_| {
                let x = keep[k];
                k += 1;
                x
            });
        }
        if let Some(v)=self.layer_provenance.as_mut(){let mut k=0;v.retain(|_|{let b=keep[k];k+=1;b});}
        r(&mut self.pos, keep);
        r(&mut self.plate, keep);
        r(&mut self.area_m2, keep);
        r(&mut self.thickness_m, keep);
        r(&mut self.sediment_m, keep);
        r(&mut self.eroded_total_m, keep);
        r(&mut self.arc_added_m, keep);
        r(&mut self.shortening_added_m, keep);
        r(&mut self.birth_t_my, keep);
        if let Some(v) = self.strata.as_mut() {
            let mut k = 0;
            v.retain(|_| {
                let x = keep[k];
                k += 1;
                x
            });
        }
        Ok(())
    }

    /// Total parcel area (m^2), crust volume (m^3) and sediment (m^3).
    pub fn totals(&self) -> (f64, f64, f64) {
        let (mut a, mut v, mut s) = (0.0, 0.0, 0.0);
        for k in 0..self.len() {
            a += self.area_m2[k];
            v += self.area_m2[k] * self.thickness_m[k];
            s += self.area_m2[k] * self.sediment_m[k];
        }
        (a, v, s)
    }
}

/// What the advance's projection wrote into continental cells (per cell;
/// `continental[c]` false elsewhere), kept for `write_back`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ParcelProjection {
    pub continental: Vec<bool>,
    pub thickness_m: Vec<f64>,
    pub sediment_m: Vec<f64>,
    pub eroded_total_m: Vec<f64>,
    pub arc_added_m: Vec<f64>,
    pub shortening_added_m: Vec<f64>,
}

/// Per-step bookkeeping of the parcel carry inside the advance.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ParcelStepReport {
    pub parcels_before: u64,
    pub area_before_m2: f64,
    pub volume_before_m3: f64,
    pub sediment_before_m3: f64,
    /// Collision consumption (the declared area sink).
    pub consumed_parcels: u64,
    pub consumed_area_m2: f64,
    pub consumed_m3: f64,
    pub consumed_sediment_m3: f64,
    pub overlap_cells: u64,
    /// The projection: continental cells, their area and grid crust volume.
    pub projected_cells: u64,
    pub projected_area_m2: f64,
    pub projected_volume_m3: f64,
    /// Oceanic columns (and births) the projection overrode.
    pub overridden_oceanic_cells: u64,
    pub overridden_births: u64,
    /// Cells the Eulerian carry left continental that the projection does
    /// not: filled from an adjacent oceanic cell of the same plate, or born.
    pub vacated_copied: u64,
    pub vacated_births: u64,
    /// Plate-id fragments folded, cells relabelled, parcels relabelled.
    pub fragments_folded: u64,
    pub fragment_cells: u64,
    pub fragment_parcels: u64,
    /// SEED-STRATIGRAPHY: consumed stack volume per lithology (m^3; empty
    /// when the parcels carry no stacks).
    pub consumed_strata_m3: Vec<f64>,
}

/// End-of-step write-back bookkeeping.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct WriteBackReport {
    /// sum(cell area x change) handed to parcels (crust m^3, sediment m^3).
    pub writeback_m3: f64,
    pub writeback_sediment_m3: f64,
    /// Rift births that ended continental: new parcels (the area source).
    pub created_parcels: u64,
    pub created_area_m2: f64,
    pub created_m3: f64,
    pub created_sediment_m3: f64,
    /// Thickness / sediment driven below zero and clamped (m^3 added back).
    pub clamped_m3: f64,
    pub clamped_sediment_m3: f64,
    /// Projected continental cells that ended the step oceanic (not
    /// expected; their change is not written back).
    pub lost_cells: u64,
    pub parcels_after: u64,
    pub area_after_m2: f64,
    pub volume_after_m3: f64,
    pub sediment_after_m3: f64,
}

/// The tent kernel of `continental-parcel-carry-v1`.
pub(crate) struct Kernel<'a> {
    grid: &'a GridSpec,
    cents: &'a [[f64; 3]],
    areas: &'a [f64],
    radius_rad: f64,
    dlon: f64,
    row_cos: Vec<f64>,
    row_lat: Vec<f64>,
    /// SEED-PERF (lane SP): squared chord beyond which the distance
    /// provably exceeds `radius_rad` (for near-unit parcel positions).
    skip_chord2: f64,
}

impl<'a> Kernel<'a> {
    pub(crate) fn new(grid: &'a GridSpec, cents: &'a [[f64; 3]], areas: &'a [f64]) -> Result<Self, KinematicError> {
        let dlat = std::f64::consts::PI / grid.nlat() as f64;
        let mut row_cos = Vec::with_capacity(grid.nlat());
        let mut row_lat = Vec::with_capacity(grid.nlat());
        for i in 0..grid.nlat() {
            let lat = grid
                .center_lat_deg(i)
                .map_err(|e| KinematicError::Grid(e.to_string()))?
                .to_radians();
            row_lat.push(lat);
            row_cos.push(lat.cos().max(1.0e-9));
        }
        Ok(Self {
            grid,
            cents,
            areas,
            radius_rad: dlat,
            skip_chord2: if dlat < 3.0 {
                let h = (0.5 * dlat).sin();
                4.0 * h * h * (1.0 + 1.0e-9) + 1.0e-13
            } else {
                f64::INFINITY
            },
            dlon: 2.0 * std::f64::consts::PI / grid.nlon() as f64,
            row_cos,
            row_lat,
        })
    }

    /// Normalised weights (sum 1) of the parcel at `p` over cells, appended
    /// to `out` in ascending cell order.
    pub(crate) fn weights(&self, p: [f64; 3], out: &mut Vec<(u32, f64)>) -> Result<(), KinematicError> {
        let c0 = nearest_cell(self.grid, p)?;
        let nlon = self.grid.nlon() as i64;
        let nlat = self.grid.nlat() as i64;
        let (i0, j0) = ((c0 / self.grid.nlon()) as i64, (c0 % self.grid.nlon()) as i64);
        let lat_p = p[2].clamp(-1.0, 1.0).asin();
        let start = out.len();
        // The chord skip assumes |p| = 1 to 1e-10 (else exact everywhere).
        let skip = if (p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0).abs() <= 1.0e-10 {
            self.skip_chord2
        } else {
            f64::INFINITY
        };
        let mut sum = 0.0;
        for di in -1..=1i64 {
            let i = i0 + di;
            if i < 0 || i >= nlat {
                continue;
            }
            let iu = i as usize;
            if (self.row_lat[iu] - lat_p).abs() >= self.radius_rad {
                continue;
            }
            let half = (self.radius_rad / (self.dlon * self.row_cos[iu])).ceil() as i64 + 1;
            let mut cols: Vec<i64> = if 2 * half + 1 >= nlon {
                (0..nlon).collect()
            } else {
                (-half..=half).map(|dj| (j0 + dj).rem_euclid(nlon)).collect()
            };
            cols.sort_unstable();
            for j in cols {
                let c = (i * nlon + j) as usize;
                // SEED-PERF (lane SP): a cell provably farther than the
                // radius (squared chord margin as in
                // `kinematic::MaterialSearch`) has weight 0 and is never
                // pushed; skip its exact distance.
                let q = self.cents[c];
                let (x, y, z) = (q[0] - p[0], q[1] - p[1], q[2] - p[2]);
                if x * x + y * y + z * z > skip {
                    continue;
                }
                let d = great_circle_rad(self.cents[c], p)?;
                let w = (1.0 - d / self.radius_rad).max(0.0) * self.areas[c];
                if w > 0.0 {
                    out.push((c as u32, w));
                    sum += w;
                }
            }
        }
        if !(sum > 0.0) {
            out.truncate(start);
            out.push((c0 as u32, 1.0));
            return Ok(());
        }
        out[start..].sort_unstable_by_key(|x| x.0);
        for x in &mut out[start..] {
            x.1 /= sum;
        }
        Ok(())
    }

    /// CSR weights of every parcel: (offsets, entries).
    fn all(&self, layer: &ParcelLayer) -> Result<(Vec<usize>, Vec<(u32, f64)>), KinematicError> {
        // SEED-PERF (lane SP): blocks of parcels in parallel, concatenated
        // in parcel order (the first error in parcel order is returned).
        const BLOCK: usize = 512;
        let nb = layer.len().div_ceil(BLOCK);
        let blocks = crate::kinematic::par::map_heavy(nb, |b| -> (Vec<usize>, Vec<(u32, f64)>, Option<KinematicError>) {
            let mut lens = Vec::with_capacity(BLOCK);
            let mut ent = Vec::with_capacity(BLOCK * 8);
            for k in b * BLOCK..((b + 1) * BLOCK).min(layer.len()) {
                if let Err(e) = self.weights(layer.pos[k], &mut ent) {
                    return (lens, ent, Some(e));
                }
                lens.push(ent.len());
            }
            (lens, ent, None)
        });
        let mut off = Vec::with_capacity(layer.len() + 1);
        let mut ent = Vec::with_capacity(layer.len() * 8);
        off.push(0);
        for (lens, e, err) in blocks {
            let base = ent.len();
            off.extend(lens.iter().map(|l| base + l));
            ent.extend(e);
            if let Some(err) = err {
                return Err(err);
            }
        }
        Ok((off, ent))
    }
}

fn plate_index(layer: &ParcelLayer, extra: &[PlateId]) -> BTreeMap<PlateId, usize> {
    let mut ids: Vec<PlateId> = layer.plate.clone();
    ids.extend_from_slice(extra);
    ids.sort();
    ids.dedup();
    ids.into_iter().enumerate().map(|(i, p)| (p, i)).collect()
}

/// Per-plate deposits (m^2) per cell: `dep[plate index][cell]`.
fn deposits(
    layer: &ParcelLayer,
    off: &[usize],
    ent: &[(u32, f64)],
    pidx: &BTreeMap<PlateId, usize>,
    n: usize,
) -> Vec<Vec<f64>> {
    let mut dep = vec![vec![0.0f64; n]; pidx.len()];
    for k in 0..layer.len() {
        let pi = pidx[&layer.plate[k]];
        for &(c, w) in &ent[off[k]..off[k + 1]] {
            dep[pi][c as usize] += layer.area_m2[k] * w;
        }
    }
    dep
}

/// The parcel step inside the advance (after the Eulerian claim loop):
/// move, consume collision overlaps, project into `out`. Books the
/// overrides in `ledger`, pushes the consumed volume to `collisions` (the
/// rate-limited-belt-v1 supply, `accrete`: sediment joins the crust).
#[allow(clippy::too_many_arguments)]
pub(crate) fn parcel_carry(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    layer: &mut ParcelLayer,
    params: &ParcelParams,
    poles: &BTreeMap<PlateId, EulerPoleV0>,
    dt_my: f64,
    oceanic_thickness_m: f64,
    accrete: bool,
    out: &mut CrustState,
    source: &mut [u32],
    collisions: &mut Vec<(usize, f64, f64, f64)>,
    ledger: &mut AdvanceLedger,
) -> Result<(ParcelProjection, ParcelStepReport), KinematicError> {
    let n = grid.cells();
    layer.reserve_identity_events(n, 0, layer.len())?;
    let mut rep = ParcelStepReport {
        parcels_before: layer.len() as u64,
        ..Default::default()
    };
    let (a0, v0, s0) = layer.totals();
    rep.area_before_m2 = a0;
    rep.volume_before_m3 = v0;
    rep.sediment_before_m3 = s0;
    // 1. Move (exact rigid rotation) and remember the step-start cell.
    let mut prev = Vec::with_capacity(layer.len());
    for k in 0..layer.len() {
        prev.push(nearest_cell(grid, layer.pos[k])? as u32);
        let pole = poles
            .get(&layer.plate[k])
            .ok_or(KinematicError::UnknownPlate { plate: layer.plate[k] })?;
        let q = rotate_by_euler(layer.pos[k], pole, dt_my)?;
        let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        layer.pos[k] = [q[0] / norm, q[1] / norm, q[2] / norm];
    }
    let kern = Kernel::new(grid, cents, areas)?;
    // 2. Collision consumption.
    {
        let (off, ent) = kern.all(layer)?;
        let pidx = plate_index(layer, &[]);
        let dep = deposits(layer, &off, &ent, &pidx, n);
        let np = pidx.len();
        let mut h_sum = vec![0.0f64; n];
        let mut d_sum = vec![0.0f64; n];
        for k in 0..layer.len() {
            for &(c, w) in &ent[off[k]..off[k + 1]] {
                let d = layer.area_m2[k] * w;
                h_sum[c as usize] += d * layer.thickness_m[k];
                d_sum[c as usize] += d;
            }
        }
        // Winner plate index per overlap cell.
        let mut overlap: Vec<Option<usize>> = vec![None; n];
        for c in 0..n {
            let tot: f64 = (0..np).map(|i| dep[i][c]).sum();
            if tot <= (1.0 + params.overlap_cell_fraction) * areas[c] {
                continue;
            }
            let mut big = 0;
            let mut win: Option<usize> = None;
            for i in 0..np {
                if dep[i][c] >= 0.5 * areas[c] {
                    big += 1;
                    if win.map_or(true, |w| dep[i][c] > dep[w][c]) {
                        win = Some(i);
                    }
                }
            }
            if big >= 2 {
                overlap[c] = win;
                rep.overlap_cells += 1;
            }
        }
        let mut keep = vec![true; layer.len()];
        for k in 0..layer.len() {
            let c = nearest_cell(grid, layer.pos[k])?;
            let Some(w) = overlap[c] else { continue };
            let pi = pidx[&layer.plate[k]];
            if pi == w || dep[pi][c] < 0.5 * areas[c] {
                continue;
            }
            if let Some(h) = layer.history.as_mut() {
                let id = layer.identity.as_ref().ok_or_else(||KinematicError::Grid("stack history requires carrier identity".into()))?.ids()[k];
                let stack = &layer.strata.as_ref().ok_or_else(||KinematicError::Grid("stack history requires strata".into()))?[k];
                h.retirement(id,out.t_my,layer.pos[k],layer.plate[k].0,layer.area_m2[k],stack).map_err(KinematicError::Grid)?;
            }
            if let Some(ids) = layer.identity.as_mut() {
                ids.retire(k, out.t_my, c, layer.plate[k].0).map_err(KinematicError::Grid)?;
            }
            keep[k] = false;
            let vol = layer.thickness_m[k] * layer.area_m2[k];
            let sed = layer.sediment_m[k] * layer.area_m2[k];
            let h_cell = if d_sum[c] > 0.0 { h_sum[c] / d_sum[c] } else { layer.thickness_m[k] };
            let mean_h = 0.5 * (layer.thickness_m[k] + h_cell);
            if accrete {
                collisions.push((c, vol + sed, 0.0, mean_h));
                ledger.belt_in_m3 += vol + sed;
                ledger.sediment_accreted_m3 += sed;
            } else {
                collisions.push((c, vol, sed, mean_h));
                ledger.belt_in_m3 += vol;
                ledger.belt_sediment_m3 += sed;
                ledger.sediment_stacked_m3 += sed;
            }
            ledger.collision_overlap_consumed_m3 += vol;
            if let Some(v) = layer.strata.as_ref() {
                if rep.consumed_strata_m3.is_empty() {
                    rep.consumed_strata_m3 = vec![0.0; crate::stratigraphy::NLITH];
                }
                v[k].lith_volumes(layer.area_m2[k], &mut rep.consumed_strata_m3);
            }
            rep.consumed_parcels += 1;
            rep.consumed_area_m2 += layer.area_m2[k];
            rep.consumed_m3 += vol;
            rep.consumed_sediment_m3 += sed;
        }
        if rep.consumed_parcels > 0 {
            layer.retain(&keep)?;
            let mut kp = 0;
            prev.retain(|_| {
                let x = keep[kp];
                kp += 1;
                x
            });
        }
    }
    // 3. Projection.
    let (off, ent) = kern.all(layer)?;
    let pidx = plate_index(layer, &[]);
    let ids: Vec<PlateId> = pidx.keys().cloned().collect();
    let dep = deposits(layer, &off, &ent, &pidx, n);
    let mut tot = vec![0.0f64; n];
    let mut acc = vec![[0.0f64; 5]; n];
    // Dominant parcel per cell: (deposit, parcel).
    let mut dom: Vec<Option<(f64, usize)>> = vec![None; n];
    for k in 0..layer.len() {
        for &(c, w) in &ent[off[k]..off[k + 1]] {
            let c = c as usize;
            let d = layer.area_m2[k] * w;
            tot[c] += d;
            acc[c][0] += d * layer.thickness_m[k];
            acc[c][1] += d * layer.sediment_m[k];
            acc[c][2] += d * layer.eroded_total_m[k];
            acc[c][3] += d * layer.arc_added_m[k];
            acc[c][4] += d * layer.shortening_added_m[k];
            if dom[c].map_or(true, |(bd, _)| d > bd) {
                dom[c] = Some((d, k));
            }
        }
    }
    let mut proj = ParcelProjection {
        continental: vec![false; n],
        thickness_m: vec![0.0; n],
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
    };
    for c in 0..n {
        proj.continental[c] = tot[c] >= 0.5 * areas[c] && tot[c] > 0.0;
    }
    let was_cont: Vec<bool> = out
        .kind
        .iter()
        .map(|k| *k == CrustKindCodeV0::Continental)
        .collect();
    let birth_src: Vec<bool> = source.iter().map(|&s| s == NO_SOURCE).collect();
    // Vacated cells first (they read only cells the projection leaves
    // untouched: oceanic in both the carry and the projection).
    let nlon = grid.nlon();
    for c in 0..n {
        if proj.continental[c] || !was_cont[c] {
            continue;
        }
        let (i, j) = (c / nlon, c % nlon);
        let mut pick: Option<(f64, usize)> = None;
        for (ti, tj, nb) in grid
            .neighbors(i, j)
            .map_err(|e| KinematicError::Grid(e.to_string()))?
        {
            let k = ti * nlon + tj;
            if proj.continental[k] || was_cont[k] || birth_src[k] || out.plate[k] != out.plate[c] {
                continue;
            }
            if pick.map_or(true, |(bd, bk)| nb.distance_m < bd || (nb.distance_m == bd && k < bk)) {
                pick = Some((nb.distance_m, k));
            }
        }
        out.material_pos[c] = cents[c];
        out.kind[c] = CrustKindCodeV0::Oceanic;
        out.sediment_m[c] = 0.0;
        out.eroded_total_m[c] = 0.0;
        match pick {
            Some((_, k)) => {
                out.birth_t_my[c] = out.birth_t_my[k];
                out.thickness_m[c] = out.thickness_m[k];
                out.arc_added_m[c] = out.arc_added_m[k];
                out.shortening_added_m[c] = out.shortening_added_m[k];
                source[c] = source[k];
                ledger.ridge_created_m3 += out.thickness_m[c] * areas[c];
                rep.vacated_copied += 1;
            }
            None => {
                out.birth_t_my[c] = out.t_my;
                out.thickness_m[c] = oceanic_thickness_m;
                out.arc_added_m[c] = 0.0;
                out.shortening_added_m[c] = 0.0;
                source[c] = NO_SOURCE;
                ledger.ridge_created_m3 += oceanic_thickness_m * areas[c];
                ledger.ridge_cells += 1;
                rep.vacated_births += 1;
            }
        }
    }
    for c in 0..n {
        if !proj.continental[c] {
            continue;
        }
        if !was_cont[c] {
            if birth_src[c] {
                ledger.ridge_created_m3 -= out.thickness_m[c] * areas[c];
                ledger.ridge_cells -= 1;
                rep.overridden_births += 1;
            } else {
                ledger.subducted_m3 += out.thickness_m[c] * areas[c];
                ledger.sediment_subducted_m3 += out.sediment_m[c] * areas[c];
                rep.overridden_oceanic_cells += 1;
            }
        }
        let t = tot[c];
        let (_, kd) = dom[c].expect("a continental cell has a depositor");
        // Plate with the largest deposit (ties: lowest id, ascending order).
        let mut best: Option<(f64, usize)> = None;
        for (i, d) in dep.iter().enumerate() {
            if d[c] > 0.0 && best.map_or(true, |(bd, _)| d[c] > bd) {
                best = Some((d[c], i));
            }
        }
        let plate = ids[best.expect("depositor").1];
        proj.thickness_m[c] = acc[c][0] / t;
        proj.sediment_m[c] = acc[c][1] / t;
        proj.eroded_total_m[c] = acc[c][2] / t;
        proj.arc_added_m[c] = acc[c][3] / t;
        proj.shortening_added_m[c] = acc[c][4] / t;
        out.material_pos[c] = cents[c];
        out.plate[c] = plate;
        out.kind[c] = CrustKindCodeV0::Continental;
        out.birth_t_my[c] = layer.birth_t_my[kd];
        out.thickness_m[c] = proj.thickness_m[c];
        out.sediment_m[c] = proj.sediment_m[c];
        out.eroded_total_m[c] = proj.eroded_total_m[c];
        out.arc_added_m[c] = proj.arc_added_m[c];
        out.shortening_added_m[c] = proj.shortening_added_m[c];
        source[c] = prev[kd];
        rep.projected_cells += 1;
        rep.projected_area_m2 += areas[c];
        rep.projected_volume_m3 += proj.thickness_m[c] * areas[c];
    }
    Ok((proj, rep))
}

/// Fold plate-id pieces (8-neighbour, longitude wraps) smaller than
/// `min_area_m2` into the plate they share the most neighbour links with
/// (ties: lowest id); up to 4 passes. Returns (cell, from, to) in
/// ascending cell order per pass.
pub fn fold_fragments(
    grid: &GridSpec,
    areas: &[f64],
    plate: &mut [PlateId],
    min_area_m2: f64,
) -> Result<(u64, Vec<(usize, PlateId, PlateId)>), KinematicError> {
    let n = grid.cells();
    let nlat = grid.nlat() as i64;
    let nlon = grid.nlon() as i64;
    let nbrs = |c: usize| -> Vec<usize> {
        let (i, j) = ((c as i64) / nlon, (c as i64) % nlon);
        let mut v = Vec::with_capacity(8);
        for di in -1..=1i64 {
            let ti = i + di;
            if ti < 0 || ti >= nlat {
                continue;
            }
            for dj in -1..=1i64 {
                if di == 0 && dj == 0 {
                    continue;
                }
                v.push((ti * nlon + (j + dj).rem_euclid(nlon)) as usize);
            }
        }
        v
    };
    let mut moved: Vec<(usize, PlateId, PlateId)> = Vec::new();
    let mut folded = 0u64;
    for _pass in 0..4 {
        let mut comp = vec![u32::MAX; n];
        let mut changes: Vec<(Vec<usize>, PlateId, PlateId)> = Vec::new();
        let mut ncomp = 0u32;
        for s in 0..n {
            if comp[s] != u32::MAX {
                continue;
            }
            let p = plate[s];
            let mut cells = vec![s];
            comp[s] = ncomp;
            let mut head = 0;
            while head < cells.len() {
                let c = cells[head];
                head += 1;
                for k in nbrs(c) {
                    if comp[k] == u32::MAX && plate[k] == p {
                        comp[k] = ncomp;
                        cells.push(k);
                    }
                }
            }
            let area: f64 = cells.iter().map(|&c| areas[c]).sum();
            if area < min_area_m2 {
                let mut links: BTreeMap<PlateId, u64> = BTreeMap::new();
                for &c in &cells {
                    for k in nbrs(c) {
                        if plate[k] != p {
                            *links.entry(plate[k]).or_insert(0) += 1;
                        }
                    }
                }
                let mut best: Option<(u64, PlateId)> = None;
                for (q, cnt) in links {
                    if best.map_or(true, |(bc, _)| cnt > bc) {
                        best = Some((cnt, q));
                    }
                }
                if let Some((_, q)) = best {
                    cells.sort_unstable();
                    changes.push((cells, p, q));
                }
            }
            ncomp += 1;
        }
        if changes.is_empty() {
            break;
        }
        for (cells, from, to) in changes {
            folded += 1;
            for c in cells {
                plate[c] = to;
                moved.push((c, from, to));
            }
        }
    }
    Ok((folded, moved))
}

/// Relabel the parcels whose nearest cell was folded (and whose plate is
/// the folded piece's). Returns the number relabelled.
pub(crate) fn relabel_parcels(
    grid: &GridSpec,
    layer: &mut ParcelLayer,
    moved: &[(usize, PlateId, PlateId)],
) -> Result<u64, KinematicError> {
    if moved.is_empty() {
        return Ok(0);
    }
    // Later passes override earlier ones for the same cell.
    let mut map: BTreeMap<usize, Vec<(PlateId, PlateId)>> = BTreeMap::new();
    for &(c, f, t) in moved {
        map.entry(c).or_default().push((f, t));
    }
    let mut count = 0;
    for k in 0..layer.len() {
        let c = nearest_cell(grid, layer.pos[k])?;
        if let Some(v) = map.get(&c) {
            let mut p = layer.plate[k];
            for &(f, t) in v {
                if p == f {
                    p = t;
                }
            }
            if p != layer.plate[k] {
                layer.plate[k] = p;
                count += 1;
            }
        }
    }
    Ok(count)
}

/// End-of-step write-back: every change the step made to a projected
/// continental cell goes to its parcels in proportion to their deposit;
/// cells that ended the step continental without being projected (rift
/// births) become new parcels.
pub fn write_back(
    grid: &GridSpec,
    areas: &[f64],
    layer: &mut ParcelLayer,
    proj: &ParcelProjection,
    state: &CrustState,
) -> Result<WriteBackReport, KinematicError> {
    write_back_impl(grid, areas, layer, proj, state, false)
}

/// SEED-RANGES (lane RG, 2026-10-13; part of the `cycle_ranges_v1`
/// laws, revision 1 in `unpublished-design-record`):
/// `write_back` with crust REMOVALS (a negative thickness change of a cell:
/// erosion, plateau cap, delamination) handed to the cell's parcels in
/// proportion to deposit x parcel thickness instead of deposit alone, so a
/// thin parcel mixed into a thick column is not driven below zero (clamped)
/// and left as a zero-thickness parcel that can later project a crustless
/// continental cell. Additions and every other field are as `write_back`.
pub fn write_back_thickness_weighted(
    grid: &GridSpec,
    areas: &[f64],
    layer: &mut ParcelLayer,
    proj: &ParcelProjection,
    state: &CrustState,
) -> Result<WriteBackReport, KinematicError> {
    write_back_impl(grid, areas, layer, proj, state, true)
}

fn write_back_impl(
    grid: &GridSpec,
    areas: &[f64],
    layer: &mut ParcelLayer,
    proj: &ParcelProjection,
    state: &CrustState,
    weighted: bool,
) -> Result<WriteBackReport, KinematicError> {
    Ok(write_back_core(grid, areas, layer, proj, state, weighted, 0)?.0)
}

/// HYPSO-LEVERS (2026-10-21, opt-in `write_back_residual_iterations`):
/// `write_back` (`weighted` = the thickness-weighted variant) followed by
/// `iterations` residual-correction passes (`residual_correction`).
pub fn write_back_corrected(
    grid: &GridSpec,
    areas: &[f64],
    layer: &mut ParcelLayer,
    proj: &ParcelProjection,
    state: &CrustState,
    weighted: bool,
    iterations: u32,
) -> Result<WriteBackReport, KinematicError> {
    Ok(write_back_core(grid, areas, layer, proj, state, weighted, iterations)?.0)
}

/// The write-back's kernel weights (CSR over the parcels present at entry),
/// the live-cell deposit and the live mask, kept for the stratigraphy pass.
struct WriteBackWeights {
    off: Vec<usize>,
    ent: Vec<(u32, f64)>,
    dep: Vec<f64>,
    live: Vec<bool>,
}

fn write_back_core(
    grid: &GridSpec,
    areas: &[f64],
    layer: &mut ParcelLayer,
    proj: &ParcelProjection,
    state: &CrustState,
    weighted: bool,
    residual_iterations: u32,
) -> Result<(WriteBackReport, WriteBackWeights), KinematicError> {
    let n = grid.cells();
    layer.reserve_identity_events(n, n, 0)?;
    let cents = crate::kinematic::centroids(grid)?;
    let kern = Kernel::new(grid, &cents, areas)?;
    let (off, ent) = kern.all(layer)?;
    let mut rep = WriteBackReport::default();
    let live: Vec<bool> = (0..n)
        .map(|c| proj.continental[c] && state.kind[c] == CrustKindCodeV0::Continental)
        .collect();
    for c in 0..n {
        if proj.continental[c] && !live[c] {
            rep.lost_cells += 1;
        }
    }
    let mut dep = vec![0.0f64; n];
    for k in 0..layer.len() {
        for &(c, w) in &ent[off[k]..off[k + 1]] {
            if live[c as usize] {
                dep[c as usize] += layer.area_m2[k] * w;
            }
        }
    }
    // SEED-RANGES: thickness-weighted deposit for removals.
    let mut dep_h = if weighted { vec![0.0f64; n] } else { Vec::new() };
    if weighted {
        for k in 0..layer.len() {
            for &(c, w) in &ent[off[k]..off[k + 1]] {
                if live[c as usize] {
                    dep_h[c as usize] += layer.area_m2[k] * w * layer.thickness_m[k].max(0.0);
                }
            }
        }
    }
    let mut g_neg = if weighted { vec![0.0f64; n] } else { Vec::new() };
    // Per-cell change times area over the deposit (m per unit weight).
    let mut g = vec![[0.0f64; 5]; n];
    for c in 0..n {
        if !live[c] || !(dep[c] > 0.0) {
            continue;
        }
        let d = [
            state.thickness_m[c] - proj.thickness_m[c],
            state.sediment_m[c] - proj.sediment_m[c],
            state.eroded_total_m[c] - proj.eroded_total_m[c],
            state.arc_added_m[c] - proj.arc_added_m[c],
            state.shortening_added_m[c] - proj.shortening_added_m[c],
        ];
        rep.writeback_m3 += d[0] * areas[c];
        rep.writeback_sediment_m3 += d[1] * areas[c];
        let f = areas[c] / dep[c];
        for x in 0..5 {
            g[c][x] = d[x] * f;
        }
        if weighted && d[0] < 0.0 && dep_h[c] > 0.0 {
            g[c][0] = 0.0;
            g_neg[c] = d[0] * areas[c] / dep_h[c];
        }
    }
    for k in 0..layer.len() {
        let mut dk = [0.0f64; 5];
        let hk = layer.thickness_m[k].max(0.0);
        for &(c, w) in &ent[off[k]..off[k + 1]] {
            let c = c as usize;
            if !live[c] {
                continue;
            }
            for x in 0..5 {
                dk[x] += g[c][x] * w;
            }
            if weighted {
                dk[0] += g_neg[c] * w * hk;
            }
        }
        layer.thickness_m[k] += dk[0];
        layer.sediment_m[k] += dk[1];
        layer.eroded_total_m[k] += dk[2];
        layer.arc_added_m[k] += dk[3];
        layer.shortening_added_m[k] += dk[4];
        if layer.thickness_m[k] < 0.0 {
            rep.clamped_m3 += -layer.thickness_m[k] * layer.area_m2[k];
            layer.thickness_m[k] = 0.0;
        }
        if layer.sediment_m[k] < 0.0 {
            rep.clamped_sediment_m3 += -layer.sediment_m[k] * layer.area_m2[k];
            layer.sediment_m[k] = 0.0;
        }
    }
    residual_correction(areas, layer, state, &off, &ent, &dep, &live, residual_iterations, &mut rep);
    for c in 0..n {
        if state.kind[c] == CrustKindCodeV0::Continental && !proj.continental[c] {
            layer.push_cell(state, &cents, areas, c)?;
            rep.created_parcels += 1;
            rep.created_area_m2 += areas[c];
            rep.created_m3 += state.thickness_m[c] * areas[c];
            rep.created_sediment_m3 += state.sediment_m[c] * areas[c];
        }
    }
    let (a, v, s) = layer.totals();
    rep.parcels_after = layer.len() as u64;
    rep.area_after_m2 = a;
    rep.volume_after_m3 = v;
    rep.sediment_after_m3 = s;
    Ok((rep, WriteBackWeights { off, ent, dep, live }))
}

/// HYPSO-LEVERS (2026-10-21) write-back residual correction (no-op at 0
/// passes). Each pass re-projects the parcels present at entry with the
/// write-back's kernel weights (the projection's deposit-weighted mean),
/// takes the residual r = end-of-step thickness - re-projection on live
/// cells, removes its area-weighted mean (so sum(parcel area x change) is
/// 0 and the parcel volume ledger is unchanged) and hands it to the
/// parcels: a positive residual in proportion to deposit (as `write_back`
/// hands a change), a negative one in proportion to deposit x parcel
/// thickness (as `write_back_thickness_weighted` hands a removal), so a thin
/// parcel mixed into a thick column is not driven to zero (with the plain
/// deposit rule for both signs, seeds 240907 and 7 projected a crustless
/// continental cell at 14 and 356 My). A thickness driven below zero is
/// clamped and booked in `clamped_m3`, as there.
/// Why: projection then write-back is not an identity; each step's changes
/// return smoothed over about two kernel radii, which lowers thick belts
/// (lane HYPSO-LEVERS measured -0.26 Mkm3/My of crust excess above 36 km,
/// seed 42, 40-400 My; 3 passes left -0.05).
#[allow(clippy::too_many_arguments)]
fn residual_correction(
    areas: &[f64],
    layer: &mut ParcelLayer,
    state: &CrustState,
    off: &[usize],
    ent: &[(u32, f64)],
    dep: &[f64],
    live: &[bool],
    iterations: u32,
    rep: &mut WriteBackReport,
) {
    let n = live.len();
    for _ in 0..iterations {
        let k_in = off.len() - 1;
        let mut num = vec![0.0f64; n];
        let mut tot = vec![0.0f64; n];
        for k in 0..k_in {
            for &(c, w) in &ent[off[k]..off[k + 1]] {
                let d = layer.area_m2[k] * w;
                tot[c as usize] += d;
                num[c as usize] += d * layer.thickness_m[k];
            }
        }
        let mut r = vec![0.0f64; n];
        let (mut sr, mut sa) = (0.0f64, 0.0f64);
        for c in 0..n {
            if live[c] && dep[c] > 0.0 && tot[c] > 0.0 {
                r[c] = state.thickness_m[c] - num[c] / tot[c];
                sr += r[c] * areas[c];
                sa += areas[c];
            }
        }
        let mean = if sa > 0.0 { sr / sa } else { 0.0 };
        // Thickness-weighted deposit for the negative residuals.
        let mut dep_h = vec![0.0f64; n];
        for k in 0..k_in {
            let hk = layer.thickness_m[k].max(0.0);
            for &(c, w) in &ent[off[k]..off[k + 1]] {
                dep_h[c as usize] += layer.area_m2[k] * w * hk;
            }
        }
        let mut g = vec![0.0f64; n];
        let mut g_neg = vec![0.0f64; n];
        for c in 0..n {
            if live[c] && dep[c] > 0.0 && tot[c] > 0.0 {
                let d = r[c] - mean;
                if d < 0.0 && dep_h[c] > 0.0 {
                    g_neg[c] = d * areas[c] / dep_h[c];
                } else {
                    g[c] = d * areas[c] / dep[c];
                }
            }
        }
        for k in 0..k_in {
            let mut dk = 0.0f64;
            let hk = layer.thickness_m[k].max(0.0);
            for &(c, w) in &ent[off[k]..off[k + 1]] {
                dk += g[c as usize] * w + g_neg[c as usize] * w * hk;
            }
            layer.thickness_m[k] += dk;
            if layer.thickness_m[k] < 0.0 {
                rep.clamped_m3 += -layer.thickness_m[k] * layer.area_m2[k];
                layer.thickness_m[k] = 0.0;
            }
        }
    }
}

/// SEED-STRATIGRAPHY stage A (lane ST, 2026-10-13) `parcel-stacks-v1`:
/// `write_back` (or the thickness-weighted variant), then the stratigraphy
/// pass on the same kernel weights: every parcel present at entry gets its
/// share of the tagged per-cell deltas (`tags`, plus the arc and
/// shortening lineage changes projection -> step end), applied by
/// `stratigraphy::apply_parcel`; births get `stratigraphy::birth_stack`.
/// The parcel fields are written exactly as `write_back` writes them.
#[allow(clippy::too_many_arguments)]
pub fn write_back_strata(
    grid: &GridSpec,
    areas: &[f64],
    layer: &mut ParcelLayer,
    proj: &ParcelProjection,
    state: &CrustState,
    weighted: bool,
    tags: &crate::stratigraphy::CellTags,
    sp: &crate::stratigraphy::StrataParams,
    t0_my: f64,
) -> Result<(WriteBackReport, crate::stratigraphy::StrataStepLog), KinematicError> {
    write_back_strata_corrected(grid, areas, layer, proj, state, weighted, 0, tags, sp, t0_my)
}

/// HYPSO-LEVERS (2026-10-21, opt-in `write_back_residual_iterations`):
/// `write_back_strata` with the residual correction (`residual_correction`)
/// after the write-back and before births; the stacks then take the
/// corrected parcel column (the correction is basement flow, step 9 of
/// `stratigraphy::apply_parcel`). 0 passes = `write_back_strata`.
#[allow(clippy::too_many_arguments)]
pub fn write_back_strata_corrected(
    grid: &GridSpec,
    areas: &[f64],
    layer: &mut ParcelLayer,
    proj: &ParcelProjection,
    state: &CrustState,
    weighted: bool,
    residual_iterations: u32,
    tags: &crate::stratigraphy::CellTags,
    sp: &crate::stratigraphy::StrataParams,
    t0_my: f64,
) -> Result<(WriteBackReport, crate::stratigraphy::StrataStepLog), KinematicError> {
    use crate::stratigraphy::{apply_parcel, birth_stack, ParcelAmounts, StrataStepLog};
    let n = grid.cells();
    let k0 = layer.len();
    let mut log = StrataStepLog::default();
    match layer.strata.as_ref() {
        Some(v) if v.len() == k0 => {}
        _ => {
            return Err(KinematicError::Grid(format!(
                "parcel-stacks-v1: {} stacks for {k0} parcels",
                layer.strata.as_ref().map_or(0, |v| v.len())
            )))
        }
    }
    {
        let v = layer.strata.as_ref().expect("checked");
        for k in 0..k0 {
            let col = layer.thickness_m[k] + layer.sediment_m[k];
            let tot = v[k].total_m();
            log.before_m3 += layer.area_m2[k] * tot;
            log.entry_mismatch_max_rel = log.entry_mismatch_max_rel.max((tot - col).abs() / col.abs().max(1.0));
        }
    }
    let (rep, w) = write_back_core(grid, areas, layer, proj, state, weighted, residual_iterations)?;
    // parcel-stacks-v2 (lane SA): positive shortening is the belt's placed
    // supply only; the rest of the lineage change is basement only.
    let v2 = sp.v2.is_some();
    if v2 && (tags.belt_placed_m.len() != n || tags.tect_elev_m.len() != n || tags.post_elev_m.len() != n || tags.warm_shelf.len() != n) {
        return Err(KinematicError::Grid("parcel-stacks-v2: tags without belt / elevation / shelf fields".into()));
    }
    // Per-cell amounts per unit deposit weight (m): g = change x area / dep.
    const NT: usize = 9;
    let mut g = vec![[0.0f64; NT]; n];
    for c in 0..n {
        if !w.live[c] || !(w.dep[c] > 0.0) {
            continue;
        }
        let f = areas[c] / w.dep[c];
        let ds = state.shortening_added_m[c] - proj.shortening_added_m[c];
        let (sp_pos, sp_neg, rest) = if v2 {
            // PREREG revision 1 (lane SA): the shortening applied to the
            // stack (pluton share + pure shear) is the cell's net lineage
            // change clamped to [0, placed]; what is left (inflow beyond the
            // placed supply, or outflow) is basement only.
            let placed = tags.belt_placed_m[c].max(0.0);
            let eff = ds.clamp(0.0, placed);
            (eff, 0.0, ds - eff)
        } else {
            (ds.max(0.0), (-ds).max(0.0), 0.0)
        };
        g[c] = [
            tags.erosion_m[c] * f,
            tags.deposit_m[c] * f,
            (state.arc_added_m[c] - proj.arc_added_m[c]).max(0.0) * f,
            tags.plume_m[c] * f,
            tags.thinning_m[c] * f,
            tags.basement_add_m[c] * f,
            sp_pos * f,
            sp_neg * f,
            rest * f,
        ];
    }
    let mut stacks = layer.strata.take().expect("checked");
    for k in 0..k0 {
        let selected_id = layer.history.as_ref().and_then(|h| {
            let id = layer.identity.as_ref()?.ids()[k];
            h.selected(id).then_some(id)
        });
        let before = selected_id.map(|_|stacks[k].clone());
        let mut support=selected_id.map(|_|Vec::with_capacity(w.off[k+1]-w.off[k]));
        let mut reference = None;
        let mut am = [0.0f64; NT];
        let (mut best_dep, mut class, mut warm_shelf) = (0.0f64, 0u8, false);
        let (mut best_w, mut land) = (-1.0f64, false);
        // v2: kernel-weighted elevations of the live cells (rock uplift).
        let (mut ew, mut et, mut ep) = (0.0f64, 0.0f64, 0.0f64);
        for &(c, wk) in &w.ent[w.off[k]..w.off[k + 1]] {
            let c = c as usize;
            if let Some(rows)=support.as_mut() {
                rows.push(crate::carrier_history::CellSupport {
                    cell:c as u32,weight:wk,live:w.live[c],area_m2:areas[c],deposited_area_m2:w.dep[c],
                    raw:[tags.erosion_m[c],tags.deposit_m[c],state.arc_added_m[c]-proj.arc_added_m[c],
                        tags.plume_m[c],tags.thinning_m[c],tags.basement_add_m[c],
                        state.shortening_added_m[c]-proj.shortening_added_m[c],if v2 {tags.belt_placed_m[c]}else{0.}],
                    per_unit:g[c],class:tags.class[c],land:tags.land[c],warm_shelf:v2&&tags.warm_shelf[c],
                    tect_elev_m:if v2 {tags.tect_elev_m[c]}else{0.},post_elev_m:if v2 {tags.post_elev_m[c]}else{0.},
                });
            }
            if wk > best_w {
                best_w = wk;
                land = tags.land[c];
            }
            if !w.live[c] {
                continue;
            }
            for x in 0..NT {
                am[x] += g[c][x] * wk;
            }
            let d = g[c][1] * wk;
            if d > best_dep {
                best_dep = d;
                class = tags.class[c];
                warm_shelf = v2 && tags.warm_shelf[c];
            }
            if v2 {
                ew += wk;
                et += wk * tags.tect_elev_m[c];
                ep += wk * tags.post_elev_m[c];
            }
        }
        if v2 && ew > 0.0 {
            reference = Some([et / ew, ep / ew]);
            let st = &mut stacks[k];
            st.rock_uplift_m += et / ew - st.elev_ref_m;
            st.elev_ref_m = ep / ew;
        }
        let a = ParcelAmounts {
            erosion: am[0],
            deposit: am[1],
            class,
            arc: am[2],
            plume: am[3],
            thinning: am[4],
            basement_add: am[5],
            short_pos: am[6],
            short_neg: am[7],
            land,
            belt_rest: am[8],
            warm_shelf,
        };
        let total = layer.thickness_m[k] + layer.sediment_m[k];
        let log_before = selected_id.map(|_|log.clone());
        let evolution=if let Some(provenance)=layer.layer_provenance.as_mut() {
            let retain=selected_id.is_some()&&layer.history.as_ref().is_some_and(|h|h.config().layer_evolution);
            let id=layer.identity.as_ref().expect("provenance requires identities").ids()[k];
            let (p,trace)=crate::stratigraphy::apply_parcel_provenance(&mut stacks[k],&a,total,layer.area_m2[k],sp,t0_my,state.t_my,&mut log,std::mem::take(&mut provenance[k]),id,retain).map_err(KinematicError::Grid)?;
            provenance[k]=p;retain.then_some(trace)
        }else if selected_id.is_some()&&layer.history.as_ref().is_some_and(|h|h.config().layer_evolution) {
            Some(crate::stratigraphy::apply_parcel_recorded(&mut stacks[k],&a,total,layer.area_m2[k],sp,t0_my,state.t_my,&mut log,crate::layer_evolution::byte_bound(sp.max_layers)).map_err(KinematicError::Grid)?)
        }else{apply_parcel(&mut stacks[k], &a, total, layer.area_m2[k], sp, t0_my, state.t_my, &mut log);None};
        if let Some(id)=selected_id {
            layer.history.as_mut().expect("selected history").operation(id,t0_my,state.t_my,
                layer.pos[k],layer.plate[k].0,layer.area_m2[k],before.expect("captured"),&stacks[k],a,
                reference,total,log_before.expect("captured"),&log,support.expect("captured"),evolution).map_err(KinematicError::Grid)?;
        }
    }
    // Births (pushed by the write-back in ascending cell order).
    let mut k = k0;
    for c in 0..n {
        if state.kind[c] == CrustKindCodeV0::Continental && !proj.continental[c] {
            let selected_id = layer.history.as_ref().and_then(|h| {
                let id = layer.identity.as_ref()?.ids()[k];
                h.selected(id).then_some(id)
            });
            let log_before = selected_id.map(|_|log.clone());
            let (born,mut evolution)=if layer.layer_provenance.is_some()||(selected_id.is_some()&&layer.history.as_ref().is_some_and(|h|h.config().layer_evolution)) {
                let (stack,trace)=crate::stratigraphy::birth_stack_recorded(state.thickness_m[c],state.sediment_m[c],state.arc_added_m[c],tags.class[c],sp,state.t_my,layer.area_m2[k],&mut log,crate::layer_evolution::byte_bound(sp.max_layers)).map_err(KinematicError::Grid)?;(stack,Some(trace))
            }else{(birth_stack(state.thickness_m[c],state.sediment_m[c],state.arc_added_m[c],tags.class[c],sp,state.t_my,layer.area_m2[k],&mut log),None)};
            if let Some(provenance)=layer.layer_provenance.as_mut(){
                let id=layer.identity.as_ref().expect("provenance requires identities").ids()[k];
                provenance[k]=crate::stratigraphy::birth_provenance(&born,evolution.as_ref().expect("recorded birth"),state.arc_added_m[c].max(0.).min(state.thickness_m[c].max(0.)),tags.class[c],sp,state.t_my,id);
                if !layer.history.as_ref().is_some_and(|h|h.config().layer_evolution){evolution=None;}
            }
            stacks[k]=born;
            if v2 {
                stacks[k].elev_ref_m = tags.post_elev_m[c];
            }
            if let Some(id)=selected_id {
                layer.history.as_mut().expect("selected history").birth(id,state.t_my,layer.pos[k],
                    layer.plate[k].0,layer.area_m2[k],[state.thickness_m[c],state.sediment_m[c],state.arc_added_m[c]],
                    tags.class[c],v2.then_some([0.,tags.post_elev_m[c]]),&stacks[k],log_before.expect("captured"),&log,evolution)
                    .map_err(KinematicError::Grid)?;
            }
            k += 1;
        }
    }
    let mut layers = 0u64;
    let (mut ua, mut uu) = (0.0f64, 0.0f64);
    for k in 0..layer.len() {
        let st = &stacks[k];
        if v2 {
            ua += layer.area_m2[k];
            uu += layer.area_m2[k] * st.rock_uplift_m.abs();
        }
        log.after_m3 += layer.area_m2[k] * st.total_m();
        log.column_m3 += layer.area_m2[k] * (layer.thickness_m[k] + layer.sediment_m[k]);
        layers += st.cover.len() as u64;
        if st.cover.len() + 1 >= sp.max_layers {
            log.stacks_at_cap += 1;
        }
        if st.basement_m < sp.basement_min_m {
            log.basement_below_min += 1;
        }
    }
    log.parcels = layer.len() as u64;
    log.stacks = stacks.len() as u64;
    layer.strata = Some(stacks);
    log.mean_cover_layers = if layer.is_empty() { 0.0 } else { layers as f64 / layer.len() as f64 };
    log.mean_abs_rock_uplift_m = if ua > 0.0 { uu / ua } else { 0.0 };
    let births: f64 = log.births_m3.iter().sum();
    log.closure_residual_m3 = log.after_m3 - (log.before_m3 + log.ops_m3() + births);
    log.closure_rel = log.closure_residual_m3.abs() / log.after_m3.abs().max(1.0);
    log.column_mismatch_rel = (log.after_m3 - log.column_m3).abs() / log.column_m3.abs().max(1.0);
    Ok((rep, log))
}

/// SEED-STRATIGRAPHY: the dominant parcel (largest deposit; ties: lowest
/// parcel) of every continental cell of `state` (`None` elsewhere).
pub fn dominant_parcels(grid: &GridSpec, layer: &ParcelLayer, state: &CrustState) -> Result<Vec<Option<usize>>, KinematicError> {
    let n = grid.cells();
    let cents = crate::kinematic::centroids(grid)?;
    let areas = crate::kinematic::cell_areas(grid)?;
    let kern = Kernel::new(grid, &cents, &areas)?;
    let (off, ent) = kern.all(layer)?;
    let mut dom: Vec<Option<(f64, usize)>> = vec![None; n];
    for k in 0..layer.len() {
        for &(c, w) in &ent[off[k]..off[k + 1]] {
            let d = layer.area_m2[k] * w;
            let c = c as usize;
            if dom[c].map_or(true, |(bd, _)| d > bd) {
                dom[c] = Some((d, k));
            }
        }
    }
    Ok((0..n)
        .map(|c| if state.kind[c] == CrustKindCodeV0::Continental { dom[c].map(|x| x.1) } else { None })
        .collect())
}

/// SEED-STRATIGRAPHY stage A.1 (lane SA, `parcel-stacks-v2`): the per-cell
/// cumulative rock uplift (kernel-weighted mean of the parcels' uplift over
/// each continental cell; 0 on oceanic cells) and its gradient east / north
/// (m per m): central differences over continental neighbours on the
/// sphere, one-sided where only one side is continental, 0 when isolated;
/// the polar rows use their single interior neighbour north-south.
pub fn rock_uplift_field(grid: &GridSpec, layer: &ParcelLayer, state: &CrustState) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), KinematicError> {
    let n = grid.cells();
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    let mut u = vec![0.0f64; n];
    let mut gx = vec![0.0f64; n];
    let mut gy = vec![0.0f64; n];
    let Some(stacks) = layer.strata.as_ref() else {
        return Ok((u, gx, gy));
    };
    let cents = crate::kinematic::centroids(grid)?;
    let areas = crate::kinematic::cell_areas(grid)?;
    let kern = Kernel::new(grid, &cents, &areas)?;
    let (off, ent) = kern.all(layer)?;
    let mut wsum = vec![0.0f64; n];
    for k in 0..layer.len() {
        for &(c, w) in &ent[off[k]..off[k + 1]] {
            let c = c as usize;
            let d = layer.area_m2[k] * w;
            wsum[c] += d;
            u[c] += d * stacks[k].rock_uplift_m;
        }
    }
    let cont: Vec<bool> = (0..n).map(|c| state.kind[c] == CrustKindCodeV0::Continental).collect();
    for c in 0..n {
        u[c] = if cont[c] && wsum[c] > 0.0 { u[c] / wsum[c] } else { 0.0 };
    }
    let r = grid.radius_m();
    let dlon = std::f64::consts::TAU / nlon as f64;
    let dlat = std::f64::consts::PI / nlat as f64;
    for i in 0..nlat {
        let lat = grid.center_lat_deg(i).map_err(|e| KinematicError::Grid(e.to_string()))?.to_radians();
        let dx = r * lat.cos().max(1e-6) * dlon;
        let dy = r * dlat;
        for j in 0..nlon {
            let c = i * nlon + j;
            if !cont[c] {
                continue;
            }
            let e = i * nlon + (j + 1) % nlon;
            let w = i * nlon + (j + nlon - 1) % nlon;
            gx[c] = match (cont[e], cont[w]) {
                (true, true) => (u[e] - u[w]) / (2.0 * dx),
                (true, false) => (u[e] - u[c]) / dx,
                (false, true) => (u[c] - u[w]) / dx,
                _ => 0.0,
            };
            let nn = if i + 1 < nlat { Some((i + 1) * nlon + j) } else { None }.filter(|&k| cont[k]);
            let ss = if i > 0 { Some((i - 1) * nlon + j) } else { None }.filter(|&k| cont[k]);
            gy[c] = match (nn, ss) {
                (Some(a), Some(b)) => (u[a] - u[b]) / (2.0 * dy),
                (Some(a), None) => (u[a] - u[c]) / dy,
                (None, Some(b)) => (u[c] - u[b]) / dy,
                _ => 0.0,
            };
        }
    }
    Ok((u, gx, gy))
}

/// Relabel the parcels of `parent` on the positive side of `cut_normal`
/// to `child` (the `split_plate` rule, applied to parcel positions).
pub fn split_parcels(layer: &mut ParcelLayer, parent: PlateId, child: PlateId, cut_normal: [f64; 3]) -> u64 {
    let mut count = 0;
    for k in 0..layer.len() {
        let x = layer.pos[k];
        if layer.plate[k] == parent
            && x[0] * cut_normal[0] + x[1] * cut_normal[1] + x[2] * cut_normal[2] > 0.0
        {
            layer.plate[k] = child;
            count += 1;
        }
    }
    count
}

/// Union-find state of `contact-length-pole-share-v1` (plate -> parent).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SutureState {
    pub parent: BTreeMap<PlateId, PlateId>,
}

impl SutureState {
    pub fn find(&self, p: PlateId) -> PlateId {
        let mut r = p;
        while let Some(&q) = self.parent.get(&r) {
            if q == r {
                break;
            }
            r = q;
        }
        r
    }
}

/// One suture: the pair, its contact length and the surviving pole's root.
#[derive(Debug, Clone, PartialEq)]
pub struct SutureEvent {
    pub t_my: f64,
    pub a: PlateId,
    pub b: PlateId,
    pub contact_m: f64,
    pub winner_root: PlateId,
    pub loser_root: PlateId,
}

/// `contact-length-pole-share-v1`: measure the converging continent-
/// continent contact per plate pair on `state`, suture pairs at or above
/// `threshold_m` (smaller group takes the larger group's root pole), then
/// give every plate its group root's pole. Returns the new sutures and the
/// largest contact seen (m).
pub fn suture_pass(
    grid: &GridSpec,
    state: &CrustState,
    plates: &mut [PlateRecordV0],
    sut: &mut SutureState,
    threshold_m: f64,
) -> Result<(Vec<SutureEvent>, f64), KinematicError> {
    let cents = crate::kinematic::centroids(grid)?;
    let areas = crate::kinematic::cell_areas(grid)?;
    let poles: BTreeMap<PlateId, EulerPoleV0> = plates
        .iter()
        .filter_map(|p| p.euler_pole.map(|e| (p.plate_id, e)))
        .collect();
    let radius = grid.radius_m();
    let dlat = std::f64::consts::PI / grid.nlat() as f64;
    let dlon = 2.0 * std::f64::consts::PI / grid.nlon() as f64;
    let nlon = grid.nlon();
    let mut contact: BTreeMap<(PlateId, PlateId), f64> = BTreeMap::new();
    for c in 0..state.cells() {
        if state.kind[c] != CrustKindCodeV0::Continental {
            continue;
        }
        let (i, j) = (c / nlon, c % nlon);
        let lat_bot = -std::f64::consts::FRAC_PI_2 + i as f64 * dlat;
        for (ti, tj, nb) in grid
            .neighbors(i, j)
            .map_err(|e| KinematicError::Grid(e.to_string()))?
        {
            if nb.di.abs() + nb.dj.abs() != 1 {
                continue;
            }
            let k = ti * nlon + tj;
            let (pa_id, pb_id) = (state.plate[c], state.plate[k]);
            if k <= c || pa_id == pb_id || state.kind[k] != CrustKindCodeV0::Continental {
                continue;
            }
            let (Some(pa), Some(pb)) = (poles.get(&pa_id), poles.get(&pb_id)) else {
                return Err(KinematicError::UnknownPlate { plate: pa_id });
            };
            let raw = [
                cents[c][0] + cents[k][0],
                cents[c][1] + cents[k][1],
                cents[c][2] + cents[k][2],
            ];
            let norm = (raw.iter().map(|x| x * x).sum::<f64>()).sqrt();
            if norm < 1e-12 {
                continue;
            }
            let mid = [raw[0] / norm, raw[1] / norm, raw[2] / norm];
            let vn = normal_rate_at(mid, cents[c], cents[k], pa, pb)? * radius;
            if !(vn > 0.0) {
                continue;
            }
            let edge = if nb.di == 0 {
                radius * dlat
            } else if nb.di > 0 {
                radius * (lat_bot + dlat).cos() * dlon
            } else {
                radius * lat_bot.cos() * dlon
            };
            let key = if pa_id < pb_id { (pa_id, pb_id) } else { (pb_id, pa_id) };
            *contact.entry(key).or_insert(0.0) += edge;
        }
    }
    let max_contact = contact.values().cloned().fold(0.0, f64::max);
    let mut pairs: Vec<((PlateId, PlateId), f64)> = contact
        .into_iter()
        .filter(|(_, l)| *l >= threshold_m)
        .collect();
    pairs.sort_by(|x, y| y.1.total_cmp(&x.1).then(x.0.cmp(&y.0)));
    let mut plate_area: BTreeMap<PlateId, f64> = BTreeMap::new();
    for c in 0..state.cells() {
        *plate_area.entry(state.plate[c]).or_insert(0.0) += areas[c];
    }
    let mut events = Vec::new();
    for ((a, b), l) in pairs {
        let (ra, rb) = (sut.find(a), sut.find(b));
        if ra == rb {
            continue;
        }
        let group_area = |r: PlateId, s: &SutureState| -> f64 {
            plate_area
                .iter()
                .filter(|(p, _)| s.find(**p) == r)
                .map(|(_, x)| *x)
                .sum()
        };
        let (aa, ab) = (group_area(ra, sut), group_area(rb, sut));
        let (w, lo) = if aa > ab || (aa == ab && ra < rb) { (ra, rb) } else { (rb, ra) };
        sut.parent.insert(lo, w);
        events.push(SutureEvent {
            t_my: state.t_my,
            a,
            b,
            contact_m: l,
            winner_root: w,
            loser_root: lo,
        });
    }
    // Every plate moves with its group root's pole.
    let root_pole: BTreeMap<PlateId, Option<EulerPoleV0>> =
        plates.iter().map(|p| (p.plate_id, p.euler_pole)).collect();
    for p in plates.iter_mut() {
        let r = sut.find(p.plate_id);
        if r != p.plate_id {
            if let Some(Some(pole)) = root_pole.get(&r) {
                p.euler_pole = Some(*pole);
            }
        }
    }
    Ok((events, max_contact))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rotation::lat_lon_to_unit;

    #[test]
    fn kernel_weights_are_normalised_and_cover_the_containing_cell() {
        let grid = GridSpec::new(36, 72, 6.371e6).unwrap();
        let cents = crate::kinematic::centroids(&grid).unwrap();
        let areas = crate::kinematic::cell_areas(&grid).unwrap();
        let k = Kernel::new(&grid, &cents, &areas).unwrap();
        for (lat, lon) in [(0.0, 0.0), (2.4, 3.3), (60.1, -170.0), (88.0, 10.0), (-89.9, 0.0)] {
            let p = lat_lon_to_unit(lat, lon).unwrap();
            let mut w = Vec::new();
            k.weights(p, &mut w).unwrap();
            let s: f64 = w.iter().map(|x| x.1).sum();
            assert!((s - 1.0).abs() < 1e-12, "{lat} {lon} {s}");
            let c0 = nearest_cell(&grid, p).unwrap() as u32;
            assert!(w.iter().any(|x| x.0 == c0), "{lat} {lon}");
            // Strictly ascending cells (no duplicates).
            assert!(w.windows(2).all(|x| x[0].0 < x[1].0));
            assert!(w.len()<=3*grid.nlon());
            if lat.abs()>89. {assert!(w.len()>8,"polar support is not bounded by an eight-entry allocation hint");}
        }
    }

    #[test]
    fn fold_relabels_small_pieces_only() {
        let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
        let areas = crate::kinematic::cell_areas(&grid).unwrap();
        let n = grid.cells();
        let mut plate = vec![PlateId(0); n];
        // Big plate 1: rows 0..9; a single speck of plate 2 inside plate 0.
        for c in 0..9 * 36 {
            plate[c] = PlateId(1);
        }
        plate[12 * 36 + 5] = PlateId(2);
        let total: f64 = areas.iter().sum();
        // 18 x 36 cells are ~0.15 % of the sphere each: fold below 1 %.
        let (folded, moved) = fold_fragments(&grid, &areas, &mut plate, 1e-2 * total).unwrap();
        assert_eq!(folded, 1);
        assert_eq!(moved, vec![(12 * 36 + 5, PlateId(2), PlateId(0))]);
        assert_eq!(plate.iter().filter(|p| **p == PlateId(1)).count(), 9 * 36);
    }
}
