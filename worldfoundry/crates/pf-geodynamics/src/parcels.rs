//! The Lagrangian crust-parcel ledger (brief K2/K4; spec §55 material
//! lineage — crust parcels are pf-geology's FIRST material-lineage
//! instance). THE authoritative material identity/history; see
//! `conservation.rs` for the two gates.
//!
//! PARCEL != CELL IS STRUCTURAL HERE (overnight agent-B acceptance):
//! - `ParcelLedger` carries NO spatial sampling of any kind. It knows what
//!   rock exists, when it was born, how it is related and who owns it —
//!   never where it currently sits in a grid.
//! - Where material currently sits is the DERIVED weighted projection:
//!   `project_onto_cells` consumes a driver-supplied position map (parcel
//!   -> element) and emits `CellParcelContribution`s at weight 1.0 for B1
//!   nearest-cell sampling. Two parcels may land in one cell and NEITHER
//!   disappears; a future conservative/PIC remapper may split one parcel
//!   across cells — the ledger neither knows nor cares.
//! - Conservation NEVER consults the projection (that ordering is the
//!   whole point: identity first, sampling second).
//!
//! Age is ACCRUED: `birth_t_my` is fixed at creation and age is
//! `t_now - birth` — never recomputed from current ridge distance (K2 rev
//! 2: that rewrite-what-happened failure is what this ledger exists to
//! prevent).

use crate::plates::PlateSeeding;
use pf_contracts::domain::{
    CellParcelContribution, CrustParcelRecordV0, EventId, ParcelBirthV0, ParcelId, ParcelLedgerV0,
    ParcelRetirementCauseV0, ParcelRetirementV0Detail, PlateId,
};
use pf_contracts::spatial::ElementId;
use pf_grid::GridSpec;
use std::collections::BTreeMap;

/// Loud parcel-ledger errors; each names its operation and context.
#[derive(Debug, Clone, PartialEq)]
pub enum ParcelsError {
    UnknownParcel {
        op: &'static str,
        parcel: ParcelId,
        live: u64,
    },
    AlreadyRetired {
        op: &'static str,
        parcel: ParcelId,
        retired_at: f64,
    },
    NotLive {
        op: &'static str,
        parcel: ParcelId,
    },
    DuplicateParcel {
        parcel: ParcelId,
    },
    BadArea {
        op: &'static str,
        element: ElementId,
        area_m2: f64,
    },
    BadTime {
        op: &'static str,
        t_my: f64,
        why: &'static str,
    },
    UnknownPlate {
        op: &'static str,
        plate: PlateId,
    },
    EmptyBirth {
        op: &'static str,
    },
    PositionMapIncomplete {
        missing: Vec<ParcelId>,
    },
    PositionMapExtra {
        extra: Vec<ParcelId>,
    },
    BadElement {
        op: &'static str,
        element: u64,
        cells: usize,
    },
    Conservation(String),
}

impl std::fmt::Display for ParcelsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParcelsError::UnknownParcel { op, parcel, live } => {
                write!(f, "{op}: parcel {parcel:?} not in ledger ({live} records)")
            }
            ParcelsError::AlreadyRetired {
                op,
                parcel,
                retired_at,
            } => {
                write!(
                    f,
                    "{op}: parcel {parcel:?} already retired at {retired_at:?} My"
                )
            }
            ParcelsError::NotLive { op, parcel } => {
                write!(f, "{op}: parcel {parcel:?} is not live")
            }
            ParcelsError::DuplicateParcel { parcel } => write!(f, "duplicate parcel id {parcel:?}"),
            ParcelsError::BadArea {
                op,
                element,
                area_m2,
            } => {
                write!(
                    f,
                    "{op}: element {element:?} has nonfinite/nonpositive area {area_m2:?} m2"
                )
            }
            ParcelsError::BadTime { op, t_my, why } => write!(f, "{op}: t_my {t_my:?}: {why}"),
            ParcelsError::UnknownPlate { op, plate } => write!(f, "{op}: unknown plate {plate:?}"),
            ParcelsError::EmptyBirth { op } => write!(f, "{op}: refusing empty birth batch"),
            ParcelsError::PositionMapIncomplete { missing } => {
                write!(
                    f,
                    "projection: {} live parcels missing positions (first: {:?})",
                    missing.len(),
                    missing.first()
                )
            }
            ParcelsError::PositionMapExtra { extra } => {
                write!(
                    f,
                    "projection: positions given for {} non-live parcels (first: {:?})",
                    extra.len(),
                    extra.first()
                )
            }
            ParcelsError::BadElement { op, element, cells } => {
                write!(f, "{op}: element {element} outside 0..{cells}")
            }
            ParcelsError::Conservation(why) => write!(f, "conservation: {why}"),
        }
    }
}

impl std::error::Error for ParcelsError {}

fn check_time(op: &'static str, t_my: f64) -> Result<(), ParcelsError> {
    if !t_my.is_finite() || t_my < 0.0 {
        return Err(ParcelsError::BadTime {
            op,
            t_my,
            why: "must be finite >= 0",
        });
    }
    Ok(())
}

/// Runtime parcel ledger over the durable records. Records are kept in
/// ascending `parcel_id` order (ids are allocated monotonically, so append
/// preserves order — checked, not assumed). All operations are loud.
#[derive(Debug, Clone, PartialEq)]
pub struct ParcelLedger {
    records: Vec<CrustParcelRecordV0>,
    next_id: u64,
    t_now_my: f64,
}

impl ParcelLedger {
    /// Seed the initial ledger from a plate seeding: exactly one parcel
    /// per grid element, `material_area_equivalent_m2` = the element's
    /// spherical cell area, continental typing from the seeding's
    /// continental plate set. Identity gate: survivors == cells; material
    /// gate: surviving material == grid total area (4*pi*R^2).
    pub fn seed(grid: &GridSpec, seeding: &PlateSeeding) -> Result<Self, ParcelsError> {
        if seeding.assignment.len() != grid.cells() {
            return Err(ParcelsError::Conservation(format!(
                "plate seeding covers {} of {} cells",
                seeding.assignment.len(),
                grid.cells()
            )));
        }
        let mut records = Vec::with_capacity(grid.cells());
        for element in 0..grid.cells() as u64 {
            let (i, _j) = grid.ij(element).map_err(|_| ParcelsError::BadElement {
                op: "seed",
                element,
                cells: grid.cells(),
            })?;
            let flat = element as usize;
            let area = grid
                .cell_area_m2(i)
                .map_err(|e| ParcelsError::Conservation(format!("cell_area_m2: {e}")))?;
            if !area.is_finite() || area <= 0.0 {
                return Err(ParcelsError::BadArea {
                    op: "seed",
                    element: ElementId(element),
                    area_m2: area,
                });
            }
            let plate = seeding.assignment[flat];
            let birth = if seeding.continental_plate_ids.contains(&plate) {
                ParcelBirthV0::InitialContinental
            } else {
                ParcelBirthV0::InitialOceanic
            };
            records.push(CrustParcelRecordV0 {
                parcel_id: ParcelId(element),
                birth_t_my: 0.0,
                birth,
                lineage_parent: None,
                plate,
                material_area_equivalent_m2: area,
                retirement: None,
            });
        }
        let ledger = Self {
            records,
            next_id: grid.cells() as u64,
            t_now_my: 0.0,
        };
        ledger.assert_ordered();
        Ok(ledger)
    }

    fn assert_ordered(&self) {
        for w in self.records.windows(2) {
            debug_assert!(
                w[0].parcel_id < w[1].parcel_id,
                "ledger order invariant broken"
            );
        }
    }

    pub fn t_now_my(&self) -> f64 {
        self.t_now_my
    }

    /// Advance world time. Ages are accrued against this clock; the clock
    /// never runs backward (loud).
    pub fn advance_time(&mut self, dt_my: f64) -> Result<(), ParcelsError> {
        check_time("advance_time", dt_my)?;
        let next = self.t_now_my + dt_my;
        if !next.is_finite() {
            return Err(ParcelsError::BadTime {
                op: "advance_time",
                t_my: next,
                why: "time overflow",
            });
        }
        self.t_now_my = next;
        Ok(())
    }

    /// Accrued age of a parcel at the current clock (never recomputed
    /// from geometry; K2 rev 2).
    pub fn age_my(&self, parcel: ParcelId) -> Result<f64, ParcelsError> {
        let r = self.find(parcel, "age_my")?;
        Ok(self.t_now_my - r.birth_t_my)
    }

    fn find(
        &self,
        parcel: ParcelId,
        op: &'static str,
    ) -> Result<&CrustParcelRecordV0, ParcelsError> {
        self.records
            .binary_search_by_key(&parcel, |r| r.parcel_id)
            .map(|i| &self.records[i])
            .map_err(|_| ParcelsError::UnknownParcel {
                op,
                parcel,
                live: self.records.len() as u64,
            })
    }

    fn find_mut(
        &mut self,
        parcel: ParcelId,
        op: &'static str,
    ) -> Result<&mut CrustParcelRecordV0, ParcelsError> {
        match self.records.binary_search_by_key(&parcel, |r| r.parcel_id) {
            Ok(i) => Ok(&mut self.records[i]),
            Err(_) => Err(ParcelsError::UnknownParcel {
                op,
                parcel,
                live: self.records.len() as u64,
            }),
        }
    }

    pub fn records(&self) -> &[CrustParcelRecordV0] {
        &self.records
    }

    pub fn live_count(&self) -> u64 {
        self.records
            .iter()
            .filter(|r| r.retirement.is_none())
            .count() as u64
    }

    /// New crust at a ridge (K2): one parcel per spreading cell, birth
    /// time = now, thermal/crust age zero, material = the cell areas.
    /// Returns the minted ids in cell order.
    pub fn ridge_birth(
        &mut self,
        cells: &[(ElementId, f64, PlateId)],
    ) -> Result<Vec<ParcelId>, ParcelsError> {
        if cells.is_empty() {
            return Err(ParcelsError::EmptyBirth { op: "ridge_birth" });
        }
        let t = self.t_now_my;
        let mut minted = Vec::with_capacity(cells.len());
        for &(element, area, plate) in cells {
            if !area.is_finite() || area <= 0.0 {
                return Err(ParcelsError::BadArea {
                    op: "ridge_birth",
                    element,
                    area_m2: area,
                });
            }
            let id = ParcelId(self.next_id);
            self.next_id += 1;
            self.records.push(CrustParcelRecordV0 {
                parcel_id: id,
                birth_t_my: t,
                birth: ParcelBirthV0::RidgeCreated,
                lineage_parent: None,
                plate,
                material_area_equivalent_m2: area,
                retirement: None,
            });
            minted.push(id);
        }
        Ok(minted)
    }

    /// Hotspot-created crust (K6 island chains): same shape as ridge birth
    /// with its own birth kind (new volcanic material, no parent rock).
    pub fn hotspot_birth(
        &mut self,
        cells: &[(ElementId, f64, PlateId)],
    ) -> Result<Vec<ParcelId>, ParcelsError> {
        if cells.is_empty() {
            return Err(ParcelsError::EmptyBirth {
                op: "hotspot_birth",
            });
        }
        let t = self.t_now_my;
        let mut minted = Vec::with_capacity(cells.len());
        for &(element, area, plate) in cells {
            if !area.is_finite() || area <= 0.0 {
                return Err(ParcelsError::BadArea {
                    op: "hotspot_birth",
                    element,
                    area_m2: area,
                });
            }
            let id = ParcelId(self.next_id);
            self.next_id += 1;
            self.records.push(CrustParcelRecordV0 {
                parcel_id: id,
                birth_t_my: t,
                birth: ParcelBirthV0::HotspotCreated,
                lineage_parent: None,
                plate,
                material_area_equivalent_m2: area,
                retirement: None,
            });
            minted.push(id);
        }
        Ok(minted)
    }

    /// Accretion (collision transfer to the overriding plate): new parcels
    /// whose lineage parent is the source rock. Parents must exist and be
    /// born no later than the children (checked at snapshot validation,
    /// and loudly here for immediate feedback).
    pub fn accrete(
        &mut self,
        children: &[(ElementId, f64, PlateId, ParcelId)],
    ) -> Result<Vec<ParcelId>, ParcelsError> {
        if children.is_empty() {
            return Err(ParcelsError::EmptyBirth { op: "accrete" });
        }
        let t = self.t_now_my;
        let mut minted = Vec::with_capacity(children.len());
        for &(element, area, plate, parent) in children {
            if !area.is_finite() || area <= 0.0 {
                return Err(ParcelsError::BadArea {
                    op: "accrete",
                    element,
                    area_m2: area,
                });
            }
            let pr = self.find(parent, "accrete.parent")?;
            if pr.birth_t_my > t {
                return Err(ParcelsError::BadTime {
                    op: "accrete",
                    t_my: pr.birth_t_my,
                    why: "parent born after accretion time (age monotonicity)",
                });
            }
            let id = ParcelId(self.next_id);
            self.next_id += 1;
            self.records.push(CrustParcelRecordV0 {
                parcel_id: id,
                birth_t_my: t,
                birth: ParcelBirthV0::Accreted,
                lineage_parent: Some(parent),
                plate,
                material_area_equivalent_m2: area,
                retirement: None,
            });
            minted.push(id);
        }
        Ok(minted)
    }

    /// Retire parcels into the subduction sink (K2): records WHEN, WHY and
    /// the ledger event that did it. Retiring an already-retired parcel is
    /// a loud error, never a no-op.
    pub fn retire_subducted(
        &mut self,
        parcels: &[ParcelId],
        event: EventId,
    ) -> Result<(), ParcelsError> {
        let t = self.t_now_my;
        for &p in parcels {
            let r = self.find_mut(p, "retire_subducted")?;
            if let Some(old) = r.retirement {
                return Err(ParcelsError::AlreadyRetired {
                    op: "retire_subducted",
                    parcel: p,
                    retired_at: old.t_my,
                });
            }
            r.retirement = Some(ParcelRetirementV0Detail {
                t_my: t,
                cause: ParcelRetirementCauseV0::Subducted,
                event,
            });
        }
        Ok(())
    }

    /// Boundary sweep (K2): live parcels change owning plate in place —
    /// identity, birth and lineage are untouched (the EVENT ledger records
    /// the sweep; this is the material side of it).
    pub fn reassign_plate(
        &mut self,
        parcels: &[ParcelId],
        new_plate: PlateId,
    ) -> Result<(), ParcelsError> {
        for &p in parcels {
            let r = self.find_mut(p, "reassign_plate")?;
            if r.retirement.is_some() {
                return Err(ParcelsError::NotLive {
                    op: "reassign_plate",
                    parcel: p,
                });
            }
            r.plate = new_plate;
        }
        Ok(())
    }

    /// Lineage chain: parcel back through parents to its initial/new-material
    /// ancestor (the wake-up acceptance story: "born at ridge X, transported
    /// on plate 4, accreted in collision Y, exposed here").
    pub fn lineage_chain(&self, parcel: ParcelId) -> Result<Vec<ParcelId>, ParcelsError> {
        let mut chain = vec![parcel];
        loop {
            let r = self.find(*chain.last().expect("nonempty"), "lineage_chain")?;
            match r.lineage_parent {
                Some(parent) if !chain.contains(&parent) => chain.push(parent),
                Some(parent) => {
                    return Err(ParcelsError::Conservation(format!(
                        "lineage cycle at parcel {parent:?} (chain {chain:?})"
                    )))
                }
                None => return Ok(chain),
            }
        }
    }

    /// DERIVED weighted projection (PARCEL != CELL): B1 nearest-cell — one
    /// contribution per live parcel at weight 1.0. `positions` must cover
    /// exactly the live parcels; many parcels per cell is legal; retired
    /// parcels with positions, or live parcels without, are loud errors.
    /// A later conservative/PIC remapper supplies its own split map; the
    /// ledger is not consulted about space at all beyond this projection.
    pub fn project_onto_cells(
        &self,
        positions: &BTreeMap<ParcelId, ElementId>,
        cells: usize,
    ) -> Result<Vec<CellParcelContribution>, ParcelsError> {
        let mut missing = Vec::new();
        for r in &self.records {
            if r.retirement.is_none() && !positions.contains_key(&r.parcel_id) {
                missing.push(r.parcel_id);
            }
        }
        if !missing.is_empty() {
            return Err(ParcelsError::PositionMapIncomplete {
                missing: missing[..5.min(missing.len())].to_vec(),
            });
        }
        let mut extra = Vec::new();
        for &pid in positions.keys() {
            match self.records.binary_search_by_key(&pid, |r| r.parcel_id) {
                Ok(i) if self.records[i].retirement.is_none() => {}
                Ok(_) | Err(_) => extra.push(pid),
            }
        }
        if !extra.is_empty() {
            return Err(ParcelsError::PositionMapExtra {
                extra: extra[..5.min(extra.len())].to_vec(),
            });
        }
        let mut contributions: Vec<CellParcelContribution> = positions
            .iter()
            .map(|(&pid, &el)| {
                if el.0 as usize >= cells {
                    return Err(ParcelsError::BadElement {
                        op: "project_onto_cells",
                        element: el.0,
                        cells,
                    });
                }
                Ok(CellParcelContribution {
                    element: el,
                    parcel_id: pid,
                    weight: 1.0,
                })
            })
            .collect::<Result<_, _>>()?;
        contributions.sort_by_key(|c| (c.element, c.parcel_id));
        Ok(contributions)
    }

    /// Durable snapshot with BOTH conservation gates recomputed from the
    /// records in canonical parcel_id order (the domain validates this
    /// snapshot; conservation.rs audits it).
    pub fn to_domain(&self, declared_error_m2: f64) -> Result<ParcelLedgerV0, ParcelsError> {
        crate::conservation::snapshot(&self.records, self.t_now_my, declared_error_m2)
            .map_err(ParcelsError::Conservation)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::{seed_plates, SplitMix64};
    use pf_contracts::params::{
        FidelityPreset, MantleBlock, PlanetParameterProfile, ScenarioPreset, SeedPath,
    };

    fn seeded() -> (GridSpec, PlateSeeding, ParcelLedger) {
        let grid = GridSpec::new(10, 18, 6.371e6).unwrap();
        let mut p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
        p.tectonics.plate_count = 5;
        p.tectonics.continental_seeds = 2;
        let sp = SeedPath {
            seed: 240907,
            derivation: "toy".into(),
            salt: "b-test".into(),
        };
        let seeding = seed_plates(&grid, &p, &sp).unwrap();
        let ledger = ParcelLedger::seed(&grid, &seeding).unwrap();
        (grid, seeding, ledger)
    }

    fn positions_all(ledger: &ParcelLedger, cells: usize) -> BTreeMap<ParcelId, ElementId> {
        // Deterministic stand-in for the driver's advection: parcel i sits
        // at element (i % cells).
        let mut m = BTreeMap::new();
        for (k, r) in ledger.records().iter().enumerate() {
            if r.retirement.is_none() {
                m.insert(r.parcel_id, ElementId((k as u64) % cells as u64));
            }
        }
        m
    }

    #[test]
    fn seed_closes_both_gates_against_the_sphere() {
        let (grid, _seeding, ledger) = seeded();
        let dom = ledger.to_domain(1.0e-6).unwrap();
        dom.validate().unwrap();
        assert_eq!(dom.identity.survivors, grid.cells() as u64);
        assert_eq!(dom.identity.initial_parcels, grid.cells() as u64);
        assert_eq!(dom.identity.births, 0);
        // Material = total sphere area (the grid guarantees sum = 4 pi R^2).
        let total = grid.total_area_m2();
        assert!((dom.material.surviving_material - total).abs() / total < 1.0e-9);
        assert_eq!(dom.material.quantity, "m2");
        // Continental typing actually applied.
        let continental = ledger
            .records()
            .iter()
            .filter(|r| matches!(r.birth, ParcelBirthV0::InitialContinental))
            .count();
        assert!(
            continental > 0,
            "continental seeds produced no continental parcels"
        );
        assert!(
            continental < grid.cells(),
            "everything continental — typing broken"
        );
    }

    #[test]
    fn ridge_birth_and_retirement_move_both_gates() {
        let (grid, _s, mut ledger) = seeded();
        ledger.advance_time(2.0).unwrap();
        let area = grid.cell_area_m2(3).unwrap();
        let minted = ledger
            .ridge_birth(&[(ElementId(30), area, pf_contracts::domain::PlateId(1))])
            .unwrap();
        assert_eq!(minted.len(), 1);
        let dom = ledger.to_domain(1.0e-6).unwrap();
        dom.validate().unwrap();
        assert_eq!(dom.identity.births, 1);
        assert_eq!(dom.identity.survivors, grid.cells() as u64 + 1);
        assert!((dom.material.created_material - area).abs() < 1.0e-6);

        // Retire it: identity exact, material moves out of surviving.
        ledger.advance_time(1.0).unwrap();
        ledger.retire_subducted(&minted, EventId(7)).unwrap();
        let dom = ledger.to_domain(1.0e-6).unwrap();
        dom.validate().unwrap();
        assert_eq!(dom.identity.survivors, grid.cells() as u64);
        assert_eq!(dom.identity.retirements, 1);
        assert!((dom.material.retired_exported_material - area).abs() < 1.0e-6);
        assert!(
            (dom.material.surviving_material - grid.total_area_m2()).abs() / grid.total_area_m2()
                < 1.0e-9
        );
        // Age accrued, birth fixed: the retired parcel's birth is unchanged.
        assert!((ledger.records()[minted[0].0 as usize].birth_t_my - 2.0).abs() == 0.0);

        // Double retirement is loud, not a no-op.
        let err = ledger.retire_subducted(&minted, EventId(8)).unwrap_err();
        assert!(err.to_string().contains("already retired"), "{err}");
    }

    #[test]
    fn two_parcels_one_cell_neither_disappears() {
        let (grid, _s, ledger) = seeded();
        let cells = grid.cells();
        // Fold every parcel into HALF the cells: many-to-one, all weight 1.
        let mut positions = BTreeMap::new();
        for (k, r) in ledger.records().iter().enumerate() {
            positions.insert(r.parcel_id, ElementId((k as u64) % (cells as u64 / 2)));
        }
        let contribs = ledger.project_onto_cells(&positions, cells).unwrap();
        assert_eq!(contribs.len(), cells);
        // Per-parcel weights close to 1 (they are exactly 1 here).
        for c in &contribs {
            assert_eq!(c.weight, 1.0);
        }
        // Conservation is untouched by the projection (structural claim,
        // tested by NOT involving positions in the snapshot).
        let dom = ledger.to_domain(1.0e-6).unwrap();
        assert_eq!(dom.identity.survivors, cells as u64);
    }

    #[test]
    fn projection_map_must_cover_exactly_the_live_parcels() {
        let (grid, _s, mut ledger) = seeded();
        let cells = grid.cells();
        // Incomplete: drop one live parcel's position.
        let mut bad = positions_all(&ledger, cells);
        let a_live = *bad.keys().next().unwrap();
        bad.remove(&a_live);
        let err = ledger.project_onto_cells(&bad, cells).unwrap_err();
        assert!(
            matches!(err, ParcelsError::PositionMapIncomplete { .. }),
            "{err}"
        );

        // Extra: add a retired parcel's position after retiring it.
        ledger.advance_time(1.0).unwrap();
        ledger.retire_subducted(&[a_live], EventId(1)).unwrap();
        let mut with_retired = positions_all(&ledger, cells);
        with_retired.insert(a_live, ElementId(0));
        let err = ledger.project_onto_cells(&with_retired, cells).unwrap_err();
        assert!(
            matches!(err, ParcelsError::PositionMapExtra { .. }),
            "{err}"
        );

        // Out-of-range element is loud (use a still-live parcel).
        let live_id = ledger.records()[1].parcel_id;
        let mut oob = positions_all(&ledger, cells);
        oob.insert(live_id, ElementId(cells as u64 + 5));
        let err = ledger.project_onto_cells(&oob, cells).unwrap_err();
        assert!(matches!(err, ParcelsError::BadElement { .. }), "{err}");
    }

    #[test]
    fn lineage_chain_walks_to_the_initial_ancestor() {
        let (grid, _s, mut ledger) = seeded();
        ledger.advance_time(3.0).unwrap();
        let a0 = ledger.records()[0].parcel_id;
        let child = ledger
            .accrete(&[(ElementId(1), grid.cell_area_m2(0).unwrap(), PlateId(1), a0)])
            .unwrap()[0];
        let grandchild = ledger
            .accrete(&[(
                ElementId(2),
                grid.cell_area_m2(0).unwrap(),
                PlateId(2),
                child,
            )])
            .unwrap()[0];
        let chain = ledger.lineage_chain(grandchild).unwrap();
        assert_eq!(chain, vec![grandchild, child, a0]);
    }

    #[test]
    fn bad_inputs_are_loud() {
        let (_grid, _s, mut ledger) = seeded();
        // Nonfinite area.
        let err = ledger
            .ridge_birth(&[(ElementId(1), f64::NAN, PlateId(1))])
            .unwrap_err();
        assert!(matches!(err, ParcelsError::BadArea { .. }), "{err}");
        // Empty batch refused (a no-op birth would hide driver bugs).
        assert!(matches!(
            ledger.ridge_birth(&[]).unwrap_err(),
            ParcelsError::EmptyBirth { .. }
        ));
        // Unknown parent for accretion.
        let err = ledger
            .accrete(&[(ElementId(1), 1.0, PlateId(1), ParcelId(999_999))])
            .unwrap_err();
        assert!(matches!(err, ParcelsError::UnknownParcel { .. }), "{err}");
        // Time never runs backward.
        assert!(ledger.advance_time(-1.0).is_err());
        // Reassigning a retired parcel is loud.
        ledger.advance_time(1.0).unwrap();
        let p0 = ledger.records()[0].parcel_id;
        ledger.retire_subducted(&[p0], EventId(2)).unwrap();
        let err = ledger.reassign_plate(&[p0], PlateId(3)).unwrap_err();
        assert!(matches!(err, ParcelsError::NotLive { .. }), "{err}");
        // Unknown parcel lookup names the op.
        let err = ledger.age_my(ParcelId(424_242)).unwrap_err();
        assert!(err.to_string().contains("age_my"), "{err}");
    }

    #[test]
    fn five_step_sequence_closes_per_step_and_cumulative() {
        let (grid, _s, mut ledger) = seeded();
        let cells = grid.cells() as u64;
        let area = grid.cell_area_m2(0).unwrap();
        // Whole-planet magnitudes (~5e14 m2) put f64 summation granularity
        // near 0.03 m2: a 1 m2 declared error is an honest numerical budget
        // here, not a hidden tolerance.
        let mut history = crate::conservation::ConservationHistory::new(0.0);
        let mut minted_prev = ParcelId(0); // assigned at step 0 before any use
        let mut step = 0;
        while step < 5 {
            ledger.advance_time(2.0).unwrap();
            let t = ledger.t_now_my();
            let minted = ledger
                .ridge_birth(&[
                    (ElementId(step as u64), area, PlateId(1)),
                    (ElementId(step as u64 + 1), area, PlateId(2)),
                ])
                .unwrap();
            // Retire the previous step's first mint (and one initial parcel
            // on the last step) so retirements are nonzero every step after
            // the first.
            if step >= 1 {
                ledger
                    .retire_subducted(&[minted_prev], EventId(step as u64 + 100))
                    .unwrap();
            }
            if step == 4 {
                ledger
                    .retire_subducted(&[ParcelId(0)], EventId(200))
                    .unwrap();
            }
            minted_prev = minted[0];
            let snap = ledger.to_domain(1.0).unwrap();
            history.push_row(&snap, t).unwrap();
            step += 1;
        }
        history.close(ledger.to_domain(1.0).unwrap()).unwrap();
        assert_eq!(history.rows().len(), 5);
        // Identity bookkeeping is exact integers end to end.
        let last = history.rows().last().unwrap();
        assert_eq!(
            last.survivors,
            cells + 10 - 4 - 1,
            "net: +2 births x5 steps, -1 retirement on steps 1-4, -1 extra initial on the last"
        );
        // Cumulative material closes within declared error.
        assert!(history.cumulative().closes());
        // A lying row is caught: tamper with one row's survivors.
        let mut lied = history.clone();
        if let Some(row) = lied.rows_mut().last_mut() {
            row.survivors += 3;
        }
        assert!(lied.close(ledger.to_domain(1.0e-6).unwrap()).is_err());
    }

    // Determinism cross-check: same seed path -> identical ledger history.
    #[test]
    fn ledger_is_deterministic_from_the_seed_path() {
        fn build(salt: &str) -> (u64, f64) {
            let grid = GridSpec::new(10, 18, 6.371e6).unwrap();
            let mut p: PlanetParameterProfile =
                ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
            p.tectonics.plate_count = 5;
            let sp = SeedPath {
                seed: 99,
                derivation: "det".into(),
                salt: salt.into(),
            };
            let seeding = seed_plates(&grid, &p, &sp).unwrap();
            let mut ledger = ParcelLedger::seed(&grid, &seeding).unwrap();
            ledger.advance_time(4.0).unwrap();
            let area = grid.cell_area_m2(1).unwrap();
            let ids = ledger
                .ridge_birth(&[(ElementId(5), area, PlateId(0))])
                .unwrap();
            ledger.retire_subducted(&ids, EventId(3)).unwrap();
            let dom = ledger.to_domain(1.0e-6).unwrap();
            (dom.identity.survivors, dom.material.surviving_material)
        }
        assert_eq!(build("x"), build("x"));
        // And the RNG stream itself replays (shared plumbing sanity).
        assert_eq!(SplitMix64::new(7).next_u64(), SplitMix64::new(7).next_u64());
    }
}
