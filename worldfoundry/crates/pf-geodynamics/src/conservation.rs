//! The two conservation gates (brief T1; spec §54.4-style for crust
//! parcels; overnight agent-B acceptance: per step AND cumulative, error
//! reported, never silently grown).
//!
//! GATE 1 — IDENTITY (exact, integer domain):
//!   initial + births - retirements == survivors
//! Counts are counts. There is no remap tolerance here; a mismatch is a
//! ledger bug, full stop.
//!
//! GATE 2 — MATERIAL (physical, declared numerical error):
//!   initial + created - retired/exported == surviving +- declared_error
//! in `material_area_equivalent_m2` — the conserved Lagrangian weight
//! (NOT the current geometric footprint; shortening/thickening changes
//! footprints, not this). Sums are taken in canonical parcel_id order so
//! the writer's and the verifier's floats agree bit-for-bit.
//!
//! `ConservationHistory` additionally enforces step-to-step monotonic
//! bookkeeping: births and retirements never DECREASE between snapshots
//! (a ledger that un-retires rock has rewritten history), initial counts
//! are constant, and survivor deltas close exactly against birth and
//! retirement deltas. Per-step rows + a cumulative row are the loud
//! report; nothing here silently widens a tolerance.

use pf_contracts::domain::{
    CrustParcelRecordV0, MaterialConservationV0, ParcelBirthV0, ParcelIdentityConservationV0,
    ParcelLedgerV0,
};

/// Recompute both gates from records in canonical (ascending parcel_id)
/// order and build a self-consistent durable snapshot. The snapshot cannot
/// lie: every number in it is derived from the records here.
pub fn snapshot(
    records: &[CrustParcelRecordV0],
    t_my: f64,
    declared_error_m2: f64,
) -> Result<ParcelLedgerV0, String> {
    if !t_my.is_finite() || t_my < 0.0 {
        return Err(format!("snapshot: t_my {t_my:?} must be finite >= 0"));
    }
    if !declared_error_m2.is_finite() || declared_error_m2 < 0.0 {
        return Err(format!(
            "snapshot: declared_error_m2 {declared_error_m2:?} must be finite >= 0"
        ));
    }
    let mut ordered: Vec<&CrustParcelRecordV0> = records.iter().collect();
    ordered.sort_by_key(|r| r.parcel_id);
    for w in ordered.windows(2) {
        if w[0].parcel_id == w[1].parcel_id {
            return Err(format!(
                "snapshot: duplicate parcel_id {:?}",
                w[0].parcel_id
            ));
        }
    }
    let mut initial: u64 = 0;
    let mut births: u64 = 0;
    let mut retirements: u64 = 0;
    let (mut m_initial, mut m_created, mut m_retired, mut m_surviving) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for r in &ordered {
        if !r.material_area_equivalent_m2.is_finite() || r.material_area_equivalent_m2 <= 0.0 {
            return Err(format!(
                "snapshot: parcel {:?} has nonfinite/nonpositive material {:?}",
                r.parcel_id, r.material_area_equivalent_m2
            ));
        }
        match r.birth {
            ParcelBirthV0::InitialOceanic | ParcelBirthV0::InitialContinental => {
                initial += 1;
                m_initial += r.material_area_equivalent_m2;
            }
            _ => {
                births += 1;
                m_created += r.material_area_equivalent_m2;
            }
        }
        if r.retirement.is_some() {
            retirements += 1;
            m_retired += r.material_area_equivalent_m2;
        } else {
            m_surviving += r.material_area_equivalent_m2;
        }
    }
    let survivors = ordered.len() as u64 - retirements;
    let identity = ParcelIdentityConservationV0 {
        t_my,
        initial_parcels: initial,
        births,
        retirements,
        survivors,
    };
    if !identity.identity_holds() {
        return Err(format!(
            "snapshot: identity gate failed during construction: {} + {} - {} != {survivors}",
            identity.initial_parcels, identity.births, identity.retirements
        ));
    }
    let material = MaterialConservationV0 {
        t_my,
        quantity: "m2".into(),
        initial_material: m_initial,
        created_material: m_created,
        retired_exported_material: m_retired,
        surviving_material: m_surviving,
        declared_error: declared_error_m2,
    };
    Ok(ParcelLedgerV0 {
        parcels: ordered.into_iter().cloned().collect(),
        identity,
        material,
    })
}

/// One full gate snapshot at one committed step (absolutes, not deltas:
/// each row independently re-derivable from the records at that time).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConservationRow {
    pub t_my: f64,
    pub initial_parcels: u64,
    pub births: u64,
    pub retirements: u64,
    pub survivors: u64,
    pub initial_material_m2: f64,
    pub created_material_m2: f64,
    pub retired_material_m2: f64,
    pub surviving_material_m2: f64,
    pub declared_error_m2: f64,
    /// |initial + created - retired - surviving| at this step (report
    /// artifact; must be <= declared_error — checked).
    pub material_residual_m2: f64,
}

impl ConservationRow {
    fn from_snapshot(snap: &ParcelLedgerV0) -> Self {
        let m = &snap.material;
        ConservationRow {
            t_my: m.t_my,
            initial_parcels: snap.identity.initial_parcels,
            births: snap.identity.births,
            retirements: snap.identity.retirements,
            survivors: snap.identity.survivors,
            initial_material_m2: m.initial_material,
            created_material_m2: m.created_material,
            retired_material_m2: m.retired_exported_material,
            surviving_material_m2: m.surviving_material,
            declared_error_m2: m.declared_error,
            material_residual_m2: (m.initial_material + m.created_material
                - m.retired_exported_material
                - m.surviving_material)
                .abs(),
        }
    }
}

/// Cumulative material closure (the last committed row's absolutes).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CumulativeMaterial {
    pub initial_m2: f64,
    pub created_m2: f64,
    pub retired_m2: f64,
    pub surviving_m2: f64,
    pub declared_error_m2: f64,
}

impl CumulativeMaterial {
    pub fn closes(&self) -> bool {
        let residual = self.initial_m2 + self.created_m2 - self.retired_m2 - self.surviving_m2;
        residual.is_finite() && residual.abs() <= self.declared_error_m2
    }

    pub fn residual_m2(&self) -> f64 {
        (self.initial_m2 + self.created_m2 - self.retired_m2 - self.surviving_m2).abs()
    }
}

/// Per-step rows + cumulative closure, with step-to-step monotonic
/// bookkeeping enforced (history cannot be rewritten between snapshots).
#[derive(Debug, Clone, PartialEq)]
pub struct ConservationHistory {
    rows: Vec<ConservationRow>,
}

impl ConservationHistory {
    pub fn new(_t_start_my: f64) -> Self {
        Self { rows: Vec::new() }
    }

    pub fn rows(&self) -> &[ConservationRow] {
        &self.rows
    }

    pub fn rows_mut(&mut self) -> &mut [ConservationRow] {
        &mut self.rows
    }

    /// Commit one step: the snapshot must be internally valid (the domain
    /// gate), the step's material residual within declared error, and the
    /// bookkeeping monotone against the previous row.
    pub fn push_row(&mut self, snap: &ParcelLedgerV0, t_my: f64) -> Result<(), String> {
        snap.validate()
            .map_err(|e| format!("row {}: snapshot invalid: {e}", self.rows.len()))?;
        let row = ConservationRow::from_snapshot(snap);
        if (row.t_my - t_my).abs() > 1.0e-12 {
            return Err(format!(
                "row {}: t_my {t_my:?} does not match snapshot {:?}",
                self.rows.len(),
                row.t_my
            ));
        }
        if row.material_residual_m2 > row.declared_error_m2 {
            return Err(format!(
                "row {} at {:?} My: material residual {:?} exceeds declared error {:?}",
                self.rows.len(),
                row.t_my,
                row.material_residual_m2,
                row.declared_error_m2
            ));
        }
        if let Some(prev) = self.rows.last() {
            if row.t_my <= prev.t_my {
                return Err(format!(
                    "row {}: time {:?} not after previous {:?}",
                    self.rows.len(),
                    row.t_my,
                    prev.t_my
                ));
            }
            if row.initial_parcels != prev.initial_parcels {
                return Err(format!(
                    "row {}: initial_parcels changed {} -> {} (history rewrite)",
                    self.rows.len(),
                    prev.initial_parcels,
                    row.initial_parcels
                ));
            }
            if row.births < prev.births || row.retirements < prev.retirements {
                return Err(format!(
                    "row {}: monotonicity violated (births {} -> {}, retirements {} -> {}): rock was un-born or un-retired",
                    self.rows.len(), prev.births, row.births, prev.retirements, row.retirements
                ));
            }
            let d_survivors = row.survivors as i64 - prev.survivors as i64;
            let d_expected =
                (row.births - prev.births) as i64 - (row.retirements - prev.retirements) as i64;
            if d_survivors != d_expected {
                return Err(format!(
                    "row {}: survivor delta {d_survivors} != births delta - retirements delta ({d_expected})",
                    self.rows.len()
                ));
            }
            let d_mat = row.surviving_material_m2 - prev.surviving_material_m2;
            let d_mat_expected = (row.created_material_m2 - prev.created_material_m2)
                - (row.retired_material_m2 - prev.retired_material_m2);
            if (d_mat - d_mat_expected).abs() > row.declared_error_m2 + prev.declared_error_m2 {
                return Err(format!(
                    "row {}: material delta {d_mat:?} != created/retired deltas {d_mat_expected:?} beyond declared errors",
                    self.rows.len()
                ));
            }
        }
        self.rows.push(row);
        Ok(())
    }

    /// Cumulative closure from the last committed row.
    pub fn cumulative(&self) -> CumulativeMaterial {
        let last = self
            .rows
            .last()
            .expect("cumulative on empty history is a caller bug (fail loud in tests)");
        CumulativeMaterial {
            initial_m2: last.initial_material_m2,
            created_m2: last.created_material_m2,
            retired_m2: last.retired_material_m2,
            surviving_m2: last.surviving_material_m2,
            declared_error_m2: last.declared_error_m2,
        }
    }

    /// Close the history against the FINAL ledger state: last row must
    /// equal the final snapshot's gates exactly (counts) and within
    /// declared error (material).
    pub fn close(&self, final_snap: ParcelLedgerV0) -> Result<(), String> {
        let last = self
            .rows
            .last()
            .ok_or_else(|| "close: empty history".to_string())?;
        final_snap
            .validate()
            .map_err(|e| format!("close: final snapshot invalid: {e}"))?;
        let f = &final_snap.identity;
        if (
            last.initial_parcels,
            last.births,
            last.retirements,
            last.survivors,
        ) != (f.initial_parcels, f.births, f.retirements, f.survivors)
        {
            return Err(format!(
                "close: last row ({}, {}, {}, {}) != final snapshot ({}, {}, {}, {})",
                last.initial_parcels,
                last.births,
                last.retirements,
                last.survivors,
                f.initial_parcels,
                f.births,
                f.retirements,
                f.survivors
            ));
        }
        let m = &final_snap.material;
        let drift = |a: f64, b: f64, name: &str| -> Result<(), String> {
            if (a - b).abs() > m.declared_error {
                return Err(format!(
                    "close: {name} drifted {a:?} vs {b:?} beyond {:?}",
                    m.declared_error
                ));
            }
            Ok(())
        };
        drift(
            last.initial_material_m2,
            m.initial_material,
            "initial_material",
        )?;
        drift(
            last.created_material_m2,
            m.created_material,
            "created_material",
        )?;
        drift(
            last.retired_material_m2,
            m.retired_exported_material,
            "retired_material",
        )?;
        drift(
            last.surviving_material_m2,
            m.surviving_material,
            "surviving_material",
        )?;
        if !self.cumulative().closes() {
            return Err(format!(
                "close: cumulative residual {:?} exceeds declared error {:?}",
                self.cumulative().residual_m2(),
                last.declared_error_m2
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_contracts::domain::{ParcelBirthV0, ParcelId, PlateId};

    fn rec(id: u64, birth: ParcelBirthV0, area: f64, retired: Option<f64>) -> CrustParcelRecordV0 {
        CrustParcelRecordV0 {
            parcel_id: ParcelId(id),
            birth_t_my: 0.0,
            birth,
            lineage_parent: None,
            plate: PlateId(1),
            material_area_equivalent_m2: area,
            retirement: retired.map(|t| pf_contracts::domain::ParcelRetirementV0Detail {
                t_my: t,
                cause: pf_contracts::domain::ParcelRetirementCauseV0::Subducted,
                event: pf_contracts::domain::EventId(1),
            }),
        }
    }

    #[test]
    fn snapshot_computes_both_gates_in_canonical_order() {
        // Deliberately unordered input: canonical summation must not care.
        let records = vec![
            rec(5, ParcelBirthV0::RidgeCreated, 3.0, Some(8.0)),
            rec(1, ParcelBirthV0::InitialOceanic, 10.0, None),
            rec(3, ParcelBirthV0::InitialContinental, 7.5, None),
            rec(9, ParcelBirthV0::RidgeCreated, 2.5, None),
        ];
        let snap = snapshot(&records, 10.0, 1.0e-9).unwrap();
        snap.validate().unwrap();
        assert_eq!(snap.identity.initial_parcels, 2);
        assert_eq!(snap.identity.births, 2);
        assert_eq!(snap.identity.retirements, 1);
        assert_eq!(snap.identity.survivors, 3);
        assert!((snap.material.initial_material - 17.5).abs() < 1.0e-12);
        assert!((snap.material.created_material - 5.5).abs() < 1.0e-12);
        assert!((snap.material.retired_exported_material - 3.0).abs() < 1.0e-12);
        assert!((snap.material.surviving_material - 20.0).abs() < 1.0e-12);
    }

    #[test]
    fn snapshot_refuses_bad_inputs_loudly() {
        assert!(snapshot(
            &[rec(1, ParcelBirthV0::InitialOceanic, -1.0, None)],
            0.0,
            0.0
        )
        .is_err());
        let dup = vec![
            rec(1, ParcelBirthV0::InitialOceanic, 1.0, None),
            rec(1, ParcelBirthV0::InitialOceanic, 1.0, None),
        ];
        assert!(snapshot(&dup, 0.0, 0.0).is_err());
        assert!(snapshot(&[], -1.0, 0.0).is_err());
        assert!(snapshot(&[], 0.0, f64::NAN).is_err());
    }

    #[test]
    fn history_rejects_rewrites_and_bad_deltas() {
        let mut h = ConservationHistory::new(0.0);
        let s0 = snapshot(
            &[rec(1, ParcelBirthV0::InitialOceanic, 10.0, None)],
            1.0,
            1.0e-9,
        )
        .unwrap();
        h.push_row(&s0, 1.0).unwrap();

        // Birth count went DOWN between snapshots: history rewrite.
        let s_bad = snapshot(&[], 2.0, 1.0e-9).unwrap();
        let err = h.push_row(&s_bad, 2.0).unwrap_err();
        assert!(err.contains("initial_parcels changed"), "{err}");

        // Survivor delta inconsistent with births/retirements deltas:
        // fabricate by hand (snapshot itself cannot lie).
        let mut lied = snapshot(
            &[rec(1, ParcelBirthV0::InitialOceanic, 10.0, None)],
            2.0,
            1.0e-9,
        )
        .unwrap();
        lied.identity.survivors += 4;
        let err = h.push_row(&lied, 2.0).unwrap_err();
        assert!(err.contains("snapshot invalid"), "{err}");

        // Time must advance (same t as the previous row is refused).
        let s1 = snapshot(
            &[
                rec(1, ParcelBirthV0::InitialOceanic, 10.0, None),
                rec(2, ParcelBirthV0::RidgeCreated, 1.0, None),
            ],
            1.0,
            1.0e-9,
        )
        .unwrap();
        let err = h.push_row(&s1, 1.0).unwrap_err();
        assert!(err.contains("not after previous"), "{err}");
    }
}
