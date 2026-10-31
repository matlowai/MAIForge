//! Opt-in `native-formation-age-v1` observations of the native crust solver.
//!
//! A handle identifies an accepted numerical origin cohort, NOT a unique rock,
//! conserved ocean parcel, or canonical body. Nearest-source remapping can copy
//! one cohort into several columns. Volumes below are receipts at creation, not
//! a claim that the subsequent Eulerian transport conserves cohort mass.
//! Initial thermal ages do not establish formation ages. Projection births are
//! numerical column initialization, not observations of basalt crystallization.
//! Time intervals span the actual operation; no substep chronology is inferred.
//! Coverage is the basal ocean column origin. Later arc/plume additions, layer
//! mixing, erosion and excavated rock identity require separate producer wiring;
//! a column origin must never be assigned indiscriminately to every rock in it.

use crate::kinematic::{CrustState, NO_SOURCE};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;
use serde_json::{json, Value};

pub const VERSION: &str = "native-formation-age-v1";

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimeInterval {
    pub start_my: f64,
    pub end_my: f64,
}
impl TimeInterval {
    pub fn new(start_my: f64, end_my: f64) -> Result<Self, String> {
        if !start_my.is_finite() || !end_my.is_finite() || end_my <= start_my {
            return Err("formation-age: expected finite increasing operation times".into());
        }
        Ok(Self { start_my, end_my })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BirthOrigin {
    /// Accepted gap under the selected kinematic law. Not necessarily a
    /// resolved divergent ridge: unfilled transform gaps also use that law.
    KinematicGap,
    /// Continental parcel projection vacated this column without an ocean donor.
    ProjectionInitialization,
}
/// A filled provisional gap can later be vacated by continental projection.
/// Retain the cancelled proposal independently of the accepted creation kind.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BirthCandidate {
    pub origin: BirthOrigin,
    pub provisional_gap: bool,
}

impl BirthOrigin {
    fn label(self) -> &'static str {
        match self {
            Self::KinematicGap => "kinematic-gap",
            Self::ProjectionInitialization => "projection-initialization",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Formation {
    UnknownInherited,
    NotOceanic,
    /// One-based index into this recorder's event archive. Namespace is the
    /// pinned run/product supplied by the producer, never a plate ID.
    Recorded(u64),
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColumnProvenance {
    pub formation: Formation,
    pub latest_thermal_reset: Option<u64>,
}
impl ColumnProvenance {
    fn unknown(kind: CrustKindCodeV0) -> Self {
        Self {
            formation: if kind == CrustKindCodeV0::Oceanic {
                Formation::UnknownInherited
            } else {
                Formation::NotOceanic
            },
            latest_thermal_reset: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum EventKind {
    Birth {
        origin: BirthOrigin,
        created_m3: f64,
    },
    ThermalReset {
        formation: Formation,
        previous_thermal_birth_my: f64,
        effective_thermal_birth_my: f64,
        /// Column volume affected, with zero material creation in this operation.
        affected_m3: f64,
    },
}
#[derive(Debug, Clone, PartialEq)]
pub struct FormationEvent {
    pub id: u64,
    pub interval: TimeInterval,
    pub cell: u32,
    pub plate_at_event: u32,
    pub area_m2: f64,
    pub kind: EventKind,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct StepReceipt {
    pub provisional_gap_cells: usize,
    pub cancelled_gap_cells: usize,
    pub accepted_gap_cells: usize,
    pub projection_birth_cells: usize,
    /// Sum of accepted new ocean columns, excluding copied projection fills,
    /// arc/plume additions and stretched continental rifts.
    pub created_m3: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FormationAgeRecorderV1 {
    nlat: usize,
    nlon: usize,
    radius_m: f64,
    time_my: f64,
    started_at_my: f64,
    event_limit: usize,
    columns: Vec<ColumnProvenance>,
    events: Vec<FormationEvent>,
    last_step: StepReceipt,
}
impl FormationAgeRecorderV1 {
    /// Attach to a native state without inventing its missing history. An old
    /// checkpoint can only start an explicitly partial archive at this time.
    pub fn inherited(
        grid: &GridSpec,
        state: &CrustState,
        event_limit: usize,
    ) -> Result<Self, String> {
        if state.cells() != grid.cells()
            || state.kind.len() != grid.cells()
            || !state.t_my.is_finite()
            || event_limit == 0
            || grid.cells() > u32::MAX as usize
        {
            return Err("formation-age: invalid initial state, grid or event limit".into());
        }
        Ok(Self {
            nlat: grid.nlat(),
            nlon: grid.nlon(),
            radius_m: grid.radius_m(),
            time_my: state.t_my,
            started_at_my: state.t_my,
            event_limit,
            columns: state
                .kind
                .iter()
                .map(|&k| ColumnProvenance::unknown(k))
                .collect(),
            events: Vec::new(),
            last_step: StepReceipt::default(),
        })
    }
    pub fn columns(&self) -> &[ColumnProvenance] {
        &self.columns
    }
    pub fn events(&self) -> &[FormationEvent] {
        &self.events
    }
    pub fn last_step(&self) -> &StepReceipt {
        &self.last_step
    }
    pub fn time_my(&self) -> f64 {
        self.time_my
    }

    /// Persistent representation, distinct from the earlier diagnostic export.
    /// IDs, capture interval and the last native receipt are retained verbatim.
    pub fn checkpoint_v1(&self) -> Value {
        json!({"schema":"native-formation-age-checkpoint-v1", "recorder":self.to_json(),
            "last_step": {"provisional_gap_cells":self.last_step.provisional_gap_cells,
                "cancelled_gap_cells":self.last_step.cancelled_gap_cells,
                "accepted_gap_cells":self.last_step.accepted_gap_cells,
                "projection_birth_cells":self.last_step.projection_birth_cells,
                "created_m3":self.last_step.created_m3}})
    }

    pub fn from_checkpoint_v1(
        value: &Value,
        grid: &GridSpec,
        state: &CrustState,
        expected_limit: usize,
    ) -> Result<Self, String> {
        if value["schema"] != "native-formation-age-checkpoint-v1" {
            return Err("formation-age: unsupported checkpoint schema".into());
        }
        let bad = || "formation-age: malformed checkpoint value".to_string();
        let number = |v: &Value| v.as_f64().filter(|x| x.is_finite()).ok_or_else(bad);
        let integer = |v: &Value| v.as_u64().ok_or_else(bad);
        let size = |v: &Value| usize::try_from(integer(v)?).map_err(|_| bad());
        let formation = |v: &Value| match v["kind"].as_str() {
            Some("unknown-inherited") => Ok(Formation::UnknownInherited),
            Some("not-oceanic") => Ok(Formation::NotOceanic),
            Some("recorded") => Ok(Formation::Recorded(integer(&v["event_id"])?)),
            _ => Err(bad()),
        };
        let r = &value["recorder"];
        if r["schema"] != VERSION { return Err(bad()); }
        let raw_events = r["events"].as_array().ok_or_else(bad)?;
        if raw_events.len() > expected_limit { return Err(bad()); }
        let mut events = Vec::with_capacity(raw_events.len());
        for e in raw_events {
            let k = &e["event"];
            let kind = match k["kind"].as_str() {
                Some("birth") => EventKind::Birth {
                    origin: match k["origin"].as_str() {
                        Some("kinematic-gap") => BirthOrigin::KinematicGap,
                        Some("projection-initialization") => BirthOrigin::ProjectionInitialization,
                        _ => return Err(bad()),
                    },
                    created_m3: number(&k["created_m3"])? },
                Some("thermal-reset") => {
                    if number(&k["created_m3"])? != 0.0 { return Err(bad()); }
                    EventKind::ThermalReset {
                        formation: formation(&k["formation"])?,
                        previous_thermal_birth_my: number(&k["previous_thermal_birth_my"])?,
                        effective_thermal_birth_my: number(&k["effective_thermal_birth_my"])?,
                        affected_m3: number(&k["affected_m3"])? }
                }
                _ => return Err(bad()),
            };
            events.push(FormationEvent { id: integer(&e["id"])?,
                interval: TimeInterval::new(number(&e["start_my"])?, number(&e["end_my"])?)?,
                cell: u32::try_from(integer(&e["cell"])?).map_err(|_| bad())?,
                plate_at_event: u32::try_from(integer(&e["plate_at_event"])?).map_err(|_| bad())?,
                area_m2: number(&e["area_m2"])?, kind });
        }
        let raw_columns = r["columns"].as_array().ok_or_else(bad)?;
        if raw_columns.len() != grid.cells() { return Err(bad()); }
        let mut columns = Vec::with_capacity(raw_columns.len());
        for c in raw_columns {
            let reset = c.get("latest_thermal_reset").ok_or_else(bad)?;
            columns.push(ColumnProvenance { formation: formation(&c["formation"])?,
                latest_thermal_reset: if reset.is_null() { None } else { Some(integer(reset)?) } });
        }
        let step = &value["last_step"];
        let recorder = Self {
            nlat: size(&r["nlat"])?, nlon: size(&r["nlon"])?, radius_m: number(&r["radius_m"])?,
            time_my: number(&r["time_my"])?, started_at_my: number(&r["started_at_my"])?,
            event_limit: size(&r["event_limit"])?, columns, events,
            last_step: StepReceipt { provisional_gap_cells: size(&step["provisional_gap_cells"])?,
                cancelled_gap_cells: size(&step["cancelled_gap_cells"])?,
                accepted_gap_cells: size(&step["accepted_gap_cells"])?,
                projection_birth_cells: size(&step["projection_birth_cells"])?,
                created_m3: number(&step["created_m3"])? },
        };
        recorder.validate_checkpoint(grid, state, expected_limit)?;
        Ok(recorder)
    }

    /// This does not reconstruct history from thermal age. The supplied native
    /// state is used only to verify the archive's current grid, time and kinds.
    pub fn validate_checkpoint(
        &self,
        grid: &GridSpec,
        state: &CrustState,
        expected_limit: usize,
    ) -> Result<(), String> {
        let bad = || "formation-age: invalid checkpoint grid, times, IDs, references or receipts".to_string();
        if self.nlat != grid.nlat() || self.nlon != grid.nlon()
            || self.radius_m != grid.radius_m() || self.time_my != state.t_my
            || !self.time_my.is_finite() || !self.started_at_my.is_finite()
            || self.started_at_my > self.time_my || self.event_limit != expected_limit
            || self.event_limit == 0 || self.events.len() > self.event_limit
            || self.columns.len() != grid.cells() || state.kind.len() != grid.cells()
            || state.cells() != grid.cells()
        { return Err(bad()); }
        let formation_valid = |f: Formation, before: usize| match f {
            Formation::UnknownInherited => true,
            Formation::NotOceanic => false,
            Formation::Recorded(id) => id > 0 && id <= before as u64
                && matches!(self.events[id as usize - 1].kind, EventKind::Birth { .. }),
        };
        let mut previous = (self.started_at_my, self.started_at_my);
        let mut gaps = 0;
        let mut projections = 0;
        let mut created = 0.0;
        for (index, event) in self.events.iter().enumerate() {
            let t = event.interval;
            if event.id != index as u64 + 1 || event.cell as usize >= grid.cells()
                || TimeInterval::new(t.start_my, t.end_my).is_err()
                || t.start_my < self.started_at_my || t.end_my > self.time_my
                || (t.end_my, t.start_my) < previous
                || !event.area_m2.is_finite() || event.area_m2 <= 0.0
            { return Err(bad()); }
            let area = grid.cell_area_m2(event.cell as usize / grid.nlon()).map_err(|e| e.to_string())?;
            if (event.area_m2 - area).abs() > area * 1e-12 { return Err(bad()); }
            previous = (t.end_my, t.start_my);
            match event.kind {
                EventKind::Birth { origin, created_m3 } => {
                    if !created_m3.is_finite() || created_m3 <= 0.0 { return Err(bad()); }
                    if t.end_my == self.time_my {
                        created += created_m3;
                        match origin {
                            BirthOrigin::KinematicGap => gaps += 1,
                            BirthOrigin::ProjectionInitialization => projections += 1,
                        }
                    }
                }
                EventKind::ThermalReset { formation, previous_thermal_birth_my: old,
                    effective_thermal_birth_my: new, affected_m3 } => {
                    if !formation_valid(formation, index) || !old.is_finite() || !new.is_finite()
                        || new <= old || new > t.end_my || !affected_m3.is_finite() || affected_m3 <= 0.0
                    { return Err(bad()); }
                }
            }
        }
        for (cell, column) in self.columns.iter().enumerate() {
            if state.kind[cell] != CrustKindCodeV0::Oceanic {
                if column.formation != Formation::NotOceanic || column.latest_thermal_reset.is_some() {
                    return Err(bad());
                }
            } else if !formation_valid(column.formation, self.events.len()) {
                return Err(bad());
            }
            if let Some(id) = column.latest_thermal_reset {
                let event = id.checked_sub(1).and_then(|i| usize::try_from(i).ok())
                    .and_then(|i| self.events.get(i)).ok_or_else(bad)?;
                if !matches!(event.kind, EventKind::ThermalReset { formation, .. } if formation == column.formation) {
                    return Err(bad());
                }
            }
        }
        let r = &self.last_step;
        if r.provisional_gap_cells > grid.cells() || r.cancelled_gap_cells > r.provisional_gap_cells
            || r.accepted_gap_cells != gaps || r.projection_birth_cells != projections
            || r.accepted_gap_cells != r.provisional_gap_cells - r.cancelled_gap_cells
            || !r.created_m3.is_finite() || (r.created_m3 - created).abs() > created.max(1.0) * 1e-10
        { return Err(bad()); }
        Ok(())
    }

    /// Reserve worst-case one event/cell before the native mutation. Full
    /// archives fail visibly; nothing is silently dropped or auto-pruned.
    pub(crate) fn preflight(&self, grid: &GridSpec, time_my: f64) -> Result<(), String> {
        if self.nlat != grid.nlat()
            || self.nlon != grid.nlon()
            || self.radius_m != grid.radius_m()
            || self.time_my != time_my
        {
            return Err(
                "formation-age: grid/time mismatch; carry the sidecar with every native step"
                    .into(),
            );
        }
        if self
            .events
            .len()
            .checked_add(grid.cells())
            .map_or(true, |n| n > self.event_limit)
        {
            return Err("formation-age: event limit cannot accommodate this operation; no records dropped; explicit archive capacity exhausted".into());
        }
        Ok(())
    }

    pub(crate) fn accept_step(
        &mut self,
        before: &CrustState,
        after: &CrustState,
        source: &[u32],
        candidates: &[Option<BirthCandidate>],
        areas: &[f64],
        birth_thickness_m: f64,
    ) -> Result<(), String> {
        // A future native creation path must supply its own observation hook.
        // Fail before appending rather than silently losing accepted births.
        for (c, &src) in source.iter().enumerate() {
            if src == NO_SOURCE && after.kind[c] == CrustKindCodeV0::Oceanic {
                let volume = areas[c] * birth_thickness_m;
                if candidates[c].is_none() || !volume.is_finite() || volume <= 0.0 {
                    return Err(format!("formation-age: accepted birth at cell {c} lacks a candidate or finite positive volume"));
                }
            }
        }
        let interval = TimeInterval {
            start_my: before.t_my,
            end_my: after.t_my,
        };
        let mut receipt = StepReceipt::default();
        let mut columns = Vec::with_capacity(source.len());
        for (c, &src) in source.iter().enumerate() {
            let candidate = candidates[c];
            let accepted = src == NO_SOURCE && after.kind[c] == CrustKindCodeV0::Oceanic;
            if candidate.is_some_and(|b| b.provisional_gap) {
                receipt.provisional_gap_cells += 1;
                if !accepted || candidate.unwrap().origin != BirthOrigin::KinematicGap {
                    receipt.cancelled_gap_cells += 1;
                }
            }
            let mut column = ColumnProvenance::unknown(after.kind[c]);
            if accepted {
                // Every accepted candidate was validated before archive mutation.
                if let Some(birth) = candidate {
                    let origin = birth.origin;
                    let created_m3 = areas[c] * birth_thickness_m;
                    let id = self.events.len() as u64 + 1;
                    self.events.push(FormationEvent {
                        id,
                        interval,
                        cell: c as u32,
                        plate_at_event: after.plate[c].0,
                        area_m2: areas[c],
                        kind: EventKind::Birth { origin, created_m3 },
                    });
                    column.formation = Formation::Recorded(id);
                    receipt.created_m3 += created_m3;
                    match origin {
                        BirthOrigin::KinematicGap => receipt.accepted_gap_cells += 1,
                        BirthOrigin::ProjectionInitialization => {
                            receipt.projection_birth_cells += 1
                        }
                    }
                }
            } else if src != NO_SOURCE
                && after.kind[c] == CrustKindCodeV0::Oceanic
                && before.kind[src as usize] == CrustKindCodeV0::Oceanic
            {
                column = self.columns[src as usize];
            }
            columns.push(column);
        }
        self.columns = columns;
        self.time_my = after.t_my;
        self.last_step = receipt;
        Ok(())
    }

    pub(crate) fn thermal_resets(
        &mut self,
        state: &CrustState,
        previous_birth: &[f64],
        areas: &[f64],
        interval: TimeInterval,
    ) {
        for c in 0..self.columns.len() {
            if state.kind[c] != CrustKindCodeV0::Oceanic {
                // A preceding island-arc conversion can change kind outside
                // the carry step. Never label that column as seafloor basalt.
                self.columns[c] = ColumnProvenance::unknown(state.kind[c]);
                continue;
            }
            if state.birth_t_my[c] > previous_birth[c] {
                let id = self.events.len() as u64 + 1;
                self.events.push(FormationEvent {
                    id,
                    interval,
                    cell: c as u32,
                    plate_at_event: state.plate[c].0,
                    area_m2: areas[c],
                    kind: EventKind::ThermalReset {
                        formation: self.columns[c].formation,
                        previous_thermal_birth_my: previous_birth[c],
                        effective_thermal_birth_my: state.birth_t_my[c],
                        affected_m3: areas[c] * state.thickness_m[c],
                    },
                });
                self.columns[c].latest_thermal_reset = Some(id);
            }
        }
    }

    /// Diagnostic/export value, not an admitted product or a checkpoint reader.
    /// Consumers must preserve the run/grid namespace and the full archive.
    pub fn to_json(&self) -> Value {
        let formation = |f: Formation| match f {
            Formation::UnknownInherited => json!({"kind":"unknown-inherited"}),
            Formation::NotOceanic => json!({"kind":"not-oceanic"}),
            Formation::Recorded(id) => json!({"kind":"recorded", "event_id":id}),
        };
        let events: Vec<_> = self.events.iter().map(|e| {
            let kind = match e.kind {
                EventKind::Birth { origin, created_m3 } => json!({"kind":"birth",
                    "origin":origin.label(), "created_m3":created_m3}),
                EventKind::ThermalReset { formation: f, previous_thermal_birth_my,
                    effective_thermal_birth_my, affected_m3 } => json!({"kind":"thermal-reset",
                    "formation":formation(f), "previous_thermal_birth_my":previous_thermal_birth_my,
                    "effective_thermal_birth_my":effective_thermal_birth_my,
                    "affected_m3":affected_m3, "created_m3":0.0}),
            };
            json!({"id":e.id, "start_my":e.interval.start_my, "end_my":e.interval.end_my,
                "cell":e.cell, "plate_at_event":e.plate_at_event, "area_m2":e.area_m2, "event":kind})
        }).collect();
        let columns: Vec<_> = self.columns.iter().map(|c| json!({
            "formation":formation(c.formation), "latest_thermal_reset":c.latest_thermal_reset,
        })).collect();
        json!({"schema":VERSION, "nlat":self.nlat, "nlon":self.nlon, "radius_m":self.radius_m,
            "time_my":self.time_my, "started_at_my":self.started_at_my, "event_limit":self.event_limit,
            "identity_scope":"run-local Eulerian origin cohorts, not conserved ocean parcels",
            "events":events, "columns":columns})
    }
}
