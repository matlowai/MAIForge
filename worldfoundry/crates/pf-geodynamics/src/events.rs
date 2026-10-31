//! Append-only tectonic event ledger with chunked deterministic sidecar
//! bytes (brief K7; acceptance: stable event IDs, monotonic allocation,
//! declared limit refuses LOUDLY with context, chunk writer bytes
//! deterministic).
//!
//! The GENERIC artifact layer binds sidecar bytes (path/length/digest) via
//! `SidecarDescriptor`; this module only produces the BYTES and the
//! domain-side chunk index. No second digest notion lives here (one
//! verification universe, preflight steering 2026-09-07).

use pf_contracts::domain::{EventId, EventKindV0, ParcelId, PlateId, TectonicEventV0};
use pf_contracts::spatial::ElementId;

/// Loud event-ledger errors.
#[derive(Debug, Clone, PartialEq)]
pub enum EventsError {
    LimitExceeded {
        limit: u64,
        attempted_total: u64,
        kind: EventKindV0,
    },
    NonMonotonicTime {
        last: f64,
        got: f64,
    },
    NonFinite {
        what: &'static str,
        value: f64,
    },
    BadRange {
        from: u64,
        to: u64,
        events: u64,
    },
    Serialization(String),
}

impl std::fmt::Display for EventsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EventsError::LimitExceeded { limit, attempted_total, kind } => write!(
                f,
                "event limit {limit} exceeded (would be {attempted_total} events; last kind {kind:?}); \
                 refusing loudly — raise the declared limit or chunk more aggressively, never truncate"
            ),
            EventsError::NonMonotonicTime { last, got } => {
                write!(f, "event at {got:?} My predates ledger head {last:?} My")
            }
            EventsError::NonFinite { what, value } => write!(f, "nonfinite {what}: {value:?}"),
            EventsError::BadRange { from, to, events } => {
                write!(f, "chunk range [{from}, {to}] invalid for {events} events")
            }
            EventsError::Serialization(why) => write!(f, "serialization: {why}"),
        }
    }
}

impl std::error::Error for EventsError {}

/// Declared default limit (loud refusal, never truncation — K7). Callers
/// may configure lower/higher explicitly.
pub const DEFAULT_MAX_EVENTS: u64 = 5_000_000;

#[derive(Debug, Clone, PartialEq)]
pub struct EventLedger {
    events: Vec<TectonicEventV0>,
    max_events: u64,
    last_t_my: f64,
}

impl EventLedger {
    pub fn new(max_events: u64) -> Self {
        Self {
            events: Vec::new(),
            max_events,
            last_t_my: 0.0,
        }
    }

    pub fn len(&self) -> u64 {
        self.events.len() as u64
    }

    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    pub fn events(&self) -> &[TectonicEventV0] {
        &self.events
    }

    pub fn next_event_id(&self) -> EventId {
        EventId(self.events.len() as u64)
    }

    /// Append one event: id allocated monotonically, time must not precede
    /// the ledger head, limit refusal carries limit + attempted total.
    pub fn push(
        &mut self,
        t_my: f64,
        kind: EventKindV0,
        plates: Vec<PlateId>,
        footprint_elements: Vec<ElementId>,
        material_lineage: Vec<ParcelId>,
        note: String,
    ) -> Result<EventId, EventsError> {
        if !t_my.is_finite() || t_my < 0.0 {
            return Err(EventsError::NonFinite {
                what: "t_my",
                value: t_my,
            });
        }
        if t_my < self.last_t_my {
            return Err(EventsError::NonMonotonicTime {
                last: self.last_t_my,
                got: t_my,
            });
        }
        if self.events.len() as u64 >= self.max_events {
            return Err(EventsError::LimitExceeded {
                limit: self.max_events,
                attempted_total: self.events.len() as u64 + 1,
                kind,
            });
        }
        let id = self.next_event_id();
        self.events.push(TectonicEventV0 {
            event_id: id,
            t_my,
            kind,
            plates,
            footprint_elements,
            material_lineage,
            note,
        });
        self.last_t_my = t_my;
        Ok(id)
    }

    /// Deterministic sidecar bytes for events `[from, to)` (by INDEX,
    /// which equals event id order): canonical JSON array. Two calls with
    /// the same ledger state produce identical bytes. The generic artifact
    /// layer hashes these bytes into the manifest sidecar digest.
    pub fn chunk_bytes(&self, from: u64, to: u64) -> Result<Vec<u8>, EventsError> {
        if from > to || to > self.events.len() as u64 {
            return Err(EventsError::BadRange {
                from,
                to,
                events: self.events.len() as u64,
            });
        }
        serde_json::to_vec(&self.events[from as usize..to as usize])
            .map_err(|e| EventsError::Serialization(e.to_string()))
    }

    /// Chunk index plan: split into chunks of at most `target` records,
    /// returning (start, end) ranges in order.
    pub fn chunk_plan(&self, target: u64) -> Vec<(u64, u64)> {
        let total = self.events.len() as u64;
        let mut plan = Vec::new();
        let mut start = 0u64;
        while start < total {
            let end = (start + target).min(total);
            plan.push((start, end));
            start = end;
        }
        plan
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn push_ev(l: &mut EventLedger, t: f64, kind: EventKindV0) -> EventId {
        l.push(
            t,
            kind,
            vec![PlateId(1)],
            vec![ElementId(3)],
            vec![ParcelId(9)],
            "test".into(),
        )
        .unwrap()
    }

    #[test]
    fn ids_are_monotonic_and_unique() {
        let mut l = EventLedger::new(DEFAULT_MAX_EVENTS);
        let a = push_ev(&mut l, 1.0, EventKindV0::RidgeBirth);
        let b = push_ev(&mut l, 2.0, EventKindV0::HotspotCrossing);
        let c = push_ev(&mut l, 2.0, EventKindV0::BoundaryReassignment);
        assert!(a < b && b < c, "{a:?} {b:?} {c:?}");
        let ids: std::collections::BTreeSet<_> = l.events().iter().map(|e| e.event_id).collect();
        assert_eq!(ids.len(), l.len() as usize);
    }

    #[test]
    fn limit_refusal_is_loud_and_specific() {
        let mut l = EventLedger::new(3);
        push_ev(&mut l, 1.0, EventKindV0::RidgeBirth);
        push_ev(&mut l, 2.0, EventKindV0::RidgeBirth);
        push_ev(&mut l, 3.0, EventKindV0::RidgeBirth);
        let err = l
            .push(
                4.0,
                EventKindV0::SubductionRetirement,
                vec![],
                vec![],
                vec![],
                "x".into(),
            )
            .unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("limit 3"), "{msg}");
        assert!(msg.contains("would be 4"), "{msg}");
        assert!(msg.contains("never truncate"), "{msg}");
        // The refusal did not append anything.
        assert_eq!(l.len(), 3);
    }

    #[test]
    fn time_must_not_run_backward() {
        let mut l = EventLedger::new(10);
        push_ev(&mut l, 5.0, EventKindV0::RidgeBirth);
        let err = l
            .push(
                4.9,
                EventKindV0::RidgeBirth,
                vec![],
                vec![],
                vec![],
                "x".into(),
            )
            .unwrap_err();
        assert!(matches!(err, EventsError::NonMonotonicTime { .. }), "{err}");
        // Equal time is fine (many events per step).
        push_ev(&mut l, 5.0, EventKindV0::TransformTear);
    }

    #[test]
    fn chunk_bytes_are_deterministic_and_range_checked() {
        let mut l = EventLedger::new(100);
        for k in 0..25u64 {
            push_ev(&mut l, k as f64, EventKindV0::RidgeBirth);
        }
        let a = l.chunk_bytes(0, 10).unwrap();
        let b = l.chunk_bytes(0, 10).unwrap();
        assert_eq!(a, b, "same state, same bytes");
        assert_ne!(a, l.chunk_bytes(0, 11).unwrap());
        assert!(l.chunk_bytes(10, 10).unwrap().is_empty() || l.chunk_bytes(10, 10).is_ok());
        let err = l.chunk_bytes(20, 26).unwrap_err();
        assert!(matches!(err, EventsError::BadRange { .. }), "{err}");
        // Chunk plan covers exactly, in order, within target sizes.
        let plan = l.chunk_plan(10);
        assert_eq!(plan, vec![(0, 10), (10, 20), (20, 25)]);
        // Bytes round-trip as the same records.
        let back: Vec<TectonicEventV0> = serde_json::from_slice(&a).unwrap();
        assert_eq!(back.len(), 10);
        assert_eq!(back[0].event_id, EventId(0));
    }
}
