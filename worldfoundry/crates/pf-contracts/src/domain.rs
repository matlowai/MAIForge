//! Durable domain schemas for geodynamics, hydrology and sediment transfer
//! (overnight preflight 2026-09-07; brief §K/§S5a/§S4, spec §54.4/§54.8/§57.5).
//!
//! THESE ARE CONTRACT SCHEMAS, NOT BACKEND STRUCTS. They are the durable
//! typed surface the overnight agents program against and that artifacts
//! persist; solver-internal state (remap bins, routing queues) stays in
//! its backend crate and never appears here.
//!
//! FOUR STRUCTURAL RULES (operator steering, 2026-09-07, pre-freeze):
//!
//! 1. IDENTITY CONSERVATION != MATERIAL CONSERVATION. Parcel counts close
//!    EXACTLY (integers): initial + births - retirements = survivors. Rock
//!    quantity closes only within DECLARED numerical error, in a physical
//!    unit (B1: material-area-equivalent m^2): initial material + created
//!    minus retired/exported = surviving material within declared_error.
//!    Two gates, two types, never merged.
//!
//! 2. PARCEL != CELL, AND THE PROJECTION IS WEIGHTED. The parcel ledger is
//!    the authoritative material identity/history (§55 lineage seam) and
//!    carries no spatial sampling. The derived projection is a set of
//!    `CellParcelContribution { element, parcel, weight }` with per-parcel
//!    weights summing to ~1. B1 nearest-cell emits ONE contribution at
//!    weight = 1.0; a later conservative/PIC remapper may legitimately
//!    split one parcel 0.63/0.37 across two cells without changing the
//!    artifact language. Nothing here forbids that.
//!
//! 3. UNBOUNDED LEDGERS NEVER LIVE IN artifact.json. Events, parcels and
//!    occupancy contributions are bulk data: they persist as GENERIC
//!    manifest-bound sidecars (`SemanticManifest::sidecars`) whose bytes
//!    are digest-verified by the generic artifact layer exactly like
//!    field blobs — ONE verification universe. The domain summary carries
//!    only DESCRIPTORS (counts, sidecar logical-name refs, time ranges);
//!    path/length/digest truth lives solely in the manifest. Small
//!    bounded tables (plates, boundaries, basins, rivers) stay inline.
//!    Chunk refusal at declared limits is a loud runtime error, never
//!    truncation (K7).
//!
//! 4. NO DUPLICATED TRUTH. Hydrology does not restate the terrain digest:
//!    its binding to the published terrain is the semantic manifest's
//!    dependency input (one truth, cross-validated at verify).
//!
//! All persisted counts are fixed-width (u64/u32/u8/f64, never usize).
//! Fail-loud validation names stage and field; no silent clamping.

use crate::artifact::SidecarDescriptor;
use crate::spatial::ElementId;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

#[derive(Debug, thiserror::Error)]
pub enum DomainError {
    #[error("domain validation failed at {stage}: {field}: {why}")]
    Violation {
        stage: &'static str,
        field: String,
        why: String,
    },
}

fn violate(stage: &'static str, field: &str, why: impl Into<String>) -> DomainError {
    DomainError::Violation {
        stage,
        field: field.into(),
        why: why.into(),
    }
}

/// Tolerance for declared fraction/weight splits (grain classes, phases,
/// provenance mixtures, projection weights). Loose enough to survive
/// canonical float encoding, tight enough to catch real drift.
const FRACTION_SUM_TOL: f64 = 1.0e-9;

fn closes_to_one(sum: f64, tol: f64) -> bool {
    (sum - 1.0).abs() <= tol
}

// ---------------------------------------------------------------------------
// Shared identifiers (fixed-width newtypes; stable across epochs)

/// Material identity of one crust parcel (PARCEL, never CELL).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ParcelId(pub u64);

/// Plate identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct PlateId(pub u32);

/// Stable tectonic event identity (K7: stable IDs, never renumbered).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct EventId(pub u64);

/// Drainage basin identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct BasinId(pub u64);

/// River reach identity (network edge between junctions/outlets).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ReachId(pub u64);

// ---------------------------------------------------------------------------
// Geodynamics: bounded inline tables

/// Euler rotation pole for one plate (K2). Degrees + rate; sign encodes
/// direction. The durable record of HOW a plate moved, for provenance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EulerPoleV0 {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub rate_deg_per_my: f64,
}

/// How a plate left the active set. Extensible; v0 has subduction only
/// because that is the only plate retirement K2 declares.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PlateRetirementV0 {
    Subducted,
    /// SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13; additive): the plate's
    /// cells were folded into another plate by a reorganization law (merge,
    /// fragment fold, consolidation into a group-mate) and it no longer
    /// claims cells or parcels. Emitted only by `plate-reorg-scripted-v2`.
    Folded,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PlateRetirementDetailV0 {
    pub t_my: f64,
    pub cause: PlateRetirementV0,
}

/// Bounded table (plate_count <= 32 per params): safe inline in the
/// domain summary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlateRecordV0 {
    pub plate_id: PlateId,
    /// Lineage: plate this one split from, if any.
    pub parent_plate: Option<PlateId>,
    pub birth_t_my: f64,
    /// Present while the plate moves; `None` only after retirement.
    pub euler_pole: Option<EulerPoleV0>,
    pub retirement: Option<PlateRetirementDetailV0>,
}

/// Boundary classification (K3/K5). Convergent splits into subduction vs
/// collision because their rules differ downstream.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoundaryKindV0 {
    Ridge,
    Subduction,
    Collision,
    Transform,
}

/// Stable boundary segment record (K3: stable IDs); bounded by plate
/// pairs, safe inline. The per-cell footprint is a field blob keyed by
/// `boundary_id`; the schema carries the count for manifest accounting.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoundaryRecordV0 {
    pub boundary_id: u64,
    pub kind: BoundaryKindV0,
    pub plate_a: PlateId,
    pub plate_b: PlateId,
    pub birth_t_my: f64,
    pub last_active_t_my: f64,
    pub footprint_element_count: u64,
}

// ---------------------------------------------------------------------------
// Geodynamics: bulk ledgers (chunk sidecars + descriptors)

/// Tectonic event kinds (K2/K5/K6/K7). Backends extend only via an
/// additive v1 (A2.1 policy).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EventKindV0 {
    /// New crust born at a ridge (parcel birth, thermal age zero).
    RidgeBirth,
    /// Parcel consumed at a subduction boundary.
    SubductionRetirement,
    /// Boundary swept a cell: parcel changed plate (lineage event).
    BoundaryReassignment,
    /// Plate crossed a fixed hotspot (K6; mints age-progressive chains).
    HotspotCrossing,
    /// Collision began between two plates (K5 wide-collision rule).
    CollisionOnset,
    /// Transform motion tore a field (K5 torn-field marker).
    TransformTear,
}

/// One ledger entry (K7): `{id, t_my, kind, plates, footprint cells,
/// material lineage, note}`. Append-only; persisted in deterministic
/// chunks; a declared limit refusal is loud, never truncation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TectonicEventV0 {
    pub event_id: EventId,
    pub t_my: f64,
    pub kind: EventKindV0,
    pub plates: Vec<PlateId>,
    /// Footprint cells as element ids of the surface domain (may be empty
    /// for purely ledger-kind events).
    pub footprint_elements: Vec<ElementId>,
    /// Material lineage touched by this event.
    pub material_lineage: Vec<ParcelId>,
    pub note: String,
}

/// How a parcel came to exist. Initial parcels cover the seeded sphere;
/// everything else is born by the physics the ledger records.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ParcelBirthV0 {
    /// Present at t=0, oceanic.
    InitialOceanic,
    /// Present at t=0, continental seed.
    InitialContinental,
    /// Created at a spreading center (K2; thermal + crust age zero).
    RidgeCreated,
    /// Created by hotspot volcanism (K6 island chains).
    HotspotCreated,
    /// Accreted in a collision (transfer to the overriding plate).
    Accreted,
}

/// Why a parcel left the ledger. v0: subduction only (K2).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ParcelRetirementCauseV0 {
    Subducted,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ParcelRetirementV0Detail {
    pub t_my: f64,
    pub cause: ParcelRetirementCauseV0,
    /// Event that retired this parcel (ledger cross-reference).
    pub event: EventId,
}

/// THE authoritative material record (PARCEL = rock identity/history;
/// §55 material lineage). Deliberately carries no spatial sampling: where
/// the material is currently sampled on the grid is the derived weighted
/// projection. `material_area_equivalent_m2` is a CONSERVED LAGRANGIAN
/// material weight, not the parcel's current geometric footprint: later
/// collision shortening/thickening may reduce exposed area without
/// annihilating material, and conservative remaps partition THIS quantity
/// (projection weights sum to it, not to a re-measured footprint).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CrustParcelRecordV0 {
    pub parcel_id: ParcelId,
    pub birth_t_my: f64,
    pub birth: ParcelBirthV0,
    /// Lineage parent for Accreted parcels; `None` for initial and
    /// ridge/hotspot births (new material, no parent rock).
    pub lineage_parent: Option<ParcelId>,
    /// Owning plate while alive.
    pub plate: PlateId,
    /// Conserved Lagrangian material weight of this parcel, in the
    /// ledger's declared unit (B1: material-area-equivalent, m^2). NOT the
    /// current geometric footprint — shortening/thickening changes
    /// footprint, not this. Volume/mass upgrades arrive additively (v1)
    /// without changing the conservation gates.
    pub material_area_equivalent_m2: f64,
    pub retirement: Option<ParcelRetirementV0Detail>,
}

/// GATE 1 — IDENTITY conservation (exact, integer domain):
/// `initial + births - retirements = survivors`. Counts never carry
/// numerical error; a mismatch is a ledger bug, full stop.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ParcelIdentityConservationV0 {
    pub t_my: f64,
    pub initial_parcels: u64,
    pub births: u64,
    pub retirements: u64,
    pub survivors: u64,
}

impl ParcelIdentityConservationV0 {
    pub fn identity_holds(&self) -> bool {
        let lhs = self.initial_parcels as i64 + self.births as i64 - self.retirements as i64;
        lhs == self.survivors as i64
    }
}

/// GATE 2 — MATERIAL conservation (physical quantity, declared error):
/// `initial + created - retired_or_exported = surviving +- declared_error`.
/// B1 unit: represented surface area m^2.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MaterialConservationV0 {
    pub t_my: f64,
    /// Units string from the pf-contracts units table (e.g. "m2").
    pub quantity: String,
    pub initial_material: f64,
    pub created_material: f64,
    pub retired_exported_material: f64,
    pub surviving_material: f64,
    /// Declared numerical error (summation/remap). Must be >= 0; a large
    /// value is a smell the reviewer hunts, not something to hide.
    pub declared_error: f64,
}

impl MaterialConservationV0 {
    pub fn closes(&self) -> bool {
        let residual = self.initial_material + self.created_material
            - self.retired_exported_material
            - self.surviving_material;
        residual.is_finite() && residual.abs() <= self.declared_error
    }
}

/// The loaded parcel ledger (chunk records + both gates). Writers validate
/// before sealing chunks; verify reloads and revalidates. Never embedded
/// wholesale in artifact.json — only its descriptor is.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParcelLedgerV0 {
    pub parcels: Vec<CrustParcelRecordV0>,
    pub identity: ParcelIdentityConservationV0,
    pub material: MaterialConservationV0,
}

/// Reference to one chunk of a bulk ledger, persisted as a generic
/// manifest-bound sidecar. Path, bytes, length and digest live in
/// `SemanticManifest::sidecars`; this ref carries only domain accounting.
/// Deterministic sidecar writers produce deterministic bytes; the GENERIC
/// artifact layer hashes those bytes into `SidecarDescriptor::digest` —
/// there is deliberately NO domain-level "chunk digest" notion (one
/// verification universe; preflight steering 2026-09-07).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LedgerChunkRefV0 {
    /// Logical name of a `SidecarDescriptor` in the same artifact's
    /// semantic manifest (e.g. "events.chunk-000").
    pub sidecar: String,
    /// Records in this chunk; chunks are never empty.
    pub record_count: u64,
}

impl LedgerChunkRefV0 {
    fn validate(&self, stage: &'static str) -> Result<(), DomainError> {
        if self.sidecar.is_empty() {
            return Err(violate(stage, "sidecar", "empty sidecar logical name"));
        }
        if self.record_count == 0 {
            return Err(violate(stage, "record_count", "empty chunk"));
        }
        Ok(())
    }
}

/// Descriptor for the event ledger: counts, chunk index, time range.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EventLedgerDescriptorV0 {
    pub schema: String,
    pub total_events: u64,
    /// Inclusive time range of the ledgered events, when any exist.
    pub t_start_my: Option<f64>,
    pub t_end_my: Option<f64>,
    pub chunks: Vec<LedgerChunkRefV0>,
}

/// Descriptor for the parcel ledger: counts + chunk index. Conservation
/// summaries live beside it in the domain (they describe the SAME loaded
/// records the chunks persist).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParcelLedgerDescriptorV0 {
    pub schema: String,
    pub total_parcels: u64,
    pub live_parcels: u64,
    pub chunks: Vec<LedgerChunkRefV0>,
}

/// Descriptor for the derived weighted occupancy projection (bulk data:
/// chunked like the ledgers, or a field-blob pair — the descriptor is the
/// contract either way).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProjectionDescriptorV0 {
    /// Instant the projection was taken; must equal the domain's t_my.
    pub projected_t_my: f64,
    pub contribution_count: u64,
    /// Distinct parcels carrying nonzero weight (<= live_parcels).
    pub projected_parcel_count: u64,
    pub chunks: Vec<LedgerChunkRefV0>,
}

// ---------------------------------------------------------------------------
// Geodynamics: the artifact.json domain summary

/// Relational summary published with `pf.geodynamics/v0` artifacts. Bulk
/// ledgers and the projection persist as chunk sidecars described by the
/// descriptors here; field blobs (crust_age, elevation, thermal, ...)
/// reference these ids.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GeodynamicsDomainV0 {
    pub t_my: f64,
    /// Bounded inline tables.
    pub plates: Vec<PlateRecordV0>,
    pub boundaries: Vec<BoundaryRecordV0>,
    /// Bulk ledgers + derived projection (descriptor + sidecars).
    pub event_ledger: EventLedgerDescriptorV0,
    pub parcel_ledger: ParcelLedgerDescriptorV0,
    pub occupancy: ProjectionDescriptorV0,
    /// The two conservation gates, recorded per committed epoch (and the
    /// final epoch's values here; per-step history is a sidecar table).
    pub identity_conservation: ParcelIdentityConservationV0,
    pub material_conservation: MaterialConservationV0,
}

// ---------------------------------------------------------------------------
// Hydrology domain (pf.hydrology/v0 relational schema)

/// Internal-drainage classification (A2/S3): classified, never erased.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DrainageClassV0 {
    Endorheic,
    OceanConnected,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BasinRecordV0 {
    pub basin_id: BasinId,
    pub classification: DrainageClassV0,
    /// Spill level for internal basins (depression filling, S2); the
    /// physical elevation field is never fill-adjusted.
    pub spill_level_m: Option<f64>,
    /// Lake outlet element for endorheic basins that overflow.
    pub outlet_element: Option<ElementId>,
    /// River reaches whose water terminates in this basin.
    pub inlet_reaches: Vec<ReachId>,
}

/// One river network edge (S4). Discharge proxy is rainfall-weighted
/// physical contributing area in m3/s-class units, resolution-independent.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RiverReachV0 {
    pub reach_id: ReachId,
    /// Downstream reach; `None` = terminates at an ocean or a classified
    /// lake basin (recorded via `outlet_basin`).
    pub downstream: Option<ReachId>,
    pub outlet_basin: Option<BasinId>,
    pub strahler_order: u8,
    pub mean_discharge_proxy_m3_s: f64,
}

/// Waterfall cause labels (§57.5 B1 set). Caprock/lava-edge labels arrive
/// with stratigraphy (additive v1). DETECTION labels — the detector's
/// thresholds are constants, never functions of `waterfall_wonder`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WaterfallCause {
    UpliftStep,
    DeformationBand,
    HangingValley,
}

/// One detected waterfall (S4/§57.5): a longitudinal-profile knickpoint
/// or perched tributary confluence above the declared constant thresholds,
/// labeled with its inferred cause. Detection, never decree.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterfallObservation {
    pub element: ElementId,
    /// Observed drop over the reach; strictly positive and finite.
    pub drop_m: f64,
    pub cause: WaterfallCause,
    /// Reach the fall sits on, when applicable.
    pub reach: Option<ReachId>,
}

/// Relational schema published with `pf.hydrology/v0`. Binds to exactly
/// one published terrain via the semantic manifest's dependency input
/// (ONE truth: no digest restated here — rule 4). Per-cell routing arrays
/// are field blobs; these tables are bounded by physical thresholds.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HydrologyDomainV0 {
    /// Epoch-local physical water bodies, distinct from terminal catchments.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub water_bodies: Vec<WaterBodyV0>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exchange_history: Option<ExchangeHistoryV0>,
    pub t_my: f64,
    pub basins: Vec<BasinRecordV0>,
    pub rivers: Vec<RiverReachV0>,
    pub waterfalls: Vec<WaterfallObservation>,
    /// Per-cell routing field blobs and their sentinel conventions
    /// (additive 2026-09-13, ADR 032). `None` for artifacts that publish
    /// only the relational tables.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_fields: Option<HydrologyCellFieldsV0>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WaterBodyKindV0 {
    Marine,
    Lake,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterTransferV0 {
    pub from_element: ElementId,
    pub to_element: ElementId,
    pub receiving_body: u32,
    /// fluvial_delivery, diffusive_delivery, lake_overflow or runoff_rate.
    pub kind: String,
    pub sediment_m3: f64,
    /// Steady proxy rate, not an integrated ocean-water inventory.
    pub runoff_m3_s: f64,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExchangeHistoryV0 {
    pub sidecar: String,
    pub epochs: u64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArcNormalSampleV0 {
    pub overlap_element: ElementId,
    pub winner: PlateId,
    pub loser: PlateId,
    pub input_interface: [ElementId; 2],
    pub source_epoch_my: f64,
    pub normal_rate_m_per_my: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterExchangeEpochV0 {
    pub arc_normals: Vec<ArcNormalSampleV0>,
    pub start_my: f64,
    pub end_my: f64,
    pub sediment_density_kg_m3: f64,
    pub ocean_seed: ElementId,
    pub ocean_area_m2: f64,
    pub water_bodies: Vec<WaterBodyV0>,
    pub transfers: Vec<WaterTransferV0>,
    /// Totals keyed by body and transfer kind; overflow is never new input.
    pub sediment_totals_m3: BTreeMap<String, f64>,
    pub crust_advection_residual_m3: f64,
    pub crust_advection_absolute_cumulative_m3: f64,
    pub wedge_in_m3: f64,
    pub wedge_stacked_m3: f64,
    pub wedge_delaminated_m3: f64,
    pub spread_bound_hits: u64,
    pub global_marine_layer_m3: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterBodyV0 {
    pub body_id: u32,
    pub kind: WaterBodyKindV0,
    pub seed_element: ElementId,
    pub water_surface_m: f64,
    pub area_m2: f64,
    /// Every directed routing exit, not an invented unique spill.
    pub outlet_edges: Vec<[ElementId; 2]>,
}

/// Names and sentinels of the per-cell routing fields a hydrology
/// artifact publishes on the terrain's surface domain: flow direction as
/// a receiver element id, basin membership, reach membership, the routing
/// (fill) surface, contributing area and the discharge proxy. Sentinels
/// are DECLARED here so a reader never infers "u32::MAX means terminal".
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HydrologyCellFieldsV0 {
    /// u32 field: receiver element id; `terminal_code` at ocean cells and
    /// endorheic lake floors.
    pub flow_to_field: String,
    pub terminal_code: u32,
    /// u32 field: basin id of the cell; `no_basin_code` where the routing
    /// assigned none (ocean cells that receive no land flow).
    pub basin_field: String,
    pub no_basin_code: u32,
    /// u32 field: reach id where the cell is a river cell; `no_reach_code`
    /// elsewhere.
    pub reach_field: String,
    pub no_reach_code: u32,
    /// f32 fields on the same domain.
    pub fill_level_field: String,
    pub contributing_area_field: String,
    pub discharge_field: String,
    /// Lake level convention: an endorheic basin's water surface is its
    /// `spill_level_m`; the physical elevation field is never fill-adjusted.
    pub lake_level_convention: String,
}

impl HydrologyCellFieldsV0 {
    pub fn v0() -> Self {
        Self {
            flow_to_field: "flow_to_element".into(),
            terminal_code: u32::MAX,
            basin_field: "basin_id".into(),
            no_basin_code: u32::MAX,
            reach_field: "reach_id".into(),
            no_reach_code: u32::MAX,
            fill_level_field: "fill_level_m".into(),
            contributing_area_field: "contributing_area_m2".into(),
            discharge_field: "discharge_proxy_m3_s".into(),
            lake_level_convention: "water_bodies records marine prescribed stage and lake fill-to-spill stages; fill_level_m is routing elevation on dry land; elevation_m remains physical".into(),
        }
    }
}

// ---------------------------------------------------------------------------
// Sediment transfer (S5a / §54.8) — provenance-preserving ledger

/// Coarse grain-class summary of transferred material (mass fractions).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GrainClassSplit {
    pub coarse_fraction: f64,
    pub sand_fraction: f64,
    pub silt_fraction: f64,
    pub clay_fraction: f64,
}

/// Particulate vs dissolved phase split (§54.8: the B1 erosion pass
/// PRESERVES this split so chemical/biogenic sediment routes can be added
/// without re-instrumenting provenance).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PhaseSplit {
    pub particulate_fraction: f64,
    pub dissolved_fraction: f64,
}

/// One provenance component: lineage label (parcel id or province name)
/// and its mass fraction. Mixtures close to 1.0.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProvenanceFraction {
    pub lineage: String,
    pub fraction: f64,
}

/// Transport system that moved the material (§54.1 loop).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TransportPathway {
    Fluvial,
    Marine,
    Eolian,
    Glacial,
}

/// Where deposition ended up. `destination_elements` required for sinks
/// inside the modeled surface; ocean/deep-sink export may leave it empty.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SedimentSinkV0 {
    Basin { basin: BasinId },
    Floodplain,
    Delta,
    OceanExport,
    DeepSink,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SedimentDepositionV0 {
    pub deposited_t_my: f64,
    pub sink: SedimentSinkV0,
    pub destination_elements: Vec<ElementId>,
}

/// One provenance-preserving sediment transfer record per epoch movement
/// (S5a): eroded mass never vanishes into a field name; it becomes this.
/// Bulk-ledger rule 3 applies at scale: epoch summaries inline, full
/// record streams chunked — the record schema is identical either way.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SedimentTransferRecord {
    pub record_id: u64,
    pub epoch_t_my: f64,
    /// Source region: elements (+ optional basin reference in the
    /// artifact's relational tables).
    pub source_region: Vec<ElementId>,
    /// Concrete parcel lineages the material came from.
    pub source_material_lineages: Vec<ParcelId>,
    /// Coarse provenance mixture (§54.2-style: "72% granite arc, ...").
    pub source_mixture: Vec<ProvenanceFraction>,
    /// Total mass moved by this record; strictly positive and finite.
    pub mass_kg: f64,
    pub grain_classes: GrainClassSplit,
    pub phase_split: PhaseSplit,
    pub transport_pathway: TransportPathway,
    /// `None` while still in transit at epoch end; the next epoch's record
    /// continues the story. Deposited material records where it landed.
    pub deposition: Option<SedimentDepositionV0>,
    /// Material exported out of the modeled surface system (§54.4 sink).
    pub exported_deep_sink: bool,
}

/// Per-epoch sediment mass balance (§54.4 gate: eroded = deposited +
/// exported +- declared error; per step AND cumulative).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SedimentEpochBalanceV0 {
    pub t_my: f64,
    pub eroded_kg: f64,
    pub deposited_kg: f64,
    pub exported_kg: f64,
    pub declared_error_kg: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SedimentMassBalanceV0 {
    pub epochs: Vec<SedimentEpochBalanceV0>,
    /// Declared cumulative closure; validated against the epoch rows.
    pub cumulative_eroded_kg: f64,
    pub cumulative_deposited_kg: f64,
    pub cumulative_exported_kg: f64,
    pub cumulative_declared_error_kg: f64,
}

// ---------------------------------------------------------------------------
// Validation (fail loud, name the field)

fn check_finite_positive(stage: &'static str, field: &str, v: f64) -> Result<(), DomainError> {
    if !v.is_finite() || v <= 0.0 {
        return Err(violate(
            stage,
            field,
            format!("expected finite > 0, found {v:?}"),
        ));
    }
    Ok(())
}

fn check_finite_nonneg(stage: &'static str, field: &str, v: f64) -> Result<(), DomainError> {
    if !v.is_finite() || v < 0.0 {
        return Err(violate(
            stage,
            field,
            format!("expected finite >= 0, found {v:?}"),
        ));
    }
    Ok(())
}

impl EulerPoleV0 {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "EulerPoleV0";
        if !self.lat_deg.is_finite() || !(-90.0..=90.0).contains(&self.lat_deg) {
            return Err(violate(
                stage,
                "lat_deg",
                format!("out of [-90, 90]: {:?}", self.lat_deg),
            ));
        }
        if !self.lon_deg.is_finite() || !(-180.0..=180.0).contains(&self.lon_deg) {
            return Err(violate(
                stage,
                "lon_deg",
                format!("out of [-180, 180]: {:?}", self.lon_deg),
            ));
        }
        if !self.rate_deg_per_my.is_finite() {
            return Err(violate(stage, "rate_deg_per_my", "nonfinite"));
        }
        Ok(())
    }
}

impl ParcelLedgerV0 {
    /// The loaded-ledger gate: BOTH conservation invariants recomputed
    /// from the records must match the declared summaries (identity
    /// exactly; material within its declared error, summed in canonical
    /// parcel_id order). Lineage must be complete and monotone.
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "ParcelLedgerV0";
        let mut by_id: HashMap<ParcelId, &CrustParcelRecordV0> = HashMap::new();
        for (i, p) in self.parcels.iter().enumerate() {
            if by_id.insert(p.parcel_id, p).is_some() {
                return Err(violate(
                    stage,
                    "parcels",
                    format!("duplicate parcel_id {:?} at index {i}", p.parcel_id),
                ));
            }
            check_finite_nonneg(stage, &format!("parcels[{i}].birth_t_my"), p.birth_t_my)?;
            check_finite_positive(
                stage,
                &format!("parcels[{i}].material_area_equivalent_m2"),
                p.material_area_equivalent_m2,
            )?;
            if let Some(r) = p.retirement {
                if r.t_my < p.birth_t_my {
                    return Err(violate(
                        stage,
                        &format!("parcels[{i}].retirement.t_my"),
                        format!("retirement {:?} precedes birth {:?}", r.t_my, p.birth_t_my),
                    ));
                }
            }
        }
        // Lineage is checked in a SECOND pass against the complete ledger:
        // a parent may appear anywhere (chunk order is not lineage order).
        for (i, p) in self.parcels.iter().enumerate() {
            if let Some(parent) = p.lineage_parent {
                let pp = by_id.get(&parent).ok_or_else(|| {
                    violate(
                        stage,
                        &format!("parcels[{i}].lineage_parent"),
                        format!(
                            "parent {:?} not present in ledger (lineage completeness)",
                            parent
                        ),
                    )
                })?;
                if pp.birth_t_my > p.birth_t_my {
                    return Err(violate(
                        stage,
                        &format!("parcels[{i}].lineage_parent"),
                        format!(
                            "parent {:?} born {:?} after child birth {:?} (age monotonicity)",
                            parent, pp.birth_t_my, p.birth_t_my
                        ),
                    ));
                }
            }
        }

        // Gate 1: identity, exact.
        let initial = self
            .parcels
            .iter()
            .filter(|p| {
                matches!(
                    p.birth,
                    ParcelBirthV0::InitialOceanic | ParcelBirthV0::InitialContinental
                )
            })
            .count() as u64;
        let births = self.parcels.len() as u64 - initial;
        let retirements = self
            .parcels
            .iter()
            .filter(|p| p.retirement.is_some())
            .count() as u64;
        let survivors = self.parcels.len() as u64 - retirements;
        let idc = &self.identity;
        if idc.initial_parcels != initial
            || idc.births != births
            || idc.retirements != retirements
            || idc.survivors != survivors
        {
            return Err(violate(
                stage,
                "identity",
                format!(
                    "declared ({}, {}, {}, {}) disagrees with recomputed ({initial}, {births}, {retirements}, {survivors})",
                    idc.initial_parcels, idc.births, idc.retirements, idc.survivors
                ),
            ));
        }
        if !idc.identity_holds() {
            return Err(violate(
                stage,
                "identity",
                format!(
                    "identity failed: {} + {} - {} != {}",
                    idc.initial_parcels, idc.births, idc.retirements, idc.survivors
                ),
            ));
        }

        // Gate 2: material, within declared error, canonical summation
        // order (ascending parcel_id) so verify reproduces the writer.
        let mut ordered: Vec<&CrustParcelRecordV0> = self.parcels.iter().collect();
        ordered.sort_by_key(|p| p.parcel_id);
        let (mut initial_m, mut created_m, mut retired_m, mut surviving_m) = (0.0, 0.0, 0.0, 0.0);
        for p in ordered {
            let q = p.material_area_equivalent_m2;
            match p.birth {
                ParcelBirthV0::InitialOceanic | ParcelBirthV0::InitialContinental => initial_m += q,
                _ => created_m += q,
            }
            if p.retirement.is_some() {
                retired_m += q;
            } else {
                surviving_m += q;
            }
        }
        let mat = &self.material;
        check_finite_nonneg(stage, "material.declared_error", mat.declared_error)?;
        let drift = |declared: f64, summed: f64, field: &str| -> Result<(), DomainError> {
            if !declared.is_finite() || (declared - summed).abs() > mat.declared_error {
                return Err(violate(
                    stage,
                    &format!("material.{field}"),
                    format!("declared {declared:?} vs canonical sum {summed:?} beyond declared_error {:?}", mat.declared_error),
                ));
            }
            Ok(())
        };
        drift(mat.initial_material, initial_m, "initial_material")?;
        drift(mat.created_material, created_m, "created_material")?;
        drift(
            mat.retired_exported_material,
            retired_m,
            "retired_exported_material",
        )?;
        drift(mat.surviving_material, surviving_m, "surviving_material")?;
        if !mat.closes() {
            return Err(violate(
                stage,
                "material",
                format!(
                    "material identity failed: {} + {} - {} != {} +- {}",
                    mat.initial_material,
                    mat.created_material,
                    mat.retired_exported_material,
                    mat.surviving_material,
                    mat.declared_error
                ),
            ));
        }
        if mat.quantity.is_empty() {
            return Err(violate(
                stage,
                "material.quantity",
                "units string required (units table)",
            ));
        }
        Ok(())
    }

    /// Cross-check against the plate table: live parcels sit on live
    /// plates (verify-side, with the domain summary loaded).
    pub fn validate_against_plates(&self, plates: &[PlateRecordV0]) -> Result<(), DomainError> {
        let stage = "ParcelLedgerV0.plates";
        for (i, p) in self.parcels.iter().enumerate() {
            if p.retirement.is_some() {
                continue;
            }
            match plates.iter().find(|q| q.plate_id == p.plate) {
                None => {
                    return Err(violate(
                        stage,
                        &format!("parcels[{i}].plate"),
                        format!("live parcel on unknown plate {:?}", p.plate),
                    ))
                }
                Some(q) if q.retirement.is_some() => {
                    return Err(violate(
                        stage,
                        &format!("parcels[{i}].plate"),
                        format!("live parcel on retired plate {:?}", p.plate),
                    ))
                }
                Some(_) => {}
            }
        }
        Ok(())
    }

    /// Is this parcel alive at `t`? (Projection input.)
    fn alive_at(&self, id: ParcelId, t: f64) -> bool {
        self.parcels.iter().any(|p| {
            p.parcel_id == id && p.birth_t_my <= t && p.retirement.is_none_or(|r| t < r.t_my)
        })
    }
}

/// A derived weighted projection in loaded form. B1 nearest-cell emits one
/// contribution per parcel at weight = 1.0; conservative/PIC remappers may
/// split one parcel across cells — per-parcel weights must close to 1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CellParcelContribution {
    pub element: ElementId,
    pub parcel_id: ParcelId,
    /// Fraction of the parcel's represented material sampled into this
    /// cell; finite, > 0; per-parcel sums close to 1.
    pub weight: f64,
}

impl CellParcelContribution {
    /// Validate a loaded projection against the ledger at one instant.
    /// LEGAL: many parcels per cell; one parcel across several cells with
    /// fractional weights. ILLEGAL: dead/unknown parcels, non-closing
    /// weights, duplicate (element, parcel) contributions.
    pub fn validate_projection(
        contributions: &[CellParcelContribution],
        ledger: &ParcelLedgerV0,
        projected_t_my: f64,
    ) -> Result<(), DomainError> {
        let stage = "CellParcelContribution";
        // Duplicate (element, parcel) pairs: sort u64 tuples, scan.
        let mut keys: Vec<(u64, u64)> = contributions
            .iter()
            .map(|c| (c.element.0, c.parcel_id.0))
            .collect();
        keys.sort_unstable();
        for w in keys.windows(2) {
            if w[0] == w[1] {
                return Err(violate(
                    stage,
                    "*",
                    format!("duplicate contribution for element/parcel {:?}", w[0]),
                ));
            }
        }
        let mut sums: BTreeMap<ParcelId, f64> = BTreeMap::new();
        for (i, c) in contributions.iter().enumerate() {
            check_finite_positive(stage, &format!("[{i}].weight"), c.weight)?;
            if c.weight > 1.0 {
                return Err(violate(
                    stage,
                    &format!("[{i}].weight"),
                    format!("weight > 1.0: {:?}", c.weight),
                ));
            }
            if !ledger.alive_at(c.parcel_id, projected_t_my) {
                return Err(violate(
                    stage,
                    &format!("[{i}].parcel_id"),
                    format!(
                        "parcel {:?} not alive at projected_t_my={projected_t_my:?}",
                        c.parcel_id
                    ),
                ));
            }
            *sums.entry(c.parcel_id).or_insert(0.0) += c.weight;
        }
        for (&pid, &sum) in &sums {
            if !closes_to_one(sum, FRACTION_SUM_TOL) {
                return Err(violate(
                    stage,
                    "weights",
                    format!("parcel {pid:?} weights sum to {sum:?}, not 1.0 (split projections must partition the parcel)"),
                ));
            }
        }
        Ok(())
    }
}

fn validate_chunk_list(
    chunks: &[LedgerChunkRefV0],
    stage: &'static str,
    total: u64,
) -> Result<(), DomainError> {
    for (i, c) in chunks.iter().enumerate() {
        c.validate(stage)
            .map_err(|e| violate(stage, &format!("chunks[{i}]"), e.to_string()))?;
    }
    let summed: u64 = chunks.iter().map(|c| c.record_count).sum();
    if summed != total {
        return Err(violate(
            stage,
            "chunks",
            format!("chunk record_counts sum to {summed}, descriptor says {total}"),
        ));
    }
    Ok(())
}

impl GeodynamicsDomainV0 {
    /// Summary-level validation: bounded tables internally consistent;
    /// descriptors close against their chunk indexes; the two gates hold;
    /// live-parcel counts agree everywhere. Record-level revalidation
    /// happens verify-side when chunks are loaded.
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "GeodynamicsDomainV0";
        check_finite_nonneg(stage, "t_my", self.t_my)?;

        let mut plates: HashMap<PlateId, ()> = HashMap::new();
        for (i, pl) in self.plates.iter().enumerate() {
            if plates.insert(pl.plate_id, ()).is_some() {
                return Err(violate(
                    stage,
                    "plates",
                    format!("duplicate plate_id {:?}", pl.plate_id),
                ));
            }
            check_finite_nonneg(stage, &format!("plates[{i}].birth_t_my"), pl.birth_t_my)?;
            if let Some(parent) = pl.parent_plate {
                if parent == pl.plate_id {
                    return Err(violate(
                        stage,
                        &format!("plates[{i}].parent_plate"),
                        "plate is its own parent",
                    ));
                }
                if !self.plates.iter().any(|q| q.plate_id == parent) {
                    return Err(violate(
                        stage,
                        &format!("plates[{i}].parent_plate"),
                        format!("parent {parent:?} not in table"),
                    ));
                }
            }
            if let Some(pole) = pl.euler_pole {
                pole.validate().map_err(|e| {
                    violate(stage, &format!("plates[{i}].euler_pole"), e.to_string())
                })?;
            }
            if let Some(r) = pl.retirement {
                if pl.euler_pole.is_some() {
                    return Err(violate(
                        stage,
                        &format!("plates[{i}].euler_pole"),
                        "retired plate still carries a rotation",
                    ));
                }
                if r.t_my < pl.birth_t_my {
                    return Err(violate(
                        stage,
                        &format!("plates[{i}].retirement.t_my"),
                        "retirement precedes birth",
                    ));
                }
            }
        }
        if self.plates.is_empty() {
            return Err(violate(stage, "plates", "no plates declared"));
        }

        for (i, b) in self.boundaries.iter().enumerate() {
            check_finite_nonneg(stage, &format!("boundaries[{i}].birth_t_my"), b.birth_t_my)?;
            if b.last_active_t_my < b.birth_t_my {
                return Err(violate(
                    stage,
                    &format!("boundaries[{i}].last_active_t_my"),
                    "activity window is negative",
                ));
            }
            for (name, id) in [("plate_a", b.plate_a), ("plate_b", b.plate_b)] {
                if !plates.contains_key(&id) {
                    return Err(violate(
                        stage,
                        &format!("boundaries[{i}].{name}"),
                        format!("unknown plate {id:?}"),
                    ));
                }
            }
            if b.footprint_element_count == 0 {
                return Err(violate(
                    stage,
                    &format!("boundaries[{i}].footprint_element_count"),
                    "boundary with empty footprint",
                ));
            }
        }

        let ev = &self.event_ledger;
        if ev.schema.is_empty() {
            return Err(violate(stage, "event_ledger.schema", "schema tag required"));
        }
        validate_chunk_list(&ev.chunks, stage, ev.total_events)?;
        match (ev.t_start_my, ev.t_end_my) {
            (Some(a), Some(b)) if b < a => {
                return Err(violate(
                    stage,
                    "event_ledger",
                    format!("t_end {b:?} precedes t_start {a:?}"),
                ));
            }
            _ => {}
        }

        let pl = &self.parcel_ledger;
        if pl.schema.is_empty() {
            return Err(violate(
                stage,
                "parcel_ledger.schema",
                "schema tag required",
            ));
        }
        validate_chunk_list(&pl.chunks, stage, pl.total_parcels)?;
        if pl.live_parcels > pl.total_parcels {
            return Err(violate(
                stage,
                "parcel_ledger.live_parcels",
                "live parcels exceed total parcels",
            ));
        }

        let occ = &self.occupancy;
        if occ.projected_t_my != self.t_my {
            return Err(violate(
                stage,
                "occupancy.projected_t_my",
                format!(
                    "projection at {:?} does not match domain t_my {:?} (stale projection)",
                    occ.projected_t_my, self.t_my
                ),
            ));
        }
        validate_chunk_list(&occ.chunks, stage, occ.contribution_count)?;
        if occ.projected_parcel_count > pl.live_parcels {
            return Err(violate(
                stage,
                "occupancy.projected_parcel_count",
                format!(
                    "{} projected parcels exceed {} live parcels",
                    occ.projected_parcel_count, pl.live_parcels
                ),
            ));
        }

        // The gates, and cross-agreement of survivor counts.
        if !self.identity_conservation.identity_holds() {
            return Err(violate(
                stage,
                "identity_conservation",
                format!(
                    "identity failed: {} + {} - {} != {}",
                    self.identity_conservation.initial_parcels,
                    self.identity_conservation.births,
                    self.identity_conservation.retirements,
                    self.identity_conservation.survivors
                ),
            ));
        }
        if self.identity_conservation.survivors != pl.live_parcels {
            return Err(violate(
                stage,
                "identity_conservation.survivors",
                format!(
                    "{} survivors but parcel ledger reports {} live",
                    self.identity_conservation.survivors, pl.live_parcels
                ),
            ));
        }
        if !self.material_conservation.closes() {
            return Err(violate(
                stage,
                "material_conservation",
                format!(
                    "material failed: {} + {} - {} != {} +- {}",
                    self.material_conservation.initial_material,
                    self.material_conservation.created_material,
                    self.material_conservation.retired_exported_material,
                    self.material_conservation.surviving_material,
                    self.material_conservation.declared_error
                ),
            ));
        }
        if self.material_conservation.quantity.is_empty() {
            return Err(violate(
                stage,
                "material_conservation.quantity",
                "units string required",
            ));
        }
        Ok(())
    }

    /// Bind every chunk ref to the artifact's manifest sidecars: each
    /// `LedgerChunkRefV0::sidecar` must resolve to exactly one
    /// `SidecarDescriptor` logical name, and when the descriptor carries a
    /// record count it must AGREE with the chunk ref's count. Generic
    /// verify already guarantees the sidecar BYTES; this closes the
    /// domain accounting against the same one truth.
    pub fn validate_sidecar_binding(
        &self,
        sidecars: &[SidecarDescriptor],
    ) -> Result<(), DomainError> {
        let stage = "GeodynamicsDomainV0.sidecars";
        let families: [(&str, &str, &[LedgerChunkRefV0], u64); 3] = [
            (
                "event_ledger",
                "events",
                &self.event_ledger.chunks,
                self.event_ledger.total_events,
            ),
            (
                "parcel_ledger",
                "parcels",
                &self.parcel_ledger.chunks,
                self.parcel_ledger.total_parcels,
            ),
            (
                "occupancy",
                "contributions",
                &self.occupancy.chunks,
                self.occupancy.contribution_count,
            ),
        ];
        for (family, what, chunks, total) in families {
            let mut summed: u64 = 0;
            for c in chunks {
                c.validate(stage)
                    .map_err(|e| violate(stage, &format!("{family}.chunks"), e.to_string()))?;
                let desc = sidecars
                    .iter()
                    .find(|s| s.logical_name == c.sidecar)
                    .ok_or_else(|| {
                        violate(
                            stage,
                            &format!("{family}.chunks"),
                            format!(
                                "chunk sidecar {:?} not declared in the semantic manifest",
                                c.sidecar
                            ),
                        )
                    })?;
                if let Some(n) = desc.record_count {
                    if n != c.record_count {
                        return Err(violate(
                            stage,
                            &format!("{family}.chunks"),
                            format!(
                                "sidecar {:?} declares {n} {what} but the chunk ref says {}",
                                c.sidecar, c.record_count
                            ),
                        ));
                    }
                }
                summed += c.record_count;
            }
            if summed != total {
                return Err(violate(
                    stage,
                    family,
                    format!("chunks account for {summed} {what}, descriptor says {total}"),
                ));
            }
        }
        Ok(())
    }
}

impl HydrologyDomainV0 {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "HydrologyDomainV0";
        let mut body_ids = std::collections::BTreeSet::new();
        for b in &self.water_bodies {
            if b.body_id == u32::MAX
                || !body_ids.insert(b.body_id)
                || !b.water_surface_m.is_finite()
                || !b.area_m2.is_finite()
                || b.area_m2 <= 0.0
            {
                return Err(violate(
                    stage,
                    "water_bodies",
                    "invalid/duplicate water body",
                ));
            }
        }
        check_finite_nonneg(stage, "t_my", self.t_my)?;

        let mut basins: HashMap<BasinId, &BasinRecordV0> = HashMap::new();
        for (i, b) in self.basins.iter().enumerate() {
            if basins.insert(b.basin_id, b).is_some() {
                return Err(violate(
                    stage,
                    "basins",
                    format!("duplicate basin_id {:?}", b.basin_id),
                ));
            }
            if let Some(spill) = b.spill_level_m {
                if !spill.is_finite() {
                    return Err(violate(
                        stage,
                        &format!("basins[{i}].spill_level_m"),
                        "nonfinite",
                    ));
                }
            }
        }

        let mut reaches: HashMap<ReachId, &RiverReachV0> = HashMap::new();
        for (i, r) in self.rivers.iter().enumerate() {
            if reaches.insert(r.reach_id, r).is_some() {
                return Err(violate(
                    stage,
                    "rivers",
                    format!("duplicate reach_id {:?}", r.reach_id),
                ));
            }
            if r.strahler_order < 1 {
                return Err(violate(
                    stage,
                    &format!("rivers[{i}].strahler_order"),
                    "Strahler order starts at 1",
                ));
            }
            check_finite_positive(
                stage,
                &format!("rivers[{i}].mean_discharge_proxy_m3_s"),
                r.mean_discharge_proxy_m3_s,
            )?;
            match (r.downstream, r.outlet_basin) {
                (Some(_), Some(_)) => {
                    return Err(violate(
                        stage,
                        &format!("rivers[{i}]"),
                        "reach cannot flow downstream AND terminate in a basin",
                    ));
                }
                (None, None) => {
                    return Err(violate(
                        stage,
                        &format!("rivers[{i}]"),
                        "reach must flow somewhere: downstream or outlet_basin",
                    ));
                }
                _ => {}
            }
        }
        for (i, r) in self.rivers.iter().enumerate() {
            if let Some(d) = r.downstream {
                if d == r.reach_id {
                    return Err(violate(
                        stage,
                        &format!("rivers[{i}].downstream"),
                        "reach is its own downstream (cycle)",
                    ));
                }
                if !reaches.contains_key(&d) {
                    return Err(violate(
                        stage,
                        &format!("rivers[{i}].downstream"),
                        format!("unknown reach {d:?}"),
                    ));
                }
            }
            if let Some(b) = r.outlet_basin {
                if !basins.contains_key(&b) {
                    return Err(violate(
                        stage,
                        &format!("rivers[{i}].outlet_basin"),
                        format!("unknown basin {b:?}"),
                    ));
                }
            }
        }
        for (i, b) in self.basins.iter().enumerate() {
            for &rin in &b.inlet_reaches {
                if !reaches.contains_key(&rin) {
                    return Err(violate(
                        stage,
                        &format!("basins[{i}].inlet_reaches"),
                        format!("unknown reach {rin:?}"),
                    ));
                }
            }
        }

        for (i, w) in self.waterfalls.iter().enumerate() {
            check_finite_positive(stage, &format!("waterfalls[{i}].drop_m"), w.drop_m)?;
            if let Some(r) = w.reach {
                if !reaches.contains_key(&r) {
                    return Err(violate(
                        stage,
                        &format!("waterfalls[{i}].reach"),
                        format!("unknown reach {r:?}"),
                    ));
                }
            }
        }
        if let Some(cf) = &self.cell_fields {
            for (name, v) in [
                ("flow_to_field", &cf.flow_to_field),
                ("basin_field", &cf.basin_field),
                ("reach_field", &cf.reach_field),
                ("fill_level_field", &cf.fill_level_field),
                ("contributing_area_field", &cf.contributing_area_field),
                ("discharge_field", &cf.discharge_field),
                ("lake_level_convention", &cf.lake_level_convention),
            ] {
                if v.is_empty() {
                    return Err(violate(stage, &format!("cell_fields.{name}"), "empty"));
                }
            }
        }
        Ok(())
    }
}

impl GrainClassSplit {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "GrainClassSplit";
        let parts = [
            ("coarse_fraction", self.coarse_fraction),
            ("sand_fraction", self.sand_fraction),
            ("silt_fraction", self.silt_fraction),
            ("clay_fraction", self.clay_fraction),
        ];
        for (name, v) in parts {
            if !v.is_finite() || !(0.0..=1.0).contains(&v) {
                return Err(violate(
                    stage,
                    name,
                    format!("fraction outside [0, 1]: {v:?}"),
                ));
            }
        }
        let sum =
            self.coarse_fraction + self.sand_fraction + self.silt_fraction + self.clay_fraction;
        if !closes_to_one(sum, FRACTION_SUM_TOL) {
            return Err(violate(
                stage,
                "*",
                format!("fractions sum to {sum:?}, not 1.0"),
            ));
        }
        Ok(())
    }
}

impl PhaseSplit {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "PhaseSplit";
        for (name, v) in [
            ("particulate_fraction", self.particulate_fraction),
            ("dissolved_fraction", self.dissolved_fraction),
        ] {
            if !v.is_finite() || !(0.0..=1.0).contains(&v) {
                return Err(violate(
                    stage,
                    name,
                    format!("fraction outside [0, 1]: {v:?}"),
                ));
            }
        }
        let sum = self.particulate_fraction + self.dissolved_fraction;
        if !closes_to_one(sum, FRACTION_SUM_TOL) {
            return Err(violate(
                stage,
                "*",
                format!(
                    "phases sum to {sum:?}, not 1.0 (§54.8 split must be declared, not defaulted)"
                ),
            ));
        }
        Ok(())
    }
}

impl SedimentTransferRecord {
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "SedimentTransferRecord";
        check_finite_nonneg(stage, "epoch_t_my", self.epoch_t_my)?;
        check_finite_positive(stage, "mass_kg", self.mass_kg)?;
        if self.source_region.is_empty() {
            return Err(violate(stage, "source_region", "empty source region"));
        }
        if self.source_material_lineages.is_empty() {
            return Err(violate(
                stage,
                "source_material_lineages",
                "no lineage provenance (§54.2)",
            ));
        }
        if self.source_mixture.is_empty() {
            return Err(violate(
                stage,
                "source_mixture",
                "no provenance mixture declared",
            ));
        }
        let mut sum = 0.0;
        for (i, f) in self.source_mixture.iter().enumerate() {
            if f.lineage.is_empty() {
                return Err(violate(
                    stage,
                    &format!("source_mixture[{i}].lineage"),
                    "empty lineage label",
                ));
            }
            if !f.fraction.is_finite() || f.fraction <= 0.0 {
                return Err(violate(
                    stage,
                    &format!("source_mixture[{i}].fraction"),
                    format!("expected > 0, found {:?}", f.fraction),
                ));
            }
            sum += f.fraction;
        }
        if !closes_to_one(sum, FRACTION_SUM_TOL) {
            return Err(violate(
                stage,
                "source_mixture",
                format!("mixture sums to {sum:?}, not 1.0"),
            ));
        }
        self.grain_classes.validate()?;
        self.phase_split.validate()?;
        if let Some(d) = &self.deposition {
            check_finite_nonneg(stage, "deposition.deposited_t_my", d.deposited_t_my)?;
            if d.deposited_t_my < self.epoch_t_my {
                return Err(violate(
                    stage,
                    "deposition.deposited_t_my",
                    "deposition precedes the transfer it closes",
                ));
            }
            match &d.sink {
                SedimentSinkV0::Basin { .. }
                | SedimentSinkV0::Floodplain
                | SedimentSinkV0::Delta => {
                    if d.destination_elements.is_empty() {
                        return Err(violate(
                            stage,
                            "deposition.destination_elements",
                            "surface sinks must name destination cells",
                        ));
                    }
                }
                SedimentSinkV0::OceanExport | SedimentSinkV0::DeepSink => {}
            }
        }
        Ok(())
    }
}

impl SedimentMassBalanceV0 {
    /// §54.4 gate: per-epoch and cumulative closure within declared error.
    pub fn validate(&self) -> Result<(), DomainError> {
        let stage = "SedimentMassBalanceV0";
        let (mut er, mut dep, mut exp) = (0.0, 0.0, 0.0);
        for (i, e) in self.epochs.iter().enumerate() {
            check_finite_nonneg(stage, &format!("epochs[{i}].eroded_kg"), e.eroded_kg)?;
            check_finite_nonneg(stage, &format!("epochs[{i}].deposited_kg"), e.deposited_kg)?;
            check_finite_nonneg(stage, &format!("epochs[{i}].exported_kg"), e.exported_kg)?;
            check_finite_nonneg(
                stage,
                &format!("epochs[{i}].declared_error_kg"),
                e.declared_error_kg,
            )?;
            let residual = e.eroded_kg - e.deposited_kg - e.exported_kg;
            if residual.abs() > e.declared_error_kg {
                return Err(violate(
                    stage,
                    &format!("epochs[{i}]"),
                    format!(
                        "eroded {} != deposited {} + exported {} +- {} (residual {residual:?})",
                        e.eroded_kg, e.deposited_kg, e.exported_kg, e.declared_error_kg
                    ),
                ));
            }
            er += e.eroded_kg;
            dep += e.deposited_kg;
            exp += e.exported_kg;
        }
        let close = |declared: f64, summed: f64, field: &str| -> Result<(), DomainError> {
            if (declared - summed).abs() > self.cumulative_declared_error_kg {
                return Err(violate(
                    stage,
                    field,
                    format!("declared cumulative {declared:?} disagrees with epoch sum {summed:?} beyond declared error {:?}", self.cumulative_declared_error_kg),
                ));
            }
            Ok(())
        };
        close(self.cumulative_eroded_kg, er, "cumulative_eroded_kg")?;
        close(self.cumulative_deposited_kg, dep, "cumulative_deposited_kg")?;
        close(self.cumulative_exported_kg, exp, "cumulative_exported_kg")?;
        let residual =
            self.cumulative_eroded_kg - self.cumulative_deposited_kg - self.cumulative_exported_kg;
        if residual.abs() > self.cumulative_declared_error_kg {
            return Err(violate(
                stage,
                "cumulative",
                format!(
                    "cumulative residual {residual:?} exceeds declared error {:?}",
                    self.cumulative_declared_error_kg
                ),
            ));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests

#[cfg(test)]
mod tests {
    use super::*;

    fn pole() -> EulerPoleV0 {
        EulerPoleV0 {
            lat_deg: 30.0,
            lon_deg: -45.0,
            rate_deg_per_my: 0.8,
        }
    }

    fn parcel(id: u64, birth_t: f64, birth: ParcelBirthV0, area: f64) -> CrustParcelRecordV0 {
        CrustParcelRecordV0 {
            parcel_id: ParcelId(id),
            birth_t_my: birth_t,
            birth,
            lineage_parent: None,
            plate: PlateId(1),
            material_area_equivalent_m2: area,
            retirement: None,
        }
    }

    fn ledger_of(parcels: &[CrustParcelRecordV0]) -> ParcelLedgerV0 {
        let initial = parcels
            .iter()
            .filter(|p| {
                matches!(
                    p.birth,
                    ParcelBirthV0::InitialOceanic | ParcelBirthV0::InitialContinental
                )
            })
            .count() as u64;
        let retirements = parcels.iter().filter(|p| p.retirement.is_some()).count() as u64;
        let mut ordered: Vec<&CrustParcelRecordV0> = parcels.iter().collect();
        ordered.sort_by_key(|p| p.parcel_id);
        let (mut i_m, mut c_m, mut r_m, mut s_m) = (0.0, 0.0, 0.0, 0.0);
        for p in &ordered {
            match p.birth {
                ParcelBirthV0::InitialOceanic | ParcelBirthV0::InitialContinental => {
                    i_m += p.material_area_equivalent_m2
                }
                _ => c_m += p.material_area_equivalent_m2,
            }
            if p.retirement.is_some() {
                r_m += p.material_area_equivalent_m2
            } else {
                s_m += p.material_area_equivalent_m2
            }
        }
        ParcelLedgerV0 {
            parcels: parcels.to_vec(),
            identity: ParcelIdentityConservationV0 {
                t_my: 10.0,
                initial_parcels: initial,
                births: parcels.len() as u64 - initial,
                retirements,
                survivors: parcels.len() as u64 - retirements,
            },
            material: MaterialConservationV0 {
                t_my: 10.0,
                quantity: "m2".into(),
                initial_material: i_m,
                created_material: c_m,
                retired_exported_material: r_m,
                surviving_material: s_m,
                declared_error: 1.0e-6,
            },
        }
    }

    fn chunk(logical: &str, count: u64) -> LedgerChunkRefV0 {
        LedgerChunkRefV0 {
            sidecar: logical.into(),
            record_count: count,
        }
    }

    fn geodynamics(
        t_my: f64,
        ledger: &ParcelLedgerV0,
        contributions: usize,
    ) -> GeodynamicsDomainV0 {
        let live = ledger
            .parcels
            .iter()
            .filter(|p| p.retirement.is_none())
            .count() as u64;
        GeodynamicsDomainV0 {
            t_my,
            plates: vec![PlateRecordV0 {
                plate_id: PlateId(1),
                parent_plate: None,
                birth_t_my: 0.0,
                euler_pole: Some(pole()),
                retirement: None,
            }],
            boundaries: vec![BoundaryRecordV0 {
                boundary_id: 7,
                kind: BoundaryKindV0::Ridge,
                plate_a: PlateId(1),
                plate_b: PlateId(1),
                birth_t_my: 0.0,
                last_active_t_my: t_my,
                footprint_element_count: 12,
            }],
            event_ledger: EventLedgerDescriptorV0 {
                schema: "pf.events/v0".into(),
                total_events: 1,
                t_start_my: Some(2.0),
                t_end_my: Some(2.0),
                chunks: vec![chunk("events.chunk-000", 1)],
            },
            parcel_ledger: ParcelLedgerDescriptorV0 {
                schema: "pf.parcels/v0".into(),
                total_parcels: ledger.parcels.len() as u64,
                live_parcels: live,
                chunks: vec![chunk("parcels.chunk-000", ledger.parcels.len() as u64)],
            },
            occupancy: ProjectionDescriptorV0 {
                projected_t_my: t_my,
                contribution_count: contributions as u64,
                projected_parcel_count: live.min(contributions as u64),
                chunks: vec![chunk("occupancy.chunk-000", contributions as u64)],
            },
            identity_conservation: ledger.identity,
            material_conservation: ledger.material.clone(),
        }
    }

    #[test]
    fn identity_and_material_are_separate_gates() {
        let mut retired = parcel(2, 4.0, ParcelBirthV0::RidgeCreated, 5.0e9);
        retired.retirement = Some(ParcelRetirementV0Detail {
            t_my: 8.0,
            cause: ParcelRetirementCauseV0::Subducted,
            event: EventId(1),
        });
        let parcels = vec![
            parcel(0, 0.0, ParcelBirthV0::InitialOceanic, 3.0e9),
            parcel(1, 0.0, ParcelBirthV0::InitialContinental, 4.0e9),
            retired,
        ];
        let good = ledger_of(&parcels);
        good.validate().unwrap();

        // Identity lie: inflate survivors — exact gate must trip.
        let mut lied_identity = good.clone();
        lied_identity.identity.survivors += 1;
        let err = lied_identity.validate().unwrap_err();
        assert!(err.to_string().contains("identity"), "{err}");

        // Material lie: survivors keep material that was retired — the
        // quantity gate trips even though counts are honest.
        let mut lied_material = good.clone();
        lied_material.material.surviving_material += 5.0e9;
        lied_material.material.declared_error = 1.0e-6;
        let err = lied_material.validate().unwrap_err();
        assert!(err.to_string().contains("surviving_material"), "{err}");

        // Small numerical drift within declared error is legitimate.
        let mut drifted = good.clone();
        drifted.material.declared_error = 1.0e3;
        drifted.material.surviving_material += 500.0;
        drifted.validate().unwrap();
    }

    #[test]
    fn projection_supports_nearest_cell_and_conservative_split() {
        let parcels = vec![
            parcel(0, 0.0, ParcelBirthV0::InitialOceanic, 3.0e9),
            parcel(1, 0.0, ParcelBirthV0::InitialContinental, 4.0e9),
        ];
        let ledger = ledger_of(&parcels);

        // B1 nearest-cell: one contribution per parcel, weight 1. Two
        // parcels in ONE cell: both survive (PARCEL != CELL).
        let nearest = vec![
            CellParcelContribution {
                element: ElementId(42),
                parcel_id: ParcelId(0),
                weight: 1.0,
            },
            CellParcelContribution {
                element: ElementId(42),
                parcel_id: ParcelId(1),
                weight: 1.0,
            },
        ];
        CellParcelContribution::validate_projection(&nearest, &ledger, 10.0).unwrap();
        assert_eq!(ledger.identity.survivors, 2);

        // Later conservative/PIC remapper: parcel 0 split 0.63/0.37
        // across two cells — same artifact language, still valid.
        let split = vec![
            CellParcelContribution {
                element: ElementId(41),
                parcel_id: ParcelId(0),
                weight: 0.63,
            },
            CellParcelContribution {
                element: ElementId(42),
                parcel_id: ParcelId(0),
                weight: 0.37,
            },
            CellParcelContribution {
                element: ElementId(42),
                parcel_id: ParcelId(1),
                weight: 1.0,
            },
        ];
        CellParcelContribution::validate_projection(&split, &ledger, 10.0).unwrap();

        // Weights that do not partition the parcel are refused.
        let leaky = vec![
            CellParcelContribution {
                element: ElementId(41),
                parcel_id: ParcelId(0),
                weight: 0.63,
            },
            CellParcelContribution {
                element: ElementId(42),
                parcel_id: ParcelId(0),
                weight: 0.30,
            },
        ];
        let err = CellParcelContribution::validate_projection(&leaky, &ledger, 10.0).unwrap_err();
        assert!(err.to_string().contains("not 1.0"), "{err}");

        // Dead parcel in the projection is refused.
        let mut dead = parcel(5, 0.0, ParcelBirthV0::RidgeCreated, 1.0e9);
        dead.retirement = Some(ParcelRetirementV0Detail {
            t_my: 3.0,
            cause: ParcelRetirementCauseV0::Subducted,
            event: EventId(1),
        });
        let mut with_dead = parcels.clone();
        with_dead.push(dead);
        let ledger2 = ledger_of(&with_dead);
        let ghost = vec![CellParcelContribution {
            element: ElementId(1),
            parcel_id: ParcelId(5),
            weight: 1.0,
        }];
        let err = CellParcelContribution::validate_projection(&ghost, &ledger2, 10.0).unwrap_err();
        assert!(err.to_string().contains("not alive"), "{err}");
    }

    #[test]
    fn descriptors_close_against_chunks_and_gates() {
        let parcels = vec![parcel(0, 0.0, ParcelBirthV0::InitialOceanic, 3.0e9)];
        let ledger = ledger_of(&parcels);
        let g = geodynamics(10.0, &ledger, 1);
        g.validate().unwrap();

        // Chunk counts must sum to the descriptor total.
        let mut bad_chunks = g.clone();
        bad_chunks.parcel_ledger.chunks[0].record_count = 2;
        let err = bad_chunks.validate().unwrap_err();
        assert!(err.to_string().contains("sum to"), "{err}");

        // Survivors must agree with the live-parcel descriptor (identity
        // itself still holds — this isolates the cross-check).
        let mut disagree = g.clone();
        disagree.parcel_ledger.live_parcels = 2;
        let err = disagree.validate().unwrap_err();
        assert!(err.to_string().contains("live"), "{err}");

        // Stale projection instant refused.
        let mut stale = g.clone();
        stale.occupancy.projected_t_my = 9.5;
        let err = stale.validate().unwrap_err();
        assert!(err.to_string().contains("stale"), "{err}");
    }

    #[test]
    fn chunk_refs_bind_to_manifest_sidecars() {
        use crate::artifact::SidecarDescriptor;
        let parcels = vec![parcel(0, 0.0, ParcelBirthV0::InitialOceanic, 3.0e9)];
        let ledger = ledger_of(&parcels);
        let g = geodynamics(10.0, &ledger, 1);

        let manifest_sidecars = |events: Option<u64>| {
            vec![SidecarDescriptor {
                logical_name: "events.chunk-000".into(),
                path: "ledger/events-000.json".into(),
                schema: "pf.events/v0".into(),
                byte_length: 128,
                digest: "blake3:sealed".into(),
                record_count: events,
            }]
        };
        // Manifest knows only SOME of the referenced sidecars: the parcel
        // and occupancy chunk refs dangle -> loud refusal.
        let err = g
            .validate_sidecar_binding(&manifest_sidecars(Some(1)))
            .unwrap_err();
        assert!(
            err.to_string()
                .contains("not declared in the semantic manifest"),
            "{err}"
        );

        // Fully declared manifest binds cleanly...
        let full = vec![
            SidecarDescriptor {
                logical_name: "events.chunk-000".into(),
                path: "ledger/events-000.json".into(),
                schema: "pf.events/v0".into(),
                byte_length: 128,
                digest: "blake3:a".into(),
                record_count: Some(1),
            },
            SidecarDescriptor {
                logical_name: "parcels.chunk-000".into(),
                path: "ledger/parcels-000.json".into(),
                schema: "pf.parcels/v0".into(),
                byte_length: 256,
                digest: "blake3:b".into(),
                record_count: Some(1),
            },
            SidecarDescriptor {
                logical_name: "occupancy.chunk-000".into(),
                path: "ledger/occupancy-000.json".into(),
                schema: "pf.occupancy/v0".into(),
                byte_length: 64,
                digest: "blake3:c".into(),
                record_count: Some(1),
            },
        ];
        g.validate_sidecar_binding(&full).unwrap();

        // ...but a record-count disagreement between chunk ref and the
        // manifest descriptor is a two-truths violation -> loud refusal.
        let mut conflicting = full.clone();
        conflicting[0].record_count = Some(7);
        let err = g.validate_sidecar_binding(&conflicting).unwrap_err();
        assert!(err.to_string().contains("but the chunk ref says"), "{err}");
    }

    #[test]
    fn lineage_completeness_and_monotonicity() {
        let mut orphan = parcel(5, 6.0, ParcelBirthV0::Accreted, 1.0e9);
        orphan.lineage_parent = Some(ParcelId(77));
        let err = ledger_of(&[orphan]).validate().unwrap_err();
        assert!(err.to_string().contains("not present in ledger"), "{err}");

        let parent = parcel(1, 0.0, ParcelBirthV0::InitialOceanic, 3.0e9);
        let mut child = parcel(2, 1.0, ParcelBirthV0::Accreted, 2.0e9);
        child.lineage_parent = Some(ParcelId(1));
        ledger_of(&[child.clone(), parent]).validate().unwrap();

        let parent_late = parcel(9, 5.0, ParcelBirthV0::RidgeCreated, 1.0e9);
        let mut child_early = parcel(2, 1.0, ParcelBirthV0::Accreted, 2.0e9);
        child_early.lineage_parent = Some(ParcelId(9));
        let err = ledger_of(&[child_early, parent_late])
            .validate()
            .unwrap_err();
        assert!(err.to_string().contains("age monotonicity"), "{err}");
    }

    #[test]
    fn live_parcel_on_retired_plate_refused() {
        let parcels = vec![parcel(0, 0.0, ParcelBirthV0::InitialOceanic, 3.0e9)];
        let ledger = ledger_of(&parcels);
        let mut g = geodynamics(10.0, &ledger, 1);
        g.plates[0].retirement = Some(PlateRetirementDetailV0 {
            t_my: 5.0,
            cause: PlateRetirementV0::Subducted,
        });
        g.plates[0].euler_pole = None;
        g.validate().unwrap(); // summary alone stays consistent...
        let ledger = ledger_of(&parcels);
        let err = ledger.validate_against_plates(&g.plates).unwrap_err();
        assert!(err.to_string().contains("retired plate"), "{err}");
    }

    fn hydrology() -> HydrologyDomainV0 {
        HydrologyDomainV0 {
            water_bodies: Vec::new(),
            exchange_history: None,
            t_my: 10.0,
            basins: vec![BasinRecordV0 {
                basin_id: BasinId(1),
                classification: DrainageClassV0::Endorheic,
                spill_level_m: Some(120.5),
                outlet_element: Some(ElementId(9)),
                inlet_reaches: vec![ReachId(2)],
            }],
            rivers: vec![
                RiverReachV0 {
                    reach_id: ReachId(1),
                    downstream: Some(ReachId(2)),
                    outlet_basin: None,
                    strahler_order: 3,
                    mean_discharge_proxy_m3_s: 42.0,
                },
                RiverReachV0 {
                    reach_id: ReachId(2),
                    downstream: None,
                    outlet_basin: Some(BasinId(1)),
                    strahler_order: 4,
                    mean_discharge_proxy_m3_s: 90.0,
                },
            ],
            waterfalls: vec![WaterfallObservation {
                element: ElementId(11),
                drop_m: 28.0,
                cause: WaterfallCause::UpliftStep,
                reach: Some(ReachId(1)),
            }],
            cell_fields: None,
        }
    }

    #[test]
    fn hydrology_validates_and_causes_are_snake_case() {
        hydrology().validate().unwrap();
        // Cause labels serialize as the §57.5 wire names.
        assert_eq!(
            serde_json::to_string(&WaterfallCause::DeformationBand).unwrap(),
            "\"deformation_band\""
        );
        assert_eq!(
            serde_json::to_string(&WaterfallCause::HangingValley).unwrap(),
            "\"hanging_valley\""
        );

        // Zero/negative drop is not a waterfall.
        let mut h = hydrology();
        h.waterfalls[0].drop_m = 0.0;
        let err = h.validate().unwrap_err();
        assert!(err.to_string().contains("drop_m"), "{err}");

        // A reach going nowhere refused.
        let mut h = hydrology();
        h.rivers[0].downstream = None;
        h.rivers[0].outlet_basin = None;
        let err = h.validate().unwrap_err();
        assert!(err.to_string().contains("flow somewhere"), "{err}");
    }

    fn sediment() -> SedimentTransferRecord {
        SedimentTransferRecord {
            record_id: 1,
            epoch_t_my: 12.0,
            source_region: vec![ElementId(5), ElementId(6)],
            source_material_lineages: vec![ParcelId(0), ParcelId(3)],
            source_mixture: vec![
                ProvenanceFraction {
                    lineage: "granite-arc".into(),
                    fraction: 0.72,
                },
                ProvenanceFraction {
                    lineage: "metamorphic-highlands".into(),
                    fraction: 0.28,
                },
            ],
            mass_kg: 1.0e12,
            grain_classes: GrainClassSplit {
                coarse_fraction: 0.1,
                sand_fraction: 0.5,
                silt_fraction: 0.3,
                clay_fraction: 0.1,
            },
            phase_split: PhaseSplit {
                particulate_fraction: 0.97,
                dissolved_fraction: 0.03,
            },
            transport_pathway: TransportPathway::Fluvial,
            deposition: Some(SedimentDepositionV0 {
                deposited_t_my: 13.0,
                sink: SedimentSinkV0::Basin { basin: BasinId(1) },
                destination_elements: vec![ElementId(7)],
            }),
            exported_deep_sink: false,
        }
    }

    #[test]
    fn sediment_record_validates_and_splits_must_close() {
        sediment().validate().unwrap();

        let mut s = sediment();
        s.phase_split = PhaseSplit {
            particulate_fraction: 0.6,
            dissolved_fraction: 0.5,
        };
        let err = s.validate().unwrap_err();
        assert!(err.to_string().contains("phases sum"), "{err}");

        let mut s = sediment();
        s.grain_classes.clay_fraction = 0.2;
        let err = s.validate().unwrap_err();
        assert!(err.to_string().contains("sum to"), "{err}");

        let mut s = sediment();
        s.source_mixture[1].fraction = 0.1;
        let err = s.validate().unwrap_err();
        assert!(err.to_string().contains("mixture sums"), "{err}");
        let mut s = sediment();
        s.source_material_lineages = vec![];
        let err = s.validate().unwrap_err();
        assert!(err.to_string().contains("lineage provenance"), "{err}");

        let mut s = sediment();
        s.deposition = Some(SedimentDepositionV0 {
            deposited_t_my: 13.0,
            sink: SedimentSinkV0::Delta,
            destination_elements: vec![],
        });
        let err = s.validate().unwrap_err();
        assert!(err.to_string().contains("destination cells"), "{err}");
    }

    #[test]
    fn sediment_mass_balance_gate() {
        let ok = SedimentMassBalanceV0 {
            epochs: vec![SedimentEpochBalanceV0 {
                t_my: 12.0,
                eroded_kg: 100.0,
                deposited_kg: 80.0,
                exported_kg: 20.0,
                declared_error_kg: 0.0,
            }],
            cumulative_eroded_kg: 100.0,
            cumulative_deposited_kg: 80.0,
            cumulative_exported_kg: 20.0,
            cumulative_declared_error_kg: 0.0,
        };
        ok.validate().unwrap();

        let mut bad = ok.clone();
        bad.epochs[0].deposited_kg = 70.0;
        let err = bad.validate().unwrap_err();
        assert!(err.to_string().contains("eroded"), "{err}");

        let mut bad2 = ok;
        bad2.cumulative_exported_kg = 25.0;
        let err = bad2.validate().unwrap_err();
        assert!(err.to_string().contains("cumulative_exported_kg"), "{err}");
    }

    #[test]
    fn schemas_roundtrip_and_chunk_digests_are_stable() {
        let parcels = vec![parcel(0, 0.0, ParcelBirthV0::InitialOceanic, 3.0e9)];
        let ledger = ledger_of(&parcels);
        let g = geodynamics(10.0, &ledger, 1);
        let a = serde_json::to_string(&g).unwrap();
        let b = serde_json::to_string(&g).unwrap();
        assert_eq!(a, b);
        let back: GeodynamicsDomainV0 = serde_json::from_str(&a).unwrap();
        assert_eq!(back, g);

        let h = hydrology();
        let ha = serde_json::to_string(&h).unwrap();
        assert_eq!(serde_json::from_str::<HydrologyDomainV0>(&ha).unwrap(), h);

        let s = sediment();
        let sa = serde_json::to_string(&s).unwrap();
        assert_eq!(
            serde_json::from_str::<SedimentTransferRecord>(&sa).unwrap(),
            s
        );
    }
}
