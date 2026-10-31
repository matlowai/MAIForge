//! pf-contracts: typed artifact contracts, parameter system, digests, receipts.
//!
//! B1 rev 2 (see docs/plans/PLANET_FOUNDRY_B1_BUILD_BRIEF.md). Fail-loud
//! everywhere: no silent clamping, default-value substitution, NaN scrubbing
//! or fallback profiles. Every error carries stage/field/index context.

pub mod planet;
pub mod world_admission;
pub mod column_envelope;
pub mod artifact;
pub mod emplacement;
pub mod volcanic_supply;
pub mod source_chemistry;
pub mod layer_provenance;
pub mod domain;
pub mod laws;
pub mod params;
/// WS-MODULE: the SimModule contract (ADR 039 §4, amendments C2-C4).
pub mod sim_module;
pub mod spatial;
pub mod terrain;
pub mod region;

pub use artifact::{
    verify_artifact, ArtifactHeader, ArtifactKind, ArtifactOnDisk, BackendId, Dtype,
    ExecutionReceipt, FieldDescriptor, SemanticManifest, SidecarDescriptor, Tier,
};
pub use domain::HydrologyCellFieldsV0;
pub use terrain::{
    AssumptionSourceV0, AssumptionV0, AvailabilityV0, BoundaryClassCodeV0, CrustKindCodeV0,
    CrustVolumeLedgerV0, ErosionRecordV0, RegionDomainV1, RegionFrameV1,
    RegionReachV1, RegionSourceV0, ResamplingV0, ResourceFamilyKindV0, ResourceFamilyV0,
    ResourcePotentialDomainV0, RockFamilyCodeV0, TerrainCodesV0, TerrainDomainV0,
    TerrainSummaryV0, UpliftPriorV0, RESOURCE_PREVIEW_USAGE_LIMITS,
};
pub use domain::{
    BasinId, BasinRecordV0, BoundaryKindV0, BoundaryRecordV0, CellParcelContribution,
    CrustParcelRecordV0, DrainageClassV0, EulerPoleV0, EventId, EventKindV0,
    EventLedgerDescriptorV0, GeodynamicsDomainV0, GrainClassSplit, HydrologyDomainV0,
    LedgerChunkRefV0, MaterialConservationV0, ParcelBirthV0, ParcelId,
    ParcelIdentityConservationV0, ParcelLedgerDescriptorV0, ParcelLedgerV0,
    ParcelRetirementCauseV0, ParcelRetirementV0Detail, PhaseSplit, PlateId, PlateRecordV0,
    PlateRetirementDetailV0, PlateRetirementV0, ProjectionDescriptorV0, ProvenanceFraction,
    ReachId, RiverReachV0, SedimentDepositionV0, SedimentEpochBalanceV0, SedimentMassBalanceV0,
    SedimentSinkV0, SedimentTransferRecord, TectonicEventV0, TransportPathway, WaterfallCause,
    WaterfallObservation,
};
pub use laws::{
    ConstantV1, LawSetV1, SourceClassV1, UniverseConstantsV1, UniverseError, UniverseV1, UnitV1,
};
pub use params::{FidelityPreset, ParamError, PlanetParameterProfile, ScenarioPreset, SeedPath};
pub use spatial::{
    CoordinateSystem, DatumDescriptor, DiscretizationDescriptor, DomainDeclaration, ElementId,
    NeighborLink, SpatialDiscretization, SpatialDomainKind, SpatialError,
};

/// Units used across Planet Foundry fields. Stated once, referenced by name
/// in every `FieldDescriptor::units` and validated against this table.
pub mod units {
    pub const MY: &str = "My"; // millions of years
    pub const M: &str = "m"; // meters
    pub const M_PER_MY: &str = "m/My";
    pub const KM: &str = "km";
    pub const MW_M2: &str = "mW/m2"; // surface heat flow
    pub const KELVIN: &str = "K"; // temperature anomaly uses K
    pub const DEG_MY: &str = "deg/My"; // plate rotation rate
    pub const M2: &str = "m2"; // represented surface area (parcel material)
    pub const M3_S: &str = "m3/s"; // discharge proxy
    pub const KG: &str = "kg"; // sediment mass
    pub const ID: &str = "id"; // integer identifiers (u32 fields)
    pub const CODE: &str = "code"; // declared enum codes (u8 fields)
    pub const FRACTION: &str = "1"; // dimensionless [0, 1] potentials
}

#[derive(Debug, thiserror::Error)]
pub enum ContractError {
    #[error("nonfinite value rejected at {stage}: field={field} index={index}")]
    NonFinite {
        stage: &'static str,
        field: String,
        index: usize,
    },
    #[error("serialization failed at {0}")]
    Serialization(&'static str, #[source] serde_json::Error),
    #[error("spatial domain violation: {0}")]
    DomainViolation(String),
    #[error("blob length mismatch for field={field}: expected {expected} bytes, found {found}")]
    BlobLength {
        field: String,
        expected: usize,
        found: usize,
    },
    #[error("blob hash mismatch for field={field}: recorded {recorded}, recomputed {recomputed}")]
    BlobHashMismatch {
        field: String,
        recorded: String,
        recomputed: String,
    },
    #[error("digest mismatch for {artifact}: recorded {recorded}, recomputed {recomputed}")]
    DigestMismatch {
        artifact: String,
        recorded: String,
        recomputed: String,
    },
    #[error("publish refused: target {target} already exists")]
    PublishTargetExists { target: String },
    #[error("sidecar {logical_name} ({path}): {why}")]
    Sidecar {
        logical_name: String,
        path: String,
        why: String,
    },
    #[error("io at {stage}: {source}")]
    Io {
        stage: &'static str,
        #[source]
        source: std::io::Error,
    },
}

/// Opt-in native excavation transactions.
pub mod excavation;

pub mod mining_history;

pub mod subunits;

/// Shared read-only rock-volume ABI; no automatic activation.
pub mod volume;

pub mod player_inventory;

pub mod resource_admission;
