//! The PREVIEW-tier substrate chain (ADR 032, proposed 2026-09-13):
//!
//! ```text
//! seeded plates + Euler poles (plates.rs, built)
//!   -> Eulerian plate field, ridge birth, subduction, collision (kinematic.rs)
//!   -> per-cell boundary classification (boundaries.rs kernels, built)
//!   -> thickening prior from convergence (uplift.rs)
//!   -> crust thickness -> isostatic elevation by the profile's UniverseV1
//!      isostasy law (isostasy.rs; WS-CONST: constants read through laws.rs)
//!      + oceanic thermal subsidence (GDH1) + water loading
//!   -> route (routing.rs, built) -> erosion (erosion.rs) -> sediment ledger (sediment.rs)
//!   -> hydrology (hydrology.rs, built) -> resource potential (here)
//!   -> publish pf.terrain/v0, pf.hydrology/v0, pf.resource/v0, BackendId tier Preview
//! ```
//!
//! Backend identity: `kinematic-preview-v0`, `Tier::Preview`. Everything
//! this backend does not do (solved mantle, melt, element mass, parcel
//! lineage, S5a provenance records, flexure, waterfall detection) is
//! listed in `TerrainDomainV0.approximations` and `assumptions`, so the
//! artifact carries its own limits (ADR 032 usage rule).
//!
//! Determinism: same `SeedPath` + same params => identical digests
//! (`tests/determinism.rs`). No randomness beyond `plates.rs`.
//!
//! WS-MODULE (lane SM, 2026-10-12): the former single loop of
//! `simulate_with` is split, code unchanged, into `PreviewEnv::new` (what
//! the recipe determines), `PreviewRun::init` (stage A), `PreviewRun::step`
//! (one loop iteration) and `PreviewRun::finish` (the tail). Batch
//! generation calls them in order; the genesis SimModule in
//! `preview_module.rs` drives the same calls (ADR 039 §4).

use pf_contracts::artifact::{digest_canonical, digest_params, BackendId, Tier};
use pf_contracts::domain::{
    BoundaryKindV0, BoundaryRecordV0, HydrologyCellFieldsV0, HydrologyDomainV0, PlateRecordV0,
};
use pf_contracts::params::{
    Compensation, ContinentalBodiesBlock, FidelityPreset, MantleBlock, OrogenyBlock,
    PlanetParameterProfile, RainfallMode, ScenarioPreset, SedimentLawsBlock, SeedPath,
    SurfaceBlock, TwoColumnIsostasyBlock,
};
use pf_contracts::terrain::{
    AssumptionSourceV0, AssumptionV0, CrustKindCodeV0, CrustVolumeLedgerV0, ErosionRecordV0,
    ResourceFamilyKindV0, ResourceFamilyV0, ResourcePotentialDomainV0, RockFamilyCodeV0,
    TerrainCodesV0, TerrainDomainV0, TerrainSummaryV0, RESOURCE_PREVIEW_USAGE_LIMITS,
};
use pf_contracts::{units, ArtifactKind, ContractError};
use pf_geodynamics::boundaries::BoundaryTable;
use pf_geodynamics::carry::CarryLaw;
use pf_contracts::laws::{
    CrustSeedingLawV1, ErosionLawV1, IsostasyLawV1, LandDepositionLawV1, SeaLevelInputV1,
    CollisionSedimentLawV1, PlateCarryLawV1, SedimentPlacementLawV1, TectonicThickeningLawV1,
    UniverseV1, HotspotVolcanismLawV1, ContinentalRiftLawV1, PlateReorganizationLawV1,
    RiftFlankLawV1, SubductionTrenchLawV1, PassiveMarginLawV1, PlateSutureLawV1,
};
use pf_contracts::laws::OrogenPlateauLawV1; // SEED-RANGES (lane RG)
use pf_contracts::laws::ClosedBasinLawV1; // INLAND-BASINS (lane INLAND-BASINS)
use pf_contracts::laws::PlateBoundaryLawV1; // SEED-TRANSFORMS (lane TF)
use pf_geodynamics::ranges::{
    plateau_spill, suture_indent_pass, IndentParams, IndentStepReport, Indentation, PlateauParams,
    PlateauReport,
};
use pf_geodynamics::conserve::{
    split_parcels, suture_pass, write_back, write_back_thickness_weighted, ParcelLayer, ParcelParams, SutureState,
};
// SEED-STRATIGRAPHY (lane ST): parcel-stacks-v1.
use pf_contracts::laws::StratigraphyLawV1;
use pf_geodynamics::stratigraphy::{CellTags, StrataParams, StrataStack, StrataStepLog};
use pf_geodynamics::isostasy::{
    single_column_airy_elevation_m, two_column_field_with, ColumnV0, SingleColumnConstants,
    TwoColumnConstants,
};
use pf_geodynamics::kinematic::{
    advance_with_extras, advance_with_laws, advance_with_parcels_transform, classify_cells_law, nearest_cell, AdvanceLedger,
    ContinentalBodiesReport, CrustState, ExtrasLaws, KinematicError, NO_SOURCE,
};
use pf_geodynamics::hotspots::{seed_plumes, Plume, PlumeParams};
use pf_geodynamics::rift::{split_plate, RiftParams};
use pf_geodynamics::reorg::{
    active_plumes, plume_split_law, reorg_pass_law, reorg_schedule, seed_plumes_renewal, ReorgCtx, ReorgLaw,
    ReorgParams,
    ReorgState, ReorgStepReport,
};
use pf_geodynamics::rift_flank::{flank_term, stretch_update, RiftFlankParams};
use pf_geodynamics::trench::{trench_term, TrenchParams};
use pf_contracts::laws::FloodBasaltLawV1;
use pf_geodynamics::flood_basalt::{
    cell_centroids, swell_term, FloodBasaltParams, LipParams, SwellParams,
};
use pf_geodynamics::orogen::{
    subgrid_step, subgrid_tau_my, BeltParams, CollisionSedimentLaw, SubgridInputs,
};
use pf_geodynamics::plates::{seed_plates, PlatesError};
use pf_contracts::laws::CratonSeedingLawV1; // SEED-IRREGULAR-SEEDS (lane IS)
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;
use pf_surface::erosion::{
    erode_deposit_step, DiffusionStencil, ErosionDepositionParams, ErosionError, ErosionParams,
    SurfaceLaws,
};
use pf_surface::hydrology::{assemble_hydrology_with_cells, HydrologyError};
use pf_surface::relief::{subgrid_lambda, NonlinearSubgridParams};
use pf_surface::uplands::CriticalSlopeParams;
use pf_surface::erosion::retain::{erode_retain_step_with_lakes, RetentionParams};
use pf_contracts::laws::SedimentRetentionLawV1;
use pf_contracts::laws::ShelfGradingLawV1; // SHELF law-set entry (lane CLIMATE-GENESIS-2)
use pf_contracts::laws::ClimateSubgridCellMeanV1; // I11 (lane CLIMATE-GENESIS-2)
use pf_contracts::laws::HillslopeThresholdLawV1;
use pf_surface::routing::{route, Routing, RoutingError};
// CLIMATE-GENESIS (lane CG): `epoch-climate-runoff-v1`.
use pf_contracts::params::ClimateForcingBlock;
use pf_surface::climate::seed::genesis::{genesis_runoff_step, GenesisClimateParams, GenesisClimateState};
use pf_surface::sediment::{
    CompensatedColumn, DepositionReport, SedimentError, SedimentLedger, ShelfParams,
};
use std::collections::BTreeMap;
use std::path::Path;
use std::time::Instant;

use crate::receipt::{f32_field, receipt, u32_field, u8_field, FieldBlob};

pub const BACKEND_NAME: &str = "kinematic-preview-v0";
pub const BACKEND_VERSION: &str = "0.2.0";
pub const SURFACE_DOMAIN: &str = "surface.primary";

/// Water density for the water-loading factor (kg/m^3).
pub const WATER_DENSITY_KG_M3: f64 = 1030.0;
/// Stein and Stein 1992 GDH1 plate-cooling depth-age constants (m, My).
pub const GDH1_RIDGE_DEPTH_M: f64 = 2600.0;
pub const GDH1_SQRT_COEF_M_PER_SQRT_MY: f64 = 365.0;
pub const GDH1_ASYMPTOTE_M: f64 = 5651.0;
pub const GDH1_EXP_COEF_M: f64 = 2473.0;
pub const GDH1_EXP_RATE_PER_MY: f64 = 0.0278;
pub const GDH1_SWITCH_MY: f64 = 20.0;
/// Sediment thickness at which rock_family reports SedimentaryCover (m).
pub const SEDIMENTARY_COVER_THRESHOLD_M: f64 = 200.0;
/// Erosion exposure saturation (m of cumulative bedrock removal).
pub const EXPOSURE_SATURATION_M: f64 = 5000.0;
/// Sediment ledger numerical budget.
pub const SEDIMENT_REL_TOL: f64 = 1.0e-9;
pub const SEDIMENT_ABS_FLOOR_KG: f64 = 1.0e-3;
/// Declared explicit-diffusion stability limit.
pub const DIFFUSION_STABILITY_LIMIT: f64 = 0.5;
/// Stream-power exponents (n = 1 required by the implicit scheme).
pub const STREAM_POWER_M: f64 = 0.5;
pub const STREAM_POWER_N: f64 = 1.0;
/// Unit conversions declared as assumptions (params carry per-year values).
pub const YEARS_PER_MY: f64 = 1.0e6;

#[derive(Debug)]
pub enum PreviewError {
    Params(String),
    Plates(PlatesError),
    Kinematic(KinematicError),
    Routing(RoutingError),
    Erosion(ErosionError),
    Sediment(SedimentError),
    Hydrology(HydrologyError),
    Contract(ContractError),
    Domain(String),
    Grid(String),
}

impl std::fmt::Display for PreviewError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PreviewError::Params(w) => write!(f, "params: {w}"),
            PreviewError::Plates(e) => write!(f, "plates: {e}"),
            PreviewError::Kinematic(e) => write!(f, "kinematic: {e}"),
            PreviewError::Routing(e) => write!(f, "routing: {e}"),
            PreviewError::Erosion(e) => write!(f, "erosion: {e}"),
            PreviewError::Sediment(e) => write!(f, "sediment: {e}"),
            PreviewError::Hydrology(e) => write!(f, "hydrology: {e}"),
            PreviewError::Contract(e) => write!(f, "contract: {e}"),
            PreviewError::Domain(w) => write!(f, "domain: {w}"),
            PreviewError::Grid(w) => write!(f, "grid: {w}"),
        }
    }
}

impl std::error::Error for PreviewError {}

macro_rules! from_err {
    ($t:ty, $v:ident) => {
        impl From<$t> for PreviewError {
            fn from(e: $t) -> Self {
                PreviewError::$v(e)
            }
        }
    };
}
from_err!(PlatesError, Plates);
from_err!(KinematicError, Kinematic);
from_err!(RoutingError, Routing);
from_err!(ErosionError, Erosion);
from_err!(SedimentError, Sediment);
from_err!(HydrologyError, Hydrology);
from_err!(ContractError, Contract);

/// Run configuration. `grid` overrides the fidelity ladder with an
/// off-ladder mini grid (recorded as resolution_level u32::MAX).
#[derive(Debug, Clone)]
pub struct PreviewConfig {
    pub scenario: ScenarioPreset,
    pub fidelity: FidelityPreset,
    pub grid: Option<(usize, usize)>,
    pub seed: u64,
    pub duration_my: Option<f64>,
    pub step_my: Option<f64>,
}

impl PreviewConfig {
    pub fn seed_path(&self) -> SeedPath {
        SeedPath {
            seed: self.seed,
            derivation: "kinematic-preview".into(),
            salt: self.scenario.name.clone(),
        }
    }

    /// The resolved parameter profile: scenario + fidelity, compensation
    /// forced to Airy (declared: flexure deferred at this tier).
    pub fn profile(&self) -> Result<PlanetParameterProfile, PreviewError> {
        let mut mantle = MantleBlock::default();
        mantle.compensation = Compensation::Airy;
        let mut p = self.scenario.profile(self.fidelity, mantle);
        p.mantle.compensation = Compensation::Airy;
        if let Some(d) = self.duration_my {
            p.tectonics.duration_my = d;
        }
        if let Some(s) = self.step_my {
            p.tectonics.step_my = s;
        }
        check_profile(&p)?;
        Ok(p)
    }

    pub fn grid(&self, p: &PlanetParameterProfile) -> Result<GridSpec, PreviewError> {
        grid_for(self.grid, p)
    }
}

/// The run grid: the off-ladder override, or the profile's fidelity level.
/// WS-CONST: the planet radius is read through the universe constants.
/// (WS-MODULE: lifted out of `PreviewConfig::grid` unchanged so a resolved
/// recipe can build its grid without a `PreviewConfig`.)
pub fn grid_for(
    grid_override: Option<(usize, usize)>,
    p: &PlanetParameterProfile,
) -> Result<GridSpec, PreviewError> {
    let radius_m = universe_of(p)?.constants.planet_radius_m.value;
    match grid_override {
        Some((nlat, nlon)) => GridSpec::new(nlat, nlon, radius_m)
            .map_err(|e| PreviewError::Grid(e.to_string())),
        None => GridSpec::level(p.fidelity.grid_level, radius_m)
            .map_err(|e| PreviewError::Grid(e.to_string())),
    }
}

/// OPT-IN physics flags for the preview chain (lanes F and K2,
/// 2026-10-11). Kept OUT of `PreviewConfig` so every existing config
/// literal and `simulate()` call is unchanged; `simulate(c)` is exactly
/// `simulate_with(c, &PreviewOptions::default())`, and the default options
/// leave the resolved profile byte-identical (both blocks are skipped on
/// write when absent). When a flag is set it is written into the resolved
/// `PlanetParameterProfile`, so it moves `params_digest` and is recorded in
/// every published artifact.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PreviewOptions {
    /// Q29 opt-in native evidence; absent preserves historical recipes.
    pub native_evidence: Option<pf_contracts::params::NativeProducerEvidenceV1>,
    /// Two-column isostasy + volume-primary sea level (lane F model).
    pub two_column_isostasy: Option<TwoColumnIsostasyBlock>,
    /// Continental crust bodies with a rifted-margin taper (lane K2).
    pub continental_bodies: Option<ContinentalBodiesBlock>,
    /// DIAGNOSTIC counterfactuals only (e.g. erosion off to separate
    /// isostasy from the erosion budget). Replaces the scenario's surface
    /// block in the resolved profile, so it is recorded in params_digest;
    /// never a proposed setting.
    pub surface_override: Option<SurfaceBlock>,
    /// SEED-OROGEN (lane OG): rate-limited belts and/or subgrid stream
    /// power. Written into the resolved profile (`orogeny` block), so it
    /// moves params_digest and the law set.
    pub orogeny: Option<OrogenyBlock>,
    /// SEED-EXTRAS D7 (lane SD): `erosion-deposition-v1` and/or
    /// `shelf-wave-base-v1`. Written into the resolved profile
    /// (`sediment_laws` block), so it moves params_digest and the law set.
    pub sediment_laws: Option<SedimentLawsBlock>,
    /// CLIMATE-GENESIS (lane CG): `epoch-climate-runoff-v1` (climate runoff
    /// into the genesis fluvial kernels). Written into the resolved profile
    /// (`climate_forcing` block), so it moves params_digest and the law set.
    pub climate_forcing: Option<ClimateForcingBlock>,
}

impl PreviewOptions {
    /// Legacy chain (no flags).
    pub fn legacy() -> Self {
        Self::default()
    }
    /// Two-column isostasy only (Earth-like v0 constants).
    pub fn two_column_v0() -> Self {
        Self {
            two_column_isostasy: Some(TwoColumnIsostasyBlock::earth_v0()),
            continental_bodies: None,
            surface_override: None,
            orogeny: None,
            sediment_laws: None,
            climate_forcing: None,
            native_evidence: None,
        }
    }
    /// Two-column isostasy plus continental bodies (CRUST-SPIKE v0).
    pub fn crust_bodies_v0() -> Self {
        Self {
            two_column_isostasy: Some(TwoColumnIsostasyBlock::earth_v0()),
            continental_bodies: Some(ContinentalBodiesBlock::v0()),
            surface_override: None,
            orogeny: None,
            sediment_laws: None,
            climate_forcing: None,
            native_evidence: None,
        }
    }
    /// SEED-OROGEN candidate profile: two-column isostasy, continental
    /// bodies, `rate-limited-belt-v1` and `stream-power-subgrid-v1`.
    pub fn orogen_v1() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-CARRY-CC candidate profile (lane CC): `orogen_v1` with the
    /// plate carry `front-rescue-carry-v1`.
    pub fn orogen_carry_v1() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::carry_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEDIMENT-TOWERS (lane ST) corrected carry profile: `orogen_carry_v1`
    /// with `collision_sediment = accreted-to-belt-crust-v1`.
    pub fn orogen_carry_v2() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::carry_v2()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEDIMENT-TOWERS (lane ST) corrected profile: `orogen_v1` with
    /// `collision_sediment = accreted-to-belt-crust-v1`.
    pub fn orogen_v2() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::v2()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-EXTRAS candidate profile (lane EX): `orogen_carry_v2` plus
    /// hotspot volcanism `fixed-plume-chains-v1`, continental rifting
    /// `stretched-rift-v1` and plate splits `plume-head-continental-split-v1`.
    pub fn extras_v1() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::extras_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-RIFT-SHOULDERS candidate profile (lane RS, 2026-10-13):
    /// `extras_v1` plus `rift_flank = shoulder-uplift-thermal-v1` (rift
    /// thermal anomaly and subsidence, flexural shoulders).
    pub fn extras_v2() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::extras_v2()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-TRENCH candidate profile (lane TN, 2026-10-13): `extras_v2`
    /// plus `subduction_trench = slab-flexure-trench-v1` (flexural trench
    /// and outer rise of the downgoing oceanic plate).
    pub fn extras_v3() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::extras_v3()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-FLOODBASALT candidate profile (lane FB, 2026-10-13):
    /// `extras_v3` plus `flood_basalt = plume-head-lip-swell-v1` (plume-head
    /// flood-basalt plateaus and mantle-fixed plume swells).
    pub fn extras_v3_flood() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::extras_v3_flood()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-CARRY-CONSERVE candidate profile (lane CV, 2026-10-13):
    /// `extras_v3` with the plate carry `continental-parcel-carry-v1` (in
    /// place of `front-rescue-carry-v1`) and `plate_suture =
    /// contact-length-pole-share-v1` (PREREG.txt in
    /// `unpublished-design-record`).
    pub fn conserve_v1() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::conserve_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-IRREGULAR-SEEDS + SEED-FREEBOARD start height candidate
    /// profile (lane IS, 2026-10-13): `conserve_v1` with `craton_seeding =
    /// irregular-cratons-v1` (old irregular cratons, worn sutures, rough
    /// basement, interiors starting near +0.4 km; PREREG.txt in
    /// `unpublished-design-record`).
    pub fn oldstart_v1() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            ..Self::conserve_v1()
        }
    }
    /// SEED-PLATE-REORG (lane PR, 2026-10-13): `oldstart_v1` plus lane SR's
    /// `sediment_retention = continental-retention-base-level-v1` at its v1
    /// values (no new constant; until now a driver-only composition).
    pub fn oldstart_retain_v1() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::conserve_retain_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-PLATE-REORG candidate profile (lane PR, 2026-10-13):
    /// `oldstart_retain_v1` plus `plate_reorganization =
    /// plate-reorg-scripted-v1` (scripted Wilson cycle; PREREG.txt in
    /// `unpublished-design-record`).
    pub fn cycle_v1() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::cycle_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-MAGMATISM candidate profile (lane MG, 2026-10-13): `cycle_v1`
    /// plus `magmatism = slab-arc-backarc-llsvp-v1` (slab-depth arcs, island
    /// arcs maturing into continental crust, back-arc extension, LLSVP-edge
    /// drifting plumes, volcano points, magmatic resource families;
    /// PREREG.txt in `unpublished-design-record`).
    pub fn cycle_magma_v1() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::cycle_magma_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-BREAKUP-CUT-V2 candidate profile (lane BC, 2026-10-13):
    /// `oldstart_retain_v1` plus `plate_reorganization =
    /// plate-reorg-scripted-v2` (breakup cut across continents; PREREG.txt in
    /// `unpublished-design-record`).
    pub fn cycle_v2() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::cycle_v2()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-RANGES candidate profile (lane RG, 2026-10-13): `cycle_v1` with
    /// `plate_suture = contact-length-indentation-taper-v1` (post-collision
    /// indentation) and `orogen_plateau = plateau-cap-spill-v1` (plateau
    /// growth; belt delamination backstop raised to 100 km). PREREG.txt in
    /// `unpublished-design-record`.
    pub fn cycle_ranges_v1() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::cycle_ranges_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-RANGES composition C (lane RG, 2026-10-13): `cycle_ranges_v1`
    /// with lane UP's erosion (`stream-power-subgrid-nonlinear-v1` n = 2 at
    /// RL v1 values plus `critical-slope-envelope-v1` v1); no new constant.
    pub fn cycle_ranges_uplands_v1() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::cycle_ranges_uplands_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-WORLD-COMPOSE candidate profile (lane CP, 2026-10-13): every
    /// approved piece in one opt-in profile: oldstart + retention base,
    /// `plate-reorg-scripted-v3` (v2 cut, symmetric ridges, age-driven
    /// subduction initiation), indentation, `plateau-cap-spill-v2` (inside
    /// the belt step, crust plus sediment, foothill taper; 70 km backstop),
    /// n = 2 erosion + envelope, `slab-arc-backarc-llsvp-v2` and
    /// `plume-head-lip-swell-v1` on the same LLSVP plumes. PREREG.txt in
    /// `unpublished-design-record`.
    pub fn world_v1() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::world_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-STRATIGRAPHY stage A candidate profile (lane ST, 2026-10-13):
    /// `world_v1` plus `stratigraphy = parcel-stacks-v1` (rock columns on
    /// the continental parcels; a passive record, so the geography is
    /// world_v1's). PREREG.txt in `unpublished-design-record`.
    pub fn world_strata_v1() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::world_strata_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-STRATIGRAPHY stage A.1 candidate profile (lane SA, 2026-10-13):
    /// `world_v1` plus `stratigraphy = parcel-stacks-v2` (belt-placed
    /// plutons, sub-cell shelf, t = 0 platform cover, 16-layer stacks with
    /// the unconformity-keeping merge, rock uplift and its gradient; still a
    /// passive record). PREREG.txt in `unpublished-design-record`.
    pub fn world_strata_v2() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::world_strata_v2()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-TRANSFORMS candidate profile (lane TF, 2026-10-13): `world_v1`
    /// with `plate_boundary = smoothed-normal-transform-v1` (smoothed
    /// boundary normal, transform class, no birth or subduction on
    /// transforms, no ridge birth on convergence) and
    /// `plate-reorg-scripted-v4` (slab-age and spacing limits on the
    /// age-triggered initiation). PREREG.txt in `unpublished-design-record`.
    pub fn world_v2() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::world_v2()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-OLD-OCEAN candidate profile (lane OR, 2026-10-13): `world_v2`
    /// with `plate-reorg-scripted-v5` (v4 plus passive-margin subduction
    /// initiation driven by ocean age). PREREG.txt in `unpublished-design-record`.
    pub fn world_v3() -> Self {
        Self {
            continental_bodies: Some(ContinentalBodiesBlock::cratons_v1()),
            orogeny: Some(OrogenyBlock::world_v3()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-RELIEF-LAW candidate profile (lane RL): `orogen_v1` with its
    /// erosion law `stream-power-subgrid-v1` replaced by
    /// `stream-power-subgrid-nonlinear-v1` (n = 2, k_ref = 500 m, x_c =
    /// 100 m, pre-registered). Belt law, bodies and two-column unchanged.
    pub fn orogen_relief_v1() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::relief_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-UPLANDS candidate profile (lane UP, 2026-10-13): `extras_v2`
    /// with the erosion law `stream-power-subgrid-nonlinear-v1` (lane RL v1)
    /// plus `hillslope_threshold = critical-slope-envelope-v1` (S_c 0.04,
    /// knickpoint celerity 70 km/My; PREREG.txt).
    pub fn uplands_v1() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::uplands_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-SEDIMENT-RETAIN candidate profile (lane SR, 2026-10-13):
    /// `extras_v3` plus `sediment_retention = continental-retention-base-level-v1`
    /// (H 1000 m, p 2, S_a 1e-4, wave base 100 m; PREREG.txt).
    pub fn retain_v1() -> Self {
        Self {
            orogeny: Some(OrogenyBlock::retain_v1()),
            ..Self::crust_bodies_v0()
        }
    }
    /// SEED-EXTRAS D7 candidate profile (lane SD): two-column isostasy,
    /// continental bodies, `erosion-deposition-v1` and `shelf-wave-base-v1`.
    pub fn sediment_v1() -> Self {
        Self {
            sediment_laws: Some(SedimentLawsBlock::v1()),
            ..Self::crust_bodies_v0()
        }
    }
}

impl PreviewConfig {
    /// `profile()` plus the opt-in flags, re-validated.
    pub fn profile_with(
        &self,
        options: &PreviewOptions,
    ) -> Result<PlanetParameterProfile, PreviewError> {
        let mut p = self.profile()?;
        p.two_column_isostasy = options.two_column_isostasy.clone();
        p.continental_bodies = options.continental_bodies.clone();
        p.orogeny = options.orogeny.clone();
        p.sediment_laws = options.sediment_laws.clone();
        p.climate_forcing = options.climate_forcing.clone();
        p.native_evidence = options.native_evidence.clone();
        if let Some(s) = &options.surface_override {
            p.surface = s.clone();
        }
        check_profile(&p)?;
        Ok(p)
    }
}

/// Validate a resolved profile for this backend: the contract's own dial
/// validation, then the backend's declared capability (uniform rainfall
/// only; LatitudeBands is refused, never substituted).
pub fn check_profile(p: &PlanetParameterProfile) -> Result<(), PreviewError> {
    p.validate()
        .map_err(|e| PreviewError::Params(e.to_string()))?;
    // WS-CONST: the law set and constants must resolve and validate
    // (typed OutOfEnvelope / Unsupported / Inconsistent / PinMismatch).
    universe_of(p)?;
    if p.surface.rainfall_mode != RainfallMode::Uniform {
        return Err(PreviewError::Params(
            "kinematic-preview-v0 implements RainfallMode::Uniform only (LatitudeBands is not built; refusing rather than substituting)".into(),
        ));
    }
    Ok(())
}

/// Per-step diagnostics (printed by the CLI; not part of any digest).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepLog {
    pub t_my: f64,
    pub ridge_cells: u64,
    pub subduction_cells: u64,
    pub collision_cells: u64,
    pub advection_residual_m3: f64,
    pub eroded_m3: f64,
    pub deposited_m3: f64,
    pub exported_m3: f64,
    pub diffusion_substeps: u32,
    pub deposit_aggraded_m3: f64,
    pub max_elevation_m: f64,
    pub min_elevation_m: f64,
    pub land_fraction: f64,
    /// Two-column mode only (skipped otherwise, so the legacy CLI output is
    /// unchanged): the solved sea level on the mantle free-surface gauge
    /// after this step's re-equilibration (m), the flooded area fraction
    /// (global fill) and the water-volume closure V(S) - V (m^3).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sea_level_gauge_m: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flooded_area_fraction: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ocean_volume_residual_m3: Option<f64>,
    /// Two-column mode only: compensation depth used, deepest root, and
    /// the deepest column's [crust, sediment] thickness (m).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub compensation_depth_m: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub deepest_root_m: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub deepest_root_crust_sediment_m: Option<[f64; 2]>,
    /// Flag modes only: area-weighted 10/50/90th percentile of land
    /// elevation after this step (m), to time continental planation.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub land_elevation_p10_p50_p90_m: Option<[f64; 3]>,
    /// SEED-EXTRAS D7 (lane SD) laws only: where this step's moved material
    /// came to rest (m^3) and the solver/placement counters.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub deposition: Option<StepDeposition>,
}

/// Lane SD: the per-step sediment partition (see
/// `pf_surface::sediment::DepositionReport`).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDeposition {
    pub erosion_law: &'static str,
    pub sediment_law: &'static str,
    pub land_fluvial_m3: f64,
    pub land_hillslope_m3: f64,
    pub land_floor_m3: f64,
    pub lake_m3: f64,
    pub marine_shelf_m3: f64,
    pub marine_slope_m3: f64,
    pub marine_deep_m3: f64,
    pub global_remainder_m3: f64,
    pub delivered_to_ocean_m3: f64,
    pub marine_search_bound_hits: u64,
    pub shelf_cells_filled: u64,
    pub solver_passes: u32,
    pub solver_converged: bool,
}

impl From<&DepositionReport> for StepDeposition {
    fn from(r: &DepositionReport) -> Self {
        Self {
            erosion_law: r.erosion_law,
            sediment_law: r.sediment_law,
            land_fluvial_m3: r.land_fluvial_m3,
            land_hillslope_m3: r.land_hillslope_m3,
            land_floor_m3: r.land_floor_m3,
            lake_m3: r.lake_m3,
            marine_shelf_m3: r.marine_shelf_m3,
            marine_slope_m3: r.marine_slope_m3,
            marine_deep_m3: r.marine_deep_m3,
            global_remainder_m3: r.global_remainder_m3,
            delivered_to_ocean_m3: r.delivered_to_ocean_m3,
            marine_search_bound_hits: r.marine_search_bound_hits,
            shelf_cells_filled: r.shelf_cells_filled,
            solver_passes: r.sweeps,
            solver_converged: r.converged,
        }
    }
}

/// Lane SD: the surface laws a universe selects (legacy = both `None`).
/// The compensated column uses the isostasy law's own densities.
pub fn surface_laws_of(u: &UniverseV1) -> SurfaceLaws {
    let erosion_deposition = match &u.law_set.land_deposition {
        LandDepositionLawV1::DetachmentLimitedV0 => None,
        LandDepositionLawV1::ErosionDepositionV1 {
            deposition_coefficient_g,
        } => Some(ErosionDepositionParams {
            deposition_coefficient_g: deposition_coefficient_g.value,
            ..ErosionDepositionParams::v1()
        }),
    };
    let shelf = match &u.law_set.sediment_placement {
        SedimentPlacementLawV1::LocalFillGlobalLayerV0 => None,
        SedimentPlacementLawV1::ShelfWaveBaseV1 {
            shelf_wave_base_depth_m,
        } => Some(ShelfParams {
            wave_base_depth_m: shelf_wave_base_depth_m.value,
            ..ShelfParams::v1(CompensatedColumn {
                mantle_density_kg_m3: u.constants.mantle_density_kg_m3.value,
                sediment_density_kg_m3: u.sediment_density_kg_m3(),
                water_density_kg_m3: u.constants.seawater_density_kg_m3.value,
            })
        }),
    };
    SurfaceLaws {
        erosion_deposition,
        shelf,
    }
}

/// Lane SD: assumptions and approximation names the sediment laws add or
/// replace (only called when one is selected; legacy lists untouched).
fn sediment_law_assumptions(
    laws: &SurfaceLaws,
    assumptions: &mut Vec<AssumptionV0>,
    approximations: &mut Vec<String>,
) {
    let pa = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source: AssumptionSourceV0::ProjectAssumption,
        note: note.into(),
    };
    if let Some(ed) = &laws.erosion_deposition {
        for a in approximations.iter_mut() {
            if a == "implicit-stream-power-n1" {
                *a = "implicit-erosion-deposition-n1-yuan-style-v1".into();
            }
        }
        assumptions.push(pa("A-DEPOSITION-G", "dimensionless deposition coefficient G: land deposition rate G Qs / (p A), p = 1 (uniform runoff, geometric area)", ed.deposition_coefficient_g, "1", "project assumption awaiting review; Yuan et al. 2019 (JGR Earth Surface 124) and Davy and Lague 2009 (JGR 114) recalled, not checked; declared before any run (lane SD FALSIFIERS.md)"));
        assumptions.push(pa("A-ED-NO-FLEXURE", "erosion-deposition stores sediment only where the graded river surface needs it; no flexural foreland subsidence (local Airy only), so foreland storage is the mountain-front alluvial wedge", 0.0, "1", "flexure is a SEED-EXTRAS item; lane SD measured land storage at about 0.4 % of eroded volume with G = 1 at the toy K"));
    }
    if let Some(sh) = &laws.shelf {
        assumptions.retain(|a| !matches!(a.id.as_str(), "A-MARINE-GLOBAL-LAYER" | "A-SEDIMENT-SPREAD-BOUND"));
        approximations.push("shelf-wave-base-v1-compensated-marine-accommodation".into());
        assumptions.push(pa("A-SHELF-WAVE-BASE", "depth below sea level to which the shelf top aggrades; cells at or above it pass sediment on; deeper cells are filled up to it outward from each river mouth (contiguous progradation)", sh.wave_base_depth_m, "m", "project assumption awaiting review; Earth shelf break ~130 m recalled, not checked; declared before any run"));
        assumptions.push(pa("A-COMPENSATED-ACCOMMODATION", "standing-water accommodation includes the local-Airy response to the sediment load (thickness = level change x (rho_m - rho_w)/(rho_m - rho_s) under water, rho_m/(rho_m - rho_s) in air) with the isostasy law's densities; sea level held fixed within a step", sh.column.rise_per_m_submerged(), "1", "the volume-primary sea-level response to displaced water appears at the next isostatic solve"));
        assumptions.push(pa("A-MARINE-SEARCH-BOUND", "computational bound on ocean cells settled per river-mouth placement; what it cannot place is the declared global remainder (uniform over the ocean reservoir, measured per step)", sh.search_max_cells as f64, "cells", "not a physical capacity; hits and the remainder volume are in the step log"));
        assumptions.push(pa("A-NO-DEEP-SEA-BYPASS", "no grain-size partition, turbidity currents or pelagic rain: deep-water deposition happens only where a shelf front advances into it", 0.0, "1", "a fines-bypass law would move part of the slope class to the abyssal plain"));
    }
}

/// Area-weighted percentiles of the land (>= sea) elevations; None when no
/// land. Diagnostic only.
fn land_percentiles(elevation: &[f64], areas: &[f64], sea: f64, qs: [f64; 3]) -> Option<[f64; 3]> {
    let mut land: Vec<(f64, f64)> = elevation
        .iter()
        .zip(areas)
        .filter(|(e, _)| **e >= sea)
        .map(|(e, a)| (*e, *a))
        .collect();
    if land.is_empty() {
        return None;
    }
    land.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f64 = land.iter().map(|x| x.1).sum();
    let mut out = [0.0; 3];
    for (k, q) in qs.iter().enumerate() {
        let mut acc = 0.0;
        out[k] = land.last().unwrap().0;
        for &(e, a) in &land {
            acc += a;
            if acc >= q * total {
                out[k] = e;
                break;
            }
        }
    }
    Some(out)
}

/// Everything the chain produced, in memory, before publication.
pub struct PreviewOutputs {
    /// Source-bound native evidence, published as a terrain sidecar only when selected.
    pub native_evidence: Option<serde_json::Value>,
    pub exchange_epochs: Vec<pf_contracts::domain::WaterExchangeEpochV0>,
    pub grid: GridSpec,
    pub params: PlanetParameterProfile,
    pub seed_path: SeedPath,
    pub plates: Vec<PlateRecordV0>,
    pub state: CrustState,
    pub elevation_m: Vec<f64>,
    pub uplift_rate_m_per_my: Vec<f64>,
    pub routing: Routing,
    pub terrain: TerrainDomainV0,
    pub hydrology: HydrologyDomainV0,
    pub cell_reach: Vec<Option<pf_contracts::domain::ReachId>>,
    pub resource: ResourcePotentialDomainV0,
    pub potentials: Vec<(ResourceFamilyKindV0, Vec<f32>)>,
    pub steps: Vec<StepLog>,
    /// Elevation at t = 0 from the initial crust (diagnostic; not published).
    pub elevation_t0_m: Vec<f64>,
    /// Two-column mode: the t = 0 sea-level solve (gauge m, flooded frac).
    pub sea_level_t0: Option<SeaLevelSolve>,
    /// Continental-bodies mode: seeding bookkeeping.
    pub continental_bodies: Option<ContinentalBodiesReport>,
    /// SEED-OROGEN laws only: per-step belt/arc/spreading/subgrid ledger and
    /// the final subgrid relief field.
    pub orogen: Option<OrogenRunState>,
    /// SEED-MAGMATISM (lane MG): the derived volcano record (law selected
    /// only; published as a resource sidecar).
    pub volcanoes: Option<pf_contracts::terrain::VolcanoFeaturesDomainV0>,
    /// SEED-STRATIGRAPHY (lane ST): the rock column record (law selected
    /// only; published as a terrain sidecar).
    pub stratigraphy: Option<pf_contracts::terrain::StratigraphyDomainV0>,
}

/// SEED-OROGEN (lane OG): per-step bookkeeping of the new laws (volumes in
/// m^3). Present only when a new law is selected.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct OrogenStepLog {
    pub t_my: f64,
    pub belt_in_m3: f64,
    pub collision_overlap_consumed_m3: f64,
    pub collision_supply_edges: u64,
    pub belt_placed_m3: f64,
    pub belt_sediment_m3: f64,
    pub belt_delaminated_m3: f64,
    pub spread_substeps: u32,
    pub spread_moved_m3: f64,
    pub spread_residual_m3: f64,
    pub spread_max_kappa_m2_per_my: f64,
    pub arc_created_m3: f64,
    pub arc_placed_m3: f64,
    pub max_thickness_m: f64,
    /// Subgrid erosion: sum of target cell-mean lowering * area, the realised
    /// removal (bedrock + sediment) on the same cells, the largest per-cell
    /// |realised - target| (m) over cells with a target, the count of
    /// targets capped at the receiver level, and land relief statistics.
    pub subgrid_target_m3: f64,
    pub subgrid_realised_m3: f64,
    pub subgrid_max_abs_dev_m: f64,
    pub subgrid_capped_cells: u64,
    pub subgrid_trunk_raised_cells: u64,
    pub subgrid_mean_land_relief_m: f64,
    pub subgrid_max_relief_m: f64,
    /// SEED-CARRY-CC (lane CC): carry bookkeeping under
    /// `front-rescue-carry-v1` only (skipped when absent, so every other
    /// law's log serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub carry: Option<CarryStepLog>,
    /// SEED-EXTRAS (lane EX): plume / rift / split bookkeeping (skipped when
    /// absent, so every other law's log serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub extras: Option<ExtrasStepLog>,
    /// SEED-TRENCH (lane TN): `slab-flexure-trench-v1` only (skipped when
    /// absent, so every other law's log serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trench: Option<TrenchStepLog>,
    /// SEED-CARRY-CONSERVE (lane CV): `continental-parcel-carry-v1` /
    /// `contact-length-pole-share-v1` ledger (skipped when absent, so every
    /// other law's log serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conserve: Option<ConserveStepLog>,
    /// SEED-STRATIGRAPHY (lane ST): `parcel-stacks-v1` ledger (skipped when
    /// absent, so every other law's log serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stratigraphy: Option<StratigraphyStepLog>,
}

/// SEED-STRATIGRAPHY (lane ST, 2026-10-13): the per-step ledger of
/// `parcel-stacks-v1` (volumes m^3; per-lithology arrays in code order,
/// deposits per class in `pf_geodynamics::stratigraphy::DEPOSIT_CLASSES`
/// order). Closure: after = before + every applied operation + births
/// (`closure_rel`), and after = the parcel column (`column_mismatch_rel`);
/// `before` is at write-back entry, after this step's consumption.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct StratigraphyStepLog {
    pub parcels: u64,
    pub stacks: u64,
    pub before_m3: f64,
    pub after_m3: f64,
    pub column_m3: f64,
    pub consumed_m3: Vec<f64>,
    pub births_m3: Vec<f64>,
    pub eroded_m3: Vec<f64>,
    pub deposited_m3: Vec<f64>,
    pub volcanic_m3: f64,
    pub pluton_m3: f64,
    pub plume_basement_m3: f64,
    pub basement_add_m3: f64,
    pub thinned_m3: f64,
    pub shortening_scale_m3: f64,
    pub shortening_removed_m3: f64,
    pub basement_flow_m3: f64,
    pub basement_flow_abs_m3: f64,
    /// |basement flow| over all volume the step's operations moved.
    pub untagged_share: f64,
    pub underflow_m3: f64,
    pub metamorphosed_m3: f64,
    pub merges_gap: u64,
    pub merges_cap: u64,
    pub merges_min: u64,
    pub closure_residual_m3: f64,
    pub closure_rel: f64,
    pub column_mismatch_rel: f64,
    pub entry_mismatch_max_rel: f64,
    pub mean_cover_layers: f64,
    pub stacks_at_cap: u64,
    pub basement_below_min: u64,
    /// SEED-STRATIGRAPHY stage A.1 (lane SA, `parcel-stacks-v2` only;
    /// skipped when absent, so every v1 log serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub v2: Option<StratigraphyV2StepLog>,
}

/// SEED-STRATIGRAPHY stage A.1 (lane SA): the `parcel-stacks-v2` additions
/// to the step ledger (volumes m^3): the non-placed shortening-lineage
/// change written to the basement (signed; in `closure`), warm shelf-facies
/// deposits (P3 denominator), capacity merges across an unconformity, and
/// the area-weighted mean |rock uplift| of the parcels.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct StratigraphyV2StepLog {
    pub belt_basement_m3: f64,
    pub warm_shelf_deposit_m3: f64,
    pub merges_unconformity: u64,
    pub mean_abs_rock_uplift_m: f64,
}

impl StratigraphyStepLog {
    fn of(l: &StrataStepLog, consumed: &[f64], v2: bool) -> Self {
        let moved = l.eroded_m3.iter().sum::<f64>()
            + l.deposited_m3.iter().sum::<f64>()
            + l.volcanic_m3
            + l.pluton_m3
            + l.plume_basement_m3
            + l.basement_add_m3
            + l.thinned_m3
            + l.shortening_scale_m3
            + l.shortening_removed_m3
            + l.basement_flow_abs_m3
            + l.belt_basement_m3.abs();
        Self {
            parcels: l.parcels,
            stacks: l.stacks,
            before_m3: l.before_m3,
            after_m3: l.after_m3,
            column_m3: l.column_m3,
            consumed_m3: if consumed.is_empty() { vec![0.0; pf_geodynamics::stratigraphy::NLITH] } else { consumed.to_vec() },
            births_m3: l.births_m3.to_vec(),
            eroded_m3: l.eroded_m3.to_vec(),
            deposited_m3: l.deposited_m3.to_vec(),
            volcanic_m3: l.volcanic_m3,
            pluton_m3: l.pluton_m3,
            plume_basement_m3: l.plume_basement_m3,
            basement_add_m3: l.basement_add_m3,
            thinned_m3: l.thinned_m3,
            shortening_scale_m3: l.shortening_scale_m3,
            shortening_removed_m3: l.shortening_removed_m3,
            basement_flow_m3: l.basement_flow_m3,
            basement_flow_abs_m3: l.basement_flow_abs_m3,
            untagged_share: if moved > 0.0 { l.basement_flow_abs_m3 / moved } else { 0.0 },
            underflow_m3: l.underflow_m3,
            metamorphosed_m3: l.metamorphosed_m3,
            merges_gap: l.merges_gap,
            merges_cap: l.merges_cap,
            merges_min: l.merges_min,
            closure_residual_m3: l.closure_residual_m3,
            closure_rel: l.closure_rel,
            column_mismatch_rel: l.column_mismatch_rel,
            entry_mismatch_max_rel: l.entry_mismatch_max_rel,
            mean_cover_layers: l.mean_cover_layers,
            stacks_at_cap: l.stacks_at_cap,
            basement_below_min: l.basement_below_min,
            v2: v2.then(|| StratigraphyV2StepLog {
                belt_basement_m3: l.belt_basement_m3,
                warm_shelf_deposit_m3: l.warm_shelf_deposit_m3,
                merges_unconformity: l.merges_unconformity,
                mean_abs_rock_uplift_m: l.mean_abs_rock_uplift_m,
            }),
        }
    }
}

/// SEED-STRATIGRAPHY (lane ST): the `parcel-stacks-v1` constants of a
/// universe (None = legacy).
pub fn strata_params_of(u: &UniverseV1) -> Option<StrataParams> {
    let law = &u.law_set.stratigraphy;
    if matches!(law, StratigraphyLawV1::NoneV0) {
        return None;
    }
    let cs = law.constants();
    let v = |id: &str| cs.iter().find(|(k, _)| *k == id).map(|(_, c)| c.value).expect("declared strata constant");
    Some(StrataParams {
        max_layers: v("strata_max_layers").round() as usize,
        merge_gap_my: v("strata_merge_gap_my"),
        unconformity_gap_my: v("strata_unconformity_gap_my"),
        min_layer_m: v("strata_min_layer_m"),
        trace_m: v("strata_trace_m"),
        arc_extrusive_fraction: v("strata_arc_extrusive_fraction"),
        plume_extrusive_fraction: v("strata_plume_extrusive_fraction"),
        shortening_pluton_fraction: v("strata_shortening_pluton_fraction"),
        pluton_roof_m: v("strata_pluton_roof_m"),
        metamorphic_depth_m: v("strata_metamorphic_depth_m"),
        high_grade_depth_m: v("strata_high_grade_depth_m"),
        shelf_depth_m: v("strata_shelf_depth_m"),
        carbonate_max_paleolat_deg: v("strata_carbonate_max_paleolat_deg"),
        carbonate_max_supply_m_per_my: v("strata_carbonate_max_supply_m_per_my"),
        arid_lat_min_deg: v("strata_arid_lat_min_deg"),
        arid_lat_max_deg: v("strata_arid_lat_max_deg"),
        humid_lat_max_deg: v("strata_humid_lat_max_deg"),
        humid_midlat_min_deg: v("strata_humid_midlat_min_deg"),
        humid_midlat_max_deg: v("strata_humid_midlat_max_deg"),
        coal_max_elevation_m: v("strata_coal_max_elevation_m"),
        conglomerate_relief_m: v("strata_conglomerate_relief_m"),
        river_mouth_area_m2: v("strata_river_mouth_area_m2"),
        publish_depth_m: v("strata_publish_depth_m"),
        basement_min_m: v("strata_basement_min_m"),
        // SEED-STRATIGRAPHY stage A.1 (lane SA): parcel-stacks-v2.
        v2: matches!(law, StratigraphyLawV1::ParcelStacksV2 { .. }).then(|| pf_geodynamics::stratigraphy::StrataV2Params {
            shelf_width_m: v("strata_shelf_width_m"),
            slope_foot_depth_m: v("strata_slope_foot_depth_m"),
            shelf_majority_fraction: v("strata_shelf_majority_fraction"),
            platform_cover_max_m: v("strata_platform_cover_max_m"),
            platform_full_elev_m: v("strata_platform_full_elev_m"),
            platform_zero_elev_m: v("strata_platform_zero_elev_m"),
            platform_sandstone_fraction: v("strata_platform_sandstone_fraction"),
            platform_carbonate_fraction: v("strata_platform_carbonate_fraction"),
            platform_age_base_my: v("strata_platform_age_base_my"),
            platform_age_top_my: v("strata_platform_age_top_my"),
        }),
    })
}

/// SEED-STRATIGRAPHY (lane ST): the assumptions every stratigraphy record
/// carries.
pub fn strata_assumptions(sp: &StrataParams) -> Vec<AssumptionV0> {
    let a = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source: AssumptionSourceV0::ProjectAssumption,
        note: note.into(),
    };
    let mut out = vec![
        a("A-STRATA-PARCEL-STACKS", "parcel-stacks-v1: every continental parcel carries a basement (the residual of its column) and up to max_layers - 1 cover layers written by erosion (truncation from the top), deposition by environment, arcs (volcanics on top, plutons at the basement top), plume crust (basalt on top, the rest into the basement), rift and back-arc thinning (every layer scaled), shortening (a pluton share, every layer thickened) and a basement-only residual (channel flow, plateau spill, delamination); cover below the metamorphic depth becomes metamorphic", sp.max_layers as f64, "1", "lane ST 2026-10-13, design review stratigraphy consult; pre-registered in unpublished-design-record; a passive record (writes nothing the simulation reads)"),
        a("A-STRATA-LATITUDE-CLIMATE", "no climate law runs: warm water (shelf carbonate) is the band within this latitude of the equator, the arid belt (closed-lake evaporite) and the humid bands (coal) are fixed latitude bands of the depositing cell (paleolatitude)", sp.carbonate_max_paleolat_deg, "deg", &format!("arid {}-{} deg, humid < {} and {}-{} deg; a climate law (CLIMATE-V0) would replace this fallback", sp.arid_lat_min_deg, sp.arid_lat_max_deg, sp.humid_lat_max_deg, sp.humid_midlat_min_deg, sp.humid_midlat_max_deg)),
        a("A-STRATA-METAMORPHISM", "a cover layer whose base is deeper than this becomes metamorphic (grade 1; grade 2 below the high-grade depth); plutons and oceanic crust keep their class", sp.metamorphic_depth_m, "m", &format!("high grade below {} m; greenschist onset recalled, not checked", sp.high_grade_depth_m)),
        a("A-STRATA-DOMINANT-PARCEL", "a cell's published column is its dominant parcel's stack (largest kernel deposit); oceanic cells: pelagic shale (the sediment field) over oceanic igneous crust", sp.publish_depth_m, "m", "listed to this depth; deeper layers are cut and flagged"),
    ];
    // SEED-STRATIGRAPHY stage A.1 (lane SA): parcel-stacks-v2 only (v1
    // records unchanged).
    if let Some(v) = &sp.v2 {
        out.push(a("A-STRATA-BELT-PLACED-SHORTENING", "parcel-stacks-v2: plutons and pure-shear thickening come only from the belt's placed collision supply (and arcs): the stack takes the cell's net shortening-lineage change clamped to [0, placed]; the rest (plateau cap/spill, channel-flow inflow beyond the placed supply, outflow, delamination) changes the basement only", sp.shortening_pluton_fraction, "1", "lane SA 2026-10-13 PREREG C1 with revision 1 (net clamped to the placed supply); removes granite sills under continent interiors"));
        out.push(a("A-STRATA-SUBCELL-SHELF", "parcel-stacks-v2: a marine cell deeper than the shelf depth takes the shelf facies when its shelf fraction min(f_dist, f_depth) reaches the majority fraction; f_dist = min(1, shelf width / sqrt(cell area)) on coast cells (an 8-neighbour not ocean), 0 elsewhere; f_depth = (slope foot - depth) / (slope foot - shelf depth / 2) clamped to [0, 1]", v.shelf_width_m, "m", &format!("slope foot {} m, majority {}; lane SA PREREG C2", v.slope_foot_depth_m, v.shelf_majority_fraction)));
        out.push(a("A-STRATA-PLATFORM-COVER", "parcel-stacks-v2: seeded parcels carry a t = 0 platform cover (bottom-up sandstone, shale, carbonate) of this thickness at or below the full-cover elevation, falling linearly to none at the zero-cover elevation (exposed shield), capped at half the column; a relabel of the basement top (volume neutral); one craton age in the seeding, so position enters only through the t = 0 elevation", v.platform_cover_max_m, "m", &format!("full at {} m, none at {} m; sandstone {}, carbonate {}; ages {} to {} My; lane SA PREREG C3", v.platform_full_elev_m, v.platform_zero_elev_m, v.platform_sandstone_fraction, v.platform_carbonate_fraction, v.platform_age_base_my, v.platform_age_top_my)));
        out.push(a("A-STRATA-ROCK-UPLIFT", "parcel-stacks-v2: rock_uplift_m is each parcel's cumulative change of post-tectonic elevation against the previous step's post-erosion elevation (tectonics plus isostatic rebound, not erosion or deposition) since t = 0 or birth, kernel-averaged per cell; uplift_grad_east/north its gradient; stage B reads strata dipping toward minus the gradient", 0.0, "m", "lane SA PREREG C5; a regional tilt input, not a fold or fault model"));
    }
    out
}

/// SEED-STRATIGRAPHY stage A.1 (lane SA): the `parcel-stacks-v2` rock
/// uplift field and its gradient into a published record.
fn strata_uplift_into(
    d: &mut pf_contracts::terrain::StratigraphyDomainV0,
    grid: &GridSpec,
    layer: &ParcelLayer,
    state: &CrustState,
) -> Result<(), PreviewError> {
    let (u, gx, gy) = pf_geodynamics::conserve::rock_uplift_field(grid, layer, state)?;
    let f = |v: Vec<f64>| v.into_iter().map(|x| x as f32).collect::<Vec<f32>>();
    d.rock_uplift_m = Some(f(u));
    d.uplift_grad_east = Some(f(gx));
    d.uplift_grad_north = Some(f(gy));
    Ok(())
}

/// SEED-STRATIGRAPHY (lane ST): the step's tagged per-cell deltas.
#[allow(clippy::too_many_arguments)]
fn strata_tags(
    grid: &GridSpec,
    sp: &StrataParams,
    ledger: &AdvanceLedger,
    backarc_dh: Option<&[f64]>,
    step: &pf_surface::erosion::ErosionStep,
    elev: &[f64],
    sea: f64,
    ocean: &[bool],
    routing: &Routing,
    dt: f64,
) -> Result<CellTags, PreviewError> {
    let n = grid.cells();
    let mut t = CellTags::zeros(n);
    for c in 0..n {
        let out = step.sediment_removed_m[c] + step.bedrock_removed_m[c];
        let inn = step.sediment_added_m[c];
        t.erosion_m[c] = (out - inn).max(0.0);
        t.deposit_m[c] = (inn - out).max(0.0);
    }
    for &(c, dh, _) in ledger.hotspot_cells.iter().chain(ledger.lip_cells.iter()) {
        t.plume_m[c as usize] += dh.max(0.0);
    }
    for &(c, dh) in &ledger.rift_thinning {
        t.thinning_m[c as usize] += dh.max(0.0);
    }
    if let Some(b) = backarc_dh {
        for c in 0..n {
            if b[c] < 0.0 {
                t.thinning_m[c] += -b[c];
            } else {
                t.basement_add_m[c] += b[c];
            }
        }
    }
    let endorheic: Vec<bool> = (0..n)
        .map(|c| routing.drainage_class(pf_contracts::spatial::ElementId(c as u64)) == Some(pf_contracts::domain::DrainageClassV0::Endorheic))
        .collect();
    let flow_to: Vec<Option<usize>> = routing.flow_to.iter().map(|x| x.map(|e| e.0 as usize)).collect();
    let rate: Vec<f64> = t.deposit_m.iter().map(|d| d / dt).collect();
    let (class, land, warm_shelf) = pf_geodynamics::stratigraphy::classify_cells_with_shelf(
        grid,
        &pf_geodynamics::stratigraphy::EnvInputs {
            elevation_m: elev,
            sea_m: sea,
            ocean,
            fill_level_m: &routing.fill_level_m,
            endorheic: &endorheic,
            flow_to: &flow_to,
            contributing_area_m2: &routing.contributing_area_m2,
            deposit_rate_m_per_my: &rate,
        },
        sp,
    )?;
    t.class = class;
    t.land = land;
    // SEED-STRATIGRAPHY stage A.1 (lane SA): parcel-stacks-v2 inputs.
    if sp.v2.is_some() {
        t.warm_shelf = warm_shelf;
        t.belt_placed_m = match ledger.belt_budget.as_ref() {
            Some(bb) => bb.placed_m.clone(),
            None => vec![0.0; n],
        };
        t.tect_elev_m = elev.to_vec();
        t.post_elev_m = (0..n).map(|c| elev[c] + step.dz_m[c]).collect();
    }
    Ok(t)
}

/// SEED-CARRY-CONSERVE (lane CV): the per-step continental ledger of
/// `continental-parcel-carry-v1` (areas m^2, volumes m^3) and the sutures
/// of `contact-length-pole-share-v1`.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ConserveStepLog {
    /// Parcels after the step and their area, crust volume and sediment.
    pub parcels: u64,
    pub parcel_area_m2: f64,
    pub parcel_volume_m3: f64,
    pub parcel_sediment_m3: f64,
    /// The projection: continental cells, their area and crust volume.
    pub grid_cont_cells: u64,
    pub grid_cont_area_m2: f64,
    pub grid_cont_volume_m3: f64,
    /// Area sink: collision consumption (crust and sediment to the belt).
    pub consumed_parcels: u64,
    pub consumed_area_m2: f64,
    pub consumed_m3: f64,
    pub consumed_sediment_m3: f64,
    pub overlap_cells: u64,
    /// Area source: rift births that ended continental (new parcels).
    pub created_parcels: u64,
    pub created_area_m2: f64,
    pub created_m3: f64,
    pub created_sediment_m3: f64,
    /// Changes written back to parcels (belt, arc, rift, plume, erosion,
    /// deposition), clamps, and the parcel volume / sediment residual
    /// after - (before - consumed + written back + created + clamped).
    pub writeback_m3: f64,
    pub writeback_sediment_m3: f64,
    pub clamped_m3: f64,
    pub clamped_sediment_m3: f64,
    pub volume_residual_m3: f64,
    pub sediment_residual_m3: f64,
    /// Projection overrides and fills.
    pub overridden_oceanic_cells: u64,
    pub overridden_births: u64,
    pub vacated_copied: u64,
    pub vacated_births: u64,
    pub lost_cells: u64,
    /// Plate-id fragments folded (pieces, cells, parcels relabelled).
    pub fragments_folded: u64,
    pub fragment_cells: u64,
    pub fragment_parcels: u64,
    /// Sutures made this step and the largest converging contact (m).
    pub sutures: u64,
    pub max_contact_m: f64,
}

/// SEED-CARRY-CONSERVE (lane CV): one suture of
/// `contact-length-pole-share-v1`.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SutureLog {
    pub t_my: f64,
    pub a: u32,
    pub b: u32,
    pub contact_m: f64,
    pub winner_root: u32,
    pub loser_root: u32,
}

/// SEED-CARRY-CONSERVE (lane CV): run state of the parcel carry (the
/// parcels) and the suture law (union-find parents, events).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ConserveRunState {
    #[serde(default,skip_serializing_if="Option::is_none")]
    pub fluvial_trace:Option<pf_surface::erosion::fluvial_trace::FluvialArchive>,
    #[serde(with = "parcel_layer_serde")]
    pub parcels: ParcelLayer,
    pub suture_parent: Vec<(u32, u32)>,
    pub sutures: Vec<SutureLog>,
    /// SEED-PLATE-REORG (lane PR): `plate-reorg-scripted-v1` state; skipped
    /// when absent (every other law's checkpoint is byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reorg: Option<ReorgRunState>,
    /// SEED-RANGES (lane RG): `contact-length-indentation-taper-v1`
    /// indentations; skipped when absent (every other law's checkpoint is
    /// byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub indent: Option<Vec<IndentLog>>,
    /// SEED-MAGMATISM (lane MG): `slab-arc-backarc-llsvp-v1` vent markers,
    /// arc heat and step logs; skipped when absent (every other law's
    /// checkpoint is byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub magma: Option<crate::magmatism::MagmaRunState>,
}

/// SEED-RANGES (lane RG): one indentation of
/// `contact-length-indentation-taper-v1` (run state and log).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct IndentLog {
    pub t0_my: f64,
    pub winner: u32,
    pub loser: u32,
    pub members: Vec<u32>,
    pub omega_rel: [f64; 3],
    pub shortening_m: f64,
    pub v0_m_per_my: Option<f64>,
    pub end_t_my: Option<f64>,
    pub end_cause: u8,
}

impl IndentLog {
    fn to_state(&self) -> Indentation {
        Indentation {
            t0_my: self.t0_my,
            winner: pf_contracts::domain::PlateId(self.winner),
            loser: pf_contracts::domain::PlateId(self.loser),
            members: self.members.iter().map(|&m| pf_contracts::domain::PlateId(m)).collect(),
            omega_rel: self.omega_rel,
            shortening_m: self.shortening_m,
            v0_m_per_my: self.v0_m_per_my,
            end_t_my: self.end_t_my,
            end_cause: self.end_cause,
        }
    }
    fn from_state(s: &Indentation) -> Self {
        Self {
            t0_my: s.t0_my,
            winner: s.winner.0,
            loser: s.loser.0,
            members: s.members.iter().map(|m| m.0).collect(),
            omega_rel: s.omega_rel,
            shortening_m: s.shortening_m,
            v0_m_per_my: s.v0_m_per_my,
            end_t_my: s.end_t_my,
            end_cause: s.end_cause,
        }
    }
}

/// SEED-RANGES (lane RG): per-step log of the indentation and plateau laws
/// (volumes m^3).
#[derive(Debug, Clone, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub struct RangesStepLog {
    pub indent_started: u32,
    pub indent_active: u32,
    pub indent_ended: u32,
    pub indent_contact_m: f64,
    pub indent_mean_rate_m_per_my: f64,
    pub indent_max_taper: f64,
    pub plateau_source_cells: u64,
    pub plateau_excess_m3: f64,
    pub plateau_moved_m3: f64,
    pub plateau_filled_cells: u64,
    pub plateau_delaminated_m3: f64,
    pub plateau_oceanic_delaminated_m3: f64,
    pub plateau_max_distance_m: f64,
    pub plateau_cells_at_cap: u64,
}

impl RangesStepLog {
    fn of(i: Option<&IndentStepReport>, p: Option<&PlateauReport>) -> Self {
        let mut l = RangesStepLog::default();
        if let Some(i) = i {
            l.indent_started = i.started;
            l.indent_active = i.active;
            l.indent_ended = i.ended;
            l.indent_contact_m = i.contact_m;
            l.indent_mean_rate_m_per_my = i.mean_rate_m_per_my;
            l.indent_max_taper = i.max_taper;
        }
        if let Some(p) = p {
            l.plateau_source_cells = p.source_cells;
            l.plateau_excess_m3 = p.excess_m3;
            l.plateau_moved_m3 = p.moved_m3;
            l.plateau_filled_cells = p.filled_cells;
            l.plateau_delaminated_m3 = p.delaminated_m3;
            l.plateau_oceanic_delaminated_m3 = p.oceanic_delaminated_m3;
            l.plateau_max_distance_m = p.max_distance_m;
            l.plateau_cells_at_cap = p.cells_at_cap;
        }
        l
    }
}

/// SEED-PLATE-REORG (lane PR): one relation of `plate-reorg-scripted-v1`
/// (child root = parent root + omega_rel across re-solves while young).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ReorgRelationLog {
    pub parent: u32,
    pub child: u32,
    pub omega_rel: [f64; 3],
    pub t_my: f64,
    pub kind: u8,
}

/// SEED-PLATE-REORG (lane PR): one event of `plate-reorg-scripted-v1`.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ReorgEventLog {
    pub t_my: f64,
    pub kind: String,
    pub root_a: u32,
    pub root_b: u32,
    pub cells_moved: u64,
    pub group_share: f64,
    pub side_b_share: f64,
    pub cut_length_m: f64,
    pub weak_length_m: f64,
    pub groups_counted: u32,
}

/// SEED-PLATE-REORG (lane PR): run state of `plate-reorg-scripted-v1`.
#[derive(Debug, Clone, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub struct ReorgRunState {
    pub tenure: Vec<(u32, f64)>,
    pub cascade: Vec<(u32, f64)>,
    pub relations: Vec<ReorgRelationLog>,
    pub next_event: u32,
    pub breakups: u32,
    pub initiations: u32,
    pub events: Vec<ReorgEventLog>,
}

impl ReorgRunState {
    fn to_state(&self) -> ReorgState {
        ReorgState {
            tenure: self.tenure.clone(),
            cascade: self.cascade.clone(),
            relations: self
                .relations
                .iter()
                .map(|r| pf_geodynamics::reorg::Relation {
                    parent: r.parent,
                    child: r.child,
                    omega_rel: r.omega_rel,
                    t_my: r.t_my,
                    kind: r.kind,
                })
                .collect(),
            next_event: self.next_event,
            breakups: self.breakups,
            initiations: self.initiations,
            events: self
                .events
                .iter()
                .map(|e| pf_geodynamics::reorg::ReorgEvent {
                    t_my: e.t_my,
                    kind: e.kind.clone(),
                    root_a: e.root_a,
                    root_b: e.root_b,
                    cells_moved: e.cells_moved,
                    group_share: e.group_share,
                    side_b_share: e.side_b_share,
                    cut_length_m: e.cut_length_m,
                    weak_length_m: e.weak_length_m,
                    groups_counted: e.groups_counted,
                })
                .collect(),
        }
    }
    fn from_state(s: &ReorgState) -> Self {
        Self {
            tenure: s.tenure.clone(),
            cascade: s.cascade.clone(),
            relations: s
                .relations
                .iter()
                .map(|r| ReorgRelationLog {
                    parent: r.parent,
                    child: r.child,
                    omega_rel: r.omega_rel,
                    t_my: r.t_my,
                    kind: r.kind,
                })
                .collect(),
            next_event: s.next_event,
            breakups: s.breakups,
            initiations: s.initiations,
            events: s
                .events
                .iter()
                .map(|e| ReorgEventLog {
                    t_my: e.t_my,
                    kind: e.kind.clone(),
                    root_a: e.root_a,
                    root_b: e.root_b,
                    cells_moved: e.cells_moved,
                    group_share: e.group_share,
                    side_b_share: e.side_b_share,
                    cut_length_m: e.cut_length_m,
                    weak_length_m: e.weak_length_m,
                    groups_counted: e.groups_counted,
                })
                .collect(),
        }
    }
}

impl ConserveRunState {
    fn suture_state(&self) -> SutureState {
        SutureState {
            parent: self
                .suture_parent
                .iter()
                .map(|&(a, b)| (pf_contracts::domain::PlateId(a), pf_contracts::domain::PlateId(b)))
                .collect(),
        }
    }
}

mod parcel_layer_serde {
    use pf_contracts::domain::PlateId;
    use pf_geodynamics::conserve::ParcelLayer;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    #[derive(Serialize, Deserialize)]
    struct Ser {
        pos: Vec<[f64; 3]>,
        plate: Vec<u32>,
        area_m2: Vec<f64>,
        thickness_m: Vec<f64>,
        sediment_m: Vec<f64>,
        eroded_total_m: Vec<f64>,
        arc_added_m: Vec<f64>,
        shortening_added_m: Vec<f64>,
        birth_t_my: Vec<f64>,
        /// SEED-STRATIGRAPHY (lane ST): skipped when absent (every other
        /// law's checkpoint is byte-identical).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        strata: Option<Vec<StackSer>>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        identity: Option<serde_json::Value>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        history: Option<serde_json::Value>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        layer_provenance: Option<Vec<pf_contracts::layer_provenance::StackProvenance>>,
    }

    /// One stack: basement [m, lith, age, last exposure] and the cover,
    /// bottom-up, as (thickness, age base, age top, max burial, lith, env,
    /// grade, flags).
    #[derive(Serialize, Deserialize)]
    struct StackSer {
        b: (f64, u8, f64, f64),
        cover: Vec<(f64, f64, f64, f64, u8, u8, u8, u8)>,
        /// SEED-STRATIGRAPHY stage A.1 (lane SA, parcel-stacks-v2): (rock
        /// uplift, reference elevation); skipped when both are +0.0, which
        /// is always so under v1 (v1 checkpoints byte-identical).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        u: Option<(f64, f64)>,
    }

    fn to_ser(st: &pf_geodynamics::stratigraphy::StrataStack) -> StackSer {
        StackSer {
            b: (st.basement_m, st.basement_lith, st.basement_age_my, st.last_exposure_my),
            cover: st
                .cover
                .iter()
                .map(|l| (l.thickness_m, l.age_base_my, l.age_top_my, l.max_burial_m, l.lith, l.env, l.grade, l.flags))
                .collect(),
            u: if st.rock_uplift_m.to_bits() == 0 && st.elev_ref_m.to_bits() == 0 {
                None
            } else {
                Some((st.rock_uplift_m, st.elev_ref_m))
            },
        }
    }

    fn from_ser(x: StackSer) -> pf_geodynamics::stratigraphy::StrataStack {
        pf_geodynamics::stratigraphy::StrataStack {
            basement_m: x.b.0,
            basement_lith: x.b.1,
            basement_age_my: x.b.2,
            last_exposure_my: x.b.3,
            cover: x
                .cover
                .into_iter()
                .map(|c| pf_geodynamics::stratigraphy::StrataLayer {
                    thickness_m: c.0,
                    age_base_my: c.1,
                    age_top_my: c.2,
                    max_burial_m: c.3,
                    lith: c.4,
                    env: c.5,
                    grade: c.6,
                    flags: c.7,
                })
                .collect(),
            rock_uplift_m: x.u.map_or(0.0, |u| u.0),
            elev_ref_m: x.u.map_or(0.0, |u| u.1),
        }
    }

    pub fn serialize<S: Serializer>(l: &ParcelLayer, s: S) -> Result<S::Ok, S::Error> {
        Ser {
            pos: l.pos.clone(),
            plate: l.plate.iter().map(|p| p.0).collect(),
            area_m2: l.area_m2.clone(),
            thickness_m: l.thickness_m.clone(),
            sediment_m: l.sediment_m.clone(),
            eroded_total_m: l.eroded_total_m.clone(),
            arc_added_m: l.arc_added_m.clone(),
            shortening_added_m: l.shortening_added_m.clone(),
            birth_t_my: l.birth_t_my.clone(),
            strata: l.strata.as_ref().map(|v| v.iter().map(to_ser).collect()),
            identity: l.identity.as_ref().map(|v| v.to_json()),
            history: l.history.as_ref().map(|v| v.to_json()),
            layer_provenance: l.layer_provenance.clone(),
        }
        .serialize(s)
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<ParcelLayer, D::Error> {
        let x = Ser::deserialize(d)?;
        let identity = x.identity.map(pf_geodynamics::conserve::CarrierIdentity::from_json)
            .transpose().map_err(serde::de::Error::custom)?;
        let history = x.history.map(pf_geodynamics::carrier_history::StackHistory::from_json)
            .transpose().map_err(serde::de::Error::custom)?;
        let k = x.pos.len();
        for (name, len) in [
            ("plate", x.plate.len()),
            ("area_m2", x.area_m2.len()),
            ("thickness_m", x.thickness_m.len()),
            ("sediment_m", x.sediment_m.len()),
            ("eroded_total_m", x.eroded_total_m.len()),
            ("arc_added_m", x.arc_added_m.len()),
            ("shortening_added_m", x.shortening_added_m.len()),
            ("birth_t_my", x.birth_t_my.len()),
            ("strata", x.strata.as_ref().map_or(k, |v| v.len())),
            ("layer_provenance", x.layer_provenance.as_ref().map_or(k, |v|v.len())),
            ("identity", identity.as_ref().map_or(k, |v| v.ids().len())),
        ] {
            if len != k {
                return Err(serde::de::Error::custom(format!(
                    "parcel field {name}: {len} entries for {k} parcels"
                )));
            }
        }
        let result = ParcelLayer {
            pos: x.pos,
            plate: x.plate.into_iter().map(PlateId).collect(),
            area_m2: x.area_m2,
            thickness_m: x.thickness_m,
            sediment_m: x.sediment_m,
            eroded_total_m: x.eroded_total_m,
            arc_added_m: x.arc_added_m,
            shortening_added_m: x.shortening_added_m,
            birth_t_my: x.birth_t_my,
            strata: x.strata.map(|v| v.into_iter().map(from_ser).collect()),
            identity,
            history,
            layer_provenance: x.layer_provenance,
        };
        if let Some(p)=&result.layer_provenance {
            let stacks=result.strata.as_ref().ok_or_else(||serde::de::Error::custom("provenance without strata"))?;
            if result.identity.is_none(){return Err(serde::de::Error::custom("provenance without IDs"));}
            for (p,st) in p.iter().zip(stacks){pf_geodynamics::layer_provenance::validate(p,st).map_err(serde::de::Error::custom)?;}
        }
        if let Some(h)=&result.history {
            h.validate_context(result.identity.as_ref().ok_or_else(||serde::de::Error::custom("history without carrier IDs"))?,
                h.latest_time(),result.strata.as_deref()).map_err(serde::de::Error::custom)?;
        }
        Ok(result)
    }
}

/// SEED-CARRY-CONSERVE (lane CV): the grid's continental ledger after a
/// step, computed under EVERY law (feeds no output): continental crust area
/// and volume, sediment on continental vs oceanic cells, land area.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ContinentalLedgerV1 {
    pub cont_area_m2: f64,
    pub cont_volume_m3: f64,
    pub sediment_cont_m3: f64,
    pub sediment_ocean_m3: f64,
    pub land_area_m2: f64,
}

/// SEED-TRENCH (lane TN): per-step `slab-flexure-trench-v1` diagnostics
/// (the term handed to the isostasy, water-loaded m).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TrenchStepLog {
    /// Trench (downgoing, overridden oceanic) cells and flexed cells.
    pub trench_cells: u64,
    pub flexed_cells: u64,
    /// Largest trench-axis deflection w0 and outer-rise term (m).
    pub max_w0_m: f64,
    pub max_bulge_wl_m: f64,
}

/// SEED-EXTRAS (lane EX): per-step plume, rift and split bookkeeping.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ExtrasStepLog {
    /// Magmatic crust added by plume heads and tails (m^3), of which heads.
    pub hotspot_created_m3: f64,
    pub hotspot_head_m3: f64,
    pub plume_heads: u64,
    /// Rift births (stretched continental crust), rift candidates left
    /// oceanic at breakup, candidates on converging plates (skipped), and
    /// the volume moved from flanks into rift births (m^3).
    pub rift_births: u64,
    pub rift_breakup_births: u64,
    pub rift_convergent_skipped: u64,
    pub rift_moved_m3: f64,
    /// Plate splits this step and the plate count after the step.
    pub plate_splits: u64,
    pub plates: u64,
    /// SEED-RIFT-SHOULDERS (lane RS): `shoulder-uplift-thermal-v1` only
    /// (skipped when absent, so extras_v1 logs serialize byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rift_flank: Option<RiftFlankStepLog>,
    /// SEED-FLOODBASALT (lane FB): `plume-head-lip-swell-v1` only (skipped
    /// when absent, so every other log serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flood_basalt: Option<FloodBasaltStepLog>,
}

/// SEED-FLOODBASALT (lane FB): per-step `plume-head-lip-swell-v1`
/// bookkeeping.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct FloodBasaltStepLog {
    /// LIP crust erupted this step (m^3; inside `hotspot_created_m3`), the
    /// part on continental columns, and the cells it reached.
    pub lip_created_m3: f64,
    pub lip_continental_m3: f64,
    pub lip_cells: u64,
    /// Largest swell term and the columns with more than 100 m of it
    /// (water-loaded m).
    pub max_swell_wl_m: f64,
    pub swell_cells: u64,
}

/// SEED-RIFT-SHOULDERS (lane RS): per-step `shoulder-uplift-thermal-v1`
/// diagnostics (the term handed to the isostasy, water-loaded m).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RiftFlankStepLog {
    /// Continental columns with a stretching factor above 1, and columns
    /// carrying a shoulder.
    pub stretched_cells: u64,
    pub shoulder_cells: u64,
    pub max_thermal_wl_m: f64,
    pub max_shoulder_wl_m: f64,
}

/// SEED-RIFT-SHOULDERS (lane RS): stretching markers carried with the
/// crust column (1 / -1 at births and on unstretched columns).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RiftFlankMarkers {
    /// Stretching factor beta (>= 1).
    pub beta: Vec<f64>,
    /// Equivalent start time of the column's thermal anomaly (My).
    pub thermal_t0_my: Vec<f64>,
}

/// SEED-EXTRAS (lane EX): one plume-head plate split.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PlateSplitRecord {
    pub t_my: f64,
    pub plume: u64,
    pub parent: u32,
    pub child: u32,
    pub cut_normal: [f64; 3],
    pub cells_moved: u64,
    pub child_pole_lat_deg: f64,
    pub child_pole_lon_deg: f64,
    pub child_pole_rate_deg_per_my: f64,
}

/// SEED-EXTRAS (lane EX): per-cell markers carried with the crust (by the
/// advance's source map; 0 / -1 at births) plus the split history.
/// `-1.0` means "never" in the time fields.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ExtrasRunState {
    /// Cumulative plume crust added to the column (m).
    pub volcanic_added_m: Vec<f64>,
    /// Latest plume eruption time onto the column (My).
    pub volcanic_t_my: Vec<f64>,
    /// Time the column was born as rift crust (My).
    pub rift_t_my: Vec<f64>,
    /// Cumulative rift stretching loss of the column (m).
    pub rift_thinned_m: Vec<f64>,
    pub splits: Vec<PlateSplitRecord>,
    /// SEED-RIFT-SHOULDERS (lane RS): `shoulder-uplift-thermal-v1` only.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rift_flank: Option<RiftFlankMarkers>,
    /// SEED-FLOODBASALT (lane FB): cumulative LIP crust added to the column
    /// (m; also inside `volcanic_added_m`); `plume-head-lip-swell-v1` only.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lip_added_m: Option<Vec<f64>>,
}

/// SEED-CARRY-CC (lane CC): per-step `front-rescue-carry-v1` bookkeeping.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CarryStepLog {
    /// Rescue claims admitted.
    pub rescues: u64,
    /// Continent-continent overlap losers, those still carried elsewhere,
    /// and the consumed parcels (count, source-cell area m^2, volume m^3,
    /// sediment m^3); the consumed volume is `belt_in_m3`.
    pub cc_losers: u64,
    pub cc_losers_carried: u64,
    pub consumed_parcels: u64,
    pub consumed_area_m2: f64,
    pub consumed_m3: f64,
    pub consumed_sediment_m3: f64,
    /// Continental parcels still carried nowhere (silent drops; remap
    /// residual) and their volume (m^3).
    pub continental_dropped_parcels: u64,
    pub continental_dropped_m3: f64,
    /// Ridge births between continents of two plates.
    pub sliver_births: u64,
}

/// SEED-OROGEN state carried between steps (checkpointed by the genesis
/// module only when a new law is selected).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct OrogenRunState {
    /// Subgrid relief r = cell mean - trunk level (m), carried with the
    /// crust; `None` unless `stream-power-subgrid-v1`.
    pub subgrid_relief_m: Option<Vec<f64>>,
    pub steps: Vec<OrogenStepLog>,
    /// SEED-EXTRAS (lane EX): markers and splits; `None` unless a plume or
    /// rift law is selected.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub extras: Option<ExtrasRunState>,
    /// SEED-CARRY-CONSERVE (lane CV): the parcels and sutures; `None`
    /// unless `continental-parcel-carry-v1` or `contact-length-pole-share-v1`
    /// is selected.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conserve: Option<ConserveRunState>,
}

/// The sea-level solve of one two-column isostasy evaluation.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct SeaLevelSolve {
    pub sea_level_gauge_m: f64,
    pub flooded_area_fraction: f64,
    pub ocean_volume_residual_m3: f64,
    /// Compensation depth used (m): derived (deepest root + margin) or the
    /// fixed `compensation_depth_m`.
    pub compensation_depth_m: f64,
    /// Deepest crustal root over all columns (m below the gauge zero) and
    /// that column's crust and sediment thickness (m).
    pub deepest_root_m: f64,
    pub deepest_root_crust_sediment_m: [f64; 2],
}

/// GDH1 depth below the ridge-crest datum as a function of crust age;
/// returns the thermal ELEVATION relative to old lithosphere (>= 0).
pub fn gdh1_thermal_elevation_m(age_my: f64) -> f64 {
    let age = age_my.max(0.0);
    let depth = if age < GDH1_SWITCH_MY {
        GDH1_RIDGE_DEPTH_M + GDH1_SQRT_COEF_M_PER_SQRT_MY * age.sqrt()
    } else {
        GDH1_ASYMPTOTE_M - GDH1_EXP_COEF_M * (-GDH1_EXP_RATE_PER_MY * age).exp()
    };
    (GDH1_ASYMPTOTE_M - depth).max(0.0)
}

/// WS-CONST: the law set and universe constants of a profile (typed
/// refusal mapped to `PreviewError::Params`).
pub fn universe_of(p: &PlanetParameterProfile) -> Result<UniverseV1, PreviewError> {
    p.universe()
        .map_err(|e| PreviewError::Params(format!("universe: {e}")))
}

/// Elevation of every cell from the crust state: Airy on (crust +
/// sediment) by the legacy single-column law, water-loading
/// amplification below sea level, GDH1 thermal term for oceanic crust.
/// WS-CONST: densities, the reference column and the sea level are read
/// through the profile's `UniverseV1` (same values as before).
pub fn elevation_from_state(
    state: &CrustState,
    params: &PlanetParameterProfile,
) -> Result<Vec<f64>, PreviewError> {
    elevation_single_column(state, &universe_of(params)?)
}

fn elevation_single_column(
    state: &CrustState,
    u: &UniverseV1,
) -> Result<Vec<f64>, PreviewError> {
    let k = SingleColumnConstants::from_universe(u)
        .map_err(|e| PreviewError::Domain(format!("isostasy: {e}")))?;
    if k.compensation != Compensation::Airy {
        // As the dispatch `isostatic_elevation_m(.., None)` refused before.
        return Err(PreviewError::Domain(
            "isostasy: Flexure mode requires the flexure solver (never silently Airy)".into(),
        ));
    }
    let sea = match u
        .sea_level_input()
        .map_err(|e| PreviewError::Params(format!("universe: {e}")))?
    {
        SeaLevelInputV1::Prescribed { sea_level_m } => sea_level_m,
        SeaLevelInputV1::OceanVolume { .. } => {
            return Err(PreviewError::Params(
                "universe: the single-column law needs a prescribed sea level".into(),
            ))
        }
    };
    let rho_m = k.mantle_density_kg_m3;
    let water_factor = rho_m / (rho_m - u.constants.seawater_density_kg_m3.value);
    let n = state.cells();
    let mut out = Vec::with_capacity(n);
    for c in 0..n {
        let column = state.thickness_m[c] + state.sediment_m[c];
        let airy = single_column_airy_elevation_m(column, &k)
            .map_err(|e| PreviewError::Domain(format!("isostasy at cell {c}: {e}")))?;
        let loaded = if airy < sea {
            sea + (airy - sea) * water_factor
        } else {
            airy
        };
        let thermal = if state.kind[c] == CrustKindCodeV0::Oceanic {
            gdh1_thermal_elevation_m(state.age_my(c))
        } else {
            0.0
        };
        out.push(loaded + thermal);
    }
    Ok(out)
}

/// SEED-WORLD-COMPOSE (lane CP): per cell, the shortening-budget distance
/// bin (0-3 by `SHORTENING_BUDGET_BINS_KM`, 4 = no source cell) from a
/// multi-source Dijkstra over edge neighbours (centre great-circle steps).
fn shortening_distance_bins(
    grid: &GridSpec,
    state: &CrustState,
    suture_parent: &BTreeMap<u32, u32>,
) -> Result<Vec<u8>, PreviewError> {
    let n = grid.cells();
    let nlon = grid.nlon();
    let cents = cell_centroids(grid).map_err(PreviewError::Domain)?;
    let r = grid.radius_m();
    let find = |mut p: u32| {
        let mut hops = 0;
        while let Some(&q) = suture_parent.get(&p) {
            if q == p || hops > 4096 {
                break;
            }
            p = q;
            hops += 1;
        }
        p
    };
    let sub = pf_contracts::terrain::BoundaryClassCodeV0::Subduction.code();
    let col = pf_contracts::terrain::BoundaryClassCodeV0::Collision.code();
    let nb4 = |c: usize| -> [usize; 4] {
        let (i, j) = (c / nlon, c % nlon);
        let up = if i + 1 < grid.nlat() { (i + 1) * nlon + j } else { c };
        let dn = if i > 0 { (i - 1) * nlon + j } else { c };
        [i * nlon + (j + 1) % nlon, i * nlon + (j + nlon - 1) % nlon, up, dn]
    };
    let mut dist = vec![f64::INFINITY; n];
    let mut heap = std::collections::BinaryHeap::new();
    for c in 0..n {
        let mut src = state.boundary_class[c] == sub || state.boundary_class[c] == col;
        if !src && state.kind[c] == CrustKindCodeV0::Continental {
            for k in nb4(c) {
                if k != c
                    && state.kind[k] == CrustKindCodeV0::Continental
                    && state.plate[k] != state.plate[c]
                    && find(state.plate[k].0) == find(state.plate[c].0)
                {
                    src = true;
                    break;
                }
            }
        }
        if src {
            dist[c] = 0.0;
            heap.push(std::cmp::Reverse((0u64, c)));
        }
    }
    // SEED-PERF (lane SP, 2026-10-13): the edge lengths (the same
    // expression per (c, k)) are evaluated up front in parallel; the
    // Dijkstra itself is unchanged and serial.
    let elen: Vec<[f64; 4]> = pf_geodynamics::kinematic::par::map_indexed(n, |c| {
        let mut e = [0.0f64; 4];
        for (q, k) in nb4(c).into_iter().enumerate() {
            if k != c {
                let dot = cents[c][0] * cents[k][0] + cents[c][1] * cents[k][1] + cents[c][2] * cents[k][2];
                e[q] = dot.clamp(-1.0, 1.0).acos() * r;
            }
        }
        e
    });
    while let Some(std::cmp::Reverse((dk, c))) = heap.pop() {
        let d = f64::from_bits(dk);
        if d > dist[c] {
            continue;
        }
        for (q, k) in nb4(c).into_iter().enumerate() {
            if k == c {
                continue;
            }
            let nd = d + elen[c][q];
            if nd < dist[k] {
                dist[k] = nd;
                heap.push(std::cmp::Reverse((nd.to_bits(), k)));
            }
        }
    }
    Ok(dist
        .iter()
        .map(|&d| {
            if !d.is_finite() {
                4
            } else {
                let km = d / 1000.0;
                SHORTENING_BUDGET_BINS_KM.iter().position(|&b| km < b).unwrap_or(3) as u8
            }
        })
        .collect())
}

/// SEED-WORLD-COMPOSE (lane CP) `slab-arc-backarc-llsvp-v2`: add a
/// water-loaded dynamic term after the sea-level solve. Submerged columns
/// rise by the water-loaded amount until they reach sea level; above it the
/// rest (and every land column) rises by the dry amount (`dry_factor` x).
/// The solved sea level is not changed: the term displaces no ocean water.
fn add_dynamic_term(elev: &mut [f64], term_wl: &[f64], sea: f64, dry_factor: f64) {
    for (e, &t) in elev.iter_mut().zip(term_wl) {
        if !(t > 0.0) {
            continue;
        }
        if *e >= sea {
            *e += t * dry_factor;
        } else {
            let need = sea - *e;
            if t <= need {
                *e += t;
            } else {
                *e = sea + (t - need) * dry_factor;
            }
        }
    }
}

/// SEED-FLOODBASALT (lane FB): a column term (or zeros) plus `extra`.
fn add_column_term(base: Option<&[f64]>, extra: &[f64]) -> Vec<f64> {
    match base {
        Some(b) => b.iter().zip(extra).map(|(a, e)| a + e).collect(),
        None => extra.to_vec(),
    }
}

/// Elevation on the grid by the SELECTED isostasy law of the profile's
/// `UniverseV1` (WS-CONST): the legacy single-column Airy
/// (`elevation_from_state`, byte-identical) under
/// `single-column-reference-v0`; otherwise the opt-in two-column model with
/// a volume-primary sea level (lane F FREEBOARD-SPIKE, wired by lane K2).
/// Output elevations are relative to the solved sea level placed at the
/// sea-level datum, so routing/erosion keep using `sea` unchanged.
pub fn elevation_on_grid(
    state: &CrustState,
    params: &PlanetParameterProfile,
    areas: &[f64],
) -> Result<(Vec<f64>, Option<SeaLevelSolve>), PreviewError> {
    elevation_on_grid_with(state, params, areas, None)
}

/// `elevation_on_grid` plus an optional per-cell column term for
/// CONTINENTAL columns (water-loaded m, like the oceanic GDH1 term): the
/// SEED-RIFT-SHOULDERS `shoulder-uplift-thermal-v1` rift thermal and
/// shoulder uplift. `None` is `elevation_on_grid` exactly. Refused under
/// the single-column law (the law set refuses that combination too).
pub fn elevation_on_grid_with(
    state: &CrustState,
    params: &PlanetParameterProfile,
    areas: &[f64],
    continental_term_wl_m: Option<&[f64]>,
) -> Result<(Vec<f64>, Option<SeaLevelSolve>), PreviewError> {
    elevation_on_grid_terms(state, params, areas, continental_term_wl_m, None)
}

/// `elevation_on_grid_with` plus an optional per-cell column term for
/// OCEANIC columns added to their GDH1 term (water-loaded m): the
/// SEED-TRENCH `slab-flexure-trench-v1` trench and outer rise. `None` is
/// `elevation_on_grid_with` exactly. Refused under the single-column law.
pub fn elevation_on_grid_terms(
    state: &CrustState,
    params: &PlanetParameterProfile,
    areas: &[f64],
    continental_term_wl_m: Option<&[f64]>,
    oceanic_term_wl_m: Option<&[f64]>,
) -> Result<(Vec<f64>, Option<SeaLevelSolve>), PreviewError> {
    // WS-CONST: law dispatch and every constant through the universe.
    let u = universe_of(params)?;
    if let Some(t) = oceanic_term_wl_m {
        if t.len() != state.cells() {
            return Err(PreviewError::Domain(format!(
                "trench term has {} cells, state {}",
                t.len(),
                state.cells()
            )));
        }
        if !matches!(u.law_set.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. }) {
            return Err(PreviewError::Domain(
                "trench term needs the two-column isostasy law".into(),
            ));
        }
    }
    if let Some(t) = continental_term_wl_m {
        if t.len() != state.cells() {
            return Err(PreviewError::Domain(format!(
                "rift flank term has {} cells, state {}",
                t.len(),
                state.cells()
            )));
        }
        if !matches!(u.law_set.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. }) {
            return Err(PreviewError::Domain(
                "rift flank term needs the two-column isostasy law".into(),
            ));
        }
    }
    let (fixed_depth_m, depth_margin_m) = match &u.law_set.isostasy {
        IsostasyLawV1::TwoColumnCommonDepthV0 {
            fixed_compensation_depth_m,
            derived_depth_margin_m,
        } => (
            fixed_compensation_depth_m.value,
            derived_depth_margin_m.as_ref().map(|m| m.value),
        ),
        _ => return Ok((elevation_single_column(state, &u)?, None)),
    };
    let columns: Vec<ColumnV0> = (0..state.cells())
        .map(|c| ColumnV0 {
            kind: state.kind[c],
            crust_m: state.thickness_m[c],
            sediment_m: state.sediment_m[c],
            thermal_uplift_water_loaded_m: if state.kind[c] == CrustKindCodeV0::Oceanic {
                match oceanic_term_wl_m {
                    Some(t) => gdh1_thermal_elevation_m(state.age_my(c)) + t[c],
                    None => gdh1_thermal_elevation_m(state.age_my(c)),
                }
            } else {
                match continental_term_wl_m {
                    Some(t) => t[c],
                    None => 0.0,
                }
            },
        })
        .collect();
    // Deepest crustal root (positive depth below the gauge zero), from the
    // same mass balance isostasy.rs uses: root = t_dry - m_layers / rho_m,
    // with sediment loading at the continental crust density (legacy
    // A-SEDIMENT-DENSITY) and t_dry the water-loaded thermal height times
    // (rho_m - rho_w) / rho_m. Independent of the compensation depth.
    let rho_m = u.constants.mantle_density_kg_m3.value;
    let rho_cc = u.constants.continental_crust_density_kg_m3.value;
    let rho_s = u.sediment_density_kg_m3();
    let rho_w = u.constants.seawater_density_kg_m3.value;
    let rho_oc = u
        .constants
        .oceanic_crust_density_kg_m3
        .as_ref()
        .ok_or_else(|| PreviewError::Params("universe: two-column law without oceanic crust density".into()))?
        .value;
    let (mut deepest, mut deepest_cell) = (f64::MIN, 0usize);
    for (c, col) in columns.iter().enumerate() {
        let rho_c = match col.kind {
            CrustKindCodeV0::Continental => rho_cc,
            CrustKindCodeV0::Oceanic => rho_oc,
        };
        let t_dry = col.thermal_uplift_water_loaded_m * (rho_m - rho_w) / rho_m;
        let depth = (rho_c * col.crust_m + rho_s * col.sediment_m) / rho_m - t_dry;
        if depth > deepest {
            deepest = depth;
            deepest_cell = c;
        }
    }
    let depth_m = match depth_margin_m {
        None => fixed_depth_m,
        Some(margin) => deepest.max(0.0) + margin,
    };
    let k = TwoColumnConstants::from_universe(&u, depth_m)
        .map_err(|e| PreviewError::Domain(format!("two-column isostasy: {e}")))?;
    let sol = two_column_field_with(&columns, areas, &k)
        .map_err(|e| PreviewError::Domain(format!("two-column isostasy: {e}")))?;
    let solve = SeaLevelSolve {
        sea_level_gauge_m: sol.sea_level_gauge_m,
        flooded_area_fraction: sol.flooded_area_fraction,
        ocean_volume_residual_m3: sol.volume_check_m3 - k.ocean_water_volume_m3,
        compensation_depth_m: depth_m,
        deepest_root_m: deepest,
        deepest_root_crust_sediment_m: [
            columns[deepest_cell].crust_m,
            columns[deepest_cell].sediment_m,
        ],
    };
    Ok((sol.elevation_m, Some(solve)))
}

/// Assumptions and approximation names the opt-in flags add or replace
/// (only called when a flag is on; the legacy lists are untouched).
fn flag_assumptions(
    p: &PlanetParameterProfile,
    prior: &UpliftPrior,
    assumptions: &mut Vec<AssumptionV0>,
    approximations: &mut Vec<String>,
) {
    let mk = |id: &str, statement: &str, value: f64, units: &str, source, note: &str| AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source,
        note: note.into(),
    };
    if let Some(b) = &p.two_column_isostasy {
        // These legacy entries describe the single-column model and are
        // false under the two-column model: removed, replaced below.
        assumptions.retain(|a| {
            !matches!(
                a.id.as_str(),
                "A-SEA-LEVEL-CONSTANT" | "A-GDH1-REFERENCE-OFFSET" | "A-WATER-DENSITY"
            )
        });
        for a in approximations.iter_mut() {
            if a == "airy-crust-only-isostasy-with-water-loading" {
                *a = "two-column-airy-isostasy-volume-primary-sea-level-global-fill".into();
            }
        }
        assumptions.push(mk("A-OCEAN-VOLUME", "liquid ocean water volume; sea level and flooded fraction are solved from it every isostasy evaluation (global fill)", b.ocean_water_volume_m3, "m3", AssumptionSourceV0::ProjectAssumption, "project assumption awaiting review; value recalled from Charette and Smith 2010, Oceanography 23(2):112 (1.335e9 km3), citation not yet checked"));
        assumptions.push(mk("A-OCEANIC-CRUST-DENSITY", "oceanic crust bulk density in the two-column balance", b.oceanic_crust_density_kg_m3, "kg/m3", AssumptionSourceV0::ProjectAssumption, "project assumption awaiting review; value recalled from Turcotte and Schubert, Geodynamics, and Carlson and Raskin 1984, Nature 311:555, citation not yet checked"));
        assumptions.push(mk("A-SEAWATER-DENSITY", "sea water density for water loading in the two-column balance", b.seawater_density_kg_m3, "kg/m3", AssumptionSourceV0::Literature, "Turcotte and Schubert 2014, Geodynamics, 3rd ed., section 2.2"));
        match b.compensation_depth_margin_m {
            None => assumptions.push(mk("A-COMPENSATION-DEPTH", "fixed common compensation depth below the mantle free-surface gauge; elevations are independent of it while it lies below every root", b.compensation_depth_m, "m", AssumptionSourceV0::ProjectAssumption, "drops out of the solution; checked per cell")),
            Some(m) => assumptions.push(mk("A-COMPENSATION-DEPTH-MARGIN", "compensation depth DERIVED at every isostasy evaluation as the deepest crustal root over all columns plus this margin (per-step value in the step log)", m, "m", AssumptionSourceV0::ProjectAssumption, "validity bound only: the depth cancels in the mass balance, so elevations do not depend on it; design review review, ADR 039 amendment C4/C5")),
        }
        assumptions.push(mk("A-SEA-LEVEL-VOLUME-PRIMARY", "published elevations are relative to the solved sea level, placed at the planet sea-level datum; the gauge sea level per step is in the step log", p.planet.sea_level_m, "m", AssumptionSourceV0::ProjectAssumption, "global fill: disconnected basins below sea level count as flooded; continental thermal term 0"));
        assumptions.push(mk("A-UPLIFT-RATE-KIND-DENSITY", "uplift-rate diagnostic for thickened cells uses (rho_m - rho_kind) / rho_m with the cell's crust-kind density", p.planet.mantle_density_kg_m3, "kg/m3", AssumptionSourceV0::ProjectAssumption, "diagnostic field only; elevation itself comes from the two-column solve"));
    }
    if let Some(b) = &p.continental_bodies {
        let c = b.resolved_area_fraction(&p.tectonics);
        assumptions.extend(prior.body_assumptions(
            c,
            b.continental_area_fraction.is_none(),
            b.rifted_margin_width_m,
        ));
        match &b.irregular_cratons {
            None => {
                approximations.push("geodesic-disc-continental-bodies-v0".into());
                approximations.push("linear-rifted-margin-taper-v0".into());
            }
            // SEED-IRREGULAR-SEEDS (lane IS): the craton law replaces both.
            Some(k) => {
                approximations.push("irregular-cratons-v1".into());
                approximations.push("two-segment-rifted-margin-v1".into());
                assumptions.push(mk("A-CRATON-INTERIOR-THICKNESS", "irregular-cratons-v1: initial continental assemblies (count = tectonics.continental_seeds, rank-size shares) of overlapping elliptic cratons with domain-warped, noise-roughened coasts ranked to the declared area; continental cells ride their assembly's plate; interior thickness = this base + multi-scale basement roughness + Gaussian suture roots on craton seams; two-segment rifted margin (necking wedge over 0.3 W, then a proximal shelf stretched by proximal_margin_stretch at the hinge)", k.craton_interior_thickness_m, "m", AssumptionSourceV0::ProjectAssumption, "lane IS 2026-10-13, PREREG.txt; model compensation for the missing continental lithosphere term: starts interiors near +0.4 km instead of +1.2 km; Earth cratonic crust is 35-45 km"));
                assumptions.push(mk("A-SUTURE-ROOT", "irregular-cratons-v1: old suture roots between initial cratons, crest thickening (m); Gaussian sigma in the law set", k.suture_crest_thickening_m, "m", AssumptionSourceV0::ProjectAssumption, "lane IS 2026-10-13; Appalachian/Ural-like worn sutures standing 1-2 km at t = 0; recalled, not checked"));
            }
        }
    }
}

/// SEED-OROGEN: assumptions and approximation names the new laws add or
/// replace (only called when one is selected; legacy lists untouched).
fn orogen_assumptions(
    belt: Option<&BeltParams>,
    subgrid_x_c_m: Option<f64>,
    relief_law: Option<&NonlinearSubgridParams>,
    ep: &ErosionParams,
    grid: &GridSpec,
    assumptions: &mut Vec<AssumptionV0>,
    approximations: &mut Vec<String>,
) {
    let pa = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source: AssumptionSourceV0::ProjectAssumption,
        note: note.into(),
    };
    let der = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source: AssumptionSourceV0::Derived,
        note: note.into(),
    };
    if let Some(b) = belt {
        // The legacy stacking entries describe a rule that did not run.
        assumptions.retain(|a| {
            !matches!(
                a.id.as_str(),
                "A-UPLIFT-ARC" | "A-UPLIFT-COLLISION-CAP" | "A-UPLIFT-WEDGE-PASS"
            )
        });
        for a in approximations.iter_mut() {
            if a == "kinematic-overlap-thickening-prior" {
                *a = "rate-limited-collision-belt-with-lower-crust-channel-flow-v1".into();
            }
        }
        approximations.push("arc-per-subducted-area-spread-on-overriding-plate-v1".into());
        assumptions.push(pa("A-BELT-SHEAR-DIP", "dip of the conjugate step-up shear zones of the two-sided crustal wedge; instantaneous belt width W = 2 h cot(dip), h the mean colliding crust thickness", b.shear_zone_dip_deg, "deg", "project assumption awaiting review; recalled from Willett, Beaumont and Fullsack 1993, Geology 21:371 (Coulomb shear zones ~30 deg), citation not checked; triangular weights over W, continental cells only, overlap cell alone when W/2 < cell spacing"));
        assumptions.push(pa("A-CHANNEL-VISCOSITY", "effective viscosity of the ductile lower-crustal channel; spreading flux q = rho_c g (1 - rho_c/rho_m) h_ch^3 / (12 mu) grad h (Poiseuille, isostatic surface slope)", b.channel_viscosity_pa_s, "Pa s", "project assumption awaiting review; recalled range 1e18 (Clark and Royden 2000, Geology 28:703) to 1e21 Pa s (Beaumont et al. 2001; Bird 1991), citations not checked; fixed a priori before any run"));
        assumptions.push(pa("A-CHANNEL-REFERENCE-THICKNESS", "crust thickness above which the excess forms the flowing channel, h_ch = max(0, h - h_ref); spreading only between continental cells", b.channel_reference_thickness_m, "m", "project assumption: normal-thickness continental crust does not flow laterally on these timescales"));
        assumptions.push(pa("A-DELAMINATION-THICKNESS", "crust above this thickness after spreading is removed each step (eclogitised-root delamination backstop; ledger term delaminated_m3)", b.delamination_thickness_m, "m", "backstop at the legacy 70 km cap (Zhang and Klemperer 2005; Kay and Kay 1993; recalled, not checked); its share of belt input is reported per step"));
        assumptions.push(der("A-ARC-PRODUCTIVITY", "arc crust added per unit area of subducted plate at convergent subduction samples", b.arc_productivity_m, "m", "derived from the legacy A-UPLIFT-ARC intent (0.005 x ~100 km arc = 500 m), i.e. 25 km3/km/My at 50 km/My, order of Jicha and Jagoutz 2015 (recalled, not checked)"));
        assumptions.push(pa("A-ARC-BELT-WIDTH", "arc addition spread with triangular weights over overriding-plate cells within this distance of the overlap cell (arc-trench gap not resolved)", b.arc_belt_width_m, "m", "project assumption awaiting review"));
        assumptions.push(der("A-SURFACE-GRAVITY-SPREADING", "surface gravity read by the lateral-spreading law", b.gravity_m_s2, "m/s2", "universe constant surface_gravity_m_s2 (stated, not derived from planet mass)"));
    }
    if let Some(x_c) = subgrid_x_c_m {
        for a in approximations.iter_mut() {
            if a == "implicit-stream-power-n1" {
                *a = "two-elevation-subgrid-stream-power-n1-v1".into();
            }
        }
        // tau_sub at the equator for the record.
        let dx = grid.cell_area_m2(grid.nlat() / 2).map(|a| a.sqrt()).unwrap_or(0.0);
        if let Some(rp) = relief_law {
            // SEED-RELIEF-LAW (lane RL): the n-generalised closure.
            for a in approximations.iter_mut() {
                if a == "two-elevation-subgrid-stream-power-n1-v1" {
                    *a = "two-elevation-subgrid-stream-power-nonlinear-v1".into();
                }
            }
            let lambda = subgrid_lambda(dx, x_c);
            assumptions.push(pa("A-SUBGRID-CHANNEL-HEAD", "length scale x_c of the smallest channels of the unresolved subgrid network; Lambda = integral_{x_c}^{dx} x^(-2m/n) dx = ln(dx/x_c) (m/n = 0.5), dx = sqrt(cell area)", x_c, "m", "project assumption awaiting review; recalled hillslope / channel-head lengths ~100 m (Montgomery and Dietrich 1988; Roering et al. 2007), not checked; enters through ln(dx/x_c)"));
            assumptions.push(pa("A-STREAM-POWER-N-NONLINEAR", "slope exponent n of stream-power-subgrid-nonlinear-v1, applied to the resolved trunk (dc/dt = -K_n A^(n/2) ((c - c_r)/L)^n, implicit) and to the subgrid closure (dz/dt = -K_n ((z - c)/Lambda)^n, backward Euler; quasi-steady self-similar network a(x) = x^2)", rp.slope_exponent_n, "1", "project assumption awaiting review; stochastic-threshold mean-flow n ~ 2 or more (Lague 2014; DiBiase and Whipple 2011; Harel et al. 2016; Baldwin et al. 2003), recalled, not checked; declared before any run (lane RL FALSIFIER.txt)"));
            assumptions.push(pa("A-REFERENCE-STEEPNESS", "channel steepness k_s = S A^0.5 at which the n-law erodes at the n = 1 law's rate: K_n = K k_ref^(1 - n)", rp.reference_steepness_m, "m", "project assumption awaiting review; active-range k_sn (Kirby and Whipple 2012) and K k_ref = 1 mm/yr (Montgomery and Brandon 2002), recalled, not checked; declared before any run"));
            assumptions.push(AssumptionV0 {
                id: "A-SUBGRID-LAMBDA".into(),
                statement: "Lambda = ln(dx/x_c) at the equator row (the subgrid network length factor; for n = 1 tau_sub = Lambda/K)".into(),
                value: lambda,
                units: "1".into(),
                source: AssumptionSourceV0::Derived,
                note: "derived from the cell size and A-SUBGRID-CHANNEL-HEAD".into(),
            });
            return;
        }
        let tau = subgrid_tau_my(dx, x_c, ep.k_per_my, ep.m);
        assumptions.push(pa("A-SUBGRID-CHANNEL-HEAD", "length scale x_c of the smallest channels of the unresolved subgrid network; cell-mean response time tau_sub = (1/K) integral_{x_c}^{dx} x^(-2m) dx with dx = sqrt(cell area)", x_c, "m", "project assumption awaiting review; recalled hillslope / channel-head lengths ~100 m (Montgomery and Dietrich 1988; Roering et al. 2007), not checked; enters through ln(dx/x_c)"));
        assumptions.push(pa("A-SUBGRID-NETWORK-GEOMETRY", "self-similar subgrid drainage network a(x) = x^2 (Hack exponent 2); tau_sub at the equator row in My is the value", tau, "My", "Hack's empirical exponent ~1.7 with mainstream sinuosity would lengthen tau_sub; cell mean relaxes toward a carried trunk level; effective per-cell contributing areas realise the target through the shared erosion stage"));
    }
}

fn area_weights(grid: &GridSpec) -> Result<Vec<f64>, PreviewError> {
    let mut out = Vec::with_capacity(grid.cells());
    for i in 0..grid.nlat() {
        let a = grid
            .cell_area_m2(i)
            .map_err(|e| PreviewError::Grid(e.to_string()))?;
        for _ in 0..grid.nlon() {
            out.push(a);
        }
    }
    Ok(out)
}

/// WS-CONST: K and D read through the universe constants; m, n and the
/// per-year unit of K are fixed by the erosion law version.
fn erosion_params(u: &UniverseV1) -> ErosionParams {
    ErosionParams {
        k_per_my: u.constants.stream_power_erodibility_per_yr.value * YEARS_PER_MY,
        m: STREAM_POWER_M,
        n: STREAM_POWER_N,
        diffusivity_m2_per_my: u.constants.hillslope_diffusivity_m2_per_yr.value * YEARS_PER_MY,
        stability_limit: DIFFUSION_STABILITY_LIMIT,
    }
}

/// WS-CONST: the preview uplift prior with its crust thicknesses read
/// through the universe constants (values identical to the prior's
/// defaults; equality pinned by tests/ws_const_preview.rs).
fn uplift_prior(u: &UniverseV1) -> UpliftPrior {
    let mut prior = UpliftPrior::preview_v0();
    prior.record.oceanic_crust_thickness_m = u.constants.oceanic_crust_thickness_m.value;
    prior.record.continental_crust_thickness_m = u.constants.continental_crust_thickness_m.value;
    prior
}

/// The labelled assumptions the driver itself adds (units, isostasy,
/// erosion constants, thermal model, what is missing).
fn driver_assumptions(p: &PlanetParameterProfile, ep: &ErosionParams) -> Vec<AssumptionV0> {
    let pa = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source: AssumptionSourceV0::ProjectAssumption,
        note: note.into(),
    };
    let lit = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source: AssumptionSourceV0::Literature,
        note: note.into(),
    };
    vec![
        pa("A-RAIN-UNITS", "SurfaceBlock.uniform_runoff_m_per_yr is effective runoff in m/yr, converted explicitly to m/My", p.surface.uniform_runoff_m_per_yr, "m/yr", "typed runoff units; legacy rainfall key refused"),
        pa("A-OCEAN-CONNECTIVITY", "largest-connected-subsea-component-by-area-v1; epoch-local reservoir, other depressions fill to spill; no water-inventory or salinity solution", 1.0, "reservoir", "equal-area ties choose lowest seed element; selection recorded in each step"),
        pa("A-SEA-LEVEL-CONSTANT", "prescribed sea level is constant; no sea-level history or coupled inland water loading", p.planet.sea_level_m, "m", "Airy water-loading still uses the prescribed sea datum independently of routing membership"),
        pa("A-GDH1-REFERENCE-OFFSET", "GDH1 thermal anomaly added to the project Airy reference, not anchored to observed normal ocean crust depth", -28000.0*(p.planet.mantle_density_kg_m3-p.planet.crust_density_kg_m3)/(p.planet.mantle_density_kg_m3-WATER_DENSITY_KG_M3)+GDH1_ASYMPTOTE_M, "m", "reference offset at default 7 km oceanic versus 35 km Airy column; retain as model convention"),
        pa("A-ARC-NORMAL-GEOMETRY", "beginning-of-step nearest interface edge for the exact plate pair; normal independent of velocity; nonlocal/missing interface refused", 3.0, "cell_diagonals", "locality radius = both plates swept angular distances plus three grid-cell diagonals; per-event source edge recorded in surface-exchange sidecar"),
        pa("A-SEDIMENT-SPREAD-BOUND", "maximum cells in each local receiving-water deposition search", 512.0, "cells", "the computational limit is not a physical capacity limit; hits recorded per epoch"),
        pa("A-MARINE-GLOBAL-LAYER", "remaining marine deposition is spread uniformly over the selected connected ocean reservoir after bounded local filling", 1.0, "enabled", "declared global redistribution shortcut, not marine transport; layer volume and sparse delivery locations recorded"),
        pa("A-K-UNITS", "SurfaceBlock.erodibility is interpreted as the stream-power K per YEAR (A in m^2, m = 0.5, n = 1) and converted to per My", p.surface.erodibility, "1/yr", "model-assumption-review Q2"),
        pa("A-STREAM-POWER-M", "stream-power area exponent m", STREAM_POWER_M, "1", "Whipple and Tucker 1999 range 0.3 to 0.6; model-assumption-review Q2"),
        pa("A-STREAM-POWER-N", "stream-power slope exponent n (implicit scheme requires 1)", STREAM_POWER_N, "1", "model-assumption-review Q2"),
        pa("A-HILLSLOPE-D", "linear hillslope diffusivity, applied at the global cell scale where it is a smoothing operator rather than a hillslope process", ep.diffusivity_m2_per_my, "m2/My", "model-assumption-review Q3"),
        pa("A-DIFFUSION-LIMIT", "explicit diffusion number limit per substep (convex-combination bound is 1.0; 0.5 declared)", DIFFUSION_STABILITY_LIMIT, "1", "model-assumption-review Q3"),
        lit("A-WATER-DENSITY", "sea water density for the water-loading factor rho_m / (rho_m - rho_w)", WATER_DENSITY_KG_M3, "kg/m3", "Turcotte and Schubert 2014, Geodynamics, 3rd ed., section 2.2 (Airy isostasy with water infill)"),
        lit("A-GDH1-ASYMPTOTE", "GDH1 plate-cooling asymptotic ocean depth; the thermal term is 5651 m minus GDH1 depth(age)", GDH1_ASYMPTOTE_M, "m", "Stein and Stein 1992, Nature 359:123, eq. 1 to 2 (d = 2600 + 365 sqrt(t) for t < 20 My; 5651 - 2473 exp(-0.0278 t) after)"),
        pa("A-SEDIMENT-DENSITY", "sediment loads the column at crust density (no compaction, no porosity)", p.planet.crust_density_kg_m3, "kg/m3", "single-density column keeps the built Airy function usable; model-assumption-review Q4"),
        pa("A-COMPENSATION", "Airy local compensation only; the bounded-support flexure kernel (flexure.rs) is not applied at this tier", 0.0, "1", "model-assumption-review Q4"),
        pa("A-ADVECTION", "semi-Lagrangian nearest-source carry of crust state; not conservative; residual reported in crust_volume", 0.0, "1", "model-assumption-review Q5"),
        pa("A-SED-COVER-THRESHOLD", "sediment thickness at which rock_family reports SedimentaryCover", SEDIMENTARY_COVER_THRESHOLD_M, "m", "coarse cell-scale threshold"),
        pa("A-NO-PARCEL-LEDGER", "the Lagrangian parcel ledger (parcels.rs) and S5a sediment provenance records are NOT driven at this tier; no material lineage is published", 0.0, "1", "returns with the solved-mantle tier (ADR 026)"),
        pa("A-NO-HOTSPOTS", "fixed mantle hotspots (K6) are not applied at this tier", 0.0, "1", "hotspot chains return with the solved-mantle upwellings"),
        pa("A-NO-WATERFALL-DETECTION", "waterfalls.rs is not run; the hydrology artifact carries an empty waterfall list", 0.0, "1", "ADR 026 s9 keeps the detector frozen; may be run on preview artifacts as an instrument later"),
    ]
}

/// Run the chain. Pure function of the config (no I/O). Legacy physics:
/// identical to `simulate_with(config, &PreviewOptions::default())`.
pub fn simulate(config: &PreviewConfig) -> Result<PreviewOutputs, PreviewError> {
    simulate_with(config, &PreviewOptions::default())
}

/// Run the chain with opt-in physics flags. Pure function (no I/O).
pub fn simulate_with(
    config: &PreviewConfig,
    options: &PreviewOptions,
) -> Result<PreviewOutputs, PreviewError> {
    let params = config.profile_with(options)?;
    run_profile(params, config.grid, config.seed_path())
}

/// WS-MODULE (lane SM, 2026-10-12): batch run of an already-resolved
/// profile, for example a recipe with a `pinned_universe` (WS-CONST idea B).
/// The profile is checked first (`check_profile`), then run exactly as
/// `simulate_with` runs it.
pub fn simulate_profile(
    params: PlanetParameterProfile,
    grid_override: Option<(usize, usize)>,
    seed_path: SeedPath,
) -> Result<PreviewOutputs, PreviewError> {
    check_profile(&params)?;
    run_profile(params, grid_override, seed_path)
}

/// Batch generation = init, every scheduled step, finish. The genesis
/// SimModule (`crate::preview_module`) drives the SAME three calls, so
/// batch generation is module replay with an empty intervention log
/// (ADR 039 §4, one implementation per tier).
fn run_profile(
    params: PlanetParameterProfile,
    grid_override: Option<(usize, usize)>,
    seed_path: SeedPath,
) -> Result<PreviewOutputs, PreviewError> {
    let mut run = PreviewRun::init(PreviewEnv::new(params, grid_override, seed_path)?)?;
    while run.steps_taken() < run.steps_scheduled() {
        run.step()?;
    }
    run.finish()
}

/// WS-MODULE: everything the chain derives from the recipe (resolved
/// profile, grid, seed path) before stage A. Pure function of the recipe;
/// rebuilt identically on restore, so it is never checkpointed.
pub struct PreviewEnv {
    pub(crate) grid_override: Option<(usize, usize)>,
    pub(crate) params: PlanetParameterProfile,
    pub(crate) universe: UniverseV1,
    pub(crate) grid: GridSpec,
    pub(crate) seed_path: SeedPath,
    pub(crate) prior: UpliftPrior,
    pub(crate) ep: ErosionParams,
    pub(crate) areas: Vec<f64>,
    pub(crate) sea: f64,
    pub(crate) rho_c: f64,
    pub(crate) airy_factor: f64,
    pub(crate) flags_on: bool,
    pub(crate) airy_factor_oceanic: Option<f64>,
    pub(crate) dt: f64,
    pub(crate) steps: u32,
    pub(crate) rain: Vec<f64>,
    pub(crate) stencil: DiffusionStencil,
    /// SEED-OROGEN: `rate-limited-belt-v1` parameters, if selected.
    pub(crate) belt: Option<BeltParams>,
    /// SEED-OROGEN: `stream-power-subgrid-v1` channel-head length, if selected.
    /// Also set under `stream-power-subgrid-nonlinear-v1` (lane RL), which
    /// carries the same subgrid relief state.
    pub(crate) subgrid_x_c_m: Option<f64>,
    /// SEED-RELIEF-LAW (lane RL): `stream-power-subgrid-nonlinear-v1`
    /// parameters, if selected (None = lane OG's n = 1 closure or legacy).
    pub(crate) relief_law: Option<NonlinearSubgridParams>,
    /// SEED-UPLANDS (lane UP): `critical-slope-envelope-v1` parameters and
    /// the grid's 8-neighbour adjacency, if selected.
    pub(crate) critical_slope: Option<(CriticalSlopeParams, Vec<Vec<(usize, f64)>>)>,
    /// SEED-SEDIMENT-RETAIN (lane SR): `continental-retention-base-level-v1`
    /// parameters, if selected.
    pub(crate) retention: Option<RetentionParams>,
    /// SHELF (lane SHELF, 2026-10-20): `wave-base-shelf-grading-v1`
    /// parameters, if `orogeny.sediment_retention.shelf_wave_grading`.
    pub(crate) shelf_grading: Option<pf_geodynamics::shelf_grade::ShelfGradeParams>,
    /// INLAND-BASINS (lane INLAND-BASINS, 2026-10-21):
    /// `closed-basin-water-balance-v1` parameters, if
    /// `orogeny.sediment_retention.closed_basin_water_balance`.
    pub(crate) closed_basin: Option<pf_surface::closed_basin::ClosedBasinParams>,
    /// INCISION-THRESHOLD (lane IT): `stream-power-excess-steepness-threshold-v1`
    /// parameters, if selected (needs `relief_law`).
    pub(crate) incision_threshold: Option<pf_surface::erosion::incision_threshold::IncisionThresholdParams>,
    /// GLACIAL-T0 (lane GLACIAL-T0, 2026-10-21): `glacial-stream-power-blend-v1`
    /// parameters, if selected (needs `relief_law` and `climate`).
    pub(crate) glacial: Option<pf_surface::erosion::glacial::GlacialParams>,
    /// INCISION-THRESHOLD part 1 (lane IT): the marine remainder layer law.
    pub(crate) ocean_layer: pf_surface::erosion::OceanRemainderLayer,
    /// SEED-EXTRAS D7 (lane SD): surface sediment laws (legacy = none).
    pub(crate) surface_laws: SurfaceLaws,
    /// SEED-CARRY-CC (lane CC): the plate carry law (legacy by default).
    pub(crate) carry: CarryLaw,
    /// SEDIMENT-TOWERS (lane ST): the collision-sediment law (legacy by default).
    pub(crate) collision_sediment: CollisionSedimentLaw,
    /// SEED-EXTRAS (lane EX): the seeded plumes and their parameters
    /// (`fixed-plume-chains-v1`), the rift parameters (`stretched-rift-v1`)
    /// and the split opening rate (`plume-head-continental-split-v1`).
    pub(crate) plumes: Option<(Vec<Plume>, PlumeParams)>,
    pub(crate) rift: Option<RiftParams>,
    pub(crate) split_rate_deg_per_my: Option<f64>,
    /// SEED-RIFT-SHOULDERS (lane RS): `shoulder-uplift-thermal-v1`.
    pub(crate) rift_flank: Option<RiftFlankParams>,
    /// SEED-TRENCH (lane TN): `slab-flexure-trench-v1`.
    pub(crate) trench: Option<TrenchParams>,
    /// SEED-TRENCH part 2 (lane TN): `seeded-margin-thermal-subsidence-v1`.
    pub(crate) seeded_margin_subsidence: bool,
    /// SEED-FLOODBASALT (lane FB): `plume-head-lip-swell-v1`.
    pub(crate) flood_basalt: Option<FloodBasaltParams>,
    /// SEED-CARRY-CONSERVE (lane CV): `continental-parcel-carry-v1`
    /// constants and the `contact-length-pole-share-v1` contact length (m).
    pub(crate) parcel_params: Option<ParcelParams>,
    /// HYPSO-LEVERS: write-back residual correction passes (0 = v1).
    pub(crate) parcel_residual_iterations: u32,
    pub(crate) suture_contact_m: Option<f64>,
    /// SEED-PLATE-REORG (lane PR): `plate-reorg-scripted-v1` constants and
    /// its re-solve schedule (seeded, prefix-stable in the run length).
    pub(crate) reorg: Option<(ReorgParams, Vec<f64>)>,
    /// SEED-BREAKUP-CUT-V2 (lane BC): which scripted law version `reorg`
    /// runs (ScriptedV1 unless `plate-reorg-scripted-v2` is selected).
    pub(crate) reorg_law: ReorgLaw,
    /// SEED-RANGES (lane RG): `contact-length-indentation-taper-v1` and
    /// `plateau-cap-spill-v1` constants.
    pub(crate) indent: Option<IndentParams>,
    pub(crate) plateau: Option<PlateauParams>,
    /// SEED-WORLD-COMPOSE (lane CP): `plateau-cap-spill-v2` (inside the
    /// belt step), the `plate-reorg-scripted-v3` initiation age, and under
    /// `slab-arc-backarc-llsvp-v2` the dry factor (rho_m - rho_w) / rho_m of
    /// the arc heat's dynamic (non-displacing) elevation term.
    pub(crate) plateau_in_belt: Option<pf_geodynamics::ranges::PlateauInBeltParams>,
    pub(crate) reorg_initiation_age_my: Option<f64>,
    pub(crate) magma_heat_dry_factor: Option<f64>,
    /// SEED-TRANSFORMS (lane TF): `smoothed-normal-transform-v1` constants
    /// and the `plate-reorg-scripted-v4` initiation limits (None = legacy).
    pub(crate) boundary_transform: Option<pf_geodynamics::boundaries::TransformParams>,
    pub(crate) reorg_limits: Option<pf_geodynamics::reorg::InitiationLimits>,
    /// SEED-OLD-OCEAN (lane OR): `plate-reorg-scripted-v5` margin initiation (None = v1-v4).
    pub(crate) reorg_margin: Option<pf_geodynamics::reorg::MarginInitiation>,
    /// PLATE-POLES (lane PLATE-POLES): `rigid-group-contacts-v1` release part
    /// (detached suture-group members leave before each re-solve).
    pub(crate) rigid_group_release: bool,
    /// SEED-MAGMATISM (lane MG): `slab-arc-backarc-llsvp-v1` constants, plume
    /// drifts and cell centroids.
    pub(crate) magma: Option<crate::magmatism::MagmaEnv>,
    /// SEED-STRATIGRAPHY (lane ST): `parcel-stacks-v1` constants.
    pub(crate) strata: Option<StrataParams>,
    /// CLIMATE-GENESIS (lane CG): `epoch-climate-runoff-v1`, if selected.
    pub(crate) climate: Option<GenesisClimateParams>,
}

/// SEED-OROGEN: the belt parameters a universe selects (None = legacy).
pub fn belt_params_of(u: &UniverseV1) -> Option<BeltParams> {
    match &u.law_set.tectonic_thickening {
        TectonicThickeningLawV1::KinematicOverlapStackingV0 => None,
        TectonicThickeningLawV1::RateLimitedBeltV1 {
            belt_shear_zone_dip_deg,
            lower_crust_channel_viscosity_pa_s,
            channel_reference_thickness_m,
            delamination_thickness_m,
            arc_productivity_m,
            arc_belt_width_m,
        } => Some(BeltParams {
            shear_zone_dip_deg: belt_shear_zone_dip_deg.value,
            channel_viscosity_pa_s: lower_crust_channel_viscosity_pa_s.value,
            channel_reference_thickness_m: channel_reference_thickness_m.value,
            delamination_thickness_m: delamination_thickness_m.value,
            arc_productivity_m: arc_productivity_m.value,
            arc_belt_width_m: arc_belt_width_m.value,
            gravity_m_s2: u.constants.surface_gravity_m_s2.value,
            crust_density_kg_m3: u.constants.continental_crust_density_kg_m3.value,
            mantle_density_kg_m3: u.constants.mantle_density_kg_m3.value,
        }),
    }
}

impl PreviewEnv {
    pub fn new(
        params: PlanetParameterProfile,
        grid_override: Option<(usize, usize)>,
        seed_path: SeedPath,
    ) -> Result<Self, PreviewError> {
        // WS-CONST: every physical constant and law choice below is read
        // through the resolved universe (legacy profile = the old values).
        let u = universe_of(&params)?;
        let grid = grid_for(grid_override, &params)?;
        let prior = uplift_prior(&u);
        let ep = erosion_params(&u);
        ep.validate()?;
        let n = grid.cells();
        let areas = area_weights(&grid)?;
        let sea = u.constants.sea_level_datum_m.value;
        let rho_c = u.constants.continental_crust_density_kg_m3.value;
        let rho_m = u.constants.mantle_density_kg_m3.value;
        let airy_factor = (rho_m - rho_c) / rho_m;
        let flags_on = !u.law_set.is_legacy();
        // Two-column mode: the uplift-rate diagnostic uses the crust-kind
        // density (legacy keeps the single factor above, byte-identical).
        let airy_factor_oceanic = match &u.law_set.isostasy {
            IsostasyLawV1::TwoColumnCommonDepthV0 { .. } => u
                .constants
                .oceanic_crust_density_kg_m3
                .as_ref()
                .map(|c| (rho_m - c.value) / rho_m),
            _ => None,
        };
        let dt = params.tectonics.step_my;
        let steps = (params.tectonics.duration_my / dt).round() as u32;
        let rain_m_per_my = params.surface.uniform_runoff_m_per_yr * YEARS_PER_MY;
        let rain = vec![rain_m_per_my; n];
        // CLIMATE-GENESIS (lane CG): the climate forcing law (r_ref = the
        // uniform runoff at which K is quoted).
        let climate = GenesisClimateParams::from_law(
            &u.law_set.climate_forcing,
            dt,
            &grid,
            params.surface.uniform_runoff_m_per_yr,
            ep.m,
            u.constants.surface_gravity_m_s2.value,
        )
        .map_err(PreviewError::Params)?;
        let stencil = DiffusionStencil::build(&grid)?;
        let belt = belt_params_of(&u);
        let (subgrid_x_c_m, relief_law) = match &u.law_set.erosion {
            ErosionLawV1::StreamPowerN1LinearDiffusionV0 => (None, None),
            ErosionLawV1::StreamPowerSubgridV1 {
                channel_head_length_m,
            } => (Some(channel_head_length_m.value), None),
            ErosionLawV1::StreamPowerSubgridNonlinearV1 {
                channel_head_length_m,
                slope_exponent_n,
                reference_steepness_m,
            } => {
                let rp = NonlinearSubgridParams {
                    k_per_my: ep.k_per_my,
                    slope_exponent_n: slope_exponent_n.value,
                    reference_steepness_m: reference_steepness_m.value,
                    channel_head_length_m: channel_head_length_m.value,
                };
                rp.validate().map_err(PreviewError::Domain)?;
                (Some(channel_head_length_m.value), Some(rp))
            }
        };
        let critical_slope = match &u.law_set.hillslope_threshold {
            HillslopeThresholdLawV1::NoneV0 => None,
            HillslopeThresholdLawV1::CriticalSlopeEnvelopeV1 {
                hillslope_critical_slope,
                knickpoint_celerity_m_per_my,
            } => {
                let cs = CriticalSlopeParams {
                    critical_slope: hillslope_critical_slope.value,
                    knickpoint_celerity_m_per_my: knickpoint_celerity_m_per_my.value,
                };
                cs.validate().map_err(PreviewError::Domain)?;
                let adj = pf_surface::uplands::adjacency(&grid).map_err(PreviewError::Domain)?;
                Some((cs, adj))
            }
        };
        // INCISION-THRESHOLD (lane IT; merged by lane CLIMATE-GENESIS-2):
        // the two opt-in laws.
        let ocean_layer = match u.law_set.marine_remainder_layer {
            pf_contracts::laws::MarineRemainderLayerLawV1::PerMouthGlobalLayerV0 => {
                pf_surface::erosion::OceanRemainderLayer::PerMouthV0
            }
            pf_contracts::laws::MarineRemainderLayerLawV1::OncePerStepCappedGlobalLayerV1 => {
                pf_surface::erosion::OceanRemainderLayer::OncePerStepV1
            }
        };
        let incision_threshold = match &u.law_set.fluvial_incision_threshold {
            pf_contracts::laws::FluvialIncisionThresholdLawV1::NoneV0 => None,
            pf_contracts::laws::FluvialIncisionThresholdLawV1::ExcessSteepnessThresholdV1 {
                critical_channel_steepness_m,
            }
            // LAKE-HOLD-2: the trunk-only variant (same constant).
            | pf_contracts::laws::FluvialIncisionThresholdLawV1::ExcessSteepnessThresholdTrunkOnlyV1 {
                critical_channel_steepness_m,
            } => {
                let t = pf_surface::erosion::incision_threshold::IncisionThresholdParams {
                    critical_channel_steepness_m: critical_channel_steepness_m.value,
                    trunk_only: matches!(
                        u.law_set.fluvial_incision_threshold,
                        pf_contracts::laws::FluvialIncisionThresholdLawV1::ExcessSteepnessThresholdTrunkOnlyV1 { .. }
                    ),
                    // ALLUVIAL-TRUNK: set below, once the retention law's S_a is known.
                    alluvial_grade_slope: None,
                };
                t.validate().map_err(PreviewError::Domain)?;
                Some(t)
            }
        };
        // GLACIAL-T0 (lane GLACIAL-T0): the opt-in glacial erosion law.
        let glacial = match &u.law_set.glacial_erosion {
            pf_contracts::laws::GlacialErosionLawV1::NoneV0 => None,
            pf_contracts::laws::GlacialErosionLawV1::GlacialStreamPowerBlendV1 {
                glacial_to_fluvial_erodibility_ratio,
                full_ice_height_above_ela_m,
                coldhouse_duty,
                glacial_max_cooling_k,
                glacial_max_precipitation_factor,
                warm_state_duties,
                warm_state_warmings_k,
            } => {
                let g = pf_surface::erosion::glacial::GlacialParams {
                    glacial_to_fluvial_erodibility_ratio: glacial_to_fluvial_erodibility_ratio.value,
                    full_ice_height_above_ela_m: full_ice_height_above_ela_m.value,
                    coldhouse_duty: coldhouse_duty.value,
                    glacial_max_cooling_k: glacial_max_cooling_k.value,
                    glacial_max_precipitation_factor: glacial_max_precipitation_factor.value,
                    warm_states: warm_state_duties
                        .iter()
                        .zip(warm_state_warmings_k)
                        .map(|(d, w)| pf_surface::erosion::glacial::WarmStateParams { duty: d.value, warming_k: w.value })
                        .collect(),
                };
                g.validate().map_err(PreviewError::Domain)?;
                if relief_law.is_none() || climate.is_none() {
                    return Err(PreviewError::Params(format!(
                        "{}: needs stream-power-subgrid-nonlinear-v1 and epoch-climate-runoff-v1",
                        pf_surface::erosion::glacial::GLACIAL_EROSION_LAW_ID
                    )));
                }
                Some(g)
            }
        };
        // SEED-SEDIMENT-RETAIN (lane SR): the retention law, if selected
        // (compensated column with the isostasy law's own densities).
        let retention = match &u.law_set.sediment_retention {
            SedimentRetentionLawV1::NoneV0 => None,
            SedimentRetentionLawV1::ContinentalRetentionBaseLevelV1 {
                base_level_relief_scale_m,
                base_level_taper_exponent,
                alluvial_grade_slope,
                shelf_wave_base_depth_m,
            } => {
                let rp = RetentionParams {
                    base_level_relief_scale_m: base_level_relief_scale_m.value,
                    taper_exponent: base_level_taper_exponent.value,
                    alluvial_grade_slope: alluvial_grade_slope.value,
                    shelf_wave_base_depth_m: shelf_wave_base_depth_m.value,
                    ocean_layer,
                    ..RetentionParams::v1(CompensatedColumn {
                        mantle_density_kg_m3: u.constants.mantle_density_kg_m3.value,
                        sediment_density_kg_m3: u.sediment_density_kg_m3(),
                        water_density_kg_m3: u.constants.seawater_density_kg_m3.value,
                    })
                };
                rp.validate().map_err(|e| PreviewError::Domain(e.to_string()))?;
                Some(rp)
            }
        };
        // ALLUVIAL-TRUNK (`stream-power-threshold-free-alluvial-channels-v1`):
        // S_a of the retention law releases the threshold on alluvial trunks
        // (here) and on alluvial lake outlet channels (rule 10, below).
        let alluvial_grade_slope = match &u.law_set.fluvial_alluvial_channel {
            pf_contracts::laws::FluvialAlluvialChannelLawV1::NoneV0 => None,
            pf_contracts::laws::FluvialAlluvialChannelLawV1::ThresholdFreeAlluvialChannelsV1 => match (&retention, &incision_threshold) {
                (Some(rp), Some(_)) => Some(rp.alluvial_grade_slope),
                _ => {
                    return Err(PreviewError::Domain(format!(
                        "{} needs an incision threshold and continental-retention-base-level-v1",
                        pf_surface::erosion::incision_threshold::ALLUVIAL_CHANNEL_LAW_ID
                    )))
                }
            },
        };
        let incision_threshold = match incision_threshold {
            Some(t) if alluvial_grade_slope.is_some() => {
                let t = pf_surface::erosion::incision_threshold::IncisionThresholdParams { alluvial_grade_slope, ..t };
                t.validate().map_err(PreviewError::Domain)?;
                Some(t)
            }
            t => t,
        };
        // SHELF: read from the resolved law set (law-set entry, lane
        // CLIMATE-GENESIS-2), so a pinned universe governs it like every
        // other law.
        let shelf_grading = match (&retention, &u.law_set.shelf_grading) {
            (_, ShelfGradingLawV1::NoneV0) => None,
            (Some(rp), ShelfGradingLawV1::WaveBaseShelfGradingV1) => {
                let g = pf_geodynamics::shelf_grade::ShelfGradeParams {
                    wave_base_depth_m: rp.shelf_wave_base_depth_m,
                    mantle_density_kg_m3: rp.column.mantle_density_kg_m3,
                    sediment_density_kg_m3: rp.column.sediment_density_kg_m3,
                    seawater_density_kg_m3: rp.column.water_density_kg_m3,
                    search_max_cells: rp.search_max_cells,
                };
                g.validate().map_err(PreviewError::Domain)?;
                Some(g)
            }
            (None, ShelfGradingLawV1::WaveBaseShelfGradingV1) => {
                return Err(PreviewError::Domain("shelf_wave_grading needs continental-retention-base-level-v1".into()))
            }
        };
        // INLAND-BASINS: read from the resolved law set.
        let closed_basin = match (&retention, &u.law_set.closed_basin) {
            (_, ClosedBasinLawV1::NoneV0) => None,
            (Some(_), ClosedBasinLawV1::WaterBalanceFocusedFillV1 { suspended_load_fraction, uniform_evaporation_to_runoff_ratio }) => {
                let cb = pf_surface::closed_basin::ClosedBasinParams {
                    suspended_load_fraction: suspended_load_fraction.value,
                    uniform_evaporation_to_runoff_ratio: uniform_evaporation_to_runoff_ratio.value,
                    filled_basin_grading_depth_m: None,
                    outlet: None,
                };
                cb.validate().map_err(PreviewError::Domain)?;
                Some(cb)
            }
            // LAKE-EXCESS: the graded variant (`filled-basin-alluvial-grade-v1`).
            (Some(_), ClosedBasinLawV1::WaterBalanceFocusedFillAlluvialGradeV1 { suspended_load_fraction, uniform_evaporation_to_runoff_ratio, filled_basin_grading_depth_m }) => {
                let cb = pf_surface::closed_basin::ClosedBasinParams {
                    suspended_load_fraction: suspended_load_fraction.value,
                    uniform_evaporation_to_runoff_ratio: uniform_evaporation_to_runoff_ratio.value,
                    filled_basin_grading_depth_m: Some(filled_basin_grading_depth_m.value),
                    outlet: None,
                };
                cb.validate().map_err(PreviewError::Domain)?;
                Some(cb)
            }
            (None, _) => {
                return Err(PreviewError::Domain("closed_basin_water_balance needs continental-retention-base-level-v1".into()))
            }
        };
        // LAKE-OUTLET (`lake-outlet-sill-incision-v1`): the outlet channel
        // law, from laws already selected (no constant of its own): the
        // relief law's trunk law with the incision threshold if selected,
        // else the retention step's linear law; K quoted at the climate
        // law's r_ref, else at the uniform runoff.
        let closed_basin = match (&u.law_set.lake_outlet, closed_basin) {
            (pf_contracts::laws::LakeOutletLawV1::NoneV0, cb) => cb,
            (pf_contracts::laws::LakeOutletLawV1::SillIncisionV1
                | pf_contracts::laws::LakeOutletLawV1::SillIncisionFlexuralSupportV1
                | pf_contracts::laws::LakeOutletLawV1::PathIncisionV1
                | pf_contracts::laws::LakeOutletLawV1::PathIncisionFlexuralSupportV1
                | pf_contracts::laws::LakeOutletLawV1::PathThroughLakesV1
                | pf_contracts::laws::LakeOutletLawV1::PathThroughLakesFlexuralSupportV1, Some(mut cb)) => {
                let channel = relief_law.as_ref().map(|rp| pf_surface::closed_basin::OutletChannelLaw {
                    k_n_per_my: rp.k_n_per_my(),
                    m: rp.m(),
                    n: rp.slope_exponent_n,
                    threshold_rate_m_per_my: incision_threshold.as_ref().map_or(0.0, |t| t.threshold_rate_m_per_my(rp)),
                });
                let o = pf_surface::closed_basin::OutletIncisionParams {
                    reference_runoff_m_per_yr: climate.as_ref().map_or(params.surface.uniform_runoff_m_per_yr, |c| c.r_ref_m_per_yr),
                    channel,
                    // SILL-REBOUND: rule 6 (no local rebound at sills).
                    flexural_support: matches!(
                        u.law_set.lake_outlet,
                        pf_contracts::laws::LakeOutletLawV1::SillIncisionFlexuralSupportV1
                            | pf_contracts::laws::LakeOutletLawV1::PathIncisionFlexuralSupportV1
                            | pf_contracts::laws::LakeOutletLawV1::PathThroughLakesFlexuralSupportV1
                    ),
                    // LAKE-HOLD: rule 7 (the whole overflow path is incised).
                    path_incision: matches!(
                        u.law_set.lake_outlet,
                        pf_contracts::laws::LakeOutletLawV1::PathIncisionV1
                            | pf_contracts::laws::LakeOutletLawV1::PathIncisionFlexuralSupportV1
                            | pf_contracts::laws::LakeOutletLawV1::PathThroughLakesV1
                            | pf_contracts::laws::LakeOutletLawV1::PathThroughLakesFlexuralSupportV1
                    ),
                    // LAKE-HOLD-2: rule 8 (the cut passes through lakes on
                    // the path within one step).
                    through_lakes: matches!(
                        u.law_set.lake_outlet,
                        pf_contracts::laws::LakeOutletLawV1::PathThroughLakesV1
                            | pf_contracts::laws::LakeOutletLawV1::PathThroughLakesFlexuralSupportV1
                    ),
                    // LAKE-HOLD-2: rule 9 (threshold grade on the path).
                    path_grade: matches!(
                        u.law_set.lake_outlet_grade,
                        pf_contracts::laws::LakeOutletGradeLawV1::PathThresholdGradeV1
                    ),
                    // ALLUVIAL-TRUNK: rule 10 (no threshold on alluvial outlet channels).
                    alluvial_grade_slope,
                };
                o.validate().map_err(PreviewError::Domain)?;
                cb.outlet = Some(o);
                Some(cb)
            }
            (_, None) => {
                return Err(PreviewError::Domain("lake_outlet_incision needs closed-basin-water-balance-v1".into()))
            }
        };
        let surface_laws = surface_laws_of(&u);
        let mut parcel_residual_iterations = 0u32;
        let (carry, parcel_params) = match &u.law_set.plate_carry {
            PlateCarryLawV1::SemiLagrangianNearestMaterialV0 => {
                (CarryLaw::SemiLagrangianNearestMaterialV0, None)
            }
            PlateCarryLawV1::FrontRescueCarryV1 => (CarryLaw::FrontRescueV1, None),
            PlateCarryLawV1::ContinentalParcelCarryV1 {
                parcel_overlap_cell_fraction,
                plate_fragment_min_area_fraction,
                write_back_residual_iterations,
            } => {
                parcel_residual_iterations = write_back_residual_iterations.as_ref().map_or(0, |n| n.value as u32);
                let pp = ParcelParams {
                    overlap_cell_fraction: parcel_overlap_cell_fraction.value,
                    fragment_min_area_fraction: plate_fragment_min_area_fraction.value,
                };
                pp.validate().map_err(PreviewError::Domain)?;
                (CarryLaw::ContinentalParcelsV1, Some(pp))
            }
        };
        let suture_contact_m = match &u.law_set.plate_suture {
            PlateSutureLawV1::NoneV0 => None,
            PlateSutureLawV1::ContactLengthPoleShareV1 {
                suture_contact_length_m,
            } => Some(suture_contact_length_m.value),
            // SEED-RANGES (lane RG): same contact length, plus indentation.
            PlateSutureLawV1::ContactLengthIndentationTaperV1 {
                suture_contact_length_m,
                ..
            } => Some(suture_contact_length_m.value),
        };
        let indent = match &u.law_set.plate_suture {
            PlateSutureLawV1::ContactLengthIndentationTaperV1 {
                indentation_time_my,
                indentation_shortening_budget_m,
                ..
            } => Some(IndentParams {
                time_my: indentation_time_my.value,
                budget_m: indentation_shortening_budget_m.value,
            }),
            _ => None,
        };
        let plateau = match &u.law_set.orogen_plateau {
            OrogenPlateauLawV1::NoneV0 => None,
            OrogenPlateauLawV1::PlateauCapSpillV1 {
                plateau_crust_cap_m,
                plateau_spill_reach_m,
            } => Some(PlateauParams {
                cap_m: plateau_crust_cap_m.value,
                reach_m: plateau_spill_reach_m.value,
            }),
            OrogenPlateauLawV1::PlateauCapSpillV2 { .. } => None,
        };
        // SEED-WORLD-COMPOSE (lane CP): `plateau-cap-spill-v2`.
        let plateau_in_belt = match &u.law_set.orogen_plateau {
            OrogenPlateauLawV1::PlateauCapSpillV2 {
                plateau_crust_cap_m,
                plateau_spill_reach_m,
                plateau_foothill_gradient,
                plateau_taper_floor_m,
            } => Some(pf_geodynamics::ranges::PlateauInBeltParams {
                cap_m: plateau_crust_cap_m.value,
                reach_m: plateau_spill_reach_m.value,
                foothill_gradient: plateau_foothill_gradient.value,
                taper_floor_m: plateau_taper_floor_m.value,
            }),
            _ => None,
        };
        let collision_sediment = match u.law_set.collision_sediment {
            CollisionSedimentLawV1::SurfaceLayerV0 => CollisionSedimentLaw::SurfaceLayerV0,
            CollisionSedimentLawV1::AccretedToBeltCrustV1 => {
                CollisionSedimentLaw::AccretedToBeltCrustV1
            }
        };
        // SEED-WORLD-COMPOSE (lane CP): v2 / v3 read through v1's constants.
        let reorg_view = u.law_set.plate_reorganization.shared_view();
        let plumes = match &u.law_set.hotspot_volcanism {
            HotspotVolcanismLawV1::NoneV0 => None,
            HotspotVolcanismLawV1::FixedPlumeChainsV1 {
                plume_count,
                plume_onset_lead_my,
                plume_head_volume_m3,
                plume_head_radius_m,
                plume_tail_flux_m3_per_my,
                plume_edifice_radius_m,
            } => Some((
                // SEED-PLATE-REORG (lane PR): under plate-reorg-scripted-v1
                // the onsets are a renewal process (run-length independent).
                match &reorg_view {
                    PlateReorganizationLawV1::PlateReorgScriptedV1 { plume_head_interval_my, .. }
                    | PlateReorganizationLawV1::PlateReorgScriptedV2 { plume_head_interval_my, .. } => {
                        seed_plumes_renewal(
                            plume_head_interval_my.value,
                            plume_onset_lead_my.value,
                            params.tectonics.duration_my,
                            &seed_path,
                        )
                    }
                    _ => seed_plumes(
                        plume_count.value as u32,
                        plume_onset_lead_my.value,
                        params.tectonics.duration_my,
                        &seed_path,
                    ),
                },
                PlumeParams {
                    head_volume_m3: plume_head_volume_m3.value,
                    head_radius_m: plume_head_radius_m.value,
                    tail_flux_m3_per_my: plume_tail_flux_m3_per_my.value,
                    edifice_radius_m: plume_edifice_radius_m.value,
                },
            )),
        };
        // SEED-MAGMATISM (lane MG): `slab-arc-backarc-llsvp-v1`; under it the
        // plumes are re-drawn at the LLSVP edges (ids and onsets unchanged).
        let mut plumes = plumes;
        // SEED-WORLD-COMPOSE (lane CP): v2 reads v1's constants plus its own.
        let magma_v2 = match &u.law_set.magmatism {
            pf_contracts::laws::MagmatismLawV1::SlabArcBackarcLlsvpV2 {
                backarc_episode_my,
                backarc_max_stretch_factor,
                kimberlite_plume_radius_m,
                ..
            }
            | pf_contracts::laws::MagmatismLawV1::SlabArcBackarcLlsvpV3 {
                backarc_episode_my,
                backarc_max_stretch_factor,
                kimberlite_plume_radius_m,
                ..
            } => Some(crate::magmatism::MagmaV2Params {
                backarc_episode_my: backarc_episode_my.value,
                backarc_max_stretch_factor: backarc_max_stretch_factor.value,
                kimberlite_plume_radius_m: kimberlite_plume_radius_m.value,
            }),
            _ => None,
        };
        // Only v2 keeps the historical non-displacing post-solve convention.
        // V3 keeps v2 episodes/resources but uses the column heat path below.
        let magma_heat_dry_factor = if matches!(&u.law_set.magmatism,
            pf_contracts::laws::MagmatismLawV1::SlabArcBackarcLlsvpV2 { .. }) {
            Some((rho_m - u.constants.seawater_density_kg_m3.value) / rho_m)
        } else { None };
        let magma = match &u.law_set.magmatism.shared_view() {
            pf_contracts::laws::MagmatismLawV1::NoneV0 => None,
            pf_contracts::laws::MagmatismLawV1::SlabArcBackarcLlsvpV2 { .. }
            | pf_contracts::laws::MagmatismLawV1::SlabArcBackarcLlsvpV3 { .. } => unreachable!("shared_view maps v2/v3 to v1"),
            pf_contracts::laws::MagmatismLawV1::SlabArcBackarcLlsvpV1 {
                slab_depth_below_arc_m,
                slab_dip_young_deg,
                slab_dip_old_deg,
                slab_dip_age_my,
                arc_magma_productivity_m,
                arc_rate_floor_m_per_my,
                arc_rate_cap_m_per_my,
                arc_half_width_m,
                arc_continental_thickness_m,
                arc_thermal_uplift_m,
                arc_thermal_half_width_m,
                arc_thermal_decay_my,
                backarc_slab_age_min_my,
                backarc_overriding_advance_max_m_per_my,
                backarc_extension_rate_m_per_my,
                backarc_axis_offset_m,
                backarc_half_width_m,
                backarc_min_crust_m,
                llsvp_edge_deg,
                llsvp_edge_sigma_deg,
                llsvp_antipode_jitter_deg,
                plume_background_fraction,
                plume_drift_rate_m_per_my,
                vent_active_window_my,
                volcano_spacing_m,
                arc_crust_cap_m,
            } => {
                let mp = crate::magmatism::MagmaParams {
                    arc: pf_geodynamics::arcs::ArcParams {
                        slab_depth_m: slab_depth_below_arc_m.value,
                        dip_young_deg: slab_dip_young_deg.value,
                        dip_old_deg: slab_dip_old_deg.value,
                        dip_age_my: slab_dip_age_my.value,
                        productivity_m: arc_magma_productivity_m.value,
                        rate_floor_m_per_my: arc_rate_floor_m_per_my.value,
                        rate_cap_m_per_my: arc_rate_cap_m_per_my.value,
                        half_width_m: arc_half_width_m.value,
                        continental_thickness_m: arc_continental_thickness_m.value,
                        thermal_uplift_m: arc_thermal_uplift_m.value,
                        thermal_sigma_m: arc_thermal_half_width_m.value,
                        thermal_decay_my: arc_thermal_decay_my.value,
                        crust_cap_m: arc_crust_cap_m.value,
                    },
                    backarc: pf_geodynamics::backarc::BackarcParams {
                        slab_age_min_my: backarc_slab_age_min_my.value,
                        advance_max_m_per_my: backarc_overriding_advance_max_m_per_my.value,
                        rate_m_per_my: backarc_extension_rate_m_per_my.value,
                        axis_offset_m: backarc_axis_offset_m.value,
                        half_width_m: backarc_half_width_m.value,
                        min_crust_m: backarc_min_crust_m.value,
                    },
                    source: pf_geodynamics::hotspots::PlumeSourceParams {
                        edge_deg: llsvp_edge_deg.value,
                        edge_sigma_deg: llsvp_edge_sigma_deg.value,
                        antipode_jitter_deg: llsvp_antipode_jitter_deg.value,
                        background_fraction: plume_background_fraction.value,
                        drift_rate_m_per_my: plume_drift_rate_m_per_my.value,
                    },
                    vent_active_window_my: vent_active_window_my.value,
                    volcano_spacing_m: volcano_spacing_m.value,
                    v2: magma_v2,
                };
                let head_radius_m = plumes.as_ref().map(|p| p.1.head_radius_m).unwrap_or(0.0);
                let tail_lifetime_my = match &reorg_view {
                    PlateReorganizationLawV1::PlateReorgScriptedV1 { plume_tail_lifetime_my, .. }
                    | PlateReorganizationLawV1::PlateReorgScriptedV2 { plume_tail_lifetime_my, .. } => Some(plume_tail_lifetime_my.value),
                    _ => None,
                };
                Some(
                    crate::magmatism::MagmaEnv::new(mp, &grid, &seed_path, plumes.as_mut().map(|p| &mut p.0), head_radius_m, tail_lifetime_my)
                        .map_err(PreviewError::Params)?,
                )
            }
        };
        let rift = match &u.law_set.continental_rift {
            ContinentalRiftLawV1::OceanicBirthV0 => None,
            ContinentalRiftLawV1::StretchedRiftV1 {
                rift_half_width_m,
                breakup_crust_thickness_m,
            } => Some(RiftParams {
                half_width_m: rift_half_width_m.value,
                breakup_thickness_m: breakup_crust_thickness_m.value,
            }),
        };
        let split_rate_deg_per_my = match &u.law_set.plate_reorganization {
            PlateReorganizationLawV1::FixedPlatesV0 => None,
            PlateReorganizationLawV1::PlumeHeadContinentalSplitV1 {
                rift_opening_rate_deg_per_my,
            } => Some(rift_opening_rate_deg_per_my.value),
            // SEED-PLATE-REORG: plume splits run through the reorg law.
            PlateReorganizationLawV1::PlateReorgScriptedV1 { .. }
            | PlateReorganizationLawV1::PlateReorgScriptedV2 { .. }
            | PlateReorganizationLawV1::PlateReorgScriptedV3 { .. }
            | PlateReorganizationLawV1::PlateReorgScriptedV4 { .. }
            | PlateReorganizationLawV1::PlateReorgScriptedV5 { .. } => None,
        };
        // SEED-PLATE-REORG (lane PR): `plate-reorg-scripted-v1`.
        // SEED-BREAKUP-CUT-V2 (lane BC): v1 and v2 share the constants.
        let reorg_law = match &u.law_set.plate_reorganization {
            PlateReorganizationLawV1::PlateReorgScriptedV2 { .. } => ReorgLaw::ScriptedV2,
            PlateReorganizationLawV1::PlateReorgScriptedV4 { .. } => ReorgLaw::ScriptedV4,
            PlateReorganizationLawV1::PlateReorgScriptedV5 { .. } => ReorgLaw::ScriptedV5,
            PlateReorganizationLawV1::PlateReorgScriptedV3 { .. } => ReorgLaw::ScriptedV3,
            _ => ReorgLaw::ScriptedV1,
        };
        let reorg_initiation_age_my = match &u.law_set.plate_reorganization {
            PlateReorganizationLawV1::PlateReorgScriptedV3 { subduction_initiation_age_my, .. } => Some(subduction_initiation_age_my.value),
            PlateReorganizationLawV1::PlateReorgScriptedV4 { subduction_initiation_age_my, .. } => Some(subduction_initiation_age_my.value),
            PlateReorganizationLawV1::PlateReorgScriptedV5 { subduction_initiation_age_my, .. } => Some(subduction_initiation_age_my.value),
            _ => None,
        };
        // SEED-TRANSFORMS (lane TF): plate-reorg-scripted-v4 limits and the
        // plate-boundary law (None = legacy).
        let reorg_limits = match &u.law_set.plate_reorganization {
            PlateReorganizationLawV1::PlateReorgScriptedV4 {
                initiation_slab_age_min_my,
                initiation_old_slab_fraction_min,
                initiation_trench_spacing_m,
                ..
            } => Some(pf_geodynamics::reorg::InitiationLimits {
                slab_age_min_my: initiation_slab_age_min_my.value,
                old_slab_fraction_min: initiation_old_slab_fraction_min.value,
                trench_spacing_m: initiation_trench_spacing_m.value,
            }),
            _ => None,
        };
        // SEED-OLD-OCEAN (lane OR): plate-reorg-scripted-v5 carries v4's
        // limits plus the passive-margin initiation (None = v1-v4).
        let reorg_limits = match &u.law_set.plate_reorganization {
            PlateReorganizationLawV1::PlateReorgScriptedV5 {
                initiation_slab_age_min_my,
                initiation_old_slab_fraction_min,
                initiation_trench_spacing_m,
                ..
            } => Some(pf_geodynamics::reorg::InitiationLimits {
                slab_age_min_my: initiation_slab_age_min_my.value,
                old_slab_fraction_min: initiation_old_slab_fraction_min.value,
                trench_spacing_m: initiation_trench_spacing_m.value,
            }),
            _ => reorg_limits,
        };
        let reorg_margin = match &u.law_set.plate_reorganization {
            PlateReorganizationLawV1::PlateReorgScriptedV5 {
                margin_initiation_age_my,
                margin_segment_length_m,
                margin_trench_spacing_m,
                ..
            } => Some(pf_geodynamics::reorg::MarginInitiation {
                age_my: margin_initiation_age_my.value,
                segment_length_m: margin_segment_length_m.value,
                spacing_m: margin_trench_spacing_m.value,
            }),
            _ => None,
        };
        let rigid_group_release = u.law_set.rigid_group == pf_contracts::laws::RigidGroupLawV1::RigidGroupContactsV1;
        let boundary_transform = match &u.law_set.plate_boundary {
            PlateBoundaryLawV1::NeighbourDirectionV0 => None,
            PlateBoundaryLawV1::SmoothedNormalTransformV1 {
                transform_obliquity_deg,
                boundary_normal_smoothing_m,
            } => Some(pf_geodynamics::boundaries::TransformParams {
                transform_obliquity_deg: transform_obliquity_deg.value,
                normal_smoothing_m: boundary_normal_smoothing_m.value,
                // PLATE-POLES (lane PLATE-POLES): `rigid-group-interior-v1` and
                // `rigid-group-contacts-v1` both classify welded contacts as interior.
                rigid_group_interior: !u.law_set.rigid_group.is_legacy(),
            }),
        };
        let reorg = match &reorg_view {
            PlateReorganizationLawV1::PlateReorgScriptedV1 {
                rift_opening_rate_deg_per_my,
                reorg_interval_min_my,
                reorg_interval_max_my,
                slab_pull_bias,
                slab_circumference_reference,
                reorg_rate_floor_fraction,
                supercontinent_share,
                supercontinent_tenure_my,
                breakup_piece_max_share,
                breakup_cascade_my,
                rift_persistence_my,
                weak_zone_strength_ratio,
                plate_count_min,
                plate_count_max,
                plate_count_min_area_fraction,
                subduction_initiation_rate_deg_per_my,
                plume_head_interval_my,
                plume_tail_lifetime_my,
                plume_split_min_share,
            }
            | PlateReorganizationLawV1::PlateReorgScriptedV2 {
                rift_opening_rate_deg_per_my,
                reorg_interval_min_my,
                reorg_interval_max_my,
                slab_pull_bias,
                slab_circumference_reference,
                reorg_rate_floor_fraction,
                supercontinent_share,
                supercontinent_tenure_my,
                breakup_piece_max_share,
                breakup_cascade_my,
                rift_persistence_my,
                weak_zone_strength_ratio,
                plate_count_min,
                plate_count_max,
                plate_count_min_area_fraction,
                subduction_initiation_rate_deg_per_my,
                plume_head_interval_my,
                plume_tail_lifetime_my,
                plume_split_min_share,
            } => {
                let head_radius = match &plumes {
                    Some((_, pp)) => pp.head_radius_m,
                    None => {
                        return Err(PreviewError::Domain(
                            "plate-reorg-scripted-v1 requires fixed-plume-chains-v1".into(),
                        ))
                    }
                };
                if parcel_params.is_none() || suture_contact_m.is_none() || rift.is_none() {
                    return Err(PreviewError::Domain(
                        "plate-reorg-scripted-v1 requires continental-parcel-carry-v1, contact-length-pole-share-v1 and stretched-rift-v1".into(),
                    ));
                }
                let rp = ReorgParams {
                    rift_opening_rate_deg_per_my: rift_opening_rate_deg_per_my.value,
                    interval_min_my: reorg_interval_min_my.value,
                    interval_max_my: reorg_interval_max_my.value,
                    slab_pull_bias: slab_pull_bias.value,
                    slab_circumference_reference: slab_circumference_reference.value,
                    rate_floor_fraction: reorg_rate_floor_fraction.value,
                    supercontinent_share: supercontinent_share.value,
                    supercontinent_tenure_my: supercontinent_tenure_my.value,
                    breakup_piece_max_share: breakup_piece_max_share.value,
                    breakup_cascade_my: breakup_cascade_my.value,
                    rift_persistence_my: rift_persistence_my.value,
                    weak_zone_strength_ratio: weak_zone_strength_ratio.value,
                    plate_count_min: plate_count_min.value as u32,
                    plate_count_max: plate_count_max.value as u32,
                    plate_count_min_area_fraction: plate_count_min_area_fraction.value,
                    subduction_initiation_rate_deg_per_my: subduction_initiation_rate_deg_per_my.value,
                    plume_head_interval_my: plume_head_interval_my.value,
                    plume_tail_lifetime_my: plume_tail_lifetime_my.value,
                    plume_split_min_share: plume_split_min_share.value,
                    rate_scale_deg_per_my: params.tectonics.rate_scale_deg_per_my,
                    plume_head_radius_m: head_radius,
                };
                rp.validate().map_err(PreviewError::Domain)?;
                let schedule = reorg_schedule(
                    &seed_path,
                    rp.interval_min_my,
                    rp.interval_max_my,
                    params.tectonics.duration_my,
                );
                Some((rp, schedule))
            }
            _ => None,
        };
        let rift_flank = match &u.law_set.rift_flank {
            RiftFlankLawV1::NoneV0 => None,
            RiftFlankLawV1::ShoulderUpliftThermalV1 {
                rift_flank_flexural_parameter_m,
                rift_shoulder_ratio,
            } => Some(RiftFlankParams {
                flexural_parameter_m: rift_flank_flexural_parameter_m.value,
                shoulder_ratio: rift_shoulder_ratio.value,
            }),
        };
        let trench = match &u.law_set.subduction_trench {
            SubductionTrenchLawV1::NoneV0 => None,
            SubductionTrenchLawV1::SlabFlexureTrenchV1 {
                trench_vertical_load_fraction,
                trench_slab_length_max_m,
                trench_thermal_expansivity_per_k,
                trench_mantle_temperature_c,
                trench_elastic_isotherm_c,
                trench_thermal_diffusivity_m2_per_yr,
                trench_elastic_modulus_pa,
                trench_poisson_ratio,
            } => {
                let tp = TrenchParams {
                    vertical_load_fraction: trench_vertical_load_fraction.value,
                    slab_length_max_m: trench_slab_length_max_m.value,
                    thermal_expansivity_per_k: trench_thermal_expansivity_per_k.value,
                    mantle_temperature_c: trench_mantle_temperature_c.value,
                    elastic_isotherm_c: trench_elastic_isotherm_c.value,
                    thermal_diffusivity_m2_per_yr: trench_thermal_diffusivity_m2_per_yr.value,
                    elastic_modulus_pa: trench_elastic_modulus_pa.value,
                    poisson_ratio: trench_poisson_ratio.value,
                    mantle_density_kg_m3: rho_m,
                    water_density_kg_m3: u.constants.seawater_density_kg_m3.value,
                    gravity_m_s2: u.constants.surface_gravity_m_s2.value,
                };
                tp.validate().map_err(PreviewError::Domain)?;
                Some(tp)
            }
        };
        // SEED-FLOODBASALT (lane FB): `plume-head-lip-swell-v1`.
        let flood_basalt = match &u.law_set.flood_basalt {
            FloodBasaltLawV1::NoneV0 => None,
            FloodBasaltLawV1::PlumeHeadLipSwellV1 {
                flood_basalt_extrusive_volume_m3,
                flood_basalt_intrusive_ratio,
                flood_basalt_radius_m,
                flood_basalt_taper_fraction,
                flood_basalt_eruption_duration_my,
                swell_tail_amplitude_m,
                swell_tail_radius_m,
                swell_head_amplitude_m,
                swell_head_radius_m,
                swell_head_decay_my,
            } => {
                let fb = FloodBasaltParams {
                    lip: LipParams {
                        extrusive_volume_m3: flood_basalt_extrusive_volume_m3.value,
                        intrusive_ratio: flood_basalt_intrusive_ratio.value,
                        radius_m: flood_basalt_radius_m.value,
                        taper_fraction: flood_basalt_taper_fraction.value,
                        eruption_duration_my: flood_basalt_eruption_duration_my.value,
                    },
                    swell: SwellParams {
                        tail_amplitude_m: swell_tail_amplitude_m.value,
                        tail_radius_m: swell_tail_radius_m.value,
                        head_amplitude_m: swell_head_amplitude_m.value,
                        head_radius_m: swell_head_radius_m.value,
                        head_decay_my: swell_head_decay_my.value,
                    },
                };
                fb.validate().map_err(PreviewError::Domain)?;
                if plumes.is_none() {
                    return Err(PreviewError::Domain(
                        "plume-head-lip-swell-v1 requires fixed-plume-chains-v1".into(),
                    ));
                }
                Some(fb)
            }
        };
        let seeded_margin_subsidence =
            !matches!(u.law_set.passive_margin, PassiveMarginLawV1::NoneV0);
        // SEED-STRATIGRAPHY (lane ST): parcel-stacks-v1 needs the parcels.
        let strata = strata_params_of(&u);
        if let Some(sp) = &strata {
            sp.validate().map_err(PreviewError::Domain)?;
            if parcel_params.is_none() {
                return Err(PreviewError::Domain("parcel-stacks-v1 requires continental-parcel-carry-v1".into()));
            }
        }
        if let Some(limit) = params.orogeny.as_ref().and_then(|o| o.rate_limited_belt.as_ref())
            .and_then(|b| b.continental_parcel_carry.as_ref()).and_then(|c| c.carrier_record_limit) {
            pf_geodynamics::conserve::CarrierIdentity::validate_limit(limit).map_err(PreviewError::Domain)?;
            if parcel_params.is_none() {
                return Err(PreviewError::Domain("carrier recording requires native parcel carry".into()));
            }
        }
        if let Some(c) = params.orogeny.as_ref().and_then(|o|o.rate_limited_belt.as_ref())
            .and_then(|b|b.continental_parcel_carry.as_ref()) {
            if c.layer_provenance && (c.carrier_record_limit.is_none() || strata.is_none()) {
                return Err(PreviewError::Domain("layer provenance requires carrier IDs and strata".into()));
            }
            if let Some(h)=&c.stack_history {
                pf_geodynamics::carrier_history::StackHistory::validate_config(h).map_err(PreviewError::Domain)?;
                if c.carrier_record_limit.is_none() || strata.is_none() {
                    return Err(PreviewError::Domain("stack history requires native carrier IDs and strata".into()));
                }
            }
        }
        if let Some(config)=params.orogeny.as_ref().and_then(|o|o.sediment_retention.as_ref()).and_then(|s|s.fluvial_trace.as_ref()) {
            pf_surface::erosion::fluvial_trace::FluvialArchive::validate_config(config).map_err(PreviewError::Domain)?;
            if retention.is_none() || (parcel_params.is_none()&&suture_contact_m.is_none()) || config.steps.last().is_some_and(|s|*s>steps as u64) {
                return Err(PreviewError::Domain("fluvial trace requires retention/conservation state and selected steps within the schedule".into()));
            }
        }
        Ok(Self {
            params,
            grid_override,
            universe: u,
            grid,
            seed_path,
            prior,
            ep,
            areas,
            sea,
            rho_c,
            airy_factor,
            flags_on,
            airy_factor_oceanic,
            dt,
            steps,
            rain,
            stencil,
            belt,
            subgrid_x_c_m,
            relief_law,
            critical_slope,
            retention,
            shelf_grading,
            closed_basin,
            surface_laws,
            carry,
            collision_sediment,
            plumes,
            rift,
            split_rate_deg_per_my,
            rift_flank,
            trench,
            seeded_margin_subsidence,
            parcel_params,
            parcel_residual_iterations,
            suture_contact_m,
            flood_basalt,
            reorg,
            reorg_law,
            indent,
            plateau,
            magma,
            plateau_in_belt,
            reorg_initiation_age_my,
            magma_heat_dry_factor,
            boundary_transform,
            reorg_limits,
            reorg_margin,
            rigid_group_release,
            strata,
            climate,
            incision_threshold,
            glacial,
            ocean_layer,
        })
    }

    /// SEED-FLOODBASALT (lane FB): the `plume-head-lip-swell-v1`
    /// parameters.
    pub fn flood_basalt_params(&self) -> Option<&FloodBasaltParams> {
        self.flood_basalt.as_ref()
    }

    /// SEED-FLOODBASALT (lane FB): the mantle-fixed plume swell column term
    /// (water-loaded m) at time `t_my` (`None` unless the law is selected).
    pub fn swell_at(&self, t_my: f64) -> Result<Option<Vec<f64>>, PreviewError> {
        let (Some(fb), Some((plumes, _))) = (&self.flood_basalt, &self.plumes) else {
            return Ok(None);
        };
        let cents = cell_centroids(&self.grid).map_err(PreviewError::Domain)?;
        // SEED-WORLD-COMPOSE (lane CP): one plume model. Under the scripted
        // reorganization tails end after their lifetime, and so do their
        // swells (no profile combined the two before world_v1).
        if let Some((rp, _)) = &self.reorg {
            let live: Vec<Plume> = plumes.iter().filter(|p| p.onset_my + rp.plume_tail_lifetime_my > t_my).cloned().collect();
            return Ok(Some(swell_term(&cents, self.grid.radius_m(), &live, &fb.swell, t_my)));
        }
        Ok(Some(swell_term(&cents, self.grid.radius_m(), plumes, &fb.swell, t_my)))
    }

    /// SEED-TRENCH (lane TN): the `slab-flexure-trench-v1` parameters.
    pub fn trench_params(&self) -> Option<&TrenchParams> {
        self.trench.as_ref()
    }

    /// SEED-TRENCH (lane TN): the `slab-flexure-trench-v1` column term and
    /// its step log for a crust state (`None` unless the law is selected).
    pub fn trench_term_of(
        &self,
        state: &CrustState,
    ) -> Result<Option<(Vec<f64>, TrenchStepLog)>, PreviewError> {
        let Some(tp) = &self.trench else {
            return Ok(None);
        };
        let t = trench_term(&self.grid, &self.areas, state, tp).map_err(PreviewError::Domain)?;
        let log = TrenchStepLog {
            trench_cells: t.trench_cells.len() as u64,
            flexed_cells: t.flexed_cells,
            max_w0_m: t.max_w0_m,
            max_bulge_wl_m: t.max_bulge_wl_m,
        };
        Ok(Some((t.total_wl_m, log)))
    }

    /// SEED-EXTRAS (lane EX): the seeded plumes (empty unless
    /// `fixed-plume-chains-v1`), for inspection and maps.
    pub fn plumes_for_inspection(&self) -> &[Plume] {
        self.plumes.as_ref().map(|(p, _)| p.as_slice()).unwrap_or(&[])
    }

    /// True when a SEED-EXTRAS (lane EX) plume or rift law is selected.
    pub fn extras_active(&self) -> bool {
        self.plumes.is_some() || self.rift.is_some()
    }

    /// SEED-CARRY-CONSERVE (lane CV): true when the parcel carry or the
    /// suture law is selected (both carry `ConserveRunState`).
    pub fn conserve_active(&self) -> bool {
        self.parcel_params.is_some() || self.suture_contact_m.is_some()
    }

    /// True when any SEED-OROGEN law is selected.
    pub fn orogen_active(&self) -> bool {
        self.belt.is_some() || self.subgrid_x_c_m.is_some()
    }

    pub fn params(&self) -> &PlanetParameterProfile {
        &self.params
    }
    pub fn universe(&self) -> &UniverseV1 {
        &self.universe
    }
    pub fn grid(&self) -> &GridSpec {
        &self.grid
    }
    pub fn seed_path(&self) -> &SeedPath {
        &self.seed_path
    }
    pub fn areas(&self) -> &[f64] {
        &self.areas
    }
    /// The recorded step length (My) and the scheduled step count.
    pub fn dt_my(&self) -> f64 {
        self.dt
    }
    pub fn steps_scheduled(&self) -> u32 {
        self.steps
    }
}

/// WS-MODULE: one boundary-table upsert as `step` performed it (the
/// genesis module journals these so a restore can rebuild the table by
/// replaying them; `BoundaryTable` keeps its fields private).
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct BoundaryUpsert {
    pub t_my: f64,
    pub a: pf_contracts::domain::PlateId,
    pub b: pf_contracts::domain::PlateId,
    pub kind: BoundaryKindV0,
    pub footprint: u64,
}

/// WS-MODULE: what one `PreviewRun::step` measured, for the genesis
/// module's receipts and conservation closures. Reads only; the batch path
/// ignores it.
#[derive(Debug, Clone)]
pub struct PreviewStepReport {
    pub t_my: f64,
    pub boundary_upserts: Vec<BoundaryUpsert>,
    pub advance: AdvanceLedger,
    /// Bedrock (crust column) volume before advance, after advance and
    /// after erosion (m^3, area-weighted with the run's cell areas).
    pub bedrock_before_m3: f64,
    pub bedrock_after_advance_m3: f64,
    pub bedrock_final_m3: f64,
    /// Sediment volume at the same three points (m^3).
    pub sediment_before_m3: f64,
    pub sediment_after_advance_m3: f64,
    pub sediment_final_m3: f64,
    pub bedrock_eroded_m3: f64,
    pub sediment_eroded_m3: f64,
    pub eroded_m3: f64,
    pub deposited_m3: f64,
    pub exported_m3: f64,
    pub deposit_aggraded_m3: f64,
    /// Lane SD laws only: the named partition of this step's moved material.
    pub deposition: Option<DepositionReport>,
    /// Two-column mode: the sea-level solves after advance (the tectonic
    /// surface) and after erosion (re-equilibration).
    pub tectonic_solve: Option<SeaLevelSolve>,
    pub solve: Option<SeaLevelSolve>,
    /// SEED-CARRY-CONSERVE (lane CV): the grid's continental ledger after
    /// the step (every law) and, under the parcel carry / suture laws, the
    /// parcel ledger.
    pub continental: ContinentalLedgerV1,
    pub conserve: Option<ConserveStepLog>,
    /// SEED-STRATIGRAPHY (lane ST): `parcel-stacks-v1` only.
    pub stratigraphy: Option<StratigraphyStepLog>,
    /// SEED-PLATE-REORG (lane PR): `plate-reorg-scripted-v1` only.
    pub reorg: Option<ReorgStepReport>,
    /// SEED-RANGES (lane RG): indentation / plateau laws only.
    pub ranges: Option<RangesStepLog>,
    /// SEED-WORLD-COMPOSE (lane CP; lane VL defect-4 diagnostic): where the
    /// step's shortening volume went, by distance to a convergent boundary
    /// or suture (`rate-limited-belt-v1` only; feeds no output).
    pub shortening_budget: Option<ShorteningBudgetLog>,
}

/// SEED-WORLD-COMPOSE (lane CP): distance bins (km) of the shortening
/// budget: [0, 250), [250, 500), [500, 1000), [1000, inf), no source.
pub const SHORTENING_BUDGET_BINS_KM: [f64; 3] = [250.0, 500.0, 1000.0];

/// Per-step shortening budget (m^3 per distance bin, see
/// `SHORTENING_BUDGET_BINS_KM`): belt supply placed, net change from the
/// plateau law (v2 in-belt pass, or v1 post-advance spill), net change from
/// the channel-flow spreading, delaminated by the backstop, and bedrock
/// eroded from continental columns. Distance is to the nearest cell on a
/// subduction or collision boundary or on a suture (continental edge
/// between member plates of one suture group) at the step end.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ShorteningBudgetLog {
    pub placed_m3: [f64; 5],
    pub plateau_m3: [f64; 5],
    pub spread_m3: [f64; 5],
    pub delaminated_m3: [f64; 5],
    pub eroded_m3: [f64; 5],
}

/// WS-MODULE: the chain's evolving state between steps. `init` is stage A,
/// `step` is one iteration of the former `simulate_with` loop, `finish` is
/// the former tail (final hydrology, terrain, resource). The code is the
/// former loop body moved here unchanged; only `self.` was added.
pub struct PreviewRun {
    pub(crate) native_evidence: Option<crate::preview_module::NativeEvidenceStateV1>,
    pub(crate) env: PreviewEnv,
    pub(crate) plates: Vec<PlateRecordV0>,
    pub(crate) state: CrustState,
    pub(crate) bodies_report: Option<ContinentalBodiesReport>,
    pub(crate) elevation: Vec<f64>,
    pub(crate) elevation_t0: Vec<f64>,
    pub(crate) sea_level_t0: Option<SeaLevelSolve>,
    pub(crate) uplift_rate: Vec<f64>,
    pub(crate) sediment_ledger: SedimentLedger,
    pub(crate) boundary_table: BoundaryTable,
    pub(crate) crust: CrustVolumeLedgerV0,
    pub(crate) nu_max_observed: f64,
    pub(crate) substeps_max: u32,
    pub(crate) max_uplift: f64,
    pub(crate) max_erosion_rate: f64,
    pub(crate) logs: Vec<StepLog>,
    pub(crate) exchange_epochs: Vec<pf_contracts::domain::WaterExchangeEpochV0>,
    pub(crate) absolute_drift: f64,
    /// SEED-OROGEN state; `None` under the legacy laws.
    pub(crate) orogen: Option<OrogenRunState>,
    /// CLIMATE-GENESIS (lane CG): `epoch-climate-runoff-v1` state; `None`
    /// unless the law is selected.
    pub(crate) climate: Option<GenesisClimateState>,
}

impl PreviewRun {
    /// Stage A: plates and initial crust.
    pub fn init(env: PreviewEnv) -> Result<Self, PreviewError> {
        let grid = &env.grid;
        let u = &env.universe;
        let n = grid.cells();
        let seeding = seed_plates(grid, &env.params, &env.seed_path)?;
        let plates = seeding.plates.clone();
        // SEED-IRREGULAR-SEEDS (lane IS): per-cell margin stretch of the
        // craton seeding (None on every other path).
        let mut craton_beta: Option<Vec<f64>> = None;
        let (mut state, bodies_report) = match &u.law_set.crust_seeding {
            CrustSeedingLawV1::WholePlateContinentsV0 => {
                (CrustState::seed(grid, &seeding, &env.prior)?, None)
            }
            CrustSeedingLawV1::GeodesicDiscBodiesV0 {
                continental_area_fraction,
                rifted_margin_width_m,
            } => match &u.law_set.craton_seeding {
                CratonSeedingLawV1::NoneV0 => {
                    let (s, r) = CrustState::seed_bodies(
                        grid,
                        &seeding,
                        &env.prior,
                        continental_area_fraction.value,
                        rifted_margin_width_m.value,
                    )?;
                    (s, Some(r))
                }
                CratonSeedingLawV1::IrregularCratonsV1 {
                    craton_interior_thickness_m,
                    craton_rank_size_exponent,
                    cratons_per_unit_share,
                    basement_roughness_std_m,
                    suture_crest_thickening_m,
                    suture_half_width_m,
                    proximal_margin_stretch,
                } => {
                    let cp = pf_geodynamics::cratons::CratonParams {
                        interior_thickness_m: craton_interior_thickness_m.value,
                        rank_size_exponent: craton_rank_size_exponent.value,
                        cratons_per_unit_share: cratons_per_unit_share.value,
                        basement_roughness_std_m: basement_roughness_std_m.value,
                        suture_crest_thickening_m: suture_crest_thickening_m.value,
                        suture_half_width_m: suture_half_width_m.value,
                        proximal_margin_stretch: proximal_margin_stretch.value,
                    };
                    let (s, r, x) = pf_geodynamics::cratons::seed_cratons(
                        grid,
                        &seeding,
                        &env.prior,
                        continental_area_fraction.value,
                        rifted_margin_width_m.value,
                        &cp,
                        env.params.tectonics.continental_seeds,
                        pf_geodynamics::plates::seed_from_path_tagged(
                            &env.seed_path,
                            "irregular-cratons-v1",
                        ),
                    )?;
                    craton_beta = Some(x.stretch_beta);
                    (s, Some(r))
                }
            },
        };
        let (class, rate, _) = classify_cells_law(grid, &state, &plates, env.boundary_transform.as_ref())?;
        state.boundary_class = class;
        state.normal_rate_m_per_my = rate;
        let initial_volume = state.crust_volume_m3(&env.areas);
        // SEED-TRENCH: the trench term from the start (None = the legacy call).
        let trench0 = env.trench_term_of(&state)?;
        // SEED-FLOODBASALT: the plume swell from the start (None = the legacy call).
        let swell0 = env.swell_at(state.t_my)?;
        let trench0_sw: Option<Vec<f64>> = swell0
            .as_ref()
            .map(|s| add_column_term(trench0.as_ref().map(|t| t.0.as_slice()), s));
        let (elevation, sea_level_t0) = elevation_on_grid_terms(
            &state,
            &env.params,
            &env.areas,
            swell0.as_deref(),
            trench0_sw.as_deref().or(trench0.as_ref().map(|t| t.0.as_slice())),
        )?;
        let elevation_t0 = elevation.clone();
        let uplift_rate = vec![0.0f64; n];
        let sediment_ledger =
            SedimentLedger::new(env.rho_c, SEDIMENT_REL_TOL, SEDIMENT_ABS_FLOOR_KG)?;
        let boundary_table = BoundaryTable::new(0.0);
        let crust = CrustVolumeLedgerV0 {
            initial_m3: initial_volume,
            ridge_created_m3: 0.0,
            arc_created_m3: 0.0,
            subducted_m3: 0.0,
            delaminated_m3: 0.0,
            eroded_m3: 0.0,
            advection_residual_m3: 0.0,
            closure_plug_m3: 0.0,
            final_m3: initial_volume,
        };
        let logs = Vec::with_capacity(env.steps as usize);
        // CLIMATE-GENESIS (lane CG): empty state; the first step solves.
        let climate = env.climate.as_ref().map(|_| GenesisClimateState::default());
        let orogen = if env.orogen_active() {
            Some(OrogenRunState {
                subgrid_relief_m: env.subgrid_x_c_m.map(|_| vec![0.0f64; n]),
                steps: Vec::with_capacity(env.steps as usize),
                extras: env.extras_active().then(|| ExtrasRunState {
                    volcanic_added_m: vec![0.0; n],
                    lip_added_m: env.flood_basalt.map(|_| vec![0.0; n]),
                    volcanic_t_my: vec![-1.0; n],
                    rift_t_my: vec![-1.0; n],
                    rift_thinned_m: vec![0.0; n],
                    splits: Vec::new(),
                    rift_flank: env.rift_flank.map(|_| {
                        let mut m = RiftFlankMarkers {
                            beta: vec![1.0; n],
                            thermal_t0_my: vec![-1.0; n],
                        };
                        // SEED-TRENCH part 2 (`seeded-margin-thermal-subsidence-v1`):
                        // the seeded rifted-margin taper is stretched crust
                        // (beta = h_cc / h) that broke up when the adjacent
                        // seeded ocean was born.
                        if env.seeded_margin_subsidence {
                            let h_cc = env.prior.record.continental_crust_thickness_m;
                            let t0 = -env.prior.initial_oceanic_age_my;
                            match &craton_beta {
                                // SEED-IRREGULAR-SEEDS (lane IS): the craton
                                // margins' own stretch (h_int / h), not 35 km / h.
                                Some(b) => {
                                    for c in 0..n {
                                        if state.kind[c] == CrustKindCodeV0::Continental && b[c] > 1.0 {
                                            m.beta[c] = b[c];
                                            m.thermal_t0_my[c] = t0;
                                        }
                                    }
                                }
                                None => {
                                    for c in 0..n {
                                        let h = state.thickness_m[c];
                                        if state.kind[c] == CrustKindCodeV0::Continental && h > 0.0 && h < h_cc {
                                            m.beta[c] = h_cc / h;
                                            m.thermal_t0_my[c] = t0;
                                        }
                                    }
                                }
                            }
                        }
                        m
                    }),
                }),
                conserve: if env.conserve_active() {
                    Some(ConserveRunState {
                        fluvial_trace:env.params.orogeny.as_ref().and_then(|o|o.sediment_retention.as_ref()).and_then(|s|s.fluvial_trace.as_ref())
                            .map(|c|pf_surface::erosion::fluvial_trace::FluvialArchive::new(c.clone(),grid,env.dt)).transpose().map_err(PreviewError::Domain)?,
                        parcels: if env.parcel_params.is_some() {
                            let limit = env.params.orogeny.as_ref().and_then(|o| o.rate_limited_belt.as_ref())
                                .and_then(|b| b.continental_parcel_carry.as_ref()).and_then(|c| c.carrier_record_limit);
                            let mut l = match limit {
                                Some(limit) => ParcelLayer::seed_with_identity(grid, &state, limit)?,
                                None => ParcelLayer::seed(grid, &state)?,
                            };
                            // SEED-STRATIGRAPHY (lane ST): basement-only stacks
                            // (the column at t = 0), aged by the parcel's birth.
                            if env.strata.is_some() {
                                l.strata = Some(
                                    (0..l.len())
                                        .map(|k| {
                                            StrataStack::basement_only(
                                                l.thickness_m[k] + l.sediment_m[k],
                                                pf_geodynamics::stratigraphy::LITH_BASEMENT,
                                                l.birth_t_my[k],
                                            )
                                        })
                                        .collect(),
                                );
                            }
                            // SEED-STRATIGRAPHY stage A.1 (lane SA,
                            // parcel-stacks-v2): the t = 0 platform cover by
                            // the cell's t = 0 elevation and the rock-uplift
                            // reference (parcels are seeded one per
                            // continental cell in ascending cell order).
                            if let (Some(v), Some(stacks)) = (env.strata.as_ref().and_then(|sp| sp.v2), l.strata.as_mut()) {
                                let mut k = 0usize;
                                for c in 0..n {
                                    if state.kind[c] == CrustKindCodeV0::Continental {
                                        pf_geodynamics::stratigraphy::seed_platform_cover(&mut stacks[k], elevation[c] - env.sea, &v);
                                        stacks[k].elev_ref_m = elevation[c];
                                        k += 1;
                                    }
                                }
                                if k != stacks.len() {
                                    return Err(PreviewError::Domain(format!("parcel-stacks-v2: {k} continental cells for {} seeded parcels", stacks.len())));
                                }
                            }
                            if let Some(config)=env.params.orogeny.as_ref().and_then(|o|o.rate_limited_belt.as_ref())
                                .and_then(|b|b.continental_parcel_carry.as_ref()).and_then(|c|c.stack_history.as_ref()) {
                                let mut h=pf_geodynamics::carrier_history::StackHistory::new(config.clone(),*env.strata.as_ref().expect("validated"),grid)
                                    .map_err(PreviewError::Domain)?;
                                h.reserve_initial().map_err(PreviewError::Domain)?;
                                for (k,id) in l.identity.as_ref().expect("validated").ids().iter().enumerate() {
                                    h.initial(*id,state.t_my,l.pos[k],l.plate[k].0,l.area_m2[k],&l.strata.as_ref().expect("validated")[k])
                                        .map_err(PreviewError::Domain)?;
                                }
                                l.history=Some(h);
                            }
                            if env.params.orogeny.as_ref().and_then(|o|o.rate_limited_belt.as_ref())
                                .and_then(|b|b.continental_parcel_carry.as_ref()).is_some_and(|c|c.layer_provenance) {
                                l.layer_provenance=Some(l.strata.as_ref().expect("validated strata").iter()
                                    .zip(l.identity.as_ref().expect("validated IDs").ids())
                                    .map(|(s,id)|pf_geodynamics::layer_provenance::initial(s,*id)).collect());
                            }
                            l
                        } else {
                            ParcelLayer::default()
                        },
                        suture_parent: Vec::new(),
                        sutures: Vec::new(),
                        reorg: env.reorg.as_ref().map(|_| ReorgRunState::default()),
                        indent: env.indent.as_ref().map(|_| Vec::new()),
                        magma: env.magma.as_ref().map(|m| crate::magmatism::MagmaRunState::new_for(n, m.params.v2.is_some())),
                    })
                } else {
                    None
                },
            })
        } else {
            None
        };
        let native_evidence = crate::preview_module::NativeEvidenceStateV1::init(&env, &state)?;
        Ok(Self {
            native_evidence,
            env,
            plates,
            state,
            bodies_report,
            elevation,
            elevation_t0,
            sea_level_t0,
            uplift_rate,
            sediment_ledger,
            boundary_table,
            crust,
            nu_max_observed: 0.0,
            substeps_max: 1,
            max_uplift: 0.0,
            max_erosion_rate: 0.0,
            logs,
            exchange_epochs: Vec::new(),
            absolute_drift: 0.0,
            orogen,
            climate,
        })
    }

    pub fn env(&self) -> &PreviewEnv {
        &self.env
    }
    pub fn native_evidence(&self) -> Option<&crate::preview_module::NativeEvidenceStateV1> {
        self.native_evidence.as_ref()
    }
    pub fn state(&self) -> &CrustState {
        &self.state
    }
    pub fn elevation_m(&self) -> &[f64] {
        &self.elevation
    }
    /// Actual evolving subgrid relief, exported for existing terrain consumers.
    pub fn subgrid_relief_m(&self) -> Option<&[f64]> {
        self.orogen.as_ref().and_then(|o|o.subgrid_relief_m.as_deref())
    }
    pub fn uplift_rate_m_per_my(&self) -> &[f64] {
        &self.uplift_rate
    }
    pub fn logs(&self) -> &[StepLog] {
        &self.logs
    }
    /// CLIMATE-GENESIS (lane CG): the climate forcing state (None unless
    /// `epoch-climate-runoff-v1` is selected).
    pub fn climate_forcing(&self) -> Option<&GenesisClimateState> {
        self.climate.as_ref()
    }
    pub fn crust_ledger(&self) -> &CrustVolumeLedgerV0 {
        &self.crust
    }
    pub fn steps_taken(&self) -> u32 {
        self.logs.len() as u32
    }
    pub fn steps_scheduled(&self) -> u32 {
        self.env.steps
    }
    /// SEED-PLATE-REORG (lane PR): the plate records (poles, lineage).
    pub fn plates(&self) -> &[PlateRecordV0] {
        &self.plates
    }
    /// SEED-PLATE-REORG (lane PR): the parcel / suture / reorg run state.
    pub fn conserve_state(&self) -> Option<&ConserveRunState> {
        self.orogen.as_ref().and_then(|o| o.conserve.as_ref())
    }

    /// Recording scope is carrier identity only, not material ancestry.
    pub fn carrier_identity(&self) -> Result<Option<&pf_geodynamics::conserve::CarrierIdentity>, PreviewError> {
        let limit = self.env.params.orogeny.as_ref().and_then(|o| o.rate_limited_belt.as_ref())
            .and_then(|b| b.continental_parcel_carry.as_ref()).and_then(|c| c.carrier_record_limit);
        let layer = self.conserve_state().map(|c| &c.parcels);
        let enabled=self.env.params.orogeny.as_ref().and_then(|o|o.rate_limited_belt.as_ref())
            .and_then(|b|b.continental_parcel_carry.as_ref()).is_some_and(|c|c.layer_provenance);
        let stored=layer.and_then(|l|l.layer_provenance.as_ref());
        if enabled!=stored.is_some() || stored.is_some_and(|p|Some(p.len())!=layer.map(|l|l.len())) {
            return Err(PreviewError::Domain("layer provenance differs from profile/parcel count".into()));
        }
        match (limit, layer.and_then(|l| l.identity.as_ref())) {
            (None, None) => Ok(None),
            (Some(limit), Some(ids)) if ids.record_limit() == limit => {
                ids.check_context(self.env.grid.cells(), layer.expect("identity has layer").len())
                    .map_err(PreviewError::Domain)?;
                Ok(Some(ids))
            }
            _ => Err(PreviewError::Domain("carrier recording mode/limit differs from profile".into())),
        }
    }

    pub fn carrier_history(&self)->Result<Option<&pf_geodynamics::carrier_history::StackHistory>,PreviewError> {
        let config=self.env.params.orogeny.as_ref().and_then(|o|o.rate_limited_belt.as_ref())
            .and_then(|b|b.continental_parcel_carry.as_ref()).and_then(|c|c.stack_history.as_ref());
        let history=self.conserve_state().and_then(|c|c.parcels.history.as_ref());
        match (config,history) {
            (None,None)=>Ok(None),
            (Some(c),Some(h)) if h.config()==c && Some(h.params())==self.env.strata.as_ref()
                && h.grid()==[self.env.grid.nlat(),self.env.grid.nlon()] && h.radius_m().to_bits()==self.env.grid.radius_m().to_bits()=>{
                h.validate_schedule(self.steps_taken() as u64,self.env.dt).map_err(PreviewError::Domain)?;
                Ok(Some(h))
            },
            _=>Err(PreviewError::Domain("stack history profile mismatch".into()))
        }
    }
    pub fn validate_carrier_history(&self)->Result<(),PreviewError> {
        if let Some(h)=self.carrier_history()? {
            let p=&self.conserve_state().expect("history layer").parcels;
            h.validate_context(self.carrier_identity()?.ok_or_else(||PreviewError::Domain("missing history IDs".into()))?,
                self.state.t_my,p.strata.as_deref()).map_err(PreviewError::Domain)?;
        }
        Ok(())
    }

    pub fn fluvial_trace(&self)->Result<Option<&pf_surface::erosion::fluvial_trace::FluvialArchive>,PreviewError>{
        let config=self.env.params.orogeny.as_ref().and_then(|o|o.sediment_retention.as_ref()).and_then(|s|s.fluvial_trace.as_ref());
        let trace=self.conserve_state().and_then(|c|c.fluvial_trace.as_ref());
        match (config,trace){
            (None,None)=>Ok(None),
            (Some(c),Some(t)) if t.config()==c&&t.grid()==[self.env.grid.nlat(),self.env.grid.nlon()]&&t.dt_my().to_bits()==self.env.dt.to_bits()=>{
                t.check_grid(&self.env.grid).map_err(PreviewError::Domain)?;
                t.check_schedule(self.steps_taken() as u64).map_err(PreviewError::Domain)?;Ok(Some(t))
            },
            _=>Err(PreviewError::Domain("fluvial trace profile mismatch".into()))
        }
    }
    /// SEED-MAGMATISM (lane MG): the magmatism run state (markers, logs).
    pub fn magma_state(&self) -> Option<&crate::magmatism::MagmaRunState> {
        self.conserve_state().and_then(|c| c.magma.as_ref())
    }
    /// SEED-MAGMATISM (lane MG): the derived volcano record of the current
    /// state (None unless the law is selected).
    pub fn volcano_features(
        &self,
    ) -> Result<Option<pf_contracts::terrain::VolcanoFeaturesDomainV0>, PreviewError> {
        match (self.env.magma.as_ref(), self.magma_state()) {
            (Some(me), Some(ms)) => crate::magmatism::volcano_features(
                me,
                &self.env.grid,
                &self.state,
                &self.elevation,
                self.env.sea,
                ms,
                self.env.plumes.as_ref().map(|p| p.0.as_slice()),
            )
            .map(Some)
            .map_err(PreviewError::Domain),
            _ => Ok(None),
        }
    }
    /// SEED-STRATIGRAPHY (lane ST): the rock column record of the current
    /// state to `depth_m` (the law's publish depth when None; pass
    /// `f64::INFINITY` for the full column). None unless the law is
    /// selected.
    pub fn stratigraphy(&self, depth_m: Option<f64>) -> Result<Option<pf_contracts::terrain::StratigraphyDomainV0>, PreviewError> {
        let (Some(sp), Some(cs)) = (self.env.strata.as_ref(), self.conserve_state()) else {
            return Ok(None);
        };
        let Some(stacks) = cs.parcels.strata.as_ref() else {
            return Ok(None);
        };
        let dom = pf_geodynamics::conserve::dominant_parcels(&self.env.grid, &cs.parcels, &self.state)?;
        let mut d = pf_geodynamics::stratigraphy::publish(
            &self.env.grid,
            self.state.t_my,
            stacks,
            &dom,
            &self.state.sediment_m,
            &self.state.thickness_m,
            &self.state.birth_t_my,
            sp,
            depth_m.unwrap_or(sp.publish_depth_m),
            strata_assumptions(sp),
        );
        if sp.v2.is_some() {
            strata_uplift_into(&mut d, &self.env.grid, &cs.parcels, &self.state)?;
        }
        d.validate().map_err(|e| PreviewError::Domain(e.to_string()))?;
        Ok(Some(d))
    }
    /// SEED-STRATIGRAPHY (lane ST): the parcel stacks themselves (None
    /// unless the law is selected).
    pub fn strata_stacks(&self) -> Option<&[StrataStack]> {
        self.conserve_state().and_then(|c| c.parcels.strata.as_deref())
    }
    /// SEED-MAGMATISM (lane MG): the additive resource families of the
    /// current state (None unless the law is selected).
    pub fn magmatic_potentials(
        &self,
    ) -> Result<Option<Vec<crate::magmatism::MagmaticFamily>>, PreviewError> {
        match self.magma_state() {
            Some(ms) if self.env.magma.is_some() => crate::magmatism::magmatic_potentials(
                &self.state,
                &self.elevation,
                self.env.sea,
                &self.env.areas,
                ms,
            )
            .map(Some)
            .map_err(PreviewError::Domain),
            _ => Ok(None),
        }
    }

    /// One step of the chain at the recorded `dt` (tectonics, routing,
    /// erosion and the sediment ledger, isostatic re-equilibration).
    pub fn step(&mut self) -> Result<PreviewStepReport, PreviewError> {
        if self.steps_taken() >= self.env.steps {
            return Err(PreviewError::Domain(format!(
                "step schedule complete ({} of {} steps); refusing to step past it",
                self.steps_taken(),
                self.env.steps
            )));
        }
        let step_start_my = self.state.t_my;
        if let Some(e) = &self.native_evidence { e.reserve_step(&self.env)?; }
        self.carrier_identity()?;
        self.carrier_history()?;
        // CLIMATE-GENESIS (lane CG): 0-based index of this step (epoch schedule).
        let climate_step = self.steps_taken();
        let fluvial_step=self.steps_taken() as u64+1;
        let fluvial_clock=self.fluvial_trace()?.filter(|t|t.due(fluvial_step)).map(|_|pf_surface::erosion::fluvial_trace::TraceClock{
            step:fluvial_step,start_my:self.state.t_my,end_my:self.state.t_my+self.env.dt});
        let surface_budgets=fluvial_clock.is_some()&&self.fluvial_trace()?.is_some_and(|t|t.config().stage_budgets);
        if let Some(c) = self.orogen.as_mut().and_then(|o| o.conserve.as_mut()) {
            if let Some(t)=c.fluvial_trace.as_mut(){t.reserve(fluvial_step).map_err(PreviewError::Domain)?;}
            if let Some(h)=c.parcels.history.as_mut(){
                let ids=c.parcels.identity.as_ref().ok_or_else(||PreviewError::Domain("history requires native identities".into()))?;
                let stacks=c.parcels.strata.as_deref().ok_or_else(||PreviewError::Domain("history requires native stacks".into()))?;
                h.reserve_live_step(ids,stacks).map_err(PreviewError::Domain)?;
            }
            let active = c.parcels.len();
            // Reserve before any plume, plate or physical-state mutation.
            c.parcels.reserve_identity_events(self.env.grid.cells(), self.env.grid.cells(), active)?;
        }
        // SEED-MAGMATISM (lane MG): plumes drift in the mantle frame; every
        // plume consumer of this step sees the step-midpoint position.
        if let (Some(me), Some((pl, _))) = (self.env.magma.as_ref(), self.env.plumes.as_mut()) {
            me.drift(pl, self.state.t_my + 0.5 * self.env.dt);
        }
        // SEED-PLATE-REORG `plate-reorg-scripted-v1`: scheduled re-solves,
        // merges, subduction initiation and breakups, then plume-head splits
        // along weak zones (before the suture pass and the carry).
        let mut splits_now: Vec<PlateSplitRecord> = Vec::new();
        let mut reorg_report: Option<ReorgStepReport> = None;
        if let Some((rp, schedule)) = self.env.reorg.as_ref() {
            let (t0, t1) = (self.state.t_my, self.state.t_my + self.env.dt);
            let cs = self
                .orogen
                .as_mut()
                .and_then(|o| o.conserve.as_mut())
                .ok_or_else(|| PreviewError::Domain("plate-reorg-scripted-v1 selected but no conserve state".into()))?;
            let mut rs = cs
                .reorg
                .as_ref()
                .ok_or_else(|| PreviewError::Domain("plate-reorg-scripted-v1 selected but no reorg state".into()))?
                .to_state();
            let mut sut = cs.suture_state();
            let ctx = ReorgCtx { grid: &self.env.grid, seed_path: &self.env.seed_path, params: rp, schedule };
            // PLATE-POLES (lane PLATE-POLES): `rigid-group-contacts-v1` release
            // (detached suture-group members leave before each re-solve).
            let release = self.env.rigid_group_release;
            let rep = match self.env.reorg_initiation_age_my {
                _ if release => pf_geodynamics::reorg::reorg_pass_release(&ctx, self.env.reorg_law, self.env.reorg_initiation_age_my,
                    self.env.reorg_limits.as_ref(), self.env.reorg_margin.as_ref(), true, &mut self.state, &mut self.plates, &mut cs.parcels, &mut sut, &mut rs, t0, t1),
                // SEED-WORLD-COMPOSE (lane CP): plate-reorg-scripted-v3.
                // SEED-TRANSFORMS (lane TF): plate-reorg-scripted-v4.
                // SEED-OLD-OCEAN (lane OR): plate-reorg-scripted-v5.
                Some(age) if self.env.reorg_margin.is_some() => pf_geodynamics::reorg::reorg_pass_v5(&ctx, age, self.env.reorg_limits.as_ref().expect("v5 carries the v4 limits"), self.env.reorg_margin.as_ref().expect("guarded"), &mut self.state, &mut self.plates, &mut cs.parcels, &mut sut, &mut rs, t0, t1),
                Some(age) if self.env.reorg_limits.is_some() => pf_geodynamics::reorg::reorg_pass_v4(&ctx, age, self.env.reorg_limits.as_ref().expect("guarded"), &mut self.state, &mut self.plates, &mut cs.parcels, &mut sut, &mut rs, t0, t1),
                Some(age) => pf_geodynamics::reorg::reorg_pass_v3(&ctx, age, &mut self.state, &mut self.plates, &mut cs.parcels, &mut sut, &mut rs, t0, t1),
                None => reorg_pass_law(&ctx, self.env.reorg_law, &mut self.state, &mut self.plates, &mut cs.parcels, &mut sut, &mut rs, t0, t1),
            }
            .map_err(PreviewError::Domain)?;
            if let Some((plumes, _)) = &self.env.plumes {
                for pl in plumes {
                    if !(pl.onset_my > t0 && pl.onset_my <= t1) {
                        continue;
                    }
                    if let Some((e, n, parent, child, child_pole)) = plume_split_law(
                        &ctx,
                        self.env.reorg_law,
                        &mut self.state,
                        &mut self.plates,
                        &mut cs.parcels,
                        &mut sut,
                        &mut rs,
                        pl,
                        pl.onset_my,
                    )
                    .map_err(PreviewError::Domain)?
                    {
                        splits_now.push(PlateSplitRecord {
                            t_my: e.t_my,
                            plume: pl.id as u64,
                            parent: parent.0,
                            child: child.0,
                            cut_normal: n,
                            cells_moved: e.cells_moved,
                            child_pole_lat_deg: child_pole.lat_deg,
                            child_pole_lon_deg: child_pole.lon_deg,
                            child_pole_rate_deg_per_my: child_pole.rate_deg_per_my,
                        });
                    }
                }
            }
            cs.suture_parent = sut.parent.iter().map(|(a, b)| (a.0, b.0)).collect();
            cs.reorg = Some(ReorgRunState::from_state(&rs));
            reorg_report = Some(rep);
        }
        // SEED-EXTRAS `plume-head-continental-split-v1`: a plume head arriving
        // this step under continental crust splits the plate under it (before
        // the carry, so the halves open during this step).
        if let (Some(rate), Some((plumes, _))) = (self.env.split_rate_deg_per_my, &self.env.plumes) {
            let (t0, t1) = (self.state.t_my, self.state.t_my + self.env.dt);
            for pl in plumes {
                if !(pl.onset_my > t0 && pl.onset_my <= t1) {
                    continue;
                }
                let c0 = nearest_cell(&self.env.grid, pl.position)?;
                if self.state.kind[c0] != CrustKindCodeV0::Continental {
                    continue;
                }
                let parent = self.state.plate[c0];
                if let Some(sp) = split_plate(
                    &mut self.plates,
                    &mut self.state.plate,
                    &self.state.material_pos,
                    parent,
                    pl.position,
                    pl.split_azimuth_rad,
                    rate,
                    pl.onset_my,
                    pl.id as u64,
                )
                .map_err(PreviewError::Domain)?
                {
                    // SEED-CARRY-CONSERVE: the parcels follow the split.
                    if let Some(cs) = self.orogen.as_mut().and_then(|o| o.conserve.as_mut()) {
                        split_parcels(&mut cs.parcels, sp.parent, sp.child, sp.cut_normal);
                    }
                    splits_now.push(PlateSplitRecord {
                        t_my: sp.t_my,
                        plume: sp.plume,
                        parent: sp.parent.0,
                        child: sp.child.0,
                        cut_normal: sp.cut_normal,
                        cells_moved: sp.cells_moved,
                        child_pole_lat_deg: sp.child_pole.lat_deg,
                        child_pole_lon_deg: sp.child_pole.lon_deg,
                        child_pole_rate_deg_per_my: sp.child_pole.rate_deg_per_my,
                    });
                }
            }
        }
        // SEED-CARRY-CONSERVE `contact-length-pole-share-v1`: suture plates
        // whose converging continental contact reached the declared length
        // (the smaller group takes the larger's pole), before the carry.
        let mut sutures_now: Vec<SutureLog> = Vec::new();
        let mut max_contact_m = 0.0;
        let mut indent_report: Option<IndentStepReport> = None;
        if let Some(th) = self.env.suture_contact_m {
            let cs = self
                .orogen
                .as_mut()
                .and_then(|o| o.conserve.as_mut())
                .ok_or_else(|| PreviewError::Domain("suture law selected but no conserve state".into()))?;
            let mut sut = cs.suture_state();
            let (events, mx) = if let Some(ip) = self.env.indent.as_ref() {
                // SEED-RANGES `contact-length-indentation-taper-v1`.
                let mut ind: Vec<Indentation> = cs
                    .indent
                    .as_ref()
                    .map(|v| v.iter().map(IndentLog::to_state).collect())
                    .unwrap_or_default();
                let (ev, mx, irep) = suture_indent_pass(
                    &self.env.grid,
                    &self.state,
                    &mut self.plates,
                    &mut sut,
                    th,
                    ip,
                    &mut ind,
                    self.env.dt,
                )?;
                cs.indent = Some(ind.iter().map(IndentLog::from_state).collect());
                indent_report = Some(irep);
                (ev, mx)
            } else {
                suture_pass(&self.env.grid, &self.state, &mut self.plates, &mut sut, th)?
            };
            max_contact_m = mx;
            cs.suture_parent = sut.parent.iter().map(|(a, b)| (a.0, b.0)).collect();
            for e in events {
                let l = SutureLog {
                    t_my: e.t_my,
                    a: e.a.0,
                    b: e.b.0,
                    contact_m: e.contact_m,
                    winner_root: e.winner_root.0,
                    loser_root: e.loser_root.0,
                };
                cs.sutures.push(l.clone());
                sutures_now.push(l);
            }
        }
        // The parcel layer is taken out of the run state for the step.
        let mut parcel_layer: Option<ParcelLayer> = if self.env.parcel_params.is_some() {
            Some(std::mem::take(
                &mut self
                    .orogen
                    .as_mut()
                    .and_then(|o| o.conserve.as_mut())
                    .ok_or_else(|| PreviewError::Domain("parcel carry selected but no conserve state".into()))?
                    .parcels,
            ))
        } else {
            None
        };
        let env = &self.env;
        // SEED-PLATE-REORG: plume tails end after their lifetime.
        let live_plumes: Option<Vec<Plume>> = match (&env.reorg, &env.plumes) {
            (Some((rp, _)), Some((p, _))) => Some(active_plumes(p, rp.plume_tail_lifetime_my, self.state.t_my)),
            _ => None,
        };
        let plume_arg: Option<(&[Plume], PlumeParams)> = match (&live_plumes, &env.plumes) {
            (Some(v), Some((_, pp))) => Some((v.as_slice(), *pp)),
            (None, Some((p, pp))) => Some((p.as_slice(), *pp)),
            _ => None,
        };
        let grid = &env.grid;
        let params = &env.params;
        let areas = &env.areas;
        let prior = &env.prior;
        let ep = &env.ep;
        let plates = &self.plates;
        let (n, dt, sea, rho_c) = (grid.cells(), env.dt, env.sea, env.rho_c);
        let (airy_factor, airy_factor_oceanic) = (env.airy_factor, env.airy_factor_oceanic);
        let mut upserts = Vec::new();
        let bedrock_before_m3 = self.state.crust_volume_m3(areas);
        let sediment_before_m3 = volume_m3(&self.state.sediment_m, areas);

        // 1. Tectonics: plate motion, boundaries, thickening prior.
        // SEED-OROGEN: `env.belt` None = the legacy rule, bit for bit.
        let formation = self.native_evidence.as_mut().and_then(|e| e.formation.as_mut());
        let (mut next, mut ledger, mut thickening, source) = if formation.is_some() {
            let parcels = parcel_layer.as_mut().zip(env.parcel_params);
            let with_parcels = parcels.is_some();
            let extras = ExtrasLaws { plumes: plume_arg, rift: env.rift,
                flood_basalt: env.flood_basalt.map(|f| f.lip) };
            pf_geodynamics::kinematic::advance_with_formation_age_v1(
                grid, &self.state, plates, prior, dt, env.belt.as_ref(),
                if with_parcels { pf_geodynamics::carry::CarryLaw::ContinentalParcelsV1 } else { env.carry },
                env.collision_sediment,
                (with_parcels || env.extras_active()).then_some(&extras), parcels,
                if with_parcels && env.reorg_law != ReorgLaw::ScriptedV1 {
                    pf_geodynamics::kinematic::RidgeAccretion::SymmetricV1
                } else { pf_geodynamics::kinematic::RidgeAccretion::PreviousOccupantV0 },
                if with_parcels { env.plateau_in_belt.as_ref() } else { None },
                if with_parcels { env.boundary_transform.as_ref() } else { None }, formation,
            )?
        } else if let (Some(layer), Some(pp)) =
            (parcel_layer.as_mut(), env.parcel_params)
        {
            // SEED-CARRY-CONSERVE: `continental-parcel-carry-v1`.
            advance_with_parcels_transform(
                grid,
                &self.state,
                plates,
                prior,
                dt,
                env.belt.as_ref(),
                env.collision_sediment,
                &ExtrasLaws {
                    plumes: plume_arg,
                    rift: env.rift,
                    flood_basalt: env.flood_basalt.map(|f| f.lip),
                },
                layer,
                pp,
                // SEED-BREAKUP-CUT-V2 (lane BC, PREREG revision 1): v2 splits
                // ridge births between the flanks (boundary mid-ridge).
                if env.reorg_law != ReorgLaw::ScriptedV1 {
                    pf_geodynamics::kinematic::RidgeAccretion::SymmetricV1
                } else {
                    pf_geodynamics::kinematic::RidgeAccretion::PreviousOccupantV0
                },
                // SEED-WORLD-COMPOSE (lane CP): plateau-cap-spill-v2 in the belt.
                env.plateau_in_belt.as_ref(),
                // SEED-TRANSFORMS (lane TF): smoothed-normal-transform-v1 (None = legacy).
                env.boundary_transform.as_ref(),
            )?
        } else if env.extras_active() {
            advance_with_extras(
                grid,
                &self.state,
                plates,
                prior,
                dt,
                env.belt.as_ref(),
                env.carry,
                env.collision_sediment,
                &ExtrasLaws {
                    plumes: plume_arg,
                    rift: env.rift,
                    flood_basalt: env.flood_basalt.map(|f| f.lip),
                },
            )?
        } else {
            advance_with_laws(
                grid,
                &self.state,
                plates,
                prior,
                dt,
                env.belt.as_ref(),
                env.carry,
                env.collision_sediment,
            )?
        };
        // SEED-RANGES `plateau-cap-spill-v1`: continental crust above the cap
        // spills into the same plate's hinterland (after the belt step,
        // before every consumer of the advanced state); the delaminated
        // remainder is a named advance ledger term.
        let mut post_spill_dh: Option<Vec<f64>> = None;
        let plateau_report: Option<PlateauReport> = match env.plateau.as_ref() {
            Some(pp) => {
                let h0 = next.thickness_m.clone();
                let r = plateau_spill(grid, &mut next, &mut thickening, pp)?;
                post_spill_dh = Some((0..n).map(|c| next.thickness_m[c] - h0[c]).collect());
                let d = r.delaminated_m3 + r.oceanic_delaminated_m3;
                ledger.delaminated_m3 += d;
                ledger.belt_delaminated_m3 += d;
                ledger.volume_after_m3 = next.crust_volume_m3(areas);
                ledger.max_thickness_m = next.thickness_m.iter().cloned().fold(0.0, f64::max);
                Some(r)
            }
            // SEED-WORLD-COMPOSE (lane CP): the in-belt v2 pass, logged in
            // the same record.
            None => ledger.plateau_belt.as_ref().map(|b| PlateauReport {
                source_cells: b.source_cells,
                excess_m3: b.excess_m3,
                moved_m3: b.moved_m3,
                filled_cells: b.filled_cells,
                delaminated_m3: b.delaminated_m3,
                oceanic_delaminated_m3: b.oceanic_delaminated_m3,
                max_distance_m: b.max_distance_m,
                cells_at_cap: b.cells_at_cap,
            }),
        };
        // SEED-MAGMATISM (lane MG, `slab-arc-backarc-llsvp-v1`): arcs at slab
        // depth, island-arc maturation, back-arc extension and vent markers on
        // the advanced state (before every consumer of it); the arc heat is a
        // column term of both isostasy evaluations of this step.
        // SEED-STRATIGRAPHY (lane ST): the back-arc part of the magma step's
        // thickness change (total minus arc) is a tagged delta.
        let strata_pre_magma: Option<(Vec<f64>, Vec<f64>)> = match (env.strata.as_ref(), env.magma.as_ref()) {
            (Some(_), Some(_)) => Some((next.thickness_m.clone(), next.arc_added_m.clone())),
            _ => None,
        };
        let magma_term: Option<Vec<f64>> = match env.magma.as_ref() {
            Some(me) => {
                let ms = self
                    .orogen
                    .as_mut()
                    .and_then(|o| o.conserve.as_mut())
                    .and_then(|c| c.magma.as_mut())
                    .ok_or_else(|| PreviewError::Domain("magmatism selected but no magma state".into()))?;
                Some(
                    crate::magmatism::magma_step_with_formation_age_v1(
                        me,
                        grid,
                        areas,
                        plates,
                        &mut next,
                        &mut ledger,
                        &mut thickening,
                        &source,
                        ms,
                        plume_arg.as_ref().map(|p| p.0),
                        self.state.t_my,
                        dt,
                        self.native_evidence.as_mut().and_then(|e| e.formation.as_mut()),
                    )
                    .map_err(PreviewError::Domain)?,
                )
            }
            None => None,
        };
        let strata_backarc_dh: Option<Vec<f64>> = strata_pre_magma.map(|(h0, a0)| {
            (0..n).map(|c| (next.thickness_m[c] - h0[c]) - (next.arc_added_m[c] - a0[c])).collect()
        });
        // SEED-EXTRAS markers ride with the crust column (0 / never at births),
        // then this step's plume additions, rift births and flank thinning.
        let extras_carried: Option<ExtrasRunState> = self
            .orogen
            .as_ref()
            .and_then(|o| o.extras.as_ref())
            .map(|x| {
                let pick = |v: &[f64], none: f64| -> Vec<f64> {
                    source
                        .iter()
                        .map(|&s| if s == NO_SOURCE { none } else { v[s as usize] })
                        .collect()
                };
                let mut y = ExtrasRunState {
                    volcanic_added_m: pick(&x.volcanic_added_m, 0.0),
                    volcanic_t_my: pick(&x.volcanic_t_my, -1.0),
                    rift_t_my: pick(&x.rift_t_my, -1.0),
                    rift_thinned_m: pick(&x.rift_thinned_m, 0.0),
                    splits: x.splits.clone(),
                    rift_flank: x.rift_flank.as_ref().map(|m| RiftFlankMarkers {
                        beta: pick(&m.beta, 1.0),
                        thermal_t0_my: pick(&m.thermal_t0_my, -1.0),
                    }),
                    lip_added_m: x.lip_added_m.as_ref().map(|v| pick(v, 0.0)),
                };
                for &(c, dh, t) in &ledger.hotspot_cells {
                    y.volcanic_added_m[c as usize] += dh;
                    y.volcanic_t_my[c as usize] = y.volcanic_t_my[c as usize].max(t);
                }
                // SEED-FLOODBASALT: LIP crust is plume crust too.
                for &(c, dh, t) in &ledger.lip_cells {
                    y.volcanic_added_m[c as usize] += dh;
                    y.volcanic_t_my[c as usize] = y.volcanic_t_my[c as usize].max(t);
                    if let Some(l) = y.lip_added_m.as_mut() {
                        l[c as usize] += dh;
                    }
                }
                for &b in &ledger.rift_births {
                    y.rift_t_my[b as usize] = next.t_my;
                }
                for &(c, dh) in &ledger.rift_thinning {
                    y.rift_thinned_m[c as usize] += dh;
                }
                y.splits.extend(splits_now.iter().cloned());
                y
            });
        // SEED-RIFT-SHOULDERS (`shoulder-uplift-thermal-v1`): update the
        // stretching markers with this step's rift births and flank
        // thinning, then evaluate the rift thermal + shoulder column term
        // (used by both isostasy evaluations of this step).
        let mut extras_carried = extras_carried;
        let mut flank: Option<(Vec<f64>, RiftFlankStepLog)> = None;
        if let (Some(fp), Some(m)) = (
            env.rift_flank.as_ref(),
            extras_carried.as_mut().and_then(|x| x.rift_flank.as_mut()),
        ) {
            stretch_update(
                &mut m.beta,
                &mut m.thermal_t0_my,
                &ledger.rift_births,
                &ledger.rift_birth_flanks,
                &ledger.rift_thinning,
                &next.thickness_m,
                next.t_my,
                &gdh1_thermal_elevation_m,
            )
            .map_err(PreviewError::Domain)?;
            let term = flank_term(
                grid,
                areas,
                &next.kind,
                &m.beta,
                &m.thermal_t0_my,
                next.t_my,
                fp,
                &gdh1_thermal_elevation_m,
            )
            .map_err(PreviewError::Domain)?;
            let log = RiftFlankStepLog {
                stretched_cells: term.stretched_cells,
                shoulder_cells: term.shoulder_cells,
                max_thermal_wl_m: term.thermal_wl_m.iter().fold(0.0, |a: f64, &b| a.max(b)),
                max_shoulder_wl_m: term.shoulder_wl_m.iter().fold(0.0, |a: f64, &b| a.max(b)),
            };
            flank = Some((term.total_wl_m, log));
        }
        let flank_wl: Option<&[f64]> = flank.as_ref().map(|f| f.0.as_slice());
        // Subgrid relief rides with the crust column (0 at ridge births).
        let relief_carried: Option<Vec<f64>> = self
            .orogen
            .as_ref()
            .and_then(|o| o.subgrid_relief_m.as_ref())
            .map(|r| {
                source
                    .iter()
                    .map(|&s| if s == NO_SOURCE { 0.0 } else { r[s as usize] })
                    .collect()
            });
        self.crust.ridge_created_m3 += ledger.ridge_created_m3;
        self.crust.arc_created_m3 += ledger.arc_created_m3;
        // SEED-EXTRAS (fixed-plume-chains-v1): plume crust is magmatic
        // addition like arc crust; CrustVolumeLedgerV0 has no plume field, so
        // it is booked here (named per step in the module control volume).
        // Exactly 0 under every other law (and then not added at all).
        if ledger.hotspot_created_m3 != 0.0 {
            self.crust.arc_created_m3 += ledger.hotspot_created_m3;
        }
        self.crust.subducted_m3 += ledger.subducted_m3;
        self.crust.delaminated_m3 += ledger.delaminated_m3;
        // SEDIMENT-TOWERS (accreted-to-belt-crust-v1): collision sediment
        // turned into crust is the reverse of the erosion move (crust ->
        // sediment), booked as NET erosion (the contract field is "net crust
        // removed by erosion"); exactly 0 under every other law.
        self.crust.eroded_m3 -= ledger.sediment_accreted_m3;
        self.crust.advection_residual_m3 += ledger.residual_m3();
        let t_new = next.t_my;
        // Boundary table: majority kind per plate pair this step.
        let (_, _, footprint) = classify_cells_law(grid, &next, plates, env.boundary_transform.as_ref())?;
        let mut per_pair: BTreeMap<
            (pf_contracts::domain::PlateId, pf_contracts::domain::PlateId),
            (BoundaryKindV0, u64, u64),
        > = BTreeMap::new();
        for (&(a, b, kind), &count) in &footprint {
            let e = per_pair.entry((a, b)).or_insert((kind, 0, 0));
            e.2 += count;
            if count > e.1 || (count == e.1 && kind_rank(kind) < kind_rank(e.0)) {
                e.0 = kind;
                e.1 = count;
            }
        }
        for ((a, b), (kind, _, total)) in per_pair {
            upserts.push(BoundaryUpsert {
                t_my: t_new,
                a,
                b,
                kind,
                footprint: total,
            });
            self.boundary_table
                .upsert(t_new, a, b, kind, total)
                .map_err(|e| PreviewError::Domain(format!("boundary table: {e}")))?;
        }
        // SEED-TRENCH (`slab-flexure-trench-v1`): the trench and outer-rise
        // term of the advanced state (plates, kinds, ages and boundary
        // classes do not change in the surface stage, so both isostasy
        // evaluations of this step use it).
        let trench = env.trench_term_of(&next)?;
        let trench_wl: Option<&[f64]> = trench.as_ref().map(|t| t.0.as_slice());
        // SEED-FLOODBASALT (`plume-head-lip-swell-v1`): the mantle-fixed plume
        // swell at the step end on both column terms (both isostasy
        // evaluations of this step); None = the legacy terms.
        let swell = env.swell_at(t_new)?;
        let flank_sw: Option<Vec<f64>> = swell.as_ref().map(|s| add_column_term(flank_wl, s));
        let trench_sw: Option<Vec<f64>> = swell.as_ref().map(|s| add_column_term(trench_wl, s));
        let flank_wl: Option<&[f64]> = flank_sw.as_deref().or(flank_wl);
        let trench_wl: Option<&[f64]> = trench_sw.as_deref().or(trench_wl);
        // SEED-MAGMATISM (lane MG): the arc heat on both column terms.
        // SEED-WORLD-COMPOSE (lane CP) v2: not in the column terms; it is a
        // dynamic term added after each sea-level solve (displaces no water).
        let heat_dynamic: Option<(&Vec<f64>, f64)> = match (magma_term.as_ref(), env.magma_heat_dry_factor) {
            (Some(t), Some(f)) => Some((t, f)),
            _ => None,
        };
        let magma_col: Option<&Vec<f64>> = if heat_dynamic.is_some() { None } else { magma_term.as_ref() };
        let flank_mg: Option<Vec<f64>> = magma_col.map(|s| add_column_term(flank_wl, s));
        let trench_mg: Option<Vec<f64>> = magma_col.map(|s| add_column_term(trench_wl, s));
        let flank_wl: Option<&[f64]> = flank_mg.as_deref().or(flank_wl);
        let trench_wl: Option<&[f64]> = trench_mg.as_deref().or(trench_wl);
        let (mut elev_tect, tectonic_solve) =
            elevation_on_grid_terms(&next, params, areas, flank_wl, trench_wl)?;
        if let Some((t, f)) = heat_dynamic {
            add_dynamic_term(&mut elev_tect, t, sea, f);
        }
        for c in 0..n {
            self.uplift_rate[c] = if next.birth_t_my[c] == t_new {
                (elev_tect[c] - self.elevation[c]) / dt
            } else {
                let f = match airy_factor_oceanic {
                    Some(fo) if next.kind[c] == CrustKindCodeV0::Oceanic => fo,
                    _ => airy_factor,
                };
                f * thickening[c] / dt
            };
            self.max_uplift = self.max_uplift.max(self.uplift_rate[c]);
        }
        self.state = next;
        let bedrock_after_advance_m3 = self.state.crust_volume_m3(areas);
        let sediment_after_advance_m3 = volume_m3(&self.state.sediment_m, areas);

        // 2. Routing on the tectonic surface.
        let ocean = pf_surface::water::largest_ocean(grid, &elev_tect, sea)?;
        let routing = route(grid, &elev_tect, &env.rain, &ocean, sea)?;
        // CLIMATE-GENESIS (lane CG, `epoch-climate-runoff-v1`): every fluvial
        // kernel below reads the runoff-weighted area A_eff = sum a R / r_ref
        // (flow graph, fill, basins, water bodies, discharge and the
        // stratigraphy inputs stay those of `routing`). None = legacy.
        // I11 (lane CLIMATE-GENESIS-2): with `cell-runoff-subgrid-erosivity-v1`
        // the per-cell ratio R / r_ref also reaches the subgrid cell-mean
        // closure (None = v1 default, the reference runoff).
        let mut subgrid_runoff_ratio: Option<Vec<f64>> = None;
        // INLAND-BASINS: the climate fields the closed-basin balance reads
        // (runoff m/yr, PET and precipitation mm/yr); None without the law.
        let mut lake_climate: Option<(Vec<f64>, Vec<f64>, Vec<f64>)> = None;
        let forced_routing: Option<Routing> = match (&env.climate, self.climate.as_mut()) {
            (Some(cp), Some(cs)) => {
                let runoff = genesis_runoff_step(
                    cp,
                    cs,
                    grid,
                    climate_step,
                    t_new,
                    &elev_tect,
                    &routing.fill_level_m,
                    &ocean,
                    sea,
                    areas,
                )
                .map_err(|e| PreviewError::Domain(format!("climate forcing at t = {t_new} My: {e}")))?;
                if env.closed_basin.is_some() {
                    lake_climate = Some((runoff.clone(), cs.last_pet_mm_per_yr.clone(), cs.last_precipitation_mm_per_yr.clone()));
                }
                // HYPSO-LEVERS (opt-in `ErosiveRunoffLawV1`): the fluvial
                // kernels read the eroding runoff R_e; the lake balance above
                // keeps the modelled runoff. None = v1 (R_e = R, same bytes).
                let erosive = cp.erosive_runoff_m_per_yr(&runoff);
                let r_e: &[f64] = erosive.as_deref().unwrap_or(&runoff);
                let a_eff = pf_surface::routing::runoff_weighted_area_m2(grid, &routing, r_e, cp.r_ref_m_per_yr)?;
                if cp.subgrid_cell_mean == ClimateSubgridCellMeanV1::CellRunoffSubgridErosivityV1 {
                    subgrid_runoff_ratio = Some(r_e.iter().map(|r| r / cp.r_ref_m_per_yr).collect());
                }
                let mut r = routing.clone();
                r.contributing_area_m2 = a_eff;
                Some(r)
            }
            (None, None) => None,
            _ => return Err(PreviewError::Domain("climate forcing law and run state disagree".into())),
        };
        let fluvial_routing: &Routing = forced_routing.as_ref().unwrap_or(&routing);

        // 3. Erosion and the sediment ledger.
        // Lane E (SEED-EROSION-BUG, 218bad8) per-step available-material
        // bound: bedrock removal is limited to the crust column present at
        // the start of erosion. Output-identical when the bound never binds.
        // SEED-OROGEN `stream-power-subgrid-v1`: trunk incision and the
        // cell-mean relaxation target; the erosion stage realises the target
        // through per-cell effective contributing areas (orogen.rs).
        let subgrid = match (env.subgrid_x_c_m, &relief_carried) {
            (Some(x_c), Some(relief)) => {
                let wet: Vec<bool> = (0..n)
                    .map(|i| ocean[i] || routing.fill_level_m[i] > elev_tect[i] + 1.0e-9)
                    .collect();
                let flow_to: Vec<Option<usize>> =
                    routing.flow_to.iter().map(|t| t.map(|e| e.0 as usize)).collect();
                let mut reach = vec![1.0f64; n];
                for i in 0..n {
                    if let Some(j) = flow_to[i] {
                        let (ii, ij) = grid.ij(i as u64).map_err(|e| PreviewError::Grid(e.to_string()))?;
                        let (ji, jj) = grid.ij(j as u64).map_err(|e| PreviewError::Grid(e.to_string()))?;
                        reach[i] = grid
                            .center_distance_m(ii, ij, ji, jj)
                            .map_err(|e| PreviewError::Grid(e.to_string()))?
                            .max(1.0);
                    }
                }
                let cell_size: Vec<f64> = areas.iter().map(|a| a.sqrt()).collect();
                // SEED-UPLANDS (lane UP): the threshold floor on the
                // cell-mean target; with OG's n = 1 law the relief.rs
                // closure at n = 1 is used (bit-identical to OG's, lane RL
                // tests/relief_law.rs).
                let floor = match &env.critical_slope {
                    Some((cs, adj)) => Some(
                        pf_surface::uplands::threshold_lowering_isostatic(
                            adj,
                            &elev_tect,
                            &wet,
                            &routing.fill_level_m,
                            &cell_size,
                            cs,
                            env.airy_factor,
                            dt,
                        )
                        .map_err(PreviewError::Domain)?,
                    ),
                    None => None,
                };
                let linear_rp = NonlinearSubgridParams {
                    k_per_my: ep.k_per_my,
                    slope_exponent_n: 1.0,
                    reference_steepness_m: 1.0,
                    channel_head_length_m: x_c,
                };
                // I11: lane OG's n = 1 path has no runoff input, so a forced
                // cell mean runs the relief.rs closure at n = 1 (bit-identical
                // to OG's at ratio 1).
                let rp_used = env
                    .relief_law
                    .as_ref()
                    .or(floor.as_ref().map(|_| &linear_rp))
                    .or(subgrid_runoff_ratio.as_ref().map(|_| &linear_rp));
                // GLACIAL-T0 (lane GLACIAL-T0, `glacial-stream-power-blend-v1`):
                // per-cell erodibility factors from the genesis climate's
                // temperature and this step's precipitation (None = off).
                let glacial_k = match (&env.glacial, &env.climate, self.climate.as_ref()) {
                    (Some(gp), Some(cp), Some(cs)) => {
                        let column = cs.column.as_ref().ok_or_else(|| {
                            PreviewError::Domain(format!("glacial erosion at t = {t_new} My: no climate column"))
                        })?;
                        let (t_sl, half, lapse) =
                            pf_surface::climate::seed::genesis::glacial_temperature_inputs(cp, column, grid, &ocean)
                                .map_err(|e| PreviewError::Domain(format!("glacial erosion at t = {t_new} My: {e}")))?;
                        let rp = env.relief_law.as_ref().ok_or_else(|| {
                            PreviewError::Domain("glacial erosion needs the nonlinear subgrid law".into())
                        })?;
                        Some(
                            pf_surface::erosion::glacial::glacial_erosivity(
                                gp,
                                rp.m(),
                                &pf_surface::erosion::glacial::GlacialInputs {
                                    receiver_first: &routing.receiver_first,
                                    flow_to: &flow_to,
                                    ocean: &ocean,
                                    wet: &wet,
                                    fill_level_m: &routing.fill_level_m,
                                    area_m2: areas,
                                    z_mean_m: &elev_tect,
                                    relief_m: relief,
                                    sea_level_temperature_c: &t_sl,
                                    half_range_k: &half,
                                    precipitation_mm_per_yr: &cs.last_precipitation_mm_per_yr,
                                    lapse_k_per_m: lapse,
                                    datum_m: sea,
                                },
                            )
                            .map_err(PreviewError::Domain)?,
                        )
                    }
                    (None, _, _) => None,
                    _ => return Err(PreviewError::Domain("glacial erosion law without climate state".into())),
                };
                let zero_threshold = pf_surface::erosion::incision_threshold::IncisionThresholdParams {
                    critical_channel_steepness_m: 0.0,
                    trunk_only: false,
                    alluvial_grade_slope: None,
                };
                let sg = if let (Some(rp), Some(gk)) = (env.relief_law.as_ref(), glacial_k.as_ref()) {
                    // GLACIAL-T0: the threshold step carries the glacial
                    // factors; without the incision threshold it runs at
                    // k_c = 0 (relief.rs's n != 1 law bit for bit, lane IT
                    // test). Validation refuses I11 with this law.
                    let r = pf_surface::erosion::incision_threshold::subgrid_step_threshold_glacial(
                        &pf_surface::relief::SubgridInputs {
                            receiver_first: &routing.receiver_first,
                            flow_to: &flow_to,
                            contributing_area_m2: &fluvial_routing.contributing_area_m2,
                            reach_length_m: &reach,
                            wet: &wet,
                            fill_level_m: &routing.fill_level_m,
                            cell_size_m: &cell_size,
                        },
                        &elev_tect,
                        relief,
                        rp,
                        env.incision_threshold.as_ref().unwrap_or(&zero_threshold),
                        ep.k_per_my,
                        ep.m,
                        dt,
                        floor.as_ref().map(|f| f.lowering_m.as_slice()),
                        Some(gk),
                    )
                    .map_err(PreviewError::Domain)?;
                    pf_geodynamics::orogen::SubgridStep {
                        trunk_m: r.trunk_m,
                        target_lowering_m: r.target_lowering_m,
                        effective_area_m2: r.effective_area_m2,
                        capped_cells: r.capped_cells,
                        trunk_raised_cells: r.trunk_raised_cells,
                    }
                } else if let (Some(rp), Some(thr)) = (env.relief_law.as_ref(), env.incision_threshold.as_ref()) {
                    // INCISION-THRESHOLD (lane IT): the nonlinear subgrid law
                    // with the excess-steepness threshold on trunk and mean
                    // (validation refuses it together with I11).
                    let r = pf_surface::erosion::incision_threshold::subgrid_step_threshold(
                        &pf_surface::relief::SubgridInputs {
                            receiver_first: &routing.receiver_first,
                            flow_to: &flow_to,
                            contributing_area_m2: &fluvial_routing.contributing_area_m2,
                            reach_length_m: &reach,
                            wet: &wet,
                            fill_level_m: &routing.fill_level_m,
                            cell_size_m: &cell_size,
                        },
                        &elev_tect,
                        relief,
                        rp,
                        thr,
                        ep.k_per_my,
                        ep.m,
                        dt,
                        floor.as_ref().map(|f| f.lowering_m.as_slice()),
                    )
                    .map_err(PreviewError::Domain)?;
                    pf_geodynamics::orogen::SubgridStep {
                        trunk_m: r.trunk_m,
                        target_lowering_m: r.target_lowering_m,
                        effective_area_m2: r.effective_area_m2,
                        capped_cells: r.capped_cells,
                        trunk_raised_cells: r.trunk_raised_cells,
                    }
                } else if let Some(rp) = rp_used {
                    // SEED-RELIEF-LAW (lane RL): the n-generalised closure;
                    // same inputs, same result type, same realisation.
                    let r = pf_surface::relief::subgrid_step_with_floor_and_runoff(
                        &pf_surface::relief::SubgridInputs {
                            receiver_first: &routing.receiver_first,
                            flow_to: &flow_to,
                            contributing_area_m2: &fluvial_routing.contributing_area_m2,
                            reach_length_m: &reach,
                            wet: &wet,
                            fill_level_m: &routing.fill_level_m,
                            cell_size_m: &cell_size,
                        },
                        &elev_tect,
                        relief,
                        rp,
                        ep.k_per_my,
                        ep.m,
                        dt,
                        floor.as_ref().map(|f| f.lowering_m.as_slice()),
                        subgrid_runoff_ratio.as_deref(),
                    )
                    .map_err(PreviewError::Domain)?;
                    pf_geodynamics::orogen::SubgridStep {
                        trunk_m: r.trunk_m,
                        target_lowering_m: r.target_lowering_m,
                        effective_area_m2: r.effective_area_m2,
                        capped_cells: r.capped_cells,
                        trunk_raised_cells: r.trunk_raised_cells,
                    }
                } else {
                    subgrid_step(
                    &SubgridInputs {
                        receiver_first: &routing.receiver_first,
                        flow_to: &flow_to,
                        contributing_area_m2: &fluvial_routing.contributing_area_m2,
                        reach_length_m: &reach,
                        wet: &wet,
                        fill_level_m: &routing.fill_level_m,
                        cell_size_m: &cell_size,
                    },
                    &elev_tect,
                    relief,
                    ep.k_per_my,
                    ep.m,
                    x_c,
                    dt,
                )
                .map_err(PreviewError::Domain)?
                };
                let mut eff = fluvial_routing.clone();
                for i in 0..n {
                    if !wet[i] && flow_to[i].is_some() {
                        eff.contributing_area_m2[i] = sg.effective_area_m2[i];
                    }
                }
                Some((sg, eff, wet))
            }
            _ => None,
        };
        let erosion_routing = subgrid.as_ref().map(|s| &s.1).unwrap_or(fluvial_routing);
        // Lane SD: `erode_deposit_step` with the legacy surface laws IS
        // `erode_step_bounded` (bit for bit, pf-surface tests/deposition_laws.rs).
        let step = match &env.retention {
            // SEED-SEDIMENT-RETAIN (lane SR): `continental-retention-base-level-v1`
            // replaces the legacy step (validated: legacy land laws only).
            Some(rp) => {
                let continental: Vec<bool> = self
                    .state
                    .kind
                    .iter()
                    .map(|k| *k == CrustKindCodeV0::Continental)
                    .collect();
                // INLAND-BASINS (`closed-basin-water-balance-v1`): lake
                // water balance on the step's routing, then the lake-aware
                // retention step; None = the retention step exactly.
                let lake_balance = match &env.closed_basin {
                    Some(cb) => {
                        let uniform: Vec<f64> = env.rain.iter().map(|r| r / 1.0e6).collect();
                        let (r_in, e_net) = match &lake_climate {
                            Some((r, pet, p)) => pf_surface::closed_basin::balance_inputs(n, 0.0, Some((r, pet, p)), cb),
                            None => {
                                let e: Vec<f64> = uniform.iter().map(|r| cb.uniform_evaporation_to_runoff_ratio * r).collect();
                                Ok((uniform.clone(), e))
                            }
                        }
                        .map_err(PreviewError::Domain)?;
                        let lb = pf_surface::closed_basin::lake_balance(&routing, &elev_tect, areas, &r_in, &e_net)
                            .map_err(PreviewError::Domain)?;
                        pf_surface::closed_basin::record_observation(pf_surface::closed_basin::LakeObservation {
                            water_level_m: pf_surface::closed_basin::cell_water_level_m(&routing, &lb, &elev_tect, sea),
                            full_lakes: lb.full_lakes,
                            terminal_lakes: lb.terminal_lakes,
                            wet_area_m2: lb.wet_area_m2,
                            exposed_area_m2: lb.exposed_area_m2,
                        });
                        Some((lb, *cb))
                    }
                    None => None,
                };
                let fines = env.universe.law_set.marine_fines.as_ref().map(|law| {
                    let (fraction, min_depth_m) = law.parameters();
                    pf_surface::erosion::marine_fines::MarineFinesParams { fraction, min_depth_m }
                });
                let (step, trace) = if let Some(fines) = fines.as_ref() {
                    let (step, _, trace, report) = pf_surface::erosion::retain::erode_retain_step_with_marine_fines(
                    grid,
                    &env.stencil,
                    erosion_routing,
                    &elev_tect,
                    &self.state.sediment_m,
                    &ocean,
                    sea,
                    dt,
                    ep,
                    Some(&self.state.thickness_m),
                    &continental,
                    rp,
                    fluvial_clock,
                    surface_budgets,
                    lake_balance.as_ref().map(|(lb, cb)| (lb, cb)),
                    Some(fines),
                    )?;
                    self.native_evidence.as_mut().ok_or_else(|| PreviewError::Domain("missing fines evidence state".into()))?
                        .record_fines(self.logs.len() as u64 + 1, step_start_my, self.state.t_my, report)?;
                    (step, trace)
                } else {
                    let (step, _, trace) = erode_retain_step_with_lakes(
                    grid,
                    &env.stencil,
                    erosion_routing,
                    &elev_tect,
                    &self.state.sediment_m,
                    &ocean,
                    sea,
                    dt,
                    ep,
                    Some(&self.state.thickness_m),
                    &continental,
                    rp,
                    fluvial_clock,
                    surface_budgets,
                    lake_balance.as_ref().map(|(lb, cb)| (lb, cb)),
                    )?;
                    (step, trace)
                };
                if let Some(frame)=trace {
                    self.orogen.as_mut().and_then(|o|o.conserve.as_mut()).and_then(|c|c.fluvial_trace.as_mut())
                        .ok_or_else(||PreviewError::Domain("missing requested fluvial recorder".into()))?
                        .record(frame).map_err(PreviewError::Domain)?;
                }
                step
            }
            // INCISION-THRESHOLD part 1 (lane IT): the capped once-per-step
            // marine layer on the legacy step (validated: legacy laws only).
            None if env.ocean_layer != pf_surface::erosion::OceanRemainderLayer::PerMouthV0
                && env.surface_laws.is_legacy() =>
            {
                pf_surface::erosion::erode_step_bounded_layer(
                    grid,
                    &env.stencil,
                    erosion_routing,
                    &elev_tect,
                    &self.state.sediment_m,
                    &ocean,
                    sea,
                    dt,
                    ep,
                    Some(&self.state.thickness_m),
                    env.ocean_layer,
                )?
            }
            None => erode_deposit_step(
                grid,
                &env.stencil,
                erosion_routing,
                &elev_tect,
                &self.state.sediment_m,
                &ocean,
                sea,
                dt,
                ep,
                &env.surface_laws,
                Some(&self.state.thickness_m),
            )?,
        };
        // SHELF (lane SHELF): wave-base grading folded into the step as an
        // internal sediment transfer (eroded/deposited volumes unchanged).
        let step = match &env.shelf_grading {
            Some(g) => {
                let mut step = step;
                let avail: Vec<f64> = (0..n)
                    .map(|c| self.state.sediment_m[c] + step.sediment_added_m[c] - step.sediment_removed_m[c])
                    .collect();
                let surf: Vec<f64> = (0..n).map(|c| elev_tect[c] + step.dz_m[c]).collect();
                pf_geodynamics::shelf_grade::grade_to_wave_base(
                    grid,
                    areas,
                    &surf,
                    &ocean,
                    &avail,
                    sea,
                    g,
                    &mut step.sediment_removed_m,
                    &mut step.sediment_added_m,
                    &mut step.dz_m,
                )
                .map_err(PreviewError::Domain)?;
                step
            }
            None => step,
        };
        self.nu_max_observed = self.nu_max_observed.max(step.diffusion_number_max);
        self.substeps_max = self.substeps_max.max(step.diffusion_substeps);
        self.max_erosion_rate = self.max_erosion_rate.max(step.max_lowering_rate_m_per_my);
        let state = &mut self.state;
        let mut eroded_rock_m3 = 0.0;
        for c in 0..n {
            state.thickness_m[c] -= step.bedrock_removed_m[c];
            state.eroded_total_m[c] += step.bedrock_removed_m[c];
            state.sediment_m[c] += step.sediment_added_m[c] - step.sediment_removed_m[c];
            if state.sediment_m[c] < 0.0 {
                // Rounding below zero is a bug signal, not something to clamp silently.
                if state.sediment_m[c] < -1.0e-9 {
                    return Err(PreviewError::Domain(format!(
                        "sediment column negative at cell {c}: {}",
                        state.sediment_m[c]
                    )));
                }
                state.sediment_m[c] = 0.0;
            }
            eroded_rock_m3 += step.bedrock_removed_m[c] * areas[c];
        }
        self.crust.eroded_m3 += eroded_rock_m3;
        // SEED-WORLD-COMPOSE (lane CP): the shortening budget diagnostic.
        let shortening_budget: Option<ShorteningBudgetLog> = match ledger.belt_budget.as_ref() {
            Some(bb) => {
                let parent: BTreeMap<u32, u32> = self
                    .orogen
                    .as_ref()
                    .and_then(|o| o.conserve.as_ref())
                    .map(|c| c.suture_parent.iter().cloned().collect())
                    .unwrap_or_default();
                let bins = shortening_distance_bins(grid, &self.state, &parent)?;
                let mut log = ShorteningBudgetLog::default();
                for c in 0..n {
                    let b = bins[c] as usize;
                    let a = areas[c];
                    log.placed_m3[b] += bb.placed_m[c] * a;
                    log.plateau_m3[b] += bb.plateau_dh_m[c] * a + post_spill_dh.as_ref().map_or(0.0, |d| d[c]) * a;
                    log.spread_m3[b] += bb.spread_dh_m[c] * a;
                    log.delaminated_m3[b] += bb.delaminated_m[c] * a;
                    if self.state.kind[c] == CrustKindCodeV0::Continental {
                        log.eroded_m3[b] += step.bedrock_removed_m[c] * a;
                    }
                }
                Some(log)
            }
            None => None,
        };
        self.sediment_ledger.record_epoch_m3(
            t_new,
            step.eroded_m3,
            step.deposited_m3,
            step.exported_m3,
        )?;
        self.absolute_drift += ledger.residual_m3().abs();
        let marine = routing
            .water_bodies
            .iter()
            .find(|b| b.kind == pf_contracts::domain::WaterBodyKindV0::Marine)
            .expect("classified marine reservoir");
        let mut transfers = step.transfers.clone();
        for (i, to) in routing.flow_to.iter().enumerate() {
            if let Some(j) = to {
                let j = j.0 as usize;
                if routing.water_body_id[j] != u32::MAX
                    && routing.water_body_id[j] != routing.water_body_id[i]
                {
                    transfers.push(pf_contracts::domain::WaterTransferV0 {
                        from_element: pf_contracts::spatial::ElementId(i as u64),
                        to_element: pf_contracts::spatial::ElementId(j as u64),
                        receiving_body: routing.water_body_id[j],
                        kind: "runoff_rate".into(),
                        sediment_m3: 0.,
                        runoff_m3_s: routing.discharge_proxy_m3_s[i],
                    });
                }
            }
        }
        let mut totals = BTreeMap::new();
        for t in &transfers {
            *totals
                .entry(format!("{}:{}", t.receiving_body, t.kind))
                .or_insert(0.0) += t.sediment_m3;
        }
        self.exchange_epochs
            .push(pf_contracts::domain::WaterExchangeEpochV0 {
                arc_normals: ledger.arc_normals.clone(),
                start_my: t_new - dt,
                end_my: t_new,
                sediment_density_kg_m3: rho_c,
                ocean_seed: marine.seed_element,
                ocean_area_m2: marine.area_m2,
                water_bodies: routing.water_bodies.clone(),
                transfers,
                sediment_totals_m3: totals,
                crust_advection_residual_m3: ledger.residual_m3(),
                crust_advection_absolute_cumulative_m3: self.absolute_drift,
                wedge_in_m3: ledger.wedge_in_m3,
                wedge_stacked_m3: ledger.wedge_stacked_m3,
                wedge_delaminated_m3: ledger.wedge_delaminated_m3,
                spread_bound_hits: step.spread_bound_hits,
                global_marine_layer_m3: step.deposit_aggraded_m3,
            });
        let bedrock_final_m3 = self.state.crust_volume_m3(areas);
        let sediment_final_m3 = volume_m3(&self.state.sediment_m, areas);

        // SEED-CARRY-CONSERVE: hand the step's changes of continental cells
        // back to the parcels (conservative), births become parcels; the
        // ledger of this step.
        let parcel_projection = ledger.parcel_projection.take();
        let mut strata_log: Option<StratigraphyStepLog> = None;
        let conserve_log: Option<ConserveStepLog> = if env.conserve_active() {
            let mut log = ConserveStepLog {
                parcels: 0,
                parcel_area_m2: 0.0,
                parcel_volume_m3: 0.0,
                parcel_sediment_m3: 0.0,
                grid_cont_cells: 0,
                grid_cont_area_m2: 0.0,
                grid_cont_volume_m3: 0.0,
                consumed_parcels: 0,
                consumed_area_m2: 0.0,
                consumed_m3: 0.0,
                consumed_sediment_m3: 0.0,
                overlap_cells: 0,
                created_parcels: 0,
                created_area_m2: 0.0,
                created_m3: 0.0,
                created_sediment_m3: 0.0,
                writeback_m3: 0.0,
                writeback_sediment_m3: 0.0,
                clamped_m3: 0.0,
                clamped_sediment_m3: 0.0,
                volume_residual_m3: 0.0,
                sediment_residual_m3: 0.0,
                overridden_oceanic_cells: 0,
                overridden_births: 0,
                vacated_copied: 0,
                vacated_births: 0,
                lost_cells: 0,
                fragments_folded: 0,
                fragment_cells: 0,
                fragment_parcels: 0,
                sutures: sutures_now.len() as u64,
                max_contact_m,
            };
            if let (Some(layer), Some(proj), Some(pr)) =
                (parcel_layer.as_mut(), parcel_projection.as_ref(), ledger.parcels.as_ref())
            {
                // SEED-RANGES (revision 1): removals by parcel thickness under
                // the plateau law (cycle_ranges profiles); write_back otherwise.
                let weighted = env.plateau.is_some() || env.plateau_in_belt.is_some();
                let wb = if let Some(sp) = env.strata.as_ref() {
                    // SEED-STRATIGRAPHY (lane ST): the same write-back, then
                    // the stacks on its kernel weights.
                    let tags = strata_tags(grid, sp, &ledger, strata_backarc_dh.as_deref(), &step, &elev_tect, sea, &ocean, &routing, dt)?;
                    let (wb, sl) = if env.parcel_residual_iterations > 0 {
                        // HYPSO-LEVERS (opt-in): the same write-back plus the residual correction.
                        pf_geodynamics::conserve::write_back_strata_corrected(grid, areas, layer, proj, &self.state, weighted, env.parcel_residual_iterations, &tags, sp, t_new - dt)?
                    } else {
                        pf_geodynamics::conserve::write_back_strata(grid, areas, layer, proj, &self.state, weighted, &tags, sp, t_new - dt)?
                    };
                    if let Some(h)=layer.history.as_mut() {
                        h.ledger(self.state.t_my,&sl,&pr.consumed_strata_m3).map_err(PreviewError::Domain)?;
                    }
                    strata_log = Some(StratigraphyStepLog::of(&sl, &pr.consumed_strata_m3, sp.v2.is_some()));
                    wb
                } else if env.parcel_residual_iterations > 0 {
                    pf_geodynamics::conserve::write_back_corrected(grid, areas, layer, proj, &self.state, weighted, env.parcel_residual_iterations)?
                } else if weighted {
                    write_back_thickness_weighted(grid, areas, layer, proj, &self.state)?
                } else {
                    write_back(grid, areas, layer, proj, &self.state)?
                };
                log.parcels = wb.parcels_after;
                log.parcel_area_m2 = wb.area_after_m2;
                log.parcel_volume_m3 = wb.volume_after_m3;
                log.parcel_sediment_m3 = wb.sediment_after_m3;
                log.grid_cont_cells = pr.projected_cells;
                log.grid_cont_area_m2 = pr.projected_area_m2;
                log.grid_cont_volume_m3 = pr.projected_volume_m3;
                log.consumed_parcels = pr.consumed_parcels;
                log.consumed_area_m2 = pr.consumed_area_m2;
                log.consumed_m3 = pr.consumed_m3;
                log.consumed_sediment_m3 = pr.consumed_sediment_m3;
                log.overlap_cells = pr.overlap_cells;
                log.created_parcels = wb.created_parcels;
                log.created_area_m2 = wb.created_area_m2;
                log.created_m3 = wb.created_m3;
                log.created_sediment_m3 = wb.created_sediment_m3;
                log.writeback_m3 = wb.writeback_m3;
                log.writeback_sediment_m3 = wb.writeback_sediment_m3;
                log.clamped_m3 = wb.clamped_m3;
                log.clamped_sediment_m3 = wb.clamped_sediment_m3;
                log.volume_residual_m3 = wb.volume_after_m3
                    - (pr.volume_before_m3 - pr.consumed_m3 + wb.writeback_m3 + wb.created_m3
                        + wb.clamped_m3);
                log.sediment_residual_m3 = wb.sediment_after_m3
                    - (pr.sediment_before_m3 - pr.consumed_sediment_m3
                        + wb.writeback_sediment_m3
                        + wb.created_sediment_m3
                        + wb.clamped_sediment_m3);
                log.overridden_oceanic_cells = pr.overridden_oceanic_cells;
                log.overridden_births = pr.overridden_births;
                log.vacated_copied = pr.vacated_copied;
                log.vacated_births = pr.vacated_births;
                log.lost_cells = wb.lost_cells;
                log.fragments_folded = pr.fragments_folded;
                log.fragment_cells = pr.fragment_cells;
                log.fragment_parcels = pr.fragment_parcels;
            }
            Some(log)
        } else {
            None
        };
        if let Some(layer) = parcel_layer.take() {
            if let Some(cs) = self.orogen.as_mut().and_then(|o| o.conserve.as_mut()) {
                cs.parcels = layer;
            }
        }
        // The grid's continental ledger (every law; feeds no output).
        let mut continental = ContinentalLedgerV1::default();
        for c in 0..n {
            let a = areas[c];
            if self.state.kind[c] == CrustKindCodeV0::Continental {
                continental.cont_area_m2 += a;
                continental.cont_volume_m3 += a * self.state.thickness_m[c];
                continental.sediment_cont_m3 += a * self.state.sediment_m[c];
            } else {
                continental.sediment_ocean_m3 += a * self.state.sediment_m[c];
            }
        }

        // SEED-OROGEN bookkeeping and the new subgrid relief.
        if let Some(o) = self.orogen.as_mut() {
            let mut log = OrogenStepLog {
                t_my: t_new,
                belt_in_m3: ledger.belt_in_m3,
                collision_overlap_consumed_m3: ledger.collision_overlap_consumed_m3,
                collision_supply_edges: ledger.collision_supply_edges,
                belt_placed_m3: ledger.belt_placed_m3,
                belt_sediment_m3: ledger.belt_sediment_m3,
                belt_delaminated_m3: ledger.belt_delaminated_m3,
                spread_substeps: ledger.spread.as_ref().map(|r| r.substeps).unwrap_or(0),
                spread_moved_m3: ledger.spread.as_ref().map(|r| r.moved_m3).unwrap_or(0.0),
                spread_residual_m3: ledger.spread.as_ref().map(|r| r.residual_m3).unwrap_or(0.0),
                spread_max_kappa_m2_per_my: ledger
                    .spread
                    .as_ref()
                    .map(|r| r.max_kappa_m2_per_my)
                    .unwrap_or(0.0),
                arc_created_m3: ledger.arc_created_m3,
                arc_placed_m3: ledger.arc_placed_m3,
                max_thickness_m: self.state.thickness_m.iter().cloned().fold(0.0, f64::max),
                subgrid_target_m3: 0.0,
                subgrid_realised_m3: 0.0,
                subgrid_max_abs_dev_m: 0.0,
                subgrid_capped_cells: 0,
                subgrid_trunk_raised_cells: 0,
                subgrid_mean_land_relief_m: 0.0,
                subgrid_max_relief_m: 0.0,
                carry: (self.env.carry == CarryLaw::FrontRescueV1).then(|| CarryStepLog {
                    rescues: ledger.carry.rescues_admitted,
                    cc_losers: ledger.carry.cc_overlap_losers,
                    cc_losers_carried: ledger.carry.cc_losers_carried,
                    consumed_parcels: ledger.carry.cc_consumed_parcels,
                    consumed_area_m2: ledger.carry.cc_consumed_area_m2,
                    consumed_m3: ledger.carry.cc_consumed_m3,
                    consumed_sediment_m3: ledger.carry.cc_consumed_sediment_m3,
                    continental_dropped_parcels: ledger.carry.continental_dropped_parcels,
                    continental_dropped_m3: ledger.carry.continental_dropped_m3,
                    sliver_births: ledger.carry.sliver_births,
                }),
                extras: extras_carried.as_ref().map(|_| ExtrasStepLog {
                    hotspot_created_m3: ledger.hotspot_created_m3,
                    hotspot_head_m3: ledger.hotspot_head_m3,
                    plume_heads: ledger.hotspot_heads,
                    rift_births: ledger.rift_births.len() as u64,
                    rift_breakup_births: ledger.rift_breakup_births,
                    rift_convergent_skipped: ledger.rift_convergent_skipped,
                    rift_moved_m3: ledger.rift_moved_m3,
                    plate_splits: splits_now.len() as u64,
                    plates: self.plates.len() as u64,
                    rift_flank: flank.as_ref().map(|f| f.1.clone()),
                    flood_basalt: swell.as_ref().map(|s| FloodBasaltStepLog {
                        lip_created_m3: ledger.lip_created_m3,
                        lip_continental_m3: ledger.lip_continental_m3,
                        lip_cells: ledger.lip_cells.len() as u64,
                        max_swell_wl_m: s.iter().fold(0.0, |a: f64, &b| a.max(b)),
                        swell_cells: s.iter().filter(|&&v| v > 100.0).count() as u64,
                    }),
                }),
                trench: trench.as_ref().map(|t| t.1.clone()),
                conserve: conserve_log.clone(),
                stratigraphy: strata_log.clone(),
            };
            if let Some((sg, _, wet)) = &subgrid {
                let mut relief = vec![0.0f64; n];
                let (mut land_a, mut land_r) = (0.0, 0.0);
                for c in 0..n {
                    let t = sg.target_lowering_m[c];
                    if t > 0.0 {
                        let realised = step.bedrock_removed_m[c] + step.sediment_removed_m[c];
                        log.subgrid_target_m3 += t * areas[c];
                        log.subgrid_realised_m3 += realised * areas[c];
                        log.subgrid_max_abs_dev_m = log.subgrid_max_abs_dev_m.max((realised - t).abs());
                    }
                    if !wet[c] {
                        relief[c] = ((elev_tect[c] + step.dz_m[c]) - sg.trunk_m[c]).max(0.0);
                        land_a += areas[c];
                        land_r += relief[c] * areas[c];
                        log.subgrid_max_relief_m = log.subgrid_max_relief_m.max(relief[c]);
                    }
                }
                log.subgrid_capped_cells = sg.capped_cells;
                log.subgrid_trunk_raised_cells = sg.trunk_raised_cells;
                log.subgrid_mean_land_relief_m = if land_a > 0.0 { land_r / land_a } else { 0.0 };
                o.subgrid_relief_m = Some(relief);
            }
            if extras_carried.is_some() {
                o.extras = extras_carried;
            }
            o.steps.push(log);
        }

        // 4. Isostatic re-equilibration after erosion.
        let (mut e_after, solve) =
            elevation_on_grid_terms(&self.state, params, areas, flank_wl, trench_wl)?;
        if let Some((t, f)) = heat_dynamic {
            add_dynamic_term(&mut e_after, t, sea, f);
        }
        self.elevation = e_after;
        let elevation = &self.elevation;
        let (mut mx, mut mn, mut land) = (f64::MIN, f64::MAX, 0.0);
        let total_area: f64 = areas.iter().sum();
        for c in 0..n {
            mx = mx.max(elevation[c]);
            mn = mn.min(elevation[c]);
            if elevation[c] >= sea {
                land += areas[c];
            }
        }
        continental.land_area_m2 = land;
        self.logs.push(StepLog {
            t_my: t_new,
            ridge_cells: ledger.ridge_cells,
            subduction_cells: ledger.subduction_cells,
            collision_cells: ledger.collision_cells,
            advection_residual_m3: ledger.residual_m3(),
            eroded_m3: step.eroded_m3,
            deposited_m3: step.deposited_m3,
            exported_m3: step.exported_m3,
            diffusion_substeps: step.diffusion_substeps,
            deposit_aggraded_m3: step.deposit_aggraded_m3,
            max_elevation_m: mx,
            min_elevation_m: mn,
            land_fraction: land / total_area,
            sea_level_gauge_m: solve.map(|s| s.sea_level_gauge_m),
            flooded_area_fraction: solve.map(|s| s.flooded_area_fraction),
            ocean_volume_residual_m3: solve.map(|s| s.ocean_volume_residual_m3),
            compensation_depth_m: solve.map(|s| s.compensation_depth_m),
            deepest_root_m: solve.map(|s| s.deepest_root_m),
            deepest_root_crust_sediment_m: solve.map(|s| s.deepest_root_crust_sediment_m),
            land_elevation_p10_p50_p90_m: if env.flags_on {
                land_percentiles(elevation, areas, sea, [0.1, 0.5, 0.9])
            } else {
                None
            },
            deposition: step.deposition.as_ref().map(StepDeposition::from),
        });
        Ok(PreviewStepReport {
            t_my: t_new,
            boundary_upserts: upserts,
            bedrock_before_m3,
            bedrock_after_advance_m3,
            bedrock_final_m3,
            sediment_before_m3,
            sediment_after_advance_m3,
            sediment_final_m3,
            bedrock_eroded_m3: step.bedrock_eroded_m3,
            sediment_eroded_m3: step.sediment_eroded_m3,
            eroded_m3: step.eroded_m3,
            deposited_m3: step.deposited_m3,
            exported_m3: step.exported_m3,
            deposit_aggraded_m3: step.deposit_aggraded_m3,
            deposition: step.deposition.clone(),
            tectonic_solve,
            solve,
            advance: ledger,
            continental,
            conserve: conserve_log,
            stratigraphy: strata_log,
            reorg: reorg_report,
            ranges: if indent_report.is_some() || plateau_report.is_some() {
                Some(RangesStepLog::of(indent_report.as_ref(), plateau_report.as_ref()))
            } else {
                None
            },
            shortening_budget,
        })
    }

    /// The tail of the chain: crust ledger closure, final hydrology,
    /// terrain and resource domains. Refuses an incomplete schedule.
    pub fn finish(self) -> Result<PreviewOutputs, PreviewError> {
        let native_product = self.native_evidence.as_ref().map(|e| e.product(&self)).transpose()?;
        self.carrier_identity()?;
        self.validate_carrier_history()?;
        if let Some(trace)=self.fluvial_trace()?{trace.validate_all().map_err(PreviewError::Domain)?;}
        if self.steps_taken() != self.env.steps {
            return Err(PreviewError::Domain(format!(
                "finish needs the whole step schedule: {} of {} steps taken",
                self.steps_taken(),
                self.env.steps
            )));
        }
        let PreviewRun {
            native_evidence: _,
            env,
            plates,
            state,
            bodies_report,
            elevation,
            elevation_t0,
            sea_level_t0,
            uplift_rate,
            sediment_ledger,
            boundary_table,
            mut crust,
            nu_max_observed,
            substeps_max,
            max_uplift,
            max_erosion_rate,
            logs,
            exchange_epochs,
            absolute_drift: _,
            orogen,
            climate: _,
        } = self;
        let PreviewEnv {
            params,
            grid,
            seed_path,
            prior,
            ep,
            areas,
            sea,
            rho_c,
            flags_on,
            dt,
            steps,
            rain,
            belt,
            subgrid_x_c_m,
            relief_law,
            critical_slope,
            retention,
            surface_laws,
            carry,
            collision_sediment,
            plumes,
            rift,
            split_rate_deg_per_my,
            rift_flank,
            trench,
            seeded_margin_subsidence,
            parcel_params,
            parcel_residual_iterations: _,
            suture_contact_m,
            flood_basalt,
            reorg,
            reorg_law,
            indent,
            plateau,
            magma,
            plateau_in_belt,
            reorg_initiation_age_my,
            magma_heat_dry_factor,
            boundary_transform,
            reorg_limits,
            reorg_margin,
            strata,
            climate: climate_law,
            ..
        } = env;
        let n = grid.cells();
        crust.final_m3 = state.crust_volume_m3(&areas);
        // advection_residual_m3 is the sum of the per-step measured residuals
        // (accumulated above); whatever the named terms still fail to explain
        // is the closure plug, kept apart so it is visible (F6).
        crust.closure_plug_m3 = crust.final_m3
            - (crust.initial_m3 + crust.ridge_created_m3 + crust.arc_created_m3
                - crust.subducted_m3
                - crust.delaminated_m3
                - crust.eroded_m3
                + crust.advection_residual_m3);

        // --- final hydrology on the final surface ------------------------------
        let ocean = pf_surface::water::largest_ocean(&grid, &elevation, sea)?;
        let routing = route(&grid, &elevation, &rain, &ocean, sea)?;
        let (mut hydrology, cell_reach) = assemble_hydrology_with_cells(&grid, &routing, state.t_my)?;
        hydrology.cell_fields = Some(HydrologyCellFieldsV0::v0());
        hydrology.exchange_history = Some(pf_contracts::domain::ExchangeHistoryV0 {
            sidecar: "surface-exchange".into(),
            epochs: exchange_epochs.len() as u64,
        });
        hydrology
            .validate()
            .map_err(|e| PreviewError::Domain(e.to_string()))?;

        // --- terrain domain -----------------------------------------------------
        let total_area: f64 = areas.iter().sum();
        let (mut mx, mut mn, mut land, mut thick_sum) = (f64::MIN, f64::MAX, 0.0, 0.0);
        for c in 0..n {
            mx = mx.max(elevation[c]);
            mn = mn.min(elevation[c]);
            if elevation[c] >= sea {
                land += areas[c];
            }
            thick_sum += state.thickness_m[c] * areas[c];
        }
        let mut boundaries: Vec<BoundaryRecordV0> = boundary_table.active().cloned().collect();
        boundaries.extend(boundary_table.retired().iter().cloned());
        boundaries.sort_by_key(|b| b.boundary_id);
        let mut assumptions = prior.assumptions();
        assumptions.extend(driver_assumptions(&params, &ep));
        let mut approximations: Vec<String> = vec![
                "rigid-plates-seeded-euler-poles-never-coupled".into(),
                "eulerian-semi-lagrangian-crust-carry".into(),
                "kinematic-overlap-thickening-prior".into(),
                "airy-crust-only-isostasy-with-water-loading".into(),
                "gdh1-oceanic-thermal-subsidence".into(),
                "implicit-stream-power-n1".into(),
                "explicit-linear-hillslope-diffusion".into(),
                "uniform-rainfall".into(),
                "no-mantle-no-melt-no-element-mass".into(),
            ];
        if flags_on {
            flag_assumptions(&params, &prior, &mut assumptions, &mut approximations);
        }
        if !surface_laws.is_legacy() {
            sediment_law_assumptions(&surface_laws, &mut assumptions, &mut approximations);
        }
        let mut uplift_record = prior.record.clone();
        let mut stream_power_scheme = "implicit-receiver-first-braun-willett-2013".to_string();
        if belt.is_some() || subgrid_x_c_m.is_some() {
            orogen_assumptions(
                belt.as_ref(),
                subgrid_x_c_m,
                relief_law.as_ref(),
                &ep,
                &grid,
                &mut assumptions,
                &mut approximations,
            );
        }
        if belt.is_some() {
            uplift_record.rule = "rate-limited-belt-v1".into();
        }
        if carry == CarryLaw::FrontRescueV1 {
            // SEED-CARRY-CC (lane CC): the carry the run used.
            for a in approximations.iter_mut() {
                if a == "eulerian-semi-lagrangian-crust-carry" {
                    *a = "eulerian-semi-lagrangian-crust-carry-with-continental-front-rescue-v1".into();
                }
            }
            assumptions.push(AssumptionV0 {
                id: "A-ADVECTION-FRONT-RESCUE".into(),
                statement: "front-rescue-carry-v1: a continental parcel the semi-Lagrangian carry would carry to no cell is moved (never copied) into the one cell whose back-traced point lies in its cell-shaped patch; continental losers of continent-continent overlaps carried nowhere are consumed once and their volume and sediment are the belt's collision supply (replacing the kinematic-flux supply)".into(),
                value: 0.0,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane CC 2026-10-12, no constant; pre-registered in unpublished-design-record; the interior raster churn of the legacy carry (drops and duplicates) is unchanged and still reported in the advection residual".into(),
            });
        }
        if let Some(pp) = parcel_params {
            // SEED-CARRY-CONSERVE (lane CV): the carry the run used.
            for a in approximations.iter_mut() {
                if a == "eulerian-semi-lagrangian-crust-carry" {
                    *a = "lagrangian-continental-parcels-projected-on-eulerian-oceanic-carry-v1".into();
                }
            }
            assumptions.push(AssumptionV0 {
                id: "A-ADVECTION-CONTINENTAL-PARCELS".into(),
                statement: "continental-parcel-carry-v1: continental crust is a set of Lagrangian parcels (area, thickness, sediment, plate) rotated rigidly every step; the grid is their tent-kernel projection (continental where half a cell is covered); every change a step makes to a continental cell is written back to its parcels in proportion to their deposit; area changes only by rift births (new parcels) and collision consumption (two plates each covering half a cell with total above 1 + this value: the smaller depositors' parcels become the belt's collision supply); plate-id pieces below the fragment share of the sphere fold into the plate they touch most; oceanic crust keeps the semi-Lagrangian carry".into(),
                value: pp.overlap_cell_fraction,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: format!("lane CV 2026-10-13, pre-registered in unpublished-design-record; fragment fold below {} of the sphere; the projection averages columns over about one cell", pp.fragment_min_area_fraction),
            });
        }
        if let Some(l) = suture_contact_m {
            // SEED-CARRY-CONSERVE (lane CV): the suture law the run used.
            for a in approximations.iter_mut() {
                if a == "rigid-plates-seeded-euler-poles-never-coupled" {
                    *a = "rigid-plates-seeded-euler-poles-sutured-by-continental-contact-v1".into();
                }
            }
            assumptions.push(AssumptionV0 {
                id: "A-SUTURE-CONTACT-LENGTH".into(),
                statement: "contact-length-pole-share-v1: when the converging continent-continent contact between two plates reaches this length, the plate group with the smaller area takes the larger group's Euler pole for the rest of the run (convergence across the suture stops)".into(),
                value: l,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane CV 2026-10-13, pre-registered; recalled Himalayan front ~2400-2500 km, Alpine ~1200 km (not checked); no pole re-solve, no slab pull".into(),
            });
        }
        if let Some(ip) = &indent {
            // SEED-RANGES (lane RG): the suture law is the indentation version.
            for a in assumptions.iter_mut() {
                if a.id == "A-SUTURE-CONTACT-LENGTH" {
                    a.statement = "contact-length-indentation-taper-v1: when the converging continent-continent contact between two plates reaches this length, the two plate groups suture (union-find) and the smaller group's member plates keep their relative rotation at the suture, scaled by a linear taper to zero over the indentation time or until the shortening budget at their contact is spent (post-collision indentation; convergence tapers instead of stopping)".into();
                    a.note = format!("lane RG 2026-10-13, pre-registered in unpublished-design-record; contact length from lane CV (recalled Himalayan front ~2400-2500 km); indentation time {} My, shortening budget {} m (India-Asia ~50 My and ~2000 km, Alps ~500 km; recalled, not checked)", ip.time_my, ip.budget_m);
                }
            }
            assumptions.push(AssumptionV0 {
                id: "A-INDENTATION-TIME".into(),
                statement: "contact-length-indentation-taper-v1: time over which the smaller group's relative rotation at the suture tapers linearly to zero".into(),
                value: ip.time_my,
                units: "My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane RG 2026-10-13, pre-registered; India-Asia collision ~50 Ma and still converging, Alpine convergence slowing over ~30 My (recalled, not checked); design review review 'try 40'".into(),
            });
            assumptions.push(AssumptionV0 {
                id: "A-INDENTATION-BUDGET".into(),
                statement: "contact-length-indentation-taper-v1: cumulative shortening (length-weighted mean converging normal rate across the indentor's contact with the rest of its group, times dt) that ends an indentation early".into(),
                value: ip.budget_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane RG 2026-10-13, pre-registered; Alps ~500 km, India-Asia ~2000 km including subducted Greater India (recalled, not checked); midpoint of design review's 500-1000 km".into(),
            });
        }
        if let Some(pp) = &plateau {
            // SEED-RANGES (lane RG): the plateau law.
            assumptions.push(AssumptionV0 {
                id: "A-PLATEAU-CAP".into(),
                statement: "plateau-cap-spill-v1: continental crust above this thickness after the belt step gives the excess to the nearest same-plate continental cells below it (breadth-first rings of edge neighbours through cells already at the cap, each filled up to the cap, ascending cell order), within the reach; the remainder and oceanic crust above it are delaminated (named ledger term)".into(),
                value: pp.cap_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: format!("lane RG 2026-10-13, pre-registered; set by plateau height: Tibetan crust ~65-75 km stands at ~5 km (Zhang and Klemperer 2005; Molnar and Lyon-Caen 1988; recalled, not checked) and this model's Airy column puts 65 km at ~5.2 km; reach {} m (Tibet ~1000 km N-S)", pp.reach_m),
            });
        }
        // SEED-WORLD-COMPOSE (lane CP): the v2 plateau law, the v3 initiation
        // and the v2 magmatism changes the run used.
        if let Some(pp) = &plateau_in_belt {
            assumptions.push(AssumptionV0 {
                id: "A-PLATEAU-CAP-V2".into(),
                statement: "plateau-cap-spill-v2: inside the belt step, after the collision supply is placed and before the channel-flow spreading (which then never sees crust above the cap): sediment of continental columns whose crust plus sediment exceeds this cap is accreted into the crust, crust above it spills to the nearest same-plate continental cells as in v1 (remainder delaminated), then the foothill taper; the spreading diffusivity is scaled per cell by min(1, (2 dx/dy)^2) (polar limiter, numerical: slows channel flow poleward of 60 deg only)".into(),
                value: pp.cap_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: format!("lane CP 2026-10-13, pre-registered in unpublished-design-record (lane RG ideas 1-3); cap and reach {} m as v1", pp.reach_m),
            });
            assumptions.push(AssumptionV0 {
                id: "A-PLATEAU-FOOTHILL-GRADIENT".into(),
                statement: format!("plateau-cap-spill-v2 foothill taper: the crust thickness step between edge-neighbouring same-plate continental cells is at most this x their centre distance wherever the thicker cell is above {} m (volume moved down-step, conservative, ascending Gauss-Seidel sweeps)", pp.taper_floor_m),
                value: pp.foothill_gradient,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane CP 2026-10-13; Tibet's crust 70 -> 40 km over ~200-300 km at its margins, Subandean ramp ~200-300 km (recalled, not checked)".into(),
            });
        }
        // SEED-TRANSFORMS (lane TF): the plate-boundary law and the v4 limits.
        if let Some(tp) = boundary_transform {
            assumptions.push(AssumptionV0 {
                id: "A-TRANSFORM-OBLIQUITY".into(),
                statement: "smoothed-normal-transform-v1: a moving plate boundary whose relative velocity is within this angle of the boundary trace is a transform: gaps on it are filled from the same side's carried material and oceanic overlaps on it (or on a divergent boundary) are retired without subduction, both named remap terms; convergent gaps are filled the same way (no ridge birth where the plates converge)".into(),
                value: tp.transform_obliquity_deg,
                units: "deg".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane TF 2026-10-13, pre-registered in unpublished-design-record (lane VL defect 1; VL validated at 30 deg, 20 and 45 as sensitivity)".into(),
            });
            assumptions.push(AssumptionV0 {
                id: "A-BOUNDARY-NORMAL-SMOOTHING".into(),
                statement: "smoothed-normal-transform-v1: Gaussian sigma of the plate-indicator window whose first moment gives the boundary normal (window cut at 3 sigma, cell-area weighted; the staircase direction where degenerate or pointing away)".into(),
                value: tp.normal_smoothing_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane TF 2026-10-13; ~1.5 L0 cells, as lane VL's validation normal".into(),
            });
        }
        if let Some(l) = reorg_limits {
            assumptions.push(AssumptionV0 {
                id: "A-INITIATION-SLAB-AGE".into(),
                statement: format!("plate-reorg-scripted-v4: an age-triggered initiation cut is rejected unless at least {} of its new boundary length runs where the older flanking cell is ocean at least this old (My)", l.old_slab_fraction_min),
                value: l.slab_age_min_my,
                units: "My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane TF 2026-10-13 (lane CP: too many intra-oceanic arcs); spontaneous initiation needs old dense lithosphere (Stern 2004, recalled, not checked)".into(),
            });
            assumptions.push(AssumptionV0 {
                id: "A-INITIATION-SPACING".into(),
                statement: "plate-reorg-scripted-v4: an age-triggered initiation cut is rejected when more than half of its new boundary runs within this distance of a subduction-class cell or of a boundary started earlier in the same pass".into(),
                value: l.trench_spacing_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane TF 2026-10-13; arc systems ~1000 km or more apart on Earth (recalled, not checked)".into(),
            });
        }
        // SEED-OLD-OCEAN (lane OR): the v5 margin initiation.
        if let Some(mi) = reorg_margin {
            assumptions.push(AssumptionV0 {
                id: "A-MARGIN-INITIATION-AGE".into(),
                statement: "plate-reorg-scripted-v5: at each scheduled event, while fewer than the maximum groups, the passive-margin ocean cell with the most ocean at least this old within 1000 km (oceanic, same suture group as a continental neighbour, touching no other group) at least this old and away from existing subduction seeds a trench along the margin; the ocean whose nearest coast is that stretch becomes a plate converging on the continent at the subduction initiation rate (at most 3 per event, one per group, side B at least 0.1% of the sphere)".into(),
                value: mi.age_my,
                units: "My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane OR 2026-10-13, pre-registered in unpublished-design-record (lane VL defect 3); ocean floor older than ~160-180 My is almost absent on Earth (Mueller et al. 2008, recalled, not checked)".into(),
            });
            assumptions.push(AssumptionV0 {
                id: "A-MARGIN-SEGMENT-LENGTH".into(),
                statement: "plate-reorg-scripted-v5: the new trench runs along at most this length of margin (margin cells at least 100 My old reachable from the seed along the margin)".into(),
                value: mi.segment_length_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane OR 2026-10-13; Lesser Antilles ~850 km, Izu-Bonin-Mariana ~2800 km (recalled, not checked)".into(),
            });
            assumptions.push(AssumptionV0 {
                id: "A-MARGIN-TRENCH-SPACING".into(),
                statement: "plate-reorg-scripted-v5: the seed and the segment keep at least this distance from subduction-class cells and from trenches started earlier in the same pass".into(),
                value: mi.spacing_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane OR 2026-10-13; as A-INITIATION-SPACING".into(),
            });
        }
        if let Some(age) = reorg_initiation_age_my {
            assumptions.push(AssumptionV0 {
                id: "A-SUBDUCTION-INITIATION-AGE".into(),
                statement: "plate-reorg-scripted-v3: at each scheduled event a group whose ocean older than this covers at least the counted-plate area starts a trench through that old ocean (oldest-area group first, while below the plate-count maximum); ocean weakens with age from 100 to 180 My, the coastal halo is at least craton strength and continents 10x stronger, so trenches open on the ocean side of passive margins and in enclosed basins".into(),
                value: age,
                units: "My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane CP 2026-10-13 (lane VL defect 3, lane BC idea 1); oldest in-situ ocean floor ~180 My, ocean-age p95 ~150-180 My (recalled, not checked)".into(),
            });
        }
        if matches!(universe_of(&params)?.law_set.magmatism,
            pf_contracts::laws::MagmatismLawV1::SlabArcBackarcLlsvpV3 { .. }) {
            assumptions.push(AssumptionV0 {
                id: "A-ARC-HEAT-WATER-CLOSURE-V3".into(),
                statement: "Arc thermal buoyancy enters the local-Airy water-loaded column before both sea-level solves; the prescribed water inventory includes its displacement. V2 back-arc episode/stretch and head-only resource rules remain.".into(),
                value: 1.0,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "Teaching opt-in coupling policy. No new observational calibration claim; numerical inventory closure must pass the teaching water gates.".into(),
            });
        }
        if let (Some(f), Some(me)) = (magma_heat_dry_factor, magma.as_ref()) {
            assumptions.push(AssumptionV0 {
                id: "A-ARC-HEAT-DYNAMIC".into(),
                statement: "slab-arc-backarc-llsvp-v2: the arc heat is a dynamic elevation term added after the sea-level solve (water-loaded under the sea, times this dry factor above it); it displaces no ocean water because the ocean volume constant is calibrated on a present-day hypsometry that already contains arc and back-arc thermal highs (v1 added it to the column, raising sea level ~80 m); back-arc extension is one episode per cell and the kimberlite marker head-only".into(),
                value: f,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: format!("lane CP 2026-10-13 (lane MG ideas 1, 2, 4); back-arc episode {} My, stretching factor {}, kimberlite radius {} m", me.params.v2.map_or(0.0, |v| v.backarc_episode_my), me.params.v2.map_or(0.0, |v| v.backarc_max_stretch_factor), me.params.v2.map_or(0.0, |v| v.kimberlite_plume_radius_m)),
            });
        }
        if let Some((rp, schedule)) = &reorg {
            // SEED-PLATE-REORG (lane PR): the scripted Wilson cycle the run used.
            for a in approximations.iter_mut() {
                if a == "rigid-plates-seeded-euler-poles-sutured-by-continental-contact-v1" {
                    *a = "rigid-plates-scripted-wilson-cycle-v1".into();
                }
            }
            assumptions.push(AssumptionV0 {
                id: "A-PLATE-REORG-SCRIPTED".into(),
                statement: if reorg_law == ReorgLaw::ScriptedV5 { "plate-reorg-scripted-v5: plate-reorg-scripted-v4 plus passive-margin subduction initiation driven by ocean age: the passive-margin ocean with the most old ocean nearby seeds a trench along a limited stretch of that margin, the frontal ocean of that margin becoming a plate that converges on the continent; re-solves keep, per group, the candidate axis (8, two sweeps) closing the least ridge or young-young boundary (OR-WORDING-R1) (A-MARGIN-INITIATION-AGE, A-MARGIN-SEGMENT-LENGTH, A-MARGIN-TRENCH-SPACING)" } else if reorg_law == ReorgLaw::ScriptedV4 { "plate-reorg-scripted-v4: plate-reorg-scripted-v3 with the age-triggered initiation limited to cuts whose new boundary runs mostly through old slab and stays away from existing subduction (A-INITIATION-SLAB-AGE, A-INITIATION-SPACING)" } else if reorg_law == ReorgLaw::ScriptedV3 { "plate-reorg-scripted-v3: plate-reorg-scripted-v2 (continental cut domain, peel rejection, symmetric ridge accretion, folded retirement) plus age-driven subduction initiation on the v2 integer cut (A-SUBDUCTION-INITIATION-AGE)" } else if reorg_law == ReorgLaw::ScriptedV2 { "plate-reorg-scripted-v2: scripted kinematics, not a force balance. As plate-reorg-scripted-v1 (seeded re-solves with slab-torque axes, supercontinent breakup after the tenure, plume-head splits, subduction initiation, folds, plume renewal), but breakups and plume splits cut only the group's continent plus a one-cell halo at least as strong as craton, weakened by cumulative shortening, old sutures and the plume head (no thickness term); the ocean boundary continues from the coast exit points to the nearest-coast equidistance line; cuts leaving under 15% of the continent on a side are rejected; folded plates retire as Folded" } else { "plate-reorg-scripted-v1: scripted kinematics, not a force balance. Suture groups get new poles at seeded re-solve times (axis half random, half slab torque; rate grows with the subducting share of the boundary); a group holding the supercontinent share of the continental crust for the tenure time breaks up along a min cut through weak zones (crust far from the group's median thickness, old shortening, old sutures); plume heads split large continents along the same cut; subduction initiates at passive margins when too few plates remain; cut leftovers fold or become plates; plume onsets are a renewal process with finite tails" }.into(),
                value: rp.supercontinent_tenure_my,
                units: "My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: format!("{}lane PR 2026-10-13, pre-registered in unpublished-design-record; {} re-solves scheduled; constants recalled (Forsyth and Uyeda 1975; Bird 2003; Coffin and Eldholm 1994), not checked", if reorg_law != ReorgLaw::ScriptedV1 { "lane BC 2026-10-13 v2 cut, pre-registered in unpublished-design-record; " } else { "" }, schedule.len()),
            });
        }
        if collision_sediment == CollisionSedimentLaw::AccretedToBeltCrustV1 {
            // SEDIMENT-TOWERS (lane ST): the collision-sediment law the run used.
            assumptions.push(AssumptionV0 {
                id: "A-COLLISION-SEDIMENT-ACCRETED".into(),
                statement: "accreted-to-belt-crust-v1: sediment of a consumed continental collision loser joins the rate-limited belt's crust supply (placed, spread by channel flow and capped by the delamination backstop like crust) instead of the surface sediment layer; a loser left to the remap keeps its sediment as it keeps its crust; the accreted volume is booked as net erosion in the crust ledger".into(),
                value: 0.0,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane ST 2026-10-12, no constant; accretion of passive-margin and trench sediments into fold-thrust belts recalled, not checked; legacy surface-layer-v0 built 35-57 km sediment towers on suture cells (unpublished-design-record)".into(),
            });
        }
        // SEED-EXTRAS (lane EX): the plume, rift and split laws the run used.
        if let Some((pl, pp)) = &plumes {
            assumptions.retain(|a| a.id != "A-NO-HOTSPOTS");
            approximations.push("fixed-mantle-plumes-head-plus-age-progressive-tail-magmatic-crust-v1".into());
            assumptions.push(AssumptionV0 {
                id: "A-HOTSPOT-PLUMES".into(),
                statement: format!(
                    "fixed-plume-chains-v1: {} plumes fixed in the mantle (datum) frame (K6 hotspot stream positions, onsets uniform in [-lead, duration)); a head adds {:.3e} m3 of crust over a {:.0} m disc at onset (any plate); the tail adds {:.3e} m3/My along the step's age-progressive track on the overlying plate ({:.0} m discs). Plume crust is booked as arc_created in the crust ledger (no contract field) and named hotspot_created per step in the module control volume",
                    pl.len(), pp.head_volume_m3, pp.head_radius_m, pp.tail_flux_m3_per_my, pp.edifice_radius_m
                ),
                value: pl.len() as f64,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane EX 2026-10-12; values fixed before any run in unpublished-design-record; no thermal swell, no flexural moat, plumes never die".into(),
            });
        }
        if let Some(rp) = &rift {
            approximations.push("stretched-continental-rift-births-v1".into());
            assumptions.push(AssumptionV0 {
                id: "A-CONTINENTAL-RIFT".into(),
                statement: format!(
                    "stretched-rift-v1: a ridge birth whose nearest older cells of two diverging plates are continental is filled with continental crust stretched out of those plates' continental cells within max({:.0} m, 1.5 cell spacings) (triangular weights, volume moved exactly); a rift birth thinner than {:.0} m is oceanic instead (breakup)",
                    rp.half_width_m, rp.breakup_thickness_m
                ),
                value: rp.breakup_thickness_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: if rift_flank.is_some() {
                    "lane EX 2026-10-12; rift thermal anomaly, thermal subsidence and shoulders: A-RIFT-FLANK (lane RS)".into()
                } else {
                    "lane EX 2026-10-12; no rift-flank uplift, no thermal subsidence of stretched crust (continental columns carry no thermal term)".into()
                },
            });
        }
        if let Some(fp) = &rift_flank {
            approximations.push("rift-thermal-anomaly-plus-broken-plate-shoulders-v1".into());
            assumptions.push(AssumptionV0 {
                id: "A-RIFT-FLANK".into(),
                statement: format!(
                    "shoulder-uplift-thermal-v1: a continental column stretched by beta (stretched-rift-v1 births: flanks' weighted mean unstretched thickness / birth thickness; flank thinning multiplies beta) carries the thermal uplift (1 - 1/beta) x GDH1 thermal elevation at its thermal age (McKenzie 1978 uniform stretching; repeated stretching keeps one continuous equivalent age), so the rift floor subsides thermally with time. Each continental column gets a shoulder of {:.2} x the thermal uplift of a stretched column within reach x the cell-mean broken-plate profile exp(-x/a) cos(x/a), a = {:.0} m, x from the source cell's edge, cut at pi a / 2; the largest such shoulder counts where it exceeds the column's own thermal uplift. Elevation term only (water-loaded column term of the two-column isostasy); no crust volume moved",
                    fp.shoulder_ratio, fp.flexural_parameter_m
                ),
                value: fp.flexural_parameter_m,
                units: "m".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane RS 2026-10-13; values fixed before any run in unpublished-design-record; crust stays locally (Airy) compensated (no flexural support of the floor's crustal deficit), shoulders are a rule not a flexure solve, no plume swell".into(),
            });
        }
        if seeded_margin_subsidence {
            assumptions.push(AssumptionV0 {
                id: "A-PASSIVE-MARGIN".into(),
                statement: "seeded-margin-thermal-subsidence-v1: at t = 0 every seeded continental column thinner than the interior thickness h_cc carries the rift-flank stretching markers beta = h_cc / h and thermal start time = -initial oceanic age (breakup when the adjacent seeded ocean was born), so shoulder-uplift-thermal-v1 gives the seeded passive margins their post-rift thermal subsidence (and landward shoulders)".into(),
                value: prior.initial_oceanic_age_my,
                units: "My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane TN 2026-10-13; no new constant (PREREG.txt part 2); uniform stretching of the seeded linear taper; the t = 0 elevation does not include the term (first applied at the first step)".into(),
            });
        }
        if let Some(tp) = &trench {
            approximations.push("broken-plate-slab-pull-trench-flexure-v1".into());
            assumptions.push(AssumptionV0 {
                id: "A-TRENCH".into(),
                statement: format!(
                    "slab-flexure-trench-v1: at a Subduction-class boundary the oceanic cell overridden by a continental or younger oceanic neighbour (the kinematic overlap rule) is a trench cell; its plate (age A, normal convergence v) bends as a broken elastic plate under the end load V0 = {:.2} x slab pull, slab pull = 2 rho_m alpha_v T_m g sqrt(kappa A / pi) min({:.0} m, v 2.32^2 A / pi^2) (alpha_v {:.1e} /K, T_m {:.0} C, kappa {:.4} m2/yr); Te = 2 erfinv({:.0} / T_m) sqrt(kappa A), E {:.2e} Pa, nu {:.2}; w0 = V0 / ((rho_m - rho_w) g alpha); every oceanic cell of the downgoing plate takes the cell-mean of -w0 exp(-x/alpha) cos(x/alpha) (x from the nearest trench axis, cut at 3 pi alpha / 2) as a water-loaded column term. Elevation term only; no crust or sediment moved",
                    tp.vertical_load_fraction,
                    tp.slab_length_max_m,
                    tp.thermal_expansivity_per_k,
                    tp.mantle_temperature_c,
                    tp.thermal_diffusivity_m2_per_yr,
                    tp.elastic_isotherm_c,
                    tp.elastic_modulus_pa,
                    tp.poisson_ratio
                ),
                value: tp.vertical_load_fraction,
                units: "1".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane TN 2026-10-13; values fixed before any run in unpublished-design-record; no bending moment or yielding, no accretionary wedge or forearc response, continental cells of the downgoing plate not flexed, cell-mean at grid resolution".into(),
            });
        }
        if let Some(fb) = &flood_basalt {
            approximations.push("plume-head-lip-plateau-plus-mantle-fixed-swell-v1".into());
            assumptions.push(AssumptionV0 {
                id: "A-FLOOD-BASALT".into(),
                statement: format!(
                    "plume-head-lip-swell-v1: a plume whose onset is inside the run erupts {:.3e} m3 extrusive x (1 + {:.1}) intrusive/underplated = {:.3e} m3 of crust over {:.1} My, placed (each step's share) where the plate under the plume has carried it, over a flat-topped {:.0} m disc tapering over the outer {:.0} % (any plate); it is real crust (carried, eroded, compensated), booked with the plume crust (hotspot_created / arc_created), and supersedes the fixed-plume-chains-v1 head (run with head volume 0). Every active plume adds the mantle-fixed water-loaded swell {:.0} m exp(-(r/{:.0} m)^2) + {:.0} m exp(-(t - onset)/{:.1} My) exp(-(r/{:.0} m)^2) (cut at 3 R) to every column term. Plumes active before t = 0 erupt no LIP in the run",
                    fb.lip.extrusive_volume_m3,
                    fb.lip.intrusive_ratio,
                    fb.lip.total_volume_m3(),
                    fb.lip.eruption_duration_my,
                    fb.lip.radius_m,
                    100.0 * fb.lip.taper_fraction,
                    fb.swell.tail_amplitude_m,
                    fb.swell.tail_radius_m,
                    fb.swell.head_amplitude_m,
                    fb.swell.head_decay_my,
                    fb.swell.head_radius_m
                ),
                value: fb.lip.total_volume_m3(),
                units: "m3".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane FB 2026-10-13; values fixed before any run in unpublished-design-record; LIP crust loads at the column's crust density (no separate basalt or underplate density, no depleted-residue buoyancy), one volume for continental and oceanic heads, plumes never die so tail swells persist, no flexural moat".into(),
            });
        }
        if let Some(rate) = split_rate_deg_per_my {
            for a in approximations.iter_mut() {
                if a == "rigid-plates-seeded-euler-poles-never-coupled" {
                    *a = "rigid-plates-seeded-euler-poles-plus-plume-head-splits-never-coupled".into();
                }
            }
            assumptions.push(AssumptionV0 {
                id: "A-PLUME-HEAD-SPLIT".into(),
                statement: "plume-head-continental-split-v1: a plume head arriving under continental crust splits the plate under it along a great circle through the plume (seeded azimuth); the new plate's angular velocity is the parent's plus the opening rate about the axis on the cut 90 deg from the plume".into(),
                value: rate,
                units: "deg/My".into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: "lane EX 2026-10-12; parent poles never change; the far end of the cut (beyond 90 deg) converges".into(),
            });
        }
        if subgrid_x_c_m.is_some() {
            stream_power_scheme =
                "implicit-receiver-first-braun-willett-2013-trunk-plus-subgrid-relaxation-v1".into();
        }
        if let Some(rp) = &relief_law {
            if rp.slope_exponent_n != 1.0 {
                stream_power_scheme =
                    "implicit-receiver-first-power-root-trunk-plus-subgrid-relaxation-nonlinear-v1"
                        .into();
            }
        }
        // SEED-SEDIMENT-RETAIN (lane SR): the retention law, if selected.
        if let Some(rp) = &retention {
            approximations.push("continental-retention-base-level-v1".into());
            let pa = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
                id: id.into(),
                statement: statement.into(),
                value,
                units: units.into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: note.into(),
            };
            assumptions.push(pa("A-BASE-LEVEL-RELIEF-SCALE", "stream-power lowering of a land cell is multiplied by phi = min(1, (r/H)^p), r = local relief (max - min level over the cell and its 8 neighbours, water cells at their water level); no change at r >= H (rule 2)", rp.base_level_relief_scale_m, "m", "project assumption awaiting review; recalled cratonic denudation ~1-10 m/My, old uplands ~20-30 m/My, nonlinear rise with relief (von Blanckenburg 2005; Bierman and Caffee 2002; Matmon et al. 2003; Montgomery and Brandon 2002), not checked; declared before any run (lane SR PREREG.txt)"));
            assumptions.push(pa("A-BASE-LEVEL-TAPER-EXPONENT", "exponent p of the base-level erosion taper", rp.taper_exponent, "1", "project assumption awaiting review; declared before any run (lane SR PREREG.txt)"));
            assumptions.push(pa("A-ALLUVIAL-GRADE-SLOPE", "material passing a continental land cell is stored on it up to the graded level receiver + S_a L (compensated thickness), before lakes, flooded continental cells and the legacy placement", rp.alluvial_grade_slope, "1", "project assumption awaiting review; recalled lowland alluvial gradients ~2e-5 to ~2e-4, not checked; declared before any run"));
            assumptions.push(pa("A-CONTINENTAL-SHELF-WAVE-BASE", "flooded continental cells at a river mouth are filled to sea level minus this depth (compensated, nearest first) before any sediment reaches oceanic cells", rp.shelf_wave_base_depth_m, "m", "lane SD's shelf wave base value, unchanged; Earth shelf break ~130 m recalled, not checked"));
        }
        // SEED-UPLANDS (lane UP): the threshold hillslope law, if selected.
        if let Some((cs, _)) = &critical_slope {
            stream_power_scheme = format!("{stream_power_scheme}+critical-slope-envelope-v1");
            approximations.push("critical-slope-envelope-threshold-hillslopes-v1".into());
            let pa = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
                id: id.into(),
                statement: statement.into(),
                value,
                units: units.into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: note.into(),
            };
            assumptions.push(pa("A-HILLSLOPE-CRITICAL-SLOPE", "critical slope S_c of cell means: the envelope z_env(i) = min_j (b_j + S_c path(i, j)) (b = land mean or water fill level, 8-neighbour paths) caps how far a cell mean may stand above the lowest ground around it; land cells above it erode at (z - z_env)/tau (crust removal; isostasy returns most of it), the cell-mean target being max(stream-power target, this)", cs.critical_slope, "1", "project assumption awaiting review; recalled range-front geometry (Burbank et al. 1996; Montgomery and Brandon 2002), not checked; Earth p99.5 steepest cell-mean gradient at 0.625 deg is 0.0375 (ETOPO 2022, lane UP); declared before any run"));
            assumptions.push(pa("A-KNICKPOINT-CELERITY", "celerity at which an over-steepened front is cut back; threshold response time tau = sqrt(cell area) / celerity, backward Euler with the envelope held for the step; crust removal sized for the post-isostasy surface, d = (z - z_env)/a dt/(tau + dt), a = (rho_m - rho_c)/rho_m (rule 2)", cs.knickpoint_celerity_m_per_my, "m/My", "project assumption awaiting review; recalled knickpoint retreat 1 to >100 km/My (Crosby and Whipple 2006; Whittaker and Boulton 2012), not checked; declared before any run"));
        }
        // CLIMATE-GENESIS (lane CG): `epoch-climate-runoff-v1`, if selected.
        if let Some(cp) = &climate_law {
            stream_power_scheme = format!("{stream_power_scheme}+epoch-climate-runoff-v1");
            for (id, _) in pf_surface::climate::seed::genesis::GENESIS_CLIMATE_APPROXIMATIONS_V1 {
                approximations.push((*id).into());
            }
            if cp.subgrid_cell_mean == ClimateSubgridCellMeanV1::CellRunoffSubgridErosivityV1 {
                approximations.push(pf_surface::climate::seed::genesis::GENESIS_CLIMATE_SUBGRID_APPROXIMATION_V1.0.into());
            }
            let pa = |id: &str, statement: &str, value: f64, units: &str, note: &str| AssumptionV0 {
                id: id.into(),
                statement: statement.into(),
                value,
                units: units.into(),
                source: AssumptionSourceV0::ProjectAssumption,
                note: note.into(),
            };
            assumptions.push(pa("A-CLIMATE-GENESIS-R-REF", "stream power on discharge: every fluvial kernel reads A_eff = sum over the catchment of a R / r_ref; K is the erodibility at the reference runoff r_ref (surface.uniform_runoff_m_per_yr)", cp.r_ref_m_per_yr, "m/yr", "project assumption awaiting review (lane CG DESIGN.md 4.1-4.2); never stacked with a uniform-runoff-equivalent K (validation refuses K < 1.5e-6 /yr)"));
            assumptions.push(pa("A-CLIMATE-GENESIS-EPOCH", "a full climate solve every epoch on the coarsened grid, frozen-column local rainout between solves", cp.epoch_steps as f64 * dt, "My", "project assumption awaiting review; S0 timing on C6 (lane CG)"));
            // HYPSO-LEVERS (opt-in): the erosive-runoff response.
            match cp.erosive_runoff.as_ref() {
                Some(pf_contracts::laws::ErosiveRunoffLawV1::SublinearAboveReferenceV1 { erosive_runoff_exponent }) => {
                    stream_power_scheme = format!("{stream_power_scheme}+erosive-runoff-sublinear-above-reference-v1");
                    assumptions.push(pa("A-EROSIVE-RUNOFF-EXPONENT", "the fluvial kernels read the eroding runoff R_e = R for R <= r_ref and r_ref (R / r_ref)^p above (sublinear response of long-term erosion to mean runoff); the lake balance and published fields keep R", erosive_runoff_exponent.value, "1", "project assumption awaiting review (lane HYPSO-LEVERS 2026-10-21); stochastic-threshold literature recalled, not checked"));
                }
                Some(pf_contracts::laws::ErosiveRunoffLawV1::CapV1 { erosive_runoff_cap_m_per_yr }) => {
                    stream_power_scheme = format!("{stream_power_scheme}+erosive-runoff-cap-v1");
                    assumptions.push(pa("A-EROSIVE-RUNOFF-CAP", "the fluvial kernels read the eroding runoff R_e = min(R, cap); the lake balance and published fields keep R", erosive_runoff_cap_m_per_yr.value, "m/yr", "project assumption awaiting review (lane HYPSO-LEVERS 2026-10-21); recalled, not checked"));
                }
                None => {}
            }
            assumptions.push(pa("A-CLIMATE-GENESIS-RUNOFF-FLOOR", "land runoff R = max(P - E, f P); lakes max(P - PET, 0); published hydrology and water-exchange discharge stay uniform-runoff in v1", cp.floor_fraction, "1", "lane CLIMATE-RUNOFF's land_p_minus_e_floor_v1 constant; integrator ruling 2026-10-20 (published T0 discharge stays uniform in v1)"));
        }
        if let Some(evidence) = &params.native_evidence {
            if let Some(fines) = &evidence.marine_fines {
                let (fraction, depth) = fines.law.parameters();
                assumptions.retain(|a| a.id != "A-NO-DEEP-SEA-BYPASS");
                approximations.push(pf_contracts::laws::MarineFinesLawV1::ID.into());
                assumptions.push(AssumptionV0 {
                    id: "A-NATIVE-TERRIGENOUS-FINES".into(),
                    statement: "selected share of finite river-mouth supply mixes instantly within each connected ocean and settles on its deep floor".into(),
                    value: fraction, units: "1".into(), source: AssumptionSourceV0::ProjectAssumption,
                    note: format!("floor depth threshold {depth} m; actual debits/credits in native-producer-evidence; no currents, pelagic production or new mass"),
                });
            }
            if evidence.formation_age.is_some() {
                approximations.push("native-formation-age-v1-basal-ocean-cohorts".into());
            }
            if let Some(a) = assumptions.iter_mut().find(|a| a.id == "A-NO-PARCEL-LEDGER") {
                a.note = "native-producer-evidence publishes basal origin and/or fines-transfer observations only; no finite ocean-parcel, excavated layer/body association or canonical admission".into();
                a.statement = "published native evidence has explicitly limited coverage; a rock's full material lineage remains unassigned".into();
            }
        }
        let terrain = TerrainDomainV0 {
            t_my: state.t_my,
            sea_level_m: sea,
            approximations,
            assumptions,
            plates: plates.clone(),
            boundaries,
            codes: TerrainCodesV0::v0(SEDIMENTARY_COVER_THRESHOLD_M),
            uplift_prior: uplift_record,
            erosion: ErosionRecordV0 {
                stream_power_k_per_my: ep.k_per_my,
                stream_power_m: ep.m,
                stream_power_n: ep.n,
                stream_power_scheme,
                hillslope_diffusivity_m2_per_my: ep.diffusivity_m2_per_my,
                diffusion_stability_limit: ep.stability_limit,
                diffusion_number_max_observed: nu_max_observed,
                diffusion_substeps_max: substeps_max,
                rock_density_kg_m3: rho_c,
            },
            sediment_balance: sediment_ledger.balance()?,
            crust_volume: crust,
            summary: TerrainSummaryV0 {
                steps,
                step_my: dt,
                min_elevation_m: mn,
                max_elevation_m: mx,
                land_area_fraction: land / total_area,
                mean_crust_thickness_m: thick_sum / total_area,
                max_uplift_rate_m_per_my: max_uplift.max(0.0),
                max_erosion_rate_m_per_my: max_erosion_rate,
            },
        };
        terrain
            .validate()
            .map_err(|e| PreviewError::Domain(e.to_string()))?;

        // --- resource potential (interpretation of conditions, spec 2.3) --------
        let (mut resource, mut potentials) = resource_potential(&state, &elevation, &routing, &areas, sea)?;
        // SEED-MAGMATISM (lane MG): the additive magmatic / hydrothermal
        // families and the derived volcano record (law selected only).
        let mut volcanoes: Option<pf_contracts::terrain::VolcanoFeaturesDomainV0> = None;
        if let (Some(me), Some(ms)) = (
            magma.as_ref(),
            orogen.as_ref().and_then(|o| o.conserve.as_ref()).and_then(|c| c.magma.as_ref()),
        ) {
            for f in crate::magmatism::magmatic_potentials(&state, &elevation, sea, &areas, ms)
                .map_err(PreviewError::Domain)?
            {
                let kind = ResourceFamilyKindV0::MAGMATIC_V1
                    .iter()
                    .copied()
                    .find(|k| k.field_name() == f.field)
                    .ok_or_else(|| PreviewError::Domain(format!("magmatism: unknown family {}", f.field)))?;
                resource.families.push(ResourceFamilyV0 {
                    kind,
                    field: f.field.into(),
                    basis: f.basis.iter().map(|s| s.to_string()).collect(),
                    formula: f.formula,
                    mean_potential: f.mean_potential,
                });
                potentials.push((kind, f.values));
            }
            resource.assumptions.extend(crate::magmatism::magmatic_assumptions_for(magma.as_ref().and_then(|m| m.params.v2.as_ref())));
            resource
                .validate()
                .map_err(|e| PreviewError::Domain(e.to_string()))?;
            volcanoes = Some(
                crate::magmatism::volcano_features(
                    me,
                    &grid,
                    &state,
                    &elevation,
                    sea,
                    ms,
                    plumes.as_ref().map(|p| p.0.as_slice()),
                )
                .map_err(PreviewError::Domain)?,
            );
        }

        // SEED-STRATIGRAPHY (lane ST): the rock column record (law selected
        // only).
        let mut stratigraphy: Option<pf_contracts::terrain::StratigraphyDomainV0> = None;
        if let (Some(sp), Some(cs)) = (strata.as_ref(), orogen.as_ref().and_then(|o| o.conserve.as_ref())) {
            if let Some(stacks) = cs.parcels.strata.as_ref() {
                let dom = pf_geodynamics::conserve::dominant_parcels(&grid, &cs.parcels, &state)?;
                let mut d = pf_geodynamics::stratigraphy::publish(
                    &grid,
                    state.t_my,
                    stacks,
                    &dom,
                    &state.sediment_m,
                    &state.thickness_m,
                    &state.birth_t_my,
                    sp,
                    sp.publish_depth_m,
                    strata_assumptions(sp),
                );
                if sp.v2.is_some() {
                    strata_uplift_into(&mut d, &grid, &cs.parcels, &state)?;
                }
                d.validate().map_err(|e| PreviewError::Domain(e.to_string()))?;
                stratigraphy = Some(d);
            }
        }

        Ok(PreviewOutputs {
            native_evidence: native_product,
            exchange_epochs,
            grid,
            params,
            seed_path,
            plates,
            state,
            elevation_m: elevation,
            uplift_rate_m_per_my: uplift_rate,
            routing,
            terrain,
            hydrology,
            cell_reach,
            resource,
            potentials,
            steps: logs,
            elevation_t0_m: elevation_t0,
            sea_level_t0,
            continental_bodies: bodies_report,
            volcanoes,
            stratigraphy,
            orogen,
        })
    }
}

/// Area-weighted volume of a per-cell thickness (m^3), summed in cell
/// order (the same reduction `CrustState::crust_volume_m3` uses).
fn volume_m3(thickness_m: &[f64], areas: &[f64]) -> f64 {
    thickness_m.iter().zip(areas).map(|(h, a)| h * a).sum()
}

fn kind_rank(k: BoundaryKindV0) -> u8 {
    match k {
        BoundaryKindV0::Ridge => 0,
        BoundaryKindV0::Subduction => 1,
        BoundaryKindV0::Collision => 2,
        BoundaryKindV0::Transform => 3,
    }
}

fn sat(x: f64, scale: f64) -> f64 {
    (x / scale).clamp(0.0, 1.0)
}

/// Resource potential families from the final state. Each formula is
/// recorded in the artifact next to the fields it reads.
fn resource_potential(
    state: &CrustState,
    elevation: &[f64],
    routing: &Routing,
    areas: &[f64],
    sea: f64,
) -> Result<
    (
        ResourcePotentialDomainV0,
        Vec<(ResourceFamilyKindV0, Vec<f32>)>,
    ),
    PreviewError,
> {
    let n = state.cells();
    let total_area: f64 = areas.iter().sum();
    const ARC_SAT_M: f64 = 2000.0;
    const SHORTENING_SAT_M: f64 = 10_000.0;
    const RIDGE_AGE_SCALE_MY: f64 = 10.0;
    const CRATON_AGE_MY: f64 = 2000.0;
    const PLACER_SED_M: f64 = 100.0;
    const PLACER_AREA_M2: f64 = 1.0e11;
    const BASIN_SED_M: f64 = 1000.0;
    let mut out: Vec<(ResourceFamilyKindV0, Vec<f32>)> = Vec::new();
    let mut families = Vec::new();
    for kind in ResourceFamilyKindV0::ALL {
        let mut field = Vec::with_capacity(n);
        let (basis, formula): (Vec<&str>, String) = match kind {
            ResourceFamilyKindV0::ConvergentMagmatic => (
                vec!["arc_added (subduction history)", "cumulative_erosion_m"],
                format!("sat(arc_added/{ARC_SAT_M}) * (0.25 + 0.75 * exposure)"),
            ),
            ResourceFamilyKindV0::DivergentHydrothermal => (
                vec!["crust_kind", "crust_age_my"],
                format!("oceanic ? exp(-age/{RIDGE_AGE_SCALE_MY}) : 0"),
            ),
            ResourceFamilyKindV0::CollisionalMetamorphic => (
                vec!["shortening_added (collision history)", "cumulative_erosion_m"],
                format!("sat(shortening/{SHORTENING_SAT_M}) * (0.25 + 0.75 * exposure)"),
            ),
            ResourceFamilyKindV0::CratonicShield => (
                vec!["crust_kind", "crust_age_my", "arc_added", "shortening_added"],
                format!("continental ? sat(age/{CRATON_AGE_MY}) * (1 - sat(arc/{ARC_SAT_M})) * (1 - sat(shortening/{SHORTENING_SAT_M})) : 0"),
            ),
            ResourceFamilyKindV0::AlluvialPlacer => (
                vec!["sediment_thickness_m", "contributing_area_m2", "elevation_m"],
                format!("land ? sat(sediment/{PLACER_SED_M}) * sat(area/{PLACER_AREA_M2}) : 0"),
            ),
            // SEED-MAGMATISM: the additive families are never in ALL.
            ResourceFamilyKindV0::PorphyryCuMoAu | ResourceFamilyKindV0::EpithermalAuAg | ResourceFamilyKindV0::VolcanogenicMassiveSulfide | ResourceFamilyKindV0::MagmaticNiCuPge | ResourceFamilyKindV0::KimberliteDiamond | ResourceFamilyKindV0::PegmatiteCrystal | ResourceFamilyKindV0::GraniteSnW => {
                unreachable!("magmatic families are appended by magmatism.rs (basis)")
            }
            ResourceFamilyKindV0::BasinSedimentary => (
                vec!["sediment_thickness_m"],
                format!("sat(sediment/{BASIN_SED_M})"),
            ),
        };
        let mut weighted = 0.0;
        for c in 0..n {
            let exposure = sat(state.eroded_total_m[c], EXPOSURE_SATURATION_M);
            let oceanic = state.kind[c] == CrustKindCodeV0::Oceanic;
            let land = elevation[c] >= sea;
            let v = match kind {
                ResourceFamilyKindV0::ConvergentMagmatic => {
                    sat(state.arc_added_m[c], ARC_SAT_M) * (0.25 + 0.75 * exposure)
                }
                ResourceFamilyKindV0::DivergentHydrothermal => {
                    if oceanic {
                        (-state.age_my(c).max(0.0) / RIDGE_AGE_SCALE_MY).exp()
                    } else {
                        0.0
                    }
                }
                ResourceFamilyKindV0::CollisionalMetamorphic => {
                    sat(state.shortening_added_m[c], SHORTENING_SAT_M) * (0.25 + 0.75 * exposure)
                }
                ResourceFamilyKindV0::CratonicShield => {
                    if oceanic {
                        0.0
                    } else {
                        sat(state.age_my(c), CRATON_AGE_MY)
                            * (1.0 - sat(state.arc_added_m[c], ARC_SAT_M))
                            * (1.0 - sat(state.shortening_added_m[c], SHORTENING_SAT_M))
                    }
                }
                ResourceFamilyKindV0::AlluvialPlacer => {
                    if land {
                        sat(state.sediment_m[c], PLACER_SED_M)
                            * sat(routing.contributing_area_m2[c], PLACER_AREA_M2)
                    } else {
                        0.0
                    }
                }
                ResourceFamilyKindV0::BasinSedimentary => sat(state.sediment_m[c], BASIN_SED_M),
                ResourceFamilyKindV0::PorphyryCuMoAu | ResourceFamilyKindV0::EpithermalAuAg | ResourceFamilyKindV0::VolcanogenicMassiveSulfide | ResourceFamilyKindV0::MagmaticNiCuPge | ResourceFamilyKindV0::KimberliteDiamond | ResourceFamilyKindV0::PegmatiteCrystal | ResourceFamilyKindV0::GraniteSnW => {
                    unreachable!("magmatic families are appended by magmatism.rs (value)")
                }
            };
            let v = v.clamp(0.0, 1.0);
            weighted += v * areas[c];
            field.push(v as f32);
        }
        families.push(ResourceFamilyV0 {
            kind,
            field: kind.field_name().into(),
            basis: basis.into_iter().map(String::from).collect(),
            formula,
            mean_potential: weighted / total_area,
        });
        out.push((kind, field));
    }
    let pa = |id: &str, statement: &str, value: f64, units: &str| AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source: AssumptionSourceV0::ProjectAssumption,
        note: "coarse potential scale; model-assumption-review Q6".into(),
    };
    let resource = ResourcePotentialDomainV0 {
        t_my: state.t_my,
        usage_limits: RESOURCE_PREVIEW_USAGE_LIMITS.into(),
        families,
        exposure_saturation_m: EXPOSURE_SATURATION_M,
        assumptions: vec![
            pa(
                "A-RES-EXPOSURE",
                "cumulative bedrock erosion at which exposure saturates",
                EXPOSURE_SATURATION_M,
                "m",
            ),
            pa(
                "A-RES-ARC-SAT",
                "arc addition at which convergent-magmatic potential saturates",
                ARC_SAT_M,
                "m",
            ),
            pa(
                "A-RES-SHORTENING-SAT",
                "collisional thickening at which metamorphic potential saturates",
                SHORTENING_SAT_M,
                "m",
            ),
            pa(
                "A-RES-RIDGE-AGE",
                "e-folding crust age of divergent-hydrothermal potential",
                RIDGE_AGE_SCALE_MY,
                "My",
            ),
            pa(
                "A-RES-CRATON-AGE",
                "crust age at which cratonic potential saturates",
                CRATON_AGE_MY,
                "My",
            ),
            pa(
                "A-RES-PLACER",
                "sediment thickness and contributing area scales for alluvial potential",
                PLACER_SED_M,
                "m",
            ),
            pa(
                "A-RES-BASIN",
                "sediment thickness at which basin potential saturates",
                BASIN_SED_M,
                "m",
            ),
        ],
    };
    resource
        .validate()
        .map_err(|e| PreviewError::Domain(e.to_string()))?;
    Ok((resource, out))
}

/// Digests of the three published artifacts.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct PublishedDigests {
    pub terrain: String,
    pub hydrology: String,
    pub resource: String,
}

fn backend() -> BackendId {
    BackendId {
        name: BACKEND_NAME.into(),
        tier: Tier::Preview,
        version: BACKEND_VERSION.into(),
    }
}

fn rock_family(state: &CrustState, c: usize) -> u8 {
    if state.sediment_m[c] >= SEDIMENTARY_COVER_THRESHOLD_M {
        RockFamilyCodeV0::SedimentaryCover.code()
    } else if state.shortening_added_m[c] > 0.0
        && state.shortening_added_m[c] >= state.arc_added_m[c]
    {
        RockFamilyCodeV0::CollisionalThickened.code()
    } else if state.arc_added_m[c] > 0.0 {
        RockFamilyCodeV0::ArcMagmatic.code()
    } else if state.kind[c] == CrustKindCodeV0::Oceanic {
        RockFamilyCodeV0::OceanicIgneous.code()
    } else {
        RockFamilyCodeV0::ContinentalUndifferentiated.code()
    }
}

/// Publish terrain, hydrology and resource artifacts under `out_root`
/// (artifact ids "terrain", "hydrology", "resource"). Refuses existing
/// targets (publish is atomic and never overwrites).
pub fn publish(
    out: &PreviewOutputs,
    out_root: &Path,
    started: Instant,
) -> Result<PublishedDigests, PreviewError> {
    publish_with_exchange_sidecars(out, out_root, started, None)
}

/// Opt-in caller-provided complete exchange serialization (e.g. bounded chunks).
/// None preserves the legacy single-sidecar path and its explicit size refusal.
pub fn publish_with_exchange_sidecars(
    out: &PreviewOutputs,
    out_root: &Path,
    started: Instant,
    exchange_sidecars: Option<Vec<(pf_contracts::artifact::SidecarDescriptor, Vec<u8>)>>,
) -> Result<PublishedDigests, PreviewError> {
    let grid = &out.grid;
    let n = grid.cells();
    let shape = [grid.nlat() as u64, grid.nlon() as u64];
    let domain = grid.domain_declaration(SURFACE_DOMAIN)?;
    let params_digest = digest_params(&out.params)?;
    let constraints_digest = digest_canonical(&serde_json::json!({"constraints": []}))?;
    let st = &out.state;
    let f32v = |v: &[f64]| -> Vec<f32> { v.iter().map(|&x| x as f32).collect() };

    // Terrain fields.
    let age: Vec<f64> = (0..n).map(|c| st.age_my(c)).collect();
    let plate: Vec<u32> = st.plate.iter().map(|p| p.0).collect();
    let crust_kind: Vec<u8> = st.kind.iter().map(|k| k.code()).collect();
    let rock: Vec<u8> = (0..n).map(|c| rock_family(st, c)).collect();
    let terrain_fields: Vec<FieldBlob> = vec![
        f32_field(
            "elevation_m",
            SURFACE_DOMAIN,
            &shape,
            units::M,
            &f32v(&out.elevation_m),
        )?,
        f32_field(
            "crust_thickness_m",
            SURFACE_DOMAIN,
            &shape,
            units::M,
            &f32v(&st.thickness_m),
        )?,
        f32_field(
            "crust_age_my",
            SURFACE_DOMAIN,
            &shape,
            units::MY,
            &f32v(&age),
        )?,
        u32_field("plate_id", SURFACE_DOMAIN, &shape, units::ID, &plate)?,
        u8_field(
            "boundary_class",
            SURFACE_DOMAIN,
            &shape,
            units::CODE,
            &st.boundary_class,
        )?,
        f32_field(
            "uplift_rate_m_per_my",
            SURFACE_DOMAIN,
            &shape,
            units::M_PER_MY,
            &f32v(&out.uplift_rate_m_per_my),
        )?,
        f32_field(
            "sediment_thickness_m",
            SURFACE_DOMAIN,
            &shape,
            units::M,
            &f32v(&st.sediment_m),
        )?,
        f32_field(
            "cumulative_erosion_m",
            SURFACE_DOMAIN,
            &shape,
            units::M,
            &f32v(&st.eroded_total_m),
        )?,
        f32_field(
            "normal_rate_m_per_my",
            SURFACE_DOMAIN,
            &shape,
            units::M_PER_MY,
            &f32v(&st.normal_rate_m_per_my),
        )?,
        u8_field(
            "crust_kind",
            SURFACE_DOMAIN,
            &shape,
            units::CODE,
            &crust_kind,
        )?,
        u8_field("rock_family", SURFACE_DOMAIN, &shape, units::CODE, &rock)?,
    ];
    // SEED-STRATIGRAPHY (lane ST): the rock column record rides as a
    // sidecar of the terrain artifact (none on every other profile:
    // publish_typed is publish_typed_with_sidecars with no sidecar).
    let mut terrain_sidecars = Vec::new();
    if let Some(sidecar) = crate::preview_module::native_evidence_sidecar(out, &terrain_fields)? {
        terrain_sidecars.push(sidecar);
    }
    if let Some(sd) = &out.stratigraphy {
        let bytes = serde_json::to_vec(sd).map_err(|e| PreviewError::Domain(e.to_string()))?;
        terrain_sidecars.push((
            pf_contracts::artifact::SidecarDescriptor {
                logical_name: "stratigraphy".into(),
                path: "derived/stratigraphy.json".into(),
                schema: pf_contracts::terrain::STRATIGRAPHY_SCHEMA_V0.into(),
                byte_length: 0,
                digest: String::new(),
                record_count: Some(sd.lith.len() as u64),
            },
            bytes,
        ));
    }
    let terrain_digest = crate::receipt::publish_typed_with_sidecars(
        out_root,
        ArtifactKind::Terrain,
        "terrain",
        &out.seed_path,
        &params_digest,
        &constraints_digest,
        vec![],
        backend(),
        vec![domain.clone()],
        terrain_fields,
        &out.terrain,
        receipt(started),
        terrain_sidecars,
    )?;

    // Hydrology fields.
    let cf = HydrologyCellFieldsV0::v0();
    let flow_to: Vec<u32> = out
        .routing
        .flow_to
        .iter()
        .map(|t| t.map(|e| e.0 as u32).unwrap_or(cf.terminal_code))
        .collect();
    let basin: Vec<u32> = (0..n)
        .map(|c| {
            out.routing
                .basin_index(pf_contracts::spatial::ElementId(c as u64))
                .map(|i| out.routing.basins[i].basin_id.0 as u32)
                .unwrap_or(cf.no_basin_code)
        })
        .collect();
    let reach: Vec<u32> = out
        .cell_reach
        .iter()
        .map(|r| r.map(|x| x.0 as u32).unwrap_or(cf.no_reach_code))
        .collect();
    let hydrology_fields: Vec<FieldBlob> = vec![
        u32_field(
            "water_body_id",
            SURFACE_DOMAIN,
            &shape,
            units::ID,
            &out.routing.water_body_id,
        )?,
        u32_field(
            &cf.flow_to_field,
            SURFACE_DOMAIN,
            &shape,
            units::ID,
            &flow_to,
        )?,
        u32_field(&cf.basin_field, SURFACE_DOMAIN, &shape, units::ID, &basin)?,
        u32_field(&cf.reach_field, SURFACE_DOMAIN, &shape, units::ID, &reach)?,
        f32_field(
            &cf.fill_level_field,
            SURFACE_DOMAIN,
            &shape,
            units::M,
            &f32v(&out.routing.fill_level_m),
        )?,
        f32_field(
            &cf.contributing_area_field,
            SURFACE_DOMAIN,
            &shape,
            units::M2,
            &f32v(&out.routing.contributing_area_m2),
        )?,
        f32_field(
            &cf.discharge_field,
            SURFACE_DOMAIN,
            &shape,
            units::M3_S,
            &f32v(&out.routing.discharge_proxy_m3_s),
        )?,
    ];
    let exchange_sidecars = if let Some(sidecars) = exchange_sidecars {
        sidecars
    } else {
        let exchange_bytes = serde_json::to_vec(&out.exchange_epochs)
            .map_err(|e| ContractError::Serialization("exchange.serialize", e))?;
        if exchange_bytes.len() > 128 * 1024 * 1024 {
            return Err(PreviewError::Domain(
                "exchange sidecar exceeds 128 MiB; chunk explicitly before publication".into(),
            ));
        }
        vec![(pf_contracts::artifact::SidecarDescriptor {
            logical_name: "surface-exchange".into(),
            path: "history/surface-exchange.json".into(),
            schema: "pf.surface-exchange/v1".into(),
            byte_length: 0, digest: String::new(),
            record_count: Some(out.exchange_epochs.len() as u64),
        }, exchange_bytes)]
    };
    let hydrology_digest = crate::receipt::publish_typed_with_sidecars(
        out_root,
        ArtifactKind::Hydrology,
        "hydrology",
        &out.seed_path,
        &params_digest,
        &constraints_digest,
        vec![terrain_digest.clone()],
        backend(),
        vec![domain.clone()],
        hydrology_fields,
        &out.hydrology,
        receipt(started),
        exchange_sidecars,
    )?;

    // Resource fields.
    let mut resource_fields = Vec::new();
    for (kind, data) in &out.potentials {
        resource_fields.push(f32_field(
            kind.field_name(),
            SURFACE_DOMAIN,
            &shape,
            units::FRACTION,
            data,
        )?);
    }
    // SEED-MAGMATISM (lane MG): the derived volcano record rides as a sidecar
    // of the resource artifact (none on every other profile: publish_typed is
    // publish_typed_with_sidecars with no sidecar).
    let mut resource_sidecars = Vec::new();
    if let Some(v) = &out.volcanoes {
        let bytes = serde_json::to_vec(v).map_err(|e| PreviewError::Domain(e.to_string()))?;
        resource_sidecars.push((
            pf_contracts::artifact::SidecarDescriptor {
                logical_name: "volcano-features".into(),
                path: "derived/volcano-features.json".into(),
                schema: "pf.volcano-features/v0".into(),
                byte_length: 0,
                digest: String::new(),
                record_count: Some(v.features.len() as u64),
            },
            bytes,
        ));
    }
    let resource_digest = crate::receipt::publish_typed_with_sidecars(

        out_root,
        ArtifactKind::Resource,
        "resource",
        &out.seed_path,
        &params_digest,
        &constraints_digest,
        vec![terrain_digest.clone(), hydrology_digest.clone()],
        backend(),
        vec![domain],
        resource_fields,
        &out.resource,
        receipt(started),
        resource_sidecars,
    )?;

    Ok(PublishedDigests {
        terrain: terrain_digest,
        hydrology: hydrology_digest,
        resource: resource_digest,
    })
}

// ---------------------------------------------------------------------------
// CLIMATE-V0 (lane CL, 2026-10-12): OPT-IN separate entry points.
//
// Neither `simulate` nor `publish` calls anything below, so default output
// bytes are unchanged. The climate tier is a pure function of PUBLISHED
// fields (pf-surface `climate.rs`): f32 elevation, the terrain sea-level
// datum, the hydrology fill level and the marine water body of the
// hydrology routing. Both entry
// points read exactly those values, so they agree bit for bit (pinned by
// tests/climate_v0_preview.rs).

/// Artifact id of the published climate tier.
pub const CLIMATE_ARTIFACT_ID: &str = "climate";
/// Backend version string of the `climate-coarse-v0` producer.
pub const CLIMATE_BACKEND_VERSION: &str = "0.1.0";

/// The climate tier from in-memory preview outputs, reading what `publish`
/// writes: elevation and the routing fill level rounded to f32,
/// `terrain.sea_level_m`, and the routing ocean (the cells of the marine
/// water body).
pub fn climate_v0_from_outputs(
    out: &PreviewOutputs,
    recipe: &pf_contracts::laws::ClimateRecipeV1,
) -> Result<pf_surface::climate::ClimateOutputV0, PreviewError> {
    let elevation: Vec<f64> = out.elevation_m.iter().map(|&z| z as f32 as f64).collect();
    let fill: Vec<f64> = out.routing.fill_level_m.iter().map(|&z| z as f32 as f64).collect();
    let inputs = pf_surface::climate::ClimateInputs {
        grid: &out.grid,
        elevation_m: &elevation,
        fill_level_m: &fill,
        sea_level_m: out.terrain.sea_level_m,
        ocean: &out.routing.ocean,
    };
    // SEED-CLIMATE (lane CM request): dispatch on the recipe's law;
    // climate-coarse-v0 still runs derive_climate_coarse_v0 unchanged.
    pf_surface::climate::derive_climate(&inputs, recipe)
        .map_err(|e| PreviewError::Domain(format!("climate: {e}")))
}

/// What `publish_climate_v0` produced.
pub struct ClimatePublished {
    pub digest: String,
    pub terrain_digest: String,
    pub hydrology_digest: String,
    pub output: pf_surface::climate::ClimateOutputV0,
}

/// Published-run entry point: verify `run_root` (terrain, hydrology,
/// resource; `region::SourceRun`), derive the climate tier from the
/// published fields, and publish `<out_root>/climate` (`pf.climate/v0`)
/// citing the terrain and hydrology digests. Refuses an existing target.
pub fn publish_climate_v0(
    run_root: &Path,
    out_root: &Path,
    recipe: &pf_contracts::laws::ClimateRecipeV1,
    started: Instant,
) -> Result<ClimatePublished, PreviewError> {
    use pf_contracts::artifact::{
        publish_artifact, read_artifact_header, ArtifactHeader, SemanticManifest,
    };
    use pf_contracts::domain::WaterBodyKindV0;
    use pf_surface::climate::{
        ClimateDomainV0, ClimateInputRecordV0, ClimateInputs, CLIMATE_FIELDS_V0,
        CLIMATE_SCHEMA_V0,
    };
    let src = crate::region::SourceRun::load(run_root)
        .map_err(|e| PreviewError::Domain(format!("climate source run: {e}")))?;
    let t = read_artifact_header::<TerrainDomainV0>(&run_root.join("terrain"))?;
    if t.digest != src.terrain_digest {
        return Err(PreviewError::Domain(
            "terrain header changed after verification".into(),
        ));
    }
    let grid = &src.source.grid;
    let n = grid.cells();
    let elevation: Vec<f64> = src
        .fields_f32
        .get("elevation_m")
        .ok_or_else(|| PreviewError::Domain("terrain publishes no elevation_m".into()))?
        .iter()
        .map(|&z| z as f64)
        .collect();
    let fill: Vec<f64> = src
        .fields_f32
        .get("fill_level_m")
        .ok_or_else(|| PreviewError::Domain("hydrology publishes no fill_level_m".into()))?
        .iter()
        .map(|&z| z as f64)
        .collect();
    let body_id = src
        .fields_u32
        .get("water_body_id")
        .ok_or_else(|| PreviewError::Domain("hydrology publishes no water_body_id".into()))?;
    let marine: BTreeMap<u32, bool> = src
        .hydrology
        .water_bodies
        .iter()
        .map(|b| (b.body_id, b.kind == WaterBodyKindV0::Marine))
        .collect();
    let mut ocean = vec![false; n];
    for c in 0..n {
        if body_id[c] != pf_surface::water::DRY {
            ocean[c] = *marine.get(&body_id[c]).ok_or_else(|| {
                PreviewError::Domain(format!(
                    "cell {c} names water body {} that the hydrology payload does not declare",
                    body_id[c]
                ))
            })?;
        }
    }
    let inputs = ClimateInputs {
        grid,
        elevation_m: &elevation,
        fill_level_m: &fill,
        sea_level_m: src.terrain.sea_level_m,
        ocean: &ocean,
    };
    // SEED-CLIMATE (lane CM request): dispatch on the recipe's law.
    let output = pf_surface::climate::derive_climate(&inputs, recipe)
        .map_err(|e| PreviewError::Domain(format!("climate: {e}")))?;
    let domain = ClimateDomainV0::build(
        &output,
        recipe,
        ClimateInputRecordV0 {
            elevation_field: "terrain/elevation_m".into(),
            fill_level_field: "hydrology/fill_level_m".into(),
            sea_level_m: src.terrain.sea_level_m,
            ocean_source: "hydrology/water_body_id cells whose body kind is Marine (the routing ocean: largest connected sub-sea component)".into(),
        },
    )?;
    let shape = [grid.nlat() as u64, grid.nlon() as u64];
    let surface = t
        .header
        .semantic
        .domains
        .first()
        .ok_or_else(|| PreviewError::Domain("terrain declares no domain".into()))?
        .clone();
    let mut descriptors = Vec::new();
    let mut blobs = Vec::new();
    for ((name, values), spec) in output.fields.named().into_iter().zip(CLIMATE_FIELDS_V0) {
        debug_assert_eq!(name, spec.name);
        let data: Vec<f32> = values.iter().map(|&x| x as f32).collect();
        let blob = f32_field(name, &surface.domain_id, &shape, spec.units, &data)?;
        blobs.push((name.to_string(), blob.bytes));
        descriptors.push(blob.descriptor);
    }
    let header = ArtifactHeader {
        schema: CLIMATE_SCHEMA_V0.into(),
        artifact_id: CLIMATE_ARTIFACT_ID.into(),
        semantic: SemanticManifest {
            schema: CLIMATE_SCHEMA_V0.into(),
            artifact_id: CLIMATE_ARTIFACT_ID.into(),
            seed_path: t.header.semantic.seed_path.clone(),
            params_digest: t.header.semantic.params_digest.clone(),
            constraints_digest: t.header.semantic.constraints_digest.clone(),
            inputs: vec![src.terrain_digest.clone(), src.hydrology_digest.clone()],
            backend: BackendId {
                // SEED-CLIMATE (lane CM request): the tier the law selects
                // (climate-coarse-v0 for the v0 recipe, so v0 bytes hold).
                name: pf_surface::climate::climate_tier_of(recipe.law).0.into(),
                tier: Tier::Preview,
                version: CLIMATE_BACKEND_VERSION.into(),
            },
            domains: vec![surface],
            fields: descriptors,
            blob_hashes: Default::default(),
            sidecars: vec![],
        },
        domain,
    };
    let digest = publish_artifact(out_root, header, blobs, vec![], receipt(started))?;
    Ok(ClimatePublished {
        digest,
        terrain_digest: src.terrain_digest,
        hydrology_digest: src.hydrology_digest,
        output,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gdh1_is_continuous_monotone_and_bounded() {
        let a = gdh1_thermal_elevation_m(0.0);
        assert!(
            (a - (GDH1_ASYMPTOTE_M - GDH1_RIDGE_DEPTH_M)).abs() < 1e-9,
            "{a}"
        );
        let left = gdh1_thermal_elevation_m(GDH1_SWITCH_MY - 1e-9);
        let right = gdh1_thermal_elevation_m(GDH1_SWITCH_MY + 1e-9);
        assert!(
            (left - right).abs() < 2.0,
            "GDH1 branches disagree at 20 My: {left} vs {right}"
        );
        let mut prev = a;
        for t in 1..200 {
            let v = gdh1_thermal_elevation_m(t as f64);
            assert!(v <= prev + 1e-9, "not monotone at {t}");
            assert!(v >= 0.0);
            prev = v;
        }
        assert!(
            gdh1_thermal_elevation_m(-5.0) == a,
            "negative age clamps to zero"
        );
    }

    #[test]
    fn preview_profile_forces_airy_and_refuses_latitude_bands() {
        let cfg = PreviewConfig {
            scenario: ScenarioPreset::benchmark_a(),
            fidelity: FidelityPreset::Preview,
            grid: Some((18, 36)),
            seed: 1,
            duration_my: Some(4.0),
            step_my: Some(2.0),
        };
        let p = cfg.profile().unwrap();
        assert_eq!(
            p.mantle.compensation,
            Compensation::Airy,
            "benchmark-a asks for flexure; preview declares Airy"
        );
        // Negative case (F5): a LatitudeBands profile is refused, loudly,
        // with the backend named, never silently run as Uniform.
        let mut bands = p.clone();
        bands.surface.rainfall_mode = RainfallMode::LatitudeBands;
        let err = check_profile(&bands).unwrap_err().to_string();
        assert!(
            err.contains("LatitudeBands") && err.contains("kinematic-preview-v0"),
            "{err}"
        );
        check_profile(&p).unwrap();
    }
}

/// PLATE-POLES (lane PLATE-POLES, 2026-10-21): the opt-in key
/// `orogeny.rate_limited_belt.plate_boundary.rigid_group_contacts` reaches
/// the universe law set (`rigid-group-contacts-v1`), key off keeps every
/// legacy law-set byte, and a short world_v3 run steps with the key on.
#[cfg(test)]
mod rigid_group_law_set_tests {
    use super::*;
    use pf_contracts::laws::RigidGroupLawV1;
    use pf_contracts::params::{FidelityPreset, PlanetParameterProfile, ScenarioPreset};

    fn cfg() -> PreviewConfig {
        PreviewConfig { scenario: ScenarioPreset::toy(), fidelity: FidelityPreset::Preview, grid: Some((36, 72)), seed: 42, duration_my: Some(4.0), step_my: Some(2.0) }
    }

    fn set_key(p: &mut PlanetParameterProfile, on: bool) {
        p.orogeny.as_mut().unwrap().rate_limited_belt.as_mut().unwrap().plate_boundary.as_mut().unwrap().rigid_group_contacts = on;
    }
    fn set_release(p: &mut PlanetParameterProfile, on: bool) {
        p.orogeny.as_mut().unwrap().rate_limited_belt.as_mut().unwrap().plate_boundary.as_mut().unwrap().rigid_group_release = on;
    }

    #[test]
    fn rigid_group_key_enters_the_law_set() {
        let off = cfg().profile_with(&PreviewOptions::world_v3()).unwrap();
        let u_off = off.universe().unwrap();
        assert_eq!(u_off.law_set.rigid_group, RigidGroupLawV1::SharedPoleV0);
        let j_off = serde_json::to_string(&u_off).unwrap();
        assert!(!j_off.contains("rigid_group") && !j_off.contains("rigid-group-contacts"), "legacy law set must not name the slot");
        assert!(!serde_json::to_string(&off).unwrap().contains("rigid_group_contacts"), "key off is skipped on write");
        let mut on = off.clone();
        set_key(&mut on, true);
        let u_on = on.universe().unwrap();
        assert_eq!(u_on.law_set.rigid_group, RigidGroupLawV1::RigidGroupInteriorV1);
        assert!(serde_json::to_string(&u_on).unwrap().contains("rigid-group-interior-v1"));
        assert_ne!(u_off.digest().unwrap(), u_on.digest().unwrap());
        let env = PreviewEnv::new(on.clone(), cfg().grid, cfg().seed_path()).unwrap();
        assert!(env.boundary_transform.map_or(false, |t| t.rigid_group_interior) && !env.rigid_group_release);
        // Release without the interior key is refused; with it, the contacts law.
        let mut rel = off.clone();
        set_release(&mut rel, true);
        assert!(rel.universe().is_err());
        set_key(&mut on, true);
        set_release(&mut on, true);
        let u_rel = on.universe().unwrap();
        assert_eq!(u_rel.law_set.rigid_group, RigidGroupLawV1::RigidGroupContactsV1);
        assert!(serde_json::to_string(&u_rel).unwrap().contains("rigid-group-contacts-v1"));
        let back: pf_contracts::laws::UniverseV1 = serde_json::from_str(&j_off).unwrap();
        assert_eq!(back, u_off);
        let mut pinned_off = off.clone();
        pinned_off.pin_universe().unwrap();
        set_key(&mut pinned_off, true);
        assert!(pinned_off.universe().is_err());
        // The key reaches the classifier; a short run steps.
        let env = PreviewEnv::new(on, cfg().grid, cfg().seed_path()).unwrap();
        assert!(env.boundary_transform.map_or(false, |t| t.rigid_group_interior) && env.rigid_group_release);
        let mut run = PreviewRun::init(env).unwrap();
        while run.steps_taken() < run.steps_scheduled() {
            run.step().unwrap();
        }
    }
}
