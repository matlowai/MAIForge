//! Parameter system: scenario blocks vs fidelity blocks (B1 rev 2 §P).
//!
//! Dial validation lives HERE so CLI and backends share one rule set (A6).
//! Resolved values and declared numerical substeps are recorded in every
//! semantic manifest; invalid combinations fail loudly with field names.

use crate::laws::{
    legacy, ConstantV1, CrustSeedingLawV1, ErosionLawV1, IsostasyLawV1, LawSetV1,
    LandDepositionLawV1, OceanicThermalLawV1, SedimentPlacementLawV1, SeaLevelLawV1,
    SourceClassV1, TectonicThickeningLawV1, PlateCarryLawV1, CollisionSedimentLawV1,
    HotspotVolcanismLawV1, ContinentalRiftLawV1, PlateReorganizationLawV1, RiftFlankLawV1,
    HillslopeThresholdLawV1,
    FluvialIncisionThresholdLawV1, MarineRemainderLayerLawV1,
    SubductionTrenchLawV1, PassiveMarginLawV1, PlateSutureLawV1,
    UniverseConstantsV1, UniverseError, UniverseV1, UnitV1, LAW_SET_SCHEMA_V1, UNIVERSE_CONSTANTS_SCHEMA_V1, UNIVERSE_SCHEMA_V1,
};
use crate::laws::FloodBasaltLawV1; // SEED-FLOODBASALT (lane FB)
use crate::laws::SedimentRetentionLawV1; // SEED-SEDIMENT-RETAIN (lane SR)
use crate::laws::CratonSeedingLawV1; // SEED-IRREGULAR-SEEDS (lane IS)
use crate::laws::OrogenPlateauLawV1; // SEED-RANGES (lane RG)
use crate::laws::PlateBoundaryLawV1; // SEED-TRANSFORMS (lane TF)
use crate::laws::MagmatismLawV1; // SEED-MAGMATISM (lane MG)
use crate::laws::StratigraphyLawV1; // SEED-STRATIGRAPHY (lane ST)
use crate::laws::{ClimateBetweenEpochsV1, ClimateErosivityModeV1, ClimateForcingLawV1, ClimateRecipeV1, ClimateSubgridCellMeanV1, ErosiveRunoffLawV1}; // CLIMATE-GENESIS (lane CG); HYPSO-LEVERS
use crate::laws::ShelfGradingLawV1; // SHELF law-set entry (lane CLIMATE-GENESIS-2)
use crate::laws::RigidGroupLawV1; // PLATE-POLES (lane PLATE-POLES)
use crate::laws::ClosedBasinLawV1; // INLAND-BASINS law-set entry (lane INLAND-BASINS)
use crate::laws::LakeOutletLawV1; // LAKE-OUTLET law-set entry (lane LAKE-OUTLET)
use crate::laws::GlacialErosionLawV1; // GLACIAL-T0 (lane GLACIAL-T0)
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

/// Deterministic seed path: seed + derivation label + salt. Every artifact
/// records it in its semantic manifest (A3).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SeedPath {
    pub seed: u64,
    pub derivation: String,
    pub salt: String,
}

#[derive(Debug, thiserror::Error)]
pub enum ParamError {
    #[error("parameter {field} = {value} outside valid range ({lo}, {hi}): {why}")]
    Range {
        field: String,
        value: String,
        lo: String,
        hi: String,
        why: &'static str,
    },
    #[error("parameter combination invalid: {0}")]
    Combination(&'static str),
    /// WS-CONST: the law set / universe constants refused the profile
    /// (typed: OutOfEnvelope, Unsupported, Inconsistent, PinMismatch, ...).
    #[error("universe: {0}")]
    Universe(#[from] UniverseError),
}

fn check(
    field: &'static str,
    value: f64,
    lo: f64,
    hi: f64,
    why: &'static str,
) -> Result<(), ParamError> {
    if !value.is_finite() || !(lo..=hi).contains(&value) {
        return Err(ParamError::Range {
            field: field.into(),
            value: format!("{value:?}"),
            lo: format!("{lo}"),
            hi: format!("{hi}"),
            why,
        });
    }
    Ok(())
}

/// Planet shape and reference densities. Sea level is the bathymetric datum.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PlanetBlock {
    pub radius_m: f64,
    pub gravity_m_s2: f64,
    pub crust_density_kg_m3: f64,
    pub mantle_density_kg_m3: f64,
    pub sea_level_m: f64,
}

/// Tectonic scenario: the world's history parameters. These are held
/// CONSTANT across fidelity comparisons (review #7).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TectonicsBlock {
    pub plate_count: u32,
    pub rate_scale_deg_per_my: f64,
    pub ocean_fraction: f64,
    pub continental_seeds: u32,
    pub hotspot_count: u32,
    pub duration_my: f64,
    pub step_my: f64,
}

/// Mantle physics switches. Radial layer COUNT is a fidelity parameter
/// (FidelityBlock); model choice and compensation mode live here.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ThermalModel {
    HalfSpace,
    Plate,
}

/// Isostatic compensation: `Airy` XOR `Flexure` — never both (M5, review #5).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Compensation {
    Airy,
    Flexure,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MantleBlock {
    pub thermal_model: ThermalModel,
    pub compensation: Compensation,
    /// Required (and only used) when compensation = Flexure.
    pub effective_elastic_thickness_km: f64,
    /// Heat injected per plume per step, K of layer-mean anomaly (M1 source term).
    pub plume_heat_k_per_step: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum RainfallMode {
    Uniform,
    LatitudeBands,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SurfaceBlock {
    pub rainfall_mode: RainfallMode,
    pub uniform_runoff_m_per_yr: f64,
    /// Stream-power K in 1/yr for area in m² and m=0.5, n=1.
    pub erodibility: f64,
    /// Linear hillslope diffusivity (m^2/yr) for the preview erosion pass
    /// (ADR 032, additive 2026-09-13). Defaulted on read and SKIPPED on
    /// write at the default, so every pre-existing profile keeps its
    /// params_digest (pinned by `toy_preview_params_digest_is_unchanged`).
    #[serde(
        default = "default_hillslope_diffusivity",
        skip_serializing_if = "is_default_hillslope_diffusivity"
    )]
    pub hillslope_diffusivity_m2_per_yr: f64,
}

fn default_hillslope_diffusivity() -> f64 {
    0.1
}

fn is_default_hillslope_diffusivity(v: &f64) -> bool {
    *v == default_hillslope_diffusivity()
}

/// Fidelity: discretization ONLY. Resolution, layers, substeps — never
/// duration or scenario (review #7: fidelity comparisons hold history const).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct FidelityBlock {
    /// 0 => 288x144, +1 doubles each axis. Hard cap 3 (2304x1152).
    pub grid_level: u8,
    /// Mantle radial layers. 1..=24. Storage is 3D from day one (A2.1).
    pub radial_layers: u32,
    /// Declared numerical substeps per tectonic step (explicit algorithm,
    /// recorded; never a silent fallback profile).
    pub substeps: u32,
}

/// World-character intent dials (spec §57): aesthetic bias resolved into
/// causal parameter bundles + candidate-scoring objectives — NEVER direct
/// feature insertion. Waterfall Wonder is the first dial; others (Mountain
/// Grandeur, River Drama, Karst Weirdness, ...) join the same pattern.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorldCharacterBlock {
    /// 0 = natural baseline .. 1 = wonderful/ridiculous waterfall worlds.
    /// A deep-history dial: biases differential-erodibility stratigraphy
    /// (major mechanism), uplift stepping and knickpoint persistence so
    /// the world's HISTORY is more likely to create durable dramatic falls
    /// that rivers then discover (§57.2). Detection is downstream; nothing
    /// inserts falls.
    #[serde(default)]
    pub waterfall_wonder: f64,
}

impl Default for WorldCharacterBlock {
    fn default() -> Self {
        Self {
            waterfall_wonder: 0.0,
        }
    }
}

/// Resolved causal channel weights (deterministic function of the dial;
/// recorded via params_digest). Each channel is consumed only by its owning
/// subsystem and only within geologically appropriate settings — no zebra
/// cake (spec §57.2).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorldCharacterBias {
    /// B1 geodynamics: sharpness of uplift-band edges (1.0 = neutral).
    pub uplift_band_sharpness: f64,
    /// B1 geodynamics: persistence of tectonic uplift steps (1.0 = neutral).
    /// A CAUSE bias — waterfall DETECTION thresholds never depend on the
    /// dial: causal-comparison runs score both worlds with the identical
    /// detector (§57.5; review preflight — no ruler-changing).
    pub uplift_step_persistence: f64,
    /// pf-stratigraphy (parked until §54): how often neighboring beds
    /// differ meaningfully in stream-power erodibility K.
    pub resistance_contrast_frequency: f64,
    /// pf-stratigraphy (parked): K ratio between resistant caps and weak
    /// underliers.
    pub resistance_contrast_magnitude: f64,
    /// pf-stratigraphy (parked): lateral persistence of resistant beds
    /// across valleys.
    pub resistant_bed_persistence: f64,
    /// pf-stratigraphy (parked): thickness bias for resistant beds.
    pub resistant_bed_thickness_bias: f64,
}

impl WorldCharacterBlock {
    /// Deterministic, monotone resolution into causal channel weights.
    /// Neutral at wonder = 0; every output finite and bounded (validated
    /// in tests). The stratigraphy channels stay capped so sedimentary
    /// basins never become alternating steel and pudding.
    pub fn resolve(&self) -> WorldCharacterBias {
        let w = self.waterfall_wonder;
        WorldCharacterBias {
            uplift_band_sharpness: 1.0 + 1.5 * w,
            uplift_step_persistence: 1.0 + w,
            resistance_contrast_frequency: 0.12 + 0.45 * w,
            resistance_contrast_magnitude: 1.1 + 1.4 * w,
            resistant_bed_persistence: 0.35 + 0.45 * w,
            resistant_bed_thickness_bias: 0.35 + 0.40 * w,
        }
    }
}

/// OPT-IN two-column isostasy with a volume-primary ocean (FREEBOARD-SPIKE,
/// lane F, 2026-10-11; WORLD_SIMULATION_EXECPLAN §5). When present:
///
/// - continental columns use `planet.crust_density_kg_m3`, oceanic columns
///   use `oceanic_crust_density_kg_m3`, sediment keeps the legacy
///   assumption (loads at `planet.crust_density_kg_m3`), and every column
///   is balanced against the same mass above `compensation_depth_m` in a
///   mantle of `planet.mantle_density_kg_m3` (no 35 km reference column);
/// - the ocean WATER VOLUME is the input; sea level and the flooded
///   fraction are OUTPUTS solved by filling the water-loaded hypsometry
///   with true spherical cell areas. `tectonics.ocean_fraction` is NOT a
///   target in this mode (it is only compared against, never tuned to).
///
/// Absent (the default) = the legacy single-column Airy model with the
/// preset sea level. The field is skipped on write when absent, so every
/// existing profile keeps its params_digest.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct TwoColumnIsostasyBlock {
    /// Oceanic crust bulk density (kg/m^3).
    pub oceanic_crust_density_kg_m3: f64,
    /// Sea-water density for water loading (kg/m^3).
    pub seawater_density_kg_m3: f64,
    /// Common compensation depth below the mantle free-surface datum (m).
    /// Must lie below the deepest crustal root (checked loudly per cell).
    /// FIXED-DEPTH mode (lane F): used only when
    /// `compensation_depth_margin_m` is None.
    pub compensation_depth_m: f64,
    /// Total liquid ocean water volume (m^3): the primary ocean input.
    pub ocean_water_volume_m3: f64,
    /// DERIVED-DEPTH mode (lane K2, 2026-10-11; design review review, ADR 039
    /// amendment C4/C5): when Some, the pf-lab preview derives the
    /// compensation depth at EVERY isostasy evaluation as (deepest crustal
    /// root over all columns) + this margin, and `compensation_depth_m` is
    /// not used. The root of a column, t_dry - (rho_c h_c + rho_s h_s) /
    /// rho_m, does not depend on the depth, and elevations are independent
    /// of the depth up to rounding (it cancels in the mass balance), so the
    /// margin is a validity bound, not a physical parameter. Skipped on
    /// write when None.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub compensation_depth_margin_m: Option<f64>,
}

impl TwoColumnIsostasyBlock {
    /// Earth-like v0 values (stated, not tuned). Integrator ruling
    /// 2026-10-11: values recalled from memory stay labelled PROJECT
    /// ASSUMPTION AWAITING REVIEW until a citation is checked.
    /// - oceanic crust 2900 kg/m^3: project assumption awaiting review
    ///   (recalled: the value Turcotte and Schubert, Geodynamics, use for
    ///   oceanic crust; Carlson and Raskin 1984, Nature 311:555, ~2860-2890).
    /// - sea water 1030 kg/m^3: same value as pf-lab WATER_DENSITY_KG_M3.
    /// - compensation depth 100 km: the fixed-depth fallback (any depth
    ///   below the deepest root gives the same elevations); not used while
    ///   the margin below is set.
    /// - compensation depth margin 10 km: PROJECT ASSUMPTION; a validity
    ///   bound below the deepest root, with no physical effect on
    ///   elevations (the depth cancels; pinned by
    ///   pf-lab tests/crust_bodies.rs derived-depth known answer).
    /// - ocean volume 1.335e18 m^3: project assumption awaiting review
    ///   (recalled: Charette and Smith 2010, Oceanography 23(2):112,
    ///   1.335e9 km^3 for Earth).
    pub fn earth_v0() -> Self {
        Self {
            oceanic_crust_density_kg_m3: 2900.0,
            seawater_density_kg_m3: 1030.0,
            compensation_depth_m: 100_000.0,
            ocean_water_volume_m3: 1.335e18,
            compensation_depth_margin_m: Some(10_000.0),
        }
    }
}

/// OPT-IN continental crust BODIES carried by plates (CRUST-SPIKE, lane K2,
/// 2026-10-11; SEED_WORLDS_EXECPLAN D2 / SEED-CRUST). When present, the
/// initial crust is NOT "continental plates are wholly continental"; instead
/// each continental plate (the same `tectonics.continental_seeds` plates the
/// legacy path uses) carries one continental body inside it and the rest of
/// the plate is oceanic. Named approximations (pf-geodynamics kinematic.rs
/// `seed_bodies`, uplift.rs `rifted_margin_thickness_m`):
///
/// - `geodesic-disc-body-v0`: a body is the set of its host plate's cells
///   nearest (great-circle) to the plate's own nucleation (Voronoi seed)
///   point, grown until its area reaches its share of the declared
///   continental area. Shares are proportional to host-plate area. No new
///   randomness (the plate stream is reused read-only); shape is the
///   minimum-perimeter initial condition, not a history.
/// - `linear-rifted-margin-taper-v0`: crust thickness rises linearly from
///   the oceanic crust thickness at the body coast to the interior
///   continental thickness over `rifted_margin_width_m` (distance measured
///   from the cell centroid to the nearest coast edge midpoint).
///
/// Absent (the default) = legacy whole-plate continents and byte-identical
/// serialization (skipped on write).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ContinentalBodiesBlock {
    /// Continental crust area as a fraction of the planet (0, 1). `None` =
    /// DERIVED as `1 - tectonics.ocean_fraction` (the scenario's declared
    /// ocean fraction read as the continental crust area it implies).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub continental_area_fraction: Option<f64>,
    /// Width (m) over which a rifted (passive) margin thins from interior
    /// continental thickness to oceanic thickness at the coast. 0 = no
    /// taper (plateau bodies).
    pub rifted_margin_width_m: f64,
    /// SEED-IRREGULAR-SEEDS (lane IS, 2026-10-13): selects `craton_seeding
    /// = irregular-cratons-v1` (irregular old cratons, sutures, rough
    /// basement, two-segment margins) in place of the disc shapes; the
    /// area fraction and margin width above still apply. Skipped on write
    /// when absent (every existing block serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub irregular_cratons: Option<IrregularCratonsBlock>,
}

/// Parameters of `irregular-cratons-v1` (see laws.rs; values fixed before
/// any run in `unpublished-design-record`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct IrregularCratonsBlock {
    pub craton_interior_thickness_m: f64,
    pub craton_rank_size_exponent: f64,
    pub cratons_per_unit_share: f64,
    pub basement_roughness_std_m: f64,
    pub suture_crest_thickening_m: f64,
    pub suture_half_width_m: f64,
    pub proximal_margin_stretch: f64,
}

impl IrregularCratonsBlock {
    /// v1 values fixed BEFORE any run (lane IS PREREG.txt, 2026-10-13
    /// 10:20 EDT).
    pub fn v1() -> Self {
        Self {
            craton_interior_thickness_m: 29_600.0,
            craton_rank_size_exponent: 1.4,
            cratons_per_unit_share: 8.0,
            basement_roughness_std_m: 1_500.0,
            suture_crest_thickening_m: 7_000.0,
            suture_half_width_m: 120_000.0,
            proximal_margin_stretch: 1.15,
        }
    }
}

impl ContinentalBodiesBlock {
    /// v0 values. Coverage derived from the scenario's ocean fraction.
    /// Margin width 200 km: PROJECT ASSUMPTION AWAITING REVIEW, informed by rifted-margin
    /// crustal transects in which crust thins from ~30-35 km to under
    /// 10 km over roughly 100-300 km (necking plus hyperextended domains,
    /// e.g. the Iberia-Newfoundland conjugate margins; Peron-Pinvidic and
    /// Manatschal 2009, Int J Earth Sci 98:1581; citation recalled,
    /// UNVERIFIED).
    pub fn v0() -> Self {
        Self {
            continental_area_fraction: None,
            rifted_margin_width_m: 200_000.0,
            irregular_cratons: None,
        }
    }

    /// SEED-IRREGULAR-SEEDS (lane IS, 2026-10-13): `v0()` plus
    /// `irregular-cratons-v1` at its v1 values.
    pub fn cratons_v1() -> Self {
        Self {
            irregular_cratons: Some(IrregularCratonsBlock::v1()),
            ..Self::v0()
        }
    }

    /// The continental area fraction in force: the declared value, else
    /// `1 - ocean_fraction`.
    pub fn resolved_area_fraction(&self, tectonics: &TectonicsBlock) -> f64 {
        self.continental_area_fraction
            .unwrap_or(1.0 - tectonics.ocean_fraction)
    }
}

/// SEED-OROGEN (lane OG, 2026-10-12): opt-in orogeny and coarse-cell
/// erosion laws. One optional sub-block per law family, each skipped on
/// write when absent, so the block is extensible and `None` everywhere (the
/// default) is the legacy law set with byte-identical serialization.
///
/// - `rate_limited_belt` selects `TectonicThickeningLawV1::RateLimitedBeltV1`
///   (absent: legacy `kinematic-overlap-stacking-v0`).
/// - `subgrid_erosion` selects `ErosionLawV1::StreamPowerSubgridV1`
///   (absent: legacy `stream-power-n1-linear-diffusion-v0`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(deny_unknown_fields)]
pub struct OrogenyBlock {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rate_limited_belt: Option<RateLimitedBeltBlock>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subgrid_erosion: Option<SubgridErosionBlock>,
    /// SEED-RELIEF-LAW (lane RL, 2026-10-12): selects
    /// `ErosionLawV1::StreamPowerSubgridNonlinearV1`. Skipped on write when
    /// absent (every existing block serializes as before). Refused together
    /// with `subgrid_erosion` (one erosion law per world).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nonlinear_subgrid_erosion: Option<NonlinearSubgridErosionBlock>,
    /// SEED-UPLANDS (lane UP, 2026-10-13): selects `hillslope_threshold =
    /// critical-slope-envelope-v1`. Skipped on write when absent. Needs a
    /// subgrid erosion law.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub critical_slope_erosion: Option<CriticalSlopeErosionBlock>,
    /// SEED-SEDIMENT-RETAIN (lane SR, 2026-10-13): selects
    /// `sediment_retention = continental-retention-base-level-v1`. Skipped
    /// on write when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sediment_retention: Option<SedimentRetentionBlock>,
    /// INCISION-THRESHOLD (lane IT, 2026-10-20): selects
    /// `fluvial_incision_threshold = stream-power-excess-steepness-threshold-v1`.
    /// Skipped on write when absent. Needs `nonlinear_subgrid_erosion`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub incision_threshold: Option<IncisionThresholdBlock>,
    /// LAKE-HOLD-2 (lane LAKE-HOLD-2, 2026-10-21): `true` selects
    /// `fluvial_incision_threshold =
    /// stream-power-excess-steepness-threshold-trunk-only-v1` (the threshold
    /// on the trunk channel only, the cell-mean closure at k_c = 0). Needs
    /// `incision_threshold`. Skipped on write when false (byte-identical).
    #[serde(default, skip_serializing_if = "is_false")]
    pub incision_threshold_trunk_only: bool,
    /// ALLUVIAL-TRUNK (lane ALLUVIAL-TRUNK, 2026-10-21): `true` selects
    /// `fluvial_alluvial_channel =
    /// stream-power-threshold-free-alluvial-channels-v1` (no threshold on
    /// trunks and lake outlet channels where k_c / A^(m/n) < S_a). Needs
    /// `incision_threshold` and `sediment_retention`. Skipped on write when
    /// false (byte-identical).
    #[serde(default, skip_serializing_if = "is_false")]
    pub incision_threshold_alluvial_channels: bool,
    /// INCISION-THRESHOLD part 1 (lane IT): `true` selects
    /// `marine_remainder_layer = once-per-step-capped-global-layer-v1`.
    /// Skipped on write when false.
    #[serde(default, skip_serializing_if = "is_false")]
    pub once_per_step_marine_layer: bool,
    /// GLACIAL-T0 (lane GLACIAL-T0, 2026-10-21): selects `glacial_erosion =
    /// glacial-stream-power-blend-v1`. Skipped on write when absent (every
    /// existing block serializes as before). Needs `nonlinear_subgrid_erosion`
    /// and a `climate_forcing` block.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub glacial_erosion: Option<GlacialErosionBlock>,
}

/// One warm (non-glacial) climate state of `glacial-stream-power-blend-v1`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GlacialWarmStateBlock {
    /// Fraction of genesis time in this state.
    pub duty: f64,
    /// Warming relative to the genesis climate (K).
    pub warming_k: f64,
}

/// Parameters of `glacial-stream-power-blend-v1` (see laws.rs and
/// `pf_surface::erosion::glacial`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GlacialErosionBlock {
    pub glacial_to_fluvial_erodibility_ratio: f64,
    pub full_ice_height_above_ela_m: f64,
    /// Fraction of genesis time in the coldhouse with glacial cycles.
    pub coldhouse_duty: f64,
    pub glacial_max_cooling_k: f64,
    pub glacial_max_precipitation_factor: f64,
    pub warm_states: Vec<GlacialWarmStateBlock>,
}

impl GlacialErosionBlock {
    /// Lane GLACIAL-T0 PREREG.txt (2026-10-21 11:25 EDT), declared before
    /// any run. Sources in the resolved constants' labels.
    pub fn v1() -> Self {
        Self {
            glacial_to_fluvial_erodibility_ratio: 7.6,
            full_ice_height_above_ela_m: 500.0,
            coldhouse_duty: 0.13,
            glacial_max_cooling_k: 7.0,
            glacial_max_precipitation_factor: 0.7,
            warm_states: vec![
                GlacialWarmStateBlock { duty: 0.18, warming_k: 6.0 },
                GlacialWarmStateBlock { duty: 0.28, warming_k: 9.5 },
                GlacialWarmStateBlock { duty: 0.41, warming_k: 14.0 },
            ],
        }
    }
}

/// Parameters of `stream-power-excess-steepness-threshold-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct IncisionThresholdBlock {
    pub critical_channel_steepness_m: f64,
}

impl IncisionThresholdBlock {
    /// Lane IT PREREG.txt (2026-10-20 23:42 EDT): 30 m and 60 m declared
    /// before any run; v1 is the lower value.
    pub fn v1() -> Self {
        Self {
            critical_channel_steepness_m: 30.0,
        }
    }
}

impl OrogenyBlock {
    /// Both v1 laws with their v1 values.
    pub fn v1() -> Self {
        Self {
            rate_limited_belt: Some(RateLimitedBeltBlock::v1()),
            subgrid_erosion: Some(SubgridErosionBlock::v1()),
            nonlinear_subgrid_erosion: None,
            critical_slope_erosion: None,
            sediment_retention: None,
            incision_threshold: None,
            incision_threshold_trunk_only: false,
            incision_threshold_alluvial_channels: false,
            once_per_step_marine_layer: false,
            glacial_erosion: None,
        }
    }

    /// SEED-CARRY-CC candidate (lane CC): `v1()` with the plate carry
    /// `front-rescue-carry-v1` selected (no new constant).
    pub fn carry_v1() -> Self {
        let mut o = Self::v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.front_rescue_carry = true;
        }
        o
    }

    /// SEDIMENT-TOWERS (lane ST) corrected carry profile: `carry_v1()` with
    /// `collision_sediment = accreted-to-belt-crust-v1` (no new constant).
    pub fn carry_v2() -> Self {
        let mut o = Self::carry_v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.accrete_collision_sediment = true;
        }
        o
    }

    /// SEDIMENT-TOWERS (lane ST) corrected profile: `v1()` with
    /// `collision_sediment = accreted-to-belt-crust-v1` (no new constant).
    pub fn v2() -> Self {
        let mut o = Self::v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.accrete_collision_sediment = true;
        }
        o
    }

    /// SEED-EXTRAS (lane EX) candidate: `carry_v2()` plus hotspot volcanism
    /// `fixed-plume-chains-v1`, continental rifting `stretched-rift-v1` and
    /// plate splits `plume-head-continental-split-v1`, all at their v1
    /// values (fixed before any run, `unpublished-design-record
    /// PREREG.txt`).
    pub fn extras_v1() -> Self {
        let mut o = Self::carry_v2();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.plume_volcanism = Some(PlumeVolcanismBlock::v1());
            b.continental_rift = Some(ContinentalRiftBlock::v1());
        }
        o
    }

    /// SEED-RIFT-SHOULDERS (lane RS, 2026-10-13) candidate: `extras_v1()`
    /// plus `rift_flank = shoulder-uplift-thermal-v1` at its v1 values
    /// (fixed before any run, `unpublished-design-record
    /// PREREG.txt`).
    pub fn extras_v2() -> Self {
        let mut o = Self::extras_v1();
        if let Some(r) = o
            .rate_limited_belt
            .as_mut()
            .and_then(|b| b.continental_rift.as_mut())
        {
            r.rift_flank = Some(RiftFlankBlock::v1());
        }
        o
    }

    /// SEED-TRENCH (lane TN, 2026-10-13) candidate: `extras_v2()` plus
    /// `subduction_trench = slab-flexure-trench-v1` at its v1 values (fixed
    /// before any run, `unpublished-design-record`).
    pub fn extras_v3() -> Self {
        let mut o = Self::extras_v2();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.subduction_trench = Some(SubductionTrenchBlock::v1());
            if let Some(f) = b.continental_rift.as_mut().and_then(|r| r.rift_flank.as_mut()) {
                f.seeded_margin_subsidence = true;
            }
        }
        o
    }

    /// SEED-CARRY-CONSERVE (lane CV, 2026-10-13) candidate: `extras_v3()`
    /// with the plate carry `continental-parcel-carry-v1` (in place of
    /// `front-rescue-carry-v1`) and `plate_suture =
    /// contact-length-pole-share-v1`, at their v1 values (fixed before any
    /// run, `unpublished-design-record`).
    pub fn conserve_v1() -> Self {
        let mut o = Self::extras_v3();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.front_rescue_carry = false;
            b.continental_parcel_carry = Some(ContinentalParcelCarryBlock::v1());
            b.plate_suture = Some(PlateSutureBlock::v1());
        }
        o
    }

    /// SEED-FLOODBASALT (lane FB, 2026-10-13) candidate: `extras_v3()`
    /// plus `flood_basalt = plume-head-lip-swell-v1` at its v1 values
    /// (fixed before any run, `unpublished-design-record`).
    pub fn extras_v3_flood() -> Self {
        let mut o = Self::extras_v3();
        if let Some(p) = o
            .rate_limited_belt
            .as_mut()
            .and_then(|b| b.plume_volcanism.as_mut())
        {
            p.flood_basalt = Some(FloodBasaltBlock::v1());
        }
        o
    }

    /// SEED-RELIEF-LAW candidate: `rate-limited-belt-v1` plus
    /// `stream-power-subgrid-nonlinear-v1` (in place of
    /// `stream-power-subgrid-v1`), all at their pre-registered values.
    pub fn relief_v1() -> Self {
        Self {
            rate_limited_belt: Some(RateLimitedBeltBlock::v1()),
            subgrid_erosion: None,
            nonlinear_subgrid_erosion: Some(NonlinearSubgridErosionBlock::v1()),
            critical_slope_erosion: None,
            sediment_retention: None,
            incision_threshold: None,
            incision_threshold_trunk_only: false,
            incision_threshold_alluvial_channels: false,
            once_per_step_marine_layer: false,
            glacial_erosion: None,
        }
    }

    /// SEED-UPLANDS (lane UP, 2026-10-13) candidate: `extras_v2()` with
    /// the erosion law `stream-power-subgrid-nonlinear-v1` (lane RL v1
    /// values, in place of OG's n = 1) plus `critical-slope-envelope-v1`
    /// at its v1 values (fixed before any run,
    /// `unpublished-design-record`).
    pub fn uplands_v1() -> Self {
        let mut o = Self::extras_v2();
        o.subgrid_erosion = None;
        o.nonlinear_subgrid_erosion = Some(NonlinearSubgridErosionBlock::v1());
        o.critical_slope_erosion = Some(CriticalSlopeErosionBlock::v1());
        o
    }

    /// SEED-SEDIMENT-RETAIN (lane SR, 2026-10-13) candidate: `extras_v3()`
    /// plus `continental-retention-base-level-v1` at its v1 values (fixed
    /// before any run, `unpublished-design-record`).
    pub fn retain_v1() -> Self {
        let mut o = Self::extras_v3();
        o.sediment_retention = Some(SedimentRetentionBlock::v1());
        o
    }

    /// SEED-PLATE-REORG (lane PR, 2026-10-13): `conserve_v1()` plus
    /// `continental-retention-base-level-v1` at its v1 values (the orogeny
    /// half of `oldstart_retain_v1`; no new constant).
    pub fn conserve_retain_v1() -> Self {
        let mut o = Self::conserve_v1();
        o.sediment_retention = Some(SedimentRetentionBlock::v1());
        o
    }

    /// SEED-PLATE-REORG (lane PR, 2026-10-13) candidate: `conserve_retain_v1()`
    /// plus `plate-reorg-scripted-v1` at its v1 values (fixed before any
    /// run, `unpublished-design-record`).
    pub fn cycle_v1() -> Self {
        let mut o = Self::conserve_retain_v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.plate_reorg = Some(PlateReorgBlock::v1());
        }
        o
    }

    /// SEED-MAGMATISM (lane MG, 2026-10-13) candidate: `cycle_v1()` plus
    /// `magmatism = slab-arc-backarc-llsvp-v1` at its v1 values; the belt's
    /// contact arc (`arc_productivity_m` 500 m) is set to 0 because the law
    /// places the arc crust at slab depth instead (fixed before any run,
    /// `unpublished-design-record`).
    pub fn cycle_magma_v1() -> Self {
        let mut o = Self::cycle_v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.arc_productivity_m = 0.0;
            b.magmatism = Some(MagmatismBlock::v1());
        }
        o
    }

    /// SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13) candidate:
    /// `conserve_retain_v1()` plus `plate-reorg-scripted-v2` at its v2 values
    /// (`unpublished-design-record`).
    pub fn cycle_v2() -> Self {
        let mut o = Self::conserve_retain_v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.plate_reorg = Some(PlateReorgBlock::v2());
        }
        o
    }

    /// SEED-RANGES (lane RG, 2026-10-13) candidate: `cycle_v1()` with the
    /// suture law `contact-length-indentation-taper-v1` (post-collision
    /// indentation) and `orogen_plateau = plateau-cap-spill-v1` at their v1
    /// values; the belt's delamination backstop rises to 100 km so the
    /// plateau law handles continental crust above its 65 km cap (fixed
    /// before any run, `unpublished-design-record`).
    pub fn cycle_ranges_v1() -> Self {
        let mut o = Self::cycle_v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.delamination_thickness_m = 100_000.0;
            b.plateau = Some(OrogenPlateauBlock::v1());
            if let Some(s) = b.plate_suture.as_mut() {
                s.indentation = Some(SutureIndentationBlock::v1());
            }
        }
        o
    }

    /// SEED-RANGES (lane RG, 2026-10-13) composition C: `cycle_ranges_v1()`
    /// with lane UP's erosion (`stream-power-subgrid-nonlinear-v1` at RL v1
    /// values in place of OG's n = 1, plus `critical-slope-envelope-v1` at
    /// its v1 values); no new constant.
    pub fn cycle_ranges_uplands_v1() -> Self {
        let mut o = Self::cycle_ranges_v1();
        o.subgrid_erosion = None;
        o.nonlinear_subgrid_erosion = Some(NonlinearSubgridErosionBlock::v1());
        o.critical_slope_erosion = Some(CriticalSlopeErosionBlock::v1());
        o
    }

    /// SEED-TRANSFORMS (lane TF, 2026-10-13) candidate: `world_v1()` with
    /// `plate_boundary = smoothed-normal-transform-v1` and
    /// `plate-reorg-scripted-v4`, both at their v1 values
    /// (`unpublished-design-record`).
    pub fn world_v2() -> Self {
        let mut o = Self::world_v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.plate_reorg = Some(PlateReorgBlock::v4());
            b.plate_boundary = Some(PlateBoundaryBlock::v1());
        }
        o
    }
    /// SEED-OLD-OCEAN (lane OR, 2026-10-13) candidate: `world_v2()` with
    /// `plate-reorg-scripted-v5` (passive-margin subduction initiation
    /// driven by ocean age) at its v5 values
    /// (`unpublished-design-record`).
    pub fn world_v3() -> Self {
        let mut o = Self::world_v2();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.plate_reorg = Some(PlateReorgBlock::v5());
        }
        o
    }
    /// SEED-WORLD-COMPOSE (lane CP, 2026-10-13) candidate: every approved
    /// piece together. `conserve_retain_v1()` (oldstart base half) plus
    /// `plate-reorg-scripted-v3` (v2 cut, symmetric ridges, age-driven
    /// initiation), the indentation suture law v1, `plateau-cap-spill-v2`
    /// (inside the belt step; the delamination backstop back at the v1
    /// 70 km), lane UP's erosion (n = 2 + envelope), `slab-arc-backarc-
    /// llsvp-v2` (belt contact arc 0 as in cycle_magma_v1) and
    /// `plume-head-lip-swell-v1` on the same (LLSVP-placed, drifting)
    /// plumes (`unpublished-design-record`).
    pub fn world_v1() -> Self {
        let mut o = Self::conserve_retain_v1();
        o.subgrid_erosion = None;
        o.nonlinear_subgrid_erosion = Some(NonlinearSubgridErosionBlock::v1());
        o.critical_slope_erosion = Some(CriticalSlopeErosionBlock::v1());
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.plate_reorg = Some(PlateReorgBlock::v3());
            b.plateau = Some(OrogenPlateauBlock::v2());
            if let Some(s) = b.plate_suture.as_mut() {
                s.indentation = Some(SutureIndentationBlock::v1());
            }
            b.arc_productivity_m = 0.0;
            b.magmatism = Some(MagmatismBlock::v2());
            if let Some(p) = b.plume_volcanism.as_mut() {
                p.flood_basalt = Some(FloodBasaltBlock::v1());
            }
        }
        o
    }

    /// SEED-STRATIGRAPHY stage A (lane ST, 2026-10-13) candidate:
    /// `world_v1()` plus `stratigraphy = parcel-stacks-v1` at its v1 values
    /// (`unpublished-design-record`). A passive record:
    /// the geography is world_v1's.
    pub fn world_strata_v1() -> Self {
        let mut o = Self::world_v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.stratigraphy = Some(StratigraphyBlock::v1());
        }
        o
    }

    /// SEED-STRATIGRAPHY stage A.1 (lane SA, 2026-10-13) candidate:
    /// `world_v1()` plus `stratigraphy = parcel-stacks-v2`
    /// (`unpublished-design-record`). A passive record:
    /// the geography is world_v1's.
    pub fn world_strata_v2() -> Self {
        let mut o = Self::world_v1();
        if let Some(b) = o.rate_limited_belt.as_mut() {
            b.stratigraphy = Some(StratigraphyBlock::v2());
        }
        o
    }
}

/// Observation of complete native donor-first passes at explicitly selected steps.
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct FluvialTraceBlock {
    #[serde(default,skip_serializing_if="std::ops::Not::not")]
    pub stage_budgets:bool,
    pub steps:Vec<u64>,
    pub record_limit:u64,
    pub byte_limit:u64,
}

/// Parameters of `continental-retention-base-level-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SedimentRetentionBlock {
    #[serde(default,skip_serializing_if="Option::is_none")]
    pub fluvial_trace:Option<FluvialTraceBlock>,
    /// SHELF (lane SHELF, 2026-10-20): opt-in `wave-base-shelf-grading-v1`
    /// (pf-geodynamics shelf_grade.rs): after each erosion step, connected
    /// ocean shallower than `shelf_wave_base_depth_m` loses sediment to wave
    /// base and the volume progrades onto deeper ocean cells. No new
    /// constant. Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub shelf_wave_grading: bool,
    /// INLAND-BASINS (lane INLAND-BASINS, 2026-10-21): opt-in
    /// `closed-basin-water-balance-v1` (pf-surface closed_basin.rs): lake
    /// water levels from the closed-basin water balance and focused floor
    /// fill of lake sediment. Constants are the law's v1 values (law set).
    /// Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub closed_basin_water_balance: bool,
    /// LAKE-EXCESS (lane LAKE-EXCESS, 2026-10-21): opt-in
    /// `filled-basin-alluvial-grade-v1` on top of
    /// `closed_basin_water_balance` (required): a lake filled to within
    /// d_c of its spill level grades to its outlet as an alluvial plain.
    /// Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub filled_basin_alluvial_grading: bool,
    /// LAKE-OUTLET (lane LAKE-OUTLET, 2026-10-21): opt-in
    /// `lake-outlet-sill-incision-v1` on top of `closed_basin_water_balance`
    /// (required): each full lake's outlet incises its sill by channel
    /// stream power on the lake's overflow. Absent/false = byte-identical
    /// to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub lake_outlet_incision: bool,
    /// SILL-REBOUND (lane SILL-REBOUND, 2026-10-21): opt-in
    /// `lake-outlet-sill-incision-flexural-support-v1` on top of
    /// `lake_outlet_incision` (required): no local isostatic rebound at an
    /// outlet sill (removal sized for the post-isostasy surface).
    /// Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub lake_outlet_flexural_support: bool,
    /// LAKE-HOLD (lane LAKE-HOLD, 2026-10-21): opt-in
    /// `lake-outlet-path-incision-v1` on top of `lake_outlet_incision`
    /// (required): the outlet channel law acts on the whole overflow path
    /// from the sill to the next water body, not on the sill alone (with
    /// `lake_outlet_flexural_support`: on every incised path cell).
    /// Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub lake_outlet_path_incision: bool,
    /// LAKE-HOLD-2 (lane LAKE-HOLD-2, 2026-10-21): opt-in
    /// `lake-outlet-path-through-lakes-v1` on top of
    /// `lake_outlet_path_incision` (required): a river entering a full lake
    /// whose sill the same step cuts is graded to that lake's post-incision
    /// spill, so the cut passes down a chain of lakes within one step (with
    /// `lake_outlet_flexural_support`: the flexural-support variant).
    /// Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub lake_outlet_path_through_lakes: bool,
    /// LAKE-HOLD-2 (lane LAKE-HOLD-2, 2026-10-21): opt-in
    /// `lake-outlet-path-threshold-grade-v1` on top of
    /// `lake_outlet_path_incision` and the incision threshold (both
    /// required): rule 2a's grade on the outlet path is min(S_a, the outlet
    /// channel's threshold slope). Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub lake_outlet_path_grade: bool,
    pub base_level_relief_scale_m: f64,
    pub base_level_taper_exponent: f64,
    pub alluvial_grade_slope: f64,
    pub shelf_wave_base_depth_m: f64,
}

impl SedimentRetentionBlock {
    /// v1 values fixed BEFORE any run (lane SR PREREG.txt, 2026-10-13
    /// 09:47 EDT).
    pub fn v1() -> Self {
        Self {
            fluvial_trace: None,
            shelf_wave_grading: false,
            closed_basin_water_balance: false,
            filled_basin_alluvial_grading: false,
            lake_outlet_incision: false,
            lake_outlet_flexural_support: false,
            lake_outlet_path_incision: false,
            lake_outlet_path_through_lakes: false,
            lake_outlet_path_grade: false,
            base_level_relief_scale_m: 1000.0,
            base_level_taper_exponent: 2.0,
            alluvial_grade_slope: 1.0e-4,
            shelf_wave_base_depth_m: 100.0,
        }
    }
}

/// Parameters of `critical-slope-envelope-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CriticalSlopeErosionBlock {
    pub hillslope_critical_slope: f64,
    pub knickpoint_celerity_m_per_my: f64,
}

impl CriticalSlopeErosionBlock {
    /// v1 values fixed BEFORE any run (lane UP PREREG.txt, 2026-10-13
    /// 08:55 EDT).
    pub fn v1() -> Self {
        Self {
            hillslope_critical_slope: 0.04,
            knickpoint_celerity_m_per_my: 70_000.0,
        }
    }
}

/// Parameters of `stream-power-subgrid-nonlinear-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct NonlinearSubgridErosionBlock {
    pub channel_head_length_m: f64,
    pub slope_exponent_n: f64,
    pub reference_steepness_m: f64,
}

impl NonlinearSubgridErosionBlock {
    /// v1 values fixed BEFORE any run (lane RL
    /// `unpublished-design-record`, 2026-10-12 18:40 EDT).
    pub fn v1() -> Self {
        Self {
            channel_head_length_m: 100.0,
            slope_exponent_n: 2.0,
            reference_steepness_m: 500.0,
        }
    }
}

/// Parameters of `rate-limited-belt-v1` (see laws.rs for the model).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RateLimitedBeltBlock {
    pub belt_shear_zone_dip_deg: f64,
    pub lower_crust_channel_viscosity_pa_s: f64,
    pub channel_reference_thickness_m: f64,
    pub delamination_thickness_m: f64,
    pub arc_productivity_m: f64,
    pub arc_belt_width_m: f64,
    /// SEED-CARRY-CC (lane CC, 2026-10-12): `true` selects the law-set entry
    /// `plate_carry = front-rescue-carry-v1` (whose consumed collision volume
    /// is this belt's supply, so it lives with the belt). Absent / `false` =
    /// the legacy carry; skipped on write when false, so every existing block
    /// serializes byte-identically.
    #[serde(default, skip_serializing_if = "is_false")]
    pub front_rescue_carry: bool,
    /// SEDIMENT-TOWERS (lane ST, 2026-10-12): `true` selects the law-set
    /// entry `collision_sediment = accreted-to-belt-crust-v1` (collision
    /// sediment joins this belt's crust supply, so it lives with the belt).
    /// Absent / `false` = legacy `surface-layer-v0`; skipped on write when
    /// false, so every existing block serializes byte-identically.
    #[serde(default, skip_serializing_if = "is_false")]
    pub accrete_collision_sediment: bool,
    /// SEED-EXTRAS (lane EX, 2026-10-12): selects `hotspot_volcanism =
    /// fixed-plume-chains-v1`. Lives with the belt because the plume pass
    /// runs inside the belt-law advance (refused with the legacy stacking
    /// rule). Skipped on write when absent (every existing block serializes
    /// byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plume_volcanism: Option<PlumeVolcanismBlock>,
    /// SEED-EXTRAS (lane EX, 2026-10-12): selects `continental_rift =
    /// stretched-rift-v1` (and, with its `plume_head_split`, `plate_
    /// reorganization = plume-head-continental-split-v1`). Skipped on write
    /// when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub continental_rift: Option<ContinentalRiftBlock>,
    /// SEED-TRENCH (lane TN, 2026-10-13): selects `subduction_trench =
    /// slab-flexure-trench-v1`. Lives with the belt like the other
    /// boundary laws of the candidate profiles (its trench is the oceanic
    /// side of the belt's subduction boundaries). Skipped on write when
    /// absent (every existing block serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subduction_trench: Option<SubductionTrenchBlock>,
    /// SEED-CARRY-CONSERVE (lane CV, 2026-10-13): selects `plate_carry =
    /// continental-parcel-carry-v1` (its collision consumption is this
    /// belt's supply, so it lives with the belt; refused together with
    /// `front_rescue_carry`). Skipped on write when absent (every existing
    /// block serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub continental_parcel_carry: Option<ContinentalParcelCarryBlock>,
    /// SEED-CARRY-CONSERVE (lane CV, 2026-10-13): selects `plate_suture =
    /// contact-length-pole-share-v1`. Skipped on write when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plate_suture: Option<PlateSutureBlock>,
    /// SEED-PLATE-REORG (lane PR, 2026-10-13): selects
    /// `plate_reorganization = plate-reorg-scripted-v1` (supersedes the
    /// plume-head split, whose opening rate it keeps; needs
    /// `continental_rift.plume_head_split`). Skipped on write when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plate_reorg: Option<PlateReorgBlock>,
    /// SEED-TRANSFORMS (lane TF, 2026-10-13): selects `plate_boundary =
    /// smoothed-normal-transform-v1`. Skipped on write when absent (every
    /// existing block serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plate_boundary: Option<PlateBoundaryBlock>,
    /// SEED-RANGES (lane RG, 2026-10-13): selects `orogen_plateau =
    /// plateau-cap-spill-v1` (acts on this belt's thickened crust).
    /// Skipped on write when absent (every existing block serializes
    /// byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plateau: Option<OrogenPlateauBlock>,
    /// SEED-MAGMATISM (lane MG, 2026-10-13): selects `magmatism =
    /// slab-arc-backarc-llsvp-v1` (its arcs replace this belt's contact arc in
    /// the candidate profile). Skipped on write when absent (every existing
    /// block serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub magmatism: Option<MagmatismBlock>,
    /// SEED-STRATIGRAPHY (lane ST, 2026-10-13): selects `stratigraphy =
    /// parcel-stacks-v1` (the stacks live on this belt's continental
    /// parcels). Skipped on write when absent (every existing block
    /// serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stratigraphy: Option<StratigraphyBlock>,
}

/// Parameters of `parcel-stacks-v1` (see laws.rs; SEED-STRATIGRAPHY stage
/// A, lane ST, 2026-10-13).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct StratigraphyBlock {
    pub strata_max_layers: f64,
    pub strata_merge_gap_my: f64,
    pub strata_unconformity_gap_my: f64,
    pub strata_min_layer_m: f64,
    pub strata_trace_m: f64,
    pub strata_arc_extrusive_fraction: f64,
    pub strata_plume_extrusive_fraction: f64,
    pub strata_shortening_pluton_fraction: f64,
    /// PREREG revision 1 (lane ST, 2026-10-13 16:25 EDT).
    pub strata_pluton_roof_m: f64,
    pub strata_metamorphic_depth_m: f64,
    pub strata_high_grade_depth_m: f64,
    pub strata_shelf_depth_m: f64,
    pub strata_carbonate_max_paleolat_deg: f64,
    pub strata_carbonate_max_supply_m_per_my: f64,
    pub strata_arid_lat_min_deg: f64,
    pub strata_arid_lat_max_deg: f64,
    pub strata_humid_lat_max_deg: f64,
    pub strata_humid_midlat_min_deg: f64,
    pub strata_humid_midlat_max_deg: f64,
    pub strata_coal_max_elevation_m: f64,
    pub strata_conglomerate_relief_m: f64,
    pub strata_river_mouth_area_m2: f64,
    pub strata_publish_depth_m: f64,
    pub strata_basement_min_m: f64,
    /// SEED-STRATIGRAPHY stage A.1 (lane SA, 2026-10-13): present = law
    /// `parcel-stacks-v2` (these constants on top of the v1 ones); skipped
    /// on write when absent, so every v1 block serializes byte-identically.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub v2: Option<StratigraphyV2Block>,
}

/// The constants `parcel-stacks-v2` adds (lane SA PREREG.txt,
/// `unpublished-design-record`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct StratigraphyV2Block {
    pub strata_shelf_width_m: f64,
    pub strata_slope_foot_depth_m: f64,
    pub strata_shelf_majority_fraction: f64,
    pub strata_platform_cover_max_m: f64,
    pub strata_platform_full_elev_m: f64,
    pub strata_platform_zero_elev_m: f64,
    pub strata_platform_sandstone_fraction: f64,
    pub strata_platform_carbonate_fraction: f64,
    pub strata_platform_age_base_my: f64,
    pub strata_platform_age_top_my: f64,
}

impl StratigraphyV2Block {
    /// v2 values fixed BEFORE any run (lane SA PREREG.txt, 2026-10-13
    /// 17:40 EDT).
    pub fn v2() -> Self {
        Self {
            strata_shelf_width_m: 80000.0,
            strata_slope_foot_depth_m: 4000.0,
            strata_shelf_majority_fraction: 0.5,
            strata_platform_cover_max_m: 2000.0,
            strata_platform_full_elev_m: 0.0,
            strata_platform_zero_elev_m: 600.0,
            strata_platform_sandstone_fraction: 0.3,
            strata_platform_carbonate_fraction: 0.3,
            strata_platform_age_base_my: -600.0,
            strata_platform_age_top_my: -50.0,
        }
    }
}

impl StratigraphyBlock {
    /// v1 values fixed BEFORE any run (lane ST PREREG.txt, 2026-10-13
    /// 16:00 EDT).
    pub fn v1() -> Self {
        Self {
            strata_max_layers: 12.0,
            strata_merge_gap_my: 2.0,
            strata_unconformity_gap_my: 5.0,
            strata_min_layer_m: 10.0,
            strata_trace_m: 0.01,
            strata_arc_extrusive_fraction: 0.3,
            strata_plume_extrusive_fraction: 0.2,
            strata_shortening_pluton_fraction: 0.1,
            strata_pluton_roof_m: 4000.0,
            strata_metamorphic_depth_m: 10000.0,
            strata_high_grade_depth_m: 20000.0,
            strata_shelf_depth_m: 200.0,
            strata_carbonate_max_paleolat_deg: 30.0,
            strata_carbonate_max_supply_m_per_my: 50.0,
            strata_arid_lat_min_deg: 15.0,
            strata_arid_lat_max_deg: 35.0,
            strata_humid_lat_max_deg: 15.0,
            strata_humid_midlat_min_deg: 40.0,
            strata_humid_midlat_max_deg: 65.0,
            strata_coal_max_elevation_m: 200.0,
            strata_conglomerate_relief_m: 1000.0,
            strata_river_mouth_area_m2: 1.0e11,
            strata_publish_depth_m: 3000.0,
            strata_basement_min_m: 5000.0,
            v2: None,
        }
    }
    /// `parcel-stacks-v2` (lane SA PREREG.txt, 2026-10-13 17:40 EDT): the v1
    /// values with a 16-layer cap, plus the v2 constants.
    pub fn v2() -> Self {
        Self {
            strata_max_layers: 16.0,
            v2: Some(StratigraphyV2Block::v2()),
            ..Self::v1()
        }
    }
}

/// Parameters of `slab-arc-backarc-llsvp-v1` (see laws.rs; SEED-MAGMATISM,
/// lane MG, 2026-10-13).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct MagmatismBlock {
    pub slab_depth_below_arc_m: f64,
    pub slab_dip_young_deg: f64,
    pub slab_dip_old_deg: f64,
    pub slab_dip_age_my: f64,
    pub arc_magma_productivity_m: f64,
    pub arc_rate_floor_m_per_my: f64,
    pub arc_rate_cap_m_per_my: f64,
    pub arc_half_width_m: f64,
    pub arc_continental_thickness_m: f64,
    pub arc_thermal_uplift_m: f64,
    pub arc_thermal_half_width_m: f64,
    pub arc_thermal_decay_my: f64,
    pub backarc_slab_age_min_my: f64,
    pub backarc_overriding_advance_max_m_per_my: f64,
    pub backarc_extension_rate_m_per_my: f64,
    pub backarc_axis_offset_m: f64,
    pub backarc_half_width_m: f64,
    pub backarc_min_crust_m: f64,
    pub llsvp_edge_deg: f64,
    pub llsvp_edge_sigma_deg: f64,
    pub llsvp_antipode_jitter_deg: f64,
    pub plume_background_fraction: f64,
    pub plume_drift_rate_m_per_my: f64,
    pub vent_active_window_my: f64,
    pub volcano_spacing_m: f64,
    pub arc_crust_cap_m: f64,
    /// SEED-WORLD-COMPOSE (lane CP, 2026-10-13): 1 = `slab-arc-backarc-llsvp-v1`,
    /// 2 = `slab-arc-backarc-llsvp-v2` (needs the three v2 constants below).
    /// Defaults to 1 and is skipped on write when 1; the v2 constants are
    /// skipped when absent, so every v1 block serializes byte-identically.
    #[serde(default = "law_version_one", skip_serializing_if = "is_law_version_one")]
    pub law_version: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub backarc_episode_my: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub backarc_max_stretch_factor: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kimberlite_plume_radius_m: Option<f64>,
}

fn law_version_one() -> u32 {
    1
}
fn is_law_version_one(v: &u32) -> bool {
    *v == 1
}

impl MagmatismBlock {
    /// v1 values fixed BEFORE any run (lane MG PREREG.txt, 2026-10-13
    /// 13:20 EDT).
    pub fn v1() -> Self {
        Self {
            slab_depth_below_arc_m: 105000.0,
            slab_dip_young_deg: 30.0,
            slab_dip_old_deg: 45.0,
            slab_dip_age_my: 150.0,
            arc_magma_productivity_m: 600.0,
            arc_rate_floor_m_per_my: 20000.0,
            arc_rate_cap_m_per_my: 100000.0,
            arc_half_width_m: 75000.0,
            arc_continental_thickness_m: 20000.0,
            arc_thermal_uplift_m: 1000.0,
            arc_thermal_half_width_m: 150000.0,
            arc_thermal_decay_my: 15.0,
            backarc_slab_age_min_my: 50.0,
            backarc_overriding_advance_max_m_per_my: 0.0,
            backarc_extension_rate_m_per_my: 30000.0,
            backarc_axis_offset_m: 150000.0,
            backarc_half_width_m: 150000.0,
            backarc_min_crust_m: 12000.0,
            llsvp_edge_deg: 35.0,
            llsvp_edge_sigma_deg: 5.0,
            llsvp_antipode_jitter_deg: 15.0,
            plume_background_fraction: 0.15,
            plume_drift_rate_m_per_my: 10000.0,
            vent_active_window_my: 3.0,
            volcano_spacing_m: 70000.0,
            arc_crust_cap_m: 50000.0,
            law_version: 1,
            backarc_episode_my: None,
            backarc_max_stretch_factor: None,
            kimberlite_plume_radius_m: None,
        }
    }

    /// SEED-WORLD-COMPOSE (lane CP): `slab-arc-backarc-llsvp-v2` values
    /// fixed BEFORE any run (lane CP PREREG.txt, 2026-10-13): v1 values plus
    /// a 20 My back-arc episode, a stretching factor of 2 and a 300 km
    /// head-only kimberlite radius (lane MG ideas 2 and 4).
    pub fn v2() -> Self {
        Self {
            law_version: 2,
            backarc_episode_my: Some(20.0),
            backarc_max_stretch_factor: Some(2.0),
            kimberlite_plume_radius_m: Some(300_000.0),
            ..Self::v1()
        }
    }
    /// Teaching opt-in: v2 parameters with volume-conserving arc heat.
    pub fn v3() -> Self {
        Self { law_version: 3, ..Self::v2() }
    }

}

/// Parameters of `plateau-cap-spill-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct OrogenPlateauBlock {
    pub plateau_crust_cap_m: f64,
    pub plateau_spill_reach_m: f64,
    /// SEED-WORLD-COMPOSE (lane CP, 2026-10-13): 1 = `plateau-cap-spill-v1`
    /// (after the advance), 2 = `plateau-cap-spill-v2` (inside the belt
    /// step, crust plus sediment, foothill taper). Defaults to 1 and is
    /// skipped on write when 1; the v2 constants are skipped when absent.
    #[serde(default = "law_version_one", skip_serializing_if = "is_law_version_one")]
    pub law_version: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plateau_foothill_gradient: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plateau_taper_floor_m: Option<f64>,
}

impl OrogenPlateauBlock {
    /// v1 values fixed BEFORE any run (lane RG PREREG.txt, 2026-10-13
    /// 12:19 EDT).
    pub fn v1() -> Self {
        Self {
            plateau_crust_cap_m: 65_000.0,
            plateau_spill_reach_m: 1.0e6,
            law_version: 1,
            plateau_foothill_gradient: None,
            plateau_taper_floor_m: None,
        }
    }

    /// SEED-WORLD-COMPOSE (lane CP): `plateau-cap-spill-v2` values fixed
    /// BEFORE any run (lane CP PREREG.txt, 2026-10-13): v1's cap and reach,
    /// a foothill gradient of 0.1 (10 km of crust per 100 km) and a 40 km
    /// taper floor.
    pub fn v2() -> Self {
        Self {
            law_version: 2,
            plateau_foothill_gradient: Some(0.1),
            plateau_taper_floor_m: Some(40_000.0),
            ..Self::v1()
        }
    }
}

/// Optional observation budget, not a change to native geological laws.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CarrierHistoryBlock {
    #[serde(default,skip_serializing_if="std::ops::Not::not")]
    pub layer_evolution:bool,
    pub selected: Vec<u64>,
    pub event_limit: u64,
    pub byte_limit: u64,
}

/// Parameters of `continental-parcel-carry-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ContinentalParcelCarryBlock {
    /// Opt-in finite body manifest; absent leaves the native profile unchanged.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub native_body_source: Option<crate::emplacement::BodyInput>,
    /// Passive layer formation recorder, omitted when disabled. Requires carrier IDs.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub layer_provenance: bool,
    /// Optional recording quota. None keeps legacy serialized profiles unchanged.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub carrier_record_limit: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stack_history: Option<CarrierHistoryBlock>,
    pub parcel_overlap_cell_fraction: f64,
    pub plate_fragment_min_area_fraction: f64,
    /// HYPSO-LEVERS (2026-10-21, opt-in): re-projection passes of the
    /// write-back residual correction (see `PlateCarryLawV1::
    /// ContinentalParcelCarryV1`). Absent = the v1 write-back.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub write_back_residual_iterations: Option<f64>,
}

impl ContinentalParcelCarryBlock {
    /// v1 values fixed BEFORE any run (lane CV PREREG.txt, 2026-10-13
    /// 09:31 EDT).
    pub fn v1() -> Self {
        Self {
            native_body_source: None,
            layer_provenance: false,
            carrier_record_limit: None,
            stack_history: None,
            parcel_overlap_cell_fraction: 0.5,
            plate_fragment_min_area_fraction: 1.0e-3,
            write_back_residual_iterations: None,
        }
    }
}

/// Parameters of `contact-length-pole-share-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PlateSutureBlock {
    pub suture_contact_length_m: f64,
    /// SEED-RANGES (lane RG, 2026-10-13): `Some` selects
    /// `contact-length-indentation-taper-v1` in place of
    /// `contact-length-pole-share-v1`. Skipped on write when absent (every
    /// existing block serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub indentation: Option<SutureIndentationBlock>,
}

/// Parameters of `contact-length-indentation-taper-v1` beyond the contact
/// length (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SutureIndentationBlock {
    pub indentation_time_my: f64,
    pub indentation_shortening_budget_m: f64,
}

impl SutureIndentationBlock {
    /// v1 values fixed BEFORE any run (lane RG PREREG.txt, 2026-10-13
    /// 12:19 EDT).
    pub fn v1() -> Self {
        Self {
            indentation_time_my: 40.0,
            indentation_shortening_budget_m: 7.5e5,
        }
    }
}

impl PlateSutureBlock {
    /// v1 value fixed BEFORE any run (lane CV PREREG.txt, 2026-10-13
    /// 09:31 EDT).
    pub fn v1() -> Self {
        Self {
            suture_contact_length_m: 2.5e6,
            indentation: None,
        }
    }
}

/// Parameters of `plate-reorg-scripted-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PlateReorgBlock {
    pub reorg_interval_min_my: f64,
    pub reorg_interval_max_my: f64,
    pub slab_pull_bias: f64,
    pub slab_circumference_reference: f64,
    pub reorg_rate_floor_fraction: f64,
    pub supercontinent_share: f64,
    pub supercontinent_tenure_my: f64,
    pub breakup_piece_max_share: f64,
    pub breakup_cascade_my: f64,
    pub rift_persistence_my: f64,
    pub weak_zone_strength_ratio: f64,
    pub plate_count_min: u32,
    pub plate_count_max: u32,
    pub plate_count_min_area_fraction: f64,
    pub subduction_initiation_rate_deg_per_my: f64,
    pub plume_head_interval_my: f64,
    pub plume_tail_lifetime_my: f64,
    pub plume_split_min_share: f64,
    /// SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13): 1 = `plate-reorg-scripted-v1`,
    /// 2 = `plate-reorg-scripted-v2` (same constants). Defaults to 1 and is
    /// skipped on write when 1, so every v1 block serializes byte-identically.
    #[serde(default = "PlateReorgBlock::law_version_v1", skip_serializing_if = "PlateReorgBlock::is_law_version_v1")]
    pub law_version: u32,
    /// SEED-WORLD-COMPOSE (lane CP, 2026-10-13): law_version 3
    /// (`plate-reorg-scripted-v3`) only; skipped on write when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subduction_initiation_age_my: Option<f64>,
    /// SEED-TRANSFORMS (lane TF, 2026-10-13): law_version 4
    /// (`plate-reorg-scripted-v4`) only; skipped on write when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initiation_slab_age_min_my: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initiation_old_slab_fraction_min: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initiation_trench_spacing_m: Option<f64>,
    /// SEED-OLD-OCEAN (lane OR, 2026-10-13): law_version 5
    /// (`plate-reorg-scripted-v5`) only; skipped on write when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub margin_initiation_age_my: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub margin_segment_length_m: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub margin_trench_spacing_m: Option<f64>,
}

impl PlateReorgBlock {
    fn law_version_v1() -> u32 {
        1
    }
    fn is_law_version_v1(v: &u32) -> bool {
        *v == 1
    }
    /// v1 values fixed BEFORE any run (lane PR PREREG.txt, 2026-10-13
    /// 10:55 EDT).
    pub fn v1() -> Self {
        Self {
            reorg_interval_min_my: 50.0,
            reorg_interval_max_my: 100.0,
            slab_pull_bias: 0.5,
            slab_circumference_reference: 0.3,
            reorg_rate_floor_fraction: 0.2,
            supercontinent_share: 0.6,
            supercontinent_tenure_my: 100.0,
            breakup_piece_max_share: 0.4,
            breakup_cascade_my: 20.0,
            rift_persistence_my: 100.0,
            weak_zone_strength_ratio: 0.2,
            plate_count_min: 7,
            plate_count_max: 12,
            plate_count_min_area_fraction: 0.005,
            subduction_initiation_rate_deg_per_my: 0.3,
            plume_head_interval_my: 25.0,
            plume_tail_lifetime_my: 150.0,
            plume_split_min_share: 0.2,
            law_version: 1,
            subduction_initiation_age_my: None,
            initiation_slab_age_min_my: None,
            initiation_old_slab_fraction_min: None,
            initiation_trench_spacing_m: None,
            margin_initiation_age_my: None,
            margin_segment_length_m: None,
            margin_trench_spacing_m: None,
        }
    }

    /// SEED-BREAKUP-CUT-V2 (lane BC): `plate-reorg-scripted-v2` values fixed
    /// BEFORE any run (lane BC PREREG.txt, 2026-10-13 12:20 EDT): v1 values
    /// except the plume-split continental share floor 0.35 (design review review).
    pub fn v2() -> Self {
        Self {
            plume_split_min_share: 0.35,
            law_version: 2,
            ..Self::v1()
        }
    }

    /// SEED-WORLD-COMPOSE (lane CP): `plate-reorg-scripted-v3` values fixed
    /// BEFORE any run (lane CP PREREG.txt, 2026-10-13): v2 values plus a
    /// 160 My subduction-initiation age (lane VL defect 3: Earth's oldest
    /// in-situ ocean floor is ~180 My, p95 ~150-180 My).
    pub fn v3() -> Self {
        Self {
            law_version: 3,
            subduction_initiation_age_my: Some(160.0),
            ..Self::v2()
        }
    }

    /// SEED-TRANSFORMS (lane TF): `plate-reorg-scripted-v4` values fixed
    /// BEFORE any run (lane TF PREREG.txt, 2026-10-13): v3 values plus the
    /// age-triggered initiation limits (old slab >= 100 My along at least
    /// half the new boundary; at most half of it within 1000 km of
    /// existing subduction).
    pub fn v4() -> Self {
        Self {
            law_version: 4,
            initiation_slab_age_min_my: Some(100.0),
            initiation_old_slab_fraction_min: Some(0.5),
            initiation_trench_spacing_m: Some(1_000_000.0),
            ..Self::v3()
        }
    }

    /// SEED-OLD-OCEAN (lane OR): `plate-reorg-scripted-v5` values fixed
    /// BEFORE any run (lane OR PREREG.txt, 2026-10-13): v4 values plus the
    /// passive-margin initiation (ocean >= 160 My at the margin, segments
    /// <= 2000 km, >= 1000 km from existing subduction).
    pub fn v5() -> Self {
        Self {
            law_version: 5,
            margin_initiation_age_my: Some(160.0),
            margin_segment_length_m: Some(2_000_000.0),
            margin_trench_spacing_m: Some(1_000_000.0),
            ..Self::v4()
        }
    }
}

/// Parameters of `smoothed-normal-transform-v1` (see laws.rs; SEED-TRANSFORMS,
/// lane TF, 2026-10-13).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PlateBoundaryBlock {
    pub transform_obliquity_deg: f64,
    pub boundary_normal_smoothing_m: f64,
    /// PLATE-POLES (lane PLATE-POLES, 2026-10-21): opt-in
    /// `rigid-group-interior-v1` (pf-geodynamics boundaries.rs): contacts
    /// between plates of one rigid suture group (bit-identical poles) are
    /// plate interior, not transforms (labels and boundary records only).
    /// No new constant. Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub rigid_group_contacts: bool,
    /// PLATE-POLES: opt-in `rigid-group-contacts-v1` on top of
    /// `rigid_group_contacts` (required): at each reorganization event,
    /// suture-group members detached from the group's main body leave the
    /// group (the event's re-solve gives them their own pole). No new
    /// constant. Absent/false = byte-identical to before.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub rigid_group_release: bool,
}

impl PlateBoundaryBlock {
    /// v1 values fixed BEFORE any run (lane TF PREREG.txt, 2026-10-13).
    pub fn v1() -> Self {
        Self {
            transform_obliquity_deg: 30.0,
            boundary_normal_smoothing_m: 200_000.0,
            rigid_group_contacts: false,
            rigid_group_release: false,
        }
    }
}
/// Parameters of `slab-flexure-trench-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SubductionTrenchBlock {
    pub trench_vertical_load_fraction: f64,
    pub trench_slab_length_max_m: f64,
    pub trench_thermal_expansivity_per_k: f64,
    pub trench_mantle_temperature_c: f64,
    pub trench_elastic_isotherm_c: f64,
    pub trench_thermal_diffusivity_m2_per_yr: f64,
    pub trench_elastic_modulus_pa: f64,
    pub trench_poisson_ratio: f64,
}

impl SubductionTrenchBlock {
    /// v1 values fixed BEFORE any run (lane TN PREREG.txt, 2026-10-13).
    pub fn v1() -> Self {
        Self {
            trench_vertical_load_fraction: 0.1,
            trench_slab_length_max_m: 660_000.0,
            trench_thermal_expansivity_per_k: 3.0e-5,
            trench_mantle_temperature_c: 1300.0,
            trench_elastic_isotherm_c: 450.0,
            trench_thermal_diffusivity_m2_per_yr: 31.5576,
            trench_elastic_modulus_pa: 7.0e10,
            trench_poisson_ratio: 0.25,
        }
    }
}

/// Parameters of `fixed-plume-chains-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PlumeVolcanismBlock {
    pub plume_count: u32,
    pub plume_onset_lead_my: f64,
    pub plume_head_volume_m3: f64,
    pub plume_head_radius_m: f64,
    pub plume_tail_flux_m3_per_my: f64,
    pub plume_edifice_radius_m: f64,
    /// SEED-FLOODBASALT (lane FB, 2026-10-13): selects `flood_basalt =
    /// plume-head-lip-swell-v1` (lives with the plumes it needs). Skipped on
    /// write when absent (every existing block serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flood_basalt: Option<FloodBasaltBlock>,
}

/// Parameters of `plume-head-lip-swell-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct FloodBasaltBlock {
    pub flood_basalt_extrusive_volume_m3: f64,
    pub flood_basalt_intrusive_ratio: f64,
    pub flood_basalt_radius_m: f64,
    pub flood_basalt_taper_fraction: f64,
    pub flood_basalt_eruption_duration_my: f64,
    pub swell_tail_amplitude_m: f64,
    pub swell_tail_radius_m: f64,
    pub swell_head_amplitude_m: f64,
    pub swell_head_radius_m: f64,
    pub swell_head_decay_my: f64,
}

impl FloodBasaltBlock {
    /// v1 values fixed BEFORE any run (lane FB PREREG.txt, 2026-10-13
    /// 09:27 EDT).
    pub fn v1() -> Self {
        Self {
            flood_basalt_extrusive_volume_m3: 1.5e15,
            flood_basalt_intrusive_ratio: 4.0,
            flood_basalt_radius_m: 800_000.0,
            flood_basalt_taper_fraction: 0.3,
            flood_basalt_eruption_duration_my: 2.0,
            swell_tail_amplitude_m: 800.0,
            swell_tail_radius_m: 500_000.0,
            swell_head_amplitude_m: 700.0,
            swell_head_radius_m: 800_000.0,
            swell_head_decay_my: 10.0,
        }
    }
}

impl PlumeVolcanismBlock {
    /// v1 values fixed BEFORE any run (lane EX PREREG.txt, 2026-10-12
    /// 23:37 EDT).
    pub fn v1() -> Self {
        Self {
            plume_count: 12,
            plume_onset_lead_my: 50.0,
            plume_head_volume_m3: 2.0e15,
            plume_head_radius_m: 500_000.0,
            plume_tail_flux_m3_per_my: 1.0e14,
            plume_edifice_radius_m: 100_000.0,
            flood_basalt: None, // SEED-FLOODBASALT (lane FB): opt-in
        }
    }
}

/// Parameters of `stretched-rift-v1` and, optionally,
/// `plume-head-continental-split-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ContinentalRiftBlock {
    pub rift_half_width_m: f64,
    pub breakup_crust_thickness_m: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plume_head_split: Option<PlumeHeadSplitBlock>,
    /// SEED-RIFT-SHOULDERS (lane RS, 2026-10-13): selects `rift_flank =
    /// shoulder-uplift-thermal-v1`. Skipped on write when absent (every
    /// existing block serializes byte-identically).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rift_flank: Option<RiftFlankBlock>,
}

/// Parameters of `shoulder-uplift-thermal-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RiftFlankBlock {
    pub rift_flank_flexural_parameter_m: f64,
    pub rift_shoulder_ratio: f64,
    /// SEED-TRENCH part 2 (lane TN, 2026-10-13): `true` selects
    /// `passive_margin = seeded-margin-thermal-subsidence-v1` (the seeded
    /// rifted margins get this law's stretching markers at t = 0). Skipped
    /// on write when false (every existing block serializes identically).
    #[serde(default, skip_serializing_if = "is_false")]
    pub seeded_margin_subsidence: bool,
}

impl RiftFlankBlock {
    /// v1 values fixed BEFORE any run (lane RS PREREG.txt, 2026-10-13).
    pub fn v1() -> Self {
        Self {
            rift_flank_flexural_parameter_m: 60_000.0,
            rift_shoulder_ratio: 1.0,
            seeded_margin_subsidence: false,
        }
    }
}

/// Parameters of `plume-head-continental-split-v1` (see laws.rs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PlumeHeadSplitBlock {
    pub rift_opening_rate_deg_per_my: f64,
}

impl ContinentalRiftBlock {
    /// v1 values fixed BEFORE any run (lane EX PREREG.txt), with the
    /// plume-head split selected.
    pub fn v1() -> Self {
        Self {
            rift_half_width_m: 100_000.0,
            breakup_crust_thickness_m: 10_000.0,
            plume_head_split: Some(PlumeHeadSplitBlock {
                rift_opening_rate_deg_per_my: 0.25,
            }),
            rift_flank: None,
        }
    }
}

fn is_false(b: &bool) -> bool {
    !*b
}

impl RateLimitedBeltBlock {
    /// v1 values, fixed BEFORE any run (lane OG falsifier file,
    /// `unpublished-design-record`). Every value is recalled
    /// or chosen and awaiting review (ADR 039 C5); see `resolve_universe`.
    pub fn v1() -> Self {
        Self {
            belt_shear_zone_dip_deg: 30.0,
            lower_crust_channel_viscosity_pa_s: 1.0e19,
            channel_reference_thickness_m: 35_000.0,
            delamination_thickness_m: 70_000.0,
            arc_productivity_m: 500.0,
            arc_belt_width_m: 100_000.0,
            front_rescue_carry: false,
            accrete_collision_sediment: false,
            plume_volcanism: None,
            continental_rift: None,
            subduction_trench: None,
            continental_parcel_carry: None,
            plate_suture: None,
            plate_reorg: None,
            plate_boundary: None,
            plateau: None,
            magmatism: None,
            stratigraphy: None,
        }
    }
}

/// Parameters of `stream-power-subgrid-v1` (see laws.rs for the model).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SubgridErosionBlock {
    pub channel_head_length_m: f64,
}

impl SubgridErosionBlock {
    pub fn v1() -> Self {
        Self {
            channel_head_length_m: 100.0,
        }
    }
}

/// SEED-EXTRAS D7 (lane SD, 2026-10-12): opt-in sediment laws. Each
/// sub-block is skipped on write when absent; `None` everywhere is the
/// legacy law set with byte-identical serialization.
///
/// - `erosion_deposition` selects `LandDepositionLawV1::ErosionDepositionV1`
///   (absent: legacy detachment-limited).
/// - `shelf` selects `SedimentPlacementLawV1::ShelfWaveBaseV1` (absent:
///   legacy local fill + global layer).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(deny_unknown_fields)]
pub struct SedimentLawsBlock {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub erosion_deposition: Option<ErosionDepositionBlock>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub shelf: Option<ShelfBlock>,
}

impl SedimentLawsBlock {
    /// Both v1 laws with their v1 values.
    pub fn v1() -> Self {
        Self {
            erosion_deposition: Some(ErosionDepositionBlock::v1()),
            shelf: Some(ShelfBlock::v1()),
        }
    }
}

/// Parameters of `erosion-deposition-v1`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ErosionDepositionBlock {
    pub deposition_coefficient_g: f64,
}

impl ErosionDepositionBlock {
    /// v1 value fixed BEFORE any run (lane SD
    /// `unpublished-design-record`); equals
    /// `pf_surface::erosion::DEPOSITION_COEFFICIENT_G`.
    pub fn v1() -> Self {
        Self {
            deposition_coefficient_g: 1.0,
        }
    }
}

/// Parameters of `shelf-wave-base-v1`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ShelfBlock {
    pub wave_base_depth_m: f64,
}

impl ShelfBlock {
    /// v1 value fixed BEFORE any run (lane SD FALSIFIERS.md); equals
    /// `pf_surface::sediment::SHELF_WAVE_BASE_DEPTH_M`.
    pub fn v1() -> Self {
        Self {
            wave_base_depth_m: 100.0,
        }
    }
}

/// CLIMATE-GENESIS (lane CG, 2026-10-20): selects `climate_forcing =
/// epoch-climate-runoff-v1`. The climate recipe is carried whole; the
/// extra constants of the tier (`SeedClimateConstantsV1`, and for
/// seed-climate-v2 `RotationClimateConstantsV2`) are pinned by digest and
/// verified by the generator. `pf_surface::climate::seed::genesis::
/// climate_forcing_block_v1` builds the v1 block.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ClimateForcingBlock {
    pub climate: ClimateRecipeV1,
    pub seed_constants_digest: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rotation_constants_digest: Option<String>,
    /// Interval between climate solves (My); a whole number of steps.
    pub epoch_my: f64,
    /// Climate solve grid = genesis grid / divisor in each axis (1, 2, 4).
    pub solve_grid_divisor: f64,
    /// Line Gauss-Seidel stop tolerance of a genesis solve.
    pub sweep_rel_tol: f64,
    /// Land runoff floor R >= f P.
    pub runoff_floor_fraction: f64,
    pub between_epochs: ClimateBetweenEpochsV1,
    pub erosivity: ClimateErosivityModeV1,
    /// I11 (lane CLIMATE-GENESIS-2): climate runoff in the subgrid
    /// cell-mean relaxation; absent = the v1 default (reference runoff).
    #[serde(default, skip_serializing_if = "ClimateSubgridCellMeanV1::is_default")]
    pub subgrid_cell_mean: ClimateSubgridCellMeanV1,
    /// HYPSO-LEVERS (2026-10-21, opt-in): exponent p of the erosive runoff
    /// above the reference runoff (law
    /// `erosive-runoff-sublinear-above-reference-v1`). Absent = linear v1.
    /// Refused together with `erosive_runoff_cap_m_per_yr`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub erosive_runoff_exponent: Option<f64>,
    /// HYPSO-LEVERS (2026-10-21, opt-in): cap of the erosive runoff in m/yr
    /// (law `erosive-runoff-cap-v1`). Absent = no cap.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub erosive_runoff_cap_m_per_yr: Option<f64>,
}

pub const NATIVE_PRODUCER_EVIDENCE_V1: &str = "native-producer-evidence-config-v1";

/// Q29: explicitly selected observations and finite terrigenous redistribution.
/// There is no calibrated default and no pelagic source in this configuration.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct NativeProducerEvidenceV1 {
    pub schema: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub marine_fines: Option<MarineFinesSelectionV1>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub formation_age: Option<FormationRecordingV1>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct MarineFinesSelectionV1 {
    pub law: crate::laws::MarineFinesLawV1,
    /// Maximum retained operation receipts. Exhaustion fails before stepping.
    pub receipt_limit: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct FormationRecordingV1 {
    /// Maximum retained birth/reset events. No archive pruning is implicit.
    pub event_limit: u64,
}

impl NativeProducerEvidenceV1 {
    pub fn validate(&self) -> Result<(), ParamError> {
        if self.schema != NATIVE_PRODUCER_EVIDENCE_V1
            || (self.marine_fines.is_none() && self.formation_age.is_none())
        {
            return Err(ParamError::Combination("native producer evidence requires its v1 schema and at least one selection"));
        }
        for limit in self.marine_fines.iter().map(|f| f.receipt_limit)
            .chain(self.formation_age.iter().map(|f| f.event_limit))
        {
            if limit == 0 || usize::try_from(limit).is_err() {
                return Err(ParamError::Combination("native producer evidence needs positive addressable archive limits"));
            }
        }
        if let Some(f) = &self.marine_fines { f.law.validate()?; }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PlanetParameterProfile {
    pub planet: PlanetBlock,
    pub tectonics: TectonicsBlock,
    pub mantle: MantleBlock,
    pub surface: SurfaceBlock,
    pub fidelity: FidelityBlock,
    /// Intent-space dials above raw physics (§57); neutral default.
    #[serde(default)]
    pub world_character: WorldCharacterBlock,
    /// Opt-in two-column isostasy + ocean-volume sea level (lane F
    /// FREEBOARD-SPIKE). `None` = legacy behaviour and byte-identical
    /// serialization (skipped on write).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub two_column_isostasy: Option<TwoColumnIsostasyBlock>,
    /// Opt-in continental crust bodies (lane K2 CRUST-SPIKE). `None` =
    /// legacy whole-plate continents and byte-identical serialization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub continental_bodies: Option<ContinentalBodiesBlock>,
    /// SEED-OROGEN (lane OG): opt-in orogeny / coarse-cell erosion laws.
    /// `None` = legacy laws and byte-identical serialization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orogeny: Option<OrogenyBlock>,
    /// SEED-EXTRAS D7 (lane SD): opt-in sediment laws. `None` = legacy
    /// laws and byte-identical serialization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sediment_laws: Option<SedimentLawsBlock>,
    /// CLIMATE-GENESIS (lane CG): opt-in climate forcing of genesis
    /// erosion. `None` = legacy and byte-identical serialization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub climate_forcing: Option<ClimateForcingBlock>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub native_evidence: Option<NativeProducerEvidenceV1>,
    /// WS-CONST (ADR 039 §2): OPTIONAL pin of the law set and universe
    /// constants this recipe resolves to (`universe()`). When present, the
    /// profile refuses to validate if this build resolves different physics
    /// (law ids, law parameters, constant values or units; source labels are
    /// not compared). `None` (every existing recipe) = not pinned, and the
    /// field is skipped on write so every existing params_digest is unchanged.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pinned_universe: Option<UniverseV1>,
}

impl Default for PlanetBlock {
    fn default() -> Self {
        Self {
            radius_m: 6.371e6,
            gravity_m_s2: 9.81,
            crust_density_kg_m3: 2800.0,
            mantle_density_kg_m3: 3300.0,
            sea_level_m: 0.0,
        }
    }
}

impl Default for MantleBlock {
    fn default() -> Self {
        Self {
            thermal_model: ThermalModel::Plate,
            compensation: Compensation::Flexure,
            effective_elastic_thickness_km: 40.0,
            plume_heat_k_per_step: 12.0,
        }
    }
}

impl Default for SurfaceBlock {
    fn default() -> Self {
        Self {
            rainfall_mode: RainfallMode::Uniform,
            uniform_runoff_m_per_yr: 1.0,
            erodibility: 2.0e-6,
            hillslope_diffusivity_m2_per_yr: default_hillslope_diffusivity(),
        }
    }
}

impl PlanetParameterProfile {
    /// Fail-loud dial validation (A6). Every rule names its field.
    pub fn validate(&self) -> Result<(), ParamError> {
        if let Some(e) = &self.native_evidence { e.validate()?; }
        let p = &self.planet;
        check(
            "planet.radius_m",
            p.radius_m,
            1.0e5,
            1.0e8,
            "planetary radius",
        )?;
        check(
            "planet.gravity_m_s2",
            p.gravity_m_s2,
            0.05,
            100.0,
            "surface gravity",
        )?;
        check(
            "planet.crust_density_kg_m3",
            p.crust_density_kg_m3,
            2000.0,
            3500.0,
            "crust density below mantle",
        )?;
        check(
            "planet.mantle_density_kg_m3",
            p.mantle_density_kg_m3,
            3000.0,
            6000.0,
            "mantle density",
        )?;
        if p.crust_density_kg_m3 >= p.mantle_density_kg_m3 {
            return Err(ParamError::Combination(
                "crust_density_kg_m3 must be < mantle_density_kg_m3 (stable layering)",
            ));
        }
        check("planet.sea_level_m", p.sea_level_m, -1.0e4, 1.0e4, "datum")?;

        let t = &self.tectonics;
        if !(2..=32).contains(&t.plate_count) {
            return Err(ParamError::Range {
                field: "tectonics.plate_count".into(),
                value: t.plate_count.to_string(),
                lo: "2".into(),
                hi: "32".into(),
                why: "need at least two plates for boundaries",
            });
        }
        check(
            "tectonics.rate_scale_deg_per_my",
            t.rate_scale_deg_per_my,
            0.001,
            10.0,
            "plate rate scale",
        )?;
        check(
            "tectonics.ocean_fraction",
            t.ocean_fraction,
            0.10,
            0.90,
            "ocean fraction",
        )?;
        if t.continental_seeds == 0 {
            return Err(ParamError::Range {
                field: "tectonics.continental_seeds".into(),
                value: "0".into(),
                lo: "1".into(),
                hi: "16".into(),
                why: "at least one continental seed",
            });
        }
        if t.hotspot_count == 0 {
            return Err(ParamError::Range {
                field: "tectonics.hotspot_count".into(),
                value: "0".into(),
                lo: "1".into(),
                hi: "8".into(),
                why: "at least one hotspot (K6)",
            });
        }
        check(
            "tectonics.duration_my",
            t.duration_my,
            1.0,
            2000.0,
            "history length",
        )?;
        check("tectonics.step_my", t.step_my, 0.05, 50.0, "time step")?;
        if t.step_my > t.duration_my {
            return Err(ParamError::Combination(
                "tectonics.step_my > duration_my: empty history",
            ));
        }
        let steps = t.duration_my / t.step_my;
        if (steps - steps.round()).abs() > 1.0e-9 {
            return Err(ParamError::Combination(
                "tectonics.duration_my must be an integral multiple of step_my \
                 (fail-loud; no implicit partial step)",
            ));
        }

        let m = &self.mantle;
        if m.compensation == Compensation::Flexure {
            check(
                "mantle.effective_elastic_thickness_km",
                m.effective_elastic_thickness_km,
                5.0,
                150.0,
                "flexure requires an effective elastic thickness",
            )?;
        }
        check(
            "mantle.plume_heat_k_per_step",
            m.plume_heat_k_per_step,
            0.0,
            500.0,
            "plume source term",
        )?;

        let s = &self.surface;
        check(
            "surface.uniform_runoff_m_per_yr",
            s.uniform_runoff_m_per_yr,
            0.0,
            50.0,
            "rainfall",
        )?;
        check(
            "surface.erodibility",
            s.erodibility,
            0.0,
            1.0e-3,
            "erodibility",
        )?;
        check(
            "surface.hillslope_diffusivity_m2_per_yr",
            s.hillslope_diffusivity_m2_per_yr,
            0.0,
            1.0e3,
            "hillslope diffusivity",
        )?;

        check(
            "world_character.waterfall_wonder",
            self.world_character.waterfall_wonder,
            0.0,
            1.0,
            "intent dial range (spec §57)",
        )?;
        let bias = self.world_character.resolve();
        check(
            "world_character.bias.uplift_band_sharpness",
            bias.uplift_band_sharpness,
            1.0,
            2.5,
            "resolved channel",
        )?;
        check(
            "world_character.bias.uplift_step_persistence",
            bias.uplift_step_persistence,
            1.0,
            2.0,
            "resolved channel (cause, never detector)",
        )?;
        check(
            "world_character.bias.resistance_contrast_magnitude",
            bias.resistance_contrast_magnitude,
            1.1,
            2.5,
            "resolved channel (no zebra cake)",
        )?;

        let f = &self.fidelity;
        if f.grid_level > 3 {
            return Err(ParamError::Range {
                field: "fidelity.grid_level".into(),
                value: f.grid_level.to_string(),
                lo: "0".into(),
                hi: "3".into(),
                why: "grid level cap (2304x1152 at L3)",
            });
        }
        if !(1..=24).contains(&f.radial_layers) {
            return Err(ParamError::Range {
                field: "fidelity.radial_layers".into(),
                value: f.radial_layers.to_string(),
                lo: "1".into(),
                hi: "24".into(),
                why: "mantle radial layers",
            });
        }
        if f.substeps == 0 || f.substeps > 64 {
            return Err(ParamError::Range {
                field: "fidelity.substeps".into(),
                value: f.substeps.to_string(),
                lo: "1".into(),
                hi: "64".into(),
                why: "declared numerical substeps",
            });
        }
        if let Some(b) = &self.two_column_isostasy {
            check(
                "two_column_isostasy.oceanic_crust_density_kg_m3",
                b.oceanic_crust_density_kg_m3,
                2000.0,
                3500.0,
                "oceanic crust density",
            )?;
            if b.oceanic_crust_density_kg_m3 >= p.mantle_density_kg_m3 {
                return Err(ParamError::Combination(
                    "two_column_isostasy.oceanic_crust_density_kg_m3 must be < \
                     planet.mantle_density_kg_m3 (stable layering)",
                ));
            }
            check(
                "two_column_isostasy.seawater_density_kg_m3",
                b.seawater_density_kg_m3,
                900.0,
                1300.0,
                "liquid water density",
            )?;
            if b.seawater_density_kg_m3 >= p.crust_density_kg_m3
                || b.seawater_density_kg_m3 >= b.oceanic_crust_density_kg_m3
            {
                return Err(ParamError::Combination(
                    "two_column_isostasy.seawater_density_kg_m3 must be below both crust densities",
                ));
            }
            check(
                "two_column_isostasy.compensation_depth_m",
                b.compensation_depth_m,
                1.0e4,
                1.0e6,
                "compensation depth below the mantle free-surface datum",
            )?;
            check(
                "two_column_isostasy.ocean_water_volume_m3",
                b.ocean_water_volume_m3,
                0.0,
                1.0e22,
                "liquid ocean water volume",
            )?;
            if let Some(m) = b.compensation_depth_margin_m {
                check(
                    "two_column_isostasy.compensation_depth_margin_m",
                    m,
                    1.0,
                    1.0e6,
                    "margin of the derived compensation depth below the deepest root",
                )?;
            }
        }
        if let Some(b) = &self.continental_bodies {
            if let Some(c) = b.continental_area_fraction {
                check(
                    "continental_bodies.continental_area_fraction",
                    c,
                    0.01,
                    0.95,
                    "continental crust area fraction",
                )?;
            }
            check(
                "continental_bodies.resolved_area_fraction",
                b.resolved_area_fraction(t),
                0.01,
                0.95,
                "continental crust area fraction in force (declared or 1 - ocean_fraction)",
            )?;
            check(
                "continental_bodies.rifted_margin_width_m",
                b.rifted_margin_width_m,
                0.0,
                2.0e6,
                "rifted margin width",
            )?;
        }
        // WS-CONST: a pinned universe must still be what this build resolves.
        // Unpinned profiles are resolved and checked by their consumers via
        // `universe()`; validate() stays cheap because the legacy isostasy
        // helpers call it per evaluation.
        // SEED-OROGEN: an orogeny block is checked against the code-owned
        // envelopes of its law constants (typed UniverseError).
        if self.pinned_universe.is_some() || self.orogeny.is_some() || self.sediment_laws.is_some()
            || self.native_evidence.is_some()
        {
            self.universe()?;
        }
        Ok(())
    }

    /// WS-CONST: the law set and universe constants this profile runs
    /// under, resolved from the v0 blocks (ADR 039 §2 legacy mapping),
    /// validated (typed `UniverseError`), and checked against
    /// `pinned_universe` when one is present. Pure and allocation-light; the
    /// generator reads its physical constants through this.
    pub fn universe(&self) -> Result<UniverseV1, UniverseError> {
        let u = self.resolve_universe()?;
        u.validate()?;
        if let Some(pin) = &self.pinned_universe {
            pin.validate()?;
            pin.same_physics(&u)?;
        }
        Ok(u)
    }

    /// Pin the current resolution into the recipe (moves params_digest:
    /// a pinned recipe is a new recipe). Replaces any existing pin.
    pub fn pin_universe(&mut self) -> Result<(), UniverseError> {
        self.pinned_universe = None;
        let u = self.universe()?;
        self.pinned_universe = Some(u);
        Ok(())
    }

    /// The legacy mapping (no validation): v0 blocks + lifted code constants.
    fn resolve_universe(&self) -> Result<UniverseV1, UniverseError> {
        let p = &self.planet;
        let s = &self.surface;
        let d = PlanetBlock::default();
        let sd = SurfaceBlock::default();
        let leg = |o: &'static str| SourceClassV1::LegacyPresetValue {
            origin: Cow::Borrowed(o),
        };
        let pa = |b: &'static str| SourceClassV1::ProjectAssumptionAwaitingReview {
            basis: Cow::Borrowed(b),
        };
        let cited = |c: &'static str| SourceClassV1::Cited {
            citation: Cow::Borrowed(c),
        };
        let derived = |f: &'static str| SourceClassV1::Derived {
            from: Cow::Borrowed(f),
        };
        // A value equal to the shipped default keeps the default's label;
        // any other value was declared by the recipe and has no citation.
        let label = |v: f64, default: f64, if_default: SourceClassV1, declared: &'static str| {
            if v.to_bits() == default.to_bits() {
                if_default
            } else {
                pa(declared)
            }
        };
        let c = ConstantV1::new;
        let two = self.two_column_isostasy.as_ref();
        let (isostasy, sea_level) = match two {
            None => (
                IsostasyLawV1::SingleColumnReferenceV0 {
                    reference_crust_thickness_m: c(
                        legacy::REFERENCE_CRUST_THICKNESS_M,
                        UnitV1::Meter,
                        leg("pf-geodynamics isostasy.rs REFERENCE_CRUST_M: the Airy reference column at elevation 0 (the D1 cause, ADR 039 Problem)"),
                    ),
                    compensation: self.mantle.compensation,
                },
                SeaLevelLawV1::PrescribedDatumV0,
            ),
            Some(b) => {
                if self.mantle.compensation != Compensation::Airy {
                    return Err(UniverseError::Unsupported {
                        law: "isostasy/two-column-common-depth-v0 with mantle.compensation = Flexure",
                        why: "the two-column law is local (Airy) compensation; flexure is not built for it",
                    });
                }
                (
                    IsostasyLawV1::TwoColumnCommonDepthV0 {
                        fixed_compensation_depth_m: c(
                            b.compensation_depth_m,
                            UnitV1::Meter,
                            pa("fixed-depth mode compensation depth; elevations do not depend on it while it lies below every root (A-COMPENSATION-DEPTH)"),
                        ),
                        derived_depth_margin_m: b.compensation_depth_margin_m.map(|m| {
                            c(
                                m,
                                UnitV1::Meter,
                                pa("validity bound below the deepest root; the depth cancels in the mass balance (A-COMPENSATION-DEPTH-MARGIN, ADR 039 C4/C5)"),
                            )
                        }),
                    },
                    SeaLevelLawV1::VolumePrimaryGlobalFillV0,
                )
            }
        };
        let crust_seeding = match &self.continental_bodies {
            None => CrustSeedingLawV1::WholePlateContinentsV0,
            Some(b) => CrustSeedingLawV1::GeodesicDiscBodiesV0 {
                continental_area_fraction: c(
                    b.resolved_area_fraction(&self.tectonics),
                    UnitV1::Dimensionless,
                    if b.continental_area_fraction.is_none() {
                        derived("1 - tectonics.ocean_fraction (the scenario's declared ocean fraction read as the continental crust area it implies; A-CRUST-BODY-AREA-FRACTION)")
                    } else {
                        pa("declared in continental_bodies.continental_area_fraction")
                    },
                ),
                rifted_margin_width_m: c(
                    b.rifted_margin_width_m,
                    UnitV1::Meter,
                    label(
                        b.rifted_margin_width_m,
                        ContinentalBodiesBlock::v0().rifted_margin_width_m,
                        pa("recalled from Peron-Pinvidic and Manatschal 2009, Int J Earth Sci 98:1581; citation not checked (A-RIFTED-MARGIN-WIDTH)"),
                        "declared in continental_bodies.rifted_margin_width_m; no citation",
                    ),
                ),
            },
        };
        // SEED-OROGEN (lane OG): law entries from the optional block; absent
        // = the legacy variants (the thickening entry is then skipped on
        // write, so legacy law sets serialize exactly as before).
        let orogeny = self.orogeny.as_ref();
        let tectonic_thickening = match orogeny.and_then(|o| o.rate_limited_belt.as_ref()) {
            None => TectonicThickeningLawV1::KinematicOverlapStackingV0,
            Some(b) => {
                let v1 = RateLimitedBeltBlock::v1();
                TectonicThickeningLawV1::RateLimitedBeltV1 {
                    belt_shear_zone_dip_deg: c(
                        b.belt_shear_zone_dip_deg,
                        UnitV1::Degree,
                        label(b.belt_shear_zone_dip_deg, v1.belt_shear_zone_dip_deg, pa("recalled: crustal-scale two-sided orogenic wedges grow above conjugate step-up shear zones rooted at the mantle singularity point, dipping ~30 deg (Coulomb, 45 - phi/2 with phi ~ 30 deg; Willett, Beaumont and Fullsack 1993, Geology 21:371); belt width W = 2 h cot(dip); citation not checked (A-BELT-SHEAR-DIP)"), "declared in orogeny.rate_limited_belt.belt_shear_zone_dip_deg; no citation"),
                    ),
                    lower_crust_channel_viscosity_pa_s: c(
                        b.lower_crust_channel_viscosity_pa_s,
                        UnitV1::PascalSecond,
                        label(b.lower_crust_channel_viscosity_pa_s, v1.lower_crust_channel_viscosity_pa_s, pa("recalled: lower-crustal channel-flow models of Tibet use ~1e18 Pa s (Clark and Royden 2000, Geology 28:703) to 1e19-1e21 Pa s (Beaumont et al. 2001, Nature 414:738; Bird 1991, JGR 96:10275); 1e19 chosen a priori mid-range; citations not checked (A-CHANNEL-VISCOSITY)"), "declared in orogeny.rate_limited_belt.lower_crust_channel_viscosity_pa_s; no citation"),
                    ),
                    channel_reference_thickness_m: c(
                        b.channel_reference_thickness_m,
                        UnitV1::Meter,
                        label(b.channel_reference_thickness_m, v1.channel_reference_thickness_m, pa("project assumption: crust of the 35 km reference column does not flow laterally on these timescales (stable continental interiors); only the excess thickness of thickened crust forms the ductile channel (A-CHANNEL-REFERENCE-THICKNESS)"), "declared in orogeny.rate_limited_belt.channel_reference_thickness_m; no citation"),
                    ),
                    delamination_thickness_m: c(
                        b.delamination_thickness_m,
                        UnitV1::Meter,
                        label(b.delamination_thickness_m, v1.delamination_thickness_m, pa("backstop at the legacy 70 km cap (A-UPLIFT-COLLISION-CAP, Tibetan crust 65-75 km, Zhang and Klemperer 2005, JGR 110 B01409); excess root removal read as eclogitised lower-crust delamination (Kay and Kay 1993, Tectonophysics 219:177); citations not checked; its share of belt input is reported (falsifier F1) (A-DELAMINATION-THICKNESS)"), "declared in orogeny.rate_limited_belt.delamination_thickness_m; no citation"),
                    ),
                    arc_productivity_m: c(
                        b.arc_productivity_m,
                        UnitV1::Meter,
                        label(b.arc_productivity_m, v1.arc_productivity_m, derived("legacy A-UPLIFT-ARC intent restated per subducted area: 0.005 x ~100 km arc width = 500 m of arc crust per unit area of subducted plate, i.e. 25 km3/km/My at 50 km/My (order of Jicha and Jagoutz 2015, Elements 11:105, recalled, not checked) (A-ARC-PRODUCTIVITY)"),
                        "declared in orogeny.rate_limited_belt.arc_productivity_m; no citation"),
                    ),
                    arc_belt_width_m: c(
                        b.arc_belt_width_m,
                        UnitV1::Meter,
                        label(b.arc_belt_width_m, v1.arc_belt_width_m, pa("project assumption: the ~100 km arc width named in the legacy A-UPLIFT-ARC note; arc-trench gap not resolved (A-ARC-BELT-WIDTH)"), "declared in orogeny.rate_limited_belt.arc_belt_width_m; no citation"),
                    ),
                }
            }
        };
        // SEED-CARRY-CC (lane CC): the plate carry entry (absent = legacy,
        // skipped on write).
        let plate_carry = match orogeny.and_then(|o| o.rate_limited_belt.as_ref()) {
            // SEED-CARRY-CONSERVE (lane CV): one carry law per world.
            Some(b) if b.front_rescue_carry && b.continental_parcel_carry.is_some() => {
                return Err(UniverseError::Unsupported {
                    law: "plate-carry/front-rescue-carry-v1 with plate-carry/continental-parcel-carry-v1",
                    why: "one plate carry law per world; continental-parcel-carry-v1 replaces the front rescue",
                });
            }
            Some(b) if b.front_rescue_carry => PlateCarryLawV1::FrontRescueCarryV1,
            Some(b) if b.continental_parcel_carry.is_some() => {
                let k = b.continental_parcel_carry.as_ref().expect("checked");
                let v1 = ContinentalParcelCarryBlock::v1();
                const D: &str = "declared in orogeny.rate_limited_belt.continental_parcel_carry; no citation";
                PlateCarryLawV1::ContinentalParcelCarryV1 {
                    parcel_overlap_cell_fraction: c(k.parcel_overlap_cell_fraction, UnitV1::Dimensionless, label(k.parcel_overlap_cell_fraction, v1.parcel_overlap_cell_fraction, pa("project assumption, grid-scale numerical choice: collision consumption starts once two continents overlap by half a cell; touching coasts share boundary cells through the one-row tent kernel without overlapping (lane CV PREREG.txt) (A-PARCEL-OVERLAP)"), D)),
                    plate_fragment_min_area_fraction: c(k.plate_fragment_min_area_fraction, UnitV1::Dimensionless, label(k.plate_fragment_min_area_fraction, v1.plate_fragment_min_area_fraction, pa("project assumption: 0.1 % of the sphere (0.51 Mkm2 on Earth). Recalled (Bird 2003, G-cubed 4:1027, not checked): Earth's smaller named plates are Juan de Fuca ~0.25, Scotia ~1.6, Cocos ~2.9 Mkm2; the toy seeds ~7 plates of tens of Mkm2, so smaller plate-id pieces are carry speckle (A-PLATE-FRAGMENT-FOLD)"), D)),
                    write_back_residual_iterations: k.write_back_residual_iterations.map(|v| c(v, UnitV1::Dimensionless, pa("project assumption (lane HYPSO-LEVERS 2026-10-21): numerical correction count; 3 passes cut the write-back mismatch of crust excess above 36 km from -0.26 to -0.05 Mkm3/My on seed 42 (A-PARCEL-WRITE-BACK-RESIDUAL)"))),
                }
            }
            _ => PlateCarryLawV1::SemiLagrangianNearestMaterialV0,
        };
        // SEED-CARRY-CONSERVE (lane CV): the suture entry (absent = legacy,
        // skipped on write).
        let plate_suture = match orogeny
            .and_then(|o| o.rate_limited_belt.as_ref())
            .and_then(|b| b.plate_suture.as_ref())
        {
            None => PlateSutureLawV1::NoneV0,
            Some(s) if s.indentation.is_none() => PlateSutureLawV1::ContactLengthPoleShareV1 {
                suture_contact_length_m: c(s.suture_contact_length_m, UnitV1::Meter, label(s.suture_contact_length_m, PlateSutureBlock::v1().suture_contact_length_m, pa("project assumption: recalled Himalayan collision front ~2400-2500 km, Alpine ~1200 km (not checked); on the toy's ~3900 km-radius disc continents a 2500 km chord is ~400 km of convergence after first contact, so belts still form before the suture (lane CV PREREG.txt) (A-SUTURE-CONTACT-LENGTH)"), "declared in orogeny.rate_limited_belt.plate_suture; no citation")),
            },
            // SEED-RANGES (lane RG): the indentation suture law.
            Some(s) => {
                let ind = s.indentation.as_ref().expect("indentation block");
                let v1 = SutureIndentationBlock::v1();
                const D: &str = "declared in orogeny.rate_limited_belt.plate_suture.indentation; no citation";
                PlateSutureLawV1::ContactLengthIndentationTaperV1 {
                    suture_contact_length_m: c(s.suture_contact_length_m, UnitV1::Meter, label(s.suture_contact_length_m, PlateSutureBlock::v1().suture_contact_length_m, pa("project assumption: recalled Himalayan collision front ~2400-2500 km, Alpine ~1200 km (not checked); on the toy's ~3900 km-radius disc continents a 2500 km chord is ~400 km of convergence after first contact, so belts still form before the suture (lane CV PREREG.txt) (A-SUTURE-CONTACT-LENGTH)"), "declared in orogeny.rate_limited_belt.plate_suture; no citation")),
                    indentation_time_my: c(ind.indentation_time_my, UnitV1::Megayear, label(ind.indentation_time_my, v1.indentation_time_my, pa("project assumption: convergence continues for tens of My after continental collision; recalled India-Asia collision ~50 Ma and still converging at ~4-5 cm/yr, Alpine collision ~35-30 Ma with slowing convergence (not checked); design review review 2026-10-13 'Earth 40-50; try 40' (lane RG PREREG.txt) (A-INDENTATION-TIME)"), D)),
                    indentation_shortening_budget_m: c(ind.indentation_shortening_budget_m, UnitV1::Meter, label(ind.indentation_shortening_budget_m, v1.indentation_shortening_budget_m, pa("project assumption: post-collision shortening recalled ~500 km in the Alps and ~2000 km for India-Asia including subducted Greater India (the extreme case; not checked); midpoint of the design review review's 500-1000 km (lane RG PREREG.txt) (A-INDENTATION-BUDGET)"), D)),
                }
            }
        };
        // SEED-RANGES (lane RG): the plateau entry (absent = legacy, skipped
        // on write).
        let orogen_plateau = match orogeny
            .and_then(|o| o.rate_limited_belt.as_ref())
            .and_then(|b| b.plateau.as_ref())
        {
            None => OrogenPlateauLawV1::NoneV0,
            Some(p) => {
                let v1 = OrogenPlateauBlock::v1();
                const D: &str = "declared in orogeny.rate_limited_belt.plateau; no citation";
                let built = OrogenPlateauLawV1::PlateauCapSpillV1 {
                    plateau_crust_cap_m: c(p.plateau_crust_cap_m, UnitV1::Meter, label(p.plateau_crust_cap_m, v1.plateau_crust_cap_m, pa("project assumption set by plateau HEIGHT: Tibetan crust ~65-75 km (Zhang and Klemperer 2005) stands at ~5 km, and plateau height is limited near ~5 km by the force balance of the range front (Molnar and Lyon-Caen 1988); recalled, not checked. This model's Airy column has no lithospheric root term and puts 60-65 km crust at ~5.2 km and 70 km at ~6.0 km (measured, lane RG PREREG.txt), so 65 km (A-PLATEAU-CAP)"), D)),
                    plateau_spill_reach_m: c(p.plateau_spill_reach_m, UnitV1::Meter, label(p.plateau_spill_reach_m, v1.plateau_spill_reach_m, pa("project assumption: the Tibetan plateau is ~1000 km N-S by ~2500 km E-W, the Altiplano 300-400 km wide (recalled, not checked); thickened crust farther than this from its source is delaminated (counted) (A-PLATEAU-REACH)"), D)),
                };
                // SEED-WORLD-COMPOSE (lane CP): law_version 2.
                match (p.law_version, p.plateau_foothill_gradient, p.plateau_taper_floor_m) {
                    (1, None, None) => built,
                    (2, Some(g), Some(f)) => {
                        let v2 = OrogenPlateauBlock::v2();
                        let OrogenPlateauLawV1::PlateauCapSpillV1 { plateau_crust_cap_m, plateau_spill_reach_m } = built else { unreachable!() };
                        OrogenPlateauLawV1::PlateauCapSpillV2 {
                            plateau_crust_cap_m,
                            plateau_spill_reach_m,
                            plateau_foothill_gradient: c(g, UnitV1::Dimensionless, label(g, v2.plateau_foothill_gradient.unwrap_or(f64::NAN), pa("project assumption (lane CP PREREG.txt): plateau margins step down over a few hundred km; recalled Tibet's crust 70 -> 40 km over ~200-300 km at its eastern and northern margins and the Subandean ramp of the Altiplano ~200-300 km wide (not checked), i.e. ~10 km of crust per 100 km (A-PLATEAU-FOOTHILL-GRADIENT)"), D)),
                            plateau_taper_floor_m: c(f, UnitV1::Meter, label(f, v2.plateau_taper_floor_m.unwrap_or(f64::NAN), pa("project assumption (lane CP PREREG.txt): the taper acts only where the thicker cell of a pair is above 40 km, i.e. on thickened crust, so normal crust and rifted margins are untouched (A-PLATEAU-TAPER-FLOOR)"), D)),
                        }
                    }
                    _ => {
                        return Err(UniverseError::Unsupported {
                            law: "orogen-plateau/plateau-cap-spill law_version other than 1 (no v2 constants) or 2 (both v2 constants)",
                            why: "only plateau-cap-spill-v1 and plateau-cap-spill-v2 exist",
                        });
                    }
                }
            }
        };
        // SEED-STRATIGRAPHY (lane ST): the parcel-stacks entry (absent =
        // legacy, skipped on write).
        let stratigraphy = match orogeny
            .and_then(|o| o.rate_limited_belt.as_ref())
            .and_then(|b| b.stratigraphy.as_ref())
        {
            None => StratigraphyLawV1::NoneV0,
            Some(g) => {
                let v1 = StratigraphyBlock::v1();
                const D: &str = "declared in orogeny.rate_limited_belt.stratigraphy; no citation";
                let law = StratigraphyLawV1::ParcelStacksV1 {
                    strata_max_layers: c(g.strata_max_layers, UnitV1::Dimensionless, label(g.strata_max_layers, v1.strata_max_layers, pa("design review stratigraphy consult 2026-10-13 section 1: hard cap of 12 layers (basement included) bounds memory; enough for a game (A-STRATA-CAPACITY)"), D)),
                    strata_merge_gap_my: c(g.strata_merge_gap_my, UnitV1::Megayear, label(g.strata_merge_gap_my, v1.strata_merge_gap_my, pa("design review consult section 1: same lithology and environment appended within 2 My extends the top layer (one step at dt 2 My) (A-STRATA-MERGE)"), D)),
                    strata_unconformity_gap_my: c(g.strata_unconformity_gap_my, UnitV1::Megayear, label(g.strata_unconformity_gap_my, v1.strata_unconformity_gap_my, pa("design review consult section 1: an age gap above 5 My between layers reads as an unconformity (A-STRATA-UNCONFORMITY)"), D)),
                    strata_min_layer_m: c(g.strata_min_layer_m, UnitV1::Meter, label(g.strata_min_layer_m, v1.strata_min_layer_m, pa("design review consult section 1: cover layers thinner than 10 m after erosion merge downward (A-STRATA-MERGE)"), D)),
                    strata_trace_m: c(g.strata_trace_m, UnitV1::Meter, label(g.strata_trace_m, v1.strata_trace_m, pa("project assumption: additions under 1 cm extend the layer they land on instead of opening a layer (numerical floor) (A-STRATA-MERGE)"), D)),
                    strata_arc_extrusive_fraction: c(g.strata_arc_extrusive_fraction, UnitV1::Dimensionless, label(g.strata_arc_extrusive_fraction, v1.strata_arc_extrusive_fraction, pa("design review consult section 2: 30 % of arc addition erupts as volcanics, 70 % stays as plutons (recalled intrusive:extrusive ratios of arcs ~2:1 to 10:1; not checked) (A-STRATA-ARC)"), D)),
                    strata_plume_extrusive_fraction: c(g.strata_plume_extrusive_fraction, UnitV1::Dimensionless, label(g.strata_plume_extrusive_fraction, v1.strata_plume_extrusive_fraction, pa("derived: 1 / (1 + r) with lane FB's flood_basalt_intrusive_ratio r = 4 (LIP extrusive share); hotspot tail crust uses the same split (A-STRATA-PLUME)"), D)),
                    strata_shortening_pluton_fraction: c(g.strata_shortening_pluton_fraction, UnitV1::Dimensionless, label(g.strata_shortening_pluton_fraction, v1.strata_shortening_pluton_fraction, pa("design review consult section 2: a declared 10 % of collisional thickening becomes granite plutons (crustal melting in thickened orogens; not checked) (A-STRATA-COLLISION)"), D)),
                    strata_pluton_roof_m: c(g.strata_pluton_roof_m, UnitV1::Meter, label(g.strata_pluton_roof_m, v1.strata_pluton_roof_m, pa("PREREG revision 1: plutons are emplaced below a crystalline roof (recalled emplacement depths of granite plutons ~2-10 km; 4 km chosen; not checked), so a pluton reaches the surface only after its roof is eroded (A-STRATA-ARC)"), D)),
                    strata_metamorphic_depth_m: c(g.strata_metamorphic_depth_m, UnitV1::Meter, label(g.strata_metamorphic_depth_m, v1.strata_metamorphic_depth_m, pa("design review consult section 2: cover buried below 10 km becomes metamorphic (greenschist onset ~300 C at ~10 km on a ~30 C/km geotherm; recalled, not checked) (A-STRATA-METAMORPHISM)"), D)),
                    strata_high_grade_depth_m: c(g.strata_high_grade_depth_m, UnitV1::Meter, label(g.strata_high_grade_depth_m, v1.strata_high_grade_depth_m, pa("project assumption: amphibolite and above below ~20 km (~600 C; recalled, not checked) (A-STRATA-METAMORPHISM)"), D)),
                    strata_shelf_depth_m: c(g.strata_shelf_depth_m, UnitV1::Meter, label(g.strata_shelf_depth_m, v1.strata_shelf_depth_m, pa("design review consult section 2: deeper than 200 m (shelf break) a marine deposit is deep-marine shale (recalled typical shelf-break depth ~130-200 m; not checked) (A-STRATA-ENVIRONMENT)"), D)),
                    strata_carbonate_max_paleolat_deg: c(g.strata_carbonate_max_paleolat_deg, UnitV1::Degree, label(g.strata_carbonate_max_paleolat_deg, v1.strata_carbonate_max_paleolat_deg, pa("LATITUDE FALLBACK (A-STRATA-LATITUDE-CLIMATE): no climate law runs (uniform rainfall); modern shallow-water carbonate platforms lie mostly within ~30 deg of the equator (recalled, not checked)"), D)),
                    strata_carbonate_max_supply_m_per_my: c(g.strata_carbonate_max_supply_m_per_my, UnitV1::MeterPerMegayear, label(g.strata_carbonate_max_supply_m_per_my, v1.strata_carbonate_max_supply_m_per_my, pa("project assumption: carbonate factories are smothered by clastic supply; a shelf cell receiving more than 50 m/My of sediment is clastic (A-STRATA-ENVIRONMENT)"), D)),
                    strata_arid_lat_min_deg: c(g.strata_arid_lat_min_deg, UnitV1::Degree, label(g.strata_arid_lat_min_deg, v1.strata_arid_lat_min_deg, pa("LATITUDE FALLBACK (A-STRATA-LATITUDE-CLIMATE): the subtropical arid belt ~15-35 deg (Hadley descent; recalled)"), D)),
                    strata_arid_lat_max_deg: c(g.strata_arid_lat_max_deg, UnitV1::Degree, label(g.strata_arid_lat_max_deg, v1.strata_arid_lat_max_deg, pa("LATITUDE FALLBACK (A-STRATA-LATITUDE-CLIMATE): as strata_arid_lat_min_deg"), D)),
                    strata_humid_lat_max_deg: c(g.strata_humid_lat_max_deg, UnitV1::Degree, label(g.strata_humid_lat_max_deg, v1.strata_humid_lat_max_deg, pa("LATITUDE FALLBACK (A-STRATA-LATITUDE-CLIMATE): humid tropics within ~15 deg (ITCZ; recalled)"), D)),
                    strata_humid_midlat_min_deg: c(g.strata_humid_midlat_min_deg, UnitV1::Degree, label(g.strata_humid_midlat_min_deg, v1.strata_humid_midlat_min_deg, pa("LATITUDE FALLBACK (A-STRATA-LATITUDE-CLIMATE): humid mid-latitude storm belt ~40-65 deg (recalled)"), D)),
                    strata_humid_midlat_max_deg: c(g.strata_humid_midlat_max_deg, UnitV1::Degree, label(g.strata_humid_midlat_max_deg, v1.strata_humid_midlat_max_deg, pa("LATITUDE FALLBACK (A-STRATA-LATITUDE-CLIMATE): as strata_humid_midlat_min_deg"), D)),
                    strata_coal_max_elevation_m: c(g.strata_coal_max_elevation_m, UnitV1::Meter, label(g.strata_coal_max_elevation_m, v1.strata_coal_max_elevation_m, pa("project assumption: coal swamps are humid coastal plains and delta tops within ~200 m of sea level (A-STRATA-ENVIRONMENT)"), D)),
                    strata_conglomerate_relief_m: c(g.strata_conglomerate_relief_m, UnitV1::Meter, label(g.strata_conglomerate_relief_m, v1.strata_conglomerate_relief_m, pa("design review consult section 2: land deposits where the local relief within one cell exceeds 1 km are conglomerate (A-STRATA-ENVIRONMENT)"), D)),
                    strata_river_mouth_area_m2: c(g.strata_river_mouth_area_m2, UnitV1::SquareMeter, label(g.strata_river_mouth_area_m2, v1.strata_river_mouth_area_m2, pa("project assumption: a marine cell fed directly by a land neighbour draining at least 1e5 km2 is a delta (deltaic sandstone) (A-STRATA-ENVIRONMENT)"), D)),
                    strata_publish_depth_m: c(g.strata_publish_depth_m, UnitV1::Meter, label(g.strata_publish_depth_m, v1.strata_publish_depth_m, pa("design review consult section 2: the published record lists the column to 3 km (deeper layers never reach a T1b surface) (A-STRATA-PUBLISH)"), D)),
                    strata_basement_min_m: c(g.strata_basement_min_m, UnitV1::Meter, label(g.strata_basement_min_m, v1.strata_basement_min_m, pa("design review consult section 1: parcels whose basement falls below 5 km are counted (count only, never enforced) (A-STRATA-CAPACITY)"), D)),
                };
                // SEED-STRATIGRAPHY stage A.1 (lane SA): parcel-stacks-v2 when the
                // v2 constants are present.
                match &g.v2 {
                    None => law,
                    Some(x) => {
                        let w = StratigraphyV2Block::v2();
                        let vb = StratigraphyBlock::v2();
                        match law {
                            StratigraphyLawV1::ParcelStacksV1 {
                            strata_max_layers,
                            strata_merge_gap_my,
                            strata_unconformity_gap_my,
                            strata_min_layer_m,
                            strata_trace_m,
                            strata_arc_extrusive_fraction,
                            strata_plume_extrusive_fraction,
                            strata_shortening_pluton_fraction,
                            strata_pluton_roof_m,
                            strata_metamorphic_depth_m,
                            strata_high_grade_depth_m,
                            strata_shelf_depth_m,
                            strata_carbonate_max_paleolat_deg,
                            strata_carbonate_max_supply_m_per_my,
                            strata_arid_lat_min_deg,
                            strata_arid_lat_max_deg,
                            strata_humid_lat_max_deg,
                            strata_humid_midlat_min_deg,
                            strata_humid_midlat_max_deg,
                            strata_coal_max_elevation_m,
                            strata_conglomerate_relief_m,
                            strata_river_mouth_area_m2,
                            strata_publish_depth_m,
                            strata_basement_min_m,
                            } => StratigraphyLawV1::ParcelStacksV2 {
                                strata_max_layers: if g.strata_max_layers == vb.strata_max_layers {
                                    c(g.strata_max_layers, UnitV1::Dimensionless, pa("lane SA PREREG 2026-10-13 (C4): cap raised 12 -> 16 layers (basement included) after half the v1 stacks sat at the cap by 1000 My; ~40 B per layer (A-STRATA-CAPACITY)"))
                                } else {
                                    strata_max_layers
                                },
                                strata_merge_gap_my,
                                strata_unconformity_gap_my,
                                strata_min_layer_m,
                                strata_trace_m,
                                strata_arc_extrusive_fraction,
                                strata_plume_extrusive_fraction,
                                strata_shortening_pluton_fraction,
                                strata_pluton_roof_m,
                                strata_metamorphic_depth_m,
                                strata_high_grade_depth_m,
                                strata_shelf_depth_m,
                                strata_carbonate_max_paleolat_deg,
                                strata_carbonate_max_supply_m_per_my,
                                strata_arid_lat_min_deg,
                                strata_arid_lat_max_deg,
                                strata_humid_lat_max_deg,
                                strata_humid_midlat_min_deg,
                                strata_humid_midlat_max_deg,
                                strata_coal_max_elevation_m,
                                strata_conglomerate_relief_m,
                                strata_river_mouth_area_m2,
                                strata_publish_depth_m,
                                strata_basement_min_m,
                                strata_shelf_width_m: c(x.strata_shelf_width_m, UnitV1::Meter, label(x.strata_shelf_width_m, w.strata_shelf_width_m, pa("recalled global mean continental-shelf width ~80 km (not checked); the declared sub-cell shelf-slope ramp: a flat shelf of this width seaward of a coast (A-STRATA-SUBCELL-SHELF)"), D)),
                                strata_slope_foot_depth_m: c(x.strata_slope_foot_depth_m, UnitV1::Meter, label(x.strata_slope_foot_depth_m, w.strata_slope_foot_depth_m, pa("recalled continental rise / abyssal plain depth ~4 km (not checked): the non-shelf part of a marine cell lies no deeper, which bounds the shelf share of a deep cell (A-STRATA-SUBCELL-SHELF)"), D)),
                                strata_shelf_majority_fraction: c(x.strata_shelf_majority_fraction, UnitV1::Dimensionless, label(x.strata_shelf_majority_fraction, w.strata_shelf_majority_fraction, pa("project assumption: a marine cell's deposit takes the shelf facies when at least half of its area is shelf (majority facies, one class per cell per step) (A-STRATA-SUBCELL-SHELF)"), D)),
                                strata_platform_cover_max_m: c(x.strata_platform_cover_max_m, UnitV1::Meter, label(x.strata_platform_cover_max_m, w.strata_platform_cover_max_m, pa("brief 2026-10-13 (0.5-2 km); recalled cratonic platform covers ~1-3 km in intracratonic basins, absent on shields (not checked) (A-STRATA-PLATFORM-COVER)"), D)),
                                strata_platform_full_elev_m: c(x.strata_platform_full_elev_m, UnitV1::Meter, label(x.strata_platform_full_elev_m, w.strata_platform_full_elev_m, pa("project assumption: parcels at or below sea level at t = 0 carry the full platform cover (A-STRATA-PLATFORM-COVER)"), D)),
                                strata_platform_zero_elev_m: c(x.strata_platform_zero_elev_m, UnitV1::Meter, label(x.strata_platform_zero_elev_m, w.strata_platform_zero_elev_m, pa("project assumption: about the 80th percentile of t = 0 land elevation (L1 s42 609 m, L0 s42 577 m): the highest fifth of the old land (suture uplands, craton highs) is exposed shield (A-STRATA-PLATFORM-COVER)"), D)),
                                strata_platform_sandstone_fraction: c(x.strata_platform_sandstone_fraction, UnitV1::Dimensionless, label(x.strata_platform_sandstone_fraction, w.strata_platform_sandstone_fraction, pa("project assumption: transgressive platform sequence, basal sandstone 30 % (A-STRATA-PLATFORM-COVER)"), D)),
                                strata_platform_carbonate_fraction: c(x.strata_platform_carbonate_fraction, UnitV1::Dimensionless, label(x.strata_platform_carbonate_fraction, w.strata_platform_carbonate_fraction, pa("project assumption: transgressive platform sequence, carbonate cap 30 % (shale 40 %) (A-STRATA-PLATFORM-COVER)"), D)),
                                strata_platform_age_base_my: c(x.strata_platform_age_base_my, UnitV1::Megayear, label(x.strata_platform_age_base_my, w.strata_platform_age_base_my, pa("project assumption: the platform cover was laid down after the -2000 My cratons and before t = 0 (A-STRATA-PLATFORM-COVER)"), D)),
                                strata_platform_age_top_my: c(x.strata_platform_age_top_my, UnitV1::Megayear, label(x.strata_platform_age_top_my, w.strata_platform_age_top_my, pa("project assumption: as strata_platform_age_base_my (A-STRATA-PLATFORM-COVER)"), D)),
                            },
                            other => other,
                        }
                    }
                }
            }
        };
        // SEED-MAGMATISM (lane MG): the magmatism entry (absent = legacy,
        // skipped on write).
        let magmatism = match orogeny
            .and_then(|o| o.rate_limited_belt.as_ref())
            .and_then(|b| b.magmatism.as_ref())
        {
            None => MagmatismLawV1::NoneV0,
            Some(m) => {
                let v1 = MagmatismBlock::v1();
                const D: &str = "declared in orogeny.rate_limited_belt.magmatism; no citation";
                let built = MagmatismLawV1::SlabArcBackarcLlsvpV1 {
                    slab_depth_below_arc_m: c(m.slab_depth_below_arc_m, UnitV1::Meter, label(m.slab_depth_below_arc_m, v1.slab_depth_below_arc_m, pa("recalled, not checked: depth to the slab under arc fronts 65-175 km, mean ~105 km (England, Engdahl and Thatcher 2004 GJI; Syracuse and Abers 2006 G3) (lane MG PREREG.txt)"), D)),
                    slab_dip_young_deg: c(m.slab_dip_young_deg, UnitV1::Degree, label(m.slab_dip_young_deg, v1.slab_dip_young_deg, pa("recalled, not checked: mean trench-to-100-km slab dip ~25-45 deg, steeper for older slabs (Uyeda and Kanamori 1979; Lallemand et al. 2005 G3); arc-trench gaps ~100-300 km (Gill 1981) (lane MG PREREG.txt)"), D)),
                    slab_dip_old_deg: c(m.slab_dip_old_deg, UnitV1::Degree, label(m.slab_dip_old_deg, v1.slab_dip_old_deg, pa("recalled, not checked: as slab_dip_young_deg (A-SLAB-DIP) (lane MG PREREG.txt)"), D)),
                    slab_dip_age_my: c(m.slab_dip_age_my, UnitV1::Megayear, label(m.slab_dip_age_my, v1.slab_dip_age_my, pa("recalled, not checked: as slab_dip_young_deg (A-SLAB-DIP) (lane MG PREREG.txt)"), D)),
                    arc_magma_productivity_m: c(m.arc_magma_productivity_m, UnitV1::Meter, label(m.arc_magma_productivity_m, v1.arc_magma_productivity_m, pa("recalled, not checked: 30 km3/km/My at 50 km/My; arc growth 20-40 km3/km/My (Reymer and Schubert 1984), intra-oceanic arcs 30-95 (Dimalanta et al. 2002), net of root foundering (Kay and Kay 1993; Jagoutz and Behn 2013) (lane MG PREREG.txt)"), D)),
                    arc_rate_floor_m_per_my: c(m.arc_rate_floor_m_per_my, UnitV1::MeterPerMegayear, label(m.arc_rate_floor_m_per_my, v1.arc_rate_floor_m_per_my, pa("project assumption: arc output depends only weakly on convergence rate (A-ARC-RATE-CLAMP) (lane MG PREREG.txt)"), D)),
                    arc_rate_cap_m_per_my: c(m.arc_rate_cap_m_per_my, UnitV1::MeterPerMegayear, label(m.arc_rate_cap_m_per_my, v1.arc_rate_cap_m_per_my, pa("project assumption (A-ARC-RATE-CLAMP) (lane MG PREREG.txt)"), D)),
                    arc_half_width_m: c(m.arc_half_width_m, UnitV1::Meter, label(m.arc_half_width_m, v1.arc_half_width_m, pa("recalled, not checked: volcanic front plus plutonic belt 50-150 km wide (Sierra Nevada batholith ~100 km) (lane MG PREREG.txt)"), D)),
                    arc_continental_thickness_m: c(m.arc_continental_thickness_m, UnitV1::Meter, label(m.arc_continental_thickness_m, v1.arc_continental_thickness_m, pa("recalled, not checked: intra-oceanic arcs above ~20 km crust carry an intermediate/felsic middle crust (Suyehiro et al. 1996; Kodaira et al. 2007) (lane MG PREREG.txt)"), D)),
                    arc_thermal_uplift_m: c(m.arc_thermal_uplift_m, UnitV1::Meter, label(m.arc_thermal_uplift_m, v1.arc_thermal_uplift_m, pa("recalled, not checked: hot thin arc/backarc lithosphere stands ~1-1.6 km above cold lithosphere with the same crust (Hyndman, Currie and Mazzotti 2005 GSA Today); taken low (lane MG PREREG.txt)"), D)),
                    arc_thermal_half_width_m: c(m.arc_thermal_half_width_m, UnitV1::Meter, label(m.arc_thermal_half_width_m, v1.arc_thermal_half_width_m, pa("project assumption (A-ARC-HEAT) (lane MG PREREG.txt)"), D)),
                    arc_thermal_decay_my: c(m.arc_thermal_decay_my, UnitV1::Megayear, label(m.arc_thermal_decay_my, v1.arc_thermal_decay_my, pa("recalled, not checked: conductive time L^2/(pi^2 kappa) of ~60 km backarc lithosphere ~12 My (lane MG PREREG.txt)"), D)),
                    backarc_slab_age_min_my: c(m.backarc_slab_age_min_my, UnitV1::Megayear, label(m.backarc_slab_age_min_my, v1.backarc_slab_age_min_my, pa("recalled, not checked: back-arc spreading sits behind old slabs (Mariana ~150, Tonga ~100, Japan ~130 My; Molnar and Atwater 1978) (lane MG PREREG.txt)"), D)),
                    backarc_overriding_advance_max_m_per_my: c(m.backarc_overriding_advance_max_m_per_my, UnitV1::MeterPerMegayear, label(m.backarc_overriding_advance_max_m_per_my, v1.backarc_overriding_advance_max_m_per_my, pa("recalled, not checked: back-arc extension where the upper plate does not advance toward the trench in the mantle frame (Heuret and Lallemand 2005 PEPI) (lane MG PREREG.txt)"), D)),
                    backarc_extension_rate_m_per_my: c(m.backarc_extension_rate_m_per_my, UnitV1::MeterPerMegayear, label(m.backarc_extension_rate_m_per_my, v1.backarc_extension_rate_m_per_my, pa("recalled, not checked: Mariana Trough 3-4 cm/yr, Lau up to 10-15, Japan Sea ~2-3 cm/yr during 25-15 Ma (lane MG PREREG.txt)"), D)),
                    backarc_axis_offset_m: c(m.backarc_axis_offset_m, UnitV1::Meter, label(m.backarc_axis_offset_m, v1.backarc_axis_offset_m, pa("recalled, not checked: Mariana Trough axis ~100-200 km behind the volcanic front (lane MG PREREG.txt)"), D)),
                    backarc_half_width_m: c(m.backarc_half_width_m, UnitV1::Meter, label(m.backarc_half_width_m, v1.backarc_half_width_m, pa("recalled, not checked: Mariana Trough ~250 km wide; Japan Sea 600-900 km after ~10 My (lane MG PREREG.txt)"), D)),
                    backarc_min_crust_m: c(m.backarc_min_crust_m, UnitV1::Meter, label(m.backarc_min_crust_m, v1.backarc_min_crust_m, pa("recalled, not checked: hyperextended continental crust of the Japan Sea Yamato Basin ~12-15 km (lane MG PREREG.txt)"), D)),
                    llsvp_edge_deg: c(m.llsvp_edge_deg, UnitV1::Degree, label(m.llsvp_edge_deg, v1.llsvp_edge_deg, pa("recalled, not checked: the two LLSVPs cover ~20-25% of the core-mantle boundary; hotspots and reconstructed LIPs lie near their margins (Burke and Torsvik 2004; Torsvik et al. 2006; Burke et al. 2008) (lane MG PREREG.txt)"), D)),
                    llsvp_edge_sigma_deg: c(m.llsvp_edge_sigma_deg, UnitV1::Degree, label(m.llsvp_edge_sigma_deg, v1.llsvp_edge_sigma_deg, pa("recalled, not checked: ~80% of reconstructed LIPs within +-10 deg of the plume generation zones (same sources) (lane MG PREREG.txt)"), D)),
                    llsvp_antipode_jitter_deg: c(m.llsvp_antipode_jitter_deg, UnitV1::Degree, label(m.llsvp_antipode_jitter_deg, v1.llsvp_antipode_jitter_deg, pa("recalled, not checked: Tuzo (Africa) and Jason (Pacific) are roughly antipodal on the equator (lane MG PREREG.txt)"), D)),
                    plume_background_fraction: c(m.plume_background_fraction, UnitV1::Dimensionless, label(m.plume_background_fraction, v1.plume_background_fraction, pa("project assumption: plumes away from the generation zones (A-PLUME-BACKGROUND) (lane MG PREREG.txt)"), D)),
                    plume_drift_rate_m_per_my: c(m.plume_drift_rate_m_per_my, UnitV1::MeterPerMegayear, label(m.plume_drift_rate_m_per_my, v1.plume_drift_rate_m_per_my, pa("recalled, not checked: modelled hotspot drift 0.5-2 cm/yr (Steinberger and O'Connell 1998); Hawaii ~4 cm/yr 80-47 Ma (Tarduno et al. 2003) is the extreme (lane MG PREREG.txt)"), D)),
                    vent_active_window_my: c(m.vent_active_window_my, UnitV1::Megayear, label(m.vent_active_window_my, v1.vent_active_window_my, pa("project assumption (A-VOLCANO-ACTIVE-WINDOW) (lane MG PREREG.txt)"), D)),
                    volcano_spacing_m: c(m.volcano_spacing_m, UnitV1::Meter, label(m.volcano_spacing_m, v1.volcano_spacing_m, pa("recalled, not checked: arc volcano spacing ~50-70 km (de Bremond d'Ars et al. 1995) (lane MG PREREG.txt)"), D)),
                    arc_crust_cap_m: c(m.arc_crust_cap_m, UnitV1::Meter, label(m.arc_crust_cap_m, v1.arc_crust_cap_m, pa("recalled, not checked: mature arc crust 30-45 km (Aleutians ~35-38, Japan ~30-35, Cascades ~40-45); thicker arc roots founder (Kay and Kay 1993; Jull and Kelemen 2001); thicker crust needs shortening (lane MG PREREG.txt revision 1)"), D)),
                };
                // SEED-WORLD-COMPOSE (lane CP): law_version 2.
                match (m.law_version, m.backarc_episode_my, m.backarc_max_stretch_factor, m.kimberlite_plume_radius_m) {
                    (1, None, None, None) => built,
                    (2 | 3, Some(e), Some(f), Some(r)) => {
                        let v2 = MagmatismBlock::v2();
                        let law = built.into_v2(
                            c(e, UnitV1::Megayear, label(e, v2.backarc_episode_my.unwrap_or(f64::NAN), pa("recalled, not checked: back-arc basins spread for ~10-30 My before stopping (Japan Sea ~25-15 Ma, Parece Vela ~30-15 Ma, Shikoku ~25-15 Ma) (lane MG idea 2; lane CP PREREG.txt) (A-BACKARC-EPISODE)"), D)),
                            c(f, UnitV1::Dimensionless, label(f, v2.backarc_max_stretch_factor.unwrap_or(f64::NAN), pa("project assumption: one back-arc episode thins continental crust by at most a factor of 2 (lane MG idea 2; lane CP PREREG.txt) (A-BACKARC-STRETCH)"), D)),
                            c(r, UnitV1::Meter, label(r, v2.kimberlite_plume_radius_m.unwrap_or(f64::NAN), pa("project assumption: kimberlite fields cluster within a few hundred km of a plume-head track (lane MG idea 4; lane CP PREREG.txt) (A-RES-KIMBERLITE-PLUME-V2)"), D)),
                        );
                        if m.law_version == 3 { law.into_v3() } else { law }
                    }
                    _ => {
                        return Err(UniverseError::Unsupported {
                            law: "magmatism/slab-arc-backarc-llsvp law_version other than 1 (no v2 constants) or 2/3 (all v2 constants)",
                            why: "only slab-arc-backarc-llsvp-v1, v2 and v3 exist",
                        });
                    }
                }
            }
        };
        // SEDIMENT-TOWERS (lane ST): the collision-sediment entry (absent =
        // legacy, skipped on write).
        let collision_sediment = match orogeny.and_then(|o| o.rate_limited_belt.as_ref()) {
            Some(b) if b.accrete_collision_sediment => {
                CollisionSedimentLawV1::AccretedToBeltCrustV1
            }
            _ => CollisionSedimentLawV1::SurfaceLayerV0,
        };
        // SEED-EXTRAS (lane EX): plume, rift and plate-split entries (absent
        // = legacy, skipped on write).
        let belt_block = orogeny.and_then(|o| o.rate_limited_belt.as_ref());
        let hotspot_volcanism = match belt_block.and_then(|b| b.plume_volcanism.as_ref()) {
            None => HotspotVolcanismLawV1::NoneV0,
            Some(h) => {
                let v1 = PlumeVolcanismBlock::v1();
                HotspotVolcanismLawV1::FixedPlumeChainsV1 {
                    plume_count: c(
                        h.plume_count as f64,
                        UnitV1::Dimensionless,
                        label(h.plume_count as f64, v1.plume_count as f64, pa("project assumption: a dozen long-lived deep plumes; recalled Courtillot et al. 2003 (EPSL 205:295) list 7 primary plumes, Montelli et al. 2004 (Science 303:338) at least 9 deep ones, French and Romanowicz 2015 (Nature 525:95) about 20 broad plumes; citations not checked (A-PLUME-COUNT)"), "declared in orogeny.rate_limited_belt.plume_volcanism.plume_count; no citation"),
                    ),
                    plume_onset_lead_my: c(
                        h.plume_onset_lead_my,
                        UnitV1::Megayear,
                        label(h.plume_onset_lead_my, v1.plume_onset_lead_my, pa("project assumption: onsets uniform in [-50 My, duration); about a third of a 100 My run's plumes are already active at t = 0 and have no head in the run (A-PLUME-ONSET)"), "declared in orogeny.rate_limited_belt.plume_volcanism.plume_onset_lead_my; no citation"),
                    ),
                    plume_head_volume_m3: c(
                        h.plume_head_volume_m3,
                        UnitV1::CubicMeter,
                        label(h.plume_head_volume_m3, v1.plume_head_volume_m3, pa("recalled: continental large igneous provinces of order 1-3 million km3 (Deccan, CAMP), oceanic plateaus larger (Ontong Java); 2e15 m3 chosen a priori; citations not checked (A-PLUME-HEAD-VOLUME)"), "declared in orogeny.rate_limited_belt.plume_volcanism.plume_head_volume_m3; no citation"),
                    ),
                    plume_head_radius_m: c(
                        h.plume_head_radius_m,
                        UnitV1::Meter,
                        label(h.plume_head_radius_m, v1.plume_head_radius_m, pa("recalled: plume heads flatten to about 1000-2000 km across beneath the lithosphere (Campbell 2005, Elements 1:265); citation not checked (A-PLUME-HEAD-RADIUS)"), "declared in orogeny.rate_limited_belt.plume_volcanism.plume_head_radius_m; no citation"),
                    ),
                    plume_tail_flux_m3_per_my: c(
                        h.plume_tail_flux_m3_per_my,
                        UnitV1::CubicMeterPerMegayear,
                        label(h.plume_tail_flux_m3_per_my, v1.plume_tail_flux_m3_per_my, pa("recalled: Hawaii's recent magma supply about 0.2 km3/yr, most plumes less; 0.1 km3/yr chosen a priori; citation not checked (A-PLUME-TAIL-FLUX)"), "declared in orogeny.rate_limited_belt.plume_volcanism.plume_tail_flux_m3_per_my; no citation"),
                    ),
                    plume_edifice_radius_m: c(
                        h.plume_edifice_radius_m,
                        UnitV1::Meter,
                        label(h.plume_edifice_radius_m, v1.plume_edifice_radius_m, pa("project assumption: a volcanic chain with its flexural moat and swell is about 200 km across (Hawaiian chain, recalled) (A-PLUME-EDIFICE-RADIUS)"), "declared in orogeny.rate_limited_belt.plume_volcanism.plume_edifice_radius_m; no citation"),
                    ),
                }
            }
        };
        let rift_block = belt_block.and_then(|b| b.continental_rift.as_ref());
        let continental_rift = match rift_block {
            None => ContinentalRiftLawV1::OceanicBirthV0,
            Some(r) => {
                let v1 = ContinentalRiftBlock::v1();
                ContinentalRiftLawV1::StretchedRiftV1 {
                    rift_half_width_m: c(
                        r.rift_half_width_m,
                        UnitV1::Meter,
                        label(r.rift_half_width_m, v1.rift_half_width_m, pa("recalled: narrow continental rifts are under ~100 km wide (Buck 1991, JGR 96:20161) and each rifted margin is 100-200 km wide (Peron-Pinvidic and Manatschal 2009); citations not checked (A-RIFT-HALF-WIDTH)"), "declared in orogeny.rate_limited_belt.continental_rift.rift_half_width_m; no citation"),
                    ),
                    breakup_crust_thickness_m: c(
                        r.breakup_crust_thickness_m,
                        UnitV1::Meter,
                        label(r.breakup_crust_thickness_m, v1.breakup_crust_thickness_m, pa("recalled: continental crust thinned below ~10 km becomes entirely brittle and the mantle is exhumed or seafloor spreading starts (Perez-Gussinye and Reston 2001, JGR 106:3961); citation not checked (A-BREAKUP-THICKNESS)"), "declared in orogeny.rate_limited_belt.continental_rift.breakup_crust_thickness_m; no citation"),
                    ),
                }
            }
        };
        let plate_reorganization_split = match rift_block.and_then(|r| r.plume_head_split.as_ref()) {
            None => PlateReorganizationLawV1::FixedPlatesV0,
            Some(sp) => PlateReorganizationLawV1::PlumeHeadContinentalSplitV1 {
                rift_opening_rate_deg_per_my: c(
                    sp.rift_opening_rate_deg_per_my,
                    UnitV1::DegreePerMegayear,
                    label(sp.rift_opening_rate_deg_per_my, 0.25, pa("recalled: young oceans open at about 1.5-4 cm/yr full rate (Red Sea, early South and Central Atlantic); 0.25 deg/My = 2.8 cm/yr at the plume on an Earth-radius planet, chosen a priori; plume-head-triggered breakup recalled from Courtillot et al. 1999 (EPSL 166:177); citations not checked (A-RIFT-OPENING-RATE)"), "declared in orogeny.rate_limited_belt.continental_rift.plume_head_split.rift_opening_rate_deg_per_my; no citation"),
                ),
            },
        };
        // SEED-PLATE-REORG (lane PR): `plate-reorg-scripted-v1` supersedes
        // the plume-head split and keeps its opening rate (absent = the
        // split law above, unchanged).
        let reorg_block = belt_block.and_then(|b| b.plate_reorg.as_ref());
        let plate_reorganization = match (reorg_block, &plate_reorganization_split) {
            (None, _) => plate_reorganization_split,
            (Some(_), PlateReorganizationLawV1::FixedPlatesV0) => {
                return Err(UniverseError::Unsupported {
                    law: "plate-reorganization/plate-reorg-scripted-v1 without continental_rift.plume_head_split",
                    why: "the scripted reorganization takes its rift opening rate from the plume-head split block and supersedes it",
                });
            }
            (Some(r), PlateReorganizationLawV1::PlumeHeadContinentalSplitV1 { rift_opening_rate_deg_per_my }) => {
                // SEED-BREAKUP-CUT-V2 (lane BC): law_version 1 or 2; labels
                // compare against that version's shipped values.
                // SEED-TRANSFORMS (lane TF): law_version 4 carries the
                // initiation age and all three initiation limits; 1-3 none.
                let v4_any = r.initiation_slab_age_min_my.is_some() || r.initiation_old_slab_fraction_min.is_some() || r.initiation_trench_spacing_m.is_some();
                let v4_all = r.initiation_slab_age_min_my.is_some() && r.initiation_old_slab_fraction_min.is_some() && r.initiation_trench_spacing_m.is_some();
                // SEED-OLD-OCEAN (lane OR): law_version 5 carries v4's
                // limits and all three margin constants; 1-4 none.
                let v5_any = r.margin_initiation_age_my.is_some() || r.margin_segment_length_m.is_some() || r.margin_trench_spacing_m.is_some();
                let v5_all = r.margin_initiation_age_my.is_some() && r.margin_segment_length_m.is_some() && r.margin_trench_spacing_m.is_some();
                if v5_any != (r.law_version == 5) || (v5_any && !v5_all) {
                    return Err(UniverseError::Unsupported {
                        law: "plate-reorganization/plate-reorg-scripted law_version 5 without all three margin-initiation constants, or 1-4 with any",
                        why: "only plate-reorg-scripted-v5 carries the margin initiation",
                    });
                }
                if !matches!((r.law_version, v4_any, v4_all), (1, false, _) | (2, false, _) | (3, false, _) | (4, true, true) | (5, true, true)) {
                    return Err(UniverseError::Unsupported {
                        law: "plate-reorganization/plate-reorg-scripted law_version 4 without all three initiation limits, or 1-3 with any",
                        why: "only plate-reorg-scripted-v4 carries the initiation limits",
                    });
                }
                if !matches!((r.law_version, r.subduction_initiation_age_my.is_some()), (1, false) | (2, false) | (3, true) | (4, true) | (5, true)) {
                    return Err(UniverseError::Unsupported {
                        law: "plate-reorganization/plate-reorg-scripted law_version other than 1, 2 or 3 (3 with its initiation age, 1-2 without)",
                        why: "only plate-reorg-scripted-v1, -v2 and -v3 exist",
                    });
                }
                let v1 = match r.law_version { 5 => PlateReorgBlock::v5(), 4 => PlateReorgBlock::v4(), 3 => PlateReorgBlock::v3(), 2 => PlateReorgBlock::v2(), _ => PlateReorgBlock::v1() };
                const D: &str = "declared in orogeny.rate_limited_belt.plate_reorg; no citation";
                let dl = UnitV1::Dimensionless;
                let my = UnitV1::Megayear;
                let built = PlateReorganizationLawV1::PlateReorgScriptedV1 {
                    rift_opening_rate_deg_per_my: rift_opening_rate_deg_per_my.clone(),
                    reorg_interval_min_my: c(r.reorg_interval_min_my, my, label(r.reorg_interval_min_my, v1.reorg_interval_min_my, pa("project assumption: pole re-solves every 50-100 My (design review review 2026-10-13); global plate reorganisations recalled near ~100 Ma and ~50 Ma (Hawaiian-Emperor bend ~47 Ma; Matthews et al. 2012, not checked) (A-REORG-INTERVAL)"), D)),
                    reorg_interval_max_my: c(r.reorg_interval_max_my, my, label(r.reorg_interval_max_my, v1.reorg_interval_max_my, pa("project assumption: see A-REORG-INTERVAL"), D)),
                    slab_pull_bias: c(r.slab_pull_bias, dl, label(r.slab_pull_bias, v1.slab_pull_bias, pa("project assumption: slab pull is the dominant plate driving force (Forsyth and Uyeda 1975; Conrad and Lithgow-Bertelloni 2002, recalled, not checked); half of a re-solved rotation axis follows the group's slab torque, half is random (A-SLAB-PULL-BIAS)"), D)),
                    slab_circumference_reference: c(r.slab_circumference_reference, dl, label(r.slab_circumference_reference, v1.slab_circumference_reference, pa("recalled: plates with more than ~30% of their circumference attached to subducting slabs move fast, 6-9 cm/yr (Forsyth and Uyeda 1975, not checked) (A-SLAB-CIRCUMFERENCE)"), D)),
                    reorg_rate_floor_fraction: c(r.reorg_rate_floor_fraction, dl, label(r.reorg_rate_floor_fraction, v1.reorg_rate_floor_fraction, pa("recalled: slab-free plates move ~1-2 cm/yr against 6-9 cm/yr for slab-attached ones (Forsyth and Uyeda 1975, not checked) (A-REORG-RATE-FLOOR)"), D)),
                    supercontinent_share: c(r.supercontinent_share, dl, label(r.supercontinent_share, v1.supercontinent_share, pa("project assumption: a suture group holding 60% of the continental crust is a supercontinent (Pangaea ~all, Gondwana ~0.6, recalled) (A-SUPERCONTINENT-SHARE)"), D)),
                    supercontinent_tenure_my: c(r.supercontinent_tenure_my, my, label(r.supercontinent_tenure_my, v1.supercontinent_tenure_my, pa("recalled: Pangaea assembled by ~300 Ma and broke up from ~200-180 Ma; mantle insulation under a supercontinent acts over ~100 My (Anderson 1982; Gurnis 1988, not checked) (A-SUPERCONTINENT-TENURE)"), D)),
                    breakup_piece_max_share: c(r.breakup_piece_max_share, dl, label(r.breakup_piece_max_share, v1.breakup_piece_max_share, pa("recalled: the largest post-Pangaea pieces (Laurasia, then Eurasia) hold ~0.35-0.4 of the continental area (A-BREAKUP-PIECE-SHARE)"), D)),
                    breakup_cascade_my: c(r.breakup_cascade_my, my, label(r.breakup_cascade_my, v1.breakup_cascade_my, pa("recalled: Central Atlantic opening ~195 Ma, then the East-West Gondwana split ~180 Ma (A-BREAKUP-CASCADE)"), D)),
                    rift_persistence_my: c(r.rift_persistence_my, my, label(r.rift_persistence_my, v1.rift_persistence_my, pa("project assumption: a young ocean or trench keeps its relative motion for 100 My across re-solves (the Atlantic has opened continuously for ~180 My, recalled) (A-RIFT-PERSISTENCE)"), D)),
                    weak_zone_strength_ratio: c(r.weak_zone_strength_ratio, dl, label(r.weak_zone_strength_ratio, v1.weak_zone_strength_ratio, pa("project assumption, factor 5 a priori: rifts reopen old sutures and orogenic belts (Wilson 1966; Tommasi and Vauchez 2001; Audet and Buergmann 2011, recalled, not checked) (A-WEAK-ZONE-STRENGTH)"), D)),
                    plate_count_min: c(r.plate_count_min as f64, dl, label(r.plate_count_min as f64, v1.plate_count_min as f64, pa("project assumption: the toy seeds 7 plates; Earth has ~7 major and ~8 medium plates (Bird 2003, recalled) (A-PLATE-COUNT)"), D)),
                    plate_count_max: c(r.plate_count_max as f64, dl, label(r.plate_count_max as f64, v1.plate_count_max as f64, pa("project assumption: see A-PLATE-COUNT"), D)),
                    plate_count_min_area_fraction: c(r.plate_count_min_area_fraction, dl, label(r.plate_count_min_area_fraction, v1.plate_count_min_area_fraction, pa("project assumption: 0.5% of the sphere, 2.55 Mkm2 on Earth (Cocos ~2.9, Caribbean ~3.3 Mkm2; Bird 2003, recalled) (A-PLATE-COUNT-AREA)"), D)),
                    subduction_initiation_rate_deg_per_my: c(r.subduction_initiation_rate_deg_per_my, UnitV1::DegreePerMegayear, label(r.subduction_initiation_rate_deg_per_my, v1.subduction_initiation_rate_deg_per_my, pa("project assumption: ~3.3 cm/yr young convergence after subduction initiation (Izu-Bonin type; Stern 2004, recalled, not checked) (A-SUBDUCTION-INITIATION-RATE)"), D)),
                    plume_head_interval_my: c(r.plume_head_interval_my, my, label(r.plume_head_interval_my, v1.plume_head_interval_my, pa("recalled: large igneous provinces about every 20-30 My over the Phanerozoic (Coffin and Eldholm 1994; Ernst 2014, not checked) (A-PLUME-HEAD-INTERVAL)"), D)),
                    plume_tail_lifetime_my: c(r.plume_tail_lifetime_my, my, label(r.plume_tail_lifetime_my, v1.plume_tail_lifetime_my, pa("recalled: hotspot tracks Hawaii-Emperor >80 My, Tristan-Walvis ~130 My, Kerguelen ~120 My (not checked) (A-PLUME-TAIL-LIFETIME)"), D)),
                    plume_split_min_share: c(r.plume_split_min_share, dl, label(r.plume_split_min_share, v1.plume_split_min_share, pa(if r.law_version >= 2 { "project assumption (lane BC, design review review 2026-10-13): a plume head under a continent holding less than 35% of the continental crust does not split it, so plumes break supercontinents and their large pieces, not small continents (A-PLUME-SPLIT-MIN-SHARE-V2)" } else { "project assumption: a plume head under a continent holding less than 20% of the continental crust does not split it (the Siberian traps did not split Siberia, recalled) (A-PLUME-SPLIT-MIN-SHARE)" }), D)),
                };
                match r.law_version {
                    2 => built.into_scripted_v2(),
                    // SEED-WORLD-COMPOSE (lane CP): v3.
                    3 => {
                        let a = r.subduction_initiation_age_my.expect("checked above");
                        built.into_scripted_v2().into_scripted_v3(c(a, my, label(a, v1.subduction_initiation_age_my.unwrap_or(f64::NAN), pa("recalled, not checked: the oldest in-situ ocean floor is ~180 My (Pigafetta Basin) and ocean-age p95 ~150-180 My (Mueller et al. 2008); passive margins and enclosed basins with older floor are where subduction must start (lane VL defect 3; lane CP PREREG.txt) (A-SUBDUCTION-INITIATION-AGE)"), D)))
                    }
                    // SEED-TRANSFORMS (lane TF): v4 = v3 plus the
                    // age-triggered initiation limits.
                    4 => {
                        let a = r.subduction_initiation_age_my.expect("checked above");
                        let sa = r.initiation_slab_age_min_my.expect("checked above");
                        let of = r.initiation_old_slab_fraction_min.expect("checked above");
                        let sp = r.initiation_trench_spacing_m.expect("checked above");
                        built
                            .into_scripted_v2()
                            .into_scripted_v3(c(a, my, label(a, v1.subduction_initiation_age_my.unwrap_or(f64::NAN), pa("recalled, not checked: the oldest in-situ ocean floor is ~180 My (Pigafetta Basin) and ocean-age p95 ~150-180 My (Mueller et al. 2008); passive margins and enclosed basins with older floor are where subduction must start (lane VL defect 3; lane CP PREREG.txt) (A-SUBDUCTION-INITIATION-AGE)"), D)))
                            .into_scripted_v4(
                                c(sa, my, label(sa, v1.initiation_slab_age_min_my.unwrap_or(f64::NAN), pa("project assumption (lane TF PREREG.txt): spontaneous subduction initiation needs old, dense lithosphere; the v3 weakening starts at 100 My; Izu-Bonin-Mariana started ~52 Ma on Cretaceous Pacific floor (Stern 2004; Arculus et al. 2015, recalled, not checked) (A-INITIATION-SLAB-AGE)"), D)),
                                c(of, dl, label(of, v1.initiation_old_slab_fraction_min.unwrap_or(f64::NAN), pa("project assumption (lane TF PREREG.txt): at least half of a new age-triggered trench runs in old slab, so cuts do not wander through young ocean (lane CP: v3 cuts were 3,000-46,000 km long with 0-60% weak length) (A-INITIATION-SLAB-AGE)"), D)),
                                c(sp, UnitV1::Meter, label(sp, v1.initiation_trench_spacing_m.unwrap_or(f64::NAN), pa("project assumption (lane TF PREREG.txt): arc systems on Earth sit ~1000 km or more apart (Izu-Bonin to Ryukyu ~1000 km, Tonga to Vanuatu ~1500 km, recalled, not checked); at most half of a new trench may run within this of existing subduction (A-INITIATION-SPACING)"), D)),
                            )
                    }
                    // SEED-OLD-OCEAN (lane OR): v5 = v4 plus the
                    // passive-margin initiation.
                    5 => {
                        let a = r.subduction_initiation_age_my.expect("checked above");
                        let sa = r.initiation_slab_age_min_my.expect("checked above");
                        let of = r.initiation_old_slab_fraction_min.expect("checked above");
                        let sp = r.initiation_trench_spacing_m.expect("checked above");
                        let ma = r.margin_initiation_age_my.expect("checked above");
                        let ml = r.margin_segment_length_m.expect("checked above");
                        let ms = r.margin_trench_spacing_m.expect("checked above");
                        built
                            .into_scripted_v2()
                            .into_scripted_v3(c(a, my, label(a, v1.subduction_initiation_age_my.unwrap_or(f64::NAN), pa("recalled, not checked: the oldest in-situ ocean floor is ~180 My (Pigafetta Basin) and ocean-age p95 ~150-180 My (Mueller et al. 2008); passive margins and enclosed basins with older floor are where subduction must start (lane VL defect 3; lane CP PREREG.txt) (A-SUBDUCTION-INITIATION-AGE)"), D)))
                            .into_scripted_v4(
                                c(sa, my, label(sa, v1.initiation_slab_age_min_my.unwrap_or(f64::NAN), pa("project assumption (lane TF PREREG.txt): spontaneous subduction initiation needs old, dense lithosphere; the v3 weakening starts at 100 My; Izu-Bonin-Mariana started ~52 Ma on Cretaceous Pacific floor (Stern 2004; Arculus et al. 2015, recalled, not checked) (A-INITIATION-SLAB-AGE)"), D)),
                                c(of, dl, label(of, v1.initiation_old_slab_fraction_min.unwrap_or(f64::NAN), pa("project assumption (lane TF PREREG.txt): at least half of a new age-triggered trench runs in old slab, so cuts do not wander through young ocean (lane CP: v3 cuts were 3,000-46,000 km long with 0-60% weak length) (A-INITIATION-SLAB-AGE)"), D)),
                                c(sp, UnitV1::Meter, label(sp, v1.initiation_trench_spacing_m.unwrap_or(f64::NAN), pa("project assumption (lane TF PREREG.txt): arc systems on Earth sit ~1000 km or more apart (Izu-Bonin to Ryukyu ~1000 km, Tonga to Vanuatu ~1500 km, recalled, not checked); at most half of a new trench may run within this of existing subduction (A-INITIATION-SPACING)"), D)),
                            )
                            .into_scripted_v5(
                                c(ma, my, label(ma, v1.margin_initiation_age_my.unwrap_or(f64::NAN), pa("recalled, not checked: ocean floor older than ~160-180 My is almost absent on Earth (Mueller et al. 2008); Atlantic-type passive margins are expected to fail and start subducting as they age (Stern 2004; Nikolaeva et al. 2010) (lane OR PREREG.txt) (A-MARGIN-INITIATION-AGE)"), D)),
                                c(ml, UnitV1::Meter, label(ml, v1.margin_segment_length_m.unwrap_or(f64::NAN), pa("project assumption (lane OR PREREG.txt): a new margin trench starts along at most ~2000 km of coast (Lesser Antilles ~850 km, Scotia ~700 km, Izu-Bonin-Mariana ~2800 km, recalled, not checked), so a single initiation does not outline a whole continent (A-MARGIN-SEGMENT-LENGTH)"), D)),
                                c(ms, UnitV1::Meter, label(ms, v1.margin_trench_spacing_m.unwrap_or(f64::NAN), pa("project assumption (lane OR PREREG.txt): as A-INITIATION-SPACING, arc systems ~1000 km or more apart (A-MARGIN-TRENCH-SPACING)"), D)),
                            )
                    }
                    _ => built,
                }
            }
            (Some(_), _) => {
                return Err(UniverseError::Unsupported {
                    law: "plate-reorganization/plate-reorg-scripted-v1",
                    why: "unexpected plate reorganization law under the scripted reorganization",
                });
            }
        };
        // SEED-RIFT-SHOULDERS (lane RS): the rift flank entry (absent =
        // legacy, skipped on write).
        let rift_flank = match rift_block.and_then(|r| r.rift_flank.as_ref()) {
            None => RiftFlankLawV1::NoneV0,
            Some(f) => {
                let v1 = RiftFlankBlock::v1();
                RiftFlankLawV1::ShoulderUpliftThermalV1 {
                    rift_flank_flexural_parameter_m: c(
                        f.rift_flank_flexural_parameter_m,
                        UnitV1::Meter,
                        label(f.rift_flank_flexural_parameter_m, v1.rift_flank_flexural_parameter_m, pa("recalled: rift flanks are flexed by a broken plate of effective elastic thickness about 10-30 km (Weissel and Karner 1989, JGR 94:13919; Ebinger et al. 1991 for East Africa); alpha = (4D / (rho_m g))^(1/4) is about 50-70 km for Te 20-25 km, so shoulders are about 50-100 km wide; 60 km chosen a priori; citations not checked (A-RIFT-FLANK-ALPHA)"), "declared in orogeny.rate_limited_belt.continental_rift.rift_flank.rift_flank_flexural_parameter_m; no citation"),
                    ),
                    rift_shoulder_ratio: c(
                        f.rift_shoulder_ratio,
                        UnitV1::Dimensionless,
                        label(f.rift_shoulder_ratio, v1.rift_shoulder_ratio, pa("project assumption: the footwall at the border fault rises with the buoyant rift floor next to it (coupled at the fault, ratio 1), so a young rift has shoulders of order the floor's thermal uplift (recalled: young rift shoulders 0.5-2 km, Red Sea, Gulf of Suez, Rhine graben); chosen a priori; citations not checked (A-RIFT-SHOULDER-RATIO)"), "declared in orogeny.rate_limited_belt.continental_rift.rift_flank.rift_shoulder_ratio; no citation"),
                    ),
                }
            }
        };
        // SEED-FLOODBASALT (lane FB): the LIP and swell entry (absent =
        // legacy, skipped on write).
        let flood_basalt = match belt_block
            .and_then(|b| b.plume_volcanism.as_ref())
            .and_then(|p| p.flood_basalt.as_ref())
        {
            None => FloodBasaltLawV1::NoneV0,
            Some(f) => {
                let v1 = FloodBasaltBlock::v1();
                const D: &str = "declared in orogeny.rate_limited_belt.plume_volcanism.flood_basalt; no citation";
                FloodBasaltLawV1::PlumeHeadLipSwellV1 {
                    flood_basalt_extrusive_volume_m3: c(f.flood_basalt_extrusive_volume_m3, UnitV1::CubicMeter, label(f.flood_basalt_extrusive_volume_m3, v1.flood_basalt_extrusive_volume_m3, pa("recalled extrusive volumes: Columbia River ~2.1e5 km3 (Reidel et al. 2013, GSA SP 497), Ethiopia-Yemen ~3.5e5 km3, Deccan ~1e6+, Parana-Etendeka ~1.5-2e6, CAMP ~2-3e6, Siberian ~2-4e6 km3 (Courtillot and Renne 2003, C. R. Geosci. 335:113; Reichow et al. 2009, EPSL 277:9); one value for all heads, 1.5e6 km3 chosen a priori; citations not checked (A-LIP-VOLUME)"), D)),
                    flood_basalt_intrusive_ratio: c(f.flood_basalt_intrusive_ratio, UnitV1::Dimensionless, label(f.flood_basalt_intrusive_ratio, v1.flood_basalt_intrusive_ratio, pa("recalled: intrusive:extrusive ~5:1 to 10:1 (Crisp 1984, JVGR 20:177); 5-15 km of high-velocity underplate beneath LIPs and volcanic margins (White and McKenzie 1989, JGR 94:7685); 4 chosen at the low end because all continental crust loads at 2800 kg/m3 here; citations not checked (A-LIP-VOLUME)"), D)),
                    flood_basalt_radius_m: c(f.flood_basalt_radius_m, UnitV1::Meter, label(f.flood_basalt_radius_m, v1.flood_basalt_radius_m, pa("recalled: plume heads flatten to 1000-2000 km across (Campbell 2005, Elements 1:265); Deccan's original extent ~1.5e6 km2; citations not checked (A-LIP-RADIUS)"), D)),
                    flood_basalt_taper_fraction: c(f.flood_basalt_taper_fraction, UnitV1::Dimensionless, label(f.flood_basalt_taper_fraction, v1.flood_basalt_taper_fraction, pa("project assumption: flat lava-pile plateau thinning over the outer 30 % of the radius (A-LIP-SHAPE)"), D)),
                    flood_basalt_eruption_duration_my: c(f.flood_basalt_eruption_duration_my, UnitV1::Megayear, label(f.flood_basalt_eruption_duration_my, v1.flood_basalt_eruption_duration_my, pa("recalled: main phases ~1 My (Deccan, Siberian, CAMP: Schoene et al. 2019, Burgess and Bowring 2015, Blackburn et al. 2013), Ethiopia 1-2 My, Ontong Java a few My; 2 My chosen; citations not checked (A-LIP-DURATION)"), D)),
                    swell_tail_amplitude_m: c(f.swell_tail_amplitude_m, UnitV1::Meter, label(f.swell_tail_amplitude_m, v1.swell_tail_amplitude_m, pa("recalled: hotspot swells 0.5-1.5 km (Crough 1983, Annu. Rev. Earth Planet. Sci. 11:165), Hawaiian swell ~1-1.5 km; 800 m water-loaded chosen; citation not checked (A-PLUME-SWELL)"), D)),
                    swell_tail_radius_m: c(f.swell_tail_radius_m, UnitV1::Meter, label(f.swell_tail_radius_m, v1.swell_tail_radius_m, pa("recalled: Hawaiian swell ~1000-1200 km across (Sleep 1990, JGR 95:6715); Gaussian e-fold 500 km; citation not checked (A-PLUME-SWELL)"), D)),
                    swell_head_amplitude_m: c(f.swell_head_amplitude_m, UnitV1::Meter, label(f.swell_head_amplitude_m, v1.swell_head_amplitude_m, pa("recalled: ~0.5-1 km of domal uplift over ~2000 km before and during flood volcanism (Campbell and Griffiths 1990, EPSL 99:79; Campbell 2005); East African dynamic support ~1-1.5 km; 700 m water-loaded above the tail chosen; citations not checked (A-PLUME-SWELL)"), D)),
                    swell_head_radius_m: c(f.swell_head_radius_m, UnitV1::Meter, label(f.swell_head_radius_m, v1.swell_head_radius_m, pa("project assumption: head swell e-fold = the province radius, 800 km (A-PLUME-SWELL)"), D)),
                    swell_head_decay_my: c(f.swell_head_decay_my, UnitV1::Megayear, label(f.swell_head_decay_my, v1.swell_head_decay_my, pa("project assumption: the flattened head spreads and cools over ~10-20 My; e-fold 10 My (A-PLUME-SWELL)"), D)),
                }
            }
        };
        // SEED-TRENCH part 2 (lane TN): the seeded passive-margin entry.
        let passive_margin = match rift_block.and_then(|r| r.rift_flank.as_ref()) {
            Some(f) if f.seeded_margin_subsidence => PassiveMarginLawV1::SeededMarginThermalSubsidenceV1,
            _ => PassiveMarginLawV1::NoneV0,
        };
        // SEED-TRENCH (lane TN): the subduction trench entry (absent =
        // legacy, skipped on write).
        let subduction_trench = match belt_block.and_then(|b| b.subduction_trench.as_ref()) {
            None => SubductionTrenchLawV1::NoneV0,
            Some(t) => {
                let v1 = SubductionTrenchBlock::v1();
                const D: &str = "declared in orogeny.rate_limited_belt.subduction_trench; no citation";
                SubductionTrenchLawV1::SlabFlexureTrenchV1 {
                    trench_vertical_load_fraction: c(t.trench_vertical_load_fraction, UnitV1::Dimensionless, label(t.trench_vertical_load_fraction, v1.trench_vertical_load_fraction, pa("project assumption: share of the slab pull (recalled ~2-5e13 N/m) carried as the vertical end load of the bending plate; trench-profile fits give end loads ~1e12-5e12 N/m (recalled: Caldwell et al. 1976, EPSL 31:239; Turcotte and Schubert 2014 ch. 3); 0.1 chosen a priori; citations not checked (A-TRENCH-LOAD)"), D)),
                    trench_slab_length_max_m: c(t.trench_slab_length_max_m, UnitV1::Meter, label(t.trench_slab_length_max_m, v1.trench_slab_length_max_m, pa("recalled: the upper-mantle slab down to the 660 km discontinuity carries the slab pull; citation not checked (A-TRENCH-LOAD)"), D)),
                    trench_thermal_expansivity_per_k: c(t.trench_thermal_expansivity_per_k, UnitV1::PerKelvin, label(t.trench_thermal_expansivity_per_k, v1.trench_thermal_expansivity_per_k, pa("recalled: Turcotte and Schubert 2014, Geodynamics, mantle rock 3e-5 /K; citation not checked (A-TRENCH-LOAD)"), D)),
                    trench_mantle_temperature_c: c(t.trench_mantle_temperature_c, UnitV1::DegreeCelsius, label(t.trench_mantle_temperature_c, v1.trench_mantle_temperature_c, pa("recalled: half-space cooling fits with T_m - T_0 ~1300-1350 C (Parsons and Sclater 1977; Turcotte and Schubert 2014); citation not checked (A-TRENCH-PLATE)"), D)),
                    trench_elastic_isotherm_c: c(t.trench_elastic_isotherm_c, UnitV1::DegreeCelsius, label(t.trench_elastic_isotherm_c, v1.trench_elastic_isotherm_c, pa("recalled: oceanic elastic thickness follows the 450 +- 150 C isotherm (Watts 2001, Isostasy and Flexure of the Lithosphere); citation not checked (A-TRENCH-PLATE)"), D)),
                    trench_thermal_diffusivity_m2_per_yr: c(t.trench_thermal_diffusivity_m2_per_yr, UnitV1::SquareMeterPerYear, label(t.trench_thermal_diffusivity_m2_per_yr, v1.trench_thermal_diffusivity_m2_per_yr, pa("recalled: kappa = 1e-6 m2/s (Turcotte and Schubert 2014); citation not checked (A-TRENCH-PLATE)"), D)),
                    trench_elastic_modulus_pa: c(t.trench_elastic_modulus_pa, UnitV1::Pascal, label(t.trench_elastic_modulus_pa, v1.trench_elastic_modulus_pa, pa("recalled: E = 70 GPa in the Turcotte and Schubert 2014 lithosphere flexure examples; citation not checked (A-TRENCH-PLATE)"), D)),
                    trench_poisson_ratio: c(t.trench_poisson_ratio, UnitV1::Dimensionless, label(t.trench_poisson_ratio, v1.trench_poisson_ratio, pa("recalled: nu = 0.25 (Turcotte and Schubert 2014); citation not checked (A-TRENCH-PLATE)"), D)),
                }
            }
        };
        // SEED-UPLANDS (lane UP): the threshold hillslope entry (absent =
        // legacy, skipped on write).
        let hillslope_threshold = match orogeny.and_then(|o| o.critical_slope_erosion.as_ref()) {
            None => HillslopeThresholdLawV1::NoneV0,
            Some(b) => {
                let v1 = CriticalSlopeErosionBlock::v1();
                HillslopeThresholdLawV1::CriticalSlopeEnvelopeV1 {
                    hillslope_critical_slope: c(
                        b.hillslope_critical_slope,
                        UnitV1::Dimensionless,
                        label(b.hillslope_critical_slope, v1.hillslope_critical_slope, pa("critical slope of 0.625 deg (~70 km) cell means. Recalled range-front geometry: Himalayan front ~4.5-5 km over 120-150 km, central Andes western flank ~4 km over 100-150 km; threshold hillslopes ~0.6-0.7 at the hillslope scale (Burbank et al. 1996, Nature 379:505; Montgomery and Brandon 2002, EPSL 201:481; Roering et al. 1999, WRR 35:853); citations not checked. Measured check: p99.5 of Earth's steepest cell-mean gradients at this grid (ETOPO 2022, land below 60 deg) is 0.0375; 0.04 declared before any run (A-HILLSLOPE-CRITICAL-SLOPE)"), "declared in orogeny.critical_slope_erosion.hillslope_critical_slope; no citation"),
                    ),
                    knickpoint_celerity_m_per_my: c(
                        b.knickpoint_celerity_m_per_my,
                        UnitV1::MeterPerMegayear,
                        label(b.knickpoint_celerity_m_per_my, v1.knickpoint_celerity_m_per_my, pa("rate at which an over-steepened front is cut back, tau = cell size / celerity (1 My at 70 km). Recalled bedrock knickpoint retreat ~1 to >100 km/My depending on drainage area (Crosby and Whipple 2006, Geomorphology 82:16; Whittaker and Boulton 2012, JGR 117) and active threshold-range erosion of 1-10 mm/yr (Burbank et al. 1996; Dadson et al. 2003); citations not checked; 70 km/My declared before any run (A-KNICKPOINT-CELERITY)"), "declared in orogeny.critical_slope_erosion.knickpoint_celerity_m_per_my; no citation"),
                    ),
                }
            }
        };
        // INCISION-THRESHOLD (lane IT): absent = legacy, skipped on write.
        // LAKE-HOLD-2: `incision_threshold_trunk_only` selects the trunk-only
        // variant (same constant) and needs the threshold block.
        let trunk_only = orogeny.is_some_and(|o| o.incision_threshold_trunk_only);
        if trunk_only && orogeny.and_then(|o| o.incision_threshold.as_ref()).is_none() {
            return Err(UniverseError::Unsupported {
                law: "fluvial-incision-threshold/stream-power-excess-steepness-threshold-trunk-only-v1 without an incision_threshold block",
                why: "orogeny.incision_threshold_trunk_only restricts the incision threshold to the trunk; set orogeny.incision_threshold too",
            });
        }
        let fluvial_incision_threshold = match orogeny.and_then(|o| o.incision_threshold.as_ref()) {
            None => FluvialIncisionThresholdLawV1::NoneV0,
            Some(b) => {
                let critical_channel_steepness_m = c(
                    b.critical_channel_steepness_m,
                    UnitV1::Meter,
                    label(b.critical_channel_steepness_m, IncisionThresholdBlock::v1().critical_channel_steepness_m, pa("project assumption awaiting review: excess-stream-power incision threshold k_c (theta = 0.5 units). San Gabriel channels with k_s 30-180 m^0.9 (theta 0.45, ~70-100 m and up here) still incise at 35-1100 m/Ma in a threshold-dominated regime with tau_c 45 Pa (DiBiase and Whipple 2011, JGR 116 F04036, text read); k_c is kept below that; rate equivalent K_n k_c^2 = 3.6 m/My at K 2e-6/yr, k_ref 500 m (cratons 1-10 m/My, recalled). Threshold extends post-orogenic decay ~x20 (Baldwin, Whipple and Tucker 2003, JGR 108(B3), abstract read by lane EARTH-NUMBERS). 30 m declared before any run (lane IT PREREG.txt) (A-INCISION-THRESHOLD)"), "declared in orogeny.incision_threshold.critical_channel_steepness_m; no citation"),
                );
                if trunk_only {
                    FluvialIncisionThresholdLawV1::ExcessSteepnessThresholdTrunkOnlyV1 { critical_channel_steepness_m }
                } else {
                    FluvialIncisionThresholdLawV1::ExcessSteepnessThresholdV1 { critical_channel_steepness_m }
                }
            }
        };
        // ALLUVIAL-TRUNK: `incision_threshold_alluvial_channels` selects the
        // alluvial-channel release; it reads k_c and the retention law's S_a.
        let fluvial_alluvial_channel = if orogeny.is_some_and(|o| o.incision_threshold_alluvial_channels) {
            if orogeny.and_then(|o| o.incision_threshold.as_ref()).is_none() || orogeny.and_then(|o| o.sediment_retention.as_ref()).is_none() {
                return Err(UniverseError::Unsupported {
                    law: "fluvial-alluvial-channel/stream-power-threshold-free-alluvial-channels-v1 without an incision_threshold or sediment_retention block",
                    why: "orogeny.incision_threshold_alluvial_channels releases the incision threshold (k_c) where k_c / sqrt(A) is below the retention law's alluvial grade S_a; set orogeny.incision_threshold and orogeny.sediment_retention too",
                });
            }
            crate::laws::FluvialAlluvialChannelLawV1::ThresholdFreeAlluvialChannelsV1
        } else {
            crate::laws::FluvialAlluvialChannelLawV1::NoneV0
        };
        // GLACIAL-T0 (lane GLACIAL-T0): absent = legacy, skipped on write.
        let glacial_erosion = match orogeny.and_then(|o| o.glacial_erosion.as_ref()) {
            None => GlacialErosionLawV1::NoneV0,
            Some(b) => {
                let v1 = GlacialErosionBlock::v1();
                let v1w = |i: usize| v1.warm_states.get(i).cloned();
                GlacialErosionLawV1::GlacialStreamPowerBlendV1 {
                    glacial_to_fluvial_erodibility_ratio: c(
                        b.glacial_to_fluvial_erodibility_ratio,
                        UnitV1::Dimensionless,
                        label(b.glacial_to_fluvial_erodibility_ratio, v1.glacial_to_fluvial_erodibility_ratio, pa("project assumption awaiting review: K_g/K_f = 19/2.5 = 7.6 from Liebl et al. 2023 (GMD 16:1315, text read): OpenLEM glacial erodibility K_g = 19 /Ma calibrated to iSOSIA (same eroded volume over 100 kyr, m/n = 0.5, n = 1, reference precipitation 1 m/a) and their fluvial steady state k_sn = 250 m at U = 0.625 mm/a (K_f = U/k_sn = 2.5 /Ma, the Robl et al. 2017 value Hergarten 2021 quotes); the n = 2 calibration (K_g 24 m^-1 Ma^-1, same volume) was not converted; Hergarten 2021 (ESurf 9:937) used 2 only for display (A-GLACIAL-ERODIBILITY-RATIO)"), "declared in orogeny.glacial_erosion.glacial_to_fluvial_erodibility_ratio; no citation"),
                    ),
                    full_ice_height_above_ela_m: c(
                        b.full_ice_height_above_ela_m,
                        UnitV1::Meter,
                        label(b.full_ice_height_above_ela_m, v1.full_ice_height_above_ela_m, pa("project assumption awaiting review: Liebl et al. 2023 (GMD 16:1315, text read) iSOSIA PDD at 1 m/a precipitation gives ELA 1330 m and full-ice altitude 1830 m (500 m), the calibration of OpenLEM's linear ramp (their eq 10); Hergarten 2021 (ESurf 9:937) example He 2000 m, Hf 2500 m; one value for all climates (mass-balance gradient not precipitation-dependent) (A-GLACIAL-FULL-ICE-HEIGHT)"), "declared in orogeny.glacial_erosion.full_ice_height_above_ela_m; no citation"),
                    ),
                    coldhouse_duty: c(
                        b.coldhouse_duty,
                        UnitV1::Dimensionless,
                        label(b.coldhouse_duty, v1.coldhouse_duty, pa("project assumption awaiting review: coldhouse (GMST 11-18 C, like today) 13 % of Phanerozoic time, Judd et al. 2024 (Science 385 eadk3705, PhanDA, text read); glacial cycles run in this state only (A-GLACIAL-COLDHOUSE-DUTY)"), "declared in orogeny.glacial_erosion.coldhouse_duty; no citation"),
                    ),
                    glacial_max_cooling_k: c(
                        b.glacial_max_cooling_k,
                        UnitV1::Kelvin,
                        label(b.glacial_max_cooling_k, v1.glacial_max_cooling_k, pa("project assumption awaiting review: glacial-maximum mountain cooling 7 K with f_P 0.7 per lane GLACIAL-CHECK (unpublished-design-record, recommendation 1a): global LGM -6.1 C (Tierney et al. 2020), NZ Southern Alps >= 6-6.5 C, observed ELA depressions ~900 m (Uinta) to 1000-1250 m (central Andes); 154 m per K at lapse 6.5 K/km; cycle sinusoidal in time (Liebl et al. 2023 sec 3.5) (A-GLACIAL-MAX-COOLING)"), "declared in orogeny.glacial_erosion.glacial_max_cooling_k; no citation"),
                    ),
                    glacial_max_precipitation_factor: c(
                        b.glacial_max_precipitation_factor,
                        UnitV1::Dimensionless,
                        label(b.glacial_max_precipitation_factor, v1.glacial_max_precipitation_factor, pa("project assumption awaiting review: LGM precipitation 0.7 of modern, lane GLACIAL-CHECK (NZ 0.75-1.0, Seguinot 2018 Alps ~0.66 at -6 K; central Chile 33 S wetter: not latitude-aware here) (A-GLACIAL-MAX-PRECIPITATION)"), "declared in orogeny.glacial_erosion.glacial_max_precipitation_factor; no citation"),
                    ),
                    warm_state_duties: b
                        .warm_states
                        .iter()
                        .enumerate()
                        .map(|(i, w)| {
                            c(
                                w.duty,
                                UnitV1::Dimensionless,
                                match v1w(i) {
                                    Some(d) if d.duty.to_bits() == w.duty.to_bits() => pa("project assumption awaiting review: Judd et al. 2024 (Science 385 eadk3705, text read) time shares: coolhouse 18 % (icehouse 31 % minus coldhouse 13 %), transitional 28 % (100 - 41 - 31), greenhouse (warmhouse + hothouse) 41 % (A-GLACIAL-WARM-DUTY)"),
                                    _ => pa("declared in orogeny.glacial_erosion.warm_states[].duty; no citation"),
                                },
                            )
                        })
                        .collect(),
                    warm_state_warmings_k: b
                        .warm_states
                        .iter()
                        .enumerate()
                        .map(|(i, w)| {
                            c(
                                w.warming_k,
                                UnitV1::Kelvin,
                                match v1w(i) {
                                    Some(d) if d.warming_k.to_bits() == w.warming_k.to_bits() => pa("project assumption awaiting review: GMST midpoints of the Judd et al. 2024 states (coolhouse 18-22 C -> 20, transitional 22-25 -> 23.5, greenhouse 25-36 -> 28) minus pre-industrial ~14 C (recalled); applied as a uniform offset (no polar amplification, which Judd et al. report) (A-GLACIAL-WARM-OFFSET)"),
                                    _ => pa("declared in orogeny.glacial_erosion.warm_states[].warming_k; no citation"),
                                },
                            )
                        })
                        .collect(),
                }
            }
        };
        let marine_remainder_layer = if orogeny.map_or(false, |o| o.once_per_step_marine_layer) {
            MarineRemainderLayerLawV1::OncePerStepCappedGlobalLayerV1
        } else {
            MarineRemainderLayerLawV1::PerMouthGlobalLayerV0
        };
        // SEED-SEDIMENT-RETAIN (lane SR): the retention entry (absent =
        // legacy, skipped on write).
        let sediment_retention = match orogeny.and_then(|o| o.sediment_retention.as_ref()) {
            None => SedimentRetentionLawV1::NoneV0,
            Some(b) => {
                let v1 = SedimentRetentionBlock::v1();
                const D: &str = "declared in orogeny.sediment_retention; no citation";
                SedimentRetentionLawV1::ContinentalRetentionBaseLevelV1 {
                    base_level_relief_scale_m: c(b.base_level_relief_scale_m, UnitV1::Meter, label(b.base_level_relief_scale_m, v1.base_level_relief_scale_m, pa("project assumption: local relief (max - min level over the cell and its 8 neighbours) below which stream-power lowering is tapered, phi = min(1, (r/H)^p) (rule 2; denudation ~ local relief, Ahnert 1970, Am J Sci 268:243, recalled). Recalled: low-relief cratons/shields denude ~1-10 m/My (von Blanckenburg 2005, EPSL 237:462; Bierman and Caffee 2002, Am J Sci 301:326; Portenga and Bierman 2011, GSA Today 21(8)), old uplands ~20-30 m/My (Matmon et al. 2003, Geology 31:155), erosion rising nonlinearly with relief (Montgomery and Brandon 2002, EPSL 201:481); citations not checked; 1000 m declared before any run (A-BASE-LEVEL-RELIEF-SCALE)"), D)),
                    base_level_taper_exponent: c(b.base_level_taper_exponent, UnitV1::Dimensionless, label(b.base_level_taper_exponent, v1.base_level_taper_exponent, pa("project assumption: taper exponent p; with the baseline's ~linear lowland lowering it gives ~2-5 m/My at 300 m, ~10-25 m/My at 500 m, no change above H (same recalled sources as A-BASE-LEVEL-RELIEF-SCALE); 2 declared before any run (A-BASE-LEVEL-TAPER-EXPONENT)"), D)),
                    alluvial_grade_slope: c(b.alluvial_grade_slope, UnitV1::Dimensionless, label(b.alluvial_grade_slope, v1.alluvial_grade_slope, pa("project assumption: graded slope of lowland alluvial cell means; recalled lowland river gradients ~2e-5 (Amazon) to ~1-2e-4 (lower Mississippi, Ganga plain), fans and piedmonts steeper but sub-cell at 70 km; not checked; 1e-4 declared before any run (A-ALLUVIAL-GRADE-SLOPE)"), D)),
                    shelf_wave_base_depth_m: c(b.shelf_wave_base_depth_m, UnitV1::Meter, label(b.shelf_wave_base_depth_m, v1.shelf_wave_base_depth_m, pa("lane SD's shelf wave base, unchanged: storm wave base / cell-mean shelf-top depth, Earth's shelf break ~130 m recalled, not checked (A-CONTINENTAL-SHELF-WAVE-BASE)"), D)),
                }
            }
        };
        if orogeny.map_or(false, |o| {
            o.subgrid_erosion.is_some() && o.nonlinear_subgrid_erosion.is_some()
        }) {
            return Err(UniverseError::Inconsistent(
                "orogeny.subgrid_erosion and orogeny.nonlinear_subgrid_erosion both select the erosion law; choose one",
            ));
        }
        let erosion = match (
            orogeny.and_then(|o| o.subgrid_erosion.as_ref()),
            orogeny.and_then(|o| o.nonlinear_subgrid_erosion.as_ref()),
        ) {
            (None, Some(b)) => {
                let v1 = NonlinearSubgridErosionBlock::v1();
                ErosionLawV1::StreamPowerSubgridNonlinearV1 {
                    channel_head_length_m: c(
                        b.channel_head_length_m,
                        UnitV1::Meter,
                        label(b.channel_head_length_m, v1.channel_head_length_m, pa("lane OG's subgrid channel-head length, unchanged: recalled hillslope lengths / channel-head source distances of order 100 m (Montgomery and Dietrich 1988, Nature 336:232; Roering et al. 2007, Geology 35:571); citations not checked; enters only through ln(dx/x_c) (A-SUBGRID-CHANNEL-HEAD)"), "declared in orogeny.nonlinear_subgrid_erosion.channel_head_length_m; no citation"),
                    ),
                    slope_exponent_n: c(
                        b.slope_exponent_n,
                        UnitV1::Dimensionless,
                        label(b.slope_exponent_n, v1.slope_exponent_n, pa("recalled: mean-flow erosion rises faster than linearly with channel steepness where floods are stochastic and incision has a threshold, effective n ~ 2 or more (Lague 2014, ESPL 39:38; DiBiase and Whipple 2011, JGR 116 F04036; Ouimet et al. 2009, Geology 37:579); a global compilation favours n > 1 (Harel, Mudd and Attal 2016, Geomorphology 268:184); n > 1 lets post-orogenic relief persist (Baldwin, Whipple and Tucker 2003, JGR 108(B3)); citations not checked; 2 declared before any run (A-STREAM-POWER-N-NONLINEAR)"), "declared in orogeny.nonlinear_subgrid_erosion.slope_exponent_n; no citation"),
                    ),
                    reference_steepness_m: c(
                        b.reference_steepness_m,
                        UnitV1::Meter,
                        label(b.reference_steepness_m, v1.reference_steepness_m, pa("channel steepness (theta = 0.5 units) at which the n-law equals the n = 1 law, read as where the declared K holds (active ranges). Recalled anchors: active-range k_sn ~ 100-500 m^0.9 at theta_ref 0.45 (Kirby and Whipple 2012, J Struct Geol 44:54; Wobus et al. 2006, GSA SP 398) x A^0.05 ~ 2-3 = ~250-1500 m; K k_ref = 2e-6/yr x 500 m = 1 mm/yr, an active-orogen erosion rate (Montgomery and Brandon 2002, EPSL 201:481); citations not checked; 500 m declared before any run (A-REFERENCE-STEEPNESS)"), "declared in orogeny.nonlinear_subgrid_erosion.reference_steepness_m; no citation"),
                    ),
                }
            }
            (None, None) => ErosionLawV1::StreamPowerN1LinearDiffusionV0,
            (Some(b), _) => ErosionLawV1::StreamPowerSubgridV1 {
                channel_head_length_m: c(
                    b.channel_head_length_m,
                    UnitV1::Meter,
                    label(b.channel_head_length_m, SubgridErosionBlock::v1().channel_head_length_m, pa("recalled: hillslope lengths / channel-head source distances of order 100 m in soil-mantled landscapes (Montgomery and Dietrich 1988, Nature 336:232; Roering et al. 2007, Geology 35:571); citations not checked; enters only through ln(dx/x_c) (A-SUBGRID-CHANNEL-HEAD)"), "declared in orogeny.subgrid_erosion.channel_head_length_m; no citation"),
                ),
            },
        };
        // SEED-EXTRAS D7 (lane SD): sediment law entries; absent = legacy
        // variants, skipped on write.
        let sed = self.sediment_laws.as_ref();
        let land_deposition = match sed.and_then(|b| b.erosion_deposition.as_ref()) {
            None => LandDepositionLawV1::DetachmentLimitedV0,
            Some(b) => LandDepositionLawV1::ErosionDepositionV1 {
                deposition_coefficient_g: c(
                    b.deposition_coefficient_g,
                    UnitV1::Dimensionless,
                    label(b.deposition_coefficient_g, ErosionDepositionBlock::v1().deposition_coefficient_g, pa("recalled: Yuan et al. 2019 (JGR Earth Surface 124) work with a dimensionless deposition coefficient G of order 1 (xi-q model of Davy and Lague 2009, JGR 114); 1.0 declared before any run; citations not checked (A-DEPOSITION-G)"), "declared in sediment_laws.erosion_deposition.deposition_coefficient_g; no citation"),
                ),
            },
        };
        let sediment_placement = match sed.and_then(|b| b.shelf.as_ref()) {
            None => SedimentPlacementLawV1::LocalFillGlobalLayerV0,
            Some(b) => SedimentPlacementLawV1::ShelfWaveBaseV1 {
                shelf_wave_base_depth_m: c(
                    b.wave_base_depth_m,
                    UnitV1::Meter,
                    label(b.wave_base_depth_m, ShelfBlock::v1().wave_base_depth_m, pa("storm wave base / cell-mean shelf-top depth; Earth's shelf break averages ~130 m (recalled, citation not checked; partly a glacial-lowstand inheritance) and a graded shelf's mean depth is shallower; 100 m declared before any run (A-SHELF-WAVE-BASE)"), "declared in sediment_laws.shelf.wave_base_depth_m; no citation"),
                ),
            },
        };
        let earth = TwoColumnIsostasyBlock::earth_v0();
        let constants = UniverseConstantsV1 {
            schema: Cow::Borrowed(UNIVERSE_CONSTANTS_SCHEMA_V1),
            planet_radius_m: c(
                p.radius_m,
                UnitV1::Meter,
                label(p.radius_m, d.radius_m, leg("PlanetBlock default (Earth mean radius; citation not checked)"), "declared in planet.radius_m; no citation"),
            ),
            surface_gravity_m_s2: c(
                p.gravity_m_s2,
                UnitV1::MeterPerSecondSquared,
                label(p.gravity_m_s2, d.gravity_m_s2, leg("PlanetBlock default; STATED, not derived from planet mass (mass is not a v1 constant); read only by the flexure kernel, which the preview tier does not run"), "declared in planet.gravity_m_s2; no citation"),
            ),
            mantle_density_kg_m3: c(
                p.mantle_density_kg_m3,
                UnitV1::KilogramPerCubicMeter,
                label(p.mantle_density_kg_m3, d.mantle_density_kg_m3, leg("PlanetBlock default"), "declared in planet.mantle_density_kg_m3; no citation"),
            ),
            continental_crust_density_kg_m3: c(
                p.crust_density_kg_m3,
                UnitV1::KilogramPerCubicMeter,
                label(p.crust_density_kg_m3, d.crust_density_kg_m3, leg("PlanetBlock default crust_density_kg_m3 (one density for every crust kind under the single-column law)"), "declared in planet.crust_density_kg_m3; no citation"),
            ),
            oceanic_crust_density_kg_m3: two.map(|b| {
                c(
                    b.oceanic_crust_density_kg_m3,
                    UnitV1::KilogramPerCubicMeter,
                    label(b.oceanic_crust_density_kg_m3, earth.oceanic_crust_density_kg_m3, pa("recalled from Turcotte and Schubert, Geodynamics, and Carlson and Raskin 1984, Nature 311:555; citation not checked (A-OCEANIC-CRUST-DENSITY, ADR 039 C5)"), "declared in two_column_isostasy.oceanic_crust_density_kg_m3; no citation"),
                )
            }),
            seawater_density_kg_m3: match two {
                None => c(
                    legacy::SEAWATER_DENSITY_KG_M3,
                    UnitV1::KilogramPerCubicMeter,
                    cited("Turcotte and Schubert 2014, Geodynamics, 3rd ed., section 2.2 (as labelled in pf-lab A-WATER-DENSITY; pf-lab WATER_DENSITY_KG_M3; not re-checked by WS-CONST)"),
                ),
                Some(b) => c(
                    b.seawater_density_kg_m3,
                    UnitV1::KilogramPerCubicMeter,
                    label(b.seawater_density_kg_m3, earth.seawater_density_kg_m3, cited("Turcotte and Schubert 2014, Geodynamics, 3rd ed., section 2.2 (as labelled in pf-lab A-SEAWATER-DENSITY; not re-checked by WS-CONST)"), "declared in two_column_isostasy.seawater_density_kg_m3; no citation"),
                ),
            },
            continental_crust_thickness_m: c(
                legacy::CONTINENTAL_CRUST_THICKNESS_M,
                UnitV1::Meter,
                derived("equals the 35 km single-column reference column (pf-geodynamics uplift.rs CONTINENTAL_CRUST_THICKNESS_M, A-CRUST-CONTINENTAL-THICKNESS)"),
            ),
            oceanic_crust_thickness_m: c(
                legacy::OCEANIC_CRUST_THICKNESS_M,
                UnitV1::Meter,
                cited("White, McKenzie and O'Nions 1992, JGR 97:19683 (7.1 +- 0.8 km), as labelled in A-CRUST-OCEANIC-THICKNESS (pf-geodynamics uplift.rs OCEANIC_CRUST_THICKNESS_M); not re-checked by WS-CONST"),
            ),
            stream_power_erodibility_per_yr: c(
                s.erodibility,
                UnitV1::PerYear,
                label(s.erodibility, sd.erodibility, leg("SurfaceBlock default erodibility, read as K per year with m = 0.5, n = 1 (A-K-UNITS)"), "declared in surface.erodibility; no citation"),
            ),
            hillslope_diffusivity_m2_per_yr: c(
                s.hillslope_diffusivity_m2_per_yr,
                UnitV1::SquareMeterPerYear,
                label(s.hillslope_diffusivity_m2_per_yr, sd.hillslope_diffusivity_m2_per_yr, leg("SurfaceBlock default hillslope diffusivity (A-HILLSLOPE-D)"), "declared in surface.hillslope_diffusivity_m2_per_yr; no citation"),
            ),
            ocean_water_volume_m3: two.map(|b| {
                c(
                    b.ocean_water_volume_m3,
                    UnitV1::CubicMeter,
                    label(b.ocean_water_volume_m3, earth.ocean_water_volume_m3, pa("recalled from Charette and Smith 2010, Oceanography 23(2):112 (1.335e9 km3); citation not checked (A-OCEAN-VOLUME, ADR 039 C5)"), "declared in two_column_isostasy.ocean_water_volume_m3; no citation"),
                )
            }),
            sea_level_datum_m: c(
                p.sea_level_m,
                UnitV1::Meter,
                match two {
                    None => label(p.sea_level_m, d.sea_level_m, leg("PlanetBlock default prescribed sea level, constant in time (A-SEA-LEVEL-CONSTANT)"), "declared prescribed sea level planet.sea_level_m; no citation"),
                    Some(_) => derived("output datum only: the sea level solved from the ocean volume is placed here (A-SEA-LEVEL-VOLUME-PRIMARY); not a physical input"),
                },
            ),
        };
        // SEED-IRREGULAR-SEEDS (lane IS): from the bodies block.
        let craton_seeding = match self
            .continental_bodies
            .as_ref()
            .and_then(|b| b.irregular_cratons.as_ref())
        {
            None => CratonSeedingLawV1::NoneV0,
            Some(k) => {
                let v1 = IrregularCratonsBlock::v1();
                const D: &str = "declared in continental_bodies.irregular_cratons; no citation";
                CratonSeedingLawV1::IrregularCratonsV1 {
                    craton_interior_thickness_m: c(k.craton_interior_thickness_m, UnitV1::Meter, label(k.craton_interior_thickness_m, v1.craton_interior_thickness_m, pa("project assumption, model compensation: two-column isostasy puts 35 km continental crust at +1209 m (lane IS measured on conserve_v1 t = 0) with slope (rho_m - rho_c)/rho_m = 0.1515, so 29.6 km starts interiors at about +0.39 km (model specification 0.3-0.5 km; Earth bedrock land p50 381 m). Earth cratonic crust is 35-45 km (Christensen and Mooney 1995, recalled): the thinner value stands in for the missing continental lithosphere term (lane IS PREREG.txt) (A-CRATON-INTERIOR-THICKNESS)"), D)),
                    craton_rank_size_exponent: c(k.craton_rank_size_exponent, UnitV1::Dimensionless, label(k.craton_rank_size_exponent, v1.craton_rank_size_exponent, pa("project assumption: Earth landmass rank sizes recalled (Afro-Eurasia ~84, Americas ~42, Antarctica ~14, Australia ~8.6 Mkm2) give exponents 1.0-1.6 (lane IS PREREG.txt) (A-CRATON-RANK-SIZE)"), D)),
                    cratons_per_unit_share: c(k.cratons_per_unit_share, UnitV1::Dimensionless, label(k.cratons_per_unit_share, v1.cratons_per_unit_share, pa("project assumption: ~35 Archean cratons grouped into ~10 major ones (Bleeker 2003, Lithos 71:99, recalled, not checked) (lane IS PREREG.txt) (A-CRATON-COUNT)"), D)),
                    basement_roughness_std_m: c(k.basement_roughness_std_m, UnitV1::Meter, label(k.basement_roughness_std_m, v1.basement_roughness_std_m, pa("project assumption: about 230 m elevation std over interiors; intracratonic basins and arches with 0.2-1 km relief (Michigan, Williston, Hudson Bay, Congo; recalled) (lane IS PREREG.txt) (A-BASEMENT-ROUGHNESS)"), D)),
                    suture_crest_thickening_m: c(k.suture_crest_thickening_m, UnitV1::Meter, label(k.suture_crest_thickening_m, v1.suture_crest_thickening_m, pa("project assumption: old orogen roots of 40-50 km (Appalachians, Urals) against platform crust and 1-2 km bedrock crests (recalled, not checked) (lane IS PREREG.txt) (A-SUTURE-ROOT)"), D)),
                    suture_half_width_m: c(k.suture_half_width_m, UnitV1::Meter, label(k.suture_half_width_m, v1.suture_half_width_m, pa("project assumption: Appalachian and Ural belts 200-600 km wide (recalled); sigma 120 km is ~280 km FWHM (lane IS PREREG.txt) (A-SUTURE-ROOT)"), D)),
                    proximal_margin_stretch: c(k.proximal_margin_stretch, UnitV1::Dimensionless, label(k.proximal_margin_stretch, v1.proximal_margin_stretch, pa("project assumption: proximal rifted-margin domain mildly stretched (beta < ~1.5), necking zone tens of km wide (Peron-Pinvidic and Manatschal 2009, Int J Earth Sci 98:1581; Sutra et al. 2013, G3; recalled, not checked) (lane IS PREREG.txt) (A-PROXIMAL-MARGIN)"), D)),
                }
            }
        };
        // SEED-TRANSFORMS (lane TF): the plate-boundary law (absent = legacy,
        // skipped on write).
        let plate_boundary = match belt_block.and_then(|b| b.plate_boundary.as_ref()) {
            None => PlateBoundaryLawV1::NeighbourDirectionV0,
            Some(t) => {
                let v1 = PlateBoundaryBlock::v1();
                const D: &str = "declared in orogeny.rate_limited_belt.plate_boundary; no citation";
                PlateBoundaryLawV1::SmoothedNormalTransformV1 {
                    transform_obliquity_deg: c(t.transform_obliquity_deg, UnitV1::Degree, label(t.transform_obliquity_deg, v1.transform_obliquity_deg, pa("project assumption (lane TF PREREG.txt; lane VL's validation used 30 deg, 20 and 45 as sensitivity): boundaries within 30 deg of pure strike-slip behave as transforms; Earth's ridge-transform systems are mostly within a few degrees of the spreading direction (recalled) (A-TRANSFORM-OBLIQUITY)"), D)),
                    boundary_normal_smoothing_m: c(t.boundary_normal_smoothing_m, UnitV1::Meter, label(t.boundary_normal_smoothing_m, v1.boundary_normal_smoothing_m, pa("project assumption (lane TF PREREG.txt): sigma 200 km ~ 1.5 L0 cells (lane VL used sigma 1.5 L0 cells), wide enough to remove the grid staircase, narrow against plate size (A-BOUNDARY-NORMAL-SMOOTHING)"), D)),
                }
            }
        };
        // PLATE-POLES (lane PLATE-POLES): the opt-in key
        // `orogeny.rate_limited_belt.plate_boundary.rigid_group_contacts`
        // (absent / false = legacy, skipped on write).
        let rigid_group = match belt_block.and_then(|b| b.plate_boundary.as_ref()) {
            Some(t) if t.rigid_group_release && !t.rigid_group_contacts => {
                return Err(UniverseError::Unsupported {
                    law: "rigid-group/rigid-group-contacts-v1 without rigid-group-interior-v1",
                    why: "orogeny.rate_limited_belt.plate_boundary.rigid_group_release extends rigid_group_contacts; set rigid_group_contacts too",
                });
            }
            Some(t) if t.rigid_group_release => RigidGroupLawV1::RigidGroupContactsV1,
            Some(t) if t.rigid_group_contacts => RigidGroupLawV1::RigidGroupInteriorV1,
            _ => RigidGroupLawV1::SharedPoleV0,
        };
        // CLIMATE-GENESIS (lane CG): the climate forcing entry (absent =
        // legacy, skipped on write).
        let climate_forcing = match &self.climate_forcing {
            None => ClimateForcingLawV1::NoneV0,
            Some(b) if b.erosive_runoff_exponent.is_some() && b.erosive_runoff_cap_m_per_yr.is_some() => {
                return Err(UniverseError::Unsupported {
                    law: "climate-forcing erosive_runoff_exponent with erosive_runoff_cap_m_per_yr",
                    why: "the erosive-runoff response is one law: choose the sublinear exponent or the cap, not both",
                });
            }
            Some(b) => ClimateForcingLawV1::EpochClimateRunoffV1 {
                climate: b.climate.clone(),
                seed_constants_digest: b.seed_constants_digest.clone(),
                rotation_constants_digest: b.rotation_constants_digest.clone(),
                climate_epoch_my: c(b.epoch_my, UnitV1::Megayear, pa("project assumption awaiting review (lane CG DESIGN.md 4.3, S0 timing): genesis climate solve every few My, well below the 20-30 My rise-and-fall time of ranges")),
                climate_solve_grid_divisor: c(b.solve_grid_divisor, UnitV1::Dimensionless, pa("project assumption awaiting review (lane CG DESIGN.md 4.3): numerics; the 300 km fallout length is about 2 cells of the coarsened L1 grid")),
                climate_sweep_rel_tol: c(b.sweep_rel_tol, UnitV1::Dimensionless, pa("project assumption awaiting review (lane CG DESIGN.md 4.3): numerics; warm-started solve to this tolerance, water closure still checked at 1e-6")),
                climate_runoff_floor_fraction: c(b.runoff_floor_fraction, UnitV1::Dimensionless, pa("project assumption awaiting review: lane CLIMATE-RUNOFF's land_p_minus_e_floor_v1 constant f = 0.1 (same rule as the region chain)")),
                between_epochs: b.between_epochs,
                erosivity: b.erosivity,
                subgrid_cell_mean: b.subgrid_cell_mean,
                erosive_runoff: match (b.erosive_runoff_exponent, b.erosive_runoff_cap_m_per_yr) {
                    (Some(p), _) => Some(ErosiveRunoffLawV1::SublinearAboveReferenceV1 {
                        erosive_runoff_exponent: c(p, UnitV1::Dimensionless, pa("project assumption awaiting review (lane HYPSO-LEVERS 2026-10-21): sublinear response of long-term bedrock erosion to mean runoff above the reference runoff (stochastic-threshold models; Tucker 2004, Lague et al. 2005, DiBiase and Whipple 2011, recalled, not checked)")),
                    }),
                    (None, Some(cap)) => Some(ErosiveRunoffLawV1::CapV1 {
                        erosive_runoff_cap_m_per_yr: c(cap, UnitV1::MeterPerYear, pa("project assumption awaiting review (lane HYPSO-LEVERS 2026-10-21): runoff above the cap adds no erosive power (saturating response to mean runoff; recalled, not checked)")),
                    }),
                    (None, None) => None,
                },
            },
        };
        // SHELF (law-set entry, lane CLIMATE-GENESIS-2): the opt-in key
        // `orogeny.sediment_retention.shelf_wave_grading` (absent / false =
        // legacy, skipped on write).
        let shelf_grading = match orogeny.and_then(|o| o.sediment_retention.as_ref()) {
            Some(b) if b.shelf_wave_grading => ShelfGradingLawV1::WaveBaseShelfGradingV1,
            _ => ShelfGradingLawV1::NoneV0,
        };
        // INLAND-BASINS (lane INLAND-BASINS): the opt-in key
        // `orogeny.sediment_retention.closed_basin_water_balance` (absent /
        // false = legacy, skipped on write).
        // LAKE-EXCESS (lane LAKE-EXCESS): `filled_basin_alluvial_grading`
        // selects the graded variant and needs the water-balance key.
        const CB_SUSPENDED: &str = "project assumption awaiting review (lane INLAND-BASINS): suspended share of river load reaching a lake; bedload about 10 % of total load, Milliman and Meade 1983 and Turowski et al. 2010 recalled, not checked; the suspended share is focused into the deepest water (Blais and Kalff 1995 recalled)";
        const CB_EVAP: &str = "project assumption awaiting review (lane INLAND-BASINS): without a climate law, net open-water evaporation = 1.5 x the uniform runoff; Earth land means (lake evaporation ~1.2 m/yr minus land precipitation ~0.8 m/yr) / land runoff ~0.27 m/yr, Oki and Kanae 2006 and Trenberth et al. 2007 recalled, not checked";
        let closed_basin = match orogeny.and_then(|o| o.sediment_retention.as_ref()) {
            Some(b) if b.filled_basin_alluvial_grading && !b.closed_basin_water_balance => {
                return Err(UniverseError::Unsupported {
                    law: "closed-basin/filled-basin-alluvial-grade-v1 without closed-basin-water-balance-v1",
                    why: "orogeny.sediment_retention.filled_basin_alluvial_grading extends the water-balance focused fill; set closed_basin_water_balance too",
                });
            }
            Some(b) if b.filled_basin_alluvial_grading => ClosedBasinLawV1::WaterBalanceFocusedFillAlluvialGradeV1 {
                suspended_load_fraction: c(0.9, UnitV1::Dimensionless, pa(CB_SUSPENDED)),
                uniform_evaporation_to_runoff_ratio: c(1.5, UnitV1::Dimensionless, pa(CB_EVAP)),
                filled_basin_grading_depth_m: c(10.0, UnitV1::Meter, pa("project assumption awaiting review (lane LAKE-EXCESS): residual water depth over a lake's sediment fill below which the inflowing river crosses the fill as a channel (its bankfull flow depth) and builds an alluvial plain graded to the outlet instead of settling in standing water; bankfull depths of large lowland rivers roughly 5-20 m (Leopold and Maddock 1953 hydraulic geometry; lower Mississippi ~15-20 m, Rhine ~5-10 m; recalled, not checked)")),
            },
            Some(b) if b.closed_basin_water_balance => ClosedBasinLawV1::WaterBalanceFocusedFillV1 {
                suspended_load_fraction: c(0.9, UnitV1::Dimensionless, pa("project assumption awaiting review (lane INLAND-BASINS): suspended share of river load reaching a lake; bedload about 10 % of total load, Milliman and Meade 1983 and Turowski et al. 2010 recalled, not checked; the suspended share is focused into the deepest water (Blais and Kalff 1995 recalled)")),
                uniform_evaporation_to_runoff_ratio: c(1.5, UnitV1::Dimensionless, pa("project assumption awaiting review (lane INLAND-BASINS): without a climate law, net open-water evaporation = 1.5 x the uniform runoff; Earth land means (lake evaporation ~1.2 m/yr minus land precipitation ~0.8 m/yr) / land runoff ~0.27 m/yr, Oki and Kanae 2006 and Trenberth et al. 2007 recalled, not checked")),
            },
            _ => ClosedBasinLawV1::NoneV0,
        };
        // LAKE-OUTLET (lane LAKE-OUTLET): `lake_outlet_incision` selects
        // `lake-outlet-sill-incision-v1` and needs the water-balance key.
        let lake_outlet = match orogeny.and_then(|o| o.sediment_retention.as_ref()) {
            Some(b) if b.lake_outlet_incision && !b.closed_basin_water_balance => {
                return Err(UniverseError::Unsupported {
                    law: "lake-outlet/lake-outlet-sill-incision-v1 without closed-basin-water-balance-v1",
                    why: "orogeny.sediment_retention.lake_outlet_incision incises with the water balance's overflow; set closed_basin_water_balance too",
                });
            }
            // SILL-REBOUND (lane SILL-REBOUND): `lake_outlet_flexural_support`
            // selects `lake-outlet-sill-incision-flexural-support-v1` and
            // needs the outlet key.
            Some(b) if b.lake_outlet_flexural_support && !b.lake_outlet_incision => {
                return Err(UniverseError::Unsupported {
                    law: "lake-outlet/lake-outlet-sill-incision-flexural-support-v1 without lake-outlet-sill-incision-v1",
                    why: "orogeny.sediment_retention.lake_outlet_flexural_support changes how the outlet law's sill lowering is compensated; set lake_outlet_incision too",
                });
            }
            // LAKE-HOLD (lane LAKE-HOLD): `lake_outlet_path_incision` selects
            // the path variants and needs the outlet key.
            Some(b) if b.lake_outlet_path_incision && !b.lake_outlet_incision => {
                return Err(UniverseError::Unsupported {
                    law: "lake-outlet/lake-outlet-path-incision-v1 without lake-outlet-sill-incision-v1",
                    why: "orogeny.sediment_retention.lake_outlet_path_incision extends the outlet law's sill incision along the overflow path; set lake_outlet_incision too",
                });
            }
            // LAKE-HOLD-2 (lane LAKE-HOLD-2): `lake_outlet_path_through_lakes`
            // selects the through-lakes variants and needs the path key.
            Some(b) if b.lake_outlet_path_through_lakes && !(b.lake_outlet_incision && b.lake_outlet_path_incision) => {
                return Err(UniverseError::Unsupported {
                    law: "lake-outlet/lake-outlet-path-through-lakes-v1 without lake-outlet-path-incision-v1",
                    why: "orogeny.sediment_retention.lake_outlet_path_through_lakes passes the outlet path's cut through lakes on the path; set lake_outlet_incision and lake_outlet_path_incision too",
                });
            }
            Some(b) if b.lake_outlet_path_through_lakes && b.lake_outlet_flexural_support => LakeOutletLawV1::PathThroughLakesFlexuralSupportV1,
            Some(b) if b.lake_outlet_path_through_lakes => LakeOutletLawV1::PathThroughLakesV1,
            Some(b) if b.lake_outlet_incision && b.lake_outlet_path_incision && b.lake_outlet_flexural_support => LakeOutletLawV1::PathIncisionFlexuralSupportV1,
            Some(b) if b.lake_outlet_incision && b.lake_outlet_path_incision => LakeOutletLawV1::PathIncisionV1,
            Some(b) if b.lake_outlet_incision && b.lake_outlet_flexural_support => LakeOutletLawV1::SillIncisionFlexuralSupportV1,
            Some(b) if b.lake_outlet_incision => LakeOutletLawV1::SillIncisionV1,
            _ => LakeOutletLawV1::NoneV0,
        };
        // LAKE-HOLD-2: `lake_outlet_path_grade` selects the outlet-path grade
        // (validated against the outlet and threshold laws in laws.rs).
        let lake_outlet_grade = match orogeny.and_then(|o| o.sediment_retention.as_ref()) {
            Some(b) if b.lake_outlet_path_grade && !(b.lake_outlet_incision && b.lake_outlet_path_incision) => {
                return Err(UniverseError::Unsupported {
                    law: "lake-outlet-grade/lake-outlet-path-threshold-grade-v1 without lake-outlet-path-incision-v1",
                    why: "orogeny.sediment_retention.lake_outlet_path_grade sets rule 2a's grade on the overflow path; set lake_outlet_incision and lake_outlet_path_incision too",
                });
            }
            Some(b) if b.lake_outlet_path_grade => crate::laws::LakeOutletGradeLawV1::PathThresholdGradeV1,
            _ => crate::laws::LakeOutletGradeLawV1::NoneV0,
        };
        Ok(UniverseV1 {
            schema: Cow::Borrowed(UNIVERSE_SCHEMA_V1),
            law_set: LawSetV1 {
                schema: Cow::Borrowed(LAW_SET_SCHEMA_V1),
                marine_fines: self.native_evidence.as_ref().and_then(|e| e.marine_fines.as_ref()).map(|f| f.law.clone()),
                isostasy,
                sea_level,
                crust_seeding,
                oceanic_thermal: OceanicThermalLawV1::Gdh1V0,
                erosion,
                tectonic_thickening,
                land_deposition,
                sediment_placement,
                plate_carry,
                collision_sediment,
                hotspot_volcanism,
                continental_rift,
                plate_reorganization,
                plate_boundary,
                rift_flank,
                hillslope_threshold,
                fluvial_incision_threshold,
                fluvial_alluvial_channel,
                glacial_erosion,
                marine_remainder_layer,
                sediment_retention,
                subduction_trench,
                passive_margin,
                plate_suture,
                flood_basalt,
                craton_seeding,
                orogen_plateau,
                magmatism,
                stratigraphy,
                climate_forcing,
                shelf_grading,
                closed_basin,
                lake_outlet,
                lake_outlet_grade,
                rigid_group,
            },
            constants,
        })
    }
}

/// Fidelity presets: discretization only (P, review #7).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FidelityPreset {
    Preview,
    Standard,
    High,
}

impl FidelityPreset {
    pub fn block(self) -> FidelityBlock {
        match self {
            // L0 288x144, 4 layers, 1 substep
            FidelityPreset::Preview => FidelityBlock {
                grid_level: 0,
                radial_layers: 4,
                substeps: 1,
            },
            // L1 576x288, 8 layers, 2 substeps
            FidelityPreset::Standard => FidelityBlock {
                grid_level: 1,
                radial_layers: 8,
                substeps: 2,
            },
            // L2 1152x576, 16 layers, 4 substeps
            FidelityPreset::High => FidelityBlock {
                grid_level: 2,
                radial_layers: 16,
                substeps: 4,
            },
        }
    }
}

/// Scenario presets: the world's history. Held constant across fidelity axes.
/// Benchmark A/B/C map to spec §44; starting tuning, revisable with evidence.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ScenarioPreset {
    pub name: String,
    pub tectonics: TectonicsBlock,
}

fn tect(
    plate_count: u32,
    rate: f64,
    ocean: f64,
    seeds: u32,
    hotspots: u32,
    dur: f64,
    step: f64,
) -> TectonicsBlock {
    TectonicsBlock {
        plate_count,
        rate_scale_deg_per_my: rate,
        ocean_fraction: ocean,
        continental_seeds: seeds,
        hotspot_count: hotspots,
        duration_my: dur,
        step_my: step,
    }
}

impl ScenarioPreset {
    pub fn toy() -> Self {
        Self {
            name: "toy".into(),
            tectonics: tect(7, 1.0, 0.62, 4, 2, 100.0, 2.0),
        }
    }
    /// §44 Benchmark A — extreme collision: big rates, long history, land-rich.
    pub fn benchmark_a() -> Self {
        Self {
            name: "benchmark-a".into(),
            tectonics: tect(6, 4.0, 0.35, 5, 1, 300.0, 2.0),
        }
    }
    /// §44 Benchmark B — oceanic arc/island world: fast rates, water-rich.
    pub fn benchmark_b() -> Self {
        Self {
            name: "benchmark-b".into(),
            tectonics: tect(8, 3.0, 0.80, 2, 3, 200.0, 2.0),
        }
    }
    /// §44 Benchmark C — ancient eroded world: slow rates, very long history.
    pub fn benchmark_c() -> Self {
        Self {
            name: "benchmark-c".into(),
            tectonics: tect(5, 0.5, 0.55, 4, 1, 900.0, 5.0),
        }
    }

    pub fn profile(&self, fidelity: FidelityPreset, mantle: MantleBlock) -> PlanetParameterProfile {
        let mut profile = PlanetParameterProfile {
            planet: PlanetBlock::default(),
            tectonics: self.tectonics.clone(),
            mantle,
            surface: SurfaceBlock::default(),
            fidelity: fidelity.block(),
            world_character: WorldCharacterBlock::default(),
            two_column_isostasy: None,
            continental_bodies: None,
            orogeny: None,
            sediment_laws: None,
            climate_forcing: None,
            native_evidence: None,
            pinned_universe: None,
        };
        // Benchmarks A requires flexure (foreland); keep operator-visible.
        if self.name == "benchmark-a" {
            profile.mantle.compensation = Compensation::Flexure;
        }
        profile
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn valid_profile() -> PlanetParameterProfile {
        ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default())
    }

    #[test]
    fn toy_profile_validates() {
        valid_profile().validate().unwrap();
    }

    #[test]
    fn typed_runoff_profile_digest_is_pinned() {
        // 2026-09-13: typed m/yr runoff intentionally changes the profile
        // contract/digest. The additive default diffusivity stays omitted.
        let p = valid_profile();
        let json = serde_json::to_string(&p).unwrap();
        assert!(
            !json.contains("hillslope"),
            "default diffusivity leaked into the profile JSON"
        );
        let d = crate::artifact::digest_params(&p).unwrap();
        assert_eq!(
            d, "blake3:9c112dfe000114c9ae36d76e4de8ce7c8695bc174061d0a378337a0fdd86b543",
            "toy/preview params digest moved"
        );
        let mut q = valid_profile();
        q.surface.hillslope_diffusivity_m2_per_yr = 0.2;
        assert!(serde_json::to_string(&q).unwrap().contains("hillslope"));
        assert_ne!(crate::artifact::digest_params(&q).unwrap(), d);
    }

    #[test]
    fn two_column_isostasy_is_opt_in_skipped_when_absent_and_validated_when_present() {
        // Absent: not serialized, and an old JSON without the key loads to None.
        let p = valid_profile();
        assert!(p.two_column_isostasy.is_none());
        let json = serde_json::to_string(&p).unwrap();
        assert!(!json.contains("two_column"), "absent flag leaked into JSON");
        let back: PlanetParameterProfile = serde_json::from_str(&json).unwrap();
        assert_eq!(back, p);
        // Present: serialized, round-trips, moves the digest, validates.
        let mut q = valid_profile();
        q.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
        q.validate().unwrap();
        let qj = serde_json::to_string(&q).unwrap();
        assert!(qj.contains("ocean_water_volume_m3"));
        assert_eq!(serde_json::from_str::<PlanetParameterProfile>(&qj).unwrap(), q);
        assert_ne!(
            crate::artifact::digest_params(&q).unwrap(),
            crate::artifact::digest_params(&p).unwrap()
        );
        // Unknown keys inside the block are refused.
        let bad = qj.replace("\"ocean_water_volume_m3\"", "\"ocean_fraction_target\"");
        assert!(serde_json::from_str::<PlanetParameterProfile>(&bad).is_err());
        // Loud range / layering failures name their field.
        let mut r = q.clone();
        r.two_column_isostasy.as_mut().unwrap().oceanic_crust_density_kg_m3 = 3400.0;
        assert!(r.validate().is_err());
        let mut r = q.clone();
        r.two_column_isostasy.as_mut().unwrap().ocean_water_volume_m3 = -1.0;
        let e = r.validate().unwrap_err().to_string();
        assert!(e.contains("ocean_water_volume_m3"), "{e}");
        let mut r = q.clone();
        r.two_column_isostasy.as_mut().unwrap().compensation_depth_m = f64::NAN;
        assert!(r.validate().is_err());
    }

    #[test]
    fn continental_bodies_is_opt_in_skipped_when_absent_and_validated_when_present() {
        let p = valid_profile();
        assert!(p.continental_bodies.is_none());
        let json = serde_json::to_string(&p).unwrap();
        assert!(!json.contains("continental_bodies"), "absent flag leaked into JSON");
        let mut q = valid_profile();
        q.continental_bodies = Some(ContinentalBodiesBlock::v0());
        q.validate().unwrap();
        // Derived coverage: 1 - ocean_fraction (toy 0.62 -> 0.38).
        let b = q.continental_bodies.as_ref().unwrap();
        assert!((b.resolved_area_fraction(&q.tectonics) - 0.38).abs() < 1e-12);
        let qj = serde_json::to_string(&q).unwrap();
        assert!(qj.contains("rifted_margin_width_m"));
        assert!(!qj.contains("continental_area_fraction"), "derived coverage must not serialize a number");
        assert_eq!(serde_json::from_str::<PlanetParameterProfile>(&qj).unwrap(), q);
        assert_ne!(
            crate::artifact::digest_params(&q).unwrap(),
            crate::artifact::digest_params(&p).unwrap()
        );
        let bad = qj.replace("\"rifted_margin_width_m\"", "\"coast_noise\"");
        assert!(serde_json::from_str::<PlanetParameterProfile>(&bad).is_err());
        let mut r = q.clone();
        r.continental_bodies.as_mut().unwrap().continental_area_fraction = Some(1.2);
        let e = r.validate().unwrap_err().to_string();
        assert!(e.contains("continental_area_fraction"), "{e}");
        let mut r = q.clone();
        r.continental_bodies.as_mut().unwrap().rifted_margin_width_m = -1.0;
        assert!(r.validate().is_err());
        let mut r = q.clone();
        r.continental_bodies.as_mut().unwrap().continental_area_fraction = Some(0.35);
        r.validate().unwrap();
        assert_eq!(
            r.continental_bodies.as_ref().unwrap().resolved_area_fraction(&r.tectonics),
            0.35
        );
    }

    #[test]
    fn duration_must_divide_step_loudly() {
        let mut p = valid_profile();
        p.tectonics.duration_my = 101.0;
        let err = p.validate().unwrap_err();
        assert!(err.to_string().contains("integral multiple"), "{err}");
    }

    #[test]
    fn flexure_requires_thickness() {
        let mut p = valid_profile();
        p.mantle.compensation = Compensation::Flexure;
        p.mantle.effective_elastic_thickness_km = 0.0;
        assert!(p.validate().is_err());
    }

    #[test]
    fn nonfinite_rejected() {
        let mut p = valid_profile();
        p.planet.radius_m = f64::NAN;
        let err = p.validate().unwrap_err();
        assert!(err.to_string().contains("planet.radius_m"), "{err}");
    }

    #[test]
    fn benchmark_profiles_validate() {
        for s in [
            ScenarioPreset::benchmark_a(),
            ScenarioPreset::benchmark_b(),
            ScenarioPreset::benchmark_c(),
        ] {
            for f in [
                FidelityPreset::Preview,
                FidelityPreset::Standard,
                FidelityPreset::High,
            ] {
                s.profile(f, MantleBlock::default())
                    .validate()
                    .unwrap_or_else(|e| {
                        panic!("{}/{} failed validation: {e}", s.name, format!("{f:?}"))
                    });
            }
        }
    }

    #[test]
    fn waterfall_wonder_validates_and_resolves_causally() {
        let mut p = valid_profile();
        p.world_character.waterfall_wonder = 1.5;
        assert!(p.validate().is_err(), "dial above 1.0 must fail loudly");
        p.world_character.waterfall_wonder = 0.8;
        p.validate().unwrap();

        let neutral = WorldCharacterBlock::default().resolve();
        let biased = p.world_character.resolve();
        // Neutral baseline is exactly neutral.
        assert_eq!(neutral.uplift_band_sharpness, 1.0);
        assert_eq!(neutral.uplift_step_persistence, 1.0);
        // Resolution is monotone in the dial (deep-history bias grows).
        assert!(biased.resistance_contrast_magnitude > neutral.resistance_contrast_magnitude);
        assert!(biased.resistant_bed_persistence > neutral.resistant_bed_persistence);
        assert!(biased.uplift_band_sharpness > neutral.uplift_band_sharpness);
        assert!(biased.uplift_step_persistence > neutral.uplift_step_persistence);
        // Zebra-cake cap: contrast frequency stays modest even at wonder=1.
        let max = WorldCharacterBlock {
            waterfall_wonder: 1.0,
        }
        .resolve();
        assert!(max.resistance_contrast_frequency <= 0.6);
        assert!(max.resistance_contrast_magnitude <= 2.5);
        // Determinism: same dial, same resolution.
        assert_eq!(p.world_character.resolve(), p.world_character.resolve());
    }

    /// SEED-CARRY-CC (lane CC): the carry selector is absent from every
    /// existing block's bytes, selects the law-set entry when set, and is
    /// refused without the belt law.
    #[test]
    fn front_rescue_carry_is_opt_in_and_needs_the_belt() {
        let v1 = serde_json::to_string(&OrogenyBlock::v1()).unwrap();
        assert!(!v1.contains("front_rescue_carry"), "{v1}");
        let mut p = valid_profile();
        p.orogeny = Some(OrogenyBlock::v1());
        let u = p.universe().unwrap();
        assert!(u.law_set.plate_carry.is_legacy());
        assert!(!serde_json::to_string(&u.law_set).unwrap().contains("plate_carry"));
        p.orogeny = Some(OrogenyBlock::carry_v1());
        let u = p.universe().unwrap();
        assert_eq!(u.law_set.plate_carry, PlateCarryLawV1::FrontRescueCarryV1);
        let js = serde_json::to_string(&u.law_set).unwrap();
        assert!(js.contains("\"plate_carry\":\"front-rescue-carry-v1\""), "{js}");
        let back: LawSetV1 = serde_json::from_str(&js).unwrap();
        assert_eq!(back, u.law_set);
        // The law entry with the legacy stacking rule is refused.
        let mut bad = u.clone();
        bad.law_set.tectonic_thickening = crate::laws::TectonicThickeningLawV1::KinematicOverlapStackingV0;
        assert!(bad.validate().is_err());
    }

    /// SEDIMENT-TOWERS (lane ST): the collision-sediment selector is absent
    /// from every existing block's bytes, selects the law-set entry when
    /// set, round-trips, and is refused without the belt law.
    #[test]
    fn accreted_collision_sediment_is_opt_in_and_needs_the_belt() {
        for o in [OrogenyBlock::v1(), OrogenyBlock::carry_v1(), OrogenyBlock::relief_v1()] {
            let js = serde_json::to_string(&o).unwrap();
            assert!(!js.contains("accrete_collision_sediment"), "{js}");
            let mut p = valid_profile();
            p.orogeny = Some(o);
            let u = p.universe().unwrap();
            assert!(u.law_set.collision_sediment.is_legacy());
            assert!(!serde_json::to_string(&u.law_set).unwrap().contains("collision_sediment"));
        }
        for o in [OrogenyBlock::v2(), OrogenyBlock::carry_v2()] {
            let mut p = valid_profile();
            p.orogeny = Some(o);
            let u = p.universe().unwrap();
            assert_eq!(
                u.law_set.collision_sediment,
                CollisionSedimentLawV1::AccretedToBeltCrustV1
            );
            let js = serde_json::to_string(&u.law_set).unwrap();
            assert!(
                js.contains("\"collision_sediment\":\"accreted-to-belt-crust-v1\""),
                "{js}"
            );
            let back: LawSetV1 = serde_json::from_str(&js).unwrap();
            assert_eq!(back, u.law_set);
            let mut bad = u.clone();
            bad.law_set.tectonic_thickening =
                crate::laws::TectonicThickeningLawV1::KinematicOverlapStackingV0;
            bad.law_set.plate_carry = PlateCarryLawV1::SemiLagrangianNearestMaterialV0;
            assert!(bad.validate().is_err());
        }
        assert_eq!(
            OrogenyBlock::carry_v2().rate_limited_belt.unwrap().front_rescue_carry,
            true
        );
    }

    /// SEED-EXTRAS (lane EX): the plume / rift / split selectors are opt-in
    /// (absent from every existing block and law set), round-trip, and are
    /// refused in an inconsistent combination.
    #[test]
    fn seed_extras_laws_are_opt_in_and_consistent() {
        for o in [OrogenyBlock::v1(), OrogenyBlock::carry_v2(), OrogenyBlock::relief_v1()] {
            let js = serde_json::to_string(&o).unwrap();
            assert!(!js.contains("plume") && !js.contains("rift"), "{js}");
            let mut p = valid_profile();
            p.orogeny = Some(o);
            let u = p.universe().unwrap();
            let ls = serde_json::to_string(&u.law_set).unwrap();
            assert!(
                !ls.contains("hotspot_volcanism")
                    && !ls.contains("continental_rift")
                    && !ls.contains("plate_reorganization"),
                "{ls}"
            );
        }
        let mut p = valid_profile();
        p.orogeny = Some(OrogenyBlock::extras_v1());
        let u = p.universe().unwrap();
        u.validate().unwrap();
        let js = serde_json::to_string(&u.law_set).unwrap();
        assert!(js.contains("fixed-plume-chains-v1"), "{js}");
        assert!(js.contains("stretched-rift-v1"), "{js}");
        assert!(js.contains("plume-head-continental-split-v1"), "{js}");
        let back: LawSetV1 = serde_json::from_str(&js).unwrap();
        assert_eq!(back, u.law_set);
        // The split needs the plume and rift laws.
        let mut bad = u.clone();
        bad.law_set.hotspot_volcanism = HotspotVolcanismLawV1::NoneV0;
        assert!(bad.validate().is_err());
        // Plumes / rifts need the belt advance.
        let mut bad = u.clone();
        bad.law_set.plate_reorganization = PlateReorganizationLawV1::FixedPlatesV0;
        bad.law_set.tectonic_thickening =
            crate::laws::TectonicThickeningLawV1::KinematicOverlapStackingV0;
        bad.law_set.plate_carry = PlateCarryLawV1::SemiLagrangianNearestMaterialV0;
        bad.law_set.collision_sediment = CollisionSedimentLawV1::SurfaceLayerV0;
        assert!(bad.validate().is_err());
        // Breakup thickness must be below the continental thickness.
        let mut q = valid_profile();
        let mut o = OrogenyBlock::extras_v1();
        o.rate_limited_belt.as_mut().unwrap().continental_rift.as_mut().unwrap().breakup_crust_thickness_m = 30_000.0;
        q.orogeny = Some(o);
        q.universe().unwrap().validate().unwrap();
        let mut o = OrogenyBlock::extras_v1();
        o.rate_limited_belt.as_mut().unwrap().plume_volcanism.as_mut().unwrap().plume_count = 0;
        q.orogeny = Some(o);
        assert!(q.universe().and_then(|u| u.validate()).is_err());
    }

    /// SEED-RIFT-SHOULDERS (lane RS): `rift_flank` is opt-in (absent from
    /// every existing block and law set, extras_v1 included), round-trips,
    /// and is refused without the rift law or the two-column isostasy.
    #[test]
    fn rift_flank_law_is_opt_in_and_consistent() {
        for o in [OrogenyBlock::v1(), OrogenyBlock::carry_v2(), OrogenyBlock::extras_v1()] {
            let js = serde_json::to_string(&o).unwrap();
            assert!(!js.contains("rift_flank") && !js.contains("shoulder"), "{js}");
            let mut p = valid_profile();
            p.orogeny = Some(o);
            let ls = serde_json::to_string(&p.universe().unwrap().law_set).unwrap();
            assert!(!ls.contains("rift_flank"), "{ls}");
        }
        let mut p = valid_profile();
        p.mantle.compensation = Compensation::Airy;
        p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
        p.orogeny = Some(OrogenyBlock::extras_v2());
        let u = p.universe().unwrap();
        u.validate().unwrap();
        let js = serde_json::to_string(&u.law_set).unwrap();
        assert!(js.contains("shoulder-uplift-thermal-v1"), "{js}");
        let back: LawSetV1 = serde_json::from_str(&js).unwrap();
        assert_eq!(back, u.law_set);
        // Needs the rift law.
        let mut bad = u.clone();
        bad.law_set.plate_reorganization = PlateReorganizationLawV1::FixedPlatesV0;
        bad.law_set.continental_rift = ContinentalRiftLawV1::OceanicBirthV0;
        assert!(bad.validate().is_err());
        // Needs the two-column isostasy.
        let mut q = valid_profile();
        q.mantle.compensation = Compensation::Airy;
        q.orogeny = Some(OrogenyBlock::extras_v2());
        assert!(q.universe().and_then(|u| u.validate()).is_err());
    }

    /// SEED-UPLANDS (lane UP): `hillslope_threshold` is opt-in (absent from
    /// every existing block and law set) and needs a subgrid erosion law.
    #[test]
    fn hillslope_threshold_law_is_opt_in_and_consistent() {
        for o in [OrogenyBlock::v1(), OrogenyBlock::carry_v2(), OrogenyBlock::extras_v2(), OrogenyBlock::relief_v1()] {
            let js = serde_json::to_string(&o).unwrap();
            assert!(!js.contains("critical_slope"), "{js}");
            let mut p = valid_profile();
            p.mantle.compensation = Compensation::Airy;
            p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
            p.orogeny = Some(o);
            let ls = serde_json::to_string(&p.universe().unwrap().law_set).unwrap();
            assert!(!ls.contains("hillslope_threshold"), "{ls}");
        }
        let mut p = valid_profile();
        p.mantle.compensation = Compensation::Airy;
        p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
        p.orogeny = Some(OrogenyBlock::uplands_v1());
        let u = p.universe().unwrap();
        u.validate().unwrap();
        let js = serde_json::to_string(&u.law_set).unwrap();
        assert!(js.contains("critical-slope-envelope-v1"), "{js}");
        assert!(js.contains("stream-power-subgrid-nonlinear-v1"), "{js}");
        let back: LawSetV1 = serde_json::from_str(&js).unwrap();
        assert_eq!(back, u.law_set);
        // Needs a subgrid erosion law.
        let mut bad = u.clone();
        bad.law_set.erosion = ErosionLawV1::StreamPowerN1LinearDiffusionV0;
        assert!(bad.validate().is_err());
    }

    /// SEED-FLOODBASALT (lane FB): `flood_basalt` is opt-in (absent from
    /// every existing block and law set, extras_v3 included), round-trips,
    /// and is refused without the plumes or the two-column isostasy.
    #[test]
    fn flood_basalt_law_is_opt_in_and_consistent() {
        for o in [
            OrogenyBlock::v1(),
            OrogenyBlock::carry_v2(),
            OrogenyBlock::extras_v1(),
            OrogenyBlock::extras_v2(),
            OrogenyBlock::extras_v3(),
        ] {
            let js = serde_json::to_string(&o).unwrap();
            assert!(!js.contains("flood_basalt"), "{js}");
            let mut p = valid_profile();
            p.mantle.compensation = Compensation::Airy;
            p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
            p.continental_bodies = Some(ContinentalBodiesBlock::v0());
            p.orogeny = Some(o);
            let ls = serde_json::to_string(&p.universe().unwrap().law_set).unwrap();
            assert!(!ls.contains("flood_basalt") && !ls.contains("lip-swell"), "{ls}");
        }
        let mut p = valid_profile();
        p.mantle.compensation = Compensation::Airy;
        p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
        p.continental_bodies = Some(ContinentalBodiesBlock::v0());
        p.orogeny = Some(OrogenyBlock::extras_v3_flood());
        let u = p.universe().unwrap();
        u.validate().unwrap();
        let js = serde_json::to_string(&u.law_set).unwrap();
        assert!(js.contains("plume-head-lip-swell-v1"), "{js}");
        let back: LawSetV1 = serde_json::from_str(&js).unwrap();
        assert_eq!(back, u.law_set);
        // Refused without the plumes.
        let mut bad = u.clone();
        bad.law_set.hotspot_volcanism = HotspotVolcanismLawV1::NoneV0;
        bad.law_set.plate_reorganization = PlateReorganizationLawV1::FixedPlatesV0;
        assert!(bad.validate().is_err());
        // Refused under the single-column isostasy.
        let mut q = valid_profile();
        q.mantle.compensation = Compensation::Airy;
        q.continental_bodies = Some(ContinentalBodiesBlock::v0());
        q.orogeny = Some(OrogenyBlock::extras_v3_flood());
        assert!(q.universe().and_then(|u| u.validate()).is_err());
    }

    /// SEED-TRENCH (lane TN): `subduction_trench` is opt-in (absent from
    /// every existing block and law set, extras_v2 included), round-trips,
    /// and is refused without the two-column isostasy.
    #[test]
    fn subduction_trench_law_is_opt_in_and_consistent() {
        for (o, two) in [
            (OrogenyBlock::v1(), false),
            (OrogenyBlock::carry_v2(), false),
            (OrogenyBlock::extras_v1(), false),
            (OrogenyBlock::extras_v2(), true),
        ] {
            let js = serde_json::to_string(&o).unwrap();
            assert!(!js.contains("subduction_trench"), "{js}");
            let mut p = valid_profile();
            if two {
                p.mantle.compensation = Compensation::Airy;
                p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
            }
            p.orogeny = Some(o);
            let ls = serde_json::to_string(&p.universe().unwrap().law_set).unwrap();
            assert!(!ls.contains("subduction_trench") && !ls.contains("slab-flexure"), "{ls}");
            assert!(!ls.contains("passive_margin") && !js.contains("seeded_margin"), "{ls}");
        }
        let mut p = valid_profile();
        p.mantle.compensation = Compensation::Airy;
        p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
        p.orogeny = Some(OrogenyBlock::extras_v3());
        // The seeded passive-margin entry needs the continental bodies.
        assert!(p.universe().and_then(|u| u.validate()).is_err());
        p.continental_bodies = Some(ContinentalBodiesBlock::v0());
        let u = p.universe().unwrap();
        u.validate().unwrap();
        let js = serde_json::to_string(&u.law_set).unwrap();
        assert!(js.contains("slab-flexure-trench-v1"), "{js}");
        assert!(js.contains("seeded-margin-thermal-subsidence-v1"), "{js}");
        let back: LawSetV1 = serde_json::from_str(&js).unwrap();
        assert_eq!(back, u.law_set);
        let mut q = valid_profile();
        q.mantle.compensation = Compensation::Airy;
        q.orogeny = Some(OrogenyBlock::extras_v3());
        assert!(q.universe().and_then(|u| u.validate()).is_err());
    }
}
