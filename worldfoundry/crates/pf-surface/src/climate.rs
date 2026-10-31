//! CLIMATE-V0 (lane CL, 2026-10-12): the declared coarse climate tier.
//!
//! Roadmap node CLIMATE-V0: "Declared coarse climate tier: latitude
//! insolation, lapse rate, orographic lift from relief, upwind distance to
//! coast". Exit: "Named approximations published as a tiered field;
//! SEED-CLIMATE later upgrades it under the same schema."
//!
//! WHAT THIS IS: a pure, deterministic function from fields an earlier
//! stage PUBLISHED (final elevation, the sea-level datum and the marine
//! mask of the hydrology routing, on the lat/lon surface grid) plus a
//! pinned `ClimateRecipeV1` (law id `climate-coarse-v0` and constants with
//! unit, source and a code-owned envelope, `pf_contracts::laws`) to a set of
//! climate fields with a fixed schema (`pf.climate/v0`, `CLIMATE_FIELDS_V0`).
//!
//! ADMISSIBILITY (ADR 039 amendment C1, ADR 026 §1):
//! (a) every approximation is named in code (`CLIMATE_APPROXIMATIONS_V0`)
//!     and recorded in the artifact with the constants it reads;
//! (b) the only spatial inputs are fields an earlier stage produced
//!     (elevation, sea level, the routing's fill level and marine mask;
//!     grid geometry);
//! (c) there is no noise, no RNG and no hand-placed field: every
//!     latitude dependence is an analytic function of latitude with
//!     declared constants (insolation geometry, the three-cell bands), and
//!     every longitude dependence comes from the inputs (proved by
//!     `climate_v0_*` tests: zonally symmetric inputs give zonally
//!     symmetric outputs; perturbing one cell changes only its latitude
//!     ring);
//! (d) the schema names fields and their status (modelled or declared
//!     zero) so a later tier (SEED-CLIMATE) fills the same names.
//!
//! CONSERVED QUANTITIES (ADR 039 C4), declared in `CLIMATE_CONSERVED_V0`
//! and closed on every call (a breach is a typed error, never a warning):
//! - `atmospheric-water`: sum over cells of (E - P) * area is zero at the
//!   steady state (every latitude ring returns its column water to the
//!   start value; the periodic steady state is solved exactly per ring);
//! - `toa-energy`: the zonal energy balance absorbs what it emits, globally;
//! - `meridional-heat-transport`: the relaxation transport sums to zero.
//!
//! WHY A PURE FUNCTION AND NOT A SimModule: CLIMATE-V0 has no state that
//! evolves in SimTime, takes no interventions and does not feed back into
//! the genesis loop (runoff stays the constant `SurfaceBlock` value; making
//! precipitation drive erosion would change default bytes and is a
//! SEED-CLIMATE decision). It is a derived observation of the FINAL genesis
//! fields with its own law id, constants, envelope and tier identity, so it
//! is published as its own artifact citing the terrain and hydrology
//! digests. The module-shaped parts of ADR 039 §4 are still here:
//! identity (`CLIMATE_TIER_ID`, `CLIMATE_TIER_VERSION`), validity
//! (`ClimateRecipeV1::validate`, typed `OutOfEnvelope` / `Unsupported`),
//! conserved quantities and determinism. A SEED-CLIMATE that steps
//! seasons or couples into runoff should become a SimModule.
//!
//! DETERMINISM: sequential loops in a fixed order, no hash maps, no clock,
//! no randomness. Transcendental functions (sin, cos, exp, acos, powf) are
//! platform-library dependent, so byte identity is a guarantee-1 claim on a
//! pinned build only (ADR 039 §6).

use pf_contracts::laws::{ClimateLawV1, ClimateRecipeV1, UniverseError};
use pf_contracts::sim_module::{ClosureCheck, ConservedQuantity};
use pf_grid::GridSpec;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// SEED-CLIMATE (lane CM, 2026-10-12): the `seed-climate-v1` tier lives in
// `climate/seed.rs`. Nothing in `derive_climate_coarse_v0` changed; the
// additions below are the dispatcher `derive_climate`, the optional
// `seed_v1` record on the output and payload (absent, and not serialized,
// for climate-coarse-v0, so v0 artifacts are byte-identical), and the
// law-dependent tier/approximation records in `ClimateDomainV0::build`.
pub mod seed;
// ROTATION-CLIMATE (lane RC, 2026-10-20): the opt-in `seed-climate-v2` tier
// (rotation-scaled circulation bands) lives in `climate/rotation.rs`. The
// additions here are its dispatch arms and the optional `seed_v2` record
// (absent and not serialized for v0 and v1, so their bytes are unchanged).
pub mod rotation;
// CLIMATE-ARIDITY (lane CLIMATE-ARIDITY, 2026-10-21): the opt-in
// `crh-rainout-v1` moisture option (column-relative-humidity rainout
// efficiency), selected by the seed-constants set `earth_crh_v1()`.
pub mod aridity;

/// Artifact schema of a published climate tier (any tier).
pub const CLIMATE_SCHEMA_V0: &str = "pf.climate/v0";
/// Tier id of this implementation (= the law id it implements).
pub const CLIMATE_TIER_ID: &str = "climate-coarse-v0";
pub const CLIMATE_TIER_VERSION: u32 = 1;

/// Julian year (s): rates are published per year.
pub const SECONDS_PER_YEAR: f64 = 31_557_600.0;
/// Orbit samples for the annual-mean insolation (numerics, fixed by law).
pub const INSOLATION_ORBIT_SAMPLES: usize = 720;
/// Bound on the periodic fixed-point rounds per latitude ring.
pub const MAX_RING_ROUNDS: u32 = 64;
/// Relative tolerance of the ring fixed point.
pub const RING_FIXED_POINT_REL_TOL: f64 = 1.0e-12;
/// Sentinel of `aridity_index` where PET is zero (no month above 0 degC).
pub const ARIDITY_SENTINEL_PET_ZERO: f64 = -1.0;
/// Sentinel of `upwind_coast_distance_m` where the latitude circle has no
/// marine cell (no upwind ocean exists along the transport path).
pub const DISTANCE_SENTINEL_NO_UPWIND_OCEAN: f64 = -1.0;
/// UNEP (1992) aridity-index class bounds (recalled, awaiting review):
/// hyper-arid < 0.05 <= arid < 0.2 <= semi-arid < 0.5 <= dry sub-humid
/// < 0.65 <= humid.
pub const UNEP_ARIDITY_BOUNDS: [f64; 4] = [0.05, 0.2, 0.5, 0.65];
/// Thornthwaite (1948) and Willmott et al. (1985) constants, fixed by the
/// law version (recalled, awaiting review).
pub const THORNTHWAITE_MM_PER_MONTH: f64 = 16.0;
pub const THORNTHWAITE_HEAT_EXPONENT: f64 = 1.514;
pub const THORNTHWAITE_HOT_BRANCH_C: f64 = 26.5;

/// One named approximation: id, statement, and the constants it reads.
#[derive(Debug, Clone, Copy)]
pub struct ApproximationV0 {
    pub id: &'static str,
    pub statement: &'static str,
    pub constants: &'static [&'static str],
}

/// Every approximation `climate-coarse-v0` makes (ADR 039 C1 (a)).
pub const CLIMATE_APPROXIMATIONS_V0: &[ApproximationV0] = &[
    ApproximationV0 {
        id: "circular-orbit-annual-mean-insolation-v0",
        statement: "Top-of-atmosphere daily-mean insolation from latitude and solar declination, averaged uniformly over solar longitude (circular orbit: no eccentricity or precession); the first annual harmonic of the same series drives seasonality. S0 is derived as L / (4 pi d^2).",
        constants: &["star_luminosity_w", "orbital_distance_m", "axial_tilt_deg"],
    },
    ApproximationV0 {
        id: "budyko-sellers-zonal-ebm-v0",
        statement: "Zonal annual-mean energy balance at a sea-level reference: (1 - albedo) Q = A + B T + C (T - global mean T), T in degC; step ice albedo where the latitude's temperature evaluated with the mean of the two albedos is below the threshold (Budyko's ice-line condition, discrete analogue); solved from an ice-free start (ice only grows, so the first, most poleward consistent ice edge is selected among the multiple equilibria); one planetary albedo per latitude (no land/ocean contrast); elevation does not feed back into the zonal balance.",
        constants: &["olr_intercept_w_m2", "olr_slope_w_m2_k", "meridional_transport_w_m2_k", "albedo_ice_free", "albedo_ice", "ice_albedo_threshold_c"],
    },
    ApproximationV0 {
        id: "constant-environmental-lapse-rate-v0",
        statement: "Land surface temperature = zonal sea-level temperature - lapse rate * (surface - sea level), where surface = max(elevation, hydrology fill level), so lake cells sit at their water surface and only dry cells below the datum are warmed by the lapse rate; marine cells sit at the zonal sea-level temperature.",
        constants: &["lapse_rate_k_per_m"],
    },
    ApproximationV0 {
        id: "annual-mean-no-maritime-correction-v0",
        statement: "The annual-mean temperature has no land/ocean or coastal term; the ocean moderates only the seasonal amplitude.",
        constants: &[],
    },
    ApproximationV0 {
        id: "three-cell-zonal-wind-bands-v0",
        statement: "Prevailing surface wind u = -U_band sin(6 |lat|): easterly trades below 30 deg, westerlies 30-60 deg, polar easterlies above 60 deg, calm at 0, 30, 60 and 90 deg; the meridional component is declared zero; no monsoon, seasonal shift or topographic steering. Valid only for an Earth-like rotation rate (envelope on rotation_period_s).",
        constants: &["trade_wind_speed_m_s", "westerly_wind_speed_m_s", "polar_easterly_speed_m_s", "rotation_period_s"],
    },
    ApproximationV0 {
        id: "upwind-coast-distance-along-zonal-path-v0",
        statement: "Distance from a cell centre to the nearest marine cell upwind along its latitude circle (half a cell for the first land cell); marine = the hydrology routing ocean (largest connected sub-sea component); lakes and disconnected sub-sea basins are not moisture sources; sentinel -1 when the circle has no marine cell.",
        constants: &[],
    },
    ApproximationV0 {
        id: "zonal-moisture-transport-v0",
        statement: "Column water is advected along the latitude circle in the prevailing wind direction at one effective speed; over the ocean it relaxes toward rh * saturation at the zonal sea-level temperature; background rainout is column / residence time; the periodic steady state is solved exactly per ring. No meridional moisture transport: each latitude ring closes its own budget (no ITCZ convergence, no export from subtropical oceans, no poleward supply to continental interiors or polar caps).",
        constants: &["moisture_transport_speed_m_s", "rainout_residence_s", "ocean_evaporation_relaxation_s", "ocean_column_relative_humidity", "column_water_saturation_ref_kg_m2", "column_water_reference_temperature_c", "clausius_clapeyron_rate_per_k"],
    },
    ApproximationV0 {
        id: "three-cell-rainout-modulation-v0",
        statement: "Background rainout rate scaled by 1 + k cos(6 lat): enhanced at the rising branches (equator, 60 deg), suppressed at the sinking branches (30 deg, poles) of the same three-cell idealization as the wind bands.",
        constants: &["circulation_rainout_modulation"],
    },
    ApproximationV0 {
        id: "saturation-rainout-orographic-v0",
        statement: "On entering a cell, column water above saturation at the cell's surface temperature rains out at once. Saturation scales as exp(cc (T - T_ref)), and the lapse rate cools elevated cells, so lift over relief precipitates on the windward side and the depleted column continues downwind (rain shadow).",
        constants: &["lapse_rate_k_per_m", "clausius_clapeyron_rate_per_k", "column_water_saturation_ref_kg_m2", "column_water_reference_temperature_c"],
    },
    ApproximationV0 {
        id: "fixed-land-recycling-v0",
        statement: "Land evapotranspiration = recycling ratio * background rainout; no soil-moisture store, no orographic recycling.",
        constants: &["land_recycling_ratio"],
    },
    ApproximationV0 {
        id: "continentality-damped-seasonal-response-v0",
        statement: "Seasonal half range = m * A1 * (1 - albedo) / (B + C): the quasi-equilibrium response to the first annual insolation harmonic A1, with m = d + (1 - d) (1 - exp(-x / L)) for upwind coast distance x (m = d on marine cells, 1 on circles without ocean).",
        constants: &["ocean_seasonal_damping", "continentality_length_m", "olr_slope_w_m2_k", "meridional_transport_w_m2_k"],
    },
    ApproximationV0 {
        id: "thornthwaite-1948-unadjusted-pet-v0",
        statement: "Potential evapotranspiration from twelve monthly temperatures sampled from mean +- half range * cos; Thornthwaite heat-index formula, the Willmott et al. (1985) branch above 26.5 degC, no day-length or month-length adjustment. Known to underestimate PET in arid climates. Computed on every cell; on marine cells it is not a land aridity.",
        constants: &[],
    },
    ApproximationV0 {
        id: "unep-aridity-index-v0",
        statement: "Aridity index = P / PET (UNEP classes 0.05 / 0.2 / 0.5 / 0.65); sentinel -1 where PET is zero.",
        constants: &[],
    },
];

/// Status of a schema field in a given tier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum FieldStatusV0 {
    /// Computed by this tier's named approximations.
    Modelled,
    /// Present for schema stability; this tier declares it identically zero
    /// (a later tier fills it under the same name).
    DeclaredZero,
}

/// One field of the `pf.climate/v0` schema (f32, row-major, on the
/// terrain's surface domain).
#[derive(Debug, Clone, Copy)]
pub struct ClimateFieldSpecV0 {
    pub name: &'static str,
    pub units: &'static str,
    pub meaning: &'static str,
    pub status: FieldStatusV0,
    pub approximation: &'static str,
    pub sentinel: Option<(f64, &'static str)>,
}

/// The schema: names, units and meanings are what a client reads. A later
/// tier keeps every name and unit; it may change a status from
/// declared-zero to modelled and add fields, never rename or re-unit one.
pub const CLIMATE_FIELDS_V0: &[ClimateFieldSpecV0] = &[
    ClimateFieldSpecV0 { name: "temperature_mean_c", units: "degC", meaning: "annual-mean near-surface air temperature", status: FieldStatusV0::Modelled, approximation: "budyko-sellers-zonal-ebm-v0 + constant-environmental-lapse-rate-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "temperature_half_range_k", units: "K", meaning: "seasonal half range (warmest-month mean = mean + half range, coldest = mean - half range)", status: FieldStatusV0::Modelled, approximation: "continentality-damped-seasonal-response-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "precipitation_mm_per_yr", units: "mm/yr", meaning: "annual precipitation (kg/m2/yr of water)", status: FieldStatusV0::Modelled, approximation: "zonal-moisture-transport-v0 + saturation-rainout-orographic-v0 + three-cell-rainout-modulation-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "precipitation_orographic_mm_per_yr", units: "mm/yr", meaning: "part of precipitation_mm_per_yr from saturation rainout (lift over relief)", status: FieldStatusV0::Modelled, approximation: "saturation-rainout-orographic-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "evaporation_mm_per_yr", units: "mm/yr", meaning: "ocean evaporation and land evapotranspiration into the column", status: FieldStatusV0::Modelled, approximation: "zonal-moisture-transport-v0 + fixed-land-recycling-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "pet_mm_per_yr", units: "mm/yr", meaning: "potential evapotranspiration", status: FieldStatusV0::Modelled, approximation: "thornthwaite-1948-unadjusted-pet-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "aridity_index", units: "1", meaning: "P / PET", status: FieldStatusV0::Modelled, approximation: "unep-aridity-index-v0", sentinel: Some((ARIDITY_SENTINEL_PET_ZERO, "PET is zero (no month above 0 degC): aridity undefined")) },
    ClimateFieldSpecV0 { name: "wind_u_m_s", units: "m/s", meaning: "prevailing surface wind, eastward component", status: FieldStatusV0::Modelled, approximation: "three-cell-zonal-wind-bands-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "wind_v_m_s", units: "m/s", meaning: "prevailing surface wind, northward component", status: FieldStatusV0::DeclaredZero, approximation: "three-cell-zonal-wind-bands-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "upwind_coast_distance_m", units: "m", meaning: "distance to the nearest ocean upwind along the moisture transport path", status: FieldStatusV0::Modelled, approximation: "upwind-coast-distance-along-zonal-path-v0", sentinel: Some((DISTANCE_SENTINEL_NO_UPWIND_OCEAN, "no ocean upwind along the transport path")) },
    ClimateFieldSpecV0 { name: "insolation_toa_w_m2", units: "W/m2", meaning: "annual-mean top-of-atmosphere insolation", status: FieldStatusV0::Modelled, approximation: "circular-orbit-annual-mean-insolation-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "column_water_kg_m2", units: "kg/m2", meaning: "annual-mean column water vapour along the transport path", status: FieldStatusV0::Modelled, approximation: "zonal-moisture-transport-v0", sentinel: None },
];

/// Conserved quantities (ADR 039 C4).
pub const CLIMATE_CONSERVED_V0: &[ConservedQuantity] = &[
    ConservedQuantity {
        id: "atmospheric-water",
        unit: "kg/yr",
        ledger: "sum over cells of (evaporation - precipitation) * cell area; zero when every latitude ring's column water returns to its start value",
        rel_tol: 1.0e-9,
        abs_tol: 1.0,
        gross: "sum over cells of precipitation * cell area",
        applies: None,
        note: "closes per latitude ring (no meridional moisture transport); the ring steady state is solved exactly, so the residual is rounding",
    },
    ConservedQuantity {
        id: "toa-energy",
        unit: "W/m2",
        ledger: "area-weighted global mean of absorbed (1 - albedo) Q minus emitted A + B T",
        rel_tol: 1.0e-9,
        abs_tol: 1.0e-9,
        gross: "global mean absorbed insolation",
        applies: None,
        note: "zonal sea-level-reference balance only; elevation and the seasonal cycle do not enter",
    },
    ConservedQuantity {
        id: "meridional-heat-transport",
        unit: "W/m2",
        ledger: "area-weighted global mean of the relaxation transport C (global mean T - T)",
        rel_tol: 1.0e-9,
        abs_tol: 1.0e-9,
        gross: "area-weighted global mean of |C (global mean T - T)|",
        applies: None,
        note: "zero-sum by construction of the global mean; checked numerically",
    },
];

/// Typed failures. Recipe refusals carry the WS-CONST typed error.
#[derive(Debug, Clone, PartialEq)]
pub enum ClimateError {
    Recipe(UniverseError),
    Inputs(String),
    NoConvergence { stage: &'static str, row: usize, rounds: u32 },
    NonFinite { field: &'static str, index: usize },
    ClosureFailed(ClosureCheck),
}

impl std::fmt::Display for ClimateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ClimateError::Recipe(e) => write!(f, "climate recipe refused: {e}"),
            ClimateError::Inputs(w) => write!(f, "climate inputs refused: {w}"),
            ClimateError::NoConvergence { stage, row, rounds } => {
                write!(f, "{stage} did not converge on row {row} after {rounds} rounds")
            }
            ClimateError::NonFinite { field, index } => {
                write!(f, "nonfinite {field} at cell {index}")
            }
            ClimateError::ClosureFailed(c) => write!(
                f,
                "conserved quantity {} does not close at {}: residual {} {} > allowed {}",
                c.quantity, c.stage, c.residual, c.unit, c.allowed
            ),
        }
    }
}

impl std::error::Error for ClimateError {}

impl From<UniverseError> for ClimateError {
    fn from(e: UniverseError) -> Self {
        ClimateError::Recipe(e)
    }
}

/// The published fields an earlier stage produced, and nothing else.
pub struct ClimateInputs<'a> {
    pub grid: &'a GridSpec,
    /// Final surface elevation (m), as published (f32 widened to f64).
    pub elevation_m: &'a [f64],
    /// Water surface of the hydrology routing (m): the depression-filled
    /// level (`fill_level_m`, never below elevation). Lake cells sit at
    /// their water surface, not at their floor.
    pub fill_level_m: &'a [f64],
    /// The sea-level datum the terrain was published against (m).
    pub sea_level_m: f64,
    /// Marine mask: the hydrology routing ocean (largest connected
    /// sub-sea component), as published through the marine water body.
    pub ocean: &'a [bool],
}

/// The climate fields, one value per cell, in schema order.
#[derive(Debug, Clone, PartialEq)]
pub struct ClimateFieldsV0 {
    pub temperature_mean_c: Vec<f64>,
    pub temperature_half_range_k: Vec<f64>,
    pub precipitation_mm_per_yr: Vec<f64>,
    pub precipitation_orographic_mm_per_yr: Vec<f64>,
    pub evaporation_mm_per_yr: Vec<f64>,
    pub pet_mm_per_yr: Vec<f64>,
    pub aridity_index: Vec<f64>,
    pub wind_u_m_s: Vec<f64>,
    pub wind_v_m_s: Vec<f64>,
    pub upwind_coast_distance_m: Vec<f64>,
    pub insolation_toa_w_m2: Vec<f64>,
    pub column_water_kg_m2: Vec<f64>,
}

impl ClimateFieldsV0 {
    /// (name, values) in `CLIMATE_FIELDS_V0` order.
    pub fn named(&self) -> [(&'static str, &[f64]); 12] {
        [
            ("temperature_mean_c", &self.temperature_mean_c),
            ("temperature_half_range_k", &self.temperature_half_range_k),
            ("precipitation_mm_per_yr", &self.precipitation_mm_per_yr),
            ("precipitation_orographic_mm_per_yr", &self.precipitation_orographic_mm_per_yr),
            ("evaporation_mm_per_yr", &self.evaporation_mm_per_yr),
            ("pet_mm_per_yr", &self.pet_mm_per_yr),
            ("aridity_index", &self.aridity_index),
            ("wind_u_m_s", &self.wind_u_m_s),
            ("wind_v_m_s", &self.wind_v_m_s),
            ("upwind_coast_distance_m", &self.upwind_coast_distance_m),
            ("insolation_toa_w_m2", &self.insolation_toa_w_m2),
            ("column_water_kg_m2", &self.column_water_kg_m2),
        ]
    }
}

/// Per-latitude-row record (published in the domain payload).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ZonalRowV0 {
    pub lat_deg: f64,
    pub insolation_annual_w_m2: f64,
    pub insolation_first_harmonic_w_m2: f64,
    pub albedo: f64,
    pub sea_level_temperature_c: f64,
    pub wind_u_m_s: f64,
    pub rainout_factor: f64,
    pub ocean_cells: u32,
    pub ring_rounds: u32,
    /// Column water start-minus-end of the final ring pass (kg/m2).
    pub ring_residual_kg_m2: f64,
}

/// Land fraction per UNEP aridity class (area-weighted, land cells only).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AridityClassesV0 {
    pub hyper_arid: f64,
    pub arid: f64,
    pub semi_arid: f64,
    pub dry_subhumid: f64,
    pub humid: f64,
    pub pet_zero: f64,
}

/// Global summary numbers (area-weighted with true cell areas).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClimateSummaryV0 {
    pub solar_constant_w_m2: f64,
    pub global_mean_insolation_w_m2: f64,
    pub ebm_global_mean_sea_level_temperature_c: f64,
    pub ebm_iterations: u32,
    /// Equatorward-most ice-covered row centre per hemisphere (None: no ice).
    pub ice_edge_lat_deg_north: Option<f64>,
    pub ice_edge_lat_deg_south: Option<f64>,
    pub global_mean_temperature_c: f64,
    pub land_mean_temperature_c: Option<f64>,
    pub ocean_mean_temperature_c: Option<f64>,
    pub land_area_fraction: f64,
    pub global_mean_precipitation_mm_per_yr: f64,
    pub land_mean_precipitation_mm_per_yr: Option<f64>,
    pub ocean_mean_precipitation_mm_per_yr: Option<f64>,
    pub max_precipitation_mm_per_yr: f64,
    /// Orographic share of land precipitation (None: no land precipitation).
    pub land_orographic_share: Option<f64>,
    pub rows_without_ocean: u32,
    pub land_aridity_classes: AridityClassesV0,
}

/// Everything `derive_climate_coarse_v0` returns.
#[derive(Debug, Clone, PartialEq)]
pub struct ClimateOutputV0 {
    pub fields: ClimateFieldsV0,
    pub rows: Vec<ZonalRowV0>,
    pub summary: ClimateSummaryV0,
    pub closures: Vec<ClosureCheck>,
    /// SEED-CLIMATE: the extra constants and solver record of a
    /// `seed-climate-v1` run (None for climate-coarse-v0).
    pub seed_v1: Option<seed::SeedClimateRecordV1>,
    /// ROTATION-CLIMATE: the circulation record of a `seed-climate-v2` run
    /// (None for the other tiers).
    pub seed_v2: Option<rotation::RotationClimateRecordV2>,
}

/// SEED-CLIMATE (lane CM): run the tier the recipe's law names, with the
/// default `seed-climate-v1` extra constants. climate-coarse-v0 goes to
/// `derive_climate_coarse_v0` unchanged.
pub fn derive_climate(
    inputs: &ClimateInputs<'_>,
    recipe: &ClimateRecipeV1,
) -> Result<ClimateOutputV0, ClimateError> {
    match recipe.law {
        ClimateLawV1::ClimateCoarseV0 => derive_climate_coarse_v0(inputs, recipe),
        ClimateLawV1::SeedClimateV1 => seed::derive_seed_climate_v1(
            inputs,
            recipe,
            &seed::SeedClimateConstantsV1::earth_v1(),
        ),
        ClimateLawV1::SeedClimateV2 => rotation::derive_seed_climate_v2(
            inputs,
            recipe,
            &seed::SeedClimateConstantsV1::earth_v1(),
            &rotation::RotationClimateConstantsV2::earth_v2(),
        ),
    }
}

/// Tier id and version of the implementation a law selects.
pub fn climate_tier_of(law: ClimateLawV1) -> (&'static str, u32) {
    match law {
        ClimateLawV1::ClimateCoarseV0 => (CLIMATE_TIER_ID, CLIMATE_TIER_VERSION),
        ClimateLawV1::SeedClimateV1 => (seed::SEED_CLIMATE_TIER_ID, seed::SEED_CLIMATE_TIER_VERSION),
        ClimateLawV1::SeedClimateV2 => (rotation::SEED_CLIMATE_V2_TIER_ID, rotation::SEED_CLIMATE_V2_TIER_VERSION),
    }
}

// ---------------------------------------------------------------------------
// Physics pieces (public for known-answer tests).

/// Daily-mean top-of-atmosphere insolation (W/m2) at latitude `lat` for
/// solar declination `decl` (radians), solar constant `s0`.
pub fn daily_mean_insolation_w_m2(s0: f64, lat: f64, decl: f64) -> f64 {
    let x = -lat.tan() * decl.tan();
    let h0 = if x >= 1.0 {
        0.0
    } else if x <= -1.0 {
        PI
    } else {
        x.acos()
    };
    (s0 / PI) * (h0 * lat.sin() * decl.sin() + lat.cos() * decl.cos() * h0.sin())
}

/// Annual mean and first-harmonic amplitude of daily-mean insolation over
/// a circular orbit (solar longitude sampled uniformly).
pub fn annual_insolation_w_m2(s0: f64, lat_deg: f64, tilt_deg: f64) -> (f64, f64) {
    let lat = lat_deg.to_radians();
    let eps = tilt_deg.to_radians();
    let n = INSOLATION_ORBIT_SAMPLES;
    let (mut sum, mut c, mut s) = (0.0, 0.0, 0.0);
    for k in 0..n {
        let lam = 2.0 * PI * (k as f64 + 0.5) / n as f64;
        let decl = (eps.sin() * lam.sin()).asin();
        let q = daily_mean_insolation_w_m2(s0, lat, decl);
        sum += q;
        c += q * lam.cos();
        s += q * lam.sin();
    }
    let mean = sum / n as f64;
    let a1 = 2.0 * (c * c + s * s).sqrt() / n as f64;
    (mean, a1)
}

/// Prevailing zonal wind (m/s, eastward positive): `-U_band sin(6 |lat|)`.
pub fn zonal_wind_m_s(lat_deg: f64, trade: f64, westerly: f64, polar: f64) -> f64 {
    let a = lat_deg.abs();
    let amp = if a < 30.0 {
        trade
    } else if a < 60.0 {
        westerly
    } else {
        polar
    };
    -amp * (6.0 * a.to_radians()).sin()
}

/// Background rainout factor of the three-cell circulation.
pub fn rainout_factor(lat_deg: f64, k: f64) -> f64 {
    1.0 + k * (6.0 * lat_deg.to_radians()).cos()
}

/// Saturated column water (kg/m2) at temperature `t_c`.
pub fn saturation_column_kg_m2(w_ref: f64, t_ref_c: f64, cc: f64, t_c: f64) -> f64 {
    w_ref * (cc * (t_c - t_ref_c)).exp()
}

/// Thornthwaite (1948) PET (mm/yr) from twelve monthly means sampled from
/// `mean +- half_range cos`, unadjusted for day and month length, with the
/// Willmott et al. (1985) branch at and above 26.5 degC.
pub fn thornthwaite_pet_mm_per_yr(mean_c: f64, half_range_k: f64) -> f64 {
    let months: [f64; 12] = std::array::from_fn(|m| {
        mean_c + half_range_k * (2.0 * PI * (m as f64 + 0.5) / 12.0).cos()
    });
    let heat: f64 = months
        .iter()
        .filter(|&&t| t > 0.0)
        .map(|&t| (t / 5.0).powf(THORNTHWAITE_HEAT_EXPONENT))
        .sum();
    if heat <= 0.0 {
        return 0.0;
    }
    let a = 6.75e-7 * heat.powi(3) - 7.71e-5 * heat.powi(2) + 1.792e-2 * heat + 0.49239;
    months
        .iter()
        .map(|&t| {
            if t <= 0.0 {
                0.0
            } else if t < THORNTHWAITE_HOT_BRANCH_C {
                THORNTHWAITE_MM_PER_MONTH * (10.0 * t / heat).powf(a)
            } else {
                -415.85 + 32.24 * t - 0.43 * t * t
            }
        })
        .sum()
}

/// Per-cell transit coefficients of the moisture ring.
#[derive(Clone, Copy)]
struct CellTransit {
    /// Saturated column at the cell's surface temperature.
    w_sat: f64,
    /// dW/dt = a - b W during transit.
    a: f64,
    b: f64,
    dt: f64,
    ocean: bool,
}

/// One pass around the ring from `g`: returns (W_end, alpha, beta) with
/// W_end = alpha + beta * g for the clip pattern this pass took.
fn ring_lap(cells: &[CellTransit], g: f64) -> (f64, f64, f64) {
    let (mut w, mut alpha, mut beta) = (g, 0.0, 1.0);
    for c in cells {
        if w > c.w_sat {
            w = c.w_sat;
            alpha = c.w_sat;
            beta = 0.0;
        }
        let e = (-c.b * c.dt).exp();
        let ws = c.a / c.b;
        w = ws + (w - ws) * e;
        alpha = ws + (alpha - ws) * e;
        beta *= e;
    }
    (w, alpha, beta)
}

/// Area-weighted mean over a mask (None if the mask has no area).
fn masked_mean(v: &[f64], areas: &[f64], mask: impl Fn(usize) -> bool) -> Option<f64> {
    let (mut s, mut a) = (0.0, 0.0);
    for c in 0..v.len() {
        if mask(c) {
            s += v[c] * areas[c];
            a += areas[c];
        }
    }
    (a > 0.0).then(|| s / a)
}

/// The declared coarse climate tier. Pure: identical inputs give
/// identical outputs; the only spatial inputs are `inputs`.
pub fn derive_climate_coarse_v0(
    inputs: &ClimateInputs<'_>,
    recipe: &ClimateRecipeV1,
) -> Result<ClimateOutputV0, ClimateError> {
    recipe.validate()?;
    match recipe.law {
        ClimateLawV1::ClimateCoarseV0 => {}
        ClimateLawV1::SeedClimateV1 => {
            return Err(ClimateError::Recipe(UniverseError::Unsupported {
                law: "climate/seed-climate-v1",
                why: "not implemented by climate-coarse-v0",
            }))
        }
        ClimateLawV1::SeedClimateV2 => {
            return Err(ClimateError::Recipe(UniverseError::Unsupported {
                law: "climate/seed-climate-v2",
                why: "not implemented by climate-coarse-v0",
            }))
        }
    }
    let grid = inputs.grid;
    let (nlat, nlon, n) = (grid.nlat(), grid.nlon(), grid.cells());
    if inputs.elevation_m.len() != n || inputs.ocean.len() != n || inputs.fill_level_m.len() != n {
        return Err(ClimateError::Inputs(format!(
            "elevation ({}), fill level ({}) and ocean ({}) must all have {n} cells",
            inputs.elevation_m.len(),
            inputs.fill_level_m.len(),
            inputs.ocean.len()
        )));
    }
    if !inputs.sea_level_m.is_finite() {
        return Err(ClimateError::Inputs("sea level is not finite".into()));
    }
    if let Some(c) = inputs.elevation_m.iter().position(|z| !z.is_finite()) {
        return Err(ClimateError::NonFinite { field: "elevation_m", index: c });
    }
    if let Some(c) = inputs.fill_level_m.iter().position(|z| !z.is_finite()) {
        return Err(ClimateError::NonFinite { field: "fill_level_m", index: c });
    }
    let k = &recipe.constants;
    let s0 = k.solar_constant_w_m2();
    let (aa, bb, cc_t) = (
        k.olr_intercept_w_m2.value,
        k.olr_slope_w_m2_k.value,
        k.meridional_transport_w_m2_k.value,
    );
    let (alb_free, alb_ice, t_ice) = (
        k.albedo_ice_free.value,
        k.albedo_ice.value,
        k.ice_albedo_threshold_c.value,
    );
    let lapse = k.lapse_rate_k_per_m.value;
    let (w_ref, t_ref, cc) = (
        k.column_water_saturation_ref_kg_m2.value,
        k.column_water_reference_temperature_c.value,
        k.clausius_clapeyron_rate_per_k.value,
    );
    let rh = k.ocean_column_relative_humidity.value;
    let tau_e = k.ocean_evaporation_relaxation_s.value;
    let tau_r = k.rainout_residence_s.value;
    let u_eff = k.moisture_transport_speed_m_s.value;
    let recyc = k.land_recycling_ratio.value;
    let k_circ = k.circulation_rainout_modulation.value;
    let damp = k.ocean_seasonal_damping.value;
    let l_cont = k.continentality_length_m.value;

    // --- geometry: true cell areas; transit length so length * row width
    // equals the cell area exactly (budgets close in true areas).
    let radius = grid.radius_m();
    let row_width = radius * PI / nlat as f64;
    let mut lat = Vec::with_capacity(nlat);
    let mut row_area = Vec::with_capacity(nlat);
    for i in 0..nlat {
        lat.push(grid.center_lat_deg(i).map_err(|e| ClimateError::Inputs(e.to_string()))?);
        row_area.push(grid.cell_area_m2(i).map_err(|e| ClimateError::Inputs(e.to_string()))?);
    }
    let total_area: f64 = row_area.iter().map(|a| a * nlon as f64).sum();
    let weight: Vec<f64> = row_area.iter().map(|a| a * nlon as f64 / total_area).collect();
    let areas: Vec<f64> = (0..n).map(|c| row_area[c / nlon]).collect();

    // --- insolation (circular-orbit-annual-mean-insolation-v0).
    let (mut q, mut a1) = (Vec::with_capacity(nlat), Vec::with_capacity(nlat));
    for &l in &lat {
        let (m, h) = annual_insolation_w_m2(s0, l, k.axial_tilt_deg.value);
        q.push(m);
        a1.push(h);
    }

    // --- zonal EBM (budyko-sellers-zonal-ebm-v0), warm start, ice grows.
    let mut alb = vec![alb_free; nlat];
    let mut t_sl = vec![0.0; nlat];
    let mut t_bar: f64;
    let mut ebm_iterations = 0u32;
    loop {
        ebm_iterations += 1;
        if ebm_iterations as usize > nlat + 2 {
            return Err(ClimateError::NoConvergence {
                stage: "zonal energy balance (ice-albedo pattern)",
                row: 0,
                rounds: ebm_iterations,
            });
        }
        t_bar = (0..nlat)
            .map(|i| weight[i] * ((1.0 - alb[i]) * q[i] - aa))
            .sum::<f64>()
            / bb;
        for i in 0..nlat {
            t_sl[i] = ((1.0 - alb[i]) * q[i] - aa + cc_t * t_bar) / (bb + cc_t);
        }
        // Budyko's ice-line condition: a latitude is ice-covered when its
        // temperature evaluated with the MEAN of the two albedos is below
        // the threshold (the edge carries the mean albedo). Depends on the
        // albedo pattern only through the global mean, so ice only grows
        // from the ice-free start and the loop is monotone.
        let alb_edge = 0.5 * (alb_free + alb_ice);
        let mut changed = false;
        for i in 0..nlat {
            let t_edge = ((1.0 - alb_edge) * q[i] - aa + cc_t * t_bar) / (bb + cc_t);
            let want = if t_edge < t_ice { alb_ice } else { alb_free };
            if want != alb[i] {
                alb[i] = want;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let absorbed: f64 = (0..nlat).map(|i| weight[i] * (1.0 - alb[i]) * q[i]).sum();
    let emitted: f64 = (0..nlat).map(|i| weight[i] * (aa + bb * t_sl[i])).sum();
    let transport: f64 = (0..nlat).map(|i| weight[i] * cc_t * (t_bar - t_sl[i])).sum();
    let transport_gross: f64 = (0..nlat)
        .map(|i| weight[i] * (cc_t * (t_bar - t_sl[i])).abs())
        .sum();

    // --- winds and rainout factor (three-cell bands).
    let (trade, west, polar) = (
        k.trade_wind_speed_m_s.value,
        k.westerly_wind_speed_m_s.value,
        k.polar_easterly_speed_m_s.value,
    );
    let u_row: Vec<f64> = lat.iter().map(|&l| zonal_wind_m_s(l, trade, west, polar)).collect();
    let f_row: Vec<f64> = lat.iter().map(|&l| rainout_factor(l, k_circ)).collect();

    // --- surface temperature (lapse rate; marine cells at the zonal T).
    let sea = inputs.sea_level_m;
    let mut t_mean = vec![0.0; n];
    for c in 0..n {
        let i = c / nlon;
        t_mean[c] = if inputs.ocean[c] {
            t_sl[i]
        } else {
            let surface = inputs.elevation_m[c].max(inputs.fill_level_m[c]);
            t_sl[i] - lapse * (surface - sea)
        };
    }

    let mut fields = ClimateFieldsV0 {
        temperature_mean_c: t_mean,
        temperature_half_range_k: vec![0.0; n],
        precipitation_mm_per_yr: vec![0.0; n],
        precipitation_orographic_mm_per_yr: vec![0.0; n],
        evaporation_mm_per_yr: vec![0.0; n],
        pet_mm_per_yr: vec![0.0; n],
        aridity_index: vec![0.0; n],
        wind_u_m_s: vec![0.0; n],
        wind_v_m_s: vec![0.0; n],
        upwind_coast_distance_m: vec![DISTANCE_SENTINEL_NO_UPWIND_OCEAN; n],
        insolation_toa_w_m2: vec![0.0; n],
        column_water_kg_m2: vec![0.0; n],
    };
    let mut rows = Vec::with_capacity(nlat);
    let mut rows_without_ocean = 0u32;

    for i in 0..nlat {
        let ds = row_area[i] / row_width;
        let dt = ds / u_eff;
        let u = u_row[i];
        if u == 0.0 || !u.is_finite() {
            return Err(ClimateError::Inputs(format!(
                "row {i} (lat {}) has no defined wind direction",
                lat[i]
            )));
        }
        // Downwind order around the ring, starting at column 0.
        let order: Vec<usize> = (0..nlon)
            .map(|s| {
                if u > 0.0 {
                    i * nlon + s
                } else {
                    i * nlon + (nlon - s) % nlon
                }
            })
            .collect();
        let ocean_cells = order.iter().filter(|&&c| inputs.ocean[c]).count() as u32;
        if ocean_cells == 0 {
            rows_without_ocean += 1;
        }

        // Upwind coast distance: two downwind laps.
        let mut run: Option<f64> = None;
        for _ in 0..2 {
            for &c in &order {
                if inputs.ocean[c] {
                    fields.upwind_coast_distance_m[c] = 0.0;
                    run = Some(0.0);
                } else if let Some(r) = run {
                    fields.upwind_coast_distance_m[c] = r + 0.5 * ds;
                    run = Some(r + ds);
                }
            }
        }

        // Moisture transit coefficients.
        let w_eq = rh * saturation_column_kg_m2(w_ref, t_ref, cc, t_sl[i]);
        let f = f_row[i];
        let cells: Vec<CellTransit> = order
            .iter()
            .map(|&c| {
                let w_sat = saturation_column_kg_m2(w_ref, t_ref, cc, fields.temperature_mean_c[c]);
                if inputs.ocean[c] {
                    CellTransit { w_sat, a: w_eq / tau_e, b: 1.0 / tau_e + f / tau_r, dt, ocean: true }
                } else {
                    CellTransit { w_sat, a: 0.0, b: (1.0 - recyc) * f / tau_r, dt, ocean: false }
                }
            })
            .collect();

        // Periodic steady state: exact affine fixed point of the lap map
        // for the clip pattern it takes; iterate until self-consistent.
        let mut g = 0.0f64;
        let mut rounds = 0u32;
        loop {
            rounds += 1;
            if rounds > MAX_RING_ROUNDS {
                return Err(ClimateError::NoConvergence {
                    stage: "moisture ring fixed point",
                    row: i,
                    rounds,
                });
            }
            let (_, alpha, beta) = ring_lap(&cells, g);
            let g_new = if beta < 1.0 { alpha / (1.0 - beta) } else { g };
            let (w_end, _, _) = ring_lap(&cells, g_new);
            let tol = RING_FIXED_POINT_REL_TOL * g_new.abs().max(f64::MIN_POSITIVE);
            g = g_new;
            if (w_end - g_new).abs() <= tol {
                break;
            }
        }

        // Final pass: record E, P, orographic P and the mean column.
        let mut w = g;
        for (s, cell) in cells.iter().enumerate() {
            let c = order[s];
            let p_oro = (w - cell.w_sat).max(0.0);
            w = w.min(cell.w_sat);
            let e = (-cell.b * cell.dt).exp();
            let ws = cell.a / cell.b;
            let w_int = ws * cell.dt + (w - ws) * (1.0 - e) / cell.b;
            let p_bg = (f / tau_r) * w_int;
            let evap = if cell.ocean {
                (w_eq * cell.dt - w_int) / tau_e
            } else {
                recyc * p_bg
            };
            w = ws + (w - ws) * e;
            let per_yr = SECONDS_PER_YEAR / cell.dt;
            fields.precipitation_mm_per_yr[c] = (p_oro + p_bg) * per_yr;
            fields.precipitation_orographic_mm_per_yr[c] = p_oro * per_yr;
            fields.evaporation_mm_per_yr[c] = evap * per_yr;
            fields.column_water_kg_m2[c] = w_int / cell.dt;
        }
        let ring_residual = g - w;

        // Seasonal half range, wind, insolation, PET, aridity.
        let resp = a1[i] * (1.0 - alb[i]) / (bb + cc_t);
        for &c in &order {
            let m = if inputs.ocean[c] {
                damp
            } else if ocean_cells == 0 {
                1.0
            } else {
                let x = fields.upwind_coast_distance_m[c];
                damp + (1.0 - damp) * (1.0 - (-x / l_cont).exp())
            };
            let h = m * resp;
            fields.temperature_half_range_k[c] = h;
            fields.wind_u_m_s[c] = u;
            fields.insolation_toa_w_m2[c] = q[i];
            let pet = thornthwaite_pet_mm_per_yr(fields.temperature_mean_c[c], h);
            fields.pet_mm_per_yr[c] = pet;
            fields.aridity_index[c] = if pet > 0.0 {
                fields.precipitation_mm_per_yr[c] / pet
            } else {
                ARIDITY_SENTINEL_PET_ZERO
            };
        }
        rows.push(ZonalRowV0 {
            lat_deg: lat[i],
            insolation_annual_w_m2: q[i],
            insolation_first_harmonic_w_m2: a1[i],
            albedo: alb[i],
            sea_level_temperature_c: t_sl[i],
            wind_u_m_s: u,
            rainout_factor: f,
            ocean_cells,
            ring_rounds: rounds,
            ring_residual_kg_m2: ring_residual,
        });
    }

    for (name, v) in fields.named() {
        if let Some(c) = v.iter().position(|x| !x.is_finite()) {
            return Err(ClimateError::NonFinite { field: name, index: c });
        }
    }

    // --- closures (ADR 039 C4).
    let (mut e_sum, mut p_sum) = (0.0, 0.0);
    for c in 0..n {
        e_sum += fields.evaporation_mm_per_yr[c] * areas[c];
        p_sum += fields.precipitation_mm_per_yr[c] * areas[c];
    }
    let closures = vec![
        CLIMATE_CONSERVED_V0[0].closure("climate-coarse-v0/moisture-rings", e_sum - p_sum, p_sum),
        CLIMATE_CONSERVED_V0[1].closure("climate-coarse-v0/zonal-ebm", absorbed - emitted, absorbed),
        CLIMATE_CONSERVED_V0[2].closure("climate-coarse-v0/zonal-ebm", transport, transport_gross),
    ];
    for cl in &closures {
        if !cl.closes() {
            return Err(ClimateError::ClosureFailed(cl.clone()));
        }
    }

    // --- summary.
    let ocean = inputs.ocean;
    let land = |c: usize| !ocean[c];
    let ice_edge = |north: bool| -> Option<f64> {
        let mut best: Option<f64> = None;
        for i in 0..nlat {
            if alb[i] == alb_ice && (lat[i] > 0.0) == north {
                let l = lat[i].abs();
                best = Some(best.map_or(l, |b: f64| b.min(l)));
            }
        }
        best.map(|l| if north { l } else { -l })
    };
    let land_area: f64 = (0..n).filter(|&c| land(c)).map(|c| areas[c]).sum();
    let mut classes = [0.0f64; 6];
    for c in (0..n).filter(|&c| land(c)) {
        let ai = fields.aridity_index[c];
        let slot = if ai == ARIDITY_SENTINEL_PET_ZERO {
            5
        } else {
            UNEP_ARIDITY_BOUNDS.iter().take_while(|&&b| ai >= b).count()
        };
        classes[slot] += areas[c];
    }
    let frac = |x: f64| if land_area > 0.0 { x / land_area } else { 0.0 };
    let land_p: f64 = (0..n).filter(|&c| land(c)).map(|c| fields.precipitation_mm_per_yr[c] * areas[c]).sum();
    let land_oro: f64 = (0..n).filter(|&c| land(c)).map(|c| fields.precipitation_orographic_mm_per_yr[c] * areas[c]).sum();
    let summary = ClimateSummaryV0 {
        solar_constant_w_m2: s0,
        global_mean_insolation_w_m2: (0..nlat).map(|i| weight[i] * q[i]).sum(),
        ebm_global_mean_sea_level_temperature_c: t_bar,
        ebm_iterations,
        ice_edge_lat_deg_north: ice_edge(true),
        ice_edge_lat_deg_south: ice_edge(false),
        global_mean_temperature_c: masked_mean(&fields.temperature_mean_c, &areas, |_| true).unwrap_or(0.0),
        land_mean_temperature_c: masked_mean(&fields.temperature_mean_c, &areas, land),
        ocean_mean_temperature_c: masked_mean(&fields.temperature_mean_c, &areas, |c| ocean[c]),
        land_area_fraction: land_area / total_area,
        global_mean_precipitation_mm_per_yr: p_sum / total_area,
        land_mean_precipitation_mm_per_yr: masked_mean(&fields.precipitation_mm_per_yr, &areas, land),
        ocean_mean_precipitation_mm_per_yr: masked_mean(&fields.precipitation_mm_per_yr, &areas, |c| ocean[c]),
        max_precipitation_mm_per_yr: fields.precipitation_mm_per_yr.iter().cloned().fold(0.0, f64::max),
        land_orographic_share: (land_p > 0.0).then(|| land_oro / land_p),
        rows_without_ocean,
        land_aridity_classes: AridityClassesV0 {
            hyper_arid: frac(classes[0]),
            arid: frac(classes[1]),
            semi_arid: frac(classes[2]),
            dry_subhumid: frac(classes[3]),
            humid: frac(classes[4]),
            pet_zero: frac(classes[5]),
        },
    };
    Ok(ClimateOutputV0 { fields, rows, summary, closures, seed_v1: None, seed_v2: None })
}

// ---------------------------------------------------------------------------
// The published domain payload (`pf.climate/v0`).

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClimateTierV0 {
    /// Tier id: `climate-coarse-v0` here; SEED-CLIMATE publishes its own.
    pub id: String,
    pub version: u32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ApproximationRecordV0 {
    pub id: String,
    pub statement: String,
    pub constants: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SentinelV0 {
    pub value: f64,
    pub meaning: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClimateFieldRecordV0 {
    pub name: String,
    pub units: String,
    pub meaning: String,
    pub status: FieldStatusV0,
    pub approximation: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sentinel: Option<SentinelV0>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConservedRecordV0 {
    pub id: String,
    pub unit: String,
    pub ledger: String,
    pub rel_tol: f64,
    pub abs_tol: f64,
    pub gross: String,
    pub note: String,
}

/// Which upstream fields the tier read (by name; the artifact's `inputs`
/// carry the digests of the artifacts that published them).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClimateInputRecordV0 {
    pub elevation_field: String,
    pub fill_level_field: String,
    pub sea_level_m: f64,
    pub ocean_source: String,
}

/// Domain payload of a `pf.climate/v0` artifact. Every tier fills the same
/// structure; `tier`, `approximations` and field statuses say what it did.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClimateDomainV0 {
    /// Optional diagnostics; absent preserves every v1 serialized byte.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seasonal_extension: Option<seed::seasons::Record>,
    pub tier: ClimateTierV0,
    pub recipe: ClimateRecipeV1,
    pub recipe_digest: String,
    pub inputs: ClimateInputRecordV0,
    pub fields: Vec<ClimateFieldRecordV0>,
    pub approximations: Vec<ApproximationRecordV0>,
    pub conserved: Vec<ConservedRecordV0>,
    pub closures: Vec<ClosureCheck>,
    pub rows: Vec<ZonalRowV0>,
    pub summary: ClimateSummaryV0,
    pub usage_limits: Vec<String>,
    /// SEED-CLIMATE: present only for `seed-climate-v1` (not serialized
    /// when absent, so climate-coarse-v0 payloads keep their bytes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_v1: Option<seed::SeedClimateRecordV1>,
    /// ROTATION-CLIMATE: present only for `seed-climate-v2` (not serialized
    /// when absent, so v0 and v1 payloads keep their bytes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_v2: Option<rotation::RotationClimateRecordV2>,
}

/// Usage limits a consumer must respect (ADR 032 usage rule).
pub const CLIMATE_USAGE_LIMITS_V0: &[&str] = &[
    "coarse tier: annual means and one seasonal harmonic; no weather, no interannual variability",
    "no meridional moisture or heat transport beyond the zonal relaxation term: no ITCZ convergence, no monsoons, no subtropical ocean export",
    "every constant is a project assumption awaiting review (ADR 039 C5); do not report against observed norms as a qualification",
    "wind_v_m_s is declared zero; wind_u_m_s is a band idealization, not a flow solution",
];

impl ClimateDomainV0 {
    pub fn build(
        out: &ClimateOutputV0,
        recipe: &ClimateRecipeV1,
        inputs: ClimateInputRecordV0,
    ) -> Result<Self, pf_contracts::ContractError> {
        // SEED-CLIMATE (lane CM): the records follow the recipe's law.
        let seed_law = recipe.law == ClimateLawV1::SeedClimateV1;
        let v2_law = recipe.law == ClimateLawV1::SeedClimateV2;
        let (field_specs, approximations, conserved, limits) = if v2_law {
            (
                rotation::SEED_CLIMATE_FIELDS_V2,
                rotation::SEED_CLIMATE_APPROXIMATIONS_V2,
                seed::SEED_CLIMATE_CONSERVED_V1,
                rotation::SEED_CLIMATE_USAGE_LIMITS_V2,
            )
        } else if seed_law {
            (
                seed::SEED_CLIMATE_FIELDS_V1,
                seed::SEED_CLIMATE_APPROXIMATIONS_V1,
                seed::SEED_CLIMATE_CONSERVED_V1,
                seed::SEED_CLIMATE_USAGE_LIMITS_V1,
            )
        } else {
            (CLIMATE_FIELDS_V0, CLIMATE_APPROXIMATIONS_V0, CLIMATE_CONSERVED_V0, CLIMATE_USAGE_LIMITS_V0)
        };
        let (tier_id, tier_version) = climate_tier_of(recipe.law);
        // CLIMATE-ARIDITY (2026-10-21): a run with the opt-in crh-rainout-v1
        // seed constants records that approximation and limit too (runs
        // without it are unchanged).
        let crh = (seed_law || v2_law) && out.seed_v1.as_ref().is_some_and(|r| r.constants.crh_rainout.is_some());
        let extra_approx: &[ApproximationV0] = if crh { std::slice::from_ref(&aridity::CRH_RAINOUT_APPROXIMATION_V1) } else { &[] };
        let extra_limit: &[&str] = if crh { &[aridity::CRH_RAINOUT_USAGE_LIMIT_V1] } else { &[] };
        Ok(Self {
            tier: ClimateTierV0 { id: tier_id.into(), version: tier_version },
            recipe: recipe.clone(),
            recipe_digest: recipe.digest()?,
            inputs,
            fields: field_specs
                .iter()
                .map(|f| ClimateFieldRecordV0 {
                    name: f.name.into(),
                    units: f.units.into(),
                    meaning: f.meaning.into(),
                    status: f.status,
                    approximation: f.approximation.into(),
                    sentinel: f.sentinel.map(|(value, meaning)| SentinelV0 { value, meaning: meaning.into() }),
                })
                .collect(),
            approximations: approximations
                .iter()
                .chain(extra_approx.iter())
                .map(|a| ApproximationRecordV0 {
                    id: a.id.into(),
                    statement: a.statement.into(),
                    constants: a.constants.iter().map(|s| s.to_string()).collect(),
                })
                .collect(),
            conserved: conserved
                .iter()
                .map(|q| ConservedRecordV0 {
                    id: q.id.into(),
                    unit: q.unit.into(),
                    ledger: q.ledger.into(),
                    rel_tol: q.rel_tol,
                    abs_tol: q.abs_tol,
                    gross: q.gross.into(),
                    note: q.note.into(),
                })
                .collect(),
            closures: out.closures.clone(),
            rows: out.rows.clone(),
            summary: out.summary.clone(),
            usage_limits: limits.iter().chain(extra_limit.iter()).map(|s| s.to_string()).collect(),
            seasonal_extension: None,
            seed_v1: if seed_law || v2_law { out.seed_v1.clone() } else { None },
            seed_v2: if v2_law { out.seed_v2.clone() } else { None },
        })
    }

    /// Structural checks a consumer can run on a loaded payload: the field
    /// list must be exactly the v0 schema names and units, in order.
    pub fn validate_schema(&self) -> Result<(), String> {
        if self.fields.len() != CLIMATE_FIELDS_V0.len() {
            return Err(format!(
                "{} fields declared, schema has {}",
                self.fields.len(),
                CLIMATE_FIELDS_V0.len()
            ));
        }
        for (f, s) in self.fields.iter().zip(CLIMATE_FIELDS_V0) {
            if f.name != s.name || f.units != s.units {
                return Err(format!(
                    "field {} [{}] does not match schema {} [{}]",
                    f.name, f.units, s.name, s.units
                ));
            }
        }
        let ids: Vec<&str> = self.approximations.iter().map(|a| a.id.as_str()).collect();
        for f in &self.fields {
            for part in f.approximation.split(" + ") {
                if !ids.contains(&part) {
                    return Err(format!("field {} names unknown approximation {part}", f.name));
                }
            }
        }
        Ok(())
    }
}
