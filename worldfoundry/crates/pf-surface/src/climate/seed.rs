//! SEED-CLIMATE (lane CM, 2026-10-12): the `seed-climate-v1` tier.
//!
//! Roadmap node SEED-CLIMATE: "Preview climate: latitude bands, winds, rain
//! shadow, lapse rate, biome classes". Policy (design 2026-10-12, LANES.md
//! "Efficiency policy"): coarse, loosely correct mechanics, "90% not 99.9%".
//!
//! WHAT CHANGES AGAINST `climate-coarse-v0` (which stays untouched and
//! byte-identical; this tier runs only when the recipe names
//! `seed-climate-v1`):
//! - moisture moves on the 2-D sphere, not per latitude ring: zonal band
//!   wind plus a declared meridional band wind (Hadley trades equatorward,
//!   Ferrel flow poleward, polar flow equatorward), so the column converges
//!   at the equator (ITCZ) and near 60 deg and diverges at 30 deg; plus an
//!   eddy diffusivity that carries water across band edges and into rings
//!   with no ocean (v0's zero-rain land rings, hard stripes);
//! - saturation excess no longer falls in the first wall cell: column water
//!   above saturation at the cell's lapse-cooled temperature is removed at
//!   the rate excess / tau_f with tau_f = fallout length / effective
//!   transport speed, so the windward maximum is spread over a declared
//!   distance and converges with grid size instead of scaling as 1/cell;
//!   the depleted column descending in the lee is the rain shadow;
//! - land evapotranspiration follows Budyko's curve E/P = B(PET/P) (capped)
//!   instead of a fixed recycling ratio, so humid interiors recycle more
//!   (v0's dry bias) and the column still loses water over land;
//! - continentality (seasonal amplitude) reads an advected and diffused
//!   ocean-influence tracer instead of a zonal fetch, so rings without ocean
//!   get a finite, smooth value.
//! The zonal energy balance, lapse rate, insolation, Thornthwaite PET and
//! the aridity index are the v0 approximations, reused unchanged.
//!
//! The steady state is solved by symmetric line Gauss-Seidel: each latitude
//! row is solved exactly (cyclic tridiagonal: zonal upwind advection and
//! diffusion), rows are coupled through the latest neighbour values, and the
//! nonlinear terms (saturation fallout, ocean evaporation that stops at the
//! equilibrium column, the Budyko ratio) are lagged by one sweep. Flux-form
//! upwind advection is conservative, so the global water ledger closes to
//! the solver tolerance (declared below, checked on every call).
//!
//! DETERMINISM: fixed sweep order, no hash maps, no clock, no randomness.

use super::{
    annual_insolation_w_m2, masked_mean, rainout_factor, saturation_column_kg_m2,
    thornthwaite_pet_mm_per_yr, zonal_wind_m_s, ApproximationV0, AridityClassesV0,
    ClimateError, ClimateFieldSpecV0, ClimateFieldsV0, ClimateInputs, ClimateOutputV0,
    ClimateSummaryV0, FieldStatusV0, ZonalRowV0, ARIDITY_SENTINEL_PET_ZERO,
    DISTANCE_SENTINEL_NO_UPWIND_OCEAN, SECONDS_PER_YEAR, UNEP_ARIDITY_BOUNDS,
};
use pf_contracts::laws::{
    ClimateLawV1, ClimateRecipeV1, ConstantSpecV1, ConstantV1, SourceClassV1, UnitV1,
    UniverseError, CLIMATE_RECIPE_SCHEMA_V1,
};
use pf_contracts::sim_module::ConservedQuantity;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

/// CLIMATE-GENESIS (lane CG, 2026-10-20): climate forcing of T0 genesis
/// erosion (`epoch-climate-runoff-v1`). Declared here rather than in
/// `climate.rs` to keep the lane's write set (path: `climate/genesis.rs`).
#[path = "genesis.rs"]
pub mod genesis;

/// Opt-in seasonal and subcell climate diagnostics.
#[path = "seasons.rs"]
pub mod seasons;

/// Tier id (= the law id it implements) and version.
pub const SEED_CLIMATE_TIER_ID: &str = "seed-climate-v1";
pub const SEED_CLIMATE_TIER_VERSION: u32 = 1;
/// Schema of the extra constants record this tier reads.
pub const SEED_CLIMATE_CONSTANTS_SCHEMA_V1: &str = "pf.seed-climate-constants/v1";
/// CLIMATE-ARIDITY (2026-10-21): schema of a constants set that carries the
/// `crh-rainout-v1` option (`SeedClimateConstantsV1::earth_crh_v1`). A set
/// without the option keeps the v1 schema, serialization and digest.
pub const SEED_CLIMATE_CONSTANTS_SCHEMA_CRH1: &str = "pf.seed-climate-constants/v1+crh-rainout-v1";

/// Solver numerics (fixed by the law version).
pub const SEED_MAX_SWEEPS: u32 = 40_000;
/// Stop when the largest change of column water in one sweep is below this
/// fraction of the largest column.
pub const SEED_SWEEP_REL_TOL: f64 = 1.0e-11;
/// Declared closure tolerance of the global water ledger (iterative solve).
pub const SEED_WATER_REL_TOL: f64 = 1.0e-6;
/// Continentality tracer below this counts as "no ocean influence".
pub const SEED_TRACER_FLOOR: f64 = 1.0e-30;

/// The constants `seed-climate-v1` reads beyond `ClimateConstantsV1`.
/// Until the integrator folds them into pf-contracts they live here with a
/// code-owned envelope and are recorded (with their digest) in the artifact.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SeedClimateConstantsV1 {
    pub schema: Cow<'static, str>,
    /// Meridional surface wind amplitude of the Hadley cell (equatorward).
    pub hadley_meridional_speed_m_s: ConstantV1,
    /// Meridional surface wind amplitude of the Ferrel cell (poleward).
    pub ferrel_meridional_speed_m_s: ConstantV1,
    /// Meridional surface wind amplitude of the polar cell (equatorward).
    pub polar_meridional_speed_m_s: ConstantV1,
    /// Eddy diffusivity of column water (transient eddies, band-edge mixing).
    pub eddy_moisture_diffusivity_m2_per_yr: ConstantV1,
    /// Distance over which saturation excess falls out, at the effective
    /// transport speed `moisture_transport_speed_m_s` of the base recipe.
    pub orographic_fallout_length_m: ConstantV1,
    /// Cap on the Budyko land evapotranspiration ratio E/P (below 1 so a
    /// land column always loses water).
    pub land_et_ratio_cap: ConstantV1,
    /// CLIMATE-ARIDITY (2026-10-21): opt-in column-relative-humidity
    /// rainout efficiency (`super::aridity`). Absent (and not serialized)
    /// in `earth_v1()`, so the v1 bytes and digest are unchanged.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub crh_rainout: Option<super::aridity::CrhRainoutV1>,
}

/// Code-owned envelopes (inclusive).
pub const SEED_CLIMATE_CONSTANT_SPECS_V1: &[ConstantSpecV1] = &[
    ConstantSpecV1 { id: "hadley_meridional_speed_m_s", unit: UnitV1::MeterPerSecond, lo: 0.0, hi: 5.0, why: "surface meridional wind of a three-cell idealization; above a few m/s the column flux would exceed any observed Hadley moisture transport" },
    ConstantSpecV1 { id: "ferrel_meridional_speed_m_s", unit: UnitV1::MeterPerSecond, lo: 0.0, hi: 5.0, why: "as above, Ferrel cell" },
    ConstantSpecV1 { id: "polar_meridional_speed_m_s", unit: UnitV1::MeterPerSecond, lo: 0.0, hi: 5.0, why: "as above, polar cell" },
    ConstantSpecV1 { id: "eddy_moisture_diffusivity_m2_per_yr", unit: UnitV1::SquareMeterPerYear, lo: 3.0e11, hi: 3.0e14, why: "eddy diffusivity of order 1e4 to 1e7 m2/s; positive keeps the line solves strictly dominant" },
    ConstantSpecV1 { id: "orographic_fallout_length_m", unit: UnitV1::Meter, lo: 1.0e4, hi: 1.0e6, why: "hydrometeor drift and windward-slope scale: tens to hundreds of km" },
    ConstantSpecV1 { id: "land_et_ratio_cap", unit: UnitV1::Dimensionless, lo: 0.0, hi: 0.95, why: "below 1 so every land column loses water and the moisture operator stays nonsingular" },
];

fn awaiting(basis: &'static str) -> SourceClassV1 {
    SourceClassV1::ProjectAssumptionAwaitingReview { basis: Cow::Borrowed(basis) }
}

impl SeedClimateConstantsV1 {
    /// Earth-like v1 values: every one recalled or chosen, not re-checked
    /// (ADR 039 C5).
    pub fn earth_v1() -> Self {
        let k = |value, unit, basis| ConstantV1::new(value, unit, awaiting(basis));
        Self {
            schema: Cow::Borrowed(SEED_CLIMATE_CONSTANTS_SCHEMA_V1),
            hadley_meridional_speed_m_s: k(1.5, UnitV1::MeterPerSecond, "project assumption: order of the zonal-mean surface meridional wind in the trades (~1-2 m/s, recalled, not checked against a reanalysis)"),
            ferrel_meridional_speed_m_s: k(0.5, UnitV1::MeterPerSecond, "project assumption: Ferrel-cell surface poleward flow is weaker than the Hadley return flow (recalled order)"),
            polar_meridional_speed_m_s: k(0.5, UnitV1::MeterPerSecond, "project assumption: weak polar-cell equatorward surface flow"),
            eddy_moisture_diffusivity_m2_per_yr: k(3.156e13, UnitV1::SquareMeterPerYear, "recalled: large-scale eddy diffusivity of order 1e6 m2/s used in moist diffusive energy-balance models (e.g. Frierson et al. 2007); not re-checked"),
            orographic_fallout_length_m: k(3.0e5, UnitV1::Meter, "project assumption: e-folding removal distance of saturation excess, lumping hydrometeor drift with a precipitation efficiency below one (the rest passes the crest as vapour); chosen by lane CM on seed 42 L1 so the steepest coastal wall peaks near 11,000 mm/yr (Earth record ~12,000, recalled) while lee shadows stay; 100 km gave 35,700 and 600 km 8,700 mm/yr"),
            land_et_ratio_cap: k(0.9, UnitV1::Dimensionless, "project assumption: cap of the Budyko E/P ratio so land columns always lose water"),
            crh_rainout: None,
        }
    }

    /// CLIMATE-ARIDITY (2026-10-21): `earth_v1()` plus the opt-in
    /// `crh-rainout-v1` option (`super::aridity`), under its own schema
    /// (so its digest differs from the v1 set's).
    pub fn earth_crh_v1() -> Self {
        let mut s = Self::earth_v1();
        s.schema = Cow::Borrowed(SEED_CLIMATE_CONSTANTS_SCHEMA_CRH1);
        s.crh_rainout = Some(super::aridity::CrhRainoutV1::earth_v1());
        s
    }

    fn entries(&self) -> [(&'static str, &ConstantV1); 6] {
        [
            ("hadley_meridional_speed_m_s", &self.hadley_meridional_speed_m_s),
            ("ferrel_meridional_speed_m_s", &self.ferrel_meridional_speed_m_s),
            ("polar_meridional_speed_m_s", &self.polar_meridional_speed_m_s),
            ("eddy_moisture_diffusivity_m2_per_yr", &self.eddy_moisture_diffusivity_m2_per_yr),
            ("orographic_fallout_length_m", &self.orographic_fallout_length_m),
            ("land_et_ratio_cap", &self.land_et_ratio_cap),
        ]
    }

    /// Fail-loud validation (schema, units, code-owned envelopes).
    pub fn validate(&self) -> Result<(), UniverseError> {
        // CLIMATE-ARIDITY: the schema names whether the crh option is
        // present (a set may not carry the option under the v1 schema).
        let expected = if self.crh_rainout.is_some() { SEED_CLIMATE_CONSTANTS_SCHEMA_CRH1 } else { SEED_CLIMATE_CONSTANTS_SCHEMA_V1 };
        if self.schema != expected {
            return Err(UniverseError::Schema {
                what: "seed climate constants",
                found: self.schema.to_string(),
                expected,
            });
        }
        if let Some(c) = &self.crh_rainout {
            c.validate()?;
        }
        for (id, c) in self.entries() {
            let spec = SEED_CLIMATE_CONSTANT_SPECS_V1
                .iter()
                .find(|s| s.id == id)
                .unwrap_or_else(|| panic!("no seed climate spec for {id}"));
            if c.unit != spec.unit {
                return Err(UniverseError::UnitMismatch {
                    constant: spec.id,
                    expected: spec.unit.as_str(),
                    found: c.unit.as_str(),
                });
            }
            if !c.value.is_finite() || !(spec.lo..=spec.hi).contains(&c.value) {
                return Err(UniverseError::OutOfEnvelope {
                    constant: spec.id,
                    value: c.value,
                    unit: spec.unit.as_str(),
                    lo: spec.lo,
                    hi: spec.hi,
                    why: spec.why,
                });
            }
        }
        Ok(())
    }

    /// blake3 over the canonical JSON (source labels included).
    pub fn digest(&self) -> Result<String, pf_contracts::ContractError> {
        pf_contracts::artifact::digest_canonical(self)
    }
}

/// The `seed-climate-v1` recipe: the v0 Earth constants under the new law
/// id, except the ice-free albedo (0.31, below); the extra constants are
/// `SeedClimateConstantsV1`.
pub fn seed_climate_v1_earth_recipe() -> ClimateRecipeV1 {
    let mut r = ClimateRecipeV1::coarse_v0_earth();
    r.schema = Cow::Borrowed(CLIMATE_RECIPE_SCHEMA_V1);
    r.law = ClimateLawV1::SeedClimateV1;
    // v0's ~3 K cool bias (global mean 11-12 degC, equator 24 degC) dries
    // every column through Clausius-Clapeyron. Lane CM measured on seed 42
    // L1: ice-free albedo 0.32 -> global mean 11.1 degC (ice edge 64.7),
    // 0.315 -> 12.5, 0.31 -> 14.0 (ice edge 70.3), 0.30 -> 17.3 (80.9).
    r.constants.albedo_ice_free = ConstantV1::new(
        0.31,
        UnitV1::Dimensionless,
        awaiting("project assumption: ice-free planetary albedo 0.31, between v0's 0.32 and the recalled Earth Bond albedo ~0.29-0.30; chosen by lane CM so the zonal EBM gives a global mean near 14 degC (recalled Earth ~14-15) with the ice edge near 70 deg; not re-checked"),
    );
    r
}

/// Recorded in the artifact payload of a `seed-climate-v1` run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SeedClimateRecordV1 {
    pub constants: SeedClimateConstantsV1,
    pub constants_digest: String,
    pub sweeps: u32,
    pub sweep_max_change_kg_m2: f64,
}

/// Every approximation `seed-climate-v1` makes. The first five and the last
/// three are the v0 approximations, reused unchanged.
pub const SEED_CLIMATE_APPROXIMATIONS_V1: &[ApproximationV0] = &[
    super::CLIMATE_APPROXIMATIONS_V0[0], // circular-orbit-annual-mean-insolation-v0
    super::CLIMATE_APPROXIMATIONS_V0[1], // budyko-sellers-zonal-ebm-v0
    super::CLIMATE_APPROXIMATIONS_V0[2], // constant-environmental-lapse-rate-v0
    super::CLIMATE_APPROXIMATIONS_V0[3], // annual-mean-no-maritime-correction-v0
    super::CLIMATE_APPROXIMATIONS_V0[4], // three-cell-zonal-wind-bands-v0
    ApproximationV0 {
        id: "three-cell-meridional-wind-bands-v1",
        statement: "Prevailing surface meridional wind v = -sign(lat) V_band sin(6 |lat|): equatorward in the Hadley (0-30 deg) and polar (60-90 deg) cells, poleward in the Ferrel cell (30-60 deg); zero at 0, 30, 60, 90 deg, so the column converges at the equator and near 60 deg and diverges at 30 deg and the poles. No seasonal migration, no monsoon, no topographic steering.",
        constants: &["hadley_meridional_speed_m_s", "ferrel_meridional_speed_m_s", "polar_meridional_speed_m_s"],
    },
    ApproximationV0 {
        id: "two-dimensional-moisture-transport-v1",
        statement: "Steady column water on the sphere: flux-form upwind advection by the band wind (u, v) plus eddy diffusion, ocean evaporation relaxing toward rh * saturation at the zonal sea-level temperature (stopping at that equilibrium), background rainout column / residence time scaled by the three-cell factor; solved by symmetric line Gauss-Seidel to a declared tolerance; the band wind is used as the column moisture-flux velocity.",
        constants: &["eddy_moisture_diffusivity_m2_per_yr", "rainout_residence_s", "ocean_evaporation_relaxation_s", "ocean_column_relative_humidity", "column_water_saturation_ref_kg_m2", "column_water_reference_temperature_c", "clausius_clapeyron_rate_per_k", "trade_wind_speed_m_s", "westerly_wind_speed_m_s", "polar_easterly_speed_m_s"],
    },
    super::CLIMATE_APPROXIMATIONS_V0[7], // three-cell-rainout-modulation-v0
    ApproximationV0 {
        id: "saturation-excess-fallout-length-v1",
        statement: "Column water above saturation at the cell's lapse-cooled surface temperature falls out at the rate excess / tau_f, tau_f = fallout length / effective transport speed, so the excess is removed over a declared distance downwind instead of in the first cell; the depleted column continues downwind (rain shadow in the lee). The orographic part is the excess above saturation at the cell's surface temperature minus the excess above saturation at the zonal sea-level temperature (the share caused by elevation). No flow blocking, no precipitation efficiency below one.",
        constants: &["orographic_fallout_length_m", "moisture_transport_speed_m_s", "lapse_rate_k_per_m", "clausius_clapeyron_rate_per_k", "column_water_saturation_ref_kg_m2", "column_water_reference_temperature_c"],
    },
    ApproximationV0 {
        id: "budyko-curve-land-evapotranspiration-v1",
        statement: "Land evapotranspiration E = P * min(cap, sqrt(phi tanh(1/phi) (1 - exp(-phi)))), phi = PET / P (Budyko 1974 curve, parameter free), returned to the column; no soil-moisture store, no seasonality of E.",
        constants: &["land_et_ratio_cap"],
    },
    ApproximationV0 {
        id: "advected-ocean-influence-continentality-v1",
        statement: "Ocean influence c is 1 on marine cells and, on land, is advected by the band wind (advective form, upwind), eddy-diffused and decays at rate U_eff / L_cont; seasonal half range = m * A1 (1 - albedo) / (B + C) with m = d + (1 - d)(1 - c); the published upwind coast distance is the effective distance -L_cont ln c (sentinel -1 where c is below 1e-30).",
        constants: &["ocean_seasonal_damping", "continentality_length_m", "moisture_transport_speed_m_s", "eddy_moisture_diffusivity_m2_per_yr", "olr_slope_w_m2_k", "meridional_transport_w_m2_k"],
    },
    super::CLIMATE_APPROXIMATIONS_V0[11], // thornthwaite-1948-unadjusted-pet-v0
    super::CLIMATE_APPROXIMATIONS_V0[12], // unep-aridity-index-v0
];

/// The same schema names and units as v0; statuses and approximations of
/// this tier (wind_v_m_s flips from declared-zero to modelled).
pub const SEED_CLIMATE_FIELDS_V1: &[ClimateFieldSpecV0] = &[
    ClimateFieldSpecV0 { name: "temperature_mean_c", units: "degC", meaning: "annual-mean near-surface air temperature", status: FieldStatusV0::Modelled, approximation: "budyko-sellers-zonal-ebm-v0 + constant-environmental-lapse-rate-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "temperature_half_range_k", units: "K", meaning: "seasonal half range (warmest-month mean = mean + half range, coldest = mean - half range)", status: FieldStatusV0::Modelled, approximation: "advected-ocean-influence-continentality-v1", sentinel: None },
    ClimateFieldSpecV0 { name: "precipitation_mm_per_yr", units: "mm/yr", meaning: "annual precipitation (kg/m2/yr of water)", status: FieldStatusV0::Modelled, approximation: "two-dimensional-moisture-transport-v1 + saturation-excess-fallout-length-v1 + three-cell-rainout-modulation-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "precipitation_orographic_mm_per_yr", units: "mm/yr", meaning: "part of precipitation_mm_per_yr from saturation rainout (lift over relief)", status: FieldStatusV0::Modelled, approximation: "saturation-excess-fallout-length-v1", sentinel: None },
    ClimateFieldSpecV0 { name: "evaporation_mm_per_yr", units: "mm/yr", meaning: "ocean evaporation and land evapotranspiration into the column", status: FieldStatusV0::Modelled, approximation: "two-dimensional-moisture-transport-v1 + budyko-curve-land-evapotranspiration-v1", sentinel: None },
    ClimateFieldSpecV0 { name: "pet_mm_per_yr", units: "mm/yr", meaning: "potential evapotranspiration", status: FieldStatusV0::Modelled, approximation: "thornthwaite-1948-unadjusted-pet-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "aridity_index", units: "1", meaning: "P / PET", status: FieldStatusV0::Modelled, approximation: "unep-aridity-index-v0", sentinel: Some((ARIDITY_SENTINEL_PET_ZERO, "PET is zero (no month above 0 degC): aridity undefined")) },
    ClimateFieldSpecV0 { name: "wind_u_m_s", units: "m/s", meaning: "prevailing surface wind, eastward component", status: FieldStatusV0::Modelled, approximation: "three-cell-zonal-wind-bands-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "wind_v_m_s", units: "m/s", meaning: "prevailing surface wind, northward component", status: FieldStatusV0::Modelled, approximation: "three-cell-meridional-wind-bands-v1", sentinel: None },
    ClimateFieldSpecV0 { name: "upwind_coast_distance_m", units: "m", meaning: "distance to the nearest ocean upwind along the moisture transport path", status: FieldStatusV0::Modelled, approximation: "advected-ocean-influence-continentality-v1", sentinel: Some((DISTANCE_SENTINEL_NO_UPWIND_OCEAN, "no ocean influence reaches the cell")) },
    ClimateFieldSpecV0 { name: "insolation_toa_w_m2", units: "W/m2", meaning: "annual-mean top-of-atmosphere insolation", status: FieldStatusV0::Modelled, approximation: "circular-orbit-annual-mean-insolation-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "column_water_kg_m2", units: "kg/m2", meaning: "annual-mean column water vapour along the transport path", status: FieldStatusV0::Modelled, approximation: "two-dimensional-moisture-transport-v1", sentinel: None },
];

/// Conserved quantities of this tier.
pub const SEED_CLIMATE_CONSERVED_V1: &[ConservedQuantity] = &[
    ConservedQuantity {
        id: "atmospheric-water",
        unit: "kg/yr",
        ledger: "sum over cells of (evaporation - precipitation) * cell area; zero at the steady state (flux-form transport telescopes)",
        rel_tol: SEED_WATER_REL_TOL,
        abs_tol: 1.0,
        gross: "sum over cells of precipitation * cell area",
        applies: None,
        note: "global closure; the residual is the iterative solver's, bounded by the declared sweep tolerance",
    },
    super::CLIMATE_CONSERVED_V0[1],
    super::CLIMATE_CONSERVED_V0[2],
];

pub const SEED_CLIMATE_USAGE_LIMITS_V1: &[&str] = &[
    "coarse tier: annual means and one seasonal harmonic; no weather, no interannual variability, no monsoon or seasonal ITCZ shift",
    "winds are three-cell band idealizations used as the moisture-flux velocity, not a flow solution; no flow blocking around high ranges",
    "every constant is a project assumption awaiting review (ADR 039 C5); do not report against observed norms as a qualification",
    "rows[].ring_rounds is the global sweep count and rows[].ring_residual_kg_m2 the row's largest column change in the last sweep (no per-ring budget in this tier)",
];

/// Meridional surface wind (m/s, northward positive) at `lat_deg`.
pub fn meridional_wind_m_s(lat_deg: f64, hadley: f64, ferrel: f64, polar: f64) -> f64 {
    let a = lat_deg.abs();
    if a == 0.0 {
        return 0.0;
    }
    let amp = if a < 30.0 {
        hadley
    } else if a < 60.0 {
        ferrel
    } else {
        polar
    };
    -lat_deg.signum() * amp * (6.0 * a.to_radians()).sin()
}

/// Budyko (1974) evapotranspiration ratio E/P for dryness index
/// phi = PET / P (no cap).
pub fn budyko_ratio(phi: f64) -> f64 {
    if phi <= 0.0 {
        return 0.0;
    }
    if !phi.is_finite() {
        return 1.0;
    }
    (phi * (1.0 / phi).tanh() * (1.0 - (-phi).exp())).sqrt()
}

/// Solve the cyclic tridiagonal system
/// diag[j] x[j] - west[j] x[j-1] - east[j] x[j+1] = rhs[j] (periodic j).
/// The rows here are strictly diagonally dominant, so no pivoting is needed.
fn solve_cyclic(diag: &[f64], west: &[f64], east: &[f64], rhs: &[f64], x: &mut [f64], scratch: &mut [Vec<f64>; 3]) {
    let n = diag.len();
    // Sub/super diagonals: a[j] multiplies x[j-1], c[j] multiplies x[j+1].
    // Corners: A[0][n-1] = -west[0] (beta), A[n-1][0] = -east[n-1] (alpha).
    let alpha = -east[n - 1];
    let beta = -west[0];
    let gamma = -diag[0];
    let [bb, z, cp] = scratch;
    bb.clear();
    bb.extend_from_slice(diag);
    bb[0] = diag[0] - gamma;
    bb[n - 1] = diag[n - 1] - alpha * beta / gamma;
    // Thomas on (a = -west, b = bb, c = -east) for rhs and for u.
    let thomas = |r: &dyn Fn(usize) -> f64, out: &mut [f64], cp: &mut Vec<f64>| {
        cp.clear();
        cp.resize(n, 0.0);
        let mut den = bb[0];
        cp[0] = -east[0] / den;
        out[0] = r(0) / den;
        for j in 1..n {
            den = bb[j] - (-west[j]) * cp[j - 1];
            cp[j] = if j + 1 < n { -east[j] / den } else { 0.0 };
            out[j] = (r(j) - (-west[j]) * out[j - 1]) / den;
        }
        for j in (0..n - 1).rev() {
            out[j] -= cp[j] * out[j + 1];
        }
    };
    thomas(&|j| rhs[j], x, cp);
    z.clear();
    z.resize(n, 0.0);
    thomas(
        &|j| {
            if j == 0 {
                gamma
            } else if j == n - 1 {
                alpha
            } else {
                0.0
            }
        },
        z,
        cp,
    );
    let fact = (x[0] + beta * x[n - 1] / gamma) / (1.0 + z[0] + beta * z[n - 1] / gamma);
    for j in 0..n {
        x[j] -= fact * z[j];
    }
}

/// Grid geometry for the finite-volume operators.
struct Geo {
    nlat: usize,
    nlon: usize,
    /// Cell area per row (m2).
    area: Vec<f64>,
    /// Zonal centre spacing per row (m).
    dx: Vec<f64>,
    /// Zonal face length (m) = meridional centre spacing.
    dy: f64,
    /// Meridional face length at edge k = 0..=nlat (m); 0 at the poles.
    lx: Vec<f64>,
    /// Latitude of edge k (deg).
    edge_lat: Vec<f64>,
}

/// Per-cell operator coefficients (already divided by the cell area, 1/s):
/// diag, west, east, south, north.
#[derive(Clone, Copy, Default)]
struct Coef {
    d: f64,
    w: f64,
    e: f64,
    s: f64,
    n: f64,
}

/// Transport operator rows. `flux_form`: conservative upwind (moisture);
/// otherwise advective upwind (tracer).
fn transport_coefs(g: &Geo, u_row: &[f64], v_edge: &[f64], k_diff: f64, flux_form: bool) -> Vec<Coef> {
    let mut out = vec![Coef::default(); g.nlat * g.nlon];
    for i in 0..g.nlat {
        let a = g.area[i];
        let u = u_row[i];
        let dzon = k_diff * g.dy / g.dx[i] / a;
        let mut row = Coef::default();
        // East face: outward normal velocity u; west face: -u.
        let (ue, uw) = (u * g.dy / a, -u * g.dy / a);
        // South face edge i: outward -v_edge[i]; north face edge i+1: +v.
        let (vs, vn) = (-v_edge[i] * g.lx[i] / a, v_edge[i + 1] * g.lx[i + 1] / a);
        let ds = k_diff * g.lx[i] / g.dy / a;
        let dn = k_diff * g.lx[i + 1] / g.dy / a;
        if flux_form {
            row.d = ue.max(0.0) + uw.max(0.0) + vs.max(0.0) + vn.max(0.0);
            row.e = (-ue).max(0.0);
            row.w = (-uw).max(0.0);
            row.s = (-vs).max(0.0);
            row.n = (-vn).max(0.0);
        } else {
            row.e = (-ue).max(0.0);
            row.w = (-uw).max(0.0);
            row.s = (-vs).max(0.0);
            row.n = (-vn).max(0.0);
            row.d = row.e + row.w + row.s + row.n;
        }
        row.d += 2.0 * dzon + ds + dn;
        row.e += dzon;
        row.w += dzon;
        row.s += ds;
        row.n += dn;
        for j in 0..g.nlon {
            out[i * g.nlon + j] = row;
        }
    }
    out
}

/// Symmetric line Gauss-Seidel on rows. `local(c, x_c)` returns the extra
/// (diag, rhs) of cell c given its current value. Returns (sweeps, last
/// max change, per-row last max change).
fn line_gauss_seidel(
    g: &Geo,
    coef: &[Coef],
    x: &mut [f64],
    mut local: impl FnMut(usize, f64) -> (f64, f64),
    stage: &'static str,
    sweep_rel_tol: f64,
) -> Result<(u32, f64, Vec<f64>), ClimateError> {
    let (nlat, nlon) = (g.nlat, g.nlon);
    let mut diag = vec![0.0; nlon];
    let mut west = vec![0.0; nlon];
    let mut east = vec![0.0; nlon];
    let mut rhs = vec![0.0; nlon];
    let mut sol = vec![0.0; nlon];
    let mut scratch = [Vec::new(), Vec::new(), Vec::new()];
    let mut row_change = vec![0.0; nlat];
    let mut sweeps = 0u32;
    loop {
        sweeps += 1;
        if sweeps > SEED_MAX_SWEEPS {
            return Err(ClimateError::NoConvergence { stage, row: 0, rounds: sweeps });
        }
        let mut max_change = 0.0f64;
        for step in 0..nlat {
            let i = if sweeps % 2 == 1 { step } else { nlat - 1 - step };
            for j in 0..nlon {
                let c = i * nlon + j;
                let k = coef[c];
                let (ld, lr) = local(c, x[c]);
                diag[j] = k.d + ld;
                west[j] = k.w;
                east[j] = k.e;
                let mut r = lr;
                if i > 0 {
                    r += k.s * x[c - nlon];
                }
                if i + 1 < nlat {
                    r += k.n * x[c + nlon];
                }
                rhs[j] = r;
            }
            solve_cyclic(&diag, &west, &east, &rhs, &mut sol, &mut scratch);
            let mut rc = 0.0f64;
            for j in 0..nlon {
                let c = i * nlon + j;
                let v = sol[j].max(0.0);
                rc = rc.max((v - x[c]).abs());
                x[c] = v;
            }
            row_change[i] = rc;
            max_change = max_change.max(rc);
        }
        let scale = x.iter().cloned().fold(0.0f64, f64::max).max(f64::MIN_POSITIVE);
        if max_change <= sweep_rel_tol * scale {
            return Ok((sweeps, max_change, row_change));
        }
    }
}

/// The `seed-climate-v1` tier. Pure: identical inputs give identical
/// outputs; the only spatial inputs are `inputs`.
pub fn derive_seed_climate_v1(
    inputs: &ClimateInputs<'_>,
    recipe: &ClimateRecipeV1,
    extra: &SeedClimateConstantsV1,
) -> Result<ClimateOutputV0, ClimateError> {
    if recipe.law != ClimateLawV1::SeedClimateV1 {
        return Err(ClimateError::Recipe(UniverseError::Unsupported {
            law: "climate/climate-coarse-v0",
            why: "derive_seed_climate_v1 implements seed-climate-v1 only",
        }));
    }
    // The base constants carry the v0 envelopes. `ClimateRecipeV1::validate`
    // refuses the seed law id until the integrator lands the request patch,
    // so validate the same constants under the v0 id (same checks).
    let mut base = recipe.clone();
    base.law = ClimateLawV1::ClimateCoarseV0;
    base.validate()?;
    extra.validate()?;
    let k = &recipe.constants;
    let (trade, west, polar) = (k.trade_wind_speed_m_s.value, k.westerly_wind_speed_m_s.value, k.polar_easterly_speed_m_s.value);
    let (vh, vf, vp) = (
        extra.hadley_meridional_speed_m_s.value,
        extra.ferrel_meridional_speed_m_s.value,
        extra.polar_meridional_speed_m_s.value,
    );
    let k_circ = k.circulation_rainout_modulation.value;
    seed_climate_core(
        inputs,
        recipe,
        extra,
        &|l| zonal_wind_m_s(l, trade, west, polar),
        &|l| meridional_wind_m_s(l, vh, vf, vp),
        &|l| rainout_factor(l, k_circ),
    )
}

/// CLIMATE-GENESIS (lane CG, 2026-10-20): `derive_seed_climate_v1` with
/// solve options, also returning the converged column state. With
/// `SeedSolveOptions::default()` the output equals `derive_seed_climate_v1`.
pub fn derive_seed_climate_v1_with(
    inputs: &ClimateInputs<'_>,
    recipe: &ClimateRecipeV1,
    extra: &SeedClimateConstantsV1,
    opts: &SeedSolveOptions<'_>,
) -> Result<(ClimateOutputV0, SeedColumnStateV1), ClimateError> {
    if recipe.law != ClimateLawV1::SeedClimateV1 {
        return Err(ClimateError::Recipe(UniverseError::Unsupported {
            law: "climate/climate-coarse-v0",
            why: "derive_seed_climate_v1 implements seed-climate-v1 only",
        }));
    }
    let mut base = recipe.clone();
    base.law = ClimateLawV1::ClimateCoarseV0;
    base.validate()?;
    extra.validate()?;
    let k = &recipe.constants;
    let (trade, west, polar) = (k.trade_wind_speed_m_s.value, k.westerly_wind_speed_m_s.value, k.polar_easterly_speed_m_s.value);
    let (vh, vf, vp) = (
        extra.hadley_meridional_speed_m_s.value,
        extra.ferrel_meridional_speed_m_s.value,
        extra.polar_meridional_speed_m_s.value,
    );
    let k_circ = k.circulation_rainout_modulation.value;
    seed_climate_core_with(
        inputs,
        recipe,
        extra,
        &|l| zonal_wind_m_s(l, trade, west, polar),
        &|l| meridional_wind_m_s(l, vh, vf, vp),
        &|l| rainout_factor(l, k_circ),
        opts,
    )
}

/// ROTATION-CLIMATE (lane RC, 2026-10-20): the `seed-climate-v1` solver
/// with the circulation geometry passed in. `u_of`, `v_of` and `f_of` give
/// the prevailing zonal wind, meridional wind and background rainout factor
/// at a latitude (deg). `seed-climate-v1` passes the fixed three-cell bands
/// (same arithmetic as before the split, so its bytes are unchanged);
/// `seed-climate-v2` (`climate/rotation.rs`) passes rotation-scaled bands.
/// The caller validates the recipe and the extra constants.
pub(super) fn seed_climate_core(
    inputs: &ClimateInputs<'_>,
    recipe: &ClimateRecipeV1,
    extra: &SeedClimateConstantsV1,
    u_of: &dyn Fn(f64) -> f64,
    v_of: &dyn Fn(f64) -> f64,
    f_of: &dyn Fn(f64) -> f64,
) -> Result<ClimateOutputV0, ClimateError> {
    seed_climate_core_with(inputs, recipe, extra, u_of, v_of, f_of, &SeedSolveOptions::default()).map(|(out, _)| out)
}

/// CLIMATE-GENESIS (lane CG, 2026-10-20): options of one `seed-climate`
/// solve. `SeedSolveOptions::default()` is the post-hoc tier exactly (cold
/// start from the equilibrium column and the ocean mask, tolerance
/// `SEED_SWEEP_REL_TOL`), so every existing artifact keeps its bytes.
/// Genesis epochs warm-start from the previous epoch's converged state and
/// declare their own (looser) tolerance; the water closure at
/// `SEED_WATER_REL_TOL` is still checked on every solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeedSolveOptions<'a> {
    /// Initial column water per cell (kg/m2); None = the v1 cold start.
    pub initial_column_water_kg_m2: Option<&'a [f64]>,
    /// Initial continentality tracer per cell; None = the v1 cold start.
    pub initial_tracer: Option<&'a [f64]>,
    /// Stop when the largest change of one sweep is below this fraction of
    /// the largest value (both line solves).
    pub sweep_rel_tol: f64,
}

impl Default for SeedSolveOptions<'_> {
    fn default() -> Self {
        Self { initial_column_water_kg_m2: None, initial_tracer: None, sweep_rel_tol: SEED_SWEEP_REL_TOL }
    }
}

/// CLIMATE-GENESIS: the converged solver state of one solve and the zonal
/// rows a local (pointwise) rainout needs. Rows are indexed by latitude
/// row of the solve grid; cells row-major on the solve grid.
#[derive(Debug, Clone, PartialEq)]
pub struct SeedColumnStateV1 {
    pub nlat: usize,
    pub nlon: usize,
    /// Converged column water (kg/m2).
    pub column_water_kg_m2: Vec<f64>,
    /// Converged continentality tracer (ocean influence, 0..1).
    pub tracer: Vec<f64>,
    /// Zonal sea-level temperature per row (degC).
    pub sea_level_temperature_c: Vec<f64>,
    /// Seasonal half-range response per row, A1 (1 - albedo) / (B + C) (K).
    pub half_range_response_k: Vec<f64>,
    /// Background rainout factor per row (circulation modulation).
    pub rainout_factor: Vec<f64>,
    pub tracer_sweeps: u32,
    pub moisture_sweeps: u32,
}

/// `seed_climate_core` with solve options; also returns the converged
/// column state (CLIMATE-GENESIS). Identical arithmetic: with
/// `SeedSolveOptions::default()` the output is the post-hoc tier's.
pub(super) fn seed_climate_core_with(
    inputs: &ClimateInputs<'_>,
    recipe: &ClimateRecipeV1,
    extra: &SeedClimateConstantsV1,
    u_of: &dyn Fn(f64) -> f64,
    v_of: &dyn Fn(f64) -> f64,
    f_of: &dyn Fn(f64) -> f64,
    opts: &SeedSolveOptions<'_>,
) -> Result<(ClimateOutputV0, SeedColumnStateV1), ClimateError> {
    if !(opts.sweep_rel_tol.is_finite() && opts.sweep_rel_tol > 0.0) {
        return Err(ClimateError::Inputs(format!("sweep tolerance {} must be finite and positive", opts.sweep_rel_tol)));
    }
    let grid = inputs.grid;
    let (nlat, nlon, n) = (grid.nlat(), grid.nlon(), grid.cells());
    if nlon < 3 || nlat < 2 {
        return Err(ClimateError::Inputs(format!("grid {nlat} x {nlon} too small for the line solver")));
    }
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
    let (aa, bb, cc_t) = (k.olr_intercept_w_m2.value, k.olr_slope_w_m2_k.value, k.meridional_transport_w_m2_k.value);
    let (alb_free, alb_ice, t_ice) = (k.albedo_ice_free.value, k.albedo_ice.value, k.ice_albedo_threshold_c.value);
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
    let damp = k.ocean_seasonal_damping.value;
    let l_cont = k.continentality_length_m.value;
    let k_diff = extra.eddy_moisture_diffusivity_m2_per_yr.value / SECONDS_PER_YEAR;
    let tau_f = extra.orographic_fallout_length_m.value / u_eff;
    let et_cap = extra.land_et_ratio_cap.value;
    // CLIMATE-ARIDITY: the opt-in crh rainout (None = v1 arithmetic).
    let crh = extra.crh_rainout.as_ref().map(|c| c.params());

    // --- geometry.
    let radius = grid.radius_m();
    let dlat = std::f64::consts::PI / nlat as f64;
    let dlon = 2.0 * std::f64::consts::PI / nlon as f64;
    let mut lat = Vec::with_capacity(nlat);
    let mut row_area = Vec::with_capacity(nlat);
    for i in 0..nlat {
        lat.push(grid.center_lat_deg(i).map_err(|e| ClimateError::Inputs(e.to_string()))?);
        row_area.push(grid.cell_area_m2(i).map_err(|e| ClimateError::Inputs(e.to_string()))?);
    }
    let edge_lat: Vec<f64> = (0..=nlat).map(|kk| -90.0 + kk as f64 * 180.0 / nlat as f64).collect();
    let lx: Vec<f64> = edge_lat
        .iter()
        .enumerate()
        .map(|(kk, &l)| if kk == 0 || kk == nlat { 0.0 } else { radius * l.to_radians().cos() * dlon })
        .collect();
    let geo = Geo {
        nlat,
        nlon,
        area: row_area.clone(),
        dx: lat.iter().map(|l| radius * l.to_radians().cos() * dlon).collect(),
        dy: radius * dlat,
        lx,
        edge_lat,
    };
    let total_area: f64 = row_area.iter().map(|a| a * nlon as f64).sum();
    let weight: Vec<f64> = row_area.iter().map(|a| a * nlon as f64 / total_area).collect();
    let areas: Vec<f64> = (0..n).map(|c| row_area[c / nlon]).collect();

    // --- insolation and zonal EBM (v0 approximations, same algorithm).
    let (mut q, mut a1) = (Vec::with_capacity(nlat), Vec::with_capacity(nlat));
    for &l in &lat {
        let (m, h) = annual_insolation_w_m2(s0, l, k.axial_tilt_deg.value);
        q.push(m);
        a1.push(h);
    }
    let mut alb = vec![alb_free; nlat];
    let mut t_sl = vec![0.0; nlat];
    let mut t_bar: f64;
    let mut ebm_iterations = 0u32;
    loop {
        ebm_iterations += 1;
        if ebm_iterations as usize > nlat + 2 {
            return Err(ClimateError::NoConvergence { stage: "zonal energy balance (ice-albedo pattern)", row: 0, rounds: ebm_iterations });
        }
        t_bar = (0..nlat).map(|i| weight[i] * ((1.0 - alb[i]) * q[i] - aa)).sum::<f64>() / bb;
        for i in 0..nlat {
            t_sl[i] = ((1.0 - alb[i]) * q[i] - aa + cc_t * t_bar) / (bb + cc_t);
        }
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
    let transport_gross: f64 = (0..nlat).map(|i| weight[i] * (cc_t * (t_bar - t_sl[i])).abs()).sum();

    // --- winds (circulation geometry from the caller).
    let u_row: Vec<f64> = lat.iter().map(|&l| u_of(l)).collect();
    let v_row: Vec<f64> = lat.iter().map(|&l| v_of(l)).collect();
    let v_edge: Vec<f64> = geo
        .edge_lat
        .iter()
        .enumerate()
        .map(|(kk, &l)| if kk == 0 || kk == nlat { 0.0 } else { v_of(l) })
        .collect();
    let f_row: Vec<f64> = lat.iter().map(|&l| f_of(l)).collect();

    // --- surface temperature (lapse rate; marine cells at the zonal T).
    let sea = inputs.sea_level_m;
    let ocean = inputs.ocean;
    let t_mean: Vec<f64> = (0..n)
        .map(|c| {
            let i = c / nlon;
            if ocean[c] {
                t_sl[i]
            } else {
                let surface = inputs.elevation_m[c].max(inputs.fill_level_m[c]);
                t_sl[i] - lapse * (surface - sea)
            }
        })
        .collect();

    // --- continentality tracer (advective form; ocean Dirichlet 1).
    let tracer_coef = transport_coefs(&geo, &u_row, &v_edge, k_diff, false);
    let decay = u_eff / l_cont;
    let mut tracer: Vec<f64> = match opts.initial_tracer {
        None => (0..n).map(|c| if ocean[c] { 1.0 } else { 0.0 }).collect(),
        Some(t) => warm_start(t, n, "initial tracer")?,
    };
    let big = 1.0e30;
    let (tracer_sweeps, _, _) = line_gauss_seidel(
        &geo,
        &tracer_coef,
        &mut tracer,
        |c, _| if ocean[c] { (big, big) } else { (decay, 0.0) },
        "continentality tracer",
        opts.sweep_rel_tol,
    )?;
    let mut half_range = vec![0.0; n];
    let mut coast_dist = vec![DISTANCE_SENTINEL_NO_UPWIND_OCEAN; n];
    let mut pet = vec![0.0; n];
    for c in 0..n {
        let i = c / nlon;
        let resp = a1[i] * (1.0 - alb[i]) / (bb + cc_t);
        let (m, x) = if ocean[c] {
            (damp, 0.0)
        } else {
            let t = tracer[c].clamp(0.0, 1.0);
            let x = if t > SEED_TRACER_FLOOR { -l_cont * t.ln() } else { DISTANCE_SENTINEL_NO_UPWIND_OCEAN };
            (damp + (1.0 - damp) * (1.0 - t), x)
        };
        half_range[c] = m * resp;
        coast_dist[c] = x;
        pet[c] = thornthwaite_pet_mm_per_yr(t_mean[c], half_range[c]) / SECONDS_PER_YEAR; // kg/m2/s
    }

    // --- moisture (flux form).
    let w_sat: Vec<f64> = (0..n).map(|c| saturation_column_kg_m2(w_ref, t_ref, cc, t_mean[c])).collect();
    let w_sat_sl: Vec<f64> = (0..n).map(|c| saturation_column_kg_m2(w_ref, t_ref, cc, t_sl[c / nlon])).collect();
    let w_eq_row: Vec<f64> = t_sl.iter().map(|&t| rh * saturation_column_kg_m2(w_ref, t_ref, cc, t)).collect();
    let moist_coef = transport_coefs(&geo, &u_row, &v_edge, k_diff, true);
    let precip_rate = |c: usize, w: f64| -> f64 {
        match crh {
            None => f_row[c / nlon] * w / tau_r + (w - w_sat[c]).max(0.0) / tau_f,
            Some(p) => p.rate(f_row[c / nlon] / tau_r, w, w_sat[c]) + (w - w_sat[c]).max(0.0) / tau_f,
        }
    };
    let et_ratio = |c: usize, w: f64| -> f64 {
        let p = precip_rate(c, w);
        if p > 0.0 {
            budyko_ratio(pet[c] / p).min(et_cap)
        } else {
            et_cap
        }
    };
    let mut wv: Vec<f64> = match opts.initial_column_water_kg_m2 {
        None => (0..n).map(|c| w_eq_row[c / nlon]).collect(),
        Some(w) => warm_start(w, n, "initial column water")?,
    };
    let (sweeps, last_change, row_change) = line_gauss_seidel(
        &geo,
        &moist_coef,
        &mut wv,
        |c, w| {
            let f = f_row[c / nlon];
            let sat_on = w > w_sat[c];
            if let Some(p) = crh {
                // CLIMATE-ARIDITY: Newton linearization of the steep
                // background sink s(w) = (f / tau_r) w G(w / w_sat):
                // s(w) ~ s'(w0) w - (s'(w0) w0 - s(w0)); the constant part
                // k G a r (1 - G) w0 is >= 0 (logistic G), so the row stays
                // an M-matrix with a non-negative right-hand side.
                let (s0, ds) = p.rate_and_slope(f / tau_r, w, w_sat[c]);
                let lin_r = ds * w - s0;
                if ocean[c] {
                    let w_eq = w_eq_row[c / nlon];
                    let (mut d, mut r) = (ds, lin_r);
                    if w < w_eq {
                        d += 1.0 / tau_e;
                        r += w_eq / tau_e;
                    }
                    if sat_on {
                        d += 1.0 / tau_f;
                        r += w_sat[c] / tau_f;
                    }
                    return (d, r);
                }
                let keep = 1.0 - et_ratio(c, w);
                let (mut d, mut r) = (keep * ds, keep * lin_r);
                if sat_on {
                    d += keep / tau_f;
                    r += keep * w_sat[c] / tau_f;
                }
                return (d, r);
            }
            if ocean[c] {
                let w_eq = w_eq_row[c / nlon];
                let (mut d, mut r) = (f / tau_r, 0.0);
                if w < w_eq {
                    d += 1.0 / tau_e;
                    r += w_eq / tau_e;
                }
                if sat_on {
                    d += 1.0 / tau_f;
                    r += w_sat[c] / tau_f;
                }
                (d, r)
            } else {
                let keep = 1.0 - et_ratio(c, w);
                let (mut d, mut r) = (keep * f / tau_r, 0.0);
                if sat_on {
                    d += keep / tau_f;
                    r += keep * w_sat[c] / tau_f;
                }
                (d, r)
            }
        },
        "two-dimensional moisture steady state",
        opts.sweep_rel_tol,
    )?;

    let mut fields = ClimateFieldsV0 {
        temperature_mean_c: t_mean,
        temperature_half_range_k: half_range,
        precipitation_mm_per_yr: vec![0.0; n],
        precipitation_orographic_mm_per_yr: vec![0.0; n],
        evaporation_mm_per_yr: vec![0.0; n],
        pet_mm_per_yr: vec![0.0; n],
        aridity_index: vec![0.0; n],
        wind_u_m_s: (0..n).map(|c| u_row[c / nlon]).collect(),
        wind_v_m_s: (0..n).map(|c| v_row[c / nlon]).collect(),
        upwind_coast_distance_m: coast_dist,
        insolation_toa_w_m2: (0..n).map(|c| q[c / nlon]).collect(),
        column_water_kg_m2: wv.clone(),
    };
    let (mut e_sum, mut p_sum) = (0.0, 0.0);
    for c in 0..n {
        let w = wv[c];
        let p = precip_rate(c, w);
        let p_oro = ((w - w_sat[c]).max(0.0) - (w - w_sat_sl[c]).max(0.0)).max(0.0) / tau_f;
        let e = if ocean[c] {
            (w_eq_row[c / nlon] - w).max(0.0) / tau_e
        } else {
            et_ratio(c, w) * p
        };
        e_sum += e * SECONDS_PER_YEAR * areas[c];
        p_sum += p * SECONDS_PER_YEAR * areas[c];
        fields.precipitation_mm_per_yr[c] = p * SECONDS_PER_YEAR;
        fields.precipitation_orographic_mm_per_yr[c] = p_oro * SECONDS_PER_YEAR;
        fields.evaporation_mm_per_yr[c] = e * SECONDS_PER_YEAR;
        let pet_yr = pet[c] * SECONDS_PER_YEAR;
        fields.pet_mm_per_yr[c] = pet_yr;
        fields.aridity_index[c] = if pet_yr > 0.0 { p * SECONDS_PER_YEAR / pet_yr } else { ARIDITY_SENTINEL_PET_ZERO };
    }
    for (name, v) in fields.named() {
        if let Some(c) = v.iter().position(|x| !x.is_finite()) {
            return Err(ClimateError::NonFinite { field: name, index: c });
        }
    }
    let closures = vec![
        SEED_CLIMATE_CONSERVED_V1[0].closure("seed-climate-v1/moisture-2d", e_sum - p_sum, p_sum),
        SEED_CLIMATE_CONSERVED_V1[1].closure("seed-climate-v1/zonal-ebm", absorbed - emitted, absorbed),
        SEED_CLIMATE_CONSERVED_V1[2].closure("seed-climate-v1/zonal-ebm", transport, transport_gross),
    ];
    for cl in &closures {
        if !cl.closes() {
            return Err(ClimateError::ClosureFailed(cl.clone()));
        }
    }

    // --- rows and summary.
    let mut rows = Vec::with_capacity(nlat);
    let mut rows_without_ocean = 0u32;
    for i in 0..nlat {
        let ocean_cells = (0..nlon).filter(|&j| ocean[i * nlon + j]).count() as u32;
        if ocean_cells == 0 {
            rows_without_ocean += 1;
        }
        rows.push(ZonalRowV0 {
            lat_deg: lat[i],
            insolation_annual_w_m2: q[i],
            insolation_first_harmonic_w_m2: a1[i],
            albedo: alb[i],
            sea_level_temperature_c: t_sl[i],
            wind_u_m_s: u_row[i],
            rainout_factor: f_row[i],
            ocean_cells,
            ring_rounds: sweeps,
            ring_residual_kg_m2: row_change[i],
        });
    }
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
        let slot = if ai == ARIDITY_SENTINEL_PET_ZERO { 5 } else { UNEP_ARIDITY_BOUNDS.iter().take_while(|&&b| ai >= b).count() };
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
    let record = SeedClimateRecordV1 {
        constants_digest: extra
            .digest()
            .map_err(|e| ClimateError::Inputs(format!("seed climate constants digest: {e}")))?,
        constants: extra.clone(),
        sweeps,
        sweep_max_change_kg_m2: last_change,
    };
    let state = SeedColumnStateV1 {
        nlat,
        nlon,
        column_water_kg_m2: wv,
        tracer,
        sea_level_temperature_c: t_sl,
        half_range_response_k: (0..nlat).map(|i| a1[i] * (1.0 - alb[i]) / (bb + cc_t)).collect(),
        rainout_factor: f_row,
        tracer_sweeps,
        moisture_sweeps: sweeps,
    };
    Ok((ClimateOutputV0 { fields, rows, summary, closures, seed_v1: Some(record), seed_v2: None }, state))
}

/// A warm-start vector: right length, finite, non-negative (fail loud).
fn warm_start(v: &[f64], n: usize, what: &'static str) -> Result<Vec<f64>, ClimateError> {
    if v.len() != n {
        return Err(ClimateError::Inputs(format!("{what}: {} values for {n} cells", v.len())));
    }
    if let Some(c) = v.iter().position(|x| !(x.is_finite() && *x >= 0.0)) {
        return Err(ClimateError::Inputs(format!("{what}: value {} at cell {c} is not finite and non-negative", v[c])));
    }
    Ok(v.to_vec())
}
