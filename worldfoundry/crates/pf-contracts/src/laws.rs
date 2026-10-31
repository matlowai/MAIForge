//! WS-CONST (lane WC, 2026-10-12): thin, versioned `LawSetV1` and
//! `UniverseConstantsV1` (ADR 039 §1-§3, design review 5.1 disposition C1-C5).
//!
//! WHAT THIS IS: the typed record of the laws a world runs under and the
//! physical constants the current generator path reads, each constant with
//! its value, unit, source class and a code-owned validity envelope. It is
//! deliberately thin (ADR 039 §1, design review: "WS-CONST stays thin: it versions
//! only the inputs the D1 calculation actually needs, with units, source
//! class and validity, plus legacy recipe preservation. No general law
//! interpreter is built first."):
//!
//! - A law id names a supported, versioned executable implementation. Law
//!   ids are serde enum variants, so an unknown id is refused at load; a
//!   string in data can never make an unknown law executable. `PrattV0` is
//!   NAMED (ADR 039 §1 "Airy or Pratt choice") and refused as `Unsupported`.
//! - Constants carry `value` + `unit` + `source`. The validity envelope is
//!   CODE-owned (`constant_spec`), never data: a recipe cannot widen it.
//! - Ocean water volume is the primary sea-level input under the two-column
//!   law (`SeaLevelLawV1::VolumePrimaryGlobalFillV0`); the sea-level constant
//!   is then only the output datum. A prescribed sea level and an ocean
//!   volume are never both inputs (ADR 039 §3, fresh-eyes F15).
//! - Recalled-but-unchecked citations are `ProjectAssumptionAwaitingReview`
//!   (design review C5), never `Cited`.
//!
//! LEGACY PRESERVATION: the v0 recipe encoding (the `PlanetParameterProfile`
//! blocks) is unchanged. `PlanetParameterProfile::universe()` resolves it to
//! a `UniverseV1`; with no physics flag set this is the `legacy` profile
//! (single-column Airy against the 35 km reference, prescribed 0 m sea level,
//! stated 9.81 m/s^2 gravity, whole-plate continents) and the generator reads
//! the same f64 values it read before, so outputs are byte-identical. A
//! recipe MAY pin its resolution (`PlanetParameterProfile.pinned_universe`,
//! skipped on write when absent so every existing params digest is
//! unchanged); a pinned recipe then refuses to run if this build resolves
//! different physics (law ids, law parameters, constant values or units).
//! Source labels are not part of that comparison: checking a citation later
//! must not invalidate pinned worlds.
//!
//! NOT BUILT HERE (by design): a law interpreter, the SimModule trait, the
//! LawContextQuery of ADR 041 §5 (in-world domains read an effective local
//! context; that consumer is future work and nothing here blocks it: it can
//! wrap `UniverseV1` as the base context), planet mass / rotation / star /
//! heat constants (no executable consumer yet), and overrides for code-fixed
//! constants (crust thicknesses, the reference column, the legacy water
//! density): in v1 they are lifted and recorded, not settable.

use crate::params::Compensation;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

/// Schema ids (versioned; a different string is refused loudly).
pub const UNIVERSE_SCHEMA_V1: &str = "pf.universe/v1";
pub const LAW_SET_SCHEMA_V1: &str = "pf.law-set/v1";
pub const UNIVERSE_CONSTANTS_SCHEMA_V1: &str = "pf.universe-constants/v1";

/// Values the legacy preview chain hard-coded outside the parameter
/// profile, LIFTED here so the legacy profile states them once. Each equals
/// the code constant named in its comment (pinned by
/// pf-lab tests/ws_const_preview.rs; pf-contracts cannot depend on the
/// crates that define the originals).
pub mod legacy {
    /// = pf_geodynamics::isostasy::REFERENCE_CRUST_M (Airy column at 0 m).
    pub const REFERENCE_CRUST_THICKNESS_M: f64 = 35_000.0;
    /// = pf_geodynamics::uplift::CONTINENTAL_CRUST_THICKNESS_M.
    pub const CONTINENTAL_CRUST_THICKNESS_M: f64 = 35_000.0;
    /// = pf_geodynamics::uplift::OCEANIC_CRUST_THICKNESS_M.
    pub const OCEANIC_CRUST_THICKNESS_M: f64 = 7_000.0;
    /// = pf_lab::preview::WATER_DENSITY_KG_M3 (single-column water loading).
    pub const SEAWATER_DENSITY_KG_M3: f64 = 1030.0;
}

/// Units a v1 constant may carry. Serialized as the same unit strings the
/// artifact assumption entries use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum UnitV1 {
    #[serde(rename = "m")]
    Meter,
    #[serde(rename = "m/s2")]
    MeterPerSecondSquared,
    #[serde(rename = "kg/m3")]
    KilogramPerCubicMeter,
    #[serde(rename = "m3")]
    CubicMeter,
    #[serde(rename = "1/yr")]
    PerYear,
    #[serde(rename = "m2/yr")]
    SquareMeterPerYear,
    #[serde(rename = "1")]
    Dimensionless,
    // CLIMATE-V0 (lane CL, 2026-10-12): additive units for the climate
    // constants below. No existing constant uses them, so every v1
    // universe serializes exactly as before.
    #[serde(rename = "W")]
    Watt,
    #[serde(rename = "s")]
    Second,
    #[serde(rename = "deg")]
    Degree,
    #[serde(rename = "degC")]
    DegreeCelsius,
    #[serde(rename = "W/m2")]
    WattPerSquareMeter,
    #[serde(rename = "W/m2/K")]
    WattPerSquareMeterKelvin,
    #[serde(rename = "K/m")]
    KelvinPerMeter,
    #[serde(rename = "1/K")]
    PerKelvin,
    #[serde(rename = "kg/m2")]
    KilogramPerSquareMeter,
    #[serde(rename = "m/s")]
    MeterPerSecond,
    // SEED-OROGEN (lane OG, 2026-10-12): additive unit for the lower-crust
    // channel viscosity. No pre-existing constant uses it.
    #[serde(rename = "Pa s")]
    PascalSecond,
    // SEED-EXTRAS (lane EX, 2026-10-12): additive units for the plume and
    // rift constants. No pre-existing constant uses them.
    #[serde(rename = "m3/My")]
    CubicMeterPerMegayear,
    #[serde(rename = "My")]
    Megayear,
    #[serde(rename = "deg/My")]
    DegreePerMegayear,
    /// SEED-UPLANDS (lane UP): knickpoint celerity.
    MeterPerMegayear,
    // SEED-TRENCH (lane TN, 2026-10-13): additive unit for the plate
    // elastic modulus. No pre-existing constant uses it.
    #[serde(rename = "Pa")]
    Pascal,
    // SEED-STRATIGRAPHY (lane ST, 2026-10-13): additive unit for the river
    // mouth drainage-area threshold. No pre-existing constant uses it.
    #[serde(rename = "m2")]
    SquareMeter,
    // GLACIAL-T0 (lane GLACIAL-T0, 2026-10-21): additive unit for
    // temperature DIFFERENCES (climate-state offsets). No pre-existing
    // constant uses it.
    #[serde(rename = "K")]
    Kelvin,
    // HYPSO-LEVERS (2026-10-21): additive unit for the erosive-runoff cap.
    // No pre-existing constant uses it.
    #[serde(rename = "m/yr")]
    MeterPerYear,
}

impl UnitV1 {
    pub fn as_str(self) -> &'static str {
        match self {
            UnitV1::Meter => "m",
            UnitV1::MeterPerSecondSquared => "m/s2",
            UnitV1::KilogramPerCubicMeter => "kg/m3",
            UnitV1::CubicMeter => "m3",
            UnitV1::PerYear => "1/yr",
            UnitV1::SquareMeterPerYear => "m2/yr",
            UnitV1::Dimensionless => "1",
            UnitV1::Watt => "W",
            UnitV1::Second => "s",
            UnitV1::Degree => "deg",
            UnitV1::DegreeCelsius => "degC",
            UnitV1::WattPerSquareMeter => "W/m2",
            UnitV1::WattPerSquareMeterKelvin => "W/m2/K",
            UnitV1::KelvinPerMeter => "K/m",
            UnitV1::PerKelvin => "1/K",
            UnitV1::KilogramPerSquareMeter => "kg/m2",
            UnitV1::MeterPerSecond => "m/s",
            UnitV1::PascalSecond => "Pa s",
            UnitV1::CubicMeterPerMegayear => "m3/My",
            UnitV1::Megayear => "My",
            UnitV1::DegreePerMegayear => "deg/My",
            UnitV1::MeterPerMegayear => "m/My",
            UnitV1::Pascal => "Pa",
            UnitV1::SquareMeter => "m2",
            UnitV1::Kelvin => "K",
            UnitV1::MeterPerYear => "m/yr",
        }
    }
}

/// Where a value comes from (ADR 026 §7, ADR 032 §3, ADR 039 §2, C5).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "class", rename_all = "snake_case", deny_unknown_fields)]
pub enum SourceClassV1 {
    /// Literature value with a full citation (as recorded in code).
    Cited { citation: Cow<'static, str> },
    /// Derived from other inputs by a stated relation.
    Derived { from: Cow<'static, str> },
    /// Chosen by the project, or recalled without a checked citation;
    /// awaiting review (ADR 039 C5).
    ProjectAssumptionAwaitingReview { basis: Cow<'static, str> },
    /// A value the legacy presets or code carried before WS-CONST, kept so
    /// old recipes keep their meaning (ADR 039 §2).
    LegacyPresetValue { origin: Cow<'static, str> },
}

/// One constant: value, unit, source. Validity comes from `constant_spec`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConstantV1 {
    pub value: f64,
    pub unit: UnitV1,
    pub source: SourceClassV1,
}

impl ConstantV1 {
    pub fn new(value: f64, unit: UnitV1, source: SourceClassV1) -> Self {
        Self {
            value,
            unit,
            source,
        }
    }
}

/// Code-owned validity envelope of one constant id (inclusive bounds).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantSpecV1 {
    pub id: &'static str,
    pub unit: UnitV1,
    pub lo: f64,
    pub hi: f64,
    pub why: &'static str,
}

/// The v1 envelope table. Bounds that existed in
/// `PlanetParameterProfile::validate` are repeated exactly; the crust
/// thickness bounds are new (positive, not above the 70 km collision cap
/// A-UPLIFT-COLLISION-CAP).
pub const CONSTANT_SPECS_V1: &[ConstantSpecV1] = &[
    // SEED-MAGMATISM (lane MG) revision 1: arc crust cap.
    ConstantSpecV1 {
        id: "arc_crust_cap_m",
        unit: UnitV1::Meter,
        lo: 3.0e4,
        hi: 1.0e5,
        why: "thickness up to which arc magma builds a column (thicker arc roots founder)",
    },
    // SEED-MAGMATISM (lane MG, 2026-10-13): slab-arc-backarc-llsvp-v1.
    ConstantSpecV1 {
        id: "slab_depth_below_arc_m",
        unit: UnitV1::Meter,
        lo: 50000.0,
        hi: 250000.0,
        why: "depth of the slab top beneath the volcanic arc front",
    },
    ConstantSpecV1 {
        id: "slab_dip_young_deg",
        unit: UnitV1::Degree,
        lo: 5.0,
        hi: 80.0,
        why: "effective trench-to-arc slab dip of a young slab",
    },
    ConstantSpecV1 {
        id: "slab_dip_old_deg",
        unit: UnitV1::Degree,
        lo: 5.0,
        hi: 80.0,
        why: "effective trench-to-arc slab dip of an old slab",
    },
    ConstantSpecV1 {
        id: "slab_dip_age_my",
        unit: UnitV1::Megayear,
        lo: 1.0,
        hi: 1000.0,
        why: "slab age at which the dip reaches its old value",
    },
    ConstantSpecV1 {
        id: "arc_magma_productivity_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 7000.0,
        why: "arc crust per unit area of subducted plate (flux = this x clamped convergence rate)",
    },
    ConstantSpecV1 {
        id: "arc_rate_floor_m_per_my",
        unit: UnitV1::MeterPerMegayear,
        lo: 0.0,
        hi: 300000.0,
        why: "convergence rate floor of the arc flux",
    },
    ConstantSpecV1 {
        id: "arc_rate_cap_m_per_my",
        unit: UnitV1::MeterPerMegayear,
        lo: 1000.0,
        hi: 300000.0,
        why: "convergence rate cap of the arc flux",
    },
    ConstantSpecV1 {
        id: "arc_half_width_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 500000.0,
        why: "radius over which arc crust is added around the arc point",
    },
    ConstantSpecV1 {
        id: "arc_continental_thickness_m",
        unit: UnitV1::Meter,
        lo: 10000.0,
        hi: 40000.0,
        why: "thickness at which an oceanic column built mostly of arc crust becomes continental",
    },
    ConstantSpecV1 {
        id: "arc_thermal_uplift_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 3000.0,
        why: "water-loaded uplift of hot arc lithosphere at full arc heat",
    },
    ConstantSpecV1 {
        id: "arc_thermal_half_width_m",
        unit: UnitV1::Meter,
        lo: 10000.0,
        hi: 1000000.0,
        why: "gaussian width of the arc heat on the overriding plate",
    },
    ConstantSpecV1 {
        id: "arc_thermal_decay_my",
        unit: UnitV1::Megayear,
        lo: 1.0,
        hi: 200.0,
        why: "e-folding time of arc heat after the arc moves or stops",
    },
    ConstantSpecV1 {
        id: "backarc_slab_age_min_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 500.0,
        why: "minimum slab age for rollback and back-arc extension",
    },
    ConstantSpecV1 {
        id: "backarc_overriding_advance_max_m_per_my",
        unit: UnitV1::MeterPerMegayear,
        lo: -100000.0,
        hi: 100000.0,
        why: "largest overriding-plate advance toward the trench that still allows back-arc extension",
    },
    ConstantSpecV1 {
        id: "backarc_extension_rate_m_per_my",
        unit: UnitV1::MeterPerMegayear,
        lo: 1000.0,
        hi: 200000.0,
        why: "full back-arc opening rate",
    },
    ConstantSpecV1 {
        id: "backarc_axis_offset_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 1000000.0,
        why: "distance of the back-arc axis behind the arc",
    },
    ConstantSpecV1 {
        id: "backarc_half_width_m",
        unit: UnitV1::Meter,
        lo: 10000.0,
        hi: 1000000.0,
        why: "half width of the extending back-arc band",
    },
    ConstantSpecV1 {
        id: "backarc_min_crust_m",
        unit: UnitV1::Meter,
        lo: 5000.0,
        hi: 30000.0,
        why: "thinnest stretched continental back-arc crust",
    },
    ConstantSpecV1 {
        id: "llsvp_edge_deg",
        unit: UnitV1::Degree,
        lo: 5.0,
        hi: 90.0,
        why: "angular radius of the LLSVP plume generation zone",
    },
    ConstantSpecV1 {
        id: "llsvp_edge_sigma_deg",
        unit: UnitV1::Degree,
        lo: 0.0,
        hi: 45.0,
        why: "spread of plume sources about the LLSVP edge",
    },
    ConstantSpecV1 {
        id: "llsvp_antipode_jitter_deg",
        unit: UnitV1::Degree,
        lo: 0.0,
        hi: 90.0,
        why: "largest offset of the second LLSVP from the first one's antipode",
    },
    ConstantSpecV1 {
        id: "plume_background_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "share of plumes placed uniformly instead of at an LLSVP edge",
    },
    ConstantSpecV1 {
        id: "plume_drift_rate_m_per_my",
        unit: UnitV1::MeterPerMegayear,
        lo: 0.0,
        hi: 100000.0,
        why: "mean plume drift rate in the mantle frame",
    },
    ConstantSpecV1 {
        id: "vent_active_window_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 100.0,
        why: "a vent is active if it erupted within this time of the snapshot",
    },
    ConstantSpecV1 {
        id: "volcano_spacing_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 1000000.0,
        why: "minimum spacing of arc and hotspot volcano points",
    },
    ConstantSpecV1 {
        id: "planet_radius_m",
        unit: UnitV1::Meter,
        lo: 1.0e5,
        hi: 1.0e8,
        why: "planetary radius",
    },
    ConstantSpecV1 {
        id: "surface_gravity_m_s2",
        unit: UnitV1::MeterPerSecondSquared,
        lo: 0.05,
        hi: 100.0,
        why: "surface gravity",
    },
    ConstantSpecV1 {
        id: "mantle_density_kg_m3",
        unit: UnitV1::KilogramPerCubicMeter,
        lo: 3000.0,
        hi: 6000.0,
        why: "mantle density",
    },
    ConstantSpecV1 {
        id: "continental_crust_density_kg_m3",
        unit: UnitV1::KilogramPerCubicMeter,
        lo: 2000.0,
        hi: 3500.0,
        why: "crust density below mantle",
    },
    ConstantSpecV1 {
        id: "oceanic_crust_density_kg_m3",
        unit: UnitV1::KilogramPerCubicMeter,
        lo: 2000.0,
        hi: 3500.0,
        why: "oceanic crust density",
    },
    ConstantSpecV1 {
        id: "seawater_density_kg_m3",
        unit: UnitV1::KilogramPerCubicMeter,
        lo: 900.0,
        hi: 1300.0,
        why: "liquid water density",
    },
    ConstantSpecV1 {
        id: "continental_crust_thickness_m",
        unit: UnitV1::Meter,
        lo: 1.0e3,
        hi: 7.0e4,
        why: "positive and not above the 70 km collision cap",
    },
    ConstantSpecV1 {
        id: "oceanic_crust_thickness_m",
        unit: UnitV1::Meter,
        lo: 1.0e3,
        hi: 7.0e4,
        why: "positive and not above the 70 km collision cap",
    },
    ConstantSpecV1 {
        id: "reference_crust_thickness_m",
        unit: UnitV1::Meter,
        lo: 1.0e3,
        hi: 7.0e4,
        why: "single-column Airy reference column",
    },
    ConstantSpecV1 {
        id: "stream_power_erodibility_per_yr",
        unit: UnitV1::PerYear,
        lo: 0.0,
        hi: 1.0e-3,
        why: "erodibility",
    },
    ConstantSpecV1 {
        id: "hillslope_diffusivity_m2_per_yr",
        unit: UnitV1::SquareMeterPerYear,
        lo: 0.0,
        hi: 1.0e3,
        why: "hillslope diffusivity",
    },
    ConstantSpecV1 {
        id: "ocean_water_volume_m3",
        unit: UnitV1::CubicMeter,
        lo: 0.0,
        hi: 1.0e22,
        why: "liquid ocean water volume",
    },
    ConstantSpecV1 {
        id: "sea_level_datum_m",
        unit: UnitV1::Meter,
        lo: -1.0e4,
        hi: 1.0e4,
        why: "datum",
    },
    ConstantSpecV1 {
        id: "fixed_compensation_depth_m",
        unit: UnitV1::Meter,
        lo: 1.0e4,
        hi: 1.0e6,
        why: "compensation depth below the mantle free-surface datum",
    },
    ConstantSpecV1 {
        id: "derived_depth_margin_m",
        unit: UnitV1::Meter,
        lo: 1.0,
        hi: 1.0e6,
        why: "margin of the derived compensation depth below the deepest root",
    },
    ConstantSpecV1 {
        id: "continental_area_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.01,
        hi: 0.95,
        why: "continental crust area fraction",
    },
    ConstantSpecV1 {
        id: "rifted_margin_width_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 2.0e6,
        why: "rifted margin width",
    },
    // SEED-OROGEN (lane OG, 2026-10-12): rate-limited-belt-v1 and
    // stream-power-subgrid-v1 parameters. Additive rows; envelopes are
    // project assumptions about where each named approximation stays
    // meaningful (the `why` names it).
    ConstantSpecV1 {
        id: "belt_shear_zone_dip_deg",
        unit: UnitV1::Degree,
        lo: 10.0,
        hi: 60.0,
        why: "dip of the conjugate crustal step-up shear zones of a two-sided orogenic wedge (Coulomb shear zones; outside this range the wedge-width relation W = 2 h cot(dip) is not a crustal-scale geometry)",
    },
    ConstantSpecV1 {
        id: "lower_crust_channel_viscosity_pa_s",
        unit: UnitV1::PascalSecond,
        lo: 1.0e17,
        hi: 1.0e23,
        why: "effective viscosity of a ductile lower-crustal channel (Poiseuille channel-flow spreading)",
    },
    ConstantSpecV1 {
        id: "channel_reference_thickness_m",
        unit: UnitV1::Meter,
        lo: 1.0e4,
        hi: 7.0e4,
        why: "crust thickness above which the excess forms a flowing lower-crustal channel",
    },
    ConstantSpecV1 {
        id: "delamination_thickness_m",
        unit: UnitV1::Meter,
        lo: 4.0e4,
        hi: 1.5e5,
        why: "crust thickness above which the excess root is removed (eclogitised-root delamination backstop)",
    },
    ConstantSpecV1 {
        id: "arc_productivity_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 7.0e3,
        why: "arc crust volume per unit area of subducted plate (at most the oceanic crust thickness)",
    },
    ConstantSpecV1 {
        id: "arc_belt_width_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 1.0e6,
        why: "width over which arc magmatic addition is spread on the overriding plate",
    },
    ConstantSpecV1 {
        id: "subgrid_channel_head_length_m",
        unit: UnitV1::Meter,
        lo: 1.0,
        hi: 1.0e5,
        why: "length scale of the smallest channels of the unresolved (subgrid) drainage network (hillslope length / channel head)",
    },
    // SEED-RELIEF-LAW (lane RL, 2026-10-12): stream-power-subgrid-
    // nonlinear-v1 parameters. Additive rows.
    ConstantSpecV1 {
        id: "stream_power_slope_exponent_n",
        unit: UnitV1::Dimensionless,
        lo: 1.0,
        hi: 4.0,
        why: "slope exponent n of mean-flow stream power (m = n/2); 1 is the linear legacy law, values above ~4 are outside the stochastic-threshold range the law stands for",
    },
    ConstantSpecV1 {
        id: "stream_power_reference_steepness_m",
        unit: UnitV1::Meter,
        lo: 1.0,
        hi: 1.0e5,
        why: "channel steepness k_s = S A^0.5 at which the n-law erodes at the n = 1 law's rate (K_n = K k_ref^(1 - n))",
    },
    // SEED-EXTRAS D7 (lane SD, 2026-10-12): erosion-deposition-v1 and
    // shelf-wave-base-v1 parameters. Additive rows.
    ConstantSpecV1 {
        id: "deposition_coefficient_g",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 100.0,
        why: "dimensionless deposition coefficient G of the xi-q / Yuan et al. 2019 erosion-deposition law (deposition rate = G Qs / (p A))",
    },
    ConstantSpecV1 {
        id: "shelf_wave_base_depth_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 500.0,
        why: "depth below sea level to which the shelf top aggrades (storm wave base / cell-mean shelf-top depth)",
    },
    // SEED-SEDIMENT-RETAIN (lane SR, 2026-10-13):
    // continental-retention-base-level-v1 parameters. Additive rows.
    ConstantSpecV1 {
        id: "base_level_relief_scale_m",
        unit: UnitV1::Meter,
        lo: 10.0,
        hi: 1.0e4,
        why: "local relief (max - min level over a cell and its 8 neighbours) at and above which stream-power lowering is not tapered (phi = min(1, (r/H)^p))",
    },
    ConstantSpecV1 {
        id: "base_level_taper_exponent",
        unit: UnitV1::Dimensionless,
        lo: 0.25,
        hi: 8.0,
        why: "exponent p of the base-level erosion taper phi = min(1, (h/H)^p)",
    },
    ConstantSpecV1 {
        id: "alluvial_grade_slope",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0e-2,
        why: "graded alluvial slope of cell means up to which continental lowlands store passing sediment",
    },
    // SEED-EXTRAS (lane EX, 2026-10-12): fixed-plume-chains-v1,
    // stretched-rift-v1 and plume-head-continental-split-v1 parameters.
    // Additive rows.
    ConstantSpecV1 {
        id: "plume_count",
        unit: UnitV1::Dimensionless,
        lo: 1.0,
        hi: 100.0,
        why: "number of fixed mantle plumes (whole number)",
    },
    ConstantSpecV1 {
        id: "plume_onset_lead_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 1.0e4,
        why: "plume onsets are drawn uniformly in [-lead, duration); an onset before t = 0 is a plume already active, without a head in the run",
    },
    ConstantSpecV1 {
        id: "plume_head_volume_m3",
        unit: UnitV1::CubicMeter,
        lo: 0.0,
        hi: 1.0e18,
        why: "magmatic crust volume of a plume head (large igneous province) added once at onset",
    },
    ConstantSpecV1 {
        id: "plume_head_radius_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 3.0e6,
        why: "radius of the disc over which a plume head's volume is spread",
    },
    ConstantSpecV1 {
        id: "plume_tail_flux_m3_per_my",
        unit: UnitV1::CubicMeterPerMegayear,
        lo: 0.0,
        hi: 1.0e16,
        why: "steady magmatic crust production of a plume tail",
    },
    ConstantSpecV1 {
        id: "plume_edifice_radius_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 1.0e6,
        why: "radius of the disc over which each tail sample is spread on the overlying plate (edifice plus swell footprint)",
    },
    ConstantSpecV1 {
        id: "rift_half_width_m",
        unit: UnitV1::Meter,
        lo: 1.0e3,
        hi: 1.0e6,
        why: "half width of the continental zone that stretches to fill a rift birth",
    },
    ConstantSpecV1 {
        id: "breakup_crust_thickness_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 3.0e4,
        why: "a rift birth thinner than this is an oceanic birth instead (continental breakup)",
    },
    ConstantSpecV1 {
        id: "rift_opening_rate_deg_per_my",
        unit: UnitV1::DegreePerMegayear,
        lo: 0.0,
        hi: 10.0,
        why: "relative angular rate of a plate split off by a plume head (opening rate at the plume = rate x planet radius)",
    },
    // SEED-RIFT-SHOULDERS (lane RS, 2026-10-13): shoulder-uplift-thermal-v1
    // parameters. Additive rows.
    ConstantSpecV1 {
        id: "rift_flank_flexural_parameter_m",
        unit: UnitV1::Meter,
        lo: 1.0e3,
        hi: 1.0e6,
        why: "flexural parameter alpha of the broken plate at a rift border fault: the shoulder profile is exp(-x/alpha) cos(x/alpha) out to x = pi alpha / 2 from the rift edge",
    },
    ConstantSpecV1 {
        id: "rift_shoulder_ratio",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 2.0,
        why: "shoulder uplift at the rift edge as a fraction of the adjacent rift floor's thermal uplift",
    },
    // SEED-UPLANDS (lane UP, 2026-10-13): critical-slope-envelope-v1
    // parameters. Additive rows.
    ConstantSpecV1 {
        id: "hillslope_critical_slope",
        unit: UnitV1::Dimensionless,
        lo: 1.0e-3,
        hi: 2.0,
        why: "critical slope of CELL MEANS at the grid's scale: the steepest cell-mean gradient threshold landscapes sustain (about 0.6-0.7 at the hillslope scale, a few percent for 0.5-1 deg means)",
    },
    // INCISION-THRESHOLD (lane IT, 2026-10-20; merged by lane
    // CLIMATE-GENESIS-2). Additive row.
    ConstantSpecV1 {
        id: "critical_channel_steepness_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 500.0,
        why: "channel steepness k_s = S A^(1/2) (theta = 0.5 units) below which floods no longer detach bedrock (excess-stream-power incision threshold); kept below the steepness of still-incising channels",
    },
    // GLACIAL-T0 (lane GLACIAL-T0, 2026-10-21): glacial-stream-power-blend-v1.
    // Additive rows.
    ConstantSpecV1 {
        id: "glacial_to_fluvial_erodibility_ratio",
        unit: UnitV1::Dimensionless,
        lo: 0.1,
        hi: 100.0,
        why: "K_g / K_f of the blended glacial + subglacial-fluvial stream-power law (Hergarten 2021 eq 27): erosion of a glacier relative to a river at the same flux (catchment-size equivalent) and slope",
    },
    ConstantSpecV1 {
        id: "full_ice_height_above_ela_m",
        unit: UnitV1::Meter,
        lo: 10.0,
        hi: 5000.0,
        why: "height above the equilibrium-line altitude at which all precipitation becomes ice (linear ice-production ramp, Hergarten 2021 eq 23)",
    },
    ConstantSpecV1 {
        id: "glacial_state_duty",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "fraction of genesis time spent in one declared climate state",
    },
    ConstantSpecV1 {
        id: "glacial_state_offset_k",
        unit: UnitV1::Kelvin,
        lo: 0.0,
        hi: 30.0,
        why: "temperature offset of a climate state from the genesis climate (cooling of the glacial maximum, or warming of a warm state)",
    },
    ConstantSpecV1 {
        id: "glacial_max_precipitation_factor",
        unit: UnitV1::Dimensionless,
        lo: 0.1,
        hi: 2.0,
        why: "precipitation at the glacial maximum relative to the genesis climate",
    },
    ConstantSpecV1 {
        id: "knickpoint_celerity_m_per_my",
        unit: UnitV1::MeterPerMegayear,
        lo: 1.0e2,
        hi: 1.0e7,
        why: "rate at which an over-steepened front is cut back; the threshold response time of a cell is its linear size over this celerity",
    },
    // CLIMATE-GENESIS (lane CG, 2026-10-20): epoch-climate-runoff-v1
    // numerics and runoff rule. Additive rows.
    ConstantSpecV1 {
        id: "climate_epoch_my",
        unit: UnitV1::Megayear,
        lo: 2.0,
        hi: 20.0,
        why: "interval between genesis climate solves; ranges rise and fall over 20-30 My, so the epoch stays well below that (lag oscillation) and at or above one 2 My step",
    },
    ConstantSpecV1 {
        id: "climate_solve_grid_divisor",
        unit: UnitV1::Dimensionless,
        lo: 1.0,
        hi: 4.0,
        why: "the climate solve grid is the genesis grid coarsened by this factor in each axis (1, 2 or 4); numerics, not physics",
    },
    ConstantSpecV1 {
        id: "climate_sweep_rel_tol",
        unit: UnitV1::Dimensionless,
        lo: 1.0e-11,
        hi: 1.0e-5,
        why: "line Gauss-Seidel stop tolerance of a genesis climate solve; the water closure at 1e-6 is still checked on every solve",
    },
    // HYPSO-LEVERS (2026-10-21): erosive-runoff response of the climate
    // forcing law (opt-in). Additive rows.
    ConstantSpecV1 {
        id: "erosive_runoff_exponent",
        unit: UnitV1::Dimensionless,
        lo: 0.05,
        hi: 1.0,
        why: "exponent p of the erosive runoff above the reference runoff, R_e = r_ref (R / r_ref)^p for R > r_ref: 1 is the linear v1 response, smaller values make orographic runoff less erosive",
    },
    ConstantSpecV1 {
        id: "erosive_runoff_cap_m_per_yr",
        unit: UnitV1::MeterPerYear,
        lo: 0.01,
        hi: 100.0,
        why: "cap of the erosive runoff, R_e = min(R, cap): runoff above it adds no erosive power",
    },
    ConstantSpecV1 {
        id: "climate_runoff_floor_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "land runoff floor R >= f P (land_p_minus_e_floor_v1): the share of precipitation that runs off even where evapotranspiration could take it all",
    },
    // SEED-TRENCH (lane TN, 2026-10-13): slab-flexure-trench-v1
    // parameters. Additive rows.
    ConstantSpecV1 {
        id: "trench_vertical_load_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "share of the slab pull carried as the vertical end load of the bending downgoing plate at the trench",
    },
    ConstantSpecV1 {
        id: "trench_slab_length_max_m",
        unit: UnitV1::Meter,
        lo: 1.0e4,
        hi: 3.0e6,
        why: "longest slab whose negative thermal buoyancy pulls on the plate (upper-mantle slab)",
    },
    ConstantSpecV1 {
        id: "trench_thermal_expansivity_per_k",
        unit: UnitV1::PerKelvin,
        lo: 1.0e-6,
        hi: 1.0e-4,
        why: "volumetric thermal expansivity of the oceanic lithosphere (slab negative buoyancy)",
    },
    ConstantSpecV1 {
        id: "trench_mantle_temperature_c",
        unit: UnitV1::DegreeCelsius,
        lo: 500.0,
        hi: 2500.0,
        why: "plate base minus surface temperature of the half-space cooling plate",
    },
    ConstantSpecV1 {
        id: "trench_elastic_isotherm_c",
        unit: UnitV1::DegreeCelsius,
        lo: 100.0,
        hi: 1000.0,
        why: "isotherm whose depth is the plate's effective elastic thickness",
    },
    ConstantSpecV1 {
        id: "trench_thermal_diffusivity_m2_per_yr",
        unit: UnitV1::SquareMeterPerYear,
        lo: 1.0,
        hi: 1000.0,
        why: "thermal diffusivity of the oceanic lithosphere (31.5576 m2/yr = 1e-6 m2/s)",
    },
    ConstantSpecV1 {
        id: "trench_elastic_modulus_pa",
        unit: UnitV1::Pascal,
        lo: 1.0e9,
        hi: 1.0e12,
        why: "Young's modulus of the elastic plate",
    },
    ConstantSpecV1 {
        id: "trench_poisson_ratio",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 0.49,
        why: "Poisson's ratio of the elastic plate",
    },
    // SEED-FLOODBASALT (lane FB, 2026-10-13): plume-head-lip-swell-v1
    // parameters. Additive rows.
    ConstantSpecV1 {
        id: "flood_basalt_extrusive_volume_m3",
        unit: UnitV1::CubicMeter,
        lo: 0.0,
        hi: 1.0e18,
        why: "erupted (extrusive) flood-basalt volume of one plume-head large igneous province",
    },
    ConstantSpecV1 {
        id: "flood_basalt_intrusive_ratio",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 20.0,
        why: "intrusive plus underplated magma per unit extrusive volume of a large igneous province",
    },
    ConstantSpecV1 {
        id: "flood_basalt_radius_m",
        unit: UnitV1::Meter,
        lo: 1.0e5,
        hi: 3.0e6,
        why: "radius of the flood-basalt province (the flattened plume head)",
    },
    ConstantSpecV1 {
        id: "flood_basalt_taper_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.01,
        hi: 1.0,
        why: "outer share of the province radius over which the lava pile thins to zero",
    },
    ConstantSpecV1 {
        id: "flood_basalt_eruption_duration_my",
        unit: UnitV1::Megayear,
        lo: 0.1,
        hi: 50.0,
        why: "duration of the main flood-basalt eruption phase",
    },
    ConstantSpecV1 {
        id: "swell_tail_amplitude_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 5000.0,
        why: "peak water-loaded dynamic topography of a plume tail's swell",
    },
    ConstantSpecV1 {
        id: "swell_tail_radius_m",
        unit: UnitV1::Meter,
        lo: 5.0e4,
        hi: 3.0e6,
        why: "e-folding radius of the plume tail swell (Gaussian)",
    },
    ConstantSpecV1 {
        id: "swell_head_amplitude_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 5000.0,
        why: "extra peak water-loaded dynamic topography of a newly arrived plume head",
    },
    ConstantSpecV1 {
        id: "swell_head_radius_m",
        unit: UnitV1::Meter,
        lo: 5.0e4,
        hi: 3.0e6,
        why: "e-folding radius of the plume head swell (Gaussian)",
    },
    ConstantSpecV1 {
        id: "swell_head_decay_my",
        unit: UnitV1::Megayear,
        lo: 0.5,
        hi: 500.0,
        why: "e-folding time of the plume head swell after onset",
    },
    // SEED-CARRY-CONSERVE (lane CV, 2026-10-13): continental-parcel-carry-v1
    // and contact-length-pole-share-v1 parameters. Additive rows.
    ConstantSpecV1 {
        id: "parcel_overlap_cell_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.05,
        hi: 1.0,
        why: "overlap beyond a full cell at which the smaller continental depositors' parcels in the cell are consumed by collision",
    },
    ConstantSpecV1 {
        id: "plate_fragment_min_area_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 0.05,
        why: "connected plate-id pieces smaller than this share of the sphere fold into the plate they touch most",
    },
    // HYPSO-LEVERS (2026-10-21): opt-in write-back residual correction of
    // continental-parcel-carry-v1. Additive row.
    ConstantSpecV1 {
        id: "write_back_residual_iterations",
        unit: UnitV1::Dimensionless,
        lo: 1.0,
        hi: 16.0,
        why: "whole number of re-projection passes after the parcel write-back that hand the residual (end-of-step cell thickness minus the parcels' re-projection, area-weighted mean removed) back to the parcels",
    },
    ConstantSpecV1 {
        id: "suture_contact_length_m",
        unit: UnitV1::Meter,
        lo: 1.0e4,
        hi: 2.0e7,
        why: "converging continent-continent contact length at which two plates suture (the smaller takes the larger's pole)",
    },
    // SEED-PLATE-REORG (lane PR, 2026-10-13): plate-reorg-scripted-v1
    // parameters. Additive rows.
    ConstantSpecV1 {
        id: "reorg_interval_min_my",
        unit: UnitV1::Megayear,
        lo: 1.0,
        hi: 1000.0,
        why: "shortest gap between scripted pole re-solves",
    },
    ConstantSpecV1 {
        id: "reorg_interval_max_my",
        unit: UnitV1::Megayear,
        lo: 1.0,
        hi: 1000.0,
        why: "longest gap between scripted pole re-solves",
    },
    ConstantSpecV1 {
        id: "slab_pull_bias",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "weight of the slab-torque axis against a random axis at a pole re-solve",
    },
    ConstantSpecV1 {
        id: "slab_circumference_reference",
        unit: UnitV1::Dimensionless,
        lo: 0.01,
        hi: 1.0,
        why: "subducting share of a group's boundary at which its re-solved rate saturates",
    },
    ConstantSpecV1 {
        id: "reorg_rate_floor_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "re-solved rate of a slab-free group as a fraction of the rate scale",
    },
    ConstantSpecV1 {
        id: "supercontinent_share",
        unit: UnitV1::Dimensionless,
        lo: 0.05,
        hi: 1.0,
        why: "share of the continental crust in one suture group that makes it a supercontinent",
    },
    ConstantSpecV1 {
        id: "supercontinent_tenure_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 2000.0,
        why: "time a supercontinent holds together before it breaks up",
    },
    ConstantSpecV1 {
        id: "breakup_piece_max_share",
        unit: UnitV1::Dimensionless,
        lo: 0.05,
        hi: 1.0,
        why: "a breakup piece holding more than this share of the continental crust breaks again",
    },
    ConstantSpecV1 {
        id: "breakup_cascade_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 500.0,
        why: "delay before an oversized breakup piece breaks again",
    },
    ConstantSpecV1 {
        id: "rift_persistence_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 2000.0,
        why: "age below which a breakup or initiation keeps its relative rotation across re-solves",
    },
    ConstantSpecV1 {
        id: "weak_zone_strength_ratio",
        unit: UnitV1::Dimensionless,
        lo: 0.01,
        hi: 1.0,
        why: "strength of the weakest zones (old sutures, belts, thinned crust) relative to normal crust in the breakup cut",
    },
    ConstantSpecV1 {
        id: "plate_count_min",
        unit: UnitV1::Dimensionless,
        lo: 1.0,
        hi: 100.0,
        why: "fewest plate groups above the counting area before subduction initiation splits one",
    },
    ConstantSpecV1 {
        id: "plate_count_max",
        unit: UnitV1::Dimensionless,
        lo: 1.0,
        hi: 100.0,
        why: "most plate groups above the counting area before the smallest oceanic one folds",
    },
    ConstantSpecV1 {
        id: "plate_count_min_area_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 0.2,
        why: "share of the sphere above which a plate group counts (and below which cut leftovers fold)",
    },
    ConstantSpecV1 {
        id: "subduction_initiation_rate_deg_per_my",
        unit: UnitV1::DegreePerMegayear,
        lo: 0.0,
        hi: 5.0,
        why: "relative rotation rate at which a new subduction zone converges",
    },
    ConstantSpecV1 {
        id: "plume_head_interval_my",
        unit: UnitV1::Megayear,
        lo: 1.0,
        hi: 1000.0,
        why: "mean gap between plume head onsets (renewal process)",
    },
    ConstantSpecV1 {
        id: "plume_tail_lifetime_my",
        unit: UnitV1::Megayear,
        lo: 1.0,
        hi: 5000.0,
        why: "time after onset for which a plume tail keeps erupting",
    },
    ConstantSpecV1 {
        id: "plume_split_min_share",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "share of the continental crust a plume's group must hold for its head to split it",
    },
    // SEED-IRREGULAR-SEEDS (lane IS, 2026-10-13): irregular-cratons-v1
    // parameters. Additive rows.
    ConstantSpecV1 {
        id: "craton_interior_thickness_m",
        unit: UnitV1::Meter,
        lo: 15_000.0,
        hi: 60_000.0,
        why: "untapered interior continental crust thickness of the initial cratons",
    },
    ConstantSpecV1 {
        id: "craton_rank_size_exponent",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 4.0,
        why: "rank-size exponent of continental assembly and craton area shares (share of rank k ~ k^-a)",
    },
    ConstantSpecV1 {
        id: "cratons_per_unit_share",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 100.0,
        why: "cratons in an assembly = 1 + floor(this x the assembly's share of the continental area)",
    },
    ConstantSpecV1 {
        id: "basement_roughness_std_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 10_000.0,
        why: "standard deviation of the multi-scale basement thickness roughness (basins and arches)",
    },
    ConstantSpecV1 {
        id: "suture_crest_thickening_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 30_000.0,
        why: "crest crustal thickening of an old suture root between initial cratons",
    },
    ConstantSpecV1 {
        id: "suture_half_width_m",
        unit: UnitV1::Meter,
        lo: 1.0e4,
        hi: 1.0e6,
        why: "Gaussian half-width (sigma) of an old suture root",
    },
    ConstantSpecV1 {
        id: "proximal_margin_stretch",
        unit: UnitV1::Dimensionless,
        lo: 1.0,
        hi: 3.0,
        why: "stretch factor of the proximal rifted margin at the necking hinge (shelf break)",
    },
    // SEED-RANGES (lane RG, 2026-10-13).
    ConstantSpecV1 {
        id: "indentation_time_my",
        unit: UnitV1::Megayear,
        lo: 1.0,
        hi: 200.0,
        why: "time over which post-collision indentation tapers to zero",
    },
    ConstantSpecV1 {
        id: "indentation_shortening_budget_m",
        unit: UnitV1::Meter,
        lo: 1.0e4,
        hi: 5.0e6,
        why: "cumulative post-collision shortening at a suture contact that ends the indentation",
    },
    ConstantSpecV1 {
        id: "plateau_crust_cap_m",
        unit: UnitV1::Meter,
        lo: 4.0e4,
        hi: 1.0e5,
        why: "continental crust thickness above which thickened crust spills into the hinterland",
    },
    ConstantSpecV1 {
        id: "plateau_spill_reach_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 5.0e6,
        why: "largest distance from its source that spilled plateau crust may travel",
    },
    // SEED-WORLD-COMPOSE (lane CP, 2026-10-13).
    ConstantSpecV1 {
        id: "plateau_foothill_gradient",
        unit: UnitV1::Dimensionless,
        lo: 0.01,
        hi: 1.0,
        why: "largest crust thickness step per unit horizontal distance at a plateau margin (foothill taper)",
    },
    ConstantSpecV1 {
        id: "plateau_taper_floor_m",
        unit: UnitV1::Meter,
        lo: 2.0e4,
        hi: 1.0e5,
        why: "crust thickness above which the foothill taper acts on the thicker cell of a pair",
    },
    ConstantSpecV1 {
        id: "backarc_episode_my",
        unit: UnitV1::Megayear,
        lo: 1.0,
        hi: 200.0,
        why: "time after a cell's first back-arc extension during which it may still be extended",
    },
    ConstantSpecV1 {
        id: "backarc_max_stretch_factor",
        unit: UnitV1::Dimensionless,
        lo: 1.0,
        hi: 10.0,
        why: "largest total stretching factor of a continental cell in one back-arc episode",
    },
    ConstantSpecV1 {
        id: "kimberlite_plume_radius_m",
        unit: UnitV1::Meter,
        lo: 1.0e4,
        hi: 2.0e6,
        why: "distance from an arriving plume head within which continental columns get the kimberlite plume-pass marker",
    },
    ConstantSpecV1 {
        id: "subduction_initiation_age_my",
        unit: UnitV1::Megayear,
        lo: 50.0,
        hi: 500.0,
        why: "ocean-crust age above which old ocean triggers subduction initiation",
    },
    // SEED-TRANSFORMS (lane TF, 2026-10-13): smoothed-normal-transform-v1 and
    // plate-reorg-scripted-v4.
    ConstantSpecV1 {
        id: "transform_obliquity_deg",
        unit: UnitV1::Degree,
        lo: 1.0,
        hi: 89.0,
        why: "angle between the relative plate velocity and the boundary trace below which a moving plate boundary is a transform",
    },
    ConstantSpecV1 {
        id: "boundary_normal_smoothing_m",
        unit: UnitV1::Meter,
        lo: 1.0e4,
        hi: 2.0e6,
        why: "Gaussian sigma of the plate-indicator window whose first moment gives the boundary normal",
    },
    ConstantSpecV1 {
        id: "initiation_slab_age_min_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 500.0,
        why: "ocean age that counts as old slab along a new age-triggered subduction boundary",
    },
    ConstantSpecV1 {
        id: "initiation_old_slab_fraction_min",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "least length fraction of a new age-triggered subduction boundary that runs through old slab",
    },
    ConstantSpecV1 {
        id: "initiation_trench_spacing_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 5.0e6,
        why: "distance from existing subduction within which a new age-triggered boundary may run for at most half of its length",
    },
    // SEED-OLD-OCEAN (lane OR, 2026-10-13): plate-reorg-scripted-v5 margin initiation.
    ConstantSpecV1 {
        id: "margin_initiation_age_my",
        unit: UnitV1::Megayear,
        lo: 50.0,
        hi: 500.0,
        why: "passive-margin ocean at least this old seeds a new trench along the margin (oldest first)",
    },
    ConstantSpecV1 {
        id: "margin_segment_length_m",
        unit: UnitV1::Meter,
        lo: 1.0e5,
        hi: 1.0e7,
        why: "longest along-margin length of a passive-margin subduction initiation",
    },
    ConstantSpecV1 {
        id: "margin_trench_spacing_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 5.0e6,
        why: "a passive-margin initiation keeps at least this distance from existing subduction and trenches started in the same pass",
    },
    // SEED-STRATIGRAPHY (lane ST, 2026-10-13): parcel-stacks-v1. Additive rows.
    ConstantSpecV1 {
        id: "strata_max_layers",
        unit: UnitV1::Dimensionless,
        lo: 2.0,
        hi: 64.0,
        why: "layers per continental parcel stack, the basement included",
    },
    ConstantSpecV1 {
        id: "strata_merge_gap_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 50.0,
        why: "age gap within which an appended layer of the same lithology and environment extends the top layer",
    },
    ConstantSpecV1 {
        id: "strata_unconformity_gap_my",
        unit: UnitV1::Megayear,
        lo: 0.0,
        hi: 1000.0,
        why: "age gap between a layer base and the layer below it read as an unconformity",
    },
    ConstantSpecV1 {
        id: "strata_min_layer_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 1000.0,
        why: "cover layer thickness below which a truncated top layer merges downward",
    },
    ConstantSpecV1 {
        id: "strata_trace_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 10.0,
        why: "addition thickness below which it extends the layer it lands on",
    },
    ConstantSpecV1 {
        id: "strata_arc_extrusive_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "share of arc crust added as volcanics on top (the rest is a granite pluton at the basement top)",
    },
    ConstantSpecV1 {
        id: "strata_plume_extrusive_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "share of plume (LIP and hotspot) crust added as basalt on top (the rest joins the basement)",
    },
    ConstantSpecV1 {
        id: "strata_shortening_pluton_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "share of positive shortening thickening added as a granite pluton (the rest thickens every layer)",
    },
    ConstantSpecV1 {
        id: "strata_pluton_roof_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 30000.0,
        why: "emplacement depth of plutons: crystalline basement roof kept above a new pluton (revision 1)",
    },
    ConstantSpecV1 {
        id: "strata_metamorphic_depth_m",
        unit: UnitV1::Meter,
        lo: 1000.0,
        hi: 60000.0,
        why: "base depth below which a cover layer becomes metamorphic (grade 1)",
    },
    ConstantSpecV1 {
        id: "strata_high_grade_depth_m",
        unit: UnitV1::Meter,
        lo: 1000.0,
        hi: 80000.0,
        why: "base depth below which a metamorphic layer is high grade",
    },
    ConstantSpecV1 {
        id: "strata_shelf_depth_m",
        unit: UnitV1::Meter,
        lo: 10.0,
        hi: 2000.0,
        why: "water depth below which a marine deposit is deep-marine shale",
    },
    ConstantSpecV1 {
        id: "strata_carbonate_max_paleolat_deg",
        unit: UnitV1::Degree,
        lo: 0.0,
        hi: 90.0,
        why: "latitude fallback for warm water: shelf carbonate only equatorward of this",
    },
    ConstantSpecV1 {
        id: "strata_carbonate_max_supply_m_per_my",
        unit: UnitV1::MeterPerMegayear,
        lo: 0.0,
        hi: 100000.0,
        why: "shelf deposition rate above which the shelf deposit is clastic shale, not carbonate",
    },
    ConstantSpecV1 {
        id: "strata_arid_lat_min_deg",
        unit: UnitV1::Degree,
        lo: 0.0,
        hi: 90.0,
        why: "latitude fallback: arid belt equatorward edge (closed lakes deposit evaporite)",
    },
    ConstantSpecV1 {
        id: "strata_arid_lat_max_deg",
        unit: UnitV1::Degree,
        lo: 0.0,
        hi: 90.0,
        why: "latitude fallback: arid belt poleward edge",
    },
    ConstantSpecV1 {
        id: "strata_humid_lat_max_deg",
        unit: UnitV1::Degree,
        lo: 0.0,
        hi: 90.0,
        why: "latitude fallback: humid tropics poleward edge",
    },
    ConstantSpecV1 {
        id: "strata_humid_midlat_min_deg",
        unit: UnitV1::Degree,
        lo: 0.0,
        hi: 90.0,
        why: "latitude fallback: humid mid-latitude band equatorward edge",
    },
    ConstantSpecV1 {
        id: "strata_humid_midlat_max_deg",
        unit: UnitV1::Degree,
        lo: 0.0,
        hi: 90.0,
        why: "latitude fallback: humid mid-latitude band poleward edge",
    },
    ConstantSpecV1 {
        id: "strata_coal_max_elevation_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 5000.0,
        why: "highest elevation above sea level of a humid coastal land deposit that is coal-bearing",
    },
    ConstantSpecV1 {
        id: "strata_conglomerate_relief_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 10000.0,
        why: "3x3 local relief above which a land deposit is conglomerate",
    },
    ConstantSpecV1 {
        id: "strata_river_mouth_area_m2",
        unit: UnitV1::SquareMeter,
        lo: 0.0,
        hi: 100000000000000.0,
        why: "drainage area of a land neighbour that makes a marine cell a river mouth (deltaic sandstone)",
    },
    ConstantSpecV1 {
        id: "strata_publish_depth_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 100000.0,
        why: "depth to which the published stratigraphy record lists layers",
    },
    ConstantSpecV1 {
        id: "strata_basement_min_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 100000.0,
        why: "basement thickness below which a parcel is counted (count only)",
    },
    // SEED-STRATIGRAPHY stage A.1 (lane SA, 2026-10-13): parcel-stacks-v2.
    // Additive rows.
    ConstantSpecV1 {
        id: "strata_shelf_width_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 2.0e6,
        why: "width of the flat shelf seaward of a coast in the declared sub-cell shelf-slope ramp",
    },
    ConstantSpecV1 {
        id: "strata_slope_foot_depth_m",
        unit: UnitV1::Meter,
        lo: 100.0,
        hi: 11000.0,
        why: "deepest depth of the non-shelf part of a marine cell in the sub-cell shelf-slope ramp (continental rise)",
    },
    ConstantSpecV1 {
        id: "strata_shelf_majority_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "sub-cell shelf fraction at or above which a marine cell's deposit takes the shelf facies",
    },
    ConstantSpecV1 {
        id: "strata_platform_cover_max_m",
        unit: UnitV1::Meter,
        lo: 0.0,
        hi: 20000.0,
        why: "platform cover seeded at t = 0 on parcels at or below the full-cover elevation",
    },
    ConstantSpecV1 {
        id: "strata_platform_full_elev_m",
        unit: UnitV1::Meter,
        lo: -10000.0,
        hi: 10000.0,
        why: "t = 0 elevation above the sea datum at or below which a parcel gets the full platform cover",
    },
    ConstantSpecV1 {
        id: "strata_platform_zero_elev_m",
        unit: UnitV1::Meter,
        lo: -10000.0,
        hi: 10000.0,
        why: "t = 0 elevation above the sea datum at or above which a parcel gets no platform cover (exposed shield)",
    },
    ConstantSpecV1 {
        id: "strata_platform_sandstone_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "basal sandstone share of the t = 0 platform cover",
    },
    ConstantSpecV1 {
        id: "strata_platform_carbonate_fraction",
        unit: UnitV1::Dimensionless,
        lo: 0.0,
        hi: 1.0,
        why: "top carbonate share of the t = 0 platform cover (shale is the rest)",
    },
    ConstantSpecV1 {
        id: "strata_platform_age_base_my",
        unit: UnitV1::Megayear,
        lo: -5000.0,
        hi: 0.0,
        why: "age of the base of the t = 0 platform cover",
    },
    ConstantSpecV1 {
        id: "strata_platform_age_top_my",
        unit: UnitV1::Megayear,
        lo: -5000.0,
        hi: 0.0,
        why: "age of the top of the t = 0 platform cover",
    },
];

/// Look up the code-owned spec of a constant id (panics on an id that is
/// not in the table: that is a programming error, never data).
pub fn constant_spec(id: &str) -> &'static ConstantSpecV1 {
    CONSTANT_SPECS_V1
        .iter()
        .find(|s| s.id == id)
        .unwrap_or_else(|| panic!("no v1 constant spec for {id}"))
}

/// Typed refusals (ADR 039 §4 "validity refuses, naming the quantity and
/// the range"; §8 acceptance 8 `OutOfEnvelope` / `Unsupported`).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum UniverseError {
    #[error("{what} schema {found:?} is not the supported {expected:?}")]
    Schema {
        what: &'static str,
        found: String,
        expected: &'static str,
    },
    #[error("constant {constant}: unit {found:?} is not the declared {expected:?}")]
    UnitMismatch {
        constant: &'static str,
        expected: &'static str,
        found: &'static str,
    },
    #[error("constant {constant} = {value:?} {unit} is outside the v1 validity envelope [{lo}, {hi}] ({why})")]
    OutOfEnvelope {
        constant: &'static str,
        value: f64,
        unit: &'static str,
        lo: f64,
        hi: f64,
        why: &'static str,
    },
    #[error("law {law} is not supported by this build: {why}")]
    Unsupported {
        law: &'static str,
        why: &'static str,
    },
    #[error("laws and constants are inconsistent: {0}")]
    Inconsistent(&'static str),
    #[error("pinned universe differs from what this build resolves at {path}: pinned {pinned}, resolved {resolved}")]
    PinMismatch {
        path: String,
        pinned: String,
        resolved: String,
    },
    #[error("universe serialization failed: {0}")]
    Serialization(String),
}

fn check(id: &'static str, c: &ConstantV1) -> Result<(), UniverseError> {
    let spec = constant_spec(id);
    if c.unit != spec.unit {
        return Err(UniverseError::UnitMismatch {
            constant: id,
            expected: spec.unit.as_str(),
            found: c.unit.as_str(),
        });
    }
    if !c.value.is_finite() || !(spec.lo..=spec.hi).contains(&c.value) {
        return Err(UniverseError::OutOfEnvelope {
            constant: id,
            value: c.value,
            unit: spec.unit.as_str(),
            lo: spec.lo,
            hi: spec.hi,
            why: spec.why,
        });
    }
    Ok(())
}

/// Isostasy law (ADR 039 §1 "isostasy: Airy or Pratt choice").
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum IsostasyLawV1 {
    /// LEGACY: every column against one reference column of
    /// `reference_crust_thickness_m` at 0 m; one crust density for every
    /// crust kind; sediment loads at the continental crust density
    /// (A-SEDIMENT-DENSITY); water loading against the prescribed datum.
    /// `compensation` is the existing Airy XOR Flexure dispatch (flexure
    /// needs the bounded-support solver, which the preview tier never runs).
    SingleColumnReferenceV0 {
        reference_crust_thickness_m: ConstantV1,
        compensation: Compensation,
    },
    /// Lane F / K2: continental and oceanic columns balanced against a
    /// common compensation depth (local Airy only). The depth is
    /// `fixed_compensation_depth_m`, or, in the preview chain when
    /// `derived_depth_margin_m` is set, the deepest crustal root plus that
    /// margin at every evaluation (the depth cancels; a validity bound).
    TwoColumnCommonDepthV0 {
        fixed_compensation_depth_m: ConstantV1,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        derived_depth_margin_m: Option<ConstantV1>,
    },
    /// Pratt compensation: named so the structure can carry it; not
    /// implemented, refused as `Unsupported`.
    PrattV0,
}

/// Sea-level law (ADR 039 §3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SeaLevelLawV1 {
    /// LEGACY: sea level is the prescribed constant `sea_level_datum_m`.
    PrescribedDatumV0,
    /// Ocean water volume is the input; sea level and flooded fraction are
    /// solved every isostasy evaluation by filling the water-loaded
    /// hypsometry with true cell areas (global fill: disconnected basins
    /// count as flooded). `sea_level_datum_m` is only the output datum.
    VolumePrimaryGlobalFillV0,
}

/// Initial crust seeding.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum CrustSeedingLawV1 {
    /// LEGACY: continental plates are wholly continental.
    WholePlateContinentsV0,
    /// Lane K2: one geodesic-disc body per continental plate with the
    /// `linear-rifted-margin-taper-v0` thickness profile.
    GeodesicDiscBodiesV0 {
        continental_area_fraction: ConstantV1,
        rifted_margin_width_m: ConstantV1,
    },
}

/// Oceanic thermal subsidence. The depth-age constants are fixed by the
/// law version (Stein and Stein 1992 GDH1, as cited in pf-lab A-GDH1-*).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum OceanicThermalLawV1 {
    Gdh1V0,
}

/// Erosion. Exponents m = 0.5, n = 1 and the unit of K (per year, area in
/// m^2) are fixed by the law version (ADR 039 §2 "erodibility: unit fixed
/// by the pinned erosion-law version").
///
/// SEED-OROGEN (lane OG, 2026-10-12): `Copy`/`Eq` dropped so a variant can
/// carry a `ConstantV1`; the legacy variant serializes exactly as before.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ErosionLawV1 {
    StreamPowerN1LinearDiffusionV0,
    /// Two-elevation subgrid stream power (lane OG). Each coarse cell
    /// carries a TRUNK-channel level and its cell-mean surface. The trunk
    /// obeys the same implicit stream power as the legacy law on the coarse
    /// network; the cell mean relaxes toward its trunk with the response
    /// time of the unresolved (subgrid) channel network,
    /// `tau_sub = (1/K) * integral_{x_c}^{dx} x^(-2m) dx` (self-similar
    /// network a(x) = x^2; = ln(dx/x_c)/K for m = 0.5), the travel time of an
    /// incision wave from the cell scale `dx` up to channel heads `x_c` under
    /// the SAME law. Material leaves the cell only at the cell-mean
    /// lowering rate. Reduces to the resolved law as dx -> x_c.
    StreamPowerSubgridV1 { channel_head_length_m: ConstantV1 },
    /// SEED-RELIEF-LAW (lane RL, 2026-10-12): `stream-power-subgrid-v1`
    /// generalised to a slope exponent n >= 1, applied consistently to the
    /// resolved trunk and to the subgrid network closure (implementation
    /// `pf_surface::relief`). Concavity m/n = 0.5 is fixed by the law
    /// version (m = n/2). Erodibility K_n = K * k_ref^(1 - n), with K the
    /// universe's `stream_power_erodibility_per_yr` (n = 1 units) and k_ref
    /// the `reference_steepness_m` at which the n-law equals the n = 1 law.
    /// Trunk: dc/dt = -K_n A^m ((c - c_r)/L)^n; cell mean: dz/dt =
    /// -K_n ((z - c)/ln(dx/x_c))^n (quasi-steady self-similar subgrid
    /// network). With n = 1 it is `stream-power-subgrid-v1` exactly.
    StreamPowerSubgridNonlinearV1 {
        channel_head_length_m: ConstantV1,
        slope_exponent_n: ConstantV1,
        reference_steepness_m: ConstantV1,
    },
}

impl ErosionLawV1 {
    /// True for the legacy single-elevation law.
    pub fn is_legacy(&self) -> bool {
        matches!(self, ErosionLawV1::StreamPowerN1LinearDiffusionV0)
    }
}

/// Tectonic crust thickening at convergent boundaries (SEED-OROGEN, lane
/// OG, 2026-10-12). The legacy rule is serialized by ABSENCE (the field is
/// skipped when legacy), so every existing law set, pinned universe and
/// params digest is byte-identical.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum TectonicThickeningLawV1 {
    /// LEGACY `kinematic-overlap-thickening-v0` (pf-geodynamics uplift.rs):
    /// collision stacks the loser column onto the winner up to the 70 km
    /// cap, passes the excess up to 8 cells outward, delaminates the rest;
    /// arc thickening 0.005 * v_n * dt at the overlap cell.
    KinematicOverlapStackingV0,
    /// Rate-limited belt (lane OG). Collision: the loser column consumed in
    /// an overlap (its volume, which is the shortening flux h * v_n * dt per
    /// unit boundary length) is spread over a two-sided wedge of width
    /// `W = 2 h cot(dip)`, h the mean colliding crust thickness, rate-wise
    /// over the step; thickened crust spreads laterally by lower-crustal
    /// channel flow (flux = rho_c g (1 - rho_c/rho_m) h_ch^3 / (12 mu) *
    /// grad h, h_ch = max(0, h - h_ref)); crust above
    /// `delamination_thickness_m` after spreading is delaminated (named
    /// ledger term). Arc: `arc_productivity_m` of crust per unit area of
    /// subducted plate (convergent samples only), spread over overriding-
    /// plate cells within `arc_belt_width_m`.
    RateLimitedBeltV1 {
        belt_shear_zone_dip_deg: ConstantV1,
        lower_crust_channel_viscosity_pa_s: ConstantV1,
        channel_reference_thickness_m: ConstantV1,
        delamination_thickness_m: ConstantV1,
        arc_productivity_m: ConstantV1,
        arc_belt_width_m: ConstantV1,
    },
}

impl TectonicThickeningLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, TectonicThickeningLawV1::KinematicOverlapStackingV0)
    }
    fn legacy() -> Self {
        TectonicThickeningLawV1::KinematicOverlapStackingV0
    }
}

/// How plate motion carries the per-cell crust state (SEED-CARRY-CC, lane
/// CC, 2026-10-12). The legacy carry is serialized by ABSENCE (the field is
/// skipped when legacy), so every existing law set, pinned universe and
/// params digest is byte-identical.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum PlateCarryLawV1 {
    /// LEGACY `semi-lagrangian-nearest-material-v0` (A-ADVECTION,
    /// pf-geodynamics kinematic.rs): each plate claims a cell when the
    /// globally nearest material parcel to the back-traced centre is its
    /// own; not conservative (residual reported).
    SemiLagrangianNearestMaterialV0,
    /// `front-rescue-carry-v1` (pf-geodynamics carry.rs): legacy claims,
    /// plus a rescue claim for every continental parcel the legacy rule
    /// would carry nowhere (moved, never copied, into the one cell whose
    /// back-traced point lies in its cell-shaped patch); continental losers
    /// of continent-continent overlaps that are carried nowhere are consumed
    /// once and exactly their volume and sediment are the collision supply
    /// of `rate-limited-belt-v1` (replacing its kinematic-flux supply). No
    /// constant. Requires `tectonic_thickening = rate-limited-belt-v1`.
    FrontRescueCarryV1,
    /// `continental-parcel-carry-v1` (SEED-CARRY-CONSERVE, lane CV,
    /// 2026-10-13; pf-geodynamics conserve.rs): continental crust is a set
    /// of Lagrangian parcels carrying area, thickness, sediment and plate,
    /// rotated rigidly every step; the grid is their projection (a cell is
    /// continental where the parcels' tent-kernel area covers half of it).
    /// Oceanic crust keeps the legacy Eulerian carry. Every change a step
    /// makes to a continental cell is written back to its parcels in
    /// proportion to their deposit (conservative). Area changes only by
    /// rift births (new parcels) and collision consumption (two plates each
    /// covering half a cell with a total above `1 +
    /// parcel_overlap_cell_fraction`: the smaller depositors' parcels there
    /// become the belt's collision supply). Plate-id pieces smaller than
    /// `plate_fragment_min_area_fraction` of the sphere fold into the plate
    /// they touch most. Requires `tectonic_thickening =
    /// rate-limited-belt-v1`.
    ///
    /// HYPSO-LEVERS (2026-10-21, opt-in
    /// `write_back_residual_iterations`; absent = the v1 write-back, skipped
    /// on write so every existing law set keeps its bytes): after the
    /// write-back, N passes re-project the parcels at their current positions
    /// with the same tent-kernel weights and hand the residual (end-of-step
    /// cell thickness minus that re-projection, its area-weighted mean
    /// removed so the parcel volume is unchanged) back by deposit. Why: the
    /// projection then write-back pair is not an identity (each step's
    /// changes come back smoothed over about two kernel radii), which
    /// removed 0.26 Mkm3/My of crust excess above 36 km (seed 42 climate
    /// recipe, 40-400 My; 0.45 %/My of that stock), independently of the
    /// thickness-weighted removal rule (the unweighted write-back lost the
    /// same). 3 passes (negative residuals handed by deposit x thickness)
    /// left -0.05 Mkm3/My. Measured effect (lane HYPSO-LEVERS, one run per
    /// seed, controlled windows before plate ids diverge): land > 1 km
    /// +0.21 / +0.11 / +0.21 points of land (seeds 42 / 240907 / 7), erosion
    /// +6 to +7 %, lakes mixed. Small; not a default.
    ContinentalParcelCarryV1 {
        parcel_overlap_cell_fraction: ConstantV1,
        plate_fragment_min_area_fraction: ConstantV1,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        write_back_residual_iterations: Option<ConstantV1>,
    },
}

impl PlateCarryLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, PlateCarryLawV1::SemiLagrangianNearestMaterialV0)
    }
    fn legacy() -> Self {
        PlateCarryLawV1::SemiLagrangianNearestMaterialV0
    }
}

/// What happens to the sediment of a continental collision loser under
/// `rate-limited-belt-v1` (SEDIMENT-TOWERS, lane ST, 2026-10-12;
/// implementation `pf_geodynamics::orogen::CollisionSedimentLaw`). The
/// legacy entry is serialized by ABSENCE (every existing law set, pinned
/// universe and params digest is byte-identical).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum CollisionSedimentLawV1 {
    /// LEGACY `surface-layer-v0`: the loser's sediment is added to the belt
    /// cells' surface sediment layer, outside the crust budget (no
    /// spreading, no delamination backstop); it piles up into sediment
    /// towers tens of km thick on suture cells (lane ST measurement).
    SurfaceLayerV0,
    /// `accreted-to-belt-crust-v1`: sediment caught in a collision is
    /// accreted into the orogenic wedge. A consumed loser's sediment
    /// (`front-rescue-carry-v1`) joins the belt's crust supply (placed,
    /// spread and delamination-capped like crust); a loser left to the
    /// remap keeps its sediment, as it keeps its crust (no copy). No
    /// constant. Requires `tectonic_thickening = rate-limited-belt-v1`.
    AccretedToBeltCrustV1,
}

impl CollisionSedimentLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, CollisionSedimentLawV1::SurfaceLayerV0)
    }
    fn legacy() -> Self {
        CollisionSedimentLawV1::SurfaceLayerV0
    }
}

/// Hotspot volcanism (SEED-EXTRAS, lane EX, 2026-10-12; implementation
/// `pf_geodynamics::hotspots::plume_pass`). The legacy entry (no hotspot
/// volcanism, assumption A-NO-HOTSPOTS) is serialized by ABSENCE.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum HotspotVolcanismLawV1 {
    /// LEGACY: fixed mantle hotspots are not applied at this tier.
    NoneV0,
    /// `fixed-plume-chains-v1`: `plume_count` plumes fixed in the mantle
    /// (datum) frame, positions from the K6 hotspot seed stream, onsets
    /// uniform in [-`plume_onset_lead_my`, duration). At onset a head adds
    /// `plume_head_volume_m3` of crust over a disc of `plume_head_radius_m`
    /// (any plate). From onset on, the tail adds `plume_tail_flux_m3_per_my`
    /// along the step's age-progressive track (the material that was over the
    /// plume at time t is placed where the overlying plate has carried it by
    /// the step end), each sample over a disc of `plume_edifice_radius_m` on
    /// that plate. Oceanic plates get seamount chains, continental plates
    /// volcanic tracks and plateaus. Volume is a named ledger term.
    FixedPlumeChainsV1 {
        plume_count: ConstantV1,
        plume_onset_lead_my: ConstantV1,
        plume_head_volume_m3: ConstantV1,
        plume_head_radius_m: ConstantV1,
        plume_tail_flux_m3_per_my: ConstantV1,
        plume_edifice_radius_m: ConstantV1,
    },
}

impl HotspotVolcanismLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, HotspotVolcanismLawV1::NoneV0)
    }
    fn legacy() -> Self {
        HotspotVolcanismLawV1::NoneV0
    }
}

/// What fills a gap that opens inside continental crust (SEED-EXTRAS, lane
/// EX, 2026-10-12; implementation `pf_geodynamics::rift`). Serialized by
/// ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ContinentalRiftLawV1 {
    /// LEGACY: every ridge birth is new oceanic crust, also between two
    /// diverging continents (instant breakup).
    OceanicBirthV0,
    /// `stretched-rift-v1`: a ridge birth whose nearest older cells of two
    /// diverging plates are continental is a rift birth. The continental
    /// cells of those plates within max(`rift_half_width_m`, 1.5 cell
    /// spacings) stretch (triangular weights, uniform-stretch share) and
    /// exactly their lost volume fills the birth as thin continental crust:
    /// a rift valley that deepens as it thins. A rift birth that would be
    /// thinner than `breakup_crust_thickness_m` is an oceanic birth instead
    /// (breakup), so the continent splits with thinned conjugate margins.
    StretchedRiftV1 {
        rift_half_width_m: ConstantV1,
        breakup_crust_thickness_m: ConstantV1,
    },
}

impl ContinentalRiftLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, ContinentalRiftLawV1::OceanicBirthV0)
    }
    fn legacy() -> Self {
        ContinentalRiftLawV1::OceanicBirthV0
    }
}

/// Changes to the plate set during a run (SEED-EXTRAS, lane EX,
/// 2026-10-12; implementation `pf_geodynamics::rift::split_plate`).
/// Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum PlateReorganizationLawV1 {
    /// LEGACY: the seeded plates and Euler poles never change.
    FixedPlatesV0,
    /// `plume-head-continental-split-v1`: when a plume head arrives under
    /// continental crust, the plate under it splits along a great circle
    /// through the plume (seeded azimuth). The side the cut normal points to
    /// becomes a new plate (lineage `parent_plate`) whose angular velocity is
    /// the parent's plus `rift_opening_rate_deg_per_my` about the axis on
    /// the cut 90 degrees from the plume, so the two halves open at the
    /// plume. Requires `fixed-plume-chains-v1` and `stretched-rift-v1`.
    PlumeHeadContinentalSplitV1 { rift_opening_rate_deg_per_my: ConstantV1 },
    /// `plate-reorg-scripted-v1` (SEED-PLATE-REORG, lane PR, 2026-10-13;
    /// implementation `pf_geodynamics::reorg`): a scripted Wilson cycle,
    /// not a force balance. Groups are suture groups of
    /// `contact-length-pole-share-v1`. (1) At seeded times (gaps uniform in
    /// [reorg_interval_min_my, reorg_interval_max_my], one stream, so the
    /// schedule does not depend on the run length) every group's pole is
    /// re-drawn in ascending root order: axis = (1 - slab_pull_bias) random
    /// + slab_pull_bias x the group's slab torque; rate = rate scale x
    /// (floor + (1 - floor) min(1, slab share / slab_circumference_reference))
    /// x (0.75 + 0.5 u); relations younger than rift_persistence_my keep
    /// their relative rotation. (2) A group holding supercontinent_share of
    /// the continental crust for supercontinent_tenure_my breaks up along a
    /// min cut over its cells whose capacity is edge length x strength
    /// (weak zones: crust far from the group's median thickness, old
    /// shortening, old sutures between member plates, down to
    /// weak_zone_strength_ratio); a piece still above breakup_piece_max_share
    /// breaks again breakup_cascade_my later; the new group opens at
    /// rift_opening_rate_deg_per_my. (3) Plume heads under groups holding
    /// plume_split_min_share split them along the same weak-zone cut
    /// (weakened within the plume head radius). (4) While fewer than
    /// plate_count_min groups exceed plate_count_min_area_fraction of the
    /// sphere, the group with the most ocean splits along its passive margin
    /// (or a random axis) and its ocean side converges at
    /// subduction_initiation_rate_deg_per_my; above plate_count_max the
    /// smallest group without continent folds into its neighbour; cut
    /// leftovers fold or become plates. (5) Plume onsets are a renewal
    /// process (mean gap plume_head_interval_my) with tails lasting
    /// plume_tail_lifetime_my, independent of the run length. Requires
    /// fixed-plume-chains-v1, stretched-rift-v1, continental-parcel-carry-v1
    /// and contact-length-pole-share-v1.
    PlateReorgScriptedV1 {
        rift_opening_rate_deg_per_my: ConstantV1,
        reorg_interval_min_my: ConstantV1,
        reorg_interval_max_my: ConstantV1,
        slab_pull_bias: ConstantV1,
        slab_circumference_reference: ConstantV1,
        reorg_rate_floor_fraction: ConstantV1,
        supercontinent_share: ConstantV1,
        supercontinent_tenure_my: ConstantV1,
        breakup_piece_max_share: ConstantV1,
        breakup_cascade_my: ConstantV1,
        rift_persistence_my: ConstantV1,
        weak_zone_strength_ratio: ConstantV1,
        plate_count_min: ConstantV1,
        plate_count_max: ConstantV1,
        plate_count_min_area_fraction: ConstantV1,
        subduction_initiation_rate_deg_per_my: ConstantV1,
        plume_head_interval_my: ConstantV1,
        plume_tail_lifetime_my: ConstantV1,
        plume_split_min_share: ConstantV1,
    },
    /// `plate-reorg-scripted-v2` (SEED-BREAKUP-CUT-V2, lane BC, 2026-10-13;
    /// implementation `pf_geodynamics::reorg` with `ReorgLaw::ScriptedV2`;
    /// pre-registered in `unpublished-design-record`):
    /// `plate-reorg-scripted-v1` with the same constants and a new breakup
    /// and plume-split cut. Strength comes from cumulative shortening, old
    /// sutures (continental edges between member plates) and the plume head
    /// disc only (no thickness term); the cut runs over the group's
    /// continental cells plus a one-cell halo at least as strong as craton,
    /// and the rest of the group takes the side of the nearest domain cell,
    /// so the new boundary leaves the coast at the two exit points and
    /// continues through the ocean; cuts leaving less than 15% of the
    /// group's continent on either side are rejected (breakups retry the
    /// axis); the max flow is on integer capacities (canonical minimal
    /// source side). Plates emptied by folds retire as `Folded`. Same
    /// requirements as v1.
    PlateReorgScriptedV2 {
        rift_opening_rate_deg_per_my: ConstantV1,
        reorg_interval_min_my: ConstantV1,
        reorg_interval_max_my: ConstantV1,
        slab_pull_bias: ConstantV1,
        slab_circumference_reference: ConstantV1,
        reorg_rate_floor_fraction: ConstantV1,
        supercontinent_share: ConstantV1,
        supercontinent_tenure_my: ConstantV1,
        breakup_piece_max_share: ConstantV1,
        breakup_cascade_my: ConstantV1,
        rift_persistence_my: ConstantV1,
        weak_zone_strength_ratio: ConstantV1,
        plate_count_min: ConstantV1,
        plate_count_max: ConstantV1,
        plate_count_min_area_fraction: ConstantV1,
        subduction_initiation_rate_deg_per_my: ConstantV1,
        plume_head_interval_my: ConstantV1,
        plume_tail_lifetime_my: ConstantV1,
        plume_split_min_share: ConstantV1,
    },
    /// `plate-reorg-scripted-v3` (SEED-WORLD-COMPOSE, lane CP, 2026-10-13;
    /// implementation `pf_geodynamics::reorg` with `ReorgLaw::ScriptedV3`;
    /// pre-registered in `unpublished-design-record`):
    /// `plate-reorg-scripted-v2` with subduction initiation on the v2 cut
    /// machinery (integer max flow) and driven by ocean-crust age. At every
    /// scheduled event, besides the v1/v2 "too few groups" trigger, a group
    /// whose ocean older than `subduction_initiation_age_my` covers at least
    /// the counted-plate area starts a trench (oldest-crust group first,
    /// while fewer than `plate_count_max` groups): the cut separates the
    /// group's continent (and the ocean behind it) from the far ocean in
    /// the direction of the old crust; ocean cells weaken linearly with age
    /// from 100 to 180 My (to `weak_zone_strength_ratio`), the one-cell
    /// coastal halo is at least craton strength and continental cells are
    /// 10x stronger, so trenches open in old ocean on the ocean side of
    /// passive margins and in enclosed basins, not along whole coasts; cuts
    /// that split a continent (5-95% of it on one side) are rejected
    /// (axis retried). Same requirements as v1.
    PlateReorgScriptedV3 {
        rift_opening_rate_deg_per_my: ConstantV1,
        reorg_interval_min_my: ConstantV1,
        reorg_interval_max_my: ConstantV1,
        slab_pull_bias: ConstantV1,
        slab_circumference_reference: ConstantV1,
        reorg_rate_floor_fraction: ConstantV1,
        supercontinent_share: ConstantV1,
        supercontinent_tenure_my: ConstantV1,
        breakup_piece_max_share: ConstantV1,
        breakup_cascade_my: ConstantV1,
        rift_persistence_my: ConstantV1,
        weak_zone_strength_ratio: ConstantV1,
        plate_count_min: ConstantV1,
        plate_count_max: ConstantV1,
        plate_count_min_area_fraction: ConstantV1,
        subduction_initiation_rate_deg_per_my: ConstantV1,
        plume_head_interval_my: ConstantV1,
        plume_tail_lifetime_my: ConstantV1,
        plume_split_min_share: ConstantV1,
        subduction_initiation_age_my: ConstantV1,
    },
    /// `plate-reorg-scripted-v4` (SEED-TRANSFORMS, lane TF, 2026-10-13;
    /// implementation `pf_geodynamics::reorg` with `ReorgLaw::ScriptedV4`;
    /// pre-registered in `unpublished-design-record`):
    /// `plate-reorg-scripted-v3` with two limits on the AGE-triggered
    /// initiation only. A candidate cut is rejected (next axis tried) when
    /// less than `initiation_old_slab_fraction_min` of the new boundary's
    /// length (4-neighbour edges between the two sides inside the group)
    /// runs where the older flanking cell is ocean at least
    /// `initiation_slab_age_min_my` old, or when more than half of it
    /// (fixed by version) runs within `initiation_trench_spacing_m` of a
    /// subduction-class cell or of a boundary started earlier in the same
    /// pass. Same requirements as v1.
    PlateReorgScriptedV4 {
        rift_opening_rate_deg_per_my: ConstantV1,
        reorg_interval_min_my: ConstantV1,
        reorg_interval_max_my: ConstantV1,
        slab_pull_bias: ConstantV1,
        slab_circumference_reference: ConstantV1,
        reorg_rate_floor_fraction: ConstantV1,
        supercontinent_share: ConstantV1,
        supercontinent_tenure_my: ConstantV1,
        breakup_piece_max_share: ConstantV1,
        breakup_cascade_my: ConstantV1,
        rift_persistence_my: ConstantV1,
        weak_zone_strength_ratio: ConstantV1,
        plate_count_min: ConstantV1,
        plate_count_max: ConstantV1,
        plate_count_min_area_fraction: ConstantV1,
        subduction_initiation_rate_deg_per_my: ConstantV1,
        plume_head_interval_my: ConstantV1,
        plume_tail_lifetime_my: ConstantV1,
        plume_split_min_share: ConstantV1,
        subduction_initiation_age_my: ConstantV1,
        initiation_slab_age_min_my: ConstantV1,
        initiation_old_slab_fraction_min: ConstantV1,
        initiation_trench_spacing_m: ConstantV1,
    },
    /// `plate-reorg-scripted-v5` (SEED-OLD-OCEAN (lane OR, 2026-10-13);
    /// implementation `pf_geodynamics::reorg` with `ReorgLaw::ScriptedV5`;
    /// pre-registered in `unpublished-design-record`):
    /// `plate-reorg-scripted-v4` plus passive-margin subduction initiation
    /// driven by ocean age (lane VL defect 3). At each scheduled event,
    /// while fewer than the maximum groups, the passive-margin ocean
    /// cell at least `margin_initiation_age_my` old with the most such ocean within 1000 km, farther than
    /// `margin_trench_spacing_m` from existing subduction seeds a trench
    /// along at most `margin_segment_length_m` of that margin; the margin's
    /// frontal ocean (nearest to that stretch of coast) becomes a new plate
    /// converging on the continent (at most 3 per event, one per group); its re-solves keep, per group, the candidate axis (8, two sweeps) closing the least ridge or young-young boundary (OR-WORDING-R1).
    /// Same requirements as v1.
    PlateReorgScriptedV5 {
        rift_opening_rate_deg_per_my: ConstantV1,
        reorg_interval_min_my: ConstantV1,
        reorg_interval_max_my: ConstantV1,
        slab_pull_bias: ConstantV1,
        slab_circumference_reference: ConstantV1,
        reorg_rate_floor_fraction: ConstantV1,
        supercontinent_share: ConstantV1,
        supercontinent_tenure_my: ConstantV1,
        breakup_piece_max_share: ConstantV1,
        breakup_cascade_my: ConstantV1,
        rift_persistence_my: ConstantV1,
        weak_zone_strength_ratio: ConstantV1,
        plate_count_min: ConstantV1,
        plate_count_max: ConstantV1,
        plate_count_min_area_fraction: ConstantV1,
        subduction_initiation_rate_deg_per_my: ConstantV1,
        plume_head_interval_my: ConstantV1,
        plume_tail_lifetime_my: ConstantV1,
        plume_split_min_share: ConstantV1,
        subduction_initiation_age_my: ConstantV1,
        initiation_slab_age_min_my: ConstantV1,
        initiation_old_slab_fraction_min: ConstantV1,
        initiation_trench_spacing_m: ConstantV1,
        margin_initiation_age_my: ConstantV1,
        margin_segment_length_m: ConstantV1,
        margin_trench_spacing_m: ConstantV1,
    },
}

impl PlateReorganizationLawV1 {
    /// SEED-WORLD-COMPOSE (lane CP): the v3 law with a v2 law's constants
    /// plus the initiation age; any other law unchanged.
    pub fn into_scripted_v3(self, subduction_initiation_age_my: ConstantV1) -> Self {
        match self {
            PlateReorganizationLawV1::PlateReorgScriptedV2 {
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
            } => PlateReorganizationLawV1::PlateReorgScriptedV3 {
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
                subduction_initiation_age_my,
            },
            other => other,
        }
    }
    /// SEED-TRANSFORMS (lane TF): the v4 law from a v3 law's constants plus
    /// the three initiation limits; any other law unchanged.
    pub fn into_scripted_v4(
        self,
        initiation_slab_age_min_my: ConstantV1,
        initiation_old_slab_fraction_min: ConstantV1,
        initiation_trench_spacing_m: ConstantV1,
    ) -> Self {
        match self {
            PlateReorganizationLawV1::PlateReorgScriptedV3 {
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
                subduction_initiation_age_my,
            } => PlateReorganizationLawV1::PlateReorgScriptedV4 {
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
                subduction_initiation_age_my,
                initiation_slab_age_min_my,
                initiation_old_slab_fraction_min,
                initiation_trench_spacing_m,
            },
            other => other,
        }
    }
    /// SEED-OLD-OCEAN (lane OR): the v5 law from a v4 law's constants plus
    /// the three margin-initiation constants; any other law unchanged.
    pub fn into_scripted_v5(
        self,
        margin_initiation_age_my: ConstantV1,
        margin_segment_length_m: ConstantV1,
        margin_trench_spacing_m: ConstantV1,
    ) -> Self {
        match self {
            PlateReorganizationLawV1::PlateReorgScriptedV4 {
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
                subduction_initiation_age_my,
                initiation_slab_age_min_my,
                initiation_old_slab_fraction_min,
                initiation_trench_spacing_m,
            } => PlateReorganizationLawV1::PlateReorgScriptedV5 {
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
                subduction_initiation_age_my,
                initiation_slab_age_min_my,
                initiation_old_slab_fraction_min,
                initiation_trench_spacing_m,
                margin_initiation_age_my,
                margin_segment_length_m,
                margin_trench_spacing_m,
            },
            other => other,
        }
    }
    /// SEED-WORLD-COMPOSE (lane CP): v2 and v3 viewed as v1 (same shared
    /// constants) for the shared checks; any other law unchanged.
    pub fn shared_view(&self) -> Self {
        match self.clone() {
            PlateReorganizationLawV1::PlateReorgScriptedV2 {
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
            | PlateReorganizationLawV1::PlateReorgScriptedV5 {
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
                ..
            }
            | PlateReorganizationLawV1::PlateReorgScriptedV4 {
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
                ..
            }
            | PlateReorganizationLawV1::PlateReorgScriptedV3 {
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
                ..
            } => PlateReorganizationLawV1::PlateReorgScriptedV1 {
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
            },
            other => other,
        }
    }
    /// SEED-BREAKUP-CUT-V2 (lane BC): the v2 law with a v1 law's constants
    /// (same ids, values and labels); any other law unchanged.
    pub fn into_scripted_v2(self) -> Self {
        match self {
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
            } => PlateReorganizationLawV1::PlateReorgScriptedV2 {
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
            },
            other => other,
        }
    }
    pub fn is_legacy(&self) -> bool {
        matches!(self, PlateReorganizationLawV1::FixedPlatesV0)
    }
    fn legacy() -> Self {
        PlateReorganizationLawV1::FixedPlatesV0
    }
}

/// Thermal and flexural response of stretched continental crust
/// (SEED-RIFT-SHOULDERS, lane RS, 2026-10-13; implementation
/// `pf_geodynamics::rift_flank`). Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum RiftFlankLawV1 {
    /// LEGACY: stretched continental columns carry no thermal term and
    /// their flanks no shoulder (local Airy on the crust only).
    NoneV0,
    /// `shoulder-uplift-thermal-v1`: (b) a stretched continental column
    /// (stretching factor beta from `stretched-rift-v1`) carries the
    /// uniform-stretching thermal uplift (1 - 1/beta) x the oceanic GDH1
    /// thermal elevation at its thermal age (McKenzie 1978 class; the
    /// sqrt-age half-space branch early, exponential late), so the rift
    /// floor subsides thermally with time; repeated stretching adds to the
    /// anomaly with a continuous equivalent age. (a) Every continental
    /// column gets a shoulder from each stretched column within reach:
    /// `rift_shoulder_ratio` x that column's thermal uplift x the
    /// cell-mean broken-plate profile exp(-x/alpha) cos(x/alpha)
    /// (x from the source cell's edge, zero beyond pi alpha / 2), applied
    /// where it exceeds the column's own thermal uplift. Elevation only:
    /// no crust volume is created or moved. Requires `stretched-rift-v1`
    /// and the two-column isostasy law.
    ShoulderUpliftThermalV1 {
        rift_flank_flexural_parameter_m: ConstantV1,
        rift_shoulder_ratio: ConstantV1,
    },
}

impl RiftFlankLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, RiftFlankLawV1::NoneV0)
    }
    fn legacy() -> Self {
        RiftFlankLawV1::NoneV0
    }
}

/// Threshold hillslopes at the cell scale (SEED-UPLANDS, lane UP,
/// 2026-10-13; implementation `pf_surface::uplands`). Serialized by
/// ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum HillslopeThresholdLawV1 {
    /// LEGACY: no slope threshold beyond the erosion law itself.
    NoneV0,
    /// `critical-slope-envelope-v1`: cell means may not stand above the
    /// envelope z_env(i) = min_j (b_j + S_c path(i, j)) (b = land mean or
    /// water fill level; 8-neighbour Dijkstra). A land cell above it erodes
    /// at (z - z_env) / tau, tau = cell size / `knickpoint_celerity_m_per_my`
    /// (backward Euler, envelope held for the step); the cell-mean target
    /// is max(subgrid stream-power target, this), realised through the
    /// shared erosion stage. Isostasy returns most of each eroded metre,
    /// so ranges taper through foothill cells and lose height gradually.
    /// Requires a subgrid erosion law (`stream-power-subgrid-v1` or
    /// `stream-power-subgrid-nonlinear-v1`).
    CriticalSlopeEnvelopeV1 {
        hillslope_critical_slope: ConstantV1,
        knickpoint_celerity_m_per_my: ConstantV1,
    },
}

impl HillslopeThresholdLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, HillslopeThresholdLawV1::NoneV0)
    }
    fn legacy() -> Self {
        HillslopeThresholdLawV1::NoneV0
    }
}

/// Fluvial incision threshold (lane INCISION-THRESHOLD, 2026-10-20;
/// implementation `pf_surface::erosion::incision_threshold`; wiring merged
/// by lane CLIMATE-GENESIS-2). Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum FluvialIncisionThresholdLawV1 {
    /// LEGACY: no threshold beyond the erosion law itself.
    NoneV0,
    /// `stream-power-excess-steepness-threshold-v1`: E = K_n max(0, k_s^n -
    /// k_c^n) on the trunk and on the subgrid cell-mean closure of
    /// `stream-power-subgrid-nonlinear-v1` (k_s = S A^(1/2), theta = 0.5
    /// units); backward Euler with the threshold term constant over the
    /// step; never raises a level. Requires that erosion law.
    ExcessSteepnessThresholdV1 { critical_channel_steepness_m: ConstantV1 },
    /// LAKE-HOLD-2 (lane LAKE-HOLD-2, 2026-10-21):
    /// `stream-power-excess-steepness-threshold-trunk-only-v1`: the same
    /// threshold on the trunk channel only; the subgrid cell-mean closure
    /// runs at k_c = 0 (the detachment threshold is a property of bedrock
    /// channels, not of the soil-mantled subgrid hillslopes it otherwise
    /// freezes below Lambda k_c of mean relief). Same constant, same
    /// requirements.
    ExcessSteepnessThresholdTrunkOnlyV1 { critical_channel_steepness_m: ConstantV1 },
}

impl FluvialIncisionThresholdLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, FluvialIncisionThresholdLawV1::NoneV0)
    }
    fn legacy() -> Self {
        FluvialIncisionThresholdLawV1::NoneV0
    }
}

/// ALLUVIAL-TRUNK (lane ALLUVIAL-TRUNK, 2026-10-21; implementation
/// `pf_surface::erosion::incision_threshold` and pf-surface closed_basin.rs
/// rule 10): where the fluvial incision threshold is released on alluvial
/// channels. Serialized by ABSENCE when legacy, so every existing law set
/// keeps its bytes. Requires a fluvial incision threshold (its k_c) and the
/// sediment retention law (its graded alluvial slope S_a). No constants of
/// its own.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum FluvialAlluvialChannelLawV1 {
    /// LEGACY: the threshold acts on every channel it is selected for.
    NoneV0,
    /// `stream-power-threshold-free-alluvial-channels-v1`: on the subgrid
    /// trunk and on the lake outlet channel (rules 5-9), the threshold term
    /// is dropped where the channel is alluvial: its bedrock threshold slope
    /// k_c / A^(m/n) is below the retention law's graded alluvial slope S_a
    /// (drainage area above (k_c / S_a)^(n/m)). Such reaches are graded by
    /// the retention law's alluvial grade (rule 2a), not held at k_c.
    ThresholdFreeAlluvialChannelsV1,
}

impl FluvialAlluvialChannelLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, FluvialAlluvialChannelLawV1::NoneV0)
    }
    fn legacy() -> Self {
        FluvialAlluvialChannelLawV1::NoneV0
    }
}

/// Glacial erosion of T0 genesis (lane GLACIAL-T0, 2026-10-21;
/// implementation `pf_surface::erosion::glacial`). Serialized by ABSENCE
/// when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum GlacialErosionLawV1 {
    /// LEGACY: no glacial erosion.
    NoneV0,
    /// `glacial-stream-power-blend-v1`: Hergarten 2021 (ESurf 9:937) blended
    /// glacial + subglacial-fluvial stream power E = K_eff A^m S^n,
    /// K_eff = gamma^m K_g + (1 - gamma)^m K_f, gamma = ice share of the
    /// precipitation flux, as a per-cell erodibility factor on the nonlinear
    /// subgrid trunk and cell-mean closure; ice from an Ohmura 1992 ELA on
    /// the genesis climate's temperature and precipitation; climate states
    /// as time-averaged duties (a coldhouse with sinusoidal glacial cycles
    /// plus warm states). Requires `stream-power-subgrid-nonlinear-v1` and
    /// `epoch-climate-runoff-v1`.
    GlacialStreamPowerBlendV1 {
        glacial_to_fluvial_erodibility_ratio: ConstantV1,
        full_ice_height_above_ela_m: ConstantV1,
        coldhouse_duty: ConstantV1,
        glacial_max_cooling_k: ConstantV1,
        glacial_max_precipitation_factor: ConstantV1,
        /// Warm states: (duty, warming K) pairs, same length.
        warm_state_duties: Vec<ConstantV1>,
        warm_state_warmings_k: Vec<ConstantV1>,
    },
}

impl GlacialErosionLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, GlacialErosionLawV1::NoneV0)
    }
    fn legacy() -> Self {
        GlacialErosionLawV1::NoneV0
    }
}

/// How the marine remainder (river-mouth volume left after the bounded
/// local fill to sea level) is laid as the declared global layer (lane
/// INCISION-THRESHOLD part 1, 2026-10-20; `pf_surface::erosion::
/// OceanRemainderLayer`). Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum MarineRemainderLayerLawV1 {
    /// LEGACY `per-mouth-global-layer-v0`: each saturated mouth raises every
    /// ocean cell (cost mouths x ocean cells; lifts sea-level cells above it).
    PerMouthGlobalLayerV0,
    /// `once-per-step-capped-global-layer-v1`: the step's remainders laid
    /// once over ocean cells below sea level, each capped at sea level;
    /// local fills land exactly on their water level.
    OncePerStepCappedGlobalLayerV1,
}

impl MarineRemainderLayerLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, MarineRemainderLayerLawV1::PerMouthGlobalLayerV0)
    }
    fn legacy() -> Self {
        MarineRemainderLayerLawV1::PerMouthGlobalLayerV0
    }
}

/// Flexure of a subducting oceanic plate (SEED-TRENCH, lane TN,
/// 2026-10-13; implementation `pf_geodynamics::trench`). Serialized by
/// ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum SubductionTrenchLawV1 {
    /// LEGACY: the downgoing plate floats at its GDH1 depth up to the
    /// boundary (no trench, no outer rise).
    NoneV0,
    /// `slab-flexure-trench-v1`: at a Subduction-class boundary the oceanic
    /// cell overridden by a continental or younger oceanic neighbour (the
    /// kinematic overlap rule) is a trench cell. Its plate (age A, normal
    /// convergence v) bends as a broken elastic plate under the end load
    /// V0 = f x slab pull, slab pull = 2 rho_m alpha_v T_m g sqrt(kappa A /
    /// pi) min(L_max, v 2.32^2 A / pi^2); Te = 2 erfinv(T_e / T_m)
    /// sqrt(kappa A), D = E Te^3 / (12 (1 - nu^2)), alpha = (4 D / ((rho_m -
    /// rho_w) g))^(1/4), w0 = V0 / ((rho_m - rho_w) g alpha). Every oceanic
    /// cell of the downgoing plate takes the cell-mean of w0 exp(-x/alpha)
    /// cos(x/alpha) (x from the trench axis, cut at 3 pi alpha / 2) of its
    /// nearest trench cell, as a negative water-loaded column term:
    /// trench at the axis, outer rise seaward. Elevation only. Requires the
    /// two-column isostasy law.
    SlabFlexureTrenchV1 {
        trench_vertical_load_fraction: ConstantV1,
        trench_slab_length_max_m: ConstantV1,
        trench_thermal_expansivity_per_k: ConstantV1,
        trench_mantle_temperature_c: ConstantV1,
        trench_elastic_isotherm_c: ConstantV1,
        trench_thermal_diffusivity_m2_per_yr: ConstantV1,
        trench_elastic_modulus_pa: ConstantV1,
        trench_poisson_ratio: ConstantV1,
    },
}

impl SubductionTrenchLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, SubductionTrenchLawV1::NoneV0)
    }
    fn legacy() -> Self {
        SubductionTrenchLawV1::NoneV0
    }
}

/// Continent-continent suturing (SEED-CARRY-CONSERVE, lane CV,
/// 2026-10-13; implementation `pf_geodynamics::conserve::suture_pass`).
/// Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum PlateSutureLawV1 {
    /// LEGACY: plates keep their poles; colliding continents converge
    /// for ever.
    NoneV0,
    /// `contact-length-pole-share-v1`: at the start of every step, when the
    /// converging continent-continent contact length between two plates
    /// (shared cell edges between their continental cells with convergent
    /// normal rate) reaches `suture_contact_length_m`, the plate group with
    /// the smaller total area takes the Euler pole of the larger
    /// (union-find; permanent), so convergence across the suture stops.
    ContactLengthPoleShareV1 { suture_contact_length_m: ConstantV1 },
    /// `contact-length-indentation-taper-v1` (SEED-RANGES, lane RG,
    /// 2026-10-13; implementation `pf_geodynamics::ranges::
    /// suture_indent_pass`): suturing exactly as
    /// `contact-length-pole-share-v1`, but the smaller group's member plates
    /// keep their relative rotation at the suture, scaled by the linear
    /// taper `1 - tau / indentation_time_my`, until that time or until the
    /// cumulative mean convergence across their contact with the rest of
    /// the group reaches `indentation_shortening_budget_m` (post-collision
    /// indentation; convergence tapers instead of stopping).
    ContactLengthIndentationTaperV1 {
        suture_contact_length_m: ConstantV1,
        indentation_time_my: ConstantV1,
        indentation_shortening_budget_m: ConstantV1,
    },
}

impl PlateSutureLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, PlateSutureLawV1::NoneV0)
    }
    fn legacy() -> Self {
        PlateSutureLawV1::NoneV0
    }
}

/// Orogenic plateau growth (SEED-RANGES, lane RG, 2026-10-13;
/// implementation `pf_geodynamics::ranges::plateau_spill`). Requires
/// `rate-limited-belt-v1`. Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum OrogenPlateauLawV1 {
    /// LEGACY: the belt's delamination backstop alone limits thickness.
    NoneV0,
    /// `plateau-cap-spill-v1`: after the belt step, continental crust above
    /// `plateau_crust_cap_m` gives its excess to the nearest same-plate
    /// continental cells below the cap (breadth-first over edge
    /// neighbours through cells already at the cap, each filled up to the
    /// cap), within `plateau_spill_reach_m` of the source; what finds no
    /// room is delaminated (counted). Oceanic crust above the cap is
    /// delaminated. Builds flat-topped plateaus that widen into the
    /// hinterland instead of one-cell spikes.
    PlateauCapSpillV1 {
        plateau_crust_cap_m: ConstantV1,
        plateau_spill_reach_m: ConstantV1,
    },
    /// `plateau-cap-spill-v2` (SEED-WORLD-COMPOSE, lane CP, 2026-10-13;
    /// implementation `pf_geodynamics::ranges::plateau_in_belt`, run inside
    /// the belt step): the step's collision supply is placed, then on every
    /// continental column whose crust plus sediment exceeds the cap the
    /// sediment is accreted into the crust (booked like
    /// accreted-to-belt-crust-v1), crust above the cap spills as in v1
    /// (same-plate breadth-first rings within the reach, remainder
    /// delaminated), then a foothill taper limits the thickness step between
    /// edge-neighbouring same-plate continental cells to
    /// `plateau_foothill_gradient` x their centre distance wherever the
    /// thicker one exceeds `plateau_taper_floor_m` (volume moved down-step,
    /// conservative, fixed sweep order), and only then does the channel flow
    /// spread the field (no source, so it never sees crust above the cap),
    /// with the polar limiter: per-cell diffusivity x min(1, (2 dx/dy)^2)
    /// (lane CP PREREG revision 1; keeps converging-meridian cells from
    /// needing 1e5+ explicit substeps).
    PlateauCapSpillV2 {
        plateau_crust_cap_m: ConstantV1,
        plateau_spill_reach_m: ConstantV1,
        plateau_foothill_gradient: ConstantV1,
        plateau_taper_floor_m: ConstantV1,
    },
}

impl OrogenPlateauLawV1 {
    /// SEED-WORLD-COMPOSE (lane CP): v2 viewed as v1 for the shared checks.
    pub fn shared_view(&self) -> Self {
        match self.clone() {
            OrogenPlateauLawV1::PlateauCapSpillV2 { plateau_crust_cap_m, plateau_spill_reach_m, .. } => {
                OrogenPlateauLawV1::PlateauCapSpillV1 { plateau_crust_cap_m, plateau_spill_reach_m }
            }
            other => other,
        }
    }
}

impl OrogenPlateauLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, OrogenPlateauLawV1::NoneV0)
    }
    fn legacy() -> Self {
        OrogenPlateauLawV1::NoneV0
    }
}

/// Subduction and plume magmatism (SEED-MAGMATISM, lane MG, 2026-10-13;
/// implementation `pf_geodynamics::{arcs, backarc, volcanoes}`, LLSVP plume
/// sources in `pf_geodynamics::hotspots`, glue `pf_lab::magmatism`;
/// pre-registered in `unpublished-design-record`). Requires
/// the two-column isostasy, `rate-limited-belt-v1` and the parcel carry or
/// a suture law (its run state lives with theirs). Serialized by ABSENCE
/// when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum MagmatismLawV1 {
    /// LEGACY: no slab-depth arcs, back-arcs or LLSVP plume sources.
    NoneV0,
    /// `slab-arc-backarc-llsvp-v1`: after the tectonic advance, arc crust is
    /// added on the overriding plate at d = slab depth / tan(dip) from every
    /// trench (dip from the slab age; flux = productivity x clamped
    /// convergence rate x trench length); oceanic columns built mostly of
    /// arc crust become continental at the threshold thickness; old slabs
    /// under a non-advancing upper plate open back-arcs (oceanic crust
    /// renewed, continental crust stretched toward the arc); a decaying
    /// arc heat lifts the overriding plate; plume sources sit at the edges
    /// of two antipodal LLSVPs and drift; vent markers give derived volcano
    /// points and seven additive resource families.
    SlabArcBackarcLlsvpV1 {
        slab_depth_below_arc_m: ConstantV1,
        slab_dip_young_deg: ConstantV1,
        slab_dip_old_deg: ConstantV1,
        slab_dip_age_my: ConstantV1,
        arc_magma_productivity_m: ConstantV1,
        arc_rate_floor_m_per_my: ConstantV1,
        arc_rate_cap_m_per_my: ConstantV1,
        arc_half_width_m: ConstantV1,
        arc_continental_thickness_m: ConstantV1,
        arc_thermal_uplift_m: ConstantV1,
        arc_thermal_half_width_m: ConstantV1,
        arc_thermal_decay_my: ConstantV1,
        backarc_slab_age_min_my: ConstantV1,
        backarc_overriding_advance_max_m_per_my: ConstantV1,
        backarc_extension_rate_m_per_my: ConstantV1,
        backarc_axis_offset_m: ConstantV1,
        backarc_half_width_m: ConstantV1,
        backarc_min_crust_m: ConstantV1,
        llsvp_edge_deg: ConstantV1,
        llsvp_edge_sigma_deg: ConstantV1,
        llsvp_antipode_jitter_deg: ConstantV1,
        plume_background_fraction: ConstantV1,
        plume_drift_rate_m_per_my: ConstantV1,
        vent_active_window_my: ConstantV1,
        volcano_spacing_m: ConstantV1,
        arc_crust_cap_m: ConstantV1,
    },
    /// `slab-arc-backarc-llsvp-v2` (SEED-WORLD-COMPOSE, lane CP,
    /// 2026-10-13): v1 with (a) the arc heat applied as a dynamic elevation
    /// term AFTER the sea-level solve (it lifts columns but displaces no
    /// ocean water; the ocean volume constant is calibrated on a hypsometry
    /// that already contains arc and back-arc thermal highs), (b) back-arc
    /// extension one episode per cell: a cell is renewed or stretched only
    /// within `backarc_episode_my` of its first back-arc extension, and
    /// continental cells never below 1/`backarc_max_stretch_factor` of their
    /// thickness at the episode start, (c) the kimberlite plume-pass marker
    /// set only by an arriving plume head, within
    /// `kimberlite_plume_radius_m` of it.
    SlabArcBackarcLlsvpV2 {
        slab_depth_below_arc_m: ConstantV1,
        slab_dip_young_deg: ConstantV1,
        slab_dip_old_deg: ConstantV1,
        slab_dip_age_my: ConstantV1,
        arc_magma_productivity_m: ConstantV1,
        arc_rate_floor_m_per_my: ConstantV1,
        arc_rate_cap_m_per_my: ConstantV1,
        arc_half_width_m: ConstantV1,
        arc_continental_thickness_m: ConstantV1,
        arc_thermal_uplift_m: ConstantV1,
        arc_thermal_half_width_m: ConstantV1,
        arc_thermal_decay_my: ConstantV1,
        backarc_slab_age_min_my: ConstantV1,
        backarc_overriding_advance_max_m_per_my: ConstantV1,
        backarc_extension_rate_m_per_my: ConstantV1,
        backarc_axis_offset_m: ConstantV1,
        backarc_half_width_m: ConstantV1,
        backarc_min_crust_m: ConstantV1,
        llsvp_edge_deg: ConstantV1,
        llsvp_edge_sigma_deg: ConstantV1,
        llsvp_antipode_jitter_deg: ConstantV1,
        plume_background_fraction: ConstantV1,
        plume_drift_rate_m_per_my: ConstantV1,
        vent_active_window_my: ConstantV1,
        volcano_spacing_m: ConstantV1,
        arc_crust_cap_m: ConstantV1,
        backarc_episode_my: ConstantV1,
        backarc_max_stretch_factor: ConstantV1,
        kimberlite_plume_radius_m: ConstantV1,
    },
    /// Teaching opt-in: v2 episodes/resources with arc heat inside the
    /// water-loaded column solve, so the prescribed water inventory closes.
    SlabArcBackarcLlsvpV3 {
        slab_depth_below_arc_m: ConstantV1,
        slab_dip_young_deg: ConstantV1,
        slab_dip_old_deg: ConstantV1,
        slab_dip_age_my: ConstantV1,
        arc_magma_productivity_m: ConstantV1,
        arc_rate_floor_m_per_my: ConstantV1,
        arc_rate_cap_m_per_my: ConstantV1,
        arc_half_width_m: ConstantV1,
        arc_continental_thickness_m: ConstantV1,
        arc_thermal_uplift_m: ConstantV1,
        arc_thermal_half_width_m: ConstantV1,
        arc_thermal_decay_my: ConstantV1,
        backarc_slab_age_min_my: ConstantV1,
        backarc_overriding_advance_max_m_per_my: ConstantV1,
        backarc_extension_rate_m_per_my: ConstantV1,
        backarc_axis_offset_m: ConstantV1,
        backarc_half_width_m: ConstantV1,
        backarc_min_crust_m: ConstantV1,
        llsvp_edge_deg: ConstantV1,
        llsvp_edge_sigma_deg: ConstantV1,
        llsvp_antipode_jitter_deg: ConstantV1,
        plume_background_fraction: ConstantV1,
        plume_drift_rate_m_per_my: ConstantV1,
        vent_active_window_my: ConstantV1,
        volcano_spacing_m: ConstantV1,
        arc_crust_cap_m: ConstantV1,
        backarc_episode_my: ConstantV1,
        backarc_max_stretch_factor: ConstantV1,
        kimberlite_plume_radius_m: ConstantV1,
    },
}

impl MagmatismLawV1 {
    /// Preserve v2 parameters while selecting the volume-conserving v3 policy.
    pub fn into_v3(self) -> Self {
        match self {
            MagmatismLawV1::SlabArcBackarcLlsvpV2 {
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
                backarc_episode_my,
                backarc_max_stretch_factor,
                kimberlite_plume_radius_m,
            } => MagmatismLawV1::SlabArcBackarcLlsvpV3 {
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
                backarc_episode_my,
                backarc_max_stretch_factor,
                kimberlite_plume_radius_m,
            },
            other => other,
        }
    }

    /// SEED-WORLD-COMPOSE (lane CP): the v2 law from a v1 law's constants
    /// plus the three v2 constants; any other law unchanged.
    pub fn into_v2(self, backarc_episode_my: ConstantV1, backarc_max_stretch_factor: ConstantV1, kimberlite_plume_radius_m: ConstantV1) -> Self {
        match self {
            MagmatismLawV1::SlabArcBackarcLlsvpV1 {
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
            } => MagmatismLawV1::SlabArcBackarcLlsvpV2 {
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
                backarc_episode_my,
                backarc_max_stretch_factor,
                kimberlite_plume_radius_m,
            },
            other => other,
        }
    }
    /// SEED-WORLD-COMPOSE (lane CP): v2 viewed as v1 for the shared checks
    /// and the shared constants; any other law unchanged.
    pub fn shared_view(&self) -> Self {
        match self.clone() {
            MagmatismLawV1::SlabArcBackarcLlsvpV2 {
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
                ..
            }
            | MagmatismLawV1::SlabArcBackarcLlsvpV3 {
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
                ..
            } => MagmatismLawV1::SlabArcBackarcLlsvpV1 {
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
            },
            other => other,
        }
    }
    pub fn is_legacy(&self) -> bool {
        matches!(self, MagmatismLawV1::NoneV0)
    }
    fn legacy() -> Self {
        MagmatismLawV1::NoneV0
    }
}

/// Initial craton geometry and thickness (SEED-IRREGULAR-SEEDS, lane IS,
/// 2026-10-13; implementation `pf_geodynamics::cratons::seed_cratons`).
/// Refines `crust_seeding = geodesic-disc-bodies-v0`, whose continental
/// area fraction and rifted margin width it keeps. Serialized by ABSENCE
/// when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum CratonSeedingLawV1 {
    /// LEGACY: the crust seeding law's own shapes (discs or whole plates).
    NoneV0,
    /// `irregular-cratons-v1`: `tectonics.continental_seeds` assemblies of
    /// rank-size areas, each of 1 + floor(N_c share) overlapping elliptic
    /// cratons; coasts from a domain-warped score plus multi-scale edge
    /// noise, ranked to the declared area; continental cells ride their
    /// assembly's plate; interior thickness = base + basement roughness +
    /// suture roots on the seams; two-segment rifted margin (necking wedge,
    /// then a proximal shelf stretched by `proximal_margin_stretch` at the
    /// hinge). Noise wavelengths and shape numbers are fixed by the version.
    IrregularCratonsV1 {
        craton_interior_thickness_m: ConstantV1,
        craton_rank_size_exponent: ConstantV1,
        cratons_per_unit_share: ConstantV1,
        basement_roughness_std_m: ConstantV1,
        suture_crest_thickening_m: ConstantV1,
        suture_half_width_m: ConstantV1,
        proximal_margin_stretch: ConstantV1,
    },
}

impl CratonSeedingLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, CratonSeedingLawV1::NoneV0)
    }
    fn legacy() -> Self {
        CratonSeedingLawV1::NoneV0
    }
}

/// Large igneous provinces and plume swells (SEED-FLOODBASALT, lane FB,
/// 2026-10-13; implementation `pf_geodynamics::flood_basalt`). Serialized
/// by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum FloodBasaltLawV1 {
    /// LEGACY: a plume head is the `fixed-plume-chains-v1` head cone only;
    /// no swell.
    NoneV0,
    /// `plume-head-lip-swell-v1`: (1) a plume whose onset is inside the
    /// run (onset >= 0) erupts V = V_extrusive (1 + r_intrusive) of crust
    /// over [onset, onset + duration]; each step's share is placed where
    /// the plate now under the plume has carried the material erupted at
    /// the overlap midpoint, over a flat-topped disc of radius R (weight 1
    /// inside (1 - taper) R, linear to 0 at R, any plate, volume exact).
    /// Real crust: carried, eroded and compensated like any crust; booked
    /// with the plume crust. It supersedes the v1 head cone (plume_pass
    /// runs with head volume 0). (2) Every active plume lifts the columns
    /// around its FIXED mantle position by the water-loaded column term
    /// A_tail exp(-(r/R_tail)^2) + A_head exp(-(t - onset)/tau)
    /// exp(-(r/R_head)^2) (cut at 3 R; elevation only). Requires
    /// `fixed-plume-chains-v1` and the two-column isostasy law.
    PlumeHeadLipSwellV1 {
        flood_basalt_extrusive_volume_m3: ConstantV1,
        flood_basalt_intrusive_ratio: ConstantV1,
        flood_basalt_radius_m: ConstantV1,
        flood_basalt_taper_fraction: ConstantV1,
        flood_basalt_eruption_duration_my: ConstantV1,
        swell_tail_amplitude_m: ConstantV1,
        swell_tail_radius_m: ConstantV1,
        swell_head_amplitude_m: ConstantV1,
        swell_head_radius_m: ConstantV1,
        swell_head_decay_my: ConstantV1,
    },
}

impl FloodBasaltLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, FloodBasaltLawV1::NoneV0)
    }
    fn legacy() -> Self {
        FloodBasaltLawV1::NoneV0
    }
}

/// Post-rift state of the SEEDED rifted margins (SEED-TRENCH part 2, lane
/// TN, 2026-10-13; implementation `pf_lab::preview` init + the rift-flank
/// law). Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum PassiveMarginLawV1 {
    /// LEGACY: the seeded rifted-margin taper is an initial thickness only
    /// (no thermal state, no post-rift subsidence).
    NoneV0,
    /// `seeded-margin-thermal-subsidence-v1`: at t = 0 every seeded
    /// continental column thinner than the interior thickness h_cc carries
    /// the rift-flank stretching markers beta = h_cc / h and thermal start
    /// time = -initial oceanic age (breakup when the adjacent seeded ocean
    /// was born), so `shoulder-uplift-thermal-v1` gives the seeded passive
    /// margins their post-rift thermal subsidence. No constant. Requires
    /// `shoulder-uplift-thermal-v1` and `geodesic-disc-bodies-v0`.
    SeededMarginThermalSubsidenceV1,
}

impl PassiveMarginLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, PassiveMarginLawV1::NoneV0)
    }
    fn legacy() -> Self {
        PassiveMarginLawV1::NoneV0
    }
}

/// Sediment deposition on land (SEED-EXTRAS D7, lane SD, 2026-10-12;
/// implementation `pf_surface::erosion::erode_deposit_step`). Serialized by
/// ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum LandDepositionLawV1 {
    /// LEGACY: detachment-limited stream power; eroded material reaches
    /// standing water without depositing on land.
    DetachmentLimitedV0,
    /// `erosion-deposition-v1`: dh/dt = -K A^m S + G Qs/(p A) (Yuan et al.
    /// 2019; Davy and Lague 2009; recalled, not checked), solved implicitly
    /// on the drainage tree; p = 1. Composes with
    /// `stream-power-n1-linear-diffusion-v0` only (refused with the subgrid
    /// law until a combined implementation exists).
    ErosionDepositionV1 { deposition_coefficient_g: ConstantV1 },
}

impl LandDepositionLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, LandDepositionLawV1::DetachmentLimitedV0)
    }
    fn legacy() -> Self {
        LandDepositionLawV1::DetachmentLimitedV0
    }
}

/// Sediment placement in standing water (lakes and the ocean) (lane SD).
/// Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum SedimentPlacementLawV1 {
    /// LEGACY `local-fill-global-layer-v0`: uncompensated fill to the water
    /// level, 512-cell bounded spread, remainder as one uniform ocean layer
    /// (A-SEDIMENT-SPREAD-BOUND, A-MARINE-GLOBAL-LAYER).
    LocalFillGlobalLayerV0,
    /// `shelf-wave-base-v1`: isostatically compensated accommodation (the
    /// isostasy law's own densities), shelf top at sea level minus the wave
    /// base, contiguous progradation from each river mouth, declared and
    /// measured global remainder only beyond the bounded search.
    ShelfWaveBaseV1 { shelf_wave_base_depth_m: ConstantV1 },
}

impl SedimentPlacementLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, SedimentPlacementLawV1::LocalFillGlobalLayerV0)
    }
    fn legacy() -> Self {
        SedimentPlacementLawV1::LocalFillGlobalLayerV0
    }
}

/// Continental sediment retention and the base-level erosion taper
/// (SEED-SEDIMENT-RETAIN, lane SR, 2026-10-13; implementation
/// `pf_surface::erosion::retain`). Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum SedimentRetentionLawV1 {
    /// LEGACY: the erosion step of the selected land laws, unchanged.
    NoneV0,
    /// `continental-retention-base-level-v1`: (1) a land cell's stream-power
    /// lowering is multiplied by phi = min(1, (r / H)^p), r its local relief
    /// (max - min level over the cell and its 8 neighbours, water at its
    /// level; no change at r >= H; rule 2, lane SR PREREG.txt); (2)
    /// detached material is stored first on continental land below the
    /// graded level (receiver + S_a L), then in lakes, then on flooded
    /// continental cells up to wave base (all with compensated
    /// thickness), and only the excess goes to the legacy placement
    /// (oceanic cells). Requires legacy land deposition and sediment
    /// placement (it replaces both).
    ContinentalRetentionBaseLevelV1 {
        base_level_relief_scale_m: ConstantV1,
        base_level_taper_exponent: ConstantV1,
        alluvial_grade_slope: ConstantV1,
        shelf_wave_base_depth_m: ConstantV1,
    },
}

impl SedimentRetentionLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, SedimentRetentionLawV1::NoneV0)
    }
    fn legacy() -> Self {
        SedimentRetentionLawV1::NoneV0
    }
}

/// Stratified rock columns on the continental parcels (SEED-STRATIGRAPHY
/// stage A, lane ST, 2026-10-13; implementation
/// `pf_geodynamics::stratigraphy`; pre-registered in
/// `unpublished-design-record`). Requires
/// `continental-parcel-carry-v1`. A passive record: it reads the step's
/// recorded changes and writes nothing the simulation reads. Serialized by
/// ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum StratigraphyLawV1 {
    /// LEGACY: no rock column (rock_family v0 only).
    NoneV0,
    /// `parcel-stacks-v1`: every continental parcel carries a basement
    /// (the residual of its column) and up to max_layers - 1 cover layers;
    /// the step's per-cell changes are handed to the parcels as tagged
    /// deltas (erosion, deposition by environment, arc, plume crust,
    /// thinning, back-arc moves, shortening, basement residual) and applied
    /// per parcel in a declared order (removals, then additions bottom-up);
    /// cover deeper than the metamorphic depth becomes metamorphic; merges
    /// keep the stack within its capacity.
    ParcelStacksV1 {
        strata_max_layers: ConstantV1,
        strata_merge_gap_my: ConstantV1,
        strata_unconformity_gap_my: ConstantV1,
        strata_min_layer_m: ConstantV1,
        strata_trace_m: ConstantV1,
        strata_arc_extrusive_fraction: ConstantV1,
        strata_plume_extrusive_fraction: ConstantV1,
        strata_shortening_pluton_fraction: ConstantV1,
        strata_pluton_roof_m: ConstantV1,
        strata_metamorphic_depth_m: ConstantV1,
        strata_high_grade_depth_m: ConstantV1,
        strata_shelf_depth_m: ConstantV1,
        strata_carbonate_max_paleolat_deg: ConstantV1,
        strata_carbonate_max_supply_m_per_my: ConstantV1,
        strata_arid_lat_min_deg: ConstantV1,
        strata_arid_lat_max_deg: ConstantV1,
        strata_humid_lat_max_deg: ConstantV1,
        strata_humid_midlat_min_deg: ConstantV1,
        strata_humid_midlat_max_deg: ConstantV1,
        strata_coal_max_elevation_m: ConstantV1,
        strata_conglomerate_relief_m: ConstantV1,
        strata_river_mouth_area_m2: ConstantV1,
        strata_publish_depth_m: ConstantV1,
        strata_basement_min_m: ConstantV1,
    },
    /// `parcel-stacks-v2` (SEED-STRATIGRAPHY stage A.1, lane SA, 2026-10-13;
    /// `unpublished-design-record`): v1 plus plutons and
    /// pure-shear thickening only from the belt's placed collision supply
    /// (plateau, channel-flow and delamination moves go to the basement), a
    /// declared sub-cell shelf-slope ramp choosing the shelf facies, a t = 0
    /// platform cover by elevation, a capacity merge that keeps
    /// unconformities and merges old, similar, thin layers first, and a
    /// per-parcel cumulative rock uplift published with its gradient. Still a
    /// passive record.
    ParcelStacksV2 {
        strata_max_layers: ConstantV1,
        strata_merge_gap_my: ConstantV1,
        strata_unconformity_gap_my: ConstantV1,
        strata_min_layer_m: ConstantV1,
        strata_trace_m: ConstantV1,
        strata_arc_extrusive_fraction: ConstantV1,
        strata_plume_extrusive_fraction: ConstantV1,
        strata_shortening_pluton_fraction: ConstantV1,
        strata_pluton_roof_m: ConstantV1,
        strata_metamorphic_depth_m: ConstantV1,
        strata_high_grade_depth_m: ConstantV1,
        strata_shelf_depth_m: ConstantV1,
        strata_carbonate_max_paleolat_deg: ConstantV1,
        strata_carbonate_max_supply_m_per_my: ConstantV1,
        strata_arid_lat_min_deg: ConstantV1,
        strata_arid_lat_max_deg: ConstantV1,
        strata_humid_lat_max_deg: ConstantV1,
        strata_humid_midlat_min_deg: ConstantV1,
        strata_humid_midlat_max_deg: ConstantV1,
        strata_coal_max_elevation_m: ConstantV1,
        strata_conglomerate_relief_m: ConstantV1,
        strata_river_mouth_area_m2: ConstantV1,
        strata_publish_depth_m: ConstantV1,
        strata_basement_min_m: ConstantV1,
        strata_shelf_width_m: ConstantV1,
        strata_slope_foot_depth_m: ConstantV1,
        strata_shelf_majority_fraction: ConstantV1,
        strata_platform_cover_max_m: ConstantV1,
        strata_platform_full_elev_m: ConstantV1,
        strata_platform_zero_elev_m: ConstantV1,
        strata_platform_sandstone_fraction: ConstantV1,
        strata_platform_carbonate_fraction: ConstantV1,
        strata_platform_age_base_my: ConstantV1,
        strata_platform_age_top_my: ConstantV1,
    },
}

impl StratigraphyLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, StratigraphyLawV1::NoneV0)
    }
    fn legacy() -> Self {
        StratigraphyLawV1::NoneV0
    }
    /// The (id, constant) pairs of `parcel-stacks-v1` in declaration order.
    pub fn constants(&self) -> Vec<(&'static str, &ConstantV1)> {
        match self {
            StratigraphyLawV1::NoneV0 => Vec::new(),
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
            } => vec![
                ("strata_max_layers", strata_max_layers),
                ("strata_merge_gap_my", strata_merge_gap_my),
                ("strata_unconformity_gap_my", strata_unconformity_gap_my),
                ("strata_min_layer_m", strata_min_layer_m),
                ("strata_trace_m", strata_trace_m),
                ("strata_arc_extrusive_fraction", strata_arc_extrusive_fraction),
                ("strata_plume_extrusive_fraction", strata_plume_extrusive_fraction),
                ("strata_shortening_pluton_fraction", strata_shortening_pluton_fraction),
                ("strata_pluton_roof_m", strata_pluton_roof_m),
                ("strata_metamorphic_depth_m", strata_metamorphic_depth_m),
                ("strata_high_grade_depth_m", strata_high_grade_depth_m),
                ("strata_shelf_depth_m", strata_shelf_depth_m),
                ("strata_carbonate_max_paleolat_deg", strata_carbonate_max_paleolat_deg),
                ("strata_carbonate_max_supply_m_per_my", strata_carbonate_max_supply_m_per_my),
                ("strata_arid_lat_min_deg", strata_arid_lat_min_deg),
                ("strata_arid_lat_max_deg", strata_arid_lat_max_deg),
                ("strata_humid_lat_max_deg", strata_humid_lat_max_deg),
                ("strata_humid_midlat_min_deg", strata_humid_midlat_min_deg),
                ("strata_humid_midlat_max_deg", strata_humid_midlat_max_deg),
                ("strata_coal_max_elevation_m", strata_coal_max_elevation_m),
                ("strata_conglomerate_relief_m", strata_conglomerate_relief_m),
                ("strata_river_mouth_area_m2", strata_river_mouth_area_m2),
                ("strata_publish_depth_m", strata_publish_depth_m),
                ("strata_basement_min_m", strata_basement_min_m),
            ],
            StratigraphyLawV1::ParcelStacksV2 {
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
                strata_shelf_width_m,
                strata_slope_foot_depth_m,
                strata_shelf_majority_fraction,
                strata_platform_cover_max_m,
                strata_platform_full_elev_m,
                strata_platform_zero_elev_m,
                strata_platform_sandstone_fraction,
                strata_platform_carbonate_fraction,
                strata_platform_age_base_my,
                strata_platform_age_top_my,
            } => vec![
                ("strata_max_layers", strata_max_layers),
                ("strata_merge_gap_my", strata_merge_gap_my),
                ("strata_unconformity_gap_my", strata_unconformity_gap_my),
                ("strata_min_layer_m", strata_min_layer_m),
                ("strata_trace_m", strata_trace_m),
                ("strata_arc_extrusive_fraction", strata_arc_extrusive_fraction),
                ("strata_plume_extrusive_fraction", strata_plume_extrusive_fraction),
                ("strata_shortening_pluton_fraction", strata_shortening_pluton_fraction),
                ("strata_pluton_roof_m", strata_pluton_roof_m),
                ("strata_metamorphic_depth_m", strata_metamorphic_depth_m),
                ("strata_high_grade_depth_m", strata_high_grade_depth_m),
                ("strata_shelf_depth_m", strata_shelf_depth_m),
                ("strata_carbonate_max_paleolat_deg", strata_carbonate_max_paleolat_deg),
                ("strata_carbonate_max_supply_m_per_my", strata_carbonate_max_supply_m_per_my),
                ("strata_arid_lat_min_deg", strata_arid_lat_min_deg),
                ("strata_arid_lat_max_deg", strata_arid_lat_max_deg),
                ("strata_humid_lat_max_deg", strata_humid_lat_max_deg),
                ("strata_humid_midlat_min_deg", strata_humid_midlat_min_deg),
                ("strata_humid_midlat_max_deg", strata_humid_midlat_max_deg),
                ("strata_coal_max_elevation_m", strata_coal_max_elevation_m),
                ("strata_conglomerate_relief_m", strata_conglomerate_relief_m),
                ("strata_river_mouth_area_m2", strata_river_mouth_area_m2),
                ("strata_publish_depth_m", strata_publish_depth_m),
                ("strata_basement_min_m", strata_basement_min_m),
                ("strata_shelf_width_m", strata_shelf_width_m),
                ("strata_slope_foot_depth_m", strata_slope_foot_depth_m),
                ("strata_shelf_majority_fraction", strata_shelf_majority_fraction),
                ("strata_platform_cover_max_m", strata_platform_cover_max_m),
                ("strata_platform_full_elev_m", strata_platform_full_elev_m),
                ("strata_platform_zero_elev_m", strata_platform_zero_elev_m),
                ("strata_platform_sandstone_fraction", strata_platform_sandstone_fraction),
                ("strata_platform_carbonate_fraction", strata_platform_carbonate_fraction),
                ("strata_platform_age_base_my", strata_platform_age_base_my),
                ("strata_platform_age_top_my", strata_platform_age_top_my),
            ],
        }
    }
}

/// Plate-boundary classification and the carry's boundary rules
/// (SEED-TRANSFORMS, lane TF, 2026-10-13; implementation
/// `pf_geodynamics::kinematic::{classify_cells_v2,
/// advance_with_parcels_transform}`; pre-registered in
/// `unpublished-design-record`). Serialized by ABSENCE when
/// legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum PlateBoundaryLawV1 {
    /// LEGACY: each boundary cell is classified by the normal rate along the
    /// direction to the neighbouring plate's cell (the grid staircase),
    /// transform only below 1e-4 rad/My; every gap is a ridge birth and
    /// every oceanic overlap subducts.
    NeighbourDirectionV0,
    /// `smoothed-normal-transform-v1`: the boundary normal is the first
    /// moment of the two plates' indicators over a Gaussian window (sigma
    /// `boundary_normal_smoothing_m`, cut at 3 sigma, area weighted; the
    /// staircase direction where degenerate); the relative velocity (w_a -
    /// w_b) x r is a TRANSFORM when its obliquity to the trace is below
    /// `transform_obliquity_deg`, else convergent / divergent by the sign
    /// of its normal part. In the carry, gaps on a transform or convergent
    /// step-start class take the same side's carried material (2 layers;
    /// no ridge birth), and oceanic overlaps on a transform or divergent
    /// class retire the loser without subduction (no arc sample); both are
    /// named remap terms. Requires `continental-parcel-carry-v1`.
    SmoothedNormalTransformV1 {
        transform_obliquity_deg: ConstantV1,
        boundary_normal_smoothing_m: ConstantV1,
    },
}

impl PlateBoundaryLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, PlateBoundaryLawV1::NeighbourDirectionV0)
    }
    fn legacy() -> Self {
        PlateBoundaryLawV1::NeighbourDirectionV0
    }
}

/// PLATE-POLES (lane PLATE-POLES, 2026-10-21; implementation
/// `pf_geodynamics::boundaries::LAW_RIGID_GROUP_CONTACTS_V1`,
/// `kinematic::classify_cells_v2`, `reorg::reorg_pass_release`): what a
/// suture group is as a rigid body. Serialized by ABSENCE when legacy, so
/// every existing law set keeps its bytes; when selected it enters the
/// universe digest and `same_physics`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum RigidGroupLawV1 {
    /// LEGACY: plates of one suture group share one pole wherever they are;
    /// their mutual contacts (zero relative velocity) are classified as
    /// transforms; a member detached from the group keeps moving with it.
    SharedPoleV0,
    /// `rigid-group-interior-v1`: a contact between plates with
    /// bit-identical angular velocity is plate interior (no boundary class,
    /// no boundary record; the carry's remap class keeps the same-body
    /// remap rule there, so the motion and material are unchanged). No new
    /// constant. Requires `smoothed-normal-transform-v1`.
    RigidGroupInteriorV1,
    /// `rigid-group-contacts-v1`: `rigid-group-interior-v1` plus release: at
    /// each scheduled reorganization event, after the fragment pass and
    /// before the re-solve, group members that touch no other member of the
    /// group's main body (largest-area connected set) leave the group, so
    /// the re-solve gives them their own pole (a weld needs contact). No
    /// new constant. The release part acts only under a plate-reorganization
    /// law.
    RigidGroupContactsV1,
}

impl RigidGroupLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, RigidGroupLawV1::SharedPoleV0)
    }
    fn legacy() -> Self {
        RigidGroupLawV1::SharedPoleV0
    }
}
/// CLIMATE-GENESIS (lane CG, 2026-10-20): how genesis runoff is carried
/// between two climate solve epochs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ClimateBetweenEpochsV1 {
    /// `frozen-column-local-rainout-v1`: keep the epoch's column water,
    /// ocean-influence tracer and zonal rows; every step recompute
    /// precipitation, evapotranspiration and runoff pointwise from the
    /// current surface (orographic rain follows the ranges as plates move;
    /// the rain shadow's upwind depletion lags until the next solve).
    FrozenColumnLocalRainoutV1,
    /// `hold-runoff-v1`: keep the epoch's runoff field in the geographic
    /// frame until the next solve.
    HoldRunoffV1,
}

/// CLIMATE-GENESIS: what the fluvial kernels receive as runoff.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ClimateErosivityModeV1 {
    /// `climate-runoff-pattern-v1`: the modelled runoff (the law).
    ClimateRunoffPatternV1,
    /// DIAGNOSTIC: runoff = the reference runoff everywhere, so the
    /// runoff-weighted area equals the geometric area bit for bit (proves
    /// the path is exact; legacy elevation bytes).
    UniformReferenceDiag,
    /// DIAGNOSTIC: every step's land runoff replaced by one uniform value
    /// with the same area-weighted land mean of R^m (same mean erosivity,
    /// no spatial pattern).
    LandMeanEquivalentDiag,
    /// DIAGNOSTIC (lane CLIMATE-GENESIS-2, 2026-10-21): every step's land
    /// runoff replaced by its area-weighted land MEAN R (no pattern). The
    /// mean-matched control when the fluvial terms are linear in R: the
    /// n = 2 subgrid law has m = n/2 = 1 on the trunk A_eff and (with
    /// cell-runoff-subgrid-erosivity-v1) on the cell mean, so
    /// `land-mean-equivalent-diag` (matched in R^0.5) is not mean-matched
    /// there.
    LandMeanRunoffDiag,
}

/// I11 (lane CLIMATE-GENESIS-2, 2026-10-20): whether climate runoff also
/// reaches the subgrid cell-mean (hillslope-channel) relaxation of the
/// subgrid erosion laws. Serialized by ABSENCE at the v1 default, so every
/// existing climate law keeps its bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ClimateSubgridCellMeanV1 {
    /// v1 default: the cell-mean relaxation runs at the reference runoff
    /// (climate reaches the subgrid law only through the trunk's A_eff).
    #[default]
    ReferenceRunoffV1,
    /// `cell-runoff-subgrid-erosivity-v1`: the unresolved network inside a
    /// cell drains that cell's own area, so its discharge is a(x) R_i and
    /// the closure becomes dz/dt = -K_n (R_i / r_ref)^m ((z - c)/Lambda)^n
    /// (m = n/2, the law's own area exponent), R_i the cell's runoff fed to
    /// erosion. The critical-slope floor and hillslope diffusion stay
    /// climate-independent.
    CellRunoffSubgridErosivityV1,
}

/// HYPSO-LEVERS (2026-10-21): how the land runoff R of
/// `epoch-climate-runoff-v1` becomes the ERODING runoff R_e that the
/// fluvial kernels read (A_eff = sum a R_e / r_ref and, with
/// `cell-runoff-subgrid-erosivity-v1`, the cell-mean ratio R_e / r_ref).
/// Absent = the v1 law (R_e = R, linear). The lake water balance and every
/// published climate field keep the modelled runoff R.
///
/// Why (lane HYPSO-CENSUS, seed 42): runoff on land above 2 km averaged
/// 3.1-4.1 m/yr against 0.24 m/yr on lowland, and with the n = 2 subgrid
/// law (m = 1) it entered erosivity linearly, so the climate recipe eroded
/// crust thicker than 36 km at 2.7 %/My against 1.4 %/My at uniform 1 m/yr.
/// Long-term bedrock erosion is expected to respond sublinearly to mean
/// runoff where an incision threshold and discharge variability matter
/// (stochastic-threshold models, Tucker 2004, Lague et al. 2005, and
/// DiBiase and Whipple 2011; and cosmogenic compilations correlate erosion
/// more with relief than with precipitation, Portenga and Bierman 2011):
/// all RECALLED, not checked; project assumption awaiting review.
///
/// Measured (lane HYPSO-LEVERS, climate-v1 + lake keys, one run per seed,
/// controlled windows before plate ids diverge): p = 0.3 raised land > 1 km
/// by only +0.13 / +0.09 / +0.09 points of land (seeds 42 / 240907 / 7),
/// erosion -1 to -3 %, drylands unchanged, lakes -1 to -6 %; p = 0.5 and a
/// 1.5 m/yr cap did the same or less on seed 42. Runoff reaches only the
/// trunk A_eff (the cell-mean closure runs at r_ref unless
/// `cell-runoff-subgrid-erosivity-v1`), and the uniform-reference diagnostic
/// had LESS high land (6.8 vs 7.9 %): the climate recipe's high-land
/// shortfall against village-c7-v1 comes from its K 2e-6 (against 1e-6),
/// not from orographic runoff. Not a default.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ErosiveRunoffLawV1 {
    /// `erosive-runoff-sublinear-above-reference-v1`: R_e = R for R <=
    /// r_ref, R_e = r_ref (R / r_ref)^p above (continuous at r_ref; dry and
    /// lowland runoff below r_ref unchanged).
    SublinearAboveReferenceV1 { erosive_runoff_exponent: ConstantV1 },
    /// `erosive-runoff-cap-v1`: R_e = min(R, cap).
    CapV1 { erosive_runoff_cap_m_per_yr: ConstantV1 },
}

impl ErosiveRunoffLawV1 {
    pub fn law_id(&self) -> &'static str {
        match self {
            ErosiveRunoffLawV1::SublinearAboveReferenceV1 { .. } => "erosive-runoff-sublinear-above-reference-v1",
            ErosiveRunoffLawV1::CapV1 { .. } => "erosive-runoff-cap-v1",
        }
    }
    pub fn validate(&self) -> Result<(), UniverseError> {
        match self {
            ErosiveRunoffLawV1::SublinearAboveReferenceV1 { erosive_runoff_exponent } => check("erosive_runoff_exponent", erosive_runoff_exponent),
            ErosiveRunoffLawV1::CapV1 { erosive_runoff_cap_m_per_yr } => check("erosive_runoff_cap_m_per_yr", erosive_runoff_cap_m_per_yr),
        }
    }
}

impl ClimateSubgridCellMeanV1 {
    pub fn is_default(&self) -> bool {
        *self == ClimateSubgridCellMeanV1::ReferenceRunoffV1
    }
}

/// CLIMATE-GENESIS (lane CG, 2026-10-20; design
/// `unpublished-design-record`): climate forcing of the
/// genesis fluvial kernels. Serialized by ABSENCE when legacy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ClimateForcingLawV1 {
    /// LEGACY: stream power on geometric contributing area (uniform runoff
    /// enters discharge only).
    NoneV0,
    /// `epoch-climate-runoff-v1`: a `seed-climate-v1`/`-v2` solve every
    /// `climate_epoch_my` (on a grid coarsened by
    /// `climate_solve_grid_divisor`, warm-started, to
    /// `climate_sweep_rel_tol`); land runoff R = max(P - E, f P) (f =
    /// `climate_runoff_floor_fraction`), lakes max(P - PET, 0); every
    /// fluvial kernel reads A_eff = sum over the catchment of a R / r_ref
    /// (approximation `runoff-weighted-area-v1`, r_ref =
    /// `surface.uniform_runoff_m_per_yr`), so E = K A_eff^m S^n with K quoted
    /// at r_ref. The climate recipe and the digests of the tier's extra
    /// constants are part of the law (ADR 044: a changed climate constant is
    /// a different world).
    EpochClimateRunoffV1 {
        climate: ClimateRecipeV1,
        seed_constants_digest: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        rotation_constants_digest: Option<String>,
        climate_epoch_my: ConstantV1,
        climate_solve_grid_divisor: ConstantV1,
        climate_sweep_rel_tol: ConstantV1,
        climate_runoff_floor_fraction: ConstantV1,
        between_epochs: ClimateBetweenEpochsV1,
        erosivity: ClimateErosivityModeV1,
        /// I11 (lane CLIMATE-GENESIS-2): skipped on write at the v1 default.
        #[serde(default, skip_serializing_if = "ClimateSubgridCellMeanV1::is_default")]
        subgrid_cell_mean: ClimateSubgridCellMeanV1,
        /// HYPSO-LEVERS: erosive-runoff response; absent = linear (v1),
        /// skipped on write so every existing climate law keeps its bytes.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        erosive_runoff: Option<ErosiveRunoffLawV1>,
    },
}

impl ClimateForcingLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, ClimateForcingLawV1::NoneV0)
    }
    fn legacy() -> Self {
        ClimateForcingLawV1::NoneV0
    }
}

/// SHELF law-set entry (lane CLIMATE-GENESIS-2, 2026-10-20; physics from
/// lane SHELF, implementation `pf_geodynamics::shelf_grade`): marine
/// reworking of shallow ocean floor to wave base. Serialized by ABSENCE
/// when legacy, so every existing law set keeps its bytes; when selected it
/// enters the universe digest and `same_physics` (a pinned universe with
/// the law off refuses a profile with it on, and the reverse).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ShelfGradingLawV1 {
    /// LEGACY: no grading (sediment on shallow ocean floor stays where the
    /// retention / placement laws put it).
    NoneV0,
    /// `wave-base-shelf-grading-v1`: after each erosion step, connected
    /// open ocean shallower than the retention law's
    /// `shelf_wave_base_depth_m` loses sediment (never bedrock) down to wave
    /// base; the volume goes to the nearest deeper ocean cells, each filled
    /// at most to wave base (the shelf edge progrades); the surface follows
    /// through the retention law's Airy sediment column. Mass conserved (an
    /// internal ledger transfer). No new constant: the wave base, densities
    /// and search bound are the retention law's. Grading completing within
    /// one step is a project assumption awaiting review (lane SHELF).
    /// Requires `continental-retention-base-level-v1`.
    WaveBaseShelfGradingV1,
}

impl ShelfGradingLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, ShelfGradingLawV1::NoneV0)
    }
    fn legacy() -> Self {
        ShelfGradingLawV1::NoneV0
    }
}

/// INLAND-BASINS (lane INLAND-BASINS, 2026-10-21): closed-basin hydrology
/// and lake infill (pf-surface closed_basin.rs). Serialized by ABSENCE when
/// legacy, so every existing law set keeps its bytes; when selected it
/// enters the universe digest and `same_physics`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ClosedBasinLawV1 {
    /// LEGACY: every closed depression is a lake full to its spill level
    /// (routing); lake sediment builds inlet deltas only.
    NoneV0,
    /// `closed-basin-water-balance-v1`: each routing lake's water level from
    /// its balance (runoff inflow + upstream overflow against net open-water
    /// evaporation; full at the spill level or terminal below it), and the
    /// suspended share of lake sediment focused into the deepest water as a
    /// level floor fill (the bedload share keeps the retention law's inlet
    /// delta). Requires `continental-retention-base-level-v1`.
    WaterBalanceFocusedFillV1 {
        suspended_load_fraction: ConstantV1,
        uniform_evaporation_to_runoff_ratio: ConstantV1,
    },
    /// LAKE-EXCESS (lane LAKE-EXCESS, 2026-10-21):
    /// `closed-basin-water-balance-v1` plus `filled-basin-alluvial-grade-v1`:
    /// a lake whose floor fill reaches within `filled_basin_grading_depth_m`
    /// of its spill level is a filled basin, and its sediment builds an
    /// alluvial plain graded to the outlet (spill + S_a x flow distance,
    /// S_a the retention law's graded alluvial slope) instead of a level
    /// floor. A separate variant, so the v1 variant keeps its bytes.
    WaterBalanceFocusedFillAlluvialGradeV1 {
        suspended_load_fraction: ConstantV1,
        uniform_evaporation_to_runoff_ratio: ConstantV1,
        filled_basin_grading_depth_m: ConstantV1,
    },
}

impl ClosedBasinLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, ClosedBasinLawV1::NoneV0)
    }
    fn legacy() -> Self {
        ClosedBasinLawV1::NoneV0
    }
}

/// LAKE-OUTLET (lane LAKE-OUTLET, 2026-10-21): outlet incision of full
/// lakes (pf-surface closed_basin.rs rule 5). Serialized by ABSENCE when
/// legacy, so every existing law set keeps its bytes. Requires a
/// `closed-basin-water-balance-v1` variant (the overflow comes from its
/// balance). No constants of its own: the channel law is the selected
/// relief law's trunk law (with the incision threshold when selected).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum LakeOutletLawV1 {
    /// LEGACY: a lake's pour cell is lowered by the step's ordinary
    /// (subgrid-closure, tapered) erosion only.
    NoneV0,
    /// `lake-outlet-sill-incision-v1`: each full lake's pour cell is
    /// lowered by channel stream power on the lake's overflow discharge
    /// (no taper, no subgrid closure; the larger of that and the step's own
    /// lowering), so the spill level and the lake level follow.
    SillIncisionV1,
    /// SILL-REBOUND (lane SILL-REBOUND, 2026-10-21):
    /// `lake-outlet-sill-incision-flexural-support-v1`: `SillIncisionV1`,
    /// and a sill's lowering gets no local isostatic rebound (an outlet
    /// channel is far narrower than the flexural wavelength, so the plate
    /// carries it): the material removed at a sill cell is sized so the
    /// post-isostasy surface drops by the lowering (the retention law's
    /// compensated column, as the critical-slope rule sizes its removal).
    /// No constants of its own.
    SillIncisionFlexuralSupportV1,
    /// LAKE-HOLD (lane LAKE-HOLD, 2026-10-21):
    /// `lake-outlet-path-incision-v1`: `SillIncisionV1`, and the outlet
    /// river's channel law acts on every dry cell of the overflow path from
    /// the sill to the next water body (ocean, another lake or a terminal
    /// floor), not on the sill alone: each such cell is lowered by channel
    /// stream power on the overflow discharge (no taper, no subgrid
    /// closure; the larger of that and the step's own lowering). No
    /// constants of its own.
    PathIncisionV1,
    /// LAKE-HOLD: `lake-outlet-path-incision-flexural-support-v1`:
    /// `PathIncisionV1` with `SillIncisionFlexuralSupportV1`'s
    /// compensation on every incised path cell. No constants of its own.
    PathIncisionFlexuralSupportV1,
    /// LAKE-HOLD-2 (lane LAKE-HOLD-2, 2026-10-21):
    /// `lake-outlet-path-through-lakes-v1`: `PathIncisionV1`, and a cell
    /// whose receiver is a full lake whose sill the same step cuts is graded
    /// to that lake's post-incision spill (the cut passes down a chain of
    /// lakes within one step). No constants of its own.
    PathThroughLakesV1,
    /// LAKE-HOLD-2: `lake-outlet-path-through-lakes-flexural-support-v1`:
    /// `PathThroughLakesV1` with `PathIncisionFlexuralSupportV1`'s
    /// compensation. No constants of its own.
    PathThroughLakesFlexuralSupportV1,
}

impl LakeOutletLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, LakeOutletLawV1::NoneV0)
    }
    fn legacy() -> Self {
        LakeOutletLawV1::NoneV0
    }
    /// LAKE-HOLD-2: the variants that incise the whole overflow path.
    pub fn incises_path(&self) -> bool {
        matches!(
            self,
            LakeOutletLawV1::PathIncisionV1
                | LakeOutletLawV1::PathIncisionFlexuralSupportV1
                | LakeOutletLawV1::PathThroughLakesV1
                | LakeOutletLawV1::PathThroughLakesFlexuralSupportV1
        )
    }
}

/// LAKE-HOLD-2 (lane LAKE-HOLD-2, 2026-10-21): the graded-plain slope of
/// rule 2a on a lake's outlet path (pf-surface closed_basin.rs rule 9).
/// Serialized by ABSENCE when legacy, so every existing law set keeps its
/// bytes. Requires a path-incising `lake_outlet` variant and the fluvial
/// incision threshold (its k_c). No constants of its own.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum LakeOutletGradeLawV1 {
    /// LEGACY: rule 2a grades every continental cell to receiver + S_a L.
    NoneV0,
    /// `lake-outlet-path-threshold-grade-v1`: on the sill and the overflow
    /// path of a full lake (cells with outlet discharge area A_q > 0) the
    /// grade is min(S_a, k_c / A_q^(m/n)), the outlet channel's own
    /// threshold slope: a plain built steeper than that is cut again by the
    /// same channel law, so storage stops at the slope the overflow river can
    /// hold. No constants of its own.
    PathThresholdGradeV1,
}

impl LakeOutletGradeLawV1 {
    pub fn is_legacy(&self) -> bool {
        matches!(self, LakeOutletGradeLawV1::NoneV0)
    }
    fn legacy() -> Self {
        LakeOutletGradeLawV1::NoneV0
    }
}

/// CLIMATE-GENESIS: the climate law quotes K at the reference runoff
/// r_ref; RELIEF-DECAY's uniform-runoff recipe K 1.0e-6 /yr (= 2e-6 x
/// (0.27 m/yr / 1 m/yr)^0.5, Earth's land runoff folded into K) already
/// carries a runoff factor. Under `epoch-climate-runoff-v1` a K below this
/// value is refused so the runoff factor is never applied twice
/// (integrator ruling 2026-10-20: "alternatives, never stacked").
/// Project assumption awaiting review: the boundary is the geometric mean
/// of the two recipes' K (1.0e-6 and 2e-6), rounded up.
pub const CLIMATE_FORCING_MIN_REFERENCE_K_PER_YR: f64 = 1.5e-6;

/// The law set a world runs under (v1: only the entries the current
/// generator path executes).
/// Q25/Q29: perfectly mixed terrigenous fines in each connected ocean.
/// Parameters are explicit dimensionless share and positive-down depth in metres.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum MarineFinesLawV1 {
    #[serde(rename = "connected-ocean-mixed-terrigenous-fines-v1")]
    ConnectedOceanMixedTerrigenousV1 { fraction: f64, min_depth_m: f64 },
}
impl MarineFinesLawV1 {
    pub const ID: &'static str = "connected-ocean-mixed-terrigenous-fines-v1";
    pub fn parameters(&self) -> (f64, f64) {
        match *self { Self::ConnectedOceanMixedTerrigenousV1 { fraction, min_depth_m } => (fraction, min_depth_m) }
    }
    pub fn validate(&self) -> Result<(), UniverseError> {
        let (fraction, depth) = self.parameters();
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction)
            || !depth.is_finite() || depth < 0.0
        {
            return Err(UniverseError::Unsupported {
                law: Self::ID,
                why: "requires finite fraction in [0,1] and finite min_depth_m >= 0",
            });
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LawSetV1 {
    pub schema: Cow<'static, str>,
    /// Q29: absent on all historical profiles, including pinned universes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub marine_fines: Option<MarineFinesLawV1>,
    pub isostasy: IsostasyLawV1,
    pub sea_level: SeaLevelLawV1,
    pub crust_seeding: CrustSeedingLawV1,
    pub oceanic_thermal: OceanicThermalLawV1,
    pub erosion: ErosionLawV1,
    /// SEED-OROGEN: skipped on write when legacy, defaulted on read, so
    /// every pre-existing law set round-trips to identical bytes.
    #[serde(
        default = "TectonicThickeningLawV1::legacy",
        skip_serializing_if = "TectonicThickeningLawV1::is_legacy"
    )]
    pub tectonic_thickening: TectonicThickeningLawV1,
    /// SEED-EXTRAS D7 (lane SD): skipped on write when legacy, defaulted
    /// on read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "LandDepositionLawV1::legacy",
        skip_serializing_if = "LandDepositionLawV1::is_legacy"
    )]
    pub land_deposition: LandDepositionLawV1,
    #[serde(
        default = "SedimentPlacementLawV1::legacy",
        skip_serializing_if = "SedimentPlacementLawV1::is_legacy"
    )]
    pub sediment_placement: SedimentPlacementLawV1,
    /// SEED-CARRY-CC (lane CC): skipped on write when legacy, defaulted on
    /// read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "PlateCarryLawV1::legacy",
        skip_serializing_if = "PlateCarryLawV1::is_legacy"
    )]
    pub plate_carry: PlateCarryLawV1,
    /// SEDIMENT-TOWERS (lane ST): skipped on write when legacy, defaulted
    /// on read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "CollisionSedimentLawV1::legacy",
        skip_serializing_if = "CollisionSedimentLawV1::is_legacy"
    )]
    pub collision_sediment: CollisionSedimentLawV1,
    /// SEED-EXTRAS (lane EX): skipped on write when legacy, defaulted on
    /// read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "HotspotVolcanismLawV1::legacy",
        skip_serializing_if = "HotspotVolcanismLawV1::is_legacy"
    )]
    pub hotspot_volcanism: HotspotVolcanismLawV1,
    #[serde(
        default = "ContinentalRiftLawV1::legacy",
        skip_serializing_if = "ContinentalRiftLawV1::is_legacy"
    )]
    pub continental_rift: ContinentalRiftLawV1,
    #[serde(
        default = "PlateReorganizationLawV1::legacy",
        skip_serializing_if = "PlateReorganizationLawV1::is_legacy"
    )]
    pub plate_reorganization: PlateReorganizationLawV1,
    /// SEED-TRANSFORMS (lane TF): skipped on write when legacy, defaulted on
    /// read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "PlateBoundaryLawV1::legacy",
        skip_serializing_if = "PlateBoundaryLawV1::is_legacy"
    )]
    pub plate_boundary: PlateBoundaryLawV1,
    /// SEED-RIFT-SHOULDERS (lane RS): skipped on write when legacy,
    /// defaulted on read (every pre-existing law set round-trips to
    /// identical bytes).
    #[serde(
        default = "RiftFlankLawV1::legacy",
        skip_serializing_if = "RiftFlankLawV1::is_legacy"
    )]
    pub rift_flank: RiftFlankLawV1,
    /// SEED-UPLANDS (lane UP): skipped on write when legacy, defaulted on
    /// read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "HillslopeThresholdLawV1::legacy",
        skip_serializing_if = "HillslopeThresholdLawV1::is_legacy"
    )]
    pub hillslope_threshold: HillslopeThresholdLawV1,
    /// INCISION-THRESHOLD (lane IT): skipped on write when legacy,
    /// defaulted on read (every pre-existing law set round-trips to
    /// identical bytes).
    #[serde(
        default = "FluvialIncisionThresholdLawV1::legacy",
        skip_serializing_if = "FluvialIncisionThresholdLawV1::is_legacy"
    )]
    pub fluvial_incision_threshold: FluvialIncisionThresholdLawV1,
    /// ALLUVIAL-TRUNK (lane ALLUVIAL-TRUNK): skipped on write when legacy,
    /// defaulted on read (every pre-existing law set keeps its bytes).
    #[serde(
        default = "FluvialAlluvialChannelLawV1::legacy",
        skip_serializing_if = "FluvialAlluvialChannelLawV1::is_legacy"
    )]
    pub fluvial_alluvial_channel: FluvialAlluvialChannelLawV1,
    /// GLACIAL-T0 (lane GLACIAL-T0, 2026-10-21): skipped on write when
    /// legacy, defaulted on read (every pre-existing law set round-trips to
    /// identical bytes).
    #[serde(
        default = "GlacialErosionLawV1::legacy",
        skip_serializing_if = "GlacialErosionLawV1::is_legacy"
    )]
    pub glacial_erosion: GlacialErosionLawV1,
    /// INCISION-THRESHOLD part 1 (lane IT): skipped on write when legacy,
    /// defaulted on read.
    #[serde(
        default = "MarineRemainderLayerLawV1::legacy",
        skip_serializing_if = "MarineRemainderLayerLawV1::is_legacy"
    )]
    pub marine_remainder_layer: MarineRemainderLayerLawV1,
    /// SEED-SEDIMENT-RETAIN (lane SR): skipped on write when legacy,
    /// defaulted on read (every pre-existing law set round-trips to
    /// identical bytes).
    #[serde(
        default = "SedimentRetentionLawV1::legacy",
        skip_serializing_if = "SedimentRetentionLawV1::is_legacy"
    )]
    pub sediment_retention: SedimentRetentionLawV1,
    /// SEED-TRENCH (lane TN): skipped on write when legacy, defaulted on
    /// read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "SubductionTrenchLawV1::legacy",
        skip_serializing_if = "SubductionTrenchLawV1::is_legacy"
    )]
    pub subduction_trench: SubductionTrenchLawV1,
    /// SEED-TRENCH part 2 (lane TN): skipped on write when legacy,
    /// defaulted on read (every pre-existing law set round-trips to
    /// identical bytes).
    #[serde(
        default = "PassiveMarginLawV1::legacy",
        skip_serializing_if = "PassiveMarginLawV1::is_legacy"
    )]
    pub passive_margin: PassiveMarginLawV1,
    /// SEED-FLOODBASALT (lane FB): skipped on write when legacy, defaulted
    /// on read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "FloodBasaltLawV1::legacy",
        skip_serializing_if = "FloodBasaltLawV1::is_legacy"
    )]
    pub flood_basalt: FloodBasaltLawV1,
    /// SEED-CARRY-CONSERVE (lane CV): skipped on write when legacy,
    /// defaulted on read (every pre-existing law set round-trips to
    /// identical bytes).
    #[serde(
        default = "PlateSutureLawV1::legacy",
        skip_serializing_if = "PlateSutureLawV1::is_legacy"
    )]
    pub plate_suture: PlateSutureLawV1,
    /// SEED-IRREGULAR-SEEDS (lane IS): skipped on write when legacy,
    /// defaulted on read (every pre-existing law set round-trips to
    /// identical bytes).
    #[serde(
        default = "CratonSeedingLawV1::legacy",
        skip_serializing_if = "CratonSeedingLawV1::is_legacy"
    )]
    pub craton_seeding: CratonSeedingLawV1,
    /// SEED-RANGES (lane RG): skipped on write when legacy, defaulted on
    /// read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "OrogenPlateauLawV1::legacy",
        skip_serializing_if = "OrogenPlateauLawV1::is_legacy"
    )]
    pub orogen_plateau: OrogenPlateauLawV1,
    /// SEED-MAGMATISM (lane MG): skipped on write when legacy, defaulted on
    /// read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "MagmatismLawV1::legacy",
        skip_serializing_if = "MagmatismLawV1::is_legacy"
    )]
    pub magmatism: MagmatismLawV1,
    /// SEED-STRATIGRAPHY (lane ST): skipped on write when legacy, defaulted
    /// on read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "StratigraphyLawV1::legacy",
        skip_serializing_if = "StratigraphyLawV1::is_legacy"
    )]
    pub stratigraphy: StratigraphyLawV1,
    /// CLIMATE-GENESIS (lane CG): skipped on write when legacy, defaulted
    /// on read (every pre-existing law set round-trips to identical bytes).
    #[serde(
        default = "ClimateForcingLawV1::legacy",
        skip_serializing_if = "ClimateForcingLawV1::is_legacy"
    )]
    pub climate_forcing: ClimateForcingLawV1,
    /// SHELF (lane SHELF; law-set entry lane CLIMATE-GENESIS-2): skipped on
    /// write when legacy, defaulted on read (every pre-existing law set
    /// round-trips to identical bytes).
    #[serde(
        default = "ShelfGradingLawV1::legacy",
        skip_serializing_if = "ShelfGradingLawV1::is_legacy"
    )]
    pub shelf_grading: ShelfGradingLawV1,
    /// INLAND-BASINS (lane INLAND-BASINS): skipped on write when legacy,
    /// defaulted on read.
    #[serde(
        default = "ClosedBasinLawV1::legacy",
        skip_serializing_if = "ClosedBasinLawV1::is_legacy"
    )]
    pub closed_basin: ClosedBasinLawV1,
    /// LAKE-OUTLET (lane LAKE-OUTLET): skipped on write when legacy,
    /// defaulted on read.
    #[serde(
        default = "LakeOutletLawV1::legacy",
        skip_serializing_if = "LakeOutletLawV1::is_legacy"
    )]
    pub lake_outlet: LakeOutletLawV1,
    /// LAKE-HOLD-2 (lane LAKE-HOLD-2): skipped on write when legacy,
    /// defaulted on read.
    #[serde(
        default = "LakeOutletGradeLawV1::legacy",
        skip_serializing_if = "LakeOutletGradeLawV1::is_legacy"
    )]
    pub lake_outlet_grade: LakeOutletGradeLawV1,
    /// PLATE-POLES (lane PLATE-POLES): skipped on write when legacy,
    /// defaulted on read.
    #[serde(
        default = "RigidGroupLawV1::legacy",
        skip_serializing_if = "RigidGroupLawV1::is_legacy"
    )]
    pub rigid_group: RigidGroupLawV1,
}

impl LawSetV1 {
    /// True for the legacy preview law profile (no physics flag).
    pub fn is_legacy(&self) -> bool {
        self.marine_fines.is_none()
            && matches!(self.isostasy, IsostasyLawV1::SingleColumnReferenceV0 { .. })
            && self.sea_level == SeaLevelLawV1::PrescribedDatumV0
            && self.crust_seeding == CrustSeedingLawV1::WholePlateContinentsV0
            && self.tectonic_thickening.is_legacy()
            && self.erosion.is_legacy()
            && self.land_deposition.is_legacy()
            && self.sediment_placement.is_legacy()
            && self.plate_carry.is_legacy()
            && self.collision_sediment.is_legacy()
            && self.hotspot_volcanism.is_legacy()
            && self.continental_rift.is_legacy()
            && self.plate_reorganization.is_legacy()
            && self.plate_boundary.is_legacy()
            && self.rift_flank.is_legacy()
            && self.hillslope_threshold.is_legacy()
            && self.fluvial_incision_threshold.is_legacy()
            && self.fluvial_alluvial_channel.is_legacy()
            && self.glacial_erosion.is_legacy()
            && self.marine_remainder_layer.is_legacy()
            && self.sediment_retention.is_legacy()
            && self.subduction_trench.is_legacy()
            && self.passive_margin.is_legacy()
            && self.flood_basalt.is_legacy()
            && self.plate_suture.is_legacy()
            && self.craton_seeding.is_legacy()
            && self.orogen_plateau.is_legacy()
            && self.magmatism.is_legacy()
            && self.stratigraphy.is_legacy()
            && self.climate_forcing.is_legacy()
            && self.shelf_grading.is_legacy()
            && self.closed_basin.is_legacy()
            && self.lake_outlet.is_legacy()
            && self.lake_outlet_grade.is_legacy()
            && self.rigid_group.is_legacy()
    }
}

/// Physical constants pinned per world (v1: the generator's current inputs).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UniverseConstantsV1 {
    pub schema: Cow<'static, str>,
    // Planet
    pub planet_radius_m: ConstantV1,
    /// STATED, not derived: planet mass is not a v1 constant. Read only by
    /// the flexure kernel, which the preview tier does not run.
    pub surface_gravity_m_s2: ConstantV1,
    // Materials
    pub mantle_density_kg_m3: ConstantV1,
    pub continental_crust_density_kg_m3: ConstantV1,
    /// Two-column law only (the single-column law uses one crust density).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oceanic_crust_density_kg_m3: Option<ConstantV1>,
    pub seawater_density_kg_m3: ConstantV1,
    pub continental_crust_thickness_m: ConstantV1,
    pub oceanic_crust_thickness_m: ConstantV1,
    pub stream_power_erodibility_per_yr: ConstantV1,
    pub hillslope_diffusivity_m2_per_yr: ConstantV1,
    // Surfaces
    /// Volume-primary sea-level law only: THE sea-level input.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ocean_water_volume_m3: Option<ConstantV1>,
    /// Prescribed sea level (legacy law) or the output datum the solved sea
    /// level is placed at (volume-primary law).
    pub sea_level_datum_m: ConstantV1,
}

/// What sets sea level, as the generator should read it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SeaLevelInputV1 {
    Prescribed { sea_level_m: f64 },
    OceanVolume { volume_m3: f64, output_datum_m: f64 },
}

/// Law set + constants: what a world pins at creation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UniverseV1 {
    pub schema: Cow<'static, str>,
    pub law_set: LawSetV1,
    pub constants: UniverseConstantsV1,
}

impl UniverseV1 {
    /// Fail-loud validation: schemas, units, envelopes, layering, and the
    /// law/constant combinations this build supports.
    pub fn validate(&self) -> Result<(), UniverseError> {
        let schema = |what, found: &str, expected| {
            if found == expected {
                Ok(())
            } else {
                Err(UniverseError::Schema {
                    what,
                    found: found.to_string(),
                    expected,
                })
            }
        };
        schema("universe", &self.schema, UNIVERSE_SCHEMA_V1)?;
        schema("law set", &self.law_set.schema, LAW_SET_SCHEMA_V1)?;
        schema(
            "universe constants",
            &self.constants.schema,
            UNIVERSE_CONSTANTS_SCHEMA_V1,
        )?;
        let k = &self.constants;
        check("planet_radius_m", &k.planet_radius_m)?;
        check("surface_gravity_m_s2", &k.surface_gravity_m_s2)?;
        check("mantle_density_kg_m3", &k.mantle_density_kg_m3)?;
        check(
            "continental_crust_density_kg_m3",
            &k.continental_crust_density_kg_m3,
        )?;
        if let Some(c) = &k.oceanic_crust_density_kg_m3 {
            check("oceanic_crust_density_kg_m3", c)?;
        }
        check("seawater_density_kg_m3", &k.seawater_density_kg_m3)?;
        check(
            "continental_crust_thickness_m",
            &k.continental_crust_thickness_m,
        )?;
        check("oceanic_crust_thickness_m", &k.oceanic_crust_thickness_m)?;
        check(
            "stream_power_erodibility_per_yr",
            &k.stream_power_erodibility_per_yr,
        )?;
        check(
            "hillslope_diffusivity_m2_per_yr",
            &k.hillslope_diffusivity_m2_per_yr,
        )?;
        if let Some(c) = &k.ocean_water_volume_m3 {
            check("ocean_water_volume_m3", c)?;
        }
        check("sea_level_datum_m", &k.sea_level_datum_m)?;

        // Layering (the same rules PlanetParameterProfile::validate states).
        let rho_m = k.mantle_density_kg_m3.value;
        let rho_cc = k.continental_crust_density_kg_m3.value;
        let rho_w = k.seawater_density_kg_m3.value;
        if rho_cc >= rho_m {
            return Err(UniverseError::Inconsistent(
                "continental_crust_density_kg_m3 must be < mantle_density_kg_m3 (stable layering)",
            ));
        }
        if let Some(c) = &k.oceanic_crust_density_kg_m3 {
            if c.value >= rho_m {
                return Err(UniverseError::Inconsistent(
                    "oceanic_crust_density_kg_m3 must be < mantle_density_kg_m3 (stable layering)",
                ));
            }
            if rho_w >= c.value {
                return Err(UniverseError::Inconsistent(
                    "seawater_density_kg_m3 must be below the oceanic crust density",
                ));
            }
        }
        if rho_w >= rho_cc {
            return Err(UniverseError::Inconsistent(
                "seawater_density_kg_m3 must be below the continental crust density",
            ));
        }
        if k.oceanic_crust_thickness_m.value > k.continental_crust_thickness_m.value {
            return Err(UniverseError::Inconsistent(
                "oceanic_crust_thickness_m must not exceed continental_crust_thickness_m \
                 (the rifted-margin taper thins toward the coast)",
            ));
        }

        // Law entries and the constants each one requires.
        let l = &self.law_set;
        if let Some(fines) = &l.marine_fines {
            fines.validate()?;
            if l.sediment_retention.is_legacy() {
                return Err(UniverseError::Unsupported {
                    law: MarineFinesLawV1::ID,
                    why: "marine fines requires the native continental retention law",
                });
            }
        }
        match &l.isostasy {
            IsostasyLawV1::PrattV0 => {
                return Err(UniverseError::Unsupported {
                    law: "isostasy/pratt-v0",
                    why:
                        "named by ADR 039 §1; no implementation with a declared model scope exists",
                })
            }
            IsostasyLawV1::SingleColumnReferenceV0 {
                reference_crust_thickness_m,
                ..
            } => {
                check("reference_crust_thickness_m", reference_crust_thickness_m)?;
                if l.sea_level != SeaLevelLawV1::PrescribedDatumV0 {
                    return Err(UniverseError::Unsupported {
                        law: "isostasy/single-column-reference-v0 + sea-level/volume-primary-global-fill-v0",
                        why: "the single-column law water-loads against a prescribed datum; a volume-primary sea level is only implemented for the two-column law",
                    });
                }
                if k.oceanic_crust_density_kg_m3.is_some() {
                    return Err(UniverseError::Inconsistent(
                        "oceanic_crust_density_kg_m3 is set but the single-column law uses one crust density",
                    ));
                }
            }
            IsostasyLawV1::TwoColumnCommonDepthV0 {
                fixed_compensation_depth_m,
                derived_depth_margin_m,
            } => {
                check("fixed_compensation_depth_m", fixed_compensation_depth_m)?;
                if let Some(m) = derived_depth_margin_m {
                    check("derived_depth_margin_m", m)?;
                }
                if l.sea_level != SeaLevelLawV1::VolumePrimaryGlobalFillV0 {
                    return Err(UniverseError::Unsupported {
                        law: "isostasy/two-column-common-depth-v0 + sea-level/prescribed-datum-v0",
                        why: "the two-column law solves sea level from the ocean volume; a prescribed sea level with it is not implemented",
                    });
                }
                if k.oceanic_crust_density_kg_m3.is_none() {
                    return Err(UniverseError::Inconsistent(
                        "the two-column law requires oceanic_crust_density_kg_m3",
                    ));
                }
            }
        }
        match l.sea_level {
            SeaLevelLawV1::PrescribedDatumV0 => {
                if k.ocean_water_volume_m3.is_some() {
                    return Err(UniverseError::Inconsistent(
                        "a prescribed sea level and an ocean water volume are two inputs for one \
                         quantity (ADR 039 §3); the prescribed-datum law takes no volume",
                    ));
                }
            }
            SeaLevelLawV1::VolumePrimaryGlobalFillV0 => {
                if k.ocean_water_volume_m3.is_none() {
                    return Err(UniverseError::Inconsistent(
                        "the volume-primary sea-level law requires ocean_water_volume_m3",
                    ));
                }
            }
        }
        if let CrustSeedingLawV1::GeodesicDiscBodiesV0 {
            continental_area_fraction,
            rifted_margin_width_m,
        } = &l.crust_seeding
        {
            check("continental_area_fraction", continental_area_fraction)?;
            check("rifted_margin_width_m", rifted_margin_width_m)?;
        }
        if let TectonicThickeningLawV1::RateLimitedBeltV1 {
            belt_shear_zone_dip_deg,
            lower_crust_channel_viscosity_pa_s,
            channel_reference_thickness_m,
            delamination_thickness_m,
            arc_productivity_m,
            arc_belt_width_m,
        } = &l.tectonic_thickening
        {
            check("belt_shear_zone_dip_deg", belt_shear_zone_dip_deg)?;
            check(
                "lower_crust_channel_viscosity_pa_s",
                lower_crust_channel_viscosity_pa_s,
            )?;
            check("channel_reference_thickness_m", channel_reference_thickness_m)?;
            check("delamination_thickness_m", delamination_thickness_m)?;
            check("arc_productivity_m", arc_productivity_m)?;
            check("arc_belt_width_m", arc_belt_width_m)?;
            if delamination_thickness_m.value <= channel_reference_thickness_m.value {
                return Err(UniverseError::Inconsistent(
                    "delamination_thickness_m must exceed channel_reference_thickness_m \
                     (the backstop acts on crust already thick enough to flow)",
                ));
            }
            if delamination_thickness_m.value <= k.continental_crust_thickness_m.value {
                return Err(UniverseError::Inconsistent(
                    "delamination_thickness_m must exceed the continental crust thickness",
                ));
            }
        }
        if !l.plate_carry.is_legacy() && l.tectonic_thickening.is_legacy() {
            return Err(UniverseError::Unsupported {
                law: "plate-carry/front-rescue-carry-v1 with tectonic-thickening/kinematic-overlap-stacking-v0",
                why: "front-rescue-carry-v1 places its consumed collision volume on the rate-limited-belt-v1 belt; a combination with the legacy stacking rule is not implemented",
            });
        }
        // SEED-CARRY-CONSERVE (lane CV): the parcel carry's constants (its
        // belt requirement is the check above).
        if let PlateCarryLawV1::ContinentalParcelCarryV1 {
            parcel_overlap_cell_fraction,
            plate_fragment_min_area_fraction,
            write_back_residual_iterations,
        } = &l.plate_carry
        {
            check("parcel_overlap_cell_fraction", parcel_overlap_cell_fraction)?;
            check("plate_fragment_min_area_fraction", plate_fragment_min_area_fraction)?;
            if let Some(n) = write_back_residual_iterations {
                check("write_back_residual_iterations", n)?;
                if n.value.fract() != 0.0 {
                    return Err(UniverseError::Inconsistent("continental-parcel-carry-v1: write_back_residual_iterations must be a whole number"));
                }
            }
        }
        if let PlateSutureLawV1::ContactLengthPoleShareV1 {
            suture_contact_length_m,
        } = &l.plate_suture
        {
            check("suture_contact_length_m", suture_contact_length_m)?;
        }
        // SEED-RANGES (lane RG): the indentation suture law and the plateau
        // law.
        if let PlateSutureLawV1::ContactLengthIndentationTaperV1 {
            suture_contact_length_m,
            indentation_time_my,
            indentation_shortening_budget_m,
        } = &l.plate_suture
        {
            check("suture_contact_length_m", suture_contact_length_m)?;
            check("indentation_time_my", indentation_time_my)?;
            check("indentation_shortening_budget_m", indentation_shortening_budget_m)?;
            if !(indentation_time_my.value > 0.0) || !(indentation_shortening_budget_m.value > 0.0) {
                return Err(UniverseError::Inconsistent(
                    "contact-length-indentation-taper-v1: indentation time and shortening budget must be positive",
                ));
            }
        }
        // SEED-WORLD-COMPOSE (lane CP): v2's extra constants, then the
        // shared v1 checks on its shared constants.
        if let OrogenPlateauLawV1::PlateauCapSpillV2 {
            plateau_foothill_gradient,
            plateau_taper_floor_m,
            ..
        } = &l.orogen_plateau
        {
            check("plateau_foothill_gradient", plateau_foothill_gradient)?;
            check("plateau_taper_floor_m", plateau_taper_floor_m)?;
        }
        let plateau_view = l.orogen_plateau.shared_view();
        if let OrogenPlateauLawV1::PlateauCapSpillV1 {
            plateau_crust_cap_m,
            plateau_spill_reach_m,
        } = &plateau_view
        {
            check("plateau_crust_cap_m", plateau_crust_cap_m)?;
            check("plateau_spill_reach_m", plateau_spill_reach_m)?;
            match &l.tectonic_thickening {
                TectonicThickeningLawV1::RateLimitedBeltV1 {
                    delamination_thickness_m,
                    ..
                } => {
                    if !(plateau_crust_cap_m.value > k.continental_crust_thickness_m.value)
                        || !(delamination_thickness_m.value > plateau_crust_cap_m.value)
                    {
                        return Err(UniverseError::Inconsistent(
                            "plateau-cap-spill-v1: the cap must exceed the continental crust thickness and the belt's delamination backstop must exceed the cap",
                        ));
                    }
                }
                _ => {
                    return Err(UniverseError::Unsupported {
                        law: "orogen-plateau/plateau-cap-spill-v1 without rate-limited-belt-v1",
                        why: "the plateau law acts on the rate-limited belt's thickened crust",
                    })
                }
            }
        }
        // SEED-MAGMATISM (lane MG): slab-arc-backarc-llsvp-v1 constants and
        // requirements.
        // SEED-WORLD-COMPOSE (lane CP): v2's extra constants, then the
        // shared v1 checks.
        if let MagmatismLawV1::SlabArcBackarcLlsvpV2 {
            backarc_episode_my,
            backarc_max_stretch_factor,
            kimberlite_plume_radius_m,
            ..
        } | MagmatismLawV1::SlabArcBackarcLlsvpV3 {
            backarc_episode_my, backarc_max_stretch_factor, kimberlite_plume_radius_m, ..
        } = &l.magmatism
        {
            check("backarc_episode_my", backarc_episode_my)?;
            check("backarc_max_stretch_factor", backarc_max_stretch_factor)?;
            check("kimberlite_plume_radius_m", kimberlite_plume_radius_m)?;
        }
        let magmatism_view = l.magmatism.shared_view();
        if let MagmatismLawV1::SlabArcBackarcLlsvpV1 {
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
        } = &magmatism_view
        {
            check("slab_depth_below_arc_m", slab_depth_below_arc_m)?;
            check("slab_dip_young_deg", slab_dip_young_deg)?;
            check("slab_dip_old_deg", slab_dip_old_deg)?;
            check("slab_dip_age_my", slab_dip_age_my)?;
            check("arc_magma_productivity_m", arc_magma_productivity_m)?;
            check("arc_rate_floor_m_per_my", arc_rate_floor_m_per_my)?;
            check("arc_rate_cap_m_per_my", arc_rate_cap_m_per_my)?;
            check("arc_half_width_m", arc_half_width_m)?;
            check("arc_continental_thickness_m", arc_continental_thickness_m)?;
            check("arc_thermal_uplift_m", arc_thermal_uplift_m)?;
            check("arc_thermal_half_width_m", arc_thermal_half_width_m)?;
            check("arc_thermal_decay_my", arc_thermal_decay_my)?;
            check("backarc_slab_age_min_my", backarc_slab_age_min_my)?;
            check("backarc_overriding_advance_max_m_per_my", backarc_overriding_advance_max_m_per_my)?;
            check("backarc_extension_rate_m_per_my", backarc_extension_rate_m_per_my)?;
            check("backarc_axis_offset_m", backarc_axis_offset_m)?;
            check("backarc_half_width_m", backarc_half_width_m)?;
            check("backarc_min_crust_m", backarc_min_crust_m)?;
            check("llsvp_edge_deg", llsvp_edge_deg)?;
            check("llsvp_edge_sigma_deg", llsvp_edge_sigma_deg)?;
            check("llsvp_antipode_jitter_deg", llsvp_antipode_jitter_deg)?;
            check("plume_background_fraction", plume_background_fraction)?;
            check("plume_drift_rate_m_per_my", plume_drift_rate_m_per_my)?;
            check("vent_active_window_my", vent_active_window_my)?;
            check("volcano_spacing_m", volcano_spacing_m)?;
            check("arc_crust_cap_m", arc_crust_cap_m)?;
            if !(slab_dip_old_deg.value >= slab_dip_young_deg.value)
                || !(arc_rate_cap_m_per_my.value >= arc_rate_floor_m_per_my.value)
            {
                return Err(UniverseError::Inconsistent(
                    "slab-arc-backarc-llsvp-v1: the old-slab dip must be at least the young-slab dip and the rate cap at least the floor",
                ));
            }
            if !matches!(l.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. })
                || l.tectonic_thickening.is_legacy()
                || (l.plate_suture.is_legacy()
                    && !matches!(l.plate_carry, PlateCarryLawV1::ContinentalParcelCarryV1 { .. }))
            {
                return Err(UniverseError::Unsupported {
                    law: "magmatism/slab-arc-backarc-llsvp-v1 without two-column isostasy, rate-limited-belt-v1 and the parcel carry or a suture law",
                    why: "its arc heat is a column term of the two-column isostasy, it lives in the belt block and its run state is checkpointed with the parcel/suture state",
                });
            }
        }
        // SEED-IRREGULAR-SEEDS (lane IS): craton constants; the law refines
        // the bodies law (area fraction and margin width come from it).
        if let CratonSeedingLawV1::IrregularCratonsV1 {
            craton_interior_thickness_m,
            craton_rank_size_exponent,
            cratons_per_unit_share,
            basement_roughness_std_m,
            suture_crest_thickening_m,
            suture_half_width_m,
            proximal_margin_stretch,
        } = &l.craton_seeding
        {
            check("craton_interior_thickness_m", craton_interior_thickness_m)?;
            check("craton_rank_size_exponent", craton_rank_size_exponent)?;
            check("cratons_per_unit_share", cratons_per_unit_share)?;
            check("basement_roughness_std_m", basement_roughness_std_m)?;
            check("suture_crest_thickening_m", suture_crest_thickening_m)?;
            check("suture_half_width_m", suture_half_width_m)?;
            check("proximal_margin_stretch", proximal_margin_stretch)?;
            if !matches!(l.crust_seeding, CrustSeedingLawV1::GeodesicDiscBodiesV0 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "craton-seeding/irregular-cratons-v1 without geodesic-disc-bodies-v0",
                    why: "the craton law takes its continental area fraction and rifted margin width from the bodies law",
                });
            }
        }
        if !l.collision_sediment.is_legacy() && l.tectonic_thickening.is_legacy() {
            return Err(UniverseError::Unsupported {
                law: "collision-sediment/accreted-to-belt-crust-v1 with tectonic-thickening/kinematic-overlap-stacking-v0",
                why: "accreted-to-belt-crust-v1 adds collision sediment to the rate-limited-belt-v1 crust supply; a combination with the legacy stacking rule is not implemented",
            });
        }
        // SEED-EXTRAS (lane EX): plume, rift and plate-split laws.
        if let HotspotVolcanismLawV1::FixedPlumeChainsV1 {
            plume_count,
            plume_onset_lead_my,
            plume_head_volume_m3,
            plume_head_radius_m,
            plume_tail_flux_m3_per_my,
            plume_edifice_radius_m,
        } = &l.hotspot_volcanism
        {
            check("plume_count", plume_count)?;
            if plume_count.value.fract() != 0.0 {
                return Err(UniverseError::Inconsistent(
                    "plume_count must be a whole number",
                ));
            }
            check("plume_onset_lead_my", plume_onset_lead_my)?;
            check("plume_head_volume_m3", plume_head_volume_m3)?;
            check("plume_head_radius_m", plume_head_radius_m)?;
            check("plume_tail_flux_m3_per_my", plume_tail_flux_m3_per_my)?;
            check("plume_edifice_radius_m", plume_edifice_radius_m)?;
        }
        if let ContinentalRiftLawV1::StretchedRiftV1 {
            rift_half_width_m,
            breakup_crust_thickness_m,
        } = &l.continental_rift
        {
            check("rift_half_width_m", rift_half_width_m)?;
            check("breakup_crust_thickness_m", breakup_crust_thickness_m)?;
            if breakup_crust_thickness_m.value >= k.continental_crust_thickness_m.value {
                return Err(UniverseError::Inconsistent(
                    "breakup_crust_thickness_m must be below the continental crust thickness",
                ));
            }
        }
        if let PlateReorganizationLawV1::PlumeHeadContinentalSplitV1 {
            rift_opening_rate_deg_per_my,
        } = &l.plate_reorganization
        {
            check("rift_opening_rate_deg_per_my", rift_opening_rate_deg_per_my)?;
            if l.hotspot_volcanism.is_legacy() || l.continental_rift.is_legacy() {
                return Err(UniverseError::Unsupported {
                    law: "plate-reorganization/plume-head-continental-split-v1 without fixed-plume-chains-v1 and stretched-rift-v1",
                    why: "the split is triggered by a plume head and its opening is filled by the rift law",
                });
            }
        }
        // SEED-PLATE-REORG (lane PR): the scripted reorganization.
        // SEED-WORLD-COMPOSE (lane CP): v3's initiation age, then the
        // shared checks on v1's constants.
        if let PlateReorganizationLawV1::PlateReorgScriptedV3 { subduction_initiation_age_my, .. } = &l.plate_reorganization {
            check("subduction_initiation_age_my", subduction_initiation_age_my)?;
        }
        // SEED-TRANSFORMS (lane TF): v4's initiation age and limits.
        if let PlateReorganizationLawV1::PlateReorgScriptedV4 {
            subduction_initiation_age_my,
            initiation_slab_age_min_my,
            initiation_old_slab_fraction_min,
            initiation_trench_spacing_m,
            ..
        } = &l.plate_reorganization
        {
            check("subduction_initiation_age_my", subduction_initiation_age_my)?;
            check("initiation_slab_age_min_my", initiation_slab_age_min_my)?;
            check("initiation_old_slab_fraction_min", initiation_old_slab_fraction_min)?;
            check("initiation_trench_spacing_m", initiation_trench_spacing_m)?;
        }
        // SEED-OLD-OCEAN (lane OR): v5's initiation age, v4 limits and the
        // margin initiation.
        if let PlateReorganizationLawV1::PlateReorgScriptedV5 {
            subduction_initiation_age_my,
            initiation_slab_age_min_my,
            initiation_old_slab_fraction_min,
            initiation_trench_spacing_m,
            margin_initiation_age_my,
            margin_segment_length_m,
            margin_trench_spacing_m,
            ..
        } = &l.plate_reorganization
        {
            check("subduction_initiation_age_my", subduction_initiation_age_my)?;
            check("initiation_slab_age_min_my", initiation_slab_age_min_my)?;
            check("initiation_old_slab_fraction_min", initiation_old_slab_fraction_min)?;
            check("initiation_trench_spacing_m", initiation_trench_spacing_m)?;
            check("margin_initiation_age_my", margin_initiation_age_my)?;
            check("margin_segment_length_m", margin_segment_length_m)?;
            check("margin_trench_spacing_m", margin_trench_spacing_m)?;
        }
        // SEED-TRANSFORMS (lane TF): the boundary law runs inside the parcel
        // carry's advance.
        if let PlateBoundaryLawV1::SmoothedNormalTransformV1 {
            transform_obliquity_deg,
            boundary_normal_smoothing_m,
        } = &l.plate_boundary
        {
            check("transform_obliquity_deg", transform_obliquity_deg)?;
            check("boundary_normal_smoothing_m", boundary_normal_smoothing_m)?;
            if !matches!(l.plate_carry, PlateCarryLawV1::ContinentalParcelCarryV1 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "plate-boundary/smoothed-normal-transform-v1 without continental-parcel-carry-v1",
                    why: "the transform gap fill and overlap rule run in the parcel carry's advance (advance_with_parcels_transform)",
                });
            }
        }
        let reorg_view = l.plate_reorganization.shared_view();
        // SEED-BREAKUP-CUT-V2 (lane BC): v2 has v1's constants and checks.
        if let PlateReorganizationLawV1::PlateReorgScriptedV1 {
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
        } = &reorg_view
        {
            check("rift_opening_rate_deg_per_my", rift_opening_rate_deg_per_my)?;
            check("reorg_interval_min_my", reorg_interval_min_my)?;
            check("reorg_interval_max_my", reorg_interval_max_my)?;
            check("slab_pull_bias", slab_pull_bias)?;
            check("slab_circumference_reference", slab_circumference_reference)?;
            check("reorg_rate_floor_fraction", reorg_rate_floor_fraction)?;
            check("supercontinent_share", supercontinent_share)?;
            check("supercontinent_tenure_my", supercontinent_tenure_my)?;
            check("breakup_piece_max_share", breakup_piece_max_share)?;
            check("breakup_cascade_my", breakup_cascade_my)?;
            check("rift_persistence_my", rift_persistence_my)?;
            check("weak_zone_strength_ratio", weak_zone_strength_ratio)?;
            check("plate_count_min", plate_count_min)?;
            check("plate_count_max", plate_count_max)?;
            check("plate_count_min_area_fraction", plate_count_min_area_fraction)?;
            check("subduction_initiation_rate_deg_per_my", subduction_initiation_rate_deg_per_my)?;
            check("plume_head_interval_my", plume_head_interval_my)?;
            check("plume_tail_lifetime_my", plume_tail_lifetime_my)?;
            check("plume_split_min_share", plume_split_min_share)?;
            if reorg_interval_min_my.value > reorg_interval_max_my.value
                || plate_count_min.value > plate_count_max.value
                || plate_count_min.value.fract() != 0.0
                || plate_count_max.value.fract() != 0.0
            {
                return Err(UniverseError::Inconsistent(
                    "plate-reorg-scripted-v1: interval min must not exceed max, plate counts must be whole and min <= max",
                ));
            }
            if l.hotspot_volcanism.is_legacy()
                || l.continental_rift.is_legacy()
                || !matches!(l.plate_carry, PlateCarryLawV1::ContinentalParcelCarryV1 { .. })
                || l.plate_suture.is_legacy()
            {
                return Err(UniverseError::Unsupported {
                    law: "plate-reorganization/plate-reorg-scripted-v1 without fixed-plume-chains-v1, stretched-rift-v1, continental-parcel-carry-v1 and contact-length-pole-share-v1",
                    why: "its plume splits need the plume law, its openings the rift law, its cuts move continental parcels, and its groups are suture groups",
                });
            }
        }
        // SEED-RIFT-SHOULDERS (lane RS): the rift flank law.
        if let RiftFlankLawV1::ShoulderUpliftThermalV1 {
            rift_flank_flexural_parameter_m,
            rift_shoulder_ratio,
        } = &l.rift_flank
        {
            check("rift_flank_flexural_parameter_m", rift_flank_flexural_parameter_m)?;
            check("rift_shoulder_ratio", rift_shoulder_ratio)?;
            if l.continental_rift.is_legacy() {
                return Err(UniverseError::Unsupported {
                    law: "rift-flank/shoulder-uplift-thermal-v1 without stretched-rift-v1",
                    why: "the stretching factor and thermal age come from the stretched-rift-v1 births and flank thinning",
                });
            }
            if !matches!(l.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "rift-flank/shoulder-uplift-thermal-v1 with the single-column isostasy law",
                    why: "the rift thermal and shoulder term is a water-loaded column term of the two-column isostasy law",
                });
            }
        }
        // SEED-TRENCH part 2 (lane TN): the seeded passive-margin law.
        if !l.passive_margin.is_legacy() {
            if l.rift_flank.is_legacy() {
                return Err(UniverseError::Unsupported {
                    law: "passive-margin/seeded-margin-thermal-subsidence-v1 without shoulder-uplift-thermal-v1",
                    why: "the seeded margins' thermal subsidence is evaluated by the rift-flank law from its stretching markers",
                });
            }
            if !matches!(l.crust_seeding, CrustSeedingLawV1::GeodesicDiscBodiesV0 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "passive-margin/seeded-margin-thermal-subsidence-v1 without geodesic-disc-bodies-v0",
                    why: "the stretching factor comes from the seeded rifted-margin taper of the continental bodies",
                });
            }
        }
        // SEED-STRATIGRAPHY (lane ST): the parcel stacks.
        if !l.stratigraphy.is_legacy() {
            for (id, c) in l.stratigraphy.constants() {
                check(id, c)?;
            }
            if !matches!(l.plate_carry, PlateCarryLawV1::ContinentalParcelCarryV1 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "stratigraphy/parcel-stacks-v1 without continental-parcel-carry-v1",
                    why: "the stacks live on the continental parcels",
                });
            }
        }
        // SEED-FLOODBASALT (lane FB): the LIP and swell law.
        if let FloodBasaltLawV1::PlumeHeadLipSwellV1 {
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
        } = &l.flood_basalt
        {
            check("flood_basalt_extrusive_volume_m3", flood_basalt_extrusive_volume_m3)?;
            check("flood_basalt_intrusive_ratio", flood_basalt_intrusive_ratio)?;
            check("flood_basalt_radius_m", flood_basalt_radius_m)?;
            check("flood_basalt_taper_fraction", flood_basalt_taper_fraction)?;
            check("flood_basalt_eruption_duration_my", flood_basalt_eruption_duration_my)?;
            check("swell_tail_amplitude_m", swell_tail_amplitude_m)?;
            check("swell_tail_radius_m", swell_tail_radius_m)?;
            check("swell_head_amplitude_m", swell_head_amplitude_m)?;
            check("swell_head_radius_m", swell_head_radius_m)?;
            check("swell_head_decay_my", swell_head_decay_my)?;
            if l.hotspot_volcanism.is_legacy() {
                return Err(UniverseError::Unsupported {
                    law: "flood-basalt/plume-head-lip-swell-v1 without fixed-plume-chains-v1",
                    why: "the provinces and swells belong to the seeded plumes and their onsets",
                });
            }
            if !matches!(l.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "flood-basalt/plume-head-lip-swell-v1 with the single-column isostasy law",
                    why: "the plume swell is a water-loaded column term of the two-column isostasy law",
                });
            }
        }
        // SEED-TRENCH (lane TN): the subduction trench law.
        if let SubductionTrenchLawV1::SlabFlexureTrenchV1 {
            trench_vertical_load_fraction,
            trench_slab_length_max_m,
            trench_thermal_expansivity_per_k,
            trench_mantle_temperature_c,
            trench_elastic_isotherm_c,
            trench_thermal_diffusivity_m2_per_yr,
            trench_elastic_modulus_pa,
            trench_poisson_ratio,
        } = &l.subduction_trench
        {
            check("trench_vertical_load_fraction", trench_vertical_load_fraction)?;
            check("trench_slab_length_max_m", trench_slab_length_max_m)?;
            check("trench_thermal_expansivity_per_k", trench_thermal_expansivity_per_k)?;
            check("trench_mantle_temperature_c", trench_mantle_temperature_c)?;
            check("trench_elastic_isotherm_c", trench_elastic_isotherm_c)?;
            check("trench_thermal_diffusivity_m2_per_yr", trench_thermal_diffusivity_m2_per_yr)?;
            check("trench_elastic_modulus_pa", trench_elastic_modulus_pa)?;
            check("trench_poisson_ratio", trench_poisson_ratio)?;
            if trench_elastic_isotherm_c.value >= trench_mantle_temperature_c.value {
                return Err(UniverseError::Inconsistent(
                    "trench_elastic_isotherm_c must be below trench_mantle_temperature_c",
                ));
            }
            if !matches!(l.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "subduction-trench/slab-flexure-trench-v1 with the single-column isostasy law",
                    why: "the trench deflection is a water-loaded column term of the two-column isostasy law",
                });
            }
        }
        // SEED-UPLANDS (lane UP): the threshold hillslope law.
        if let HillslopeThresholdLawV1::CriticalSlopeEnvelopeV1 {
            hillslope_critical_slope,
            knickpoint_celerity_m_per_my,
        } = &l.hillslope_threshold
        {
            check("hillslope_critical_slope", hillslope_critical_slope)?;
            check("knickpoint_celerity_m_per_my", knickpoint_celerity_m_per_my)?;
            if !matches!(
                l.erosion,
                ErosionLawV1::StreamPowerSubgridV1 { .. }
                    | ErosionLawV1::StreamPowerSubgridNonlinearV1 { .. }
            ) {
                return Err(UniverseError::Unsupported {
                    law: "hillslope-threshold/critical-slope-envelope-v1 without a subgrid erosion law",
                    why: "the threshold target is realised through the subgrid cell-mean target and effective areas",
                });
            }
        }
        // INCISION-THRESHOLD (lane IT): the excess-steepness threshold is
        // realised inside the nonlinear subgrid step.
        if let FluvialIncisionThresholdLawV1::ExcessSteepnessThresholdV1 {
            critical_channel_steepness_m,
        }
        | FluvialIncisionThresholdLawV1::ExcessSteepnessThresholdTrunkOnlyV1 {
            critical_channel_steepness_m,
        } = &l.fluvial_incision_threshold
        {
            check("critical_channel_steepness_m", critical_channel_steepness_m)?;
            if !matches!(l.erosion, ErosionLawV1::StreamPowerSubgridNonlinearV1 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "fluvial-incision-threshold/excess-steepness-threshold-v1 without stream-power-subgrid-nonlinear-v1",
                    why: "the threshold is implemented on the nonlinear subgrid trunk and cell-mean closure only",
                });
            }
            // Merge rule (lane CLIMATE-GENESIS-2): the threshold step has no
            // per-cell runoff input, so I11 would be silently dropped.
            if let ClimateForcingLawV1::EpochClimateRunoffV1 { subgrid_cell_mean, .. } = &l.climate_forcing {
                if !subgrid_cell_mean.is_default() {
                    return Err(UniverseError::Unsupported {
                        law: "fluvial-incision-threshold/excess-steepness-threshold-v1 with climate-forcing subgrid_cell_mean = cell-runoff-subgrid-erosivity-v1",
                        why: "the threshold subgrid step does not take the per-cell runoff ratio yet",
                    });
                }
            }
        }
        // GLACIAL-T0 (lane GLACIAL-T0): the glacial factors act inside the
        // nonlinear subgrid step and read the genesis climate fields.
        if let GlacialErosionLawV1::GlacialStreamPowerBlendV1 {
            glacial_to_fluvial_erodibility_ratio,
            full_ice_height_above_ela_m,
            coldhouse_duty,
            glacial_max_cooling_k,
            glacial_max_precipitation_factor,
            warm_state_duties,
            warm_state_warmings_k,
        } = &l.glacial_erosion
        {
            check("glacial_to_fluvial_erodibility_ratio", glacial_to_fluvial_erodibility_ratio)?;
            check("full_ice_height_above_ela_m", full_ice_height_above_ela_m)?;
            check("glacial_state_duty", coldhouse_duty)?;
            check("glacial_state_offset_k", glacial_max_cooling_k)?;
            check("glacial_max_precipitation_factor", glacial_max_precipitation_factor)?;
            for d in warm_state_duties {
                check("glacial_state_duty", d)?;
            }
            for t in warm_state_warmings_k {
                check("glacial_state_offset_k", t)?;
            }
            let total = coldhouse_duty.value + warm_state_duties.iter().map(|d| d.value).sum::<f64>();
            if warm_state_duties.len() != warm_state_warmings_k.len() || total > 1.0 + 1e-12 {
                return Err(UniverseError::Unsupported {
                    law: "glacial-erosion/glacial-stream-power-blend-v1 with unpaired warm states or duties summing above 1",
                    why: "each warm state needs one duty and one warming, and the states share one timeline",
                });
            }
            if !matches!(l.erosion, ErosionLawV1::StreamPowerSubgridNonlinearV1 { .. }) {
                return Err(UniverseError::Unsupported {
                    law: "glacial-erosion/glacial-stream-power-blend-v1 without stream-power-subgrid-nonlinear-v1",
                    why: "the glacial factors act on the nonlinear subgrid trunk and cell-mean closure only",
                });
            }
            match &l.climate_forcing {
                ClimateForcingLawV1::EpochClimateRunoffV1 { subgrid_cell_mean, .. } => {
                    if !subgrid_cell_mean.is_default() {
                        return Err(UniverseError::Unsupported {
                            law: "glacial-erosion/glacial-stream-power-blend-v1 with climate-forcing subgrid_cell_mean = cell-runoff-subgrid-erosivity-v1",
                            why: "the glacial subgrid step does not take the per-cell runoff ratio",
                        });
                    }
                }
                _ => {
                    return Err(UniverseError::Unsupported {
                        law: "glacial-erosion/glacial-stream-power-blend-v1 without epoch-climate-runoff-v1",
                        why: "the equilibrium-line altitude needs the genesis climate's temperature and precipitation fields",
                    });
                }
            }
        }
        if !l.marine_remainder_layer.is_legacy()
            && l.sediment_retention.is_legacy()
            && !(l.land_deposition.is_legacy() && l.sediment_placement.is_legacy())
        {
            return Err(UniverseError::Unsupported {
                law: "marine-remainder-layer/once-per-step-capped-global-layer-v1 with erosion-deposition or shelf-wave-base",
                why: "wired into the legacy step and continental-retention-base-level-v1 only",
            });
        }
        // SHELF (law-set entry, lane CLIMATE-GENESIS-2): grading reuses the
        // retention law's wave base and Airy sediment column.
        // INLAND-BASINS: the lake placement lives in the retention step.
        if !l.closed_basin.is_legacy() && l.sediment_retention.is_legacy() {
            return Err(UniverseError::Inconsistent(
                "closed-basin/closed-basin-water-balance-v1 requires sediment-retention/continental-retention-base-level-v1 (its lake placement is part of that step)",
            ));
        }
        // ALLUVIAL-TRUNK: the release reads k_c (threshold law) and S_a
        // (retention law).
        if !l.fluvial_alluvial_channel.is_legacy() && (l.fluvial_incision_threshold.is_legacy() || l.sediment_retention.is_legacy()) {
            return Err(UniverseError::Inconsistent(
                "fluvial-alluvial-channel/stream-power-threshold-free-alluvial-channels-v1 requires a fluvial incision threshold (its k_c) and sediment-retention/continental-retention-base-level-v1 (its graded alluvial slope S_a)",
            ));
        }
        // LAKE-OUTLET: the outlet incision reads the water balance.
        // LAKE-HOLD-2: the path grade acts on the overflow path and reads k_c.
        if !l.lake_outlet_grade.is_legacy() && (!l.lake_outlet.incises_path() || l.fluvial_incision_threshold.is_legacy()) {
            return Err(UniverseError::Inconsistent(
                "lake-outlet-grade/lake-outlet-path-threshold-grade-v1 requires a path-incising lake-outlet law and a fluvial incision threshold (its k_c)",
            ));
        }
        if !l.lake_outlet.is_legacy() && l.closed_basin.is_legacy() {
            return Err(UniverseError::Inconsistent(
                "lake-outlet/lake-outlet-sill-incision-v1 requires closed-basin/closed-basin-water-balance-v1 (it incises with that balance's overflow)",
            ));
        }
        // PLATE-POLES: the rigid-group classification lives in the
        // smoothed-normal classifier.
        if !l.rigid_group.is_legacy() && l.plate_boundary.is_legacy() {
            return Err(UniverseError::Inconsistent(
                "rigid-group/rigid-group-interior-v1 and rigid-group-contacts-v1 require plate-boundary/smoothed-normal-transform-v1 (the classification is part of that classifier)",
            ));
        }
        if !l.shelf_grading.is_legacy() && l.sediment_retention.is_legacy() {
            return Err(UniverseError::Inconsistent(
                "shelf-grading/wave-base-shelf-grading-v1 requires sediment-retention/continental-retention-base-level-v1 (it reuses its wave base and sediment column)",
            ));
        }
        // CLIMATE-GENESIS (lane CG): the climate forcing law.
        if let ClimateForcingLawV1::EpochClimateRunoffV1 {
            climate,
            seed_constants_digest,
            rotation_constants_digest,
            climate_epoch_my,
            climate_solve_grid_divisor,
            climate_sweep_rel_tol,
            climate_runoff_floor_fraction,
            subgrid_cell_mean,
            erosive_runoff,
            ..
        } = &l.climate_forcing
        {
            if let Some(e) = erosive_runoff {
                e.validate()?;
            }
            if !subgrid_cell_mean.is_default() && matches!(l.erosion, ErosionLawV1::StreamPowerN1LinearDiffusionV0) {
                return Err(UniverseError::Inconsistent(
                    "epoch-climate-runoff-v1: subgrid_cell_mean = cell-runoff-subgrid-erosivity-v1 needs a subgrid erosion law (there is no cell-mean relaxation to force)",
                ));
            }
            check("climate_epoch_my", climate_epoch_my)?;
            check("climate_solve_grid_divisor", climate_solve_grid_divisor)?;
            check("climate_sweep_rel_tol", climate_sweep_rel_tol)?;
            check("climate_runoff_floor_fraction", climate_runoff_floor_fraction)?;
            if ![1.0, 2.0, 4.0].contains(&climate_solve_grid_divisor.value) {
                return Err(UniverseError::Inconsistent(
                    "climate_solve_grid_divisor must be 1, 2 or 4 (whole 2x2 / 4x4 blocks)",
                ));
            }
            climate.validate()?;
            match (climate.law, rotation_constants_digest.is_some()) {
                (ClimateLawV1::SeedClimateV1, false) | (ClimateLawV1::SeedClimateV2, true) => {}
                (ClimateLawV1::SeedClimateV1, true) | (ClimateLawV1::SeedClimateV2, false) => {
                    return Err(UniverseError::Inconsistent(
                        "epoch-climate-runoff-v1: rotation_constants_digest is required with seed-climate-v2 and refused with seed-climate-v1",
                    ));
                }
                (ClimateLawV1::ClimateCoarseV0, _) => {
                    return Err(UniverseError::Unsupported {
                        law: "climate-forcing/epoch-climate-runoff-v1 with climate-coarse-v0",
                        why: "the genesis forcing needs the 2-D column state of seed-climate-v1 or -v2",
                    });
                }
            }
            if !seed_constants_digest.starts_with("blake3:") {
                return Err(UniverseError::Inconsistent(
                    "epoch-climate-runoff-v1: seed_constants_digest must be a blake3: digest",
                ));
            }
            if self.constants.stream_power_erodibility_per_yr.value < CLIMATE_FORCING_MIN_REFERENCE_K_PER_YR {
                return Err(UniverseError::Inconsistent(
                    "epoch-climate-runoff-v1 quotes K at the reference runoff; a K below 1.5e-6 /yr is the uniform-runoff-equivalent recipe (RELIEF-DECAY K 1.0e-6) and would apply the runoff factor twice (alternatives, never stacked)",
                ));
            }
        }
        if (!l.hotspot_volcanism.is_legacy() || !l.continental_rift.is_legacy())
            && l.tectonic_thickening.is_legacy()
        {
            return Err(UniverseError::Unsupported {
                law: "hotspot-volcanism / continental-rift with tectonic-thickening/kinematic-overlap-stacking-v0",
                why: "the plume and rift passes are wired through the rate-limited-belt-v1 advance and its run state; a combination with the legacy stacking rule is not implemented",
            });
        }
        if let ErosionLawV1::StreamPowerSubgridV1 {
            channel_head_length_m,
        } = &l.erosion
        {
            check("subgrid_channel_head_length_m", channel_head_length_m)?;
        }
        if let ErosionLawV1::StreamPowerSubgridNonlinearV1 {
            channel_head_length_m,
            slope_exponent_n,
            reference_steepness_m,
        } = &l.erosion
        {
            check("subgrid_channel_head_length_m", channel_head_length_m)?;
            check("stream_power_slope_exponent_n", slope_exponent_n)?;
            check("stream_power_reference_steepness_m", reference_steepness_m)?;
        }
        if let LandDepositionLawV1::ErosionDepositionV1 {
            deposition_coefficient_g,
        } = &l.land_deposition
        {
            check("deposition_coefficient_g", deposition_coefficient_g)?;
            if !l.erosion.is_legacy() {
                return Err(UniverseError::Unsupported {
                    law: "land-deposition/erosion-deposition-v1 with a non-legacy erosion law",
                    why: "erosion-deposition-v1 is implemented on the single-elevation implicit stream power only (the subgrid law rewrites contributing areas, which the deposition term also reads)",
                });
            }
        }
        // SEED-SEDIMENT-RETAIN (lane SR): continental retention + taper.
        if let SedimentRetentionLawV1::ContinentalRetentionBaseLevelV1 {
            base_level_relief_scale_m,
            base_level_taper_exponent,
            alluvial_grade_slope,
            shelf_wave_base_depth_m,
        } = &l.sediment_retention
        {
            check("base_level_relief_scale_m", base_level_relief_scale_m)?;
            check("base_level_taper_exponent", base_level_taper_exponent)?;
            check("alluvial_grade_slope", alluvial_grade_slope)?;
            check("shelf_wave_base_depth_m", shelf_wave_base_depth_m)?;
            if !l.land_deposition.is_legacy() || !l.sediment_placement.is_legacy() {
                return Err(UniverseError::Unsupported {
                    law: "sediment-retention/continental-retention-base-level-v1 with erosion-deposition-v1 or shelf-wave-base-v1",
                    why: "the retention step replaces the legacy erosion step and placement; a combination with lane SD's laws is not implemented",
                });
            }
        }
        if let SedimentPlacementLawV1::ShelfWaveBaseV1 {
            shelf_wave_base_depth_m,
        } = &l.sediment_placement
        {
            check("shelf_wave_base_depth_m", shelf_wave_base_depth_m)?;
        }
        Ok(())
    }

    /// The sea-level input the active law reads (validated universes only).
    pub fn sea_level_input(&self) -> Result<SeaLevelInputV1, UniverseError> {
        let k = &self.constants;
        match self.law_set.sea_level {
            SeaLevelLawV1::PrescribedDatumV0 => Ok(SeaLevelInputV1::Prescribed {
                sea_level_m: k.sea_level_datum_m.value,
            }),
            SeaLevelLawV1::VolumePrimaryGlobalFillV0 => Ok(SeaLevelInputV1::OceanVolume {
                volume_m3: k
                    .ocean_water_volume_m3
                    .as_ref()
                    .ok_or(UniverseError::Inconsistent(
                        "the volume-primary sea-level law requires ocean_water_volume_m3",
                    ))?
                    .value,
                output_datum_m: k.sea_level_datum_m.value,
            }),
        }
    }

    /// Sediment loads at the continental crust density under every v1
    /// isostasy law (A-SEDIMENT-DENSITY; a law rule, not a stored input).
    pub fn sediment_density_kg_m3(&self) -> f64 {
        self.constants.continental_crust_density_kg_m3.value
    }

    /// blake3 over the canonical JSON (source labels included). For
    /// recipes; the published params digest covers a pinned universe as
    /// part of the profile.
    pub fn digest(&self) -> Result<String, crate::ContractError> {
        crate::artifact::digest_canonical(self)
    }

    /// Compare the physics of two universes: law ids, law parameters,
    /// constant values and units. Source labels are ignored. Returns the
    /// first differing JSON path.
    pub fn same_physics(&self, other: &UniverseV1) -> Result<(), UniverseError> {
        let a = physics_value(self)?;
        let b = physics_value(other)?;
        match first_difference(&a, &b, String::new()) {
            None => Ok(()),
            Some((path, x, y)) => Err(UniverseError::PinMismatch {
                path,
                pinned: x,
                resolved: y,
            }),
        }
    }
}

fn physics_value(u: &UniverseV1) -> Result<serde_json::Value, UniverseError> {
    let mut v = serde_json::to_value(u).map_err(|e| UniverseError::Serialization(e.to_string()))?;
    strip_sources(&mut v);
    Ok(v)
}

fn strip_sources(v: &mut serde_json::Value) {
    match v {
        serde_json::Value::Object(m) => {
            m.remove("source");
            for x in m.values_mut() {
                strip_sources(x);
            }
        }
        serde_json::Value::Array(a) => {
            for x in a {
                strip_sources(x);
            }
        }
        _ => {}
    }
}

fn first_difference(
    a: &serde_json::Value,
    b: &serde_json::Value,
    path: String,
) -> Option<(String, String, String)> {
    use serde_json::Value::Object;
    match (a, b) {
        (Object(x), Object(y)) => {
            let mut keys: Vec<&String> = x.keys().chain(y.keys()).collect();
            keys.sort();
            keys.dedup();
            for k in keys {
                let p = format!("{path}/{k}");
                match (x.get(k), y.get(k)) {
                    (Some(xa), Some(yb)) => {
                        if let Some(d) = first_difference(xa, yb, p) {
                            return Some(d);
                        }
                    }
                    (xa, yb) => {
                        return Some((
                            p,
                            xa.map_or("absent".into(), |v| v.to_string()),
                            yb.map_or("absent".into(), |v| v.to_string()),
                        ))
                    }
                }
            }
            None
        }
        _ if a == b => None,
        _ => Some((
            if path.is_empty() { "/".into() } else { path },
            a.to_string(),
            b.to_string(),
        )),
    }
}

// ===========================================================================
// CLIMATE-V0 (lane CL, 2026-10-12): the coarse climate law and its constants.
//
// ADDITIVE ONLY. Nothing above this line changed except the new `UnitV1`
// variants, which no existing constant uses: every `LawSetV1`,
// `UniverseConstantsV1` and pinned profile serializes byte-identically.
//
// Why a separate record and not new `LawSetV1` / `UniverseConstantsV1`
// fields: those structs are built by struct literal in params.rs
// (`resolve_universe`, owned elsewhere), and the climate tier is not part of
// the genesis recipe that the preview generator replays. `ClimateRecipeV1`
// is what a climate artifact pins (its digest is in the artifact). Folding
// the planet/star rows (luminosity, orbital distance, tilt, rotation) into
// `UniverseConstantsV1` as optional fields is a later, separate change
// (ADR 039 §2 lists them; WS-CONST deferred them for lack of a consumer).
//
// Admissibility (ADR 039 C1, ADR 026 §1): every approximation the law
// implements is named by the implementation (pf-surface `climate.rs`,
// `CLIMATE_APPROXIMATIONS_V0`) and recorded in the artifact; every constant
// carries value, unit and source; the envelope below is CODE-owned and a
// recipe cannot widen it. Every value here was recalled or chosen, never
// checked against a source in this session, so every source is
// `ProjectAssumptionAwaitingReview` (ADR 039 C5).
// ===========================================================================

/// Schema id of a climate recipe (law + constants).
pub const CLIMATE_RECIPE_SCHEMA_V1: &str = "pf.climate-recipe/v1";
/// Schema id of the climate constants record.
pub const CLIMATE_CONSTANTS_SCHEMA_V1: &str = "pf.climate-constants/v1";

/// The climate law a climate artifact was produced under. An unknown id is
/// a load error (serde enum), exactly like the `LawSetV1` entries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ClimateLawV1 {
    /// Declared coarse tier (roadmap CLIMATE-V0): zonal Budyko-Sellers
    /// energy balance with ice albedo, lapse rate, three-cell zonal wind
    /// bands, zonal moisture transport with orographic (saturation)
    /// rainout, Thornthwaite PET. Implemented by
    /// `pf_surface::climate::derive_climate_coarse_v0`.
    ClimateCoarseV0,
    /// SEED-CLIMATE tier (lane CM, 2026-10-12): 2-D moisture transport
    /// with meridional band wind and eddy diffusion, saturation excess
    /// removed over a declared fallout length, Budyko land
    /// evapotranspiration. Implemented by
    /// `pf_surface::climate::seed::derive_seed_climate_v1`.
    SeedClimateV1,
    /// ROTATION-CLIMATE tier (lane RC, 2026-10-20): `seed-climate-v1`'s
    /// moisture solver with circulation geometry that follows the planet:
    /// Hadley edge from an Earth-calibrated Held-style scaling in rotation
    /// rate, radius and gravity, eddy-driven band count from a Rhines-scale
    /// rule, a global Hadley cell for slow rotators. Reads
    /// `rotation_period_s` (with its own envelope,
    /// `CLIMATE_V2_ROTATION_SPEC`). Implemented by
    /// `pf_surface::climate::rotation::derive_seed_climate_v2`.
    SeedClimateV2,
}

impl ClimateLawV1 {
    pub fn id(self) -> &'static str {
        match self {
            ClimateLawV1::ClimateCoarseV0 => "climate-coarse-v0",
            ClimateLawV1::SeedClimateV1 => "seed-climate-v1",
            ClimateLawV1::SeedClimateV2 => "seed-climate-v2",
        }
    }
}

/// Code-owned envelopes of the climate constants (inclusive bounds). The
/// bounds are project assumptions about where the named approximations
/// stay meaningful; the `why` names the approximation that limits them.
pub const CLIMATE_CONSTANT_SPECS_V1: &[ConstantSpecV1] = &[
    ConstantSpecV1 { id: "star_luminosity_w", unit: UnitV1::Watt, lo: 1.0e25, hi: 1.0e28, why: "stellar luminosity (the energy-balance envelope is enforced on the derived solar constant)" },
    ConstantSpecV1 { id: "orbital_distance_m", unit: UnitV1::Meter, lo: 1.0e10, hi: 1.0e13, why: "circular-orbit distance" },
    ConstantSpecV1 { id: "axial_tilt_deg", unit: UnitV1::Degree, lo: 0.0, hi: 40.0, why: "annual-mean zonal EBM with ice albedo and three-cell bands: Earth-like obliquity only (high obliquity reverses the pole-equator insolation contrast)" },
    ConstantSpecV1 { id: "rotation_period_s", unit: UnitV1::Second, lo: 43_200.0, hi: 172_800.0, why: "three-cell (Hadley/Ferrel/polar) wind bands assume an Earth-like rotation rate (0.5 to 2 days)" },
    ConstantSpecV1 { id: "olr_intercept_w_m2", unit: UnitV1::WattPerSquareMeter, lo: 150.0, hi: 260.0, why: "linear OLR = A + B T (T in degC); A is the greenhouse proxy" },
    ConstantSpecV1 { id: "olr_slope_w_m2_k", unit: UnitV1::WattPerSquareMeterKelvin, lo: 1.0, hi: 3.0, why: "linear OLR slope" },
    ConstantSpecV1 { id: "meridional_transport_w_m2_k", unit: UnitV1::WattPerSquareMeterKelvin, lo: 0.5, hi: 10.0, why: "Budyko relaxation heat transport coefficient" },
    ConstantSpecV1 { id: "albedo_ice_free", unit: UnitV1::Dimensionless, lo: 0.05, hi: 0.6, why: "planetary albedo without ice" },
    ConstantSpecV1 { id: "albedo_ice", unit: UnitV1::Dimensionless, lo: 0.3, hi: 0.9, why: "planetary albedo of ice-covered latitudes" },
    ConstantSpecV1 { id: "ice_albedo_threshold_c", unit: UnitV1::DegreeCelsius, lo: -30.0, hi: 0.0, why: "zonal sea-level temperature below which a latitude is ice-covered" },
    ConstantSpecV1 { id: "lapse_rate_k_per_m", unit: UnitV1::KelvinPerMeter, lo: 0.003, hi: 0.0098, why: "mean environmental lapse rate, at most the dry adiabatic 9.8 K/km" },
    ConstantSpecV1 { id: "trade_wind_speed_m_s", unit: UnitV1::MeterPerSecond, lo: 0.1, hi: 30.0, why: "band amplitude; positive so every row has a defined wind direction" },
    ConstantSpecV1 { id: "westerly_wind_speed_m_s", unit: UnitV1::MeterPerSecond, lo: 0.1, hi: 30.0, why: "band amplitude; positive so every row has a defined wind direction" },
    ConstantSpecV1 { id: "polar_easterly_speed_m_s", unit: UnitV1::MeterPerSecond, lo: 0.1, hi: 30.0, why: "band amplitude; positive so every row has a defined wind direction" },
    ConstantSpecV1 { id: "column_water_saturation_ref_kg_m2", unit: UnitV1::KilogramPerSquareMeter, lo: 1.0, hi: 100.0, why: "saturated column water at the reference temperature" },
    ConstantSpecV1 { id: "column_water_reference_temperature_c", unit: UnitV1::DegreeCelsius, lo: -20.0, hi: 40.0, why: "reference temperature of the saturated column" },
    ConstantSpecV1 { id: "clausius_clapeyron_rate_per_k", unit: UnitV1::PerKelvin, lo: 0.05, hi: 0.08, why: "exponential Clausius-Clapeyron scaling of column water near Earth surface temperatures" },
    ConstantSpecV1 { id: "ocean_column_relative_humidity", unit: UnitV1::Dimensionless, lo: 0.3, hi: 1.0, why: "ocean equilibrium column relative to saturation" },
    ConstantSpecV1 { id: "ocean_evaporation_relaxation_s", unit: UnitV1::Second, lo: 3_600.0, hi: 3.0e7, why: "relaxation time of the column toward the ocean equilibrium" },
    ConstantSpecV1 { id: "rainout_residence_s", unit: UnitV1::Second, lo: 86_400.0, hi: 3.0e7, why: "background rainout residence time of column water" },
    ConstantSpecV1 { id: "moisture_transport_speed_m_s", unit: UnitV1::MeterPerSecond, lo: 0.5, hi: 30.0, why: "effective column moisture transport speed along the prevailing wind" },
    ConstantSpecV1 { id: "land_recycling_ratio", unit: UnitV1::Dimensionless, lo: 0.0, hi: 0.9, why: "fraction of background land rainout re-evaporated into the column; below 1 so land columns lose water" },
    ConstantSpecV1 { id: "circulation_rainout_modulation", unit: UnitV1::Dimensionless, lo: 0.0, hi: 0.9, why: "amplitude of the three-cell rising/sinking rainout factor 1 + k cos(6 lat); below 1 keeps rainout positive" },
    ConstantSpecV1 { id: "ocean_seasonal_damping", unit: UnitV1::Dimensionless, lo: 0.0, hi: 1.0, why: "ocean seasonal amplitude as a fraction of the land equilibrium response" },
    ConstantSpecV1 { id: "continentality_length_m", unit: UnitV1::Meter, lo: 1.0e4, hi: 1.0e7, why: "e-folding upwind distance from the coast for the seasonal amplitude" },
];

/// ROTATION-CLIMATE (lane RC, 2026-10-20): the rotation envelope of
/// `seed-climate-v2`, replacing the three-cell 0.5-2 day envelope for that
/// law only (the other laws keep `CLIMATE_CONSTANT_SPECS_V1`). Project
/// assumption awaiting review: below 4 h the Rhines band rule asks for more
/// eddy-driven bands than a coarse lat/lon grid resolves and the planet nears
/// rotational break-up for Earth-like densities; above 100 days the
/// daily-mean insolation the energy balance uses stops being a fair average
/// (the solar day approaches the year and a day-night circulation, which no
/// tier models, takes over). Between, the Hadley edge saturates to a global
/// cell, so slow rotators are inside the law, not refused.
pub const CLIMATE_V2_ROTATION_SPEC: ConstantSpecV1 = ConstantSpecV1 {
    id: "rotation_period_s",
    unit: UnitV1::Second,
    lo: 14_400.0,
    hi: 8_640_000.0,
    why: "seed-climate-v2 rotation-scaled bands: 4 h (band count vs grid, break-up) to 100 days (daily-mean insolation, no day-night circulation)",
};

/// Derived-quantity envelope: the solar constant S0 = L / (4 pi d^2).
pub const CLIMATE_SOLAR_CONSTANT_LO_W_M2: f64 = 1000.0;
pub const CLIMATE_SOLAR_CONSTANT_HI_W_M2: f64 = 1800.0;

/// Look up a climate constant spec (panics on an unknown id: a programming
/// error, never data).
pub fn climate_constant_spec(id: &str) -> &'static ConstantSpecV1 {
    CLIMATE_CONSTANT_SPECS_V1
        .iter()
        .find(|s| s.id == id)
        .unwrap_or_else(|| panic!("no climate v1 constant spec for {id}"))
}

fn check_climate(id: &'static str, c: &ConstantV1) -> Result<(), UniverseError> {
    check_climate_spec(id, climate_constant_spec(id), c)
}

fn check_climate_spec(id: &'static str, spec: &ConstantSpecV1, c: &ConstantV1) -> Result<(), UniverseError> {
    if c.unit != spec.unit {
        return Err(UniverseError::UnitMismatch {
            constant: id,
            expected: spec.unit.as_str(),
            found: c.unit.as_str(),
        });
    }
    if !c.value.is_finite() || !(spec.lo..=spec.hi).contains(&c.value) {
        return Err(UniverseError::OutOfEnvelope {
            constant: id,
            value: c.value,
            unit: spec.unit.as_str(),
            lo: spec.lo,
            hi: spec.hi,
            why: spec.why,
        });
    }
    Ok(())
}

/// The constants the coarse climate tier reads. Planet radius is NOT here:
/// it comes from the grid the upstream terrain was published on.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClimateConstantsV1 {
    pub schema: Cow<'static, str>,
    // Star and orbit (ADR 039 §2 "Star" rows). S0 is DERIVED from these.
    pub star_luminosity_w: ConstantV1,
    pub orbital_distance_m: ConstantV1,
    // Planet (ADR 039 §2 "Planet" rows).
    pub axial_tilt_deg: ConstantV1,
    /// Read only by the validity envelope (three-cell band structure).
    pub rotation_period_s: ConstantV1,
    // Zonal energy balance (approximation `budyko-sellers-zonal-ebm-v0`).
    pub olr_intercept_w_m2: ConstantV1,
    pub olr_slope_w_m2_k: ConstantV1,
    pub meridional_transport_w_m2_k: ConstantV1,
    pub albedo_ice_free: ConstantV1,
    pub albedo_ice: ConstantV1,
    pub ice_albedo_threshold_c: ConstantV1,
    // Lapse rate (approximation `constant-environmental-lapse-rate-v0`).
    pub lapse_rate_k_per_m: ConstantV1,
    // Wind bands (approximation `three-cell-zonal-wind-bands-v0`).
    pub trade_wind_speed_m_s: ConstantV1,
    pub westerly_wind_speed_m_s: ConstantV1,
    pub polar_easterly_speed_m_s: ConstantV1,
    // Moisture (approximations `zonal-moisture-transport-v0`,
    // `saturation-rainout-orographic-v0`, `three-cell-rainout-modulation-v0`).
    pub column_water_saturation_ref_kg_m2: ConstantV1,
    pub column_water_reference_temperature_c: ConstantV1,
    pub clausius_clapeyron_rate_per_k: ConstantV1,
    pub ocean_column_relative_humidity: ConstantV1,
    pub ocean_evaporation_relaxation_s: ConstantV1,
    pub rainout_residence_s: ConstantV1,
    pub moisture_transport_speed_m_s: ConstantV1,
    pub land_recycling_ratio: ConstantV1,
    pub circulation_rainout_modulation: ConstantV1,
    // Seasonality (approximation `continentality-damped-seasonal-response-v0`).
    pub ocean_seasonal_damping: ConstantV1,
    pub continentality_length_m: ConstantV1,
}

fn awaiting(basis: &'static str) -> SourceClassV1 {
    SourceClassV1::ProjectAssumptionAwaitingReview {
        basis: Cow::Borrowed(basis),
    }
}

impl ClimateConstantsV1 {
    /// Earth-like v0 values. EVERY value is recalled or chosen and was not
    /// checked against its source in this session (ADR 039 C5).
    pub fn earth_v0() -> Self {
        let k = |value, unit, basis| ConstantV1::new(value, unit, awaiting(basis));
        Self {
            schema: Cow::Borrowed(CLIMATE_CONSTANTS_SCHEMA_V1),
            star_luminosity_w: k(3.828e26, UnitV1::Watt, "recalled: IAU 2015 Resolution B3 nominal solar luminosity; not re-checked"),
            orbital_distance_m: k(1.495_978_707e11, UnitV1::Meter, "recalled: IAU 2012 Resolution B2 astronomical unit; circular orbit assumed; not re-checked"),
            axial_tilt_deg: k(23.44, UnitV1::Degree, "recalled: present Earth obliquity; not re-checked"),
            rotation_period_s: k(86_164.1, UnitV1::Second, "recalled: Earth sidereal day; read only by the wind-band validity envelope; not re-checked"),
            olr_intercept_w_m2: k(202.0, UnitV1::WattPerSquareMeter, "recalled: Budyko (1969) linear OLR intercept as tabulated in textbook Budyko-Sellers treatments (T in degC); exact source and value not re-checked"),
            olr_slope_w_m2_k: k(1.90, UnitV1::WattPerSquareMeterKelvin, "recalled: Budyko (1969) linear OLR slope as tabulated in textbook treatments; not re-checked"),
            meridional_transport_w_m2_k: k(3.04, UnitV1::WattPerSquareMeterKelvin, "recalled: Budyko relaxation transport coefficient as tabulated in textbook treatments; not re-checked"),
            albedo_ice_free: k(0.32, UnitV1::Dimensionless, "recalled: Budyko-Sellers ice-free planetary albedo; not re-checked"),
            albedo_ice: k(0.62, UnitV1::Dimensionless, "recalled: Budyko-Sellers ice albedo; not re-checked"),
            ice_albedo_threshold_c: k(-10.0, UnitV1::DegreeCelsius, "recalled: Budyko ice-line temperature; not re-checked"),
            lapse_rate_k_per_m: k(0.0065, UnitV1::KelvinPerMeter, "recalled: ICAO / ISO 2533 standard-atmosphere tropospheric lapse rate, used as a mean environmental lapse rate; not re-checked"),
            trade_wind_speed_m_s: k(6.0, UnitV1::MeterPerSecond, "project assumption: order of Earth's zonal-mean surface trade winds (recalled, not checked against a reanalysis)"),
            westerly_wind_speed_m_s: k(7.0, UnitV1::MeterPerSecond, "project assumption: order of Earth's zonal-mean surface westerlies (recalled, not checked against a reanalysis)"),
            polar_easterly_speed_m_s: k(3.0, UnitV1::MeterPerSecond, "project assumption: order of Earth's polar easterlies (recalled, not checked against a reanalysis)"),
            column_water_saturation_ref_kg_m2: k(30.0, UnitV1::KilogramPerSquareMeter, "project assumption: chosen so the ocean equilibrium column (x 0.8) is ~24 kg/m2, near the recalled global-mean precipitable water of ~25 kg/m2; not re-checked"),
            column_water_reference_temperature_c: k(15.0, UnitV1::DegreeCelsius, "project assumption: recalled Earth global-mean surface temperature"),
            clausius_clapeyron_rate_per_k: k(0.07, UnitV1::PerKelvin, "recalled: ~7 %/K Clausius-Clapeyron scaling of near-surface water vapour (e.g. Held and Soden 2006); not re-checked"),
            ocean_column_relative_humidity: k(0.8, UnitV1::Dimensionless, "recalled: marine boundary-layer relative humidity of order 0.8; not re-checked"),
            ocean_evaporation_relaxation_s: k(172_800.0, UnitV1::Second, "project assumption: 2-day relaxation of the column toward the ocean equilibrium"),
            rainout_residence_s: k(734_400.0, UnitV1::Second, "recalled: atmospheric water-vapour residence time of ~8.5 days (Trenberth 1998; van der Ent and Tuinenburg 2017); not re-checked"),
            moisture_transport_speed_m_s: k(5.0, UnitV1::MeterPerSecond, "project assumption: effective column moisture transport speed"),
            land_recycling_ratio: k(0.3, UnitV1::Dimensionless, "recalled: continental moisture recycling ratios of order 0.3-0.4 (Eltahir and Bras 1996; van der Ent et al. 2010); not re-checked"),
            circulation_rainout_modulation: k(0.6, UnitV1::Dimensionless, "project assumption: amplitude of the rising (0, 60 deg) / sinking (30, 90 deg) rainout factor"),
            ocean_seasonal_damping: k(0.2, UnitV1::Dimensionless, "project assumption: ocean seasonal amplitude ~ one fifth of the land equilibrium response (recalled order, not checked)"),
            continentality_length_m: k(1.0e6, UnitV1::Meter, "project assumption: e-folding upwind distance from the coast for the seasonal amplitude"),
        }
    }

    /// Derived solar constant S0 = L / (4 pi d^2) in W/m^2 (ADR 039 §2:
    /// derived quantities are never stored as inputs).
    pub fn solar_constant_w_m2(&self) -> f64 {
        let d = self.orbital_distance_m.value;
        self.star_luminosity_w.value / (4.0 * std::f64::consts::PI * d * d)
    }
}

/// What a climate artifact pins: the law and its constants.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClimateRecipeV1 {
    pub schema: Cow<'static, str>,
    pub law: ClimateLawV1,
    pub constants: ClimateConstantsV1,
}

impl ClimateRecipeV1 {
    /// The declared coarse tier with Earth-like v0 constants.
    pub fn coarse_v0_earth() -> Self {
        Self {
            schema: Cow::Borrowed(CLIMATE_RECIPE_SCHEMA_V1),
            law: ClimateLawV1::ClimateCoarseV0,
            constants: ClimateConstantsV1::earth_v0(),
        }
    }

    /// Fail-loud validation with typed refusals (schemas, units, code-owned
    /// envelopes, the derived solar constant, ordering, law support).
    pub fn validate(&self) -> Result<(), UniverseError> {
        if self.schema != CLIMATE_RECIPE_SCHEMA_V1 {
            return Err(UniverseError::Schema {
                what: "climate recipe",
                found: self.schema.to_string(),
                expected: CLIMATE_RECIPE_SCHEMA_V1,
            });
        }
        let k = &self.constants;
        if k.schema != CLIMATE_CONSTANTS_SCHEMA_V1 {
            return Err(UniverseError::Schema {
                what: "climate constants",
                found: k.schema.to_string(),
                expected: CLIMATE_CONSTANTS_SCHEMA_V1,
            });
        }
        // SEED-CLIMATE (lane CM, 2026-10-12): seed-climate-v1 is implemented
        // by pf_surface::climate::seed::derive_seed_climate_v1 (same base
        // constants and envelopes; its extra constants carry their own
        // code-owned envelope in pf-surface). climate-coarse-v0's
        // derive_climate_coarse_v0 still refuses it as Unsupported.
        for (id, c) in [
            ("star_luminosity_w", &k.star_luminosity_w),
            ("orbital_distance_m", &k.orbital_distance_m),
            ("axial_tilt_deg", &k.axial_tilt_deg),
            ("rotation_period_s", &k.rotation_period_s),
            ("olr_intercept_w_m2", &k.olr_intercept_w_m2),
            ("olr_slope_w_m2_k", &k.olr_slope_w_m2_k),
            ("meridional_transport_w_m2_k", &k.meridional_transport_w_m2_k),
            ("albedo_ice_free", &k.albedo_ice_free),
            ("albedo_ice", &k.albedo_ice),
            ("ice_albedo_threshold_c", &k.ice_albedo_threshold_c),
            ("lapse_rate_k_per_m", &k.lapse_rate_k_per_m),
            ("trade_wind_speed_m_s", &k.trade_wind_speed_m_s),
            ("westerly_wind_speed_m_s", &k.westerly_wind_speed_m_s),
            ("polar_easterly_speed_m_s", &k.polar_easterly_speed_m_s),
            ("column_water_saturation_ref_kg_m2", &k.column_water_saturation_ref_kg_m2),
            ("column_water_reference_temperature_c", &k.column_water_reference_temperature_c),
            ("clausius_clapeyron_rate_per_k", &k.clausius_clapeyron_rate_per_k),
            ("ocean_column_relative_humidity", &k.ocean_column_relative_humidity),
            ("ocean_evaporation_relaxation_s", &k.ocean_evaporation_relaxation_s),
            ("rainout_residence_s", &k.rainout_residence_s),
            ("moisture_transport_speed_m_s", &k.moisture_transport_speed_m_s),
            ("land_recycling_ratio", &k.land_recycling_ratio),
            ("circulation_rainout_modulation", &k.circulation_rainout_modulation),
            ("ocean_seasonal_damping", &k.ocean_seasonal_damping),
            ("continentality_length_m", &k.continentality_length_m),
        ] {
            // ROTATION-CLIMATE (lane RC): seed-climate-v2 reads the rotation
            // period through its own law, so it carries its own envelope.
            if id == "rotation_period_s" && self.law == ClimateLawV1::SeedClimateV2 {
                check_climate_spec(id, &CLIMATE_V2_ROTATION_SPEC, c)?;
            } else {
                check_climate(id, c)?;
            }
        }
        let s0 = k.solar_constant_w_m2();
        if !s0.is_finite()
            || !(CLIMATE_SOLAR_CONSTANT_LO_W_M2..=CLIMATE_SOLAR_CONSTANT_HI_W_M2).contains(&s0)
        {
            return Err(UniverseError::OutOfEnvelope {
                constant: "solar_constant_w_m2 (derived: star_luminosity_w / (4 pi orbital_distance_m^2))",
                value: s0,
                unit: "W/m2",
                lo: CLIMATE_SOLAR_CONSTANT_LO_W_M2,
                hi: CLIMATE_SOLAR_CONSTANT_HI_W_M2,
                why: "the linear-OLR energy balance (A, B) is calibrated near present Earth; far from it the outgoing longwave is not linear in temperature",
            });
        }
        if k.albedo_ice.value <= k.albedo_ice_free.value {
            return Err(UniverseError::Inconsistent(
                "albedo_ice must exceed albedo_ice_free (ice-albedo feedback sign)",
            ));
        }
        Ok(())
    }

    /// blake3 over the canonical JSON (source labels included).
    pub fn digest(&self) -> Result<String, crate::ContractError> {
        crate::artifact::digest_canonical(self)
    }
}
