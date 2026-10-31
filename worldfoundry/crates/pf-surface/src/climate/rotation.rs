//! ROTATION-CLIMATE (lane RC, 2026-10-20): the `seed-climate-v2` tier.
//!
//! Problem it fixes: `climate-coarse-v0` and `seed-climate-v1` hard-wire
//! Earth's circulation (band edges at 30 and 60 deg, `sin(6 lat)` winds,
//! `cos(6 lat)` rainout), and the recipe refuses any rotation period outside
//! 0.5-2 days. A parametric planet with another rotation rate or radius got
//! Earth's wind belts, ITCZ, subtropical dry belts and rain shadows.
//!
//! WHAT CHANGES AGAINST `seed-climate-v1` (which stays untouched and
//! byte-identical; this tier runs only when the recipe names
//! `seed-climate-v2`): only the circulation GEOMETRY. The moisture solver,
//! energy balance, fallout, Budyko evapotranspiration, continentality, PET
//! and aridity are `seed-climate-v1`'s, called through the same core
//! (`seed::seed_climate_core`) with different band functions.
//!
//! THE LAW (every constant a project assumption awaiting review, ADR 039 C5):
//! 1. Hadley edge (HADLEY-P, integrator decision 2026-10-20 23:35; lane
//!    CLIMATE-GENESIS-2): eddy-limited scaling on the SINE of the edge,
//!    calibrated to Earth: sin(phi_H) = sin(phi_E) * S^(-p) * G^(p/2), with
//!    S = (Omega a) / (Omega_E a_E) and G = (g H) / (g_E H_E); phi_H = 90
//!    where the right side reaches 1. Default p = 1/3: Mitchell and Hill
//!    (2021, J. Atmos. Sci. 78:1445, arXiv 1911.05860; primary text read by
//!    lane EARTH-NUMBERS) find sin(phi_H) ~ Omega^(-1/3) in all three
//!    idealized GCMs over Omega* 1/16-2; lane TW-ORACLE's fit to the
//!    ThousandWorlds GCM set gives p 0.25-0.4 best, 1/2 close, 1 worst
//!    (tidally locked runs; shape only). Envelope 0.25-0.5: Held (2000)'s
//!    two-layer edge (p = 1/2) is the upper edge; the axisymmetric
//!    Held-Hou (1980) p = 1 is an upper bound, not a plausible value, and
//!    is outside the envelope. Circulation depth H scales as g^(-e)
//!    (e = 1: the tropopause depth follows the scale height R T / g, so
//!    gravity cancels in g H; e = 0: fixed depth, gravity enters as
//!    g^(p/2)). The pole-equator contrast dtheta is held at Earth's value.
//! 2. Eddy-driven bands. Poleward of phi_H the extratropics hold m pairs of
//!    (westerly, easterly) bands. A pair is one eddy-driven jet spacing,
//!    which follows the Rhines scale L ~ sqrt(U / beta), beta = 2 Omega
//!    cos(lat) / a, so its angular width scales as S^(-q), q = 1/2, with U
//!    held at Earth's value; Earth's pair (Ferrel 30-60 + polar 60-90) is
//!    60 deg wide. With x = W / P (extratropical width W = 90 - phi_H,
//!    pair width P), the band pattern is a BLEND of the m = max(1,
//!    floor(x)) pair geometry and the m + 1 pair geometry (each splitting W
//!    evenly) with weight w = x - m (0 below x = 1): winds and rainout
//!    modulation change continuously with rotation instead of jumping when
//!    round(x) steps (HADLEY-P "smooth the band-count change"). At an
//!    integer x (Earth: x = 1) it is the single geometry exactly. Faster
//!    rotation gives more, narrower bands (the multiple-jet regime of
//!    fast-rotating GCMs, e.g. Kaspi and Showman 2015, recalled).
//! 3. Global cell. If W < f * P the Hadley cell covers the whole hemisphere
//!    (slow-rotator regime, Venus/Titan-like): trades and equatorward flow to
//!    the pole. Default f = 0.13 puts the switch at S = 0.15 with p = 1/3
//!    (TW-ORACLE: bounded cells in 4-6 % of runs at S < 0.1, 33 % at
//!    0.1-0.2, 96-100 % above 0.2).
//! 4. Band profiles generalize v1's sin(6 lat) / cos(6 lat) per band k with
//!    edges [lo, hi]: x = pi (|lat| - lo) / (hi - lo); thermally direct
//!    bands (Hadley, extratropical easterly) have u = -U sin x, v =
//!    -sign(lat) V sin x (equatorward), rainout 1 + k_c cos x (rising at the
//!    equatorward edge); eddy-driven westerly bands flip all three signs.
//!    On Earth (edges 0, 30, 60, 90) this is v1's geometry exactly, up to
//!    floating-point rounding of the sine arguments.
//!
//! DETERMINISM: as `seed-climate-v1` (fixed sweep order, no RNG, no clock).

use super::seed::{seed_climate_core_with, seed_climate_v1_earth_recipe, SeedClimateConstantsV1, SeedColumnStateV1, SeedSolveOptions, SEED_CLIMATE_APPROXIMATIONS_V1};
use super::{ApproximationV0, ClimateError, ClimateFieldSpecV0, ClimateInputs, ClimateOutputV0, FieldStatusV0, ARIDITY_SENTINEL_PET_ZERO, DISTANCE_SENTINEL_NO_UPWIND_OCEAN};
use pf_contracts::laws::{ClimateLawV1, ClimateRecipeV1, ConstantSpecV1, ConstantV1, SourceClassV1, UnitV1, UniverseError};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::f64::consts::PI;

/// Tier id (= the law id it implements) and version.
pub const SEED_CLIMATE_V2_TIER_ID: &str = "seed-climate-v2";
/// 2: HADLEY-P (lane CLIMATE-GENESIS-2, 2026-10-20): sine-form edge with
/// p = 1/3, global cell at S < 0.15, blended band count. Version 1 (lane RC:
/// phi_H = phi_E S^-p, p = 1/2, f = 1/2, rounded band count) is not kept.
pub const SEED_CLIMATE_V2_TIER_VERSION: u32 = 2;
/// Schema of the rotation-scaling constants record.
pub const ROTATION_CLIMATE_CONSTANTS_SCHEMA_V2: &str = "pf.rotation-climate-constants/v2";

/// The constants `seed-climate-v2` reads beyond `ClimateConstantsV1` and
/// `SeedClimateConstantsV1`. Code-owned envelope below; recorded with its
/// digest in the artifact.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RotationClimateConstantsV2 {
    pub schema: Cow<'static, str>,
    /// Calibration planet: rotation period, radius and gravity at which the
    /// reference edge and pair width hold (Earth).
    pub reference_rotation_period_s: ConstantV1,
    pub reference_planet_radius_m: ConstantV1,
    pub reference_surface_gravity_m_s2: ConstantV1,
    /// The planet's surface gravity. The climate inputs carry no universe
    /// constants, so this is a recipe value; it should equal the world's
    /// `UniverseConstantsV1::surface_gravity_m_s2` (9.81 m/s2 for the
    /// current profiles).
    pub surface_gravity_m_s2: ConstantV1,
    /// Hadley edge latitude of the calibration planet.
    pub reference_hadley_edge_deg: ConstantV1,
    /// p in sin(phi_H) ~ S^(-p) G^(p/2) (1/3 = Mitchell and Hill 2021,
    /// 1/2 = Held 2000).
    pub hadley_scaling_exponent: ConstantV1,
    /// e in H ~ g^(-e) (1: depth follows the scale height; 0: fixed depth).
    pub circulation_depth_gravity_exponent: ConstantV1,
    /// Angular width of one (westerly, easterly) extratropical band pair on
    /// the calibration planet.
    pub reference_eddy_pair_width_deg: ConstantV1,
    /// q in pair width ~ S^(-q) (Rhines scale: 1/2).
    pub eddy_pair_scaling_exponent: ConstantV1,
    /// f: the Hadley cell becomes global when the extratropics are narrower
    /// than f pair widths.
    pub global_cell_pair_fraction: ConstantV1,
}

/// Code-owned envelopes (inclusive).
pub const ROTATION_CLIMATE_CONSTANT_SPECS_V2: &[ConstantSpecV1] = &[
    ConstantSpecV1 { id: "reference_rotation_period_s", unit: UnitV1::Second, lo: 43_200.0, hi: 172_800.0, why: "calibration planet must be Earth-like (the reference edge and pair width are Earth's)" },
    ConstantSpecV1 { id: "reference_planet_radius_m", unit: UnitV1::Meter, lo: 3.0e6, hi: 1.2e7, why: "calibration planet must be Earth-like" },
    ConstantSpecV1 { id: "reference_surface_gravity_m_s2", unit: UnitV1::MeterPerSecondSquared, lo: 5.0, hi: 20.0, why: "calibration planet must be Earth-like" },
    ConstantSpecV1 { id: "surface_gravity_m_s2", unit: UnitV1::MeterPerSecondSquared, lo: 0.5, hi: 50.0, why: "terrestrial planets with a thin troposphere (Titan 1.35 to super-Earths); the scaling is a power law in g H" },
    ConstantSpecV1 { id: "reference_hadley_edge_deg", unit: UnitV1::Degree, lo: 20.0, hi: 40.0, why: "Earth's annual-mean Hadley terminus / subtropical-high latitude" },
    ConstantSpecV1 { id: "hadley_scaling_exponent", unit: UnitV1::Dimensionless, lo: 0.25, hi: 0.5, why: "eddy-limited GCM scalings: lane TW-ORACLE best fit 0.25-0.4 (ThousandWorlds, tidally locked), Mitchell and Hill 2021 1/3, Held (2000) 1/2; the axisymmetric Held-Hou p = 1 is an upper bound and fits worst (HADLEY-P, 2026-10-20)" },
    ConstantSpecV1 { id: "circulation_depth_gravity_exponent", unit: UnitV1::Dimensionless, lo: 0.0, hi: 1.0, why: "circulation depth between fixed (0) and scale-height-like R T / g (1)" },
    ConstantSpecV1 { id: "reference_eddy_pair_width_deg", unit: UnitV1::Degree, lo: 30.0, hi: 90.0, why: "Earth's extratropics hold one Ferrel + polar pair over about 60 deg" },
    ConstantSpecV1 { id: "eddy_pair_scaling_exponent", unit: UnitV1::Dimensionless, lo: 0.25, hi: 1.0, why: "Rhines-scale jet spacing gives 1/2" },
    ConstantSpecV1 { id: "global_cell_pair_fraction", unit: UnitV1::Dimensionless, lo: 0.0, hi: 1.0, why: "fraction of a jet spacing the extratropics must hold before an eddy-driven band forms" },
];

fn awaiting(basis: &'static str) -> SourceClassV1 {
    SourceClassV1::ProjectAssumptionAwaitingReview { basis: Cow::Borrowed(basis) }
}

impl RotationClimateConstantsV2 {
    /// Earth-calibrated v2 values (every one recalled or chosen; ADR 039 C5).
    pub fn earth_v2() -> Self {
        let k = |value, unit, basis| ConstantV1::new(value, unit, awaiting(basis));
        Self {
            schema: Cow::Borrowed(ROTATION_CLIMATE_CONSTANTS_SCHEMA_V2),
            reference_rotation_period_s: k(86_164.1, UnitV1::Second, "recalled: Earth sidereal day (same value as the climate recipe's Earth rotation_period_s); not re-checked"),
            reference_planet_radius_m: k(6_371_000.0, UnitV1::Meter, "recalled: Earth mean radius (same value as the world profiles' planet_radius_m); not re-checked"),
            reference_surface_gravity_m_s2: k(9.81, UnitV1::MeterPerSecondSquared, "recalled: Earth surface gravity (same stated value as the world profiles); not re-checked"),
            surface_gravity_m_s2: k(9.81, UnitV1::MeterPerSecondSquared, "project assumption: the planet's surface gravity, set equal to the world's stated UniverseConstantsV1 surface_gravity_m_s2 (9.81 in the current profiles); the climate inputs do not carry it"),
            reference_hadley_edge_deg: k(30.0, UnitV1::Degree, "recalled: Earth's annual-mean Hadley cell terminus and subtropical dry belt near 30 deg (the v0/v1 band edge); not re-checked against a reanalysis"),
            hadley_scaling_exponent: k(1.0 / 3.0, UnitV1::Dimensionless, "project assumption awaiting review (HADLEY-P, integrator decision 2026-10-20 23:35): sin(phi_H) ~ Omega^-1/3 in all three idealized GCMs of Mitchell and Hill (2021, J. Atmos. Sci. 78:1445; primary text read by lane EARTH-NUMBERS); lane TW-ORACLE ThousandWorlds fit p = 0.265, 95 % CI 0.21-1.3, rms best at 0.25-0.4 (tidally locked runs, shape only)"),
            circulation_depth_gravity_exponent: k(1.0, UnitV1::Dimensionless, "project assumption: tropopause depth scales as the scale height R T / g (the dry adiabatic lapse rate is g / c_p), so g H, and the Hadley edge, do not depend on gravity at Earth-like temperatures; 0 would hold the depth fixed"),
            reference_eddy_pair_width_deg: k(60.0, UnitV1::Degree, "project assumption: Earth's extratropics (30-90 deg) hold one eddy-driven westerly band and one polar easterly band, the v0/v1 Ferrel and polar cells"),
            eddy_pair_scaling_exponent: k(0.5, UnitV1::Dimensionless, "recalled: Rhines (1975) scale L ~ pi sqrt(2 U / beta), beta = 2 Omega cos(lat) / a, sets eddy-driven jet spacing; with U fixed the angular spacing scales as (Omega a)^-1/2; fast-rotating GCMs show more, narrower jets (e.g. Kaspi and Showman 2015, recalled, not re-checked)"),
            global_cell_pair_fraction: k(0.13, UnitV1::Dimensionless, "project assumption awaiting review (HADLEY-P, 2026-10-20): no eddy-driven band forms unless the extratropics hold 0.13 jet spacings; with p = 1/3 on sin(phi_H) this puts the switch to a global Hadley cell at S = 0.15 (lane TW-ORACLE: bounded cells in 4-6 % of ThousandWorlds runs at S < 0.1, 33 % at 0.1-0.2, 96-100 % above 0.2)"),
        }
    }

    fn entries(&self) -> [(&'static str, &ConstantV1); 10] {
        [
            ("reference_rotation_period_s", &self.reference_rotation_period_s),
            ("reference_planet_radius_m", &self.reference_planet_radius_m),
            ("reference_surface_gravity_m_s2", &self.reference_surface_gravity_m_s2),
            ("surface_gravity_m_s2", &self.surface_gravity_m_s2),
            ("reference_hadley_edge_deg", &self.reference_hadley_edge_deg),
            ("hadley_scaling_exponent", &self.hadley_scaling_exponent),
            ("circulation_depth_gravity_exponent", &self.circulation_depth_gravity_exponent),
            ("reference_eddy_pair_width_deg", &self.reference_eddy_pair_width_deg),
            ("eddy_pair_scaling_exponent", &self.eddy_pair_scaling_exponent),
            ("global_cell_pair_fraction", &self.global_cell_pair_fraction),
        ]
    }

    /// Fail-loud validation (schema, units, code-owned envelopes).
    pub fn validate(&self) -> Result<(), UniverseError> {
        if self.schema != ROTATION_CLIMATE_CONSTANTS_SCHEMA_V2 {
            return Err(UniverseError::Schema {
                what: "rotation climate constants",
                found: self.schema.to_string(),
                expected: ROTATION_CLIMATE_CONSTANTS_SCHEMA_V2,
            });
        }
        for (id, c) in self.entries() {
            let spec = ROTATION_CLIMATE_CONSTANT_SPECS_V2
                .iter()
                .find(|s| s.id == id)
                .unwrap_or_else(|| panic!("no rotation climate spec for {id}"));
            if c.unit != spec.unit {
                return Err(UniverseError::UnitMismatch { constant: spec.id, expected: spec.unit.as_str(), found: c.unit.as_str() });
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

/// The `seed-climate-v2` recipe with Earth constants: `seed-climate-v1`'s
/// recipe under the new law id; only the rotation period's source label
/// changes (v2 reads it through the law, not just an envelope). Select
/// another planet by setting `constants.rotation_period_s` (4 h to 100 d).
pub fn seed_climate_v2_earth_recipe() -> ClimateRecipeV1 {
    let mut r = seed_climate_v1_earth_recipe();
    r.law = ClimateLawV1::SeedClimateV2;
    r.constants.rotation_period_s = ConstantV1::new(
        86_164.1,
        UnitV1::Second,
        awaiting("recalled: Earth sidereal day; read by the seed-climate-v2 circulation law (Hadley edge and eddy band count); not re-checked"),
    );
    r
}

/// Kind of one circulation band (per hemisphere, equator to pole).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BandKindV2 {
    /// Thermally direct tropical cell: easterly trades, equatorward flow.
    Hadley,
    /// Eddy-driven band: surface westerlies, poleward flow (Ferrel-like).
    Westerly,
    /// Thermally direct extratropical band: easterlies, equatorward flow.
    Easterly,
}

impl BandKindV2 {
    /// +1 for thermally direct bands, -1 for eddy-driven westerly bands.
    fn sign(self) -> f64 {
        match self {
            BandKindV2::Hadley | BandKindV2::Easterly => 1.0,
            BandKindV2::Westerly => -1.0,
        }
    }
}

/// The circulation geometry the law derived for one planet (published).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CirculationV2 {
    pub rotation_period_s: f64,
    pub planet_radius_m: f64,
    pub surface_gravity_m_s2: f64,
    /// S = (Omega a) / (Omega_E a_E).
    pub rotation_speed_ratio: f64,
    /// G = (g H) / (g_E H_E) under the depth-gravity exponent.
    pub depth_gravity_ratio: f64,
    /// asin(min(1, sin(phi_E) S^-p G^(p/2))) before the global-cell rule.
    pub hadley_edge_scaled_deg: f64,
    /// The Hadley edge used (90 when the cell is global).
    pub hadley_edge_deg: f64,
    pub global_hadley_cell: bool,
    /// Angular width of one eddy (westerly, easterly) pair.
    pub eddy_pair_width_deg: f64,
    /// m = max(1, floor(W / P)) (0 when global): pairs of the primary
    /// geometry poleward of the Hadley edge.
    pub eddy_pairs: u32,
    /// Band edges from the equator (0) to the pole (90), per hemisphere
    /// (primary geometry, m pairs).
    pub band_edges_deg: Vec<f64>,
    pub band_kinds: Vec<BandKindV2>,
    /// Weight w = W / P - m of the m + 1 pair geometry in the blend (0 =
    /// primary geometry only; skipped on write when 0).
    #[serde(default, skip_serializing_if = "is_zero_f64")]
    pub blend_weight: f64,
    /// The m + 1 pair geometry (present when `blend_weight` > 0).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub blend_band_edges_deg: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub blend_band_kinds: Vec<BandKindV2>,
}

fn is_zero_f64(v: &f64) -> bool {
    *v == 0.0
}

/// Hadley band plus m (westerly, easterly) pairs splitting [phi, 90] evenly.
fn pair_geometry(phi: f64, m: u32) -> (Vec<f64>, Vec<BandKindV2>) {
    let w = 90.0 - phi;
    let nb = 2 * m as usize;
    let mut edges = vec![0.0, phi];
    let mut kinds = vec![BandKindV2::Hadley];
    for b in 0..nb {
        edges.push(if b + 1 == nb { 90.0 } else { phi + (b + 1) as f64 * w / nb as f64 });
        kinds.push(if b % 2 == 0 { BandKindV2::Westerly } else { BandKindV2::Easterly });
    }
    (edges, kinds)
}

/// (kind, phase x in [0, pi]) of |lat| in one geometry.
fn band_of(edges: &[f64], kinds: &[BandKindV2], lat_deg: f64) -> (BandKindV2, f64) {
    let a = lat_deg.abs().min(90.0);
    let nb = kinds.len();
    let mut b = nb - 1;
    for k in 0..nb {
        if a < edges[k + 1] {
            b = k;
            break;
        }
    }
    let (lo, hi) = (edges[b], edges[b + 1]);
    (kinds[b], PI * (a - lo) / (hi - lo))
}

impl CirculationV2 {
    /// Derive the band geometry from planet rotation, radius and gravity.
    pub fn derive(rotation_period_s: f64, planet_radius_m: f64, k: &RotationClimateConstantsV2) -> Self {
        let g = k.surface_gravity_m_s2.value;
        // Omega / Omega_E = T_E / T; a / a_E.
        let s = (k.reference_rotation_period_s.value / rotation_period_s) * (planet_radius_m / k.reference_planet_radius_m.value);
        // g H / (g_E H_E) with H ~ g^-e: (g / g_E)^(1 - e).
        let gr = (g / k.reference_surface_gravity_m_s2.value).powf(1.0 - k.circulation_depth_gravity_exponent.value);
        let p = k.hadley_scaling_exponent.value;
        let factor = s.powf(-p) * gr.powf(0.5 * p);
        // Sine form (Mitchell and Hill 2021); the calibration planet keeps
        // its reference edge exactly.
        let phi_ref = k.reference_hadley_edge_deg.value;
        let phi_scaled = if factor == 1.0 {
            phi_ref
        } else {
            let x = phi_ref.to_radians().sin() * factor;
            if x >= 1.0 {
                90.0
            } else {
                x.asin().to_degrees()
            }
        };
        let pair = k.reference_eddy_pair_width_deg.value * s.powf(-k.eddy_pair_scaling_exponent.value);
        let w = 90.0 - phi_scaled;
        let global = w < k.global_cell_pair_fraction.value * pair;
        let (phi, pairs, (edges, kinds), blend_weight, (be, bk)) = if global {
            (90.0, 0u32, (vec![0.0, 90.0], vec![BandKindV2::Hadley]), 0.0, (Vec::new(), Vec::new()))
        } else {
            let x = w / pair;
            let m = (x.floor() as u32).max(1);
            let bw = (x - m as f64).clamp(0.0, 1.0);
            let blend = if bw > 0.0 { pair_geometry(phi_scaled, m + 1) } else { (Vec::new(), Vec::new()) };
            (phi_scaled, m, pair_geometry(phi_scaled, m), bw, blend)
        };
        Self {
            rotation_period_s,
            planet_radius_m,
            surface_gravity_m_s2: g,
            rotation_speed_ratio: s,
            depth_gravity_ratio: gr,
            hadley_edge_scaled_deg: phi_scaled,
            hadley_edge_deg: phi,
            global_hadley_cell: global,
            eddy_pair_width_deg: pair,
            eddy_pairs: pairs,
            band_edges_deg: edges,
            band_kinds: kinds,
            blend_weight,
            blend_band_edges_deg: be,
            blend_band_kinds: bk,
        }
    }

    /// (1 - w) f(primary) + w f(m + 1 geometry); the primary alone at w = 0.
    fn blended(&self, lat_deg: f64, f: impl Fn(BandKindV2, f64) -> f64) -> f64 {
        let (k0, x0) = band_of(&self.band_edges_deg, &self.band_kinds, lat_deg);
        let a = f(k0, x0);
        if self.blend_weight == 0.0 {
            return a;
        }
        let (k1, x1) = band_of(&self.blend_band_edges_deg, &self.blend_band_kinds, lat_deg);
        (1.0 - self.blend_weight) * a + self.blend_weight * f(k1, x1)
    }

    /// Prevailing zonal wind (m/s, eastward positive).
    pub fn zonal_wind_m_s(&self, lat_deg: f64, trade: f64, westerly: f64, polar: f64) -> f64 {
        self.blended(lat_deg, |kind, x| {
            let amp = match kind {
                BandKindV2::Hadley => trade,
                BandKindV2::Westerly => westerly,
                BandKindV2::Easterly => polar,
            };
            -kind.sign() * amp * x.sin()
        })
    }

    /// Prevailing meridional wind (m/s, northward positive).
    pub fn meridional_wind_m_s(&self, lat_deg: f64, hadley: f64, ferrel: f64, polar: f64) -> f64 {
        if lat_deg == 0.0 {
            return 0.0;
        }
        self.blended(lat_deg, |kind, x| {
            let amp = match kind {
                BandKindV2::Hadley => hadley,
                BandKindV2::Westerly => ferrel,
                BandKindV2::Easterly => polar,
            };
            -lat_deg.signum() * kind.sign() * amp * x.sin()
        })
    }

    /// Background rainout factor (rising at convergent band edges).
    pub fn rainout_factor(&self, lat_deg: f64, k: f64) -> f64 {
        self.blended(lat_deg, |kind, x| 1.0 + k * kind.sign() * x.cos())
    }
}

/// Recorded in the artifact payload of a `seed-climate-v2` run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RotationClimateRecordV2 {
    pub constants: RotationClimateConstantsV2,
    pub constants_digest: String,
    pub circulation: CirculationV2,
}

const V1: &[ApproximationV0] = SEED_CLIMATE_APPROXIMATIONS_V1;

/// Every approximation `seed-climate-v2` makes: v1's, with the three
/// fixed-band ones replaced by the rotation-scaled law.
pub const SEED_CLIMATE_APPROXIMATIONS_V2: &[ApproximationV0] = &[
    V1[0], // circular-orbit-annual-mean-insolation-v0
    V1[1], // budyko-sellers-zonal-ebm-v0
    V1[2], // constant-environmental-lapse-rate-v0
    V1[3], // annual-mean-no-maritime-correction-v0
    ApproximationV0 {
        id: "sine-scaled-hadley-edge-v2",
        statement: "Hadley edge sin(phi_H) = sin(phi_E) S^-p G^(p/2) (phi_H = 90 where this reaches 1), S = (Omega a) / (Omega_E a_E) from the recipe's rotation period and the grid's planet radius, G = (g / g_E)^(1 - e) (circulation depth H ~ g^-e); p = 1/3 is Mitchell and Hill (2021), p = 1/2 Held (2000) (envelope 0.25-0.5); the pole-equator temperature contrast is held at Earth's value (not read from the energy balance). Annual mean: no seasonal migration, no cross-equatorial cell.",
        constants: &["rotation_period_s", "reference_rotation_period_s", "reference_planet_radius_m", "reference_surface_gravity_m_s2", "surface_gravity_m_s2", "reference_hadley_edge_deg", "hadley_scaling_exponent", "circulation_depth_gravity_exponent"],
    },
    ApproximationV0 {
        id: "rhines-scaled-blended-eddy-bands-v2",
        statement: "Poleward of phi_H, with x = W / P (W = 90 - phi_H, pair width P = P_E S^-q: Rhines scale with Earth's eddy velocity), winds and rainout modulation are the blend (1 - w) of the m = max(1, floor(x)) pair geometry and w of the m + 1 pair geometry, w = x - m (each geometry: pairs of eddy-driven westerly and thermally direct easterly bands splitting W evenly), so the pattern changes continuously with rotation; if W < f P the Hadley cell is global (reaches the pole). Bands are hemispherically symmetric.",
        constants: &["reference_eddy_pair_width_deg", "eddy_pair_scaling_exponent", "global_cell_pair_fraction", "rotation_period_s"],
    },
    ApproximationV0 {
        id: "rotation-scaled-zonal-wind-bands-v2",
        statement: "Prevailing surface wind u = -s U_band sin(pi (|lat| - lo) / (hi - lo)) in each band [lo, hi] of the rotation-scaled geometry, s = +1 for Hadley and easterly bands (amplitudes: trade, polar easterly), -1 for westerly bands (westerly amplitude); calm at every band edge; reduces to v0's sin(6 |lat|) for edges 0, 30, 60, 90. No monsoon, seasonal shift or topographic steering.",
        constants: &["trade_wind_speed_m_s", "westerly_wind_speed_m_s", "polar_easterly_speed_m_s"],
    },
    ApproximationV0 {
        id: "rotation-scaled-meridional-wind-bands-v2",
        statement: "Prevailing surface meridional wind v = -sign(lat) s V_band sin(pi (|lat| - lo) / (hi - lo)): equatorward in Hadley and easterly bands, poleward in westerly bands, so the column converges at the equator and at each westerly-to-easterly edge and diverges at the Hadley edge, each easterly-to-westerly edge and the poles; reduces to v1's three-cell meridional bands for edges 0, 30, 60, 90.",
        constants: &["hadley_meridional_speed_m_s", "ferrel_meridional_speed_m_s", "polar_meridional_speed_m_s"],
    },
    V1[6], // two-dimensional-moisture-transport-v1 (with the v2 band wind and rainout factor)
    ApproximationV0 {
        id: "rotation-scaled-rainout-modulation-v2",
        statement: "Background rainout rate scaled by 1 + k s cos(pi (|lat| - lo) / (hi - lo)) per band: enhanced at rising branches (equator, poleward edges of westerly bands), suppressed at sinking branches (Hadley edge, poleward edges of easterly bands, poles); replaces v0's 1 + k cos(6 lat), which it equals for edges 0, 30, 60, 90.",
        constants: &["circulation_rainout_modulation"],
    },
    V1[8],  // saturation-excess-fallout-length-v1
    V1[9],  // budyko-curve-land-evapotranspiration-v1
    V1[10], // advected-ocean-influence-continentality-v1
    V1[11], // thornthwaite-1948-unadjusted-pet-v0
    V1[12], // unep-aridity-index-v0
];

/// Same schema names and units as v0/v1; v2's wind and rainout approximations.
pub const SEED_CLIMATE_FIELDS_V2: &[ClimateFieldSpecV0] = &[
    ClimateFieldSpecV0 { name: "temperature_mean_c", units: "degC", meaning: "annual-mean near-surface air temperature", status: FieldStatusV0::Modelled, approximation: "budyko-sellers-zonal-ebm-v0 + constant-environmental-lapse-rate-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "temperature_half_range_k", units: "K", meaning: "seasonal half range (warmest-month mean = mean + half range, coldest = mean - half range)", status: FieldStatusV0::Modelled, approximation: "advected-ocean-influence-continentality-v1", sentinel: None },
    ClimateFieldSpecV0 { name: "precipitation_mm_per_yr", units: "mm/yr", meaning: "annual precipitation (kg/m2/yr of water)", status: FieldStatusV0::Modelled, approximation: "two-dimensional-moisture-transport-v1 + saturation-excess-fallout-length-v1 + rotation-scaled-rainout-modulation-v2", sentinel: None },
    ClimateFieldSpecV0 { name: "precipitation_orographic_mm_per_yr", units: "mm/yr", meaning: "part of precipitation_mm_per_yr from saturation rainout (lift over relief)", status: FieldStatusV0::Modelled, approximation: "saturation-excess-fallout-length-v1", sentinel: None },
    ClimateFieldSpecV0 { name: "evaporation_mm_per_yr", units: "mm/yr", meaning: "ocean evaporation and land evapotranspiration into the column", status: FieldStatusV0::Modelled, approximation: "two-dimensional-moisture-transport-v1 + budyko-curve-land-evapotranspiration-v1", sentinel: None },
    ClimateFieldSpecV0 { name: "pet_mm_per_yr", units: "mm/yr", meaning: "potential evapotranspiration", status: FieldStatusV0::Modelled, approximation: "thornthwaite-1948-unadjusted-pet-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "aridity_index", units: "1", meaning: "P / PET", status: FieldStatusV0::Modelled, approximation: "unep-aridity-index-v0", sentinel: Some((ARIDITY_SENTINEL_PET_ZERO, "PET is zero (no month above 0 degC): aridity undefined")) },
    ClimateFieldSpecV0 { name: "wind_u_m_s", units: "m/s", meaning: "prevailing surface wind, eastward component", status: FieldStatusV0::Modelled, approximation: "rotation-scaled-zonal-wind-bands-v2 + sine-scaled-hadley-edge-v2 + rhines-scaled-blended-eddy-bands-v2", sentinel: None },
    ClimateFieldSpecV0 { name: "wind_v_m_s", units: "m/s", meaning: "prevailing surface wind, northward component", status: FieldStatusV0::Modelled, approximation: "rotation-scaled-meridional-wind-bands-v2 + sine-scaled-hadley-edge-v2 + rhines-scaled-blended-eddy-bands-v2", sentinel: None },
    ClimateFieldSpecV0 { name: "upwind_coast_distance_m", units: "m", meaning: "distance to the nearest ocean upwind along the moisture transport path", status: FieldStatusV0::Modelled, approximation: "advected-ocean-influence-continentality-v1", sentinel: Some((DISTANCE_SENTINEL_NO_UPWIND_OCEAN, "no ocean influence reaches the cell")) },
    ClimateFieldSpecV0 { name: "insolation_toa_w_m2", units: "W/m2", meaning: "annual-mean top-of-atmosphere insolation", status: FieldStatusV0::Modelled, approximation: "circular-orbit-annual-mean-insolation-v0", sentinel: None },
    ClimateFieldSpecV0 { name: "column_water_kg_m2", units: "kg/m2", meaning: "annual-mean column water vapour along the transport path", status: FieldStatusV0::Modelled, approximation: "two-dimensional-moisture-transport-v1", sentinel: None },
];

pub const SEED_CLIMATE_USAGE_LIMITS_V2: &[&str] = &[
    "coarse tier: annual means and one seasonal harmonic; no weather, no interannual variability, no monsoon or seasonal ITCZ shift",
    "band geometry is an Earth-calibrated scaling law (Hadley edge, Rhines band count), not a flow solution; band wind amplitudes stay Earth constants at every rotation rate; no flow blocking",
    "the pole-equator temperature contrast and eddy velocity in the scalings are held at Earth's values; obliquity, luminosity and atmospheric mass do not move the bands",
    "every constant is a project assumption awaiting review (ADR 039 C5); do not report against observed norms as a qualification",
    "rows[].ring_rounds is the global sweep count and rows[].ring_residual_kg_m2 the row's largest column change in the last sweep (no per-ring budget in this tier)",
];

/// The `seed-climate-v2` tier. Pure: identical inputs give identical
/// outputs; the only spatial inputs are `inputs`.
pub fn derive_seed_climate_v2(
    inputs: &ClimateInputs<'_>,
    recipe: &ClimateRecipeV1,
    extra: &SeedClimateConstantsV1,
    rot: &RotationClimateConstantsV2,
) -> Result<ClimateOutputV0, ClimateError> {
    derive_seed_climate_v2_with(inputs, recipe, extra, rot, &SeedSolveOptions::default()).map(|(out, _)| out)
}

/// CLIMATE-GENESIS (lane CG, 2026-10-20): pass-through of the solve options
/// (warm start, tolerance) and the converged column state. With
/// `SeedSolveOptions::default()` the output equals `derive_seed_climate_v2`.
pub fn derive_seed_climate_v2_with(
    inputs: &ClimateInputs<'_>,
    recipe: &ClimateRecipeV1,
    extra: &SeedClimateConstantsV1,
    rot: &RotationClimateConstantsV2,
    opts: &SeedSolveOptions<'_>,
) -> Result<(ClimateOutputV0, SeedColumnStateV1), ClimateError> {
    if recipe.law != ClimateLawV1::SeedClimateV2 {
        return Err(ClimateError::Recipe(UniverseError::Unsupported {
            law: "climate/not-seed-climate-v2",
            why: "derive_seed_climate_v2 implements seed-climate-v2 only",
        }));
    }
    recipe.validate()?;
    extra.validate()?;
    rot.validate()?;
    let circ = CirculationV2::derive(recipe.constants.rotation_period_s.value, inputs.grid.radius_m(), rot);
    let k = &recipe.constants;
    let (trade, west, polar) = (k.trade_wind_speed_m_s.value, k.westerly_wind_speed_m_s.value, k.polar_easterly_speed_m_s.value);
    let (vh, vf, vp) = (
        extra.hadley_meridional_speed_m_s.value,
        extra.ferrel_meridional_speed_m_s.value,
        extra.polar_meridional_speed_m_s.value,
    );
    let k_circ = k.circulation_rainout_modulation.value;
    let (mut out, state) = seed_climate_core_with(
        inputs,
        recipe,
        extra,
        &|l| circ.zonal_wind_m_s(l, trade, west, polar),
        &|l| circ.meridional_wind_m_s(l, vh, vf, vp),
        &|l| circ.rainout_factor(l, k_circ),
        opts,
    )?;
    for c in &mut out.closures {
        c.stage = c.stage.replace("seed-climate-v1/", "seed-climate-v2/");
    }
    out.seed_v2 = Some(RotationClimateRecordV2 {
        constants_digest: rot
            .digest()
            .map_err(|e| ClimateError::Inputs(format!("rotation climate constants digest: {e}")))?,
        constants: rot.clone(),
        circulation: circ,
    });
    Ok((out, state))
}

#[cfg(test)]
mod tests {
    //! Every map here is a SYNTHETIC FIXTURE (ADR 026 §1: legal as instrument
    //! calibration in unit tests, never world input).
    use super::*;
    use crate::climate::seed::seed_climate_v1_earth_recipe;
    use crate::climate::{derive_climate, ClimateInputs, ClimateOutputV0};
    use pf_grid::GridSpec;

    const R: f64 = 6.371e6;
    const DAY: f64 = 86_164.1;

    fn run(g: &GridSpec, z: &[f64], ocean: &[bool], recipe: &ClimateRecipeV1) -> ClimateOutputV0 {
        let inputs = ClimateInputs { grid: g, elevation_m: z, fill_level_m: z, sea_level_m: 0.0, ocean };
        derive_climate(&inputs, recipe).unwrap()
    }

    fn v2(period_s: f64) -> ClimateRecipeV1 {
        let mut r = seed_climate_v2_earth_recipe();
        r.constants.rotation_period_s.value = period_s;
        r
    }

    /// Continent from lon_lo to lon_hi below 70 deg with an optional N-S
    /// wall 6 deg wide (same fixture as tests/seed_climate_v1.rs).
    fn band_world(g: &GridSpec, lon_lo: f64, lon_hi: f64, wall: Option<f64>, wall_m: f64) -> (Vec<f64>, Vec<bool>) {
        let (nlat, nlon) = (g.nlat(), g.nlon());
        let mut z = vec![-4000.0; g.cells()];
        let mut ocean = vec![true; g.cells()];
        for i in 0..nlat {
            let lat = g.center_lat_deg(i).unwrap();
            for j in 0..nlon {
                let lon = g.center_lon_deg(j).unwrap();
                if lon >= lon_lo && lon < lon_hi && lat.abs() < 70.0 {
                    let c = i * nlon + j;
                    ocean[c] = false;
                    z[c] = 200.0;
                    if let Some(w) = wall {
                        if lon >= w && lon < w + 6.0 {
                            z[c] = 200.0 + wall_m;
                        }
                    }
                }
            }
        }
        (z, ocean)
    }

    #[test]
    fn v2_band_geometry_follows_rotation() {
        let k = RotationClimateConstantsV2::earth_v2();
        let earth = CirculationV2::derive(DAY, R, &k);
        assert_eq!(earth.band_edges_deg, vec![0.0, 30.0, 60.0, 90.0]);
        assert_eq!(earth.band_kinds, vec![BandKindV2::Hadley, BandKindV2::Westerly, BandKindV2::Easterly]);
        let slow = CirculationV2::derive(8.0 * 86_400.0, R, &k);
        assert!(slow.global_hadley_cell && slow.band_edges_deg == vec![0.0, 90.0]);
        let fast = CirculationV2::derive(8.0 * 3_600.0, R, &k);
        assert_eq!(fast.eddy_pairs, 2);
        // HADLEY-P: sin(phi) = 0.5 * 2.99^-1/3 -> about 20.3 deg (lane
        // EARTH-NUMBERS: "the 8 h world edge is ~20 deg").
        assert!(fast.hadley_edge_deg < 21.0 && fast.hadley_edge_deg > 19.5, "{}", fast.hadley_edge_deg);
        // Monotone: faster rotation never widens the Hadley cell or removes bands.
        let mut last = (91.0, 0usize);
        for hours in [2400.0, 400.0, 192.0, 96.0, 48.0, 24.0, 12.0, 8.0, 6.0, 4.0] {
            let c = CirculationV2::derive(hours * 3600.0, R, &k);
            assert!(c.hadley_edge_deg <= last.0 && c.band_kinds.len() >= last.1, "{hours} h: {c:?}");
            last = (c.hadley_edge_deg, c.band_kinds.len());
        }
        // A bigger planet at Earth's rotation rate behaves like faster rotation.
        let big = CirculationV2::derive(DAY, 2.0 * R, &k);
        assert!(big.hadley_edge_deg < 30.0);
        // Gravity cancels with the default depth exponent (H ~ 1 / g)...
        let mut kg = k.clone();
        kg.surface_gravity_m_s2.value = 20.0;
        assert_eq!(CirculationV2::derive(DAY, R, &kg).hadley_edge_deg, 30.0);
        // ...and moves the edge with a fixed depth.
        kg.circulation_depth_gravity_exponent.value = 0.0;
        assert!(CirculationV2::derive(DAY, R, &kg).hadley_edge_deg > 30.0);
    }

    /// HADLEY-P (2026-10-20): p = 1/3 default, p = 1 outside the envelope,
    /// Earth stays at 30 deg with one pair and no blend, the global cell
    /// switches at S = 0.15, and the band pattern is continuous in S.
    #[test]
    fn hadley_p_defaults_switch_and_smooth_bands() {
        let k = RotationClimateConstantsV2::earth_v2();
        assert_eq!(k.hadley_scaling_exponent.value, 1.0 / 3.0);
        assert!(k.validate().is_ok());
        let mut k1 = k.clone();
        k1.hadley_scaling_exponent.value = 1.0;
        assert!(matches!(k1.validate(), Err(UniverseError::OutOfEnvelope { .. })));
        let earth = CirculationV2::derive(DAY, R, &k);
        assert_eq!((earth.hadley_edge_deg, earth.eddy_pairs, earth.blend_weight), (30.0, 1, 0.0));
        assert!(earth.blend_band_edges_deg.is_empty());
        // Global-cell switch at S = 0.15 (rotation period / 0.15).
        assert!(CirculationV2::derive(DAY / 0.149, R, &k).global_hadley_cell);
        assert!(!CirculationV2::derive(DAY / 0.152, R, &k).global_hadley_cell);
        assert!(!CirculationV2::derive(DAY / 0.2, R, &k).global_hadley_cell);
        // Continuity across the m = 1 -> 2 step (x = W / P crosses 2): the
        // zonal wind at every latitude moves by a small amount for a small
        // change of S, where rounding jumped.
        let mut s_lo = 1.0;
        let mut s_hi = 6.0;
        for _ in 0..60 {
            let mid = 0.5 * (s_lo + s_hi);
            let c = CirculationV2::derive(DAY / mid, R, &k);
            if c.eddy_pairs >= 2 {
                s_hi = mid;
            } else {
                s_lo = mid;
            }
        }
        let a = CirculationV2::derive(DAY / s_lo, R, &k);
        let b = CirculationV2::derive(DAY / s_hi, R, &k);
        assert_eq!((a.eddy_pairs, b.eddy_pairs), (1, 2));
        assert!(a.blend_weight > 0.999 && b.blend_weight < 1e-3, "{} {}", a.blend_weight, b.blend_weight);
        for l in 1..90 {
            let lat = l as f64;
            let (ua, ub) = (a.zonal_wind_m_s(lat, -7.0, 8.0, -4.0), b.zonal_wind_m_s(lat, -7.0, 8.0, -4.0));
            assert!((ua - ub).abs() < 0.05, "lat {lat}: {ua} vs {ub}");
            let (ra, rb) = (a.rainout_factor(lat, 0.5), b.rainout_factor(lat, 0.5));
            assert!((ra - rb).abs() < 0.01, "lat {lat}: {ra} vs {rb}");
        }
    }

    /// The Earth case is seed-climate-v1, field for field.
    #[test]
    fn v2_earth_equals_v1() {
        let g = GridSpec::new(36, 72, R).unwrap();
        let (z, o) = band_world(&g, -60.0, 40.0, Some(-20.0), 3000.0);
        let a = run(&g, &z, &o, &seed_climate_v1_earth_recipe());
        let b = run(&g, &z, &o, &seed_climate_v2_earth_recipe());
        for ((n, x), (_, y)) in a.fields.named().into_iter().zip(b.fields.named()) {
            for c in 0..x.len() {
                assert!((x[c] - y[c]).abs() <= 1e-9 * x[c].abs().max(1.0), "{n}[{c}] v1 {} v2 {}", x[c], y[c]);
            }
        }
        assert!(b.closures.iter().all(|c| c.closes() && c.stage.starts_with("seed-climate-v2/")));
        assert!(b.seed_v2.is_some() && a.seed_v2.is_none());
    }

    /// Rain shadow follows the wind: at 45 deg Earth has westerlies (wet west
    /// flank), an 8 h day puts 45 deg in an easterly band (wet east flank).
    #[test]
    fn v2_wet_flank_follows_band_wind() {
        let g = GridSpec::new(90, 180, R).unwrap();
        let (z, o) = band_world(&g, -60.0, 40.0, Some(-20.0), 3000.0);
        let i = (0..g.nlat()).find(|&i| (g.center_lat_deg(i).unwrap() - 45.0).abs() < 1.1).unwrap();
        let col = |lon: f64| (0..g.nlon()).find(|&j| (g.center_lon_deg(j).unwrap() - lon).abs() < 1.1).unwrap();
        let (west, east) = (i * g.nlon() + col(-19.0), i * g.nlon() + col(-15.0));
        let earth = run(&g, &z, &o, &v2(DAY));
        let fast = run(&g, &z, &o, &v2(8.0 * 3600.0));
        let p = |o: &ClimateOutputV0, c: usize| o.fields.precipitation_mm_per_yr[c];
        assert!(earth.fields.wind_u_m_s[west] > 0.0 && fast.fields.wind_u_m_s[west] < 0.0);
        assert!(p(&earth, west) > 1.2 * p(&earth, east), "earth west {} east {}", p(&earth, west), p(&earth, east));
        assert!(p(&fast, east) > 1.2 * p(&fast, west), "fast west {} east {}", p(&fast, west), p(&fast, east));
        for out in [&earth, &fast] {
            assert!(out.closures.iter().all(|c| c.closes()));
        }
    }

    /// Envelopes: v2 accepts 8 h and 8 d; v0 and v1 still refuse 8 d; v2
    /// refuses 1 h; v0's entry point refuses the v2 law.
    #[test]
    fn v2_rotation_envelope() {
        assert!(v2(8.0 * 3600.0).validate().is_ok());
        assert!(v2(8.0 * 86_400.0).validate().is_ok());
        assert!(matches!(v2(3600.0).validate(), Err(UniverseError::OutOfEnvelope { .. })));
        let mut r1 = seed_climate_v1_earth_recipe();
        r1.constants.rotation_period_s.value = 8.0 * 86_400.0;
        assert!(matches!(r1.validate(), Err(UniverseError::OutOfEnvelope { .. })));
        let mut r0 = ClimateRecipeV1::coarse_v0_earth();
        r0.constants.rotation_period_s.value = 8.0 * 86_400.0;
        assert!(matches!(r0.validate(), Err(UniverseError::OutOfEnvelope { .. })));
        let g = GridSpec::new(18, 36, R).unwrap();
        let z = vec![-4000.0; g.cells()];
        let o = vec![true; g.cells()];
        let inputs = ClimateInputs { grid: &g, elevation_m: &z, fill_level_m: &z, sea_level_m: 0.0, ocean: &o };
        assert!(crate::climate::derive_climate_coarse_v0(&inputs, &v2(DAY)).is_err());
    }
}
