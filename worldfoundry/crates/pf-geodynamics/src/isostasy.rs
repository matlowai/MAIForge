//! Airy isostasy + compensation dispatch (brief M2). Crust thickness ->
//! elevation against a reference column; `Airy` XOR `Flexure`, never both
//! (review #5). Flexure lives in `flexure.rs`; this module dispatches.

use crate::flexure::FlexureSolver;
use pf_contracts::laws::{IsostasyLawV1, SeaLevelInputV1, UniverseV1};
use pf_contracts::params::{Compensation, PlanetParameterProfile};

#[derive(Debug, Clone, PartialEq)]
pub enum IsostasyError {
    Params(String),
    NonFinite {
        what: &'static str,
        value: f64,
    },
    CompensationConflict {
        mode: Compensation,
        why: &'static str,
    },
}

impl std::fmt::Display for IsostasyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IsostasyError::Params(why) => write!(f, "parameters: {why}"),
            IsostasyError::NonFinite { what, value } => write!(f, "nonfinite {what}: {value:?}"),
            IsostasyError::CompensationConflict { mode, why } => {
                write!(f, "compensation mode {mode:?}: {why}")
            }
        }
    }
}

impl std::error::Error for IsostasyError {}

/// Reference continental crust thickness (m): the column with elevation 0.
/// WS-CONST: the legacy law reads this value from
/// `UniverseV1` (`pf_contracts::laws::legacy::REFERENCE_CRUST_THICKNESS_M`,
/// equality pinned by pf-lab tests/ws_const_preview.rs); this constant is
/// kept for the flexure kernel and existing callers.
pub const REFERENCE_CRUST_M: f64 = 35_000.0;

/// WS-CONST: resolved constants of the legacy single-column law, read from
/// the `UniverseV1` (never from the parameter blocks directly).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SingleColumnConstants {
    pub mantle_density_kg_m3: f64,
    /// One density for every crust kind under this law.
    pub crust_density_kg_m3: f64,
    pub reference_crust_thickness_m: f64,
    pub compensation: Compensation,
}

impl SingleColumnConstants {
    /// Refuses unless the universe's isostasy law is the single-column one.
    pub fn from_universe(u: &UniverseV1) -> Result<Self, IsostasyError> {
        match &u.law_set.isostasy {
            IsostasyLawV1::SingleColumnReferenceV0 {
                reference_crust_thickness_m,
                compensation,
            } => Ok(Self {
                mantle_density_kg_m3: u.constants.mantle_density_kg_m3.value,
                crust_density_kg_m3: u.constants.continental_crust_density_kg_m3.value,
                reference_crust_thickness_m: reference_crust_thickness_m.value,
                compensation: *compensation,
            }),
            other => Err(IsostasyError::Params(format!(
                "single-column isostasy refused: the law set's isostasy law is {other:?}"
            ))),
        }
    }
}

/// Single-column Airy elevation from resolved constants (the arithmetic of
/// `airy_elevation_m`, unchanged): `(rho_m - rho_c)/rho_m * (h - h_ref)`.
pub fn single_column_airy_elevation_m(
    crust_thickness_m: f64,
    k: &SingleColumnConstants,
) -> Result<f64, IsostasyError> {
    if !crust_thickness_m.is_finite() || crust_thickness_m <= 0.0 {
        return Err(IsostasyError::NonFinite {
            what: "crust_thickness_m",
            value: crust_thickness_m,
        });
    }
    let rho_m = k.mantle_density_kg_m3;
    let rho_c = k.crust_density_kg_m3;
    Ok((rho_m - rho_c) / rho_m * (crust_thickness_m - k.reference_crust_thickness_m))
}

/// Under the opt-in two-column mode the single-column reference functions
/// are refused loudly: they cannot know crust kind or the solved sea
/// level, and silently running the legacy model would substitute physics.
fn refuse_if_two_column(params: &PlanetParameterProfile) -> Result<(), IsostasyError> {
    if params.two_column_isostasy.is_some() {
        return Err(IsostasyError::CompensationConflict {
            mode: params.mantle.compensation,
            why: "two_column_isostasy is set: per-cell single-column isostasy against the \
                  35 km reference is refused; use isostasy::two_column_field",
        });
    }
    Ok(())
}

/// Airy elevation of a crust column vs the reference column (M2):
/// `e = (rho_m - rho_c)/rho_m * (h_c - h_ref)` — thicker crust rides
/// higher; zero at the reference thickness; linear in thickness. Densities
/// come from the parameter profile (fail-loud validated).
pub fn airy_elevation_m(
    crust_thickness_m: f64,
    params: &PlanetParameterProfile,
) -> Result<f64, IsostasyError> {
    params
        .validate()
        .map_err(|e| IsostasyError::Params(e.to_string()))?;
    refuse_if_two_column(params)?;
    let u = params
        .universe()
        .map_err(|e| IsostasyError::Params(e.to_string()))?;
    single_column_airy_elevation_m(crust_thickness_m, &SingleColumnConstants::from_universe(&u)?)
}

/// Compensation dispatch (M2, review #5): Airy computes directly; Flexure
/// REQUIRES a solver — Airy with a solver present, or Flexure without one,
/// are both loud errors (never silently double-compensate or substitute).
pub fn isostatic_elevation_m(
    crust_thickness_m: f64,
    params: &PlanetParameterProfile,
    flexure: Option<&FlexureSolver>,
) -> Result<f64, IsostasyError> {
    params
        .validate()
        .map_err(|e| IsostasyError::Params(e.to_string()))?;
    refuse_if_two_column(params)?;
    match params.mantle.compensation {
        Compensation::Airy => {
            if flexure.is_some() {
                return Err(IsostasyError::CompensationConflict {
                    mode: Compensation::Airy,
                    why: "Airy mode with a flexure solver supplied (double compensation?)",
                });
            }
            airy_elevation_m(crust_thickness_m, params)
        }
        Compensation::Flexure => {
            let solver = flexure.ok_or(IsostasyError::CompensationConflict {
                mode: Compensation::Flexure,
                why: "Flexure mode requires the flexure solver (never silently Airy)",
            })?;
            solver
                .airy_equivalent_elevation_m(crust_thickness_m, params)
                .map_err(|e| IsostasyError::Params(e.to_string()))
        }
    }
}

// ---------------------------------------------------------------------------
// OPT-IN two-column isostasy + volume-primary sea level (FREEBOARD-SPIKE,
// lane F, 2026-10-11; WORLD_SIMULATION_EXECPLAN §5). Active only when
// `params.two_column_isostasy` is `Some`; nothing above changes otherwise.
//
// Model (local Airy compensation, SI units throughout):
//
// - Datum ("gauge"): the surface height a crust-free, water-free column of
//   pure mantle would have. Elevations below are on this gauge until the
//   final re-reference to sea level.
// - Every column has the same mass per unit area above the common
//   compensation depth D (m below the gauge zero) as that pure-mantle
//   column: M_ref = rho_m * D.
// - A dry column (crust h_c at rho_c by kind, sediment h_s at rho_s,
//   mantle filling from the crust base down to -D) with surface e and a
//   thermal buoyancy t_dry (a mantle-lithosphere density deficit written
//   as an equivalent mantle thickness) satisfies
//       rho_c h_c + rho_s h_s + rho_m (D + e - h_c - h_s) - rho_m t_dry = rho_m D.
// - A submerged column adds a water layer rho_w (S - e). Solving gives
//       e_wet = (rho_m e_dry - rho_w S) / (rho_m - rho_w),
//   i.e. the water depth is phi (S - e_dry) with phi = rho_m/(rho_m - rho_w).
// - The ocean volume is the INPUT: V(S) = sum_i A_i max(S - e_i(S), 0)
//   = phi sum_i A_i max(S - e_dry_i, 0), with A_i true spherical cell
//   areas. V is continuous, piecewise linear and non-decreasing in S, so S
//   is solved exactly segment by segment. Sea level and flooded fraction
//   are OUTPUTS.
// - GDH1 thermal subsidence is observed under water, so the caller passes
//   it as a water-loaded height; it is converted to dry-equivalent with
//   (rho_m - rho_w)/rho_m, which makes a submerged cell's thermal relief
//   exactly the GDH1 value again.
//
// Declared limits: global fill (every cell below S is flooded, even a
// basin not connected to the main ocean; the review's connected-ocean
// caveat is MEASURED by callers, not solved here); no sea-level-dependent
// lithosphere thermal state for continents (continental thermal term 0,
// isopycnic assumption); local compensation only (no flexure).
// ---------------------------------------------------------------------------

use pf_contracts::terrain::CrustKindCodeV0;

/// One column's inputs (SI units).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ColumnV0 {
    pub kind: CrustKindCodeV0,
    /// Crystalline crust thickness (m), > 0.
    pub crust_m: f64,
    /// Sediment thickness on top of the crust (m), >= 0.
    pub sediment_m: f64,
    /// Thermal uplift relative to old lithosphere, expressed as a
    /// WATER-LOADED height (m) as GDH1 reports it; 0 for continents.
    pub thermal_uplift_water_loaded_m: f64,
}

/// Resolved constants for the two-column model (all SI).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TwoColumnConstants {
    pub mantle_density_kg_m3: f64,
    pub continental_crust_density_kg_m3: f64,
    pub oceanic_crust_density_kg_m3: f64,
    pub sediment_density_kg_m3: f64,
    pub seawater_density_kg_m3: f64,
    pub compensation_depth_m: f64,
    pub ocean_water_volume_m3: f64,
    /// Downstream datum the solved sea level is placed at (m).
    pub output_sea_level_datum_m: f64,
}

impl TwoColumnConstants {
    /// Resolve from a validated profile. Refuses when the flag is absent
    /// (never silently legacy) or when compensation is Flexure (the
    /// two-column model is local; flexure + two-column is not built).
    pub fn from_params(params: &PlanetParameterProfile) -> Result<Self, IsostasyError> {
        params
            .validate()
            .map_err(|e| IsostasyError::Params(e.to_string()))?;
        // Flag-absent and flexure refusals keep their pre-WS-CONST messages.
        params
            .two_column_isostasy
            .as_ref()
            .ok_or(IsostasyError::CompensationConflict {
                mode: params.mantle.compensation,
                why: "two_column_field called without params.two_column_isostasy (opt-in flag absent)",
            })?;
        if params.mantle.compensation != Compensation::Airy {
            return Err(IsostasyError::CompensationConflict {
                mode: params.mantle.compensation,
                why: "two-column isostasy is local (Airy) compensation; flexure is not built for it",
            });
        }
        // WS-CONST: every value is read through the resolved universe; the
        // fixed-depth mode uses the law's fixed compensation depth.
        let u = params
            .universe()
            .map_err(|e| IsostasyError::Params(e.to_string()))?;
        let depth = match &u.law_set.isostasy {
            IsostasyLawV1::TwoColumnCommonDepthV0 {
                fixed_compensation_depth_m,
                ..
            } => fixed_compensation_depth_m.value,
            other => {
                return Err(IsostasyError::Params(format!(
                    "two-column constants refused: the law set's isostasy law is {other:?}"
                )))
            }
        };
        Self::from_universe(&u, depth)
    }

    /// WS-CONST: resolve from a validated `UniverseV1` with an explicit
    /// compensation depth (the law's fixed depth, or the preview chain's
    /// derived deepest-root-plus-margin depth). Refuses unless the law set
    /// is the two-column law with the volume-primary sea level.
    pub fn from_universe(
        u: &UniverseV1,
        compensation_depth_m: f64,
    ) -> Result<Self, IsostasyError> {
        if !matches!(u.law_set.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. }) {
            return Err(IsostasyError::Params(format!(
                "two-column constants refused: the law set's isostasy law is {:?}",
                u.law_set.isostasy
            )));
        }
        let k = &u.constants;
        let oceanic = k.oceanic_crust_density_kg_m3.as_ref().ok_or_else(|| {
            IsostasyError::Params("two-column law without oceanic_crust_density_kg_m3".into())
        })?;
        let volume = match u
            .sea_level_input()
            .map_err(|e| IsostasyError::Params(e.to_string()))?
        {
            SeaLevelInputV1::OceanVolume { volume_m3, .. } => volume_m3,
            SeaLevelInputV1::Prescribed { .. } => {
                return Err(IsostasyError::Params(
                    "two-column law requires the volume-primary sea-level law".into(),
                ))
            }
        };
        Ok(Self {
            mantle_density_kg_m3: k.mantle_density_kg_m3.value,
            continental_crust_density_kg_m3: k.continental_crust_density_kg_m3.value,
            oceanic_crust_density_kg_m3: oceanic.value,
            // Legacy assumption A-SEDIMENT-DENSITY kept: sediment loads at
            // the (continental) crust density (a law rule in UniverseV1).
            sediment_density_kg_m3: u.sediment_density_kg_m3(),
            seawater_density_kg_m3: k.seawater_density_kg_m3.value,
            compensation_depth_m,
            ocean_water_volume_m3: volume,
            output_sea_level_datum_m: k.sea_level_datum_m.value,
        })
    }

    /// Water-loading amplification phi = rho_m / (rho_m - rho_w).
    pub fn water_load_factor(&self) -> f64 {
        self.mantle_density_kg_m3 / (self.mantle_density_kg_m3 - self.seawater_density_kg_m3)
    }

    fn crust_density(&self, kind: CrustKindCodeV0) -> f64 {
        match kind {
            CrustKindCodeV0::Continental => self.continental_crust_density_kg_m3,
            CrustKindCodeV0::Oceanic => self.oceanic_crust_density_kg_m3,
        }
    }
}

/// Unloaded (dry) surface elevation of one column on the mantle
/// free-surface gauge, by mass balance above the compensation depth.
pub fn dry_elevation_m(col: &ColumnV0, k: &TwoColumnConstants) -> Result<f64, IsostasyError> {
    if !col.crust_m.is_finite() || col.crust_m <= 0.0 {
        return Err(IsostasyError::NonFinite {
            what: "crust_m",
            value: col.crust_m,
        });
    }
    if !col.sediment_m.is_finite() || col.sediment_m < 0.0 {
        return Err(IsostasyError::NonFinite {
            what: "sediment_m",
            value: col.sediment_m,
        });
    }
    if !col.thermal_uplift_water_loaded_m.is_finite() {
        return Err(IsostasyError::NonFinite {
            what: "thermal_uplift_water_loaded_m",
            value: col.thermal_uplift_water_loaded_m,
        });
    }
    let rho_m = k.mantle_density_kg_m3;
    let d = k.compensation_depth_m;
    let h = col.crust_m + col.sediment_m;
    let m_ref = rho_m * d;
    let m_layers =
        k.crust_density(col.kind) * col.crust_m + k.sediment_density_kg_m3 * col.sediment_m;
    let thermal_dry = col.thermal_uplift_water_loaded_m * (rho_m - k.seawater_density_kg_m3) / rho_m;
    // rho_m (D + e - h) + m_layers - rho_m t_dry = m_ref
    let e = (m_ref - m_layers + rho_m * thermal_dry) / rho_m - d + h;
    let root = e - h;
    if !(root > -d) {
        return Err(IsostasyError::Params(format!(
            "crustal root at {root} m lies at or below the compensation depth -{d} m"
        )));
    }
    Ok(e)
}

/// Water-loaded surface elevation on the gauge for sea level `s` (gauge).
pub fn loaded_elevation_m(e_dry: f64, s: f64, k: &TwoColumnConstants) -> f64 {
    if e_dry >= s {
        e_dry
    } else {
        let rho_m = k.mantle_density_kg_m3;
        let rho_w = k.seawater_density_kg_m3;
        (rho_m * e_dry - rho_w * s) / (rho_m - rho_w)
    }
}

/// Water volume (m^3) held below sea level `s` by the water-loaded surface:
/// V(S) = phi * sum_i A_i max(S - e_dry_i, 0).
pub fn ocean_volume_m3(e_dry: &[f64], areas_m2: &[f64], s: f64, k: &TwoColumnConstants) -> f64 {
    let phi = k.water_load_factor();
    e_dry
        .iter()
        .zip(areas_m2)
        .map(|(&e, &a)| if e < s { a * phi * (s - e) } else { 0.0 })
        .sum()
}

/// Solve the sea level S (gauge) that holds exactly `volume_m3` of water.
/// Exact on each linear segment of V(S); zero volume gives the deepest
/// dry elevation; a volume above every column floods the whole planet.
pub fn solve_sea_level_m(
    e_dry: &[f64],
    areas_m2: &[f64],
    volume_m3: f64,
    k: &TwoColumnConstants,
) -> Result<f64, IsostasyError> {
    if e_dry.is_empty() || e_dry.len() != areas_m2.len() {
        return Err(IsostasyError::Params(format!(
            "sea-level solve needs equal non-empty inputs ({} elevations, {} areas)",
            e_dry.len(),
            areas_m2.len()
        )));
    }
    if !volume_m3.is_finite() || volume_m3 < 0.0 {
        return Err(IsostasyError::NonFinite {
            what: "ocean_water_volume_m3",
            value: volume_m3,
        });
    }
    for (&e, &a) in e_dry.iter().zip(areas_m2) {
        if !e.is_finite() {
            return Err(IsostasyError::NonFinite {
                what: "dry elevation",
                value: e,
            });
        }
        if !a.is_finite() || a <= 0.0 {
            return Err(IsostasyError::NonFinite {
                what: "cell area",
                value: a,
            });
        }
    }
    let phi = k.water_load_factor();
    let mut idx: Vec<usize> = (0..e_dry.len()).collect();
    idx.sort_by(|&a, &b| e_dry[a].total_cmp(&e_dry[b]).then(a.cmp(&b)));
    let mut e_prev = e_dry[idx[0]];
    if volume_m3 == 0.0 {
        return Ok(e_prev);
    }
    let mut v_prev = 0.0;
    let mut a_wet = 0.0;
    for &i in &idx {
        let e_i = e_dry[i];
        if a_wet > 0.0 {
            let v_i = v_prev + phi * a_wet * (e_i - e_prev);
            if v_i >= volume_m3 {
                return Ok(e_prev + (volume_m3 - v_prev) / (phi * a_wet));
            }
            v_prev = v_i;
        }
        e_prev = e_i;
        a_wet += areas_m2[i];
    }
    Ok(e_prev + (volume_m3 - v_prev) / (phi * a_wet))
}

/// Result of the two-column field solve.
#[derive(Debug, Clone, PartialEq)]
pub struct TwoColumnSolution {
    /// Water-loaded elevation relative to the SOLVED sea level, placed at
    /// `planet.sea_level_m` (the datum every downstream stage uses), m.
    pub elevation_m: Vec<f64>,
    /// Dry elevation on the mantle free-surface gauge (m).
    pub dry_elevation_gauge_m: Vec<f64>,
    /// Solved sea level on the mantle free-surface gauge (m).
    pub sea_level_gauge_m: f64,
    /// Area fraction with dry elevation below the solved sea level (global
    /// fill; includes basins not connected to the main ocean).
    pub flooded_area_fraction: f64,
    /// V(S) recomputed at the solved S (m^3), for a closure check.
    pub volume_check_m3: f64,
}

/// Two-column isostasy over a whole field, then the volume-primary sea
/// level. Requires `params.two_column_isostasy` (opt-in) and Airy.
pub fn two_column_field(
    columns: &[ColumnV0],
    areas_m2: &[f64],
    params: &PlanetParameterProfile,
) -> Result<TwoColumnSolution, IsostasyError> {
    let k = TwoColumnConstants::from_params(params)?;
    two_column_field_with(columns, areas_m2, &k)
}

/// WS-CONST: `two_column_field` with already-resolved constants (the
/// preview chain resolves them from `UniverseV1` with its derived depth).
pub fn two_column_field_with(
    columns: &[ColumnV0],
    areas_m2: &[f64],
    k: &TwoColumnConstants,
) -> Result<TwoColumnSolution, IsostasyError> {
    let k = *k;
    if columns.len() != areas_m2.len() {
        return Err(IsostasyError::Params(format!(
            "{} columns but {} areas",
            columns.len(),
            areas_m2.len()
        )));
    }
    let mut e_dry = Vec::with_capacity(columns.len());
    for (c, col) in columns.iter().enumerate() {
        e_dry.push(dry_elevation_m(col, &k).map_err(|e| {
            IsostasyError::Params(format!("two-column isostasy at cell {c}: {e}"))
        })?);
    }
    let s = solve_sea_level_m(&e_dry, areas_m2, k.ocean_water_volume_m3, &k)?;
    let phi = k.water_load_factor();
    let datum = k.output_sea_level_datum_m;
    let (mut flooded, mut total) = (0.0, 0.0);
    let mut elevation = Vec::with_capacity(columns.len());
    for (&e, &a) in e_dry.iter().zip(areas_m2) {
        total += a;
        if e < s {
            flooded += a;
            elevation.push(datum - phi * (s - e));
        } else {
            elevation.push(datum + (e - s));
        }
    }
    let volume_check_m3 = ocean_volume_m3(&e_dry, areas_m2, s, &k);
    Ok(TwoColumnSolution {
        elevation_m: elevation,
        dry_elevation_gauge_m: e_dry,
        sea_level_gauge_m: s,
        flooded_area_fraction: flooded / total,
        volume_check_m3,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_contracts::params::{Compensation, FidelityPreset, MantleBlock, ScenarioPreset};

    fn params(mode: Compensation) -> PlanetParameterProfile {
        let mut p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
        p.mantle.compensation = mode;
        p
    }

    #[test]
    fn airy_reference_column_is_zero_and_thickness_is_monotone() {
        let p = params(Compensation::Airy);
        assert!(airy_elevation_m(REFERENCE_CRUST_M, &p).unwrap().abs() < 1e-9);
        let thin = airy_elevation_m(REFERENCE_CRUST_M * 0.5, &p).unwrap();
        let thick = airy_elevation_m(REFERENCE_CRUST_M * 1.5, &p).unwrap();
        assert!(thin < 0.0 && thick > 0.0);
        // Value pinned: (3300-2800)/3300 * 17500 m.
        assert!((thick - (500.0 / 3300.0) * 17_500.0).abs() < 1e-6);
        assert!(airy_elevation_m(-1.0, &p).is_err());
        assert!(airy_elevation_m(f64::NAN, &p).is_err());
    }

    #[test]
    fn dispatch_is_xor_and_loud() {
        let solver = FlexureSolver::reference().unwrap();
        // Airy + a solver present = suspicious double-compensation -> loud.
        assert!(isostatic_elevation_m(
            REFERENCE_CRUST_M,
            &params(Compensation::Airy),
            Some(&solver)
        )
        .is_err());
        assert!(
            isostatic_elevation_m(REFERENCE_CRUST_M, &params(Compensation::Airy), None)
                .unwrap()
                .abs()
                < 1e-9
        );
        // Flexure without a solver -> loud (never silently Airy).
        assert!(
            isostatic_elevation_m(REFERENCE_CRUST_M, &params(Compensation::Flexure), None).is_err()
        );
        // Flexure WITH a solver: near-reference thickness -> near-zero
        // elevation (the solver's compensation volume normalizes to Airy).
        let e = isostatic_elevation_m(
            REFERENCE_CRUST_M,
            &params(Compensation::Flexure),
            Some(&solver),
        )
        .unwrap();
        assert!(
            e.abs() < 1.0,
            "reference column through flexure ~ sea level, got {e}"
        );
    }
}
