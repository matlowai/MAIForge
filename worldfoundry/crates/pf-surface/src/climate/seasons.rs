//! Opt-in diagnostic downscaling of seed-climate-v1, CX-CLIMATE-SEASONS.
//! Existing fields and the annual water solve are never modified.
//! Sources and applicability are persisted by the publisher, not implied validation.
use super::{ClimateOutputV0, SECONDS_PER_YEAR};
use pf_contracts::laws::ClimateRecipeV1;
use serde::{Deserialize, Serialize};
use std::f64::consts::{PI, TAU};

pub const SCHEMA: &str = "pf.climate-seasons/v1";
pub const MONTH_DAYS: [usize; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub maritime: bool,
    pub precipitation_seasons: bool,
    pub orographic: bool,
}
impl Options {
    pub fn any(self) -> bool {
        self.maritime || self.precipitation_seasons || self.orographic
    }
}

/// Source-inspired reduced order model, not McKinnon et al.'s trajectory solve.
/// C dT/dt + lambda T = F1 cos(omega t). C = atmosphere + mixed surface.
/// Land: 2 m water equivalent. Ocean: 50 m slab. Constant depths are assumptions.
/// Ocean fraction is the seed generator's advected/diffused exposure tracer.
/// A phase diagnostic is returned; old weather phase is NOT silently changed.
pub fn thermal_response(
    forcing_w_m2: f64,
    lambda: f64,
    coast_m: f64,
    length_m: f64,
) -> (f64, f64, f64) {
    let ocean = if coast_m < 0.0 {
        0.0
    } else {
        (-coast_m / length_m).exp()
    };
    let capacity = 1.0e7 + 4.1e6 * (2.0 + 48.0 * ocean); // J m^-2 K^-1
    let omega_c = TAU / SECONDS_PER_YEAR * capacity;
    let half = forcing_w_m2 / lambda.hypot(omega_c);
    let lag_days = omega_c.atan2(lambda) / TAU * 365.0;
    (half, lag_days, ocean)
}
fn smooth(a: f64, b: f64, x: f64) -> f64 {
    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Regime closure, NOT a dynamical monthly moisture solver or a Koppen classifier.
/// Summer convection on continents; winter storms on maritime midlatitudes;
/// subtropical summer drying strongest where westerly marine exposure exists.
/// Coefficients are declared assumptions informed by station normals and Koppen
/// season definitions. No monsoon, ENSO, storm track longitude or rain/snow drift.
pub fn winter_target(lat_deg: f64, ocean: f64) -> f64 {
    let lat = lat_deg.abs();
    let tropical = 0.5 - 0.25 * smooth(0.0, 15.0, lat);
    let mediterranean = smooth(25.0, 32.0, lat) * (1.0 - smooth(38.0, 48.0, lat));
    let extratropical = (1.0 - ocean) * 0.30 + ocean * (0.60 + 0.20 * mediterranean);
    let w = smooth(20.0, 35.0, lat);
    (1.0 - w) * tropical + w * extratropical
}
pub fn winter_month(north: bool, m: usize) -> bool {
    if north {
        m >= 9 || m <= 2
    } else {
        (3..=8).contains(&m)
    }
}

/// Normalize each half-year separately to make the published cold-half fraction
/// exact (Oct-Mar N, Apr-Sep S). Month shares are annual totals, not daily rates.
pub fn monthly_shares(lat_deg: f64, winter: f64) -> [f64; 12] {
    let mut a = [0.0; 12];
    let mut sums = [0.0; 2];
    let mut d0 = 0;
    let peak = if lat_deg >= 0.0 { 200.0 } else { 17.5 };
    for m in 0..12 {
        for d in d0..d0 + MONTH_DAYS[m] {
            a[m] += 1.0 - PI * (winter - 0.5) * (TAU * (d as f64 - peak) / 365.0).cos();
        }
        sums[winter_month(lat_deg >= 0.0, m) as usize] += a[m];
        d0 += MONTH_DAYS[m];
    }
    for m in 0..12 {
        let cold = winter_month(lat_deg >= 0.0, m);
        a[m] *= if cold {
            winter / sums[1]
        } else {
            (1.0 - winter) / sums[0]
        };
    }
    a
}

/// Frozen-column extension of the EXACT seed saturation-excess fallout law.
/// P(z)=P_ref+r{[q-s_ref exp(-k(z-z_ref))]+-[q-s_ref]+}.
/// r = seconds_per_year / fallout_time, k = CC rate * lapse.
/// Reference precipitation is at max(T0 ground, filled level), NOT sea level.
/// Finite moisture supply bounds P. No exponential precipitation growth.
/// This is a column diagnostic, not area-conservative subgrid hydrology:
/// hydrology adoption requires area normalization with actual terrain hypsometry.
pub fn precipitation_at_height(
    p_ref: f64,
    q: f64,
    s_ref: f64,
    r: f64,
    k: f64,
    z_ref: f64,
    z: f64,
) -> f64 {
    p_ref + r * ((q - s_ref * (-k * (z - z_ref)).exp()).max(0.0) - (q - s_ref).max(0.0))
}

#[derive(Debug)]
pub struct Field {
    pub name: String,
    pub units: &'static str,
    pub meaning: &'static str,
    pub data: Vec<f64>,
}

/// Derived new fields only. Off returns no fields and does not touch output.
/// Temperatures and saturation are evaluated from the same unquantized solve.
pub fn derive(
    out: &ClimateOutputV0,
    recipe: &ClimateRecipeV1,
    nlon: usize,
    z_ref: &[f64],
    opt: Options,
) -> Result<Vec<Field>, String> {
    if !opt.any() {
        return Ok(vec![]);
    }
    let n = out.fields.temperature_mean_c.len();
    if nlon == 0
        || out.rows.len() * nlon != n
        || z_ref.len() != n
        || z_ref.iter().any(|x| !x.is_finite())
    {
        return Err("seasons grid/reference-height mismatch".into());
    }
    if recipe.law != pf_contracts::laws::ClimateLawV1::SeedClimateV1 {
        return Err("seasons requires seed-climate-v1".into());
    }
    let k = &recipe.constants;
    let extra = out.seed_v1.as_ref().ok_or("missing seed constants")?;
    if extra.constants.crh_rainout.is_some() {
        return Err("seasons v1 not validated with crh rainout".into());
    }
    let rate = SECONDS_PER_YEAR * k.moisture_transport_speed_m_s.value
        / extra.constants.orographic_fallout_length_m.value;
    let cc = k.clausius_clapeyron_rate_per_k.value;
    let lapse = k.lapse_rate_k_per_m.value;
    let mut fields = Vec::new();
    let mut add = |name: &str, units, meaning, data| {
        fields.push(Field {
            name: name.into(),
            units,
            meaning,
            data,
        })
    };
    let thermal: Vec<_> = (0..n)
        .map(|c| {
            let row = &out.rows[c / nlon];
            thermal_response(
                (1.0 - row.albedo) * row.insolation_first_harmonic_w_m2,
                k.olr_slope_w_m2_k.value + k.meridional_transport_w_m2_k.value,
                out.fields.upwind_coast_distance_m[c],
                k.continentality_length_m.value,
            )
        })
        .collect();
    if opt.maritime {
        add("temperature_half_range_maritime_k","K","annual first harmonic with mixed heat capacity; old temperature_half_range_k unchanged",thermal.iter().map(|x|x.0).collect());
        add(
            "temperature_response_lag_days",
            "day",
            "phase lag after insolation maximum; diagnostic only; existing weather phase unchanged",
            thermal.iter().map(|x| x.1).collect(),
        );
        add(
            "maritime_exposure_fraction",
            "1",
            "advected and diffused ocean tracer recovered from effective upwind distance",
            thermal.iter().map(|x| x.2).collect(),
        );
    }
    if opt.precipitation_seasons {
        let winter: Vec<_> = (0..n)
            .map(|c| winter_target(out.rows[c / nlon].lat_deg, thermal[c].2))
            .collect();
        let shares: Vec<_> = (0..n)
            .map(|c| monthly_shares(out.rows[c / nlon].lat_deg, winter[c]))
            .collect();
        add(
            "precipitation_winter_fraction",
            "1",
            "Oct-Mar fraction north, Apr-Sep south; annual P unchanged",
            winter,
        );
        for m in 0..12 {
            add(&format!("precipitation_month_{:02}_fraction",m+1),"1","fraction of annual precipitation in calendar month; divide by days in month for daily forcing",shares.iter().map(|s|s[m]).collect());
        }
    }
    if opt.orographic {
        let sat: Vec<_> = out
            .fields
            .temperature_mean_c
            .iter()
            .map(|&t| {
                k.column_water_saturation_ref_kg_m2.value
                    * (cc * (t - k.column_water_reference_temperature_c.value)).exp()
            })
            .collect();
        let grad: Vec<_> = (0..n)
            .map(|c| {
                if out.fields.column_water_kg_m2[c] >= sat[c] {
                    rate * cc * lapse * sat[c]
                } else {
                    0.0
                }
            })
            .collect();
        add(
            "precipitation_reference_height_m",
            "m",
            "max(T0 elevation,filled water level) minus sea; marine reference zero",
            z_ref.to_vec(),
        );
        add(
            "precipitation_reference_saturation_kg_m2",
            "kg/m2",
            "saturation column at the T0 temperature used in the annual solve",
            sat,
        );
        add(
            "precipitation_fallout_per_year",
            "1/yr",
            "seconds per year divided by fallout residence time",
            vec![rate; n],
        );
        add(
            "precipitation_saturation_decay_per_m",
            "1/m",
            "Clausius Clapeyron rate times environmental lapse",
            vec![cc * lapse; n],
        );
        add("precipitation_height_gradient_mm_yr_per_m","mm/yr/m","local right derivative at reference height; use full bounded saturation law away from reference",grad);
    }
    for f in &fields {
        if f.data.iter().any(|x| !x.is_finite()) {
            return Err(format!("non-finite seasons field {}", f.name));
        }
    }
    Ok(fields)
}

/// Typed persisted extension. Old domain payloads omit this record exactly.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThermalConstants {
    pub atmosphere_j_m2_k: f64,
    pub water_volumetric_capacity_j_m3_k: f64,
    pub land_water_equivalent_depth_m: f64,
    pub ocean_mixed_layer_depth_m: f64,
    pub source_class: String,
    pub ocean_depth_sensitivity_m: [f64; 2],
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub schema: String,
    pub options: Options,
    pub annual_base_digest: String,
    pub sea_level_m: f64,
    pub legacy_fields: String,
    pub thermal_law: String,
    pub thermal_constants: ThermalConstants,
    pub thermal_sources: Vec<String>,
    pub rain_season_law: String,
    pub rain_season_sources: Vec<String>,
    pub orographic_law: String,
    pub orographic_sources: Vec<String>,
    pub limits: Vec<String>,
}
