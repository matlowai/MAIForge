//! CLIMATE-ARIDITY (lane CLIMATE-ARIDITY, 2026-10-21): the opt-in
//! `crh-rainout-v1` option of the `seed-climate-v1` / `seed-climate-v2`
//! moisture solver, selected by the seed-constants set
//! `SeedClimateConstantsV1::earth_crh_v1()` (schema
//! `SEED_CLIMATE_CONSTANTS_SCHEMA_CRH1`). The `earth_v1()` set, its digest
//! and every artifact made with it are unchanged.
//!
//! PROBLEM (measured, `unpublished-design-record`): with the v1
//! rainout P = f(lat) W / tau_r a column's precipitation depends on its
//! water W only, so a continental column that still holds about half the
//! ocean's water rains about half the ocean rate. Land P was nearly zonal;
//! subtropical land had P / W ~ 28 /yr, where Earth's subtropical deserts
//! hold W ~ 15-25 kg/m2 yet get P < 100 mm/yr (P / W of order 1-5 /yr;
//! recalled). Arid + hyper-arid land (UNEP AI < 0.2) was 2.6 % on C6 and
//! 2.9 % on a C7 genesis snapshot; Earth (UNEP 1992, Thornthwaite PET)
//! 19.6 %. Uniform PET scaling cannot fix it (PET x2 gives drylands 63-71 %
//! with arid still 15-18 %): the AI distribution is too narrow because the
//! P field lacks the dry tail.
//!
//! LAW. Observed precipitation rises exponentially with the column relative
//! humidity r = W / W_sat (Bretherton, Peters and Back 2004, J. Climate
//! 17:1517-1528: monthly means over all tropical oceans,
//! P = exp(11.4 (r - 0.522)) mm/day; daily means a = 15.6, r = 0.603; the
//! slope falls as the averaging period grows). Large-scale subsidence (the
//! subtropical highs, interiors downstream of rain-out) holds r low and P
//! collapses although W is not small. The option multiplies the v1
//! background rainout by a logistic efficiency whose dry tail is that
//! observed exponential and which tends to the v1 rate in moist columns:
//!   P_bg = f(lat) W / tau_r * G(r),  G(r) = 1 / (1 + exp(-a (r - r_0))).
//! a = 11.4 (the monthly fit, the longest averaging published). r_0 = 0.55:
//! for r << r_0, P_bg ~ (f W_sat r / tau_r) exp(a (r - r_0)); equating this
//! to the monthly fit at a dry column (r 0.45, f 1) gives r_0 = 0.56 at
//! 15 degC to 0.63 at 25 degC, so the dry-tail rate is the observed one to
//! within a factor ~2; inside that bracket 0.55 was chosen from measured
//! aridity (0.6 put genesis arid + hyper-arid at 25.7 % against Earth's
//! 19.6 %; see the constant's source label). The wet side is not amplified beyond the v1 rate: the
//! coarse grid cannot resolve the convective pick-up (an amplified variant
//! put ~1e5 mm/yr in single coastal cells). A first variant with
//! min(1, exp(a (r - r_0))) had a slope kink at r_0 and failed to converge
//! (40,000 sweeps) on a 24 x 48 fixture; the logistic is smooth.
//! Saturation-excess (orographic) fallout is forced ascent and unchanged.
//!
//! CONSERVATION. The rainout is a sink of the same transported column; the
//! global water ledger closes to `SEED_WATER_REL_TOL` on every solve (typed
//! failure otherwise). The steep sink is Newton-linearized inside the line
//! Gauss-Seidel solve: lagging it like v1's other nonlinear terms would
//! oscillate (the lagged map has slope a r ~ 7 > 1).
//!
//! LIMITS: the fit is for tropical oceans; it is applied at all latitudes
//! and over land (project assumption awaiting review; extratropical frontal
//! rain is probably less sensitive to r). r uses the saturation column of
//! the annual-mean surface temperature, not a humidity profile. Moist
//! columns rain at most at the v1 rate, so global (ocean) precipitation
//! falls: C6 global P 1049 -> 887 mm/yr (Earth ~1000, recalled), land P
//! 1042 -> 957. No subtropical-high zonal asymmetry: west-coast deserts
//! (Atacama, Namib) and humid east coasts are still missing.

use pf_contracts::laws::{ConstantSpecV1, ConstantV1, SourceClassV1, UniverseError, UnitV1};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

/// Law id of the option (recorded inside the constants set).
pub const CRH_RAINOUT_LAW_V1: &str = "crh-rainout-v1";

/// The approximation a run with the option records (appended to its tier's
/// list; runs without the option record nothing new).
pub const CRH_RAINOUT_APPROXIMATION_V1: super::ApproximationV0 = super::ApproximationV0 {
    id: "crh-rainout-v1",
    statement: "Background rainout f(lat) W / tau_r is multiplied by the logistic efficiency 1 / (1 + exp(-a (r - r_0))), r = W / W_sat at the cell's surface temperature: its dry tail is the observed exponential collapse of precipitation with falling column relative humidity (Bretherton, Peters and Back 2004, tropical-ocean monthly fit slope and, through r_0, its absolute rate), and it tends to the v1 rate in moist columns; applied at all latitudes and over land; saturation-excess fallout unchanged; Newton-linearized in the line solve.",
    constants: &["crh_rainout_exponent", "crh_reference", "rainout_residence_s", "circulation_rainout_modulation", "column_water_saturation_ref_kg_m2", "column_water_reference_temperature_c", "clausius_clapeyron_rate_per_k"],
};

/// Usage limit a run with the option records (appended).
pub const CRH_RAINOUT_USAGE_LIMIT_V1: &str = "crh-rainout-v1: tropical-ocean precipitation/humidity fit applied everywhere; no subtropical-high zonal asymmetry or monsoon, so west-coast deserts and humid east coasts are missing";

/// Constants of `crh-rainout-v1`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CrhRainoutV1 {
    pub law: Cow<'static, str>,
    /// a in G(r) = 1 / (1 + exp(-a (r - r_0))).
    pub crh_rainout_exponent: ConstantV1,
    /// r_0: column relative humidity at which G = 1/2.
    pub crh_reference: ConstantV1,
}

/// Code-owned envelopes (inclusive).
pub const CRH_RAINOUT_CONSTANT_SPECS_V1: &[ConstantSpecV1] = &[
    ConstantSpecV1 { id: "crh_rainout_exponent", unit: UnitV1::Dimensionless, lo: 5.0, hi: 20.0, why: "observed fits of ln P against column relative humidity have slopes 11.4 (monthly) to 15.6 (daily) (Bretherton et al. 2004); longer averaging lowers the slope, so the envelope extends below the monthly value" },
    ConstantSpecV1 { id: "crh_reference", unit: UnitV1::Dimensionless, lo: 0.3, hi: 1.0, why: "a column relative humidity between subsidence-dry (~0.3) and saturated; observed tropical r spans about 0.25-0.85" },
];

fn awaiting(basis: &'static str) -> SourceClassV1 {
    SourceClassV1::ProjectAssumptionAwaitingReview { basis: Cow::Borrowed(basis) }
}

impl CrhRainoutV1 {
    /// Earth values (project assumptions awaiting review, ADR 039 C5).
    pub fn earth_v1() -> Self {
        Self {
            law: Cow::Borrowed(CRH_RAINOUT_LAW_V1),
            crh_rainout_exponent: ConstantV1::new(
                11.4,
                UnitV1::Dimensionless,
                awaiting("Bretherton, Peters and Back (2004, J. Climate 17:1517-1528; paper read 2026-10-21): monthly-mean fit P = exp(11.4 (r - 0.522)) mm/day over all tropical oceans, r = SSM/I water vapour path / ECMWF saturation path (daily fit a = 15.6, r = 0.603); the slope of the longest published averaging is used for annual means; applied at all latitudes and over land (project assumption)"),
            ),
            crh_reference: ConstantV1::new(
                0.55,
                UnitV1::Dimensionless,
                awaiting("project assumption (lane CLIMATE-ARIDITY, 2026-10-21): r at which the efficiency is 1/2. Physical bracket: matching the dry tail (f W_sat r / tau_r) exp(a (r - r_0)) to the Bretherton et al. (2004) monthly fit at r 0.45 gives r_0 0.56 (T 15 degC) to 0.63 (T 25 degC). Measured (UNEP AI classes, % of land; Earth UNEP 1992 arid + hyper-arid 19.6, drylands 47.2): r_0 0.55 gives C6 post hoc 18.3 / 37.3 and genesis (3 seeds, 120-400 My means) 21.3 / 42.3; r_0 0.6 gives 20.5 / 39.1 and 25.7 / 45.5; 0.55 chosen as the value inside both targets; C6 global P 887 mm/yr at 0.55 (Earth ~1000, recalled)"),
            ),
        }
    }

    pub fn validate(&self) -> Result<(), UniverseError> {
        if self.law != CRH_RAINOUT_LAW_V1 {
            return Err(UniverseError::Schema { what: "crh rainout law", found: self.law.to_string(), expected: CRH_RAINOUT_LAW_V1 });
        }
        for (id, c) in [("crh_rainout_exponent", &self.crh_rainout_exponent), ("crh_reference", &self.crh_reference)] {
            let spec = CRH_RAINOUT_CONSTANT_SPECS_V1.iter().find(|s| s.id == id).unwrap_or_else(|| panic!("no crh spec for {id}"));
            if c.unit != spec.unit {
                return Err(UniverseError::UnitMismatch { constant: spec.id, expected: spec.unit.as_str(), found: c.unit.as_str() });
            }
            if !c.value.is_finite() || !(spec.lo..=spec.hi).contains(&c.value) {
                return Err(UniverseError::OutOfEnvelope { constant: spec.id, value: c.value, unit: spec.unit.as_str(), lo: spec.lo, hi: spec.hi, why: spec.why });
            }
        }
        Ok(())
    }

    /// Pointwise parameters (a, r_0).
    pub fn params(&self) -> CrhParams {
        CrhParams { a: self.crh_rainout_exponent.value, r0: self.crh_reference.value }
    }
}

/// Resolved (a, r_0).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrhParams {
    pub a: f64,
    pub r0: f64,
}

impl CrhParams {
    /// Background rainout rate s(w) = k w G(w / w_sat) with the logistic
    /// efficiency G(r) = 1 / (1 + exp(-a (r - r_0))), and its derivative
    /// ds/dw = k G (1 + a r (1 - G)). `k` = f(lat) / tau_r (1/s), `w_sat` > 0.
    pub fn rate_and_slope(&self, k: f64, w: f64, w_sat: f64) -> (f64, f64) {
        let r = w / w_sat;
        let g = 1.0 / (1.0 + (-self.a * (r - self.r0)).exp());
        (k * w * g, k * g * (1.0 + self.a * r * (1.0 - g)))
    }

    /// The rate alone (pointwise evaluation of the same law).
    pub fn rate(&self, k: f64, w: f64, w_sat: f64) -> f64 {
        self.rate_and_slope(k, w, w_sat).0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crh_rate_known_answers() {
        let c = CrhRainoutV1::earth_v1();
        let p = c.params();
        let k = 1.0 / 734_400.0;
        let (a, r0) = (c.crh_rainout_exponent.value, c.crh_reference.value);
        // r = r_0: half the v1 rate.
        let w = r0 * 30.0;
        assert!((p.rate(k, w, 30.0) - 0.5 * k * w).abs() <= 1e-15 * k * w);
        // Dry tail: r_0 - 0.3 -> 1 / (1 + exp(0.3 a)) of the v1 rate
        // (a = 11.4: about 1/31).
        let w = (r0 - 0.3) * 30.0;
        let want = k * w / (1.0 + (0.3 * a).exp());
        assert!((p.rate(k, w, 30.0) - want).abs() <= 1e-12 * want);
        // Moist column: tends to the v1 rate (r = 1: 1 / (1 + exp(-0.4 a)),
        // 0.990 of it).
        let want = k * 30.0 / (1.0 + (-(1.0 - r0) * a).exp());
        assert!((p.rate(k, 30.0, 30.0) - want).abs() <= 1e-12 * want && want > 0.98 * k * 30.0);
    }

    #[test]
    fn crh_slope_is_the_derivative_and_newton_constant_is_non_negative() {
        let p = CrhRainoutV1::earth_v1().params();
        let (k, ws) = (1.3e-6, 40.0);
        for &w in &[5.0, 12.0, 20.0, 23.0, 30.0, 39.0] {
            let h = 1e-5;
            let (s1, _) = p.rate_and_slope(k, w + h, ws);
            let (s0, _) = p.rate_and_slope(k, w - h, ws);
            let (s, d) = p.rate_and_slope(k, w, ws);
            assert!(((s1 - s0) / (2.0 * h) - d).abs() <= 1e-6 * d.abs(), "w {w}");
            // Newton constant part s'(w) w - s(w) >= 0 (keeps the line
            // solve an M-matrix with a non-negative right-hand side).
            assert!(d * w - s >= -1e-18, "w {w}");
        }
    }

    #[test]
    fn crh_set_has_its_own_schema_and_digest_and_v1_is_unchanged() {
        use crate::climate::seed::{SeedClimateConstantsV1, SEED_CLIMATE_CONSTANTS_SCHEMA_CRH1, SEED_CLIMATE_CONSTANTS_SCHEMA_V1};
        let v1 = SeedClimateConstantsV1::earth_v1();
        let crh = SeedClimateConstantsV1::earth_crh_v1();
        assert!(v1.validate().is_ok() && crh.validate().is_ok());
        assert_eq!(v1.schema, SEED_CLIMATE_CONSTANTS_SCHEMA_V1);
        assert_eq!(crh.schema, SEED_CLIMATE_CONSTANTS_SCHEMA_CRH1);
        assert_ne!(v1.digest().unwrap(), crh.digest().unwrap());
        // The v1 set serializes without the option key (its bytes and digest
        // are the pre-option ones).
        let j = serde_json::to_string(&v1).unwrap();
        assert!(!j.contains("crh_rainout"), "{j}");
        // The option under the v1 schema, or the crh schema without the
        // option, is refused.
        let mut bad = v1.clone();
        bad.crh_rainout = Some(CrhRainoutV1::earth_v1());
        assert!(matches!(bad.validate(), Err(UniverseError::Schema { .. })));
        let mut bad = crh.clone();
        bad.crh_rainout = None;
        assert!(matches!(bad.validate(), Err(UniverseError::Schema { .. })));
        // Round trip.
        let back: SeedClimateConstantsV1 = serde_json::from_str(&serde_json::to_string(&crh).unwrap()).unwrap();
        assert_eq!(back, crh);
    }

    /// Synthetic fixture (instrument calibration, ADR 026 §1): a
    /// 100-degree-wide low continent from 50 S to 50 N. With the option the
    /// subtropical interior is drier than under v1 (subsidence suppression),
    /// the water ledger closes, and the deep tropics stay wet.
    #[test]
    fn crh_dries_the_subtropical_interior_and_closes_water() {
        use crate::climate::seed::{derive_seed_climate_v1, seed_climate_v1_earth_recipe, SeedClimateConstantsV1};
        use crate::climate::ClimateInputs;
        let g = pf_grid::GridSpec::new(36, 72, 6.371e6).unwrap();
        let (nlat, nlon) = (g.nlat(), g.nlon());
        let mut z = vec![-4000.0; g.cells()];
        for i in 0..nlat {
            let lat = g.center_lat_deg(i).unwrap();
            for j in 0..nlon {
                let lon = g.center_lon_deg(j).unwrap();
                if lat.abs() < 50.0 && (-50.0..50.0).contains(&lon) {
                    z[i * nlon + j] = 300.0;
                }
            }
        }
        let ocean: Vec<bool> = z.iter().map(|v| *v < 0.0).collect();
        let inputs = ClimateInputs { grid: &g, elevation_m: &z, fill_level_m: &z, sea_level_m: 0.0, ocean: &ocean };
        let r = seed_climate_v1_earth_recipe();
        let a = derive_seed_climate_v1(&inputs, &r, &SeedClimateConstantsV1::earth_v1()).unwrap();
        let b = derive_seed_climate_v1(&inputs, &r, &SeedClimateConstantsV1::earth_crh_v1()).unwrap();
        for c in &b.closures {
            assert!(c.closes(), "{c:?}");
        }
        let mean = |o: &crate::climate::ClimateOutputV0, lo: f64, hi: f64| -> f64 {
            let (mut s, mut n) = (0.0, 0.0);
            for i in 0..nlat {
                let lat = g.center_lat_deg(i).unwrap().abs();
                for j in 0..nlon {
                    let lon = g.center_lon_deg(j).unwrap();
                    if lat >= lo && lat < hi && (-20.0..20.0).contains(&lon) {
                        s += o.fields.precipitation_mm_per_yr[i * nlon + j];
                        n += 1.0;
                    }
                }
            }
            s / n
        };
        let (sub_v1, sub_crh) = (mean(&a, 20.0, 35.0), mean(&b, 20.0, 35.0));
        let (eq_v1, eq_crh) = (mean(&a, 0.0, 10.0), mean(&b, 0.0, 10.0));
        assert!(sub_crh < 0.5 * sub_v1, "subtropical interior: crh {sub_crh} vs v1 {sub_v1}");
        assert!(eq_crh > 0.7 * eq_v1, "deep tropics: crh {eq_crh} vs v1 {eq_v1}");
        assert!(b.summary.land_aridity_classes.hyper_arid + b.summary.land_aridity_classes.arid > a.summary.land_aridity_classes.hyper_arid + a.summary.land_aridity_classes.arid);
    }

    #[test]
    fn crh_validate_refuses_out_of_envelope() {
        let mut c = CrhRainoutV1::earth_v1();
        assert!(c.validate().is_ok());
        c.crh_rainout_exponent.value = 25.0;
        assert!(matches!(c.validate(), Err(UniverseError::OutOfEnvelope { .. })));
        let mut c = CrhRainoutV1::earth_v1();
        c.crh_reference.value = 0.1;
        assert!(matches!(c.validate(), Err(UniverseError::OutOfEnvelope { .. })));
        let mut c = CrhRainoutV1::earth_v1();
        c.law = Cow::Borrowed("other");
        assert!(matches!(c.validate(), Err(UniverseError::Schema { .. })));
    }
}
