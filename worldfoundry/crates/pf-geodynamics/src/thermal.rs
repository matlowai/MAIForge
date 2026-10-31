//! Layered thermal state (brief M0/M1): layer geometry, analytic cooling
//! references, and the explicit per-step anomaly update with ridge reset,
//! trench cold injection and plume warmth. No driver integration (Stage A
//! first); kernels only.
//!
//! UNITS: depth m, time My, temperature anomaly K, diffusivity m^2/My
//! (declared default ~1e-6 m^2/s = 31.5576 m^2/My; documented, overridable).
//!
//! What EVOLVES vs what is DIAGNOSTIC (M1 requirement):
//! - EVOLVES: `ThermalState::anomaly_k` (per element per layer), by the
//!   explicit relaxation + source-term update in `step`.
//! - DIAGNOSTIC (pure functions of inputs, never stored here): the
//!   half-space / plate cooling reference profiles, and the surface
//!   heat-flow proxy consumers derive from state + crust thickness (M3,
//!   computed downstream).

use pf_contracts::params::{Compensation, PlanetParameterProfile, ThermalModel};

#[derive(Debug, Clone, PartialEq)]
pub enum ThermalError {
    BadGeometry {
        why: String,
    },
    NonFinite {
        what: &'static str,
        value: f64,
        index: usize,
    },
    SizeMismatch {
        what: &'static str,
        expected: usize,
        found: usize,
    },
    BadMask {
        index: usize,
        cells: usize,
    },
    BadLayer {
        layer: usize,
        layers: usize,
    },
    NegativeDt {
        dt_my: f64,
    },
}

impl std::fmt::Display for ThermalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ThermalError::BadGeometry { why } => write!(f, "layer geometry: {why}"),
            ThermalError::NonFinite { what, value, index } => {
                write!(f, "nonfinite {what} {value:?} at flat index {index}")
            }
            ThermalError::SizeMismatch {
                what,
                expected,
                found,
            } => {
                write!(f, "{what}: expected {expected} entries, found {found}")
            }
            ThermalError::BadMask { index, cells } => {
                write!(f, "mask index {index} outside 0..{cells}")
            }
            ThermalError::BadLayer { layer, layers } => {
                write!(f, "layer {layer} outside 0..{layers}")
            }
            ThermalError::NegativeDt { dt_my } => write!(f, "dt_my {dt_my:?} must be > 0"),
        }
    }
}

impl std::error::Error for ThermalError {}

/// Default thermal diffusivity: 1e-6 m^2/s expressed in m^2/My
/// (1 My = 31,557,600 s, so MULTIPLY: 1e-6 * 3.15576e13 = 3.15576e7).
/// Sanity anchor: sqrt(kappa * 100 My) ~ 56 km (real lithosphere).
/// Declared, overridable per call.
pub const DEFAULT_KAPPA_M2_PER_MY: f64 = 1.0e-6 * 31_557_600_000_000.0;

/// Cold-anomaly injected down the layer column at trenches per step (K).
/// Coarse B1 source-term constant; declared here so tests and receipts can
/// pin it (a later cycle couples it to convergence rate).
pub const TRENCH_COLD_INJECTION_K_PER_STEP: f64 = 1.5;

/// M0 layer geometry: strictly increasing interface depths, surface first.
#[derive(Debug, Clone, PartialEq)]
pub struct LayerGeometry {
    /// n_layers + 1 interface depths [m]; interfaces[0] = 0 (surface).
    pub interfaces_m: Vec<f64>,
}

impl LayerGeometry {
    /// Default 4 layers (lithosphere, upper, transition, lower) over the
    /// declared mantle depth (M0 default 2.89e6 m to the core-mantle
    /// boundary-ish; parameter, not physics claim).
    pub fn default_four(mantle_depth_m: f64) -> Result<Self, ThermalError> {
        Self::uniform(4, mantle_depth_m)
    }

    /// n equal-thickness layers over the depth range.
    pub fn uniform(n_layers: usize, mantle_depth_m: f64) -> Result<Self, ThermalError> {
        if n_layers == 0 || n_layers > 64 {
            return Err(ThermalError::BadGeometry {
                why: format!("{n_layers} layers (1..=64)"),
            });
        }
        if !mantle_depth_m.is_finite() || mantle_depth_m <= 0.0 {
            return Err(ThermalError::BadGeometry {
                why: format!("depth {mantle_depth_m:?} m"),
            });
        }
        let mut interfaces_m = Vec::with_capacity(n_layers + 1);
        for k in 0..=n_layers {
            interfaces_m.push(mantle_depth_m * k as f64 / n_layers as f64);
        }
        Ok(Self { interfaces_m })
    }

    pub fn from_interfaces(interfaces_m: Vec<f64>) -> Result<Self, ThermalError> {
        if interfaces_m.len() < 2 {
            return Err(ThermalError::BadGeometry {
                why: "need at least 2 interfaces".into(),
            });
        }
        if interfaces_m[0] != 0.0 {
            return Err(ThermalError::BadGeometry {
                why: "first interface must be the surface (0 m)".into(),
            });
        }
        for w in interfaces_m.windows(2) {
            if w[1] <= w[0] || !w[1].is_finite() {
                return Err(ThermalError::BadGeometry {
                    why: format!("interfaces not strictly increasing: {w:?}"),
                });
            }
        }
        Ok(Self { interfaces_m })
    }

    pub fn n_layers(&self) -> usize {
        self.interfaces_m.len() - 1
    }

    pub fn layer_midpoint_m(&self, layer: usize) -> Result<f64, ThermalError> {
        if layer >= self.n_layers() {
            return Err(ThermalError::BadLayer {
                layer,
                layers: self.n_layers(),
            });
        }
        Ok(0.5 * (self.interfaces_m[layer] + self.interfaces_m[layer + 1]))
    }
}

/// Abramowitz & Stegun 7.1.26 erf approximation (|err| <= 1.5e-7) — no
/// external dependency, precision declared.
pub fn erf_approx(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * x);
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-x * x).exp();
    sign * y
}

/// Half-space cooling reference temperature (diagnostic, M1): T = Ts +
/// (Tm - Ts) erf(z / (2 sqrt(kappa t))).
pub fn half_space_temperature_c(
    depth_m: f64,
    age_my: f64,
    t_surface_c: f64,
    t_mantle_c: f64,
    kappa_m2_per_my: f64,
) -> Result<f64, ThermalError> {
    if !depth_m.is_finite() || depth_m < 0.0 || !age_my.is_finite() || age_my <= 0.0 {
        return Err(ThermalError::NonFinite {
            what: "half_space inputs",
            value: age_my,
            index: 0,
        });
    }
    let arg = depth_m / (2.0 * (kappa_m2_per_my * age_my).sqrt());
    Ok(t_surface_c + (t_mantle_c - t_surface_c) * erf_approx(arg))
}

/// Plate cooling reference temperature (diagnostic, M1) with a truncated
/// Fourier series; `terms` declared (error ~ exp(-pi^2 n^2 ...) tail).
pub fn plate_temperature_c(
    depth_m: f64,
    age_my: f64,
    plate_thickness_m: f64,
    t_surface_c: f64,
    t_mantle_c: f64,
    kappa_m2_per_my: f64,
    terms: usize,
) -> Result<f64, ThermalError> {
    if !depth_m.is_finite()
        || !(0.0..=plate_thickness_m).contains(&depth_m)
        || !age_my.is_finite()
        || age_my <= 0.0
        || plate_thickness_m <= 0.0
        || terms == 0
    {
        return Err(ThermalError::NonFinite {
            what: "plate inputs",
            value: depth_m,
            index: 0,
        });
    }
    let l = plate_thickness_m;
    let z = depth_m;
    let mut series = 0.0f64;
    for n in 1..=terms {
        let nf = n as f64;
        let sin = (std::f64::consts::PI * nf * z / l).sin();
        let decay =
            (-std::f64::consts::PI * std::f64::consts::PI * nf * nf * kappa_m2_per_my * age_my
                / (l * l))
                .exp();
        let before = series;
        series += sin / nf * decay;
        if series == before {
            break; // tail below f64 resolution — declared convergence stop
        }
    }
    Ok(t_surface_c + (t_mantle_c - t_surface_c) * (z / l + 2.0 / std::f64::consts::PI * series))
}

/// EVOLVING state: per-element per-layer temperature anomaly (K), flat
/// row-major [element * n_layers + layer]. What evolves: this vector, by
/// `step` only. Nonfinite entries are refused everywhere, loudly.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalState {
    pub anomaly_k: Vec<f64>,
    pub n_layers: usize,
    /// Per-layer relaxation timescale toward ambient (My); declared.
    pub relaxation_tau_my: Vec<f64>,
}

impl ThermalState {
    pub fn zeros(
        cells: usize,
        geometry: &LayerGeometry,
        tau_my: Vec<f64>,
    ) -> Result<Self, ThermalError> {
        let n = geometry.n_layers();
        if tau_my.len() != n {
            return Err(ThermalError::SizeMismatch {
                what: "relaxation_tau_my",
                expected: n,
                found: tau_my.len(),
            });
        }
        for &t in &tau_my {
            if !t.is_finite() || t <= 0.0 {
                return Err(ThermalError::NonFinite {
                    what: "tau_my",
                    value: t,
                    index: 0,
                });
            }
        }
        Ok(Self {
            anomaly_k: vec![0.0; cells * n],
            n_layers: n,
            relaxation_tau_my: tau_my,
        })
    }

    pub fn cells(&self) -> usize {
        self.anomaly_k.len() / self.n_layers
    }

    pub fn get(&self, element: usize, layer: usize) -> Result<f64, ThermalError> {
        self.index(element, layer).map(|i| self.anomaly_k[i])
    }

    fn index(&self, element: usize, layer: usize) -> Result<usize, ThermalError> {
        if layer >= self.n_layers {
            return Err(ThermalError::BadLayer {
                layer,
                layers: self.n_layers,
            });
        }
        let i = element * self.n_layers + layer;
        if i >= self.anomaly_k.len() {
            return Err(ThermalError::BadMask {
                index: element,
                cells: self.cells(),
            });
        }
        Ok(i)
    }

    /// One explicit thermal step (M1): aging relaxation toward ambient,
    /// ridge reset (anomaly -> 0 across the column), trench cold injection
    /// (down the column), plume warmth (across the column). Masks are
    /// sorted element indices; plumes are element indices.
    pub fn step(
        &mut self,
        dt_my: f64,
        ridge_mask: &[usize],
        trench_mask: &[usize],
        plume_elements: &[usize],
        plume_heat_k_per_step: f64,
    ) -> Result<(), ThermalError> {
        if !dt_my.is_finite() || dt_my <= 0.0 {
            return Err(ThermalError::NegativeDt { dt_my });
        }
        if !plume_heat_k_per_step.is_finite() || plume_heat_k_per_step < 0.0 {
            return Err(ThermalError::NonFinite {
                what: "plume_heat_k_per_step",
                value: plume_heat_k_per_step,
                index: 0,
            });
        }
        let cells = self.cells();
        for &e in ridge_mask
            .iter()
            .chain(trench_mask.iter())
            .chain(plume_elements.iter())
        {
            if e >= cells {
                return Err(ThermalError::BadMask { index: e, cells });
            }
        }
        // 1) aging relaxation per layer (element-parallel, layer-seriatim).
        for layer in 0..self.n_layers {
            let tau = self.relaxation_tau_my[layer];
            let decay = (-dt_my / tau).exp();
            for e in 0..cells {
                let i = e * self.n_layers + layer;
                self.anomaly_k[i] *= decay;
            }
        }
        // 2) ridge reset: spreading centers reset to ambient.
        for &e in ridge_mask {
            for layer in 0..self.n_layers {
                self.anomaly_k[e * self.n_layers + layer] = 0.0;
            }
        }
        // 3) trench cold injection down the column (declared constant).
        for &e in trench_mask {
            for layer in 0..self.n_layers {
                let i = e * self.n_layers + layer;
                self.anomaly_k[i] -= TRENCH_COLD_INJECTION_K_PER_STEP;
            }
        }
        // 4) plume warm injection across the column (per-step heat).
        for &e in plume_elements {
            for layer in 0..self.n_layers {
                let i = e * self.n_layers + layer;
                self.anomaly_k[i] += plume_heat_k_per_step;
            }
        }
        // 5) fail loud on any nonfinite result.
        for (i, &v) in self.anomaly_k.iter().enumerate() {
            if !v.is_finite() {
                return Err(ThermalError::NonFinite {
                    what: "anomaly_k after step",
                    value: v,
                    index: i,
                });
            }
        }
        Ok(())
    }
}

/// M2 compensation dispatch helper: Airy XOR Flexure, loudly (never both).
pub fn compensation_mode(params: &PlanetParameterProfile) -> Result<Compensation, ThermalError> {
    params.validate().map_err(|e| ThermalError::BadGeometry {
        why: format!("params: {e}"),
    })?;
    Ok(params.mantle.compensation)
}

/// Which analytic reference applies for a scenario (M1 diagnostic choice).
pub fn reference_model(params: &PlanetParameterProfile) -> ThermalModel {
    params.mantle.thermal_model
}

#[cfg(test)]
mod tests {
    use super::*;

    const K: f64 = DEFAULT_KAPPA_M2_PER_MY;
    const TS: f64 = 0.0;
    const TM: f64 = 1350.0;

    #[test]
    fn layer_geometry_validates_and_rejects() {
        let g = LayerGeometry::default_four(2.89e6).unwrap();
        assert_eq!(g.n_layers(), 4);
        assert_eq!(g.interfaces_m[0], 0.0);
        assert!((g.layer_midpoint_m(0).unwrap() - 2.89e6 / 8.0).abs() < 1.0);
        assert!(LayerGeometry::from_interfaces(vec![0.0, 1.0, 0.5]).is_err());
        assert!(LayerGeometry::from_interfaces(vec![1.0, 2.0]).is_err());
        assert!(LayerGeometry::uniform(0, 1.0).is_err());
        assert!(LayerGeometry::uniform(4, f64::NAN).is_err());
    }

    #[test]
    fn half_space_properties() {
        // Surface is at Ts at any age — within the A&S erf approximation's
        // declared 1.5e-7 error (|T error| <= 1350 * 1.5e-7 ~ 2e-4).
        assert!((half_space_temperature_c(0.0, 10.0, TS, TM, K).unwrap() - TS).abs() < 1e-3);
        // Diffusive-length sanity with the corrected kappa: sqrt(kappa*100
        // My) ~ 56 km, so 200 km depth at 100 My is ~1.8 diffusive lengths
        // below the surface -> >= 98% of the mantle temperature.
        let deep = half_space_temperature_c(2.0e5, 100.0, TS, TM, K).unwrap();
        assert!(
            deep > 0.98 * TM,
            "deep old lithosphere ~ mantle, got {deep}"
        );
        // Scaling invariance of the half-space solution: T(z, t) depends
        // only on z/sqrt(t): quadrupling both gives the same temperature.
        let a = half_space_temperature_c(50_000.0, 10.0, TS, TM, K).unwrap();
        let b = half_space_temperature_c(200_000.0, 160.0, TS, TM, K).unwrap();
        assert!((a - b).abs() < 1e-6, "erf scaling invariance: {a} vs {b}");
        // Monotone in depth.
        let t1 = half_space_temperature_c(30_000.0, 20.0, TS, TM, K).unwrap();
        let t2 = half_space_temperature_c(60_000.0, 20.0, TS, TM, K).unwrap();
        assert!(t2 > t1);
        assert!(half_space_temperature_c(-1.0, 10.0, TS, TM, K).is_err());
    }

    #[test]
    fn plate_model_steady_state_is_linear_and_ends_are_pinned() {
        let l = 125_000.0;
        // Steady state (huge age): series decayed, linear gradient. The
        // approach is asymptotic: at t=1e9 My the n=1 term still leaves a
        // ~2e-6 K residual — tolerance honest to that.
        let mid = plate_temperature_c(l / 2.0, 1.0e9, l, TS, TM, K, 64).unwrap();
        assert!(
            (mid - TM / 2.0).abs() < 1e-4,
            "steady mid-plate ~ (Ts+Tm)/2, got {mid}"
        );
        let top = plate_temperature_c(0.0, 1.0e9, l, TS, TM, K, 64).unwrap();
        let bot = plate_temperature_c(l, 1.0e9, l, TS, TM, K, 64).unwrap();
        assert!((top - TS).abs() < 1e-3 && (bot - TM).abs() < 1e-3);
        // Young plate at depth is HOTTER than the steady profile: the
        // plate model cools from a mantle-temperature initial state
        // (Fourier series -> T(t=0) = Tm everywhere), approaching the
        // linear steady gradient from above.
        let young = plate_temperature_c(l / 2.0, 1.0, l, TS, TM, K, 64).unwrap();
        assert!(
            young > mid,
            "cooling from mantle-initial: young {young} should exceed steady {mid}"
        );
        // Boundary inputs validated.
        assert!(plate_temperature_c(l + 1.0, 10.0, l, TS, TM, K, 8).is_err());
    }

    #[test]
    fn thermal_step_sources_and_relaxation() {
        let g = LayerGeometry::default_four(2.89e6).unwrap();
        let mut s = ThermalState::zeros(6, &g, vec![50.0; 4]).unwrap();
        // Plume warmth lands across the column, exactly.
        s.step(1.0, &[], &[], &[2], 12.0).unwrap();
        for layer in 0..4 {
            assert!((s.get(2, layer).unwrap() - 12.0).abs() < 1e-12);
        }
        // Trench cold injection after one more step: previous warmth decayed
        // one step then cold added.
        s.step(1.0, &[], &[5], &[], 0.0).unwrap();
        for layer in 0..4 {
            let v = s.get(5, layer).unwrap();
            assert!(
                (v + TRENCH_COLD_INJECTION_K_PER_STEP).abs() < 1e-12,
                "got {v}"
            );
        }
        // Ridge reset zeroes the column.
        s.step(1.0, &[2], &[], &[], 0.0).unwrap();
        for layer in 0..4 {
            assert_eq!(s.get(2, layer).unwrap(), 0.0);
        }
        // Relaxation decays exponentially: value after two steps of dt on
        // untouched element 1... set one explicitly then step twice.
        let mut s2 = ThermalState::zeros(2, &g, vec![10.0; 4]).unwrap();
        s2.anomaly_k[0] = 100.0;
        s2.step(10.0, &[], &[], &[], 0.0).unwrap();
        let expected = 100.0 * (-10.0f64 / 10.0).exp();
        assert!((s2.anomaly_k[0] - expected).abs() < 1e-9);
        // Loud: bad mask index, negative dt, nonfinite heat.
        assert!(s2.step(1.0, &[99], &[], &[], 0.0).is_err());
        assert!(s2.step(-1.0, &[], &[], &[], 0.0).is_err());
        assert!(s2.step(1.0, &[], &[], &[], f64::NAN).is_err());
    }

    #[test]
    fn step_is_deterministic() {
        let g = LayerGeometry::default_four(2.89e6).unwrap();
        let mut a = ThermalState::zeros(8, &g, vec![30.0; 4]).unwrap();
        let mut b = ThermalState::zeros(8, &g, vec![30.0; 4]).unwrap();
        for s in [&mut a, &mut b] {
            s.anomaly_k[3] = 5.0;
            s.step(2.0, &[1], &[4], &[6], 9.0).unwrap();
            s.step(3.0, &[], &[4], &[6], 9.0).unwrap();
        }
        assert_eq!(a.anomaly_k, b.anomaly_k);
    }

    #[test]
    fn compensation_mode_is_validated_and_xor() {
        use pf_contracts::params::{FidelityPreset, MantleBlock, ScenarioPreset};
        let mut p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
        assert_eq!(compensation_mode(&p).unwrap(), p.mantle.compensation);
        // Breaking the profile (invalid Te with flexure) is loud.
        p.mantle.effective_elastic_thickness_km = 0.0;
        assert!(compensation_mode(&p).is_err());
    }
}
