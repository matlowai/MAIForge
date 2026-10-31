//! GLACIAL-T0 (lane GLACIAL-T0, 2026-10-21): law id
//! `glacial-stream-power-blend-v1`, an opt-in glacial erosion law for T0
//! genesis on top of `stream-power-subgrid-nonlinear-v1` (`crate::relief`,
//! with or without `crate::erosion::incision_threshold`). Nothing here runs
//! unless the law is selected.
//!
//! PHYSICS (Hergarten 2021, Earth Surf. Dynam. 9:937-952, doi:10.5194/
//! esurf-9-937-2021, text read; benchmarked against the iSOSIA ice-flow model
//! by Liebl et al. 2023, Geosci. Model Dev. 16:1315-1343, doi:10.5194/
//! gmd-16-1315-2023, text read). In the sliding-dominated, thin-ice limit the
//! sliding law v_s ~ h^(psi-1) S^psi, Hallet-type erosion E ~ v_s^l and the
//! glacier width scaling w ~ q_i^0.3 (Bahr 1997; Prasicek et al. 2020) give a
//! stream-power law in the ice flux, E = K_g A_i^m S^n (Hergarten eq 19), with
//! the glacial concavity 0.47 ~ the fluvial one, so m and n equal the fluvial
//! exponents. A_i is the ice flux as a catchment-size equivalent (eq 17,
//! steady balance eq 21). Adding subglacial meltwater erosion with the fluvial
//! erodibility (eq 24-25) gives one law for the whole landscape (eq 26-28):
//!     E = K_eff A^m S^n,  K_eff = gamma^m K_g + (1 - gamma)^m K_f,
//!     gamma = A_i / A  (ice share of the total precipitation flux).
//! With K_f = K the profile never stands above the fluvial one and the
//! fluvial-to-glacial transition is smooth (Hergarten sec 4).
//!
//! T0 FORM. A T0 cell (~70 km) holds unresolved valleys, so the law acts as a
//! per-cell erodibility factor F = gamma^m r + (1 - gamma)^m (r = K_g / K_f) on
//! the two levels of the subgrid law:
//! - trunk (the resolved cell-to-cell channel): gamma from the routed ice and
//!   precipitation fluxes (Hergarten eq 20-21; ice flux clamped at 0, the
//!   glacier terminus);
//! - cell mean (the subgrid network): gamma at each height of the cell's
//!   hypsometry is the local accumulation share clamp((h - He)/(Hf - He), 0, 1)
//!   (declared: tongues below the ELA are not credited to the cell mean);
//!   F is averaged over the hypsometry.
//! The incision threshold (if selected) is a fluvial property: it is weighted
//! by the fluvial share (1 - gamma)^m. With no ice anywhere F = 1 and the
//! weight = 1 exactly, so the step is bit-identical to the law being off.
//!
//! ICE PRODUCTION (Hergarten eq 23; Liebl 2023 eq 10-11): above the ELA He,
//! p_i = P min((H - He)/(Hf - He), 1); below it the melt rate of Liebl eq 11,
//! p_i = -P (x + x^2/3), x = (He - H)/(Hf - He) (fits iSOSIA's PDD melt).
//! Cell hypsometry: uniform on [z - R, z + R] (z cell mean, R the carried
//! subgrid relief, so the trunk level is the cell's low end), midpoint
//! quadrature with `HYPSOMETRY_POINTS` points.
//!
//! ELA (Ohmura, Kasser and Funk 1992, J. Glaciol. 38(130):397-411, 70 glaciers
//! worldwide, standard error 200 mm; verified by lane GLACIAL-CHECK): at the
//! ELA the annual precipitation P (mm w.e.) and the summer (JJA north, DJF
//! south; the warmest 3 months) free-air temperature T satisfy
//!     P = 645 + 296 T + 9 T^2,  so  T_ela(P) = (-296 + sqrt(296^2 - 36 (645 - P))) / 18.
//! Genesis climate: free-air T(z) = t_sl - dT - lapse (z - datum), seasonal cycle
//! one harmonic of half range A, so the warmest-3-month mean is
//! T + A sin(pi/4)/(pi/4). Hence He = datum + (t_sl - dT + 0.9003 A - T_ela(f_P P)) / lapse.
//! P is the cell's precipitation from the genesis climate at the cell surface
//! (declared: no precipitation lapse inside a cell).
//!
//! CLIMATE STATES (time-averaged duty; the T0 step is 2 My, glacial cycles
//! ~0.1 My, so the law uses the duty-weighted mean rate). Each state is a
//! uniform temperature offset dT (no polar amplification, declared) and a
//! precipitation factor f_P on the genesis (modern-Earth-calibrated) climate:
//! - coldhouse with glacial cycles: dT varies sinusoidally in time between 0
//!   (interglacial = the genesis climate) and the glacial-maximum cooling, f_P
//!   between 1 and its glacial value (Liebl 2023 sec 3.5 used a sinusoidal ELA
//!   in time); `CYCLE_PHASES`-point midpoint quadrature in phase;
//! - warmer states (coolhouse, transitional, greenhouse) with dT < 0, f_P 1.
//! Since E is linear in K_eff at fixed topography, the duty-weighted mean
//! factor is exact for the rate at the step's topography.
//!
//! DETERMINISM: pure functions of the inputs; fixed loops, no threads, no
//! randomness. The diagnostics are a thread-local observation copy (never
//! run state, never in an output digest), as `incision_threshold`'s tally.

/// Law id (the typed law set lives in pf-contracts laws.rs).
pub const GLACIAL_EROSION_LAW_ID: &str = "glacial-stream-power-blend-v1";

/// Midpoints of the uniform cell hypsometry (numerics).
pub const HYPSOMETRY_POINTS: usize = 16;
/// Phase points of one glacial cycle (numerics).
pub const CYCLE_PHASES: usize = 8;
/// Mean of cos over the warmest quarter of a one-harmonic seasonal cycle:
/// sin(pi/4) / (pi/4).
pub const WARMEST_QUARTER_FACTOR: f64 = 0.900_316_316_157_106_1;

/// One warm (non-cycling) climate state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WarmStateParams {
    /// Fraction of time (0..1).
    pub duty: f64,
    /// Warming relative to the genesis climate (K, >= 0).
    pub warming_k: f64,
}

/// Resolved parameters of `glacial-stream-power-blend-v1`.
#[derive(Debug, Clone, PartialEq)]
pub struct GlacialParams {
    /// r = K_g / K_f.
    pub glacial_to_fluvial_erodibility_ratio: f64,
    /// Hf - He (m): height above the ELA where all precipitation is ice.
    pub full_ice_height_above_ela_m: f64,
    /// Fraction of time in the coldhouse with glacial cycles.
    pub coldhouse_duty: f64,
    /// Glacial-maximum cooling of the cycle (K, >= 0).
    pub glacial_max_cooling_k: f64,
    /// Precipitation factor at the glacial maximum (0..].
    pub glacial_max_precipitation_factor: f64,
    /// The warmer states (their duties plus the coldhouse duty sum to <= 1;
    /// any remainder is time without ice anywhere, i.e. contributes F = 1).
    pub warm_states: Vec<WarmStateParams>,
}

impl GlacialParams {
    pub fn validate(&self) -> Result<(), String> {
        let id = GLACIAL_EROSION_LAW_ID;
        let fin = |v: f64| v.is_finite();
        if !fin(self.glacial_to_fluvial_erodibility_ratio) || self.glacial_to_fluvial_erodibility_ratio <= 0.0 {
            return Err(format!("{id}: erodibility ratio {} must be finite and > 0", self.glacial_to_fluvial_erodibility_ratio));
        }
        if !fin(self.full_ice_height_above_ela_m) || self.full_ice_height_above_ela_m <= 0.0 {
            return Err(format!("{id}: full-ice height {} must be finite and > 0", self.full_ice_height_above_ela_m));
        }
        if !fin(self.glacial_max_cooling_k) || self.glacial_max_cooling_k < 0.0 {
            return Err(format!("{id}: glacial-maximum cooling {} must be finite and >= 0", self.glacial_max_cooling_k));
        }
        if !fin(self.glacial_max_precipitation_factor) || self.glacial_max_precipitation_factor <= 0.0 {
            return Err(format!("{id}: glacial-maximum precipitation factor {} must be finite and > 0", self.glacial_max_precipitation_factor));
        }
        if !fin(self.coldhouse_duty) || !(0.0..=1.0).contains(&self.coldhouse_duty) {
            return Err(format!("{id}: coldhouse duty {} outside [0, 1]", self.coldhouse_duty));
        }
        let mut total = self.coldhouse_duty;
        for w in &self.warm_states {
            if !fin(w.duty) || !(0.0..=1.0).contains(&w.duty) || !fin(w.warming_k) || w.warming_k < 0.0 {
                return Err(format!("{id}: warm state {w:?} needs duty in [0, 1] and warming >= 0"));
            }
            total += w.duty;
        }
        if total > 1.0 + 1e-12 {
            return Err(format!("{id}: state duties sum to {total} > 1"));
        }
        Ok(())
    }

    /// The climate states as (weight, cooling K, precipitation factor):
    /// the cycle phases first, then the warm states. Weights sum to the
    /// declared duties; the remainder (ice-free time) is implicit.
    pub fn states(&self) -> Vec<(f64, f64, f64)> {
        let mut s = Vec::with_capacity(CYCLE_PHASES + self.warm_states.len());
        if self.coldhouse_duty > 0.0 {
            let w = self.coldhouse_duty / CYCLE_PHASES as f64;
            for k in 0..CYCLE_PHASES {
                let phi = 2.0 * std::f64::consts::PI * (k as f64 + 0.5) / CYCLE_PHASES as f64;
                let c = 0.5 * (1.0 - phi.cos()); // 0 interglacial .. 1 glacial maximum
                s.push((w, c * self.glacial_max_cooling_k, 1.0 - c * (1.0 - self.glacial_max_precipitation_factor)));
            }
        }
        for ws in &self.warm_states {
            if ws.duty > 0.0 {
                s.push((ws.duty, -ws.warming_k, 1.0));
            }
        }
        s
    }

    /// Ice-free time fraction (contributes F = 1, weight 1).
    fn ice_free_duty(&self) -> f64 {
        (1.0 - self.coldhouse_duty - self.warm_states.iter().map(|w| w.duty).sum::<f64>()).max(0.0)
    }
}

/// Ohmura et al. 1992: summer free-air temperature (deg C) at the ELA for
/// annual precipitation `p_mm` (mm w.e. per yr, >= 0).
pub fn ohmura_ela_summer_temperature_c(p_mm: f64) -> f64 {
    let p = p_mm.max(0.0);
    (-296.0 + (296.0f64 * 296.0 - 36.0 * (645.0 - p)).sqrt()) / 18.0
}

/// ELA (m, same datum as the elevations) of one cell in one state.
pub fn ela_m(
    sea_level_temperature_c: f64,
    half_range_k: f64,
    precipitation_mm_per_yr: f64,
    cooling_k: f64,
    precipitation_factor: f64,
    lapse_k_per_m: f64,
    datum_m: f64,
) -> f64 {
    let t_ela = ohmura_ela_summer_temperature_c(precipitation_mm_per_yr * precipitation_factor);
    datum_m + (sea_level_temperature_c - cooling_k + WARMEST_QUARTER_FACTOR * half_range_k - t_ela) / lapse_k_per_m
}

/// Ice production share p_i / P at height h (Hergarten eq 23 above the ELA,
/// Liebl 2023 eq 11 melt below).
pub fn ice_share(h: f64, ela: f64, full: f64) -> f64 {
    let x = (h - ela) / full;
    if x >= 0.0 {
        x.min(1.0)
    } else {
        let y = -x;
        -(y + y * y / 3.0)
    }
}

/// Per-cell inputs of one step.
pub struct GlacialInputs<'a> {
    /// Receivers before donors (`Routing::receiver_first`).
    pub receiver_first: &'a [usize],
    pub flow_to: &'a [Option<usize>],
    /// Ocean cells (ice flux ends there).
    pub ocean: &'a [bool],
    /// Ocean or lake cells (hypsometry collapses to the water surface).
    pub wet: &'a [bool],
    pub fill_level_m: &'a [f64],
    pub area_m2: &'a [f64],
    pub z_mean_m: &'a [f64],
    /// Carried subgrid relief R (cell mean above the trunk), m.
    pub relief_m: &'a [f64],
    /// Genesis climate per cell: sea-level (zonal) temperature, seasonal
    /// half range, precipitation at the cell surface.
    pub sea_level_temperature_c: &'a [f64],
    pub half_range_k: &'a [f64],
    pub precipitation_mm_per_yr: &'a [f64],
    pub lapse_k_per_m: f64,
    pub datum_m: f64,
}

/// Erodibility factors of one step (per cell).
#[derive(Debug, Clone, PartialEq)]
pub struct GlacialErosivity {
    /// F on the trunk (resolved channel) erodibility.
    pub trunk_factor: Vec<f64>,
    /// Weight of the incision threshold on the trunk ((1 - gamma)^m mean).
    pub trunk_threshold_weight: Vec<f64>,
    /// F on the cell-mean (subgrid network) closure.
    pub mean_factor: Vec<f64>,
    pub mean_threshold_weight: Vec<f64>,
}

/// DIAGNOSTIC ONLY (never run state, never in an output or digest).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GlacialDiagnostics {
    /// ELA of the interglacial (genesis climate) and glacial-maximum states (m).
    pub ela_interglacial_m: Vec<f64>,
    pub ela_glacial_max_m: Vec<f64>,
    /// Duty-weighted mean trunk ice share gamma.
    pub trunk_gamma: Vec<f64>,
    /// Duty-weighted mean share of the cell's hypsometry above the ELA
    /// (accumulation area fraction), and the same at the glacial maximum.
    pub accumulation_fraction: Vec<f64>,
    pub accumulation_fraction_glacial_max: Vec<f64>,
    pub trunk_factor: Vec<f64>,
    pub mean_factor: Vec<f64>,
    /// GLACIAL-T1 (lane GLACIAL-T1, 2026-10-21): the ELA of each coldhouse
    /// cycle phase (`CYCLE_PHASES` fields, phase k at cooling
    /// c_k x glacial-maximum cooling, c_k = (1 - cos(2 pi (k + 0.5) / 8)) / 2;
    /// empty when the coldhouse duty is 0) and their mean (the cycle-mean
    /// ELA, empty likewise), plus the climate inputs of the step, so a
    /// published `t0.glacial` layer can be re-evaluated by T1b without
    /// re-deriving the climate. Observation only, like every field here.
    pub ela_cycle_phase_m: Vec<Vec<f64>>,
    pub ela_cycle_mean_m: Vec<f64>,
    pub sea_level_temperature_c: Vec<f64>,
    pub half_range_k: Vec<f64>,
    pub precipitation_mm_per_yr: Vec<f64>,
    pub lapse_k_per_m: f64,
    pub datum_m: f64,
}

thread_local! {
    static LAST: std::cell::RefCell<Option<GlacialDiagnostics>> = const { std::cell::RefCell::new(None) };
}

/// Diagnostic: the latest step's fields on this thread (taken, not copied).
pub fn take_last_diagnostics() -> Option<GlacialDiagnostics> {
    LAST.with(|t| t.borrow_mut().take())
}

fn hypsometry(z: f64, r: f64, k: usize) -> f64 {
    let r = r.max(0.0);
    z - r + 2.0 * r * (k as f64 + 0.5) / HYPSOMETRY_POINTS as f64
}

/// The per-cell factors of one step (and the thread-local diagnostics).
pub fn glacial_erosivity(p: &GlacialParams, m: f64, inp: &GlacialInputs<'_>) -> Result<GlacialErosivity, String> {
    p.validate()?;
    let n = inp.z_mean_m.len();
    for (what, len) in [
        ("flow_to", inp.flow_to.len()),
        ("ocean", inp.ocean.len()),
        ("wet", inp.wet.len()),
        ("fill", inp.fill_level_m.len()),
        ("area", inp.area_m2.len()),
        ("relief", inp.relief_m.len()),
        ("t_sl", inp.sea_level_temperature_c.len()),
        ("half_range", inp.half_range_k.len()),
        ("precipitation", inp.precipitation_mm_per_yr.len()),
    ] {
        if len != n {
            return Err(format!("{GLACIAL_EROSION_LAW_ID}: input {what}: {len} != {n}"));
        }
    }
    if !(m > 0.0 && m.is_finite()) || !(inp.lapse_k_per_m > 0.0 && inp.lapse_k_per_m.is_finite()) {
        return Err(format!("{GLACIAL_EROSION_LAW_ID}: m {m} and lapse {} must be positive", inp.lapse_k_per_m));
    }
    let r = p.glacial_to_fluvial_erodibility_ratio;
    let full = p.full_ice_height_above_ela_m;
    let states = p.states();
    let free = p.ice_free_duty();
    let blend = |g: f64| -> (f64, f64) {
        if g <= 0.0 {
            (1.0, 1.0)
        } else {
            let w = (1.0 - g).powf(m);
            (g.powf(m) * r + w, w)
        }
    };
    let mut tf = vec![0.0f64; n];
    let mut tw = vec![0.0f64; n];
    let mut mf = vec![0.0f64; n];
    let mut mw = vec![0.0f64; n];
    let mut diag = GlacialDiagnostics {
        ela_interglacial_m: vec![0.0; n],
        ela_glacial_max_m: vec![0.0; n],
        trunk_gamma: vec![0.0; n],
        accumulation_fraction: vec![0.0; n],
        accumulation_fraction_glacial_max: vec![0.0; n],
        trunk_factor: Vec::new(),
        mean_factor: Vec::new(),
        ela_cycle_phase_m: Vec::new(),
        ela_cycle_mean_m: Vec::new(),
        sea_level_temperature_c: inp.sea_level_temperature_c.to_vec(),
        half_range_k: inp.half_range_k.to_vec(),
        precipitation_mm_per_yr: inp.precipitation_mm_per_yr.to_vec(),
        lapse_k_per_m: inp.lapse_k_per_m,
        datum_m: inp.datum_m,
    };
    // GLACIAL-T1: per-phase and cycle-mean ELA of the coldhouse cycle (diagnostic).
    if p.coldhouse_duty > 0.0 {
        let mut mean = vec![0.0f64; n];
        for &(_, cool, fp) in states.iter().take(CYCLE_PHASES) {
            let e: Vec<f64> = (0..n)
                .map(|c| ela_m(inp.sea_level_temperature_c[c], inp.half_range_k[c], inp.precipitation_mm_per_yr[c], cool, fp, inp.lapse_k_per_m, inp.datum_m))
                .collect();
            for c in 0..n {
                mean[c] += e[c] / CYCLE_PHASES as f64;
            }
            diag.ela_cycle_phase_m.push(e);
        }
        diag.ela_cycle_mean_m = mean;
    }
    for c in 0..n {
        diag.ela_interglacial_m[c] = ela_m(inp.sea_level_temperature_c[c], inp.half_range_k[c], inp.precipitation_mm_per_yr[c], 0.0, 1.0, inp.lapse_k_per_m, inp.datum_m);
        diag.ela_glacial_max_m[c] = ela_m(
            inp.sea_level_temperature_c[c],
            inp.half_range_k[c],
            inp.precipitation_mm_per_yr[c],
            p.glacial_max_cooling_k,
            p.glacial_max_precipitation_factor,
            inp.lapse_k_per_m,
            inp.datum_m,
        );
    }
    // Glacial-maximum accumulation fraction (diagnostic).
    for c in 0..n {
        if inp.ocean[c] {
            continue;
        }
        let (z, rr) = if inp.wet[c] { (inp.fill_level_m[c], 0.0) } else { (inp.z_mean_m[c], inp.relief_m[c]) };
        let e = diag.ela_glacial_max_m[c];
        let above = (0..HYPSOMETRY_POINTS).filter(|&k| hypsometry(z, rr, k) > e).count();
        diag.accumulation_fraction_glacial_max[c] = above as f64 / HYPSOMETRY_POINTS as f64;
    }
    let mut ice = vec![0.0f64; n];
    let mut flux = vec![0.0f64; n];
    let mut prod = vec![0.0f64; n];
    for &(weight, cool, fp) in &states {
        // Per cell: ELA, net ice production share (hypsometry mean), the
        // cell-mean factor and threshold weight.
        for c in 0..n {
            if inp.ocean[c] {
                // No erosion there; the factors are 1 (no glacial change).
                prod[c] = 0.0;
                mf[c] += weight;
                mw[c] += weight;
                continue;
            }
            let e = ela_m(inp.sea_level_temperature_c[c], inp.half_range_k[c], inp.precipitation_mm_per_yr[c], cool, fp, inp.lapse_k_per_m, inp.datum_m);
            let (z, rr) = if inp.wet[c] { (inp.fill_level_m[c], 0.0) } else { (inp.z_mean_m[c], inp.relief_m[c]) };
            let (mut net, mut fsum, mut wsum, mut above) = (0.0, 0.0, 0.0, 0usize);
            for k in 0..HYPSOMETRY_POINTS {
                let h = hypsometry(z, rr, k);
                let s = ice_share(h, e, full);
                net += s;
                let (f, w) = blend(s.clamp(0.0, 1.0));
                fsum += f;
                wsum += w;
                if h > e {
                    above += 1;
                }
            }
            let q = HYPSOMETRY_POINTS as f64;
            prod[c] = net / q;
            mf[c] += weight * fsum / q;
            mw[c] += weight * wsum / q;
            diag.accumulation_fraction[c] += weight * above as f64 / q;
        }
        // Routed ice and precipitation fluxes (donors before receivers).
        for c in 0..n {
            let a = inp.area_m2[c] * inp.precipitation_mm_per_yr[c].max(0.0) * fp;
            flux[c] = a;
            ice[c] = a * prod[c];
        }
        for &c in inp.receiver_first.iter().rev() {
            if inp.ocean[c] {
                ice[c] = 0.0;
                continue;
            }
            ice[c] = ice[c].max(0.0);
            if let Some(j) = inp.flow_to[c] {
                flux[j] += flux[c];
                ice[j] += ice[c];
            }
        }
        for c in 0..n {
            let g = if flux[c] > 0.0 { (ice[c] / flux[c]).clamp(0.0, 1.0) } else { 0.0 };
            let (f, w) = blend(g);
            tf[c] += weight * f;
            tw[c] += weight * w;
            diag.trunk_gamma[c] += weight * g;
        }
    }
    // Ice-free time: F = 1, weight 1. The weights sum to 1 up to rounding,
    // so each field is divided by the summed weight; a field whose duty-
    // weighted sum equals that total to 1e-12 (no ice in any state) is set
    // to exactly 1, which keeps the zero-ice identity bit-exact.
    let duty: f64 = states.iter().map(|s| s.0).sum::<f64>() + free;
    for v in tf.iter_mut().chain(tw.iter_mut()).chain(mf.iter_mut()).chain(mw.iter_mut()) {
        *v += free;
        if (*v - duty).abs() <= 1e-12 * duty {
            *v = 1.0;
        } else {
            *v /= duty;
        }
    }
    diag.trunk_factor = tf.clone();
    diag.mean_factor = mf.clone();
    LAST.with(|t| *t.borrow_mut() = Some(diag));
    Ok(GlacialErosivity {
        trunk_factor: tf,
        trunk_threshold_weight: tw,
        mean_factor: mf,
        mean_threshold_weight: mw,
    })
}

#[cfg(test)]
mod tests {
    //! SYNTHETIC FIXTURES only (ADR 026 §1: instrument calibration).
    use super::*;

    fn params() -> GlacialParams {
        GlacialParams {
            glacial_to_fluvial_erodibility_ratio: 7.6,
            full_ice_height_above_ela_m: 500.0,
            coldhouse_duty: 0.13,
            glacial_max_cooling_k: 7.0,
            glacial_max_precipitation_factor: 0.7,
            warm_states: vec![
                WarmStateParams { duty: 0.18, warming_k: 6.0 },
                WarmStateParams { duty: 0.28, warming_k: 9.5 },
                WarmStateParams { duty: 0.41, warming_k: 14.0 },
            ],
        }
    }

    fn chain(z: [f64; 4], t_sl: f64) -> GlacialErosivity {
        let flow_to = [Some(1usize), Some(2), Some(3), None];
        let inp = GlacialInputs {
            receiver_first: &[3, 2, 1, 0],
            flow_to: &flow_to,
            ocean: &[false, false, false, true],
            wet: &[false, false, false, true],
            fill_level_m: &[z[0], z[1], z[2], 0.0],
            area_m2: &[4.9e9; 4],
            z_mean_m: &z,
            relief_m: &[800.0, 400.0, 100.0, 0.0],
            sea_level_temperature_c: &[t_sl; 4],
            half_range_k: &[10.0; 4],
            precipitation_mm_per_yr: &[1000.0; 4],
            lapse_k_per_m: 0.0065,
            datum_m: 0.0,
        };
        glacial_erosivity(&params(), 1.0, &inp).unwrap()
    }

    #[test]
    fn ohmura_inverts_the_regression() {
        for p in [0.0, 300.0, 1000.0, 2500.0, 4000.0] {
            let t = ohmura_ela_summer_temperature_c(p);
            let back = 645.0 + 296.0 * t + 9.0 * t * t;
            if p >= 0.0 {
                assert!((back - p).abs() < 1e-9 * p.max(1.0), "{p} {back}");
            }
        }
        // Wetter glaciers have warmer ELAs (lower ELA at the same climate).
        assert!(ohmura_ela_summer_temperature_c(3000.0) > ohmura_ela_summer_temperature_c(500.0));
    }

    #[test]
    fn hot_world_is_exactly_ice_free() {
        let g = chain([3000.0, 1500.0, 400.0, -4000.0], 60.0);
        for c in 0..4 {
            assert_eq!(g.trunk_factor[c], 1.0);
            assert_eq!(g.trunk_threshold_weight[c], 1.0);
            assert_eq!(g.mean_factor[c], 1.0);
            assert_eq!(g.mean_threshold_weight[c], 1.0);
        }
    }

    #[test]
    fn cold_high_cells_erode_faster_and_factor_is_bounded() {
        let g = chain([4500.0, 1500.0, 400.0, -4000.0], 15.0);
        let d = take_last_diagnostics().unwrap();
        // Interglacial ELA for t_sl 15, A 10, P 1000: T_ela = (-296 +
        // sqrt(296^2 + 36 x 355))/18 = 1.158 C -> (15 + 9.003 - 1.158)/0.0065
        // = 3514 m.
        assert!((d.ela_interglacial_m[0] - 3514.6).abs() < 1.0, "{}", d.ela_interglacial_m[0]);
        assert!(d.ela_glacial_max_m[0] < d.ela_interglacial_m[0] - 900.0);
        // GLACIAL-T1: 8 cycle phases between the two end states; the
        // symmetric sinusoid puts the cycle mean near the half-way ELA.
        assert_eq!(d.ela_cycle_phase_m.len(), CYCLE_PHASES);
        for e in &d.ela_cycle_phase_m {
            assert!(e[0] < d.ela_interglacial_m[0] && e[0] > d.ela_glacial_max_m[0]);
        }
        let mid = 0.5 * (d.ela_interglacial_m[0] + d.ela_glacial_max_m[0]);
        assert!((d.ela_cycle_mean_m[0] - mid).abs() < 60.0, "{} {}", d.ela_cycle_mean_m[0], mid);
        assert!(g.mean_factor[0] > 1.5, "{}", g.mean_factor[0]);
        assert!(g.mean_factor[0] <= 7.6 + 1e-12);
        assert!(g.mean_factor[2] == 1.0 || g.mean_factor[2] < g.mean_factor[0]);
        // Threshold weight drops where ice dominates.
        assert!(g.mean_threshold_weight[0] < 1.0);
        // Ice reaches the trunk of the highest cell only.
        assert!(d.trunk_gamma[0] > 0.0);
        assert!(g.trunk_factor[0] >= 1.0);
    }

    #[test]
    fn full_glaciation_gives_the_ratio() {
        let mut p = params();
        p.coldhouse_duty = 1.0;
        p.warm_states.clear();
        p.glacial_max_cooling_k = 0.0;
        let flow_to = [None];
        let inp = GlacialInputs {
            receiver_first: &[0],
            flow_to: &flow_to,
            ocean: &[false],
            wet: &[false],
            fill_level_m: &[5000.0],
            area_m2: &[1.0e9],
            z_mean_m: &[5000.0],
            relief_m: &[0.0],
            sea_level_temperature_c: &[-20.0],
            half_range_k: &[5.0],
            precipitation_mm_per_yr: &[500.0],
            lapse_k_per_m: 0.0065,
            datum_m: 0.0,
        };
        let g = glacial_erosivity(&p, 1.0, &inp).unwrap();
        assert!((g.mean_factor[0] - 7.6).abs() < 1e-12);
        assert!((g.trunk_factor[0] - 7.6).abs() < 1e-12);
        assert_eq!(g.mean_threshold_weight[0], 0.0);
    }
}
