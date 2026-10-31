//! SEED-TRENCH (lane TN, 2026-10-13): `slab-flexure-trench-v1`, the
//! flexural depression of a subducting oceanic plate (trench) and its
//! outer rise, from the slab's pull on a broken elastic plate.
//!
//! Without it the downgoing plate floats at its GDH1 depth right up to the
//! boundary, so the deepest ocean is the oldest abyssal plain and no
//! trench forms. This law adds an ELEVATION term only (no crust or
//! sediment is created or moved), as a water-loaded height like the
//! oceanic GDH1 term, for the two-column isostasy:
//!
//! 1. Polarity, the kinematic overlap rule: at a cell classified
//!    Subduction (convergent and not both continental), an OCEANIC cell is
//!    a trench cell when a neighbour of another plate overrides it: the
//!    neighbour is continental, or both are oceanic and the neighbour is
//!    younger (equal age: the lower plate id overrides). Continental crust
//!    never subducts here, so there is no trench at continent-continent,
//!    divergent or transform boundaries.
//! 2. Load: the slab pull per unit trench length of a half-space plate of
//!    age A, F = 2 rho_m alpha_v T_m g sqrt(kappa A / pi) L (Turcotte and
//!    Schubert class), with slab length L = min(L_max, v tau), tau =
//!    h_L^2 / (pi^2 kappa), h_L = 2.32 sqrt(kappa A): a slab stays cold
//!    for about v h^2 / kappa (McKenzie 1969 class), so young, slow slabs
//!    pull little. A fraction f of it is the vertical end load V0 = f F.
//! 3. Plate: Te = depth of the T_e isotherm in the half-space plate,
//!    2 erfinv(T_e / T_m) sqrt(kappa A) (Watts 2001 class); D = E Te^3 /
//!    (12 (1 - nu^2)); alpha = (4 D / ((rho_m - rho_w) g))^(1/4).
//! 4. Broken-plate deflection w(x) = w0 exp(-x/alpha) cos(x/alpha), w0 =
//!    V0 alpha^3 / (4 D) = V0 / ((rho_m - rho_w) g alpha), x from the
//!    trench axis (the edge of the trench cell toward the overriding
//!    plate), cut at the second zero 3 pi alpha / 2: trench at x = 0,
//!    outer rise of 0.067 w0 at 3 pi alpha / 4. Every oceanic cell of the
//!    downgoing plate takes the profile of its nearest trench cell (axis
//!    distance; ties to the lowest cell id), averaged over its own extent.
//!
//! Recomputed from the current state at every isostasy evaluation, so the
//! trench migrates with the boundary. Declared limits: no bending moment
//! (pure end load), no yielding, no accretionary wedge or forearc
//! response, continental cells of the downgoing plate are not flexed, a
//! cell-mean at grid resolution (a 69 km cell averages the trench V).
//!
//! Determinism: ascending cell order, strict improvements only.

use crate::kinematic::centroids;
use crate::rift::disc_with_distance;
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_grid::GridSpec;

/// Seconds per My (Julian year), for kappa in m^2/yr x ages in My.
const SECONDS_PER_MY: f64 = 3.15576e13;
/// h_L = 2.32 sqrt(kappa A): thermal plate thickness of the half-space
/// (the 0.9 T_m isotherm, Turcotte and Schubert class).
const PLATE_THICKNESS_COEF: f64 = 2.32;

/// Parameters of `slab-flexure-trench-v1` plus the universe densities and
/// gravity it reads.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrenchParams {
    pub vertical_load_fraction: f64,
    pub slab_length_max_m: f64,
    pub thermal_expansivity_per_k: f64,
    pub mantle_temperature_c: f64,
    pub elastic_isotherm_c: f64,
    pub thermal_diffusivity_m2_per_yr: f64,
    pub elastic_modulus_pa: f64,
    pub poisson_ratio: f64,
    pub mantle_density_kg_m3: f64,
    pub water_density_kg_m3: f64,
    pub gravity_m_s2: f64,
}

impl TrenchParams {
    pub fn validate(&self) -> Result<(), String> {
        let pos = [
            self.vertical_load_fraction,
            self.slab_length_max_m,
            self.thermal_expansivity_per_k,
            self.mantle_temperature_c,
            self.elastic_isotherm_c,
            self.thermal_diffusivity_m2_per_yr,
            self.elastic_modulus_pa,
            self.mantle_density_kg_m3,
            self.water_density_kg_m3,
            self.gravity_m_s2,
        ];
        if pos.iter().any(|v| !(v.is_finite() && *v > 0.0))
            || !(self.poisson_ratio >= 0.0 && self.poisson_ratio < 0.5)
            || !(self.elastic_isotherm_c < self.mantle_temperature_c)
            || !(self.water_density_kg_m3 < self.mantle_density_kg_m3)
        {
            return Err(format!("trench: bad parameters {self:?}"));
        }
        Ok(())
    }
}

/// erf (Abramowitz and Stegun 7.1.26, |error| < 1.5e-7).
fn erf(x: f64) -> f64 {
    let s = x.signum();
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * x);
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-x * x).exp();
    s * y
}

/// erfinv on [0, 1) by bisection of `erf`.
fn erfinv(y: f64) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 6.0f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if erf(mid) < y {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Trench-axis deflection w0 (m, positive down, water-loaded) and the
/// flexural parameter alpha (m) of a plate of age `age_my` subducting at
/// `rate_m_per_my`. (0, 0) for a non-positive age or rate.
pub fn trench_w0_alpha_m(age_my: f64, rate_m_per_my: f64, p: &TrenchParams) -> (f64, f64) {
    if !(age_my > 0.0) || !(rate_m_per_my > 0.0) {
        return (0.0, 0.0);
    }
    let kappa_s = p.thermal_diffusivity_m2_per_yr / 3.15576e7; // m^2/s
    let kt = kappa_s * age_my * SECONDS_PER_MY; // m^2
    let g = p.gravity_m_s2;
    let drho = p.mantle_density_kg_m3 - p.water_density_kg_m3;
    // tau = h_L^2 / (pi^2 kappa) = 2.32^2 A / pi^2 (kappa cancels), in My.
    let tau_my = PLATE_THICKNESS_COEF * PLATE_THICKNESS_COEF * age_my
        / (std::f64::consts::PI * std::f64::consts::PI);
    let slab_len = p.slab_length_max_m.min(rate_m_per_my * tau_my);
    let pull = 2.0
        * p.mantle_density_kg_m3
        * p.thermal_expansivity_per_k
        * p.mantle_temperature_c
        * g
        * (kt / std::f64::consts::PI).sqrt()
        * slab_len;
    let v0 = p.vertical_load_fraction * pull;
    let te = 2.0 * erfinv(p.elastic_isotherm_c / p.mantle_temperature_c) * kt.sqrt();
    let d = p.elastic_modulus_pa * te.powi(3) / (12.0 * (1.0 - p.poisson_ratio * p.poisson_ratio));
    let alpha = (4.0 * d / (drho * g)).powf(0.25);
    (v0 / (drho * g * alpha), alpha)
}

/// Antiderivative of exp(-u) cos(u).
fn profile_integral(u: f64) -> f64 {
    0.5 * (-u).exp() * (u.sin() - u.cos())
}

/// Mean of the broken-plate profile exp(-u) cos(u) (u in [0, 3 pi / 2],
/// 0 beyond) over [u1, u2] (alpha units); u below 0 is clamped to 0.
pub fn trench_profile_mean(u1: f64, u2: f64) -> f64 {
    let (u1, u2) = (u1.max(0.0), u2.max(0.0));
    if !(u2 > u1) {
        return if u1 <= 0.0 { 1.0 } else { 0.0 };
    }
    let end = 1.5 * std::f64::consts::PI;
    let (a, b) = (u1.min(end), u2.min(end));
    let integral = if b > a {
        profile_integral(b) - profile_integral(a)
    } else {
        0.0
    };
    integral / (u2 - u1)
}

/// The law's elevation term (water-loaded m per cell; negative = down)
/// and its diagnostics.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TrenchTerm {
    pub total_wl_m: Vec<f64>,
    /// Per trench cell: (cell, w0 m, alpha m, plate age My, rate m/My).
    pub trench_cells: Vec<(u32, f64, f64, f64, f64)>,
    pub flexed_cells: u64,
    pub max_w0_m: f64,
    pub max_bulge_wl_m: f64,
}

/// True when `nb` overrides `c` at a subduction contact (the kinematic
/// overlap rule): c oceanic, nb continental, or both oceanic and nb
/// younger (equal birth: the lower plate id overrides).
fn overrides(
    kind: &[CrustKindCodeV0],
    birth_t_my: &[f64],
    plate: &[pf_contracts::domain::PlateId],
    c: usize,
    nb: usize,
) -> bool {
    if kind[c] != CrustKindCodeV0::Oceanic {
        return false;
    }
    match kind[nb] {
        CrustKindCodeV0::Continental => true,
        CrustKindCodeV0::Oceanic => {
            birth_t_my[nb] > birth_t_my[c] || (birth_t_my[nb] == birth_t_my[c] && plate[nb] < plate[c])
        }
    }
}

/// Evaluate `slab-flexure-trench-v1` on the grid for a crust state.
#[allow(clippy::too_many_arguments)]
pub fn trench_term(
    grid: &GridSpec,
    areas: &[f64],
    state: &crate::kinematic::CrustState,
    params: &TrenchParams,
) -> Result<TrenchTerm, String> {
    params.validate()?;
    let n = grid.cells();
    if state.cells() != n || areas.len() != n {
        return Err("trench_term: field sizes differ from the grid".into());
    }
    let sub = BoundaryClassCodeV0::Subduction.code();
    // SEED-PERF (lane SP, 2026-10-13): the per-cell source test in parallel
    // (neighbour ids only; the distances were unused), collected in cell
    // order (bit for bit; first error in cell order).
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    let found = crate::kinematic::par::map_indexed(n, |c| -> Result<Option<(u32, f64, f64, f64, f64)>, String> {
        if state.boundary_class[c] != sub
            || state.kind[c] != CrustKindCodeV0::Oceanic
            || !(state.normal_rate_m_per_my[c] > 0.0)
        {
            return Ok(None);
        }
        grid.ij(c as u64).map_err(|e| e.to_string())?;
        let mut nbs = [0usize; 8];
        let nn = crate::kinematic::neighbor_ids(nlat, nlon, c, &mut nbs);
        let mut downgoing = false;
        for &nb in &nbs[..nn] {
            if state.plate[nb] != state.plate[c]
                && overrides(&state.kind, &state.birth_t_my, &state.plate, c, nb)
            {
                downgoing = true;
                break;
            }
        }
        if !downgoing {
            return Ok(None);
        }
        let age = state.age_my(c);
        let rate = state.normal_rate_m_per_my[c];
        let (w0, alpha) = trench_w0_alpha_m(age, rate, params);
        Ok((w0 > 0.0 && alpha > 0.0).then_some((c as u32, w0, alpha, age, rate)))
    });
    let mut sources: Vec<(u32, f64, f64, f64, f64)> = Vec::new();
    for f in found {
        if let Some(x) = f? {
            sources.push(x);
        }
    }
    let mut total = vec![0.0f64; n];
    let mut best_x = vec![f64::INFINITY; n];
    let mut flexed = 0u64;
    let mut max_w0 = 0.0f64;
    let mut max_bulge = 0.0f64;
    if !sources.is_empty() {
        let cents = centroids(grid).map_err(|e| e.to_string())?;
        let max_half = areas.iter().fold(0.0f64, |m, a| m.max(a.sqrt() * 0.5));
        // SEED-PERF (lane SP): the sources' discs in parallel (pure), then
        // the min-distance fold serially in source order (bit for bit).
        // Per source in parallel: its disc's eligible targets with their
        // axis distance and deflection (pure); then the nearest-axis fold
        // serially in source order (bit for bit).
        let per_src = crate::kinematic::par::map_heavy(sources.len(), |k| -> Result<Vec<(usize, f64, f64)>, String> {
            let (s, w0, alpha, _, _) = sources[k];
            let s = s as usize;
            let s_half = areas[s].sqrt() * 0.5;
            let reach = 1.5 * std::f64::consts::PI * alpha + s_half + max_half;
            let mut out = Vec::new();
            for (t, d) in disc_with_distance(grid, &cents, s, reach)? {
                if state.plate[t] != state.plate[s] || state.kind[t] != CrustKindCodeV0::Oceanic {
                    continue;
                }
                // Axis distance of the target's centre: the axis lies half a
                // cell behind the trench cell's centre.
                let x = d + s_half;
                let t_half = areas[t].sqrt() * 0.5;
                let kk = trench_profile_mean((x - t_half) / alpha, (x + t_half) / alpha);
                out.push((t, x, -w0 * kk));
            }
            Ok(out)
        });
        for (&(_, w0, _, _, _), cand) in sources.iter().zip(per_src) {
            max_w0 = max_w0.max(w0);
            for (t, x, v) in cand? {
                if !(x < best_x[t]) {
                    continue;
                }
                best_x[t] = x;
                total[t] = v;
            }
        }
        for c in 0..n {
            if best_x[c].is_finite() && total[c] != 0.0 {
                flexed += 1;
            }
            max_bulge = max_bulge.max(total[c]);
        }
    }
    Ok(TrenchTerm {
        total_wl_m: total,
        trench_cells: sources,
        flexed_cells: flexed,
        max_w0_m: max_w0,
        max_bulge_wl_m: max_bulge,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p() -> TrenchParams {
        TrenchParams {
            vertical_load_fraction: 0.1,
            slab_length_max_m: 660_000.0,
            thermal_expansivity_per_k: 3.0e-5,
            mantle_temperature_c: 1300.0,
            elastic_isotherm_c: 450.0,
            thermal_diffusivity_m2_per_yr: 31.5576,
            elastic_modulus_pa: 7.0e10,
            poisson_ratio: 0.25,
            mantle_density_kg_m3: 3300.0,
            water_density_kg_m3: 1030.0,
            gravity_m_s2: 9.81,
        }
    }

    #[test]
    fn erf_and_inverse_are_consistent() {
        assert!((erf(0.5) - 0.520_499_877_8).abs() < 2e-7);
        assert!((erf(erfinv(1.0 / 3.0)) - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn deflection_grows_with_age_and_rate_and_saturates() {
        let q = p();
        let (w_old, a_old) = trench_w0_alpha_m(100.0, 50_000.0, &q);
        // Pre-registered: about 2.8 km relief, alpha about 84 km.
        assert!((w_old - 2813.0).abs() < 30.0, "{w_old}");
        assert!((a_old - 84_000.0).abs() < 1_500.0, "{a_old}");
        let (w_young, _) = trench_w0_alpha_m(10.0, 30_000.0, &q);
        assert!(w_young < 0.25 * w_old, "{w_young}");
        let (w_slow, _) = trench_w0_alpha_m(20.0, 10_000.0, &q);
        let (w_fast, _) = trench_w0_alpha_m(20.0, 60_000.0, &q);
        assert!(w_slow < w_fast);
        // Saturated slab length: rate no longer matters.
        let (w1, _) = trench_w0_alpha_m(150.0, 30_000.0, &q);
        let (w2, _) = trench_w0_alpha_m(150.0, 100_000.0, &q);
        assert_eq!(w1, w2);
        assert_eq!(trench_w0_alpha_m(0.0, 1.0, &q), (0.0, 0.0));
        assert_eq!(trench_w0_alpha_m(50.0, -1.0, &q), (0.0, 0.0));
    }

    #[test]
    fn profile_mean_has_trench_and_outer_rise() {
        let half = std::f64::consts::FRAC_PI_2;
        // Point values through narrow intervals.
        let at = |u: f64| trench_profile_mean(u - 1e-6, u + 1e-6);
        assert!((at(0.5) - (-0.5f64).exp() * 0.5f64.cos()).abs() < 1e-6);
        assert!(at(half + 0.1) < 0.0);
        let bulge = at(0.75 * std::f64::consts::PI);
        assert!((bulge + 0.0670).abs() < 1e-3, "{bulge}");
        assert_eq!(trench_profile_mean(5.0, 6.0), 0.0);
    }

    #[test]
    fn trench_forms_on_the_oceanic_side_only() {
        use pf_contracts::domain::PlateId;
        // 1-degree grid: plate 0 continental west of column 180, plate 1
        // oceanic east of it, converging (class/rate set directly).
        let grid = GridSpec::new(180, 360, 6_371_000.0).unwrap();
        let n = grid.cells();
        let areas: Vec<f64> = (0..n).map(|c| grid.cell_area_m2(c / grid.nlon()).unwrap()).collect();
        let mut st = crate::kinematic::CrustState {
            t_my: 100.0,
            plate: vec![PlateId(0); n],
            kind: vec![CrustKindCodeV0::Continental; n],
            birth_t_my: vec![-1000.0; n],
            thickness_m: vec![35_000.0; n],
            sediment_m: vec![0.0; n],
            eroded_total_m: vec![0.0; n],
            arc_added_m: vec![0.0; n],
            shortening_added_m: vec![0.0; n],
            boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
            normal_rate_m_per_my: vec![0.0; n],
            material_pos: vec![[0.0, 0.0, 1.0]; n],
        };
        for i in 0..180 {
            for j in 180..360 {
                let c = i * 360 + j;
                st.plate[c] = PlateId(1);
                st.kind[c] = CrustKindCodeV0::Oceanic;
                st.birth_t_my[c] = 0.0;
                st.thickness_m[c] = 7_000.0;
            }
            for j in [179usize, 180] {
                st.boundary_class[i * 360 + j] = BoundaryClassCodeV0::Subduction.code();
                st.normal_rate_m_per_my[i * 360 + j] = 50_000.0;
            }
        }
        let t = trench_term(&grid, &areas, &st, &p()).unwrap();
        let row = 90 * 360;
        assert!(t.trench_cells.iter().all(|s| (s.0 as usize) % 360 == 180));
        assert!(!t.trench_cells.is_empty());
        assert_eq!(t.total_wl_m[row + 179], 0.0, "continental overriding side untouched");
        let trench = t.total_wl_m[row + 180];
        assert!(trench < -1000.0, "{trench}");
        // Shallows seaward, then an outer rise, then nothing.
        assert!(t.total_wl_m[row + 181] > trench);
        let bulge = (181..190).map(|j| t.total_wl_m[row + j]).fold(f64::MIN, f64::max);
        assert!(bulge > 50.0 && bulge < 0.1 * -trench, "{bulge} {trench}");
        assert_eq!(t.total_wl_m[row + 200], 0.0);
        // Continent-continent: no trench.
        let mut cc = st.clone();
        cc.kind = vec![CrustKindCodeV0::Continental; n];
        let t2 = trench_term(&grid, &areas, &cc, &p()).unwrap();
        assert!(t2.trench_cells.is_empty() && t2.total_wl_m.iter().all(|&v| v == 0.0));
        // Divergent (ridge class): no trench.
        let mut rd = st.clone();
        rd.boundary_class = vec![BoundaryClassCodeV0::Ridge.code(); n];
        assert!(trench_term(&grid, &areas, &rd, &p()).unwrap().trench_cells.is_empty());
    }
}
