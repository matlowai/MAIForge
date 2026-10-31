//! SEED-RIFT-SHOULDERS (lane RS, 2026-10-13): `shoulder-uplift-thermal-v1`,
//! the thermal and flexural response of continental crust stretched by
//! `stretched-rift-v1` (rift.rs).
//!
//! Without it a stretched column only floats lower (local Airy on the
//! thinner crust), so every rift is a step-sided trough that goes marine at
//! once and its flanks never rise. This law adds an ELEVATION term only (no
//! crust volume is created or moved), expressed as a water-loaded height
//! like the oceanic GDH1 term so the two-column isostasy converts it for
//! columns above sea level:
//!
//! (b) Thermal anomaly and thermal subsidence of the rift floor. A column
//!     stretched by a factor beta carries the uniform-stretching thermal
//!     uplift (1 - 1/beta) x C(age), with C the oceanic thermal elevation
//!     curve (GDH1: sqrt-age half-space branch early, exponential plate
//!     branch late) and age the time since stretching (McKenzie 1978
//!     class: the stretched lithosphere is hot in proportion to 1 - 1/beta
//!     and cools like oceanic plate). The floor therefore stands high while
//!     young and subsides thermally with time. Repeated stretching adds
//!     (f_new - f_old) C(0) to the current anomaly; the column keeps one
//!     equivalent age that reproduces the summed anomaly (continuous).
//! (a) Shoulders. Every continental column gets a shoulder from each
//!     stretched column within reach: `shoulder_ratio` x that column's
//!     thermal uplift x the cell-mean of the broken-plate flexural profile
//!     exp(-x/alpha) cos(x/alpha), with x measured from the source cell's
//!     edge and the profile cut at its first zero x = pi alpha / 2 (no
//!     negative lobe). The largest contribution counts (no summing along
//!     the rift), and it counts only where it exceeds the column's own
//!     thermal uplift. So the border-fault footwalls rise with the buoyant
//!     floor and the shoulders fade with the floor's thermal age.
//!
//! Declared limits: crust stays locally (Airy) compensated, so the floor's
//! crustal deficit is not flexurally supported and the flanks are not
//! down-warped by it; the shoulder is not volume-conserving buoyancy
//! (a rule, not a flexure solve); no plume swell; erosion acts on the
//! uplifted surface through the normal erosion law.
//!
//! Determinism: ascending cell order; maxima are order-independent.

use crate::kinematic::centroids;
use crate::rift::disc_with_distance;
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;

/// Parameters of `shoulder-uplift-thermal-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RiftFlankParams {
    /// Broken-plate flexural parameter alpha (m).
    pub flexural_parameter_m: f64,
    /// Shoulder at the rift edge / adjacent floor thermal uplift.
    pub shoulder_ratio: f64,
}

/// Age horizon (My) of the inverse thermal curve: beyond it the anomaly
/// is treated as fully decayed.
const AGE_HORIZON_MY: f64 = 1.0e4;

/// The equivalent age a with `curve(a) = y` for a non-increasing curve
/// (bisection; 0 when y is at or above curve(0), the horizon when y is at
/// or below curve(horizon)).
fn inverse_age(curve: &dyn Fn(f64) -> f64, y: f64) -> f64 {
    if y >= curve(0.0) {
        return 0.0;
    }
    if y <= curve(AGE_HORIZON_MY) {
        return AGE_HORIZON_MY;
    }
    let (mut lo, mut hi) = (0.0, AGE_HORIZON_MY);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if curve(mid) > y {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Thermal uplift (water-loaded m) of a column with stretching factor
/// `beta` and thermal start time `t0_my`, at `t_my`.
pub fn thermal_uplift_wl_m(beta: f64, t0_my: f64, t_my: f64, curve: &dyn Fn(f64) -> f64) -> f64 {
    if !(beta > 1.0) {
        return 0.0;
    }
    (1.0 - 1.0 / beta) * curve((t_my - t0_my).max(0.0))
}

/// Update the stretching markers after one step's rift pass. `beta` and
/// `t0_my` are the carried per-cell markers (1 and any value for an
/// unstretched column), `births` / `birth_flanks` the step's rift births
/// and their flank cells with weights, `thinning` the per-cell flank
/// thinning (m) and `thickness_m` the post-advance crust.
///
/// A birth's stretching factor is its flanks' weighted mean UNSTRETCHED
/// thickness (step-start thickness x step-start beta) over its own
/// thickness, so a birth between already-thinned flanks is not counted as
/// barely stretched. Births are set first, from the step-start markers;
/// then each flank's beta is multiplied by (h + dh) / h.
#[allow(clippy::too_many_arguments)]
pub fn stretch_update(
    beta: &mut [f64],
    t0_my: &mut [f64],
    births: &[u32],
    birth_flanks: &[Vec<(u32, f64)>],
    thinning: &[(u32, f64)],
    thickness_m: &[f64],
    t_now_my: f64,
    curve: &dyn Fn(f64) -> f64,
) -> Result<(), String> {
    let n = thickness_m.len();
    if beta.len() != n || t0_my.len() != n || births.len() != birth_flanks.len() {
        return Err("stretch_update: field sizes differ".into());
    }
    let c0 = curve(0.0);
    let mut step_dh = vec![0.0f64; 0];
    if !births.is_empty() {
        step_dh = vec![0.0; n];
        for &(c, dh) in thinning {
            step_dh[c as usize] += dh;
        }
    }
    let mut birth_beta = Vec::with_capacity(births.len());
    for (&b, flanks) in births.iter().zip(birth_flanks) {
        let b = b as usize;
        let h = thickness_m[b];
        let (mut num, mut den) = (0.0, 0.0);
        for &(c, q) in flanks {
            let c = c as usize;
            num += q * (thickness_m[c] + step_dh[c]) * beta[c].max(1.0);
            den += q;
        }
        if !(h > 0.0) || !(den > 0.0) || !num.is_finite() {
            return Err(format!("stretch_update: birth {b} thickness {h}, flank weight {den}"));
        }
        birth_beta.push((b, (num / den / h).max(1.0)));
    }
    for (b, bb) in birth_beta {
        beta[b] = bb;
        t0_my[b] = t_now_my;
    }
    for &(c, dh) in thinning {
        let c = c as usize;
        let h = thickness_m[c];
        if !(h > 0.0) || !(dh >= 0.0) {
            return Err(format!("stretch_update: flank {c} thickness {h}, thinning {dh}"));
        }
        let b_old = beta[c].max(1.0);
        let b_new = b_old * (h + dh) / h;
        let (f_old, f_new) = (1.0 - 1.0 / b_old, 1.0 - 1.0 / b_new);
        if !(f_new > 0.0) {
            continue;
        }
        let t_old = if f_old > 0.0 {
            f_old * curve((t_now_my - t0_my[c]).max(0.0))
        } else {
            0.0
        };
        let t_new = t_old + (f_new - f_old) * c0;
        beta[c] = b_new;
        t0_my[c] = t_now_my - inverse_age(curve, t_new / f_new);
    }
    Ok(())
}

/// The law's elevation term and its parts (water-loaded m per cell).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RiftFlankTerm {
    /// Thermal uplift of stretched continental columns.
    pub thermal_wl_m: Vec<f64>,
    /// Shoulder uplift beyond a column's own thermal uplift.
    pub shoulder_wl_m: Vec<f64>,
    /// thermal + shoulder: the column term handed to the isostasy.
    pub total_wl_m: Vec<f64>,
    pub stretched_cells: u64,
    pub shoulder_cells: u64,
}

/// Antiderivative of exp(-u) cos(u).
fn profile_integral(u: f64) -> f64 {
    0.5 * (-u).exp() * (u.sin() - u.cos())
}

/// Mean of the broken-plate profile p(u) = exp(-u) cos(u) for u in
/// [0, pi/2] (1 for u < 0, 0 beyond pi/2) over [u1, u2] (alpha units).
pub fn profile_cell_mean(u1: f64, u2: f64) -> f64 {
    let half_pi = std::f64::consts::FRAC_PI_2;
    if !(u2 > u1) {
        return 0.0;
    }
    let mut integral = 0.0;
    if u1 < 0.0 {
        integral += u2.min(0.0) - u1;
    }
    let (a, b) = (u1.max(0.0), u2.min(half_pi));
    if b > a {
        integral += profile_integral(b) - profile_integral(a);
    }
    integral / (u2 - u1)
}

/// Evaluate `shoulder-uplift-thermal-v1` on the grid at `t_my`.
#[allow(clippy::too_many_arguments)]
pub fn flank_term(
    grid: &GridSpec,
    areas: &[f64],
    kind: &[CrustKindCodeV0],
    beta: &[f64],
    t0_my: &[f64],
    t_my: f64,
    params: &RiftFlankParams,
    curve: &dyn Fn(f64) -> f64,
) -> Result<RiftFlankTerm, String> {
    let n = grid.cells();
    if [areas.len(), kind.len(), beta.len(), t0_my.len()]
        .iter()
        .any(|&l| l != n)
    {
        return Err("flank_term: field sizes differ from the grid".into());
    }
    let alpha = params.flexural_parameter_m;
    if !(alpha > 0.0) || !(params.shoulder_ratio >= 0.0) {
        return Err(format!("flank_term: bad parameters {params:?}"));
    }
    let mut thermal = vec![0.0; n];
    let mut stretched_cells = 0u64;
    for c in 0..n {
        if kind[c] == CrustKindCodeV0::Continental && beta[c] > 1.0 {
            thermal[c] = thermal_uplift_wl_m(beta[c], t0_my[c], t_my, curve);
            stretched_cells += 1;
        }
    }
    let mut shoulder_src = vec![0.0f64; n];
    if params.shoulder_ratio > 0.0 && stretched_cells > 0 {
        let cents = centroids(grid).map_err(|e| e.to_string())?;
        let max_half = areas.iter().fold(0.0f64, |m, a| m.max(a.sqrt() * 0.5));
        // SEED-PERF (lane SP): the sources' discs in parallel (pure), then
        // the max fold serially in source order (bit for bit).
        let srcs: Vec<usize> = (0..n).filter(|&r| thermal[r] > 0.0).collect();
        let per_src = crate::kinematic::par::map_heavy(srcs.len(), |k| -> Result<Vec<(usize, f64)>, String> {
            let r = srcs[k];
            let s_r = areas[r].sqrt() * 0.5;
            let reach = std::f64::consts::FRAC_PI_2 * alpha + s_r + max_half;
            let mut out = Vec::new();
            for (c, d) in disc_with_distance(grid, &cents, r, reach)? {
                if c == r || kind[c] != CrustKindCodeV0::Continental {
                    continue;
                }
                let s_c = areas[c].sqrt() * 0.5;
                let k = profile_cell_mean((d - s_r - s_c) / alpha, (d - s_r + s_c) / alpha);
                out.push((c, params.shoulder_ratio * thermal[r] * k));
            }
            Ok(out)
        });
        for cand in per_src {
            for (c, v) in cand? {
                if v > shoulder_src[c] {
                    shoulder_src[c] = v;
                }
            }
        }
    }
    let mut shoulder = vec![0.0; n];
    let mut total = vec![0.0; n];
    let mut shoulder_cells = 0u64;
    for c in 0..n {
        if shoulder_src[c] > thermal[c] {
            shoulder[c] = shoulder_src[c] - thermal[c];
            shoulder_cells += 1;
        }
        total[c] = thermal[c] + shoulder[c];
    }
    Ok(RiftFlankTerm {
        thermal_wl_m: thermal,
        shoulder_wl_m: shoulder,
        total_wl_m: total,
        stretched_cells,
        shoulder_cells,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// GDH1-shaped test curve (water-loaded thermal elevation, m).
    fn gdh1(age: f64) -> f64 {
        let a = age.max(0.0);
        let depth = if a < 20.0 {
            2600.0 + 365.0 * a.sqrt()
        } else {
            5651.0 - 2473.0 * (-0.0278 * a).exp()
        };
        (5651.0 - depth).max(0.0)
    }

    #[test]
    fn profile_mean_is_bounded_and_zero_beyond_first_zero() {
        assert!((profile_cell_mean(-1.0, -0.5) - 1.0).abs() < 1e-15);
        assert_eq!(profile_cell_mean(1.6, 3.0), 0.0);
        let m = profile_cell_mean(0.0, 1.0);
        // int_0^1 e^-u cos u du = 0.5 + e^-1 (sin 1 - cos 1) / 2
        let exact = 0.5 + 0.5 * (-1.0f64).exp() * (1.0f64.sin() - 1.0f64.cos());
        assert!((m - exact).abs() < 1e-12, "{m} {exact}");
        let mut prev = 1.0;
        for k in 0..20 {
            let u = k as f64 * 0.1;
            let v = profile_cell_mean(u, u + 0.1);
            assert!(v <= prev + 1e-12 && v >= 0.0);
            prev = v;
        }
    }

    #[test]
    fn repeated_stretching_is_continuous_and_matches_one_shot_when_simultaneous() {
        // One column thinned 30 -> 20 km in one step equals 30 -> 25 -> 20
        // within the same step time.
        let curve = |a: f64| gdh1(a);
        let mut b1 = vec![1.0];
        let mut t1 = vec![-1.0];
        stretch_update(&mut b1, &mut t1, &[], &[], &[(0, 10_000.0)], &[20_000.0], 10.0, &curve).unwrap();
        let mut b2 = vec![1.0];
        let mut t2 = vec![-1.0];
        stretch_update(&mut b2, &mut t2, &[], &[], &[(0, 5_000.0)], &[25_000.0], 10.0, &curve).unwrap();
        stretch_update(&mut b2, &mut t2, &[], &[], &[(0, 5_000.0)], &[20_000.0], 10.0, &curve).unwrap();
        assert!((b1[0] - 1.5).abs() < 1e-12 && (b2[0] - 1.5).abs() < 1e-12);
        assert!((t1[0] - 10.0).abs() < 1e-9 && (t2[0] - 10.0).abs() < 1e-6, "{t1:?} {t2:?}");
        // A later increment keeps the anomaly continuous: just after it the
        // uplift is the decayed old one plus the new fresh part.
        let before = thermal_uplift_wl_m(b1[0], t1[0], 40.0, &curve);
        let f_old = 1.0 - 1.0 / b1[0];
        stretch_update(&mut b1, &mut t1, &[], &[], &[(0, 2_000.0)], &[18_000.0], 40.0, &curve).unwrap();
        let f_new = 1.0 - 1.0 / b1[0];
        let after = thermal_uplift_wl_m(b1[0], t1[0], 40.0, &curve);
        assert!((after - (before + (f_new - f_old) * gdh1(0.0))).abs() < 1e-6, "{before} {after}");
        assert!(t1[0] > 10.0 && t1[0] < 40.0);
        // Thermal subsidence: the anomaly decays with time.
        let later = thermal_uplift_wl_m(b1[0], t1[0], 80.0, &curve);
        assert!(later < after);
    }

    #[test]
    fn a_birth_between_thinned_flanks_counts_their_earlier_stretching() {
        // Flanks 0 and 2 were 30 km, already stretched to 20 km (beta 1.5);
        // this step each loses 2 km into birth 1, which ends 10 km thick.
        let curve = |a: f64| gdh1(a);
        let mut beta = vec![1.5, 1.0, 1.5];
        let mut t0 = vec![0.0, -1.0, 0.0];
        let h = [18_000.0, 10_000.0, 18_000.0];
        let flanks = vec![vec![(0u32, 1.0), (2u32, 1.0)]];
        stretch_update(&mut beta, &mut t0, &[1], &flanks, &[(0, 2_000.0), (2, 2_000.0)], &h, 10.0, &curve)
            .unwrap();
        assert!((beta[1] - 3.0).abs() < 1e-12, "{}", beta[1]);
        assert_eq!(t0[1], 10.0);
        assert!((beta[0] - 1.5 * 20.0 / 18.0).abs() < 1e-12);
    }

    #[test]
    fn a_rift_line_gets_raised_shoulders_that_fade_outward() {
        // 1-degree grid; a north-south line of stretched continental cells.
        let grid = GridSpec::new(180, 360, 6_371_000.0).unwrap();
        let n = grid.cells();
        let areas: Vec<f64> = (0..n)
            .map(|c| grid.cell_area_m2(c / grid.nlon()).unwrap())
            .collect();
        let kind = vec![CrustKindCodeV0::Continental; n];
        let mut beta = vec![1.0; n];
        let t0 = vec![0.0; n];
        let j0 = 180;
        for i in 80..100 {
            beta[i * 360 + j0] = 2.0;
        }
        let p = RiftFlankParams {
            flexural_parameter_m: 60_000.0,
            shoulder_ratio: 1.0,
        };
        let curve = |a: f64| gdh1(a);
        let t = flank_term(&grid, &areas, &kind, &beta, &t0, 0.0, &p, &curve).unwrap();
        let floor = t.total_wl_m[90 * 360 + j0];
        assert!((floor - 0.5 * gdh1(0.0)).abs() < 1e-9);
        let s1 = t.total_wl_m[90 * 360 + j0 + 1];
        let s2 = t.total_wl_m[90 * 360 + j0 - 1];
        assert!(s1 > 0.0 && (s1 - s2).abs() < 1e-9 && s1 < floor, "{s1} {s2} {floor}");
        let mut prev = s1;
        for dj in 2..4 {
            let v = t.total_wl_m[90 * 360 + j0 + dj];
            assert!(v <= prev, "{dj}: {v} > {prev}");
            prev = v;
        }
        // Zero beyond reach.
        assert_eq!(t.total_wl_m[90 * 360 + j0 + 5], 0.0);
        // Oceanic cells get nothing.
        let mut k2 = kind.clone();
        k2[90 * 360 + j0 + 1] = CrustKindCodeV0::Oceanic;
        let t2 = flank_term(&grid, &areas, &k2, &beta, &t0, 0.0, &p, &curve).unwrap();
        assert_eq!(t2.total_wl_m[90 * 360 + j0 + 1], 0.0);
        // The shoulders fade with the floor's thermal age.
        let t3 = flank_term(&grid, &areas, &kind, &beta, &t0, 60.0, &p, &curve).unwrap();
        assert!(t3.total_wl_m[90 * 360 + j0 + 1] < s1);
    }
}
