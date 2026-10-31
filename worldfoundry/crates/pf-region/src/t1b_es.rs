//! Lane E-SOLVE-GENERALITY (2026-10-21): opt-in E-solve generality options for the chart-wide T1b rivernet E solve
//! (`t1b_rivernet::solve_e`). Key `es=1` (layer AND windows; windows only read the table, the key keeps the recipes equal).
//! Absent = the 0.9.x solve unchanged, byte for byte. Never seen by the keyed draws (draw_digest strips it).
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon; a construction solve (no SimTime).
//!
//! Diagnosis (GT1-L1-20261021 site 1 and the C6 site, VZ_SOLVE_TRACE, `unpublished-design-record`): the default damped
//! Jacobi assumes the cell relief above the trunk level responds as M - L ~ E (elasticity g = d ln(M - L) / d ln E = 1, step
//! 0.6 ln((T - L)/(M - L))). Measured g in the converging phase is 0.1 to 0.6 (bed slopes reach the cap x Sc, hillslopes
//! saturate at Sc f(beta) -> Sc, and the hat-blended E shares posts with the neighbours), so each update closes only
//! 0.6 g of the log error and 20 updates are not enough. Three further classes cannot meet the target by any E:
//! (1) chart-rim cells whose centre lies outside the chart lattice (no or few own posts; their in-chart hat mass is the tail
//! over neighbour cells, whose E controls those posts), (2) cells whose target lies below the anchor mean (T < L: sea-level
//! floored trunk levels on marine / coastal cells) and (3) cells pinned at an E bound whose residual points outward.
//!
//! Options (EsParamsV1, one switch each so they can be A/B'd):
//! - update 1 (secant_diagonal_v1): per-cell secant estimate of g from the cell's own last two evaluations, step = ln((T - L)
//!   / (M - L)) / g, g geometric-smoothed and clamped to [g_min, g_max]; the first update uses g = 1/0.6 (= the default
//!   Jacobi step); a residual sign change halves that one step (Illinois); fallback +-ln 2 by sign when T or M is within
//!   1 m of L (as the default) times a multiplier halved on every sign change there; step clamped to +-step_max.
//! - coverage_min (chart_rim_v1): a cell whose hat mass inside the chart lattice is below this fraction of its full hat mass
//!   (= T0 cell area / lattice spacing^2; 500 permille <=> the cell centre is outside the chart) is a rim cell: its F6 mean
//!   is not evaluable in this chart. rim_mode 2 (v1): solved with the default Jacobi step, outside the stop test (residuals
//!   reported); rim_mode 1: not solved, E = the solved 8-neighbours' mean (measured and rejected: rim strips +600..+1400 m).
//! - kkt 1 (bound_kkt_stop_v1): the loop also stops when every solved cell is within tolerance OR at an E bound with its
//!   residual pointing out of the box (E = e_min and M > T, or E = e_max and M < T), and every coastal cell meets F4 or has
//!   its flooding offset F pinned at a bound the same way: the bound-constrained fixed point.
//!   The strict `converged` flag is unchanged; the record lists every bound-limited cell with its residual.
//! - target 1 (plain_land_target_v1, v1 default): resolvable T0 Land cells solve the PLAIN whole-cell mean (ADR 043 F6
//!   literal) instead of the hat-weighted mean; the hat-mean solve left plain cell means off by up to 289 m (site 2: 5/13
//!   Land cells within F6; C6 4/6); with target 1: 13/13 and 6/6 (GT1 8/9: 44299 is bound-limited, T below the trunk mean).
//! The layer header gets `lane_es` with the per-cell record and the ADR 043 F6 literal (plain whole-cell mean over the
//! cell's own chart posts vs the bilinear mean, tolerance max(50 m, 10 % of T0 subgrid relief)) on T0 Land cells.

use serde::{Deserialize, Serialize};

pub const ES_RULE: &str = "e_solve_generality_v1";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct EsParamsV1 {
    /// 0 the default damped Jacobi, 1 secant_diagonal_v1.
    pub update: u32,
    pub g_min_permille: i64,
    pub g_max_permille: i64,
    /// Step clamp in ln E (permille of 1).
    pub step_max_permille: i64,
    /// chart_rim_v1 threshold (permille of the full hat mass; 0 off).
    pub coverage_min_permille: i64,
    /// chart_rim_v1 treatment: 1 not solved, neighbour-mean E (MEASURED AND REJECTED on GT1 / C6: the rim strips' plain
    /// means move +600 to +1400 m off the bilinear); 2 solved with the default Jacobi step (no secant) and left out of the
    /// stop test (their residuals are reported).
    pub rim_mode: u32,
    /// 1 bound_kkt_stop_v1.
    pub kkt: u32,
    /// 0 the default hat-weighted target; 1 plain_land_target_v1 (ADR 043 F6 literal on resolvable Land cells; v1).
    pub target: u32,
    pub update_rule: String,
    pub coverage_rule: String,
    pub target_rule: String,
    pub kkt_rule: String,
    /// Lane E-SOLVE-CONVERGE (opt-in esc=1): the convergent solve driver (t1b_esc.rs); None = unchanged (not serialized).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conv: Option<crate::t1b_esc::EscParamsV1>,
}

impl EsParamsV1 {
    pub fn v1() -> Self {
        EsParamsV1 {
            update: 1,
            g_min_permille: 100,
            g_max_permille: 2000,
            step_max_permille: 1386,
            coverage_min_permille: 500,
            rim_mode: 2,
            kkt: 1,
            target: 1,
            update_rule: "secant_diagonal_v1: x = ln E, y = ln(M - L) per solved cell; after an update with |dx| > 0.01 the observed g = dy/dx (own cell; the neighbours move too, so an effective value); g_new = sqrt(g_old x clamp(g_obs, g_min, g_max)) when g_obs > 0, else g_old; g starts at 1/0.6 (the default Jacobi step); step = ln((T - L)/(M - L)) / g, halved once when the residual sign changed since the previous evaluation; fallback +-ln 2 by the sign of M - T when T or M is within 1 m of L, times a multiplier halved on every sign change in the fallback (floor 1/32); step clamped to +-step_max; E clamped to [e_min, e_max]".into(),
            coverage_rule: "chart_rim_v1: coverage = (hat mass of the cell on the chart solve lattice) / (T0 cell area / spacing^2); 500 permille is the cell centre on the lattice edge; a cell the default rule solves with coverage < coverage_min is a rim cell: the chart holds less than half of its hat, so its ADR 043 F6 mean is not evaluable here and its posts are mostly other cells' (measured elasticity 0 to 0.1); rim_mode 2: rim cells are solved with the default damped Jacobi step (never the secant, which drove weakly controlled rim cells to E 1286 m/My and broke F4 on GT1) and do not enter the stop test; rim_mode 1: not solved, the mean E of the solved 8-neighbours (recomputed every update; E0 if none)".into(),
            target_rule: "plain_land_target_v1 (target 1): on T0 Land cells with coverage >= coverage_min (not rim) and land posts on the solve lattice, M, T and L are the PLAIN means over the cell's own land posts (edifices excluded; the ADR 043 F6 literal, per T0 cell) instead of the hat-weighted means; other cells keep the hat-weighted means".into(),
            kkt_rule: "bound_kkt_stop_v1: stop when every solved non-rim cell has |M - T| <= tolerance, or E <= e_min with M - T > 0, or E >= e_max with M - T < 0, AND every solved coastal cell has |land - bilinear land| <= the F4 tolerance, or its flooding offset F = flood_max with land in excess, or F = 0 with land short (projected fixed point of the bound-constrained E and F solve); converged (strict) is still max |M - T| <= tolerance and the F4 rule".into(),
            conv: None,
        }
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.update > 1 || self.kkt > 1 { return Err("es: update and kkt are 0 or 1".into()); }
        if !(1..=2).contains(&self.rim_mode) { return Err("es: rim_mode is 1 or 2".into()); }
        if self.target > 1 { return Err("es: target is 0 or 1".into()); }
        if self.g_min_permille <= 0 || self.g_max_permille < self.g_min_permille { return Err("es: g bounds".into()); }
        if self.step_max_permille <= 0 || !(0..=1000).contains(&self.coverage_min_permille) { return Err("es: step_max / coverage_min".into()); }
        Ok(())
    }
}

/// Per-cell secant state.
#[derive(Clone, Copy, Debug)]
pub struct SecantState { pub x: f64, pub y: f64, pub g: f64, pub sign: f64,
    /// Fallback-step multiplier (halved on each residual sign change while in the +-ln 2 fallback, floor 1/32).
    pub fb: f64 }

/// One secant_diagonal_v1 step for one cell. `num` = T - L, `den` = M - L, `diff` = M - T, `x` = ln E.
/// Returns (step in ln E, the new state's g).
pub fn secant_step(p: &EsParamsV1, st: Option<&SecantState>, x: f64, num: f64, den: f64, diff: f64) -> (f64, Option<SecantState>) {
    let (gmin, gmax) = (p.g_min_permille as f64 / 1000.0, p.g_max_permille as f64 / 1000.0);
    let smax = p.step_max_permille as f64 / 1000.0;
    let sign = if diff.abs() < 1e-9 { 0.0 } else { diff.signum() };
    if !(num > 1.0 && den > 1.0) {
        // Fallback (T or M within 1 m of L): +-ln 2 by the sign of M - T, the multiplier halved on every sign change
        // (the default Jacobi's CF sign damping, so a cell straddling L cannot flip between two E values forever).
        let l2 = 2f64.ln();
        let mut fb = st.map(|s0| s0.fb).unwrap_or(1.0);
        if let Some(s0) = st { if s0.sign * sign < 0.0 { fb = (fb * 0.5).max(1.0 / 32.0); } }
        let s = (if diff < 0.0 { l2 } else { -l2 }) * fb;
        // Keep g; y is undefined (den <= 1), so the next secant restarts from this point.
        return (s.clamp(-smax, smax), Some(SecantState { x, y: f64::NAN, g: st.map(|s0| s0.g).unwrap_or(1.0 / 0.6), sign, fb }));
    }
    let y = den.ln();
    let mut g = st.map(|s| s.g).unwrap_or(1.0 / 0.6);
    if let Some(s0) = st {
        let dx = x - s0.x;
        if dx.abs() > 0.01 && s0.y.is_finite() {
            let gobs = (y - s0.y) / dx;
            if gobs > 0.0 { g = (s0.g * gobs.clamp(gmin, gmax)).sqrt().clamp(gmin, gmax); }
        }
    }
    let mut step = (num / den).ln() / g;
    if let Some(s0) = st { if s0.sign * sign < 0.0 { step *= 0.5; } }
    (step.clamp(-smax, smax), Some(SecantState { x, y, g, sign, fb: st.map(|s0| s0.fb).unwrap_or(1.0) }))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cell whose relief responds with elasticity 0.25 (M - L = c E^0.25) converges in a few secant steps, where the
    /// default 0.6 ln step needs far more (the GT1 / C6 slow class).
    #[test]
    fn secant_converges_on_low_elasticity_cell() {
        let p = EsParamsV1::v1();
        let (l, t, c, g) = (500.0f64, 1500.0f64, 120.0f64, 0.25f64);
        let m_of = |e: f64| l + c * e.powf(g);
        let tol = 20.0;
        let mut e = 400.0f64;
        let mut st: Option<SecantState> = None;
        let mut n_sec = 0;
        for k in 0..40 {
            let m = m_of(e);
            if (m - t).abs() <= tol { n_sec = k; break; }
            let (s, ns) = secant_step(&p, st.as_ref(), e.ln(), t - l, m - l, m - t);
            st = ns; e = (e.ln() + s).exp().clamp(5.0, 20_000.0);
            n_sec = k + 1;
        }
        assert!((m_of(e) - t).abs() <= tol, "secant did not converge: E {e} M {}", m_of(e));
        let mut e2 = 400.0f64;
        let mut n_jac = 40;
        for k in 0..40 {
            let m = m_of(e2);
            if (m - t).abs() <= tol { n_jac = k; break; }
            let s = (0.6 * ((t - l) / (m - l)).ln()).clamp(-(2f64.ln()), 2f64.ln());
            e2 = (e2.ln() + s).exp().clamp(5.0, 20_000.0);
        }
        assert!(n_sec < n_jac, "secant {n_sec} vs jacobi {n_jac}");
        assert!(n_sec <= 8, "secant took {n_sec}");
    }

    /// A cell straddling L (fallback branch) with alternating residual sign does not flip between two E values forever.
    #[test]
    fn fallback_sign_changes_damp() {
        let p = EsParamsV1::v1();
        let mut st: Option<SecantState> = None;
        let mut steps = Vec::new();
        for k in 0..8 {
            let diff = if k % 2 == 0 { 30.0 } else { -30.0 };
            let (s, ns) = secant_step(&p, st.as_ref(), 5.0, 0.5, 0.5, diff);
            st = ns; steps.push(s.abs());
        }
        assert!(steps[7] < steps[0] / 8.0, "{steps:?}");
    }

    /// Falsifier: an unreachable target (M - L never below 300 m) drives E to the step clamp direction without NaN.
    #[test]
    fn secant_unreachable_stays_finite() {
        let p = EsParamsV1::v1();
        let (l, t) = (1000.0f64, 1100.0f64);
        let m_of = |e: f64| l + 300.0 + e * 0.01;
        let mut e = 300.0f64;
        let mut st: Option<SecantState> = None;
        for _ in 0..30 {
            let m = m_of(e);
            let (s, ns) = secant_step(&p, st.as_ref(), e.ln(), t - l, m - l, m - t);
            assert!(s.is_finite());
            st = ns; e = (e.ln() + s).exp().clamp(5.0, 20_000.0);
        }
        assert_eq!(e, 5.0);
    }
}
