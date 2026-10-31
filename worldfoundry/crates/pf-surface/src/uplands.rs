//! SEED-UPLANDS (lane UP, 2026-10-13; design: old ranges wear down into
//! 1-2 km uplands with isostatic rebound, and ranges taper through
//! foothills instead of standing as one step): law id
//! `critical-slope-envelope-v1`, a coarse threshold-hillslope rule. Opt-in
//! only; no legacy path calls this module. Constants and predictions were
//! fixed before any run in `unpublished-design-record`.
//!
//! THE LAW. Cell means may not stand higher above the lowest ground around
//! them than a critical slope allows:
//!     z_env(i) = min_j ( b_j + S_c * path_distance(i, j) ),
//! with b_j the cell-mean elevation of land cells and the fill (water)
//! level of wet cells, path distances along the 8-neighbour grid (Dijkstra
//! from the lowest cells). A land cell above its envelope erodes at
//!     E_th = (z - z_env) / tau,   tau = dx / v_kp,
//! (dx the linear cell size, v_kp a knickpoint celerity: the time for an
//! over-steepened front to be cut back across one cell). The threshold is
//! a SURFACE condition, so the crust removal is sized for the surface
//! AFTER isostatic rebound: with a = (rho_m - rho_c)/rho_m the surface
//! lowering per metre of crust removed (local Airy), backward Euler with
//! z_env held for the step,
//!     d_th = (z - z_env) / a * dt / (tau + dt).
//! (Rule 2, 2026-10-13 09:02 EDT: rule 1 omitted 1/a and was inert,
//! PREREG.txt.)
//! The cell-mean target handed to the shared erosion stage is
//! max(stream-power target, d_th) (`relief::subgrid_step_with_floor`): the
//! faster process sets the rate, nothing is double counted, and transport,
//! deposition and the sediment ledgers stay with the erosion stage.
//! Isostasy (the two-column solve after erosion) returns most of every
//! eroded metre as rock uplift, so a front loses its excess over the
//! envelope gradually (e-folding ~ tau / (1 - rho_c/rho_m)), a high young
//! belt erodes fast at its edges, and a belt ends as a ridge tapering
//! through foothill cells instead of a plateau behind one step.
//!
//! SCALE. Nothing here assumes 70 km cells: S_c is the critical slope OF
//! CELL MEANS at the grid's scale (about 0.6-0.7 for hillslopes, about 0.04
//! for 0.625 deg means on Earth; relief does not scale linearly with the
//! window, so each scale declares its own S_c), and tau shrinks with dx.
//! The same functions can be reused on a finer regional grid.
//!
//! LIMITS (declared): one critical slope everywhere (no lithology, climate
//! or glacial dependence); the envelope is recomputed each step from the
//! step-start surface; wet cells are never lowered here.

use pf_grid::GridSpec;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Law id (the typed law set lives in pf-contracts laws.rs).
pub const CRITICAL_SLOPE_LAW_ID: &str = "critical-slope-envelope-v1";

/// Resolved parameters of `critical-slope-envelope-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CriticalSlopeParams {
    /// Critical slope of cell means at the grid's scale (dimensionless).
    pub critical_slope: f64,
    /// Knickpoint celerity v_kp (m/My); tau = dx / v_kp.
    pub knickpoint_celerity_m_per_my: f64,
}

impl CriticalSlopeParams {
    pub fn validate(&self) -> Result<(), String> {
        if !(self.critical_slope.is_finite() && self.critical_slope > 0.0 && self.critical_slope <= 2.0) {
            return Err(format!(
                "critical_slope {} outside (0, 2]",
                self.critical_slope
            ));
        }
        if !(self.knickpoint_celerity_m_per_my.is_finite() && self.knickpoint_celerity_m_per_my > 0.0) {
            return Err(format!(
                "knickpoint_celerity_m_per_my {} must be finite and > 0",
                self.knickpoint_celerity_m_per_my
            ));
        }
        Ok(())
    }

    /// Response time tau = dx / v_kp (My) for a cell of linear size dx (m).
    pub fn response_time_my(&self, dx_m: f64) -> f64 {
        dx_m / self.knickpoint_celerity_m_per_my
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Key(f64, usize);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    /// Reversed so `BinaryHeap` pops the LOWEST level first; ties break to
    /// the lowest cell id (deterministic).
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .0
            .partial_cmp(&self.0)
            .expect("envelope levels are finite")
            .then_with(|| other.1.cmp(&self.1))
    }
}

/// 8-neighbour adjacency with centre distances (m), from the grid.
pub fn adjacency(grid: &GridSpec) -> Result<Vec<Vec<(usize, f64)>>, String> {
    let n = grid.cells();
    let mut adj = Vec::with_capacity(n);
    for c in 0..n {
        let (i, j) = grid.ij(c as u64).map_err(|e| e.to_string())?;
        let nb = grid.neighbors(i, j).map_err(|e| e.to_string())?;
        adj.push(
            nb.into_iter()
                .map(|(ti, tj, l)| (grid.flat(ti, tj) as usize, l.distance_m))
                .collect(),
        );
    }
    Ok(adj)
}

/// The critical-slope envelope z_env(i) = min_j (b_j + S_c * path(i, j))
/// (Dijkstra from the lowest cells). `z_env <= base` everywhere.
pub fn critical_slope_envelope(
    adj: &[Vec<(usize, f64)>],
    base_m: &[f64],
    critical_slope: f64,
) -> Result<Vec<f64>, String> {
    let n = base_m.len();
    if adj.len() != n {
        return Err(format!("envelope adjacency {} != {n}", adj.len()));
    }
    if let Some(i) = base_m.iter().position(|v| !v.is_finite()) {
        return Err(format!("envelope base[{i}] is not finite"));
    }
    let mut env = base_m.to_vec();
    let mut done = vec![false; n];
    let mut heap: BinaryHeap<Key> = (0..n).map(|i| Key(env[i], i)).collect();
    while let Some(Key(level, i)) = heap.pop() {
        if done[i] || level > env[i] {
            continue;
        }
        done[i] = true;
        for &(j, d) in &adj[i] {
            let cand = level + critical_slope * d;
            if cand < env[j] {
                env[j] = cand;
                heap.push(Key(cand, j));
            }
        }
    }
    Ok(env)
}

/// One step of the threshold term.
#[derive(Debug, Clone, PartialEq)]
pub struct ThresholdStep {
    /// Envelope used this step (m).
    pub envelope_m: Vec<f64>,
    /// Backward-Euler threshold lowering d_th per cell (m, >= 0; 0 on wet
    /// cells). This is the FLOOR handed to `subgrid_step_with_floor`.
    pub lowering_m: Vec<f64>,
    /// Land cells above their envelope.
    pub cells_above: u64,
    /// Largest excess z - z_env over land cells (m).
    pub max_excess_m: f64,
}

/// RULE 1 (inert, superseded 2026-10-13 09:02 EDT): `threshold_lowering_isostatic`
/// with a = 1. Kept only so the applied wiring request 01 builds until
/// request 02 switches the call; not used by any profile after that.
pub fn threshold_lowering(
    adj: &[Vec<(usize, f64)>],
    z_mean_m: &[f64],
    wet: &[bool],
    fill_level_m: &[f64],
    cell_size_m: &[f64],
    p: &CriticalSlopeParams,
    dt_my: f64,
) -> Result<ThresholdStep, String> {
    threshold_lowering_isostatic(adj, z_mean_m, wet, fill_level_m, cell_size_m, p, 1.0, dt_my)
}

/// The threshold lowering of `critical-slope-envelope-v1` for one step.
/// `z_mean_m` is the step-start cell-mean surface, `wet` marks ocean and
/// lake cells whose base is `fill_level_m`, `cell_size_m` = sqrt(area),
/// `isostatic_lowering_per_m` = a = (rho_m - rho_c)/rho_m in (0, 1].
#[allow(clippy::too_many_arguments)]
pub fn threshold_lowering_isostatic(
    adj: &[Vec<(usize, f64)>],
    z_mean_m: &[f64],
    wet: &[bool],
    fill_level_m: &[f64],
    cell_size_m: &[f64],
    p: &CriticalSlopeParams,
    isostatic_lowering_per_m: f64,
    dt_my: f64,
) -> Result<ThresholdStep, String> {
    p.validate()?;
    let a = isostatic_lowering_per_m;
    if !(a.is_finite() && a > 0.0 && a <= 1.0) {
        return Err(format!("isostatic_lowering_per_m {a} outside (0, 1]"));
    }
    let n = z_mean_m.len();
    for (what, len) in [
        ("wet", wet.len()),
        ("fill", fill_level_m.len()),
        ("cell_size", cell_size_m.len()),
    ] {
        if len != n {
            return Err(format!("threshold input {what}: {len} != {n}"));
        }
    }
    if !(dt_my > 0.0) {
        return Err("threshold step: dt must be positive".into());
    }
    let base: Vec<f64> = (0..n)
        .map(|i| if wet[i] { fill_level_m[i] } else { z_mean_m[i] })
        .collect();
    let env = critical_slope_envelope(adj, &base, p.critical_slope)?;
    let mut low = vec![0.0f64; n];
    let mut cells = 0u64;
    let mut max_x = 0.0f64;
    for i in 0..n {
        if wet[i] {
            continue;
        }
        let x = z_mean_m[i] - env[i];
        if x > 0.0 {
            cells += 1;
            max_x = max_x.max(x);
            let tau = p.response_time_my(cell_size_m[i]);
            low[i] = x / a * dt_my / (tau + dt_my);
        }
    }
    Ok(ThresholdStep {
        envelope_m: env,
        lowering_m: low,
        cells_above: cells,
        max_excess_m: max_x,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::relief::{subgrid_step, subgrid_step_with_floor, NonlinearSubgridParams, SubgridInputs};

    fn line_adj(n: usize, dx: f64) -> Vec<Vec<(usize, f64)>> {
        (0..n)
            .map(|i| {
                let mut v = Vec::new();
                if i > 0 {
                    v.push((i - 1, dx));
                }
                if i + 1 < n {
                    v.push((i + 1, dx));
                }
                v
            })
            .collect()
    }

    fn params() -> CriticalSlopeParams {
        CriticalSlopeParams {
            critical_slope: 0.04,
            knickpoint_celerity_m_per_my: 70_000.0,
        }
    }

    #[test]
    fn envelope_tapers_a_plateau_through_foothills() {
        // 9 cells of 70 km: lowland 200 m, a 5-cell plateau at 6 km.
        let dx = 70_000.0;
        let z = [200.0, 200.0, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 200.0, 200.0];
        let env = critical_slope_envelope(&line_adj(9, dx), &z, 0.04).unwrap();
        let s = 0.04 * dx; // 2800 m per cell
        assert_eq!(env[0], 200.0);
        assert!((env[2] - (200.0 + s)).abs() < 1e-9);
        assert!((env[3] - (200.0 + 2.0 * s)).abs() < 1e-9);
        assert_eq!(env[4], 6000.0); // 3 cells in: 8600 m allowed, not binding
        assert!((env[6] - (200.0 + s)).abs() < 1e-9);
        for i in 0..9 {
            assert!(env[i] <= z[i]);
        }
    }

    #[test]
    fn threshold_lowering_is_backward_euler_toward_the_envelope() {
        let dx = 70_000.0;
        let z = [0.0, 5000.0, 300.0];
        let wet = [true, false, false];
        let fill = [0.0, 5000.0, 300.0];
        let size = [dx; 3];
        let p = params();
        let t = threshold_lowering_isostatic(&line_adj(3, dx), &z, &wet, &fill, &size, &p, 0.5, 2.0).unwrap();
        // Envelope of the middle cell: min(0 + 2800, 300 + 2800) = 2800.
        assert!((t.envelope_m[1] - 2800.0).abs() < 1e-9);
        let tau = 1.0; // 70 km / 70 km per My
        assert!((t.lowering_m[1] - 2200.0 / 0.5 * 2.0 / (tau + 2.0)).abs() < 1e-9);
        assert_eq!(t.lowering_m[0], 0.0);
        assert_eq!(t.lowering_m[2], 0.0);
        assert_eq!(t.cells_above, 1);
        assert!((t.max_excess_m - 2200.0).abs() < 1e-9);
        // Loud on bad input.
        let mut bad = p;
        bad.critical_slope = 0.0;
        assert!(threshold_lowering_isostatic(&line_adj(3, dx), &z, &wet, &fill, &size, &bad, 0.5, 2.0).is_err());
        assert!(threshold_lowering_isostatic(&line_adj(3, dx), &z, &wet, &fill, &size, &p, 0.5, 0.0).is_err());
        assert!(threshold_lowering_isostatic(&line_adj(3, dx), &z, &wet, &fill, &size, &p, 0.0, 2.0).is_err());
        assert!(threshold_lowering_isostatic(&line_adj(3, dx), &z, &wet, &fill, &size, &p, 1.5, 2.0).is_err());
    }

    /// Two land cells and an ocean in a line: cell 0 -> 1 -> 2 (ocean).
    fn chain<'a>(
        flow: &'a [Option<usize>],
        area: &'a [f64],
        reach: &'a [f64],
        wet: &'a [bool],
        fill: &'a [f64],
        size: &'a [f64],
    ) -> SubgridInputs<'a> {
        SubgridInputs {
            receiver_first: &[2, 1, 0],
            flow_to: flow,
            contributing_area_m2: area,
            reach_length_m: reach,
            wet,
            fill_level_m: fill,
            cell_size_m: size,
        }
    }

    #[test]
    fn floor_none_is_subgrid_step_and_floor_is_realised() {
        let dx = 70_000.0;
        let flow = [Some(1usize), Some(2), None];
        let area = [5e9, 1e10, 1e10];
        let reach = [dx, dx, 1.0];
        let wet = [false, false, true];
        let z = [6000.0, 2500.0, 0.0];
        let fill = [6000.0, 2500.0, 0.0];
        let size = [dx; 3];
        let relief = [100.0, 100.0, 0.0];
        for n in [1.0, 2.0] {
            let rp = NonlinearSubgridParams {
                k_per_my: 2.0,
                slope_exponent_n: n,
                reference_steepness_m: 500.0,
                channel_head_length_m: 100.0,
            };
            let inp = chain(&flow, &area, &reach, &wet, &fill, &size);
            let a = subgrid_step(&inp, &z, &relief, &rp, 2.0, 0.5, 2.0).unwrap();
            let b = subgrid_step_with_floor(&inp, &z, &relief, &rp, 2.0, 0.5, 2.0, None).unwrap();
            assert_eq!(a, b);
            let zero = [0.0; 3];
            let c = subgrid_step_with_floor(&inp, &z, &relief, &rp, 2.0, 0.5, 2.0, Some(&zero)).unwrap();
            assert_eq!(a, c);
            // The threshold floor raises the target where it is larger.
            let t = threshold_lowering_isostatic(&line_adj(3, dx), &z, &wet, &fill, &size, &params(), 1.0, 2.0).unwrap();
            let f = subgrid_step_with_floor(&inp, &z, &relief, &rp, 2.0, 0.5, 2.0, Some(&t.lowering_m)).unwrap();
            for i in 0..2 {
                assert!(f.target_lowering_m[i] >= a.target_lowering_m[i] - 1e-12);
            }
            assert!(t.lowering_m[0] > 0.0);
            let gap0 = z[0] - (z[1] - f.target_lowering_m[1]);
            assert!(f.target_lowering_m[0] >= t.lowering_m[0].min(gap0 * (1.0 - 1e-9)) - 1e-9);
            // Realised through the stage's n = 1 update: z' = (z + f b)/(1 + f).
            let fe = 2.0 * f.effective_area_m2[1].sqrt() * 2.0 / dx;
            let realised = z[1] - (z[1] + fe * 0.0) / (1.0 + fe);
            assert!((realised - f.target_lowering_m[1]).abs() < 1e-6, "{realised} vs {}", f.target_lowering_m[1]);
        }
        let bad = [-1.0, 0.0, 0.0];
        let rp = NonlinearSubgridParams {
            k_per_my: 2.0,
            slope_exponent_n: 2.0,
            reference_steepness_m: 500.0,
            channel_head_length_m: 100.0,
        };
        let inp = chain(&flow, &area, &reach, &wet, &fill, &size);
        assert!(subgrid_step_with_floor(&inp, &z, &relief, &rp, 2.0, 0.5, 2.0, Some(&bad)).is_err());
    }

    #[test]
    fn envelope_is_deterministic_and_order_free() {
        let adj = line_adj(6, 1000.0);
        let z = [5.0, 900.0, 3.0, 3.0, 800.0, 1.0];
        let a = critical_slope_envelope(&adj, &z, 0.1).unwrap();
        let b = critical_slope_envelope(&adj, &z, 0.1).unwrap();
        assert_eq!(a, b);
        assert_eq!(a, vec![5.0, 103.0, 3.0, 3.0, 101.0, 1.0]);
    }
}
