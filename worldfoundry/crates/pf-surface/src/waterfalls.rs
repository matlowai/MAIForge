//! Waterfall detection (brief S4, spec §57.5): knickpoints and perched
//! confluences with CONSTANT thresholds and measured-evidence cause
//! labels. Detection, never decree: `waterfall_wonder` biases causes
//! upstream in geodynamics; this ruler never changes (same detector
//! scores both worlds in the Gate 4 comparison).
//!
//! CONSTANT THRESHOLDS (compile-time, tested as such):
//! - KNICKPOINT_DROP_M: minimum fill-surface drop across a flow step.
//! - KNICKPOINT_SLOPE_RATIO: drop / mean downstream slope must exceed
//!   this multiple for a SHARP profile break (not just a steep world).
//! - PERCH_HEIGHT_M: tributary junction height above the trunk's local
//!   fill level for a perched (hanging valley) confluence.
//!
//! CAUSE RULES (measured evidence only, documented per §57.5):
//! - UpliftStep: sharp gradient STEP with planar continuity — the drop is
//!   concentrated in ONE flow step while upstream/downstream steps are
//!   much gentler (a knickpoint, not a sustained slope).
//! - DeformationBand: a gradient break ALIGNED across multiple adjacent
//!   reaches/cells at similar elevations (a band, not a point): at least
//!   MIN_BAND_CELLS drops within BAND_WINDOW cells at comparable heights.
//! - HangingValley: tributary discharge junction far above the trunk
//!   local base level (perched confluence).

use crate::routing::Routing;
use pf_contracts::domain::{WaterfallCause, WaterfallObservation};
use pf_contracts::spatial::ElementId;
use pf_grid::GridSpec;

fn distance_between(grid: &GridSpec, a: usize, b: usize) -> Result<f64, WaterfallError> {
    let (ai, aj) = grid
        .ij(a as u64)
        .map_err(|e| WaterfallError::Routing(e.to_string()))?;
    let (bi, bj) = grid
        .ij(b as u64)
        .map_err(|e| WaterfallError::Routing(e.to_string()))?;
    grid.center_distance_m(ai, aj, bi, bj)
        .map_err(|e| WaterfallError::Routing(e.to_string()))
}

pub const KNICKPOINT_DROP_M: f64 = 15.0;
pub const KNICKPOINT_SLOPE_RATIO: f64 = 5.0;
pub const PERCH_HEIGHT_M: f64 = 30.0;
/// Channel-continuity tolerance: a hanging valley's own upstream sits at
/// (nearly) its height; a knickpoint's upstream sits near the DOWNSTREAM
/// level (it is a step, not a perched valley).
pub const PERCH_CONTINUITY_M: f64 = 15.0;
pub const MIN_BAND_CELLS: usize = 3;
pub const BAND_WINDOW_CELLS: usize = 12;
pub const BAND_HEIGHT_TOL_M: f64 = 50.0;

#[derive(Debug, Clone, PartialEq)]
pub enum WaterfallError {
    Routing(String),
    NonFinite { index: usize },
}

impl std::fmt::Display for WaterfallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WaterfallError::Routing(why) => write!(f, "routing: {why}"),
            WaterfallError::NonFinite { index } => write!(f, "nonfinite fill level at {index}"),
        }
    }
}

impl std::error::Error for WaterfallError {}

/// Detect waterfalls on a routed terrain. `flow_step_drops` records the
/// fill-surface drop per flow step (index-aligned with cells); cause
/// labels come from the measured rules above.
pub fn detect_waterfalls(
    grid: &GridSpec,
    routing: &Routing,
) -> Result<Vec<WaterfallObservation>, WaterfallError> {
    let n = grid.cells();
    for (i, &f) in routing.fill_level_m.iter().enumerate() {
        if !f.is_finite() {
            return Err(WaterfallError::NonFinite { index: i });
        }
    }
    // Inflow counts (for confluence evidence).
    let mut inflow: Vec<u32> = vec![0; n];
    for t in routing.flow_to.iter().flatten() {
        inflow[t.0 as usize] += 1;
    }
    // Mean downstream slope per cell (fill surface over flow distance).
    let mut mean_slope = vec![0.0f64; n];
    for i in 0..n {
        if let Some(ElementId(j)) = routing.flow_to[i] {
            let d = distance_between(grid, i, j as usize)?.max(1.0);
            mean_slope[i] = (routing.fill_level_m[i] - routing.fill_level_m[j as usize]) / d;
        }
    }
    // Drop per flow step.
    let drop_at: Vec<Option<(usize, f64)>> = (0..n)
        .map(|i| match routing.flow_to[i] {
            Some(ElementId(j)) => {
                let cell = j as usize;
                Some((cell, routing.fill_level_m[i] - routing.fill_level_m[cell]))
            }
            None => None,
        })
        .collect();

    let mut out: Vec<WaterfallObservation> = Vec::new();

    // Pass 1 — perched confluences (CAUSE EVIDENCE DOMINATES: a junction
    // where a tributary enters far above the trunk's local base level is
    // a HangingValley even when the entry drop is sharp; in-channel steps
    // without confluence are the knickpoint cases below).
    let mut hanging: Vec<bool> = vec![false; n];
    let mut inflow_source_level = vec![f64::MIN; n];
    for k in 0..n {
        if let Some(ElementId(t)) = routing.flow_to[k] {
            let t = t as usize;
            inflow_source_level[t] = inflow_source_level[t].max(routing.fill_level_m[k]);
        }
    }
    for i in 0..n {
        let Some(ElementId(j)) = routing.flow_to[i] else {
            continue;
        };
        let j = j as usize;
        if inflow[j] < 2 {
            continue; // no confluence
        }
        let perch = routing.fill_level_m[i] - routing.fill_level_m[j];
        if perch >= PERCH_HEIGHT_M {
            // Channel continuity: a perched VALLEY continues upstream at
            // the tributary's height; a bare step does not.
            let upstream_at_height =
                inflow_source_level[i] >= routing.fill_level_m[i] - PERCH_CONTINUITY_M;
            if upstream_at_height {
                hanging[i] = true;
                out.push(WaterfallObservation {
                    element: ElementId(i as u64),
                    drop_m: perch,
                    cause: WaterfallCause::HangingValley,
                    reach: None,
                });
            }
        }
    }

    // Pass 2 — knickpoints on non-confluence cells.
    let spacing = grid_mean_spacing_m(grid);
    for i in 0..n {
        if hanging[i] {
            continue;
        }
        let Some((j, drop)) = drop_at[i] else {
            continue;
        };
        if drop < KNICKPOINT_DROP_M {
            continue;
        }
        // Sharpness: this step's slope vs the SURROUNDING slopes, EXCLUDING
        // the step pair itself — upstream cells flowing into i, plus the
        // cell downstream of j. A step, not a sustained ramp.
        let mut local: Vec<f64> = Vec::new();
        for k in 0..n {
            if let Some(ElementId(t)) = routing.flow_to[k] {
                if t as usize == i && k != i {
                    local.push(mean_slope[k].abs());
                }
            }
        }
        if let Some((k2, _)) = drop_at[j] {
            local.push(mean_slope[k2].abs());
        }
        let local: Vec<f64> = local.into_iter().filter(|s| *s > 1.0e-9).collect();
        let neighborhood = if local.is_empty() {
            0.0
        } else {
            local.iter().sum::<f64>() / local.len() as f64
        };
        let this_slope = drop / distance_between(grid, i, j)?.max(1.0);
        if neighborhood > 1.0e-9 && this_slope / neighborhood >= KNICKPOINT_SLOPE_RATIO {
            // Band evidence: nearby drops at comparable heights -> band;
            // isolated -> uplift step.
            let mut band_cells = 0usize;
            for k in 0..n {
                if k == i {
                    continue;
                }
                let dist = distance_between(grid, i, k)?;
                if dist > BAND_WINDOW_CELLS as f64 * spacing {
                    continue;
                }
                if let Some((_, d2)) = drop_at[k] {
                    if d2 >= KNICKPOINT_DROP_M
                        && (routing.fill_level_m[k] - routing.fill_level_m[i]).abs()
                            <= BAND_HEIGHT_TOL_M
                    {
                        band_cells += 1;
                    }
                }
            }
            let cause = if band_cells + 1 >= MIN_BAND_CELLS {
                WaterfallCause::DeformationBand
            } else {
                WaterfallCause::UpliftStep
            };
            out.push(WaterfallObservation {
                element: ElementId(i as u64),
                drop_m: drop,
                cause,
                reach: None,
            });
        }
    }
    out.sort_by_key(|w| (w.element, w.drop_m.to_bits()));
    out.dedup_by_key(|w| w.element);
    Ok(out)
}

fn grid_mean_spacing_m(grid: &GridSpec) -> f64 {
    // Mean cell spacing from total circumference-ish: sqrt(total area / n).
    let n = grid.cells() as f64;
    (grid.total_area_m2() / n).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::routing::route;

    #[test]
    fn thresholds_are_compile_time_constants() {
        // The ruler is constant: these are consts, not parameters; the
        // detector signature takes no dial (structural check by reading).
        assert_eq!(KNICKPOINT_DROP_M, 15.0);
        assert_eq!(PERCH_HEIGHT_M, 30.0);
        assert_eq!(KNICKPOINT_SLOPE_RATIO, 5.0);
        // The public API: detect_waterfalls(grid, routing) — no dial
        // parameter exists to pass (compile-time property).
        let f: fn(&GridSpec, &Routing) -> Result<Vec<WaterfallObservation>, WaterfallError> =
            detect_waterfalls;
        let _ = f;
    }

    #[test]
    fn synthetic_knickpoint_detected_uplift_step_smooth_not() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        // Gentle ramp falling toward element 127 with ONE sharp step at 40.
        let mut elev: Vec<f64> = (0..n).map(|i| 1000.0 - (i as f64) * 0.5).collect();
        elev[40] += 60.0; // a step on the ramp: upstream cell much higher
        let mut ocean = vec![false; n];
        ocean[127] = true; // base level at the BOTTOM of the ramp
        let rain = vec![50.0; n];
        let r = route(&grid, &elev, &rain, &ocean, elev[127]).unwrap();
        let falls = detect_waterfalls(&grid, &r).unwrap();
        assert!(
            falls
                .iter()
                .any(|w| w.element.0 == 40 && w.cause == WaterfallCause::UpliftStep),
            "synthetic step at 40 must be an UpliftStep, got {falls:?}"
        );
        // A perfectly smooth ramp produces no falls.
        let smooth: Vec<f64> = (0..n).map(|i| 1000.0 - (i as f64) * 0.5).collect();
        let mut ocean2 = vec![false; n];
        ocean2[127] = true;
        let r2 = route(&grid, &smooth, &rain, &ocean2, smooth[127]).unwrap();
        assert!(
            detect_waterfalls(&grid, &r2).unwrap().is_empty(),
            "smooth world must have no falls"
        );
    }

    #[test]
    fn perched_confluence_detected_hanging_valley() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        // Trunk ramp falling toward 127; a perched tributary head at 41
        // feeds cell 40, which is a real confluence (trunk cell 39 also
        // flows in from upstream).
        // A high plateau (elements 8,9,10 / 24,25,26 / 41,42 in rows
        // 0-2, cols 8-10) drains internally 9 -> 25 -> 41 and pours over
        // its edge into trunk cell 40, which also receives the trunk's
        // own continuation from 39: a perched confluence.
        let mut elev: Vec<f64> = (0..n).map(|i| 500.0 - (i as f64) * 0.3).collect();
        elev[8] = 604.0;
        elev[9] = 604.0;
        elev[10] = 604.0;
        elev[24] = 603.0;
        elev[25] = 602.0;
        elev[26] = 603.0;
        elev[41] = 600.0; // perched edge of the plateau
        elev[42] = 601.0; // closes the eastern side so 41 drains to 40
        let mut ocean = vec![false; n];
        ocean[127] = true;
        let rain = vec![50.0; n];
        let r = route(&grid, &elev, &rain, &ocean, elev[127]).unwrap();
        let falls = detect_waterfalls(&grid, &r).unwrap();
        // The plateau edge pours into the trunk at a perched confluence
        // with upstream channel continuity: SOME cell must carry that
        // causal label at the plateau's magnitude (the exact pour cell is
        // a D8 steepest-descent property, not asserted).
        assert!(
            falls
                .iter()
                .any(|w| w.cause == WaterfallCause::HangingValley && w.drop_m >= 80.0),
            "plateau edge must yield a HangingValley of ~112 m, got {falls:?}"
        );
        // Structural: every HangingValley carries at least the perch height.
        for w in falls
            .iter()
            .filter(|w| w.cause == WaterfallCause::HangingValley)
        {
            assert!(w.drop_m >= PERCH_HEIGHT_M);
        }
    }
}
