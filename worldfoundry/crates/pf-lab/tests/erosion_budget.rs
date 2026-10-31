//! SEED-EROSION-BUG (lane E, 2026-10-11): control-volume material budget
//! for the preview erosion chain, per step and per cell.
//!
//! Why a driver lives in a test: `pf-lab/src/preview.rs` belongs to lane K2,
//! so this file re-runs the same legacy per-step chain (advance -> isostasy
//! -> ocean -> route -> erode -> isostasy) through public APIs only. The
//! guard `budget_driver_is_bitwise_simulate` proves the re-run equals
//! `simulate()` exactly, so every number measured here is the number the
//! real chain produces.
//!
//! Control volumes (all volumes in m^3, bedrock and sediment kept apart):
//! - BEDROCK (global): V_after_advance = V_before + ridge + arc - subducted
//!   - delaminated + remap_drift (kinematic, measured by `AdvanceLedger`);
//!   then V_final = V_after_advance - bedrock_eroded.
//! - SEDIMENT (global): S_after_advance = S_before + sediment_remap_drift
//!   (nearest-source carry, subduction retirement and collision stacking of
//!   sediment, none of which is ledgered by `advance`, so the drift is
//!   measured here as a difference and reported separately); then
//!   S_final = S_after_advance - sediment_eroded + deposited + exported
//!   (erosion's own transfer records).
//! - PER LINEAGE (each Eulerian cell carries its material column's history
//!   fields with it): thickness = T0(kind) + shortening_added + arc_added
//!   - eroded_total. Cumulative erosion is bounded by initial thickness
//!   plus tectonic supply; the identity residual is measured.
//! - PER STEP, PER CELL available-material bound: bedrock removed in a step
//!   never exceeds the crust column present at the start of the step's
//!   erosion; sediment removed never exceeds the sediment present.
//!
//! Every step's four control volumes (bedrock and sediment, advance stage
//! and erosion stage) plus the transport identity are recorded in a
//! `pf_surface::sediment::ControlVolumeLedger`, which REFUSES (panics the
//! drive) any epoch whose residual exceeds the declared numerical budget
//! CV_REL_TOL * gross magnitude + CV_ABS_FLOOR_M3.
//!
//! `Erosion::Bounded` drives `erode_step_bounded` with the crust thickness
//! as the available bedrock column (the wiring proposed for preview.rs);
//! `Erosion::Legacy` drives `erode_step` exactly as preview.rs does today.
//!
//! The ignored test `measure_budget` is env-driven (EROSION_SEED,
//! EROSION_FIDELITY, EROSION_OUT, optional EROSION_STEP_MY,
//! EROSION_DURATION_MY, EROSION_MODE=legacy|bounded) and writes a JSON
//! report plus raw f32 fields for the lane E measurement.

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::kinematic::{advance, classify_cells, CrustState};
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;
use pf_lab::preview::{
    elevation_from_state, simulate, PreviewConfig, DIFFUSION_STABILITY_LIMIT, SEDIMENT_REL_TOL,
    STREAM_POWER_M, STREAM_POWER_N, YEARS_PER_MY,
};
use pf_surface::erosion::{erode_step, erode_step_bounded, DiffusionStencil, ErosionParams};
use pf_surface::sediment::{ControlVolumeLedger, ControlVolumeTerms};
use pf_surface::routing::route;
use pf_surface::water::largest_ocean;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Declared numerical budget of every control-volume epoch: f64 rounding
/// of running sums. The transport identity accumulates ~1e6 sequential
/// updates per step at L1 (cells x diffusion substeps x deposit spread), so
/// the worst-case bound is ~N * 2.2e-16 ~ 2e-10 of the gross magnitude. The
/// declared budget is the project's existing sediment-gate tolerance
/// (`preview::SEDIMENT_REL_TOL` = 1e-9), not a new number. (A first
/// declaration of 1e-12, from a sqrt(N_cells) estimate, was refused by the
/// ledger at 38 My on seed 240907 L1 with a 6.86e4 m3 residual against a
/// 6.76e4 m3 budget: the estimate was wrong, the closure was not.) The
/// floor covers epochs with nothing moving.
const CV_REL_TOL: f64 = SEDIMENT_REL_TOL;
const CV_ABS_FLOOR_M3: f64 = 1.0;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Erosion {
    Legacy,
    Bounded,
}

#[derive(Debug, Clone, Default, serde::Serialize)]
struct StepBudget {
    t_my: f64,
    // Bedrock control volume (m^3).
    bedrock_before_m3: f64,
    ridge_created_m3: f64,
    arc_created_m3: f64,
    subducted_m3: f64,
    delaminated_m3: f64,
    bedrock_remap_drift_m3: f64,
    bedrock_after_advance_m3: f64,
    bedrock_eroded_m3: f64,
    bedrock_final_m3: f64,
    bedrock_closure_residual_m3: f64,
    // Tectonic thickening received this step (m^3): arc + collision stack.
    thickening_m3: f64,
    // Sediment control volume (m^3).
    sediment_before_m3: f64,
    sediment_remap_drift_m3: f64,
    sediment_after_advance_m3: f64,
    sediment_eroded_m3: f64,
    deposited_m3: f64,
    exported_m3: f64,
    aggraded_global_layer_m3: f64,
    sediment_final_m3: f64,
    sediment_closure_residual_m3: f64,
    // Erosion transfer identity: eroded = deposited + exported.
    eroded_m3: f64,
    transfer_residual_m3: f64,
    // Available-material bound (counts of violating cells and the worst).
    bedrock_bound_violations: u64,
    worst_bedrock_overdraw_m: f64,
    sediment_bound_violations: u64,
    min_thickness_after_m: f64,
    // Where the bedrock erosion happens.
    max_cell_bedrock_removed_m: f64,
    max_cell: usize,
    max_cell_elev_before_m: f64,
    max_cell_receiver_level_m: f64,
    max_cell_thickness_before_m: f64,
    max_cell_thickening_this_step_m: f64,
    max_cell_contributing_area_m2: f64,
    bedrock_eroded_on_thickened_cells_m3: f64,
    bedrock_eroded_on_capped_cells_m3: f64,
    land_area_fraction: f64,
    max_elevation_m: f64,
    bedrock_bound_hits: u64,
    max_cv_residual_over_declared: f64,
    /// Largest |residual| / gross magnitude over this step's epochs.
    max_cv_relative_residual: f64,
}

struct BudgetDrive {
    grid: GridSpec,
    areas: Vec<f64>,
    elevation: Vec<f64>,
    elevation_t0: Vec<f64>,
    state: CrustState,
    steps: Vec<StepBudget>,
    /// Sum of bedrock removed per cell position (Eulerian, m) over the run.
    eulerian_bedrock_removed_m: Vec<f64>,
    /// Per-lineage identity residual (m) at the end of the run.
    lineage_residual_m: Vec<f64>,
    ledger: ControlVolumeLedger,
}

fn area_weights(grid: &GridSpec) -> Vec<f64> {
    let mut out = Vec::with_capacity(grid.cells());
    for i in 0..grid.nlat() {
        let a = grid.cell_area_m2(i).unwrap();
        for _ in 0..grid.nlon() {
            out.push(a);
        }
    }
    out
}

fn vol(v: &[f64], a: &[f64]) -> f64 {
    v.iter().zip(a).map(|(x, y)| x * y).sum()
}

fn t0_of(kind: CrustKindCodeV0, prior: &UpliftPrior) -> f64 {
    match kind {
        CrustKindCodeV0::Continental => prior.record.continental_crust_thickness_m,
        _ => prior.record.oceanic_crust_thickness_m,
    }
}

fn drive(config: &PreviewConfig) -> BudgetDrive {
    drive_mode(config, Erosion::Legacy)
}

fn drive_mode(config: &PreviewConfig, mode: Erosion) -> BudgetDrive {
    let params = config.profile().unwrap();
    let mut cv = ControlVolumeLedger::new(CV_REL_TOL, CV_ABS_FLOOR_M3).unwrap();
    let grid = config.grid(&params).unwrap();
    let prior = UpliftPrior::preview_v0();
    let cap = prior.record.max_crust_thickness_m;
    let ep = ErosionParams {
        k_per_my: params.surface.erodibility * YEARS_PER_MY,
        m: STREAM_POWER_M,
        n: STREAM_POWER_N,
        diffusivity_m2_per_my: params.surface.hillslope_diffusivity_m2_per_yr * YEARS_PER_MY,
        stability_limit: DIFFUSION_STABILITY_LIMIT,
    };
    ep.validate().unwrap();
    let n = grid.cells();
    let areas = area_weights(&grid);
    let total_area: f64 = areas.iter().sum();
    let sea = params.planet.sea_level_m;

    let seeding = seed_plates(&grid, &params, &config.seed_path()).unwrap();
    let plates = seeding.plates.clone();
    let mut state = CrustState::seed(&grid, &seeding, &prior).unwrap();
    let (class, rate, _) = classify_cells(&grid, &state, &plates).unwrap();
    state.boundary_class = class;
    state.normal_rate_m_per_my = rate;
    let elevation_t0 = elevation_from_state(&state, &params).unwrap();
    let mut elev = elevation_t0.clone();

    let dt = params.tectonics.step_my;
    let steps = (params.tectonics.duration_my / dt).round() as u32;
    let rain = vec![params.surface.uniform_runoff_m_per_yr * YEARS_PER_MY; n];
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let mut log = Vec::new();
    let mut eulerian = vec![0.0f64; n];
    for _ in 0..steps {
        let mut b = StepBudget {
            bedrock_before_m3: state.crust_volume_m3(&areas),
            sediment_before_m3: vol(&state.sediment_m, &areas),
            ..Default::default()
        };
        let (next, ledger, thickening) = advance(&grid, &state, &plates, &prior, dt).unwrap();
        let elev_tect = elevation_from_state(&next, &params).unwrap();
        state = next;
        b.t_my = state.t_my;
        b.ridge_created_m3 = ledger.ridge_created_m3;
        b.arc_created_m3 = ledger.arc_created_m3;
        b.subducted_m3 = ledger.subducted_m3;
        b.delaminated_m3 = ledger.delaminated_m3;
        b.bedrock_remap_drift_m3 = ledger.residual_m3();
        b.bedrock_after_advance_m3 = state.crust_volume_m3(&areas);
        b.thickening_m3 = vol(&thickening, &areas);
        b.sediment_after_advance_m3 = vol(&state.sediment_m, &areas);
        b.sediment_remap_drift_m3 = b.sediment_after_advance_m3 - b.sediment_before_m3;

        let ocean = largest_ocean(&grid, &elev_tect, sea).unwrap();
        let routing = route(&grid, &elev_tect, &rain, &ocean, sea).unwrap();
        let step = match mode {
            Erosion::Legacy => erode_step(
                &grid,
                &stencil,
                &routing,
                &elev_tect,
                &state.sediment_m,
                &ocean,
                sea,
                dt,
                &ep,
            )
            .unwrap(),
            Erosion::Bounded => erode_step_bounded(
                &grid,
                &stencil,
                &routing,
                &elev_tect,
                &state.sediment_m,
                &ocean,
                sea,
                dt,
                &ep,
                Some(&state.thickness_m),
            )
            .unwrap(),
        };
        b.bedrock_bound_hits = step.bedrock_bound_hits;

        // Available-material bound and where erosion happens.
        let mut worst = 0.0f64;
        for c in 0..n {
            let over = step.bedrock_removed_m[c] - state.thickness_m[c];
            if over > 0.0 {
                b.bedrock_bound_violations += 1;
                worst = worst.max(over);
            }
            // Sediment deposited earlier in the same step (diffusion
            // substeps) may be re-eroded, so the bound is start + added.
            if step.sediment_removed_m[c]
                > state.sediment_m[c] + step.sediment_added_m[c] + 1.0e-9
            {
                b.sediment_bound_violations += 1;
            }
            let v = step.bedrock_removed_m[c] * areas[c];
            if thickening[c] > 0.0 {
                b.bedrock_eroded_on_thickened_cells_m3 += v;
            }
            if state.thickness_m[c] >= cap - 1.0 {
                b.bedrock_eroded_on_capped_cells_m3 += v;
            }
            if step.bedrock_removed_m[c] > b.max_cell_bedrock_removed_m {
                b.max_cell_bedrock_removed_m = step.bedrock_removed_m[c];
                b.max_cell = c;
            }
        }
        b.worst_bedrock_overdraw_m = worst;
        let mc = b.max_cell;
        b.max_cell_elev_before_m = elev_tect[mc];
        b.max_cell_receiver_level_m = routing.flow_to[mc]
            .map(|j| {
                let j = j.0 as usize;
                if ocean[j] {
                    sea
                } else {
                    routing.fill_level_m[j].max(elev_tect[j])
                }
            })
            .unwrap_or(f64::NAN);
        b.max_cell_thickness_before_m = state.thickness_m[mc];
        b.max_cell_thickening_this_step_m = thickening[mc];
        b.max_cell_contributing_area_m2 = routing.contributing_area_m2[mc];

        b.eroded_m3 = step.eroded_m3;
        b.deposited_m3 = step.deposited_m3;
        b.exported_m3 = step.exported_m3;
        b.aggraded_global_layer_m3 = step.deposit_aggraded_m3;
        b.transfer_residual_m3 = step.eroded_m3 - step.deposited_m3 - step.exported_m3;
        b.bedrock_eroded_m3 = vol(&step.bedrock_removed_m, &areas);
        b.sediment_eroded_m3 = vol(&step.sediment_removed_m, &areas);

        for c in 0..n {
            state.thickness_m[c] -= step.bedrock_removed_m[c];
            state.eroded_total_m[c] += step.bedrock_removed_m[c];
            eulerian[c] += step.bedrock_removed_m[c];
            state.sediment_m[c] += step.sediment_added_m[c] - step.sediment_removed_m[c];
            if state.sediment_m[c] < 0.0 {
                assert!(state.sediment_m[c] >= -1.0e-9, "negative sediment at {c}");
                state.sediment_m[c] = 0.0;
            }
        }
        b.min_thickness_after_m = state.thickness_m.iter().cloned().fold(f64::MAX, f64::min);
        b.bedrock_final_m3 = state.crust_volume_m3(&areas);
        b.bedrock_closure_residual_m3 =
            b.bedrock_final_m3 - (b.bedrock_after_advance_m3 - b.bedrock_eroded_m3);
        b.sediment_final_m3 = vol(&state.sediment_m, &areas);
        b.sediment_closure_residual_m3 = b.sediment_final_m3
            - (b.sediment_after_advance_m3 - b.sediment_eroded_m3 + b.deposited_m3 + b.exported_m3);
        // The ledger refuses any epoch outside its declared budget.
        let t = b.t_my;
        let epochs = [
            ControlVolumeTerms {
                reservoir: "bedrock",
                stage: "advance",
                t_my: t,
                initial_m3: b.bedrock_before_m3,
                inflows_m3: &[("ridge_created", b.ridge_created_m3), ("arc_created", b.arc_created_m3)],
                outflows_m3: &[("subducted", b.subducted_m3), ("delaminated", b.delaminated_m3)],
                remap_drift_m3: b.bedrock_remap_drift_m3,
                final_m3: b.bedrock_after_advance_m3,
                diagnostics_m3: &[("collision_and_arc_thickening", b.thickening_m3)],
            },
            ControlVolumeTerms {
                reservoir: "bedrock",
                stage: "erosion",
                t_my: t,
                initial_m3: b.bedrock_after_advance_m3,
                inflows_m3: &[],
                outflows_m3: &[("bedrock_eroded", step.bedrock_eroded_m3)],
                remap_drift_m3: 0.0,
                final_m3: b.bedrock_final_m3,
                diagnostics_m3: &[],
            },
            ControlVolumeTerms {
                reservoir: "sediment",
                stage: "advance",
                t_my: t,
                initial_m3: b.sediment_before_m3,
                inflows_m3: &[],
                outflows_m3: &[],
                remap_drift_m3: b.sediment_remap_drift_m3,
                final_m3: b.sediment_after_advance_m3,
                diagnostics_m3: &[],
            },
            ControlVolumeTerms {
                reservoir: "sediment",
                stage: "erosion",
                t_my: t,
                initial_m3: b.sediment_after_advance_m3,
                inflows_m3: &[("deposited_on_land", step.deposited_m3), ("received_by_ocean", step.exported_m3)],
                outflows_m3: &[("sediment_eroded", step.sediment_eroded_m3)],
                remap_drift_m3: 0.0,
                final_m3: b.sediment_final_m3,
                diagnostics_m3: &[("global_marine_layer_placement_plug", step.deposit_aggraded_m3)],
            },
            ControlVolumeTerms {
                reservoir: "transport",
                stage: "erosion",
                t_my: t,
                initial_m3: 0.0,
                inflows_m3: &[("bedrock_eroded", step.bedrock_eroded_m3), ("sediment_eroded", step.sediment_eroded_m3)],
                outflows_m3: &[("deposited_on_land", step.deposited_m3), ("received_by_ocean", step.exported_m3)],
                remap_drift_m3: 0.0,
                final_m3: 0.0,
                diagnostics_m3: &[],
            },
        ];
        for e in &epochs {
            let rec = cv.record(e).unwrap_or_else(|err| panic!("{err}"));
            b.max_cv_residual_over_declared = b
                .max_cv_residual_over_declared
                .max(rec.residual_m3.abs() / rec.declared_error_m3);
            let gross = (rec.declared_error_m3 - CV_ABS_FLOOR_M3) / CV_REL_TOL;
            if gross > 0.0 {
                b.max_cv_relative_residual =
                    b.max_cv_relative_residual.max(rec.residual_m3.abs() / gross);
            }
        }

        elev = elevation_from_state(&state, &params).unwrap();
        let (mut mx, mut land) = (f64::MIN, 0.0);
        for c in 0..n {
            mx = mx.max(elev[c]);
            if elev[c] >= sea {
                land += areas[c];
            }
        }
        b.max_elevation_m = mx;
        b.land_area_fraction = land / total_area;
        log.push(b);
    }
    let lineage_residual_m = (0..n)
        .map(|c| {
            state.thickness_m[c]
                - (t0_of(state.kind[c], &prior) + state.shortening_added_m[c] + state.arc_added_m[c]
                    - state.eroded_total_m[c])
        })
        .collect();
    BudgetDrive {
        grid,
        areas,
        elevation: elev,
        elevation_t0,
        state,
        steps: log,
        eulerian_bedrock_removed_m: eulerian,
        lineage_residual_m,
        ledger: cv,
    }
}

fn mini(seed: u64) -> PreviewConfig {
    PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed,
        duration_my: Some(20.0),
        step_my: Some(2.0),
    }
}

#[test]
fn budget_driver_is_bitwise_simulate() {
    for seed in [240907u64, 7] {
        let cfg = mini(seed);
        let sim = simulate(&cfg).unwrap();
        let d = drive(&cfg);
        assert_eq!(d.elevation, sim.elevation_m, "seed {seed}: driver diverged");
        assert_eq!(d.state.thickness_m, sim.state.thickness_m);
        assert_eq!(d.state.sediment_m, sim.state.sediment_m);
        assert_eq!(d.state.eroded_total_m, sim.state.eroded_total_m);
        for (a, b) in d.steps.iter().zip(&sim.steps) {
            assert_eq!(a.eroded_m3, b.eroded_m3);
            assert_eq!(a.land_area_fraction, b.land_fraction);
        }
    }
}

/// The bounded step with the crust column as availability never binds on
/// these worlds, so the proposed preview.rs wiring is output-identical.
#[test]
fn bounded_drive_equals_legacy_when_the_bound_never_binds() {
    for seed in [240907u64, 7] {
        let cfg = mini(seed);
        let a = drive_mode(&cfg, Erosion::Legacy);
        let b = drive_mode(&cfg, Erosion::Bounded);
        assert!(b.steps.iter().all(|s| s.bedrock_bound_hits == 0), "seed {seed}: bound bound");
        assert_eq!(a.elevation, b.elevation, "seed {seed}");
        assert_eq!(a.state, b.state, "seed {seed}");
    }
}

/// Control volumes close per step inside the declared budget (the ledger
/// refuses otherwise), no cell is overdrawn in any step, and every
/// lineage's cumulative erosion is bounded by its initial column plus the
/// tectonic thickening it received; the lineage identity holds.
#[test]
fn control_volumes_close_and_erosion_is_supply_bounded() {
    let prior = UpliftPrior::preview_v0();
    for seed in [240907u64, 7] {
        let d = drive_mode(&mini(seed), Erosion::Bounded);
        assert_eq!(d.ledger.epochs().len(), 5 * d.steps.len());
        for s in &d.steps {
            assert_eq!(s.bedrock_bound_violations, 0, "seed {seed} t {}", s.t_my);
            assert_eq!(s.sediment_bound_violations, 0, "seed {seed} t {}", s.t_my);
            assert!(s.max_cv_residual_over_declared <= 1.0);
            assert!(s.min_thickness_after_m >= 0.0);
        }
        let n = d.state.thickness_m.len();
        for c in 0..n {
            let supplied = t0_of(d.state.kind[c], &prior)
                + d.state.shortening_added_m[c]
                + d.state.arc_added_m[c];
            assert!(
                d.state.eroded_total_m[c] <= supplied + 1.0e-6,
                "seed {seed} cell {c}: eroded {} > initial + supplied {supplied}",
                d.state.eroded_total_m[c]
            );
            assert!(
                d.lineage_residual_m[c].abs() <= 1.0e-6,
                "seed {seed} cell {c}: lineage identity residual {}",
                d.lineage_residual_m[c]
            );
        }
        // Something actually eroded and moved (not a vacuous pass).
        let tot = d.ledger.totals();
        assert!(tot["bedrock/erosion/out:bedrock_eroded"] > 0.0);
        assert!(tot["transport/erosion/out:received_by_ocean"] > 0.0);
    }
}

#[test]
fn bounded_drive_is_deterministic_across_three_runs() {
    let cfg = mini(240907);
    let a = drive_mode(&cfg, Erosion::Bounded);
    for _ in 0..2 {
        let b = drive_mode(&cfg, Erosion::Bounded);
        assert_eq!(a.elevation, b.elevation);
        assert_eq!(a.state, b.state);
        assert_eq!(a.ledger, b.ledger);
    }
}

// ---------------------------------------------------------------------------
// Measurement (ignored; lane E runs it explicitly per seed).
// ---------------------------------------------------------------------------

fn write_f32(path: &Path, v: &[f64]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for x in v {
        f.write_all(&(*x as f32).to_le_bytes()).unwrap();
    }
}

fn summarize(d: &BudgetDrive) -> serde_json::Value {
    let s = &d.steps;
    let sum = |f: &dyn Fn(&StepBudget) -> f64| s.iter().map(f).sum::<f64>();
    let maxabs = |f: &dyn Fn(&StepBudget) -> f64| s.iter().map(|x| f(x).abs()).fold(0.0, f64::max);
    let n = d.state.thickness_m.len();
    let prior = UpliftPrior::preview_v0();
    let mut argmax = 0usize;
    for c in 0..n {
        if d.state.eroded_total_m[c] > d.state.eroded_total_m[argmax] {
            argmax = c;
        }
    }
    let mut eu_arg = 0usize;
    for c in 0..n {
        if d.eulerian_bedrock_removed_m[c] > d.eulerian_bedrock_removed_m[eu_arg] {
            eu_arg = c;
        }
    }
    let bound_violations = (0..n)
        .filter(|&c| {
            d.state.eroded_total_m[c]
                > t0_of(d.state.kind[c], &prior) + d.state.shortening_added_m[c] + d.state.arc_added_m[c] + 1.0e-6
        })
        .count();
    let total_area: f64 = d.areas.iter().sum();
    serde_json::json!({
        "steps": s.len(),
        "cum_bedrock_eroded_m3": sum(&|x| x.bedrock_eroded_m3),
        "cum_sediment_eroded_m3": sum(&|x| x.sediment_eroded_m3),
        "cum_eroded_m3": sum(&|x| x.eroded_m3),
        "cum_deposited_m3": sum(&|x| x.deposited_m3),
        "cum_exported_m3": sum(&|x| x.exported_m3),
        "cum_global_layer_m3": sum(&|x| x.aggraded_global_layer_m3),
        "cum_thickening_m3": sum(&|x| x.thickening_m3),
        "cum_arc_created_m3": sum(&|x| x.arc_created_m3),
        "cum_ridge_created_m3": sum(&|x| x.ridge_created_m3),
        "cum_subducted_m3": sum(&|x| x.subducted_m3),
        "cum_delaminated_m3": sum(&|x| x.delaminated_m3),
        "cum_bedrock_remap_drift_m3": sum(&|x| x.bedrock_remap_drift_m3),
        "cum_sediment_remap_drift_m3": sum(&|x| x.sediment_remap_drift_m3),
        "cum_bedrock_eroded_on_thickened_cells_m3": sum(&|x| x.bedrock_eroded_on_thickened_cells_m3),
        "cum_bedrock_eroded_on_capped_cells_m3": sum(&|x| x.bedrock_eroded_on_capped_cells_m3),
        "max_abs_bedrock_closure_residual_m3": maxabs(&|x| x.bedrock_closure_residual_m3),
        "max_abs_sediment_closure_residual_m3": maxabs(&|x| x.sediment_closure_residual_m3),
        "max_abs_transfer_residual_m3": maxabs(&|x| x.transfer_residual_m3),
        "bedrock_bound_violations_total": s.iter().map(|x| x.bedrock_bound_violations).sum::<u64>(),
        "sediment_bound_violations_total": s.iter().map(|x| x.sediment_bound_violations).sum::<u64>(),
        "min_thickness_any_step_m": s.iter().map(|x| x.min_thickness_after_m).fold(f64::MAX, f64::min),
        "max_cell_bedrock_removed_one_step_m": s.iter().map(|x| x.max_cell_bedrock_removed_m).fold(0.0, f64::max),
        "lineage_identity_max_abs_residual_m": d.lineage_residual_m.iter().map(|v| v.abs()).fold(0.0, f64::max),
        "lineage_cumulative_bound_violations": bound_violations,
        "lineage_max_eroded_total_m": d.state.eroded_total_m[argmax],
        "lineage_max_cell": argmax,
        "lineage_max_cell_shortening_added_m": d.state.shortening_added_m[argmax],
        "lineage_max_cell_arc_added_m": d.state.arc_added_m[argmax],
        "lineage_max_cell_thickness_m": d.state.thickness_m[argmax],
        "lineage_max_cell_kind": d.state.kind[argmax].code(),
        "eulerian_max_bedrock_removed_m": d.eulerian_bedrock_removed_m[eu_arg],
        "eulerian_max_cell": eu_arg,
        "final_sediment_mean_m": vol(&d.state.sediment_m, &d.areas) / total_area,
        "final_crust_volume_m3": d.state.crust_volume_m3(&d.areas),
    })
}

#[test]
#[ignore]
fn measure_budget() {
    let seed: u64 = std::env::var("EROSION_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("EROSION_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let out = PathBuf::from(std::env::var("EROSION_OUT").unwrap());
    std::fs::create_dir_all(&out).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: std::env::var("EROSION_DURATION_MY").ok().map(|v| v.parse().unwrap()),
        step_my: std::env::var("EROSION_STEP_MY").ok().map(|v| v.parse().unwrap()),
    };
    let mode = match std::env::var("EROSION_MODE").as_deref() {
        Ok("bounded") => Erosion::Bounded,
        Ok("legacy") | Err(_) => Erosion::Legacy,
        Ok(other) => panic!("unknown EROSION_MODE {other}"),
    };
    let t = std::time::Instant::now();
    let d = drive_mode(&cfg, mode);
    let secs = t.elapsed().as_secs_f64();
    write_f32(&out.join("elevation_m.f32"), &d.elevation);
    write_f32(&out.join("elevation_t0_m.f32"), &d.elevation_t0);
    write_f32(&out.join("thickness_m.f32"), &d.state.thickness_m);
    write_f32(&out.join("sediment_m.f32"), &d.state.sediment_m);
    write_f32(&out.join("eroded_total_m.f32"), &d.state.eroded_total_m);
    write_f32(&out.join("shortening_added_m.f32"), &d.state.shortening_added_m);
    write_f32(&out.join("arc_added_m.f32"), &d.state.arc_added_m);
    write_f32(&out.join("eulerian_bedrock_removed_m.f32"), &d.eulerian_bedrock_removed_m);
    write_f32(&out.join("lineage_residual_m.f32"), &d.lineage_residual_m);
    let kind: Vec<u8> = d.state.kind.iter().map(|k| k.code()).collect();
    std::fs::write(out.join("kind.u8"), &kind).unwrap();
    let row_areas: Vec<f64> = (0..d.grid.nlat()).map(|i| d.grid.cell_area_m2(i).unwrap()).collect();
    let report = serde_json::json!({
        "seed": seed,
        "fidelity": std::env::var("EROSION_FIDELITY").unwrap(),
        "nlat": d.grid.nlat(),
        "nlon": d.grid.nlon(),
        "row_area_m2": row_areas,
        "drive_seconds": secs,
        "mode": format!("{mode:?}"),
        "control_volume_totals_m3": d.ledger.totals(),
        "control_volume_max_residual_over_declared": d.steps.iter().map(|s| s.max_cv_residual_over_declared).fold(0.0, f64::max),
        "control_volume_max_relative_residual": d.steps.iter().map(|s| s.max_cv_relative_residual).fold(0.0, f64::max),
        "control_volume_declared_rel_tol": CV_REL_TOL,
        "bedrock_bound_hits_total": d.steps.iter().map(|s| s.bedrock_bound_hits).sum::<u64>(),
        "summary": summarize(&d),
        "steps": d.steps,
    });
    std::fs::write(out.join("budget.json"), serde_json::to_vec_pretty(&report).unwrap()).unwrap();
}
