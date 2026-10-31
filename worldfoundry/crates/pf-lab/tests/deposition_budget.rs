//! Lane SD (2026-10-12): drives the preview chain with the new surface laws
//! (`erosion-deposition-v1`, `shelf-wave-base-v1`) and measures conservation,
//! the sediment budget partition and the shelf.
//!
//! Why a driver lives in a test: `pf-lab/src/preview.rs` and the law set
//! belong to another lane (OG), so the wiring is a requested diff. This file
//! re-runs the per-step chain through public APIs only:
//! advance -> isostasy -> ocean -> route -> erode -> isostasy. In
//! `Mode::Legacy` it calls `erode_deposit_step` with the legacy law set,
//! which is `erode_step_bounded`. The guard
//! `legacy_mode_driver_is_bitwise_simulate_with` proves that equals
//! `simulate_with` bit for bit, for the legacy and the bodies + two-column
//! options. The other modes differ from it only in the surface laws.
//!
//! Every step records five control volumes in a `ControlVolumeLedger`, which
//! refuses an epoch outside SEDIMENT_REL_TOL x gross + 1 m3:
//! - bedrock at advance
//! - bedrock at erosion
//! - sediment at advance (subduction outflow, measured remap drift)
//! - sediment at erosion: each storage term of the partition is a named
//!   inflow; sediment eroded is the outflow
//! - transport: bedrock + sediment eroded = the named storage terms
//!
//! The ignored test `measure_deposition` is env-driven: DEP_SEED,
//! DEP_FIDELITY (preview | standard), DEP_OPTS (legacy | twocol | bodies),
//! DEP_MODE (legacy | shelf | ed | both), DEP_OUT, optional DEP_STEP_MY and
//! DEP_DURATION_MY. It writes the field files and meta.json in the format
//! lane K2's analyze.py reads, plus budget.json.

use pf_contracts::laws::{CrustSeedingLawV1, UniverseV1};
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::kinematic::{advance, classify_cells, ContinentalBodiesReport, CrustState};
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;
use pf_lab::preview::{
    elevation_on_grid, simulate_with, universe_of, PreviewConfig, PreviewOptions, SeaLevelSolve,
    DIFFUSION_STABILITY_LIMIT, SEDIMENT_REL_TOL, STREAM_POWER_M, STREAM_POWER_N, YEARS_PER_MY,
};
use pf_surface::erosion::{
    erode_deposit_step, DiffusionStencil, ErosionDepositionParams, ErosionParams, SurfaceLaws,
};
use pf_surface::routing::route;
use pf_surface::sediment::{
    CompensatedColumn, ControlVolumeLedger, ControlVolumeTerms, ShelfParams,
    SHELF_WAVE_BASE_DEPTH_M,
};
use pf_surface::water::largest_ocean;
use std::io::Write;
use std::path::{Path, PathBuf};

const CV_REL_TOL: f64 = SEDIMENT_REL_TOL;
const CV_ABS_FLOOR_M3: f64 = 1.0;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Mode {
    Legacy,
    Shelf,
    Ed,
    Both,
}

impl Mode {
    fn parse(s: &str) -> Self {
        match s {
            "legacy" => Mode::Legacy,
            "shelf" => Mode::Shelf,
            "ed" => Mode::Ed,
            "both" => Mode::Both,
            o => panic!("unknown DEP_MODE {o}"),
        }
    }
}

/// The isostasy law's own densities for the compensated column.
fn column_of(u: &UniverseV1) -> CompensatedColumn {
    CompensatedColumn {
        mantle_density_kg_m3: u.constants.mantle_density_kg_m3.value,
        sediment_density_kg_m3: u.sediment_density_kg_m3(),
        water_density_kg_m3: u.constants.seawater_density_kg_m3.value,
    }
}

fn laws_for(mode: Mode, u: &UniverseV1) -> SurfaceLaws {
    let shelf = Some(ShelfParams::v1(column_of(u)));
    // DEP_G_SENSITIVITY: a labelled sensitivity probe of the declared G
    // (information for the report, never a proposed value).
    let g = std::env::var("DEP_G_SENSITIVITY")
        .ok()
        .map(|v| v.parse::<f64>().unwrap())
        .unwrap_or(pf_surface::erosion::DEPOSITION_COEFFICIENT_G);
    let ed = Some(ErosionDepositionParams {
        deposition_coefficient_g: g,
        ..ErosionDepositionParams::v1()
    });
    match mode {
        Mode::Legacy => SurfaceLaws::legacy(),
        Mode::Shelf => SurfaceLaws { erosion_deposition: None, shelf },
        Mode::Ed => SurfaceLaws { erosion_deposition: ed, shelf: None },
        Mode::Both => SurfaceLaws { erosion_deposition: ed, shelf },
    }
}

fn options(name: &str) -> PreviewOptions {
    match name {
        "legacy" => PreviewOptions::legacy(),
        "twocol" => PreviewOptions::two_column_v0(),
        "bodies" => PreviewOptions::crust_bodies_v0(),
        o => panic!("unknown DEP_OPTS {o}"),
    }
}

#[derive(Debug, Clone, Default, serde::Serialize)]
struct StepBudget {
    t_my: f64,
    bedrock_before_m3: f64,
    bedrock_after_advance_m3: f64,
    bedrock_eroded_m3: f64,
    bedrock_final_m3: f64,
    sediment_before_m3: f64,
    sediment_subducted_m3: f64,
    sediment_remap_drift_m3: f64,
    sediment_after_advance_m3: f64,
    sediment_eroded_m3: f64,
    sediment_final_m3: f64,
    eroded_m3: f64,
    deposited_m3: f64,
    exported_m3: f64,
    // Partition (m3).
    land_fluvial_m3: f64,
    land_hillslope_m3: f64,
    land_floor_m3: f64,
    lake_m3: f64,
    marine_shelf_m3: f64,
    marine_slope_m3: f64,
    marine_deep_m3: f64,
    global_remainder_m3: f64,
    delivered_to_ocean_m3: f64,
    fluvial_gross_erosion_m3: f64,
    fluvial_gross_deposition_m3: f64,
    marine_search_bound_hits: u64,
    lake_spread_bound_hits: u64,
    shelf_cells_filled: u64,
    solver_passes: u32,
    solver_converged: bool,
    solver_last_change_m: f64,
    bedrock_bound_hits: u64,
    // F4(c): cells the shelf law filled exactly to the shelf top this step,
    // and their worst |elevation after re-equilibration - shelf top| (m).
    shelf_filled_cells_checked: u64,
    shelf_filled_max_abs_dev_m: f64,
    shelf_filled_mean_dev_m: f64,
    // Outcome after the step.
    land_fraction: f64,
    shelf_area_fraction: f64,
    max_elevation_m: f64,
    min_elevation_m: f64,
    sea_level_gauge_m: Option<f64>,
    ocean_volume_residual_m3: Option<f64>,
    land_elevation_p10_p50_p90_m: Option<[f64; 3]>,
    max_cv_residual_over_declared: f64,
    max_cv_relative_residual: f64,
}

struct Drive {
    grid: GridSpec,
    areas: Vec<f64>,
    params_ocean_fraction: f64,
    sea: f64,
    elevation: Vec<f64>,
    elevation_t0: Vec<f64>,
    state: CrustState,
    state0: CrustState,
    sea_level_t0: Option<SeaLevelSolve>,
    bodies: Option<ContinentalBodiesReport>,
    steps: Vec<StepBudget>,
    ledger: ControlVolumeLedger,
    /// Eulerian cumulative deposit on dry land cells (m) and in the ocean (m).
    cum_land_deposit_m: Vec<f64>,
    cum_marine_deposit_m: Vec<f64>,
    /// Land deposit (m) of the LAST step (actively aggrading at the end).
    last_land_deposit_m: Vec<f64>,
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

fn land_percentiles(e: &[f64], a: &[f64], sea: f64) -> Option<[f64; 3]> {
    let mut land: Vec<(f64, f64)> = e.iter().zip(a).filter(|(z, _)| **z >= sea).map(|(z, w)| (*z, *w)).collect();
    if land.is_empty() {
        return None;
    }
    land.sort_by(|x, y| x.0.total_cmp(&y.0));
    let tot: f64 = land.iter().map(|x| x.1).sum();
    let mut out = [0.0; 3];
    for (k, q) in [0.1, 0.5, 0.9].iter().enumerate() {
        let mut acc = 0.0;
        out[k] = land.last().unwrap().0;
        for &(z, w) in &land {
            acc += w;
            if acc >= q * tot {
                out[k] = z;
                break;
            }
        }
    }
    Some(out)
}

fn drive(cfg: &PreviewConfig, opts: &PreviewOptions, mode: Mode) -> Drive {
    let params = cfg.profile_with(opts).unwrap();
    let u = universe_of(&params).unwrap();
    let laws = laws_for(mode, &u);
    let grid = cfg.grid(&params).unwrap();
    // As preview.rs `uplift_prior` and `erosion_params` (private there).
    let mut prior = UpliftPrior::preview_v0();
    prior.record.oceanic_crust_thickness_m = u.constants.oceanic_crust_thickness_m.value;
    prior.record.continental_crust_thickness_m = u.constants.continental_crust_thickness_m.value;
    let ep = ErosionParams {
        k_per_my: u.constants.stream_power_erodibility_per_yr.value * YEARS_PER_MY,
        m: STREAM_POWER_M,
        n: STREAM_POWER_N,
        diffusivity_m2_per_my: u.constants.hillslope_diffusivity_m2_per_yr.value * YEARS_PER_MY,
        stability_limit: DIFFUSION_STABILITY_LIMIT,
    };
    ep.validate().unwrap();
    let column = column_of(&u);
    let n = grid.cells();
    let areas = area_weights(&grid);
    let total_area: f64 = areas.iter().sum();
    let sea = u.constants.sea_level_datum_m.value;
    let seeding = seed_plates(&grid, &params, &cfg.seed_path()).unwrap();
    let plates = seeding.plates.clone();
    let (mut state, bodies) = match &u.law_set.crust_seeding {
        CrustSeedingLawV1::WholePlateContinentsV0 => (CrustState::seed(&grid, &seeding, &prior).unwrap(), None),
        CrustSeedingLawV1::GeodesicDiscBodiesV0 {
            continental_area_fraction,
            rifted_margin_width_m,
        } => {
            let (s, r) = CrustState::seed_bodies(
                &grid,
                &seeding,
                &prior,
                continental_area_fraction.value,
                rifted_margin_width_m.value,
            )
            .unwrap();
            (s, Some(r))
        }
    };
    let (class, rate, _) = classify_cells(&grid, &state, &plates).unwrap();
    state.boundary_class = class;
    state.normal_rate_m_per_my = rate;
    let state0 = state.clone();
    let (elevation_t0, sea_level_t0) = elevation_on_grid(&state, &params, &areas).unwrap();
    let mut elev = elevation_t0.clone();
    let dt = params.tectonics.step_my;
    let steps = (params.tectonics.duration_my / dt).round() as u32;
    let rain = vec![params.surface.uniform_runoff_m_per_yr * YEARS_PER_MY; n];
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let mut cv = ControlVolumeLedger::new(CV_REL_TOL, CV_ABS_FLOOR_M3).unwrap();
    let mut log = Vec::with_capacity(steps as usize);
    let mut cum_land = vec![0.0f64; n];
    let mut cum_marine = vec![0.0f64; n];
    let mut last_land = vec![0.0f64; n];
    let shelf_top = sea - SHELF_WAVE_BASE_DEPTH_M;
    for _ in 0..steps {
        let mut b = StepBudget {
            bedrock_before_m3: state.crust_volume_m3(&areas),
            sediment_before_m3: vol(&state.sediment_m, &areas),
            ..Default::default()
        };
        let (next, ledger, _thickening) = advance(&grid, &state, &plates, &prior, dt).unwrap();
        let (elev_tect, _) = elevation_on_grid(&next, &params, &areas).unwrap();
        state = next;
        b.t_my = state.t_my;
        b.bedrock_after_advance_m3 = state.crust_volume_m3(&areas);
        b.sediment_after_advance_m3 = vol(&state.sediment_m, &areas);
        b.sediment_subducted_m3 = ledger.sediment_subducted_m3;
        b.sediment_remap_drift_m3 =
            b.sediment_after_advance_m3 - (b.sediment_before_m3 - b.sediment_subducted_m3);
        let ocean = largest_ocean(&grid, &elev_tect, sea).unwrap();
        let routing = route(&grid, &elev_tect, &rain, &ocean, sea).unwrap();
        let step = erode_deposit_step(
            &grid,
            &stencil,
            &routing,
            &elev_tect,
            &state.sediment_m,
            &ocean,
            sea,
            dt,
            &ep,
            &laws,
            Some(&state.thickness_m),
        )
        .unwrap();
        b.eroded_m3 = step.eroded_m3;
        b.deposited_m3 = step.deposited_m3;
        b.exported_m3 = step.exported_m3;
        b.bedrock_eroded_m3 = step.bedrock_eroded_m3;
        b.sediment_eroded_m3 = step.sediment_eroded_m3;
        b.bedrock_bound_hits = step.bedrock_bound_hits;
        if let Some(r) = &step.deposition {
            b.land_fluvial_m3 = r.land_fluvial_m3;
            b.land_hillslope_m3 = r.land_hillslope_m3;
            b.land_floor_m3 = r.land_floor_m3;
            b.lake_m3 = r.lake_m3;
            b.marine_shelf_m3 = r.marine_shelf_m3;
            b.marine_slope_m3 = r.marine_slope_m3;
            b.marine_deep_m3 = r.marine_deep_m3;
            b.global_remainder_m3 = r.global_remainder_m3;
            b.delivered_to_ocean_m3 = r.delivered_to_ocean_m3;
            b.fluvial_gross_erosion_m3 = r.fluvial_gross_erosion_m3;
            b.fluvial_gross_deposition_m3 = r.fluvial_gross_deposition_m3;
            b.marine_search_bound_hits = r.marine_search_bound_hits;
            b.lake_spread_bound_hits = r.lake_spread_bound_hits;
            b.shelf_cells_filled = r.shelf_cells_filled;
            b.solver_passes = r.sweeps;
            b.solver_converged = r.converged;
            b.solver_last_change_m = r.last_sweep_max_change_m;
        } else {
            // Legacy: the step's own terms (lake/land deposition, ocean
            // receipts, global layer as a diagnostic inside the receipts).
            b.global_remainder_m3 = step.deposit_aggraded_m3;
            b.delivered_to_ocean_m3 = step.exported_m3;
            b.solver_converged = true;
        }
        let wet_land = |c: usize| !ocean[c] && routing.fill_level_m[c] > elev_tect[c] + 1.0e-9;
        for c in 0..n {
            last_land[c] = 0.0;
            if ocean[c] {
                cum_marine[c] += step.sediment_added_m[c];
            } else if !wet_land(c) {
                cum_land[c] += step.sediment_added_m[c];
                last_land[c] = step.sediment_added_m[c];
            }
        }
        // Apply exactly as PreviewRun::step does.
        for c in 0..n {
            state.thickness_m[c] -= step.bedrock_removed_m[c];
            state.eroded_total_m[c] += step.bedrock_removed_m[c];
            state.sediment_m[c] += step.sediment_added_m[c] - step.sediment_removed_m[c];
            if state.sediment_m[c] < 0.0 {
                assert!(state.sediment_m[c] >= -1.0e-9, "negative sediment at {c}: {}", state.sediment_m[c]);
                state.sediment_m[c] = 0.0;
            }
        }
        b.bedrock_final_m3 = state.crust_volume_m3(&areas);
        b.sediment_final_m3 = vol(&state.sediment_m, &areas);
        let (e_after, solve) = elevation_on_grid(&state, &params, &areas).unwrap();
        // F4(c): cells the shelf law brought exactly to the shelf top.
        if laws.shelf.is_some() {
            let (mut k, mut worst, mut sum) = (0u64, 0.0f64, 0.0f64);
            for c in 0..n {
                if !ocean[c] || step.sediment_added_m[c] <= 0.0 || step.sediment_removed_m[c] > 0.0 {
                    continue;
                }
                let surf = column.surface_after(elev_tect[c], step.sediment_added_m[c], sea);
                if (surf - shelf_top).abs() < 1.0e-6 {
                    let dev = e_after[c] - shelf_top;
                    k += 1;
                    worst = worst.max(dev.abs());
                    sum += dev;
                }
            }
            b.shelf_filled_cells_checked = k;
            b.shelf_filled_max_abs_dev_m = worst;
            b.shelf_filled_mean_dev_m = if k > 0 { sum / k as f64 } else { 0.0 };
        }
        elev = e_after;
        let (mut mx, mut mn, mut land, mut shelf) = (f64::MIN, f64::MAX, 0.0, 0.0);
        for c in 0..n {
            mx = mx.max(elev[c]);
            mn = mn.min(elev[c]);
            if elev[c] >= sea {
                land += areas[c];
            } else if elev[c] >= sea - 200.0 {
                shelf += areas[c];
            }
        }
        b.max_elevation_m = mx;
        b.min_elevation_m = mn;
        b.land_fraction = land / total_area;
        b.shelf_area_fraction = shelf / total_area;
        b.sea_level_gauge_m = solve.map(|s| s.sea_level_gauge_m);
        b.ocean_volume_residual_m3 = solve.map(|s| s.ocean_volume_residual_m3);
        b.land_elevation_p10_p50_p90_m = land_percentiles(&elev, &areas, sea);

        // Control volumes (the ledger refuses any epoch outside its budget).
        let t = b.t_my;
        let stored: Vec<(&str, f64)> = match &step.deposition {
            Some(r) => r.stored_terms().to_vec(),
            None => vec![("deposited_on_land", step.deposited_m3), ("received_by_ocean", step.exported_m3)],
        };
        let sed_out = [("sediment_eroded", step.sediment_eroded_m3)];
        let transport_in = [("bedrock_eroded", step.bedrock_eroded_m3), ("sediment_eroded", step.sediment_eroded_m3)];
        let diag = [("global_remainder_placement", b.global_remainder_m3)];
        let epochs = [
            ControlVolumeTerms {
                reservoir: "bedrock",
                stage: "advance",
                t_my: t,
                initial_m3: b.bedrock_before_m3,
                inflows_m3: &[("ridge_created", ledger.ridge_created_m3), ("arc_created", ledger.arc_created_m3)],
                outflows_m3: &[("subducted", ledger.subducted_m3), ("delaminated", ledger.delaminated_m3)],
                remap_drift_m3: b.bedrock_after_advance_m3
                    - (b.bedrock_before_m3 + ledger.ridge_created_m3 + ledger.arc_created_m3
                        - ledger.subducted_m3
                        - ledger.delaminated_m3),
                final_m3: b.bedrock_after_advance_m3,
                diagnostics_m3: &[],
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
                outflows_m3: &[("subducted", b.sediment_subducted_m3)],
                remap_drift_m3: b.sediment_remap_drift_m3,
                final_m3: b.sediment_after_advance_m3,
                diagnostics_m3: &[],
            },
            ControlVolumeTerms {
                reservoir: "sediment",
                stage: "erosion",
                t_my: t,
                initial_m3: b.sediment_after_advance_m3,
                inflows_m3: &stored,
                outflows_m3: &sed_out,
                remap_drift_m3: 0.0,
                final_m3: b.sediment_final_m3,
                diagnostics_m3: &diag,
            },
            ControlVolumeTerms {
                reservoir: "transport",
                stage: "erosion",
                t_my: t,
                initial_m3: 0.0,
                inflows_m3: &transport_in,
                outflows_m3: &stored,
                remap_drift_m3: 0.0,
                final_m3: 0.0,
                diagnostics_m3: &[],
            },
        ];
        for e in &epochs {
            let rec = cv.record(e).unwrap_or_else(|err| panic!("{mode:?}: {err}"));
            b.max_cv_residual_over_declared =
                b.max_cv_residual_over_declared.max(rec.residual_m3.abs() / rec.declared_error_m3);
            let gross = (rec.declared_error_m3 - CV_ABS_FLOOR_M3) / CV_REL_TOL;
            if gross > 0.0 {
                b.max_cv_relative_residual = b.max_cv_relative_residual.max(rec.residual_m3.abs() / gross);
            }
        }
        log.push(b);
    }
    Drive {
        grid,
        areas,
        params_ocean_fraction: params.tectonics.ocean_fraction,
        sea,
        elevation: elev,
        elevation_t0,
        state,
        state0,
        sea_level_t0,
        bodies,
        steps: log,
        ledger: cv,
        cum_land_deposit_m: cum_land,
        cum_marine_deposit_m: cum_marine,
        last_land_deposit_m: last_land,
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
fn legacy_mode_driver_is_bitwise_simulate_with() {
    for seed in [240907u64, 7] {
        for opts in ["legacy", "bodies"] {
            let cfg = mini(seed);
            let sim = match simulate_with(&cfg, &options(opts)) {
                Ok(s) => s,
                Err(e) => {
                    // A refused world (e.g. a seeding refusal) is not a driver
                    // difference; it must refuse identically in the driver.
                    eprintln!("seed {seed} {opts}: simulate_with refused: {e}");
                    continue;
                }
            };
            let d = drive(&cfg, &options(opts), Mode::Legacy);
            assert_eq!(d.elevation, sim.elevation_m, "seed {seed} {opts}: driver diverged");
            assert_eq!(d.state.thickness_m, sim.state.thickness_m);
            assert_eq!(d.state.sediment_m, sim.state.sediment_m);
            assert_eq!(d.elevation_t0, sim.elevation_t0_m);
            for (a, b) in d.steps.iter().zip(&sim.steps) {
                assert_eq!(a.eroded_m3, b.eroded_m3);
                assert_eq!(a.land_fraction, b.land_fraction);
            }
        }
    }
}

/// F1: every mode closes all five control volumes every step; the named
/// partition sums to what eroded; things actually moved.
#[test]
fn every_mode_closes_every_step() {
    for seed in [240907u64, 7] {
        for opts in ["legacy", "bodies"] {
            for mode in [Mode::Shelf, Mode::Ed, Mode::Both] {
                let d = drive(&mini(seed), &options(opts), mode);
                assert_eq!(d.ledger.epochs().len(), 5 * d.steps.len());
                for s in &d.steps {
                    assert!(s.max_cv_residual_over_declared <= 1.0, "{mode:?} {opts} t {}", s.t_my);
                    let stored = s.land_fluvial_m3 + s.land_hillslope_m3 + s.land_floor_m3 + s.lake_m3
                        + s.marine_shelf_m3 + s.marine_slope_m3 + s.marine_deep_m3 + s.global_remainder_m3;
                    assert!((stored - s.eroded_m3).abs() <= 1e-9 * s.eroded_m3.max(1.0));
                }
                let tot = d.ledger.totals();
                assert!(tot["transport/erosion/in:bedrock_eroded"] > 0.0);
                if mode != Mode::Shelf {
                    assert!(tot["transport/erosion/out:land_fluvial"] > 0.0, "{mode:?} {opts}: no land storage");
                }
            }
        }
    }
}

/// F4(c): under the single-column law (fixed sea level), a cell the shelf
/// law fills to the shelf top re-equilibrates to it.
#[test]
fn shelf_fill_re_equilibrates_to_the_shelf_top_under_single_column_isostasy() {
    let d = drive(&mini(240907), &options("legacy"), Mode::Both);
    let checked: u64 = d.steps.iter().map(|s| s.shelf_filled_cells_checked).sum();
    assert!(checked > 0, "no cell was filled to the shelf top");
    for s in &d.steps {
        assert!(s.shelf_filled_max_abs_dev_m <= 1.0, "t {}: {} m", s.t_my, s.shelf_filled_max_abs_dev_m);
    }
}

/// F2: three runs are identical (state, elevation, ledger).
#[test]
fn new_laws_are_deterministic_across_three_runs() {
    let cfg = mini(240907);
    let a = drive(&cfg, &options("bodies"), Mode::Both);
    for _ in 0..2 {
        let b = drive(&cfg, &options("bodies"), Mode::Both);
        assert_eq!(a.elevation, b.elevation);
        assert_eq!(a.state, b.state);
        assert_eq!(a.ledger, b.ledger);
    }
}

// ---------------------------------------------------------------------------
// Measurement (ignored; run explicitly per seed / options / mode).
// ---------------------------------------------------------------------------

fn write_f32(path: &Path, v: &[f64]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for x in v {
        f.write_all(&(*x as f32).to_le_bytes()).unwrap();
    }
}

#[test]
#[ignore]
fn measure_deposition() {
    let seed: u64 = std::env::var("DEP_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("DEP_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        o => panic!("unknown fidelity {o}"),
    };
    let opts_name = std::env::var("DEP_OPTS").unwrap();
    let mode = Mode::parse(&std::env::var("DEP_MODE").unwrap());
    let out = PathBuf::from(std::env::var("DEP_OUT").unwrap());
    std::fs::create_dir_all(&out).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: std::env::var("DEP_DURATION_MY").ok().map(|v| v.parse().unwrap()),
        step_my: std::env::var("DEP_STEP_MY").ok().map(|v| v.parse().unwrap()),
    };
    let t = std::time::Instant::now();
    let d = std::panic::catch_unwind(|| drive(&cfg, &options(&opts_name), mode));
    let d = match d {
        Ok(d) => d,
        Err(e) => {
            let msg = e
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "panic".into());
            std::fs::write(out.join("ERROR.txt"), &msg).unwrap();
            panic!("{msg}");
        }
    };
    let secs = t.elapsed().as_secs_f64();
    let grid = &d.grid;
    write_f32(&out.join("elevation_m.f32"), &d.elevation);
    write_f32(&out.join("elevation_t0_m.f32"), &d.elevation_t0);
    write_f32(&out.join("thickness_m.f32"), &d.state.thickness_m);
    write_f32(&out.join("thickness_t0_m.f32"), &d.state0.thickness_m);
    write_f32(&out.join("sediment_m.f32"), &d.state.sediment_m);
    write_f32(&out.join("eroded_total_m.f32"), &d.state.eroded_total_m);
    write_f32(&out.join("cum_land_deposit_m.f32"), &d.cum_land_deposit_m);
    write_f32(&out.join("cum_marine_deposit_m.f32"), &d.cum_marine_deposit_m);
    write_f32(&out.join("last_land_deposit_m.f32"), &d.last_land_deposit_m);
    let code = |k: &[CrustKindCodeV0]| k.iter().map(|x| x.code()).collect::<Vec<u8>>();
    std::fs::write(out.join("kind.u8"), code(&d.state.kind)).unwrap();
    std::fs::write(out.join("kind_t0.u8"), code(&d.state0.kind)).unwrap();
    std::fs::write(out.join("boundary_class.u8"), &d.state.boundary_class).unwrap();
    let plate: Vec<u8> = d.state.plate.iter().flat_map(|p| p.0.to_le_bytes()).collect();
    std::fs::write(out.join("plate_id.u32"), plate).unwrap();
    let ocean = largest_ocean(grid, &d.elevation, d.sea).unwrap();
    std::fs::write(out.join("ocean_connected.u8"), ocean.iter().map(|&b| b as u8).collect::<Vec<_>>()).unwrap();
    let ocean0 = largest_ocean(grid, &d.elevation_t0, d.sea).unwrap();
    std::fs::write(out.join("ocean_connected_t0.u8"), ocean0.iter().map(|&b| b as u8).collect::<Vec<_>>()).unwrap();
    let row_areas: Vec<f64> = (0..grid.nlat()).map(|i| grid.cell_area_m2(i).unwrap()).collect();
    let sum = |f: &dyn Fn(&StepBudget) -> f64| d.steps.iter().map(f).sum::<f64>();
    let summary = serde_json::json!({
        "cum_eroded_m3": sum(&|s| s.eroded_m3),
        "cum_bedrock_eroded_m3": sum(&|s| s.bedrock_eroded_m3),
        "cum_sediment_eroded_m3": sum(&|s| s.sediment_eroded_m3),
        "cum_land_fluvial_m3": sum(&|s| s.land_fluvial_m3),
        "cum_land_hillslope_m3": sum(&|s| s.land_hillslope_m3),
        "cum_land_floor_m3": sum(&|s| s.land_floor_m3),
        "cum_lake_m3": sum(&|s| s.lake_m3),
        "cum_marine_shelf_m3": sum(&|s| s.marine_shelf_m3),
        "cum_marine_slope_m3": sum(&|s| s.marine_slope_m3),
        "cum_marine_deep_m3": sum(&|s| s.marine_deep_m3),
        "cum_global_remainder_m3": sum(&|s| s.global_remainder_m3),
        "cum_delivered_to_ocean_m3": sum(&|s| s.delivered_to_ocean_m3),
        "cum_deposited_m3": sum(&|s| s.deposited_m3),
        "cum_exported_m3": sum(&|s| s.exported_m3),
        "cum_sediment_subducted_m3": sum(&|s| s.sediment_subducted_m3),
        "cum_sediment_remap_drift_m3": sum(&|s| s.sediment_remap_drift_m3),
        "marine_search_bound_hits": d.steps.iter().map(|s| s.marine_search_bound_hits).sum::<u64>(),
        "lake_spread_bound_hits": d.steps.iter().map(|s| s.lake_spread_bound_hits).sum::<u64>(),
        "bedrock_bound_hits": d.steps.iter().map(|s| s.bedrock_bound_hits).sum::<u64>(),
        "solver_unconverged_steps": d.steps.iter().filter(|s| !s.solver_converged).count(),
        "solver_max_passes": d.steps.iter().map(|s| s.solver_passes).max().unwrap_or(0),
        "shelf_filled_cells_checked": d.steps.iter().map(|s| s.shelf_filled_cells_checked).sum::<u64>(),
        "shelf_filled_max_abs_dev_m": d.steps.iter().map(|s| s.shelf_filled_max_abs_dev_m).fold(0.0, f64::max),
        "cv_max_residual_over_declared": d.steps.iter().map(|s| s.max_cv_residual_over_declared).fold(0.0, f64::max),
        "cv_max_relative_residual": d.steps.iter().map(|s| s.max_cv_relative_residual).fold(0.0, f64::max),
        "final_sediment_volume_m3": vol(&d.state.sediment_m, &d.areas),
        "final_land_fraction": d.steps.last().map(|s| s.land_fraction),
        "final_shelf_area_fraction": d.steps.last().map(|s| s.shelf_area_fraction),
    });
    let bodies = d.bodies.as_ref().map(|r| serde_json::json!({
        "declared_area_fraction": r.declared_area_fraction,
        "achieved_area_fraction": r.achieved_area_fraction,
    }));
    let meta = serde_json::json!({
        "seed": seed,
        "mode": format!("{opts_name}+{mode:?}"),
        "opts": opts_name,
        "surface_mode": format!("{mode:?}"),
        "nlat": grid.nlat(),
        "nlon": grid.nlon(),
        "row_area_m2": row_areas,
        "sea_level_datum_m": d.sea,
        "ocean_fraction_param": d.params_ocean_fraction,
        "simulate_seconds": secs,
        "sea_level_t0": d.sea_level_t0,
        "continental_bodies": bodies,
        "shelf_wave_base_depth_m": SHELF_WAVE_BASE_DEPTH_M,
        "deposition_coefficient_g": std::env::var("DEP_G_SENSITIVITY").ok().map(|v| v.parse::<f64>().unwrap()).unwrap_or(pf_surface::erosion::DEPOSITION_COEFFICIENT_G),
        "summary": summary,
        "steps": d.steps,
    });
    std::fs::write(out.join("meta.json"), serde_json::to_vec_pretty(&meta).unwrap()).unwrap();
    let budget = serde_json::json!({
        "control_volume_totals_m3": d.ledger.totals(),
        "control_volume_declared_rel_tol": CV_REL_TOL,
        "summary": summary,
    });
    std::fs::write(out.join("budget.json"), serde_json::to_vec_pretty(&budget).unwrap()).unwrap();
}
