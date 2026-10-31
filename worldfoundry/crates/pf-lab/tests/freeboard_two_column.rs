//! FREEBOARD-SPIKE (lane F, 2026-10-11): drive the preview chain with the
//! OPT-IN two-column isostasy + volume-primary sea level, and measure it
//! against the legacy chain.
//!
//! Why a driver lives in a test: `pf-lab/src/preview.rs` is not in lane F's
//! owned files, so the flag is not wired into `simulate()` yet (with the
//! flag set, `simulate()` refuses loudly via `isostatic_elevation_m`). This
//! file re-runs the same per-step chain (advance -> isostasy -> ocean ->
//! route -> erode -> isostasy) using only public APIs. The guard
//! `driver_with_flag_off_is_bitwise_simulate` proves the re-run equals
//! `simulate()` exactly when the flag is off, so the only difference in a
//! flag-on run is the isostasy/sea-level step.
//!
//! The ignored test `measure_seed` writes raw fields for the lane F
//! measurement; it is env-driven (FREEBOARD_SEED, FREEBOARD_FIDELITY,
//! FREEBOARD_OUT, FREEBOARD_CHECK_OFF) and is not part of the normal suite.

use pf_contracts::params::{PlanetParameterProfile, TwoColumnIsostasyBlock};
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::isostasy::{two_column_field, ColumnV0, TwoColumnSolution};
use pf_geodynamics::kinematic::{advance, classify_cells, CrustState};
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;
use pf_lab::preview::{
    elevation_from_state, gdh1_thermal_elevation_m, simulate, PreviewConfig,
    DIFFUSION_STABILITY_LIMIT, STREAM_POWER_M, STREAM_POWER_N, YEARS_PER_MY,
};
use pf_surface::erosion::{erode_step, DiffusionStencil, ErosionParams};
use pf_surface::routing::route;
use pf_surface::water::largest_ocean;
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, serde::Serialize)]
struct DriveStep {
    t_my: f64,
    sea_level_gauge_m: Option<f64>,
    flooded_area_fraction: Option<f64>,
    land_area_fraction: f64,
    eroded_m3: f64,
    deposited_m3: f64,
    exported_m3: f64,
    max_elevation_m: f64,
    min_elevation_m: f64,
}

struct Drive {
    grid: GridSpec,
    areas: Vec<f64>,
    elevation_t0: Vec<f64>,
    elevation: Vec<f64>,
    state: CrustState,
    steps: Vec<DriveStep>,
    t0_solution: Option<TwoColumnSolution>,
    final_solution: Option<TwoColumnSolution>,
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

fn columns(state: &CrustState) -> Vec<ColumnV0> {
    (0..state.cells())
        .map(|c| ColumnV0 {
            kind: state.kind[c],
            crust_m: state.thickness_m[c],
            sediment_m: state.sediment_m[c],
            thermal_uplift_water_loaded_m: if state.kind[c] == CrustKindCodeV0::Oceanic {
                gdh1_thermal_elevation_m(state.age_my(c))
            } else {
                0.0
            },
        })
        .collect()
}

/// Elevation by the selected isostasy: legacy (`elevation_from_state`, the
/// exact preview.rs function) when the flag is absent, two-column otherwise.
fn elevation(
    state: &CrustState,
    params: &PlanetParameterProfile,
    areas: &[f64],
) -> (Vec<f64>, Option<TwoColumnSolution>) {
    if params.two_column_isostasy.is_some() {
        let sol = two_column_field(&columns(state), areas, params).unwrap();
        (sol.elevation_m.clone(), Some(sol))
    } else {
        (elevation_from_state(state, params).unwrap(), None)
    }
}

fn drive(config: &PreviewConfig, flag: Option<TwoColumnIsostasyBlock>) -> Drive {
    drive_with(config, flag, false)
}

/// `no_erosion` is a DIAGNOSTIC counterfactual only (stream-power K and
/// hillslope D set to 0) used to separate isostasy from the erosion budget
/// (D4); it is never a proposed setting.
fn drive_with(config: &PreviewConfig, flag: Option<TwoColumnIsostasyBlock>, no_erosion: bool) -> Drive {
    let mut params = config.profile().unwrap();
    params.two_column_isostasy = flag;
    if no_erosion {
        params.surface.erodibility = 0.0;
        params.surface.hillslope_diffusivity_m2_per_yr = 0.0;
    }
    params.validate().unwrap();
    let grid = config.grid(&params).unwrap();
    let seed_path = config.seed_path();
    let prior = UpliftPrior::preview_v0();
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

    let seeding = seed_plates(&grid, &params, &seed_path).unwrap();
    let plates = seeding.plates.clone();
    let mut state = CrustState::seed(&grid, &seeding, &prior).unwrap();
    let (class, rate, _) = classify_cells(&grid, &state, &plates).unwrap();
    state.boundary_class = class;
    state.normal_rate_m_per_my = rate;
    let (elevation_t0, t0_solution) = elevation(&state, &params, &areas);
    let mut elev = elevation_t0.clone();
    let mut final_solution = t0_solution.clone();

    let dt = params.tectonics.step_my;
    let steps = (params.tectonics.duration_my / dt).round() as u32;
    let rain = vec![params.surface.uniform_runoff_m_per_yr * YEARS_PER_MY; n];
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let mut log = Vec::new();
    for _ in 0..steps {
        let (next, _ledger, _thickening) = advance(&grid, &state, &plates, &prior, dt).unwrap();
        let (elev_tect, _) = elevation(&next, &params, &areas);
        state = next;
        let ocean = largest_ocean(&grid, &elev_tect, sea).unwrap();
        let routing = route(&grid, &elev_tect, &rain, &ocean, sea).unwrap();
        let step = erode_step(
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
        .unwrap();
        for c in 0..n {
            state.thickness_m[c] -= step.bedrock_removed_m[c];
            state.eroded_total_m[c] += step.bedrock_removed_m[c];
            state.sediment_m[c] += step.sediment_added_m[c] - step.sediment_removed_m[c];
            if state.sediment_m[c] < 0.0 {
                assert!(state.sediment_m[c] >= -1.0e-9, "negative sediment at {c}");
                state.sediment_m[c] = 0.0;
            }
        }
        let (e, sol) = elevation(&state, &params, &areas);
        elev = e;
        let (mut mx, mut mn, mut land) = (f64::MIN, f64::MAX, 0.0);
        for c in 0..n {
            mx = mx.max(elev[c]);
            mn = mn.min(elev[c]);
            if elev[c] >= sea {
                land += areas[c];
            }
        }
        log.push(DriveStep {
            t_my: state.t_my,
            sea_level_gauge_m: sol.as_ref().map(|s| s.sea_level_gauge_m),
            flooded_area_fraction: sol.as_ref().map(|s| s.flooded_area_fraction),
            land_area_fraction: land / total_area,
            eroded_m3: step.eroded_m3,
            deposited_m3: step.deposited_m3,
            exported_m3: step.exported_m3,
            max_elevation_m: mx,
            min_elevation_m: mn,
        });
        final_solution = sol;
    }
    Drive {
        grid,
        areas,
        elevation_t0,
        elevation: elev,
        state,
        steps: log,
        t0_solution,
        final_solution,
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
fn driver_with_flag_off_is_bitwise_simulate() {
    for seed in [240907u64, 7] {
        let cfg = mini(seed);
        let sim = simulate(&cfg).unwrap();
        let d = drive(&cfg, None);
        assert_eq!(d.elevation, sim.elevation_m, "seed {seed}: driver diverged");
        assert_eq!(d.state.thickness_m, sim.state.thickness_m);
        assert_eq!(d.state.sediment_m, sim.state.sediment_m);
        for (a, b) in d.steps.iter().zip(&sim.steps) {
            assert_eq!(a.land_area_fraction, b.land_fraction);
            assert_eq!(a.eroded_m3, b.eroded_m3);
        }
    }
}

#[test]
fn flag_on_is_deterministic_conserves_water_and_differs_from_legacy() {
    let cfg = mini(240907);
    let a = drive(&cfg, Some(TwoColumnIsostasyBlock::earth_v0()));
    let b = drive(&cfg, Some(TwoColumnIsostasyBlock::earth_v0()));
    assert_eq!(a.elevation, b.elevation, "flag-on run not deterministic");
    let off = drive(&cfg, None);
    assert_ne!(a.elevation, off.elevation, "flag had no effect");
    let v = TwoColumnIsostasyBlock::earth_v0().ocean_water_volume_m3;
    for sol in [a.t0_solution.as_ref().unwrap(), a.final_solution.as_ref().unwrap()] {
        assert!(((sol.volume_check_m3 - v) / v).abs() < 1e-9, "water not conserved");
    }
}

#[test]
fn simulate_refuses_the_flag_instead_of_ignoring_it() {
    // simulate() builds its profile from the preset, so the flag cannot be
    // injected there; but the per-cell function it calls refuses it.
    let cfg = mini(240907);
    let mut p = cfg.profile().unwrap();
    p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
    let grid = cfg.grid(&p).unwrap();
    let seeding = seed_plates(&grid, &p, &cfg.seed_path()).unwrap();
    let state = CrustState::seed(&grid, &seeding, &UpliftPrior::preview_v0()).unwrap();
    assert!(elevation_from_state(&state, &p).is_err());
}

// ---------------------------------------------------------------------------
// Measurement (ignored; lane F runs it explicitly per seed).
// ---------------------------------------------------------------------------

fn write_f32(path: &Path, v: &[f64]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for x in v {
        f.write_all(&(*x as f32).to_le_bytes()).unwrap();
    }
}

fn write_u8(path: &Path, v: &[u8]) {
    std::fs::write(path, v).unwrap();
}

fn dump(dir: &Path, d: &Drive, sea: f64, extra: serde_json::Value) {
    std::fs::create_dir_all(dir).unwrap();
    write_f32(&dir.join("elevation_t0_m.f32"), &d.elevation_t0);
    write_f32(&dir.join("elevation_m.f32"), &d.elevation);
    write_f32(&dir.join("thickness_m.f32"), &d.state.thickness_m);
    write_f32(&dir.join("sediment_m.f32"), &d.state.sediment_m);
    write_f32(&dir.join("eroded_total_m.f32"), &d.state.eroded_total_m);
    let kind: Vec<u8> = d.state.kind.iter().map(|k| k.code()).collect();
    write_u8(&dir.join("kind.u8"), &kind);
    let ocean = largest_ocean(&d.grid, &d.elevation, sea).unwrap();
    write_u8(
        &dir.join("ocean_connected.u8"),
        &ocean.iter().map(|&b| b as u8).collect::<Vec<_>>(),
    );
    let ocean0 = largest_ocean(&d.grid, &d.elevation_t0, sea).unwrap();
    write_u8(
        &dir.join("ocean_connected_t0.u8"),
        &ocean0.iter().map(|&b| b as u8).collect::<Vec<_>>(),
    );
    let row_areas: Vec<f64> = (0..d.grid.nlat()).map(|i| d.grid.cell_area_m2(i).unwrap()).collect();
    let meta = serde_json::json!({
        "nlat": d.grid.nlat(),
        "nlon": d.grid.nlon(),
        "row_area_m2": row_areas,
        "total_area_m2": d.areas.iter().sum::<f64>(),
        "sea_level_datum_m": sea,
        "steps": d.steps,
        "t0_sea_level_gauge_m": d.t0_solution.as_ref().map(|s| s.sea_level_gauge_m),
        "t0_flooded_area_fraction": d.t0_solution.as_ref().map(|s| s.flooded_area_fraction),
        "t0_volume_check_m3": d.t0_solution.as_ref().map(|s| s.volume_check_m3),
        "final_sea_level_gauge_m": d.final_solution.as_ref().map(|s| s.sea_level_gauge_m),
        "final_flooded_area_fraction": d.final_solution.as_ref().map(|s| s.flooded_area_fraction),
        "final_volume_check_m3": d.final_solution.as_ref().map(|s| s.volume_check_m3),
        "extra": extra,
    });
    std::fs::write(dir.join("meta.json"), serde_json::to_vec_pretty(&meta).unwrap()).unwrap();
}

#[test]
#[ignore]
fn measure_seed() {
    let seed: u64 = std::env::var("FREEBOARD_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("FREEBOARD_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let out = PathBuf::from(std::env::var("FREEBOARD_OUT").unwrap());
    let check_off = std::env::var("FREEBOARD_CHECK_OFF").map(|v| v == "1").unwrap_or(false);
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: None,
        step_my: None,
    };
    let sea = cfg.profile().unwrap().planet.sea_level_m;
    if std::env::var("FREEBOARD_NO_EROSION").map(|v| v == "1").unwrap_or(false) {
        for (name, flag) in [("off", None), ("on", Some(TwoColumnIsostasyBlock::earth_v0()))] {
            let t = std::time::Instant::now();
            let r = std::panic::catch_unwind(|| drive_with(&cfg, flag.clone(), true));
            match r {
                Ok(d) => dump(&out.join(name), &d, sea, serde_json::json!({
                    "counterfactual": "no erosion: surface.erodibility = 0 and hillslope diffusivity = 0 (diagnostic only)",
                    "drive_seconds": t.elapsed().as_secs_f64()})),
                Err(_) => {
                    std::fs::create_dir_all(out.join(name)).unwrap();
                    std::fs::write(out.join(name).join("ERROR.txt"), "drive panicked; see test stderr").unwrap();
                }
            }
        }
        return;
    }

    // OFF: the real simulate() is the baseline; the driver's own flag-off
    // run is compared bitwise when FREEBOARD_CHECK_OFF=1.
    let t = std::time::Instant::now();
    let sim = simulate(&cfg);
    let sim_s = t.elapsed().as_secs_f64();
    match &sim {
        Ok(sim) => {
            let mut off = drive_from_sim(&cfg, sim);
            let mut extra = serde_json::json!({"simulate_seconds": sim_s, "source": "pf_lab::preview::simulate"});
            if check_off {
                let d = drive(&cfg, None);
                let same = d.elevation == sim.elevation_m;
                extra["driver_off_bitwise_equal_to_simulate"] = serde_json::json!(same);
                assert!(same, "driver diverged from simulate at seed {seed}");
                off.steps = d.steps;
            }
            dump(&out.join("off"), &off, sea, extra);
        }
        Err(e) => {
            std::fs::create_dir_all(out.join("off")).unwrap();
            std::fs::write(out.join("off/ERROR.txt"), format!("{e}")).unwrap();
        }
    }
    // ON.
    let t = std::time::Instant::now();
    let on = std::panic::catch_unwind(|| drive(&cfg, Some(TwoColumnIsostasyBlock::earth_v0())));
    let on_s = t.elapsed().as_secs_f64();
    match on {
        Ok(d) => dump(
            &out.join("on"),
            &d,
            sea,
            serde_json::json!({"drive_seconds": on_s, "block": TwoColumnIsostasyBlock::earth_v0(), "source": "freeboard_two_column::drive"}),
        ),
        Err(_) => {
            std::fs::create_dir_all(out.join("on")).unwrap();
            std::fs::write(out.join("on/ERROR.txt"), "flag-on drive panicked; see test stderr").unwrap();
        }
    }
}

/// Package a real simulate() result in the Drive shape; t0 comes from the
/// identical seeding + preview.rs elevation_from_state.
fn drive_from_sim(cfg: &PreviewConfig, sim: &pf_lab::preview::PreviewOutputs) -> Drive {
    let params = &sim.params;
    let grid = sim.grid.clone();
    let seeding = seed_plates(&grid, params, &cfg.seed_path()).unwrap();
    let state0 = CrustState::seed(&grid, &seeding, &UpliftPrior::preview_v0()).unwrap();
    let elevation_t0 = elevation_from_state(&state0, params).unwrap();
    let areas = area_weights(&grid);
    let steps = sim
        .steps
        .iter()
        .map(|s| DriveStep {
            t_my: s.t_my,
            sea_level_gauge_m: None,
            flooded_area_fraction: None,
            land_area_fraction: s.land_fraction,
            eroded_m3: s.eroded_m3,
            deposited_m3: s.deposited_m3,
            exported_m3: s.exported_m3,
            max_elevation_m: s.max_elevation_m,
            min_elevation_m: s.min_elevation_m,
        })
        .collect();
    Drive {
        grid,
        areas,
        elevation_t0,
        elevation: sim.elevation_m.clone(),
        state: sim.state.clone(),
        steps,
        t0_solution: None,
        final_solution: None,
    }
}
