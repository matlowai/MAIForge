//! SEED-RIFT-SHOULDERS (lane RS, 2026-10-13): the `extras_v2` profile
//! (`extras_v1` + `rift_flank = shoulder-uplift-thermal-v1`).
//!
//! The normal test runs the profile on lane OG's mini grid (36 x 72) through
//! the batch path, the stepper and the genesis module (whose control-volume
//! closures are checked every step): the three agree, the crust ledger
//! closes, the assumption is recorded, the markers are carried and the step
//! log reports a nonzero rift term.
//!
//! The ignored driver `rift_shoulders_trace` runs one (seed, fidelity, mode)
//! and writes the final fields, markers and the rift term parts for maps
//! and cross-sections. Env: RS_SEED, RS_FIDELITY (preview | standard),
//! RS_MODE (extras_v1 | extras_v2), RS_OUT, optional RS_STEP_MY.

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{ModuleRun, RunBinding};
use pf_geodynamics::rift_flank::{flank_term, RiftFlankParams};
use pf_lab::preview::{
    gdh1_thermal_elevation_m, simulate_with, PreviewConfig, PreviewEnv, PreviewOptions, PreviewRun,
};
use pf_lab::preview_module::{outputs_digest, preview_recipe, PreviewGenesisModule as M};
use std::io::Write;
use std::path::{Path, PathBuf};

fn write_f32(path: &Path, v: &[f64]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for x in v {
        f.write_all(&(*x as f32).to_le_bytes()).unwrap();
    }
}

fn options(mode: &str) -> PreviewOptions {
    match mode {
        "extras_v1" => PreviewOptions::extras_v1(),
        "extras_v2" => PreviewOptions::extras_v2(),
        other => panic!("unknown mode {other}"),
    }
}

#[test]
fn rift_flank_profile_closes_replays_and_repeats() {
    let opts = PreviewOptions::extras_v2();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed: 240907,
        duration_my: Some(30.0),
        step_my: Some(2.0),
    };
    let a = simulate_with(&cfg, &opts).unwrap();
    let da = outputs_digest(&a).unwrap();
    let b = simulate_with(&cfg, &opts).unwrap();
    assert_eq!(outputs_digest(&b).unwrap(), da);
    let cv = &a.terrain.crust_volume;
    assert!(
        cv.closure_plug_m3.abs() <= 1e-9 * cv.initial_m3,
        "closure plug {} of {}",
        cv.closure_plug_m3,
        cv.initial_m3
    );
    let ids: Vec<&str> = a.terrain.assumptions.iter().map(|x| x.id.as_str()).collect();
    for want in ["A-RIFT-FLANK", "A-CONTINENTAL-RIFT", "A-HOTSPOT-PLUMES"] {
        assert!(ids.contains(&want), "{want} missing");
    }
    let og = a.orogen.as_ref().unwrap();
    let m = og.extras.as_ref().unwrap().rift_flank.as_ref().unwrap();
    assert!(m.beta.iter().all(|&b| b >= 1.0));
    let logs: Vec<_> = og
        .steps
        .iter()
        .filter_map(|s| s.extras.as_ref().and_then(|x| x.rift_flank.clone()))
        .collect();
    assert_eq!(logs.len(), og.steps.len());
    assert!(logs.iter().any(|l| l.stretched_cells > 0 && l.max_thermal_wl_m > 0.0));
    // The extras_v1 profile carries no flank state and logs none.
    let v1 = simulate_with(&cfg, &PreviewOptions::extras_v1()).unwrap();
    let og1 = v1.orogen.as_ref().unwrap();
    assert!(og1.extras.as_ref().unwrap().rift_flank.is_none());
    assert!(og1.steps.iter().all(|s| s.extras.as_ref().unwrap().rift_flank.is_none()));
    assert!(!v1.terrain.assumptions.iter().any(|x| x.id == "A-RIFT-FLANK"));
    // Stepper and module replay equal the batch.
    let env = PreviewEnv::new(cfg.profile_with(&opts).unwrap(), cfg.grid, cfg.seed_path()).unwrap();
    let mut stepper = PreviewRun::init(env).unwrap();
    while stepper.steps_taken() < stepper.steps_scheduled() {
        stepper.step().unwrap();
    }
    let stepped = stepper.finish().unwrap();
    assert_eq!(outputs_digest(&stepped).unwrap(), da);
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut run = ModuleRun::<M>::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
    run.run_to_end(None).unwrap();
    let out = M::finish(run.into_state()).unwrap();
    assert_eq!(outputs_digest(&out).unwrap(), da, "module != batch");
    let last = logs.last().unwrap();
    eprintln!(
        "rift flank mini: stretched {} shoulder {} max thermal {:.0} m max shoulder {:.0} m (water-loaded)",
        last.stretched_cells, last.shoulder_cells, last.max_thermal_wl_m, last.max_shoulder_wl_m
    );
}

#[test]
#[ignore]
fn rift_shoulders_trace() {
    let seed: u64 = std::env::var("RS_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("RS_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("RS_MODE").unwrap();
    let step_my: Option<f64> = std::env::var("RS_STEP_MY").ok().map(|s| s.parse().unwrap());
    let out_dir = PathBuf::from(std::env::var("RS_OUT").unwrap());
    std::fs::create_dir_all(&out_dir).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: None,
        step_my,
    };
    let params = cfg.profile_with(&options(&mode)).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let grid = run.env().grid().clone();
    let areas = run.env().areas().to_vec();
    let mut rows: Vec<(f64, usize, u64, f64, f64)> = Vec::new();
    while run.steps_taken() < run.steps_scheduled() {
        let r = run.step().unwrap();
        let log = run.logs().last().unwrap();
        rows.push((
            r.t_my,
            r.advance.rift_births.len(),
            r.advance.rift_breakup_births,
            log.max_elevation_m,
            log.land_fraction,
        ));
    }
    let out = run.finish().unwrap();
    let mut csv = std::io::BufWriter::new(std::fs::File::create(out_dir.join("steps.csv")).unwrap());
    writeln!(
        csv,
        "t_my,rift_births,breakup_births,stretched_cells,shoulder_cells,max_thermal_wl_m,max_shoulder_wl_m,max_elev_m,land_frac"
    )
    .unwrap();
    let og_steps = &out.orogen.as_ref().unwrap().steps;
    for (k, (t, rb, bb, mx, lf)) in rows.iter().enumerate() {
        let (sc, shc, mt, ms) = og_steps[k]
            .extras
            .as_ref()
            .and_then(|x| x.rift_flank.clone())
            .map(|f| (f.stretched_cells, f.shoulder_cells, f.max_thermal_wl_m, f.max_shoulder_wl_m))
            .unwrap_or((0, 0, 0.0, 0.0));
        writeln!(csv, "{t},{rb},{bb},{sc},{shc},{mt:.1},{ms:.1},{mx:.1},{lf:.5}").unwrap();
    }
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    let kind: Vec<u8> = out.state.kind.iter().map(|k| *k as u8).collect();
    std::fs::write(out_dir.join("kind.u8"), kind).unwrap();
    let plate: Vec<f64> = out.state.plate.iter().map(|p| p.0 as f64).collect();
    write_f32(&out_dir.join("plate.f32"), &plate);
    if let Some(x) = out.orogen.as_ref().and_then(|o| o.extras.as_ref()) {
        write_f32(&out_dir.join("rift_t_my.f32"), &x.rift_t_my);
        write_f32(&out_dir.join("rift_thinned_m.f32"), &x.rift_thinned_m);
        write_f32(&out_dir.join("volcanic_added_m.f32"), &x.volcanic_added_m);
        if let Some(m) = &x.rift_flank {
            write_f32(&out_dir.join("rift_beta.f32"), &m.beta);
            write_f32(&out_dir.join("rift_thermal_t0_my.f32"), &m.thermal_t0_my);
            let p = RiftFlankParams {
                flexural_parameter_m: 60_000.0,
                shoulder_ratio: 1.0,
            };
            let t = flank_term(
                &grid,
                &areas,
                &out.state.kind,
                &m.beta,
                &m.thermal_t0_my,
                out.state.t_my,
                &p,
                &gdh1_thermal_elevation_m,
            )
            .unwrap();
            write_f32(&out_dir.join("rift_thermal_wl_m.f32"), &t.thermal_wl_m);
            write_f32(&out_dir.join("rift_shoulder_wl_m.f32"), &t.shoulder_wl_m);
        }
    }
    let cv = &out.terrain.crust_volume;
    let meta = serde_json::json!({
        "seed": seed, "mode": mode, "nlat": grid.nlat(), "nlon": grid.nlon(),
        "t_my": out.state.t_my, "step_my": step_my,
        "row_area_m2": (0..grid.nlat()).map(|i| areas[i * grid.nlon()]).collect::<Vec<_>>(),
        "closure_plug_m3": cv.closure_plug_m3, "initial_m3": cv.initial_m3,
        "land_area_fraction": out.terrain.summary.land_area_fraction,
        "max_elevation_m": out.terrain.summary.max_elevation_m,
        "outputs_digest": outputs_digest(&out).unwrap(),
    });
    std::fs::write(out_dir.join("meta.json"), serde_json::to_string_pretty(&meta).unwrap()).unwrap();
}
