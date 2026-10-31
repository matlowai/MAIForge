//! SEED-UPLANDS (lane UP, 2026-10-13): the `uplands_v1` profile
//! (`extras_v2` + `stream-power-subgrid-nonlinear-v1` +
//! `critical-slope-envelope-v1`) and labelled variants. A quick mini-grid
//! test (deterministic, closes, stepper = batch = module) and an ignored
//! trace driver. Env: UP_SEED, UP_FIDELITY (preview | standard), UP_MODE
//! (extras_v2 | uplands_v1 | env_n1 | n2_only), UP_OUT, optional
//! UP_DURATION_MY, UP_STEP_MY, UP_SNAP_EVERY_MY (elevation + thickness
//! snapshots for the aging analysis).

use pf_contracts::params::{
    CriticalSlopeErosionBlock, FidelityPreset, NonlinearSubgridErosionBlock, OrogenyBlock,
    ScenarioPreset,
};
use pf_contracts::sim_module::{ModuleRun, RunBinding};
use pf_lab::preview::{simulate_with, PreviewConfig, PreviewEnv, PreviewOptions, PreviewRun};
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
        "extras_v2" => PreviewOptions::extras_v2(),
        "uplands_v1" => PreviewOptions::uplands_v1(),
        // Labelled variants (same constants): the envelope alone on OG's
        // n = 1 law, and RL's n = 2 law alone, both on extras_v2.
        "env_n1" => PreviewOptions {
            orogeny: Some(OrogenyBlock {
                critical_slope_erosion: Some(CriticalSlopeErosionBlock::v1()),
                ..OrogenyBlock::extras_v2()
            }),
            ..PreviewOptions::extras_v2()
        },
        "n2_only" => PreviewOptions {
            orogeny: Some(OrogenyBlock {
                subgrid_erosion: None,
                nonlinear_subgrid_erosion: Some(NonlinearSubgridErosionBlock::v1()),
                ..OrogenyBlock::extras_v2()
            }),
            ..PreviewOptions::extras_v2()
        },
        // LABELLED SENSITIVITY (lane UP, not a proposal): uplands_v1 with
        // RL's reference steepness halved (K_n doubled at n = 2).
        "uplands_kref250" => {
            let mut b = NonlinearSubgridErosionBlock::v1();
            b.reference_steepness_m = 250.0;
            PreviewOptions {
                orogeny: Some(OrogenyBlock {
                    nonlinear_subgrid_erosion: Some(b),
                    ..OrogenyBlock::uplands_v1()
                }),
                ..PreviewOptions::extras_v2()
            }
        }
        other => panic!("unknown mode {other}"),
    }
}

#[test]
fn uplands_profile_closes_replays_and_repeats() {
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed: 240907,
        duration_my: Some(30.0),
        step_my: Some(2.0),
    };
    for mode in ["uplands_v1", "env_n1"] {
        let opts = options(mode);
        let a = simulate_with(&cfg, &opts).unwrap();
        let da = outputs_digest(&a).unwrap();
        let b = simulate_with(&cfg, &opts).unwrap();
        assert_eq!(outputs_digest(&b).unwrap(), da);
        let cv = &a.terrain.crust_volume;
        assert!(
            cv.closure_plug_m3.abs() <= 1e-9 * cv.initial_m3,
            "{mode}: closure plug {} of {}",
            cv.closure_plug_m3,
            cv.initial_m3
        );
        let ids: Vec<&str> = a.terrain.assumptions.iter().map(|x| x.id.as_str()).collect();
        for want in ["A-HILLSLOPE-CRITICAL-SLOPE", "A-KNICKPOINT-CELERITY", "A-RIFT-FLANK"] {
            assert!(ids.contains(&want), "{mode}: {want} missing");
        }
        // The envelope must change the result against its base profile.
        let base = if mode == "env_n1" { options("extras_v2") } else { options("n2_only") };
        let o = simulate_with(&cfg, &base).unwrap();
        assert_ne!(outputs_digest(&o).unwrap(), da, "{mode}: envelope had no effect");
        assert!(!o.terrain.assumptions.iter().any(|x| x.id == "A-HILLSLOPE-CRITICAL-SLOPE"));
        // Stepper and module replay equal the batch.
        let env = PreviewEnv::new(cfg.profile_with(&opts).unwrap(), cfg.grid, cfg.seed_path()).unwrap();
        let mut stepper = PreviewRun::init(env).unwrap();
        while stepper.steps_taken() < stepper.steps_scheduled() {
            stepper.step().unwrap();
        }
        assert_eq!(outputs_digest(&stepper.finish().unwrap()).unwrap(), da);
        let recipe = preview_recipe(&cfg, &opts, false).unwrap();
        let mut run = ModuleRun::<M>::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
        run.run_to_end(None).unwrap();
        let out = M::finish(run.into_state()).unwrap();
        assert_eq!(outputs_digest(&out).unwrap(), da, "{mode}: module != batch");
    }
}

#[test]
#[ignore]
fn uplands_trace() {
    let seed: u64 = std::env::var("UP_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("UP_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("UP_MODE").unwrap();
    let num = |k: &str| std::env::var(k).ok().map(|s| s.parse::<f64>().unwrap());
    let out_dir = PathBuf::from(std::env::var("UP_OUT").unwrap());
    std::fs::create_dir_all(&out_dir).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: num("UP_DURATION_MY"),
        step_my: num("UP_STEP_MY"),
    };
    let snap = num("UP_SNAP_EVERY_MY");
    let env = PreviewEnv::new(cfg.profile_with(&options(&mode)).unwrap(), cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let grid = run.env().grid().clone();
    // Start state for the continental ledger (integrator request 2026-10-13).
    write_f32(&out_dir.join("thick_t0000.f32"), &run.state().thickness_m);
    write_f32(&out_dir.join("sed_t0000.f32"), &run.state().sediment_m);
    write_f32(&out_dir.join("elev_t0000.f32"), run.elevation_m());
    let kind0: Vec<u8> = run.state().kind.iter().map(|k| *k as u8).collect();
    std::fs::write(out_dir.join("kind_t0000.u8"), kind0).unwrap();
    let mut csv = String::from("t_my,max_elev_m,land_frac\n");
    while run.steps_taken() < run.steps_scheduled() {
        let r = run.step().unwrap();
        let log = run.logs().last().unwrap();
        csv += &format!("{},{:.1},{:.5}\n", r.t_my, log.max_elevation_m, log.land_fraction);
        if let Some(every) = snap {
            let k = (r.t_my / every).round();
            if (r.t_my - k * every).abs() < 1e-9 {
                write_f32(&out_dir.join(format!("elev_t{:04}.f32", r.t_my as i64)), run.elevation_m());
                write_f32(&out_dir.join(format!("thick_t{:04}.f32", r.t_my as i64)), &run.state().thickness_m);
                let plate: Vec<f64> = run.state().plate.iter().map(|p| p.0 as f64).collect();
                write_f32(&out_dir.join(format!("plate_t{:04}.f32", r.t_my as i64)), &plate);
            }
        }
    }
    std::fs::write(out_dir.join("steps.csv"), csv).unwrap();
    let out = run.finish().unwrap();
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    let kind: Vec<u8> = out.state.kind.iter().map(|k| *k as u8).collect();
    std::fs::write(out_dir.join("kind.u8"), kind).unwrap();
    if let Some(r) = out.orogen.as_ref().and_then(|o| o.subgrid_relief_m.as_ref()) {
        write_f32(&out_dir.join("subgrid_relief_m.f32"), r);
    }
    let cv = &out.terrain.crust_volume;
    let meta = serde_json::json!({
        "seed": seed, "mode": mode, "nlat": grid.nlat(), "nlon": grid.nlon(),
        "t_my": out.state.t_my, "duration_my": cfg.duration_my, "step_my": cfg.step_my,
        "closure_plug_m3": cv.closure_plug_m3, "initial_m3": cv.initial_m3,
        "eroded_m3": cv.eroded_m3,
        "land_area_fraction": out.terrain.summary.land_area_fraction,
        "max_elevation_m": out.terrain.summary.max_elevation_m,
        "outputs_digest": outputs_digest(&out).unwrap(),
    });
    std::fs::write(out_dir.join("meta.json"), serde_json::to_string_pretty(&meta).unwrap()).unwrap();
}
