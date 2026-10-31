//! SEED-SEDIMENT-RETAIN (lane SR, 2026-10-13): the `retain_v1` profile
//! (`extras_v3` + `continental-retention-base-level-v1`). A quick mini-grid
//! test (opt-in law set, deterministic, closes, stepper = batch = module)
//! and an ignored trace driver. Env: SR_SEED, SR_FIDELITY (preview |
//! standard), SR_MODE (extras_v3 | retain_v1), SR_OUT, optional
//! SR_DURATION_MY, SR_STEP_MY, SR_SNAP_EVERY_MY (elevation, thickness,
//! sediment, kind snapshots and the retention partition per interval).

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{ModuleRun, RunBinding};
use pf_lab::preview::{simulate_with, PreviewConfig, PreviewEnv, PreviewOptions, PreviewRun};
use pf_lab::preview_module::{outputs_digest, preview_recipe, PreviewGenesisModule as M};
use pf_surface::erosion::retain::{take_thread_totals, RetentionReport};
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
        "extras_v3" => PreviewOptions::extras_v3(),
        "retain_v1" => PreviewOptions::retain_v1(),
        // SENSITIVITY ONLY (labelled; never a proposal): H = 500 m.
        "retain_h500_sens" => {
            let mut o = pf_contracts::params::OrogenyBlock::retain_v1();
            if let Some(b) = o.sediment_retention.as_mut() {
                b.base_level_relief_scale_m = 500.0;
            }
            PreviewOptions { orogeny: Some(o), ..PreviewOptions::retain_v1() }
        }
        other => panic!("unknown mode {other}"),
    }
}

fn report_json(r: &RetentionReport) -> serde_json::Value {
    serde_json::json!({
        "continental_grade_m3": r.continental_grade_m3,
        "lake_m3": r.lake_m3,
        "land_floor_m3": r.land_floor_m3,
        "continental_shelf_m3": r.continental_shelf_m3,
        "legacy_marine_m3": r.legacy_marine_m3,
        "global_layer_m3": r.global_layer_m3,
        "hillslope_to_water_m3": r.hillslope_to_water_m3,
        "added_on_continental_m3": r.added_on_continental_m3,
        "added_on_oceanic_m3": r.added_on_oceanic_m3,
        "taper_withheld_m3": r.taper_withheld_m3,
        "tapered_cells": r.tapered_cells,
        "mouths": r.mouths,
        "mouths_without_continental_shelf": r.mouths_without_continental_shelf,
        "shelf_search_bound_hits": r.shelf_search_bound_hits,
    })
}

#[test]
fn retain_profile_is_opt_in_closes_replays_and_repeats() {
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed: 240907,
        duration_my: Some(30.0),
        step_my: Some(2.0),
    };
    // Opt-in: absent from extras_v3's law set, present in retain_v1's.
    let base = cfg.profile_with(&options("extras_v3")).unwrap();
    let ls = serde_json::to_string(&base.universe().unwrap().law_set).unwrap();
    assert!(!ls.contains("sediment_retention"), "{ls}");
    let opts = options("retain_v1");
    let p = cfg.profile_with(&opts).unwrap();
    let u = p.universe().unwrap();
    u.validate().unwrap();
    let js = serde_json::to_string(&u.law_set).unwrap();
    assert!(js.contains("continental-retention-base-level-v1"), "{js}");
    let back: pf_contracts::laws::LawSetV1 = serde_json::from_str(&js).unwrap();
    assert_eq!(back, u.law_set);

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
    for want in [
        "A-BASE-LEVEL-RELIEF-SCALE",
        "A-BASE-LEVEL-TAPER-EXPONENT",
        "A-ALLUVIAL-GRADE-SLOPE",
        "A-CONTINENTAL-SHELF-WAVE-BASE",
    ] {
        assert!(ids.contains(&want), "{want} missing");
    }
    let o = simulate_with(&cfg, &options("extras_v3")).unwrap();
    assert_ne!(outputs_digest(&o).unwrap(), da, "retention had no effect");
    assert!(!o.terrain.assumptions.iter().any(|x| x.id == "A-BASE-LEVEL-RELIEF-SCALE"));
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
    assert_eq!(outputs_digest(&out).unwrap(), da, "module != batch");
}

#[test]
#[ignore]
fn sediment_retain_trace() {
    let seed: u64 = std::env::var("SR_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("SR_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("SR_MODE").unwrap();
    let num = |k: &str| std::env::var(k).ok().map(|s| s.parse::<f64>().unwrap());
    let out_dir = PathBuf::from(std::env::var("SR_OUT").unwrap());
    std::fs::create_dir_all(&out_dir).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: num("SR_DURATION_MY"),
        step_my: num("SR_STEP_MY"),
    };
    let snap = num("SR_SNAP_EVERY_MY");
    let env = PreviewEnv::new(cfg.profile_with(&options(&mode)).unwrap(), cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let grid = run.env().grid().clone();
    let dump = |run: &PreviewRun, tag: &str| {
        write_f32(&out_dir.join(format!("elev_{tag}.f32")), run.elevation_m());
        write_f32(&out_dir.join(format!("thick_{tag}.f32")), &run.state().thickness_m);
        write_f32(&out_dir.join(format!("sed_{tag}.f32")), &run.state().sediment_m);
        let kind: Vec<u8> = run.state().kind.iter().map(|k| *k as u8).collect();
        std::fs::write(out_dir.join(format!("kind_{tag}.u8")), kind).unwrap();
    };
    dump(&run, "t0000");
    let _ = take_thread_totals();
    let mut total = RetentionReport::default();
    let mut intervals = Vec::new();
    let mut csv = String::from("t_my,max_elev_m,land_frac,sea_level_gauge_m,eroded_m3,exported_m3,deposited_m3\n");
    while run.steps_taken() < run.steps_scheduled() {
        let r = run.step().unwrap();
        let log = run.logs().last().unwrap();
        csv += &format!(
            "{},{:.1},{:.5},{:.2},{:.6e},{:.6e},{:.6e}\n",
            r.t_my,
            log.max_elevation_m,
            log.land_fraction,
            log.sea_level_gauge_m.unwrap_or(f64::NAN),
            log.eroded_m3,
            log.exported_m3,
            log.deposited_m3
        );
        if let Some(every) = snap {
            let k = (r.t_my / every).round();
            if (r.t_my - k * every).abs() < 1e-9 {
                dump(&run, &format!("t{:04}", r.t_my as i64));
                let iv = take_thread_totals();
                total.accumulate(&iv);
                intervals.push(serde_json::json!({"t_my": r.t_my, "retention": report_json(&iv)}));
            }
        }
    }
    total.accumulate(&take_thread_totals());
    std::fs::write(out_dir.join("steps.csv"), csv).unwrap();
    let out = run.finish().unwrap();
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    let kind: Vec<u8> = out.state.kind.iter().map(|k| *k as u8).collect();
    std::fs::write(out_dir.join("kind.u8"), kind).unwrap();
    let cv = &out.terrain.crust_volume;
    let meta = serde_json::json!({
        "seed": seed, "mode": mode, "nlat": grid.nlat(), "nlon": grid.nlon(),
        "t_my": out.state.t_my, "duration_my": cfg.duration_my, "step_my": cfg.step_my,
        "closure_plug_m3": cv.closure_plug_m3, "initial_m3": cv.initial_m3,
        "final_m3": cv.final_m3, "eroded_m3": cv.eroded_m3, "subducted_m3": cv.subducted_m3,
        "land_area_fraction": out.terrain.summary.land_area_fraction,
        "max_elevation_m": out.terrain.summary.max_elevation_m,
        "sea_level_m": out.terrain.sea_level_m,
        "retention_total": report_json(&total),
        "retention_intervals": intervals,
        "outputs_digest": outputs_digest(&out).unwrap(),
    });
    std::fs::write(out_dir.join("meta.json"), serde_json::to_string_pretty(&meta).unwrap()).unwrap();
}
