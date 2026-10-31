//! SEED-TRENCH (lane TN, 2026-10-13): the `extras_v3` profile (`extras_v2`
//! + `subduction_trench = slab-flexure-trench-v1` + the passive-margin
//! shelf law, when present).
//!
//! The normal test runs the profile on lane OG's mini grid (36 x 72) through
//! the batch path, the stepper and the genesis module: the three agree, the
//! crust ledger closes, the assumption is recorded, the step log reports a
//! trench, and extras_v2 logs none.
//!
//! The ignored driver `trenches_trace` runs one (seed, fidelity, mode) and
//! writes the final fields plus the trench term for maps and
//! cross-sections. Env: TN_SEED, TN_FIDELITY (preview | standard), TN_MODE
//! (extras_v2 | extras_v3), TN_OUT, optional TN_STEP_MY.

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
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
        "extras_v3" => PreviewOptions::extras_v3(),
        other => panic!("unknown mode {other}"),
    }
}

#[test]
fn trench_profile_closes_replays_and_repeats() {
    let opts = PreviewOptions::extras_v3();
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
    for want in ["A-TRENCH", "A-RIFT-FLANK", "A-PASSIVE-MARGIN"] {
        assert!(ids.contains(&want), "{want} missing");
    }
    let og = a.orogen.as_ref().unwrap();
    let logs: Vec<_> = og.steps.iter().filter_map(|s| s.trench.clone()).collect();
    assert_eq!(logs.len(), og.steps.len());
    assert!(logs.iter().any(|l| l.trench_cells > 0 && l.max_w0_m > 0.0));
    // extras_v2 logs no trench and records no assumption.
    let v2 = simulate_with(&cfg, &PreviewOptions::extras_v2()).unwrap();
    assert!(v2.orogen.as_ref().unwrap().steps.iter().all(|s| s.trench.is_none()));
    assert!(!v2.terrain.assumptions.iter().any(|x| x.id == "A-TRENCH" || x.id == "A-PASSIVE-MARGIN"));
    // Seeded margins carry stretching markers under extras_v3 only.
    let beta = |o: &pf_lab::preview::PreviewOutputs| {
        o.orogen.as_ref().unwrap().extras.as_ref().unwrap().rift_flank.as_ref().unwrap().beta.iter().filter(|&&b| b > 1.0).count()
    };
    assert!(beta(&a) > beta(&v2), "{} {}", beta(&a), beta(&v2));
    // The trench deepens the ocean floor.
    let min = |e: &[f64]| e.iter().cloned().fold(f64::MAX, f64::min);
    assert!(min(&a.elevation_m) < min(&v2.elevation_m));
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
        "trench mini: trench cells {} flexed {} max w0 {:.0} m max bulge {:.0} m; min elevation {:.0} vs {:.0} m",
        last.trench_cells,
        last.flexed_cells,
        last.max_w0_m,
        last.max_bulge_wl_m,
        min(&a.elevation_m),
        min(&v2.elevation_m)
    );
}

#[test]
#[ignore]
fn trenches_trace() {
    let seed: u64 = std::env::var("TN_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("TN_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("TN_MODE").unwrap();
    let step_my: Option<f64> = std::env::var("TN_STEP_MY").ok().map(|s| s.parse().unwrap());
    let out_dir = PathBuf::from(std::env::var("TN_OUT").unwrap());
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
    let env = PreviewEnv::new(params.clone(), cfg.grid, cfg.seed_path()).unwrap();
    // A trench-law env evaluates the term for diagnostics in every mode.
    let probe = PreviewEnv::new(
        cfg.profile_with(&PreviewOptions::extras_v3()).unwrap(),
        cfg.grid,
        cfg.seed_path(),
    )
    .unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let grid = run.env().grid().clone();
    let areas = run.env().areas().to_vec();
    let total: f64 = areas.iter().sum();
    let mut csv = std::io::BufWriter::new(std::fs::File::create(out_dir.join("steps.csv")).unwrap());
    writeln!(
        csv,
        "t_my,sea_level_gauge_m,land_frac,min_elev_m,shelf_0_200_frac,below_6km_frac,trench_cells,max_w0_m,max_bulge_wl_m,deposit_aggraded_m3,exported_m3"
    )
    .unwrap();
    // Optional snapshots: TN_DUMP_AT = comma-separated times (My).
    let dump_at: Vec<f64> = std::env::var("TN_DUMP_AT")
        .map(|v| v.split(',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect())
        .unwrap_or_default();
    let dump = |run: &PreviewRun, t: f64| {
        let d = out_dir.join(format!("snap_{:.0}", t));
        std::fs::create_dir_all(&d).unwrap();
        write_f32(&d.join("elevation_m.f32"), run.elevation_m());
        write_f32(&d.join("sediment_m.f32"), &run.state().sediment_m);
        write_f32(&d.join("thickness_m.f32"), &run.state().thickness_m);
        let kind: Vec<u8> = run.state().kind.iter().map(|k| *k as u8).collect();
        std::fs::write(d.join("kind.u8"), kind).unwrap();
        let plate: Vec<f64> = run.state().plate.iter().map(|p| p.0 as f64).collect();
        write_f32(&d.join("plate.f32"), &plate);
    };
    if dump_at.contains(&0.0) {
        dump(&run, 0.0);
    }
    while run.steps_taken() < run.steps_scheduled() {
        let r = run.step().unwrap();
        if dump_at.iter().any(|&t| (t - r.t_my).abs() < 1e-9) {
            dump(&run, r.t_my);
        }
        let log = run.logs().last().unwrap().clone();
        let e = run.elevation_m();
        let (mut shelf, mut deep) = (0.0, 0.0);
        for c in 0..e.len() {
            if e[c] < 0.0 && e[c] >= -200.0 {
                shelf += areas[c];
            }
            if e[c] < -6000.0 {
                deep += areas[c];
            }
        }
        let (tc, w0, bu) = probe
            .trench_term_of(run.state())
            .unwrap()
            .map(|t| (t.1.trench_cells, t.1.max_w0_m, t.1.max_bulge_wl_m))
            .unwrap();
        writeln!(
            csv,
            "{},{:.2},{:.5},{:.1},{:.5},{:.6},{},{:.1},{:.1},{:.4e},{:.4e}",
            r.t_my,
            log.sea_level_gauge_m.unwrap_or(f64::NAN),
            log.land_fraction,
            log.min_elevation_m,
            shelf / total,
            deep / total,
            tc,
            w0,
            bu,
            r.deposit_aggraded_m3,
            r.exported_m3
        )
        .unwrap();
    }
    let st = run.state().clone();
    let term = probe.trench_term_of(&st).unwrap().unwrap();
    let tt = pf_geodynamics::trench::trench_term(
        &grid,
        &areas,
        &st,
        probe.trench_params().unwrap(),
    )
    .unwrap();
    let mut is_trench = vec![0u8; st.cells()];
    let mut w0 = vec![0.0f64; st.cells()];
    for &(c, w, _, _, _) in &tt.trench_cells {
        is_trench[c as usize] = 1;
        w0[c as usize] = w;
    }
    let out = run.finish().unwrap();
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    let kind: Vec<u8> = out.state.kind.iter().map(|k| *k as u8).collect();
    std::fs::write(out_dir.join("kind.u8"), kind).unwrap();
    std::fs::write(out_dir.join("boundary_class.u8"), &out.state.boundary_class).unwrap();
    std::fs::write(out_dir.join("trench_cell.u8"), is_trench).unwrap();
    let plate: Vec<f64> = out.state.plate.iter().map(|p| p.0 as f64).collect();
    write_f32(&out_dir.join("plate.f32"), &plate);
    let age: Vec<f64> = (0..out.state.cells()).map(|c| out.state.age_my(c)).collect();
    write_f32(&out_dir.join("age_my.f32"), &age);
    write_f32(&out_dir.join("normal_rate_m_per_my.f32"), &out.state.normal_rate_m_per_my);
    write_f32(&out_dir.join("trench_term_wl_m.f32"), &term.0);
    write_f32(&out_dir.join("trench_w0_m.f32"), &w0);
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
