//! SEED-EXTRAS (lane EX, 2026-10-12): the `extras_v1` profile
//! (`orogen_carry_v2` + `fixed-plume-chains-v1` + `stretched-rift-v1` +
//! `plume-head-continental-split-v1`).
//!
//! The normal test runs the profile on lane OG's mini grid (36 x 72) through
//! the batch path, the stepper and the genesis module (whose control-volume
//! closures are checked every step): three agree, the crust ledger closes,
//! the assumptions are recorded, and plume crust is added.
//!
//! The ignored driver `seed_extras_trace` runs one (seed, fidelity, mode)
//! and writes a per-step CSV plus the final fields and markers for maps.
//! Env: EX_SEED, EX_FIDELITY (preview | standard), EX_MODE (orogen_carry_v2 |
//! extras_v1), EX_OUT, optional EX_STEP_MY (dt sanity runs).

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
        "orogen_carry_v2" => PreviewOptions::orogen_carry_v2(),
        "extras_v1" => PreviewOptions::extras_v1(),
        other => panic!("unknown mode {other}"),
    }
}

#[test]
fn extras_profile_closes_replays_and_repeats() {
    let opts = PreviewOptions::extras_v1();
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
    for want in ["A-HOTSPOT-PLUMES", "A-CONTINENTAL-RIFT", "A-PLUME-HEAD-SPLIT", "A-COLLISION-SEDIMENT-ACCRETED"] {
        assert!(ids.contains(&want), "{want} missing");
    }
    assert!(!ids.contains(&"A-NO-HOTSPOTS"));
    let x = a.orogen.as_ref().unwrap().extras.as_ref().unwrap();
    assert!(x.volcanic_added_m.iter().any(|&h| h > 0.0));
    // Stepper and module replay equal the batch.
    let env = PreviewEnv::new(cfg.profile_with(&opts).unwrap(), cfg.grid, cfg.seed_path()).unwrap();
    let mut stepper = PreviewRun::init(env).unwrap();
    let mut plume = 0.0;
    while stepper.steps_taken() < stepper.steps_scheduled() {
        plume += stepper.step().unwrap().advance.hotspot_created_m3;
    }
    assert!(plume > 0.0);
    let stepped = stepper.finish().unwrap();
    assert_eq!(outputs_digest(&stepped).unwrap(), da);
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut run = ModuleRun::<M>::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
    run.run_to_end(None).unwrap();
    let out = M::finish(run.into_state()).unwrap();
    assert_eq!(outputs_digest(&out).unwrap(), da, "module != batch");
    eprintln!(
        "extras mini: plume crust {plume:.3e} m3, splits {}, plates {}",
        x.splits.len(),
        a.plates.len()
    );
}

#[test]
#[ignore]
fn seed_extras_trace() {
    let seed: u64 = std::env::var("EX_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("EX_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("EX_MODE").unwrap();
    let step_my: Option<f64> = std::env::var("EX_STEP_MY").ok().map(|s| s.parse().unwrap());
    let out_dir = PathBuf::from(std::env::var("EX_OUT").unwrap());
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
    let mut csv = std::io::BufWriter::new(std::fs::File::create(out_dir.join("steps.csv")).unwrap());
    writeln!(
        csv,
        "t_my,hotspot_m3,head_m3,heads,rift_births,breakup_births,convergent_skipped,rift_moved_m3,ridge_cells,adv_residual_m3,max_elev_m,land_frac,boundary_pairs"
    )
    .unwrap();
    while run.steps_taken() < run.steps_scheduled() {
        let r = run.step().unwrap();
        let a = &r.advance;
        let log = run.logs().last().unwrap();
        writeln!(
            csv,
            "{},{:.6e},{:.6e},{},{},{},{},{:.6e},{},{:.6e},{:.1},{:.5},{}",
            r.t_my,
            a.hotspot_created_m3,
            a.hotspot_head_m3,
            a.hotspot_heads,
            a.rift_births.len(),
            a.rift_breakup_births,
            a.rift_convergent_skipped,
            a.rift_moved_m3,
            a.ridge_cells,
            a.residual_m3(),
            log.max_elevation_m,
            log.land_fraction,
            r.boundary_upserts.len(),
        )
        .unwrap();
    }
    csv.flush().unwrap();
    let plumes: Vec<serde_json::Value> = run
        .env()
        .plumes_for_inspection()
        .iter()
        .map(|p| {
            let (lat, lon) = pf_geodynamics::rotation::unit_to_lat_lon(p.position).unwrap();
            serde_json::json!({"id": p.id, "lat": lat, "lon": lon, "onset_my": p.onset_my})
        })
        .collect();
    let out = run.finish().unwrap();
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    let kind: Vec<u8> = out.state.kind.iter().map(|k| *k as u8).collect();
    std::fs::write(out_dir.join("kind.u8"), kind).unwrap();
    let plate: Vec<f64> = out.state.plate.iter().map(|p| p.0 as f64).collect();
    write_f32(&out_dir.join("plate.f32"), &plate);
    let age: Vec<f64> = (0..out.state.cells()).map(|c| out.state.age_my(c)).collect();
    write_f32(&out_dir.join("age_my.f32"), &age);
    let mut splits = serde_json::Value::Null;
    if let Some(x) = out.orogen.as_ref().and_then(|o| o.extras.as_ref()) {
        write_f32(&out_dir.join("volcanic_added_m.f32"), &x.volcanic_added_m);
        write_f32(&out_dir.join("volcanic_t_my.f32"), &x.volcanic_t_my);
        write_f32(&out_dir.join("rift_t_my.f32"), &x.rift_t_my);
        write_f32(&out_dir.join("rift_thinned_m.f32"), &x.rift_thinned_m);
        splits = serde_json::to_value(&x.splits).unwrap();
    }
    let cv = &out.terrain.crust_volume;
    let meta = serde_json::json!({
        "seed": seed, "mode": mode, "nlat": grid.nlat(), "nlon": grid.nlon(),
        "t_my": out.state.t_my,
        "row_area_m2": (0..grid.nlat()).map(|i| areas[i * grid.nlon()]).collect::<Vec<_>>(),
        "plumes": plumes, "splits": splits,
        "plates": out.plates.iter().map(|p| serde_json::json!({"id": p.plate_id.0, "parent": p.parent_plate.map(|q| q.0), "birth_t_my": p.birth_t_my})).collect::<Vec<_>>(),
        "closure_plug_m3": cv.closure_plug_m3, "initial_m3": cv.initial_m3,
        "arc_created_m3": cv.arc_created_m3, "ridge_created_m3": cv.ridge_created_m3,
        "land_area_fraction": out.terrain.summary.land_area_fraction,
        "max_elevation_m": out.terrain.summary.max_elevation_m,
        "outputs_digest": outputs_digest(&out).unwrap(),
    });
    std::fs::write(out_dir.join("meta.json"), serde_json::to_string_pretty(&meta).unwrap()).unwrap();
}
