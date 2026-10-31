//! SEED-FLOODBASALT (lane FB, 2026-10-13): the `extras_v3_flood` profile
//! (`extras_v3` + `flood_basalt = plume-head-lip-swell-v1`).
//!
//! The normal test runs the profile on lane OG's mini grid (36 x 72) through
//! the batch path, the stepper and the genesis module: the three agree, the
//! crust ledger closes, the LIP volume logged is the declared volume of the
//! heads that erupted, the LIP marker holds it, the swell is logged, the
//! assumption is recorded, and extras_v3 logs none of it.
//!
//! The ignored driver `flood_basalt_trace` runs one (seed, fidelity, mode)
//! and writes the final fields, the LIP marker, the swell term and the
//! plumes for maps and cross-sections. Env: FB_SEED, FB_FIDELITY (preview |
//! standard), FB_MODE (extras_v3 | extras_v3_flood), FB_OUT, optional
//! FB_STEP_MY, FB_DURATION_MY, FB_DUMP_AT.

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{ModuleRun, RunBinding};
use pf_contracts::terrain::CrustKindCodeV0;
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
        "extras_v3" => PreviewOptions::extras_v3(),
        "extras_v3_flood" => PreviewOptions::extras_v3_flood(),
        other => panic!("unknown mode {other}"),
    }
}

#[test]
fn flood_basalt_profile_closes_replays_and_books_the_lip() {
    let opts = PreviewOptions::extras_v3_flood();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed: 240907,
        duration_my: Some(40.0),
        step_my: Some(2.0),
    };
    let a = simulate_with(&cfg, &opts).unwrap();
    let da = outputs_digest(&a).unwrap();
    let cv = &a.terrain.crust_volume;
    assert!(
        cv.closure_plug_m3.abs() <= 1e-9 * cv.initial_m3,
        "closure plug {} of {}",
        cv.closure_plug_m3,
        cv.initial_m3
    );
    let ids: Vec<&str> = a.terrain.assumptions.iter().map(|x| x.id.as_str()).collect();
    for want in ["A-FLOOD-BASALT", "A-HOTSPOT-PLUMES", "A-TRENCH"] {
        assert!(ids.contains(&want), "{want} missing");
    }
    let og = a.orogen.as_ref().unwrap();
    let logs: Vec<_> = og
        .steps
        .iter()
        .map(|s| s.extras.as_ref().unwrap().flood_basalt.clone().unwrap())
        .collect();
    assert_eq!(logs.len(), og.steps.len());
    // The LIP volume logged is whole provinces of the declared volume
    // (eruptions straddling the run end aside), all inside the plume crust.
    let env = PreviewEnv::new(cfg.profile_with(&opts).unwrap(), cfg.grid, cfg.seed_path()).unwrap();
    let fb = *env.flood_basalt_params().unwrap();
    let v = fb.lip.total_volume_m3();
    let lip: f64 = logs.iter().map(|l| l.lip_created_m3).sum();
    let in_run = env
        .plumes_for_inspection()
        .iter()
        .filter(|p| p.onset_my >= 0.0 && p.onset_my + fb.lip.eruption_duration_my <= 40.0)
        .count();
    assert!(in_run > 0, "no head erupted inside the mini run");
    assert!(lip >= in_run as f64 * v * (1.0 - 1e-9), "{lip} vs {in_run} x {v}");
    for s in &og.steps {
        let x = s.extras.as_ref().unwrap();
        let f = x.flood_basalt.as_ref().unwrap();
        assert!(f.lip_created_m3 <= x.hotspot_created_m3 * (1.0 + 1e-12));
        assert!(f.lip_continental_m3 <= f.lip_created_m3 * (1.0 + 1e-12));
    }
    assert!(logs.iter().all(|l| l.max_swell_wl_m > 0.0), "a plume is active from t = 0 in this seed");
    // The carried LIP marker holds crust. Not bounded by the volume erupted:
    // the semi-Lagrangian carry duplicates or drops source columns (lane CV's
    // SEED-CARRY-CONSERVE), so only presence and order of magnitude are checked.
    let marker: f64 = og
        .extras
        .as_ref()
        .unwrap()
        .lip_added_m
        .as_ref()
        .unwrap()
        .iter()
        .zip(env.areas())
        .map(|(h, a)| h * a)
        .sum();
    assert!(marker > 0.2 * lip && marker < 2.0 * lip, "{marker} {lip}");
    // extras_v3 carries none of it.
    let v3 = simulate_with(&cfg, &PreviewOptions::extras_v3()).unwrap();
    let og3 = v3.orogen.as_ref().unwrap();
    assert!(og3.steps.iter().all(|s| s.extras.as_ref().unwrap().flood_basalt.is_none()));
    assert!(og3.extras.as_ref().unwrap().lip_added_m.is_none());
    assert!(!v3.terrain.assumptions.iter().any(|x| x.id == "A-FLOOD-BASALT"));
    // Stepper and module replay equal the batch.
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
    eprintln!(
        "flood mini: LIP {:.3e} m3 ({} in-run heads x {:.3e}), marker {:.3e} m3, max swell {:.0} m",
        lip,
        in_run,
        v,
        marker,
        logs.iter().map(|l| l.max_swell_wl_m).fold(0.0, f64::max)
    );
}

#[test]
#[ignore]
fn flood_basalt_trace() {
    let seed: u64 = std::env::var("FB_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("FB_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("FB_MODE").unwrap();
    let step_my: Option<f64> = std::env::var("FB_STEP_MY").ok().map(|s| s.parse().unwrap());
    let out_dir = PathBuf::from(std::env::var("FB_OUT").unwrap());
    std::fs::create_dir_all(&out_dir).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        // Optional shorter run (snapshots of the carried markers): FB_DURATION_MY.
        duration_my: std::env::var("FB_DURATION_MY").ok().filter(|s| !s.is_empty()).map(|s| s.parse().unwrap()),
        step_my,
    };
    let env = PreviewEnv::new(cfg.profile_with(&options(&mode)).unwrap(), cfg.grid, cfg.seed_path()).unwrap();
    let plumes: Vec<serde_json::Value> = env
        .plumes_for_inspection()
        .iter()
        .map(|p| serde_json::json!({"id": p.id, "position": p.position, "onset_my": p.onset_my}))
        .collect();
    let mut run = PreviewRun::init(env).unwrap();
    let grid = run.env().grid().clone();
    let areas = run.env().areas().to_vec();
    let total: f64 = areas.iter().sum();
    let mut csv = std::io::BufWriter::new(std::fs::File::create(out_dir.join("steps.csv")).unwrap());
    writeln!(
        csv,
        "t_my,land_frac,land_above_1km_frac_of_land,land_above_2km_frac_of_land,cont_area_frac,cont_crust_m3,max_elev_m"
    )
    .unwrap();
    // Optional snapshots: FB_DUMP_AT = comma-separated times (My).
    let dump_at: Vec<f64> = std::env::var("FB_DUMP_AT")
        .map(|v| v.split(',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect())
        .unwrap_or_default();
    while run.steps_taken() < run.steps_scheduled() {
        let r = run.step().unwrap();
        if dump_at.iter().any(|&t| (t - r.t_my).abs() < 1e-9) {
            let d = out_dir.join(format!("snap_{:.0}", r.t_my));
            std::fs::create_dir_all(&d).unwrap();
            write_f32(&d.join("elevation_m.f32"), run.elevation_m());
            write_f32(&d.join("thickness_m.f32"), &run.state().thickness_m);
            let kind: Vec<u8> = run.state().kind.iter().map(|k| *k as u8).collect();
            std::fs::write(d.join("kind.u8"), kind).unwrap();
        }
        let e = run.elevation_m();
        let st = run.state();
        let (mut land, mut a1, mut a2, mut ca, mut cv) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for c in 0..e.len() {
            if e[c] > 0.0 {
                land += areas[c];
                if e[c] > 1000.0 {
                    a1 += areas[c];
                }
                if e[c] > 2000.0 {
                    a2 += areas[c];
                }
            }
            if st.kind[c] == CrustKindCodeV0::Continental {
                ca += areas[c];
                cv += st.thickness_m[c] * areas[c];
            }
        }
        writeln!(
            csv,
            "{},{:.5},{:.5},{:.5},{:.5},{:.6e},{:.1}",
            r.t_my,
            land / total,
            a1 / land.max(1.0),
            a2 / land.max(1.0),
            ca / total,
            cv,
            e.iter().cloned().fold(f64::MIN, f64::max)
        )
        .unwrap();
    }
    let t_end = run.state().t_my;
    let swell = run.env().swell_at(t_end).unwrap();
    let out = run.finish().unwrap();
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    let kind: Vec<u8> = out.state.kind.iter().map(|k| *k as u8).collect();
    std::fs::write(out_dir.join("kind.u8"), kind).unwrap();
    std::fs::write(out_dir.join("boundary_class.u8"), &out.state.boundary_class).unwrap();
    let plate: Vec<f64> = out.state.plate.iter().map(|p| p.0 as f64).collect();
    write_f32(&out_dir.join("plate.f32"), &plate);
    let og = out.orogen.as_ref().unwrap();
    let x = og.extras.as_ref().unwrap();
    write_f32(&out_dir.join("volcanic_added_m.f32"), &x.volcanic_added_m);
    write_f32(&out_dir.join("volcanic_t_my.f32"), &x.volcanic_t_my);
    let n = out.state.cells();
    write_f32(&out_dir.join("lip_added_m.f32"), x.lip_added_m.as_deref().unwrap_or(&vec![0.0; n]));
    write_f32(&out_dir.join("swell_wl_m.f32"), swell.as_deref().unwrap_or(&vec![0.0; n]));
    let steps: Vec<serde_json::Value> = og
        .steps
        .iter()
        .map(|s| {
            let e = s.extras.as_ref().unwrap();
            serde_json::json!({
                "hotspot_created_m3": e.hotspot_created_m3,
                "hotspot_head_m3": e.hotspot_head_m3,
                "plume_heads": e.plume_heads,
                "plate_splits": e.plate_splits,
                "flood_basalt": e.flood_basalt,
            })
        })
        .collect();
    let cv = &out.terrain.crust_volume;
    let meta = serde_json::json!({
        "seed": seed, "mode": mode, "nlat": grid.nlat(), "nlon": grid.nlon(),
        "t_my": out.state.t_my, "step_my": step_my,
        "closure_plug_m3": cv.closure_plug_m3, "initial_m3": cv.initial_m3,
        "arc_created_m3": cv.arc_created_m3, "final_m3": cv.final_m3,
        "land_area_fraction": out.terrain.summary.land_area_fraction,
        "max_elevation_m": out.terrain.summary.max_elevation_m,
        "plumes": plumes,
        "steps": steps,
        "outputs_digest": outputs_digest(&out).unwrap(),
    });
    std::fs::write(out_dir.join("meta.json"), serde_json::to_string_pretty(&meta).unwrap()).unwrap();
}
