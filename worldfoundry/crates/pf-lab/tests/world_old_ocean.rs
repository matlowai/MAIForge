//! SEED-OLD-OCEAN (lane OR, 2026-10-13): `PreviewOptions::world_v3` on the
//! mini grid (36 x 72, 2 My steps):
//! - opt-in: the law set carries `plate-reorg-scripted-v5`; no existing
//!   profile (world_v2 included) serializes its keys;
//! - the run is deterministic and prefix-stable, passive-margin initiations
//!   happen, the parcel ledger closes;
//! - finish, genesis-module replay and checkpoint / restore / continue equal
//!   the batch across v5 events.

use pf_contracts::laws::PlateReorganizationLawV1;
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding};
use pf_lab::preview::{simulate_with, PreviewConfig, PreviewEnv, PreviewOptions, PreviewRun};
use pf_lab::preview_module::{outputs_digest, preview_recipe, PreviewGenesisModule as M};

fn mini(seed: u64, duration: f64) -> PreviewConfig {
    PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed,
        duration_my: Some(duration),
        step_my: Some(2.0),
    }
}

#[test]
fn world_v3_is_opt_in_and_selects_v5() {
    let cfg = mini(42, 20.0);
    let u = cfg.profile_with(&PreviewOptions::world_v3()).unwrap().universe().unwrap();
    assert!(matches!(u.law_set.plate_reorganization, PlateReorganizationLawV1::PlateReorgScriptedV5 { .. }));
    let s3 = serde_json::to_string(&u).unwrap();
    assert!(s3.contains("plate-reorg-scripted-v5") && s3.contains("margin_initiation_age_my"));
    for o in [
        PreviewOptions::cycle_v2(),
        PreviewOptions::cycle_magma_v1(),
        PreviewOptions::world_v1(),
        PreviewOptions::world_v2(),
    ] {
        let p = cfg.profile_with(&o).unwrap();
        let s = serde_json::to_string(&p).unwrap();
        let us = serde_json::to_string(&p.universe().unwrap()).unwrap();
        for key in ["plate-reorg-scripted-v5", "margin_initiation_age_my", "margin_segment_length_m", "margin_trench_spacing_m", "law_version\":5"] {
            assert!(!s.contains(key) && !us.contains(key), "{key} in an existing profile");
        }
    }
}

fn run_to(cfg: &PreviewConfig, stop_steps: u32) -> (PreviewRun, usize) {
    let params = cfg.profile_with(&PreviewOptions::world_v3()).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let mut margins = 0;
    while run.steps_taken() < stop_steps {
        let r = run.step().unwrap();
        let c = r.conserve.as_ref().expect("conserve ledger");
        assert!(c.volume_residual_m3.abs() <= 1e-9 * c.parcel_volume_m3.max(1.0), "{c:?}");
        if let Some(x) = &r.reorg {
            margins += x.margin_seeds.len();
            for (age, _, area, len) in &x.margin_seeds {
                assert!(*age >= 160.0 && *area > 0.0 && *len > 0.0, "{age} {area} {len}");
            }
        }
    }
    (run, margins)
}

#[test]
fn world_v3_is_deterministic_prefix_stable_and_initiates_at_margins() {
    let (a, ma) = run_to(&mini(42, 300.0), 120);
    let (b, mb) = run_to(&mini(42, 400.0), 120);
    assert!(ma > 0, "no passive-margin initiation in 240 My");
    assert_eq!(ma, mb);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(a.state().plate, b.state().plate);
    assert_eq!(a.state().boundary_class, b.state().boundary_class);
    assert_eq!(bits(&a.state().birth_t_my), bits(&b.state().birth_t_my));
    assert_eq!(bits(a.elevation_m()), bits(b.elevation_m()));
}

#[test]
fn world_v3_finish_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::world_v3();
    let cfg = mini(42, 240.0);
    let batch = simulate_with(&cfg, &opts).unwrap();
    batch.resource.validate().unwrap();
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    full.run_to_end(None).unwrap();
    let full_cp = full.checkpoint().unwrap().digest().unwrap();
    let full_out = M::finish(full.into_state()).unwrap();
    assert_eq!(outputs_digest(&full_out).unwrap(), outputs_digest(&batch).unwrap());
    let mut a = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    a.run_until(61, None).unwrap();
    let bytes = a.checkpoint().unwrap().to_bytes().unwrap();
    drop(a);
    let cp = CheckpointV1::from_bytes(&bytes).unwrap();
    let mut b = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    b.run_to_end(None).unwrap();
    assert_eq!(b.checkpoint().unwrap().digest().unwrap(), full_cp);
    let resumed = M::finish(b.into_state()).unwrap();
    assert_eq!(outputs_digest(&resumed).unwrap(), outputs_digest(&full_out).unwrap());
}
