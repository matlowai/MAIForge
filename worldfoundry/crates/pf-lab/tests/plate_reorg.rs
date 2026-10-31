//! SEED-PLATE-REORG (lane PR, 2026-10-13): `PreviewOptions::cycle_v1`
//! (`oldstart_retain_v1` + `plate-reorg-scripted-v1`) on the mini grid
//! (36 x 72, 2 My steps):
//! - the run is independent of its scheduled length: the first 60 My of a
//!   120 My run equal a 60 My run bit for bit (re-solve schedule and plume
//!   renewal are prefix-stable);
//! - at least one scheduled re-solve happens by 120 My and the parcel
//!   ledger still closes;
//! - the genesis module replays the batch bit for bit and a checkpoint /
//!   restore / continue equals the uninterrupted run (the reorg state is
//!   checkpointed).

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

fn run_to(cfg: &PreviewConfig, stop_steps: u32) -> (PreviewRun, u32) {
    let params = cfg.profile_with(&PreviewOptions::cycle_v1()).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let mut resolves = 0;
    while run.steps_taken() < stop_steps {
        let r = run.step().unwrap();
        let c = r.conserve.expect("conserve ledger under cycle_v1");
        assert!(c.volume_residual_m3.abs() <= 1e-9 * c.parcel_volume_m3.max(1.0), "{c:?}");
        resolves += r.reorg.expect("reorg report under cycle_v1").resolves;
    }
    (run, resolves)
}

#[test]
fn cycle_v1_is_independent_of_run_length() {
    let (a, _) = run_to(&mini(42, 60.0), 30);
    let (b, resolves) = run_to(&mini(42, 120.0), 30);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(a.state().plate, b.state().plate);
    assert_eq!(bits(&a.state().thickness_m), bits(&b.state().thickness_m));
    assert_eq!(bits(a.elevation_m()), bits(b.elevation_m()));
    assert_eq!(a.plates(), b.plates());
    let _ = resolves;
    let (_, resolves) = run_to(&mini(42, 120.0), 60);
    assert!(resolves >= 1, "a re-solve by 120 My (gaps 50-100 My)");
}

#[test]
fn cycle_v1_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::cycle_v1();
    let cfg = mini(240907, 120.0);
    let batch = simulate_with(&cfg, &opts).unwrap();
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    full.run_to_end(None).unwrap();
    let full_cp = full.checkpoint().unwrap().digest().unwrap();
    let full_out = M::finish(full.into_state()).unwrap();
    assert_eq!(outputs_digest(&full_out).unwrap(), outputs_digest(&batch).unwrap());

    let mut a = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    a.run_until(37, None).unwrap();
    let bytes = a.checkpoint().unwrap().to_bytes().unwrap();
    drop(a);
    let cp = CheckpointV1::from_bytes(&bytes).unwrap();
    let mut b = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    b.run_to_end(None).unwrap();
    assert_eq!(b.checkpoint().unwrap().digest().unwrap(), full_cp);
    let resumed = M::finish(b.into_state()).unwrap();
    assert_eq!(outputs_digest(&resumed).unwrap(), outputs_digest(&full_out).unwrap());
}
