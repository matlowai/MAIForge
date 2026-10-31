//! SEED-CARRY-CONSERVE (lane CV, 2026-10-13): `PreviewOptions::conserve_v1`
//! (`continental-parcel-carry-v1` + `contact-length-pole-share-v1`) on the
//! mini grid (36 x 72, 20 My in 2 My steps):
//! - the per-step parcel ledger closes (area changes only by the declared
//!   rift source and collision sink; volume residual at rounding);
//! - the genesis module replays the batch bit for bit and a checkpoint /
//!   restore / continue equals the uninterrupted run (the parcels and the
//!   suture state are checkpointed).

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding};
use pf_lab::preview::{simulate_with, PreviewConfig, PreviewEnv, PreviewOptions, PreviewRun};
use pf_lab::preview_module::{outputs_digest, preview_recipe, PreviewGenesisModule as M};

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
fn conserve_v1_parcel_ledger_closes_every_step() {
    for seed in [42u64, 240907] {
        let cfg = mini(seed);
        let params = cfg.profile_with(&PreviewOptions::conserve_v1()).unwrap();
        let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
        let mut run = PreviewRun::init(env).unwrap();
        let mut area_prev: Option<f64> = None;
        while run.steps_taken() < run.steps_scheduled() {
            let r = run.step().unwrap();
            let c = r.conserve.expect("conserve ledger under conserve_v1");
            let scale = c.parcel_volume_m3.max(1.0);
            assert!(c.volume_residual_m3.abs() <= 1e-9 * scale, "seed {seed}: {c:?}");
            assert!(c.parcel_sediment_m3 >= 0.0 && c.lost_cells == 0, "seed {seed}: {c:?}");
            if let Some(a) = area_prev {
                let expect = a - c.consumed_area_m2 + c.created_area_m2;
                assert!((c.parcel_area_m2 - expect).abs() <= 1e-9 * a, "seed {seed}: area {c:?}");
            }
            // The grid projection tracks the parcels (coast rasterisation).
            assert!(
                (c.grid_cont_area_m2 - c.parcel_area_m2).abs() <= 0.05 * c.parcel_area_m2,
                "seed {seed}: grid {} vs parcels {}",
                c.grid_cont_area_m2,
                c.parcel_area_m2
            );
            assert!(r.continental.cont_area_m2 > 0.0);
            area_prev = Some(c.parcel_area_m2);
        }
    }
}

#[test]
fn conserve_v1_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::conserve_v1();
    let cfg = mini(42);
    let batch = simulate_with(&cfg, &opts).unwrap();
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    full.run_to_end(None).unwrap();
    let full_cp = full.checkpoint().unwrap().digest().unwrap();
    let full_out = M::finish(full.into_state()).unwrap();
    assert_eq!(outputs_digest(&full_out).unwrap(), outputs_digest(&batch).unwrap());

    let mut a = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    a.run_until(4, None).unwrap();
    let bytes = a.checkpoint().unwrap().to_bytes().unwrap();
    drop(a);
    let cp = CheckpointV1::from_bytes(&bytes).unwrap();
    let mut b = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    b.run_to_end(None).unwrap();
    assert_eq!(b.checkpoint().unwrap().digest().unwrap(), full_cp);
    let resumed = M::finish(b.into_state()).unwrap();
    assert_eq!(outputs_digest(&resumed).unwrap(), outputs_digest(&full_out).unwrap());
}
