//! WS-MODULE (lane SM, 2026-10-12): the preview chain as a genesis
//! SimModule (ADR 039 §4, acceptance 4, 5 and 10; design review C2-C4), on the
//! off-ladder mini grid the other preview tests use (36 x 72, 20 My in
//! 2 My steps). Full-size proofs (L0/L1, separate processes) are the
//! ignored drivers in `ws_module_capture.rs`.
//!
//! - batch equals module replay, bit for bit (legacy and both flag modes);
//! - checkpoint at step k, restore through bytes, continue = uninterrupted
//!   (outputs and the receipts after k);
//! - state bytes are canonical (decode then encode is the identity);
//! - one pinned recipe replays three times byte-identically, equals the
//!   pinned batch, and pinning changes no physics;
//! - every declared conserved quantity closes every step (C4);
//! - typed refusals: schedule, envelope, binding, law source, tampering.

use pf_contracts::params::{Compensation, FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{
    AuthorityMode, CheckpointV1, ModuleRun, RunBinding, SimModule, SimModuleError, StepDt, Tempo,
};
use pf_lab::preview::{
    publish, simulate_profile, simulate_with, PreviewConfig, PreviewOptions, PreviewOutputs,
};
use pf_lab::preview_module::{
    outputs_digest, preview_recipe, PreviewGenesisModule as M, PreviewObservation,
    PreviewObserveRequest, PreviewRecipe, CONSERVED,
};
use std::path::PathBuf;
use std::time::Instant;

type Run = ModuleRun<M>;

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

fn modes() -> Vec<(&'static str, PreviewOptions, u64)> {
    vec![
        ("legacy", PreviewOptions::legacy(), 240907),
        ("legacy", PreviewOptions::legacy(), 7),
        ("two_column", PreviewOptions::two_column_v0(), 7),
        ("two_column", PreviewOptions::two_column_v0(), 240907),
        ("bodies", PreviewOptions::crust_bodies_v0(), 3),
        ("bodies", PreviewOptions::crust_bodies_v0(), 240907),
    ]
}

fn scratch(tag: &str) -> PathBuf {
    let base = std::env::var_os("WS_MODULE_SCRATCH")
        .map(PathBuf::from)
        .unwrap_or_else(std::env::temp_dir);
    let dir = base.join(format!("ws-module-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    dir
}

fn replay(
    recipe: PreviewRecipe,
) -> (
    PreviewOutputs,
    Vec<pf_contracts::sim_module::StepReceiptV1>,
    String,
) {
    let mut run = Run::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
    run.run_to_end(None).unwrap();
    let final_cp = run.checkpoint().unwrap().digest().unwrap();
    let receipts = run.receipts().to_vec();
    (M::finish(run.into_state()).unwrap(), receipts, final_cp)
}

fn published(o: &PreviewOutputs, tag: &str) -> pf_lab::preview::PublishedDigests {
    let root = scratch(tag);
    let d = publish(o, &root, Instant::now()).unwrap();
    let _ = std::fs::remove_dir_all(&root);
    d
}

#[test]
fn ws_module_batch_equals_module_replay_bitwise() {
    for (mode, opts, seed) in modes() {
        let cfg = mini(seed);
        let batch = simulate_with(&cfg, &opts).unwrap();
        let (module, receipts, _) = replay(preview_recipe(&cfg, &opts, false).unwrap());
        assert_eq!(receipts.len(), 10, "{mode} {seed}");
        assert_eq!(module.elevation_m, batch.elevation_m, "{mode} {seed}");
        assert_eq!(module.state, batch.state, "{mode} {seed}");
        assert_eq!(module.terrain, batch.terrain, "{mode} {seed}");
        assert_eq!(module.hydrology, batch.hydrology, "{mode} {seed}");
        assert_eq!(
            module.exchange_epochs, batch.exchange_epochs,
            "{mode} {seed}"
        );
        assert_eq!(
            outputs_digest(&module).unwrap(),
            outputs_digest(&batch).unwrap(),
            "{mode} seed {seed}: module replay differs from batch"
        );
        assert_eq!(
            published(&module, "m"),
            published(&batch, "b"),
            "{mode} seed {seed}: published digests differ"
        );
    }
}

#[test]
fn ws_module_checkpoint_restore_continue_equals_uninterrupted() {
    for ((_, opts, seed), k) in modes().into_iter().zip([4u64, 1, 5, 9, 9, 3]) {
        let cfg = mini(seed);
        let recipe = preview_recipe(&cfg, &opts, false).unwrap();
        let (full, full_receipts, full_cp) = replay(recipe.clone());

        let mut a = Run::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
        a.run_until(k, None).unwrap();
        let bytes = a.checkpoint().unwrap().to_bytes().unwrap();
        drop(a);
        let cp = CheckpointV1::from_bytes(&bytes).unwrap();
        assert_eq!(cp.header.clock.step_index(), k);
        let mut b = Run::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
        b.run_to_end(None).unwrap();
        assert_eq!(
            b.checkpoint().unwrap().digest().unwrap(),
            full_cp,
            "seed {seed} k {k}"
        );
        assert_eq!(
            b.receipts(),
            &full_receipts[k as usize..],
            "seed {seed} k {k}: receipts after the restore differ"
        );
        let resumed = M::finish(b.into_state()).unwrap();
        assert_eq!(
            outputs_digest(&resumed).unwrap(),
            outputs_digest(&full).unwrap(),
            "seed {seed}: continuation differs from the uninterrupted run"
        );
    }
}

#[test]
fn ws_module_state_bytes_are_canonical() {
    for (_, opts, seed) in modes().into_iter().take(5).skip(2) {
        let recipe = preview_recipe(&mini(seed), &opts, false).unwrap();
        let mut run = Run::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
        for k in [0u64, 3, 10] {
            run.run_until(k, None).unwrap();
            let bytes = M::encode_state(run.state()).unwrap();
            let again = M::encode_state(&M::decode_state(&bytes).unwrap()).unwrap();
            assert!(
                bytes == again,
                "seed {seed} k {k}: decode then encode changed the bytes"
            );
            assert_eq!(
                run.state_digest(),
                pf_contracts::sim_module::digest_bytes(&bytes)
            );
        }
    }
}

#[test]
fn ws_module_pinned_recipe_replays_three_times_identically() {
    for (opts, seed) in [
        (PreviewOptions::legacy(), 240907u64),
        (PreviewOptions::crust_bodies_v0(), 7),
    ] {
        let cfg = mini(seed);
        let recipe = preview_recipe(&cfg, &opts, true).unwrap();
        assert!(recipe.inputs.profile.pinned_universe.is_some());
        let first = replay(recipe.clone());
        let d0 = outputs_digest(&first.0).unwrap();
        for _ in 0..2 {
            let again = replay(recipe.clone());
            assert_eq!(outputs_digest(&again.0).unwrap(), d0);
            assert_eq!(again.1, first.1, "receipts differ between replays");
            assert_eq!(again.2, first.2, "final checkpoint differs between replays");
        }
        // Pinned batch = pinned module replay.
        let p = &recipe.inputs.profile;
        let batch = simulate_profile(p.clone(), cfg.grid, cfg.seed_path()).unwrap();
        assert_eq!(outputs_digest(&batch).unwrap(), d0);
        // Pinning moves the recipe (params digest) but not the physics.
        let unpinned = simulate_with(&cfg, &opts).unwrap();
        assert_eq!(unpinned.elevation_m, first.0.elevation_m);
        assert_eq!(unpinned.state, first.0.state);
        assert_eq!(unpinned.exchange_epochs, first.0.exchange_epochs);
        assert_ne!(unpinned.params, first.0.params);
        let mut stripped = first.0.params.clone();
        stripped.pinned_universe = None;
        assert_eq!(stripped, unpinned.params);
    }
}

#[test]
fn ws_module_conservation_closes_every_step() {
    for (mode, opts, seed) in modes() {
        let recipe = preview_recipe(&mini(seed), &opts, false).unwrap();
        let (out, receipts, _) = replay(recipe);
        let two_column = mode != "legacy";
        for r in &receipts {
            for q in CONSERVED {
                let n = r.closures.iter().filter(|c| c.quantity == q.id).count();
                let expect = match q.id {
                    "bedrock-volume" | "sediment-volume" => 2,
                    "ocean-water-volume" => {
                        if two_column {
                            2
                        } else {
                            0
                        }
                    }
                    _ => 1,
                };
                assert_eq!(n, expect, "{mode} {seed} {} at {:?}", q.id, r.clock_after);
            }
            for c in &r.closures {
                assert!(c.closes(), "{mode} {seed}: {c:?}");
                assert!(c.allowed > 0.0 && c.gross >= 0.0);
            }
        }
        // Not vacuous: material moved, and the published crust ledger's
        // plug sits inside the declared cumulative tolerance.
        let moved: f64 = receipts
            .iter()
            .flat_map(|r| &r.closures)
            .filter(|c| c.quantity == "erosion-transport")
            .map(|c| c.gross)
            .sum();
        assert!(moved > 0.0, "{mode} {seed}: nothing eroded");
        let cv = &out.terrain.crust_volume;
        let last = receipts.last().unwrap();
        let cum = last
            .closures
            .iter()
            .find(|c| c.quantity == "crust-volume-cumulative")
            .unwrap();
        assert_eq!(
            cum.residual, cv.closure_plug_m3,
            "{mode} {seed}: plug differs from the artifact"
        );
        assert!(cv.closure_plug_m3.abs() <= cum.allowed);
    }
}

#[test]
fn ws_module_observe_is_read_only() {
    let recipe = preview_recipe(&mini(7), &PreviewOptions::two_column_v0(), false).unwrap();
    let mut run = Run::start(recipe, RunBinding::LAB_SIM_TIME).unwrap();
    run.run_until(3, None).unwrap();
    let before = run.state_digest().to_string();
    for req in [
        PreviewObserveRequest::Summary,
        PreviewObserveRequest::Ledgers,
        PreviewObserveRequest::Field(pf_lab::preview_module::PreviewField::ElevationM),
    ] {
        let _ = run.observe(&req).unwrap();
    }
    match run.observe(&PreviewObserveRequest::Summary).unwrap() {
        PreviewObservation::Summary {
            step_index, t_my, ..
        } => {
            assert_eq!(step_index, 3);
            assert_eq!(t_my, 6.0);
        }
        other => panic!("{other:?}"),
    }
    let after = pf_contracts::sim_module::digest_bytes(&M::encode_state(run.state()).unwrap());
    assert_eq!(before, after);
    // And the checkpoint still matches the runner's digest (no hidden mutation).
    run.checkpoint().unwrap();
}

#[test]
fn ws_module_typed_refusals() {
    let cfg = mini(7);
    let recipe = preview_recipe(&cfg, &PreviewOptions::legacy(), false).unwrap();

    // Another dt partition is guarantee 3, not replay: refused.
    let mut r2 = recipe.clone();
    r2.schedule.dt = StepDt::SimTime {
        amount: 1.0,
        unit: pf_contracts::sim_module::SimTimeUnit::MillionYears,
    };
    r2.schedule.steps = 20;
    assert!(matches!(
        Run::start(r2, RunBinding::GENESIS_SIM_TIME),
        Err(SimModuleError::Schedule(_))
    ));
    // Flexure is outside the preview envelope (never silently Airy).
    let mut r3 = recipe.clone();
    r3.inputs.profile.mantle.compensation = Compensation::Flexure;
    let e = Run::start(r3, RunBinding::GENESIS_SIM_TIME).err().unwrap();
    assert!(matches!(e, SimModuleError::Unsupported(_)), "{e}");
    // An out-of-envelope constant is the universe's typed refusal.
    let mut r4 = recipe.clone();
    r4.inputs.profile.planet.mantle_density_kg_m3 = 1.0e6;
    assert!(Run::start(r4, RunBinding::GENESIS_SIM_TIME).is_err());
    // Bindings: SimTime never in-world; this module has no WorldTime.
    let e = Run::start(
        recipe.clone(),
        RunBinding {
            tempo: Tempo::SimTime,
            mode: AuthorityMode::InWorld,
        },
    )
    .err()
    .unwrap();
    assert!(matches!(e, SimModuleError::Binding(_)), "{e}");
    let e = Run::start(
        recipe.clone(),
        RunBinding {
            tempo: Tempo::WorldTime,
            mode: AuthorityMode::Lab,
        },
    )
    .err()
    .unwrap();
    assert!(matches!(e, SimModuleError::Binding(_)), "{e}");

    // Tampered checkpoint, wrong recipe, wrong module, step past schedule.
    let mut run = Run::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    run.run_until(2, None).unwrap();
    let mut bytes = run.checkpoint().unwrap().to_bytes().unwrap();
    let last = bytes.len() - 1;
    bytes[last] ^= 1;
    assert!(matches!(
        CheckpointV1::from_bytes(&bytes),
        Err(SimModuleError::Checkpoint(_))
    ));
    let cp = run.checkpoint().unwrap();
    let other = preview_recipe(&mini(8), &PreviewOptions::legacy(), false).unwrap();
    assert!(matches!(
        Run::resume(other, RunBinding::GENESIS_SIM_TIME, &cp),
        Err(SimModuleError::Checkpoint(_))
    ));
    let mut wrong = cp.clone();
    wrong.header.module_id = "pf.genesis.something-else".into();
    assert!(matches!(
        M::restore(&wrong),
        Err(SimModuleError::ModuleMismatch { .. })
    ));
    run.run_to_end(None).unwrap();
    assert!(matches!(
        run.step_next(None),
        Err(SimModuleError::Schedule(_))
    ));
    // Finish before the schedule ends is refused.
    let mut early = Run::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
    early.run_until(3, None).unwrap();
    assert!(M::finish(early.into_state()).is_err());
}
