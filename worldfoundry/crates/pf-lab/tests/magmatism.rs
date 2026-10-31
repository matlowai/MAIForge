//! SEED-MAGMATISM (lane MG, 2026-10-13): `PreviewOptions::cycle_magma_v1`
//! (`cycle_v1` + `magmatism = slab-arc-backarc-llsvp-v1`) on the mini grid
//! (36 x 72, 2 My steps):
//! - the law is opt-in: cycle_v1's universe carries no magmatism entry and
//!   its plumes are not moved;
//! - arcs are built (arc crust > 0, booked in the crust ledger), the parcel
//!   ledger closes, the run is deterministic and prefix-stable;
//! - finish appends the seven magmatic families and a valid volcano record;
//! - the genesis module replays the batch and checkpoint / restore / continue
//!   equals the uninterrupted run (the magma state is checkpointed).

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding};
use pf_contracts::terrain::ResourceFamilyKindV0;
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

fn run_to(cfg: &PreviewConfig, opts: &PreviewOptions, stop_steps: u32) -> (PreviewRun, f64) {
    let params = cfg.profile_with(opts).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let mut arc = 0.0;
    while run.steps_taken() < stop_steps {
        let r = run.step().unwrap();
        let c = r.conserve.expect("conserve ledger");
        assert!(c.volume_residual_m3.abs() <= 1e-9 * c.parcel_volume_m3.max(1.0), "{c:?}");
        arc += r.advance.arc_created_m3;
    }
    (run, arc)
}

#[test]
fn magmatism_is_opt_in() {
    let cfg = mini(42, 20.0);
    let old = cfg.profile_with(&PreviewOptions::cycle_v1()).unwrap();
    let u = old.universe().unwrap();
    assert!(u.law_set.magmatism.is_legacy());
    assert!(!serde_json::to_string(&u.law_set).unwrap().contains("magmatism"));
    assert!(!serde_json::to_string(&old).unwrap().contains("magmatism"));
    let new = cfg.profile_with(&PreviewOptions::cycle_magma_v1()).unwrap();
    assert!(!new.universe().unwrap().law_set.magmatism.is_legacy());
    // cycle_v1 plumes are where they always were; the law moves them.
    let e_old = PreviewEnv::new(old, cfg.grid, cfg.seed_path()).unwrap();
    let e_new = PreviewEnv::new(new, cfg.grid, cfg.seed_path()).unwrap();
    let (a, b) = (e_old.plumes_for_inspection(), e_new.plumes_for_inspection());
    assert_eq!(a.len(), b.len());
    assert!(a.iter().zip(b).all(|(x, y)| x.id == y.id && x.onset_my == y.onset_my));
    assert!(a.iter().zip(b).any(|(x, y)| x.position != y.position));
}

#[test]
fn magmatism_builds_arcs_deterministically_and_prefix_stable() {
    let opts = PreviewOptions::cycle_magma_v1();
    let (a, arc_a) = run_to(&mini(42, 40.0), &opts, 20);
    let (b, arc_b) = run_to(&mini(42, 80.0), &opts, 20);
    assert!(arc_a > 0.0, "arc crust is created");
    assert_eq!(arc_a.to_bits(), arc_b.to_bits());
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(a.state().plate, b.state().plate);
    assert_eq!(bits(&a.state().thickness_m), bits(&b.state().thickness_m));
    assert_eq!(bits(a.elevation_m()), bits(b.elevation_m()));
    assert_eq!(a.magma_state(), b.magma_state());
    let ms = a.magma_state().expect("magma state");
    assert_eq!(ms.logs.len(), 20);
    assert!(ms.logs.iter().any(|l| l.fronts > 0 && l.arc_m3 > 0.0));
    let v = a.volcano_features().unwrap().expect("volcano record");
    v.validate().unwrap();
    assert!(!v.features.is_empty());
}

#[test]
fn magmatism_finish_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::cycle_magma_v1();
    let cfg = mini(240907, 60.0);
    let batch = simulate_with(&cfg, &opts).unwrap();
    for k in ResourceFamilyKindV0::MAGMATIC_V1 {
        assert!(batch.resource.families.iter().any(|f| f.kind == k), "{k:?}");
        assert!(batch.potentials.iter().any(|(x, v)| *x == k && v.len() == 36 * 72));
    }
    batch.resource.validate().unwrap();
    batch.volcanoes.as_ref().expect("volcano record").validate().unwrap();
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    full.run_to_end(None).unwrap();
    let full_cp = full.checkpoint().unwrap().digest().unwrap();
    let full_out = M::finish(full.into_state()).unwrap();
    assert_eq!(outputs_digest(&full_out).unwrap(), outputs_digest(&batch).unwrap());
    let mut a = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    a.run_until(17, None).unwrap();
    let bytes = a.checkpoint().unwrap().to_bytes().unwrap();
    drop(a);
    let cp = CheckpointV1::from_bytes(&bytes).unwrap();
    let mut b = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    b.run_to_end(None).unwrap();
    assert_eq!(b.checkpoint().unwrap().digest().unwrap(), full_cp);
    let resumed = M::finish(b.into_state()).unwrap();
    assert_eq!(outputs_digest(&resumed).unwrap(), outputs_digest(&full_out).unwrap());
    assert_eq!(resumed.volcanoes, full_out.volcanoes);
}
