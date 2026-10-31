//! SEED-WORLD-COMPOSE (lane CP, 2026-10-13): `PreviewOptions::world_v1` on the
//! mini grid (36 x 72, 2 My steps):
//! - opt-in: the composed law set carries `plate-reorg-scripted-v3`,
//!   `plateau-cap-spill-v2` and `slab-arc-backarc-llsvp-v2`, and none of the
//!   existing profiles does (their law sets serialize without the new keys);
//! - the run is deterministic and prefix-stable, the parcel ledger closes,
//!   the in-belt plateau pass runs and nothing exceeds the cap after it;
//! - finish, genesis-module replay and checkpoint / restore / continue equal
//!   the batch (the v2 back-arc episode markers are checkpointed).

use pf_contracts::laws::{MagmatismLawV1, OrogenPlateauLawV1, PlateReorganizationLawV1};
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
fn world_v1_is_opt_in_and_selects_the_new_versions() {
    let cfg = mini(42, 20.0);
    let u = cfg.profile_with(&PreviewOptions::world_v1()).unwrap().universe().unwrap();
    assert!(matches!(u.law_set.plate_reorganization, PlateReorganizationLawV1::PlateReorgScriptedV3 { .. }));
    assert!(matches!(u.law_set.orogen_plateau, OrogenPlateauLawV1::PlateauCapSpillV2 { .. }));
    assert!(matches!(u.law_set.magmatism, MagmatismLawV1::SlabArcBackarcLlsvpV2 { .. }));
    assert!(!u.law_set.flood_basalt.is_legacy());
    for o in [
        PreviewOptions::cycle_v1(),
        PreviewOptions::cycle_v2(),
        PreviewOptions::cycle_ranges_uplands_v1(),
        PreviewOptions::cycle_magma_v1(),
        PreviewOptions::extras_v3_flood(),
    ] {
        let p = cfg.profile_with(&o).unwrap();
        let s = serde_json::to_string(&p).unwrap();
        for key in ["law_version\":3", "plateau_foothill_gradient", "backarc_episode_my", "kimberlite_plume_radius_m", "subduction_initiation_age_my", "plateau-cap-spill-v2", "slab-arc-backarc-llsvp-v2", "plate-reorg-scripted-v3"] {
            assert!(!s.contains(key), "{key} in an existing profile");
        }
    }
}

fn run_to(cfg: &PreviewConfig, stop_steps: u32) -> (PreviewRun, u32) {
    let params = cfg.profile_with(&PreviewOptions::world_v1()).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let mut plateau_steps = 0;
    while run.steps_taken() < stop_steps {
        let r = run.step().unwrap();
        let c = r.conserve.as_ref().expect("conserve ledger");
        assert!(c.volume_residual_m3.abs() <= 1e-9 * c.parcel_volume_m3.max(1.0), "{c:?}");
        if let Some(p) = r.advance.plateau_belt.as_ref() {
            plateau_steps += 1;
            assert!(p.taper_residual_m < 1.0e3, "{p:?}");
        }
        assert!(r.shortening_budget.is_some());
    }
    (run, plateau_steps)
}

#[test]
fn world_v1_is_deterministic_and_prefix_stable() {
    let (a, pa) = run_to(&mini(42, 40.0), 20);
    let (b, pb) = run_to(&mini(42, 80.0), 20);
    assert_eq!(pa, 20);
    assert_eq!(pa, pb);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(a.state().plate, b.state().plate);
    assert_eq!(bits(&a.state().thickness_m), bits(&b.state().thickness_m));
    assert_eq!(bits(a.elevation_m()), bits(b.elevation_m()));
    assert_eq!(a.magma_state(), b.magma_state());
    let ms = a.magma_state().expect("magma state");
    assert_eq!(ms.backarc_first_t_my.len(), 36 * 72);
}

#[test]
fn world_v1_finish_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::world_v1();
    let cfg = mini(240907, 60.0);
    let batch = simulate_with(&cfg, &opts).unwrap();
    batch.resource.validate().unwrap();
    assert!(batch.resource.assumptions.iter().any(|a| a.id == "A-RES-KIMBERLITE-PLUME-V2"));
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
}
