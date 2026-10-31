//! SEED-TRANSFORMS (lane TF, 2026-10-13): `PreviewOptions::world_v2` on the
//! mini grid (36 x 72, 2 My steps):
//! - opt-in: the law set carries `smoothed-normal-transform-v1` and
//!   `plate-reorg-scripted-v4`; no existing profile serializes either key;
//!   the boundary law is refused without the parcel carry;
//! - the run is deterministic and prefix-stable, the parcel ledger closes,
//!   the transform rules act (gap fills, non-convergent overlaps);
//! - finish, genesis-module replay and checkpoint / restore / continue equal
//!   the batch.

use pf_contracts::laws::{PlateBoundaryLawV1, PlateReorganizationLawV1};
use pf_contracts::params::{FidelityPreset, PlateBoundaryBlock, ScenarioPreset};
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
fn world_v2_is_opt_in_and_selects_the_new_versions() {
    let cfg = mini(42, 20.0);
    let u = cfg.profile_with(&PreviewOptions::world_v2()).unwrap().universe().unwrap();
    assert!(matches!(u.law_set.plate_reorganization, PlateReorganizationLawV1::PlateReorgScriptedV4 { .. }));
    assert!(matches!(u.law_set.plate_boundary, PlateBoundaryLawV1::SmoothedNormalTransformV1 { .. }));
    for o in [
        PreviewOptions::cycle_v1(),
        PreviewOptions::cycle_v2(),
        PreviewOptions::cycle_ranges_uplands_v1(),
        PreviewOptions::cycle_magma_v1(),
        PreviewOptions::world_v1(),
    ] {
        let p = cfg.profile_with(&o).unwrap();
        let s = serde_json::to_string(&p).unwrap();
        let us = serde_json::to_string(&p.universe().unwrap()).unwrap();
        for key in ["plate_boundary", "smoothed-normal-transform-v1", "plate-reorg-scripted-v4", "initiation_slab_age_min_my", "law_version\":4", "transform_obliquity_deg"] {
            assert!(!s.contains(key) && !us.contains(key), "{key} in an existing profile");
        }
    }
    // The boundary law needs the parcel carry (its advance runs there).
    let mut p = cfg.profile_with(&PreviewOptions::orogen_v1()).unwrap();
    assert!(p.universe().is_ok());
    let belt = p.orogeny.as_mut().and_then(|o| o.rate_limited_belt.as_mut()).expect("belt");
    assert!(belt.continental_parcel_carry.is_none());
    belt.plate_boundary = Some(PlateBoundaryBlock::v1());
    let e = format!("{:?}", p.universe().unwrap_err());
    assert!(e.contains("smoothed-normal-transform-v1"), "{e}");
}

fn run_to(cfg: &PreviewConfig, stop_steps: u32) -> (PreviewRun, u64, u64) {
    let params = cfg.profile_with(&PreviewOptions::world_v2()).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let (mut filled, mut overlaps) = (0, 0);
    while run.steps_taken() < stop_steps {
        let r = run.step().unwrap();
        let c = r.conserve.as_ref().expect("conserve ledger");
        assert!(c.volume_residual_m3.abs() <= 1e-9 * c.parcel_volume_m3.max(1.0), "{c:?}");
        let t = &r.advance.transform;
        assert!(t.boundary_cells > 0);
        filled += t.filled_transform_cells + t.filled_convergent_cells;
        overlaps += t.nonconvergent_overlaps;
    }
    (run, filled, overlaps)
}

#[test]
fn world_v2_is_deterministic_and_prefix_stable() {
    let (a, fa, oa) = run_to(&mini(42, 40.0), 20);
    let (b, fb, ob) = run_to(&mini(42, 80.0), 20);
    assert!(fa > 0 && oa > 0, "fills {fa} overlaps {oa}");
    assert_eq!((fa, oa), (fb, ob));
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(a.state().plate, b.state().plate);
    assert_eq!(a.state().boundary_class, b.state().boundary_class);
    assert_eq!(bits(&a.state().birth_t_my), bits(&b.state().birth_t_my));
    assert_eq!(bits(a.elevation_m()), bits(b.elevation_m()));
    // Transform cells exist at the step end.
    let tr = pf_contracts::terrain::BoundaryClassCodeV0::Transform.code();
    assert!(a.state().boundary_class.iter().any(|&k| k == tr));
}

#[test]
fn world_v2_finish_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::world_v2();
    let cfg = mini(240907, 60.0);
    let batch = simulate_with(&cfg, &opts).unwrap();
    batch.resource.validate().unwrap();
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
