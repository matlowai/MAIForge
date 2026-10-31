//! SEED-STRATIGRAPHY stage A.1 (lane SA, 2026-10-13): `PreviewOptions::
//! world_strata_v2` (world_v1 + `stratigraphy = parcel-stacks-v2`) on the
//! mini grid (36 x 72, 2 My steps):
//! - opt-in: only the new profile carries parcel-stacks-v2; world_strata_v1
//!   keeps parcel-stacks-v1 and its serialization carries no v2 field;
//! - passive record: the geography is world_v1's bit for bit;
//! - t = 0 platform cover present and volume neutral; the ledger closes every
//!   step; one stack per parcel; prefix-stable;
//! - the published record carries the rock-uplift field and its gradient
//!   (v1 records do not); checkpoint / restore / continue and module replay
//!   give the same stratigraphy.

use pf_contracts::laws::StratigraphyLawV1;
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

fn start(cfg: &PreviewConfig, o: &PreviewOptions) -> PreviewRun {
    let params = cfg.profile_with(o).unwrap();
    PreviewRun::init(PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap()).unwrap()
}

#[test]
fn world_strata_v2_is_opt_in_and_v1_unchanged() {
    let cfg = mini(42, 20.0);
    let u = cfg.profile_with(&PreviewOptions::world_strata_v2()).unwrap().universe().unwrap();
    assert!(matches!(u.law_set.stratigraphy, StratigraphyLawV1::ParcelStacksV2 { .. }));
    u.validate().unwrap();
    let p1 = cfg.profile_with(&PreviewOptions::world_strata_v1()).unwrap();
    let u1 = p1.universe().unwrap();
    assert!(matches!(u1.law_set.stratigraphy, StratigraphyLawV1::ParcelStacksV1 { .. }));
    let s1 = serde_json::to_string(&p1).unwrap() + &serde_json::to_string(&u1.law_set).unwrap();
    assert!(!s1.contains("strata_shelf_width_m") && !s1.contains("parcel-stacks-v2") && !s1.contains("\"v2\""), "v2 in v1");
    // The v1 record carries no uplift field.
    let mut r1 = start(&cfg, &PreviewOptions::world_strata_v1());
    r1.step().unwrap();
    let d1 = r1.stratigraphy(None).unwrap().unwrap();
    assert!(d1.rock_uplift_m.is_none());
    let j1 = serde_json::to_string(&d1).unwrap();
    assert!(!j1.contains("rock_uplift_m") && !j1.contains("uplift_grad"));
}

#[test]
fn world_strata_v2_is_passive_closes_and_publishes_uplift() {
    let cfg = mini(42, 60.0);
    let mut a = start(&cfg, &PreviewOptions::world_strata_v2());
    let mut w = start(&cfg, &PreviewOptions::world_v1());
    // t = 0: platform cover on the seeded parcels, volume neutral.
    let st0 = a.strata_stacks().unwrap();
    let covered = st0.iter().filter(|s| s.cover_m() > 0.0).count();
    assert!(covered * 4 > st0.len(), "platform cover on {covered} of {} parcels", st0.len());
    assert!(st0.iter().all(|s| s.cover_m() <= 2000.0 + 1e-6));
    let d0 = a.stratigraphy(None).unwrap().unwrap();
    d0.validate().unwrap();
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    let mut uplift_seen = 0.0f64;
    while a.steps_taken() < 25 {
        let r = a.step().unwrap();
        w.step().unwrap();
        let g = r.stratigraphy.as_ref().expect("strata ledger");
        assert!(g.closure_rel < 1e-9, "{g:?}");
        assert!(g.column_mismatch_rel < 1e-9, "{g:?}");
        assert!(g.entry_mismatch_max_rel < 1e-9, "{g:?}");
        assert_eq!(g.parcels, g.stacks);
        let v2 = g.v2.as_ref().expect("v2 ledger");
        uplift_seen = uplift_seen.max(v2.mean_abs_rock_uplift_m);
        assert_eq!(a.state().plate, w.state().plate);
        assert_eq!(a.state().kind, w.state().kind);
        assert_eq!(bits(&a.state().thickness_m), bits(&w.state().thickness_m));
        assert_eq!(bits(&a.state().sediment_m), bits(&w.state().sediment_m));
        assert_eq!(bits(a.elevation_m()), bits(w.elevation_m()));
    }
    assert!(uplift_seen > 0.0, "no rock uplift recorded");
    let d = a.stratigraphy(None).unwrap().expect("record");
    d.validate().unwrap();
    let n = (d.nlat * d.nlon) as usize;
    assert_eq!(d.rock_uplift_m.as_ref().map(|v| v.len()), Some(n));
    assert!(d.uplift_grad_east.as_ref().unwrap().iter().any(|x| *x != 0.0));
    assert!(d.uplift_grad_north.as_ref().unwrap().iter().any(|x| *x != 0.0));
    assert!(a.strata_stacks().unwrap().iter().all(|s| s.cover.len() + 1 <= 16));
    // Prefix-stable and deterministic.
    let mut b = start(&mini(42, 120.0), &PreviewOptions::world_strata_v2());
    while b.steps_taken() < 25 {
        b.step().unwrap();
    }
    assert_eq!(a.strata_stacks(), b.strata_stacks());
    assert_eq!(a.stratigraphy(None).unwrap(), b.stratigraphy(None).unwrap());
}

#[test]
fn world_strata_v2_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::world_strata_v2();
    let cfg = mini(240907, 40.0);
    let batch = simulate_with(&cfg, &opts).unwrap();
    let sd = batch.stratigraphy.as_ref().expect("stratigraphy record");
    sd.validate().unwrap();
    assert!(sd.rock_uplift_m.is_some());
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    full.run_to_end(None).unwrap();
    let full_cp = full.checkpoint().unwrap().digest().unwrap();
    let full_out = M::finish(full.into_state()).unwrap();
    assert_eq!(outputs_digest(&full_out).unwrap(), outputs_digest(&batch).unwrap());
    let mut a = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    a.run_until(11, None).unwrap();
    let bytes = a.checkpoint().unwrap().to_bytes().unwrap();
    drop(a);
    let cp = CheckpointV1::from_bytes(&bytes).unwrap();
    let mut b = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    b.run_to_end(None).unwrap();
    assert_eq!(b.checkpoint().unwrap().digest().unwrap(), full_cp);
    let resumed = M::finish(b.into_state()).unwrap();
    assert_eq!(resumed.stratigraphy, batch.stratigraphy);
    assert_eq!(outputs_digest(&resumed).unwrap(), outputs_digest(&full_out).unwrap());
}
