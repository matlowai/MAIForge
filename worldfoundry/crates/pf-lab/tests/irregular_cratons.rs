//! SEED-IRREGULAR-SEEDS (lane IS, 2026-10-13): `PreviewOptions::oldstart_v1`
//! (`conserve_v1` + `craton_seeding = irregular-cratons-v1`) on the mini
//! grid (36 x 72, 20 My in 2 My steps):
//! - the law set carries the craton law only when selected (conserve_v1's
//!   resolved profile does not mention it);
//! - init is deterministic and starts continental interiors low (median
//!   continental elevation 0.2-0.6 km, against conserve_v1's ~1.2 km);
//! - the genesis module replays the batch bit for bit and a checkpoint /
//!   restore / continue equals the uninterrupted run.

use pf_contracts::laws::CratonSeedingLawV1;
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding};
use pf_contracts::terrain::CrustKindCodeV0;
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

fn median_cont_elev(opts: &PreviewOptions, seed: u64) -> f64 {
    let cfg = mini(seed);
    let params = cfg.profile_with(opts).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let run = PreviewRun::init(env).unwrap();
    let st = run.state();
    let mut e: Vec<f64> = (0..st.kind.len())
        .filter(|&c| st.kind[c] == CrustKindCodeV0::Continental)
        .map(|c| run.elevation_m()[c])
        .collect();
    e.sort_by(|a, b| a.total_cmp(b));
    e[e.len() / 2]
}

#[test]
fn oldstart_v1_law_set_and_low_start() {
    let cfg = mini(42);
    let p = cfg.profile_with(&PreviewOptions::oldstart_v1()).unwrap();
    let u = p.universe().unwrap();
    assert!(matches!(u.law_set.craton_seeding, CratonSeedingLawV1::IrregularCratonsV1 { .. }));
    let q = cfg.profile_with(&PreviewOptions::conserve_v1()).unwrap();
    assert!(q.universe().unwrap().law_set.craton_seeding.is_legacy());
    let qj = serde_json::to_string(&q).unwrap();
    assert!(!qj.contains("irregular_cratons") && !qj.contains("craton_seeding"));
    for seed in [42u64, 240907] {
        let low = median_cont_elev(&PreviewOptions::oldstart_v1(), seed);
        let high = median_cont_elev(&PreviewOptions::conserve_v1(), seed);
        assert!(low > 200.0 && low < 600.0, "seed {seed}: oldstart median {low}");
        assert!(high > 1000.0, "seed {seed}: conserve median {high}");
    }
}

#[test]
fn oldstart_v1_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::oldstart_v1();
    let cfg = mini(42);
    let batch = simulate_with(&cfg, &opts).unwrap();
    let again = simulate_with(&cfg, &opts).unwrap();
    assert_eq!(outputs_digest(&batch).unwrap(), outputs_digest(&again).unwrap());
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
