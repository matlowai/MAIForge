//! SEED-STRATIGRAPHY stage A (lane ST, 2026-10-13): `PreviewOptions::
//! world_strata_v1` (world_v1 + `stratigraphy = parcel-stacks-v1`) on the
//! mini grid (36 x 72, 2 My steps):
//! - opt-in: only the new profile's law set carries the stratigraphy law;
//!   world_v1 and the earlier profiles serialize without it;
//! - passive record: the geography (plates, kinds, thickness, sediment,
//!   elevation) is world_v1's bit for bit;
//! - every step the stack ledger closes, there is one stack per parcel and
//!   the stacks equal the parcel columns; deterministic and prefix-stable;
//! - the published record validates; finish, module replay and
//!   checkpoint / restore / continue give the same stratigraphy.

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
fn world_strata_v1_is_opt_in() {
    let cfg = mini(42, 20.0);
    let u = cfg.profile_with(&PreviewOptions::world_strata_v1()).unwrap().universe().unwrap();
    assert!(matches!(u.law_set.stratigraphy, StratigraphyLawV1::ParcelStacksV1 { .. }));
    u.validate().unwrap();
    for o in [PreviewOptions::world_v1(), PreviewOptions::cycle_magma_v1(), PreviewOptions::conserve_v1(), PreviewOptions::legacy()] {
        let p = cfg.profile_with(&o).unwrap();
        let s = serde_json::to_string(&p).unwrap();
        assert!(!s.contains("strata_") && !s.contains("stratigraphy"), "stratigraphy in an existing profile");
        let ls = serde_json::to_string(&p.universe().unwrap().law_set).unwrap();
        assert!(!ls.contains("stratigraphy") && !ls.contains("parcel-stacks"), "{ls}");
    }
}

#[test]
fn world_strata_v1_is_passive_closes_and_is_deterministic() {
    let cfg = mini(42, 60.0);
    let mut a = start(&cfg, &PreviewOptions::world_strata_v1());
    let mut w = start(&cfg, &PreviewOptions::world_v1());
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    let mut layers_seen = 0.0f64;
    while a.steps_taken() < 25 {
        let r = a.step().unwrap();
        w.step().unwrap();
        let g = r.stratigraphy.as_ref().expect("strata ledger");
        assert!(g.closure_rel < 1e-9, "{g:?}");
        assert!(g.column_mismatch_rel < 1e-9, "{g:?}");
        assert!(g.entry_mismatch_max_rel < 1e-9, "{g:?}");
        assert_eq!(g.parcels, g.stacks);
        layers_seen = layers_seen.max(g.mean_cover_layers);
        // Passive: world_v1's geography.
        assert_eq!(a.state().plate, w.state().plate);
        assert_eq!(a.state().kind, w.state().kind);
        assert_eq!(bits(&a.state().thickness_m), bits(&w.state().thickness_m));
        assert_eq!(bits(&a.state().sediment_m), bits(&w.state().sediment_m));
        assert_eq!(bits(a.elevation_m()), bits(w.elevation_m()));
    }
    assert!(layers_seen > 0.1, "no cover written: {layers_seen}");
    let d = a.stratigraphy(None).unwrap().expect("record");
    d.validate().unwrap();
    assert!(d.thickness_m.iter().all(|t| *t <= 3000.0 + 1e-3));
    let full = a.stratigraphy(Some(f64::INFINITY)).unwrap().unwrap();
    full.validate().unwrap();
    // Prefix-stable and deterministic.
    let mut b = start(&mini(42, 120.0), &PreviewOptions::world_strata_v1());
    while b.steps_taken() < 25 {
        b.step().unwrap();
    }
    assert_eq!(a.strata_stacks(), b.strata_stacks());
    assert_eq!(a.stratigraphy(None).unwrap(), b.stratigraphy(None).unwrap());
}

#[test]
fn world_strata_v1_finish_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::world_strata_v1();
    let cfg = mini(240907, 40.0);
    let batch = simulate_with(&cfg, &opts).unwrap();
    let sd = batch.stratigraphy.as_ref().expect("stratigraphy record");
    sd.validate().unwrap();
    // Published as a sidecar of the terrain artifact (target-dir scratch).
    let root = std::path::PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "world_strata_publish_{}_{}",
        std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()
    ));
    pf_lab::preview::publish(&batch, &root, std::time::Instant::now()).unwrap();
    let mut found = None;
    for e in walk(&root) {
        if e.ends_with("derived/stratigraphy.json") {
            found = Some(e);
        }
    }
    let side: pf_contracts::terrain::StratigraphyDomainV0 =
        serde_json::from_slice(&std::fs::read(found.expect("stratigraphy sidecar")).unwrap()).unwrap();
    assert_eq!(&side, sd);
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

fn walk(d: &std::path::Path) -> Vec<std::path::PathBuf> {
    let mut out = Vec::new();
    for e in std::fs::read_dir(d).unwrap() {
        let p = e.unwrap().path();
        if p.is_dir() {
            out.extend(walk(&p));
        } else {
            out.push(p);
        }
    }
    out.sort();
    out
}
