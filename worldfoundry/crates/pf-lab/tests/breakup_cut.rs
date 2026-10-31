//! SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13): `PreviewOptions::cycle_v2`
//! (`oldstart_retain_v1` + `plate-reorg-scripted-v2`) on the mini grid
//! (36 x 72, 2 My steps). Pre-registered in
//! `generated-test-data/PREREG.txt`.
//!
//! Prefix test ACROSS A BREAKUP STEP (design review review 2026-10-13 item 3.3): the
//! supercontinent share and tenure are lowered (test-only values) so a
//! scripted breakup, with its plate-id reassignment, happens inside the
//! first 40 My; a 40 My run must equal the first 40 My of a 100 My run bit
//! for bit (cells, plate records, parcels, suture and reorg state).

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_lab::preview::{PreviewConfig, PreviewEnv, PreviewOptions, PreviewRun};

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

fn forced_breakup() -> PreviewOptions {
    let mut o = PreviewOptions::cycle_v2();
    let r = o
        .orogeny
        .as_mut()
        .and_then(|g| g.rate_limited_belt.as_mut())
        .and_then(|b| b.plate_reorg.as_mut())
        .expect("cycle_v2 carries the reorg block");
    assert_eq!(r.law_version, 2);
    r.supercontinent_share = 0.3;
    r.supercontinent_tenure_my = 10.0;
    o
}

fn run_to(cfg: &PreviewConfig, stop_steps: u32) -> (PreviewRun, u32) {
    let params = cfg.profile_with(&forced_breakup()).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let mut breakups = 0;
    while run.steps_taken() < stop_steps {
        let r = run.step().unwrap();
        let c = r.conserve.expect("conserve ledger under cycle_v2");
        assert!(c.volume_residual_m3.abs() <= 1e-9 * c.parcel_volume_m3.max(1.0), "{c:?}");
        breakups += r.reorg.expect("reorg report under cycle_v2").breakups;
    }
    (run, breakups)
}

#[test]
fn cycle_v2_prefix_holds_across_a_breakup_step() {
    let (a, breakups) = run_to(&mini(42, 40.0), 20);
    assert!(breakups >= 1, "a scripted breakup inside the first 40 My");
    let (b, _) = run_to(&mini(42, 100.0), 20);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(a.state().plate, b.state().plate);
    assert_eq!(bits(&a.state().thickness_m), bits(&b.state().thickness_m));
    assert_eq!(bits(a.elevation_m()), bits(b.elevation_m()));
    assert_eq!(a.plates(), b.plates());
    assert_eq!(a.conserve_state(), b.conserve_state());
    let ev = &a.conserve_state().unwrap().reorg.as_ref().unwrap().events;
    for e in ev.iter().filter(|e| e.kind == "breakup" || e.kind == "cascade") {
        assert!((0.15..=0.85).contains(&e.side_b_share), "peel accepted: {e:?}");
    }
}
