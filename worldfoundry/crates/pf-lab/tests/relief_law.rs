//! SEED-RELIEF-LAW (lane RL, 2026-10-12): `stream-power-subgrid-nonlinear-v1`
//! in the preview chain, on the off-ladder mini grid (36 x 72, 20 My in
//! 2 My steps). Full-size evidence comes from the ignored driver
//! `relief_measure.rs` and the lane folder `generated-test-data/`.
//!
//! - selection: absent block adds no key; the candidate profile selects the
//!   law id; both erosion blocks together, erosion-deposition with it, and
//!   out-of-envelope constants are refused;
//! - falsifier F3: with n = 1 the new law reproduces lane OG's
//!   `stream-power-subgrid-v1` bit for bit on every state field;
//! - determinism (3 runs), crust closure, F4 realisation (D = 0), and the
//!   genesis module checkpoint / restore / continue = uninterrupted = batch.

use pf_contracts::params::{
    FidelityPreset, NonlinearSubgridErosionBlock, OrogenyBlock, ScenarioPreset, SedimentLawsBlock,
};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding};
use pf_lab::preview::{simulate_with, PreviewConfig, PreviewOptions, PreviewOutputs};
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

fn relief_with(n: f64, k_ref: f64) -> PreviewOptions {
    PreviewOptions {
        orogeny: Some(OrogenyBlock {
            rate_limited_belt: OrogenyBlock::v1().rate_limited_belt,
            subgrid_erosion: None,
            nonlinear_subgrid_erosion: Some(NonlinearSubgridErosionBlock {
                slope_exponent_n: n,
                reference_steepness_m: k_ref,
                ..NonlinearSubgridErosionBlock::v1()
            }),
            critical_slope_erosion: None,
            sediment_retention: None,
            // lane IT (2026-10-20): later opt-in OrogenyBlock fields stay absent.
            ..Default::default()
        }),
        ..PreviewOptions::crust_bodies_v0()
    }
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn relief_selection_keys_and_refusals() {
    let cfg = mini(7);
    // Absent: no key anywhere (orogen_v1 serializes as before).
    let q = cfg.profile_with(&PreviewOptions::orogen_v1()).unwrap();
    let json = serde_json::to_string(&q).unwrap();
    assert!(!json.contains("nonlinear_subgrid_erosion"), "{json}");
    let laws = serde_json::to_string(&q.universe().unwrap().law_set).unwrap();
    assert!(!laws.contains("nonlinear"), "{laws}");
    // Present: the law id and its constants.
    let r = cfg
        .profile_with(&PreviewOptions::orogen_relief_v1())
        .unwrap();
    let u = r.universe().unwrap();
    let ids = pf_contracts::sim_module::law_ids(&u.law_set).unwrap();
    assert_eq!(ids["erosion"], "stream-power-subgrid-nonlinear-v1");
    assert_eq!(ids["tectonic_thickening"], "rate-limited-belt-v1");
    let mut pinned = r.clone();
    pinned.pin_universe().unwrap();
    let back: pf_contracts::params::PlanetParameterProfile =
        serde_json::from_str(&serde_json::to_string(&pinned).unwrap()).unwrap();
    assert_eq!(back, pinned);
    back.validate().unwrap();
    // Both erosion blocks: refused.
    let mut o = OrogenyBlock::v1();
    o.nonlinear_subgrid_erosion = Some(NonlinearSubgridErosionBlock::v1());
    let both = PreviewOptions {
        orogeny: Some(o),
        ..PreviewOptions::crust_bodies_v0()
    };
    let e = cfg.profile_with(&both).unwrap_err().to_string();
    assert!(e.contains("nonlinear_subgrid_erosion"), "{e}");
    // Envelope refusals.
    for (n, k) in [(0.5, 500.0), (5.0, 500.0), (2.0, 0.0), (2.0, 1.0e6)] {
        let e = cfg
            .profile_with(&relief_with(n, k))
            .unwrap_err()
            .to_string();
        assert!(
            e.contains("outside the v1 validity envelope"),
            "n {n} k {k}: {e}"
        );
    }
    // Erosion-deposition with it: refused (as with lane OG's subgrid law).
    let ed = PreviewOptions {
        sediment_laws: Some(SedimentLawsBlock::v1()),
        ..PreviewOptions::orogen_relief_v1()
    };
    assert!(cfg.profile_with(&ed).is_err());
}

fn assert_state_fields_equal(a: &PreviewOutputs, b: &PreviewOutputs, what: &str) {
    assert_eq!(
        bits(&a.elevation_m),
        bits(&b.elevation_m),
        "{what}: elevation"
    );
    assert_eq!(
        bits(&a.state.thickness_m),
        bits(&b.state.thickness_m),
        "{what}: thickness"
    );
    assert_eq!(
        bits(&a.state.sediment_m),
        bits(&b.state.sediment_m),
        "{what}: sediment"
    );
    assert_eq!(
        bits(&a.state.eroded_total_m),
        bits(&b.state.eroded_total_m),
        "{what}: eroded"
    );
    let ra = a
        .orogen
        .as_ref()
        .unwrap()
        .subgrid_relief_m
        .as_ref()
        .unwrap();
    let rb = b
        .orogen
        .as_ref()
        .unwrap()
        .subgrid_relief_m
        .as_ref()
        .unwrap();
    assert_eq!(bits(ra), bits(rb), "{what}: subgrid relief");
    let sa = &a.orogen.as_ref().unwrap().steps;
    let sb = &b.orogen.as_ref().unwrap().steps;
    assert_eq!(sa.len(), sb.len());
    for (x, y) in sa.iter().zip(sb) {
        assert_eq!(
            x.subgrid_target_m3.to_bits(),
            y.subgrid_target_m3.to_bits(),
            "{what}"
        );
        assert_eq!(
            x.subgrid_realised_m3.to_bits(),
            y.subgrid_realised_m3.to_bits(),
            "{what}"
        );
        assert_eq!(x.subgrid_capped_cells, y.subgrid_capped_cells, "{what}");
        assert_eq!(
            x.subgrid_trunk_raised_cells, y.subgrid_trunk_raised_cells,
            "{what}"
        );
    }
}

/// Falsifier F3: n = 1 is lane OG's law, bit for bit (any k_ref).
#[test]
fn relief_n1_is_stream_power_subgrid_v1_bitwise() {
    for seed in [240907u64, 3, 7] {
        let cfg = mini(seed);
        let og = simulate_with(&cfg, &PreviewOptions::orogen_v1()).unwrap();
        for k_ref in [500.0, 37.0] {
            let rl = simulate_with(&cfg, &relief_with(1.0, k_ref)).unwrap();
            assert_state_fields_equal(&og, &rl, &format!("seed {seed} k_ref {k_ref}"));
        }
    }
}

#[test]
fn relief_n2_deterministic_closing_and_different_from_n1() {
    for seed in [240907u64, 3, 7] {
        let cfg = mini(seed);
        let opts = PreviewOptions::orogen_relief_v1();
        let a = simulate_with(&cfg, &opts).unwrap();
        let da = outputs_digest(&a).unwrap();
        for _ in 0..2 {
            let b = simulate_with(&cfg, &opts).unwrap();
            assert_eq!(outputs_digest(&b).unwrap(), da, "seed {seed}");
        }
        let cv = &a.terrain.crust_volume;
        assert!(
            cv.closure_plug_m3.abs() <= 1e-9 * cv.initial_m3,
            "seed {seed}: closure plug {} of {}",
            cv.closure_plug_m3,
            cv.initial_m3
        );
        let og = a.orogen.as_ref().unwrap();
        assert!(og
            .subgrid_relief_m
            .as_ref()
            .unwrap()
            .iter()
            .all(|r| *r >= 0.0));
        let (t, r): (f64, f64) = og.steps.iter().fold((0.0, 0.0), |(t, r), s| {
            (t + s.subgrid_target_m3, r + s.subgrid_realised_m3)
        });
        assert!(t > 0.0, "seed {seed}: the law must erode");
        assert!(
            (r - t).abs() <= 1e-2 * t,
            "seed {seed}: target {t} realised {r}"
        );
        assert!(a
            .terrain
            .approximations
            .iter()
            .any(|x| x == "two-elevation-subgrid-stream-power-nonlinear-v1"));
        assert!(a
            .terrain
            .assumptions
            .iter()
            .any(|x| x.id == "A-REFERENCE-STEEPNESS"));
        // n = 2 at k_ref = 500 m erodes the (low-relief) mini world less.
        let n1 = simulate_with(&cfg, &PreviewOptions::orogen_v1()).unwrap();
        let eroded = |o: &PreviewOutputs| o.steps.iter().map(|s| s.eroded_m3).sum::<f64>();
        eprintln!(
            "seed {seed}: eroded n2 {:.4e} vs n1 {:.4e} m3",
            eroded(&a),
            eroded(&n1)
        );
        assert!(eroded(&a) < eroded(&n1), "seed {seed}");
    }
}

/// Falsifier F4 (realisation): with hillslope diffusion off the erosion
/// stage removes exactly the n = 2 target.
#[test]
fn relief_target_is_realised_exactly_without_diffusion() {
    for seed in [240907u64, 3] {
        let mut s = pf_contracts::params::SurfaceBlock::default();
        s.hillslope_diffusivity_m2_per_yr = 0.0;
        let opts = PreviewOptions {
            surface_override: Some(s),
            ..PreviewOptions::orogen_relief_v1()
        };
        let a = simulate_with(&mini(seed), &opts).unwrap();
        let og = a.orogen.as_ref().unwrap();
        let mut worst = 0.0f64;
        for st in &og.steps {
            if st.subgrid_target_m3 > 0.0 {
                worst = worst.max(
                    (st.subgrid_realised_m3 - st.subgrid_target_m3).abs() / st.subgrid_target_m3,
                );
            }
        }
        eprintln!("seed {seed}: D=0 worst relative volume dev {worst:.3e}");
        assert!(worst < 1e-6, "seed {seed}: {worst}");
    }
}

#[test]
fn relief_module_checkpoint_restore_continue_equals_uninterrupted() {
    let cfg = mini(240907);
    let opts = PreviewOptions::orogen_relief_v1();
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    full.run_to_end(None).unwrap();
    let full_cp = full.checkpoint().unwrap().digest().unwrap();
    let full_out = M::finish(full.into_state()).unwrap();
    let batch = simulate_with(&cfg, &opts).unwrap();
    assert_eq!(
        outputs_digest(&full_out).unwrap(),
        outputs_digest(&batch).unwrap(),
        "module replay != batch under the relief law"
    );
    for k in [3u64, 6] {
        let mut a = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
        a.run_until(k, None).unwrap();
        let bytes = a.checkpoint().unwrap().to_bytes().unwrap();
        drop(a);
        let cp = CheckpointV1::from_bytes(&bytes).unwrap();
        let mut b =
            ModuleRun::<M>::resume(recipe.clone(), RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
        b.run_to_end(None).unwrap();
        assert_eq!(b.checkpoint().unwrap().digest().unwrap(), full_cp, "k {k}");
        let out = M::finish(b.into_state()).unwrap();
        assert_eq!(
            outputs_digest(&out).unwrap(),
            outputs_digest(&full_out).unwrap(),
            "k {k}"
        );
    }
}
