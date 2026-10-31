//! SEED-OROGEN (lane OG, 2026-10-12): `rate-limited-belt-v1` and
//! `stream-power-subgrid-v1` in the preview chain, on the off-ladder mini
//! grid (36 x 72, 20 My in 2 My steps). Full-size evidence is produced by
//! the ignored driver in `orogen_measure.rs`.
//!
//! - legacy selection is untouched: `advance_with(None)` equals `advance`
//!   bit for bit over a multi-step history, and an absent orogeny block
//!   leaves the profile JSON and the law set JSON free of new keys;
//! - the new laws: three runs byte-identical; crust ledger closes; belt and
//!   arc volumes are placed exactly; spreading conserves; the subgrid target
//!   is realised by the shared erosion stage;
//! - genesis module: checkpoint at k, restore, continue = uninterrupted;
//! - envelope refusal for an out-of-range law constant.

use pf_contracts::params::{FidelityPreset, OrogenyBlock, ScenarioPreset};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding};
use pf_geodynamics::kinematic::{advance, advance_with, classify_cells, CrustState};
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::uplift::UpliftPrior;
use pf_lab::preview::{simulate_with, PreviewConfig, PreviewOptions};
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
fn orogen_legacy_advance_with_none_is_advance_bitwise() {
    for seed in [240907u64, 3, 7] {
        let cfg = mini(seed);
        let p = cfg.profile().unwrap();
        let grid = cfg.grid(&p).unwrap();
        let seeding = seed_plates(&grid, &p, &cfg.seed_path()).unwrap();
        let prior = UpliftPrior::preview_v0();
        let mut a = CrustState::seed(&grid, &seeding, &prior).unwrap();
        let (class, rate, _) = classify_cells(&grid, &a, &seeding.plates).unwrap();
        a.boundary_class = class;
        a.normal_rate_m_per_my = rate;
        let mut b = a.clone();
        for _ in 0..12 {
            let (na, la, ta) = advance(&grid, &a, &seeding.plates, &prior, 2.0).unwrap();
            let (nb, lb, tb, src) =
                advance_with(&grid, &b, &seeding.plates, &prior, 2.0, None).unwrap();
            assert_eq!(na, nb, "seed {seed}");
            assert_eq!(la, lb, "seed {seed}");
            assert_eq!(
                ta.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                tb.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
            assert_eq!(src.len(), grid.cells());
            a = na;
            b = nb;
        }
    }
}

#[test]
fn orogen_absent_block_adds_no_keys_and_present_block_selects_laws() {
    let cfg = mini(7);
    let p = cfg.profile().unwrap();
    let json = serde_json::to_string(&p).unwrap();
    assert!(!json.contains("orogeny"));
    let u = p.universe().unwrap();
    let laws = serde_json::to_string(&u.law_set).unwrap();
    assert!(!laws.contains("tectonic_thickening"), "{laws}");
    assert!(laws.contains("stream-power-n1-linear-diffusion-v0"));
    assert!(u.law_set.is_legacy());

    let q = cfg.profile_with(&PreviewOptions::orogen_v1()).unwrap();
    let uq = q.universe().unwrap();
    assert!(!uq.law_set.is_legacy());
    let ids = pf_contracts::sim_module::law_ids(&uq.law_set).unwrap();
    assert_eq!(ids["tectonic_thickening"], "rate-limited-belt-v1");
    assert_eq!(ids["erosion"], "stream-power-subgrid-v1");
    // A pinned universe with the new laws round-trips and validates.
    let mut r = q.clone();
    r.pin_universe().unwrap();
    let back: pf_contracts::params::PlanetParameterProfile =
        serde_json::from_str(&serde_json::to_string(&r).unwrap()).unwrap();
    assert_eq!(back, r);
    back.validate().unwrap();
}

#[test]
fn orogen_envelope_refuses_out_of_range_constants() {
    let cfg = mini(7);
    let mut o = OrogenyBlock::v1();
    o.rate_limited_belt
        .as_mut()
        .unwrap()
        .lower_crust_channel_viscosity_pa_s = 1.0e30;
    let opts = PreviewOptions {
        orogeny: Some(o),
        ..PreviewOptions::crust_bodies_v0()
    };
    let e = cfg.profile_with(&opts).unwrap_err().to_string();
    assert!(e.contains("lower_crust_channel_viscosity_pa_s"), "{e}");
    let mut o = OrogenyBlock::v1();
    o.rate_limited_belt
        .as_mut()
        .unwrap()
        .delamination_thickness_m = 30_000.0;
    let opts = PreviewOptions {
        orogeny: Some(o),
        ..PreviewOptions::crust_bodies_v0()
    };
    assert!(cfg.profile_with(&opts).is_err());
}

#[test]
fn orogen_new_laws_deterministic_and_closing() {
    for (seed, opts) in [
        (240907u64, PreviewOptions::orogen_v1()),
        (3, PreviewOptions::orogen_v1()),
        (
            7,
            PreviewOptions {
                orogeny: Some(OrogenyBlock::v1()),
                ..PreviewOptions::legacy()
            },
        ),
    ] {
        let cfg = mini(seed);
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
        let og = a.orogen.as_ref().expect("orogen state");
        assert_eq!(og.steps.len(), 10);
        let mut belt_in = 0.0;
        let mut target = 0.0;
        let mut realised = 0.0;
        for s in &og.steps {
            belt_in += s.belt_in_m3;
            assert!(
                (s.belt_placed_m3 - s.belt_in_m3).abs() <= 1e-9 * s.belt_in_m3.max(1.0),
                "seed {seed} t {}: belt placed {} of {}",
                s.t_my,
                s.belt_placed_m3,
                s.belt_in_m3
            );
            assert!((s.arc_placed_m3 - s.arc_created_m3).abs() <= 1e-9 * s.arc_created_m3.max(1.0));
            assert!(s.max_thickness_m <= 70_000.0 + 1e-6);
            target += s.subgrid_target_m3;
            realised += s.subgrid_realised_m3;
        }
        eprintln!(
            "seed {seed}: belt in {belt_in:.3e} m3, subgrid target {target:.4e} realised {realised:.4e} m3, max dev {:.3e} m",
            og.steps.iter().map(|s| s.subgrid_max_abs_dev_m).fold(0.0, f64::max)
        );
        if target > 0.0 {
            // Realised includes the (tiny) hillslope-diffusion removal.
            assert!((realised - target).abs() <= 1e-2 * target, "seed {seed}");
        }
        assert!(og
            .subgrid_relief_m
            .as_ref()
            .unwrap()
            .iter()
            .all(|r| *r >= 0.0));
    }
}

#[test]
fn orogen_module_checkpoint_restore_continue_equals_uninterrupted() {
    let cfg = mini(240907);
    let opts = PreviewOptions::orogen_v1();
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    full.run_to_end(None).unwrap();
    let full_cp = full.checkpoint().unwrap().digest().unwrap();
    let full_out = M::finish(full.into_state()).unwrap();
    let batch = simulate_with(&cfg, &opts).unwrap();
    assert_eq!(
        outputs_digest(&full_out).unwrap(),
        outputs_digest(&batch).unwrap(),
        "module replay != batch under the new laws"
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

/// Falsifier F4 (injection half): with hillslope diffusion off, the
/// erosion stage removes exactly the subgrid target on every cell (the only
/// difference is floating-point rounding of the effective-area round trip).
#[test]
fn orogen_subgrid_target_is_realised_exactly_without_diffusion() {
    for seed in [240907u64, 3] {
        let mut s = pf_contracts::params::SurfaceBlock::default();
        s.hillslope_diffusivity_m2_per_yr = 0.0;
        let opts = PreviewOptions {
            surface_override: Some(s),
            ..PreviewOptions::orogen_v1()
        };
        let a = simulate_with(&mini(seed), &opts).unwrap();
        let og = a.orogen.as_ref().unwrap();
        let mut worst_rel = 0.0f64;
        for st in &og.steps {
            if st.subgrid_target_m3 > 0.0 {
                worst_rel = worst_rel.max(
                    (st.subgrid_realised_m3 - st.subgrid_target_m3).abs() / st.subgrid_target_m3,
                );
            }
        }
        let dev = og
            .steps
            .iter()
            .map(|s| s.subgrid_max_abs_dev_m)
            .fold(0.0, f64::max);
        eprintln!("seed {seed}: D=0 worst relative volume dev {worst_rel:.3e}, max per-cell dev {dev:.3e} m");
        assert!(worst_rel < 1e-6, "seed {seed}: {worst_rel}");
    }
}

/// Merge guard (lane SD wiring, applied by lane OG): the deposition
/// partition rides in the module checkpoint only under non-legacy surface
/// laws; checkpoint, restore and continue equals the uninterrupted run and
/// the batch, with lane SD's `sediment_v1` and with SD + the belt law.
#[test]
fn orogen_sd_merge_module_checkpoint_with_deposition_laws() {
    let belt_plus_sd = PreviewOptions {
        orogeny: Some(OrogenyBlock {
            subgrid_erosion: None,
            ..OrogenyBlock::v1()
        }),
        ..PreviewOptions::sediment_v1()
    };
    for opts in [PreviewOptions::sediment_v1(), belt_plus_sd] {
        let cfg = mini(3);
        let recipe = preview_recipe(&cfg, &opts, false).unwrap();
        let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
        full.run_to_end(None).unwrap();
        let full_cp = full.checkpoint().unwrap().digest().unwrap();
        let full_out = M::finish(full.into_state()).unwrap();
        assert!(full_out.steps.iter().all(|s| s.deposition.is_some()));
        assert_eq!(
            outputs_digest(&full_out).unwrap(),
            outputs_digest(&simulate_with(&cfg, &opts).unwrap()).unwrap()
        );
        let mut a = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
        a.run_until(5, None).unwrap();
        let cp = CheckpointV1::from_bytes(&a.checkpoint().unwrap().to_bytes().unwrap()).unwrap();
        drop(a);
        let mut b = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
        b.run_to_end(None).unwrap();
        assert_eq!(b.checkpoint().unwrap().digest().unwrap(), full_cp);
    }
    // Erosion-deposition with the subgrid law is refused (no combined law).
    let both = PreviewOptions {
        orogeny: Some(OrogenyBlock::v1()),
        ..PreviewOptions::sediment_v1()
    };
    assert!(mini(3).profile_with(&both).is_err());
}
