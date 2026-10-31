//! SEED-RANGES (lane RG, 2026-10-13): `PreviewOptions::cycle_ranges_v1`
//! (`cycle_v1` + `contact-length-indentation-taper-v1` +
//! `plateau-cap-spill-v1`) on the mini grid (36 x 72, 2 My steps):
//! - the law set names the two new laws only in the new profiles; cycle_v1
//!   serializes without them (byte-identity of existing law sets);
//! - indentations start at sutures, keep the contact converging, and end
//!   by time or budget; the parcel ledger closes every step and the plateau
//!   law conserves volume up to its named delamination;
//! - the genesis module replays the batch bit for bit and a checkpoint /
//!   restore / continue equals the uninterrupted run (the indentation state
//!   is checkpointed).

use pf_contracts::laws::{OrogenPlateauLawV1, PlateSutureLawV1};
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::sim_module::{CheckpointV1, ModuleRun, RunBinding};
use pf_contracts::terrain::CrustKindCodeV0;
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
fn law_set_names_the_new_laws_only_in_the_new_profiles() {
    let cfg = mini(42, 40.0);
    let old = cfg.profile_with(&PreviewOptions::cycle_v1()).unwrap().universe().unwrap();
    assert!(matches!(old.law_set.plate_suture, PlateSutureLawV1::ContactLengthPoleShareV1 { .. }));
    assert!(old.law_set.orogen_plateau.is_legacy());
    let oj = serde_json::to_string(&old).unwrap();
    assert!(!oj.contains("orogen_plateau") && !oj.contains("indentation"));
    let pj = serde_json::to_string(&cfg.profile_with(&PreviewOptions::cycle_v1()).unwrap()).unwrap();
    assert!(!pj.contains("plateau") && !pj.contains("indentation"));
    for opts in [PreviewOptions::cycle_ranges_v1(), PreviewOptions::cycle_ranges_uplands_v1()] {
        let u = cfg.profile_with(&opts).unwrap().universe().unwrap();
        assert!(matches!(
            u.law_set.plate_suture,
            PlateSutureLawV1::ContactLengthIndentationTaperV1 { .. }
        ));
        assert!(matches!(u.law_set.orogen_plateau, OrogenPlateauLawV1::PlateauCapSpillV1 { .. }));
        // Round trip.
        let j = serde_json::to_string(&u).unwrap();
        let back: pf_contracts::laws::UniverseV1 = serde_json::from_str(&j).unwrap();
        assert_eq!(serde_json::to_string(&back).unwrap(), j);
    }
}

#[test]
fn indentation_keeps_sutures_converging_and_ledgers_close() {
    let cfg = mini(42, 200.0);
    let params = cfg.profile_with(&PreviewOptions::cycle_ranges_v1()).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let (mut started, mut converging_steps, mut moved) = (0u32, 0u32, 0.0f64);
    while run.steps_taken() < run.steps_scheduled() {
        let r = run.step().unwrap();
        let c = r.conserve.expect("conserve ledger");
        assert!(c.volume_residual_m3.abs() <= 1e-9 * c.parcel_volume_m3.max(1.0), "{c:?}");
        let g = r.ranges.expect("ranges log under cycle_ranges_v1");
        started += g.indent_started;
        if g.indent_active > 0 && g.indent_mean_rate_m_per_my > 0.0 {
            converging_steps += 1;
        }
        moved += g.plateau_moved_m3;
        // The plateau pass conserves the spilled excess up to its named
        // delamination.
        assert!(
            (g.plateau_excess_m3 - g.plateau_moved_m3 - g.plateau_delaminated_m3).abs()
                <= 1e-9 * g.plateau_excess_m3.max(1.0),
            "{g:?}"
        );
        // No continental cell leaves the advance above the cap.
        assert!(r.advance.max_thickness_m <= 65_000.0 + 1e-6 || {
            let st = run.state();
            (0..st.cells()).all(|i| st.kind[i] != CrustKindCodeV0::Continental || st.thickness_m[i] <= 65_000.0 + 1e-6)
        });
    }
    assert!(started >= 1, "at least one suture (indentation) by 200 My");
    assert!(converging_steps >= 2, "the contact keeps converging after the suture");
    let cs = run.conserve_state().unwrap();
    let ind = cs.indent.as_ref().expect("indentation state");
    assert_eq!(ind.len() as u32, started);
    for i in ind {
        assert!(i.end_cause <= 3);
        if i.end_cause == 0 {
            assert!(run.state().t_my - i.t0_my < 40.0 + 1e-9);
        }
    }
    let _ = moved;
}

#[test]
fn cycle_ranges_v1_module_replay_and_checkpoint_restore() {
    let opts = PreviewOptions::cycle_ranges_v1();
    let cfg = mini(240907, 120.0);
    let batch = simulate_with(&cfg, &opts).unwrap();
    let recipe = preview_recipe(&cfg, &opts, false).unwrap();
    let mut full = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    full.run_to_end(None).unwrap();
    let full_cp = full.checkpoint().unwrap().digest().unwrap();
    let full_out = M::finish(full.into_state()).unwrap();
    assert_eq!(outputs_digest(&full_out).unwrap(), outputs_digest(&batch).unwrap());

    let mut a = ModuleRun::<M>::start(recipe.clone(), RunBinding::GENESIS_SIM_TIME).unwrap();
    a.run_until(37, None).unwrap();
    let bytes = a.checkpoint().unwrap().to_bytes().unwrap();
    drop(a);
    let cp = CheckpointV1::from_bytes(&bytes).unwrap();
    let mut b = ModuleRun::<M>::resume(recipe, RunBinding::GENESIS_SIM_TIME, &cp).unwrap();
    b.run_to_end(None).unwrap();
    assert_eq!(b.checkpoint().unwrap().digest().unwrap(), full_cp);
    let resumed = M::finish(b.into_state()).unwrap();
    assert_eq!(outputs_digest(&resumed).unwrap(), outputs_digest(&full_out).unwrap());
}
