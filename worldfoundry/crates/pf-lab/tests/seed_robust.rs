//! SEED-ROBUST (lane R5, 2026-10-12): the full preview chain (tectonics,
//! routing, erosion, ledgers) runs through a step that used to refuse.
//!
//! Toy preset, seed 3, L0 (144 x 288), legacy crust: before this node
//! `pf-lab simulate` refused at the 16 -> 18 My step with "arc overlap 11334:
//! no input interface for plates 3 / 6; event-time geometry required". Plate
//! 1's tongue, one cell wide, separated plates 3 and 6; both reached its tip
//! in that step, so their contact formed during the step. The kernel now
//! samples that new contact across the material that met.

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_lab::preview::{simulate, PreviewConfig};

fn cfg() -> PreviewConfig {
    PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: None,
        seed: 3,
        duration_my: Some(18.0),
        step_my: None,
    }
}

#[test]
fn seed3_l0_runs_through_the_former_refusal_and_closes_its_ledgers() {
    let out = simulate(&cfg()).expect("the 16 -> 18 My step no longer refuses");
    assert_eq!(out.state.t_my, 18.0);
    // The new contact is published with its source pair (plate 3's source
    // 11335 first, then plate 6's 11332; not grid neighbours, by design).
    let epoch = out
        .exchange_epochs
        .iter()
        .find(|e| e.start_my == 16.0)
        .expect("epoch 16 -> 18");
    let s = epoch
        .arc_normals
        .iter()
        .find(|s| s.overlap_element.0 == 11334)
        .expect("the former refusal is sampled");
    assert_eq!((s.winner.0, s.loser.0), (3, 6));
    assert_eq!(
        [s.input_interface[0].0, s.input_interface[1].0],
        [11335, 11332]
    );
    assert_eq!(s.source_epoch_my, 16.0);
    assert!(s.normal_rate_m_per_my > 0.0 && s.normal_rate_m_per_my.is_finite());
    // Crust volume: whatever the named terms and the measured advection
    // residual do not explain stays at rounding level.
    let cv = &out.terrain.crust_volume;
    assert!(
        cv.closure_plug_m3.abs() <= 1e-9 * cv.initial_m3,
        "closure plug {} of {}",
        cv.closure_plug_m3,
        cv.initial_m3
    );
    // Deterministic: a second run is identical in state, surface and steps.
    let again = simulate(&cfg()).unwrap();
    assert_eq!(out.state, again.state);
    assert_eq!(out.elevation_m, again.elevation_m);
    assert_eq!(out.exchange_epochs, again.exchange_epochs);
}
