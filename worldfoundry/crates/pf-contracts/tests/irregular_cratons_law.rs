//! SEED-IRREGULAR-SEEDS (lane IS, 2026-10-13): the `craton_seeding` law
//! family is opt-in (absent block -> legacy law, nothing serialized), the
//! cratons_v1 bodies block resolves to `irregular-cratons-v1`, round-trips
//! through JSON, and out-of-envelope constants are refused.

use pf_contracts::laws::CratonSeedingLawV1;
use pf_contracts::params::{
    Compensation, ContinentalBodiesBlock, FidelityPreset, MantleBlock, PlanetParameterProfile, ScenarioPreset,
    TwoColumnIsostasyBlock,
};

fn base() -> PlanetParameterProfile {
    let mut p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
    p.mantle.compensation = Compensation::Airy;
    p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
    p
}

#[test]
fn craton_law_is_opt_in_and_round_trips() {
    let mut p = base();
    p.continental_bodies = Some(ContinentalBodiesBlock::v0());
    let j = serde_json::to_string(&p).unwrap();
    assert!(!j.contains("irregular_cratons"));
    let u = p.universe().unwrap();
    assert!(u.law_set.craton_seeding.is_legacy());
    assert!(!serde_json::to_string(&u).unwrap().contains("craton"));

    let mut q = base();
    q.continental_bodies = Some(ContinentalBodiesBlock::cratons_v1());
    let uq = q.universe().unwrap();
    match &uq.law_set.craton_seeding {
        CratonSeedingLawV1::IrregularCratonsV1 { craton_interior_thickness_m, proximal_margin_stretch, .. } => {
            assert_eq!(craton_interior_thickness_m.value, 29_600.0);
            assert_eq!(proximal_margin_stretch.value, 1.15);
        }
        other => panic!("{other:?}"),
    }
    let qj = serde_json::to_string(&q).unwrap();
    let back: PlanetParameterProfile = serde_json::from_str(&qj).unwrap();
    assert_eq!(back, q);
    let uj = serde_json::to_string(&uq).unwrap();
    assert!(uj.contains("irregular-cratons-v1"));
}

#[test]
fn craton_constants_out_of_envelope_are_refused() {
    let mut q = base();
    let mut b = ContinentalBodiesBlock::cratons_v1();
    b.irregular_cratons.as_mut().unwrap().proximal_margin_stretch = 0.5;
    q.continental_bodies = Some(b);
    let e = q.universe().unwrap_err().to_string();
    assert!(e.contains("proximal_margin_stretch"), "{e}");
}
