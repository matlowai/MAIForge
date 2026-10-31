use pf_contracts::laws::MarineFinesLawV1;
use pf_contracts::params::{
    FidelityPreset, FormationRecordingV1, MantleBlock, MarineFinesSelectionV1,
    NativeProducerEvidenceV1, ScenarioPreset, NATIVE_PRODUCER_EVIDENCE_V1,
};

#[test]
fn evidence_controls_are_explicit_versioned_and_pinned() {
    let mut profile =
        ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
    let legacy = serde_json::to_value(&profile).unwrap();
    assert!(legacy.get("native_evidence").is_none());
    profile.native_evidence = Some(NativeProducerEvidenceV1 {
        schema: NATIVE_PRODUCER_EVIDENCE_V1.into(),
        formation_age: Some(FormationRecordingV1 {
            event_limit: 10_000,
        }),
        marine_fines: None,
    });
    profile.validate().unwrap();
    // Recording is an observation; it does not replace the selected physics.
    assert!(profile.universe().unwrap().law_set.marine_fines.is_none());
    let c = profile.native_evidence.as_mut().unwrap();
    c.schema = "unknown-config-v2".into();
    assert!(c.validate().is_err());
    c.schema = NATIVE_PRODUCER_EVIDENCE_V1.into();
    c.formation_age.as_mut().unwrap().event_limit = 0;
    assert!(c.validate().is_err());
    c.formation_age.as_mut().unwrap().event_limit = 10_000;
    c.marine_fines = Some(MarineFinesSelectionV1 {
        law: MarineFinesLawV1::ConnectedOceanMixedTerrigenousV1 {
            fraction: 0.25,
            min_depth_m: 1000.0,
        },
        receipt_limit: 30,
    });
    // Missing physical retention cannot be silently replaced by another solver.
    assert!(profile.validate().is_err());
    profile.orogeny = Some(pf_contracts::params::OrogenyBlock::world_strata_v2());
    profile.mantle.compensation = pf_contracts::params::Compensation::Airy;
    profile.two_column_isostasy = Some(pf_contracts::params::TwoColumnIsostasyBlock::earth_v0());
    profile.continental_bodies = Some(pf_contracts::params::ContinentalBodiesBlock::cratons_v1());
    profile.pin_universe().unwrap();
    profile.validate().unwrap();
    profile
        .native_evidence
        .as_mut()
        .unwrap()
        .marine_fines
        .as_mut()
        .unwrap()
        .law = MarineFinesLawV1::ConnectedOceanMixedTerrigenousV1 {
        fraction: 0.5,
        min_depth_m: 1000.0,
    };
    assert!(
        profile.validate().is_err(),
        "a pinned recipe must reject changed fines physics"
    );
}
