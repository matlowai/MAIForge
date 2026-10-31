use pf_rockspec::{redbed_v1, v3_1, v3_2c, SurfaceState};
#[test]
fn additive_catalogue_preserves_every_legacy_realization_and_closed_coating_mass() {
    let old = v3_1::catalogue().unwrap();
    let new = redbed_v1::catalogue().unwrap();
    assert_eq!(new.realizations.len(), old.realizations.len() + 1);
    for r in &old.realizations {
        assert_eq!(
            serde_json::to_vec(r).unwrap(),
            serde_json::to_vec(new.realization(&r.id).unwrap()).unwrap()
        );
    }
    let r = new.realization(redbed_v1::REALIZATION).unwrap();
    let sum: f64 = r.spec.mode.iter().map(|m| m.mass_fraction.value).sum();
    assert!((sum - 1.).abs() < 1e-12);
    assert_eq!(
        r.parameters.sedimentary.as_ref().unwrap().grain_coatings[0].mass_fraction,
        0.007
    );
}
#[test]
fn legacy_samples_are_numerically_identical_with_additive_catalogue() {
    let old = v3_1::catalogue().unwrap();
    let new = redbed_v1::catalogue().unwrap();
    for f in v3_1::FixtureSet::swatches().unwrap().samples {
        // Reuses the crate's explicitly manufactured comparison set. No world
        // body is assigned one of these references by this test.
        let a = f.resolve_fixture(&old).unwrap();
        let b = f.resolve_fixture(&new).unwrap();
        let s = SurfaceState::fresh(
            "same-face",
            "same-event",
            [0.031, 0.017, 0.029],
            [0., 0., 1.],
        );
        assert_eq!(
            v3_2c::Evaluator::new(&a)
                .unwrap()
                .surface_sample(&s, 2.)
                .unwrap(),
            v3_2c::Evaluator::new(&b)
                .unwrap()
                .surface_sample(&s, 2.)
                .unwrap()
        );
    }
}
