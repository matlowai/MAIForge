//! WS-CONST (lane WC, 2026-10-12): `LawSetV1` / `UniverseConstantsV1`
//! types, the legacy mapping, validation (typed refusals), serde
//! conventions and the optional recipe pin. ADR 039 §1-§3, acceptance 1
//! and 8, design review C5.

use pf_contracts::artifact::digest_params;
use pf_contracts::laws::{
    constant_spec, legacy, ConstantV1, CrustSeedingLawV1, IsostasyLawV1, SeaLevelInputV1,
    SeaLevelLawV1, SourceClassV1, UnitV1, UniverseError, UniverseV1, CONSTANT_SPECS_V1,
};
use pf_contracts::params::{
    Compensation, ContinentalBodiesBlock, FidelityPreset, MantleBlock, ParamError,
    PlanetParameterProfile, ScenarioPreset, TwoColumnIsostasyBlock,
};

fn toy_airy() -> PlanetParameterProfile {
    let mut m = MantleBlock::default();
    m.compensation = Compensation::Airy;
    ScenarioPreset::toy().profile(FidelityPreset::Preview, m)
}

fn two_column() -> PlanetParameterProfile {
    let mut p = toy_airy();
    p.two_column_isostasy = Some(TwoColumnIsostasyBlock::earth_v0());
    p
}

fn bodies() -> PlanetParameterProfile {
    let mut p = two_column();
    p.continental_bodies = Some(ContinentalBodiesBlock::v0());
    p
}

fn is_class(s: &SourceClassV1, class: &str) -> bool {
    matches!(
        (s, class),
        (SourceClassV1::Cited { .. }, "cited")
            | (SourceClassV1::Derived { .. }, "derived")
            | (
                SourceClassV1::ProjectAssumptionAwaitingReview { .. },
                "assumption"
            )
            | (SourceClassV1::LegacyPresetValue { .. }, "legacy")
    )
}

#[test]
fn ws_const_legacy_profile_reproduces_the_stated_values() {
    let u = toy_airy().universe().unwrap();
    assert!(u.law_set.is_legacy());
    match &u.law_set.isostasy {
        IsostasyLawV1::SingleColumnReferenceV0 {
            reference_crust_thickness_m,
            compensation,
        } => {
            assert_eq!(reference_crust_thickness_m.value, 35_000.0);
            assert_eq!(reference_crust_thickness_m.unit, UnitV1::Meter);
            assert!(is_class(&reference_crust_thickness_m.source, "legacy"));
            assert_eq!(*compensation, Compensation::Airy);
        }
        other => panic!("{other:?}"),
    }
    let k = &u.constants;
    // Stated gravity and 0 m sea level (ADR 039 §2 "legacy presets keep
    // their meaning").
    assert_eq!(k.surface_gravity_m_s2.value, 9.81);
    assert!(is_class(&k.surface_gravity_m_s2.source, "legacy"));
    assert_eq!(k.sea_level_datum_m.value, 0.0);
    assert_eq!(
        u.sea_level_input().unwrap(),
        SeaLevelInputV1::Prescribed { sea_level_m: 0.0 }
    );
    assert!(k.ocean_water_volume_m3.is_none());
    assert!(k.oceanic_crust_density_kg_m3.is_none());
    assert_eq!(k.planet_radius_m.value, 6.371e6);
    assert_eq!(k.mantle_density_kg_m3.value, 3300.0);
    assert_eq!(k.continental_crust_density_kg_m3.value, 2800.0);
    assert_eq!(u.sediment_density_kg_m3(), 2800.0);
    assert_eq!(
        k.seawater_density_kg_m3.value,
        legacy::SEAWATER_DENSITY_KG_M3
    );
    assert_eq!(k.continental_crust_thickness_m.value, 35_000.0);
    assert!(is_class(&k.continental_crust_thickness_m.source, "derived"));
    assert_eq!(k.oceanic_crust_thickness_m.value, 7_000.0);
    assert_eq!(k.stream_power_erodibility_per_yr.value, 2.0e-6);
    assert_eq!(k.hillslope_diffusivity_m2_per_yr.value, 0.1);
    // Contract-level presets (Flexure default) resolve too.
    for s in [
        ScenarioPreset::toy(),
        ScenarioPreset::benchmark_a(),
        ScenarioPreset::benchmark_b(),
        ScenarioPreset::benchmark_c(),
    ] {
        let u = s
            .profile(FidelityPreset::Standard, MantleBlock::default())
            .universe()
            .unwrap();
        assert!(matches!(
            u.law_set.isostasy,
            IsostasyLawV1::SingleColumnReferenceV0 {
                compensation: Compensation::Flexure,
                ..
            }
        ));
    }
}

#[test]
fn ws_const_two_column_makes_ocean_volume_the_sea_level_input() {
    let u = two_column().universe().unwrap();
    assert_eq!(
        u.law_set.sea_level,
        SeaLevelLawV1::VolumePrimaryGlobalFillV0
    );
    let v = u.constants.ocean_water_volume_m3.as_ref().unwrap();
    assert_eq!(v.value, 1.335e18);
    assert_eq!(v.unit, UnitV1::CubicMeter);
    // design review C5: recalled, unchecked -> project assumption awaiting review.
    assert!(is_class(&v.source, "assumption"));
    let rho_oc = u.constants.oceanic_crust_density_kg_m3.as_ref().unwrap();
    assert_eq!(rho_oc.value, 2900.0);
    assert!(is_class(&rho_oc.source, "assumption"));
    // The sea-level constant is only the output datum here.
    assert!(is_class(&u.constants.sea_level_datum_m.source, "derived"));
    assert_eq!(
        u.sea_level_input().unwrap(),
        SeaLevelInputV1::OceanVolume {
            volume_m3: 1.335e18,
            output_datum_m: 0.0
        }
    );
    match &u.law_set.isostasy {
        IsostasyLawV1::TwoColumnCommonDepthV0 {
            fixed_compensation_depth_m,
            derived_depth_margin_m,
        } => {
            assert_eq!(fixed_compensation_depth_m.value, 100_000.0);
            assert_eq!(derived_depth_margin_m.as_ref().unwrap().value, 10_000.0);
        }
        other => panic!("{other:?}"),
    }
    // A declared (non-default) volume is labelled as declared, uncited.
    let mut p = two_column();
    p.two_column_isostasy
        .as_mut()
        .unwrap()
        .ocean_water_volume_m3 = 1.0e18;
    let v = p
        .universe()
        .unwrap()
        .constants
        .ocean_water_volume_m3
        .unwrap();
    assert_eq!(v.value, 1.0e18);
    match v.source {
        SourceClassV1::ProjectAssumptionAwaitingReview { basis } => {
            assert!(basis.contains("declared"), "{basis}")
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn ws_const_bodies_area_fraction_is_derived_from_ocean_fraction() {
    let u = bodies().universe().unwrap();
    assert!(!u.law_set.is_legacy());
    match &u.law_set.crust_seeding {
        CrustSeedingLawV1::GeodesicDiscBodiesV0 {
            continental_area_fraction,
            rifted_margin_width_m,
        } => {
            assert_eq!(continental_area_fraction.value, 1.0 - 0.62);
            assert!(is_class(&continental_area_fraction.source, "derived"));
            assert_eq!(rifted_margin_width_m.value, 200_000.0);
            assert!(is_class(&rifted_margin_width_m.source, "assumption"));
        }
        other => panic!("{other:?}"),
    }
    // Bodies under the legacy isostasy are a valid combination too.
    let mut p = toy_airy();
    p.continental_bodies = Some(ContinentalBodiesBlock::v0());
    let u = p.universe().unwrap();
    assert!(!u.law_set.is_legacy());
    assert_eq!(u.law_set.sea_level, SeaLevelLawV1::PrescribedDatumV0);
}

#[test]
fn ws_const_serde_round_trips_and_names_law_ids() {
    for p in [toy_airy(), two_column(), bodies()] {
        let u = p.universe().unwrap();
        let j = serde_json::to_string(&u).unwrap();
        let back: UniverseV1 = serde_json::from_str(&j).unwrap();
        assert_eq!(back, u);
        back.validate().unwrap();
        assert!(j.contains("\"schema\":\"pf.universe/v1\""), "{j}");
        assert!(j.contains("pf.law-set/v1") && j.contains("pf.universe-constants/v1"));
        assert_eq!(back.digest().unwrap(), u.digest().unwrap());
    }
    let j = serde_json::to_string(&toy_airy().universe().unwrap()).unwrap();
    assert!(j.contains("\"single-column-reference-v0\""), "{j}");
    assert!(j.contains("\"prescribed-datum-v0\"") && j.contains("\"whole-plate-continents-v0\""));
    assert!(j.contains("\"gdh1-v0\"") && j.contains("\"stream-power-n1-linear-diffusion-v0\""));
    assert!(j.contains("\"unit\":\"kg/m3\"") && j.contains("\"class\":\"legacy_preset_value\""));
    let j = serde_json::to_string(&bodies().universe().unwrap()).unwrap();
    assert!(
        j.contains("\"two-column-common-depth-v0\"")
            && j.contains("\"volume-primary-global-fill-v0\"")
    );
    assert!(j.contains("\"geodesic-disc-bodies-v0\"") && j.contains("\"ocean_water_volume_m3\""));
}

#[test]
fn ws_const_unknown_or_malformed_data_is_refused_at_load() {
    let j = serde_json::to_string(&two_column().universe().unwrap()).unwrap();
    // An unknown law id can never become executable.
    for bad in [
        j.replace("\"two-column-common-depth-v0\"", "\"light-near-sound-v0\""),
        j.replace(
            "\"volume-primary-global-fill-v0\"",
            "\"fraction-targeted-v0\"",
        ),
        j.replace("\"gdh1-v0\"", "\"gdh2-v0\""),
        // Unknown keys: in a law entry, a constant, the constants, the root.
        j.replace("\"fixed_compensation_depth_m\"", "\"fixed_depth\""),
        j.replacen("\"value\"", "\"val\"", 1),
        j.replace("\"seawater_density_kg_m3\"", "\"water_density\""),
        j.replace("\"law_set\"", "\"laws\""),
        // Unknown unit string and unknown source class.
        j.replacen("\"unit\":\"m\"", "\"unit\":\"km\"", 1),
        j.replacen(
            "\"class\":\"project_assumption_awaiting_review\"",
            "\"class\":\"vibes\"",
            1,
        ),
    ] {
        assert_ne!(bad, j, "replacement did not apply");
        let e = serde_json::from_str::<UniverseV1>(&bad)
            .expect_err("malformed universe accepted")
            .to_string();
        // Refused for the right reason: the unknown id / key is named.
        assert!(
            e.contains("unknown variant") || e.contains("unknown field"),
            "{e}"
        );
    }
}

#[test]
fn ws_const_validation_refusals_are_typed() {
    let base = two_column().universe().unwrap();
    // Unit mismatch (well-formed JSON, wrong unit for the id).
    let mut u = base.clone();
    u.constants.mantle_density_kg_m3.unit = UnitV1::Meter;
    assert!(matches!(
        u.validate(),
        Err(UniverseError::UnitMismatch {
            constant: "mantle_density_kg_m3",
            ..
        })
    ));
    // Out of envelope, naming the quantity and the range.
    let mut u = base.clone();
    u.constants.mantle_density_kg_m3.value = 7000.0;
    match u.validate() {
        Err(UniverseError::OutOfEnvelope {
            constant, lo, hi, ..
        }) => {
            assert_eq!(constant, "mantle_density_kg_m3");
            assert_eq!((lo, hi), (3000.0, 6000.0));
        }
        other => panic!("{other:?}"),
    }
    let mut u = base.clone();
    u.constants.ocean_water_volume_m3.as_mut().unwrap().value = f64::NAN;
    assert!(matches!(
        u.validate(),
        Err(UniverseError::OutOfEnvelope {
            constant: "ocean_water_volume_m3",
            ..
        })
    ));
    let mut u = base.clone();
    u.constants.oceanic_crust_thickness_m.value = 90_000.0;
    assert!(matches!(
        u.validate(),
        Err(UniverseError::OutOfEnvelope {
            constant: "oceanic_crust_thickness_m",
            ..
        })
    ));
    // Layering.
    let mut u = base.clone();
    u.constants
        .oceanic_crust_density_kg_m3
        .as_mut()
        .unwrap()
        .value = 3290.0;
    u.constants.mantle_density_kg_m3.value = 3200.0;
    assert!(matches!(u.validate(), Err(UniverseError::Inconsistent(_))));
    // Exotic / unimplemented law: named, refused as Unsupported.
    let mut u = base.clone();
    u.law_set.isostasy = IsostasyLawV1::PrattV0;
    assert!(matches!(
        u.validate(),
        Err(UniverseError::Unsupported { .. })
    ));
    // Two inputs for one quantity: a prescribed sea level plus a volume.
    let mut u = toy_airy().universe().unwrap();
    u.constants.ocean_water_volume_m3 = Some(ConstantV1::new(
        1.0e18,
        UnitV1::CubicMeter,
        SourceClassV1::ProjectAssumptionAwaitingReview {
            basis: "test".into(),
        },
    ));
    assert!(matches!(u.validate(), Err(UniverseError::Inconsistent(_))));
    // Volume-primary without a volume.
    let mut u = base.clone();
    u.constants.ocean_water_volume_m3 = None;
    assert!(matches!(u.validate(), Err(UniverseError::Inconsistent(_))));
    // Two-column with a prescribed sea level, single-column with volume
    // primary: not implemented combinations.
    let mut u = base.clone();
    u.law_set.sea_level = SeaLevelLawV1::PrescribedDatumV0;
    assert!(matches!(
        u.validate(),
        Err(UniverseError::Unsupported { .. })
    ));
    let mut u = toy_airy().universe().unwrap();
    u.law_set.sea_level = SeaLevelLawV1::VolumePrimaryGlobalFillV0;
    assert!(matches!(
        u.validate(),
        Err(UniverseError::Unsupported { .. })
    ));
    // Wrong schema version.
    let mut u = base.clone();
    u.law_set.schema = "pf.law-set/v2".into();
    assert!(matches!(u.validate(), Err(UniverseError::Schema { .. })));
    // Profile level: two-column + flexure is refused as Unsupported.
    let mut p = two_column();
    p.mantle.compensation = Compensation::Flexure;
    assert!(matches!(
        p.universe(),
        Err(UniverseError::Unsupported { .. })
    ));
    // Every spec id is unique and every bound finite.
    for (i, s) in CONSTANT_SPECS_V1.iter().enumerate() {
        assert!(
            s.lo.is_finite() && s.hi.is_finite() && s.lo <= s.hi,
            "{}",
            s.id
        );
        assert_eq!(constant_spec(s.id).id, s.id);
        assert!(
            CONSTANT_SPECS_V1[i + 1..].iter().all(|t| t.id != s.id),
            "{}",
            s.id
        );
    }
}

#[test]
fn ws_const_pin_is_opt_in_digest_neutral_when_absent_and_refuses_drift() {
    for p in [toy_airy(), two_column(), bodies()] {
        let json = serde_json::to_string(&p).unwrap();
        assert!(
            !json.contains("pinned_universe"),
            "absent pin leaked into JSON"
        );
        let d0 = digest_params(&p).unwrap();
        let mut q = p.clone();
        q.pin_universe().unwrap();
        q.validate().unwrap();
        let qj = serde_json::to_string(&q).unwrap();
        assert!(qj.contains("\"pinned_universe\""));
        let back: PlanetParameterProfile = serde_json::from_str(&qj).unwrap();
        assert_eq!(back, q);
        // A pinned recipe is a new recipe.
        assert_ne!(digest_params(&q).unwrap(), d0);
        // Same physics as the unpinned resolution.
        assert_eq!(q.universe().unwrap(), p.universe().unwrap());
        // Re-pinning is idempotent.
        let mut r = q.clone();
        r.pin_universe().unwrap();
        assert_eq!(r, q);
    }
    // Drift: the recipe's blocks now resolve to other physics than pinned.
    let mut q = two_column();
    q.pin_universe().unwrap();
    q.two_column_isostasy
        .as_mut()
        .unwrap()
        .ocean_water_volume_m3 = 1.2e18;
    match q.universe() {
        Err(UniverseError::PinMismatch { path, .. }) => {
            assert_eq!(path, "/constants/ocean_water_volume_m3/value")
        }
        other => panic!("{other:?}"),
    }
    assert!(matches!(
        q.validate(),
        Err(ParamError::Universe(UniverseError::PinMismatch { .. }))
    ));
    // A law change is drift too.
    let mut q = toy_airy();
    q.pin_universe().unwrap();
    q.continental_bodies = Some(ContinentalBodiesBlock::v0());
    assert!(matches!(
        q.universe(),
        Err(UniverseError::PinMismatch { .. })
    ));
    // A tampered pin (edited value) is drift as well.
    let mut q = toy_airy();
    q.pin_universe().unwrap();
    q.pinned_universe
        .as_mut()
        .unwrap()
        .constants
        .surface_gravity_m_s2
        .value = 3.71;
    assert!(matches!(
        q.universe(),
        Err(UniverseError::PinMismatch { .. })
    ));
    // An invalid pin is refused with its own typed error before comparison.
    let mut q = toy_airy();
    q.pin_universe().unwrap();
    q.pinned_universe.as_mut().unwrap().law_set.isostasy = IsostasyLawV1::PrattV0;
    assert!(matches!(
        q.universe(),
        Err(UniverseError::Unsupported { .. })
    ));
    // Source labels are NOT physics: re-labelling a pinned constant (e.g. a
    // citation later checked) does not invalidate the pinned world.
    let mut q = two_column();
    q.pin_universe().unwrap();
    q.pinned_universe
        .as_mut()
        .unwrap()
        .constants
        .ocean_water_volume_m3
        .as_mut()
        .unwrap()
        .source = SourceClassV1::Cited {
        citation: "a checked citation".into(),
    };
    q.universe().unwrap();
    q.validate().unwrap();
}

#[test]
fn ws_const_existing_pinned_params_digest_is_unchanged() {
    // The digest pinned by params.rs typed_runoff_profile_digest_is_pinned
    // (toy / preview, contract-level MantleBlock default).
    let p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
    assert_eq!(
        digest_params(&p).unwrap(),
        "blake3:9c112dfe000114c9ae36d76e4de8ce7c8695bc174061d0a378337a0fdd86b543"
    );
    // And the pf-lab CLI toy / preview profile digest recorded in every
    // pre-WS-CONST artifact (seed-robust runs/after, params_digest).
    assert_eq!(
        digest_params(&toy_airy()).unwrap(),
        "blake3:5ce5cfc8eaa0ba74d785746ca4ca78a98bc66ea144fba44cc0e80a64e5ee4767"
    );
}
