use pf_rockspec::{surface_v1::*, *};
use serde_json::json;

fn estimate(value: f64, unit: &str) -> Estimate {
    Estimate {value,interval:[value;2],unit:unit.into(),interval_kind:"declared_scenario_value".into(),
        distribution:"fixed".into(),reference:"surface_v1 test scenario, not world truth".into(),
        evidence_kind:"manufactured_diagnostic".into(),support_scale:"one exposed test patch".into(),
        conditions:"explicit scenario for model qualification; no claim of measured Moab age or dust budget".into(),
        covariance_group:"test_scenario".into()}
}
fn body(id: &str) -> v3::ResolvedBody {
    let c = v3_1::catalogue().unwrap();
    v3_1::SurfaceFixture {
        id: "surface-test".into(),
        realization_id: id.into(),
        history: "fresh_dry".into(),
        weathering_family: "test".into(),
        surface: fresh(),
        footprint_m: 2.,
    }
    .resolve_fixture(&c)
    .unwrap()
}
fn fresh() -> SurfaceState {
    SurfaceState::fresh("test-face", "test-collapse", [0.; 3], [0., 0., 1.])
}
fn input(age: f64, ai: f64, p: f64, t: f64, dust: f64) -> Input {
    Input {
        exposure: Exposure::SinceEvent {
            years: estimate(age, "yr"),
            event_ref: "declared-test-collapse".into(),
        },
        climate: Some(Climate {
            artifact_ref: "test-canonical-field-semantics-not-real-artifact".into(),
            sample_ref: "test-point".into(),
            aridity_index: Some(estimate(ai, "1")),
            temperature_mean_c: estimate(t, "degC"),
            precipitation_mm_per_yr: estimate(p, "mm/yr"),
            stationary_history_assumption: "declared constant forcing over test exposure".into(),
        }),
        varnish_supply: Some(VarnishSupply {
            dry_mass_kg_m2: estimate(dust, "kg/m2"),
            source_ref: "test-external-dust-account".into(),
        }),
    }
}
fn predicted(m: &Model, b: &v3::ResolvedBody, i: &Input) -> (Prediction, OpticalRow) {
    match m.evaluate(b, &fresh(), 2., Some(i)).unwrap() {
        Evaluation::Predicted { prediction, row } => (prediction, row),
        x => panic!("expected prediction, got {x:?}"),
    }
}
fn y(rgb: [f64; 3]) -> f64 {
    rgb.iter()
        .zip([0.2126, 0.7152, 0.0722])
        .map(|(r, w)| r * w)
        .sum()
}

#[test]
fn all_fresh_rows_key_off_are_byte_identical_and_inputs_immutable() {
    let c = v3_1::catalogue().unwrap();
    let m = Model::embedded().unwrap();
    let cases = v3_1::FixtureSet::comparisons(&c).unwrap();
    let mut n = 0;
    for case in cases.samples.iter().filter(|s| s.history == "fresh_dry") {
        let b = case.resolve_fixture(&c).unwrap();
        let before = serde_json::to_vec(&b).unwrap();
        let expected = v3_2c::Evaluator::new(&b)
            .unwrap()
            .surface_sample(&case.surface, case.footprint_m)
            .unwrap();
        let Evaluation::Disabled { sample } = m
            .evaluate(&b, &case.surface, case.footprint_m, None)
            .unwrap()
        else {
            panic!()
        };
        assert_eq!(
            serde_json::to_vec(&expected).unwrap(),
            serde_json::to_vec(&sample).unwrap()
        );
        assert_eq!(before, serde_json::to_vec(&b).unwrap());
        n += 1;
    }
    assert!(n >= 67, "only {n} rows compared");
    println!("KEY_OFF_FRESH_ROWS={n}");
}
#[test]
fn known_zero_age_preserves_fresh_optics_and_unknown_is_not_fresh() {
    let m = Model::embedded().unwrap();
    let b = body("sandstone_quartz.reference");
    let mut i = input(0., 0.2, 200., 15., 1.);
    i.climate = None;
    i.varnish_supply = None;
    let (p, row) = predicted(&m, &b, &i);
    assert_eq!(p.weights.fresh, 1.);
    assert_eq!(
        row.base_reflectance_linear,
        row.substrate.base_reflectance_linear
    );
    assert_eq!(row.roughness.to_bits(), row.substrate.roughness.to_bits());
    assert!(p.weights.lichen.is_none());
    assert!(p.recognized_oxide_fe_mass_fraction > 0.);
    assert!(p.available_fe_kg_m2 > 0.);
    i.exposure = Exposure::Unknown {
        reason: "no T1b exposure event".into(),
    };
    let x = m.evaluate(&b, &fresh(), 2., Some(&i)).unwrap();
    let j = serde_json::to_value(x).unwrap();
    assert_eq!(j["status"], "unknown");
    assert!(j.get("row").is_none());
}
#[test]
fn age_climate_iron_and_external_supply_have_distinct_effects() {
    let m = Model::embedded().unwrap();
    let b = body("sandstone_quartz.reference");
    let (_, fresh) = predicted(&m, &b, &input(0., 0.2, 200., 15., 1.));
    let (red, r) = predicted(&m, &b, &input(10000., 0.2, 200., 15., 0.));
    let (old, v) = predicted(&m, &b, &input(100000., 0.2, 200., 15., 1.));
    let (humid, h) = predicted(&m, &b, &input(10000., 1.5, 1500., 18., 0.));
    assert!(red.weights.hematite > red.weights.rind);
    assert!(r.base_srgb[0] > r.base_srgb[1] && r.base_srgb[1] > r.base_srgb[2]);
    assert!(old.weights.varnish > 0.9);
    assert!(y(v.base_reflectance_linear) < y(r.base_reflectance_linear));
    assert!(y(fresh.base_reflectance_linear) > y(r.base_reflectance_linear));
    assert!(humid.weights.rind > humid.weights.hematite);
    assert_ne!(r.base_reflectance_linear, h.base_reflectance_linear);
    for &(rain, t) in &[(0., 15.), (200., -10.)] {
        let (p, _) = predicted(&m, &b, &input(10000., 0.2, rain, t, 1.));
        assert_eq!(
            p.weights.fresh, 1.,
            "no rain or frozen mean must not generate coating"
        );
    }
    // Chert has no recognized oxide Fe. Dust still may varnish it, independently.
    let nofe = body("chert.white");
    let (p, _) = predicted(&m, &nofe, &input(10000., 0.2, 200., 15., 0.));
    assert_eq!(p.recognized_oxide_fe_mass_fraction, 0.);
    assert_eq!(p.weights.fresh, 1.);
    let (p, _) = predicted(&m, &nofe, &input(100000., 0.2, 200., 15., 1.));
    assert!(p.weights.varnish > 0.9);
}
#[test]
fn deterministic_grid_inventory_optical_energy_and_body_preservation() {
    let m = Model::embedded().unwrap();
    let mut n = 0;
    for id in [
        "sandstone_quartz.reference",
        "arkose.reference",
        "basalt_tholeiitic.reference",
        "chert.white",
    ] {
        let b = body(id);
        let before = serde_json::to_vec(&b).unwrap();
        for age in [0., 1., 1000., 100000., 1e6] {
            for (ai, p, t) in [
                (0., 0., 15.),
                (0.2, 200., 15.),
                (1.5, 1500., 18.),
                (0.2, 200., -10.),
            ] {
                for dust in [0., 0.001, 1.] {
                    let i = input(age, ai, p, t, dust);
                    let (pred, row) = predicted(&m, &b, &i);
                    let twice = predicted(&m, &b, &i);
                    assert_eq!(
                        serde_json::to_vec(&(&pred, &row)).unwrap(),
                        serde_json::to_vec(&twice).unwrap()
                    );
                    let w = &pred.weights;
                    assert!((w.fresh + w.hematite + w.rind + w.varnish - 1.).abs() < 1e-12);
                    assert!([w.fresh, w.hematite, w.rind, w.varnish]
                        .iter()
                        .all(|&x| (0. ..=1.).contains(&x)));
                    assert!(pred.redistributed_fe_kg_m2 <= pred.available_fe_kg_m2 * (1. + 1e-12));
                    assert!(pred.external_varnish_used_kg_m2 <= dust * (1. + 1e-12));
                    let p = m.parameters();
                    let recovered = pred.hematite_thickness_m
                        * p["hematite_density"].value
                        * p["hematite_fe"].value
                        + pred.goethite_rind_thickness_m
                            * p["goethite_density"].value
                            * p["goethite_fe"].value;
                    assert!((recovered - pred.redistributed_fe_kg_m2).abs() < 1e-12);
                    for j in 0..3 {
                        assert!(row.diffuse_energy_linear[j] + row.f0_linear[j] <= 1. + 1e-12);
                    }
                    assert_eq!(before, serde_json::to_vec(&b).unwrap());
                    n += 1;
                }
            }
        }
    }
    println!("DETERMINISTIC_GRID_ROWS={n}");
    assert_eq!(n, 240);
}
#[test]
fn missing_inputs_and_erosion_cannot_silently_become_known_age() {
    let m = Model::embedded().unwrap();
    let b = body("sandstone_quartz.reference");
    let mut i = input(10000., 0.2, 200., 15., 1.);
    i.climate.as_mut().unwrap().aridity_index = None;
    assert!(matches!(
        m.evaluate(&b, &fresh(), 2., Some(&i)).unwrap(),
        Evaluation::Unknown { .. }
    ));
    i = input(10000., 0.2, 200., 15., 1.);
    i.varnish_supply = None;
    assert!(matches!(
        m.evaluate(&b, &fresh(), 2., Some(&i)).unwrap(),
        Evaluation::Unknown { .. }
    ));
    i = input(10000., 0.2, 200., 15., 1.);
    i.exposure = Exposure::ErosionResidence {
        erosion_m_per_yr: estimate(0., "m/yr"),
        renewal_depth_m: estimate(0.1, "m"),
        record_ref: "test".into(),
    };
    assert!(matches!(
        m.evaluate(&b, &fresh(), 2., Some(&i)).unwrap(),
        Evaluation::Unknown { .. }
    ));
    if let Exposure::ErosionResidence {
        erosion_m_per_yr, ..
    } = &mut i.exposure
    {
        *erosion_m_per_yr = estimate(1e-5, "m/yr");
    }
    let (p, _) = predicted(&m, &b, &i);
    assert_eq!(p.exposure_years, 10000.);
    assert!(p.exposure_basis.contains("not_last_event"));
    i.exposure = Exposure::SinceEvent {
        years: estimate(1e8, "yr"),
        event_ref: "test".into(),
    };
    assert!(matches!(
        m.evaluate(&b, &fresh(), 2., Some(&i)).unwrap(),
        Evaluation::Unknown { .. }
    ));
}
#[test]
fn rejects_sentinels_negative_nonfinite_unprovenanced_and_double_coatings() {
    let m = Model::embedded().unwrap();
    let b = body("sandstone_quartz.reference");
    for bad in [-1., f64::NAN, f64::INFINITY] {
        let i = input(1000., bad, 200., 15., 1.);
        assert!(m.evaluate(&b, &fresh(), 2., Some(&i)).is_err());
        let i = input(bad, 0.2, 200., 15., 1.);
        assert!(m.evaluate(&b, &fresh(), 2., Some(&i)).is_err());
        let i = input(1000., 0.2, bad, 15., 1.);
        assert!(m.evaluate(&b, &fresh(), 2., Some(&i)).is_err());
        let i = input(1000., 0.2, 200., 15., bad);
        assert!(m.evaluate(&b, &fresh(), 2., Some(&i)).is_err());
    }
    let mut i = input(1000., 0.2, 200., 15., 1.);
    i.climate.as_mut().unwrap().artifact_ref.clear();
    assert!(m.evaluate(&b, &fresh(), 2., Some(&i)).is_err());
    let i = input(1000., 0.2, 200., 15., 1.);
    let mut s = fresh();
    s.wet_coverage = 0.5;
    assert!(m.evaluate(&b, &s, 2., Some(&i)).is_err());
    let mut coated = fresh();
    coated.coating = Some(CoatingLayer {
        source_event: "test-old-joint".into(),
        mineral_mode: vec![ConstituentFraction {
            mineral: "goethite".into(),
            mass_fraction: 1.,
        }],
        grain_diameter_m: 1e-5,
        thickness_m: 1e-5,
        coverage: 1.,
        extinction_m_inv: 1e4,
    });
    assert!(m.evaluate(&b, &coated, 2., Some(&i)).is_err());
    let mut d: serde_json::Value =
        serde_json::from_str(include_str!("../data/surface_v1.json")).unwrap();
    d["parameters"]["hematite_extinction"]["value"] = json!(-1.);
    assert!(Model::from_json(&d.to_string()).is_err());
}
#[test]
fn spectral_data_are_valid_and_hold_the_existing_white_point() {
    let m = Model::embedded().unwrap();
    for v in m.spectrum_to_linear_rgb(&vec![1.; 81]).unwrap() {
        assert!((v - 1.).abs() < 1e-5);
    }
    assert!(m.spectrum_to_linear_rgb(&[0.; 80]).is_err());
    for (name, s) in m.spectra() {
        let vals: Vec<_> = s.reflectance.iter().map(|e| e.value).collect();
        assert_eq!(
            m.endmember_linear_rgb(name).unwrap(),
            m.spectrum_to_linear_rgb(&vals).unwrap()
        );
    }
}
#[test]
fn optional_evidence_export() {
    let Ok(out) = std::env::var("SURFACE_V1_EVIDENCE_DIR") else {
        eprintln!("SURFACE_V1_EVIDENCE_DIR unset; no artifact export requested");
        return;
    };
    let dir = std::path::Path::new(&out);
    std::fs::create_dir_all(dir).unwrap();
    let m = Model::embedded().unwrap();
    let b = body("sandstone_quartz.reference");
    let mut rows = Vec::new();
    for (name, i) in [
        ("fresh_fall", input(0., 0.2, 200., 15., 1.)),
        ("arid_1k_no_dust", input(1000., 0.2, 200., 15., 0.)),
        ("arid_10k_no_dust", input(10000., 0.2, 200., 15., 0.)),
        ("arid_10k_with_dust", input(10000., 0.2, 200., 15., 1.)),
        ("arid_100k_varnish", input(100000., 0.2, 200., 15., 1.)),
        ("humid_10k_rind", input(10000., 1.5, 1500., 18., 0.)),
        ("frozen_10k", input(10000., 0.2, 200., -10., 1.)),
        ("no_rain_10k", input(10000., 0.2, 0., 15., 1.)),
    ] {
        rows.push(json!({"name":name,"realization":b.reference(),"evaluation":m.evaluate(&b,&fresh(),2.,Some(&i)).unwrap()}));
    }
    let mut i = input(10000., 0.2, 200., 15., 1.);
    i.exposure = Exposure::Unknown {
        reason: "T1b collapse time not published".into(),
    };
    rows.push(json!({"name":"unknown_exposure","evaluation":m.evaluate(&b,&fresh(),2.,Some(&i)).unwrap()}));
    let colours: std::collections::BTreeMap<_, _> = m
        .spectra()
        .keys()
        .map(|k| (k, m.endmember_linear_rgb(k).unwrap()))
        .collect();
    std::fs::write(dir.join("samples.json"),serde_json::to_vec_pretty(&json!({"qualification":"manufactured climate/age/dust scenarios on actual catalog body; not generated Moab evidence","rows":rows,"endmember_linear_rgb":colours})).unwrap()).unwrap();
}

#[test]
fn fresh_surface_does_not_whiten_inherited_iron_rich_bulk() {
    let m = Model::embedded().unwrap();
    let b = body("gossan.reference");
    let (p, row) = predicted(&m, &b, &input(0., 0.2, 200., 15., 0.));
    let original = v3_2c::Evaluator::new(&b)
        .unwrap()
        .surface_sample(&fresh(), 2.)
        .unwrap();
    assert_eq!(
        row.base_reflectance_linear,
        original.base_reflectance_linear
    );
    assert_eq!(p.weights.fresh, 1.);
    assert!(p.recognized_oxide_fe_mass_fraction > 0.1);
    let (_, pale) = predicted(
        &m,
        &body("sandstone_quartz.reference"),
        &input(0., 0.2, 200., 15., 0.),
    );
    assert_ne!(row.base_reflectance_linear, pale.base_reflectance_linear);
}
