use pf_rockspec::{v3, v3_1, v4, *};
use serde_json::{json, Value};

fn fixture(id: &str) -> v3_1::SurfaceFixture {
    let c = v3_1::catalogue().unwrap();
    v3_1::FixtureSet::comparisons(&c)
        .unwrap()
        .samples
        .into_iter()
        .find(|s| s.realization_id == id && s.history == "fresh_dry" && s.footprint_m == 2.)
        .unwrap()
}
fn optical_equal(a: SurfaceSample, b: SurfaceSample) {
    // Versioned identity changes; every optical/density/history output must match.
    let normalize = |s: SurfaceSample| {
        let mut v = serde_json::to_value(s).unwrap();
        v["contract_version"] = json!("comparison");
        v["numerical_model"] = json!("comparison");
        v["material_ref"]["spec_version"] = json!("comparison");
        v
    };
    assert!(
        normalize(a) == normalize(b),
        "physical output differs on unchanged inputs"
    );
}
fn energy(s: &SurfaceSample) {
    for j in 0..3 {
        assert!((0. ..=1.).contains(&s.base_reflectance_linear[j]));
        assert!(s.diffuse_energy_linear[j] + s.f0_linear[j] <= 1. + 1e-12);
    }
}
fn y(s: &SurfaceSample) -> f64 {
    s.base_reflectance_linear
        .iter()
        .zip([0.2126, 0.7152, 0.0722])
        .map(|(a, b)| a * b)
        .sum()
}

#[test]
fn opt_in_versions_reject_rejected_v31_and_wrong_identity() {
    let c = v3_1::catalogue().unwrap();
    assert_eq!(c.core.schema_version, "pf-rockspec/3.1b");
    assert_eq!(c.core.numerical_model, "assemblage+connected-km/3.1b");
    assert_eq!(
        Catalogue::embedded().unwrap().schema_version,
        "pf-rockspec/1"
    );
    assert_eq!(
        Catalogue::embedded_v2().unwrap().schema_version,
        "pf-rockspec/2"
    );
    assert_eq!(
        v3::Catalogue::embedded().unwrap().core.schema_version,
        "pf-rockspec/3"
    );
    let f = fixture("granite.i_type_pink");
    let b = f.resolve_fixture(&c).unwrap();
    assert!(c
        .core
        .realize_prior("granite", b.reference().clone())
        .is_err());
    assert!(v4::Optics::embedded()
        .unwrap()
        .surface_sample(&b, &f.surface, 2.)
        .is_err());
    let mut bad = serde_json::to_value(&c).unwrap();
    bad["schema_version"] = json!("pf-rockspec-assemblage/3.1");
    bad["core"]["schema_version"] = json!("pf-rockspec/3.1");
    bad["core"]["numerical_model"] = json!("assemblage+grain-transport/3.1");
    assert!(v3::Catalogue::from_json(&bad.to_string()).is_err());
    let mut bad = serde_json::to_value(&b).unwrap();
    bad["reference"]["spec_version"] = json!("pf-rockspec/3");
    assert!(v3::ResolvedBody::from_json(&bad.to_string()).is_err());
}

#[test]
fn unchanged_65_realizations_including_all_granites_match_v3_exactly() {
    let c = v3_1::catalogue().unwrap();
    let old = v3::Catalogue::embedded().unwrap();
    let cases = v3_1::FixtureSet::comparisons(&c).unwrap();
    let mut n = 0;
    let mut granite = 0;
    for r in &c.realizations {
        if ["basalt_tholeiitic.lip", "basalt_tholeiitic.island_shield"].contains(&r.id.as_str()) {
            continue;
        }
        let f = cases
            .samples
            .iter()
            .find(|s| s.realization_id == r.id)
            .unwrap();
        let a = f.resolve_fixture(&old).unwrap();
        let b = f.resolve_fixture(&c).unwrap();
        for f in cases.samples.iter().filter(|s| s.realization_id == r.id) {
            optical_equal(
                a.surface_sample(&f.surface, f.footprint_m).unwrap(),
                b.surface_sample(&f.surface, f.footprint_m).unwrap(),
            );
            n += 1;
            if r.lithofacies == "granite" {
                granite += 1;
            }
        }
    }
    assert_eq!(n, 585);
    assert_eq!(granite, 54);
}

#[test]
fn retained_basalt_modes_use_v3_optics_and_keep_requested_dark_range() {
    let c = v3_1::catalogue().unwrap();
    let old = v3::Catalogue::embedded().unwrap();
    let mut count = 0;
    for r in c
        .realizations
        .iter()
        .filter(|r| r.lithofacies.starts_with("basalt_"))
    {
        let original = old.realization(&r.id).unwrap();
        assert_eq!(
            serde_json::to_value(&r.parameters).unwrap(),
            serde_json::to_value(&original.parameters).unwrap()
        );
        for key in [
            "grain",
            "optical_domain_diameter_m",
            "structure_v2",
            "intrinsic_porosity",
            "vesicularity",
        ] {
            assert_eq!(
                serde_json::to_value(&r.spec).unwrap()[key],
                serde_json::to_value(&original.spec).unwrap()[key]
            );
        }
        let f = fixture(&r.id);
        let b = f.resolve_fixture(&c).unwrap();
        // Independent old-version catalogue evaluates the same modal change.
        let mut v3_mode = old.clone();
        v3_mode
            .core
            .references
            .insert("RS31".into(), c.core.references["RS31"].clone());
        v3_mode
            .realizations
            .iter_mut()
            .find(|x| x.id == r.id)
            .unwrap()
            .spec = r.spec.clone();
        v3_mode.validate().unwrap();
        let expected = f.resolve_fixture(&v3_mode).unwrap();
        for scale in [0.1, 2., 50.] {
            let out = b.surface_sample(&f.surface, scale).unwrap();
            assert!(
                (0.05..=0.10).contains(&y(&out)),
                "{} {scale}: Y={}",
                r.id,
                y(&out)
            );
            optical_equal(out, expected.surface_sample(&f.surface, scale).unwrap());
            count += 1;
        }
    }
    assert_eq!(count, 15);
    // Same mode under a different setting has no hidden tint.
    let mut same = c.clone();
    let reference = same
        .realization("basalt_tholeiitic.reference")
        .unwrap()
        .clone();
    let lip = same
        .realizations
        .iter_mut()
        .find(|r| r.id == "basalt_tholeiitic.lip")
        .unwrap();
    lip.spec = reference.spec;
    lip.parameters = reference.parameters;
    let f = fixture("basalt_tholeiitic.lip");
    let a = f.resolve_fixture(&same).unwrap();
    let b = fixture("basalt_tholeiitic.reference")
        .resolve_fixture(&same)
        .unwrap();
    assert_eq!(
        a.surface_sample(&f.surface, 2.)
            .unwrap()
            .base_reflectance_linear,
        b.surface_sample(&f.surface, 2.)
            .unwrap()
            .base_reflectance_linear
    );
}

#[test]
fn exactly_54_frozen_surfaces_use_v3_optics_without_bulk_reaction() {
    let c = v3_1::catalogue().unwrap();
    let old = v3::Catalogue::embedded().unwrap();
    let set = v3_1::FixtureSet::diagnostics().unwrap();
    let swatches = v3_1::FixtureSet::swatches().unwrap();
    assert_eq!(set.samples.len(), 54);
    assert_eq!(
        serde_json::to_value(&set).unwrap(),
        serde_json::to_value(&swatches).unwrap()
    );
    let mut weather = 0;
    for f in set.samples {
        let b = f.resolve_fixture(&c).unwrap();
        let a = f.resolve_fixture(&old).unwrap();
        let before = serde_json::to_value(&b).unwrap();
        let sample = b.surface_sample(&f.surface, f.footprint_m).unwrap();
        energy(&sample);
        optical_equal(sample, a.surface_sample(&f.surface, f.footprint_m).unwrap());
        assert_eq!(f.weathering_family, b.spec().weathering_family);
        if f.history == "weathered_dry" {
            weather += 1;
            assert!(f.surface.exposure_elapsed_s.is_none());
            assert!(matches!(f.surface.mechanism, FaceMechanism::Inherited));
            let l = f.surface.coating.as_ref().unwrap();
            assert!(
                (l.mineral_mode.iter().map(|x| x.mass_fraction).sum::<f64>() - 1.).abs() < 1e-12
            );
            for m in &l.mineral_mode {
                assert!(b.spec().weathering_products.contains(&m.mineral));
            }
            if b.reference().spec_id == "banded_iron_formation" {
                assert_eq!(f.weathering_family, "oxidation_residual");
                assert!(l
                    .mineral_mode
                    .iter()
                    .all(|x| ["quartz", "hematite", "goethite"].contains(&x.mineral.as_str())));
            }
            let fresh = fixture(&f.realization_id).surface;
            assert!(b
                .surface_sample(&fresh, f.footprint_m)
                .unwrap()
                .coating
                .is_none());
        } else {
            assert!(f.surface.coating.is_none());
        }
        assert_eq!(before, serde_json::to_value(&b).unwrap());
    }
    assert_eq!(weather, 18);
}

#[test]
fn all_603_comparison_cases_are_passive_deterministic_and_pure() {
    let c = v3_1::catalogue().unwrap();
    let set = v3_1::FixtureSet::comparisons(&c).unwrap();
    assert_eq!(set.samples.len(), 603);
    for r in &c.realizations {
        let cases: Vec<_> = set
            .samples
            .iter()
            .filter(|s| s.realization_id == r.id)
            .collect();
        let b = cases[0].resolve_fixture(&c).unwrap();
        let before = serde_json::to_value(&b).unwrap();
        for f in cases {
            let a = b.surface_sample(&f.surface, f.footprint_m).unwrap();
            energy(&a);
            assert_eq!(a, b.surface_sample(&f.surface, f.footprint_m).unwrap());
            assert_eq!(a.contract_version, v3_1::VERSION);
            assert_eq!(a.numerical_model, v3_1::MODEL);
            assert_eq!(a.material_ref, *b.reference());
        }
        assert_eq!(before, serde_json::to_value(&b).unwrap());
    }
}

#[test]
fn fixture_priors_match_actual_parameters_and_modes() {
    let c = v3_1::catalogue().unwrap();
    let w: Value =
        serde_json::from_str(include_str!("../data/weathering_priors_v3_1.json")).unwrap();
    let p: Value =
        serde_json::from_str(include_str!("../data/realization_priors_v3_1.json")).unwrap();
    assert_eq!(w["families"].as_object().unwrap().len(), 6);
    assert_eq!(p["realizations"].as_object().unwrap().len(), 2);
    for (id, q) in p["realizations"].as_object().unwrap() {
        assert_eq!(
            q["mode"],
            serde_json::to_value(&c.realization(id).unwrap().spec.mode).unwrap()
        );
    }
    for f in v3_1::FixtureSet::diagnostics().unwrap().samples {
        if f.history != "weathered_dry" {
            continue;
        }
        let q = &w["families"][f.realization_id.split('.').next().unwrap()];
        let l = f.surface.coating.unwrap();
        assert_eq!(q["family"], f.weathering_family);
        for (name, v) in [
            ("coverage", l.coverage),
            ("thickness_m", l.thickness_m),
            ("grain_diameter_m", l.grain_diameter_m),
            ("rms_slope_history_rad", f.surface.rms_slope_history_rad),
            ("abrasion", f.surface.abrasion),
        ] {
            assert_eq!(q[name]["value"].as_f64().unwrap(), v);
            assert_eq!(q[name]["reference"], "RS31");
            assert_eq!(q[name]["evidence_kind"], "project_prior");
            assert!(
                q[name]["interval"][0].as_f64().unwrap() <= v
                    && v <= q[name]["interval"][1].as_f64().unwrap()
            );
        }
        for m in l.mineral_mode {
            let r = q["residue_mode"]
                .as_array()
                .unwrap()
                .iter()
                .find(|x| x["mineral"] == m.mineral)
                .unwrap();
            assert_eq!(
                r["mass_fraction"]["value"].as_f64().unwrap(),
                m.mass_fraction
            );
        }
    }
}

#[test]
fn malformed_versions_and_duplicate_fixture_ids_fail() {
    let mut c = v3_1::catalogue().unwrap();
    c.schema_version = "pf-rockspec-assemblage/3".into();
    assert!(c.validate().is_err());
    let mut s = v3_1::FixtureSet::swatches().unwrap();
    s.samples.push(s.samples[0].clone());
    assert!(v3_1::FixtureSet::from_json(&serde_json::to_string(&s).unwrap()).is_err());
}
