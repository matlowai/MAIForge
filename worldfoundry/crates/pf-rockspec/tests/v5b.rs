use pf_rockspec::{v3::*, v5b::*, *};
fn body(c: &v3::Catalogue, id: &str) -> ResolvedBody {
    let r = c.realization(id).unwrap();
    let reference = MaterialRef {
        spec_id: r.lithofacies.clone(),
        spec_version: "pf-rockspec/3".into(),
        realization_id: id.into(),
        body_id: "test-optics-body".into(),
    };
    let a = BodyAdmission {
        material_ref: reference.clone(),
        kind: BodyKind::Host,
        event_id: "fixture".into(),
        snapshot_id: "fixture".into(),
        geometry_ref: "fixture".into(),
        stock_account_ref: "fixture".into(),
        material_frame_ref: "fixture".into(),
        host_body_id: None,
    };
    c.resolve(
        reference,
        FormationContext {
            setting: r.settings[0].clone(),
            protolith: None,
            metamorphic_grade: None,
            burial_depth_m: None,
            provenance_ref: "test-only".into(),
        },
        a,
    )
    .unwrap()
}
fn fresh() -> SurfaceState {
    SurfaceState::fresh("same-face", "test-cut", [0.; 3], [0., 0., 1.])
}
fn y(s: &SurfaceSample) -> f64 {
    s.base_reflectance_linear
        .iter()
        .zip([0.2126, 0.7152, 0.0722])
        .map(|(x, w)| x * w)
        .sum()
}
fn energy(s: &SurfaceSample) {
    for j in 0..3 {
        assert!((0. ..=1.).contains(&s.base_reflectance_linear[j]));
        assert!(s.diffuse_energy_linear[j] + s.f0_linear[j] <= 1. + 1e-12);
    }
}

fn cut() -> SurfaceForm {
    SurfaceForm::Coherent {
        finish: SolidFinish::Cut,
        additional_rms_slope: 0.,
    }
}
fn powder(d: f64) -> SurfaceForm {
    SurfaceForm::Particulate {
        particle_diameter_m: d,
        additional_rms_slope: 0.,
    }
}
#[test]
fn form_does_not_mutate_stock_identity_or_legacy_outputs() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    let old = v5::Optics::embedded().unwrap();
    for r in &c.realizations {
        let b = body(&c, &r.id);
        let before = serde_json::to_vec(&b).unwrap();
        let v5 = old
            .prepare(&b, None)
            .unwrap()
            .surface_sample(&fresh(), 2.)
            .unwrap();
        for form in [
            cut(),
            SurfaceForm::Coherent {
                finish: SolidFinish::Broken,
                additional_rms_slope: 0.2,
            },
            SurfaceForm::Coherent {
                finish: SolidFinish::Outcrop,
                additional_rms_slope: 0.4,
            },
            powder(1e-4),
        ] {
            for seed in [None, Some(31)] {
                let p = o.prepare(&b, form, seed).unwrap();
                for f in [0.0001, 0.1, 50.] {
                    let sample = p.surface_sample(&fresh(), f).unwrap();
                    energy(&sample);
                    assert_eq!(sample.material_ref, b.reference().clone());
                    assert_eq!(sample.contract_version, VERSION);
                    assert_eq!(sample.dry_bulk_density_kg_m3, v5.dry_bulk_density_kg_m3);
                    assert_eq!(
                        sample.dry_deposit_density_kg_m3,
                        v5.dry_deposit_density_kg_m3
                    );
                    assert_eq!(
                        serde_json::to_vec(&sample).unwrap(),
                        serde_json::to_vec(&p.surface_sample(&fresh(), f).unwrap()).unwrap()
                    );
                }
            }
        }
        assert_eq!(before, serde_json::to_vec(&b).unwrap());
        assert_eq!(
            serde_json::to_vec(&v5).unwrap(),
            serde_json::to_vec(
                &old.prepare(&b, None)
                    .unwrap()
                    .surface_sample(&fresh(), 2.)
                    .unwrap()
            )
            .unwrap()
        );
    }
}
#[test]
fn form_changes_interfaces_not_absorption_or_rock_key_tints() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    for r in &c.realizations {
        let b = body(&c, &r.id);
        let solid = o.prepare(&b, cut(), Some(4)).unwrap();
        let grains = o.prepare(&b, powder(1e-6), Some(4)).unwrap();
        for (a, z) in solid.phases.iter().zip(&grains.phases) {
            assert_eq!(a.absorption_m_inv, z.absorption_m_inv);
            assert!(a.scattering_after_m_inv <= z.scattering_after_m_inv + 1e-9);
            assert!(a.chord_m > 0. && a.scattering_after_m_inv > 0.);
            if a.mineral.contains("glass") {
                assert_eq!(a.internal_scattering_retention, 0.);
            }
            for j in 0..3 {
                assert!(a.reflectance[j] <= z.reflectance[j] + 1e-12);
            }
        }
    }
}
#[test]
fn roughness_is_passive_escape_not_absorption_or_finish_name_magic() {
    let c = v3::Catalogue::embedded().unwrap();
    let b = body(&c, "granite.reference");
    let o = Optics::embedded().unwrap();
    let flat = o
        .prepare(&b, cut(), None)
        .unwrap()
        .surface_sample(&fresh(), 2.)
        .unwrap();
    let mut previous = y(&flat);
    for slope in [0., 0.2, 0.5, 1.] {
        let a = o
            .prepare(
                &b,
                SurfaceForm::Coherent {
                    finish: SolidFinish::Cut,
                    additional_rms_slope: slope,
                },
                None,
            )
            .unwrap()
            .surface_sample(&fresh(), 2.)
            .unwrap();
        let z = o
            .prepare(
                &b,
                SurfaceForm::Coherent {
                    finish: SolidFinish::Outcrop,
                    additional_rms_slope: slope,
                },
                None,
            )
            .unwrap()
            .surface_sample(&fresh(), 2.)
            .unwrap();
        assert_eq!(a.base_reflectance_linear, z.base_reflectance_linear);
        assert!(y(&a) <= previous + 1e-12);
        previous = y(&a);
        for j in 0..3 {
            let g = 1. / (1. + slope * slope);
            let r = flat.base_reflectance_linear[j];
            assert!((a.base_reflectance_linear[j] - r * g / (1. - r * (1. - g))).abs() < 1e-12);
        }
        assert!(a.roughness >= flat.roughness);
        energy(&a);
    }
}
#[test]
fn invalid_form_and_configuration_fail_without_unknown_to_powder_default() {
    for d in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(powder(d).validate().is_err());
    }
    for s in [-1., 1.1, f64::NAN] {
        assert!(SurfaceForm::Coherent {
            finish: SolidFinish::Broken,
            additional_rms_slope: s
        }
        .validate()
        .is_err());
    }
    let mut c = Config::embedded().unwrap();
    c.coherence_length_m.value = 0.;
    assert!(c.validate().is_err());
    assert!(serde_json::from_str::<SurfaceForm>(r#"{"form":"unknown"}"#).is_err());
    assert!(serde_json::from_str::<SurfaceForm>(
        r#"{"form":"particulate","additional_rms_slope":0}"#
    )
    .is_err());
}
#[test]
fn all_diagnostic_face_histories_are_preserved() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    let cases = DiagnosticSet::frozen_v2().unwrap();
    assert_eq!(cases.samples.len(), 128);
    for original in cases.samples {
        let (input, out) = o.diagnostic(&c, &original, cut()).unwrap();
        energy(&out);
        assert_eq!(
            serde_json::to_value(&input.surface).unwrap(),
            serde_json::to_value(&original.surface).unwrap()
        );
        if original.surface.surface_id.contains("altered_fresh") {
            assert_eq!(input.bulk, original.bulk);
        }
    }
}
#[test]
fn usgs_known_solid_regression_is_separate_from_validation() {
    let c = v3::Catalogue::embedded().unwrap();
    let b = body(&c, "basalt_tholeiitic.reference");
    let o = Optics::embedded().unwrap();
    let cal = v5::Calibration::embedded().unwrap();
    let refdata = cal.rocks["basalt_tholeiitic"]
        .observations
        .iter()
        .find(|r| r.record == 11848)
        .unwrap();
    let out = o
        .prepare(&b, cut(), None)
        .unwrap()
        .surface_sample(&fresh(), 50.)
        .unwrap();
    let rmse = (out
        .base_reflectance_linear
        .iter()
        .zip(refdata.rgb)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        / 3.)
        .sqrt();
    assert!(rmse < 0.015);
    println!("training target USGS11848 RGB RMSE={rmse}, not independent validation");
}
