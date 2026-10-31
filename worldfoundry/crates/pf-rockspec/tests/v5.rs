use pf_rockspec::{v3::*, v5::*, *};
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
#[test]
fn calibrated_bodies_are_passive_repeatable_and_leave_admitted_state_untouched() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    let mut changed = 0;
    for r in &c.realizations {
        let b = body(&c, &r.id);
        let before = serde_json::to_vec(&b).unwrap();
        let old = b.surface_sample(&fresh(), 2.).unwrap();
        for seed in [None, Some(0), Some(1), Some(u64::MAX)] {
            let p = o.prepare(&b, seed).unwrap();
            for f in [0.00001, 0.1, 2., 50.] {
                let a = p.surface_sample(&fresh(), f).unwrap();
                energy(&a);
                assert_eq!(a.material_ref, b.reference().clone());
                assert_eq!(a.contract_version, VERSION);
                assert_eq!(a.dry_bulk_density_kg_m3, old.dry_bulk_density_kg_m3);
                assert_eq!(
                    serde_json::to_vec(&a).unwrap(),
                    serde_json::to_vec(&p.surface_sample(&fresh(), f).unwrap()).unwrap()
                );
            }
        }
        if (y(&o
            .prepare(&b, None)
            .unwrap()
            .surface_sample(&fresh(), 2.)
            .unwrap())
            - y(&old))
        .abs()
            > 1e-6
        {
            changed += 1;
        }
        assert_eq!(before, serde_json::to_vec(&b).unwrap());
        assert_eq!(
            serde_json::to_vec(&old).unwrap(),
            serde_json::to_vec(&b.surface_sample(&fresh(), 2.).unwrap()).unwrap()
        );
    }
    assert!(changed > 50);
    println!("67 bodies x4 seeds x4 footprints passive; {changed} opt-in differences; old body bytes exact");
}
#[test]
fn migration_preserves_all_128_face_histories_and_altered_modes() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    let cases = DiagnosticSet::frozen_v2().unwrap();
    assert_eq!(cases.samples.len(), 128);
    for original in cases.samples {
        let (input, output) = o.diagnostic(&c, &original).unwrap();
        energy(&output);
        assert_eq!(
            serde_json::to_value(&input.surface).unwrap(),
            serde_json::to_value(&original.surface).unwrap()
        );
        assert_eq!(input.footprint_m, original.footprint_m);
        if original.surface.surface_id.contains("altered_fresh") {
            assert_eq!(input.bulk, original.bulk);
        }
    }
}
#[test]
fn malformed_calibration_is_rejected_and_unmeasured_minerals_keep_priors() {
    let c = Calibration::embedded().unwrap();
    assert_eq!(c.minerals.len(), 56);
    assert_eq!(
        c.minerals
            .values()
            .filter(|p| !p.observations.is_empty())
            .count(),
        41
    );
    for k in &c.missing_minerals {
        assert!(c.minerals[k].observations.is_empty());
    }
    let mut bad = c.clone();
    bad.minerals.get_mut("quartz").unwrap().observations[0].rgb[1] = -0.1;
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.transport_transition_lengths = [5., 1.];
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.minerals.remove("quartz");
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.oxide_scattering_multiplier = 0.;
    assert!(bad.validate().is_err());
    let cat = v3::Catalogue::embedded().unwrap();
    let b = body(&cat, "obsidian.reference");
    let o = Optics::new(c).unwrap();
    let p = o.prepare(&b, None).unwrap();
    let glass = p
        .receipts
        .iter()
        .find(|r| r.mineral == "volcanic_glass")
        .unwrap();
    assert!(glass.records.is_empty());
    assert_eq!(
        glass.mean,
        cat.core
            .minerals
            .iter()
            .find(|m| m.id == "volcanic_glass")
            .unwrap()
            .diffuse_reflectance
            .clone()
            .map(|q| q.value)
    );
}
#[test]
fn oxide_cross_section_changes_mixing_without_changing_pure_grain_reflectance() {
    let cat = v3::Catalogue::embedded().unwrap();
    let b = body(&cat, "basalt_tholeiitic.reference");
    let mut c = Calibration::embedded().unwrap();
    c.oxide_scattering_multiplier = 1.;
    let o1 = Optics::new(c.clone()).unwrap();
    c.oxide_scattering_multiplier = 32.;
    let o2 = Optics::new(c).unwrap();
    for mineral in ["hematite", "goethite"] {
        let mut bulk = b.bulk().state.clone();
        bulk.mineral_mode = vec![ConstituentFraction {
            mineral: mineral.into(),
            mass_fraction: 1.,
        }];
        for d in [1e-6, 1e-4, 1e-2] {
            let m = Material {
                catalogue: &cat.core,
                spec: b.spec(),
                reference: b.reference().clone(),
                bulk: bulk.clone(),
                grain_diameter_m: d,
            };
            let a = o1
                .prepare_material(&m, b.parameters(), None)
                .unwrap()
                .mixing_receipt()
                .unwrap();
            let z = o2
                .prepare_material(&m, b.parameters(), None)
                .unwrap()
                .mixing_receipt()
                .unwrap();
            for j in 0..3 {
                assert!((a.reflectance[j] - z.reflectance[j]).abs() < 1e-12);
            }
        }
    }
}
#[test]
fn empirical_draws_are_correlated_per_body_and_not_per_face() {
    let c = v3::Catalogue::embedded().unwrap();
    let b = body(&c, "granite.i_type_pink");
    let o = Optics::embedded().unwrap();
    let a = o.prepare(&b, Some(17)).unwrap();
    let repeat = o.prepare(&b, Some(17)).unwrap();
    let other = o.prepare(&b, Some(18)).unwrap();
    let x = a.surface_sample(&fresh(), 50.).unwrap();
    let z = repeat.surface_sample(&fresh(), 50.).unwrap();
    assert_eq!(x.base_reflectance_linear, z.base_reflectance_linear);
    assert_ne!(
        x.base_reflectance_linear,
        other
            .surface_sample(&fresh(), 50.)
            .unwrap()
            .base_reflectance_linear
    );
    assert_eq!(
        x.material_ref,
        other.surface_sample(&fresh(), 50.).unwrap().material_ref
    );
}
