use pf_rockspec::{v3::*, v4::*, *};
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
fn transition_is_monotone_has_exact_limits_and_is_not_a_rock_key() {
    for length in [1e-6, 1e-4, 1e-2] {
        let mut previous = 0.;
        for i in 1..=1000 {
            let grain = i as f64 * length / 100.;
            let w = transition(grain, length, [1., 5.]).unwrap();
            assert!(w >= previous);
            assert!((w - transition(grain * 7., length * 7., [1., 5.]).unwrap()).abs() < 1e-14);
            previous = w;
        }
        assert_eq!(transition(length, length, [1., 5.]).unwrap(), 0.);
        assert_eq!(transition(5. * length, length, [1., 5.]).unwrap(), 1.);
    }
    for bad in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(transition(bad, 1., [1., 5.]).is_err());
    }
}
#[test]
fn every_mineral_is_passive_and_single_grain_darkens_with_size() {
    let c = v3::Catalogue::embedded().unwrap();
    for mineral in &c.core.minerals {
        let modes = vec![ConstituentFraction {
            mineral: mineral.id.clone(),
            mass_fraction: 1.,
        }];
        let mut previous = [1.; 3];
        for d in [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1] {
            let r = mixing(&c.core, &modes, d, d, [1., 5.]).unwrap();
            for j in 0..3 {
                assert!((r.areal[j] - r.intimate[j]).abs() < 1e-12);
                assert!(r.reflectance[j] <= previous[j] + 1e-12);
                assert!((0. ..=1.).contains(&r.reflectance[j]));
            }
            previous = r.reflectance;
        }
    }
}
#[test]
fn heterogeneous_transition_brightens_when_coefficients_are_fixed() {
    let c = v3::Catalogue::embedded().unwrap();
    let modes = vec![
        ConstituentFraction {
            mineral: "quartz".into(),
            mass_fraction: 0.7,
        },
        ConstituentFraction {
            mineral: "magnetite".into(),
            mass_fraction: 0.3,
        },
    ];
    let mut prev = [0.; 3];
    for i in 1..=200 {
        // Hold K,S fixed to isolate the closure. Changing grain also changes S,
        // so total mixed-rock reflectance has no universal monotonicity theorem.
        let r = mixing(&c.core, &modes, i as f64 * 1e-5, 0.001, [1., 5.]).unwrap();
        for j in 0..3 {
            assert!(r.areal[j] >= r.intimate[j]);
            assert!(r.reflectance[j] >= prev[j] - 1e-12);
            assert!(
                r.reflectance[j] >= r.intimate[j] - 1e-12 && r.reflectance[j] <= r.areal[j] + 1e-12
            );
        }
        prev = r.reflectance;
    }
}
#[test]
fn slab_passivity_limits_and_bright_substrate_darkening() {
    let mut count = 0;
    for k in [0., 1e-3, 1e3, 1e6, 1e9] {
        for s in [0., 1e-3, 1e3, 1e6, 1e9] {
            for t in [0., 1e-9, 1e-6, 1e-3, 1.] {
                let (r, tr) = slab(k, s, t).unwrap();
                assert!(r + tr <= 1. + 1e-12);
                for sub in [0., 0.1, 0.5, 1.] {
                    assert!(layer_over(sub, r, tr).unwrap() <= 1. + 1e-12);
                    count += 1;
                }
            }
        }
    }
    assert_eq!(count, 500);
    assert_eq!(slab(1., 1., 0.).unwrap(), (0., 1.));
    let (r, t) = slab(0., 1e6, 1e-3).unwrap();
    assert!((r + t - 1.).abs() < 1e-12);
    assert!((layer_over(1., r, t).unwrap() - 1.).abs() < 1e-12);
    let (r, t) = slab(2e6, 1e6, 1.).unwrap();
    assert!((r - km_reflectance(2.)).abs() < 1e-12 && t == 0.);
    let mut previous = 0.8;
    for i in 0..=100 {
        let (r, t) = slab(2e6, 1e6, i as f64 * 1e-7).unwrap();
        let current = layer_over(0.8, r, t).unwrap();
        assert!(current <= previous + 1e-12);
        previous = current;
    }
    assert!(slab(-1., 1., 1.).is_err());
    assert!(layer_over(0.5, 0.8, 0.8).is_err());
}
#[test]
fn morphology_closes_stock_and_geometry_and_coats_darken_all_six() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    assert_eq!(o.sediments.len(), 6);
    for id in o.sediments.keys() {
        let b = body(&c, id);
        let before = serde_json::to_vec(&b).unwrap();
        let r = o.morphology_receipt(&b).unwrap().unwrap();
        assert!(
            (r.core_mass_fraction + r.coating_mass_fraction + r.cement_mass_fraction - 1.).abs()
                < 1e-12
        );
        assert!(
            (r.reconstructed_coating_volume_m3_per_kg / r.coating_volume_m3_per_kg - 1.).abs()
                < 1e-12
        );
        assert!(r.thickness_within_prior_support);
        assert!(
            y(&o.surface_sample(&b, &fresh(), 2.).unwrap())
                < y(&o.without_coatings(&b, &fresh(), 2.).unwrap())
        );
        assert_eq!(before, serde_json::to_vec(&b).unwrap());
        if id.starts_with("sand.") || id.starts_with("regolith.") {
            assert_eq!(r.cement_mass_fraction, 0.);
        }
    }
}
#[test]
fn malformed_prior_and_unfunded_coating_are_rejected() {
    let mut o = Optics::embedded().unwrap();
    o.sediments
        .get_mut("sand.reference")
        .unwrap()
        .coverage
        .value = 1.5;
    assert!(Optics::from_json(&serde_json::to_string(&o).unwrap()).is_err());
    let mut o = Optics::embedded().unwrap();
    o.sediments.get_mut("sand.reference").unwrap().coating_mode[0].mass_fraction = 0.99;
    let c = v3::Catalogue::embedded().unwrap();
    let b = body(&c, "sand.reference");
    assert!(o.surface_sample(&b, &fresh(), 2.).is_err());
    let mut o = Optics::embedded().unwrap();
    o.sediments
        .get_mut("sand.reference")
        .unwrap()
        .pigment_scattering_length_m
        .unit = "nm".into();
    assert!(o.validate().is_err());
}
#[test]
fn all_realizations_three_scales_dry_wet_are_passive_and_pure() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    let mut n = 0;
    for r in &c.realizations {
        let b = body(&c, &r.id);
        let before = serde_json::to_vec(&b).unwrap();
        for footprint in [0.1, 2., 50.] {
            for wet in [0., 1.] {
                let mut s = fresh();
                s.wet_coverage = wet;
                s.water_film_m = wet * 1e-4;
                let a = o.surface_sample(&b, &s, footprint).unwrap();
                energy(&a);
                assert_eq!(a, o.surface_sample(&b, &s, footprint).unwrap());
                assert_eq!(a.material_ref, b.reference().clone());
                assert_eq!(a.contract_version, "pf-rockspec-optics/4");
                n += 1;
            }
        }
        assert_eq!(before, serde_json::to_vec(&b).unwrap());
    }
    assert_eq!(n, 402);
}
#[test]
fn frozen_diagnostics_migrate_all_128_without_losing_surface_history() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    let set = DiagnosticSet::frozen_v2().unwrap();
    assert_eq!(set.samples.len(), 128);
    for old in set.samples {
        let (input, out) = o.diagnostic(&c, &old).unwrap();
        energy(&out);
        assert_eq!(
            serde_json::to_vec(&input.surface).unwrap(),
            serde_json::to_vec(&old.surface).unwrap()
        );
        assert_eq!(input.footprint_m, old.footprint_m);
        if old.surface.surface_id.contains("altered_fresh") {
            assert_eq!(input.bulk.mineral_mode, old.bulk.mineral_mode);
        }
    }
}
#[test]
fn protected_fine_and_single_phase_cases_keep_their_physical_limits() {
    let c = v3::Catalogue::embedded().unwrap();
    let o = Optics::embedded().unwrap();
    for id in ["chalk.reference", "shale.organic_rich"] {
        let b = body(&c, id);
        let old = b.surface_sample(&fresh(), 2.).unwrap();
        let new = o.surface_sample(&b, &fresh(), 2.).unwrap();
        assert!((y(&new) - y(&old)).abs() < 1e-12);
        assert_eq!(o.mixing_receipt(&b).unwrap().areal_weight, [0.; 3]);
    }
    for id in ["marble.reference", "obsidian.reference"] {
        let b = body(&c, id);
        let old = b.surface_sample(&fresh(), 2.).unwrap();
        let new = o.surface_sample(&b, &fresh(), 2.).unwrap();
        assert!((y(&new) - y(&old)).abs() < 0.001);
    }
    // Basalt has a small but nonzero crossover. Do not disguise it as key-off equality.
    let b = body(&c, "basalt_tholeiitic.reference");
    let w = o.mixing_receipt(&b).unwrap().areal_weight;
    assert!(w.iter().all(|x| *x > 0. && *x < 0.1));
}
