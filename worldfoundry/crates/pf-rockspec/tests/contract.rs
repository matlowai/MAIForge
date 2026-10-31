use pf_rockspec::*;
fn material<'a>(c: &'a Catalogue, key: &str) -> Material<'a> {
    c.realize_prior(
        key,
        MaterialRef {
            spec_id: key.into(),
            spec_version: c.schema_version.clone(),
            realization_id: "test-prior-v1".into(),
            body_id: "test-body".into(),
        },
    )
    .unwrap()
}
fn fresh() -> SurfaceState {
    SurfaceState::fresh("test-face", "test-cut", [0.0; 3], [0.0, 0.0, 1.0])
}
#[test]
fn complete_keys_and_geological_distinctions() {
    let c = Catalogue::embedded().unwrap();
    let expected="basalt_tholeiitic basalt_alkali andesite dacite rhyolite obsidian scoria pumice ignimbrite tuff_ash granite granodiorite diorite gabbro peridotite sandstone_quartz arkose greywacke conglomerate siltstone mudstone shale limestone dolostone chalk chert rock_salt gypsum_rock coal slate phyllite schist gneiss marble quartzite serpentinite amphibolite hornfels gravel sand silt_clay till regolith";
    assert_eq!(c.rocks.len(), 43);
    for key in expected.split_whitespace() {
        c.rock(key).unwrap();
    }
    assert_eq!(c.minerals.len(), 32);
    assert_eq!(c.mineral("quartz").unwrap().cleavage.directions.value, 0.0);
    assert_eq!(c.mineral("calcite").unwrap().cleavage.directions.value, 3.0);
    assert_eq!(
        c.mineral("muscovite").unwrap().cleavage.directions.value,
        1.0
    );
    assert_eq!(
        c.mineral("hornblende").unwrap().cleavage.angles_deg[0].value,
        56.0
    );
    assert!(matches!(
        c.rock("sand").unwrap().mechanical["ucs_pa"],
        Property::Unknown { .. }
    ));
    assert!(c.rock("sand").unwrap().joint_sets.is_empty());
    assert!(c
        .rock("basalt_tholeiitic")
        .unwrap()
        .joint_sets
        .iter()
        .any(|j| j.relation == "cooling_column_side"));
}
#[test]
fn frozen_all_finite_closed_energy_and_exact_determinism() {
    let c = Catalogue::embedded().unwrap();
    let cases = DiagnosticSet::frozen().unwrap();
    assert_eq!(cases.samples.len(), 128);
    let first: Vec<_> = cases
        .samples
        .iter()
        .map(|case| case.evaluate(&c).unwrap())
        .collect();
    for (i, case) in cases.samples.iter().enumerate().rev() {
        let v = case.evaluate(&c).unwrap();
        assert_eq!(v, first[i]);
        assert_eq!(
            serde_json::to_vec(&v).unwrap(),
            serde_json::to_vec(&first[i]).unwrap()
        );
        for j in 0..3 {
            assert!(v.diffuse_energy_linear[j] + v.f0_linear[j] <= 1.0 + 1e-12);
        }
        assert!((0.0..=1.0).contains(&v.roughness));
    }
}
#[test]
fn all_387_swatch_states_evaluate() {
    let c = Catalogue::embedded().unwrap();
    let set = DiagnosticSet::swatches().unwrap();
    assert_eq!(set.samples.len(), 387);
    for s in set.samples {
        s.evaluate(&c).unwrap();
    }
}
#[test]
fn fresnel_and_transport_analytic_limits() {
    assert!((fresnel_f0(1.5, 0., 1.) - 0.04).abs() < 1e-15);
    assert!((fresnel_f0(1.5, 2., 1.) - 4.25 / 10.25).abs() < 1e-15);
    assert_eq!(km_reflectance(0.0), 1.0);
    assert_eq!(wet_diffuse(0.0, 1.333), 0.0);
    assert_eq!(wet_diffuse(1.0, 1.333), 1.0);
    for r in [0.01, 0.1, 0.3, 0.8, 0.99] {
        assert!(wet_diffuse(r, 1.333) < r);
        assert!((wet_diffuse(r, 1.) - r).abs() < 1e-14);
    }
    assert!(km_reflectance(100.) < km_reflectance(1.));
}
#[test]
fn wet_darker_smoother_without_bulk_mutation() {
    let c = Catalogue::embedded().unwrap();
    for r in &c.rocks {
        let m = material(&c, &r.id);
        let before = serde_json::to_vec(&m.bulk).unwrap();
        let mut s = fresh();
        let dry = surface_sample(&m, &s, 0.1).unwrap();
        s.wet_coverage = 1.;
        s.water_film_m = 0.0001;
        let wet = surface_sample(&m, &s, 0.1).unwrap();
        assert!(wet.roughness < dry.roughness);
        for j in 0..3 {
            assert!(wet.base_reflectance_linear[j] < dry.base_reflectance_linear[j]);
        }
        assert_eq!(serde_json::to_vec(&m.bulk).unwrap(), before);
    }
}
#[test]
fn footprint_filters_detail_not_material_or_mean_albedo() {
    let c = Catalogue::embedded().unwrap();
    let m = material(&c, "granite");
    let s = fresh();
    let a = surface_sample(&m, &s, 1e-5).unwrap();
    let b = surface_sample(&m, &s, 50.).unwrap();
    assert_eq!(a.base_reflectance_linear, b.base_reflectance_linear);
    assert_eq!(a.material_ref, b.material_ref);
    assert_eq!(a.normal_detail.wavelength_m, b.normal_detail.wavelength_m);
    assert!(a.normal_detail.resolved_rms_slope_rad > b.normal_detail.resolved_rms_slope_rad);
    assert!(a.roughness < b.roughness);
    assert!(a.sparkle.unwrap().resolved_coverage > b.sparkle.unwrap().resolved_coverage);
    for x in [a.normal_detail, b.normal_detail] {
        assert!(
            (x.resolved_rms_slope_rad.powi(2) + x.unresolved_slope_variance
                - x.unfiltered_rms_slope_rad.powi(2))
            .abs()
                < 1e-12
        );
    }
}
#[test]
fn grain_size_changes_transport_without_changing_mass_mode() {
    let c = Catalogue::embedded().unwrap();
    let mut m = material(&c, "granite");
    let s = fresh();
    m.grain_diameter_m = 0.00001;
    let fine = surface_sample(&m, &s, 0.1).unwrap();
    m.grain_diameter_m = 0.01;
    let coarse = surface_sample(&m, &s, 0.1).unwrap();
    for j in 0..3 {
        assert!(fine.base_reflectance_linear[j] > coarse.base_reflectance_linear[j]);
    }
    assert_eq!(fine.dry_bulk_density_kg_m3, coarse.dry_bulk_density_kg_m3);
}
#[test]
fn fresh_cut_retains_altered_bulk_and_old_face_retains_coating() {
    let c = Catalogue::embedded().unwrap();
    let mut m = material(&c, "basalt_tholeiitic");
    let f = fresh();
    let original = surface_sample(&m, &f, 0.1).unwrap();
    m.bulk.mineral_mode = vec![ConstituentFraction {
        mineral: "chlorite".into(),
        mass_fraction: 1.,
    }];
    m.bulk.alteration_record = "buried replacement admitted by test fixture".into();
    let altered = surface_sample(&m, &f, 0.1).unwrap();
    assert_ne!(
        altered.base_reflectance_linear,
        original.base_reflectance_linear
    );
    let mut old = f.clone();
    old.coating = Some(CoatingLayer {
        source_event: "old-rind".into(),
        mineral_mode: vec![ConstituentFraction {
            mineral: "goethite".into(),
            mass_fraction: 1.,
        }],
        grain_diameter_m: 0.00001,
        thickness_m: 0.001,
        coverage: 1.,
        extinction_m_inv: 10000.,
    });
    let old_result = surface_sample(&m, &old, 0.1).unwrap();
    assert_ne!(
        old_result.base_reflectance_linear,
        altered.base_reflectance_linear
    );
    assert_eq!(old_result.substrate, altered.substrate);
    assert!(old_result.coating.unwrap().optical_weight > 0.99999);
    old.coating.as_mut().unwrap().thickness_m = 0.;
    let zero = surface_sample(&m, &old, 0.1).unwrap();
    assert_eq!(
        zero.base_reflectance_linear,
        altered.base_reflectance_linear
    );
    assert_eq!(zero.f0_linear, altered.f0_linear);
    assert_eq!(zero.roughness, altered.roughness);
}
#[test]
fn exposure_time_alone_never_selects_a_colour() {
    let c = Catalogue::embedded().unwrap();
    let m = material(&c, "granite");
    let mut s = fresh();
    let a = surface_sample(&m, &s, 2.).unwrap();
    s.exposure_elapsed_s = Some(1e12);
    let b = surface_sample(&m, &s, 2.).unwrap();
    assert_eq!(a, b);
    s.exposure_elapsed_s = None;
    assert_eq!(a, surface_sample(&m, &s, 2.).unwrap());
}
#[test]
fn fabrication_frame_controls_anisotropy_and_world_motion_does_not_swim() {
    let c = Catalogue::embedded().unwrap();
    let mut m = material(&c, "schist");
    let mut s = fresh();
    s.mechanism = FaceMechanism::Cleavage;
    let a = surface_sample(&m, &s, 0.0001).unwrap();
    s.normal_material = [1., 0., 0.];
    let b = surface_sample(&m, &s, 0.0001).unwrap();
    assert!(a.normal_detail.unfiltered_rms_slope_rad < b.normal_detail.unfiltered_rms_slope_rad);
    m.bulk.fabric_to_material = [0.5_f64.sqrt(), 0., 0.5_f64.sqrt(), 0.];
    let rotated = surface_sample(&m, &s, 0.0001).unwrap();
    assert!(
        (rotated.normal_detail.unfiltered_rms_slope_rad - a.normal_detail.unfiltered_rms_slope_rad)
            .abs()
            < 1e-14
    );
    // Rest coordinates are metadata for homogeneous reference; no world/tile RNG exists.
    s.material_coordinate_m = [1e9, -2e8, 3.];
    assert_eq!(rotated, surface_sample(&m, &s, 0.0001).unwrap());
}
#[test]
fn bad_inputs_refused_without_default_rock_or_air() {
    let c = Catalogue::embedded().unwrap();
    assert!(matches!(c.rock("no_rock"), Err(Error::UnknownRock(_))));
    let mut m = material(&c, "granite");
    let s = fresh();
    for f in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(surface_sample(&m, &s, f).is_err());
    }
    m.bulk.mineral_mode[0].mass_fraction += 0.1;
    assert!(surface_sample(&m, &s, 1.).is_err());
    let m = material(&c, "granite");
    let mut s = fresh();
    s.normal_material = [0.; 3];
    assert!(surface_sample(&m, &s, 1.).is_err());
    s = fresh();
    s.water_film_m = 0.001;
    assert!(surface_sample(&m, &s, 1.).is_err());
    s = fresh();
    s.exposure_elapsed_s = Some(-1.);
    assert!(surface_sample(&m, &s, 1.).is_err());
}
#[test]
fn corrupted_catalogue_provenance_and_mass_fail_loudly() {
    let c = Catalogue::embedded().unwrap();
    let mut bad = c.clone();
    bad.minerals[0].density_kg_m3.interval = [3000., 1000.];
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.minerals[0].density_kg_m3.reference = "missing".into();
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.rocks[0].dry_bulk_density_kg_m3.value *= 0.9;
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.rocks[0].mode[0].mass_fraction.value -= 0.1;
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.rocks[0].vesicularity.value = 0.8;
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.minerals.push(bad.minerals[0].clone());
    assert!(bad.validate().is_err());
}

#[test]
fn glass_response_uses_constituent_and_fracture_slope_not_rock_key() {
    let c = Catalogue::embedded().unwrap();
    let s = fresh();
    let obsidian = surface_sample(&material(&c, "obsidian"), &s, 0.1).unwrap();
    let basalt = surface_sample(&material(&c, "basalt_tholeiitic"), &s, 0.1).unwrap();
    let pumice = surface_sample(&material(&c, "pumice"), &s, 0.1).unwrap();
    assert!(obsidian.roughness < basalt.roughness);
    assert!(pumice.base_reflectance_linear[1] > obsidian.base_reflectance_linear[1]);
    let mut m = material(&c, "basalt_tholeiitic");
    m.bulk.mineral_mode = material(&c, "obsidian").bulk.mineral_mode;
    m.grain_diameter_m = c.rock("obsidian").unwrap().optical_domain_diameter_m.value;
    let same = surface_sample(&m, &s, 0.1).unwrap();
    assert_eq!(
        same.base_reflectance_linear,
        obsidian.base_reflectance_linear
    );
    assert_eq!(same.roughness, obsidian.roughness);
}
#[test]
fn two_layer_detail_keeps_two_actual_wavelengths_and_variances() {
    let c = Catalogue::embedded().unwrap();
    let m = material(&c, "granite");
    let mut s = fresh();
    s.coating = Some(CoatingLayer {
        source_event: "two-band-test".into(),
        mineral_mode: vec![ConstituentFraction {
            mineral: "muscovite".into(),
            mass_fraction: 1.0,
        }],
        grain_diameter_m: 0.00002,
        thickness_m: 0.0001,
        coverage: 0.5,
        extinction_m_inv: 10000.0,
    });
    s.wet_coverage = 0.4;
    s.water_film_m = 0.00001;
    let a = surface_sample(&m, &s, 0.0001).unwrap();
    assert_eq!(a.normal_detail.wavelength_m, m.grain_diameter_m);
    assert_eq!(
        a.coating_normal_detail.as_ref().unwrap().wavelength_m,
        0.00002
    );
    for d in [&a.normal_detail, a.coating_normal_detail.as_ref().unwrap()] {
        assert!(
            (d.resolved_rms_slope_rad.powi(2) + d.unresolved_slope_variance
                - d.unfiltered_rms_slope_rad.powi(2))
            .abs()
                < 1e-12
        );
    }
    assert!(
        a.sparkle.unwrap().specular_energy_fraction
            + a.coating_sparkle.unwrap().specular_energy_fraction
            <= c.model.sparkle_fraction.value
    );
}
#[test]
fn intrinsic_and_deposit_pores_are_not_counted_twice() {
    let c = Catalogue::embedded().unwrap();
    let mut m = material(&c, "sand");
    let s = fresh();
    let a = surface_sample(&m, &s, 0.1).unwrap();
    assert!(
        (a.dry_deposit_density_kg_m3 / a.dry_bulk_density_kg_m3
            - (1.0 - m.bulk.interparticle_porosity))
            .abs()
            < 1e-12
    );
    m.bulk.interparticle_porosity = 0.0;
    let b = surface_sample(&m, &s, 0.1).unwrap();
    assert_eq!(a.dry_bulk_density_kg_m3, b.dry_bulk_density_kg_m3);
    assert_eq!(b.dry_bulk_density_kg_m3, b.dry_deposit_density_kg_m3);
    assert_eq!(a.base_reflectance_linear, b.base_reflectance_linear);
    assert_ne!(
        c.rock("gravel").unwrap().grain.d50_m.value,
        c.rock("gravel").unwrap().optical_domain_diameter_m.value
    );
}
