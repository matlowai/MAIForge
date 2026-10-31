use pf_rockspec::v2::*;
use pf_rockspec::*;
fn mat<'a>(c: &'a Catalogue, k: &str) -> Material<'a> {
    c.realize_prior(
        k,
        MaterialRef {
            spec_id: k.into(),
            spec_version: c.schema_version.clone(),
            realization_id: "falsifier-v2".into(),
            body_id: "body".into(),
        },
    )
    .unwrap()
}
fn fresh() -> SurfaceState {
    SurfaceState::fresh("face", "cut", [0.123, 0.456, 0.789], [0., 0., 1.])
}
fn y(o: &SurfaceSample) -> f64 {
    o.base_reflectance_linear
        .iter()
        .zip([0.2126, 0.7152, 0.0722])
        .map(|(a, b)| a * b)
        .sum()
}
#[test]
fn v2_catalogue_and_all_frozen_states_are_valid_and_reproducible() {
    let c = Catalogue::embedded_v2().unwrap();
    assert_eq!(c.rocks.len(), 43);
    assert_eq!(c.minerals.len(), 33);
    let cases = DiagnosticSet::frozen_v2().unwrap();
    assert_eq!(cases.samples.len(), 128);
    for x in &cases.samples {
        let a = x.evaluate(&c).unwrap();
        assert_eq!(a, x.evaluate(&c).unwrap());
        assert!(a.structure_v2.is_some());
        for j in 0..3 {
            assert!(a.f0_linear[j] + a.diffuse_energy_linear[j] <= 1. + 1e-12);
        }
    }
    let sw = DiagnosticSet::swatches_v2().unwrap();
    assert_eq!(sw.samples.len(), 387);
    for x in sw.samples {
        x.evaluate(&c).unwrap();
    }
}
#[test]
fn protected_fresh_materials_stay_in_arc1_bands_at_all_three_footprints() {
    let c = Catalogue::embedded_v2().unwrap();
    for (key, lo, hi) in [
        ("basalt_tholeiitic", 0.05, 0.10),
        ("granite", 0.25, 0.35),
        ("obsidian", 0.008, 0.03),
        ("marble", 0.6, 0.9),
    ] {
        for f in [0.1, 2., 50.] {
            let value = y(&surface_sample(&mat(&c, key), &fresh(), f).unwrap());
            assert!((lo..=hi).contains(&value), "{key} {f}m Y={value}");
        }
    }
}
#[test]
fn sediments_do_not_return_to_clean_powder_brightness() {
    let c = Catalogue::embedded_v2().unwrap();
    for (key, hi) in [
        ("shale", 0.20),
        ("mudstone", 0.22),
        ("siltstone", 0.30),
        ("sandstone_quartz", 0.50),
        ("arkose", 0.40),
        ("conglomerate", 0.40),
        ("chert", 0.35),
        ("sand", 0.50),
        ("silt_clay", 0.35),
        ("till", 0.30),
        ("regolith", 0.35),
        ("quartzite", 0.60),
        ("pumice", 0.60),
    ] {
        let value = y(&surface_sample(&mat(&c, key), &fresh(), 2.).unwrap());
        assert!(value < hi, "{key} Y={value}");
    }
    // Original shale already had carbon. This falsifies the "just add carbon" explanation.
    let c1 = Catalogue::embedded().unwrap();
    assert_eq!(
        c1.rock("shale")
            .unwrap()
            .mode
            .iter()
            .find(|x| x.mineral == "organic_matter")
            .unwrap()
            .mass_fraction
            .value,
        0.05
    );
    let mut powder = c.clone();
    powder
        .rocks
        .iter_mut()
        .find(|r| r.id == "shale")
        .unwrap()
        .structure_v2
        .as_mut()
        .unwrap()
        .optical_boundary_efficiency
        .value = 1.;
    assert!(
        y(&surface_sample(&mat(&powder, "shale"), &fresh(), 2.).unwrap())
            > 2. * y(&surface_sample(&mat(&c, "shale"), &fresh(), 2.).unwrap())
    );
}
#[test]
fn exposure_is_explicit_family_specific_and_keeps_coal_dark() {
    let c = Catalogue::embedded_v2().unwrap();
    for r in &c.rocks {
        let m = mat(&c, &r.id);
        let before = serde_json::to_vec(&m.bulk).unwrap();
        let f = fresh();
        let dry = surface_sample(&m, &f, 2.).unwrap();
        let wet_state =
            diagnostic_weathering(&m, f.clone(), ExposureScenario::WetTemperate).unwrap();
        let arid_state =
            diagnostic_weathering(&m, f.clone(), ExposureScenario::AridVarnish).unwrap();
        let wet = surface_sample(&m, &wet_state, 2.).unwrap();
        let arid = surface_sample(&m, &arid_state, 2.).unwrap();
        assert_eq!(serde_json::to_vec(&m.bulk).unwrap(), before);
        assert_eq!(
            wet.substrate.diffuse_reflectance_linear,
            dry.substrate.diffuse_reflectance_linear
        );
        if r.id == "coal" {
            assert!(y(&wet) < 0.06 && y(&arid) < 0.06);
        } else {
            assert_ne!(wet_state.coating, arid_state.coating);
        }
        if r.id == "granite" || r.id == "gneiss" {
            assert!(y(&wet) - y(&dry) < 0.10);
            assert!(y(&arid) < y(&dry));
        }
        if r.id == "basalt_tholeiitic" {
            assert!(y(&wet) - y(&dry) < 0.08);
            assert!(
                wet.base_reflectance_linear[0] / wet.base_reflectance_linear[2]
                    > dry.base_reflectance_linear[0] / dry.base_reflectance_linear[2]
            );
        }
        let mut aged = f.clone();
        aged.exposure_elapsed_s = Some(1e12);
        assert_eq!(dry, surface_sample(&m, &aged, 2.).unwrap());
        if wet_state.coating.is_some() {
            let mut zero = wet_state;
            zero.coating.as_mut().unwrap().thickness_m = 0.;
            assert_eq!(
                surface_sample(&m, &zero, 2.)
                    .unwrap()
                    .base_reflectance_linear,
                dry.base_reflectance_linear
            );
        }
    }
}
#[test]
fn scale_changes_come_from_joint_descriptors_and_not_ids() {
    let c = Catalogue::embedded_v2().unwrap();
    let m = mat(&c, "granite");
    let a = surface_sample(&m, &fresh(), 0.1).unwrap();
    let b = surface_sample(&m, &fresh(), 2.).unwrap();
    let far = surface_sample(&m, &fresh(), 50.).unwrap();
    assert!(y(&a) > y(&b) && y(&b) > y(&far));
    assert!(a.roughness < b.roughness);
    assert_eq!(a.dry_bulk_density_kg_m3, far.dry_bulk_density_kg_m3);
    let mut no_joints = c.clone();
    no_joints
        .rocks
        .iter_mut()
        .find(|r| r.id == "granite")
        .unwrap()
        .joint_sets
        .clear();
    let flat = surface_sample(&mat(&no_joints, "granite"), &fresh(), 50.).unwrap();
    assert!(y(&flat) > y(&far));
    let mut renamed = c.clone();
    renamed
        .rocks
        .iter_mut()
        .find(|r| r.id == "granite")
        .unwrap()
        .id = "same-material-new-label".into();
    let same = surface_sample(&mat(&renamed, "same-material-new-label"), &fresh(), 50.).unwrap();
    assert_eq!(same.base_reflectance_linear, far.base_reflectance_linear);
    for r in [0., 0.02, 0.5, 1.] {
        assert!((cavity_reflectance(r, 1.) - r).abs() < 1e-12);
        assert!(cavity_reflectance(r, 0.4) <= r + 1e-12);
    }
}
#[test]
fn texture_is_material_addressed_and_far_field_is_homogenized() {
    let c = Catalogue::embedded_v2().unwrap();
    for key in [
        "gneiss",
        "conglomerate",
        "gravel",
        "till",
        "schist",
        "scoria",
        "pumice",
        "andesite",
        "dacite",
        "shale",
        "siltstone",
        "granite",
    ] {
        let m = mat(&c, key);
        let mut s = fresh();
        let d = m.spec.structure_v2.as_ref().unwrap().wavelength_m.value;
        let mut values = Vec::new();
        let mut normal_values = Vec::new();
        for i in 0..128 {
            s.material_coordinate_m = [
                i as f64 * d * 0.13,
                i as f64 * d * 0.31,
                i as f64 * d * 0.57,
            ];
            let v = surface_sample(&m, &s, d / 8.).unwrap();
            values.push(y(&v));
            normal_values.push(v.structure_v2.unwrap().normal_material);
        }
        let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            hi - lo > 0.001 || normal_values.windows(2).any(|w| w[0] != w[1]),
            "{key} missing structure"
        );
        let a = surface_sample(&m, &s, d / 8.).unwrap();
        s.surface_id = "different-tile-or-fragment-face".into();
        s.creation_event = "new-cut-same-bulk".into();
        let b = surface_sample(&m, &s, d / 8.).unwrap();
        assert_eq!(a.base_reflectance_linear, b.base_reflectance_linear);
        assert_eq!(a.structure_v2, b.structure_v2);
        let far = surface_sample(&m, &s, 50.).unwrap();
        s.material_coordinate_m = [99., -15., 30.];
        assert_eq!(
            far.base_reflectance_linear,
            surface_sample(&m, &s, 50.).unwrap().base_reflectance_linear
        );
        assert_eq!(
            far.base_reflectance_linear,
            far.structure_v2.unwrap().mean_reflectance_linear
        );
    }
}
#[test]
fn texture_rotation_and_filtered_area_mean() {
    let c = Catalogue::embedded_v2().unwrap();
    for key in [
        "granite",
        "gneiss",
        "conglomerate",
        "schist",
        "pumice",
        "shale",
    ] {
        let mut m = mat(&c, key);
        let d = m.spec.structure_v2.as_ref().unwrap().wavelength_m.value;
        let mut s = fresh();
        s.material_coordinate_m = [0.3 * d, 0.6 * d, 0.8 * d];
        let a = surface_sample(&m, &s, d / 10.).unwrap();
        // Exact 180 degree fabric/body rotation, avoids coordinate-boundary roundoff.
        m.bulk.fabric_to_material = [0., 0., 1., 0.];
        s.material_coordinate_m = [
            -s.material_coordinate_m[0],
            s.material_coordinate_m[1],
            -s.material_coordinate_m[2],
        ];
        s.normal_material = [0., 0., -1.];
        let b = surface_sample(&m, &s, d / 10.).unwrap();
        for j in 0..3 {
            assert!((a.base_reflectance_linear[j] - b.base_reflectance_linear[j]).abs() < 1e-12);
        }
        m.bulk.fabric_to_material = [1., 0., 0., 0.];
        s.normal_material = [0., 0., 1.];
        let mut observed = 0.;
        let mut expected = 0.;
        // Independent locations across 32^3 grain domains, not same-code replay.
        for i in 0..4096 {
            s.material_coordinate_m = [
                (i % 16) as f64 * d * 2.13,
                ((i / 16) % 16) as f64 * d * 2.31,
                (i / 256) as f64 * d * 2.57,
            ];
            let out = surface_sample(&m, &s, d / 16.).unwrap();
            observed += y(&out);
            let st = out.structure_v2.unwrap();
            expected += st
                .mean_reflectance_linear
                .iter()
                .zip([0.2126, 0.7152, 0.0722])
                .map(|(a, b)| a * b)
                .sum::<f64>();
        }
        let error = ((observed - expected) / 4096.).abs();
        assert!(error < 0.02, "{key} area mean error={error}");
    }
}
#[test]
fn v2_wetness_energy_and_malformed_structure_falsifiers() {
    let c = Catalogue::embedded_v2().unwrap();
    for r in &c.rocks {
        let m = mat(&c, &r.id);
        let s = fresh();
        let a = surface_sample(&m, &s, 2.).unwrap();
        let mut wet = s;
        wet.wet_coverage = 1.;
        wet.water_film_m = 0.0001;
        let b = surface_sample(&m, &wet, 2.).unwrap();
        assert!(y(&b) < y(&a));
        assert!(b.roughness < a.roughness);
    }
    assert!(
        surface_sample(&mat(&c, "pumice"), &fresh(), 2.)
            .unwrap()
            .roughness
            > 0.8
    );
    let mut bad = c.clone();
    bad.rocks[0]
        .structure_v2
        .as_mut()
        .unwrap()
        .optical_boundary_efficiency
        .value = 0.;
    assert!(bad.validate().is_err());
    let mut bad = c.clone();
    bad.rocks[0].structure_v2 = None;
    assert!(bad.validate().is_err());
    let mut bad = c;
    bad.rocks[0]
        .structure_v2
        .as_mut()
        .unwrap()
        .wavelength_m
        .value = f64::NAN;
    assert!(bad.validate().is_err());
}
