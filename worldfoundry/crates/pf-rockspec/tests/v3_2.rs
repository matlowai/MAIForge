use pf_rockspec::{
    v3::{self, BodyAdmission, BodyKind, FormationContext, ResolvedBody},
    v3_2::{Evaluator, VERSION},
    MaterialRef, SurfaceState,
};
fn body(c: &v3::Catalogue, id: &str, address: &str) -> ResolvedBody {
    let r = c.realization(id).unwrap();
    let reference = MaterialRef {
        spec_id: r.lithofacies.clone(),
        spec_version: c.core.schema_version.clone(),
        realization_id: id.into(),
        body_id: address.into(),
    };
    let admission = BodyAdmission {
        material_ref: reference.clone(),
        kind: BodyKind::Host,
        event_id: "manufactured-test".into(),
        snapshot_id: "fixture".into(),
        geometry_ref: "fixture".into(),
        stock_account_ref: "fixture-no-world-stock".into(),
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
            provenance_ref: "manufactured-test".into(),
        },
        admission,
    )
    .unwrap()
}
fn surface(p: [f64; 3]) -> SurfaceState {
    SurfaceState::fresh("face", "cut", p, [0., 1., 0.])
}
fn maxdiff(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter()
        .zip(b)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max)
}
fn rot(v: [f64; 3]) -> [f64; 3] {
    [-v[1], v[0], v[2]]
}
#[test]
fn all_realizations_all_scales_bounded_and_far_bulk_exact() {
    let c = v3::Catalogue::embedded().unwrap();
    let mut count = 0;
    let mut far = 0;
    for r in &c.realizations {
        let b = body(&c, &r.id, "test-body");
        let e = Evaluator::new(&b).unwrap();
        let before = serde_json::to_vec(&b).unwrap();
        for footprint in [
            1e-6,
            0.001 / 128.,
            0.1 / 128.,
            2. / 128.,
            10. / 128.,
            0.1,
            2.,
            100.,
        ] {
            for i in 0..16 {
                let mut s = surface([0.137 * i as f64 - 0.4, 0.017, 0.093 * i as f64 - 0.2]);
                if i % 4 == 1 {
                    s.normal_material = [0., 0., 1.];
                }
                if i % 4 == 2 {
                    s.wet_coverage = 0.6;
                    s.water_film_m = 0.0001;
                }
                if i % 4 == 3 {
                    s.coating = Some(pf_rockspec::CoatingLayer {
                        source_event: "test-coating".into(),
                        mineral_mode: vec![pf_rockspec::ConstituentFraction {
                            mineral: "kaolinite".into(),
                            mass_fraction: 1.,
                        }],
                        grain_diameter_m: 0.0001,
                        thickness_m: 0.0002,
                        coverage: 0.8,
                        extinction_m_inv: 12000.,
                    });
                }
                let out = e.surface_sample(&s, footprint).unwrap();
                let st = out.structure_v2.as_ref().unwrap();
                assert_eq!(out.contract_version, VERSION);
                assert_eq!(out.material_ref, *b.reference());
                for j in 0..3 {
                    assert!((0. ..=1.).contains(&out.base_reflectance_linear[j]));
                    assert!(out.diffuse_energy_linear[j] + out.f0_linear[j] <= 1. + 1e-14);
                }
                assert!(out.roughness.is_finite() && (0. ..=1.).contains(&out.roughness));
                assert!((st.normal_material.iter().map(|x| x * x).sum::<f64>() - 1.).abs() < 1e-12);
                assert!(
                    st.normal_material
                        .iter()
                        .zip(s.normal_material)
                        .map(|(a, b)| a * b)
                        .sum::<f64>()
                        > 0.
                );
                if footprint == 100. {
                    let old = b.surface_sample(&s, footprint).unwrap();
                    assert_eq!(
                        out.base_reflectance_linear, old.base_reflectance_linear,
                        "{}",
                        r.id
                    );
                    assert_eq!(out.roughness, old.roughness);
                    assert_eq!(out.f0_linear, old.f0_linear);
                    assert_eq!(
                        st.normal_material,
                        old.structure_v2.unwrap().normal_material
                    );
                    far += 1;
                }
                count += 1;
            }
        }
        assert_eq!(before, serde_json::to_vec(&b).unwrap());
    }
    println!("catalogue={} optical samples={count}, exact far samples={far}; dry/wet/coated and two face directions",c.realizations.len());
}
#[test]
fn material_address_and_frame_survive_fragment_and_view_changes() {
    let c = v3::Catalogue::embedded().unwrap();
    let mut tested = 0;
    let mut changed = 0;
    for id in [
        "schist.reference",
        "phyllite.reference",
        "gneiss.reference",
        "pumice.reference",
        "ignimbrite.welded",
        "granite.reference",
    ] {
        let b = body(&c, id, "parent-body");
        let e = Evaluator::new(&b).unwrap();
        let clone = ResolvedBody::from_json(&serde_json::to_string(&b).unwrap()).unwrap();
        let ce = Evaluator::new(&clone).unwrap();
        let other = body(&c, id, "other-body");
        let oe = Evaluator::new(&other).unwrap();
        let mut encoded = serde_json::to_value(&b).unwrap();
        encoded["bulk"]["state"]["fabric_to_material"] =
            serde_json::json!([0.5_f64.sqrt(), 0., 0., 0.5_f64.sqrt()]);
        let rotated: ResolvedBody = serde_json::from_value(encoded).unwrap();
        let re = Evaluator::new(&rotated).unwrap();
        for i in 0..64 {
            let s = surface([i as f64 * 0.00131 - 0.04, 0.017, 0.00197 * i as f64 - 0.05]);
            let a = e.surface_sample(&s, 1e-5).unwrap();
            let same = ce.surface_sample(&s, 1e-5).unwrap();
            assert_eq!(a, same);
            let mut moved = s.clone();
            moved.surface_id = "detached-fragment-face-in-another-tile".into();
            moved.creation_event = "fracture".into();
            let v = e.surface_sample(&moved, 1e-5).unwrap();
            assert_eq!(a.base_reflectance_linear, v.base_reflectance_linear);
            assert_eq!(a.structure_v2, v.structure_v2);
            if maxdiff(
                a.base_reflectance_linear,
                oe.surface_sample(&s, 1e-5).unwrap().base_reflectance_linear,
            ) > 1e-8
            {
                changed += 1;
            }
            moved.material_coordinate_m = rot(s.material_coordinate_m);
            moved.normal_material = rot(s.normal_material);
            let r = re.surface_sample(&moved, 1e-5).unwrap();
            assert!(
                maxdiff(a.base_reflectance_linear, r.base_reflectance_linear) < 1e-11,
                "{id}"
            );
            assert!(
                maxdiff(
                    rot(a.structure_v2.unwrap().normal_material),
                    r.structure_v2.unwrap().normal_material
                ) < 1e-9,
                "{id}"
            );
            tested += 1;
        }
    }
    assert!(changed > 100);
    println!("frame/reload/fragment cases={tested}; changed-address colour differences={changed}");
}
#[test]
fn foliation_lenses_and_outcrop_have_physical_scale() {
    let c = v3::Catalogue::embedded().unwrap();
    for id in ["schist.reference", "phyllite.reference"] {
        let b = body(&c, id, "test-body");
        let e = Evaluator::new(&b).unwrap();
        let r = e.record();
        assert!(r.foliated_cells);
        assert!(
            (r.cell_axes_m[0] / r.cell_axes_m[2] - b.spec().grain.aspect_ratio.value).abs() < 1e-12
        );
        let mut horizontal = 0.;
        let mut vertical = 0.;
        let d = b.spec().grain.d50_m.value;
        for i in 0..4096 {
            let p = [
                (i % 64) as f64 * d * 0.337,
                0.073,
                (i / 64) as f64 * d * 0.419,
            ];
            let a = e
                .surface_sample(&surface(p), d / 100.)
                .unwrap()
                .base_reflectance_linear;
            let mut px = p;
            let mut pz = p;
            px[0] += d * 0.08;
            pz[2] += d * 0.08;
            horizontal += maxdiff(
                a,
                e.surface_sample(&surface(px), d / 100.)
                    .unwrap()
                    .base_reflectance_linear,
            );
            vertical += maxdiff(
                a,
                e.surface_sample(&surface(pz), d / 100.)
                    .unwrap()
                    .base_reflectance_linear,
            );
        }
        assert!(
            vertical > 1.7 * horizontal,
            "{id} x={horizontal} z={vertical}"
        );
        println!(
            "{id} foliation vertical/horizontal finite-difference variation={:.5}",
            vertical / horizontal
        );
    }
    let b = body(&c, "ignimbrite.welded", "test-body");
    let e = Evaluator::new(&b).unwrap();
    assert!(e.record().welded_fiamme);
    assert_eq!(e.record().fiamme_flattening.value, 12.);
    let b = body(&c, "ignimbrite.unwelded", "test-body");
    assert!(!Evaluator::new(&b).unwrap().record().welded_fiamme);
    for id in [
        "granite.reference",
        "schist.reference",
        "gneiss.reference",
        "shale.reference",
    ] {
        let b = body(&c, id, "test-body");
        let e = Evaluator::new(&b).unwrap();
        let mut variation: f64 = 0.;
        for i in 0..1024 {
            let s = surface([i as f64 * 0.017, 0.17, i as f64 * 0.021]);
            let a = e.surface_sample(&s, 10. / 128.).unwrap();
            variation = variation.max(maxdiff(
                a.base_reflectance_linear,
                a.structure_v2.unwrap().mean_reflectance_linear,
            ));
        }
        assert!(variation > 1e-5, "outcrop missing {id}");
        println!("{id} 10m/128 footprint outcrop max reflectance departure={variation:.8}");
    }
}
#[test]
fn v31_explicit_opt_in_invalid_inputs_and_single_normal_sparkle() {
    let c = pf_rockspec::v3_1::catalogue().unwrap();
    // v3.1 catalogue is accessed through its body resolver; no schema mutation.
    let r = c.realization("schist.reference").unwrap();
    let b = body(&c, &r.id, "test-body");
    let e = Evaluator::new(&b).unwrap();
    let mut selected = 0;
    let mut varied = 0.;
    for i in 0..2048 {
        let s = surface([i as f64 * 0.000317, 0.013, 0.000173 * i as f64]);
        let o = e.surface_sample(&s, 1e-6).unwrap();
        let st = o.structure_v2.unwrap();
        assert_eq!(st.flake_specular_fraction, 0.);
        if st.qualification.contains("selected=true") {
            selected += 1;
            varied += maxdiff(st.normal_material, s.normal_material);
        }
    }
    assert!(selected > 5 && selected < 600 && varied > 1.);
    for f in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(e.surface_sample(&surface([0.; 3]), f).is_err());
    }
    assert!(e
        .surface_sample(&surface([f64::INFINITY, 0., 0.]), 1e-5)
        .is_err());
    assert!(e.surface_sample(&surface([1e20, 0., 0.]), 1e-5).is_err());
    println!("v3.1b schist single-normal flakes={selected}/2048; aggregate normal variation={varied:.5}; invalid footprint/coordinate refused");
}
fn random(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*state >> 11) as f64 / ((1_u64 << 53) as f64)
}
#[test]
fn resolved_substrate_mean_is_preserved_statistically() {
    let c = v3::Catalogue::embedded().unwrap();
    let n = 32768;
    for id in [
        "schist.reference",
        "phyllite.reference",
        "gneiss.reference",
        "banded_iron_formation.reference",
        "pumice.reference",
        "scoria.reference",
        "ignimbrite.welded",
        "granite.reference",
    ] {
        let b = body(&c, id, "test-body");
        let e = Evaluator::new(&b).unwrap();
        let mut rng = 20261025;
        let mut sum = [0.; 3];
        let mut sum2 = [0.; 3];
        for _ in 0..n {
            // Wide volume supplies many independent band IDs. A 64 m box repeats
            // only ~1067 gneiss bands, so 32768 points are NOT 32768 independent
            // phase draws. Preserve small-patch sampling variability in the figures.
            let p = std::array::from_fn(|_| 64000. * random(&mut rng) - 32000.);
            let out = e.surface_sample(&surface(p), 1e-5).unwrap();
            let st = out.structure_v2.unwrap();
            for j in 0..3 {
                let d = out.base_reflectance_linear[j] - st.mean_reflectance_linear[j];
                sum[j] += d;
                sum2[j] += d * d;
            }
        }
        let means = sum.map(|v| v / n as f64);
        let se = std::array::from_fn::<_, 3, _>(|j| {
            ((sum2[j] / n as f64 - means[j] * means[j]) / n as f64).sqrt()
        });
        for j in 0..3 {
            assert!(
                means[j].abs() < 6. * se[j] + 1e-6,
                "{id} channel {j} mean error={} SE={}",
                means[j],
                se[j]
            );
        }
        println!("{id} mean closure n={n} error={means:?} SE={se:?}");
    }
}
#[test]
fn welded_structure_is_lenticular_not_infinite_bands() {
    let c = v3::Catalogue::embedded().unwrap();
    let b = body(&c, "ignimbrite.welded", "test-body");
    let e = Evaluator::new(&b).unwrap();
    let d = b.spec().structure_v2.as_ref().unwrap().wavelength_m.value;
    let mut rng = 41;
    let mut xdiff = 0.;
    let mut zdiff = 0.;
    let mut nonzero = 0;
    for _ in 0..8192 {
        let p = [
            d * (30. * random(&mut rng) - 15.),
            d * (30. * random(&mut rng) - 15.),
            d * (30. * random(&mut rng) - 15.),
        ];
        let a = e
            .surface_sample(&surface(p), 1e-6)
            .unwrap()
            .base_reflectance_linear;
        let mut px = p;
        let mut pz = p;
        px[0] += d * 0.01;
        pz[2] += d * 0.01;
        let dx = maxdiff(
            a,
            e.surface_sample(&surface(px), 1e-6)
                .unwrap()
                .base_reflectance_linear,
        );
        let dz = maxdiff(
            a,
            e.surface_sample(&surface(pz), 1e-6)
                .unwrap()
                .base_reflectance_linear,
        );
        xdiff += dx;
        zdiff += dz;
        if dx > 1e-10 {
            nonzero += 1;
        }
    }
    assert!(
        nonzero > 100,
        "infinite straight bands have no lateral lens termination"
    );
    assert!(
        zdiff > 4. * xdiff,
        "flattened lenses not expressed: {xdiff} {zdiff}"
    );
    println!(
        "welded lens n=8192 lateral changes={nonzero}, vertical/horizontal variation={:.5}",
        zdiff / xdiff
    );
}
#[test]
fn thin_film_does_not_erase_outcrop_when_microstructure_is_unresolved() {
    let c = v3::Catalogue::embedded().unwrap();
    let b = body(&c, "granite.reference", "film-body");
    let e = Evaluator::new(&b).unwrap();
    let mut strongest: f64 = 0.;
    let mut relative_error: f64 = 0.;
    for i in 0..512 {
        let s = surface([i as f64 * 0.013, 0.017, 0.023]);
        let dry = e
            .surface_sample(&s, 0.04)
            .unwrap()
            .structure_v2
            .unwrap()
            .normal_material;
        let mut wet = s.clone();
        wet.wet_coverage = 1.;
        wet.water_film_m = 0.0001;
        let n = e
            .surface_sample(&wet, 0.04)
            .unwrap()
            .structure_v2
            .unwrap()
            .normal_material;
        let strength = maxdiff(dry, s.normal_material);
        if strength > 0.01 {
            strongest = strongest.max(strength);
            relative_error = relative_error.max(maxdiff(dry, n) / strength);
        }
    }
    assert!(strongest > 0.01);
    assert!(relative_error < 0.001);
    println!("thin film at 40mm footprint: maximum outcrop normal departure={strongest:.7}, relative change={relative_error:.7}");
}
#[test]
fn outcrop_uses_declared_spacing_and_requires_a_plane_intersection() {
    let c = v3::Catalogue::embedded().unwrap();
    let original = body(&c, "granite.reference", "joint-fixture");
    // Manufactured single-plane fixture, derived from the real supplied joint.
    let mut encoded = serde_json::to_value(&original).unwrap();
    let joint = encoded["catalogue"]["rocks"][0]["joint_sets"][0].clone();
    encoded["catalogue"]["rocks"][0]["joint_sets"] = serde_json::json!([joint]);
    let b: ResolvedBody = serde_json::from_value(encoded.clone()).unwrap();
    let e = Evaluator::new(&b).unwrap();
    for field in ["spacing_m", "persistence_m", "aperture_m"] {
        let estimate = &mut encoded["catalogue"]["rocks"][0]["joint_sets"][0][field];
        estimate["value"] = serde_json::json!(estimate["value"].as_f64().unwrap() * 4.);
        for j in 0..2 {
            estimate["interval"][j] =
                serde_json::json!(estimate["interval"][j].as_f64().unwrap() * 4.);
        }
        estimate["conditions"] =
            serde_json::json!("manufactured 4x scale-covariance test, not admitted field data");
    }
    let scaled: ResolvedBody = serde_json::from_value(encoded).unwrap();
    let se = Evaluator::new(&scaled).unwrap();
    let mut nonflat = 0;
    for i in 0..256 {
        let s = surface([i as f64 * 0.027 - 2., 0.017, 0.031]);
        let a = e.surface_sample(&s, 0.04).unwrap();
        let st = a.structure_v2.as_ref().unwrap();
        if maxdiff(st.normal_material, s.normal_material) > 0.01 {
            nonflat += 1;
        }
        let mut ss = s.clone();
        ss.material_coordinate_m = ss.material_coordinate_m.map(|v| 4. * v);
        let v = se.surface_sample(&ss, 0.16).unwrap();
        assert!(maxdiff(st.normal_material, v.structure_v2.unwrap().normal_material) < 1e-12);
        let mut parallel = s;
        parallel.normal_material = [1., 0., 0.];
        let p = e.surface_sample(&parallel, 0.04).unwrap();
        let pst = p.structure_v2.unwrap();
        assert_eq!(p.base_reflectance_linear, pst.mean_reflectance_linear);
        assert_eq!(pst.normal_material, parallel.normal_material);
    }
    assert!(nonflat > 10);
    println!("single-joint fixture: 256/256 parallel traces absent, 256/256 4x scale normals equal; nonflat={nonflat}");
}
