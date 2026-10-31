use pf_rockspec::{
    v3::{self, BodyAdmission, BodyKind, FormationContext, ResolvedBody},
    v3_2c::{Evaluator, VERSION},
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
fn all_bodies_passivity_far_equality_and_frozen_body_state() {
    let c = v3::Catalogue::embedded().unwrap();
    let mut samples = 0;
    let mut far = 0;
    for r in &c.realizations {
        let b = body(&c, &r.id, "test-body");
        let before = serde_json::to_vec(&b).unwrap();
        let e = Evaluator::new(&b).unwrap();
        let old = pf_rockspec::v3_2b::Evaluator::new(&b).unwrap();
        for f in [1e-6, 0.1 / 128., 2. / 128., 10. / 128., 1., 100.] {
            for i in 0..16 {
                let mut s = surface([i as f64 * 0.137 - 0.7, 0.017, 0.193 * i as f64 - 1.]);
                if i % 4 == 1 {
                    s.normal_material = [0., 0., 1.];
                }
                if i % 4 == 2 {
                    s.wet_coverage = 0.6;
                    s.water_film_m = 0.0001;
                }
                if i % 4 == 3 {
                    s.coating = Some(pf_rockspec::CoatingLayer {
                        source_event: "manufactured-coating".into(),
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
                let out = e.surface_sample(&s, f).unwrap();
                let st = out.structure_v2.as_ref().unwrap();
                assert_eq!(out.contract_version, VERSION);
                assert_eq!(out.material_ref, *b.reference());
                for j in 0..3 {
                    assert!((0. ..=1.).contains(&out.base_reflectance_linear[j]));
                    assert!(out.diffuse_energy_linear[j] + out.f0_linear[j] <= 1. + 1e-14);
                }
                assert!((0. ..=1.).contains(&out.roughness));
                assert!((st.normal_material.iter().map(|x| x * x).sum::<f64>() - 1.).abs() < 1e-12);
                if f == 100. {
                    let o = old.surface_sample(&s, f).unwrap();
                    assert_eq!(out.base_reflectance_linear, o.base_reflectance_linear);
                    assert_eq!(out.f0_linear, o.f0_linear);
                    assert_eq!(out.roughness, o.roughness);
                    assert_eq!(st.normal_material, o.structure_v2.unwrap().normal_material);
                    far += 1;
                }
                samples += 1;
            }
        }
        assert_eq!(before, serde_json::to_vec(&b).unwrap());
    }
    println!("v3.2c passivity/state samples={samples}; exact v3.2b far samples={far}");
}
fn peak_ratio(values: &[f64]) -> (f64, f64) {
    let n = values.len();
    let mean = values.iter().sum::<f64>() / n as f64;
    let mut total = 0.;
    let mut peak: f64 = 0.;
    for k in 1..=n / 2 {
        let angle = -std::f64::consts::TAU * k as f64 / n as f64;
        let (s, c) = angle.sin_cos();
        let (mut cs, mut sn) = (1., 0.);
        let (mut re, mut im) = (0., 0.);
        for &v in values {
            re += (v - mean) * cs;
            im += (v - mean) * sn;
            let next = cs * c - sn * s;
            sn = sn * c + cs * s;
            cs = next;
        }
        let power = re * re + im * im;
        total += power;
        peak = peak.max(power);
    }
    (
        peak / total,
        (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64).sqrt(),
    )
}
#[test]
fn outcrop_periodicity_falsifier_and_nonzero_structure() {
    let c = v3::Catalogue::embedded().unwrap();
    let mut worst: f64 = 0.;
    let mut min_std = f64::INFINITY;
    let mut tested = 0;
    for id in [
        "schist.reference",
        "phyllite.reference",
        "ignimbrite.welded",
        "ignimbrite.unwelded",
        "pumice.reference",
        "granite.reference",
        "basalt_tholeiitic.reference",
        "sandstone_quartz.reference",
    ] {
        for address in ["spectral-body-a", "spectral-body-b", "spectral-body-c"] {
            let b = body(&c, id, address);
            let e = Evaluator::new(&b).unwrap();
            let (index, law) = e
                .record()
                .outcrop_distributions
                .iter()
                .enumerate()
                .find(|(_, v)| v.kind == "joint")
                .unwrap();
            let mut values = Vec::new();
            for i in 0..2048 {
                let x = (i as f64 + 0.5) * 128. * law.spacing_m.value / 2048.;
                let p = std::array::from_fn(|j| {
                    law.normal[j] * x + law.tangent_u[j] * 0.371 * law.persistence_m.value
                        - law.tangent_v[j] * 0.217 * law.persistence_m.value
                });
                values.push(
                    e.trace_probe(index, p, law.tangent_v, law.spacing_m.value * 0.02)
                        .unwrap()
                        .value,
                );
            }
            let (ratio, std) = peak_ratio(&values);
            worst = worst.max(ratio);
            min_std = min_std.min(std);
            tested += 1;
            println!("spectrum {id} {address}: peak/total={ratio:.7}, std={std:.7}");
            assert!(ratio < 0.10, "periodic trace field {id} {address}: {ratio}");
            assert!(
                std > 0.005,
                "cannot pass by erasing trace field {id} {address}: {std}"
            );
        }
    }
    let periodic = (0..2048)
        .map(|i| (std::f64::consts::TAU * i as f64 / 16.).sin())
        .collect::<Vec<_>>();
    let (control, _) = peak_ratio(&periodic);
    assert!(control > 0.99);
    println!("spectral cases={tested}, maximum ratio={worst:.7}, minimum std={min_std:.7}; periodic control={control:.7}");
}
#[test]
fn frame_address_and_reload_are_material_bound() {
    let c = v3::Catalogue::embedded().unwrap();
    let mut count = 0;
    let mut changed = 0;
    for id in [
        "schist.reference",
        "granite.reference",
        "ignimbrite.welded",
        "sandstone_quartz.reference",
    ] {
        let b = body(&c, id, "frame-body");
        let e = Evaluator::new(&b).unwrap();
        let mut encoded = serde_json::to_value(&b).unwrap();
        encoded["bulk"]["state"]["fabric_to_material"] =
            serde_json::json!([0.5_f64.sqrt(), 0., 0., 0.5_f64.sqrt()]);
        let rotated: ResolvedBody = serde_json::from_value(encoded).unwrap();
        let re = Evaluator::new(&rotated).unwrap();
        let other = body(&c, id, "other-body");
        let oe = Evaluator::new(&other).unwrap();
        let reload = ResolvedBody::from_json(&serde_json::to_string(&b).unwrap()).unwrap();
        let le = Evaluator::new(&reload).unwrap();
        for i in 0..64 {
            let s = surface([i as f64 * 0.021 - 0.4, 0.017, 0.031 * i as f64 - 0.8]);
            let a = e.surface_sample(&s, 0.0001).unwrap();
            assert_eq!(a, le.surface_sample(&s, 0.0001).unwrap());
            let mut moved = s.clone();
            moved.surface_id = "fragment-in-another-tile".into();
            assert_eq!(
                a.structure_v2,
                e.surface_sample(&moved, 0.0001).unwrap().structure_v2
            );
            moved.material_coordinate_m = rot(s.material_coordinate_m);
            moved.normal_material = rot(s.normal_material);
            let r = re.surface_sample(&moved, 0.0001).unwrap();
            assert!(maxdiff(a.base_reflectance_linear, r.base_reflectance_linear) < 1e-10);
            assert!(
                maxdiff(
                    rot(a.structure_v2.unwrap().normal_material),
                    r.structure_v2.unwrap().normal_material
                ) < 1e-8
            );
            if maxdiff(
                a.base_reflectance_linear,
                oe.surface_sample(&s, 0.0001)
                    .unwrap()
                    .base_reflectance_linear,
            ) > 1e-8
            {
                changed += 1;
            }
            count += 1;
        }
    }
    assert!(changed > 30);
    println!("frame/reload/tile cases={count}, changed-body colours={changed}");
}
#[test]
fn fiamme_constituents_close_and_contrast_survives_dealignment() {
    let c = v3::Catalogue::embedded().unwrap();
    let b = body(&c, "ignimbrite.welded", "declared-comparison-body");
    let e = Evaluator::new(&b).unwrap();
    let old = pf_rockspec::v3_2b::Evaluator::new(&b).unwrap();
    let record = e.record().fiamme.as_ref().unwrap();
    let oldrecord = old.portable_structure();
    let mut error: f64 = 0.;
    for i in 0..record.lens_volume_mode.len() {
        let v = record.mean_lens_volume_fraction * record.lens_volume_mode[i]
            + (1. - record.mean_lens_volume_fraction) * record.matrix_volume_mode[i];
        error = error.max((v - oldrecord.phases[i].area_fraction).abs());
    }
    assert!(error < 1e-14);
    assert!(record.matrix_volume_mode.iter().all(|&v| v >= 0.));
    let mut values = [Vec::new(), Vec::new()];
    for y in 0..128 {
        for x in 0..128 {
            let s = surface([
                ((x as f64 + 0.5) / 128. - 0.5) * 0.1,
                0.017,
                ((y as f64 + 0.5) / 128. - 0.5) * 0.1,
            ]);
            for (i, out) in [
                old.surface_sample(&s, 0.1 / 128.).unwrap(),
                e.surface_sample(&s, 0.1 / 128.).unwrap(),
            ]
            .iter()
            .enumerate()
            {
                values[i].push(
                    out.base_reflectance_linear
                        .iter()
                        .zip([0.2126, 0.7152, 0.0722])
                        .map(|(v, w)| v * w)
                        .sum::<f64>(),
                );
            }
        }
    }
    let stats = |v: &mut Vec<f64>| {
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        let std = (v.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / v.len() as f64).sqrt();
        v.sort_by(f64::total_cmp);
        (mean, std, v[v.len() * 99 / 100] - v[v.len() / 100])
    };
    let a = stats(&mut values[0]);
    let n = stats(&mut values[1]);
    println!("fiamme mode volume closure max error={error:.3e}; lens fraction={}; optics lens={:?} matrix={:?}; old(mean,std,p99-p1)={a:?}; new={n:?}; std gain={}",record.mean_lens_volume_fraction,record.lens_reflectance,record.matrix_reflectance,n.1/a.1);
    assert!(
        n.1 > 0.8 * a.1,
        "jittered lenses retain optical contrast, without tint"
    );
    assert!(
        n.2 > 0.035,
        "resolved lens contrast must exceed a few linear encoding steps"
    );
}
fn random(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*state >> 11) as f64 / ((1_u64 << 53) as f64)
}
#[test]
fn statistical_bulk_mean_and_visible_outcrop_normal_survive() {
    let c = v3::Catalogue::embedded().unwrap();
    let n = 16384;
    for id in [
        "schist.reference",
        "phyllite.reference",
        "ignimbrite.welded",
        "ignimbrite.unwelded",
        "pumice.reference",
        "granite.reference",
        "basalt_tholeiitic.reference",
        "sandstone_quartz.reference",
    ] {
        let b = body(&c, id, "mean-body");
        let e = Evaluator::new(&b).unwrap();
        let mut rng = 20261025;
        let mut sum = [0.; 3];
        let mut sum2 = [0.; 3];
        for _ in 0..n {
            let p = std::array::from_fn(|_| 64000. * random(&mut rng) - 32000.);
            let o = e.surface_sample(&surface(p), 1e-5).unwrap();
            let mean = o.structure_v2.unwrap().mean_reflectance_linear;
            for j in 0..3 {
                let d = o.base_reflectance_linear[j] - mean[j];
                sum[j] += d;
                sum2[j] += d * d;
            }
        }
        let errors = sum.map(|v| v / n as f64);
        let se = std::array::from_fn::<_, 3, _>(|j| {
            ((sum2[j] / n as f64 - errors[j] * errors[j]) / n as f64).sqrt()
        });
        for j in 0..3 {
            assert!(
                errors[j].abs() < 6. * se[j] + 1e-6,
                "{id} channel {j}: mean={} SE={}",
                errors[j],
                se[j]
            );
        }
        let mut slope_sum = 0.;
        let mut slope_max: f64 = 0.;
        for i in 0..1024 {
            let s = surface([i as f64 * 0.037 - 10., 0.017, i as f64 * 0.029 - 10.]);
            let o = e.surface_sample(&s, 10. / 128.).unwrap();
            let diff = maxdiff(o.structure_v2.unwrap().normal_material, s.normal_material);
            slope_sum += diff * diff;
            slope_max = slope_max.max(diff);
        }
        let rms = (slope_sum / 1024.).sqrt();
        assert!(
            slope_max > 0.01 && rms > 0.001,
            "can't remove outcrop normal: {id} max={slope_max} rms={rms}"
        );
        println!("mean {id}: n={n}, error={errors:?}, SE={se:?}; 10m/128 normal max={slope_max:.6}, RMS={rms:.6}");
    }
}
#[test]
fn declared_distributions_track_body_descriptors_and_reject_invalid_queries() {
    let c = v3::Catalogue::embedded().unwrap();
    let mut sets = 0;
    for r in &c.realizations {
        let b = body(&c, &r.id, "law-body");
        let e = Evaluator::new(&b).unwrap();
        let laws = e
            .record()
            .outcrop_distributions
            .iter()
            .filter(|v| v.kind == "joint");
        for (law, j) in laws.zip(&b.spec().joint_sets) {
            assert_eq!(law.spacing_m.value, j.spacing_m.value);
            assert_eq!(law.persistence_m.value, j.persistence_m.value);
            assert_eq!(law.aperture_m.value, j.aperture_m.value);
            for (support, prior) in [
                (law.persistence_support_m, &j.persistence_m),
                (law.aperture_support_m, &j.aperture_m),
            ] {
                assert!(support[0] >= prior.interval[0] && support[1] <= prior.interval[1]);
                assert!(((support[0] + support[1]) / 2. - prior.value).abs() < 1e-10);
            }
            assert!(law.site_probability > 0. && law.site_probability <= 1.);
            sets += 1;
        }
    }
    let b = body(&c, "granite.reference", "invalid-body");
    let e = Evaluator::new(&b).unwrap();
    for f in [0., -1., f64::INFINITY, f64::NAN] {
        assert!(e.surface_sample(&surface([0.; 3]), f).is_err());
    }
    assert!(e
        .trace_probe(0, [f64::NAN, 0., 0.], [0., 1., 0.], 1.)
        .is_err());
    assert!(e.trace_probe(0, [0.; 3], [0.; 3], 1.).is_err());
    println!("{sets} joint distributions retain supplied means/supports; invalid queries refused");
}
