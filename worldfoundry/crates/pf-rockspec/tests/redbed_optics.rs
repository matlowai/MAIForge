use pf_rockspec::v3::{BodyAdmission, BodyKind};
use pf_rockspec::{redbed_v1, v3, v3_1, v3_2c, *};
use serde_json::{json, Value};

fn body(c: &v3::Catalogue) -> v3::ResolvedBody {
    let r = c.realization(redbed_v1::REALIZATION).unwrap();
    let reference = MaterialRef {
        spec_id: r.lithofacies.clone(),
        spec_version: c.core.schema_version.clone(),
        realization_id: r.id.clone(),
        body_id: "manufactured-coating-diagnostic".into(),
    };
    c.resolve(
        reference.clone(),
        v3::FormationContext {
            setting: r.settings[0].clone(),
            protolith: None,
            metamorphic_grade: None,
            burial_depth_m: None,
            provenance_ref: "manufactured-test-truth".into(),
        },
        BodyAdmission {
            material_ref: reference,
            kind: BodyKind::Host,
            event_id: "test".into(),
            snapshot_id: "test".into(),
            geometry_ref: "test".into(),
            stock_account_ref: "test".into(),
            material_frame_ref: "test".into(),
            host_body_id: None,
        },
    )
    .unwrap()
}
fn fresh() -> SurfaceState {
    SurfaceState::fresh("diagnostic", "manufactured-cut", [0.; 3], [0., 0., 1.])
}

/// Explicit output path only. Existing frozen comparison fixtures and optional
/// real saved native book. Never loads project-local files implicitly.
#[test]
fn export_optics_evidence() {
    let Ok(out) = std::env::var("VZ_REDBED_EVIDENCE_OUT") else {
        eprintln!("SKIP optional evidence export: VZ_REDBED_EVIDENCE_OUT unset");
        return;
    };
    let c = redbed_v1::catalogue().unwrap();
    let mut rows = vec![];
    for f in v3_1::FixtureSet::swatches().unwrap().samples {
        let b = f.resolve_fixture(&c).unwrap();
        let e = v3_2c::Evaluator::new(&b).unwrap();
        rows.push(
            json!({"fixture":f,"fresh_dry":e.surface_sample(&fresh(),2.).unwrap(),
            "transport":e.transport_inputs().unwrap(),"structure":e.portable_structure()}),
        );
    }
    let b = body(&c);
    let e = v3_2c::Evaluator::new(&b).unwrap();
    let mut native = vec![];
    if let Ok(path) = std::env::var("VZ_REDBED_NATIVE_BOOK") {
        let file = std::fs::File::open(path).unwrap();
        assert!(file.metadata().unwrap().len() < 200_000_000);
        let book: Value = serde_json::from_reader(file).unwrap();
        for (id, saved) in book["bodies"].as_object().unwrap() {
            if saved["resolved"].is_null() {
                continue;
            }
            let b: v3::ResolvedBody = serde_json::from_value(saved["resolved"].clone()).unwrap();
            let e = v3_2c::Evaluator::new(&b).unwrap();
            native.push(
                json!({"body":id,"sample":e.surface_sample(&fresh(),2.).unwrap(),
                "transport":e.transport_inputs().unwrap()}),
            );
        }
    }
    let result = json!({"legacy":rows,"redbed":e.surface_sample(&fresh(),2.).unwrap(),
        "redbed_transport":e.transport_inputs().unwrap(),"native":native});
    std::fs::write(out, serde_json::to_vec_pretty(&result).unwrap()).unwrap();
}

fn variant(hematite: f64, clay_rim: f64, diameter: f64) -> v3::ResolvedBody {
    let mut c = redbed_v1::catalogue().unwrap();
    let r = c
        .realizations
        .iter_mut()
        .find(|r| r.id == redbed_v1::REALIZATION)
        .unwrap();
    for p in &mut r.spec.mode {
        if p.mineral == "hematite" {
            p.mass_fraction.value = hematite;
            p.mass_fraction.interval = [0., 0.03];
        }
        if p.mineral == "quartz" {
            p.mass_fraction.value = 0.96 - hematite;
        }
    }
    let coat = &mut r.parameters.sedimentary.as_mut().unwrap().grain_coatings;
    coat[0].mass_fraction = hematite;
    coat[1].mass_fraction = clay_rim;
    r.spec.optical_domain_diameter_m.value = diameter;
    body(&c)
}
#[test]
fn mass_energy_determinism_and_parameter_sensitivity() {
    let mut n = 0;
    let mut evidence = vec![];
    for w in [0., 0.0001, 0.0005, 0.001, 0.002, 0.007, 0.01, 0.02, 0.03] {
        for clay in [0., 0.0125, 0.025] {
            for diameter in [0.00015, 0.0003, 0.0006] {
                let b = variant(w, clay, diameter);
                let before = serde_json::to_vec(&b).unwrap();
                let c = redbed_v1::CoatingReceipt::for_body(&b).unwrap().unwrap();
                assert!(
                    (c.coating_volume_m3_per_kg - c.reconstructed_coating_volume_m3_per_kg).abs()
                        < 1e-18
                );
                let e = v3_2c::Evaluator::new(&b).unwrap();
                for footprint in [0.00001, 0.002, 2.] {
                    let s = e.surface_sample(&fresh(), footprint).unwrap();
                    assert_eq!(s, e.surface_sample(&fresh(), footprint).unwrap());
                    for j in 0..3 {
                        assert!(c.layer_reflectance[j] >= 0. && c.layer_transmittance[j] >= 0.);
                        assert!(c.layer_reflectance[j] + c.layer_transmittance[j] <= 1. + 1e-12);
                        assert!(s.diffuse_energy_linear[j] + s.f0_linear[j] <= 1. + 1e-12);
                    }
                }
                assert_eq!(before, serde_json::to_vec(&b).unwrap());
                evidence.push(
                    json!({"hematite":w,"clay_rim":clay,"diameter_m":diameter,"receipt":c,
                    "sample":e.surface_sample(&fresh(),2.).unwrap()}),
                );
                n += 1;
            }
        }
    }
    assert_eq!(n, 81);
    let b = variant(0., 0., 0.0003);
    let c = redbed_v1::CoatingReceipt::for_body(&b).unwrap().unwrap();
    assert_eq!(c.thickness_m, 0.);
    assert_eq!(c.coated_reflectance, c.uncoated_host_reflectance);
    let a = redbed_v1::CoatingReceipt::for_body(&variant(0.007, 0.025, 0.0003))
        .unwrap()
        .unwrap();
    let z = redbed_v1::CoatingReceipt::for_body(&variant(0.007, 0.025, 0.0006))
        .unwrap()
        .unwrap();
    assert!((z.thickness_m / a.thickness_m - 2.).abs() < 1e-12);
    let mut sensitivity = vec![];
    for l in [1e-6, 2e-6, 2.5e-6, 3e-6, 10e-6] {
        sensitivity.push(
            redbed_v1::CoatingReceipt::with_scattering_length(&variant(0.007, 0.025, 0.0003), l)
                .unwrap(),
        );
    }
    if let Ok(out) = std::env::var("VZ_REDBED_SWEEP_OUT") {
        std::fs::write(
            out,
            serde_json::to_vec_pretty(
                &json!({"scenarios":evidence,"scattering_sensitivity":sensitivity}),
            )
            .unwrap(),
        )
        .unwrap();
    }
    println!("COATING_CASES={n}; ENERGY_SURFACES={}", 3 * n);
}
#[test]
fn invalid_stock_and_lengths_refused_and_old_key_off() {
    let mut b = serde_json::to_value(variant(0.007, 0.025, 0.0003)).unwrap();
    // Saved body is immutable through the public API; falsify the actual saved
    // parameter record via its wire format to exercise the trust boundary.
    let before = b.clone();
    let obj = b.as_object_mut().unwrap();
    let params = obj.get_mut("parameters").expect("saved body parameters");
    params["sedimentary"]["grain_coatings"][0]["mass_fraction"] = json!(0.1);
    let bad: v3::ResolvedBody = serde_json::from_value(b).unwrap();
    assert!(redbed_v1::CoatingReceipt::for_body(&bad).is_err());
    let good: v3::ResolvedBody = serde_json::from_value(before).unwrap();
    for length in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(redbed_v1::CoatingReceipt::with_scattering_length(&good, length).is_err());
    }
    let old = v3_1::catalogue().unwrap();
    for f in v3_1::FixtureSet::comparisons(&old).unwrap().samples {
        let b = f.resolve_fixture(&old).unwrap();
        assert!(redbed_v1::CoatingReceipt::for_body(&b).unwrap().is_none());
    }
}
