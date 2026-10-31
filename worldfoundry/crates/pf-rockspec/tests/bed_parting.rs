use pf_rockspec::redbed_v1::bed_parting_v1 as law;
use pf_rockspec::{v3, v3_2c, SurfaceState};
use serde_json::{json, Value};
fn source() -> Option<(Value, Value)> {
    let Ok(path) = std::env::var("VZ_BED_PROPERTY_BOOK") else {
        eprintln!("SKIP VZ_BED_PROPERTY_BOOK unset");
        return None;
    };
    assert!(std::fs::metadata(&path).unwrap().len() < 200_000_000);
    let v: Value = serde_json::from_reader(std::fs::File::open(path).unwrap()).unwrap();
    let b = v["bodies"]
        .as_array()
        .unwrap()
        .iter()
        .find(|b| b["stack"]["parent_prov_id"] == 15607552)
        .unwrap()
        .clone();
    Some((v, b))
}
#[test]
fn real_beds_material_energy_retention_and_export() {
    let Some((v, b)) = source() else {
        return;
    };
    let mut rows = Vec::new();
    let fresh = SurfaceState::fresh(
        "bed-diagnostic",
        "conditional-substrate",
        [0.; 3],
        [0., 0., 1.],
    );
    let mut covered = 0.;
    for ((band, beds), props) in b["stack"]["bands"]
        .as_array()
        .unwrap()
        .iter()
        .zip(b["bed_contacts"]["bands"].as_array().unwrap())
        .zip(b["bed_properties"]["bands"].as_array().unwrap())
    {
        let id = band["id"].as_str().unwrap();
        let base: v3::ResolvedBody =
            serde_json::from_value(v["materials"][id]["resolved"].clone()).unwrap();
        let before = serde_json::to_vec(&base).unwrap();
        for (bed, p) in beds["beds"]
            .as_array()
            .unwrap()
            .iter()
            .zip(props.as_array().unwrap())
        {
            let delta = law::Delta {
                hematite_mass_fraction: p["hematite_mass_fraction"].as_f64(),
                ..Default::default()
            };
            let h = bed["thickness_m"].as_f64().unwrap();
            let original = p["original_deposit_thickness_m"].as_f64().unwrap();
            let r =
                law::realize(&base, &delta, Some(law::Environment::Fluvial), original, h).unwrap();
            assert_eq!(
                json!(r.receipt.retained_thickness_m),
                p["retained_parting_thickness_m"]
            );
            assert!(r.receipt.max_mode_closure_error < 1e-12);
            assert!(r.receipt.relative_volume_closure_error < 1e-12);
            let packed = law::compact_optics(&r.host).unwrap();
            assert_eq!(
                packed["material_seed_words_le"],
                json!(
                    v3_2c::Evaluator::new(&r.host)
                        .unwrap()
                        .portable_structure()
                        .material_seed_words
                )
            );
            let host = v3_2c::Evaluator::new(&r.host)
                .unwrap()
                .surface_sample(&fresh, 2.)
                .unwrap();
            let parting = r.parting.as_ref().map(|b| {
                v3_2c::Evaluator::new(b)
                    .unwrap()
                    .surface_sample(&fresh, 2.)
                    .unwrap()
            });
            for sample in std::iter::once(&host).chain(parting.iter()) {
                for i in 0..3 {
                    assert!(
                        sample.diffuse_energy_linear[i] >= 0.
                            && sample.diffuse_energy_linear[i] + sample.f0_linear[i] <= 1. + 1e-12
                    );
                }
                assert!((0. ..=1.).contains(&sample.roughness));
            }
            let a =
                law::realize(&base, &delta, Some(law::Environment::Fluvial), original, h).unwrap();
            assert_eq!(
                serde_json::to_vec(&r.receipt).unwrap(),
                serde_json::to_vec(&a.receipt).unwrap()
            );
            assert_eq!(
                law::compact_optics(&r.host).unwrap(),
                law::compact_optics(&a.host).unwrap()
            );
            let full = r.sample_interval(&fresh, 2., 0., h, 0., h).unwrap();
            assert!((full.host_weight + full.parting_weight - 1.).abs() < 1e-12);
            let mut filters = Vec::new();
            for width in [0.001, 0.002, 0.005, 0.01, 0.1, 1., 10.] {
                let lo = h - r.receipt.retained_thickness_m * 0.5 - width * 0.5;
                let s = r
                    .sample_interval(&fresh, width, 0., h, lo, lo + width)
                    .unwrap();
                assert!(s.host_weight + s.parting_weight <= 1. + 1e-12);
                assert!(s.parting_weight <= r.receipt.retained_thickness_m / width + 1e-12);
                filters.push(json!({"vertical_pixel_width_m":width,"host_weight":s.host_weight,"parting_weight":s.parting_weight}));
            }
            rows.push(json!({"bed":bed,"properties":p,"band_base_m":band["offsets_m"][0],"red":band["redbed"],"receipt":r.receipt,"host":host,"parting":parting,"filters":filters}));
            covered += h;
        }
        assert_eq!(before, serde_json::to_vec(&base).unwrap());
    }
    assert_eq!(rows.len(), 51);
    println!("REAL_BEDS={} HEIGHT_M={covered}", rows.len());
    if let Ok(out) = std::env::var("VZ_BED_OPTICS_OUT") {
        std::fs::write(
            out,
            serde_json::to_vec_pretty(&json!({"book_pin":v["pin"],"rows":rows})).unwrap(),
        )
        .unwrap();
    }
}
#[test]
fn real_base_falsifiers_and_eolian_sensitivity() {
    let Some((v, b)) = source() else {
        return;
    };
    let id = b["stack"]["bands"][0]["id"].as_str().unwrap();
    let base: v3::ResolvedBody =
        serde_json::from_value(v["materials"][id]["resolved"].clone()).unwrap();
    let bed = &b["bed_contacts"]["bands"][0]["beds"][0];
    let h = bed["thickness_m"].as_f64().unwrap();
    for lat in [f64::NAN, f64::INFINITY, -91., 91.] {
        assert!(law::aridity_strength(lat).is_err());
    }
    assert_eq!(law::aridity_strength(30.).unwrap(), 1.);
    assert_eq!(law::aridity_strength(-30.).unwrap(), 1.);
    for w in [-0.1, 1.1, f64::NAN] {
        assert!(law::realize(
            &base,
            &law::Delta {
                hematite_mass_fraction: Some(w),
                ..Default::default()
            },
            None,
            h,
            h
        )
        .is_err());
    }
    let d = law::Delta::default();
    assert!(law::realize(&base, &d, None, h, -1.).is_err());
    assert!(law::realize(&base, &d, None, h, 2. * h).is_err());
    let eroded = law::realize(&base, &d, Some(law::Environment::Fluvial), h + 0.01, h).unwrap();
    assert_eq!(eroded.receipt.retained_thickness_m, 0.);
    assert!(eroded.parting.is_none());
    // Sensitivity uses a real sandstone base; explicitly not native eolian data.
    let eolian = law::realize(&base, &d, Some(law::Environment::DryEolian), h, h).unwrap();
    assert_eq!(eolian.receipt.whole_bed_mode, eolian.receipt.parting_mode);
    assert_eq!(eolian.receipt.retained_thickness_m, 0.003);
    let none = law::realize(&base, &d, None, h, h).unwrap();
    assert!(none.parting.is_none());
    assert!(law::overlap_m(2., 1., 0., 3.).is_err());
}
