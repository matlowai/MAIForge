//! Reuse the frozen existing comparison set for the explicit key-off gate.
use pf_rockspec::{v3_1, v3_2c};
#[test]
fn all_existing_rows_byte_identical() {
    let Ok(path) = std::env::var("VZ_BED_LEGACY_BASELINE") else {
        eprintln!("SKIP VZ_BED_LEGACY_BASELINE unset");
        return;
    };
    assert!(std::fs::metadata(&path).unwrap().len() < 200_000_000);
    let before = std::fs::read(path).unwrap();
    let c = v3_1::catalogue().unwrap();
    let mut rows = vec![];
    for f in v3_1::FixtureSet::comparisons(&c).unwrap().samples {
        let b = f.resolve_fixture(&c).unwrap();
        let e = v3_2c::Evaluator::new(&b).unwrap();
        rows.push(serde_json::json!({"id":f.id,"history":f.history,"sample":e.surface_sample(&f.surface,f.footprint_m).unwrap(),"portable":e.portable_structure()}));
    }
    let bytes = serde_json::to_vec_pretty(&rows).unwrap();
    if let Ok(out) = std::env::var("VZ_BED_LEGACY_OUT") {
        std::fs::write(out, &bytes).unwrap();
    }
    assert!(
        before == bytes,
        "legacy sample/portable-structure bytes differ"
    );
    assert_eq!(rows.len(), 603);
    println!("LEGACY_ROWS={}; BYTE_IDENTICAL=true", rows.len());
}
