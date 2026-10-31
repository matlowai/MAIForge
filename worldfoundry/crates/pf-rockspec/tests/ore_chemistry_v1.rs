#[test]
fn reference_data_loads_without_mutating_legacy_catalogue() {
    let data: serde_json::Value = pf_rockspec::ore_chemistry_v1::load().unwrap();
    assert_eq!(data["schema"], "pf.ore-chemistry-table/v1");
    assert_eq!(data["endmembers"].as_array().unwrap().len(), 18);
    assert_eq!(data["groups"].as_array().unwrap().len(), 7);
    assert!(data["groups"]
        .as_array()
        .unwrap()
        .iter()
        .any(|g| g["id"] == "organic_matter"));
}
