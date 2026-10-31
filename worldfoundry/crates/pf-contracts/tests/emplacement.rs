use pf_contracts::{emplacement::*,params::ContinentalParcelCarryBlock};
#[test]
fn absent_key_has_identical_legacy_bytes_and_future_kind_refuses() {
    let c=ContinentalParcelCarryBlock::v1();
    assert_eq!(serde_json::to_string(&c).unwrap(),"{\"parcel_overlap_cell_fraction\":0.5,\"plate_fragment_min_area_fraction\":0.001}");
    assert!(serde_json::from_str::<BodyKind>("\"batholith_v99\"").is_err());
    for k in ["sheet","bed","sill","dyke","pluton","edifice","vein_set","ore_shell","unresolved"] {
        assert_eq!(serde_json::to_value(serde_json::from_value::<BodyKind>(serde_json::json!(k)).unwrap()).unwrap(),k);
    }
}
#[test]
fn finite_tet_and_prism_validate_degenerate_and_uncovered_refuse() {
    let mut g=FiniteGeometry{frame_id:"test_truth".into(),vertices_m:vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[1.,0.,1.],[0.,1.,1.]],elements:vec![Element::Tetrahedron([0,1,2,3]),Element::TriangularPrism([0,1,2,3,4,5])],coverage_min_m:[0.;3],coverage_max_m:[1.;3],coordinate_error_m:0.};
    g.validate().unwrap();g.elements=vec![Element::Tetrahedron([0,1,2,2])];assert!(g.validate().is_err());
    g.elements=vec![Element::Tetrahedron([0,1,2,3])];g.coverage_max_m[2]=0.5;assert!(g.validate().is_err());
    assert!(BodyInput{schema:SCHEMA.into(),manifest_sha256:"not-a-pin".into()}.validate().is_err());
}
