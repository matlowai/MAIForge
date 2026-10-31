use pf_region::{t1b_litho::EnvelopeConfig,t1b_ve::material::{self,CompactInput}};
use std::path::PathBuf;
#[test]
fn actual_compact_source_ownership_and_falsifiers(){
    let (Some(tile),Some(env),Some(out))=(std::env::var_os("VZ_ARC4_MATERIAL"),std::env::var_os("VZ_ARC4_ENVELOPE"),std::env::var_os("VZ_ARC4_REPORT")) else {eprintln!("SKIP: supply VZ_ARC4_MATERIAL, VZ_ARC4_ENVELOPE, VZ_ARC4_REPORT");return};
    let tile=PathBuf::from(tile);let out=PathBuf::from(out);assert!(!out.exists(),"preserve evidence");
    let cfg=EnvelopeConfig::from_path(env.to_str().unwrap()).unwrap();
    let input=CompactInput::load(&tile).unwrap();
    let legacy=std::env::var_os("VZ_ARC4_AUDIT_LEGACY").is_some();
    let owned=material::bind(input,&cfg,legacy).unwrap();
    assert!(owned.report.proxy_positive_posts>0);
    assert!(!owned.report.units.is_empty());
    assert!(owned.report.solid_closure_residual_m3<=owned.report.quadrature_and_roundoff_bound_m3);
    assert!(owned.report.units.iter().all(|u|u.unallocated_parent_m3>=0.0));
    assert!(owned.report.source_mobile_m3<=owned.report.requested_proxy_m3);
    for i in 0..owned.support_m2.len(){if owned.excess_depth_m[i]>0.0{assert!(owned.mobile_columns[i].is_empty());}}
    let status=material::require_evolution_stock(&owned);
    let report=owned.report.clone();
    let mut bad=CompactInput::load(&tile).unwrap();bad.header["source_envelope"]["digest"]="wrong-envelope".into();
    assert!(material::bind(bad,&cfg,false).is_err());
    let twin=material::bind(CompactInput::load(&tile).unwrap(),&cfg,legacy).unwrap();
    let one=material::partition(std::slice::from_ref(&owned)).unwrap();
    assert_eq!(one.finite_occupancy_m3,report.finite_occupancy_m3);
    let duplicate=material::partition(&[owned,twin]).unwrap_err();assert!(duplicate.contains("duplicate source"));
    std::fs::write(out,serde_json::to_vec_pretty(&serde_json::json!({"input":tile,"ownership":report,"process_admission":status,"falsifiers":{"wrong_envelope_rejected":true,"duplicate_owner_rejected":true}})).unwrap()).unwrap();
}
