//! Reservation falsifiers for actual native mutation functions.
use pf_geodynamics::{carrier_history::StackHistory,conserve::CarrierIdentity,stratigraphy::*,layer_evolution};
use pf_contracts::params::CarrierHistoryBlock;
use pf_grid::GridSpec;
    fn params() -> StrataParams {
        StrataParams {
            max_layers: 64,
            merge_gap_my: 2.0,
            unconformity_gap_my: 5.0,
            min_layer_m: 10.0,
            trace_m: 0.01,
            arc_extrusive_fraction: 0.3,
            plume_extrusive_fraction: 0.2,
            shortening_pluton_fraction: 0.1,
            pluton_roof_m: 4_000.0,
            metamorphic_depth_m: 10_000.0,
            high_grade_depth_m: 20_000.0,
            shelf_depth_m: 200.0,
            carbonate_max_paleolat_deg: 30.0,
            carbonate_max_supply_m_per_my: 50.0,
            arid_lat_min_deg: 15.0,
            arid_lat_max_deg: 35.0,
            humid_lat_max_deg: 15.0,
            humid_midlat_min_deg: 40.0,
            humid_midlat_max_deg: 65.0,
            coal_max_elevation_m: 200.0,
            conglomerate_relief_m: 1000.0,
            river_mouth_area_m2: 1.0e11,
            publish_depth_m: 3000.0,
            basement_min_m: 5000.0,
            v2: None,
        }
    }


fn history(selected:Vec<u64>,bytes:u64)->StackHistory {
 StackHistory::new(CarrierHistoryBlock{selected,layer_evolution:true,event_limit:100000,byte_limit:bytes},params(),&GridSpec::new(288,576,6371000.).unwrap()).unwrap()
}
#[test]
fn maximum_insertions_and_erosion_respect_shape_and_change_bounds() {
 let p=params();let mut max_add=0;let mut max_merge=0;let mut max_changes=0;
 for n in 0..64 {
  for erosion in [0.,1.,9.5,10.,n as f64*10.+1.] {
   for total in [0.,100.,100000.] {
    let mut st=StrataStack::basement_only(50000.,0,-2000.);
    st.cover=(0..n).map(|i|StrataLayer{thickness_m:10.,age_base_my:-500.,age_top_my:-400.,lith:if i%2==0{5}else{6},env:6,..Default::default()}).collect();
    let input=st.clone();let mut log=StrataStepLog::default();let mut plain_log=log.clone();let mut plain=st.clone();
    let a=ParcelAmounts{erosion,deposit:10.,arc:100.,plume:100.,thinning:1.,basement_add:1.,short_pos:1.,short_neg:1.,land:true,class:5,belt_rest:-1.,warm_shelf:false};
    apply_parcel(&mut plain,&a,total,1.,&p,0.,2.,&mut plain_log);
    let trace=apply_parcel_recorded(&mut st,&a,total,1.,&p,0.,2.,&mut log,layer_evolution::byte_bound(64)).unwrap();
    assert_eq!(plain,st);assert_eq!(plain_log,log);
    assert!(trace.changes.len()<=n+32);
    for c in &trace.changes {assert!(c.before.cover.len()<=n+4 && c.after.cover.len()<=n+4);max_add=max_add.max(c.after.cover.len().saturating_sub(n));}
    max_merge=max_merge.max(trace.changes.iter().filter(|c|c.kind==17).count());max_changes=max_changes.max(trace.changes.len());
    let h=history(vec![0],1<<30);let ids=CarrierIdentity::initial(288*576,0.,&[(0,0)],1000000).unwrap();
    assert!(h.live_step_budget(&ids,&[input]).unwrap().1<(1<<30));
   }
  }
 }
 assert_eq!(max_add,4);assert_eq!(max_merge,4);eprintln!("960 native operation cases; max changes={max_changes}; transient additions={max_add}; capacity merges={max_merge}");
}
#[test]
fn future_birth_and_maximum_shapes_reserve_and_oversize_refuses() {
 let ids=CarrierIdentity::initial(288*576,0.,&(0..70).map(|i|(i,0)).collect::<Vec<_>>(),1000000).unwrap();
 let mut s=StrataStack::basement_only(1000.,0,0.);s.cover=vec![StrataLayer{thickness_m:1.,..Default::default()};63];
 let mut h=history((0..70).collect(),1<<30);h.reserve_live_step(&ids,&vec![s.clone();70]).unwrap();
 let (count,bytes)=h.live_step_budget(&ids,&vec![s.clone();70]).unwrap();assert_eq!(count,71);assert_eq!(bytes,599695360);
 let before=h.clone();assert!(h.reserve_live_step(&ids,&[]).is_err());assert_eq!(h,before);
 s.cover.push(StrataLayer::default());assert!(h.reserve_live_step(&ids,&vec![s;70]).is_err());assert_eq!(h,before);
 let mut births=history(vec![99999],1<<30);births.reserve_live_step(&ids,&vec![StrataStack::basement_only(1000.,0,0.);70]).unwrap();
 assert!(StackHistory::validate_config(&CarrierHistoryBlock{selected:(0..257).collect(),layer_evolution:true,event_limit:100000,byte_limit:1<<30}).is_err());
}
#[test]
fn actual_retained_selected_stack_reservation() {
 let Ok(path)=std::env::var("VZ_BODY_SELECTED_STACKS") else {eprintln!("SKIP retained shape gate: set VZ_BODY_SELECTED_STACKS");return};
 let v:serde_json::Value=serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
 let selected:Vec<u64>=serde_json::from_value(v["selected"].clone()).unwrap();
 let mut stacks=vec![StrataStack::basement_only(1.,0,0.);selected.last().copied().unwrap() as usize+1];
 for &id in &selected {
  let r=&v["strata"][id.to_string()];let b=r["b"].as_array().unwrap();
  let mut s=StrataStack::basement_only(b[0].as_f64().unwrap(),b[1].as_u64().unwrap() as u8,b[2].as_f64().unwrap());
  s.last_exposure_my=b[3].as_f64().unwrap();s.rock_uplift_m=r["u"][0].as_f64().unwrap();s.elev_ref_m=r["u"][1].as_f64().unwrap();
  for l in r["cover"].as_array().unwrap() {s.cover.push(StrataLayer{thickness_m:l[0].as_f64().unwrap(),age_base_my:l[1].as_f64().unwrap(),age_top_my:l[2].as_f64().unwrap(),max_burial_m:l[3].as_f64().unwrap(),lith:l[4].as_u64().unwrap() as u8,env:l[5].as_u64().unwrap() as u8,grade:l[6].as_u64().unwrap() as u8,flags:l[7].as_u64().unwrap() as u8});}
  stacks[id as usize]=s;
 }
 let ids=CarrierIdentity::initial(288*576,0.,&(0..stacks.len()).map(|i|(i,0)).collect::<Vec<_>>(),1000000).unwrap();
 let mut h=history(selected,1<<30);let (_,bytes)=h.live_step_budget(&ids,&stacks).unwrap();h.reserve_live_step(&ids,&stacks).unwrap();
 assert_eq!(bytes,155103232);eprintln!("retained 70 selected shapes: {bytes} reserved bytes; synthetic identity carrier positions are irrelevant to the polar bound");
}
#[test]
fn captured_native_final_coverage_does_not_imply_historical_coverage() {
 let Ok(path)=std::env::var("VZ_BODY_CAPTURE") else {eprintln!("SKIP real capture: set VZ_BODY_CAPTURE");return};
 let h=StackHistory::from_json(serde_json::from_reader(std::io::BufReader::new(std::fs::File::open(path).unwrap())).unwrap()).unwrap();
 let a=pf_geodynamics::body_history::audit(&h,200,2.).unwrap();
 let ops=a["operations"].as_array().unwrap();
 for (id,time,cell,fraction) in [(44000,332.,145322,0.1379915612661005),(34123,166.,141352,0.3874213253430756)] {
  let past=ops.iter().find(|o|o["carrier"]==id && o["end_my"]==time).unwrap();
  let gap=past["cell_coverage_gaps"].as_array().unwrap().iter().find(|g|g["cell"]==cell).unwrap();
  assert!((gap["unrecorded_fraction"].as_f64().unwrap()-fraction).abs()<1e-12);
  let final_row=ops.iter().find(|o|o["carrier"]==id && o["end_my"]==400.).unwrap();
  assert!(final_row["cell_coverage_gaps"].as_array().unwrap().is_empty());
 }
}
#[test]
fn actual_source_accounts_close_without_changing_retained_provenance() {
    let Ok(path)=std::env::var("VZ_BODY_CAPTURE") else {eprintln!("SKIP source accounting: set VZ_BODY_CAPTURE");return};
    let h=StackHistory::from_json(serde_json::from_reader(std::io::BufReader::new(std::fs::File::open(path).unwrap())).unwrap()).unwrap();
    let targets=vec![(44000,"carrier:44000:plume-volcanic:4074a00000000000".into()),(25323,"carrier:25323:initial-cover-0:c082c00000000000".into())];
    let books=pf_geodynamics::body_history::source_books(&h,&targets).unwrap();
    for (b,expected) in books.iter().zip([76.34752871064003,226.8946150451091]) {
        assert!((b.report["retained_m"].as_f64().unwrap()-expected).abs()<1e-10);
        assert!(b.report["mixed_target_boundaries"].as_array().unwrap().is_empty());
        assert_eq!(b.report["final_top_face_creation_interval_my"],serde_json::json!([398.,400.]));
        assert!(b.report["balance_residual_m"].as_f64().unwrap().abs()<=b.report["roundoff_bound_m"].as_f64().unwrap());
    }
    assert_eq!(books[0].report["emplacement_phases"].as_array().unwrap().len(),2);
    assert!(!books[0].report["within_episode_temporal_differences"].as_array().unwrap().is_empty());
}
