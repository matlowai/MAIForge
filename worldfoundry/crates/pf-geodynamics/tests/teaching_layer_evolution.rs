use pf_geodynamics::stratigraphy::*;
use pf_geodynamics::layer_evolution::{byte_bound,change_bound,BASEMENT};
fn params(capacity: usize) -> StrataParams {
        StrataParams {
            max_layers: capacity,
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
            v2: Some(StrataV2Params {
                shelf_width_m: 80000.0, slope_foot_depth_m: 4000.0,
                shelf_majority_fraction: 0.5, platform_cover_max_m: 2000.0,
                platform_full_elev_m: 0.0, platform_zero_elev_m: 600.0,
                platform_sandstone_fraction: 0.3, platform_carbonate_fraction: 0.3,
                platform_age_base_my: -600.0, platform_age_top_my: -50.0,
            }),
        }
    }

fn layer(thickness:f64,lith:u8,age:f64)->StrataLayer{StrataLayer{thickness_m:thickness,lith,env:ENV_FLUVIAL,age_base_my:age,age_top_my:age+1.,..Default::default()}}
fn check(st:StrataStack,a:ParcelAmounts,total:f64,p:&StrataParams)->pf_geodynamics::layer_evolution::LayerEvolution{
 let mut plain=st.clone();let mut actual=st;let mut plain_log=StrataStepLog::default();let mut log=plain_log.clone();
 apply_parcel(&mut plain,&a,total,1e6,p,10.,12.,&mut plain_log);
 let trace=apply_parcel_recorded(&mut actual,&a,total,1e6,p,10.,12.,&mut log,byte_bound(p.max_layers)).unwrap();
 assert_eq!(actual,plain);assert_eq!(log,plain_log);assert!(trace.changes.len()<=change_bound(p.max_layers));
 for c in &trace.changes {
  assert!(c.links.iter().all(|[i,j]|(*i==BASEMENT||(*i as usize)<c.before.cover.len())&&(*j==BASEMENT||(*j as usize)<c.after.cover.len())));
  if c.kind==17||c.kind==3 {
   for target in 0..c.after.cover.len(){let parents:Vec<_>=c.links.iter().filter(|[i,j]|*i!=BASEMENT&&*j==target as u32).map(|[i,_]|*i as usize).collect();if parents.len()>1{let sum:f64=parents.iter().map(|i|c.before.cover[*i].thickness_m).sum();assert_eq!(sum,c.after.cover[target].thickness_m,"native merge links must name the actual source layers");}}
   assert!((c.before.total_m()-c.after.total_m()).abs()<1e-7,"merge/absorption must preserve local column thickness");
  }
  if c.kind==7{assert!((c.before.total_m()-c.after.total_m()).abs()<1e-7,"roof is transferred from basement, not created");}
 }

 trace
}
#[test]
fn actual_capacity_merge_and_bottom_up_deficit_are_recorded(){
 for cap in [16,64]{for v2 in [false,true]{let mut p=params(cap);if !v2{p.v2=None;}
 let mut st=StrataStack::basement_only(10000.,LITH_BASEMENT,0.);st.cover=(0..cap-1).map(|i|layer(100.+i as f64,(i%4+3) as u8,-100.+i as f64*3.)).collect();
 let mut a=ParcelAmounts::default();a.deposit=200.;a.class=3;
 let trace=check(st.clone(),a,st.total_m()+200.,&p);assert!(trace.changes.iter().any(|c|c.kind==17));
 let trace=check(st.clone(),ParcelAmounts::default(),10.,&p);assert!(trace.changes.iter().any(|c|c.kind==15));assert!(trace.changes.iter().any(|c|c.kind==16));
 }}
}
#[test]
fn thin_absorption_roof_birth_and_budget_refusal_are_explicit(){
 let p=params(16);let mut st=StrataStack::basement_only(20000.,LITH_BASEMENT,-100.);st.cover=vec![layer(100.,LITH_SHALE,-20.),layer(11.,LITH_SANDSTONE,-10.)];
 let mut a=ParcelAmounts::default();a.erosion=2.;let trace=check(st.clone(),a,st.total_m()-2.,&p);assert!(trace.changes.iter().any(|c|c.kind==3));
 let mut a=ParcelAmounts::default();a.arc=1000.;let trace=check(st.clone(),a,st.total_m()+1000.,&p);assert!(trace.changes.iter().any(|c|c.kind==7&&c.links.iter().any(|[i,j]|*i==BASEMENT&&*j!=BASEMENT)));
 let before=st.clone();let mut log=StrataStepLog::default();let old_log=log.clone();assert!(apply_parcel_recorded(&mut st,&ParcelAmounts::default(),20000.,1e6,&p,0.,2.,&mut log,byte_bound(p.max_layers)-1).is_err());assert_eq!(st,before);assert_eq!(log,old_log);
 let mut plain_log=StrataStepLog::default();let mut log=plain_log.clone();let plain=birth_stack(30000.,100.,2000.,3,&p,20.,1e6,&mut plain_log);let (actual,trace)=birth_stack_recorded(30000.,100.,2000.,3,&p,20.,1e6,&mut log,byte_bound(p.max_layers)).unwrap();assert_eq!(plain,actual);assert_eq!(log,plain_log);assert!(!trace.changes.is_empty());
}
