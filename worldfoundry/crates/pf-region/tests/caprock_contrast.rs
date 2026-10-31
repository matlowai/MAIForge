//! Real E3 Moab surface and envelope. The uniform-lithology counterfactual is
//! the explicitly required falsifier, never a substitute for the native run.
use pf_region::t1b_ct::{cap_contact,contrast_kernel,ContrastParams};
use pf_region::t1b_litho::{EnvelopeConfig,LithologyV1};
use pf_region::t1b_rivernet::{StrataCols,RockPosts};
use pf_region::t0::{SiteRecord,T0};
use std::path::Path;

#[test]
fn real_moab_contact_process_and_uniform_lithology_falsifier() {
    let Ok(path)=std::env::var("VZ_E5_INPUT")else{eprintln!("SKIP real E5: set VZ_E5_INPUT to input manifest");return};
    let j:serde_json::Value=serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    let get=|k:&str|Path::new(j[k].as_str().unwrap());
    let site=SiteRecord::load(get("site")).unwrap();let t0=T0::load(get("world"),get("subgrid"),&site,0.).unwrap();
    let cols=StrataCols::load(get("strata")).unwrap();let cfg=EnvelopeConfig::from_path(j["envelope"].as_str().unwrap()).unwrap();
    let mut table=LithologyV1::standin_v1();table.envelope=Some(cfg.clone());
    let w=get("window");let wh:serde_json::Value=serde_json::from_slice(&std::fs::read(w.join("window.json")).unwrap()).unwrap();
    let h=&wh["header"];let lat=&h["lattice"];let n=lat["n"].as_u64().unwrap()as usize;
    let gx0=lat["gx0"].as_i64().unwrap();let gy0=lat["gy0"].as_i64().unwrap();
    let chart=&h["recipe"]["chart"];let r=chart["datum_radius_mm"].as_f64().unwrap()/1000.;
    let la0=(chart["origin_lat_e7"].as_f64().unwrap()/1e7).to_radians();let lo0=(chart["origin_lon_e7"].as_f64().unwrap()/1e7).to_radians();
    let bytes=std::fs::read(w.join("fields/height_mm.i32")).unwrap();
    let initial:Vec<i64>=bytes.chunks_exact(4).map(|b|i32::from_le_bytes(b.try_into().unwrap())as i64).collect();assert_eq!(initial.len(),n*n);
    let pos:Vec<_>=(0..n*n).map(|i|{let x=(gx0+(i%n)as i64)as f64*32.;let y=(gy0+(i/n)as i64)as f64*32.;
        let a=x/(2.*r);let b=y/(2.*r);let q=a*a+b*b;let ve=2.*a/(1.+q);let vn=2.*b/(1.+q);let vu=(1.-q)/(1.+q);
        let px=la0.cos()*vu-la0.sin()*vn;let pz=la0.sin()*vu+la0.cos()*vn;
        (pz.atan2(px.hypot(ve)).to_degrees(),(lo0+ve.atan2(px)).to_degrees(),x,y)}).collect();
    let mut rock=pf_region::t1b_lcont::envelope_rock(&cfg,&table,&t0,&cols,&pos).unwrap();let p=ContrastParams::v1();
    let contacts:Vec<_>=(0..n*n).map(|i|cap_contact(&rock,i,initial[i]as f64/1000.,&p.weathering_table)).collect();
    let count=contacts.iter().filter(|c|c.is_some()).count();assert!(count>0,"real contact coverage");
    // Kernel test uses real heights; pin only outer ring here. Full pipeline
    // additionally pins actual channels/floors/banks/overwritten terrain.
    let pinned:Vec<_>=(0..n*n).map(|i|i%n==0||i/n==0||i%n+1==n||i/n+1==n).collect();
    let mut z=initial.clone();let stats=contrast_kernel(n,32000,&pinned,&contacts,&p,&mut z).unwrap();
    assert_eq!(stats.local_ledger_residual_mm,0);assert!(z.iter().zip(&initial).all(|(a,b)|a<=b));assert!(z!=initial,"real contrast changes surface");
    for i in 0..n*n{if pinned[i]{assert_eq!(z[i],initial[i]);}}
    // A real cap geometry with equal effective rates is NOT a contrast.
    let mut equal = p.weathering_table.clone();
    equal.classes[6].bare_weathering_mm_per_kyr=equal.classes[5].bare_weathering_mm_per_kyr;
    assert!((0..n*n).all(|i|cap_contact(&rock,i,initial[i]as f64/1000.,&equal).is_none()));
    // Uniform sandstone control. All bodies retain actual source geometry,
    // positions and thickness but have the SAME manufactured lithology.
    for code in [6,5] {
        let RockPosts::Strata{cont:Some(cp),..}=&mut rock else{panic!("expected native envelope")};
        for s in &mut cp.stencils {for b in &mut s.bodies{b.code=code;}}
        let uniform:Vec<_>=(0..n*n).map(|i|cap_contact(&rock,i,initial[i]as f64/1000.,&p.weathering_table)).collect();
        assert!(uniform.iter().all(|c|c.is_none()));let mut u=initial.clone();let control=contrast_kernel(n,32000,&pinned,&uniform,&p,&mut u).unwrap();assert_eq!(u,initial);assert_eq!(control.failed_slabs,0);assert_eq!(control.cap_edge_posts,0);
    }
    println!("E5 real Moab contacts={count}; slabs={}; vertical_edges={}; uniform_posts={} changed=0; ledger_mm=0",stats.failed_slabs,stats.vertical_cap_edge_posts,n*n);
}
#[test]
fn defaults_keep_legacy_table_and_hillslope_recipe_exact() {
    let t=serde_json::to_value(LithologyV1::standin_v1()).unwrap();
    for c in t["classes"].as_array().unwrap(){assert!(c.get("bare_weathering_mm_per_kyr").is_none());}
    let h=serde_json::to_value(pf_region::t1b_hs::HsParamsV1::v0()).unwrap();assert!(h.get("contrast").is_none());
    let mut p=ContrastParams::v1();p.steps=0;assert!(p.validate().is_err());p.steps=257;assert!(p.validate().is_err());
}
