//! Real retained C7 input only. Set VZ_E3_INPUT to the lane's input manifest.
use pf_region::t1b_litho::{EnvelopeConfig, ExposedLithology, LithologyV1};
use pf_region::t1b_rivernet::StrataCols;
use pf_region::t0::{SiteRecord,T0};
use std::path::Path;
fn inputs() -> Option<(EnvelopeConfig,T0,StrataCols)> {
    let Ok(p)=std::env::var("VZ_E3_INPUT") else {eprintln!("SKIP real E3: set VZ_E3_INPUT");return None};
    let j:serde_json::Value=serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap();
    let path=|k:&str|Path::new(j[k].as_str().unwrap());
    let site=SiteRecord::load(path("site")).unwrap();
    Some((EnvelopeConfig::from_path(j["envelope"].as_str().unwrap()).unwrap(),T0::load(path("world"),path("subgrid"),&site,0.).unwrap(),StrataCols::load(path("strata")).unwrap()))
}
#[test]
fn every_native_layer_including_thin_has_its_original_identity_and_rock() {
    let Some((cfg,t0,cols))=inputs() else{return};
    let env=cfg.load().unwrap();env.validate_inputs(&t0,&cols).unwrap();
    let cells:Vec<u32>=env.cells.keys().copied().collect();
    let pos:Vec<_>=cells.iter().map(|&c|{let(a,b)=t0.centre_deg(c);(a,b,0.,0.)}).collect();
    let table=LithologyV1::standin_v1();
    let rock=pf_region::t1b_lcont::envelope_rock(&cfg,&table,&t0,&cols,&pos).unwrap();
    let(mut thin,mut positive,mut zero,mut unresolved)=(0,0,0,0);
    for(i,&c)in cells.iter().enumerate(){let e=&env.cells[&c];for l in &e.layers{
        if l.h==0.{zero+=1;continue;}
        // All positive source thicknesses survive the adapter, including those
        // smaller than f64 absolute-height resolution. Do not claim such a layer
        // has a representable interior surface or discard its source budget.
        if let pf_region::t1b_rivernet::RockPosts::Strata {cont:Some(cp),..}=&rock {
            let b=cp.stencils[cp.st[i] as usize].bodies.iter().find(|b|b.src[0]==Some((c,l.index))).unwrap();
            assert_eq!(b.th[0],l.h); assert_eq!(b.code,l.code);
        }else{panic!("envelope adapter missing")}
        positive+=1;if l.h<1.{thin+=1;}
        let surface=l.base+l.h*0.5;
        if surface<=l.base || surface>=l.base+l.h {unresolved+=1;continue;}
        assert_eq!(rock.layer_at(i,surface),Some((c,l.index)),"{}",l.id);
        assert_eq!(rock.at(i,surface).2,l.code,"{}",l.id);
        assert_eq!(e.exposed(0.,surface).unwrap().index,l.index,"{} h={} base={} surface={} mean={}",l.id,l.h,l.base,surface,e.mean_m);
    }}
    println!("real C7 cells={} positive={} thin={} zero={} unrepresentable_interiors={}",cells.len(),positive,thin,zero,unresolved);
    assert_eq!(cells.len(),39363);assert_eq!(thin,20104);assert_eq!(zero,22);
}
#[test]
fn source_mismatch_and_digest_mutation_are_refused() {
    let Some((mut cfg,mut t0,cols))=inputs() else{return};
    let e=cfg.load().unwrap();let c=*e.cells.keys().next().unwrap()as usize;
    t0.elevation_m[c]+=1.;assert!(e.validate_inputs(&t0,&cols).is_err());
    cfg.digest="blake3:invalid".into();assert!(cfg.load().is_err());
}
#[test]
fn finite_occupancy_does_not_mint_top_extension_or_below_datum_stock() {
    let Some((cfg,_,_))=inputs()else{return};let env=cfg.load().unwrap();
    for e in env.cells.values(){
        for surface in [e.mean_m-e.relief_m,e.mean_m,e.mean_m+e.relief_m]{
            let v:f64=e.layers.iter().map(|l|e.occupancy(l,0.,surface)).sum();
            let extra=(surface-e.mean_m.max(e.sea_m)).max(0.);
            let solid=(surface-e.sea_m).max(0.);
            let bottom=e.layers.last().unwrap().base;
            let uncovered=(surface.min(bottom)-e.sea_m).max(0.);
            assert!((v+extra+uncovered-solid).abs()<1e-7);
            assert!(v<=e.layers.iter().map(|l|l.h).sum::<f64>()+1e-7);
        }
    }
}
#[test]
fn legacy_table_omits_the_absent_key_and_spherical_area_converges() {
    let table=serde_json::to_value(LithologyV1::standin_v1()).unwrap();assert!(table.get("envelope").is_none());
    // Coordinates derived from the actual C7 chart and radius, not a mock world.
    let Some((_,t0,_))=inputs()else{return};
    for x in [0.,16384.,131056.]{let(a,err)=pf_region::t1b_checks::envelope_area(x,0.,16.,16.,t0.radius_m);
        let mut refined=0.;for dx in [-8.,8.]{for dy in [-8.,8.]{refined+=pf_region::t1b_checks::envelope_area(x+dx,dy,8.,8.,t0.radius_m).0;}}
        assert!((a-refined).abs()<=err);assert!(a<=1024.&&a>1023.);
    }
}
#[test]
fn retained_legacy_keysets_reserialize_byte_equivalent_without_env() {
    let Ok(p)=std::env::var("VZ_E3_INPUT")else{eprintln!("SKIP legacy recipes: set VZ_E3_INPUT");return};
    let j:serde_json::Value=serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap();
    let mut count=0;
    let paths:Vec<std::path::PathBuf>=if let Some(v)=j["legacy_recipes"].as_array(){v.iter().map(|p|p.as_str().unwrap().into()).collect()}else{j["sites"].as_array().unwrap().iter().map(|s|Path::new(s["layer"].as_str().unwrap()).join("rivernet.json")).collect()};
    let mut unique=std::collections::BTreeSet::new();
    for path in paths {
        let h:serde_json::Value=serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        let params=&h["header"]["recipe"]["params"];
        let original=params.get("rivernet").unwrap_or(params);let mut v=original.clone();
        unique.insert(serde_json::to_string(original).unwrap());
        if v.get("f6_reference").is_none(){v["f6_reference"]=0.into();}
        let mut p:pf_region::t1b_rivernet::RiverNetParamsV1=serde_json::from_value(v).unwrap();
        assert!(p.lithology.envelope.is_none());
        assert_eq!(serde_json::to_vec(&serde_json::to_value(&p).unwrap()).unwrap(),serde_json::to_vec(original).unwrap());
        // env=0 has exactly the same representation on each retained recipe.
        p.lithology.envelope=None;
        assert_eq!(serde_json::to_value(&p).unwrap(),*original);
        count+=1;
    }
    println!("retained legacy recipes exactly preserved: {count}; distinct keysets: {}",unique.len());assert!(count>=3);
}
#[test]
fn real_constructed_lake_has_water_exposure_and_bad_rock_is_rejected() {
    let Ok(p)=std::env::var("VZ_E3_INPUT")else{eprintln!("SKIP real lake: set VZ_E3_INPUT");return};
    let j:serde_json::Value=serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap();
    let Some(w)=j["lake_window"].as_str()else{eprintln!("SKIP real lake: lake_window not supplied");return};
    let(raw,_)=pf_region::t1b_window::load(Path::new(w)).unwrap();
    let cfg=EnvelopeConfig::from_path(j["envelope"].as_str().unwrap()).unwrap();let env=cfg.load().unwrap();
    let(mut lakes,mut rocks)=(0,0);
    for i in 0..raw.lat.nodes(){if !raw.lat.in_core(i){continue;}let Some(e)=env.cells.get(&raw.t0_cell[i])else{continue};let code=raw.exposed_litho.as_ref().unwrap()[i];
        if raw.fixed[i]==pf_region::t1b::LAKE {
            assert!(pf_region::t1b_checks::envelope_exposure_matches(raw.fixed[i],Some(code),e.layers[0].code));
            assert!(!pf_region::t1b_checks::envelope_exposure_matches(raw.fixed[i],Some(e.layers[0].code),e.layers[0].code));lakes+=1;
        }else if code!=255 {
            assert!(!pf_region::t1b_checks::envelope_exposure_matches(raw.fixed[i],Some((code+1)%12),code));rocks+=1;
        }
    }
    println!("real lake regression: {lakes} water posts, {rocks} rock-post falsifiers");assert_eq!(lakes,378);assert!(rocks>0);
}
