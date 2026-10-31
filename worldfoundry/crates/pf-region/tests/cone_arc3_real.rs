use pf_region::t1b_surface_completion::*;
use pf_region::t1b_surface_domain::PostId;
use std::path::{Path,PathBuf};
fn load(root:&Path,name:&str)->(Rect,Vec<f64>,Vec<u8>,String){
    let dir=std::fs::read_dir(root.join(name)).unwrap().map(|r|r.unwrap().path()).find(|p|p.file_name().unwrap().to_str().unwrap().starts_with("coast_final_")).unwrap();
    let meta:serde_json::Value=serde_json::from_slice(&std::fs::read(dir.join("capture.json")).unwrap()).unwrap();
    let bytes=std::fs::read(dir.join("ground.f64")).unwrap();assert_eq!(blake3::hash(&bytes).to_hex().as_str(),meta["fields"]["ground.f64"].as_str().unwrap());
    let fixed=std::fs::read(dir.join("fixed.u8")).unwrap();assert_eq!(blake3::hash(&fixed).to_hex().as_str(),meta["fields"]["fixed.u8"].as_str().unwrap());
    let rect=Rect{gx0:meta["gx0"].as_i64().unwrap(),gy0:meta["gy0"].as_i64().unwrap(),width:meta["n"].as_u64().unwrap() as usize,height:meta["n"].as_u64().unwrap() as usize};
    let z=bytes.chunks_exact(8).map(|c|f64::from_le_bytes(c.try_into().unwrap())).collect();
    (rect,z,fixed,format!("captured-stage:{}",meta["fields"]["ground.f64"].as_str().unwrap()))
}
fn root()->Option<PathBuf>{let p=std::env::var_os("VZ_SURFACE_WITNESS");if p.is_none(){eprintln!("SKIP C7 selective witness: set VZ_SURFACE_WITNESS to native capture root");}p.map(PathBuf::from)}
fn packet(pin:&str,samples:Vec<Sample>)->Fetch{Fetch{predecessor:pin.into(),content_digest:pf_region::t1b_surface_domain::content_digest(&samples).unwrap(),samples}}
#[test]
fn selective_real_fill_completes_witness_without_ring_terminal(){
    let Some(root)=root()else{return};let(rect,z,fixed,pin)=load(&root,"capture_1_5");
    let targets:Vec<usize>=(2496..3137).flat_map(|gy|(960..1089).map(move|gx|PostId{gy,gx})).map(|id|rect.index(id).unwrap()).collect();assert_eq!(targets.len(),82689);
    // One pinned native predecessor dataset, initially cropped at the historical
    // bad window's west edge. This is a source-bound kernel test, not an assertion
    // that the diagnostic capture has certified native upstream stages.
    let mut working=WorkingSet::new(rect,pin.clone(),2_000_000).unwrap();
    let initial=(0..z.len()).filter(|&i|rect.id(i).gx>=704).map(|i|Sample{id:rect.id(i),ground_m:z[i],fixed:fixed[i]}).collect();working.ingest(packet(&pin,initial)).unwrap();
    let before=pf_region::t1b_ls::condition_detail_bounds(&working.patch()).unwrap();let cert0=before.certify(&working.patch(),&targets,1024.0).unwrap();assert!(!cert0.pass);let wi=rect.index(PostId{gx:960,gy:2520}).unwrap();assert!(!before.height_exact(wi));
    // The actual chart contract must not close a cropped window inside it.
    working.native_chart=Some(pin.clone());
    assert!(matches!(before.certify(&working.patch(),&targets,1024.0),Err(pf_region::t1b_surface_domain::DomainError::Stale{..})));
    let crop_check=pf_region::t1b_ls::condition_detail_bounds(&working.patch()).unwrap().certify(&working.patch(),&targets,1024.0).unwrap();
    assert_eq!(crop_check.unresolved_posts,cert0.unresolved_posts);assert_eq!(crop_check.missing,cert0.missing);
    working.native_chart=None;
    let mut requests=Vec::new();
    let mut fetch=|ids:&[PostId]|{
        let tiles:std::collections::BTreeSet<_>=ids.iter().map(|id|(id.gy.div_euclid(128),id.gx.div_euclid(128))).collect();requests.push(tiles.iter().copied().collect::<Vec<_>>());
        for id in ids{if rect.index(*id).is_none(){return Err(pf_region::t1b_surface_domain::DomainError::NeedsDomain{stage:"outside retained exact predecessor".into(),frontier:ids.to_vec()});}}
        let samples=(0..z.len()).filter(|&i|{let id=rect.id(i);tiles.contains(&(id.gy.div_euclid(128),id.gx.div_euclid(128)))}).map(|i|Sample{id:rect.id(i),ground_m:z[i],fixed:fixed[i]}).collect();Ok(packet(&pin,samples))
    };
    let(after,report)=working.complete(&targets,1024.0,32,&mut fetch).unwrap();
    let p=Patch{rect,ground:&z,fixed:&fixed,present:&[],predecessor:&pin,native_chart:None};let oracle=pf_region::t1b_ls::condition_detail_bounds(&p).unwrap();let full=oracle.certify(&p,&targets,1024.0).unwrap();assert!(full.pass);
    for &i in &targets{assert_eq!(after.upper[i].to_bits(),oracle.upper[i].to_bits());}
    assert_eq!(after.upper[wi].to_bits(),z[wi].to_bits());let path=after.path(&working.patch(),wi).unwrap();assert!(path.iter().any(|id|id.gx<704));
    let end=rect.index(*path.last().unwrap()).unwrap();assert_ne!(fixed[end],pf_region::t1b::FREE);assert_ne!(fixed[end],pf_region::t1b::RING);
    let native=root.join("probe_1_5/w_1_5/fields/height_mm.i32");let native=std::fs::read(native).unwrap();let native:Vec<i32>=native.chunks_exact(4).map(|b|i32::from_le_bytes(b.try_into().unwrap())).collect();
    let mut mm_bad=0;let mut max_mm=0i64;
    for &i in &targets{let id=rect.id(i);let j=((id.gy-2496)*641+(id.gx-448)) as usize;let actual=vz_tilekey::chart::quantize_mm(after.upper[i]*1000.0);let delta=(actual-native[j] as i64).abs();if delta!=0{mm_bad+=1;max_mm=max_mm.max(delta);}}
    let record=serde_json::json!({"scope":"conditional native pre-coast-fill input test; no full world acceptance","before":cert0,"completion":report,"oracle":full,"witness_ground_m":z[wi],"witness_after_m":after.upper[wi],"known_terminal":path.last(),"path_posts":path.len(),"path":path,"requests_128post_tiles":requests,"published_reference_mm_mismatches":mm_bad,"published_reference_max_mm":max_mm});
    if let Some(dir)=std::env::var_os("VZ_SURFACE_WITNESS_OUTPUT"){std::fs::create_dir_all(&dir).unwrap();std::fs::write(PathBuf::from(dir).join("selective-witness.json"),serde_json::to_vec_pretty(&record).unwrap()).unwrap();}
    eprintln!("selective witness: {} exact posts, {} fetch iterations, {} fetched posts; published reference {} differences max {} mm",targets.len(),report.iterations,report.fetched_posts,mm_bad,max_mm);
    assert_eq!(mm_bad,0);
}

#[test]
fn actual_chart_boundary_is_not_an_arbitrary_crop_outlet(){
    let Some(root)=std::env::var_os("VZ_SURFACE_NATIVE_CACHE").map(PathBuf::from)else{eprintln!("SKIP actual chart boundary: set VZ_SURFACE_NATIVE_CACHE to the failed native run's exact source-pin directory");return};
    let rect=Rect{gx0:0,gy0:2560,width:1536,height:2048};let pin=root.file_name().unwrap().to_str().unwrap().to_string();let mut working=WorkingSet::new(rect,pin.clone(),4_000_000).unwrap();
    let half=vz_tilekey::CHART_HALF_EXTENT_MM/pf_region::t1b::T1B_SPACING_MM;let mut physical_boundaries=0;
    for ty in 5..=8{for tx in 0..=2{
        let stem=format!("main_{tx}_{ty}");let meta:serde_json::Value=serde_json::from_slice(&std::fs::read(root.join(format!("{stem}.json"))).unwrap()).unwrap();let bytes=std::fs::read(root.join(format!("{stem}.bin"))).unwrap();
        assert_eq!(meta["source_pin"],pin);assert_eq!(blake3::hash(&bytes).to_hex().as_str(),meta["digest"].as_str().unwrap());assert_eq!(bytes.len(),512*512*13);
        for y in 0..512{for x in 0..512{let id=PostId{gx:tx*512+x as i64,gy:ty*512+y as i64};if id.gy>half{continue;}let i=rect.index(id).unwrap();let k=y*512+x;working.ground[i]=f64::from_le_bytes(bytes[k*8..k*8+8].try_into().unwrap());working.fixed[i]=bytes[512*512*8+k];working.present[i]=true;if working.fixed[i]!=0{physical_boundaries+=1;}}}
    }}
    let targets:Vec<usize>=(2944..=4096).flat_map(|gy|(192..=1344).map(move|gx|PostId{gy,gx})).map(|id|rect.index(id).unwrap()).collect();
    let open=pf_region::t1b_ls::condition_detail_bounds(&working.patch()).unwrap();let before=open.certify(&working.patch(),&targets,1024.0).unwrap();assert!(!before.pass);assert!(before.missing.iter().all(|id|id.gy>half||id.gy< -half||id.gx>half||id.gx< -half));
    working.native_chart=Some(pin.clone());
    assert!(matches!(open.certify(&working.patch(),&targets,1024.0),Err(pf_region::t1b_surface_domain::DomainError::Stale{..})));
    let closed=pf_region::t1b_ls::condition_detail_bounds(&working.patch()).unwrap();let after=closed.certify(&working.patch(),&targets,1024.0).unwrap();assert!(after.pass);
    for &i in &targets{assert_eq!(open.upper[i].to_bits(),closed.upper[i].to_bits());}
    assert_eq!(working.fixed.iter().enumerate().filter(|(i,f)|working.present[*i]&&**f!=0).count(),physical_boundaries);
    let record=serde_json::json!({"source_pin":pin,"old_unknown_exterior":before,"native_chart_contract":after,"physical_boundaries_unchanged":physical_boundaries,"new_terminals":0,"upper_height_bit_changes":0,"contract":"existing native LS finite chart D8 graph, no exterior land edge; native T1a/water seeds unchanged"});
    if let Some(dir)=std::env::var_os("VZ_SURFACE_WITNESS_OUTPUT"){std::fs::create_dir_all(&dir).unwrap();std::fs::write(PathBuf::from(dir).join("native-chart-boundary.json"),serde_json::to_vec_pretty(&record).unwrap()).unwrap();}
    eprintln!("native chart bound: {} old unresolved -> {}, no height or terminal change",before.unresolved_posts,after.unresolved_posts);
}
