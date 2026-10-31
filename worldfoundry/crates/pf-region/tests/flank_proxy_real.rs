//! Native proxy binding audit on actual registered construction output.
//! This does not run or claim positive-stock evolution or geology qualification.
use std::path::PathBuf;
#[test]
fn bind_real_native_proxy_and_preserve_boundary_stock_scope(){
    let Some(path)=std::env::var_os("VZ_NATIVE_PROXY_WINDOW").map(PathBuf::from)else{eprintln!("SKIP native proxy binding: set VZ_NATIVE_PROXY_WINDOW");return};
    let input=pf_region::t1b_ve::flank3_load_native_window(&path).unwrap();
    let lat=input.lat;let n=lat.n;let z=input.ground_m;let fixed=input.geometry_fixed;
    let proxy=input.proxy;
    assert!(proxy.report["positive_mobile_posts"].as_u64().unwrap()>0);assert!(!proxy.source_names.is_empty());
    assert!(proxy.stock_depth_m.iter().all(|x|x.is_finite()&&*x>=0.0));assert!(proxy.material_m.iter().all(|x|x.is_finite()&&*x>=0.0));
    for i in 0..n*n{if proxy.stock_depth_m[i]>0.0{assert_ne!(proxy.source_id[i],u32::MAX);assert!(proxy.source_names.contains_key(&proxy.source_id[i]));assert!(proxy.floor_m[i]<=z[i]+1e-9);}if fixed[i]!=pf_region::t1b::FREE{assert_eq!(proxy.stock_depth_m[i],0.0);}}
    if let Some(dir)=std::env::var_os("VZ_NATIVE_PROXY_REPORT"){let p=PathBuf::from(dir);std::fs::create_dir_all(p.parent().unwrap()).unwrap();assert!(!p.exists(),"preserve prior proxy evidence");std::fs::write(p,serde_json::to_vec_pretty(&serde_json::json!({"input":path,"input_digest":proxy.input_digest,"report":proxy.report,"sources":proxy.source_names,"scope":"actual native proxy binding only; not positive-stock process or geological stock"})).unwrap()).unwrap();}
    eprintln!("native proxy input: {}",proxy.report);
}


#[test]
fn bind_full_native_matrix_and_record_actual_missing_support(){
    let Some(root)=std::env::var_os("VZ_NATIVE_PROXY_ROOT").map(PathBuf::from)else{eprintln!("SKIP full native binding: set VZ_NATIVE_PROXY_ROOT and VZ_NATIVE_PROXY_GROUP_REPORT");return};
    let output=PathBuf::from(std::env::var_os("VZ_NATIVE_PROXY_GROUP_REPORT").expect("explicit persistent report path required"));
    assert!(!output.exists(),"preserve prior group evidence");
    let mut groups=Vec::new();let mut windows=Vec::new();
    for (label,x0,y0,side) in [("central",-1,-1,3),("east",4,0,2),("north",1,4,2)]{
        let n=(side*512+129)as usize;let count=n*n;
        let mut seen=vec![false;count];let mut ground=vec![0.0;count];let mut fixed=vec![255;count];
        let mut stock=vec![0.0;count];let mut material=vec![0.0;count];let mut exposure=vec![0.0;count];let mut source=vec![u32::MAX;count];let mut own=vec![0.0;count];
        let mut names=std::collections::BTreeMap::new();let mut overlaps=0usize;
        for ry in y0..y0+side{for rx in x0..x0+side{
            let w=pf_region::t1b_ve::flank3_load_native_window(&root.join(format!("w_{rx}_{ry}"))).unwrap();
            windows.push(serde_json::json!({"window":[rx,ry],"report":w.report,"input_digest":w.proxy.input_digest}));
            for (&id,name)in &w.proxy.source_names{if let Some(old)=names.insert(id,name.clone()){assert_eq!(old,*name);}}
            for i in 0..w.lat.nodes(){let(gx,gy)=w.lat.g(i);let j=((gy-(y0*512-64))*n as i64+gx-(x0*512-64))as usize;
                if seen[j]{overlaps+=1;assert_eq!(ground[j],w.ground_m[i]);assert_eq!(fixed[j],w.geometry_fixed[i]);assert_eq!(stock[j],w.proxy.stock_depth_m[i]);assert_eq!(material[j],w.proxy.material_m[i]);assert_eq!(exposure[j],w.proxy.exposure_fraction[i]);assert_eq!(source[j],w.proxy.source_id[i]);assert_eq!(own[j],w.own_source_m2[i]);}
                else{seen[j]=true;ground[j]=w.ground_m[i];fixed[j]=w.geometry_fixed[i];stock[j]=w.proxy.stock_depth_m[i];material[j]=w.proxy.material_m[i];exposure[j]=w.proxy.exposure_fraction[i];source[j]=w.proxy.source_id[i];own[j]=w.own_source_m2[i];}
            }
        }}
        assert!(seen.iter().all(|x|*x));
        let mut missing=std::collections::BTreeSet::new();let mut edge_stock=Vec::new();let mut volume=std::collections::BTreeMap::<u32,f64>::new();
        for i in 0..count{
            if stock[i]>0.0{*volume.entry(source[i]).or_default()+=stock[i]*1024.0;}
            let x=i%n;let y=i/n;let gx=x0*512-64+x as i64;let gy=y0*512-64+y as i64;
            if x==0||y==0||x+1==n||y+1==n{
                if stock[i]>0.0{edge_stock.push(serde_json::json!({"gx":gx,"gy":gy,"depth_m":stock[i],"source":source[i]}));}
                if fixed[i]==pf_region::t1b::FREE{for(ox,oy)in pf_region::t1b_surface_domain::D8{
                    let xx=x as i64+ox;let yy=y as i64+oy;if xx<0||yy<0||xx>=n as i64||yy>=n as i64{missing.insert((gy+oy,gx+ox));}
                }}
            }
        }
        let row=serde_json::json!({"group":label,"gx0":x0*512-64,"gy0":y0*512-64,"n":n,"posts":count,"exact_duplicate_reads":overlaps,"positive_stock_posts":stock.iter().filter(|&&x|x>0.0).count(),"unique_proxy_volume_m3":volume.values().sum::<f64>(),"source_volume_m3":volume,"source_names":names,"positive_stock_edge_posts":edge_stock.len(),"edge_stock":edge_stock,"unsupplied_d8_neighbours":missing.len(),"unsupplied_neighbours":missing,"current_complete_graph_api_closed":missing.is_empty(),"note":"Missing full-D8 support is a conservative API limitation, not proof every outside post physically influences this process. Mutable selective MFD/donor/flat/bank/deposition closure must replace that criterion before this can run. Positive edge stock is real supplied input and was not zeroed."});
        eprintln!("native group {label}: {} positive stock, {} positive edge posts, {} unsupplied D8 neighbours",row["positive_stock_posts"],row["positive_stock_edge_posts"],row["unsupplied_d8_neighbours"]);groups.push(row);
    }
    std::fs::create_dir_all(output.parent().unwrap()).unwrap();
    std::fs::write(output,serde_json::to_vec_pretty(&serde_json::json!({"windows":windows,"groups":groups,"scope":"actual full17 native proxy binding and source support; no evolution or mutable-domain qualification"})).unwrap()).unwrap();
}

#[test]
fn native_input_rejects_perturbed_real_own_source(){
    let(Some(path),Some(root))=(std::env::var_os("VZ_NATIVE_PROXY_WINDOW").map(PathBuf::from),std::env::var_os("VZ_NATIVE_PROXY_FALSIFIERS").map(PathBuf::from))else{eprintln!("SKIP native input falsifier: set explicit native window and persistent evidence root");return};
    assert!(!root.exists(),"preserve prior falsifier evidence");std::fs::create_dir_all(root.join("w_bad/fields")).unwrap();
    let target=root.join("w_bad");std::fs::copy(path.join("window.json"),target.join("window.json")).unwrap();
    for entry in std::fs::read_dir(path.join("fields")).unwrap(){let entry=entry.unwrap();std::os::unix::fs::symlink(entry.path(),target.join("fields").join(entry.file_name())).unwrap();}
    // Use a separate changed file and a correspondingly changed locator. The
    // canonical header is rehashed, but the pinned actual source bytes are not.
    let mut bytes=std::fs::read(path.join("fields/runoff_weight.f32")).unwrap();let value=f32::from_le_bytes(bytes[..4].try_into().unwrap());bytes[..4].copy_from_slice(&(value+1.0).to_le_bytes());
    std::fs::write(target.join("fields/perturbed_runoff.f32"),bytes).unwrap();
    let mut j:serde_json::Value=serde_json::from_slice(&std::fs::read(target.join("window.json")).unwrap()).unwrap();
    for f in j["header"]["fields"].as_array_mut().unwrap(){if f["name"]=="runoff_weight"{f["file"]="fields/perturbed_runoff.f32".into();}}
    j["digest"]=pf_contracts::artifact::digest_canonical(&j["header"]).unwrap().into();std::fs::write(target.join("window.json"),serde_json::to_vec(&j).unwrap()).unwrap();
    let error=match pf_region::t1b_ve::flank3_load_native_window(&target){Ok(_)=>panic!("mutated native own source admitted"),Err(e)=>e};assert!(error.contains("runoff_weight size/hash mismatch"),"{error}");
    std::fs::write(root.join("result.json"),serde_json::to_vec_pretty(&serde_json::json!({"actual_native_input":path,"mutation":"one actual runoff sample +1; same expected field content pin, rehashed locator header","rejected":true,"error":error})).unwrap()).unwrap();
}

#[test]
fn mutable_support_uses_native_physical_boundaries_and_preserves_stock(){
    let(Some(root),Some(output))=(std::env::var_os("VZ_NATIVE_SUPPORT_ROOT"),std::env::var_os("VZ_NATIVE_SUPPORT_REPORT"))else{eprintln!("SKIP mutable native support: set VZ_NATIVE_SUPPORT_ROOT and VZ_NATIVE_SUPPORT_REPORT");return};
    let root=PathBuf::from(root);let output=PathBuf::from(output);assert!(!output.exists(),"preserve prior support evidence");
    let mut paths=Vec::new();for(label,x0,y0,side)in [("central",-1,-1,3),("east",4,0,2),("north",1,4,2)]{let _=label;for y in y0..y0+side{for x in x0..x0+side{paths.push(root.join(format!("w_{x}_{y}")));}}}
    if let Some(list)=std::env::var_os("VZ_NATIVE_SUPPORT_PATHS") {paths=serde_json::from_slice(&std::fs::read(list).unwrap()).unwrap();}
    let limit=std::env::var("VZ_NATIVE_SUPPORT_MAX_POSTS").map(|v|v.parse().unwrap()).unwrap_or(8_000_000);
    let report=pf_region::t1b_ve::flank3_native_support(&paths,limit).unwrap();
    assert!(report.positive_stock_posts>0);assert_eq!(report.complete,report.missing.is_empty());assert_eq!(report.complete,report.required_windows.is_empty());
    assert!(report.visited_free_posts>=report.positive_stock_posts);assert!(report.touched_fixed_posts>0);
    std::fs::write(output,serde_json::to_vec_pretty(&report).unwrap()).unwrap();
    eprintln!("mutable native support: {} unique posts; {} positive stock; {} connected free; {} native fixed support; {} missing posts in {} owner windows",report.unique_posts,report.positive_stock_posts,report.visited_free_posts,report.touched_fixed_posts,report.missing.len(),report.required_windows.len());
}

#[test]
fn compact_canonical_material_matches_actual_native_fields_and_proxy(){
    let (Some(tile),Some(reference),Some(world),Some(report))=(std::env::var_os("VZ_MATERIAL_TILE"),std::env::var_os("VZ_NATIVE_PROXY_WINDOW"),std::env::var_os("VZ_MATERIAL_WORLD_JSON"),std::env::var_os("VZ_MATERIAL_REPORT"))else{eprintln!("SKIP canonical material identity: set explicit material/native/world/report paths");return};
    let tile=PathBuf::from(tile);let reference=PathBuf::from(reference);let report=PathBuf::from(report);assert!(!report.exists());
    let mut h:serde_json::Value=serde_json::from_slice(&std::fs::read(tile.join("material.json")).unwrap()).unwrap();let digest=h.as_object_mut().unwrap().remove("digest").unwrap();assert_eq!(digest,pf_region::t1b_surface_domain::content_digest(&h).unwrap());
    let rect:pf_region::t1b_surface_completion::Rect=serde_json::from_value(h["rect"].clone()).unwrap();let count=rect.count().unwrap();
    let w=pf_region::t1b_ve::flank3_load_native_window(&reference).unwrap();assert_eq!(h["runoff_reference_m_s"],w.runoff_reference_m_s);assert_eq!(h["threshold_m3_s"],w.threshold_m3_s);let map:Vec<usize>=(0..count).map(|i|{let p=rect.id(i);w.lat.idx(p.gx,p.gy).unwrap()}).collect();
    let mut fields=std::collections::BTreeMap::new();for d in h["fields"].as_array().unwrap(){let name=d["file"].as_str().unwrap();let b=std::fs::read(tile.join(name)).unwrap();assert_eq!(d["bytes"],b.len());assert_eq!(d["digest"],blake3::hash(&b).to_hex().to_string());fields.insert(name.to_string(),b);}
    let old:serde_json::Value=serde_json::from_slice(&std::fs::read(reference.join("window.json")).unwrap()).unwrap();
    for name in ["height_mm","process_geometry_fixed","process_native_link","volcano_id","process_edifice_style","process_edifice_unit","runoff_weight"]{
        let d=old["header"]["fields"].as_array().unwrap().iter().find(|d|d["name"]==name).unwrap();let file=d["file"].as_str().unwrap();let bytes=std::fs::read(reference.join(file)).unwrap();let name=std::path::Path::new(file).file_name().unwrap().to_str().unwrap();let size=bytes.len()/w.lat.nodes();let expected:Vec<u8>=map.iter().flat_map(|&i|bytes[i*size..(i+1)*size].iter().copied()).collect();assert!(fields[name]==expected,"native material field {name}");
    }
    for(name,v)in [("proxy_stock_depth_m.f64",&w.proxy.stock_depth_m),("proxy_material_m.f64",&w.proxy.material_m),("proxy_exposure_fraction.f64",&w.proxy.exposure_fraction),("proxy_floor_m.f64",&w.proxy.floor_m)]{let expected:Vec<u8>=map.iter().flat_map(|&i|v[i].to_le_bytes()).collect();assert!(fields[name]==expected,"proxy field {name}");}
    let expected:Vec<u8>=map.iter().flat_map(|&i|w.proxy.source_id[i].to_le_bytes()).collect();assert!(fields["proxy_source_id.u32"]==expected,"proxy source IDs");
    let cfg:serde_json::Value=serde_json::from_slice(&std::fs::read(world).unwrap()).unwrap();let path=|k:&str|std::path::Path::new(cfg[k].as_str().unwrap());
    let site=pf_region::t0::SiteRecord::load(path("site")).unwrap();let t0=pf_region::t0::T0::load(path("run"),path("subgrid"),&site,10.0).unwrap();
    let chart=vz_tilekey::chart::Chart::new(serde_json::from_value(h["chart"].clone()).unwrap()).unwrap();let mut cells=Vec::new();let mut basins=Vec::new();let mut bm=std::collections::BTreeMap::new();
    for i in 0..count{let p=rect.id(i);let(la,lo)=chart.inverse_mm(p.gx*32000,p.gy*32000).unwrap();let c=t0.cell_of_deg(la*vz_tilekey::chart::RAD2DEG,lo*vz_tilekey::chart::RAD2DEG);cells.extend_from_slice(&c.to_le_bytes());let basin=*bm.entry(c).or_insert_with(||t0.final_basin(c).unwrap());basins.extend_from_slice(&basin.to_le_bytes());}
    assert!(fields["native_t0_cell.u32"]==cells,"native T0 cell IDs");assert!(fields["native_t0_basin.u32"]==basins,"native raw T0 basin IDs");
    let positive=map.iter().filter(|&&i|w.proxy.stock_depth_m[i]>0.0).count();assert!(positive>0);
    let result=serde_json::json!({"material_tile":tile,"native_reference":reference,"posts":count,"positive_stock_posts":positive,"raw_fields_exact":7,"proxy_fields_exact":5,"raw_t0_fields_exact":2,"field_bytes":fields.values().map(Vec::len).sum::<usize>(),"mutable_support_closed":false,"positive_evolution_run":false});
    std::fs::write(report,serde_json::to_vec_pretty(&result).unwrap()).unwrap();eprintln!("canonical material identity {result}");
}
