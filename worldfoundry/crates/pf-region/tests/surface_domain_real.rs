//! Real native pre-fill component through an explicit caller-supplied artifact.
//! No cached project paths or synthetic terrain gates.
use pf_region::t1b_surface_domain::*;
use pf_region::t1b_ls::route_detail_surface;
fn input()->Option<Tile> {
    let Some(path)=std::env::var_os("VZ_SURFACE_REAL_TILE") else {eprintln!("SKIP real surface domain: set VZ_SURFACE_REAL_TILE to a qualified native component tile");return None;};
    let mut tile:Tile=serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    // The fixture is extracted from native capture, not an accepted upstream
    // stage. Pin its exact local fields here for API/falsifier tests only.
    let (g,s,b)=component_digests(&tile.posts).unwrap();
    tile.pin.ground=g;tile.pin.sources=s;tile.pin.boundary=b;
    Some(tile)
}
fn domain(tile:Tile)->Domain {Domain::assemble(&tile.pin.clone(),tile.policy,vec![tile],5_000_000).unwrap()}
#[test]
fn real_partition_permutation_and_legacy_minimax() {
    let Some(tile)=input() else{return};
    let d=domain(tile.clone());
    let r=route_detail_surface(&d,Mode::SingleFlowV1).unwrap();
    for split in [7,11] {
        let mut tiles=vec![Tile {pin:tile.pin.clone(),policy:tile.policy,posts:Vec::new()};split];
        for (i,p) in tile.posts.iter().rev().enumerate() {tiles[i%split].posts.push(p.clone());}
        tiles.reverse();
        let other=Domain::assemble(&tile.pin,tile.policy,tiles,5_000_000).unwrap();
        assert_eq!(d.digest(),other.digest());
        assert_eq!(serde_json::to_vec(&r).unwrap(),serde_json::to_vec(&route_detail_surface(&other,Mode::SingleFlowV1).unwrap()).unwrap());
    }
    // Re-embedding the complete real component in its original square produces
    // exactly the legacy pinned minimax heights at all component posts. Outside
    // posts are isolated at +infinity; this is an oracle domain extension, not
    // manufactured terrain claimed as game data.
    let x0=tile.posts.iter().map(|p|p.id.gx).min().unwrap();
    let y0=tile.posts.iter().map(|p|p.id.gy).min().unwrap();
    let n=((tile.posts.iter().map(|p|p.id.gx).max().unwrap()-x0).max(tile.posts.iter().map(|p|p.id.gy).max().unwrap()-y0)+1) as usize;
    let lat=pf_region::t1b_rivernet::EvalLat {gx0:x0,gy0:y0,n,dx_mm:tile.pin.spacing_mm};
    let mut z=vec![f64::INFINITY;n*n];let mut fixed=vec![false;n*n];
    for p in &tile.posts {let i=lat.idx(p.id.gx,p.id.gy).unwrap();z[i]=p.ground_m;fixed[i]=p.boundary.is_some();}
    let mut sm=pf_region::t1b_seam::SeamTaint::new(n,0);
    pf_region::t1b_seam::pinned_fill(&lat,&mut z,&fixed,&mut sm,"real-oracle");
    for (i,p) in d.posts().iter().enumerate() {assert_eq!(r.routing_m[i].to_bits(),z[lat.idx(p.id.gx,p.id.gy).unwrap()].to_bits());}
    eprintln!("real component: {} posts, two partitions, minimax exact",d.posts().len());
}
#[test]
fn real_missing_duplicate_stale_boundary_and_resource_falsifiers() {
    let Some(tile)=input() else{return};
    let active=tile.posts.iter().position(|p|p.boundary.is_none()).unwrap();
    let p=tile.posts[active].id;
    let neighbour=tile.posts.iter().position(|q|(0..8).any(|k|p.neighbour(k)==q.id)).unwrap();
    let mut missing=tile.clone();missing.posts.remove(neighbour);
    assert!(matches!(Domain::assemble(&tile.pin,tile.policy,vec![missing],5_000_000),Err(DomainError::NeedsDomain{..})));
    let mut duplicate=tile.clone();duplicate.posts.push(duplicate.posts[active].clone());
    assert!(matches!(Domain::assemble(&tile.pin,tile.policy,vec![duplicate],5_000_000),Err(DomainError::DuplicateSource{..})));
    let mut stale=tile.clone();stale.pin.predecessor.push('x');
    assert!(matches!(Domain::assemble(&tile.pin,tile.policy,vec![stale],5_000_000),Err(DomainError::Stale{..})));
    let mut doubled=tile.clone();doubled.posts[active].source_m2+=doubled.posts[active].support_m2;
    assert!(matches!(Domain::assemble(&tile.pin,tile.policy,vec![doubled],5_000_000),Err(DomainError::Stale{field}) if field=="sources"));
    let mut bad=tile.clone();bad.posts.iter_mut().find(|p|p.boundary.is_some()).unwrap().boundary.as_mut().unwrap().kind=pf_region::t1b::RING;
    assert!(matches!(Domain::assemble(&tile.pin,tile.policy,vec![bad],5_000_000),Err(DomainError::Invalid{..})));
    assert!(matches!(Domain::assemble(&tile.pin,tile.policy,vec![tile.clone()],tile.posts.len()-1),Err(DomainError::Resource{..})));
}
#[test]
fn real_mfd_flux_and_potential() {
    let Some(mut tile)=input() else{return};
    // This test requires a real single-basin, unconduited FREE component; the
    // export script checks that precondition, not an invented basin assignment.
    assert!(tile.posts.iter().all(|p|p.basin==tile.posts[0].basin));
    tile.policy=Policy::NativeDetailV1;
    let d=domain(tile);
    let r=route_detail_surface(&d,Mode::MultipleFlowSlope11V1).unwrap();
    let mut geometric:Vec<f64>=d.posts().iter().map(|p|if p.boundary.is_none(){p.support_m2}else{0.0}).collect();
    let mut effective:Vec<f64>=d.posts().iter().map(|p|p.source_m2+p.external_m2).collect();
    for &ii in r.downstream_first.iter().rev() {
        let i=ii as usize;
        assert!((r.edges[i].iter().map(|e|e.weight).sum::<f64>()-1.0).abs()<1e-12);
        for e in &r.edges[i] {
            let j=e.receiver as usize;
            assert!(r.routing_m[i]>r.routing_m[j] || (r.routing_m[i]==r.routing_m[j] && r.rank[i]>r.rank[j]));
            assert_eq!(d.posts()[i].basin,d.posts()[j].basin);
            geometric[j]+=geometric[i]*e.weight;effective[j]+=effective[i]*e.weight;
        }
    }
    assert_eq!(geometric,r.geometric_area_m2);assert_eq!(effective,r.effective_area_m2);
    for b in &r.balances {
        assert!((b.geometric_source_m2-b.geometric_export_m2).abs()<=1e-9*b.geometric_source_m2.max(1.0));
        assert!((b.effective_source_m2-b.effective_export_m2).abs()<=1e-9*b.effective_source_m2.max(1.0));
    }
    eprintln!("real MFD: {} posts, capacity {} m3, balances {:?}",d.posts().len(),r.fill_capacity_m3,r.balances);
}

#[test]
fn real_artifact_roundtrip_and_changed_content_fails() {
    let Some(tile)=input() else{return};
    let path=std::path::PathBuf::from(std::env::var_os("VZ_SURFACE_REAL_TILE").unwrap());
    // Sandbox PIDs repeat across invocations. Allocate an unused ordinal under
    // that PID rather than overwriting or deleting earlier test evidence.
    let mut ordinal=0u64;
    let out=loop {
        let candidate=path.parent().unwrap().join(format!("artifact-test-{}-{ordinal}",std::process::id()));
        match std::fs::create_dir(&candidate) {
            Ok(())=>break candidate,
            Err(e) if e.kind()==std::io::ErrorKind::AlreadyExists=>ordinal+=1,
            Err(e)=>panic!("cannot allocate preserved test evidence: {e}"),
        }
    };
    let bytes=serde_json::to_vec(&tile).unwrap();
    let digest=blake3::hash(&bytes).to_hex().to_string();
    std::fs::write(out.join("tile.json"),&bytes).unwrap();
    let request=serde_json::json!({"pin":tile.pin,"policy":tile.policy,"mode":"SingleFlowV1","max_posts":5000000,"tiles":[{"path":"tile.json","digest":digest}]});
    std::fs::write(out.join("request.json"),serde_json::to_vec_pretty(&request).unwrap()).unwrap();
    build_from_request(&out.join("request.json"),&out.join("candidate")).unwrap();
    let certificate:serde_json::Value=serde_json::from_slice(&std::fs::read(out.join("candidate/certificate.json")).unwrap()).unwrap();
    assert_eq!(certificate["closed"],true);assert_eq!(certificate["accepted_for_publication"],false);
    let mut changed=bytes;changed.push(b' ');std::fs::write(out.join("tile.json"),changed).unwrap();
    assert!(matches!(build_from_request(&out.join("request.json"),&out.join("rejected")),Err(DomainError::Stale{..})));
    assert!(out.join("rejected/failure.json").exists());
    eprintln!("real artifact and stale-byte rejection: {}",out.display());
}

#[test]
fn real_owner_registry_coverage_and_graph_falsifiers(){
    let Some(mut tile)=input()else{return};tile.policy=Policy::NativeDetailV1;let d=domain(tile);
    let r=route_detail_surface(&d,Mode::MultipleFlowSlope11V1).unwrap();
    // Deliberately low API-test threshold exercises graph registration on the
    // actual small component. It is not the native F2 threshold or a world gate.
    let reg=pf_region::t1b_ls::register_detail_channels(&d,&r,1e-8,1e-12).unwrap();
    let a=pf_region::t1b_surface_channels::audit(&d,&r,&reg).unwrap();
    assert_eq!(a.supra_threshold_posts,d.posts().iter().filter(|p|p.boundary.is_none()).count());
    assert!(a.connected_edges>0);assert!(a.pass);
    let mut gap=reg.clone();let p=gap.posts.iter_mut().find(|p|!p.edges.is_empty()).unwrap();p.edges.pop();gap.digest.clear();gap.digest=content_digest(&gap).unwrap();
    assert!(pf_region::t1b_surface_channels::audit(&d,&r,&gap).is_err());
    let mut false_outlet=reg.clone();let boundary=false_outlet.posts.iter().find_map(|p|p.boundary.clone()).unwrap();
    false_outlet.posts.iter_mut().find(|p|p.boundary.is_none()).unwrap().boundary=Some(boundary);false_outlet.digest.clear();false_outlet.digest=content_digest(&false_outlet).unwrap();
    assert!(pf_region::t1b_surface_channels::audit(&d,&r,&false_outlet).is_err());
    let mut stale=reg.clone();stale.domain_digest.push('x');stale.digest.clear();stale.digest=content_digest(&stale).unwrap();
    assert!(matches!(pf_region::t1b_surface_channels::audit(&d,&r,&stale),Err(DomainError::Stale{..})));
    let mut broken=r.clone();broken.edges.iter_mut().find(|e|!e.is_empty()).unwrap()[0].weight*=0.5;broken.digest.clear();broken.digest=content_digest(&broken).unwrap();
    assert!(pf_region::t1b_surface_channels::audit_route(&d,&broken).is_err());
    let projection=pf_region::t1b_surface_channels::guide_projection(&d,&r,&reg);
    if a.ambiguous_legacy_attachments>0{assert!(projection.is_err());}else{let p=projection.unwrap();assert_eq!(p.len(),a.supra_threshold_posts);assert!(p.values().all(|v|v.0==pf_region::t1b::GUIDE));}
    eprintln!("real owner registry (API-test threshold only): {:?}",a);
}

#[test]
fn owner_flank_dry_stock_keeps_real_routes_and_rejects_unidentified_stock(){
    let Some(mut tile)=input()else{return};tile.policy=Policy::NativeDetailV1;
    // Native millimetre representation of the captured geometry. The explicit
    // zero-removable-stock condition is an API falsifier, not inferred geology
    // and not the real C7 120/240-step morphology/refinement experiment.
    for p in &mut tile.posts{p.ground_m=vz_tilekey::chart::quantize_mm(p.ground_m*1000.0)as f64/1000.0;}
    let(g,s,b)=component_digests(&tile.posts).unwrap();tile.pin.ground=g;tile.pin.sources=s;tile.pin.boundary=b;let d=domain(tile);let n=d.posts().len();
    let zero=vec![0.0;n];let exposure=vec![1.0;n];let floors:Vec<f64>=d.posts().iter().map(|p|p.ground_m).collect();let source=vec![u32::MAX;n];let names=std::collections::BTreeMap::new();
    let proxy=||pf_region::t1b_ve::Flank3Proxy{stock_depth_m:&zero,material_m:&zero,exposure_fraction:&exposure,floor_m:&floors,source_id:&source,source_names:&names};
    let mut reports=0;
    let a=pf_region::t1b_ve::flank3_evolve_owner(&d,proxy(),120,1e-8,1e-12,&mut |_|{reports+=1;Ok(())}).unwrap();assert_eq!(reports,120);
    for pair in a.steps.windows(2){assert_eq!(pair[0].route_digest,pair[1].input_route_digest);}
    assert_eq!(a.steps.last().unwrap().accepted_domain_digest,a.report["physical_pre_projection_revision"].as_str().unwrap());
    assert_eq!(a.physical_ground_m,floors);assert_eq!(a.report["exported_m3"],0.0);assert!(a.sources.is_empty());
    let text=include_str!("../src/t1b_ve.rs");let begin=text.find("pub const FLANK3_OWNER_RULE").unwrap()+"pub const FLANK3_OWNER_RULE".len();
    assert_eq!(a.report["kernel_digest"],blake3::hash(text[begin..].as_bytes()).to_hex().to_string());
    let b=pf_region::t1b_ve::flank3_evolve_owner(&d,proxy(),240,1e-8,1e-12,&mut |_|Ok(())).unwrap();assert_eq!(a.published_height_mm,b.published_height_mm);assert_eq!(a.route.effective_area_m2,b.route.effective_area_m2);
    let mut unidentified=zero.clone();unidentified[d.posts().iter().position(|p|p.boundary.is_none()).unwrap()]=1.0;
    let bad=pf_region::t1b_ve::Flank3Proxy{stock_depth_m:&unidentified,material_m:&zero,exposure_fraction:&exposure,floor_m:&floors,source_id:&source,source_names:&names};
    assert!(pf_region::t1b_ve::flank3_evolve_owner(&d,bad,120,1e-8,1e-12,&mut |_|Ok(())).is_err());
    eprintln!("owner flank dry-stock falsifier: {n} real geometry posts; 120/240 owner reroutes, zero geometry/sediment change; not native volcano refinement acceptance");
}

#[test]
fn physical_policy_keeps_basin_metadata_without_inventing_walls(){
    let Some(mut tile)=input()else{return};tile.policy=Policy::NativePhysicalV1;
    let original=domain(tile.clone());let a=route_detail_surface(&original,Mode::MultipleFlowSlope11V1).unwrap();
    // Mutate provenance labels of actual native geometry. These labels must not
    // change the solve under the physical policy; transfer accounting must see
    // them. No changed/fabricated terrain is used as an acceptance fixture.
    for p in &mut tile.posts {if p.boundary.is_none(){p.basin=p.basin.wrapping_add(1);}}
    let(g,s,b)=component_digests(&tile.posts).unwrap();tile.pin.ground=g;tile.pin.sources=s;tile.pin.boundary=b;
    let changed=domain(tile);let b=route_detail_surface(&changed,Mode::MultipleFlowSlope11V1).unwrap();
    assert_eq!(a.routing_m,b.routing_m);assert_eq!(a.effective_area_m2,b.effective_area_m2);
    assert_eq!(serde_json::to_value(&a.edges).unwrap(),serde_json::to_value(&b.edges).unwrap());assert!(!b.basin_transfers.is_empty());
    let reg=pf_region::t1b_ls::register_detail_channels(&changed,&b,1e-8,1e-12).unwrap();
    let product=pf_region::t1b_surface_channels::native_graph_product(&changed,&b,&reg).unwrap();
    pf_region::t1b_surface_channels::audit_native_graph_product(&changed,&b,&reg,&product).unwrap();
    assert!(!product.coarse_network_composed);assert!(!product.ports.is_empty());
    assert!(product.registry.posts.iter().any(|p|p.native_attachments.len()>1));
    let mut false_transfer=b.clone();false_transfer.basin_transfers[0].effective_m2+=1.0;false_transfer.digest.clear();false_transfer.digest=content_digest(&false_transfer).unwrap();
    assert!(pf_region::t1b_surface_channels::audit_route(&changed,&false_transfer).is_err());
    assert!((product.local_effective_source_m2-product.ports.iter().map(|p|p.effective_export_m2).sum::<f64>()).abs()<1e-9*product.local_effective_source_m2.max(1.0));
    let mut false_balance=b.clone();false_balance.balances[0].effective_source_m2+=1.0;false_balance.digest.clear();false_balance.digest=content_digest(&false_balance).unwrap();
    assert!(pf_region::t1b_surface_channels::audit_route(&changed,&false_balance).is_err());
    let mut bad=product.clone();bad.ports[0].effective_export_m2+=1.0;bad.digest.clear();bad.digest=content_digest(&bad).unwrap();
    assert!(pf_region::t1b_surface_channels::audit_native_graph_product(&changed,&b,&reg,&bad).is_err());
    let mut lost=product.clone();let row=lost.registry.posts.iter_mut().find(|p|!p.edges.is_empty()).unwrap();row.edges.pop();lost.registry.digest.clear();lost.registry.digest=content_digest(&lost.registry).unwrap();lost.digest.clear();lost.digest=content_digest(&lost).unwrap();
    assert!(pf_region::t1b_surface_channels::audit_native_graph_product(&changed,&b,&reg,&lost).is_err());
    eprintln!("physical-label mutation: {} actual native posts retain exact routing; {} explicit basin transfers; {} ports; changed port and missing branch rejected",changed.posts().len(),b.basin_transfers.len(),product.ports.len());
}
