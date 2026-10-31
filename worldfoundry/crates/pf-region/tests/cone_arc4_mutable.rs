use pf_region::{t1b_surface_area::{Grid,Sample,qualify_mutable},t1b_ve::material::CompactInput};
use std::path::PathBuf;
#[test]
fn actual_mutable_dependencies_and_unknown_support_falsifier(){
    let (Some(tile),Some(out))=(std::env::var_os("VZ_ARC4_MATERIAL"),std::env::var_os("VZ_ARC4_MUTABLE_REPORT"))else{eprintln!("SKIP: VZ_ARC4_MATERIAL and VZ_ARC4_MUTABLE_REPORT");return};
    let input=CompactInput::load(&PathBuf::from(tile)).unwrap();let out=PathBuf::from(out);assert!(!out.exists());
    let mut g=Grid::new(input.rect,input.header["native_context_pin"].as_str().unwrap().into(),input.height_mm.len()).unwrap();
    let samples:Vec<_>=(0..input.height_mm.len()).map(|i|Sample{id:input.rect.id(i),height_mm:input.height_mm[i],fixed:input.fixed[i],link:input.links[i],own_m2:if input.fixed[i]==0{1024.0*input.runoff[i] as f64}else{0.0}}).collect();g.ingest(&samples).unwrap();
    let mobile:Vec<_>=(0..g.z.len()).filter(|&i|input.proxy_depth[i]>0.0).collect();assert!(!mobile.is_empty());
    let r=qualify_mutable(&g,&mobile).unwrap();
    let mut edges=0;for i in 0..g.z.len(){let sum:f64=r.weights[i].iter().sum();if g.fixed[i]==0&&r.potential_m[i].is_finite(){assert!((sum-1.0).abs()<1e-12);edges+=r.weights[i].iter().filter(|&&w|w>0.0).count();}else{assert_eq!(sum,0.0);}}
    // Find an actual complete subproblem inside the supplied native tile. This
    // tests positive qualification as well as refusals; no synthetic terrain.
    let valid=mobile.iter().copied().find(|&start|{
        if (0..8).any(|k|g.rect.neighbour(start,k).is_none_or(|j|r.bad_root[j]!=u32::MAX)){return false;}
        let mut todo=vec![start];let mut seen=std::collections::BTreeSet::new();
        while let Some(i)=todo.pop(){if !seen.insert(i){continue;}if r.bad_root[i]!=u32::MAX{return false;}
            for k in 0..8{if r.weights[i][k]>0.0{let Some(j)=g.rect.neighbour(i,k)else{return false};todo.push(j);}}
        }true
    }).expect("actual tile should contain at least one complete mutable target");
    let valid_route=qualify_mutable(&g,&[valid]).unwrap();assert!(valid_route.certificate.pass);
    // Remove a real immediate neighbour of an actual positive-stock post.
    let (&target,absent)=mobile.iter().find_map(|i|(0..8).filter_map(|k|g.rect.neighbour(*i,k)).find(|&j|j!=*i).map(|j|(i,j))).unwrap();
    g.present[absent]=false;
    let missing=qualify_mutable(&g,&[target]).unwrap();assert!(!missing.certificate.pass);assert!(missing.certificate.missing.contains(&g.rect.id(absent)));
    // The accepted proof must change with an actual proposed geometry revision.
    g.present[absent]=true;g.z[target]+=0.001;
    let changed=qualify_mutable(&g,&mobile).unwrap();assert_ne!(r.certificate.geometry_digest,changed.certificate.geometry_digest);assert_ne!(r.certificate.route_digest,changed.certificate.route_digest);
    std::fs::write(out,serde_json::to_vec_pretty(&serde_json::json!({"certificate":r.certificate,"positive_branches":edges,"complete_actual_target":g.rect.id(valid),"complete_target_certificate":valid_route.certificate,"missing_actual_neighbour_rejected":true,"changed_geometry_requalified":true,"scope":"Actual material support qualification, no physical step acceptance or coarse composition claimed"})).unwrap()).unwrap();
}
