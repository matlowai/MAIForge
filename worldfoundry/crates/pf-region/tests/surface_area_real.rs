use pf_region::t1b_surface_area::{Grid,Sample,solve};
use pf_region::t1b_surface_completion::Rect;
use std::path::PathBuf;
#[test]
fn owner_sfd_matches_native_law_and_keeps_missing_donors_visible(){
    let Some(root)=std::env::var_os("VZ_SURFACE_WITNESS").map(PathBuf::from)else{eprintln!("SKIP native area real test: set VZ_SURFACE_WITNESS");return};
    let directory=root.join("capture_1_5");let dir=std::fs::read_dir(&directory).unwrap().map(|p|p.unwrap().path()).find(|p|p.file_name().unwrap().to_str().unwrap().starts_with("coast_final_")).unwrap();
    let native=std::fs::read_dir(&directory).unwrap().map(|p|p.unwrap().path()).find(|p|p.file_name().unwrap().to_str().unwrap().starts_with("native_")).unwrap();
    let meta:serde_json::Value=serde_json::from_slice(&std::fs::read(dir.join("capture.json")).unwrap()).unwrap();let n=meta["n"].as_u64().unwrap()as usize;
    let rect=Rect{gx0:meta["gx0"].as_i64().unwrap(),gy0:meta["gy0"].as_i64().unwrap(),width:n,height:n};
    let bytes=std::fs::read(dir.join("ground.f64")).unwrap();assert_eq!(blake3::hash(&bytes).to_hex().as_str(),meta["fields"]["ground.f64"].as_str().unwrap());let z:Vec<f64>=bytes.chunks_exact(8).map(|x|f64::from_le_bytes(x.try_into().unwrap())).collect();
    let fixed=std::fs::read(dir.join("fixed.u8")).unwrap();let links=std::fs::read(native.join("chan_link.u32")).unwrap();let own=std::fs::read(native.join("source.f64")).unwrap();
    let samples:Vec<Sample>=(0..n*n).map(|i|Sample{id:rect.id(i),height_mm:vz_tilekey::chart::quantize_mm(z[i]*1000.0)as i32,fixed:fixed[i],link:u32::from_le_bytes(links[i*4..i*4+4].try_into().unwrap()),own_m2:f64::from_le_bytes(own[i*8..i*8+8].try_into().unwrap())}).collect();
    let mut grid=Grid::new(rect,format!("native-real:{}",meta["fields"]["ground.f64"].as_str().unwrap()),2_000_000).unwrap();grid.ingest(&samples).unwrap();
    let targets:Vec<usize>=(256..n-256).flat_map(|y|(256..n-256).map(move|x|y*n+x)).collect();let actual=solve(&grid,&targets).unwrap();
    let legacy=pf_region::t1b_lem::route(n,&grid.z,&grid.fixed,32.0);let mut area:Vec<f64>=samples.iter().map(|p|p.own_m2).collect();for &iu in legacy.order.iter().rev(){let i=iu as usize;let r=legacy.receiver[i];if r!=pf_region::NONE{area[r as usize]+=area[i];}}
    let mut compared=0;for i in 0..n*n{if fixed[i]==pf_region::t1b::FREE{assert_eq!(actual.receiver[i],legacy.receiver[i]);assert_eq!(actual.area[i].to_bits(),area[i].to_bits());compared+=1;}}
    assert!(actual.certificate.unresolved_requested_posts>0);assert!(!actual.certificate.missing.is_empty());
    // Repeated reads do not double sources; conflicting known bytes are rejected.
    grid.ingest(&samples[..100]).unwrap();let first=samples.iter().find(|p|p.fixed==pf_region::t1b::FREE).unwrap();let mut conflict=*first;conflict.own_m2+=1024.0;assert!(grid.ingest(&[conflict]).is_err());
    eprintln!("real native SFD: {compared} free receivers and area values exact; {} requested posts still require real outside inputs",actual.certificate.unresolved_requested_posts);
}
