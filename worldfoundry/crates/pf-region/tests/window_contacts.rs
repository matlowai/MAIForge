//! Real rnwin/T2 regression. Supply a small JSON path manifest using
//! PF_WINDOW_CONTACTS_FIXTURE; no private project paths or synthetic worlds.
use std::{fs, path::{Path,PathBuf}};
use serde_json::Value;
fn bytes(p:&Path,f:&str)->Vec<u8>{fs::read(p.join(f)).unwrap()}
fn i32s(p:&Path,f:&str)->Vec<i32>{bytes(p,f).chunks_exact(4).map(|c|i32::from_le_bytes(c.try_into().unwrap())).collect()}
fn floats(p:&Path,f:&str)->Vec<f32>{bytes(p,f).chunks_exact(4).map(|c|f32::from_le_bytes(c.try_into().unwrap())).collect()}
#[test]
fn real_window_contacts_match_detail_seams_and_legacy() {
    let Some(path)=std::env::var_os("PF_WINDOW_CONTACTS_FIXTURE") else {
        eprintln!("SKIP real contact acceptance: set PF_WINDOW_CONTACTS_FIXTURE to the small path manifest"); return;
    };
    let cfg:Value=serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    let path=|k:&str|PathBuf::from(cfg[k].as_str().unwrap());
    let on=path("on"); let off=path("off"); let old=path("legacy"); let served=path("served"); let neighbor=path("neighbor");
    let (w,h)=pf_region::t1b_window::load(&on).unwrap();
    let (_,off_h)=pf_region::t1b_window::load(&off).unwrap();
    let (_,old_h)=pf_region::t1b_window::load(&old).unwrap();
    let (_,served_h)=pf_region::t1b_window::load(&served).unwrap();
    assert_eq!(off_h,old_h,"key-off header identity");
    assert_eq!(h,served_h,"served/batch header identity");
    assert_eq!(h["recipe_digest"],off_h["recipe_digest"],"publication cannot re-key physics");
    let fields=off_h["fields"].as_array().unwrap();
    for f in fields {let file=f["file"].as_str().unwrap();
        assert_eq!(bytes(&old,file),bytes(&off,file),"key-off {file}");
        assert_eq!(bytes(&off,file),bytes(&on,file),"opt-in changed terrain {file}");
    }
    assert_eq!(h["fields"].as_array().unwrap().len(),fields.len()+3);
    for f in h["fields"].as_array().unwrap() { let file=f["file"].as_str().unwrap();assert_eq!(bytes(&on,file),bytes(&served,file),"served {file}"); }
    let base=i32s(&on,"fields/unit_base_mm.i32");let top=i32s(&on,"fields/unit_top_mm.i32");let layer=bytes(&on,"fields/layer.u8");
    assert_eq!(base.len(),w.lat.n*w.lat.n);assert_eq!(top.len(),base.len());assert_eq!(layer.len(),base.len());
    let mut water=0;
    for i in 0..base.len(){if w.fixed[i]==2 || w.fixed[i]==4 {water+=1;assert_eq!((base[i],top[i],layer[i]),(i32::MIN,i32::MIN,255));}}
    let(mut compared,mut worst_base,mut worst_top)=(0usize,0.0_f64,0.0_f64);
    for t in cfg["tiles"].as_array().unwrap(){
        let p=PathBuf::from(t["path"].as_str().unwrap());let tx=t["x"].as_i64().unwrap();let ty=t["y"].as_i64().unwrap();
        let th=i32s(&p,"height_mm.i32");let hab=floats(&p,"height_above_unit_base_m.f32");let diu=floats(&p,"depth_in_unit_m.f32");let tl=bytes(&p,"layer.u8");
        for j in (0..129).step_by(16){for i in (0..129).step_by(16){
            let k=j*129+i;let x=tx*8+i as i64/16;let y=ty*8+j as i64/16;
            let q=((y-w.lat.gy0) as usize)*w.lat.n+(x-w.lat.gx0) as usize;
            assert_eq!(w.height_mm[q],th[k],"same-post height");assert_eq!(layer[q],tl[k],"same-post layer");
            assert!(hab[k].is_finite() && diu[k].is_finite());assert_ne!(base[q],i32::MIN);assert_ne!(top[q],i32::MIN);
            worst_base=worst_base.max((base[q] as f64-(th[k] as f64-hab[k] as f64*1000.)).abs());
            worst_top=worst_top.max((top[q] as f64-(th[k] as f64+diu[k] as f64*1000.)).abs());compared+=1;
        }}
    }
    // Half-mm quantization plus f32 relative-depth rounding at these measured contacts.
    assert!(compared>0);assert!(worst_base<=0.55,"base error {worst_base} mm");assert!(worst_top<=0.55,"top error {worst_top} mm");
    let(n,_) = pf_region::t1b_window::load(&neighbor).unwrap();
    let mut seam=0;
    for (file,size) in [("fields/unit_base_mm.i32",4),("fields/unit_top_mm.i32",4),("fields/layer.u8",1)]{
        let a=bytes(&on,file);let b=bytes(&neighbor,file);let mut count=0;
        for y in w.lat.gy0.max(n.lat.gy0)..(w.lat.gy0+w.lat.n as i64).min(n.lat.gy0+n.lat.n as i64){
            for x in w.lat.gx0.max(n.lat.gx0)..(w.lat.gx0+w.lat.n as i64).min(n.lat.gx0+n.lat.n as i64){
                let ai=((y-w.lat.gy0)as usize*w.lat.n+(x-w.lat.gx0)as usize)*size;
                let bi=((y-n.lat.gy0)as usize*n.lat.n+(x-n.lat.gx0)as usize)*size;
                assert_eq!(&a[ai..ai+size],&b[bi..bi+size],"seam {file} at {x},{y}");count+=1;
            }
        }assert!(count>0);seam=count;
    }
    println!("CONTACT_ACCEPTANCE old_fields={} water_posts={water} t2_posts={compared} max_base_error_mm={worst_base} max_top_error_mm={worst_top} seam_posts_per_field={seam}",fields.len());
}
