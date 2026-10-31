use pf_detail::input::T1bWindow;
use std::path::PathBuf;
fn ints(p:PathBuf)->Vec<i32>{std::fs::read(p).unwrap().chunks_exact(4).map(|b|i32::from_le_bytes(b.try_into().unwrap())).collect()}
#[test]
fn real_moab_parent_posts_seams_cliff_and_provenance() {
    let Ok(p)=std::env::var("VZ_E5B_INPUT")else{eprintln!("SKIP E5B native detail: set VZ_E5B_INPUT");return};
    let m:serde_json::Value=serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap();let w=T1bWindow::load(&PathBuf::from(m["windows"][0].as_str().unwrap())).unwrap();assert!(w.evolved.is_some());
    let root=PathBuf::from(m["detail"].as_str().unwrap());let mut nparent=0;let mut ncliff=0;
    for ix in [25,26]{let d=root.join(format!("tiles/t_{ix}_33"));let z=ints(d.join("height_mm.i32"));let rg=std::fs::read(d.join("ct_regime.u8")).unwrap();let l=std::fs::read(d.join("litho.u8")).unwrap();let ids=ints(d.join("prov_id.u32"));
        for y in 0..129{for x in 0..129{let i=y*129+x;let (xx,yy)=((ix*256+x*2) as f64,(33*256+y*2)as f64);
            if x%16==0&&y%16==0{assert_eq!(z[i],pf_detail::to_mm(w.h((xx/32.)as i64,(yy/32.)as i64)));nparent+=1;}
            assert_eq!(z[i],pf_detail::contrast::height_mm(&w,(xx/2.)as i64,(yy/2.)as i64));
            assert_ne!(ids[i],0,"real provenance");
            if x>0&&x<128&&y>0&&y<128&&l[i]==6 {let mut slope=0.0_f64;for dy in -1isize..=1{for dx in -1isize..=1{if dx==0&&dy==0{continue;}let j=(y as isize+dy)as usize*129+(x as isize+dx)as usize;slope=slope.max((z[i]as f64-z[j]as f64).abs()/2000./(dx as f64).hypot(dy as f64));}}
                if slope>=3.732{assert_eq!(rg[i],2);ncliff+=1;}
            }
        }}
        if ix==25{let a=z[48*129+80];let b=z[64*129+80];assert_eq!(b-a,276527);assert!(f64::atan((b-a)as f64/32000.).to_degrees()>75.);}
    }
    let a=ints(root.join("tiles/t_25_33/height_mm.i32"));let b=ints(root.join("tiles/t_26_33/height_mm.i32"));for y in 0..129{assert_eq!(a[y*129+128],b[y*129]);}
    assert!(ncliff>0);println!("E5B native fine parent_matches={nparent} seam=129 cliffs={ncliff} missing=0");
}
