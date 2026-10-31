//! Native output acceptance, enabled explicitly by VZ_E5B_INPUT. No substitute terrain.
use std::path::PathBuf;
fn ints(p:PathBuf)->Vec<i32>{std::fs::read(p).unwrap().chunks_exact(4).map(|b|i32::from_le_bytes(b.try_into().unwrap())).collect()}
#[test]
fn final_cliffs_and_source_export_accounts() {
    let Ok(p)=std::env::var("VZ_E5B_INPUT")else{eprintln!("SKIP E5B native publication: set VZ_E5B_INPUT");return};
    let m:serde_json::Value=serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap();
    let mut total=0;
    for dir in m["windows"].as_array().unwrap() {
        let w=PathBuf::from(dir.as_str().unwrap());let doc:serde_json::Value=serde_json::from_slice(&std::fs::read(w.join("window.json")).unwrap()).unwrap();let h=&doc["header"];
        let n=h["lattice"]["n"].as_u64().unwrap() as usize;let hp=h["lattice"]["halo_posts"].as_u64().unwrap() as usize;
        let z=ints(w.join("fields/height_mm.i32"));let rd=|s:&str|std::fs::read(w.join(format!("fields/{s}.u8"))).unwrap();
        let (rock,rg,fixed,over)=(rd("exposed_litho"),rd("ct_regime"),rd("fixed"),rd("over"));
        let c=pf_region::t1b_ct::CtParamsV1::v1();let mut count=0;
        for y in hp..n-hp {for x in hp..n-hp {let i=y*n+x;if fixed[i]!=0||over[i]!=0{continue;}
            let Some(k)=c.classes.iter().find(|k|k.code==rock[i]&&k.mode!=0)else{continue};
            let mut slope=0.0_f64;for dy in -1isize..=1{for dx in -1isize..=1{if dx==0&&dy==0{continue;}let j=(y as isize+dy) as usize*n+(x as isize+dx)as usize;
                slope=slope.max((z[i] as f64-z[j]as f64).abs()/32000./(dx as f64).hypot(dy as f64));}}
            if slope>=k.cliff_tan_permille as f64/1000. {assert_eq!(rg[i],2,"cliff at {x},{y}");count+=1;}
        }}
        let hs=&h["prep"]["net"]["hs"];let cs=&hs["contrast"];assert_eq!(cs["local_ledger_residual_mm"],0);assert_eq!(cs["talus_stored_m3"],0.);
        let mut sum=[0i128;3];for a in cs["export_accounts"].as_array().unwrap(){for(k,key)in ["cap_weathering_mm3","weak_weathering_mm3","collapsed_cap_mm3"].iter().enumerate(){let v=a[key].as_str().unwrap().parse::<i128>().unwrap();assert!(v>=0);sum[k]+=v;}}
        for(k,key)in ["cap_weathering_m3","weak_weathering_m3","collapse_export_m3"].iter().enumerate(){assert!((sum[k]as f64/1e9-cs[key].as_f64().unwrap()).abs()<1e-5);}
        assert!(cs["export_reservoir"].as_str().unwrap().contains("unresolved_external_sediment"));
        println!("E5B {} eligible_cliffs={} missing=0 exported_mm3={sum:?}",w.display(),count);total+=count;
    }
    assert!(total>0,"actual Moab cliffs required");
}

#[test]
fn existing_recipe_keysets_are_unchanged() {
    let Ok(p)=std::env::var("VZ_E5B_INPUT")else{eprintln!("SKIP E5B recipes: set VZ_E5B_INPUT");return};
    let m:serde_json::Value=serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap();let mut keys=std::collections::BTreeSet::new();let mut n=0;
    for path in m["legacy_recipes"].as_array().unwrap(){let v:serde_json::Value=serde_json::from_slice(&std::fs::read(path.as_str().unwrap()).unwrap()).unwrap();
        let params=&v["header"]["recipe"]["params"];let p=params.get("rivernet").unwrap_or(params);
        // Existing recipe adapter: zero f6_reference is omitted by Serialize,
        // but legacy Deserialize requires it. Re-serialization must omit it again.
        let mut decoded=p.clone();if decoded.get("f6_reference").is_none(){decoded["f6_reference"]=0.into();}
        let q:pf_region::t1b_rivernet::RiverNetParamsV1=serde_json::from_value(decoded).unwrap();let r=serde_json::to_value(q).unwrap();assert_eq!(&r,p);
        keys.insert(serde_json::to_string(p).unwrap());n+=1;
    }
    assert!(n>0);println!("E5B legacy recipes={n} distinct_param_sets={} exact",keys.len());
}
