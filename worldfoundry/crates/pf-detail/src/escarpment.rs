//! The native sparse ESCARP-1 fine surface, consumed by terrain AND collision.
//! All removal happened in T1b. Detail only samples the retained geometry.
use crate::R;
use std::collections::HashMap;
pub struct FineSurface {
    pub meta: serde_json::Value,
    gx0:i64, gy0:i64, n:usize, base:Vec<i32>, delta:HashMap<(i64,i64),i32>,
}
impl FineSurface {
    pub fn load(meta:&serde_json::Value,gx0:i64,gy0:i64,n:usize,base:&[u8],delta:&[u8])->R<Self> {
        if meta["schema"]!="escarp-sparse-e5-2m/v1"||meta["spacing_mm"]!=2000||base.len()!=n*n*4||delta.len()%12!=0||meta["published_changed_fine_posts"].as_u64()!=Some((delta.len()/12)as u64) {
            return Err("ESCARP-1 corrupt or unsupported fine field".into());
        }
        let rd=|v:&[u8]|i32::from_le_bytes(v.try_into().unwrap());
        let mut map=HashMap::new();let mut prev=None;
        for v in delta.chunks_exact(12) {
            let x=rd(&v[..4])as i64;let y=rd(&v[4..8])as i64;let d=rd(&v[8..]);
            if x<gx0*16||x>(gx0+n as i64-1)*16||y<gy0*16||y>(gy0+n as i64-1)*16||d>=0||prev.is_some_and(|p|p>=(y,x)) {
                return Err("ESCARP-1 out of domain, non-removal, duplicate or unsorted fine record".into());
            }
            prev=Some((y,x));map.insert((x,y),d);
        }
        Ok(Self{meta:meta.clone(),gx0,gy0,n,base:base.chunks_exact(4).map(rd).collect(),delta:map})
    }
    pub fn height_mm(&self,px:i64,py:i64)->i32 {
        let gx=px.div_euclid(16)-self.gx0;let gy=py.div_euclid(16)-self.gy0;
        assert!(gx>=0&&gy>=0&&gx+1<self.n as i64&&gy+1<self.n as i64,"ESCARP-1 sample outside retained parent halo");
        let i=gy as usize*self.n+gx as usize;let u=px.rem_euclid(16);let v=py.rem_euclid(16);
        let h=|j|self.base[j]as i64;
        let num=(h(i)*(16-u)+h(i+1)*u)*(16-v)+(h(i+self.n)*(16-u)+h(i+self.n+1)*u)*v;
        let z=if num>=0{(num+128)/256}else{-((-num+128)/256)};
        i32::try_from(z+self.delta.get(&(px,py)).copied().unwrap_or(0)as i64).expect("ESCARP-1 fine height overflow")
    }
}

/// Exact (to bisection precision) native-contact intersections of the same
/// anti-diagonal triangles used by the current 2m Unreal/collision mesh.
/// Rim/toe segments are the actual corners of those faces, not smoothed lines.
/// They are local crossing-face corners. They coincide with the full lip/foot
/// on one-cell collapsed faces (the measured cliff). General multi-cell or
/// pinned faces need a further rim traversal; no global rim coverage is claimed.
#[derive(serde::Serialize)]
pub struct Front {
    pub source_cell:u32,
    pub cap_layer:u16,
    pub triangle:[i64;3],
    pub basal_xyz_m:[[f64;3];2],
    pub rim_xyz_m:[[f64;3];2],
    pub toe_xyz_m:[[f64;3];2],
    pub above_recorded_top:bool,
    pub max_edge_tangent:f64,
}
pub fn fronts(e:&crate::contrast::Evolved,h:&[i32],ox:i64,oy:i64)->Vec<Front> {
    let n=crate::TILE_POSTS;assert_eq!(h.len(),n*n);let mut out=Vec::new();
    for y in 0..n-1 {for x in 0..n-1 {
        let vertices=[y*n+x,y*n+x+1,(y+1)*n+x,(y+1)*n+x+1];
        for (ti,tri) in [[0usize,1,2],[1,3,2]].into_iter().enumerate() {
            let xyz:[[f64;3];3]=std::array::from_fn(|j|{let k=vertices[tri[j]];[(ox+(k%n)as i64)as f64*2.,(oy+(k/n)as i64)as f64*2.,h[k]as f64/1000.]});
            let cols:Vec<_>=xyz.iter().map(|p|e.column(p[0],p[1])).collect();
            if cols.iter().any(|c|c.0.cell!=cols[0].0.cell) {continue;}
            let col=cols[0].0;
            for k in 0..col.layers.len().saturating_sub(1) {
                let cap=&col.layers[k];if cap.code!=6||col.layers[k+1].code!=5 {continue;}
                let phi:[f64;3]=std::array::from_fn(|i|xyz[i][2]-(cols[i].1-cap.depth));
                let hard=phi.iter().filter(|&&v|v>0.).count();if hard==0||hard==3 {continue;}
                let mut basal=Vec::new();let mut rim=Vec::new();let mut toe=Vec::new();let mut slope=0.0_f64;
                for (a,b) in [(0,1),(1,2),(2,0)] {
                    let pa=xyz[a];let pb=xyz[b];let d=(pa[0]-pb[0]).hypot(pa[1]-pb[1]);
                    slope=slope.max((pa[2]-pb[2]).abs()/d);
                    if (phi[a]>0.)==(phi[b]>0.) {continue;}
                    // Canonical endpoint ordering makes shared-edge roots bit
                    // identical across tiles, independent of triangle winding.
                    let (lo,hi)=if (pa[1],pa[0])<(pb[1],pb[0]) {(pa,pb)}else{(pb,pa)};
                    let f=|t:f64| {let p=std::array::from_fn::<_,3,_>(|j|lo[j]+t*(hi[j]-lo[j]));let (c,z)=e.column(p[0],p[1]);assert_eq!(c.cell,col.cell);(p,p[2]-(z-cap.depth))};
                    let mut l=0.;let mut r=1.;let sign=f(l).1>0.;
                    for _ in 0..48 {let m=(l+r)/2.;if (f(m).1>0.)==sign {l=m;} else {r=m;}}
                    basal.push(f((l+r)/2.).0);
                    let (up,dn)=if phi[a]>0. {(pa,pb)}else{(pb,pa)};rim.push(up);toe.push(dn);
                }
                assert_eq!(basal.len(),2,"triangle contact must have exactly two crossings");
                out.push(Front{source_cell:col.cell,cap_layer:cap.index,triangle:[ox+x as i64,oy+y as i64,ti as i64],
                    basal_xyz_m:[basal[0],basal[1]],rim_xyz_m:[rim[0],rim[1]],toe_xyz_m:[toe[0],toe[1]],
                    above_recorded_top:xyz.iter().enumerate().any(|(i,p)|p[2]>cols[i].1-(cap.depth-cap.h)),max_edge_tangent:slope});
            }
        }
    }}
    out
}

/// ESCARP-2 is a conservative projection of the calibrated native solution.
/// Cache keys are global parent posts; neither tile request order nor tile
/// boundaries can change the reconstruction or its exact per-post ledger.
pub struct Reconstructed {
    antialias:bool,
    toes:bool,
    grid:pf_region::t1b_escarp::Reconstruction,
    band:Option<pf_region::t1b_escarp::BandReconstruction>,
    cache:std::sync::Mutex<HashMap<(i64,i64),std::sync::Arc<pf_region::t1b_escarp::ReconstructedCell>>>,
}
impl Reconstructed {
    pub fn from_window(w:&crate::input::T1bWindow)->R<Option<Self>> {
        let mode=w.ct.as_ref().and_then(|c|c.params.get("escarp_reconstruct"));
        if mode.is_none() {return Ok(None);}
        if !matches!(mode.and_then(|v|v.as_u64()),Some(2|3|4|6)) {return Err("unsupported ESCARP reconstruction".into());}
        let e=w.evolved.as_ref().ok_or("ESCARP-2 needs a calibrated process-evolved window")?;
        let mut contacts=Vec::with_capacity(w.n*w.n);
        for y in 0..w.n {for x in 0..w.n {
            let (c,datum)=e.column((w.gx0+x as i64)as f64*32.,(w.gy0+y as i64)as f64*32.);
            let h=w.height_m[y*w.n+x];let mut best:Option<(f64,pf_region::t1b_ct::CapContact)>=None;
            for pair in c.layers.windows(2) {let a=&pair[0];let b=&pair[1];if a.code!=6||b.code!=5 {continue;}
                let base=datum-a.depth;let top=base+a.h;let bottom=datum-b.depth;
                let distance=(h-top).max(bottom-h).max(0.)+(h-base).abs()*1e-9;
                if best.as_ref().is_none_or(|v|distance<v.0) {
                    let q=|v:f64|(v*1000.).round()as i64;
                    best=Some((distance,pf_region::t1b_ct::CapContact{source:(c.cell,a.index),base_mm:q(base),top_mm:q(top),weak_bottom_mm:q(bottom),cap_rate:0,weak_rate:0}));
                }
            }
            contacts.push(best.map(|v|v.1));
        }}
        let h=w.height_m.iter().map(|v|(v*1000.).round()as i64).collect();let fixed=w.fixed.iter().map(|v|*v!=0).collect();
        let grid=pf_region::t1b_escarp::Reconstruction::new(w.n,[w.gx0,w.gy0],h,contacts,fixed)?;
        let toes=mode.and_then(|v|v.as_u64())==Some(6);
        let antialias=matches!(mode.and_then(|v|v.as_u64()),Some(4|6));
        let mut band=match mode.and_then(|v|v.as_u64()) {
            Some(4|6)=>Some(pf_region::t1b_escarp::BandReconstruction::build_antialiased(&grid)?),
            Some(3)=>Some(pf_region::t1b_escarp::BandReconstruction::build(&grid)?),
            _=>None,
        };
        if toes {
            let tangent=e.classes.iter().find(|c|c.code==6&&c.mode!=pf_region::t1b_ct::CT_NONE)
                .ok_or("ESCARP6 missing native sandstone class")?.talus_tan_permille as f64/1000.;
            band.as_mut().unwrap().soften_toes(&grid,tangent)?;
        }
        Ok(Some(Self{antialias,toes,grid,band,cache:std::sync::Mutex::new(HashMap::new())}))
    }
    fn cell(&self,gx:i64,gy:i64)->std::sync::Arc<pf_region::t1b_escarp::ReconstructedCell> {
        if let Some(c)=self.cache.lock().unwrap().get(&(gx,gy)).cloned() {return c;}
        let c=std::sync::Arc::new(if let Some(b)=&self.band {b.cell(&self.grid,gx,gy)}else{self.grid.cell(gx,gy)});self.cache.lock().unwrap().insert((gx,gy),c.clone());c
    }
    pub fn height_mm(&self,px:i64,py:i64)->i32 {
        let gx=(px+8).div_euclid(16);let gy=(py+8).div_euclid(16);
        let u=(px+8).rem_euclid(16)as usize;let v=(py+8).rem_euclid(16)as usize;
        self.cell(gx,gy).heights[v*16+u]
    }
    /// ESCARP-5 opt-in composition. Keep every conserved nearest-post prism
    /// exactly. Outside its union, a C1 partition blends over 64m to bounded
    /// C1 base interpolation. The product of compact rectangle complements
    /// avoids nearest-distance medial-axis creases. Coordinates are global
    /// 2m posts, so adjacent tiles use the identical mask.
    pub fn floor_weight(&self,px:i64,py:i64)->f64 {
        assert!(self.antialias,"floor C1 composition requires ESCARP-4");
        let band=self.band.as_ref().unwrap();
        let (gx,gy)=((px+8).div_euclid(16),(py+8).div_euclid(16));
        if band.cells.contains_key(&(gx,gy)) {return 0.;}
        let bump=|p:i64,g:i64| {
            let distance=((p-g*16).abs()-8).max(0) as f64;
            let t=(distance/32.).min(1.);
            1.-t*t*(3.-2.*t)
        };
        let mut weight=1.;
        for y in gy-3..=gy+3 {for x in gx-3..=gx+3 {
            if band.cells.contains_key(&(x,y)) {
                weight*=1.-bump(px,x)*bump(py,y);
            }
        }}
        weight
    }
    pub fn height_with_floor_mm(&self,px:i64,py:i64,c1:&crate::base_interp::C1Surface)->i32 {
        let weight=self.floor_weight(px,py);
        if weight==0. {return self.height_mm(px,py);}
        let smooth=c1.height_mm(px,py);
        if weight==1. {return smooth;}
        let old=self.height_mm(px,py);
        (old as f64+weight*(smooth as f64-old as f64)).round()as i32
    }
    pub fn rule(&self)->&'static str {
        if self.toes {pf_region::t1b_escarp::TOE_RULE}else if self.antialias {pf_region::t1b_escarp::ANTIALIASED_BAND_RULE}else if self.band.is_some(){pf_region::t1b_escarp::BAND_RECONSTRUCTION_RULE}else{pf_region::t1b_escarp::RECONSTRUCTION_RULE}
    }
    pub fn has_toes(&self)->bool {self.toes}
    pub fn toe_ledger(&self)->Option<serde_json::Value> {
        self.toes.then(||serde_json::json!({"rule":pf_region::t1b_escarp::TOE_RULE,"accounts":self.band.as_ref().unwrap().toe_accounts,
            "account":"incremental signed volume relative to ESCARP4, balanced per source/cap; original band certificate plus these zero-net changes. Do not treat original gross debit/credit as the modified gross distribution."}))
    }
    pub fn is_band(&self)->bool {self.band.is_some()}
    pub fn is_antialiased(&self)->bool {self.antialias}
    pub fn ledger(&self,ox:i64,oy:i64)->serde_json::Value {self.ledger_with_floor(ox,oy,None)}
    pub fn ledger_with_floor(&self,ox:i64,oy:i64,c1:Option<&crate::base_interp::C1Surface>)->serde_json::Value {
        let mut debit=0i128;let mut credit=0i128;let mut max_correction=0;let mut fronts=0;let mut posts=0;
        // Canonical ownership: a tile reports only parent centres in its
        // half-open core. Adjacent tile ledgers never duplicate a post.
        for gy in oy.div_euclid(16)..(oy+128).div_euclid(16) {for gx in ox.div_euclid(16)..(ox+128).div_euclid(16) {
            let c=self.cell(gx,gy);
            if let Some(c1)=c1 {
                let target=self.grid.parent_height_mm(gx,gy);
                for v in 0..16 {for u in 0..16 {
                    let h=self.height_with_floor_mm(gx*16+u-8,gy*16+v-8,c1);
                    let dv=(h as i128-target as i128)*4_000_000;
                    if dv>=0 {credit+=dv;}else{debit-=dv;}
                }}
            }else{debit+=c.debit_mm3;credit+=c.credit_mm3;}
            max_correction=max_correction.max(c.max_correction_mm);fronts+=c.front_samples;posts+=1;
        }}
        if let Some(band)=&self.band {
            return serde_json::json!({"schema":"escarp-reconstruct-calibrated32-band-2m/v1","rule":self.rule(),"new_simulated_time_yr":0,
                "parent_posts":posts,"fine_prisms":posts*256,"front_samples":fronts,
                "debit_mm3":debit.to_string(),"credit_mm3":credit.to_string(),"residual_mm3":(credit-debit).to_string(),
                "max_conservative_correction_m":max_correction as f64/1000.,
                "account":"tile slice, generally nonzero; exact conservation belongs to the complete source bands below, not individual posts or tiles",
                "window_band_origin":self.grid.origin,"window_band_parent_n":self.grid.n,"window_band_accounts":band.accounts,
                "band_account_aggregation":"one certificate per parent window; repeated in each tile for audit, DO NOT sum duplicate window certificates",
                "support":"owner centre in tile half-open core; full nearest-post support [-8,7]; outside contact band C1 cubic reference has no per-post mean constraint"});
        }
        serde_json::json!({"schema":"escarp-reconstruct-calibrated32-2m/v1","rule":pf_region::t1b_escarp::RECONSTRUCTION_RULE,"new_simulated_time_yr":0,"parent_posts":posts,"fine_prisms":posts*256,"front_samples":fronts,"debit_mm3":debit.to_string(),"credit_mm3":credit.to_string(),"residual_mm3":(credit-debit).to_string(),"max_mean_error_mm":0,"max_conservative_correction_m":max_correction as f64/1000.,"account":"within each calibrated coarse post prism; not new process removal or collectible stock","support":"owner centre in tile half-open core, full nearest-post fine support [-8,7] included even across tile boundary"})
    }
}

#[cfg(test)]
mod reconstruction_tests {
    use super::*;
    #[test]
    fn native_calibrated_post_means_and_mode_guards() {
        let Ok(path)=std::env::var("VZ_ESCARP2_NATIVE_WINDOW")else{eprintln!("SKIP: set VZ_ESCARP2_NATIVE_WINDOW to the calibrated native Moab window");return};
        let mut w=crate::input::T1bWindow::load(std::path::Path::new(&path)).unwrap();
        assert!(Reconstructed::from_window(&w).unwrap().is_none(),"key-off unexpectedly enabled reconstruction");
        w.ct.as_mut().unwrap().params["escarp_reconstruct"]=2.into();
        let fine=Reconstructed::from_window(&w).unwrap().unwrap();
        let mut checked=0;let mut debit=0i128;let mut credit=0i128;let mut front=0;
        // Actual requested cliff and three bead tiles, plus one cell of support.
        for (x0,y0) in [(200,264),(240,248),(248,248),(248,256)] {
            for gy in y0..y0+9 {for gx in x0..x0+9 {
                let c=fine.cell(gx,gy);let sum:i64=c.heights.iter().map(|v|*v as i64).sum();
                let parent=(w.height_m[w.idx(gx,gy).unwrap()]*1000.).round()as i64;
                assert_eq!(sum,256*parent,"post mean changed at {gx},{gy}");
                assert_eq!(c.debit_mm3,c.credit_mm3);
                debit+=c.debit_mm3;credit+=c.credit_mm3;front+=c.front_samples;checked+=1;
            }}
        }
        assert!(front>0,"vacuous contact-front reconstruction");
        w.ct.as_mut().unwrap().params["escarp_reconstruct"]=99.into();assert!(Reconstructed::from_window(&w).is_err());
        w.ct.as_mut().unwrap().params["escarp_reconstruct"]=2.into();w.evolved=None;assert!(Reconstructed::from_window(&w).is_err());
        let report=serde_json::json!({"input":path,"source_cell":60967,"parent_post_means_checked":checked,"fine_prisms_checked":checked*256,"front_samples":front,"debit_mm3":debit.to_string(),"credit_mm3":credit.to_string(),"mean_residual_mm":0,"ledger_residual_mm3":(credit-debit).to_string(),"shape_acceptance":"not asserted by accounting test; measured on actual T2 outputs separately"});
        println!("{}",serde_json::to_string_pretty(&report).unwrap());
        if let Ok(out)=std::env::var("VZ_ESCARP2_ACCOUNT_REPORT") {std::fs::write(out,serde_json::to_vec_pretty(&report).unwrap()).unwrap();}
    }
}

#[cfg(test)]
mod band_tests {
    use super::*;
    #[test]
    fn native_band_conserves_by_source_without_per_post_constraint() {native_band(3);}
    #[test]
    fn native_antialiased_band_conserves_by_source() {native_band(4);}
    fn native_band(mode:u32) {
        let key=format!("VZ_ESCARP{mode}_NATIVE_WINDOW");
        let Ok(path)=std::env::var(&key)else{eprintln!("SKIP: set {key} for the native source-band test");return};
        let mut w=crate::input::T1bWindow::load(std::path::Path::new(&path)).unwrap();
        assert!(Reconstructed::from_window(&w).unwrap().is_none());
        w.ct.as_mut().unwrap().params["escarp_reconstruct"]=mode.into();
        let f=Reconstructed::from_window(&w).unwrap().unwrap();let b=f.band.as_ref().unwrap();
        assert_eq!(f.is_antialiased(),mode==4);
        assert!(!b.accounts.is_empty());assert!(!b.cells.is_empty());
        let mut totals=std::collections::BTreeMap::<u32,(i128,i128)>::new();let mut changed_means=0;
        for (&(gx,gy),c) in &b.cells {
            let target=(w.height_m[w.idx(gx,gy).unwrap()]*1000.).round()as i64;
            let sum:i64=c.heights.iter().map(|v|*v as i64).sum();changed_means+=usize::from(sum!=target*256);
            let (col,_)=w.evolved.as_ref().unwrap().column(gx as f64*32.,gy as f64*32.);
            let t=totals.entry(col.cell).or_default();t.0+=(sum-target*256)as i128;t.1+=c.credit_mm3-c.debit_mm3;
            assert_eq!((sum-target*256)as i128*4_000_000,c.credit_mm3-c.debit_mm3);
        }
        assert!(changed_means>0,"per-post constraint was accidentally retained");
        for (_, (sum,volume)) in totals {assert_eq!(sum,0);assert_eq!(volume,0);}
        for a in &b.accounts {assert_eq!(a.debit_mm3,a.credit_mm3);assert_eq!(a.residual_mm3,"0");assert!(a.rms_post_mean_deviation_m>0.);}
        // Adjacent tile requests see exactly the same global post, independent
        // of query order and cache population. This is NOT a cross-window test.
        for py in 4224..4353 {let a=f.height_mm(4096,py);let _=f.height_mm(4097,py);assert_eq!(a,f.height_mm(4096,py));}
        let report=serde_json::json!({"input":path,"band_parent_posts":b.cells.len(),"changed_post_means":changed_means,"accounts":b.accounts,"ledger_residual_mm3":"0"});
        println!("{}",serde_json::to_string_pretty(&report).unwrap());
        if let Ok(out)=std::env::var(format!("VZ_ESCARP{mode}_ACCOUNT_REPORT")){std::fs::write(out,serde_json::to_vec_pretty(&report).unwrap()).unwrap();}
    }
}

#[cfg(test)]
mod floor_tests {
    use super::*;
    #[test]
    fn native_c1_floor_keeps_all_band_prisms_and_parent_posts() {
        let Ok(path)=std::env::var("VZ_ESCARP5_NATIVE_WINDOW")else{
            eprintln!("SKIP: set VZ_ESCARP5_NATIVE_WINDOW for actual ESCARP4+C1 composition");return
        };
        let w=crate::input::T1bWindow::load(std::path::Path::new(&path)).unwrap();
        let f=Reconstructed::from_window(&w).unwrap().unwrap();assert!(f.is_antialiased());
        let c1=crate::base_interp::C1Surface::new(&w).unwrap();let b=f.band.as_ref().unwrap();
        let mut checked=0;let mut volume=0i128;
        for (&(gx,gy),cell) in &b.cells {
            let target=f.grid.parent_height_mm(gx,gy);
            for v in 0..16 {for u in 0..16 {
                let (px,py)=(gx*16+u-8,gy*16+v-8);
                assert_eq!(f.floor_weight(px,py),0.);
                let h=f.height_with_floor_mm(px,py,&c1);
                assert_eq!(h,cell.heights[(v*16+u)as usize]);
                volume+=(h as i128-target as i128)*4_000_000;checked+=1;
            }}
        }
        assert!(checked>0);assert_eq!(volume,0);
        let mut exact_posts=0;let mut fixed_posts=0;let mut changed=0;let mut pure=0;
        for (tx,ty) in [(28,30),(27,29),(29,30)] {
            for gy in ty*8..=ty*8+8 {for gx in tx*8..=tx*8+8 {
                if f.floor_weight(gx*16,gy*16)==1. {
                    assert_eq!(f.height_with_floor_mm(gx*16,gy*16,&c1)as i64,f.grid.parent_height_mm(gx,gy));
                    exact_posts+=1;fixed_posts+=usize::from(w.fixed[w.idx(gx,gy).unwrap()]!=0);
                }
            }}
            for py in ty*128..=ty*128+128 {for px in tx*128..=tx*128+128 {
                let weight=f.floor_weight(px,py);assert!((0.0..=1.0).contains(&weight));
                if weight==1. {
                    let h=f.height_with_floor_mm(px,py,&c1);assert_eq!(h,c1.height_mm(px,py));pure+=1;
                    changed+=usize::from(h!=f.height_mm(px,py));
                }
            }}
        }
        assert!(fixed_posts>0&&exact_posts>0&&changed>0&&pure>0);
        let report=serde_json::json!({"input":path,"band_prisms_unchanged":checked,"band_residual_mm3":volume.to_string(),"floor_exact_posts":exact_posts,"fixed_exact_posts":fixed_posts,"pure_c1_samples":pure,"changed_floor_samples":changed});
        println!("{}",report);
        if let Ok(out)=std::env::var("VZ_ESCARP5_ACCOUNT_REPORT") {std::fs::write(out,serde_json::to_vec_pretty(&report).unwrap()).unwrap();}
    }
}

#[cfg(test)]
mod toe_tests {
    use super::*;
    #[test]
    fn native_toe_increment_conserves_and_keeps_resolved_cliff() {
        let Ok(path)=std::env::var("VZ_ESCARP6_NATIVE_WINDOW")else{eprintln!("SKIP: set VZ_ESCARP6_NATIVE_WINDOW to the actual mode4 Moab parent");return};
        let mut w=crate::input::T1bWindow::load(std::path::Path::new(&path)).unwrap();
        assert_eq!(w.ct.as_ref().unwrap().params["escarp_reconstruct"],4);
        let before=Reconstructed::from_window(&w).unwrap().unwrap();
        let c1=crate::base_interp::C1Surface::new(&w).unwrap();
        w.ct.as_mut().unwrap().params["escarp_reconstruct"]=6.into();
        let after=Reconstructed::from_window(&w).unwrap().unwrap();
        let (a,b)=(before.band.as_ref().unwrap(),after.band.as_ref().unwrap());
        assert_eq!(a.cells.len(),b.cells.len());
        let mut changed=0usize;let mut delta=0i128;let mut fixed=0usize;
        for (&xy,cell) in &a.cells {
            for k in 0..256 {let d=b.cells[&xy].heights[k]as i128-cell.heights[k]as i128;delta+=d;changed+=usize::from(d!=0);}
        }
        assert!(changed>0);assert_eq!(delta,0);
        for a in &b.toe_accounts {assert_eq!(a.debit_mm3,a.credit_mm3);assert_eq!(a.residual_mm3,"0");}
        let mut cliff=0;
        for py in 4269..=4291 {for px in 3200..=3328 {
            assert_eq!(before.height_with_floor_mm(px,py,&c1),after.height_with_floor_mm(px,py,&c1),"main cliff changed at {},{}",px*2,py*2);cliff+=1;
        }}
        // Real window fixed locations are absent from both contact bands.
        for gy in 200..280 {for gx in 190..250 {
            if w.fixed[w.idx(gx,gy).unwrap()]!=0 {
                let (px,py)=(gx*16,gy*16);
                assert_eq!(before.height_with_floor_mm(px,py,&c1),after.height_with_floor_mm(px,py,&c1));fixed+=1;
            }
        }}
        assert!(fixed>0);
        // A real-parent parameter overlay for served-path tests. This does
        // not pretend to be a new regional solve or silently reuse its digest.
        if let Ok(dest)=std::env::var("VZ_ESCARP6_PARENT_OUTPUT") {
            let original=std::path::Path::new(&path);
            let mut doc:serde_json::Value=serde_json::from_slice(&std::fs::read(original.join("window.json")).unwrap()).unwrap();
            let source_digest=doc["digest"].clone();
            let h=&mut doc["header"];
            let digest=|v:&serde_json::Value|pf_contracts::artifact::digest_canonical(v).unwrap();
            assert_eq!(h["recipe"]["params_digest"],digest(&h["recipe"]["params"]));
            assert_eq!(h["recipe_digest"],digest(&h["recipe"]));
            h["recipe"]["params"]["rivernet"]["ls"]["tb"]["tc"]["ct"]["escarp_reconstruct"]=6.into();
            h["recipe"]["params_digest"]=digest(&h["recipe"]["params"]).into();
            h["recipe_digest"]=digest(&h["recipe"]).into();
            h["refinement_test_parent"]=serde_json::json!({"source_header_digest":source_digest,"source_directory":path,
                "kind":"native mode4 parent fields retained exactly; test-only mode6 refinement selector; no regional re-solve claimed"});
            doc["digest"]=digest(&doc["header"]).into();
            let dest=std::path::Path::new(&dest);assert!(!dest.exists(),"refuse to overwrite native test carrier");
            std::fs::create_dir_all(dest).unwrap();
            std::os::unix::fs::symlink(original.join("fields"),dest.join("fields")).unwrap();
            std::fs::write(dest.join("window.json"),serde_json::to_vec_pretty(&doc).unwrap()).unwrap();
            let loaded=crate::input::T1bWindow::load(dest).unwrap();
            assert_eq!(loaded.ct.as_ref().unwrap().params["escarp_reconstruct"],6);
            assert_eq!(loaded.height_m,w.height_m);assert_eq!(loaded.fixed,w.fixed);
        }
        let out=std::env::var("VZ_ESCARP6_TEST_OUTPUT").ok();
        if let Some(out)=out {
            let out=std::path::Path::new(&out);std::fs::create_dir_all(out).unwrap();
            for (name,f) in [("before",&before),("after",&after)] {
                for (tx,ty) in [(29,27),(29,30),(25,33)] {
                    let mut bytes=Vec::new();
                    for py in ty*128..=ty*128+128 {for px in tx*128..=tx*128+128 {bytes.extend(f.height_with_floor_mm(px,py,&c1).to_le_bytes());}}
                    std::fs::write(out.join(format!("{name}_{tx}_{ty}.i32")),bytes).unwrap();
                }
            }
            let report=serde_json::json!({"input":path,"band_parent_posts":a.cells.len(),"compared_fine_prisms":a.cells.len()*256,"changed_fine_prisms":changed,"incremental_residual_mm3":(delta*4_000_000).to_string(),"cliff_samples_exact":cliff,"fixed_samples_exact":fixed,"toe":after.toe_ledger()});
            std::fs::write(out.join("native_accounts.json"),serde_json::to_vec_pretty(&report).unwrap()).unwrap();
            println!("{}",report);
        }
    }
}

#[cfg(test)]
mod toe_served_tests {
    use super::*;
    #[test]
    fn served_toe_uses_exact_native_geological_columns() {
        let Ok(path)=std::env::var("VZ_ESCARP6_NATIVE_WINDOW")else{eprintln!("SKIP: set VZ_ESCARP6_NATIVE_WINDOW");return};
        let Ok(out)=std::env::var("VZ_ESCARP6_SERVED_OUTPUT")else{eprintln!("SKIP: set VZ_ESCARP6_SERVED_OUTPUT to actual pf-serve mode6 export");return};
        let w=crate::input::T1bWindow::load(std::path::Path::new(&path)).unwrap();let e=w.evolved.as_ref().unwrap();
        let mut posts=0;let mut ids=std::collections::BTreeSet::new();
        for (tx,ty) in [(29,27),(29,30),(25,33),(28,30)] {
            let dir=std::path::Path::new(&out).join(format!("tiles/t_{tx}_{ty}"));
            let read=|f:&str|std::fs::read(dir.join(f)).unwrap();
            let h=read("height_mm.i32");let pv=read("prov_id.u32");let layer=read("layer.u8");let litho=read("litho.u8");
            let above=read("height_above_unit_base_m.f32");let depth=read("depth_in_unit_m.f32");
            for j in 0..129 {for i in 0..129 {
                let k=j*129+i;let z=i32::from_le_bytes(h[k*4..k*4+4].try_into().unwrap())as f64/1000.;
                let (c,datum)=e.column((tx*256+i*2)as f64,(ty*256+j*2)as f64);
                let d=datum-z;let l=c.layers.iter().find(|l|d<l.depth).unwrap();
                let id=u32::from_le_bytes(pv[k*4..k*4+4].try_into().unwrap());
                assert_eq!(id,crate::strata_prov::strata_id(c.cell,l.index));
                assert_eq!(layer[k],l.index as u8);assert_eq!(litho[k],l.code);
                assert_eq!(f32::from_le_bytes(above[k*4..k*4+4].try_into().unwrap()),(l.depth-d)as f32);
                assert_eq!(f32::from_le_bytes(depth[k*4..k*4+4].try_into().unwrap()),(d-(l.depth-l.h))as f32);
                ids.insert(id);posts+=1;
            }}
        }
        assert!(ids.len()>1,"vacuous contact check");
        let report=serde_json::json!({"parent":path,"actual_served_export":out,"posts":posts,"height_identity_layer_and_depth_mismatches":0,"provenance_ids":ids});
        println!("{}",report);
        if let Ok(p)=std::env::var("VZ_ESCARP6_PROVENANCE_REPORT") {std::fs::write(p,serde_json::to_vec_pretty(&report).unwrap()).unwrap();}
    }
}
