//! ESCARP-1 candidate. Sparse 2 m execution of E5, with Jacobi updates and
//! the same finite interval, contacts, pinning and explicit export reservoir.
//! Coarse samples restrict removal; the pre-process surface plus fine removals
//! is the authoritative detailed geometry. This is not smoothing or extra time.
use crate::{R, t1b_ct::{CapContact, ContrastParams, ContrastStats, ContrastExportAccount}};
use serde::Serialize;
use std::collections::{HashMap, HashSet, BTreeMap};
pub const EVOLVED_SCHEMA: &str = "e5b-evolved-surface/escarp-v1";
const Q: usize = 16;
const AREA_MM2: i128 = 4_000_000;

#[derive(Clone, Debug, Serialize)]
pub struct FineSurface {
    pub schema: String,
    pub spacing_mm: i64,
    pub evaluation_origin_parent: [i64;2],
    pub evaluation_parent_n: usize,
    pub dt_yr: i64,
    pub steps: u32,
    pub changed_fine_posts: usize,
    pub published_changed_fine_posts: usize,
    pub peak_active_posts: usize,
    pub coarse_restriction_rounding_m3: f64,
    pub rule: String,
    #[serde(skip)] pub base_mm: Vec<i32>,
    /// Global 2 m lattice x,y and signed removal in mm. Canonical row order.
    #[serde(skip)] pub delta_xyz: Vec<[i32;3]>,
}
struct Grid<'a> {
    n: usize, nf: usize, h: &'a [i64], c: &'a [Option<CapContact>], pin: &'a [bool],
}
impl Grid<'_> {
    fn interp(&self, k: u64, f: impl Fn(usize)->i64) -> i64 {
        let x=k as usize%self.nf; let y=k as usize/self.nf;
        let gx=(x/Q).min(self.n-2);let gy=(y/Q).min(self.n-2);
        let u=(x-gx*Q) as i64;let v=(y-gy*Q) as i64;let q=Q as i64;
        let i=gy*self.n+gx;
        let a=(f(i)*(q-u)+f(i+1)*u)*(q-v)+(f(i+self.n)*(q-u)+f(i+self.n+1)*u)*v;
        if a>=0 {(a+128)/256}else{-((-a+128)/256)}
    }
    fn height(&self,k:u64,d:&HashMap<u64,i64>)->i64 {self.interp(k,|i|self.h[i])+d.get(&k).copied().unwrap_or(0)}
    fn contact(&self,k:u64)->Option<CapContact> {
        let x=k as usize%self.nf;let y=k as usize/self.nf;
        if x==0||y==0||x+1==self.nf||y+1==self.nf{return None;}
        let i=(y/Q)*self.n+x/Q;
        let c=self.c[i]?;
        for j in [i+1,i+self.n,i+self.n+1] {if self.c[j]?.source!=c.source{return None;}}
        // Interpolation is only inside a single real source contact. No
        // renumbering or interpolation through absent/different columns.
        let out=CapContact {source:c.source,
            base_mm:self.interp(k,|j|self.c[j].unwrap().base_mm),
            top_mm:self.interp(k,|j|self.c[j].unwrap().top_mm),
            weak_bottom_mm:self.interp(k,|j|self.c[j].unwrap().weak_bottom_mm),
            cap_rate:c.cap_rate,weak_rate:c.weak_rate};
        let h0=self.interp(k,|j|self.h[j]);
        if h0>out.top_mm||h0<=out.weak_bottom_mm {None} else {Some(out)}
    }
    fn pinned(&self,k:u64)->bool {
        let x=k as usize%self.nf;let y=k as usize/self.nf;
        x==0||y==0||x+1==self.nf||y+1==self.nf||self.pin[((y+8)/16)*self.n+(x+8)/16]
    }
    fn neighbours(&self,k:u64)->[u64;4] {[k-1,k+1,k-self.nf as u64,k+self.nf as u64]}
}

/// The returned coarse removals are sums of 256 equal-area fine post prisms,
/// NOT rounded heights. Keep their fractional mm until coarse publication.
pub fn run(n:usize,h0:&[i64],contacts:&[Option<CapContact>],pinned:&[bool],p:&ContrastParams,
    global_origin:[i64;2],pub_box:(usize,usize,usize,usize)) -> R<(Vec<i64>,ContrastStats)> {
    p.validate()?;
    if n<3||h0.len()!=n*n||contacts.len()!=n*n||pinned.len()!=n*n||p.dt_yr%1000!=0 {
        return Err("ESCARP-1 invalid lattice or interval not divisible by 1000yr".into());
    }
    let dt=1000i64;let steps=p.steps.checked_mul((p.dt_yr/dt)as u32).ok_or("ESCARP-1 step overflow")?;
    let retreat=p.retreat_mm_per_kyr*dt/1000;
    if retreat>2000{return Err("ESCARP-1 retreat exceeds 2m; reduce the fine time step".into());}
    let nf=(n-1).checked_mul(Q).and_then(|x|x.checked_add(1)).ok_or("ESCARP-1 extent overflow")?;
    let g=Grid{n,nf,h:h0,c:contacts,pin:pinned};
    let mut active=HashSet::new();
    // Search coarse mixed cells, then test the actual fine interface. The law
    // acts ONLY next to a matching hard/weak contact. Everything else is static.
    // One fine-node ring includes nodes on the other side of a coarse boundary.
    for y in 0..n-1 {for x in 0..n-1 {
        let ids=[y*n+x,y*n+x+1,(y+1)*n+x,(y+1)*n+x+1];
        let mut hard=false;let mut weak=false;
        for i in ids {if let Some(c)=contacts[i] {if h0[i]>c.base_mm {hard=true}else{weak=true}}}
        if !hard||!weak {continue;}
        for fy in (y*Q).saturating_sub(1).max(1)..=((y+1)*Q+1).min(nf-2) {
            for fx in (x*Q).saturating_sub(1).max(1)..=((x+1)*Q+1).min(nf-2) {active.insert((fy*nf+fx)as u64);}
        }
    }}
    let mut delta=HashMap::<u64,i64>::new();let mut debt=HashMap::<u64,i64>::new();
    let mut accounts=BTreeMap::<(u32,u16),[i64;3]>::new();
    let mut st=ContrastStats{conditional:true,contact_posts:contacts.iter().filter(|v|v.is_some()).count()as u64,..Default::default()};
    let mut peak=active.len();
    for step in 0..steps {
        let mut changes=Vec::new();let mut next_active=HashSet::with_capacity(active.len());
        for &k in &active {
            if g.pinned(k) {continue;}
            let Some(c)=g.contact(k)else{continue};let h=g.height(k,&delta);
            let nb=g.neighbours(k);let low=nb.iter().map(|&j|g.height(j,&delta)).min().unwrap();
            let floor=(low+1).max(c.weak_bottom_mm+1);
            if h<=floor {debt.remove(&k);continue;}
            let hard=h>c.base_mm;let mut toe=false;let mut unsupported=false;
            for j in nb {let Some(cj)=g.contact(j)else{continue};if cj.source!=c.source{continue;}
                let hj=g.height(j,&delta);
                if hard&&hj<=cj.base_mm {toe=true;if cj.base_mm-hj>p.support_tan_permille*2000/1000 {unsupported=true;}}
                if !hard&&hj>cj.base_mm {toe=true;}
            }
            if !toe {debt.remove(&k);continue;}
            let a=accounts.entry(c.source).or_default();let mut dh=0;
            if hard {
                let w=(c.cap_rate*dt/1000).min(h-c.base_mm).min(h-floor);dh-=w;a[0]+=w;
                if unsupported {
                    st.unsupported_post_steps+=1;let d=debt.entry(k).or_default();*d+=retreat;
                    if *d>=2000 {
                        let f=(h+dh-c.base_mm.max(floor)).max(0);dh-=f;a[2]+=f;
                        if f>0 {st.failed_slabs+=1;}*d=0;
                    }
                } else {debt.remove(&k);}
            } else {let w=(c.weak_rate*dt/1000).min(h-floor);dh-=w;a[1]+=w;}
            if dh!=0 {changes.push((k,dh));}
            // Debt may advance without a height change. A still exposed
            // contact must remain active; newly exposed neighbours enter below.
            next_active.insert(k);
        }
        // Jacobi: not one height changes until every old-state flux is read.
        for (k,dh) in changes {
            *delta.entry(k).or_default()+=dh;
            for j in g.neighbours(k) {next_active.insert(j);}
        }
        active=next_active;peak=peak.max(active.len());
        if step%24==0 {println!("  escarp fine step {step}/{steps}, active {}, removed {}",active.len(),delta.len());}
    }
    let mut coarse=vec![0i64;n*n];let mut entries:Vec<_>=delta.into_iter().collect();entries.sort_unstable_by_key(|v|v.0);
    let mut surface=FineSurface{schema:"escarp-sparse-e5-2m/v1".into(),spacing_mm:2000,evaluation_origin_parent:global_origin,evaluation_parent_n:n,dt_yr:dt,steps,
        changed_fine_posts:entries.len(),published_changed_fine_posts:0,peak_active_posts:peak,coarse_restriction_rounding_m3:0.,
        rule:"Candidate: E5 at 2m, 1000yr Jacobi steps for the same total interval; all removal credited to the E5 external reservoir; 32m means restrict 16x16 fine post-prism removals. Initial geometry is bilinear pre-E5 32m terrain. Contact planes interpolate only the same retained source unit; fixed masks remain nearest native-post masks. No talus, joints, new geology, mean compensation or time duplication. Four-direction law still has possible directional bias at 2m; not a convergence claim.".into(),base_mm:Vec::new(),delta_xyz:Vec::new()};
    let mut sum_delta=0i64;
    for (k,dh) in entries {
        let x=k as usize%nf;let y=k as usize/nf;let owner=((y+8)/16)*n+(x+8)/16;
        coarse[owner]+=dh;sum_delta+=dh;
        if x>=pub_box.0*Q&&x<=pub_box.1*Q&&y>=pub_box.2*Q&&y<=pub_box.3*Q {
            let gx=global_origin[0]*16+x as i64;let gy=global_origin[1]*16+y as i64;
            surface.delta_xyz.push([i32::try_from(gx).map_err(|_|"fine x overflow")?,i32::try_from(gy).map_err(|_|"fine y overflow")?,i32::try_from(dh).map_err(|_|"fine removal overflow")?]);
        }
    }
    surface.published_changed_fine_posts=surface.delta_xyz.len();
    for y in pub_box.2..=pub_box.3 {for x in pub_box.0..=pub_box.1 {surface.base_mm.push(i32::try_from(h0[y*n+x]).map_err(|_|"base overflow")?);}}
    let mut sums=[0i64;3];
    for ((cell,layer),v) in accounts {if v==[0;3]{continue;}
        for k in 0..3 {sums[k]+=v[k];}
        st.export_accounts.push(ContrastExportAccount{source_cell:cell,cap_layer:layer,
            cap_weathering_mm3:(v[0]as i128*AREA_MM2).to_string(),weak_weathering_mm3:(v[1]as i128*AREA_MM2).to_string(),collapsed_cap_mm3:(v[2]as i128*AREA_MM2).to_string()});
    }
    st.local_ledger_residual_mm=sum_delta+sums.iter().sum::<i64>();
    if st.local_ledger_residual_mm!=0||coarse.iter().sum::<i64>()!=sum_delta {return Err("ESCARP-1 fine-to-coarse material ledger failed".into());}
    st.cap_weathering_m3=sums[0]as f64*0.004;st.weak_weathering_m3=sums[1]as f64*0.004;st.collapse_export_m3=sums[2]as f64*0.004;
    st.mean_shift_m=sum_delta as f64/256./1000./(n*n)as f64;
    st.export_reservoir="e5_unresolved_external_sediment: fine removed in-situ volume, efficient evacuation conditional; not collectible stock".into();
    st.subgrid=Some(surface);
    Ok((coarse,st))
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Real-data differential certificate: the optimized sparse solver must
    /// equal the existing full Jacobi kernel, with the very same 2m inputs.
    #[test]
    fn native_sparse_equals_dense_and_restricts_every_removed_mm() {
        let Ok(path)=std::env::var("VZ_ESCARP_NATIVE_WINDOW")else{eprintln!("SKIP: set VZ_ESCARP_NATIVE_WINDOW to a published E5B native window");return};
        let dir=std::path::Path::new(&path);
        let doc:serde_json::Value=serde_json::from_slice(&std::fs::read(dir.join("window.json")).unwrap()).unwrap();
        let header=&doc["header"];let lat=&header["lattice"];
        let wn=lat["n"].as_u64().unwrap()as usize;let gx0=lat["gx0"].as_i64().unwrap();let gy0=lat["gy0"].as_i64().unwrap();
        let read=|name:&str| {let b=std::fs::read(dir.join(format!("fields/{name}.i32"))).unwrap();b.chunks_exact(4).map(|v|i32::from_le_bytes(v.try_into().unwrap())).collect::<Vec<_>>()};
        let heights=read("height_mm");let bases=read("unit_base_mm");let tops=read("unit_top_mm");let layers=std::fs::read(dir.join("fields/layer.u8")).unwrap();
        let fixed=std::fs::read(dir.join("fields/fixed.u8")).unwrap();
        let col=header["prep"]["net"]["evolved_surface"]["columns"].as_array().unwrap().iter().find(|c|c["cell"]==60967).unwrap();
        assert_eq!(col["layers"][0]["code"],6);assert_eq!(col["layers"][1]["code"],5);
        let cap=(col["layers"][0]["h"].as_f64().unwrap()*1000.).round()as i64;
        let weak=(col["layers"][1]["h"].as_f64().unwrap()*1000.).round()as i64;
        let (n,ox,oy)=(17usize,197i64,259i64);let mut h0=Vec::new();let mut c=Vec::new();let mut pin=Vec::new();
        for y in 0..n {for x in 0..n {
            let j=(oy+y as i64-gy0)as usize*wn+(ox+x as i64-gx0)as usize;let h=heights[j]as i64;h0.push(h);
            let base=match layers[j] {0=>Some(bases[j]as i64),1=>Some(tops[j]as i64),_=>None};
            c.push(base.filter(|b|h>*b-weak&&h<=*b+cap).map(|base|CapContact{source:(60967,0),base_mm:base,top_mm:base+cap,weak_bottom_mm:base-weak,cap_rate:10,weak_rate:500}));
            pin.push(x==0||y==0||x+1==n||y+1==n||fixed[j]!=0);
        }}
        let nf=(n-1)*16+1;let g=Grid{n,nf,h:&h0,c:&c,pin:&pin};let none=HashMap::new();
        let fine_c:Vec<_>=(0..nf*nf).map(|k|g.contact(k as u64)).collect();
        let fine_p:Vec<_>=(0..nf*nf).map(|k|g.pinned(k as u64)).collect();
        let mut dense:Vec<_>=(0..nf*nf).map(|k|g.height(k as u64,&none)).collect();let initial=dense.clone();
        let mut p=ContrastParams::v1();p.steps=96;p.dt_yr=1000;p.subgrid_spacing_mm=Some(2000);
        let (coarse,stats)=run(n,&h0,&c,&pin,&p,[ox,oy],(0,n-1,0,n-1)).unwrap();
        let reference=crate::t1b_ct::contrast_kernel(nf,2000,&fine_p,&fine_c,&p,&mut dense).unwrap();
        let mut sparse=initial.clone();let surface=stats.subgrid.as_ref().unwrap();
        for &[x,y,d] in &surface.delta_xyz {let i=(y as i64-oy*16)as usize*nf+(x as i64-ox*16)as usize;sparse[i]+=d as i64;}
        assert_eq!(sparse,dense,"sparse scheduler changed the native E5 law");
        let removed:i64=dense.iter().zip(&initial).map(|(a,b)|a-b).sum();
        assert!(removed<0&&surface.changed_fine_posts>0,"vacuous native contact test");
        assert_eq!(coarse.iter().sum::<i64>(),removed);
        assert_eq!(stats.local_ledger_residual_mm,0);
        assert_eq!(stats.cap_weathering_m3,reference.cap_weathering_m3);
        assert_eq!(stats.weak_weathering_m3,reference.weak_weathering_m3);
        assert_eq!(stats.collapse_export_m3,reference.collapse_export_m3);
        let absent=vec![None;n*n];
        let (zero,zst)=run(n,&h0,&absent,&pin,&p,[ox,oy],(0,n-1,0,n-1)).unwrap();
        assert!(zero.iter().all(|v|*v==0)&&zst.export_accounts.is_empty()&&zst.subgrid.as_ref().unwrap().delta_xyz.is_empty(),"no-contact falsifier changed real terrain");
        let report=serde_json::json!({"no_contact_falsifier_exact_fine_posts":nf*nf,"input":path,"coarse_crop":[ox,oy,n as i64],"fine_posts_compared":nf*nf,"steps":96,"dt_yr":1000,"removed_fine_mm":-removed,"changed_fine_posts":surface.changed_fine_posts,"fine_height_mismatches":0,"restriction_residual_mm":0,"ledger_residual_mm":0,"export_m3":stats.cap_weathering_m3+stats.weak_weathering_m3+stats.collapse_export_m3});
        println!("{}",serde_json::to_string_pretty(&report).unwrap());
        if let Ok(out)=std::env::var("VZ_ESCARP_DIFFERENTIAL_REPORT") {std::fs::write(out,serde_json::to_vec_pretty(&report).unwrap()).unwrap();}
    }
}

/// Explicit diagnostic capture of real native inputs, before any fine update.
/// No window is published by this path. Source components cannot exchange
/// evolving neighbours: four matching coarse IDs are required for each fine
/// contact, leaving an unsupported gap between different source components.
pub fn capture_native(inp:&crate::t1b_hs::HsIn,h0:&[i64],c:&[Option<CapContact>],pin:&[bool],p:&ContrastParams,path:&std::path::Path,target:u32)->R<()> {
    use std::io::{Write,BufWriter};
    let crate::t1b_rivernet::RockPosts::Strata{cell,table,..}=inp.rock else{return Err("capture needs native columns".into());};
    let env=table.envelope.as_ref().ok_or("capture needs E3 envelope")?;
    let data=env.load()?;let source=data.cells.get(&target).ok_or("capture target absent")?;
    let n=inp.lat.n;let mut x0=n;let mut y0=n;let mut x1=0;let mut y1=0;let mut count=0;
    for i in 0..cell.len(){if cell[i]==target{let x=i%n;let y=i/n;x0=x0.min(x);y0=y0.min(y);x1=x1.max(x);y1=y1.max(y);count+=1;}}
    if count==0||x0<2||y0<2||x1+2>=n||y1+2>=n{return Err("capture source cell not enclosed".into());}
    let sn=(x1-x0).max(y1-y0)+5;x0=(x0-2).min(n-sn);y0=(y0-2).min(n-sn);
    std::fs::create_dir(path).map_err(|e|e.to_string())?;
    let writer=|name:&str|->R<BufWriter<std::fs::File>> {Ok(BufWriter::new(std::fs::File::create(path.join(name)).map_err(|e|e.to_string())?))};
    let mut wh=writer("height_mm.i32")?;let mut wp=writer("pinned.u8")?;let mut wl=writer("cap_layer.u16")?;
    let mut wb=writer("base_mm.i32")?;let mut wt=writer("top_mm.i32")?;let mut ww=writer("weak_bottom_mm.i32")?;
    let mut rates=None;let mut contact_count=0;
    for y in y0..y0+sn {for x in x0..x0+sn {let i=y*n+x;
        wh.write_all(&i32::try_from(h0[i]).map_err(|_|"capture height overflow")?.to_le_bytes()).map_err(|e|e.to_string())?;
        wp.write_all(&[pin[i]as u8]).map_err(|e|e.to_string())?;
        let cc=c[i].filter(|v|v.source.0==target);
        if let Some(v)=cc {if v.source.1==u16::MAX{return Err("capture sentinel collision".into());}let r=(v.cap_rate,v.weak_rate);if rates.is_some_and(|q|q!=r){return Err("capture rates vary".into());}rates=Some(r);contact_count+=1;}
        wl.write_all(&cc.map(|v|v.source.1).unwrap_or(u16::MAX).to_le_bytes()).map_err(|e|e.to_string())?;
        for (f,v) in [(&mut wb,cc.map(|v|v.base_mm)),(&mut wt,cc.map(|v|v.top_mm)),(&mut ww,cc.map(|v|v.weak_bottom_mm))] {f.write_all(&i32::try_from(v.unwrap_or(0)).map_err(|_|"capture contact overflow")?.to_le_bytes()).map_err(|e|e.to_string())?;}
    }}
    for f in [&mut wh,&mut wp,&mut wl,&mut wb,&mut wt,&mut ww] {f.flush().map_err(|e|e.to_string())?;}
    let meta=serde_json::json!({"schema":"escarp-native-source-input/v1","source_cell":target,"n":sn,"origin":[inp.lat.gx0+x0 as i64,inp.lat.gy0+y0 as i64],"native_source_posts":count,"contact_posts":contact_count,"rates":rates,"params":p,"envelope_digest":env.digest,"source_area_m2":source.area_m2,"source_mean_m":source.mean_m,"source_shape":source.shape,"source_radius_m":source.radius_m,"original_eval_n":n,"capture_kind":"actual native pre-process heights, pinned mask and contact geometry; other source contacts are absent, their initial heights remain as immutable floor neighbours"});
    std::fs::write(path.join("input.json"),serde_json::to_vec_pretty(&meta).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod whole_native {
    use super::*;
    #[test]
    fn captured_full_source_cell_fine_restriction() {
        let Ok(path)=std::env::var("VZ_ESCARP_CAPTURE_TEST")else{eprintln!("SKIP: set VZ_ESCARP_CAPTURE_TEST for whole-source native execution");return};
        let dir=std::path::Path::new(&path);let meta:serde_json::Value=serde_json::from_slice(&std::fs::read(dir.join("input.json")).unwrap()).unwrap();
        assert_eq!(meta["schema"],"escarp-native-source-input/v1");
        let n=meta["n"].as_u64().unwrap()as usize;let cell=meta["source_cell"].as_u64().unwrap()as u32;
        let origin=[meta["origin"][0].as_i64().unwrap(),meta["origin"][1].as_i64().unwrap()];
        let read=|name:&str| {let b=std::fs::read(dir.join(name)).unwrap();let a:Vec<i64>=b.chunks_exact(4).map(|q|i32::from_le_bytes(q.try_into().unwrap())as i64).collect();assert_eq!(a.len(),n*n);a};
        let h=read("height_mm.i32");let b=read("base_mm.i32");let t=read("top_mm.i32");let w=read("weak_bottom_mm.i32");
        let layer=std::fs::read(dir.join("cap_layer.u16")).unwrap();assert_eq!(layer.len(),n*n*2);
        let pin:Vec<bool>=std::fs::read(dir.join("pinned.u8")).unwrap().into_iter().map(|v|v!=0).collect();assert_eq!(pin.len(),n*n);
        let cr=meta["rates"][0].as_i64().unwrap();let wr=meta["rates"][1].as_i64().unwrap();
        let contacts:Vec<_>=layer.chunks_exact(2).enumerate().map(|(i,v)| {let l=u16::from_le_bytes(v.try_into().unwrap());if l==u16::MAX{None}else{Some(CapContact{source:(cell,l),base_mm:b[i],top_mm:t[i],weak_bottom_mm:w[i],cap_rate:cr,weak_rate:wr})}}).collect();
        drop(b);drop(t);drop(w);drop(layer);
        let mut p=ContrastParams::v1();p.subgrid_spacing_mm=Some(2000);assert_eq!(serde_json::to_value(&p).unwrap(),meta["params"],"captured scenario differs; never substitute default rates");
        let (coarse,st)=run(n,&h,&contacts,&pin,&p,origin,(0,n-1,0,n-1)).unwrap();
        assert_eq!(st.local_ledger_residual_mm,0);
        let mut volume=0i128;for a in &st.export_accounts {for s in [&a.cap_weathering_mm3,&a.weak_weathering_mm3,&a.collapsed_cap_mm3] {volume+=s.parse::<i128>().unwrap();}}
        assert_eq!(volume,-coarse.iter().map(|v|*v as i128).sum::<i128>()*AREA_MM2);
        let mut witnesses=Vec::new();let surface=st.subgrid.as_ref().unwrap();
        // Full fine deltas are retained for independent overlap and geometric
        // area audit, not hidden behind an aggregate pass.
        let mut bytes=Vec::with_capacity(surface.delta_xyz.len()*12);
        for v in &surface.delta_xyz {for a in v {bytes.extend_from_slice(&a.to_le_bytes());}}
        std::fs::write(dir.join("fine_delta_xyz.i32"),bytes).unwrap();
        let mut count=0;for (i,c) in contacts.iter().enumerate(){if c.is_some(){count+=1;if witnesses.len()<5 {witnesses.push([origin[0]+(i%n)as i64,origin[1]+(i/n)as i64,h[i]]);}}}
        let report=serde_json::json!({"input":meta,"fine":surface,"source_contact_posts":count,"export_accounts":st.export_accounts,"fine_removed_volume_mm3":volume.to_string(),"fine_removed_volume_m3":volume as f64/1e9,"restriction_residual_mm3":0,"fine_ledger_residual_mm":0,"input_witnesses":witnesses,"result_semantics":"complete actual native source-cell contact component, same sparse2m kernel; geometric whole-cell F6 integration and overlap checked separately"});
        std::fs::write(dir.join("result.json"),serde_json::to_vec_pretty(&report).unwrap()).unwrap();
        println!("whole source {cell}: {} changed fine posts, {}m3 removed, exact ledger",surface.changed_fine_posts,volume as f64/1e9);
    }
}

/// This mode leaves the calibrated E5 32 m solution untouched. Detail reconstructs
/// the same finite-volume degrees of freedom, rather than evolving another grid.
pub const RECONSTRUCTION_RULE: &str = "ESCARP-2 candidate: calibrated 32m E5 post-prism means reconstructed on 2m posts. Same-source cap/weak crossings locate the front; compact moving least-squares contact-relative branch planes and weighted crossing normals produce a continuous front. A one-fine-spacing face resolves its geometric discontinuity. A smooth cell-interior correction preserves each 16x16 fine-post sum exactly (integer mm); debit and credit relative to the owning coarse prism balance per post. No new erosion, rock identity, simulated time or collectible stock. Cell support is nearest-post [-8,7] fine indices, as ESCARP-1 restriction; its quadrature centre is one metre southwest of the native post. Fit radius 96m is a geometric reconstruction stencil, not a geological length. Outside fronts use a C1 cubic reference and the same conservative correction. This candidate may overshoot; its shape gates must be measured.";

#[derive(Clone)]
struct ContactEdge { xy:[f64;2], normal:[f64;2], source:(u32,u16) }
/// Small immutable native grid. Fine cells are evaluated on demand by pf-detail,
/// with the same evaluator used for geometry, collision and per-post accounting.
pub struct Reconstruction {
    pub n:usize,
    pub origin:[i64;2],
    h:Vec<i64>,
    contacts:Vec<Option<CapContact>>,
    fixed:Vec<bool>,
    edges:Vec<Vec<ContactEdge>>,
}
#[derive(Clone)]
pub struct ReconstructedCell {
    pub heights:[i32;256],
    pub debit_mm3:i128,
    pub credit_mm3:i128,
    pub max_correction_mm:i64,
    pub front_samples:usize,
}
impl Reconstruction {
    pub fn new(n:usize, origin:[i64;2], h:Vec<i64>, contacts:Vec<Option<CapContact>>, fixed:Vec<bool>)->R<Self> {
        if n<9 || h.len()!=n*n || contacts.len()!=n*n || fixed.len()!=n*n {return Err("ESCARP-2 bad native reconstruction lattice".into());}
        let mut edges=vec![Vec::new();n*n];
        for y in 0..n-1 {for x in 0..n-1 {let i=y*n+x;let Some(a)=contacts[i]else{continue};
            if fixed[i] {continue;}
            for (j,axis) in [(i+1,0),(i+n,1)] {let Some(b)=contacts[j]else{continue};
                if fixed[j]||a.source!=b.source||(h[i]>a.base_mm)==(h[j]>b.base_mm) {continue;}
                // Only actual retained cap/weak interfaces, not construction
                // excess above the column or a different source's rock.
                if h[i]>a.top_mm||h[j]>b.top_mm||h[i]<=a.weak_bottom_mm||h[j]<=b.weak_bottom_mm {continue;}
                let mut xy=[(origin[0]+x as i64)as f64*32.-1.,(origin[1]+y as i64)as f64*32.-1.];
                xy[axis]+=16.;let mut normal=[0.,0.];normal[axis]=if h[j]>b.base_mm {1.}else{-1.};
                edges[i].push(ContactEdge{xy,normal,source:a.source});
            }
        }}
        Ok(Self{n,origin,h,contacts,fixed,edges})
    }
    pub fn parent_height_mm(&self,gx:i64,gy:i64)->i64 {
        self.h[self.idx(gx-self.origin[0],gy-self.origin[1])]
    }
    fn idx(&self,x:i64,y:i64)->usize {assert!(x>=0&&y>=0&&x<self.n as i64&&y<self.n as i64,"ESCARP-2 outside native halo");y as usize*self.n+x as usize}
    fn cubic(&self,x:f64,y:f64)->f64 {
        let gx=x/32.-self.origin[0]as f64;let gy=y/32.-self.origin[1]as f64;
        let (ix,iy)=(gx.floor()as i64,gy.floor()as i64);let weights=|t:f64|[-0.5*t+t*t-0.5*t*t*t,1.-2.5*t*t+1.5*t*t*t,0.5*t+2.*t*t-1.5*t*t*t,-0.5*t*t+0.5*t*t*t];
        let a=weights(gx-ix as f64);let b=weights(gy-iy as f64);let mut z=0.;
        for v in 0..4 {for u in 0..4 {z+=a[u]*b[v]*self.h[self.idx(ix+u as i64-1,iy+v as i64-1)]as f64;}}z
    }
    fn raw(&self,x:f64,y:f64,source:Option<(u32,u16)>)->(f64,bool) {
        let reference=self.cubic(x,y);let Some(source)=source else{return(reference,false)};
        let gx=(x/32.).floor()as i64-self.origin[0];let gy=(y/32.).floor()as i64-self.origin[1];
        let mut centre=[0.,0.];let mut normal=[0.,0.];let mut ew=0.;let mut nearest=f64::INFINITY;
        let weight=|r2:f64| if r2<9216. {(1.-r2/9216.).powi(3)}else{0.};
        for yy in gy-3..=gy+3 {for xx in gx-3..=gx+3 {
            if xx<0||yy<0||xx>=self.n as i64||yy>=self.n as i64 {continue;}
            for e in &self.edges[self.idx(xx,yy)] {if e.source!=source {continue;}
                let r2=(x-e.xy[0]).powi(2)+(y-e.xy[1]).powi(2);let w=weight(r2);if w==0. {continue;}
                nearest=nearest.min(r2.sqrt());ew+=w;for k in 0..2 {centre[k]+=w*e.xy[k];normal[k]+=w*e.normal[k];}
            }
        }}
        let norm=normal[0].hypot(normal[1]);if ew==0.||norm<1e-12||nearest>=64. {return(reference,false);}
        let distance=((x-centre[0]/ew)*normal[0]+(y-centre[1]/ew)*normal[1])/norm;
        let mut matrices=[[[0.;4];3];2];let mut base_sum=0.;let mut base_w=0.;
        for yy in gy-3..=gy+3 {for xx in gx-3..=gx+3 {
            if xx<0||yy<0||xx>=self.n as i64||yy>=self.n as i64 {continue;}
            let i=self.idx(xx,yy);let Some(c)=self.contacts[i]else{continue};if c.source!=source||self.fixed[i]||self.h[i]>c.top_mm||self.h[i]<=c.weak_bottom_mm {continue;}
            let dx=(self.origin[0]+xx)as f64*32.-1.-x;let dy=(self.origin[1]+yy)as f64*32.-1.-y;
            let w=weight(dx*dx+dy*dy);if w==0. {continue;}
            let basis=[1.,dx/32.,dy/32.];let k=usize::from(self.h[i]>c.base_mm);let value=(self.h[i]-c.base_mm)as f64;
            for row in 0..3 {for col in 0..3 {matrices[k][row][col]+=w*basis[row]*basis[col];}matrices[k][row][3]+=w*basis[row]*value;}
            // A native contact plane, sampled at the shifted prism centres;
            // recover its local value by fitting it separately below.
            base_sum+=w*c.base_mm as f64;base_w+=w;
        }}
        let (Some(low),Some(high))=(solve_plane(matrices[0]),solve_plane(matrices[1]))else{return(reference,false)};
        // Contact-plane evaluation is bilinear only within the same source.
        let fx=(x+1.)/32.-self.origin[0]as f64;let fy=(y+1.)/32.-self.origin[1]as f64;
        let (ix,iy)=(fx.floor()as i64,fy.floor()as i64);let (u,v)=(fx-ix as f64,fy-iy as f64);
        let mut basal=0.;
        for (xx,yy,w) in [(ix,iy,(1.-u)*(1.-v)),(ix+1,iy,u*(1.-v)),(ix,iy+1,(1.-u)*v),(ix+1,iy+1,u*v)] {
            let Some(c)=self.contacts[self.idx(xx,yy)].filter(|c|c.source==source)else{return(reference,false)};basal+=w*c.base_mm as f64;
        }
        let _=(base_sum,base_w);
        if high<=low {return(reference,false);}
        let phase=(distance/2.+0.5).clamp(0.,1.);
        let face=basal+low+(high-low)*phase;
        let t=((64.-nearest)/32.).clamp(0.,1.);let blend=t*t*(3.-2.*t);
        (reference+(face-reference)*blend,true)
    }
    fn raw_band(&self,x:f64,y:f64,source:Option<(u32,u16)>,antialias:bool)->(f64,bool,f64) {
        let (h,f,w,_,_)=self.raw_band_toe(x,y,source,antialias,None);
        (h,f,w)
    }
    // The legacy arithmetic is retained verbatim before the optional toe delta.
    fn raw_band_toe(&self,x:f64,y:f64,source:Option<(u32,u16)>,antialias:bool,toe_tangent:Option<f64>)->(f64,bool,f64,f64,f64) {
        let reference=self.cubic(x,y);let Some(source)=source else{return(reference,false,0.,0.,0.)};
        let gx=(x/32.).floor()as i64-self.origin[0];let gy=(y/32.).floor()as i64-self.origin[1];
        let mut centre=[0.,0.];let mut normal=[0.,0.];let mut ew=0.;let mut nearest=f64::INFINITY;
        let weight=|r2:f64| if r2<25600. {(1.-r2/25600.).powi(3)}else{0.};
        for yy in gy-6..=gy+6 {for xx in gx-6..=gx+6 {
            if xx<0||yy<0||xx>=self.n as i64||yy>=self.n as i64 {continue;}
            for e in &self.edges[self.idx(xx,yy)] {if e.source!=source {continue;}
                let r2=(x-e.xy[0]).powi(2)+(y-e.xy[1]).powi(2);let w=weight(r2);if w==0. {continue;}
                nearest=nearest.min(r2.sqrt());ew+=w;for k in 0..2 {centre[k]+=w*e.xy[k];normal[k]+=w*e.normal[k];}
            }
        }}
        let norm=normal[0].hypot(normal[1]);if ew==0.||norm<1e-12||nearest>=128. {return(reference,false,0.,0.,0.);}
        let distance=((x-centre[0]/ew)*normal[0]+(y-centre[1]/ew)*normal[1])/norm;
        let mut matrices=[[[0.;4];3];2];let mut base_sum=0.;let mut base_w=0.;
        for yy in gy-6..=gy+6 {for xx in gx-6..=gx+6 {
            if xx<0||yy<0||xx>=self.n as i64||yy>=self.n as i64 {continue;}
            let i=self.idx(xx,yy);let Some(c)=self.contacts[i]else{continue};if c.source!=source||self.fixed[i]||self.h[i]>c.top_mm||self.h[i]<=c.weak_bottom_mm {continue;}
            let dx=(self.origin[0]+xx)as f64*32.-1.-x;let dy=(self.origin[1]+yy)as f64*32.-1.-y;
            let w=weight(dx*dx+dy*dy);if w==0. {continue;}
            let basis=[1.,dx/32.,dy/32.];let k=usize::from(self.h[i]>c.base_mm);let value=(self.h[i]-c.base_mm)as f64;
            for row in 0..3 {for col in 0..3 {matrices[k][row][col]+=w*basis[row]*basis[col];}matrices[k][row][3]+=w*basis[row]*value;}
            // A native contact plane, sampled at the shifted prism centres;
            // recover its local value by fitting it separately below.
            base_sum+=w*c.base_mm as f64;base_w+=w;
        }}
        let (Some(low),Some(high))=(solve_plane(matrices[0]),solve_plane(matrices[1]))else{return(reference,false,0.,0.,0.)};
        // Contact-plane evaluation is bilinear only within the same source.
        let fx=(x+1.)/32.-self.origin[0]as f64;let fy=(y+1.)/32.-self.origin[1]as f64;
        let (ix,iy)=(fx.floor()as i64,fy.floor()as i64);let (u,v)=(fx-ix as f64,fy-iy as f64);
        let mut basal=0.;
        for (xx,yy,w) in [(ix,iy,(1.-u)*(1.-v)),(ix+1,iy,u*(1.-v)),(ix,iy+1,(1.-u)*v),(ix+1,iy+1,u*v)] {
            let Some(c)=self.contacts[self.idx(xx,yy)].filter(|c|c.source==source)else{return(reference,false,0.,0.,0.)};basal+=w*c.base_mm as f64;
        }
        let _=(base_sum,base_w);
        if high<=low {return(reference,false,0.,0.,0.);}
        // ESCARP-4 resolves the sub-cell front over three 2m spacings.
        // Keep the upper edge at d=8m (1m below mode3's edge), spread
        // toward the foot to d=2m. A linear interior gives a phase-independent
        // maximum cell rise, unlike a one-cell step or a steep sigmoid.
        // This is a mesh-resolution closure, not a new erosion process.
        let phase=if antialias {((distance-2.)/6.).clamp(0.,1.)}
            else {((distance-8.)/2.+0.5).clamp(0.,1.)};
        let face=basal+low+(high-low)*phase;
        let t=((128.-nearest)/64.).clamp(0.,1.);let blend=t*t*(3.-2.*t);
        let original=reference+(face-reference)*blend;
        let Some(tangent)=toe_tangent else{return(original,true,blend,0.,0.)};
        let rise=(high-low)/1000.;
        // One parent spacing at the inherited class talus tangent separates
        // a small unresolved toe from a resolved cliff. Fade to zero over a
        // second parent spacing; large faces retain every original bit.
        // This is an explicit reconstruction closure, not a talus deposit.
        let q=(rise/(32.*tangent)-1.).clamp(0.,1.);
        let eligibility=1.-q*q*(3.-2.*q);
        if eligibility==0. {return(original,true,blend,0.,0.);}
        // C1 sub-cell profile, centred on the SAME fitted front as mode4.
        // The 1.5 factor is smoothstep's maximum derivative. Six fine
        // spacings give small toes a resolved shoulder and foot as well.
        let width=(1.5*rise/tangent).max(12.);
        let t=((distance-5.)/width+0.5).clamp(0.,1.);
        let phase_toe=t*t*(3.-2.*t);
        let delta=(high-low)*(phase_toe-phase)*blend*eligibility;
        // Scale balancing support with the actual step height. Without
        // this factor an arbitrarily small fitted step could receive a finite
        // source-wide correction and grow an artificial bump as rise -> 0.
        let weight=blend*eligibility*(rise/(32.*tangent))*(4.*t*(1.-t)).powi(2);
        (original,true,blend,delta,weight)
    }
    /// Finite-volume prolongation. Each integer fine-post prism has area4m².
    /// The smooth zero-edge bubble changes the shape, never the coarse mean.
    pub fn cell(&self,gx:i64,gy:i64)->ReconstructedCell {
        let (x,y)=(gx-self.origin[0],gy-self.origin[1]);let i=self.idx(x,y);
        assert!(x>=5&&y>=5&&x+5<self.n as i64&&y+5<self.n as i64,"ESCARP-2 needs five parent posts of halo");
        let source=self.contacts[i].map(|c|c.source);let target=self.h[i];let mut raw=[0.;256];let mut weights=[0.;256];let mut sum=0.;let mut ws=0.;let mut front_samples=0;
        for v in 0..16 {for u in 0..16 {let k=v*16+u;let px=gx as f64*32.+(u as f64-8.)*2.;let py=gy as f64*32.+(v as f64-8.)*2.;
            let (h,front)=if self.fixed[i] {(target as f64,false)}else{self.raw(px,py,source)};raw[k]=h;sum+=h;front_samples+=usize::from(front);
            let bubble=|a:usize| (std::f64::consts::PI*(a as f64+0.5)/16.).sin().powi(2);
            weights[k]=bubble(u)*bubble(v);ws+=weights[k];
        }}
        let correction=(target as f64*256.-sum)/ws;let mut heights=[0i32;256];let mut isum=0i64;let mut max_correction_mm=0;
        for k in 0..256 {let z=(raw[k]+correction*weights[k]).round()as i64;heights[k]=i32::try_from(z).expect("ESCARP-2 height overflow");isum+=z;max_correction_mm=max_correction_mm.max((z-raw[k].round()as i64).abs());}
        let residual=target*256-isum;
        heights[8*16+8]=i32::try_from(heights[8*16+8]as i64+residual).expect("ESCARP-2 rounding overflow");
        let mut debit_mm3=0;let mut credit_mm3=0;
        for &h in &heights {let v=(h as i64-target)as i128*AREA_MM2;if v>=0 {credit_mm3+=v;}else{debit_mm3-=v;}}
        assert_eq!(debit_mm3,credit_mm3,"ESCARP-2 per-post reconstruction volume failed");
        ReconstructedCell{heights,debit_mm3,credit_mm3,max_correction_mm,front_samples}
    }
}
fn solve_plane(mut a:[[f64;4];3])->Option<f64> {
    for k in 0..3 {let row=(k..3).max_by(|&i,&j|a[i][k].abs().total_cmp(&a[j][k].abs())).unwrap();a.swap(k,row);
        let d=a[k][k];if d.abs()<1e-10{return None;}for j in k..4 {a[k][j]/=d;}
        for i in 0..3 {if i==k {continue;}let f=a[i][k];for j in k..4 {a[i][j]-=f*a[k][j];}}
    }
    Some(a[0][3])
}

#[cfg(test)]
mod reconstruction_mode_tests {
    #[test]
    fn missing_native_parameter_carrier_is_refused() {
        let mut p=crate::t1b_ct::CtParamsV1::v1();
        assert!(serde_json::to_value(&p).unwrap().get("escarp_reconstruct").is_none());
        p.escarp_reconstruct=Some(2);
        assert!(p.validate().unwrap_err().contains("ct_base=1"));
        p.escarp_reconstruct=Some(4);
        assert!(p.validate().unwrap_err().contains("ct_base=1"));
        p.escarp_reconstruct=Some(6);
        assert!(p.validate().unwrap_err().contains("ct_base=1"));
        p.escarp_reconstruct=Some(99);
        assert!(p.validate().unwrap_err().contains("unsupported escarp"));
    }
}

pub const BAND_RECONSTRUCTION_RULE: &str = "ESCARP-3 candidate: zero-time reconstruction of the calibrated 32m E5 surface. Same-source cap/weak crossing normals and branch planes use a160m compact fit. The face is2m wide, centred8m toward the cap from the fitted crossing to represent the upper part of the coarse transition. Blend to the C1 cubic reference over64..128m. Redistribute the signed volume discrepancy with this smooth band weight, separately for each native source/cap identity in the retained window. No per-post mean constraint or cell-interior bubble. Exact integer-mm post-prism volume for the complete band; bounded1mm balanced rounding. Nonzero per-post deviations are reported. No new erosion, simulated time or collectible stock. Window-scoped band coefficients are not certified across distinct parent windows. Outside the band use uncorrected C1 cubic interpolation. Exact post-prism accounting is not triangle-integrated volume.";

pub const ANTIALIASED_BAND_RULE: &str = "ESCARP-4 opt-in: ESCARP-3 same-source contact fit, calibrated 32m E5 surface, and exact source/cap band post-prism accounting. The sub-cell signed-distance face is linear from d=2m to d=8m toward the cap (6m = three fine spacings), with constant branch values outside. Upper edge is1m downslope of mode3's edge; widen toward the foot, avoiding shoulder encroachment. This is a2m mesh-resolution closure, not a physical erosion length or new simulated time. Exact integer-mm band volume and bounded1mm rounding; per-post means may differ. Geometry/collision share this surface. Existing window-boundary and triangle-integrated-volume limitations remain. No new rock, stock, or geology.";

pub const TOE_RULE: &str = "ESCARP6 opt-in: ESCARP4 baseline plus a C1 sub-cell small-contact-step profile centred at the same signed distance d=5m. Width max(12m,1.5*branch_rise/native_sandstone_talus_tangent); eligibility is1 below32m*tangent, smoothly0 above64m*tangent. This uses the inherited conditional class angle as a reconstruction scale, NOT evidence of a deposit or additional erosion. Exact integer-mm incremental volume balance per source/cap on compact, step-height-weighted toe support only; resolved high cliffs receive no change. Native planes and final-height provenance lookup remain exact. Requires base_interp=c1. No new stock or simulated time. Inherits window-scoped and post-prism quadrature limits.";

#[derive(Default, Serialize)]
pub struct ToeAccount {
    pub source_cell:u32,
    pub cap_layer:u16,
    pub parent_posts:usize,
    pub changed_fine_prisms:usize,
    pub debit_mm3:String,
    pub credit_mm3:String,
    pub residual_mm3:String,
    pub max_delta_mm:i64,
    pub max_balance_correction_mm:i64,
}

/// A window-scoped geometric band, independent of tile request order. Its
/// source/cap accounts cannot exchange volume across geological identities.
/// All included parent prisms count, including the zero-weight fringe samples.
pub struct BandReconstruction {
    pub cells: HashMap<(i64,i64),ReconstructedCell>,
    pub accounts: Vec<BandAccount>,
    pub toe_accounts: Vec<ToeAccount>,
}
#[derive(Default, Serialize)]
pub struct BandAccount {
    pub source_cell:u32,
    pub cap_layer:u16,
    pub parent_posts:usize,
    pub fine_prisms:usize,
    pub debit_mm3:String,
    pub credit_mm3:String,
    pub residual_mm3:String,
    pub raw_volume_difference_mm3:String,
    pub max_post_mean_deviation_m:f64,
    pub rms_post_mean_deviation_m:f64,
    pub max_correction_m:f64,
}
impl BandReconstruction {
    pub fn build(g:&Reconstruction)->R<Self> {Self::build_profile(g,false)}
    pub fn build_antialiased(g:&Reconstruction)->R<Self> {Self::build_profile(g,true)}
    fn build_profile(g:&Reconstruction,antialias:bool)->R<Self> {
        use std::collections::BTreeSet;
        // Expand only real same-source contacts. The seven-post halo covers
        // the fit plus cubic stencil; no extrapolated or invented contacts.
        let mut candidates=BTreeSet::new();
        for (i,edges) in g.edges.iter().enumerate() {
            if edges.is_empty(){continue;}
            let (x,y)=((i%g.n)as i64,(i/g.n)as i64);
            for yy in (y-5).max(7)..=(y+6).min(g.n as i64-8) {
                for xx in (x-5).max(7)..=(x+6).min(g.n as i64-8) {
                    let j=g.idx(xx,yy);if g.fixed[j]{continue;}
                    if g.contacts[j].is_some_and(|c|edges.iter().any(|e|e.source==c.source)) {candidates.insert((yy,xx));}
                }
            }
        }
        struct Pending {xy:(i64,i64), source:(u32,u16), target:i64, cell:ReconstructedCell, weights:[u32;256]}
        let mut pending=Vec::new();
        // Sum in integer millimetres and integer weight units. The second pass
        // is exact, including rounding: each correction is a difference of two
        // cumulative quotients. Every sample differs from its ideal by <1mm.
        let mut totals=BTreeMap::<(u32,u16),(i128,i128)>::new();
        for (yy,xx) in candidates {
            let i=g.idx(xx,yy);let source=g.contacts[i].unwrap().source;
            let (gx,gy)=(g.origin[0]+xx,g.origin[1]+yy);
            let mut cell=ReconstructedCell{heights:[0;256],debit_mm3:0,credit_mm3:0,max_correction_mm:0,front_samples:0};
            let mut weights=[0u32;256];let mut sum=0i128;let mut ws=0i128;
            for v in 0..16 {for u in 0..16 {let k=v*16+u;
                let (h,front,w)=g.raw_band(gx as f64*32.+(u as f64-8.)*2.,gy as f64*32.+(v as f64-8.)*2.,Some(source),antialias);
                if !h.is_finite(){return Err("ESCARP-3 nonfinite native fit".into());}
                cell.heights[k]=i32::try_from(h.round()as i64).map_err(|_|"ESCARP-3 native fit overflow")?;
                weights[k]=(w*1_000_000_000.).round()as u32;
                sum+=cell.heights[k]as i128;ws+=weights[k]as i128;cell.front_samples+=usize::from(front);
            }}
            if ws==0 {continue;}
            let t=totals.entry(source).or_default();t.0+=g.h[i]as i128*256-sum;t.1+=ws;
            pending.push(Pending{xy:(gx,gy),source,target:g.h[i],cell,weights});
        }
        let mut cumulative=BTreeMap::<(u32,u16),(i128,i128)>::new();
        let mut accounts=BTreeMap::<(u32,u16),(BandAccount,i128,i128,f64)>::new();
        let mut cells=HashMap::new();
        for mut p in pending {
            let (residual,total_weight)=totals[&p.source];
            let prefix=cumulative.entry(p.source).or_default();let mut sum=0i64;
            for k in 0..256 {
                prefix.0+=p.weights[k]as i128;
                let q=residual*prefix.0/total_weight;let correction=q-prefix.1;prefix.1=q;
                p.cell.heights[k]=i32::try_from(p.cell.heights[k]as i128+correction).map_err(|_|"ESCARP-3 correction overflow")?;
                p.cell.max_correction_mm=p.cell.max_correction_mm.max(correction.abs()as i64);
                sum+=p.cell.heights[k]as i64;
                let dv=(p.cell.heights[k]as i128-p.target as i128)*AREA_MM2;
                if dv>=0 {p.cell.credit_mm3+=dv;}else{p.cell.debit_mm3-=dv;}
            }
            let (a,debit,credit,squared)=accounts.entry(p.source).or_default();
            a.source_cell=p.source.0;a.cap_layer=p.source.1;a.parent_posts+=1;a.fine_prisms+=256;
            let mean=(sum as f64/256.-p.target as f64)/1000.;
            a.max_post_mean_deviation_m=a.max_post_mean_deviation_m.max(mean.abs());*squared+=mean*mean;
            a.max_correction_m=a.max_correction_m.max(p.cell.max_correction_mm as f64/1000.);
            *debit+=p.cell.debit_mm3;*credit+=p.cell.credit_mm3;
            cells.insert(p.xy,p.cell);
        }
        let mut output=Vec::new();
        for (source,(mut a,debit,credit,squared)) in accounts {
            if debit!=credit||cumulative[&source].1!=totals[&source].0 {return Err("ESCARP-3 source-band volume residual".into());}
            a.debit_mm3=debit.to_string();a.credit_mm3=credit.to_string();a.residual_mm3=(credit-debit).to_string();
            a.raw_volume_difference_mm3=(-totals[&source].0*AREA_MM2).to_string();
            a.rms_post_mean_deviation_m=(squared/a.parent_posts as f64).sqrt();output.push(a);
        }
        println!("  escarp3 band: {} parent posts, {} native source accounts",cells.len(),output.len());
        Ok(Self{cells,accounts:output,toe_accounts:Vec::new()})
    }
    /// Apply an incremental, zero-volume toe reconstruction to an existing
    /// ESCARP4 band. All balancing weights vanish on the resolved cliff.
    /// The original source-band certificate remains valid; the incremental
    /// certificate below accounts for the changed distribution separately.
    pub fn soften_toes(&mut self,g:&Reconstruction,tangent:f64)->R<()> {
        if !tangent.is_finite()||tangent<=0. {return Err("ESCARP6 needs a positive native class talus tangent".into());}
        struct Pending {xy:(i64,i64),source:(u32,u16),delta:[i64;256],weights:[u32;256]}
        let mut keys:Vec<_>=self.cells.keys().copied().collect();keys.sort_by_key(|&(x,y)|(y,x));
        let mut pending=Vec::new();let mut totals=BTreeMap::<(u32,u16),(i128,i128)>::new();
        for (gx,gy) in keys {
            let i=g.idx(gx-g.origin[0],gy-g.origin[1]);let source=g.contacts[i].unwrap().source;
            let mut delta=[0i64;256];let mut weights=[0u32;256];let mut sum=0i128;let mut ws=0i128;
            for v in 0..16 {for u in 0..16 {let k=v*16+u;
                let (_,_,_,d,w)=g.raw_band_toe(gx as f64*32.+(u as f64-8.)*2.,gy as f64*32.+(v as f64-8.)*2.,Some(source),true,Some(tangent));
                if !d.is_finite()||!w.is_finite(){return Err("ESCARP6 nonfinite toe fit".into());}
                delta[k]=d.round()as i64;weights[k]=(w*1_000_000_000.).round()as u32;
                sum+=delta[k]as i128;ws+=weights[k]as i128;
            }}
            if ws==0 {if sum!=0{return Err("ESCARP6 toe delta without balancing support".into());}continue;}
            let total=totals.entry(source).or_default();total.0+=sum;total.1+=ws;
            pending.push(Pending{xy:(gx,gy),source,delta,weights});
        }
        let mut cumulative=BTreeMap::<(u32,u16),(i128,i128)>::new();
        let mut accounts=BTreeMap::<(u32,u16),(ToeAccount,i128,i128)>::new();
        for p in pending {
            let (sum,ws)=totals[&p.source];let prefix=cumulative.entry(p.source).or_default();
            let (account,debit,credit)=accounts.entry(p.source).or_default();
            account.source_cell=p.source.0;account.cap_layer=p.source.1;account.parent_posts+=1;
            let c=self.cells.get_mut(&p.xy).unwrap();let target=g.parent_height_mm(p.xy.0,p.xy.1);
            c.debit_mm3=0;c.credit_mm3=0;
            for k in 0..256 {
                prefix.0+=p.weights[k]as i128;let q=-sum*prefix.0/ws;let correction=q-prefix.1;prefix.1=q;
                let d=p.delta[k]as i128+correction;
                c.heights[k]=i32::try_from(c.heights[k]as i128+d).map_err(|_|"ESCARP6 toe height overflow")?;
                account.changed_fine_prisms+=usize::from(d!=0);
                account.max_delta_mm=account.max_delta_mm.max(d.abs()as i64);
                account.max_balance_correction_mm=account.max_balance_correction_mm.max(correction.abs()as i64);
                if d>=0{*credit+=d*AREA_MM2}else{*debit-=d*AREA_MM2}
                let dv=(c.heights[k]as i128-target as i128)*AREA_MM2;
                if dv>=0{c.credit_mm3+=dv}else{c.debit_mm3-=dv}
            }
        }
        for (source,(mut account,debit,credit)) in accounts {
            if debit!=credit||cumulative[&source].1!=-totals[&source].0{return Err("ESCARP6 toe volume residual".into());}
            account.debit_mm3=debit.to_string();account.credit_mm3=credit.to_string();account.residual_mm3=(credit-debit).to_string();
            self.toe_accounts.push(account);
        }
        Ok(())
    }
    pub fn cell(&self,g:&Reconstruction,gx:i64,gy:i64)->ReconstructedCell {
        if let Some(c)=self.cells.get(&(gx,gy)){return c.clone();}
        let i=g.idx(gx-g.origin[0],gy-g.origin[1]);let target=g.h[i];
        let mut c=ReconstructedCell{heights:[0;256],debit_mm3:0,credit_mm3:0,max_correction_mm:0,front_samples:0};
        for v in 0..16 {for u in 0..16 {
            let h=if g.fixed[i]{target}else{g.cubic(gx as f64*32.+(u as f64-8.)*2.,gy as f64*32.+(v as f64-8.)*2.).round()as i64};
            c.heights[v*16+u]=i32::try_from(h).expect("ESCARP-3 cubic height overflow");
            let dv=(h-target)as i128*AREA_MM2;if dv>=0{c.credit_mm3+=dv;}else{c.debit_mm3-=dv;}
        }}c
    }
}
