//! Owner's selective native SFD area closure. Same construction routing law as
//! t1b_lem::route, with global addresses, exact ground and explicit missing data.
//! A candidate is published only when every donor, flat and capture neighbour
//! of every requested free post is known. Missing sources are never estimated.
use crate::t1b_surface_completion::{Patch,Rect};
use crate::t1b_surface_domain::{PostId,D8,content_digest};
use crate::t1b_surface_stages::{NativeOwner,Stage};
use crate::t1b_rivernet::EvalLat;
use crate::t1b::{FREE,T1B_SPACING_MM,T1B_SPACING_M};
use serde::{Serialize,Deserialize};
use std::collections::{BTreeMap,BTreeSet};
use std::cell::RefCell;
use std::sync::Arc;
type R<T>=Result<T,String>;
const NONE:u32=u32::MAX;
const BLOCK:i64=128;
pub const RULE:&str="ls_native_area_sfd_v1: native construction single-flow steepest descent on exact mm ground's fill-to-spill potential, flats to canonical global-address flood parent; immutable 32m own area times native f32 runoff weight; fixed native boundaries receive local contributions but are not fresh sources; complete height/flat/donor/capture dependencies before publishing requested free-post areas; no crop-ring outlet, accumulated-area source or unknown inflow estimate; native coarse channel totals remain separately owned";
#[derive(Clone,Copy,Debug)]
pub struct Sample{pub id:PostId,pub height_mm:i32,pub fixed:u8,pub link:u32,pub own_m2:f64}
pub struct Grid{pub rect:Rect,pub z:Vec<f64>,pub fixed:Vec<u8>,pub links:Vec<u32>,pub own:Vec<f64>,pub present:Vec<bool>,pub source_pin:String}
impl Grid{
    pub fn new(rect:Rect,pin:String,limit:usize)->R<Self>{let n=rect.count().map_err(|e|e.to_string())?;if n>limit{return Err(format!("NeedsDomain resource: {n} working posts exceed {limit}"));}Ok(Self{rect,z:vec![f64::NAN;n],fixed:vec![FREE;n],links:vec![NONE;n],own:vec![0.0;n],present:vec![false;n],source_pin:pin})}
    pub fn ingest(&mut self,samples:&[Sample])->R<()>{for p in samples{let i=self.rect.index(p.id).ok_or("area sample outside working rectangle")?;let z=p.height_mm as f64/1000.0;
        if !p.own_m2.is_finite()||p.own_m2<0.0||![0,1,2,4,5].contains(&p.fixed){return Err("invalid native area sample".into());}
        if self.present[i]&&(self.z[i].to_bits()!=z.to_bits()||self.fixed[i]!=p.fixed||self.links[i]!=p.link||self.own[i].to_bits()!=p.own_m2.to_bits()){return Err(format!("stale canonical area input at {:?}",p.id));}
        self.z[i]=z;self.fixed[i]=p.fixed;self.links[i]=p.link;self.own[i]=p.own_m2;self.present[i]=true;
    }Ok(())}
    fn patch(&self)->Patch<'_>{Patch{rect:self.rect,ground:&self.z,fixed:&self.fixed,present:&self.present,predecessor:&self.source_pin,native_chart:Some(&self.source_pin)}}
    fn inside_chart(id:PostId)->bool{let h=vz_tilekey::CHART_HALF_EXTENT_MM/T1B_SPACING_MM;id.gx.abs()<=h&&id.gy.abs()<=h}
    fn missing_neighbours(&self,i:usize)->Vec<PostId>{(0..8).map(|k|self.rect.id(i).neighbour(k)).filter(|id|Self::inside_chart(*id)&&self.rect.index(*id).map(|j|!self.present[j]).unwrap_or(true)).collect()}
    fn neighbour(&self,i:usize,k:usize)->Option<usize>{self.rect.neighbour(i,k).filter(|&j|self.present[j])}
}
#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct Certificate{pub rule:String,pub source_pin:String,pub input_digest:String,pub route_digest:String,pub requested_free_posts:usize,pub unresolved_requested_posts:usize,pub known_posts:usize,pub known_source_m2:f64,pub known_boundary_export_m2:f64,pub unresolved_working_sources:usize,pub iterations:usize,pub fetched_blocks:usize,pub missing:Vec<PostId>,pub pass:bool}
pub struct Routed{pub area:Vec<f64>,pub receiver:Vec<u32>,pub terminal:Vec<u32>,pub bad_root:Vec<u32>,pub order:Vec<u32>,pub certificate:Certificate}
/// Exposed for real-data/falsifier tests. The complete atlas remains owned by
/// NativeAreaOwner; this provisional solve cannot claim whole-domain balance
/// when unclosed working sources remain outside the requested component.
pub fn solve(grid:&Grid,targets:&[usize])->R<Routed>{
    let p=grid.patch();let b=crate::t1b_ls::condition_detail_bounds(&p).map_err(|e|e.to_string())?;let n=grid.z.len();
    let mut pos=vec![NONE;n];for(k,&i)in b.order.iter().enumerate(){pos[i as usize]=k as u32;}
    let mut receiver=vec![NONE;n];
    for i in 0..n{if !grid.present[i]||grid.fixed[i]!=FREE||!b.upper[i].is_finite(){continue;}
        let mut best:Option<(f64,PostId,usize)>=None;
        for k in 0..8{if let Some(j)=grid.neighbour(i,k){if b.upper[j]<b.upper[i]&&pos[j]<pos[i]{let len=T1B_SPACING_M*if D8[k].0!=0&&D8[k].1!=0{2f64.sqrt()}else{1.0};let slope=(b.upper[i]-b.upper[j])/len;let id=grid.rect.id(j);if best.map(|(s,q,_)|slope>s||(slope==s&&id<q)).unwrap_or(true){best=Some((slope,id,j));}}}}
        receiver[i]=best.map(|(_,_,j)|j as u32).unwrap_or(b.parent[i]);
    }
    let mut bad=vec![NONE;n];let mut frontier=vec![false;n];let mut uncertain=vec![false;n];
    for i in 0..n{if grid.present[i]&&grid.fixed[i]==FREE{frontier[i]=!grid.missing_neighbours(i).is_empty();if frontier[i]{bad[i]=i as u32;}uncertain[i]=!b.height_exact(i);}}
    // A lower<upper component has a lower-bound path to an unknown frontier.
    // That path cannot pass an exact post: its known outlet would also make
    // the upstream minimax exact. Assign a real frontier witness per component.
    let mut seen=vec![false;n];let mut component=Vec::new();let mut stack=Vec::new();
    for start in 0..n{if !uncertain[start]||seen[start]{continue;}component.clear();stack.push(start);seen[start]=true;let mut root=NONE;
        while let Some(i)=stack.pop(){component.push(i);if frontier[i]{root=root.min(i as u32);}for k in 0..8{if let Some(j)=grid.neighbour(i,k){if uncertain[j]&&!seen[j]{seen[j]=true;stack.push(j);}}}}
        if root==NONE{return Err("uncertain area height has no actual missing-frontier witness".into());}for &i in &component{bad[i]=root;}
    }
    // Complete equal-potential components, including flat portions whose rank
    // can change when another seed arrives. Actual fixed boundaries stop taint.
    seen.fill(false);
    for start in 0..n{if !grid.present[start]||grid.fixed[start]!=FREE||seen[start]{continue;}component.clear();stack.push(start);seen[start]=true;let mut root=NONE;
        while let Some(i)=stack.pop(){component.push(i);root=root.min(bad[i]);for k in 0..8{if let Some(j)=grid.neighbour(i,k){if grid.fixed[j]==FREE&&!seen[j]&&b.upper[j]==b.upper[i]{seen[j]=true;stack.push(j);}}}}
        if root!=NONE{for &i in &component{bad[i]=root;}}
    }
    // A neighbouring uncertain post can capture a donor even if its current
    // provisional receiver points elsewhere. Preserve that dependency.
    let seeds=bad.clone();for i in 0..n{if seeds[i]!=NONE{for k in 0..8{if let Some(j)=grid.neighbour(i,k){if grid.fixed[j]==FREE{bad[j]=bad[j].min(seeds[i]);}}}}}
    let geometry_bad=bad.clone();
    let mut area=vec![0.0;n];let mut input=0.0;let mut unresolved_sources=0;
    for i in 0..n{if grid.present[i]&&grid.fixed[i]==FREE{area[i]=grid.own[i];input+=area[i];if receiver[i]==NONE{unresolved_sources+=1;}}}
    for &iu in b.order.iter().rev(){let i=iu as usize;let r=receiver[i];if r!=NONE{let j=r as usize;area[j]+=area[i];bad[j]=bad[j].min(bad[i]);}}
    let mut terminal=vec![NONE;n];for &iu in &b.order{let i=iu as usize;terminal[i]=if grid.fixed[i]!=FREE{i as u32}else{let r=receiver[i];if r==NONE{NONE}else{terminal[r as usize]}};}
    let mut path_bad=geometry_bad.clone();for &iu in &b.order{let i=iu as usize;let r=receiver[i];if r!=NONE{path_bad[i]=path_bad[i].min(path_bad[r as usize]);}}
    let mut missing=BTreeSet::new();let mut count=0;let mut unresolved=0;
    for &i in targets{if i>=n||!grid.present[i]{return Err("requested area post is absent".into());}if grid.fixed[i]!=FREE{continue;}count+=1;
        let root=bad[i].min(path_bad[i]);if root!=NONE{unresolved+=1;for id in grid.missing_neighbours(root as usize){missing.insert(id);}}
        else if terminal[i]==NONE||!area[i].is_finite(){return Err("certified area has no supported native terminal or finite flux".into());}
    }
    let export:f64=(0..n).filter(|&i|grid.present[i]&&grid.fixed[i]!=FREE).map(|i|area[i]).sum();
    if !input.is_finite()||!export.is_finite()||area.iter().any(|x|!x.is_finite()){return Err("nonfinite native area accounting".into());}
    if unresolved_sources==0&&(input-export).abs()>1e-9*input.max(1.0){return Err("native area known-source accounting failure".into());}
    let mut h=blake3::Hasher::new();h.update(RULE.as_bytes());h.update(b.input_digest.as_bytes());
    for i in 0..n{if grid.present[i]{h.update(&grid.own[i].to_le_bytes());h.update(&grid.links[i].to_le_bytes());h.update(&receiver[i].to_le_bytes());h.update(&area[i].to_le_bytes());}}
    let certificate=Certificate{rule:RULE.into(),source_pin:grid.source_pin.clone(),input_digest:b.input_digest,route_digest:h.finalize().to_hex().to_string(),requested_free_posts:count,unresolved_requested_posts:unresolved,known_posts:grid.present.iter().filter(|v|**v).count(),known_source_m2:input,known_boundary_export_m2:export,unresolved_working_sources:unresolved_sources,iterations:0,fetched_blocks:0,missing:missing.into_iter().collect(),pass:unresolved==0};
    Ok(Routed{area,receiver,terminal,bad_root:bad,order:b.order,certificate})
}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct OutputPost{pub id:PostId,pub area_m2:f64,pub receiver:Option<PostId>,pub terminal:PostId,pub terminal_fixed:u8,pub native_attachment:u32}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Output{pub certificate:Certificate,pub posts:Vec<OutputPost>,#[serde(default,skip_serializing_if="Option::is_none")]pub registry:Option<SfdRegistry>,pub digest:String}
pub struct NativeAreaOwner<'a,'b,'c>{publish_channels:bool,root:std::path::PathBuf,ground:&'a NativeOwner<'b,'c>,blocks:RefCell<BTreeMap<(i64,i64),Arc<Vec<Sample>>>>}
impl<'a,'b,'c> NativeAreaOwner<'a,'b,'c>{
    pub fn new(ground:&'a NativeOwner<'b,'c>,root:&std::path::Path,publish_channels:bool)->R<Self>{std::fs::create_dir_all(root).map_err(|e|e.to_string())?;Ok(Self{publish_channels,root:root.into(),ground,blocks:RefCell::new(BTreeMap::new())})}
    fn block(&self,tx:i64,ty:i64)->R<Arc<Vec<Sample>>>{
        if let Some(b)=self.blocks.borrow().get(&(ty,tx)){return Ok(b.clone());}
        let half=vz_tilekey::CHART_HALF_EXTENT_MM/T1B_SPACING_MM;let x0=tx*BLOCK;let y0=ty*BLOCK;
        if x0>half||y0>half||x0+BLOCK-1 < -half||y0+BLOCK-1 < -half{return Err("area fetch outside native chart".into());}
        let lat=EvalLat{gx0:x0.clamp(-half,half-BLOCK+1),gy0:y0.clamp(-half,half-BLOCK+1),n:BLOCK as usize,dx_mm:T1B_SPACING_MM};
        let(z,fixed,links,_)=self.ground.final_ground(&lat)?;let ctx=self.ground.native_context();let mut out=Vec::new();
        for i in 0..lat.nodes(){let(gx,gy)=lat.g(i);if gx.div_euclid(BLOCK)!=tx||gy.div_euclid(BLOCK)!=ty{continue;}
            let weight=if ctx.params.climate_runoff.is_some(){let(la,lo)=ctx.chart.inverse_mm(gx*T1B_SPACING_MM,gy*T1B_SPACING_MM).map_err(|e|e.to_string())?;ctx.t0.climate.as_ref().ok_or("native area climate missing")?.weight_at_deg(ctx.t0,la*vz_tilekey::chart::RAD2DEG,lo*vz_tilekey::chart::RAD2DEG)as f32}else{1.0};
            let height_mm=i32::try_from(vz_tilekey::chart::quantize_mm(z[i]*1000.0)).map_err(|_|"area ground overflow")?;
            out.push(Sample{id:PostId{gx,gy},height_mm,fixed:fixed[i],link:links[i],own_m2:T1B_SPACING_M*T1B_SPACING_M*weight as f64});
        }
        let out=Arc::new(out);self.blocks.borrow_mut().insert((ty,tx),out.clone());Ok(out)
    }
    /// Start from the already owned evaluation surface, then fetch only actual
    /// donor/capture/flat dependencies. Publication targets include the halo.
    pub fn query(&self,lat:&EvalLat,z:&[f64],fixed:&[u8],rw:&[f32],targets:&[PostId])->R<Output>{
        let(_,_,links)=self.ground.source_samples(Stage::VfDelta,lat)?;
        let mut initial=Vec::with_capacity(lat.nodes());for i in 0..lat.nodes(){let(gx,gy)=lat.g(i);initial.push(Sample{id:PostId{gx,gy},height_mm:i32::try_from(vz_tilekey::chart::quantize_mm(z[i]*1000.0)).map_err(|_|"initial area ground overflow")?,fixed:fixed[i],link:links[i],own_m2:T1B_SPACING_M*T1B_SPACING_M*if rw.is_empty(){1.0}else{rw[i]as f64}});}
        let mut keys=BTreeSet::<(i64,i64)>::new();let mut required:BTreeSet<PostId>=targets.iter().copied().collect();
        let threshold=crate::t1b_rivernet::window_params(self.ground.native_context().params,0).threshold_m3_s();
        for iteration in 0..64{
            let mut x0=lat.gx0;let mut y0=lat.gy0;let mut x1=lat.gx0+lat.n as i64-1;let mut y1=lat.gy0+lat.n as i64-1;
            for &(ty,tx)in &keys{x0=x0.min(tx*BLOCK);y0=y0.min(ty*BLOCK);x1=x1.max((tx+1)*BLOCK-1);y1=y1.max((ty+1)*BLOCK-1);}
            let rect=Rect{gx0:x0,gy0:y0,width:(x1-x0+1)as usize,height:(y1-y0+1)as usize};let mut grid=Grid::new(rect,self.ground.source_pin().into(),8_000_000)?;grid.ingest(&initial)?;
            for &(ty,tx)in &keys{grid.ingest(&self.block(tx,ty)?)?;}
            let indices:Vec<usize>=required.iter().map(|id|rect.index(*id).ok_or("area target out of domain")).collect::<Result<_,_>>()?;let mut r=solve(&grid,&indices)?;r.certificate.iterations=iteration;r.certificate.fetched_blocks=keys.len();
            eprintln!("ls area closure round {iteration}: {} requested free, {} unresolved, {} loaded, {} fetched blocks",r.certificate.requested_free_posts,r.certificate.unresolved_requested_posts,r.certificate.known_posts,keys.len());
            if r.certificate.pass{
                let registered=if self.publish_channels{channel_set(&grid,&r,&indices,self.ground.native_context().t0.runoff_m_per_yr,threshold)?}else{BTreeSet::new()};
                let old_required=required.len();for &i in &registered{required.insert(rect.id(i));}if required.len()!=old_required{continue;}
                let registry=if self.publish_channels{Some(make_registry(self.ground.native_context(),&grid,&r,&registered,threshold)?)}else{None};
                if let Some(reg)=&registry{audit_registry(self.ground.native_context(),&grid,&r,reg,&indices)?;}
                if let Some(reg)=&registry{exercise_native_falsifiers(self.ground.native_context(),&grid,&r,reg,&indices)?;}
                let mut posts=Vec::new();for &id in targets{let i=rect.index(id).ok_or("area publication target missing")?;if grid.fixed[i]!=FREE{continue;}let t=r.terminal[i]as usize;posts.push(OutputPost{id,area_m2:r.area[i],receiver:if r.receiver[i]==NONE{None}else{Some(rect.id(r.receiver[i]as usize))},terminal:rect.id(t),terminal_fixed:grid.fixed[t],native_attachment:grid.links[t]});}
                let mut out=Output{certificate:r.certificate,posts,registry,digest:String::new()};out.digest=content_digest(&out).map_err(|e|e.to_string())?;
                let file=self.root.join(format!("{}.json",out.digest));
                let bytes=serde_json::to_vec(&out).map_err(|e|e.to_string())?;
                if file.exists(){if std::fs::read(&file).map_err(|e|e.to_string())?!=bytes{return Err("stale native area artifact at content address".into());}}else{std::fs::write(file,bytes).map_err(|e|e.to_string())?;}
                return Ok(out);
            }
            let n=keys.len();for id in r.certificate.missing{keys.insert((id.gy.div_euclid(BLOCK),id.gx.div_euclid(BLOCK)));}if keys.len()==n{return Err("native area completion made no progress".into());}
        }
        Err("NeedsDomain: native area exceeded64 fetch rounds".into())
    }
}

#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Port{pub donor:PostId,pub receiver:PostId,pub area_m2:f64}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Channel{pub id:PostId,pub reach_id:String,pub receiver:PostId,pub area_m2:f64,pub own_source_m2:f64,pub native_attachment:u32,pub terminal:PostId,pub terminal_fixed:u8,pub source_t0_basin:u32,pub attachment_basin:u32,pub geometry_fixed:u8}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct SfdRegistry{pub rule:String,pub source_pin:String,pub route_digest:String,pub threshold_m3_s:f64,pub runoff_m_per_yr:f64,pub channels:Vec<Channel>,pub incoming_owner_ports:Vec<Port>,pub outgoing_native_ports:Vec<Port>,pub supra_threshold_posts:usize,pub wrong_basin_posts:usize,pub max_node_balance_relative:f64,pub local_source_m2:f64,pub incoming_owner_m2:f64,pub outgoing_native_m2:f64,pub relative_residual:f64,pub digest:String}
const REGISTRY_RULE:&str="ls_native_sfd_channels_v1: actual certified construction SFD graph; every reporting-threshold post plus connected downstream support; immutable own sources and explicit incoming/outgoing owner ports, never accumulated area relabelled as local rainfall; unique native attachment from actual terminal; geometry mask stays in canonical owner source; GUIDE is publication role only; coarse native channel totals are not reinjected or added to this graph";
fn channel_set(grid:&Grid,r:&Routed,targets:&[usize],runoff:f64,threshold:f64)->R<BTreeSet<usize>>{
    let mut registered=BTreeSet::new();
    for &start in targets{if grid.fixed[start]!=FREE||(r.area[start]as f32)as f64*runoff/crate::SECONDS_PER_YEAR<threshold{continue;}
        let mut i=start;let mut guard=0;
        while grid.fixed[i]==FREE{if !registered.insert(i){break;}let j=r.receiver[i];if j==NONE{return Err("registered source has no supported downstream path".into());}i=j as usize;guard+=1;if guard>grid.z.len(){return Err("registered source cycle".into());}}
    }
    Ok(registered)
}
fn make_registry(ctx:&crate::t1b_rivernet::RnContext<'_>,grid:&Grid,r:&Routed,registered:&BTreeSet<usize>,threshold:f64)->R<SfdRegistry>{
    let mut incoming=BTreeMap::<usize,Vec<usize>>::new();let mut pos=vec![NONE;grid.z.len()];for(k,&i)in r.order.iter().enumerate(){pos[i as usize]=k as u32;}
    for i in 0..grid.z.len(){let j=r.receiver[i];if j!=NONE&&registered.contains(&(j as usize)){incoming.entry(j as usize).or_default().push(i);}}
    let mut channels=Vec::new();let mut inputs=Vec::new();let mut outputs=Vec::new();let mut local=0.0;let mut incoming_total=0.0;let mut outgoing_total=0.0;let mut max_error=0f64;let mut supra=0;let mut wrong=0;
    for &i in registered{
        if r.bad_root[i]!=NONE{return Err("cannot register an unclosed area".into());}
        let t=r.terminal[i]as usize;let attachment=grid.links[t];
        if attachment==NONE||attachment as usize>=ctx.net.kind.len(){return Err(format!("native attachment adapter missing for terminal {:?} kind {}",grid.rect.id(t),grid.fixed[t]));}
        let native_kind=ctx.net.kind[attachment as usize];
        if ![crate::t1b_rivernet::K_T1A,crate::t1b_rivernet::K_MARINE,crate::t1b_rivernet::K_LAKE].contains(&native_kind)&&!ctx.net.conduit[attachment as usize]{return Err("registered terminal has no supported native graph attachment".into());}
        let id=grid.rect.id(i);let(la,lo)=ctx.chart.inverse_mm(id.gx*T1B_SPACING_MM,id.gy*T1B_SPACING_MM).map_err(|e|e.to_string())?;let cell=ctx.t0.cell_of_deg(la*vz_tilekey::chart::RAD2DEG,lo*vz_tilekey::chart::RAD2DEG);let basin=ctx.t0.final_basin(cell)?;
        let attachment_basin=ctx.skel.outlet_basin[attachment as usize];if attachment_basin!=basin{wrong+=1;}
        let mut expected=grid.own[i];let donors=incoming.entry(i).or_default();donors.sort_by_key(|&j|std::cmp::Reverse(pos[j]));
        for &j in donors.iter(){if r.bad_root[j]!=NONE{return Err("registered inflow port has unclosed source dependencies".into());}expected+=r.area[j];if !registered.contains(&j){incoming_total+=r.area[j];inputs.push(Port{donor:grid.rect.id(j),receiver:id,area_m2:r.area[j]});}}
        let error=(expected-r.area[i]).abs()/r.area[i].max(1.0);max_error=max_error.max(error);if error>1e-12{return Err("registered node source/inflow/outflow mismatch".into());}
        let j=r.receiver[i]as usize;if !registered.contains(&j){if grid.fixed[j]==FREE{return Err("registered downstream edge ends on unregistered free ground".into());}outgoing_total+=r.area[i];outputs.push(Port{donor:id,receiver:grid.rect.id(j),area_m2:r.area[i]});}
        local+=grid.own[i];if (r.area[i]as f32)as f64*ctx.t0.runoff_m_per_yr/crate::SECONDS_PER_YEAR>=threshold{supra+=1;}
        channels.push(Channel{id,reach_id:content_digest(&(REGISTRY_RULE,&ctx.rn_recipe_digest,id)).map_err(|e|e.to_string())?,receiver:grid.rect.id(j),area_m2:r.area[i],own_source_m2:grid.own[i],native_attachment:attachment,terminal:grid.rect.id(t),terminal_fixed:grid.fixed[t],source_t0_basin:basin,attachment_basin,geometry_fixed:FREE});
    }
    let residual=(local+incoming_total-outgoing_total).abs()/(local+incoming_total).max(1.0);
    if !residual.is_finite()||residual>1e-9{return Err("registered graph port/source ledger fails".into());}
    let mut out=SfdRegistry{rule:REGISTRY_RULE.into(),source_pin:grid.source_pin.clone(),route_digest:r.certificate.route_digest.clone(),threshold_m3_s:threshold,runoff_m_per_yr:ctx.t0.runoff_m_per_yr,channels,incoming_owner_ports:inputs,outgoing_native_ports:outputs,supra_threshold_posts:supra,wrong_basin_posts:wrong,max_node_balance_relative:max_error,local_source_m2:local,incoming_owner_m2:incoming_total,outgoing_native_m2:outgoing_total,relative_residual:residual,digest:String::new()};out.digest=content_digest(&out).map_err(|e|e.to_string())?;Ok(out)
}

/// Independent admission audit. Membership is derived from pre-registration
/// flow, ports from the actual receiver graph, and identities from native data.
fn audit_registry(ctx:&crate::t1b_rivernet::RnContext<'_>,grid:&Grid,r:&Routed,reg:&SfdRegistry,targets:&[usize])->R<()>{
    let threshold=crate::t1b_rivernet::window_params(ctx.params,0).threshold_m3_s();
    if reg.rule!=REGISTRY_RULE||reg.source_pin!=grid.source_pin||reg.route_digest!=r.certificate.route_digest||reg.threshold_m3_s!=threshold||reg.runoff_m_per_yr!=ctx.t0.runoff_m_per_yr{return Err("stale native registration binding".into());}
    let mut copy=reg.clone();copy.digest.clear();if content_digest(&copy).map_err(|e|e.to_string())?!=reg.digest{return Err("stale native registration content".into());}
    let mut expected=BTreeSet::new();let mut todo:Vec<usize>=targets.iter().copied().filter(|&i|grid.fixed[i]==FREE&&(r.area[i]as f32)as f64*ctx.t0.runoff_m_per_yr/crate::SECONDS_PER_YEAR>=threshold).collect();
    while let Some(i)=todo.pop(){if grid.fixed[i]!=FREE||!expected.insert(i){continue;}let next=r.receiver[i];if next==NONE{return Err("unconnected threshold flow".into());}todo.push(next as usize);}
    let mut actual=BTreeMap::new();for node in &reg.channels{let i=grid.rect.index(node.id).ok_or("registered post outside owner graph")?;if actual.insert(i,node).is_some(){return Err("duplicated registered post".into());}}
    if expected.len()!=actual.len()||expected.iter().any(|i|!actual.contains_key(i)){return Err("native registration coverage/orphan failure".into());}
    for (&i,node)in &actual{
        let j=r.receiver[i]as usize;let t=r.terminal[i]as usize;
        if grid.fixed[i]!=FREE||r.bad_root[i]!=NONE||node.geometry_fixed!=FREE||node.receiver!=grid.rect.id(j)||node.terminal!=grid.rect.id(t)||node.terminal_fixed!=grid.fixed[t]||node.native_attachment!=grid.links[t]||node.area_m2!=r.area[i]||node.own_source_m2!=grid.own[i]{return Err("native registration changes geometry, source or edge".into());}
        if node.reach_id!=content_digest(&(REGISTRY_RULE,&ctx.rn_recipe_digest,node.id)).map_err(|e|e.to_string())?{return Err("native registration reach identity".into());}
        let(la,lo)=ctx.chart.inverse_mm(node.id.gx*T1B_SPACING_MM,node.id.gy*T1B_SPACING_MM).map_err(|e|e.to_string())?;
        let cell=ctx.t0.cell_of_deg(la*vz_tilekey::chart::RAD2DEG,lo*vz_tilekey::chart::RAD2DEG);
        if node.source_t0_basin!=ctx.t0.final_basin(cell)?||node.attachment_basin!=*ctx.skel.outlet_basin.get(node.native_attachment as usize).ok_or("native attachment index")?{return Err("registered basin identities differ from actual native/reference data".into());}
    }
    let mut expected_in=BTreeMap::new();let mut expected_out=BTreeMap::new();let mut local=0.0;
    for i in 0..grid.z.len(){let j=r.receiver[i];if j==NONE{continue;}let j=j as usize;
        if expected.contains(&j)&&!expected.contains(&i){if r.bad_root[i]!=NONE{return Err("unclosed source port".into());}expected_in.insert((grid.rect.id(i),grid.rect.id(j)),r.area[i]);}
        if expected.contains(&i){local+=grid.own[i];if !expected.contains(&j){if grid.fixed[j]==FREE{return Err("channel gap before native boundary".into());}expected_out.insert((grid.rect.id(i),grid.rect.id(j)),r.area[i]);}}
    }
    let actual_in:BTreeMap<_,_>=reg.incoming_owner_ports.iter().map(|p|((p.donor,p.receiver),p.area_m2)).collect();let actual_out:BTreeMap<_,_>=reg.outgoing_native_ports.iter().map(|p|((p.donor,p.receiver),p.area_m2)).collect();
    if actual_in.len()!=reg.incoming_owner_ports.len()||actual_out.len()!=reg.outgoing_native_ports.len()||actual_in!=expected_in||actual_out!=expected_out{return Err("duplicated, omitted or invented native flow port".into());}
    let input:f64=expected_in.values().sum();let output:f64=expected_out.values().sum();if (local+input-output).abs()>1e-9*(local+input).max(1.0){return Err("independent native source/port balance failure".into());}
    let same=|a:f64,b:f64|a.is_finite()&&b.is_finite()&&(a-b).abs()<=1e-12*a.abs().max(b.abs()).max(1.0);
    if !same(local,reg.local_source_m2)||!same(input,reg.incoming_owner_m2)||!same(output,reg.outgoing_native_m2)||!reg.relative_residual.is_finite()||reg.relative_residual>1e-9||!reg.max_node_balance_relative.is_finite()||reg.max_node_balance_relative>1e-12{return Err("native registration summary ledger disagrees with graph".into());}
    let wrong=reg.channels.iter().filter(|c|c.source_t0_basin!=c.attachment_basin).count();
    let supra=reg.channels.iter().filter(|c|(c.area_m2 as f32)as f64*ctx.t0.runoff_m_per_yr/crate::SECONDS_PER_YEAR>=threshold).count();
    if reg.wrong_basin_posts!=wrong||reg.supra_threshold_posts!=supra{return Err("native registration summary hides actual roles/basin disagreements".into());}
    Ok(())
}

/// Read-only actual-native route witness, usable on legacy parent and owner
/// candidates with the same explicitly supplied global post list.
pub fn audit_requested_routes(ctx:&crate::t1b_rivernet::RnContext<'_>,lat:&EvalLat,z:&[f64],fixed:&[u8],links:&[u32],rx:i64,ry:i64)->R<()>{
    let Some(dir)=std::env::var_os("VZ_NATIVE_ROUTE_AUDIT")else{return Ok(())};
    let file=std::env::var_os("VZ_NATIVE_ROUTE_TARGETS").ok_or("VZ_NATIVE_ROUTE_AUDIT needs explicit VZ_NATIVE_ROUTE_TARGETS")?;
    let ids:Vec<PostId>=serde_json::from_slice(&std::fs::read(file).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
    let zq:Vec<f64>=z.iter().map(|v|vz_tilekey::chart::quantize_mm(*v*1000.0)as f64/1000.0).collect();let route=crate::t1b_lem::route(lat.n,&zq,fixed,lat.dx());let term=crate::t1b_lem::terminals(&route);let mut rows=Vec::new();
    for id in ids{let Some(i)=lat.idx(id.gx,id.gy)else{continue};let t=term[i]as usize;if t>=lat.nodes(){return Err("native audit route missing terminal".into());}let link=links[t];let(gx,gy)=lat.g(t);
        let(la,lo)=ctx.chart.inverse_mm(id.gx*T1B_SPACING_MM,id.gy*T1B_SPACING_MM).map_err(|e|e.to_string())?;let cell=ctx.t0.cell_of_deg(la*vz_tilekey::chart::RAD2DEG,lo*vz_tilekey::chart::RAD2DEG);
        rows.push(serde_json::json!({"post":id,"height_mm":vz_tilekey::chart::quantize_mm(z[i]*1000.0),"fixed":fixed[i],"source_t0_cell":cell,"source_t0_basin":ctx.t0.final_basin(cell)?,"terminal":PostId{gx,gy},"terminal_fixed":fixed[t],"native_attachment":link,"attachment_basin":ctx.skel.outlet_basin.get(link as usize).copied(),"attachment_kind":ctx.net.kind.get(link as usize).copied()}));
    }
    let dir=std::path::PathBuf::from(dir);std::fs::create_dir_all(&dir).map_err(|e|e.to_string())?;let out=dir.join(format!("w_{rx}_{ry}.json"));if out.exists(){return Err("native audit refuses evidence overwrite".into());}
    std::fs::write(out,serde_json::to_vec_pretty(&serde_json::json!({"source":"actual native evaluated parent/candidate, local-route diagnostic","layer_digest":ctx.layer_digest,"recipe_digest":ctx.rn_recipe_digest,"rows":rows})).unwrap()).map_err(|e|e.to_string())
}


fn exercise_native_falsifiers(ctx:&crate::t1b_rivernet::RnContext<'_>,grid:&Grid,r:&Routed,reg:&SfdRegistry,targets:&[usize])->R<()>{
    let Some(root)=std::env::var_os("VZ_NATIVE_REGISTRY_FALSIFIERS")else{return Ok(())};if reg.channels.is_empty(){return Ok(());}
    let mut cases=Vec::new();
    let mut check=|name:&str,mut candidate:SfdRegistry|->R<()>{candidate.digest.clear();candidate.digest=content_digest(&candidate).map_err(|e|e.to_string())?;let result=audit_registry(ctx,grid,r,&candidate,targets);if result.is_ok(){return Err(format!("native falsifier {name} was incorrectly admitted"));}cases.push(serde_json::json!({"case":name,"rejected":true,"reason":result.unwrap_err()}));Ok(())};
    let mut bad=reg.clone();bad.channels.pop();check("deleted-channel",bad)?;
    let mut bad=reg.clone();bad.channels[0].receiver=bad.channels[0].id;check("false-receiver",bad)?;
    let mut bad=reg.clone();bad.channels[0].native_attachment=NONE;check("false-native-attachment",bad)?;
    let mut bad=reg.clone();bad.source_pin.push('x');check("stale-surface",bad)?;
    if !reg.incoming_owner_ports.is_empty(){let mut bad=reg.clone();bad.incoming_owner_ports[0].area_m2*=2.0;check("doubled-owner-inflow",bad)?;}
    let mut bad=reg.clone();bad.local_source_m2+=1024.0;check("invented-local-source",bad)?;
    if reg.channels.iter().any(|c|c.source_t0_basin!=c.attachment_basin){let mut bad=reg.clone();for c in &mut bad.channels{c.attachment_basin=c.source_t0_basin;}bad.wrong_basin_posts=0;check("hidden-reference-basin-disagreement",bad)?;}
    let root=std::path::PathBuf::from(root);std::fs::create_dir_all(&root).map_err(|e|e.to_string())?;let path=root.join(format!("{}.json",reg.digest));
    let bytes=serde_json::to_vec_pretty(&serde_json::json!({"registry":reg.digest,"native_channels":reg.channels.len(),"source":"actual native candidate, mutations only in independent audit copies","cases":cases})).unwrap();
    if path.exists(){if std::fs::read(&path).map_err(|e|e.to_string())?!=bytes{return Err("native falsifier evidence changed at the same registry digest".into());}}else{std::fs::write(path,bytes).map_err(|e|e.to_string())?;}Ok(())
}

/// Mutable physical MFD support. Optimistic missing-frontier seeds occur only
/// in the lower-bound proof, never in this route or its export ledger.
pub const MUTABLE_RULE:&str="native_mutable_mfd_support_v1: exact minimax bounds; complete equal-potential flat and neighbour-capture dependencies; all downhill slope^1.1/distance branches; native global-address flat parent; donor taint propagated through every positive branch; all sediment recipients to real fixed terminals; immediate bank stencils for supplied mobile posts; certificates bind the actual changed geometry and own sources; downstream deposition and bank activation require requalification before accepting a revision";
#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct MutableCertificate {
    pub rule:String,pub source_pin:String,pub geometry_digest:String,pub route_digest:String,
    pub requested_mobile_posts:usize,pub sediment_recipient_posts:usize,pub bank_neighbour_posts:usize,
    pub qualified_required_posts:usize,pub unresolved_required_posts:usize,
    pub known_source_m2:f64,pub known_boundary_export_m2:f64,pub unresolved_working_sources:usize,
    pub missing:Vec<PostId>,pub pass:bool,
}
pub struct MutableRoute {
    pub potential_m:Vec<f64>,pub parent:Vec<u32>,pub order:Vec<u32>,
    /// D8 weights. On a flat the only positive entry is the known-outlet parent.
    pub weights:Vec<[f64;8]>,pub effective_area_m2:Vec<f64>,pub bad_root:Vec<u32>,
    pub required:Vec<bool>,pub certificate:MutableCertificate,
}
/// Does not accept a physical update. Call again on the proposed geometry with
/// every post that now carries movable stock, including deposited/bank stock.
/// The process must reject the proposal if that certificate is not complete.
pub fn qualify_mutable(grid:&Grid,mobile:&[usize])->R<MutableRoute>{
    let patch=grid.patch();let bounds=crate::t1b_ls::condition_detail_bounds(&patch).map_err(|e|e.to_string())?;
    let n=grid.z.len();let mut mobile_seen=BTreeSet::new();
    for &i in mobile {if i>=n||!grid.present[i]||grid.fixed[i]!=FREE||!mobile_seen.insert(i){return Err("mutable support has absent, fixed or duplicate mobile target".into());}}
    let mut weights=vec![[0.0;8];n];let mut bad=vec![NONE;n];let mut frontier=vec![false;n];let mut uncertain=vec![false;n];
    for i in 0..n {if !grid.present[i]||grid.fixed[i]!=FREE {continue;}
        frontier[i]=!grid.missing_neighbours(i).is_empty();if frontier[i]{bad[i]=i as u32;}
        uncertain[i]=!bounds.height_exact(i);
        if !bounds.upper[i].is_finite(){continue;}
        for k in 0..8 {if let Some(j)=grid.neighbour(i,k) {if bounds.upper[i]>bounds.upper[j] {
            let d=32.0*if D8[k].0!=0&&D8[k].1!=0{2f64.sqrt()}else{1.0};
            weights[i][k]=((bounds.upper[i]-bounds.upper[j])/d).powf(1.1)/d;
        }}}
        let mut sum: f64=weights[i].iter().sum();
        if sum==0.0 {let parent=bounds.parent[i];if parent==NONE{return Err("known mutable potential has no native escape parent".into());}
            let k=(0..8).find(|&k|grid.neighbour(i,k)==Some(parent as usize)).ok_or("mutable flat parent is not a known neighbour")?;weights[i][k]=1.0;sum=1.0;
        }
        if !sum.is_finite()||sum<=0.0{return Err("nonfinite mutable MFD normalization".into());}
        for w in &mut weights[i]{*w/=sum;}
    }
    // Unknown escape paths can affect an entire minimax component. Associate
    // each with an actual fetchable frontier, not a fabricated basin wall.
    let mut seen=vec![false;n];let mut component=Vec::new();let mut stack=Vec::new();
    for start in 0..n {if !uncertain[start]||seen[start]{continue;}component.clear();stack.push(start);seen[start]=true;let mut root=NONE;
        while let Some(i)=stack.pop(){component.push(i);if frontier[i]{root=root.min(i as u32);}
            for k in 0..8{if let Some(j)=grid.neighbour(i,k){if uncertain[j]&&!seen[j]{seen[j]=true;stack.push(j);}}}}
        if root==NONE{return Err("mutable uncertain component lacks a fetchable frontier".into());}
        for &i in &component{bad[i]=root;}
    }
    // A newly discovered outlet can change the global parent on an entire flat,
    // even when the minimax height itself is already exact.
    seen.fill(false);
    for start in 0..n {if !grid.present[start]||grid.fixed[start]!=FREE||seen[start]{continue;}component.clear();stack.push(start);seen[start]=true;let mut root=NONE;
        while let Some(i)=stack.pop(){component.push(i);root=root.min(bad[i]);for k in 0..8{if let Some(j)=grid.neighbour(i,k){if grid.fixed[j]==FREE&&!seen[j]&&bounds.upper[i]==bounds.upper[j]{seen[j]=true;stack.push(j);}}}}
        if root!=NONE{for &i in &component{bad[i]=root;}}
    }
    // Potential changes next door can capture a previously differently routed
    // donor. This is a geometry dependency, independent of its present branch.
    let seeds=bad.clone();for i in 0..n{if seeds[i]!=NONE{for k in 0..8{if let Some(j)=grid.neighbour(i,k){if grid.fixed[j]==FREE{bad[j]=bad[j].min(seeds[i]);}}}}}
    let mut area=vec![0.0;n];let mut input=0.0;let mut unresolved_sources=0;
    for i in 0..n {if grid.present[i]&&grid.fixed[i]==FREE{if !grid.own[i].is_finite()||grid.own[i]<0.0{return Err("invalid mutable own source".into());}area[i]=grid.own[i];input+=area[i];if !bounds.upper[i].is_finite(){unresolved_sources+=1;}}}
    for &iu in bounds.order.iter().rev(){let i=iu as usize;
        for k in 0..8 {let w=weights[i][k];if w<=0.0{continue;}let j=grid.neighbour(i,k).ok_or("mutable branch receiver absent")?;
            area[j]+=area[i]*w;bad[j]=bad[j].min(bad[i]);
        }
    }
    // Required targets include every possible sediment recipient. No branch is
    // pruned by a visual guide, accumulated-area threshold or T0 basin label.
    let mut required=vec![false;n];stack.extend(mobile.iter().copied());let mut descendants=0;
    while let Some(i)=stack.pop(){if required[i]{continue;}required[i]=true;descendants+=1;
        for k in 0..8{if weights[i][k]>0.0{stack.push(grid.neighbour(i,k).ok_or("sediment recipient absent")?);}}
    }
    let mut missing=BTreeSet::new();let mut bank=vec![false;n];
    for &i in mobile {for k in 0..8 {let id=grid.rect.id(i).neighbour(k);if !Grid::inside_chart(id){continue;}
        if let Some(j)=grid.neighbour(i,k){bank[j]=true;required[j]=true;}else{missing.insert(id);}
    }}
    let mut qualified=0;let mut unresolved=0;
    for i in 0..n {if !required[i]{continue;}
        if bad[i]!=NONE {unresolved+=1;for id in grid.missing_neighbours(bad[i] as usize){missing.insert(id);}}
        else if !bounds.upper[i].is_finite(){return Err("mutable required post lacks a supported terminal".into());}
        else {qualified+=1;}
    }
    let export:f64=(0..n).filter(|&i|grid.present[i]&&grid.fixed[i]!=FREE).map(|i|area[i]).sum();
    if !input.is_finite()||!export.is_finite()||area.iter().any(|a|!a.is_finite()){return Err("nonfinite mutable water ledger".into());}
    if unresolved_sources==0&&(input-export).abs()>1e-9*input.max(1.0){return Err("mutable known source/export mismatch".into());}
    let mut hash=blake3::Hasher::new();hash.update(MUTABLE_RULE.as_bytes());hash.update(bounds.input_digest.as_bytes());
    for i in 0..n{if grid.present[i]{hash.update(&grid.own[i].to_le_bytes());for w in &weights[i]{hash.update(&w.to_le_bytes());}hash.update(&area[i].to_le_bytes());}}
    let certificate=MutableCertificate{rule:MUTABLE_RULE.into(),source_pin:grid.source_pin.clone(),geometry_digest:bounds.input_digest,route_digest:hash.finalize().to_hex().to_string(),requested_mobile_posts:mobile.len(),sediment_recipient_posts:descendants,bank_neighbour_posts:bank.iter().filter(|&&b|b).count(),qualified_required_posts:qualified,unresolved_required_posts:unresolved,known_source_m2:input,known_boundary_export_m2:export,unresolved_working_sources:unresolved_sources,pass:unresolved==0&&missing.is_empty(),missing:missing.into_iter().collect()};
    Ok(MutableRoute{potential_m:bounds.upper,parent:bounds.parent,order:bounds.order,weights,effective_area_m2:area,bad_root:bad,required,certificate})
}
