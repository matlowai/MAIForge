//! Native staged construction interface. The legacy path passes no owner.
//! Canonical stages are initial-condition construction, not process sediment.
use crate::t1b_rivernet::{EvalLat,Surface};
use serde::{Serialize,Deserialize};
pub type R<T>=std::result::Result<T,String>;
#[derive(Clone,Copy,Debug,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
pub enum Stage{Main,CoastFinal,VfDelta}
impl Stage{pub fn name(self)->&'static str{match self{Self::Main=>"main",Self::CoastFinal=>"coast_final",Self::VfDelta=>"vf_delta"}}}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Proof{
    pub stage:Stage,pub source_pin:String,pub input_digest:String,pub output_digest:String,
    pub requested_posts:usize,pub fetched_tiles:usize,pub loaded_posts:usize,pub completion_iterations:usize,
    pub canonical_raised_posts:u32,pub initial_raise_m3:f64,pub max_raise_m:f64,
    pub local_input_adoption_net_m3:f64,pub local_input_adoption_max_m:f64,
    pub read_set:Vec<String>,pub rule:String,
}
pub trait Conditioner{
    fn condition(&self,stage:Stage,lat:&EvalLat,z:&mut[f64],fixed:&mut[u8])->R<Proof>;
}
#[derive(Clone,Debug)]
pub struct Capture{pub z:Vec<f64>,pub fixed:Vec<u8>,pub links:Vec<u32>,pub dependencies:Vec<Proof>}
pub enum Evaluated{Surface(Surface),Captured(Capture)}
impl Evaluated{pub fn surface(self)->R<Surface>{match self{Self::Surface(s)=>Ok(s),Self::Captured(_)=>Err("unexpected stage capture in window evaluation".into())}}}
pub fn capture(z:&[f64],fixed:&[u8],sm:Option<&crate::t1b_seam::SeamTaint>)->Capture{
    Capture{z:z.to_vec(),fixed:fixed.to_vec(),links:Vec::new(),dependencies:sm.map(|s|s.owner.clone()).unwrap_or_default()}
}
pub fn apply(owner:&dyn Conditioner,stage:Stage,lat:&EvalLat,z:&mut[f64],fixed:&mut[u8],sm:Option<&mut crate::t1b_seam::SeamTaint>)->R<(u32,f64)>{
    let proof=owner.condition(stage,lat,z,fixed)?;let result=(proof.canonical_raised_posts,proof.max_raise_m);
    let sm=sm.ok_or("owner construction requires seam_cert=1 for predecessor proofs")?;
    // The owner replaced the entire requested stage from its canonical, closed
    // predecessor. This is not removal of taint from a locally filled surface.
    sm.taint.fill(false);sm.ext.fill(0);sm.owner.push(proof);Ok(result)
}

use std::{cell::RefCell,collections::{BTreeMap,BTreeSet},path::{Path,PathBuf},sync::Arc};
use crate::t1b_surface_completion::{Rect,WorkingSet};
use crate::t1b_surface_domain::{content_digest,PostId};
const TILE:i64=512;
const PAD:i64=320;
pub const NATIVE_RULE:&str="ls_native_stages_v1: unique 512-post half-open owner tiles at global integer addresses; initial-condition pre-main sources use the existing native 1153-square evaluation with 256-post margin plus 64-post halo; at the physical chart edge that evaluation is shifted wholly inside the chart, with explicit one-sided support and no invented external fields; subsequent stage sources use canonical closed predecessor fills; original unrestricted D8 construction minimax, exact missing-frontier lower bounds, no ring outlet; stage output replaces local construction input explicitly; this is initial construction, not process mass or water; actual chart edges inherit the native LS finite D8 graph: no exterior land edge, no ring terminal; native T1a/water exports remain explicit; global eikonal/other non-fill dependency qualification remains a separate gate";
#[derive(Clone,Copy,Debug,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
struct TileKey{stage:Stage,ty:i64,tx:i64}
#[derive(Clone,Serialize,Deserialize)]
struct TileMeta{key:TileKey,source_pin:String,digest:String,dependencies:Vec<Proof>,rule:String}
struct NativeTile{meta:TileMeta,z:Vec<f64>,fixed:Vec<u8>,links:Vec<u32>}
/// All caches belong to this native context. Cache hit/miss does not change a
/// working domain or its read set. There is no process-global current window.
pub struct NativeOwner<'a,'world>{ctx:&'a crate::t1b_rivernet::RnContext<'world>,root:PathBuf,pin:String,tiles:RefCell<BTreeMap<TileKey,Arc<NativeTile>>>,building:RefCell<BTreeSet<TileKey>>,max_working_posts:usize,legacy_root:Option<(PathBuf,String,String)>}
impl<'a,'world> NativeOwner<'a,'world>{
    pub fn new(ctx:&'a crate::t1b_rivernet::RnContext<'world>,root:&Path)->R<Self>{
        if ctx.params.seam_cert.is_none(){return Err("ls native owner requires seam_cert=1".into());}
        if ctx.params.hs.as_ref().is_some_and(|h|h.runs()){return Err("ls native owner: transient HS is not yet a qualified canonical downstream modifier".into());}
        if crate::t1b_ve::ve(ctx.params).and_then(|v|v.flank2.as_ref()).is_some_and(|v|v.phase!=0){return Err("ls native owner rejects legacy flank2 evolution; use owner-driven evolution".into());}
        let rn=include_str!("t1b_rivernet.rs");let rn=rn.split("pub fn evaluate_owned(").nth(1).ok_or("owner source anchor")?.split("// ------------------------------------------------------------------ lane CO inputs").next().unwrap();
        let pin=content_digest(&(NATIVE_RULE,&ctx.rn_recipe_digest,&ctx.layer_digest,ctx.params,rn,include_str!("t1b_coast.rs"),include_str!("t1b_surface_completion.rs"),include_str!("t1b_surface_stages.rs"))).map_err(|e|e.to_string())?;
        // A complete unknown-exterior proof is stronger than the inherited
        // closed-chart-land-edge contract. Reuse only this exact reviewed old
        // producer, under the same native context and unchanged non-fill kernels.
        let profile_path=root.join("pre-native-chart-contract.json");
        let legacy_root=if profile_path.exists(){
            let bytes=std::fs::read(&profile_path).map_err(|e|e.to_string())?;
            if blake3::hash(&bytes).to_hex().as_str()!="45ce4ca79ed6b530d482cd45efb2bf50136274a8d0ce455dcd97e8bef5b7ab6e"{return Err("stale chart-refinement source receipt".into());}
            let old:serde_json::Value=serde_json::from_slice(&bytes).map_err(|e|e.to_string())?;
            let text=|key:&str|old[key].as_str().ok_or_else(||format!("chart refinement missing {key}"));
            if text("evaluation_source")?!=rn||text("coast_source")?!=include_str!("t1b_coast.rs"){return Err("chart refinement cannot reuse changed non-fill kernels".into());}
            let old_pin=content_digest(&(text("rule")?,&ctx.rn_recipe_digest,&ctx.layer_digest,ctx.params,text("evaluation_source")?,text("coast_source")?,text("completion_source")?,text("stages_source")?)).map_err(|e|e.to_string())?;
            let directory=root.join(&old_pin);if directory.is_dir(){eprintln!("ls chart boundary refinement reuses strictly certified source {old_pin}");Some((directory,old_pin,text("rule")?.to_owned()))}else{None}
        }else{None};
        let root=root.join(&pin);std::fs::create_dir_all(&root).map_err(|e|e.to_string())?;
        Ok(Self{ctx,root,pin,tiles:RefCell::new(BTreeMap::new()),building:RefCell::new(BTreeSet::new()),max_working_posts:8_000_000,legacy_root})
    }
    pub fn source_pin(&self)->&str{&self.pin}
    pub fn native_context(&self)->&crate::t1b_rivernet::RnContext<'world>{self.ctx}
    /// Exact canonical pre-stage samples and original native attachments. Locators
    /// and cache population are not part of the sampled values.
    pub fn source_samples(&self,stage:Stage,lat:&EvalLat)->R<(Vec<f64>,Vec<u8>,Vec<u32>)>{
        let mut tiles=BTreeMap::new();let mut z=Vec::with_capacity(lat.nodes());let mut fixed=Vec::with_capacity(lat.nodes());let mut links=Vec::with_capacity(lat.nodes());
        for i in 0..lat.nodes(){let(gx,gy)=lat.g(i);let key=TileKey{stage,ty:gy.div_euclid(TILE),tx:gx.div_euclid(TILE)};
            if !tiles.contains_key(&key){tiles.insert(key,self.tile(key)?);}
            let t=&tiles[&key];let j=(gy.rem_euclid(TILE)*TILE+gx.rem_euclid(TILE))as usize;
            if !t.z[j].is_finite(){return Err("NeedsBoundary: sample outside supported native chart".into());}
            z.push(t.z[j]);fixed.push(t.fixed[j]);links.push(t.links[j]);
        }
        Ok((z,fixed,links))
    }
    pub fn final_ground(&self,lat:&EvalLat)->R<(Vec<f64>,Vec<u8>,Vec<u32>,Proof)>{
        let(mut z,mut fixed,links)=self.source_samples(Stage::VfDelta,lat)?;
        let proof=self.condition(Stage::VfDelta,lat,&mut z,&mut fixed)?;Ok((z,fixed,links,proof))
    }

    fn tile(&self,key:TileKey)->R<Arc<NativeTile>>{
        if let Some(t)=self.tiles.borrow().get(&key){return Ok(t.clone());}
        if !self.building.borrow_mut().insert(key){return Err(format!("canonical stage dependency cycle: {key:?}"));}
        let result=self.load_or_build(key);self.building.borrow_mut().remove(&key);
        let t=Arc::new(result?);let mut cache=self.tiles.borrow_mut();if cache.len()>=96{cache.pop_first();}cache.insert(key,t.clone());Ok(t)
    }
    fn load_or_build(&self,key:TileKey)->R<NativeTile>{
        let stem=format!("{}_{}_{}",key.stage.name(),key.tx,key.ty);let mut meta_path=self.root.join(format!("{stem}.json"));let mut data_path=self.root.join(format!("{stem}.bin"));
        let mut expected_pin=self.pin.as_str();let mut expected_rule=NATIVE_RULE;
        if !meta_path.exists(){if let Some((directory,pin,rule))=&self.legacy_root{let old=directory.join(format!("{stem}.json"));if old.exists(){meta_path=old;data_path=directory.join(format!("{stem}.bin"));expected_pin=pin;expected_rule=rule;}}}
        if meta_path.exists(){
            let meta:TileMeta=serde_json::from_slice(&std::fs::read(&meta_path).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
            if meta.key!=key||meta.source_pin!=expected_pin||meta.rule!=expected_rule{return Err("stale canonical tile metadata".into());}
            let b=std::fs::read(&data_path).map_err(|e|e.to_string())?;let n=(TILE*TILE)as usize;
            if b.len()!=n*13||blake3::hash(&b).to_hex().as_str()!=meta.digest{return Err("stale canonical tile bytes".into());}
            let z=b[..n*8].chunks_exact(8).map(|c|f64::from_le_bytes(c.try_into().unwrap())).collect();let fixed=b[n*8..n*9].to_vec();let links=b[n*9..].chunks_exact(4).map(|c|u32::from_le_bytes(c.try_into().unwrap())).collect();return Ok(NativeTile{meta,z,fixed,links});
        }
        if data_path.exists(){return Err(format!("partial tile retained without metadata: {}",data_path.display()));}
        let half=vz_tilekey::CHART_HALF_EXTENT_MM/crate::t1b::T1B_SPACING_MM;
        if key.tx*TILE>half||key.ty*TILE>half||(key.tx+1)*TILE-1 < -half||(key.ty+1)*TILE-1 < -half{return Err(format!("NeedsBoundary: outside native chart at {key:?}; no invented outlet"));}
        let n=(TILE+2*PAD+1)as usize;
        let lat=EvalLat{gx0:(key.tx*TILE-PAD).clamp(-half,half-n as i64+1),gy0:(key.ty*TILE-PAD).clamp(-half,half-n as i64+1),n,dx_mm:crate::t1b::T1B_SPACING_MM};
        eprintln!("ls canonical source {:?} {} {}",key.stage,key.tx,key.ty);
        let c=self.ctx;let mut lp=crate::t1b_rivernet::lat_prep(lat,&c.chart,c.t0,c.extras,c.params,c.strata,c.anchor)?;lp.publish_margin_posts=256;
        let result=crate::t1b_rivernet::evaluate_owned(&lp,&c.net,&c.curves,c.t0,&c.extras.laws,&c.table,c.params,&c.coast,Some(self),Some(key.stage))?;
        let capture=match result{Evaluated::Captured(c)=>c,Evaluated::Surface(_)=>return Err("canonical source missed requested stage".into())};
        let n=(TILE*TILE)as usize;let mut z=Vec::with_capacity(n);let mut fixed=Vec::with_capacity(n);let mut links=Vec::with_capacity(n);
        for y in 0..TILE as usize{for x in 0..TILE as usize{
            let gx=key.tx*TILE+x as i64;let gy=key.ty*TILE+y as i64;
            if gx.abs()>half||gy.abs()>half{z.push(f64::NAN);fixed.push(crate::t1b::FREE);links.push(crate::NONE);continue;}
            let i=lat.idx(gx,gy).ok_or("canonical owner core outside its evaluation")?;z.push(capture.z[i]);fixed.push(capture.fixed[i]);links.push(capture.links[i]);
        }}
        let mut b=Vec::with_capacity(n*13);for v in &z{b.extend_from_slice(&v.to_le_bytes());}b.extend_from_slice(&fixed);for v in &links{b.extend_from_slice(&v.to_le_bytes());}
        let meta=TileMeta{key,source_pin:self.pin.clone(),digest:blake3::hash(&b).to_hex().to_string(),dependencies:capture.dependencies,rule:NATIVE_RULE.into()};
        use std::io::Write;
        std::fs::OpenOptions::new().create_new(true).write(true).open(&data_path).and_then(|mut f|f.write_all(&b)).map_err(|e|e.to_string())?;
        std::fs::OpenOptions::new().create_new(true).write(true).open(&meta_path).and_then(|mut f|f.write_all(&serde_json::to_vec(&meta).unwrap())).map_err(|e|e.to_string())?;
        Ok(NativeTile{meta,z,fixed,links})
    }
}
impl Conditioner for NativeOwner<'_,'_>{
    fn condition(&self,stage:Stage,lat:&EvalLat,z:&mut[f64],fixed:&mut[u8])->R<Proof>{
        if lat.dx_mm!=crate::t1b::T1B_SPACING_MM{return Err("canonical construction currently requires native 32m posts".into());}
        let x0=lat.gx0.div_euclid(TILE);let y0=lat.gy0.div_euclid(TILE);let x1=(lat.gx0+lat.n as i64-1).div_euclid(TILE);let y1=(lat.gy0+lat.n as i64-1).div_euclid(TILE);
        let mut keys=BTreeSet::new();for ty in y0..=y1{for tx in x0..=x1{keys.insert(TileKey{stage,ty,tx});}}
        for iteration in 0..32{
            let minx=keys.iter().map(|k|k.tx).min().unwrap();let maxx=keys.iter().map(|k|k.tx).max().unwrap();let miny=keys.iter().map(|k|k.ty).min().unwrap();let maxy=keys.iter().map(|k|k.ty).max().unwrap();
            let rect=Rect{gx0:minx*TILE,gy0:miny*TILE,width:((maxx-minx+1)*TILE)as usize,height:((maxy-miny+1)*TILE)as usize};
            let mut working=WorkingSet::new(rect,format!("{}:{}",self.pin,stage.name()),self.max_working_posts).map_err(|e|e.to_string())?;
            working.native_chart=Some(self.pin.clone());
            let mut read_set=Vec::new();
            for &key in &keys{
                let tile=self.tile(key)?;read_set.push(tile.meta.digest.clone());
                // Tile identity and its on-disk content have already been checked.
                // Core ownership is half-open; no shared post is counted twice.
                for y in 0..TILE as usize{for x in 0..TILE as usize{let i=y*TILE as usize+x;let id=PostId{gy:key.ty*TILE+y as i64,gx:key.tx*TILE+x as i64};let j=rect.index(id).unwrap();let half=vz_tilekey::CHART_HALF_EXTENT_MM/crate::t1b::T1B_SPACING_MM;if id.gx.abs()>half||id.gy.abs()>half{continue;}working.present[j]=true;working.ground[j]=tile.z[i];working.fixed[j]=tile.fixed[i];}}
            }
            let p=working.patch();let bounds=crate::t1b_ls::condition_detail_bounds(&p).map_err(|e|e.to_string())?;
            let targets:Vec<usize>=(0..lat.nodes()).map(|i|{let(gx,gy)=lat.g(i);rect.index(PostId{gx,gy}).unwrap()}).collect();
            let cert=bounds.certify(&p,&targets,lat.dx()*lat.dx()).map_err(|e|e.to_string())?;
            if cert.pass{
                let mut adoption=0.0f64;let mut max_adoption=0.0f64;let mut out=blake3::Hasher::new();out.update(self.pin.as_bytes());out.update(stage.name().as_bytes());
                for (i,&j) in targets.iter().enumerate(){
                    if (fixed[i]==crate::t1b::CHANNEL)!=(working.fixed[j]==crate::t1b::CHANNEL){return Err(format!("canonical stage would change immutable T1a classification at {:?}",rect.id(j)));}
                    let d=working.ground[j]-z[i];adoption+=d*lat.dx()*lat.dx();max_adoption=max_adoption.max(d.abs());z[i]=bounds.upper[j];fixed[i]=working.fixed[j];out.update(&z[i].to_le_bytes());out.update(&[fixed[i]]);
                }
                if !adoption.is_finite(){return Err("nonfinite construction input adoption".into());}
                return Ok(Proof{stage,source_pin:self.pin.clone(),input_digest:bounds.input_digest,output_digest:out.finalize().to_hex().to_string(),requested_posts:targets.len(),fetched_tiles:keys.len(),loaded_posts:working.present.iter().filter(|x|**x).count(),completion_iterations:iteration,canonical_raised_posts:cert.raised_posts as u32,initial_raise_m3:cert.initial_raise_m3,max_raise_m:cert.max_raise_m,local_input_adoption_net_m3:adoption,local_input_adoption_max_m:max_adoption,read_set,rule:NATIVE_RULE.into()});
            }
            let before=keys.len();for id in cert.missing{keys.insert(TileKey{stage,ty:id.gy.div_euclid(TILE),tx:id.gx.div_euclid(TILE)});}
            if keys.len()==before{return Err(format!("NeedsBoundary: no new native tile for {stage:?} at [{},{},{}]; chart contract or missing interior input must be reconciled",lat.gx0,lat.gy0,lat.n));}
        }
        Err("NeedsDomain: canonical completion exceeded 32 fetch rounds; no publication".into())
    }
}
