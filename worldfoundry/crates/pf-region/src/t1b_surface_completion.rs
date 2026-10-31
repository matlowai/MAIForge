//! Selective dependency certificates for the LS owner. A frontier is an
//! optimistic *bound*, never a physical outlet. Known-native-outlet minimax H+
//! is an upper bound. Seeding missing-domain neighbours at their own ground
//! gives H- <= true H <= H+. Equality certifies the requested height exactly.
//! Routing/area needs the additional donor/capture/flat certificate below.
use crate::t1b_surface_domain::{invalid, DomainError, PostId, Result};
use serde::{Deserialize,Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeSet,BinaryHeap};

pub const RULE:&str="ls_selective_bounds_v1: exact native-outlet upper minimax and optimistic missing-frontier lower minimax; equal finite bounds certify heights; frontier seeds exist only in the lower-bound solve and never become physical terminals; actual routing separately requires complete flat/donor/capture dependencies";
const NONE:u32=u32::MAX;

#[derive(Clone,Copy,Debug,Serialize,Deserialize,PartialEq,Eq)]
pub struct Rect {pub gx0:i64,pub gy0:i64,pub width:usize,pub height:usize}
impl Rect {
    pub fn count(self)->Result<usize>{self.width.checked_mul(self.height).filter(|&n|n>0 && n<u32::MAX as usize).ok_or_else(||invalid("invalid detail rectangle"))}
    pub fn id(self,i:usize)->PostId{PostId{gy:self.gy0+(i/self.width) as i64,gx:self.gx0+(i%self.width) as i64}}
    pub fn index(self,id:PostId)->Option<usize>{let(x,y)=(id.gx-self.gx0,id.gy-self.gy0);if x<0||y<0||x>=self.width as i64||y>=self.height as i64{None}else{Some(y as usize*self.width+x as usize)}}
    pub fn neighbour(self,i:usize,k:usize)->Option<usize>{self.index(self.id(i).neighbour(k))}
}
/// Missing posts have present=false; they are not ocean, void or boundary data.
/// The immutable pin names the supplied predecessor surface, not a locator.
pub struct Patch<'a>{pub rect:Rect,pub ground:&'a[f64],pub fixed:&'a[u8],pub present:&'a[bool],pub predecessor:&'a str, pub native_chart:Option<&'a str>}
impl Patch<'_>{
    fn valid(&self)->Result<usize>{
        let n=self.rect.count()?;
        if self.ground.len()!=n||self.fixed.len()!=n||(!self.present.is_empty()&&self.present.len()!=n)||self.predecessor.is_empty(){return Err(invalid("patch shape or predecessor pin"));}
        // Protect global-neighbour arithmetic and native chart-sized indexing.
        if self.rect.gx0.abs_diff(0)>1_000_000_000||self.rect.gy0.abs_diff(0)>1_000_000_000||self.rect.width>1_000_000||self.rect.height>1_000_000{return Err(invalid("patch address range"));}
        if self.native_chart.is_some_and(|s|s.is_empty()){return Err(invalid("missing native chart contract pin"));}
        for i in 0..n{if self.has(i)&&self.outside_native_chart(self.rect.id(i)){return Err(invalid("known sample outside native chart contract"));}if self.has(i)&&(!self.ground[i].is_finite()||![crate::t1b::FREE,crate::t1b::CHANNEL,crate::t1b::GUIDE,crate::t1b::MARINE,crate::t1b::LAKE].contains(&self.fixed[i])){return Err(invalid(format!("nonfinite ground or unsupported boundary at {:?}",self.rect.id(i))));}}
        Ok(n)
    }
    pub fn has(&self,i:usize)->bool{self.present.is_empty()||self.present[i]}
    fn outside_native_chart(&self,id:PostId)->bool{let half=vz_tilekey::CHART_HALF_EXTENT_MM/crate::t1b::T1B_SPACING_MM;self.native_chart.is_some()&&(id.gx.abs()>half||id.gy.abs()>half)}
    fn frontier(&self,i:usize)->bool{(0..8).any(|k|{let id=self.rect.id(i).neighbour(k);!self.outside_native_chart(id)&&self.rect.index(id).map(|j|!self.has(j)).unwrap_or(true)})}
}
#[derive(Clone,Copy)]
struct Key{z:f64,id:PostId,i:u32}
impl PartialEq for Key{fn eq(&self,o:&Self)->bool{self.cmp(o)==Ordering::Equal}}
impl Eq for Key{}
impl PartialOrd for Key{fn partial_cmp(&self,o:&Self)->Option<Ordering>{Some(self.cmp(o))}}
impl Ord for Key{fn cmp(&self,o:&Self)->Ordering{o.z.total_cmp(&self.z).then(o.id.cmp(&self.id))}}

pub struct Bounds {
    pub upper:Vec<f64>,pub lower:Vec<f64>,
    /// Actual known-outlet parent, never the optimistic-frontier parent.
    pub parent:Vec<u32>,pub order:Vec<u32>,
    /// Global-index frontier representative for a possibly lower escape.
    lower_frontier:Vec<u32>,
    pub input_digest:String,
    rect:Rect, predecessor:String,
}
#[derive(Debug,Serialize,Deserialize)]
pub struct FillCertificate{
    pub rule:String,pub input_digest:String,pub predecessor:String,
    pub requested_posts:usize,pub exact_posts:usize,pub unresolved_posts:usize,
    pub raised_posts:usize,pub initial_raise_m3:f64,pub max_raise_m:f64,
    pub missing:Vec<PostId>,pub pass:bool,
}

fn flood(p:&Patch<'_>,optimistic:bool)->(Vec<f64>,Vec<u32>,Vec<u32>,Vec<u32>){
    let n=p.ground.len();let mut h=vec![f64::INFINITY;n];let mut parent=vec![NONE;n];let mut origin=vec![NONE;n];let mut seen=vec![false;n];let mut heap=BinaryHeap::new();let mut order=Vec::new();
    for i in 0..n{if p.has(i)&&(p.fixed[i]!=crate::t1b::FREE||(optimistic&&p.frontier(i))){seen[i]=true;h[i]=p.ground[i];if p.fixed[i]==crate::t1b::FREE{origin[i]=i as u32;}heap.push(Key{z:h[i],id:p.rect.id(i),i:i as u32});}}
    while let Some(Key{z,i,..})=heap.pop(){let i=i as usize;order.push(i as u32);
        for k in 0..8{if let Some(j)=p.rect.neighbour(i,k){if p.has(j)&&!seen[j]{seen[j]=true;h[j]=p.ground[j].max(z);parent[j]=i as u32;origin[j]=origin[i];heap.push(Key{z:h[j],id:p.rect.id(j),i:j as u32});}}}
    }
    (h,parent,origin,order)
}
/// Owner-only implementation; use t1b_ls::condition_detail_bounds.
pub(crate) fn solve_bounds(p:&Patch<'_>)->Result<Bounds>{
    p.valid()?;
    let input_digest=patch_digest(p)?;
    let (upper,parent,_,order)=flood(p,false);let(lower,_,lower_frontier,_)=flood(p,true);
    Ok(Bounds{upper,lower,parent,order,lower_frontier,input_digest,rect:p.rect,predecessor:p.predecessor.into()})
}
fn patch_digest(p:&Patch<'_>)->Result<String>{
    p.valid()?;
    let mut hash=blake3::Hasher::new();hash.update(RULE.as_bytes());hash.update(p.predecessor.as_bytes());if let Some(pin)=p.native_chart{hash.update(b"native_chart_d8_no_exterior_land_edges_v1");hash.update(pin.as_bytes());}hash.update(&serde_json::to_vec(&p.rect).map_err(|e|invalid(e.to_string()))?);
    for i in 0..p.ground.len(){hash.update(&[p.has(i) as u8]);if p.has(i){hash.update(&p.ground[i].to_bits().to_le_bytes());hash.update(&[p.fixed[i]]);}}
    Ok(hash.finalize().to_hex().to_string())
}
impl Bounds{
    pub fn height_exact(&self,i:usize)->bool{self.upper[i].is_finite()&&self.upper[i]==self.lower[i]}
    fn add_frontier(p:&Patch<'_>,root:usize,missing:&mut BTreeSet<PostId>){for k in 0..8{let id=p.rect.id(root).neighbour(k);if !p.outside_native_chart(id)&&p.rect.index(id).map(|j|!p.has(j)).unwrap_or(true){missing.insert(id);}}}
    pub fn certify(&self,p:&Patch<'_>,targets:&[usize],support_m2:f64)->Result<FillCertificate>{
        if self.rect!=p.rect||self.predecessor!=p.predecessor||self.input_digest!=patch_digest(p)?{return Err(DomainError::Stale{field:"fill proof input".into()});}
        if !support_m2.is_finite()||support_m2<=0.0{return Err(invalid("positive finite support required"));}
        let mut c=FillCertificate{rule:RULE.into(),input_digest:self.input_digest.clone(),predecessor:p.predecessor.into(),requested_posts:targets.len(),exact_posts:0,unresolved_posts:0,raised_posts:0,initial_raise_m3:0.0,max_raise_m:0.0,missing:Vec::new(),pass:false};let mut missing=BTreeSet::new();let mut unique=BTreeSet::new();
        for &i in targets{
            if i>=p.ground.len()||!p.has(i)||!unique.insert(i){return Err(invalid("missing or duplicate requested post"));}
            if self.height_exact(i){c.exact_posts+=1;let dh=self.upper[i]-p.ground[i];if dh>0.0{c.raised_posts+=1;c.initial_raise_m3+=dh*support_m2;c.max_raise_m=c.max_raise_m.max(dh);}}
            else {c.unresolved_posts+=1;let root=self.lower_frontier[i];if root==NONE{return Err(invalid("unreachable component has no supported outlet or fetchable frontier"));}Self::add_frontier(p,root as usize,&mut missing);}
        }
        if !c.initial_raise_m3.is_finite(){return Err(invalid("nonfinite initial-condition raise volume"));}
        c.missing=missing.into_iter().collect();c.pass=c.unresolved_posts==0;Ok(c)
    }
    pub fn require(&self,p:&Patch<'_>,targets:&[usize],support_m2:f64,stage:&str)->Result<FillCertificate>{let c=self.certify(p,targets,support_m2)?;if c.pass{Ok(c)}else{Err(DomainError::NeedsDomain{stage:stage.into(),frontier:c.missing})}}
    /// Reconstruct a known physical escape witness. Optimistic seeds cannot enter it.
    pub fn path(&self,p:&Patch<'_>,start:usize)->Result<Vec<PostId>>{
        if start>=p.ground.len()||!self.upper[start].is_finite(){return Err(invalid("no known outlet path"));}
        let mut path=Vec::new();let mut i=start;
        loop{path.push(p.rect.id(i));if p.fixed[i]!=crate::t1b::FREE{return Ok(path);}let j=self.parent[i];if j==NONE||path.len()>p.ground.len(){return Err(invalid("broken known-outlet path"));}i=j as usize;}
    }
}

#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Sample{pub id:PostId,pub ground_m:f64,pub fixed:u8}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Fetch{pub predecessor:String,pub content_digest:String,pub samples:Vec<Sample>}
#[derive(Debug,Serialize,Deserialize)]
pub struct CompletionReport{pub iterations:usize,pub loaded_posts:usize,pub fetched_posts:usize,pub read_set:Vec<String>,pub certificate:FillCertificate}
/// A bounded working rectangle with explicit holes, not an invented rectangular
/// hydrological boundary. Fetches come from the canonical predecessor producer.
pub struct WorkingSet{pub rect:Rect,pub ground:Vec<f64>,pub fixed:Vec<u8>,pub present:Vec<bool>,pub predecessor:String,pub native_chart:Option<String>,pub read_set:BTreeSet<String>}
impl WorkingSet{
    pub fn new(rect:Rect,predecessor:String,max_posts:usize)->Result<Self>{let n=rect.count()?;if n>max_posts{return Err(DomainError::Resource{posts:n,limit:max_posts});}Ok(Self{rect,ground:vec![f64::NAN;n],fixed:vec![crate::t1b::FREE;n],present:vec![false;n],predecessor,native_chart:None,read_set:BTreeSet::new()})}
    pub fn ingest(&mut self,fetch:Fetch)->Result<usize>{
        if fetch.predecessor!=self.predecessor{return Err(DomainError::Stale{field:"fetched predecessor".into()});}
        if crate::t1b_surface_domain::content_digest(&fetch.samples)?!=fetch.content_digest{return Err(DomainError::Stale{field:"fetched content".into()});}
        let mut seen=BTreeSet::new();let mut added=0;
        // Validate the complete fetch before mutating the cache.
        for p in &fetch.samples{
            if !seen.insert(p.id){return Err(DomainError::DuplicateSource{post:p.id});}
            let i=self.rect.index(p.id).ok_or_else(||DomainError::NeedsDomain{stage:"working rectangle expansion".into(),frontier:vec![p.id]})?;
            if !p.ground_m.is_finite()||![0,1,2,4,5].contains(&p.fixed){return Err(invalid("invalid fetched ground or physical boundary"));}
            if self.present[i]&&(self.ground[i].to_bits()!=p.ground_m.to_bits()||self.fixed[i]!=p.fixed){return Err(DomainError::Stale{field:format!("canonical input disagreement at {:?}",p.id)});}
        }
        for p in fetch.samples{let i=self.rect.index(p.id).unwrap();if !self.present[i]{added+=1;self.present[i]=true;self.ground[i]=p.ground_m;self.fixed[i]=p.fixed;}}
        self.read_set.insert(fetch.content_digest);Ok(added)
    }
    pub fn patch(&self)->Patch<'_>{Patch{rect:self.rect,ground:&self.ground,fixed:&self.fixed,present:&self.present,predecessor:&self.predecessor,native_chart:self.native_chart.as_deref()}}
    /// Recompute bounds after each exact predecessor fetch. Never fill an
    /// uncertain output or advance an accepted pointer. A finite iteration or
    /// working-rectangle cap returns NeedsDomain with the actual next request.
    pub fn complete(&mut self,targets:&[usize],support_m2:f64,max_iterations:usize,fetch:&mut dyn FnMut(&[PostId])->Result<Fetch>)->Result<(Bounds,CompletionReport)>{
        let mut fetched=0;
        for iteration in 0..=max_iterations{
            let p=self.patch();let bounds=crate::t1b_ls::condition_detail_bounds(&p)?;let c=bounds.certify(&p,targets,support_m2)?;
            if c.pass{return Ok((bounds,CompletionReport{iterations:iteration,loaded_posts:self.present.iter().filter(|x|**x).count(),fetched_posts:fetched,read_set:self.read_set.iter().cloned().collect(),certificate:c}));}
            if iteration==max_iterations{return Err(DomainError::NeedsDomain{stage:"selective fill completion".into(),frontier:c.missing});}
            let packet=fetch(&c.missing)?;let n=self.ingest(packet)?;
            if n==0{return Err(invalid("dependency fetch made no progress"));}fetched+=n;
        }
        unreachable!()
    }
}
