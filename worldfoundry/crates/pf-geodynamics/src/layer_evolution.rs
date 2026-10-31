//! Optional observation of actual native layer mutations. Not grain provenance.
use crate::stratigraphy::{StrataStack,StrataParams};
pub const BASEMENT:u32=u32::MAX;
#[derive(Debug,Clone,PartialEq)]
pub struct LayerChange {
    pub kind:u8,
    /// Actual pre/post cover-slot relationships; BASEMENT names the basement.
    /// These are parent relationships, not source mass fractions.
    pub links:Vec<[u32;2]>,
    pub before:StrataStack,
    pub after:StrataStack,
}
#[derive(Debug,Clone,PartialEq,Default)]
pub struct LayerEvolution { pub changes:Vec<LayerChange> }
/// Numbered native mutation kinds; kept source-pinned in archived histories.
pub const NAMES:[&str;22]=[
 "erosion_remove","erosion_truncate","erosion_basement","thin_remainder_absorb",
 "thinning","basement_outflow","belt_basement","pluton_roof","pluton_insert",
 "plume_basement","basement_add","shortening_scale","volcanic_add","sediment_add",
 "basement_reconciliation","basement_underflow","zero_layer_filter","capacity_merge",
 "burial_metamorphism","exposure","top_add","bottom_add",
];
/// Conservative encoded JSON allowance, including both snapshots and all links.
/// apply_parcel starts with <= max_layers cover slots, adds <= four, then removes
/// them in bounded erosion/underflow/capacity loops. Most branches emit one row.
pub fn change_bound(max_layers:usize)->usize {8*max_layers+64}
pub fn byte_bound(max_layers:usize)->u64 {change_bound(max_layers) as u64*(2048+1024*(max_layers as u64+4))}
pub(crate) enum Mapping {
 Identity,
 Insert(usize,bool), // true: transferred from basement; false: unresolved external input
 Remove(usize,Option<usize>), // optional surviving pre-removal cover slot; BASEMENT via BaseAbsorb
 Merge(usize,usize),
 BaseAbsorb(usize),
 BaseTo(usize),
 Filter(Vec<usize>), // surviving pre-filter indices, in native order
}
pub(crate) struct Observer {pub trace:LayerEvolution,last:StrataStack,max_layers:usize, pub provenance:Option<crate::layer_provenance::Recorder>, pub retain_trace:bool}
impl Observer {
 pub fn new(st:&StrataStack,p:&StrataParams,byte_limit:u64)->Result<Self,String>{
  if !(2..=256).contains(&p.max_layers)||st.cover.len()>p.max_layers||byte_limit<byte_bound(p.max_layers){return Err("layer evolution budget/shape refused before mutation".into())}
  let mut changes=Vec::new();changes.try_reserve(change_bound(p.max_layers)).map_err(|_|"layer evolution allocation before mutation")?;
  Ok(Self{trace:LayerEvolution{changes},last:st.clone(),max_layers:p.max_layers,provenance:None,retain_trace:true})
 }
 pub fn capture(&mut self,kind:u8,mapping:Mapping,st:&StrataStack,b:f64){
  let mut after=st.clone();after.basement_m=b;
  if self.last==after{return}
  let n=self.last.cover.len();let mut links=vec![[BASEMENT,BASEMENT]];
  for i in 0..n {
   let target=match &mapping {
    Mapping::Identity|Mapping::BaseTo(_)=>Some(i as u32),
    Mapping::Insert(j,_)=>Some((i+usize::from(i>=*j)) as u32),
    Mapping::Remove(j,into)=>if i==*j{into.map(|k|(k-usize::from(k>*j)) as u32)}else{Some((i-usize::from(i>*j)) as u32)},
    Mapping::Merge(keep,drop)=>Some((if i==*drop{*keep}else{i}-usize::from((if i==*drop{*keep}else{i})>*drop)) as u32),
    Mapping::BaseAbsorb(j)=>if i==*j{Some(BASEMENT)}else{Some((i-usize::from(i>*j)) as u32)},
    Mapping::Filter(surviving)=>surviving.iter().position(|x|*x==i).map(|j|j as u32),
   };
   if let Some(j)=target{links.push([i as u32,j]);}
  }
  match mapping {Mapping::Insert(i,true)|Mapping::BaseTo(i)=>links.push([BASEMENT,i as u32]),_=>{}}
  assert!((kind as usize)<NAMES.len()&&after.cover.len()<=self.max_layers+4,"native layer evolution shape invariant");
  assert!(links.iter().all(|[i,j]|(*i==BASEMENT||(*i as usize)<n)&&(*j==BASEMENT||(*j as usize)<after.cover.len())),"native layer evolution links");
  assert!(self.trace.changes.len()<change_bound(self.max_layers),"native layer evolution reservation invariant");
  let change=LayerChange{kind,links,before:self.last.clone(),after:after.clone()};
  if let Some(p)=self.provenance.as_mut(){p.capture(&change);}
  if self.retain_trace {self.trace.changes.push(change);}
  self.last=after;
 }
 pub fn finish(self,st:&StrataStack)->LayerEvolution {assert_eq!(&self.last,st,"unobserved native stack mutation");self.trace}
}
