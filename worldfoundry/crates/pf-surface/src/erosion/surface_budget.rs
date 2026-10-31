//! Cumulative native stage counters; no inferred source mixtures.
use super::{ErosionStep,Lowering};
use super::fluvial_trace::{bits,FluvialFrame};
use pf_grid::GridSpec;
use serde::{Deserialize,Serialize};
type Result<T> = std::result::Result<T,String>;
pub const STAGES:[&str;6]=["fluvial_detachment","graded_retention","lake_floor_placement","continental_shelf","legacy_marine","hillslope_diffusion"];
// 23 exact arrays (19 JSON bytes/value), plus counter records and headroom.
pub const BYTES_PER_CELL:u64=512;
pub const BASE_BYTES:u64=16384;
mod array_bits {
 use serde::{Deserialize,Deserializer,Serializer,ser::SerializeSeq};
 pub fn serialize<S:Serializer>(v:&[f64],s:S)->Result<S::Ok,S::Error>{
  let mut seq=s.serialize_seq(Some(v.len()))?;for x in v{if !x.is_finite(){return Err(serde::ser::Error::custom("nonfinite surface counter"))}seq.serialize_element(&format!("{:016x}",x.to_bits()))?;}seq.end()
 }
 pub fn deserialize<'de,D:Deserializer<'de>>(d:D)->Result<Vec<f64>,D::Error>{
  #[derive(Deserialize)]struct Number(#[serde(with="super::bits")]f64);
  Ok(Vec::<Number>::deserialize(d)?.into_iter().map(|n|n.0).collect())
 }
}
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Counters {
 #[serde(with="array_bits")]pub bedrock_removed_m:Vec<f64>,
 #[serde(with="array_bits")]pub sediment_removed_m:Vec<f64>,
 #[serde(with="array_bits")]pub sediment_added_m:Vec<f64>,
 #[serde(with="bits")]pub eroded_m3:f64,
 #[serde(with="bits")]pub bedrock_eroded_m3:f64,
 #[serde(with="bits")]pub sediment_eroded_m3:f64,
 #[serde(with="bits")]pub deposited_m3:f64,
 #[serde(with="bits")]pub exported_m3:f64,
}
impl Counters {
 fn capture(s:&Lowering<'_>)->Self{Self{bedrock_removed_m:s.bedrock_removed.clone(),sediment_removed_m:s.sediment_removed.clone(),sediment_added_m:s.sediment_added.clone(),eroded_m3:s.eroded_m3,bedrock_eroded_m3:s.bedrock_eroded_m3,sediment_eroded_m3:s.sediment_eroded_m3,deposited_m3:s.deposited_m3,exported_m3:s.exported_m3}}
 fn final_step(s:&ErosionStep)->Self{Self{bedrock_removed_m:s.bedrock_removed_m.clone(),sediment_removed_m:s.sediment_removed_m.clone(),sediment_added_m:s.sediment_added_m.clone(),eroded_m3:s.eroded_m3,bedrock_eroded_m3:s.bedrock_eroded_m3,sediment_eroded_m3:s.sediment_eroded_m3,deposited_m3:s.deposited_m3,exported_m3:s.exported_m3}}
 fn arrays(&self)->[&[f64];3]{[&self.bedrock_removed_m,&self.sediment_removed_m,&self.sediment_added_m]}
 fn scalars(&self)->[f64;5]{[self.eroded_m3,self.bedrock_eroded_m3,self.sediment_eroded_m3,self.deposited_m3,self.exported_m3]}
 fn validate(&self,areas:&[f64])->Result<()>{
  for a in self.arrays(){if a.len()!=areas.len()||a.iter().any(|v|!v.is_finite()||*v<0.){return Err("invalid surface counter array".into())}}
  if self.scalars().iter().any(|v|!v.is_finite()||*v<0.){return Err("invalid surface counter volume".into())}
  let vol=|a:&[f64]|->f64{a.iter().zip(areas).map(|(x,a)|x*a).sum()};
  close(vol(&self.bedrock_removed_m),self.bedrock_eroded_m3)?;close(vol(&self.sediment_removed_m),self.sediment_eroded_m3)?;
  close(self.bedrock_eroded_m3+self.sediment_eroded_m3,self.eroded_m3)?;
  // exported_m3 is marine receipt already included in sediment_added.
  close(vol(&self.sediment_added_m),self.deposited_m3+self.exported_m3)
 }
}
fn close(a:f64,b:f64)->Result<()>{if !a.is_finite()||!b.is_finite()||(a-b).abs()>1.+1e-9*a.abs().max(b.abs()){Err("surface volume reconciliation".into())}else{Ok(())}}
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Stage{pub name:String,pub cumulative:Counters}
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SurfaceBudget{
 pub schema:String,
 #[serde(with="bits")]pub radius_m:f64,
 #[serde(with="array_bits")]pub areas_m2:Vec<f64>,
 pub stages:Vec<Stage>,pub final_counters:Option<Counters>,
 #[serde(with="array_bits")]pub final_dz_m:Vec<f64>,
}
impl SurfaceBudget{
 pub(super) fn new(grid:&GridSpec,areas:&[f64])->Self{Self{schema:"vz.native-surface-budget/v1".into(),radius_m:grid.radius_m(),areas_m2:areas.to_vec(),stages:Vec::with_capacity(6),final_counters:None,final_dz_m:Vec::new()}}
 pub(super) fn capture(&mut self,s:&Lowering<'_>){self.stages.push(Stage{name:STAGES[self.stages.len()].into(),cumulative:Counters::capture(s)});}
 pub(super) fn finish(&mut self,s:&ErosionStep){self.final_counters=Some(Counters::final_step(s));self.final_dz_m=s.dz_m.clone();}
 pub fn validate(&self,frame:&FluvialFrame,shape:[usize;2])->Result<()>{
  let grid=GridSpec::new(shape[0],shape[1],self.radius_m).map_err(|e|e.to_string())?;let n=grid.cells();
  if self.schema!="vz.native-surface-budget/v1"||self.areas_m2.len()!=n||self.final_dz_m.len()!=n||self.stages.len()!=6{return Err("surface budget shape/schema".into())}
  for i in 0..shape[0]{let area=grid.cell_area_m2(i).map_err(|e|e.to_string())?;for j in 0..shape[1]{if self.areas_m2[i*shape[1]+j].to_bits()!=area.to_bits(){return Err("surface cell area mismatch".into())}}}
  for (k,s) in self.stages.iter().enumerate(){
   if s.name!=STAGES[k]{return Err("surface stage order".into())}s.cumulative.validate(&self.areas_m2)?;
   if k>0{let prev=&self.stages[k-1].cumulative;
    for (a,b) in s.cumulative.arrays().into_iter().zip(prev.arrays()){if a.iter().zip(b).any(|(a,b)|a<b){return Err("nonmonotone cumulative surface counter".into())}}
    if s.cumulative.scalars().into_iter().zip(prev.scalars()).any(|(a,b)|a<b){return Err("nonmonotone cumulative surface volume".into())}
    if k<5&&(s.cumulative.bedrock_removed_m!=prev.bedrock_removed_m||s.cumulative.sediment_removed_m!=prev.sediment_removed_m){return Err("removal outside detachment/diffusion".into())}
   }
  }
  let first=&self.stages[0].cumulative;let graded=&self.stages[1].cumulative;
  if first.sediment_added_m.iter().any(|x|*x!=0.){return Err("deposition before routing".into())}
  if frame.rows.iter().any(|r|r.cell as usize>=n){return Err("surface row index".into())}
  let mut rows=vec![None;n];for r in &frame.rows{rows[r.cell as usize]=Some(r);}
  for i in 0..n{let (detached,retained)=rows[i].map_or((0.,0.),|r|(r.detached_m3,r.retained_m3));close((first.bedrock_removed_m[i]+first.sediment_removed_m[i])*self.areas_m2[i],detached)?;close(graded.sediment_added_m[i]*self.areas_m2[i],retained)?;}
  for k in [4,5]{let c=&self.stages[k].cumulative;close(c.eroded_m3,c.deposited_m3+c.exported_m3)?;}
  let c=self.final_counters.as_ref().ok_or("missing final surface counters")?;if c!=&self.stages[5].cumulative{return Err("final surface arrays differ from last stage".into())}
  for (i,&dz) in self.final_dz_m.iter().enumerate(){let b=c.bedrock_removed_m[i];let s=c.sediment_removed_m[i];let a=c.sediment_added_m[i];if !dz.is_finite()||(dz-(a-b-s)).abs()>1e-7+1e-9*(a+b+s){return Err("surface net elevation reconciliation".into())}}
  Ok(())
 }
}
