//! Optional observation of the native donor-first fluvial pass, not material ancestry.
use pf_contracts::params::FluvialTraceBlock;
use pf_grid::GridSpec;
use serde::{Deserialize,Serialize};
type Result<T> = std::result::Result<T,String>;
pub const HEADER_BYTES:u64=65536;
pub const ROW_BYTES:u64=512;
pub const TOPOLOGY_BYTES_PER_CELL:u64=36;
pub const FRAME_BYTES:u64=1024;
const MAX_CELLS:usize=4_194_304;
pub(super) mod bits {
    use serde::{Deserialize,Deserializer,Serializer};
    pub fn serialize<S:Serializer>(value:&f64,s:S)->Result<S::Ok,S::Error>{
        if !value.is_finite(){return Err(serde::ser::Error::custom("nonfinite trace number"))}
        s.serialize_str(&format!("{:016x}",value.to_bits()))
    }
    pub fn deserialize<'de,D:Deserializer<'de>>(d:D)->Result<f64,D::Error>{
        let s=String::deserialize(d)?;
        if s.len()!=16 || !s.bytes().all(|c|c.is_ascii_digit()||(b'a'..=b'f').contains(&c)){return Err(serde::de::Error::custom("invalid trace binary64 bits"))}
        let n=f64::from_bits(u64::from_str_radix(&s,16).map_err(serde::de::Error::custom)?);
        if !n.is_finite(){return Err(serde::de::Error::custom("nonfinite trace number"))}Ok(n)
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TraceClock {
    pub step:u64,
    #[serde(with="bits")] pub start_my:f64,
    #[serde(with="bits")] pub end_my:f64,
}
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FluvialRow {
    pub cell:u32,pub receiver:Option<u32>,pub continental:bool,pub grade_evaluated:bool,
    #[serde(with="bits")] pub incoming_m3:f64,
    #[serde(with="bits")] pub detached_m3:f64,
    #[serde(with="bits")] pub grade_need_m3:f64,
    #[serde(with="bits")] pub retained_m3:f64,
    #[serde(with="bits")] pub outgoing_m3:f64,
}
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FluvialFrame {
    pub clock:TraceClock,
    pub wet:Vec<bool>,pub ocean:Vec<bool>,pub receivers:Vec<Option<u32>>,pub donor_order:Vec<u32>,
    pub rows:Vec<FluvialRow>,
    #[serde(default,skip_serializing_if="Option::is_none")]
    pub surface_budget:Option<super::surface_budget::SurfaceBudget>,
}
#[derive(Debug,Clone,PartialEq,Default,Serialize)]
pub struct FluvialSummary {
    pub detached_m3:f64,pub retained_m3:f64,pub marine_delivery_m3:f64,pub lake_delivery_m3:f64,
    pub terminal_pending_m3:f64,pub closure_residual_m3:f64,pub rounding_budget_m3:f64,
    pub land_edges:u64,pub water_edges:u64,pub terminals:u64,
}
impl FluvialFrame {
    pub fn new(clock:TraceClock,wet:Vec<bool>,ocean:Vec<bool>,receivers:Vec<Option<u32>>,donor_order:Vec<u32>)->Self {
        let rows=Vec::with_capacity(wet.iter().filter(|&&w|!w).count());
        Self{clock,wet,ocean,receivers,donor_order,rows,surface_budget:None}
    }
    pub fn validate(&self,cells:usize,dt:f64)->Result<FluvialSummary>{
        if cells==0||cells>MAX_CELLS||self.wet.len()!=cells||self.ocean.len()!=cells||self.receivers.len()!=cells
            ||self.donor_order.len()!=cells||self.rows.len()>cells{return Err("fluvial trace topology shape".into())}
        let c=self.clock;
        let eps=32.*f64::EPSILON*c.end_my.abs().max(dt).max(1.);
        if c.step==0||!c.start_my.is_finite()||!c.end_my.is_finite()||!dt.is_finite()||dt<=0.
            ||(c.start_my-(c.step-1) as f64*dt).abs()>eps||(c.end_my-c.step as f64*dt).abs()>eps {
            return Err("fluvial trace clock mismatch".into());
        }
        let mut order_index=vec![usize::MAX;cells];
        for (rank,&i) in self.donor_order.iter().enumerate(){let i=i as usize;if i>=cells||order_index[i]!=usize::MAX{return Err("fluvial trace order is not a permutation".into())}order_index[i]=rank;}
        for i in 0..cells {
            if self.ocean[i]&&!self.wet[i]{return Err("ocean trace cell must be wet".into())}
            if let Some(j)=self.receivers[i]{if j as usize>=cells||j as usize==i{return Err("invalid fluvial receiver".into())}}
        }
        let expected:Vec<_>=self.donor_order.iter().copied().filter(|&i|!self.wet[i as usize]).collect();
        if expected.len()!=self.rows.len()||expected.iter().zip(&self.rows).any(|(i,r)|*i!=r.cell){return Err("missing or reordered native land-cell trace".into())}
        let mut incoming=vec![0.0f64;cells];let mut summary=FluvialSummary::default();
        for r in &self.rows {
            let i=r.cell as usize;
            if [r.incoming_m3,r.detached_m3,r.grade_need_m3,r.retained_m3,r.outgoing_m3].iter().any(|x|!x.is_finite()||*x<0.) {
                return Err("invalid fluvial trace volume".into());
            }
            if r.receiver!=self.receivers[i]||r.incoming_m3.to_bits()!=incoming[i].to_bits(){return Err("native fluvial receiver/incoming replay differs".into())}
            if !r.grade_evaluated&&(r.grade_need_m3!=0.||r.retained_m3!=0.){return Err("retention recorded outside evaluated branch".into())}
            if r.grade_evaluated&&(!r.continental||r.receiver.is_none()||r.incoming_m3<=0.){return Err("invalid evaluated grade branch".into())}
            let take=if r.grade_evaluated{r.incoming_m3.min(r.grade_need_m3)}else{0.};
            if take.to_bits()!=r.retained_m3.to_bits(){return Err("native graded retention replay differs".into())}
            let mut outgoing=incoming[i];if take>0.{outgoing-=take;}outgoing+=r.detached_m3;
            if outgoing.to_bits()!=r.outgoing_m3.to_bits(){return Err("native fluvial outgoing replay differs".into())}
            summary.detached_m3+=r.detached_m3;summary.retained_m3+=take;
            match r.receiver {
                None=>{summary.terminals+=1;summary.terminal_pending_m3+=outgoing;},
                Some(j)=>{let j=j as usize;
                    if self.wet[j]{summary.water_edges+=1;if self.ocean[j]{summary.marine_delivery_m3+=outgoing}else{summary.lake_delivery_m3+=outgoing}}
                    else {if order_index[j]<=order_index[i]{return Err("fluvial edge violates donor-first order".into())}
                        summary.land_edges+=1;if outgoing>0.{incoming[j]+=outgoing;}}
                }
            }
        }
        summary.closure_residual_m3=summary.detached_m3-summary.retained_m3-summary.marine_delivery_m3-summary.lake_delivery_m3-summary.terminal_pending_m3;
        summary.rounding_budget_m3=1.+1e-9*summary.detached_m3.abs().max(1.);
        if !summary.closure_residual_m3.is_finite()||summary.closure_residual_m3.abs()>summary.rounding_budget_m3{return Err("native fluvial pass volume closure".into())}
        Ok(summary)
    }
}
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
struct ArchiveWire {
    schema:String,scope:String,float_encoding:String,sources:[String;4],
    config:FluvialTraceBlock,grid:[usize;2],
    #[serde(with="bits")] dt_my:f64,
    frames:Vec<FluvialFrame>,rows_used:u64,charged_bytes:u64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct FluvialArchive {wire:ArchiveWire}
fn sources()->[String;4]{
    use pf_contracts::sim_module::digest_bytes;
    [digest_bytes(include_bytes!("retain.rs")),digest_bytes(include_bytes!("../erosion.rs")),digest_bytes(include_bytes!("fluvial_trace.rs")),digest_bytes(include_bytes!("surface_budget.rs"))]
}
impl Serialize for FluvialArchive {fn serialize<S:serde::Serializer>(&self,s:S)->std::result::Result<S::Ok,S::Error>{self.wire.serialize(s)}}
impl<'de> Deserialize<'de> for FluvialArchive {
    fn deserialize<D:serde::Deserializer<'de>>(d:D)->std::result::Result<Self,D::Error>{
        let archive=Self{wire:ArchiveWire::deserialize(d)?};archive.validate_all().map_err(serde::de::Error::custom)?;Ok(archive)
    }
}
impl FluvialArchive {
    pub fn validate_config(config:&FluvialTraceBlock)->Result<()> {
        if config.steps.is_empty()||config.steps.len()>32||config.steps[0]==0||config.steps.windows(2).any(|w|w[0]>=w[1])
            ||!(1..=10_000_000).contains(&config.record_limit)||!(HEADER_BYTES..=1024*1024*1024).contains(&config.byte_limit){return Err("fluvial trace needs 1..32 sorted positive steps, 1..10000000 records and 64KiB..1GiB bytes".into())}Ok(())
    }
    pub fn new(config:FluvialTraceBlock,grid:&GridSpec,dt_my:f64)->Result<Self>{
        Self::validate_config(&config)?;
        if grid.nlat().checked_mul(grid.nlon()).is_none_or(|n|n==0||n>MAX_CELLS)||!dt_my.is_finite()||dt_my<=0.{return Err("invalid fluvial trace grid/clock".into())}
        Ok(Self{wire:ArchiveWire{schema:"vz.native-fluvial-trace/v2".into(),scope:"native_donor_first_pass_not_material_lineage".into(),float_encoding:"ieee754-binary64-lowercase-hex".into(),sources:sources(),config,grid:[grid.nlat(),grid.nlon()],dt_my,frames:Vec::new(),rows_used:0,charged_bytes:HEADER_BYTES}})
    }
    pub fn config(&self)->&FluvialTraceBlock{&self.wire.config}
    pub fn grid(&self)->[usize;2]{self.wire.grid}
    pub fn dt_my(&self)->f64{self.wire.dt_my}
    pub fn frames(&self)->&[FluvialFrame]{&self.wire.frames}
    pub fn charged_bytes(&self)->u64{self.wire.charged_bytes}
    pub fn due(&self,step:u64)->bool{self.wire.config.steps.binary_search(&step).is_ok()}
    fn stage_bytes(&self)->u64 {if self.config().stage_budgets {super::surface_budget::BASE_BYTES+super::surface_budget::BYTES_PER_CELL*(self.wire.grid[0]*self.wire.grid[1]) as u64}else{0}}
    fn validate_surface(&self,frame:&FluvialFrame)->Result<()> {
        if frame.surface_budget.is_some()!=self.config().stage_budgets{return Err("surface budget profile mismatch".into())}
        if let Some(b)=&frame.surface_budget{b.validate(frame,self.wire.grid)?;}Ok(())
    }
    pub fn step_reservation_bytes(&self)->u64{self.stage_bytes()+FRAME_BYTES+(ROW_BYTES+TOPOLOGY_BYTES_PER_CELL)*(self.wire.grid[0]*self.wire.grid[1]) as u64}
    pub fn preflight(&self,step:u64)->Result<()> {
        if !self.due(step){return Ok(())}
        let cells=(self.wire.grid[0]*self.wire.grid[1]) as u64;
        if self.wire.rows_used+cells>self.wire.config.record_limit||self.wire.charged_bytes+self.step_reservation_bytes()>self.wire.config.byte_limit{return Err("fluvial trace budget exhausted before mutation".into())}Ok(())
    }
    pub fn reserve(&mut self,step:u64)->Result<()> {
        self.preflight(step)?;if self.due(step){self.wire.frames.try_reserve(1).map_err(|_|"fluvial frame allocation")?;}Ok(())
    }
    pub fn record(&mut self,frame:FluvialFrame)->Result<()> {
        self.preflight(frame.clock.step)?;
        if !self.due(frame.clock.step)||self.wire.frames.last().is_some_and(|f|f.clock.step>=frame.clock.step){return Err("unexpected fluvial frame".into())}
        let cells=self.wire.grid[0]*self.wire.grid[1];frame.validate(cells,self.wire.dt_my)?;self.validate_surface(&frame)?;
        self.wire.rows_used+=frame.rows.len() as u64;
        self.wire.charged_bytes+=self.stage_bytes()+FRAME_BYTES+TOPOLOGY_BYTES_PER_CELL*cells as u64+ROW_BYTES*frame.rows.len() as u64;
        self.wire.frames.push(frame);Ok(())
    }
    pub fn check_grid(&self,grid:&GridSpec)->Result<()> {
        if self.grid()!=[grid.nlat(),grid.nlon()]||self.frames().iter().any(|f|f.surface_budget.as_ref().is_some_and(|b|b.radius_m.to_bits()!=grid.radius_m().to_bits())){return Err("surface/fluvial grid profile mismatch".into())}Ok(())
    }
    pub fn check_schedule(&self,done:u64)->Result<()> {
        if !self.wire.config.steps.iter().copied().filter(|s|*s<=done).eq(self.wire.frames.iter().map(|f|f.clock.step)) {return Err("missing or unexpected selected fluvial step".into())}Ok(())
    }
    pub fn validate_all(&self)->Result<()> {
        Self::validate_config(&self.wire.config)?;
        let cells=self.wire.grid[0].checked_mul(self.wire.grid[1]).filter(|n|*n>0&&*n<=MAX_CELLS).ok_or("fluvial grid bound")?;
        if !self.wire.dt_my.is_finite()||self.wire.dt_my<=0.||self.wire.schema!="vz.native-fluvial-trace/v2"||self.wire.scope!="native_donor_first_pass_not_material_lineage"
            ||self.wire.float_encoding!="ieee754-binary64-lowercase-hex"||self.wire.sources!=sources(){return Err("fluvial trace schema/source/clock mismatch".into())}
        let mut rows=0u64;let mut charge=HEADER_BYTES;
        for frame in &self.wire.frames {frame.validate(cells,self.wire.dt_my)?;self.validate_surface(&frame)?;rows+=frame.rows.len() as u64;charge+=self.stage_bytes()+FRAME_BYTES+TOPOLOGY_BYTES_PER_CELL*cells as u64+ROW_BYTES*frame.rows.len() as u64;}
        if rows!=self.wire.rows_used||charge!=self.wire.charged_bytes||rows>self.wire.config.record_limit||charge>self.wire.config.byte_limit{return Err("fluvial trace budget accounting".into())}
        self.check_schedule(self.wire.frames.last().map_or(0,|f|f.clock.step))?;Ok(())
    }
}
