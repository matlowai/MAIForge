//! SEED-MAGMATISM (lane MG, 2026-10-13): per-cell vent markers and the
//! volcano feature points of `slab-arc-backarc-llsvp-v1`
//! (`unpublished-design-record`).
//!
//! Markers ride with the crust column through the advance's source map
//! (reset at births, like the extras markers). Volcano points are DERIVED
//! (a map projection, never canon): vent cells ordered by last eruption, then
//! magmatic volume, then cell id, thinned greedily to a minimum spacing per
//! kind and deduplicated by origin; plus one point per active plume.
//!
//! Determinism: ascending cell order, explicit sort keys, no hashed maps.

use crate::kinematic::NO_SOURCE;
use pf_grid::GridSpec;

/// Vent kind codes (0 = none).
pub const VENT_NONE: u8 = 0;
pub const VENT_ISLAND_ARC: u8 = 1;
pub const VENT_CONTINENTAL_ARC: u8 = 2;
pub const VENT_HOTSPOT: u8 = 3;
pub const VENT_RIFT: u8 = 4;
pub const VENT_LIP: u8 = 5;

/// Per-cell magmatic markers (all `grid.cells()` long). Times are My,
/// `-1.0` = never.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct VentMarkers {
    /// Kind of the latest eruption onto the column.
    pub kind: Vec<u8>,
    pub first_t_my: Vec<f64>,
    pub last_t_my: Vec<f64>,
    /// Magmatic crust added by the recorded processes (m).
    pub added_m: Vec<f64>,
    /// Cell where the column's first recorded eruption happened.
    pub origin: Vec<u32>,
    /// Arc heat in [0, 1] (arcs.rs).
    pub heat: Vec<f64>,
    /// Latest back-arc extension time (My).
    pub backarc_t_my: Vec<f64>,
    /// Plume-head (LIP) crust added (m).
    pub head_m: Vec<f64>,
    /// 1 once an active plume passed within the plume-pass radius.
    pub plume_pass: Vec<u8>,
    /// SEED-WORLD-COMPOSE (lane CP) `slab-arc-backarc-llsvp-v2` only (empty
    /// otherwise): first back-arc extension time of the current episode
    /// (My, -1 never) and the continental thickness at its start (m).
    pub backarc_first_t_my: Vec<f64>,
    pub backarc_start_h_m: Vec<f64>,
}

impl VentMarkers {
    pub fn new(n: usize) -> Self {
        Self {
            kind: vec![VENT_NONE; n],
            first_t_my: vec![-1.0; n],
            last_t_my: vec![-1.0; n],
            added_m: vec![0.0; n],
            origin: vec![u32::MAX; n],
            heat: vec![0.0; n],
            backarc_t_my: vec![-1.0; n],
            head_m: vec![0.0; n],
            plume_pass: vec![0; n],
            backarc_first_t_my: Vec::new(),
            backarc_start_h_m: Vec::new(),
        }
    }

    /// SEED-WORLD-COMPOSE (lane CP): markers with the v2 back-arc episode
    /// fields allocated.
    pub fn with_episodes(mut self) -> Self {
        let n = self.kind.len();
        if self.backarc_first_t_my.is_empty() {
            self.backarc_first_t_my = vec![-1.0; n];
            self.backarc_start_h_m = vec![0.0; n];
        }
        self
    }

    pub fn len(&self) -> usize {
        self.kind.len()
    }

    pub fn is_empty(&self) -> bool {
        self.kind.is_empty()
    }

    pub fn check(&self, n: usize) -> Result<(), String> {
        let ok = self.kind.len() == n
            && self.first_t_my.len() == n
            && self.last_t_my.len() == n
            && self.added_m.len() == n
            && self.origin.len() == n
            && self.heat.len() == n
            && self.backarc_t_my.len() == n
            && self.head_m.len() == n
            && self.plume_pass.len() == n
            && (self.backarc_first_t_my.is_empty() || self.backarc_first_t_my.len() == n)
            && self.backarc_start_h_m.len() == self.backarc_first_t_my.len();
        if ok {
            Ok(())
        } else {
            Err(format!("vent markers do not have {n} cells"))
        }
    }

    /// The markers carried to the advanced grid (`NO_SOURCE` = birth: empty).
    pub fn carry(&self, source: &[u32]) -> Self {
        let mut out = Self::new(source.len());
        if !self.backarc_first_t_my.is_empty() {
            out = out.with_episodes();
        }
        for (dst, &s) in source.iter().enumerate() {
            if s == NO_SOURCE {
                continue;
            }
            let s = s as usize;
            out.kind[dst] = self.kind[s];
            out.first_t_my[dst] = self.first_t_my[s];
            out.last_t_my[dst] = self.last_t_my[s];
            out.added_m[dst] = self.added_m[s];
            out.origin[dst] = self.origin[s];
            out.heat[dst] = self.heat[s];
            out.backarc_t_my[dst] = self.backarc_t_my[s];
            out.head_m[dst] = self.head_m[s];
            out.plume_pass[dst] = self.plume_pass[s];
            if !self.backarc_first_t_my.is_empty() {
                out.backarc_first_t_my[dst] = self.backarc_first_t_my[s];
                out.backarc_start_h_m[dst] = self.backarc_start_h_m[s];
            }
        }
        out
    }

    /// Record an eruption of `kind` adding `dh` m at time `t`.
    pub fn record(&mut self, cell: usize, kind: u8, t: f64, dh: f64) {
        if self.first_t_my[cell] < 0.0 {
            self.first_t_my[cell] = t;
            self.origin[cell] = cell as u32;
        }
        if t >= self.last_t_my[cell] {
            self.last_t_my[cell] = t;
            self.kind[cell] = kind;
        }
        self.added_m[cell] += dh;
    }
}

/// One derived volcano point.
#[derive(Debug, Clone, PartialEq)]
pub struct VolcanoPoint {
    pub kind: u8,
    pub cell: u32,
    pub origin: u32,
    pub first_t_my: f64,
    pub last_t_my: f64,
    pub added_m: f64,
    pub active: bool,
    /// Active plume id for plume points (None for marker points).
    pub plume: Option<usize>,
}

/// Spacing (m) per kind given the arc/hotspot base spacing.
pub fn spacing_for(kind: u8, base_m: f64) -> f64 {
    match kind {
        VENT_RIFT => 150_000.0f64.max(base_m),
        VENT_LIP => 500_000.0f64.max(base_m),
        _ => base_m,
    }
}

/// Minimum magmatic crust (m) for a marker cell to be a vent candidate
/// (rift vents are rift-born and need none).
pub const MIN_VENT_ADDED_M: f64 = 500.0;

/// Volcano points of a marker field at `t_my`. `active_plumes` are (plume
/// id, current cell) of the plumes active at `t_my`, in id order.
pub fn volcano_points(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    m: &VentMarkers,
    t_my: f64,
    active_window_my: f64,
    spacing_m: f64,
    active_plumes: &[(usize, usize)],
) -> Result<Vec<VolcanoPoint>, String> {
    let n = grid.cells();
    m.check(n)?;
    let mut cand: Vec<usize> = (0..n)
        .filter(|&c| {
            m.kind[c] != VENT_NONE
                && m.last_t_my[c] >= 0.0
                && (m.kind[c] == VENT_RIFT || m.added_m[c] >= MIN_VENT_ADDED_M)
        })
        .collect();
    cand.sort_by(|&a, &b| {
        m.last_t_my[b]
            .partial_cmp(&m.last_t_my[a])
            .expect("finite times")
            .then(m.added_m[b].partial_cmp(&m.added_m[a]).expect("finite volumes"))
            .then(a.cmp(&b))
    });
    let mut out = Vec::new();
    // Blocked cells per kind (1..=5) and taken origins per kind.
    let mut blocked: Vec<Vec<bool>> = vec![vec![false; n]; 6];
    let mut origins: Vec<std::collections::BTreeSet<u32>> = vec![Default::default(); 6];
    let block = |blocked: &mut Vec<Vec<bool>>, kind: u8, c: usize| -> Result<(), String> {
        let w = crate::orogen::disc_weights(grid, cents, c, spacing_for(kind, spacing_m), &|_| true)?;
        for (k, _) in w {
            blocked[kind as usize][k] = true;
        }
        Ok(())
    };
    // Active plumes first: they are the live hotspots.
    for &(id, c) in active_plumes {
        out.push(VolcanoPoint {
            kind: VENT_HOTSPOT,
            cell: c as u32,
            origin: c as u32,
            first_t_my: m.first_t_my[c],
            last_t_my: t_my,
            added_m: m.added_m[c],
            active: true,
            plume: Some(id),
        });
        block(&mut blocked, VENT_HOTSPOT, c)?;
    }
    for c in cand {
        let kind = m.kind[c];
        if kind as usize >= blocked.len() || blocked[kind as usize][c] {
            continue;
        }
        if !origins[kind as usize].insert(m.origin[c]) {
            continue;
        }
        out.push(VolcanoPoint {
            kind,
            cell: c as u32,
            origin: m.origin[c],
            first_t_my: m.first_t_my[c],
            last_t_my: m.last_t_my[c],
            added_m: m.added_m[c],
            active: t_my - m.last_t_my[c] <= active_window_my,
            plume: None,
        });
        block(&mut blocked, kind, c)?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn carry_follows_source_and_resets_births() {
        let mut m = VentMarkers::new(3);
        m.record(0, VENT_ISLAND_ARC, 4.0, 700.0);
        m.record(0, VENT_HOTSPOT, 6.0, 100.0);
        assert_eq!(m.kind[0], VENT_HOTSPOT);
        assert_eq!(m.first_t_my[0], 4.0);
        assert_eq!(m.added_m[0], 800.0);
        let c = m.carry(&[NO_SOURCE, 0, 0]);
        assert_eq!(c.kind, vec![VENT_NONE, VENT_HOTSPOT, VENT_HOTSPOT]);
        assert_eq!(c.origin[1], 0);
        assert_eq!(c.first_t_my[0], -1.0);
    }
}

/// VS2 passive native source observation. The recorder is explicitly scoped to
/// the stepping thread and is not simulation/checkpoint state. Disabled calls
/// allocate nothing. A failed observer returns an error; it never changes fluxes.
pub mod supply_support;

pub mod supply_capture {
    use pf_contracts::planet::{PlanetKey,PlanetLattice,WorldScope,PLANET_KIND};
    use pf_contracts::volcanic_supply::*;
    use serde_json::{json,Value};
    use std::{cell::RefCell,collections::{BTreeMap,BTreeSet},fs::{File,OpenOptions},io::{BufWriter,Write},path::Path};
    use super::supply_support::{SupplySupport,SupplyCell};
    type Result<T> = std::result::Result<T,String>;
    const MAX_STEP_SAMPLES:usize=100_000;
    const MAX_STEP_CELLS:usize=2_000_000;
    const MAX_SIDECAR_BYTES:u64=3_000_000_000;
    #[derive(Debug,Clone)]
    pub struct CellVolume {
        pub cell:u32,
        pub kernel_weight:f64,
        pub whole_cell_area_m2:f64,
        pub gross_m3:f64,
        pub admitted_m3:f64,
    }
    #[derive(Debug,Clone)]
    pub struct Sample {
        pub branch:String,
        pub lineage:String,
        pub anchor:u32,
        pub start_my:f64,
        pub end_my:f64,
        pub gross_m3:f64,
        pub placement_denominator_m2:f64,
        pub cells:Vec<CellVolume>,
    }
    impl Sample {
        pub fn to_json(&self)->Value {
            json!({"branch":self.branch,"lineage":self.lineage,"anchor":self.anchor,"start_my":self.start_my,
                "end_my":self.end_my,"gross_m3":self.gross_m3,"placement_denominator_m2":self.placement_denominator_m2,
                "cells":self.cells.iter().map(|c|json!({"cell":c.cell,"kernel_weight":c.kernel_weight,
                    "whole_cell_area_m2":c.whole_cell_area_m2,"gross_m3":c.gross_m3,"admitted_m3":c.admitted_m3})).collect::<Vec<_>>()})
        }
        pub fn placed(branch:&str,lineage:String,anchor:u32,start_my:f64,end_my:f64,
            gross_m3:f64,weights:&[(usize,f64)],areas:&[f64])->Result<Self> {
            let norm:f64=weights.iter().map(|(c,w)|w*areas[*c]).sum();
            if !norm.is_finite() || norm<=0. || !gross_m3.is_finite() || gross_m3<0. {
                return Err("VS2 invalid source placement".into());
            }
            let cells=weights.iter().filter(|(_,w)|*w>0.).map(|&(c,w)| {
                // Identical arithmetic to native place/put/lip, including the
                // conversion back from per-cell thickness to volume.
                let v=(gross_m3*w/norm)*areas[c];
                CellVolume{cell:c as u32,kernel_weight:w,whole_cell_area_m2:areas[c],gross_m3:v,admitted_m3:v}
            }).collect();
            Ok(Self{branch:branch.into(),lineage,anchor,start_my,end_my,gross_m3,
                placement_denominator_m2:norm,cells})
        }
    }
    #[derive(Default)]
    struct Buffer { samples:Vec<Sample>,cells:usize, arc_control:Option<[f64;2]>, plume_control:Option<[f64;3]> }
    thread_local! { static ACTIVE:RefCell<Option<Buffer>>=const{RefCell::new(None)}; }
    pub fn enabled()->bool {ACTIVE.with(|v|v.borrow().is_some())}
    pub fn record(mut samples:Vec<Sample>)->Result<()> {
        samples.retain(|s|s.gross_m3!=0.);
        ACTIVE.with(|v| {
            let mut b=v.borrow_mut();let Some(b)=b.as_mut() else{return Ok(())};
            let count:usize=samples.iter().map(|s|s.cells.len()).sum();
            if b.samples.len()+samples.len()>MAX_STEP_SAMPLES || b.cells+count>MAX_STEP_CELLS {
                return Err("VS2 step observation budget exceeded".into());
            }
            b.cells+=count;b.samples.extend(samples);Ok(())
        })
    }
    pub fn arc_control(admitted:f64,capped:f64)->Result<()> {
        ACTIVE.with(|v| {
            let mut b=v.borrow_mut();if let Some(b)=b.as_mut() {
                if b.arc_control.replace([admitted,capped]).is_some(){return Err("VS2 duplicate arc pass".into());}
            }Ok(())
        })
    }
    pub fn plume_control(ordinary_head:f64,tail:f64,lip:f64)->Result<()> {
        ACTIVE.with(|v| {
            let mut b=v.borrow_mut();if let Some(b)=b.as_mut() {
                if b.plume_control.replace([ordinary_head,tail,lip]).is_some(){return Err("VS2 duplicate plume pass".into());}
                if lip>0. && ordinary_head!=0. {return Err("VS2 duplicated LIP/ordinary head".into());}
            }Ok(())
        })
    }
    fn close(a:f64,b:f64)->Result<f64> {
        let r=(a-b).abs();if !r.is_finite() || r>0.001f64.max(a.abs().max(b.abs())*1e-9) {
            Err(format!("VS2 volume closure: {a} vs {b}"))
        }else{Ok(r)}
    }
    fn known<T>(value:T)->Knowledge<T>{Knowledge::Known{value}}
    fn unknown<T>(reason:&str)->Knowledge<T>{Knowledge::Unknown{reason:reason.into()}}
    fn evidence(value:f64,unit:&str,reference:&str)->EvidenceValue {
        EvidenceValue{value,unit:unit.into(),interval:EvidenceInterval::NotReported{reason:"native model constant; no inferred physical confidence interval".into()},
            conditions:"C7 legacy uniform accounting; not composition-specific density or measured eruptive chemistry".into(),
            reference:reference.into(),evidence_kind:"pinned_native_recipe".into()}
    }
    fn mass(volume:f64,rho:f64)->Result<MassUg> {
        let x=volume*rho*1e9;
        if !x.is_finite() || x<0. || x>=i128::MAX as f64 {return Err("VS2 mass conversion range".into());}
        MassUg::new(x.round_ties_even() as i128)
    }
    fn source(s:&Sample,world:&WorldScope,lattice:&PlanetLattice,radius:f64,rho:f64)->Result<MeltSupply> {
        let identity=SourceIdentity{world:world.clone(),anchor:PlanetKey::Cell{cell:s.anchor},lineage:s.lineage.clone(),
            interval:TimeInterval{start_my:s.start_my,end_my:s.end_my}};
        let id=identity.id()?;
        let sum:f64=s.cells.iter().map(|c|c.gross_m3).sum();close(s.gross_m3,sum)?;
        let admitted:f64=s.cells.iter().map(|c|c.admitted_m3).sum();
        let rejected:f64=s.cells.iter().map(|c|c.gross_m3-c.admitted_m3).sum();close(sum,admitted+rejected)?;
        // Quantize the two partition leaves independently. Define the parent as
        // their exact integer sum and retain its difference from the native
        // potential in the adjacent quantity receipt. Never clip a cap failure.
        let am=mass(admitted,rho)?;let rm=mass(rejected,rho)?;
        let denominator=1u64<<52;
        let mut weights:Vec<u64>=s.cells.iter().map(|c|((c.gross_m3/sum)*denominator as f64).floor().max(1.) as u64).collect();
        let used:u64=weights.iter().sum();
        let largest=s.cells.iter().enumerate().max_by(|a,b|a.1.gross_m3.total_cmp(&b.1.gross_m3)).ok_or("empty supply support")?.0;
        if used<=denominator {weights[largest]+=denominator-used;}
        else {weights[largest]=weights[largest].checked_sub(used-denominator).ok_or("VS2 support quantization")?;}
        let support=NativeSupport{radius_m:radius,denominator,
            cells:s.cells.iter().zip(weights).map(|(c,w)|SupportCell{cell:PlanetKey::Cell{cell:c.cell},whole_cell_area_m2:c.whole_cell_area_m2,
                represented_area_m2:c.whole_cell_area_m2,weight:w}).collect(),
            selection_rule:"complete native placement support; integer volume shares, binary64 originals in adjacent sample; whole spherical cells".into(),
            evidence_ref:format!("sources.jsonl#{id}/sample")};
        let density=DensityBasis{id:format!("c7-legacy-crust-{}",rho.to_bits()),kind:DensityKind::LegacyUniformDry,
            density_kg_m3:evidence(rho,"kg/m3","recipe.inputs.profile.planet.crust_density_kg_m3"),
            intrinsic_pore_fraction:unknown("not modeled"),interparticle_pore_fraction:unknown("not modeled")};
        let v=MeltSupply{id:id.clone(),identity,model_version:"c7-passive-source-observation-v1".into(),
            origin:SupplyOrigin::BoundaryInflux{law_version:s.branch.clone(),evidence_ref:"native arc/plume/LIP production law; no finite mantle debit".into()},support,
            dry_rock:ProductionPartition{gross_potential:known(am.checked_add(rm)?),admitted:known(am),capped_or_rejected:known(rm),unextracted:known(MassUg::ZERO)},
            volatile:ProductionPartition{gross_potential:unknown("not modeled"),admitted:unknown("not modeled"),capped_or_rejected:unknown("not modeled"),unextracted:unknown("not modeled")},
            density_basis:known(density),material:unknown("source branch does not determine chemistry or mineral modes"),
            provenance_refs:vec![s.lineage.clone()],quantity_evidence_refs:vec![format!("sources.jsonl#{id}/quantity_receipt")],constants:vec![],composition:None};
        v.validate(world,lattice)?;Ok(v)
    }
    struct Stream { writer:BufWriter<File>,bytes:u64,rows:u64 }
    impl Stream {
        fn new(path:&Path)->Result<Self> {Ok(Self{writer:BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path).map_err(|e|e.to_string())?),bytes:0,rows:0})}
        fn row(&mut self,value:&Value)->Result<()> {
            let bytes=serde_json::to_vec(value).map_err(|e|e.to_string())?;
            if self.bytes+bytes.len() as u64+1>MAX_SIDECAR_BYTES {return Err("VS2 sidecar disk budget exceeded".into());}
            self.writer.write_all(&bytes).and_then(|_|self.writer.write_all(b"\n")).map_err(|e|e.to_string())?;
            self.bytes+=bytes.len() as u64+1;self.rows+=1;Ok(())
        }
        fn flush(&mut self)->Result<()> {self.writer.flush().map_err(|e|e.to_string())}
    }
    /// Scoped process-local observer. Not Send: native stepping and hooks must
    /// execute on the owner thread. Source hooks themselves are serial.
    pub struct Capture {
        sources:Stream,supports:Stream,index:Vec<Value>,world:WorldScope,lattice:PlanetLattice,
        radius:f64,rho:f64,arc_fraction:f64,plume_fraction:f64,
        seen:BTreeMap<String,[u64;2]>,totals:BTreeMap<String,[f64;3]>,targets:BTreeMap<u64,Value>,
        max_residual:f64,checked:usize,root:std::path::PathBuf,
        _thread:std::marker::PhantomData<std::rc::Rc<()>>,
    }
    impl Drop for Capture {fn drop(&mut self){ACTIVE.with(|v|{v.borrow_mut().take();});}}
    impl Capture {
        pub fn start(root:&Path,world:WorldScope,grid:&pf_grid::GridSpec,rho:f64,arc_fraction:f64,plume_fraction:f64)->Result<Self> {
            if enabled(){return Err("VS2 observer already active on thread".into());}
            world.validate()?;
            if !rho.is_finite() || rho<=0. || [arc_fraction,plume_fraction].iter().any(|x|!x.is_finite() || !(0.0..=1.0).contains(x)) {
                return Err("VS2 invalid density/partition".into());
            }
            std::fs::create_dir(root).map_err(|e|e.to_string())?;
            let sources=Stream::new(&root.join("sources.jsonl"))?;
            let supports=Stream::new(&root.join("supports.jsonl"))?;
            ACTIVE.with(|v|*v.borrow_mut()=Some(Buffer::default()));
            Ok(Self{sources,supports,index:Vec::new(),world,lattice:PlanetLattice{kind:PLANET_KIND.into(),rows:grid.nlat() as u32,cols:grid.nlon() as u32,subdivisions:1},
                radius:grid.radius_m(),rho,arc_fraction,plume_fraction,seen:BTreeMap::new(),totals:BTreeMap::new(),targets:BTreeMap::new(),
                max_residual:0.,checked:0,root:root.into(),_thread:std::marker::PhantomData})
        }
        pub fn step(&mut self,grid:&pf_grid::GridSpec,layer:&crate::conserve::ParcelLayer,state:&crate::kinematic::CrustState)->Result<()> {
            let b=ACTIVE.with(|v|v.borrow_mut().replace(Buffer::default())).ok_or("VS2 inactive observer")?;
            let mut wanted=BTreeSet::new();for s in &b.samples {for c in &s.cells {wanted.insert(c.cell);}}
            let h=layer.history.as_ref().ok_or("VS2 requires selected native history")?;
            let support=super::supply_support::observe(h,grid,layer,state,&wanted)?;
            self.write_step(state.t_my,b,support)
        }
        fn write_step(&mut self,time:f64,b:Buffer,support:SupplySupport)->Result<()> {
            let cell_map:BTreeMap<u32,&SupplyCell>=support.cells.iter().map(|c|(c.cell,c)).collect();
            let source_offset=self.sources.bytes;let source_row=self.sources.rows;let support_offset=self.supports.bytes;
            let mut control=[0.;2];
            let mut plume_control=[0.;3];
            let mut shares:BTreeMap<u64,BTreeMap<String,f64>>=BTreeMap::new();
            let targets:BTreeSet<_>=support.targets.iter().map(|t|t.carrier).collect();
            for s in b.samples {
                if s.end_my!=time {return Err("VS2 source step mismatch".into());}
                let contract=source(&s,&self.world,&self.lattice,self.radius,self.rho)?;
                if self.seen.contains_key(&contract.id){return Err("VS2 duplicate source event".into());}
                let source_id=contract.id.clone();
                let offset=self.sources.bytes;
                let fraction=if s.branch=="arc" {self.arc_fraction}else{self.plume_fraction};
                let admitted:f64=s.cells.iter().map(|c|c.admitted_m3).sum();
                let capped:f64=s.cells.iter().map(|c|c.gross_m3-c.admitted_m3).sum();
                if s.branch=="arc" {control[0]+=admitted;control[1]+=capped;}
                match s.branch.as_str() {"plume_head"=>plume_control[0]+=admitted,"plume_tail"=>plume_control[1]+=admitted,"lip_head"=>plume_control[2]+=admitted,_=>{}}
                let totals=self.totals.entry(s.branch.clone()).or_insert([0.;3]);totals[0]+=s.gross_m3;totals[1]+=admitted;totals[2]+=capped;
                let mut routes=Vec::new();
                for c in &s.cells {
                    let native=cell_map.get(&c.cell).ok_or("VS2 missing source cell support")?;
                    let mut allocated=0.;let mut target_routes=Vec::new();
                    if native.live && native.deposited_area_m2>0. {
                        for n in &native.neighbours {
                            let volume=c.admitted_m3*(n.area_m2*n.weight/native.deposited_area_m2);
                            allocated+=volume;
                            if targets.contains(&n.carrier) {
                                *shares.entry(n.carrier).or_default().entry(s.branch.clone()).or_default()+=volume;
                                target_routes.push(json!({"carrier":n.carrier,"volume_m3":volume}));
                            }
                        }
                        self.max_residual=self.max_residual.max(close(c.admitted_m3,allocated)?);
                    }
                    routes.push(json!({"cell":c.cell,"existing_carriers_m3":allocated,
                        "unrepresented_or_birth_import_m3":if native.live && native.deposited_area_m2>0. {0.}else{c.admitted_m3},
                        "target_carriers":target_routes,"support_ref":format!("supports.jsonl@{support_offset}/{}",c.cell)}));
                }
                let native_gross_mass=s.gross_m3*self.rho;
                let parent_mass=contract.dry_rock.gross_potential.require()?.micrograms() as f64*1e-9;
                self.sources.row(&json!({"kind":"source","supply":contract,"sample":s.to_json(),"cell_routes":routes,
                    "partition":{"native_extrusive_fraction":fraction,"extrusive_m3":admitted*fraction,"intrusive_m3":admitted*(1.-fraction),
                        "qualification":"fraction is native only for existing carrier operations; unresolved/birth recipients are not admitted erupted stock",
                        "lava_vs_pyroclastic":"unknown; native model has no such partition","timing_below_step":"conditional, not recovered"},
                    "quantity_receipt":{"basis":"legacy uniform dry mass proxy using pinned crust density; no chemical dry/volatile solve",
                        "native_gross_kg":native_gross_mass,"exact_partition_parent_kg":parent_mass,
                        "parent_rounding_difference_kg":parent_mass-native_gross_mass,
                        "binary64_bound_kg":native_gross_mass.abs()*1e-12+1e-9,
                        "integer_support_denominator":1u64<<52,"integer_share_absolute_error_bound":(s.cells.len()+1) as f64/(1u64<<52) as f64}}))?;
                self.seen.insert(source_id,[offset,self.sources.bytes-offset]);
            }
            if let Some(expected)=b.arc_control {for i in 0..2 {self.max_residual=self.max_residual.max(close(control[i],expected[i])?);}}
            else {return Err("VS2 missing native arc control".into());}
            if let Some(expected)=b.plume_control {for i in 0..3 {self.max_residual=self.max_residual.max(close(plume_control[i],expected[i])?);}}
            else if plume_control.iter().any(|v|*v!=0.) {return Err("VS2 missing native plume control".into());}
            self.checked+=support.checked_native_denominators;
            for t in &support.targets {
                let row=self.targets.entry(t.carrier).or_insert_with(||json!({"carrier":t.carrier,"first_kind":t.kind,"birth_my":time,
                    "events":0,"operations":0,"native_arc_m3":0.,"native_plume_m3":0.,"birth_import_m3":t.birth_import_m3,
                    "source_branches_m3":{},"arc_tag_minus_fresh_source_m3":0.}));
                row["events"]=json!(row["events"].as_u64().unwrap()+1);
                if t.kind==1 {row["operations"]=json!(row["operations"].as_u64().unwrap()+1);}
                for (key,v) in [("native_arc_m3",t.arc_m3),("native_plume_m3",t.plume_m3)] {
                    row[key]=json!(row[key].as_f64().unwrap()+v);
                }
                let branch=shares.get(&t.carrier).cloned().unwrap_or_default();
                for (key,v) in &branch {let prev=row["source_branches_m3"][key].as_f64().unwrap_or(0.);row["source_branches_m3"][key]=json!(prev+v);}
                let residual=t.arc_m3-branch.get("arc").copied().unwrap_or(0.);
                row["arc_tag_minus_fresh_source_m3"]=json!(row["arc_tag_minus_fresh_source_m3"].as_f64().unwrap()+residual);
            }
            self.supports.row(&json!({"time_my":time,"support":support.to_json(),"target_fresh_sources_m3":shares,
                "route_basis":"observed native cell-to-carrier kernel; births/import and noncarrier cells are explicit unresolved recipients; stack net lineage is separately recorded"}))?;
            self.sources.flush()?;self.supports.flush()?;
            self.index.push(json!({"end_my":time,"sources_byte_start":source_offset,"sources_byte_end":self.sources.bytes,
                "source_row_start":source_row,"source_row_count":self.sources.rows-source_row,
                "supports_byte_start":support_offset,"supports_byte_end":self.supports.bytes}));
            Ok(())
        }
        pub fn finish(mut self)->Result<Value> {
            if ACTIVE.with(|v|v.borrow().as_ref().is_some_and(|b|!b.samples.is_empty())) {return Err("VS2 unflushed source events".into());}
            self.sources.flush()?;self.supports.flush()?;
            let mut source_index=OpenOptions::new().write(true).create_new(true).open(self.root.join("source-index.json")).map_err(|e|e.to_string())?;
            serde_json::to_writer(&mut source_index,&self.seen).map_err(|e|e.to_string())?;
            source_index.flush().map_err(|e|e.to_string())?;
            let summary=json!({"schema":"pf.native-volcanic-capture/v1","world":self.world,"lattice":self.lattice,
                "steps":self.index,"source_rows":self.sources.rows,"source_bytes":self.sources.bytes,"support_bytes":self.supports.bytes,
                "branch_gross_admitted_capped_m3":self.totals,"targets":self.targets.values().collect::<Vec<_>>(),
                "checked_native_denominators":self.checked,"missing_required_neighbours":0,"max_volume_closure_residual_m3":self.max_residual,
                "density_kg_m3":self.rho,"arc_extrusive_fraction":self.arc_fraction,"plume_extrusive_fraction":self.plume_fraction,
                "limits":"observation only; fresh source projection is distinct from net stack lineage and retained stock; births imported; substep timing conditional; lava/pyroclastic and volatiles unknown; no cone funding or extraction authority"});
            let mut f=OpenOptions::new().write(true).create_new(true).open(self.root.join("index.json")).map_err(|e|e.to_string())?;
            serde_json::to_writer(&mut f,&summary).map_err(|e|e.to_string())?;f.flush().map_err(|e|e.to_string())?;
            Ok(summary)
        }
    }
}
