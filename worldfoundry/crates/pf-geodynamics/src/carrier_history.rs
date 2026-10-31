//! Optional native selected-carrier stack history. Does not assert material ancestry.
use crate::conserve::CarrierIdentity;
use crate::stratigraphy::{apply_parcel,birth_stack,ParcelAmounts,StrataLayer,StrataParams,StrataStack,StrataStepLog,StrataV2Params,NLITH,NCLASS};
use pf_contracts::params::CarrierHistoryBlock;
use serde_json::{json,Value};
use std::collections::BTreeMap;
type Result<T> = std::result::Result<T,String>;
trait Wire: Sized {
    fn wire(&self) -> Value;
    fn unwire(value: &Value) -> Result<Self>;
}
macro_rules! integer_wire {
    ($t:ty) => {
        impl Wire for $t {
            fn wire(&self) -> Value {
                json!(*self)
            }
            fn unwire(v: &Value) -> Result<Self> {
                v.as_u64()
                    .and_then(|x| <$t>::try_from(x).ok())
                    .ok_or("invalid unsigned integer".into())
            }
        }
    };
}
integer_wire!(u8);
integer_wire!(u32);
integer_wire!(u64);
integer_wire!(usize);
impl Wire for bool {
    fn wire(&self) -> Value {
        json!(*self)
    }
    fn unwire(v: &Value) -> Result<Self> {
        v.as_bool().ok_or("invalid bool".into())
    }
}
impl Wire for f64 {
    fn wire(&self) -> Value {
        json!(format!("{:016x}", self.to_bits()))
    }
    fn unwire(v: &Value) -> Result<Self> {
        let s = v.as_str().ok_or("float must be bit string")?;
        if s.len() != 16
            || !s
                .bytes()
                .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
        {
            return Err("invalid binary64 hex".into());
        }
        let f = Self::from_bits(u64::from_str_radix(s, 16).map_err(|_| "invalid float bits")?);
        if !f.is_finite() {
            return Err("nonfinite float".into());
        }
        Ok(f)
    }
}
impl<T: Wire> Wire for Vec<T> {
    fn wire(&self) -> Value {
        Value::Array(self.iter().map(Wire::wire).collect())
    }
    fn unwire(v: &Value) -> Result<Self> {
        v.as_array()
            .ok_or("expected array")?
            .iter()
            .map(T::unwire)
            .collect()
    }
}
impl<T: Wire, const N: usize> Wire for [T; N] {
    fn wire(&self) -> Value {
        Value::Array(self.iter().map(Wire::wire).collect())
    }
    fn unwire(v: &Value) -> Result<Self> {
        let a = v.as_array().ok_or("expected fixed array")?;
        if a.len() != N {
            return Err("wrong array length".into());
        }
        Vec::<T>::unwire(v)?
            .try_into()
            .map_err(|_| "wrong array length".into())
    }
}
impl<T: Wire> Wire for Option<T> {
    fn wire(&self) -> Value {
        self.as_ref().map_or(Value::Null, Wire::wire)
    }
    fn unwire(v: &Value) -> Result<Self> {
        if v.is_null() {
            Ok(None)
        } else {
            Ok(Some(T::unwire(v)?))
        }
    }
}
macro_rules! struct_wire {
    ($t:ident { $($f:ident),* $(,)? }) => { impl Wire for $t {
        fn wire(&self) -> Value {
            let mut object=serde_json::Map::new();
            $(object.insert(stringify!($f).into(), self.$f.wire());)*
            Value::Object(object)
        }
        fn unwire(v: &Value) -> Result<Self> {
            let object=v.as_object().ok_or("expected object")?;
            let fields=[$(stringify!($f)),*];
            if object.len()!=fields.len() || !fields.iter().all(|f| object.contains_key(*f)) {
                return Err(concat!("wrong fields: ", stringify!($t)).into());
            }
            Ok(Self { $($f: Wire::unwire(&object[stringify!($f)])?),* })
        }
    }};
}

use crate::layer_evolution::{LayerChange,LayerEvolution};
struct_wire!(LayerChange {kind,links,before,after});
struct_wire!(LayerEvolution {changes});

struct_wire!(ParcelAmounts {
    erosion,
    deposit,
    class,
    arc,
    plume,
    thinning,
    basement_add,
    short_pos,
    short_neg,
    land,
    belt_rest,
    warm_shelf
});

struct_wire!(StrataLayer {
    thickness_m,
    age_base_my,
    age_top_my,
    max_burial_m,
    lith,
    env,
    grade,
    flags
});

struct_wire!(StrataStack {
    basement_m,
    basement_lith,
    basement_age_my,
    last_exposure_my,
    cover,
    rock_uplift_m,
    elev_ref_m
});

struct_wire!(StrataParams {
    max_layers,
    merge_gap_my,
    unconformity_gap_my,
    min_layer_m,
    trace_m,
    arc_extrusive_fraction,
    plume_extrusive_fraction,
    shortening_pluton_fraction,
    pluton_roof_m,
    metamorphic_depth_m,
    high_grade_depth_m,
    shelf_depth_m,
    carbonate_max_paleolat_deg,
    carbonate_max_supply_m_per_my,
    arid_lat_min_deg,
    arid_lat_max_deg,
    humid_lat_max_deg,
    humid_midlat_min_deg,
    humid_midlat_max_deg,
    coal_max_elevation_m,
    conglomerate_relief_m,
    river_mouth_area_m2,
    publish_depth_m,
    basement_min_m,
    v2
});

struct_wire!(StrataV2Params {
    shelf_width_m,
    slope_foot_depth_m,
    shelf_majority_fraction,
    platform_cover_max_m,
    platform_full_elev_m,
    platform_zero_elev_m,
    platform_sandstone_fraction,
    platform_carbonate_fraction,
    platform_age_base_my,
    platform_age_top_my
});

struct_wire!(StrataStepLog {
    parcels,
    stacks,
    before_m3,
    after_m3,
    column_m3,
    consumed_m3,
    births_m3,
    eroded_m3,
    deposited_m3,
    volcanic_m3,
    pluton_m3,
    plume_basement_m3,
    basement_add_m3,
    thinned_m3,
    shortening_scale_m3,
    shortening_removed_m3,
    basement_flow_m3,
    basement_flow_abs_m3,
    underflow_m3,
    metamorphosed_m3,
    merges_gap,
    merges_cap,
    merges_min,
    closure_residual_m3,
    closure_rel,
    column_mismatch_rel,
    entry_mismatch_max_rel,
    mean_cover_layers,
    stacks_at_cap,
    basement_below_min,
    belt_basement_m3,
    warm_shelf_deposit_m3,
    merges_unconformity,
    mean_abs_rock_uplift_m
});

fn valid_stack(s: &StrataStack) -> bool {
    [
        s.basement_m,
        s.basement_age_my,
        s.last_exposure_my,
        s.rock_uplift_m,
        s.elev_ref_m,
    ]
    .iter()
    .all(|x| x.is_finite())
        && s.basement_m >= 0.0
        && s.total_m().is_finite()
        && (s.basement_lith as usize) < NLITH
        && s.cover.iter().all(|l| {
            [l.thickness_m, l.age_base_my, l.age_top_my, l.max_burial_m]
                .iter()
                .all(|x| x.is_finite())
                && l.thickness_m >= 0.0
                && l.max_burial_m >= 0.0
                && l.age_top_my >= l.age_base_my
                && (l.lith as usize) < NLITH
        })
}

fn valid_params(p: &StrataParams) -> bool {
    p.validate().is_ok()
        && p.max_layers <= 64
        && [
            p.merge_gap_my,
            p.unconformity_gap_my,
            p.min_layer_m,
            p.trace_m,
            p.arc_extrusive_fraction,
            p.plume_extrusive_fraction,
            p.shortening_pluton_fraction,
            p.pluton_roof_m,
            p.metamorphic_depth_m,
            p.high_grade_depth_m,
            p.shelf_depth_m,
            p.carbonate_max_paleolat_deg,
            p.carbonate_max_supply_m_per_my,
            p.arid_lat_min_deg,
            p.arid_lat_max_deg,
            p.humid_lat_max_deg,
            p.humid_midlat_min_deg,
            p.humid_midlat_max_deg,
            p.coal_max_elevation_m,
            p.conglomerate_relief_m,
            p.river_mouth_area_m2,
            p.publish_depth_m,
            p.basement_min_m,
        ]
        .iter()
        .all(|x| x.is_finite())
        && p.v2.as_ref().is_none_or(|v| {
            [
                v.shelf_width_m,
                v.slope_foot_depth_m,
                v.shelf_majority_fraction,
                v.platform_cover_max_m,
                v.platform_full_elev_m,
                v.platform_zero_elev_m,
                v.platform_sandstone_fraction,
                v.platform_carbonate_fraction,
                v.platform_age_base_my,
                v.platform_age_top_my,
            ]
            .iter()
            .all(|x| x.is_finite())
        })
}

fn valid_log(l: &StrataStepLog) -> bool {
    [
        l.before_m3,
        l.after_m3,
        l.column_m3,
        l.volcanic_m3,
        l.pluton_m3,
        l.plume_basement_m3,
        l.basement_add_m3,
        l.thinned_m3,
        l.shortening_scale_m3,
        l.shortening_removed_m3,
        l.basement_flow_m3,
        l.basement_flow_abs_m3,
        l.underflow_m3,
        l.metamorphosed_m3,
        l.closure_residual_m3,
        l.closure_rel,
        l.column_mismatch_rel,
        l.entry_mismatch_max_rel,
        l.mean_cover_layers,
        l.belt_basement_m3,
        l.warm_shelf_deposit_m3,
        l.mean_abs_rock_uplift_m,
    ]
    .iter()
    .all(|x| x.is_finite())
        && l.consumed_m3
            .iter()
            .chain(&l.births_m3)
            .chain(&l.eroded_m3)
            .chain(&l.deposited_m3)
            .all(|x| x.is_finite())
}

/// Payload byte quota includes a conservative fixed allowance for metadata.
const HEADER_CHARGE: u64 = 65536;
/// Fixed fields and binary64 values encode below this conservative bound.
const SUPPORT_BYTES: u64 = 1024;
/// Bounds the geometry cache used to validate recorded cell IDs and weights.
const MAX_HISTORY_GRID_CELLS: usize = 4_194_304;
#[derive(Debug,Clone,PartialEq)]
pub struct CellSupport {
    pub cell:u32, pub weight:f64, pub live:bool,
    pub area_m2:f64, pub deposited_area_m2:f64,
    /// erosion, deposit, signed arc delta, plume, thinning, basement add,
    /// signed shortening delta, placed belt thickness (v2 only).
    pub raw:[f64;8],
    /// Exact producer g: erosion, deposit, arc, plume, thinning, basement add,
    /// positive shortening, negative shortening, signed belt remainder.
    pub per_unit:[f64;9],
    pub class:u8, pub land:bool, pub warm_shelf:bool,
    pub tect_elev_m:f64, pub post_elev_m:f64,
}
struct_wire!(CellSupport {cell,weight,live,area_m2,deposited_area_m2,raw,per_unit,class,land,warm_shelf,tect_elev_m,post_elev_m});
fn aggregate_support(rows:&[CellSupport],v2:bool)->Result<(ParcelAmounts,Option<[f64;2]>)> {
    let mut am=[0.0;9];let(mut best_dep,mut class,mut warm)=(0.,0u8,false);
    let(mut best_w,mut land)=(-1.,false);let(mut ew,mut et,mut ep)=(0.,0.,0.);
    for row in rows {
        if row.raw.iter().chain(row.per_unit.iter()).copied()
            .chain([row.weight,row.area_m2,row.deposited_area_m2,row.tect_elev_m,row.post_elev_m])
            .any(|x|!x.is_finite()) || row.weight<=0. || row.weight>1. || row.area_m2<=0.
            || row.deposited_area_m2<0. || row.class as usize>=NCLASS {
            return Err("invalid native cell support".into());
        }
        let mut expected=[0.0;9];
        if row.live && row.deposited_area_m2>0. {
            let f=row.area_m2/row.deposited_area_m2;let ds=row.raw[6];
            let(pos,neg,rest)=if v2 {let eff=ds.clamp(0.,row.raw[7].max(0.));(eff,0.,ds-eff)} else {(ds.max(0.),(-ds).max(0.),0.)};
            expected=[row.raw[0]*f,row.raw[1]*f,row.raw[2].max(0.)*f,row.raw[3]*f,
                row.raw[4]*f,row.raw[5]*f,pos*f,neg*f,rest*f];
        }
        if expected.wire()!=row.per_unit.wire(){return Err("native cell tag/weight conversion differs".into())}
        if row.weight>best_w {best_w=row.weight;land=row.land;}
        if !row.live {continue}
        for i in 0..9 {am[i]+=row.per_unit[i]*row.weight;}
        let d=row.per_unit[1]*row.weight;
        if d>best_dep {best_dep=d;class=row.class;warm=v2&&row.warm_shelf;}
        if v2 {ew+=row.weight;et+=row.weight*row.tect_elev_m;ep+=row.weight*row.post_elev_m;}
    }
    Ok((ParcelAmounts{erosion:am[0],deposit:am[1],class,arc:am[2],plume:am[3],thinning:am[4],basement_add:am[5],short_pos:am[6],short_neg:am[7],land,belt_rest:am[8],warm_shelf:warm},
        if v2&&ew>0. {Some([et/ew,ep/ew])}else{None}))
}
#[derive(Debug,Clone,PartialEq)]
struct Event {
    // 0 initial, 1 operation, 2 birth, 3 retirement, 4 global step ledger.
    kind:u8, id:Option<u64>, start_my:f64, end_my:f64,
    position:[f64;3], plate:u32, area_m2:f64,
    before:Option<StrataStack>, after:Option<StrataStack>,
    amounts:Option<ParcelAmounts>, reference:Option<[f64;2]>,
    birth_inputs:Option<[f64;3]>, class:u8, total_after_m:f64,
    log_before:Option<StrataStepLog>, log_after:Option<StrataStepLog>, consumed:Vec<f64>,
    support:Vec<CellSupport>,
    evolution:Option<LayerEvolution>,
}
struct_wire!(Event {kind,id,start_my,end_my,position,plate,area_m2,before,after,amounts,reference,birth_inputs,class,total_after_m,log_before,log_after,consumed,support,evolution});
impl Event {
    fn new(kind:u8,id:Option<u64>,start:f64,end:f64,pos:[f64;3],plate:u32,area:f64)->Self {
        Self{kind,id,start_my:start,end_my:end,position:pos,plate,area_m2:area,
            before:None,after:None,amounts:None,reference:None,birth_inputs:None,
            class:0,total_after_m:0.,log_before:None,log_after:None,consumed:Vec::new(),support:Vec::new(),evolution:None}
    }
}
#[derive(Debug,Clone,PartialEq)]
pub struct StackHistory {
    config:CarrierHistoryBlock, params:StrataParams, events:Vec<Event>, charged_bytes:u64,
    grid:[usize;2], radius_m:f64,
    // A compatible imported v3 archive keeps its original producer attribution.
    // Reading/re-exporting it does not silently stamp it with this build's pins.
    retained_sources:Option<Value>,
}
fn sources()->Value {
    use pf_contracts::sim_module::digest_bytes;
    json!({"stratigraphy":digest_bytes(include_bytes!("stratigraphy.rs")),
        "producer":digest_bytes(include_bytes!("conserve.rs")),
        "codec":digest_bytes(include_bytes!("carrier_history.rs")),"layer_observer":digest_bytes(include_bytes!("layer_evolution.rs"))})
}

// D-0155 / Q27: version 1 of the v3 read-compatibility policy. These are complete
// source bundles, never independently interchangeable producer/codec digests.
// 59109e2 matches retained C7 rec1; 562087f matches the corrected 2026-10-25
// native-body-source capture. Later opt-in observers and reservation bounds do
// not replace replay: from_json still verifies all native operations and charges.
const V3_SOURCE_COMPATIBILITY_V1: [[&str; 4]; 2] = [
    [
        "blake3:1bcb4069251626492af2396b0fc84d099d0b9306663f52b457aac6ea86cf866d",
        "blake3:71c56c3fd77a450cbbedf19cad14d249e5a7b376e03e2759dbd6a60248911919",
        "blake3:604937d8de1531763d2f85813e611e200706dd7236850a3ea9603bee63e0df64",
        "blake3:3d03e7d3326147509359920c3eeeeab743b18696c50e0f0da75f48b574894222",
    ],
    [
        "blake3:eef9b56b51450efc642925cf5bf4386f686496e1e70c6504e7c1e161943c6b94",
        "blake3:21caf2213292eb16d937bc3c3a32722a34ab18b512a8ff484b06cb8ec82bb942",
        "blake3:ec53ebbeca37c4082429e1811c25f8ea00d88821093a4c7def097242048f3e32",
        "blake3:4971e3599296ac1e320a816dbb04beb6aa4b0c2339704b05b0e8ffedffa68319",
    ],
];

fn compatible_v3_sources_v1(value: &Value) -> bool {
    let keys = ["producer", "codec", "stratigraphy", "layer_observer"];
    value.as_object().is_some_and(|o| o.len() == keys.len())
        && V3_SOURCE_COMPATIBILITY_V1.iter().any(|bundle| {
            keys.iter().zip(bundle).all(|(key, digest)| value[*key].as_str() == Some(*digest))
        })
}

impl StackHistory {
    pub fn validate_config(c:&CarrierHistoryBlock)->Result<()> {
        if c.selected.is_empty() || c.selected.len()>256 || c.selected.windows(2).any(|w|w[0]>=w[1])
            || c.event_limit==0 || c.event_limit>1_000_000 || !(HEADER_CHARGE..=1024*1024*1024).contains(&c.byte_limit) {
            return Err("stack history requires 1..256 sorted unique IDs, 1..1000000 events and 64KiB..1GiB bytes".into());
        }
        Ok(())
    }
    pub fn new(config:CarrierHistoryBlock,params:StrataParams,grid:&pf_grid::GridSpec)->Result<Self> {
        Self::validate_config(&config)?;
        if !valid_params(&params) {return Err("invalid stack history parameters".into())}
        if grid.nlat().checked_mul(grid.nlon()).is_none_or(|n|n>MAX_HISTORY_GRID_CELLS) {return Err("native history grid exceeds 4194304-cell verification bound".into())}
        Ok(Self{config,params,events:Vec::new(),charged_bytes:HEADER_CHARGE,
            grid:[grid.nlat(),grid.nlon()],radius_m:grid.radius_m(),retained_sources:None})
    }
    pub fn config(&self)->&CarrierHistoryBlock {&self.config}
    pub fn params(&self)->&StrataParams {&self.params}
    pub fn grid(&self)->[usize;2] {self.grid}
    pub fn radius_m(&self)->f64 {self.radius_m}
    fn support_limit(&self)->usize {self.grid[0].min(3)*self.grid[1]}
    pub fn selected(&self,id:u64)->bool {self.config.selected.binary_search(&id).is_ok()}
    pub fn event_count(&self)->usize {self.events.len()}
    pub fn latest_time(&self)->f64 {self.events.last().map_or(0.,|e|e.end_my)}
    fn require_current_producer(&self)->Result<()> {
        if self.retained_sources.is_some() {
            return Err("compatible historical v3 stack history is read-only; mixed-producer append requires a new versioned capture".into());
        }
        Ok(())
    }
    pub fn validate_schedule(&self,steps:u64,dt:f64)->Result<()> {
        let mut found=0u64;
        for e in self.events.iter().filter(|e|e.kind==4) {
            found+=1;
            if e.start_my!=e.end_my || e.end_my!=found as f64*dt {return Err("stack global ledger schedule mismatch".into())}
        }
        if found!=steps || self.latest_time()!=steps as f64*dt {return Err("missing stack step ledger".into())}
        Ok(())
    }
    fn per_event_bound(&self,support:bool)->u64 {
        self.shaped_event_bound(self.params.max_layers, support)
    }
    // apply_parcel can add at most four cover slots (roof, pluton, volcanic,
    // sediment). Erosion emits <= n removals, one truncation, one basement
    // change and one thin-remainder absorption. Capacity merging emits <= 4
    // changes for an admitted input n < max_layers. All other branches emit
    // <= 17 changes, including both zero filters. Thus n+32 is conservative,
    // even allowing n == max_layers in the context-free bound. Each change
    // has two stacks of <= n+4 slots and <= n+6 lineage links. Keep the old
    // conservative per-change wire allowance; tighten only the operation
    // count and live input shape. Birth starts empty and adds <= three slots.
    fn shaped_event_bound(&self,n:usize,support:bool)->u64 {
        16_384+1024*(n as u64+4)+if support {
            self.support_limit() as u64*SUPPORT_BYTES+if self.config.layer_evolution {
                (n as u64+32)*(2048+1024*(n as u64+4))
            }else{0}
        }else{0}
    }
    /// Whole-step reservation from the native live stacks, before mutation.
    /// A selected live ID either retires or receives one operation. An as-yet
    /// unallocated ID can be born once; a retired ID cannot be reused. The
    /// global ledger receives its own fixed record. Polar support remains
    /// bounded by three complete longitude rows, regardless of motion.
    pub fn live_step_budget(&self,ids:&CarrierIdentity,stacks:&[StrataStack])->Result<(usize,u64)> {
        if ids.ids().len()!=stacks.len() {return Err("stack history live shape mismatch".into())}
        let mut shapes=BTreeMap::new();
        for (&id,stack) in ids.ids().iter().zip(stacks) {
            if self.selected(id) {
                if stack.cover.len()>=self.params.max_layers || !valid_stack(stack)
                    || shapes.insert(id,stack.cover.len()).is_some() {
                    return Err("invalid selected live stack before mutation".into());
                }
            }
        }
        let mut count=1usize;
        let mut bytes=self.shaped_event_bound(0,false);
        for id in &self.config.selected {
            let n=if let Some(n)=shapes.get(id) {*n}
                else if *id>=ids.births().len() as u64 {0}
                else {continue};
            count+=1;
            bytes=bytes.checked_add(self.shaped_event_bound(n,true)).ok_or("stack reservation overflow")?;
        }
        Ok((count,bytes))
    }
    pub fn reserve_live_step(&mut self,ids:&CarrierIdentity,stacks:&[StrataStack])->Result<()> {
        self.require_current_producer()?;
        let (count,bytes)=self.live_step_budget(ids,stacks)?;
        if self.events.len().checked_add(count).is_none_or(|n|n as u64>self.config.event_limit)
            || self.charged_bytes.checked_add(bytes).is_none_or(|n|n>self.config.byte_limit) {
            return Err("stack history budget exhausted before mutation".into());
        }
        self.events.try_reserve(count).map_err(|_|"stack history event allocation".into())
    }
    fn preflight(&self,count:usize,support:bool)->Result<()> {
        self.require_current_producer()?;
        if self.events.len().checked_add(count).is_none_or(|n|n as u64>self.config.event_limit)
            || self.charged_bytes.checked_add(self.per_event_bound(support)*count as u64).is_none_or(|n|n>self.config.byte_limit) {
            return Err("stack history budget exhausted before mutation".into());
        }
        Ok(())
    }
    pub fn budget_estimate(&self)->Value {
        json!({"initial_records_bound":self.config.selected.len(),
            "initial_encoded_bytes_bound":HEADER_CHARGE+self.per_event_bound(false)*self.config.selected.len() as u64,
            "step_records_bound":self.config.selected.len()+1,
            "step_encoded_reservation_bytes":self.per_event_bound(true)*(self.config.selected.len()+1) as u64,
            "support_cells_per_operation_bound":self.support_limit(),
            "configured_byte_limit":self.config.byte_limit,"configured_event_limit":self.config.event_limit,
            "layer_evolution":self.config.layer_evolution,"layer_changes_per_operation_bound":self.params.max_layers+32,
            "scope":"context-free maximum; native step uses selected live shapes, four insertion slots and n+32 changes; polar support retained; not an RSS estimate"})
    }
    pub fn preflight_step(&self)->Result<()> {self.preflight(self.config.selected.len()+1,true)}
    pub fn reserve_step(&mut self)->Result<()> {
        self.preflight_step()?;
        self.events.try_reserve(self.config.selected.len()+1).map_err(|_|"stack history event allocation".into())
    }
    pub fn reserve_initial(&mut self)->Result<()> {
        self.preflight(self.config.selected.len(),false)?;
        self.events.try_reserve(self.config.selected.len()).map_err(|_|"stack history initial allocation".into())
    }
    fn push(&mut self,e:Event)->Result<()> {
        self.require_current_producer()?;
        if e.evolution.is_some()!=(self.config.layer_evolution&&(e.kind==1||e.kind==2)){return Err("layer evolution profile mismatch".into())}
        if e.evolution.as_ref().is_some_and(|t|t.changes.len()>crate::layer_evolution::change_bound(self.params.max_layers)){return Err("layer evolution change bound".into())}
        let bytes=serde_json::to_vec(&e.wire()).map_err(|e|e.to_string())?.len() as u64+1;
        let n=if e.kind==2 {0}else{e.before.as_ref().or(e.after.as_ref()).map_or(0,|s|s.cover.len())};
        if e.evolution.as_ref().is_some_and(|t|t.changes.len()>n+32
            || t.changes.iter().any(|c|c.before.cover.len()>n+4 || c.after.cover.len()>n+4))
            || bytes>self.shaped_event_bound(n,e.kind==1||e.kind==2) {
            return Err("stack history live reservation invariant failed".into());
        }
        if e.support.len()>self.support_limit() || e.support.iter().any(|c|serde_json::to_vec(&c.wire()).map_or(true,|v|v.len() as u64+1>SUPPORT_BYTES))
            || bytes>self.per_event_bound(e.kind==1||e.kind==2) || self.events.len() as u64>=self.config.event_limit
            || self.charged_bytes+bytes>self.config.byte_limit {
            return Err("stack history reservation invariant failed".into());
        }
        self.charged_bytes+=bytes;self.events.push(e);Ok(())
    }
    pub fn initial(&mut self,id:u64,time:f64,pos:[f64;3],plate:u32,area:f64,stack:&StrataStack)->Result<()> {
        if !self.selected(id){return Ok(())}
        let mut e=Event::new(0,Some(id),time,time,pos,plate,area);e.after=Some(stack.clone());self.push(e)
    }
    pub fn retirement(&mut self,id:u64,time:f64,pos:[f64;3],plate:u32,area:f64,stack:&StrataStack)->Result<()> {
        if !self.selected(id){return Ok(())}
        let mut e=Event::new(3,Some(id),time,time,pos,plate,area);e.before=Some(stack.clone());self.push(e)
    }
    #[allow(clippy::too_many_arguments)]
    pub fn operation(&mut self,id:u64,start:f64,end:f64,pos:[f64;3],plate:u32,area:f64,
        before:StrataStack,after:&StrataStack,amounts:ParcelAmounts,reference:Option<[f64;2]>,total:f64,
        log_before:StrataStepLog,log_after:&StrataStepLog,support:Vec<CellSupport>,evolution:Option<LayerEvolution>)->Result<()> {
        let mut e=Event::new(1,Some(id),start,end,pos,plate,area);
        e.before=Some(before);e.after=Some(after.clone());e.amounts=Some(amounts);e.reference=reference;
        e.total_after_m=total;e.log_before=Some(log_before);e.log_after=Some(log_after.clone());e.support=support;e.evolution=evolution;self.push(e)
    }
    #[allow(clippy::too_many_arguments)]
    pub fn birth(&mut self,id:u64,time:f64,pos:[f64;3],plate:u32,area:f64,inputs:[f64;3],class:u8,
        reference:Option<[f64;2]>,after:&StrataStack,log_before:StrataStepLog,log_after:&StrataStepLog,evolution:Option<LayerEvolution>)->Result<()> {
        let mut e=Event::new(2,Some(id),time,time,pos,plate,area);e.birth_inputs=Some(inputs);e.class=class;
        e.evolution=evolution;e.reference=reference;e.after=Some(after.clone());e.log_before=Some(log_before);e.log_after=Some(log_after.clone());self.push(e)
    }
    pub fn ledger(&mut self,time:f64,log:&StrataStepLog,consumed:&[f64])->Result<()> {
        let mut e=Event::new(4,None,time,time,[0.,0.,0.],0,0.);e.log_after=Some(log.clone());e.consumed=consumed.to_vec();self.push(e)
    }
    pub fn to_json(&self)->Value {
        json!({"schema":"vz.native-stack-history/v3","scope":"selected_carrier_operations_not_material_ancestry",
            "float_encoding":"ieee754-binary64-lowercase-hex","sources":self.retained_sources.clone().unwrap_or_else(sources),"config":self.config,
            "params":self.params.wire(),"grid":self.grid,"radius_m":self.radius_m.wire(),"charged_bytes":self.charged_bytes,"events":self.events.iter().map(Wire::wire).collect::<Vec<_>>()})
    }
    pub fn from_json(value:Value)->Result<Self> {
        let current_sources = value["sources"] == sources();
        if value.as_object().is_none_or(|o|o.len()!=10)
            || value["schema"]!="vz.native-stack-history/v3"
            || value["scope"]!="selected_carrier_operations_not_material_ancestry"
            || value["float_encoding"]!="ieee754-binary64-lowercase-hex"
            || (!current_sources && !compatible_v3_sources_v1(&value["sources"])) {
            return Err("stack history schema/source mismatch".into());
        }
        let config:CarrierHistoryBlock=serde_json::from_value(value["config"].clone()).map_err(|e|e.to_string())?;
        let params=StrataParams::unwire(&value["params"])?;
        let shape=<[usize;2]>::unwire(&value["grid"])?;
        if shape[0].checked_mul(shape[1]).is_none_or(|n|n>MAX_HISTORY_GRID_CELLS) {return Err("native history grid verification bound".into())}
        let grid=pf_grid::GridSpec::new(shape[0],shape[1],f64::unwire(&value["radius_m"])?).map_err(|e|e.to_string())?;
        let mut h=Self::new(config,params,&grid)?;
        let events=value["events"].as_array().ok_or("missing stack events")?;
        if events.len() as u64>h.config.event_limit {return Err("stack event count limit".into())}
        for v in events {
            let e=Event::unwire(v)?;
            h.push(e)?;
        }
        if value["charged_bytes"].as_u64()!=Some(h.charged_bytes) {return Err("stack byte accounting mismatch".into())}
        h.replay()?;
        if !current_sources {h.retained_sources=Some(value["sources"].clone());}
        Ok(h)
    }
    fn replay(&self)->Result<BTreeMap<u64,(StrataStack,f64,bool)>> {
        let mut known:BTreeMap<u64,(StrataStack,f64,bool)>=BTreeMap::new();
        let mut time=f64::NEG_INFINITY;
        let grid=pf_grid::GridSpec::new(self.grid[0],self.grid[1],self.radius_m).map_err(|e|e.to_string())?;
        let centers=crate::kinematic::centroids(&grid).map_err(|e|e.to_string())?;
        let areas=crate::kinematic::cell_areas(&grid).map_err(|e|e.to_string())?;
        let kernel=crate::conserve::Kernel::new(&grid,&centers,&areas).map_err(|e|e.to_string())?;
        for e in &self.events {
            if !e.start_my.is_finite() || !e.end_my.is_finite() || e.start_my>e.end_my || e.end_my<time {
                return Err("stack history clock order".into());
            }
            time=e.end_my;
            if e.kind!=1 && !e.support.is_empty(){return Err("unexpected cell support on non-operation".into())}
            if e.before.as_ref().is_some_and(|s|!valid_stack(s)||s.cover.len()>=self.params.max_layers)
                || e.after.as_ref().is_some_and(|s|!valid_stack(s)||s.cover.len()>=self.params.max_layers)
                || e.log_before.as_ref().is_some_and(|l|!valid_log(l)) || e.log_after.as_ref().is_some_and(|l|!valid_log(l)) {
                return Err("invalid captured native stack/log".into());
            }
            if e.kind==4 {
                if e.id.is_some() || e.log_after.is_none() || (!e.consumed.is_empty() && e.consumed.len()!=NLITH) || e.consumed.iter().any(|x|!x.is_finite()||*x<0.) {return Err("invalid global step ledger".into())}
                continue;
            }
            let id=e.id.ok_or("missing carrier event ID")?;
            if !self.selected(id) || !e.area_m2.is_finite() || e.area_m2<=0.
                || e.position.iter().any(|x|!x.is_finite()) || (e.position.iter().map(|x|x*x).sum::<f64>()-1.).abs()>1e-8 {
                return Err("invalid carrier event support".into());
            }
            match e.kind {
                0|2 => {
                    if known.contains_key(&id) || e.before.is_some() {return Err("duplicate carrier history birth".into())}
                    let after=e.after.as_ref().ok_or("missing initial/birth stack")?;
                    if e.kind==2 {
                        let x=e.birth_inputs.ok_or("missing birth inputs")?;
                        if x.iter().any(|x|!x.is_finite()||*x<0.) || e.class as usize>=NCLASS {return Err("invalid birth inputs".into())}
                        let mut log=e.log_before.clone().ok_or("missing birth log")?;
                        let mut actual=if let Some(expected)=&e.evolution {
                            let (stack,observed)=crate::stratigraphy::birth_stack_recorded(x[0],x[1],x[2],e.class,&self.params,e.end_my,e.area_m2,&mut log,crate::layer_evolution::byte_bound(self.params.max_layers))?;
                            if observed.wire()!=expected.wire(){return Err("native birth layer evolution replay differs".into())}stack
                        }else{birth_stack(x[0],x[1],x[2],e.class,&self.params,e.end_my,e.area_m2,&mut log)};
                        if let Some(r)=e.reference {actual.elev_ref_m=r[1]}
                        if actual.wire()!=after.wire() || Some(log.wire())!=e.log_after.as_ref().map(Wire::wire) {return Err("native birth replay differs".into())}
                    }
                    known.insert(id,(after.clone(),e.end_my,false));
                }
                1|3 => {
                    let previous=known.get_mut(&id).ok_or("carrier operation before initialization")?;
                    let before=e.before.as_ref().ok_or("missing input stack")?;
                    if previous.2 || before.wire()!=previous.0.wire() {return Err("carrier stack continuity differs".into())}
                    if e.kind==3 {previous.2=true;previous.1=e.end_my;continue}
                    if previous.1!=e.start_my {return Err("carrier operation interval gap".into())}
                    let mut actual=before.clone();
                    if let Some(r)=e.reference {actual.rock_uplift_m+=r[0]-actual.elev_ref_m;actual.elev_ref_m=r[1];}
                    let amounts=e.amounts.as_ref().ok_or("missing operation inputs")?;
                    let mut weights=Vec::new();kernel.weights(e.position,&mut weights).map_err(|e|e.to_string())?;
                    if weights.len()!=e.support.len() || e.support.windows(2).any(|w|w[0].cell>=w[1].cell)
                        || weights.iter().zip(&e.support).any(|((c,w),s)|*c!=s.cell || w.to_bits()!=s.weight.to_bits()
                            || areas[*c as usize].to_bits()!=s.area_m2.to_bits()) {
                        return Err("native kernel cell/weight/area support differs".into());
                    }
                    let (replayed,reference)=aggregate_support(&e.support,self.params.v2.is_some())?;
                    if replayed.wire()!=amounts.wire() || reference.wire()!=e.reference.wire() {
                        return Err("native cell support aggregation differs".into());
                    }
                    if [amounts.erosion,amounts.deposit,amounts.arc,amounts.plume,amounts.thinning,amounts.basement_add,amounts.short_pos,amounts.short_neg].iter().any(|x|!x.is_finite()||*x<0.)
                        || !amounts.belt_rest.is_finite() || amounts.class as usize>=NCLASS || !e.total_after_m.is_finite() || e.total_after_m<0. {return Err("invalid operation inputs".into())}
                    let mut log=e.log_before.clone().ok_or("missing input log")?;
                    if let Some(expected)=&e.evolution {
                        let observed=crate::stratigraphy::apply_parcel_recorded(&mut actual,amounts,e.total_after_m,e.area_m2,&self.params,e.start_my,e.end_my,&mut log,crate::layer_evolution::byte_bound(self.params.max_layers))?;
                        if observed.wire()!=expected.wire(){return Err("native layer evolution replay differs".into())}
                    }else{apply_parcel(&mut actual,amounts,e.total_after_m,e.area_m2,&self.params,e.start_my,e.end_my,&mut log);}
                    if Some(actual.wire())!=e.after.as_ref().map(Wire::wire) || Some(log.wire())!=e.log_after.as_ref().map(Wire::wire) {return Err("native stack operation replay differs".into())}
                    *previous=(actual,e.end_my,false);
                }
                _=>return Err("unknown native stack event".into()),
            }
        }
        Ok(known)
    }
    pub fn validate_context(&self,ids:&CarrierIdentity,time:f64,stacks:Option<&[StrataStack]>)->Result<()> {
        ids.check_context(self.grid[0]*self.grid[1],ids.ids().len())?;
        let known=self.replay()?;
        if self.events.iter().any(|e|e.end_my>time) {return Err("stack event after snapshot".into())}
        for id in &self.config.selected {
            let birth=ids.births().get(*id as usize);
            match (birth,known.get(id)) {
                (None,None)=>{},
                (Some(b),Some((last,t,retired)))=>{
                    let first=self.events.iter().find(|e|e.id==Some(*id)).ok_or("missing carrier initialization")?;
                    if first.end_my!=b.model_time_my || first.kind!=if b.initial {0}else{2} {return Err("stack/carrier allocation mismatch".into())}
                    let death=ids.retirements().iter().find(|r|r.id==*id);
                    if *retired!=death.is_some() || (*retired && *t!=death.unwrap().model_time_my) || (!*retired && *t!=time) {return Err("stack/carrier lifecycle mismatch".into())}
                    if let Some(stacks)=stacks {
                        if !*retired {
                            let k=ids.ids().iter().position(|x|x==id).ok_or("active stack handle absent")?;
                            if stacks.get(k).is_none_or(|s|s.wire()!=last.wire()) {return Err("stack history differs from native final stack".into())}
                        }
                    }
                }
                _=>return Err("selected allocated carrier lacks history".into()),
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod layer_observer_native_archive_tests {
 use super::*;
 #[test]
 #[ignore="requires retained source-matched native C5 verifier and archive"]
 fn earlier_actual_world_operations_preserve_exact_stacks_and_logs(){
  use crate::stratigraphy::apply_parcel_recorded;
  let root=std::path::PathBuf::from(std::env::var("GEO_C6_SOURCE_RUN").expect("source run"));
  let verifier=std::env::var("GEO_C6_SOURCE_VERIFIER").expect("source-matched verifier");
  let output=std::process::Command::new(verifier).arg("--verify-run").arg(&root).output().unwrap();assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
  let data:Value=serde_json::from_slice(&std::fs::read(root.join("carrier-identities.json")).unwrap()).unwrap();
  let h=&data["stack_history"];let p=StrataParams::unwire(&h["params"]).unwrap();let mut count=0;let mut rows=0;let mut kinds=std::collections::BTreeSet::new();let mut max_encoded=0;
  for value in h["events"].as_array().unwrap(){
   let mut value=value.clone();value.as_object_mut().unwrap().entry("evolution").or_insert(Value::Null);
   let e=Event::unwire(&value).unwrap();if e.kind!=1{continue}
   let mut observed=e.before.clone().unwrap();if let Some(r)=e.reference{observed.rock_uplift_m+=r[0]-observed.elev_ref_m;observed.elev_ref_m=r[1];}
   let mut plain=observed.clone();let mut log=e.log_before.clone().unwrap();let mut plain_log=log.clone();let a=e.amounts.as_ref().unwrap();
   apply_parcel(&mut plain,a,e.total_after_m,e.area_m2,&p,e.start_my,e.end_my,&mut plain_log);
   let trace=apply_parcel_recorded(&mut observed,a,e.total_after_m,e.area_m2,&p,e.start_my,e.end_my,&mut log,crate::layer_evolution::byte_bound(p.max_layers)).unwrap();
   assert_eq!(plain.wire(),e.after.as_ref().unwrap().wire());assert_eq!(observed.wire(),plain.wire());assert_eq!(log.wire(),plain_log.wire());assert_eq!(log.wire(),e.log_after.as_ref().unwrap().wire());
   for pair in trace.changes.windows(2){assert_eq!(pair[0].after.wire(),pair[1].before.wire());}
   let encoded=serde_json::to_vec(&trace.wire()).unwrap();max_encoded=max_encoded.max(encoded.len());assert!(encoded.len() as u64<=crate::layer_evolution::byte_bound(p.max_layers));
   for c in &trace.changes{kinds.insert(c.kind);}rows+=trace.changes.len();count+=1;
  }
  assert!(count>500&&rows>count);eprintln!("actual native calls={count} changes={rows} max_encoded_operation={max_encoded} kinds={kinds:?}");
 }
}

#[cfg(test)]
mod formation_record_tests {
    use super::*;
    #[test]
    fn actual_c7_operations_keep_native_bytes_and_transport_sources() {
        let Some(path)=std::env::var_os("VZ_LAYER_PROVENANCE_HISTORY") else {
            eprintln!("SKIP: VZ_LAYER_PROVENANCE_HISTORY unset");return;
        };
        let data:Value=serde_json::from_reader(std::io::BufReader::new(std::fs::File::open(path).unwrap())).unwrap();
        let h=&data["stack_history"];let p=StrataParams::unwire(&h["params"]).unwrap();
        let mut count=0;let mut merges=0;let mut transformations=0;
        for value in h["events"].as_array().unwrap(){
            let mut value=value.clone();value.as_object_mut().unwrap().entry("evolution").or_insert(Value::Null);
            let e=Event::unwire(&value).unwrap();if e.kind!=1{continue}
            let mut native=e.before.clone().unwrap();
            if let Some(r)=e.reference{native.rock_uplift_m+=r[0]-native.elev_ref_m;native.elev_ref_m=r[1];}
            let provenance=crate::layer_provenance::initial(&native,e.id.unwrap());
            let mut plain=native.clone();let mut native_log=e.log_before.clone().unwrap();let mut plain_log=native_log.clone();
            let a=e.amounts.as_ref().unwrap();
            apply_parcel(&mut plain,a,e.total_after_m,e.area_m2,&p,e.start_my,e.end_my,&mut plain_log);
            let (result,trace)=crate::stratigraphy::apply_parcel_provenance(&mut native,a,e.total_after_m,e.area_m2,&p,e.start_my,e.end_my,&mut native_log,provenance,e.id.unwrap(),true).unwrap();
            assert_eq!(plain.wire(),native.wire());assert_eq!(plain_log.wire(),native_log.wire());
            assert_eq!(native.wire(),e.after.as_ref().unwrap().wire());assert_eq!(native_log.wire(),e.log_after.as_ref().unwrap().wire());
            crate::layer_provenance::validate(&result,&native).unwrap();
            merges+=trace.changes.iter().filter(|c|c.kind==17||c.kind==3).count();
            transformations+=trace.changes.iter().filter(|c|c.kind==18).count();count+=1;
        }
        assert!(count>500&&merges>0&&transformations>0);
        eprintln!("actual C7 operations={count}, native merges={merges}, burial/metamorphism observations={transformations}; exact stack/log equality");
    }
}
