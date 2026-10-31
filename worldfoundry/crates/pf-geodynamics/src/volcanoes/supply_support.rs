//! Passive VS2 view over the unchanged authenticated native history codec.
//! No checkpoint source digest, schema, wire encoder or decoder is modified.
use crate::carrier_history::StackHistory;
use serde_json::{json,Value};
use std::collections::BTreeMap;
type Result<T> = std::result::Result<T,String>;
fn number(v:&Value)->Result<f64> {
    let s=v.as_str().ok_or("VS2 expected native binary64")?;
    if s.len()!=16 || !s.bytes().all(|c|c.is_ascii_digit()||(b'a'..=b'f').contains(&c)) {
        return Err("VS2 invalid native binary64".into());
    }
    let f=f64::from_bits(u64::from_str_radix(s,16).map_err(|_|"VS2 native binary64")?);
    if !f.is_finite(){return Err("VS2 nonfinite history view".into());} Ok(f)
}
fn unsigned(v:&Value)->Result<u64>{v.as_u64().ok_or("VS2 native unsigned field".into())}
fn triple(v:&Value)->Result<[f64;3]> {
    let a=v.as_array().ok_or("VS2 native triple")?;
    if a.len()!=3{return Err("VS2 native triple size".into());}
    Ok([number(&a[0])?,number(&a[1])?,number(&a[2])?])
}
struct SupportView {cell:u32,weight:f64,live:bool,deposited_area_m2:f64}
struct EventView {
    id:u64,kind:u8,area_m2:f64,position:[f64;3],support:Vec<SupportView>,
    amounts:Option<[f64;2]>,birth_inputs:Option<[f64;3]>,
}
impl EventView {
    fn read(e:&Value)->Result<Self> {
        let support=e["support"].as_array().ok_or("VS2 native support")?.iter().map(|s|Ok(SupportView {
            cell:u32::try_from(unsigned(&s["cell"])?).map_err(|_|"VS2 cell range")?,
            weight:number(&s["weight"])?,live:s["live"].as_bool().ok_or("VS2 native live flag")?,
            deposited_area_m2:number(&s["deposited_area_m2"])?,
        })).collect::<Result<Vec<_>>>()?;
        Ok(Self{id:unsigned(&e["id"])?,kind:u8::try_from(unsigned(&e["kind"])?).map_err(|_|"VS2 event kind range")?,
            area_m2:number(&e["area_m2"])?,position:triple(&e["position"])?,support,
            amounts:if e["amounts"].is_null(){None}else{Some([number(&e["amounts"]["arc"])?,number(&e["amounts"]["plume"]) ?])},
            birth_inputs:if e["birth_inputs"].is_null(){None}else{Some(triple(&e["birth_inputs"])?)}})
    }
}

/// VS2 observation of the actual write-back support. Recompute the same kernel
/// from the completed native parcels, excluding births which did not participate
/// in this step's write-back. All neighbours are retained, not only selected IDs.
#[derive(Debug, Clone)]
pub struct SupplyNeighbour {
    pub carrier: u64,
    pub weight: f64,
    pub area_m2: f64,
}
#[derive(Debug, Clone)]
pub struct SupplyCell {
    pub cell: u32,
    pub live: bool,
    pub area_m2: f64,
    pub deposited_area_m2: f64,
    pub neighbours: Vec<SupplyNeighbour>,
}
#[derive(Debug, Clone)]
pub struct SupplyTarget {
    pub carrier: u64,
    pub kind: u8,
    pub time_my: f64,
    pub area_m2: f64,
    pub position: [f64;3],
    pub arc_m3: f64,
    pub plume_m3: f64,
    pub birth_import_m3: Option<[f64;3]>,
    pub support_cells: Vec<u32>,
}
#[derive(Debug, Clone)]
pub struct SupplySupport {
    pub cells: Vec<SupplyCell>,
    pub targets: Vec<SupplyTarget>,
    pub checked_native_denominators: usize,
    pub missing_required_neighbours: usize,
}
pub fn observe(history:&StackHistory, grid:&pf_grid::GridSpec, layer:&crate::conserve::ParcelLayer,
        state:&crate::kinematic::CrustState, source_cells:&std::collections::BTreeSet<u32>) -> Result<SupplySupport> {
        use pf_contracts::terrain::CrustKindCodeV0;
        let time=state.t_my;
        let ids=layer.identity.as_ref().ok_or("supply support requires identities")?;
        let centers=crate::kinematic::centroids(grid).map_err(|e|e.to_string())?;
        let areas=crate::kinematic::cell_areas(grid).map_err(|e|e.to_string())?;
        let kernel=crate::conserve::Kernel::new(grid,&centers,&areas).map_err(|e|e.to_string())?;
        // The existing public archive view is sufficient. Do not add an
        // accessor to the authenticated legacy codec source file.
        let archive=history.to_json();
        if archive["schema"]!="vz.native-stack-history/v3" {return Err("VS2 unsupported history view".into());}
        let mut events=Vec::new();
        for e in archive["events"].as_array().ok_or("VS2 missing events")? {
            if number(&e["end_my"])?==time && !e["id"].is_null() {events.push(EventView::read(e)?);}
        }
        drop(archive);
        let mut wanted=source_cells.clone();
        for e in &events {for s in &e.support {wanted.insert(s.cell);}}
        let birth_cells:std::collections::BTreeSet<_>=ids.births().iter()
            .filter(|b|b.model_time_my==time && !b.initial).map(|b|b.event_cell as u32).collect();
        let mut cells:BTreeMap<u32,SupplyCell>=wanted.into_iter().map(|c| (c,SupplyCell {
            cell:c,live:state.kind[c as usize]==CrustKindCodeV0::Continental && !birth_cells.contains(&c),
            area_m2:areas[c as usize],deposited_area_m2:0.,neighbours:Vec::new(),
        })).collect();
        let mut weights=Vec::new();
        for k in 0..layer.len() {
            // Native write-back has already appended these births. They did not
            // contribute to its denominator and own imported initial columns.
            if ids.births()[ids.ids()[k] as usize].model_time_my==time {continue;}
            weights.clear();kernel.weights(layer.pos[k],&mut weights).map_err(|e|e.to_string())?;
            for &(c,w) in &weights {
                if let Some(s)=cells.get_mut(&c) {
                    if s.live {s.deposited_area_m2+=layer.area_m2[k]*w;}
                    s.neighbours.push(SupplyNeighbour{carrier:ids.ids()[k],weight:w,area_m2:layer.area_m2[k]});
                }
            }
        }
        let mut checked=0;
        let mut targets=Vec::new();
        for e in events {
            for s in &e.support {
                let actual=&cells[&s.cell];
                if actual.live!=s.live || actual.deposited_area_m2.to_bits()!=s.deposited_area_m2.to_bits() {
                    return Err(format!("VS2 support denominator differs at carrier {:?} cell {} time {time}: {} vs {}",e.id,s.cell,actual.deposited_area_m2,s.deposited_area_m2));
                }
                let id=e.id;
                let own=actual.neighbours.iter().find(|v|v.carrier==id).ok_or("VS2 missing required neighbour")?;
                if own.weight.to_bits()!=s.weight.to_bits() || own.area_m2.to_bits()!=e.area_m2.to_bits() {
                    return Err("VS2 support weight differs".into());
                }
                checked+=1;
            }
            targets.push(SupplyTarget{carrier:e.id,kind:e.kind,time_my:time,area_m2:e.area_m2,position:e.position,
                arc_m3:e.amounts.map_or(0.,|a|a[0]*e.area_m2),
                plume_m3:e.amounts.map_or(0.,|a|a[1]*e.area_m2),
                birth_import_m3:e.birth_inputs.map(|a|a.map(|x|x*e.area_m2)),
                support_cells:e.support.iter().map(|s|s.cell).collect()});
        }
        Ok(SupplySupport{cells:cells.into_values().collect(),targets,
            checked_native_denominators:checked,missing_required_neighbours:0})
    }

impl SupplySupport {
    pub fn to_json(&self)->Value {
        json!({"cells":self.cells.iter().map(|c|json!({"cell":c.cell,"live":c.live,"area_m2":c.area_m2,
            "deposited_area_m2":c.deposited_area_m2,"neighbours":c.neighbours.iter().map(|n|json!({"carrier":n.carrier,
                "weight":n.weight,"area_m2":n.area_m2})).collect::<Vec<_>>()})).collect::<Vec<_>>(),
            "targets":self.targets.iter().map(|t|json!({"carrier":t.carrier,"kind":t.kind,"time_my":t.time_my,
                "area_m2":t.area_m2,"position":t.position,"arc_m3":t.arc_m3,"plume_m3":t.plume_m3,
                "birth_import_m3":t.birth_import_m3,"support_cells":t.support_cells})).collect::<Vec<_>>(),
            "checked_native_denominators":self.checked_native_denominators,"missing_required_neighbours":self.missing_required_neighbours})
    }
}
