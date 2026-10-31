//! Adapter over the native selected-carrier recorder. No new simulation.
//! Cell forcing coverage is a necessary diagnostic, not proof of historical
//! material/contact ancestry. Missing deposited area stays explicitly unknown.
use crate::carrier_history::StackHistory;
use serde_json::{json,Value};
use std::collections::BTreeMap;
type Result<T> = std::result::Result<T,String>;
fn number(v:&Value)->Result<f64> {
    let s=v.as_str().ok_or("history float must be encoded binary64")?;
    let x=f64::from_bits(u64::from_str_radix(s,16).map_err(|_|"history float encoding")?);
    if !x.is_finite() {return Err("nonfinite history value".into())} Ok(x)
}
/// Independently sum the selected parcels' deposited area at every captured
/// cell/step and compare to the recorder's full native projection denominator.
/// Does not declare a body or assign the unknown remainder to any material.
pub fn audit(history:&StackHistory,steps:u64,dt_my:f64)->Result<Value> {
    history.validate_schedule(steps,dt_my)?;
    let archive=history.to_json();
    let events=archive["events"].as_array().ok_or("missing history events")?;
    // Keys use positive model-time bits; equality is exact, not a fuzzy join.
    let mut cells:BTreeMap<(u64,u64),(f64,f64)>=BTreeMap::new();
    let mut kinds=BTreeMap::<u64,usize>::new();
    for e in events {
        let kind=e["kind"].as_u64().ok_or("event kind")?;*kinds.entry(kind).or_default()+=1;
        if kind!=1 {continue}
        let time=number(&e["end_my"])?;let area=number(&e["area_m2"])?;
        for c in e["support"].as_array().ok_or("operation support")? {
            let cell=c["cell"].as_u64().ok_or("cell ID")?;
            let total=number(&c["deposited_area_m2"])?;
            let row=cells.entry((time.to_bits(),cell)).or_insert((0.,total));
            if row.1.to_bits()!=total.to_bits() {return Err("inconsistent native cell area denominator".into())}
            row.0+=area*number(&c["weight"])?;
        }
    }
    let mut operations=Vec::new();let mut missing_operations=0usize;
    for e in events {
        if e["kind"]!=1 {continue}
        let end=number(&e["end_my"])?;let mut gaps=Vec::new();let mut max_missing=0.0f64;
        for c in e["support"].as_array().ok_or("operation support")? {
            if c["live"]!=true {continue}
            let cell=c["cell"].as_u64().ok_or("cell ID")?;
            let (selected,total)=cells[&(end.to_bits(),cell)];
            if total<=0. {return Err("live support has no native deposited area".into())}
            // At most 256 nonnegative terms, binary64. This allowance exceeds
            // their summation error and cannot hide a geological fraction.
            let tolerance=1e-9*total;
            if selected>total+tolerance {return Err("selected deposited area exceeds native total".into())}
            let missing=(total-selected).max(0.);
            if missing>tolerance {
                max_missing=max_missing.max(missing/total);
                gaps.push(json!({"cell":cell,"native_deposited_area_m2":total,
                    "selected_deposited_area_m2":selected,"unrecorded_area_m2":missing,
                    "unrecorded_fraction":missing/total}));
            }
        }
        if !gaps.is_empty() {missing_operations+=1;}
        operations.push(json!({"carrier":e["id"],"start_my":number(&e["start_my"])?,"end_my":end,
            "before_cover_count":e["before"]["cover"].as_array().ok_or("before stack")?.len(),
            "after_cover_count":e["after"]["cover"].as_array().ok_or("after stack")?.len(),
            "changes":e["evolution"]["changes"].as_array().map(|a|a.len()),
            "max_unrecorded_cell_area_fraction":max_missing,"cell_coverage_gaps":gaps}));
    }
    Ok(json!({"schema":"vz.body-history-audit/v1","steps":steps,"dt_my":dt_my,
        "selected":history.config().selected,"event_kinds":kinds,"charged_bytes":archive["charged_bytes"],
        "operations":operations,"operations_with_cell_coverage_gaps":missing_operations,
        "coverage_basis":"selected parcel area times native tent weight versus recorded full deposited area",
        "relative_numerical_tolerance":1e-9,
        "material_neighbour_lineage":"not established by cell area equality",
        "terrain_join":"pending","body_qualification":"not performed"}))
}

use crate::stratigraphy::{StrataStack,StrataLayer,ParcelAmounts};
use crate::layer_provenance::{StackProvenance,Recorder};
use crate::layer_evolution::LayerChange;
fn unsigned(v:&Value)->Result<u64>{v.as_u64().ok_or("expected native unsigned field".into())}
fn boolean(v:&Value)->Result<bool>{v.as_bool().ok_or("expected native boolean field".into())}
fn stack(v:&Value)->Result<StrataStack>{
    Ok(StrataStack{basement_m:number(&v["basement_m"])?,basement_lith:unsigned(&v["basement_lith"])? as u8,
        basement_age_my:number(&v["basement_age_my"])?,last_exposure_my:number(&v["last_exposure_my"])?,
        rock_uplift_m:number(&v["rock_uplift_m"])?,elev_ref_m:number(&v["elev_ref_m"])?,
        cover:v["cover"].as_array().ok_or("native cover array")?.iter().map(|l|Ok(StrataLayer{
            thickness_m:number(&l["thickness_m"])?,age_base_my:number(&l["age_base_my"])?,age_top_my:number(&l["age_top_my"])?,
            max_burial_m:number(&l["max_burial_m"])?,lith:unsigned(&l["lith"])? as u8,env:unsigned(&l["env"])? as u8,
            grade:unsigned(&l["grade"])? as u8,flags:unsigned(&l["flags"])? as u8})).collect::<Result<Vec<_>>>()?})
}
fn amounts(v:&Value)->Result<ParcelAmounts>{Ok(ParcelAmounts{
    erosion:number(&v["erosion"])?,deposit:number(&v["deposit"])?,class:unsigned(&v["class"])? as u8,
    arc:number(&v["arc"])?,plume:number(&v["plume"])?,thinning:number(&v["thinning"])?,basement_add:number(&v["basement_add"])?,
    short_pos:number(&v["short_pos"])?,short_neg:number(&v["short_neg"])?,belt_rest:number(&v["belt_rest"])?,
    land:boolean(&v["land"])?,warm_shelf:boolean(&v["warm_shelf"])?})}
// Native reconciliation may temporarily carry a negative basement deficit.
// The unchanged provenance observer represents that as zero material. Preserve
// the signed native stack; validate only its material quantity against max(0,b).
fn validate_observed_material(p:&StackProvenance,st:&StrataStack)->Result<()> {
    p.basement.validate(st.basement_m.max(0.))?;
    if p.cover.len()!=st.cover.len() {return Err("observed cover shape mismatch".into())}
    for (a,b) in p.cover.iter().zip(&st.cover) {a.validate(b.thickness_m.max(0.))?;}
    Ok(())
}
fn source_height(p:&StackProvenance,episode:&str)->f64 {
    std::iter::once(&p.basement).chain(&p.cover).flat_map(|l|&l.contributions).filter(|c|c.episode==episode).map(|c|c.thickness_m).sum()
}
fn mixed(p:&StackProvenance,episode:&str)->bool {
    std::iter::once(&p.basement).chain(&p.cover).any(|l|l.contributions.iter().any(|c|c.episode==episode&&c.thickness_m>0.)
        &&l.contributions.iter().any(|c|c.episode!=episode&&c.thickness_m>0.))
}
fn top(p:&StackProvenance,episode:&str)->bool {p.cover.last().is_some_and(|l|l.contributions.iter().any(|c|c.episode==episode&&c.thickness_m>0.))}
#[derive(Debug,Clone)]
struct PhaseStock { interval_my:[f64;2], retained_m:f64 }
fn scale_phases(p:&mut [PhaseStock],before:f64,after:f64)->Result<()> {
    if before<=0. {if after==0. {return Ok(())}return Err("source phase appears without emplacement".into())}
    for s in p {s.retained_m*=after/before;}Ok(())
}
fn truncate_phases(p:&mut [PhaseStock],amount:f64,from_top:bool)->Result<()> {
    let mut need=amount;
    let indices:Vec<usize>=if from_top {(0..p.len()).rev().collect()}else{(0..p.len()).collect()};
    for i in indices {let take=need.min(p[i].retained_m);p[i].retained_m-=take;need-=take;}
    if need>1e-9*amount.max(1.) {return Err("ordered source phase removal does not close".into())}Ok(())
}
#[derive(Debug,Clone)]
pub struct SourceBook {
    pub carrier:u64,
    pub episode:String,
    pub area_m2:f64,
    pub position:[f64;3],
    pub final_stack:StrataStack,
    pub provenance:StackProvenance,
    pub report:Value,
}
/// Reapply the existing provenance observer to actual recorded native mutations.
/// This supplies forward source accounts before geometric subdivision. It does
/// not substitute a different erosion or provenance law. Mixed target layers
/// explicitly withhold ordered-removal compatibility rather than relabel them.
pub fn source_books(history:&StackHistory,targets:&[(u64,String)])->Result<Vec<SourceBook>> {
    if !history.config().layer_evolution {return Err("source book requires evolution".into())}
    let archive=history.to_json();let events=archive["events"].as_array().ok_or("events")?;
    let mut output=Vec::new();
    for (id,episode) in targets {
        if !history.selected(*id) {return Err("target was not captured".into())}
        let first=events.iter().find(|e|e["id"]==*id).ok_or("missing target initialization")?;
        if first["kind"]!=0 {return Err("source book birth adapter not qualified for this target".into())}
        let initial=stack(&first["after"])?;let mut provenance=crate::layer_provenance::initial(&initial,*id);
        let initial_h=source_height(&provenance,episode);
        let mut ordered=if initial_h>0. {let c=std::iter::once(&provenance.basement).chain(&provenance.cover).flat_map(|p|&p.contributions).find(|c|c.episode==*episode).ok_or("initial source")?;vec![PhaseStock{interval_my:[c.start_my,c.end_my],retained_m:initial_h}]}else{vec![]};
        let mut homogeneous=ordered.clone();let mut temporal_differences=Vec::new();
        let mut final_stack=initial;
        let area=number(&first["area_m2"])?;let mut position=[0.;3];let mut terms=BTreeMap::<String,[f64;2]>::new();
        let mut changes=Vec::new();let mut mixtures=Vec::new();let mut top_creation:Option<[f64;2]>=None;
        let mut max_base_depth=0.0f64;let mut phases=Vec::new();let mut final_land=None;let mut final_support_land=None;let mut initial_exposed=top(&provenance,episode);
        if mixed(&provenance,episode) {mixtures.push(json!({"time_my":0.,"kind":"initial"}));}
        for e in events.iter().filter(|e|e["id"]==*id) {
            let end=number(&e["end_my"])?;let start=number(&e["start_my"])?;
            position=[number(&e["position"][0])?,number(&e["position"][1])?,number(&e["position"][2])?];
            if number(&e["area_m2"])?!=area {return Err("carrier area changed in history".into())}
            if e["kind"]==3 {return Err("requested source owner retired".into())}
            if e["kind"]!=1 {continue}
            let input=amounts(&e["amounts"])?;
            final_land=Some(input.land);final_support_land=Some(e["support"].as_array().ok_or("native support")?.iter().all(|c|c["live"]!=true||c["land"]==true));
            let mut recorder=Recorder{state:provenance,carrier:*id,amounts:input,params:*history.params(),t0:start,t1:end};
            for (index,c) in e["evolution"]["changes"].as_array().ok_or("mutation trace")?.iter().enumerate() {
                let change=LayerChange{kind:unsigned(&c["kind"])? as u8,links:serde_json::from_value(c["links"].clone()).map_err(|e|e.to_string())?,before:stack(&c["before"])?,after:stack(&c["after"])?};
                validate_observed_material(&recorder.state,&change.before)?;
                let before=source_height(&recorder.state,episode);let was_top=top(&recorder.state,episode);
                if mixed(&recorder.state,episode) {mixtures.push(json!({"interval_my":[start,end],"mutation":index,"side":"before"}));}
                recorder.capture(&change);validate_observed_material(&recorder.state,&change.after)?;
                let after=source_height(&recorder.state,episode);let is_top=top(&recorder.state,episode);let delta=after-before;
                if delta>0. && matches!(change.kind,8|12|13) {
                    let p=PhaseStock{interval_my:[start,end],retained_m:delta};ordered.push(p.clone());homogeneous.push(p);
                } else if delta!=0. {
                    scale_phases(&mut homogeneous,before,after)?;
                    if delta<0. && matches!(change.kind,0|1|2|15) {truncate_phases(&mut ordered,-delta,change.kind!=15)?;}
                    else {scale_phases(&mut ordered,before,after)?;}
                }
                let disagreement=ordered.iter().zip(&homogeneous).map(|(a,b)|(a.retained_m-b.retained_m).abs()).fold(0f64,f64::max);
                if disagreement>1e-9*after.max(1.) {temporal_differences.push(json!({"interval_my":[start,end],"mutation":index,"max_phase_retention_difference_m":disagreement}));}

                if mixed(&recorder.state,episode) {mixtures.push(json!({"interval_my":[start,end],"mutation":index,"side":"after"}));}
                if is_top && (!was_top || (delta<0. && matches!(change.kind,0|1|2|15)) || (delta>0. && matches!(change.kind,8|12|13))) {top_creation=Some([start,end]);initial_exposed=false;}
                let name=crate::layer_evolution::NAMES[change.kind as usize];let t=terms.entry(name.into()).or_default();
                if delta>=0. {t[0]+=delta;} else {t[1]-=delta;}
                if delta!=0. {
                    changes.push(json!({"interval_my":[start,end],"mutation":index,"kind":name,"before_m":before,"after_m":after,"delta_m":delta,
                        "native_amounts":{"erosion_m":input.erosion,"deposit_m":input.deposit,"arc_m":input.arc,"plume_m":input.plume,"thinning_m":input.thinning,"shortening_m":input.short_pos}}));
                    if delta>0. && matches!(change.kind,8|12|13) {phases.push(json!({"interval_my":[start,end],"kind":name,"added_thickness_m":delta,"not_individual_eruptions":true}));}
                }
                for (i,l) in recorder.state.cover.iter().enumerate() {
                    if l.contributions.iter().any(|c|c.episode==*episode&&c.thickness_m>0.) {
                        max_base_depth=max_base_depth.max(change.after.cover[i..].iter().map(|l|l.thickness_m).sum());
                    }
                }
            }
            provenance=recorder.state;final_stack=stack(&e["after"])?;crate::layer_provenance::validate(&provenance,&final_stack)?;
        }
        let retained=source_height(&provenance,episode);if retained<=0. {return Err("target source is absent at final state".into())}
        let additions:f64=terms.values().map(|v|v[0]).sum();let removals:f64=terms.values().map(|v|v[1]).sum();
        let residual=initial_h+additions-removals-retained;let gross=initial_h+additions+removals;
        let bound=32.*f64::EPSILON*(changes.len()+1) as f64*gross.max(1.);
        if residual.abs()>bound {return Err("forward source account does not close".into())}
        let report=json!({"schema":"vz.native-source-account/v1","carrier":id,"episode":episode,"area_m2":area,
            "quantity_basis":"native thin-shell area times recorded contribution thickness",
            "initial_m":initial_h,"positive_operation_total_m":additions,"negative_operation_total_m":removals,"retained_m":retained,
            "initial_m3":initial_h*area,"retained_m3":retained*area,"balance_residual_m":residual,"roundoff_bound_m":bound,
            "operation_totals_m":terms,"changes":changes,"emplacement_phases":phases,
            "mixed_target_boundaries":mixtures,"ordered_removal_compatibility":if mixtures.is_empty(){"single_episode_at_all_observed_boundaries"}else{"unresolved_mixed_layer_order"},
            "intra_episode_pulses":"compressed model intervals are not resolved eruptions; no thermal scenario admitted",
            "ordered_native_addition_phase_retention":ordered.iter().map(|p|json!({"interval_my":p.interval_my,"retained_m":p.retained_m})).collect::<Vec<_>>(),"homogeneous_extension_phase_retention":homogeneous.iter().map(|p|json!({"interval_my":p.interval_my,"retained_m":p.retained_m})).collect::<Vec<_>>(),
            "within_episode_temporal_differences":temporal_differences,
            "temporal_comparison_scope":"native episode totals checked exactly; homogeneous subinterval fractions are a diagnostic extension, not a field present in the retained compressed episode; no subepisode ancestry is qualified",

            "max_observed_layer_base_depth_m":max_base_depth,"final_top":top(&provenance,episode),
            "final_top_face_creation_interval_my":top_creation,"initial_exposed_with_unknown_pre_t0":initial_exposed,
            "final_native_land":final_land,"final_live_support_all_land":final_support_land,
            "pre_t0":"declared initial condition only; no simulated earlier history",
            "thermal":"unresolved","terrain_join":"pending"});
        output.push(SourceBook{carrier:*id,episode:episode.clone(),area_m2:area,position,final_stack,provenance,report});
    }
    Ok(output)
}

#[derive(Debug,Clone)]
pub struct CarrierEndpoints {
    pub id:u64,
    pub initial_position:[f64;3],
    pub final_position:[f64;3],
    pub area_m2:f64,
    pub initial_stack:StrataStack,
    pub final_stack:StrataStack,
}
/// Native endpoints for the conditional transported-initial-cell partition.
/// V1 refuses converted/rift births rather than inventing their earlier domain.
pub fn initial_carrier_endpoints(history:&StackHistory)->Result<BTreeMap<u64,CarrierEndpoints>> {
    let archive=history.to_json();let events=archive["events"].as_array().ok_or("events")?;
    let mut first=BTreeMap::new();let mut last=BTreeMap::new();
    for e in events {if let Some(id)=e["id"].as_u64(){first.entry(id).or_insert(e);last.insert(id,e);}}
    let mut result=BTreeMap::new();
    for (id,a) in first {
        let b=last[&id];
        if a["kind"]!=0 || b["kind"]==3 {return Err("initial-cell partition requires surviving initial carriers".into())}
        let pos=|e:&Value|->Result<[f64;3]>{Ok([number(&e["position"][0])?,number(&e["position"][1])?,number(&e["position"][2])?])};
        let area=number(&a["area_m2"])?;if area!=number(&b["area_m2"])? {return Err("partition area changed".into())}
        result.insert(id,CarrierEndpoints{id,initial_position:pos(a)?,final_position:pos(b)?,area_m2:area,initial_stack:stack(&a["after"])?,final_stack:stack(&b["after"])?});
    }
    Ok(result)
}

/// Adapter identity, separate from the unchanged native recorder codec.
pub fn source_fingerprint()->String {pf_contracts::sim_module::digest_bytes(include_bytes!("body_history.rs"))}

/// Quantifies the native projection/area convention without promoting grid means
/// to conserved parents. Outer incomplete selected support remains explicit.
pub fn cell_mean_diagnostic(history:&StackHistory)->Result<Value>{
    let archive=history.to_json();let events=archive["events"].as_array().ok_or("events")?;
    let time=history.latest_time();let mut cells=BTreeMap::<u64,(f64,f64,f64,f64,usize)>::new();
    for e in events.iter().filter(|e|e["kind"]==1) {
        if number(&e["end_my"])?!=time{continue}
        let area=number(&e["area_m2"])?;let height=number(&e["total_after_m"])?;
        for c in e["support"].as_array().ok_or("support")? {
            let id=unsigned(&c["cell"])?;let cell_area=number(&c["area_m2"])?;let total=number(&c["deposited_area_m2"])?;
            let r=cells.entry(id).or_insert((cell_area,total,0.,0.,0));
            if r.0!=cell_area||r.1!=total{return Err("cell mean diagnostic input mismatch".into())}
            let d=area*number(&c["weight"])?;r.2+=d;r.3+=d*height;r.4+=1;
        }
    }
    let rows:Vec<_>=cells.into_iter().map(|(cell,(a,d,selected,volume,n))|{
        let complete=d>0.&&(selected-d).abs()<=1e-9*d;
        let mean=complete.then(||volume/d);
        json!({"cell":cell,"cell_area_m2":a,"native_deposited_area_m2":d,"captured_deposited_area_m2":selected,
            "captured_parcels":n,"complete_material_support":complete,"deposited_area_over_cell_area":d/a,
            "native_tent_mean_total_column_m":mean,"captured_native_deposited_volume_m3":volume,
            "cell_area_times_mean_minus_parcel_volume_m3":mean.map(|h|a*h-volume),
            "equivalent_thickness_gap_m":mean.map(|h|(a-d)*h/a),
            "interpretation":"projection diagnostic only; cell area times kernel mean is not another conserved parent"})
    }).collect();
    Ok(json!({"schema":"vz.native-cell-mean-diagnostic/v1","time_my":time,"rows":rows,"native_parent":"Lagrangian parcel ledger","quantity_basis":"native thin-shell total crust+sediment column, not mineral mass","outer_incomplete_support":"unknown mean; retained, not silently omitted"}))
}
