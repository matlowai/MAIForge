//! Stage-input surface dependency certificates. Final heights never enter the certificate.
//! A dirty computation domain conservatively influences its entire published rectangle.
//! No finite eikonal/fill attenuation radius or final-height difference is used as a mask.
use std::{cell::RefCell, collections::{BTreeMap, BTreeSet}, path::Path};
use serde_json::{json, Value};
use crate::{t1b_rivernet::{Curve, LatPrep, Network, RiverNetParamsV1}, t1b_coast::CoastCtx};
type R<T> = Result<T,String>;
pub const RULE: &str = "esc_window_input_domain_v3";
static CODE: std::sync::OnceLock<String> = std::sync::OnceLock::new();
fn code_digest() -> R<String> {
    if let Some(s) = CODE.get() { return Ok(s.clone()); }
    let path = std::env::current_exe().map_err(|e|e.to_string())?;
    let bytes = std::fs::read(path).map_err(|e|e.to_string())?;
    let s = format!("blake3:{}",blake3::hash(&bytes).to_hex()); let _ = CODE.set(s.clone()); Ok(s)
}
struct WindowKey { x:i64, y:i64, root:String, layer:String }
struct Binding { layer:String, reference:Option<String> }
struct Geometry { h:blake3::Hasher, counts:[u64;4] }
thread_local! {
    static GEOMETRY: RefCell<Option<Geometry>> = const { RefCell::new(None) };
    static BINDING: RefCell<Option<Binding>> = const { RefCell::new(None) };
    static WINDOW: RefCell<Option<WindowKey>> = const { RefCell::new(None) };
    static RAW: RefCell<BTreeMap<(i64,i64),crate::t1b_window::RawWindow>> = const { RefCell::new(BTreeMap::new()) };
}
pub fn bind_layer(dir:&Path, layer:&str, p:&RiverNetParamsV1)->R<()> {
    if std::env::var_os("VZ_ESC_DEPS_OUT").is_none() {return Ok(());}
    let reference=if std::env::var_os("VZ_ESC_DEPS_REF").is_some() {
        let r=crate::t1b_esc_local::reference(p).ok_or("dependency comparison needs esc_ref=<reference layer>")?;
        let (rd,_,_,_)=crate::t1b_rivernet::load_layer_co(r)?;
        let v:Value=serde_json::from_slice(&std::fs::read(dir.join("rivernet.json")).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
        if let Some(expected)=v.pointer("/header/lane_es/incremental/reference_digest").and_then(Value::as_str) {
            if expected!=rd {return Err("dependency reference is not the incremental layer's parent".into());}
        }
        Some(rd)
    }else{None};
    BINDING.with(|b|*b.borrow_mut()=Some(Binding{layer:layer.into(),reference}));Ok(())
}
pub fn begin_window(x:i64,y:i64,root:&str,layer:&str) {
    let on=std::env::var_os("VZ_ESC_DEPS_OUT").is_some();
    WINDOW.with(|w| *w.borrow_mut() = on.then(||WindowKey{x,y,root:root.into(),layer:layer.into()}));
    GEOMETRY.with(|g| *g.borrow_mut() = on.then(||Geometry{h:blake3::Hasher::new(),counts:[0;4]}));
}
fn floats(h:&mut blake3::Hasher,v:&[f64]) { for x in v { h.update(&x.to_bits().to_le_bytes()); } }
fn f32s(h:&mut blake3::Hasher,v:&[f32]) { for x in v { h.update(&x.to_bits().to_le_bytes()); } }
fn u32s(h:&mut blake3::Hasher,v:&[u32]) { for x in v { h.update(&x.to_le_bytes()); } }
fn bools(h:&mut blake3::Hasher,v:&[bool]) { for &x in v { h.update(&[x as u8]); } }
fn digest(h:blake3::Hasher)->String { format!("blake3:{}",h.finalize().to_hex()) }
fn overlaps(pts:&[[f64;4]],r:[f64;4])->bool {
    let mut b=[f64::INFINITY,f64::INFINITY,f64::NEG_INFINITY,f64::NEG_INFINITY];
    for p in pts {b[0]=b[0].min(p[0]);b[1]=b[1].min(p[1]);b[2]=b[2].max(p[0]);b[3]=b[3].max(p[1]);}
    b[2]>=r[0] && b[0]<=r[2] && b[3]>=r[1] && b[1]<=r[3]
}
#[allow(clippy::too_many_arguments)]
pub fn surface_inputs(lp:&LatPrep, net:&Network, curves:&[Curve], table:&BTreeMap<u32,f64>, p:&RiverNetParamsV1, cx:&CoastCtx,
    e:&[f32], z:&[f64], fixed:&[u8], floor:&[bool], known:&[bool], labels:&[u32], areas:&[f32], slopes:&[f32]) -> R<Option<Value>> {
    WINDOW.with(|w| {
        let key=w.borrow_mut().take(); let Some(k)=key.as_ref() else {return Ok(None)};
        if cx.gl.is_some() || curves.iter().any(|c|c.mq.is_some()) || cx.cf.as_ref().and_then(|c|c.ring_plain.as_ref()).is_some() {
            return Err("esc dependencies: glacier, swept MQ and ring-plain domains are not tracked".into());
        }
        let lat=lp.lat; let dx=lat.dx();
        // Same conservative raster box as evaluate(), including its upper rounding cell.
        let r=[lat.gx0 as f64*dx,lat.gy0 as f64*dx,(lat.gx0+lat.n as i64)as f64*dx,(lat.gy0+lat.n as i64)as f64*dx];
        // Carrier and fan-apex scans use rounded lattice indices without the raster's bbox
        // rejection. One post of support conservatively includes every half-post rounding case.
        let cr=[r[0]-dx,r[1]-dx,r[2]+dx,r[3]+dx];
        let mut parts=BTreeMap::<String,String>::new();
        let mut h=blake3::Hasher::new();
        h.update(k.root.as_bytes());h.update(code_digest()?.as_bytes());
        let env:BTreeMap<String,String> = std::env::vars().filter(|(k,_)| k.starts_with("VZ_") && !["VZ_ESC_DEPS_OUT","VZ_ESC_DEPS_REF","VZ_ESC_EVIDENCE","VZ_ESC_EVIDENCE_REUSE","VZ_ESC_PERTURB","VZ_ESC_LOCAL_DIAGNOSTIC"].contains(&k.as_str())).collect();
        h.update(&serde_json::to_vec(&env).map_err(|e|e.to_string())?);
        h.update(&serde_json::to_vec(p).map_err(|e|e.to_string())?);
        floats(&mut h,&cx.ct_steps);floats(&mut h,&cx.ct_beds.glob);floats(&mut h,&cx.ct_beds.glob_eta);
        for row in &cx.ct_beds.pk {h.update(&(row.len()as u64).to_le_bytes());for &(a,b) in row {floats(&mut h,&[a,b]);}}
        parts.insert("immutable_recipe_and_closure".into(),digest(h));
        let mut h=blake3::Hasher::new();let mut hc=BTreeSet::new();
        for i in 0..lat.nodes() {for j in 0..4 {if lp.hat_w[i][j]>0.0 {hc.insert(lp.hat_c[i][j]);}}}
        for c in hc {h.update(&c.to_le_bytes());floats(&mut h,&[table.get(&c).copied().unwrap_or(0.0),cx.flood.get(&c).copied().unwrap_or(0.0)]);}
        f32s(&mut h,e);parts.insert("local_E_F_hats".into(),digest(h));
        let mut h=blake3::Hasher::new();floats(&mut h,z);h.update(fixed);bools(&mut h,floor);bools(&mut h,known);u32s(&mut h,labels);f32s(&mut h,areas);f32s(&mut h,slopes);
        parts.insert("raster_floor_seeds".into(),digest(h));
        let mut h=blake3::Hasher::new();let mut ids:BTreeSet<u32>=labels.iter().copied().filter(|&x|x!=u32::MAX).collect();let mut ncurve=0;
        // Full relevant curves conservatively cover raster samples and tributary-fan apex scans.
        // The global maximum carrier belt is an actual input to the floor rule.
        let pad=curves.iter().filter_map(|c|c.carrier.as_ref().map(|x|x.1)).fold(0.0,f64::max);
        floats(&mut h,&[pad]);h.update(&[curves.iter().any(|c|c.carrier.is_some()) as u8]);
        for c in curves {if overlaps(&c.pts,cr) || c.carrier.as_ref().map(|a|overlaps(&a.0,cr)).unwrap_or(false) {
            ncurve+=1;h.update(&[c.code]);u32s(&mut h,&[c.link,c.label]);h.update(&(c.pts.len()as u64).to_le_bytes());
            for q in &c.pts {floats(&mut h,q);}
            if let Some((v,b))=&c.carrier {h.update(&[1]);floats(&mut h,&[*b]);for q in v {floats(&mut h,q);}}else{h.update(&[0]);}
            if c.label!=u32::MAX {ids.insert(c.label);}
        }}
        parts.insert("relevant_curves_and_carriers".into(),digest(h));
        let mut h=blake3::Hasher::new();for id in ids {let n=id as usize;if n>=net.kind.len(){return Err("esc dependencies: invalid source label".into());}u32s(&mut h,&[id,net.term[n]]);h.update(&[net.kind[n]]);floats(&mut h,&[net.rockmass[n]]);}
        parts.insert("source_terminal_and_rockmass".into(),digest(h));
        let mut h=blake3::Hasher::new();let mut nf=0;
        for f in &cx.fans {if f.x_m+f.radius_m>=r[0] && f.x_m-f.radius_m<=r[2] && f.y_m+f.radius_m>=r[1] && f.y_m-f.radius_m<=r[3] {
            nf+=1;h.update(&serde_json::to_vec(f).map_err(|e|e.to_string())?);
        }}
        parts.insert("relevant_coastal_fans".into(),digest(h));
        let mut h=blake3::Hasher::new();let mut nv=0;
        for (i,v) in cx.volcanoes.iter().enumerate() {
            let rad=v.r_influence_m.max(if cx.ve_units {v.ve.as_ref().map(|x|x.r_units_m).unwrap_or(0.0)}else{0.0});
            if v.vent_xy_m[0]+rad>=r[0] && v.vent_xy_m[0]-rad<=r[2] && v.vent_xy_m[1]+rad>=r[1] && v.vent_xy_m[1]-rad<=r[3] {
                nv+=1;h.update(&(i as u32).to_le_bytes());h.update(&serde_json::to_vec(v).map_err(|e|e.to_string())?);
                if let Some(&v)=cx.tf_vent_c1.get(i){floats(&mut h,&[v]);}
            }
        }
        h.update(&[cx.ve_units as u8]);parts.insert("relevant_volcanoes".into(),digest(h));
        let geometry=json!({"gx0":lat.gx0,"gy0":lat.gy0,"n":lat.n,"dx_mm":lat.dx_mm});
        let input_digest=crate::t1b_rivernet::RN_SCHEMA.to_string()+":"+&pf_contracts::artifact::digest_canonical(&json!({"rule":RULE,"geometry":geometry,"parts":parts})).map_err(|e|e.to_string())?;
        let mut record=json!({"rule":RULE,"layer_digest":k.layer,"window":[k.x,k.y],"geometry":geometry,"input_digest":input_digest,"parts":parts,
            "relevant_curves":ncurve,"relevant_fans":nf,"relevant_volcanoes":nv,"decision_stage":"prefix inputs plus native cap/shadow operand trace; no final-height comparisons",
            "scope":"entire published rectangle of every input-different computation domain; union on world posts"});
        if let Ok(dir)=std::env::var("VZ_ESC_DEPS_REF") {
            let bytes=std::fs::read(Path::new(&dir).join(format!("w_{}_{}.json",k.x,k.y))).map_err(|e|format!("reference dependency manifest: {e}"))?;
            let old:Value=serde_json::from_slice(&bytes).map_err(|e|e.to_string())?;
            let expected_digest=crate::t1b_rivernet::RN_SCHEMA.to_string()+":"+&pf_contracts::artifact::digest_canonical(&json!({"rule":old["rule"],"geometry":old["geometry"],"parts":old["parts"]})).map_err(|e|e.to_string())?;
            if old["input_digest"]!=expected_digest {return Err("reference input certificate digest mismatch".into());}
            BINDING.with(|b|->R<()> {let bb=b.borrow();let binding=bb.as_ref().ok_or("dependency comparison needs native layer binding")?;
                if binding.layer!=k.layer || old["layer_digest"].as_str()!=binding.reference.as_deref() {return Err("input certificate is not bound to the selected reference layer".into());}Ok(())})?;
            if old["rule"]!=RULE || old["geometry"]!=geometry {return Err("reference dependency rule/grid differs".into());}
            if old["parts"]["immutable_recipe_and_closure"]!=parts["immutable_recipe_and_closure"] {return Err("reference dependency immutable recipe, code or environment differs".into());}
            record["prefix_inputs_equal"]=json!(parts.iter().all(|(key,value)| old["parts"][key.as_str()]==*value));
            record["_reference"]=old;
        }
        Ok(Some(record))
    })
}

/// Native validation without retaining every material-field variant. The package deliberately
/// is not named window.json and cannot be loaded as a publishable window. Seam comparisons use
/// the actual complete RawWindow objects, not reconstructed or mocked fields.
pub fn publish_evidence(wr:&crate::t1b_window::WindowRun,out:&Path)->R<bool> {
    if std::env::var("VZ_ESC_EVIDENCE").as_deref()!=Ok("1"){return Ok(false);}
    std::fs::create_dir_all(out.join("fields")).map_err(|e|e.to_string())?;
    for (name,b) in &wr.blobs {if ["fields/height_mm.i32","fields/t0_cell.u32","fields/fixed.u8"].contains(&name.as_str()) {
        if let Ok(root)=std::env::var("VZ_ESC_EVIDENCE_REUSE") {
            let old=Path::new(&root).join(out.file_name().ok_or("evidence window name")?).join(name);
            let bytes=std::fs::read(&old).map_err(|e|format!("evidence reuse {}: {e}",old.display()))?;
            if bytes!=*b {return Err(format!("evidence reuse mismatch for {}: newly computed bytes differ",old.display()));}
            std::fs::hard_link(&old,out.join(name)).map_err(|e|format!("evidence hard link: {e}"))?;
        } else { std::fs::write(out.join(name),b).map_err(|e|e.to_string())?; }
    }}
    let fields:Vec<Value>=wr.blobs.iter().map(|(n,b)|json!({"file":n,"bytes":b.len(),"digest":format!("blake3:{}",blake3::hash(b).to_hex())})).collect();
    let v=json!({"kind":"nonpublishable native validation evidence","header":wr.header,"digest":wr.raw.digest,"complete_field_hashes":fields});
    std::fs::write(out.join("window_evidence.json"),serde_json::to_vec_pretty(&v).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
    RAW.with(|cell|->R<()> {
        let mut map=cell.borrow_mut();let (x,y)=(wr.raw.lat.rx,wr.raw.lat.ry);
        let seams=out.parent().ok_or("evidence parent")?.join("seams");std::fs::create_dir_all(&seams).map_err(|e|e.to_string())?;
        for (dx,dy) in [(-1,0),(1,0),(0,-1),(0,1)] {if let Some(other)=map.get(&(x+dx,y+dy)) {
            let (s,_)=crate::t1b_checks::seam_compare(other,&wr.raw,wr.header.recipe.params.seam_half_width_mm/crate::t1b::T1B_SPACING_MM)?;
            std::fs::write(seams.join(format!("s_{}_{}__{}_{}.json",x+dx,y+dy,x,y)),serde_json::to_vec_pretty(&s).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
        }}
        map.insert((x,y),wr.raw.clone());Ok(())
    })?;
    Ok(true)
}

#[cfg(test)] mod tests {
    #[test] fn rounded_fan_points_are_inside_the_dependency_envelope() {
        let lat=crate::t1b_rivernet::EvalLat{gx0:0,gy0:0,n:4,dx_mm:32_000};
        let pts=[[-8.0,32.0,100.0,1.0],[-8.0,64.0,99.0,2.0]];
        assert!(lat.idx((pts[0][0]/32.0_f64).round()as i64,1).is_some());
        assert!(!super::overlaps(&pts,[0.0,0.0,128.0,128.0]));
        assert!(super::overlaps(&pts,[-32.0,-32.0,160.0,160.0]));
        assert!(!super::overlaps(&[[-40.0,32.0,0.0,0.0]],[-32.0,-32.0,160.0,160.0]));
    }
}


/// Observe exact native cap operands, after native bbox rejection and before min/clamp.
pub fn tracing()->bool {GEOMETRY.with(|g|g.borrow().is_some())}
pub fn cap_begin(stage:u8,x:f64,y:f64,dneed:f64) {GEOMETRY.with(|g|if let Some(g)=g.borrow_mut().as_mut(){g.h.update(&[0,stage]);floats(&mut g.h,&[x,y,dneed]);g.counts[0]+=1;});}
pub fn cap_segment(s:&[f64;4]) {GEOMETRY.with(|g|if let Some(g)=g.borrow_mut().as_mut(){g.h.update(&[1]);floats(&mut g.h,s);g.counts[1]+=1;});}
pub fn geometry_end(){GEOMETRY.with(|g|if let Some(g)=g.borrow_mut().as_mut(){g.h.update(&[2]);});}
pub fn shadow_begin(a:[f64;2],b:[f64;2]){GEOMETRY.with(|g|if let Some(g)=g.borrow_mut().as_mut(){g.h.update(&[3]);floats(&mut g.h,&[a[0],a[1],b[0],b[1]]);g.counts[2]+=1;});}
#[inline(always)]
fn shadow_possible(a:[f64;2],b:[f64;2],s:&[f64;4])->bool {
    let pad=32.0;
    s[0].max(s[2])>=a[0].min(b[0])-pad && s[0].min(s[2])<=a[0].max(b[0])+pad && s[1].max(s[3])>=a[1].min(b[1])-pad && s[1].min(s[3])<=a[1].max(b[1])+pad
}
#[inline(always)]
pub fn shadow_segment(a:[f64;2],b:[f64;2],s:&[f64;4]) {if shadow_possible(a,b,s){GEOMETRY.with(|g|if let Some(g)=g.borrow_mut().as_mut(){g.h.update(&[4]);floats(&mut g.h,s);g.counts[3]+=1;});}}
pub fn finish_inputs(mut record:Value)->R<Value> {
    let g=GEOMETRY.with(|g|g.borrow_mut().take()).ok_or("missing native geometry input trace")?;
    record["parts"]["native_corridor_operands"]=json!(digest(g.h));
    record["corridor_trace_counts"]=json!({"cap_queries":g.counts[0],"cap_segment_operands":g.counts[1],"shadow_queries":g.counts[2],"shadow_segment_operands":g.counts[3]});
    let d=crate::t1b_rivernet::RN_SCHEMA.to_string()+":"+&pf_contracts::artifact::digest_canonical(&json!({"rule":RULE,"geometry":record["geometry"],"parts":record["parts"]})).map_err(|e|e.to_string())?;
    record["input_digest"]=json!(d);
    if let Some(old)=record.as_object_mut().unwrap().remove("_reference") {
        let changed:Vec<String>=record["parts"].as_object().unwrap().iter().filter(|(k,v)|old["parts"][k.as_str()]!=**v).map(|(k,_)|k.clone()).collect();
        record["influenced"]=json!(old["input_digest"]!=record["input_digest"]);
        record["changed_input_groups"]=json!(changed);record["reference_input_digest"]=old["input_digest"].clone();record["reference_layer_digest"]=old["layer_digest"].clone();
    }
    let dir=std::env::var("VZ_ESC_DEPS_OUT").map_err(|e|e.to_string())?;std::fs::create_dir_all(&dir).map_err(|e|e.to_string())?;
    let x=record["window"][0].as_i64().ok_or("window x")?;let y=record["window"][1].as_i64().ok_or("window y")?;
    std::fs::write(Path::new(&dir).join(format!("w_{x}_{y}.json")),serde_json::to_vec_pretty(&record).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
    Ok(record)
}
#[cfg(test)] mod geometry_tests {
    #[test] fn far_segments_cannot_shadow_but_crossing_segments_are_kept() {
        assert!(!super::shadow_possible([0.0,0.0],[32.0,32.0],&[10000.0,0.0,10000.0,32.0]));
        assert!(super::shadow_possible([0.0,16.0],[1000.0,16.0],&[500.0,0.0,500.0,32.0]));
    }
}
