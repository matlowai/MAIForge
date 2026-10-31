//! E1/E3 source-solid ownership for the native cone process.
//!
//! A compact tile is a content-bound input, not a budget. Disjoint half-open
//! post owners restrict the existing envelope. Geometry above its retained top
//! is a named zero-stock account. It cannot be detached as invented rock.
use crate::t1b_surface_completion::Rect;
use crate::t1b_surface_domain::content_digest;
use crate::t1b_litho::EnvelopeConfig;
use serde::{Deserialize, Serialize};
use std::{collections::{BTreeMap, BTreeSet}, path::Path};
type R<T> = Result<T, String>;

pub const RULE: &str = "native_envelope_stock_v1: content-bound half-open native material posts; spherical post-prism support using E3 quadrature; native E1 cell/record owns each retained interval; finite occupancy restricted above source sea datum; below-datum and unoccupied parent budget remain reserve; construction excess is zero stock; proxy requested removal never admits stock; no source renaming, budget fitting or inferred accumulated water source";

pub struct CompactInput {
    pub header: serde_json::Value,
    pub rect: Rect,
    pub digest: String,
    pub height_mm: Vec<i32>,
    pub fixed: Vec<u8>,
    pub links: Vec<u32>,
    pub cells: Vec<u32>,
    pub basins: Vec<u32>,
    pub runoff: Vec<f32>,
    pub proxy_depth: Vec<f64>,
    pub material: Vec<f64>,
    pub exposure: Vec<f64>,
    pub proxy_floor: Vec<f64>,
}
impl CompactInput {
    pub fn load(path: &Path) -> R<Self> {
        let mut header: serde_json::Value = serde_json::from_slice(&std::fs::read(path.join("material.json")).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
        let digest = header.as_object_mut().ok_or("material header object")?.remove("digest").and_then(|v|v.as_str().map(str::to_owned)).ok_or("material digest missing")?;
        if content_digest(&header).map_err(|e|e.to_string())? != digest { return Err("material header digest mismatch".into()); }
        if header["schema"] != "vz.native-physical-material/v1" || header["spacing_mm"] != 32000 || header["support_m2"] != 1024.0 { return Err("unsupported material sampling contract".into()); }
        for key in ["native_context_pin", "source_digest", "rn_recipe_digest", "layer_digest"] {
            if header[key].as_str().is_none_or(str::is_empty) { return Err(format!("material missing {key}")); }
        }
        let rect: Rect = serde_json::from_value(header["rect"].clone()).map_err(|e|e.to_string())?;
        let n=rect.count().map_err(|e|e.to_string())?;
        let half=vz_tilekey::CHART_HALF_EXTENT_MM/32000;
        if rect.gx0%512!=0 || rect.gy0%512!=0 || rect.gx0 < -half || rect.gy0 < -half || rect.gx0>half || rect.gy0>half
            || rect.width!=(half+1-rect.gx0).min(512) as usize || rect.height!=(half+1-rect.gy0).min(512) as usize { return Err("material is not a canonical half-open owner tile".into()); }
        let mut fields=BTreeMap::new();
        for d in header["fields"].as_array().ok_or("material fields")? {
            let name=d["file"].as_str().ok_or("material field name")?;
            if Path::new(name).components().count()!=1 || !matches!(Path::new(name).components().next(),Some(std::path::Component::Normal(_))) { return Err("material field path escapes tile".into()); }
            let bytes=std::fs::read(path.join(name)).map_err(|e|e.to_string())?;
            if d["bytes"].as_u64()!=Some(bytes.len() as u64) || d["digest"]!=blake3::hash(&bytes).to_hex().as_str() { return Err(format!("material {name} size/hash mismatch")); }
            if fields.insert(name.to_owned(),bytes).is_some(){return Err("duplicate material field".into());}
        }
        let take=|name:&str,size:usize|->R<&[u8]> {let b=fields.get(name).ok_or_else(||format!("missing material {name}"))?;if b.len()!=n*size{return Err(format!("material {name} dimensions"));}Ok(b)};
        let u32s=|name|->R<Vec<u32>>{Ok(take(name,4)?.chunks_exact(4).map(|c|u32::from_le_bytes(c.try_into().unwrap())).collect())};
        let f64s=|name|->R<Vec<f64>>{Ok(take(name,8)?.chunks_exact(8).map(|c|f64::from_le_bytes(c.try_into().unwrap())).collect())};
        let out=Self{rect,digest,height_mm:take("height_mm.i32",4)?.chunks_exact(4).map(|c|i32::from_le_bytes(c.try_into().unwrap())).collect(),
            fixed:take("process_geometry_fixed.u8",1)?.to_vec(),links:u32s("process_native_link.u32")?,cells:u32s("native_t0_cell.u32")?,basins:u32s("native_t0_basin.u32")?,
            runoff:take("runoff_weight.f32",4)?.chunks_exact(4).map(|c|f32::from_le_bytes(c.try_into().unwrap())).collect(),
            proxy_depth:f64s("proxy_stock_depth_m.f64")?,material:f64s("proxy_material_m.f64")?,exposure:f64s("proxy_exposure_fraction.f64")?,proxy_floor:f64s("proxy_floor_m.f64")?,header};
        for i in 0..n {
            if ![0,1,2,4,5].contains(&out.fixed[i]) || !out.runoff[i].is_finite() || out.runoff[i]<0.0
                || !out.proxy_depth[i].is_finite() || out.proxy_depth[i]<0.0 || !out.material[i].is_finite() || out.material[i]<0.0
                || !out.exposure[i].is_finite() || !(0.0..=1.0).contains(&out.exposure[i]) || !out.proxy_floor[i].is_finite() || out.proxy_floor[i]>out.height_mm[i] as f64/1000.0+1e-9 {
                return Err(format!("invalid compact physical sample at {:?}",rect.id(i)));
            }
            if out.fixed[i]!=crate::t1b::FREE && out.proxy_depth[i]!=0.0 {return Err("mobile proxy on immutable boundary".into());}
        }
        Ok(out)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UnitOwnership {
    pub cell:u32, pub index:u16, pub source_csr:u32, pub record_id:String,
    pub native_parent_budget_m3:f64, pub owned_support_m2:f64,
    pub finite_occupancy_m3:f64, pub requested_mobile_m3:f64,
    pub unallocated_parent_m3:f64,
}
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct OwnershipReport {
    pub rule:String, pub envelope_digest:String, pub tile_digests:Vec<String>,
    pub posts:usize, pub land_posts:usize, pub proxy_positive_posts:usize,
    pub source_mobile_posts:usize, pub proxy_positive_on_excess_posts:usize,
    pub total_solid_m3:f64, pub finite_occupancy_m3:f64,
    pub zero_stock_construction_excess_m3:f64, pub requested_proxy_m3:f64,
    pub source_mobile_m3:f64, pub solid_closure_residual_m3:f64,
    pub quadrature_and_roundoff_bound_m3:f64, pub units:Vec<UnitOwnership>,
    pub native_context_bound:bool,
}
/// A prepared column's mobile layers are ordered top to bottom. Each tuple
/// names the native E1 CSR record and actual available volume, not a proxy ID.
pub struct OwnedTile {
    pub input:CompactInput, pub support_m2:Vec<f64>,
    pub mobile_columns:Vec<Vec<(u32,f64)>>, pub excess_depth_m:Vec<f64>,
    pub report:OwnershipReport,
}
/// `audit_legacy` permits measurement of old proxy evidence against E1. Such a
/// result explicitly cannot enter evolution or claim an E3-native binding.
pub fn bind(input:CompactInput,cfg:&EnvelopeConfig,audit_legacy:bool)->R<OwnedTile>{
    let bound=input.header["source_envelope"]["digest"]==cfg.digest;
    if !bound && !audit_legacy {return Err("native material lacks matching E3 source-envelope binding; regenerate with env and conf=terrain".into());}
    let env=cfg.load()?;
    let decl:vz_tilekey::chart::ChartDeclV1=serde_json::from_value(input.header["chart"].clone()).map_err(|e|e.to_string())?;
    let chart=vz_tilekey::chart::Chart::new(decl).map_err(|e|e.to_string())?;
    let n=input.height_mm.len();let mut support=vec![0.0;n];let mut mobile=vec![Vec::new();n];let mut excess=vec![0.0;n];
    let mut report=OwnershipReport{rule:RULE.into(),envelope_digest:cfg.digest.clone(),tile_digests:vec![input.digest.clone()],posts:n,native_context_bound:bound,..Default::default()};
    let mut units:BTreeMap<(u32,u16),UnitOwnership>=BTreeMap::new();
    let half=vz_tilekey::CHART_HALF_EXTENT_MM/32000;
    for i in 0..n {
        let p=input.rect.id(i);let z=input.height_mm[i] as f64/1000.0;
        // The finite chart's end posts own half a dual square, as in E3.
        let (xa,xb)=((p.gx as f64-0.5).max(-half as f64),(p.gx as f64+0.5).min(half as f64));
        let (ya,yb)=((p.gy as f64-0.5).max(-half as f64),(p.gy as f64+0.5).min(half as f64));
        let (area,error)=crate::t1b_checks::envelope_area((xa+xb)*16.0,(ya+yb)*16.0,(xb-xa)*16.0,(yb-ya)*16.0,chart.radius_mm()/1000.0);
        support[i]=area;
        let positive=input.proxy_depth[i]>0.0; if positive {report.proxy_positive_posts+=1;report.requested_proxy_m3+=area*input.proxy_depth[i];}
        let cell=input.cells[i];
        let Some(e)=env.cells.get(&cell) else {
            if positive{return Err(format!("positive requested stock has no E1 land-cell source: {cell}"));}continue;
        };
        report.land_posts+=1;
        let (lat,lon)=chart.inverse_mm(p.gx*32000,p.gy*32000).map_err(|e|e.to_string())?;
        let tilt=e.tilt(cell,lat.to_degrees(),lon.to_degrees());
        let height=(z-e.sea_m).max(0.0);let over=(z-(e.mean_m+tilt).max(e.sea_m)).max(0.0);excess[i]=over;
        report.total_solid_m3+=area*height;report.zero_stock_construction_excess_m3+=area*over;
        report.quadrature_and_roundoff_bound_m3+=error*height;
        if positive && over>0.0 {report.proxy_positive_on_excess_posts+=1;}
        let mut accounted=over;
        for l in &e.layers {
            let h=e.occupancy(l,tilt,z);accounted+=h;
            let u=units.entry((cell,l.index)).or_insert_with(||UnitOwnership{cell,index:l.index,source_csr:l.csr,record_id:l.id.clone(),native_parent_budget_m3:l.h*e.area_m2,owned_support_m2:0.0,finite_occupancy_m3:0.0,requested_mobile_m3:0.0,unallocated_parent_m3:0.0});
            u.owned_support_m2+=area;u.finite_occupancy_m3+=area*h;report.finite_occupancy_m3+=area*h;
            // Do not tunnel through unowned construction excess to erase a
            // legitimate source underneath. Preserve and report the blocker.
            if positive && over==0.0 {
                let bottom=(l.base+tilt).max(e.sea_m).max(z-input.proxy_depth[i]).max(input.proxy_floor[i]);
                let depth=(z.min(l.base+tilt+l.h)-bottom).max(0.0);
                if depth>0.0 {let v=area*depth;mobile[i].push((l.csr,v));u.requested_mobile_m3+=v;report.source_mobile_m3+=v;}
            }
        }
        if (height-accounted).abs()>1e-7*height.max(1.0) {return Err(format!("E3 source-solid gap at {:?}: {} m",p,height-accounted));}
        if !mobile[i].is_empty(){report.source_mobile_posts+=1;}
    }
    report.solid_closure_residual_m3=(report.total_solid_m3-report.finite_occupancy_m3-report.zero_stock_construction_excess_m3).abs();
    report.quadrature_and_roundoff_bound_m3+=64.0*f64::EPSILON*n as f64*report.total_solid_m3.abs();
    for u in units.values_mut(){
        if u.finite_occupancy_m3>u.native_parent_budget_m3+1e-9*u.native_parent_budget_m3.max(1.0){return Err(format!("native parent budget exceeded {}",u.record_id));}
        u.unallocated_parent_m3=u.native_parent_budget_m3-u.finite_occupancy_m3;
    }
    report.units=units.into_values().collect();
    Ok(OwnedTile{input,support_m2:support,mobile_columns:mobile,excess_depth_m:excess,report})
}

/// Combine reports only after disjoint spatial ownership and common native
/// context are proved. A duplicate tile is a duplicate source, not a cache hit.
pub fn partition(tiles:&[OwnedTile])->R<OwnershipReport>{
    let first=tiles.first().ok_or("empty native material partition")?;
    let mut rects=BTreeSet::new();let mut units:BTreeMap<(u32,u16),UnitOwnership>=BTreeMap::new();
    let mut out=OwnershipReport{rule:RULE.into(),envelope_digest:first.report.envelope_digest.clone(),native_context_bound:true,..Default::default()};
    for t in tiles {
        if !rects.insert((t.input.rect.gy0,t.input.rect.gx0)){return Err("duplicate source owner tile".into());}
        if t.report.envelope_digest!=out.envelope_digest || t.input.header["native_context_pin"]!=first.input.header["native_context_pin"] {return Err("mixed native material context/envelope".into());}
        out.native_context_bound &= t.report.native_context_bound;
        out.tile_digests.extend(t.report.tile_digests.iter().cloned());
        out.posts+=t.report.posts;out.land_posts+=t.report.land_posts;out.proxy_positive_posts+=t.report.proxy_positive_posts;out.source_mobile_posts+=t.report.source_mobile_posts;out.proxy_positive_on_excess_posts+=t.report.proxy_positive_on_excess_posts;
        out.total_solid_m3+=t.report.total_solid_m3;out.finite_occupancy_m3+=t.report.finite_occupancy_m3;out.zero_stock_construction_excess_m3+=t.report.zero_stock_construction_excess_m3;out.requested_proxy_m3+=t.report.requested_proxy_m3;out.source_mobile_m3+=t.report.source_mobile_m3;out.quadrature_and_roundoff_bound_m3+=t.report.quadrature_and_roundoff_bound_m3;
        for u in &t.report.units {
            if let Some(v)=units.get_mut(&(u.cell,u.index)) {
                if v.record_id!=u.record_id || v.source_csr!=u.source_csr || v.native_parent_budget_m3!=u.native_parent_budget_m3{return Err("source record conflict".into());}
                v.owned_support_m2+=u.owned_support_m2;v.finite_occupancy_m3+=u.finite_occupancy_m3;v.requested_mobile_m3+=u.requested_mobile_m3;
            }else{units.insert((u.cell,u.index),u.clone());}
        }
    }
    for u in units.values_mut(){if u.finite_occupancy_m3>u.native_parent_budget_m3+1e-9*u.native_parent_budget_m3.max(1.0){return Err(format!("partition overdraws native source {}",u.record_id));}u.unallocated_parent_m3=u.native_parent_budget_m3-u.finite_occupancy_m3;}
    out.solid_closure_residual_m3=(out.total_solid_m3-out.finite_occupancy_m3-out.zero_stock_construction_excess_m3).abs();out.units=units.into_values().collect();out.tile_digests.sort();Ok(out)
}

pub fn require_evolution_stock(t:&OwnedTile)->R<()> {
    if !t.report.native_context_bound {return Err("historical material audit is not an E3-native process input".into());}
    if t.report.proxy_positive_on_excess_posts!=0 {return Err(format!("{} positive proxy posts lie on zero-stock construction excess; source construction must allocate or evolve these explicitly",t.report.proxy_positive_on_excess_posts));}
    if t.report.solid_closure_residual_m3>t.report.quadrature_and_roundoff_bound_m3{return Err("source-solid closure fails".into());}Ok(())
}

// Same-crate VS4 candidate gate; no key-off evaluator or stock is changed.
#[path = "t1b_emplacement.rs"]
pub mod emplacement;
