//! E5B: refine a process-evolved parent without running its construction twice.
//! This is geometric prolongation, not additional erosion or new material stock.
use crate::{input::T1bWindow, tile::{Ctx,TileOut}, *};
use serde::Deserialize;
use std::collections::{BTreeMap};

#[derive(Deserialize)]
pub struct Evolved {
    pub schema:String,
    pub rule:String,
    pub envelope_digest:String,
    pub columns:Vec<Column>,
    pub chart:Chart,
    pub classes:Vec<pf_region::t1b_ct::CtClassV1>,
}
#[derive(Deserialize)]
pub struct Chart { origin_lat_e7:i64, origin_lon_e7:i64, datum_radius_mm:i64 }
#[derive(Deserialize)]
pub struct Column { pub(crate) cell:u32, shape:[u32;2], radius_m:f64, mean_m:f64, grad:[f64;2], pub(crate) layers:Vec<Layer> }
#[derive(Deserialize)]
pub struct Layer { pub(crate) index:u16, pub(crate) code:u8, pub(crate) depth:f64, pub(crate) h:f64 }
impl Evolved {
    pub fn validate(&self)->R<()> {
        if (self.schema!="e5b-evolved-surface/v1" && self.schema!=pf_region::t1b_escarp::EVOLVED_SCHEMA) || self.columns.is_empty() || self.chart.datum_radius_mm<=0 { return Err("invalid E5B evolved surface".into()); }
        let shape=self.columns[0].shape;
        if shape.contains(&0) || self.columns.iter().any(|c|c.shape!=shape || c.layers.is_empty() || c.layers.iter().any(|l|l.index>255)) { return Err("E5B invalid column or unsupported >255 layer identity".into()); }
        Ok(())
    }
    pub(crate) fn column(&self,x:f64,y:f64)->(&Column,f64) {
        let r=self.chart.datum_radius_mm as f64/1000.;
        let (a,b)=(x/(2.*r),y/(2.*r)); let q=a*a+b*b;
        let (ve,vn,vu)=(2.*a/(1.+q),2.*b/(1.+q),(1.-q)/(1.+q));
        let la0=(self.chart.origin_lat_e7 as f64/1e7).to_radians();
        let lo0=(self.chart.origin_lon_e7 as f64/1e7).to_radians();
        let px=la0.cos()*vu-la0.sin()*vn;let pz=la0.sin()*vu+la0.cos()*vn;
        let (la,lo)=(pz.atan2(px.hypot(ve)).to_degrees(),(lo0+ve.atan2(px)).to_degrees());
        let shape=self.columns[0].shape;
        let row=(((la+90.)/180.*shape[0] as f64).floor() as i64).clamp(0,shape[0] as i64-1);
        let col=(((lo+180.).rem_euclid(360.)/360.*shape[1] as f64).floor() as i64).rem_euclid(shape[1] as i64);
        let cell=(row*shape[1] as i64+col) as u32;
        let c=self.columns.iter().find(|c|c.cell==cell).expect("E5B source column missing, refusing invented rock");
        let cla=-90.+(row as f64+0.5)*180./shape[0] as f64;
        let clo=-180.+(col as f64+0.5)*360./shape[1] as f64;
        let east=c.radius_m*cla.to_radians().cos()*((lo-clo+540.).rem_euclid(360.)-180.).to_radians();
        let north=c.radius_m*(la-cla).to_radians();
        (c,c.mean_m+c.grad[0]*east+c.grad[1]*north)
    }
}
/// Bilinear prolongation on the 2 m integer post lattice. Weights are exact
/// sixteenths of the 32 m parent cell. Round the integer numerator once, ties
/// away from zero, like to_mm, without a metres/mm floating-point tie reversal.
/// This is convex, shared-edge exact after quantization, and cannot overshoot.
pub fn height_mm(w:&T1bWindow,px:i64,py:i64)->i32 {
    let (gx,gy)=(px.div_euclid(16),py.div_euclid(16));
    let (u,v)=(px.rem_euclid(16),py.rem_euclid(16));
    let h=|dx,dy|(w.height_m[w.idx(gx+dx,gy+dy).expect("E5B refinement outside parent halo")]*1000.).round() as i64;
    let numerator=(h(0,0)*(16-u)+h(1,0)*u)*(16-v)+(h(0,1)*(16-u)+h(1,1)*u)*v;
    (if numerator>=0 {(numerator+128)/256} else {-((-numerator+128)/256)}) as i32
}
fn refined_height_mm(ctx:&Ctx,px:i64,py:i64)->i32 {
    if let Some(fine)=&ctx.escarp2 {
        return match &ctx.base_c1 {Some(c1)=>fine.height_with_floor_mm(px,py,c1),None=>fine.height_mm(px,py)};
    }
    if let Some(fine)=&ctx.w.escarpment {return fine.height_mm(px,py);}
    match &ctx.base_c1 { Some(c1)=>c1.height_mm(px,py), None=>height_mm(ctx.w,px,py) }
}
pub fn assemble(ctx:&Ctx,ix:i64,iy:i64)->TileOut {
    let t=std::time::Instant::now();let e=ctx.w.evolved.as_ref().expect("E5B parent");
    let np=TILE_POSTS*TILE_POSTS;let (ox,oy)=crate::tile::tile_origin_post(ix,iy);
    let mut hm=Vec::with_capacity(np);let mut litho=Vec::with_capacity(np);let mut layer=Vec::with_capacity(np);
    let mut fixed=Vec::with_capacity(np);let mut regime=Vec::with_capacity(np);let mut ids=Vec::with_capacity(np);
    let mut diu=Vec::with_capacity(np);let mut hab=Vec::with_capacity(np);let mut used=BTreeMap::new();
    for j in 0..TILE_POSTS {for i in 0..TILE_POSTS {
        let (x,y)=((ox+i as i64) as f64*POST_M,(oy+j as i64) as f64*POST_M);
        let h=refined_height_mm(ctx,ox+i as i64,oy+j as i64);let z=h as f64/1000.;hm.push(h);
        let (c,datum)=e.column(x,y);let depth=datum-z;
        let l=c.layers.iter().find(|l|depth<l.depth).expect("E5B surface below retained envelope, refusing invented rock");
        litho.push(l.code);layer.push(l.index as u8);
        let parent=ctx.w.nearest(x,y).expect("E5B parent support");
        let f=match ctx.w.fixed[parent] {0=>0,1|5|7=>1,2|4=>3,_=>2};fixed.push(f);
        let mut slope=0.0_f64;
        for dy in -1i32..=1 {for dx in -1i32..=1 {if dx==0 && dy==0 {continue;}
            let dz=(h-refined_height_mm(ctx,ox+i as i64+dx as i64,oy+j as i64+dy as i64)) as f64/1000.;
            slope=slope.max(dz.abs()/(POST_M*(dx as f64).hypot(dy as f64)));
        }}
        let k=e.classes.iter().find(|k|k.code==l.code && k.mode!=pf_region::t1b_ct::CT_NONE);
        regime.push(if f!=0 {0} else if let Some(k)=k {if slope>=k.cliff_tan_permille as f64/1000. {2} else {3}} else {0});
        // These contacts and indices are actual retained envelope layers, not
        // procedural bedding. No minimum-thickness filtering or renumbering.
        diu.push((depth-(l.depth-l.h)) as f32);hab.push((l.depth-depth) as f32);
        let id=crate::strata_prov::strata_id(c.cell,l.index);ids.push(id);used.insert(id,(c,l));
    }}
    let rows=if let Some(pv)=&ctx.prov {used.iter().map(|(_id,(c,l))|{
        let a=ctx.t0.off[c.cell as usize] as usize;let b=ctx.t0.off[c.cell as usize+1] as usize;
        assert_eq!(ctx.t0.lith[a+l.index as usize],l.code,"E5B source lithology mismatch");
        pv.strata.strata_row(c.cell,l.index,l.code,l.h as f32,b-a)
    }).collect()}else{Vec::new()};
    if ctx.prov.is_none(){ids.clear();}
    let raw_mean=(0..TILE_POSTS).flat_map(|j|(0..TILE_POSTS).map(move |i|ctx.w.bicubic((ox+i as i64)as f64*POST_M,(oy+j as i64)as f64*POST_M))).sum::<f64>()/np as f64;
    let mean=hm.iter().map(|&h|h as f64/1000.).sum::<f64>()/np as f64;
    let counts:Vec<_>=(0..4).map(|r|regime.iter().filter(|&&v|v==r).count()).collect();
    let mut out = TileOut {ix,iy,height_mm:hm.clone(),base_mm:hm,litho,layer,fixed,closure_mm:vec![0;np],cut_minus_base_m:vec![0.;np],
        base_mean_m:mean,raw_t1b_mean_m:raw_mean,mean_shift_c_m:0.,mean_deficit_m:0.,shift_max_m:0.,lip_posts:0,lip_raise_mm:0,
        fill_raised_posts:0,fill_max_raise_mm:0,relax_seconds_own:0.,assemble_seconds:t.elapsed().as_secs_f64(),
        bed:Vec::new(),over:Vec::new(),remap_mm:Vec::new(),prov_id:ids,prov_rows:rows,depth_in_unit_m:diu,height_above_unit_base_m:hab,
        bed_style_id:Vec::new(),flow_unit:Vec::new(),pyro_sub:Vec::new(),pyro_stats:None,de_mm:Vec::new(),de_talus:Vec::new(),de_stats:None,
        sharp_mm:vec![0;np],ct_regime:regime,sharp_stats:Some(serde_json::json!({"stage":"evolved_parent_refinement","rule":e.rule,
            "envelope_digest":e.envelope_digest,"regime_posts":counts,"new_simulated_time_yr":0,"new_material_volume_m3":0,
            "base_mm_semantics":"evolved parent bilinear reference, supersedes the legacy generic B1 bicubic field note",
            "raw_bicubic_mean_m":raw_mean,"bilinear_minus_bicubic_mean_m":mean-raw_mean,
            "regime_semantics":"0 non-cliff soft rock or fixed domain; 2 resolved bare cliff; 3 exposed non-cliff hard rock; no talus deposit is implied",
            "limitation":"2m sampling of the evolved32m surface; no additional fine erosion, construction, mantle or procedural bedding"})),
        levee_mean_m:None,carve_mean_m:None,b8_cap_infeasible:None};
    if let Some(fine)=&ctx.w.escarpment {
        let fronts=crate::escarpment::fronts(e,&out.height_mm,ox,oy);
        let stats=out.sharp_stats.as_mut().unwrap();
        stats["base_mm_semantics"]="pre-E5 bilinear parent plus already ledgered 2m removal; final coarse heights are restricted means".into();
        stats["escarpment"]=fine.meta.clone();
        stats["escarpment_fronts"]=serde_json::to_value(fronts).expect("finite contact geometry");
        stats["front_semantics"]="Triangle-mesh intersections with actual E3 cap basal surfaces. Basal endpoints are sub-post roots; rim/toe endpoints are the actual high/low corners of the same collision triangles. Equal endpoints denote a triangular fan apex. Adjacent elements weld by global coordinates. Above-recorded-top corners are explicitly flagged, not invented rock.".into();
        stats["limitation"]="2m four-face E5 process; no joints, talus, overhangs or convergence claim; 32m pinning masks remain".into();
        stats.as_object_mut().unwrap().remove("bilinear_minus_bicubic_mean_m");
    }
    if let Some(fine)=&ctx.escarp2 {
        let stats=out.sharp_stats.as_mut().unwrap();
        stats["rule"]=fine.rule().into();
        stats["base_mm_semantics"]=if fine.is_antialiased(){"ESCARP-4 anti-aliased band-conservative 2m reconstruction; per-post means may differ"}else if fine.is_band(){"ESCARP-3 band-conservative 2m reconstruction; per-post means may differ"}else{"ESCARP-2 conservative 2m reconstruction of calibrated 32m post-prism means"}.into();
        stats["escarpment"]=fine.ledger_with_floor(ox,oy,ctx.base_c1.as_ref());
        stats["escarpment_fronts"]=serde_json::to_value(crate::escarpment::fronts(e,&out.height_mm,ox,oy)).expect("finite contact geometry");
        stats["limitation"]=if fine.is_band(){"Window-scoped source-band conservation. Individual posts/tiles can have signed residuals. Across-parent-window agreement is unverified. Outside the contact band uses uncorrected C1 cubic interpolation. Post-prism volume is not triangle-integrated volume. No new erosion time or collectible stock."}else{"Geometric reconstruction only; no new erosion time, talus, joints or collectible stock. Exact post-prism means do not imply exact triangle-integrated volume. Coarse constraints can introduce local extrema; measured acceptance is required."}.into();
        stats.as_object_mut().unwrap().remove("bilinear_minus_bicubic_mean_m");
    }
    if ctx.base_c1.is_some() {
        let stats=out.sharp_stats.as_mut().unwrap();
        stats["base_mm_semantics"]="evolved parent C1 bounded bicubic reference; supersedes legacy generic B1 field note".into();
        stats["base_interp_rule"]=crate::base_interp::RULE.into();
        stats.as_object_mut().unwrap().remove("bilinear_minus_bicubic_mean_m");
        stats["c1_minus_legacy_bicubic_mean_m"]=(mean-raw_mean).into();
        if ctx.escarp2.is_some() {
            stats["base_mm_semantics"]="ESCARP-4 exact conserved band; bounded C1 floor outside with 64m C1 exterior partition".into();
            stats["base_interp_scope"]="base_interp=c1 with ct_escarp=4 only: all existing band prisms unchanged; outside blend fixed parent heights remain exact at posts, not flattened over 32m prisms. Monotone edges and boundedness apply to the pure C1 area, not the ESCARP band or blend.".into();
            stats["escarpment"]["support"]="owner centre in tile half-open core; full nearest-post support [-8,7]; tile accounting uses final composed heights; unchanged window certificates cover the ESCARP band only".into();
            stats["limitation"]="Outside-band C1 reconstruction is geometry, with no per-post volume constraint or new material stock. The blend preserves every band sample and its exact volume. Across-parent-window agreement remains unverified. Integer-mm output approximates continuous C1; no new erosion time.".into();
        }
    }
    if let Some(toe)=ctx.escarp2.as_ref().and_then(|f|f.toe_ledger()) {
        let stats=out.sharp_stats.as_mut().unwrap();
        stats["escarpment_toes"]=toe;
        stats["escarpment"]["window_band_accounts_semantics"]="ESCARP4 baseline certificate, plus the separately published zero-net toe increment. Original gross debit/credit describes the baseline distribution; tile debit/credit uses final heights.".into();
        stats["escarpment"]["support"]="Owner centre in half-open tile core, full nearest-post support [-8,7]; final composed tile heights. The complete window has the original source/cap net volume plus exact zero-net toe increments.".into();
        stats["base_mm_semantics"]="ESCARP4 plus compact C1 small-contact-step reconstruction, exact incremental source/cap volume balance, bounded C1 outside band".into();
        stats["base_interp_scope"]="ct_escarp=6 requires base_interp=c1; high-cliff samples unchanged; low-relief toe prisms modified with a separate exact incremental ledger".into();
        stats["limitation"]="Conditional geometric reconstruction, not a talus or erosion process. Window-scoped post-prism balance is not triangle-integrated volume. No new material stock; exposed identity is looked up at the final height in unchanged native columns.".into();
    }
    out
}
