//! Channel registration owned by LS on an audited supplied surface graph.
//! Classification is separate from physical protection and never changes routing.
//! Every positive-weight downstream branch is retained through a native boundary.
use crate::t1b_surface_domain::{Domain,PostId,Boundary,DomainError,Result,invalid,content_digest};
use crate::t1b_ls::DetailRoute;
use serde::{Serialize,Deserialize};
use std::collections::{BTreeMap,BTreeSet};
pub const RULE:&str="ls_channel_graph_v1: native reporting discharge from once-quantized effective area; register every above-threshold post and every downstream positive-weight edge to a declared native boundary; stable point/reach identity from generation, chart and global post; no change to physical protection, source flux or terrain; audit against pre-registration geometry; legacy single-attachment projection rejects ambiguous attachment";
#[derive(Clone,Debug,Serialize,Deserialize,PartialEq)]
pub struct Edge{pub receiver:PostId,pub weight:f64,pub effective_flux_m2:f64}
#[derive(Clone,Debug,Serialize,Deserialize,PartialEq)]
pub struct ChannelPost{
    pub id:PostId,pub reach_id:String,pub above_threshold:bool,
    pub ground_m:f64,pub routing_m:f64,pub effective_area_m2:f64,pub published_area_m2:f32,pub discharge_m3_s:f64,
    pub edges:Vec<Edge>,pub boundary:Option<Boundary>,pub native_attachments:Vec<u32>,
}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Registry{
    pub rule:String,pub generation:String,pub chart:String,pub domain_digest:String,pub route_digest:String,
    pub runoff_reference_m_s:f64,pub threshold_m3_s:f64,pub posts:Vec<ChannelPost>,pub digest:String,
}
#[derive(Debug,Serialize,Deserialize)]
pub struct Audit{
    pub domain_posts:usize,pub registered_posts:usize,pub supra_threshold_posts:usize,pub connected_edges:usize,
    pub boundary_posts:usize,pub ambiguous_legacy_attachments:usize,
    pub geometric_source_m2:f64,pub geometric_export_m2:f64,pub effective_source_m2:f64,pub effective_export_m2:f64,
    pub max_effective_field_error_m2:f64,pub pass:bool,
}
/// Rebuild source accumulation and traverse actual edges. Do not trust a mask,
/// terminal label, certificate flag, serialized order, or the route's own ledger.
pub fn audit_route(domain:&Domain,route:&DetailRoute)->Result<Audit>{
    let p=domain.posts();let n=p.len();
    if route.domain_digest!=domain.digest(){return Err(DomainError::Stale{field:"channel surface/domain".into()});}
    if [&route.routing_m,&route.geometric_area_m2,&route.effective_area_m2].iter().any(|v|v.len()!=n)||route.edges.len()!=n||route.rank.len()!=n{return Err(invalid("channel route dimensions"));}
    let mut copy=route.clone();copy.digest.clear();if content_digest(&copy)?!=route.digest{return Err(DomainError::Stale{field:"channel route content".into()});}
    let mut geometric=vec![0.0;n];let mut effective=vec![0.0;n];let mut source_g=0.0;let mut source_q=0.0;let mut edge_count=0;let mut capacity=0.0f64;
    for i in 0..n{
        if !route.routing_m[i].is_finite()||route.routing_m[i]<p[i].ground_m{return Err(invalid("invalid routing potential"));}
        if p[i].boundary.is_some(){if !route.edges[i].is_empty()||route.routing_m[i]!=p[i].ground_m{return Err(invalid("changed physical boundary"));}continue;}
        capacity+=(route.routing_m[i]-p[i].ground_m)*p[i].support_m2;
        geometric[i]=p[i].support_m2;effective[i]=p[i].source_m2+p[i].external_m2;source_g+=geometric[i];source_q+=effective[i];
        let mut sum=0.0;let mut seen=BTreeSet::new();
        for e in &route.edges[i]{
            let j=e.receiver as usize;if j>=n||!seen.insert(j)||!e.weight.is_finite()||e.weight<=0.0{return Err(invalid("invalid or duplicate channel edge"));}
            if !(0..8).any(|k|domain.edge(i,k)==Some(j)){return Err(invalid("channel edge violates native legal/protection policy"));}
            if !(route.routing_m[i]>route.routing_m[j]||(route.routing_m[i]==route.routing_m[j]&&route.rank[i]>route.rank[j])){return Err(invalid("uphill or cyclic channel potential"));}
            sum+=e.weight;edge_count+=1;
        }
        if (sum-1.0).abs()>1e-12{return Err(invalid("missing channel edge or nonconservative partition"));}
    }
    if !capacity.is_finite()||capacity!=route.fill_capacity_m3{return Err(invalid("route fill capacity disagrees with supplied ground/support"));}
    let active=p.iter().filter(|p|p.boundary.is_none()).count();
    if route.downstream_first.len()!=active{return Err(invalid("route order loses active sources"));}
    let mut positions=vec![usize::MAX;n];
    for (k,&iu) in route.downstream_first.iter().enumerate(){let i=iu as usize;if i>=n||p[i].boundary.is_some()||positions[i]!=usize::MAX{return Err(invalid("invalid or duplicated route order"));}positions[i]=k;}
    for i in 0..n{for e in &route.edges[i]{let j=e.receiver as usize;if p[j].boundary.is_none()&&positions[j]>=positions[i]{return Err(invalid("route order is not downstream-first"));}}}
    let mut order:Vec<usize>=(0..n).collect();order.sort_by(|&i,&j|route.routing_m[i].total_cmp(&route.routing_m[j]).then(route.rank[i].cmp(&route.rank[j])).then(p[i].id.cmp(&p[j].id)));
    for &i in order.iter().rev(){for e in &route.edges[i]{let j=e.receiver as usize;geometric[j]+=geometric[i]*e.weight;effective[j]+=effective[i]*e.weight;}}
    let mut export_g=0.0;let mut export_q=0.0;let mut err=0f64;
    for i in 0..n{
        if !effective[i].is_finite()||!geometric[i].is_finite(){return Err(invalid("nonfinite channel accumulation"));}
        err=err.max((effective[i]-route.effective_area_m2[i]).abs());
        if geometric[i]!=route.geometric_area_m2[i]||effective[i]!=route.effective_area_m2[i]{return Err(invalid(format!("published area does not equal actual source accumulation at {:?}",p[i].id)));}
        if p[i].boundary.is_some(){export_g+=geometric[i];export_q+=effective[i];}
    }
    if (source_g-export_g).abs()>1e-9*source_g.max(1.0)||(source_q-export_q).abs()>1e-9*source_q.max(1.0){return Err(invalid("channel source/export balance"));}
    let mut balances=BTreeMap::<u32,(f64,f64,f64,f64)>::new();
    for (i,post) in p.iter().enumerate() {
        let b=balances.entry(post.basin).or_default();
        if post.boundary.is_none(){b.0+=post.support_m2;b.1+=post.source_m2+post.external_m2;}
        else{b.2+=geometric[i];b.3+=effective[i];}
    }
    if route.balances.len()!=balances.len(){return Err(invalid("missing or extraneous basin source ledger"));}
    for (actual,(&basin,&(gs,qs,ge,qe))) in route.balances.iter().zip(&balances) {
        if actual.basin!=basin||actual.geometric_source_m2!=gs||actual.effective_source_m2!=qs||actual.geometric_export_m2!=ge||actual.effective_export_m2!=qe{return Err(invalid("basin ledger differs from actual native source accumulation"));}
    }
    // Independently reconstruct provenance-label transfers from actual source
    // accumulation. Rehashing a false transfer receipt cannot certify it.
    let mut transfers=BTreeMap::<(u32,u32),(f64,f64)>::new();
    if domain.policy==crate::t1b_surface_domain::Policy::NativePhysicalV1 {
        for i in 0..n {for e in &route.edges[i] {let j=e.receiver as usize;
            if p[i].basin!=p[j].basin {let v=transfers.entry((p[i].basin,p[j].basin)).or_default();v.0+=geometric[i]*e.weight;v.1+=effective[i]*e.weight;}
        }}
    }
    if route.basin_transfers.len()!=transfers.len(){return Err(invalid("missing or extraneous physical basin transfer"));}
    for (actual,(&(from,to),&(g,q))) in route.basin_transfers.iter().zip(&transfers) {
        if actual.from_basin!=from||actual.to_basin!=to||actual.geometric_m2!=g||actual.effective_m2!=q{return Err(invalid("physical basin transfer differs from actual graph/source accumulation"));}
    }
    Ok(Audit{domain_posts:n,registered_posts:0,supra_threshold_posts:0,connected_edges:edge_count,boundary_posts:p.iter().filter(|p|p.boundary.is_some()).count(),ambiguous_legacy_attachments:0,geometric_source_m2:source_g,geometric_export_m2:export_g,effective_source_m2:source_q,effective_export_m2:export_q,max_effective_field_error_m2:err,pass:true})
}
pub fn register(domain:&Domain,route:&DetailRoute,runoff_reference_m_s:f64,threshold_m3_s:f64)->Result<Registry>{
    audit_route(domain,route)?;
    if !matches!(domain.policy,crate::t1b_surface_domain::Policy::NativeDetailV1|crate::t1b_surface_domain::Policy::NativePhysicalV1){return Err(invalid("channel registration requires the native routing policy, not construction minimax"));}
    if !runoff_reference_m_s.is_finite()||runoff_reference_m_s<=0.0||!threshold_m3_s.is_finite()||threshold_m3_s<=0.0{return Err(invalid("positive native runoff reference and literal discharge threshold required"));}
    let p=domain.posts();let n=p.len();let mut member=vec![false;n];let mut above=vec![false;n];
    for i in 0..n{let a=route.effective_area_m2[i]as f32;let q=a as f64*runoff_reference_m_s;if !a.is_finite()||!q.is_finite(){return Err(invalid("published area/discharge overflows representation"));}if p[i].boundary.is_none(){above[i]=q>=threshold_m3_s;member[i]=above[i];}}
    let mut stack:Vec<_>=(0..n).filter(|&i|member[i]).collect();
    while let Some(i)=stack.pop(){for e in &route.edges[i]{let j=e.receiver as usize;if !member[j]{member[j]=true;stack.push(j);}}}
    let mut order:Vec<_>=(0..n).collect();order.sort_by(|&i,&j|route.routing_m[i].total_cmp(&route.routing_m[j]).then(route.rank[i].cmp(&route.rank[j])).then(p[i].id.cmp(&p[j].id)));
    let mut attach=vec![BTreeSet::new();n];
    for &i in &order{if member[i]{if let Some(b)=&p[i].boundary{if let Some(a)=b.attachment{attach[i].insert(a);}}else{for e in &route.edges[i]{let a=attach[e.receiver as usize].clone();attach[i].extend(a);}}}}
    let mut posts=Vec::new();
    for i in 0..n{if member[i]{
        let reach_id=content_digest(&(RULE,&domain.pin().generation,&domain.pin().chart,p[i].id))?;
        posts.push(ChannelPost{id:p[i].id,reach_id,above_threshold:above[i],ground_m:p[i].ground_m,routing_m:route.routing_m[i],effective_area_m2:route.effective_area_m2[i],published_area_m2:route.effective_area_m2[i]as f32,discharge_m3_s:(route.effective_area_m2[i]as f32)as f64*runoff_reference_m_s,edges:route.edges[i].iter().map(|e|Edge{receiver:p[e.receiver as usize].id,weight:e.weight,effective_flux_m2:route.effective_area_m2[i]*e.weight}).collect(),boundary:p[i].boundary.clone(),native_attachments:attach[i].iter().copied().collect()});
    }}
    let mut reg=Registry{rule:RULE.into(),generation:domain.pin().generation.clone(),chart:domain.pin().chart.clone(),domain_digest:domain.digest().into(),route_digest:route.digest.clone(),runoff_reference_m_s,threshold_m3_s,posts,digest:String::new()};reg.digest=content_digest(&reg)?;Ok(reg)
}
/// Complete independent coverage check against pre-registration geometry. Rebuild
/// the required closure from supra-threshold flow, not from the submitted mask.
pub fn audit(domain:&Domain,route:&DetailRoute,reg:&Registry)->Result<Audit>{
    let mut result=audit_route(domain,route)?;
    if reg.domain_digest!=domain.digest()||reg.route_digest!=route.digest||reg.rule!=RULE{return Err(DomainError::Stale{field:"channel registry revision".into()});}
    let mut copy=reg.clone();copy.digest.clear();if content_digest(&copy)?!=reg.digest{return Err(DomainError::Stale{field:"channel registry content".into()});}
    if !matches!(domain.policy,crate::t1b_surface_domain::Policy::NativeDetailV1|crate::t1b_surface_domain::Policy::NativePhysicalV1)||!reg.runoff_reference_m_s.is_finite()||reg.runoff_reference_m_s<=0.0||!reg.threshold_m3_s.is_finite()||reg.threshold_m3_s<=0.0{return Err(invalid("invalid registry policy/runoff/threshold"));}
    // Traverse the authoritative graph independently of register() and of the
    // submitted membership/attachment lists. Every submitted edge is compared
    // directly to the owner; every required node must be reached from a source.
    let p=domain.posts();let mut submitted=BTreeMap::new();
    for row in &reg.posts{if submitted.insert(row.id,row).is_some(){return Err(invalid("duplicate registered post"));}}
    let supra:Vec<usize>=(0..p.len()).filter(|&i|p[i].boundary.is_none()&&(route.effective_area_m2[i]as f32)as f64*reg.runoff_reference_m_s>=reg.threshold_m3_s).collect();
    let mut required=BTreeSet::new();let mut stack=supra.clone();
    while let Some(i)=stack.pop(){if !required.insert(i){continue;}for e in &route.edges[i]{stack.push(e.receiver as usize);}}
    if submitted.len()!=required.len(){return Err(invalid("orphaned or missing registered nodes"));}
    for &i in &required{
        let row=submitted.get(&p[i].id).ok_or_else(||invalid("supra-threshold flow or downstream support is unregistered"))?;
        let area=route.effective_area_m2[i]as f32;let q=area as f64*reg.runoff_reference_m_s;
        let above=p[i].boundary.is_none()&&q>=reg.threshold_m3_s;
        if !q.is_finite()||row.above_threshold!=above||row.ground_m!=p[i].ground_m||row.routing_m!=route.routing_m[i]||row.effective_area_m2!=route.effective_area_m2[i]||row.published_area_m2.to_bits()!=area.to_bits()||row.discharge_m3_s!=q||row.boundary!=p[i].boundary||row.edges.len()!=route.edges[i].len(){return Err(invalid("registered geometry, flux or role differs from the owner"));}
        if row.reach_id!=content_digest(&(RULE,&domain.pin().generation,&domain.pin().chart,p[i].id))?{return Err(invalid("unstable registered reach identity"));}
        let mut edges=BTreeMap::new();for e in &row.edges{if edges.insert(e.receiver,e).is_some(){return Err(invalid("duplicated registered edge"));}}
        for e in &route.edges[i]{let id=p[e.receiver as usize].id;let actual=edges.get(&id).ok_or_else(||invalid("deleted registered channel edge"))?;if actual.weight!=e.weight||actual.effective_flux_m2!=route.effective_area_m2[i]*e.weight||!submitted.contains_key(&id){return Err(invalid("registered edge loses flux or downstream connection"));}}
        // A separate endpoint walk verifies the submitted attachment set. This
        // is deliberately independent of the producer's reverse-order union.
        let mut endpoints=BTreeSet::new();let mut visited=BTreeSet::new();let mut todo=vec![i];
        while let Some(j)=todo.pop(){if !visited.insert(j){continue;}if let Some(b)=&p[j].boundary{if let Some(a)=b.attachment{endpoints.insert(a);}}else{for e in &route.edges[j]{todo.push(e.receiver as usize);}}}
        if row.native_attachments!=endpoints.into_iter().collect::<Vec<_>>(){return Err(invalid("false or missing native attachment"));}
    }
    result.registered_posts=reg.posts.len();result.supra_threshold_posts=reg.posts.iter().filter(|p|p.above_threshold).count();result.ambiguous_legacy_attachments=reg.posts.iter().filter(|p|p.boundary.is_none()&&p.native_attachments.len()!=1).count();Ok(result)
}
/// Legacy GUIDE publication is allowed only when it has a unique actual native
/// attachment. Physical fixed/protection fields in the owner domain are untouched.
/// A multi-terminal MFD branch must be represented by a native graph adapter;
/// never choose an arbitrary attachment merely to make the old checker pass.
pub fn guide_projection(domain:&Domain,route:&DetailRoute,reg:&Registry)->Result<BTreeMap<PostId,(u8,u32,f32)>>{
    audit(domain,route,reg)?;let mut out=BTreeMap::new();
    for p in &reg.posts{if p.boundary.is_none(){if p.native_attachments.len()!=1{return Err(invalid(format!("native GUIDE projection cannot represent {} attachments at {:?}",p.native_attachments.len(),p.id)));}out.insert(p.id,(crate::t1b::GUIDE,p.native_attachments[0],p.published_area_m2));}}
    Ok(out)
}

/// A graph product carries every branch and its native terminal. It is not a
/// GUIDE raster and does not change the physical fixed/protection field.
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct NativeGraphProduct {
    pub schema:String,pub domain_digest:String,pub route_digest:String,pub registry:Registry,
    pub local_geometric_source_m2:f64,pub local_effective_source_m2:f64,
    pub ports:Vec<NativePort>,pub coarse_network_composed:bool,pub digest:String,
}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct NativePort {
    pub boundary:Boundary,pub geometric_export_m2:f64,pub effective_export_m2:f64,
    pub native_posts:Vec<PostId>,
}
/// Compose once-owned native sources into actual boundary ports. Coarse fixed
/// channel discharge is deliberately NOT added: it contains overlapping source
/// area. Replacing coarse contributions requires an explicit spatial partition
/// and remains a separate acceptance gate, recorded in the product.
pub fn native_graph_product(domain:&Domain,route:&DetailRoute,reg:&Registry)->Result<NativeGraphProduct>{
    let audit=audit(domain,route,reg)?;
    let mut ports=BTreeMap::<String,NativePort>::new();
    for(i,p)in domain.posts().iter().enumerate(){if let Some(b)=&p.boundary{
        let port=ports.entry(b.id.clone()).or_insert_with(||NativePort{boundary:b.clone(),geometric_export_m2:0.0,effective_export_m2:0.0,native_posts:Vec::new()});
        if port.boundary!=*b{return Err(invalid("native port identity conflict"));}
        port.geometric_export_m2+=route.geometric_area_m2[i];port.effective_export_m2+=route.effective_area_m2[i];port.native_posts.push(p.id);
    }}
    let mut out=NativeGraphProduct{schema:"vz.native-detail-graph/v1".into(),domain_digest:domain.digest().into(),route_digest:route.digest.clone(),registry:reg.clone(),local_geometric_source_m2:audit.geometric_source_m2,local_effective_source_m2:audit.effective_source_m2,ports:ports.into_values().collect(),coarse_network_composed:false,digest:String::new()};
    out.digest=content_digest(&out)?;Ok(out)
}
pub fn audit_native_graph_product(domain:&Domain,route:&DetailRoute,reg:&Registry,product:&NativeGraphProduct)->Result<()>{
    let expected=native_graph_product(domain,route,reg)?;
    if content_digest(product)?!=content_digest(&expected)?{return Err(invalid("native graph product loses a branch, source, terminal or ownership boundary"));}Ok(())
}

/// Own-source replacement against the ACTUAL published T1a source quadrature.
/// A fine post removes its footprint's coarse own source before its routed
/// export is injected at the verified native attachment. Accumulated discharge
/// is a baseline output, never a second local source.
pub struct CoarseWater {
    pub skeleton_digest:String,pub chart:vz_tilekey::chart::ChartDeclV1,
    own_effective_m2:Vec<f64>,baseline_effective_m2:Vec<f64>,receiver:Vec<u32>,
}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct SourceReplacement {
    pub coarse_node:u32,pub replaced_fraction:f64,pub removed_coarse_own_m2:f64,
}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct WaterPartition {
    pub coarse_digest:String,pub fine_posts:usize,pub fine_own_m2:f64,
    pub removed_coarse_own_m2:f64,pub refinement_delta_m2:f64,
    pub sources:Vec<SourceReplacement>,pub digest:String,
}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct WaterComposition {
    pub partition:WaterPartition,pub graph_digest:String,pub coarse_effective_m2:Vec<f64>,
    pub direct_native_water_exports_m2:f64,pub coarse_terminal_delta_m2:f64,
    pub source_export_residual_m2:f64,pub coarse_network_composed:bool,pub digest:String,
}
impl CoarseWater {
    pub fn load(path:&std::path::Path,runoff_reference_m_s:f64)->std::result::Result<Self,String>{
        if !runoff_reference_m_s.is_finite()||runoff_reference_m_s<=0.0{return Err("coarse runoff reference must be positive".into());}
        let view=crate::t1b_view::SkeletonView::load(path)?;
        let j:serde_json::Value=serde_json::from_slice(&std::fs::read(path.join("skeleton.json")).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
        let h=&j["header"];
        // SkeletonView verifies all source fields, including these two. Verify
        // these re-read bytes as well so a concurrent mutation is not admitted.
        let read=|name:&str|->std::result::Result<Vec<f64>,String>{
            let f=h["fields"].as_array().ok_or("coarse fields")?.iter().find(|f|f["name"]==name).ok_or_else(||format!("coarse source field {name} absent"))?;
            let file=f["file"].as_str().ok_or("coarse source locator")?;
            let b=std::fs::read(path.join(file)).map_err(|e|e.to_string())?;
            if f["digest"]!=format!("blake3:{}",blake3::hash(&b).to_hex())||b.len()!=crate::NODES*4{return Err("coarse source field digest/dimensions".into());}
            Ok(b.chunks_exact(4).map(|b|f32::from_le_bytes(b.try_into().unwrap()) as f64).collect())
        };
        let area=read("area_m2")?;
        let runoff=if view.climate_digest.is_some(){read("runoff_m_per_yr")?}else{vec![runoff_reference_m_s*crate::SECONDS_PER_YEAR;crate::NODES]};
        let mut own=vec![0.0;crate::NODES];
        for i in 0..crate::NODES{if [0,1,4].contains(&view.kind[i]){own[i]=area[i]*runoff[i]/(runoff_reference_m_s*crate::SECONDS_PER_YEAR);}
            if !own[i].is_finite()||own[i]<0.0||!view.discharge_m3_s[i].is_finite()||view.discharge_m3_s[i]<0.0{return Err("invalid actual coarse water source".into());}}
        Ok(Self{skeleton_digest:view.digest,chart:serde_json::from_value(h["recipe"]["chart"].clone()).map_err(|e|e.to_string())?,own_effective_m2:own,baseline_effective_m2:view.discharge_m3_s.iter().map(|q|*q as f64/runoff_reference_m_s).collect(),receiver:view.receiver})
    }
    /// Each fine dual square is intersected with coarse dual squares, including
    /// exact half shares at node boundaries. The original coarse source is a
    /// node quadrature, so its density is constant on its own support. The fine
    /// interpolation difference is a declared refinement delta, never fitted.
    fn partition(&self,samples:&[(PostId,f64)])->Result<WaterPartition>{
        let mut fractions=BTreeMap::<u32,f64>::new();let mut seen=BTreeSet::new();let mut fine=0.0;
        for &(p,own) in samples {
            if !seen.insert(p){return Err(DomainError::DuplicateSource{post:p});}
            if !own.is_finite()||own<0.0{return Err(invalid("invalid fine own source"));}fine+=own;
            let h=vz_tilekey::CHART_HALF_EXTENT_MM/32000;
            if p.gx.abs()>h||p.gy.abs()>h{return Err(invalid("fine source outside actual coarse chart"));}
            let (xa,xb,ya,yb)=((p.gx as f64-0.5).max(-h as f64),(p.gx as f64+0.5).min(h as f64),(p.gy as f64-0.5).max(-h as f64),(p.gy as f64+0.5).min(h as f64));
            for gy in ((ya+8.0)/16.0).floor() as i64..=((yb+8.0)/16.0).floor() as i64 {
                for gx in ((xa+8.0)/16.0).floor() as i64..=((xb+8.0)/16.0).floor() as i64 {
                    let x=(xb.min(gx as f64*16.0+8.0)-xa.max(gx as f64*16.0-8.0)).max(0.0);
                    let y=(yb.min(gy as f64*16.0+8.0)-ya.max(gy as f64*16.0-8.0)).max(0.0);
                    if x*y==0.0{continue;}
                    let node=crate::node_index(gx,gy).ok_or_else(||invalid("fine ownership beyond native coarse node support"))?;
                    *fractions.entry(node as u32).or_default()+=x*y/256.0;
                }
            }
        }
        let mut removed=0.0;let mut sources=Vec::new();
        for (node,f) in fractions {if f>1.0{return Err(invalid("fine ownership overlaps more than one coarse source support"));}let q=self.own_effective_m2[node as usize]*f;removed+=q;sources.push(SourceReplacement{coarse_node:node,replaced_fraction:f,removed_coarse_own_m2:q});}
        let mut p=WaterPartition{coarse_digest:self.skeleton_digest.clone(),fine_posts:samples.len(),fine_own_m2:fine,removed_coarse_own_m2:removed,refinement_delta_m2:fine-removed,sources,digest:String::new()};p.digest=content_digest(&p)?;Ok(p)
    }
    pub fn material_partition(&self,tiles:&[crate::t1b_ve::material::CompactInput])->Result<WaterPartition>{
        let mut samples=Vec::new();let mut pin=None;
        for t in tiles{
            if serde_json::to_value(&self.chart).map_err(|e|invalid(e.to_string()))?!=t.header["chart"]{return Err(invalid("coarse/fine chart mismatch"));}
            let tp=t.header["native_context_pin"].as_str().ok_or_else(||invalid("material context"))?;
            if pin.is_some_and(|p|p!=tp){return Err(invalid("mixed fine source contexts"));}pin=Some(tp);
            for i in 0..t.height_mm.len(){if t.fixed[i]==crate::t1b::FREE{samples.push((t.rect.id(i),1024.0*t.runoff[i] as f64));}}
        }
        self.partition(&samples)
    }
    pub fn compose(&self,domain:&Domain,route:&DetailRoute,reg:&Registry)->Result<WaterComposition>{
        if domain.pin().spacing_mm!=32000{return Err(invalid("coarse replacement requires native32m posts"));}
        if domain.pin().chart!=content_digest(&self.chart)? {return Err(invalid("fine domain must bind the actual coarse chart declaration digest"));}
        let graph=native_graph_product(domain,route,reg)?;
        let samples:Vec<_>=domain.posts().iter().filter(|p|p.boundary.is_none()).map(|p|(p.id,p.source_m2)).collect();
        if domain.posts().iter().any(|p|p.external_m2!=0.0){return Err(invalid("external fine inflow needs a separately disjoint source owner"));}
        let partition=self.partition(&samples)?;
        if (partition.fine_own_m2-graph.local_effective_source_m2).abs()>1e-9*partition.fine_own_m2.max(1.0){return Err(invalid("fine graph and source partition disagree"));}
        let mut delta=vec![0.0;crate::NODES];for s in &partition.sources{delta[s.coarse_node as usize]-=s.removed_coarse_own_m2;}
        let mut direct=0.0;
        for p in &graph.ports {match p.boundary.attachment {
            Some(node) if [crate::t1b::CHANNEL,crate::t1b::GUIDE].contains(&p.boundary.kind) => {let d=delta.get_mut(node as usize).ok_or_else(||invalid("fine port attachment outside coarse graph"))?;*d+=p.effective_export_m2;},
            None if [crate::t1b::MARINE,crate::t1b::LAKE].contains(&p.boundary.kind) => {direct+=p.effective_export_m2;},
            _=>return Err(invalid("unsupported coarse/fine port")),
        }}
        // Kahn ordering of actual receiver edges. No hierarchy or basin label
        // replaces the graph; terminals include the declared off-chart receiver.
        let mut indegree=vec![0usize;crate::NODES];for &r in &self.receiver{if (r as usize)<crate::NODES{indegree[r as usize]+=1;}}
        let mut ready:BTreeSet<_>=(0..crate::NODES).filter(|&i|indegree[i]==0).collect();let mut visited=0;let mut terminal=0.0;
        while let Some(i)=ready.pop_first(){visited+=1;let r=self.receiver[i] as usize;if r<crate::NODES{delta[r]+=delta[i];indegree[r]-=1;if indegree[r]==0{ready.insert(r);}}else{terminal+=delta[i];}}
        if visited!=crate::NODES{return Err(invalid("cycle in actual coarse receiver graph"));}
        let residual=(terminal+direct-partition.refinement_delta_m2).abs();if residual>1e-9*partition.fine_own_m2.max(partition.removed_coarse_own_m2).max(1.0){return Err(invalid("coarse/fine replacement source/export failure"));}
        let total:Vec<_>=self.baseline_effective_m2.iter().zip(&delta).map(|(a,b)|a+b).collect();
        if total.iter().any(|a|!a.is_finite()||*a<0.0){return Err(invalid("replacement exceeds a published coarse discharge; baseline quadrature/precision needs explicit reconciliation"));}
        let mut out=WaterComposition{partition,graph_digest:graph.digest,coarse_effective_m2:total,direct_native_water_exports_m2:direct,coarse_terminal_delta_m2:terminal,source_export_residual_m2:residual,coarse_network_composed:true,digest:String::new()};out.digest=content_digest(&out)?;Ok(out)
    }
}
