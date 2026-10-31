//! VS4 native source replacement and whole-cell restriction. These are candidate
//! transactions, not permission to infer rock from an edifice height. Failure
//! leaves the existing E3 construction excess at zero stock.
use pf_contracts::column_envelope::{SourceColumnEnvelope, RetainedVolcanicBatch,
    ReplacementLedger, RetainedReplacement};
use pf_contracts::volcanic_supply::MassUg;
use serde::{Serialize,Deserialize};
use std::collections::BTreeMap;
type R<T> = Result<T,String>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UnitMeasure {
    pub record_id: String,
    pub native_alias: Option<String>,
    pub occupied_full_depth_m3: f64,
    pub occupied_above_sea_m3: f64,
    pub published_h_times_area_m3: f64,
    pub dry_mass_kg: f64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WholeCellMeasure {
    pub cell:u32,
    pub subdivisions:usize,
    pub spherical_area_m2:f64,
    pub sampled_area_m2:f64,
    pub mean_surface_m:f64,
    pub original_mean_m:f64,
    pub mean_residual_m:f64,
    pub positive_excess_m3:f64,
    pub surface_removed_from_envelope_m3:f64,
    pub units:Vec<UnitMeasure>,
    pub geometry_kind:String,
}
/// Integrate on equal-area patches in (sin latitude, longitude). Every patch
/// lies in the native cell; marine areas must be supplied rather than clipped.
/// Volume uses C7's declared thin-shell reference area, not a silently different
/// radial/density law. The sampler declares its quantity in geometry_kind.
/// Water stages may be measured only as diagnostic upper bounds on solid bed;
/// they cannot supply a successful CellRestriction or a positive stock claim.
pub fn restrict_whole_cell(
    e:&SourceColumnEnvelope, n:usize, density:f64, geometry_kind:&str,
    surface:&impl Fn(f64,f64)->R<f64>,
)->R<WholeCellMeasure>{
    e.validate()?;
    if !(2..=8192).contains(&n) || !density.is_finite() || density<=0. || geometry_kind.is_empty(){return Err("invalid spherical restriction inputs".into());}
    let cell=match e.planet_key{pf_contracts::planet::PlanetKey::Cell{cell}=>cell,_=>return Err("not a cell".into())};
    let [la0,la1]=e.support.latitude_bounds_deg.map(f64::to_radians);
    let [lo0,lo1]=e.support.longitude_bounds_deg.map(f64::to_radians);
    let (s0,s1)=(la0.sin(),la1.sin());let radius=e.support.radius_m;
    let area=radius*radius*(lo1-lo0)*(s1-s0);
    if (area-e.support.whole_cell_area_m2).abs()>1e-12*area{return Err("native spherical area disagreement".into());}
    let da=area/(n*n) as f64;let mut sum=0.;let mut excess=0.;let mut removed=0.;
    let mut units:Vec<_>=e.layers.iter().map(|l|UnitMeasure{record_id:l.record_id.clone(),native_alias:l.native_unit_key.clone(),occupied_full_depth_m3:0.,occupied_above_sea_m3:0.,published_h_times_area_m3:l.h_retained_m as f64*area,dry_mass_kg:0.}).collect();
    // Integrate removal from the analytic h*A baseline. Summing the full deep
    // crust at every post would introduce error unrelated to this local edit.
    let mut removed_units=vec![0.;units.len()];let mut removed_correction=vec![0.;units.len()];
    // Compensated sum avoids a mesh-size-dependent height accumulation error.
    let mut correction=0.;
    for y in 0..n {let lat=(s0+(y as f64+0.5)*(s1-s0)/n as f64).asin();
        for x in 0..n {let lon=lo0+(x as f64+0.5)*(lo1-lo0)/n as f64;
            let tilt=envelope_tilt(e,lat.to_degrees(),lon.to_degrees());
            let roof=e.moments.mean_z_m+tilt;
            let z=surface(lat.to_degrees(),lon.to_degrees())?;
            if !z.is_finite(){return Err("unknown/nonfinite solid surface".into());}
            let v=z-correction;let next=sum+v;correction=(next-sum)-v;sum=next;
            excess+=(z-roof).max(0.)*da;removed+=(roof-z).max(0.)*da;
            for (index,(l,u)) in e.layers.iter().zip(&mut units).enumerate(){
                let base=l.base_at_cell_centre_m+tilt;let h=l.h_retained_m as f64;let top=base+h;
                let loss=(top-z).clamp(0.,h)*da-removed_correction[index];
                let next=removed_units[index]+loss;removed_correction[index]=(next-removed_units[index])-loss;removed_units[index]=next;
                u.occupied_above_sea_m3+=(z.min(top)-base.max(e.support.sea_level_m)).max(0.)*da;
            }
        }
    }
    for (u,loss) in units.iter_mut().zip(removed_units) {
        u.occupied_full_depth_m3=u.published_h_times_area_m3-loss;
        u.dry_mass_kg=u.occupied_full_depth_m3*density;
    }
    let mean=sum/(n*n) as f64;
    Ok(WholeCellMeasure{cell,subdivisions:n,spherical_area_m2:area,sampled_area_m2:area,mean_surface_m:mean,original_mean_m:e.moments.mean_z_m,mean_residual_m:mean-e.moments.mean_z_m,positive_excess_m3:excess,surface_removed_from_envelope_m3:removed,units,geometry_kind:geometry_kind.into()})
}
pub fn envelope_tilt(e:&SourceColumnEnvelope,lat_deg:f64,lon_deg:f64)->f64{
    let la=(e.support.latitude_bounds_deg[0]+e.support.latitude_bounds_deg[1])*0.5;
    let lo=(e.support.longitude_bounds_deg[0]+e.support.longitude_bounds_deg[1])*0.5;
    let east=e.support.radius_m*la.to_radians().cos()*((lon_deg-lo+540.).rem_euclid(360.)-180.).to_radians();
    let north=e.support.radius_m*(lat_deg-la).to_radians();
    e.moments.uplift_gradient_east_north[0]*east+e.moments.uplift_gradient_east_north[1]*north
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DiffuseSupportMeasure {
    pub carrier:u64,
    pub batch_id:String,
    pub source_mass:MassUg,
    pub represented_mass:MassUg,
    pub native_carrier_area_m2:f64,
    pub represented_area_m2:f64,
    pub relative_residual:f64,
    pub residual_kg:f64,
    pub tolerance_kg:f64,
    pub quantization_bound_kg:f64,
    pub cells:Vec<u32>,
    pub qualified:bool,
    pub model:String,
}
/// Test the inherited diffuse-copy hypothesis WITHOUT normalizing excess away.
/// This is a falsifiable spatial proposal under the VS3 homogeneous-retention
/// abstraction. A dominant-cell map is not a measured carrier mixing weight.
/// It therefore cannot be used as a successful source admission by itself.
pub fn measure_diffuse_copy(
    b:&RetainedVolcanicBatch, native_area:f64, cells:&[(u32,f64)],
)->R<DiffuseSupportMeasure>{
    b.validate()?;
    if b.remaining == MassUg::ZERO {return Err("zero retained batch cannot fund replacement".into());}
    if !native_area.is_finite() || native_area<=0. || cells.is_empty(){return Err("missing native diffuse support".into());}
    let mut seen=std::collections::BTreeSet::new();let mut area=0.;let mut represented=MassUg::ZERO;
    for &(cell,a) in cells {
        if !seen.insert(cell) || !a.is_finite() || a<=0.{return Err("duplicate/invalid diffuse cell support".into());}
        area+=a;
        let ug=b.remaining.micrograms() as f64*(a/native_area);
        if !ug.is_finite() || ug>=i128::MAX as f64{return Err("diffuse mass overflow".into());}
        represented=represented.checked_add(MassUg::new(ug.round_ties_even() as i128)?)?;
    }
    let kg=b.remaining.micrograms() as f64*1e-9;
    let residual=(represented.micrograms()-b.remaining.micrograms()) as f64*1e-9;
    let quantization=(cells.len() as f64*0.5+16.*f64::EPSILON*represented.micrograms() as f64)*1e-9;
    let tolerance=0.001f64.max(1e-9*kg);
    Ok(DiffuseSupportMeasure{carrier:b.carrier,batch_id:b.batch_id.clone(),source_mass:b.remaining,represented_mass:represented,
        native_carrier_area_m2:native_area,represented_area_m2:area,relative_residual:residual/kg,
        residual_kg:residual,tolerance_kg:tolerance,quantization_bound_kg:quantization,cells:cells.iter().map(|x|x.0).collect(),
        qualified:residual.abs()<=tolerance+quantization,
        model:"candidate homogeneous retained concentration over each dominant-carrier cell; unnormalized full native support; not a recovered carrier-to-cell mixing map".into()})
}

/// Called from the native material owner, never by interpreting a mesh as stock.
/// Registration verifies unique diffuse funding; application is atomic on error.
/// Geometry qualification remains mandatory even when a finite batch is ample.
pub fn prepare_replacement(
    ledger:&ReplacementLedger, proposal:&RetainedReplacement,
)->R<ReplacementLedger>{let mut candidate=ledger.clone();candidate.apply(proposal)?;Ok(candidate)}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_contracts::column_envelope::DiffuseWithdrawal;
    use std::path::Path;
    use serde_json::json;
    #[test]
    fn actual_c7_whole_cell_restriction_and_diffuse_falsifiers(){
        let (Ok(input),Ok(output))=(std::env::var("VZ_VS4_INPUT"),std::env::var("VZ_VS4_OUTPUT"))else{
            eprintln!("SKIP VS4 actual C7: set VZ_VS4_INPUT and VZ_VS4_OUTPUT");return;
        };
        let manifest:serde_json::Value=serde_json::from_slice(&std::fs::read(input).unwrap()).unwrap();
        let read=|key:&str|std::fs::read(manifest[key].as_str().unwrap()).unwrap();
        let envelopes:Vec<SourceColumnEnvelope>=serde_json::from_slice(&read("envelopes")).unwrap();
        let batches:Vec<RetainedVolcanicBatch>=serde_json::from_slice(&read("batches")).unwrap();
        let areas:BTreeMap<u64,f64>=serde_json::from_slice(&read("native_areas")).unwrap();
        assert!(envelopes.iter().all(|e|batches.iter().all(|b|e.stratigraphy_sha256==b.stratigraphy_sha256)));
        let mut geometry=Vec::new();let mut supports:BTreeMap<u64,Vec<(u32,f64)>>=BTreeMap::new();
        for e in &envelopes {
            let cell=match e.planet_key{pf_contracts::planet::PlanetKey::Cell{cell}=>cell,_=>panic!()};
            let roof=|lat,lon|Ok(e.moments.mean_z_m+envelope_tilt(e,lat,lon));
            let mut mesh=Vec::new();
            for n in [32,64,128,256] {mesh.push(restrict_whole_cell(e,n,2800.,"actual E1 retained source planes; not final T1b terrain",&roof).unwrap());}
            for m in &mesh {assert!((m.sampled_area_m2-m.spherical_area_m2).abs()<1e-6);}
            let a=&mesh[2];let b=&mesh[3];
            let max_volume_delta=a.units.iter().zip(&b.units).map(|(a,b)|(a.occupied_full_depth_m3-b.occupied_full_depth_m3).abs()).fold(0.,f64::max);
            let max_mean_delta=(a.mean_surface_m-b.mean_surface_m).abs();
            let [lo,hi]=e.support.latitude_bounds_deg.map(f64::to_radians);
            let spherical_centroid_lat=((hi*hi.sin()+hi.cos())-(lo*lo.sin()+lo.cos()))/(hi.sin()-lo.sin());
            let analytic_mean=e.moments.mean_z_m+e.moments.uplift_gradient_east_north[1]*e.support.radius_m*(spherical_centroid_lat-(lo+hi)*0.5);
            geometry.push(json!({"cell":cell,"refinements":mesh,"last_volume_change_m3":max_volume_delta,"last_mean_change_m":max_mean_delta,
                "analytic_source_plane_mean_m":analytic_mean,"analytic_source_plane_residual_m":analytic_mean-e.moments.mean_z_m}));
            supports.entry(e.native_datum.as_ref().unwrap().carrier).or_default().push((cell,e.support.whole_cell_area_m2));
        }
        let mut measures=Vec::new();let mut missing=Vec::new();let mut falsifiers=0;
        for b in &batches {
            let Some(cells)=supports.get(&b.carrier)else{missing.push(b.batch_id.clone());continue;};
            let m=measure_diffuse_copy(b,areas[&b.carrier],cells).unwrap();
            // Try the unnormalized proposal in the actual exact ledger. Its
            // rejection is evidence, not a substitute geometry pass.
            let owners:Vec<_>=cells.iter().map(|&(cell,a)|DiffuseWithdrawal{cell,record_id:format!("candidate-homogeneous:{}:{cell}",b.batch_id),mass:MassUg::new((b.remaining.micrograms() as f64*a/areas[&b.carrier]).round_ties_even() as i128).unwrap()}).collect();
            let mut ledger=ReplacementLedger::default();let before=ledger.clone();
            assert!(ledger.register(b.clone(),owners).is_err());assert_eq!(ledger,before);falsifiers+=1;
            let mut dup=cells.clone();dup.push(cells[0]);assert!(measure_diffuse_copy(b,areas[&b.carrier],&dup).is_err());falsifiers+=1;
            assert!(!m.qualified);measures.push(m);
        }
        let cfg=crate::t1b_litho::EnvelopeConfig::from_path(manifest["envelope_jsonl"].as_str().unwrap()).unwrap();
        let owned=super::super::bind(super::super::CompactInput::load(Path::new(manifest["native_tile"].as_str().unwrap())).unwrap(),&cfg,false).unwrap();
        assert!(owned.report.zero_stock_construction_excess_m3>0.);
        assert!(super::super::require_evolution_stock(&owned).is_err());falsifiers+=1;
        // A real finite positive proposal, using the actual legacy tile excess.
        // Its declared cost fraction is NOT adopted as a physical placement law.
        // Missing full-cell geometry must still refuse after exact mass closure.
        use pf_contracts::{column_envelope::FundedDeposit,emplacement::{TimingStatus,EruptionBranch},volcanic_supply::TimeInterval};
        let batch=batches.iter().filter(|b|b.carrier==71200).max_by(|a,b|a.interval_my[1].total_cmp(&b.interval_my[1])).unwrap();
        let donor=envelopes.iter().find(|e|e.planet_key==(pf_contracts::planet::PlanetKey::Cell{cell:32652})).unwrap().layers[0].record_id.clone();
        let amount=MassUg::new((owned.report.zero_stock_construction_excess_m3*batch.dry_density_kg_m3*1e9).round_ties_even() as i128).unwrap();
        assert!(amount < batch.remaining);
        let mut ledger=ReplacementLedger::default();
        let opening=vec![DiffuseWithdrawal{cell:32652,record_id:donor.clone(),mass:batch.remaining}];
        ledger.register(batch.clone(),opening.clone()).unwrap();
        let opening_state=ledger.clone();
        assert!(ledger.register(batch.clone(),opening).is_err());assert_eq!(ledger,opening_state);falsifiers+=1;
        let mut proposal=RetainedReplacement{
            id:"vs4-native-tile-excess-proposal/1".into(),event_id:String::new(),batch_id:batch.batch_id.clone(),source_book_sha256:batch.source_book_sha256.clone(),
            observed_interval_my:batch.interval_my,event_timing:TimingStatus::Conditional{
                observed_interval:TimeInterval{start_my:batch.interval_my[0],end_my:batch.interval_my[1]},
                realization_id:"vs4-candidate-only/1".into(),assumption_ref:"conditional reconstruction; individual eruption ancestry unresolved".into()},
            placement_efficiency:amount.micrograms() as f64/batch.remaining.micrograms() as f64,
            placement_efficiency_evidence:"diagnostic demand / finite retained mass; not an adopted physical efficiency or emplacement law".into(),gross:amount,
            diffuse:vec![DiffuseWithdrawal{cell:32652,record_id:donor,mass:amount}],
            deposits:vec![FundedDeposit{id:"vs4-unadmitted-construction-tile/1".into(),cell:32652,record_id:"vs4-proposed-refined-lot/1".into(),
                eruption_branch:EruptionBranch::Lava,mass:amount,occupied_m3:owned.report.zero_stock_construction_excess_m3,
                integration_error_m3:owned.report.quadrature_and_roundoff_bound_m3+0.0005*owned.support_m2.iter().sum::<f64>()}],restrictions:vec![]};
        proposal.event_id=proposal.expected_event_id().unwrap();
        let refusal=ledger.apply(&proposal).unwrap_err();assert_eq!(refusal,"missing or unrelated affected cell restriction");
        assert_eq!(ledger,opening_state);falsifiers+=1;
        let check=|q:&RetainedReplacement,fragment:&str|{
            let mut l=opening_state.clone();let error=l.apply(q).unwrap_err();assert!(error.contains(fragment),"{error}");assert_eq!(l,opening_state);
        };
        let mut q=proposal.clone();q.source_book_sha256="0".repeat(64);check(&q,"unbound");falsifiers+=1;
        let mut q=proposal.clone();q.batch_id="unknown-source".into();check(&q,"unknown finite batch");falsifiers+=1;
        let mut q=proposal.clone();q.gross=batch.remaining.checked_add(MassUg::new(1).unwrap()).unwrap();check(&q,"overdrawn");falsifiers+=1;
        let mut q=proposal.clone();q.deposits.clear();check(&q,"missing sink or source");falsifiers+=1;
        let mut q=proposal.clone();q.deposits[0].occupied_m3*=1.1;check(&q,"occupied volume and dry mass disagree");falsifiers+=1;
        let mut q=proposal.clone();q.deposits[0].eruption_branch=EruptionBranch::Intrusion;check(&q,"invalid or duplicate positive deposit");falsifiers+=1;
        let mut q=proposal.clone();q.event_timing=TimingStatus::Captured{event_ref:batch.native_episode.clone()};q.event_id=q.expected_event_id().unwrap();check(&q,"cannot assert a recovered eruption");falsifiers+=1;
        let mut q=proposal.clone();q.diffuse.push(q.diffuse[0].clone());check(&q,"duplicate or overdrawn diffuse");falsifiers+=1;
        let mut q=proposal.clone();q.id="different-transfer-id".into();assert_eq!(q.expected_event_id().unwrap(),proposal.event_id);falsifiers+=1;
        let mut section=Vec::new();
        let peak=owned.input.height_mm.iter().enumerate().max_by_key(|(_,z)|*z).unwrap().0;
        let row=peak/owned.input.rect.width;
        let env=cfg.load().unwrap();let chart:vz_tilekey::chart::ChartDeclV1=serde_json::from_value(owned.input.header["chart"].clone()).unwrap();
        let chart=vz_tilekey::chart::Chart::new(chart).unwrap();
        for x in 0..owned.input.rect.width{
            let i=row*owned.input.rect.width+x;let p=owned.input.rect.id(i);let e=&env.cells[&owned.input.cells[i]];
            let (la,lo)=chart.inverse_mm(p.gx*32000,p.gy*32000).unwrap();
            section.push(json!({"x_m":p.gx*32,"solid_m":owned.input.height_mm[i] as f64/1000.,"retained_roof_m":e.mean_m+e.tilt(owned.input.cells[i],la.to_degrees(),lo.to_degrees()),"funded_construction_m":0.,"zero_stock_excess_m":owned.excess_depth_m[i]}));
        }
        let out=json!({"schema":"pf.vs4-restriction-audit/1","whole_cell_source_geometry":geometry,"diffuse_proposals":measures,
            "unknown_marine_batches":missing,"native_32m_ownership":owned.report,"section_north_m":(owned.input.rect.gy0+row as i64)*32,"section":section,"falsifiers_passed":falsifiers,
            "funded_cone_volume_m3":0.,"committed_transfers":0,"positive_proposal":proposal,"positive_proposal_refusal":refusal,
            "proposal_opening_spatial_status":"conditional isolated scalar owner for refusal test; no full-cell native allocation qualified",
            "proposal_integer_selection_rounding_bound_kg":16.*f64::EPSILON*amount.micrograms() as f64*1e-9+0.5e-9,
            "proposal_ledger_residual_ug":"0","proposal_ledger_tolerance_kg":0.001f64.max(1e-9*amount.micrograms() as f64*1e-9),
            "status":"not qualified: inherited diffuse allocation overdraws finite native support; no normalization or datum fit; final whole-cell T1b terrain and displacement remain unknown"});
        let f=std::fs::OpenOptions::new().write(true).create_new(true).open(output).unwrap();serde_json::to_writer_pretty(f,&out).unwrap();
        println!("VS4 actual C7: {} whole source cells, {} failed finite allocations, {} marine unknown batches, {} real falsifier/invariant checks; funded cone volume zero",envelopes.len(),out["diffuse_proposals"].as_array().unwrap().len(),out["unknown_marine_batches"].as_array().unwrap().len(),falsifiers);
    }
}

/// Authenticated native 32 m construction, sampled as the same half-open post
/// prisms used by E3. All tiles must share their context and envelope pin.
/// This compact export precedes the window dam-lake classifier and does not
/// retain its pre-fill beds. All heights are conservative construction-top
/// upper bounds on solid, even when the explicit water flag is absent.
pub struct NativeMosaic {
    tiles:BTreeMap<(i64,i64),super::CompactInput>,
    chart:vz_tilekey::chart::Chart,
    pub context_pin:String,
    pub envelope_digest:String,
}
impl NativeMosaic {
    pub fn load(root:&std::path::Path)->R<Self>{
        let mut tiles=BTreeMap::new();let mut chart=None;let mut context=None;let mut envelope=None;
        for entry in std::fs::read_dir(root).map_err(|e|e.to_string())? {
            let entry=entry.map_err(|e|e.to_string())?;
            if !entry.path().join("material.json").is_file(){continue;}
            let t=super::CompactInput::load(&entry.path())?;
            let p=t.header["native_context_pin"].as_str().ok_or("missing native context")?.to_owned();
            let e=t.header["source_envelope"]["digest"].as_str().ok_or("missing native source envelope")?.to_owned();
            let c:vz_tilekey::chart::ChartDeclV1=serde_json::from_value(t.header["chart"].clone()).map_err(|e|e.to_string())?;
            if context.as_ref().is_some_and(|old| old != &p) || envelope.as_ref().is_some_and(|old| old != &e)
                || chart.as_ref().is_some_and(|old| old != &c){return Err("mixed native material contexts".into());}
            context=Some(p);envelope=Some(e);chart=Some(c);
            let key=(t.rect.gx0.div_euclid(512),t.rect.gy0.div_euclid(512));
            if tiles.insert(key,t).is_some(){return Err("duplicate native material owner".into());}
        }
        Ok(Self{tiles,chart:vz_tilekey::chart::Chart::new(chart.ok_or("empty native material")?).map_err(|e|e.to_string())?,context_pin:context.unwrap(),envelope_digest:envelope.unwrap()})
    }
    pub fn sample(&self,lat:f64,lon:f64)->R<Option<(f64,u8,u32)>>{
        let xy=self.chart.forward_f64(lat.to_radians(),lon.to_radians()).map_err(|e|e.to_string())?;
        // forward_f64 returns millimetres. Boundary posts own only their
        // in-chart half prism; do not extend them into unknown exterior.
        let half=vz_tilekey::CHART_HALF_EXTENT_MM as f64;
        if xy[0].abs()>half || xy[1].abs()>half{return Ok(None);}
        let gx=(xy[0]/32000.+0.5).floor() as i64;let gy=(xy[1]/32000.+0.5).floor() as i64;
        let Some(t)=self.tiles.get(&(gx.div_euclid(512),gy.div_euclid(512)))else{return Ok(None);};
        let Some(i)=t.rect.index(crate::t1b_surface_domain::PostId{gx,gy})else{return Ok(None);};
        Ok(Some((t.height_mm[i] as f64/1000.,t.fixed[i],t.cells[i])))
    }
    pub fn tile_count(&self)->usize{self.tiles.len()}
}

#[cfg(test)]
mod native_measure_tests {
    use super::*;
    use serde_json::json;
    #[test]
    fn actual_native_whole_cells_include_cones_and_off_cone_terrain(){
        let (Ok(input),Ok(root),Ok(output))=(std::env::var("VZ_VS4_INPUT"),std::env::var("VZ_VS4_MOSAIC"),std::env::var("VZ_VS4_MOSAIC_OUTPUT"))else{
            eprintln!("SKIP VS4 native whole cells: set VZ_VS4_INPUT, VZ_VS4_MOSAIC, VZ_VS4_MOSAIC_OUTPUT");return;
        };
        let manifest:serde_json::Value=serde_json::from_slice(&std::fs::read(input).unwrap()).unwrap();
        let envelopes:Vec<SourceColumnEnvelope>=serde_json::from_slice(&std::fs::read(manifest["envelopes"].as_str().unwrap()).unwrap()).unwrap();
        let mosaic=NativeMosaic::load(std::path::Path::new(&root)).unwrap();
        let cfg=crate::t1b_litho::EnvelopeConfig::from_path(manifest["envelope_jsonl"].as_str().unwrap()).unwrap();
        assert_eq!(mosaic.envelope_digest,cfg.digest);
        let mut cells=Vec::new();
        for e in &envelopes {
            let mut levels=Vec::new();
            for n in [128,256,512,1024,2048,4096,8192] {
                let water=std::cell::Cell::new(0usize);let mismatch=std::cell::Cell::new(0usize);
                let own=match e.planet_key{pf_contracts::planet::PlanetKey::Cell{cell}=>cell,_=>panic!()};
                let sample=|la,lo|->R<f64>{
                    let Some((z,f,c))=mosaic.sample(la,lo)? else {return Err("native chart/exterior solid missing".into());};
                    if f==crate::t1b::MARINE || f==crate::t1b::LAKE {water.set(water.get()+1);}
                    if c != own {mismatch.set(mismatch.get()+1);}
                    Ok(z)
                };
                let measure=restrict_whole_cell(e,n,2800.,"native32m post-prism construction including cones; upper bound on solid because compact export omits dam-lake pre-fill beds",&sample);
                match measure {
                    Ok(m)=>levels.push(json!({"measure":m,"water_stage_samples":water.get(),"post_owner_boundary_samples":mismatch.get(),
                        "mean_interpretation":"upper bound on solid mean: explicit water stages and possible later-classified dam lakes; compact export lacks pre-fill beds; never a solid-volume pass"})),
                    Err(reason)=>{
                        // Retain partial support quantitatively without using a
                        // land-only denominator as a whole-cell mean.
                        let [la0,la1]=e.support.latitude_bounds_deg.map(f64::to_radians);
                        let [lo0,lo1]=e.support.longitude_bounds_deg.map(f64::to_radians);
                        let (s0,s1)=(la0.sin(),la1.sin());let da=e.support.whole_cell_area_m2/(n*n) as f64;
                        let(mut known,mut water_count,mut z_area,mut excess)=(0usize,0usize,0.,0.);
                        for y in 0..n {let la=(s0+(y as f64+0.5)*(s1-s0)/n as f64).asin().to_degrees();
                            for x in 0..n {let lo=(lo0+(x as f64+0.5)*(lo1-lo0)/n as f64).to_degrees();
                                if let Some((z,f,_))=mosaic.sample(la,lo).unwrap(){
                                    known+=1;z_area+=z*da;
                                    if f==crate::t1b::MARINE || f==crate::t1b::LAKE {water_count+=1;}
                                    else {excess+=(z-e.moments.mean_z_m-envelope_tilt(e,la,lo)).max(0.)*da;}
                                }
                            }
                        }
                        levels.push(json!({"subdivisions":n,"unknown":reason,"whole_cell_area_m2":e.support.whole_cell_area_m2,
                            "known_area_m2":known as f64*da,"unknown_area_m2":(n*n-known) as f64*da,
                            "water_stage_area_m2":water_count as f64*da,"published_surface_integral_on_known_support_m3":z_area,
                            "land_construction_excess_on_known_support_m3":excess,"whole_cell_solid_mean_m":null}));
                    }
                }
            }
            cells.push(json!({"cell":e.planet_key,"levels":levels}));
        }
        let plan:serde_json::Value=serde_json::from_slice(&std::fs::read(manifest["native_plan"].as_str().unwrap()).unwrap()).unwrap();
        let mut marine=Vec::new();
        for c in plan["cells"].as_array().unwrap().iter().filter(|c|c["has_land_envelope"]==false) {
            let cell=c["cell"].as_u64().unwrap() as u32;let grid=envelopes[0].support.shape;let radius=envelopes[0].support.radius_m;
            let la0=(-90.+180.*(cell/grid[1]) as f64/grid[0] as f64).to_radians();
            let la1=la0+std::f64::consts::PI/grid[0] as f64;
            let lo0=(-180.+360.*(cell%grid[1]) as f64/grid[1] as f64).to_radians();
            let dlo=2.*std::f64::consts::PI/grid[1] as f64;
            let area=radius*radius*dlo*(la1.sin()-la0.sin());let mut levels=Vec::new();
            for n in [128,256,512,1024] {
                let(mut known,mut water,mut sum)=(0usize,0usize,0.);
                for y in 0..n {let la=(la0.sin()+(y as f64+0.5)*(la1.sin()-la0.sin())/n as f64).asin().to_degrees();
                    for x in 0..n {let lo=(lo0+(x as f64+0.5)*dlo/n as f64).to_degrees();
                        if let Some((z,f,_))=mosaic.sample(la,lo).unwrap(){known+=1;sum+=z;if f==crate::t1b::MARINE || f==crate::t1b::LAKE {water+=1;}}
                    }
                }
                levels.push(json!({"subdivisions":n,"whole_cell_area_m2":area,"known_area_m2":area*known as f64/(n*n) as f64,
                    "water_stage_area_m2":area*water as f64/(n*n) as f64,"published_height_mean_on_known_support_m":if known>0{Some(sum/known as f64)}else{None},
                    "whole_cell_solid_mean_m":null,"source_unit_occupancy":null}));
            }
            marine.push(json!({"cell":cell,"carrier":c["carrier"],"levels":levels,
                "status":"unknown: no E1 marine source envelope; regional water stage is not a seabed; exterior chart support missing"}));
        }
        let out=json!({"schema":"pf.vs4-native-whole-cells/1","context_pin":mosaic.context_pin,"envelope":mosaic.envelope_digest,
            "tiles":mosaic.tile_count(),"cells":cells,"marine_cells":marine,"note":"No datum fit, no edifice mask, no land-only denominator. Every construction-top mean is an upper bound on solid: compact material precedes the window dam-lake classifier and omits pre-fill beds. Explicit water flags do not exhaust lake support. Missing chart exterior is unknown. Negative complete-cell upper-bound residuals falsify the original mean without inventing beds; positive upper bounds do not establish a solid mean. This is restriction evidence, not a source debit or terrain acceptance."});
        let f=std::fs::OpenOptions::new().write(true).create_new(true).open(output).unwrap();serde_json::to_writer_pretty(f,&out).unwrap();
        println!("VS4 native mosaic: {} authenticated tiles; {} source cells measured/refused with explicit support",mosaic.tile_count(),envelopes.len());
    }
}
