//! Export a bounded 5x5 region using one persistent native sampler.
//! No source fields, terrain math or canonical process-node ownership are changed.
use pf_region::planet_relief::WorldRelief;
use serde_json::{json,Value};
use std::{collections::BTreeMap,path::Path,time::Instant};
fn hash(b:&[u8])->String {pf_hydro_catchment::input::sha256_bytes(b)}
fn blake(b:&[u8])->String {format!("blake3:{}",blake3::hash(b))}
fn sample_map(nodes:&BTreeMap<(usize,usize),f64>,w:&WorldRelief,lon:f64,lat:f64)->Result<f64,String> {
    let r=pf_region::planet_catchment::RESOLUTION;let nx=w.input.cols*r;let ny=w.input.rows*r;
    let x=(lon+180.).rem_euclid(360.)/360.*nx as f64-0.5;let y=(lat+90.)/180.*ny as f64-0.5;
    let i=x.floor() as isize;let j=y.floor() as isize;let fx=x-x.floor();let fy=y-y.floor();
    let at=|i:isize,j:isize|nodes.get(&(i.rem_euclid(nx as isize) as usize,j.clamp(0,ny as isize-1) as usize)).copied().ok_or("boundary node absent".to_string());
    Ok((1.-fy)*((1.-fx)*at(i,j)?+fx*at(i+1,j)?)+fy*((1.-fx)*at(i,j+1)?+fx*at(i+1,j+1)?))
}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let a:Vec<String>=std::env::args().collect();if !(4..=5).contains(&a.len()){return Err("usage: maiforge_native_region WORLD_ROOT CENTRE_CELL OUTPUT_DIR [MAX_WORK_NODES_PER_BLOCK=8000000]".into());}
    let root=Path::new(&a[1]);let centre=a[2].parse::<u32>()?;let out=Path::new(&a[3]);
    if out.exists(){return Err("output exists; choose a fresh directory".into());}
    let parent=out.parent().filter(|p|!p.as_os_str().is_empty()).unwrap_or_else(||Path::new(".")).canonicalize()?;
    if parent.join(out.file_name().ok_or("output directory")?).starts_with(root.canonicalize()?){return Err("output must remain outside source world".into());}
    let budget=a.get(4).map(|s|s.parse::<usize>()).transpose()?.unwrap_or(8_000_000);let started=Instant::now();
    let world=WorldRelief::from_published(root).map_err(std::io::Error::other)?;let load_s=started.elapsed().as_secs_f64();
    let rows=world.input.rows;let cols=world.input.cols;let r=pf_region::planet_catchment::RESOLUTION;
    if centre as usize>=rows*cols||rows<5||cols<10 {return Err("unsupported region source grid or centre cell".into());}
    let row=centre as usize/cols;let col=centre as usize%cols;
    if row<2||row+2>=rows {return Err("complete 5x5 region would exceed source latitude bounds".into());}
    let row_start=row-2;let col_start=col as isize-2;
    let mut requests=Vec::new();let mut preflight=Vec::new();let mut bound_failure=false;
    for y in 0..5 {for x in 0..5 {let cell=((row_start+y)*cols+(col_start+x).rem_euclid(cols as isize) as usize) as u32;let work=world.process_block_work(cell).map_err(std::io::Error::other)?;let count:usize=work.iter().map(|w|w.physical_nodes).sum();if count>budget {bound_failure=true;}
        requests.push(cell);preflight.push(json!({"cell":cell,"workNodes":count,"withinBudget":count<=budget,"work":work}));}}
    std::fs::create_dir(out)?;
    std::fs::write(out.join("preflight.json"),serde_json::to_vec_pretty(&json!({"budgetPerBlock":budget,"completeCoverage":25,"blocks":preflight,"loadSeconds":load_s}))?)?;
    if bound_failure{return Err("region work bound exceeded; no blocks solved and no region manifest published".into());}
    let mut nodes=BTreeMap::<(usize,usize),f64>::new();let mut blocks=Vec::new();let mut shared=0usize;let mut max_shared=0.0f64;let mut source:Option<Value>=None;
    let mut maximum_quantization=0.0f64;let mut minimum=f64::INFINITY;let mut maximum=f64::NEG_INFINITY;let mut dry_max:Option<(f64,f64,f64)>=None;
    for (ordinal,cell) in requests.into_iter().enumerate() {
        let block=world.process_block(cell,budget).map_err(std::io::Error::other)?;
        if let Some(s)=&source {if *s!=block.source{return Err("block source changed during region export".into());}} else {source=Some(block.source.clone());}
        let mut canonical_error=0.0f64;let probes=[(0,0),(129,0),(0,129),(129,129),(1,1),(64,64),(128,128),(0,64),(129,64),(64,0),(64,129)];
        for (x,y) in probes {let gx=(block.origin_post[0]+x as isize).rem_euclid(block.process_grid[1] as isize);let gy=(block.origin_post[1]+y as isize).clamp(0,block.process_grid[0] as isize-1);let lon=-180.+(gx as f64+0.5)*360./block.process_grid[1] as f64;let lat=-90.+(gy as f64+0.5)*180./block.process_grid[0] as f64;
            let z=world.sample(lat,lon).map_err(std::io::Error::other)?.height_m;canonical_error=canonical_error.max((z-block.height_m[y*block.n+x]).abs());}
        if canonical_error>1e-6{return Err("native canonical sample check failed".into());}
        let mut bytes=Vec::with_capacity(block.height_m.len()*4);let canonical:Vec<u8>=block.height_m.iter().flat_map(|x|x.to_le_bytes()).collect();let mut quant=0.0f64;
        for (index,&z) in block.height_m.iter().enumerate() {
            let x=index%block.n;let y=index/block.n;let gx=(block.origin_post[0]+x as isize).rem_euclid(block.process_grid[1] as isize) as usize;let gy=(block.origin_post[1]+y as isize).clamp(0,block.process_grid[0] as isize-1) as usize;
            if let Some(previous)=nodes.insert((gx,gy),z){shared+=1;max_shared=max_shared.max((z-previous).abs());if z.to_bits()!=previous.to_bits(){return Err("shared canonical apron changed between blocks".into());}}
            let q=z as f32;if !q.is_finite(){return Err("nonfinite f32 transfer".into());}quant=quant.max((z-q as f64).abs());bytes.extend(q.to_le_bytes());minimum=minimum.min(z);maximum=maximum.max(z);
            if x>=1&&x<=r&&y>=1&&y<=r&&world.input.land[cell as usize] {let lat=-90.+(gy as f64+0.5)*180./block.process_grid[0] as f64;let lon=-180.+(gx as f64+0.5)*360./block.process_grid[1] as f64;if dry_max.map_or(true,|v|z>v.0){dry_max=Some((z,lat,lon));}}
        }
        maximum_quantization=maximum_quantization.max(quant);
        let summaries:Vec<_>=block.process_receipts.as_array().ok_or("receipts")?.iter().map(|v|json!({"outletCell":v["outlet_cell"],"coreCells":v["core_cells"],"haloCells":v["halo_cells"],"nativeNodes":v["physical_nodes"],"waterRelativeResidual":v["final_routing"]["relative_residual"],"transportRelativeResidual":v["final_transport"]["relative_residual"],"bulkRelativeResidual":v["geometric_bulk_relative_residual"],"diffusionRelativeResidual":v["maximum_diffusion_balance_residual"],"maximumDryDiagnosticFillM":v["final_maximum_dry_fill_m"]})).collect();
        let mut manifest=json!({"schema":block.schema,"model":block.model,"cell":cell,"nativeGrid":block.native_grid,"processGrid":block.process_grid,"resolution":r,"apron":1,"shape":[block.n,block.n],"originPost":block.origin_post,
          "sampling":"global cell-centred process nodes; south-first rows, longitude periodic","interpolation":"bilinear; polar half-row blends to native poleHeightM; apron y clamps to valid node rows","poleHeightM":block.pole_height_m,
          "sourceResolutionM":world.input.radius_m*std::f64::consts::PI/block.process_grid[0] as f64,"radiusM":world.input.radius_m,"grid":{"rows":block.n,"cols":block.n,"order":"south-first","sampling":"process-cell-centres"},
          "fields":{"elevation":{"file":"height.f32","type":"f32le","unit":"m","bytes":bytes.len(),"sha256":hash(&bytes),"digest":blake(&bytes),"canonicalF64Digest":blake(&canonical),"maxQuantizationErrorM":quant}},
          "source":block.source,"processReceiptSummaries":summaries,"work":block.work,"limits":block.limits,"canonicalSampleChecks":{"count":probes.len(),"maximumErrorM":canonical_error},"materialQualification":"Height-only rendering data; no finite-stock admission, exposed-layer mapping or gameplay authority is supplied."});
        manifest["id"]=blake(&serde_json::to_vec(&manifest)?).into();
        let dir=out.join("blocks").join(format!("cell-{cell}"));std::fs::create_dir_all(&dir)?;
        std::fs::write(dir.join("height.f32"),bytes)?;let manifest_bytes=serde_json::to_vec_pretty(&manifest)?;std::fs::write(dir.join("manifest.json"),&manifest_bytes)?;
        blocks.push(json!({"cell":cell,"file":format!("blocks/cell-{cell}/manifest.json"),"sha256":hash(&manifest_bytes)}));
        eprintln!("block {}/25 cell={} elapsed={:.3}s",ordinal+1,cell,started.elapsed().as_secs_f64());
    }
    let source=source.ok_or("no region source")?;let south=-90.+row_start as f64*180./rows as f64;let north=-90.+(row_start+5) as f64*180./rows as f64;let west=-180.+col_start as f64*360./cols as f64;let east=-180.+(col_start+5) as f64*360./cols as f64;
    let mut boundary=Vec::new();let mut max_delta=0.0f64;let mut max_dry_delta=0.0f64;let mut count=0usize;let mut dry_count=0usize;
    for edge in 0..4 {for k in 0..=5*r {let t=k as f64/(5*r) as f64;let (lon,lat)=match edge {0=>(west+(east-west)*t,south),1=>(west+(east-west)*t,north),2=>(west,south+(north-south)*t),_=>(east,south+(north-south)*t)};let lon=(lon+180.).rem_euclid(360.)-180.;let fine=sample_map(&nodes,&world,lon,lat).map_err(std::io::Error::other)?;let base=world.datum(lat,lon);let delta=fine-base;let dry=world.input.land[world.cell(lat,lon) as usize];count+=1;max_delta=max_delta.max(delta.abs());if dry{dry_count+=1;max_dry_delta=max_dry_delta.max(delta.abs());}boundary.push(json!({"longitude":lon,"latitude":lat,"nativeSolidM":fine,"baseSolidM":base,"differenceM":delta,"coarseDry":dry}));}}
    let (clat,clon)=world.centre(centre as usize);let centre_height=sample_map(&nodes,&world,clon,clat).map_err(std::io::Error::other)?;
    let audit=json!({"sharedNodeComparisons":shared,"maximumSharedNodeDifferenceM":max_shared,"heightRangeM":[minimum,maximum],"maxQuantizationErrorM":maximum_quantization,"boundarySamples":count,"dryBoundarySamples":dry_count,"maximumSolidBoundaryDifferenceM":max_delta,"maximumDrySolidBoundaryDifferenceM":max_dry_delta,"centreRawHeightM":centre_height,"dryMaximum":dry_max.map(|(z,lat,lon)|json!({"heightM":z,"latitude":lat,"longitude":lon})),"comparison":"Native solid versus original bilinear-pole-mean solid datum; neither water-surface substitution nor rendering guard is applied.","boundary":boundary});
    std::fs::write(out.join("boundary-audit.json"),serde_json::to_vec_pretty(&audit)?)?;
    let mut region=json!({"schema":"maiforge.viewer-native-region/v1","model":pf_region::planet_catchment::MODEL,"seed":world.input.seed,"epochMy":source["input"]["epochMy"],"radiusM":world.input.radius_m,"nativeGrid":[rows,cols],"processGrid":[rows*r,cols*r],"resolution":r,"apron":1,"sourceResolutionM":world.input.radius_m*std::f64::consts::PI/(rows*r) as f64,
        "coarseRectangle":{"rowStart":row_start,"colStart":col_start,"rowCount":5,"colCount":5},"sourceArtifacts":{"terrain":source["input"]["terrainArtifact"],"hydrology":source["input"]["hydrologyArtifact"]},"blocks":blocks,"source":source,
        "quality":{"maximumSolidBoundaryDifferenceM":max_delta,"maximumDrySolidBoundaryDifferenceM":max_dry_delta,"sharedNodeComparisons":shared,"maximumSharedNodeDifferenceM":max_shared},
        "limits":"Raw conditional native solid-height region. No coarse/fine blending, water-surface replacement, finite-material admission or physical-authority claim. Rendering transition is a separate labelled proxy."});
    region["id"]=blake(&serde_json::to_vec(&region)?).into();std::fs::write(out.join("region.json"),serde_json::to_vec_pretty(&region)?)?;
    std::fs::write(out.join("execution.json"),serde_json::to_vec_pretty(&json!({"sourceLoadSeconds":load_s,"totalSeconds":started.elapsed().as_secs_f64(),"blockCount":25,"note":"Timings excluded from source/data identity."}))?)?;
    println!("{}",serde_json::to_string_pretty(&json!({"id":region["id"],"blocks":25,"boundaryMaxM":max_delta,"dryBoundaryMaxM":max_dry_delta,"sharedMaxM":max_shared,"centreRawHeightM":centre_height,"totalSeconds":started.elapsed().as_secs_f64()}))?);Ok(())
}
