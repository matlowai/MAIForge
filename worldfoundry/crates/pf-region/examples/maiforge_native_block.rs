//! Export one native process-grid block, not an oversampled display mesh.
//! CLI: maiforge_native_block WORLD_ROOT CELL OUTPUT_DIR [MAX_WORK_NODES]
use pf_region::planet_relief::WorldRelief;
use serde_json::json;
use std::{path::Path, time::Instant};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a: Vec<String> = std::env::args().collect();
    if !(4..=5).contains(&a.len()) {
        return Err(
            "usage: maiforge_native_block WORLD_ROOT CELL OUTPUT_DIR [MAX_WORK_NODES=2000000]"
                .into(),
        );
    }
    let root = Path::new(&a[1]);
    let cell = a[2].parse::<u32>()?;
    let out = Path::new(&a[3]);
    if out.exists() {
        return Err("output exists; choose a fresh directory".into());
    }
    let parent = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
        .canonicalize()?;
    if parent
        .join(out.file_name().ok_or("output directory")?)
        .starts_with(root.canonicalize()?)
    {
        return Err("output must remain outside source world".into());
    }
    let budget = a
        .get(4)
        .map(|x| x.parse::<usize>())
        .transpose()?
        .unwrap_or(2_000_000);
    let start = Instant::now();
    let world = WorldRelief::from_published(root).map_err(std::io::Error::other)?;
    let load_s = start.elapsed().as_secs_f64();
    let block = world
        .process_block(cell, budget)
        .map_err(std::io::Error::other)?;
    let mut sample_error = 0.0f64;
    let probes = [
        (0, 0),
        (129, 0),
        (0, 129),
        (129, 129),
        (1, 1),
        (64, 64),
        (128, 128),
        (0, 64),
        (129, 64),
        (64, 0),
        (64, 129),
    ];
    for (x, y) in probes {
        let gx = (block.origin_post[0] + x as isize).rem_euclid(block.process_grid[1] as isize);
        let gy = (block.origin_post[1] + y as isize).clamp(0, block.process_grid[0] as isize - 1);
        let lon = -180. + (gx as f64 + 0.5) * 360. / block.process_grid[1] as f64;
        let lat = -90. + (gy as f64 + 0.5) * 180. / block.process_grid[0] as f64;
        let z = world
            .sample(lat, lon)
            .map_err(std::io::Error::other)?
            .height_m;
        sample_error = sample_error.max((z - block.height_m[y * block.n + x]).abs());
    }
    if sample_error > 1e-6 {
        return Err("native block disagrees with canonical owner sampling".into());
    }
    let canonical: Vec<u8> = block
        .height_m
        .iter()
        .flat_map(|x| x.to_le_bytes())
        .collect();
    let mut max_quant = 0.0f64;
    let mut bytes = Vec::with_capacity(block.height_m.len() * 4);
    for &x in &block.height_m {
        let y = x as f32;
        if !y.is_finite() {
            return Err("height cannot be transferred as f32".into());
        }
        max_quant = max_quant.max((x - y as f64).abs());
        bytes.extend(y.to_le_bytes());
    }
    let digest = |b: &[u8]| format!("blake3:{}", blake3::hash(b));
    let summaries:Vec<_>=block.process_receipts.as_array().ok_or("receipt array")?.iter().map(|r|json!({
        "outletCell":r["outlet_cell"],"coreCells":r["core_cells"],"haloCells":r["halo_cells"],"nativeNodes":r["physical_nodes"],
        "waterRelativeResidual":r["final_routing"]["relative_residual"],"transportRelativeResidual":r["final_transport"]["relative_residual"],
        "bulkRelativeResidual":r["geometric_bulk_relative_residual"],"diffusionRelativeResidual":r["maximum_diffusion_balance_residual"],
        "maximumDryDiagnosticFillM":r["final_maximum_dry_fill_m"]})).collect();
    let mut m = json!({"schema":block.schema,"model":block.model,"cell":cell,"nativeGrid":block.native_grid,"processGrid":block.process_grid,
        "resolution":block.resolution,"apron":block.apron,"shape":[block.n,block.n],"originPost":block.origin_post,
        "sampling":"global cell-centred process nodes; south-first rows, longitude periodic",
        "interpolation":"bilinear; polar half-row blends to native poleHeightM; apron y clamps to valid node rows",
        "poleHeightM":block.pole_height_m,"sourceResolutionM":world.input.radius_m*std::f64::consts::PI/block.process_grid[0] as f64,
        "radiusM":world.input.radius_m,"grid":{"rows":block.n,"cols":block.n,"order":"south-first","sampling":"process-cell-centres"},
        "fields":{"elevation":{"file":"height.f32","type":"f32le","unit":"m","bytes":bytes.len(),"sha256":pf_hydro_catchment::input::sha256_bytes(&bytes),"digest":digest(&bytes),"canonicalF64Digest":digest(&canonical),"maxQuantizationErrorM":max_quant}},
        "source":block.source,"processReceiptSummaries":summaries,"work":block.work,"limits":block.limits,
        "canonicalSampleChecks":{"count":probes.len(),"maximumErrorM":sample_error},
        "materialQualification":"Height-only rendering data; no finite-stock admission, exposed-layer mapping or gameplay authority is supplied."});
    m["id"] = digest(&serde_json::to_vec(&m)?).into();
    std::fs::create_dir(out)?;
    std::fs::write(out.join("height.f32"), bytes)?;
    std::fs::write(
        out.join("process-receipts.json"),
        serde_json::to_vec_pretty(&block.process_receipts)?,
    )?;
    std::fs::write(
        out.join("execution.json"),
        serde_json::to_vec_pretty(
            &json!({"loadSeconds":load_s,"totalSeconds":start.elapsed().as_secs_f64(),"identity":"Execution timing is excluded from the data identity."}),
        )?,
    )?;
    std::fs::write(out.join("manifest.json"), serde_json::to_vec_pretty(&m)?)?;
    println!(
        "{}",
        serde_json::to_string_pretty(
            &json!({"id":m["id"],"cell":cell,"shape":m["shape"],"fields":m["fields"],"work":m["work"],"totalSeconds":start.elapsed().as_secs_f64()})
        )?
    );
    Ok(())
}
