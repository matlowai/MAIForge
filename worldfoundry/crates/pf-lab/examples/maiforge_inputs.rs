//! Generate native upstream planet products and terrain-consumer sidecars.
//! The small default recipe is an interface example, not Earth calibration.
use pf_contracts::artifact::digest_canonical;
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_lab::preview::{self, PreviewConfig, PreviewEnv, PreviewOptions, PreviewRun};
use serde_json::json;
use std::{path::Path, time::Instant};
type R<T> = Result<T, Box<dyn std::error::Error>>;
fn f32_file(p: &Path, values: &[f32]) -> R<()> {
    std::fs::write(p, values.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<_>>())?; Ok(())
}
fn main() -> R<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 { return Err("usage: maiforge_inputs OUTPUT_DIRECTORY".into()); }
    let root = Path::new(&args[1]);
    if root.exists() { return Err("output exists; choose a new directory to preserve prior results".into()); }
    let started=Instant::now();
    let config=PreviewConfig { scenario:ScenarioPreset::toy(), fidelity:FidelityPreset::Preview,
        grid:Some((36,72)), seed:240907, duration_my:Some(4.0), step_my:Some(2.0) };
    let options=PreviewOptions::world_strata_v2();
    let profile=config.profile_with(&options)?;
    let universe_digest=digest_canonical(&preview::universe_of(&profile)?)?;
    let mut run=PreviewRun::init(PreviewEnv::new(profile, config.grid, config.seed_path())?)?;
    while run.steps_taken()<run.steps_scheduled() {run.step()?;}
    let depth=run.state().thickness_m.iter().zip(&run.state().sediment_m)
        .map(|(r,s)|r+s).fold(0.0f64,f64::max)+1.0;
    let full=run.stratigraphy(Some(depth))?.ok_or("selected profile produced no stratigraphy")?;
    let relief:Vec<f32>=run.subgrid_relief_m().ok_or("selected profile produced no subgrid relief")?
        .iter().map(|v|*v as f32).collect();
    let out=run.finish()?;
    std::fs::create_dir_all(root.join("final"))?;
    let digests=preview::publish(&out,&root.join("run"),started)?;
    f32_file(&root.join("subgrid_relief_m.f32"),&relief)?;
    let final_dir=root.join("final");
    std::fs::write(final_dir.join("stratafull_off.u32"), full.layer_offset.iter().flat_map(|v|v.to_le_bytes()).collect::<Vec<_>>())?;
    std::fs::write(final_dir.join("stratafull_lith.u8"),&full.lith)?;
    std::fs::write(final_dir.join("stratafull_env.u8"),&full.env)?;
    f32_file(&final_dir.join("stratafull_th.f32"),&full.thickness_m)?;
    f32_file(&final_dir.join("stratafull_age_base.f32"),&full.age_base_my)?;
    f32_file(&final_dir.join("stratafull_age_top.f32"),&full.age_top_my)?;
    f32_file(&final_dir.join("strata_grad_e.f32"),full.uplift_grad_east.as_ref().ok_or("missing native east gradient")?)?;
    f32_file(&final_dir.join("strata_grad_n.f32"),full.uplift_grad_north.as_ref().ok_or("missing native north gradient")?)?;
    std::fs::write(final_dir.join("stratigraphy.json"),serde_json::to_vec_pretty(&full)?)?;
    // A deterministic example request chosen from generated land, not a production
    // geographic restriction or a hand-authored terrain location.
    let cell=(0..out.grid.cells()).filter(|&c|out.elevation_m[c]>out.terrain.sea_level_m)
        .filter(|&c|out.routing.water_body_id[c]==u32::MAX)
        .filter(|&c|out.grid.center_lat_deg(c/out.grid.nlon()).unwrap().abs()<60.0)
        .max_by(|&a,&b|out.routing.discharge_proxy_m3_s[a].total_cmp(&out.routing.discharge_proxy_m3_s[b]).then(b.cmp(&a)))
        .ok_or("recipe has no suitable generated land for the example")?;
    let lat=out.grid.center_lat_deg(cell/out.grid.nlon())?;
    let lon=out.grid.center_lon_deg(cell%out.grid.nlon())?;
    let header:serde_json::Value=serde_json::from_slice(&std::fs::read(root.join("run/terrain/artifact.json"))?)?;
    let hydro:serde_json::Value=serde_json::from_slice(&std::fs::read(root.join("run/hydrology/artifact.json"))?)?;
    let site=json!({"seed":config.seed,"level":"example-36x72","label":"generated-land-example","physical_keying":true,
        "chart":{"proposed_chart_centre_latlon_deg":[lat,lon],"site_in_chart_mm":[0,0],
        "walk_tile":{"lod":0,"ix":0,"iy":0},"region":{"lod":6,"ix":0,"iy":0}},
        "site":{"site_point_latlon_deg":[lat,lon],"t0_cell":cell},
        "reach":{"reach_id":out.cell_reach[cell].map(|r|r.0).unwrap_or(0)},
        "recipe":{"universe_digest":universe_digest,"artifact_digests":{"terrain":header["digest"],"hydrology":hydro["digest"]}}});
    std::fs::write(root.join("site.json"),serde_json::to_vec_pretty(&site)?)?;
    let report=json!({"schema":"maiforge.example-inputs/v1","grid":[36,72],"seed":config.seed,
        "duration_my":4.0,"step_my":2.0,"profile":"world_strata_v2","cell":cell,"position":[lat,lon],
        "full_column_depth_bound_m":depth,"layers":full.lith.len(),"digests":digests,
        "seconds":started.elapsed().as_secs_f64(),"limit":"Coarse demonstration recipe, not calibrated terrain or a gameplay acceptance."});
    std::fs::write(root.join("generation.json"),serde_json::to_vec_pretty(&report)?)?;
    println!("{}",serde_json::to_string_pretty(&report)?); Ok(())
}
