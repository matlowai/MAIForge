//! CLIMATE-V0 (lane CL, 2026-10-12): IGNORED evidence driver. Run one
//! process per capture (memory-capped by the caller):
//!
//! ```text
//! CL_RUN=<verified run root> CL_OUT=<new dir> \
//!   climate_v0_capture --ignored --exact climate_v0_capture_run
//! # optional: simulate + publish the source run first
//! CL_SIM=1 CL_SEED=7 CL_FIDELITY=standard CL_MODE=bodies CL_RUN=<new dir> CL_OUT=<new dir> ...
//! ```
//!
//! Writes `<CL_OUT>/climate` (the published `pf.climate/v0` artifact) and
//! `<CL_OUT>/summary.json` (digests, summary, closures, zonal rows).

use pf_contracts::laws::ClimateRecipeV1;
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_lab::preview::{publish, publish_climate_v0, simulate_with, PreviewConfig, PreviewOptions};
use std::path::PathBuf;
use std::time::Instant;

fn env(k: &str) -> Option<String> {
    std::env::var(k).ok().filter(|v| !v.is_empty())
}

#[test]
#[ignore]
fn climate_v0_capture_run() {
    let run = PathBuf::from(env("CL_RUN").expect("CL_RUN"));
    let out = PathBuf::from(env("CL_OUT").expect("CL_OUT"));
    let started = Instant::now();
    let mut sim = serde_json::Value::Null;
    if env("CL_SIM").as_deref() == Some("1") {
        let seed: u64 = env("CL_SEED").expect("CL_SEED").parse().unwrap();
        let fidelity = match env("CL_FIDELITY").as_deref() {
            Some("preview") => FidelityPreset::Preview,
            Some("standard") => FidelityPreset::Standard,
            other => panic!("CL_FIDELITY {other:?}"),
        };
        let opts = match env("CL_MODE").as_deref() {
            Some("legacy") => PreviewOptions::legacy(),
            Some("bodies") => PreviewOptions::crust_bodies_v0(),
            Some("two_column") => PreviewOptions::two_column_v0(),
            other => panic!("CL_MODE {other:?}"),
        };
        let cfg = PreviewConfig {
            scenario: ScenarioPreset::toy(),
            fidelity,
            grid: None,
            seed,
            duration_my: None,
            step_my: None,
        };
        let o = simulate_with(&cfg, &opts).unwrap();
        std::fs::create_dir_all(&run).unwrap();
        let d = publish(&o, &run, started).unwrap();
        sim = serde_json::json!({
            "seed": seed, "mode": env("CL_MODE"), "fidelity": env("CL_FIDELITY"),
            "terrain": d.terrain, "hydrology": d.hydrology, "resource": d.resource,
            "simulate_seconds": started.elapsed().as_secs_f64(),
        });
    }
    std::fs::create_dir_all(&out).unwrap();
    let t0 = Instant::now();
    let recipe = ClimateRecipeV1::coarse_v0_earth();
    let p = publish_climate_v0(&run, &out, &recipe, t0).unwrap();
    let summary = serde_json::json!({
        "source_run": run,
        "simulated": sim,
        "climate_digest": p.digest,
        "terrain_digest": p.terrain_digest,
        "hydrology_digest": p.hydrology_digest,
        "recipe_digest": recipe.digest().unwrap(),
        "climate_seconds": t0.elapsed().as_secs_f64(),
        "summary": p.output.summary,
        "closures": p.output.closures,
        "rows": p.output.rows,
    });
    std::fs::write(
        out.join("summary.json"),
        serde_json::to_string_pretty(&summary).unwrap(),
    )
    .unwrap();
    println!("climate {} in {:.2}s", p.digest, t0.elapsed().as_secs_f64());
}
