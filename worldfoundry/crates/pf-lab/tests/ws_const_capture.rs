//! WS-CONST (lane WC, 2026-10-12): capture drivers for the before/after
//! byte-identity proof. Both tests are `#[ignore]`d measurement drivers that
//! use only APIs that existed before WS-CONST, so the SAME file compiles and
//! runs against the pre-change and post-change source. Lane WC runs them
//! explicitly (env-driven) and compares the outputs byte for byte.
//!
//! - `ws_const_capture_profiles`: every scenario preset x fidelity x
//!   physics-option profile, as the canonical JSON the params digest is taken
//!   over, plus that digest, written to `$WS_OUT/profiles.json`.
//! - `ws_const_capture_flag_run`: one `simulate_with` run in a flag mode
//!   (`WS_MODE` = two_column | bodies | legacy), `WS_SEED`, `WS_FIDELITY`
//!   (preview | standard), published to `$WS_OUT` with the step log, the
//!   t = 0 elevation, the t = 0 sea-level solve and the resolved params.

use pf_contracts::artifact::digest_params;
use pf_contracts::params::{FidelityPreset, MantleBlock, ScenarioPreset};
use pf_lab::preview::{publish, simulate_with, PreviewConfig, PreviewOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

fn env(k: &str) -> String {
    std::env::var(k).unwrap_or_else(|_| panic!("set {k}"))
}

fn scenarios() -> Vec<ScenarioPreset> {
    vec![
        ScenarioPreset::toy(),
        ScenarioPreset::benchmark_a(),
        ScenarioPreset::benchmark_b(),
        ScenarioPreset::benchmark_c(),
    ]
}

const FIDELITIES: [FidelityPreset; 3] = [
    FidelityPreset::Preview,
    FidelityPreset::Standard,
    FidelityPreset::High,
];

fn options() -> Vec<(&'static str, PreviewOptions)> {
    vec![
        ("legacy", PreviewOptions::legacy()),
        ("two_column", PreviewOptions::two_column_v0()),
        ("bodies", PreviewOptions::crust_bodies_v0()),
    ]
}

#[test]
#[ignore]
fn ws_const_capture_profiles() {
    let out = PathBuf::from(env("WS_OUT"));
    std::fs::create_dir_all(&out).unwrap();
    let mut rows = Vec::new();
    for s in scenarios() {
        for f in FIDELITIES {
            // Contract-level profile (MantleBlock default: Flexure).
            let p = s.profile(f, MantleBlock::default());
            rows.push(serde_json::json!({
                "name": format!("contract/{}/{:?}", s.name, f),
                "json": serde_json::to_string(&p).unwrap(),
                "digest": digest_params(&p).unwrap(),
            }));
            // pf-lab resolved profiles (Airy forced) per physics option.
            for (tag, o) in options() {
                let cfg = PreviewConfig {
                    scenario: s.clone(),
                    fidelity: f,
                    grid: None,
                    seed: 1,
                    duration_my: None,
                    step_my: None,
                };
                let p = cfg.profile_with(&o).unwrap();
                rows.push(serde_json::json!({
                    "name": format!("preview/{}/{:?}/{tag}", s.name, f),
                    "json": serde_json::to_string(&p).unwrap(),
                    "digest": digest_params(&p).unwrap(),
                }));
            }
        }
    }
    let mut file = std::fs::File::create(out.join("profiles.json")).unwrap();
    file.write_all(serde_json::to_string_pretty(&rows).unwrap().as_bytes())
        .unwrap();
}

#[test]
#[ignore]
fn ws_const_capture_flag_run() {
    let out = PathBuf::from(env("WS_OUT"));
    let seed: u64 = env("WS_SEED").parse().unwrap();
    let fidelity = match env("WS_FIDELITY").as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("fidelity {other}"),
    };
    let opts = match env("WS_MODE").as_str() {
        "legacy" => PreviewOptions::legacy(),
        "two_column" => PreviewOptions::two_column_v0(),
        "bodies" => PreviewOptions::crust_bodies_v0(),
        other => panic!("mode {other}"),
    };
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: None,
        step_my: None,
    };
    let started = Instant::now();
    let o = simulate_with(&cfg, &opts).unwrap();
    let digests = publish(&o, &out, started).unwrap();
    let extra = |name: &str, bytes: Vec<u8>| {
        std::fs::write(out.join(name), bytes).unwrap();
    };
    extra("steps.json", serde_json::to_vec(&o.steps).unwrap());
    extra("params.json", serde_json::to_vec(&o.params).unwrap());
    extra(
        "sea_level_t0.json",
        serde_json::to_vec(&o.sea_level_t0).unwrap(),
    );
    extra(
        "elevation_t0_m.f64",
        o.elevation_t0_m
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect(),
    );
    if let Some(r) = &o.continental_bodies {
        extra(
            "bodies_coast_distance_m.f64",
            r.coast_distance_m
                .iter()
                .flat_map(|v| v.to_le_bytes())
                .collect(),
        );
        extra(
            "bodies_summary.txt",
            format!(
                "{:?}\n{:?}\n{:?}\n",
                r.bodies, r.achieved_area_fraction, r.continental_crust_volume_m3
            )
            .into_bytes(),
        );
    }
    extra("digests.json", serde_json::to_vec(&digests).unwrap());
}
