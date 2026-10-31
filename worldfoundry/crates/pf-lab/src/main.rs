//! pf-lab: the Planet Foundry driver.
//!
//! ```text
//! pf-lab simulate --preset toy --seed 240907 --output DIR
//!     [--fidelity preview|standard|high] [--mini NLATxNLON]
//!     [--duration-my X] [--step-my Y]
//! pf-lab extract-region --run DIR --lat L --lon L --half-extent-mm N
//!     --spacing-mm R --output DIR [--id region]
//!     (metre flags --half-extent-m / --resolution-m are accepted and
//!     converted; the frame is integer millimetres per ADR 029 section 3)
//! ```
//!
//! No argument parser dependency (offline build); errors are loud and
//! name the flag.

use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_lab::preview::{publish, simulate, PreviewConfig};
use pf_lab::region::{extract, RegionRequest, SourceRun};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

fn usage() -> ! {
    eprintln!(
        "usage:\n  pf-lab simulate --preset <toy|benchmark-a|benchmark-b|benchmark-c> --seed <u64> --output <dir> [--fidelity preview|standard|high] [--mini NLATxNLON] [--duration-my X] [--step-my Y]\n  pf-lab extract-region --run <dir> --lat <deg> --lon <deg> --half-extent-mm <i64> --spacing-mm <u32> --output <dir> [--id <artifact id>]  (or --half-extent-m / --resolution-m, converted x1000)"
    );
    std::process::exit(2);
}

fn parse_flags(args: &[String]) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if let Some(key) = a.strip_prefix("--") {
            let Some(v) = args.get(i + 1) else {
                eprintln!("flag --{key} needs a value");
                usage();
            };
            out.insert(key.to_string(), v.clone());
            i += 2;
        } else {
            eprintln!("unexpected argument {a:?}");
            usage();
        }
    }
    out
}

fn required<'a>(flags: &'a BTreeMap<String, String>, key: &str) -> &'a str {
    match flags.get(key) {
        Some(v) => v,
        None => {
            eprintln!("missing --{key}");
            usage();
        }
    }
}

fn parse_num<T: std::str::FromStr>(flags: &BTreeMap<String, String>, key: &str) -> Option<T> {
    flags.get(key).map(|v| {
        v.parse::<T>().unwrap_or_else(|_| {
            eprintln!("--{key}: cannot parse {v:?}");
            usage();
        })
    })
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(cmd) = args.first() else { usage() };
    let flags = parse_flags(&args[1..]);
    let started = Instant::now();
    match cmd.as_str() {
        "simulate" => {
            let scenario = match required(&flags, "preset") {
                "toy" => ScenarioPreset::toy(),
                "benchmark-a" => ScenarioPreset::benchmark_a(),
                "benchmark-b" => ScenarioPreset::benchmark_b(),
                "benchmark-c" => ScenarioPreset::benchmark_c(),
                other => {
                    eprintln!("unknown preset {other:?}");
                    usage();
                }
            };
            let fidelity = match flags
                .get("fidelity")
                .map(String::as_str)
                .unwrap_or("preview")
            {
                "preview" => FidelityPreset::Preview,
                "standard" => FidelityPreset::Standard,
                "high" => FidelityPreset::High,
                other => {
                    eprintln!("unknown fidelity {other:?}");
                    usage();
                }
            };
            let grid = flags.get("mini").map(|s| {
                let (a, b) = s.split_once('x').unwrap_or_else(|| {
                    eprintln!("--mini expects NLATxNLON");
                    usage();
                });
                (
                    a.parse::<usize>().unwrap_or_else(|_| usage()),
                    b.parse::<usize>().unwrap_or_else(|_| usage()),
                )
            });
            let seed: u64 = parse_num(&flags, "seed").unwrap_or_else(|| {
                eprintln!("missing --seed");
                usage()
            });
            let output = PathBuf::from(required(&flags, "output"));
            let cfg = PreviewConfig {
                scenario,
                fidelity,
                grid,
                seed,
                duration_my: parse_num(&flags, "duration-my"),
                step_my: parse_num(&flags, "step-my"),
            };
            let out = match simulate(&cfg) {
                Ok(o) => o,
                Err(e) => {
                    eprintln!("simulate failed: {e}");
                    std::process::exit(1);
                }
            };
            let sim_seconds = started.elapsed().as_secs_f64();
            if let Err(e) = std::fs::create_dir_all(&output) {
                eprintln!("cannot create {}: {e}", output.display());
                std::process::exit(1);
            }
            let digests = match publish(&out, &output, started) {
                Ok(d) => d,
                Err(e) => {
                    eprintln!("publish failed: {e}");
                    std::process::exit(1);
                }
            };
            let summary = serde_json::json!({
                "backend": {"name": pf_lab::preview::BACKEND_NAME, "tier": "Preview", "version": pf_lab::preview::BACKEND_VERSION},
                "grid": {"nlat": out.grid.nlat(), "nlon": out.grid.nlon(), "cells": out.grid.cells(), "ladder_level": out.grid.ladder_level()},
                "seed_path": out.seed_path,
                "steps": out.terrain.summary.steps,
                "step_my": out.terrain.summary.step_my,
                "t_my": out.terrain.t_my,
                "digests": digests,
                "summary": out.terrain.summary,
                "crust_volume": out.terrain.crust_volume,
                "crust_volume_residual_fraction": out.terrain.crust_volume.residual_fraction(),
                "sediment_cumulative_kg": {
                    "eroded": out.terrain.sediment_balance.cumulative_eroded_kg,
                    "deposited": out.terrain.sediment_balance.cumulative_deposited_kg,
                    "exported": out.terrain.sediment_balance.cumulative_exported_kg,
                    "declared_error": out.terrain.sediment_balance.cumulative_declared_error_kg,
                },
                "erosion": out.terrain.erosion,
                "hydrology": {"basins": out.hydrology.basins.len(), "reaches": out.hydrology.rivers.len(), "waterfalls": out.hydrology.waterfalls.len()},
                "boundaries": out.terrain.boundaries.len(),
                "resource_mean_potential": out.resource.families.iter().map(|f| (format!("{:?}", f.kind), f.mean_potential)).collect::<BTreeMap<_, _>>(),
                "timing": {"simulate_seconds": sim_seconds, "total_seconds": started.elapsed().as_secs_f64(), "peak_rss_kib": pf_lab::receipt::peak_rss_kib()},
                "step_log_tail": out.steps.iter().rev().take(3).collect::<Vec<_>>(),
            });
            println!("{}", serde_json::to_string_pretty(&summary).expect("json"));
        }
        "extract-region" => {
            let run = PathBuf::from(required(&flags, "run"));
            let half_extent_mm: i64 = parse_num::<i64>(&flags, "half-extent-mm")
                .or_else(|| {
                    parse_num::<i64>(&flags, "half-extent-m").map(|m| {
                        m.checked_mul(1000).unwrap_or_else(|| {
                            eprintln!("metre to millimetre overflow");
                            usage()
                        })
                    })
                })
                .unwrap_or_else(|| {
                    eprintln!("missing --half-extent-mm (or --half-extent-m)");
                    usage()
                });
            let spacing_mm: u32 = parse_num::<u32>(&flags, "spacing-mm")
                .or_else(|| {
                    parse_num::<u32>(&flags, "resolution-m").map(|m| {
                        m.checked_mul(1000).unwrap_or_else(|| {
                            eprintln!("metre to millimetre overflow");
                            usage()
                        })
                    })
                })
                .unwrap_or_else(|| {
                    eprintln!("missing --spacing-mm (or --resolution-m)");
                    usage()
                });
            let req = RegionRequest {
                lat_deg: parse_num(&flags, "lat").unwrap_or_else(|| usage()),
                lon_deg: parse_num(&flags, "lon").unwrap_or_else(|| usage()),
                half_extent_mm,
                spacing_mm,
            };
            let output = PathBuf::from(required(&flags, "output"));
            let id = flags.get("id").cloned().unwrap_or_else(|| "region".into());
            let source = match SourceRun::load(&run) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("cannot load run {}: {e}", run.display());
                    std::process::exit(1);
                }
            };
            let out = match extract(&source, &req) {
                Ok(o) => o,
                Err(e) => {
                    eprintln!("extract failed: {e}");
                    std::process::exit(1);
                }
            };
            if let Err(e) = std::fs::create_dir_all(&output) {
                eprintln!("cannot create {}: {e}", output.display());
                std::process::exit(1);
            }
            let digest = match pf_lab::region::publish(&out, &source, &output, &id, started) {
                Ok(d) => d,
                Err(e) => {
                    eprintln!("publish failed: {e}");
                    std::process::exit(1);
                }
            };
            if let Some(tile_path) = flags.get("tile-output") {
                let verified = pf_contracts::region::read_region_v1_fields(&output.join(&id))
                    .unwrap_or_else(|e| {
                        eprintln!("region verify failed: {e}");
                        std::process::exit(1)
                    });
                assert_eq!(verified.artifact().digest, digest, "published region digest");
                let tile = pf_lab::region::producer_tile(&verified, &source, &id)
                    .unwrap_or_else(|e| {
                        eprintln!("tile export failed: {e}");
                        std::process::exit(1)
                    });
                let file = std::fs::OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(tile_path)
                    .unwrap_or_else(|e| {
                        eprintln!("tile output {tile_path}: {e}");
                        std::process::exit(1)
                    });
                serde_json::to_writer(file, &tile).unwrap_or_else(|e| {
                    eprintln!("tile write failed: {e}");
                    std::process::exit(1)
                });
            }
            let summary = serde_json::json!({
                "artifact_id": id,
                "digest": digest,
                "frame": out.domain.frame,
                "source": out.domain.source,
                "resampling": out.domain.resampling,
                "reaches_in_region": out.domain.reaches.len(),
                "lakes_in_region": out.domain.water_bodies.iter().filter(|b|b.kind==pf_contracts::domain::WaterBodyKindV0::Lake).count(),
                "availability": out.domain.availability,
                "state": out.domain.state,
                "timing": {"total_seconds": started.elapsed().as_secs_f64(), "peak_rss_kib": pf_lab::receipt::peak_rss_kib()},
            });
            println!("{}", serde_json::to_string_pretty(&summary).expect("json"));
        }
        _ => usage(),
    }
}
