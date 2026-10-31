//! Build the T1a skeleton for a picked site and write it with its
//! falsifier report. Labelled pre-signoff.
//!
//! Usage: pf-region-skeleton <t0_run_dir> <subgrid_relief.f32> <site_record.json> <out_dir> [amplitude_permille]
//!        [climate=<climate_artifact_dir> [climate_floor_permille=N]]
//!   climate= (lane CLIMATE-RUNOFF, opt-in): discharge carries the seed-climate-v1 runoff field (pf_region::climate_runoff);
//!   without it every byte is unchanged.
//!   Lane T0-IMPRINT-PROD (opt-in, ADR 045; absent = every byte unchanged):
//!     junction=0|1|2      params.imprint: junction mode (0 cell centre, 1 keyed offset, 2 C1 low ground + keyed jitter); F3 per T0 outlet
//!     junction_offset_permille=N, junction_jitter_permille=N   (after junction=; defaults 600 and 0 / 150)
//!     corridor_mm=N, control_spacing_mm=N                       planform corridor half width and control spacing (existing params)

use pf_region::skeleton::{build_with_twin, Inputs};
use pf_region::t0::{SiteRecord, T0};
use pf_region::SkeletonParamsV1;
use std::path::PathBuf;
use std::time::Instant;

fn peak_rss_kib() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    s.lines().find(|l| l.starts_with("VmHWM:"))?.split_whitespace().nth(1)?.parse().ok()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 5 {
        eprintln!("usage: {} <t0_run_dir> <subgrid_relief.f32> <site_record.json> <out_dir> [amplitude_permille]", a[0]);
        std::process::exit(2);
    }
    let (run, sg, site, out) = (PathBuf::from(&a[1]), PathBuf::from(&a[2]), PathBuf::from(&a[3]), PathBuf::from(&a[4]));
    let mut params = SkeletonParamsV1::default();
    let mut climate: Option<PathBuf> = None;
    let mut floor: Option<i64> = None;
    for s in &a[5..] {
        if let Some(d) = s.strip_prefix("climate=") { climate = Some(PathBuf::from(d)); }
        else if let Some(v) = s.strip_prefix("climate_floor_permille=") { floor = Some(v.parse().expect("climate_floor_permille")); }
        else if let Some(v) = s.strip_prefix("junction=") { params.imprint = Some(pf_region::skeleton::ImprintParamsV1::v1(v.parse().expect("junction mode"))); }
        else if let Some(v) = s.strip_prefix("junction_offset_permille=") { params.imprint.as_mut().expect("junction_offset_permille needs junction= before it").junction_offset_permille = v.parse().expect("permille"); }
        else if let Some(v) = s.strip_prefix("junction_jitter_permille=") { params.imprint.as_mut().expect("junction_jitter_permille needs junction= before it").junction_jitter_permille = v.parse().expect("permille"); }
        else if let Some(v) = s.strip_prefix("corridor_mm=") { params.corridor_half_width_mm = v.parse().expect("corridor_mm"); }
        else if let Some(v) = s.strip_prefix("control_spacing_mm=") { params.planform_control_spacing_mm = v.parse().expect("control_spacing_mm"); }
        else { params.planform_amplitude_permille = s.parse().expect("amplitude permille"); }
    }
    if floor.is_some() && climate.is_none() { eprintln!("climate_floor_permille needs climate=<dir>"); std::process::exit(2); }
    let t_all = Instant::now();
    let site = SiteRecord::load(&site).unwrap_or_else(|e| { eprintln!("site: {e}"); std::process::exit(1) });
    let t_load = Instant::now();
    let mut t0 = T0::load(&run, &sg, &site, params.threshold_m3_s()).unwrap_or_else(|e| { eprintln!("t0: {e}"); std::process::exit(1) });
    if let Some(d) = &climate {
        let mut cp = pf_region::climate_runoff::ClimateRunoffParamsV1::v1();
        if let Some(f) = floor { cp.floor_permille = f; }
        t0.attach_climate(d, cp.clone()).unwrap_or_else(|e| { eprintln!("climate: {e}"); std::process::exit(1) });
        params.runoff_rule = format!("climate_runoff_v1: R(x) x node area (params.climate_runoff; upstream t0.climate); r_ref = A-RAIN-UNITS {} m/yr", t0.runoff_m_per_yr);
        params.climate_runoff = Some(cp);
    }
    let load_s = t_load.elapsed().as_secs_f64();
    let t_build = Instant::now();
    let sk = build_with_twin(&Inputs { t0: &t0, site: &site, params }).unwrap_or_else(|e| { eprintln!("build: {e}"); std::process::exit(1) });
    let build_s = t_build.elapsed().as_secs_f64();
    let t_check = Instant::now();
    let checks = pf_region::checks::run(&sk, &t0).unwrap_or_else(|e| { eprintln!("checks: {e}"); std::process::exit(1) });
    let check_s = t_check.elapsed().as_secs_f64();
    let t_pub = Instant::now();
    let p = pf_region::publish::publish(&sk, &out).unwrap_or_else(|e| { eprintln!("publish: {e}"); std::process::exit(1) });
    let pub_s = t_pub.elapsed().as_secs_f64();
    std::fs::write(out.join("checks.json"), serde_json::to_vec_pretty(&checks).unwrap()).unwrap();
    let exec = serde_json::json!({
        "load_t0_verify_seconds": load_s, "build_seconds": build_s, "checks_seconds": check_s, "publish_seconds": pub_s,
        "total_seconds": t_all.elapsed().as_secs_f64(), "peak_rss_kib": peak_rss_kib(),
        "not_digested": "execution receipt; not part of the skeleton digest",
    });
    std::fs::write(out.join("execution.json"), serde_json::to_vec_pretty(&exec).unwrap()).unwrap();
    println!("skeleton {} recipe {}", p.digest, sk.recipe_digest);
    println!("segments t0={} head={} terminal={} retried={} straight={} forced={}", sk.stats.segments_t0_flow, sk.stats.segments_head,
        sk.stats.segments_terminal, sk.stats.segments_retried, sk.stats.segments_straight_fallback, sk.stats.segments_forced);
    if let Some(f) = &checks.f3_outlet {
        println!("F3 per T0 outlet (ADR 045) pass {} (per-cell diagnostic pass {:?}) | cells {} | basins {}", f.pass, checks.f3_cell_pass,
            f.cells.iter().map(|c| format!("{}:{:.0}", c.cell, c.share_permille)).collect::<Vec<_>>().join(" "),
            f.basins.iter().map(|b| format!("{}:{:+.0}", b.basin, b.rel_err_permille)).collect::<Vec<_>>().join(" "));
    }
    println!("F1 {} ({} reaches)  F2 {} (strict {} lenient {} nonriver {})  F3 {} (min gated {:.1} permille)  B4 {} ({} nodes, {:.0} m in tile, min Q {:.1})  corridor {} (max {:.0} m)",
        checks.f1_pass, checks.f1_reaches_checked, checks.f2_pass, checks.f2_strict_hillslope_nodes_ge_threshold,
        checks.f2_lenient_leaving_own_cell, checks.f2_axis_ge_threshold_on_non_river_cells, checks.f3_pass, checks.f3_min_gated_share_permille,
        checks.b4_pass, checks.b4_channel_nodes_in_tile, checks.b4_channel_length_in_tile_m, checks.b4_min_q_in_tile_m3_s,
        checks.corridor_pass, checks.corridor_max_dist_m);
    if let Some(r) = &sk.climate_runoff_m_per_yr {
        let c = t0.climate.as_ref().unwrap();
        let land: Vec<f64> = (0..r.len()).filter(|&n| sk.kind[n] != 2 && sk.kind[n] != 3).map(|n| r[n]).collect();
        let mut s = land.clone(); s.sort_by(|a, b| a.total_cmp(b));
        let pct = |q: f64| s.get(((s.len() as f64 - 1.0) * q).round() as usize).copied().unwrap_or(f64::NAN);
        println!("climate {} tier {} r_ref {} m/yr: chart land R min {:.3} p10 {:.3} p50 {:.3} p90 {:.3} max {:.3} m/yr; max Q {:.1} m3/s; twin {}",
            c.digest, c.tier, c.r_ref_m_per_yr, pct(0.0), pct(0.1), pct(0.5), pct(0.9), pct(1.0),
            sk.discharge_m3_s.iter().cloned().fold(0.0, f64::max), sk.climate_twin_digest.as_deref().unwrap_or("-"));
        for l in &sk.climate_lakes { println!("  lake {}", serde_json::to_string(l).unwrap()); }
    }
    println!("{}", exec);
}
