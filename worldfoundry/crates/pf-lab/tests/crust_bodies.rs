//! CRUST-SPIKE (lane K2, 2026-10-11): the preview chain with the opt-in
//! two-column isostasy wired into `simulate_with` and the opt-in
//! continental crust bodies.
//!
//! - `simulate()` equals `simulate_with(legacy options)` bitwise;
//! - flags on: three runs give identical state and published digests;
//! - two-column mode closes the ocean water volume at every step and logs
//!   the solved sea level; the flags move params_digest and declare their
//!   assumptions and approximations in the terrain artifact;
//! - the ignored `measure_crust` writes raw fields for the lane K2
//!   measurement (env-driven: CRUST_SEED, CRUST_FIDELITY, CRUST_MODE,
//!   CRUST_OUT). It is not part of the normal suite.

use pf_contracts::artifact::{read_artifact_header, verify_artifact};
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::terrain::{CrustKindCodeV0, TerrainDomainV0};
use pf_geodynamics::kinematic::CrustState;
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::uplift::UpliftPrior;
use pf_lab::preview::{
    elevation_on_grid, publish, simulate, simulate_with, PreviewConfig, PreviewOptions,
};
use pf_surface::water::largest_ocean;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn mini(seed: u64) -> PreviewConfig {
    PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed,
        duration_my: Some(20.0),
        step_my: Some(2.0),
    }
}

/// Persistent scratch root: std::env::temp_dir() honours TMPDIR, which the
/// lane K2 runner points at `generated-test-data/tmp` (never /tmp).
fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("pf-lab-k2-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn simulate_equals_simulate_with_legacy_options_bitwise() {
    for seed in [240907u64, 7] {
        let cfg = mini(seed);
        let a = simulate(&cfg).unwrap();
        let b = simulate_with(&cfg, &PreviewOptions::legacy()).unwrap();
        assert_eq!(a.elevation_m, b.elevation_m);
        assert_eq!(a.state, b.state);
        assert_eq!(a.terrain, b.terrain);
        assert_eq!(a.params, b.params);
        assert!(a.sea_level_t0.is_none() && a.continental_bodies.is_none());
        assert!(a.steps.iter().all(|s| s.sea_level_gauge_m.is_none()));
        // The step log serializes exactly as before (no new keys).
        let j = serde_json::to_string(&a.steps).unwrap();
        assert!(!j.contains("sea_level_gauge_m") && !j.contains("flooded"));
    }
}

#[test]
fn flags_on_three_runs_are_byte_identical() {
    for opts in [PreviewOptions::two_column_v0(), PreviewOptions::crust_bodies_v0()] {
        let cfg = mini(240907);
        let mut digests = Vec::new();
        let mut first: Option<pf_lab::preview::PreviewOutputs> = None;
        for k in 0..3 {
            let out = simulate_with(&cfg, &opts).unwrap();
            let root = scratch(&format!("det{k}"));
            digests.push(publish(&out, &root, Instant::now()).unwrap());
            if let Some(f) = &first {
                assert_eq!(f.elevation_m, out.elevation_m);
                assert_eq!(f.state, out.state);
                assert_eq!(f.terrain, out.terrain);
                assert_eq!(f.elevation_t0_m, out.elevation_t0_m);
                // NaN marks non-body cells in coast_distance_m, so the
                // report is compared field by field with bit patterns.
                match (&f.continental_bodies, &out.continental_bodies) {
                    (None, None) => {}
                    (Some(x), Some(y)) => {
                        assert_eq!(x.bodies, y.bodies);
                        assert_eq!(x.achieved_area_fraction, y.achieved_area_fraction);
                        assert_eq!(x.continental_crust_volume_m3, y.continental_crust_volume_m3);
                        let bits = |v: &[f64]| v.iter().map(|d| d.to_bits()).collect::<Vec<_>>();
                        assert_eq!(bits(&x.coast_distance_m), bits(&y.coast_distance_m));
                    }
                    _ => panic!("bodies report present in one run only"),
                }
                let sa = serde_json::to_string(&f.steps).unwrap();
                let sb = serde_json::to_string(&out.steps).unwrap();
                assert_eq!(sa, sb);
            } else {
                first = Some(out);
            }
            let _ = std::fs::remove_dir_all(&root);
        }
        assert_eq!(digests[0], digests[1]);
        assert_eq!(digests[1], digests[2]);
        let legacy = simulate(&cfg).unwrap();
        let root = scratch("legacy");
        let dl = publish(&legacy, &root, Instant::now()).unwrap();
        let _ = std::fs::remove_dir_all(&root);
        assert_ne!(dl.terrain, digests[0].terrain, "flag had no effect");
    }
}

#[test]
fn two_column_wiring_closes_water_and_declares_itself() {
    let cfg = mini(7);
    let out = simulate_with(&cfg, &PreviewOptions::crust_bodies_v0()).unwrap();
    let v = out.params.two_column_isostasy.as_ref().unwrap().ocean_water_volume_m3;
    let t0 = out.sea_level_t0.unwrap();
    assert!(t0.ocean_volume_residual_m3.abs() <= 1e-9 * v);
    assert_eq!(out.steps.len(), 10);
    for s in &out.steps {
        let r = s.ocean_volume_residual_m3.expect("logged per step");
        assert!(r.abs() <= 1e-9 * v, "step {} residual {r}", s.t_my);
        assert!(s.sea_level_gauge_m.unwrap().is_finite());
        let d = s.compensation_depth_m.unwrap();
        assert!((d - (s.deepest_root_m.unwrap() + 10_000.0)).abs() < 1e-9, "derived depth");
        let f = s.flooded_area_fraction.unwrap();
        assert!((0.0..=1.0).contains(&f));
    }
    let ids: Vec<&str> = out.terrain.assumptions.iter().map(|a| a.id.as_str()).collect();
    for id in [
        "A-OCEAN-VOLUME",
        "A-OCEANIC-CRUST-DENSITY",
        "A-COMPENSATION-DEPTH-MARGIN",
        "A-SEA-LEVEL-VOLUME-PRIMARY",
        "A-CRUST-BODY-AREA-FRACTION",
        "A-RIFTED-MARGIN-WIDTH",
    ] {
        assert!(ids.contains(&id), "missing {id}");
    }
    for gone in ["A-SEA-LEVEL-CONSTANT", "A-GDH1-REFERENCE-OFFSET"] {
        assert!(!ids.contains(&gone), "{gone} is false under two-column");
    }
    let ap = &out.terrain.approximations;
    assert!(ap.iter().any(|a| a.starts_with("two-column")));
    assert!(!ap.iter().any(|a| a == "airy-crust-only-isostasy-with-water-loading"));
    assert!(ap.iter().any(|a| a == "linear-rifted-margin-taper-v0"));
    // Published, verified, and the params digest differs from legacy.
    let root = scratch("verify");
    let d = publish(&out, &root, Instant::now()).unwrap();
    assert_eq!(
        verify_artifact::<TerrainDomainV0>(&root.join("terrain")).unwrap(),
        d.terrain
    );
    let h_on = read_artifact_header::<TerrainDomainV0>(&root.join("terrain")).unwrap();
    let legacy = simulate(&cfg).unwrap();
    let root2 = scratch("verify-legacy");
    publish(&legacy, &root2, Instant::now()).unwrap();
    let h_off = read_artifact_header::<TerrainDomainV0>(&root2.join("terrain")).unwrap();
    assert_ne!(
        h_on.header.semantic.params_digest,
        h_off.header.semantic.params_digest
    );
    let _ = std::fs::remove_dir_all(&root);
    let _ = std::fs::remove_dir_all(&root2);
    // Body bookkeeping on the real chain: coverage near the derived 0.38.
    let rep = out.continental_bodies.as_ref().unwrap();
    assert!((rep.declared_area_fraction - 0.38).abs() < 1e-12);
    assert!((rep.achieved_area_fraction - 0.38).abs() < 0.02, "{}", rep.achieved_area_fraction);
}

#[test]
fn derived_compensation_depth_known_answer() {
    // One column with a root deeper than the fixed 100 km: 70 km
    // continental crust under 60 km of sediment, both at 2800 kg/m3, no
    // thermal term. Hand value: root depth = 2800 * 130000 / 3300
    // = 110303.0303 m; derived depth = 110303.0303 + 10000 (margin).
    let cfg = mini(240907);
    let mut fixed = PreviewOptions::crust_bodies_v0();
    fixed.two_column_isostasy.as_mut().unwrap().compensation_depth_margin_m = None;
    let p_fixed = cfg.profile_with(&fixed).unwrap();
    let p_derived = cfg.profile_with(&PreviewOptions::crust_bodies_v0()).unwrap();
    let grid = cfg.grid(&p_derived).unwrap();
    let seeding = seed_plates(&grid, &p_derived, &cfg.seed_path()).unwrap();
    let (mut state, _) =
        CrustState::seed_bodies(&grid, &seeding, &UpliftPrior::preview_v0(), 0.38, 200_000.0)
            .unwrap();
    let mut areas = Vec::new();
    for i in 0..grid.nlat() {
        let a = grid.cell_area_m2(i).unwrap();
        for _ in 0..grid.nlon() {
            areas.push(a);
        }
    }
    // Before the deep column: fixed 100 km and derived depth give the
    // same elevations up to rounding (the depth cancels).
    let (e_fixed, s_fixed) = elevation_on_grid(&state, &p_fixed, &areas).unwrap();
    let (e_der, s_der) = elevation_on_grid(&state, &p_derived, &areas).unwrap();
    let max_diff = e_fixed
        .iter()
        .zip(&e_der)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(max_diff < 1e-6, "elevations moved with the depth: {max_diff} m");
    assert_eq!(s_fixed.unwrap().compensation_depth_m, 100_000.0);
    let sd = s_der.unwrap();
    // Deepest root at t0 is a 35 km continental interior: 2800*35000/3300.
    assert!((sd.deepest_root_m - 2800.0 * 35_000.0 / 3300.0).abs() < 1e-6);
    assert!((sd.compensation_depth_m - (sd.deepest_root_m + 10_000.0)).abs() < 1e-9);
    // Add the deep column: fixed mode refuses loudly, derived mode solves.
    let c = state.kind.iter().position(|k| *k == CrustKindCodeV0::Continental).unwrap();
    state.thickness_m[c] = 70_000.0;
    state.sediment_m[c] = 60_000.0;
    let err = elevation_on_grid(&state, &p_fixed, &areas).unwrap_err().to_string();
    assert!(err.contains("compensation depth"), "{err}");
    let (e2, s2) = elevation_on_grid(&state, &p_derived, &areas).unwrap();
    let s2 = s2.unwrap();
    let hand_root: f64 = 2800.0 * 130_000.0 / 3300.0;
    assert!((hand_root - 110_303.030_303).abs() < 1e-3);
    assert!((s2.deepest_root_m - hand_root).abs() < 1e-6, "{}", s2.deepest_root_m);
    assert!((s2.compensation_depth_m - (hand_root + 10_000.0)).abs() < 1e-6);
    assert_eq!(s2.deepest_root_crust_sediment_m, [70_000.0, 60_000.0]);
    // Same field solved with a hand-set fixed depth of 200 km agrees.
    let mut p200 = p_fixed.clone();
    p200.two_column_isostasy.as_mut().unwrap().compensation_depth_m = 200_000.0;
    let (e200, _) = elevation_on_grid(&state, &p200, &areas).unwrap();
    let d = e200.iter().zip(&e2).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    assert!(d < 1e-6, "{d}");
    // Water still closes.
    let v = p_derived.two_column_isostasy.as_ref().unwrap().ocean_water_volume_m3;
    assert!(s2.ocean_volume_residual_m3.abs() <= 1e-9 * v);
}

// ---------------------------------------------------------------------------
// Measurement (ignored; lane K2 runs it explicitly per seed and mode).
// ---------------------------------------------------------------------------

fn write_f32(path: &Path, v: &[f64]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for x in v {
        f.write_all(&(*x as f32).to_le_bytes()).unwrap();
    }
}

#[test]
#[ignore]
fn measure_crust() {
    let seed: u64 = std::env::var("CRUST_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("CRUST_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("CRUST_MODE").unwrap();
    let opts = match mode.as_str() {
        "legacy" => PreviewOptions::legacy(),
        "twocol" => PreviewOptions::two_column_v0(),
        "bodies" => PreviewOptions::crust_bodies_v0(),
        // DIAGNOSTIC counterfactual (never a proposed setting): bodies with
        // stream-power K = 0 and hillslope D = 0, to separate isostasy and
        // seeding from the erosion budget (D4).
        "bodies_noerosion" => {
            let mut o = PreviewOptions::crust_bodies_v0();
            let mut s = pf_contracts::params::SurfaceBlock::default();
            s.erodibility = 0.0;
            s.hillslope_diffusivity_m2_per_yr = 0.0;
            o.surface_override = Some(s);
            o
        }
        // DIAGNOSTIC: bodies with the legacy single-column isostasy (which
        // has no root-depth check), to inspect columns where the
        // two-column mode refuses.
        "bodies_airy" => PreviewOptions {
            continental_bodies: Some(pf_contracts::params::ContinentalBodiesBlock::v0()),
            ..PreviewOptions::default()
        },
        other => panic!("unknown mode {other}"),
    };
    let out_dir = PathBuf::from(std::env::var("CRUST_OUT").unwrap());
    std::fs::create_dir_all(&out_dir).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: None,
        step_my: None,
    };
    // CRUST_T0_ONLY=1: the t = 0 state only (seeding + isostasy, no
    // kinematics), for seeds whose full run aborts (D5 / root refusals).
    if std::env::var("CRUST_T0_ONLY").map(|v| v == "1").unwrap_or(false) {
        dump_t0_only(&cfg, &opts, &out_dir, seed, &mode);
        return;
    }
    let t = Instant::now();
    let out = match simulate_with(&cfg, &opts) {
        Ok(o) => o,
        Err(e) => {
            std::fs::write(out_dir.join("ERROR.txt"), format!("{e}")).unwrap();
            return;
        }
    };
    let secs = t.elapsed().as_secs_f64();
    let grid = &out.grid;
    let sea = out.params.planet.sea_level_m;
    // Initial crust kind (t = 0), rebuilt from the same public seeding.
    let seeding = seed_plates(grid, &out.params, &cfg.seed_path()).unwrap();
    let prior = UpliftPrior::preview_v0();
    let state0 = match &out.params.continental_bodies {
        None => CrustState::seed(grid, &seeding, &prior).unwrap(),
        Some(b) => {
            CrustState::seed_bodies(
                grid,
                &seeding,
                &prior,
                b.resolved_area_fraction(&out.params.tectonics),
                b.rifted_margin_width_m,
            )
            .unwrap()
            .0
        }
    };
    write_f32(&out_dir.join("elevation_t0_m.f32"), &out.elevation_t0_m);
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_t0_m.f32"), &state0.thickness_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    write_f32(&out_dir.join("eroded_total_m.f32"), &out.state.eroded_total_m);
    let code = |k: &[CrustKindCodeV0]| k.iter().map(|x| x.code()).collect::<Vec<u8>>();
    std::fs::write(out_dir.join("kind.u8"), code(&out.state.kind)).unwrap();
    std::fs::write(out_dir.join("kind_t0.u8"), code(&state0.kind)).unwrap();
    std::fs::write(out_dir.join("boundary_class.u8"), &out.state.boundary_class).unwrap();
    let plate: Vec<u8> = out
        .state
        .plate
        .iter()
        .flat_map(|p| p.0.to_le_bytes())
        .collect();
    std::fs::write(out_dir.join("plate_id.u32"), plate).unwrap();
    let ocean = largest_ocean(grid, &out.elevation_m, sea).unwrap();
    std::fs::write(
        out_dir.join("ocean_connected.u8"),
        ocean.iter().map(|&b| b as u8).collect::<Vec<_>>(),
    )
    .unwrap();
    let ocean0 = largest_ocean(grid, &out.elevation_t0_m, sea).unwrap();
    std::fs::write(
        out_dir.join("ocean_connected_t0.u8"),
        ocean0.iter().map(|&b| b as u8).collect::<Vec<_>>(),
    )
    .unwrap();
    let row_areas: Vec<f64> = (0..grid.nlat()).map(|i| grid.cell_area_m2(i).unwrap()).collect();
    let bodies = out.continental_bodies.as_ref().map(|r| {
        serde_json::json!({
            "declared_area_fraction": r.declared_area_fraction,
            "achieved_area_fraction": r.achieved_area_fraction,
            "rifted_margin_width_m": r.rifted_margin_width_m,
            "interior_thickness_m": r.interior_thickness_m,
            "edge_thickness_m": r.edge_thickness_m,
            "continental_crust_volume_m3": r.continental_crust_volume_m3,
            "plateau_equivalent_volume_m3": r.plateau_equivalent_volume_m3,
            "coast_midpoints": r.coast_midpoints,
            "bodies": r.bodies.iter().map(|b| serde_json::json!({
                "plate": b.plate.0, "target_area_m2": b.target_area_m2,
                "achieved_area_m2": b.achieved_area_m2, "cells": b.cells,
                "crust_volume_m3": b.crust_volume_m3, "tapered_cells": b.tapered_cells,
            })).collect::<Vec<_>>(),
        })
    });
    let meta = serde_json::json!({
        "seed": seed,
        "mode": mode,
        "nlat": grid.nlat(),
        "nlon": grid.nlon(),
        "row_area_m2": row_areas,
        "total_area_m2": row_areas.iter().sum::<f64>() * grid.nlon() as f64,
        "sea_level_datum_m": sea,
        "ocean_fraction_param": out.params.tectonics.ocean_fraction,
        "simulate_seconds": secs,
        "sea_level_t0": out.sea_level_t0,
        "steps": out.steps,
        "continental_bodies": bodies,
        "summary": out.terrain.summary,
        "crust_volume": out.terrain.crust_volume,
        "approximations": out.terrain.approximations,
    });
    std::fs::write(out_dir.join("meta.json"), serde_json::to_vec_pretty(&meta).unwrap()).unwrap();
}

/// t = 0 only: the same public seeding and the same `elevation_on_grid`
/// that `simulate_with` uses for its initial elevation.
fn dump_t0_only(cfg: &PreviewConfig, opts: &PreviewOptions, out_dir: &Path, seed: u64, mode: &str) {
    let params = cfg.profile_with(opts).unwrap();
    let grid = cfg.grid(&params).unwrap();
    let seeding = seed_plates(&grid, &params, &cfg.seed_path()).unwrap();
    let prior = UpliftPrior::preview_v0();
    let (state0, bodies) = match &params.continental_bodies {
        None => (CrustState::seed(&grid, &seeding, &prior).unwrap(), None),
        Some(b) => {
            let (s, r) = CrustState::seed_bodies(
                &grid,
                &seeding,
                &prior,
                b.resolved_area_fraction(&params.tectonics),
                b.rifted_margin_width_m,
            )
            .unwrap();
            (s, Some(r))
        }
    };
    let mut areas = Vec::with_capacity(grid.cells());
    for i in 0..grid.nlat() {
        let a = grid.cell_area_m2(i).unwrap();
        for _ in 0..grid.nlon() {
            areas.push(a);
        }
    }
    let (e0, solve) = elevation_on_grid(&state0, &params, &areas).unwrap();
    let sea = params.planet.sea_level_m;
    write_f32(&out_dir.join("elevation_t0_m.f32"), &e0);
    write_f32(&out_dir.join("thickness_t0_m.f32"), &state0.thickness_m);
    std::fs::write(
        out_dir.join("kind_t0.u8"),
        state0.kind.iter().map(|k| k.code()).collect::<Vec<u8>>(),
    )
    .unwrap();
    let plate: Vec<u8> = state0.plate.iter().flat_map(|p| p.0.to_le_bytes()).collect();
    std::fs::write(out_dir.join("plate_id_t0.u32"), plate).unwrap();
    let ocean0 = largest_ocean(&grid, &e0, sea).unwrap();
    std::fs::write(
        out_dir.join("ocean_connected_t0.u8"),
        ocean0.iter().map(|&b| b as u8).collect::<Vec<_>>(),
    )
    .unwrap();
    let row_areas: Vec<f64> = (0..grid.nlat()).map(|i| grid.cell_area_m2(i).unwrap()).collect();
    let meta = serde_json::json!({
        "seed": seed, "mode": mode, "t0_only": true,
        "nlat": grid.nlat(), "nlon": grid.nlon(), "row_area_m2": row_areas,
        "sea_level_datum_m": sea,
        "ocean_fraction_param": params.tectonics.ocean_fraction,
        "sea_level_t0": solve,
        "continental_bodies": bodies.as_ref().map(|r| serde_json::json!({
            "declared_area_fraction": r.declared_area_fraction,
            "achieved_area_fraction": r.achieved_area_fraction,
            "continental_crust_volume_m3": r.continental_crust_volume_m3,
            "plateau_equivalent_volume_m3": r.plateau_equivalent_volume_m3,
            "tapered_cells": r.bodies.iter().map(|b| b.tapered_cells).sum::<u64>(),
            "cells": r.bodies.iter().map(|b| b.cells).sum::<u64>(),
        })),
    });
    std::fs::write(out_dir.join("meta_t0.json"), serde_json::to_vec_pretty(&meta).unwrap()).unwrap();
}
