//! Lane SR (SEED-SEDIMENT-RETAIN, 2026-10-13): `continental-retention-base-level-v1`
//! (`pf_surface::erosion::retain`) on small synthetic worlds. Pre-registered
//! falsifiers (`generated-test-data/PREREG.txt`):
//! - F1 mass closes: eroded = deposited + exported, per cell dz = added - removed;
//! - F3 the taper (rule 3, h = max(local relief, height above base level)):
//!   steep ground and a flat plateau >= H above the sea lower exactly as in
//!   the legacy step; a low coastal cell lowers by phi(h) times the legacy
//!   lowering;
//! - F4 nothing reaches oceanic cells while reachable continental
//!   accommodation (graded plain, flooded continental shelf) is unfilled;
//! - F5 determinism;
//! - with the taper off (tiny H), no continent and no lakes, the step equals
//!   the legacy step to rounding (only the order of volume sums differs).

use pf_grid::GridSpec;
use pf_surface::erosion::retain::{erode_retain_step, local_relief, RetentionParams, RetentionReport};
use pf_surface::erosion::{erode_step_bounded, DiffusionStencil, ErosionParams, ErosionStep};
use pf_surface::routing::route;
use pf_surface::sediment::CompensatedColumn;

const COLUMN: CompensatedColumn = CompensatedColumn {
    mantle_density_kg_m3: 3300.0,
    sediment_density_kg_m3: 2800.0,
    water_density_kg_m3: 1030.0,
};

fn params(k: f64, d: f64) -> ErosionParams {
    ErosionParams {
        k_per_my: k,
        m: 0.5,
        n: 1.0,
        diffusivity_m2_per_my: d,
        stability_limit: 0.5,
    }
}

/// Rows from the south: `ocean_rows` of ocean at `ocean_depth`, then land
/// rising by `row_rise` m per row from `coast` m (with a small zonal
/// texture so flow directions are well defined).
fn ramp(nlat: usize, nlon: usize, ocean_rows: usize, ocean_depth: f64, coast: f64, row_rise: f64) -> (GridSpec, Vec<f64>, Vec<bool>) {
    let grid = GridSpec::new(nlat, nlon, 6.371e6).unwrap();
    let n = grid.cells();
    let mut ocean = vec![false; n];
    let mut elev = vec![0.0; n];
    for i in 0..nlat {
        for j in 0..nlon {
            let c = i * nlon + j;
            if i < ocean_rows {
                ocean[c] = true;
                elev[c] = -ocean_depth;
            } else {
                elev[c] = coast + row_rise * (i - ocean_rows) as f64 + (j % 3) as f64 * 0.5;
            }
        }
    }
    (grid, elev, ocean)
}

fn run(
    grid: &GridSpec,
    elev: &[f64],
    ocean: &[bool],
    cont: &[bool],
    p: &ErosionParams,
    rp: &RetentionParams,
) -> (ErosionStep, RetentionReport, ErosionStep) {
    let n = grid.cells();
    let r = route(grid, elev, &vec![1.0; n], ocean, 0.0).unwrap();
    let stencil = DiffusionStencil::build(grid).unwrap();
    let sed = vec![0.0; n];
    let rock = vec![30_000.0; n];
    let (a, rep) = erode_retain_step(grid, &stencil, &r, elev, &sed, ocean, 0.0, 2.0, p, Some(&rock), cont, rp).unwrap();
    let legacy = erode_step_bounded(grid, &stencil, &r, elev, &sed, ocean, 0.0, 2.0, p, Some(&rock)).unwrap();
    (a, rep, legacy)
}

fn assert_closed(name: &str, s: &ErosionStep, rr: &RetentionReport, areas: &[f64]) {
    let tol = 1.0e-9 * s.eroded_m3.max(1.0);
    assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: eroded {} vs {} + {}", s.eroded_m3, s.deposited_m3, s.exported_m3);
    assert!((s.eroded_m3 - s.bedrock_eroded_m3 - s.sediment_eroded_m3).abs() <= tol, "{name}: split");
    let added: f64 = s.sediment_added_m.iter().zip(areas).map(|(h, a)| h * a).sum();
    assert!((added - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: added volume");
    assert!((added - rr.added_on_continental_m3 - rr.added_on_oceanic_m3).abs() <= tol, "{name}: kind split");
    let named = rr.continental_grade_m3 + rr.lake_m3 + rr.land_floor_m3 + rr.continental_shelf_m3 + rr.legacy_marine_m3;
    // Hillslope moves land-to-land are not in the named terms; with D = 0
    // the named terms are everything.
    assert!(named <= s.eroded_m3 + tol, "{name}: named {named} > eroded {}", s.eroded_m3);
    for c in 0..s.dz_m.len() {
        let expect = s.sediment_added_m[c] - s.sediment_removed_m[c] - s.bedrock_removed_m[c];
        assert!((s.dz_m[c] - expect).abs() < 1.0e-6, "{name}: cell {c}");
    }
    assert!(s.deposition.is_none());
}

#[test]
fn taper_leaves_high_ground_alone_and_scales_low_coastal_cells() {
    let p = params(2.0, 0.0);
    let rp = RetentionParams::v1(COLUMN);
    // Steep: 1100 m per row from 1500 m at the coast (local relief >= H everywhere).
    let (grid, elev, ocean) = ramp(12, 24, 3, 3000.0, 1500.0, 1100.0);
    let n = grid.cells();
    let cont = vec![false; n];
    let (a, rr, legacy) = run(&grid, &elev, &ocean, &cont, &p, &rp);
    assert_eq!(rr.tapered_cells, 0);
    assert_eq!(a.bedrock_removed_m, legacy.bedrock_removed_m, "high ground must erode exactly as legacy");
    assert!(a.eroded_m3 > 0.0);
    assert_closed("high", &a, &rr, &DiffusionStencil::build(&grid).unwrap().areas_m2);
    // Flat high plateau (1 cm per row, 1500 m above the sea): low local
    // relief inside, but >= H above base level, so unchanged (rule 3).
    let (grid, elev, ocean) = ramp(12, 24, 3, 3000.0, 1500.0, 0.01);
    let (a, rr, legacy) = run(&grid, &elev, &ocean, &cont, &p, &rp);
    assert_eq!(rr.tapered_cells, 0);
    assert_eq!(a.bedrock_removed_m, legacy.bedrock_removed_m, "a high flat plateau must erode exactly as legacy");
    // Low: a coastal lowland from 200 m rising 150 m per row; cells draining
    // straight to the sea (same receiver level) lower by exactly phi(r) x legacy.
    let (grid, elev, ocean) = ramp(12, 24, 3, 3000.0, 200.0, 150.0);
    let (a, rr, legacy) = run(&grid, &elev, &ocean, &cont, &p, &rp);
    assert!(rr.tapered_cells > 0);
    let r = route(&grid, &elev, &vec![1.0; n], &ocean, 0.0).unwrap();
    let relief = local_relief(&grid, &r, &elev, &ocean, 0.0).unwrap();
    let mut checked = 0;
    for c in 0..n {
        let Some(to) = r.flow_to[c] else { continue };
        if ocean[c] || !ocean[to.0 as usize] || legacy.bedrock_removed_m[c] <= 0.0 {
            continue;
        }
        let phi = rp.taper(relief[c]);
        assert!(relief[c] >= elev[c] + 149.0 && relief[c] < elev[c] + 152.0, "relief {} at {}", relief[c], elev[c]);
        assert!(phi < 1.0);
        let want = phi * legacy.bedrock_removed_m[c];
        assert!((a.bedrock_removed_m[c] - want).abs() <= 1e-9 * want.max(1.0), "cell {c}: {} vs {want}", a.bedrock_removed_m[c]);
        checked += 1;
    }
    assert!(checked >= 10, "checked {checked}");
    assert!(a.eroded_m3 < legacy.eroded_m3);
    assert!(rr.taper_withheld_m3 > 0.0);
    // No cell is lowered more than the legacy step lowers it.
    for c in 0..n {
        assert!(a.bedrock_removed_m[c] <= legacy.bedrock_removed_m[c] + 1e-9);
    }
}

#[test]
fn without_taper_continent_or_lakes_it_is_the_legacy_step_to_rounding() {
    let p = params(2.0, 1.0e6);
    let mut rp = RetentionParams::v1(COLUMN);
    rp.base_level_relief_scale_m = 1.0e-9;
    let (grid, elev, ocean) = ramp(12, 24, 3, 3000.0, 20.0, 250.0);
    let cont = vec![false; grid.cells()];
    let (a, rr, legacy) = run(&grid, &elev, &ocean, &cont, &p, &rp);
    let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(1.0);
    assert_eq!(rr.lake_m3, 0.0, "fixture must have no lakes");
    assert!(rel(a.eroded_m3, legacy.eroded_m3) < 1e-12);
    assert!(rel(a.exported_m3, legacy.exported_m3) < 1e-12);
    assert!(rel(a.deposit_aggraded_m3, legacy.deposit_aggraded_m3) < 1e-9);
    for c in 0..grid.cells() {
        assert!((a.dz_m[c] - legacy.dz_m[c]).abs() < 1e-6, "cell {c}: {} vs {}", a.dz_m[c], legacy.dz_m[c]);
    }
    assert_closed("no-taper", &a, &rr, &DiffusionStencil::build(&grid).unwrap().areas_m2);
}

/// A continental world: a range, a nearly flat continental lowland below
/// grade, a flooded continental shelf (two rows at 600 m), then oceanic
/// abyss. Sediment fills the lowland toward grade and the flooded shelf
/// toward wave base before any reaches oceanic cells.
#[test]
fn continental_accommodation_fills_before_oceanic_cells() {
    let nlat = 16;
    let nlon = 24;
    let grid = GridSpec::new(nlat, nlon, 6.371e6).unwrap();
    let n = grid.cells();
    let mut ocean = vec![false; n];
    let mut cont = vec![true; n];
    let mut elev = vec![0.0; n];
    for i in 0..nlat {
        for j in 0..nlon {
            let c = i * nlon + j;
            let t = (j % 3) as f64 * 0.01;
            if i < 3 {
                ocean[c] = true;
                cont[c] = false;
                elev[c] = -4000.0;
            } else if i < 5 {
                ocean[c] = true;
                elev[c] = -600.0;
            } else if i < 9 {
                // lowland, 1 cm per row: far below the S_a = 1e-4 grade
                elev[c] = 5.0 + 0.01 * (i - 5) as f64 + t;
            } else {
                elev[c] = 1500.0 + 800.0 * (i - 9) as f64 + t * 100.0;
            }
        }
    }
    // A small K keeps the delivered volume below the continental capacity,
    // so the "nothing on oceanic cells" clause below is exercised.
    let p = params(0.05, 0.0);
    let rp = RetentionParams::v1(COLUMN);
    let (a, rr, _legacy) = run(&grid, &elev, &ocean, &cont, &p, &rp);
    let areas = DiffusionStencil::build(&grid).unwrap().areas_m2;
    assert_closed("continental", &a, &rr, &areas);
    assert!(rr.continental_grade_m3 > 0.0, "{rr:?}");
    assert!(rr.mouths > 0 && rr.mouths_without_continental_shelf == 0, "{rr:?}");
    let shelf_cap: f64 = (0..n)
        .filter(|&c| ocean[c] && cont[c])
        .map(|c| COLUMN.thickness_to_raise(elev[c], -100.0, 0.0) * areas[c])
        .sum();
    let oceanic: f64 = (0..n).filter(|&c| !cont[c]).map(|c| a.sediment_added_m[c] * areas[c]).sum();
    assert!(rr.continental_shelf_m3 < shelf_cap * (1.0 - 1e-9), "fixture must leave shelf room: {rr:?}");
    {
        assert_eq!(oceanic, 0.0, "oceanic cells received sediment while the flooded shelf had room: {rr:?}");
        assert_eq!(rr.legacy_marine_m3, 0.0);
    }
    assert!(rr.continental_shelf_m3 > 0.0, "{rr:?}");
    // The flooded shelf never builds above wave base (compensated surface).
    for c in 0..n {
        if ocean[c] && cont[c] && a.sediment_added_m[c] > 0.0 {
            let s = COLUMN.surface_after(elev[c], a.sediment_added_m[c], 0.0);
            assert!(s <= -100.0 + 1e-6, "cell {c} shelf surface {s}");
        }
    }
    // Lowland deposits never raise a cell's compensated surface above grade:
    // receiver level before its own fill (<= 5.1 m) + S_a L (L < 2000 km here).
    for c in 0..n {
        if !ocean[c] && cont[c] && a.sediment_added_m[c] > 0.0 {
            let s = COLUMN.surface_after(elev[c] - a.bedrock_removed_m[c], a.sediment_added_m[c], 0.0);
            assert!(s <= 5.1 + rp.alluvial_grade_slope * 2.0e6, "cell {c} lowland surface {s}");
        }
    }
    // Same world, oceanic coast: everything passes to the legacy placement.
    let none = vec![false; n];
    let (b, rrb, legacy) = run(&grid, &elev, &ocean, &none, &p, &rp);
    assert_eq!(rrb.continental_grade_m3 + rrb.continental_shelf_m3, 0.0);
    assert!(rrb.legacy_marine_m3 > 0.0);
    assert!(b.eroded_m3 <= legacy.eroded_m3 + 1e-6);
    assert_closed("oceanic coast", &b, &rrb, &areas);
}

#[test]
fn determinism_and_validation() {
    let p = params(2.0, 1.0e6);
    let rp = RetentionParams::v1(COLUMN);
    let (grid, elev, ocean) = ramp(12, 24, 3, 3000.0, 150.0, 200.0);
    let cont: Vec<bool> = (0..grid.cells()).map(|c| c % 5 != 0).collect();
    let (a, ra, _) = run(&grid, &elev, &ocean, &cont, &p, &rp);
    let (b, rb, _) = run(&grid, &elev, &ocean, &cont, &p, &rp);
    assert_eq!(a, b);
    assert_eq!(ra, rb);
    let mut bad = rp;
    bad.taper_exponent = 0.0;
    assert!(bad.validate().is_err());
    let mut bad = rp;
    bad.alluvial_grade_slope = -1.0;
    assert!(bad.validate().is_err());
    assert_eq!(rp.taper(1000.0), 1.0);
    assert_eq!(rp.taper(5000.0), 1.0);
    assert_eq!(rp.taper(-5.0), 0.0);
    assert!((rp.taper(300.0) - 0.09).abs() < 1e-12);
}
