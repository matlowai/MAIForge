//! Lane SD (2026-10-12): `erosion-deposition-v1` and `shelf-wave-base-v1`
//! through `erode_deposit_step`, on small synthetic worlds with known
//! structure. What must hold:
//! - the legacy law set is `erode_step_bounded`, bit for bit;
//! - G = 0 with legacy placement reproduces the legacy step to rounding
//!   (the new fluvial solver reduces to the old one);
//! - mass closes and the named partition sums to what eroded, for every law
//!   combination, and per cell dz = added - removed;
//! - the shelf law never builds a cell's compensated surface above the shelf
//!   top, builds a shelf body connected to the river mouth, and leaves no
//!   global remainder when the search bound is not reached;
//! - erosion-deposition stores sediment on the low-slope plain, not on
//!   the steep range;
//! - repeated runs are identical.

use pf_contracts::spatial::ElementId;
use pf_grid::GridSpec;
use pf_surface::erosion::{
    erode_deposit_step, erode_step_bounded, DiffusionStencil, ErosionDepositionParams,
    ErosionParams, ErosionStep, SurfaceLaws,
};
use pf_surface::routing::route;
use pf_surface::sediment::{CompensatedColumn, ShelfParams, SHELF_WAVE_BASE_DEPTH_M};

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

/// Ocean in the southern rows (depth `ocean_depth`), a nearly flat lowland
/// (1 m per row) north of it, then a steep range (400 m per row).
fn world(nlat: usize, nlon: usize, ocean_rows: usize, lowland_rows: usize, ocean_depth: f64) -> (GridSpec, Vec<f64>, Vec<bool>) {
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
            } else if i < ocean_rows + lowland_rows {
                elev[c] = 5.0 + (i - ocean_rows) as f64 + (j % 3) as f64 * 0.1;
            } else {
                let r = (i - ocean_rows - lowland_rows + 1) as f64;
                elev[c] = 5.0 + lowland_rows as f64 + 400.0 * r + (j % 3) as f64 * 3.0;
            }
        }
    }
    (grid, elev, ocean)
}

fn run(
    grid: &GridSpec,
    elev: &[f64],
    ocean: &[bool],
    sed: &[f64],
    p: &ErosionParams,
    laws: &SurfaceLaws,
) -> ErosionStep {
    let r = route(grid, elev, &vec![1.0; grid.cells()], ocean, 0.0).unwrap();
    let stencil = DiffusionStencil::build(grid).unwrap();
    erode_deposit_step(grid, &stencil, &r, elev, sed, ocean, 0.0, 2.0, p, laws, None).unwrap()
}

fn all_laws() -> Vec<(&'static str, SurfaceLaws)> {
    let shelf = Some(ShelfParams::v1(COLUMN));
    vec![
        ("shelf only", SurfaceLaws { erosion_deposition: None, shelf }),
        ("erosion-deposition only", SurfaceLaws { erosion_deposition: Some(ErosionDepositionParams::v1()), shelf: None }),
        ("both", SurfaceLaws { erosion_deposition: Some(ErosionDepositionParams::v1()), shelf }),
    ]
}

fn assert_closed(name: &str, s: &ErosionStep, areas: &[f64]) {
    let rep = s.deposition.as_ref().expect("non-legacy laws report a partition");
    let tol = 1.0e-9 * s.eroded_m3.max(1.0);
    let stored = rep.total_stored_m3();
    assert!((s.eroded_m3 - stored).abs() <= tol, "{name}: eroded {} vs stored {stored}", s.eroded_m3);
    assert!((s.deposited_m3 - rep.land_and_lake_m3()).abs() <= tol, "{name}: deposited");
    assert!((s.exported_m3 - rep.marine_m3()).abs() <= tol, "{name}: exported");
    assert!(
        (s.eroded_m3 - s.bedrock_eroded_m3 - s.sediment_eroded_m3).abs() <= tol,
        "{name}: bedrock + sediment split"
    );
    let added: f64 = s.sediment_added_m.iter().zip(areas).map(|(h, a)| h * a).sum();
    assert!((added - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: added volume");
    for c in 0..s.dz_m.len() {
        let expect = s.sediment_added_m[c] - s.sediment_removed_m[c] - s.bedrock_removed_m[c];
        assert!((s.dz_m[c] - expect).abs() < 1.0e-6, "{name}: cell {c}");
    }
}

#[test]
fn legacy_law_set_is_erode_step_bounded_bit_for_bit() {
    let (grid, elev, ocean) = world(12, 24, 3, 2, 3000.0);
    let n = grid.cells();
    let r = route(&grid, &elev, &vec![1.0; n], &ocean, 0.0).unwrap();
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let sed = vec![5.0; n];
    let rock = vec![30_000.0; n];
    let p = params(2.0, 1.0e6);
    for col in [None, Some(&rock[..])] {
        let a = erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, col).unwrap();
        let b = erode_deposit_step(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, &SurfaceLaws::legacy(), col).unwrap();
        assert_eq!(a, b);
        assert!(b.deposition.is_none());
    }
}

#[test]
fn g_zero_with_legacy_placement_reproduces_the_legacy_step_to_rounding() {
    let (grid, elev, ocean) = world(12, 24, 3, 2, 3000.0);
    let n = grid.cells();
    let r = route(&grid, &elev, &vec![1.0; n], &ocean, 0.0).unwrap();
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let sed = vec![5.0; n];
    let p = params(2.0, 1.0e6);
    let legacy = erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, None).unwrap();
    let laws = SurfaceLaws {
        erosion_deposition: Some(ErosionDepositionParams { deposition_coefficient_g: 0.0, ..ErosionDepositionParams::v1() }),
        shelf: None,
    };
    let g0 = erode_deposit_step(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, &laws, None).unwrap();
    let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(1.0);
    assert!(rel(legacy.eroded_m3, g0.eroded_m3) < 1e-12, "{} vs {}", legacy.eroded_m3, g0.eroded_m3);
    assert!(rel(legacy.exported_m3, g0.exported_m3) < 1e-12);
    assert!(rel(legacy.deposit_aggraded_m3, g0.deposit_aggraded_m3) < 1e-9);
    let worst = legacy.dz_m.iter().zip(&g0.dz_m).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    assert!(worst < 1.0e-6, "max |dz difference| {worst} m");
    assert!(legacy.exported_m3 > 0.0);
    let rep = g0.deposition.unwrap();
    assert_eq!(rep.land_fluvial_m3, 0.0, "G = 0 deposits nothing on land");
}

#[test]
fn every_law_combination_closes_mass_with_named_terms() {
    for (ocean_depth, k) in [(3000.0, 2.0), (300.0, 50.0)] {
        let (grid, elev, ocean) = world(12, 24, 3, 2, ocean_depth);
        let n = grid.cells();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let sed = vec![5.0; n];
        for (name, laws) in all_laws() {
            let s = run(&grid, &elev, &ocean, &sed, &params(k, 1.0e9), &laws);
            assert!(s.eroded_m3 > 0.0, "{name}: nothing eroded");
            assert_closed(name, &s, &stencil.areas_m2);
        }
    }
}

#[test]
fn shelf_law_builds_to_wave_base_from_the_mouth_and_never_above_it() {
    // A 300 m sea of five rows and a big range: supply fills some cells
    // but not the whole sea (at 3000 m and three rows each mouth only
    // part-fills its 15-degree entry cell; at 600 m it fills all of it).
    let (grid, elev, ocean) = world(12, 24, 5, 2, 300.0);
    let n = grid.cells();
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let sed = vec![0.0; n];
    let laws = SurfaceLaws { erosion_deposition: None, shelf: Some(ShelfParams::v1(COLUMN)) };
    let s = run(&grid, &elev, &ocean, &sed, &params(5.0, 0.0), &laws);
    let rep = s.deposition.clone().unwrap();
    assert_closed("shelf", &s, &stencil.areas_m2);
    assert!(rep.marine_search_bound_hits == 0 && rep.global_remainder_m3 == 0.0, "{rep:?}");
    let top = -SHELF_WAVE_BASE_DEPTH_M;
    let mut filled = vec![false; n];
    for c in 0..n {
        if !ocean[c] {
            continue;
        }
        let surf = COLUMN.surface_after(elev[c], s.sediment_added_m[c], 0.0);
        assert!(surf <= top + 1.0e-6, "ocean cell {c} built to {surf} m, above the shelf top");
        if s.sediment_added_m[c] > 0.0 && (surf - top).abs() < 1.0e-6 {
            filled[c] = true;
        }
    }
    let nfilled = filled.iter().filter(|&&f| f).count();
    assert!(nfilled > 0, "no cell reached the shelf top; {rep:?}");
    assert_eq!(nfilled as u64, rep.shelf_cells_filled);
    // Every filled cell touches land or another filled cell (a shelf body
    // grown from the coast, not scattered fills).
    for c in 0..n {
        if !filled[c] {
            continue;
        }
        let (i, j) = grid.ij(c as u64).unwrap();
        let ok = grid.neighbors(i, j).unwrap().into_iter().any(|(ti, tj, _)| {
            let nb = ti * grid.nlon() + tj;
            !ocean[nb] || filled[nb]
        });
        assert!(ok, "filled cell {c} is detached from the coast and the shelf");
    }
    // Compensated accommodation: the thickness laid in a filled cell is the
    // depth times (rho_m - rho_w) / (rho_m - rho_s), not the depth.
    let c = (0..n).find(|&c| filled[c]).unwrap();
    let expect = (300.0 - SHELF_WAVE_BASE_DEPTH_M) * (3300.0 - 1030.0) / (3300.0 - 2800.0);
    assert!((s.sediment_added_m[c] - expect).abs() < 1.0e-6, "{} vs {expect}", s.sediment_added_m[c]);
}

#[test]
fn shallow_sea_is_crossed_not_filled_above_wave_base() {
    // Ocean only 50 m deep: already above the shelf top, so nothing may be
    // laid until the search reaches deeper water; with none, it is the
    // declared remainder.
    let (grid, elev, ocean) = world(12, 24, 3, 2, 50.0);
    let n = grid.cells();
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let laws = SurfaceLaws { erosion_deposition: None, shelf: Some(ShelfParams::v1(COLUMN)) };
    let s = run(&grid, &elev, &ocean, &vec![0.0; n], &params(5.0, 0.0), &laws);
    let rep = s.deposition.clone().unwrap();
    assert_closed("shallow", &s, &stencil.areas_m2);
    assert_eq!(rep.marine_shelf_m3 + rep.marine_slope_m3 + rep.marine_deep_m3, 0.0);
    assert!(rep.global_remainder_m3 > 0.0);
}

#[test]
fn erosion_deposition_stores_sediment_on_the_plain_not_the_range() {
    let (grid, elev, ocean) = world(24, 48, 4, 6, 3000.0);
    let n = grid.cells();
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let sed = vec![0.0; n];
    let p = params(2.0, 0.0);
    let legacy = run(&grid, &elev, &ocean, &sed, &p, &SurfaceLaws { erosion_deposition: None, shelf: Some(ShelfParams::v1(COLUMN)) });
    let ed = run(&grid, &elev, &ocean, &sed, &p, &SurfaceLaws { erosion_deposition: Some(ErosionDepositionParams::v1()), shelf: Some(ShelfParams::v1(COLUMN)) });
    let lr = legacy.deposition.unwrap();
    let er = ed.deposition.clone().unwrap();
    assert_closed("ed", &ed, &stencil.areas_m2);
    assert_eq!(lr.land_fluvial_m3, 0.0);
    assert!(er.land_fluvial_m3 > 0.0, "{er:?}");
    assert!(er.converged, "sweeps {} change {}", er.sweeps, er.last_sweep_max_change_m);
    // Where did land deposition land: area-weighted by row class.
    let (mut plain, mut range) = (0.0, 0.0);
    for c in 0..n {
        if ocean[c] || ed.sediment_added_m[c] <= 0.0 {
            continue;
        }
        let row = c / grid.nlon();
        let v = ed.sediment_added_m[c] * stencil.areas_m2[c];
        if row < 4 + 6 {
            plain += v;
        } else {
            range += v;
        }
    }
    assert!(plain > range, "plain {plain} vs range {range}");
    // Less reaches the sea than without deposition (storage on land).
    assert!(er.delivered_to_ocean_m3 < lr.delivered_to_ocean_m3);
}

#[test]
fn laws_are_deterministic() {
    let (grid, elev, ocean) = world(12, 24, 3, 2, 3000.0);
    let sed = vec![3.0; grid.cells()];
    for (_, laws) in all_laws() {
        let a = run(&grid, &elev, &ocean, &sed, &params(2.0, 1.0e8), &laws);
        let b = run(&grid, &elev, &ocean, &sed, &params(2.0, 1.0e8), &laws);
        let c = run(&grid, &elev, &ocean, &sed, &params(2.0, 1.0e8), &laws);
        assert_eq!(a, b);
        assert_eq!(b, c);
    }
}

#[test]
fn routing_shape_sanity_for_the_fixture() {
    // Guard on the fixture itself: every land cell reaches the ocean.
    let (grid, elev, ocean) = world(12, 24, 3, 2, 3000.0);
    let r = route(&grid, &elev, &vec![1.0; grid.cells()], &ocean, 0.0).unwrap();
    for c in 0..grid.cells() {
        if ocean[c] {
            continue;
        }
        let mut cur = c;
        let mut hops = 0;
        while let Some(ElementId(j)) = r.flow_to[cur] {
            cur = j as usize;
            hops += 1;
            assert!(hops <= grid.cells());
        }
        assert!(ocean[cur], "cell {c} ends at {cur}");
    }
}

#[test]
fn headwater_cells_satisfy_the_erosion_deposition_equation() {
    // Known answer on the real solver: a cell with no donors has A = a, so
    // G a / (p A) = G and no inflow, and its solved level must satisfy
    //     h = (h0 (1 + G) + f h_r) / (1 + f + G),  f = K A^m dt / L,
    // against its receiver's SOLVED level h_r (or fill level if wet).
    let (grid, elev, ocean) = world(24, 48, 4, 6, 3000.0);
    let n = grid.cells();
    let r = route(&grid, &elev, &vec![1.0; n], &ocean, 0.0).unwrap();
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let (k, g, dt) = (2.0, 1.0, 2.0);
    let laws = SurfaceLaws { erosion_deposition: Some(ErosionDepositionParams::v1()), shelf: None };
    let s = erode_deposit_step(&grid, &stencil, &r, &elev, &vec![0.0; n], &ocean, 0.0, dt, &params(k, 0.0), &laws, None).unwrap();
    let h: Vec<f64> = elev.iter().zip(&s.dz_m).map(|(a, b)| a + b).collect();
    let mut has_donor = vec![false; n];
    for c in 0..n {
        if let Some(ElementId(j)) = r.flow_to[c] {
            has_donor[j as usize] = true;
        }
    }
    let mut checked = 0;
    for c in 0..n {
        if ocean[c] || has_donor[c] || r.fill_level_m[c] > elev[c] + 1e-9 {
            continue;
        }
        let Some(ElementId(j)) = r.flow_to[c] else { continue };
        let j = j as usize;
        let hr = if ocean[j] || r.fill_level_m[j] > elev[j] + 1e-9 { r.fill_level_m[j] } else { h[j] };
        let (ci, cj) = grid.ij(c as u64).unwrap();
        let (ji, jj) = grid.ij(j as u64).unwrap();
        let len = grid.center_distance_m(ci, cj, ji, jj).unwrap().max(1.0);
        let f = k * r.contributing_area_m2[c].powf(0.5) * dt / len;
        let expect = if elev[c] > hr {
            (elev[c] * (1.0 + g) + f * hr) / (1.0 + f + g)
        } else {
            elev[c]
        };
        assert!((h[c] - expect).abs() <= 1e-9 * expect.abs().max(1.0), "cell {c}: {} vs {expect}", h[c]);
        checked += 1;
    }
    assert!(checked > 10, "only {checked} headwater cells");
}
