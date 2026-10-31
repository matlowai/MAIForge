//! SEED-EROSION-BUG (lane E, 2026-10-11): the per-step available-material
//! bound and the bedrock/sediment split of `erode_step_bounded`.
//!
//! What must hold:
//! - with no column supplied (and with a column too thick to bind) the
//!   step equals the legacy `erode_step` exactly (default path unchanged);
//! - with a thin column the bound binds, no cell loses more bedrock than
//!   it holds, and mass still closes (eroded = deposited + exported);
//! - the bedrock and sediment shares sum to `eroded_m3`;
//! - per cell, dz = added - removed (no material from nowhere).

use pf_grid::GridSpec;
use pf_surface::erosion::{erode_step, erode_step_bounded, DiffusionStencil, ErosionParams, ErosionStep};
use pf_surface::routing::{route, Routing};

fn params(k: f64, d: f64) -> ErosionParams {
    ErosionParams {
        k_per_my: k,
        m: 0.5,
        n: 1.0,
        diffusivity_m2_per_my: d,
        stability_limit: 0.5,
    }
}

/// Southern rows ocean at -3000 m; land rising north by 200 m per row.
fn ramp_world() -> (GridSpec, Vec<f64>, Vec<bool>, Routing, DiffusionStencil) {
    let grid = GridSpec::new(12, 24, 6.371e6).unwrap();
    let n = grid.cells();
    let mut ocean = vec![false; n];
    let mut elev = vec![0.0; n];
    for i in 0..grid.nlat() {
        for j in 0..grid.nlon() {
            let c = i * grid.nlon() + j;
            if i < 3 {
                ocean[c] = true;
                elev[c] = -3000.0;
            } else {
                elev[c] = 200.0 * (i as f64 - 2.0) + (j % 3) as f64 * 5.0;
            }
        }
    }
    let r = route(&grid, &elev, &vec![1.0e6; n], &ocean, 0.0).unwrap();
    let stencil = DiffusionStencil::build(&grid).unwrap();
    (grid, elev, ocean, r, stencil)
}

fn closes(step: &ErosionStep) {
    let residual = step.eroded_m3 - step.deposited_m3 - step.exported_m3;
    assert!(
        residual.abs() <= 1.0e-9 * step.eroded_m3.max(1.0),
        "eroded {} != deposited {} + exported {} (residual {residual})",
        step.eroded_m3,
        step.deposited_m3,
        step.exported_m3
    );
    let split = step.bedrock_eroded_m3 + step.sediment_eroded_m3 - step.eroded_m3;
    assert!(
        split.abs() <= 1.0e-12 * step.eroded_m3.max(1.0),
        "bedrock {} + sediment {} != eroded {}",
        step.bedrock_eroded_m3,
        step.sediment_eroded_m3,
        step.eroded_m3
    );
    for c in 0..step.dz_m.len() {
        let expect = step.sediment_added_m[c] - step.sediment_removed_m[c] - step.bedrock_removed_m[c];
        assert!(
            (step.dz_m[c] - expect).abs() < 1.0e-9,
            "cell {c}: dz {} vs added - removed {expect}",
            step.dz_m[c]
        );
    }
}

#[test]
fn no_column_and_a_non_binding_column_equal_the_legacy_step() {
    let (grid, elev, ocean, r, stencil) = ramp_world();
    let n = grid.cells();
    let sed = vec![10.0; n];
    for p in [params(2.0, 1.0e6), params(1.0e3, 0.0), params(0.0, 1.0e12)] {
        let legacy = erode_step(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p).unwrap();
        let none =
            erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, None).unwrap();
        assert_eq!(legacy, none, "None must be the legacy step");
        let thick = vec![1.0e9; n];
        let wide =
            erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, Some(&thick))
                .unwrap();
        assert_eq!(wide.bedrock_bound_hits, 0);
        assert_eq!(legacy, wide, "a column that never binds must not change the step");
        closes(&legacy);
    }
}

#[test]
fn thin_column_binds_the_stream_power_lowering_and_mass_closes() {
    let (grid, elev, ocean, r, stencil) = ramp_world();
    let n = grid.cells();
    let sed = vec![0.0; n];
    let rock = vec![50.0; n];
    let p = params(1.0e3, 0.0);
    let free = erode_step(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 50.0, &p).unwrap();
    assert!(
        free.bedrock_removed_m.iter().any(|&v| v > 50.0),
        "fixture must overdraw a 50 m column without the bound"
    );
    let b = erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 50.0, &p, Some(&rock))
        .unwrap();
    assert!(b.bedrock_bound_hits > 0, "bound never engaged");
    for c in 0..n {
        assert!(
            b.bedrock_removed_m[c] <= rock[c] + 1.0e-9,
            "cell {c} lost {} m of a {} m column",
            b.bedrock_removed_m[c],
            rock[c]
        );
    }
    closes(&b);
    assert!(b.eroded_m3 < free.eroded_m3);
}

#[test]
fn sediment_cover_counts_toward_what_may_be_removed() {
    let (grid, elev, ocean, r, stencil) = ramp_world();
    let n = grid.cells();
    let sed = vec![30.0; n];
    let rock = vec![20.0; n];
    let p = params(1.0e3, 0.0);
    let b = erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 50.0, &p, Some(&rock))
        .unwrap();
    for c in 0..n {
        if ocean[c] {
            continue;
        }
        let removed = b.bedrock_removed_m[c] + b.sediment_removed_m[c];
        assert!(removed <= 50.0 + 1.0e-9, "cell {c} removed {removed} of 50 m");
        assert!(b.bedrock_removed_m[c] <= 20.0 + 1.0e-9);
    }
    assert!(b.sediment_eroded_m3 > 0.0 && b.bedrock_eroded_m3 > 0.0);
    closes(&b);
}

#[test]
fn diffusion_donors_are_scaled_to_what_they_hold() {
    let (grid, elev, ocean, r, stencil) = ramp_world();
    let n = grid.cells();
    let sed = vec![0.0; n];
    let rock = vec![1.0; n];
    // K = 0 isolates diffusion; D large enough to move far more than 1 m.
    let p = params(0.0, 1.0e12);
    let free = erode_step(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p).unwrap();
    assert!(free.bedrock_removed_m.iter().any(|&v| v > 1.0));
    let b = erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, Some(&rock))
        .unwrap();
    assert!(b.bedrock_bound_hits > 0);
    for c in 0..n {
        assert!(
            b.bedrock_removed_m[c] <= 1.0 + 1.0e-9,
            "cell {c} lost {} m of a 1 m column",
            b.bedrock_removed_m[c]
        );
    }
    closes(&b);
}

#[test]
fn bounded_step_is_deterministic_and_refuses_bad_columns() {
    let (grid, elev, ocean, r, stencil) = ramp_world();
    let n = grid.cells();
    let sed = vec![5.0; n];
    let rock = vec![80.0; n];
    let p = params(5.0, 1.0e9);
    let a = erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, Some(&rock))
        .unwrap();
    let b = erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, Some(&rock))
        .unwrap();
    assert_eq!(a, b);
    let mut bad = rock.clone();
    bad[7] = -1.0;
    assert!(
        erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, Some(&bad)).is_err()
    );
    bad[7] = f64::NAN;
    assert!(
        erode_step_bounded(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, Some(&bad)).is_err()
    );
    assert!(erode_step_bounded(
        &grid,
        &stencil,
        &r,
        &elev,
        &sed,
        &ocean,
        0.0,
        2.0,
        &p,
        Some(&rock[..3])
    )
    .is_err());
}
