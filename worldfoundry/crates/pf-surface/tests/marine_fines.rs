//! Q25 integration-test draft. Apply as pf-surface/tests/marine_fines.rs.
//! Exercises the native retention step; no copied world or replacement solver.
use pf_grid::GridSpec;
use pf_surface::erosion::fluvial_trace::TraceClock;
use pf_surface::erosion::marine_fines::{MarineFinesParams, MarineFinesReport};
use pf_surface::erosion::retain::{
    erode_retain_step_with_lakes, erode_retain_step_with_marine_fines, RetentionParams,
};
use pf_surface::erosion::{DiffusionStencil, ErosionParams, ErosionStep};
use pf_surface::routing::route;
use pf_surface::sediment::CompensatedColumn;

const COLUMN: CompensatedColumn = CompensatedColumn {
    mantle_density_kg_m3: 3300.0,
    sediment_density_kg_m3: 2800.0,
    water_density_kg_m3: 1030.0,
};

struct Fixture {
    grid: GridSpec,
    elevation: Vec<f64>,
    sediment: Vec<f64>,
    ocean: Vec<bool>,
    continental: Vec<bool>,
    available: Vec<f64>,
}

impl Fixture {
    fn new(two_oceans: bool) -> Self {
        let grid = GridSpec::new(16, 24, 6_371_000.0).unwrap();
        let n = grid.cells();
        let mut f = Self {
            grid,
            elevation: vec![0.0; n],
            sediment: vec![0.0; n],
            ocean: vec![false; n],
            continental: vec![true; n],
            available: vec![0.0; n],
        };
        for row in 0..16 {
            for lon in 0..24 {
                let c = row * 24 + lon;
                let texture = (lon % 3) as f64 * 0.01;
                match row {
                    0..=2 => {
                        f.elevation[c] = -4000.0;
                        f.ocean[c] = true;
                        f.continental[c] = false;
                    }
                    3..=4 => {
                        f.elevation[c] = -600.0;
                        f.ocean[c] = true;
                    }
                    13 if two_oceans => {
                        f.elevation[c] = -600.0;
                        f.ocean[c] = true;
                    }
                    14..=15 if two_oceans => {
                        f.elevation[c] = -4000.0;
                        f.ocean[c] = true;
                        f.continental[c] = false;
                    }
                    _ => {
                        let rise = if two_oceans && row > 8 {
                            12 - row
                        } else {
                            row - 5
                        };
                        f.elevation[c] = 100.0 + 600.0 * rise as f64 + texture;
                        // Only the southern catchment supplies finite erodible rock.
                        if !two_oceans || row <= 8 {
                            f.available[c] = 10.0;
                        }
                    }
                }
            }
        }
        f
    }
}

fn params() -> ErosionParams {
    ErosionParams {
        k_per_my: 0.05,
        m: 0.5,
        n: 1.0,
        diffusivity_m2_per_my: 0.0,
        stability_limit: 0.5,
    }
}

fn retention() -> RetentionParams {
    let mut p = RetentionParams::v1(COLUMN);
    p.base_level_relief_scale_m = 1e-9;
    p.alluvial_grade_slope = 0.0;
    p
}

fn assert_closed(step: &ErosionStep, stencil: &DiffusionStencil, available: &[f64]) {
    let tol = 1.0 + 1e-10 * step.eroded_m3;
    assert!((step.eroded_m3 - step.deposited_m3 - step.exported_m3).abs() <= tol);
    let added: f64 = step
        .sediment_added_m
        .iter()
        .zip(&stencil.areas_m2)
        .map(|(h, a)| h * a)
        .sum();
    assert!((added - step.eroded_m3).abs() <= tol);
    for c in 0..available.len() {
        assert!(step.bedrock_removed_m[c] <= available[c] + 1e-9);
        assert!(
            (step.dz_m[c]
                - (step.sediment_added_m[c]
                    - step.bedrock_removed_m[c]
                    - step.sediment_removed_m[c]))
                .abs()
                < 1e-7
        );
    }
}

fn assert_receipt(report: &MarineFinesReport) {
    let tolerance = 1.0 + 1e-10 * report.requested_m3;
    assert!(
        (report.requested_m3 - report.deposited_m3 - report.returned_to_shelf_m3).abs()
            <= tolerance
    );
    let mut received = 0.0;
    for component in &report.components {
        let source: f64 = component.sources.iter().map(|s| s.volume_m3).sum();
        let receiver: f64 = component.receivers.iter().map(|s| s.volume_m3).sum();
        assert!((source - receiver).abs() <= tolerance);
        assert!((receiver - component.deposited_m3).abs() <= tolerance);
        assert!(
            (component.requested_m3 - source - component.returned_to_shelf_m3).abs() <= tolerance
        );
        assert!(component
            .sources
            .iter()
            .all(|s| s.volume_m3 >= 0.0 && s.volume_m3.is_finite()));
        assert!(component
            .receivers
            .iter()
            .all(|s| s.volume_m3 > 0.0 && s.volume_m3.is_finite()));
        received += receiver;
    }
    assert!((received - report.deposited_m3).abs() <= tolerance);
}

#[test]
fn disabled_and_zero_fraction_are_exact_existing_api() {
    let f = Fixture::new(false);
    let r = route(
        &f.grid,
        &f.elevation,
        &vec![1.0; f.grid.cells()],
        &f.ocean,
        0.0,
    )
    .unwrap();
    let s = DiffusionStencil::build(&f.grid).unwrap();
    let old = erode_retain_step_with_lakes(
        &f.grid,
        &s,
        &r,
        &f.elevation,
        &f.sediment,
        &f.ocean,
        0.0,
        1.0,
        &params(),
        Some(&f.available),
        &f.continental,
        &retention(),
        None,
        false,
        None,
    )
    .unwrap();
    let zero = MarineFinesParams {
        fraction: 0.0,
        min_depth_m: 2500.0,
    };
    for setting in [None, Some(&zero)] {
        let new = erode_retain_step_with_marine_fines(
            &f.grid,
            &s,
            &r,
            &f.elevation,
            &f.sediment,
            &f.ocean,
            0.0,
            1.0,
            &params(),
            Some(&f.available),
            &f.continental,
            &retention(),
            None,
            false,
            None,
            setting,
        )
        .unwrap();
        assert_eq!(old.0, new.0);
        assert_eq!(old.1, new.1);
        assert_eq!(old.2, new.2);
    }
}

#[test]
fn native_bypass_reaches_abyss_without_waiting_for_shelf_fill_and_closes() {
    let f = Fixture::new(false);
    let r = route(
        &f.grid,
        &f.elevation,
        &vec![1.0; f.grid.cells()],
        &f.ocean,
        0.0,
    )
    .unwrap();
    let s = DiffusionStencil::build(&f.grid).unwrap();
    let old = erode_retain_step_with_lakes(
        &f.grid,
        &s,
        &r,
        &f.elevation,
        &f.sediment,
        &f.ocean,
        0.0,
        1.0,
        &params(),
        Some(&f.available),
        &f.continental,
        &retention(),
        None,
        false,
        None,
    )
    .unwrap();
    let p = MarineFinesParams {
        fraction: 0.25,
        min_depth_m: 2500.0,
    };
    let new = erode_retain_step_with_marine_fines(
        &f.grid,
        &s,
        &r,
        &f.elevation,
        &f.sediment,
        &f.ocean,
        0.0,
        1.0,
        &params(),
        Some(&f.available),
        &f.continental,
        &retention(),
        Some(TraceClock {
            step: 1,
            start_my: 0.0,
            end_my: 1.0,
        }),
        true,
        None,
        Some(&p),
    )
    .unwrap();
    assert_closed(&old.0, &s, &f.available);
    assert_closed(&new.0, &s, &f.available);
    assert!(new.0.eroded_m3 > 0.0);
    assert_eq!(new.0.bedrock_removed_m, old.0.bedrock_removed_m);
    assert_receipt(&new.3);
    assert!(new.3.deposited_m3 > 0.0 && new.3.deposited_m3 <= new.0.eroded_m3);
    let deep: Vec<usize> = (0..f.grid.cells())
        .filter(|&c| f.elevation[c] < -2500.0)
        .collect();
    assert!(
        deep.iter().all(|&c| old.0.sediment_added_m[c] == 0.0),
        "fixture must leave old abyss bare"
    );
    assert!(deep.iter().all(|&c| new.0.sediment_added_m[c] > 0.0));
    // Well-mixed deposition is thickness-uniform, not equal volume per latitude.
    let h = new.0.sediment_added_m[deep[0]];
    assert!(deep
        .iter()
        .all(|&c| (new.0.sediment_added_m[c] - h).abs() < 1e-9 * h.max(1.0)));
    assert!(
        new.2.is_some(),
        "native trace/stage budget must survive the new path"
    );
}

#[test]
fn disconnected_unsupplied_ocean_receives_no_fines() {
    let f = Fixture::new(true);
    let r = route(
        &f.grid,
        &f.elevation,
        &vec![1.0; f.grid.cells()],
        &f.ocean,
        0.0,
    )
    .unwrap();
    let s = DiffusionStencil::build(&f.grid).unwrap();
    let p = MarineFinesParams {
        fraction: 0.25,
        min_depth_m: 2500.0,
    };
    let new = erode_retain_step_with_marine_fines(
        &f.grid,
        &s,
        &r,
        &f.elevation,
        &f.sediment,
        &f.ocean,
        0.0,
        1.0,
        &params(),
        Some(&f.available),
        &f.continental,
        &retention(),
        None,
        false,
        None,
        Some(&p),
    )
    .unwrap();
    assert_closed(&new.0, &s, &f.available);
    assert_receipt(&new.3);
    assert!(new
        .3
        .components
        .iter()
        .flat_map(|c| &c.receivers)
        .all(|r| r.cell < 3 * 24));
    assert!((0..3 * 24).all(|c| new.0.sediment_added_m[c] > 0.0));
    assert!(
        (13 * 24..16 * 24).all(|c| new.0.sediment_added_m[c] == 0.0),
        "no source or mixing across the land barrier"
    );
}

#[test]
fn no_deep_receiver_falls_back_to_exact_old_placement() {
    let f = Fixture::new(false);
    let r = route(
        &f.grid,
        &f.elevation,
        &vec![1.0; f.grid.cells()],
        &f.ocean,
        0.0,
    )
    .unwrap();
    let s = DiffusionStencil::build(&f.grid).unwrap();
    let old = erode_retain_step_with_lakes(
        &f.grid,
        &s,
        &r,
        &f.elevation,
        &f.sediment,
        &f.ocean,
        0.0,
        1.0,
        &params(),
        Some(&f.available),
        &f.continental,
        &retention(),
        None,
        false,
        None,
    )
    .unwrap();
    let p = MarineFinesParams {
        fraction: 0.25,
        min_depth_m: 5000.0,
    };
    let new = erode_retain_step_with_marine_fines(
        &f.grid,
        &s,
        &r,
        &f.elevation,
        &f.sediment,
        &f.ocean,
        0.0,
        1.0,
        &params(),
        Some(&f.available),
        &f.continental,
        &retention(),
        None,
        false,
        None,
        Some(&p),
    )
    .unwrap();
    assert_eq!(old.0, new.0);
    assert_eq!(old.1, new.1);
    assert_receipt(&new.3);
    assert_eq!(new.3.deposited_m3, 0.0);
    assert!(new.3.returned_to_shelf_m3 > 0.0);
}

#[test]
fn absent_erodible_source_creates_no_deposit() {
    let mut f = Fixture::new(false);
    f.available.fill(0.0);
    let r = route(
        &f.grid,
        &f.elevation,
        &vec![1.0; f.grid.cells()],
        &f.ocean,
        0.0,
    )
    .unwrap();
    let s = DiffusionStencil::build(&f.grid).unwrap();
    let p = MarineFinesParams {
        fraction: 0.25,
        min_depth_m: 2500.0,
    };
    let new = erode_retain_step_with_marine_fines(
        &f.grid,
        &s,
        &r,
        &f.elevation,
        &f.sediment,
        &f.ocean,
        0.0,
        1.0,
        &params(),
        Some(&f.available),
        &f.continental,
        &retention(),
        None,
        false,
        None,
        Some(&p),
    )
    .unwrap();
    assert_eq!(new.0.eroded_m3, 0.0);
    assert!(new.0.sediment_added_m.iter().all(|&h| h == 0.0));
    assert_receipt(&new.3);
    assert_eq!(new.3.deposited_m3, 0.0);
}

#[test]
fn accommodation_limit_returns_unsettled_mass_to_original_shelf_path() {
    let mut f = Fixture::new(false);
    for c in 0..f.grid.cells() {
        if f.ocean[c] {
            f.elevation[c] = -600.0;
        } else {
            f.available[c] = 10_000.0;
        }
    }
    let r = route(
        &f.grid,
        &f.elevation,
        &vec![1.0; f.grid.cells()],
        &f.ocean,
        0.0,
    )
    .unwrap();
    let s = DiffusionStencil::build(&f.grid).unwrap();
    let p = MarineFinesParams {
        fraction: 1.0,
        min_depth_m: 599.99,
    };
    let new = erode_retain_step_with_marine_fines(
        &f.grid,
        &s,
        &r,
        &f.elevation,
        &f.sediment,
        &f.ocean,
        0.0,
        1.0,
        &params(),
        Some(&f.available),
        &f.continental,
        &retention(),
        None,
        false,
        None,
        Some(&p),
    )
    .unwrap();
    assert_closed(&new.0, &s, &f.available);
    assert_receipt(&new.3);
    assert!(new.3.deposited_m3 > 0.0);
    assert!(
        new.3.returned_to_shelf_m3 > 0.0,
        "fixture must exhaust fine-receiver capacity"
    );
    for receiver in new.3.components.iter().flat_map(|c| &c.receivers) {
        let h = receiver.volume_m3 / s.areas_m2[receiver.cell];
        let surface = COLUMN.surface_after(f.elevation[receiver.cell], h, 0.0);
        assert!(surface <= -p.min_depth_m + 1e-7);
    }
}

#[test]
fn invalid_fine_parameters_are_rejected() {
    let f = Fixture::new(false);
    let r = route(
        &f.grid,
        &f.elevation,
        &vec![1.0; f.grid.cells()],
        &f.ocean,
        0.0,
    )
    .unwrap();
    let s = DiffusionStencil::build(&f.grid).unwrap();
    for (fraction, depth) in [
        (-0.1, 2500.0),
        (1.1, 2500.0),
        (f64::NAN, 2500.0),
        (0.2, -1.0),
        (0.2, f64::INFINITY),
    ] {
        let p = MarineFinesParams {
            fraction,
            min_depth_m: depth,
        };
        assert!(erode_retain_step_with_marine_fines(
            &f.grid,
            &s,
            &r,
            &f.elevation,
            &f.sediment,
            &f.ocean,
            0.0,
            1.0,
            &params(),
            Some(&f.available),
            &f.continental,
            &retention(),
            None,
            false,
            None,
            Some(&p)
        )
        .is_err());
    }
}
