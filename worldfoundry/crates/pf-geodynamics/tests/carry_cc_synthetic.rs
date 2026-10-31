//! SEED-CARRY-CC (lane CC, 2026-10-12): synthetic checks of
//! `front-rescue-carry-v1` against the legacy semi-Lagrangian carry
//! (pre-registered in `generated-test-data/FALSIFIER.txt`, C1-C3).
//!
//! Two all-continental hemispheres (35 km crust) rotate in opposite senses.
//! - HEAD-ON: both rotate about a common axis g and are split by a plane
//!   containing g, so they converge head-on along one half of the boundary
//!   great circle and diverge along the other. The exact convergent area
//!   closure is 4 r R^2 per My (r in rad/My, each plate at r).
//!   Geometries: `polar` (g = the polar axis: lane OG's degenerate lockstep
//!   case) and `tilted` (g at lat 60, lon 30: the boundary is oblique to the
//!   raster).
//! - TRANSFORM: the plates are the two hemispheres of a great circle G and
//!   each rotates about G's own pole, so they slide past each other along G
//!   (no convergence). `tilted` only (G is a raster staircase).
//!
//! The ignored driver prints one row per (case, grid, law, dt); the normal
//! test pins the conservation identities and a coarse area check.

use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_geodynamics::carry::CarryLaw;
use pf_geodynamics::kinematic::{advance_with_carry, CrustState};
use pf_geodynamics::orogen::BeltParams;
use pf_geodynamics::rotation::{lat_lon_to_unit, unit_to_lat_lon};
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;

const R: f64 = 6.371e6;
const RATE: f64 = 0.5; // deg/My per plate

fn belt() -> BeltParams {
    BeltParams {
        shear_zone_dip_deg: 30.0,
        channel_viscosity_pa_s: 1e19,
        channel_reference_thickness_m: 35_000.0,
        delamination_thickness_m: 70_000.0,
        arc_productivity_m: 500.0,
        arc_belt_width_m: 100_000.0,
        gravity_m_s2: 9.81,
        crust_density_kg_m3: 2800.0,
        mantle_density_kg_m3: 3300.0,
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    [a[0] / l, a[1] / l, a[2] / l]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// (rotation axis, split-plane normal) for a case.
fn geometry(case: &str) -> ([f64; 3], [f64; 3]) {
    match case {
        "headon_polar" => ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        "headon_tilted" => {
            let g = lat_lon_to_unit(60.0, 30.0).unwrap();
            (g, unit(cross(g, [0.0, 0.0, 1.0])))
        }
        "transform_tilted" => {
            let h = lat_lon_to_unit(60.0, 30.0).unwrap();
            (h, h)
        }
        other => panic!("unknown case {other}"),
    }
}

#[derive(Debug, Default, Clone)]
struct Stats {
    consumed_area_per_my: f64,
    consumed_m3: f64,
    belt_in_m3: f64,
    belt_placed_m3: f64,
    dropped_area_per_my: f64,
    cc_overlap_losers: u64,
    cc_losers_carried: u64,
    rescues: u64,
    edges_per_step: f64,
    ridge_area_per_my: f64,
    sliver_births: u64,
    resid_m3: f64,
    volume_change_m3: f64,
    delaminated_m3: f64,
    v0_m3: f64,
    cont_area_loss_per_my: f64,
}

fn run(case: &str, nlat: usize, nlon: usize, law: CarryLaw, dt: f64, duration: f64) -> Stats {
    let grid = GridSpec::new(nlat, nlon, R).unwrap();
    let n = grid.cells();
    let prior = UpliftPrior::preview_v0();
    let areas: Vec<f64> = (0..n)
        .map(|c| grid.cell_area_m2(c / nlon).unwrap())
        .collect();
    let (axis, m) = geometry(case);
    let (alat, alon) = unit_to_lat_lon(axis).unwrap();
    let rec = |id: u32, r: f64| PlateRecordV0 {
        plate_id: PlateId(id),
        parent_plate: None,
        birth_t_my: 0.0,
        euler_pole: Some(EulerPoleV0 {
            lat_deg: alat,
            lon_deg: alon,
            rate_deg_per_my: r,
        }),
        retirement: None,
    };
    let plates = vec![rec(0, RATE), rec(1, -RATE)];
    let cents: Vec<[f64; 3]> = (0..n)
        .map(|c| {
            lat_lon_to_unit(
                grid.center_lat_deg(c / nlon).unwrap(),
                grid.center_lon_deg(c % nlon).unwrap(),
            )
            .unwrap()
        })
        .collect();
    let mut st = CrustState {
        t_my: 0.0,
        plate: cents
            .iter()
            .map(|p| PlateId(if dot(*p, m) >= 0.0 { 0 } else { 1 }))
            .collect(),
        kind: vec![CrustKindCodeV0::Continental; n],
        birth_t_my: vec![-2000.0; n],
        thickness_m: vec![35_000.0; n],
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: cents.clone(),
    };
    let b = belt();
    let steps = (duration / dt).round() as usize;
    let mut s = Stats {
        v0_m3: st.thickness_m.iter().zip(&areas).map(|(h, a)| h * a).sum(),
        ..Default::default()
    };
    let cont_area = |st: &CrustState| -> f64 {
        (0..n)
            .filter(|&c| st.kind[c] == CrustKindCodeV0::Continental)
            .map(|c| areas[c])
            .sum()
    };
    let a0 = cont_area(&st);
    let mut edges = 0u64;
    for _ in 0..steps {
        let (next, l, _, _) =
            advance_with_carry(&grid, &st, &plates, &prior, dt, Some(&b), law).unwrap();
        s.consumed_m3 += l.carry.cc_consumed_m3;
        s.consumed_area_per_my += l.carry.cc_consumed_area_m2;
        s.belt_in_m3 += l.belt_in_m3;
        s.belt_placed_m3 += l.belt_placed_m3;
        s.dropped_area_per_my += l.carry.continental_dropped_area_m2;
        s.cc_overlap_losers += l.carry.cc_overlap_losers;
        s.cc_losers_carried += l.carry.cc_losers_carried;
        s.rescues += l.carry.rescues_admitted;
        s.ridge_area_per_my += l.ridge_created_m3 / prior.record.oceanic_crust_thickness_m;
        s.sliver_births += l.carry.sliver_births;
        s.resid_m3 += l.residual_m3();
        s.delaminated_m3 += l.delaminated_m3;
        edges += l.collision_supply_edges;
        st = next;
    }
    s.edges_per_step = edges as f64 / steps as f64;
    s.consumed_area_per_my /= duration;
    s.dropped_area_per_my /= duration;
    s.ridge_area_per_my /= duration;
    s.cont_area_loss_per_my = (a0 - cont_area(&st)) / duration;
    let v1: f64 = st.thickness_m.iter().zip(&areas).map(|(h, a)| h * a).sum();
    s.volume_change_m3 = v1 - s.v0_m3;
    s
}

fn exact_area_per_my() -> f64 {
    4.0 * RATE.to_radians() * R * R
}

#[test]
fn front_rescue_head_on_consumes_and_conserves_volume() {
    // Regression pin (36x72, 5-degree cells, 40 My, dt 1), NOT the
    // pre-registered C1 band: C1 (+-25% at 36x72, +-15% at 72x144) FIRED
    // for this law (measured 0.61-1.00 of the exact closure, falling with
    // dt; FALSIFIER.txt). Pinned here: the belt receives exactly the
    // consumed volume, the consumed area stays above half the exact
    // closure, and the legacy carry consumes nothing (it drops the front
    // parcels silently instead).
    let exact = exact_area_per_my();
    for case in ["headon_polar", "headon_tilted"] {
        let s = run(case, 36, 72, CarryLaw::FrontRescueV1, 1.0, 40.0);
        let rel = s.consumed_area_per_my / exact;
        assert!((0.5..=1.25).contains(&rel), "{case}: consumed/exact {rel}");
        assert!(
            (s.belt_in_m3 - s.consumed_m3).abs() <= 1e-9 * s.consumed_m3,
            "{case}: belt input {} vs consumed {}",
            s.belt_in_m3,
            s.consumed_m3
        );
        assert!(
            (s.belt_placed_m3 - s.belt_in_m3).abs() <= 1e-9 * s.belt_in_m3,
            "{case}: placed {} vs input {}",
            s.belt_placed_m3,
            s.belt_in_m3
        );
        let legacy = run(
            case,
            36,
            72,
            CarryLaw::SemiLagrangianNearestMaterialV0,
            1.0,
            40.0,
        );
        assert_eq!(legacy.consumed_m3, 0.0);
        assert_eq!(legacy.rescues, 0);
    }
}

#[test]
fn front_rescue_requires_the_belt_law() {
    let grid = GridSpec::new(18, 36, R).unwrap();
    let n = grid.cells();
    let st = CrustState {
        t_my: 0.0,
        plate: vec![PlateId(0); n],
        kind: vec![CrustKindCodeV0::Continental; n],
        birth_t_my: vec![-2000.0; n],
        thickness_m: vec![35_000.0; n],
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![0; n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: (0..n)
            .map(|c| {
                lat_lon_to_unit(
                    grid.center_lat_deg(c / 36).unwrap(),
                    grid.center_lon_deg(c % 36).unwrap(),
                )
                .unwrap()
            })
            .collect(),
    };
    let plates = vec![PlateRecordV0 {
        plate_id: PlateId(0),
        parent_plate: None,
        birth_t_my: 0.0,
        euler_pole: Some(EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: 0.5,
        }),
        retirement: None,
    }];
    let prior = UpliftPrior::preview_v0();
    assert!(advance_with_carry(
        &grid,
        &st,
        &plates,
        &prior,
        1.0,
        None,
        CarryLaw::FrontRescueV1
    )
    .is_err());
    // Determinism and order: the same inputs give the same bytes.
    let b = belt();
    let a = advance_with_carry(
        &grid,
        &st,
        &plates,
        &prior,
        1.0,
        Some(&b),
        CarryLaw::FrontRescueV1,
    )
    .unwrap();
    let c = advance_with_carry(
        &grid,
        &st,
        &plates,
        &prior,
        1.0,
        Some(&b),
        CarryLaw::FrontRescueV1,
    )
    .unwrap();
    assert_eq!(a.0, c.0);
    assert_eq!(a.3, c.3);
}

#[test]
#[ignore]
fn carry_cc_synthetic_table() {
    let exact = exact_area_per_my();
    let dur: f64 = std::env::var("CC_DURATION")
        .map(|v| v.parse().unwrap())
        .unwrap_or(40.0);
    println!(
        "exact convergent closure {exact:.4e} m2/My; volume {:.4e} m3/My at 35 km; duration {dur} My",
        exact * 35_000.0
    );
    println!("case grid law dt consumed_area/My rel_exact consumed_m3 belt_in_m3 belt_placed_m3 dropped_area/My cc_losers losers_carried rescues edges/step ridge_area/My cont_area_loss/My sliver_births resid_m3 dV_m3 delam_m3");
    for case in ["headon_polar", "headon_tilted", "transform_tilted"] {
        for (nlat, nlon) in [(36usize, 72usize), (72, 144)] {
            for (lname, law) in [
                ("legacy", CarryLaw::SemiLagrangianNearestMaterialV0),
                ("rescue_v1", CarryLaw::FrontRescueV1),
            ] {
                for dt in [2.0, 1.0, 0.5, 0.25] {
                    let s = run(case, nlat, nlon, law, dt, dur);
                    println!(
                        "{case} {nlat}x{nlon} {lname} {dt} {:.4e} {:.3} {:.4e} {:.4e} {:.4e} {:.4e} {} {} {} {:.1} {:.4e} {:.4e} {} {:.4e} {:.4e} {:.4e}",
                        s.consumed_area_per_my,
                        s.consumed_area_per_my / exact,
                        s.consumed_m3,
                        s.belt_in_m3,
                        s.belt_placed_m3,
                        s.dropped_area_per_my,
                        s.cc_overlap_losers,
                        s.cc_losers_carried,
                        s.rescues,
                        s.edges_per_step,
                        s.ridge_area_per_my,
                        s.cont_area_loss_per_my,
                        s.sliver_births,
                        s.resid_m3,
                        s.volume_change_m3,
                        s.delaminated_m3
                    );
                }
            }
        }
    }
}

/// Diagnostic: where the convergent closure goes. Per run, continental
/// parcel area carried nowhere (drops) and extra copies (dups), split into
/// the contact zone (step-start cells within 2 rings of a different-plate
/// cell) on the convergent and divergent halves, plus consumed area.
#[test]
#[ignore]
fn carry_cc_budget_diag() {
    let exact = exact_area_per_my();
    let case = std::env::var("CC_CASE").unwrap_or("headon_tilted".into());
    let (nlat, nlon) = (72usize, 144usize);
    for (lname, law) in [
        ("legacy", CarryLaw::SemiLagrangianNearestMaterialV0),
        ("rescue_v1", CarryLaw::FrontRescueV1),
    ] {
        for dt in [2.0, 1.0, 0.5, 0.25] {
            let grid = GridSpec::new(nlat, nlon, R).unwrap();
            let n = grid.cells();
            let prior = UpliftPrior::preview_v0();
            let areas: Vec<f64> = (0..n)
                .map(|c| grid.cell_area_m2(c / nlon).unwrap())
                .collect();
            let (axis, m) = geometry(&case);
            let conv_dir = unit(cross(axis, m));
            let (alat, alon) = unit_to_lat_lon(axis).unwrap();
            let rec = |id: u32, r: f64| PlateRecordV0 {
                plate_id: PlateId(id),
                parent_plate: None,
                birth_t_my: 0.0,
                euler_pole: Some(EulerPoleV0 {
                    lat_deg: alat,
                    lon_deg: alon,
                    rate_deg_per_my: r,
                }),
                retirement: None,
            };
            let plates = vec![rec(0, RATE), rec(1, -RATE)];
            let cents: Vec<[f64; 3]> = (0..n)
                .map(|c| {
                    lat_lon_to_unit(
                        grid.center_lat_deg(c / nlon).unwrap(),
                        grid.center_lon_deg(c % nlon).unwrap(),
                    )
                    .unwrap()
                })
                .collect();
            let mut st = CrustState {
                t_my: 0.0,
                plate: cents
                    .iter()
                    .map(|p| PlateId(if dot(*p, m) >= 0.0 { 0 } else { 1 }))
                    .collect(),
                kind: vec![CrustKindCodeV0::Continental; n],
                birth_t_my: vec![-2000.0; n],
                thickness_m: vec![35_000.0; n],
                sediment_m: vec![0.0; n],
                eroded_total_m: vec![0.0; n],
                arc_added_m: vec![0.0; n],
                shortening_added_m: vec![0.0; n],
                boundary_class: vec![0; n],
                normal_rate_m_per_my: vec![0.0; n],
                material_pos: cents.clone(),
            };
            let b = belt();
            let steps = (40.0 / dt) as usize;
            // [conv drops, conv dups, div drops, div dups, far drops, far dups]
            let mut acc = [0.0f64; 6];
            let mut consumed = 0.0;
            for _ in 0..steps {
                // zone: 2 rings
                let mut zone = vec![false; n];
                for c in 0..n {
                    let (i, j) = grid.ij(c as u64).unwrap();
                    for (ti, tj, _) in grid.neighbors(i, j).unwrap() {
                        let k = ti * nlon + tj;
                        if st.plate[k] != st.plate[c] {
                            zone[c] = true;
                        }
                    }
                }
                let z1 = zone.clone();
                for c in 0..n {
                    if z1[c] {
                        let (i, j) = grid.ij(c as u64).unwrap();
                        for (ti, tj, _) in grid.neighbors(i, j).unwrap() {
                            zone[ti * nlon + tj] = true;
                        }
                    }
                }
                let (next, l, _, src) =
                    advance_with_carry(&grid, &st, &plates, &prior, dt, Some(&b), law).unwrap();
                consumed += l.carry.cc_consumed_area_m2;
                let mut k = vec![0u32; n];
                for &s in &src {
                    if s != u32::MAX {
                        k[s as usize] += 1;
                    }
                }
                for s in 0..n {
                    let slot = if zone[s] {
                        if dot(cents[s], conv_dir) > 0.0 {
                            0
                        } else {
                            2
                        }
                    } else {
                        4
                    };
                    if k[s] == 0 {
                        acc[slot] += areas[s];
                    } else if k[s] > 1 {
                        acc[slot + 1] += (k[s] - 1) as f64 * areas[s];
                    }
                }
                st = next;
            }
            let f = |x: f64| x / 40.0 / exact;
            println!(
                "{case} {lname} dt {dt}: per My / exact: consumed {:.3} | conv drops {:.3} dups {:.3} net {:.3} | div drops {:.3} dups {:.3} | far drops {:.3} dups {:.3}",
                f(consumed), f(acc[0]), f(acc[1]), f(acc[0] - acc[1]), f(acc[2]), f(acc[3]), f(acc[4]), f(acc[5])
            );
        }
    }
}
