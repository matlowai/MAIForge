//! SEED-EXTRAS (lane EX, 2026-10-12): synthetic checks of the plume and
//! rift passes behind `advance_with_extras`.
//!
//! - With no extras selected, `advance_with_extras` is `advance_with_laws`
//!   bit for bit.
//! - `stretched-rift-v1` on two all-continental hemispheres pulled apart:
//!   rift births are continental, the crust they hold is exactly the flank
//!   volume moved (one step against the legacy step from the same state),
//!   and continued opening reaches breakup (oceanic births).
//! - `fixed-plume-chains-v1` on one rotating oceanic plate: the step's tail
//!   volume and a head arriving in the step are placed exactly, and the
//!   tail track is age-progressive (older eruptions farther downstream).

use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::params::SeedPath;
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_geodynamics::carry::CarryLaw;
use pf_geodynamics::hotspots::{seed_plumes, Plume, PlumeParams};
use pf_geodynamics::kinematic::{advance_with_extras, advance_with_laws, CrustState, ExtrasLaws};
use pf_geodynamics::orogen::{BeltParams, CollisionSedimentLaw};
use pf_geodynamics::rift::RiftParams;
use pf_geodynamics::rotation::{great_circle_rad, lat_lon_to_unit, unit_to_lat_lon};
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;

const R: f64 = 6.371e6;

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

fn cents(grid: &GridSpec) -> Vec<[f64; 3]> {
    let nlon = grid.nlon();
    (0..grid.cells())
        .map(|c| {
            lat_lon_to_unit(
                grid.center_lat_deg(c / nlon).unwrap(),
                grid.center_lon_deg(c % nlon).unwrap(),
            )
            .unwrap()
        })
        .collect()
}

fn state(
    grid: &GridSpec,
    plate: Vec<PlateId>,
    kind: CrustKindCodeV0,
    h: f64,
) -> CrustState {
    let n = grid.cells();
    CrustState {
        t_my: 0.0,
        plate,
        kind: vec![kind; n],
        birth_t_my: vec![-200.0; n],
        thickness_m: vec![h; n],
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: cents(grid),
    }
}

/// Two all-continental hemispheres split by a tilted great circle, rotating
/// in opposite senses about an axis in that circle's plane, so they open.
fn diverging() -> (GridSpec, Vec<PlateRecordV0>, CrustState, Vec<f64>) {
    let (nlat, nlon) = (36, 72);
    let grid = GridSpec::new(nlat, nlon, R).unwrap();
    let n = grid.cells();
    let areas: Vec<f64> = (0..n).map(|c| grid.cell_area_m2(c / nlon).unwrap()).collect();
    let g = lat_lon_to_unit(60.0, 30.0).unwrap();
    let z = [0.0, 0.0, 1.0];
    let m = [g[1] * z[2] - g[2] * z[1], g[2] * z[0] - g[0] * z[2], g[0] * z[1] - g[1] * z[0]];
    let (alat, alon) = unit_to_lat_lon(g).unwrap();
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
    // Lane CC's head-on synthetic with the senses swapped (opening).
    let plates = vec![rec(0, -1.0), rec(1, 1.0)];
    let plate = cents(&grid)
        .iter()
        .map(|p| PlateId(if p[0] * m[0] + p[1] * m[1] + p[2] * m[2] >= 0.0 { 0 } else { 1 }))
        .collect();
    let st = state(&grid, plate, CrustKindCodeV0::Continental, 35_000.0);
    (grid, plates, st, areas)
}

fn vol(v: &[f64], areas: &[f64]) -> f64 {
    v.iter().zip(areas).map(|(h, a)| h * a).sum()
}

#[test]
fn no_extras_is_advance_with_laws_bit_for_bit() {
    let (grid, plates, st0, _) = diverging();
    let prior = UpliftPrior::preview_v0();
    let b = belt();
    let (mut a, mut x) = (st0.clone(), st0);
    for _ in 0..6 {
        let (na, la, ta, sa) = advance_with_laws(
            &grid,
            &a,
            &plates,
            &prior,
            1.0,
            Some(&b),
            CarryLaw::FrontRescueV1,
            CollisionSedimentLaw::AccretedToBeltCrustV1,
        )
        .unwrap();
        let (nx, lx, tx, sx) = advance_with_extras(
            &grid,
            &x,
            &plates,
            &prior,
            1.0,
            Some(&b),
            CarryLaw::FrontRescueV1,
            CollisionSedimentLaw::AccretedToBeltCrustV1,
            &ExtrasLaws::default(),
        )
        .unwrap();
        assert_eq!(na, nx);
        assert_eq!(la, lx);
        assert_eq!(ta, tx);
        assert_eq!(sa, sx);
        a = na;
        x = nx;
    }
    // Extras without the belt are refused.
    let rift = ExtrasLaws {
        plumes: None,
        rift: Some(RiftParams {
            half_width_m: 100_000.0,
            breakup_thickness_m: 10_000.0,
        }),
        flood_basalt: None, // SEED-FLOODBASALT literal 1
    };
    assert!(advance_with_extras(
        &grid,
        &a,
        &plates,
        &prior,
        1.0,
        None,
        CarryLaw::SemiLagrangianNearestMaterialV0,
        CollisionSedimentLaw::SurfaceLayerV0,
        &rift
    )
    .is_err());
}

#[test]
fn rift_births_are_stretched_continent_and_reach_breakup() {
    let (grid, plates, st0, areas) = diverging();
    let prior = UpliftPrior::preview_v0();
    let b = belt();
    let rp = RiftParams {
        half_width_m: 100_000.0,
        breakup_thickness_m: 10_000.0,
    };
    let x = ExtrasLaws {
        plumes: None,
        rift: Some(rp),
        flood_basalt: None, // SEED-FLOODBASALT literal 2
    };
    let step = |st: &CrustState, ex: &ExtrasLaws| {
        advance_with_extras(
            &grid,
            st,
            &plates,
            &prior,
            1.0,
            Some(&b),
            CarryLaw::FrontRescueV1,
            CollisionSedimentLaw::AccretedToBeltCrustV1,
            ex,
        )
        .unwrap()
    };
    // One step from the same state: the same carry, so the only difference
    // is the conversion of rift births. Total crust differs by exactly the
    // oceanic crust the rift births no longer receive.
    let mut st = st0.clone();
    let mut found = false;
    for _ in 0..10 {
        let (lg, ll, _, _) = step(&st, &ExtrasLaws::default());
        let (rf, rl, _, _) = step(&st, &x);
        if !rl.rift_births.is_empty() {
            found = true;
            let ocean: f64 = rl
                .rift_births
                .iter()
                .map(|&c| prior.record.oceanic_crust_thickness_m * areas[c as usize])
                .sum();
            let d = vol(&rf.thickness_m, &areas) - (vol(&lg.thickness_m, &areas) - ocean);
            assert!(d.abs() <= 1e-9 * vol(&lg.thickness_m, &areas), "volume drift {d}");
            let moved: f64 = rl.rift_births.iter().map(|&c| rf.thickness_m[c as usize] * areas[c as usize]).sum();
            assert!((moved - rl.rift_moved_m3).abs() <= 1e-9 * moved);
            let thin: f64 = rl.rift_thinning.iter().map(|&(c, dh)| dh * areas[c as usize]).sum();
            assert!((thin - moved).abs() <= 1e-9 * moved);
            for &c in &rl.rift_births {
                assert_eq!(rf.kind[c as usize], CrustKindCodeV0::Continental);
                assert!(rf.thickness_m[c as usize] >= rp.breakup_thickness_m);
                assert!(rf.thickness_m[c as usize] < 35_000.0);
            }
            assert_eq!(ll.ridge_cells, rl.ridge_cells + rl.rift_births.len() as u64);
            // The step's residual is unchanged by the conversion (to rounding).
            assert!((ll.residual_m3() - rl.residual_m3()).abs() <= 1e-9 * ll.volume_before_m3);
        }
        st = rf;
    }
    assert!(found, "no rift birth in 10 My of opening");
    // Keep opening: the rift thins until births break up into ocean.
    let mut breakup = 0;
    for _ in 0..60 {
        let (rf, rl, _, _) = step(&st, &x);
        breakup += rl.rift_breakup_births;
        st = rf;
    }
    assert!(breakup > 0, "no breakup in 70 My at 2 deg/My relative");
    assert!(st.kind.iter().any(|k| *k == CrustKindCodeV0::Oceanic));
}

#[test]
fn plume_tail_is_age_progressive_and_volumes_are_exact() {
    let (nlat, nlon) = (72, 144);
    let grid = GridSpec::new(nlat, nlon, R).unwrap();
    let n = grid.cells();
    let areas: Vec<f64> = (0..n).map(|c| grid.cell_area_m2(c / nlon).unwrap()).collect();
    let plates = vec![PlateRecordV0 {
        plate_id: PlateId(0),
        parent_plate: None,
        birth_t_my: 0.0,
        euler_pole: Some(EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: 1.0,
        }),
        retirement: None,
    }];
    let st = state(&grid, vec![PlateId(0); n], CrustKindCodeV0::Oceanic, 7_000.0);
    let p = lat_lon_to_unit(1.25, 0.0).unwrap();
    let pp = PlumeParams {
        head_volume_m3: 2.0e15,
        head_radius_m: 500_000.0,
        tail_flux_m3_per_my: 1.0e14,
        edifice_radius_m: 100_000.0,
    };
    let prior = UpliftPrior::preview_v0();
    let b = belt();
    let run = |plumes: &[Plume]| {
        advance_with_extras(
            &grid,
            &st,
            &plates,
            &prior,
            8.0,
            Some(&b),
            CarryLaw::FrontRescueV1,
            CollisionSedimentLaw::AccretedToBeltCrustV1,
            &ExtrasLaws {
                plumes: Some((plumes, pp)),
                rift: None,
                flood_basalt: None, // SEED-FLOODBASALT literal 3
            },
        )
        .unwrap()
    };
    // Tail only (onset before the step).
    let tail = [Plume {
        id: 0,
        position: p,
        onset_my: -10.0,
        split_azimuth_rad: 0.0,
    }];
    let (_, l, _, _) = run(&tail);
    assert_eq!(l.hotspot_heads, 0);
    assert!((l.hotspot_created_m3 - 8.0e14).abs() <= 1e-9 * 8.0e14, "{}", l.hotspot_created_m3);
    let placed: f64 = l.hotspot_cells.iter().map(|&(c, dh, _)| dh * areas[c as usize]).sum();
    assert!((placed - l.hotspot_created_m3).abs() <= 1e-9 * placed);
    // Age progression: eastward plate motion (positive rate about +z) carries
    // older eruptions farther east of the plume.
    let mut by_t: Vec<(f64, f64)> = l
        .hotspot_cells
        .iter()
        .map(|&(c, _, t)| {
            let x = [
                grid.center_lat_deg(c as usize / nlon).unwrap(),
                grid.center_lon_deg(c as usize % nlon).unwrap(),
            ];
            let d = great_circle_rad(p, lat_lon_to_unit(x[0], x[1]).unwrap()).unwrap();
            (t, d)
        })
        .collect();
    by_t.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (oldest, youngest) = (by_t.first().unwrap(), by_t.last().unwrap());
    assert!(oldest.0 < youngest.0);
    assert!(oldest.1 > youngest.1, "oldest {oldest:?} youngest {youngest:?}");
    // A head arriving inside the step adds its volume once.
    let head = [Plume {
        id: 0,
        position: p,
        onset_my: 3.0,
        split_azimuth_rad: 0.0,
    }];
    let (_, l, _, _) = run(&head);
    assert_eq!(l.hotspot_heads, 1);
    assert!((l.hotspot_head_m3 - 2.0e15).abs() <= 1e-9 * 2.0e15);
    assert!((l.hotspot_created_m3 - (2.0e15 + 5.0e14)).abs() <= 1e-9 * 2.5e15);
    // The step's residual names the plume volume (a rigid single plate: the
    // carry itself is a permutation here, so the residual is rounding).
    assert!(l.residual_m3().abs() <= 1e-6 * l.volume_before_m3, "{}", l.residual_m3());
}

#[test]
fn plumes_seed_deterministically_with_onsets_in_range() {
    let sp = SeedPath {
        seed: 240907,
        derivation: "d".into(),
        salt: "s".into(),
    };
    let a = seed_plumes(12, 50.0, 100.0, &sp);
    assert_eq!(a, seed_plumes(12, 50.0, 100.0, &sp));
    assert_eq!(a.len(), 12);
    let k6 = pf_geodynamics::hotspots::seed_hotspots(12, &sp);
    for (p, h) in a.iter().zip(&k6) {
        assert_eq!(p.position, h.position);
        assert!(p.onset_my >= -50.0 && p.onset_my < 100.0);
        assert!(p.split_azimuth_rad >= 0.0 && p.split_azimuth_rad < std::f64::consts::TAU);
    }
}
