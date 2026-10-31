//! SEED-PLATE-REORG (lane PR, 2026-10-13): synthetic checks of
//! `plate-reorg-scripted-v1` (pre-registered in
//! `generated-test-data/PREREG.txt`).
//!
//! 1. A plume split of an all-continental plate follows an oblique weak belt
//!    (thick crust) instead of the plume's great circle.
//! 2. The split's bookkeeping: new ids on side B, parcels follow their cells,
//!    two suture groups, and side B opens at the rift rate at the plume.
//! 3. A scheduled re-solve gives every suture group one new pole (members
//!    share it) and is deterministic.

use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::params::SeedPath;
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_geodynamics::boundaries::relative_velocity_at;
use pf_geodynamics::conserve::{ParcelLayer, SutureState};
use pf_geodynamics::hotspots::Plume;
use pf_geodynamics::kinematic::{nearest_cell, CrustState};
use pf_geodynamics::reorg::{plume_split, reorg_pass, ReorgCtx, ReorgParams, ReorgState};
use pf_geodynamics::rotation::{lat_lon_to_unit, unit_to_lat_lon};
use pf_grid::GridSpec;

const R: f64 = 6.371e6;

fn params() -> ReorgParams {
    ReorgParams {
        rift_opening_rate_deg_per_my: 0.25,
        interval_min_my: 50.0,
        interval_max_my: 100.0,
        slab_pull_bias: 0.5,
        slab_circumference_reference: 0.3,
        rate_floor_fraction: 0.2,
        supercontinent_share: 0.6,
        supercontinent_tenure_my: 100.0,
        breakup_piece_max_share: 0.4,
        breakup_cascade_my: 20.0,
        rift_persistence_my: 100.0,
        weak_zone_strength_ratio: 0.2,
        plate_count_min: 1,
        plate_count_max: 12,
        plate_count_min_area_fraction: 0.005,
        subduction_initiation_rate_deg_per_my: 0.3,
        plume_head_interval_my: 25.0,
        plume_tail_lifetime_my: 150.0,
        plume_split_min_share: 0.2,
        rate_scale_deg_per_my: 1.0,
        plume_head_radius_m: 500_000.0,
    }
}

fn rec(id: u32, lat: f64, lon: f64, r: f64) -> PlateRecordV0 {
    PlateRecordV0 {
        plate_id: PlateId(id),
        parent_plate: None,
        birth_t_my: 0.0,
        euler_pole: Some(EulerPoleV0 { lat_deg: lat, lon_deg: lon, rate_deg_per_my: r }),
        retirement: None,
    }
}

/// Plate 0 = an all-continental block (lat -40..40, lon -50..50), plate 1 =
/// the rest (ocean). With `belt`, crust along lon = 0.4 lat + 12 is 10 km
/// thicker (a worn suture root / belt).
fn world(belt: bool) -> (GridSpec, CrustState) {
    let (nlat, nlon) = (72, 144);
    let grid = GridSpec::new(nlat, nlon, R).unwrap();
    let n = grid.cells();
    let mut plate = vec![PlateId(1); n];
    let mut kind = vec![CrustKindCodeV0::Oceanic; n];
    let mut th = vec![7_000.0; n];
    let mut cents = Vec::with_capacity(n);
    for c in 0..n {
        let lat = grid.center_lat_deg(c / nlon).unwrap();
        let lon = grid.center_lon_deg(c % nlon).unwrap();
        let lon = if lon > 180.0 { lon - 360.0 } else { lon };
        cents.push(lat_lon_to_unit(lat, lon).unwrap());
        if lat.abs() < 40.0 && lon.abs() < 50.0 {
            plate[c] = PlateId(0);
            kind[c] = CrustKindCodeV0::Continental;
            th[c] = 35_000.0;
            if belt && (lon - (0.4 * lat + 12.0)).abs() < 3.0 {
                th[c] = 45_000.0;
            }
        }
    }
    let st = CrustState {
        t_my: 0.0,
        plate,
        kind,
        birth_t_my: vec![-100.0; n],
        thickness_m: th,
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: cents,
    };
    (grid, st)
}

fn sp() -> SeedPath {
    SeedPath { seed: 42, derivation: "test".into(), salt: "plate-reorg".into() }
}

fn split(belt: bool) -> (GridSpec, CrustState, Vec<PlateRecordV0>, ParcelLayer, SutureState, PlateId, PlateId) {
    let (grid, mut st) = world(belt);
    let mut plates = vec![rec(0, 10.0, 20.0, 0.3), rec(1, -30.0, 100.0, 0.5)];
    let mut layer = ParcelLayer::seed(&grid, &st).unwrap();
    let mut sut = SutureState::default();
    let mut rs = ReorgState::default();
    let p = params();
    let ctx = ReorgCtx { grid: &grid, seed_path: &sp(), params: &p, schedule: &[] };
    // Plume at the block centre; azimuth 0 -> the cut runs north-south.
    let plume = Plume { id: 3, position: lat_lon_to_unit(0.0, 0.0).unwrap(), onset_my: 1.0, split_azimuth_rad: 0.0 };
    let (_, _, a, b, _) = plume_split(&ctx, &mut st, &mut plates, &mut layer, &mut sut, &mut rs, &plume, 1.0)
        .unwrap()
        .expect("the plume splits the block");
    (grid, st, plates, layer, sut, a, b)
}

/// Fractions of continental cells west / east of the belt line on side A.
fn belt_sides(belt: bool) -> (f64, f64) {
    let (grid, st, _, _, sut, a, _) = split(belt);
    let nlon = grid.nlon();
    // Continental cells well west / east of the belt line lie on opposite
    // sides (by suture group).
    let (mut west_a, mut west, mut east_a, mut east) = (0, 0, 0, 0);
    for c in 0..grid.cells() {
        if st.kind[c] != CrustKindCodeV0::Continental {
            continue;
        }
        let lat = grid.center_lat_deg(c / nlon).unwrap();
        let lon = grid.center_lon_deg(c % nlon).unwrap();
        let lon = if lon > 180.0 { lon - 360.0 } else { lon };
        let d = lon - (0.4 * lat + 12.0);
        let on_a = sut.find(st.plate[c]) == a;
        if d < -6.0 {
            west += 1;
            west_a += on_a as u32;
        } else if d > 6.0 {
            east += 1;
            east_a += on_a as u32;
        }
    }
    (west_a as f64 / west as f64, east_a as f64 / east as f64)
}

#[test]
fn plume_split_follows_the_weak_belt() {
    let (wf, ef) = belt_sides(true);
    assert!((wf > 0.9 && ef < 0.1) || (wf < 0.1 && ef > 0.9), "west on A {wf}, east on A {ef}");
    // Control: without the belt the cut stays near the plume's meridian, so
    // the belt-line partition does not hold.
    let (wf0, ef0) = belt_sides(false);
    assert!(!((wf0 > 0.9 && ef0 < 0.1) || (wf0 < 0.1 && ef0 > 0.9)), "control west {wf0} east {ef0}");
}

#[test]
fn split_bookkeeping_and_opening() {
    let (grid, st, plates, layer, sut, a, b) = split(false);
    assert_ne!(sut.find(a), sut.find(b));
    // Every parcel's plate is its nearest cell's plate.
    for k in 0..layer.len() {
        let c = nearest_cell(&grid, layer.pos[k]).unwrap();
        assert_eq!(layer.plate[k], st.plate[c], "parcel {k}");
    }
    // Both sides hold continent.
    let mut sides = [0u32; 2];
    for c in 0..grid.cells() {
        if st.kind[c] == CrustKindCodeV0::Continental {
            if sut.find(st.plate[c]) == sut.find(a) {
                sides[0] += 1
            } else if sut.find(st.plate[c]) == sut.find(b) {
                sides[1] += 1
            }
        }
    }
    assert!(sides[0] > 50 && sides[1] > 50, "{sides:?}");
    // B opens away from A at the plume at the rift rate (azimuth 0: n is
    // east or west).
    let pa = plates.iter().find(|r| r.plate_id == sut.find(a)).unwrap().euler_pole.unwrap();
    let pb = plates.iter().find(|r| r.plate_id == sut.find(b)).unwrap().euler_pole.unwrap();
    let p = lat_lon_to_unit(0.0, 0.0).unwrap();
    let v = relative_velocity_at(p, &pb, &pa).unwrap();
    let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    assert!((speed - 0.25f64.to_radians()).abs() < 1e-9, "{speed}");
    let (_, _) = unit_to_lat_lon(v.map(|x| x / speed)).unwrap();
}

#[test]
fn resolve_gives_groups_one_new_pole_and_is_deterministic() {
    let go = || {
        let (grid, mut st, mut plates, mut layer, mut sut, _, _) = split(true);
        let p = params();
        let schedule = [2.0];
        let ctx = ReorgCtx { grid: &grid, seed_path: &sp(), params: &p, schedule: &schedule };
        let mut rs = ReorgState::default();
        // Suture plate 1 (ocean) to the first group so a group has 2 members.
        let r0 = sut.find(PlateId(0));
        sut.parent.insert(PlateId(1), r0);
        let before: Vec<_> = plates.iter().map(|r| r.euler_pole).collect();
        let rep = reorg_pass(&ctx, &mut st, &mut plates, &mut layer, &mut sut, &mut rs, 1.0, 3.0).unwrap();
        assert_eq!(rep.resolves, 1);
        assert_eq!(rs.next_event, 1);
        let after: Vec<_> = plates.iter().map(|r| r.euler_pole).collect();
        assert!(before.iter().zip(&after).any(|(x, y)| x != y));
        for r in &plates {
            let root = sut.find(r.plate_id);
            let rp = plates.iter().find(|q| q.plate_id == root).unwrap().euler_pole;
            assert_eq!(r.euler_pole, rp, "plate {} shares its root's pole", r.plate_id.0);
        }
        (after, st.plate.clone())
    };
    assert_eq!(go(), go());
}
