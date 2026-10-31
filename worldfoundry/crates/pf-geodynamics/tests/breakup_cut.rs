//! SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13): synthetic checks of
//! `plate-reorg-scripted-v2` against v1 (pre-registered in
//! `generated-test-data/PREREG.txt`).
//!
//! 1. A scripted breakup of a continent with thin (rifted) margins, sutured to
//!    its ocean: the v2 cut crosses the continent (both sides hold a quarter
//!    or more of it) and the new boundary meets the coast only at its exit
//!    points (little boundary on continent/ocean edges). v1 is run on the
//!    same world as the reference (its coast share is printed and must be
//!    larger).
//! 2. A member plate folded by consolidation is retired as `Folded` within
//!    the same pass under v2; v1 leaves it to the next pass and `Subducted`.
//! 3. The v2 breakup is deterministic (run twice).

use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0, PlateRetirementV0};
use pf_contracts::params::SeedPath;
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_geodynamics::conserve::{ParcelLayer, SutureState};
use pf_geodynamics::kinematic::CrustState;
use pf_geodynamics::reorg::{reorg_pass_law, ReorgCtx, ReorgLaw, ReorgParams, ReorgState};
use pf_geodynamics::rotation::lat_lon_to_unit;
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
        plume_split_min_share: 0.35,
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

fn sp() -> SeedPath {
    SeedPath { seed: 42, derivation: "test".into(), salt: "breakup-cut".into() }
}

/// Plate 0 = a continent (lat -40..40, lon -50..50, 35 km) whose outer 5
/// degrees are thinned to 20 km (a rifted margin, the v1 weak zone); plate 1
/// = ocean. With `member2`, ocean at lon 50..90 is plate 2.
fn world(member2: bool) -> (GridSpec, CrustState) {
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
            th[c] = if lat.abs() > 35.0 || lon.abs() > 45.0 { 20_000.0 } else { 35_000.0 };
        } else if member2 && lat.abs() < 40.0 && lon >= 50.0 && lon < 90.0 {
            plate[c] = PlateId(2);
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

struct Out {
    plate: Vec<PlateId>,
    plates: Vec<PlateRecordV0>,
    group: Vec<PlateId>,
    side_b_share: f64,
    kind: Vec<CrustKindCodeV0>,
}

/// One pass at t 200 -> 202 with the supercontinent tenure already served:
/// a scripted breakup of the group {0, 1} (continent sutured to its ocean).
fn breakup(law: ReorgLaw) -> (GridSpec, Out) {
    let (grid, mut st) = world(false);
    let mut plates = vec![rec(0, 10.0, 20.0, 0.3), rec(1, 10.0, 20.0, 0.3)];
    let mut layer = ParcelLayer::seed(&grid, &st).unwrap();
    let mut sut = SutureState::default();
    sut.parent.insert(PlateId(1), PlateId(0));
    let mut rs = ReorgState { tenure: vec![(0, 0.0)], ..Default::default() };
    let p = params();
    let ctx = ReorgCtx { grid: &grid, seed_path: &sp(), params: &p, schedule: &[] };
    let rep = reorg_pass_law(&ctx, law, &mut st, &mut plates, &mut layer, &mut sut, &mut rs, 200.0, 202.0).unwrap();
    assert_eq!(rep.breakups, 1, "{law:?}: the breakup happens");
    let e = rs.events.iter().find(|e| e.kind == "breakup").unwrap();
    let group = st.plate.iter().map(|q| sut.find(*q)).collect();
    (grid, Out { plate: st.plate.clone(), plates, group, side_b_share: e.side_b_share, kind: st.kind.clone() })
}

/// Share of the new group-boundary length (4-neighbour edges) on
/// continent/ocean edges, and the continental-cell shares of the two sides.
fn coast_share(grid: &GridSpec, o: &Out) -> (f64, f64) {
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    let (mut tot, mut coast) = (0.0f64, 0.0f64);
    for i in 0..nlat {
        for j in 0..nlon {
            let c = i * nlon + j;
            let mut nb = vec![i * nlon + (j + 1) % nlon];
            if i + 1 < nlat {
                nb.push((i + 1) * nlon + j);
            }
            for k in nb {
                if o.group[c] != o.group[k] {
                    tot += 1.0;
                    if o.kind[c] != o.kind[k] {
                        coast += 1.0;
                    }
                }
            }
        }
    }
    let root_a = o.group[(nlat / 2) * nlon];
    let (mut ca, mut cn) = (0.0, 0.0);
    for c in 0..grid.cells() {
        if o.kind[c] == CrustKindCodeV0::Continental {
            cn += 1.0;
            if o.group[c] == root_a {
                ca += 1.0;
            }
        }
    }
    (coast / tot.max(1.0), ca / cn)
}

#[test]
fn v2_breakup_crosses_the_continent_not_its_margin() {
    let (grid, v2) = breakup(ReorgLaw::ScriptedV2);
    let (grid1, v1) = breakup(ReorgLaw::ScriptedV1);
    let (c2, a2) = coast_share(&grid, &v2);
    let (c1, a1) = coast_share(&grid1, &v1);
    eprintln!("coast-boundary share v1 {c1:.3} v2 {c2:.3}; continental share of one side v1 {a1:.3} v2 {a2:.3}; side_b v1 {:.3} v2 {:.3}", v1.side_b_share, v2.side_b_share);
    assert!((0.25..=0.75).contains(&v2.side_b_share), "v2 side B continental share {}", v2.side_b_share);
    assert!((0.25..=0.75).contains(&a2), "v2 sides {a2}");
    assert!(c2 < 0.15, "v2 boundary on coast edges {c2}");
    assert!(c1 > c2, "v1 reference traces more coast ({c1}) than v2 ({c2})");
}

#[test]
fn v2_breakup_is_deterministic() {
    let (_, a) = breakup(ReorgLaw::ScriptedV2);
    let (_, b) = breakup(ReorgLaw::ScriptedV2);
    assert_eq!(a.plate, b.plate);
    assert_eq!(a.plates, b.plates);
}

/// Plate 2 (ocean, member of plate 0's group) is consolidated into plate 0
/// at a scheduled event.
fn fold(law: ReorgLaw, passes: u32) -> Option<PlateRetirementV0> {
    let (grid, mut st) = world(true);
    let mut plates = vec![rec(0, 10.0, 20.0, 0.3), rec(1, -30.0, 100.0, 0.5), rec(2, 10.0, 20.0, 0.3)];
    let mut layer = ParcelLayer::seed(&grid, &st).unwrap();
    let mut sut = SutureState::default();
    sut.parent.insert(PlateId(2), PlateId(0));
    let mut rs = ReorgState::default();
    let p = params();
    let schedule = [2.0];
    let ctx = ReorgCtx { grid: &grid, seed_path: &sp(), params: &p, schedule: &schedule };
    for k in 0..passes {
        let t0 = 1.0 + 2.0 * k as f64;
        reorg_pass_law(&ctx, law, &mut st, &mut plates, &mut layer, &mut sut, &mut rs, t0, t0 + 2.0).unwrap();
    }
    assert!(st.plate.iter().all(|q| *q != PlateId(2)), "plate 2 folded");
    plates.iter().find(|r| r.plate_id == PlateId(2)).unwrap().retirement.map(|r| r.cause)
}

#[test]
fn folded_plates_retire_as_folded_under_v2_only() {
    assert_eq!(fold(ReorgLaw::ScriptedV2, 1), Some(PlateRetirementV0::Folded));
    assert_eq!(fold(ReorgLaw::ScriptedV1, 1), None);
    assert_eq!(fold(ReorgLaw::ScriptedV1, 2), Some(PlateRetirementV0::Subducted));
}
