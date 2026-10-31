//! SEED-CARRY-CONSERVE (lane CV, 2026-10-13): synthetic checks of
//! `continental-parcel-carry-v1` and `contact-length-pole-share-v1`
//! (pre-registered in `generated-test-data/PREREG.txt`).
//!
//! 1. A continental disc drifting with its plate over the other plate's
//!    ocean keeps its area (parcels exactly, the grid projection within a
//!    few percent) at dt 1, 2 and 4 My.
//! 2. Two continental hemispheres converging head-on: the consumed area per
//!    My is near the exact closure 4 r R^2 and the belt receives exactly the
//!    consumed volume.
//! 3. The write-back hands exactly sum(cell area x change) to the parcels.
//! 4. A long converging continental contact sutures the two plates (the
//!    smaller takes the larger's pole; equal areas: the lower id survives).

use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_geodynamics::conserve::{suture_pass, write_back, write_back_corrected, ParcelLayer, ParcelParams, SutureState};
use pf_geodynamics::kinematic::{advance_with_parcels, CrustState, ExtrasLaws};
use pf_geodynamics::orogen::{BeltParams, CollisionSedimentLaw};
use pf_geodynamics::rotation::{lat_lon_to_unit, unit_to_lat_lon};
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;

const R: f64 = 6.371e6;
const RATE: f64 = 0.5; // deg/My

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

fn params() -> ParcelParams {
    ParcelParams {
        overlap_cell_fraction: 0.5,
        fragment_min_area_fraction: 1.0e-3,
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn rec(id: u32, axis: [f64; 3], r: f64) -> PlateRecordV0 {
    let (lat, lon) = unit_to_lat_lon(axis).unwrap();
    PlateRecordV0 {
        plate_id: PlateId(id),
        parent_plate: None,
        birth_t_my: 0.0,
        euler_pole: Some(EulerPoleV0 {
            lat_deg: lat,
            lon_deg: lon,
            rate_deg_per_my: r,
        }),
        retirement: None,
    }
}

struct World {
    grid: GridSpec,
    areas: Vec<f64>,
    cents: Vec<[f64; 3]>,
}

fn world(nlat: usize, nlon: usize) -> World {
    let grid = GridSpec::new(nlat, nlon, R).unwrap();
    let n = grid.cells();
    let areas = (0..n).map(|c| grid.cell_area_m2(c / nlon).unwrap()).collect();
    let cents = (0..n)
        .map(|c| {
            lat_lon_to_unit(
                grid.center_lat_deg(c / nlon).unwrap(),
                grid.center_lon_deg(c % nlon).unwrap(),
            )
            .unwrap()
        })
        .collect();
    World { grid, areas, cents }
}

fn state(w: &World, plate: Vec<PlateId>, cont: Vec<bool>) -> CrustState {
    let n = w.grid.cells();
    CrustState {
        t_my: 0.0,
        plate,
        kind: cont
            .iter()
            .map(|&c| if c { CrustKindCodeV0::Continental } else { CrustKindCodeV0::Oceanic })
            .collect(),
        birth_t_my: vec![-100.0; n],
        thickness_m: cont.iter().map(|&c| if c { 35_000.0 } else { 7_000.0 }).collect(),
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: w.cents.clone(),
    }
}

fn cont_area(w: &World, st: &CrustState) -> f64 {
    (0..w.grid.cells())
        .filter(|&c| st.kind[c] == CrustKindCodeV0::Continental)
        .map(|c| w.areas[c])
        .sum()
}

/// Run `steps` of the parcel carry with the write-back after each step.
/// Returns (final state, layer, consumed area, belt_in, consumed volume,
/// max |parcel volume residual| / volume).
fn run(
    w: &World,
    mut st: CrustState,
    plates: &[PlateRecordV0],
    dt: f64,
    duration: f64,
) -> (CrustState, ParcelLayer, f64, f64, f64, f64) {
    let prior = UpliftPrior::preview_v0();
    let b = belt();
    let mut layer = ParcelLayer::seed(&w.grid, &st).unwrap();
    let steps = (duration / dt).round() as usize;
    let (mut consumed_a, mut belt_in, mut consumed_v, mut worst) = (0.0, 0.0, 0.0, 0.0f64);
    for _ in 0..steps {
        let (next, mut l, _, _) = advance_with_parcels(
            &w.grid,
            &st,
            plates,
            &prior,
            dt,
            Some(&b),
            CollisionSedimentLaw::SurfaceLayerV0,
            &ExtrasLaws::default(),
            &mut layer,
            params(),
        )
        .unwrap();
        let pr = l.parcels.clone().unwrap();
        let proj = l.parcel_projection.take().unwrap();
        let wb = write_back(&w.grid, &w.areas, &mut layer, &proj, &next).unwrap();
        let resid = wb.volume_after_m3
            - (pr.volume_before_m3 - pr.consumed_m3 + wb.writeback_m3 + wb.created_m3 + wb.clamped_m3);
        worst = worst.max(resid.abs() / pr.volume_before_m3.max(1.0));
        consumed_a += pr.consumed_area_m2;
        consumed_v += pr.consumed_m3;
        belt_in += l.belt_in_m3;
        st = next;
    }
    (st, layer, consumed_a, belt_in, consumed_v, worst)
}

#[test]
fn drifting_disc_keeps_its_area_at_any_dt() {
    let w = world(36, 72);
    let n = w.grid.cells();
    // Plate 0: the hemisphere around (0, 0) with a 30-degree continental
    // disc at its centre; plate 1: the other hemisphere, all ocean, fixed.
    let c0 = lat_lon_to_unit(0.0, 0.0).unwrap();
    let plate: Vec<PlateId> = (0..n)
        .map(|c| PlateId(if dot(w.cents[c], c0) >= 0.0 { 0 } else { 1 }))
        .collect();
    let cont: Vec<bool> = (0..n)
        .map(|c| dot(w.cents[c], c0) >= 30f64.to_radians().cos())
        .collect();
    let axis = lat_lon_to_unit(40.0, 120.0).unwrap();
    let plates = vec![rec(0, axis, RATE), rec(1, axis, 0.0)];
    let a0 = {
        let st = state(&w, plate.clone(), cont.clone());
        cont_area(&w, &st)
    };
    let mut finals = Vec::new();
    for dt in [1.0, 2.0, 4.0] {
        let st = state(&w, plate.clone(), cont.clone());
        let (fin, layer, consumed, _, _, worst) = run(&w, st, &plates, dt, 40.0);
        let (pa, _, _) = layer.totals();
        assert_eq!(consumed, 0.0, "no other continent, nothing consumed");
        assert!((pa - a0).abs() <= 1e-9 * a0, "parcel area {pa} vs {a0}");
        assert!(worst < 1e-9, "parcel volume residual {worst}");
        let ga = cont_area(&w, &fin);
        assert!((ga - a0).abs() <= 0.05 * a0, "dt {dt}: grid area {ga} vs {a0}");
        finals.push(ga);
    }
    let lo = finals.iter().cloned().fold(f64::MAX, f64::min);
    let hi = finals.iter().cloned().fold(f64::MIN, f64::max);
    assert!(hi - lo <= 0.03 * a0, "dt spread {finals:?} of {a0}");
}

#[test]
fn head_on_collision_consumes_the_closure_and_feeds_the_belt() {
    let w = world(36, 72);
    let n = w.grid.cells();
    // Two continental hemispheres split by the plane x = 0, rotating about
    // the polar axis in opposite senses: head-on along one half of the
    // boundary (lane CC's `headon_polar`).
    let m = [1.0, 0.0, 0.0];
    let plate: Vec<PlateId> = (0..n)
        .map(|c| PlateId(if dot(w.cents[c], m) >= 0.0 { 0 } else { 1 }))
        .collect();
    let axis = [0.0, 0.0, 1.0];
    let plates = vec![rec(0, axis, RATE), rec(1, axis, -RATE)];
    let st = state(&w, plate, vec![true; n]);
    let duration = 20.0;
    let (_, _, consumed, belt_in, consumed_v, worst) = run(&w, st, &plates, 1.0, duration);
    let exact = 4.0 * RATE.to_radians() * R * R;
    let per_my = consumed / duration;
    assert!(worst < 1e-9, "parcel volume residual {worst}");
    assert!(
        per_my > 0.6 * exact && per_my < 1.4 * exact,
        "consumed {per_my:.4e} m2/My vs exact closure {exact:.4e}"
    );
    assert!((belt_in - consumed_v).abs() <= 1e-9 * consumed_v, "belt {belt_in} vs {consumed_v}");
}

#[test]
fn write_back_hands_exactly_the_cell_changes_to_the_parcels() {
    let w = world(36, 72);
    let n = w.grid.cells();
    let c0 = lat_lon_to_unit(20.0, 50.0).unwrap();
    let cont: Vec<bool> = (0..n)
        .map(|c| dot(w.cents[c], c0) >= 40f64.to_radians().cos())
        .collect();
    let st = state(&w, vec![PlateId(0); n], cont);
    let plates = vec![rec(0, lat_lon_to_unit(-10.0, 0.0).unwrap(), RATE)];
    let mut layer = ParcelLayer::seed(&w.grid, &st).unwrap();
    let (mut next, mut l, _, _) = advance_with_parcels(
        &w.grid,
        &st,
        &plates,
        &UpliftPrior::preview_v0(),
        2.0,
        Some(&belt()),
        CollisionSedimentLaw::SurfaceLayerV0,
        &ExtrasLaws::default(),
        &mut layer,
        params(),
    )
    .unwrap();
    let proj = l.parcel_projection.take().unwrap();
    let (_, v0, s0) = layer.totals();
    // An "erosion" step: every continental cell loses 100 m of crust and
    // gains 30 m of sediment.
    let (mut dv, mut ds) = (0.0, 0.0);
    for c in 0..n {
        if next.kind[c] == CrustKindCodeV0::Continental {
            next.thickness_m[c] -= 100.0;
            next.sediment_m[c] += 30.0;
            dv -= 100.0 * w.areas[c];
            ds += 30.0 * w.areas[c];
        }
    }
    let wb = write_back(&w.grid, &w.areas, &mut layer, &proj, &next).unwrap();
    assert_eq!(wb.created_parcels, 0);
    assert!((wb.writeback_m3 - dv).abs() <= 1e-9 * dv.abs());
    assert!((wb.volume_after_m3 - (v0 + dv)).abs() <= 1e-9 * v0, "{} vs {}", wb.volume_after_m3, v0 + dv);
    assert!((wb.sediment_after_m3 - (s0 + ds)).abs() <= 1e-9 * ds.abs());
}

/// HYPSO-LEVERS (2026-10-21): the write-back residual correction. A narrow
/// belt (+3 km on a 2-cell-wide band) written back plainly comes back
/// smoothed when the parcels are re-projected where they are; 3 correction
/// passes bring the re-projection much closer to the end-of-step cells and
/// leave the parcel volume exactly as the plain write-back leaves it.
#[test]
fn write_back_residual_correction_reprojects_the_end_of_step_cells() {
    let w = world(36, 72);
    let n = w.grid.cells();
    let c0 = lat_lon_to_unit(20.0, 50.0).unwrap();
    let cont: Vec<bool> = (0..n)
        .map(|c| dot(w.cents[c], c0) >= 40f64.to_radians().cos())
        .collect();
    let st = state(&w, vec![PlateId(0); n], cont);
    let plates = vec![rec(0, lat_lon_to_unit(-10.0, 0.0).unwrap(), RATE)];
    let mut layer = ParcelLayer::seed(&w.grid, &st).unwrap();
    let adv = |s: &CrustState, layer: &mut ParcelLayer, plates: &[PlateRecordV0]| {
        advance_with_parcels(&w.grid, s, plates, &UpliftPrior::preview_v0(), 2.0, Some(&belt()), CollisionSedimentLaw::SurfaceLayerV0, &ExtrasLaws::default(), layer, params())
            .unwrap()
    };
    let (mut next, mut l, _, _) = adv(&st, &mut layer, &plates);
    let proj = l.parcel_projection.take().unwrap();
    let mut band = vec![false; n];
    for c in 0..n {
        let lon = unit_to_lat_lon(w.cents[c]).unwrap().1;
        if next.kind[c] == CrustKindCodeV0::Continental && proj.continental[c] && (lon - 50.0).abs() < 5.0 {
            next.thickness_m[c] += 3000.0;
            band[c] = true;
        }
    }
    assert!(band.iter().filter(|b| **b).count() >= 6);
    let mut plain = layer.clone();
    let mut fixed = layer.clone();
    let wp = write_back(&w.grid, &w.areas, &mut plain, &proj, &next).unwrap();
    let wf = write_back_corrected(&w.grid, &w.areas, &mut fixed, &proj, &next, false, 3).unwrap();
    assert!((wf.volume_after_m3 - wp.volume_after_m3).abs() <= 1e-9 * wp.volume_after_m3, "{} vs {}", wf.volume_after_m3, wp.volume_after_m3);
    assert_eq!(wf.clamped_m3, 0.0);
    // Re-project where the parcels are (a zero-rate pole: no motion).
    let still = vec![rec(0, lat_lon_to_unit(-10.0, 0.0).unwrap(), 0.0)];
    let resid = |lay: &mut ParcelLayer| {
        let (_, mut l2, _, _) = adv(&next, lay, &still);
        let p2 = l2.parcel_projection.take().unwrap();
        let (mut s2, mut a) = (0.0, 0.0);
        for c in 0..n {
            if band[c] && p2.continental[c] {
                s2 += w.areas[c] * (p2.thickness_m[c] - next.thickness_m[c]).powi(2);
                a += w.areas[c];
            }
        }
        (s2 / a).sqrt()
    };
    let (rp, rf) = (resid(&mut plain), resid(&mut fixed));
    eprintln!("band rms re-projection residual: plain {rp:.1} m, 3 passes {rf:.1} m");
    assert!(rp > 300.0, "the plain write-back must smooth the band: rms {rp} m");
    assert!(rf < 0.5 * rp, "3 passes must at least halve the band residual: {rf} vs {rp} m");
}

#[test]
fn long_converging_contact_sutures_the_plates() {
    let w = world(36, 72);
    let n = w.grid.cells();
    let m = [1.0, 0.0, 0.0];
    let plate: Vec<PlateId> = (0..n)
        .map(|c| PlateId(if dot(w.cents[c], m) >= 0.0 { 0 } else { 1 }))
        .collect();
    let axis = [0.0, 0.0, 1.0];
    let mut plates = vec![rec(0, axis, RATE), rec(1, axis, -RATE)];
    let st = state(&w, plate, vec![true; n]);
    let mut sut = SutureState::default();
    // Threshold far above the contact: nothing happens.
    let (ev, max_contact) = suture_pass(&w.grid, &st, &mut plates, &mut sut, 1.0e8).unwrap();
    assert!(ev.is_empty() && max_contact > 2.5e6, "{max_contact}");
    assert_ne!(plates[0].euler_pole, plates[1].euler_pole);
    let (ev, _) = suture_pass(&w.grid, &st, &mut plates, &mut sut, 2.5e6).unwrap();
    assert_eq!(ev.len(), 1);
    assert_eq!(ev[0].winner_root, PlateId(0));
    assert_eq!(plates[0].euler_pole, plates[1].euler_pole);
    // Sutured: no converging contact remains.
    let (ev, max_contact) = suture_pass(&w.grid, &st, &mut plates, &mut sut, 2.5e6).unwrap();
    assert!(ev.is_empty() && max_contact == 0.0, "{max_contact}");
}
