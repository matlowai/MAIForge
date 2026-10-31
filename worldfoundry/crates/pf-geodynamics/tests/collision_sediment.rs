//! SEDIMENT-TOWERS (lane ST, 2026-10-12): `accreted-to-belt-crust-v1`
//! against the legacy `surface-layer-v0` collision-sediment rule, on lane
//! CC's head-on synthetic (two all-continental hemispheres converging along
//! a tilted boundary, 36x72 grid), with 2 km of sediment everywhere.
//!
//! - Legacy: collision-loser sediment is added to the belt's surface layer,
//!   so suture cells end far above the initial 2 km (the tower mechanism).
//! - Accreted: no collision move adds surface sediment (the carry only
//!   copies columns), so no cell exceeds 2 km; the sediment leaves the
//!   sediment field and enters the crust in one step exactly as named.
//! - `SurfaceLayerV0` through `advance_with_laws` is `advance_with_carry`
//!   bit for bit; the accreted law without the belt is refused.

use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_geodynamics::carry::CarryLaw;
use pf_geodynamics::kinematic::{advance_with_carry, advance_with_laws, CrustState};
use pf_geodynamics::orogen::{BeltParams, CollisionSedimentLaw};
use pf_geodynamics::rotation::{lat_lon_to_unit, unit_to_lat_lon};
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;

const R: f64 = 6.371e6;
const RATE: f64 = 0.5; // deg/My per plate
const SED0: f64 = 2_000.0;

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

fn setup() -> (GridSpec, Vec<PlateRecordV0>, CrustState, Vec<f64>) {
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
    let st = CrustState {
        t_my: 0.0,
        plate: cents
            .iter()
            .map(|p| PlateId(if p[0] * m[0] + p[1] * m[1] + p[2] * m[2] >= 0.0 { 0 } else { 1 }))
            .collect(),
        kind: vec![CrustKindCodeV0::Continental; n],
        birth_t_my: vec![-2000.0; n],
        thickness_m: vec![35_000.0; n],
        sediment_m: vec![SED0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: cents,
    };
    (grid, plates, st, areas)
}

fn vol(v: &[f64], areas: &[f64]) -> f64 {
    v.iter().zip(areas).map(|(h, a)| h * a).sum()
}

#[test]
fn accreted_collision_sediment_builds_no_towers_and_closes() {
    let (grid, plates, st0, areas) = setup();
    let prior = UpliftPrior::preview_v0();
    let b = belt();
    let max_after = |law: CollisionSedimentLaw| -> (f64, f64) {
        let mut st = st0.clone();
        let mut accreted = 0.0;
        for _ in 0..20 {
            let (next, l, _, _) = advance_with_laws(
                &grid,
                &st,
                &plates,
                &prior,
                1.0,
                Some(&b),
                CarryLaw::FrontRescueV1,
                law,
            )
            .unwrap();
            accreted += l.sediment_accreted_m3;
            st = next;
        }
        (st.sediment_m.iter().cloned().fold(0.0, f64::max), accreted)
    };
    let (legacy_max, legacy_acc) = max_after(CollisionSedimentLaw::SurfaceLayerV0);
    let (acc_max, acc_acc) = max_after(CollisionSedimentLaw::AccretedToBeltCrustV1);
    // Legacy: the belt's surface layer piles loser sediment on suture cells.
    assert!(legacy_max > 2.0 * SED0, "legacy max sediment {legacy_max}");
    assert_eq!(legacy_acc, 0.0);
    // Accreted: no cell gains surface sediment in a collision.
    assert!(acc_max <= SED0, "accreted max sediment {acc_max}");
    assert!(acc_acc > 0.0);

    // One step from the same state: identical claims, so the accreted volume
    // is exactly what leaves the sediment field and enters the crust
    // (delamination counted), and the belt places all of its input.
    let run1 = |law| {
        advance_with_laws(&grid, &st0, &plates, &prior, 1.0, Some(&b), CarryLaw::FrontRescueV1, law)
            .unwrap()
    };
    let (sl, ll, _, _) = run1(CollisionSedimentLaw::SurfaceLayerV0);
    let (sa, la, _, _) = run1(CollisionSedimentLaw::AccretedToBeltCrustV1);
    let acc = la.sediment_accreted_m3;
    assert!(acc > 0.0);
    assert!((ll.belt_sediment_m3 - acc).abs() <= 1e-9 * acc);
    let d_sed = vol(&sl.sediment_m, &areas) - vol(&sa.sediment_m, &areas);
    let d_crust = (vol(&sa.thickness_m, &areas) + la.delaminated_m3)
        - (vol(&sl.thickness_m, &areas) + ll.delaminated_m3);
    let scale = vol(&st0.thickness_m, &areas);
    assert!((d_sed - acc).abs() <= 1e-12 * scale, "sediment {d_sed} vs accreted {acc}");
    assert!((d_crust - acc).abs() <= 1e-12 * scale, "crust {d_crust} vs accreted {acc}");
    assert!((la.belt_in_m3 - la.belt_placed_m3).abs() <= 1e-9 * la.belt_in_m3);
    // Both residuals name the move: equal to the legacy residuals.
    assert!((la.residual_m3() - ll.residual_m3()).abs() <= 1e-12 * scale);
    assert!((la.sediment_residual_m3() - ll.sediment_residual_m3()).abs() <= 1e-12 * scale);
}

#[test]
fn surface_layer_is_the_legacy_carry_bit_for_bit_and_accretion_needs_the_belt() {
    let (grid, plates, st0, _) = setup();
    let prior = UpliftPrior::preview_v0();
    let b = belt();
    for (belt, carry) in [
        (None, CarryLaw::SemiLagrangianNearestMaterialV0),
        (Some(&b), CarryLaw::SemiLagrangianNearestMaterialV0),
        (Some(&b), CarryLaw::FrontRescueV1),
    ] {
        let old = advance_with_carry(&grid, &st0, &plates, &prior, 1.0, belt, carry).unwrap();
        let new = advance_with_laws(
            &grid,
            &st0,
            &plates,
            &prior,
            1.0,
            belt,
            carry,
            CollisionSedimentLaw::SurfaceLayerV0,
        )
        .unwrap();
        assert_eq!(old.0, new.0);
        assert_eq!(old.1, new.1);
        assert_eq!(old.2, new.2);
        assert_eq!(old.3, new.3);
    }
    assert!(advance_with_laws(
        &grid,
        &st0,
        &plates,
        &prior,
        1.0,
        None,
        CarryLaw::SemiLagrangianNearestMaterialV0,
        CollisionSedimentLaw::AccretedToBeltCrustV1,
    )
    .is_err());
}
