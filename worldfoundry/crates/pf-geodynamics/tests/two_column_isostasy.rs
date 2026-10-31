//! Known-answer tests for the OPT-IN two-column isostasy and the
//! volume-primary sea-level solve (FREEBOARD-SPIKE, lane F, 2026-10-11).
//!
//! Expected values are computed BY HAND in this file (closed forms written
//! out with their arithmetic), never by calling the producer to get them.
//! Mass balance above the compensation depth is re-derived here
//! independently from the producer's outputs.

use pf_contracts::params::{
    Compensation, FidelityPreset, MantleBlock, PlanetParameterProfile, ScenarioPreset,
    TwoColumnIsostasyBlock,
};
use pf_contracts::terrain::CrustKindCodeV0::{Continental, Oceanic};
use pf_geodynamics::isostasy::{
    airy_elevation_m, dry_elevation_m, isostatic_elevation_m, loaded_elevation_m,
    ocean_volume_m3, solve_sea_level_m, two_column_field, ColumnV0, TwoColumnConstants,
    REFERENCE_CRUST_M,
};

// Densities used by every hand calculation below (kg/m^3).
const RHO_M: f64 = 3300.0; // planet.mantle_density_kg_m3 default
const RHO_CC: f64 = 2800.0; // planet.crust_density_kg_m3 default (continental)
const RHO_OC: f64 = 2900.0; // TwoColumnIsostasyBlock::earth_v0 oceanic crust
const RHO_W: f64 = 1030.0; // sea water
const D: f64 = 100_000.0; // compensation depth (m)

fn params(volume_m3: f64) -> PlanetParameterProfile {
    let mut p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
    p.mantle.compensation = Compensation::Airy;
    let mut b = TwoColumnIsostasyBlock::earth_v0();
    b.ocean_water_volume_m3 = volume_m3;
    p.two_column_isostasy = Some(b);
    p
}

fn k(volume_m3: f64) -> TwoColumnConstants {
    TwoColumnConstants::from_params(&params(volume_m3)).unwrap()
}

fn col(kind: pf_contracts::terrain::CrustKindCodeV0, crust_m: f64) -> ColumnV0 {
    ColumnV0 {
        kind,
        crust_m,
        sediment_m: 0.0,
        thermal_uplift_water_loaded_m: 0.0,
    }
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

#[test]
fn constants_resolve_with_stated_values() {
    let k = k(1.0e18);
    assert_eq!(k.mantle_density_kg_m3, RHO_M);
    assert_eq!(k.continental_crust_density_kg_m3, RHO_CC);
    assert_eq!(k.oceanic_crust_density_kg_m3, RHO_OC);
    assert_eq!(k.seawater_density_kg_m3, RHO_W);
    assert_eq!(k.compensation_depth_m, D);
    // phi = 3300 / 2270
    assert!(close(k.water_load_factor(), 3300.0 / 2270.0, 1e-15));
}

#[test]
fn equal_columns_stand_at_equal_elevation() {
    let k = k(0.0);
    let a = dry_elevation_m(&col(Continental, 35_000.0), &k).unwrap();
    let b = dry_elevation_m(&col(Continental, 35_000.0), &k).unwrap();
    assert_eq!(a, b);
    // And in a field of identical columns with unequal areas, all share
    // one elevation whatever the water volume.
    let cols = vec![col(Oceanic, 7_000.0); 5];
    let areas = [1.0e10, 2.0e10, 3.0e10, 4.0e10, 5.0e10];
    let sol = two_column_field(&cols, &areas, &params(1.0e13)).unwrap();
    for e in &sol.elevation_m {
        assert_eq!(*e, sol.elevation_m[0]);
    }
}

#[test]
fn lighter_or_thicker_crust_stands_higher_by_the_hand_amount() {
    let k = k(0.0);
    // Same thickness 10 km, continental (2800) vs oceanic (2900):
    // difference = 10000 * (2900 - 2800) / 3300 = 303.0303... m.
    let light = dry_elevation_m(&col(Continental, 10_000.0), &k).unwrap();
    let dense = dry_elevation_m(&col(Oceanic, 10_000.0), &k).unwrap();
    assert!(light > dense);
    assert!(close(light - dense, 1_000_000.0 / 3300.0, 1e-6));
    // Thicker continental crust: 40 km vs 35 km at 2800:
    // difference = 5000 * 500 / 3300 = 757.5757... m.
    let thick = dry_elevation_m(&col(Continental, 40_000.0), &k).unwrap();
    let base = dry_elevation_m(&col(Continental, 35_000.0), &k).unwrap();
    assert!(close(thick - base, 2_500_000.0 / 3300.0, 1e-6));
    // Absolute value on the mantle free-surface gauge: 35 km at 2800
    // stands at 35000 * 500 / 3300 = 5303.0303... m.
    assert!(close(base, 17_500_000.0 / 3300.0, 1e-6));
}

#[test]
fn compensation_depth_drops_out_and_a_root_below_it_is_loud() {
    let mut p2 = params(0.0);
    p2.two_column_isostasy.as_mut().unwrap().compensation_depth_m = 200_000.0;
    let k1 = k(0.0);
    let k2 = TwoColumnConstants::from_params(&p2).unwrap();
    for c in [col(Continental, 35_000.0), col(Oceanic, 7_000.0), col(Continental, 70_000.0)] {
        let a = dry_elevation_m(&c, &k1).unwrap();
        let b = dry_elevation_m(&c, &k2).unwrap();
        assert!(close(a, b, 1e-6), "D changed the elevation: {a} vs {b}");
    }
    let mut p3 = params(0.0);
    p3.two_column_isostasy.as_mut().unwrap().compensation_depth_m = 30_000.0;
    let k3 = TwoColumnConstants::from_params(&p3).unwrap();
    // 35 km of crust: root at 5303 - 35000 = -29697 m, above -30000: ok.
    assert!(dry_elevation_m(&col(Continental, 35_000.0), &k3).is_ok());
    // 40 km: root at 6060.6 - 40000 = -33939 m, below -30000: refused.
    assert!(dry_elevation_m(&col(Continental, 40_000.0), &k3).is_err());
}

#[test]
fn zero_water_puts_sea_level_at_the_deepest_point() {
    let k = k(0.0);
    let e = [10.0, -250.0, 3.0, 400.0];
    let a = [1.0, 2.0, 3.0, 4.0];
    let s = solve_sea_level_m(&e, &a, 0.0, &k).unwrap();
    assert_eq!(s, -250.0);
    // Field form: the deepest column sits exactly at the datum, every
    // other column is dry above it by its dry-elevation difference.
    let cols = [col(Oceanic, 7_000.0), col(Continental, 35_000.0)];
    let sol = two_column_field(&cols, &[1.0e12, 1.0e12], &params(0.0)).unwrap();
    assert_eq!(sol.elevation_m[0], 0.0);
    // 35000*500/3300 - 7000*400/3300 = 14_700_000/3300 = 4454.5454... m
    assert!(close(sol.elevation_m[1], 14_700_000.0 / 3300.0, 1e-6));
    assert_eq!(sol.flooded_area_fraction, 0.0);
}

#[test]
fn ocean_volume_is_monotone_and_the_solve_inverts_it() {
    let k = k(0.0);
    let e: Vec<f64> = (0..50).map(|i| ((i * 37) % 50) as f64 * 113.0 - 3000.0).collect();
    let a: Vec<f64> = (0..50).map(|i| 1.0e9 * (1.0 + (i % 7) as f64)).collect();
    let mut last = -1.0;
    let mut s = -4000.0;
    while s < 4000.0 {
        let v = ocean_volume_m3(&e, &a, s, &k);
        assert!(v >= last, "V(S) decreased at S = {s}");
        if s > -3000.0 {
            assert!(v > last, "V(S) not strictly increasing above the deepest point at {s}");
        }
        last = v;
        s += 17.0;
    }
    for target in [1.0e9, 3.3e12, 7.7e13, 1.0e15, 5.0e15] {
        let s = solve_sea_level_m(&e, &a, target, &k).unwrap();
        let v = ocean_volume_m3(&e, &a, s, &k);
        assert!(((v - target) / target).abs() < 1e-12, "{v} vs {target}");
    }
}

#[test]
fn analytic_two_column_world_by_hand() {
    // Two cells of equal area A: one continent (35 km at 2800), one ocean
    // (7 km at 2900). Water volume chosen so the ocean is 5000 m deep:
    // V = A * 5000.
    //   dry continent  e_c = 35000 * 500 / 3300 = 17_500_000 / 3300
    //   dry ocean      e_o =  7000 * 400 / 3300 =  2_800_000 / 3300
    //   water depth    d = phi (S - e_o) = 5000  =>  S - e_o = 5000 * 2270 / 3300
    //   freeboard      f = e_c - S
    //                    = (17_500_000 - 2_800_000 - 11_350_000) / 3300
    //                    = 3_350_000 / 3300 = 1015.1515... m
    let a = 2.0e12;
    let cols = [col(Continental, 35_000.0), col(Oceanic, 7_000.0)];
    let sol = two_column_field(&cols, &[a, a], &params(a * 5000.0)).unwrap();
    assert!(close(sol.elevation_m[0], 3_350_000.0 / 3300.0, 1e-6), "{:?}", sol.elevation_m);
    assert!(close(sol.elevation_m[1], -5000.0, 1e-6), "{:?}", sol.elevation_m);
    // S on the gauge: e_o + 5000*2270/3300 = (2_800_000 + 11_350_000)/3300
    assert!(close(sol.sea_level_gauge_m, 14_150_000.0 / 3300.0, 1e-6));
    assert!(close(sol.flooded_area_fraction, 0.5, 1e-15));
    assert!(((sol.volume_check_m3 - a * 5000.0) / (a * 5000.0)).abs() < 1e-12);

    // Independent mass balance above -D from the producer's surfaces.
    let s = sol.sea_level_gauge_m;
    let e_c = sol.elevation_m[0] + s;
    let e_o = sol.elevation_m[1] + s;
    let m_cont = RHO_CC * 35_000.0 + RHO_M * (D + e_c - 35_000.0);
    let m_ocean = RHO_W * (s - e_o) + RHO_OC * 7_000.0 + RHO_M * (D + e_o - 7_000.0);
    let m_ref = RHO_M * D;
    assert!(((m_cont - m_ref) / m_ref).abs() < 1e-12, "{m_cont} vs {m_ref}");
    assert!(((m_ocean - m_ref) / m_ref).abs() < 1e-12, "{m_ocean} vs {m_ref}");
}

#[test]
fn thermal_uplift_is_water_loaded_gdh1_relief_when_submerged() {
    // A submerged ocean column 1000 m of GDH1 (water-loaded) uplift is
    // exactly 1000 m shallower at fixed S than the same column without it.
    let k = k(0.0);
    let cold = col(Oceanic, 7_000.0);
    let mut hot = cold;
    hot.thermal_uplift_water_loaded_m = 1000.0;
    let s = 5000.0; // well above both dry surfaces (848.48 m and below)
    let e_cold = loaded_elevation_m(dry_elevation_m(&cold, &k).unwrap(), s, &k);
    let e_hot = loaded_elevation_m(dry_elevation_m(&hot, &k).unwrap(), s, &k);
    assert!(close(e_hot - e_cold, 1000.0, 1e-9), "{}", e_hot - e_cold);
    // Dry, the same thermal buoyancy is 1000 * 2270/3300 = 687.8787... m.
    let d_hot = dry_elevation_m(&hot, &k).unwrap();
    let d_cold = dry_elevation_m(&cold, &k).unwrap();
    assert!(close(d_hot - d_cold, 2_270_000.0 / 3300.0, 1e-9));
}

#[test]
fn sediment_loads_at_crust_density() {
    // 1000 m of sediment at 2800 on oceanic crust raises the dry surface by
    // 1000 * 500 / 3300 = 151.5151... m.
    let k = k(0.0);
    let bare = col(Oceanic, 7_000.0);
    let mut sed = bare;
    sed.sediment_m = 1000.0;
    let d = dry_elevation_m(&sed, &k).unwrap() - dry_elevation_m(&bare, &k).unwrap();
    assert!(close(d, 500_000.0 / 3300.0, 1e-9));
}

#[test]
fn flag_gates_both_paths_loudly() {
    // Flag on: the legacy single-column functions refuse.
    let on = params(1.0e18);
    assert!(airy_elevation_m(REFERENCE_CRUST_M, &on).is_err());
    assert!(isostatic_elevation_m(REFERENCE_CRUST_M, &on, None).is_err());
    // Flag off: the two-column solve refuses, legacy is unchanged (0 m at
    // the reference column).
    let mut off = on.clone();
    off.two_column_isostasy = None;
    assert!(TwoColumnConstants::from_params(&off).is_err());
    assert!(two_column_field(&[col(Oceanic, 7_000.0)], &[1.0], &off).is_err());
    assert!(airy_elevation_m(REFERENCE_CRUST_M, &off).unwrap().abs() < 1e-9);
    // Flexure with the flag: refused (two-column is local compensation).
    let mut flex = on.clone();
    flex.mantle.compensation = Compensation::Flexure;
    assert!(TwoColumnConstants::from_params(&flex).is_err());
    // Bad inputs are loud.
    let k = k(0.0);
    assert!(dry_elevation_m(&col(Oceanic, 0.0), &k).is_err());
    assert!(dry_elevation_m(&col(Oceanic, f64::NAN), &k).is_err());
    assert!(solve_sea_level_m(&[0.0], &[1.0], -1.0, &k).is_err());
    assert!(solve_sea_level_m(&[0.0, 1.0], &[1.0], 1.0, &k).is_err());
}

#[test]
fn excess_water_floods_everything() {
    // Two cells, areas 1 and 3 m^2, dry at 0 and 100 m. Capacity up to
    // S = 100 m: phi * 1 * 100. A further 400 phi m^3 spreads over 4 m^2,
    // so S = 100 + 100 = 200 m.
    let k = k(0.0);
    let phi = 3300.0 / 2270.0;
    let s = solve_sea_level_m(&[0.0, 100.0], &[1.0, 3.0], phi * 500.0, &k).unwrap();
    assert!(close(s, 200.0, 1e-9), "{s}");
}
