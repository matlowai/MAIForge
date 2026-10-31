//! WS-CONST (lane WC, 2026-10-12): the preview chain reads its physical
//! constants through `UniverseV1` and computes exactly what it computed
//! before. The pre-WS-CONST arithmetic is RESTATED here from the old
//! source (backed up in `generated-test-data/`) and compared
//! bitwise on evolved states; the full-run byte identity is proven by the
//! lane's CLI and flag-mode before/after runs (REPORT.md).

use pf_contracts::laws::{legacy, IsostasyLawV1};
use pf_contracts::params::{FidelityPreset, ScenarioPreset};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::isostasy::{
    two_column_field_with, ColumnV0, TwoColumnConstants, REFERENCE_CRUST_M,
};
use pf_geodynamics::kinematic::CrustState;
use pf_geodynamics::uplift::{
    UpliftPrior, CONTINENTAL_CRUST_THICKNESS_M, OCEANIC_CRUST_THICKNESS_M,
};
use pf_lab::preview::{
    check_profile, elevation_from_state, elevation_on_grid, gdh1_thermal_elevation_m, simulate,
    simulate_with, PreviewConfig, PreviewOptions, WATER_DENSITY_KG_M3,
};

fn mini(seed: u64) -> PreviewConfig {
    PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity: FidelityPreset::Preview,
        grid: Some((36, 72)),
        seed,
        duration_my: Some(20.0),
        step_my: Some(2.0),
    }
}

fn areas(cfg: &PreviewConfig, p: &pf_contracts::params::PlanetParameterProfile) -> Vec<f64> {
    let grid = cfg.grid(p).unwrap();
    let mut out = Vec::new();
    for i in 0..grid.nlat() {
        let a = grid.cell_area_m2(i).unwrap();
        for _ in 0..grid.nlon() {
            out.push(a);
        }
    }
    out
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn ws_const_lifted_legacy_constants_equal_the_code_constants() {
    assert_eq!(
        legacy::REFERENCE_CRUST_THICKNESS_M.to_bits(),
        REFERENCE_CRUST_M.to_bits()
    );
    assert_eq!(
        legacy::CONTINENTAL_CRUST_THICKNESS_M.to_bits(),
        CONTINENTAL_CRUST_THICKNESS_M.to_bits()
    );
    assert_eq!(
        legacy::OCEANIC_CRUST_THICKNESS_M.to_bits(),
        OCEANIC_CRUST_THICKNESS_M.to_bits()
    );
    assert_eq!(
        legacy::SEAWATER_DENSITY_KG_M3.to_bits(),
        WATER_DENSITY_KG_M3.to_bits()
    );
    // The prior the chain builds from the universe equals the default prior.
    let p = mini(1).profile().unwrap();
    let u = p.universe().unwrap();
    let prior = UpliftPrior::preview_v0();
    assert_eq!(
        u.constants.oceanic_crust_thickness_m.value,
        prior.record.oceanic_crust_thickness_m
    );
    assert_eq!(
        u.constants.continental_crust_thickness_m.value,
        prior.record.continental_crust_thickness_m
    );
}

/// The pre-WS-CONST `elevation_from_state`, restated verbatim in its
/// arithmetic (planet block + REFERENCE_CRUST_M + WATER_DENSITY_KG_M3).
fn old_elevation_from_state(
    state: &CrustState,
    params: &pf_contracts::params::PlanetParameterProfile,
) -> Vec<f64> {
    let rho_m = params.planet.mantle_density_kg_m3;
    let rho_c = params.planet.crust_density_kg_m3;
    let water_factor = rho_m / (rho_m - WATER_DENSITY_KG_M3);
    let sea = params.planet.sea_level_m;
    (0..state.cells())
        .map(|c| {
            let column = state.thickness_m[c] + state.sediment_m[c];
            let airy = (rho_m - rho_c) / rho_m * (column - REFERENCE_CRUST_M);
            let loaded = if airy < sea {
                sea + (airy - sea) * water_factor
            } else {
                airy
            };
            let thermal = if state.kind[c] == CrustKindCodeV0::Oceanic {
                gdh1_thermal_elevation_m(state.age_my(c))
            } else {
                0.0
            };
            loaded + thermal
        })
        .collect()
}

#[test]
fn ws_const_single_column_path_is_bitwise_the_old_arithmetic() {
    for seed in [240907u64, 7, 3] {
        let cfg = mini(seed);
        let out = simulate(&cfg).unwrap();
        let p = &out.params;
        assert!(p.universe().unwrap().law_set.is_legacy());
        let new = elevation_from_state(&out.state, p).unwrap();
        let old = old_elevation_from_state(&out.state, p);
        assert_eq!(bits(&new), bits(&old), "seed {seed}");
        let (grid_e, solve) = elevation_on_grid(&out.state, p, &areas(&cfg, p)).unwrap();
        assert!(solve.is_none());
        assert_eq!(bits(&grid_e), bits(&old), "seed {seed}");
        // The evolved state exercises sediment and oceanic columns.
        assert!(out.state.sediment_m.iter().any(|s| *s > 0.0));
        assert!(out
            .state
            .kind
            .iter()
            .any(|k| *k == CrustKindCodeV0::Oceanic));
    }
}

#[test]
fn ws_const_two_column_path_is_bitwise_the_old_construction() {
    for (seed, opts) in [
        (240907u64, PreviewOptions::two_column_v0()),
        (7, PreviewOptions::crust_bodies_v0()),
    ] {
        let cfg = mini(seed);
        let out = simulate_with(&cfg, &opts).unwrap();
        let p = &out.params;
        let a = areas(&cfg, p);
        let (new_e, new_s) = elevation_on_grid(&out.state, p, &a).unwrap();
        let new_s = new_s.unwrap();
        // Old construction: constants straight from the blocks, depth =
        // deepest root + margin (lane K2), sediment at the crust density.
        let b = p.two_column_isostasy.as_ref().unwrap();
        let rho_m = p.planet.mantle_density_kg_m3;
        let rho_cc = p.planet.crust_density_kg_m3;
        let cols: Vec<ColumnV0> = (0..out.state.cells())
            .map(|c| ColumnV0 {
                kind: out.state.kind[c],
                crust_m: out.state.thickness_m[c],
                sediment_m: out.state.sediment_m[c],
                thermal_uplift_water_loaded_m: if out.state.kind[c] == CrustKindCodeV0::Oceanic {
                    gdh1_thermal_elevation_m(out.state.age_my(c))
                } else {
                    0.0
                },
            })
            .collect();
        let mut deepest = f64::MIN;
        for col in &cols {
            let rho_c = match col.kind {
                CrustKindCodeV0::Continental => rho_cc,
                CrustKindCodeV0::Oceanic => b.oceanic_crust_density_kg_m3,
            };
            let t_dry =
                col.thermal_uplift_water_loaded_m * (rho_m - b.seawater_density_kg_m3) / rho_m;
            deepest = deepest.max((rho_c * col.crust_m + rho_cc * col.sediment_m) / rho_m - t_dry);
        }
        let depth = deepest.max(0.0) + b.compensation_depth_margin_m.unwrap();
        let k_old = TwoColumnConstants {
            mantle_density_kg_m3: rho_m,
            continental_crust_density_kg_m3: rho_cc,
            oceanic_crust_density_kg_m3: b.oceanic_crust_density_kg_m3,
            sediment_density_kg_m3: rho_cc,
            seawater_density_kg_m3: b.seawater_density_kg_m3,
            compensation_depth_m: depth,
            ocean_water_volume_m3: b.ocean_water_volume_m3,
            output_sea_level_datum_m: p.planet.sea_level_m,
        };
        // The universe-resolved constants are exactly the old ones.
        let u = p.universe().unwrap();
        assert_eq!(TwoColumnConstants::from_universe(&u, depth).unwrap(), k_old);
        let old = two_column_field_with(&cols, &a, &k_old).unwrap();
        assert_eq!(bits(&new_e), bits(&old.elevation_m), "seed {seed}");
        assert_eq!(
            new_s.sea_level_gauge_m.to_bits(),
            old.sea_level_gauge_m.to_bits()
        );
        assert_eq!(new_s.compensation_depth_m.to_bits(), depth.to_bits());
        // Fixed-depth mode: from_params is the old block-derived constants.
        let mut q = p.clone();
        q.two_column_isostasy
            .as_mut()
            .unwrap()
            .compensation_depth_margin_m = None;
        let k_fixed = TwoColumnConstants::from_params(&q).unwrap();
        assert_eq!(
            k_fixed,
            TwoColumnConstants {
                compensation_depth_m: b.compensation_depth_m,
                ..k_old
            }
        );
        assert!(matches!(
            q.universe().unwrap().law_set.isostasy,
            IsostasyLawV1::TwoColumnCommonDepthV0 {
                derived_depth_margin_m: None,
                ..
            }
        ));
    }
}

#[test]
fn ws_const_pinned_profile_runs_the_same_physics_and_drift_is_refused() {
    for opts in [
        PreviewOptions::legacy(),
        PreviewOptions::two_column_v0(),
        PreviewOptions::crust_bodies_v0(),
    ] {
        let cfg = mini(7);
        let out = simulate_with(&cfg, &opts).unwrap();
        let p = out.params.clone();
        let mut pinned = p.clone();
        pinned.pin_universe().unwrap();
        check_profile(&pinned).unwrap();
        let a = areas(&cfg, &p);
        let (e0, s0) = elevation_on_grid(&out.state, &p, &a).unwrap();
        let (e1, s1) = elevation_on_grid(&out.state, &pinned, &a).unwrap();
        assert_eq!(bits(&e0), bits(&e1));
        assert_eq!(s0, s1);
        // Drift between the pinned universe and the blocks is refused
        // loudly before any physics runs.
        let mut drift = pinned.clone();
        drift.planet.mantle_density_kg_m3 = 3310.0;
        let e = check_profile(&drift).unwrap_err().to_string();
        assert!(e.contains("pinned universe differs"), "{e}");
        assert!(e.contains("mantle_density_kg_m3"), "{e}");
        assert!(elevation_on_grid(&out.state, &drift, &a).is_err());
    }
}
