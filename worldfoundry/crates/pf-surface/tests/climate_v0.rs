//! CLIMATE-V0 (lane CL, 2026-10-12): tests of the declared coarse climate
//! tier. Every map here is a SYNTHETIC FIXTURE (ADR 026 §1: legal as
//! instrument calibration in unit tests, never world input).

use pf_contracts::laws::{
    ClimateConstantsV1, ClimateLawV1, ClimateRecipeV1, UnitV1, UniverseError,
    CLIMATE_CONSTANT_SPECS_V1,
};
use pf_grid::GridSpec;
use pf_surface::climate::*;

const R: f64 = 6.371e6;

fn recipe() -> ClimateRecipeV1 {
    ClimateRecipeV1::coarse_v0_earth()
}

fn grid() -> GridSpec {
    GridSpec::new(36, 72, R).unwrap()
}

/// Two synthetic continents with a north-south mountain belt, so rows see
/// windward and leeward relief and both wind directions.
fn world(g: &GridSpec) -> (Vec<f64>, Vec<bool>) {
    let (nlat, nlon) = (g.nlat(), g.nlon());
    let mut z = vec![-4000.0; g.cells()];
    let mut ocean = vec![true; g.cells()];
    for i in 0..nlat {
        for j in 0..nlon {
            let c = i * nlon + j;
            let lat = g.center_lat_deg(i).unwrap();
            let a = (20..40).contains(&j) && lat.abs() < 70.0;
            let b = (50..62).contains(&j) && lat > -40.0 && lat < 50.0;
            if a || b {
                ocean[c] = false;
                z[c] = 200.0 + if j == 25 || j == 26 { 3000.0 } else { 0.0 };
            }
        }
    }
    (z, ocean)
}

fn run(g: &GridSpec, z: &[f64], ocean: &[bool], r: &ClimateRecipeV1) -> ClimateOutputV0 {
    derive_climate_coarse_v0(
        &ClimateInputs { grid: g, elevation_m: z, fill_level_m: z, sea_level_m: 0.0, ocean },
        r,
    )
    .unwrap()
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn climate_v0_determinism_identical_inputs_identical_bits() {
    let g = grid();
    let (z, o) = world(&g);
    let a = run(&g, &z, &o, &recipe());
    let b = run(&g, &z.clone(), &o.clone(), &recipe());
    for ((na, va), (nb, vb)) in a.fields.named().iter().zip(b.fields.named().iter()) {
        assert_eq!(na, nb);
        assert_eq!(bits(va), bits(vb), "{na} differs between identical runs");
    }
    assert_eq!(a, b);
}

#[test]
fn climate_v0_no_noise_zonally_symmetric_inputs_give_zonally_symmetric_outputs() {
    // Inputs that depend on latitude only (land poleward of 55 deg at 500 m,
    // ocean elsewhere): every output must be identical along each row, bit
    // for bit. A hidden longitude-dependent field (noise, hand placement)
    // would break this.
    let g = grid();
    let (nlat, nlon) = (g.nlat(), g.nlon());
    let mut z = vec![-3000.0; g.cells()];
    let mut o = vec![true; g.cells()];
    for i in 0..nlat {
        if g.center_lat_deg(i).unwrap().abs() > 55.0 {
            for j in 0..nlon {
                z[i * nlon + j] = 500.0;
                o[i * nlon + j] = false;
            }
        }
    }
    let out = run(&g, &z, &o, &recipe());
    for (name, v) in out.fields.named() {
        for i in 0..nlat {
            let first = v[i * nlon].to_bits();
            assert!(
                (0..nlon).all(|j| v[i * nlon + j].to_bits() == first),
                "{name} varies along row {i} under zonally symmetric inputs"
            );
        }
    }
    // All-land rows have no moisture source: zero precipitation, sentinel
    // distance, and they are counted.
    assert!(out.summary.rows_without_ocean > 0);
    for i in 0..nlat {
        if out.rows[i].ocean_cells == 0 {
            assert_eq!(out.fields.precipitation_mm_per_yr[i * nlon], 0.0);
            assert_eq!(
                out.fields.upwind_coast_distance_m[i * nlon],
                DISTANCE_SENTINEL_NO_UPWIND_OCEAN
            );
        }
    }
}

#[test]
fn climate_v0_no_noise_changes_only_through_inputs_and_stay_in_their_ring() {
    // Raise one land cell by 2 km: only its latitude ring may change (no
    // meridional coupling of moisture; elevation does not feed the zonal
    // energy balance). Every other row stays bit-identical.
    let g = grid();
    let (z, o) = world(&g);
    let base = run(&g, &z, &o, &recipe());
    let nlon = g.nlon();
    let i0 = 22usize;
    let c0 = i0 * nlon + 30;
    assert!(!o[c0]);
    let mut z2 = z.clone();
    z2[c0] += 2000.0;
    let pert = run(&g, &z2, &o, &recipe());
    let mut changed_in_ring = false;
    for (name, v) in base.fields.named() {
        let w = pert.fields.named().into_iter().find(|(n, _)| *n == name).unwrap().1;
        for c in 0..g.cells() {
            if c / nlon != i0 {
                assert_eq!(v[c].to_bits(), w[c].to_bits(), "{name} changed outside the ring at {c}");
            } else if v[c].to_bits() != w[c].to_bits() {
                changed_in_ring = true;
            }
        }
    }
    assert!(changed_in_ring, "the perturbation had no effect at all");
    assert!(pert.fields.temperature_mean_c[c0] < base.fields.temperature_mean_c[c0]);
    // The ocean mask is an input too: flipping one marine cell to land
    // changes outputs, again only in its ring.
    let mut o3 = o.clone();
    let c1 = 10 * nlon + 5;
    assert!(o3[c1]);
    o3[c1] = false;
    let mut z3 = z.clone();
    z3[c1] = 10.0;
    let flip = run(&g, &z3, &o3, &recipe());
    for c in 0..g.cells() {
        if c / nlon != 10 {
            assert_eq!(
                base.fields.precipitation_mm_per_yr[c].to_bits(),
                flip.fields.precipitation_mm_per_yr[c].to_bits()
            );
        }
    }
    assert_ne!(
        bits(&base.fields.precipitation_mm_per_yr),
        bits(&flip.fields.precipitation_mm_per_yr)
    );
}

#[test]
fn climate_v0_moisture_budget_closes_globally_and_per_ring() {
    let g = grid();
    let (z, o) = world(&g);
    let out = run(&g, &z, &o, &recipe());
    let nlon = g.nlon();
    for cl in &out.closures {
        assert!(cl.closes(), "{cl:?}");
    }
    let water = out.closures.iter().find(|c| c.quantity == "atmospheric-water").unwrap();
    assert!(water.gross > 0.0);
    // Per ring: E and P balance within rounding.
    for i in 0..g.nlat() {
        let a = g.cell_area_m2(i).unwrap();
        let (mut e, mut p) = (0.0, 0.0);
        for j in 0..nlon {
            e += out.fields.evaporation_mm_per_yr[i * nlon + j] * a;
            p += out.fields.precipitation_mm_per_yr[i * nlon + j] * a;
        }
        assert!((e - p).abs() <= 1e-9 * p.max(1.0), "row {i}: E {e} P {p}");
        assert!(out.rows[i].ring_residual_kg_m2.abs() <= 1e-9, "row {i}");
    }
    // Orographic rain is part of total rain, never more.
    for c in 0..g.cells() {
        assert!(out.fields.precipitation_orographic_mm_per_yr[c] <= out.fields.precipitation_mm_per_yr[c]);
        assert!(out.fields.precipitation_mm_per_yr[c] >= 0.0);
    }
}

#[test]
fn climate_v0_rain_shadow_and_windward_enhancement_known_answer() {
    // One westerly row (about 45 deg): ocean for columns 0..30, land after.
    // With a 3 km ridge at columns 40-41, the windward side gains
    // orographic rain and the lee gets less than the same cell without the
    // ridge.
    let g = grid();
    let nlon = g.nlon();
    let i = (0..g.nlat()).find(|&i| (g.center_lat_deg(i).unwrap() - 42.5).abs() < 1e-9).unwrap();
    assert!(zonal_wind_m_s(g.center_lat_deg(i).unwrap(), 6.0, 7.0, 3.0) > 0.0, "westerly row");
    let mut z = vec![-4000.0; g.cells()];
    let mut o = vec![true; g.cells()];
    for j in 30..nlon {
        z[i * nlon + j] = 100.0;
        o[i * nlon + j] = false;
    }
    let flat = run(&g, &z, &o, &recipe());
    let mut zr = z.clone();
    zr[i * nlon + 40] = 3000.0;
    zr[i * nlon + 41] = 3000.0;
    let ridge = run(&g, &zr, &o, &recipe());
    let p = |o: &ClimateOutputV0, j: usize| o.fields.precipitation_mm_per_yr[i * nlon + j];
    assert!(p(&ridge, 40) > p(&flat, 40) * 1.5, "windward {} vs {}", p(&ridge, 40), p(&flat, 40));
    assert!(ridge.fields.precipitation_orographic_mm_per_yr[i * nlon + 40] > 0.0);
    for j in 42..50 {
        assert!(p(&ridge, j) < p(&flat, j) * 0.8, "lee col {j}: {} vs {}", p(&ridge, j), p(&flat, j));
    }
    // Interior drying with distance on the flat case (continental interior).
    assert!(p(&flat, 31) > p(&flat, 50) && p(&flat, 50) > p(&flat, 70));
    // Upwind coast distance grows along the wind.
    let d = |j: usize| flat.fields.upwind_coast_distance_m[i * nlon + j];
    assert_eq!(d(29), 0.0);
    assert!(d(30) > 0.0 && d(31) > d(30));
}

#[test]
fn climate_v0_insolation_and_energy_known_answers() {
    let r = recipe();
    let s0 = r.constants.solar_constant_w_m2();
    // Derived S0 from L and d is ~1361 W/m2 for the recalled IAU values.
    assert!((s0 - 1361.0).abs() < 2.0, "{s0}");
    // Annual mean at the pole for a circular orbit: S0 sin(eps) / pi.
    let eps = r.constants.axial_tilt_deg.value;
    let (q_pole, _) = annual_insolation_w_m2(s0, 89.999, eps);
    let exact = s0 * eps.to_radians().sin() / std::f64::consts::PI;
    assert!(((q_pole - exact) / exact).abs() < 2e-3, "{q_pole} vs {exact}");
    // Global area-weighted mean of the annual mean is S0 / 4.
    let g = GridSpec::new(180, 8, R).unwrap();
    let mut s = 0.0;
    for i in 0..g.nlat() {
        let w = g.cell_area_m2(i).unwrap() * g.nlon() as f64 / g.total_area_m2();
        s += w * annual_insolation_w_m2(s0, g.center_lat_deg(i).unwrap(), eps).0;
    }
    assert!(((s - s0 / 4.0) / (s0 / 4.0)).abs() < 1e-3, "{s} vs {}", s0 / 4.0);
    // Zero tilt: no seasonal harmonic.
    assert!(annual_insolation_w_m2(s0, 45.0, 0.0).1 < 1e-9);
    // Winds: bands and calms.
    assert!(zonal_wind_m_s(15.0, 6.0, 7.0, 3.0) < 0.0);
    assert!(zonal_wind_m_s(-45.0, 6.0, 7.0, 3.0) > 0.0);
    assert!(zonal_wind_m_s(75.0, 6.0, 7.0, 3.0) < 0.0);
    // Thornthwaite: frozen gives zero, warmer gives more.
    assert_eq!(thornthwaite_pet_mm_per_yr(-20.0, 5.0), 0.0);
    assert!(thornthwaite_pet_mm_per_yr(25.0, 2.0) > thornthwaite_pet_mm_per_yr(10.0, 2.0));
}

#[test]
fn climate_v0_aquaplanet_sanity_against_recalled_earth_numbers() {
    // SANITY, not qualification: recalled Earth numbers (not re-checked):
    // equatorial annual mean ~26-27 degC, global mean surface ~14-15 degC,
    // global mean precipitation ~1000 mm/yr.
    let g = GridSpec::new(72, 8, R).unwrap();
    let z = vec![-4000.0; g.cells()];
    let o = vec![true; g.cells()];
    let out = run(&g, &z, &o, &recipe());
    let eq = out.rows[g.nlat() / 2].sea_level_temperature_c;
    assert!((20.0..32.0).contains(&eq), "equator {eq}");
    let gm = out.summary.global_mean_temperature_c;
    assert!((5.0..20.0).contains(&gm), "global {gm}");
    let p = out.summary.global_mean_precipitation_mm_per_yr;
    assert!((400.0..2000.0).contains(&p), "precip {p}");
    eprintln!(
        "aquaplanet: T_eq {eq:.2} C, global T {gm:.2} C, P {p:.0} mm/yr, ice edge N {:?}",
        out.summary.ice_edge_lat_deg_north
    );
    assert!(out.summary.ice_edge_lat_deg_north.is_some());
}

#[test]
fn climate_v0_envelope_refusals_are_typed_and_produce_no_output() {
    let g = grid();
    let (z, o) = world(&g);
    let go = |r: &ClimateRecipeV1| {
        derive_climate_coarse_v0(
            &ClimateInputs { grid: &g, elevation_m: &z, fill_level_m: &z, sea_level_m: 0.0, ocean: &o },
            r,
        )
    };
    let oob = |e: ClimateError, want: &str| match e {
        ClimateError::Recipe(UniverseError::OutOfEnvelope { constant, .. }) => {
            assert!(constant.contains(want), "{constant}")
        }
        other => panic!("expected OutOfEnvelope({want}), got {other:?}"),
    };
    let mut r = recipe();
    r.constants.axial_tilt_deg.value = 60.0;
    oob(go(&r).unwrap_err(), "axial_tilt_deg");
    let mut r = recipe();
    r.constants.rotation_period_s.value = 10.0 * 3600.0;
    oob(go(&r).unwrap_err(), "rotation_period_s");
    let mut r = recipe();
    r.constants.star_luminosity_w.value *= 2.0;
    oob(go(&r).unwrap_err(), "solar_constant_w_m2 (derived");
    let mut r = recipe();
    r.constants.lapse_rate_k_per_m.value = 0.012;
    oob(go(&r).unwrap_err(), "lapse_rate_k_per_m");
    let mut r = recipe();
    r.constants.land_recycling_ratio.value = 1.0;
    oob(go(&r).unwrap_err(), "land_recycling_ratio");
    let mut r = recipe();
    r.constants.trade_wind_speed_m_s.value = f64::NAN;
    oob(go(&r).unwrap_err(), "trade_wind_speed_m_s");
    // Unit mismatch.
    let mut r = recipe();
    r.constants.lapse_rate_k_per_m.unit = UnitV1::Meter;
    assert!(matches!(go(&r), Err(ClimateError::Recipe(UniverseError::UnitMismatch { .. }))));
    // Ordering.
    let mut r = recipe();
    r.constants.albedo_ice.value = 0.31;
    assert!(matches!(go(&r), Err(ClimateError::Recipe(UniverseError::Inconsistent(_)))));
    // Named but unimplemented law.
    let mut r = recipe();
    r.law = ClimateLawV1::SeedClimateV1;
    assert!(matches!(go(&r), Err(ClimateError::Recipe(UniverseError::Unsupported { .. }))));
    // Schema.
    let mut r = recipe();
    r.schema = "pf.climate-recipe/v2".into();
    assert!(matches!(go(&r), Err(ClimateError::Recipe(UniverseError::Schema { .. }))));
    // Inputs.
    let short = &z[..10];
    assert!(matches!(
        derive_climate_coarse_v0(
            &ClimateInputs { grid: &g, elevation_m: short, fill_level_m: &z, sea_level_m: 0.0, ocean: &o },
            &recipe()
        ),
        Err(ClimateError::Inputs(_))
    ));
    let mut zn = z.clone();
    zn[7] = f64::INFINITY;
    assert!(matches!(
        derive_climate_coarse_v0(
            &ClimateInputs { grid: &g, elevation_m: &zn, fill_level_m: &z, sea_level_m: 0.0, ocean: &o },
            &recipe()
        ),
        Err(ClimateError::NonFinite { field: "elevation_m", index: 7 })
    ));
    // Every spec id unique, bounds finite and ordered; earth_v0 validates.
    for (i, s) in CLIMATE_CONSTANT_SPECS_V1.iter().enumerate() {
        assert!(s.lo.is_finite() && s.hi.is_finite() && s.lo <= s.hi, "{}", s.id);
        assert!(CLIMATE_CONSTANT_SPECS_V1[i + 1..].iter().all(|t| t.id != s.id), "{}", s.id);
    }
    recipe().validate().unwrap();
}

#[test]
fn climate_v0_schema_round_trip_and_refusal_of_unknown_ids() {
    let r = recipe();
    let j = serde_json::to_string(&r).unwrap();
    assert!(j.contains("\"law\":\"climate-coarse-v0\""), "{j}");
    let back: ClimateRecipeV1 = serde_json::from_str(&j).unwrap();
    assert_eq!(back, r);
    assert_eq!(back.digest().unwrap(), r.digest().unwrap());
    for bad in [
        j.replace("climate-coarse-v0", "climate-magic-v9"),
        j.replace("\"lapse_rate_k_per_m\"", "\"lapse\""),
        j.replacen("\"unit\":\"K/m\"", "\"unit\":\"K/km\"", 1),
        j.replacen("\"class\":\"project_assumption_awaiting_review\"", "\"class\":\"vibes\"", 1),
    ] {
        assert_ne!(bad, j);
        let e = serde_json::from_str::<ClimateRecipeV1>(&bad).expect_err("accepted").to_string();
        assert!(e.contains("unknown variant") || e.contains("unknown field"), "{e}");
    }
    // Every constant of the v0 record is awaiting review (ADR 039 C5).
    let v = serde_json::to_value(ClimateConstantsV1::earth_v0()).unwrap();
    let n = v.as_object().unwrap().values().filter(|x| x.get("source").is_some()).count();
    assert_eq!(n, CLIMATE_CONSTANT_SPECS_V1.len());
    assert!(!j.contains("\"cited\""));

    // Domain payload: build, round-trip, schema check.
    let g = grid();
    let (z, o) = world(&g);
    let out = run(&g, &z, &o, &r);
    let d = ClimateDomainV0::build(
        &out,
        &r,
        ClimateInputRecordV0 {
            elevation_field: "terrain/elevation_m".into(),
            fill_level_field: "fixture".into(),
            sea_level_m: 0.0,
            ocean_source: "fixture".into(),
        },
    )
    .unwrap();
    d.validate_schema().unwrap();
    assert_eq!(d.tier.id, CLIMATE_TIER_ID);
    let dj = serde_json::to_string(&d).unwrap();
    let dback: ClimateDomainV0 = serde_json::from_str(&dj).unwrap();
    assert_eq!(dback, d);
    assert_eq!(serde_json::to_string(&dback).unwrap(), dj, "re-serialization is byte-identical");
    let e = serde_json::from_str::<ClimateDomainV0>(&dj.replacen("\"tier\"", "\"tear\"", 1))
        .expect_err("accepted");
    assert!(e.to_string().contains("unknown field") || e.to_string().contains("missing field"));
    // A later tier fills the same names: the schema check rejects a rename.
    let mut renamed = d.clone();
    renamed.fields[2].name = "rain_mm".into();
    assert!(renamed.validate_schema().is_err());
    // Approximation ids are unique and every field cites known ones.
    let ids: Vec<&str> = CLIMATE_APPROXIMATIONS_V0.iter().map(|a| a.id).collect();
    for (i, a) in ids.iter().enumerate() {
        assert!(!ids[i + 1..].contains(a), "{a}");
    }
    // Every constant read by an approximation exists in the spec table.
    for a in CLIMATE_APPROXIMATIONS_V0 {
        for c in a.constants {
            assert!(CLIMATE_CONSTANT_SPECS_V1.iter().any(|s| s.id == *c), "{} reads unknown {c}", a.id);
        }
    }
}

#[test]
fn climate_v0_lake_cells_sit_at_their_water_surface() {
    // A disconnected sub-sea basin (floor -3000 m) filled to a lake surface
    // at +50 m: its temperature follows the lake surface, not the floor.
    let g = grid();
    let (mut z, o) = world(&g);
    let nlon = g.nlon();
    let c = 22 * nlon + 30;
    assert!(!o[c]);
    let mut fill = z.clone();
    z[c] = -3000.0;
    fill[c] = 50.0;
    let lake = derive_climate_coarse_v0(
        &ClimateInputs { grid: &g, elevation_m: &z, fill_level_m: &fill, sea_level_m: 0.0, ocean: &o },
        &recipe(),
    )
    .unwrap();
    let t_sl = lake.rows[22].sea_level_temperature_c;
    let lapse = recipe().constants.lapse_rate_k_per_m.value;
    assert_eq!(lake.fields.temperature_mean_c[c], t_sl - lapse * 50.0);
    // Without the fill level the floor would be 19.5 K warmer.
    let dry = derive_climate_coarse_v0(
        &ClimateInputs { grid: &g, elevation_m: &z, fill_level_m: &z, sea_level_m: 0.0, ocean: &o },
        &recipe(),
    )
    .unwrap();
    assert!(dry.fields.temperature_mean_c[c] > lake.fields.temperature_mean_c[c] + 19.0);
    // Fill-level length is checked.
    assert!(matches!(
        derive_climate_coarse_v0(
            &ClimateInputs { grid: &g, elevation_m: &z, fill_level_m: &fill[..5], sea_level_m: 0.0, ocean: &o },
            &recipe()
        ),
        Err(ClimateError::Inputs(_))
    ));
}
