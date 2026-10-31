//! SEED-CLIMATE (lane CM, 2026-10-12): tests of the `seed-climate-v1` tier.
//! Every map here is a SYNTHETIC FIXTURE (ADR 026 §1: legal as instrument
//! calibration in unit tests, never world input).

use pf_contracts::laws::{ClimateLawV1, ClimateRecipeV1, UniverseError};
use pf_grid::GridSpec;
use pf_surface::climate::seed::*;
use pf_surface::climate::*;

const R: f64 = 6.371e6;

fn inputs<'a>(g: &'a GridSpec, z: &'a [f64], ocean: &'a [bool]) -> ClimateInputs<'a> {
    ClimateInputs { grid: g, elevation_m: z, fill_level_m: z, sea_level_m: 0.0, ocean }
}

fn v1(g: &GridSpec, z: &[f64], ocean: &[bool]) -> ClimateOutputV0 {
    derive_climate(&inputs(g, z, ocean), &seed_climate_v1_earth_recipe()).unwrap()
}

fn v0(g: &GridSpec, z: &[f64], ocean: &[bool]) -> ClimateOutputV0 {
    derive_climate_coarse_v0(&inputs(g, z, ocean), &ClimateRecipeV1::coarse_v0_earth()).unwrap()
}

/// A continent spanning `lon_cols` at every latitude below 70 deg, so many
/// latitude rings have no ocean; a N-S ridge `wall_m` high from `wall_lon`
/// to `wall_lon` + 6 deg (its windward rise is one cell at any resolution
/// used here: a step, not a resolved ramp).
fn band_world(g: &GridSpec, lon_lo_deg: f64, lon_hi_deg: f64, wall_lon_deg: Option<f64>, wall_m: f64) -> (Vec<f64>, Vec<bool>) {
    let (nlat, nlon) = (g.nlat(), g.nlon());
    let mut z = vec![-4000.0; g.cells()];
    let mut ocean = vec![true; g.cells()];
    for i in 0..nlat {
        let lat = g.center_lat_deg(i).unwrap();
        for j in 0..nlon {
            let lon = g.center_lon_deg(j).unwrap();
            if lon >= lon_lo_deg && lon < lon_hi_deg && lat.abs() < 70.0 {
                let c = i * nlon + j;
                ocean[c] = false;
                z[c] = 200.0;
                if let Some(w) = wall_lon_deg {
                    if lon >= w && lon < w + 6.0 {
                        z[c] = 200.0 + wall_m;
                    }
                }
            }
        }
    }
    (z, ocean)
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn seed_v1_determinism_and_closure() {
    let g = GridSpec::new(36, 72, R).unwrap();
    let (z, o) = band_world(&g, -60.0, 40.0, Some(-20.0), 3000.0);
    let a = v1(&g, &z, &o);
    let b = v1(&g, &z.clone(), &o.clone());
    for ((n, x), (_, y)) in a.fields.named().iter().zip(b.fields.named().iter()) {
        assert_eq!(bits(x), bits(y), "field {n} not bit-identical");
    }
    for c in &a.closures {
        assert!(c.closes(), "{c:?}");
    }
    let w = &a.closures[0];
    assert!(w.residual.abs() <= SEED_WATER_REL_TOL * w.gross.abs() + 1.0, "{w:?}");
    let rec = a.seed_v1.as_ref().unwrap();
    assert!(rec.sweeps > 1 && rec.sweeps < SEED_MAX_SWEEPS);
}

/// v0's defect: rings with no ocean get zero rain. v1 brings water in
/// meridionally, so a 100-deg-wide continent has rain in every land row.
#[test]
fn seed_v1_land_rings_without_ocean_get_rain() {
    let g = GridSpec::new(36, 72, R).unwrap();
    // A continent covering all longitudes between 20 and 50 deg N: those
    // rings have no ocean at all.
    let (nlat, nlon) = (g.nlat(), g.nlon());
    let mut z = vec![-4000.0; g.cells()];
    let mut o = vec![true; g.cells()];
    for i in 0..nlat {
        let lat = g.center_lat_deg(i).unwrap();
        if (20.0..50.0).contains(&lat) {
            for j in 0..nlon {
                z[i * nlon + j] = 300.0;
                o[i * nlon + j] = false;
            }
        }
    }
    let old = v0(&g, &z, &o);
    let new = v1(&g, &z, &o);
    assert!(old.summary.rows_without_ocean > 0);
    let mut v0_zero_rows = 0;
    for i in 0..nlat {
        let lat = g.center_lat_deg(i).unwrap();
        if !(20.0..50.0).contains(&lat) {
            continue;
        }
        let p0 = old.fields.precipitation_mm_per_yr[i * nlon];
        let p1 = new.fields.precipitation_mm_per_yr[i * nlon];
        if p0 < 1.0 {
            v0_zero_rows += 1;
        }
        assert!(p1 > 100.0, "lat {lat}: seed-climate-v1 P {p1} mm/yr on a ring without ocean");
        assert!(new.fields.upwind_coast_distance_m[i * nlon] > 0.0);
    }
    assert!(v0_zero_rows > 0, "fixture no longer shows the v0 defect");
}

/// v0's defect: the saturation excess falls in the first wall cell, so the
/// windward maximum scales as 1 / cell size. v1 removes it over a declared
/// fallout length: halving the cell leaves the maximum nearly unchanged.
#[test]
fn seed_v1_windward_maximum_converges_with_resolution() {
    let max_land_p = |g: &GridSpec, f: &dyn Fn(&GridSpec, &[f64], &[bool]) -> ClimateOutputV0| {
        let (z, o) = band_world(g, -60.0, 40.0, Some(-20.0), 3000.0);
        let out = f(g, &z, &o);
        (0..g.cells()).filter(|&c| !o[c]).map(|c| out.fields.precipitation_mm_per_yr[c]).fold(0.0, f64::max)
    };
    let coarse = GridSpec::new(90, 180, R).unwrap(); // 2 deg, ~222 km
    let fine = GridSpec::new(180, 360, R).unwrap(); // 1 deg, ~111 km
    let r1 = max_land_p(&fine, &v1) / max_land_p(&coarse, &v1);
    let r0 = max_land_p(&fine, &v0) / max_land_p(&coarse, &v0);
    eprintln!("max windward P fine/coarse: v1 {r1:.3}, v0 {r0:.3}");
    assert!((0.67..1.5).contains(&r1), "seed-climate-v1 fine/coarse max P ratio {r1}");
    assert!(r0 > 1.6, "fixture no longer shows the v0 1/cell scaling (ratio {r0})");
}

/// Rain shadow: in the westerlies the windward (west) flank of a wall is
/// wetter than the same longitude without the wall, the lee is drier.
#[test]
fn seed_v1_rain_shadow_known_answer() {
    let g = GridSpec::new(90, 180, R).unwrap();
    let (zw, o) = band_world(&g, -60.0, 40.0, Some(-20.0), 3000.0);
    let (zf, _) = band_world(&g, -60.0, 40.0, None, 0.0);
    let wall = v1(&g, &zw, &o);
    let flat = v1(&g, &zf, &o);
    let i = (0..g.nlat()).find(|&i| (g.center_lat_deg(i).unwrap() - 45.0).abs() < 1.1).unwrap();
    let col = |lon: f64| (0..g.nlon()).find(|&j| (g.center_lon_deg(j).unwrap() - lon).abs() < 1.1).unwrap();
    let windward = i * g.nlon() + col(-19.0);
    let lee = i * g.nlon() + col(15.0);
    let p = |o: &ClimateOutputV0, c: usize| o.fields.precipitation_mm_per_yr[c];
    assert!(p(&wall, windward) > 1.5 * p(&flat, windward), "windward {} vs flat {}", p(&wall, windward), p(&flat, windward));
    assert!(p(&wall, lee) < 0.8 * p(&flat, lee), "lee {} vs flat {}", p(&wall, lee), p(&flat, lee));
    assert!(wall.fields.precipitation_orographic_mm_per_yr[windward] > 0.0);
}

/// Meridional transport: on an aquaplanet the equator (convergence) is
/// wetter than 30 deg (divergence), and every field is zonally uniform.
#[test]
fn seed_v1_aquaplanet_itcz_and_zonal_symmetry() {
    let g = GridSpec::new(36, 72, R).unwrap();
    let z = vec![-4000.0; g.cells()];
    let o = vec![true; g.cells()];
    let out = v1(&g, &z, &o);
    let nlon = g.nlon();
    let row_p = |lat: f64| {
        let i = (0..g.nlat()).min_by(|&a, &b| {
            let da = (g.center_lat_deg(a).unwrap() - lat).abs();
            let db = (g.center_lat_deg(b).unwrap() - lat).abs();
            da.partial_cmp(&db).unwrap()
        }).unwrap();
        out.fields.precipitation_mm_per_yr[i * nlon]
    };
    assert!(row_p(2.5) > 1.5 * row_p(27.5), "equator {} vs 30 deg {}", row_p(2.5), row_p(27.5));
    assert!(row_p(-2.5) > 1.5 * row_p(-27.5));
    for i in 0..g.nlat() {
        let row = &out.fields.precipitation_mm_per_yr[i * nlon..(i + 1) * nlon];
        let (lo, hi) = row.iter().fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
        assert!(hi - lo <= 1.0e-6 * hi.max(1.0), "row {i} not zonally uniform: {lo}..{hi}");
    }
    let gp = out.summary.global_mean_precipitation_mm_per_yr;
    assert!((500.0..2000.0).contains(&gp), "aquaplanet global P {gp}");
}

#[test]
fn seed_v1_dispatch_and_typed_refusals() {
    let g = GridSpec::new(18, 36, R).unwrap();
    let (z, o) = band_world(&g, -60.0, 40.0, Some(-20.0), 3000.0);
    // v0 through the dispatcher is bit-identical to the v0 entry point.
    let a = derive_climate(&inputs(&g, &z, &o), &ClimateRecipeV1::coarse_v0_earth()).unwrap();
    let b = v0(&g, &z, &o);
    assert_eq!(a, b);
    assert!(a.seed_v1.is_none());
    // The v0 entry point still refuses the seed law.
    assert!(matches!(
        derive_climate_coarse_v0(&inputs(&g, &z, &o), &seed_climate_v1_earth_recipe()),
        Err(ClimateError::Recipe(UniverseError::Unsupported { .. }))
    ));
    // The seed entry point refuses the v0 law.
    assert!(matches!(
        derive_seed_climate_v1(&inputs(&g, &z, &o), &ClimateRecipeV1::coarse_v0_earth(), &SeedClimateConstantsV1::earth_v1()),
        Err(ClimateError::Recipe(UniverseError::Unsupported { .. }))
    ));
    // Extra constants: envelope, unit, schema.
    let r = seed_climate_v1_earth_recipe();
    let mut e = SeedClimateConstantsV1::earth_v1();
    e.land_et_ratio_cap.value = 1.0;
    assert!(matches!(derive_seed_climate_v1(&inputs(&g, &z, &o), &r, &e), Err(ClimateError::Recipe(UniverseError::OutOfEnvelope { .. }))));
    let mut e = SeedClimateConstantsV1::earth_v1();
    e.orographic_fallout_length_m.value = f64::NAN;
    assert!(matches!(derive_seed_climate_v1(&inputs(&g, &z, &o), &r, &e), Err(ClimateError::Recipe(UniverseError::OutOfEnvelope { .. }))));
    let mut e = SeedClimateConstantsV1::earth_v1();
    e.eddy_moisture_diffusivity_m2_per_yr.unit = pf_contracts::laws::UnitV1::Meter;
    assert!(matches!(derive_seed_climate_v1(&inputs(&g, &z, &o), &r, &e), Err(ClimateError::Recipe(UniverseError::UnitMismatch { .. }))));
    let mut e = SeedClimateConstantsV1::earth_v1();
    e.schema = "pf.seed-climate-constants/v0".into();
    assert!(matches!(derive_seed_climate_v1(&inputs(&g, &z, &o), &r, &e), Err(ClimateError::Recipe(UniverseError::Schema { .. }))));
    // Base constants keep the v0 envelopes.
    let mut r2 = seed_climate_v1_earth_recipe();
    r2.constants.axial_tilt_deg.value = 60.0;
    assert!(matches!(derive_climate(&inputs(&g, &z, &o), &r2), Err(ClimateError::Recipe(UniverseError::OutOfEnvelope { .. }))));
    assert_eq!(r.law, ClimateLawV1::SeedClimateV1);
}

#[test]
fn seed_v1_payload_schema_round_trip() {
    let g = GridSpec::new(18, 36, R).unwrap();
    let (z, o) = band_world(&g, -60.0, 40.0, Some(-20.0), 3000.0);
    let rec = || ClimateInputRecordV0 {
        elevation_field: "terrain/elevation_m".into(),
        fill_level_field: "hydrology/fill_level_m".into(),
        sea_level_m: 0.0,
        ocean_source: "fixture".into(),
    };
    let r = seed_climate_v1_earth_recipe();
    let out = v1(&g, &z, &o);
    let dom = ClimateDomainV0::build(&out, &r, rec()).unwrap();
    dom.validate_schema().unwrap();
    assert_eq!(dom.tier.id, SEED_CLIMATE_TIER_ID);
    let wv = dom.fields.iter().find(|f| f.name == "wind_v_m_s").unwrap();
    assert_eq!(wv.status, FieldStatusV0::Modelled);
    let s = dom.seed_v1.as_ref().unwrap();
    assert_eq!(s.constants_digest, SeedClimateConstantsV1::earth_v1().digest().unwrap());
    let js = serde_json::to_string(&dom).unwrap();
    let back: ClimateDomainV0 = serde_json::from_str(&js).unwrap();
    assert_eq!(back, dom);
    assert_eq!(serde_json::to_string(&back).unwrap(), js);
    // A v0 payload carries no seed record and no `seed_v1` key at all.
    let r0 = ClimateRecipeV1::coarse_v0_earth();
    let d0 = ClimateDomainV0::build(&v0(&g, &z, &o), &r0, rec()).unwrap();
    assert_eq!(d0.tier.id, CLIMATE_TIER_ID);
    assert!(!serde_json::to_string(&d0).unwrap().contains("seed_v1"));
}
