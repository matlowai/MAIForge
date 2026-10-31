use pf_surface::climate::seed::seasons::*;
#[test]
fn thermal_coast_inertia_and_continental_limit() {
    let coast = thermal_response(152.0, 4.94, 0.0, 1e6);
    let inland = thermal_response(152.0, 4.94, 3e6, 1e6);
    let dry = thermal_response(152.0, 4.94, -1.0, 1e6);
    assert!(coast.0 > 3.0 && coast.0 < 4.0);
    assert!(coast.0 < inland.0 && inland.0 < dry.0);
    assert!(coast.1 > inland.1 && inland.1 > dry.1);
    assert_eq!(thermal_response(0.0, 4.94, 0.0, 1e6).0, 0.0);
}
#[test]
fn monthly_budget_and_hemisphere_falsifiers() {
    for lat in [-89., -55., -35., 0., 35., 55., 89.] {
        for ocean in [0., 0.5, 1.] {
            let f = winter_target(lat, ocean);
            let a = monthly_shares(lat, f);
            assert!(a.iter().all(|x| *x > 0.0));
            assert!((a.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            let cold: f64 = a
                .iter()
                .enumerate()
                .filter(|(m, _)| winter_month(lat >= 0.0, *m))
                .map(|(_, x)| x)
                .sum();
            assert!((cold - f).abs() < 1e-14);
        }
    }
    assert!(winter_target(35., 1.) > 0.7);
    assert!(winter_target(55., 0.) < 0.5);
    assert!(monthly_shares(55., 0.6)[0] > monthly_shares(55., 0.6)[6]);
    assert!(monthly_shares(-55., 0.6)[0] < monthly_shares(-55., 0.6)[6]);
}
#[test]
fn orographic_reference_supply_and_dry_falsifiers() {
    let f = |z| precipitation_at_height(1300., 10., 8., 500., 0.000455, 1600., z);
    assert_eq!(f(1600.), 1300.);
    assert!(f(800.) < f(1600.));
    assert!(f(2000.) > f(1600.));
    assert!(f(100000.) <= 5300. && f(0.) >= 0.);
    assert_eq!(
        precipitation_at_height(300., 1., 8., 500., 0.000455, 1600., 800.),
        300.
    );
    let derivative = (f(1600.001) - f(1600. - 0.001)) / 0.002;
    assert!((derivative - 500. * 0.000455 * 8.).abs() < 1e-6);
}

#[test]
fn independent_flags_and_typed_legacy_roundtrip() {
    use pf_grid::GridSpec;
    use pf_surface::climate::seed::seed_climate_v1_earth_recipe;
    use pf_surface::climate::{
        derive_climate, ClimateDomainV0, ClimateInputRecordV0, ClimateInputs,
    };
    let g = GridSpec::new(12, 24, 6.371e6).unwrap();
    let z = vec![0.; g.cells()];
    let ocean = vec![true; g.cells()];
    let recipe = seed_climate_v1_earth_recipe();
    let out = derive_climate(
        &ClimateInputs {
            grid: &g,
            elevation_m: &z,
            fill_level_m: &z,
            sea_level_m: 0.,
            ocean: &ocean,
        },
        &recipe,
    )
    .unwrap();
    for (opt, count) in [
        (Options::default(), 0),
        (
            Options {
                maritime: true,
                ..Options::default()
            },
            3,
        ),
        (
            Options {
                precipitation_seasons: true,
                ..Options::default()
            },
            13,
        ),
        (
            Options {
                orographic: true,
                ..Options::default()
            },
            5,
        ),
    ] {
        assert_eq!(
            derive(&out, &recipe, g.nlon(), &z, opt).unwrap().len(),
            count
        );
    }
    let dom = ClimateDomainV0::build(
        &out,
        &recipe,
        ClimateInputRecordV0 {
            elevation_field: "terrain/elevation_m".into(),
            fill_level_field: "hydrology/fill_level_m".into(),
            sea_level_m: 0.,
            ocean_source: "test instrument".into(),
        },
    )
    .unwrap();
    let serialized = serde_json::to_value(&dom).unwrap();
    assert!(serialized.get("seasonal_extension").is_none());
    let back: ClimateDomainV0 = serde_json::from_value(serialized).unwrap();
    assert_eq!(dom, back);
}
