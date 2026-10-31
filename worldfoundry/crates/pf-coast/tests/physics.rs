use pf_coast::*;
fn c() -> Config {
    Config::new(2.0).unwrap()
}
fn p() -> Point {
    Point {
        height_m: 0.5,
        slope: 0.02,
        litho: 10,
        shore_m: 1.0,
        fixed: 0,
        shore_path_m: 0.0,
        shore_connected: true,
        tide_connected: true,
        touches_marine: true,
        mouth_area_m2: 0.0,
        wave_height_m: 1.0,
        wave_period_s: 5.0,
    }
}
#[test]
fn wave_limits_and_monotonicity() {
    assert_eq!(wind_wave(0.0, 1e5), (0.0, 0.0));
    assert_eq!(wind_wave(10.0, 0.0), (0.0, 0.0));
    let mut old = 0.0;
    for f in [32.0, 256.0, 2048.0, 16384.0, 128000.0] {
        let (h, t) = wind_wave(10.0, f);
        assert!(h > old && t > 0.0);
        old = h;
    }
    assert!(wind_wave(20.0, 10000.0).0 > wind_wave(10.0, 10000.0).0);
    let ua = 0.71 * 10_f64.powf(1.23);
    assert!((wind_wave(10.0, 1e20).0 - 0.283 * ua * ua / G).abs() < 1e-12);
}
#[test]
fn zero_wave_runup_and_dimensional_scaling() {
    assert_eq!(runup_m(0.0, 4.0, 0.05), 0.0);
    assert!((runup_m(4.0, 8.0, 0.05) / runup_m(1.0, 4.0, 0.05) - 4.0).abs() < 1e-12);
}
#[test]
fn supply_zero_monotonic_and_bounded() {
    assert_eq!(supply_index(0.0, 3.0, &c()), 0.0);
    for a in [1.0, 1e3, 1e6, 1e9] {
        let s = supply_index(a, 0.5, &c());
        assert!((0.0..1.0).contains(&s));
        assert!(supply_index(a, 3.0, &c()) > s);
    }
}
#[test]
fn meaningful_class_falsifiers() {
    let base = p();
    assert_eq!(classify(base, &c()).unwrap().class, 4);
    assert_eq!(
        classify(Point { slope: 1.0, ..base }, &c()).unwrap().class,
        1
    );
    assert_eq!(classify(Point { litho: 3, ..base }, &c()).unwrap().class, 2);
    assert_eq!(
        classify(Point { slope: 0.2, ..base }, &c()).unwrap().class,
        3
    );
    assert_eq!(
        classify(
            Point {
                wave_height_m: 0.01,
                height_m: 0.0,
                ..base
            },
            &c()
        )
        .unwrap()
        .class,
        5
    );
    let m = classify(
        Point {
            wave_height_m: 0.01,
            ..base
        },
        &c(),
    )
    .unwrap();
    assert_eq!(m.class, 6);
    assert_ne!(m.flags & MARSH_POTENTIAL, 0);
    let d = classify(
        Point {
            mouth_area_m2: 2e6,
            ..base
        },
        &c(),
    )
    .unwrap();
    assert_eq!(d.class, 7);
    assert_ne!(d.flags & DELTA_POTENTIAL, 0);
    assert_ne!(
        classify(
            Point {
                mouth_area_m2: 2e6,
                slope: 0.2,
                ..base
            },
            &c()
        )
        .unwrap()
        .class,
        7
    );
}
#[test]
fn water_and_inland_are_not_fabricated_beaches() {
    for fixed in [2, 4] {
        assert_eq!(
            classify(
                Point {
                    fixed,
                    litho: 255,
                    ..p()
                },
                &c()
            )
            .unwrap()
            .class,
            0
        );
    }
    assert_eq!(
        classify(
            Point {
                shore_m: 513.0,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        0
    );
    assert_eq!(
        classify(
            Point {
                height_m: 20.0,
                shore_m: 40.0,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        0
    );
}
#[test]
fn invalid_inputs_fail() {
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(Config::new(bad).is_err());
        assert!(classify(Point { slope: bad, ..p() }, &c()).is_err());
    }
    assert!(classify(Point { litho: 255, ..p() }, &c()).is_err());
}
#[test]
fn flat_width_is_censored_not_infinite() {
    let o = classify(Point { slope: 0.0, ..p() }, &c()).unwrap();
    assert_eq!(o.width_m, 512.0);
    assert_ne!(o.flags & WIDTH_CENSORED, 0);
}
#[test]
fn thresholds_are_not_site_tuned() {
    let threshold = 30_f64.to_radians().tan();
    assert_eq!(
        classify(
            Point {
                slope: threshold,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        1
    );
    assert_ne!(
        classify(
            Point {
                slope: threshold - 1e-8,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        1
    );
    let base = Point {
        wave_height_m: 2.0 / 15.0,
        height_m: 0.0,
        ..p()
    };
    assert_ne!(classify(base, &c()).unwrap().class, 5);
    assert_eq!(
        classify(
            Point {
                wave_height_m: base.wave_height_m - 1e-8,
                ..base
            },
            &c()
        )
        .unwrap()
        .class,
        5
    );
}

#[test]
fn marsh_uses_pioneer_elevation_and_flat_never_extends_above_high_water() {
    let lo = marsh_lower_m(&c());
    assert!((lo - 0.14169523569287323).abs() < 1e-12);
    let p = Point {
        wave_height_m: 0.01,
        height_m: lo - 0.001,
        ..p()
    };
    assert_eq!(classify(p, &c()).unwrap().class, 5);
    assert_eq!(
        classify(Point { height_m: lo, ..p }, &c()).unwrap().class,
        6
    );
    assert_ne!(
        classify(Point { height_m: 1.1, ..p }, &c()).unwrap().class,
        5
    );
    assert_ne!(
        classify(Point { height_m: 1.1, ..p }, &c()).unwrap().class,
        6
    );
}
#[test]
fn water_channels_detached_and_subtidal_posts_cannot_be_coast() {
    for fixed in [1, 2, 4, 5] {
        assert_eq!(classify(Point { fixed, ..p() }, &c()).unwrap().class, 0);
    }
    assert_eq!(
        classify(
            Point {
                shore_connected: false,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        0
    );
    assert_eq!(
        classify(
            Point {
                height_m: -1.01,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        0
    );
    assert_eq!(
        classify(
            Point {
                shore_path_m: 513.0,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        0
    );
}
#[test]
fn sediment_labels_have_no_first_post_elevation_exception() {
    assert_eq!(
        classify(
            Point {
                height_m: 8.0,
                slope: 0.2,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        0
    );
    assert_eq!(
        classify(
            Point {
                height_m: 8.0,
                slope: 0.8,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        1
    );
    assert_eq!(
        classify(
            Point {
                height_m: 8.0,
                slope: 0.8,
                touches_marine: false,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        0
    );
}

#[test]
fn low_ground_behind_barrier_is_not_intertidal() {
    assert_eq!(
        classify(
            Point {
                tide_connected: false,
                ..p()
            },
            &c()
        )
        .unwrap()
        .class,
        0
    );
    assert_ne!(classify(p(), &c()).unwrap().class, 0);
}
