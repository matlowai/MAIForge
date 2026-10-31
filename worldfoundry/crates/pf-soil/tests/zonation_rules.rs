//! Rule-level tests of pf-soil (lane SOIL-VEG-1): exact maths of the seasonal cycle, the Ohmura inversion, the
//! insolation geometry, monotone zonation with height, determinism of the threaded map.

use pf_soil::params::SoilParams;
use pf_soil::tier32::par_map;
use pf_soil::zonation::*;

#[test]
fn season_matches_closed_form() {
    // T = 0 + 10 cos: the half of the year above 0 C, mean 10 / (pi / 2)
    let (f, m) = season(0.0, 10.0, 0.0);
    assert!((f - 0.5).abs() < 1e-12);
    assert!((m - 20.0 / std::f64::consts::PI).abs() < 1e-12);
    assert_eq!(season(5.0, 1.0, 0.9).0, 1.0);
    assert_eq!(season(-5.0, 1.0, 0.9).0, 0.0);
}

#[test]
fn treeline_threshold_is_the_boundary() {
    let p = SoilParams::v1();
    for a in [0.0, 3.0, 6.6, 12.0, 20.0] {
        let t = tm_treeline(&p, a);
        assert!(tree_ok(&p, t + 1e-6, a) && !tree_ok(&p, t - 1e-3, a), "a {a} t {t}");
    }
    // no seasonality: the threshold is the 6.4 C growing-season mean itself
    assert!((tm_treeline(&p, 0.0) - 6.4).abs() < 1e-9);
}

#[test]
fn ohmura_inverse_round_trips() {
    let p = SoilParams::v1();
    for pmm in [300.0, 645.0, 1306.0, 3019.0] {
        let t = ohmura_t_ela(&p, pmm);
        let back = p.ohm_a + p.ohm_b * t + p.ohm_c * t * t;
        assert!((back - pmm).abs() < 1e-6, "{pmm} -> {t} -> {back}");
    }
    assert!(ohmura_t_ela(&p, 645.0).abs() < 1e-9);
}

#[test]
fn insolation_geometry() {
    let p = SoilParams::v1();
    let lat = 45f64.to_radians();
    assert!((insolation_ratio(&p, lat, 0.0, 0.0) - 1.0).abs() < 1e-12);
    let t = 30f64.to_radians().tan();
    // northern hemisphere, warm half year: a 30 deg south-facing slope (z rises to the north, zy > 0) gains a little
    // (sun already high: measured 1.040), a north-facing one loses much more (0.713)
    let south = insolation_ratio(&p, lat, 0.0, t);
    let north = insolation_ratio(&p, lat, 0.0, -t);
    assert!(south > 1.02 && north < 0.8, "south {south} north {north}");
    // the southern hemisphere mirrors it
    let s2 = insolation_ratio(&p, -lat, 0.0, -t);
    assert!((s2 - south).abs() < 1e-9);
    // east and west faces are equal
    assert!((insolation_ratio(&p, lat, t, 0.0) - insolation_ratio(&p, lat, -t, 0.0)).abs() < 1e-9);
}

#[test]
fn cavm_curve_is_continuous_and_monotone() {
    let p = SoilParams::v1();
    let mut last = -1.0;
    for k in 0..500 {
        let s = k as f64 * 0.1;
        let c = cavm_cover(&p, s);
        assert!(c >= last - 1e-12 && (0.0..=1.0).contains(&c));
        last = c;
    }
    assert_eq!(cavm_cover(&p, 35.0), 1.0);
}

fn post(h: f64) -> PostIn {
    PostIn { h, sea: 0.0, t_sl: 5.19, half: 6.6, p: 1306.0, ai: 5.0, lapse: 0.0065, r_sun: 1.0, slope_deg: 10.0, regime: 0, twi: 6.0,
        f_age: 1.0, f_kill: 1.0, soil_age_f: 1.0, flat_slope_deg: 10.0, ..Default::default() }
}

#[test]
fn cover_falls_with_height_through_the_zones() {
    let p = SoilParams::v1();
    let mut last = 2.0;
    let mut classes = vec![];
    for k in 0..90 {
        let o = zone(&p, &post(50.0 * k as f64));
        assert!(o.cover <= last + 1e-12, "cover rose at {} m", 50 * k);
        last = o.cover;
        if classes.last() != Some(&o.class) { classes.push(o.class); }
    }
    assert_eq!(classes.first(), Some(&C_CLOSED));
    assert!(classes.contains(&C_ECOTONE) && classes.contains(&C_NIVAL));
    assert_eq!(classes.last(), Some(&C_NIVAL));
}

#[test]
fn par_map_is_order_and_thread_invariant() {
    let f = |k: usize| ((k as f64).sin() * 1e6).round() as i64;
    let a = par_map(10_007, 1, f);
    let b = par_map(10_007, 7, f);
    let c = par_map(10_007, 32, f);
    assert_eq!(a, b);
    assert_eq!(a, c);
}
