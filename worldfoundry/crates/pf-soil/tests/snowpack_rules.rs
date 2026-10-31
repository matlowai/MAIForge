//! Rule-level tests of the SNOWPACK-1 module (lane SNOWPACK-1): exact maths of the normal-distribution helpers, the
//! irradiance geometry, the pack's monotone response to height and aspect, periodic spin-up, determinism.

use pf_soil::snowpack::*;

#[test]
fn erfc_and_normal_helpers() {
    // reference values of erfc (Abramowitz and Stegun table values)
    for (x, r) in [(0.0, 1.0), (0.5, 0.4795001222), (1.0, 0.1572992071), (2.0, 0.004677734981), (-1.0, 1.842700793)] {
        assert!((erfc(x) - r).abs() < 2e-7 * r.max(1e-3), "erfc({x}) = {} vs {r}", erfc(x));
    }
    assert!((ncdf(0.0) - 0.5).abs() < 1e-7); // erfcc fractional error < 1.2e-7
    assert!((ncdf(1.96) - 0.9750021).abs() < 1e-6);
    // PDD: s phi(0) at t = 0; -> t for t >> s; -> 0 for t << -s; even-odd identity pdd(t) - pdd(-t) = t
    let s = 3.5;
    assert!((pdd(0.0, s) - s / (std::f64::consts::TAU).sqrt()).abs() < 1e-7);
    assert!((pdd(30.0, s) - 30.0).abs() < 1e-6);
    assert!(pdd(-30.0, s) < 1e-9);
    for t in [-4.0, -1.0, 0.5, 2.0, 7.0] { assert!((pdd(t, s) - pdd(-t, s) - t).abs() < 1e-6, "t {t}"); }
}

#[test]
fn snow_share_limits() {
    let p = SnowParams::v1();
    assert!(snow_share(&p, -30.0) > 0.999999);
    assert!(snow_share(&p, 30.0) < 1e-6);
    // monotone decreasing
    let mut prev = 2.0;
    for k in -40..40 { let s = snow_share(&p, k as f64 * 0.5); assert!(s <= prev + 1e-12); prev = s; }
}

#[test]
fn irradiance_geometry() {
    let p = SnowParams::v1();
    let lat = 60f64.to_radians();
    let dec = 23.44f64.to_radians();
    let flat = daily_direct(&p, lat, 0.0, 0.0, dec, 1.0);
    let (nx, ny) = facing(30.0, 0.0);
    let (sx, sy) = facing(30.0, 180.0);
    let north = daily_direct(&p, lat, nx, ny, dec, 1.0);
    let south = daily_direct(&p, lat, sx, sy, dec, 1.0);
    assert!(south > flat && flat > north, "N hemisphere: south-facing gets more ({south} {flat} {north})");
    // mirror: southern hemisphere at the opposite declination, north-facing equals the northern south-facing
    let (mx, my) = facing(30.0, 0.0);
    let m = daily_direct(&p, -lat, mx, my, -dec, 1.0);
    assert!((m - south).abs() < 1e-9 * south.max(1.0));
    // polar night: nothing
    assert_eq!(daily_direct(&p, 80f64.to_radians(), 0.0, 0.0, -dec, 1.0), 0.0);
    // higher ground (lower pressure) gets more direct beam
    assert!(daily_direct(&p, lat, 0.0, 0.0, dec, 0.6) > flat);
}

#[test]
fn pack_responds_to_height_and_aspect() {
    let p = SnowParams::v1();
    let cal = Calendar::new(&p);
    let lat = -54.75f64.to_radians();
    let days: Vec<usize> = (0..DAYS).collect();
    let at = |h: f64, zx: f64, zy: f64| post(&p, &cal, &SnowIn { lat_rad: lat, t_sl: 5.15, half: 6.57, p_mm: 1257.0, lapse: 0.0065, h, sea: 0.0, zx, zy }, &days);
    let mut prev_b = f64::NEG_INFINITY;
    let mut prev_days = 0u16;
    for k in 0..30 {
        let o = at(k as f64 * 100.0, 0.0, 0.0);
        assert!(o.balance_mm > prev_b && o.snow_days >= prev_days, "monotone in height at {k}");
        prev_b = o.balance_mm;
        prev_days = o.snow_days;
    }
    // southern hemisphere: a south-facing (poleward) slope keeps more snow than a north-facing one
    let (sx, sy) = facing(30.0, 180.0);
    let (nx, ny) = facing(30.0, 0.0);
    let s = at(1200.0, sx, sy);
    let n = at(1200.0, nx, ny);
    assert!(s.balance_mm > n.balance_mm && s.snow_days > n.snow_days && s.irr_ratio < 1.0 && n.irr_ratio > 1.0);
    // determinism
    assert_eq!(at(1300.0, 0.1, -0.2), at(1300.0, 0.1, -0.2));
}

#[test]
fn spin_up_is_periodic_below_the_line() {
    let mut p = SnowParams::v1();
    let cal = Calendar::new(&p);
    let days: Vec<usize> = (0..DAYS).collect();
    let x = SnowIn { lat_rad: 45f64.to_radians(), t_sl: 8.0, half: 9.0, p_mm: 1200.0, lapse: 0.0065, h: 800.0, sea: 0.0, zx: 0.0, zy: 0.0 };
    let a = post(&p, &cal, &x, &days);
    assert!(!a.perennial && a.min_swe == 0.0 && a.max_swe > 0.0);
    p.spin_years = 6.0;
    let b = post(&p, &cal, &x, &days);
    assert_eq!(a.swe, b.swe);
}

#[test]
fn ela_brackets_the_sign_change() {
    let p = SnowParams::v1();
    let cal = Calendar::new(&p);
    let lat = 45f64.to_radians();
    let e = ela(&p, &cal, lat, 6.0, 8.0, 1500.0, 0.0065, 0.0, 0.0, -6000.0, 12000.0).expect("an equilibrium line");
    assert!(balance_at(&p, &cal, lat, 6.0, 8.0, 1500.0, 0.0065, e - 5.0, 0.0, 0.0) < 0.0);
    assert!(balance_at(&p, &cal, lat, 6.0, 8.0, 1500.0, 0.0065, e + 5.0, 0.0, 0.0) > 0.0);
}

#[test]
fn ohmura_copies_match_pf_soil() {
    let s = pf_soil::params::SoilParams::v1();
    let p = SnowParams::v1();
    assert_eq!((p.ohm_a, p.ohm_b, p.ohm_c, p.jja_factor), (s.ohm_a, s.ohm_b, s.ohm_c, s.jja_factor));
}

#[test]
fn anchored_line_is_the_ohmura_line() {
    let p = SnowParams::v1();
    let cal = Calendar::new(&p);
    let s = pf_soil::params::SoilParams::v1();
    for (lat, t_sl, half, pmm) in [(-54.75, 5.15, 6.57, 1257.0), (44.1, 12.06, 5.44, 2836.0), (62.0, 2.0, 12.0, 600.0)] {
        let la = f64::to_radians(lat);
        let ohm = (t_sl + s.jja_factor * half - pf_soil::zonation::ohmura_t_ela(&s, pmm)) / 0.0065;
        let e = ela(&p, &cal, la, t_sl, half, pmm, 0.0065, 0.0, 0.0, -6000.0, 12000.0).unwrap();
        assert!((e - ohm).abs() < 0.5, "lat {lat}: anchored {e} vs Ohmura {ohm}");
    }
}
