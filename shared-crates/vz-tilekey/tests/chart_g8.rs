//! ADR 044 G8 (chart precision): chart -> geodetic -> chart within 1 mm over
//! the whole chart, plus geodetic -> chart -> geodetic, the corner scale
//! error, and agreement with an independent std-libm Snyder implementation.
use vz_tilekey::chart::*;
use vz_tilekey::golden::{origins, Mix};
use vz_tilekey::CHART_HALF_EXTENT_MM as H;

fn charts() -> Vec<Chart> {
    let mut v: Vec<Chart> = origins().into_iter().map(|d| Chart::new(d).unwrap()).collect();
    // Extremes the declaration allows: 89 degrees (chart crosses the pole) and the antimeridian.
    v.push(Chart::new(ChartDeclV1::new(890_000_000, 1_800_000_000, 6_371_000_000)).unwrap());
    v.push(Chart::new(ChartDeclV1::new(-890_000_000, -1_800_000_000, 6_371_000_000)).unwrap());
    v
}

fn sweep_points() -> Vec<(i64, i64)> {
    let mut p = Vec::new();
    let n = 256i64;
    for i in 0..=n { for j in 0..=n { p.push((-H + (2 * H) * i / n, -H + (2 * H) * j / n)); } }
    let mut m = Mix(0x6008);
    for _ in 0..50_000 { p.push((m.int(-H, H), m.int(-H, H))); }
    p
}

#[test]
fn g8_chart_geodetic_chart_within_1mm_whole_chart() {
    let pts = sweep_points();
    let mut report = Vec::new();
    for c in charts() {
        let mut max_f = 0f64;
        let mut int_mismatch = 0usize;
        for (x, y) in &pts {
            let (la, lo) = c.inverse_mm(*x, *y).unwrap();
            let [fx, fy] = c.forward_f64(la, lo).unwrap();
            let e = ((fx - *x as f64).powi(2) + (fy - *y as f64).powi(2)).sqrt();
            max_f = max_f.max(e);
            if c.forward_mm(la, lo).unwrap() != [*x, *y] { int_mismatch += 1; }
        }
        assert!(max_f <= 1.0, "G8 failed: {max_f} mm for {:?}", c.decl);
        assert_eq!(int_mismatch, 0, "integer round trip mismatches for {:?}", c.decl);
        report.push(format!("origin ({}, {}) R {}: points {} max |chart->geo->chart| {:.3e} mm, integer mismatches {}",
            c.decl.origin_lat_e7, c.decl.origin_lon_e7, c.decl.datum_radius_mm, pts.len(), max_f, int_mismatch));
    }
    eprintln!("G8\n{}", report.join("\n"));
}

/// Great-circle distance in mm (std libm, independent of detmath).
fn gc(r: f64, a: (f64, f64), b: (f64, f64)) -> f64 {
    let (p1, l1, p2, l2) = (a.0, a.1, b.0, b.1);
    let h = ((p2 - p1) / 2.0).sin().powi(2) + p1.cos() * p2.cos() * ((l2 - l1) / 2.0).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

#[test]
fn geodetic_chart_geodetic_within_1mm_on_the_sphere() {
    let mut worst = 0f64;
    for c in charts() {
        let mut m = Mix(0x6009);
        for _ in 0..20_000 {
            let (x, y) = (m.int(-H, H) as f64 + 0.25, m.int(-H, H) as f64 - 0.25);
            let g = c.inverse_f64(x.clamp(-H as f64, H as f64), y.clamp(-H as f64, H as f64)).unwrap();
            let [fx, fy] = c.forward_f64(g.0, g.1).unwrap();
            let g2 = c.inverse_f64(fx.clamp(-H as f64, H as f64), fy.clamp(-H as f64, H as f64)).unwrap();
            worst = worst.max(gc(c.radius_mm(), g, g2));
        }
    }
    assert!(worst <= 1.0, "{worst} mm");
    eprintln!("geo->chart->geo worst great-circle error {worst:.3e} mm");
}

#[test]
fn matches_independent_std_libm_snyder_within_a_micrometre() {
    // Snyder (1987) eq. 21-2 to 21-4 with std trig; inverse eq. 20-14, 21-15.
    let mut worst = 0f64;
    for c in charts() {
        let (p1, l0) = c.origin_rad();
        let r = c.radius_mm();
        let mut m = Mix(0x600a);
        for _ in 0..20_000 {
            let (x, y) = (m.int(-H, H) as f64, m.int(-H, H) as f64);
            let (la, lo) = c.inverse_f64(x, y).unwrap();
            let k = 2.0 / (1.0 + p1.sin() * la.sin() + p1.cos() * la.cos() * (lo - l0).cos());
            let sx = r * k * la.cos() * (lo - l0).sin();
            let sy = r * k * (p1.cos() * la.sin() - p1.sin() * la.cos() * (lo - l0).cos());
            worst = worst.max(((sx - x).powi(2) + (sy - y).powi(2)).sqrt());
        }
    }
    assert!(worst < 1e-3, "{worst} mm");
    eprintln!("detmath chart vs std-libm Snyder: worst {worst:.3e} mm");
}

#[test]
fn corner_scale_error_about_2e_4() {
    let c = Chart::new(origins()[1].clone()).unwrap();
    let k = c.scale_factor(H as f64, H as f64);
    assert!((k - 1.0 - 2.116e-4).abs() < 1e-6, "{k}");
    // Measured, not just the formula: a 1 m chart step near the corner vs the
    // great-circle distance between the two inverse points.
    let (x, y) = ((H - 2000) as f64, (H - 2000) as f64);
    let a = c.inverse_f64(x, y).unwrap();
    let b = c.inverse_f64(x + 1000.0, y).unwrap();
    let ground = gc(c.radius_mm(), a, b);
    let measured = 1000.0 / ground;
    assert!((measured - k).abs() < 1e-6, "measured {measured} formula {k}");
    eprintln!("corner scale factor formula {:.6e} measured {:.6e}", k - 1.0, measured - 1.0);
}

#[test]
fn refusals() {
    let c = Chart::new(origins()[0].clone()).unwrap();
    assert_eq!(c.inverse_mm(H + 1, 0).unwrap_err().0, "outside_chart");
    assert_eq!(c.forward_mm(e9_to_rad(10_000_000_000), 0.0).unwrap_err().0, "outside_chart");
    assert_eq!(c.forward_mm(0.0, std::f64::consts::PI).unwrap_err().0, "outside_chart");
    assert_eq!(c.forward_f64(2.0, 0.0).unwrap_err().0, "geodetic_out_of_range");
    assert_eq!(Chart::new(ChartDeclV1::new(900_000_000, 0, 6_371_000_000)).unwrap_err().0, "chart_origin_out_of_range");
    let mut d = ChartDeclV1::new(0, 0, 6_371_000_000);
    d.kind = "local_equirectangular_v1".into();
    assert_eq!(Chart::new(d).unwrap_err().0, "unsupported_chart_kind");
    assert_eq!(Chart::new(ChartDeclV1::new(0, 0, 1_000)).unwrap_err().0, "chart_radius_out_of_range");
    assert_eq!(c.inverse_mm(0, 0).unwrap(), c.origin_rad());
    assert_eq!(c.forward_mm(0.0, 0.0).unwrap(), [0, 0]);
}

#[test]
fn chart_decl_has_ground_projection_shape() {
    let d = ChartDeclV1::new(123_456_789, -45_678_901, 6_371_000_000);
    let s = serde_json::to_string(&d).unwrap();
    assert_eq!(s, r#"{"kind":"local_oblique_stereographic_v1","origin_lat_e7":123456789,"origin_lon_e7":-45678901,"datum_radius_mm":6371000000}"#);
    let back: ChartDeclV1 = serde_json::from_str(&s).unwrap();
    assert_eq!(back, d);
}
