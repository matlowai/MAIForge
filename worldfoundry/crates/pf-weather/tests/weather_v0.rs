//! WEATHER-V0 properties on a synthetic, uniform climate (no world files).
//! Run with `cargo test --release -p pf-weather` (the calibration test
//! evaluates ~10^5 states).

use pf_weather::{ClimateGrid, WeatherModel};
use std::path::PathBuf;

fn uniform(p_mm: f32, t_c: f32, ai: f32, wind_u: f32) -> ClimateGrid {
    let (nlat, nlon) = (36, 72);
    let n = nlat * nlon;
    ClimateGrid {
        nlat,
        nlon,
        temperature_mean_c: vec![t_c; n],
        temperature_half_range_k: vec![8.0; n],
        precipitation_mm_per_yr: vec![p_mm; n],
        precipitation_orographic_mm_per_yr: vec![0.0; n],
        aridity_index: vec![ai; n],
        wind_u_m_s: vec![wind_u; n],
        wind_v_m_s: vec![0.0; n],
        elevation_m: vec![100.0; n],
        climate_digest: "synthetic".into(),
        terrain_digest: "synthetic".into(),
        climate_dir: PathBuf::new(),
        terrain_dir: PathBuf::new(),
        world_seed: 7,
        axial_tilt_deg: Some(23.44),
        lapse_rate_k_per_m: Some(0.0065),
    }
}

fn bits(s: &pf_weather::WeatherState) -> Vec<u64> {
    let e = &s.extinction_per_km;
    [s.cloud_cover, s.precip_rate_mm_h, s.wind_u_m_s, s.wind_v_m_s, s.air_temperature_c, s.dewpoint_c, s.relative_humidity,
     s.visibility_km, e.rayleigh, e.haze, e.precip, e.dust, e.fog, s.radiative_cooling]
        .iter().map(|v| v.to_bits()).collect()
}

#[test]
fn deterministic_and_random_access() {
    let c = uniform(1200.0, 12.0, 1.0, 5.0);
    let a = WeatherModel::new(&c, 42, None).unwrap();
    let series: Vec<Vec<u64>> = (0..500).map(|h| bits(&a.state_at(44.0, -50.3, 1000.0 + h as f64, None))).collect();
    // A fresh model asked for single hours out of order gives the same bits.
    for &h in &[499usize, 0, 250, 17] {
        let b = WeatherModel::new(&c, 42, None).unwrap();
        assert_eq!(bits(&b.state_at(44.0, -50.3, 1000.0 + h as f64, None)), series[h], "hour {h}");
    }
    // A different seed gives a different record.
    let d = WeatherModel::new(&c, 43, None).unwrap();
    let other: Vec<Vec<u64>> = (0..500).map(|h| bits(&d.state_at(44.0, -50.3, 1000.0 + h as f64, None))).collect();
    assert_ne!(series, other);
}

#[test]
fn wet_fraction_and_mean_rate_are_calibrated() {
    // Uniform climate: every place has the same target, so pooling points and
    // hours estimates the expectation. Points spread over both hemispheres.
    let c = uniform(1500.0, 14.0, 1.2, 5.0);
    let m = WeatherModel::new(&c, 42, None).unwrap();
    let mut wet = 0usize;
    let mut sum = 0.0;
    let mut n = 0usize;
    let mut target = 0.0;
    for p in 0..24 {
        let lat = -55.0 + 110.0 * (p as f64 + 0.5) / 24.0;
        let lon = -170.0 + 337.0 * ((p * 7) % 24) as f64 / 24.0;
        for h in (0..8760).step_by(3) {
            let s = m.state_at(lat, lon, h as f64, None);
            if s.precip_rate_mm_h > 0.0 { wet += 1; }
            sum += s.precip_rate_mm_h;
            target = s.wet_hour_fraction;
            n += 1;
        }
        m.field.clear_memo();
    }
    let f = wet as f64 / n as f64;
    let annual = sum / n as f64 * 8760.0;
    eprintln!("wet fraction {f:.4} (target {target:.4}); mean annual {annual:.0} mm (target 1500); {n} samples");
    assert!((f / target - 1.0).abs() < 0.12, "wet fraction {f} vs {target}");
    assert!((annual / 1500.0 - 1.0).abs() < 0.12, "annual {annual}");
}
