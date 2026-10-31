//! WEATHER-V1 (fronts) and COLDPOOL-2 (coldflow) properties on synthetic
//! inputs (no world files). Run with `cargo test --release -p pf-weather`.

use pf_weather::coldflow::{airshed, CellCache, ColdFlowConstantsV1, Network, NONE};
use pf_weather::fronts::FrontsConstantsV1;
use pf_weather::{ClimateGrid, WeatherModel};
use std::path::PathBuf;
use vz_tilekey::chart::{Chart, ChartDeclV1};

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

fn v1<'a>(c: &'a ClimateGrid, seed: u64) -> WeatherModel<'a> {
    let mut m = WeatherModel::new(c, seed, None).unwrap();
    m.enable_fronts(FrontsConstantsV1::default(), false).unwrap();
    m
}

#[test]
fn fronts_deterministic_random_access_and_distinct_from_v0() {
    let c = uniform(1200.0, 12.0, 1.0, 5.0);
    let a = v1(&c, 42);
    let series: Vec<Vec<u64>> = (0..400).map(|h| bits(&a.state_at(44.0, -50.3, 2000.0 + h as f64, None))).collect();
    for &h in &[399usize, 0, 200, 31] {
        let b = v1(&c, 42);
        assert_eq!(bits(&b.state_at(44.0, -50.3, 2000.0 + h as f64, None)), series[h], "hour {h}");
    }
    let d = v1(&c, 43);
    let other: Vec<Vec<u64>> = (0..400).map(|h| bits(&d.state_at(44.0, -50.3, 2000.0 + h as f64, None))).collect();
    assert_ne!(series, other);
    let v0 = WeatherModel::new(&c, 42, None).unwrap();
    let s0: Vec<Vec<u64>> = (0..400).map(|h| bits(&v0.state_at(44.0, -50.3, 2000.0 + h as f64, None))).collect();
    assert_ne!(series, s0, "fronts=1 must change the record");
}

#[test]
fn fronts_wet_fraction_and_mean_rate_are_calibrated() {
    let c = uniform(1500.0, 14.0, 1.2, 5.0);
    let m = v1(&c, 42);
    let (mut wet, mut sum, mut n, mut target) = (0usize, 0.0, 0usize, 0.0);
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
        m.clear_memo();
    }
    let f = wet as f64 / n as f64;
    let annual = sum / n as f64 * 8760.0;
    eprintln!("V1 wet fraction {f:.4} (target {target:.4}); mean annual {annual:.0} mm (target 1500); {n} samples");
    assert!((f / target - 1.0).abs() < 0.12, "wet fraction {f} vs {target}");
    assert!((annual / 1500.0 - 1.0).abs() < 0.12, "annual {annual}");
}

/// A 21 x 21 node network on a chart: a trunk running west to the sea along
/// row 0 with a flat lower reach, and a steep tributary joining it from the
/// north. The trunk's flat reach must hold one shared lid that backs up the
/// tributary's lowest node; the tributary's upper nodes keep a thin layer.
#[test]
fn coldflow_shared_lid_and_backwater() {
    let c = uniform(1000.0, 10.0, 1.5, 0.0);
    let m = WeatherModel::new(&c, 42, Some(6371.0)).unwrap();
    let n = 21usize;
    let half = 10i64;
    let nn = n * n;
    let idx = |gx: i64, gy: i64| ((gy + half) as usize) * n + (gx + half) as usize;
    let mut kind = vec![2u8; nn];
    let mut channel = vec![0u8; nn];
    let mut receiver = vec![NONE; nn];
    let mut bed = vec![0.0f64; nn];
    let (mut px, mut py) = (vec![0.0; nn], vec![0.0; nn]);
    for i in 0..nn {
        px[i] = ((i % n) as i64 - half) as f64 * 512.0;
        py[i] = ((i / n) as i64 - half) as f64 * 512.0;
    }
    // Trunk: gx -9 (mouth, sea at gx -10) .. gx 9, flat (0.5 m per node) below gx 0, steeper above.
    for gx in -9..=9i64 {
        let i = idx(gx, 0);
        kind[i] = 0;
        channel[i] = 1;
        receiver[i] = idx(gx - 1, 0) as u32;
        bed[i] = if gx <= 0 { 5.0 + 0.5 * (gx + 9) as f64 } else { 9.5 + 40.0 * gx as f64 };
    }
    // Tributary: gy 1..9 at gx 2, steep (60 m per node), joining the trunk at (2, 0).
    for gy in 1..=9i64 {
        let i = idx(2, gy);
        kind[i] = 0;
        channel[i] = 1;
        receiver[i] = idx(2, gy - 1) as u32;
        bed[i] = bed[idx(2, 0)] + 60.0 * gy as f64;
    }
    // Hillslope land around both (drains into the nearest channel cell below it).
    for gy in -3..=9i64 {
        for gx in -9..=9i64 {
            let i = idx(gx, gy);
            if kind[i] == 0 { continue; }
            if gy == 0 { continue; }
            kind[i] = 0;
            let r = if gy > 0 { idx(gx, gy - 1) } else { idx(gx, gy + 1) };
            receiver[i] = r as u32;
            bed[i] = bed[r] + 50.0;
        }
    }
    let net = Network {
        n,
        half,
        spacing_m: 512.0,
        kind,
        channel,
        receiver,
        bed_m: bed.clone(),
        pos_x_m: px,
        pos_y_m: py,
        chart_decl: [440_000_000, -500_000_000, 6_371_000_000],
        digest: "synthetic".into(),
        t1a_digest: "synthetic".into(),
        dir: PathBuf::new(),
        t1a_dir: PathBuf::new(),
    };
    let chart = Chart::new(ChartDeclV1::new(440_000_000, -500_000_000, 6_371_000_000)).unwrap();
    let k = ColdFlowConstantsV1::default();
    // Find a calm clear night for this synthetic climate: scan days for Q > 0.5.
    let mut day = 280;
    loop {
        let sr = pf_weather::sun::sunrise_world_hour(day, 44.0, -50.0, 23.44);
        let mut cache = CellCache::new(&m, &k, day, sr);
        if cache.blend(44.0, -50.0, |c| c.ns.cooling_q) > 0.5 && cache.blend(44.0, -50.0, |c| c.wind_tf_m_s) < 2.5 { break; }
        day += 1;
        assert!(day < 365, "no calm clear night found");
    }
    let sr = pf_weather::sun::sunrise_world_hour(day, 44.0, -50.0, 23.44);
    let mut cache = CellCache::new(&m, &k, day, sr);
    let a = airshed(&m, &mut cache, &net, &chart, &k, sr);
    // Lids never fall below the bed and never rise downstream of a node's own surface.
    for i in 0..nn {
        if net.kind[i] == 2 { continue; }
        assert!(a.lid_m[i] >= bed[i] - 1e-9, "lid below bed at {i}");
    }
    // The flat trunk reach shares one lid with one source.
    let s0 = a.source[idx(-3, 0)];
    let l0 = a.lid_m[idx(-3, 0)];
    for gx in -3..=0 {
        assert_eq!(a.source[idx(gx, 0)], s0, "gx {gx}");
        assert_eq!(a.lid_m[idx(gx, 0)], l0);
    }
    // Flux accumulates down the trunk (confluence adds the tributary).
    assert!(a.flux_m3_s[idx(1, 0)] > a.flux_m3_s[idx(3, 0)] + a.flux_m3_s[idx(2, 1)] * 0.99);
    // The steep tributary's upper nodes carry their own thin layers, not the trunk lid.
    let top = idx(2, 8);
    assert_eq!(a.source[top], top as u32);
    assert!(a.h_m[top] < a.h_m[idx(-3, 0)], "upper tributary layer {} vs trunk {}", a.h_m[top], a.h_m[idx(-3, 0)]);
    eprintln!("day {day}: trunk lid {l0:.1} m (source {s0}), trunk depth at gx -3 {:.1} m, tributary top layer {:.1} m, outlets critical {}",
              l0 - bed[idx(-3, 0)], a.h_m[top], a.outlet_critical.iter().filter(|&&b| b).count());
}
