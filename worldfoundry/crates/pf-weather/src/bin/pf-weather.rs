//! pf-weather CLI (lane WEATHER-COLDPOOL, 2026-10-21).
//!
//! usage (key=value arguments; outputs are never overwritten):
//!   pf-weather constants
//!   pf-weather series  climate=DIR terrain=DIR (cell=ID | lat=DEG lon=DEG) [elev=M] [year=0] [seed=N] out=PREFIX
//!   pf-weather vmap    climate=DIR terrain=DIR hour=H (window=DIR [half_km=131] [px_km=1] | lat0=D lon0=D [span_deg=12] [px_deg=0.05]) [seed=N] out=PREFIX
//!   pf-weather nights  climate=DIR terrain=DIR cell=ID [year=0] [seed=N] out=FILE.tsv
//!   pf-weather coldpool climate=DIR terrain=DIR windows=DIR day=D [times=-1.5,0,2,4] [only=w_0_0,..] [seed=N] out=DIR
//!       [coldflow=1 net=RIVERNET_DIR t1a=T1A_DIR]   (lane WEATHER-V1: COLDPOOL-2 cold-air drainage)
//!
//! Every weather command also takes `fronts=1` (lane WEATHER-V1: frontal band
//! families and convective clusters) and, with it, `oro=1` (orographic
//! modulation of frontal rain). With these keys absent the outputs are the
//! WEATHER-V0 / COLDPOOL-1 outputs byte for byte.
//!
//! Hours are world hours (hour 0 = 00:00 at longitude 0, day 0, year 0, 365-day
//! calendar). `day` is a local calendar day index of the morning in question.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use pf_weather::coldflow::{airshed, coldflow_window, CellCache, ColdFlowConstantsV1, Network, MODEL_ID as COLDFLOW_MODEL_ID, SCHEMA as COLDFLOW_SCHEMA};
use pf_weather::coldpool::{coldpool_window, Window, SCHEMA as COLDPOOL_SCHEMA};
use pf_weather::fronts::{FrontsConstantsV1, FRONTS_MODEL_ID};
use pf_weather::state::SCHEMA as STATE_SCHEMA;
use pf_weather::{ClimateGrid, WeatherModel, WeatherState, LAW_ID};
use serde_json::json;
use vz_tilekey::chart::Chart;

fn args() -> (String, BTreeMap<String, String>) {
    let mut it = std::env::args().skip(1);
    let cmd = it.next().unwrap_or_default();
    let mut m = BTreeMap::new();
    for a in it {
        if let Some((k, v)) = a.split_once('=') {
            m.insert(k.to_string(), v.to_string());
        } else {
            die(&format!("argument {a} is not key=value"));
        }
    }
    (cmd, m)
}

fn die(msg: &str) -> ! {
    eprintln!("pf-weather: {msg}");
    std::process::exit(2)
}

fn get<'a>(m: &'a BTreeMap<String, String>, k: &str) -> &'a str {
    m.get(k).map(|s| s.as_str()).unwrap_or_else(|| die(&format!("missing {k}=")))
}
fn getf(m: &BTreeMap<String, String>, k: &str, d: Option<f64>) -> f64 {
    match m.get(k) {
        Some(v) => v.parse().unwrap_or_else(|_| die(&format!("{k} not a number"))),
        None => d.unwrap_or_else(|| die(&format!("missing {k}="))),
    }
}

fn fresh(p: &Path) {
    if p.exists() {
        die(&format!("{} exists; outputs are never overwritten", p.display()));
    }
}

fn write_new(p: &Path, bytes: &[u8]) {
    fresh(p);
    if let Some(d) = p.parent() { std::fs::create_dir_all(d).unwrap_or_else(|e| die(&e.to_string())); }
    std::fs::write(p, bytes).unwrap_or_else(|e| die(&format!("{}: {e}", p.display())));
}

fn digest(bytes: &[u8]) -> String { format!("blake3:{}", blake3::hash(bytes).to_hex()) }

fn f32_bytes(v: &[f32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }

fn load(m: &BTreeMap<String, String>) -> ClimateGrid {
    ClimateGrid::load(Path::new(get(m, "climate")), Path::new(get(m, "terrain"))).unwrap_or_else(|e| die(&e))
}

fn seed(m: &BTreeMap<String, String>, c: &ClimateGrid) -> u64 {
    m.get("seed").map(|s| s.parse().unwrap_or_else(|_| die("seed"))).unwrap_or(c.world_seed)
}

fn flag(m: &BTreeMap<String, String>, k: &str) -> bool {
    match m.get(k).map(|s| s.as_str()) {
        None | Some("0") => false,
        Some("1") => true,
        Some(v) => die(&format!("{k}={v}: use 0 or 1")),
    }
}

/// Apply the opt-in WEATHER-V1 keys to a model.
fn apply_v1(model: &mut WeatherModel, m: &BTreeMap<String, String>) {
    let oro = flag(m, "oro");
    if flag(m, "fronts") {
        model.enable_fronts(FrontsConstantsV1::default(), oro).unwrap_or_else(|e| die(&e));
    } else if oro {
        die("oro=1 needs fronts=1");
    }
}

fn inputs_json(c: &ClimateGrid, model: &WeatherModel, seed: u64) -> serde_json::Value {
    let mut v = inputs_json_v0(c, model, seed);
    if let Some(f) = &model.fronts {
        v["model"] = json!({"id": FRONTS_MODEL_ID, "fronts_constants_digest": f.k.digest(), "oro": model.oro,
                            "fronts_constants": f.k.table(),
                            "note": "WEATHER-V1 (lane WEATHER-V1): frontal band families + convective clusters replace the V0 discs; V0 law for climatology, calibration, temperature, humidity, wind and extinction"});
    }
    v
}

fn inputs_json_v0(c: &ClimateGrid, model: &WeatherModel, seed: u64) -> serde_json::Value {
    json!({
        "law": LAW_ID,
        "constants_digest": model.k.digest(),
        "world_seed": seed,
        "keyed_rng": vz_keyed_rng::KEYED_RNG_V1,
        "climate_digest": c.climate_digest,
        "terrain_digest": c.terrain_digest,
        "climate_dir": c.climate_dir.display().to_string(),
        "terrain_dir": c.terrain_dir.display().to_string(),
        "grid": format!("surface.latlon.{}x{}", c.nlat, c.nlon),
        "axial_tilt_deg": model.k.axial_tilt_deg,
        "lapse_rate_k_per_m": model.k.lapse_rate_k_per_m,
        "planet_radius_km": model.k.planet_radius_km,
    })
}

const COLS: [&str; 23] = [
    "t_h", "cloud_cover", "precip_type", "precip_rate_mm_h", "convective", "wind_u_m_s", "wind_v_m_s", "wind_speed_m_s",
    "air_temperature_c", "dewpoint_c", "relative_humidity", "visibility_km", "visibility_open_ground_km", "beta_rayleigh", "beta_haze", "beta_precip",
    "beta_dust", "beta_fog", "radiative_cooling", "flat_fog", "synoptic_activity", "sun_elevation_deg", "solar_hour",
];

fn row(s: &WeatherState) -> [f32; 23] {
    let e = &s.extinction_per_km;
    [
        s.t_h as f32, s.cloud_cover as f32, s.precip_type.code() as f32, s.precip_rate_mm_h as f32, s.convective as u8 as f32,
        s.wind_u_m_s as f32, s.wind_v_m_s as f32, s.wind_speed_m_s as f32, s.air_temperature_c as f32, s.dewpoint_c as f32,
        s.relative_humidity as f32, s.visibility_km as f32, s.visibility_open_ground_km as f32, e.rayleigh as f32, e.haze as f32, e.precip as f32, e.dust as f32,
        e.fog as f32, s.radiative_cooling as f32, s.flat_fog as f32, s.synoptic_activity as f32, s.sun_elevation_deg as f32,
        s.solar_hour as f32,
    ]
}

fn pct(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() { return f64::NAN; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = ((v.len() - 1) as f64 * p).round() as usize;
    v[i]
}

fn cmd_series(m: &BTreeMap<String, String>) {
    let c = load(m);
    let sd = seed(m, &c);
    let mut model = WeatherModel::new(&c, sd, None).unwrap_or_else(|e| die(&e));
    apply_v1(&mut model, m);
    let (lat, lon, cell) = if let Some(id) = m.get("cell") {
        let id: u32 = id.parse().unwrap_or_else(|_| die("cell"));
        let (la, lo) = c.cell_centre(id);
        (la, lo, Some(id))
    } else {
        (getf(m, "lat", None), getf(m, "lon", None), None)
    };
    let elev = m.get("elev").map(|v| v.parse::<f64>().unwrap_or_else(|_| die("elev")));
    let year = getf(m, "year", Some(0.0)) as u64;
    let out = PathBuf::from(get(m, "out"));
    let h0 = year as f64 * 8760.0;
    let mut data: Vec<f32> = Vec::with_capacity(8760 * COLS.len());
    let mut states = Vec::with_capacity(8760);
    for h in 0..8760u64 {
        // The memo only caches keyed draws; clearing it changes no output.
        if model.fronts.is_some() && h % 720 == 0 { model.clear_memo(); }
        let s = model.state_at(lat, lon, h0 + h as f64, elev);
        data.extend_from_slice(&row(&s));
        states.push(s);
    }
    // Summary (every number names its metric; one year, one seed).
    let n = states.len() as f64;
    let precip: f64 = states.iter().map(|s| s.precip_rate_mm_h).sum();
    let wet = states.iter().filter(|s| s.precip_rate_mm_h > 0.0).count() as f64;
    let snow: f64 = states.iter().filter(|s| s.precip_type == pf_weather::PrecipType::Snow).map(|s| s.precip_rate_mm_h).sum();
    let mut vis: Vec<f64> = states.iter().map(|s| s.visibility_km).collect();
    let fog_h = states.iter().filter(|s| s.visibility_km < 1.0).count();
    let fog_open_h = states.iter().filter(|s| s.visibility_open_ground_km < 1.0).count();
    let dust_h = states.iter().filter(|s| s.extinction_per_km.dust > 0.0).count();
    let calm_clear_night_h = states.iter().filter(|s| s.radiative_cooling > 0.5).count();
    let mut temps: Vec<f64> = states.iter().map(|s| s.air_temperature_c).collect();
    let mean_cloud = states.iter().map(|s| s.cloud_cover).sum::<f64>() / n;
    let mean_wind = states.iter().map(|s| s.wind_speed_m_s).sum::<f64>() / n;
    let mean_rh = states.iter().map(|s| s.relative_humidity).sum::<f64>() / n;
    let cs = c.sample(lat, lon);
    let summary = json!({
        "hours": states.len(),
        "precip_realized_mm": precip,
        "precip_climatology_mm_per_yr": cs.precipitation_mm_per_yr,
        "precip_ratio_realized_over_climatology": if cs.precipitation_mm_per_yr > 0.0 { precip / cs.precipitation_mm_per_yr } else { f64::NAN },
        "wet_hour_fraction_realized": wet / n,
        "wet_hour_fraction_target": states[0].wet_hour_fraction,
        "snow_share_of_precip": if precip > 0.0 { snow / precip } else { 0.0 },
        "visibility_km_p05_p50_p95": [pct(&mut vis.clone(), 0.05), pct(&mut vis.clone(), 0.5), pct(&mut vis, 0.95)],
        "hours_visibility_below_1km": fog_h,
        "hours_open_ground_visibility_below_1km": fog_open_h,
        "hours_with_dust": dust_h,
        "hours_radiative_cooling_above_0_5": calm_clear_night_h,
        "air_temperature_c_p05_p50_p95": [pct(&mut temps.clone(), 0.05), pct(&mut temps.clone(), 0.5), pct(&mut temps, 0.95)],
        "temperature_mean_c_climatology_at_t0_elevation": cs.temperature_mean_c,
        "mean_cloud_cover": mean_cloud,
        "mean_wind_speed_m_s": mean_wind,
        "mean_relative_humidity": mean_rh,
        "memo_entries": model.field.memo_len(),
    });
    let bin = f32_bytes(&data);
    let header = json!({
        "schema": "pf.weather-series/v0",
        "state_schema": STATE_SCHEMA,
        "label": "derived weather record (pure function of canonical inputs); presentation and perception read it; not canon by itself",
        "inputs": inputs_json(&c, &model, sd),
        "point": {"lat_deg": lat, "lon_deg": lon, "cell": cell, "t0_cell": c.cell_of(lat, lon), "elevation_m": states[0].elevation_m, "climate": {
            "temperature_mean_c": cs.temperature_mean_c, "temperature_half_range_k": cs.temperature_half_range_k,
            "precipitation_mm_per_yr": cs.precipitation_mm_per_yr, "precipitation_orographic_mm_per_yr": cs.precipitation_orographic_mm_per_yr,
            "aridity_index": cs.aridity_index, "wind_u_m_s": cs.wind_u_m_s, "wind_v_m_s": cs.wind_v_m_s, "t0_elevation_m": cs.t0_elevation_m}},
        "year": year,
        "first_hour": h0,
        "columns": COLS,
        "precip_type_codes": {"0": "none", "1": "drizzle", "2": "rain", "3": "mixed", "4": "snow"},
        "rows": states.len(),
        "file": format!("{}.f32", out.file_name().unwrap().to_string_lossy()),
        "file_digest": digest(&bin),
        "summary": summary,
    });
    write_new(&PathBuf::from(format!("{}.f32", out.display())), &bin);
    write_new(&PathBuf::from(format!("{}.json", out.display())), serde_json::to_string_pretty(&header).unwrap().as_bytes());
    println!("{}", serde_json::to_string_pretty(&header["summary"]).unwrap());
}

const MAP_LAYERS: [&str; 8] = ["visibility_km", "visibility_open_ground_km", "precip_rate_mm_h", "precip_type", "cloud_cover", "air_temperature_c", "wind_speed_m_s", "t0_elevation_m"];

fn cmd_vmap(m: &BTreeMap<String, String>) {
    let c = load(m);
    let sd = seed(m, &c);
    let hour = getf(m, "hour", None);
    let out = PathBuf::from(get(m, "out"));
    // Two placements of the same planet-keyed samples: a chart grid (window=DIR,
    // chart km, the chart is +-131.072 km) or a lat/lon grid (lat0, lon0,
    // span_deg, px_deg) for regional context beyond one chart.
    let window = m.get("window").map(|w| Window::load(Path::new(w)).unwrap_or_else(|e| die(&e)));
    let radius = window.as_ref().map(|w| w.chart.radius_mm() / 1e6);
    let mut model = WeatherModel::new(&c, sd, radius).unwrap_or_else(|e| die(&e));
    apply_v1(&mut model, m);
    let (nx, grid_json, pts): (usize, serde_json::Value, Vec<(f64, f64)>) = if let Some(w) = &window {
        let chart: &Chart = &w.chart;
        let half = getf(m, "half_km", Some(131.0));
        let px = getf(m, "px_km", Some(1.0));
        let nx = (2.0 * half / px).round() as usize;
        let mut pts = Vec::with_capacity(nx * nx);
        for j in 0..nx {
            for i in 0..nx {
                let x = (-half + (i as f64 + 0.5) * px) * 1e6;
                let y = (-half + (j as f64 + 0.5) * px) * 1e6;
                let (la, lo) = chart.inverse_f64(x, y).unwrap_or_else(|e| die(&format!("chart: {e:?}")));
                pts.push((la / pf_weather::sun::DEG, lo / pf_weather::sun::DEG));
            }
        }
        let (la0, lo0) = chart.origin_rad();
        (nx, json!({"kind": "chart", "nx": nx, "ny": nx, "px_km": px, "x0_km": -half, "y0_km": -half,
                    "origin_lat_deg": la0 / pf_weather::sun::DEG, "origin_lon_deg": lo0 / pf_weather::sun::DEG,
                    "from_window": get(m, "window"), "order": "layer, row (y north), column (x east)"}), pts)
    } else {
        let (lat0, lon0) = (getf(m, "lat0", None), getf(m, "lon0", None));
        let span = getf(m, "span_deg", Some(12.0));
        let px = getf(m, "px_deg", Some(0.05));
        let nx = (span / px).round() as usize;
        let mut pts = Vec::with_capacity(nx * nx);
        for j in 0..nx {
            for i in 0..nx {
                pts.push((lat0 - 0.5 * span + (j as f64 + 0.5) * px, lon0 - 0.5 * span + (i as f64 + 0.5) * px));
            }
        }
        (nx, json!({"kind": "latlon", "nx": nx, "ny": nx, "px_deg": px, "lat_min": lat0 - 0.5 * span, "lon_min": lon0 - 0.5 * span,
                    "order": "layer, row (lat north), column (lon east)"}), pts)
    };
    let mut layers = vec![0f32; MAP_LAYERS.len() * nx * nx];
    for (k, &(la, lo)) in pts.iter().enumerate() {
        let s = model.state_at(la, lo, hour, None);
        let vals = [s.visibility_km, s.visibility_open_ground_km, s.precip_rate_mm_h, s.precip_type.code() as f64, s.cloud_cover, s.air_temperature_c, s.wind_speed_m_s, s.elevation_m];
        for (l, v) in vals.iter().enumerate() {
            layers[l * nx * nx + k] = *v as f32;
        }
    }
    let bin = f32_bytes(&layers);
    let header = json!({
        "schema": "pf.weather-map/v0",
        "label": "derived weather record sampled on a grid (planet-keyed; the grid only places the samples)",
        "inputs": inputs_json(&c, &model, sd),
        "hour": hour,
        "grid": grid_json,
        "layers": MAP_LAYERS,
        "file": format!("{}.f32", out.file_name().unwrap().to_string_lossy()),
        "file_digest": digest(&bin),
        "elevation_note": "samples at the T0 surface elevation (bilinear), not T1b heights",
    });
    write_new(&PathBuf::from(format!("{}.f32", out.display())), &bin);
    write_new(&PathBuf::from(format!("{}.json", out.display())), serde_json::to_string_pretty(&header).unwrap().as_bytes());
    let mut v: Vec<f64> = (0..nx * nx).map(|k| layers[k] as f64).collect();
    println!("vmap hour {hour}: V km p05 {:.2} p50 {:.2} p95 {:.2}", pct(&mut v.clone(), 0.05), pct(&mut v.clone(), 0.5), pct(&mut v, 0.95));
}

fn cmd_nights(m: &BTreeMap<String, String>) {
    let c = load(m);
    let sd = seed(m, &c);
    let mut model = WeatherModel::new(&c, sd, None).unwrap_or_else(|e| die(&e));
    apply_v1(&mut model, m);
    let cell: u32 = get(m, "cell").parse().unwrap_or_else(|_| die("cell"));
    let year = getf(m, "year", Some(0.0)) as i64;
    let out = PathBuf::from(get(m, "out"));
    let mut s = String::from("local_day\tdoy\tcooling_q\train_hours\tmean_wind\tmean_cloud\tt_sunrise_c\ttd_sunrise_c\trh_sunrise\tsunrise_h\tburn_m_h\n");
    for d in 0..365i64 {
        if model.fronts.is_some() && d % 30 == 0 { model.clear_memo(); }
        let day = year * 365 + d;
        let ns = pf_weather::coldpool::night_summary(&model, cell, day, f64::INFINITY);
        s.push_str(&format!(
            "{}\t{}\t{:.4}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.4}\t{:.4}\t{:.2}\n",
            day, d, ns.cooling_q, ns.night_hours_raining, ns.mean_wind_m_s, ns.mean_cloud, ns.t_ref_c, ns.td_ref_c, ns.rh_ref, ns.sunrise_h, ns.burn_rate_m_h
        ));
    }
    write_new(&out, s.as_bytes());
    println!("wrote {}", out.display());
}

fn cmd_coldpool(m: &BTreeMap<String, String>) {
    let c = load(m);
    let sd = seed(m, &c);
    let dir = PathBuf::from(get(m, "windows"));
    let day = getf(m, "day", None) as i64;
    let times: Vec<f64> = m.get("times").map(|s| s.as_str()).unwrap_or("-1.5,0,2,4").split(',').map(|v| v.parse().unwrap_or_else(|_| die("times"))).collect();
    let out = PathBuf::from(get(m, "out"));
    let mut names: Vec<String> = match m.get("only") {
        Some(o) => o.split(',').map(|s| s.to_string()).collect(),
        None => std::fs::read_dir(&dir).unwrap_or_else(|e| die(&e.to_string()))
            .filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().to_string())
            .filter(|n| n.starts_with("w_")).collect(),
    };
    names.sort();
    if flag(m, "coldflow") {
        return cmd_coldflow(m, &c, sd, &dir, day, &times, &out, &names);
    }
    let mut summary = Vec::new();
    for name in &names {
        let w = Window::load(&dir.join(name)).unwrap_or_else(|e| die(&e));
        let mut model = WeatherModel::new(&c, sd, Some(w.chart.radius_mm() / 1e6)).unwrap_or_else(|e| die(&e));
        apply_v1(&mut model, m);
        // Sunrise at the chart origin (shared by every window of a chart)
        // defines the time offsets, so all windows are evaluated at the same
        // world hour and agree on shared posts.
        let (la, lo) = w.chart.origin_rad();
        let (la, lo) = (la / pf_weather::sun::DEG, lo / pf_weather::sun::DEG);
        let sr = pf_weather::sun::sunrise_world_hour(day, la, lo, model.k.axial_tilt_deg);
        for &dt in &times {
            let t = sr + dt;
            let tag = if dt == 0.0 { "sunrise".to_string() } else if dt < 0.0 { format!("sr_minus_{}h", -dt) } else { format!("sr_plus_{}h", dt) };
            let p = coldpool_window(&model, &w, day, t);
            let od = out.join(&tag).join(name);
            let fb = f32_bytes(&p.fog_top_m);
            let pb = f32_bytes(&p.pool_top_m);
            let header = json!({
                "schema": COLDPOOL_SCHEMA,
                "label": "derived presentation/perception layer (dawn cold-air pool and radiation-fog top); pure function of the window heights, planet location and the WEATHER-V0 record; never canon by itself",
                "inputs": inputs_json(&c, &model, sd),
                "window": {"dir": dir.join(name).display().to_string(), "digest": w.digest, "key": w.window_key, "rx": w.rx, "ry": w.ry, "n": w.n, "halo_posts": w.halo, "gx0": w.gx0, "gy0": w.gy0,
                           "order": "post = (gy - gy0) * n + (gx - gx0); x east, y north; chart mm = g * 32000 (same lattice as the T1b window)"},
                "local_day": day,
                "world_hour": t,
                "hours_from_sunrise": dt,
                "sunrise_world_hour_at_chart_origin": sr,
                "cells": p.cells,
                "fields": [
                    {"name": "fog_top_m", "dtype": "f32", "file": "fog_top_m.f32", "units": "m", "note": "fog-top height above the datum where the post is inside fog, NaN elsewhere", "digest": digest(&fb)},
                    {"name": "pool_top_m", "dtype": "f32", "file": "pool_top_m.f32", "units": "m", "note": "cold-air pool top (inversion lid) where a pool is published, NaN elsewhere", "digest": digest(&pb)},
                    {"name": "mask", "dtype": "u8", "file": "mask.u8", "units": "code", "note": "0 none, 1 fog, 2 cold pool without fog, 3 marine (no radiation fog)", "digest": digest(&p.mask)},
                ],
                "core_fog_fraction": p.fog_core_fraction,
                "core_pool_fraction": p.pool_core_fraction,
                "core_fog_top_min_m": if p.fog_top_core_min_m.is_finite() { json!(p.fog_top_core_min_m) } else { json!(null) },
                "core_fog_top_max_m": if p.fog_top_core_max_m.is_finite() { json!(p.fog_top_core_max_m) } else { json!(null) },
            });
            write_new(&od.join("fog_top_m.f32"), &fb);
            write_new(&od.join("pool_top_m.f32"), &pb);
            write_new(&od.join("mask.u8"), &p.mask);
            write_new(&od.join("coldpool.json"), serde_json::to_string_pretty(&header).unwrap().as_bytes());
            println!("{tag} {name}: fog {:.3} pool {:.3} of core; fog top {:.0}..{:.0} m; Q {:?}", p.fog_core_fraction, p.pool_core_fraction,
                     p.fog_top_core_min_m, p.fog_top_core_max_m, p.cells.iter().map(|s| (s.cell, (s.cooling_q * 1000.0).round() / 1000.0)).collect::<Vec<_>>());
            summary.push(json!({"tag": tag, "window": name, "fog": p.fog_core_fraction, "pool": p.pool_core_fraction}));
        }
    }
    write_new(&out.join("summary.json"), serde_json::to_string_pretty(&summary).unwrap().as_bytes());
}

fn time_tag(dt: f64) -> String {
    if dt == 0.0 { "sunrise".to_string() } else if dt < 0.0 { format!("sr_minus_{}h", -dt) } else { format!("sr_plus_{}h", dt) }
}

/// COLDPOOL-2 (coldflow=1): chart-wide airshed per time, then every window.
#[allow(clippy::too_many_arguments)]
fn cmd_coldflow(m: &BTreeMap<String, String>, c: &ClimateGrid, sd: u64, dir: &Path, day: i64, times: &[f64], out: &Path, names: &[String]) {
    let net = Network::load(Path::new(get(m, "net")), Path::new(get(m, "t1a"))).unwrap_or_else(|e| die(&e));
    let chart = net.chart().unwrap_or_else(|e| die(&e));
    let mut model = WeatherModel::new(c, sd, Some(chart.radius_mm() / 1e6)).unwrap_or_else(|e| die(&e));
    apply_v1(&mut model, m);
    let kc = ColdFlowConstantsV1::default();
    let (la, lo) = chart.origin_rad();
    let (la, lo) = (la / pf_weather::sun::DEG, lo / pf_weather::sun::DEG);
    let sr = pf_weather::sun::sunrise_world_hour(day, la, lo, model.k.axial_tilt_deg);
    let windows: Vec<Window> = names.iter().map(|n| Window::load(&dir.join(n)).unwrap_or_else(|e| die(&e))).collect();
    for (w, name) in windows.iter().zip(names) {
        let d = [w.chart.origin_rad().0, w.chart.origin_rad().1, w.chart.radius_mm()];
        let e = [chart.origin_rad().0, chart.origin_rad().1, chart.radius_mm()];
        if d != e { die(&format!("window {name} is not on the network's chart")); }
    }
    let net_json = json!({"rivernet_dir": net.dir.display().to_string(), "rivernet_digest": net.digest,
                          "t1a_dir": net.t1a_dir.display().to_string(), "t1a_digest": net.t1a_digest,
                          "lattice": format!("{} x {} nodes, {} m", net.n, net.n, net.spacing_m)});
    let model_json = json!({"id": COLDFLOW_MODEL_ID, "constants_digest": kc.digest(), "constants": kc.table()});
    let mut summary = Vec::new();
    for &dt in times {
        let t = sr + dt;
        let tag = time_tag(dt);
        let mut cache = CellCache::new(&model, &kc, day, t);
        let a = airshed(&model, &mut cache, &net, &chart, &kc, t);
        // Chart-wide airshed layer (lod-8 node lattice).
        let cd = out.join(&tag).join("chart");
        let f32v = |v: &[f64]| f32_bytes(&v.iter().map(|&x| x as f32).collect::<Vec<_>>());
        let files: Vec<(&str, Vec<u8>, &str)> = vec![
            ("lid_m.f32", f32v(&a.lid_m), "cold-air lid (m) at the node at this hour (NaN on sea nodes)"),
            ("fog_top_m.f32", f32v(&a.fog_top_m), "fog top (m) of the node's pool (NaN: no fog)"),
            ("h_m.f32", f32v(&a.h_m), "own drainage depth (m) at this hour"),
            ("flux_m3_s.f32", f32v(&a.flux_m3_s), "accumulated cold-air volume flux (m3/s)"),
            ("q.f32", f32v(&a.q), "night cooling factor Q"),
            ("speed_m_s.f32", f32v(&a.speed_m_s), "own steady drainage speed (m/s)"),
            ("slope.f32", f32v(&a.slope), "bed slope to the receiver used in the flow law (floored)"),
            ("pooled.u8", a.pooled.iter().map(|&b| b as u8).collect(), "1 where the node holds stagnant pooled air (fog possible)"),
            ("source.u32", a.source.iter().flat_map(|x| x.to_le_bytes()).collect(), "pool source node (airshed id); 4294967295 none"),
        ];
        let mut fl = Vec::new();
        for (fname, bytes, note) in &files {
            write_new(&cd.join(fname), bytes);
            fl.push(json!({"file": fname, "note": note, "digest": digest(bytes)}));
        }
        let crit = a.outlet_critical.iter().filter(|&&b| b).count();
        let chead = json!({"schema": "pf.coldflow-airshed/v1", "label": "derived presentation/perception layer; never canon by itself",
            "inputs": inputs_json(c, &model, sd), "model": model_json, "network": net_json, "local_day": day, "world_hour": t, "hours_from_sunrise": dt,
            "order": "node = (gy + half) * n + (gx + half), the rivernet lattice", "fields": fl,
            "land_nodes": a.land_nodes, "pooled_nodes": a.pooled_nodes, "fog_nodes": a.fog_nodes, "outlets_at_critical_depth": crit,
            "cells": cache.cells.values().collect::<Vec<_>>()});
        write_new(&cd.join("airshed.json"), serde_json::to_string_pretty(&chead).unwrap().as_bytes());
        println!("{tag} chart: land nodes {} pooled {} fog {} critical outlets {}", a.land_nodes, a.pooled_nodes, a.fog_nodes, crit);
        for (w, name) in windows.iter().zip(names) {
            let p = coldflow_window(&model, &mut cache, &net, &a, &kc, w, t);
            let od = out.join(&tag).join(name);
            let fb = f32_bytes(&p.fog_top_m);
            let pb = f32_bytes(&p.pool_top_m);
            let ib: Vec<u8> = p.pool_id.iter().flat_map(|x| x.to_le_bytes()).collect();
            let header = json!({
                "schema": COLDFLOW_SCHEMA,
                "label": "derived presentation/perception layer (dawn cold-air drainage, pools per airshed, radiation and steam fog); pure function of the window heights, the chart-wide rivernet, planet location and the WEATHER record; never canon by itself",
                "inputs": inputs_json(c, &model, sd),
                "model": model_json,
                "network": net_json,
                "window": {"dir": dir.join(name).display().to_string(), "digest": w.digest, "key": w.window_key, "rx": w.rx, "ry": w.ry, "n": w.n, "halo_posts": w.halo, "gx0": w.gx0, "gy0": w.gy0,
                           "order": "post = (gy - gy0) * n + (gx - gx0); x east, y north; chart mm = g * 32000 (same lattice as the T1b window)"},
                "local_day": day,
                "world_hour": t,
                "hours_from_sunrise": dt,
                "sunrise_world_hour_at_chart_origin": sr,
                "fields": [
                    {"name": "fog_top_m", "dtype": "f32", "file": "fog_top_m.f32", "units": "m", "note": "fog top above the datum where the post is in radiation fog (mask 1) or steam fog (mask 5), NaN elsewhere", "digest": digest(&fb)},
                    {"name": "pool_top_m", "dtype": "f32", "file": "pool_top_m.f32", "units": "m", "note": "cold-air lid where the post is below it (masks 1, 2, 4), NaN elsewhere", "digest": digest(&pb)},
                    {"name": "mask", "dtype": "u8", "file": "mask.u8", "units": "code", "note": "0 none, 1 radiation fog in a pool, 2 cold pool without fog, 3 water (no fog), 4 cold drainage layer (moving, thinner than pool_depth_min, no fog), 5 steam fog over enclosed water", "digest": digest(&p.mask)},
                    {"name": "pool_id", "dtype": "u32", "file": "pool_id.u32", "units": "node", "note": "airshed id = the rivernet node whose bed + depth set the lid (4294967295 none)", "digest": digest(&ib)},
                ],
                "core_fog_fraction": p.fog_core_fraction,
                "core_pool_fraction": p.pool_core_fraction,
                "core_drainage_fraction": p.drainage_core_fraction,
                "core_steam_fraction": p.steam_core_fraction,
                "core_unmatched_fraction": p.unmatched_core_fraction,
                "core_pool_ids": p.pool_ids_core,
                "core_fog_top_min_m": if p.fog_top_core_min_m.is_finite() { json!(p.fog_top_core_min_m) } else { json!(null) },
                "core_fog_top_max_m": if p.fog_top_core_max_m.is_finite() { json!(p.fog_top_core_max_m) } else { json!(null) },
            });
            write_new(&od.join("fog_top_m.f32"), &fb);
            write_new(&od.join("pool_top_m.f32"), &pb);
            write_new(&od.join("mask.u8"), &p.mask);
            write_new(&od.join("pool_id.u32"), &ib);
            write_new(&od.join("coldpool.json"), serde_json::to_string_pretty(&header).unwrap().as_bytes());
            println!("{tag} {name}: fog {:.3} pool {:.3} drainage {:.3} steam {:.3} unmatched {:.3} of core; fog top {:.0}..{:.0} m; airsheds {}",
                     p.fog_core_fraction, p.pool_core_fraction, p.drainage_core_fraction, p.steam_core_fraction, p.unmatched_core_fraction,
                     p.fog_top_core_min_m, p.fog_top_core_max_m, p.pool_ids_core);
            summary.push(json!({"tag": tag, "window": name, "fog": p.fog_core_fraction, "pool": p.pool_core_fraction, "drainage": p.drainage_core_fraction, "steam": p.steam_core_fraction}));
        }
    }
    write_new(&out.join("summary.json"), serde_json::to_string_pretty(&summary).unwrap().as_bytes());
}

fn main() {
    let (cmd, m) = args();
    match cmd.as_str() {
        "constants" => {
            let k = pf_weather::WeatherConstantsV0::default();
            println!("{}", serde_json::to_string_pretty(&json!({"law": LAW_ID, "digest": k.digest(), "constants": k.table()})).unwrap());
        }
        "constants-v1" => {
            let f = FrontsConstantsV1::default();
            let c = ColdFlowConstantsV1::default();
            println!("{}", serde_json::to_string_pretty(&json!({
                "fronts": {"id": FRONTS_MODEL_ID, "digest": f.digest(), "constants": f.table(),
                           "derived": {"mu_front_per_cube_slot": f.mu_front(6371.0), "mu_cluster_per_cube_slot": f.mu_cluster(6371.0),
                                       "front_reach_max_km": f.front_reach_max_km(), "cluster_reach_max_km": f.cluster_reach_max_km(),
                                       "front_cover_per_system_km2_h": f.front_cover_per_system(), "cluster_cover_per_system_km2_h": f.cluster_cover_per_system()}},
                "coldflow": {"id": COLDFLOW_MODEL_ID, "digest": c.digest(), "constants": c.table()}})).unwrap());
        }
        "series" => cmd_series(&m),
        "vmap" => cmd_vmap(&m),
        "nights" => cmd_nights(&m),
        "coldpool" => cmd_coldpool(&m),
        _ => die("usage: pf-weather constants | series | vmap | nights | coldpool (see source header)"),
    }
}
