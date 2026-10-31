//! pf-snow CLI (lane SNOWPACK-1): the seasonal snowpack of `pf_soil::snowpack` on pf-soil's lattices.
//!   pf-snow window <t1b_window_dir> <out_dir> climate=<seed-climate-v1 dir> t0=<T0 run dir> days=32,91,152,244 [full=1] [threads=N] [<param>=<v> ...]
//!   pf-snow block <pf-detail block dir> <t1b_window_dir> <out_dir> climate= t0= days=... [origin=x_m,y_m] [crop=i0,j0,ni,nj] [threads=N] [<param>=<v> ...]
//!   pf-snow column climate= t0= window=<t1b_window_dir> | (lat= t_sl= half= p= lapse=) [slopes=30] [series=z1,z2,..] [<param>=<v> ...]
//! `days` are 1-based calendar days (1 = 1 Jan, 32 = 1 Feb, 365 = 31 Dec); fields are named _dNNN with that number.
//! window: the window core (edge-inclusive) as pf-soil publishes it, or with full=1 the whole lattice. block: the block
//! lattice with a no-data ring of insol_baseline_2m_posts. column: flat-ground and aspect equilibrium lines (gates).

use pf_soil::climate::Climate;
use pf_soil::params::SoilParams;
use pf_soil::snowpack::*;
use pf_soil::tier2::Block;
use pf_soil::tier32::{chart_of, climate_at, par_map};
use pf_soil::window::{Window, POST_M};
use pf_soil::zonation::ohmura_t_ela;
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn main() {
    if let Err(e) = run() { eprintln!("pf-snow: {e}"); std::process::exit(1); }
}

struct Kv {
    climate: Option<PathBuf>, t0: Option<PathBuf>, window: Option<PathBuf>, origin: Option<(f64, f64)>, crop: Option<(usize, usize, usize, usize)>,
    full: bool, threads: usize, days: Vec<usize>, params: SnowParams, overrides: Vec<(String, f64)>, col: Vec<(String, f64)>, slopes: Vec<f64>, series: Vec<f64>,
}

fn nums(v: &str) -> Result<Vec<f64>, String> { v.split(',').map(|s| s.trim().parse::<f64>().map_err(|_| format!("bad number {s}"))).collect() }

fn parse_kv(args: &[String]) -> Result<Kv, String> {
    let mut kv = Kv { climate: None, t0: None, window: None, origin: None, crop: None, full: false, threads: 1, days: vec![], params: SnowParams::v1(),
        overrides: vec![], col: vec![], slopes: vec![30.0], series: vec![] };
    for a in args {
        let (k, v) = a.split_once('=').ok_or_else(|| format!("argument {a}: expected k=v"))?;
        match k {
            "climate" => kv.climate = Some(PathBuf::from(v)),
            "t0" => kv.t0 = Some(PathBuf::from(v)),
            "window" => kv.window = Some(PathBuf::from(v)),
            "full" => kv.full = v == "1",
            "threads" => kv.threads = v.parse().map_err(|_| "threads")?,
            "days" => { kv.days = nums(v)?.into_iter().map(|x| x as usize).collect(); if kv.days.iter().any(|&d| d < 1 || d > DAYS) { return Err("days are 1..365".into()); } }
            "slopes" => kv.slopes = nums(v)?,
            "series" => kv.series = nums(v)?,
            "origin" => { let p = nums(v)?; if p.len() != 2 { return Err("origin=x,y".into()); } kv.origin = Some((p[0], p[1])); }
            "crop" => { let p = nums(v)?; if p.len() != 4 { return Err("crop=i0,j0,ni,nj".into()); } kv.crop = Some((p[0] as usize, p[1] as usize, p[2] as usize, p[3] as usize)); }
            "lat" | "t_sl" | "half" | "p" | "lapse" => kv.col.push((k.to_string(), v.parse().map_err(|_| format!("{k}: not a number"))?)),
            _ => { let x: f64 = v.parse().map_err(|_| format!("param {k}: not a number"))?; kv.params.set(k, x)?; kv.overrides.push((k.to_string(), x)); }
        }
    }
    Ok(kv)
}

struct Fields { swe: Vec<Vec<u16>>, cover: Vec<Vec<u16>>, balance: Vec<i32>, snow_days: Vec<u16>, perennial: Vec<u8>, irr: Vec<u16>, ms: Vec<u16>, ms_clamped: usize }

fn quantize(p: &SnowParams, n_out: usize, ndays: usize, pick: &dyn Fn(usize) -> Option<usize>, out: &[Option<SnowOut>]) -> Fields {
    let mut f = Fields { swe: vec![vec![u16::MAX; n_out]; ndays], cover: vec![vec![u16::MAX; n_out]; ndays], balance: vec![i32::MIN; n_out],
        snow_days: vec![u16::MAX; n_out], perennial: vec![255; n_out], irr: vec![u16::MAX; n_out], ms: vec![u16::MAX; n_out], ms_clamped: 0 };
    for o in 0..n_out {
        let Some(k) = pick(o) else { continue };
        let Some(s) = &out[k] else { continue };
        for (i, &w) in s.swe.iter().enumerate() {
            f.swe[i][o] = w.round().clamp(0.0, 65534.0) as u16;
            f.cover[i][o] = (fsc(p, w) * 1000.0).round().clamp(0.0, 1000.0) as u16;
        }
        f.balance[o] = s.balance_mm.round().clamp(-2.0e9, 2.0e9) as i32;
        f.snow_days[o] = s.snow_days;
        f.perennial[o] = s.perennial as u8;
        f.irr[o] = (s.irr_ratio * 1000.0).round().clamp(0.0, 65534.0) as u16;
        f.ms[o] = (s.ms * 1000.0).round().clamp(0.0, 65534.0) as u16;
        if s.ms_clamped { f.ms_clamped += 1; }
    }
    f
}

fn write(dir: &Path, f: &Fields, days: &[usize]) -> Result<Vec<Value>, String> {
    std::fs::create_dir_all(dir.join("fields")).map_err(|e| e.to_string())?;
    let mut out = Vec::new();
    let mut put = |name: String, dt: &str, unit: &str, note: String, bytes: Vec<u8>| -> Result<(), String> {
        let file = format!("fields/{name}.{dt}");
        std::fs::write(dir.join(&file), &bytes).map_err(|e| format!("{file}: {e}"))?;
        out.push(json!({"name": name, "file": file, "dtype": dt, "units": unit, "note": note, "digest": format!("blake3:{}", blake3::hash(&bytes).to_hex())}));
        Ok(())
    };
    let u16b = |v: &[u16]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    for (i, d) in days.iter().enumerate() {
        put(format!("swe_mm_d{d:03}"), "u16", "mm w.e.", format!("snow water equivalent on calendar day {d} (1-based) of the last spin-up year; 65535 no data"), u16b(&f.swe[i]))?;
        put(format!("snow_cover_permille_d{d:03}"), "u16", "permille", format!("snow-cover fraction tanh(SWE / swe_cover_mm) x 1000 on calendar day {d}; 65535 no data"), u16b(&f.cover[i]))?;
    }
    put("balance_mm".into(), "i32", "mm w.e./yr", "annual climatic balance of a snow surface (snowfall - melt potential); > 0 above the model equilibrium line; i32::MIN no data".into(),
        f.balance.iter().flat_map(|x| x.to_le_bytes()).collect())?;
    put("snow_days".into(), "u16", "d", "days of the last spin-up year with fsc >= 0.5; 65535 no data".into(), u16b(&f.snow_days))?;
    put("perennial".into(), "u8", "flag", "1 = SWE > 0 on every day of the last spin-up year (perennial snow); 255 no data".into(), f.perennial.clone())?;
    put("melt_multiplier_permille".into(), "u16", "permille", "melt multiplier used x 1000 (melt_scale, or the Ohmura-anchored value of the post's climate); 65535 no data".into(), u16b(&f.ms))?;
    put("irradiance_ratio_permille".into(), "u16", "permille", "annual mean clear-sky direct irradiance slope / flat x 1000 (diagnostic); 65535 no data".into(), u16b(&f.irr))?;
    Ok(out)
}

fn stats(f: &Fields, days: &[usize], land: &dyn Fn(usize) -> bool) -> Value {
    let n = f.balance.len();
    let mut per = Vec::new();
    for (i, d) in days.iter().enumerate() {
        let (mut c, mut s, mut cn, mut half) = (0f64, 0f64, 0u64, 0u64);
        for o in 0..n { if f.cover[i][o] != u16::MAX && land(o) { c += f.cover[i][o] as f64; s += f.swe[i][o] as f64; cn += 1; if f.cover[i][o] >= 500 { half += 1; } } }
        per.push(json!({"day": d, "land_posts": cn, "cover_mean_permille": if cn > 0 { c / cn as f64 } else { 0.0 }, "swe_mean_mm": if cn > 0 { s / cn as f64 } else { 0.0 },
            "cover_ge_half_share": if cn > 0 { half as f64 / cn as f64 } else { 0.0 }}));
    }
    let pc = (0..n).filter(|&o| f.perennial[o] == 1 && land(o)).count();
    let ln = (0..n).filter(|&o| f.perennial[o] != 255 && land(o)).count();
    let msv: Vec<f64> = (0..n).filter(|&o| f.ms[o] != u16::MAX).map(|o| f.ms[o] as f64 / 1000.0).collect();
    let (mmin, mmax) = msv.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &v| (a.0.min(v), a.1.max(v)));
    json!({"days": per, "perennial_land_share": if ln > 0 { pc as f64 / ln as f64 } else { 0.0 }, "land_posts": ln,
        "melt_multiplier_min": mmin, "melt_multiplier_max": mmax, "melt_multiplier_clamped_posts": f.ms_clamped})
}

fn header(tier: &str, layout: Value, inputs: Value, kv: &Kv) -> Value {
    let recipe = json!({"schema": SCHEMA, "generator": {"name": GENERATOR, "version": GENERATOR_VERSION}, "tier": tier, "params": kv.params.table(),
        "rules": rules_json(), "days": kv.days, "inputs": inputs, "layout": layout});
    json!({"schema": SCHEMA, "label": LABEL, "recipe_digest": pf_soil::recipe_digest(&recipe), "recipe": recipe,
        "param_overrides": kv.overrides.iter().map(|(k, v)| json!({k.clone(): v})).collect::<Vec<_>>()})
}

fn inputs_json(w: &Window, cl: &Climate) -> Value {
    json!({"t1b_window": {"digest": w.digest, "recipe_digest": w.recipe_digest, "window": w.window_key,
            "fields_read": w.fields_read.iter().filter(|(n, _)| n == "height_mm" || n == "fixed").map(|(n, d)| json!({"name": n, "digest": d})).collect::<Vec<_>>()},
          "climate": {"digest": cl.climate_digest, "tier": cl.tier, "lapse_rate_k_per_m": cl.lapse, "artifact_axial_tilt_deg": cl.axial_tilt_deg},
          "t0": {"terrain_digest": cl.terrain_digest, "hydrology_digest": cl.hydrology_digest, "sea_level_m": cl.sea}})
}

fn load_cl(kv: &Kv) -> Result<Climate, String> {
    let cl = Climate::load(kv.climate.as_ref().ok_or("climate= required")?, kv.t0.as_ref().ok_or("t0= required")?)?;
    if let Some(t) = cl.axial_tilt_deg { if (t - kv.params.tilt_deg).abs() > 1e-9 { return Err(format!("climate axial tilt {t} != tilt_deg {}; pass tilt_deg={t}", kv.params.tilt_deg)); } }
    Ok(cl)
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 { return Err("usage: pf-snow window|block|column ...".into()); }
    let t_all = Instant::now();
    match args[1].as_str() {
        "window" => {
            if args.len() < 4 { return Err("pf-snow window <window_dir> <out_dir> k=v...".into()); }
            let kv = parse_kv(&args[4..])?;
            if kv.days.is_empty() { return Err("days= required".into()); }
            let w = Window::load(Path::new(&args[2]))?;
            let cl = load_cl(&kv)?;
            let p = &kv.params;
            let cal = Calendar::new(p);
            let n = w.n;
            let nn = n * n;
            let chart = chart_of(&w)?;
            let pm = (POST_M * 1000.0) as i64;
            let z = &w.height_m;
            let d0: Vec<usize> = kv.days.iter().map(|d| d - 1).collect();
            let t = Instant::now();
            let res: Vec<Option<SnowOut>> = par_map(nn, kv.threads, |k| {
                let (i, j) = (k % n, k / n);
                if i < 1 || j < 1 || i + 1 >= n || j + 1 >= n { return None; }
                let (gx, gy) = w.g(k);
                let (la, t_sl, half, pp, _) = climate_at(&cl, &chart, gx * pm, gy * pm).ok()?;
                let zx = (z[k + 1] - z[k - 1]) / (2.0 * POST_M);
                let zy = (z[k + n] - z[k - n]) / (2.0 * POST_M);
                let x = SnowIn { lat_rad: la, t_sl, half, p_mm: pp, lapse: cl.lapse, h: z[k], sea: cl.sea, zx, zy };
                Some(post(p, &cal, &x, &d0))
            });
            let eval_s = t.elapsed().as_secs_f64();
            let (o0, nc) = if kv.full { (0usize, n) } else { (w.halo, n - 2 * w.halo) };
            if !kv.full && res.iter().enumerate().any(|(k, r)| r.is_none() && { let (i, j) = (k % n, k / n); i >= o0 && j >= o0 && i < o0 + nc && j < o0 + nc }) {
                return Err("a core post has no snowpack (chart inverse failed)".into());
            }
            let pick = |o: usize| -> Option<usize> { let (i, j) = (o % nc + o0, o / nc + o0); Some(j * n + i) };
            let f = quantize(p, nc * nc, kv.days.len(), &pick, &res);
            let out = PathBuf::from(&args[3]);
            let fields = write(&out, &f, &kv.days)?;
            let layout = json!({"lattice": "T1b lod-4 posts, 32 m", "order": "row-major, row 0 south, x east", "size": [nc, nc],
                "post_0_global": [w.gx0 + o0 as i64, w.gy0 + o0 as i64], "post_0_chart_m": [(w.gx0 + o0 as i64) * 32, (w.gy0 + o0 as i64) * 32],
                "extent": if kv.full { "whole window lattice (outer ring = no data)" } else { "window core, edge-inclusive (shared edge rows with the neighbour windows)" }});
            let mut h = header("t1b_32m", layout, inputs_json(&w, &cl), &kv);
            h["fields"] = json!(fields);
            let land = |o: usize| { let (i, j) = (o % nc + o0, o / nc + o0); let k = j * n + i; w.fixed[k] != 2 && w.fixed[k] != 4 && z[k] > cl.sea };
            h["stats"] = stats(&f, &kv.days, &land);
            h["not_digested"] = json!({"paths": {"window": args[2], "climate": kv.climate, "t0": kv.t0}, "eval_seconds": eval_s, "total_seconds": t_all.elapsed().as_secs_f64(), "threads": kv.threads});
            std::fs::write(out.join("snow.json"), serde_json::to_vec_pretty(&h).unwrap()).map_err(|e| e.to_string())?;
            println!("pf-snow window {} -> {} posts {}x{} recipe {} eval {:.1}s", args[2], out.display(), nc, nc, h["recipe_digest"], eval_s);
        }
        "block" => {
            if args.len() < 5 { return Err("pf-snow block <block_dir> <window_dir> <out_dir> k=v...".into()); }
            let kv = parse_kv(&args[5..])?;
            if kv.days.is_empty() { return Err("days= required".into()); }
            let w = Window::load(Path::new(&args[3]))?;
            let cl = load_cl(&kv)?;
            let b = Block::load(Path::new(&args[2]), kv.origin, kv.crop)?;
            let p = &kv.params;
            let cal = Calendar::new(p);
            let chart = chart_of(&w)?;
            let (nx, ny) = (b.n_x, b.n_y);
            let bl = p.insol_baseline_2m_posts.max(1.0) as usize;
            let z = &b.height_m;
            let d0: Vec<usize> = kv.days.iter().map(|d| d - 1).collect();
            let t = Instant::now();
            let res: Vec<Option<SnowOut>> = par_map(nx * ny, kv.threads, |k| {
                let (i, j) = (k % nx, k / nx);
                if i < bl || j < bl || i + bl >= nx || j + bl >= ny { return None; }
                let x = b.x0_m + 2.0 * i as f64;
                let y = b.y0_m + 2.0 * j as f64;
                let (la, t_sl, half, pp, _) = climate_at(&cl, &chart, (x * 1000.0).round() as i64, (y * 1000.0).round() as i64).ok()?;
                let zx = (z[k + bl] - z[k - bl]) / (4.0 * bl as f64);
                let zy = (z[k + bl * nx] - z[k - bl * nx]) / (4.0 * bl as f64);
                let xin = SnowIn { lat_rad: la, t_sl, half, p_mm: pp, lapse: cl.lapse, h: z[k], sea: cl.sea, zx, zy };
                Some(post(p, &cal, &xin, &d0))
            });
            let eval_s = t.elapsed().as_secs_f64();
            let pick = |o: usize| -> Option<usize> { Some(o) };
            let f = quantize(p, nx * ny, kv.days.len(), &pick, &res);
            let out = PathBuf::from(&args[4]);
            let fields = write(&out, &f, &kv.days)?;
            let layout = json!({"lattice": "pf-detail 2 m posts", "order": "row-major, row 0 south, x east", "size": [nx, ny],
                "post_0_chart_m": [b.x0_m, b.y0_m], "ring_no_data_posts": bl, "crop": kv.crop.map(|c| vec![c.0, c.1, c.2, c.3])});
            let mut inputs = inputs_json(&w, &cl);
            inputs["block"] = json!({"files": b.digests.iter().map(|(n, d)| json!({"file": n, "digest": d})).collect::<Vec<_>>()});
            let mut h = header("detail_2m", layout, inputs, &kv);
            h["fields"] = json!(fields);
            let land = |o: usize| b.fixed[o] != 1 && b.fixed[o] != 3 && z[o] > cl.sea;
            h["stats"] = stats(&f, &kv.days, &land);
            h["not_digested"] = json!({"paths": {"block": args[2], "window": args[3], "climate": kv.climate, "t0": kv.t0}, "eval_seconds": eval_s, "total_seconds": t_all.elapsed().as_secs_f64(), "threads": kv.threads});
            std::fs::write(out.join("snow.json"), serde_json::to_vec_pretty(&h).unwrap()).map_err(|e| e.to_string())?;
            println!("pf-snow block {} -> {} posts {}x{} recipe {} eval {:.1}s", args[2], out.display(), nx, ny, h["recipe_digest"], eval_s);
        }
        "column" => {
            let kv = parse_kv(&args[2..])?;
            let p = &kv.params;
            let cal = Calendar::new(p);
            let getc = |k: &str| kv.col.iter().find(|(a, _)| a == k).map(|x| x.1);
            let (lat_deg, t_sl, half, pmm, lapse, src) = if let Some(wd) = &kv.window {
                let w = Window::load(wd)?;
                let cl = load_cl(&kv)?;
                let chart = chart_of(&w)?;
                let c = w.gx0 + (w.n as i64) / 2;
                let r = w.gy0 + (w.n as i64) / 2;
                let (la, t_sl, half, pp, _) = climate_at(&cl, &chart, c * 32_000, r * 32_000)?;
                (la.to_degrees(), t_sl, half, pp, cl.lapse, json!({"window": wd, "centre_post": [c, r]}))
            } else {
                (getc("lat").ok_or("lat=")?, getc("t_sl").ok_or("t_sl=")?, getc("half").ok_or("half=")?, getc("p").ok_or("p=")?, getc("lapse").unwrap_or(0.0065), json!("explicit"))
            };
            let la = lat_deg.to_radians();
            // Ohmura line exactly as pf-soil computes nival_margin = 0 on flat ground
            let sp = SoilParams::v1();
            let t_jja = t_sl + sp.jja_factor * half;
            let ohm = (t_jja - ohmura_t_ela(&sp, pmm)) / lapse;
            let flat = ela(p, &cal, la, t_sl, half, pmm, lapse, 0.0, 0.0, -6000.0, 12000.0);
            let mut pf = p.clone();
            pf.anchor_ohmura = 0.0;
            let free = ela(&pf, &cal, la, t_sl, half, pmm, lapse, 0.0, 0.0, -6000.0, 12000.0);
            let (msv, msc) = melt_multiplier(p, &cal, la, t_sl, half, pmm, lapse);
            let mut asp = Vec::new();
            for &s in &kv.slopes {
                let mut row = serde_json::Map::new();
                for (nm, a) in [("N", 0.0), ("NE", 45.0), ("E", 90.0), ("SE", 135.0), ("S", 180.0), ("SW", 225.0), ("W", 270.0), ("NW", 315.0)] {
                    let (zx, zy) = facing(s, a);
                    row.insert(nm.to_string(), json!(ela(p, &cal, la, t_sl, half, pmm, lapse, zx, zy, -6000.0, 12000.0)));
                }
                asp.push(json!({"slope_deg": s, "ela_m_by_facing": row}));
            }
            let mut ser = Vec::new();
            for &zz in &kv.series {
                let x = SnowIn { lat_rad: la, t_sl: t_sl - lapse * zz.min(0.0), half, p_mm: pmm, lapse, h: zz.max(0.0), sea: 0.0, zx: 0.0, zy: 0.0 };
                let all: Vec<usize> = (0..DAYS).collect();
                let o = post(p, &cal, &x, &all);
                ser.push(json!({"z_m": zz, "swe_mm": o.swe.iter().map(|v| (v * 10.0).round() / 10.0).collect::<Vec<_>>(), "balance_mm": o.balance_mm, "perennial": o.perennial, "snow_days": o.snow_days}));
            }
            let outj = json!({"source": src, "lat_deg": lat_deg, "t_sl_c": t_sl, "half_k": half, "p_mm": pmm, "lapse_k_m": lapse, "t_jja_sl_c": t_jja,
                "ohmura_ela_m": ohm, "model_ela_flat_m": flat, "model_minus_ohmura_m": flat.map(|f| f - ohm),
                "free_model_ela_flat_m": free, "free_minus_ohmura_m": free.map(|f| f - ohm), "free_melt_scale": p.melt_scale, "melt_multiplier": msv, "melt_multiplier_clamped": msc, "aspect": asp, "series": ser,
                "params": p.table(), "overrides": kv.overrides.iter().map(|(k, v)| json!({k.clone(): v})).collect::<Vec<_>>()});
            println!("{}", serde_json::to_string(&outj).unwrap());
        }
        m => return Err(format!("unknown mode {m}")),
    }
    Ok(())
}
