//! pf-soil CLI (lane SOIL-VEG-1).
//!   pf-soil window <t1b_window_dir> <out_dir> climate=<seed-climate-v1 dir> t0=<T0 run dir> [full=1] [threads=N] [<param>=<value> ...]
//!   pf-soil block <pf-detail block dir> <t1b_window_dir> <out_dir> climate=<dir> t0=<dir> [talus=<pf-talus out dir>]
//!                 [origin=x_m,y_m] [crop=i0,j0,ni,nj] [threads=N] [<param>=<value> ...]
//! window: publishes the window core (edge-inclusive 513 x 513 posts), or with full=1 the whole lattice (stencil-invalid
//! posts as no-data). block: the block lattice with a no-data ring of ring_2m_posts.

use pf_soil::output::{quantize, stats, write, FIELD_NOTES};
use pf_soil::params::SoilParams;
use pf_soil::{classes_json, recipe_digest, rules_json, tier2, tier32, window::Window, zonation::LIMIT_NAMES, GENERATOR, GENERATOR_VERSION, LABEL, SCHEMA};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn main() {
    if let Err(e) = run() { eprintln!("pf-soil: {e}"); std::process::exit(1); }
}

struct Kv { climate: Option<PathBuf>, t0: Option<PathBuf>, talus: Option<PathBuf>, origin: Option<(f64, f64)>, crop: Option<(usize, usize, usize, usize)>, full: bool, threads: usize, params: SoilParams, overrides: Vec<(String, f64)> }

fn parse_kv(args: &[String]) -> Result<Kv, String> {
    let mut kv = Kv { climate: None, t0: None, talus: None, origin: None, crop: None, full: false, threads: 1, params: SoilParams::v1(), overrides: vec![] };
    for a in args {
        let (k, v) = a.split_once('=').ok_or_else(|| format!("argument {a}: expected k=v"))?;
        match k {
            "climate" => kv.climate = Some(PathBuf::from(v)),
            "t0" => kv.t0 = Some(PathBuf::from(v)),
            "talus" => kv.talus = Some(PathBuf::from(v)),
            "full" => kv.full = v == "1",
            "threads" => kv.threads = v.parse().map_err(|_| "threads")?,
            "origin" => { let p: Vec<f64> = v.split(',').map(|s| s.parse().map_err(|_| "origin")).collect::<Result<_, _>>()?; if p.len() != 2 { return Err("origin=x,y".into()); } kv.origin = Some((p[0], p[1])); }
            "crop" => { let p: Vec<usize> = v.split(',').map(|s| s.parse().map_err(|_| "crop")).collect::<Result<_, _>>()?; if p.len() != 4 { return Err("crop=i0,j0,ni,nj".into()); } kv.crop = Some((p[0], p[1], p[2], p[3])); }
            _ => { let x: f64 = v.parse().map_err(|_| format!("param {k}: not a number"))?; kv.params.set(k, x)?; kv.overrides.push((k.to_string(), x)); }
        }
    }
    Ok(kv)
}

fn inputs_json(w: &Window, cl: &pf_soil::climate::Climate) -> Value {
    json!({"t1b_window": {"digest": w.digest, "recipe_digest": w.recipe_digest, "window": w.window_key,
            "fields_read": w.fields_read.iter().map(|(n, d)| json!({"name": n, "digest": d})).collect::<Vec<_>>(),
            "volcano_records": w.volcanoes.iter().map(|v| json!({"id": v.id, "units": v.units.len(), "youngest_unit_yr": if v.youngest_yr().is_finite() { json!(v.youngest_yr()) } else { Value::Null }, "age_my": v.age_my, "stage": v.stage})).collect::<Vec<_>>()},
          "climate": {"digest": cl.climate_digest, "tier": cl.tier, "lapse_rate_k_per_m": cl.lapse, "artifact_axial_tilt_deg": cl.axial_tilt_deg},
          "t0": {"terrain_digest": cl.terrain_digest, "hydrology_digest": cl.hydrology_digest, "sea_level_m": cl.sea}})
}

fn header(tier: &str, layout: Value, inputs: Value, kv: &Kv) -> Value {
    let recipe = json!({"schema": SCHEMA, "generator": {"name": GENERATOR, "version": GENERATOR_VERSION}, "tier": tier, "params": kv.params.table(),
        "rules": rules_json(), "inputs": inputs, "layout": layout});
    json!({"schema": SCHEMA, "label": LABEL, "recipe_digest": recipe_digest(&recipe), "recipe": recipe, "classes": classes_json(), "limit_codes": LIMIT_NAMES,
        "param_overrides": kv.overrides.iter().map(|(k, v)| json!({k.clone(): v})).collect::<Vec<_>>(),
        "field_notes": FIELD_NOTES.iter().map(|(n, d, u, t)| json!({"name": n, "dtype": d, "units": u, "note": t})).collect::<Vec<_>>()})
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 { return Err("usage: pf-soil window <window> <out> climate= t0= ... | pf-soil block <block> <window> <out> climate= t0= ...".into()); }
    let t_all = Instant::now();
    match args[1].as_str() {
        "window" => {
            if args.len() < 4 { return Err("pf-soil window <window_dir> <out_dir> k=v...".into()); }
            let kv = parse_kv(&args[4..])?;
            let w = Window::load(Path::new(&args[2]))?;
            let cl = pf_soil::climate::Climate::load(kv.climate.as_ref().ok_or("climate= required")?, kv.t0.as_ref().ok_or("t0= required")?)?;
            let t = Instant::now();
            let t32 = tier32::run(&kv.params, &w, &cl, kv.threads)?;
            let eval_s = t.elapsed().as_secs_f64();
            let n = w.n;
            let (o0, nc) = if kv.full { (0usize, n) } else { (w.halo, n - 2 * w.halo) };
            let pick = |o: usize| -> Option<usize> { let (i, j) = (o % nc + o0, o / nc + o0); let k = j * n + i; if t32.valid[k] { Some(k) } else { None } };
            let f = quantize(nc * nc, &pick, &t32.out, &t32.twi, &t32.r_sun);
            let out = PathBuf::from(&args[3]);
            let fields = write(&out, &f)?;
            let layout = json!({"lattice": "T1b lod-4 posts, 32 m", "order": "row-major, row 0 south, x east", "size": [nc, nc],
                "post_0_global": [w.gx0 + o0 as i64, w.gy0 + o0 as i64], "post_0_chart_m": [(w.gx0 + o0 as i64) * 32, (w.gy0 + o0 as i64) * 32],
                "extent": if kv.full { "whole window lattice (stencil-invalid ring = no data)" } else { "window core, edge-inclusive (shared edge rows with the neighbour windows)" }});
            let mut h = header("t1b_32m", layout, inputs_json(&w, &cl), &kv);
            h["fields"] = json!(fields);
            h["stats"] = stats(&f);
            h["not_digested"] = json!({"paths": {"window": args[2], "climate": kv.climate, "t0": kv.t0}, "eval_seconds": eval_s, "total_seconds": t_all.elapsed().as_secs_f64(), "threads": kv.threads});
            std::fs::write(out.join("soil.json"), serde_json::to_vec_pretty(&h).unwrap()).map_err(|e| e.to_string())?;
            println!("pf-soil window {} -> {} posts {}x{} recipe {} cover_mean {:.1}", args[2], out.display(), nc, nc, h["recipe_digest"], h["stats"]["cover_mean_permille"]);
        }
        "block" => {
            if args.len() < 5 { return Err("pf-soil block <block_dir> <window_dir> <out_dir> k=v...".into()); }
            let kv = parse_kv(&args[5..])?;
            let w = Window::load(Path::new(&args[3]))?;
            let cl = pf_soil::climate::Climate::load(kv.climate.as_ref().ok_or("climate= required")?, kv.t0.as_ref().ok_or("t0= required")?)?;
            let b = tier2::Block::load(Path::new(&args[2]), kv.origin, kv.crop)?;
            let talus = match &kv.talus { Some(d) => Some(tier2::Talus::load(d)?), None => None };
            let t = Instant::now();
            let t32 = tier32::run(&kv.params, &w, &cl, kv.threads)?;
            let s32 = t.elapsed().as_secs_f64();
            let t = Instant::now();
            let t2 = tier2::run(&kv.params, &w, &t32, &cl, &b, talus.as_ref(), kv.threads)?;
            let s2 = t.elapsed().as_secs_f64();
            let nn = t2.n_x * t2.n_y;
            let pick = |o: usize| -> Option<usize> { if t2.valid[o] { Some(o) } else { None } };
            let f = quantize(nn, &pick, &t2.out, &t2.twi, &t2.r_sun);
            let out = PathBuf::from(&args[4]);
            let fields = write(&out, &f)?;
            let layout = json!({"lattice": "pf-detail 2 m posts", "order": "row-major, row 0 south, x east", "size": [t2.n_x, t2.n_y],
                "post_0_chart_m": [t2.x0_m, t2.y0_m], "ring_no_data_posts": kv.params.ring_2m_posts, "crop": kv.crop.map(|c| vec![c.0, c.1, c.2, c.3])});
            let mut inputs = inputs_json(&w, &cl);
            inputs["block"] = json!({"files": b.digests.iter().map(|(n, d)| json!({"file": n, "digest": d})).collect::<Vec<_>>()});
            inputs["talus"] = match &talus { Some(tt) => json!({"debris_cover_permille": tt.digest}), None => json!("absent: debris factor 0 (no pf-talus run for this block)") };
            let mut h = header("detail_2m", layout, inputs, &kv);
            h["fields"] = json!(fields);
            h["stats"] = stats(&f);
            h["stats"]["talus_covered_posts"] = json!(t2.talus_posts);
            h["not_digested"] = json!({"paths": {"block": args[2], "window": args[3], "climate": kv.climate, "t0": kv.t0, "talus": kv.talus}, "tier32_seconds": s32, "tier2_seconds": s2, "total_seconds": t_all.elapsed().as_secs_f64(), "threads": kv.threads});
            std::fs::write(out.join("soil.json"), serde_json::to_vec_pretty(&h).unwrap()).map_err(|e| e.to_string())?;
            println!("pf-soil block {} -> {} posts {}x{} recipe {} cover_mean {:.1}", args[2], out.display(), t2.n_x, t2.n_y, h["recipe_digest"], h["stats"]["cover_mean_permille"]);
        }
        m => return Err(format!("unknown mode {m}")),
    }
    Ok(())
}
