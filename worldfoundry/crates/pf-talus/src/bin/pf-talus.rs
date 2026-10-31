//! pf-talus CLI (lane TALUS-1).
//!
//!   pf-talus block <block_dir> <x0_m> <y0_m> <out_dir> [facies=<u8 raster>] [window=<T1b window dir>]
//!                  [core=gx0,gy0,gx1,gy1] [halo_m=<m>|auto] [k=v ...]
//!   pf-talus window <window_dir> <gx0> <gy0> <gx1> <gy1> <out_dir> [k=v ...]      (32 m posts, coarse mode)
//!   pf-talus seam <block_dir> <x0_m> <y0_m> <out_dir> tiles=<nx>x<ny> [halo_m=..] [facies=..] [window=..] [k=v ...]
//!
//! `core` is an inclusive global post rectangle on the 2 m lattice (default: the whole block); the run reads the
//! block over core + halo (default `TalusParams::halo_m`, the exact-identity halo).

use pf_talus::model::{run, Rect, RunOut};
use pf_talus::output::{block_record, field_bytes, write_run};
use pf_talus::params::TalusParams;
use pf_talus::surface::{load_block, load_window_surface, Surface};
use pf_talus::{domain_for, full_rect};
use serde_json::json;
use std::path::{Path, PathBuf};

fn die(m: String) -> ! { eprintln!("{m}"); std::process::exit(1) }

struct Opts {
    p: TalusParams,
    facies: Option<PathBuf>,
    window: Option<PathBuf>,
    core: Option<Rect>,
    halo: Option<f64>,
    tiles: Option<(usize, usize)>,
}

fn opts(a: &[String]) -> Opts {
    let mut o = Opts { p: TalusParams::default(), facies: None, window: None, core: None, halo: None, tiles: None };
    for kv in a {
        let (k, v) = kv.split_once('=').unwrap_or_else(|| die(format!("bad arg {kv}")));
        match k {
            "facies" => o.facies = Some(PathBuf::from(v)),
            "window" => o.window = Some(PathBuf::from(v)),
            "core" => {
                let n: Vec<i64> = v.split(',').map(|x| x.parse().unwrap_or_else(|e| die(format!("core: {e}")))).collect();
                o.core = Some(Rect { gx0: n[0], gy0: n[1], gx1: n[2], gy1: n[3] });
            }
            "halo_m" => o.halo = if v == "auto" { None } else { Some(v.parse().unwrap_or_else(|e| die(format!("halo: {e}")))) },
            "tiles" => {
                let (a, b) = v.split_once('x').unwrap_or_else(|| die("tiles=NxM".into()));
                o.tiles = Some((a.parse().unwrap(), b.parse().unwrap()));
            }
            _ => o.p.set(k, v).unwrap_or_else(|e| die(e)),
        }
    }
    o
}

fn load(a: &[String], o: &Opts) -> Surface {
    let x0: f64 = a[1].parse().unwrap_or_else(|e| die(format!("x0: {e}")));
    let y0: f64 = a[2].parse().unwrap_or_else(|e| die(format!("y0: {e}")));
    load_block(Path::new(&a[0]), x0, y0, o.facies.as_deref(), o.window.as_deref(), o.p.fluvial_area_m2).unwrap_or_else(|e| die(e))
}

fn one(s: &Surface, o: &Opts, core: Rect) -> RunOut {
    let halo = o.halo.unwrap_or_else(|| o.p.halo_m(s.dx));
    let d = domain_for(s, core, halo);
    run(&d, &o.p, core).unwrap_or_else(|e| die(e))
}

fn summary(r: &RunOut) {
    let st = &r.stats;
    println!("domain {:?} core {:?} sources {} chutes {} blocks(core) {} visible {} relax lost {:.1} m3 residual {}/{} steps field {} block {} timing {:?}",
        st.domain, st.core, st.sources, st.chutes, r.blocks.len(), r.blocks.iter().filter(|b| b.visible).count(), st.relax_lost_m3,
        st.repose_residual_over, st.repose_residual_posts, st.field_traces_steps, st.block_traces_steps, st.timing_s);
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.is_empty() { die("usage: pf-talus block|window|seam ...".into()); }
    match a[0].as_str() {
        "block" => {
            if a.len() < 5 { die("block <block_dir> <x0_m> <y0_m> <out_dir> [k=v]".into()); }
            let o = opts(&a[5..]);
            let s = load(&a[1..4], &o);
            let core = o.core.unwrap_or_else(|| full_rect(&s));
            let r = one(&s, &o, core);
            summary(&r);
            let halo = o.halo.unwrap_or_else(|| o.p.halo_m(s.dx));
            write_run(Path::new(&a[4]), &r, &o.p, &s.digests, json!({"mode": "block", "halo_m": halo, "block": a[1], "origin_m": [a[2], a[3]]})).unwrap_or_else(|e| die(e));
        }
        "window" => {
            if a.len() < 7 { die("window <window_dir> <gx0> <gy0> <gx1> <gy1> <out_dir> [k=v]".into()); }
            let o = opts(&a[7..]);
            let g: Vec<i64> = a[2..6].iter().map(|x| x.parse().unwrap_or_else(|e| die(format!("{e}")))).collect();
            let s = load_window_surface(Path::new(&a[1]), g[0], g[1], g[2], g[3], o.p.fluvial_area_m2).unwrap_or_else(|e| die(e));
            let core = o.core.unwrap_or_else(|| full_rect(&s));
            let r = one(&s, &o, core);
            summary(&r);
            write_run(Path::new(&a[6]), &r, &o.p, &s.digests, json!({"mode": "window-32m"})).unwrap_or_else(|e| die(e));
        }
        "seam" => {
            if a.len() < 5 { die("seam <block_dir> <x0_m> <y0_m> <out_dir> tiles=NxM [k=v]".into()); }
            let o = opts(&a[5..]);
            let s = load(&a[1..4], &o);
            let (tx, ty) = o.tiles.unwrap_or((2, 2));
            let halo = o.halo.unwrap_or_else(|| o.p.halo_m(s.dx));
            let out = PathBuf::from(&a[4]);
            std::fs::create_dir_all(&out).unwrap();
            let full = full_rect(&s);
            let fr = one(&s, &o, full);
            summary(&fr);
            let ffields = field_bytes(&fr);
            // tile cores: edge-inclusive, adjacent tiles share one row / column of posts
            let cuts = |n: usize, t: usize| -> Vec<usize> { (0..=t).map(|k| k * (n - 1) / t).collect() };
            let (cx, cy) = (cuts(s.nx, tx), cuts(s.ny, ty));
            let mut report = vec![];
            let mut all_ok = true;
            let mut tile_fields: Vec<Vec<Vec<u8>>> = vec![];
            let mut cores = vec![];
            for j in 0..ty {
                for i in 0..tx {
                    let core = Rect { gx0: s.gx0 + cx[i] as i64, gy0: s.gy0 + cy[j] as i64, gx1: s.gx0 + cx[i + 1] as i64, gy1: s.gy0 + cy[j + 1] as i64 };
                    let tr = one(&s, &o, core);
                    let tf = field_bytes(&tr);
                    let mut fdiff = vec![];
                    // compare every field with the same posts of the full run
                    for (f_full, f_tile) in ffields.iter().zip(tf.iter()) {
                        let bpp = if f_full.1 == "i32" { 4 } else { 2 };
                        let mut nd = 0u64;
                        let mut maxd = 0i64;
                        for r in 0..tr.ny {
                            for c in 0..tr.nx {
                                let kf = ((core.gy0 - s.gy0) as usize + r) * fr.nx + (core.gx0 - s.gx0) as usize + c;
                                let kt = r * tr.nx + c;
                                let a = &f_full.2[kf * bpp..kf * bpp + bpp];
                                let b = &f_tile.2[kt * bpp..kt * bpp + bpp];
                                if a != b {
                                    nd += 1;
                                    let va = if bpp == 4 { i32::from_le_bytes([a[0], a[1], a[2], a[3]]) as i64 } else { u16::from_le_bytes([a[0], a[1]]) as i64 };
                                    let vb = if bpp == 4 { i32::from_le_bytes([b[0], b[1], b[2], b[3]]) as i64 } else { u16::from_le_bytes([b[0], b[1]]) as i64 };
                                    maxd = maxd.max((va - vb).abs());
                                }
                            }
                        }
                        if nd > 0 { all_ok = false; }
                        fdiff.push(json!({"field": f_full.0, "posts_differ": nd, "max_abs_diff": maxd}));
                    }
                    let (ox0, ox1, oy0, oy1) = (core.gx0 as f64 * s.dx, core.gx1 as f64 * s.dx, core.gy0 as f64 * s.dx, core.gy1 as f64 * s.dx);
                    let want: Vec<Vec<u8>> = fr.blocks.iter().filter(|b| b.x >= ox0 && b.x < ox1 && b.y >= oy0 && b.y < oy1).map(block_record).collect();
                    let got: Vec<Vec<u8>> = tr.blocks.iter().map(block_record).collect();
                    let blocks_equal = want == got;
                    let mut first_diff = serde_json::Value::Null;
                    if !blocks_equal {
                        let fb: Vec<&pf_talus::model::Block> = fr.blocks.iter().filter(|b| b.x >= ox0 && b.x < ox1 && b.y >= oy0 && b.y < oy1).collect();
                        if let Some(m) = (0..want.len().min(got.len())).find(|&m| want[m] != got[m]) {
                            first_diff = json!({"full": format!("{:?}", fb[m]), "tile": format!("{:?}", tr.blocks[m])});
                            eprintln!("first block diff: full {:?}\n tile {:?}", fb[m], tr.blocks[m]);
                        }
                    }
                    if !blocks_equal { all_ok = false; }
                    println!("tile {i},{j} core {:?} domain {:?} fields {} blocks full {} tile {} equal {}", core, tr.stats.domain, fdiff.iter().map(|d| d["posts_differ"].as_u64().unwrap()).sum::<u64>(), want.len(), got.len(), blocks_equal);
                    report.push(json!({"tile": [i, j], "core": core, "domain": tr.stats.domain, "fields": fdiff, "blocks_full": want.len(), "blocks_tile": got.len(), "blocks_byte_equal": blocks_equal, "first_block_diff": first_diff}));
                    tile_fields.push(tf.into_iter().map(|f| f.2).collect());
                    cores.push(core);
                }
            }
            // shared edges between adjacent tiles (thickness field): the same posts computed by two different domains
            let mut edge_posts = 0u64;
            let mut edge_diff = 0u64;
            for a_i in 0..cores.len() {
                for b_i in 0..cores.len() {
                    let (ca, cb) = (cores[a_i], cores[b_i]);
                    let right = ca.gx1 == cb.gx0 && ca.gy0 == cb.gy0;
                    let up = ca.gy1 == cb.gy0 && ca.gx0 == cb.gx0;
                    if !(right || up) { continue; }
                    let (fa, fb) = (&tile_fields[a_i][0], &tile_fields[b_i][0]);
                    let (wa, wb) = ((ca.gx1 - ca.gx0 + 1) as usize, (cb.gx1 - cb.gx0 + 1) as usize);
                    let pairs: Vec<(usize, usize)> = if right {
                        (0..(ca.gy1 - ca.gy0 + 1) as usize).map(|r| (r * wa + wa - 1, r * wb)).collect()
                    } else {
                        let ha = (ca.gy1 - ca.gy0) as usize;
                        (0..wa).map(|c| (ha * wa + c, c)).collect()
                    };
                    for (ka, kb) in pairs {
                        edge_posts += 1;
                        if fa[ka * 4..ka * 4 + 4] != fb[kb * 4..kb * 4 + 4] { edge_diff += 1; }
                    }
                }
            }
            if edge_diff > 0 { all_ok = false; }
            println!("seam: tiles {tx}x{ty} halo {halo} m: all identical {all_ok}; shared-edge posts {edge_posts} differ {edge_diff}");
            let rep = json!({"block": a[1], "tiles": [tx, ty], "halo_m": halo, "halo_exact_m": o.p.halo_m(s.dx), "all_identical": all_ok, "shared_edge_posts": edge_posts, "shared_edge_differ": edge_diff, "full_blocks": fr.blocks.len(), "tiles_report": report});
            std::fs::write(out.join("seam.json"), serde_json::to_vec_pretty(&rep).unwrap()).unwrap();
        }
        _ => die(format!("unknown command {}", a[0])),
    }
}
