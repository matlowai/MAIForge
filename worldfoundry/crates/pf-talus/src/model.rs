//! One talus run over a domain (a Surface already cropped to core + halo), outputs restricted to the core.

use crate::chute::{chute_field, pack};
use crate::params::{class_fractions, rock_of, Rock, TalusParams, LAW_SET};
use crate::relax::{relax_multilevel, residual};
use crate::runout::{trace, Kin, Stop, Terrain, TraceCfg};
use crate::surface::{Surface, SINK_NONE, SINK_WATER};
use serde::Serialize;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicI64, Ordering};
use std::time::Instant;
use vz_keyed_rng::Scope;

/// Volume unit of every budget integer: 1e-6 m3.
pub const VU: f64 = 1.0e-6;
/// Distance histogram bin (m) and bins.
pub const HBIN_M: f64 = 2.0;
pub const HBINS: usize = 400;

/// Inclusive global post rectangle.
#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
pub struct Rect {
    pub gx0: i64,
    pub gy0: i64,
    pub gx1: i64,
    pub gy1: i64,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct Block {
    /// Stable id: keyed draw of (world seed, law set, lattice, source post world index, size class, event).
    pub id: u64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub d_m: f64,
    pub src_gx: i64,
    pub src_gy: i64,
    pub class: u8,
    pub event: u32,
    pub litho: u8,
    pub facies: u8,
    pub unit: u32,
    pub d_src_m: f64,
    pub d_foot_m: f64,
    pub path_m: f64,
    pub burial_m: f64,
    pub visible: bool,
    /// 0 rest, 1 sink (see sink), 2 off-domain, 3 capped.
    pub stop: u8,
    pub sink: u8,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct UnitBudget {
    pub litho: u8,
    pub facies: u8,
    pub sources: u64,
    /// Released rock (solid), weathering loss, rock left after weathering (all 1e-6 m3).
    pub rock_u: i64,
    pub weathered_u: i64,
    /// Field classes (< d_point): bulk (= rock / (1 - n)) volume released, and where it went.
    pub field_bulk_u: i64,
    pub field_dep_u: i64,
    pub field_water_u: i64,
    pub field_off_u: i64,
    pub field_cap_u: i64,
    /// Block classes: rock volume the size distribution assigns, rock volume actually emitted as blocks.
    pub block_class_rock_u: i64,
    pub block_rock_u: i64,
    pub blocks: u64,
    pub blocks_rest: u64,
    pub blocks_sink: u64,
    pub blocks_off: u64,
    pub blocks_cap: u64,
    pub blocks_visible: u64,
}

impl UnitBudget {
    fn add(&mut self, o: &UnitBudget) {
        self.sources += o.sources;
        self.rock_u += o.rock_u;
        self.weathered_u += o.weathered_u;
        self.field_bulk_u += o.field_bulk_u;
        self.field_dep_u += o.field_dep_u;
        self.field_water_u += o.field_water_u;
        self.field_off_u += o.field_off_u;
        self.field_cap_u += o.field_cap_u;
        self.block_class_rock_u += o.block_class_rock_u;
        self.block_rock_u += o.block_rock_u;
        self.blocks += o.blocks;
        self.blocks_rest += o.blocks_rest;
        self.blocks_sink += o.blocks_sink;
        self.blocks_off += o.blocks_off;
        self.blocks_cap += o.blocks_cap;
        self.blocks_visible += o.blocks_visible;
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ClassStats {
    pub lo_m: f64,
    pub hi_m: f64,
    pub block_class: bool,
    pub parcels: u64,
    pub bulk_u: i64,
    pub rest_u: i64,
    /// Volume-weighted histograms (units) of horizontal distance from the face foot and from the source, HBIN_M bins.
    pub h_foot: Vec<i64>,
    pub h_src: Vec<i64>,
}

#[derive(Debug, Serialize)]
pub struct RunStats {
    pub domain: Rect,
    pub core: Rect,
    pub dx: f64,
    pub sources: u64,
    pub chutes: u64,
    pub chute_posts: u64,
    pub relax_lost_m3: f64,
    /// Relaxation throughput (residual against the gross moved volume and against the stock).
    pub relax_throughput: crate::budget::RelaxThroughput,
    pub relax_iters: u32,
    pub repose_residual_over: u64,
    pub repose_residual_posts: u64,
    pub field_traces_steps: u64,
    pub block_traces_steps: u64,
    pub timing_s: BTreeMap<String, f64>,
}

pub struct RunOut {
    pub core: Rect,
    pub nx: usize,
    pub ny: usize,
    pub thickness_mm: Vec<i32>,
    pub chute_mm: Vec<i32>,
    pub source_mm: Vec<i32>,
    pub d50_mm: Vec<u16>,
    pub cover_pm: Vec<u16>,
    pub cobble_per_100m2: Vec<u16>,
    pub blocks: Vec<Block>,
    /// Budget per source unit: key = (exposed-unit id << 8) | volcanic facies, so a split unit reports lava, ash and
    /// riser lava separately.
    pub budget: BTreeMap<u64, UnitBudget>,
    pub classes: Vec<ClassStats>,
    pub stats: RunStats,
}

fn grad(nx: usize, ny: usize, dx: f64, h: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = nx * ny;
    let mut gx = vec![0.0; n];
    let mut gy = vec![0.0; n];
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            gx[k] = if i == 0 { (h[k + 1] - h[k]) / dx } else if i == nx - 1 { (h[k] - h[k - 1]) / dx } else { (h[k + 1] - h[k - 1]) / (2.0 * dx) };
            gy[k] = if j == 0 { (h[k + nx] - h[k]) / dx } else if j == ny - 1 { (h[k] - h[k - nx]) / dx } else { (h[k + nx] - h[k - nx]) / (2.0 * dx) };
        }
    }
    (gx, gy)
}

fn u16s(w: u64, slot: u32) -> f64 { ((w >> (16 * slot)) & 0xffff) as f64 / 65536.0 }

/// Cell and weights of a position given as source post (oi, oj) + offset (px, py) posts (domain independent weights).
fn cell_rel(nx: usize, ny: usize, oi: i64, oj: i64, px: f64, py: f64) -> (usize, usize, f64, f64) {
    let (fpx, fpy) = (px.floor(), py.floor());
    let mut i = oi + fpx as i64;
    let mut j = oj + fpy as i64;
    let (mut a, mut b) = (px - fpx, py - fpy);
    if i > nx as i64 - 2 { i = nx as i64 - 2; a = 1.0; }
    if j > ny as i64 - 2 { j = ny as i64 - 2; b = 1.0; }
    if i < 0 { i = 0; a = 0.0; }
    if j < 0 { j = 0; b = 0.0; }
    (i as usize, j as usize, a, b)
}

/// Integer bilinear split of `u` units at a source-relative position; remainder to the nearest corner. Exact.
fn split4(nx: usize, ny: usize, oi: i64, oj: i64, px: f64, py: f64, u: i64) -> [(usize, i64); 4] {
    let (i, j, a, b) = cell_rel(nx, ny, oi, oj, px, py);
    let k = j * nx + i;
    let w = [(1.0 - a) * (1.0 - b), a * (1.0 - b), (1.0 - a) * b, a * b];
    let ks = [k, k + 1, k + nx, k + nx + 1];
    let mut parts = [0i64; 4];
    let mut s = 0i64;
    for m in 0..4 { parts[m] = (u as f64 * w[m]).floor() as i64; s += parts[m]; }
    let near = if a < 0.5 { if b < 0.5 { 0 } else { 2 } } else if b < 0.5 { 1 } else { 3 };
    parts[near] += u - s;
    [(ks[0], parts[0]), (ks[1], parts[1]), (ks[2], parts[2]), (ks[3], parts[3])]
}

fn bil(nx: usize, ny: usize, v: &[f64], oi: i64, oj: i64, px: f64, py: f64) -> f64 {
    let (i, j, a, b) = cell_rel(nx, ny, oi, oj, px, py);
    let k = j * nx + i;
    (1.0 - a) * (1.0 - b) * v[k] + a * (1.0 - b) * v[k + 1] + (1.0 - a) * b * v[k + nx] + a * b * v[k + nx + 1]
}

struct Source {
    k: usize,
    rock_m3: f64,
    rock: Rock,
    fr: Vec<f64>,
}

#[derive(Default)]
struct ChunkOut {
    budget: BTreeMap<u64, UnitBudget>,
    cls_parcels: Vec<u64>,
    cls_bulk: Vec<i64>,
    cls_rest: Vec<i64>,
    h_foot: Vec<Vec<i64>>,
    h_src: Vec<Vec<i64>>,
    blocks: Vec<Block>,
    steps: u64,
}

impl ChunkOut {
    fn new(nc: usize) -> ChunkOut {
        ChunkOut { cls_parcels: vec![0; nc], cls_bulk: vec![0; nc], cls_rest: vec![0; nc], h_foot: vec![vec![0; HBINS]; nc], h_src: vec![vec![0; HBINS]; nc], ..Default::default() }
    }
}

fn hbin(d: f64) -> usize { ((d / HBIN_M) as usize).min(HBINS - 1) }

pub fn scope_for(p: &TalusParams, dx: f64) -> Result<Scope, String> {
    Scope::root(p.seed)
        .law_set(LAW_SET).map_err(|e| e.to_string())?
        .chart(0).map_err(|e| e.to_string())?
        .address(&format!("lattice-{}mm", (dx * 1000.0).round() as i64)).map_err(|e| e.to_string())
}

pub fn run(s: &Surface, p: &TalusParams, core: Rect) -> Result<RunOut, String> {
    let t_all = Instant::now();
    let mut timing = BTreeMap::new();
    let (nx, ny, dx) = (s.nx, s.ny, s.dx);
    let n = s.n();
    if nx < 3 || ny < 3 { return Err("domain too small".into()); }
    let scope = scope_for(p, dx)?;
    let nc = p.n_classes();
    let tan_lo = p.src_lo_deg.to_radians().tan();
    let threads = p.threads.max(1);

    // 1. input slope, producing posts, chutes, routing surface
    let t0 = Instant::now();
    let (gx0, gy0) = grad(nx, ny, dx, &s.z);
    let tan_in: Vec<f64> = gx0.iter().zip(gy0.iter()).map(|(a, b)| (a * a + b * b).sqrt()).collect();
    let tan_talus_rock = p.talus_rock_deg.to_radians().tan();
    let producing: Vec<bool> = (0..n)
        .map(|k| s.sink[k] == SINK_NONE && rock_of(s.litho[k], s.facies[k]).is_some() && !(s.regime[k] == 1 && tan_in[k] < tan_talus_rock))
        .collect();
    let (chute, n_chutes) = chute_field(s, &tan_in, &producing, p, &scope)?;
    let zr: Vec<f64> = s.z.iter().zip(chute.iter()).map(|(z, c)| z - c).collect();
    let (gxr, gyr) = grad(nx, ny, dx, &zr);
    let rock_tan: Vec<f64> = gxr.iter().zip(gyr.iter()).map(|(a, b)| (a * a + b * b).sqrt()).collect();
    let soil = |k: usize| -> bool { s.litho[k] != 255 && (s.regime[k] == 0 || s.regime[k] == 3 || s.litho[k] == 10) && !(s.regime[k] == 0 && rock_tan[k] > tan_lo) };
    let ground0: Vec<u8> = (0..n).map(|k| soil(k) as u8).collect();
    timing.insert("prep".into(), t0.elapsed().as_secs_f64());

    // 2. sources (input slope; production weight ramps between src_lo and src_hi)
    let mut sources: Vec<Source> = Vec::new();
    let mut source_depth = vec![0.0f64; n];
    for k in 0..n {
        if !producing[k] || tan_in[k] < tan_lo { continue; }
        let th = tan_in[k].atan().to_degrees();
        let w = ((th - p.src_lo_deg) / (p.src_hi_deg - p.src_lo_deg)).clamp(0.0, 1.0);
        if w <= 0.0 { continue; }
        let rock = rock_of(s.litho[k], s.facies[k]).unwrap();
        let depth = rock.retreat_mm_yr * 1.0e-3 * p.exposure_yr * w;
        let cos = (1.0 / (1.0 + tan_in[k] * tan_in[k])).sqrt().max(p.cos_floor);
        let rock_m3 = depth * dx * dx / cos;
        source_depth[k] = depth;
        let fr = class_fractions(p, &rock);
        sources.push(Source { k, rock_m3, rock, fr });
    }
    let chunk = 2048usize;
    let n_chunks = sources.len().div_ceil(chunk);
    let kin_of = |d: f64, rock: &Rock, wbits: u64| -> Kin {
        let mult = (p.mu_spread_ln * (2.0 * u16s(wbits, 2) - 1.0)).exp();
        let phi = p.phi_d_deg(d, rock.phi_extra_deg);
        let tn = |deg: f64| deg.min(80.0).to_radians().tan() * mult;
        let mu_d = [tn(phi), tn(phi + p.soil_phi_extra_deg)];
        let mu_s = [tn(phi + p.phi_static_extra_deg), tn(phi + p.soil_phi_extra_deg + p.phi_static_extra_deg)];
        let rt = (p.rt(d) * (1.0 + 0.1 * (2.0 * u16s(wbits, 3) - 1.0))).clamp(0.0, 1.0);
        Kin { mu_d, mu_s, rt, seed: wbits }
    };
    let cfg = TraceCfg {
        tan_ff: p.ff_deg.to_radians().tan(),
        tan_foot: tan_lo,
        mu_face: p.mu_face,
        g: p.g,
        l_max_m: p.l_max_m,
        dt_max_s: p.dt_max_s,
        step_posts: p.step_posts,
        max_steps: p.max_steps,
        lat_impact_rad: p.lat_impact_deg.to_radians(),
        lat_roll_rad: p.lat_roll_deg.to_radians(),
        turn_every_m: p.turn_every_m,
    };

    // 3. field classes on the routing surface
    let t1 = Instant::now();
    let acc_v: Vec<AtomicI64> = (0..n).map(|_| AtomicI64::new(0)).collect();
    let acc_m: Vec<AtomicI64> = (0..n).map(|_| AtomicI64::new(0)).collect();
    let acc_c: Vec<AtomicI64> = (0..n).map(|_| AtomicI64::new(0)).collect();
    let terr = Terrain { nx, ny, dx, gx: &gxr, gy: &gyr, sink: &s.sink, ground: &ground0, rock_tan: &rock_tan };
    let field_out: Vec<ChunkOut> = {
        let slots: Vec<std::sync::Mutex<Option<ChunkOut>>> = (0..n_chunks).map(|_| std::sync::Mutex::new(None)).collect();
        std::thread::scope(|sc| {
            for t in 0..threads {
                let (sources, slots, terr, acc_v, acc_m, acc_c, scope, cfg) = (&sources, &slots, &terr, &acc_v, &acc_m, &acc_c, &scope, &cfg);
                let kin_of = &kin_of;
                sc.spawn(move || {
                    let mut ci = t;
                    while ci < n_chunks {
                        let mut o = ChunkOut::new(nc);
                        for src in &sources[ci * chunk..((ci + 1) * chunk).min(sources.len())] {
                            let (i, j) = (src.k % nx, src.k / nx);
                            let (sgx, sgy) = (s.gx0 + i as i64, s.gy0 + j as i64);
                            let ub = o.budget.entry(((s.unit[src.k] as u64) << 8) | s.facies[src.k] as u64).or_default();
                            ub.litho = s.litho[src.k];
                            ub.facies = s.facies[src.k];
                            ub.sources += 1;
                            let rock_u = (src.rock_m3 / VU).round() as i64;
                            let weath_u = (src.rock_m3 * p.weathering_loss / VU).round() as i64;
                            ub.rock_u += rock_u;
                            ub.weathered_u += weath_u;
                            let rem = src.rock_m3 * (1.0 - p.weathering_loss);
                            for c in 0..nc {
                                if src.fr[c] <= 0.0 { continue; }
                                if p.is_block_class(c) { ub.block_class_rock_u += (rem * src.fr[c] / VU).round() as i64; continue; }
                                let bulk = rem * src.fr[c] / (1.0 - p.porosity);
                                let dmid = p.class_mid(c).min(src.rock.dmax_m);
                                // tracer ln(D / 1 mm) (positive for every class), 1e-3 units
                                let lnd_milli = ((dmid / 1.0e-3).ln() * 1000.0).round() as i64;
                                let cob = dmid >= p.cobble_lo_m && dmid < p.cobble_hi_m;
                                for m in 0..p.sub_parcels {
                                    let u = (bulk / p.sub_parcels as f64 / VU).round() as i64;
                                    if u <= 0 { continue; }
                                    let w = scope.draw("talus.field", pack(sgx, sgy), ((c as u32) << 8) | m).unwrap().u64();
                                    let fx0 = u16s(w, 0) - 0.5;
                                    let fy0 = u16s(w, 1) - 0.5;
                                    let kin = kin_of(dmid, &src.rock, w);
                                    let tr = trace(terr, i as i64, j as i64, fx0, fy0, &kin, cfg);
                                    o.steps += tr.steps as u64;
                                    ub.field_bulk_u += u;
                                    o.cls_parcels[c] += 1;
                                    o.cls_bulk[c] += u;
                                    match tr.stop {
                                        Stop::Rest | Stop::Capped => {
                                            if tr.stop == Stop::Capped { ub.field_cap_u += u; }
                                            ub.field_dep_u += u;
                                            for (kk, uu) in split4(nx, ny, i as i64, j as i64, tr.fx, tr.fy, u) {
                                                if uu == 0 { continue; }
                                                acc_v[kk].fetch_add(uu, Ordering::Relaxed);
                                                acc_m[kk].fetch_add(uu * lnd_milli, Ordering::Relaxed);
                                                if cob { acc_c[kk].fetch_add(uu, Ordering::Relaxed); }
                                            }
                                            o.cls_rest[c] += u;
                                            let dsrc = (tr.fx * tr.fx + tr.fy * tr.fy).sqrt() * dx;
                                            let dfoot = tr.foot.map(|(a, b)| ((tr.fx - a).powi(2) + (tr.fy - b).powi(2)).sqrt() * dx).unwrap_or(0.0);
                                            o.h_src[c][hbin(dsrc)] += u;
                                            o.h_foot[c][hbin(dfoot)] += u;
                                        }
                                        Stop::Sink => ub.field_water_u += u,
                                        Stop::OffDomain => ub.field_off_u += u,
                                    }
                                }
                            }
                        }
                        *slots[ci].lock().unwrap() = Some(o);
                        ci += threads;
                    }
                });
            }
        });
        slots.into_iter().map(|sl| sl.into_inner().unwrap().unwrap()).collect()
    };
    timing.insert("field_trace".into(), t1.elapsed().as_secs_f64());

    // 4. thickness, tracers, repose relaxation
    let t2 = Instant::now();
    let a2 = dx * dx;
    let mut thick: Vec<f64> = acc_v.iter().map(|a| a.load(Ordering::Relaxed) as f64 * VU / a2).collect();
    let mut tracers = vec![
        acc_m.iter().map(|a| a.load(Ordering::Relaxed) as f64 * VU * 1.0e-3 / a2).collect::<Vec<f64>>(),
        acc_c.iter().map(|a| a.load(Ordering::Relaxed) as f64 * VU / a2).collect::<Vec<f64>>(),
    ];
    let mut pre_lost = 0.0;
    for k in 0..n {
        if s.sink[k] != SINK_NONE && thick[k] > 0.0 {
            pre_lost += thick[k] * a2;
            thick[k] = 0.0;
            tracers[0][k] = 0.0;
            tracers[1][k] = 0.0;
        }
    }
    let tan_r = p.repose_deg.to_radians().tan();
    let mut levels = p.relax_levels.clone();
    levels.push((1, p.relax_iters));
    let tan_acc = (p.repose_deg + 5.0).to_radians().tan();
    let accept: Vec<bool> = rock_tan.iter().map(|&v| v <= tan_acc).collect();
    let stock_in: f64 = thick.iter().sum::<f64>() * a2;
    // keyed per-post repose jitter (world post index): breaks the axis-aligned straight toe fronts of the 8-neighbour
    // sandpile rule; stays inside the 33-37 deg band
    let tan_post: Vec<f64> = if p.repose_jitter_deg <= 0.0 { Vec::new() } else { (0..n).map(|k| {
        let (gi, gj) = (s.gx0 + (k % nx) as i64, s.gy0 + (k / nx) as i64);
        let u = scope.draw("talus.repose", pack(gi, gj), 0).map(|d| d.uniform_f64()).unwrap_or(0.5);
        (p.repose_deg + p.repose_jitter_deg * (2.0 * u - 1.0)).to_radians().tan()
    }).collect() };
    let ro = relax_multilevel(s.gx0, s.gy0, nx, ny, dx, &zr, &accept, if p.repose_jitter_deg > 0.0 { Some(&tan_post) } else { None }, &mut thick, &mut tracers, &s.sink, tan_r, &levels, p.relax_beta, threads);
    let stock_out: f64 = thick.iter().sum::<f64>() * a2;
    let relax_tp = crate::budget::relax_throughput(stock_in, stock_out, ro.lost_m3, ro.moved_m3);
    let tan_res = (p.repose_deg + 2.0).to_radians().tan();
    let (r_over, r_all) = residual(nx, ny, dx, &zr, &thick, tan_res, 0.1);
    let h: Vec<f64> = zr.iter().zip(thick.iter()).map(|(a, b)| a + b).collect();
    // talus grows through the exposure window: a block released at time fraction f runs over rock + f x talus.
    // Levels f_l = (l + 1/2) / L; each block takes a keyed level (uniform release time).
    let lv: Vec<f64> = (0..p.block_time_levels).map(|l| (l as f64 + 0.5) / p.block_time_levels as f64).collect();
    let lv_grad: Vec<(Vec<f64>, Vec<f64>, Vec<u8>)> = lv.iter().map(|&f| {
        let hl: Vec<f64> = zr.iter().zip(thick.iter()).map(|(a, b)| a + f * b).collect();
        let (gx, gy) = grad(nx, ny, dx, &hl);
        let gr: Vec<u8> = (0..n).map(|k| if f * thick[k] > p.debris_ground_m { 0 } else { ground0[k] }).collect();
        (gx, gy, gr)
    }).collect();
    timing.insert("relax".into(), t2.elapsed().as_secs_f64());

    // 5. blocks over the growing talus (time levels)
    let t3 = Instant::now();
    let terr_b: Vec<Terrain> = lv_grad.iter().map(|(gx, gy, gr)| Terrain { nx, ny, dx, gx, gy, sink: &s.sink, ground: gr, rock_tan: &rock_tan }).collect();
    let block_out: Vec<ChunkOut> = {
        let slots: Vec<std::sync::Mutex<Option<ChunkOut>>> = (0..n_chunks).map(|_| std::sync::Mutex::new(None)).collect();
        std::thread::scope(|sc| {
            for t in 0..threads {
                let (sources, slots, terrs, scope, cfg, thick, h, lv) = (&sources, &slots, &terr_b, &scope, &cfg, &thick, &h, &lv);
                let kin_of = &kin_of;
                sc.spawn(move || {
                    let mut ci = t;
                    while ci < n_chunks {
                        let mut o = ChunkOut::new(nc);
                        for src in &sources[ci * chunk..((ci + 1) * chunk).min(sources.len())] {
                            let (i, j) = (src.k % nx, src.k / nx);
                            let (sgx, sgy) = (s.gx0 + i as i64, s.gy0 + j as i64);
                            let ub = o.budget.entry(((s.unit[src.k] as u64) << 8) | s.facies[src.k] as u64).or_default();
                            ub.litho = s.litho[src.k];
                            ub.facies = s.facies[src.k];
                            let rem = src.rock_m3 * (1.0 - p.weathering_loss);
                            for c in 0..nc {
                                if src.fr[c] <= 0.0 || !p.is_block_class(c) { continue; }
                                let lo = p.class_edges_m[c];
                                let hi = p.class_edges_m[c + 1].min(src.rock.dmax_m);
                                if hi <= lo { continue; }
                                // expected block volume for D log-uniform in [lo, hi]: shape (hi^3 - lo^3) / (3 ln(hi / lo))
                                let vmid = p.block_shape * (hi.powi(3) - lo.powi(3)) / (3.0 * (hi / lo).ln());
                                let nexp = rem * src.fr[c] / vmid;
                                let wn = scope.draw("talus.block-n", pack(sgx, sgy), c as u32).unwrap().uniform_f64();
                                let nb = nexp.floor() as u32 + if wn < nexp - nexp.floor() { 1 } else { 0 };
                                for e in 0..nb {
                                    let att = ((c as u32) << 20) | e;
                                    let w = scope.draw("talus.block", pack(sgx, sgy), att).unwrap().u64();
                                    let w2 = scope.draw("talus.block-vis", pack(sgx, sgy), att).unwrap().uniform_f64();
                                    let id = scope.draw("talus.block-id", pack(sgx, sgy), att).unwrap().u64();
                                    let wsz = scope.draw("talus.block-size", pack(sgx, sgy), att).unwrap().uniform_f64();
                                    let d = (lo.ln() + wsz * (hi.ln() - lo.ln())).exp();
                                    let fx0 = u16s(w, 0) - 0.5;
                                    let fy0 = u16s(w, 1) - 0.5;
                                    let kin = kin_of(d, &src.rock, w);
                                    let wl = scope.draw("talus.block-time", pack(sgx, sgy), att).unwrap().below(lv.len() as u64).unwrap() as usize;
                                    let tr = trace(&terrs[wl], i as i64, j as i64, fx0, fy0, &kin, cfg);
                                    o.steps += tr.steps as u64;
                                    let vol = p.block_shape * d.powi(3);
                                    let vu = (vol / VU).round() as i64;
                                    ub.block_rock_u += vu;
                                    ub.blocks += 1;
                                    o.cls_parcels[c] += 1;
                                    o.cls_bulk[c] += vu;
                                    let (ex, ey) = ((sgx as f64 + tr.fx) * dx, (sgy as f64 + tr.fy) * dx);
                                    let dsrc = (tr.fx * tr.fx + tr.fy * tr.fy).sqrt() * dx;
                                    let dfoot = tr.foot.map(|(a, b)| ((tr.fx - a).powi(2) + (tr.fy - b).powi(2)).sqrt() * dx).unwrap_or(0.0);
                                    let inside = tr.stop != Stop::OffDomain;
                                    let (zz, burial) = if inside { (bil(nx, ny, h, i as i64, j as i64, tr.fx, tr.fy), (1.0 - lv[wl]) * bil(nx, ny, thick, i as i64, j as i64, tr.fx, tr.fy)) } else { (f64::NAN, 0.0) };
                                    let mut visible = false;
                                    match tr.stop {
                                        Stop::Rest | Stop::Capped => {
                                            if tr.stop == Stop::Rest { ub.blocks_rest += 1 } else { ub.blocks_cap += 1 }
                                            let pv = if burial > 1e-9 { (p.vis_k * d / burial).min(1.0) } else { 1.0 };
                                            visible = w2 < pv;
                                            o.cls_rest[c] += vu;
                                            o.h_src[c][hbin(dsrc)] += vu;
                                            o.h_foot[c][hbin(dfoot)] += vu;
                                        }
                                        Stop::Sink => {
                                            ub.blocks_sink += 1;
                                            visible = tr.sink != SINK_WATER;
                                        }
                                        Stop::OffDomain => ub.blocks_off += 1,
                                    }
                                    if visible { ub.blocks_visible += 1; }
                                    if inside {
                                        o.blocks.push(Block {
                                            id, x: ex, y: ey, z: zz, d_m: d, src_gx: sgx, src_gy: sgy, class: c as u8, event: e,
                                            litho: s.litho[src.k], facies: s.facies[src.k], unit: s.unit[src.k],
                                            d_src_m: dsrc, d_foot_m: dfoot, path_m: tr.path_m, burial_m: burial, visible,
                                            stop: tr.stop as u8, sink: tr.sink,
                                        });
                                    }
                                }
                            }
                        }
                        *slots[ci].lock().unwrap() = Some(o);
                        ci += threads;
                    }
                });
            }
        });
        slots.into_iter().map(|sl| sl.into_inner().unwrap().unwrap()).collect()
    };
    timing.insert("block_trace".into(), t3.elapsed().as_secs_f64());

    // 6. merge (chunk order) and cut the core
    let mut budget: BTreeMap<u64, UnitBudget> = BTreeMap::new();
    let mut classes: Vec<ClassStats> = (0..nc)
        .map(|c| ClassStats { lo_m: p.class_edges_m[c], hi_m: p.class_edges_m[c + 1], block_class: p.is_block_class(c), parcels: 0, bulk_u: 0, rest_u: 0, h_foot: vec![0; HBINS], h_src: vec![0; HBINS] })
        .collect();
    let (mut fsteps, mut bsteps) = (0u64, 0u64);
    let mut blocks_all: Vec<Block> = Vec::new();
    for (which, outs) in [(0, field_out), (1, block_out)] {
        for o in outs {
            if which == 0 { fsteps += o.steps } else { bsteps += o.steps }
            for (u, b) in &o.budget {
                let e = budget.entry(*u).or_default();
                e.litho = b.litho;
                e.facies = b.facies;
                e.add(b);
            }
            for c in 0..nc {
                classes[c].parcels += o.cls_parcels[c];
                classes[c].bulk_u += o.cls_bulk[c];
                classes[c].rest_u += o.cls_rest[c];
                for b in 0..HBINS { classes[c].h_foot[b] += o.h_foot[c][b]; classes[c].h_src[b] += o.h_src[c][b]; }
            }
            blocks_all.extend(o.blocks);
        }
    }
    // a block belongs to the tile whose cells [gx0, gx1) x [gy0, gy1) (chart m) contain it: no duplicates on shared edges
    let (ox0, ox1, oy0, oy1) = (core.gx0 as f64 * dx, core.gx1 as f64 * dx, core.gy0 as f64 * dx, core.gy1 as f64 * dx);
    let mut blocks: Vec<Block> = blocks_all.into_iter().filter(|b| b.x >= ox0 && b.x < ox1 && b.y >= oy0 && b.y < oy1).collect();
    blocks.sort_by(|a, b| (a.src_gy, a.src_gx, a.class, a.event).cmp(&(b.src_gy, b.src_gx, b.class, b.event)));

    let ci0 = (core.gx0 - s.gx0) as usize;
    let cj0 = (core.gy0 - s.gy0) as usize;
    let cnx = (core.gx1 - core.gx0 + 1) as usize;
    let cny = (core.gy1 - core.gy0 + 1) as usize;
    if core.gx0 < s.gx0 || core.gy0 < s.gy0 || ci0 + cnx > nx || cj0 + cny > ny { return Err("core outside the domain".into()); }
    let cn = cnx * cny;
    let mut out = RunOut {
        core, nx: cnx, ny: cny,
        thickness_mm: Vec::with_capacity(cn), chute_mm: Vec::with_capacity(cn), source_mm: Vec::with_capacity(cn),
        d50_mm: Vec::with_capacity(cn), cover_pm: Vec::with_capacity(cn), cobble_per_100m2: Vec::with_capacity(cn),
        blocks, budget, classes,
        stats: RunStats {
            domain: Rect { gx0: s.gx0, gy0: s.gy0, gx1: s.gx0 + nx as i64 - 1, gy1: s.gy0 + ny as i64 - 1 }, core, dx,
            sources: sources.len() as u64, chutes: n_chutes, chute_posts: chute.iter().filter(|&&c| c > 0.0).count() as u64,
            relax_lost_m3: ro.lost_m3 + pre_lost, relax_throughput: relax_tp, relax_iters: ro.iters, repose_residual_over: r_over, repose_residual_posts: r_all,
            field_traces_steps: fsteps, block_traces_steps: bsteps, timing_s: BTreeMap::new(),
        },
    };
    let cob_mid = (p.cobble_lo_m * p.cobble_hi_m).sqrt();
    for j in cj0..cj0 + cny {
        for i in ci0..ci0 + cnx {
            let k = j * nx + i;
            let t = thick[k];
            out.thickness_mm.push((t * 1000.0).round() as i32);
            out.chute_mm.push((chute[k] * 1000.0).round() as i32);
            out.source_mm.push((source_depth[k] * 1000.0).round() as i32);
            let d50 = if t > 1e-3 { 1.0e-3 * (tracers[0][k] / t).exp() } else { 0.0 };
            out.d50_mm.push((d50 * 1000.0).round().clamp(0.0, 65535.0) as u16);
            let cover = if t > 0.0 && d50 > 0.0 { (t / (p.cover_k * d50)).min(1.0) } else { 0.0 };
            out.cover_pm.push((cover * 1000.0).round() as u16);
            let fc = if t > 1e-3 { (tracers[1][k] / t).clamp(0.0, 1.0) } else { 0.0 };
            // visible cobbles per 100 m2: covered area share x cobble volume share / clast footprint (pi/4 D^2)
            let dens = 100.0 * cover * fc / (std::f64::consts::FRAC_PI_4 * cob_mid * cob_mid);
            out.cobble_per_100m2.push(dens.round().clamp(0.0, 65535.0) as u16);
        }
    }
    timing.insert("total".into(), t_all.elapsed().as_secs_f64());
    out.stats.timing_s = timing;
    Ok(out)
}
