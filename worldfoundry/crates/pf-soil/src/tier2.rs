//! The 2 m tier: zonation on a pf-detail block. Exact at 2 m: height (temperature lapse, snow line, tree line), the
//! chart position (climate), the stability slope (3 x 3 central differences), the 2 m regime (cliff / talus / soil /
//! tread), carved channels (water) and banks (riparian), lake and sea fringes from the 2 m height, the vent and tephra kill
//! zones, pf-talus debris cover. Sampled bilinearly from the 32 m tier (same window): topographic wetness, the volcanic
//! colonisation factor and the soil-age factor (the published 2 m block carries no VE unit id), the beach code, the
//! shore distance. Insolation uses a +-`insol_baseline_2m_posts` baseline (the canopy-scale aspect, not 2 m facets).
//! Posts within `ring_2m_posts` of the block edge are no-data: a consumer runs a block one ring larger than it needs.

use crate::climate::Climate;
use crate::params::SoilParams;
use crate::tier32::{climate_at, kill_factor, par_map, Tier32};
use crate::window::Window;
use crate::zonation::{insolation_ratio, post_t_gs, zone, PostIn, PostOut};
use serde_json::Value;
use std::path::Path;

pub const POST2_M: f64 = 2.0;

pub struct Block {
    pub n_x: usize,
    pub n_y: usize,
    pub x0_m: f64,
    pub y0_m: f64,
    pub height_m: Vec<f64>,
    pub fixed: Vec<u8>,
    pub regime: Vec<u8>,
    pub digests: Vec<(String, String)>,
}

fn rd(dir: &Path, name: &str) -> Result<(Vec<u8>, String), String> {
    let b = std::fs::read(dir.join(name)).map_err(|e| format!("{}/{name}: {e}", dir.display()))?;
    let d = format!("blake3:{}", blake3::hash(&b).to_hex());
    Ok((b, d))
}

impl Block {
    /// Load a pf-detail block (block.json + block_*.{i32,u8}); origin from the lod-0 tile list (256 m tiles) unless given.
    /// `crop` = (i0, j0, ni, nj) in block posts: read only that sub-array (seam tests).
    pub fn load(dir: &Path, origin: Option<(f64, f64)>, crop: Option<(usize, usize, usize, usize)>) -> Result<Block, String> {
        let bj: Value = serde_json::from_slice(&std::fs::read(dir.join("block.json")).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let sz = bj["block_size_posts"].as_array().ok_or("block_size_posts")?;
        let (nx, ny) = (sz[0].as_u64().ok_or("nx")? as usize, sz[1].as_u64().ok_or("ny")? as usize);
        let (x0, y0) = match origin {
            Some(o) => o,
            None => {
                let t = bj["tiles"].as_array().ok_or("block tiles")?;
                let tx = t.iter().filter_map(|a| a[0].as_i64()).min().ok_or("tiles")?;
                let ty = t.iter().filter_map(|a| a[1].as_i64()).min().ok_or("tiles")?;
                (tx as f64 * 256.0, ty as f64 * 256.0)
            }
        };
        let mut digests = Vec::new();
        let (hb, hd) = rd(dir, "block_height_mm.i32")?; digests.push(("block_height_mm.i32".to_string(), hd));
        let (fb, fd) = rd(dir, "block_fixed.u8")?; digests.push(("block_fixed.u8".to_string(), fd));
        let (rb, rdg) = rd(dir, "block_ct_regime.u8")?; digests.push(("block_ct_regime.u8".to_string(), rdg));
        if hb.len() != nx * ny * 4 || fb.len() != nx * ny || rb.len() != nx * ny { return Err("block field sizes".into()); }
        let (ci, cj, cn, cm) = crop.unwrap_or((0, 0, nx, ny));
        if ci + cn > nx || cj + cm > ny { return Err("crop outside the block".into()); }
        let mut height_m = Vec::with_capacity(cn * cm);
        let mut fixed = Vec::with_capacity(cn * cm);
        let mut regime = Vec::with_capacity(cn * cm);
        for j in cj..cj + cm { for i in ci..ci + cn {
            let k = j * nx + i;
            height_m.push(i32::from_le_bytes(hb[4 * k..4 * k + 4].try_into().unwrap()) as f64 / 1000.0);
            fixed.push(fb[k]);
            regime.push(rb[k]);
        } }
        Ok(Block { n_x: cn, n_y: cm, x0_m: x0 + POST2_M * ci as f64, y0_m: y0 + POST2_M * cj as f64, height_m, fixed, regime, digests })
    }
}

/// pf-talus debris cover on its own 2 m core lattice.
pub struct Talus { pub n: usize, pub x0_m: f64, pub y0_m: f64, pub cover: Vec<f64>, pub digest: String }

impl Talus {
    pub fn load(dir: &Path) -> Result<Talus, String> {
        let tj: Value = serde_json::from_slice(&std::fs::read(dir.join("talus.json")).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let c = &tj["recipe"]["core"];
        let (gx0, gx1, gy0, gy1) = (c["gx0"].as_i64().ok_or("core")?, c["gx1"].as_i64().ok_or("core")?, c["gy0"].as_i64().ok_or("core")?, c["gy1"].as_i64().ok_or("core")?);
        if gx1 - gx0 != gy1 - gy0 { return Err("talus core not square".into()); }
        let n = (gx1 - gx0 + 1) as usize;
        let (b, d) = rd(dir, "fields/debris_cover_permille.u16")?;
        if b.len() != n * n * 2 { return Err(format!("talus debris_cover size {} != {}", b.len(), n * n * 2)); }
        let cover = b.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]]) as f64 / 1000.0).collect();
        Ok(Talus { n, x0_m: gx0 as f64 * POST2_M, y0_m: gy0 as f64 * POST2_M, cover, digest: d })
    }
    pub fn at(&self, x: f64, y: f64) -> Option<f64> {
        let (i, j) = (((x - self.x0_m) / POST2_M).round(), ((y - self.y0_m) / POST2_M).round());
        if i < 0.0 || j < 0.0 || i >= self.n as f64 || j >= self.n as f64 { return None; }
        Some(self.cover[j as usize * self.n + i as usize])
    }
}

pub struct Tier2 {
    pub n_x: usize,
    pub n_y: usize,
    pub x0_m: f64,
    pub y0_m: f64,
    pub valid: Vec<bool>,
    pub twi: Vec<f64>,
    pub r_sun: Vec<f64>,
    pub out: Vec<PostOut>,
    pub talus_posts: usize,
}

pub fn run(p: &SoilParams, w: &Window, t32: &Tier32, cl: &Climate, b: &Block, talus: Option<&Talus>, threads: usize) -> Result<Tier2, String> {
    let (nx, ny) = (b.n_x, b.n_y);
    let nn = nx * ny;
    let ring = p.ring_2m_posts.max(1.0) as usize;
    let bl = p.insol_baseline_2m_posts.max(1.0) as usize;
    let rip = (p.riparian_2m_m / POST2_M).ceil() as usize;
    let rs = p.substrate_smooth_2m_posts.max(0.0) as usize;
    let rd = p.debris_smooth_2m_posts.max(0.0) as usize;
    if bl > ring || rip > ring || rs > ring || rd > ring { return Err(format!("2 m stencils (insolation {bl}, riparian {rip}, substrate {rs}, debris {rd}) exceed the ring {ring}")); }
    // per-post regime factor and debris cover, then box means (fixed summation order: deterministic)
    let regf: Vec<f64> = b.regime.iter().map(|&r| crate::zonation::regime_factor(p, r)).collect();
    let deb: Vec<f64> = (0..nn).map(|k| { let (i, j) = (k % nx, k / nx); talus.and_then(|t| t.at(b.x0_m + POST2_M * i as f64, b.y0_m + POST2_M * j as f64)).unwrap_or(0.0) }).collect();
    let boxm = |f: &[f64], k: usize, r: usize| -> f64 {
        let (i, j) = ((k % nx) as i64, (k / nx) as i64);
        let r = r as i64;
        let mut s = 0.0;
        for dj in -r..=r { for di in -r..=r { s += f[((j + dj) as usize) * nx + (i + di) as usize]; } }
        s / ((2 * r + 1) * (2 * r + 1)) as f64
    };
    let chart = crate::tier32::chart_of(w)?;
    let z = &b.height_m;
    let valid: Vec<bool> = (0..nn).map(|k| { let (i, j) = (k % nx, k / nx); i >= ring && j >= ring && i + ring < nx && j + ring < ny }).collect();
    let res: Vec<Option<(f64, f64, PostOut, bool)>> = par_map(nn, threads, |k| {
        if !valid[k] { return None; }
        let (i, j) = (k % nx, k / nx);
        let x = b.x0_m + POST2_M * i as f64;
        let y = b.y0_m + POST2_M * j as f64;
        let xm = (x * 1000.0).round() as i64;
        let ym = (y * 1000.0).round() as i64;
        let (la, t_sl, half, pp, ai) = match climate_at(cl, &chart, xm, ym) { Ok(c) => c, Err(_) => return Some((f64::NAN, f64::NAN, PostOut::default(), false)) };
        // stability slope (3 x 3 central differences) and canopy-scale gradient (+-bl posts)
        let sx = (z[k + 1] - z[k - 1]) / (2.0 * POST2_M);
        let sy = (z[k + nx] - z[k - nx]) / (2.0 * POST2_M);
        let slope = sx.hypot(sy).atan().to_degrees();
        let gx = (z[k + bl] - z[k - bl]) / (2.0 * POST2_M * bl as f64);
        let gy = (z[k + bl * nx] - z[k - bl * nx]) / (2.0 * POST2_M * bl as f64);
        let flat = gx.hypot(gy).atan().to_degrees();
        let r_sun = insolation_ratio(p, la, gx, gy);
        let s = |f: &[f64]| t32.sample(f, x, y);
        let (Some(twi), Some(fa), Some(sa), Some(bch)) = (s(&t32.twi), s(&t32.f_age), s(&t32.soil_age_f), s(&t32.beach)) else { return Some((f64::NAN, f64::NAN, PostOut::default(), false)); };
        let shore = match &t32.shore_m { Some(sf) => s(sf), None => None };
        let stage = t32.corner_max(&t32.lake_stage, x, y);
        let fx = b.fixed[k];
        let h = z[k];
        let water = fx == 1 || fx == 3;
        let mut riparian = fx == 2;
        if !water && !riparian {
            'o: for dj in -(rip as i64)..=(rip as i64) { for di in -(rip as i64)..=(rip as i64) {
                if ((di * di + dj * dj) as f64).sqrt() * POST2_M > p.riparian_2m_m { continue; }
                let q = ((j as i64 + dj) as usize) * nx + (i as i64 + di) as usize;
                if b.fixed[q] == 1 { riparian = true; break 'o; }
            } }
        }
        let lake_fringe = stage.is_finite() && !water && h - stage > 0.0 && h - stage <= p.lake_marsh_m;
        let debris = talus.and_then(|t| t.at(x, y));
        let mut xin = PostIn { h, sea: cl.sea, water, t_sl, half, p: pp, ai, lapse: cl.lapse, r_sun, slope_deg: slope, regime: b.regime[k],
            reg_factor: Some(boxm(&regf, k, rs)), twi, f_age: fa, f_kill: 1.0, soil_age_f: sa, debris: boxm(&deb, k, rd), beach: bch >= 0.5, riparian,
            lake_fringe, flat_slope_deg: flat, shore_m: shore };
        let t_gs = post_t_gs(p, &xin);
        xin.f_kill = kill_factor(p, &w.volcanoes, x, y, t_gs, ai);
        Some((twi, r_sun, zone(p, &xin), debris.is_some()))
    });
    let mut twi = vec![f64::NAN; nn];
    let mut r_sun = vec![f64::NAN; nn];
    let mut out = vec![PostOut { class: crate::zonation::C_NODATA, ..Default::default() }; nn];
    let mut vv = valid.clone();
    let mut talus_posts = 0;
    for (k, r) in res.into_iter().enumerate() {
        if let Some((t, rs, o, tp)) = r {
            if !t.is_finite() { return Err(format!("2 m post {k}: climate or 32 m sample unavailable (window does not cover the block plus one 32 m post)")); }
            twi[k] = t; r_sun[k] = rs; out[k] = o; if tp { talus_posts += 1; }
        } else { vv[k] = false; }
    }
    Ok(Tier2 { n_x: nx, n_y: ny, x0_m: b.x0_m, y0_m: b.y0_m, valid: vv, twi, r_sun, out, talus_posts })
}
