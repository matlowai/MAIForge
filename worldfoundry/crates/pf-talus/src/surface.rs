//! Input surface: a regular grid of posts on the global chart lattice (post (i, j) of a surface sits at chart
//! ((gx0 + i) dx, (gy0 + j) dx), row 0 south), with the fields the talus process reads. Loaders for a pf-detail
//! 2 m block (+ optional VOLCANIC-SPLIT facies raster, + optional T1b window for fluvial drainage) and for a T1b
//! window alone (32 m posts).

use crate::R;
use serde_json::Value;
use std::path::Path;

/// Sink codes: what happens to a clast that stops on the post.
pub const SINK_NONE: u8 = 0;
pub const SINK_CHANNEL: u8 = 1;
pub const SINK_WATER: u8 = 2;
pub const SINK_FLUVIAL: u8 = 3;

#[derive(Clone)]
pub struct Surface {
    pub dx: f64,
    pub gx0: i64,
    pub gy0: i64,
    pub nx: usize,
    pub ny: usize,
    /// Height, m.
    pub z: Vec<f64>,
    /// LithologyCodeV1; 255 water / none.
    pub litho: Vec<u8>,
    /// VOLCANIC-SPLIT facies (1 thin lava, 2 ash, 3 riser lava), 0 none.
    pub facies: Vec<u8>,
    /// ct regime (0 soil, 1 talus, 2 cliff, 3 tread), 255 unknown.
    pub regime: Vec<u8>,
    pub sink: Vec<u8>,
    /// Source unit id for the budget (exposed-unit provenance id when published, else litho << 8 | facies).
    pub unit: Vec<u32>,
    /// Input digests (name, blake3) for the recipe.
    pub digests: Vec<(String, String)>,
}

fn rd(p: &Path) -> R<Vec<u8>> { std::fs::read(p).map_err(|e| format!("{}: {e}", p.display())) }
pub fn blake(b: &[u8]) -> String { format!("blake3:{}", blake3::hash(b).to_hex()) }
fn i32s(b: &[u8]) -> Vec<i32> { b.chunks_exact(4).map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
fn u32s(b: &[u8]) -> Vec<u32> { b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
fn f32s(b: &[u8]) -> Vec<f32> { b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }

impl Surface {
    pub fn n(&self) -> usize { self.nx * self.ny }
    pub fn idx(&self, i: usize, j: usize) -> usize { j * self.nx + i }
    pub fn x_of(&self, i: usize) -> f64 { (self.gx0 + i as i64) as f64 * self.dx }
    pub fn y_of(&self, j: usize) -> f64 { (self.gy0 + j as i64) as f64 * self.dx }

    /// Sub-surface over global post rectangle [gx0, gx1] x [gy0, gy1] (inclusive), clipped to this surface.
    pub fn crop(&self, gx0: i64, gy0: i64, gx1: i64, gy1: i64) -> Surface {
        let a0 = gx0.max(self.gx0);
        let b0 = gy0.max(self.gy0);
        let a1 = gx1.min(self.gx0 + self.nx as i64 - 1);
        let b1 = gy1.min(self.gy0 + self.ny as i64 - 1);
        let (nx, ny) = ((a1 - a0 + 1) as usize, (b1 - b0 + 1) as usize);
        let (oi, oj) = ((a0 - self.gx0) as usize, (b0 - self.gy0) as usize);
        let mut s = Surface {
            dx: self.dx, gx0: a0, gy0: b0, nx, ny,
            z: Vec::with_capacity(nx * ny), litho: Vec::with_capacity(nx * ny), facies: Vec::with_capacity(nx * ny),
            regime: Vec::with_capacity(nx * ny), sink: Vec::with_capacity(nx * ny), unit: Vec::with_capacity(nx * ny),
            digests: self.digests.clone(),
        };
        for j in 0..ny {
            let r0 = (oj + j) * self.nx + oi;
            s.z.extend_from_slice(&self.z[r0..r0 + nx]);
            s.litho.extend_from_slice(&self.litho[r0..r0 + nx]);
            s.facies.extend_from_slice(&self.facies[r0..r0 + nx]);
            s.regime.extend_from_slice(&self.regime[r0..r0 + nx]);
            s.sink.extend_from_slice(&self.sink[r0..r0 + nx]);
            s.unit.extend_from_slice(&self.unit[r0..r0 + nx]);
        }
        s
    }
}

/// A pf-detail block (`pf-detail block` output dir: block.json + block_*.{i32,u8,u32}, row-major, row 0 south,
/// edge-inclusive (nx 128 + 1) posts at 2 m). `x0_m, y0_m` = chart position of post (0, 0) (tile ix0 x 256 m).
pub fn load_block(dir: &Path, x0_m: f64, y0_m: f64, facies: Option<&Path>, window: Option<&Path>, fluvial_area_m2: f64) -> R<Surface> {
    let hb = rd(&dir.join("block.json"))?;
    let h: Value = serde_json::from_slice(&hb).map_err(|e| format!("block.json: {e}"))?;
    let sz = h["block_size_posts"].as_array().ok_or("block_size_posts")?;
    let (nx, ny) = (sz[0].as_u64().ok_or("nx")? as usize, sz[1].as_u64().ok_or("ny")? as usize);
    let n = nx * ny;
    let dx = 2.0;
    let gx0 = (x0_m / dx).round() as i64;
    let gy0 = (y0_m / dx).round() as i64;
    if ((gx0 as f64) * dx - x0_m).abs() > 1e-9 || ((gy0 as f64) * dx - y0_m).abs() > 1e-9 { return Err("block origin not on the 2 m lattice".into()); }
    let mut digests = vec![("block.json".to_string(), blake(&hb))];
    let zb = rd(&dir.join("block_height_mm.i32"))?;
    digests.push(("block_height_mm".into(), blake(&zb)));
    let z: Vec<f64> = i32s(&zb).into_iter().map(|v| v as f64 / 1000.0).collect();
    if z.len() != n { return Err(format!("height has {} posts, want {n}", z.len())); }
    let litho = rd(&dir.join("block_litho.u8"))?;
    digests.push(("block_litho".into(), blake(&litho)));
    let regime = match std::fs::read(dir.join("block_ct_regime.u8")) {
        Ok(b) => { digests.push(("block_ct_regime".into(), blake(&b))); b }
        Err(_) => vec![255u8; n],
    };
    let fixed = rd(&dir.join("block_fixed.u8"))?;
    digests.push(("block_fixed".into(), blake(&fixed)));
    let fac = match facies {
        Some(p) => { let b = rd(p)?; digests.push(("facies".into(), blake(&b))); if b.len() != n { return Err("facies size".into()); } b }
        None => vec![0u8; n],
    };
    let prov = match std::fs::read(dir.join("block_prov_id.u32")) {
        Ok(b) => { digests.push(("block_prov_id".into(), blake(&b))); Some(u32s(&b)) }
        Err(_) => None,
    };
    let mut sink: Vec<u8> = fixed.iter().zip(litho.iter()).map(|(&f, &l)| match f {
        1 | 5 => SINK_CHANNEL,
        2..=4 => SINK_WATER,
        _ => if l == 255 { SINK_WATER } else { SINK_NONE },
    }).collect();
    if let Some(wd) = window {
        let w = Window::load(wd)?;
        digests.push(("t1b.window".into(), w.digest.clone()));
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                if sink[k] != SINK_NONE { continue; }
                let a = w.bilinear_area((gx0 + i as i64) as f64 * dx, (gy0 + j as i64) as f64 * dx);
                if a > fluvial_area_m2 { sink[k] = SINK_FLUVIAL; }
            }
        }
    }
    let unit = (0..n).map(|k| match &prov { Some(p) => p[k], None => ((litho[k] as u32) << 8) | fac[k] as u32 }).collect();
    Ok(Surface { dx, gx0, gy0, nx, ny, z, litho, facies: fac, regime, sink, unit, digests })
}

/// The subset of a published T1b window (pf-region t1b_window layout) the talus process reads.
pub struct Window {
    pub digest: String,
    pub gx0: i64,
    pub gy0: i64,
    pub n: usize,
    pub height_mm: Vec<i32>,
    pub fixed: Vec<u8>,
    pub drainage: Vec<f32>,
    pub litho: Option<Vec<u8>>,
    pub regime: Option<Vec<u8>>,
    pub prov: Option<Vec<u32>>,
}

impl Window {
    pub fn load(dir: &Path) -> R<Window> {
        let hb = rd(&dir.join("window.json"))?;
        let v: Value = serde_json::from_slice(&hb).map_err(|e| format!("window.json: {e}"))?;
        let h = &v["header"];
        let lat = &h["lattice"];
        let (gx0, gy0, n) = (lat["gx0"].as_i64().ok_or("gx0")?, lat["gy0"].as_i64().ok_or("gy0")?, lat["n"].as_u64().ok_or("n")? as usize);
        let file = |name: &str| -> Option<std::path::PathBuf> {
            h["fields"].as_array()?.iter().find(|f| f["name"] == name).and_then(|f| f["file"].as_str()).map(|s| dir.join(s))
        };
        let req = |name: &str| file(name).ok_or(format!("window field {name} missing"));
        let height_mm = i32s(&rd(&req("height_mm")?)?);
        let fixed = rd(&req("fixed")?)?;
        let drainage = f32s(&rd(&req("drainage_area_m2")?)?);
        let litho = file("exposed_litho").map(|p| rd(&p)).transpose()?;
        let regime = file("ct_regime").map(|p| rd(&p)).transpose()?;
        let prov = file("prov_id").map(|p| rd(&p).map(|b| u32s(&b))).transpose()?;
        let digest = v["digest"].as_str().unwrap_or("").to_string();
        Ok(Window { digest, gx0, gy0, n, height_mm, fixed, drainage, litho, regime, prov })
    }

    /// Drainage area bilinear between 32 m posts at chart (x, y) m; 0 outside.
    pub fn bilinear_area(&self, x: f64, y: f64) -> f64 {
        let fx = x / 32.0 - self.gx0 as f64;
        let fy = y / 32.0 - self.gy0 as f64;
        if fx < 0.0 || fy < 0.0 || fx >= (self.n - 1) as f64 || fy >= (self.n - 1) as f64 { return 0.0; }
        let (i, j) = (fx.floor() as usize, fy.floor() as usize);
        let (a, b) = (fx - i as f64, fy - j as f64);
        let at = |i: usize, j: usize| self.drainage[j * self.n + i] as f64;
        (1.0 - a) * (1.0 - b) * at(i, j) + a * (1.0 - b) * at(i + 1, j) + (1.0 - a) * b * at(i, j + 1) + a * b * at(i + 1, j + 1)
    }
}

/// T1b window alone (32 m posts) over global 32 m post rectangle [gx0, gx1] x [gy0, gy1] (inclusive).
pub fn load_window_surface(dir: &Path, gx0: i64, gy0: i64, gx1: i64, gy1: i64, fluvial_area_m2: f64) -> R<Surface> {
    let w = Window::load(dir)?;
    let dx = 32.0;
    let a0 = gx0.max(w.gx0);
    let b0 = gy0.max(w.gy0);
    let a1 = gx1.min(w.gx0 + w.n as i64 - 1);
    let b1 = gy1.min(w.gy0 + w.n as i64 - 1);
    let (nx, ny) = ((a1 - a0 + 1) as usize, (b1 - b0 + 1) as usize);
    let mut s = Surface { dx, gx0: a0, gy0: b0, nx, ny, z: vec![], litho: vec![], facies: vec![0; nx * ny], regime: vec![], sink: vec![], unit: vec![], digests: vec![("t1b.window".into(), w.digest.clone())] };
    for j in 0..ny {
        for i in 0..nx {
            let k = ((b0 - w.gy0) as usize + j) * w.n + (a0 - w.gx0) as usize + i;
            let l = w.litho.as_ref().map(|v| v[k]).unwrap_or(255);
            s.z.push(w.height_mm[k] as f64 / 1000.0);
            s.litho.push(l);
            s.regime.push(w.regime.as_ref().map(|v| v[k]).unwrap_or(255));
            s.sink.push(match w.fixed[k] { 1 | 5 => SINK_CHANNEL, 2 | 4 => SINK_WATER, _ => if w.drainage[k] as f64 > fluvial_area_m2 { SINK_FLUVIAL } else if l == 255 { SINK_WATER } else { SINK_NONE } });
            s.unit.push(w.prov.as_ref().map(|v| v[k]).unwrap_or((l as u32) << 8));
        }
    }
    Ok(s)
}
