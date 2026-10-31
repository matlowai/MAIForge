//! Upstream inputs: a published T1b window, the T1a skeleton, T0 terrain fields and
//! the full stratigraphy column dumps (lane SA format). All read-only; every blob is
//! digested so the T2 recipe can name its upstream bytes.

use crate::R;
use serde_json::Value;
use std::path::{Path, PathBuf};

fn rd(p: &Path) -> R<Vec<u8>> { std::fs::read(p).map_err(|e| format!("{}: {e}", p.display())) }
pub fn blake(b: &[u8]) -> String { format!("blake3:{}", blake3::hash(b).to_hex()) }
fn f32s(b: &[u8]) -> Vec<f32> { b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
fn i32s(b: &[u8]) -> Vec<i32> { b.chunks_exact(4).map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
fn u32s(b: &[u8]) -> Vec<u32> { b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }

/// A published T1b window (pf-region `t1b_window::publish` layout), verified.
pub struct T1bWindow {
    pub escarpment: Option<crate::escarpment::FineSurface>,
    pub evolved: Option<crate::contrast::Evolved>,
    pub key: String,
    pub digest: String,
    pub recipe_digest: String,
    pub generator_version: String,
    pub gx0: i64,
    pub gy0: i64,
    pub n: usize,
    pub height_m: Vec<f64>,
    pub fixed: Vec<u8>,
    pub drainage_area_m2: Vec<f32>,
    pub t0_cell: Vec<u32>,
    /// v2 (lane T2B): the eikonal surface before crest rounding, fill and the coast/volcano construct (edifice test).
    pub base_m: Option<Vec<f64>>,
    /// v2 (lane T2B): T1b exposed lithology with the over-codes (3 edifice, 10 fan/beach), 255 water.
    pub exposed_litho: Option<Vec<u8>>,
    /// Lane PV (read only when present; nothing else reads them): T1b over code, volcano index and the header's volcano list.
    pub over: Option<Vec<u8>>,
    pub volcano_id: Option<Vec<u8>>,
    pub volcanoes: Vec<crate::strata_prov::VolcInfo>,
    /// Lane BS (read only by beds=2): per listed volcano (vent x, vent y chart m, r_influence_m, cone height h_m); NaN = absent.
    pub volcano_geo: Vec<[f64; 4]>,
    /// Lane DE (read only by de_lobes=1): per listed volcano, T1b's recorded lava-flow lobes (cf.flows: angle rad, r0 m, length m,
    /// half-width m, thickness m, wander amplitude m, wander length m, wander phase rad); empty when absent.
    pub volcano_flows: Vec<Vec<[f64; 8]>>,
    /// Lane DE: lane VE's per-post edifice_style (u8) when the window publishes it (read only by de_lobes=1).
    pub edifice_style: Option<Vec<u8>>,
    /// Lane DE2 (read only by de_ve=1): per listed volcano, lane VE's unit list (header prep.net.volcanoes[k].ve.units) when present.
    pub volcano_ve_units: Vec<Option<Value>>,
    /// Lane PYRO (read only by pyro_sub=1): per listed volcano, lane VE's event counts (header ve.counts) when present, and the
    /// main cone radius r_main_m when recorded.
    pub volcano_ve_counts: Vec<Option<Value>>,
    pub volcano_r_body: Vec<Option<f64>>,
    /// Lane VH-I (read only by pyro_wind=1): per listed volcano, the T0 downwind azimuth (header ve.wind.chart_azimuth_rad).
    pub volcano_ve_wind: Vec<Option<f64>>,
    /// Lane DATUM (read only by ed_datum=1): T1b's surface before the wave cut (lane TH uncut_mm) where published, else the
    /// published height; None when the window has no uncut_mm field.
    pub pre_cut_m: Option<Vec<f64>>,
    /// Lane SHARP (read only by ct_sharp=1): T1b's cliff-and-talus fields when the window publishes them (pf-region ct=1 with
    /// ct_base=1); None otherwise.
    pub ct: Option<CtWindow>,
    /// Lane CLIMATE-RUNOFF (read only by climate_runoff=1): the window's runoff_weight (R / r_ref per T1b post) when published.
    pub runoff_weight: Option<Vec<f32>>,
}

/// Lane SHARP: the published T1b ct fields and header params of a window.
pub struct CtWindow {
    /// ct_base_mm (m): the surface without the staircase shift.
    pub base_m: Vec<f64>,
    pub regime: Vec<u8>,
    pub cliff_permille: Vec<u16>,
    /// Header params.ls.tb.tc.ct (class table, step period) and the rivernet rockmass bounds (min, max multiplier).
    pub params: Value,
    pub rockmass: (f64, f64),
    /// Lane SHARP2: pf-region ct_tread's 3-limb fields (ct_tread_permille ft, ct_tread_tan_permille tt) when the window
    /// publishes them; None otherwise (the 2-limb re-profile, every byte as lane SHARP / TREAD).
    pub tread: Option<(Vec<u16>, Vec<u16>)>,
    /// Lane PERIOD-ENDS-2: pf-region ct_ends's per-post period ends, knots, gradient and kind when the window publishes them
    /// (read only by ct_sharp_ends=1); None otherwise.
    pub ends: Option<CtEnds>,
}

/// Lane PERIOD-ENDS-2: pf-region ct_ends fields (mm as published; i32::MIN = none).
pub struct CtEnds {
    pub s0: Vec<i32>,
    pub s1: Vec<i32>,
    pub s_dn: Vec<i32>,
    pub s_up: Vec<i32>,
    pub h1: Vec<i32>,
    pub z1: Vec<i32>,
    pub h2: Vec<i32>,
    pub z2: Vec<i32>,
    pub grad_permille: Vec<u16>,
    pub kind: Vec<u8>,
}

impl T1bWindow {
    pub fn load(dir: &Path) -> R<T1bWindow> {
        let bytes = rd(&dir.join("window.json"))?;
        let v: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
        let h = v.get("header").ok_or("window.json: no header")?.clone();
        let digest = v.get("digest").and_then(Value::as_str).ok_or("window.json: no digest")?.to_string();
        let got = pf_contracts::artifact::digest_canonical(&h).map_err(|e| e.to_string())?;
        if got != digest { return Err(format!("T1b window header digest {got} != stated {digest}; refusing")); }
        let mut blobs = std::collections::BTreeMap::new();
        for f in h.get("fields").and_then(Value::as_array).ok_or("fields")? {
            let name = f.get("name").and_then(Value::as_str).ok_or("field name")?.to_string();
            let file = f.get("file").and_then(Value::as_str).ok_or("field file")?;
            let b = rd(&dir.join(file))?;
            if Some(blake(&b).as_str()) != f.get("digest").and_then(Value::as_str) { return Err(format!("T1b field {name} digest mismatch; refusing")); }
            blobs.insert(name, b);
        }
        let l = h.get("lattice").ok_or("lattice")?;
        let gi = |k: &str| l.get(k).and_then(Value::as_i64).ok_or_else(|| format!("lattice.{k}"));
        let n = gi("n")? as usize;
        let b = |k: &str| blobs.get(k).cloned().ok_or_else(|| format!("T1b field {k} missing"));
        let evolved:Option<crate::contrast::Evolved>=h.pointer("/prep/net/evolved_surface").map(|v|serde_json::from_value(v.clone()).map_err(|e|e.to_string())).transpose()?;
        if let Some(e)=&evolved {e.validate()?;}
        let fine_meta=h.pointer("/prep/net/evolved_surface/escarpment");
        let fine_schema=evolved.as_ref().is_some_and(|e|e.schema==pf_region::t1b_escarp::EVOLVED_SCHEMA);
        if fine_schema!=fine_meta.is_some() {return Err("ESCARP-1 descriptor/schema mismatch; refusing coarse-mean fallback or obsolete prototype carrier".into());}
        let escarpment=if let Some(meta)=fine_meta {
            Some(crate::escarpment::FineSurface::load(meta,gi("gx0")?,gi("gy0")?,n,&b("escarp_base_mm")?,&b("escarp_delta_xyz")?)?)
        } else { None };
        let w = T1bWindow {
            escarpment,
            evolved,
            key: format!("{}.{}", gi("rx")?, gi("ry")?),
            digest,
            recipe_digest: h.get("recipe_digest").and_then(Value::as_str).unwrap_or("?").to_string(),
            generator_version: h.pointer("/generator/version").and_then(Value::as_str).unwrap_or("?").to_string(),
            gx0: gi("gx0")?, gy0: gi("gy0")?, n,
            height_m: i32s(&b("height_mm")?).into_iter().map(|v| v as f64 / 1000.0).collect(),
            fixed: b("fixed")?,
            drainage_area_m2: f32s(&b("drainage_area_m2")?),
            t0_cell: u32s(&b("t0_cell")?),
            base_m: blobs.get("base_mm").map(|v| i32s(v).into_iter().map(|v| v as f64 / 1000.0).collect()),
            exposed_litho: blobs.get("exposed_litho").cloned(),
            over: blobs.get("over").cloned(),
            volcano_id: blobs.get("volcano_id").cloned(),
            volcanoes: h.pointer("/prep/net/volcanoes").and_then(Value::as_array).map(|a| a.iter().map(|v| {
                let s = |k: &str| v.get(k).and_then(Value::as_str).unwrap_or("").to_string();
                let f = |k: &str| v.get(k).and_then(Value::as_f64).unwrap_or(f64::NAN);
                crate::strata_prov::VolcInfo { id: s("id"), kind: s("kind"), shape: s("shape"), age_my: f("age_my"), stage: s("stage"), denudation_m: f("denudation_m") }
            }).collect()).unwrap_or_default(),
            volcano_geo: h.pointer("/prep/net/volcanoes").and_then(Value::as_array).map(|a| a.iter().map(|v| {
                let xy = v.get("vent_xy_m").and_then(Value::as_array);
                let g = |i: usize| xy.and_then(|q| q.get(i)).and_then(Value::as_f64).unwrap_or(f64::NAN);
                let f = |k: &str| v.get(k).and_then(Value::as_f64).unwrap_or(f64::NAN);
                [g(0), g(1), f("r_influence_m"), f("h_m")]
            }).collect()).unwrap_or_default(),
            volcano_flows: h.pointer("/prep/net/volcanoes").and_then(Value::as_array).map(|a| a.iter().map(|v| {
                v.pointer("/cf/flows").and_then(Value::as_array).map(|fl| fl.iter().filter_map(|f| {
                    let q = f.as_array()?;
                    if q.len() < 8 { return None; }
                    let mut o = [0.0f64; 8];
                    for (i, x) in q.iter().take(8).enumerate() { o[i] = x.as_f64()?; }
                    Some(o)
                }).collect()).unwrap_or_default()
            }).collect()).unwrap_or_default(),
            edifice_style: blobs.get("edifice_style").cloned(),
            volcano_ve_units: h.pointer("/prep/net/volcanoes").and_then(Value::as_array).map(|a| a.iter().map(|v| v.pointer("/ve/units").cloned()).collect()).unwrap_or_default(),
            volcano_ve_counts: h.pointer("/prep/net/volcanoes").and_then(Value::as_array).map(|a| a.iter().map(|v| v.pointer("/ve/counts").cloned()).collect()).unwrap_or_default(),
            volcano_r_body: h.pointer("/prep/net/volcanoes").and_then(Value::as_array).map(|a| a.iter().map(|v| v.get("r_main_m").and_then(Value::as_f64)).collect()).unwrap_or_default(),
            volcano_ve_wind: h.pointer("/prep/net/volcanoes").and_then(Value::as_array).map(|a| a.iter().map(|v| v.pointer("/ve/wind/chart_azimuth_rad").and_then(Value::as_f64)).collect()).unwrap_or_default(),
            pre_cut_m: None,
            ct: None,
            runoff_weight: blobs.get("runoff_weight").map(|b| f32s(b)),
        };
        let mut w = w;
        if let Some(u) = blobs.get("uncut_mm") {
            let u = i32s(u);
            if u.len() == w.height_m.len() { w.pre_cut_m = Some(u.iter().zip(w.height_m.iter()).map(|(&q, &h)| if q == i32::MIN { h } else { q as f64 / 1000.0 }).collect()); }
        }
        if let (Some(b0), Some(r0), Some(c0)) = (blobs.get("ct_base_mm"), blobs.get("ct_regime"), blobs.get("ct_cliff_permille")) {
            let rn = h.pointer("/recipe/params/rivernet");
            let pm = |k: &str, d: f64| rn.and_then(|v| v.get(k)).and_then(Value::as_i64).map(|v| v as f64 / 1000.0).unwrap_or(d);
            w.ct = Some(CtWindow {
                base_m: i32s(b0).into_iter().map(|v| v as f64 / 1000.0).collect(),
                regime: r0.clone(),
                cliff_permille: c0.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect(),
                params: rn.and_then(|v| v.pointer("/ls/tb/tc/ct")).cloned().unwrap_or(Value::Null),
                rockmass: (pm("rockmass_min_permille", 0.6), pm("rockmass_max_permille", 1.5)),
                tread: match (blobs.get("ct_tread_permille"), blobs.get("ct_tread_tan_permille")) {
                    (Some(a), Some(b)) => Some((a.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect(), b.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect())),
                    _ => None,
                },
                ends: {
                    let gi = |k: &str| blobs.get(k).map(|v| i32s(v));
                    match (gi("ct_s0_mm"), gi("ct_s1_mm"), gi("ct_s_dn_mm"), gi("ct_s_up_mm"), gi("ct_knot_h1_mm"), gi("ct_knot_z1_mm"), gi("ct_knot_h2_mm"), gi("ct_knot_z2_mm"), blobs.get("ct_grad_permille"), blobs.get("ct_period_kind")) {
                        (Some(s0), Some(s1), Some(s_dn), Some(s_up), Some(h1), Some(z1), Some(h2), Some(z2), Some(g), Some(k)) => Some(CtEnds { s0, s1, s_dn, s_up, h1, z1, h2, z2,
                            grad_permille: g.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect(), kind: k.clone() }),
                        _ => None,
                    }
                },
            });
        }
        if w.height_m.len() != n * n || w.fixed.len() != n * n { return Err("T1b field sizes".into()); }
        Ok(w)
    }
    #[inline]
    pub fn idx(&self, gx: i64, gy: i64) -> Option<usize> {
        let (x, y) = (gx - self.gx0, gy - self.gy0);
        if x < 0 || y < 0 || x >= self.n as i64 || y >= self.n as i64 { return None; }
        Some(y as usize * self.n + x as usize)
    }
    /// Height at a T1b post, clamped to the lattice.
    #[inline]
    pub fn h(&self, gx: i64, gy: i64) -> f64 {
        let x = (gx - self.gx0).clamp(0, self.n as i64 - 1) as usize;
        let y = (gy - self.gy0).clamp(0, self.n as i64 - 1) as usize;
        self.height_m[y * self.n + x]
    }
    /// B1: Catmull-Rom bicubic at chart position (x, y) metres.
    pub fn bicubic(&self, x: f64, y: f64) -> f64 {
        let fx = x / crate::T1B_M;
        let fy = y / crate::T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let wx = cr(tx);
        let wy = cr(ty);
        let mut acc = 0.0;
        for (j, wyj) in wy.iter().enumerate() {
            let mut row = 0.0;
            for (i, wxi) in wx.iter().enumerate() {
                row += wxi * self.h(ix - 1 + i as i64, iy - 1 + j as i64);
            }
            acc += wyj * row;
        }
        acc
    }
    /// Lane BS: Catmull-Rom bicubic of another full-window lattice field (same layout as height_m), clamped at the edge.
    pub fn bicubic_of(&self, f: &[f64], x: f64, y: f64) -> f64 {
        let fx = x / crate::T1B_M;
        let fy = y / crate::T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (wx, wy) = (cr(fx - ix as f64), cr(fy - iy as f64));
        let g = |gx: i64, gy: i64| { let x = (gx - self.gx0).clamp(0, self.n as i64 - 1) as usize; let y = (gy - self.gy0).clamp(0, self.n as i64 - 1) as usize; f[y * self.n + x] };
        let mut acc = 0.0;
        for (j, wyj) in wy.iter().enumerate() {
            let mut row = 0.0;
            for (i, wxi) in wx.iter().enumerate() { row += wxi * g(ix - 1 + i as i64, iy - 1 + j as i64); }
            acc += wyj * row;
        }
        acc
    }
    /// Nearest T1b post of chart position (x, y) metres.
    #[inline]
    pub fn nearest(&self, x: f64, y: f64) -> Option<usize> {
        self.idx((x / crate::T1B_M).round() as i64, (y / crate::T1B_M).round() as i64)
    }
    /// Lane CLIMATE-RUNOFF: the runoff weight R / r_ref at (x, y) chart metres, bilinear between T1b posts (clamped to the
    /// lattice); None when the window publishes no runoff_weight.
    pub fn runoff_weight_at(&self, x: f64, y: f64) -> Option<f64> {
        let rw = self.runoff_weight.as_ref()?;
        let fx = x / crate::T1B_M;
        let fy = y / crate::T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let g = |gx: i64, gy: i64| { let x = (gx - self.gx0).clamp(0, self.n as i64 - 1) as usize; let y = (gy - self.gy0).clamp(0, self.n as i64 - 1) as usize; rw[y * self.n + x] as f64 };
        Some((g(ix, iy) * (1.0 - tx) + g(ix + 1, iy) * tx) * (1.0 - ty) + (g(ix, iy + 1) * (1.0 - tx) + g(ix + 1, iy + 1) * tx) * ty)
    }
    /// Bilinear drainage area (m2) at (x, y).
    /// Lane INTS2n4 ringflow_ch: the bilinear drainage area over the free (fixed 0) T1b posts only, weights renormalised over
    /// those corners (0 if none).
    pub fn area_bilinear_free(&self, x: f64, y: f64) -> f64 {
        let fx = x / crate::T1B_M;
        let fy = y / crate::T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let (mut s, mut ws) = (0.0f64, 0.0f64);
        for (gx, gy, wt) in [(ix, iy, (1.0 - tx) * (1.0 - ty)), (ix + 1, iy, tx * (1.0 - ty)), (ix, iy + 1, (1.0 - tx) * ty), (ix + 1, iy + 1, tx * ty)] {
            if let Some(i) = self.idx(gx, gy) { if self.fixed[i] == 0 { s += wt * self.drainage_area_m2[i] as f64; ws += wt; } }
        }
        if ws > 0.0 { s / ws } else { 0.0 }
    }
    pub fn area_bilinear(&self, x: f64, y: f64) -> f64 {
        let fx = x / crate::T1B_M;
        let fy = y / crate::T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let a = |gx: i64, gy: i64| self.idx(gx, gy).map(|i| self.drainage_area_m2[i] as f64).unwrap_or(0.0);
        (a(ix, iy) * (1.0 - tx) + a(ix + 1, iy) * tx) * (1.0 - ty) + (a(ix, iy + 1) * (1.0 - tx) + a(ix + 1, iy + 1) * tx) * ty
    }
}

#[inline]
fn cr(t: f64) -> [f64; 4] {
    let (t2, t3) = (t * t, t * t * t);
    [(-t3 + 2.0 * t2 - t) * 0.5, (3.0 * t3 - 5.0 * t2 + 2.0) * 0.5, (-3.0 * t3 + 4.0 * t2 + t) * 0.5, (t3 - t2) * 0.5]
}

/// The T1a skeleton (lod-8 lattice, 513 x 513 nodes, gx/gy in -256..=256).
pub struct T1a {
    pub digest: String,
    pub channel: Vec<u8>,
    pub discharge_m3_s: Vec<f32>,
}
pub const T1A_GMIN: i64 = -256;
pub const T1A_N: usize = 513;

impl T1a {
    pub fn load(dir: &Path) -> R<T1a> {
        let sk = rd(&dir.join("skeleton.json"))?;
        let v: Value = serde_json::from_slice(&sk).map_err(|e| e.to_string())?;
        let ch = rd(&dir.join("fields/channel.u8"))?;
        let q = rd(&dir.join("fields/discharge_m3_s.f32"))?;
        for (name, b) in [("channel", &ch), ("discharge_m3_s", &q)] {
            let want = v.pointer("/header/fields").and_then(Value::as_array).and_then(|a| a.iter().find(|f| f.get("name").and_then(Value::as_str) == Some(name)))
                .and_then(|f| f.get("digest")).and_then(Value::as_str).ok_or_else(|| format!("skeleton field {name} not listed"))?;
            if blake(b) != want { return Err(format!("T1a field {name} digest mismatch; refusing")); }
        }
        let digest = v.get("digest").and_then(Value::as_str).unwrap_or("?").to_string();
        Ok(T1a { digest, channel: ch, discharge_m3_s: f32s(&q) })
    }
    pub fn node(&self, gx8: i64, gy8: i64) -> Option<usize> {
        let (x, y) = (gx8 - T1A_GMIN, gy8 - T1A_GMIN);
        if x < 0 || y < 0 || x >= T1A_N as i64 || y >= T1A_N as i64 { return None; }
        Some(y as usize * T1A_N + x as usize)
    }
    /// Discharge of the nearest T1a channel node (channel == 1) to the T1b post (gx, gy), searching 5 x 5 nodes.
    pub fn channel_q_near(&self, gx: i64, gy: i64) -> Option<f64> {
        let (cx, cy) = (gx as f64 / crate::L8_TO_L4 as f64, gy as f64 / crate::L8_TO_L4 as f64);
        let (rx, ry) = (cx.round() as i64, cy.round() as i64);
        let mut best: Option<(f64, f64)> = None;
        for dy in -2..=2 { for dx in -2..=2 {
            if let Some(n) = self.node(rx + dx, ry + dy) {
                if self.channel[n] == 1 {
                    let d = ((rx + dx) as f64 - cx).hypot((ry + dy) as f64 - cy);
                    if best.map(|b| d < b.0).unwrap_or(true) { best = Some((d, self.discharge_m3_s[n] as f64)); }
                }
            }
        } }
        best.map(|b| b.1)
    }
}

/// T0 per-cell data the T2 cut needs: elevation, the stratigraphy CSR column, the uplift gradient.
pub struct T0Strata {
    pub elevation_m: Vec<f32>,
    pub off: Vec<u32>,
    pub lith: Vec<u8>,
    pub th: Vec<f32>,
    pub grad_e: Vec<f32>,
    pub grad_n: Vec<f32>,
    pub sea_level_m: f64,
    pub digests: Vec<(String, String)>,
}

impl T0Strata {
    /// `run_root` = lane T2 driver output: `run/terrain/...` (published) and `final/stratafull_*`, `final/strata_grad_*`.
    pub fn load(run_root: &Path) -> R<T0Strata> {
        let p = |s: &str| -> PathBuf { run_root.join(s) };
        let mut digests = Vec::new();
        let mut get = |role: &str, s: &str| -> R<Vec<u8>> { let b = rd(&p(s))?; digests.push((role.to_string(), blake(&b))); Ok(b) };
        let el = get("t0.terrain.elevation_m", "run/terrain/fields/elevation_m.f32")?;
        let off = get("t0.stratigraphy.full.offset", "final/stratafull_off.u32")?;
        let lith = get("t0.stratigraphy.full.lith", "final/stratafull_lith.u8")?;
        let th = get("t0.stratigraphy.full.thickness", "final/stratafull_th.f32")?;
        let ge = get("t0.stratigraphy.uplift_grad_east", "final/strata_grad_e.f32")?;
        let gn = get("t0.stratigraphy.uplift_grad_north", "final/strata_grad_n.f32")?;
        let tj = rd(&p("run/terrain/artifact.json"))?;
        let tv: Value = serde_json::from_slice(&tj).map_err(|e| e.to_string())?;
        let sea = tv.pointer("/header/domain/sea_level_m").and_then(Value::as_f64).ok_or("terrain sea_level_m")?;
        digests.push(("t0.terrain.artifact".into(), tv.get("digest").and_then(Value::as_str).unwrap_or("?").to_string()));
        digests.sort();
        Ok(T0Strata { elevation_m: f32s(&el), off: u32s(&off), lith, th: f32s(&th), grad_e: f32s(&ge), grad_n: f32s(&gn), sea_level_m: sea, digests })
    }
}

/// Run configuration file (JSON): upstream paths and the anchor.
#[derive(Clone, Debug, serde::Deserialize, serde::Serialize)]
pub struct RunConfig {
    pub t1b_window_dir: String,
    pub t1a_dir: String,
    pub t0_run_root: String,
    pub site_record: String,
    pub world_seed: u64,
    pub chart: u32,
}

/// The site point in chart metres (the record's `chart.site_in_chart_mm`).
pub fn site_anchor_m(site_record: &Path) -> R<((f64, f64), String)> {
    let b = rd(site_record)?;
    let v: Value = serde_json::from_slice(&b).map_err(|e| e.to_string())?;
    let a = v.pointer("/chart/site_in_chart_mm").and_then(Value::as_array).ok_or("site_in_chart_mm")?;
    let x = a[0].as_i64().ok_or("x")? as f64 / 1000.0;
    let y = a[1].as_i64().ok_or("y")? as f64 / 1000.0;
    Ok(((x, y), blake(&b)))
}
