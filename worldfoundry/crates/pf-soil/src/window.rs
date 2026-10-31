//! Reader of a published T1b window (`vz.t1b-window/v1`): the lattice, chart, the fields this lane reads and the
//! volcano records (lane VE units with ages) from the header.

use serde_json::Value;
use std::path::{Path, PathBuf};

pub const POST_M: f64 = 32.0;

#[derive(Clone, Debug)]
pub struct Unit { pub age_yr: f64, pub style: u8, pub vei: u8, pub kind: String }

#[derive(Clone, Debug)]
pub struct Volc {
    pub id: String,
    pub vent: [f64; 2],
    pub crater_r_m: f64,
    pub age_my: f64,
    pub stage: String,
    /// History window T_h (yr): the composite body is older than every mappable unit.
    pub t_h_yr: f64,
    pub units: Vec<Unit>,
}

impl Volc {
    pub fn youngest_yr(&self) -> f64 { self.units.iter().map(|u| u.age_yr).fold(f64::INFINITY, f64::min) }
}

pub struct Window {
    pub dir: PathBuf,
    pub n: usize,
    pub gx0: i64,
    pub gy0: i64,
    pub halo: usize,
    pub origin_lat_e7: i64,
    pub origin_lon_e7: i64,
    pub radius_mm: i64,
    pub digest: String,
    pub recipe_digest: String,
    pub window_key: String,
    pub height_m: Vec<f64>,
    pub fixed: Vec<u8>,
    pub area: Vec<f64>,
    pub over: Vec<u8>,
    pub volcano_id: Vec<u8>,
    pub edifice_unit: Vec<u16>,
    pub edifice_style: Vec<u8>,
    pub ct_regime: Vec<u8>,
    pub shore_m: Option<Vec<f64>>,
    pub volcanoes: Vec<Volc>,
    pub fields_read: Vec<(String, String)>,
}

fn field_file(h: &Value, name: &str) -> Option<(String, String, String)> {
    h.get("fields")?.as_array()?.iter().find(|f| f.get("name").and_then(Value::as_str) == Some(name))
        .map(|f| (f["file"].as_str().unwrap_or("").to_string(), f["dtype"].as_str().unwrap_or("").to_string(), f["digest"].as_str().unwrap_or("").to_string()))
}

impl Window {
    pub fn load(dir: &Path) -> Result<Window, String> {
        let raw = std::fs::read(dir.join("window.json")).map_err(|e| format!("{}: {e}", dir.display()))?;
        let w: Value = serde_json::from_slice(&raw).map_err(|e| e.to_string())?;
        let h = w.get("header").ok_or("window.json: no header")?;
        let lat = &h["lattice"];
        let n = lat["n"].as_u64().ok_or("lattice n")? as usize;
        let gx0 = lat["gx0"].as_i64().ok_or("lattice gx0")?;
        let gy0 = lat["gy0"].as_i64().ok_or("lattice gy0")?;
        let halo = lat["halo_posts"].as_u64().ok_or("lattice halo_posts")? as usize;
        let sp = h["recipe"]["lattice"]["spacing_mm"].as_i64().ok_or("recipe lattice spacing_mm")?;
        if sp != 32_000 { return Err(format!("window spacing {sp} mm; pf-soil reads the 32 m T1b lattice")); }
        let ch = &h["recipe"]["chart"];
        let (origin_lat_e7, origin_lon_e7, radius_mm) = (ch["origin_lat_e7"].as_i64().ok_or("chart lat")?, ch["origin_lon_e7"].as_i64().ok_or("chart lon")?, ch["datum_radius_mm"].as_i64().ok_or("chart radius")?);
        let nn = n * n;
        let mut fields_read = Vec::new();
        let mut rd = |name: &str, need: bool| -> Result<Option<(Vec<u8>, String)>, String> {
            match field_file(h, name) {
                None => if need { Err(format!("window field {name} missing")) } else { Ok(None) },
                Some((file, dt, dg)) => {
                    let b = std::fs::read(dir.join(&file)).map_err(|e| format!("{file}: {e}"))?;
                    let dg_now = format!("blake3:{}", blake3::hash(&b).to_hex());
                    if dg_now != dg { return Err(format!("window field {name}: digest {dg_now} != header {dg}")); }
                    fields_read.push((name.to_string(), dg));
                    Ok(Some((b, dt)))
                }
            }
        };
        let chk = |b: &[u8], w: usize, name: &str| -> Result<(), String> { if b.len() != nn * w { Err(format!("{name}: {} bytes, expected {}", b.len(), nn * w)) } else { Ok(()) } };
        let (b, _) = rd("height_mm", true)?.unwrap(); chk(&b, 4, "height_mm")?;
        let height_m: Vec<f64> = b.chunks_exact(4).map(|c| i32::from_le_bytes(c.try_into().unwrap()) as f64 / 1000.0).collect();
        let (fixed, _) = rd("fixed", true)?.unwrap(); chk(&fixed, 1, "fixed")?;
        let (b, _) = rd("drainage_area_m2", true)?.unwrap(); chk(&b, 4, "drainage_area_m2")?;
        let area: Vec<f64> = b.chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap()) as f64).collect();
        let u8f = |o: Option<(Vec<u8>, String)>, d: u8| -> Vec<u8> { o.map(|x| x.0).unwrap_or_else(|| vec![d; nn]) };
        let over = u8f(rd("over", false)?, 0);
        let volcano_id = u8f(rd("volcano_id", false)?, 255);
        let edifice_style = u8f(rd("edifice_style", false)?, 0);
        let ct_regime = u8f(rd("ct_regime", false)?, 0);
        let edifice_unit: Vec<u16> = match rd("edifice_unit", false)? { Some((b, _)) => b.chunks_exact(2).map(|c| u16::from_le_bytes(c.try_into().unwrap())).collect(), None => vec![65535; nn] };
        let shore_m = rd("shore_mm", false)?.map(|(b, _)| b.chunks_exact(4).map(|c| i32::from_le_bytes(c.try_into().unwrap()) as f64 / 1000.0).collect::<Vec<f64>>());
        for (v, nm) in [(&over, "over"), (&volcano_id, "volcano_id"), (&edifice_style, "edifice_style"), (&ct_regime, "ct_regime")] { chk(v, 1, nm)?; }
        if edifice_unit.len() != nn { return Err("edifice_unit size".into()); }
        if let Some(s) = &shore_m { if s.len() != nn { return Err("shore_mm size".into()); } }
        let mut volcanoes = Vec::new();
        if let Some(vs) = h["prep"]["net"]["volcanoes"].as_array() {
            for v in vs {
                let ve = &v["ve"];
                let units = ve["units"].as_array().map(|a| a.iter().map(|u| Unit {
                    age_yr: u["age_yr"].as_f64().unwrap_or(f64::INFINITY), style: u["style"].as_u64().unwrap_or(0) as u8,
                    vei: u["vei"].as_u64().unwrap_or(0) as u8, kind: u["kind"].as_str().unwrap_or("").to_string() }).collect()).unwrap_or_default();
                let vent = v["vent_xy_m"].as_array().map(|a| [a[0].as_f64().unwrap_or(f64::NAN), a[1].as_f64().unwrap_or(f64::NAN)]).unwrap_or([f64::NAN, f64::NAN]);
                volcanoes.push(Volc { id: v["id"].as_str().unwrap_or("").to_string(), vent, crater_r_m: v["crater_radius_m"].as_f64().unwrap_or(0.0),
                    age_my: v["age_my"].as_f64().unwrap_or(f64::INFINITY), stage: v["stage"].as_str().unwrap_or("").to_string(),
                    t_h_yr: ve["t_h_ky"].as_f64().map(|k| k * 1000.0).unwrap_or(f64::INFINITY), units });
            }
        }
        Ok(Window { dir: dir.to_path_buf(), n, gx0, gy0, halo, origin_lat_e7, origin_lon_e7, radius_mm,
            digest: w["digest"].as_str().unwrap_or("").to_string(), recipe_digest: h["recipe_digest"].as_str().unwrap_or("").to_string(),
            window_key: h["window"].as_str().unwrap_or("").to_string(),
            height_m, fixed, area, over, volcano_id, edifice_unit, edifice_style, ct_regime, shore_m, volcanoes, fields_read })
    }

    pub fn idx(&self, gx: i64, gy: i64) -> Option<usize> {
        let (i, j) = (gx - self.gx0, gy - self.gy0);
        if i < 0 || j < 0 || i >= self.n as i64 || j >= self.n as i64 { None } else { Some(j as usize * self.n + i as usize) }
    }
    pub fn g(&self, k: usize) -> (i64, i64) { (self.gx0 + (k % self.n) as i64, self.gy0 + (k / self.n) as i64) }
}
