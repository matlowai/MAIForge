//! Quantized fields and their files (little-endian, row-major, row 0 south) with blake3 digests.

use crate::zonation::{PostOut, C_NODATA};
use serde_json::{json, Value};
use std::path::Path;

pub struct Fields {
    pub veg_cover_permille: Vec<u16>,
    pub veg_class: Vec<u8>,
    pub limit_code: Vec<u8>,
    pub treeline_margin_m: Vec<i16>,
    pub nival_margin_m: Vec<i16>,
    pub wetness_q: Vec<u8>,
    pub insolation_permille: Vec<u16>,
    pub soil_depth_cm: Vec<u16>,
    pub t_gs_dc: Vec<i16>,
}

fn i16m(v: f64) -> i16 { if v.is_finite() { v.round().clamp(-32767.0, 32767.0) as i16 } else { i16::MIN } }

/// Quantize posts; `pick` maps output index -> source index (None = no data).
pub fn quantize(n_out: usize, pick: &dyn Fn(usize) -> Option<usize>, out: &[PostOut], twi: &[f64], r_sun: &[f64]) -> Fields {
    let mut f = Fields { veg_cover_permille: vec![u16::MAX; n_out], veg_class: vec![C_NODATA; n_out], limit_code: vec![255; n_out],
        treeline_margin_m: vec![i16::MIN; n_out], nival_margin_m: vec![i16::MIN; n_out], wetness_q: vec![255; n_out],
        insolation_permille: vec![u16::MAX; n_out], soil_depth_cm: vec![u16::MAX; n_out], t_gs_dc: vec![i16::MIN; n_out] };
    for o in 0..n_out {
        let Some(k) = pick(o) else { continue };
        let p = &out[k];
        if p.class == C_NODATA { continue; }
        f.veg_cover_permille[o] = (p.cover * 1000.0).round().clamp(0.0, 1000.0) as u16;
        f.veg_class[o] = p.class;
        f.limit_code[o] = p.limit;
        f.treeline_margin_m[o] = i16m(p.tl_margin_m);
        f.nival_margin_m[o] = i16m(p.niv_margin_m);
        f.wetness_q[o] = if twi[k].is_finite() { (twi[k] * 10.0).round().clamp(0.0, 254.0) as u8 } else { 255 };
        f.insolation_permille[o] = if r_sun[k].is_finite() { (r_sun[k] * 1000.0).round().clamp(0.0, 65534.0) as u16 } else { u16::MAX };
        f.soil_depth_cm[o] = (p.soil_m * 100.0).round().clamp(0.0, 65534.0) as u16;
        f.t_gs_dc[o] = if p.t_gs.is_finite() { i16m(p.t_gs * 10.0) } else { i16::MIN };
    }
    f
}

pub const FIELD_NOTES: [(&str, &str, &str, &str); 9] = [
    ("veg_cover_permille", "u16", "permille", "vegetation (vascular) cover 0..1000 = f_T x f_snow x f_W x f_S x f_age x f_kill (wet classes: max with wet_cover x f_T x f_snow x f_age x f_kill); 65535 no data"),
    ("veg_class", "u8", "code", "zonation class (table classes); 255 no data"),
    ("limit_code", "u8", "code", "factor that limits cover most when the smallest factor is < 0.95: 0 none, 1 temperature, 2 moisture, 3 substrate, 4 volcanic age, 5 vent / tephra kill; 255 no data"),
    ("treeline_margin_m", "i16", "m", "height of the post's own Paulsen-Koerner tree line minus the post height (> 0 below the tree line), aspect term included; -32768 no data"),
    ("nival_margin_m", "i16", "m", "Ohmura snow-line height minus the post height (> 0 below the snow line); -32768 no data"),
    ("wetness_q", "u8", "0.1 ln(m)", "topographic wetness index ln(a / tan b) x 10 (a = effective drainage area / 32 m contour width, the T1b climate-key area; 2 m: bilinear of the 32 m TWI); 255 no data"),
    ("insolation_permille", "u16", "permille", "warm-half-year clear-sky direct-beam ratio slope / flat (no horizon shading) x 1000; 65535 no data"),
    ("soil_depth_cm", "u16", "cm", "soil depth PROXY (not a regolith model): soil_dmax x (1 - exp(-t / tau_soil)) on volcanic units x slope taper x regime factor x (1 - debris) x (0.25 + 0.75 f_T f_snow); 65535 no data"),
    ("t_gs_dc", "i16", "0.1 degC", "mean growing-season air temperature (days above 0.9 C of the seasonal cycle), aspect term included; -32768 no season or no data"),
];

pub fn write(dir: &Path, f: &Fields) -> Result<Vec<Value>, String> {
    std::fs::create_dir_all(dir.join("fields")).map_err(|e| e.to_string())?;
    let mut out = Vec::new();
    let mut put = |name: &str, bytes: Vec<u8>| -> Result<(), String> {
        let (_, dt, unit, note) = FIELD_NOTES.iter().find(|x| x.0 == name).copied().unwrap();
        let file = format!("fields/{name}.{dt}");
        std::fs::write(dir.join(&file), &bytes).map_err(|e| format!("{file}: {e}"))?;
        out.push(json!({"name": name, "file": file, "dtype": dt, "units": unit, "note": note, "digest": format!("blake3:{}", blake3::hash(&bytes).to_hex())}));
        Ok(())
    };
    let u16b = |v: &[u16]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    let i16b = |v: &[i16]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    put("veg_cover_permille", u16b(&f.veg_cover_permille))?;
    put("veg_class", f.veg_class.clone())?;
    put("limit_code", f.limit_code.clone())?;
    put("treeline_margin_m", i16b(&f.treeline_margin_m))?;
    put("nival_margin_m", i16b(&f.nival_margin_m))?;
    put("wetness_q", f.wetness_q.clone())?;
    put("insolation_permille", u16b(&f.insolation_permille))?;
    put("soil_depth_cm", u16b(&f.soil_depth_cm))?;
    put("t_gs_dc", i16b(&f.t_gs_dc))?;
    Ok(out)
}

/// Class histogram and cover statistics over the published posts.
pub fn stats(f: &Fields) -> Value {
    let mut hist = [0u64; 256];
    let mut lim = [0u64; 256];
    let (mut cs, mut cn) = (0f64, 0u64);
    for k in 0..f.veg_class.len() {
        hist[f.veg_class[k] as usize] += 1;
        lim[f.limit_code[k] as usize] += 1;
        if f.veg_cover_permille[k] != u16::MAX { cs += f.veg_cover_permille[k] as f64; cn += 1; }
    }
    let h: serde_json::Map<String, Value> = (0..256).filter(|&c| hist[c] > 0).map(|c| (c.to_string(), json!(hist[c]))).collect();
    let l: serde_json::Map<String, Value> = (0..256).filter(|&c| lim[c] > 0).map(|c| (c.to_string(), json!(lim[c]))).collect();
    json!({"posts": f.veg_class.len(), "data_posts": cn, "class_posts": h, "limit_posts": l, "cover_mean_permille": if cn > 0 { cs / cn as f64 } else { 0.0 }})
}
