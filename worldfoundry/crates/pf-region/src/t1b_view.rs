//! The T1a skeleton as T1b reads it (lane TR, TER-REGION): the published
//! node fields, verified against the skeleton header digest, the per-field
//! digests and the lod-8 tile digests of the tiles a window touches
//! (ADR 044 section 5: a regenerated or re-read layer must match its digest).
//!
//! Two sources give the same view: the published artifact on disk
//! (`skeleton.json` plus `fields/`), or an in-memory `Skeleton` (tests),
//! which goes through the same `publish::header` bytes.

use crate::skeleton::Skeleton;
use crate::*;
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;

/// The node fields T1b uses, on the lod-8 chart lattice (513 x 513).
pub struct SkeletonView {
    /// blake3 over canonical JSON of the skeleton header (the layer digest).
    pub digest: String,
    pub recipe_digest: String,
    pub physical_digest: String,
    /// Whether the verified parent recipe explicitly carries physical draw roles.
    pub physical_keying: bool,
    pub physical_climate_twin_digest: Option<String>,
    pub kind: Vec<u8>,
    pub channel: Vec<u8>,
    pub receiver: Vec<u32>,
    pub trunk_level_m: Vec<f32>,
    pub discharge_m3_s: Vec<f32>,
    pub t0_cell: Vec<u32>,
    pub outlet_basin: Vec<u32>,
    pub segment: Vec<u32>,
    pub sea_level_m: f64,
    pub threshold_m3_s: f64,
    pub corridor_half_width_mm: i64,
    /// lod-8 tile key ("8.ix.iy") -> published digest.
    pub tile_digests: BTreeMap<String, String>,
    /// Lane CLIMATE-RUNOFF (a climate T1a only): the climate artifact digest (upstream t0.climate), the climate params,
    /// and the digest of the climate-free twin skeleton (T1b's draw twin).
    pub climate_digest: Option<String>,
    pub climate_params: Option<crate::climate_runoff::ClimateRunoffParamsV1>,
    pub climate_twin_digest: Option<String>,
    /// Lane CHANNEL-CENTRELINE: the T1a segments' keyed planform curves as published in the digested header (start, end,
    /// amplitude, control offsets); `skeleton::curve` re-evaluated without draws (`t1b_cl`). Empty when the header has none.
    pub t1a_segments: Vec<crate::t1b_cl::T1aSegCurve>,
    /// Field order and element sizes as published (for tile digests).
    field_order: Vec<(String, usize)>,
    blobs: BTreeMap<String, Vec<u8>>,
}

fn u32s(b: &[u8]) -> Vec<u32> { b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
fn f32s(b: &[u8]) -> Vec<f32> { b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }

impl SkeletonView {
    /// Build from a header value, its stated digest and the field blobs by name.
    pub fn from_parts(header: &Value, stated_digest: &str, blobs: BTreeMap<String, Vec<u8>>) -> R<SkeletonView> {
        let recomputed = pf_contracts::artifact::digest_canonical(header).map_err(|e| e.to_string())?;
        if recomputed != stated_digest {
            return Err(format!("T1a skeleton header digest {recomputed} differs from the stated {stated_digest}; refusing"));
        }
        let mut field_order = Vec::new();
        for f in header.get("fields").and_then(Value::as_array).ok_or("skeleton header: no fields")? {
            let name = f.get("name").and_then(Value::as_str).ok_or("field name")?.to_string();
            let want = f.get("digest").and_then(Value::as_str).ok_or("field digest")?;
            let bytes = blobs.get(&name).ok_or_else(|| format!("T1a field {name} missing"))?;
            let got = format!("blake3:{}", blake3::hash(bytes).to_hex());
            if got != want { return Err(format!("T1a field {name} digest {got} differs from the header's {want}; refusing")); }
            let size = bytes.len() / NODES;
            if size * NODES != bytes.len() { return Err(format!("T1a field {name}: {} bytes", bytes.len())); }
            field_order.push((name, size));
        }
        let mut tile_digests = BTreeMap::new();
        for t in header.get("tile_digests").and_then(Value::as_array).ok_or("skeleton header: no tile digests")? {
            tile_digests.insert(t.get("key").and_then(Value::as_str).ok_or("tile key")?.to_string(),
                t.get("digest").and_then(Value::as_str).ok_or("tile digest")?.to_string());
        }
        let g = |n: &str| blobs.get(n).cloned().ok_or_else(|| format!("T1a field {n} missing"));
        let params = header.pointer("/recipe/params").ok_or("skeleton header: no recipe params")?;
        let physical_keying = header.pointer("/recipe/upstream").and_then(Value::as_array)
            .is_some_and(|u| u.iter().any(|e| e.get(0).and_then(Value::as_str).is_some_and(|r| r.starts_with("physical."))));
        Ok(SkeletonView {
            digest: stated_digest.to_string(),
            // Legacy readers must not acquire extra physical-projection parsing
            // requirements, nor use an unrelated physical twin for their draws.
            physical_digest: if physical_keying { crate::skeleton::physical_header_digest(header.clone())? } else { stated_digest.to_string() },
            physical_keying,
            physical_climate_twin_digest: if physical_keying { header.get("physical_climate_twin_digest").and_then(Value::as_str).map(str::to_string) } else { None },
            recipe_digest: header.get("recipe_digest").and_then(Value::as_str).ok_or("recipe_digest")?.to_string(),
            kind: g("kind")?,
            channel: g("channel")?,
            receiver: u32s(&g("receiver")?),
            trunk_level_m: f32s(&g("trunk_level_m")?),
            discharge_m3_s: f32s(&g("discharge_m3_s")?),
            t0_cell: u32s(&g("t0_cell")?),
            outlet_basin: u32s(&g("outlet_basin")?),
            segment: u32s(&g("segment")?),
            sea_level_m: header.get("sea_level_m").and_then(Value::as_f64).ok_or("sea_level_m")?,
            threshold_m3_s: header.get("river_threshold_m3_s").and_then(Value::as_f64).ok_or("river_threshold_m3_s")?,
            corridor_half_width_mm: params.get("corridor_half_width_mm").and_then(Value::as_i64).ok_or("corridor")?,
            climate_digest: header.pointer("/recipe/upstream").and_then(Value::as_array)
                .and_then(|u| u.iter().find(|e| e.get(0).and_then(Value::as_str) == Some("t0.climate")))
                .and_then(|e| e.get(1)).and_then(Value::as_str).map(str::to_string),
            climate_params: match params.get("climate_runoff") {
                Some(v) => Some(serde_json::from_value(v.clone()).map_err(|e| format!("T1a climate_runoff params: {e}"))?),
                None => None,
            },
            climate_twin_digest: header.get("climate_twin_digest").and_then(Value::as_str).map(str::to_string),
            tile_digests,
            t1a_segments: crate::t1b_cl::T1aSegCurve::from_header(header)?,
            field_order,
            blobs,
        })
    }

    /// Read the published artifact (`skeleton.json` and `fields/`).
    pub fn load(dir: &Path) -> R<SkeletonView> {
        let bytes = std::fs::read(dir.join("skeleton.json")).map_err(|e| format!("skeleton.json in {}: {e}", dir.display()))?;
        let v: Value = serde_json::from_slice(&bytes).map_err(|e| format!("skeleton.json parse: {e}"))?;
        let header = v.get("header").ok_or("skeleton.json: no header")?;
        let digest = v.get("digest").and_then(Value::as_str).ok_or("skeleton.json: no digest")?;
        let mut blobs = BTreeMap::new();
        for f in header.get("fields").and_then(Value::as_array).ok_or("no fields")? {
            let name = f.get("name").and_then(Value::as_str).ok_or("name")?;
            let file = f.get("file").and_then(Value::as_str).ok_or("file")?;
            blobs.insert(name.to_string(), std::fs::read(dir.join(file)).map_err(|e| format!("{file}: {e}"))?);
        }
        SkeletonView::from_parts(header, digest, blobs)
    }

    /// The same view from an in-memory skeleton (the bytes `publish` writes).
    pub fn from_skeleton(sk: &Skeleton) -> R<SkeletonView> {
        let (h, blobs) = crate::publish::header(sk)?;
        let digest = pf_contracts::artifact::digest_canonical(&h).map_err(|e| e.to_string())?;
        let v = serde_json::to_value(&h).map_err(|e| e.to_string())?;
        // Round-trip through text as the disk path does.
        let v: Value = serde_json::from_str(&serde_json::to_string(&v).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let mut by_name = BTreeMap::new();
        for (file, bytes) in blobs {
            let name = file.trim_start_matches("fields/").split('.').next().unwrap_or("").to_string();
            by_name.insert(name, bytes);
        }
        SkeletonView::from_parts(&v, &digest, by_name)
    }

    /// Recompute the digest of lod-8 tile (ix, iy) from the fields and
    /// compare it with the published one.
    pub fn verify_tile(&self, ix: i64, iy: i64) -> R<()> {
        let key = format!("{LOD}.{ix}.{iy}");
        let want = self.tile_digests.get(&key).ok_or_else(|| format!("T1a tile {key} not published"))?;
        let mut h = blake3::Hasher::new();
        for gy in iy * 128..=iy * 128 + 128 { for gx in ix * 128..=ix * 128 + 128 {
            let n = node_index(gx, gy).ok_or("tile post outside lattice")?;
            for (name, size) in &self.field_order {
                let b = &self.blobs[name];
                h.update(&b[n * size..(n + 1) * size]);
            }
        } }
        let got = format!("blake3:{}", h.finalize().to_hex());
        if &got != want { return Err(format!("T1a tile {key} digest {got} differs from published {want}; refusing")); }
        Ok(())
    }

    /// Verify every lod-8 tile whose posts a chart-mm rectangle touches;
    /// returns the keys verified.
    pub fn verify_tiles_for_rect(&self, min_mm: [i64; 2], max_mm: [i64; 2]) -> R<Vec<String>> {
        let side = 128 * SPACING_MM;
        let lo = |v: i64| (v - SPACING_MM).div_euclid(side);
        let hi = |v: i64| (v + SPACING_MM).div_euclid(side);
        let (krange_lo, krange_hi) = vz_tilekey::key_range(LOD);
        let mut keys = Vec::new();
        for iy in lo(min_mm[1]).max(krange_lo)..=hi(max_mm[1]).min(krange_hi - 1) {
            for ix in lo(min_mm[0]).max(krange_lo)..=hi(max_mm[0]).min(krange_hi - 1) {
                self.verify_tile(ix, iy)?;
                keys.push(format!("{LOD}.{ix}.{iy}"));
            }
        }
        Ok(keys)
    }

    #[inline]
    pub fn idx(&self, gx: i64, gy: i64) -> Option<usize> { node_index(gx, gy) }
}
