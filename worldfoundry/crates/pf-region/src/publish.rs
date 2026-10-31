//! Writing the skeleton: raw little-endian node fields plus a header whose
//! digest is blake3 over canonical JSON (ADR 044 section 1, OPEN 3 ruling).

use crate::skeleton::{BuildStats, CellInfo, Segment, Skeleton};
use crate::*;
use serde::Serialize;
use std::path::Path;

pub const SCHEMA: &str = "vz.t1a-skeleton/v1";

#[derive(Clone, Debug, Serialize)]
pub struct FieldEntry {
    pub name: String,
    pub dtype: String,
    pub units: String,
    pub file: String,
    pub digest: String,
    pub note: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct TileDigest { pub key: String, pub digest: String }

#[derive(Clone, Debug, Serialize)]
pub struct Header {
    pub schema: String,
    pub layer: String,
    pub label: String,
    pub recipe: SkeletonRecipeV1,
    pub recipe_digest: String,
    pub site_in_chart_mm_recomputed: [i64; 2],
    pub site_in_chart_mm_record: [i64; 2],
    pub walk_tile: String,
    pub region: String,
    pub sea_level_m: f64,
    pub river_threshold_m3_s: f64,
    pub fields: Vec<FieldEntry>,
    pub tile_digests: Vec<TileDigest>,
    pub segments: Vec<Segment>,
    pub cells: Vec<CellInfo>,
    pub stats: BuildStats,
    pub preserved: Vec<String>,
    pub invented: Vec<String>,
    pub approximations: Vec<String>,
    /// Lane CLIMATE-RUNOFF (key on only): digest of the same skeleton without the key (T1b's draw twin) and the climate
    /// lake budgets of the T0 lakes touching the computation lattice.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub climate_twin_digest: Option<String>,
    /// Metadata-free identity of the climate-free twin, for downstream draws only.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub physical_climate_twin_digest: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub climate_lakes: Option<Vec<crate::climate_runoff::LakeBudget>>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Published { pub header: Header, pub digest: String }

fn le_u32(v: &[u32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }
fn le_f32(v: &[f64]) -> Vec<u8> { v.iter().flat_map(|x| (*x as f32).to_le_bytes()).collect() }

/// The node fields, in a fixed order: (name, dtype, units, note, bytes, element size).
pub fn field_bytes(sk: &Skeleton) -> Vec<(&'static str, &'static str, &'static str, &'static str, Vec<u8>, usize)> {
    let thr = sk.threshold_m3_s;
    let channel: Vec<u8> = (0..NODES).map(|n| match sk.kind[n] {
        1 if sk.discharge_m3_s[n] >= thr => 1,
        0 | 4 if sk.discharge_m3_s[n] >= thr => 2,
        _ => 0,
    }).collect();
    let fd: Vec<f64> = sk.flow_dist_m.iter().map(|v| if v.is_finite() { *v } else { 0.0 }).collect();
    let mut v = vec![
        ("t0_cell", "u32", "id", "T0 L-grid element owning the node (chart inverse, floor)", le_u32(&sk.t0_cell), 4),
        ("kind", "u8", "code", "0 hillslope, 1 valley axis (corridor path), 2 marine, 3 mouth, 4 lake hillslope", sk.kind.clone(), 1),
        ("channel", "u8", "code", "1 axis node with T1a discharge >= threshold (a skeleton channel), 2 hillslope node >= threshold (F2 violation), 0 otherwise", channel, 1),
        ("receiver", "u32", "node", "downstream node; 4294967295 terminal, 4294967294 off-chart", le_u32(&sk.receiver), 4),
        ("segment", "u32", "id", "axis nodes: owning segment id; others: the segment their flow joins (4294967295 none)",
            le_u32(&(0..NODES).map(|n| if sk.kind[n] == 1 { sk.owner[n] } else { sk.join_segment[n] }).collect::<Vec<_>>()), 4),
        ("join_node", "u32", "node", "hillslope nodes: the axis node their flow path joins (4294967295 none)", le_u32(&sk.join_node), 4),
        ("discharge_m3_s", "f32", "m3/s", if sk.climate_runoff_m_per_yr.is_some() { "accumulated T1a discharge (climate runoff R(x) x area, climate_runoff.rs, plus chart-entry inflow from the climate T0 discharge)" } else { "accumulated T1a discharge (T0 uniform runoff x area, plus chart-entry inflow)" }, le_f32(&sk.discharge_m3_s), 4),
        ("trunk_level_m", "f32", "m", "axis: carved bed level; hillslope: bed level at its join node; marine: sea level", le_f32(&sk.trunk_level_m), 4),
        ("flow_dist_to_axis_m", "f32", "m", "hillslope: D8 path length to the join node; 0 on axes, marine and off-chart exits", le_f32(&fd), 4),
        ("outlet_basin", "u32", "id", "T0 basin id of the T0 terminal the node drains to", le_u32(&sk.basin), 4),
        ("t0_reach", "u32", "id", "T0 reach id the node maps to (its axis segment's start cell reach)", le_u32(&sk.t0_reach), 4),
        ("area_m2", "f32", "m2", "true area represented by the post (trapezoid weight / k^2)", le_f32(&sk.area_m2), 4),
    ];
    // Lane CLIMATE-RUNOFF (key on only): the climate runoff at each node.
    if let Some(r) = &sk.climate_runoff_m_per_yr {
        v.push(("runoff_m_per_yr", "f32", "m/yr", "lane CLIMATE-RUNOFF: climate runoff R at the node (land_p_minus_e_floor_v1, land_only_bilinear_v1; climate_runoff.rs)", le_f32(r), 4));
    }
    v
}

pub fn header(sk: &Skeleton) -> R<(Header, Vec<(String, Vec<u8>)>)> {
    let mut fields = Vec::new();
    let mut blobs = Vec::new();
    let mut per_field = Vec::new();
    for (name, dtype, units, note, bytes, size) in field_bytes(sk) {
        if bytes.len() != NODES * size { return Err(format!("field {name}: {} bytes", bytes.len())); }
        let file = format!("fields/{name}.{dtype}");
        fields.push(FieldEntry { name: name.into(), dtype: dtype.into(), units: units.into(), file: file.clone(),
            digest: format!("blake3:{}", blake3::hash(&bytes).to_hex()), note: note.into() });
        per_field.push((bytes.clone(), size));
        blobs.push((file, bytes));
    }
    // lod-8 tile digests (129 x 129 edge-inclusive posts per tile, all fields in order).
    let mut tiles = Vec::new();
    let (lo, hi) = vz_tilekey::key_range(LOD);
    for iy in lo..hi { for ix in lo..hi {
        let mut h = blake3::Hasher::new();
        for gy in iy * 128..=iy * 128 + 128 { for gx in ix * 128..=ix * 128 + 128 {
            let n = node_index(gx, gy).ok_or("tile post outside lattice")?;
            for (bytes, size) in &per_field { h.update(&bytes[n * size..(n + 1) * size]); }
        } }
        tiles.push(TileDigest { key: format!("{LOD}.{ix}.{iy}"), digest: format!("blake3:{}", h.finalize().to_hex()) });
    } }
    let (wl, wx, wy) = sk.site.walk_tile;
    let (rl, rx, ry) = sk.site.region;
    let h = Header {
        schema: SCHEMA.into(), layer: LAYER.into(), label: LABEL.into(), recipe: sk.recipe.clone(), recipe_digest: sk.recipe_digest.clone(),
        site_in_chart_mm_recomputed: sk.stats.site_in_chart_mm_recomputed, site_in_chart_mm_record: sk.site.site_in_chart_mm,
        walk_tile: format!("{wl}.{wx}.{wy}"), region: format!("{rl}.{rx}.{ry}"),
        sea_level_m: sk.sea_level_m, river_threshold_m3_s: sk.threshold_m3_s,
        fields, tile_digests: tiles, segments: sk.segments.clone(), cells: sk.cells.values().cloned().collect(), stats: sk.stats.clone(),
        preserved: vec![
            "every T0 cell-to-receiver flow edge touching the chart as a connected valley axis to the same T0 receiver, hence the same T0 reach outlet basin (F1)".into(),
            "T0 cell partition of drainage: each node drains through its own T0 cell's outlet (F3)".into(),
            "T0 trunk level at every T0 cell centre (bed = T0 mean - subgrid relief at the junction, before carving)".into(),
            "T0 sea level at mouths; T0 runoff assumption for discharge".into(),
        ],
        invented: vec![
            "keyed planform offsets inside the corridor (vz.keyed-rng/v1, purpose t1a.planform.offset), |offset| <= amplitude < corridor half width".into(),
            "headwater extension axes (opposite_half_segment_v1) and their bed rise to the T0 cell mean".into(),
            "bed profile between T0 junction levels (linear, carved monotone)".into(),
            "hillslope flow paths (declared shortest-path rule; a partition, not a landform)".into(),
        ],
        approximations: vec![
            "T0 resolves one value per ~70 km cell plus the reach line; inside a cell nothing here is resolved by T0".into(),
            "chart-entry inflow = T0 discharge of the start cell minus the runoff of its in-chart area (edge effect, labelled)".into(),
            "subgrid relief read from lane SP's driver sidecar (not a published artifact field); its blake3 is pinned as an upstream".into(),
            if sk.climate_runoff_m_per_yr.is_some() { "lane CLIMATE-RUNOFF: discharge carries the climate runoff field of the seed-climate-v1 artifact (upstream t0.climate); T0's river network, reaches and planform are T0's (computed with uniform runoff), so climate changes discharge, not topology; dry lakes are reported, not removed".into() } else { "no climate input (ADR 043 OPEN 6); discharge uses T0's uniform runoff".into() },
            "fidelity tier Preview; pre-signoff substrate".into(),
        ],
        climate_twin_digest: sk.climate_twin_digest.clone(),
        physical_climate_twin_digest: sk.physical_climate_twin_digest.clone(),
        climate_lakes: if sk.climate_runoff_m_per_yr.is_some() { Some(sk.climate_lakes.clone()) } else { None },
    };
    Ok((h, blobs))
}

pub fn publish(sk: &Skeleton, out: &Path) -> R<Published> {
    let (h, blobs) = header(sk)?;
    let digest = pf_contracts::artifact::digest_canonical(&h).map_err(|e| e.to_string())?;
    std::fs::create_dir_all(out.join("fields")).map_err(|e| e.to_string())?;
    for (file, bytes) in &blobs { std::fs::write(out.join(file), bytes).map_err(|e| format!("{file}: {e}"))?; }
    let p = Published { header: h, digest };
    std::fs::write(out.join("skeleton.json"), serde_json::to_vec_pretty(&p).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    Ok(p)
}

/// The digest without writing anything (for determinism tests).
pub fn digest(sk: &Skeleton) -> R<String> {
    let (h, _) = header(sk)?;
    pf_contracts::artifact::digest_canonical(&h).map_err(|e| e.to_string())
}
