//! Running, publishing and re-reading one T1b window (lane TR).
//!
//! A published window is RAW (pre-blend): its posts over the computation
//! lattice (core + halo). Final region posts are the seam blend of raw
//! windows (`t1b_checks::blend_post`), a pure function of raw windows that
//! are themselves pure functions of their keys.

use crate::t1b::*;
use crate::t1b_lem::{self, LemConfig, Ledger, SeriesPoint};
use crate::*;
use serde::Serialize;
use serde_json::Value;
use std::path::Path;

pub const WINDOW_SCHEMA: &str = "vz.t1b-window/v1";

/// Lane TG G1: the shore_mm field note (definition and units; the rule is params.ls.tb.tc.tg.shore_rule in the recipe).
pub const SHORE_NOTE: &str = "lane TG (T1b 0.9.4, opt-in tg shore): signed horizontal distance in mm from the post to the sea shoreline, measured at sub-post precision; > 0 on land (every post not marine), < 0 on marine posts (fixed = 2); |value| clamped to the recipe's shore_cap_mm (a clamped value means at least that far). The shoreline is the zero contour (marching squares, linear crossings on 4-neighbour links) of phi = surface - sea level taken after the coast stage and before marine posts are published at sea level; posts the wave cut sends to sea carry the cut profile extended seaward, so a cliff foot sits at the platform-width distance, not on the 32 m post staircase. Lakes are not sea. Consumer use: interpolate shore_mm bilinearly; its zero isoline is the coastline and its value the distance from it (seaward negative)";

/// Lane TH H2: the uncut_mm field note (definition and units; the rule is params.ls.tb.tc.th.uncut_rule in the recipe).
pub const UNCUT_NOTE: &str = "lane TH (T1b 0.9.4, opt-in th uncut): on cliff-class posts, the surface height in mm immediately BEFORE the coast stage's wave cut (after volcano edifices and fans); -2147483648 (i32::MIN) on every other post (not cliff class: height_mm is the surface there). Cliff class = the posts on which the hard-rock cliff profile was evaluated (rock k_mult <= cliff_k_max at the post's height, edifice posts count as hard, cliff_profile on, within the wave-cut reach of the sea; not marine, lake, channel, T1a or fan no-cut), including the posts the cut sent to sea (now marine). T1b's cliff rule on these posts is z = min(U, P(s)) with s = shore_mm and P from the recipe's params.rivernet.cf: P(s) = talus_slope s for 0 <= s < talus_width, then talus_slope talus_width + wall_slope (s - talus_width). Consumer use at any spacing: bilinear U and s, rule = min(U, P(s)); blend to the published surface by the bilinear 0/1 class flag. Exact at a post only where the platform width varies slowly along the coast (measured in the lane TH report)";

/// Lane CLIMATE-RUNOFF: the runoff_weight field note.
pub const RUNOFF_WEIGHT_NOTE: &str = "lane CLIMATE-RUNOFF (climate key only): R(x) / r_ref at the post, R the climate runoff (climate_runoff.rs land_p_minus_e_floor_v1, land_only_bilinear_v1 between T0 cell centres), r_ref the T0 genesis runoff (A-RAIN-UNITS); every drainage_area_m2 of this window is the effective area Q / r_ref (free posts: own area a_post x runoff_weight accumulated); consumers multiply a post's own area by it (pf-detail climate_runoff=1)";

/// A raw window as published (the bytes are the contract).
#[derive(Clone)]
pub struct RawWindow {
    pub lat: WinLat,
    pub height_mm: Vec<i32>,
    pub fixed: Vec<u8>,
    pub drainage_area_m2: Vec<f32>,
    pub t0_bilinear_m: Vec<f32>,
    pub base_mm: Vec<i32>,
    pub roughness_mm: Vec<i32>,
    pub initial_mm: Vec<i32>,
    pub uplift_m_per_my: Vec<f32>,
    pub t0_cell: Vec<u32>,
    pub inflow_area_m2: Vec<f32>,
    /// Dominant exposed-lithology code on the published surface (lithology_mode 1 only; 255 on water).
    pub exposed_litho: Option<Vec<u8>>,
    /// Lane TD6 (T1b 0.9.3 with fx_publish only): coast/volcano over-code (t1b_coast::OVER_*) and volcano index (255 none).
    pub over: Option<Vec<u8>>,
    pub volcano_id: Option<Vec<u8>>,
    /// Lane TG G1 (tg shore only): signed sub-post shore distance, mm (see the field note).
    pub shore_mm: Option<Vec<i32>>,
    /// Lane TH H2 (th uncut only): height before the wave cut on cliff-class posts, mm (i32::MIN elsewhere; see the note).
    pub uncut_mm: Option<Vec<i32>>,
    /// Lane PV (prov=1 only): stable id of the exposed unit (see `strata_prov`), u32::MAX on water.
    pub prov_id: Option<Vec<u32>>,
    /// Lane CLIMATE-RUNOFF (climate key only): R / r_ref per post (the factor on the post's own drainage area).
    pub runoff_weight: Option<Vec<f32>>,
    /// Lane VE (ve only): eruption style of the exposed volcanic unit (t1b_ve::STYLE_NOTE) and its unit index (UNIT_NOTE).
    pub edifice_style: Option<Vec<u8>>,
    pub edifice_unit: Option<Vec<u16>>,
    /// Lane VE4 (ve4 front only): signed distance to the snout outline, mm (i32::MIN none; t1b_ve::FRONT_NOTE).
    pub front_mm: Option<Vec<i32>>,
    /// Lane GLACIAL-T1 (glacial key only): lowering of the post by the glacial window carve (U section + cirque), mm.
    pub gl_carve_mm: Option<Vec<i32>>,
    /// Lane STEP (ct only): cliff-and-talus regime, cliff fraction and massive talus top (t1b_ct notes).
    pub ct_regime: Option<Vec<u8>>,
    pub ct_cliff_permille: Option<Vec<u16>>,
    pub ct_top_mm: Option<Vec<i32>>,
    /// Lane STEP2 (ct_base only): the surface without the staircase shift (t1b_ct::BASE_NOTE).
    pub ct_base_mm: Option<Vec<i32>>,
    /// Lane TREAD2 (ct_tread only): tread height fraction and tread tangent (t1b_ct::TREAD_NOTE, TREAD_TAN_NOTE).
    pub ct_tread_permille: Option<Vec<u16>>,
    pub ct_tread_tan_permille: Option<Vec<u16>>,
    /// Lane PERIOD-ENDS (ct_ends only): period ends, gradient, kind and bed id (t1b_ct::ENDS_NOTE etc.).
    pub ct_ends: Option<CtEndsWin>,
    pub digest: String,
    pub recipe_digest: String,
}

/// Lane PERIOD-ENDS: the ct_ends window fields.
#[derive(Clone)]
pub struct CtEndsWin {
    pub s0_mm: Vec<i32>,
    pub s1_mm: Vec<i32>,
    pub s_dn_mm: Vec<i32>,
    pub s_up_mm: Vec<i32>,
    /// ct_knot_h1_mm, ct_knot_z1_mm, ct_knot_h2_mm, ct_knot_z2_mm.
    pub knot_mm: Vec<Vec<i32>>,
    pub grad_permille: Vec<u16>,
    pub kind: Vec<u8>,
    pub bed_id: Vec<u32>,
}

impl RawWindow {
    #[inline]
    pub fn z_m(&self, i: usize) -> f64 { self.height_mm[i] as f64 / 1000.0 }
    pub fn z_vec(&self) -> Vec<f64> { self.height_mm.iter().map(|v| *v as f64 / 1000.0).collect() }
}

#[derive(Clone, Debug, Serialize)]
pub struct FieldEntry { pub name: String, pub dtype: String, pub units: String, pub file: String, pub digest: String, pub note: String }

#[derive(Clone, Debug, Serialize)]
pub struct TileDigest { pub key: String, pub digest: String }

#[derive(Clone, Debug, Serialize)]
pub struct WindowHeader {
    pub schema: String,
    pub layer: String,
    pub label: String,
    pub simtime_label: String,
    pub closure_label: String,
    pub window: String,
    pub lattice: WinLat,
    pub core_posts: [i64; 4],
    pub recipe: ReliefRecipeV1,
    pub recipe_digest: String,
    pub prep: PrepStats,
    pub ledger: Ledger,
    pub series: Vec<SeriesPoint>,
    pub fields: Vec<FieldEntry>,
    pub raw_core_tile_digests: Vec<TileDigest>,
    pub preserved: Vec<String>,
    pub invented: Vec<String>,
    pub approximations: Vec<String>,
    /// Lane PV (prov=1 only; skipped when absent so every other window serializes byte-identically): the provenance
    /// rule, sources and the record of every prov_id present.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub provenance: Option<Value>,
}

/// Lane GLACIER-1 (gl_now=1 only): the present-ice window fields (t1b_glacial::window_ice), appended after the
/// header_for field list by `push_ice_fields`.
#[derive(Clone, Debug)]
pub struct IceWin { pub ice_now_mm: Vec<i32>, pub ice_surface_mm: Vec<i32>, pub ice_class: Vec<u8> }

/// Lane GLACIER-1: append the ice fields (same file / digest form as header_for).
pub fn push_ice_fields(h: &mut WindowHeader, blobs: &mut Vec<(String, Vec<u8>)>, ice: &IceWin) {
    for (name, dtype, units, note, bytes) in [
        ("ice_now_mm", "i32", "mm", crate::t1b_glacial::ICE_NOW_NOTE, le_i32(&ice.ice_now_mm)),
        ("ice_surface_mm", "i32", "mm", crate::t1b_glacial::ICE_SURFACE_NOTE, le_i32(&ice.ice_surface_mm)),
        ("ice_class", "u8", "code", crate::t1b_glacial::ICE_CLASS_NOTE, ice.ice_class.clone()),
    ] {
        let file = format!("fields/{name}.{dtype}");
        h.fields.push(FieldEntry { name: name.into(), dtype: dtype.into(), units: units.into(), file: file.clone(), digest: format!("blake3:{}", blake3::hash(&bytes).to_hex()), note: note.into() });
        blobs.push((file, bytes));
    }
}

/// Lane MS-CEIL (ceil only): the ceiling window fields, appended after the header_for field list by `push_ceil_fields`.
#[derive(Clone, Debug)]
pub struct CeilWin { pub ceil_m: Vec<f32>, pub surface_src: Vec<u8> }

/// Lane MS-CEIL: append ceil_m and surface_src (same file / digest form as header_for).
pub fn push_ceil_fields(h: &mut WindowHeader, blobs: &mut Vec<(String, Vec<u8>)>, c: &CeilWin) {
    for (name, dtype, units, note, bytes) in [
        ("ceil_m", "f32", "m", crate::t1b_rivernet::CEIL_NOTE, le_f32(&c.ceil_m)),
        ("surface_src", "u8", "code", crate::t1b_rivernet::CEIL_SRC_NOTE, c.surface_src.clone()),
    ] {
        let file = format!("fields/{name}.{dtype}");
        h.fields.push(FieldEntry { name: name.into(), dtype: dtype.into(), units: units.into(), file: file.clone(), digest: format!("blake3:{}", blake3::hash(&bytes).to_hex()), note: note.into() });
        blobs.push((file, bytes));
    }
}

pub struct WindowRun {
    pub raw: RawWindow,
    pub header: WindowHeader,
    pub prep: WindowInputs,
    pub blobs: Vec<(String, Vec<u8>)>,
    pub seconds_prepare: f64,
    pub seconds_relax: f64,
}

fn le_i32(v: &[i32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }
fn le_f32(v: &[f32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }
fn le_u32(v: &[u32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }
fn le_u16(v: &[u16]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }
fn mm(v: f64) -> R<i32> {
    let q = vz_tilekey::chart::quantize_mm(v * 1000.0);
    i32::try_from(q).map_err(|_| format!("height {v} m does not fit i32 mm"))
}

pub fn lem_config(p: &ReliefParamsV1, laws: &ReliefLawsV1) -> LemConfig {
    LemConfig {
        k_per_my: laws.stream_power_k_per_my,
        m: laws.stream_power_m,
        d_m2_per_my: laws.hillslope_diffusivity_m2_per_my * p.hillslope_d_ppm_of_t0.map_or(1.0, |v| v as f64 / 1e6),
        stability_limit: laws.diffusion_stability_limit,
        dt_my: p.dt_my(),
        steps: p.steps(),
        stats_every_steps: p.stats_every_steps(),
        dx_m: T1B_SPACING_M,
    }
}

/// Compute one window from scratch. Reads nothing about other windows.
pub fn run_window(inp: &ReliefInputs, recipe: &ReliefRecipeV1, recipe_digest: &str, rx: i64, ry: i64) -> R<WindowRun> {
    let t = std::time::Instant::now();
    let prep = prepare(inp, recipe, recipe_digest, rx, ry)?;
    let seconds_prepare = t.elapsed().as_secs_f64();
    let lat = prep.lat;
    let nn = lat.nodes();
    let core: Vec<bool> = (0..nn).map(|i| lat.in_core(i)).collect();
    let cfg = lem_config(&inp.params, &inp.extras.laws);
    let t = std::time::Instant::now();
    let res = match &prep.litho {
        None => t1b_lem::run(lat.n, &prep.fixed, &prep.z0, &prep.uplift_m_per_my, &prep.inflow_m2, &core, &cfg),
        Some(l) if inp.params.hillslope_mode == 1 => {
            let rc = t1b_lem::RoeringConfig { ratio_max: inp.params.creep_ratio_max_ppm.ok_or("creep_ratio_max_ppm")? as f64 / 1e6 };
            t1b_lem::run_full(lat.n, &prep.fixed, &prep.z0, &prep.uplift_m_per_my, &prep.inflow_m2, &core, &cfg, None, None, Some(l), Some(rc))
        }
        Some(l) => t1b_lem::run_litho(lat.n, &prep.fixed, &prep.z0, &prep.uplift_m_per_my, &prep.inflow_m2, &core, &cfg, None, None, Some(l)),
    };
    let seconds_relax = t.elapsed().as_secs_f64();
    let mut height_mm = Vec::with_capacity(nn);
    for i in 0..nn {
        // Fixed posts publish their fixed value exactly (identical in every window).
        height_mm.push(mm(if prep.fixed[i] == FREE { res.z[i] } else { prep.z0[i] })?);
    }
    let zq: Vec<f64> = height_mm.iter().map(|v| *v as f64 / 1000.0).collect();
    let r = t1b_lem::route(lat.n, &zq, &prep.fixed, T1B_SPACING_M);
    let acc = t1b_lem::accumulate(&r, T1B_SPACING_M * T1B_SPACING_M, Some(&prep.inflow_m2));
    let mut raw = RawWindow {
        lat,
        height_mm,
        fixed: prep.fixed.clone(),
        drainage_area_m2: acc.iter().map(|v| *v as f32).collect(),
        t0_bilinear_m: prep.bilinear_m.iter().map(|v| *v as f32).collect(),
        base_mm: prep.base_m.iter().map(|v| mm(*v)).collect::<R<_>>()?,
        roughness_mm: prep.roughness_m.iter().map(|v| mm(*v)).collect::<R<_>>()?,
        initial_mm: prep.z0.iter().map(|v| mm(*v)).collect::<R<_>>()?,
        uplift_m_per_my: prep.uplift_m_per_my.iter().map(|v| *v as f32).collect(),
        t0_cell: prep.t0_cell.clone(),
        inflow_area_m2: prep.inflow_m2.iter().map(|v| *v as f32).collect(),
        exposed_litho: prep.litho.as_ref().map(|l| {
            use crate::t1b_litho::ExposedLithology;
            (0..nn).map(|i| if prep.fixed[i] == FREE || prep.fixed[i] == GUIDE || prep.fixed[i] == CHANNEL { l.at(i, zq[i]).2 } else { 255 }).collect()
        }),
        over: None,
        volcano_id: None,
        shore_mm: None,
        uncut_mm: None,
        ct_regime: None, ct_cliff_permille: None, ct_top_mm: None, ct_base_mm: None, ct_tread_permille: None, ct_tread_tan_permille: None, ct_ends: None,
        prov_id: None,
        runoff_weight: None,
        edifice_style: None,
        edifice_unit: None,
        front_mm: None,
        gl_carve_mm: None,
        digest: String::new(),
        recipe_digest: recipe_digest.to_string(),
    };
    let (header, blobs) = header_for(&raw, recipe, recipe_digest, &prep, res.ledger, res.series)?;
    raw.digest = pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())?;
    Ok(WindowRun { raw, header, prep, blobs, seconds_prepare, seconds_relax })
}

fn field_list(raw: &RawWindow) -> Vec<(&'static str, &'static str, &'static str, &'static str, Vec<u8>)> {
    let mut v = vec![
        ("height_mm", "i32", "mm", "raw (pre-blend) T1b surface after the relaxation; fixed posts carry their fixed value", le_i32(&raw.height_mm)),
        ("fixed", "u8", "code", if raw.fixed.contains(&PLAIN) { "0 free, 1 T1a channel (fixed bed), 2 marine (sea level; includes estuary arms), 3 outer ring (guide-surface ring, 0.7.0), 4 lake (stage), 5 chart-wide tributary guide (fixed guide bed), 6 reclaimed coastal plain (fixed declared surface)" } else if raw.fixed.contains(&GUIDE) { "0 free, 1 T1a channel (fixed bed), 2 marine (sea level), 3 outer ring (initial surface), 4 lake (stage), 5 chart-wide tributary guide (fixed guide bed)" } else { "0 free, 1 T1a channel (fixed bed), 2 marine (sea level), 3 outer ring (initial surface), 4 lake (stage)" }, raw.fixed.clone()),
        ("drainage_area_m2", "f32", "m2", "window-local drainage area on the published surface (priority flood from fixed posts); no inflow across the ring", le_f32(&raw.drainage_area_m2)),
        ("t0_bilinear_m", "f32", "m", "T0 elevation bilinear between L1 cell centres (the T0 crop)", le_f32(&raw.t0_bilinear_m)),
        ("base_mm", "i32", "mm", "deterministic base: T0 bilinear floored at sea level + 1 m on land; sea level / stage on water cells", le_i32(&raw.base_mm)),
        ("roughness_mm", "i32", "mm", "keyed initial roughness (sub-resolution-closure-v1); 0 on water cells", le_i32(&raw.roughness_mm)),
        ("initial_mm", "i32", "mm", "initial surface: base + roughness on land, fixed values on fixed posts", le_i32(&raw.initial_mm)),
        ("uplift_m_per_my", "f32", "m/My", "relaxation source term U (t1b-relaxation-not-history; not T0's uplift field)", le_f32(&raw.uplift_m_per_my)),
        ("t0_cell", "u32", "id", "T0 L1 element owning the post (chart inverse, floor)", le_u32(&raw.t0_cell)),
        ("inflow_area_m2", "f32", "m2", "T1a drainage area entering across the closed ring (closed_ring_t1a_inflow_v1); 0 elsewhere", le_f32(&raw.inflow_area_m2)),
    ];
    if let Some(e) = &raw.exposed_litho {
        v.push(("exposed_litho", "u8", "code", "dominant exposed lithology on the published surface (t1b_litho codes, design review LithologyV1: 0 crystalline basement, 1 granite, 2 metamorphic, 3 basalt/volcanic, 4 carbonate, 5 shale, 6 sandstone, 7 conglomerate, 8 evaporite, 9 coal-bearing, 10 unconsolidated, 11 ophiolite; 255 water); stand-in cover_over_basement_v1", e.clone()));
    }
    if let Some(o) = &raw.over {
        v.push(("over", "u8", "code", "lane TD (T1b 0.9.3): coast / volcano over-code of the post (0 none, 1 volcanic edifice: cone, flows, gullies, crater, or relief buried by the cone; 2 fan deposit; 3 beach; 4 wave cut: platform, cliff or soft-rock bluff)", o.clone()));
    }
    if let Some(o) = &raw.volcano_id {
        v.push(("volcano_id", "u8", "index", "lane TD (T1b 0.9.3): on edifice posts (over = 1) the index into the header's prep.net.volcanoes list of the cone that built or buried the post; 255 elsewhere", o.clone()));
    }
    if let Some(s) = &raw.shore_mm {
        v.push(("shore_mm", "i32", "mm", SHORE_NOTE, le_i32(s)));
    }
    if let Some(u) = &raw.uncut_mm {
        v.push(("uncut_mm", "i32", "mm", UNCUT_NOTE, le_i32(u)));
    }
    if let Some(p) = &raw.prov_id {
        v.push(("prov_id", "u32", "id", crate::strata_prov::ID_NOTE, le_u32(p)));
    }
    if let Some(w) = &raw.runoff_weight {
        v.push(("runoff_weight", "f32", "1", RUNOFF_WEIGHT_NOTE, le_f32(w)));
    }
    if let Some(e) = &raw.edifice_style {
        v.push(("edifice_style", "u8", "code", crate::t1b_ve::STYLE_NOTE, e.clone()));
    }
    if let Some(e) = &raw.edifice_unit {
        v.push(("edifice_unit", "u16", "index", crate::t1b_ve::UNIT_NOTE, le_u16(e)));
    }
    if let Some(f) = &raw.front_mm {
        v.push(("front_mm", "i32", "mm", crate::t1b_ve::FRONT_NOTE, le_i32(f)));
    }
    if let Some(r) = &raw.ct_regime {
        v.push(("ct_regime", "u8", "code", crate::t1b_ct::REGIME_NOTE, r.clone()));
    }
    if let Some(g) = &raw.gl_carve_mm {
        v.push(("gl_carve_mm", "i32", "mm", crate::t1b_glacial::CARVE_NOTE, le_i32(g)));
    }
    if let Some(c) = &raw.ct_cliff_permille {
        v.push(("ct_cliff_permille", "u16", "permille", crate::t1b_ct::CLIFF_NOTE, le_u16(c)));
    }
    if let Some(t) = &raw.ct_top_mm {
        v.push(("ct_top_mm", "i32", "mm", crate::t1b_ct::TOP_NOTE, le_i32(t)));
    }
    if let Some(b) = &raw.ct_base_mm {
        v.push(("ct_base_mm", "i32", "mm", crate::t1b_ct::BASE_NOTE, le_i32(b)));
    }
    if let Some(t) = &raw.ct_tread_permille {
        v.push(("ct_tread_permille", "u16", "permille", crate::t1b_ct::TREAD_NOTE, le_u16(t)));
    }
    if let Some(t) = &raw.ct_tread_tan_permille {
        v.push(("ct_tread_tan_permille", "u16", "permille", crate::t1b_ct::TREAD_TAN_NOTE, le_u16(t)));
    }
    if let Some(e) = &raw.ct_ends {
        v.push(("ct_s0_mm", "i32", "mm", crate::t1b_ct::ENDS_NOTE, le_i32(&e.s0_mm)));
        v.push(("ct_s1_mm", "i32", "mm", crate::t1b_ct::ENDS_NOTE, le_i32(&e.s1_mm)));
        v.push(("ct_s_dn_mm", "i32", "mm", crate::t1b_ct::ENDS_NOTE, le_i32(&e.s_dn_mm)));
        v.push(("ct_s_up_mm", "i32", "mm", crate::t1b_ct::ENDS_NOTE, le_i32(&e.s_up_mm)));
        for (k, nm) in ["ct_knot_h1_mm", "ct_knot_z1_mm", "ct_knot_h2_mm", "ct_knot_z2_mm"].iter().enumerate() { v.push((*nm, "i32", "mm", crate::t1b_ct::KNOT_NOTE, le_i32(&e.knot_mm[k]))); }
        v.push(("ct_grad_permille", "u16", "permille", crate::t1b_ct::GRAD_NOTE, le_u16(&e.grad_permille)));
        v.push(("ct_period_kind", "u8", "code", crate::t1b_ct::KIND_NOTE, e.kind.clone()));
        v.push(("ct_bed_id", "u32", "id", crate::t1b_ct::BED_ID_NOTE, le_u32(&e.bed_id)));
    }
    v
}

pub(crate) fn header_for(raw: &RawWindow, recipe: &ReliefRecipeV1, recipe_digest: &str, prep: &WindowInputs, ledger: Ledger, series: Vec<SeriesPoint>)
    -> R<(WindowHeader, Vec<(String, Vec<u8>)>)> {
    let lat = raw.lat;
    let mut fields = Vec::new();
    let mut blobs = Vec::new();
    for (name, dtype, units, note, bytes) in field_list(raw) {
        let file = format!("fields/{name}.{dtype}");
        fields.push(FieldEntry { name: name.into(), dtype: dtype.into(), units: units.into(), file: file.clone(),
            digest: format!("blake3:{}", blake3::hash(&bytes).to_hex()), note: note.into() });
        blobs.push((file, bytes));
    }
    let (x0, x1, y0, y1) = lat.core();
    let mut tiles = Vec::new();
    for ty in (y0 / TILE_CELLS)..(y1 / TILE_CELLS) { for tx in (x0 / TILE_CELLS)..(x1 / TILE_CELLS) {
        let mut h = blake3::Hasher::new();
        for gy in ty * TILE_CELLS..=ty * TILE_CELLS + TILE_CELLS { for gx in tx * TILE_CELLS..=tx * TILE_CELLS + TILE_CELLS {
            let i = lat.idx(gx, gy).ok_or("core tile post outside the window")?;
            h.update(&raw.height_mm[i].to_le_bytes());
        } }
        tiles.push(TileDigest { key: format!("{T1B_LOD}.{tx}.{ty}"), digest: format!("blake3:{}", h.finalize().to_hex()) });
    } }
    let mut preserved = vec![
        "T1a channel==1 links as fixed posts at the T1a bed (F1): same planform, same levels, same outlet".to_string(),
        "T0 marine and lake cells as fixed water posts (T0 cell classification, T1a's coast reading)".to_string(),
        "T0 erosion law and its constants (K, m, n, D, stability limit) read from the T0 terrain header".to_string(),
    ];
    if recipe.params.guide_mode == 1 {
        preserved.push("chart-wide tributary guides (t1b.tributary_guides) as fixed posts at the guide beds: the same layout in every window and every roughness realization".to_string());
    }
    let mut h = WindowHeader {
        schema: WINDOW_SCHEMA.into(), layer: T1B_LAYER.into(), label: T1B_LABEL.into(), simtime_label: SIMTIME_LABEL.into(),
        closure_label: format!("{CLOSURE_LABEL}: initial roughness rms = min({} ppm x T0 subgrid relief, {} mm), keyed uniform, white at 32 m, realization {}",
            recipe.params.roughness_rms_fraction_ppm, recipe.params.roughness_rms_cap_mm, recipe.params.roughness_realization),
        window: lat.key(), lattice: lat, core_posts: [x0, x1, y0, y1],
        recipe: recipe.clone(), recipe_digest: recipe_digest.to_string(), prep: prep.stats.clone(), ledger, series, fields,
        raw_core_tile_digests: tiles,
        preserved,
        invented: vec![
            "relief between and around the fixed channels, by the relaxation process only".into(),
            "keyed initial roughness (sub-resolution-closure-v1), bounded per ADR 043 section 2 condition 7; separable (zero amplitude = deterministic base)".into(),
        ],
        approximations: vec![
            "SimTime is t1b-relaxation-not-history: a numerical device, not world time, not T0's history; not an age".into(),
            "relaxation source term U = gain K R_sub / Lambda (T0's subgrid-closure rate), not T0's uplift field (reported in prep)".into(),
            "closed outer ring: no diffusive flux or outflow across it; T1a hillslope inflow spread over 512 m edge strips (halo absorbs; seams blended)".into(),
            "planar post area (point scale <= 2.1e-4 over the chart, not applied)".into(),
            "no climate input (ADR 043 OPEN 6); no glacial or coastal processes (OPEN 5)".into(),
            "fidelity tier Preview; pre-signoff substrate".into(),
        ],
        provenance: None,
    };
    if recipe.params.lithology_mode == 1 {
        h.invented.push("exposed lithology by the cover_over_basement_v1 stand-in (no draws): cover where the surface stays above T0 surface minus T0 sediment thickness, basement (T0 rock_family class) below; per-post critical slope and K multiplier from the declared table".into());
        h.approximations.push("rock class from T0 cells (about 70 km): no layer stack, no structural fabric (dip, strike, foliation, cleavage); class contacts are land-weighted bilinear gradients, not mapped contacts (SEED-STRATIGRAPHY queued)".into());
        if recipe.params.hillslope_mode != 1 {
            h.approximations.push("threshold-slope failure removes material instantly and exports it (no landslide deposits, no runout)".into());
        }
    }
    let p = &recipe.params;
    if p.ring_mode == 1 {
        h.approximations.push("outer ring fixed at the guide layer's relaxed 512 m surface (bilinear): the same boundary function in every window; ring posts are outlets; the halo absorbs the 512 m vs 32 m mismatch".into());
    }
    if p.hillslope_mode == 1 {
        h.approximations.push("nonlinear creep with the T0 header D (legacy A-HILLSLOPE-D); x/y operator splitting (fully implicit along each sweep, cross gradient lagged); factor tangent-extended beyond creep_ratio_max x Sc; no landslide deposits".into());
    }
    if p.guide_curve_passes > 0 {
        h.invented.push("guide valley planform between break nodes: Chaikin corner cutting of the 512 m D8 guide polylines (geometry only, no draws)".into());
    }
    if p.coast_mode >= 1 {
        h.invented.push("coastline at 32 m: the T0 bilinear zero contour (no sub-T0 coastal process: no wave cut, no deltas, no drowned valleys)".into());
    }
    if p.coast_mode == 2 {
        h.invented.push("reclaimed coastal plain on T0 marine-cell posts above the zero contour (zero relaxation uplift, declared plain ratio) and river mouths extended down the T0 bilinear surface to the new shore".into());
    }
    Ok((h, blobs))
}

pub fn publish(run: &WindowRun, out: &Path) -> R<()> {
    std::fs::create_dir_all(out.join("fields")).map_err(|e| e.to_string())?;
    for (file, bytes) in &run.blobs { std::fs::write(out.join(file), bytes).map_err(|e| format!("{file}: {e}"))?; }
    let doc = serde_json::json!({ "header": run.header, "digest": run.raw.digest });
    std::fs::write(out.join("window.json"), serde_json::to_vec_pretty(&doc).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    Ok(())
}

fn i32s(b: &[u8]) -> Vec<i32> { b.chunks_exact(4).map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
fn f32s(b: &[u8]) -> Vec<f32> { b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }
fn u32s(b: &[u8]) -> Vec<u32> { b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() }

/// Re-read a published window, verifying the header digest and every field digest.
pub fn load(dir: &Path) -> R<(RawWindow, Value)> {
    let bytes = std::fs::read(dir.join("window.json")).map_err(|e| format!("window.json in {}: {e}", dir.display()))?;
    let v: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    let h = v.get("header").ok_or("no header")?.clone();
    let digest = v.get("digest").and_then(Value::as_str).ok_or("no digest")?.to_string();
    let got = pf_contracts::artifact::digest_canonical(&h).map_err(|e| e.to_string())?;
    if got != digest { return Err(format!("window header digest {got} differs from stated {digest}; refusing")); }
    let mut blobs = std::collections::BTreeMap::new();
    for f in h.get("fields").and_then(Value::as_array).ok_or("fields")? {
        let name = f.get("name").and_then(Value::as_str).ok_or("name")?.to_string();
        let file = f.get("file").and_then(Value::as_str).ok_or("file")?;
        let b = std::fs::read(dir.join(file)).map_err(|e| format!("{file}: {e}"))?;
        let d = format!("blake3:{}", blake3::hash(&b).to_hex());
        if Some(d.as_str()) != f.get("digest").and_then(Value::as_str) { return Err(format!("field {name} digest mismatch; refusing")); }
        blobs.insert(name, b);
    }
    let lat: WinLat = {
        let l = h.get("lattice").ok_or("lattice")?;
        let gi = |k: &str| l.get(k).and_then(Value::as_i64).ok_or_else(|| format!("lattice.{k}"));
        WinLat { rx: gi("rx")?, ry: gi("ry")?, gx0: gi("gx0")?, gy0: gi("gy0")?, n: gi("n")? as usize, halo_posts: gi("halo_posts")? }
    };
    let b = |n: &str| blobs.get(n).cloned().ok_or_else(|| format!("field {n} missing"));
    let raw = RawWindow {
        lat,
        height_mm: i32s(&b("height_mm")?),
        fixed: b("fixed")?,
        drainage_area_m2: f32s(&b("drainage_area_m2")?),
        t0_bilinear_m: f32s(&b("t0_bilinear_m")?),
        base_mm: i32s(&b("base_mm")?),
        roughness_mm: i32s(&b("roughness_mm")?),
        initial_mm: i32s(&b("initial_mm")?),
        uplift_m_per_my: f32s(&b("uplift_m_per_my")?),
        t0_cell: u32s(&b("t0_cell")?),
        inflow_area_m2: f32s(&b("inflow_area_m2")?),
        exposed_litho: blobs.get("exposed_litho").cloned(),
        over: blobs.get("over").cloned(),
        volcano_id: blobs.get("volcano_id").cloned(),
        shore_mm: blobs.get("shore_mm").map(|b| i32s(b)),
        uncut_mm: blobs.get("uncut_mm").map(|b| i32s(b)),
        prov_id: blobs.get("prov_id").map(|b| u32s(b)),
        runoff_weight: blobs.get("runoff_weight").map(|b| f32s(b)),
        edifice_style: blobs.get("edifice_style").cloned(),
        edifice_unit: blobs.get("edifice_unit").map(|b| b.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect()),
        front_mm: blobs.get("front_mm").map(|b| i32s(b)),
        gl_carve_mm: blobs.get("gl_carve_mm").map(|b| i32s(b)),
        ct_regime: blobs.get("ct_regime").cloned(),
        ct_cliff_permille: blobs.get("ct_cliff_permille").map(|b| b.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect()),
        ct_top_mm: blobs.get("ct_top_mm").map(|b| i32s(b)),
        ct_base_mm: blobs.get("ct_base_mm").map(|b| i32s(b)),
        ct_tread_permille: blobs.get("ct_tread_permille").map(|b| b.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect()),
        ct_tread_tan_permille: blobs.get("ct_tread_tan_permille").map(|b| b.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect()),
        ct_ends: match (blobs.get("ct_s0_mm"), blobs.get("ct_s1_mm"), blobs.get("ct_s_dn_mm"), blobs.get("ct_s_up_mm"), blobs.get("ct_grad_permille"), blobs.get("ct_period_kind"), blobs.get("ct_bed_id")) {
            (Some(a), Some(b), Some(c), Some(d), Some(g), Some(k), Some(bi)) => Some(CtEndsWin { s0_mm: i32s(a), s1_mm: i32s(b), s_dn_mm: i32s(c), s_up_mm: i32s(d),
                knot_mm: ["ct_knot_h1_mm", "ct_knot_z1_mm", "ct_knot_h2_mm", "ct_knot_z2_mm"].iter().map(|k| blobs.get(*k).map(|b| i32s(b)).unwrap_or_default()).collect(),
                grad_permille: g.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect(), kind: k.clone(), bed_id: u32s(bi) }),
            _ => None,
        },
        digest,
        recipe_digest: h.get("recipe_digest").and_then(Value::as_str).ok_or("recipe_digest")?.to_string(),
    };
    if raw.height_mm.len() != lat.nodes() { return Err("height field size".into()); }
    Ok((raw, h))
}
