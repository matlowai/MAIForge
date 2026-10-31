//! Chart-wide tributary guides for T1b (lane TW, 2026-10-12; ADR 043
//! amendment "chart-wide tributary network", design 2026-10-12).
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon. SimTime
//! label `t1b-relaxation-not-history`.
//!
//! Why: lane TR measured that each T1b window grew its own 5-15 km tributary
//! network, which depended on the window edge and the roughness draw (seam
//! rms 759 m; F12 mean 57 m). This layer fixes that network ONCE PER CHART:
//! one relaxation of the same law (`t1b_lem`, stream power n = 1 implicit,
//! linear diffusion, constants from the T0 terrain header) on the T1a lod-8
//! lattice itself (512 m posts, 513 x 513, the whole chart), with the T1a
//! channels and T0 water fixed exactly as T1b fixes them. Its free nodes whose
//! drainage area reaches a declared threshold are GUIDES: every T1b window
//! fixes the posts on each guide's D8 link at the guide bed (the priority-flood
//! fill of the relaxed surface, monotone along receivers), so the tributary
//! layout no longer depends on which window is computed or on the window's
//! roughness draw. T1b then relaxes only the hillslopes between guides.
//!
//! Inputs: T0 (closure rate, water classes), the T1a skeleton (verified
//! view), the site seed. One keyed draw: the initial roughness per lod-8 node
//! (purpose `t1b.guides.initial.roughness`, layer = this recipe's digest), so
//! the guide layout is a pure function of (seed, T0, T1a, recipe).

use crate::t0::{CellWater, SiteRecord, T0};
use crate::t1b::{bilinear, closure_rate_m_per_my, T0Extras, ReliefLawsV1, CHANNEL, FREE, LAKE, MARINE, SIMTIME_LABEL, CLOSURE_LABEL};
use crate::t1b_lem::{self, Ledger, LemConfig, SeriesPoint};
use crate::t1b_view::SkeletonView;
use crate::*;
use serde::Serialize;
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;
use vz_keyed_rng::Scope;
use vz_tilekey::chart::{Chart, ChartDeclV1, RAD2DEG};

pub const GUIDES_GENERATOR_NAME: &str = "t1b-tributary-guides";
pub const GUIDES_GENERATOR_VERSION: &str = "0.3.0";
/// With exposed lithology and threshold-slope failure (lithology_mode 1; lane T3, 2026-10-13).
pub const GUIDES_GENERATOR_VERSION_LITHO: &str = "0.4.0";
pub const GUIDES_RECIPE_SCHEMA: &str = "vz.t1b-guides-recipe/v1";
pub const GUIDES_SCHEMA: &str = "vz.t1b-guides/v1";
pub const GUIDES_LAYER: &str = "t1b.tributary_guides";
pub const PURPOSE_GUIDE_ROUGHNESS: &str = "t1b.guides.initial.roughness";
pub const GUIDES_LABEL: &str = "pre-signoff; T0 substrate R2 not met; fidelity tier Preview; derived, never canon; SimTime t1b-relaxation-not-history";
/// Node class in the published `fixed` field for a guide node (free in the relaxation).
pub const GUIDE_NODE: u8 = 5;

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct GuideParamsV1 {
    /// Relaxation SimTime [default 40 My, same as T1b 0.4.0].
    pub sim_time_kyr: i64,
    /// Step [default 20 kyr, same as T1b].
    pub dt_kyr: i64,
    pub stats_every_kyr: i64,
    /// Source-term gain, permille [default 1000 since 0.3.0: chosen by lane TW from runs at 440 and 1000 (guide beds in
    /// 6.0.0 718 vs 810-814 m above a 710 m channel); not fitted; F6 still short, needs a ruling].
    pub uplift_gain_permille: i64,
    /// Initial roughness rms fraction of R_sub, ppm, and cap [default 5000 ppm, 5 m: T1b's rule at 512 m].
    pub roughness_rms_fraction_ppm: i64,
    pub roughness_rms_cap_mm: i64,
    pub base_floor_above_sea_mm: i64,
    /// Guide area threshold [default 4 km2 since 0.3.0, project assumption (lane TW): 2 km2 (guides_c2) gave a dense
    /// grid-aligned maze; 8 km2 gave larger seam errors; with 4 km2 a 4096 m T1b halo passes all four seams].
    pub guide_area_threshold_m2: i64,
    pub lattice_rule: String,
    pub fixed_rule: String,
    pub uplift_rule: String,
    pub base_rule: String,
    pub roughness_rule: String,
    pub routing_rule: String,
    pub guide_rule: String,
    /// 1 = t0_cell_confined_routing_v1 [default since 0.2.0]: flow stays in the node's T0 cell or enters a fixed
    /// node of the cell's T0 receiver cell (T1a's F3 rule), so guides never carry one T0 cell's area into
    /// another's (0.1.0, unconfined, did: F3 failed in 18 of 21 cells, lane TW run guides_g2).
    pub confine_mode: u32,
    /// 1 = t1a_trunk_floor_v1 [default since 0.3.0]: a free node never sits below the T1a trunk level (the bed of
    /// the channel T1a drains it to; T0's closure reads it as the cell's lowest ground). Without it (0.2.0, run
    /// guides_c2) the 512 m network bypassed the T1a channel in parallel valleys down to the cell's outlet, so guide
    /// beds in window 6.0.0 averaged 345 m against a channel at 710 m (channel perched; T1b mean 956 m vs 2228 m).
    pub floor_mode: u32,
    /// 0 = off (0.3.0, serialized exactly as before); 1 = exposed lithology with per-node K and threshold-slope
    /// failure, the same law as T1b 0.6.0 (guides 0.4.0, lane T3).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32")]
    pub lithology_mode: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lithology: Option<crate::t1b_litho::LithologyV1>,
}

impl Default for GuideParamsV1 {
    fn default() -> Self {
        GuideParamsV1 {
            sim_time_kyr: 40_000,
            dt_kyr: 20,
            stats_every_kyr: 1_000,
            uplift_gain_permille: 1000,
            roughness_rms_fraction_ppm: 5_000,
            roughness_rms_cap_mm: 5_000,
            base_floor_above_sea_mm: 1_000,
            guide_area_threshold_m2: 4_000_000,
            lattice_rule: "t1a_lod8_chart_lattice_v1: the T1a lattice itself (lod-8, 512 m, 513 x 513 edge-inclusive posts over the whole chart); closed edge (no flux or water across it)".into(),
            fixed_rule: "t1a_channels_t0_water_v1: node fixed at sea level if its T0 cell (T1a t0_cell) is marine, at the stage if lake, else at the T1a trunk level if T1a kind == 1 and channel == 1 or kind == 3 (mouth); all other nodes free".into(),
            uplift_rule: "same as t1b uplift_rule (t0_subgrid_closure_rate_v1, gain x K x R_sub / Lambda, land-weighted bilinear) at the node".into(),
            base_rule: "the T1a trunk level at the node, floored at sea level + floor on land".into(),
            roughness_rule: "keyed_uniform_white_v1 at 512 m: vz.keyed-rng/v1 root(seed) -> layer(guides recipe digest) -> chart 0 -> tile(8, ix, iy) of the node -> draw(t1b.guides.initial.roughness, ly*128+lx, 0); (2u-1) sqrt3 rms, rms = min(fraction x R_sub, cap)".into(),
            routing_rule: "priority_flood_steepest_v1 (t1b_lem::route) on the relaxed surface quantized to mm".into(),
            confine_mode: 1,
            floor_mode: 1,
            lithology_mode: 0,
            lithology: None,
            guide_rule: "guide_area_threshold_v1: a free node is a guide iff its drainage area (512 m cells, D8, closed chart edge) >= the threshold; its bed is the priority-flood fill of the relaxed quantized surface (monotone along receivers); the guide link is the D8 step to its receiver".into(),
        }
    }
}

impl GuideParamsV1 {
    pub fn steps(&self) -> usize { (self.sim_time_kyr / self.dt_kyr) as usize }
    /// Guides 0.4.0 (lane T3, TER-REGION v3): fewer, bigger valleys (threshold 16 km2) and the same
    /// exposed-lithology / threshold-slope law as T1b 0.6.0.
    pub fn v3() -> Self {
        GuideParamsV1 { guide_area_threshold_m2: 16_000_000, lithology_mode: 1, lithology: Some(crate::t1b_litho::LithologyV1::standin_v1()), ..GuideParamsV1::default() }
    }
    pub fn generator_version(&self) -> &'static str { if self.lithology_mode == 1 { GUIDES_GENERATOR_VERSION_LITHO } else { GUIDES_GENERATOR_VERSION } }
}

#[derive(Clone, Debug, Serialize)]
pub struct GuideRecipeV1 {
    pub schema: String,
    pub layer: String,
    pub generator: GeneratorRef,
    pub params: GuideParamsV1,
    pub laws: ReliefLawsV1,
    pub universe_digest: String,
    pub upstream: Vec<(String, String)>,
    pub seed: u64,
    pub chart: ChartDeclV1,
    pub simtime_label: String,
    pub closure_label: String,
    pub label: String,
}

pub fn guide_recipe(t0: &T0, extras: &T0Extras, site: &SiteRecord, skel: &SkeletonView, params: GuideParamsV1) -> R<GuideRecipeV1> {
    if params.dt_kyr <= 0 || params.sim_time_kyr % params.dt_kyr != 0 { return Err("guides: sim time must be whole steps".into()); }
    match (params.lithology_mode, &params.lithology) {
        (0, None) => {}
        (1, Some(l)) => l.validate()?,
        _ => return Err("guides: lithology_mode 1 needs a lithology table (and 0 none)".into()),
    }
    let mut upstream = vec![
        ("site".to_string(), site.bytes_digest.clone()),
        ("t0.hydrology".to_string(), t0.hydrology_digest.clone()),
        ("t0.subgrid_relief_sidecar".to_string(), t0.subgrid_digest.clone()),
        ("t0.terrain".to_string(), t0.terrain_digest.clone()),
        ("t1a.skeleton".to_string(), skel.digest.clone()),
    ];
    upstream.sort();
    Ok(GuideRecipeV1 {
        schema: GUIDES_RECIPE_SCHEMA.into(),
        layer: GUIDES_LAYER.into(),
        generator: GeneratorRef { name: GUIDES_GENERATOR_NAME.into(), version: params.generator_version().into(), tier: "Preview".into() },
        params,
        laws: extras.laws.clone(),
        universe_digest: site.universe_digest.clone(),
        upstream,
        seed: site.seed,
        chart: ChartDeclV1::new(site.chart_origin_lat_e7, site.chart_origin_lon_e7, (t0.radius_m * 1000.0).round() as i64),
        simtime_label: SIMTIME_LABEL.into(),
        closure_label: CLOSURE_LABEL.into(),
        label: GUIDES_LABEL.into(),
    })
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct GuideStats {
    pub fixed_nodes: u32,
    pub free_nodes: u32,
    pub guide_nodes: u32,
    pub guide_threshold_km2: f64,
    pub guide_max_area_km2: f64,
    /// Guide nodes whose discharge (area x T0 runoff) reaches the river threshold: rivers T0 does not have (F2 at 512 m).
    pub guide_nodes_ge_river_threshold: u32,
    pub guide_max_discharge_m3_s: f64,
    /// F3 at 512 m: per land T0 cell, the share of its free land nodes whose flow path ends on a fixed node
    /// of the same T0 cell or of the cell's T0 receiver; the worst cell (cells with >= 100 free nodes).
    pub f3_min_share_permille: f64,
    pub f3_min_share_cell: u32,
    pub f3_cells_below_950: u32,
    pub f3_cells_checked: u32,
    pub pits_after: u32,
    /// Nodes the confined flood could not reach (flooded unconstrained), final surface.
    pub unconfined_nodes: u32,
    pub mean_free_m: f64,
}

pub struct GuideLayer {
    pub recipe: GuideRecipeV1,
    pub recipe_digest: String,
    /// Relaxed surface, mm (fixed nodes carry their fixed value).
    pub height_mm: Vec<i32>,
    /// Priority-flood fill of `height_mm`, mm: the guide bed.
    pub bed_mm: Vec<i32>,
    /// 0 free, 1 T1a channel / mouth, 2 marine, 4 lake, 5 guide (free node at or above the threshold).
    pub fixed: Vec<u8>,
    pub receiver: Vec<u32>,
    pub area_m2: Vec<f32>,
    /// The fixed node each node's flow path ends on.
    pub terminal: Vec<u32>,
    pub stats: GuideStats,
    pub ledger: Ledger,
    pub series: Vec<SeriesPoint>,
    pub digest: String,
    pub header: Value,
    pub blobs: Vec<(String, Vec<u8>)>,
}

fn mm(v: f64) -> R<i32> {
    i32::try_from(vz_tilekey::chart::quantize_mm(v * 1000.0)).map_err(|_| format!("height {v} m does not fit i32 mm"))
}

/// Build the chart-wide guide layer. Pure function of its inputs.
pub fn build_guides(t0: &T0, extras: &T0Extras, site: &SiteRecord, skel: &SkeletonView, params: GuideParamsV1) -> R<GuideLayer> {
    let recipe = guide_recipe(t0, extras, site, skel, params)?;
    let recipe_digest = pf_contracts::artifact::digest_canonical(&recipe).map_err(|e| e.to_string())?;
    let p = &recipe.params;
    let chart = Chart::new(recipe.chart.clone()).map_err(|e| e.to_string())?;
    let sea = t0.sea_level_m;
    let floor = sea + p.base_floor_above_sea_mm as f64 / 1000.0;
    let gain = p.uplift_gain_permille as f64 / 1000.0;
    let rate = |c: u32| closure_rate_m_per_my(t0, &extras.laws, c);
    let layer = Scope::root(site.seed).layer(&recipe_digest).map_err(|e| e.to_string())?.chart(0).map_err(|e| e.to_string())?;
    let mut tiles: BTreeMap<(i64, i64), Scope> = BTreeMap::new();
    let mut fixed = vec![FREE; NODES];
    let mut z0 = vec![0f64; NODES];
    let mut up = vec![0f64; NODES];
    for n in 0..NODES {
        let (gx, gy) = node_g(n);
        let c = skel.t0_cell[n];
        match t0.water(c) {
            CellWater::Marine => { fixed[n] = MARINE; z0[n] = sea; continue; }
            CellWater::Lake => { fixed[n] = LAKE; z0[n] = t0.level(c); continue; }
            CellWater::Land => {}
        }
        if (skel.kind[n] == 1 && skel.channel[n] == 1) || skel.kind[n] == 3 {
            fixed[n] = CHANNEL; z0[n] = skel.trunk_level_m[n] as f64; continue;
        }
        let (la, lo) = chart.inverse_mm(gx * SPACING_MM, gy * SPACING_MM).map_err(|e| e.to_string())?;
        up[n] = gain * bilinear(t0, &rate, la * RAD2DEG, lo * RAD2DEG, true);
        let rms = ((t0.subgrid_relief_m[c as usize] as f64).max(0.0) * p.roughness_rms_fraction_ppm as f64 / 1e6)
            .min(p.roughness_rms_cap_mm as f64 / 1000.0);
        let (tx, ty) = (gx.div_euclid(128), gy.div_euclid(128));
        let (lx, ly) = (gx - tx * 128, gy - ty * 128);
        if !tiles.contains_key(&(tx, ty)) {
            tiles.insert((tx, ty), layer.tile(LOD, tx, ty).map_err(|e| e.to_string())?);
        }
        let u = tiles[&(tx, ty)].draw(PURPOSE_GUIDE_ROUGHNESS, (ly * 128 + lx) as u64, 0).map_err(|e| e.to_string())?.uniform_f64();
        z0[n] = (skel.trunk_level_m[n] as f64).max(floor) + (2.0 * u - 1.0) * 3f64.sqrt() * rms;
    }
    let cfg = LemConfig {
        k_per_my: extras.laws.stream_power_k_per_my,
        m: extras.laws.stream_power_m,
        d_m2_per_my: extras.laws.hillslope_diffusivity_m2_per_my,
        stability_limit: extras.laws.diffusion_stability_limit,
        dt_my: p.dt_kyr as f64 / 1000.0,
        steps: p.steps(),
        stats_every_steps: ((p.stats_every_kyr / p.dt_kyr) as usize).max(1),
        dx_m: SPACING_MM as f64 / 1000.0,
    };
    let core = vec![true; NODES];
    let inflow = vec![0f64; NODES];
    let zone: Vec<u32> = skel.t0_cell.clone();
    let drain: Vec<u32> = zone.iter().map(|&c| t0.flow(c).unwrap_or(NONE)).collect();
    let zones = (p.confine_mode == 1).then_some(t1b_lem::Zones { zone: &zone, drain: &drain });
    let floor_v: Vec<f64> = (0..NODES).map(|n| (skel.trunk_level_m[n] as f64).max(floor)).collect();
    let litho = match &p.lithology {
        Some(table) if p.lithology_mode == 1 => {
            let pos = |n: usize| -> R<(f64, f64)> {
                let (gx, gy) = node_g(n);
                let (la, lo) = chart.inverse_mm(gx * SPACING_MM, gy * SPACING_MM).map_err(|e| e.to_string())?;
                Ok((la * RAD2DEG, lo * RAD2DEG))
            };
            Some(crate::t1b_litho::build_standin(table, t0, &extras.rock(), NODES, &pos)?)
        }
        _ => None,
    };
    let res = t1b_lem::run_litho(N, &fixed, &z0, &up, &inflow, &core, &cfg, zones, (p.floor_mode == 1).then_some(&floor_v[..]),
        litho.as_ref().map(|l| l as &dyn crate::t1b_litho::ExposedLithology));
    let height_mm: Vec<i32> = (0..NODES).map(|n| mm(if fixed[n] == FREE { res.z[n] } else { z0[n] })).collect::<R<_>>()?;
    let zq: Vec<f64> = height_mm.iter().map(|v| *v as f64 / 1000.0).collect();
    let r = t1b_lem::route_zoned(N, &zq, &fixed, cfg.dx_m, zones);
    let area = t1b_lem::accumulate(&r, cfg.dx_m * cfg.dx_m, None);
    let term = t1b_lem::terminals(&r);
    let thr = p.guide_area_threshold_m2 as f64;
    let mut st = GuideStats { guide_threshold_km2: thr / 1e6, unconfined_nodes: r.unconfined, ..Default::default() };
    let mut pub_fixed = fixed.clone();
    let q_of = |a: f64| a * t0.runoff_m_per_yr / SECONDS_PER_YEAR;
    let mut mean = (0.0, 0usize);
    for n in 0..NODES {
        if fixed[n] != FREE { st.fixed_nodes += 1; continue; }
        st.free_nodes += 1;
        mean.0 += zq[n]; mean.1 += 1;
        if zq[n] < r.fill[n] { st.pits_after += 1; }
        if area[n] >= thr {
            pub_fixed[n] = GUIDE_NODE;
            st.guide_nodes += 1;
            st.guide_max_area_km2 = st.guide_max_area_km2.max(area[n] / 1e6);
            let q = q_of(area[n]);
            st.guide_max_discharge_m3_s = st.guide_max_discharge_m3_s.max(q);
            if q >= 10.0 { st.guide_nodes_ge_river_threshold += 1; }
        }
    }
    st.mean_free_m = mean.0 / mean.1.max(1) as f64;
    // F3 at 512 m.
    let mut per: BTreeMap<u32, (u64, u64)> = BTreeMap::new();
    for n in 0..NODES {
        if fixed[n] != FREE { continue; }
        let c = skel.t0_cell[n];
        let tc = skel.t0_cell[term[n] as usize];
        let ok = tc == c || t0.flow(c) == Some(tc);
        let e = per.entry(c).or_default();
        e.0 += 1;
        if ok { e.1 += 1; }
    }
    st.f3_min_share_permille = 1000.0;
    for (c, (tot, ok)) in &per {
        if *tot < 100 { continue; }
        st.f3_cells_checked += 1;
        let s = 1000.0 * *ok as f64 / *tot as f64;
        if s < 950.0 { st.f3_cells_below_950 += 1; }
        if s < st.f3_min_share_permille { st.f3_min_share_permille = s; st.f3_min_share_cell = *c; }
    }
    let bed_mm: Vec<i32> = r.fill.iter().map(|v| mm(*v)).collect::<R<_>>()?;
    let area_f: Vec<f32> = area.iter().map(|v| *v as f32).collect();
    let le_i32 = |v: &[i32]| -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() };
    let le_u32 = |v: &[u32]| -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() };
    let le_f32 = |v: &[f32]| -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() };
    let fields: Vec<(&str, &str, &str, Vec<u8>)> = vec![
        ("height_mm", "i32", "relaxed 512 m surface (t1b-relaxation-not-history), mm; fixed nodes carry their fixed value", le_i32(&height_mm)),
        ("bed_mm", "i32", "priority-flood fill of height_mm, mm: the guide bed (monotone along receivers)", le_i32(&bed_mm)),
        ("fixed", "u8", "0 free, 1 T1a channel or mouth (fixed), 2 marine, 4 lake, 5 guide (free node with area >= threshold)", pub_fixed.clone()),
        ("receiver", "u32", "D8 receiver on the quantized relaxed surface (NONE on fixed nodes)", le_u32(&r.receiver)),
        ("area_m2", "f32", "drainage area, 512 m cells, closed chart edge", le_f32(&area_f)),
        ("terminal", "u32", "fixed node the flow path ends on", le_u32(&term)),
    ];
    let mut entries = Vec::new();
    let mut blobs = Vec::new();
    for (name, dtype, note, bytes) in fields {
        let file = format!("fields/{name}.{dtype}");
        entries.push(serde_json::json!({"name": name, "dtype": dtype, "file": file, "note": note,
            "digest": format!("blake3:{}", blake3::hash(&bytes).to_hex())}));
        blobs.push((file, bytes));
    }
    let header = serde_json::json!({
        "schema": GUIDES_SCHEMA, "layer": GUIDES_LAYER, "label": GUIDES_LABEL, "simtime_label": SIMTIME_LABEL,
        "lattice": {"lod": LOD, "spacing_mm": SPACING_MM, "g_min": G_MIN, "g_max": G_MAX, "order": "node = (gy - g_min) * 513 + (gx - g_min) (the T1a order)"},
        "recipe": recipe, "recipe_digest": recipe_digest, "stats": st, "ledger": res.ledger, "series": res.series, "fields": entries,
        "invented": ["the tributary layout between T1a channels, by the relaxation process only (keyed initial roughness at 512 m)"],
        "approximations": ["512 m D8 guide links are straight; T1b rasterizes them at 32 m", "same law and constants as T1b; gain reused from T1b 0.4.0, not re-fitted at 512 m", "SimTime is t1b-relaxation-not-history"],
    });
    let header = serde_json::from_str::<Value>(&header.to_string()).map_err(|e| e.to_string())?;
    let digest = pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())?;
    Ok(GuideLayer { recipe, recipe_digest, height_mm, bed_mm, fixed: pub_fixed, receiver: r.receiver, area_m2: area_f, terminal: term,
        stats: st, ledger: res.ledger, series: res.series, digest, header, blobs })
}

pub fn publish_guides(g: &GuideLayer, out: &Path) -> R<()> {
    std::fs::create_dir_all(out.join("fields")).map_err(|e| e.to_string())?;
    for (file, bytes) in &g.blobs { std::fs::write(out.join(file), bytes).map_err(|e| format!("{file}: {e}"))?; }
    let doc = serde_json::json!({"header": g.header, "digest": g.digest});
    std::fs::write(out.join("guides.json"), serde_json::to_vec_pretty(&doc).map_err(|e| e.to_string())?).map_err(|e| e.to_string())
}

/// What T1b reads from the guide layer (digest-verified when loaded from disk).
pub struct GuideView {
    pub digest: String,
    pub threshold_m2: i64,
    pub bed_mm: Vec<i32>,
    pub fixed: Vec<u8>,
    pub receiver: Vec<u32>,
    pub terminal: Vec<u32>,
    /// Relaxed 512 m surface, mm (the T1b 0.7.0 ring condition; lane T3b). Read-only view field, not digested here.
    pub height_mm: Vec<i32>,
}

impl GuideView {
    pub fn from_layer(g: &GuideLayer) -> GuideView {
        GuideView { digest: g.digest.clone(), threshold_m2: g.recipe.params.guide_area_threshold_m2, bed_mm: g.bed_mm.clone(),
            fixed: g.fixed.clone(), receiver: g.receiver.clone(), terminal: g.terminal.clone(), height_mm: g.height_mm.clone() }
    }
    pub fn load(dir: &Path) -> R<GuideView> {
        let bytes = std::fs::read(dir.join("guides.json")).map_err(|e| format!("guides.json in {}: {e}", dir.display()))?;
        let v: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
        let h = v.get("header").ok_or("guides: no header")?;
        let stated = v.get("digest").and_then(Value::as_str).ok_or("guides: no digest")?;
        let got = pf_contracts::artifact::digest_canonical(h).map_err(|e| e.to_string())?;
        if got != stated { return Err(format!("guides header digest {got} differs from stated {stated}; refusing")); }
        let mut blobs = BTreeMap::new();
        for f in h.get("fields").and_then(Value::as_array).ok_or("guides: fields")? {
            let name = f.get("name").and_then(Value::as_str).ok_or("name")?;
            let file = f.get("file").and_then(Value::as_str).ok_or("file")?;
            let b = std::fs::read(dir.join(file)).map_err(|e| format!("{file}: {e}"))?;
            let d = format!("blake3:{}", blake3::hash(&b).to_hex());
            if Some(d.as_str()) != f.get("digest").and_then(Value::as_str) { return Err(format!("guides field {name} digest mismatch; refusing")); }
            blobs.insert(name.to_string(), b);
        }
        let b = |n: &str| blobs.get(n).cloned().ok_or_else(|| format!("guides field {n} missing"));
        let i32s = |b: Vec<u8>| b.chunks_exact(4).map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect::<Vec<_>>();
        let u32s = |b: Vec<u8>| b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect::<Vec<_>>();
        let gv = GuideView {
            digest: stated.to_string(),
            threshold_m2: h.pointer("/recipe/params/guide_area_threshold_m2").and_then(Value::as_i64).ok_or("threshold")?,
            bed_mm: i32s(b("bed_mm")?), fixed: b("fixed")?, receiver: u32s(b("receiver")?), terminal: u32s(b("terminal")?),
            height_mm: i32s(b("height_mm")?),
        };
        if gv.bed_mm.len() != NODES || gv.fixed.len() != NODES || gv.height_mm.len() != NODES { return Err("guides: field sizes".into()); }
        Ok(gv)
    }
    #[inline]
    pub fn is_guide(&self, n: usize) -> bool { self.fixed[n] == GUIDE_NODE }
    #[inline]
    pub fn bed_m(&self, n: usize) -> f64 { self.bed_mm[n] as f64 / 1000.0 }
    /// The relaxed 512 m surface by Catmull-Rom bicubic interpolation between lod-8 nodes at lod-4 post (gx, gy)
    /// (C1: no slope breaks along the 512 m lines; nodes off the chart: nearest on it). T1b 0.7.0 initial base and ring.
    pub fn surface_bicubic_m(&self, gx: i64, gy: i64) -> f64 {
        let l8 = crate::t1b::L8_TO_L4;
        let (x8, y8) = (gx.div_euclid(l8), gy.div_euclid(l8));
        let (tx, ty) = (gx.rem_euclid(l8) as f64 / l8 as f64, gy.rem_euclid(l8) as f64 / l8 as f64);
        let w = |t: f64| -> [f64; 4] {
            let (t2, t3) = (t * t, t * t * t);
            [0.5 * (-t3 + 2.0 * t2 - t), 0.5 * (3.0 * t3 - 5.0 * t2 + 2.0), 0.5 * (-3.0 * t3 + 4.0 * t2 + t), 0.5 * (t3 - t2)]
        };
        let (wx, wy) = (w(tx), w(ty));
        let v = |x: i64, y: i64| self.height_mm[node_index(x.clamp(G_MIN, G_MAX), y.clamp(G_MIN, G_MAX)).unwrap()] as f64 / 1000.0;
        let mut acc = 0.0;
        for (j, wyj) in wy.iter().enumerate() {
            let mut row = 0.0;
            for (i, wxi) in wx.iter().enumerate() { row += wxi * v(x8 - 1 + i as i64, y8 - 1 + j as i64); }
            acc += wyj * row;
        }
        acc
    }
    /// The relaxed 512 m surface bilinear between lod-8 nodes at lod-4 post (gx, gy) (nodes off the chart: nearest on it).
    pub fn surface_bilinear_m(&self, gx: i64, gy: i64) -> f64 {
        let l8 = crate::t1b::L8_TO_L4;
        let (x8, y8) = (gx.div_euclid(l8), gy.div_euclid(l8));
        let (tx, ty) = (gx.rem_euclid(l8) as f64 / l8 as f64, gy.rem_euclid(l8) as f64 / l8 as f64);
        let v = |x: i64, y: i64| self.height_mm[node_index(x.clamp(G_MIN, G_MAX), y.clamp(G_MIN, G_MAX)).unwrap()] as f64 / 1000.0;
        let a = v(x8, y8) * (1.0 - tx) + v(x8 + 1, y8) * tx;
        let b = v(x8, y8 + 1) * (1.0 - tx) + v(x8 + 1, y8 + 1) * tx;
        a * (1.0 - ty) + b * ty
    }
}
