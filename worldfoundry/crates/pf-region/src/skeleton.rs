//! Building the T1a skeleton: segments, keyed planform, placement with
//! conflict retries, bed profile, hillslope routing, discharge.

use crate::t0::{CellWater, SiteRecord, T0};
use crate::*;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, VecDeque};
use vz_keyed_rng::Scope;
use vz_tilekey::chart::{Chart, ChartDeclV1, DEG2RAD, RAD2DEG};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, serde::Serialize)]
pub enum SegKind {
    /// Mouth-less end of an endorheic T0 floor: one node, terminal.
    Terminal,
    /// One T0 flow edge, cell centre to receiver centre.
    T0Flow,
    /// Headwater extension of a T0 headwater cell (declared invention).
    HeadExtension,
}

#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub enum SegEnd {
    Merge { node: u32, into_segment: u32 },
    Mouth { node: u32 },
    OffChart,
    Terminal,
    Dangling,
    Outside,
}

#[derive(Clone, Debug, serde::Serialize)]
pub struct Segment {
    pub id: u32,
    pub kind: SegKind,
    pub address: String,
    /// Start cell b(s) and entry cell e(s) (where the water goes).
    pub from_cell: u32,
    pub to_cell: u32,
    pub start_mm: [i64; 2],
    pub end_mm: [i64; 2],
    pub level_start_m: f64,
    pub level_end_m: f64,
    pub q_t0_m3_s: f64,
    pub is_t0_river: bool,
    pub t0_reach: u32,
    pub depth: u32,
    pub amplitude_mm: i64,
    pub attempts: u32,
    pub straight_fallback: bool,
    pub forced: bool,
    pub conflict_notes: Vec<String>,
    pub end: SegEnd,
    pub enters_chart: bool,
    pub enters_computation_lattice: bool,
    pub inflow_m3_s: f64,
    #[serde(skip)]
    pub path: Vec<u32>,
    /// Path on the computation lattice (chart + halo), before cropping.
    #[serde(skip)]
    pub path_computation: Vec<u32>,
    pub path_nodes: u32,
    pub path_nodes_computation: u32,
    /// Keyed control offsets actually used, permille of the amplitude.
    pub control_offsets_permille: Vec<i64>,
}

/// A built skeleton: node fields plus tables.
pub struct Skeleton {
    pub recipe: SkeletonRecipeV1,
    pub recipe_digest: String,
    pub chart: Chart,
    pub params: SkeletonParamsV1,
    pub t0_cell: Vec<u32>,
    pub area_m2: Vec<f64>,
    /// 0 hillslope, 1 axis (corridor path), 2 marine, 3 mouth, 4 lake hillslope.
    pub kind: Vec<u8>,
    pub owner: Vec<u32>,
    pub receiver: Vec<u32>,
    pub join_node: Vec<u32>,
    pub join_segment: Vec<u32>,
    /// Hillslope nodes of a cell whose own channel lies outside the chart:
    /// they drain off-chart across the chart edge inside their own cell.
    pub offchart_fallback: Vec<bool>,
    pub bed_m: Vec<f64>,
    pub trunk_level_m: Vec<f64>,
    pub flow_dist_m: Vec<f64>,
    pub discharge_m3_s: Vec<f64>,
    pub basin: Vec<u32>,
    pub t0_reach: Vec<u32>,
    pub segments: Vec<Segment>,
    pub cells: BTreeMap<u32, CellInfo>,
    pub stats: BuildStats,
    pub site: SiteRecord,
    pub sea_level_m: f64,
    pub threshold_m3_s: f64,
    /// Lane CLIMATE-RUNOFF (key on only): the climate runoff R (m/yr) at each chart node, the digest of the same skeleton
    /// built without the key (the draw twin T1b keys its draws by, so climate changes physics, not draws), and the
    /// climate T0 lake budgets.
    pub climate_runoff_m_per_yr: Option<Vec<f64>>,
    pub climate_twin_digest: Option<String>,
    pub physical_climate_twin_digest: Option<String>,
    pub climate_lakes: Vec<crate::climate_runoff::LakeBudget>,
}

#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct BuildStats {
    pub segments_t0_flow: u32,
    pub segments_head: u32,
    pub segments_terminal: u32,
    pub segments_retried: u32,
    pub segments_straight_fallback: u32,
    pub segments_forced: u32,
    pub bed_nodes_carved: u32,
    pub bed_max_carve_m: f64,
    pub hillslope_pass1: u32,
    pub hillslope_offchart_fallback: u32,
    pub hillslope_pass2_foreign: u32,
    pub hillslope_pass3_coast: u32,
    pub hillslope_unrouted: u32,
    pub hillslope_offchart_fallback_in_chart: u32,
    pub halo_nodes: i64,
    pub computation_posts_per_axis: u32,
    pub site_in_chart_mm_recomputed: [i64; 2],
}

#[derive(Clone, Debug, serde::Serialize)]
pub struct CellInfo {
    pub cell: u32,
    pub ij: [usize; 2],
    pub centre_mm: [i64; 2],
    pub water: String,
    pub elevation_m: f64,
    pub subgrid_relief_m: f64,
    pub level_m: f64,
    pub flow_to: Option<u32>,
    pub t0_reach: Option<u32>,
    pub basin: Option<u32>,
    pub final_basin: Option<u32>,
    pub q_t0_m3_s: f64,
    pub river_inflow: u32,
    pub nodes_in_chart: u32,
    pub area_in_chart_m2: f64,
    pub cell_area_m2: f64,
    pub inside_permille: i64,
    /// Lane T0-IMPRINT-PROD (params.imprint only): the valley-axis junction in chart mm (centre_mm stays the T0 cell centre).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub junction_mm: Option<[i64; 2]>,
}

fn opt(v: u32) -> Option<u32> { (v != NONE).then_some(v) }

/// Lane T0-IMPRINT-PROD keyed-draw purpose of the junction offsets (declared in the recipe seed scope only when on).
pub const PURPOSE_JUNCTION: &str = "t1a.junction.offset";

/// Lane T0-IMPRINT-PROD (ADR 045, design 2026-10-21): where a T0 cell's valley axes meet (the junction). Opt-in; absent from
/// the params JSON when off, so every byte of every tier is unchanged. Integers only (exact canonical JSON).
#[derive(Clone, Debug, serde::Serialize, PartialEq, Eq)]
pub struct ImprintParamsV1 {
    /// 0 = the T0 cell centre (today's rule, kept so the F3 reading can be switched on alone);
    /// 1 = keyed offset: centre + (u_lat, u_lon) x offset x half cell, u keyed uniforms in [-1, 1];
    /// 2 = low ground: the lowest point of the Catmull-Rom bicubic (C1) T0 elevation surface over the central region
    ///     |d| <= offset x half cell (17 x 17 samples, ties to the first), plus a keyed jitter of +-jitter x half cell,
    ///     clamped to the same region.
    pub junction_mode: u32,
    /// Half size of the junction region per axis (lat and lon), permille of the half cell (< 1000 keeps the junction in its cell).
    pub junction_offset_permille: i64,
    /// Mode 2 only: keyed jitter per axis, permille of the half cell.
    pub junction_jitter_permille: i64,
    pub junction_rule: String,
    /// F3 reading of the T1a checks: "per_t0_outlet_v1" (ADR 045; the per-cell share stays as a diagnostic).
    pub f3_rule: String,
    pub head_rule: String,
}

impl ImprintParamsV1 {
    pub fn v1(mode: u32) -> Self {
        ImprintParamsV1 {
            junction_mode: mode,
            junction_offset_permille: 600,
            junction_jitter_permille: if mode == 2 { 150 } else { 0 },
            junction_rule: match mode {
                0 => "cell_centre_v1: the junction is the T0 cell centre".into(),
                1 => "keyed_offset_v1: junction = cell centre + (u_lat, u_lon) x junction_offset_permille/1000 x half cell in (lat, lon) degrees; u = 2 w - 1, w = vz.keyed-rng/v1 draw (purpose t1a.junction.offset, address junc/<cell>, counters 0 lat and 1 lon, attempt 0) under the skeleton draw scope".into(),
                _ => "low_ground_bicubic_v1: junction = argmin of the Catmull-Rom bicubic of the T0 cell-centre elevations (lat, lon index space, latitude clamped at the poles, longitude wrapped) over 17 x 17 samples of the region |d| <= junction_offset_permille/1000 x half cell, first minimum in row order; plus keyed jitter (purpose t1a.junction.offset, address junc/<cell>, counters 0 and 1) x junction_jitter_permille/1000 x half cell, clamped to the region".into(),
            },
            f3_rule: "per_t0_outlet_v1: ADR 045 F3 (per T0 final outlet with an area tolerance; per-cell share diagnostic)".into(),
            head_rule: "head_clip_own_cell_v1: the head extension (head_extension_rule, from the junction) is shortened to the largest of 64 equal fractions whose points, and 8 points at (planform amplitude + 512 m) around each, lie in the cell; a head shorter than 1024 m is dropped".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.junction_mode > 2 { return Err("imprint: junction_mode is 0, 1 or 2".into()); }
        if !(0..1000).contains(&self.junction_offset_permille) || !(0..=1000).contains(&self.junction_jitter_permille) {
            return Err("imprint: junction_offset_permille in [0, 1000), junction_jitter_permille in [0, 1000]".into());
        }
        Ok(())
    }
}

/// Catmull-Rom weights at fraction t.
fn cr_w(t: f64) -> [f64; 4] {
    let (t2, t3) = (t * t, t * t * t);
    [-0.5 * t3 + t2 - 0.5 * t, 1.5 * t3 - 2.5 * t2 + 1.0, -1.5 * t3 + 2.0 * t2 + 0.5 * t, 0.5 * t3 - 0.5 * t2]
}

/// Catmull-Rom bicubic of the T0 cell-centre elevations at fractional cell index (fi, fj) (cell (i, j) centre = (i, j)).
fn t0_bicubic(t0: &T0, fi: f64, fj: f64) -> f64 {
    let (i0, j0) = (fi.floor() as i64, fj.floor() as i64);
    let (wi, wj) = (cr_w(fi - i0 as f64), cr_w(fj - j0 as f64));
    let mut v = 0.0;
    for (a, wa) in wi.iter().enumerate() {
        let ii = (i0 - 1 + a as i64).clamp(0, t0.nlat as i64 - 1) as usize;
        for (b, wb) in wj.iter().enumerate() {
            let jj = (j0 - 1 + b as i64).rem_euclid(t0.nlon as i64) as usize;
            v += wa * wb * t0.elevation_m[ii * t0.nlon + jj] as f64;
        }
    }
    v
}

/// Lane T0-IMPRINT-PROD head_clip_own_cell_v1: the largest fraction k in [0, 1] (64 steps) such that every point of the
/// head extension from the junction `pc` toward `h` (chart mm) up to k, and the 8 points at distance `clear_mm` around it,
/// lies in T0 cell `c`.
fn head_clip_fraction(decl: &ChartDeclV1, chart: &Chart, t0: &T0, c: u32, pc: [f64; 2], h: [f64; 2], clear_mm: f64) -> f64 {
    let cell_at = |x: f64, y: f64| -> u32 {
        let (la, lo) = match chart.inverse_mm(x.round() as i64, y.round() as i64) { Ok(v) => v, Err(_) => inverse_any(decl, x, y) };
        t0.cell_of_deg(la * RAD2DEG, lo * RAD2DEG)
    };
    const STEPS: usize = 64;
    let mut ok = 0usize;
    for q in 1..=STEPS {
        let t = q as f64 / STEPS as f64;
        let (x, y) = (pc[0] + t * (h[0] - pc[0]), pc[1] + t * (h[1] - pc[1]));
        let inside = cell_at(x, y) == c && (0..8).all(|d| {
            let a = d as f64 * std::f64::consts::FRAC_PI_4;
            cell_at(x + clear_mm * a.cos(), y + clear_mm * a.sin()) == c
        });
        if !inside { break; }
        ok = q;
    }
    ok as f64 / STEPS as f64
}

/// The junction of T0 cell `c` in (lat, lon) degrees under `ip` (lane T0-IMPRINT-PROD).
pub fn junction_deg(t0: &T0, c: u32, ip: &ImprintParamsV1, root: &Scope) -> R<(f64, f64)> {
    let (la, lo) = t0.centre_deg(c);
    if ip.junction_mode == 0 { return Ok((la, lo)); }
    let (hla, hlo) = (0.5 * 180.0 / t0.nlat as f64, 0.5 * 360.0 / t0.nlon as f64);
    let f = ip.junction_offset_permille as f64 / 1000.0;
    let scope = root.address(&format!("junc/{c}")).map_err(|e| e.to_string())?;
    let u = |k: u64| -> R<f64> { Ok(2.0 * scope.draw(PURPOSE_JUNCTION, k, 0).map_err(|e| e.to_string())?.uniform_f64() - 1.0) };
    if ip.junction_mode == 1 { return Ok((la + u(0)? * f * hla, lo + u(1)? * f * hlo)); }
    // Mode 2: lowest C1 ground over the central region (offsets in half-cell units), then keyed jitter.
    let (i, j) = t0.ij(c);
    const K: i64 = 8;
    let mut best = (f64::INFINITY, 0.0f64, 0.0f64);
    for a in -K..=K { for b in -K..=K {
        let (da, db) = (f * a as f64 / K as f64, f * b as f64 / K as f64);
        let z = t0_bicubic(t0, i as f64 + 0.5 * da, j as f64 + 0.5 * db);
        if z < best.0 { best = (z, da, db); }
    } }
    let jf = ip.junction_jitter_permille as f64 / 1000.0;
    let (da, db) = ((best.1 + u(0)? * jf).clamp(-f, f), (best.2 + u(1)? * jf).clamp(-f, f));
    Ok((la + da * hla, lo + db * hlo))
}

pub struct Inputs<'a> {
    pub t0: &'a T0,
    pub site: &'a SiteRecord,
    pub params: SkeletonParamsV1,
}

/// Content identities for the physical inputs read by the regional construction.
/// Source artifact and site-byte digests remain in `upstream` for verification;
/// these separate identities drive stochastic realization. Labels, paths, producer
/// receipts and per-layer provenance never enter them. Physical params/laws and
/// the existing climate-pairing rules remain in the draw recipe unchanged.
pub fn physical_upstream(t0: &T0, site: &SiteRecord) -> R<Vec<(String, String)>> {
    let digest = |v: serde_json::Value| pf_contracts::artifact::digest_canonical(&v).map_err(|e| e.to_string());
    let site_digest = digest(serde_json::json!({
        "schema": "vz.regional-physical-site/v1", "seed": site.seed, "level": site.level,
        "chart_origin_e7": [site.chart_origin_lat_e7, site.chart_origin_lon_e7],
        "site_point_latlon_deg": site.site_point_latlon_deg, "site_in_chart_mm": site.site_in_chart_mm,
        "walk_tile": site.walk_tile, "region": site.region, "t0_cell": site.t0_cell, "reach_id": site.reach_id
    }))?;
    // Hash IEEE bytes, not rounded JSON numbers. The framing names and lengths
    // distinguish arrays, including empty arrays, across future schema changes.
    fn fields(parts: &[(&str, Vec<u8>)]) -> String {
        let mut h = blake3::Hasher::new();
        h.update(b"vz.regional-physical-fields/v1");
        for (name, bytes) in parts {
            h.update(&(name.len() as u64).to_le_bytes()); h.update(name.as_bytes());
            h.update(&(bytes.len() as u64).to_le_bytes()); h.update(bytes);
        }
        format!("blake3:{}", h.finalize().to_hex())
    }
    let f32b = |a: &[f32]| a.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<_>>();
    let u32b = |a: &[u32]| a.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<_>>();
    let terrain = fields(&[
        ("shape", [t0.nlat as u64, t0.nlon as u64].iter().flat_map(|v| v.to_le_bytes()).collect()),
        ("radius_sea_runoff", [t0.radius_m, t0.sea_level_m, t0.runoff_m_per_yr].iter().flat_map(|v| v.to_le_bytes()).collect()),
        ("elevation_m", f32b(&t0.elevation_m)),
        ("all_verified_terrain_payloads", serde_json::to_vec(&t0.terrain_field_digests).map_err(|e| e.to_string())?),
    ]);
    let hydrology = fields(&[
        ("flow_to", u32b(&t0.flow_to)), ("reach_id", u32b(&t0.reach_id)),
        ("basin_id", u32b(&t0.basin_id)), ("water_body_id", u32b(&t0.water_body_id)),
        ("discharge_m3_s", f32b(&t0.discharge_m3_s)), ("fill_level_m", f32b(&t0.fill_level_m)),
        ("bodies", serde_json::to_vec(&t0.bodies).map_err(|e| e.to_string())?),
        ("river_inflow", u32b(&t0.river_inflow)),
        ("all_verified_hydrology_payloads", serde_json::to_vec(&t0.hydrology_field_digests).map_err(|e| e.to_string())?),
        ("water_bodies", serde_json::to_vec(&t0.hydrology.water_bodies).map_err(|e| e.to_string())?),
        ("basins", serde_json::to_vec(&t0.hydrology.basins).map_err(|e| e.to_string())?),
        ("rivers", serde_json::to_vec(&t0.hydrology.rivers).map_err(|e| e.to_string())?),
    ]);
    let mut physical = vec![("physical.site".into(), site_digest), ("physical.t0.terrain".into(), terrain), ("physical.t0.hydrology".into(), hydrology)];
    physical.sort();
    if let Some(binding) = &site.draw_realization {
        if binding.schema != "vz.physical-draw-realization/v1" || binding.physical_upstream != physical {
            return Err("draw realization physical inputs differ; refusing legacy random-seed binding".into());
        }
        let expected: Vec<_> = physical.iter().map(|(r, _)| r).collect();
        let got: Vec<_> = binding.draw_upstream.iter().map(|(r, _)| r).collect();
        if got != expected || binding.draw_upstream.iter().any(|(_, d)| !d.starts_with("blake3:") || d.len() != 71 || !d[7..].bytes().all(|c| c.is_ascii_hexdigit())) {
            return Err("draw realization requires one blake3 seed per physical role".into());
        }
        // This is an explicit retained realization. No arbitrary old digest is
        // trusted as evidence: full current sources were verified by T0::load,
        // and the binding's physical precondition was recomputed above.
        return Ok(binding.draw_upstream.clone());
    }
    Ok(physical)
}

/// Replace only source identity roles with their physical counterparts. All other
/// input roles (including physical subgrid/strata payloads) retain existing keying.
pub fn project_physical_upstream(upstream: &mut Vec<(String, String)>) {
    let physical: Vec<_> = upstream.iter().filter_map(|(r, d)| r.strip_prefix("physical.").map(|r| (r.to_string(), d.clone()))).collect();
    // Even ordering remains legacy when no physical mode was declared.
    if physical.is_empty() { return; }
    upstream.retain(|(r, _)| !r.starts_with("physical.") && !physical.iter().any(|(p, _)| p == r));
    upstream.extend(physical);
    upstream.sort();
}

/// The same published physical skeleton, with source identities projected out of
/// its recipe. The normal header digest and all payload verification stay intact.
pub fn physical_header_digest(mut header: serde_json::Value) -> R<String> {
    let recipe = header.get_mut("recipe").ok_or("physical skeleton: no recipe")?;
    let mut upstream: Vec<(String, String)> = serde_json::from_value(recipe["upstream"].clone()).map_err(|e| e.to_string())?;
    project_physical_upstream(&mut upstream);
    recipe["upstream"] = serde_json::to_value(upstream).map_err(|e| e.to_string())?;
    let rd = pf_contracts::artifact::digest_canonical(recipe).map_err(|e| e.to_string())?;
    header["recipe_digest"] = rd.into();
    if let Some(twin) = header.as_object_mut().ok_or("physical skeleton: not an object")?.remove("physical_climate_twin_digest") {
        header["climate_twin_digest"] = twin;
    }
    // Climate is intentionally paired out by the downstream climate-free twin.
    // This identity still includes every field, curve, lattice, and model parameter.
    pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())
}

pub fn recipe_for(inp: &Inputs) -> R<SkeletonRecipeV1> {
    let params_digest = pf_contracts::artifact::digest_canonical(&inp.params).map_err(|e| e.to_string())?;
    let mut upstream = vec![
        ("site".to_string(), inp.site.bytes_digest.clone()),
        ("t0.hydrology".to_string(), inp.t0.hydrology_digest.clone()),
        ("t0.subgrid_relief_sidecar".to_string(), inp.t0.subgrid_digest.clone()),
        ("t0.terrain".to_string(), inp.t0.terrain_digest.clone()),
    ];
    if inp.site.physical_keying_enabled() {
        upstream.extend(physical_upstream(inp.t0, inp.site)?);
    }
    // Lane CLIMATE-RUNOFF: the climate artifact is an upstream only when the key is on.
    if let Some(cp) = &inp.params.climate_runoff {
        let c = inp.t0.climate.as_ref().ok_or("climate_runoff params set but no climate attached to T0")?;
        if &c.params != cp { return Err("climate_runoff params differ from the attached climate's".into()); }
        upstream.push(("t0.climate".to_string(), c.digest.clone()));
    }
    upstream.sort();
    Ok(SkeletonRecipeV1 {
        schema: RECIPE_SCHEMA.into(),
        layer: LAYER.into(),
        generator: GeneratorRef { name: GENERATOR_NAME.into(), version: GENERATOR_VERSION.into(), tier: "Preview".into() },
        params: inp.params.clone(),
        params_digest,
        universe_digest: inp.site.universe_digest.clone(),
        upstream,
        seed_scope: SeedScopeV1 {
            keyed_rng: vz_keyed_rng::KEYED_RNG_V1.into(),
            seed: inp.site.seed,
            ladder: vec!["layer = this recipe's digest".into(), "chart = 0".into(),
                "address = seg/<from_cell>-<to_cell> | head/<cell> (T0 L-grid element ids)".into()],
            purposes: if inp.params.imprint.as_ref().is_some_and(|i| i.junction_mode > 0) { vec![PURPOSE_PLANFORM.into(), PURPOSE_JUNCTION.into()] } else { vec![PURPOSE_PLANFORM.into()] },
        },
        chart: ChartDeclV1::new(inp.site.chart_origin_lat_e7, inp.site.chart_origin_lon_e7, (inp.t0.radius_m * 1000.0).round() as i64),
        lattice: lattice_v1(),
        key_domain: format!("{LAYER}/c0/{LOD}.*"),
        label: LABEL.into(),
    })
}

/// Legacy recipes keep their original keys. Explicit physical-key recipes project
/// their physical upstream roles before hashing. Climate pairing in either mode
/// removes climate params/upstreams so climate changes discharge, not draws.
pub fn draw_digest(recipe: &SkeletonRecipeV1) -> R<String> {
    let mut c = recipe.clone();
    project_physical_upstream(&mut c.upstream);
    if recipe.params.climate_runoff.is_none() { return c.digest(); }
    c.params.climate_runoff = None;
    c.params.runoff_rule = SkeletonParamsV1::default().runoff_rule;
    c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
    c.upstream.retain(|(r, _)| r != "t0.climate");
    c.digest()
}

/// Lane CLIMATE-RUNOFF: build with the key; with climate on, also build the climate-free twin (same T0, params without
/// climate) and record its digest (T1b keys its draws by it). Without the key this is `build`.
pub fn build_with_twin(inp: &Inputs) -> R<Skeleton> {
    let mut sk = build(inp)?;
    if inp.params.climate_runoff.is_some() {
        let mut p0 = inp.params.clone();
        p0.climate_runoff = None;
        p0.runoff_rule = SkeletonParamsV1::default().runoff_rule;
        let twin = build(&Inputs { t0: inp.t0, site: inp.site, params: p0 })?;
        sk.climate_twin_digest = Some(crate::publish::digest(&twin)?);
        if inp.site.physical_keying_enabled() {
            let (header, _) = crate::publish::header(&twin)?;
            sk.physical_climate_twin_digest = Some(physical_header_digest(serde_json::to_value(header).map_err(|e| e.to_string())?)?);
        }
    }
    Ok(sk)
}

/// Curve sample: chart mm and curve parameter t in [0, 1].
type Sample = ([f64; 2], f64);

/// Keyed planform curve of one segment. Interior control values are keyed
/// uniforms in [-1, 1]; the ends are 0, so junctions stay at the cell
/// junctions (T0 cell centres unless params.imprint). Interpolation is smoothstep between control values, so
/// |offset| <= amplitude everywhere (no overshoot).
pub fn curve(start: [f64; 2], end: [f64; 2], amplitude_mm: f64, p: &SkeletonParamsV1, scope: &Scope, attempt: u32)
    -> R<(Vec<Sample>, Vec<i64>)> {
    let d = [end[0] - start[0], end[1] - start[1]];
    let len = (d[0] * d[0] + d[1] * d[1]).sqrt();
    if len == 0.0 { return Ok((vec![(start, 0.0)], vec![])); }
    let n_int = ((len / p.planform_control_spacing_mm as f64).round() as usize).max(1);
    let mut u = vec![0.0f64; n_int + 1];
    if amplitude_mm > 0.0 && n_int >= 2 {
        for (k, uk) in u.iter_mut().enumerate().take(n_int).skip(1) {
            let w = scope.draw(PURPOSE_PLANFORM, k as u64, attempt).map_err(|e| e.to_string())?.uniform_f64();
            *uk = 2.0 * w - 1.0;
        }
    }
    let normal = [-d[1] / len, d[0] / len];
    let ns = ((len / p.planform_sample_step_mm as f64).ceil() as usize).max(1);
    let mut out = Vec::with_capacity(ns + 1);
    for q in 0..=ns {
        let t = q as f64 / ns as f64;
        let pos = t * n_int as f64;
        let k = (pos.floor() as usize).min(n_int - 1);
        let f = pos - k as f64;
        let sm = f * f * (3.0 - 2.0 * f);
        let v = u[k] + (u[k + 1] - u[k]) * sm;
        let off = amplitude_mm * v;
        out.push(([start[0] + t * d[0] + off * normal[0], start[1] + t * d[1] + off * normal[1]], t));
    }
    let ctrl = u.iter().map(|x| (x * 1000.0).round() as i64).collect();
    Ok((out, ctrl))
}

/// Rasterize samples to an 8-connected, loop-free sequence of lattice
/// coordinates (which may lie outside the lattice), each with its t.
pub fn rasterize(samples: &[Sample]) -> Vec<((i64, i64), f64)> {
    let mut out: Vec<((i64, i64), f64)> = Vec::new();
    let mut pos: BTreeMap<(i64, i64), usize> = BTreeMap::new();
    let push = |out: &mut Vec<((i64, i64), f64)>, pos: &mut BTreeMap<(i64, i64), usize>, g: (i64, i64), t: f64| {
        if let Some(&p) = pos.get(&g) {
            for (gg, _) in out.drain(p + 1..) { pos.remove(&gg); }
            return;
        }
        pos.insert(g, out.len());
        out.push((g, t));
    };
    for &(pt, t) in samples {
        let g = ((pt[0] / SPACING_MM as f64).round() as i64, (pt[1] / SPACING_MM as f64).round() as i64);
        match out.last().copied() {
            None => push(&mut out, &mut pos, g, t),
            Some((last, lt)) => {
                if last == g { continue; }
                let (dx, dy) = (g.0 - last.0, g.1 - last.1);
                let steps = dx.abs().max(dy.abs());
                // Bresenham-style 8-connected fill (DDA on the dominant axis).
                for s in 1..=steps {
                    let fx = last.0 as f64 + dx as f64 * s as f64 / steps as f64;
                    let fy = last.1 as f64 + dy as f64 * s as f64 / steps as f64;
                    let gg = (fx.round() as i64, fy.round() as i64);
                    let tt = lt + (t - lt) * s as f64 / steps as f64;
                    if out.last().map(|x| x.0) != Some(gg) { push(&mut out, &mut pos, gg, tt); }
                }
            }
        }
    }
    // Corner cutting: drop a node whose neighbours are already 8-adjacent,
    // so an L-step becomes one diagonal step (no staircase on the lattice).
    let cheb = |a: (i64, i64), b: (i64, i64)| (a.0 - b.0).abs().max((a.1 - b.1).abs());
    let mut taut: Vec<((i64, i64), f64)> = Vec::with_capacity(out.len());
    for item in out {
        while taut.len() >= 2 && cheb(taut[taut.len() - 2].0, item.0) <= 1 { taut.pop(); }
        taut.push(item);
    }
    taut
}

/// The chart inverse (same formula and operation order as vz-tilekey
/// `Chart::inverse_f64`, `vz.detmath/v1` trig) evaluated without the domain
/// check, for computation-halo nodes outside the declared chart domain.
pub fn inverse_any(decl: &ChartDeclV1, x: f64, y: f64) -> (f64, f64) {
    use vz_tilekey::chart::{e7_to_rad, TWO_PI};
    use vz_tilekey::detmath::{atan2, cos, sin, PI};
    let lat0 = e7_to_rad(decl.origin_lat_e7);
    let lon0 = e7_to_rad(decl.origin_lon_e7);
    let (s1, c1) = (sin(lat0), cos(lat0));
    let two_r = 2.0 * decl.datum_radius_mm as f64;
    if x == 0.0 && y == 0.0 { return (lat0, lon0); }
    let a = x / two_r;
    let b = y / two_r;
    let q = a * a + b * b;
    let d = 1.0 + q;
    let ve = (2.0 * a) / d;
    let vn = (2.0 * b) / d;
    let vu = (1.0 - q) / d;
    let px = c1 * vu - s1 * vn;
    let py = ve;
    let pz = s1 * vu + c1 * vn;
    let lat = atan2(pz, (px * px + py * py).sqrt());
    let mut lon = lon0 + atan2(py, px);
    if lon > PI { lon -= TWO_PI; }
    if lon <= -PI { lon += TWO_PI; }
    (lat, lon)
}

fn dist_point_segment(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 == 0.0 { 0.0 } else { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0) };
    let (qx, qy) = (a[0] + t * dx - p[0], a[1] + t * dy - p[1]);
    (qx * qx + qy * qy).sqrt()
}

pub fn segment_distance_mm(s: &Segment, n: usize) -> f64 {
    let (x, y) = node_mm(n);
    dist_point_segment([x as f64, y as f64], [s.start_mm[0] as f64, s.start_mm[1] as f64], [s.end_mm[0] as f64, s.end_mm[1] as f64])
}

fn bbox_hits_lattice(lat: &Lat, a: [f64; 2], b: [f64; 2], margin: f64) -> bool {
    let half = lat.half_mm();
    let (x0, x1) = (a[0].min(b[0]) - margin, a[0].max(b[0]) + margin);
    let (y0, y1) = (a[1].min(b[1]) - margin, a[1].max(b[1]) + margin);
    x1 >= -half && x0 <= half && y1 >= -half && y0 <= half
}

enum Walk { Done { path: Vec<(u32, f64)>, end: SegEnd, enters: bool, mouth_t: Option<f64>, forced: Option<String> }, Conflict(String) }

pub fn build(inp: &Inputs) -> R<Skeleton> {
    let t0 = inp.t0;
    let p = &inp.params;
    let thr = p.threshold_m3_s();
    let recipe = recipe_for(inp)?;
    let recipe_digest = recipe.digest()?;
    let chart = Chart::new(recipe.chart.clone()).map_err(|e| e.to_string())?;
    let site_mm = chart.forward_mm(inp.site.site_point_latlon_deg[0] * DEG2RAD, inp.site.site_point_latlon_deg[1] * DEG2RAD)
        .map_err(|e| e.to_string())?;

    // 1. Computation lattice (chart + halo); node -> T0 cell, node area.
    let lat = Lat::new(G_MIN - p.halo_nodes, G_MAX + p.halo_nodes);
    let nn = lat.nodes();
    let chart_lat = Lat::chart();
    let in_chart = |n: usize| { let (gx, gy) = lat.g(n); chart_lat.idx(gx, gy).is_some() };
    let mut t0_cell = vec![0u32; nn];
    let mut area = vec![0f64; nn];
    // Lane CLIMATE-RUNOFF: R (m/yr) at each computation node (key on only).
    let clim = match (&p.climate_runoff, &t0.climate) {
        (None, _) => None,
        (Some(_), Some(c)) => Some(c),
        (Some(_), None) => return Err("climate_runoff params set but no climate attached to T0".into()),
    };
    let mut r_node = vec![0f64; if clim.is_some() { nn } else { 0 }];
    for n in 0..nn {
        let (x, y) = lat.mm(n);
        let (la, lo) = if in_chart(n) { chart.inverse_mm(x, y).map_err(|e| e.to_string())? } else { inverse_any(&recipe.chart, x as f64, y as f64) };
        t0_cell[n] = t0.cell_of_deg(la * RAD2DEG, lo * RAD2DEG);
        if let Some(c) = clim { r_node[n] = c.r_at_deg(t0, la * RAD2DEG, lo * RAD2DEG); }
        let k = chart.scale_factor(x as f64, y as f64);
        area[n] = (SPACING_MM as f64 / 1000.0).powi(2) / (k * k);
    }
    let mut cell_nodes: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for n in 0..nn { cell_nodes.entry(t0_cell[n]).or_default().push(n); }

    // 2. Candidate cells: cells in the chart plus a 3-cell ring.
    let mut cand: BTreeSet<u32> = BTreeSet::new();
    for &c in cell_nodes.keys() {
        let (i, j) = t0.ij(c);
        for di in -3i64..=3 { for dj in -3i64..=3 {
            let ii = i as i64 + di;
            if ii < 0 || ii >= t0.nlat as i64 { continue; }
            let jj = (j as i64 + dj).rem_euclid(t0.nlon as i64);
            cand.insert((ii as usize * t0.nlon + jj as usize) as u32);
        } }
    }
    // Lane T0-IMPRINT-PROD (params.imprint, opt-in): the valley-axis junction of a cell; without the key the T0 cell centre.
    let jroot = match &p.imprint {
        Some(ip) => { ip.validate()?; Some(Scope::root(inp.site.seed).layer(&draw_digest(&recipe)?).and_then(|s| s.chart(0)).map_err(|e| e.to_string())?) }
        None => None,
    };
    let centre = |c: u32| -> R<[f64; 2]> {
        let (la, lo) = match (&p.imprint, &jroot) { (Some(ip), Some(r)) => junction_deg(t0, c, ip, r)?, _ => t0.centre_deg(c) };
        let v = chart.forward_f64(la * DEG2RAD, lo * DEG2RAD).map_err(|e| e.to_string())?;
        Ok(v)
    };

    // 3. Segments.
    let w = p.corridor_half_width_mm as f64;
    let margin = w + SPACING_MM as f64;
    let mut segs: Vec<Segment> = Vec::new();
    let mk = |kind: SegKind, address: String, from: u32, to: u32, a: [f64; 2], b: [f64; 2], ls: f64, le: f64, depth: u32| Segment {
        id: 0, kind, address, from_cell: from, to_cell: to,
        start_mm: [a[0].round() as i64, a[1].round() as i64], end_mm: [b[0].round() as i64, b[1].round() as i64],
        level_start_m: ls, level_end_m: le, q_t0_m3_s: t0.discharge_m3_s[from as usize] as f64,
        is_t0_river: (t0.discharge_m3_s[from as usize] as f64) >= thr && !t0.is_marine(from),
        t0_reach: t0.reach_id[from as usize], depth, amplitude_mm: 0, attempts: 0, straight_fallback: false,
        forced: false, conflict_notes: vec![], end: SegEnd::Outside, enters_chart: false, enters_computation_lattice: false, inflow_m3_s: 0.0,
        path: vec![], path_computation: vec![], path_nodes: 0, path_nodes_computation: 0, control_offsets_permille: vec![],
    };
    let mut seg_points: Vec<([f64; 2], [f64; 2])> = Vec::new();
    for &c in &cand {
        if t0.is_marine(c) { continue; }
        let pc = centre(c)?;
        let depth = t0.chain_from(c)?.len() as u32;
        match t0.flow(c) {
            Some(r) => {
                let pr = centre(r)?;
                if bbox_hits_lattice(&lat, pc, pr, margin) {
                    segs.push(mk(SegKind::T0Flow, format!("seg/{c}-{r}"), c, r, pc, pr, t0.level(c), t0.level(r), depth));
                    seg_points.push((pc, pr));
                }
                if t0.river_inflow[c as usize] == 0 {
                    let f = p.head_extension_permille as f64 / 1000.0;
                    let mut h = [pc[0] - f * (pr[0] - pc[0]), pc[1] - f * (pr[1] - pc[1])];
                    // Lane T0-IMPRINT-PROD (params.imprint only): a displaced junction can push the head extension (and its
                    // meander) across the cell edge, where it captures a foreign cell's hillslopes; clip it to its own cell.
                    let mut keep = true;
                    if p.imprint.is_some() {
                        let k = head_clip_fraction(&recipe.chart, &chart, t0, c, pc, h, p.amplitude_mm() + SPACING_MM as f64);
                        h = [pc[0] + k * (h[0] - pc[0]), pc[1] + k * (h[1] - pc[1])];
                        keep = ((h[0] - pc[0]).powi(2) + (h[1] - pc[1]).powi(2)).sqrt() >= 2.0 * SPACING_MM as f64;
                    }
                    if keep && bbox_hits_lattice(&lat, h, pc, margin) {
                        let mean = t0.elevation_m[c as usize] as f64;
                        segs.push(mk(SegKind::HeadExtension, format!("head/{c}"), c, c, h, pc, mean.max(t0.level(c)), t0.level(c), depth));
                        seg_points.push((h, pc));
                    }
                }
            }
            None => {
                if bbox_hits_lattice(&lat, pc, pc, 0.0) {
                    segs.push(mk(SegKind::Terminal, format!("term/{c}"), c, c, pc, pc, t0.level(c), t0.level(c), depth));
                    seg_points.push((pc, pc));
                }
            }
        }
    }
    // Placement order: terminals, then T0 flow edges downstream-first (depth
    // ascending, larger T0 discharge first, then cell id), then heads.
    let mut order: Vec<usize> = (0..segs.len()).collect();
    order.sort_by(|&a, &b| {
        let (sa, sb) = (&segs[a], &segs[b]);
        sa.kind.cmp(&sb.kind)
            .then(sa.depth.cmp(&sb.depth))
            .then(sb.q_t0_m3_s.total_cmp(&sa.q_t0_m3_s))
            .then(sa.from_cell.cmp(&sb.from_cell))
            .then(sa.to_cell.cmp(&sb.to_cell))
    });
    let segs: Vec<Segment> = order.iter().enumerate().map(|(i, &k)| { let mut s = segs[k].clone(); s.id = i as u32; s }).collect();
    let seg_points: Vec<([f64; 2], [f64; 2])> = order.iter().map(|&k| seg_points[k]).collect();
    let mut segs = segs;
    // Entry cell e(s) and the T0 chain from it.
    let entry = |s: &Segment| if s.kind == SegKind::T0Flow { s.to_cell } else { s.from_cell };
    let mut chain_set: Vec<BTreeSet<u32>> = Vec::with_capacity(segs.len());
    for s in &segs { chain_set.push(t0.chain_from(entry(s))?.into_iter().collect()); }
    let allowed = |s: &Segment, t: &Segment, cs: &BTreeSet<u32>| -> bool {
        cs.contains(&t.from_cell) || entry(t) == entry(s) || entry(t) == s.from_cell
    };

    // 4. Placement.
    // Lane CLIMATE-RUNOFF: draws keyed by the climate-free recipe (= recipe_digest without the key).
    let ddig = draw_digest(&recipe)?;
    let root = Scope::root(inp.site.seed).layer(&ddig).and_then(|s| s.chart(0)).map_err(|e| e.to_string())?;
    let mut owner = vec![NONE; nn];
    let mut receiver = vec![NONE; nn];
    let mut seg_t = vec![f64::NAN; nn];
    let mut kind = vec![0u8; nn];
    for n in 0..nn {
        kind[n] = match t0.water(t0_cell[n]) { CellWater::Marine => 2, CellWater::Lake => 4, CellWater::Land => 0 };
    }
    let mut mouth_t: Vec<Option<f64>> = vec![None; segs.len()];
    let mut stats = BuildStats { site_in_chart_mm_recomputed: site_mm, ..Default::default() };
    let amp = p.amplitude_mm();
    for si in 0..segs.len() {
        let scope = root.address(&segs[si].address).map_err(|e| e.to_string())?;
        let (a, b) = seg_points[si];
        let max_att = p.planform_max_attempts;
        let mut placed = false;
        for attempt in 0..=max_att {
            let straight = attempt == max_att || amp == 0.0;
            let a_mm = if straight { 0.0 } else { amp };
            let (samples, ctrl) = curve(a, b, a_mm, p, &scope, attempt)?;
            let raster = rasterize(&samples);
            let last_try = attempt == max_att;
            match walk(&lat, &segs, si, &raster, &owner, &kind, &t0_cell, t0, &chain_set, &allowed, last_try) {
                Walk::Conflict(why) => {
                    segs[si].conflict_notes.push(format!("attempt {attempt}: {why}"));
                    continue;
                }
                Walk::Done { path, end, enters, mouth_t: mt, forced } => {
                    let s = &mut segs[si];
                    s.attempts = attempt + 1;
                    s.straight_fallback = straight && amp > 0.0;
                    s.amplitude_mm = a_mm.round() as i64;
                    s.control_offsets_permille = ctrl;
                    s.enters_computation_lattice = enters;
                    mouth_t[si] = mt;
                    for (k, &(n, t)) in path.iter().enumerate() {
                        owner[n as usize] = s.id;
                        seg_t[n as usize] = t;
                        kind[n as usize] = 1;
                        receiver[n as usize] = if k + 1 < path.len() { path[k + 1].0 } else {
                            match &end { SegEnd::Merge { node, .. } | SegEnd::Mouth { node } => *node,
                                SegEnd::OffChart => OFF_CHART, _ => NONE }
                        };
                    }
                    if let SegEnd::Mouth { node } = end { kind[node as usize] = 3; }
                    if let Some(why) = forced { s.forced = true; s.conflict_notes.push(format!("forced: {why}")); }
                    s.path = path.iter().map(|x| x.0).collect();
                    s.path_nodes = s.path.len() as u32;
                    s.end = end;
                    placed = true;
                    break;
                }
            }
        }
        if !placed { return Err(format!("segment {} could not be placed", segs[si].address)); }
        let s = &segs[si];
        if s.attempts > 1 { stats.segments_retried += 1; }
        if s.straight_fallback { stats.segments_straight_fallback += 1; }
        if s.forced { stats.segments_forced += 1; }
        match s.kind { SegKind::T0Flow => stats.segments_t0_flow += 1, SegKind::HeadExtension => stats.segments_head += 1, SegKind::Terminal => stats.segments_terminal += 1 }
    }

    // 5. Bed levels on axis nodes, then carve monotone along flow.
    let mut bed = vec![f64::NAN; nn];
    for (si, s) in segs.iter().enumerate() {
        for &n in &s.path {
            let t = seg_t[n as usize];
            let v = match (s.kind, mouth_t[si]) {
                (SegKind::T0Flow, Some(tm)) if tm > 0.0 => s.level_start_m + (t0.sea_level_m - s.level_start_m) * (t / tm).min(1.0),
                (SegKind::Terminal, _) => s.level_start_m,
                _ => s.level_start_m + (s.level_end_m - s.level_start_m) * t,
            };
            bed[n as usize] = v;
        }
        if let SegEnd::Mouth { node } = s.end { bed[node as usize] = t0.sea_level_m; }
    }
    let topo = topo_order(&receiver)?;
    for &n in &topo {
        let r = receiver[n];
        if r >= OFF_CHART || bed[n].is_nan() || bed[r as usize].is_nan() { continue; }
        if bed[r as usize] > bed[n] {
            let carve = bed[r as usize] - bed[n];
            stats.bed_nodes_carved += 1;
            if carve > stats.bed_max_carve_m { stats.bed_max_carve_m = carve; }
            bed[r as usize] = bed[n];
        }
    }

    // 6. Hillslope routing: per T0 cell shortest path to its own axes.
    let mut own_segs: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for s in &segs {
        own_segs.entry(s.from_cell).or_default().push(s.id);
        if entry(s) != s.from_cell { own_segs.entry(entry(s)).or_default().push(s.id); }
    }
    let mut join_node = vec![NONE; nn];
    let mut join_seg = vec![NONE; nn];
    let mut fallback = vec![false; nn];
    let mut dist = vec![u64::MAX; nn];
    let routable = |n: usize, kind: &[u8], owner: &[u32]| (kind[n] == 0 || kind[n] == 4) && owner[n] == NONE;
    let is_axis = |n: usize, owner: &[u32]| owner[n] != NONE;
    for (&c, nodes) in &cell_nodes {
        if t0.is_marine(c) { continue; }
        let mut seeds: Vec<usize> = Vec::new();
        if let Some(list) = own_segs.get(&c) {
            for &sid in list {
                let s = &segs[sid as usize];
                seeds.extend(s.path.iter().map(|&n| n as usize));
                if s.path.is_empty() { if let SegEnd::Merge { node, .. } = s.end { seeds.push(node as usize); } }
            }
        }
        let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
        let mut seed_set: BTreeSet<usize> = BTreeSet::new();
        if seeds.is_empty() {
            // The cell's own channel lies outside the computation lattice: drain across its edge.
            let out_seg = segs.iter().find(|s| s.kind == SegKind::T0Flow && s.from_cell == c).map(|s| s.id).unwrap_or(NONE);
            for &n in nodes {
                if lat.on_edge(n) && routable(n, &kind, &owner) {
                    dist[n] = 0;
                    receiver[n] = OFF_CHART;
                    join_seg[n] = out_seg;
                    fallback[n] = true;
                    stats.hillslope_offchart_fallback += 1;
                    heap.push(Reverse((0, n)));
                }
            }
        } else {
            for n in seeds { if seed_set.insert(n) { heap.push(Reverse((0, n))); } }
        }
        let mut done: BTreeSet<usize> = BTreeSet::new();
        while let Some(Reverse((d, u))) = heap.pop() {
            let is_seed = seed_set.contains(&u);
            if !is_seed && d > dist[u] { continue; }
            if !done.insert(u) { continue; }
            if !is_seed && receiver[u] != OFF_CHART {
                let par = receiver[u] as usize;
                if is_axis(par, &owner) || seed_set.contains(&par) {
                    join_node[u] = par as u32;
                    join_seg[u] = if owner[par] != NONE { owner[par] } else { NONE };
                } else {
                    join_node[u] = join_node[par];
                    join_seg[u] = join_seg[par];
                    fallback[u] = fallback[par];
                }
                stats.hillslope_pass1 += 1;
            }
            let (ux, uy) = lat.g(u);
            for (k, &(dx, dy)) in D8.iter().enumerate() {
                let Some(v) = lat.idx(ux + dx, uy + dy) else { continue };
                if t0_cell[v] != c || !routable(v, &kind, &owner) || done.contains(&v) { continue; }
                let diag = k % 2 == 1;
                if diag {
                    let a = lat.idx(ux + dx, uy).unwrap();
                    let b = lat.idx(ux, uy + dy).unwrap();
                    if is_axis(a, &owner) && is_axis(b, &owner) { continue; }
                }
                let nd = d + if diag { STEP_DIAG_MM } else { STEP_ORTHO_MM };
                if nd < dist[v] {
                    dist[v] = nd;
                    receiver[v] = u as u32;
                    heap.push(Reverse((nd, v)));
                }
            }
        }
    }
    // Pass 2: land nodes not reached inside their own cell (cut off by a
    // foreign axis): shortest path to any routed node. Counted against F3.
    let pending: Vec<bool> = (0..nn).map(|n| routable(n, &kind, &owner) && receiver[n] == NONE).collect();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
    let mut d2 = vec![u64::MAX; nn];
    for n in 0..nn {
        if !pending[n] { continue; }
        let (x, y) = lat.g(n);
        for &(dx, dy) in &D8 {
            if let Some(v) = lat.idx(x + dx, y + dy) {
                let routed = is_axis(v, &owner) || (routable(v, &kind, &owner) && !pending[v]);
                if routed && d2[v] != 0 { d2[v] = 0; heap.push(Reverse((0, v))); }
            }
        }
    }
    let mut done2 = vec![false; nn];
    while let Some(Reverse((d, u))) = heap.pop() {
        if d > d2[u] || done2[u] { continue; }
        done2[u] = true;
        if pending[u] {
            let par = receiver[u] as usize;
            if is_axis(par, &owner) { join_node[u] = par as u32; join_seg[u] = owner[par]; }
            else { join_node[u] = join_node[par]; join_seg[u] = join_seg[par]; }
            dist[u] = d;
            stats.hillslope_pass2_foreign += 1;
        }
        let (ux, uy) = lat.g(u);
        for (k, &(dx, dy)) in D8.iter().enumerate() {
            let Some(v) = lat.idx(ux + dx, uy + dy) else { continue };
            if !pending[v] || done2[v] { continue; }
            let diag = k % 2 == 1;
            if diag {
                let a = lat.idx(ux + dx, uy).unwrap();
                let b = lat.idx(ux, uy + dy).unwrap();
                if is_axis(a, &owner) && is_axis(b, &owner) { continue; }
            }
            let nd = d + if diag { STEP_DIAG_MM } else { STEP_ORTHO_MM };
            if nd < d2[v] { d2[v] = nd; receiver[v] = u as u32; heap.push(Reverse((nd, v))); }
        }
    }
    // Pass 3: still unrouted land next to the sea drains into it.
    for n in 0..nn {
        if !(routable(n, &kind, &owner) && receiver[n] == NONE) { continue; }
        let (x, y) = lat.g(n);
        let sea = D8.iter().filter_map(|&(dx, dy)| lat.idx(x + dx, y + dy)).filter(|&v| kind[v] == 2 || kind[v] == 3).min();
        if let Some(v) = sea { receiver[n] = v as u32; stats.hillslope_pass3_coast += 1; } else { stats.hillslope_unrouted += 1; }
    }

    // 7. Discharge: uniform T0 runoff x node area, plus inflow where a
    // channel enters the chart.
    let mut cells: BTreeMap<u32, CellInfo> = BTreeMap::new();
    let mut area_in: BTreeMap<u32, f64> = BTreeMap::new();
    for (&c, nodes) in &cell_nodes { area_in.insert(c, nodes.iter().map(|&n| area[n]).sum()); }
    // Lane CLIMATE-RUNOFF: chart-entry inflow = climate T0 discharge of the start cell minus the climate runoff of its
    // in-lattice nodes; T0 discharge with the lake rule (climate_runoff.rs).
    let (q_t0_clim, mut climate_lakes) = match clim { Some(c) => { let (q, l) = c.t0_discharge(t0)?; (Some(q), l) } None => (None, Vec::new()) };
    climate_lakes.retain(|l| cell_nodes.keys().any(|&c| t0.water_body_id[c as usize] == l.body_id));
    let mut runoff_in: BTreeMap<u32, f64> = BTreeMap::new();
    if clim.is_some() { for (&c, nodes) in &cell_nodes { runoff_in.insert(c, nodes.iter().map(|&n| r_node[n] * area[n]).sum()); } }
    let mut inflow = vec![0f64; nn];
    for s in segs.iter_mut() {
        if !s.enters_computation_lattice || s.kind != SegKind::T0Flow { continue; }
        let target = match (s.path.first(), &s.end) {
            (Some(&n), _) => n as usize,
            (None, SegEnd::Merge { node, .. }) => *node as usize,
            _ => continue,
        };
        if let Some(qc) = &q_t0_clim {
            let own = runoff_in.get(&s.from_cell).copied().unwrap_or(0.0) / SECONDS_PER_YEAR;
            s.inflow_m3_s = (qc[s.from_cell as usize] - own).max(0.0);
            inflow[target] += s.inflow_m3_s;
            continue;
        }
        let own = t0.runoff_m_per_yr * area_in.get(&s.from_cell).copied().unwrap_or(0.0) / SECONDS_PER_YEAR;
        s.inflow_m3_s = (s.q_t0_m3_s - own).max(0.0);
        inflow[target] += s.inflow_m3_s;
    }
    let mut q = vec![0f64; nn];
    let topo = topo_order(&receiver)?;
    for &n in &topo {
        if kind[n] != 2 { q[n] += inflow[n]; }
        if kind[n] == 0 || kind[n] == 1 || kind[n] == 4 {
            if clim.is_some() { q[n] += r_node[n] * area[n] / SECONDS_PER_YEAR; } else { q[n] += t0.runoff_m_per_yr * area[n] / SECONDS_PER_YEAR; }
        }
        let r = receiver[n];
        if r < OFF_CHART { q[r as usize] += q[n]; }
    }

    // 8. Derived per-node fields.
    let mut trunk = vec![f64::NAN; nn];
    let mut basin = vec![NONE; nn];
    let mut reach = vec![NONE; nn];
    let mut final_basin_of_seg: Vec<u32> = Vec::with_capacity(segs.len());
    for s in &segs { final_basin_of_seg.push(t0.final_basin(entry(s))?); }
    for n in 0..nn {
        let c = t0_cell[n];
        match kind[n] {
            2 | 3 => { trunk[n] = t0.sea_level_m; basin[n] = t0.basin_id[c as usize]; }
            1 => {
                let s = owner[n] as usize;
                trunk[n] = bed[n];
                basin[n] = final_basin_of_seg[s];
                reach[n] = segs[s].t0_reach;
            }
            _ => {
                let j = join_node[n];
                trunk[n] = if j != NONE { bed[j as usize] } else { t0.level(c) };
                let s = join_seg[n];
                if s != NONE { basin[n] = final_basin_of_seg[s as usize]; reach[n] = segs[s as usize].t0_reach; }
                else { basin[n] = t0.final_basin(c)?; reach[n] = t0.reach_id[c as usize]; }
            }
        }
    }
    for n in 0..nn { if kind[n] == 3 { bed[n] = t0.sea_level_m; } }

    // 9. Crop to the published chart lattice.
    let cmap = |ci: u32| -> u32 {
        if ci >= OFF_CHART { return ci; }
        let (gx, gy) = lat.g(ci as usize);
        chart_lat.idx(gx, gy).map(|k| k as u32).unwrap_or(OFF_CHART)
    };
    let cidx: Vec<usize> = (0..NODES).map(|k| { let (gx, gy) = chart_lat.g(k); lat.idx(gx, gy).unwrap() }).collect();
    let crop_u32 = |v: &[u32]| -> Vec<u32> { cidx.iter().map(|&i| v[i]).collect() };
    let crop_f64 = |v: &[f64]| -> Vec<f64> { cidx.iter().map(|&i| v[i]).collect() };
    for s in segs.iter_mut() {
        s.path_computation = std::mem::take(&mut s.path);
        s.path_nodes_computation = s.path_computation.len() as u32;
        s.path = s.path_computation.iter().map(|&n| cmap(n)).filter(|&k| k < OFF_CHART).collect();
        s.path_nodes = s.path.len() as u32;
        let first_inside = s.path_computation.iter().position(|&n| in_chart(n as usize));
        s.enters_chart = match first_inside { Some(0) => s.enters_computation_lattice, Some(_) => true, None => false };
        s.end = match s.end.clone() {
            SegEnd::Merge { node, into_segment } => SegEnd::Merge { node: cmap(node), into_segment },
            SegEnd::Mouth { node } => SegEnd::Mouth { node: cmap(node) },
            e => e,
        };
        if s.enters_chart && s.kind == SegKind::T0Flow && s.path.is_empty() {
            // Its in-chart part is a mouth or merge only.
            if let SegEnd::Mouth { node } = s.end { if node >= OFF_CHART { s.end = SegEnd::Outside; } }
        }
        if !s.enters_chart && s.path.is_empty() && !matches!(s.end, SegEnd::Mouth { node } if node < OFF_CHART) { s.end = SegEnd::Outside; }
    }
    let o_t0_cell = crop_u32(&t0_cell);
    let o_area = crop_f64(&area);
    let o_kind: Vec<u8> = cidx.iter().map(|&i| kind[i]).collect();
    let o_owner = crop_u32(&owner);
    let o_receiver: Vec<u32> = cidx.iter().map(|&i| cmap(receiver[i])).collect();
    let o_join_node: Vec<u32> = cidx.iter().map(|&i| if join_node[i] == NONE { NONE } else { cmap(join_node[i]) }).collect();
    let o_join_seg = crop_u32(&join_seg);
    let o_fallback: Vec<bool> = cidx.iter().map(|&i| fallback[i]).collect();
    let o_dist: Vec<f64> = cidx.iter().map(|&i| if dist[i] == u64::MAX { f64::NAN } else { dist[i] as f64 / 1000.0 }).collect();
    stats.halo_nodes = p.halo_nodes;
    stats.computation_posts_per_axis = lat.n as u32;
    stats.hillslope_offchart_fallback_in_chart = o_fallback.iter().filter(|&&f| f).count() as u32;

    // 10. Cell table (chart nodes only).
    let mut chart_cell_nodes: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for k in 0..NODES { chart_cell_nodes.entry(o_t0_cell[k]).or_default().push(k); }
    let mut cell_ids: BTreeSet<u32> = chart_cell_nodes.keys().copied().collect();
    for s in segs.iter().filter(|s| !s.path.is_empty()) { cell_ids.insert(s.from_cell); cell_ids.insert(s.to_cell); }
    for c in cell_ids {
        let (i, j) = t0.ij(c);
        let pj = centre(c)?;
        let pc = if p.imprint.is_some() { let (la, lo) = t0.centre_deg(c); chart.forward_f64(la * DEG2RAD, lo * DEG2RAD).map_err(|e| e.to_string())? } else { pj };
        let ai: f64 = chart_cell_nodes.get(&c).map(|v| v.iter().map(|&k| o_area[k]).sum()).unwrap_or(0.0);
        let ca = t0.cell_area_m2(c);
        cells.insert(c, CellInfo {
            cell: c, ij: [i, j], centre_mm: [pc[0].round() as i64, pc[1].round() as i64],
            water: format!("{:?}", t0.water(c)), elevation_m: t0.elevation_m[c as usize] as f64,
            subgrid_relief_m: t0.subgrid_relief_m[c as usize] as f64, level_m: t0.level(c), flow_to: t0.flow(c),
            t0_reach: opt(t0.reach_id[c as usize]), basin: opt(t0.basin_id[c as usize]), final_basin: opt(t0.final_basin(c)?),
            q_t0_m3_s: t0.discharge_m3_s[c as usize] as f64, river_inflow: t0.river_inflow[c as usize],
            nodes_in_chart: chart_cell_nodes.get(&c).map(|v| v.len() as u32).unwrap_or(0), area_in_chart_m2: ai, cell_area_m2: ca,
            inside_permille: ((ai / ca) * 1000.0).round() as i64,
            junction_mm: p.imprint.as_ref().map(|_| [pj[0].round() as i64, pj[1].round() as i64]),
        });
    }

    Ok(Skeleton {
        recipe, recipe_digest, chart, params: p.clone(), t0_cell: o_t0_cell, area_m2: o_area, kind: o_kind, owner: o_owner,
        receiver: o_receiver, join_node: o_join_node, join_segment: o_join_seg, offchart_fallback: o_fallback,
        bed_m: crop_f64(&bed), trunk_level_m: crop_f64(&trunk), flow_dist_m: o_dist,
        discharge_m3_s: crop_f64(&q), basin: crop_u32(&basin), t0_reach: crop_u32(&reach), segments: segs, cells, stats,
        site: inp.site.clone(), sea_level_m: t0.sea_level_m, threshold_m3_s: thr,
        climate_runoff_m_per_yr: if clim.is_some() { Some(crop_f64(&r_node)) } else { None },
        climate_twin_digest: None, physical_climate_twin_digest: None, climate_lakes,
    })
}

/// Kahn order over the receiver forest (upstream first). Errors on a cycle.
pub fn topo_order(receiver: &[u32]) -> R<Vec<usize>> {
    let n = receiver.len();
    let mut indeg = vec![0u32; n];
    for &r in receiver { if r < OFF_CHART { indeg[r as usize] += 1; } }
    let mut q: VecDeque<usize> = (0..n).filter(|&i| indeg[i] == 0).collect();
    let mut out = Vec::with_capacity(n);
    while let Some(i) = q.pop_front() {
        out.push(i);
        let r = receiver[i];
        if r < OFF_CHART {
            indeg[r as usize] -= 1;
            if indeg[r as usize] == 0 { q.push_back(r as usize); }
        }
    }
    if out.len() != n { return Err(format!("receiver graph has a cycle ({} of {n} nodes ordered)", out.len())); }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn walk(lat: &Lat, segs: &[Segment], si: usize, raster: &[((i64, i64), f64)], owner: &[u32], kind: &[u8], t0_cell: &[u32], t0: &T0,
    chain_set: &[BTreeSet<u32>], allowed: &dyn Fn(&Segment, &Segment, &BTreeSet<u32>) -> bool, force: bool) -> Walk {
    let s = &segs[si];
    let cs = &chain_set[si];
    let mut path: Vec<(u32, f64)> = Vec::new();
    let mut inside = false;
    let mut enters = false;
    let mut forced_note: Option<String> = None;
    let conflict = |why: String, forced_note: &mut Option<String>| -> Option<Walk> {
        if force { *forced_note = Some(why); None } else { Some(Walk::Conflict(why)) }
    };
    for (k, &((gx, gy), t)) in raster.iter().enumerate() {
        let Some(n) = lat.idx(gx, gy) else {
            if inside {
                return finish(path, SegEnd::OffChart, enters, None, forced_note);
            }
            continue;
        };
        if !inside { inside = true; enters = k > 0; }
        if s.kind == SegKind::Terminal {
            path.push((n as u32, t));
            return finish(path, SegEnd::Terminal, enters, None, forced_note);
        }
        // Sea.
        if kind[n] == 2 || kind[n] == 3 {
            let m = t0_cell[n];
            let ok = s.kind == SegKind::T0Flow && t0.is_marine(s.to_cell)
                && (m == s.to_cell || (t0.basin_id[m as usize] != NONE && t0.basin_id[m as usize] == t0.basin_id[s.to_cell as usize]));
            if !ok {
                if let Some(w) = conflict(format!("reaches the sea in cell {m} at node {n}"), &mut forced_note) { return w; }
            }
            return finish(path, SegEnd::Mouth { node: n as u32 }, enters, Some(t), forced_note);
        }
        // Occupied node.
        if owner[n] != NONE {
            let o = &segs[owner[n] as usize];
            if !allowed(s, o, cs) {
                if let Some(w) = conflict(format!("hits {} at node {n}", o.address), &mut forced_note) { return w; }
            }
            return finish(path, SegEnd::Merge { node: n as u32, into_segment: o.id }, enters, None, forced_note);
        }
        // Adjacent contacts.
        let mut ok_adj: Option<usize> = None;
        for &(dx, dy) in &D8 {
            let Some(v) = lat.idx(gx + dx, gy + dy) else { continue };
            if owner[v] == NONE { continue; }
            let o = &segs[owner[v] as usize];
            if allowed(s, o, cs) {
                if ok_adj.is_none_or(|x| v < x) { ok_adj = Some(v); }
            } else if let Some(w) = conflict(format!("touches {} at node {v}", o.address), &mut forced_note) {
                return w;
            } else {
                ok_adj = Some(ok_adj.map_or(v, |x| x.min(v)));
            }
        }
        path.push((n as u32, t));
        if let Some(v) = ok_adj {
            if path.len() > 1 || s.kind != SegKind::T0Flow || enters {
                return finish(path, SegEnd::Merge { node: v as u32, into_segment: owner[v] }, enters, None, forced_note);
            }
        }
    }
    if inside { finish(path, SegEnd::Dangling, enters, None, forced_note) } else { finish(path, SegEnd::Outside, enters, None, forced_note) }
}

fn finish(path: Vec<(u32, f64)>, end: SegEnd, enters: bool, mouth_t: Option<f64>, forced: Option<String>) -> Walk {
    Walk::Done { path, end, enters, mouth_t, forced }
}
