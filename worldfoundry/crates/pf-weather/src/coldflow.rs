//! COLDPOOL-2 "coldflow" (lane WEATHER-V1, 2026-10-21; opt-in `coldflow=1`):
//! cold-air drainage over the chart-wide river network, pools per connected
//! airshed with a shared lid, spill over sills, outlets to the sea.
//!
//! A derived presentation/perception layer (never canon by itself), like
//! COLDPOOL-1, but the pool is a flow, not a floor-fill:
//!
//! 1. Network (chart-wide, the layer the T1b windows came from): the
//!    `t1b.rivernet` lod-8 lattice (513 x 513 nodes, 512 m): receiver, bed,
//!    node position, channel flag, kind; the T1a channel seeds (kind 1) take
//!    their receivers from the T1a skeleton. Beds never rise along flow.
//! 2. Production and drainage: every land node makes cold air at
//!    p = p_max Q per unit area (Q = night cooling from WEATHER, as COLDPOOL-1);
//!    the volume flux F accumulates down the receivers (confluences add).
//! 3. Drainage depth at a node: the steady depth h of a gravity current that
//!    carries F: F = u(h) A(h), u = sqrt(g' h S / Cd) (slope-driven, bulk
//!    drag), g' = g dT / T_ref with dT = the inversion strength. A(h) is the
//!    valley cross-section: a floor of minimum width, hillslope sides at a
//!    fixed tangent, plus the network nodes of the 9 x 9 stencil whose bed is
//!    within h of this bed (each adds ~one lattice cell of floor width):
//!    broad flat floors are wide, gorges are narrow (constrictions back the
//!    cold air up). At an outlet to the sea the depth is at most the critical
//!    depth (free outflow: F^2 = g' A^3 / W), so valleys open to the sea drain.
//!    Wind gate: h is scaled by 1 - smoothstep(2.5, 3.5, wind at formation).
//! 4. Lids (airsheds): from the outlets upstream, lid(i) = max(bed(i) + h(i),
//!    lid(receiver)). A node whose own flow surface is below its receiver's lid
//!    is inside the downstream pool and shares that lid (connected valleys
//!    share a lid; a pool behind a sill/constriction spills over it and its lid
//!    propagates up every tributary until the beds rise above it). A high
//!    steep valley keeps only its thin drainage layer: its cold air goes down.
//!    The source node of a lid (where bed + h set it) is the pool's floor and
//!    identity.
//! 5. Fog: pool depth D = h(source), inversion dT = dT_max Q (as COLDPOOL-1),
//!    dewpoint depression DD at the lid elevation from WEATHER; fog top
//!    F = lid - D DD / dT (whole pool when DD <= 0, none when DD >= dT), capped
//!    at the pool floor + `fog_depth_max_m` (200 m). Only
//!    stagnant pooled air fogs: lid - bed >= pool_depth_min AND the drainage
//!    speed at the airshed's source node is at most `stagnant_speed_m_s` (the
//!    whole airshed shares the flag). Faster or thinner layers are moving,
//!    mixed drainage flow (mask 4, no fog).
//! 6. Burn-off after sunrise: lid and fog top fall at the COLDPOOL-1 rate; per
//!    post the fog top falls faster on sun-facing slopes and slower in shade
//!    (factor cos(incidence) / sin(sun elevation), clamped), so edges and sunny
//!    slopes clear first and deep shaded basins last.
//! 7. Steam fog (mask 5): over enclosed water (lakes; sea posts with at least
//!    `steam_enclosed_frac` land within 0.5 km: inlets, estuaries) when the
//!    water is warmer than the drained cold air by >= `steam_dt_min_k`.
//!
//! Posts: each post follows steepest descent on the window heights for at most
//! 32 posts (or to a channel / water post); the end point is matched to the
//! nearest network channel link within 2 lattice nodes; the lid and fog top are
//! interpolated along that link. Every per-post stencil (descent 32, aspect 1,
//! enclosure 16) fits inside the 64-post halo and the network is chart-wide,
//! so windows agree exactly at shared posts.
//!
//! Bound: the airshed is computed on the whole chart (513 x 513 nodes,
//! +-131 km). Cold air entering the chart from outside is not modelled (the
//! rivernet has no off-chart inflow for land nodes; T1a seeds leaving the chart
//! are outlets with no downstream lid).

use std::collections::{BTreeMap, VecDeque};
use std::path::{Path, PathBuf};

use serde::Serialize;
use vz_tilekey::chart::{Chart, ChartDeclV1};

use crate::constants::ConstantRow;
use crate::coldpool::{night_summary, NightSummary, Window};
use crate::state::{smoothstep, WeatherModel};
use crate::sun::{solar_clock, DEG};

pub const SCHEMA: &str = "pf.coldpool-window/v1";
pub const MODEL_ID: &str = "pf-coldflow-v1";
pub const NONE: u32 = u32::MAX;
pub const OFF_CHART: u32 = u32::MAX - 1;

#[derive(Debug, Clone, PartialEq)]
pub struct ColdFlowConstantsV1 {
    pub production_max_m_s: f64,
    pub gravity_m_s2: f64,
    pub t_ref_k: f64,
    pub dt_floor_k: f64,
    pub drag_cd: f64,
    pub slope_min: f64,
    pub side_slope: f64,
    pub floor_width_min_m: f64,
    pub width_stencil_nodes: f64,
    pub depth_cap_m: f64,
    pub stagnant_speed_m_s: f64,
    pub fog_depth_max_m: f64,
    pub wind_gate_lo_m_s: f64,
    pub wind_gate_hi_m_s: f64,
    pub descent_max_posts: f64,
    pub link_search_nodes: f64,
    pub aspect_factor_min: f64,
    pub aspect_factor_max: f64,
    pub steam_dt_min_k: f64,
    pub steam_dt_full_k: f64,
    pub steam_depth_max_m: f64,
    pub steam_water_lag_days: f64,
    pub steam_enclosed_frac: f64,
    pub steam_enclosed_r_posts: f64,
}

impl Default for ColdFlowConstantsV1 {
    fn default() -> Self {
        ColdFlowConstantsV1 {
            production_max_m_s: 0.0015,
            gravity_m_s2: 9.81,
            t_ref_k: 283.15,
            dt_floor_k: 1.0,
            drag_cd: 0.1,
            slope_min: 0.003,
            side_slope: 0.4,
            floor_width_min_m: 30.0,
            width_stencil_nodes: 4.0,
            depth_cap_m: 400.0,
            stagnant_speed_m_s: 2.0,
            fog_depth_max_m: 200.0,
            wind_gate_lo_m_s: 2.5,
            wind_gate_hi_m_s: 3.5,
            descent_max_posts: 32.0,
            link_search_nodes: 2.0,
            aspect_factor_min: 0.4,
            aspect_factor_max: 2.5,
            steam_dt_min_k: 6.0,
            steam_dt_full_k: 10.0,
            steam_depth_max_m: 40.0,
            steam_water_lag_days: 30.0,
            steam_enclosed_frac: 0.65,
            steam_enclosed_r_posts: 16.0,
        }
    }
}

const R: &str = "recalled";
const P: &str = "project_assumption";

impl ColdFlowConstantsV1 {
    pub fn table(&self) -> Vec<ConstantRow> {
        let k = self;
        let row = |name, value, unit, class, basis| ConstantRow { name, value, unit, class, basis };
        vec![
            row("production_max_m_s", k.production_max_m_s, "m/s", P, "cold-air volume exported downslope per unit area at Q = 1: surface sensible heat loss ~30 W/m2 (recalled order 20-50) over rho cp dT_ref (1.2 kg/m3, 1005 J/kg/K, 5 K) = 5e-3 m/s, times an export fraction 0.3 (the rest stays in the local slope inversion); 5e-3 gave 130-200 m deep, ~5 m/s drainage flows in site-1 valleys (lane figure cf_s1_draft1_section)"),
            row("gravity_m_s2", k.gravity_m_s2, "m/s2", R, "standard gravity"),
            row("t_ref_k", k.t_ref_k, "K", P, "reference absolute temperature in the reduced gravity g' = g dT / T_ref"),
            row("dt_floor_k", k.dt_floor_k, "K", P, "minimum temperature deficit used in g' (avoids a vanishing g' at tiny Q)"),
            row("drag_cd", k.drag_cd, "1", P, "bulk drag (with entrainment) of the drainage flow: u = sqrt(g' h S / Cd); chosen so valley drainage flows run ~1-3 m/s (recalled order of katabatic speeds)"),
            row("slope_min", k.slope_min, "1", P, "floor of the network bed slope in the flow law (flats, basins)"),
            row("side_slope", k.side_slope, "1", P, "valley side tangent bounding the cross-section above the floor (the 512 m network does not resolve hillslopes)"),
            row("floor_width_min_m", k.floor_width_min_m, "m", P, "minimum valley floor width"),
            row("width_stencil_nodes", k.width_stencil_nodes, "nodes", P, "half-width of the network stencil (9 x 9, +-2 km) whose channel nodes with bed within h of this bed widen the floor by spacing^2 / ((2r+1) spacing) each (5 x 5 under-counted the 5 km wide site-5 floor: 400 m capped pools)"),
            row("depth_cap_m", k.depth_cap_m, "m", P, "cap of the drainage/pool depth (deep valley inversions 300-500 m, recalled)"),
            row("stagnant_speed_m_s", k.stagnant_speed_m_s, "m/s", P, "an airshed holds stagnant pooled air (fog possible) when the drainage speed at its source node (where bed + h sets the lid) is at most this; every node sharing that lid shares the flag (a per-node flag made patchy fog/drainage mosaics on the site-5 floor); faster layers are mixed drainage flow (mask 4, no fog): high sloping valleys drain instead of pooling"),
            row("fog_depth_max_m", k.fog_depth_max_m, "m", "cited", "radiation/valley fog is usually <= 200 m deep (design review mechanics 1.4 [V]); the fog top is capped at the pool floor (source bed) + this, deeper pools keep a clear upper inversion"),
            row("wind_gate_lo_m_s", k.wind_gate_lo_m_s, "m/s", R, "cold pools and radiation fog do not form above ~3 m/s 10 m wind (design review image critique 7.2 [V]); depth scaled by 1 - smoothstep(lo, hi, wind at formation)"),
            row("wind_gate_hi_m_s", k.wind_gate_hi_m_s, "m/s", R, "see wind_gate_lo_m_s"),
            row("descent_max_posts", k.descent_max_posts, "posts", P, "steepest-descent walk length from a post to the valley (fits the 64-post halo with the core)"),
            row("link_search_nodes", k.link_search_nodes, "nodes", P, "the descent end point is matched to the nearest channel link within this many lattice nodes"),
            row("aspect_factor_min", k.aspect_factor_min, "1", P, "burn-off factor clamp: fog top falls at burn x cos(incidence) / sin(sun elevation) (shaded slopes slower, sun-facing faster; design review 7.3 [R])"),
            row("aspect_factor_max", k.aspect_factor_max, "1", P, "see aspect_factor_min"),
            row("steam_dt_min_k", k.steam_dt_min_k, "K", R, "steam (evaporation) fog over water warmer than the air by more than ~6-10 K (recalled; design review 7.4 [V] qualitative)"),
            row("steam_dt_full_k", k.steam_dt_full_k, "K", P, "water-air difference at which the steam fog reaches its full depth"),
            row("steam_depth_max_m", k.steam_depth_max_m, "m", R, "steam fog is shallow, ~10-50 m (recalled)"),
            row("steam_water_lag_days", k.steam_water_lag_days, "day", P, "water temperature follows the seasonal air temperature (at sea level) with this lag"),
            row("steam_enclosed_frac", k.steam_enclosed_frac, "1", P, "a sea post is an inlet/estuary when at least this land fraction lies within steam_enclosed_r_posts (0.4 lit a strip along every straight coast, where the shore post sees ~0.5 land)"),
            row("steam_enclosed_r_posts", k.steam_enclosed_r_posts, "posts", P, "half-width of the enclosure box (16 posts = 512 m)"),
        ]
    }

    pub fn digest(&self) -> String {
        let bytes = serde_json::to_vec(&self.table()).expect("constants serialize");
        format!("blake3:{}", blake3::hash(&bytes).to_hex())
    }
}

/// The chart-wide drainage network (rivernet lod-8 lattice).
pub struct Network {
    pub n: usize,
    pub half: i64,
    pub spacing_m: f64,
    pub kind: Vec<u8>,
    pub channel: Vec<u8>,
    pub receiver: Vec<u32>,
    pub bed_m: Vec<f64>,
    pub pos_x_m: Vec<f64>,
    pub pos_y_m: Vec<f64>,
    pub chart_decl: [i64; 3],
    pub digest: String,
    pub t1a_digest: String,
    pub dir: PathBuf,
    pub t1a_dir: PathBuf,
}

fn json(p: &Path) -> Result<serde_json::Value, String> {
    let s = std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
    serde_json::from_str(&s).map_err(|e| format!("{}: {e}", p.display()))
}

fn read_bytes(p: &Path, n: usize, w: usize) -> Result<Vec<u8>, String> {
    let b = std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?;
    if b.len() != n * w { return Err(format!("{}: {} bytes, expected {}", p.display(), b.len(), n * w)); }
    Ok(b)
}

impl Network {
    pub fn load(dir: &Path, t1a_dir: &Path) -> Result<Network, String> {
        let v = json(&dir.join("rivernet.json"))?;
        let h = &v["header"];
        if h["schema"].as_str() != Some("vz.t1b-rivernet/v1") {
            return Err(format!("unsupported rivernet schema {:?}", h["schema"]));
        }
        let ch = &h["recipe"]["chart"];
        let decl = [
            ch["origin_lat_e7"].as_i64().ok_or("origin_lat_e7")?,
            ch["origin_lon_e7"].as_i64().ok_or("origin_lon_e7")?,
            ch["datum_radius_mm"].as_i64().ok_or("datum_radius_mm")?,
        ];
        let t = json(&t1a_dir.join("skeleton.json"))?;
        let kb = std::fs::read(dir.join("fields/kind.u8")).map_err(|e| e.to_string())?;
        let nn = kb.len();
        let n = (nn as f64).sqrt().round() as usize;
        if n * n != nn || n % 2 == 0 { return Err("rivernet lattice is not an odd square".into()); }
        let f32s = |name: &str| -> Result<Vec<f64>, String> {
            let b = read_bytes(&dir.join("fields").join(name), nn, 4)?;
            Ok(b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]) as f64).collect())
        };
        let u32s = |p: &Path| -> Result<Vec<u32>, String> {
            let b = read_bytes(p, nn, 4)?;
            Ok(b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect())
        };
        let channel = read_bytes(&dir.join("fields/channel.u8"), nn, 1)?;
        let mut receiver = u32s(&dir.join("fields/receiver.u32"))?;
        let trc = u32s(&t1a_dir.join("fields/receiver.u32"))?;
        let bed_m = f32s("bed_m.f32")?;
        let pos_x_m = f32s("pos_x_m.f32")?;
        let pos_y_m = f32s("pos_y_m.f32")?;
        // Seeds (kind 1) carry their downstream links in T1a.
        for i in 0..nn {
            if kb[i] == 1 && receiver[i] as usize >= nn {
                receiver[i] = trc[i];
            }
        }
        let half = (n as i64 - 1) / 2;
        let spacing_m = 512.0;
        for i in 0..nn {
            let gx = (i % n) as i64 - half;
            let gy = (i / n) as i64 - half;
            if (pos_x_m[i] - gx as f64 * spacing_m).abs() > spacing_m || (pos_y_m[i] - gy as f64 * spacing_m).abs() > spacing_m {
                return Err(format!("node {i}: position not near its 512 m lattice point"));
            }
        }
        Ok(Network {
            n,
            half,
            spacing_m,
            kind: kb,
            channel,
            receiver,
            bed_m,
            pos_x_m,
            pos_y_m,
            chart_decl: decl,
            digest: v["digest"].as_str().unwrap_or("unknown").to_string(),
            t1a_digest: t["digest"].as_str().unwrap_or("unknown").to_string(),
            dir: dir.to_path_buf(),
            t1a_dir: t1a_dir.to_path_buf(),
        })
    }

    pub fn is_land(&self, i: usize) -> bool { self.kind[i] != 2 }
    pub fn is_channel(&self, i: usize) -> bool { self.is_land(i) && (self.channel[i] == 1 || self.kind[i] == 1) }
    pub fn node(&self, gx: i64, gy: i64) -> Option<usize> {
        if gx.abs() > self.half || gy.abs() > self.half { return None; }
        Some(((gy + self.half) as usize) * self.n + (gx + self.half) as usize)
    }

    pub fn chart(&self) -> Result<Chart, String> {
        Chart::new(ChartDeclV1::new(self.chart_decl[0], self.chart_decl[1], self.chart_decl[2])).map_err(|e| format!("chart: {e:?}"))
    }
}

/// One T0 cell's night, plus what COLDPOOL-2 needs at the formation time.
#[derive(Debug, Clone, Serialize)]
pub struct CellNight {
    pub ns: NightSummary,
    pub wind_tf_m_s: f64,
    pub cloud_tf: f64,
    pub t_air_sea_level_c: f64,
    pub t_water_c: f64,
}

pub struct CellCache<'m, 'a> {
    m: &'m WeatherModel<'a>,
    k: ColdFlowConstantsV1,
    local_day: i64,
    t_h: f64,
    pub cells: BTreeMap<u32, CellNight>,
}

impl<'m, 'a> CellCache<'m, 'a> {
    pub fn new(m: &'m WeatherModel<'a>, k: &ColdFlowConstantsV1, local_day: i64, t_h: f64) -> Self {
        CellCache { m, k: k.clone(), local_day, t_h, cells: BTreeMap::new() }
    }

    pub fn get(&mut self, cell: u32) -> &CellNight {
        if !self.cells.contains_key(&cell) {
            let m = self.m;
            let ns = night_summary(m, cell, self.local_day, self.t_h);
            let (lat, lon) = (ns.lat_deg, ns.lon_deg);
            let tf = self.t_h.min(ns.sunrise_h);
            let s = m.state_at(lat, lon, tf, None);
            let c = m.climate.cell_sample(cell);
            let sc = solar_clock(tf, lat, lon, m.k.axial_tilt_deg);
            let peak = if lat >= 0.0 { m.k.season_peak_day_north } else { m.k.season_peak_day_north - 182.5 };
            let ph = 2.0 * std::f64::consts::PI * (sc.doy - peak - self.k.steam_water_lag_days) / 365.0;
            let t_sl_mean = c.temperature_mean_c + m.k.lapse_rate_k_per_m * c.t0_elevation_m;
            let t_water = (t_sl_mean + c.temperature_half_range_k * ph.cos()).max(0.0);
            let t_air_sl = ns.profile_t_td[0][0];
            self.cells.insert(cell, CellNight { ns, wind_tf_m_s: s.wind_speed_m_s, cloud_tf: s.cloud_cover, t_air_sea_level_c: t_air_sl, t_water_c: t_water });
        }
        &self.cells[&cell]
    }

    /// Bilinear blend of a per-cell quantity at a point.
    pub fn blend(&mut self, lat: f64, lon: f64, f: impl Fn(&CellNight) -> f64) -> f64 {
        let w = self.m.climate.bilinear(lat, lon);
        let mut s = 0.0;
        for (c, wt) in w {
            if wt == 0.0 { continue; }
            s += wt * f(self.get(c));
        }
        s
    }

    /// Dewpoint depression of near-surface air at elevation z, blended at a point.
    pub fn dd_at(&mut self, lat: f64, lon: f64, z: f64) -> f64 { self.blend(lat, lon, |c| c.ns.dd_at(z)) }
}

/// The chart-wide airshed state at one world hour.
pub struct Airshed {
    pub n: usize,
    pub lat_deg: Vec<f64>,
    pub lon_deg: Vec<f64>,
    pub q: Vec<f64>,
    pub flux_m3_s: Vec<f64>,
    pub slope: Vec<f64>,
    /// Steady (full-night) drainage depth after the wind gate.
    pub h_ss_m: Vec<f64>,
    /// Depth at the query time (night growth).
    pub h_m: Vec<f64>,
    pub lid_m: Vec<f64>,
    pub fog_top_m: Vec<f64>,
    pub source: Vec<u32>,
    pub burn_m: Vec<f64>,
    pub sunrise_h: Vec<f64>,
    pub outlet_critical: Vec<bool>,
    /// Own drainage speed (m/s) and whether the node holds stagnant pooled air.
    pub speed_m_s: Vec<f64>,
    pub pooled: Vec<bool>,
    pub land_nodes: usize,
    pub pooled_nodes: usize,
    pub fog_nodes: usize,
}

/// Solve an increasing g(h) = target on [0, cap] by bisection.
fn solve(target: f64, cap: f64, g: impl Fn(f64) -> f64) -> f64 {
    if target <= 0.0 { return 0.0; }
    if g(cap) <= target { return cap; }
    let (mut lo, mut hi) = (0.0f64, cap);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if g(mid) < target { lo = mid; } else { hi = mid; }
    }
    0.5 * (lo + hi)
}

pub fn airshed(m: &WeatherModel, cache: &mut CellCache, net: &Network, chart: &Chart, k: &ColdFlowConstantsV1, t_h: f64) -> Airshed {
    let kv = &m.k;
    let nn = net.n * net.n;
    let nan = f64::NAN;
    let mut a = Airshed {
        n: net.n,
        lat_deg: vec![nan; nn],
        lon_deg: vec![nan; nn],
        q: vec![0.0; nn],
        flux_m3_s: vec![0.0; nn],
        slope: vec![nan; nn],
        h_ss_m: vec![0.0; nn],
        h_m: vec![0.0; nn],
        lid_m: vec![nan; nn],
        fog_top_m: vec![nan; nn],
        source: vec![NONE; nn],
        burn_m: vec![0.0; nn],
        sunrise_h: vec![nan; nn],
        outlet_critical: vec![false; nn],
        speed_m_s: vec![0.0; nn],
        pooled: vec![false; nn],
        land_nodes: 0,
        pooled_nodes: 0,
        fog_nodes: 0,
    };
    let area = net.spacing_m * net.spacing_m;
    let mut ss = vec![nan; nn];
    let mut burn_rate = vec![0.0; nn];
    let mut wind = vec![0.0; nn];
    // 1. node location and night quantities
    for i in 0..nn {
        if !net.is_land(i) { continue; }
        a.land_nodes += 1;
        // A few edge nodes are jittered just past the chart square: their
        // planet location is taken at the clamped position (weather only).
        let lim = vz_tilekey::CHART_HALF_EXTENT_MM as f64;
        let (la, lo) = chart
            .inverse_f64((net.pos_x_m[i] * 1000.0).clamp(-lim, lim), (net.pos_y_m[i] * 1000.0).clamp(-lim, lim))
            .expect("node inside chart");
        let (la, lo) = (la / DEG, lo / DEG);
        a.lat_deg[i] = la;
        a.lon_deg[i] = lo;
        a.q[i] = cache.blend(la, lo, |c| c.ns.cooling_q);
        ss[i] = cache.blend(la, lo, |c| c.ns.sunset_h);
        a.sunrise_h[i] = cache.blend(la, lo, |c| c.ns.sunrise_h);
        burn_rate[i] = cache.blend(la, lo, |c| c.ns.burn_rate_m_h);
        wind[i] = cache.blend(la, lo, |c| c.wind_tf_m_s);
        a.flux_m3_s[i] = k.production_max_m_s * a.q[i] * area;
    }
    // 2. topological order (donors before receivers) and flux accumulation
    let land_recv = |i: usize| -> Option<usize> {
        let r = net.receiver[i] as usize;
        if r < nn && net.is_land(r) { Some(r) } else { None }
    };
    let mut ndon = vec![0u32; nn];
    for i in 0..nn {
        if net.is_land(i) { if let Some(r) = land_recv(i) { ndon[r] += 1; } }
    }
    let mut queue: VecDeque<usize> = (0..nn).filter(|&i| net.is_land(i) && ndon[i] == 0).collect();
    let mut order = Vec::with_capacity(a.land_nodes);
    while let Some(i) = queue.pop_front() {
        order.push(i);
        if let Some(r) = land_recv(i) {
            a.flux_m3_s[r] += a.flux_m3_s[i];
            ndon[r] -= 1;
            if ndon[r] == 0 { queue.push_back(r); }
        }
    }
    assert_eq!(order.len(), a.land_nodes, "the drainage network has a cycle");
    // 3. steady drainage depth per node
    let rs = k.width_stencil_nodes as i64;
    let w_node = net.spacing_m / (2 * rs + 1) as f64;
    for &i in &order {
        let r = net.receiver[i] as usize;
        let (bed_r, marine_out) = if r < nn {
            if net.is_land(r) { (net.bed_m[r], false) } else { (net.bed_m[r].max(0.0), true) }
        } else {
            (nan, false)
        };
        let s = if bed_r.is_finite() {
            let d = (net.pos_x_m[i] - net.pos_x_m[r]).hypot(net.pos_y_m[i] - net.pos_y_m[r]).max(1.0);
            ((net.bed_m[i] - bed_r) / d).max(k.slope_min)
        } else {
            k.slope_min
        };
        a.slope[i] = s;
        let gx = (i % net.n) as i64 - net.half;
        let gy = (i / net.n) as i64 - net.half;
        let mut offs: Vec<f64> = Vec::new();
        for dy in -rs..=rs {
            for dx in -rs..=rs {
                if dx == 0 && dy == 0 { continue; }
                if let Some(j) = net.node(gx + dx, gy + dy) {
                    if net.is_channel(j) { offs.push((net.bed_m[j] - net.bed_m[i]).abs()); }
                }
            }
        }
        let dt = (kv.pool_inversion_max_k * a.q[i]).max(k.dt_floor_k);
        let gp = k.gravity_m_s2 * dt / k.t_ref_k;
        let area_cs = |h: f64| k.floor_width_min_m * h + h * h / k.side_slope + w_node * offs.iter().map(|&o| (h - o).max(0.0)).sum::<f64>();
        let width = |h: f64| k.floor_width_min_m + 2.0 * h / k.side_slope + w_node * offs.iter().filter(|&&o| o < h).count() as f64;
        let f = a.flux_m3_s[i];
        let mut h = solve(f, k.depth_cap_m, |h| (gp * h * s / k.drag_cd).sqrt() * area_cs(h));
        if marine_out {
            let hc = solve(f, k.depth_cap_m, |h| (gp * area_cs(h).powi(3) / width(h)).sqrt());
            if hc < h {
                h = hc;
                a.outlet_critical[i] = true;
            }
        }
        let gate = 1.0 - smoothstep(k.wind_gate_lo_m_s, k.wind_gate_hi_m_s, wind[i]);
        a.h_ss_m[i] = h * gate;
        a.speed_m_s[i] = (gp * h * s / k.drag_cd).sqrt();
        let grow = if t_h >= a.sunrise_h[i] { 1.0 } else { ((t_h - ss[i]) / (a.sunrise_h[i] - ss[i])).clamp(0.0, 1.0).sqrt() };
        a.h_m[i] = a.h_ss_m[i] * grow;
    }
    // 4. lids from the outlets upstream (receivers first)
    let mut lid0 = vec![nan; nn];
    for &i in order.iter().rev() {
        let r = net.receiver[i] as usize;
        let lr = if r < nn { if net.is_land(r) { lid0[r] } else { 0.0 } } else { f64::NEG_INFINITY };
        let own = net.bed_m[i] + a.h_m[i];
        if own >= lr {
            lid0[i] = own;
            a.source[i] = i as u32;
        } else if r < nn && net.is_land(r) {
            lid0[i] = lr;
            a.source[i] = a.source[r];
        } else {
            // Bed below sea level next to the sea: the lid is the sea surface.
            lid0[i] = lr;
            a.source[i] = i as u32;
        }
    }
    // 5. pool properties at each source, then burn-off; copy to every node
    let mut fog0 = vec![nan; nn];
    let mut burn = vec![0.0; nn];
    for &i in &order {
        if a.source[i] != i as u32 { continue; }
        let l = lid0[i];
        let d = a.h_m[i];
        let grow = if t_h >= a.sunrise_h[i] { 1.0 } else { ((t_h - ss[i]) / (a.sunrise_h[i] - ss[i])).clamp(0.0, 1.0).sqrt() };
        let dt = kv.pool_inversion_max_k * a.q[i] * grow;
        let dd = cache.dd_at(a.lat_deg[i], a.lon_deg[i], l);
        let f = if d < kv.pool_depth_min_m { nan } else if dd <= 0.0 { l } else if dt > dd { l - d * dd / dt } else { nan };
        // (f64::min would turn a NaN into the cap; keep no-fog as NaN.)
        fog0[i] = if f.is_nan() { nan } else { f.min(net.bed_m[i] + k.fog_depth_max_m) };
        burn[i] = if t_h > a.sunrise_h[i] { burn_rate[i] * (t_h - a.sunrise_h[i]) } else { 0.0 };
    }
    for &i in &order {
        let s = a.source[i] as usize;
        a.burn_m[i] = burn[s];
        a.lid_m[i] = lid0[s] - burn[s];
        let stagnant = a.speed_m_s[s] <= k.stagnant_speed_m_s;
        a.pooled[i] = stagnant && a.lid_m[i] - net.bed_m[i] >= kv.pool_depth_min_m;
        a.fog_top_m[i] = if a.pooled[i] { fog0[s] - burn[s] } else { nan };
        if a.pooled[i] {
            a.pooled_nodes += 1;
            if a.fog_top_m[i].is_finite() && a.fog_top_m[i] > net.bed_m[i] { a.fog_nodes += 1; }
        }
    }
    a
}

pub struct ColdFlowProduct {
    pub n: usize,
    pub fog_top_m: Vec<f32>,
    pub pool_top_m: Vec<f32>,
    pub mask: Vec<u8>,
    pub pool_id: Vec<u32>,
    pub fog_core_fraction: f64,
    pub pool_core_fraction: f64,
    pub drainage_core_fraction: f64,
    pub steam_core_fraction: f64,
    pub unmatched_core_fraction: f64,
    pub fog_top_core_min_m: f64,
    pub fog_top_core_max_m: f64,
    pub pool_ids_core: usize,
}

/// Steepest descent (8-neighbour, drop / distance) from post i, at most
/// `max_steps`, stopping on channel (1, 5) or water (2, 4) posts.
fn descend(w: &Window, i: usize, max_steps: usize) -> usize {
    let n = w.n;
    let mut p = i;
    for _ in 0..max_steps {
        if matches!(w.fixed[p], 1 | 2 | 4 | 5) { break; }
        let (x, y) = ((p % n) as i64, (p / n) as i64);
        let z = w.height_m[p];
        let mut best = None;
        let mut best_s = 0.0;
        for (dx, dy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)] {
            let (qx, qy) = (x + dx, y + dy);
            if qx < 0 || qy < 0 || qx >= n as i64 || qy >= n as i64 { continue; }
            let q = qy as usize * n + qx as usize;
            let dist = if dx != 0 && dy != 0 { std::f64::consts::SQRT_2 } else { 1.0 };
            let s = (z - w.height_m[q]) / dist;
            if s > best_s { best_s = s; best = Some(q); }
        }
        match best { Some(q) => p = q, None => break }
    }
    p
}

/// Sun direction (east, north, up) at a place and hour.
fn sun_vec(t_h: f64, lat: f64, lon: f64, tilt: f64) -> [f64; 3] {
    let sc = solar_clock(t_h, lat, lon, tilt);
    let e = sc.elevation_deg * DEG;
    let hra = (sc.solar_hour - 12.0) * 15.0 * DEG;
    let phi = lat * DEG;
    let az = (-hra.sin()).atan2(sc.declination_rad.tan() * phi.cos() - phi.sin() * hra.cos());
    [e.cos() * az.sin(), e.cos() * az.cos(), e.sin()]
}

pub fn coldflow_window(m: &WeatherModel, cache: &mut CellCache, net: &Network, a: &Airshed, k: &ColdFlowConstantsV1, w: &Window, t_h: f64) -> ColdFlowProduct {
    let kv = &m.k;
    let n = w.n;
    let nn = net.n * net.n;
    let mut fog_top = vec![f32::NAN; n * n];
    let mut pool_top = vec![f32::NAN; n * n];
    let mut mask = vec![0u8; n * n];
    let mut pool_id = vec![NONE; n * n];
    let rl = k.link_search_nodes as i64;
    let max_steps = k.descent_max_posts as usize;
    // Water posts and the enclosure test (box sums of land posts).
    let water = |f: u8| f == 2 || f == 4;
    let re = k.steam_enclosed_r_posts as usize;
    let mut pre = vec![0u32; (n + 1) * (n + 1)];
    for y in 0..n {
        for x in 0..n {
            let land = if water(w.fixed[y * n + x]) { 0 } else { 1 };
            pre[(y + 1) * (n + 1) + x + 1] = pre[y * (n + 1) + x + 1] + pre[(y + 1) * (n + 1) + x] - pre[y * (n + 1) + x] + land;
        }
    }
    let (mut nc, mut nf, mut np, mut nd, mut nst, mut nun) = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut fmin, mut fmax) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut ids = std::collections::BTreeSet::new();
    for i in 0..n * n {
        let gx = (i % n) as i64 + w.gx0;
        let gy = (i / n) as i64 + w.gy0;
        let core = gx >= w.rx * 512 && gx <= w.rx * 512 + 512 && gy >= w.ry * 512 && gy <= w.ry * 512 + 512;
        if core { nc += 1; }
        let z = w.height_m[i];
        if water(w.fixed[i]) {
            mask[i] = 3;
            // Steam fog over enclosed water.
            let (x, y) = (i % n, i / n);
            let (x0, x1, y0, y1) = (x.saturating_sub(re), (x + re).min(n - 1), y.saturating_sub(re), (y + re).min(n - 1));
            let land = pre[(y1 + 1) * (n + 1) + x1 + 1] + pre[y0 * (n + 1) + x0] - pre[y0 * (n + 1) + x1 + 1] - pre[(y1 + 1) * (n + 1) + x0];
            let tot = ((x1 - x0 + 1) * (y1 - y0 + 1)) as f64;
            if w.fixed[i] == 4 || land as f64 / tot >= k.steam_enclosed_frac {
                let (la, lo) = w.post_latlon(i);
                let q = cache.blend(la, lo, |c| c.ns.cooling_q);
                let t_air = cache.blend(la, lo, |c| c.t_air_sea_level_c) - kv.pool_inversion_max_k * q;
                let t_water = cache.blend(la, lo, |c| c.t_water_c);
                let wind = cache.blend(la, lo, |c| c.wind_tf_m_s);
                let sr = cache.blend(la, lo, |c| c.ns.sunrise_h);
                let br = cache.blend(la, lo, |c| c.ns.burn_rate_m_h);
                let gate = 1.0 - smoothstep(k.wind_gate_lo_m_s, k.wind_gate_hi_m_s, wind);
                let depth = k.steam_depth_max_m * smoothstep(k.steam_dt_min_k, k.steam_dt_full_k, t_water - t_air) * gate
                    - if t_h > sr { br * (t_h - sr) } else { 0.0 };
                if depth > 0.0 {
                    mask[i] = 5;
                    fog_top[i] = (z.max(0.0) + depth) as f32;
                    if core { nst += 1; }
                }
            }
            continue;
        }
        // Land post: descend, match a channel link, interpolate the airshed.
        let e = descend(w, i, max_steps);
        let ex = ((e % n) as i64 + w.gx0) as f64 * 32.0;
        let ey = ((e / n) as i64 + w.gy0) as f64 * 32.0;
        let cgx = (ex / net.spacing_m).round() as i64;
        let cgy = (ey / net.spacing_m).round() as i64;
        let mut best: Option<(f64, usize, f64)> = None;
        for dy in -rl..=rl {
            for dx in -rl..=rl {
                let Some(j) = net.node(cgx + dx, cgy + dy) else { continue };
                if !net.is_channel(j) { continue; }
                let r = net.receiver[j] as usize;
                let (px, py) = (net.pos_x_m[j], net.pos_y_m[j]);
                let (d, tau) = if r < nn {
                    let (vx, vy) = (net.pos_x_m[r] - px, net.pos_y_m[r] - py);
                    let l2 = vx * vx + vy * vy;
                    let tau = if l2 > 0.0 { (((ex - px) * vx + (ey - py) * vy) / l2).clamp(0.0, 1.0) } else { 0.0 };
                    ((ex - px - tau * vx).hypot(ey - py - tau * vy), tau)
                } else {
                    ((ex - px).hypot(ey - py), 0.0)
                };
                if best.map_or(true, |(bd, bj, _)| d < bd || (d == bd && j < bj)) { best = Some((d, j, tau)); }
            }
        }
        let Some((_, j, tau)) = best else {
            if core { nun += 1; }
            continue;
        };
        let r = net.receiver[j] as usize;
        // Receiver end values: land node, the sea (lid 0, bed 0, no fog), or none.
        let (lr, br, fr, bmr, idr, pr) = if r < nn && net.is_land(r) {
            (a.lid_m[r], net.bed_m[r], a.fog_top_m[r], a.burn_m[r], a.source[r], a.pooled[r])
        } else if r < nn {
            (0.0, 0.0, f64::NAN, 0.0, NONE, false)
        } else {
            (a.lid_m[j], net.bed_m[j], a.fog_top_m[j], a.burn_m[j], a.source[j], a.pooled[j])
        };
        let pooled = if tau < 0.5 { a.pooled[j] } else { pr };
        let lid = (1.0 - tau) * a.lid_m[j] + tau * lr;
        let bed = (1.0 - tau) * net.bed_m[j] + tau * br;
        let burn = (1.0 - tau) * a.burn_m[j] + tau * bmr;
        let (fj, near_id) = (a.fog_top_m[j], if tau < 0.5 { a.source[j] } else { idr });
        let mut f = if fj.is_finite() && fr.is_finite() { (1.0 - tau) * fj + tau * fr } else if tau < 0.5 { fj } else { fr };
        if !(z < lid) { continue; }
        pool_top[i] = lid as f32;
        pool_id[i] = near_id;
        if core && near_id != NONE { ids.insert(near_id); }
        if !pooled || lid - bed < kv.pool_depth_min_m {
            mask[i] = 4;
            if core { nd += 1; }
            continue;
        }
        // Aspect-dependent burn-off: faster on sun-facing slopes, slower in shade.
        if burn > 0.0 && f.is_finite() {
            let (x, y) = (i % n, i / n);
            let zx = |xx: usize| w.height_m[y * n + xx];
            let zy = |yy: usize| w.height_m[yy * n + x];
            let (xa, xb) = (x.saturating_sub(1), (x + 1).min(n - 1));
            let (ya, yb) = (y.saturating_sub(1), (y + 1).min(n - 1));
            let gxs = (zx(xb) - zx(xa)) / (32.0 * (xb - xa) as f64);
            let gys = (zy(yb) - zy(ya)) / (32.0 * (yb - ya) as f64);
            let nrm = [-gxs, -gys, 1.0];
            let nl = (nrm[0] * nrm[0] + nrm[1] * nrm[1] + 1.0).sqrt();
            let (la, lo) = w.post_latlon(i);
            let s = sun_vec(t_h, la, lo, kv.axial_tilt_deg);
            if s[2] > 0.0 {
                let cosi = (nrm[0] * s[0] + nrm[1] * s[1] + nrm[2] * s[2]) / nl;
                let fa = (cosi / s[2].max(0.1)).clamp(k.aspect_factor_min, k.aspect_factor_max);
                f -= burn * (fa - 1.0);
            }
        }
        if f.is_finite() && z < f {
            mask[i] = 1;
            fog_top[i] = f as f32;
            if core {
                nf += 1;
                fmin = fmin.min(f);
                fmax = fmax.max(f);
            }
        } else {
            mask[i] = 2;
            if core { np += 1; }
        }
    }
    let c = nc.max(1) as f64;
    ColdFlowProduct {
        n,
        fog_top_m: fog_top,
        pool_top_m: pool_top,
        mask,
        pool_id,
        fog_core_fraction: nf as f64 / c,
        pool_core_fraction: (nf + np) as f64 / c,
        drainage_core_fraction: nd as f64 / c,
        steam_core_fraction: nst as f64 / c,
        unmatched_core_fraction: nun as f64 / c,
        fog_top_core_min_m: fmin,
        fog_top_core_max_m: fmax,
        pool_ids_core: ids.len(),
    }
}
