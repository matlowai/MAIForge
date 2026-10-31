//! T1b 0.8.0, NET-ONE (lane NO, 2026-10-13; design review terrain strategy Q1/Q5 lane A,
//! `unpublished-design-record`). LABEL: pre-signoff,
//! fidelity tier Preview, derived, never canon. SimTime label
//! `t1b-relaxation-not-history`.
//!
//! Why: per-window relaxations at 32 m (0.4.0-0.7.0) could not pass seams or
//! held-out F6 at any affordable halo, because a threshold landscape's height
//! depends on its whole catchment (5-15 km). Two stages instead:
//!
//! 1. `t1b.net64`: ONE relaxation per lod-8 tile (65.536 km) on the lod-5 post
//!    lattice (64 m, 1025 x 1025 core posts) plus a declared halo, with the T1a
//!    channels and T0 water (zero-contour coast) fixed; the source term is the
//!    T0 subgrid-closure rate times a PER-T0-CELL gain, blended with the plain
//!    (non-land-weighted) bilinear hat weights so it tapers to the coast and
//!    has no cell-edge steps; the gains are solved during the run so each
//!    cell's hat-weighted mean meets the T0 bilinear reference (F6); flow is
//!    confined per T0 cell; routing is MFD (Freeman 1991) below a declared
//!    channel-head area and D8 above; hillslopes are Roering creep with the
//!    exposed-lithology Sc and K.
//! 2. The 32 m window (lod-6 region + a small halo) is a LOCAL FILL over its
//!    owning tile: every 64 m channel down to the channel-head area is fixed
//!    (Chaikin-smoothed planform, analytic concave bed slope ~ A^-theta between
//!    break nodes), the T1a channels and water as before, the ring at the 64 m
//!    surface (bicubic), and only hillslopes relax (same law, same gains).
//!
//! Everything is a pure function of (seed, T0, T1a, recipe, key): a tile never
//! reads another tile, a window reads only its owning tile.

use crate::t0::{CellWater, SiteRecord, T0};
use crate::t1b::{closure_rate_m_per_my, T0Extras, ReliefLawsV1, CHANNEL, FREE, LAKE, MARINE, SIMTIME_LABEL, CLOSURE_LABEL};
use crate::t1b_lem::{Ledger, RoeringConfig, RoeringScratch, SeriesPoint, D8N};
use crate::t1b_litho::{ExposedLithology, LithologyV1};
use crate::t1b_view::SkeletonView;
use crate::*;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap, VecDeque};
use std::path::Path;
use vz_keyed_rng::Scope;
use vz_tilekey::chart::{Chart, ChartDeclV1, RAD2DEG};

pub const NET_GENERATOR_NAME: &str = "t1b-net64";
pub const NET_GENERATOR_VERSION: &str = "0.8.0";
pub const NET_RECIPE_SCHEMA: &str = "vz.t1b-net64-recipe/v1";
pub const NET_SCHEMA: &str = "vz.t1b-net64-tile/v1";
pub const NET_LAYER: &str = "t1b.net64";
pub const PURPOSE_NET_ROUGHNESS: &str = "t1b.net64.initial.roughness";
pub const NET_LABEL: &str = "pre-signoff; T0 substrate R2 not met; fidelity tier Preview; derived, never canon; SimTime t1b-relaxation-not-history";
/// lod-5 posts (64 m) on the TileKeyV2 ladder.
pub const NET_LOD: u8 = 5;
pub const NET_SPACING_MM: i64 = 64_000;
pub const NET_SPACING_M: f64 = 64.0;
/// lod-5 cells per lod-8 tile side (65.536 km / 64 m).
pub const TILE8_CELLS5: i64 = 1024;
/// lod-8 node to lod-5 index factor (512 m / 64 m).
pub const L8_TO_L5: i64 = 8;

/// Declared parameters of T1b 0.8.0 (both stages). Integers only. Values are project assumptions
/// pre-registered in `unpublished-design-record` (lane NO) until ruled.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct NetParamsV1 {
    /// 64 m stage halo around the lod-8 tile [4096 m, design review Q1].
    pub halo_mm: i64,
    pub dt_kyr: i64,
    /// Steps before the first gain update, steps between updates, most updates, steps after the last update.
    pub warm_steps: u32,
    pub gain_every_steps: u32,
    pub gain_updates_max: u32,
    pub final_steps: u32,
    pub stats_every_steps: u32,
    /// Initial per-cell gain on the T0 closure rate, permille [620: T3/T3b ruling].
    pub gain0_permille: i64,
    pub gain_min_permille: i64,
    pub gain_max_permille: i64,
    /// One update multiplies a gain by the relief ratio clamped to [1/c, c], c = this / 1000.
    pub gain_step_clamp_permille: i64,
    /// A cell is converged when |hat-weighted mean - hat-weighted T0 bilinear| <= this.
    pub gain_tolerance_mm: i64,
    /// A cell's gain is solved only if its hat mass over core land posts is >= this permille of its area.
    pub gain_min_hat_mass_permille: i64,
    /// Channel-head area: D8 routing and fixed channels at or above, MFD below [0.5 km2, design review Q1].
    pub channel_head_area_m2: i64,
    /// Freeman (1991) MFD exponent, ppm [1.1].
    pub mfd_exponent_ppm: i64,
    pub roughness_rms_fraction_ppm: i64,
    pub roughness_rms_cap_mm: i64,
    pub roughness_realization: u32,
    pub base_floor_above_sea_mm: i64,
    pub creep_ratio_max_ppm: i64,
    pub lithology: LithologyV1,
    /// 32 m fill: halo, step, steps, profile concavity theta (ppm), Chaikin passes on channel segments.
    pub fill_halo_mm: i64,
    pub fill_dt_kyr: i64,
    pub fill_steps: u32,
    pub profile_theta_ppm: i64,
    pub curve_passes: u32,
    /// Hillslope solver: 0 = t1b_lem::roering_diffuse (Newton, pre-registered); 1 = the same Newton with a warm start
    /// (`roering_newton_ws`; same backward-Euler root to the 1e-9 m stopping tolerance, fewer iterations).
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32")]
    pub creep_solver: u32,
    /// Lane NO PREREG REVISION 1 (0 = as first registered): 1 = chart-wide stage-0 gains (fixed per cell for every tile),
    /// 64 m hillslopes by implicit linear diffusion + exact threshold cap, estuary rule.
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32")]
    pub rev: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rev1: Option<NetRev1V1>,
    pub lattice_rule: String,
    pub fixed_rule: String,
    pub uplift_rule: String,
    pub gain_rule: String,
    pub base_rule: String,
    pub routing_rule: String,
    pub hillslope_rule: String,
    pub fill_rule: String,
    pub profile_rule: String,
}

impl Default for NetParamsV1 {
    fn default() -> Self {
        NetParamsV1 {
            halo_mm: 4_096_000,
            dt_kyr: 250,
            warm_steps: 80,
            gain_every_steps: 20,
            gain_updates_max: 6,
            final_steps: 40,
            stats_every_steps: 4,
            gain0_permille: 620,
            gain_min_permille: 50,
            gain_max_permille: 20_000,
            gain_step_clamp_permille: 2_000,
            gain_tolerance_mm: 50_000,
            gain_min_hat_mass_permille: 100,
            channel_head_area_m2: 500_000,
            mfd_exponent_ppm: 1_100_000,
            roughness_rms_fraction_ppm: 5_000,
            roughness_rms_cap_mm: 5_000,
            roughness_realization: 0,
            base_floor_above_sea_mm: 1_000,
            creep_ratio_max_ppm: 990_000,
            lithology: LithologyV1::standin_v1(),
            fill_halo_mm: 2_048_000,
            fill_dt_kyr: 250,
            fill_steps: 40,
            profile_theta_ppm: 500_000,
            curve_passes: 3,
            creep_solver: 0,
            rev: 0,
            rev1: None,
            lattice_rule: "lod8_tile_lod5_lattice_v1: one lattice per lod-8 tile (8, tx, ty): lod-5 global posts [1024 tx - h, 1024 tx + 1024 + h] x [1024 ty - h, 1024 ty + 1024 + h], h = halo / 64 m; closed outer edge (no flux, no water across)".into(),
            fixed_rule: "t1a_channels_zero_contour_water_v1: a post is marine (sea level) iff its T0 cell is land or marine and the T0 elevation bilinear between cell centres is <= sea level; lake cells at their stage; T1a channel==1 links (node -> receiver, D8 at 512 m) rasterized as 8 lod-5 steps on the remaining posts at the T1a trunk levels interpolated linearly (lower bed kept where two links write one post); all other posts free".into(),
            uplift_rule: "t0_closure_rate_cell_gain_hat_v1: U(p) = sum over the 4 surrounding T0 cell centres of w_c(p) g_c R_c, w_c the PLAIN bilinear weights (not land-renormalised), R_c = K R_sub,c / ln(sqrt(area_c) / x_c) on land cells and 0 on water cells (so U tapers to the coast and has no cell-edge step), g_c the per-cell gain".into(),
            gain_rule: "cell_gain_hat_mean_v1: g_c starts at gain0; after warm_steps and then every gain_every_steps (at most gain_updates_max times, stopping once every solved cell is within tolerance), for every land cell c whose hat mass W_c = sum of w_c(p) over core non-water posts is >= gain_min_hat_mass x area_c / post area: M = hat-weighted mean of z, T = hat-weighted mean of the T0 bilinear, B = hat-weighted mean of the initial base; g_c <- clamp(g_c x clamp((T - B) / max(M - B, 1 m), 1/c, c), gain_min, gain_max); other cells keep gain0; then final_steps with the gains frozen".into(),
            base_rule: "t1a_trunk_level_bicubic_v1: T1a trunk_level_m by Catmull-Rom bicubic between lod-8 nodes (clamped to the chart), floored at sea level + base floor on land; plus keyed roughness".into(),
            routing_rule: "priority_flood_fifo_mfd_v1: priority flood from all fixed posts, key (fill, index), pits filled through a FIFO (Barnes et al. 2014), confined per T0 cell (flow j -> i iff same T0 cell, or i fixed in the T0 receiver cell of j's cell; unreachable posts flooded unconfined in a second pass); area accumulated donors first: a post with area < channel-head area that is not under a fill sends its area to every allowed lower-fill neighbour popped before it with weight ~ (drop / step length)^p (Freeman 1991), else to its D8 steepest-descent receiver (flats: flood parent); stream power n = 1 implicit receivers first, z' = (z + sum f_j b_j) / (1 + sum f_j), f_j = K k_litho A^m w_j dt / L_j".into(),
            hillslope_rule: "nonlinear_creep_roering_v1 (t1b_lem::roering_diffuse, unchanged): q = D S / (1 - (S/Sc)^2), D = T0 header D, Sc = exposed lithology (cover_over_basement_v1 stand-in), factor tangent-extended beyond creep_ratio_max x Sc; x then y implicit sweeps".into(),
            fill_rule: "net_local_fill_v1: 32 m window = region core + fill halo, owned by the lod-8 tile containing its core; fixed: zero-contour water and T1a channel links (16 lod-4 steps) as 0.4.0, then every 64 m channel (free post with area >= channel-head area, D8 receiver) as Chaikin-smoothed segments between break nodes (heads, junctions with >= 2 channel donors, the fixed post a channel ends on) rasterized at 32 m (samples <= 8 m, nearest post, lowest bed kept) at the analytic profile bed (code 7), outer ring at the 64 m surface bicubic (code 3); free posts start at the 64 m surface (Catmull-Rom bicubic between lod-5 posts) and relax fill_steps with the same law, the tile's gains and no roughness".into(),
            profile_rule: "analytic_concave_profile_v1: along a 64 m channel segment n0..nk the bed falls from fill(n0) to fill(nk) in proportion to the cumulative integral of A^-theta ds (A linear between the 64 m nodes' areas; theta = m/n of the T0 law), i.e. slope ~ A^-theta between fixed nodes; then Chaikin with pinned ends (keeps beds non-rising)".into(),
        }
    }
}

/// Revision-1 constants (PREREG.txt REVISION 1).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct NetRev1V1 {
    pub gains_dt_kyr: i64,
    pub gains_warm_steps: u32,
    pub gains_settle_steps: u32,
    pub gains_updates_max: u32,
    pub gains_rule: String,
    pub hillslope_rule_64m: String,
    pub estuary_rule: String,
}

impl NetParamsV1 {
    /// PREREG REVISION 1 (lane NO, 2026-10-13).
    pub fn rev1() -> NetParamsV1 {
        NetParamsV1 {
            rev: 1, dt_kyr: 500, warm_steps: 80, gain_updates_max: 0, final_steps: 0, stats_every_steps: 4, creep_solver: 0,
            rev1: Some(NetRev1V1 {
                gains_dt_kyr: 500, gains_warm_steps: 60, gains_settle_steps: 40, gains_updates_max: 4,
                gains_rule: "chart_cell_gains_v1 (stage 0): one relaxation on the T1a lod-8 lattice (512 m, 513 x 513, closed chart edge, no roughness) with the 0.8.0 fixed rules, hat uplift, routing and the 64 m hillslope law; after warm steps and after each settle, every land cell whose hat mass over chart non-water nodes is >= gain_min_hat_mass of its area is updated g <- clamp(g x clamp((T - B) / max(M - B, 1 m), 1/c, c), bounds) (M, T, B hat-weighted means of z, T0 bilinear, initial base); stops when all solved cells are within tolerance or after the last update's settle; unsolved land cells take the hat-mass-weighted mean of the solved gains. The gains are a function of the T0 cell only: every lod-8 tile uses the same gains, fixed".into(),
                hillslope_rule_64m: "linear_implicit_threshold_cap_v1: implicit linear diffusion with the T0 header D (x then y backward Euler, one tridiagonal per run of free posts, fixed posts Dirichlet, closed lattice edge), then the exact min-plus threshold cap z_i <= z_j + Sc_i d_ij over 8 neighbours at the exposed-lithology Sc (t1b_lem::threshold_cap, T1b 0.6.0's law); the removed volume is exported (ledger failure_eroded)".into(),
                estuary_rule: "t1a_estuary_v1: a T1a channel post whose interpolated bed is <= sea level is marine (an arm of the sea at sea level)".into(),
            }),
            ..NetParamsV1::default()
        }
    }
    pub fn halo_posts(&self) -> i64 { self.halo_mm / NET_SPACING_MM }
    pub fn validate(&self) -> R<()> {
        if self.halo_mm % NET_SPACING_MM != 0 || self.halo_mm < 0 { return Err("net halo must be whole 64 m posts".into()); }
        if self.fill_halo_mm % crate::t1b::T1B_SPACING_MM != 0 || self.fill_halo_mm <= 0 || self.fill_halo_mm > self.halo_mm {
            return Err("fill halo must be whole 32 m posts, positive, and inside the 64 m halo".into());
        }
        if self.dt_kyr <= 0 || self.fill_dt_kyr <= 0 { return Err("dt must be positive".into()); }
        if self.gain_min_permille <= 0 || self.gain_max_permille < self.gain_min_permille || self.gain_step_clamp_permille < 1000 { return Err("gain bounds".into()); }
        if !(0..1_000_000).contains(&self.profile_theta_ppm) { return Err("theta in [0, 1)".into()); }
        if !(1..1_000_000).contains(&self.creep_ratio_max_ppm) { return Err("creep ratio".into()); }
        self.lithology.validate()
    }
    pub fn total_steps(&self, updates_done: u32) -> u32 { self.warm_steps + updates_done * self.gain_every_steps + self.final_steps }
}

/// One lod-8 tile's 64 m computation lattice (lod-5 global post indices).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct NetLat { pub tx: i64, pub ty: i64, pub gx0: i64, pub gy0: i64, pub n: usize, pub halo_posts: i64 }

impl NetLat {
    pub fn new(tx: i64, ty: i64, halo_posts: i64) -> NetLat {
        NetLat { tx, ty, gx0: tx * TILE8_CELLS5 - halo_posts, gy0: ty * TILE8_CELLS5 - halo_posts, n: (TILE8_CELLS5 + 2 * halo_posts + 1) as usize, halo_posts }
    }
    pub fn nodes(&self) -> usize { self.n * self.n }
    #[inline]
    pub fn idx(&self, gx: i64, gy: i64) -> Option<usize> {
        let (x, y) = (gx - self.gx0, gy - self.gy0);
        if x < 0 || y < 0 || x >= self.n as i64 || y >= self.n as i64 { return None; }
        Some(y as usize * self.n + x as usize)
    }
    #[inline]
    pub fn g(&self, i: usize) -> (i64, i64) { ((i % self.n) as i64 + self.gx0, (i / self.n) as i64 + self.gy0) }
    #[inline]
    pub fn in_core(&self, i: usize) -> bool {
        let (gx, gy) = self.g(i);
        let (x0, y0) = (self.tx * TILE8_CELLS5, self.ty * TILE8_CELLS5);
        gx >= x0 && gx <= x0 + TILE8_CELLS5 && gy >= y0 && gy <= y0 + TILE8_CELLS5
    }
    pub fn key(&self) -> String { format!("8.{}.{}", self.tx, self.ty) }
    /// The lod-8 tile owning lod-6 region (rx, ry).
    pub fn owner_of_region(rx: i64, ry: i64) -> (i64, i64) { (rx.div_euclid(4), ry.div_euclid(4)) }
}

#[derive(Clone, Debug, Serialize)]
pub struct NetRecipeV1 {
    pub schema: String,
    pub layer: String,
    pub generator: GeneratorRef,
    pub params: NetParamsV1,
    pub laws: ReliefLawsV1,
    pub universe_digest: String,
    pub upstream: Vec<(String, String)>,
    pub seed: u64,
    pub chart: ChartDeclV1,
    pub key_domain: String,
    pub simtime_label: String,
    pub closure_label: String,
    pub label: String,
}

pub fn net_recipe(t0: &T0, extras: &T0Extras, site: &SiteRecord, skel: &SkeletonView, params: NetParamsV1) -> R<NetRecipeV1> {
    params.validate()?;
    let mut upstream = vec![
        ("site".to_string(), site.bytes_digest.clone()),
        ("t0.hydrology".to_string(), t0.hydrology_digest.clone()),
        ("t0.subgrid_relief_sidecar".to_string(), t0.subgrid_digest.clone()),
        ("t0.terrain".to_string(), t0.terrain_digest.clone()),
        ("t1a.skeleton".to_string(), skel.digest.clone()),
    ];
    upstream.sort();
    Ok(NetRecipeV1 {
        schema: NET_RECIPE_SCHEMA.into(), layer: NET_LAYER.into(),
        generator: GeneratorRef { name: NET_GENERATOR_NAME.into(), version: NET_GENERATOR_VERSION.into(), tier: "Preview".into() },
        params, laws: extras.laws.clone(), universe_digest: site.universe_digest.clone(), upstream, seed: site.seed,
        chart: ChartDeclV1::new(site.chart_origin_lat_e7, site.chart_origin_lon_e7, (t0.radius_m * 1000.0).round() as i64),
        key_domain: format!("{NET_LAYER}/c0/8.*"),
        simtime_label: SIMTIME_LABEL.into(), closure_label: CLOSURE_LABEL.into(), label: NET_LABEL.into(),
    })
}

// ------------------------------------------------------------------ geometry helpers

/// The 4 T0 cells around (lat, lon) degrees with their plain bilinear weights (same indexing as `t1b::bilinear`).
pub fn hat4(t0: &T0, lat: f64, lon: f64) -> ([u32; 4], [f64; 4]) {
    let dlat = 180.0 / t0.nlat as f64;
    let dlon = 360.0 / t0.nlon as f64;
    let fi = (lat + 90.0) / dlat - 0.5;
    let fj = (lon + 180.0) / dlon - 0.5;
    let (i0, j0) = (fi.floor(), fj.floor());
    let (ti, tj) = (fi - i0, fj - j0);
    let mut c = [0u32; 4];
    let mut w = [0f64; 4];
    let mut k = 0;
    for (di, wi) in [(0i64, 1.0 - ti), (1, ti)] {
        let ii = (i0 as i64 + di).clamp(0, t0.nlat as i64 - 1) as usize;
        for (dj, wj) in [(0i64, 1.0 - tj), (1, tj)] {
            let jj = (j0 as i64 + dj).rem_euclid(t0.nlon as i64) as usize;
            c[k] = (ii * t0.nlon + jj) as u32;
            w[k] = wi * wj;
            k += 1;
        }
    }
    (c, w)
}

/// Catmull-Rom weights at t in [0, 1).
#[inline]
fn cr(t: f64) -> [f64; 4] {
    let (t2, t3) = (t * t, t * t * t);
    [0.5 * (-t3 + 2.0 * t2 - t), 0.5 * (3.0 * t3 - 5.0 * t2 + 2.0), 0.5 * (-3.0 * t3 + 4.0 * t2 + t), 0.5 * (t3 - t2)]
}

/// Catmull-Rom bicubic of node values `v(ix, iy)` at fractional node coordinates (x, y).
pub fn bicubic(v: &dyn Fn(i64, i64) -> f64, x: f64, y: f64) -> f64 {
    let (x0, y0) = (x.floor(), y.floor());
    let (wx, wy) = (cr(x - x0), cr(y - y0));
    let (x0, y0) = (x0 as i64, y0 as i64);
    let mut acc = 0.0;
    for (j, wyj) in wy.iter().enumerate() {
        if *wyj == 0.0 { continue; }
        let mut row = 0.0;
        for (i, wxi) in wx.iter().enumerate() { if *wxi != 0.0 { row += wxi * v(x0 - 1 + i as i64, y0 - 1 + j as i64); } }
        acc += wyj * row;
    }
    acc
}

/// The T1a trunk level, bicubic between lod-8 nodes, at lod-8 fractional coordinates.
pub fn trunk_bicubic(sk: &SkeletonView, x8: f64, y8: f64) -> f64 {
    bicubic(&|x, y| sk.trunk_level_m[node_index(x.clamp(G_MIN, G_MAX), y.clamp(G_MIN, G_MAX)).unwrap()] as f64, x8, y8)
}

/// Water class of a post under the zero-contour rule (T0 land/marine cells by the bilinear sign, lake cells by the cell).
pub fn water_zero_contour(t0: &T0, c: u32, bil: f64) -> CellWater {
    match t0.water(c) {
        CellWater::Lake => CellWater::Lake,
        _ => if bil > t0.sea_level_m { CellWater::Land } else { CellWater::Marine },
    }
}

// ------------------------------------------------------------------ routing (0.8.0)

#[derive(Clone, Copy, PartialEq)]
struct Key(f64, u32);
impl Eq for Key {}
impl PartialOrd for Key { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for Key { fn cmp(&self, o: &Self) -> Ordering { o.0.total_cmp(&self.0).then(o.1.cmp(&self.1)) } }

pub struct NetRouting {
    pub order: Vec<u32>,
    pub pos: Vec<u32>,
    pub receiver: Vec<u32>,
    pub dist: Vec<f64>,
    pub fill: Vec<f64>,
    pub phase2: Vec<bool>,
    pub unconfined: u32,
}

#[derive(Clone, Copy)]
pub struct Zones<'a> { pub zone: &'a [u32], pub drain: &'a [u32] }

#[inline]
fn allowed(z: Option<Zones>, fixed: &[u8], j: usize, i: usize) -> bool {
    match z { None => true, Some(zz) => zz.zone[i] == zz.zone[j] || (fixed[i] != FREE && zz.zone[i] == zz.drain[j]) }
}

/// Priority flood with a FIFO for pit posts (Barnes, Lehman and Mulla 2014), optional T0-cell confinement.
pub fn route(n: usize, z: &[f64], fixed: &[u8], dx: f64, zones: Option<Zones>) -> NetRouting {
    let nn = n * n;
    let mut fill = vec![0f64; nn];
    let mut parent = vec![NONE; nn];
    let mut seen = vec![false; nn];
    let mut heap = BinaryHeap::with_capacity(nn / 8);
    let mut pitq: VecDeque<u32> = VecDeque::new();
    for i in 0..nn { if fixed[i] != FREE { seen[i] = true; fill[i] = z[i]; heap.push(Key(z[i], i as u32)); } }
    let mut order = Vec::with_capacity(nn);
    let mut pos = vec![u32::MAX; nn];
    let mut phase2 = vec![false; nn];
    let mut unconfined = 0u32;
    for phase in 0..2 {
        if phase == 1 {
            if order.len() == nn || zones.is_none() { break; }
            for &i in &order { heap.push(Key(fill[i as usize], i)); }
        }
        loop {
            let (f, i) = if let Some(i) = pitq.pop_front() { (fill[i as usize], i as usize) } else if let Some(Key(f, i)) = heap.pop() { (f, i as usize) } else { break };
            if pos[i] == u32::MAX { pos[i] = order.len() as u32; order.push(i as u32); }
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            for (ox, oy) in D8N {
                let (xx, yy) = (x + ox, y + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let j = yy as usize * n + xx as usize;
                if seen[j] { continue; }
                if phase == 0 && !allowed(zones, fixed, j, i) { continue; }
                seen[j] = true;
                if phase == 1 { phase2[j] = true; unconfined += 1; }
                parent[j] = i as u32;
                if z[j] <= f { fill[j] = f; pitq.push_back(j as u32); } else { fill[j] = z[j]; heap.push(Key(z[j], j as u32)); }
            }
        }
    }
    let diag = dx * std::f64::consts::SQRT_2;
    let mut receiver = vec![NONE; nn];
    let mut dist = vec![0f64; nn];
    for i in 0..nn {
        if fixed[i] != FREE || pos[i] == u32::MAX { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let mut best: Option<(f64, usize, f64)> = None;
        for (ox, oy) in D8N {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if fill[j] < fill[i] && pos[j] < pos[i] && (phase2[i] || allowed(zones, fixed, i, j)) {
                let d = if ox != 0 && oy != 0 { diag } else { dx };
                let s = (fill[i] - fill[j]) / d;
                if best.map_or(true, |(bs, bj, _)| s > bs || (s == bs && j < bj)) { best = Some((s, j, d)); }
            }
        }
        match best {
            Some((_, j, d)) => { receiver[i] = j as u32; dist[i] = d; }
            None => {
                let p = parent[i] as usize;
                receiver[i] = p as u32;
                dist[i] = if (p % n) as i64 != x && (p / n) as i64 != y { diag } else { dx };
            }
        }
    }
    NetRouting { order, pos, receiver, dist, fill, phase2, unconfined }
}

/// Drainage area with MFD below `a_ch` (see routing_rule). Returns the area and, per post, the MFD weights over D8N
/// directions (all zero = D8 to `receiver`).
pub fn accumulate_mfd(n: usize, r: &NetRouting, z: &[f64], fixed: &[u8], dx: f64, a_ch: f64, p: f64, zones: Option<Zones>) -> (Vec<f64>, Vec<[f32; 8]>) {
    let nn = n * n;
    let a = dx * dx;
    let mut acc = vec![a; nn];
    let mut w8 = vec![[0f32; 8]; nn];
    let diag = dx * std::f64::consts::SQRT_2;
    for &i in r.order.iter().rev() {
        let i = i as usize;
        if fixed[i] != FREE { continue; }
        let wet = z[i] < r.fill[i];
        if !wet && acc[i] < a_ch {
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            let mut s = [0f64; 8];
            let mut tot = 0.0;
            for (k, (ox, oy)) in D8N.iter().enumerate() {
                let (xx, yy) = (x + ox, y + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let j = yy as usize * n + xx as usize;
                if r.fill[j] < r.fill[i] && r.pos[j] < r.pos[i] && (r.phase2[i] || allowed(zones, fixed, i, j)) {
                    let d = if *ox != 0 && *oy != 0 { diag } else { dx };
                    s[k] = ((r.fill[i] - r.fill[j]) / d).powf(p);
                    tot += s[k];
                }
            }
            if tot > 0.0 {
                for (k, (ox, oy)) in D8N.iter().enumerate() {
                    if s[k] == 0.0 { continue; }
                    let j = (y + oy) as usize * n + (x + ox) as usize;
                    let wk = s[k] / tot;
                    w8[i][k] = wk as f32;
                    acc[j] += wk * acc[i];
                }
                continue;
            }
        }
        let rc = r.receiver[i];
        if rc != NONE { acc[rc as usize] += acc[i]; }
    }
    (acc, w8)
}

/// D8 terminal (fixed post each post's D8 path ends on).
pub fn terminals(r: &NetRouting) -> Vec<u32> {
    let mut t = vec![NONE; r.receiver.len()];
    for &i in &r.order {
        let i = i as usize;
        t[i] = if r.receiver[i] == NONE { i as u32 } else { t[r.receiver[i] as usize] };
    }
    t
}

// ------------------------------------------------------------------ warm-started Roering creep (creep_solver 1)

/// Transport factor and derivative, identical to t1b_lem's private `creep_factor`.
#[inline]
fn creep_factor(s: f64, sc: f64, rm: f64) -> (f64, f64) {
    let r = s / sc;
    if r < rm { let q = 1.0 - r * r; (1.0 / q, 2.0 * r / sc / (q * q)) }
    else { let q = 1.0 - rm * rm; let (fm, dfm) = (1.0 / q, 2.0 * rm / sc / (q * q)); (fm + dfm * (s - rm * sc), dfm) }
}

/// Scratch for `roering_newton_ws` (same layout as t1b_lem::RoeringScratch, whose fields are private there).
pub struct WsScratch { crx: Vec<f64>, cry: Vec<f64>, scx: Vec<f64>, scy: Vec<f64>, a: Vec<f64>, b: Vec<f64>, c: Vec<f64>, d: Vec<f64>, idx: Vec<usize>, zl: Vec<f64>, z0: Vec<f64>, zb: Vec<f64> }
impl WsScratch {
    pub fn new(n: usize) -> Self {
        let v = |k: usize| vec![0.0; k];
        WsScratch { crx: v(n * n), cry: v(n * n), scx: v(n * n), scy: v(n * n), a: v(n), b: v(n), c: v(n), d: v(n), idx: Vec::with_capacity(n), zl: v(n + 2), z0: v(n + 2), zb: v(n + 2) }
    }
}

/// `t1b_lem::roering_diffuse` (same law, same LOD sweeps, same Newton with line search and the same 1e-9 m stopping rule)
/// with a WARM START: each line's Newton iteration starts from the start value plus that post's increment in the same sweep
/// of the previous step (`ws`), instead of from the start value. The root it converges to is the same backward-Euler
/// solution (to the stopping tolerance); only the iteration count changes. Counts total Newton iterations.
#[allow(clippy::too_many_arguments)]
pub fn roering_newton_ws(n: usize, z: &mut [f64], fixed: &[u8], sc: &[f64], d: f64, dt: f64, dx: f64, rc: &crate::t1b_lem::RoeringConfig, s: &mut WsScratch, area: f64, ws: &mut [Vec<f64>; 2], iters_total: &mut u64) -> crate::t1b_lem::RoeringStepLedger {
    let mut led = crate::t1b_lem::RoeringStepLedger::default();
    let rm = rc.ratio_max;
    let k_dt = dt / (dx * dx);
    for sweep in 0..2usize {
        // Lagged cross gradients and edge Sc for this sweep's edges (edge i = (i, i + stride)).
        let stride = if sweep == 0 { 1 } else { n };
        for y in 0..n { for x in 0..n {
            let i = y * n + x;
            let has = if sweep == 0 { x + 1 < n } else { y + 1 < n };
            if !has { continue; }
            let j = i + stride;
            let (fi, fj) = (fixed[i] == FREE, fixed[j] == FREE);
            if !fi && !fj { continue; }
            let up = if fi && fj { if z[i] >= z[j] { i } else { j } } else if fi { i } else { j };
            let cross_at = |k: usize| -> f64 {
                let (kx, ky) = (k % n, k / n);
                if sweep == 0 {
                    let (u, dd, h) = if ky + 1 < n && ky > 0 { (k + n, k - n, 2.0) } else if ky + 1 < n { (k + n, k, 1.0) } else { (k, k - n, 1.0) };
                    (z[u] - z[dd]) / (h * dx)
                } else {
                    let (u, dd, h) = if kx + 1 < n && kx > 0 { (k + 1, k - 1, 2.0) } else if kx + 1 < n { (k + 1, k, 1.0) } else { (k, k - 1, 1.0) };
                    (z[u] - z[dd]) / (h * dx)
                }
            };
            let cr = 0.5 * (cross_at(i) + cross_at(j));
            if sweep == 0 { s.crx[i] = cr; s.scx[i] = sc[up]; } else { s.cry[i] = cr; s.scy[i] = sc[up]; }
        } }
        for line in 0..n {
            let start = if sweep == 0 { line * n } else { line };
            let mut k = 0usize;
            while k < n {
                if fixed[start + k * stride] != FREE { k += 1; continue; }
                let mut k2 = k;
                s.idx.clear();
                while k2 < n && fixed[start + k2 * stride] == FREE { s.idx.push(start + k2 * stride); k2 += 1; }
                let m = s.idx.len();
                let left = (k > 0).then(|| start + (k - 1) * stride);
                let right = (k2 < n).then(|| start + k2 * stride);
                // Local line: zl[0] = left fixed (if any), zl[1..=m] the run, zl[m+1] = right fixed (if any).
                // Warm start: the start value plus this sweep's increment of the previous step (the converged root does not
                // depend on the start; near steady state the increments repeat, so Newton needs one or two iterations).
                for t in 0..m { s.zl[t + 1] = z[s.idx[t]] + ws[sweep][s.idx[t]]; s.z0[t + 1] = z[s.idx[t]]; }
                if let Some(j) = left { s.zl[0] = z[j]; }
                if let Some(j) = right { s.zl[m + 1] = z[j]; }
                // Edge e (0..=m) joins zl[e] and zl[e+1]; edge 0 exists iff left, edge m iff right. Its lagged data
                // lives at the lower-index post of the pair in lattice order.
                let edge_at = |s: &WsScratch, e: usize| -> Option<(f64, f64)> {
                    if (e == 0 && left.is_none()) || (e == m && right.is_none()) { return None; }
                    let lo = if e == 0 { left.unwrap() } else { s.idx[e - 1] };
                    Some(if sweep == 0 { (s.crx[lo], s.scx[lo]) } else { (s.cry[lo], s.scy[lo]) })
                };
                // flux on edge e toward zl[e] (mass from e+1 to e when positive) and dF/dG.
                let flux = |s: &WsScratch, e: usize| -> Option<(f64, f64, f64)> {
                    let (cr, scv) = edge_at(s, e)?;
                    let g = (s.zl[e + 1] - s.zl[e]) / dx;
                    let sm = (g * g + cr * cr).sqrt();
                    let (f, df) = creep_factor(sm, scv, rm);
                    let dfdg = if sm > 0.0 { f + df * g * g / sm } else { f };
                    Some((d * f * g, d * dfdg, sm / scv))
                };
                // max |R| on the current line state.
                let resid_max = |s: &WsScratch| -> f64 {
                    let mut fl = flux(s, 0).map_or(0.0, |v| v.0);
                    let mut mx = 0.0f64;
                    for t in 0..m {
                        let fr = flux(s, t + 1).map_or(0.0, |v| v.0);
                        mx = mx.max((s.zl[t + 1] - s.z0[t + 1] - dt / dx * (fr - fl)).abs());
                        fl = fr;
                    }
                    mx
                };
                let mut iters = 0u32;
                loop {
                    iters += 1;
                    let mut r0 = 0.0f64;
                    let mut fl = flux(s, 0);
                    for t in 0..m {
                        let fr = flux(s, t + 1);
                        let (fl_v, gl) = fl.map_or((0.0, 0.0), |v| (v.0, v.1));
                        let (fr_v, gr) = fr.map_or((0.0, 0.0), |v| (v.0, v.1));
                        // R = z - z0 - dt/dx (F_r - F_l)
                        s.d[t] = -(s.zl[t + 1] - s.z0[t + 1] - dt / dx * (fr_v - fl_v));
                        r0 = r0.max(s.d[t].abs());
                        s.b[t] = 1.0 + k_dt * (gl + gr);
                        s.a[t] = if t > 0 { -k_dt * gl } else { 0.0 };
                        s.c[t] = if t + 1 < m { -k_dt * gr } else { 0.0 };
                        fl = fr;
                    }
                    for t in 1..m { let w = s.a[t] / s.b[t - 1]; s.b[t] -= w * s.c[t - 1]; s.d[t] -= w * s.d[t - 1]; }
                    let mut nxt = 0.0;
                    let mut big = 0.0f64;
                    for t in (0..m).rev() {
                        let v = if t + 1 == m { s.d[t] / s.b[t] } else { (s.d[t] - s.c[t] * nxt) / s.b[t] };
                        nxt = v; s.d[t] = v; big = big.max(v.abs());
                    }
                    // Backtracking line search on max |R| (halve the step until the residual falls; at most 30 halvings).
                    for t in 0..m { s.zb[t + 1] = s.zl[t + 1]; }
                    let mut lam = 1.0f64;
                    loop {
                        for t in 0..m { s.zl[t + 1] = s.zb[t + 1] + lam * s.d[t]; }
                        if resid_max(s) < r0 || lam < 1e-9 { break; }
                        lam *= 0.5;
                    }
                    if big * lam <= 1e-9 || r0 <= 1e-10 { break; }
                    if iters >= 200 { led.lines_not_converged += 1; break; }
                }
                led.newton_iters_max = led.newton_iters_max.max(iters);
                *iters_total += iters as u64;
                for t in 0..m { ws[sweep][s.idx[t]] = s.zl[t + 1] - s.z0[t + 1]; }
                for t in 0..m {
                    let i = s.idx[t];
                    let dz = s.zl[t + 1] - z[i];
                    if dz < 0.0 { led.eroded -= dz * area; } else { led.deposited += dz * area; }
                    led.net_export -= dz * area;
                    z[i] = s.zl[t + 1];
                }
                // Diagnostics and export into the fixed ends (nonlinear flux on the final line state).
                for e in 0..=m {
                    if let Some((f, _, r)) = flux(s, e) {
                        if r >= 1.0 { led.over_sc_edges += 1; }
                        if r >= rm { led.clamped_edges += 1; }
                        // f > 0 moves mass from zl[e+1] to zl[e]: out of the run at the left end, into it at the right end.
                        if e == 0 { let j = left.unwrap(); let q = f * dt / dx * area; led.flux_export += q; led.export_by_class[fixed[j] as usize] += q; }
                        if e == m { let j = right.unwrap(); let q = -f * dt / dx * area; led.flux_export += q; led.export_by_class[fixed[j] as usize] += q; }
                    }
                }
                k = k2;
            }
        }
    }
    led
}


/// Implicit linear diffusion dz/dt = D lap z by x then y backward-Euler sweeps (one tridiagonal per run of free posts;
/// fixed posts Dirichlet with the flux into them exported; closed lattice edge). Unconditionally stable and monotone.
#[allow(clippy::too_many_arguments)]
pub fn linear_implicit(n: usize, z: &mut [f64], fixed: &[u8], d: f64, dt: f64, dx: f64, s: &mut WsScratch, area: f64) -> crate::t1b_lem::RoeringStepLedger {
    let mut led = crate::t1b_lem::RoeringStepLedger::default();
    let c = d * dt / (dx * dx);
    for sweep in 0..2usize {
        let stride = if sweep == 0 { 1 } else { n };
        for line in 0..n {
            let start = if sweep == 0 { line * n } else { line };
            let mut k = 0usize;
            while k < n {
                if fixed[start + k * stride] != FREE { k += 1; continue; }
                let mut k2 = k;
                s.idx.clear();
                while k2 < n && fixed[start + k2 * stride] == FREE { s.idx.push(start + k2 * stride); k2 += 1; }
                let m = s.idx.len();
                let left = (k > 0).then(|| start + (k - 1) * stride);
                let right = (k2 < n).then(|| start + k2 * stride);
                for t in 0..m {
                    let gl = if t > 0 || left.is_some() { c } else { 0.0 };
                    let gr = if t + 1 < m || right.is_some() { c } else { 0.0 };
                    s.b[t] = 1.0 + gl + gr;
                    s.a[t] = if t > 0 { -c } else { 0.0 };
                    s.c[t] = if t + 1 < m { -c } else { 0.0 };
                    let mut rhs = z[s.idx[t]];
                    if t == 0 { if let Some(j) = left { rhs += c * z[j]; } }
                    if t + 1 == m { if let Some(j) = right { rhs += c * z[j]; } }
                    s.d[t] = rhs;
                }
                for t in 1..m { let w = s.a[t] / s.b[t - 1]; s.b[t] -= w * s.c[t - 1]; s.d[t] -= w * s.d[t - 1]; }
                let mut nxt = 0.0;
                for t in (0..m).rev() { let v = if t + 1 == m { s.d[t] / s.b[t] } else { (s.d[t] - s.c[t] * nxt) / s.b[t] }; nxt = v; s.d[t] = v; }
                for t in 0..m {
                    let i = s.idx[t];
                    let dz = s.d[t] - z[i];
                    if dz < 0.0 { led.eroded -= dz * area; } else { led.deposited += dz * area; }
                    led.net_export -= dz * area;
                    z[i] = s.d[t];
                }
                if let Some(j) = left { let q = c * (z[s.idx[0]] - z[j]) * area; led.flux_export += q; led.export_by_class[fixed[j] as usize] += q; }
                if let Some(j) = right { let q = c * (z[s.idx[m - 1]] - z[j]) * area; led.flux_export += q; led.export_by_class[fixed[j] as usize] += q; }
                k = k2;
            }
        }
    }
    led
}

// ------------------------------------------------------------------ relaxation (both stages)

#[derive(Clone, Copy, Debug)]
pub struct RelaxCfg { pub n: usize, pub dx: f64, pub k_per_my: f64, pub m: f64, pub d_m2_per_my: f64, pub dt_my: f64, pub a_ch: f64, pub mfd_p: f64, pub ratio_max: f64, pub creep_solver: u32 }

pub struct Relax<'a> {
    pub cfg: RelaxCfg,
    pub fixed: &'a [u8],
    pub zones: Option<Zones<'a>>,
    pub litho: &'a dyn ExposedLithology,
    pub core: &'a [bool],
    pub z: Vec<f64>,
    pub ledger: Ledger,
    pub series: Vec<SeriesPoint>,
    sc: Vec<f64>,
    kmul: Vec<f64>,
    scratch: RoeringScratch,
    wscratch: WsScratch,
    ws: [Vec<f64>; 2],
    pub newton_iters_total: u64,
    z_prev: Vec<f64>,
    t_prev: f64,
    pub steps_done: u64,
    pub newton_iters_max: u32,
    pub lines_not_converged: u64,
    pub seconds_route: f64,
    pub seconds_fluvial: f64,
    pub seconds_creep: f64,
    diff_flux_export: f64,
}

impl<'a> Relax<'a> {
    pub fn new(cfg: RelaxCfg, fixed: &'a [u8], zones: Option<Zones<'a>>, litho: &'a dyn ExposedLithology, core: &'a [bool], z0: Vec<f64>) -> Self {
        let nn = cfg.n * cfg.n;
        let a = cfg.dx * cfg.dx;
        let mut ledger = Ledger::default();
        ledger.free_volume_start_m3 = (0..nn).filter(|&i| fixed[i] == FREE).map(|i| z0[i] * a).sum();
        let mut s = Relax { cfg, fixed, zones, litho, core, z_prev: z0.clone(), z: z0, ledger, series: Vec::new(), sc: vec![0.0; nn], kmul: vec![1.0; nn],
            scratch: RoeringScratch::new(if cfg.creep_solver == 1 { 1 } else { cfg.n }), wscratch: WsScratch::new(if cfg.creep_solver >= 1 { cfg.n } else { 1 }), ws: [vec![0.0; if cfg.creep_solver == 1 { nn } else { 0 }], vec![0.0; if cfg.creep_solver == 1 { nn } else { 0 }]], newton_iters_total: 0, t_prev: 0.0, steps_done: 0, newton_iters_max: 0, lines_not_converged: 0,
            seconds_route: 0.0, seconds_fluvial: 0.0, seconds_creep: 0.0, diff_flux_export: 0.0 };
        s.push_series(0);
        s
    }
    fn push_series(&mut self, pits: u32) {
        let t = self.steps_done as f64 * self.cfg.dt_my;
        let span = t - self.t_prev;
        let (mut s, mut mx, mut cnt, mut rs, mut rmx) = (0.0, f64::MIN, 0usize, 0.0, 0.0f64);
        for i in 0..self.z.len() {
            if !self.core[i] || self.fixed[i] != FREE { continue; }
            s += self.z[i]; mx = mx.max(self.z[i]); cnt += 1;
            if span > 0.0 { let r = ((self.z[i] - self.z_prev[i]) / span).abs(); rs += r; rmx = rmx.max(r); }
        }
        self.series.push(SeriesPoint { t_relax_my: t, core_free_mean_m: s / cnt.max(1) as f64, core_free_max_m: mx,
            core_mean_abs_rate_m_per_my: rs / cnt.max(1) as f64, core_max_abs_rate_m_per_my: rmx, core_pits: pits });
        self.z_prev.copy_from_slice(&self.z);
        self.t_prev = t;
    }
    /// Run `steps` steps with the source term `uplift` (m/My on free posts); a series point every `every` steps.
    pub fn run(&mut self, uplift: &[f64], steps: u32, every: u32) {
        let c = self.cfg;
        let n = c.n;
        let nn = n * n;
        let a = c.dx * c.dx;
        let dt = c.dt_my;
        let rc = RoeringConfig { ratio_max: c.ratio_max };
        for _ in 0..steps {
            for i in 0..nn {
                if self.fixed[i] == FREE { let u = uplift[i] * dt; self.z[i] += u; self.ledger.uplift_m3 += u * a; }
            }
            for i in 0..nn { if self.fixed[i] == FREE { let (s, k, _) = self.litho.at(i, self.z[i]); self.sc[i] = s; self.kmul[i] = k; } }
            let t = std::time::Instant::now();
            let r = route(n, &self.z, self.fixed, c.dx, self.zones);
            let (acc, w8) = accumulate_mfd(n, &r, &self.z, self.fixed, c.dx, c.a_ch, c.mfd_p, self.zones);
            let term = terminals(&r);
            self.seconds_route += t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            let wet: Vec<bool> = (0..nn).map(|i| self.z[i] < r.fill[i]).collect();
            let diag = c.dx * std::f64::consts::SQRT_2;
            for &i in &r.order {
                let i = i as usize;
                if self.fixed[i] != FREE { continue; }
                if wet[i] { self.ledger.submerged_post_steps += 1; continue; }
                let zi = self.z[i];
                let base = c.k_per_my * self.kmul[i] * acc[i].powf(c.m) * dt;
                let zn = if w8[i].iter().any(|w| *w > 0.0) {
                    let (x, y) = ((i % n) as i64, (i / n) as i64);
                    let (mut num, mut den) = (zi, 1.0);
                    for (k, (ox, oy)) in D8N.iter().enumerate() {
                        let wk = w8[i][k] as f64;
                        if wk == 0.0 { continue; }
                        let j = (y + oy) as usize * n + (x + ox) as usize;
                        let b = if wet[j] { r.fill[j] } else { self.z[j] };
                        if b >= zi { continue; }
                        let f = base * wk / if *ox != 0 && *oy != 0 { diag } else { c.dx };
                        num += f * b; den += f;
                    }
                    num / den
                } else {
                    let rcv = r.receiver[i] as usize;
                    let b = if wet[rcv] { r.fill[rcv] } else { self.z[rcv] };
                    if zi <= b { continue; }
                    let f = base / r.dist[i];
                    (zi + f * b) / (1.0 + f)
                };
                let v = (zi - zn) * a;
                self.ledger.stream_power_eroded_m3 += v;
                let cls = self.fixed[term[i] as usize] as usize;
                if cls < 5 { self.ledger.stream_power_export_by_class_m3[cls] += v; }
                self.z[i] = zn;
            }
            self.seconds_fluvial += t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            if c.creep_solver == 2 {
                let st = linear_implicit(n, &mut self.z, self.fixed, c.d_m2_per_my, dt, c.dx, &mut self.wscratch, a);
                self.ledger.diffusion_eroded_m3 += st.eroded;
                self.ledger.diffusion_deposited_m3 += st.deposited;
                self.ledger.diffusion_export_m3 += st.net_export;
                self.diff_flux_export += st.flux_export;
                for (k, v) in st.export_by_class.iter().enumerate() { if k < 5 { self.ledger.diffusion_export_by_class_m3[k] += v; } }
                let v = crate::t1b_lem::threshold_cap(n, &mut self.z, self.fixed, &self.sc, c.dx) * a;
                *self.ledger.failure_eroded_m3.get_or_insert(0.0) += v;
                self.seconds_creep += t.elapsed().as_secs_f64();
                self.steps_done += 1;
                if every > 0 && self.steps_done % every as u64 == 0 {
                    let pits = (0..nn).filter(|&i| self.core[i] && self.fixed[i] == FREE && self.z[i] < r.fill[i]).count() as u32;
                    self.push_series(pits);
                }
                continue;
            }
            let st = if c.creep_solver == 1 { roering_newton_ws(n, &mut self.z, self.fixed, &self.sc, c.d_m2_per_my, dt, c.dx, &rc, &mut self.wscratch, a, &mut self.ws, &mut self.newton_iters_total) }
                else { crate::t1b_lem::roering_diffuse(n, &mut self.z, self.fixed, &self.sc, c.d_m2_per_my, dt, c.dx, &rc, &mut self.scratch, a) };
            self.seconds_creep += t.elapsed().as_secs_f64();
            self.ledger.diffusion_eroded_m3 += st.eroded;
            self.ledger.diffusion_deposited_m3 += st.deposited;
            self.ledger.diffusion_export_m3 += st.net_export;
            self.diff_flux_export += st.flux_export;
            for (k, v) in st.export_by_class.iter().enumerate() { if k < 5 { self.ledger.diffusion_export_by_class_m3[k] += v; } }
            *self.ledger.creep_clamped_edge_steps.get_or_insert(0) += st.clamped_edges;
            *self.ledger.creep_over_sc_edge_steps.get_or_insert(0) += st.over_sc_edges;
            self.newton_iters_max = self.newton_iters_max.max(st.newton_iters_max);
            self.lines_not_converged += st.lines_not_converged as u64;
            self.steps_done += 1;
            if every > 0 && self.steps_done % every as u64 == 0 {
                let pits = (0..nn).filter(|&i| self.core[i] && self.fixed[i] == FREE && self.z[i] < r.fill[i]).count() as u32;
                self.push_series(pits);
            }
        }
    }
    /// Close the ledger (same definitions as `t1b_lem::run_full`).
    pub fn finish(&mut self) {
        let nn = self.z.len();
        let a = self.cfg.dx * self.cfg.dx;
        let l = &mut self.ledger;
        l.free_volume_end_m3 = (0..nn).filter(|&i| self.fixed[i] == FREE).map(|i| self.z[i] * a).sum();
        l.creep_newton_iters_max = Some(self.newton_iters_max);
        l.creep_lines_not_converged = Some(self.lines_not_converged);
        let fail = l.failure_eroded_m3.unwrap_or(0.0);
        l.eroded_m3 = l.stream_power_eroded_m3 + l.diffusion_eroded_m3 + fail;
        l.deposited_m3 = l.diffusion_deposited_m3;
        l.exported_m3 = l.stream_power_eroded_m3 + l.diffusion_export_m3 + fail;
        l.f7_relative_residual = if l.eroded_m3 > 0.0 { (l.eroded_m3 - l.deposited_m3 - l.exported_m3).abs() / l.eroded_m3 } else { 0.0 };
        let dv = l.free_volume_end_m3 - l.free_volume_start_m3;
        let denom = l.uplift_m3 + l.exported_m3.abs();
        l.volume_relative_residual = if denom > 0.0 { (dv - (l.uplift_m3 - l.exported_m3)).abs() / denom } else { 0.0 };
        l.diffusion_flux_relative_residual = if l.diffusion_export_m3.abs() > 0.0 { (self.diff_flux_export - l.diffusion_export_m3).abs() / l.diffusion_export_m3.abs() } else { (self.diff_flux_export - l.diffusion_export_m3).abs() };
        l.diffusion_substeps_per_step = 0;
        l.diffusion_number_per_substep = 4.0 * self.cfg.d_m2_per_my * self.cfg.dt_my / (a);
    }
}

// ------------------------------------------------------------------ per-cell gains

/// Per-post hat data and per-cell rates for the source term.
pub struct CellGains {
    pub cells: Vec<u32>,
    pub rate: Vec<f64>,
    pub gain: Vec<f64>,
    pub area_m2: Vec<f64>,
    pub post_cells: Vec<[u32; 4]>,
    pub post_w: Vec<[f64; 4]>,
}

impl CellGains {
    pub fn build(t0: &T0, laws: &ReliefLawsV1, latlon: &[(f64, f64)], gain0: f64, given: Option<&BTreeMap<u32, f64>>) -> CellGains {
        let mut map: BTreeMap<u32, u32> = BTreeMap::new();
        let mut post_cells = Vec::with_capacity(latlon.len());
        let mut post_w = Vec::with_capacity(latlon.len());
        for &(la, lo) in latlon {
            let (c, w) = hat4(t0, la, lo);
            post_cells.push(c);
            post_w.push(w);
            for cc in c { map.entry(cc).or_insert(0); }
        }
        let cells: Vec<u32> = map.keys().copied().collect();
        for (k, c) in cells.iter().enumerate() { map.insert(*c, k as u32); }
        for pc in post_cells.iter_mut() { for c in pc.iter_mut() { *c = map[c]; } }
        let rate: Vec<f64> = cells.iter().map(|&c| closure_rate_m_per_my(t0, laws, c)).collect();
        let gain: Vec<f64> = cells.iter().map(|c| given.and_then(|g| g.get(c).copied()).unwrap_or(gain0)).collect();
        let area_m2 = cells.iter().map(|&c| t0.cell_area_m2(c)).collect();
        CellGains { cells, rate, gain, area_m2, post_cells, post_w }
    }
    pub fn uplift(&self, fixed: &[u8]) -> Vec<f64> {
        (0..self.post_cells.len()).map(|i| if fixed[i] != FREE { 0.0 } else {
            (0..4).map(|q| { let k = self.post_cells[i][q] as usize; self.post_w[i][q] * self.gain[k] * self.rate[k] }).sum()
        }).collect()
    }
    pub fn gain_map(&self) -> BTreeMap<u32, f64> { self.cells.iter().zip(&self.gain).map(|(c, g)| (*c, *g)).collect() }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CellGainRow {
    pub cell: u32,
    pub land: bool,
    pub solved: bool,
    pub hat_mass_fraction: f64,
    pub gain: f64,
    pub hat_mean_z_m: f64,
    pub hat_mean_bilinear_m: f64,
    pub hat_mean_base_m: f64,
    pub diff_m: f64,
    pub at_bound: bool,
}

// ------------------------------------------------------------------ the 64 m tile

#[derive(Clone, Debug, Default, Serialize)]
pub struct NetStats {
    pub posts: u64,
    pub fixed_channel_posts: u64,
    pub marine_posts: u64,
    pub lake_posts: u64,
    pub free_posts: u64,
    pub channel_post_conflicts: u32,
    pub roughness_abs_max_m: f64,
    pub gain_updates_done: u32,
    pub gains_converged: bool,
    pub gain_history_max_abs_diff_m: Vec<f64>,
    pub cells: Vec<CellGainRow>,
    pub core_mean_m: f64,
    pub core_bilinear_mean_m: f64,
    pub net_channel_posts_core: u64,
    pub mfd_posts_last_step_core: u64,
    pub unconfined_posts_final: u32,
    pub f3_min_share_permille: f64,
    pub f3_cells_checked: u32,
    pub f3_cells_below_950: u32,
    pub pits_core_final: u32,
    pub creep_newton_iters_max: u32,
    pub creep_lines_not_converged: u64,
    pub steps: u64,
    pub sim_time_my: f64,
}

pub struct NetTile {
    pub recipe: NetRecipeV1,
    pub recipe_digest: String,
    pub lat: NetLat,
    pub height_mm: Vec<i32>,
    pub fill_mm: Vec<i32>,
    pub fixed: Vec<u8>,
    pub receiver: Vec<u32>,
    pub area_m2: Vec<f32>,
    pub terminal_node: Vec<u32>,
    pub uplift: Vec<f32>,
    pub stats: NetStats,
    pub header: Value,
    pub digest: String,
    pub blobs: Vec<(String, Vec<u8>)>,
    pub seconds: BTreeMap<String, f64>,
    /// Newton iterations summed over lines (creep_solver 1 only; not published).
    pub newton_iters_total: u64,
}

/// Stage 0 (REVISION 1): chart-wide per-cell gains.
pub struct GainLayer { pub gains: BTreeMap<u32, f64>, pub rows: Vec<CellGainRow>, pub header: Value, pub digest: String, pub seconds: f64 }

pub fn solve_gains(t0: &T0, extras: &T0Extras, site: &SiteRecord, skel: &SkeletonView, params: NetParamsV1) -> R<GainLayer> {
    let t_all = std::time::Instant::now();
    let recipe = net_recipe(t0, extras, site, skel, params)?;
    let recipe_digest = pf_contracts::artifact::digest_canonical(&recipe).map_err(|e| e.to_string())?;
    let p = &recipe.params;
    let r1 = p.rev1.as_ref().ok_or("solve_gains needs REVISION 1 params")?;
    let chart = Chart::new(recipe.chart.clone()).map_err(|e| e.to_string())?;
    let sea = t0.sea_level_m;
    let floor = sea + p.base_floor_above_sea_mm as f64 / 1000.0;
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let (mut fixed, mut z0, mut base, mut bil, mut latlon) = (vec![FREE; NODES], vec![0f64; NODES], vec![0f64; NODES], vec![0f64; NODES], Vec::with_capacity(NODES));
    let mut zone = vec![0u32; NODES];
    for i in 0..NODES {
        let (x, y) = node_mm(i);
        let (la, lo) = chart.inverse_mm(x, y).map_err(|e| e.to_string())?;
        let (la, lo) = (la * RAD2DEG, lo * RAD2DEG);
        latlon.push((la, lo));
        let c = t0.cell_of_deg(la, lo);
        zone[i] = c;
        bil[i] = crate::t1b::bilinear(t0, &elev, la, lo, false);
        base[i] = (skel.trunk_level_m[i] as f64).max(floor);
        match water_zero_contour(t0, c, bil[i]) {
            CellWater::Marine => { fixed[i] = MARINE; z0[i] = sea; base[i] = sea; continue; }
            CellWater::Lake => { fixed[i] = LAKE; z0[i] = t0.level(c); base[i] = z0[i]; continue; }
            CellWater::Land => {}
        }
        if (skel.kind[i] == 1 && skel.channel[i] == 1) || skel.kind[i] == 3 {
            let lev = skel.trunk_level_m[i] as f64;
            if lev <= sea { fixed[i] = MARINE; z0[i] = sea; } else { fixed[i] = CHANNEL; z0[i] = lev; }
            base[i] = z0[i];
            continue;
        }
        z0[i] = base[i];
    }
    let drain: Vec<u32> = zone.iter().map(|&c| t0.flow(c).unwrap_or(NONE)).collect();
    let pos = |i: usize| -> R<(f64, f64)> { Ok(latlon[i]) };
    let litho = crate::t1b_litho::build_standin(&p.lithology, t0, &extras.rock(), NODES, &pos)?;
    let mut cg = CellGains::build(t0, &extras.laws, &latlon, p.gain0_permille as f64 / 1000.0, None);
    let core = vec![true; NODES];
    let dx = SPACING_MM as f64 / 1000.0;
    let cfg = RelaxCfg { n: N, dx, k_per_my: extras.laws.stream_power_k_per_my, m: extras.laws.stream_power_m, d_m2_per_my: extras.laws.hillslope_diffusivity_m2_per_my,
        dt_my: r1.gains_dt_kyr as f64 / 1000.0, a_ch: p.channel_head_area_m2 as f64, mfd_p: p.mfd_exponent_ppm as f64 / 1e6, ratio_max: p.creep_ratio_max_ppm as f64 / 1e6, creep_solver: 2 };
    let mut rx = Relax::new(cfg, &fixed, Some(Zones { zone: &zone, drain: &drain }), &litho, &core, z0.clone());
    let mut up = cg.uplift(&fixed);
    rx.run(&up, r1.gains_warm_steps, 4);
    let land: Vec<usize> = (0..NODES).filter(|&i| fixed[i] != MARINE && fixed[i] != LAKE).collect();
    let rows_of = |cg: &CellGains, z: &[f64]| -> Vec<CellGainRow> {
        let k = cg.cells.len();
        let (mut w, mut mz, mut mb, mut mbase) = (vec![0.0; k], vec![0.0; k], vec![0.0; k], vec![0.0; k]);
        for &i in &land { for q in 0..4 { let c = cg.post_cells[i][q] as usize; let wq = cg.post_w[i][q]; w[c] += wq; mz[c] += wq * z[i]; mb[c] += wq * bil[i]; mbase[c] += wq * base[i]; } }
        (0..k).map(|c| {
            let landc = t0.water(cg.cells[c]) == CellWater::Land;
            let frac = w[c] * dx * dx / cg.area_m2[c];
            let (hz, hb, hbase) = if w[c] > 0.0 { (mz[c] / w[c], mb[c] / w[c], mbase[c] / w[c]) } else { (0.0, 0.0, 0.0) };
            CellGainRow { cell: cg.cells[c], land: landc, solved: landc && frac * 1000.0 >= p.gain_min_hat_mass_permille as f64, hat_mass_fraction: frac,
                gain: cg.gain[c], hat_mean_z_m: hz, hat_mean_bilinear_m: hb, hat_mean_base_m: hbase, diff_m: hz - hb,
                at_bound: (cg.gain[c] - p.gain_min_permille as f64 / 1000.0).abs() < 1e-12 || (cg.gain[c] - p.gain_max_permille as f64 / 1000.0).abs() < 1e-12 }
        }).collect()
    };
    let tol = p.gain_tolerance_mm as f64 / 1000.0;
    let cc = p.gain_step_clamp_permille as f64 / 1000.0;
    let (gmin, gmax) = (p.gain_min_permille as f64 / 1000.0, p.gain_max_permille as f64 / 1000.0);
    let mut history = Vec::new();
    let mut updates = 0u32;
    let mut converged = false;
    loop {
        let rows = rows_of(&cg, &rx.z);
        let worst = rows.iter().filter(|r| r.solved).map(|r| r.diff_m.abs()).fold(0.0, f64::max);
        history.push(worst);
        if worst <= tol { converged = true; break; }
        if updates >= r1.gains_updates_max { break; }
        for (k, r) in rows.iter().enumerate() {
            if !r.solved { continue; }
            let ratio = ((r.hat_mean_bilinear_m - r.hat_mean_base_m) / (r.hat_mean_z_m - r.hat_mean_base_m).max(1.0)).clamp(1.0 / cc, cc);
            cg.gain[k] = (cg.gain[k] * ratio).clamp(gmin, gmax);
        }
        updates += 1;
        up = cg.uplift(&fixed);
        rx.run(&up, r1.gains_settle_steps, 4);
    }
    rx.finish();
    let rows = rows_of(&cg, &rx.z);
    let (mut sw, mut sg) = (0.0, 0.0);
    for r in &rows { if r.solved { sw += r.hat_mass_fraction; sg += r.hat_mass_fraction * r.gain; } }
    let mean_solved = if sw > 0.0 { sg / sw } else { p.gain0_permille as f64 / 1000.0 };
    let gains: BTreeMap<u32, f64> = rows.iter().filter(|r| r.land).map(|r| (r.cell, if r.solved { r.gain } else { mean_solved })).collect();
    let hz: Vec<u8> = rx.z.iter().flat_map(|v| ((v * 1000.0).round() as i32).to_le_bytes()).collect();
    let header = serde_json::json!({
        "schema": "vz.t1b-net64-gains/v1", "layer": format!("{NET_LAYER}.gains"), "label": NET_LABEL, "simtime_label": SIMTIME_LABEL,
        "recipe_digest": recipe_digest, "rule": r1.gains_rule, "updates": updates, "converged": converged, "history_max_abs_diff_m": history,
        "mean_solved_gain": mean_solved, "gains": gains.iter().map(|(c, g)| serde_json::json!([c, g])).collect::<Vec<_>>(), "rows": rows,
        "ledger": rx.ledger, "series": rx.series, "height_512m_digest": format!("blake3:{}", blake3::hash(&hz).to_hex()),
    });
    let header = serde_json::from_str::<Value>(&header.to_string()).map_err(|e| e.to_string())?;
    let digest = pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())?;
    Ok(GainLayer { gains, rows, header, digest, seconds: t_all.elapsed().as_secs_f64() })
}

pub fn publish_gains(g: &GainLayer, out: &Path) -> R<()> {
    std::fs::create_dir_all(out).map_err(|e| e.to_string())?;
    std::fs::write(out.join("gains.json"), serde_json::to_vec_pretty(&serde_json::json!({"header": g.header, "digest": g.digest})).map_err(|e| e.to_string())?).map_err(|e| e.to_string())
}

/// Read stage-0 gains (digest-verified): (gains, digest, recipe digest).
pub fn load_gains(dir: &Path) -> R<(BTreeMap<u32, f64>, String, String)> {
    let v: Value = serde_json::from_slice(&std::fs::read(dir.join("gains.json")).map_err(|e| format!("gains.json: {e}"))?).map_err(|e| e.to_string())?;
    let h = v.get("header").ok_or("gains: header")?;
    let stated = v.get("digest").and_then(Value::as_str).ok_or("gains: digest")?;
    if pf_contracts::artifact::digest_canonical(h).map_err(|e| e.to_string())? != stated { return Err("gains digest mismatch; refusing".into()); }
    let g = h.get("gains").and_then(Value::as_array).ok_or("gains")?.iter().map(|v| Ok((v[0].as_u64().ok_or("cell")? as u32, v[1].as_f64().ok_or("g")?))).collect::<R<BTreeMap<u32, f64>>>()?;
    Ok((g, stated.into(), h.get("recipe_digest").and_then(Value::as_str).ok_or("recipe_digest")?.into()))
}

fn mm(v: f64) -> R<i32> { i32::try_from(vz_tilekey::chart::quantize_mm(v * 1000.0)).map_err(|_| format!("height {v} m does not fit i32 mm")) }

/// T1a channel links rasterized on a lattice with `steps` sub-steps per 512 m link. `put(gx, gy, level, n8)`.
fn t1a_links(sk: &SkeletonView, lo: (i64, i64), hi: (i64, i64), steps: i64, put: &mut dyn FnMut(i64, i64, f64, u32)) -> R<u32> {
    let mut links = 0;
    for gy8 in lo.1.div_euclid(steps) - 1..=hi.1.div_euclid(steps) + 1 { for gx8 in lo.0.div_euclid(steps) - 1..=hi.0.div_euclid(steps) + 1 {
        let Some(n8) = sk.idx(gx8, gy8) else { continue };
        if !(sk.kind[n8] == 1 && sk.channel[n8] == 1) { continue; }
        let r8 = sk.receiver[n8];
        if r8 == NONE || r8 == OFF_CHART { continue; }
        let (rgx, rgy) = node_g(r8 as usize);
        let (dx, dy) = (rgx - gx8, rgy - gy8);
        if dx.abs() > 1 || dy.abs() > 1 { return Err(format!("T1a link {n8} -> {r8} is not a D8 step")); }
        let (la, lb) = (sk.trunk_level_m[n8] as f64, sk.trunk_level_m[r8 as usize] as f64);
        for s in 0..=steps { put(gx8 * steps + dx * s, gy8 * steps + dy * s, la + (lb - la) * s as f64 / steps as f64, n8 as u32); }
        links += 1;
    } }
    Ok(links)
}

/// Build one lod-8 tile of the 64 m layer. Pure function of its inputs and the key.
pub fn build_tile(t0: &T0, extras: &T0Extras, site: &SiteRecord, skel: &SkeletonView, params: NetParamsV1, tx: i64, ty: i64) -> R<NetTile> {
    build_tile_g(t0, extras, site, skel, params, tx, ty, None)
}

/// `build_tile` with stage-0 gains (required iff params.rev == 1): (gains, gains digest, gains recipe digest).
#[allow(clippy::too_many_arguments)]
pub fn build_tile_g(t0: &T0, extras: &T0Extras, site: &SiteRecord, skel: &SkeletonView, params: NetParamsV1, tx: i64, ty: i64, gains: Option<&(BTreeMap<u32, f64>, String, String)>) -> R<NetTile> {
    let t_all = std::time::Instant::now();
    let mut seconds = BTreeMap::new();
    let recipe = net_recipe(t0, extras, site, skel, params)?;
    let recipe_digest = pf_contracts::artifact::digest_canonical(&recipe).map_err(|e| e.to_string())?;
    let p = &recipe.params;
    let lat = NetLat::new(tx, ty, p.halo_posts());
    let (n, nn) = (lat.n, lat.nodes());
    let half5 = vz_tilekey::CHART_HALF_EXTENT_MM / NET_SPACING_MM;
    if lat.gx0 < -half5 || lat.gy0 < -half5 || lat.gx0 + n as i64 - 1 > half5 || lat.gy0 + n as i64 - 1 > half5 {
        return Err(format!("tile {} with halo leaves the chart; refusing", lat.key()));
    }
    let chart = Chart::new(recipe.chart.clone()).map_err(|e| e.to_string())?;
    let sea = t0.sea_level_m;
    let floor = sea + p.base_floor_above_sea_mm as f64 / 1000.0;
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let mut latlon = Vec::with_capacity(nn);
    let mut t0_cell = vec![0u32; nn];
    let mut bil = vec![0f64; nn];
    let mut fixed = vec![FREE; nn];
    let mut base = vec![0f64; nn];
    let mut z0 = vec![0f64; nn];
    let mut st = NetStats { posts: nn as u64, ..Default::default() };
    let layer = Scope::root(site.seed).layer(&recipe_digest).map_err(|e| e.to_string())?.chart(0).map_err(|e| e.to_string())?;
    let mut tiles: BTreeMap<(i64, i64), Scope> = BTreeMap::new();
    for i in 0..nn {
        let (gx, gy) = lat.g(i);
        let (la, lo) = chart.inverse_mm(gx * NET_SPACING_MM, gy * NET_SPACING_MM).map_err(|e| e.to_string())?;
        let (la, lo) = (la * RAD2DEG, lo * RAD2DEG);
        latlon.push((la, lo));
        let c = t0.cell_of_deg(la, lo);
        t0_cell[i] = c;
        bil[i] = crate::t1b::bilinear(t0, &elev, la, lo, false);
        base[i] = trunk_bicubic(skel, gx as f64 / L8_TO_L5 as f64, gy as f64 / L8_TO_L5 as f64).max(floor);
        match water_zero_contour(t0, c, bil[i]) {
            CellWater::Marine => { fixed[i] = MARINE; z0[i] = sea; base[i] = sea; }
            CellWater::Lake => { fixed[i] = LAKE; z0[i] = t0.level(c); base[i] = z0[i]; }
            CellWater::Land => {
                let rms = ((t0.subgrid_relief_m[c as usize] as f64).max(0.0) * p.roughness_rms_fraction_ppm as f64 / 1e6).min(p.roughness_rms_cap_mm as f64 / 1000.0);
                let (ttx, tty) = (gx.div_euclid(128), gy.div_euclid(128));
                if !tiles.contains_key(&(ttx, tty)) { tiles.insert((ttx, tty), layer.tile(NET_LOD, ttx, tty).map_err(|e| e.to_string())?); }
                let u = tiles[&(ttx, tty)].draw(PURPOSE_NET_ROUGHNESS, ((gy - tty * 128) * 128 + (gx - ttx * 128)) as u64, 0).map_err(|e| e.to_string())?.uniform_f64();
                let rough = (2.0 * u - 1.0) * 3f64.sqrt() * rms;
                st.roughness_abs_max_m = st.roughness_abs_max_m.max(rough.abs());
                z0[i] = base[i] + rough;
            }
        }
    }
    // T1a channels.
    let mut chan: BTreeMap<usize, (f64, u32)> = BTreeMap::new();
    let mut conflicts = 0u32;
    t1a_links(skel, (lat.gx0, lat.gy0), (lat.gx0 + n as i64, lat.gy0 + n as i64), L8_TO_L5, &mut |gx, gy, lev, n8| {
        let Some(i) = lat.idx(gx, gy) else { return };
        if fixed[i] == MARINE || fixed[i] == LAKE { return; }
        match chan.get(&i).copied() {
            Some((old, _)) if old != lev => { conflicts += 1; if lev < old { chan.insert(i, (lev, n8)); } }
            Some(_) => {}
            None => { chan.insert(i, (lev, n8)); }
        }
    })?;
    st.channel_post_conflicts = conflicts;
    let mut chan_node = vec![NONE; nn];
    for (&i, &(lev, n8)) in &chan {
        if p.rev == 1 && lev <= sea { fixed[i] = MARINE; z0[i] = sea; base[i] = sea; continue; } // t1a_estuary_v1
        fixed[i] = CHANNEL; z0[i] = lev; base[i] = lev; chan_node[i] = n8;
    }
    for i in 0..nn { match fixed[i] { FREE => st.free_posts += 1, CHANNEL => st.fixed_channel_posts += 1, MARINE => st.marine_posts += 1, LAKE => st.lake_posts += 1, _ => {} } }
    // Lithology and gains.
    let pos = |i: usize| -> R<(f64, f64)> { Ok(latlon[i]) };
    let litho = crate::t1b_litho::build_standin(&p.lithology, t0, &extras.rock(), nn, &pos)?;
    if (p.rev == 1) != gains.is_some() { return Err("stage-0 gains are required iff rev == 1".into()); }
    if let Some(g) = gains {
        // Stage 0 reads neither the tile halo nor the roughness realization: compare with those normalised.
        let mut q = recipe.clone();
        q.params.halo_mm = NetParamsV1::rev1().halo_mm;
        q.params.roughness_realization = 0;
        let want = pf_contracts::artifact::digest_canonical(&q).map_err(|e| e.to_string())?;
        if g.2 != want { return Err("stage-0 gains were solved for a different recipe".into()); }
    }
    let mut cg = CellGains::build(t0, &extras.laws, &latlon, p.gain0_permille as f64 / 1000.0, gains.map(|g| &g.0));
    let core: Vec<bool> = (0..nn).map(|i| lat.in_core(i)).collect();
    let zone = t0_cell.clone();
    let drain: Vec<u32> = zone.iter().map(|&c| t0.flow(c).unwrap_or(NONE)).collect();
    seconds.insert("prepare".into(), t_all.elapsed().as_secs_f64());
    let cfg = RelaxCfg { n, dx: NET_SPACING_M, k_per_my: extras.laws.stream_power_k_per_my, m: extras.laws.stream_power_m,
        d_m2_per_my: extras.laws.hillslope_diffusivity_m2_per_my, dt_my: p.dt_kyr as f64 / 1000.0, a_ch: p.channel_head_area_m2 as f64,
        mfd_p: p.mfd_exponent_ppm as f64 / 1e6, ratio_max: p.creep_ratio_max_ppm as f64 / 1e6, creep_solver: if p.rev == 1 { 2 } else { p.creep_solver } };
    let t_relax = std::time::Instant::now();
    let mut rx = Relax::new(cfg, &fixed, Some(Zones { zone: &zone, drain: &drain }), &litho, &core, z0.clone());
    let mut up = cg.uplift(&fixed);
    rx.run(&up, p.warm_steps, p.stats_every_steps);
    let post_a = NET_SPACING_M * NET_SPACING_M;
    let land_core: Vec<usize> = (0..nn).filter(|&i| core[i] && fixed[i] != MARINE && fixed[i] != LAKE).collect();
    let cell_rows = |cg: &CellGains, z: &[f64]| -> Vec<CellGainRow> {
        let k = cg.cells.len();
        let (mut w, mut mz, mut mb, mut mbase) = (vec![0.0; k], vec![0.0; k], vec![0.0; k], vec![0.0; k]);
        for &i in &land_core { for q in 0..4 { let c = cg.post_cells[i][q] as usize; let wq = cg.post_w[i][q]; w[c] += wq; mz[c] += wq * z[i]; mb[c] += wq * bil[i]; mbase[c] += wq * base[i]; } }
        (0..k).map(|c| {
            let land = t0.water(cg.cells[c]) == CellWater::Land;
            let frac = w[c] * post_a / cg.area_m2[c];
            let (hz, hb, hbase) = if w[c] > 0.0 { (mz[c] / w[c], mb[c] / w[c], mbase[c] / w[c]) } else { (0.0, 0.0, 0.0) };
            CellGainRow { cell: cg.cells[c], land, solved: land && frac * 1000.0 >= p.gain_min_hat_mass_permille as f64, hat_mass_fraction: frac,
                gain: cg.gain[c], hat_mean_z_m: hz, hat_mean_bilinear_m: hb, hat_mean_base_m: hbase, diff_m: hz - hb,
                at_bound: (cg.gain[c] - p.gain_min_permille as f64 / 1000.0).abs() < 1e-12 || (cg.gain[c] - p.gain_max_permille as f64 / 1000.0).abs() < 1e-12 }
        }).collect()
    };
    let tol = p.gain_tolerance_mm as f64 / 1000.0;
    let clamp_c = p.gain_step_clamp_permille as f64 / 1000.0;
    let (gmin, gmax) = (p.gain_min_permille as f64 / 1000.0, p.gain_max_permille as f64 / 1000.0);
    let mut updates = 0u32;
    loop {
        let rows = cell_rows(&cg, &rx.z);
        let worst = rows.iter().filter(|r| r.solved).map(|r| r.diff_m.abs()).fold(0.0, f64::max);
        st.gain_history_max_abs_diff_m.push(worst);
        if worst <= tol { st.gains_converged = true; break; }
        if updates >= p.gain_updates_max { break; }
        for (k, r) in rows.iter().enumerate() {
            if !r.solved { continue; }
            let ratio = ((r.hat_mean_bilinear_m - r.hat_mean_base_m) / (r.hat_mean_z_m - r.hat_mean_base_m).max(1.0)).clamp(1.0 / clamp_c, clamp_c);
            cg.gain[k] = (cg.gain[k] * ratio).clamp(gmin, gmax);
        }
        updates += 1;
        up = cg.uplift(&fixed);
        rx.run(&up, p.gain_every_steps, p.stats_every_steps);
    }
    st.gain_updates_done = updates;
    rx.run(&up, p.final_steps, p.stats_every_steps);
    rx.finish();
    seconds.insert("relax".into(), t_relax.elapsed().as_secs_f64());
    seconds.insert("relax_route_accumulate".into(), rx.seconds_route);
    seconds.insert("relax_fluvial".into(), rx.seconds_fluvial);
    seconds.insert("relax_creep".into(), rx.seconds_creep);
    st.cells = cell_rows(&cg, &rx.z);
    let newton_iters_total = rx.newton_iters_total;
    st.steps = rx.steps_done;
    st.sim_time_my = rx.steps_done as f64 * cfg.dt_my;
    st.creep_newton_iters_max = rx.newton_iters_max;
    st.creep_lines_not_converged = rx.lines_not_converged;
    // Published surface, routing and network.
    let t = std::time::Instant::now();
    let height_mm: Vec<i32> = (0..nn).map(|i| mm(if fixed[i] == FREE { rx.z[i] } else { z0[i] })).collect::<R<_>>()?;
    let zq: Vec<f64> = height_mm.iter().map(|v| *v as f64 / 1000.0).collect();
    let zones = Some(Zones { zone: &zone, drain: &drain });
    let r = route(n, &zq, &fixed, NET_SPACING_M, zones);
    let (acc, w8) = accumulate_mfd(n, &r, &zq, &fixed, NET_SPACING_M, cfg.a_ch, cfg.mfd_p, zones);
    let term = terminals(&r);
    st.unconfined_posts_final = r.unconfined;
    let terminal_node: Vec<u32> = (0..nn).map(|i| {
        let tp = term[i] as usize;
        if tp == NONE as usize { return NONE; }
        if chan_node[tp] != NONE { return chan_node[tp]; }
        let (gx, gy) = lat.g(tp);
        skel.idx((gx as f64 / L8_TO_L5 as f64).round() as i64, (gy as f64 / L8_TO_L5 as f64).round() as i64).map_or(NONE, |v| v as u32)
    }).collect();
    let (mut sz, mut sb, mut cnt) = (0.0, 0.0, 0u64);
    for i in 0..nn {
        if !core[i] { continue; }
        sz += zq[i]; sb += bil[i]; cnt += 1;
        if fixed[i] == FREE {
            if acc[i] >= cfg.a_ch { st.net_channel_posts_core += 1; }
            if w8[i].iter().any(|w| *w > 0.0) { st.mfd_posts_last_step_core += 1; }
            if zq[i] < r.fill[i] { st.pits_core_final += 1; }
        }
    }
    st.core_mean_m = sz / cnt as f64;
    st.core_bilinear_mean_m = sb / cnt as f64;
    // F3 at 64 m (core free land posts): the D8 path ends in the post's own T0 cell or its receiver cell.
    let mut per: BTreeMap<u32, (u64, u64)> = BTreeMap::new();
    for i in 0..nn {
        if !core[i] || fixed[i] != FREE || term[i] == NONE { continue; }
        let c = t0_cell[i];
        let tn = terminal_node[i];
        let tc = if tn == NONE { NONE } else { skel.t0_cell[tn as usize] };
        let ok = tc == c || t0.flow(c) == Some(tc) || t0_cell[term[i] as usize] == c || t0.flow(c) == Some(t0_cell[term[i] as usize]);
        let e = per.entry(c).or_default();
        e.0 += 1; if ok { e.1 += 1; }
    }
    st.f3_min_share_permille = 1000.0;
    for (tot, ok) in per.values() {
        if *tot < 1000 { continue; }
        st.f3_cells_checked += 1;
        let s = 1000.0 * *ok as f64 / *tot as f64;
        if s < 950.0 { st.f3_cells_below_950 += 1; }
        st.f3_min_share_permille = st.f3_min_share_permille.min(s);
    }
    let fill_mm: Vec<i32> = r.fill.iter().map(|v| mm(*v)).collect::<R<_>>()?;
    let area_f: Vec<f32> = acc.iter().map(|v| *v as f32).collect();
    let up_f: Vec<f32> = up.iter().map(|v| *v as f32).collect();
    let le = |v: &[i32]| -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() };
    let leu = |v: &[u32]| -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() };
    let lef = |v: &[f32]| -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() };
    let bil_f: Vec<f32> = bil.iter().map(|v| *v as f32).collect();
    let base_mm: Vec<i32> = base.iter().map(|v| mm(*v)).collect::<R<_>>()?;
    let fields: Vec<(&str, &str, &str, Vec<u8>)> = vec![
        ("height_mm", "i32", "relaxed 64 m surface (t1b-relaxation-not-history), mm; fixed posts carry their fixed value", le(&height_mm)),
        ("fill_mm", "i32", "priority-flood fill of height_mm (all fixed posts outlets, T0-cell confined), mm: channel beds", le(&fill_mm)),
        ("fixed", "u8", "0 free, 1 T1a channel (fixed bed), 2 marine (zero contour), 4 lake (stage)", fixed.clone()),
        ("receiver", "u32", "D8 receiver (local lattice index) on the quantized surface; NONE on fixed posts", leu(&r.receiver)),
        ("area_m2", "f32", "drainage area, MFD below the channel-head area, D8 above (routing_rule)", lef(&area_f)),
        ("terminal_node", "u32", "T1a lod-8 node of the fixed post the D8 path ends on (link start node for T1a channels; nearest node for water)", leu(&terminal_node)),
        ("uplift_m_per_my", "f32", "final relaxation source term (t1b-relaxation-not-history; per-cell gains)", lef(&up_f)),
        ("t0_cell", "u32", "T0 L1 cell owning the post", leu(&t0_cell)),
        ("t0_bilinear_m", "f32", "T0 elevation bilinear between cell centres", lef(&bil_f)),
        ("base_mm", "i32", "initial base (T1a trunk level bicubic, floored) before roughness", le(&base_mm)),
    ];
    let mut entries = Vec::new();
    let mut blobs = Vec::new();
    for (name, dtype, note, bytes) in fields {
        let file = format!("fields/{name}.{dtype}");
        entries.push(serde_json::json!({"name": name, "dtype": dtype, "file": file, "note": note, "digest": format!("blake3:{}", blake3::hash(&bytes).to_hex())}));
        blobs.push((file, bytes));
    }
    let header = serde_json::json!({
        "schema": NET_SCHEMA, "layer": NET_LAYER, "label": NET_LABEL, "simtime_label": SIMTIME_LABEL, "tile": lat.key(),
        "lattice": lat, "lattice_order": "post = (gy - gy0) * n + (gx - gx0); lod-5 global posts; chart mm = g * 64000",
        "recipe": recipe, "recipe_digest": recipe_digest, "gains_stage0_digest": gains.map(|g| g.1.clone()), "gains": cg.gain_map().into_iter().map(|(c, g)| serde_json::json!([c, g])).collect::<Vec<_>>(),
        "stats": st, "ledger": rx.ledger, "series": rx.series, "fields": entries,
        "invented": ["the drainage network and relief between T1a channels by the relaxation process (keyed initial roughness at 64 m)", "per-T0-cell gains on the closure rate solved to the T0 bilinear (F6)"],
        "approximations": ["SimTime is t1b-relaxation-not-history", "closed lattice edge (halo absorbs it); tiles are independent, so lod-8 tile edges are seams", "T1a mouth links on reclaimed land keep the T1a sea-level bed (inland outlet; no mouth extension or deposition)", "fidelity tier Preview; pre-signoff substrate"],
    });
    let header = serde_json::from_str::<Value>(&header.to_string()).map_err(|e| e.to_string())?;
    let digest = pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())?;
    seconds.insert("publish_prep".into(), t.elapsed().as_secs_f64());
    seconds.insert("total".into(), t_all.elapsed().as_secs_f64());
    Ok(NetTile { recipe, recipe_digest, lat, height_mm, fill_mm, fixed, receiver: r.receiver, area_m2: area_f, terminal_node, uplift: up_f, stats: st, header, digest, blobs, seconds, newton_iters_total })
}

pub fn publish_tile(t: &NetTile, out: &Path) -> R<()> {
    std::fs::create_dir_all(out.join("fields")).map_err(|e| e.to_string())?;
    for (file, bytes) in &t.blobs { std::fs::write(out.join(file), bytes).map_err(|e| format!("{file}: {e}"))?; }
    let doc = serde_json::json!({"header": t.header, "digest": t.digest});
    std::fs::write(out.join("net.json"), serde_json::to_vec_pretty(&doc).map_err(|e| e.to_string())?).map_err(|e| e.to_string())
}

/// What the 32 m fill reads from a tile (digest-verified).
pub struct NetView {
    pub digest: String,
    pub recipe_digest: String,
    pub params: NetParamsV1,
    pub lat: NetLat,
    pub height_mm: Vec<i32>,
    pub fill_mm: Vec<i32>,
    pub fixed: Vec<u8>,
    pub receiver: Vec<u32>,
    pub area_m2: Vec<f32>,
    pub terminal_node: Vec<u32>,
    pub gains: BTreeMap<u32, f64>,
}

impl NetView {
    pub fn from_tile(t: &NetTile) -> NetView {
        NetView { digest: t.digest.clone(), recipe_digest: t.recipe_digest.clone(), params: t.recipe.params.clone(), lat: t.lat, height_mm: t.height_mm.clone(),
            fill_mm: t.fill_mm.clone(), fixed: t.fixed.clone(), receiver: t.receiver.clone(), area_m2: t.area_m2.clone(), terminal_node: t.terminal_node.clone(),
            gains: t.header["gains"].as_array().map(|a| a.iter().map(|v| (v[0].as_u64().unwrap() as u32, v[1].as_f64().unwrap())).collect()).unwrap_or_default() }
    }
    pub fn load(dir: &Path) -> R<NetView> {
        let bytes = std::fs::read(dir.join("net.json")).map_err(|e| format!("net.json in {}: {e}", dir.display()))?;
        let v: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
        let h = v.get("header").ok_or("net: no header")?;
        let stated = v.get("digest").and_then(Value::as_str).ok_or("net: no digest")?;
        let got = pf_contracts::artifact::digest_canonical(h).map_err(|e| e.to_string())?;
        if got != stated { return Err(format!("net header digest {got} differs from stated {stated}; refusing")); }
        let mut blobs = BTreeMap::new();
        for f in h.get("fields").and_then(Value::as_array).ok_or("net: fields")? {
            let name = f.get("name").and_then(Value::as_str).ok_or("name")?;
            let file = f.get("file").and_then(Value::as_str).ok_or("file")?;
            let b = std::fs::read(dir.join(file)).map_err(|e| format!("{file}: {e}"))?;
            if Some(format!("blake3:{}", blake3::hash(&b).to_hex()).as_str()) != f.get("digest").and_then(Value::as_str) { return Err(format!("net field {name} digest mismatch; refusing")); }
            blobs.insert(name.to_string(), b);
        }
        let b = |n: &str| blobs.get(n).cloned().ok_or_else(|| format!("net field {n} missing"));
        let i32s = |b: Vec<u8>| b.chunks_exact(4).map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect::<Vec<_>>();
        let u32s = |b: Vec<u8>| b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect::<Vec<_>>();
        let f32s = |b: Vec<u8>| b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect::<Vec<_>>();
        let l = h.get("lattice").ok_or("lattice")?;
        let gi = |k: &str| l.get(k).and_then(Value::as_i64).ok_or_else(|| format!("lattice.{k}"));
        let lat = NetLat { tx: gi("tx")?, ty: gi("ty")?, gx0: gi("gx0")?, gy0: gi("gy0")?, n: gi("n")? as usize, halo_posts: gi("halo_posts")? };
        let params: NetParamsV1 = serde_json::from_value(h.pointer("/recipe/params").ok_or("params")?.clone()).map_err(|e| e.to_string())?;
        let gains = h.get("gains").and_then(Value::as_array).ok_or("gains")?.iter()
            .map(|v| Ok((v[0].as_u64().ok_or("gain cell")? as u32, v[1].as_f64().ok_or("gain")?))).collect::<R<BTreeMap<u32, f64>>>()?;
        let nv = NetView { digest: stated.into(), recipe_digest: h.get("recipe_digest").and_then(Value::as_str).ok_or("recipe_digest")?.into(), params, lat,
            height_mm: i32s(b("height_mm")?), fill_mm: i32s(b("fill_mm")?), fixed: b("fixed")?, receiver: u32s(b("receiver")?), area_m2: f32s(b("area_m2")?),
            terminal_node: u32s(b("terminal_node")?), gains };
        if nv.height_mm.len() != lat.nodes() { return Err("net: field sizes".into()); }
        Ok(nv)
    }
    /// The 64 m surface by Catmull-Rom bicubic at fractional lod-5 global coordinates (clamped to the lattice).
    pub fn surface_m(&self, x5: f64, y5: f64) -> f64 {
        let n = self.lat.n as i64;
        bicubic(&|x, y| self.height_mm[((y - self.lat.gy0).clamp(0, n - 1) * n + (x - self.lat.gx0).clamp(0, n - 1)) as usize] as f64 / 1000.0, x5, y5)
    }
}

/// One smoothed 64 m channel segment in lod-4 units with its beds, the T1a node its path ends on, and its 64 m start post.
pub struct NetCurve { pub curve: crate::t1b_curves::GuideCurve, pub start: usize }

/// The tile's channel segments (free posts with area >= channel-head area), in start-post order.
pub fn net_curves(nv: &NetView) -> Vec<NetCurve> {
    let nn = nv.lat.nodes();
    let p = &nv.params;
    let a_ch = p.channel_head_area_m2 as f64;
    let theta = p.profile_theta_ppm as f64 / 1e6;
    let is_ch = |i: usize| nv.fixed[i] == FREE && nv.area_m2[i] as f64 >= a_ch && nv.receiver[i] != NONE;
    let mut donors = vec![0u32; nn];
    for i in 0..nn { if is_ch(i) { donors[nv.receiver[i] as usize] += 1; } }
    let n = nv.lat.n;
    let mut out = Vec::new();
    for s in 0..nn {
        if !is_ch(s) || donors[s] == 1 { continue; }
        let mut chain = vec![s];
        let mut cur = s;
        loop {
            let r = nv.receiver[cur] as usize;
            chain.push(r);
            if !is_ch(r) || donors[r] != 1 { break; }
            cur = r;
        }
        let k = chain.len() - 1;
        let (zs, ze) = (nv.fill_mm[s] as f64 / 1000.0, nv.fill_mm[chain[k]] as f64 / 1000.0);
        let area = |t: usize| -> f64 { let i = if t == k && !is_ch(chain[k]) { chain[k - 1] } else { chain[t] }; (nv.area_m2[i] as f64).max(1.0) };
        let mut cum = vec![0.0; k + 1];
        for t in 0..k {
            let (a, b) = (chain[t], chain[t + 1]);
            let len = if (a % n) != (b % n) && (a / n) != (b / n) { std::f64::consts::SQRT_2 } else { 1.0 };
            let am = 0.5 * (area(t) + area(t + 1));
            cum[t + 1] = cum[t] + len * am.powf(-theta);
        }
        let tot = cum[k].max(1e-300);
        let mut pts: Vec<[f64; 3]> = (0..=k).map(|t| {
            let (gx, gy) = nv.lat.g(chain[t]);
            [(2 * gx) as f64, (2 * gy) as f64, zs - (zs - ze) * cum[t] / tot]
        }).collect();
        for _ in 0..p.curve_passes { pts = crate::t1b_curves::chaikin_pinned(&pts); }
        out.push(NetCurve { curve: crate::t1b_curves::GuideCurve { pts, terminal: nv.terminal_node[s], end_node: chain[k] }, start: s });
    }
    out
}

// ------------------------------------------------------------------ the 32 m fill

#[derive(Clone, Debug, Default, Serialize)]
pub struct NetFillStats {
    pub net_tile: String,
    pub net_digest: String,
    pub net_segments_in_window: u32,
    pub net_posts: u32,
    pub net_post_conflicts: u32,
    pub net_posts_under_fixed: u32,
    pub ring_posts: u32,
    pub creep_newton_iters_max: u32,
    pub creep_lines_not_converged: u64,
    /// Timings: never serialized (the header is digested).
    #[serde(skip_serializing)]
    pub seconds_route: f64,
    #[serde(skip_serializing)]
    pub seconds_fluvial: f64,
    #[serde(skip_serializing)]
    pub seconds_creep: f64,
}

/// ReliefParamsV1 for a 0.8.0 window (the net block carries every 0.8.0 constant).
pub fn fill_params(np: &NetParamsV1) -> crate::t1b::ReliefParamsV1 {
    let mut p = crate::t1b::ReliefParamsV1::default().with_coast(1);
    p.halo_mm = np.fill_halo_mm;
    p.dt_kyr = np.fill_dt_kyr;
    p.sim_time_kyr = np.fill_dt_kyr * np.fill_steps as i64;
    p.stats_every_kyr = np.fill_dt_kyr * 4;
    p.roughness_rms_fraction_ppm = 0;
    p.uplift_gain_permille = np.gain0_permille;
    p.net = Some(np.clone());
    p
}

/// Compute one 32 m window from its owning tile. Reads nothing about other windows or tiles.
pub fn run_fill(inp: &crate::t1b::ReliefInputs, nv: &NetView, net_recipe_digest: &str, rx: i64, ry: i64) -> R<crate::t1b_window::WindowRun> {
    use crate::t1b::*;
    let t_all = std::time::Instant::now();
    let p = &inp.params;
    let np = p.net.as_ref().ok_or("run_fill needs params.net")?;
    if np != &nv.params { return Err("window net params differ from the tile's".into()); }
    if NetLat::owner_of_region(rx, ry) != (nv.lat.tx, nv.lat.ty) { return Err(format!("tile {} does not own region {rx},{ry}", nv.lat.key())); }
    if net_recipe_digest != nv.recipe_digest { return Err("net recipe digest mismatch".into()); }
    let mut recipe = recipe_for(inp)?;
    recipe.upstream.push((NET_LAYER.to_string(), net_recipe_digest.to_string()));
    recipe.upstream.sort();
    let recipe_digest = recipe.digest()?;
    let t0 = inp.t0;
    let lat = WinLat::new(rx, ry, p.halo_posts());
    let (n, nn) = (lat.n, lat.nodes());
    let chart = Chart::new(recipe.chart.clone()).map_err(|e| e.to_string())?;
    let sea = t0.sea_level_m;
    let floor = sea + np.base_floor_above_sea_mm as f64 / 1000.0;
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let mut stats = PrepStats::default();
    let mut fs = NetFillStats { net_tile: nv.lat.key(), net_digest: nv.digest.clone(), ..Default::default() };
    let (mut t0_cell, mut bil, mut base, mut fixed, mut z0) = (vec![0u32; nn], vec![0f64; nn], vec![0f64; nn], vec![FREE; nn], vec![0f64; nn]);
    let mut latlon = Vec::with_capacity(nn);
    let mut cells = std::collections::BTreeSet::new();
    for i in 0..nn {
        let (gx, gy) = lat.g(i);
        let (la, lo) = chart.inverse_mm(gx * T1B_SPACING_MM, gy * T1B_SPACING_MM).map_err(|e| e.to_string())?;
        let (la, lo) = (la * RAD2DEG, lo * RAD2DEG);
        latlon.push((la, lo));
        let c = t0.cell_of_deg(la, lo);
        cells.insert(c);
        t0_cell[i] = c;
        bil[i] = crate::t1b::bilinear(t0, &elev, la, lo, false);
        base[i] = nv.surface_m(gx as f64 / 2.0, gy as f64 / 2.0);
        match water_zero_contour(t0, c, bil[i]) {
            CellWater::Marine => { fixed[i] = MARINE; z0[i] = sea; }
            CellWater::Lake => { fixed[i] = LAKE; z0[i] = t0.level(c); }
            CellWater::Land => { base[i] = base[i].max(floor); z0[i] = base[i]; }
        }
    }
    stats.cells = cells.into_iter().collect();
    let mut chan_link_node = vec![NONE; nn];
    let mut chan: BTreeMap<usize, (f64, u32)> = BTreeMap::new();
    stats.channel_links = t1a_links(inp.skel, (lat.gx0, lat.gy0), (lat.gx0 + n as i64, lat.gy0 + n as i64), L8_TO_L4, &mut |gx, gy, lev, n8| {
        let Some(i) = lat.idx(gx, gy) else { return };
        if fixed[i] == MARINE || fixed[i] == LAKE { return; }
        match chan.get(&i).copied() {
            Some((old, _)) if old != lev => { stats.channel_post_conflicts += 1; stats.channel_post_conflict_max_m = stats.channel_post_conflict_max_m.max((old - lev).abs()); if lev < old { chan.insert(i, (lev, n8)); } }
            Some(_) => {}
            None => { chan.insert(i, (lev, n8)); }
        }
    })?;
    for (&i, &(lev, n8)) in &chan {
        if np.rev == 1 && lev <= sea { fixed[i] = MARINE; z0[i] = sea; continue; } // t1a_estuary_v1
        fixed[i] = CHANNEL; z0[i] = lev; chan_link_node[i] = n8;
    }
    // 64 m channels at 32 m.
    let curves = net_curves(nv);
    let (lx0, ly0, lx1, ly1) = ((lat.gx0 - 1) as f64, (lat.gy0 - 1) as f64, (lat.gx0 + n as i64) as f64, (lat.gy0 + n as i64) as f64);
    let mut net_level: BTreeMap<usize, (f64, u32)> = BTreeMap::new();
    let mut under = std::collections::BTreeSet::new();
    for cv in &curves {
        let c = &cv.curve;
        let inside = c.pts.iter().any(|q| q[0] >= lx0 && q[0] <= lx1 && q[1] >= ly0 && q[1] <= ly1);
        if !inside { continue; }
        fs.net_segments_in_window += 1;
        crate::t1b_curves::rasterize(c, &mut |gx, gy, lev| {
            let Some(i) = lat.idx(gx, gy) else { return };
            if fixed[i] != FREE { under.insert(i); return; }
            match net_level.get(&i).copied() {
                Some((old, _)) if lev >= old => {}
                Some(_) => { fs.net_post_conflicts += 1; net_level.insert(i, (lev, c.terminal)); }
                None => { net_level.insert(i, (lev, c.terminal)); }
            }
        });
    }
    fs.net_posts_under_fixed = under.len() as u32;
    for (&i, &(lev, term)) in &net_level { fixed[i] = NET; z0[i] = lev; chan_link_node[i] = term; }
    fs.net_posts = net_level.len() as u32;
    for i in 0..nn {
        if lat.on_ring(i) && fixed[i] == FREE { fixed[i] = RING; fs.ring_posts += 1; }
    }
    for i in 0..nn { match fixed[i] { FREE => stats.free_posts += 1, CHANNEL => stats.channel_posts += 1, MARINE => stats.marine_posts += 1, LAKE => stats.lake_posts += 1, _ => {} } }
    stats.ring_posts = fs.ring_posts;
    let pos = |i: usize| -> R<(f64, f64)> { Ok(latlon[i]) };
    let litho = crate::t1b_litho::build_standin(&np.lithology, t0, &inp.extras.rock(), nn, &pos)?;
    stats.lithology = Some(crate::t1b_litho::prep_stats(&litho, &|i| fixed[i] == FREE && lat.in_core(i)));
    let cg = CellGains::build(t0, &inp.extras.laws, &latlon, np.gain0_permille as f64 / 1000.0, Some(&nv.gains));
    let up = cg.uplift(&fixed);
    let core: Vec<bool> = (0..nn).map(|i| lat.in_core(i)).collect();
    let (mut us, mut uc) = (0.0, 0usize);
    for i in 0..nn { if core[i] && fixed[i] == FREE { us += up[i]; uc += 1; } }
    stats.uplift_core_mean_m_per_my = us / uc.max(1) as f64;
    let seconds_prepare = t_all.elapsed().as_secs_f64();
    let zone = t0_cell.clone();
    let drain: Vec<u32> = zone.iter().map(|&c| t0.flow(c).unwrap_or(NONE)).collect();
    let cfg = RelaxCfg { n, dx: T1B_SPACING_M, k_per_my: inp.extras.laws.stream_power_k_per_my, m: inp.extras.laws.stream_power_m,
        d_m2_per_my: inp.extras.laws.hillslope_diffusivity_m2_per_my, dt_my: np.fill_dt_kyr as f64 / 1000.0, a_ch: np.channel_head_area_m2 as f64,
        mfd_p: np.mfd_exponent_ppm as f64 / 1e6, ratio_max: np.creep_ratio_max_ppm as f64 / 1e6, creep_solver: np.creep_solver };
    let t = std::time::Instant::now();
    let mut rxl = Relax::new(cfg, &fixed, Some(Zones { zone: &zone, drain: &drain }), &litho, &core, z0.clone());
    rxl.run(&up, np.fill_steps, 4);
    rxl.finish();
    let seconds_relax = t.elapsed().as_secs_f64();
    fs.creep_newton_iters_max = rxl.newton_iters_max;
    fs.creep_lines_not_converged = rxl.lines_not_converged;
    fs.seconds_route = rxl.seconds_route;
    fs.seconds_fluvial = rxl.seconds_fluvial;
    fs.seconds_creep = rxl.seconds_creep;
    stats.net = Some(serde_json::to_value(&fs).map_err(|e| e.to_string())?);
    let height_mm: Vec<i32> = (0..nn).map(|i| mm(if fixed[i] == FREE { rxl.z[i] } else { z0[i] })).collect::<R<_>>()?;
    let zq: Vec<f64> = height_mm.iter().map(|v| *v as f64 / 1000.0).collect();
    let r = crate::t1b_lem::route(n, &zq, &fixed, T1B_SPACING_M);
    let acc = crate::t1b_lem::accumulate(&r, T1B_SPACING_M * T1B_SPACING_M, None);
    let ledger = rxl.ledger.clone();
    let series = rxl.series.clone();
    let prep = WindowInputs { lat, t0_cell: t0_cell.clone(), bilinear_m: bil.clone(), base_m: base.clone(), roughness_m: vec![0.0; nn], uplift_m_per_my: up.clone(),
        fixed: fixed.clone(), z0: z0.clone(), inflow_m2: vec![0.0; nn], chan_link_node, litho: None, stats };
    let mut raw = crate::t1b_window::RawWindow {
        lat, height_mm, fixed: fixed.clone(), drainage_area_m2: acc.iter().map(|v| *v as f32).collect(),
        t0_bilinear_m: bil.iter().map(|v| *v as f32).collect(), base_mm: base.iter().map(|v| mm(*v)).collect::<R<_>>()?,
        roughness_mm: vec![0; nn], initial_mm: z0.iter().map(|v| mm(*v)).collect::<R<_>>()?, uplift_m_per_my: up.iter().map(|v| *v as f32).collect(),
        t0_cell, inflow_area_m2: vec![0.0; nn],
        exposed_litho: Some((0..nn).map(|i| if fixed[i] == MARINE || fixed[i] == LAKE { 255 } else { litho.at(i, zq[i]).2 }).collect()),
        over: None, volcano_id: None, shore_mm: None, uncut_mm: None, ct_regime: None, ct_cliff_permille: None, ct_top_mm: None, ct_base_mm: None, ct_tread_permille: None, ct_tread_tan_permille: None, ct_ends: None, prov_id: None, runoff_weight: None, edifice_style: None, edifice_unit: None, front_mm: None, gl_carve_mm: None, digest: String::new(), recipe_digest: recipe_digest.clone(),
    };
    let (mut header, blobs) = crate::t1b_window::header_for(&raw, &recipe, &recipe_digest, &prep, ledger, series)?;
    header.preserved.push("the owning lod-8 tile's 64 m drainage network down to the channel-head area, as fixed posts (code 7) at the analytic profile beds: the same network in every window of the tile".into());
    header.invented.push("channel planform between 64 m nodes (Chaikin corner cutting, no draws) and the bed between break nodes (analytic concave profile, slope ~ A^-theta)".into());
    header.approximations.push("T1b 0.8.0 local fill: the window relaxes hillslopes only (channels >= channel-head area fixed from the 64 m tile); the ring is the 64 m surface; no roughness at 32 m".into());
    header.approximations.push("the window reads only its owning lod-8 tile: windows owned by different tiles meet at a lod-8 seam".into());
    raw.digest = pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())?;
    Ok(crate::t1b_window::WindowRun { raw, header, prep, blobs, seconds_prepare, seconds_relax })
}
