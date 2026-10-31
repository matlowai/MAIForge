//! Terrain tier T1b, regional relief (ADR 043 section 1, node TER-REGION;
//! lane TR, 2026-10-12). LABEL: pre-signoff, fidelity tier Preview, derived,
//! never canon.
//!
//! A T1b WINDOW is one region tile (TileKeyV2 lod 6, 16.384 km) plus a
//! declared halo, on the lod-4 post lattice (32 m). Each window is computed
//! lazily and independently (ADR 043 B6; ADR 044 section 6.2): its inputs are
//! only T0, the T1a skeleton and keyed draws addressed by GLOBAL lod-4 posts,
//! so a window's bytes never depend on which neighbours exist or ran first.
//!
//! The process is a relaxation run of the T0 run's own erosion law (stream
//! power n = 1, implicit receiver-first; linear hillslope diffusion, explicit
//! with substeps; constants read from the T0 terrain header), with:
//! - the T1a channels (channel == 1 links) FIXED at their T1a bed levels,
//! - marine T0 cells fixed at sea level, lake cells at their stage,
//! - the outer ring of the computation lattice fixed at the initial surface,
//! - a relaxation source term U = gain K R_sub / Lambda per T0 land cell: the
//!   rate T0's own subgrid closure (stream-power-subgrid-v1) assigns to the
//!   cell's subgrid relief, so the steady state under the same law carries
//!   the relief T0's closure assumes,
//! - keyed, bounded initial roughness (ADR 043 section 2, conditions 1-8,
//!   `sub-resolution-closure-v1`).
//! Its SimTime is labelled `t1b-relaxation-not-history` (ADR 043 B1): a
//! numerical device, not world time, not T0's history; nothing may cite it as
//! an age.

use crate::t0::{CellWater, SiteRecord, T0};
use crate::t1b_guides::GuideView;
use crate::t1b_view::SkeletonView;
use crate::*;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::Path;
use vz_keyed_rng::Scope;
use vz_tilekey::chart::{Chart, ChartDeclV1, RAD2DEG};

pub const T1B_GENERATOR_NAME: &str = "t1b-regional-relief";
pub const T1B_GENERATOR_VERSION: &str = "0.4.0";
/// Version with chart-wide tributary guides (guide_mode 1; lane TW, 2026-10-12).
pub const T1B_GENERATOR_VERSION_GUIDES: &str = "0.5.0";
/// Version with exposed lithology and threshold-slope failure (lithology_mode 1; lane T3, 2026-10-13).
pub const T1B_GENERATOR_VERSION_LITHO: &str = "0.6.0";
/// Version with the guide-surface ring, nonlinear creep, curved guides and the zero-contour coast (lane T3b, 2026-10-13).
pub const T1B_GENERATOR_VERSION_V4: &str = "0.7.0";
/// Version computed as a local 32 m fill over the one-lattice-per-lod-8-tile 64 m network (lane NO, T1B-NET-ONE, 2026-10-13).
pub const T1B_GENERATOR_VERSION_NET: &str = "0.8.0";
/// Version constructed river-network-first from the chart-wide `t1b.rivernet` layer (lane RV, T1B-RIVERNET, 2026-10-13).
pub const T1B_GENERATOR_VERSION_RIVERNET: &str = "0.9.0";
pub const T1B_RECIPE_SCHEMA: &str = "vz.t1b-relief-recipe/v1";
pub const T1B_LAYER: &str = "t1b.relief";
pub const SIMTIME_LABEL: &str = "t1b-relaxation-not-history";
pub const CLOSURE_LABEL: &str = "sub-resolution-closure-v1";
pub const PURPOSE_ROUGHNESS: &str = "t1b.initial.roughness";
pub const T1B_LABEL: &str = "pre-signoff; T0 substrate R2 not met; fidelity tier Preview; derived, never canon; SimTime t1b-relaxation-not-history";

pub const T1B_LOD: u8 = vz_tilekey::LOD_RELIEF;
pub const T1B_SPACING_MM: i64 = 32_000;
pub const T1B_SPACING_M: f64 = 32.0;
pub const REGION_LOD: u8 = vz_tilekey::LOD_REGION;
/// lod-4 cells per region side (16.384 km / 32 m).
pub const REGION_CELLS: i64 = 512;
/// lod-4 cells per lod-4 tile side.
pub const TILE_CELLS: i64 = 128;
/// lod-8 to lod-4 index factor (512 m / 32 m).
pub const L8_TO_L4: i64 = 16;

/// Fixed-node classes.
pub const FREE: u8 = 0;
pub const CHANNEL: u8 = 1;
pub const MARINE: u8 = 2;
/// Unused since 0.2.0 (the ring is free); kept so published codes stay stable.
pub const RING: u8 = 3;
pub const LAKE: u8 = 4;
/// Chart-wide tributary guide post (guide_mode 1): fixed at the guide bed of `t1b_guides`.
pub const GUIDE: u8 = 5;
/// Reclaimed coastal plain post (coast_mode 2): fixed at the declared plain surface (deposition is not modelled).
pub const PLAIN: u8 = 6;
/// T1b 0.8.0: a 64 m network channel (area >= channel-head area) fixed at its analytic profile bed (`t1b_net`).
pub const NET: u8 = 7;

pub(crate) fn is_zero_u32(v: &u32) -> bool { *v == 0 }
pub fn is_zero_u32_pub(v: &u32) -> bool { *v == 0 }

/// Declared parameters. Integers only (exact canonical JSON). [default] values
/// are project assumptions (ADR 039 C5) until ruled.
#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct ReliefParamsV1 {
    /// Halo around the region core [default 2048 m, ADR 043 section 4].
    pub halo_mm: i64,
    /// Seam blend half width on each side of a shared core edge [default 512 m, ADR 043 section 4].
    pub seam_half_width_mm: i64,
    /// Relaxation SimTime (t1b-relaxation-not-history) [default 40 My since 0.3.0: the
    /// core mean of window 6.0.0 drifts <= 4 m/My after 30 My from both initial bases
    /// (lane TR runs v3_*); 20 My (pre-registered) was still transient].
    pub sim_time_kyr: i64,
    /// Stream-power step [default 20 kyr; the implicit scheme is unconditionally stable,
    /// the n = 1 steady state does not depend on dt].
    pub dt_kyr: i64,
    /// Convergence series sampling [default 1 My].
    pub stats_every_kyr: i64,
    /// Calibration gain on the relaxation source term, permille [default 440 since 0.4.0:
    /// calibrated once on window 6.0.0 (seed 42) at SimTime 40 My so the core mean meets the
    /// T0 bilinear mean 2228 m (runs v4_g420 2128 m, v4_g460 2346 m; interpolated 438,
    /// rounded); pinned; tested on held-out windows. Pre-declared remedy in FALSIFIER.txt P5.
    /// Network-to-network scatter at a fixed gain is about +-100 m (v4_g380/420/460)].
    pub uplift_gain_permille: i64,
    /// Initial roughness rms as a fraction of the post's T0 cell subgrid relief, ppm
    /// [default 5000 = 0.5%] (ADR 043 section 2 condition 7, fraction).
    pub roughness_rms_fraction_ppm: i64,
    /// Initial roughness rms hard cap [default 5 m] (condition 7, cap).
    pub roughness_rms_cap_mm: i64,
    /// Roughness realization index (F12 uses 0 and 1).
    pub roughness_realization: u32,
    /// Base surface floor above sea level on land posts [default 1 m].
    pub base_floor_above_sea_mm: i64,
    /// Initial base: 0 = T0 bilinear (0.1.0, 0.2.0); 1 = T1a trunk level, bilinear between
    /// lod-8 nodes (relief grown under the source term from the channel-level surface) [default 1 since 0.3.0].
    pub initial_base_mode: u32,
    /// Ring inflow placement: 0 = ring posts of the 512 m strip (0.2.0); 1 = spread over the
    /// free halo-band posts of that strip (0.3.0 default), so no line is imposed.
    pub inflow_mode: u32,
    /// River threshold, milli m3/s (pf-surface RIVER_THRESHOLD_M3_S = 10; ADR 043 F1 default).
    pub river_threshold_milli_m3_s: i64,
    /// Valley-post area threshold for statistics [default 1 km2].
    pub valley_area_m2: i64,
    /// Local-relief square half width in posts [default 31 = 992 m].
    pub local_relief_half_width_posts: i64,
    pub uplift_rule: String,
    pub base_rule: String,
    pub roughness_rule: String,
    pub channel_rule: String,
    pub marine_rule: String,
    pub boundary_rule: String,
    pub routing_rule: String,
    pub area_rule: String,
    pub blend_rule: String,
    pub node_cell_rule: String,
    /// 0 = no guides (0.4.0, serialized exactly as before); 1 = chart-wide tributary guides fixed (0.5.0,
    /// ADR 043 amendment, design 2026-10-12).
    #[serde(skip_serializing_if = "is_zero_u32")]
    pub guide_mode: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub guide_rule: Option<String>,
    /// 0 = off (0.4.0 / 0.5.0, serialized exactly as before); 1 = exposed lithology with per-post K and
    /// threshold-slope failure (0.6.0, lane T3; design 2026-10-13: fewer bigger valleys + critical slope,
    /// rock-dependent).
    #[serde(skip_serializing_if = "is_zero_u32")]
    pub lithology_mode: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lithology: Option<crate::t1b_litho::LithologyV1>,
    /// Hillslope diffusivity as ppm of the T0 header D (absent = T0's D, as in 0.4.0-0.6.0 defaults).
    /// Lane T3 exploratory (not pre-registered): 50000 ppm = 0.005 m2/yr, a field-scale soil-creep value
    /// (recalled, unverified: Roering et al. 1999 WRR 35:853 ~0.003 m2/yr; Fernandes and Dietrich 1997 WRR 33:1307
    /// ~0.004-0.01 m2/yr); T0's 0.1 m2/yr is a coarse-grid legacy assumption (A-HILLSLOPE-D).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hillslope_d_ppm_of_t0: Option<i64>,
    /// T1b 0.7.0 switches (lane T3b, TER-REGION v4; all 0 = absent from the JSON, so 0.4.0-0.6.0 recipes are unchanged).
    /// 1 = the outer ring of the computation lattice is fixed at the guide layer's relaxed 512 m surface.
    #[serde(skip_serializing_if = "is_zero_u32")]
    pub ring_mode: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ring_rule: Option<String>,
    /// 1 = nonlinear (Roering-type) creep replaces linear diffusion and the threshold cap.
    #[serde(skip_serializing_if = "is_zero_u32")]
    pub hillslope_mode: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hillslope_rule: Option<String>,
    /// |grad z| / Sc beyond which the creep factor continues along its tangent, ppm (hillslope_mode 1).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub creep_ratio_max_ppm: Option<i64>,
    /// 1 = guide segments Chaikin-smoothed (`t1b_curves`), with this many passes.
    #[serde(skip_serializing_if = "is_zero_u32")]
    pub guide_curve_passes: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub guide_curve_rule: Option<String>,
    /// 1 = marine posts from the T0 bilinear zero contour instead of the T0 cell edge.
    #[serde(skip_serializing_if = "is_zero_u32")]
    pub coast_mode: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coast_rule: Option<String>,
    /// coast_mode 2: reclaimed-plain gradient as a fraction of the T0 bilinear height above sea, ppm.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coast_plain_ratio_ppm: Option<i64>,
    /// T1b 0.8.0 (lane NO): every constant of the 64 m network stage and the 32 m fill (absent = 0.4.0-0.7.0, serialized as before).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub net: Option<crate::t1b_net::NetParamsV1>,
    /// T1b 0.9.0 (lane RV): every constant of the river-network-first construction (absent = 0.4.0-0.8.0, serialized as before).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rivernet: Option<crate::t1b_rivernet::RiverNetParamsV1>,
}

impl Default for ReliefParamsV1 {
    fn default() -> Self {
        ReliefParamsV1 {
            halo_mm: 2_048_000,
            seam_half_width_mm: 512_000,
            sim_time_kyr: 40_000,
            dt_kyr: 20,
            stats_every_kyr: 1_000,
            uplift_gain_permille: 440,
            roughness_rms_fraction_ppm: 5_000,
            roughness_rms_cap_mm: 5_000,
            roughness_realization: 0,
            base_floor_above_sea_mm: 1_000,
            initial_base_mode: 1,
            inflow_mode: 1,
            river_threshold_milli_m3_s: 10_000,
            valley_area_m2: 1_000_000,
            local_relief_half_width_posts: 31,
            uplift_rule: "t0_subgrid_closure_rate_v1: U_c = gain x K x R_sub,c / ln(sqrt(cell area_c) / x_c) on T0 land cells (stream-power-subgrid-v1, n = 1); per post, bilinear between the 4 surrounding T0 cell centres with weights renormalised over land cells; 0 if none is land".into(),
            base_rule: "initial_base_mode 0 = t0_bilinear_floored_v1: T0 elevation_m bilinear between L-grid cell centres at the post's (lat, lon); 1 = t1a_trunk_level_bilinear_v1: the T1a trunk_level_m field (bed level of the channel each lod-8 node drains to) bilinear between lod-8 nodes; both floored at sea level + base_floor_above_sea on land posts".into(),
            roughness_rule: "keyed_uniform_white_v1: one draw per lod-4 post, vz.keyed-rng/v1 root(seed) -> layer(T1b recipe digest) -> chart 0 -> tile(4, owner tile of the post) -> draw(t1b.initial.roughness, counter = local index ly*128+lx, attempt 0); value = (2u - 1) x sqrt(3) x rms, rms = min(fraction x R_sub of the post's T0 cell, cap); 0 on marine and lake cells".into(),
            channel_rule: "t1a_channel_links_v1: every T1a node with channel == 1 and its receiver define a D8 link of 16 lod-4 steps; posts on it (land T0 cells only) are fixed at the T1a trunk levels interpolated linearly along the link; where two links write one post the lower bed is kept (counted)".into(),
            marine_rule: "t0_cell_classification_v1: posts whose T0 cell is marine are fixed at sea level; lake cells at the lake stage (T1a's coast reading)".into(),
            boundary_rule: "closed_ring_t1a_inflow_v1: the outer ring of the computation lattice is free (no diffusive flux and no water across the lattice edge); every T1a hillslope flow edge X -> R with X outside the lattice and R inside (lod-8) adds X's T1a accumulated area (discharge x seconds per year / runoff) as inflow, spread uniformly over the free posts of the 512 m strip it crosses (16 lod-4 posts centred on R along the side X lies beyond): inflow_mode 0 = the ring posts of that strip; 1 = all halo-band posts of that strip (the halo absorbs the inflow as extra contributing area, no line imposed). Replaces dirichlet_ring_initial_v1 of 0.1.0, which made the ring a base level (lane TR run r1, 2026-10-12)".into(),
            routing_rule: "priority_flood_steepest_v1: priority flood from all fixed posts (key (fill, index)); receiver = steepest descent on the filled surface over 8 neighbours (ties: lowest index), flats: flood parent; posts under a fill get no stream-power change (pf-surface rule)".into(),
            area_rule: "planar_post_cell_v1: 32 m x 32 m per post in the chart plane (point-scale error <= 2.1e-4 over the chart, not applied)".into(),
            blend_rule: "smoothstep_seam_v1: final post = sum_k w_k z_k / sum_k w_k over windows k in (ry, rx) order, w_k = s(dx_k) s(dy_k), s(d) = smoothstep((d + h) / 2h) of the signed distance d into window k's core (edge-inclusive), h = seam half width; fixed posts take their fixed value".into(),
            node_cell_rule: "chart_inverse_floor_v1 (same as T1a)".into(),
            guide_mode: 0,
            guide_rule: None,
            lithology_mode: 0,
            lithology: None,
            hillslope_d_ppm_of_t0: None,
            ring_mode: 0, ring_rule: None, hillslope_mode: 0, hillslope_rule: None, creep_ratio_max_ppm: None,
            guide_curve_passes: 0, guide_curve_rule: None, coast_mode: 0, coast_rule: None, coast_plain_ratio_ppm: None,
            net: None,
            rivernet: None,
        }
    }
}

impl ReliefParamsV1 {
    /// Switch on the chart-wide tributary guides (0.5.0).
    pub fn with_guides(mut self) -> Self {
        self.guide_mode = 1;
        self.guide_rule = Some("chart_wide_guides_v1: every guide node of the t1b.tributary_guides layer (upstream digest) and its D8 receiver define a link of 16 lod-4 steps; free posts on it (land T0 cells, not T1a channel) are fixed at the guide beds interpolated linearly along the link; where two guide links write one post the lower bed is kept (counted); no T1a ring inflow (the guides carry the area from outside the window; inflow_mode is ignored)".into());
        self
    }
    /// Switch on exposed lithology and threshold-slope failure (0.6.0, lane T3).
    pub fn with_lithology(mut self) -> Self {
        self.lithology_mode = 1;
        self.lithology = Some(crate::t1b_litho::LithologyV1::standin_v1());
        self
    }
    /// T1b 0.7.0 (lane T3b): guides + lithology + guide-surface ring + nonlinear creep + curved guides + zero-contour coast.
    pub fn with_v4(self) -> Self {
        let mut p = self.with_guides().with_lithology();
        p.ring_mode = 1;
        p.initial_base_mode = 2;
        p.base_rule.push_str("; 2 = guide_surface_bicubic_v1 (0.7.0): the t1b.tributary_guides relaxed 512 m surface (height_mm) by Catmull-Rom bicubic interpolation between lod-8 nodes (C1, so no slope break along the 512 m lines seeds straight valleys), floored at sea level + base floor on land");
        p.ring_rule = Some("guide_surface_ring_v1: every free post on the outer ring of the computation lattice is fixed at the t1b.tributary_guides relaxed 512 m surface (height_mm, fixed nodes at their fixed value) by Catmull-Rom bicubic interpolation between lod-8 nodes (the 0.7.0 initial-base field), floored at sea level + base floor on land; the same function of the global post in every window. Ring posts are outlets (class 3); flux into them is exported".into());
        p.hillslope_mode = 1;
        p.hillslope_rule = Some("nonlinear_creep_roering_v1: q = D S / (1 - (S / Sc)^2) on the 4-neighbour edges, D = the T0 header hillslope diffusivity, S = |grad z| at the edge (axis difference with the averaged cross-axis central difference), Sc = the exposed-lithology critical slope of the edge's upper free post (true gradient), beyond S = creep_ratio_max x Sc the factor 1/(1 - (S/Sc)^2) continues along its tangent (finite, monotone); locally one-dimensional backward Euler (x sweep then y sweep), each sweep fully implicit in the along-line gradient (Newton, tridiagonal Jacobian per run of free posts, cross-gradient and Sc lagged from the start of the sweep); replaces the explicit linear diffusion and the threshold cap of 0.6.0; flux into fixed posts exported".into());
        p.creep_ratio_max_ppm = Some(990_000);
        p.guide_curve_passes = 3;
        p.guide_curve_rule = Some("chaikin_guide_segments_v1: the guide network is cut into segments between break nodes (heads, junctions with >= 2 guide donors, and the non-guide node a guide drains into); each segment's polyline of (lod-8 node position, guide bed) is Chaikin corner-cut guide_curve_passes times with both ends pinned; rasterized at 32 m by samples every <= 8 m to the nearest post (8-connected), bed linear along each piece; posts already channel or water are skipped; where several samples write one post the lowest bed is kept".into());
        p.with_coast(2)
    }
    /// Coast rule switch (0.7.0 keeps 0 reachable for comparison): 0 = T0 cell edge (as 0.4.0-0.6.0), 1 = T0 bilinear zero
    /// contour only (lane T3b build2: reclaimed marine-cell land kept the land-weighted uplift and had no drainage), 2 = zero
    /// contour + reclaimed coastal plain + mouths extended to the new shore (integrator ruling 2026-10-13).
    pub fn with_coast(mut self, mode: u32) -> Self {
        self.coast_mode = mode;
        self.coast_plain_ratio_ppm = None;
        self.coast_rule = match mode {
            1 => Some("t0_zero_contour_coast_v1: on T0 land and marine cells a post is marine (fixed at sea level) iff the T0 elevation bilinear between cell centres (t0_bilinear_m, the surface whose zero contour F4 uses) is <= sea level, else land; lake cells keep the cell rule".into()),
            2 => {
                self.coast_plain_ratio_ppm = Some(40_000);
                Some("t0_zero_contour_coast_v2: marine iff the T0 bilinear <= sea level on T0 land and marine cells (lake cells keep the cell rule). RECLAIMED posts (T0 cell marine, bilinear > sea) are a coastal plain FIXED (code 6, an outlet) at max(sea + coast_plain_ratio x (bilinear - sea), sea + base floor), relaxation uplift 0: a depositional surface T1b does not model, held as declared. MOUTH EXTENSION: every T1a channel link whose receiver lies on a T0 marine cell or is a T1a mouth (kind 3), and every guide segment ending on a T0 marine-cell node, whose end node is land under this rule, continues from that node down the T0 bilinear surface (steepest descent, steps of 16 m, central-difference gradient over 32 m) to the first point with bilinear <= sea; its bed falls linearly with path length from the end node's level to min(that level, sea), capped at half the reclaimed-plain height above sea at each point (the river is incised into the plain, never perched on it); rasterized like the guide curves (8-connected, lowest bed kept) as fixed posts on free posts (channel code where any carried T1a channel passes, else guide code; lowest bed kept); a carried bed at or below sea level makes the post an estuary arm of the sea (fixed marine at sea level); chart-wide, so identical in every window".into())
            }
            _ => None,
        };
        self
    }
    pub fn generator_version(&self) -> &'static str {
        if let Some(r) = &self.rivernet { return r.generator_version(); }
        if self.net.is_some() { return T1B_GENERATOR_VERSION_NET; }
        if self.ring_mode == 1 || self.hillslope_mode == 1 || self.guide_curve_passes > 0 || self.coast_mode >= 1 { return T1B_GENERATOR_VERSION_V4; }
        if self.lithology_mode == 1 { T1B_GENERATOR_VERSION_LITHO } else if self.guide_mode == 1 { T1B_GENERATOR_VERSION_GUIDES } else { T1B_GENERATOR_VERSION }
    }
    pub fn halo_posts(&self) -> i64 { self.halo_mm / T1B_SPACING_MM }
    pub fn dt_my(&self) -> f64 { self.dt_kyr as f64 / 1000.0 }
    pub fn steps(&self) -> usize { (self.sim_time_kyr / self.dt_kyr) as usize }
    pub fn stats_every_steps(&self) -> usize { ((self.stats_every_kyr / self.dt_kyr) as usize).max(1) }
    pub fn threshold_m3_s(&self) -> f64 { self.river_threshold_milli_m3_s as f64 / 1000.0 }
    pub fn validate(&self) -> R<()> {
        if self.halo_mm % T1B_SPACING_MM != 0 || self.halo_mm < 0 { return Err("halo must be a whole number of lod-4 posts".into()); }
        if self.seam_half_width_mm % T1B_SPACING_MM != 0 || self.seam_half_width_mm > self.halo_mm {
            return Err("seam half width must be whole posts and inside the halo".into());
        }
        if self.dt_kyr <= 0 || self.sim_time_kyr < 0 || self.sim_time_kyr % self.dt_kyr != 0 { return Err("sim time must be a whole number of steps".into()); }
        if self.uplift_gain_permille < 0 || self.roughness_rms_fraction_ppm < 0 || self.roughness_rms_cap_mm < 0 { return Err("negative parameter".into()); }
        match (self.lithology_mode, &self.lithology) {
            (0, None) => {}
            (1, Some(l)) => l.validate()?,
            _ => return Err("lithology_mode 1 needs a lithology table (and 0 none)".into()),
        }
        if matches!(self.hillslope_d_ppm_of_t0, Some(v) if v <= 0) { return Err("hillslope_d_ppm_of_t0 must be positive".into()); }
        if self.ring_mode == 1 && self.guide_mode != 1 { return Err("ring_mode 1 needs the guide layer".into()); }
        if self.guide_curve_passes > 0 && self.guide_mode != 1 { return Err("guide curves need the guide layer".into()); }
        if let Some(r) = &self.rivernet { r.validate()?; if self.net.is_some() || self.guide_mode != 0 || self.lithology_mode != 0 || self.hillslope_mode != 0 { return Err("0.9.0 carries its own network and lithology (no net, guides or lithology_mode)".into()); } }
        if let Some(n) = &self.net { n.validate()?; if self.guide_mode != 0 || self.lithology_mode != 0 { return Err("0.8.0 carries its own lithology and no guides".into()); } }
        if self.hillslope_mode == 1 && (self.lithology_mode != 1 || !matches!(self.creep_ratio_max_ppm, Some(v) if v > 0 && v < 1_000_000)) {
            return Err("hillslope_mode 1 needs lithology (Sc) and 0 < creep_ratio_max_ppm < 1e6".into());
        }
        Ok(())
    }
}

/// The erosion law constants T1b uses, read from the T0 terrain header.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct ReliefLawsV1 {
    pub law_id: String,
    pub scheme: String,
    pub stream_power_k_per_my: f64,
    pub stream_power_m: f64,
    pub stream_power_n: f64,
    pub hillslope_diffusivity_m2_per_my: f64,
    pub diffusion_stability_limit: f64,
    pub subgrid_channel_head_length_m: f64,
    pub source: String,
}

/// T0 extras T1b needs beyond `T0`: the law constants and the uplift field
/// (reported for comparison, not used).
pub struct T0Extras {
    pub laws: ReliefLawsV1,
    pub uplift_rate_m_per_my: Vec<f32>,
    /// T0 rock_family, crust_kind and sediment_thickness_m (exposed-lithology stand-in, lane T3).
    pub rock_family: Vec<u8>,
    pub crust_kind: Vec<u8>,
    pub sediment_m: Vec<f32>,
}

impl T0Extras {
    pub fn load(run_dir: &Path, t0: &T0) -> R<T0Extras> {
        use pf_contracts::artifact::{read_artifact_header, read_field_f32};
        use serde_json::Value;
        let tdir = run_dir.join("terrain");
        let th = read_artifact_header::<Value>(&tdir).map_err(|e| e.to_string())?;
        let d = &th.header.domain;
        let e = d.get("erosion").ok_or("terrain header: no erosion block")?;
        let f = |k: &str| e.get(k).and_then(Value::as_f64).ok_or_else(|| format!("terrain erosion.{k} missing"));
        let xc = d.get("assumptions").and_then(Value::as_array)
            .and_then(|a| a.iter().find(|x| x.get("id").and_then(Value::as_str) == Some("A-SUBGRID-CHANNEL-HEAD")))
            .and_then(|x| x.get("value")).and_then(Value::as_f64).ok_or("terrain assumption A-SUBGRID-CHANNEL-HEAD missing")?;
        let laws = ReliefLawsV1 {
            law_id: "stream-power-n1-linear-diffusion-v0 (planar lod-4 re-implementation, lane TR)".into(),
            scheme: format!("stream power implicit receiver-first (Braun and Willett 2013); explicit 4-neighbour finite-volume diffusion with substeps; T0 scheme id: {}",
                e.get("stream_power_scheme").and_then(Value::as_str).unwrap_or("?")),
            stream_power_k_per_my: f("stream_power_k_per_my")?,
            stream_power_m: f("stream_power_m")?,
            stream_power_n: f("stream_power_n")?,
            hillslope_diffusivity_m2_per_my: f("hillslope_diffusivity_m2_per_my")?,
            diffusion_stability_limit: f("diffusion_stability_limit")?,
            subgrid_channel_head_length_m: xc,
            source: format!("T0 terrain artifact {} header domain.erosion and assumption A-SUBGRID-CHANNEL-HEAD (legacy project assumptions A-K-UNITS, A-HILLSLOPE-D; awaiting technical review Q2/Q3)", t0.terrain_digest),
        };
        if laws.stream_power_n != 1.0 { return Err("T1b implements n = 1 only (the T0 law's exponent)".into()); }
        let up = read_field_f32(&tdir, "uplift_rate_m_per_my").map_err(|e| e.to_string())?;
        if up.len() != t0.nlat * t0.nlon { return Err("uplift field size".into()); }
        let rf = pf_contracts::artifact::read_field_u8(&tdir, "rock_family").map_err(|e| e.to_string())?;
        let ck = pf_contracts::artifact::read_field_u8(&tdir, "crust_kind").map_err(|e| e.to_string())?;
        let sed = read_field_f32(&tdir, "sediment_thickness_m").map_err(|e| e.to_string())?;
        if rf.len() != up.len() || ck.len() != up.len() || sed.len() != up.len() { return Err("rock field sizes".into()); }
        Ok(T0Extras { laws, uplift_rate_m_per_my: up, rock_family: rf, crust_kind: ck, sediment_m: sed })
    }
    pub fn rock(&self) -> crate::t1b_litho::T0Rock<'_> {
        crate::t1b_litho::T0Rock { rock_family: &self.rock_family, crust_kind: &self.crust_kind, sediment_m: &self.sediment_m }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ReliefLatticeV1 {
    pub tilekey: String,
    pub lod: u8,
    pub spacing_mm: i64,
    pub window_rule: String,
    pub order: String,
}

/// DerivationRecipeV1 shape (ADR 044 section 1) for the `t1b.relief` layer.
/// One recipe per layer: windows are keys in its domain, not recipe fields.
#[derive(Clone, Debug, Serialize)]
pub struct ReliefRecipeV1 {
    pub schema: String,
    pub layer: String,
    pub generator: GeneratorRef,
    pub params: ReliefParamsV1,
    pub params_digest: String,
    pub laws: ReliefLawsV1,
    pub universe_digest: String,
    /// Sorted (role, digest).
    pub upstream: Vec<(String, String)>,
    pub seed_scope: SeedScopeV1,
    pub chart: ChartDeclV1,
    pub lattice: ReliefLatticeV1,
    pub key_domain: String,
    pub simtime_label: String,
    pub closure_label: String,
    pub label: String,
}

impl ReliefRecipeV1 {
    pub fn digest(&self) -> R<String> {
        pf_contracts::artifact::digest_canonical(self).map_err(|e| format!("recipe digest: {e}"))
    }
}

pub struct ReliefInputs<'a> {
    pub t0: &'a T0,
    pub extras: &'a T0Extras,
    pub site: &'a SiteRecord,
    pub skel: &'a SkeletonView,
    /// Chart-wide tributary guides (required iff params.guide_mode == 1).
    pub guides: Option<&'a GuideView>,
    pub params: ReliefParamsV1,
}

pub fn recipe_for(inp: &ReliefInputs) -> R<ReliefRecipeV1> {
    inp.params.validate()?;
    let params_digest = pf_contracts::artifact::digest_canonical(&inp.params).map_err(|e| e.to_string())?;
    let mut upstream = vec![
        ("site".to_string(), inp.site.bytes_digest.clone()),
        ("t0.hydrology".to_string(), inp.t0.hydrology_digest.clone()),
        ("t0.subgrid_relief_sidecar".to_string(), inp.t0.subgrid_digest.clone()),
        ("t0.terrain".to_string(), inp.t0.terrain_digest.clone()),
        ("t1a.skeleton".to_string(), inp.skel.digest.clone()),
    ];
    if inp.params.guide_mode == 1 {
        let g = inp.guides.ok_or("guide_mode 1 needs the t1b.tributary_guides layer")?;
        upstream.push(("t1b.tributary_guides".to_string(), g.digest.clone()));
    } else if inp.guides.is_some() {
        return Err("guides given but guide_mode is 0".into());
    }
    upstream.sort();
    Ok(ReliefRecipeV1 {
        schema: T1B_RECIPE_SCHEMA.into(),
        layer: T1B_LAYER.into(),
        generator: GeneratorRef { name: T1B_GENERATOR_NAME.into(), version: inp.params.generator_version().into(), tier: "Preview".into() },
        params: inp.params.clone(),
        params_digest,
        laws: inp.extras.laws.clone(),
        universe_digest: inp.site.universe_digest.clone(),
        upstream,
        seed_scope: SeedScopeV1 {
            keyed_rng: vz_keyed_rng::KEYED_RNG_V1.into(),
            seed: inp.site.seed,
            ladder: vec!["layer = this recipe's digest".into(), "chart = 0".into(),
                "tile = (4, ix, iy), the lod-4 tile owning the post (half-open floor)".into()],
            purposes: vec![PURPOSE_ROUGHNESS.into()],
        },
        chart: ChartDeclV1::new(inp.site.chart_origin_lat_e7, inp.site.chart_origin_lon_e7, (inp.t0.radius_m * 1000.0).round() as i64),
        lattice: ReliefLatticeV1 {
            tilekey: vz_tilekey::TILEKEY_V2.into(),
            lod: T1B_LOD,
            spacing_mm: T1B_SPACING_MM,
            window_rule: "window(rx, ry) = region lod-6 key core posts [512 rx, 512 rx + 512] x [512 ry, 512 ry + 512] (lod-4 global post indices, edge-inclusive) plus halo_mm on every side".into(),
            order: "post = (gy - gy0) * n + (gx - gx0); x east, y north; chart mm = g * 32000".into(),
        },
        key_domain: format!("{T1B_LAYER}/c0/{REGION_LOD}.*"),
        simtime_label: SIMTIME_LABEL.into(),
        closure_label: CLOSURE_LABEL.into(),
        label: T1B_LABEL.into(),
    })
}

/// The computation lattice of one window (lod-4 global post indices).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct WinLat {
    pub rx: i64,
    pub ry: i64,
    pub gx0: i64,
    pub gy0: i64,
    pub n: usize,
    pub halo_posts: i64,
}

impl WinLat {
    pub fn new(rx: i64, ry: i64, halo_posts: i64) -> WinLat {
        WinLat { rx, ry, gx0: rx * REGION_CELLS - halo_posts, gy0: ry * REGION_CELLS - halo_posts,
            n: (REGION_CELLS + 2 * halo_posts + 1) as usize, halo_posts }
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
    pub fn core(&self) -> (i64, i64, i64, i64) {
        (self.rx * REGION_CELLS, self.rx * REGION_CELLS + REGION_CELLS, self.ry * REGION_CELLS, self.ry * REGION_CELLS + REGION_CELLS)
    }
    #[inline]
    pub fn in_core(&self, i: usize) -> bool {
        let (gx, gy) = self.g(i);
        let (x0, x1, y0, y1) = self.core();
        gx >= x0 && gx <= x1 && gy >= y0 && gy <= y1
    }
    pub fn on_ring(&self, i: usize) -> bool {
        let (x, y) = (i % self.n, i / self.n);
        x == 0 || y == 0 || x == self.n - 1 || y == self.n - 1
    }
    pub fn key(&self) -> String { format!("{REGION_LOD}.{}.{}", self.rx, self.ry) }
}

/// Everything a window run starts from (all window-independent per post).
pub struct WindowInputs {
    pub lat: WinLat,
    pub t0_cell: Vec<u32>,
    pub bilinear_m: Vec<f64>,
    pub base_m: Vec<f64>,
    pub roughness_m: Vec<f64>,
    pub uplift_m_per_my: Vec<f64>,
    pub fixed: Vec<u8>,
    pub z0: Vec<f64>,
    /// T1a inflow drainage area added at ring posts (m2).
    pub inflow_m2: Vec<f64>,
    /// Channel posts: the T1a lod-8 node starting the link (else NONE).
    pub chan_link_node: Vec<u32>,
    /// Exposed-lithology stand-in (lithology_mode 1).
    pub litho: Option<crate::t1b_litho::CoverOverBasementV1>,
    pub stats: PrepStats,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct PrepStats {
    pub t1a_tiles_verified: Vec<String>,
    pub channel_links: u32,
    pub channel_posts: u32,
    pub channel_post_conflicts: u32,
    pub channel_post_conflict_max_m: f64,
    pub soft_axis_nodes_in_window: u32,
    pub marine_posts: u32,
    pub lake_posts: u32,
    pub ring_posts: u32,
    pub free_posts: u32,
    pub inflow_edges: u32,
    /// T1a inflow area by side [W, E, S, N] (km2).
    pub inflow_area_by_side_km2: [f64; 4],
    pub inflow_area_unplaced_km2: f64,
    pub roughness_rms_m_max: f64,
    pub roughness_abs_max_m: f64,
    pub uplift_core_mean_m_per_my: f64,
    pub t0_uplift_core_mean_m_per_my: f64,
    pub cells: Vec<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub guides: Option<GuidePrepStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lithology: Option<crate::t1b_litho::LithoPrepStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mouth_extensions: Option<MouthPrepStats>,
    /// T1b 0.8.0 fill statistics (absent otherwise).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub net: Option<serde_json::Value>,
}

/// Mouth extensions to the zero-contour shore (coast_mode 2).
#[derive(Clone, Debug, Default, Serialize)]
pub struct MouthPrepStats { pub mouths_chart: u32, pub mouths_in_window: u32, pub posts: u32, pub estuary_posts: u32, pub max_length_km: f64, pub not_reaching_sea: u32 }

#[derive(Clone, Debug, Default, Serialize)]
pub struct GuidePrepStats {
    pub guide_links: u32,
    pub guide_posts: u32,
    pub guide_post_conflicts: u32,
    pub guide_post_conflict_max_m: f64,
    pub guide_posts_under_channel_or_water: u32,
}

/// Bilinear T0 field at (lat, lon) degrees between L-grid cell centres.
/// `land_only`: weights renormalised over land cells (0 if none).
pub fn bilinear(t0: &T0, field: &dyn Fn(u32) -> f64, lat: f64, lon: f64, land_only: bool) -> f64 {
    let dlat = 180.0 / t0.nlat as f64;
    let dlon = 360.0 / t0.nlon as f64;
    let fi = (lat + 90.0) / dlat - 0.5;
    let fj = (lon + 180.0) / dlon - 0.5;
    let i0 = fi.floor();
    let j0 = fj.floor();
    let (ti, tj) = (fi - i0, fj - j0);
    let mut acc = 0.0;
    let mut wsum = 0.0;
    for (di, wi) in [(0i64, 1.0 - ti), (1, ti)] {
        let ii = (i0 as i64 + di).clamp(0, t0.nlat as i64 - 1) as usize;
        for (dj, wj) in [(0i64, 1.0 - tj), (1, tj)] {
            let jj = (j0 as i64 + dj).rem_euclid(t0.nlon as i64) as usize;
            let c = (ii * t0.nlon + jj) as u32;
            if land_only && t0.water(c) != CellWater::Land { continue; }
            let w = wi * wj;
            acc += w * field(c);
            wsum += w;
        }
    }
    if land_only { if wsum > 0.0 { acc / wsum } else { 0.0 } } else { acc }
}

/// The T1a trunk level bilinear between lod-8 nodes at lod-4 post (gx, gy);
/// nodes off the T1a lattice fall back to sea level.
pub fn trunk_bilinear(sk: &SkeletonView, gx: i64, gy: i64, sea: f64) -> f64 {
    let (x8, y8) = (gx.div_euclid(L8_TO_L4), gy.div_euclid(L8_TO_L4));
    let (tx, ty) = (gx.rem_euclid(L8_TO_L4) as f64 / L8_TO_L4 as f64, gy.rem_euclid(L8_TO_L4) as f64 / L8_TO_L4 as f64);
    let v = |x: i64, y: i64| sk.idx(x, y).map(|n| sk.trunk_level_m[n] as f64).unwrap_or(sea);
    let a = v(x8, y8) * (1.0 - tx) + v(x8 + 1, y8) * tx;
    let b = v(x8, y8 + 1) * (1.0 - tx) + v(x8 + 1, y8 + 1) * tx;
    a * (1.0 - ty) + b * ty
}

/// The T0 subgrid-closure rate of a land cell (m/My) before the gain.
pub fn closure_rate_m_per_my(t0: &T0, laws: &ReliefLawsV1, c: u32) -> f64 {
    if t0.water(c) != CellWater::Land { return 0.0; }
    let dx = t0.cell_area_m2(c).sqrt();
    let lambda = (dx / laws.subgrid_channel_head_length_m).ln();
    laws.stream_power_k_per_my * (t0.subgrid_relief_m[c as usize] as f64).max(0.0) / lambda
}

/// Keyed roughness for one post (window-independent).
pub struct RoughnessKeys { layer: Scope, tiles: BTreeMap<(i64, i64), Scope> }

impl RoughnessKeys {
    pub fn new(seed: u64, recipe_digest: &str) -> R<RoughnessKeys> {
        let layer = Scope::root(seed).layer(recipe_digest).map_err(|e| e.to_string())?.chart(0).map_err(|e| e.to_string())?;
        Ok(RoughnessKeys { layer, tiles: BTreeMap::new() })
    }
    /// Uniform in [-1, 1) for global lod-4 post (gx, gy).
    pub fn unit(&mut self, gx: i64, gy: i64) -> R<f64> {
        let (tx, ty) = (gx.div_euclid(TILE_CELLS), gy.div_euclid(TILE_CELLS));
        let (lx, ly) = (gx - tx * TILE_CELLS, gy - ty * TILE_CELLS);
        if !self.tiles.contains_key(&(tx, ty)) {
            let s = self.layer.tile(T1B_LOD, tx, ty).map_err(|e| e.to_string())?;
            self.tiles.insert((tx, ty), s);
        }
        let u = self.tiles[&(tx, ty)].draw(PURPOSE_ROUGHNESS, (ly * TILE_CELLS + lx) as u64, 0).map_err(|e| e.to_string())?.uniform_f64();
        Ok(2.0 * u - 1.0)
    }
}

/// Assemble the window inputs. Pure function of (inputs, recipe, window key).
pub fn prepare(inp: &ReliefInputs, recipe: &ReliefRecipeV1, recipe_digest: &str, rx: i64, ry: i64) -> R<WindowInputs> {
    let p = &inp.params;
    let t0 = inp.t0;
    let lat = WinLat::new(rx, ry, p.halo_posts());
    let nn = lat.nodes();
    let chart = Chart::new(recipe.chart.clone()).map_err(|e| e.to_string())?;
    // The window must lie inside the chart and inside the T1a lattice.
    let lo_mm = [lat.gx0 * T1B_SPACING_MM, lat.gy0 * T1B_SPACING_MM];
    let hi_mm = [(lat.gx0 + lat.n as i64 - 1) * T1B_SPACING_MM, (lat.gy0 + lat.n as i64 - 1) * T1B_SPACING_MM];
    let half = vz_tilekey::CHART_HALF_EXTENT_MM;
    if lo_mm[0] < -half || lo_mm[1] < -half || hi_mm[0] > half || hi_mm[1] > half {
        return Err(format!("window {} with halo leaves the chart; refusing", lat.key()));
    }
    let mut stats = PrepStats::default();
    stats.t1a_tiles_verified = inp.skel.verify_tiles_for_rect(lo_mm, hi_mm)?;

    let sea = t0.sea_level_m;
    let floor = sea + p.base_floor_above_sea_mm as f64 / 1000.0;
    let gain = p.uplift_gain_permille as f64 / 1000.0;
    let mut keys = RoughnessKeys::new(inp.site.seed, recipe_digest)?;
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let rate = |c: u32| closure_rate_m_per_my(t0, &inp.extras.laws, c);
    let t0_up = |c: u32| inp.extras.uplift_rate_m_per_my[c as usize] as f64;

    let mut t0_cell = vec![0u32; nn];
    let mut bil = vec![0f64; nn];
    let mut base = vec![0f64; nn];
    let mut rough = vec![0f64; nn];
    let mut up = vec![0f64; nn];
    let mut fixed = vec![FREE; nn];
    let mut z0 = vec![0f64; nn];
    let mut t0_up_core = (0.0, 0usize);
    let mut cells = std::collections::BTreeSet::new();
    for i in 0..nn {
        let (gx, gy) = lat.g(i);
        let (la, lo) = chart.inverse_mm(gx * T1B_SPACING_MM, gy * T1B_SPACING_MM).map_err(|e| e.to_string())?;
        let (la, lo) = (la * RAD2DEG, lo * RAD2DEG);
        let c = t0.cell_of_deg(la, lo);
        cells.insert(c);
        t0_cell[i] = c;
        bil[i] = bilinear(t0, &elev, la, lo, false);
        up[i] = gain * bilinear(t0, &rate, la, lo, true);
        let water = match (p.coast_mode, t0.water(c)) {
            (1 | 2, CellWater::Land) | (1 | 2, CellWater::Marine) => if bil[i] > sea { CellWater::Land } else { CellWater::Marine },
            (_, w) => w,
        };
        let reclaimed = p.coast_mode == 2 && water == CellWater::Land && t0.water(c) == CellWater::Marine;
        if reclaimed { up[i] = 0.0; }
        match water {
            CellWater::Land => {
                base[i] = if reclaimed { sea + p.coast_plain_ratio_ppm.unwrap_or(0) as f64 / 1e6 * (bil[i] - sea) } else { match p.initial_base_mode {
                    1 => trunk_bilinear(inp.skel, gx, gy, sea),
                    2 => inp.guides.ok_or("initial_base_mode 2 needs the guide layer")?.surface_bicubic_m(gx, gy),
                    _ => bil[i],
                } }.max(floor);
                let rms = (t0.subgrid_relief_m[c as usize] as f64).max(0.0) * p.roughness_rms_fraction_ppm as f64 / 1e6;
                let rms = rms.min(p.roughness_rms_cap_mm as f64 / 1000.0);
                stats.roughness_rms_m_max = stats.roughness_rms_m_max.max(rms);
                let a = 3f64.sqrt() * rms;
                rough[i] = if a > 0.0 { a * keys.unit(gx, gy)? } else { 0.0 };
                stats.roughness_abs_max_m = stats.roughness_abs_max_m.max(rough[i].abs());
                z0[i] = base[i] + rough[i];
            }
            CellWater::Marine => { base[i] = sea; fixed[i] = MARINE; z0[i] = sea; }
            CellWater::Lake => { base[i] = t0.level(c); fixed[i] = LAKE; z0[i] = base[i]; }
        }
        if lat.in_core(i) && t0.water(c) == CellWater::Land { t0_up_core.0 += t0_up(c); t0_up_core.1 += 1; }
    }
    stats.cells = cells.into_iter().collect();

    // T1a channel links, rasterized at 32 m (D8 links of 16 lod-4 steps).
    let sk = inp.skel;
    let mut chan_link_node = vec![NONE; nn];
    let g8_lo = (lat.gx0.div_euclid(L8_TO_L4) - 1, lat.gy0.div_euclid(L8_TO_L4) - 1);
    let g8_hi = ((lat.gx0 + lat.n as i64).div_euclid(L8_TO_L4) + 1, (lat.gy0 + lat.n as i64).div_euclid(L8_TO_L4) + 1);
    let mut chan_level: BTreeMap<usize, f64> = BTreeMap::new();
    for gy8 in g8_lo.1..=g8_hi.1 { for gx8 in g8_lo.0..=g8_hi.0 {
        let Some(n8) = sk.idx(gx8, gy8) else { continue };
        if sk.kind[n8] == 1 && sk.channel[n8] != 1 {
            let inside = lat.idx(gx8 * L8_TO_L4, gy8 * L8_TO_L4).is_some();
            if inside { stats.soft_axis_nodes_in_window += 1; }
        }
        if !(sk.kind[n8] == 1 && sk.channel[n8] == 1) { continue; }
        let r8 = sk.receiver[n8];
        if r8 == NONE || r8 == OFF_CHART { continue; }
        let (rgx, rgy) = node_g(r8 as usize);
        let (dx, dy) = (rgx - gx8, rgy - gy8);
        if dx.abs() > 1 || dy.abs() > 1 { return Err(format!("T1a link {n8} -> {r8} is not a D8 step")); }
        let (la, lb) = (sk.trunk_level_m[n8] as f64, sk.trunk_level_m[r8 as usize] as f64);
        let mut any = false;
        for s in 0..=L8_TO_L4 {
            let (gx, gy) = (gx8 * L8_TO_L4 + dx * s, gy8 * L8_TO_L4 + dy * s);
            let Some(i) = lat.idx(gx, gy) else { continue };
            any = true;
            if fixed[i] == MARINE || fixed[i] == LAKE { continue; }
            let lev = la + (lb - la) * s as f64 / L8_TO_L4 as f64;
            match chan_level.get(&i).copied() {
                Some(old) if old != lev => {
                    stats.channel_post_conflicts += 1;
                    stats.channel_post_conflict_max_m = stats.channel_post_conflict_max_m.max((old - lev).abs());
                    if lev < old { chan_level.insert(i, lev); chan_link_node[i] = n8 as u32; }
                }
                Some(_) => {}
                None => { chan_level.insert(i, lev); chan_link_node[i] = n8 as u32; }
            }
        }
        if any { stats.channel_links += 1; }
    } }
    for (&i, &lev) in &chan_level { fixed[i] = CHANNEL; z0[i] = lev; }
    // Mouth extensions to the zero-contour shore (t0_zero_contour_coast_v2).
    if p.coast_mode == 2 {
        let ext = crate::t1b_curves::mouth_extensions(t0, inp.skel, inp.guides, &chart, p.guide_curve_passes, p.coast_plain_ratio_ppm.unwrap_or(0) as f64 / 1e6)?;
        let mut ext_level: BTreeMap<usize, (f64, u32, u8)> = BTreeMap::new();
        let mut es = MouthPrepStats { mouths_chart: ext.len() as u32, ..Default::default() };
        for e in &ext {
            let mut touched = false;
            crate::t1b_curves::rasterize(&e.curve, &mut |gx, gy, lev| {
                let Some(i) = lat.idx(gx, gy) else { return };
                if fixed[i] != FREE { return; }
                touched = true;
                // Channel mouths are T1a channels carried on (code 1); guide mouths stay guides (code 5). Where both write a
                // post (descent paths merge) it is a channel post (lowest bed; the channel's T1a node), so the carried
                // T1a channel stays continuous.
                let code = if e.from_channel { CHANNEL } else { GUIDE };
                let next = match ext_level.get(&i).copied() {
                    None => (lev, e.curve.terminal, code),
                    Some((old, n_old, c_old)) => {
                        let (n, c) = if c_old == CHANNEL && code != CHANNEL { (n_old, c_old) } else if code == CHANNEL && c_old != CHANNEL { (e.curve.terminal, code) } else if lev < old { (e.curve.terminal, code) } else { (n_old, c_old) };
                        (old.min(lev), n, c)
                    }
                };
                ext_level.insert(i, next);
            });
            if touched { es.mouths_in_window += 1; es.max_length_km = es.max_length_km.max(e.length_m / 1000.0); }
            if !e.reached_sea { es.not_reaching_sea += 1; }
        }
        es.posts = ext_level.len() as u32;
        // A carried bed at or below sea level is an estuary arm of the sea (fixed marine at sea level): the T1a mouth sits at
        // sea level on the old coast, so the river meets the sea there and the reclaimed plain is cut by a tidal inlet.
        for (&i, &(lev, node, code)) in &ext_level {
            if lev <= sea { fixed[i] = MARINE; z0[i] = sea; es.estuary_posts += 1; } else { fixed[i] = code; z0[i] = lev; chan_link_node[i] = node; }
        }
        stats.mouth_extensions = Some(es);
    }
    // Reclaimed coastal plain (coast_mode 2): held at its declared surface. Without uplift the detachment-limited law would
    // plane it to sea level within ~1 My (lane T3b run final_g620 w0m2, build7: whole window at 0 m); real coastal plains are
    // depositional, which T1b does not model, so the plain is a fixed surface and an outlet for the land above it.
    if p.coast_mode == 2 {
        for i in 0..nn {
            if fixed[i] == FREE && t0.water(t0_cell[i]) == CellWater::Marine && bil[i] > sea { fixed[i] = PLAIN; z0[i] = base[i]; }
        }
    }
    // Chart-wide tributary guides (chart_wide_guides_v1).
    if p.guide_mode == 1 {
        let gv = inp.guides.ok_or("guide_mode 1 without guides")?;
        let mut gs = GuidePrepStats::default();
        let mut guide_level: BTreeMap<usize, (f64, u32)> = BTreeMap::new();
        if p.guide_curve_passes > 0 {
            // chaikin_guide_segments_v1 (0.7.0).
            let curves = crate::t1b_curves::guide_curves(gv, p.guide_curve_passes);
            let (lx0, ly0, lx1, ly1) = (lat.gx0 - 1, lat.gy0 - 1, lat.gx0 + lat.n as i64, lat.gy0 + lat.n as i64);
            let mut under = std::collections::BTreeSet::new();
            let mut first_writer: BTreeMap<usize, usize> = BTreeMap::new();
            let mut multi = std::collections::BTreeSet::new();
            for (ci, cv) in curves.iter().enumerate() {
                let (mut bx0, mut bx1, mut by0, mut by1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
                for q in &cv.pts { bx0 = bx0.min(q[0]); bx1 = bx1.max(q[0]); by0 = by0.min(q[1]); by1 = by1.max(q[1]); }
                let inside = bx1 >= lx0 as f64 && bx0 <= lx1 as f64 && by1 >= ly0 as f64 && by0 <= ly1 as f64;
                if !inside { continue; }
                gs.guide_links += 1;
                crate::t1b_curves::rasterize(cv, &mut |gx, gy, lev| {
                    let Some(i) = lat.idx(gx, gy) else { return };
                    if fixed[i] != FREE { under.insert(i); return; }
                    match first_writer.get(&i) { Some(&w) if w != ci => { multi.insert(i); } None => { first_writer.insert(i, ci); } _ => {} }
                    match guide_level.get(&i).copied() {
                        Some((old, _)) if lev >= old => {}
                        Some((old, _)) => { if multi.contains(&i) { gs.guide_post_conflict_max_m = gs.guide_post_conflict_max_m.max(old - lev); } guide_level.insert(i, (lev, cv.terminal)); }
                        None => { guide_level.insert(i, (lev, cv.terminal)); }
                    }
                });
            }
            gs.guide_post_conflicts = multi.len() as u32;
            gs.guide_posts_under_channel_or_water = under.len() as u32;
        }
        for gy8 in g8_lo.1..=g8_hi.1 { for gx8 in g8_lo.0..=g8_hi.0 {
            if p.guide_curve_passes > 0 { break; }
            let Some(n8) = sk.idx(gx8, gy8) else { continue };
            if !gv.is_guide(n8) { continue; }
            let r8 = gv.receiver[n8];
            if r8 == NONE { return Err(format!("guide node {n8} has no receiver")); }
            let (rgx, rgy) = node_g(r8 as usize);
            let (dx, dy) = (rgx - gx8, rgy - gy8);
            if dx.abs() > 1 || dy.abs() > 1 { return Err(format!("guide link {n8} -> {r8} is not a D8 step")); }
            let (la, lb) = (gv.bed_m(n8), gv.bed_m(r8 as usize));
            let mut any = false;
            for s in 0..=L8_TO_L4 {
                let (gx, gy) = (gx8 * L8_TO_L4 + dx * s, gy8 * L8_TO_L4 + dy * s);
                let Some(i) = lat.idx(gx, gy) else { continue };
                any = true;
                if fixed[i] != FREE { gs.guide_posts_under_channel_or_water += 1; continue; }
                let lev = la + (lb - la) * s as f64 / L8_TO_L4 as f64;
                match guide_level.get(&i).copied() {
                    Some((old, _)) if old != lev => {
                        gs.guide_post_conflicts += 1;
                        gs.guide_post_conflict_max_m = gs.guide_post_conflict_max_m.max((old - lev).abs());
                        if lev < old { guide_level.insert(i, (lev, gv.terminal[n8])); }
                    }
                    Some(_) => {}
                    None => { guide_level.insert(i, (lev, gv.terminal[n8])); }
                }
            }
            if any { gs.guide_links += 1; }
        } }
        for (&i, &(lev, term)) in &guide_level { fixed[i] = GUIDE; z0[i] = lev; chan_link_node[i] = term; }
        gs.guide_posts = guide_level.len() as u32;
        stats.guides = Some(gs);
    }
    // Closed ring with T1a inflow (closed_ring_t1a_inflow_v1).
    let mut inflow = vec![0f64; nn];
    let b8 = (lat.gx0.div_euclid(L8_TO_L4) + if lat.gx0.rem_euclid(L8_TO_L4) == 0 { 0 } else { 1 },
              (lat.gx0 + lat.n as i64 - 1).div_euclid(L8_TO_L4),
              lat.gy0.div_euclid(L8_TO_L4) + if lat.gy0.rem_euclid(L8_TO_L4) == 0 { 0 } else { 1 },
              (lat.gy0 + lat.n as i64 - 1).div_euclid(L8_TO_L4));
    let inside8 = |x: i64, y: i64| x >= b8.0 && x <= b8.1 && y >= b8.2 && y <= b8.3;
    let (xr0, xr1, yr0, yr1) = (lat.gx0, lat.gx0 + lat.n as i64 - 1, lat.gy0, lat.gy0 + lat.n as i64 - 1);
    for y8 in b8.2 - 1..=b8.3 + 1 { for x8 in b8.0 - 1..=b8.1 + 1 {
        if inside8(x8, y8) { continue; }
        let Some(nx) = sk.idx(x8, y8) else { continue };
        if !(sk.kind[nx] == 0 || sk.kind[nx] == 4) { continue; }
        let r = sk.receiver[nx];
        if r == NONE || r == OFF_CHART { continue; }
        let (rx8, ry8) = node_g(r as usize);
        if !inside8(rx8, ry8) { continue; }
        let area = sk.discharge_m3_s[nx] as f64 * SECONDS_PER_YEAR / t0.runoff_m_per_yr;
        let side = if x8 < b8.0 { 0 } else if x8 > b8.1 { 1 } else if y8 < b8.2 { 2 } else { 3 };
        let (cx0, cx1, cy0, cy1) = lat.core();
        if p.guide_mode == 1 { continue; } // guides carry the outside area; no ring inflow (chart_wide_guides_v1)
        let depth: Vec<i64> = if p.inflow_mode == 1 { (0..lat.halo_posts.max(1)).collect() } else { vec![0] };
        let posts: Vec<usize> = (-8..8).flat_map(|k| depth.iter().map(move |&d| (k, d))).filter_map(|(k, d)| {
            let (gx, gy) = match side {
                0 => (xr0 + d, ry8 * L8_TO_L4 + k), 1 => (xr1 - d, ry8 * L8_TO_L4 + k),
                2 => (rx8 * L8_TO_L4 + k, yr0 + d), _ => (rx8 * L8_TO_L4 + k, yr1 - d),
            };
            let in_band = match side { 0 => gx < cx0, 1 => gx > cx1, 2 => gy < cy0, _ => gy > cy1 };
            if !in_band { return None; }
            lat.idx(gx, gy).filter(|&i| fixed[i] == FREE)
        }).collect();
        stats.inflow_edges += 1;
        stats.inflow_area_by_side_km2[side] += area / 1e6;
        if posts.is_empty() { stats.inflow_area_unplaced_km2 += area / 1e6; continue; }
        let each = area / posts.len() as f64;
        for i in posts { inflow[i] += each; }
    } }
    if p.ring_mode == 1 {
        let gv = inp.guides.ok_or("ring_mode 1 without guides")?;
        for i in 0..nn {
            if !lat.on_ring(i) || fixed[i] != FREE { continue; }
            let (gx, gy) = lat.g(i);
            fixed[i] = RING;
            z0[i] = gv.surface_bicubic_m(gx, gy).max(floor);
        }
    }
    for i in 0..nn {
        match fixed[i] {
            FREE => stats.free_posts += 1,
            CHANNEL => stats.channel_posts += 1,
            MARINE => stats.marine_posts += 1,
            LAKE => stats.lake_posts += 1,
            _ => {}
        }
    }
    let mut ucore = (0.0, 0usize);
    for i in 0..nn { if lat.in_core(i) && fixed[i] == FREE { ucore.0 += up[i]; ucore.1 += 1; } }
    stats.uplift_core_mean_m_per_my = ucore.0 / ucore.1.max(1) as f64;
    stats.t0_uplift_core_mean_m_per_my = t0_up_core.0 / t0_up_core.1.max(1) as f64;
    for i in 0..nn { if lat.on_ring(i) && (fixed[i] == FREE || fixed[i] == RING) { stats.ring_posts += 1; } }
    let litho = match &p.lithology {
        Some(table) if p.lithology_mode == 1 => {
            let pos = |i: usize| -> R<(f64, f64)> {
                let (gx, gy) = lat.g(i);
                let (la, lo) = chart.inverse_mm(gx * T1B_SPACING_MM, gy * T1B_SPACING_MM).map_err(|e| e.to_string())?;
                Ok((la * RAD2DEG, lo * RAD2DEG))
            };
            let l = crate::t1b_litho::build_standin(table, t0, &inp.extras.rock(), nn, &pos)?;
            stats.lithology = Some(crate::t1b_litho::prep_stats(&l, &|i| fixed[i] == FREE && lat.in_core(i)));
            Some(l)
        }
        _ => None,
    };
    Ok(WindowInputs { lat, t0_cell, bilinear_m: bil, base_m: base, roughness_m: rough, uplift_m_per_my: up, fixed, z0, inflow_m2: inflow, chan_link_node, litho, stats })
}
