//! T1b 0.9.0, RIVERNET (lane RV, 2026-10-13; model specification approach A "sounds good to me":
//! river-network-first terrain, cf. Genevaux et al. 2013 "Terrain generation using procedural
//! models based on hydrology"). LABEL: pre-signoff, fidelity tier Preview, derived, never canon.
//! 0.9.0 is a CONSTRUCTION, not a relaxation: nothing is integrated in time, so there is no
//! SimTime and nothing may be cited as an age.
//!
//! Why: the erosion-window approaches (0.4.0-0.8.0) could not pass seams or held-out F6, because
//! n = 1 stream-power relief depends on grid resolution and a per-window (or per-tile) drainage
//! solve disagrees at window edges (drainage spans catchments). Here the drainage network is
//! fixed FIRST and the surface is an explicit function of position:
//!
//! 1. `t1b.rivernet` (chart layer, 512 m, once per chart): water and the T1a channels are seeds;
//!    every other lod-8 node joins a drainage tree grown by invasion percolation (Stark 1991) on
//!    keyed node weights, confined to its T0 cell or T0 receiver cell (F3). Node positions get
//!    keyed jitter; each link a keyed meander midpoint; chains are Chaikin-smoothed (512 m ->
//!    128 m -> 32 m planform). Channels are links with area >= a channel-head area set by rock.
//! 2. Channel beds rise upstream at slope k_s A^-theta (steady-state stream power n = 1,
//!    theta = m/n of the T0 law), anchored at the T1a trunk levels, lake stages and sea level.
//! 3. Hillslopes: an eikonal surface |grad z| = S from channels, valley floors and water, with S
//!    the steady-state nonlinear-creep gradient (Roering et al. 2007) of the exposed rock; ridges
//!    form where slopes from neighbouring channels meet; crests rounded to the linear-diffusion
//!    curvature E/D by a min(z, box mean) rule; a priority-flood fill removes any closed depression.
//! 4. One scalar per T0 cell, the steady-state erosion rate E (hat-blended, so no cell-edge
//!    steps), is solved once per chart on a 128 m lattice so each cell's hat-weighted land mean
//!    meets the T0 bilinear (F6). Every tile reads the same table, so tiles cannot disagree.
//!
//! Every post is a pure function of (seed, T0, T1a, recipe, E table, position); a window computes
//! its core + halo + a declared margin and publishes core + halo, so seams are exact when the
//! margin exceeds the dependency radius (measured, not assumed).

use crate::t0::{CellWater, SiteRecord, T0};
use crate::t1b::{bilinear, closure_rate_m_per_my, ReliefLawsV1, T0Extras, CHANNEL, FREE, GUIDE, LAKE, MARINE};
use crate::t1b_litho::{basement_code, ExposedLithology, LithologyV1};
use crate::t1b_net::hat4;
use crate::t1b_view::SkeletonView;
use crate::*;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};
use std::path::Path;
use vz_keyed_rng::Scope;
use vz_tilekey::chart::{Chart, ChartDeclV1, RAD2DEG};

pub const RN_GENERATOR_NAME: &str = "t1b-rivernet";
pub const RN_GENERATOR_VERSION: &str = "0.9.0";
pub const RN_RECIPE_SCHEMA: &str = "vz.t1b-rivernet-recipe/v1";
pub const RN_SCHEMA: &str = "vz.t1b-rivernet/v1";
pub const RN_LAYER: &str = "t1b.rivernet";
pub const RN_LABEL: &str = "pre-signoff; T0 substrate pre-signoff; fidelity tier Preview; derived, never canon; construction (no SimTime)";
pub const PURPOSE_IP: &str = "t1b.rivernet.route";
pub const PURPOSE_ROCKMASS: &str = "t1b.rivernet.rockmass";
pub const PURPOSE_JITTER: &str = "t1b.rivernet.node.jitter";
pub const PURPOSE_MEANDER: &str = "t1b.rivernet.meander";
pub const LITHO_SRC_STANDIN: &str = "standin_cover_over_basement_v1";
pub const LITHO_SRC_STRATA: &str = "t0_strata_full_v1";

/// Node kinds of the 512 m network.
pub const K_LAND: u8 = 0;
pub const K_T1A: u8 = 1;
pub const K_MARINE: u8 = 2;
pub const K_LAKE: u8 = 3;

/// Declared constants (integers and rule strings only). Pre-registered in
/// `unpublished-design-record` (lane RV); project assumptions until ruled.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct RiverNetParamsV1 {
    /// Draw realization (F12 uses 0 and 1): the attempt segment of every keyed draw.
    pub realization: u32,
    /// 0 invasion percolation (run r0), 1 first-passage percolation (REVISION 1, run r1), 2 random walk down the geodesic
    /// distance to the seeds (Scheidegger-type; exploratory run r2 and later).
    pub network_mode: u32,
    /// Rock-mass strength closure: per channel chain a keyed multiplier exp(sigma N) on the hillslope Sc, clamped (0 = off).
    pub rockmass_sigma_permille: i64,
    pub rockmass_min_permille: i64,
    pub rockmass_max_permille: i64,
    /// E-solve reference: 0 = T0 bilinear (F6 as read by lanes TR-NO), 1 = the T0 cell's own mean (F6 literal) at each land post,
    /// 2 (lane CO, provisional pending design) = the T0 cell mean on coastal T0 land cells (bilinear land fraction in (0.01, 0.99);
    /// REVISION 1: marine T0 cells keep the bilinear), bilinear elsewhere.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub f6_reference: u32,
    pub jitter_permille: i64,
    pub meander_permille: i64,
    pub chaikin_passes: u32,
    pub channel_head_area_m2: i64,
    pub bed_slope_cap_permille_of_sc: i64,
    pub floor_width_coeff_ppm: i64,
    pub floor_width_max_mm: i64,
    pub floor_slope_ppm: i64,
    pub hill_ref_length_mm: i64,
    pub sc_ref_permille: i64,
    pub hill_slope_min_ppm: i64,
    pub crest_radius_max_mm: i64,
    pub e0_gain_permille: i64,
    pub e_min_m_per_my: i64,
    pub e_max_m_per_my: i64,
    pub solve_spacing_mm: i64,
    pub solve_tolerance_mm: i64,
    pub solve_updates_max: u32,
    pub solve_min_land_hat_permille: i64,
    pub margin_mm: i64,
    pub dip_clamp_ppm: i64,
    pub lithology_source: String,
    pub lithology: LithologyV1,
    pub network_rule: String,
    pub planform_rule: String,
    pub channel_rule: String,
    pub bed_rule: String,
    pub floor_rule: String,
    pub hillslope_rule: String,
    pub crest_rule: String,
    pub fill_rule: String,
    pub e_rule: String,
    pub solve_rule: String,
    pub lithology_rule: String,
    pub water_rule: String,
    pub network_mode_rule: String,
    pub rockmass_rule: String,
    pub estuary_rule: String,
    pub bank_rule: String,
    /// T1b 0.9.1 (lane CO): coast cut (shelf, flooding offset per coastal cell, fans, wave cut, E taper); absent = 0.9.0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coast: Option<crate::t1b_coast::CoastParamsV1>,
    /// T1b 0.9.1 (lane CO): volcano edifices from the T0 volcano records; absent = 0.9.0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub volcano: Option<crate::t1b_coast::VolcanoParamsV1>,
    /// T1b 0.9.2 (lane CF): cone scaling and foot blend, flows, gullies, attached fans, damped E/F solve, cliff profile; absent = 0.9.1.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cf: Option<crate::t1b_coast::CfParamsV1>,
    /// T1b 0.9.3 (lane TD): perched-bed and one-post-cliff fixes, over-code publication; absent = 0.9.2.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fx: Option<crate::t1b_fix::FxParamsV1>,
    /// T1b 0.9.4 (lane TA): consistent level solve, T1a bank ramp, Broyden E update; absent = 0.9.3.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ls: Option<crate::t1b_ls::LsParamsV1>,
    /// Lane CLIMATE-RUNOFF (opt-in key climate=<dir>): land node and post own areas weighted by R(x) / r_ref, so every area
    /// is the effective area Q / r_ref (climate_runoff.rs); needs a climate T1a of the same climate. None = unchanged.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub climate_runoff: Option<crate::climate_runoff::ClimateRunoffParamsV1>,
    /// Lane VALLEY-FLOOR-FIX (opt-in keys vf_floor / vf_foot / vf_fan): levee-backswamp floor, colluvial footslope, tributary
    /// fans (t1b_vf.rs). None = unchanged; never seen by the keyed draws (draw_digest strips it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vf: Option<crate::t1b_vf::VfParamsV1>,
    /// Lane HS-PASS (opt-in keys hs=1 / hs_* / hs_volcano): window-only hillslope transient (t1b_hs.rs). None = unchanged;
    /// never seen by the keyed draws (draw_digest strips it) nor by the layer (rnwin checks the layer without it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hs: Option<crate::t1b_hs::HsParamsV1>,
    /// Lane SEAM-CERT (opt-in key seam_cert=1): window-only pinned-only fills, evaluation-lattice routing and certificate
    /// (t1b_seam.rs). None = unchanged; never seen by the keyed draws (draw_digest strips it) nor by the layer.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seam_cert: Option<crate::t1b_seam::SeamCertParamsV1>,
    /// Lane CHANNEL-CENTRELINE (opt-in keys cl=1 / cl_mn=1): window-only T1a-curve centreline and declared meander belt
    /// (t1b_cl.rs). None = unchanged; never seen by the keyed draws (draw_digest strips it) nor by the layer.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cl: Option<crate::t1b_cl::ClParamsV1>,
    /// Lane GLACIAL-T1 (opt-in key glacial=<t0.glacial dir>, layer AND windows): chart glacier table, glacial long profile
    /// of the network beds and the window U-section / cirque carve (t1b_glacial.rs). None = unchanged; never seen by the
    /// keyed draws (draw_digest strips it), so key-on and key-off share one network.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub glacial: Option<crate::t1b_glacial::GlacialParamsV1>,
    /// Lane E-SOLVE-GENERALITY (opt-in key es=1, layer AND windows): secant E update, chart-rim cells, bound-KKT stop
    /// (t1b_es.rs). None = unchanged; never seen by the keyed draws (draw_digest strips it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub es: Option<crate::t1b_es::EsParamsV1>,
    /// Lane DAMMED-VALLEY (opt-in key dam_lake=1): window-only edifice-dammed lakes (t1b_dam.rs). None = unchanged; never
    /// seen by the keyed draws (draw_digest strips it) nor by the layer.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dam_lake: Option<crate::t1b_dam::DamLakeParamsV1>,
    /// Lane EDIFICE-FOOT (opt-in key ef_cap=1): window-only composite foot drape cap (t1b_ef.rs). None = unchanged; never
    /// seen by the keyed draws (draw_digest strips it) nor by the layer.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edifice_foot: Option<crate::t1b_ef::EdificeFootParamsV1>,
    /// Lane T1B-PEDESTAL (opt-in key ped_solve=1, layer AND windows): the pedestal rules (t1b_ped.rs). None = unchanged;
    /// never seen by the keyed draws (draw_digest strips it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pedestal: Option<crate::t1b_ped::PedestalParamsV1>,
    /// Lane T0-IMPRINT-PROD (opt-in key conf=cell|basin|none, layer AND windows; ADR 045): the network confinement rule
    /// (t1b_outlet.rs). None = the cell rule (recipe unchanged); never seen by the keyed draws (draw_digest strips it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confinement: Option<crate::t1b_outlet::ConfinementParamsV1>,
    /// Lane MS-CEIL (opt-in key ceil=<c_h permille>, layer AND windows): ceiling-min v1 (MS-NEXT design), free non-floor
    /// land posts take min(z_ss, Z_c) (CEIL_RULE). None = unchanged; never seen by the keyed draws (draw_digest strips it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ceil: Option<CeilParamsV1>,
    /// Lane T0-IMPRINT-PROD (key outlet_checks=1): adds the ADR 045 per-T0-outlet F3 / whole-cell F6 gates to the checks.
    /// A checker switch, not a recipe parameter: never serialized, so no digest or published byte depends on it.
    #[serde(skip)]
    pub outlet_checks: u32,
    /// Lane DRAIN-FIX (key wall_checks=1): adds the one-post wall census (t1b_checks::wall_census) to each window's
    /// checks.json. A checker switch, not a recipe parameter: never serialized.
    #[serde(skip)]
    pub wall_checks: u32,
    /// Lane LS1-GLOBAL (opt-in keys ls_global=1 / ls_settle=1, layer AND windows): LS1 global order and settle (t1b_ls.rs).
    /// None = unchanged; never seen by the keyed draws (draw_digest strips it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ls_global: Option<crate::t1b_ls::LsGlobalParamsV1>,
}

impl RiverNetParamsV1 {
    /// 0.9.1 when the lane CO additions are present, else 0.9.0 (earlier recipes unchanged).
    /// The network is rerouted in the layer (TD2 capture or LS1) and windows rebuild it from the published receivers.
    pub fn rerouted(&self) -> bool { self.ls.as_ref().map(|l| l.solve == 1).unwrap_or(false) || self.fx.as_ref().map(|f| f.capture == 1).unwrap_or(false) }
    pub fn generator_version(&self) -> &'static str { if self.glacial.is_some() { crate::t1b_glacial::GL_GENERATOR_VERSION } else if self.ls.is_some() { crate::t1b_ls::LS_GENERATOR_VERSION } else if self.fx.is_some() { crate::t1b_fix::FX_GENERATOR_VERSION } else if self.cf.is_some() { crate::t1b_coast::CF_GENERATOR_VERSION } else if self.coast.is_some() || self.volcano.is_some() { crate::t1b_coast::CV_GENERATOR_VERSION } else { RN_GENERATOR_VERSION } }
    pub fn v1(strata: bool) -> Self {
        RiverNetParamsV1 {
            realization: 0,
            network_mode: 2,
            rockmass_sigma_permille: 250,
            rockmass_min_permille: 600,
            rockmass_max_permille: 1500,
            f6_reference: 0,
            jitter_permille: 300,
            meander_permille: 150,
            chaikin_passes: 2,
            channel_head_area_m2: 500_000,
            bed_slope_cap_permille_of_sc: 600,
            floor_width_coeff_ppm: 30_000,
            floor_width_max_mm: 1_200_000,
            floor_slope_ppm: 5_000,
            hill_ref_length_mm: 400_000,
            sc_ref_permille: 700,
            hill_slope_min_ppm: 10_000,
            crest_radius_max_mm: 384_000,
            e0_gain_permille: 1000,
            e_min_m_per_my: 5,
            e_max_m_per_my: 20_000,
            solve_spacing_mm: 128_000,
            solve_tolerance_mm: 20_000,
            solve_updates_max: 14,
            solve_min_land_hat_permille: 100,
            margin_mm: 3_072_000,
            dip_clamp_ppm: 87_489,
            lithology_source: if strata { LITHO_SRC_STRATA.into() } else { LITHO_SRC_STANDIN.into() },
            lithology: LithologyV1::standin_v1(),
            network_rule: "REVISION 1 first_passage_percolation_t0_confined_v1 (replaces invasion percolation): the tree is grown from the seeds by Dijkstra on node costs, entering node j from i costs step x (0.5 + u_j) (step 1 or sqrt 2 lattice units, u_j the keyed uniform of t1b.rivernet.fpp.cost), heap key (cost, node, from), so paths are near-geodesic with keyed wander and coalesce into trees (first-passage percolation); PRE-REVISION TEXT: invasion_percolation_t0_confined_v1: lod-8 nodes (512 m, the T1a lattice); water nodes = T0 lake cells and zero-contour marine (T0 bilinear <= sea); seeds = water + T1a valley-axis land nodes (kind 1 channel and sub-threshold axes, kind 3 mouths; level = T1a trunk level) + sea level / lake stage; every land node joins the tree grown from the seeds by Prim on node weights (weight = keyed u64, vz.keyed-rng/v1 root(seed) -> layer(rivernet recipe digest) -> chart 0 -> tile(8, owner lod-8 tile) -> draw(t1b.rivernet.ip.priority, counter = local index ly*128+lx, attempt = realization); heap key (weight, node, from)); a node may be entered from a node of its own T0 cell, of its T0 receiver cell, or from water; a diagonal link is refused while the crossing diagonal is a link; nodes left unreached are then entered without the cell rule (counted). Basis (recalled, unverified): Stark 1991 Nature 352:423, invasion percolation grows natural-looking dendritic networks".into(),
            planform_rule: "jitter_meander_chaikin_v1: land nodes move by (2u-1) x jitter x 512 m per axis (t1b.rivernet.node.jitter, counters 2k, 2k+1); T1a and water nodes stay on the lattice; each land link gets one midpoint moved sideways by (2u-1) x meander x link length (t1b.rivernet.meander); channel chains (between heads, junctions with >= 2 channel donors, and the seed or junction they end on) are Chaikin corner-cut chaikin_passes times with both ends pinned (Chaikin 1974); bed and area are carried linearly. T1a valley-axis links (kind 1; channel==1 published as T1a channels, sub-threshold axes as network channels): Catmull-Rom through the T1a nodes (previous = largest-discharge channel donor or the reflected receiver; next = the receiver's receiver or the reflected node), 16 samples, bed linear in the parameter".into(),
            channel_rule: "rock_channel_head_v1: a link is a channel iff its area >= channel_head_area / k_mult(rock at the node at its receiver's bed), or a channel drains into it (channels never stop downstream); rasterized at the lattice spacing by samples every <= 1/4 post to the nearest post; T1a links win over network links at a post, otherwise the lowest bed wins (first curve in chart order on ties); published code 1 (T1a) or 5 (network channel)".into(),
            bed_rule: "steady_state_stream_power_v1: bed(node) = bed(receiver) + L x min(k_s A^-theta, cap x Sc), k_s = E / (K_T0 x k_mult), theta = m/n of the T0 law, A = network area (m2), L = jittered path length through the meander midpoint, Sc and k_mult of the rock exposed at the node at the receiver's bed; seeds keep their levels; never rises along flow".into(),
            floor_rule: "nearest_channel_floor_v1: a free land post within W/2 of a channel post (W = floor_width_coeff x sqrt(A), at most floor_width_max) takes the bed of its NEAREST such channel post (ties: lower global post key) + max(floor_slope, that post's bed gradient) x distance (REVISION 1: the bed-gradient floor keeps a channel post's downstream channel neighbour lower than the floor beside it, F1; run s42_r0 had floor posts 1.8 m below the next T1a post)".into(),
            hillslope_rule: "REVISION 1 two passes: pass 1 as below with L_ref; then each free post's hillslope length L_H = (D8 steepest-descent flow length up to the divide, max over donors) + (flow length down to a channel, floor or water post) on the pass-1 surface, floored at one spacing and box-smoothed over 64 m (box_mean), and pass 2 repeats the eikonal with beta = 2 E L_H / (D Sc_ref). PASS-1 TEXT: roering_steady_eikonal_v1: first-order fast marching (4-neighbour upwind quadratic, heap key (z, global gy, gx)) from channels, floors and water; |grad z| = max(hill_slope_min, Sc(rock at the post, at the upwind height) x f(beta)), f(beta) = (sqrt(1 + beta^2) - 1)/beta, beta = 2 E L_ref / (D Sc_ref), D = T0 header hillslope diffusivity: the steady-state nonlinear-creep gradient at hillslope length L_ref (Roering et al. 2007 EPSL 264:245, recalled), with the rock contrast kept at every erosion rate".into(),
            crest_rule: "box_crest_v1: on posts that are not sources, z = min(z, b), b = separable box mean (x pass then y pass, each post's own half-width, indices clamped to the lattice; REVISION 1: metric half-width r with weights clamp(r/dx + 1/2 - |k|, 0, 1), continuous in r, replacing round(r / spacing) whole posts), r = min(D Sc_ref / E, crest_radius_max): the box mean of a V crest is the parabola of curvature E/D (steady-state linear diffusion), so only convex crests and noses are lowered".into(),
            fill_rule: "priority_flood_fill_v1: 8-neighbour priority flood from channels, water and the lattice ring (key (z, global gy, gx)); a post below its spill level is raised to it".into(),
            e_rule: "hat_blended_cell_rate_v1: E(x) = sum over the 4 surrounding T0 cell centres of the plain bilinear weight x E_cell (0 on lake cells and on cells without enough land, which keep E0 only if solved); E0 = e0_gain x T0 subgrid closure rate K R_sub / ln(sqrt(cell area) / x_c)".into(),
            solve_rule: "f6_hat_galerkin_v1: chart-wide lattice at solve_spacing; for each T0 cell with land hat mass >= solve_min_land_hat of its hat mass (land = zero-contour land), M = hat-weighted mean of z over land posts, T = same of the T0 bilinear, L = same of the T1a trunk level (floored at sea); Jacobi updates: damped Jacobi in ln E, step = 0.6 ln((T - L)/(M - L)) clamped to +-ln 2 (or +-ln 2 by the sign of M - T when M or T is within 1 m of L), bounds [e_min, e_max]; stop when every |M - T| <= tolerance or after solve_updates_max updates; the published table is the last evaluated one".into(),
            lithology_rule: if strata {
                "t0_strata_full_v1: per T0 cell the full column (lane SA stratafull_off/lith/th, layers >= 1 m) is placed with its top at the T0 cell elevation, tilted by the T0 uplift gradient (strata_grad_e/n) clamped to dip_clamp, anchored at the site point (chart m); the exposed code at (x, y, z) is the layer containing depth datum - z (above the top: layer 0; below the base: the last layer); empty (oceanic) columns use the stand-in basement code; Sc and k_mult from the table (pf-detail B2 rule without bedding or remap)".into()
            } else { "standin_cover_over_basement_v1 (lane T3 CoverOverBasementV1 with LithologyV1::standin_v1)".into() },
            water_rule: "t0_zero_contour_v1: marine iff the T0 bilinear <= sea level on T0 land and marine cells; lake cells at their stage".into(),
            network_mode_rule: "EXPLORATORY (beyond REVISION 1, run r2): random_walk_geodesic_v1: d = geodesic distance from the seeds over land nodes (Dijkstra, 16-neighbour stencil with Euclidean step lengths, the same T0-cell confinement, then unconfined for the rest); land nodes in increasing (d, node) choose a receiver among their 8 neighbours j with d_j < d_i that they may drain to (own cell, T0 receiver cell, water) and whose diagonal does not cross an existing link, with probability proportional to (d_i - d_j) / step, by the keyed uniform (t1b.rivernet.route); none left: the lowest-d neighbour. Paths always approach the seeds (near-geodesic lengths, so steady-state beds of neighbouring channels agree) and wander and coalesce at random (Scheidegger 1967, recalled: drainage into an intramontane trench); mode 0 = invasion percolation (r0), 1 = first-passage percolation (r1)".into(),
            rockmass_rule: "EXPLORATORY sub-resolution-closure-v1 (ADR 043 section 2): rock-mass strength (joint spacing, weathering) below the resolved rock class: every channel chain (network chain start node, or T1a node) carries a keyed multiplier m = clamp(exp(sigma x N), min, max) (t1b.rivernet.rockmass, normal_q32 table draw, counter = local node index); a hillslope post takes the m of the channel or floor post the pass-1 eikonal reached it from, and pass 2 uses Sc x m; channel beds, floors and water are not affected; sigma 0 reproduces the closure-free surface".into(),
            estuary_rule: "EXPLORATORY t1a_estuary_v1 (lane NO's rule): a channel post whose bed is at or below sea level is an estuary arm of the sea (fixed marine at sea level); the T1a mouth sits on the T0 cell-edge coast, which the zero-contour coast can leave inland".into(),
            coast: None,
            volcano: None,
            cf: None,
            fx: None,
            ls: None,
            climate_runoff: None,
            vf: None,
            hs: None,
            seam_cert: None,
            cl: None,
            glacial: None,
            es: None,
            dam_lake: None,
            edifice_foot: None,
            pedestal: None,
            confinement: None,
            ceil: None,
            outlet_checks: 0,
            wall_checks: 0,
            ls_global: None,
            bank_rule: "EXPLORATORY t1a_bank_v1: after crest rounding, a free post 8-adjacent to a T1a channel post is raised to at least the highest adjacent T1a bed + 0.01 m (the river is incised below its banks), so steepest descent from a T1a post stays in the channel (F1)".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        self.lithology.validate()?;
        if self.lithology.volcanic_split.is_some() && self.lithology_source != LITHO_SRC_STRATA { return Err("vs (volcanic split) needs the strata lithology source (strata=)".into()); }
        if self.solve_spacing_mm % 32_000 != 0 || self.margin_mm % 32_000 != 0 { return Err("rivernet spacings must be whole 32 m posts".into()); }
        if self.e_min_m_per_my <= 0 || self.e_max_m_per_my <= self.e_min_m_per_my { return Err("E bounds".into()); }
        if self.f6_reference > 2 { return Err("f6_reference: 0 bilinear, 1 T0 cell mean, 2 T0 cell mean on coastal cells".into()); }
        if let Some(c) = &self.coast { c.validate()?; if self.network_mode != 2 { return Err("coast needs network_mode 2".into()); } }
        if let Some(v) = &self.volcano { v.validate()?; }
        if let Some(c) = &self.cf { c.validate()?; if self.coast.is_none() || self.volcano.is_none() { return Err("cf needs coast=1 and volcano=1".into()); } }
        if let Some(f) = &self.fx { f.validate()?; if self.cf.is_none() { return Err("fx needs cf=1".into()); } if self.network_mode != 2 { return Err("fx needs network_mode 2".into()); } }
        if let Some(l) = &self.ls { l.validate()?; if self.fx.is_none() { return Err("ls needs fx=1".into()); } }
        if let Some(t) = crate::t1b_tf::tf(self) { t.validate()?; }
        if let Some(t) = crate::t1b_tg::tg(self) { t.validate()?; if t.sill_grade == 1 && crate::t1b_tf::tf(self).map(|f| f.bluff_v2 != 1).unwrap_or(true) { return Err("tg sill_grade needs tf bluff_v2".into()); } }
        if let Some(t) = crate::t1b_th::th(self) {
            t.validate()?;
            if t.expo_smooth == 1 && crate::t1b_tf::tf(self).map(|f| f.exposure != 1).unwrap_or(true) { return Err("th expo_smooth needs tf exposure".into()); }
            if t.uncut == 1 && self.cf.as_ref().map(|c| c.cliff_profile != 1).unwrap_or(true) { return Err("th uncut needs the cf cliff profile".into()); }
        }
        if let Some(v) = crate::t1b_ls::tc(self).and_then(|t| t.ve.as_ref()) {
            v.validate()?;
            if v.on == 1 && (self.volcano.is_none() || self.cf.is_none()) { return Err("ve needs volcano=1 and cf=1".into()); }
        }
        if let Some(c) = crate::t1b_ct::ct(self) { c.validate()?; }
        if let Some(v) = &self.vf {
            v.validate()?;
        }
        if let Some(h) = &self.hs { h.validate()?; }
        if let Some(c) = &self.seam_cert { c.validate()?; }
        if let Some(c) = &self.cl { c.validate()?; if crate::t1b_ls::tc_curve(self) { return Err("cl and ls_curve (lane TC C5 keyed midpoints) are exclusive".into()); } }
        if let Some(g) = &self.glacial { g.validate()?; if crate::t1b_ls::te_cell(self) { return Err("glacial needs the chart-wide E solve (not ls_cell)".into()); } }
        if let Some(d) = &self.dam_lake { d.validate()?; if self.volcano.is_none() { return Err("dam_lake needs volcano=1".into()); } if self.climate_runoff.is_none() { return Err("dam_lake needs the climate key (climate=<dir>): the lake water balance reads P and PET".into()); } }
        if let Some(q) = &self.pedestal { q.validate()?; if self.volcano.is_none() { return Err("ped needs volcano=1".into()); } if crate::t1b_ls::te_cell(self) { return Err("ped: the cell-local solve (te) is not wired".into()); } if q.base == 1 && self.cf.as_ref().map(|c| c.cone_hw_permille <= 0).unwrap_or(true) { return Err("ped_base needs cf=1 (the CF scaled cone)".into()); } }
        if let Some(f) = &self.edifice_foot { f.validate()?; if self.volcano.is_none() || self.cf.is_none() { return Err("ef_cap needs volcano=1 and cf=1 (the CF / VE composite cone)".into()); } }
        if let Some(c) = &self.ceil { c.validate()?; }
        crate::t1b_ls::validate_global(self)?;
        crate::t1b_outlet::validate_terrain(self)?;
        if let Some(e) = &self.es { e.validate()?; if crate::t1b_ls::te_cell(self) { return Err("es needs the chart-wide E solve (not ls_cell)".into()); } if e.update == 1 && self.ls.as_ref().map(|l| l.e_update == 1).unwrap_or(false) { return Err("es update 1 and ls_e=1 (Broyden) are exclusive".into()); } }
        Ok(())
    }
}

/// Lane MS-CEIL: the ceiling-min v1 law (MS-NEXT design section 3, unpublished-design-record).
pub const CEIL_RULE: &str = "ceiling_min_v1 (lane MS-CEIL, MS-NEXT ceiling-min v1; opt-in ceil=<c_h permille>, layer AND windows): after the crest rounding and the ct remap, before the valley-floor footslope, tributary fans, bank rule and fill, every free land post that is not a valley-floor post takes z = min(z_ss, Z_c): z_ss = the construction surface at the solved table E_v (unchanged; eikonal passes, crest rounding, ct), Z_c = B + c_h x R, B = the plain T0 bilinear at the post (the lattice bil), R = sum over the post's 4 hat cells of hat weight x the T0 subgrid relief (t0.subgrid_relief_m, the field that sets the T1a trunk level as cell mean minus R; all four cells as recorded), c_h = c_h_permille / 1000 [H: one declared constant of the stage: the un-incised upland (constructional ramp) stands c_h x the cell's subgrid relief above the T0 bilinear; prototype A/B, published]. Channels, floors and water are sources and are never ceiled. The coast / volcano stage builds the edifices afterwards over this ground (unchanged code: on the cone core the body does not read the ground; flank blending does); final edifice posts are surface_src 2. The E solve (solve_e, unchanged) sees the min in every evaluation; E stays the solved per-cell table. Two surfaces: threshold walls of the active network (z_ss) and the constructional ramp (Z_c) where the network has not cut (L2 ceiling-min)";
/// Window field notes (lane MS-CEIL).
pub const CEIL_NOTE: &str = "lane MS-CEIL ceiling_min_v1 (ceil=): Z_c = T0 bilinear + c_h x hat-blended T0 subgrid relief at the post (m, f32), the ceiling every free non-floor land post was clipped to (params.ceil, CEIL_RULE)";
/// Lane MS-CEIL-2: channels cut into the ceiling (sub-key ceil_bed=1 of ceil).
pub const CEIL_BED_RULE: &str = "ceiling_bed_grade_v1 (lane MS-CEIL-2, sub-key ceil_bed=1 of ceil; layer AND windows): in the node bed law (node_beds_core: every evaluation, LS1's final beds, the windows), a land node n (not a conduit) and a T1a node graded from its outlet take a graded bed (below): bed_r the receiver's (outlet's) bed as already graded, law bed = bed_r + S len (the 0.9.2 / TD2b stream-power rise, unchanged), Z_c(n) = the ceiling at the node (T0 bilinear + c_h x hat-blended T0 subgrid relief, as the lattice ceiling, f64), d_n = S_h(n) x L_n the valley-side relief of the node's link: L_n = the TD2b network hillslope length (link_hill_length, A_node / (2 len); len to the receiver through the meander midpoint, to the T1a receiver for T1a nodes), S_h = the TD2b valley-side gradient side_gradient(E_n, Sc of the rock at the height Z_c(n), rock-mass multiplier of n, L_n) (the steady-state hillslope law at the node's E: Sc f(beta), beta = 2 E L / (D Sc_ref)); the bed law takes bed_n = max(bed_r, min(law bed, cap(n))), cap(n) = max(Z_c(n) - d_n, present sea level) (Z_c is a subaerial constructional surface; next to marine T0 cells the T0 bilinear is pulled far below sea level and grades nothing); then a river never rises: in upstream-first order a land node still above its cap (held at its receiver where the ceiling falls upstream) is lowered to its cap, and each lowered land node lowers its land receiver to at most its own bed - floor_slope x len (floor_slope = params.floor_slope_ppm, the declared valley-floor slope 0.005) until the reach meets its own graded profile; never below the first T1a, water or conduit bed the path reaches (T1a trunks are bases here: lowering them ran low-gradient trunks down to the lowstand and the coast stage drowned them); tributaries of a lowered reach keep their beds (hanging at the junction). Reading: the ceiling is the divide of every valley cut into it, so a channel stands one steady-state hillslope relief below it (the valley cross-section of the network's own hillslope law); where the law bed is already lower the channel is unchanged (the knickpoint has passed); a bed never rises above the ceiling minus that relief and never falls below its receiver (monotone by construction, the same form as the TC C1 floor cap). LS1's routing cost is unchanged (t1b_ls.rs not edited): the routes are the key-off routes, only the beds are graded. Walls: the unchanged eikonal from the graded beds, clipped by the unchanged ceiling min";
/// Lane MS-CEIL-2: the cone's pedestal ring reads the un-ceiled ground (sub-key ceil_ped=1 of ceil).
pub const CEIL_PED_RULE: &str = "ceiling_free_pedestal_v1 (lane MS-CEIL-2, sub-key ceil_ped=1 of ceil; layer AND windows; needs ped_base to matter): the T1B-PEDESTAL ground (the heights before the coast / volcano stage, the ring the local base G is measured on) takes on every post the ceiling lowered the larger of its height there and its height before the ceiling min (z_ss), so the ceiling itself does not move the cone's base G or its footprint R_f; the published surface is unchanged (the min stays where CEIL_RULE puts it). The re-solved E table still moves z_ss and with it G (reported, not hidden)";
pub const CEIL_SRC_NOTE: &str = "lane MS-CEIL ceiling_min_v1 (ceil=): 0 active surface (z_ss <= Z_c at the hook: walls and slopes of the network), 1 ceiling (z_ss > Z_c, the post was lowered to Z_c; later stages may raise it by fill), 2 not a ceiling candidate (channel, water, valley floor, or a final volcano edifice post)";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CeilParamsV1 {
    pub rule: String,
    /// c_h in permille (the one declared constant of the stage).
    pub c_h_permille: i64,
    /// Lane MS-CEIL-2 (sub-key ceil_bed=1, layer AND windows): CEIL_BED_RULE. 0 = off (serialization unchanged).
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32")]
    pub bed: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bed_rule: Option<String>,
    /// Lane MS-CEIL-2 (sub-key ceil_ped=1, layer AND windows): CEIL_PED_RULE. 0 = off (serialization unchanged).
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32")]
    pub ped: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ped_rule: Option<String>,
}
impl CeilParamsV1 {
    pub fn v1(c_h_permille: i64) -> Self { CeilParamsV1 { rule: CEIL_RULE.into(), c_h_permille, bed: 0, bed_rule: None, ped: 0, ped_rule: None } }
    /// Lane MS-CEIL-2: switch a sub-rule on (key ceil_bed=1 / ceil_ped=1; 0 leaves it off).
    pub fn set_bed(&mut self, v: u32) { self.bed = v; self.bed_rule = if v == 1 { Some(CEIL_BED_RULE.into()) } else { None }; }
    pub fn set_ped(&mut self, v: u32) { self.ped = v; self.ped_rule = if v == 1 { Some(CEIL_PED_RULE.into()) } else { None }; }
    pub fn validate(&self) -> R<()> {
        if self.rule != CEIL_RULE { return Err("ceil: unknown rule".into()); }
        if !(0..=5000).contains(&self.c_h_permille) { return Err("ceil: c_h_permille must be in 0..=5000".into()); }
        if self.bed > 1 || self.ped > 1 { return Err("ceil_bed / ceil_ped must be 0 or 1".into()); }
        if (self.bed == 1) != (self.bed_rule.as_deref() == Some(CEIL_BED_RULE)) { return Err("ceil: unknown bed rule".into()); }
        if (self.ped == 1) != (self.ped_rule.as_deref() == Some(CEIL_PED_RULE)) { return Err("ceil: unknown ped rule".into()); }
        Ok(())
    }
    pub fn c_h(&self) -> f64 { self.c_h_permille as f64 / 1000.0 }
}

/// Lane MS-CEIL: Z_c per lattice post (f32 as published; the min uses the same value).
pub fn ceil_surface(lp: &LatPrep, t0: &T0, cp: &CeilParamsV1) -> Vec<f32> {
    let ch = cp.c_h();
    (0..lp.lat.nodes()).map(|i| {
        let mut r = 0.0f64;
        for k in 0..4 { let w = lp.hat_w[i][k]; if w > 0.0 { r += w as f64 * t0.subgrid_relief_m[lp.hat_c[i][k] as usize] as f64; } }
        (lp.bil[i] as f64 + ch * r) as f32
    }).collect()
}

/// DerivationRecipeV1 shape (ADR 044 section 1) for the `t1b.rivernet` chart layer.
#[derive(Clone, Debug, Serialize)]
pub struct RiverNetRecipeV1 {
    pub schema: String,
    pub layer: String,
    pub generator: GeneratorRef,
    pub params: RiverNetParamsV1,
    pub params_digest: String,
    pub laws: ReliefLawsV1,
    pub universe_digest: String,
    pub upstream: Vec<(String, String)>,
    pub seed_scope: SeedScopeV1,
    pub chart: ChartDeclV1,
    pub key_domain: String,
    pub label: String,
    /// Lane CLIMATE-RUNOFF (key on only; not serialized, not digested): the digest of the climate-free twin of the T1a
    /// skeleton (its header's climate_twin_digest), which replaces the T1a upstream in the draw digest.
    #[serde(skip)]
    pub t1a_draw_twin: Option<String>,
}

impl RiverNetRecipeV1 {
    pub fn digest(&self) -> R<String> { pf_contracts::artifact::digest_canonical(self).map_err(|e| format!("recipe digest: {e}")) }
    /// The digest keyed draws use (lane TB tb_paired_keying_v1): the recipe without params.ls.tb (equal to `digest` when absent).
    /// Lane CLIMATE-RUNOFF (climate_paired_keying_v1, key on only): also without params.climate_runoff and the t0.climate
    /// upstream, and with the climate-free T1a twin's digest as the t1a.skeleton upstream, so the network, meander, rock-mass
    /// and volcano draws are the same with and without climate (climate changes physics, not draws).
    pub fn draw_digest(&self) -> R<String> {
        if self.upstream.iter().any(|(r, _)| r.starts_with("physical.")) {
            let mut c = self.clone();
            // laws.source is provenance prose, not a physical law. Its legacy
            // serialization is a fixed v1 salt for compatibility; only the
            // projected physical identity enters it. Editing the actual source
            // description or its artifact digest must not re-key the world.
            if let Some((_, physical)) = self.upstream.iter().find(|(r, _)| r == "physical.t0.terrain") {
                c.laws.source = format!("T0 terrain artifact {} header domain.erosion and assumption A-SUBGRID-CHANNEL-HEAD (legacy project assumptions A-K-UNITS, A-HILLSLOPE-D; awaiting technical review Q2/Q3)", physical);
            }
            crate::skeleton::project_physical_upstream(&mut c.upstream);
            return c.draw_digest();
        }
        // Lane SPLIT-T1B: its volcanic split keys (vs_ct_bed, vs_ct_blend, vs_aa, vs_band_mm, vs_aa_dx_mm) change physics, not
        // draws; only those fields are stripped, so every vs=1 output without them keeps its draws.
        if self.params.lithology.volcanic_split.as_ref().map(|v| v.t1b_keys()).unwrap_or(false) {
            let mut c = self.clone();
            c.params.lithology.volcanic_split = c.params.lithology.volcanic_split.as_ref().map(|v| v.draw_twin());
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane LITHO-CONTINUITY: the rock-body continuity rule changes the rock cut (physics), not draws.
        if self.params.lithology.continuity.is_some() {
            let mut c = self.clone();
            c.params.lithology.continuity = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane T0-IMPRINT-PROD: the confinement rule changes which links are allowed, not draws.
        if self.params.confinement.is_some() {
            let mut c = self.clone();
            c.params.confinement = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane T1B-PEDESTAL: the pedestal params change the E solve / edifice base, not draws.
        if self.params.pedestal.is_some() {
            let mut c = self.clone();
            c.params.pedestal = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane EDIFICE-FOOT: the edifice_foot params change the window edifice surface, not draws.
        if self.params.edifice_foot.is_some() {
            let mut c = self.clone();
            c.params.edifice_foot = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane DAMMED-VALLEY: the dam_lake params change the window water mask, not draws.
        if self.params.dam_lake.is_some() {
            let mut c = self.clone();
            c.params.dam_lake = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane E-SOLVE-GENERALITY: the es params change the E solve, not draws.
        // Lane MS-CEIL: the ceil params lower hillslope posts, not draws (key-on shares the key-off network).
        if self.params.ceil.is_some() {
            let mut c = self.clone();
            c.params.ceil = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane LS1-GLOBAL: the ls_global params change LS1's order and stop, not draws (key-on shares the key-off network).
        if self.params.ls_global.is_some() {
            let mut c = self.clone();
            c.params.ls_global = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        if self.params.es.is_some() {
            let mut c = self.clone();
            c.params.es = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane GLACIAL-T1: the glacial params change beds and the carve, not draws (key-on shares the key-off network).
        if self.params.glacial.is_some() {
            // Lane GLACIER-1: present ice only changes no rock, so its draws are keyed exactly as key off, generator version
            // included (the GLACIAL-T1 version string otherwise reaches the draw digest and re-draws the network).
            let present = self.params.glacial.as_ref().map(|g| g.present_ice_only == 1).unwrap_or(false);
            let mut c = self.clone();
            c.params.glacial = None;
            if present { c.generator.version = c.params.generator_version().into(); }
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane CHANNEL-CENTRELINE: the cl params change the drawn centreline (window only), not draws.
        if self.params.cl.is_some() {
            let mut c = self.clone();
            c.params.cl = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane SEAM-CERT: the seam_cert params change the window fills and routing, not draws.
        if self.params.seam_cert.is_some() {
            let mut c = self.clone();
            c.params.seam_cert = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane HS-PASS: the hs params change physics (window only), not draws.
        if self.params.hs.is_some() {
            let mut c = self.clone();
            c.params.hs = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        // Lane VALLEY-FLOOR-FIX: the vf params change physics, not draws.
        if self.params.vf.is_some() {
            let mut c = self.clone();
            c.params.vf = None;
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            return c.draw_digest();
        }
        if self.params.climate_runoff.is_some() {
            let mut c = self.clone();
            c.params.climate_runoff = None;
            c.upstream.retain(|(r, _)| r != "t0.climate");
            let twin = self.t1a_draw_twin.clone().ok_or("climate rivernet recipe without the T1a draw twin digest")?;
            for u in c.upstream.iter_mut() { if u.0 == "t1a.skeleton" { u.1 = twin.clone(); } }
            c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
            c.t1a_draw_twin = None;
            return c.draw_digest();
        }
        if crate::t1b_ls::tb(&self.params).is_none() { return self.digest(); }
        let mut c = self.clone();
        if let Some(l) = c.params.ls.as_mut() { l.tb = None; }
        c.params_digest = pf_contracts::artifact::digest_canonical(&c.params).map_err(|e| e.to_string())?;
        c.digest()
    }
}

// ------------------------------------------------------------------ lithology sources

/// The T0 full stratigraphy columns (lane SA dumps) the strata rule reads.
pub struct StrataCols {
    pub off: Vec<u32>,
    pub lith: Vec<u8>,
    pub th: Vec<f32>,
    pub grad_e: Vec<f32>,
    pub grad_n: Vec<f32>,
    pub digests: Vec<(String, String)>,
    /// cell -> (layer base depths below the top, codes, index in the full column), layers >= 1 m only.
    /// Lane PF-SERVE: a per-cell OnceLock memo (was RefCell<BTreeMap>) so one chart context can be shared by window threads;
    /// each entry is the same pure function of the cell, so values are unchanged.
    cols: std::sync::OnceLock<Vec<std::sync::OnceLock<(Vec<f64>, Vec<u8>, Vec<u16>)>>>,
    /// Lane LITHO-CONTINUITY (lcont=1 only): the planet-wide rock-body identity, built once per process on first use.
    pub lcont: std::sync::OnceLock<crate::t1b_lcont::Identity>,
}

impl StrataCols {
    pub fn load(final_dir: &Path) -> R<StrataCols> {
        let mut digests = Vec::new();
        let mut rd = |role: &str, f: &str| -> R<Vec<u8>> {
            let b = std::fs::read(final_dir.join(f)).map_err(|e| format!("{f}: {e}"))?;
            digests.push((role.to_string(), format!("blake3:{}", blake3::hash(&b).to_hex())));
            Ok(b)
        };
        let u32s = |b: Vec<u8>| -> Vec<u32> { b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() };
        let f32s = |b: Vec<u8>| -> Vec<f32> { b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect() };
        let off = u32s(rd("t0.stratigraphy.full.offset", "stratafull_off.u32")?);
        let lith = rd("t0.stratigraphy.full.lith", "stratafull_lith.u8")?;
        let th = f32s(rd("t0.stratigraphy.full.thickness", "stratafull_th.f32")?);
        let grad_e = f32s(rd("t0.stratigraphy.uplift_grad_east", "strata_grad_e.f32")?);
        let grad_n = f32s(rd("t0.stratigraphy.uplift_grad_north", "strata_grad_n.f32")?);
        if lith.len() != th.len() || off.len() != grad_e.len() + 1 || grad_n.len() != grad_e.len() { return Err("strata field sizes".into()); }
        digests.sort();
        Ok(StrataCols { off, lith, th, grad_e, grad_n, digests, cols: Default::default(), lcont: Default::default() })
    }
    /// Lane LITHO-CONTINUITY: f(kept column of cell c) (the cut's layers >= 1 m: depth bases, codes, full indices).
    pub fn kept<T>(&self, c: u32, f: impl FnOnce(&(Vec<f64>, Vec<u8>, Vec<u16>)) -> T) -> T { self.with_col(c, f) }
    fn with_col<T>(&self, c: u32, f: impl FnOnce(&(Vec<f64>, Vec<u8>, Vec<u16>)) -> T) -> T {
        let memo = self.cols.get_or_init(|| (0..self.off.len().saturating_sub(1)).map(|_| std::sync::OnceLock::new()).collect());
        f(memo[c as usize].get_or_init(|| {
            let (a, b) = (self.off[c as usize] as usize, self.off[c as usize + 1] as usize);
            let (mut base, mut code, mut idx) = (Vec::new(), Vec::new(), Vec::new());
            let mut d = 0.0;
            for k in a..b {
                let t = self.th[k] as f64;
                if t < 1.0 { continue; }
                d += t;
                base.push(d);
                code.push(self.lith[k]);
                idx.push((k - a).min(u16::MAX as usize) as u16);
            }
            (base, code, idx)
        }))
    }
}

/// Per-post rock: either the stand-in or a strata cut.
pub enum RockPosts<'a> {
    Standin(crate::t1b_litho::CoverOverBasementV1),
    Strata { table: &'a LithologyV1, cols: &'a StrataCols, cell: Vec<u32>, datum: Vec<f64>, fallback: Vec<u8>, sc: [f64; 12], k: [f64; 12],
        /// Lane VOLCANIC-SPLIT (opt-in vs=1): lava / volcaniclastic packages inside lithology-3 layers; None = unchanged.
        volc: Option<crate::t1b_volc::VolcPosts>,
        /// Lane LITHO-CONTINUITY (opt-in lcont=1): continuous rock bodies across T0 cells (carries the volcanic split itself);
        /// None = unchanged.
        cont: Option<Box<crate::t1b_lcont::ContPosts>> },
}

impl<'a> RockPosts<'a> {
    /// Build for posts at (lat, lon, x_m, y_m).
    pub fn build(p: &'a RiverNetParamsV1, t0: &T0, extras: &T0Extras, strata: Option<&'a StrataCols>, anchor_m: (f64, f64), pos: &[(f64, f64, f64, f64)]) -> R<RockPosts<'a>> {
        if let Some(env) = &p.lithology.envelope {
            return crate::t1b_lcont::envelope_rock(env, &p.lithology, t0, strata.ok_or("env needs strata=")?, pos);
        }
        match (p.lithology_source.as_str(), strata) {
            (LITHO_SRC_STANDIN, _) => {
                if p.lithology.continuity.is_some() { return Err("lcont (rock-body continuity) needs the strata lithology source (strata=)".into()); }
                let s = crate::t1b_litho::build_standin(&p.lithology, t0, &extras.rock(), pos.len(), &|i| Ok((pos[i].0, pos[i].1)))?;
                Ok(RockPosts::Standin(s))
            }
            (LITHO_SRC_STRATA, Some(cols)) => {
                let clamp = p.dip_clamp_ppm as f64 / 1e6;
                let mut cell = Vec::with_capacity(pos.len());
                let mut datum = Vec::with_capacity(pos.len());
                let mut fallback = Vec::with_capacity(pos.len());
                for &(la, lo, x, y) in pos {
                    let c = t0.cell_of_deg(la, lo);
                    let (ge, gn) = (cols.grad_e[c as usize] as f64, cols.grad_n[c as usize] as f64);
                    let m = ge.hypot(gn);
                    let s = if m > clamp { clamp / m } else { 1.0 };
                    cell.push(c);
                    datum.push(t0.elevation_m[c as usize] as f64 + ge * s * (x - anchor_m.0) + gn * s * (y - anchor_m.1));
                    fallback.push(basement_code(extras.rock_family[c as usize], extras.crust_kind[c as usize]));
                }
                let mut sc = [0.0; 12];
                let mut k = [0.0; 12];
                for c in 0..12u8 { sc[c as usize] = p.lithology.sc(c); k[c as usize] = p.lithology.k(c); }
                // Lane LITHO-CONTINUITY: with lcont the bodies (and their volcanic split) come from t1b_lcont; datum per post from it.
                if let Some(lc) = p.lithology.continuity.as_ref() {
                    let id = cols.lcont.get_or_init(|| crate::t1b_lcont::Identity::build(cols, t0));
                    let cp = crate::t1b_lcont::ContPosts::build(lc, id, cols, t0, &extras.rock_family, clamp, p.lithology.volcanic_split.as_ref(), pos);
                    println!("  lane LITHO-CONTINUITY posts: {}", serde_json::to_string(&cp.stats).unwrap_or_default());
                    let datum = cp.datum.clone();
                    return Ok(RockPosts::Strata { table: &p.lithology, cols, cell, datum, fallback, sc, k, volc: None, cont: Some(Box::new(cp)) });
                }
                let volc = p.lithology.volcanic_split.as_ref().map(|v| crate::t1b_volc::VolcPosts::build(v, &cell, pos.iter().map(|q| (q.2, q.3)).collect(), &extras.rock_family,
                    &|c| cols.with_col(c, |(b, k, ix)| if k.is_empty() { None } else { Some((b.clone(), k.clone(), ix.clone())) })));
                Ok(RockPosts::Strata { table: &p.lithology, cols, cell, datum, fallback, sc, k, volc, cont: None })
            }
            _ => Err("lithology source needs its inputs (strata dir for t0_strata_full_v1)".into()),
        }
    }
}

impl<'a> ExposedLithology for RockPosts<'a> {
    #[inline]
    fn at(&self, i: usize, z: f64) -> (f64, f64, u8) {
        match self {
            RockPosts::Standin(s) => s.at(i, z),
            RockPosts::Strata { cont: Some(c), fallback, sc, k, .. } => c.at(i, z, sc, k, fallback[i]),
            RockPosts::Strata { cols, cell, datum, fallback, sc, k, volc, .. } => {
                let d = datum[i] - z;
                if let Some(r) = volc.as_ref().and_then(|v| v.at(i, d, &|c| { let c = (c as usize).min(11); (sc[c], k[c]) })) { return (r.0, r.1, r.2); }
                let code = cols.with_col(cell[i], |(base, code, _)| {
                    if code.is_empty() { return fallback[i]; }
                    if d <= 0.0 { return code[0]; }
                    let j = match base.binary_search_by(|b| b.partial_cmp(&d).unwrap()) { Ok(j) => (j + 1).min(code.len() - 1), Err(j) => j.min(code.len() - 1) };
                    code[j]
                });
                let c = (code as usize).min(11);
                (sc[c], k[c], code)
            }
        }
    }
}

impl<'a> RockPosts<'a> {
    /// Lane PV: (T0 cell, index in the cell's FULL column, top = 0) of the layer exposed at post i, surface z, by the same
    /// cut as `at`; None for the stand-in (no stack) and for an empty (oceanic) column.
    pub fn layer_at(&self, i: usize, z: f64) -> Option<(u32, u16)> {
        match self {
            RockPosts::Standin(_) => None,
            RockPosts::Strata { cont: Some(c), datum, .. } => c.unit_at(i, datum[i] - z).and_then(|u| c.src(i, u.k)),
            RockPosts::Strata { cols, cell, datum, volc, .. } => {
                let d = datum[i] - z;
                if let Some(c) = volc.as_ref().and_then(|v| v.col(i)) { return Some((cell[i], c.idx[c.unit_at(d)])); }
                cols.with_col(cell[i], |(base, code, idx)| {
                    if code.is_empty() { return None; }
                    if d <= 0.0 { return Some((cell[i], idx[0])); }
                    let j = match base.binary_search_by(|b| b.partial_cmp(&d).unwrap()) { Ok(j) => (j + 1).min(code.len() - 1), Err(j) => j.min(code.len() - 1) };
                    Some((cell[i], idx[j]))
                })
            }
        }
    }
    /// Recorded exposed source-layer contacts (absolute metres), not procedural beds.
    /// The first ordinary strata top is unrecorded, as in pf-detail::Ctx::unit_depth.
    /// No last-layer extension supplies a finite contact below the retained column.
    pub fn unit_contacts(&self, i: usize, z: f64, envelope: bool) -> Option<(u16, f64, Option<f64>)> {
        match self {
            Self::Standin(_) => None,
            Self::Strata { cont: Some(c), datum, .. } => {
                let u = c.unit_at(i, datum[i] - z)?;
                if datum[i] - z >= u.base { return None; }
                let (_, layer) = c.src(i, u.k)?;
                Some((layer, datum[i] - u.base, if u.first && !envelope { None } else { Some(datum[i] - u.top) }))
            }
            Self::Strata { cols, cell, datum, .. } => cols.with_col(cell[i], |(base, code, index)| {
                if code.is_empty() { return None; }
                let depth = datum[i] - z;
                let k = base.iter().position(|&b| depth < b)?;
                Some((index[k], datum[i] - base[k], if k == 0 { None } else { Some(datum[i] - base[k - 1]) }))
            }),
        }
    }
    /// Lane STEP2: f(cell, datum, unit depth bases, codes) of post i's column (strata only; None for the stand-in and an
    /// empty column).
    pub fn column<T>(&self, i: usize, f: impl FnOnce(u32, f64, &[f64], &[u8]) -> T) -> Option<T> {
        match self {
            RockPosts::Standin(_) => None,
            RockPosts::Strata { cont: Some(c), cell, datum, .. } => { let (b, k, _) = c.expanded(i); if k.is_empty() { None } else { Some(f(cell[i], datum[i], &b, &k)) } }
            RockPosts::Strata { cols, cell, datum, volc, .. } => {
                if let Some(c) = volc.as_ref().and_then(|v| v.col(i)) { return Some(f(cell[i], datum[i], &c.base, &c.code)); }
                cols.with_col(cell[i], |(base, code, _)| if code.is_empty() { None } else { Some(f(cell[i], datum[i], base, code)) })
            }
        }
    }
    /// Lane BENCH-VAR: as `column`, plus each unit's index in the cell's FULL column (top = 0) and the full column length.
    pub fn column_ix<T>(&self, i: usize, f: impl FnOnce(u32, f64, &[f64], &[u8], &[u16], usize) -> T) -> Option<T> {
        match self {
            RockPosts::Standin(_) => None,
            RockPosts::Strata { cont: Some(c), cols, cell, datum, .. } => {
                // Each unit's index counted from the bottom of ITS recording layer's column (the bed key), mapped onto the own cell's length.
                let (b, k, s) = c.expanded(i);
                if k.is_empty() { return None; }
                let full = (cols.off[cell[i] as usize + 1] - cols.off[cell[i] as usize]) as usize;
                let idx: Vec<u16> = s.iter().map(|&(sc_, kx)| {
                    let fl = (cols.off[sc_ as usize + 1] - cols.off[sc_ as usize]) as usize;
                    let fb = fl.saturating_sub(1 + kx as usize);
                    full.saturating_sub(1 + fb).min(u16::MAX as usize) as u16
                }).collect();
                Some(f(cell[i], datum[i], &b, &k, &idx, full))
            }
            RockPosts::Strata { cols, cell, datum, volc, .. } => {
                let full = (cols.off[cell[i] as usize + 1] - cols.off[cell[i] as usize]) as usize;
                if let Some(c) = volc.as_ref().and_then(|v| v.col(i)) { return Some(f(cell[i], datum[i], &c.base, &c.code, &c.idx, full)); }
                cols.with_col(cell[i], |(base, code, idx)| if code.is_empty() { None } else { Some(f(cell[i], datum[i], base, code, idx, full)) })
            }
        }
    }
    /// Lane TREAD2: the bed-top gradient (east, north; m/m) of post i's column: the datum tilt of `build` (strata only).
    pub fn bed_grad(&self, i: usize, clamp: f64) -> Option<(f64, f64)> {
        match self {
            RockPosts::Strata { cont: Some(c), .. } => {
                let (ge, gn) = c.grad[i];
                let m = ge.hypot(gn);
                let s = if m > clamp { clamp / m } else { 1.0 };
                Some((ge * s, gn * s))
            }
            RockPosts::Strata { cols, cell, .. } => {
                let c = cell[i] as usize;
                let (ge, gn) = (cols.grad_e[c] as f64, cols.grad_n[c] as f64);
                let m = ge.hypot(gn);
                let s = if m > clamp { clamp / m } else { 1.0 };
                Some((ge * s, gn * s))
            }
            _ => None,
        }
    }
    /// Lane VOLCANIC-SPLIT: facies of the package exposed at post i, surface z (t1b_volc FAC_*; 0 when the split is off,
    /// the cell is not split or the unit is not lithology 3).
    pub fn volc_facies(&self, i: usize, z: f64) -> u8 {
        match self {
            RockPosts::Strata { cont: Some(c), datum, .. } => c.volc_at(i, datum[i] - z).map(|q| q.0).unwrap_or(crate::t1b_volc::FAC_NONE),
            RockPosts::Strata { datum, volc: Some(v), .. } => v.facies(i, datum[i] - z),
            _ => crate::t1b_volc::FAC_NONE,
        }
    }
    /// Lane VOLCANIC-SPLIT volcanic_ct_v2 at post i, height z: None = the ct class rule applies; Some(None) = a split post
    /// without a riser period (no staircase); Some(Some((ladder [s_dn, s0, s1, s_up], eta))). And the split's params.
    pub fn volc_ct(&self, i: usize, z: f64) -> Option<Option<([f64; 4], f64)>> {
        match self { RockPosts::Strata { cont: Some(c), .. } => c.ct_ladder(i, z), RockPosts::Strata { datum, volc: Some(v), .. } => v.ct_ladder(i, datum[i], z), _ => None }
    }
    /// Lane SPLIT-T1B volcanic_ct_blend_v1: the weighted ladder alternatives at post i, height z (None as `volc_ct`).
    pub fn volc_ct_alts(&self, i: usize, z: f64) -> Option<Vec<crate::t1b_volc::LadderAlt>> {
        match self { RockPosts::Strata { cont: Some(c), .. } => c.ct_ladder_alts(i, z), RockPosts::Strata { datum, volc: Some(v), .. } => v.ct_ladder_alts(i, datum[i], z), _ => None }
    }
    pub fn volc_rule(&self) -> Option<&crate::t1b_volc::VolcanicSplitV1> { match self { RockPosts::Strata { cont: Some(c), .. } => c.vrule.as_ref(), RockPosts::Strata { volc: Some(v), .. } => Some(&v.rule), _ => None } }
    /// Lane LITHO-CONTINUITY: true when columns are per post (lcont), so callers must not memoize them per T0 cell.
    pub fn per_post(&self) -> bool { matches!(self, RockPosts::Strata { cont: Some(_), .. }) }
    /// Lane PV: the T0 cell of post i (strata only).
    pub fn cell_of(&self, i: usize) -> Option<u32> { match self { RockPosts::Strata { cell, .. } => Some(cell[i]), _ => None } }
}

/// Lane PV (opt-in `prov=1`): the provenance sources of a window run.
pub struct RnProv {
    pub strata: crate::strata_prov::StrataProv,
    /// T0 volcano records by id (onset, last eruption, active) and their digest (None: no volcano records given).
    pub volc: Option<(BTreeMap<String, (Option<f64>, Option<f64>, Option<bool>)>, String)>,
}

// ------------------------------------------------------------------ the 512 m network

#[derive(Clone)]
pub struct Network {
    pub kind: Vec<u8>,
    pub recv: Vec<u32>,
    /// Seeds and water (index order), then land nodes in invasion order (receivers before donors).
    pub order: Vec<u32>,
    pub area_m2: Vec<f64>,
    /// Chart metres (jittered for land nodes).
    pub pos: Vec<[f64; 2]>,
    /// Meander midpoint of the link node -> receiver (chart m).
    pub mid: Vec<[f64; 2]>,
    /// The seed (T1a or water node) a land node's path ends on.
    pub term: Vec<u32>,
    pub level: Vec<f64>,
    pub unreached_confined: u32,
    pub nodes_land: u32,
    /// Rock-mass multiplier per node (1 without the closure).
    pub rockmass: Vec<f64>,
    /// Lane CO (coast): conduit nodes continuing T1a terminal seeds to the construction water, their bed cap (the terminal
    /// level), and per terminal its path (terminal node first, then conduits to the water node).
    pub conduit: Vec<bool>,
    pub conduit_cap: Vec<f64>,
    pub conduit_paths: Vec<Vec<u32>>,
    /// Lane CF (R1d): the T1a receiver of every node (copy of the skeleton's), for the shelf grade's final-reach test.
    pub t1a_recv: Vec<u32>,
    /// Lane TD (capture routing): each node's T0 cell (the skeleton's) and its keyed meander uniform (internal, not published).
    pub t0_cell: Vec<u32>,
    pub meander_u: Vec<f64>,
    /// Lane TC (C1 floor 2): T1a valley-axis nodes with channel == 1 that are T1a seeds here (internal, not published).
    pub t1a_ch: Vec<bool>,
    /// Lane T0-IMPRINT-PROD conf=warp only: each node's confinement cell (t1b_outlet::warp_cells); empty otherwise.
    pub conf_cell: Vec<u32>,
    /// Lane DRAIN-FREE flux mode only (params.confinement.flux = 1; None otherwise): the T1a skeleton seed areas and the
    /// conf=free bound; the seeds' `area_m2` is recomputed from the T1b catchments whenever the network is rebuilt.
    pub flux: Option<crate::t1b_outlet::FluxBase>,
}

/// Dijkstra key (path cost, node, from), min-first.
#[derive(Clone, Copy)]
struct FpKey(f64, u32, u32);
impl PartialEq for FpKey { fn eq(&self, o: &Self) -> bool { self.cmp(o) == Ordering::Equal } }
impl Eq for FpKey {}
impl PartialOrd for FpKey { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for FpKey { fn cmp(&self, o: &Self) -> Ordering { o.0.total_cmp(&self.0).then(o.1.cmp(&self.1)).then(o.2.cmp(&self.2)) } }

fn tile_scope(cache: &mut BTreeMap<(i64, i64), Scope>, layer: &Scope, tx: i64, ty: i64) -> R<Scope> {
    if let Some(s) = cache.get(&(tx, ty)) { return Ok(s.clone()); }
    let s = layer.tile(8, tx, ty).map_err(|e| e.to_string())?;
    cache.insert((tx, ty), s.clone());
    Ok(s)
}

/// Draws for node n: (ip weight, jitter x, jitter y, meander) as u64 / uniforms.
fn node_draws(layer: &Scope, cache: &mut BTreeMap<(i64, i64), Scope>, n: usize, attempt: u32) -> R<(u64, f64, f64, f64)> {
    let (gx, gy) = node_g(n);
    let (tx, ty) = (gx.div_euclid(128), gy.div_euclid(128));
    let local = ((gy - ty * 128) * 128 + (gx - tx * 128)) as u64;
    let s = tile_scope(cache, layer, tx, ty)?;
    let w = s.draw(PURPOSE_IP, local, attempt).map_err(|e| e.to_string())?.u64();
    let _ = ();
    let jx = s.draw(PURPOSE_JITTER, 2 * local, attempt).map_err(|e| e.to_string())?.uniform_f64();
    let jy = s.draw(PURPOSE_JITTER, 2 * local + 1, attempt).map_err(|e| e.to_string())?.uniform_f64();
    let m = s.draw(PURPOSE_MEANDER, local, attempt).map_err(|e| e.to_string())?.uniform_f64();
    Ok((w, jx, jy, m))
}

pub fn rivernet_recipe(t0: &T0, extras: &T0Extras, site: &SiteRecord, skel: &SkeletonView, strata: Option<&StrataCols>, p: &RiverNetParamsV1) -> R<RiverNetRecipeV1> {
    p.validate()?;
    let mut upstream = vec![
        ("site".to_string(), site.bytes_digest.clone()),
        ("t0.hydrology".to_string(), t0.hydrology_digest.clone()),
        ("t0.subgrid_relief_sidecar".to_string(), t0.subgrid_digest.clone()),
        ("t0.terrain".to_string(), t0.terrain_digest.clone()),
        ("t1a.skeleton".to_string(), skel.digest.clone()),
    ];
    let physical_keying = site.physical_keying_enabled();
    if physical_keying != skel.physical_keying {
        return Err("site and T1a keying modes differ; rebuild the skeleton for this site record".into());
    }
    if physical_keying {
        upstream.extend(crate::skeleton::physical_upstream(t0, site)?);
        upstream.push(("physical.t1a.skeleton".into(), skel.physical_digest.clone()));
    }
    if p.lithology_source == LITHO_SRC_STRATA {
        for (r, d) in &strata.ok_or("strata rule needs the strata columns")?.digests { upstream.push((r.clone(), d.clone())); }
    }
    // Lane CLIMATE-RUNOFF: T0 climate, rivernet params and the T1a's climate must agree; no mixing of climate and uniform tiers.
    let mut t1a_draw_twin = None;
    match (&p.climate_runoff, &skel.climate_digest) {
        (None, None) => {}
        (None, Some(_)) => return Err("the T1a skeleton carries climate runoff; build T1b with the same climate=<dir> key".into()),
        (Some(_), None) => return Err("climate=<dir> needs a T1a skeleton built with the same climate key".into()),
        (Some(cp), Some(sd)) => {
            let c = t0.climate.as_ref().ok_or("climate_runoff params set but no climate attached to T0")?;
            if &c.params != cp { return Err("climate_runoff params differ from the attached climate's".into()); }
            if &c.digest != sd { return Err(format!("T1a climate {sd} differs from this climate {}; refusing", c.digest)); }
            if skel.climate_params.as_ref() != Some(cp) { return Err("T1a climate_runoff params differ from these; refusing".into()); }
            upstream.push(("t0.climate".to_string(), c.digest.clone()));
            t1a_draw_twin = Some(if physical_keying {
                skel.physical_climate_twin_digest.clone().ok_or("physical-key climate T1a without physical_climate_twin_digest; rebuild skeleton")?
            } else {
                skel.climate_twin_digest.clone().ok_or("climate T1a without climate_twin_digest")?
            });
        }
    }
    upstream.sort();
    Ok(RiverNetRecipeV1 {
        schema: RN_RECIPE_SCHEMA.into(), layer: RN_LAYER.into(),
        generator: GeneratorRef { name: RN_GENERATOR_NAME.into(), version: p.generator_version().into(), tier: "Preview".into() },
        params: p.clone(), params_digest: pf_contracts::artifact::digest_canonical(p).map_err(|e| e.to_string())?,
        laws: extras.laws.clone(), universe_digest: site.universe_digest.clone(), upstream,
        seed_scope: SeedScopeV1 { keyed_rng: vz_keyed_rng::KEYED_RNG_V1.into(), seed: site.seed,
            ladder: vec!["layer = this recipe's digest".into(), "chart = 0".into(), "tile = (8, ix, iy), the lod-8 tile owning the node (half-open floor)".into(), "attempt = realization".into()],
            // Lane STEP2's purpose t1b.rivernet.ct.step is declared in params.ls.tb.tc.ct.v2.var_rule (not here: this list is
            // part of the draw digest, and a ct sub-key must not re-key the network, curve and volcano draws).
            purposes: vec![PURPOSE_IP.into(), PURPOSE_JITTER.into(), PURPOSE_MEANDER.into(), PURPOSE_ROCKMASS.into()] },
        chart: ChartDeclV1::new(site.chart_origin_lat_e7, site.chart_origin_lon_e7, (t0.radius_m * 1000.0).round() as i64),
        key_domain: format!("{RN_LAYER}/c0/8.*"),
        label: RN_LABEL.into(),
        t1a_draw_twin,
    })
}

/// Build the chart-wide 512 m network (no beds; beds depend on the E table).
pub fn build_network(t0: &T0, skel: &SkeletonView, chart: &Chart, seed: u64, recipe_digest: &str, p: &RiverNetParamsV1) -> R<Network> {
    let sea = t0.sea_level_m;
    // Lane CO: the exposed shelf (bilinear in (sea - shelf, sea]) is land in the construction.
    let shelf = p.coast.as_ref().map(|c| c.shelf_m()).unwrap_or(0.0);
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let mut kind = vec![K_LAND; NODES];
    let mut level = vec![f64::NAN; NODES];
    let mut bil_c = vec![f64::NAN; if crate::t1b_ls::tc_clip(p) { NODES } else { 0 }];
    let mut clipped_count = 0u32;
    // Lane TG G3: construction water from the box-averaged bilinear.
    let sb = match crate::t1b_tg::tg(p).filter(|t| t.smooth_bil == 1) { Some(t) => Some(crate::t1b_tg::chart_smooth(t0, chart, t.smooth_half_mm, t.smooth_grid_mm)?), None => None };
    for n in 0..NODES {
        let (gx, gy) = node_g(n);
        let (la, lo) = chart.inverse_mm(gx * SPACING_MM, gy * SPACING_MM).map_err(|e| e.to_string())?;
        let (la, lo) = (la * RAD2DEG, lo * RAD2DEG);
        let c = t0.cell_of_deg(la, lo);
        let bil = bilinear(t0, &elev, la, lo, false);
        let bil_w = match &sb { Some(s) => s.at_mm((gx * SPACING_MM) as f64, (gy * SPACING_MM) as f64), None => bil };
        if t0.water(c) == CellWater::Lake { kind[n] = K_LAKE; level[n] = t0.level(c); }
        else if bil_w <= sea - shelf { kind[n] = K_MARINE; level[n] = sea; }
        else if skel.kind[n] == 1 || skel.kind[n] == 3 { kind[n] = K_T1A; level[n] = skel.trunk_level_m[n] as f64; }
        if crate::t1b_ls::tc_clip(p) { bil_c[n] = bil; }
    }
    // Lane TC C2 (t1a_shelf_head_clip_v1): the drowned head prefix of T1a axes is ordinary construction land.
    if crate::t1b_ls::tc_clip(p) {
        let cand = |t: usize, kind: &[u8]| kind[t] == K_T1A && bil_c[t] <= sea;
        let mut pending = vec![0u32; NODES];
        for d in 0..NODES {
            if kind[d] != K_T1A { continue; }
            let r = skel.receiver[d];
            if r != NONE && r != OFF_CHART && kind[r as usize] == K_T1A { pending[r as usize] += 1; }
        }
        let mut stack: Vec<usize> = (0..NODES).filter(|&t| cand(t, &kind) && pending[t] == 0).collect();
        let mut clipped = Vec::new();
        while let Some(t) = stack.pop() {
            clipped.push(t);
            let r = skel.receiver[t];
            if r == NONE || r == OFF_CHART { continue; }
            let r = r as usize;
            if kind[r] != K_T1A { continue; }
            pending[r] -= 1;
            if pending[r] == 0 && cand(r, &kind) { stack.push(r); }
        }
        clipped.sort();
        for &t in &clipped { kind[t] = K_LAND; level[t] = f64::NAN; }
        clipped_count = clipped.len() as u32;
    }
    let layer = Scope::root(seed).layer(recipe_digest).map_err(|e| e.to_string())?.chart(0).map_err(|e| e.to_string())?;
    let mut cache = BTreeMap::new();
    let mut weight = vec![0u64; NODES];
    let mut jit = vec![(0.0, 0.0, 0.0); NODES];
    let mut rockmass = vec![1.0f64; NODES];
    let (sig, mlo, mhi) = (p.rockmass_sigma_permille as f64 / 1000.0, p.rockmass_min_permille as f64 / 1000.0, p.rockmass_max_permille as f64 / 1000.0);
    for n in 0..NODES {
        if kind[n] == K_MARINE || kind[n] == K_LAKE { continue; }
        if sig > 0.0 {
            let (gx, gy) = node_g(n);
            let (tx, ty) = (gx.div_euclid(128), gy.div_euclid(128));
            let local = ((gy - ty * 128) * 128 + (gx - tx * 128)) as u64;
            let sc = tile_scope(&mut cache, &layer, tx, ty)?;
            let g = sc.draw(PURPOSE_ROCKMASS, local, p.realization).map_err(|e| e.to_string())?.normal_f64();
            rockmass[n] = (sig * g).exp().clamp(mlo, mhi);
        }
        if kind[n] != K_LAND {
            // Lane TC C5: the meander uniform a land node at this index would get, for T1a links.
            if kind[n] == K_T1A && (crate::t1b_ls::tc_curve(p) || crate::t1b_cl::cl(p).map(|c| c.meander_on()).unwrap_or(false)) { jit[n].2 = node_draws(&layer, &mut cache, n, p.realization)?.3; }
            continue;
        }
        let (w, jx, jy, m) = node_draws(&layer, &mut cache, n, p.realization)?;
        weight[n] = w;
        jit[n] = (jx, jy, m);
    }
    let cell = |n: usize| skel.t0_cell[n];
    // Lane CO conduits (coast only): T1a terminal seeds continued down the geodesic distance to construction water.
    let mut conduit = vec![false; NODES];
    let mut conduit_cap = vec![f64::NAN; NODES];
    let mut conduit_paths: Vec<Vec<u32>> = Vec::new();
    let mut conduit_recv: Vec<(u32, u32)> = Vec::new();
    if p.coast.is_some() {
        let mut dw = vec![f64::INFINITY; NODES];
        let mut heap: BinaryHeap<FpKey> = BinaryHeap::new();
        for n in 0..NODES { if kind[n] == K_MARINE { dw[n] = 0.0; heap.push(FpKey(0.0, n as u32, n as u32)); } }
        let mut done = vec![false; NODES];
        while let Some(FpKey(c, i, _)) = heap.pop() {
            let i = i as usize;
            if done[i] || c > dw[i] { continue; }
            done[i] = true;
            let (gx, gy) = node_g(i);
            for (dx, dy) in D8 {
                let Some(j) = node_index(gx + dx, gy + dy) else { continue };
                if done[j] || kind[j] != K_LAND { continue; }
                let nc = c + ((dx * dx + dy * dy) as f64).sqrt();
                if nc < dw[j] { dw[j] = nc; heap.push(FpKey(nc, j as u32, i as u32)); }
            }
        }
        for t in 0..NODES {
            if kind[t] != K_T1A { continue; }
            let r = skel.receiver[t];
            let terminal = r == NONE || r == OFF_CHART || kind[r as usize] == K_LAND;
            if !terminal { continue; }
            let mut path = vec![t as u32];
            let mut cur = t;
            loop {
                let (gx, gy) = node_g(cur);
                let mut best: Option<(f64, usize)> = None;
                for (dx, dy) in D8 {
                    let Some(j) = node_index(gx + dx, gy + dy) else { continue };
                    let dj = if kind[j] == K_MARINE { 0.0 } else if kind[j] == K_LAND { dw[j] } else { continue };
                    if !dj.is_finite() || (cur != t && !(dj < dw[cur])) { continue; }
                    if best.map(|(bd, bj)| dj < bd || (dj == bd && j < bj)).unwrap_or(true) { best = Some((dj, j)); }
                }
                let Some((_, j)) = best else { break };
                path.push(j as u32);
                if kind[j] == K_MARINE || conduit[j] { break; }
                conduit[j] = true;
                conduit_cap[j] = level[t];
                conduit_recv.push((cur as u32, j as u32));
                cur = j;
            }
            if path.len() >= 2 && (kind[*path.last().unwrap() as usize] == K_MARINE || conduit[*path.last().unwrap() as usize]) {
                // The last hop (into water or an earlier conduit) is a link too.
                let k = path.len();
                conduit_recv.push((path[k - 2], path[k - 1]));
                conduit_paths.push(path);
            }
        }
    }
    // Lane T0-IMPRINT-PROD: the one confinement rule (params.confinement; absent = the cell rule): j may drain into i.
    let conf_cell = match p.confinement.as_ref().filter(|c| c.mode == crate::t1b_outlet::CONF_WARP) {
        Some(cp) => crate::t1b_outlet::warp_cells(cp, &skel.t0_cell, seed, recipe_digest)?,
        None => Vec::new(),
    };
    // Lane DRAIN-FREE conf=free: the node zones (land: own T0 cell in the core, NONE in the edge band; T1a: its outlet basin).
    let conf_cell = match p.confinement.as_ref().filter(|c| c.mode == crate::t1b_outlet::CONF_FREE) {
        Some(cp) => crate::t1b_outlet::free_zones(cp, t0, chart, &skel.t0_cell, &kind, &skel.outlet_basin)?,
        None => conf_cell,
    };
    let conf = crate::t1b_outlet::Confinement::new(p, t0, &skel.t0_cell, &conf_cell, &skel.receiver);
    let _ = &cell;
    let allowed = |i: usize, j: usize, kind: &[u8]| -> bool {
        if kind[i] == K_MARINE || kind[i] == K_LAKE || conduit[i] { return true; }
        conf.allows_nodes(t0, j, i, kind[i] == K_T1A)
    };
    let mut recv = vec![NONE; NODES];
    let mut in_tree = vec![false; NODES];
    let mut order: Vec<u32> = Vec::with_capacity(NODES);
    for n in 0..NODES { if kind[n] != K_LAND { in_tree[n] = true; order.push(n as u32); } }
    let nodes_land = (NODES - order.len()) as u32;
    // Conduits: receivers set, in order downstream first (each path appended reversed), routing seeds.
    let mut t1a_conduit_recv: Vec<(u32, u32)> = Vec::new();
    for path in &conduit_paths {
        for k in (1..path.len()).rev() {
            let (a, b) = (path[k - 1] as usize, path[k] as usize);
            if kind[a] == K_T1A { t1a_conduit_recv.push((a as u32, b as u32)); continue; }
            if in_tree[a] { continue; }
            recv[a] = b as u32; in_tree[a] = true; order.push(a as u32);
        }
    }
    let _ = &conduit_recv;
    let unit = |j: usize| (weight[j] >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0);
    let crosses_fn = |j: usize, i: usize, recv: &[u32]| -> bool {
        let (xj, yj) = node_g(j);
        let (xi, yi) = node_g(i);
        let (dx, dy) = (xi - xj, yi - yj);
        if dx == 0 || dy == 0 { return false; }
        let (Some(a), Some(b)) = (node_index(xj + dx, yj), node_index(xj, yj + dy)) else { return false };
        recv[a] == b as u32 || recv[b] == a as u32
    };
    let mut unreached_confined = 0u32;
    if p.network_mode == 0 {
        // Run r0: invasion percolation (Prim on node weights).
        let mut heap: BinaryHeap<std::cmp::Reverse<(u64, u32, u32)>> = BinaryHeap::new();
        let push = |i: usize, heap: &mut BinaryHeap<std::cmp::Reverse<(u64, u32, u32)>>, in_tree: &[bool], relaxed: bool| {
            let (gx, gy) = node_g(i);
            for (dx, dy) in D8 {
                let Some(j) = node_index(gx + dx, gy + dy) else { continue };
                if in_tree[j] || kind[j] != K_LAND { continue; }
                if !relaxed && !allowed(i, j, &kind) { continue; }
                heap.push(std::cmp::Reverse((weight[j], j as u32, i as u32)));
            }
        };
        for relaxed in [false, true] {
            for k in 0..order.len() { let i = order[k] as usize; push(i, &mut heap, &in_tree, relaxed); }
            while let Some(std::cmp::Reverse((_, j, i))) = heap.pop() {
                let (j, i) = (j as usize, i as usize);
                if in_tree[j] || crosses_fn(j, i, &recv) { continue; }
                recv[j] = i as u32; in_tree[j] = true; order.push(j as u32);
                push(j, &mut heap, &in_tree, relaxed);
            }
            if !relaxed { unreached_confined = (0..NODES).filter(|&n| !in_tree[n]).count() as u32; }
        }
    } else if p.network_mode == 2 {
        // Exploratory r2: random walk down the geodesic distance to the seeds.
        const K16: [(i64, i64); 16] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (2, 1), (1, 2), (-1, 2), (-2, 1), (-2, -1), (-1, -2), (1, -2), (2, -1)];
        let mut dg = vec![f64::INFINITY; NODES];
        let mut done = vec![false; NODES];
        let mut heap: BinaryHeap<FpKey> = BinaryHeap::new();
        for &n in &order { dg[n as usize] = 0.0; heap.push(FpKey(0.0, n, n)); }
        for relaxed in [false, true] {
            if relaxed { for n in 0..NODES { if done[n] { heap.push(FpKey(dg[n], n as u32, n as u32)); done[n] = false; } } }
            while let Some(FpKey(c, i, _)) = heap.pop() {
                let i = i as usize;
                if done[i] || c > dg[i] { continue; }
                done[i] = true;
                let (gx, gy) = node_g(i);
                for (dx, dy) in K16 {
                    let Some(j) = node_index(gx + dx, gy + dy) else { continue };
                    if done[j] || kind[j] != K_LAND { continue; }
                    if !relaxed && !allowed(i, j, &kind) { continue; }
                    let nc = c + ((dx * dx + dy * dy) as f64).sqrt();
                    if nc < dg[j] { dg[j] = nc; heap.push(FpKey(nc, j as u32, i as u32)); }
                }
            }
            if !relaxed { unreached_confined = (0..NODES).filter(|&n| kind[n] == K_LAND && !done[n]).count() as u32; }
        }
        let mut land: Vec<u32> = (0..NODES as u32).filter(|&n| kind[n as usize] == K_LAND && !conduit[n as usize] && dg[n as usize].is_finite()).collect();
        land.sort_by(|&a, &b| dg[a as usize].total_cmp(&dg[b as usize]).then(a.cmp(&b)));
        for &i in &land {
            let i = i as usize;
            let (gx, gy) = node_g(i);
            let mut cand: Vec<(usize, f64)> = Vec::new();
            let mut lowest: Option<(f64, usize)> = None;
            for (dx, dy) in D8 {
                let Some(j) = node_index(gx + dx, gy + dy) else { continue };
                if !(dg[j] < dg[i]) { continue; }
                if lowest.map(|(v, _)| dg[j] < v).unwrap_or(true) { lowest = Some((dg[j], j)); }
                if !allowed(j, i, &kind) || crosses_fn(i, j, &recv) { continue; }
                let step = ((dx * dx + dy * dy) as f64).sqrt();
                cand.push((j, (dg[i] - dg[j]) / step));
            }
            let pick = if cand.is_empty() { lowest.map(|(_, j)| j) } else {
                let tot: f64 = cand.iter().map(|c| c.1).sum();
                let u = unit(i) * tot;
                let mut acc = 0.0;
                let mut ch = cand[cand.len() - 1].0;
                for c in &cand { acc += c.1; if u < acc { ch = c.0; break; } }
                Some(ch)
            };
            if let Some(j) = pick { recv[i] = j as u32; in_tree[i] = true; order.push(i as u32); }
        }
    } else {
        // REVISION 1: first-passage percolation (Dijkstra on keyed node costs) replaces invasion percolation: IP paths are
        // tortuous, so neighbouring nodes had path lengths (and steady-state beds) kilometres apart (perched channels, 1.9 km
        // crest cuts in run s42_r0); near-geodesic FPP paths coalesce into trees without that disequilibrium.
        let mut dist = vec![f64::INFINITY; NODES];
        for &n in &order { dist[n as usize] = 0.0; }
        let mut heap: BinaryHeap<FpKey> = BinaryHeap::new();
        let push = |i: usize, heap: &mut BinaryHeap<FpKey>, in_tree: &[bool], dist: &mut [f64], relaxed: bool| {
            let (gx, gy) = node_g(i);
            for (dx, dy) in D8 {
                let Some(j) = node_index(gx + dx, gy + dy) else { continue };
                if in_tree[j] || kind[j] != K_LAND { continue; }
                if !relaxed && !allowed(i, j, &kind) { continue; }
                let step = if dx != 0 && dy != 0 { std::f64::consts::SQRT_2 } else { 1.0 };
                let c = dist[i] + step * (0.5 + unit(j));
                if c < dist[j] || relaxed { dist[j] = dist[j].min(c); }
                heap.push(FpKey(c, j as u32, i as u32));
            }
        };
        let _crosses_unused = |j: usize, i: usize, recv: &[u32]| -> bool {
            let (xj, yj) = node_g(j);
            let (xi, yi) = node_g(i);
            let (dx, dy) = (xi - xj, yi - yj);
            if dx == 0 || dy == 0 { return false; }
            let (Some(a), Some(b)) = (node_index(xj + dx, yj), node_index(xj, yj + dy)) else { return false };
            recv[a] == b as u32 || recv[b] == a as u32
        };
        for relaxed in [false, true] {
            for k in 0..order.len() { let i = order[k] as usize; push(i, &mut heap, &in_tree, &mut dist, relaxed); }
            while let Some(FpKey(c, j, i)) = heap.pop() {
                let (j, i) = (j as usize, i as usize);
                if in_tree[j] { continue; }
                if crosses_fn(j, i, &recv) { continue; }
                recv[j] = i as u32;
                in_tree[j] = true;
                dist[j] = c;
                order.push(j as u32);
                push(j, &mut heap, &in_tree, &mut dist, relaxed);
            }
            if !relaxed { unreached_confined = (0..NODES).filter(|&n| !in_tree[n]).count() as u32; }
        }

    }
    // Any node still out (isolated by refused crossings): drain to its lowest-index tree neighbour.
    for n in 0..NODES {
        if in_tree[n] { continue; }
        let (gx, gy) = node_g(n);
        if let Some(j) = D8.iter().filter_map(|(dx, dy)| node_index(gx + dx, gy + dy)).filter(|&j| in_tree[j]).min() { recv[n] = j as u32; in_tree[n] = true; order.push(n as u32); }
    }
    let mut area = vec![0f64; NODES];
    let a_node = (SPACING_MM as f64 / 1000.0).powi(2);
    // Lane CLIMATE-RUNOFF: a land node's own area is a_node x R(x) / r_ref (effective area; key on only).
    let clim = if p.climate_runoff.is_some() { Some(t0.climate.as_ref().ok_or("climate_runoff params set but no climate attached to T0")?) } else { None };
    for n in 0..NODES {
        if kind[n] == K_LAND {
            area[n] = a_node;
            if let Some(c) = clim {
                let (gx, gy) = node_g(n);
                let (la, lo) = chart.inverse_mm(gx * SPACING_MM, gy * SPACING_MM).map_err(|e| e.to_string())?;
                area[n] = a_node * c.weight_at_deg(t0, la * RAD2DEG, lo * RAD2DEG);
            }
        }
        if kind[n] == K_T1A { area[n] = skel.discharge_m3_s[n] as f64 * SECONDS_PER_YEAR / t0.runoff_m_per_yr; }
    }
    for &(t, c) in &t1a_conduit_recv { if kind[c as usize] == K_LAND { area[c as usize] += skel.discharge_m3_s[t as usize] as f64 * SECONDS_PER_YEAR / t0.runoff_m_per_yr; } }
    for k in (0..order.len()).rev() {
        let i = order[k] as usize;
        if kind[i] != K_LAND || recv[i] == NONE { continue; }
        let r = recv[i] as usize;
        if kind[r] == K_LAND { area[r] += area[i]; }
    }
    let mut pos = vec![[0.0; 2]; NODES];
    let jm = p.jitter_permille as f64 / 1000.0 * SPACING_MM as f64 / 1000.0;
    for n in 0..NODES {
        let (gx, gy) = node_g(n);
        let (x, y) = ((gx * SPACING_MM) as f64 / 1000.0, (gy * SPACING_MM) as f64 / 1000.0);
        pos[n] = if kind[n] == K_LAND && !conduit[n] { [x + (2.0 * jit[n].0 - 1.0) * jm, y + (2.0 * jit[n].1 - 1.0) * jm] } else { [x, y] };
    }
    let mut mid = vec![[0.0; 2]; NODES];
    let mut term = vec![NONE; NODES];
    for &n in &order {
        let n = n as usize;
        if kind[n] != K_LAND || recv[n] == NONE { continue; }
        let r = recv[n] as usize;
        let (a, b) = (pos[n], pos[r]);
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let len = dx.hypot(dy).max(1e-9);
        let off = (2.0 * jit[n].2 - 1.0) * p.meander_permille as f64 / 1000.0 * len;
        mid[n] = [0.5 * (a[0] + b[0]) - dy / len * off, 0.5 * (a[1] + b[1]) + dx / len * off];
        term[n] = if kind[r] == K_LAND { term[r] } else { r as u32 };
    }
    let t1a_recv = skel.receiver.clone();
    let t0_cell = skel.t0_cell.clone();
    let meander_u = jit.iter().map(|j| j.2).collect();
    let t1a_ch = (0..NODES).map(|n| kind[n] == K_T1A && skel.kind[n] == 1 && skel.channel[n] == 1).collect();
    if clipped_count > 0 { println!("lane TC clip: {clipped_count} T1a head nodes on the shelf band are construction land"); }
    // Lane DRAIN-FREE (flux mode only): the T1a seed discharge recomputed from this network's land catchments.
    let flux = p.confinement.as_ref().filter(|c| c.flux == 1).map(|c| crate::t1b_outlet::FluxBase {
        t1a_area0: (0..NODES).map(|n| if kind[n] == K_T1A { area[n] } else { 0.0 }).collect(), bound_mm: c.bound_mm,
        outlet_basin: skel.outlet_basin.clone(), t0_q_cell: (0..NODES).map(|n| if kind[n] == K_T1A { t0.discharge_m3_s[skel.t0_cell[n] as usize] } else { 0.0 }).collect(),
        r_ref_m_per_yr: t0.runoff_m_per_yr, own: Vec::new(), levels_t1a: c.levels == 1 });
    let mut net = Network { kind, recv, order, area_m2: area, pos, mid, term, level, unreached_confined, nodes_land, rockmass, conduit, conduit_cap, conduit_paths, t1a_recv, t0_cell, meander_u, t1a_ch, conf_cell, flux };
    if net.flux.is_some() { crate::t1b_outlet::flux_apply_built(&mut net); }
    Ok(net)
}

/// Lattice helpers for a node at chart (x, y) metres.
/// Jittered nodes on the chart edge are clamped to the chart for the lookup (their positions are not changed).
fn latlon_of(chart: &Chart, x: f64, y: f64) -> R<(f64, f64)> {
    let h = vz_tilekey::CHART_HALF_EXTENT_MM;
    let (xm, ym) = (((x * 1000.0).round() as i64).clamp(-h, h), ((y * 1000.0).round() as i64).clamp(-h, h));
    let (la, lo) = chart.inverse_mm(xm, ym).map_err(|e| e.to_string())?;
    Ok((la * RAD2DEG, lo * RAD2DEG))
}

/// E at the node positions (hat blend of the cell table).
pub fn e_at(t0: &T0, table: &BTreeMap<u32, f64>, la: f64, lo: f64) -> f64 {
    let (c, w) = hat4(t0, la, lo);
    let mut e = 0.0;
    for k in 0..4 { if w[k] > 0.0 { e += w[k] * table.get(&c[k]).copied().unwrap_or(0.0); } }
    e
}

/// Node beds and channel flags for an E table.
pub struct Beds { pub bed: Vec<f64>, pub chan: Vec<bool>,
    /// Lane CF (t1a_shelf_grade): effective T1a curve levels per node (None = the T1a trunk levels).
    pub t1a_level: Option<Vec<f64>>,
    /// Lane TD (capture routing): the rerouted network these beds belong to (None = the input network) and its stats.
    pub net: Option<Network>,
    pub capture: Option<crate::t1b_fix::CaptureStats>,
    /// Lane TA (LS1): level-solve statistics.
    pub ls: Option<crate::t1b_ls::LsStats>,
    /// Lane MEANDER-Q (cl_mq=1): rock erodibility multiplier per node for the valley-wall factor (set by the window binary).
    pub mq_rock: Option<Vec<f64>> }

impl Beds {
    /// The network these beds were computed on.
    pub fn net<'a>(&'a self, base: &'a Network) -> &'a Network { self.net.as_ref().unwrap_or(base) }
}

pub struct NodeCtx<'a> { pub latlon: Vec<(f64, f64)>, pub rock: RockPosts<'a>, pub bil: Vec<f64>,
    /// Lane MS-CEIL-2 (ceil with ceil_bed=1 only; empty otherwise): the ceiling Z_c at each node (CEIL_BED_RULE).
    pub ceil_z: Vec<f64> }

pub fn node_ctx<'a>(net: &Network, chart: &Chart, t0: &T0, extras: &T0Extras, p: &'a RiverNetParamsV1, strata: Option<&'a StrataCols>, anchor: (f64, f64)) -> R<NodeCtx<'a>> {
    let mut latlon = Vec::with_capacity(NODES);
    let mut pos4 = Vec::with_capacity(NODES);
    let mut bil = Vec::with_capacity(NODES);
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    // Lane TG G3: tau, the T1a shelf grade and final reach read the box-averaged bilinear.
    let sb = match crate::t1b_tg::tg(p).filter(|t| t.smooth_bil == 1 && p.coast.is_some()) { Some(t) => Some(crate::t1b_tg::chart_smooth(t0, chart, t.smooth_half_mm, t.smooth_grid_mm)?), None => None };
    let h = vz_tilekey::CHART_HALF_EXTENT_MM;
    for n in 0..NODES {
        let (la, lo) = latlon_of(chart, net.pos[n][0], net.pos[n][1])?;
        latlon.push((la, lo));
        pos4.push((la, lo, net.pos[n][0], net.pos[n][1]));
        bil.push(match &sb {
            Some(s) => s.at_mm(((net.pos[n][0] * 1000.0).round() as i64).clamp(-h, h) as f64, ((net.pos[n][1] * 1000.0).round() as i64).clamp(-h, h) as f64),
            None => if p.coast.is_some() { bilinear(t0, &elev, la, lo, false) } else { 0.0 },
        });
    }
    let rock = RockPosts::build(p, t0, extras, strata, anchor, &pos4)?;
    // Lane MS-CEIL-2 (ceil_bed=1): Z_c at the nodes, the lattice ceiling's law (ceil_surface) at the node position.
    let ceil_z: Vec<f64> = match p.ceil.as_ref().filter(|c| c.bed == 1) {
        Some(cp) => latlon.iter().map(|&(la, lo)| {
            let (hc, hw) = hat4(t0, la, lo);
            let mut r = 0.0f64;
            for k in 0..4 { if hw[k] > 0.0 { r += hw[k] * t0.subgrid_relief_m[hc[k] as usize] as f64; } }
            bilinear(t0, &elev, la, lo, false) + cp.c_h() * r
        }).collect(),
        None => Vec::new(),
    };
    Ok(NodeCtx { latlon, rock, bil, ceil_z })
}

pub fn node_beds(net: &Network, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, p: &RiverNetParamsV1) -> Beds {
    node_beds_warm(net, None, ctx, t0, laws, table, flood, p)
}

/// As `node_beds`; lane TD2 capture routing starts from `start` (a previous evaluation's network) when given.
#[allow(clippy::too_many_arguments)]
pub fn node_beds_warm(net: &Network, start: Option<&Network>, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, p: &RiverNetParamsV1) -> Beds {
    node_beds_warm_f(net, start, ctx, t0, laws, table, flood, p, false)
}

/// As `node_beds_warm`; lane TC C4 `frozen` (needs ls and a start network): LS1 does not run, the network is `start`.
#[allow(clippy::too_many_arguments)]
pub fn node_beds_warm_f(net: &Network, start: Option<&Network>, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, p: &RiverNetParamsV1, frozen: bool) -> Beds {
    // Lane TA LS1: the consistent level solve replaces TD2 capture routing when on.
    if let Some(ls) = p.ls.as_ref().filter(|l| l.solve == 1) {
        let (n2, st) = crate::t1b_ls::level_solve_f(net, start, ctx, t0, laws, table, p, ls, &|nn: &Network| node_beds_core(nn, ctx, t0, laws, table, flood, p).bed,
            &|nn: &Network| steep_labels(nn, ctx, t0, laws, table, flood, p), frozen && start.is_some());
        let mut b = node_beds_core(&n2, ctx, t0, laws, table, flood, p);
        b.net = Some(n2);
        b.ls = Some(st);
        return b;
    }
    // Lane TD2: capture routing (the network is rerouted for these tables), then the 0.9.2 bed rule on it.
    if let Some(fx) = p.fx.as_ref().filter(|f| f.capture == 1) {
        let (n2, st) = crate::t1b_fix::capture(net, start, ctx, t0, laws, table, p, fx, &|nn: &Network| node_beds_core(nn, ctx, t0, laws, table, flood, p).bed);
        let mut b = node_beds_core(&n2, ctx, t0, laws, table, flood, p);
        b.net = Some(n2);
        b.capture = Some(st);
        return b;
    }
    node_beds_core(net, ctx, t0, laws, table, flood, p)
}

/// The 0.9.2 bed rule on a given network.
pub fn node_beds_core(net: &Network, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, p: &RiverNetParamsV1) -> Beds {
    // Lane TB S2: the TD2b cap with the chain labels of a first pass on the same network.
    if crate::t1b_ls::steep_eik(p) {
        let first = node_beds_core_impl(net, ctx, t0, laws, table, flood, p, None);
        let lab = crate::t1b_ls::chain_labels(net, &first.chan);
        return node_beds_core_impl(net, ctx, t0, laws, table, flood, p, Some(&lab));
    }
    node_beds_core_impl(net, ctx, t0, laws, table, flood, p, None)
}

/// Lane TB S2: the chain labels node_beds_core uses for this network (None without the option).
pub fn steep_labels(net: &Network, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, p: &RiverNetParamsV1) -> Option<Vec<u32>> {
    if !crate::t1b_ls::steep_eik(p) { return None; }
    let first = node_beds_core_impl(net, ctx, t0, laws, table, flood, p, None);
    Some(crate::t1b_ls::chain_labels(net, &first.chan))
}

/// Lane TD2b slope cap at node n (with lane TB S2's eikonal inputs when `lab` is given): `s` is the stream-power slope
/// already capped at cap x Sc, `sc` the Sc at the receiver's bed `zr`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn steep_capped(s: f64, n: usize, zr: f64, len: f64, e: f64, sc: f64, net: &Network, ctx: &NodeCtx, laws: &ReliefLawsV1, p: &RiverNetParamsV1, fx: &crate::t1b_fix::FxParamsV1, lab: Option<&[u32]>) -> f64 {
    let g = fx.steep_gamma_permille as f64 / 1000.0;
    match lab {
        None => s.min(g * crate::t1b_fix::side_gradient(e, sc, net.rockmass[n], crate::t1b_fix::link_hill_length(len), laws, p)),
        Some(lab) => {
            let rm = if lab[n] != NONE { net.rockmass[lab[n] as usize] } else { net.rockmass[n] };
            let s1 = s.min(g * crate::t1b_fix::side_gradient(e, sc, rm, crate::t1b_fix::link_hill_length(len), laws, p));
            let sc_top = ctx.rock.at(n, zr + s1 * len).0;
            s1.min(g * crate::t1b_fix::side_gradient(e, sc.min(sc_top), rm, crate::t1b_fix::link_hill_length(len), laws, p))
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn node_beds_core_impl(net: &Network, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, p: &RiverNetParamsV1, lab: Option<&[u32]>) -> Beds {
    let theta = laws.stream_power_m / laws.stream_power_n;
    let cap = p.bed_slope_cap_permille_of_sc as f64 / 1000.0;
    let mut bed = net.level.clone();
    let mut chan = vec![false; NODES];
    // Lane CO: construction water stands at the base level sea - F (hat blend of the per-cell flooding offsets).
    if p.coast.is_some() {
        for n in 0..NODES { if net.kind[n] == K_MARINE { let (la, lo) = ctx.latlon[n]; bed[n] = t0.sea_level_m - e_at(t0, flood, la, lo); } }
    }
    let tau = |n: usize| p.coast.as_ref().map(|c| c.tau(ctx.bil[n])).unwrap_or(1.0);
    // The bed slope of node n over a link of length len from a receiver bed zr (the 0.9.2 / 0.9.3 bed law).
    let slope = |n: usize, zr: f64, area: f64, len: f64| -> f64 {
        let (sc, km, _) = ctx.rock.at(n, zr);
        let (la, lo) = ctx.latlon[n];
        let e = e_at(t0, table, la, lo) * tau(n);
        let ks = e / (laws.stream_power_k_per_my * km);
        let mut s = (ks * area.powf(-theta)).min(cap * sc);
        // Lane TD2b (REVISION 1): the channel is never steeper than gamma x its valley sides.
        if let Some(fx) = p.fx.as_ref().filter(|f| f.steep_cap == 1) { s = steep_capped(s, n, zr, len, e, sc, net, ctx, laws, p, fx, lab); }
        s
    };
    // Lane MS-CEIL-2 ceiling_bed_grade_v1 (ceil_bed=1; CEIL_BED_RULE): cap(n) = Z_c(n) - S_h(n) L_n, the ceiling minus the
    // valley-side relief of n's link (Sc of the rock at the ceiling height), over land (not conduit) and T1a nodes; the bed law
    // takes max(bed_r, min(law bed, cap(n))); a post-pass (below) lowers the reaches a ceiling falling upstream left above it.
    let cgraded = |n: usize| (net.kind[n] == K_LAND && !net.conduit[n]) || net.kind[n] == K_T1A;
    let crecv = |n: usize| -> u32 { if net.kind[n] == K_T1A { net.t1a_recv[n] } else { net.recv[n] } };
    let cdown = |n: usize| -> Option<usize> { let r = crecv(n); if r == NONE || r == OFF_CHART || !cgraded(r as usize) { None } else { Some(r as usize) } };
    let (capl, clen, ctopo): (Vec<f64>, Vec<f64>, Vec<usize>) = if ctx.ceil_z.is_empty() { (Vec::new(), Vec::new(), Vec::new()) } else {
        let mut c = vec![f64::INFINITY; NODES];
        let mut clen = vec![0f64; NODES];
        let mut indeg = vec![0u32; NODES];
        for n in 0..NODES {
            if !cgraded(n) { continue; }
            let r = crecv(n);
            let len = if r == NONE || r == OFF_CHART { SPACING_MM as f64 / 1000.0 } else if net.kind[n] == K_T1A {
                (net.pos[n][0] - net.pos[r as usize][0]).hypot(net.pos[n][1] - net.pos[r as usize][1])
            } else {
                let (a, m, b) = (net.pos[n], net.mid[n], net.pos[r as usize]);
                (a[0] - m[0]).hypot(a[1] - m[1]) + (m[0] - b[0]).hypot(m[1] - b[1])
            };
            let sc = ctx.rock.at(n, ctx.ceil_z[n]).0;
            let (la, lo) = ctx.latlon[n];
            let e = e_at(t0, table, la, lo) * tau(n);
            let l = crate::t1b_fix::link_hill_length(len);
            // never below present sea level: Z_c is a subaerial surface; near marine T0 cells the bilinear is pulled far below it
            c[n] = (ctx.ceil_z[n] - crate::t1b_fix::side_gradient(e, sc, net.rockmass[n], l, laws, p) * l).max(t0.sea_level_m);
            clen[n] = len;
            if let Some(d) = cdown(n) { indeg[d] += 1; }
        }
        // upstream-first order of the graded nodes (Kahn on the receiver links)
        let mut topo = Vec::new();
        let mut q: std::collections::VecDeque<usize> = (0..NODES).filter(|&n| cgraded(n) && indeg[n] == 0).collect();
        while let Some(n) = q.pop_front() {
            topo.push(n);
            if let Some(d) = cdown(n) { indeg[d] -= 1; if indeg[d] == 0 { q.push_back(d); } }
        }
        (c, clen, topo)
    };
    let ceil_cap = |n: usize| -> Option<f64> { if capl.is_empty() || !capl[n].is_finite() { None } else { Some(capl[n]) } };
    let t1a_mode = crate::t1b_ls::t1a_mode(p);
    // Lane CF C6b: T1a levels seaward of the T0 zero contour graded to the lowstand (sea - F min(1, -bil / shelf)).
    let mut t1a_level = None;
    let cf_grade = p.cf.as_ref().map(|c| c.t1a_shelf_grade == 1).unwrap_or(false) && p.coast.is_some();
    let sea = t0.sea_level_m;
    let shelf = p.coast.as_ref().map(|cp| cp.shelf_m().max(1.0)).unwrap_or(1.0);
    // R1d: only nodes whose whole downstream T1a chain lies at bil <= sea (the last reach to the mouth).
    let final_reach = |n: usize| -> bool {
        let mut c = n;
        for _ in 0..NODES {
            let r = net.t1a_recv[c];
            if r == NONE || r == OFF_CHART { return true; }
            let r = r as usize;
            if net.kind[r] != K_T1A { return true; }
            if ctx.bil[r] > sea { return false; }
            c = r;
        }
        true
    };
    let shelf_cap = |n: usize| -> Option<f64> {
        if cf_grade && net.kind[n] == K_T1A && ctx.bil[n] <= sea && final_reach(n) {
            let (la, lo) = ctx.latlon[n];
            Some(sea - e_at(t0, flood, la, lo) * ((sea - ctx.bil[n]) / shelf).min(1.0))
        } else { None }
    };
    if t1a_mode == 0 {
        if cf_grade {
            for n in 0..NODES { if let Some(c) = shelf_cap(n) { if c < bed[n] { bed[n] = c; } } }
            t1a_level = Some((0..NODES).map(|n| if net.kind[n] == K_T1A || net.kind[n] == K_MARINE { bed[n] } else { f64::INFINITY }).collect::<Vec<f64>>());
        }
    } else {
        // Lane TB S1 (R2 EXPLORATORY): conduits first (bed law from their receivers, no terminal cap), then the T1a levels
        // by the bed law downstream first along the T1a receivers.
        // Lane DRAIN-FIX conf_levels=1 (t1b_outlet::LEVELS_RULE): the T1a discharge areas, not the recomputed flux areas.
        let lva = crate::t1b_outlet::level_areas(net);
        let lv = |n: usize| lva.as_ref().map(|v| v[n]).unwrap_or(net.area_m2[n]);
        for &n in &net.order {
            let n = n as usize;
            if net.kind[n] != K_LAND || !net.conduit[n] || net.recv[n] == NONE { continue; }
            let r = net.recv[n] as usize;
            let zr = bed[r];
            let (a, m, b) = (net.pos[n], net.mid[n], net.pos[r]);
            let len = (a[0] - m[0]).hypot(a[1] - m[1]) + (m[0] - b[0]).hypot(m[1] - b[1]);
            bed[n] = zr.max(zr + slope(n, zr, lv(n), len) * len);
            chan[n] = true;
        }
        let mut conduit_out = vec![NONE; NODES];
        for path in &net.conduit_paths { if path.len() >= 2 { conduit_out[path[0] as usize] = path[1]; } }
        let outlet = |t: usize| -> Option<usize> {
            let r = net.t1a_recv[t];
            if r != NONE && r != OFF_CHART && net.kind[r as usize] == K_T1A { return Some(r as usize); }
            if conduit_out[t] != NONE { return Some(conduit_out[t] as usize); }
            if r != NONE && r != OFF_CHART && (net.kind[r as usize] == K_MARINE || net.kind[r as usize] == K_LAKE) { return Some(r as usize); }
            None
        };
        let mut done = vec![false; NODES];
        for n in 0..NODES { if net.kind[n] != K_T1A { done[n] = true; } }
        for t0n in 0..NODES {
            if done[t0n] { continue; }
            // Walk down to the first settled outlet, then settle the stack upward.
            let mut stack = vec![t0n];
            loop {
                let t = *stack.last().unwrap();
                match outlet(t) { Some(o) if !done[o] && !stack.contains(&o) => stack.push(o), _ => break }
            }
            while let Some(t) = stack.pop() {
                let seed = t1a_mode == 2 && shelf_cap(t).is_some();
                let mut z = match (seed, outlet(t)) {
                    (false, Some(o)) if done[o] && bed[o].is_finite() => {
                        let len = (net.pos[t][0] - net.pos[o][0]).hypot(net.pos[t][1] - net.pos[o][1]);
                        let zl = bed[o] + slope(t, bed[o], lv(t), len) * len;
                        // Lane MS-CEIL-2 (ceil_bed=1): never above the ceiling minus the valley-side relief, never below the outlet.
                        match ceil_cap(t) { Some(c) => bed[o].max(zl.min(c)), None => zl }
                    }
                    _ => net.level[t],
                };
                if let Some(c) = shelf_cap(t) { z = z.min(c); }
                bed[t] = z;
                done[t] = true;
            }
        }
        t1a_level = Some((0..NODES).map(|n| if net.kind[n] == K_T1A || net.kind[n] == K_MARINE || net.kind[n] == K_LAKE || net.conduit[n] { bed[n] } else { f64::INFINITY }).collect::<Vec<f64>>());
    }
    // Lane CF C6b: a conduit never rises above its (graded) terminal T1a node.
    let ccap: Option<Vec<f64>> = if t1a_mode == 0 && t1a_level.is_some() {
        let mut c = net.conduit_cap.clone();
        for path in &net.conduit_paths { let t = path[0] as usize; for &j in &path[1..] { let j = j as usize; if net.conduit[j] { c[j] = c[j].min(bed[t]); } } }
        Some(c)
    } else { None };
    // Lane TC C1 floor 2: T1a floor cap per land node (from the T1a beds above).
    let fcap = crate::t1b_ls::t1a_floor_cap(net, &bed, p);
    for &n in &net.order {
        let n = n as usize;
        if net.kind[n] != K_LAND || net.recv[n] == NONE { continue; }
        if t1a_mode != 0 && net.conduit[n] { continue; }
        let r = net.recv[n] as usize;
        let zr = bed[r];
        let (a, m, b) = (net.pos[n], net.mid[n], net.pos[r]);
        let len = (a[0] - m[0]).hypot(a[1] - m[1]) + (m[0] - b[0]).hypot(m[1] - b[1]);
        let s = slope(n, zr, net.area_m2[n], len);
        let km = ctx.rock.at(n, zr).1;
        bed[n] = zr + s * len;
        if let Some(c) = &fcap { if c[n].is_finite() && !net.conduit[n] { bed[n] = zr.max(bed[n].min(c[n])); } }
        // Lane MS-CEIL-2 (ceil_bed=1): never above the ceiling minus the valley-side relief, never below the receiver.
        if !net.conduit[n] { if let Some(c) = ceil_cap(n) { bed[n] = zr.max(bed[n].min(c)); } }
        chan[n] = net.area_m2[n] >= p.channel_head_area_m2 as f64 / km;
        if net.conduit[n] { let cap = match &ccap { Some(c) => c[n], None => net.conduit_cap[n] }; bed[n] = zr.max(bed[n].min(cap)); chan[n] = true; }
    }
    // Lane MS-CEIL-2 (ceil_bed=1): a river never rises. A graded node the bed law left above its cap (held at its receiver
    // where the ceiling falls upstream) is lowered to the cap, and the reach below it falls at the floor slope until it meets
    // its own graded profile (never below the fixed water / conduit level it ends at); tributaries keep their beds.
    if !capl.is_empty() {
        let fsl = p.floor_slope_ppm as f64 / 1e6;
        // only land network nodes are lowered here: a T1a trunk (the record's trunk, graded by the bed law with the local cap
        // above) is a base like water, so the lowering never runs down a low-gradient trunk to the lowstand
        let land = |n: usize| net.kind[n] == K_LAND && !net.conduit[n];
        let ldown = |n: usize| -> Option<usize> { cdown(n).filter(|&d| land(d)) };
        let mut base = vec![f64::NEG_INFINITY; NODES];
        for &n in ctopo.iter().rev() {
            if !land(n) { continue; }
            base[n] = match ldown(n) { Some(d) => base[d], None => { let r = crecv(n); if r == NONE || r == OFF_CHART { f64::NEG_INFINITY } else { bed[r as usize] } } };
        }
        let mut low = vec![false; NODES];
        for &n in &ctopo {
            if !land(n) { continue; }
            if bed[n] > capl[n] + 1e-9 { bed[n] = capl[n].max(base[n]); low[n] = true; }
            if !low[n] { continue; }
            if let Some(d) = ldown(n) {
                let v = (bed[n] - fsl * clen[n]).max(base[d]);
                if v < bed[d] { bed[d] = v; low[d] = true; }
            }
        }
    }
    for k in (0..net.order.len()).rev() {
        let n = net.order[k] as usize;
        if net.kind[n] != K_LAND || !chan[n] || net.recv[n] == NONE { continue; }
        let r = net.recv[n] as usize;
        if net.kind[r] == K_LAND { chan[r] = true; }
    }
    Beds { bed, chan, t1a_level, net: None, capture: None, ls: None, mq_rock: None }
}

/// A channel curve in chart metres: points (x, y, bed, area).
pub struct Curve { pub pts: Vec<[f64; 4]>, pub code: u8, pub link: u32, pub label: u32,
    /// Lane TC C5 (ls_curve, fix B1): the T1a link as published without the meander (Catmull-Rom through the nodes), used for
    /// volcano siting so the drawn curve does not move vents; None = `pts`.
    pub site: Option<Vec<[f64; 4]>>,
    /// Lane CHANNEL-CENTRELINE (cl_mn=1): the unmeandered T1a-curve link (floor extent) and the link's largest belt offset (m).
    pub carrier: Option<(Vec<[f64; 4]>, f64)>,
    /// Lane MEANDER-Q (cl_mq=1): the chart's swept-envelope floodplain, shared by every migrated T1a link.
    pub mq: Option<std::sync::Arc<crate::t1b_mq::MqEnv>> }

fn chaikin4(p: &[[f64; 4]]) -> Vec<[f64; 4]> {
    if p.len() < 3 { return p.to_vec(); }
    let mut q = Vec::with_capacity(2 * p.len());
    q.push(p[0]);
    let m = p.len() - 1;
    for k in 0..m {
        let (a, b) = (p[k], p[k + 1]);
        let mix = |u: f64| [a[0] + (b[0] - a[0]) * u, a[1] + (b[1] - a[1]) * u, a[2] + (b[2] - a[2]) * u, a[3] + (b[3] - a[3]) * u];
        if k > 0 { q.push(mix(0.25)); }
        if k + 1 < m { q.push(mix(0.75)); }
    }
    q.push(p[m]);
    q
}

/// All channel curves of the chart, in a fixed order: T1a links (node order), then network chains (start node order).
pub fn curves(net: &Network, beds: &Beds, skel: &SkeletonView, t0: &T0, p: &RiverNetParamsV1) -> Vec<Curve> {
    curves_with_stats(net, beds, skel, t0, p, None).0
}

/// `curves` plus the lane CHANNEL-CENTRELINE statistics (None without params.cl); `confine`: the meander belt's
/// largest floor half-width at a point from other constructions (volcano edifices), None = unconfined.
pub fn curves_with_stats(net: &Network, beds: &Beds, skel: &SkeletonView, t0: &T0, p: &RiverNetParamsV1, confine: Option<&dyn Fn([f64; 2]) -> f64>) -> (Vec<Curve>, Option<crate::t1b_cl::ClStats>) {
    let net = beds.net(net);
    let mut out = Vec::new();
    // Lane DRAIN-FREE flux mode: a T1a seed's area is the network's (the T1a discharge recomputed from the actual T1b
    // catchments, t1b_outlet::flux_seed_areas); drowned axis nodes (not seeds) keep the T1a discharge.
    let flux = net.flux.is_some();
    let t1a_area = |n: usize| if flux && net.kind[n] == K_T1A { net.area_m2[n] } else { skel.discharge_m3_s[n] as f64 * SECONDS_PER_YEAR / t0.runoff_m_per_yr };
    let lattice = |n: usize| { let (gx, gy) = node_g(n); [(gx * SPACING_MM) as f64 / 1000.0, (gy * SPACING_MM) as f64 / 1000.0] };
    // T1a valley axes: channel==1 links are T1a channels (code 1, F1); sub-threshold axes (soft axes, head extensions) are
    // network channels (code 5) on the same T1a planform and levels.
    // Lane TC C2: clipped head nodes (construction land in the network) are not T1a.
    let clip = crate::t1b_ls::tc_clip(p);
    let is_ch = |n: usize| skel.kind[n] == 1 && !(clip && net.kind[n] == K_LAND);
    let ok_r = |r: u32| r != NONE && r != OFF_CHART;
    // Lane TC C5: the keyed meander midpoint of the T1a link a -> receiver (no offset on construction water).
    let tcurve = crate::t1b_ls::tc_curve(p);
    let mid_of = |a: usize| -> Option<[f64; 2]> {
        if !is_ch(a) || !ok_r(skel.receiver[a]) { return None; }
        let (pa, pr) = (lattice(a), lattice(skel.receiver[a] as usize));
        let (dx, dy) = (pr[0] - pa[0], pr[1] - pa[1]);
        let len = dx.hypot(dy).max(1e-9);
        let off = if net.kind[a] == K_T1A { (2.0 * net.meander_u[a] - 1.0) * p.meander_permille as f64 / 1000.0 * len } else { 0.0 };
        Some([0.5 * (pa[0] + pr[0]) - dy / len * off, 0.5 * (pa[1] + pr[1]) + dx / len * off])
    };
    // Lane CHANNEL-CENTRELINE (t1b_cl.rs): T1a links on the T1a planform curve (and the optional declared meander belt).
    let clp = crate::t1b_cl::cl(p);
    let mut cls = crate::t1b_cl::ClStats::default();
    let mut tabs: BTreeMap<u32, crate::t1b_cl::SegTab> = BTreeMap::new();
    let mut node_k: BTreeMap<usize, (u32, usize)> = BTreeMap::new();
    let mut jpos: BTreeMap<usize, [f64; 2]> = BTreeMap::new();
    // Lane MEANDER-Q (cl_mq=1, t1b_mq.rs): migrated channel per T1a segment and the chart's swept-envelope floodplain.
    let mut mq_segs: BTreeMap<u32, crate::t1b_mq::SegOut> = BTreeMap::new();
    let mut mq_env: Option<std::sync::Arc<crate::t1b_mq::MqEnv>> = None;
    if let Some(cp) = clp {
        let nseg = skel.t1a_segments.len() as u32;
        for n in 0..NODES {
            if !(skel.kind[n] == 1 || skel.kind[n] == 3) { continue; }
            let sg = skel.segment[n];
            if sg >= nseg { continue; }
            let tab = tabs.entry(sg).or_insert_with(|| crate::t1b_cl::SegTab::new(&skel.t1a_segments[sg as usize]));
            node_k.insert(n, (sg, tab.project(lattice(n))));
        }
        if cp.mq_on() {
            let unconf = |_: [f64; 2]| f64::INFINITY;
            let (segs, env, mst) = crate::t1b_mq::run_chart(&tabs, &node_k, skel, beds, t0, p, confine.unwrap_or(&unconf)).unwrap_or_else(|e| panic!("lane MEANDER-Q: {e}"));
            mq_segs = segs; mq_env = Some(std::sync::Arc::new(env)); cls.mq = Some(mst);
        } else if cp.sine_on() {
            let sgt = crate::t1b_cl::SgTables::new();
            let vfd = crate::t1b_vf::VfParamsV1::off();
            let vfp = p.vf.as_ref().unwrap_or(&vfd);
            let cw = p.floor_width_coeff_ppm as f64 / 1e6;
            let hmax = p.floor_width_max_mm as f64 / 2000.0;
            let hyd = |a: f64| -> (f64, f64) { (vfp.hydraulic(a).0, (0.5 * cw * a.max(0.0).sqrt()).min(hmax)) };
            let mut per: BTreeMap<u32, Vec<crate::t1b_cl::SegNode>> = BTreeMap::new();
            for (&n, &(sg, k)) in &node_k { per.entry(sg).or_default().push(crate::t1b_cl::SegNode { k, area_m2: t1a_area(n), u: net.meander_u[n] }); }
            for (sg, mut v) in per {
                v.sort_by(|a, b| a.k.cmp(&b.k).then(a.area_m2.total_cmp(&b.area_m2)));
                let tab = tabs.get_mut(&sg).unwrap();
                let unconf = |_: [f64; 2]| f64::INFINITY;
                let m = crate::t1b_cl::meander(tab, &v, cp, &sgt, &hyd, confine.unwrap_or(&unconf), sg);
                cls.meander.push(m);
            }
        }
        let (mut sh_sum, mut sh_n) = (0.0f64, 0u32);
        for (&n, &(sg, k)) in &node_k {
            let q = match mq_segs.get(&sg) { Some(o) => o.point_at(k as f64), None => tabs[&sg].pt(k) };
            let l = lattice(n);
            let d = (q[0] - l[0]).hypot(q[1] - l[1]);
            cls.node_shift_max_m = cls.node_shift_max_m.max(d);
            sh_sum += d; sh_n += 1;
            jpos.insert(n, q);
        }
        cls.node_shift_mean_m = if sh_n > 0 { sh_sum / sh_n as f64 } else { 0.0 };
    }
    let jp = |n: usize| -> [f64; 2] { jpos.get(&n).copied().unwrap_or_else(|| lattice(n)) };
    for n in 0..NODES {
        if !is_ch(n) || !ok_r(skel.receiver[n]) { continue; }
        let r = skel.receiver[n] as usize;
        let (pn, pr) = (lattice(n), lattice(r));
        let (gx, gy) = node_g(n);
        let mut prev: Option<(f32, usize)> = None;
        for (dx, dy) in D8 {
            if let Some(d) = node_index(gx + dx, gy + dy) {
                if is_ch(d) && skel.receiver[d] == n as u32 {
                    let q = skel.discharge_m3_s[d];
                    if prev.map(|(pq, _)| q > pq).unwrap_or(true) { prev = Some((q, d)); }
                }
            }
        }
        let pp = prev.map(|(_, d)| lattice(d)).unwrap_or([2.0 * pn[0] - pr[0], 2.0 * pn[1] - pr[1]]);
        let pq = if ok_r(skel.receiver[r]) { lattice(skel.receiver[r] as usize) } else { [2.0 * pr[0] - pn[0], 2.0 * pr[1] - pn[1]] };
        let (mut la, mut lb) = (skel.trunk_level_m[n] as f64, skel.trunk_level_m[r] as f64);
        if let Some(tl) = &beds.t1a_level {
            // Lane TB S1: the solved levels replace the fixed ones (an outlet without a solved bed keeps the link flat).
            if crate::t1b_ls::t1a_mode(p) != 0 { la = tl[n]; lb = if tl[r].is_finite() { tl[r].min(la) } else { la }; }
            else { la = la.min(tl[n]); lb = lb.min(tl[r]).min(la); }
        }
        let a = t1a_area(n);
        let mut curved: Option<Vec<[f64; 4]>> = None;
        let mut carrier: Option<(Vec<[f64; 4]>, f64)> = None;
        if clp.is_some() {
            cls.t1a_links += 1;
            // (x, y) of the channel and of the unmeandered curve, then beds linear in channel path length.
            let mut xy: Vec<([f64; 2], [f64; 2])> = Vec::new();
            let mut belt = 0.0f64;
            match node_k.get(&n) {
                Some(&(sg, kn)) if mq_segs.contains_key(&sg) => {
                    // Lane MEANDER-Q: the link follows the migrated channel between its nodes' carrier labels.
                    let o = &mq_segs[&sg];
                    let same = node_k.get(&r).map(|&(sr, _)| sr == sg).unwrap_or(false);
                    if same {
                        let kr = node_k[&r].1;
                        if kr >= kn { cls.links_on_curve += 1; for q in o.between(kn as f64, kr as f64) { xy.push((q, q)); } }
                        else { cls.nonmonotone_links_straight += 1; let (a, b) = (o.point_at(kn as f64), o.point_at(kr as f64)); xy.push((a, a)); xy.push((b, b)); }
                    } else {
                        cls.junction_links += 1;
                        let jr = jp(r);
                        let je = o.nearest_from(jr, kn as f64);
                        for q in o.to_index(kn as f64, je) { xy.push((q, q)); }
                        let last = xy.last().unwrap().0;
                        let gap = (last[0] - jr[0]).hypot(last[1] - jr[1]);
                        cls.connector_max_m = cls.connector_max_m.max(gap);
                        if gap > 1e-6 { xy.push((jr, jr)); }
                    }
                }
                Some(&(sg, kn)) => {
                    let tab = &tabs[&sg];
                    let same = node_k.get(&r).map(|&(sr, _)| sr == sg).unwrap_or(false);
                    if same {
                        let kr = node_k[&r].1;
                        if kr >= kn {
                            cls.links_on_curve += 1;
                            for k in kn..=kr { xy.push((tab.pt(k), tab.x[k])); if let Some(y) = tab.y.get(k) { belt = belt.max(y.abs()); } }
                        } else {
                            cls.nonmonotone_links_straight += 1;
                            xy.push((tab.pt(kn), tab.x[kn]));
                            xy.push((tab.pt(kr), tab.x[kr]));
                        }
                    } else {
                        cls.junction_links += 1;
                        let jr = jp(r);
                        let ke = tab.project_from(jr, kn);
                        for k in kn..=ke { xy.push((tab.pt(k), tab.x[k])); if let Some(y) = tab.y.get(k) { belt = belt.max(y.abs()); } }
                        let last = xy.last().unwrap().0;
                        let gap = (last[0] - jr[0]).hypot(last[1] - jr[1]);
                        cls.connector_max_m = cls.connector_max_m.max(gap);
                        if gap > 1e-6 { let jc = node_k.get(&r).map(|&(sr, kr)| tabs[&sr].x[kr]).unwrap_or(jr); xy.push((jr, jc)); }
                    }
                }
                None => {
                    cls.links_without_segment += 1;
                    xy.push((jp(n), jp(n)));
                    xy.push((jp(r), jp(r)));
                }
            }
            if xy.len() < 2 { let q = xy[0]; xy.push(q); }
            let mut cum = vec![0.0f64; xy.len()];
            for k in 1..xy.len() { cum[k] = cum[k - 1] + (xy[k].0[0] - xy[k - 1].0[0]).hypot(xy[k].0[1] - xy[k - 1].0[1]); }
            let tot = cum[xy.len() - 1].max(1e-9);
            let bed = |k: usize| la + (lb - la) * cum[k] / tot;
            let pts: Vec<[f64; 4]> = (0..xy.len()).map(|k| [xy[k].0[0], xy[k].0[1], bed(k), a]).collect();
            if clp.map(|c| c.sine_on()).unwrap_or(false) {
                carrier = Some(((0..xy.len()).map(|k| [xy[k].1[0], xy[k].1[1], bed(k), a]).collect(), belt));
            }
            curved = Some(pts);
        } else if tcurve {
            // Two Catmull-Rom spans: n -> m_n -> r, outer controls the neighbouring links' midpoints (else reflected).
            let mn = mid_of(n).unwrap();
            let q0 = prev.and_then(|(_, d)| mid_of(d)).unwrap_or([2.0 * pn[0] - mn[0], 2.0 * pn[1] - mn[1]]);
            let q4 = mid_of(r).unwrap_or([2.0 * pr[0] - mn[0], 2.0 * pr[1] - mn[1]]);
            let mut pts = Vec::with_capacity(33);
            for (h, c) in [[q0, pn, mn, pr], [pn, mn, pr, q4]].iter().enumerate() {
                for s in 0..=16 {
                    if h == 1 && s == 0 { continue; }
                    let t = s as f64 / 16.0;
                    let (t2, t3) = (t * t, t * t * t);
                    let w = [0.5 * (-t3 + 2.0 * t2 - t), 0.5 * (3.0 * t3 - 5.0 * t2 + 2.0), 0.5 * (-3.0 * t3 + 4.0 * t2 + t), 0.5 * (t3 - t2)];
                    let x = w[0] * c[0][0] + w[1] * c[1][0] + w[2] * c[2][0] + w[3] * c[3][0];
                    let y = w[0] * c[0][1] + w[1] * c[1][1] + w[2] * c[2][1] + w[3] * c[3][1];
                    let u = (h as f64 * 16.0 + s as f64) / 32.0;
                    pts.push([x, y, la + (lb - la) * u, a]);
                }
            }
            curved = Some(pts);
        }
        let mut pts = Vec::with_capacity(17);
        for s in 0..=16 {
            let t = s as f64 / 16.0;
            let (t2, t3) = (t * t, t * t * t);
            let w = [0.5 * (-t3 + 2.0 * t2 - t), 0.5 * (3.0 * t3 - 5.0 * t2 + 2.0), 0.5 * (-3.0 * t3 + 4.0 * t2 + t), 0.5 * (t3 - t2)];
            let x = w[0] * pp[0] + w[1] * pn[0] + w[2] * pr[0] + w[3] * pq[0];
            let y = w[0] * pp[1] + w[1] * pn[1] + w[2] * pr[1] + w[3] * pq[1];
            pts.push([x, y, la + (lb - la) * t, a]);
        }
        if let Some(cp) = curved {
            let mq = node_k.get(&n).and_then(|&(sg, _)| mq_segs.get(&sg)).and(mq_env.clone());
            out.push(Curve { pts: cp, code: if skel.channel[n] == 1 { CHANNEL } else { GUIDE }, link: n as u32, label: n as u32, site: Some(pts), carrier, mq });
            continue;
        }
        out.push(Curve { pts, code: if skel.channel[n] == 1 { CHANNEL } else { GUIDE }, link: n as u32, label: n as u32, site: None, carrier: None, mq: None });
    }
    // Network chains.
    let mut donors = vec![0u32; NODES];
    for n in 0..NODES { if net.kind[n] == K_LAND && beds.chan[n] && net.recv[n] != NONE { donors[net.recv[n] as usize] += 1; } }
    for s in 0..NODES {
        if net.kind[s] != K_LAND || !beds.chan[s] || net.recv[s] == NONE || donors[s] == 1 || net.conduit[s] { continue; }
        let mut pts = vec![[net.pos[s][0], net.pos[s][1], beds.bed[s], net.area_m2[s]]];
        let mut cur = s;
        loop {
            let r = net.recv[cur] as usize;
            let m = net.mid[cur];
            pts.push([m[0], m[1], 0.5 * (beds.bed[cur] + beds.bed[r]), net.area_m2[cur]]);
            let ar = if net.kind[r] == K_LAND { net.area_m2[r] } else { net.area_m2[cur] };
            // Lane CHANNEL-CENTRELINE: a chain ending on a T1a node ends on the node's point of the T1a curve.
            let pr = if clp.is_some() && net.kind[r] == K_T1A { jp(r) } else { net.pos[r] };
            pts.push([pr[0], pr[1], beds.bed[r], ar]);
            if net.kind[r] != K_LAND || donors[r] != 1 || net.recv[r] == NONE { break; }
            cur = r;
        }
        for _ in 0..p.chaikin_passes { pts = chaikin4(&pts); }
        out.push(Curve { pts, code: GUIDE, link: net.term[s], label: s as u32, site: None, carrier: None, mq: None });
    }
    // Lane CO conduits: terminal T1a node, then down the conduit chain (following receivers) to the water.
    for path in &net.conduit_paths {
        let t = path[0] as usize;
        let pt0 = if clp.is_some() { jp(t) } else { net.pos[t] };
        let mut pts = vec![[pt0[0], pt0[1], beds.bed[t], t1a_area(t)]];
        let mut cur = path[1] as usize;
        let mut guard = 0;
        loop {
            pts.push([net.pos[cur][0], net.pos[cur][1], beds.bed[cur], net.area_m2[cur].max(t1a_area(t))]);
            if net.kind[cur] != K_LAND || net.recv[cur] == NONE || guard > NODES { break; }
            let r = net.recv[cur] as usize;
            let m = net.mid[cur];
            pts.push([m[0], m[1], 0.5 * (beds.bed[cur] + beds.bed[r]), net.area_m2[cur].max(t1a_area(t))]);
            cur = r;
            guard += 1;
        }
        for _ in 0..p.chaikin_passes { pts = chaikin4(&pts); }
        out.push(Curve { pts, code: GUIDE, link: t as u32, label: t as u32, site: None, carrier: None, mq: None });
    }
    (out, clp.map(|_| cls))
}

// ------------------------------------------------------------------ surface evaluation on a lattice

/// A square lattice of posts at `dx_mm` spacing, global post indices (gx, gy) = chart mm / dx_mm.
#[derive(Clone, Copy, Debug)]
pub struct EvalLat { pub gx0: i64, pub gy0: i64, pub n: usize, pub dx_mm: i64 }

impl EvalLat {
    #[inline] pub fn nodes(&self) -> usize { self.n * self.n }
    #[inline] pub fn g(&self, i: usize) -> (i64, i64) { ((i % self.n) as i64 + self.gx0, (i / self.n) as i64 + self.gy0) }
    #[inline] pub fn idx(&self, gx: i64, gy: i64) -> Option<usize> {
        let (x, y) = (gx - self.gx0, gy - self.gy0);
        if x < 0 || y < 0 || x >= self.n as i64 || y >= self.n as i64 { None } else { Some(y as usize * self.n + x as usize) }
    }
    #[inline] pub fn dx(&self) -> f64 { self.dx_mm as f64 / 1000.0 }
}

/// Static per-post data of a lattice.
pub struct LatPrep<'a> {
    pub lat: EvalLat,
    pub cell: Vec<u32>,
    pub bil: Vec<f32>,
    pub hat_c: Vec<[u32; 4]>,
    pub hat_w: Vec<[f32; 4]>,
    /// FREE (land), MARINE or LAKE.
    pub water: Vec<u8>,
    pub water_level: Vec<f32>,
    pub rock: RockPosts<'a>,
    /// Lane TF T5: the C1 (Catmull-Rom) T0 elevation per post (empty unless tf cone_c1).
    pub bil_c1: Vec<f32>,
    /// Lane TF T1: margin (posts) between the eval lattice ring and the published box (0 unless a window sets it).
    pub publish_margin_posts: i64,
    /// Lane TG G3: the box-averaged bilinear per post (empty unless tg smooth_bil).
    pub bil_s: Vec<f32>,
    /// Lane CLIMATE-RUNOFF: R / r_ref per post (empty unless the climate key is on).
    pub rw: Vec<f32>,
}

pub fn lat_prep<'a>(lat: EvalLat, chart: &Chart, t0: &T0, extras: &T0Extras, p: &'a RiverNetParamsV1, strata: Option<&'a StrataCols>, anchor: (f64, f64)) -> R<LatPrep<'a>> {
    let nn = lat.nodes();
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let (mut cell, mut bil, mut hat_c, mut hat_w, mut water, mut wl) = (Vec::with_capacity(nn), Vec::with_capacity(nn), Vec::with_capacity(nn), Vec::with_capacity(nn), Vec::with_capacity(nn), Vec::with_capacity(nn));
    let mut pos = Vec::with_capacity(nn);
    let tf_c1 = crate::t1b_tf::tf(p).map(|t| t.cone_c1 == 1).unwrap_or(false);
    let mut bil_c1: Vec<f32> = Vec::new();
    // Lane TG G3: B on the global grid over the lattice (the same values as the chart-wide grid).
    let sb = match crate::t1b_tg::tg(p).filter(|t| t.smooth_bil == 1) {
        Some(t) => {
            let (x0, y0) = ((lat.gx0 * lat.dx_mm) as f64, (lat.gy0 * lat.dx_mm) as f64);
            let (x1, y1) = (((lat.gx0 + lat.n as i64 - 1) * lat.dx_mm) as f64, ((lat.gy0 + lat.n as i64 - 1) * lat.dx_mm) as f64);
            Some(crate::t1b_tg::build_smooth(&crate::t1b_tg::raw_bilinear(t0, chart), x0, x1, y0, y1, t.smooth_grid_mm, t.smooth_half_mm)?)
        }
        None => None,
    };
    let mut bil_s: Vec<f32> = Vec::new();
    // Lane CLIMATE-RUNOFF: R / r_ref per post (key on only).
    let clim = if p.climate_runoff.is_some() { Some(t0.climate.as_ref().ok_or("climate_runoff params set but no climate attached to T0")?) } else { None };
    let mut rw: Vec<f32> = Vec::new();
    for i in 0..nn {
        let (gx, gy) = lat.g(i);
        let (la, lo) = chart.inverse_mm(gx * lat.dx_mm, gy * lat.dx_mm).map_err(|e| e.to_string())?;
        let (la, lo) = (la * RAD2DEG, lo * RAD2DEG);
        if let Some(c) = clim { rw.push(c.weight_at_deg(t0, la, lo) as f32); }
        let c = t0.cell_of_deg(la, lo);
        let b = bilinear(t0, &elev, la, lo, false);
        let (hc, hw) = hat4(t0, la, lo);
        cell.push(c);
        bil.push(b as f32);
        hat_c.push(hc);
        hat_w.push([hw[0] as f32, hw[1] as f32, hw[2] as f32, hw[3] as f32]);
        let shelf = p.coast.as_ref().map(|cp| cp.shelf_m()).unwrap_or(0.0);
        let bw = match &sb { Some(s) => { let v = s.at_mm((gx * lat.dx_mm) as f64, (gy * lat.dx_mm) as f64); bil_s.push(v as f32); v } None => b };
        match crate::t1b_net::water_zero_contour(t0, c, bw + shelf) {
            CellWater::Land => { water.push(FREE); wl.push(0.0); }
            CellWater::Marine => { water.push(MARINE); wl.push(t0.sea_level_m as f32); }
            CellWater::Lake => { water.push(LAKE); wl.push(t0.level(c) as f32); }
        }
        pos.push((la, lo, (gx * lat.dx_mm) as f64 / 1000.0, (gy * lat.dx_mm) as f64 / 1000.0));
    }
    let rock = RockPosts::build(p, t0, extras, strata, anchor, &pos)?;
    // Lane TF T5 (REVISION 1): the box-averaged bilinear.
    if tf_c1 { bil_c1 = crate::t1b_tf::smooth_bil(&lat, &bil, crate::t1b_tf::tf(p).unwrap().cone_smooth_mm as f64 / 1000.0); }
    Ok(LatPrep { lat, cell, bil, hat_c, hat_w, water, water_level: wl, rock, bil_c1, publish_margin_posts: 0, bil_s, rw })
}

/// One evaluated surface.
pub struct Surface {
    pub dependency: Option<Value>,
    /// Lane T1B-PEDESTAL (ped_solve / ped_base only; empty otherwise): heights before the coast / volcano stage.
    pub ground_pre: Vec<f64>,
    pub z: Vec<f64>,
    /// Before crest rounding and fill (the eikonal surface).
    pub z_eik: Vec<f64>,
    /// FREE, CHANNEL, MARINE, LAKE, GUIDE (network channel).
    pub fixed: Vec<u8>,
    /// Floor posts (published free).
    pub floor: Vec<bool>,
    /// Channel posts: (network area m2, link node: T1a node or terminal seed).
    pub chan_area: Vec<f32>,
    pub chan_link: Vec<u32>,
    pub e: Vec<f32>,
    pub stats: SurfStats,
    /// Lane CO: per-post over-code (t1b_coast::OVER_*; all 0 without coast/volcano) and coast stats.
    pub over: Vec<u8>,
    pub coast: crate::t1b_coast::CoastStats,
    /// Lane TD diagnostics (filled only when `DIAG` is set; never published or digested).
    pub diag: SurfDiag,
    /// Lane TD6: volcano index (into the coast context's volcano list) on edifice posts, 255 elsewhere.
    pub vid: Vec<u8>,
    /// Lane TG G1: signed shore distance (m) per post (empty unless tg shore).
    pub shore: Vec<f64>,
    /// Lane TH H2: height before the wave cut on cliff-class posts (m; NaN elsewhere; empty unless th uncut).
    pub uncut: Vec<f64>,
    /// Lane VE: edifice style and unit per post (empty unless ve).
    pub ve_style: Vec<u8>,
    pub ve_unit: Vec<u16>,
    /// Lane VE4: front_mm distances (m, NaN none; empty unless ve4 front).
    pub ve_front: Vec<f64>,
    /// Lane STEP (ct only; empty otherwise): regime, cliff fraction (permille) and massive talus top (mm) per post.
    pub ct_regime: Vec<u8>,
    pub ct_cliff: Vec<u16>,
    pub ct_top: Vec<i32>,
    /// Lane STEP2 (ct_base only; empty otherwise): the staircase shift per post (m; 0 where none).
    pub ct_shift: Vec<f64>,
    /// Lane TREAD2 (ct_tread only; empty otherwise): tread height fraction and tread tangent per post (permille).
    pub ct_tread: Vec<u16>,
    pub ct_tread_tan: Vec<u16>,
    /// Lane PERIOD-ENDS (ct_ends only; empty otherwise): per post [s0, s1, s_dn, s_up] (m; NaN none), remap gradient S,
    /// period kind (t1b_ct::PK_*) and var 3 bed id (t1b_ct::BED_ID_NOTE).
    pub ct_ends: Vec<[f64; 4]>,
    /// F knots [u1, f1, u2, f2] as fractions of the period (t1b_ct::knots2 / knots3).
    pub ct_knots: Vec<[f64; 4]>,
    pub ct_grad: Vec<f64>,
    pub ct_pkind: Vec<u8>,
    pub ct_bed: Vec<u32>,
    /// Lane SEAM-CERT (seam_cert only): the fill taint of the construction fills.
    pub seam: Option<crate::t1b_seam::SeamTaint>,
    /// Lane GLACIAL-T1 (glacial only; empty otherwise): per-post lowering by the window carve (m).
    pub gl_low: Vec<f32>,
    /// Lane MS-CEIL (ceil only; empty otherwise): the ceiling Z_c (m) and surface_src (CEIL_SRC_NOTE) per post.
    pub ceil_m: Vec<f32>,
    pub ceil_src: Vec<u8>,
}

/// Lane TD: set by the binary's diag_dump=1 to keep per-stage arrays for the defect diagnosis (not published).
pub static DIAG: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
/// Lane TD: where run_window writes the diagnostic crop (`<dir>/w_<rx>_<ry>/*.bin`) when `DIAG` is set.
pub static DIAG_DIR: std::sync::Mutex<Option<std::path::PathBuf>> = std::sync::Mutex::new(None);

/// Lane TD per-post diagnostics: source label, channel bed gradient, pass-2 eikonal speed at the post's own height,
/// surface after crest rounding, after the first fill, and the fixed codes before the coast steps.
#[derive(Default)]
pub struct SurfDiag { pub lh: Vec<f32>, pub speed1: Vec<f32>, pub lab: Vec<u32>, pub chan_slope: Vec<f32>, pub speed: Vec<f32>, pub z_crest: Vec<f64>, pub z_fill: Vec<f64>, pub fixed_pre: Vec<u8> }

#[derive(Clone, Debug, Default, Serialize)]
pub struct SurfStats {
    pub curves_in_lattice: u32,
    pub hill_length_mean_m: f64,
    pub channel_posts_t1a: u32,
    pub channel_posts_net: u32,
    pub floor_posts: u32,
    pub crest_lowered_posts: u32,
    pub crest_lowered_max_m: f64,
    pub fill_raised_posts: u32,
    pub fill_raised_max_m: f64,
    pub fmm_unreached: u32,
    pub estuary_posts: u32,
    pub bank_raised_posts: u32,
    pub bank_raised_max_m: f64,
    /// Lane TD3 (only with fx): floor posts lowered by a lower valley's hillslope.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub floor_lowered_posts: u32,
    /// Lane VALLEY-FLOOR-FIX (only with vf): levee / footslope / fan stats.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vf: Option<crate::t1b_vf::VfStats>,
    /// Lane GLACIAL-T1 (only with the glacial key): window carve statistics.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gl: Option<crate::t1b_glacial::GlStats>,
    /// Lane TA LS2 (only with ls bank_ramp): posts raised by the bank ramp and the largest raise.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub bank_ramp_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b_ls::is_zero_f64")]
    pub bank_ramp_max_m: f64,
    /// Lane TB soft_channel_cap_v1 (only with ls_steep=2): network channel posts lowered and the largest lowering (pass 2).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub soft_channel_lowered_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b_ls::is_zero_f64")]
    pub soft_channel_lowered_max_m: f64,
    /// Lane GUIDE-SPIKES (only with ls_undercut=1): network channel posts demoted to free hillslope by guide_undercut_v1 and
    /// the largest bed-minus-height of those posts.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub guide_undercut_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b_ls::is_zero_f64")]
    pub guide_undercut_max_m: f64,
    /// Lane TC C1 (only with ls_floor): network channel posts lowered to a T1a floor and the largest lowering.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub tc_floor_lowered_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b_ls::is_zero_f64")]
    pub tc_floor_lowered_max_m: f64,
    /// Lane TF T1 (only with tf report): fill guard and basin report of the pre-coast fill.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tf_fill: Option<crate::t1b_tf::FillReport>,
    /// Lane STEP (only with ct): free hillslope posts in the cliff and talus regimes at their eikonal height, and bare-rock
    /// posts the crest rounding would have lowered (exempt).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_cliff_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_talus_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_crest_exempt_posts: u32,
    /// Lane STEP apply 1: posts moved by z_ct - z_v0 and the largest |shift|.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_shifted_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b_ls::is_zero_f64")]
    pub ct_shift_max_abs_m: f64,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_shift_capped_posts: u32,
    /// Lane STEP2 apply 2: remapped bare-rock posts whose period came from a strata package / the closure (or regular P),
    /// and the mean shift F_c(h) - h over remapped posts (m).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_strata_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_closure_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b_ls::is_zero_f64")]
    pub ct_shift_mean_m: f64,
    /// Lane STEP2 local: posts bare by S but mantled by the local gradient (G <= talus tangent).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_local_mantled_posts: u32,
    /// Lane TREAD2 (ct_tread): posts on the tread limb; remapped layered posts with a 3-limb period, with the 2-limb fallback
    /// (S too steep for a bench), and 3-limb or fallback posts whose hard share came from the record / the closure.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_tread_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_tread_profile_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_tread_fallback_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_tread_record_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_tread_closure_posts: u32,
    /// Lane TREAD2: remapped posts held at their base level (t1b_ct::TREAD_BASE_RULE).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_tread_base_clamped_posts: u32,
    /// Lane SPLIT-T1B (vs_ct_blend): split posts whose height is a blend of riser-threshold ladder alternatives.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_volc_blend_posts: u32,
    /// Lane TREAD2: remapped posts raised to their lowest formerly-lower D8 neighbour (t1b_ct::TREAD_DRAIN_RULE).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_tread_drain_raised_posts: u32,
    /// Lane BENCH-VAR (ct var 3): remapped posts whose period is a keyed bed inside a thick package, and 3-limb / fallback
    /// posts whose hard share came from the keyed bed (t1b_ct::BED_RULE).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_bed_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_bed_eta_posts: u32,
    /// Lane PERIOD-ENDS (ct_ends): remapped posts with published period ends, and neighbour periods not contiguous with the
    /// own period (the period below does not end at s0 / the one above does not start at s1).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_ends_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_ends_gap_posts: u32,
    /// Lane PERIOD-ENDS: of those, apron posts whose period below overlaps (the un-aproned period [s_dn, s1] holds s0 - 1 um).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_ends_gap_apron_dn_posts: u32,
    /// Lane PERIOD-ENDS: posts with ends whose height changed after the remap (bank rule, fill, coast; flagged PK_HELD).
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ct_ends_moved_posts: u32,
    /// Lane MS-CEIL (only with ceil): free non-floor land posts lowered to the ceiling Z_c and the largest lowering.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ceil_posts: u32,
    #[serde(skip_serializing_if = "crate::t1b_ls::is_zero_f64")]
    pub ceil_lowered_max_m: f64,
    #[serde(skip_serializing)]
    pub seconds: [f64; 5],
}

#[derive(Clone, Copy)]
struct HK { z: f64, gy: i64, gx: i64, i: u32 }
impl PartialEq for HK { fn eq(&self, o: &Self) -> bool { self.cmp(o) == Ordering::Equal } }
impl Eq for HK {}
impl PartialOrd for HK { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for HK { fn cmp(&self, o: &Self) -> Ordering { o.z.total_cmp(&self.z).then(o.gy.cmp(&self.gy)).then(o.gx.cmp(&self.gx)) } }

/// First-order fast marching: `known` posts are sources with values in `z`; `speed(i, upwind z)` is |grad z|.
pub fn fmm(lat: &EvalLat, z: &mut [f64], known: &mut [bool], speed: &dyn Fn(usize, f64) -> f64) -> u32 {
    let mut lab = vec![0u32; lat.nodes()];
    fmm_lab(lat, z, known, speed, &mut lab)
}

/// Fast marching that also carries a source label: each accepted post takes the label of its lowest known 4-neighbour.
pub fn fmm_lab(lat: &EvalLat, z: &mut [f64], known: &mut [bool], speed: &dyn Fn(usize, f64) -> f64, lab: &mut [u32]) -> u32 {
    fmm_lab_soft(lat, z, known, speed, lab, None)
}

/// Lane TD3: as `fmm_lab`; `soft` posts (not known) start TENTATIVE at their `z` value with their `lab` label, so each
/// ends at min(its own value, the arrival from the others).
pub fn fmm_lab_soft(lat: &EvalLat, z: &mut [f64], known: &mut [bool], speed: &dyn Fn(usize, f64) -> f64, lab: &mut [u32], soft: Option<&[bool]>) -> u32 {
    fmm_lab_soft_acc(lat, z, known, speed, lab, soft, None)
}

/// As `fmm_lab_soft`; lane TB soft_channel_cap_v1: an arrival at a soft post j from a front labelled l is ignored unless
/// `accept(j, l)` (None: every arrival counts, the TD3 behaviour).
pub fn fmm_lab_soft_acc(lat: &EvalLat, z: &mut [f64], known: &mut [bool], speed: &dyn Fn(usize, f64) -> f64, lab: &mut [u32], soft: Option<&[bool]>, accept: Option<&dyn Fn(usize, u32) -> bool>) -> u32 {
    let n = lat.n;
    let nn = lat.nodes();
    let dx = lat.dx();
    let mut heap = BinaryHeap::new();
    let mut tent = vec![f64::INFINITY; nn];
    let mut tlab = vec![0u32; nn];
    let update = |j: usize, z: &[f64], known: &[bool], lab: &[u32], tent: &mut [f64], tlab: &mut [u32], heap: &mut BinaryHeap<HK>| {
        let (x, y) = (j % n, j / n);
        let mut best = (f64::INFINITY, 0u32);
        let mut kv = |k: usize| -> f64 { if known[k] { if z[k] < best.0 { best = (z[k], lab[k]); } z[k] } else { f64::INFINITY } };
        let a = (if x > 0 { kv(j - 1) } else { f64::INFINITY }).min(if x + 1 < n { kv(j + 1) } else { f64::INFINITY });
        let b = (if y > 0 { kv(j - n) } else { f64::INFINITY }).min(if y + 1 < n { kv(j + n) } else { f64::INFINITY });
        let lo = a.min(b);
        if !lo.is_finite() { return; }
        if let Some(acc) = accept { if !acc(j, best.1) { return; } }
        let h = speed(j, lo) * dx;
        let c = if !a.is_finite() || !b.is_finite() || (a - b).abs() >= h { lo + h } else { 0.5 * (a + b + (2.0 * h * h - (a - b) * (a - b)).sqrt()) };
        if c < tent[j] {
            tent[j] = c;
            tlab[j] = best.1;
            let (gx, gy) = lat.g(j);
            heap.push(HK { z: c, gy, gx, i: j as u32 });
        }
    };
    for j in 0..nn {
        if known[j] { continue; }
        let (x, y) = (j % n, j / n);
        let near = (x > 0 && known[j - 1]) || (x + 1 < n && known[j + 1]) || (y > 0 && known[j - n]) || (y + 1 < n && known[j + n]);
        if near { update(j, z, known, lab, &mut tent, &mut tlab, &mut heap); }
    }
    if let Some(sf) = soft {
        for j in 0..nn {
            if sf[j] && !known[j] && z[j].is_finite() && z[j] < tent[j] {
                tent[j] = z[j];
                tlab[j] = lab[j];
                let (gx, gy) = lat.g(j);
                heap.push(HK { z: z[j], gy, gx, i: j as u32 });
            }
        }
    }
    while let Some(HK { z: v, i, .. }) = heap.pop() {
        let j = i as usize;
        if known[j] || v > tent[j] { continue; }
        known[j] = true;
        z[j] = v;
        lab[j] = tlab[j];
        let (x, y) = (j % n, j / n);
        let mut nb = [usize::MAX; 4];
        if x > 0 { nb[0] = j - 1; }
        if x + 1 < n { nb[1] = j + 1; }
        if y > 0 { nb[2] = j - n; }
        if y + 1 < n { nb[3] = j + n; }
        for k in nb { if k != usize::MAX && !known[k] { update(k, z, known, lab, &mut tent, &mut tlab, &mut heap); } }
    }
    let mut unreached = 0;
    for j in 0..nn { if !known[j] { unreached += 1; z[j] = f64::NAN; } }
    unreached
}

/// Priority-flood fill: raise every post below its spill level. Sources: `outlet` posts and the lattice ring.
pub fn fill(lat: &EvalLat, z: &mut [f64], outlet: &[bool]) -> (u32, f64) {
    let n = lat.n;
    let nn = lat.nodes();
    let mut done = vec![false; nn];
    let mut heap = BinaryHeap::new();
    for i in 0..nn {
        let (x, y) = (i % n, i / n);
        if outlet[i] || x == 0 || y == 0 || x == n - 1 || y == n - 1 {
            done[i] = true;
            let (gx, gy) = lat.g(i);
            heap.push(HK { z: z[i], gy, gx, i: i as u32 });
        }
    }
    let (mut raised, mut maxr) = (0u32, 0f64);
    while let Some(HK { z: v, i, .. }) = heap.pop() {
        let i = i as usize;
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if done[j] { continue; }
            done[j] = true;
            if z[j] < v { raised += 1; maxr = maxr.max(v - z[j]); z[j] = v; }
            let (gx, gy) = lat.g(j);
            heap.push(HK { z: z[j], gy, gx, i: j as u32 });
        }
    }
    (raised, maxr)
}

fn known_src(fixed: &[u8], floor: &[bool], i: usize) -> bool { fixed[i] != FREE || floor[i] }

/// Separable box mean (x pass then y pass) with a metric half-width per post: weights w_k = clamp(r/dx + 1/2 - |k|, 0, 1)
/// (a box of width 2r on the post lattice, continuous in r, so lattices of different spacing see the same box), indices
/// clamped to the lattice; computed where `on(i)`, else the input value is returned.
pub fn box_mean(lat: &EvalLat, v: &[f64], on: &dyn Fn(usize) -> bool, r_m: &dyn Fn(usize) -> f64) -> Vec<f64> {
    let (n, dx) = (lat.n, lat.dx());
    let nn = lat.nodes();
    let kern = |r: f64| -> (i64, f64) { let k = (r / dx + 0.5).ceil() as i64; (k, r / dx + 0.5) };
    let mut bx = vec![0f64; nn];
    for i in 0..nn {
        let (k, h) = kern(r_m(i));
        if k <= 0 || h <= 1.0 { bx[i] = v[i]; continue; }
        let (x, y) = ((i % n) as i64, i / n);
        let (mut s, mut ws) = (0.0, 0.0);
        for o in -k..=k { let w = (h - o.abs() as f64).clamp(0.0, 1.0); if w <= 0.0 { continue; } s += w * v[y * n + (x + o).clamp(0, n as i64 - 1) as usize]; ws += w; }
        bx[i] = s / ws;
    }
    let mut out = v.to_vec();
    for i in 0..nn {
        if !on(i) { continue; }
        let (k, h) = kern(r_m(i));
        if k <= 0 || h <= 1.0 { continue; }
        let (x, y) = (i % n, (i / n) as i64);
        let (mut s, mut ws) = (0.0, 0.0);
        for o in -k..=k { let w = (h - o.abs() as f64).clamp(0.0, 1.0); if w <= 0.0 { continue; } s += w * bx[(y + o).clamp(0, n as i64 - 1) as usize * n + x]; ws += w; }
        out[i] = s / ws;
    }
    out
}

/// Evaluate the surface on a prepared lattice for an E table.
#[allow(clippy::too_many_arguments)]
pub fn evaluate(lp: &LatPrep, net: &Network, curves: &[Curve], t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, p: &RiverNetParamsV1, cx: &crate::t1b_coast::CoastCtx) -> Surface {
    evaluate_owned(lp, net, curves, t0, laws, table, p, cx, None, None).and_then(crate::t1b_surface_stages::Evaluated::surface).expect("legacy evaluation without owner cannot yield a stage error")
}

pub fn evaluate_owned(lp: &LatPrep, net: &Network, curves: &[Curve], t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>, p: &RiverNetParamsV1, cx: &crate::t1b_coast::CoastCtx, owner: Option<&dyn crate::t1b_surface_stages::Conditioner>, stop: Option<crate::t1b_surface_stages::Stage>) -> R<crate::t1b_surface_stages::Evaluated> {
    let lat = lp.lat;
    let (n, nn, dx) = (lat.n, lat.nodes(), lat.dx());
    let mut st = SurfStats::default();
    let t = std::time::Instant::now();
    let hat = |i: usize, tab: &BTreeMap<u32, f64>| -> f64 {
        let mut v = 0.0f64;
        for k in 0..4 { let w = lp.hat_w[i][k]; if w > 0.0 { v += w as f64 * tab.get(&lp.hat_c[i][k]).copied().unwrap_or(0.0); } }
        v
    };
    // E per post (lane CO: tapered by the plain bilinear when the coast is on).
    let bil_w: &[f32] = if lp.bil_s.is_empty() { &lp.bil } else { &lp.bil_s };
    let e: Vec<f32> = (0..nn).map(|i| (hat(i, table) * p.coast.as_ref().map(|c| c.tau(bil_w[i] as f64)).unwrap_or(1.0)) as f32).collect();
    // Rasterize channels.
    let mut fixed = lp.water.clone();
    let mut z = vec![f64::NAN; nn];
    for i in 0..nn { if fixed[i] != FREE { z[i] = lp.water_level[i] as f64; } }
    // Lane CO: construction water at the base level sea - F.
    if p.coast.is_some() { for i in 0..nn { if fixed[i] == MARINE { z[i] = t0.sea_level_m - hat(i, &cx.flood); } } }
    let mut chan_area = vec![0f32; nn];
    let mut chan_link = vec![NONE; nn];
    let mut chan_slope = vec![0f32; nn];
    let mut lab = vec![u32::MAX; nn];
    // Lane CO (B1): every post on a T1a channel curve, over land or water.
    let mut t1a_post = vec![false; if p.coast.is_some() || p.volcano.is_some() { nn } else { 0 }];
    let (lx0, ly0) = ((lat.gx0 - 1) as f64 * dx, (lat.gy0 - 1) as f64 * dx);
    let (lx1, ly1) = ((lat.gx0 + n as i64) as f64 * dx, (lat.gy0 + n as i64) as f64 * dx);
    // Lane CHANNEL-CENTRELINE (cl_mn=1): posts drawn by a meandered T1a link (their floor extent is the carrier band).
    let mn_on = curves.iter().any(|c| c.carrier.is_some() || c.mq.is_some());
    // Lane MEANDER-Q (cl_mq=1): the swept-envelope floodplain replaces the carrier band; the floor stands the inset depth above the bed.
    let mq_env = curves.iter().find_map(|c| c.mq.clone());
    let mq_vfd = crate::t1b_vf::VfParamsV1::off();
    let mq_hyd = p.vf.as_ref().unwrap_or(&mq_vfd);
    let cl_on = crate::t1b_cl::cl(p).is_some();
    let mut mn_post = vec![false; if mn_on { nn } else { 0 }];
    for c in curves {
        let (mut bx0, mut bx1, mut by0, mut by1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for q in &c.pts { bx0 = bx0.min(q[0]); bx1 = bx1.max(q[0]); by0 = by0.min(q[1]); by1 = by1.max(q[1]); }
        if bx1 < lx0 || bx0 > lx1 || by1 < ly0 || by0 > ly1 { continue; }
        st.curves_in_lattice += 1;
        for w in c.pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            let len = (b[0] - a[0]).hypot(b[1] - a[1]);
            let ns = ((len / dx * 4.0).ceil() as usize).max(1);
            let sb = if len > 1e-9 { ((b[2] - a[2]).abs() / len) as f32 } else { 0.0 };
            for s in 0..=ns {
                let u = s as f64 / ns as f64;
                let (x, y) = (a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1]));
                let Some(i) = lat.idx((x / dx).round() as i64, (y / dx).round() as i64) else { continue };
                if c.code == CHANNEL && !t1a_post.is_empty() { t1a_post[i] = true; }
                // Lane CHANNEL-CENTRELINE: a T1a link over construction water records its link on the water post (F1 reads it).
                if lp.water[i] != FREE { if cl_on && c.code == CHANNEL && chan_link[i] == NONE { chan_link[i] = c.link; } continue; }
                let bed = a[2] + u * (b[2] - a[2]);
                let area = a[3] + u * (b[3] - a[3]);
                let take = match fixed[i] {
                    FREE => true,
                    CHANNEL => c.code == CHANNEL && bed < z[i],
                    _ => c.code == CHANNEL || bed < z[i],
                };
                if take { fixed[i] = c.code; z[i] = bed; chan_area[i] = area as f32; chan_link[i] = c.link; chan_slope[i] = sb; lab[i] = c.label; if mn_on { mn_post[i] = c.carrier.is_some() || c.mq.is_some(); } }
            }
        }
    }
    // Estuary rule: a channel bed at or below sea level is an arm of the sea.
    let sea = t0.sea_level_m;
    if p.coast.is_none() { for i in 0..nn { if (fixed[i] == CHANNEL || fixed[i] == GUIDE) && z[i] <= sea { fixed[i] = MARINE; z[i] = sea; st.estuary_posts += 1; } } }
    for i in 0..nn { match fixed[i] { CHANNEL => st.channel_posts_t1a += 1, GUIDE => st.channel_posts_net += 1, _ => {} } }
    st.seconds[0] = t.elapsed().as_secs_f64();
    // Valley floors: nearest channel post within W/2.
    let t = std::time::Instant::now();
    let cw = p.floor_width_coeff_ppm as f64 / 1e6;
    let hmax = p.floor_width_max_mm as f64 / 2000.0;
    let fs = p.floor_slope_ppm as f64 / 1e6;
    // Lane VALLEY-FLOOR-FIX (t1b_vf.rs): levee-backswamp floor / footslope / tributary fans (all off without params.vf).
    let vfp = p.vf.as_ref();
    let vf_floor = vfp.map(|v| v.floor_on()).unwrap_or(false);
    let vf_fan = vfp.map(|v| v.fan_on()).unwrap_or(false);
    let mut vfs = crate::t1b_vf::VfStats::default();
    // Lane CHANNEL-CENTRELINE (cl_mn=1): the floor extent is the W/2 band of the unmeandered T1a curves (carriers), so the
    // valley walls do not follow the declared bends; a meandered channel post reaches floor posts up to W/2 + its belt.
    let (mn_ext, mn_pad): (Vec<bool>, f64) = if let Some(env) = mq_env.as_ref() { crate::t1b_mq::env_ext(env, &lat, dx) } else if mn_on {
        let mut carea = vec![0f32; nn];
        let mut pad = 0.0f64;
        for c in curves {
            let Some((cp, belt)) = &c.carrier else { continue };
            pad = pad.max(*belt);
            for w in cp.windows(2) {
                let (a, b) = (w[0], w[1]);
                let len = (b[0] - a[0]).hypot(b[1] - a[1]);
                let ns = ((len / dx * 4.0).ceil() as usize).max(1);
                for q in 0..=ns {
                    let u = q as f64 / ns as f64;
                    let (x, y) = (a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1]));
                    let Some(i) = lat.idx((x / dx).round() as i64, (y / dx).round() as i64) else { continue };
                    carea[i] = carea[i].max((a[3] + u * (b[3] - a[3])) as f32);
                }
            }
        }
        let mut ext = vec![false; nn];
        for i in 0..nn {
            if carea[i] <= 0.0 { continue; }
            let h = (0.5 * cw * (carea[i] as f64).sqrt()).min(hmax);
            let r = (h / dx).floor() as i64;
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            for oy in -r..=r { for ox in -r..=r {
                if ((ox * ox + oy * oy) as f64).sqrt() * dx > h { continue; }
                let (xx, yy) = (x + ox, y + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                ext[yy as usize * n + xx as usize] = true;
            } }
        }
        (ext, pad)
    } else { (Vec::new(), 0.0) };
    // (reach (m), needs the carrier band) of a channel post's floor disc.
    let reach = |i: usize, h: f64| -> (f64, bool) { if mn_on && mn_post[i] { (h + mn_pad, true) } else { (h, false) } };
    // Lane TC C1 (t1a_floor_grade_v1): network channel posts inside a T1a floor take min(own bed, the T1a floor value).
    if crate::t1b_ls::tc_floor(p) >= 1 {
        let mut best: Vec<(f64, i64, i64)> = vec![(f64::INFINITY, 0, 0); nn];
        let mut val = vec![f64::NAN; nn];
        for i in 0..nn {
            if fixed[i] != CHANNEL { continue; }
            let h = (0.5 * cw * (chan_area[i] as f64).sqrt()).min(hmax);
            if h < 0.5 * dx { continue; }
            let (hr, need_ext) = reach(i, h);
            let r = (hr / dx).floor() as i64;
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            let (gx, gy) = lat.g(i);
            for oy in -r..=r { for ox in -r..=r {
                let d = ((ox * ox + oy * oy) as f64).sqrt() * dx;
                if d > hr || (ox == 0 && oy == 0) { continue; }
                let (xx, yy) = (x + ox, y + oy);
                if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                let j = yy as usize * n + xx as usize;
                if need_ext && !mn_ext[j] { continue; }
                if fixed[j] != GUIDE { continue; }
                let cand = (d, gy, gx);
                let cur = best[j];
                if cand.0 < cur.0 || (cand.0 == cur.0 && (cand.1, cand.2) < (cur.1, cur.2)) { best[j] = cand; val[j] = z[i] + if vf_floor { vfp.unwrap().tc_rel(chan_slope[i] as f64, d, dx) } else { fs.max(chan_slope[i] as f64) * d }; }
            } }
        }
        for j in 0..nn { if val[j].is_finite() && val[j] < z[j] { st.tc_floor_lowered_posts += 1; st.tc_floor_lowered_max_m = st.tc_floor_lowered_max_m.max(z[j] - val[j]); z[j] = val[j]; } }
    }
    let mut floor_best: Vec<(f64, i64, i64)> = vec![(f64::INFINITY, 0, 0); nn];
    let mut floor_val = vec![f64::NAN; nn];
    let mut floor_lab = vec![u32::MAX; nn];
    // Lane VALLEY-FLOOR-FIX: the nearest_channel_floor_v1 value of the same nearest post (footslope reference, stats) and,
    // for fans, the covering channel post of largest area per post (area, distance, floor value; channel posts included).
    let mut floor_old: Vec<f64> = if vf_floor { vec![f64::NAN; nn] } else { Vec::new() };
    // footslope: the floor half-width of each floor post's nearest channel post
    let vf_foot = vfp.map(|v| v.foot_on()).unwrap_or(false);
    let mut floor_h: Vec<f64> = if vf_foot { vec![0.0; nn] } else { Vec::new() };
    let (mut tr_area, mut tr_d, mut tr_val): (Vec<f64>, Vec<f64>, Vec<f64>) = if vf_fan { (vec![0.0; nn], vec![f64::INFINITY; nn], vec![f64::NAN; nn]) } else { (Vec::new(), Vec::new(), Vec::new()) };
    for i in 0..nn {
        if fixed[i] != CHANNEL && fixed[i] != GUIDE { continue; }
        let h = (0.5 * cw * (chan_area[i] as f64).sqrt()).min(hmax);
        if h < 0.5 * dx { continue; }
        let (hr, need_ext) = reach(i, h);
        let r = (hr / dx).floor() as i64;
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let (gx, gy) = lat.g(i);
        for oy in -r..=r { for ox in -r..=r {
            let d = ((ox * ox + oy * oy) as f64).sqrt() * dx;
            if d > hr || (ox == 0 && oy == 0) { continue; }
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if need_ext && !mn_ext[j] { continue; }
            let rel_old = fs.max(chan_slope[i] as f64) * d;
            let rel = if vf_floor { vfp.unwrap().floor_rel(chan_area[i] as f64, chan_slope[i] as f64, d, dx) } else { rel_old };
            let rel = if need_ext && mq_env.as_ref().map(|e| e.inset).unwrap_or(false) { rel + mq_hyd.hydraulic(chan_area[i] as f64).1 } else { rel };
            if vf_fan {
                let a = chan_area[i] as f64;
                if a > tr_area[j] || (a == tr_area[j] && d < tr_d[j]) { tr_area[j] = a; tr_d[j] = d; tr_val[j] = z[i] + rel; }
            }
            if fixed[j] != FREE { continue; }
            let cand = (d, gy, gx);
            let cur = floor_best[j];
            if cand.0 < cur.0 || (cand.0 == cur.0 && (cand.1, cand.2) < (cur.1, cur.2)) { floor_best[j] = cand; floor_val[j] = z[i] + rel; floor_lab[j] = lab[i]; if vf_floor { floor_old[j] = z[i] + rel_old; } if vf_foot { floor_h[j] = h; } }
        } }
    }
    if vf_floor {
        let (mut c, mut sd) = (0u32, 0.0f64);
        for i in 0..nn { if fixed[i] == FREE && floor_val[i].is_finite() { sd += floor_val[i] - floor_old[i]; if floor_val[i] - floor_old[i] > 0.1 { vfs.floor_posts_levee_gt_0p1m += 1; } c += 1; } }
        vfs.floor_mean_change_m = if c > 0 { sd / c as f64 } else { 0.0 };
    }
    let mut floor = vec![false; nn];
    let mut known = vec![false; nn];
    for i in 0..nn {
        if fixed[i] != FREE { known[i] = true; }
        else if floor_val[i].is_finite() { floor[i] = true; known[i] = true; z[i] = floor_val[i]; lab[i] = floor_lab[i]; st.floor_posts += 1; }
    }
    st.seconds[1] = t.elapsed().as_secs_f64();
    // Hillslopes: eikonal.
    let t = std::time::Instant::now();
    let d_m2 = laws.hillslope_diffusivity_m2_per_my;
    let sc_ref = p.sc_ref_permille as f64 / 1000.0;
    let lref = p.hill_ref_length_mm as f64 / 1000.0;
    let smin = p.hill_slope_min_ppm as f64 / 1e6;
    let fr: Vec<f64> = e.iter().map(|&ev| {
        let beta = 2.0 * ev as f64 * lref / (d_m2 * sc_ref);
        if beta < 1e-9 { 0.0 } else { ((1.0 + beta * beta).sqrt() - 1.0) / beta }
    }).collect();
    // Lane TD flags (all off without fx: 0.9.2 code paths).
    let fx_on = |f: fn(&crate::t1b_fix::FxParamsV1) -> u32| p.fx.as_ref().map(|x| f(x) == 1).unwrap_or(false);
    let (td_lh, td_floor, td_contact) = (fx_on(|x| x.lh_free_mean), fx_on(|x| x.soft_floor), fx_on(|x| x.contact_step));
    let speed = |i: usize, zu: f64| -> f64 {
        if td_contact { return crate::t1b_fix::contact_speed(&lp.rock, i, zu, dx, fr[i], smin); }
        let (sc, _, _) = lp.rock.at(i, zu); (sc * fr[i]).max(smin)
    };
    let known0 = known.clone();
    let z0 = z.clone();
    let lab0 = lab.clone();
    // Lane TB soft_channel_cap_v1 (ls_steep=2): network channel posts are tentative sources too, lowered only by a front of
    // their own basin (same terminal); T1a channels stay hard.
    let soft_ch = td_floor && crate::t1b_ls::soft_channels(p);
    let zbed: Vec<f64> = if soft_ch { z0.clone() } else { Vec::new() };
    let softch: Vec<bool> = if soft_ch { (0..nn).map(|i| fixed[i] == GUIDE).collect() } else { Vec::new() };
    let soft_set: Vec<bool> = if soft_ch { (0..nn).map(|i| floor[i] || softch[i]).collect() } else { Vec::new() };
    let term_of = |l: u32| -> u32 { if l == u32::MAX { u32::MAX } else if net.kind[l as usize] != K_LAND { l } else { net.term[l as usize] } };
    let own_term: Vec<u32> = if soft_ch { (0..nn).map(|i| if softch[i] { term_of(lab[i]) } else { u32::MAX }).collect() } else { Vec::new() };
    // Lane GUIDE-SPIKES guide_undercut_v1 (ls_undercut=1): network channel posts accept every front (demoted after pass 2).
    let undercut = soft_ch && crate::t1b_ls::undercut(p);
    let accept = |j: usize, l: u32| -> bool { !softch[j] || undercut || (l != u32::MAX && own_term[j] != u32::MAX && term_of(l) == own_term[j]) };
    let acc: Option<&dyn Fn(usize, u32) -> bool> = if soft_ch { Some(&accept) } else { None };
    // TD3: floors are tentative sources (min of the floor value and the hillslope arrival).
    let hard: Vec<bool> = if td_floor { (0..nn).map(|i| fixed[i] != FREE && !(soft_ch && softch[i])).collect() } else { Vec::new() };
    let dependency = crate::t1b_esc_deps::surface_inputs(lp, net, curves, table, p, cx, &e, &z, &fixed, &floor, &known, &lab, &chan_area, &chan_slope)
        .unwrap_or_else(|err| panic!("surface dependency certification: {err}"));
    if td_floor {
        let mut kf = hard.clone();
        st.fmm_unreached = if soft_ch { fmm_lab_soft_acc(&lat, &mut z, &mut kf, &speed, &mut lab, Some(&soft_set), acc) } else { fmm_lab_soft(&lat, &mut z, &mut kf, &speed, &mut lab, Some(&floor)) };
    } else {
        st.fmm_unreached = fmm_lab(&lat, &mut z, &mut known, &speed, &mut lab);
    }
    let rmk: Vec<f64> = (0..nn).map(|i| { let l = lab[i]; if l == u32::MAX { 1.0 } else { net.rockmass[l as usize] } }).collect();
    // REVISION 1: second pass with the steady-state gradient at each post's own hillslope length L_H = flow length up
    // to the divide + flow length down to the channel or floor, on the first-pass surface (D8 steepest descent), box
    // smoothed over 64 m; one constant gradient per hillslope length gave a 1-degree slope spike (run s42_r0).
    let mut rcv = vec![u32::MAX; nn];
    let mut rd = vec![0f64; nn];
    let mut ord: Vec<u32> = (0..nn as u32).filter(|&i| !known0[i as usize] && z[i as usize].is_finite()).collect();
    for &i in &ord {
        let i = i as usize;
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let mut best = (0.0f64, u32::MAX, 0.0f64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            let d = if ox != 0 && oy != 0 { dx * std::f64::consts::SQRT_2 } else { dx };
            let g = (z[i] - z[j]) / d;
            if g > best.0 { best = (g, j as u32, d); }
        }
        rcv[i] = best.1;
        rd[i] = best.2;
    }
    let key = |i: usize| { let (gx, gy) = lat.g(i); (gy, gx) };
    ord.sort_by(|&a, &b| z[a as usize].total_cmp(&z[b as usize]).then(key(a as usize).cmp(&key(b as usize))));
    let mut ldn = vec![0f64; nn];
    for &i in &ord { let i = i as usize; let r = rcv[i]; if r != u32::MAX { ldn[i] = ldn[r as usize] + rd[i]; } }
    let mut lup = vec![0f64; nn];
    for &i in ord.iter().rev() { let i = i as usize; let r = rcv[i]; if r != u32::MAX && !known0[r as usize] { let v = lup[i] + rd[i]; if v > lup[r as usize] { lup[r as usize] = v; } } }
    let lh_raw: Vec<f64> = (0..nn).map(|i| if known0[i] { 0.0 } else { (lup[i] + ldn[i]).max(dx) }).collect();
    // Lane STEP (ct): hillslope base (source height reached by D8 descent) and crest (highest post draining through) on the
    // pass-1 surface.
    let ctp = crate::t1b_ct::ct(p);
    let geo: Vec<crate::t1b_ct::Geo> = if ctp.is_some() {
        let mut base = vec![f64::NAN; nn];
        let mut crest = vec![f64::NAN; nn];
        for &i in &ord { let i = i as usize; let r = rcv[i]; if r != u32::MAX { let r = r as usize; base[i] = if known0[r] { z[r] } else { base[r] }; } crest[i] = z[i]; }
        for &i in ord.iter().rev() { let i = i as usize; let r = rcv[i]; if r != u32::MAX && !known0[r as usize] { let r = r as usize; if crest[i] > crest[r] { crest[r] = crest[i]; } } }
        (0..nn).map(|i| crate::t1b_ct::Geo { base: base[i], crest: crest[i] }).collect()
    } else { Vec::new() };
    let lh = if td_lh {
        // TD1: the box mean over hillslope posts only.
        let w: Vec<f64> = (0..nn).map(|i| if known0[i] { 0.0 } else { 1.0 }).collect();
        let vw: Vec<f64> = (0..nn).map(|i| lh_raw[i] * w[i]).collect();
        let a = box_mean(&lat, &vw, &|i| !known0[i], &|_| 64.0);
        let b = box_mean(&lat, &w, &|i| !known0[i], &|_| 64.0);
        // B4 (TD3): a soft floor post that the eikonal reaches as hillslope needs a hillslope length too: the same mean
        // over the hillslope posts within 64 m, else L_ref (without it its speed was hill_slope_min and a lowered floor
        // spread almost flat from the lower valley).
        let sfl = |i: usize| floor[i] || (soft_ch && softch[i]);
        let hb: Vec<f64> = if td_floor { box_mean(&lat, &vw, &|i| sfl(i), &|_| 64.0) } else { Vec::new() };
        let hw: Vec<f64> = if td_floor { box_mean(&lat, &w, &|i| sfl(i), &|_| 64.0) } else { Vec::new() };
        (0..nn).map(|i| if !known0[i] && b[i] > 1e-12 { a[i] / b[i] } else if td_floor && sfl(i) { if hw[i] > 1e-12 { hb[i] / hw[i] } else { lref } } else { lh_raw[i] }).collect::<Vec<f64>>()
    } else { box_mean(&lat, &lh_raw, &|i| !known0[i], &|_| 64.0) };
    let fr2: Vec<f64> = (0..nn).map(|i| {
        let beta = 2.0 * e[i] as f64 * lh[i] / (d_m2 * sc_ref);
        if beta < 1e-9 { 0.0 } else { ((1.0 + beta * beta).sqrt() - 1.0) / beta }
    }).collect();
    // Lane STEP: ct drives this pass only with apply 0; with apply 1 a second pass below switches it on.
    let ct_live = std::cell::Cell::new(ctp.map(|c| c.apply == 0).unwrap_or(false));
    let ct_sup = ctp.map(|c| c.apply == 1).unwrap_or(false);
    let (z0_ct, lab0_ct, floor_ct, known0_ct) = if ct_sup { (z0.clone(), lab0.clone(), floor.clone(), known0.clone()) } else { (Vec::new(), Vec::new(), Vec::new(), Vec::new()) };
    let speed2 = |i: usize, zu: f64| -> f64 {
        if ct_live.get() { if let Some(c) = ctp { if let Some(v) = crate::t1b_ct::speed(c, &lp.rock, i, zu, dx, rmk[i] * fr2[i], rmk[i], smin, geo[i]) { return v; } } }
        if td_contact { return crate::t1b_fix::contact_speed(&lp.rock, i, zu, dx, rmk[i] * fr2[i], smin); }
        let (sc, _, _) = lp.rock.at(i, zu); (sc * rmk[i] * fr2[i]).max(smin)
    };
    z = z0;
    let mut lab2 = lab0;
    if td_floor {
        let mut kf = hard.clone();
        st.fmm_unreached = if soft_ch { fmm_lab_soft_acc(&lat, &mut z, &mut kf, &speed2, &mut lab2, Some(&soft_set), acc) } else { fmm_lab_soft(&lat, &mut z, &mut kf, &speed2, &mut lab2, Some(&floor)) };
        // Floors lowered by a lower valley's slope are hillslope posts from here on.
        for i in 0..nn { if floor[i] && z[i] < floor_val[i] - 1e-9 { floor[i] = false; st.floor_lowered_posts += 1; } }
        // Lane GUIDE-SPIKES guide_undercut_v1: a network channel post another basin's front lowered below its bed is a free
        // hillslope post from here on (no channel area or link).
        if undercut {
            for i in 0..nn {
                if softch[i] && z[i] < zbed[i] - 1e-9 && term_of(lab2[i]) != own_term[i] {
                    fixed[i] = FREE; chan_area[i] = 0.0; chan_link[i] = NONE;
                    st.guide_undercut_posts += 1; st.guide_undercut_max_m = st.guide_undercut_max_m.max(zbed[i] - z[i]);
                }
            }
        }
        // Lane TB: network channel posts the soft cap lowered (they stay channels).
        if soft_ch { for i in 0..nn { if softch[i] && fixed[i] == GUIDE && z[i] < zbed[i] - 1e-9 { st.soft_channel_lowered_posts += 1; st.soft_channel_lowered_max_m = st.soft_channel_lowered_max_m.max(zbed[i] - z[i]); } } }
    } else {
        known = known0;
        st.fmm_unreached = fmm_lab(&lat, &mut z, &mut known, &speed2, &mut lab2);
    }
    st.hill_length_mean_m = { let v: Vec<f64> = (0..nn).filter(|&i| !known_src(&fixed, &floor, i)).map(|i| lh[i]).collect(); v.iter().sum::<f64>() / v.len().max(1) as f64 };
    // Lane STEP apply 1: the ct pass-2 eikonal from the same sources (floors as before the v0 pass); dz = z_ct - z_v0.
    let mut z_ct: Vec<f64> = Vec::new();
    if ct_sup {
        ct_live.set(true);
        let mut zc = z0_ct;
        let mut lc = lab0_ct;
        if td_floor {
            let mut kf = hard.clone();
            if soft_ch { fmm_lab_soft_acc(&lat, &mut zc, &mut kf, &speed2, &mut lc, Some(&soft_set), acc) } else { fmm_lab_soft(&lat, &mut zc, &mut kf, &speed2, &mut lc, Some(&floor_ct)) };
        } else {
            let mut kn = known0_ct;
            fmm_lab(&lat, &mut zc, &mut kn, &speed2, &mut lc);
        }
        ct_live.set(false);
        z_ct = zc;
    }
    let z_eik = z.clone();
    // Lane STEP (ct): regime, cliff fraction and massive talus top at each free hillslope post's eikonal height.
    let (mut ct_regime, mut ct_cliff, mut ct_top): (Vec<u8>, Vec<u16>, Vec<i32>) = (Vec::new(), Vec::new(), Vec::new());
    if let Some(c) = ctp {
        ct_regime = vec![crate::t1b_ct::RG_SOIL; nn];
        ct_cliff = vec![0u16; nn];
        ct_top = vec![crate::t1b_ct::TOP_NONE; nn];
        for i in 0..nn {
            if c.apply == 2 { break; }
            let zr = if ct_sup { z_ct[i] } else { z[i] };
            if fixed[i] != FREE || floor[i] || !zr.is_finite() { continue; }
            let (rg, cf, tp) = crate::t1b_ct::regime(c, &lp.rock, i, zr, rmk[i] * fr2[i], rmk[i], smin, geo[i]);
            if rg == crate::t1b_ct::RG_SOIL { continue; }
            ct_regime[i] = rg; ct_cliff[i] = cf; ct_top[i] = tp;
            if rg == crate::t1b_ct::RG_CLIFF { st.ct_cliff_posts += 1; } else { st.ct_talus_posts += 1; }
        }
    }
    let ct_exempt = ctp.map(|c| c.apply == 0 && c.crest_exempt == 1).unwrap_or(false);
    // Lane STEP2: the staircase shift per post for ct_base_mm.
    let ct_base = ctp.and_then(|c| c.v2.as_ref()).map(|v| v.base == 1).unwrap_or(false);
    let mut ct_shift: Vec<f64> = if ct_base { vec![0.0; nn] } else { Vec::new() };
    let (mut ct_tread, mut ct_tread_tan): (Vec<u16>, Vec<u16>) = (Vec::new(), Vec::new());
    let (mut ct_ends, mut ct_knots, mut ct_grad, mut ct_pkind, mut ct_bed): (Vec<[f64; 4]>, Vec<[f64; 4]>, Vec<f64>, Vec<u8>, Vec<u32>) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut diag = SurfDiag::default();
    if DIAG.load(std::sync::atomic::Ordering::Relaxed) {
        diag.lab = lab2.clone();
        diag.chan_slope = chan_slope.clone();
        diag.speed = (0..nn).map(|i| if z[i].is_finite() { speed2(i, z[i]) as f32 } else { f32::NAN }).collect();
        diag.speed1 = (0..nn).map(|i| if z[i].is_finite() { let (sc, _, _) = lp.rock.at(i, z[i]); (sc * rmk[i] * fr[i]).max(smin) as f32 } else { f32::NAN }).collect();
        diag.lh = lh.iter().map(|v| *v as f32).collect();
    }
    st.seconds[2] = t.elapsed().as_secs_f64();
    // Crest rounding.
    let t = std::time::Instant::now();
    let rmax = p.crest_radius_max_mm as f64 / 1000.0;
    let rr: Vec<f64> = e.iter().map(|&ev| if ev > 0.0 { (d_m2 * sc_ref / ev as f64).min(rmax) } else { rmax }).collect();
    let src = |i: usize| fixed[i] != FREE || floor[i];
    let b = box_mean(&lat, &z_eik, &|i| !src(i), &|i| rr[i]);
    for i in 0..nn {
        if src(i) { continue; }
        if ct_exempt && ct_regime[i] != crate::t1b_ct::RG_SOIL { if b[i] < z[i] { st.ct_crest_exempt_posts += 1; } continue; }
        if b[i] < z[i] { st.crest_lowered_posts += 1; st.crest_lowered_max_m = st.crest_lowered_max_m.max(z[i] - b[i]); z[i] = b[i]; }
    }
    // Lane STEP apply 1: the staircase rides on the rounded v0 surface.
    if ct_sup {
        let cap = ctp.map(|c| c.shift_cap_permille as f64 / 1000.0 * c.step_mm as f64 / 1000.0).unwrap_or(0.0);
        for i in 0..nn {
            if src(i) || !z[i].is_finite() || !z_ct[i].is_finite() { continue; }
            let mut d = z_ct[i] - z_eik[i];
            if cap > 0.0 && d.abs() > cap { st.ct_shift_capped_posts += 1; d = d.clamp(-cap, cap); }
            if d != 0.0 { st.ct_shifted_posts += 1; st.ct_shift_max_abs_m = st.ct_shift_max_abs_m.max(d.abs()); }
            z[i] += d;
            if ct_base { ct_shift[i] = d; }
        }
    }
    // Lane STEP2 apply 2: the monotone height remap F_c of the rounded v0 surface (no second eikonal).
    if let Some(c) = ctp.filter(|c| c.apply == 2) {
        let tab: &[f64] = &cx.ct_steps;
        let var = c.v2.as_ref().map(|v| v.var).unwrap_or(0);
        let bedp = c.v2.as_ref().and_then(|v| v.bed.as_ref());
        let step_m = c.step_mm as f64 / 1000.0;
        let umax_m = c.v2.as_ref().map(|v| v.unit_max_permille as f64 / 1000.0 * step_m).unwrap_or(0.0);
        assert!(var == 0 || !tab.is_empty(), "ct var >= 1 needs the keyed closure table (CoInputs::ctx)");
        let mut pk: BTreeMap<u32, Vec<(f64, f64)>> = BTreeMap::new();
        let mut sum = 0.0f64;
        // ct_remap_local_v1: the local gradient of the rounded surface (before any remap).
        let local = c.v2.as_ref().map(|v| v.local == 1).unwrap_or(false);
        // Lane TREAD2: the 3-limb profile and the uphill direction of the rounded surface (for the dip term).
        let tread = c.v2.as_ref().and_then(|v| v.tread.as_ref());
        let gv: Vec<(f64, f64)> = if tread.map(|t| t.dip == 1).unwrap_or(false) {
            let zr = &z;
            let d1 = |i: usize, a: Option<usize>, b: Option<usize>| -> f64 {
                let za = a.map(|j| zr[j]).filter(|v| v.is_finite());
                let zb = b.map(|j| zr[j]).filter(|v| v.is_finite());
                match (za, zb) { (Some(p), Some(q)) => (q - p) / (2.0 * dx), (None, Some(q)) => (q - zr[i]) / dx, (Some(p), None) => (zr[i] - p) / dx, _ => 0.0 }
            };
            (0..nn).map(|i| {
                if !zr[i].is_finite() { return (0.0, 0.0); }
                let (x, y) = (i % n, i / n);
                (d1(i, if x > 0 { Some(i - 1) } else { None }, if x + 1 < n { Some(i + 1) } else { None }), d1(i, if y > 0 { Some(i - n) } else { None }, if y + 1 < n { Some(i + n) } else { None }))
            }).collect()
        } else { Vec::new() };
        if tread.is_some() { ct_tread = vec![0u16; nn]; ct_tread_tan = vec![0u16; nn]; }
        let dip_clamp = p.dip_clamp_ppm as f64 / 1e6;
        let gl: Vec<f64> = if local {
            let zr = z.clone();
            let d1 = |i: usize, a: Option<usize>, b: Option<usize>| -> f64 {
                let za = a.map(|j| zr[j]).filter(|v| v.is_finite());
                let zb = b.map(|j| zr[j]).filter(|v| v.is_finite());
                match (za, zb) { (Some(p), Some(q)) => (q - p) / (2.0 * dx), (None, Some(q)) => (q - zr[i]) / dx, (Some(p), None) => (zr[i] - p) / dx, _ => 0.0 }
            };
            (0..nn).map(|i| {
                if !zr[i].is_finite() { return f64::NAN; }
                let (x, y) = (i % n, i / n);
                let gx = d1(i, if x > 0 { Some(i - 1) } else { None }, if x + 1 < n { Some(i + 1) } else { None });
                let gy = d1(i, if y > 0 { Some(i - n) } else { None }, if y + 1 < n { Some(i + n) } else { None });
                gx.hypot(gy)
            }).collect()
        } else { Vec::new() };
        let h_pre: Vec<f64> = if tread.is_some() { z.clone() } else { Vec::new() };
        // Lane PERIOD-ENDS (ct_ends): publish the period ends of each remapped post.
        let ends = c.v2.as_ref().map(|v| v.ends == 1).unwrap_or(false);
        if ends { ct_ends = vec![[f64::NAN; 4]; nn]; ct_knots = vec![[f64::NAN; 4]; nn]; ct_grad = vec![0.0; nn]; ct_pkind = vec![crate::t1b_ct::PK_NONE; nn]; ct_bed = vec![0u32; nn]; }
        // The package (lo, hi heights; var 2 / var 3), the var 3 bed hard share, keyed (a var 3 bed inside a thick package) and
        // the var 3 bed id of post i at height x reading a class of structure `mode`.
        let per_post = lp.rock.per_post();
        let mut pkg_at = |i: usize, x: f64, mode: u32| -> (Option<(f64, f64)>, Option<f64>, bool, u32) {
            if var == 2 && mode == crate::t1b_ct::CT_LAYERED {
                let p = lp.rock.column(i, |cell, datum, base, codes| {
                    let own;
                    let v: &Vec<(f64, f64)> = if per_post { own = crate::t1b_ct::packages(c, base, codes); &own } else { pk.entry(cell).or_insert_with(|| crate::t1b_ct::packages(c, base, codes)) };
                    let d = datum - x;
                    v.iter().find(|q| d >= q.0 && d < q.1).map(|q| (datum - q.1, datum - q.0))
                }).flatten();
                (p, None, false, 0)
            } else if let (3, Some(bp), true) = (var, bedp, mode == crate::t1b_ct::CT_LAYERED) {
                let r = lp.rock.column_ix(i, |cell, datum, base, codes, idx, full| {
                    let own;
                    let v: &Vec<(f64, f64)> = if per_post { own = crate::t1b_ct::packages(c, base, codes); &own } else { pk.entry(cell).or_insert_with(|| crate::t1b_ct::packages(c, base, codes)) };
                    let d = datum - x;
                    let q = *v.iter().find(|q| d >= q.0 && d < q.1)?;
                    if !q.0.is_finite() && !q.1.is_finite() { return Some((f64::NEG_INFINITY, f64::INFINITY, None, false, 0u32)); }
                    let j = base.iter().position(|x| *x == q.1).unwrap_or(base.len() - 1);
                    let id = full.saturating_sub(1 + idx[j] as usize);
                    let tb_ = cx.ct_beds.pk.get(id).map(|t| t.as_slice()).unwrap_or(&[]);
                    let bid = |k: usize| (((id as u32).saturating_add(1)).min(0xFFFF) << 16) | (k.min(0xFFFF) as u32);
                    if q.1 - q.0 <= umax_m {
                        return Some((datum - q.1, datum - q.0, Some(tb_.first().map(|e| e.1).unwrap_or(bp.eta_mid())), false, bid(0)));
                    }
                    let (b0, b1, e, k) = crate::t1b_ct::bed_in_package_k(bp, step_m, tb_, q.0, q.1, d)?;
                    Some((datum - b1, datum - b0, Some(e), true, bid(k)))
                }).flatten();
                match r {
                    Some((lo, hi, e, keyed, bid)) => (Some((lo, hi)), if !lo.is_finite() { crate::t1b_ct::glob_eta(&cx.ct_beds, x) } else { e }, keyed, bid),
                    None => (None, crate::t1b_ct::glob_eta(&cx.ct_beds, x), false, 0),
                }
            } else { (None, None, false, 0) }
        };
        for i in 0..nn {
            if src(i) || !z[i].is_finite() { continue; }
            let h = z[i];
            let (sc, _, code) = lp.rock.at(i, h);
            let s0 = (sc * rmk[i] * fr2[i]).max(smin);
            let s = if local { gl[i] } else { s0 };
            // Lane VOLCANIC-SPLIT volcanic_ct_v2 (vs=1): split posts take the riser ladder and the lens hard share instead.
            // Lane SPLIT-T1B (vs_ct_blend=1): the weighted mean over the riser-threshold ladder alternatives; vs_ct_bed=1: the
            // cliff limb tied to the lens's lava (t1b_volc CT_BLEND_RULE / CT_BED_RULE).
            let vblend = lp.rock.volc_rule().map(|v| v.ct_blend.is_some()).unwrap_or(false);
            let valts = if vblend { lp.rock.volc_ct_alts(i, h) } else { None };
            if let Some(vl) = if vblend { valts.as_ref().map(|_| None) } else { lp.rock.volc_ct(i, h) } {
                let Some(vr) = lp.rock.volc_rule() else { continue };
                let tt = tread.map(|t| t.tread_tan_permille as f64 / 1000.0).unwrap_or(0.14);
                let bl = if geo[i].base.is_finite() { geo[i].base.max(sea) } else { f64::NAN };
                let bed = vr.ct_bed.is_some();
                let (l, r, zb) = if let Some(alts) = valts.as_ref() {
                    let mut zb = 0.0;
                    let mut best: Option<(f64, [f64; 4], crate::t1b_volc::VolcRemap)> = None;
                    for (al, w) in alts {
                        let Some((l, eta)) = al else { zb += w * h; continue };
                        let r = crate::t1b_volc::ct_remap_bed(s, *eta, tt, vr.mass(), vr.face(), h, l[1], l[2], bl, bed);
                        zb += w * r.z;
                        if best.map(|b| *w > b.0).unwrap_or(true) { best = Some((*w, *l, r)); }
                    }
                    let Some((_, l, r)) = best else { continue };
                    if (zb - r.z).abs() > 1e-3 { st.ct_volc_blend_posts += 1; }
                    (l, r, zb)
                } else {
                    let Some((l, eta)) = vl else { continue };
                    let r = crate::t1b_volc::ct_remap_bed(s, eta, tt, vr.mass(), vr.face(), h, l[1], l[2], bl, bed);
                    (l, r, r.z)
                };
                let held = r.clamped || (zb - r.z).abs() > 1e-3;
                if r.clamped { st.ct_tread_base_clamped_posts += 1; }
                if ends { ct_ends[i] = [l[1], l[2], l[0], l[3]]; ct_knots[i] = r.knots; ct_grad[i] = s; ct_pkind[i] = crate::t1b_ct::PK_PACKAGE | if held { crate::t1b_ct::PK_HELD } else { 0 }; st.ct_ends_posts += 1; }
                if !ct_tread.is_empty() && r.ft > 0.0 { ct_tread[i] = (r.ft * 1000.0).round().clamp(0.0, 1000.0) as u16; ct_tread_tan[i] = (tt * 1000.0).round() as u16; st.ct_tread_profile_posts += 1; }
                if r.rg == crate::t1b_ct::RG_TREAD { st.ct_tread_posts += 1; } else if r.rg == crate::t1b_ct::RG_CLIFF { st.ct_cliff_posts += 1; } else { st.ct_talus_posts += 1; }
                ct_regime[i] = r.rg;
                ct_cliff[i] = (r.fc * 1000.0).round().clamp(1.0, 1000.0) as u16;
                let d = zb - h;
                if d != 0.0 { st.ct_shifted_posts += 1; st.ct_shift_max_abs_m = st.ct_shift_max_abs_m.max(d.abs()); sum += d; }
                z[i] = zb;
                if ct_base { ct_shift[i] = d; }
                continue;
            }
            let Some((tb, ta, cf)) = crate::t1b_ct::split(c, code, s, rmk[i]) else {
                if local && crate::t1b_ct::split(c, code, s0, rmk[i]).is_some() { st.ct_local_mantled_posts += 1; }
                continue
            };
            let mode = c.classes[code as usize].mode;
            // Lane BENCH-VAR (var 3): the post's keyed bed (or record package) and its hard share.
            let (pkg, bed_eta, keyed, bed_id) = pkg_at(i, h, mode);
            if keyed { st.ct_bed_posts += 1; }
            if mode == crate::t1b_ct::CT_LAYERED { if pkg.is_some() { st.ct_strata_posts += 1; } else { st.ct_closure_posts += 1; } }
            let (lo, hi) = pkg.unwrap_or((f64::NEG_INFINITY, f64::INFINITY));
            let mut kn = [f64::NAN; 4];
            let (zn, rg, cfp) = if let Some(t) = tread {
                let t0 = t.tread_tan_permille as f64 / 1000.0;
                let a = if gv.is_empty() { 0.0 } else {
                    let (gx, gy) = gv[i];
                    let m = gx.hypot(gy);
                    match lp.rock.bed_grad(i, dip_clamp) { Some((be, bn)) if m > 1e-12 => (be * gx + bn * gy) / m, _ => 0.0 }
                };
                let tt = (t0 + a).clamp(t.tread_min_tan_permille as f64 / 1000.0, t.tread_max_tan_permille as f64 / 1000.0);
                let used_bed = std::cell::Cell::new(false);
                let mut er = |s0: f64, s1: f64| -> Option<f64> {
                    let rec = lp.rock.column(i, |_, datum, base, codes| crate::t1b_ct::hard_share(c, datum, base, codes, s0, s1));
                    // Lane BENCH-VAR: the record share where it has one (mixed classes), else the bed's keyed share.
                    if bed_eta.is_some() && !rec.map(|e| e > 1e-9 && e < 1.0 - 1e-9).unwrap_or(false) { used_bed.set(true); return bed_eta; }
                    rec
                };
                let bl = if geo[i].base.is_finite() { geo[i].base.max(sea) } else { f64::NAN };
                let r = crate::t1b_ct::remap3(c, t, tab, mode, tb, ta, cf, s, geo[i], h, lo, hi, tt, bl, &mut er);
                if r.clamped { st.ct_tread_base_clamped_posts += 1; if ends { ct_pkind[i] |= crate::t1b_ct::PK_HELD; } }
                if r.limbs3 { st.ct_tread_profile_posts += 1; ct_tread[i] = (r.ft * 1000.0).round().clamp(0.0, 1000.0) as u16; ct_tread_tan[i] = (tt * 1000.0).round() as u16; }
                if r.fallback { st.ct_tread_fallback_posts += 1; }
                if r.limbs3 || r.fallback { if used_bed.get() { st.ct_bed_eta_posts += 1; } else if r.record { st.ct_tread_record_posts += 1; } else { st.ct_tread_closure_posts += 1; } }
                if r.rg == crate::t1b_ct::RG_TREAD { st.ct_tread_posts += 1; }
                if ends { kn = if r.limbs3 { crate::t1b_ct::knots3(r.ft, r.fb, tt, tb, s) } else { crate::t1b_ct::knots2(cf, tb, s) }; }
                (r.z, r.rg, r.fc)
            } else {
                let (zn, rg) = crate::t1b_ct::remap(c, tab, mode, tb, ta, cf, s, geo[i], h, lo, hi);
                if ends { kn = crate::t1b_ct::knots2(cf, tb, s); }
                (zn, rg, cf)
            };
            if ends {
                if let Some((s0, s1, k)) = crate::t1b_ct::remap_period(c, tab, mode, cf, geo[i], h, lo, hi) {
                    let k = if k == crate::t1b_ct::PK_PACKAGE && keyed { crate::t1b_ct::PK_BED } else { k };
                    let (mut sdn, mut sup) = (f64::NAN, f64::NAN);
                    let (mut gap, mut gap_dn) = (false, false);
                    if k != crate::t1b_ct::PK_MASSIVE {
                        for (up, x) in [(false, s0 - 1e-6), (true, s1 + 1e-6)] {
                            let (_, _, code2) = lp.rock.at(i, x);
                            let m2 = c.classes.get(code2 as usize).map(|q| q.mode).unwrap_or(crate::t1b_ct::CT_NONE);
                            if m2 != crate::t1b_ct::CT_LAYERED { continue; }
                            let (p2, _, _, _) = pkg_at(i, x, m2);
                            let (l2, h2) = p2.unwrap_or((f64::NEG_INFINITY, f64::INFINITY));
                            let Some((a, b, _)) = crate::t1b_ct::remap_period(c, tab, m2, 0.0, geo[i], x, l2, h2) else { continue };
                            if up { sup = b; gap |= (a - s1).abs() > 1e-3; } else { sdn = a; gap_dn = (b - s0).abs() > 1e-3; gap |= gap_dn; }
                        }
                    }
                    ct_ends[i] = [s0, s1, sdn, sup];
                    ct_knots[i] = kn;
                    ct_grad[i] = s;
                    ct_pkind[i] = (ct_pkind[i] & crate::t1b_ct::PK_HELD) | k;
                    ct_bed[i] = bed_id;
                    st.ct_ends_posts += 1;
                    if gap { st.ct_ends_gap_posts += 1; if gap_dn && k == crate::t1b_ct::PK_APRON { st.ct_ends_gap_apron_dn_posts += 1; } }
                }
            }
            let d = zn - h;
            ct_regime[i] = rg;
            ct_cliff[i] = (cfp * 1000.0).round().clamp(1.0, 1000.0) as u16;
            let g = geo[i];
            if mode == crate::t1b_ct::CT_MASSIVE && cf < 1.0 && g.base.is_finite() && g.crest.is_finite() && g.crest > g.base {
                ct_top[i] = vz_tilekey::chart::quantize_mm((g.base + (1.0 - cf) * (g.crest - g.base)) * 1000.0) as i32;
            }
            if rg == crate::t1b_ct::RG_CLIFF { st.ct_cliff_posts += 1; } else if rg == crate::t1b_ct::RG_TALUS { st.ct_talus_posts += 1; }
            if d != 0.0 { st.ct_shifted_posts += 1; st.ct_shift_max_abs_m = st.ct_shift_max_abs_m.max(d.abs()); sum += d; }
            z[i] = zn;
            if ct_base { ct_shift[i] = d; }
        }
        // Lane TREAD2 (t1b_ct::TREAD_DRAIN_RULE): the remap creates no closed depression: in ascending pre-remap height every
        // remapped post stands at least as high as the lowest of its D8 neighbours that were lower before the remap.
        if tread.is_some() {
            let mut idx: Vec<usize> = (0..nn).filter(|&i| ct_regime[i] != crate::t1b_ct::RG_SOIL && h_pre[i].is_finite()).collect();
            idx.sort_by(|&a, &b| h_pre[a].total_cmp(&h_pre[b]).then(a.cmp(&b)));
            for &i in &idx {
                let (x, y) = ((i % n) as i64, (i / n) as i64);
                let mut lo = f64::INFINITY;
                for (ox, oy) in D8 {
                    let (xx, yy) = (x + ox, y + oy);
                    if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                    let j = (yy as usize) * n + xx as usize;
                    if h_pre[j].is_finite() && z[j].is_finite() && h_pre[j] < h_pre[i] && z[j] < lo { lo = z[j]; }
                }
                if lo.is_finite() && z[i] < lo {
                    let d0 = z[i] - h_pre[i];
                    z[i] = lo;
                    let d = lo - h_pre[i];
                    sum += d - d0;
                    st.ct_tread_drain_raised_posts += 1;
                    if !ct_pkind.is_empty() { ct_pkind[i] |= crate::t1b_ct::PK_HELD; }
                    if ct_base { ct_shift[i] = d; }
                }
            }
        }
        if st.ct_shifted_posts > 0 { st.ct_shift_mean_m = sum / st.ct_shifted_posts as f64; }
    }
    st.seconds[3] = t.elapsed().as_secs_f64();
    // Lane GLACIAL-T1 (glacial key): U sections and cirque bowls from the published chart table (pointwise reads).
    let mut gl_low: Vec<f32> = Vec::new();
    if let (Some(gt), Some(gp)) = (cx.gl.as_ref(), p.glacial.as_ref()) {
        let (gs, low) = crate::t1b_glacial::carve(&lat, &mut z, &fixed, &lp.water, net, curves, gt, gp);
        st.gl = Some(gs);
        gl_low = low;
    }
    // Lane PERIOD-ENDS: the heights at the end of the remap (later raises flag the post: height != F(h)).
    let z_rm: Vec<f64> = if ct_ends.is_empty() { Vec::new() } else { z.clone() };
    if !diag.lab.is_empty() { diag.z_crest = z.clone(); }
    // Lane MS-CEIL ceiling_min_v1 (ceil only, CEIL_RULE): free non-floor land posts take min(z_ss, Z_c); after the ct remap
    // (z_rm above, so ct posts the ceiling lowers are flagged moved) and before the footslope, fans, bank rule and fill.
    let (ceil_m, mut ceil_on): (Vec<f32>, Vec<bool>) = match p.ceil.as_ref() { Some(cp) => (ceil_surface(lp, t0, cp), vec![false; nn]), None => (Vec::new(), Vec::new()) };
    // Lane MS-CEIL-2 (ceil_ped=1, CEIL_PED_RULE): the height before the min of every post the ceiling lowers.
    let mut ceil_pre: Vec<f64> = if p.ceil.as_ref().map(|c| c.ped == 1).unwrap_or(false) { vec![f64::NAN; nn] } else { Vec::new() };
    if !ceil_m.is_empty() {
        for i in 0..nn {
            if fixed[i] != FREE || floor[i] || !z[i].is_finite() { continue; }
            let c = ceil_m[i] as f64;
            if z[i] > c { st.ceil_posts += 1; st.ceil_lowered_max_m = st.ceil_lowered_max_m.max(z[i] - c); if !ceil_pre.is_empty() { ceil_pre[i] = z[i]; } z[i] = c; ceil_on[i] = true; }
        }
    }
    // Lane VALLEY-FLOOR-FIX colluvial_footslope_v2: the concave wall foot of wide floors filled on the rounded surface.
    if vf_foot {
        let lf = vfp.unwrap().foot_length_mm as f64 / 1000.0;
        let free: Vec<bool> = (0..nn).map(|i| fixed[i] == FREE).collect();
        let wide: Vec<bool> = (0..nn).map(|i| floor[i] && floor_h[i] >= lf).collect();
        crate::t1b_vf::foot_fill(n, dx, &free, &floor, &wide, &mut z, &|i, zz| lp.rock.at(i, zz).0, p.sc_ref_permille as f64 / 1000.0, t0.sea_level_m + 1.0, lf, &mut vfs);
    }
    // Lane VALLEY-FLOOR-FIX tributary_fan_v1: deposit on the rounded surface, before the bank rule and the fill.
    let mut fan_mask: Vec<bool> = Vec::new();
    if vf_fan {
        let v = vfp.unwrap();
        fan_mask = vec![false; nn];
        let idx = |x: f64, y: f64| lat.idx((x / dx).round() as i64, (y / dx).round() as i64);
        let ap = crate::t1b_vf::fan_apexes(curves, &idx, &tr_area, v.fan_area_ratio_permille as f64 / 1000.0);
        for (i, q) in ap {
            if !tr_val[i].is_finite() || fixed[i] == MARINE || fixed[i] == LAKE { continue; }
            let step = q[2] - tr_val[i];
            if step < v.fan_height_min_mm as f64 / 1000.0 { continue; }
            let a = q[3];
            let sf = v.fan_slope(a);
            let supply = e[i] as f64 * 1e-6 * a * v.fan_time_yr as f64;
            let r_mass = (6.0 * supply / (std::f64::consts::PI * sf)).cbrt();
            let r_geom = step / sf;
            let r_trim = tr_d[i];
            let rad = r_mass.min(r_geom).min(r_trim);
            if !(rad >= dx) { continue; }
            let za = tr_val[i] + sf * rad;
            let rr = (rad / dx).ceil() as i64;
            let (cx0, cy0) = ((q[0] / dx).round() as i64, (q[1] / dx).round() as i64);
            let mut vol = 0.0f64;
            for oy in -rr..=rr { for ox in -rr..=rr {
                let (gx, gy) = (cx0 + ox, cy0 + oy);
                let Some(k) = lat.idx(gx, gy) else { continue };
                if fixed[k] == CHANNEL || fixed[k] == MARINE || fixed[k] == LAKE || !z[k].is_finite() { continue; }
                let r = (gx as f64 * dx - q[0]).hypot(gy as f64 * dx - q[1]);
                if r > rad { continue; }
                let zc = za - sf * r;
                if zc > z[k] { vol += (zc - z[k]) * dx * dx; if zc - z[k] >= 0.1 { fan_mask[k] = true; } z[k] = zc; vfs.fan_raised_posts += 1; }
            } }
            vfs.fan_volume_m3 += vol;
            vfs.fans.push(crate::t1b_vf::FanRec { x_m: q[0], y_m: q[1], area_m2: a, step_m: step, slope: sf, supply_m3: supply, r_mass_m: r_mass, r_geom_m: r_geom, r_trim_m: r_trim, radius_m: rad, apex_m: za });
        }
    }
    if vfp.is_some() { st.vf = Some(std::mem::take(&mut vfs)); }
    // Bank rule (T1a channels): a free post next to a T1a channel post stands at least 1 cm above the highest adjacent bed.
    let ls_ramp = p.ls.as_ref().map(|l| l.bank_ramp == 1).unwrap_or(false);
    let mut bank_raised = vec![false; if ls_ramp { nn } else { 0 }];
    for i in 0..nn {
        if fixed[i] != FREE { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let mut hb = f64::NEG_INFINITY;
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if fixed[j] == CHANNEL { hb = hb.max(z[j]); }
        }
        if hb.is_finite() && z[i] < hb + 0.01 { st.bank_raised_posts += 1; st.bank_raised_max_m = st.bank_raised_max_m.max(hb + 0.01 - z[i]); z[i] = hb + 0.01; if ls_ramp { bank_raised[i] = true; } }
    }
    // Lane TA LS2: the bank raise spread as a ramp at the hillslope gradient (no one-post walls).
    if ls_ramp {
        let zpre = z.clone();
        let (c, m) = crate::t1b_ls::bank_ramp(n, dx, &mut z, &|i| fixed[i] == FREE, &bank_raised, &|i| speed2(i, zpre[i]));
        st.bank_ramp_posts = c;
        st.bank_ramp_max_m = m;
    }
    let t = std::time::Instant::now();
    let outlet: Vec<bool> = (0..nn).map(|i| fixed[i] != FREE).collect();
    // Lane SEAM-CERT: with seam_cert every construction fill is seeded only from fixed posts and records its lake taint.
    let mut seam: Option<crate::t1b_seam::SeamTaint> = p.seam_cert.as_ref().map(|_| crate::t1b_seam::SeamTaint::new(n, lp.publish_margin_posts.max(0) as usize));
    if stop == Some(crate::t1b_surface_stages::Stage::Main) { return Ok(crate::t1b_surface_stages::Evaluated::Captured({ let mut c=crate::t1b_surface_stages::capture(&z, &fixed, seam.as_ref());c.links=chan_link.clone();c })); }
    let (r, m) = if let Some(owner) = owner { crate::t1b_surface_stages::apply(owner, crate::t1b_surface_stages::Stage::Main, &lat, &mut z, &mut fixed, seam.as_mut())? } else if let Some(sm) = seam.as_mut() { crate::t1b_seam::pinned_fill_native(&lat, &mut z, &fixed, sm, "main") } else if crate::t1b_tf::tf(p).map(|t| t.report == 1).unwrap_or(false) {
        let (r, m, rep) = crate::t1b_tf::fill_report(&lat, &mut z, &outlet, sea, &|i| fixed[i] == FREE, lp.publish_margin_posts, 8000.0);
        st.tf_fill = Some(rep);
        (r, m)
    } else { fill(&lat, &mut z, &outlet) };
    st.fill_raised_posts = r;
    st.fill_raised_max_m = m;
    st.seconds[4] = t.elapsed().as_secs_f64();
    if !diag.lab.is_empty() { diag.z_fill = z.clone(); diag.fixed_pre = fixed.clone(); }
    // Lane T1B-PEDESTAL: measurement dump (env var VZ_PEDESTAL_DUMP, read-only, no bytes; t1b_ped.rs).
    let ped_dump = crate::t1b_ped::dump_dir();
    let mut ped_zpre: Vec<f64> = if ped_dump.is_some() || crate::t1b_ped::needs_ground(p) { z.clone() } else { Vec::new() };
    // Lane MS-CEIL-2 (ceil_ped=1, CEIL_PED_RULE): the pedestal ground on ceiled posts is the larger of now and before the min.
    if !ceil_pre.is_empty() && !ped_zpre.is_empty() { for i in 0..nn { if ceil_pre[i] > ped_zpre[i] { ped_zpre[i] = ceil_pre[i]; } } }
    // Lane CO: volcano edifices, flooding at present sea level, fans, wave cut, final fill.
    let mut cst = crate::t1b_coast::CoastStats::default();
    cst.seam = seam.take();
    let mut vid = vec![255u8; nn];
    let mut shore: Vec<f64> = Vec::new();
    let mut uncut: Vec<f64> = Vec::new();
    let (mut ve_style, mut ve_unit): (Vec<u8>, Vec<u16>) = (Vec::new(), Vec::new());
    let mut ve_front: Vec<f64> = Vec::new();
    let over = if p.coast.is_some() || p.volcano.is_some() {
        let rock_k = |i: usize, zz: f64| lp.rock.at(i, zz).1;
        let rock_sc = |i: usize, zz: f64| lp.rock.at(i, zz).0;
        let tfin = crate::t1b_tf::TfIn { bil_c1: &lp.bil_c1, margin_posts: lp.publish_margin_posts };
        let tg_shore = p.coast.is_some() && crate::t1b_tg::tg_on(p, |t| t.shore);
        if tg_shore { cst.tg_cut_ext = vec![f64::NAN; nn]; }
        let th_uncut = p.coast.is_some() && crate::t1b_th::th_on(p, |t| t.uncut);
        if th_uncut { cst.th_uncut = vec![f64::NAN; nn]; }
        // Lane VE: per-post style / unit arrays (the coast stage fills them).
        if crate::t1b_ve::ve(p).and_then(|v| v.ve4.as_ref()).map(|v| v.front == 1).unwrap_or(false) && p.volcano.is_some() && cx.ve_units { cst.ve_front = vec![f64::NAN; nn]; }
        // Lane DAMMED-VALLEY: heights before the coast / volcano final fill (the stage records them when allocated).
        if p.dam_lake.is_some() && p.volcano.is_some() { cst.dam_pre = vec![f64::NAN; nn]; }
        // Lane EDIFICE-FOOT: measurement dump (env var VZ_EDIFICE_DUMP, read-only, no bytes; t1b_ef.rs).
        let ef_dump = std::env::var_os("VZ_EDIFICE_DUMP").filter(|_| p.volcano.is_some());
        let ef_zpre: Vec<f64> = if ef_dump.is_some() { if cst.dam_pre.len() != nn { cst.dam_pre = vec![f64::NAN; nn]; } z.clone() } else { Vec::new() };
        if crate::t1b_ve::ve(p).is_some() && p.volcano.is_some() && cx.ve_units { cst.ve_style = vec![crate::t1b_ve::ST_NONE; nn]; cst.ve_unit = vec![crate::t1b_ve::UNIT_NONE; nn]; }
        // implementation review material-stage observer hook (integrated by lane TG, 2026-10-15): the binary passes an observer only when the
        // environment variable VZ_COAST_OBSERVE is set (parity measurement: read-only callbacks, stdout summary, no bytes).
        let tfin_o = if cx.tf.is_some() { Some(&tfin) } else { None };
        let (o, v) = if std::env::var_os("VZ_COAST_OBSERVE").is_some() {
            let (mut stages, mut h) = (Vec::new(), 0xcbf29ce484222325u64);
            let mut obs = |stage: crate::t1b_coast::CoastMaterialStage, zz: &[f64], ff: &[u8], oo: &[u8]| {
                stages.push(stage);
                for (k, v) in zz.iter().enumerate() { h = (h ^ v.to_bits() ^ ((ff[k] as u64) << 56) ^ ((oo[k] as u64) << 48)).wrapping_mul(0x100000001b3); }
            };
            let r = crate::t1b_coast::coast_and_volcano_tf_owned(&lat, sea, p.coast.as_ref(), p.volcano.as_ref(), cx, &rock_k, if p.fx.is_some() { Some(&rock_sc) } else { None }, &t1a_post, bil_w, &mut z, &mut fixed, &mut cst, tfin_o, Some(&mut obs), owner, stop)?;
            println!("  coast observer: stages {:?} snapshot hash {:016x}", stages, h);
            r
        } else {
            crate::t1b_coast::coast_and_volcano_tf_owned(&lat, sea, p.coast.as_ref(), p.volcano.as_ref(), cx, &rock_k, if p.fx.is_some() { Some(&rock_sc) } else { None }, &t1a_post, bil_w, &mut z, &mut fixed, &mut cst, tfin_o, None, owner, stop)?
        };
        if let Some(mut capture)=cst.surface_capture.take() { capture.links=chan_link.clone();return Ok(crate::t1b_surface_stages::Evaluated::Captured(capture)); }
        // Lane TG G1: the shore field on the surface before marine posts are published at sea level (cut posts: the cut
        // profile extended seaward).
        if tg_shore {
            let zc: Vec<f64> = (0..nn).map(|i| if cst.tg_cut_ext[i].is_finite() && fixed[i] == MARINE { cst.tg_cut_ext[i] } else { z[i] }).collect();
            let cap = crate::t1b_tg::tg(p).unwrap().shore_cap_mm as f64 / 1000.0;
            shore = crate::t1b_tg::shore_field(&lat, &zc, sea, &|i| fixed[i] == MARINE, cap);
            cst.tg_cut_ext = Vec::new();
        }
        if th_uncut { uncut = std::mem::take(&mut cst.th_uncut); }
        if let Some(d) = ef_dump.as_ref() { crate::t1b_ef::dump(d, &lat, &ef_zpre, &cst.dam_pre, &z, &o, &cst.ve_style, &fixed, &v); }
        ve_style = std::mem::take(&mut cst.ve_style);
        ve_unit = std::mem::take(&mut cst.ve_unit);
        ve_front = std::mem::take(&mut cst.ve_front);
        if p.coast.is_some() { for i in 0..nn { if fixed[i] == MARINE { z[i] = sea; } } }
        vid = v;
        o
    } else { vec![0u8; nn] };
    // Lane VALLEY-FLOOR-FIX: tributary fan posts expose the fan deposit.
    let mut over = over;
    seam = cst.seam.take();
    // Lane VALLEY-FLOOR-FIX delta_surface_v1 (vf_floor, coast on): a valley floor buried by lane CO's bayhead fan / fan delta
    // takes the same natural levee beside its (now estuary) channels and the same colluvial footslope at its walls, on the
    // fan surface; then the final fill again (outlets = fixed posts).
    let mut owner_vf_applied=false;
    if vf_floor && p.coast.is_some() {
        let v = vfp.unwrap();
        let co_fan: Vec<bool> = (0..nn).map(|i| over[i] == crate::t1b_coast::OVER_FAN && fixed[i] == FREE && z[i].is_finite()).collect();
        if co_fan.iter().any(|&b| b) {
            let mut best: Vec<(f64, i64, i64, f64)> = vec![(f64::INFINITY, 0, 0, 0.0); nn];
            for i in 0..nn {
                if chan_area[i] <= 0.0 { continue; }
                let a = chan_area[i] as f64;
                let (w, _) = v.hydraulic(a);
                let rmax = (4.0 * v.levee_decay_widths_permille as f64 / 1000.0 * w + 0.5 * w).min(600.0);
                if rmax < dx { continue; }
                let r = (rmax / dx).floor() as i64;
                let (x, y) = ((i % n) as i64, (i / n) as i64);
                let (gx, gy) = lat.g(i);
                for oy in -r..=r { for ox in -r..=r {
                    let d = ((ox * ox + oy * oy) as f64).sqrt() * dx;
                    if d > rmax || (ox == 0 && oy == 0) { continue; }
                    let (xx, yy) = (x + ox, y + oy);
                    if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
                    let j = yy as usize * n + xx as usize;
                    if !co_fan[j] { continue; }
                    let cur = best[j];
                    if d < cur.0 || (d == cur.0 && (gy, gx) < (cur.1, cur.2)) { best[j] = (d, gy, gx, a); }
                } }
            }
            let mut dsum = 0.0f64;
            for j in 0..nn { if best[j].0.is_finite() { let l = v.levee(best[j].3, best[j].0); if l > 0.0 { z[j] += l; dsum += l * dx * dx; if let Some(s) = st.vf.as_mut() { s.delta_levee_posts += 1; } } } }
            if let Some(s) = st.vf.as_mut() { s.delta_levee_volume_m3 += dsum; }
            if v.foot_on() {
                let lf = v.foot_length_mm as f64 / 1000.0;
                let free: Vec<bool> = (0..nn).map(|i| fixed[i] == FREE).collect();
                let mut fs2 = crate::t1b_vf::VfStats::default();
                crate::t1b_vf::foot_fill(n, dx, &free, &co_fan, &co_fan, &mut z, &|i, zz| lp.rock.at(i, zz).0, p.sc_ref_permille as f64 / 1000.0, t0.sea_level_m + 1.0, lf, &mut fs2);
                if let Some(s) = st.vf.as_mut() { s.delta_foot_raised_posts = fs2.foot_raised_posts; s.delta_foot_volume_m3 = fs2.foot_volume_m3; }
            }
            let outlet: Vec<bool> = (0..nn).map(|i| fixed[i] != FREE).collect();
            if stop == Some(crate::t1b_surface_stages::Stage::VfDelta) { return Ok(crate::t1b_surface_stages::Evaluated::Captured({ let mut c=crate::t1b_surface_stages::capture(&z, &fixed, seam.as_ref());c.links=chan_link.clone();c })); }
            let (r, m) = if let Some(owner) = owner { { owner_vf_applied=true; crate::t1b_surface_stages::apply(owner, crate::t1b_surface_stages::Stage::VfDelta, &lat, &mut z, &mut fixed, seam.as_mut())? } } else if let Some(sm) = seam.as_mut() { crate::t1b_seam::pinned_fill_native(&lat, &mut z, &fixed, sm, "vf_delta") } else { fill(&lat, &mut z, &outlet) };
            if let Some(s) = st.vf.as_mut() { s.delta_fill_raised_posts = r; s.delta_fill_raised_max_m = m; }
        }
    }
    if stop.is_none() && !owner_vf_applied { if let Some(owner)=owner { crate::t1b_surface_stages::apply(owner,crate::t1b_surface_stages::Stage::VfDelta,&lat,&mut z,&mut fixed,seam.as_mut())?; } }
    for k in 0..fan_mask.len() { if fan_mask[k] && over[k] == 0 { over[k] = crate::t1b_coast::OVER_FAN; } }
    if !z_rm.is_empty() {
        for i in 0..nn {
            if ct_pkind[i] & 0x7F == crate::t1b_ct::PK_NONE || ct_pkind[i] & crate::t1b_ct::PK_HELD != 0 { continue; }
            if z[i] != z_rm[i] { ct_pkind[i] |= crate::t1b_ct::PK_HELD; st.ct_ends_moved_posts += 1; }
        }
    }
    if let Some(d) = ped_dump.as_ref() { crate::t1b_ped::dump(d, &lat, &ped_zpre, &z, &over, &fixed, &lp.bil, &lp.cell, &e); }
    let ground_pre = if crate::t1b_ped::needs_ground(p) { ped_zpre } else { Vec::new() };
    // Lane MS-CEIL: surface_src on the final surface (CEIL_SRC_NOTE).
    let ceil_src: Vec<u8> = if ceil_on.is_empty() { Vec::new() } else {
        (0..nn).map(|i| if fixed[i] != FREE || floor[i] || over[i] == crate::t1b_coast::OVER_EDIFICE || !z[i].is_finite() { 2 } else if ceil_on[i] { 1 } else { 0 }).collect()
    };
    if stop.is_some() { return Ok(crate::t1b_surface_stages::Evaluated::Captured({ let mut c=crate::t1b_surface_stages::capture(&z, &fixed, seam.as_ref());c.links=chan_link.clone();c })); }
    Ok(crate::t1b_surface_stages::Evaluated::Surface(Surface { dependency, ground_pre, z, z_eik, fixed, floor, chan_area, chan_link, e, stats: st, over, coast: cst, diag, vid, shore, uncut, ve_style, ve_unit, ve_front, ct_regime, ct_cliff, ct_top, ct_shift, ct_tread, ct_tread_tan, ct_ends, ct_knots, ct_grad, ct_pkind, ct_bed, seam, gl_low, ceil_m, ceil_src }))
}

// ------------------------------------------------------------------ lane CO inputs

/// Inputs the lane CO steps need beyond 0.9.0: the volcano records and the keyed scope.
pub struct CoInputs<'a> { pub volcanoes: Option<&'a crate::t1b_coast::VolcanoRecords>, pub seed: u64, pub recipe_digest: String, pub realization: u32,
    /// Lane GLACIAL-T1: the t0.glacial layer (layer build only; windows read the published table).
    pub glacial: Option<&'a crate::t1b_glacial::T0Glacial>,
    /// Lane T1B-PEDESTAL edifice_local_base_v1: the published local ground per volcano id (windows; the layer solve
    /// computes its own). None = bilinear bases.
    pub ped_ground: Option<BTreeMap<String, f64>> }

impl<'a> CoInputs<'a> {
    pub fn none() -> CoInputs<'static> { CoInputs { volcanoes: None, seed: 0, recipe_digest: String::new(), realization: 0, glacial: None, ped_ground: None } }
    /// The coast context for these tables (empty when the lane CO params are absent).
    pub fn ctx(&self, p: &RiverNetParamsV1, t0: &T0, chart: &Chart, cv: &[Curve], table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>) -> R<crate::t1b_coast::CoastCtx> {
        self.ctx_ped(p, t0, chart, cv, table, flood, self.ped_ground.as_ref())
    }
    /// As `ctx` with an explicit local-ground map (lane T1B-PEDESTAL; ignored without ped_base).
    #[allow(clippy::too_many_arguments)]
    pub fn ctx_ped(&self, p: &RiverNetParamsV1, t0: &T0, chart: &Chart, cv: &[Curve], table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, ground: Option<&BTreeMap<String, f64>>) -> R<crate::t1b_coast::CoastCtx> {
        let ground = if crate::t1b_ped::local_base(p) { ground } else { None };
        let mut cx = crate::t1b_coast::CoastCtx { flood: flood.clone(), cf: p.cf.clone(), fx: p.fx.clone(), tf: crate::t1b_tf::tf(p).cloned(), tg: crate::t1b_tg::tg(p).cloned(), tg_floor_slope: p.floor_slope_ppm as f64 / 1e6, th: crate::t1b_th::th(p).cloned(), ve: crate::t1b_ve::ve(p).cloned(), ef: p.edifice_foot.clone(), ..Default::default() };
        // Lane STEP2: the keyed step-thickness closure (global boundaries: the same in every window).
        if let Some(c) = crate::t1b_ct::ct(p).filter(|c| c.v2.as_ref().map(|v| v.var >= 1).unwrap_or(false)) {
            let sc = Scope::root(self.seed).layer(&self.recipe_digest).map_err(|e| e.to_string())?.chart(0).map_err(|e| e.to_string())?;
            if c.v2.as_ref().map(|v| v.var == 3).unwrap_or(false) {
                // Lane BENCH-VAR: per-bed keyed thickness and hard share (purpose declared in params.ls.tb.tc.ct.v2.bed only).
                cx.ct_beds = crate::t1b_ct::bed_tables(c, &mut |k| sc.draw(crate::t1b_ct::PURPOSE_CT_BED, k, self.realization).map(|d| d.uniform_f64()).map_err(|e| e.to_string()))?;
                cx.ct_steps = cx.ct_beds.glob.clone();
            } else {
            cx.ct_steps = crate::t1b_ct::closure_table(c, &mut |k| sc.draw(crate::t1b_ct::PURPOSE_CT_STEP, k, self.realization).map(|d| d.uniform_f64()).map_err(|e| e.to_string()))?;
            }
        }
        if let Some(cp) = &p.coast { cx.fans = crate::t1b_coast::fans_for(cp, cv, t0, chart, table, flood, p.cf.as_ref())?; }
        if let Some(vp) = &p.volcano {
            let recs = self.volcanoes.ok_or("volcano params need the volcano records (volcanoes=)")?;
            // Lane TC fix B1: with ls_curve the vents are sited on the unmeandered T1a links (as without the option).
            let sited: Vec<Curve>;
            let cvs: &[Curve] = if cv.iter().any(|c| c.site.is_some()) {
                sited = cv.iter().map(|c| Curve { pts: c.site.clone().unwrap_or_else(|| c.pts.clone()), code: c.code, link: c.link, label: c.label, site: None, carrier: None, mq: None }).collect();
                &sited
            } else { cv };
            let vep = crate::t1b_ve::ve(p);
            cx.volcanoes = crate::t1b_coast::volcano_sites_ped2(vp, recs, t0, chart, cvs, table, self.seed, &self.recipe_digest, self.realization, None, p.cf.as_ref(), vep, ground, crate::t1b_ped::local_relief(p))?;
            // Lane HS-PASS hs_volcano (window only): empty the radial gully lists after their draws, before the VE units copy them.
            if let Some(h) = p.hs.as_ref().filter(|h| h.volcano >= 1) { crate::t1b_coast::hs_strip_gullies(&mut cx.volcanoes, h.volcano); }
            // Lane VE: eruption style and keyed history per volcano (paths on the drawn curves).
            if let (Some(vep), Some(cfp)) = (vep, p.cf.as_ref()) {
                let onset: BTreeMap<String, f64> = recs.features.iter().filter_map(|f| f.onset_my.map(|o| (f.id.clone(), o))).collect();
                let bi = crate::t1b_ve::BuildIn { vp, cf: cfp, ve: vep, t0, chart, curves: cv, table, seed: self.seed, recipe_digest: &self.recipe_digest, realization: self.realization, onset: &onset, t_my: recs.t_my };
                crate::t1b_ve::build_all(&bi, &mut cx.volcanoes)?;
            }
            cx.t1a_segs = crate::t1b_coast::t1a_segments_near(vp, &cx.volcanoes, cv);
            // Lane TF T5: the C1 T0 surface at each vent (chart m -> lat/lon as the bilinear base).
            if let Some(tfp) = crate::t1b_tf::tf(p).filter(|t| t.cone_c1 == 1) {
                let elev = |c: u32| t0.elevation_m[c as usize] as f64;
                // The bilinear at a global 32 m post exactly as lat_prep computes it (f32 as stored there).
                let bil_g = |gx: i64, gy: i64| -> R<f64> {
                    let (la, lo) = chart.inverse_mm(gx * 32_000, gy * 32_000).map_err(|e| e.to_string())?;
                    Ok(bilinear(t0, &elev, la * vz_tilekey::chart::RAD2DEG, lo * vz_tilekey::chart::RAD2DEG, false) as f32 as f64)
                };
                let half = tfp.cone_smooth_mm as f64 / 1000.0;
                cx.tf_vent_c1 = cx.volcanoes.iter().map(|v| crate::t1b_tf::smooth_bil_at(&bil_g, v.vent_xy_m[0], v.vent_xy_m[1], half)).collect::<R<Vec<f64>>>()?;
            }
        }
        Ok(cx)
    }
}

// ------------------------------------------------------------------ the E solve (chart-wide)

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CellRow {
    pub cell: u32,
    pub solved: bool,
    pub land_hat_fraction: f64,
    pub e0_m_per_my: f64,
    pub e_m_per_my: f64,
    pub hat_mean_z_m: f64,
    pub hat_mean_bilinear_m: f64,
    pub hat_mean_trunk_m: f64,
    pub diff_m: f64,
    pub at_bound: bool,
    /// Lane CO (only with coast / volcano / f6_reference 2): per WHOLE T0 cell, plain over the cell's chart posts.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coastal: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flood_offset_m: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub land_fraction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub land_fraction_bilinear: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_mean_z_m: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_mean_ref_m: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_mean_bilinear_m: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_land_posts: Option<u64>,
    /// Lane TB S3 (R4 measurement, only on cells it applies to): the plain-mean reference is used; the hat residual is kept.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r4_plain_reference: Option<bool>,
    /// Lane TC C3 (R4 variant measurement, only on cells it applies to): not solved, E = the solved 8-neighbours' mean.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r4n_neighbour_mean: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_mean_trunk_m: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hat_diff_m: Option<f64>,
    /// Lane E-SOLVE-GENERALITY es target 1 (only on cells it applies to): the plain whole-cell means are the solve target.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub es_plain_target: Option<bool>,
}

pub struct SolveResult { pub table: BTreeMap<u32, f64>, pub rows: Vec<CellRow>, pub history_max_abs_diff_m: Vec<f64>, pub converged: bool, pub evals: u32, pub seconds: f64,
    /// Lane CO: flooding offsets (empty without coast), the final fans and volcano sites, F4 history (max |land - bilinear| over coastal cells).
    pub flood: BTreeMap<u32, f64>, pub cx: crate::t1b_coast::CoastCtx, pub history_max_abs_f4: Vec<f64>,
    /// Lane TD2: the last evaluation's captured network (warm start of the published one).
    pub net: Option<Network>,
    /// Lane TA (only with ls): per evaluation [passes, receivers changed by the last pass, receivers changed vs the previous evaluation].
    pub ls_history: Vec<[u32; 3]>,
    /// Lane TC C4: the last evaluations ran on a frozen network (the published network is that network).
    pub frozen: bool,
    /// Lane TE (solve_cell only): per-cell solve record for the layer header.
    pub te: Option<Value>,
    /// Lane E-SOLVE-GENERALITY (es only): solve record for the layer header.
    pub es: Option<Value>,
    /// Lane T1B-PEDESTAL (ped only): record for the layer header (lane_ped).
    pub ped: Option<Value> }

#[allow(clippy::too_many_arguments)]
pub fn solve_e(t0: &T0, extras: &T0Extras, skel: &SkeletonView, chart: &Chart, net: &Network, nctx: &NodeCtx, p: &RiverNetParamsV1, strata: Option<&StrataCols>, anchor: (f64, f64), co: &CoInputs, log: &mut dyn FnMut(String)) -> R<SolveResult> {
    // Lane TE: the cell-local solve (opt-in; absent = this chart-wide solve, unchanged).
    if crate::t1b_ls::te_cell(p) { return crate::t1b_cell::solve_e_cells(t0, extras, skel, chart, net, nctx, p, strata, anchor, co, log); }
    let t_all = std::time::Instant::now();
    let dxmm = p.solve_spacing_mm;
    let half = vz_tilekey::CHART_HALF_EXTENT_MM / dxmm;
    let lat = EvalLat { gx0: -half, gy0: -half, n: (2 * half + 1) as usize, dx_mm: dxmm };
    let lp = lat_prep(lat, chart, t0, extras, p, strata, anchor)?;
    let nn = lat.nodes();
    let sea = t0.sea_level_m;
    let f = dxmm / crate::t1b::T1B_SPACING_MM;
    let trunk: Vec<f64> = (0..nn).map(|i| { let (gx, gy) = lat.g(i); crate::t1b::trunk_bilinear(skel, gx * f, gy * f, sea).max(sea) }).collect();
    // Hat masses.
    let mut mass: BTreeMap<u32, (f64, f64)> = BTreeMap::new();
    for i in 0..nn { for k in 0..4 {
        let w = lp.hat_w[i][k] as f64;
        if w <= 0.0 { continue; }
        let e = mass.entry(lp.hat_c[i][k]).or_default();
        e.1 += w;
        if lp.water[i] == FREE { e.0 += w; }
    } }
    let g0 = p.e0_gain_permille as f64 / 1000.0;
    let (emin, emax) = (p.e_min_m_per_my as f64, p.e_max_m_per_my as f64);
    let mut table: BTreeMap<u32, f64> = BTreeMap::new();
    let mut solved: BTreeMap<u32, bool> = BTreeMap::new();
    for (&c, &(lm, tm)) in &mass {
        let lake = t0.water(c) == CellWater::Lake;
        let s = !lake && lm >= p.solve_min_land_hat_permille as f64 / 1000.0 * tm;
        solved.insert(c, s);
        let e0 = if t0.water(c) == CellWater::Land { closure_rate_m_per_my(t0, &extras.laws, c) } else {
            // A marine T0 cell with zero-contour land: start from its land neighbours' mean rate.
            let (i, j) = t0.ij(c);
            let mut acc = (0.0, 0);
            for di in -1i64..=1 { for dj in -1i64..=1 {
                let ii = i as i64 + di; if ii < 0 || ii >= t0.nlat as i64 { continue; }
                let jj = (j as i64 + dj).rem_euclid(t0.nlon as i64);
                let cc = (ii as usize * t0.nlon + jj as usize) as u32;
                if t0.water(cc) == CellWater::Land { acc.0 += closure_rate_m_per_my(t0, &extras.laws, cc); acc.1 += 1; }
            } }
            if acc.1 > 0 { acc.0 / acc.1 as f64 } else { 0.0 }
        };
        table.insert(c, if s || t0.water(c) == CellWater::Land { (g0 * e0).clamp(emin, emax) } else { 0.0 });
        if lake { table.insert(c, 0.0); }
    }
    let e0_table = table.clone();
    // Lane TC C3 (r4_neighbour_mean_v1): cells the rule solves with land hat below r4n are not solved; neighbour-mean E.
    let r4n = crate::t1b_ls::tc_r4n(p);
    let mut r4n_cells: Vec<u32> = Vec::new();
    if r4n > 0 {
        for (&c, &(lm, tm)) in &mass { if solved[&c] && lm < r4n as f64 / 1000.0 * tm { r4n_cells.push(c); } }
        for c in &r4n_cells { solved.insert(*c, false); }
    }
    // Lane E-SOLVE-GENERALITY chart_rim_v1 (only with es coverage_min > 0): cells whose hat mass on the chart lattice is below
    // coverage_min of the full hat mass (T0 cell area / spacing^2) are rim cells: rim_mode 1 not solved (neighbour-mean E like
    // r4n), rim_mode 2 solved with the default Jacobi step and left out of the stop test.
    let es = p.es.as_ref();
    let dx_m = dxmm as f64 / 1000.0;
    let coverage: BTreeMap<u32, f64> = mass.iter().map(|(&c, &(_, tm))| (c, tm / (t0.cell_area_m2(c) / (dx_m * dx_m)))).collect();
    let mut rim_cells: Vec<u32> = Vec::new();
    // rim_mode 2: solved rim cells (default Jacobi step, outside the stop test).
    let mut rim_solved: Vec<u32> = Vec::new();
    if let Some(e) = es.filter(|e| e.coverage_min_permille > 0) {
        let rim: Vec<u32> = coverage.iter().filter(|(c, cv)| solved[c] && **cv < e.coverage_min_permille as f64 / 1000.0).map(|(c, _)| *c).collect();
        if e.rim_mode == 1 { rim_cells = rim; for c in &rim_cells { solved.insert(*c, false); } } else { rim_solved = rim; }
    }
    let nb_mean = |c: u32, table: &BTreeMap<u32, f64>, solved: &BTreeMap<u32, bool>| -> Option<f64> {
        let (i, j) = t0.ij(c);
        let mut acc = (0.0, 0);
        for di in -1i64..=1 { for dj in -1i64..=1 {
            if di == 0 && dj == 0 { continue; }
            let ii = i as i64 + di; if ii < 0 || ii >= t0.nlat as i64 { continue; }
            let jj = (j as i64 + dj).rem_euclid(t0.nlon as i64);
            let cc = (ii as usize * t0.nlon + jj as usize) as u32;
            if solved.get(&cc).copied().unwrap_or(false) { acc.0 += table[&cc]; acc.1 += 1; }
        } }
        if acc.1 > 0 { Some(acc.0 / acc.1 as f64) } else { None }
    };
    let apply_r4n = |table: &mut BTreeMap<u32, f64>| { for &c in &r4n_cells { if let Some(m) = nb_mean(c, table, &solved) { table.insert(c, m); } } };
    if r4n > 0 { apply_r4n(&mut table); }
    if !r4n_cells.is_empty() { log(format!("  lane TC r4n: cells not solved (neighbour-mean E) {:?} E {:?}", r4n_cells, r4n_cells.iter().map(|c| table[c]).collect::<Vec<_>>())); }
    let apply_rim = |table: &mut BTreeMap<u32, f64>| { for &c in &rim_cells { if let Some(m) = nb_mean(c, table, &solved) { table.insert(c, m); } } };
    if !rim_cells.is_empty() {
        apply_rim(&mut table);
        log(format!("  lane E-SOLVE-GENERALITY chart_rim_v1: cells not solved (neighbour-mean E) {:?} coverage {:?} E {:?}", rim_cells,
            rim_cells.iter().map(|c| (coverage[c] * 1000.0).round() / 1000.0).collect::<Vec<_>>(), rim_cells.iter().map(|c| table[c].round()).collect::<Vec<_>>()));
    }
    if !rim_solved.is_empty() {
        log(format!("  lane E-SOLVE-GENERALITY chart_rim_v1 mode 2: rim cells (Jacobi step, outside the stop test) {:?} coverage {:?}", rim_solved,
            rim_solved.iter().map(|c| (coverage[c] * 1000.0).round() / 1000.0).collect::<Vec<_>>()));
    }
    let mut sec: BTreeMap<u32, crate::t1b_es::SecantState> = BTreeMap::new();
    let mut kkt_stop = false;
    let mut kkt_hist: Vec<u32> = Vec::new();
    // Lane TC C4 pert_mode 2: the initial table perturbed.
    let (pert_ppm, pert_mode) = crate::t1b_ls::tc_pert(p);
    if pert_mode == 2 { for v in table.values_mut() { *v *= 1.0 + pert_ppm as f64 / 1e6; } }
    let freeze_from = crate::t1b_ls::tc_freeze(p);
    let mut frozen_any = false;
    // Lane CO: per-cell plain bilinear land fraction (coastal flags) and the flooding offsets.
    let lane_co = p.coast.is_some() || p.volcano.is_some() || p.f6_reference == 2;
    let cmin = p.coast.as_ref().map(|c| c.coastal_min_permille).unwrap_or(10) as f64 / 1000.0;
    let mut cnt: BTreeMap<u32, (u64, u64)> = BTreeMap::new();
    for i in 0..nn { let e = cnt.entry(lp.cell[i]).or_default(); e.0 += 1; if lp.bil[i] as f64 > sea && t0.water(lp.cell[i]) != CellWater::Lake { e.1 += 1; } }
    let coastal: BTreeMap<u32, bool> = cnt.iter().map(|(&c, &(n, l))| { let f = l as f64 / n as f64; (c, f > cmin && f < 1.0 - cmin) }).collect();
    let mut flood: BTreeMap<u32, f64> = BTreeMap::new();
    if let Some(cp) = &p.coast { for (&c, &k) in &coastal { if k { flood.insert(c, cp.flood_init_mm as f64 / 1000.0); } } }
    let mut hist_f4 = Vec::new();
    let mut cx = crate::t1b_coast::CoastCtx::default();
    let mut hist = Vec::new();
    let mut prev: BTreeMap<u32, (f64, f64)> = BTreeMap::new(); // (ln E, f)
    // Lane CF C7: per-cell step multipliers (E, F) and the previous residual signs.
    let damp = p.cf.as_ref().map(|c| c.solve_damp == 1).unwrap_or(false);
    let dfloor = p.cf.as_ref().map(|c| c.solve_damp_floor_permille as f64 / 1000.0).unwrap_or(1.0);
    let mut mult: BTreeMap<u32, (f64, f64)> = BTreeMap::new();
    let mut last_sign: BTreeMap<u32, (f64, f64)> = BTreeMap::new();
    let mut rows = Vec::new();
    let mut converged = false;
    let mut evals = 0;
    let mut warm: Option<Network> = None;
    let mut ls_history: Vec<[u32; 3]> = Vec::new();
    // Lane TA LS3 (only with ls e_update 1): Broyden state over the solved cells.
    let broyden_on = p.ls.as_ref().map(|l| l.e_update == 1).unwrap_or(false);
    let mut broyden = crate::t1b_ls::Broyden::default();
    // Lane TA diagnosis (env only; stdout, no bytes): per-evaluation rows of every solved cell, and an end-state probe.
    let trace = std::env::var("VZ_SOLVE_TRACE").is_ok();
    let ped_g = crate::t1b_ped::solve_ground(p);
    if ped_g { log(format!("  lane T1B-PEDESTAL: {}", crate::t1b_ped::SOLVE_RULE)); }
    let ped_b = crate::t1b_ped::local_base(p);
    if ped_b { log(format!("  lane T1B-PEDESTAL: {}", crate::t1b_ped::BASE_RULE)); }
    let mut ped_map: Option<BTreeMap<String, f64>> = None;
    let mut ped_hist: Vec<Vec<f64>> = Vec::new();
    let mut ped_last: Vec<crate::t1b_ped::LocalGround> = Vec::new();
    let mut ped_cv: Vec<Curve> = Vec::new();
    // Lane E-SOLVE-CONVERGE (esc=1 with es=1): the convergent driver chooses every evaluation's table and whether LS1 runs
    // (t1b_esc.rs); the network is held between reroutes and published as last evaluated.
    let mut incremental = crate::t1b_esc_local::prepare(p, net, nctx, t0, &extras.laws)?;
    if let Some(ir) = &incremental {
        if ir.table.keys().ne(table.keys()) { return Err("esc_ref cell table differs from this chart".into()); }
        for c in &ir.changed { if solved.get(c) != Some(&true) || rim_solved.contains(c) { return Err(format!("esc_ref: changed cell {c} is not a resolvable interior solve cell")); } }
        table = ir.table.clone(); flood = ir.flood.clone(); warm = Some(ir.net.clone());
        // Causal ablation on a fixed reference network, diagnosis only; no solve or publication changes.
        if let Some(c) = crate::t1b_esc_local::probe_cell()? {
            let rb = node_beds_warm_f(net, Some(&ir.net), nctx, t0, &extras.laws, &table, &flood, p, true);
            if let Some(st) = &rb.ls { log(format!("esc_local_probe_reference_reloaded_global_ls1 {}", serde_json::to_string(st).unwrap())); }
            for mode in [crate::t1b_ls::t1a_mode(p), 0] {
                let mut pp = p.clone();
                pp.ls.as_mut().unwrap().tb.as_mut().unwrap().t1a_levels = mode;
                let b0 = node_beds_core(&ir.net, nctx, t0, &extras.laws, &table, &flood, &pp);
                let mut tt = table.clone();
                if let Some(e) = tt.get_mut(&c) { *e *= 1.25; } else { return Err("esc local probe cell absent".into()); }
                let b1 = node_beds_core(&ir.net, nctx, t0, &extras.laws, &tt, &flood, &pp);
                let mut by_cell: BTreeMap<u32, (u32, u32, f64, f64)> = BTreeMap::new();
                let mut by_hat: BTreeMap<u32, (u32, f64)> = BTreeMap::new();
                for n in 0..NODES {
                    let dz = b1.bed[n] - b0.bed[n];
                    if dz.abs() <= 0.001 { continue; }
                    let v = by_cell.entry(net.t0_cell[n]).or_default();
                    v.0 += 1; if net.kind[n] == K_T1A { v.1 += 1; v.3 = v.3.max(dz.abs()); } v.2 = v.2.max(dz.abs());
                    let (hc, hw) = crate::t1b_net::hat4(t0, nctx.latlon[n].0, nctx.latlon[n].1);
                    for (cc, w) in hc.into_iter().zip(hw) { if w > 0.0 { let v = by_hat.entry(cc).or_default(); v.0 += 1; v.1 = v.1.max(dz.abs()); } }
                }
                log(format!("esc_local_probe source={c} factor=1.25 t1a_mode={mode} cells(changed_nodes,changed_T1a,max_abs_dbed_m,max_abs_T1a_dbed_m) {}", serde_json::to_string(&by_cell).unwrap()));
                log(format!("esc_local_probe_hat source={c} t1a_mode={mode} (changed_nodes_in_hat_support,max_abs_dbed_m) {}", serde_json::to_string(&by_hat).unwrap()));
            }
            return Err("esc_ref diagnosis-only probe complete; no layer published".into());
        }
        log(format!("esc_ref: changed {:?}, LS1 drainage scope {:?}", ir.changed, ir.scope));
    }
    let esc_p = es.and_then(|e| e.conv.as_ref());
    let mut esc = esc_p.map(|cp| {
        let sc: Vec<u32> = solved.iter().filter(|(c, s)| **s && incremental.as_ref().map(|ir| ir.changed.contains(c)).unwrap_or(true)).map(|(c, _)| *c).collect();
        crate::t1b_esc::Esc::new(cp, es.expect("esc rides in es"), &sc, &rim_solved, &|c| { let (i, j) = t0.ij(c); (i as i64, j as i64) }, t0.nlon as i64, p.solve_tolerance_mm as f64 / 1000.0, emin, emax)
    });
    if let Some(cp) = esc_p { log(format!("  lane E-SOLVE-CONVERGE: {}", cp.rule)); }
    let max_it = esc_p.map(|c| c.evals_max + 1).unwrap_or(p.solve_updates_max);
    // Diagnosis only (env VZ_ESC_PERTURB=<cell>:<dz m>): one cell's solve target raised by dz (locality test).
    let esc_pert = crate::t1b_esc::perturb_target();
    let mut esc_ref: Option<Network> = incremental.as_ref().map(|ir| ir.net.clone());
    if let Some((c, dz)) = esc_pert { log(format!("  lane E-SOLVE-CONVERGE diagnosis: VZ_ESC_PERTURB target of cell {c} raised by {dz} m")); }
    for it in 0..=max_it {
        // Lane TC C4: from evaluation freeze_from on, LS1 is not re-run (the previous evaluation's network).
        let frozen = (freeze_from > 0 && it >= freeze_from && warm.is_some()) || esc.as_ref().map(|e| !e.reroute_next() && warm.is_some()).unwrap_or(false);
        frozen_any |= frozen;
        // Lane E-SOLVE-CONVERGE ls_each: LS1 starts every evaluation from the fixed reference network (evaluation 0's).
        let start = match (&esc_ref, esc.as_ref().map(|e| e.ls_each()).unwrap_or(false)) { (Some(r), true) => Some(r), _ => warm.as_ref() };
        let mut beds = if let Some(ir) = incremental.as_mut() {
            ir.solve_beds(net, nctx, t0, &extras.laws, &table, &flood, p, log)?
        } else { node_beds_warm_f(net, start, nctx, t0, &extras.laws, &table, &flood, p, frozen) };
        if let (Some(l), Some(n2)) = (&beds.ls, &beds.net) {
            let prev = warm.as_ref().unwrap_or(net);
            let ch = (0..NODES).filter(|&k| n2.kind[k] == K_LAND && n2.recv[k] != prev.recv[k]).count() as u32;
            log(format!("  level solve: passes {} changed per pass {:?} unreached {:?} vs input {} vs previous evaluation {} of {} land; residual {} (max {:.2} m; crossing {} in hysteresis {}) top {:?}", l.passes, l.changed_per_pass, l.unreached_confined_per_pass, l.receivers_changed_vs_input, ch, l.land_nodes, l.residual_pairs, l.residual_max_m, l.residual_crossing, l.residual_in_hysteresis, &l.residual_top[..l.residual_top.len().min(4)]));
            if !l.last_changes.is_empty() || !l.residual_detail.is_empty() { log(format!("  level solve diag: last changes {:?} residual detail {:?}", l.last_changes, l.residual_detail)); }
            ls_history.push([l.passes, *l.changed_per_pass.last().unwrap_or(&0), ch]);
        }
        if let Some(c) = &beds.capture { log(format!("  capture: passes {} last {:?} changed {} frozen {} residual {} (max {:.0} m)", c.passes, c.captured_per_pass.last(), c.receivers_changed, c.frozen, c.residual_violations, c.residual_max_m)); }
        if beds.net.is_some() { warm = beds.net.take(); }
        if esc_ref.is_none() && esc.as_ref().map(|e| e.ls_each()).unwrap_or(false) { esc_ref = warm.clone(); }
        let mut beds = Beds { net: warm.clone(), ..beds };
        // Lane GLACIAL-T1: the glacier table of this evaluation's fluvial beds; land beds follow the glacial long profile.
        let gl = match (&p.glacial, co.glacial) {
            (Some(gp), Some(g)) => {
                let tb = gl_table_for(&lp, net, &beds, nctx, t0, extras, skel, chart, p, co, &table, &flood, g, gp)?;
                let kind = beds.net(net).kind.clone();
                tb.apply_beds(&kind, &mut beds.bed);
                if it == 0 || trace { log(format!("  glacial table eval {it}: {}", tb.stats)); }
                Some(std::sync::Arc::new(tb))
            }
            (Some(_), None) => return Err("glacial params without the t0.glacial layer (CoInputs::glacial)".into()),
            _ => None,
        };
        let cv = curves(net, &beds, skel, t0, p);
        // Lane T1B-PEDESTAL: the solve evaluations keep the bilinear bases (as lane VE keeps CF's cone), so the solved tables
        // equal the key-off ones; the local ground is measured on each evaluation and the final one published.
        cx = co.ctx_ped(p, t0, chart, &cv, &table, &flood, None)?;
        cx.gl = gl;
        // The network the curves belong to (lane TB: terminals for the soft channel cap; rock-mass multipliers are per node).
        let s = evaluate(&lp, beds.net(net), &cv, t0, &extras.laws, &table, p, &cx);
        evals += 1;
        // Lane T1B-PEDESTAL edifice_local_base_v1: the local ground of this evaluation (the final one is published).
        if ped_b {
            let lg = crate::t1b_ped::local_ground2(&lp.lat, &s.ground_pre, sea, &cx.volcanoes, p.cf.as_ref().ok_or("ped_base needs cf")?, crate::t1b_ve::ve(p).is_some(), crate::t1b_ped::local_relief(p));
            log(format!("  lane T1B-PEDESTAL eval {it}: local ground {:?}", lg.iter().map(|g| (g.id.clone(), g.ground_m.round(), g.base_bilinear_m.round(), g.h_m.round(), g.r_foot_m.round())).collect::<Vec<_>>()));
            ped_map = Some(lg.iter().map(|g| (g.id.clone(), g.ground_m)).collect());
            ped_hist.push(lg.iter().map(|g| g.ground_m).collect());
            ped_last = lg;
            ped_cv = cv;
        }
        let mut acc: BTreeMap<u32, [f64; 4]> = BTreeMap::new();
        // Plain per-cell sums: [posts, land posts, sum z, sum ref, sum bil] (land = after flooding; edifices excluded from the means).
        let mut plain: BTreeMap<u32, [f64; 6]> = BTreeMap::new();
        // Lane TB S3: plain per-cell T1a trunk-level sums over the same posts (only with the option).
        let r4 = crate::t1b_ls::r4_permille(p);
        let es_target = es.map(|e| e.target == 1).unwrap_or(false);
        let mut plain_tr: BTreeMap<u32, f64> = BTreeMap::new();
        let refv_of = |i: usize| -> f64 {
            let cm = t0.elevation_m[lp.cell[i] as usize] as f64;
            match p.f6_reference { 1 => cm, 2 => if coastal.get(&lp.cell[i]).copied().unwrap_or(false) && t0.water(lp.cell[i]) == CellWater::Land { cm } else { lp.bil[i] as f64 }, _ => lp.bil[i] as f64 }
        };
        for i in 0..nn {
            let land = if p.coast.is_some() { s.fixed[i] != MARINE && s.fixed[i] != LAKE } else { lp.water[i] == FREE };
            // Lane T1B-PEDESTAL edifice_ground_in_solve_v1: an edifice post counts at its pre-edifice ground (land ground only).
            let zg = if ped_g && s.over[i] == crate::t1b_coast::OVER_EDIFICE && s.ground_pre.len() == nn && s.ground_pre[i].is_finite() && s.ground_pre[i] > sea { Some(s.ground_pre[i]) } else { None };
            if lane_co {
                let e = plain.entry(lp.cell[i]).or_default();
                e[0] += 1.0;
                if land && s.z[i].is_finite() { e[1] += 1.0; if s.over[i] != crate::t1b_coast::OVER_EDIFICE || zg.is_some() { e[2] += zg.unwrap_or(s.z[i]); e[3] += refv_of(i); e[4] += lp.bil[i] as f64; e[5] += 1.0; if r4 > 0 || es_target { *plain_tr.entry(lp.cell[i]).or_default() += trunk[i]; } } }
            }
            if !land || !s.z[i].is_finite() || (s.over[i] == crate::t1b_coast::OVER_EDIFICE && zg.is_none()) { continue; }
            let zi = zg.unwrap_or(s.z[i]);
            for k in 0..4 {
                let w = lp.hat_w[i][k] as f64;
                if w <= 0.0 { continue; }
                let a = acc.entry(lp.hat_c[i][k]).or_default();
                let refv = refv_of(i);
                a[0] += w; a[1] += w * zi; a[2] += w * refv; a[3] += w * trunk[i];
            }
        }
        rows.clear();
        let mut maxd: f64 = 0.0;
        for (&c, &sv) in &solved {
            let a = acc.get(&c).copied().unwrap_or_default();
            let (m, tg, l) = if a[0] > 0.0 { (a[1] / a[0], a[2] / a[0], a[3] / a[0]) } else { (0.0, 0.0, 0.0) };
            let diff = m - tg;
            let (lm, tm) = mass[&c];
            let mut row = CellRow { cell: c, solved: sv, land_hat_fraction: lm / tm, e0_m_per_my: e0_table[&c], e_m_per_my: table[&c], hat_mean_z_m: m, hat_mean_bilinear_m: tg, hat_mean_trunk_m: l, diff_m: diff, at_bound: table[&c] <= emin || table[&c] >= emax,
                coastal: None, flood_offset_m: None, land_fraction: None, land_fraction_bilinear: None, cell_mean_z_m: None, cell_mean_ref_m: None, cell_mean_bilinear_m: None, cell_land_posts: None,
                r4_plain_reference: None, cell_mean_trunk_m: None, hat_diff_m: None, r4n_neighbour_mean: if r4n_cells.contains(&c) { Some(true) } else { None }, es_plain_target: None };
            if lane_co {
                let pl = plain.get(&c).copied().unwrap_or_default();
                let (np, _) = cnt.get(&c).copied().unwrap_or_default();
                row.coastal = Some(coastal.get(&c).copied().unwrap_or(false));
                row.flood_offset_m = flood.get(&c).copied().or(Some(0.0));
                row.land_fraction = Some(if pl[0] > 0.0 { pl[1] / pl[0] } else { 0.0 });
                row.land_fraction_bilinear = Some(if np > 0 { cnt[&c].1 as f64 / np as f64 } else { 0.0 });
                row.cell_mean_z_m = Some(if pl[5] > 0.0 { pl[2] / pl[5] } else { 0.0 });
                row.cell_mean_ref_m = Some(if pl[5] > 0.0 { pl[3] / pl[5] } else { 0.0 });
                row.cell_mean_bilinear_m = Some(if pl[5] > 0.0 { pl[4] / pl[5] } else { 0.0 });
                row.cell_land_posts = Some(pl[5] as u64);
                // Lane TB S3 (R4 measurement): marine cells under the land-hat threshold use the plain per-cell means.
                if r4 > 0 && sv && t0.water(c) == CellWater::Marine && lm / tm < r4 as f64 / 1000.0 && pl[5] > 0.0 {
                    let trm = plain_tr.get(&c).copied().unwrap_or(0.0) / pl[5];
                    row.r4_plain_reference = Some(true);
                    row.cell_mean_trunk_m = Some(trm);
                    row.hat_diff_m = Some(diff);
                    row.diff_m = pl[2] / pl[5] - pl[3] / pl[5];
                }
                // Lane E-SOLVE-GENERALITY plain_land_target_v1: resolvable Land cells solve the plain whole-cell means.
                if es_target && sv && t0.water(c) == CellWater::Land && !rim_solved.contains(&c) && pl[5] > 0.0
                    && coverage[&c] >= es.map(|e| e.coverage_min_permille).unwrap_or(0) as f64 / 1000.0 {
                    let trm = plain_tr.get(&c).copied().unwrap_or(0.0) / pl[5];
                    row.es_plain_target = Some(true);
                    row.cell_mean_trunk_m = Some(trm);
                    row.hat_diff_m = Some(diff);
                    row.diff_m = pl[2] / pl[5] - pl[3] / pl[5];
                }
            }
            if let Some((_, dz)) = esc_pert.filter(|q| q.0 == c) { row.diff_m -= dz; row.hat_mean_bilinear_m += dz; if let Some(v) = row.cell_mean_ref_m.as_mut() { *v += dz; } }
            if sv { maxd = maxd.max(row.diff_m.abs()); }
            rows.push(row);
        }
        if trace { for r in rows.iter().filter(|r| r.solved) { log(format!("  trace eval {it} cell {} E {:.3} M {:.2} T {:.2} L {:.2} diff {:.2} landhat {:.3} F {:.1}", r.cell, r.e_m_per_my, r.hat_mean_z_m, r.hat_mean_bilinear_m, r.hat_mean_trunk_m, r.diff_m, r.land_hat_fraction, r.flood_offset_m.unwrap_or(0.0))); } }
        let maxf4 = rows.iter().filter(|r| r.coastal == Some(true) && r.solved).map(|r| (r.land_fraction.unwrap_or(0.0) - r.land_fraction_bilinear.unwrap_or(0.0)).abs()).fold(0.0, f64::max);
        hist_f4.push(maxf4);
        if let Some(cp) = &p.coast { log(format!("  coast: max |land - bilinear land| over solved coastal cells {maxf4:.3}; fans {}; flooded {} fan-raised {} cliff-cut {} beach {} edifice {} buried {}", cx.fans.len(), s.coast.flooded_posts, s.coast.fan_raised_posts, s.coast.cliff_cut_posts, s.coast.beach_posts, s.coast.edifice_posts, s.coast.buried_channel_posts)); let _ = cp; }
        hist.push(maxd);
        let worst = rows.iter().filter(|r| r.solved).max_by(|a, b| a.diff_m.abs().total_cmp(&b.diff_m.abs()).then(b.cell.cmp(&a.cell))).map(|r| (r.cell, r.diff_m.round(), r.e_m_per_my.round())).unwrap_or_default();
        let nbad = rows.iter().filter(|r| r.solved && r.diff_m.abs() > p.solve_tolerance_mm as f64 / 1000.0).count();
        log(format!("solve eval {it}: max |M - T| {maxd:.1} m (worst cell, diff, E {worst:?}; {nbad} cells over tolerance); crest lowered {} (max {:.0} m) fill raised {} (max {:.1} m); seconds raster/floor/fmm/crest/fill {:?}; fmm unreached {}",
            s.stats.crest_lowered_posts, s.stats.crest_lowered_max_m, s.stats.fill_raised_posts, s.stats.fill_raised_max_m, s.stats.seconds.map(|v| (v * 10.0).round() / 10.0), s.stats.fmm_unreached));
        let f4_ok = p.coast.as_ref().map(|cp| maxf4 <= cp.flood_tolerance_permille as f64 / 1000.0).unwrap_or(true);
        if esc.is_none() && maxd <= p.solve_tolerance_mm as f64 / 1000.0 && f4_ok { converged = true; break; }
        let mut esc_f4_open = false;
        // Lane E-SOLVE-GENERALITY bound_kkt_stop_v1: within tolerance or pinned at a bound with an outward residual.
        if let Some(e) = es {
            let tol = p.solve_tolerance_mm as f64 / 1000.0;
            let open: Vec<&CellRow> = rows.iter().filter(|r| r.solved && !rim_solved.contains(&r.cell) && !(r.diff_m.abs() <= tol || (r.e_m_per_my <= emin && r.diff_m > 0.0) || (r.e_m_per_my >= emax && r.diff_m < 0.0))).collect();
            kkt_hist.push(open.len() as u32);
            log(format!("  lane E-SOLVE-GENERALITY eval {it}: cells open (not within tolerance, not bound-limited) {} {:?}", open.len(), open.iter().map(|r| (r.cell, r.diff_m.round(), r.e_m_per_my.round())).collect::<Vec<_>>()));
            // F4 by the same rule: a coastal cell within the F4 tolerance, or with its flooding offset F at a bound
            // (F = flood_max and land in excess, or F = 0 and land short).
            let f4_open: Vec<(u32, f64, f64)> = match &p.coast {
                Some(cp) => {
                    let (ftol, fmax) = (cp.flood_tolerance_permille as f64 / 1000.0, cp.flood_max_mm as f64 / 1000.0);
                    rows.iter().filter(|r| r.coastal == Some(true) && r.solved).filter_map(|r| {
                        let ex = r.land_fraction.unwrap_or(0.0) - r.land_fraction_bilinear.unwrap_or(0.0);
                        let fv = r.flood_offset_m.unwrap_or(0.0);
                        if ex.abs() <= ftol || (fv >= fmax && ex > 0.0) || (fv <= 0.0 && ex < 0.0) { None } else { Some((r.cell, ex, fv)) }
                    }).collect()
                }
                None => Vec::new(),
            };
            if !f4_ok { log(format!("  lane E-SOLVE-GENERALITY eval {it}: F4 open (not within tolerance, F not at a bound) {:?}", f4_open)); }
            esc_f4_open = !f4_open.is_empty();
            if esc.is_none() && e.kkt == 1 && open.is_empty() && f4_open.is_empty() { kkt_stop = true; break; }
        }
        // Lane E-SOLVE-CONVERGE: the driver's next table (or the stop; the table just evaluated is then the published one).
        let mut esc_next: Option<BTreeMap<u32, f64>> = None;
        let mut esc_flood = true;
        if let Some(ec) = esc.as_mut() {
            let obs: Vec<crate::t1b_esc::CellObs> = rows.iter().filter(|r| r.solved).map(|r| {
                let (num, den) = if r.r4_plain_reference == Some(true) || r.es_plain_target == Some(true) {
                    let trm = r.cell_mean_trunk_m.unwrap_or(0.0);
                    (r.cell_mean_ref_m.unwrap_or(0.0) - trm, r.cell_mean_z_m.unwrap_or(0.0) - trm)
                } else { (r.hat_mean_bilinear_m - r.hat_mean_trunk_m, r.hat_mean_z_m - r.hat_mean_trunk_m) };
                crate::t1b_esc::CellObs { cell: r.cell, diff_m: r.diff_m, e: r.e_m_per_my, num, den }
            }).collect();
            match ec.after_eval(&obs, &table, &flood, if incremental.is_some() { rows.iter().filter(|r| incremental.as_ref().unwrap().changed.contains(&r.cell) && r.coastal == Some(true)).all(|r| (r.land_fraction.unwrap_or(0.0) - r.land_fraction_bilinear.unwrap_or(0.0)).abs() <= p.coast.as_ref().map(|c| c.flood_tolerance_permille as f64 / 1000.0).unwrap_or(0.0)) } else { f4_ok || !esc_f4_open }, log) {
                crate::t1b_esc::Act::Stop => { converged = maxd <= p.solve_tolerance_mm as f64 / 1000.0 && f4_ok; kkt_stop = ec.converged_kkt(); break; }
                crate::t1b_esc::Act::Next { table: t, update_flood, flood: f, .. } => { esc_next = Some(t); esc_flood = update_flood; if let Some(f) = f { flood = f; } }
            }
            if it == max_it { log("  lane E-SOLVE-CONVERGE: loop cap reached before the driver stopped (published table = last evaluated)".into()); }
        }
        if it == max_it { break; }
        if let Some(cp) = p.coast.as_ref().filter(|_| esc_flood) {
            let (g, smax, fmax) = (cp.flood_gain_mm as f64 / 1000.0, cp.flood_step_max_mm as f64 / 1000.0, cp.flood_max_mm as f64 / 1000.0);
            for r in &rows {
                if r.coastal != Some(true) || incremental.as_ref().map(|ir| !ir.changed.contains(&r.cell)).unwrap_or(false) { continue; }
                let f = flood.entry(r.cell).or_insert(0.0);
                let ex = r.land_fraction.unwrap_or(0.0) - r.land_fraction_bilinear.unwrap_or(0.0);
                // Lane TG G4: only the error outside the inner band moves F.
                let ex = match crate::t1b_tg::tg(p).filter(|t| t.f_deadband == 1) { Some(t) => crate::t1b_tg::deadband(ex, t.f_deadband_permille as f64 / 1000.0), None => ex };
                let mf = if damp {
                    let m = mult.entry(r.cell).or_insert((1.0, 1.0));
                    let ls = last_sign.entry(r.cell).or_insert((0.0, 0.0));
                    let sg = if ex.abs() < 1e-9 { 0.0 } else { ex.signum() };
                    if sg * ls.1 < 0.0 { m.1 = (m.1 * 0.5).max(dfloor); }
                    if sg != 0.0 { ls.1 = sg; }
                    m.1
                } else { 1.0 };
                let step = (mf * g * ex).clamp(-smax, smax);
                *f = (*f + step).clamp(0.0, fmax);
            }
        }
        if let Some(t) = esc_next.take() {
            table = t;
            if incremental.is_none() {
                if r4n > 0 { apply_r4n(&mut table); }
                if !rim_cells.is_empty() { apply_rim(&mut table); }
            }
            continue;
        }
        let mut next = table.clone();
        if broyden_on {
            let sv: Vec<&CellRow> = rows.iter().filter(|r| r.solved).collect();
            if broyden.cells.is_empty() { broyden.cells = sv.iter().map(|r| r.cell).collect(); }
            let x: Vec<f64> = sv.iter().map(|r| table[&r.cell].ln()).collect();
            let rr: Vec<f64> = sv.iter().map(|r| r.diff_m).collect();
            let d: Vec<f64> = sv.iter().map(|r| (r.hat_mean_z_m - r.hat_mean_trunk_m).max(1.0)).collect();
            let st = broyden.step(&x, &rr, &d);
            for (k, r) in sv.iter().enumerate() { next.insert(r.cell, (x[k] + st[k]).exp().clamp(emin, emax)); }
            table = next;
            if r4n > 0 { apply_r4n(&mut table); }
            if !rim_cells.is_empty() { apply_rim(&mut table); }
            continue;
        }
        // Damped Jacobi in ln E (REVISION-free default): step = omega x ln((T - L)/(M - L)), clamped to +-ln 2;
        // omega = 0.6 (neighbouring cells share hats and receiver paths, so undamped steps oscillate: run s42_r0).
        for r in &rows {
            if !r.solved { continue; }
            let e = table[&r.cell];
            let (num, den) = if r.r4_plain_reference == Some(true) || r.es_plain_target == Some(true) {
                let trm = r.cell_mean_trunk_m.unwrap_or(0.0);
                (r.cell_mean_ref_m.unwrap_or(0.0) - trm, r.cell_mean_z_m.unwrap_or(0.0) - trm)
            } else { (r.hat_mean_bilinear_m - r.hat_mean_trunk_m, r.hat_mean_z_m - r.hat_mean_trunk_m) };
            // Lane E-SOLVE-GENERALITY secant_diagonal_v1 (only with es update 1).
            if let Some(ep) = es.filter(|x| x.update == 1 && !rim_solved.contains(&r.cell)) {
                let (st, ns) = crate::t1b_es::secant_step(ep, sec.get(&r.cell), e.ln(), num, den, r.diff_m);
                if let Some(ns) = ns { sec.insert(r.cell, ns); }
                next.insert(r.cell, (e.ln() + st).exp().clamp(emin, emax));
                prev.insert(r.cell, (e.ln(), r.diff_m));
                continue;
            }
            let step = if num > 1.0 && den > 1.0 { 0.6 * (num / den).ln() } else if r.diff_m < 0.0 { 2f64.ln() } else { -(2f64.ln()) };
            let step = step.clamp(-(2f64.ln()), 2f64.ln());
            let step = if damp {
                let m = mult.entry(r.cell).or_insert((1.0, 1.0));
                let ls = last_sign.entry(r.cell).or_insert((0.0, 0.0));
                let sg = if r.diff_m.abs() < 1e-9 { 0.0 } else { r.diff_m.signum() };
                if sg * ls.0 < 0.0 { m.0 = (m.0 * 0.5).max(dfloor); }
                if sg != 0.0 { ls.0 = sg; }
                step * m.0
            } else { step };
            next.insert(r.cell, (e.ln() + step).exp().clamp(emin, emax));
            prev.insert(r.cell, (e.ln(), r.diff_m));
        }
        table = next;
        if r4n > 0 { apply_r4n(&mut table); }
        if !rim_cells.is_empty() { apply_rim(&mut table); }
    }
    // Lane T1B-PEDESTAL: the published volcano sites stand on the final evaluation's local ground (what windows read).
    let ped_rec = if ped_b {
        let cx2 = co.ctx_ped(p, t0, chart, &ped_cv, &table, &flood, ped_map.as_ref())?;
        cx.volcanoes = cx2.volcanoes; cx.t1a_segs = cx2.t1a_segs; cx.tf_vent_c1 = cx2.tf_vent_c1;
        Some(serde_json::json!({"rule": crate::t1b_ped::BASE_RULE, "label": crate::t1b_ped::LABEL, "ground": ped_last, "ground_history_m": ped_hist,
            "note": "ground_m = G measured on the final E-solve evaluation (the published table; the solve evaluations themselves keep the bilinear bases, so the table equals the key-off one); windows read it; the volcano sites below are built with it"}))
    } else if ped_g { Some(serde_json::json!({"rule": crate::t1b_ped::SOLVE_RULE, "label": crate::t1b_ped::LABEL})) } else { None };
    // Lane T1B-PEDESTAL measurement (env VZ_PEDESTAL_DUMP only; no bytes): one extra evaluation of the published table with the
    // published (local-ground) volcano sites, dumped by evaluate's hook (whole-cell means with the new cones).
    if ped_b && crate::t1b_ped::dump_dir().is_some() {
        let _ = evaluate(&lp, warm.as_ref().unwrap_or(net), &ped_cv, t0, &extras.laws, &table, p, &cx);
        log("  lane T1B-PEDESTAL: measurement evaluation with the published sites dumped (VZ_PEDESTAL_DUMP)".into());
    }
    // Lane GLACIAL-T1: the published glacier table on the published network and beds (as the caller recomputes them).
    if let (Some(gp), Some(g)) = (&p.glacial, co.glacial) {
        let fb = node_beds_warm_f(net, warm.as_ref(), nctx, t0, &extras.laws, &table, &flood, p, frozen_any);
        let tb = gl_table_for(&lp, net, &fb, nctx, t0, extras, skel, chart, p, co, &table, &flood, g, gp)?;
        log(format!("  glacial table (published): {}", tb.stats));
        cx.gl = Some(std::sync::Arc::new(tb));
    }
    // Lane TC C4 pert_mode 1: the published table perturbed (the caller recomputes the network from the warm start).
    if pert_mode == 1 { for v in table.values_mut() { *v *= 1.0 + pert_ppm as f64 / 1e6; } log(format!("  lane TC: published table x (1 + {pert_ppm} ppm)")); }
    // Lane TA diagnosis (env VZ_SOLVE_PROBE=c1,c2,...; stdout only): response of every solved cell's hat mean to +ln 1.25 in E
    // of each listed cell, from the final table (flooding offsets held).
    if let Ok(list) = std::env::var("VZ_SOLVE_PROBE") {
        // Lane E-SOLVE-CONVERGE diagnosis: VZ_SOLVE_PROBE_FROZEN holds the network (LS1 not re-run); list "all" = every solved cell.
        let frz = std::env::var("VZ_SOLVE_PROBE_FROZEN").is_ok();
        let list = if list == "all" { solved.iter().filter(|(_, s)| **s).map(|(c, _)| c.to_string()).collect::<Vec<_>>().join(",") } else { list };
        let means = |tab: &BTreeMap<u32, f64>, warm: Option<&Network>| -> R<BTreeMap<u32, (f64, f64, f64)>> {
            let beds = node_beds_warm_f(net, warm, nctx, t0, &extras.laws, tab, &flood, p, frz && warm.is_some());
            let cv = curves(net, &beds, skel, t0, p);
            let cx = co.ctx(p, t0, chart, &cv, tab, &flood)?;
            let s = evaluate(&lp, beds.net(net), &cv, t0, &extras.laws, tab, p, &cx);
            let mut acc: BTreeMap<u32, [f64; 4]> = BTreeMap::new();
            for i in 0..nn {
                let land = if p.coast.is_some() { s.fixed[i] != MARINE && s.fixed[i] != LAKE } else { lp.water[i] == FREE };
                if !land || !s.z[i].is_finite() || s.over[i] == crate::t1b_coast::OVER_EDIFICE { continue; }
                let cm = t0.elevation_m[lp.cell[i] as usize] as f64;
                let refv = match p.f6_reference { 1 => cm, 2 => if coastal.get(&lp.cell[i]).copied().unwrap_or(false) && t0.water(lp.cell[i]) == CellWater::Land { cm } else { lp.bil[i] as f64 }, _ => lp.bil[i] as f64 };
                for k in 0..4 { let w = lp.hat_w[i][k] as f64; if w <= 0.0 { continue; } let a = acc.entry(lp.hat_c[i][k]).or_default(); a[0] += w; a[1] += w * s.z[i]; a[2] += w * refv; a[3] += w * trunk[i]; }
            }
            Ok(acc.into_iter().filter(|(c, _)| solved.get(c).copied().unwrap_or(false)).map(|(c, a)| (c, (a[1] / a[0], a[2] / a[0], a[3] / a[0]))).collect())
        };
        let base = means(&table, warm.as_ref())?;
        for (c, (m, t, l)) in &base { log(format!("  probe base cell {c} E {:.3} M {m:.2} T {t:.2} L {l:.2} (last evaluation M {:.2}; frozen {frz})", table[c], rows.iter().find(|r| r.cell == *c).map(|r| r.hat_mean_z_m).unwrap_or(f64::NAN))); }
        for c in list.split(',').filter_map(|v| v.trim().parse::<u32>().ok()) {
            let mut tab = table.clone();
            if let Some(e) = tab.get_mut(&c) { *e *= 1.25; } else { continue; }
            let m2 = means(&tab, warm.as_ref())?;
            for (k, (m, _, _)) in &m2 { log(format!("  probe dlnE({c}) = ln 1.25: cell {k} dM {:.2} m (dM/dlnE {:.1} m)", m - base[k].0, (m - base[k].0) / 1.25f64.ln())); }
        }
        // Lane E-SOLVE-CONVERGE diagnosis (env VZ_SOLVE_PROBE_SCAN=c:s1,s2,...): one cell's ln E stepped by each s; every cell's dM.
        if let Some((c, ss)) = std::env::var("VZ_SOLVE_PROBE_SCAN").ok().and_then(|v| v.split_once(':').map(|(a, b)| (a.trim().parse::<u32>().ok(), b.to_string()))).and_then(|(a, b)| a.map(|a| (a, b))) {
            for sv in ss.split(',').filter_map(|v| v.trim().parse::<f64>().ok()) {
                let mut tab = table.clone();
                if let Some(e) = tab.get_mut(&c) { *e *= sv.exp(); } else { continue; }
                let m2 = means(&tab, warm.as_ref())?;
                log(format!("  probe scan dlnE({c}) = {sv}: {}", m2.iter().map(|(k, (m, _, _))| format!("{k}:{:.1}", m - base[k].0)).collect::<Vec<_>>().join(" ")));
            }
        }
    }
    // Lane E-SOLVE-GENERALITY: the solve record (rim cells, bound-limited cells, secant g, ADR 043 F6 literal on Land cells).
    let es_rec = es.map(|e| {
        let tol = p.solve_tolerance_mm as f64 / 1000.0;
        let f_bound: Vec<Value> = match &p.coast {
            Some(cp) => {
                let (ftol, fmax) = (cp.flood_tolerance_permille as f64 / 1000.0, cp.flood_max_mm as f64 / 1000.0);
                rows.iter().filter(|r| r.coastal == Some(true) && r.solved).filter_map(|r| {
                    let ex = r.land_fraction.unwrap_or(0.0) - r.land_fraction_bilinear.unwrap_or(0.0);
                    let fv = r.flood_offset_m.unwrap_or(0.0);
                    if ex.abs() > ftol && ((fv >= fmax && ex > 0.0) || (fv <= 0.0 && ex < 0.0)) { Some(serde_json::json!({"cell": r.cell, "land_fraction_excess": ex, "flood_offset_m": fv})) } else { None }
                }).collect()
            }
            None => Vec::new(),
        };
        let rim2: Vec<Value> = rows.iter().filter(|r| rim_solved.contains(&r.cell)).map(|r| serde_json::json!({"cell": r.cell, "water": format!("{:?}", t0.water(r.cell)), "coverage": coverage[&r.cell], "e_m_per_my": r.e_m_per_my, "diff_m": r.diff_m, "within_tolerance": r.diff_m.abs() <= tol})).collect();
        let bound_limited: Vec<Value> = rows.iter().filter(|r| r.solved && !rim_solved.contains(&r.cell) && r.diff_m.abs() > tol && ((r.e_m_per_my <= emin && r.diff_m > 0.0) || (r.e_m_per_my >= emax && r.diff_m < 0.0)))
            .map(|r| serde_json::json!({"cell": r.cell, "water": format!("{:?}", t0.water(r.cell)), "e_m_per_my": r.e_m_per_my, "diff_m": r.diff_m, "hat_mean_z_m": r.hat_mean_z_m, "hat_mean_bilinear_m": r.hat_mean_bilinear_m, "hat_mean_trunk_m": r.hat_mean_trunk_m, "plain_target": r.es_plain_target == Some(true),
                "target_below_trunk": if r.es_plain_target == Some(true) { r.cell_mean_ref_m.unwrap_or(0.0) < r.cell_mean_trunk_m.unwrap_or(0.0) } else { r.hat_mean_bilinear_m < r.hat_mean_trunk_m }, "coverage": coverage[&r.cell]})).collect();
        let rim: Vec<Value> = rim_cells.iter().map(|c| serde_json::json!({"cell": c, "water": format!("{:?}", t0.water(*c)), "coverage": coverage[c], "e_m_per_my": table[c]})).collect();
        let g: BTreeMap<String, f64> = sec.iter().map(|(c, s)| (c.to_string(), s.g)).collect();
        let mut f6 = Vec::new();
        let (mut f6_n, mut f6_pass) = (0u32, 0u32);
        for r in &rows {
            if t0.water(r.cell) != CellWater::Land { continue; }
            let posts = r.cell_land_posts.unwrap_or(0);
            if posts == 0 { continue; }
            let cov = coverage[&r.cell];
            let ftol = (0.1 * t0.subgrid_relief_m[r.cell as usize] as f64).max(50.0);
            let d = r.cell_mean_z_m.unwrap_or(0.0) - r.cell_mean_bilinear_m.unwrap_or(0.0);
            let pass = d.abs() <= ftol;
            let full = cov >= 0.5;
            if full { f6_n += 1; if pass { f6_pass += 1; } }
            f6.push(serde_json::json!({"cell": r.cell, "coverage": cov, "centre_in_chart": full, "land_posts": posts, "cell_mean_z_m": r.cell_mean_z_m, "cell_mean_bilinear_m": r.cell_mean_bilinear_m, "t0_mean_m": t0.elevation_m[r.cell as usize], "diff_m": d, "tolerance_m": ftol, "pass": pass, "solved": r.solved}));
        }
        let mut v = serde_json::json!({"rule": crate::t1b_es::ES_RULE, "update": e.update, "coverage_min_permille": e.coverage_min_permille, "kkt": e.kkt,
            "converged_strict": converged, "kkt_stop": kkt_stop, "evals": evals, "open_cells_history": kkt_hist,
            "rim_mode": e.rim_mode, "target": e.target, "rim_cells": rim, "rim_cells_solved": rim2, "bound_limited": bound_limited, "f4_flood_bound_limited": f_bound, "secant_g": g,
            "f6_cell_literal": {"rule": "ADR 043 F6 per T0 cell: plain mean of the layer surface over the cell's own land posts on the solve lattice (edifices excluded) vs the T0 bilinear over the same posts, tolerance max(50 m, 10 % of T0 subgrid relief); Land cells with land posts; counted when the cell centre is on the chart (coverage >= 0.5); the layer surface is before window-only keys", "counted": f6_n, "pass": f6_pass, "cells": f6},
            "note": "EXPLORATORY opt-in (lane E-SOLVE-GENERALITY); the strict converged flag is unchanged; bound_limited cells are reported, never hidden"});
        // Lane E-SOLVE-CONVERGE: the driver record (only with esc=1; es-only headers unchanged).
        if let Some(ec) = &esc { v["conv"] = ec.record(); }
        if let Some(ir) = &incremental { v["incremental"] = ir.record(warm.as_ref().unwrap_or(net), &table); }
        v
    });
    if let Some(ec) = &esc { let v = ec.record(); log(format!("  lane E-SOLVE-CONVERGE: stop {} evals {} rounds {} converged_kkt {} held {}", v["stop_reason"], v["evals"], v["rounds"], v["converged_kkt"], v["held_cells"])); }
    if let Some(r) = &es_rec { log(format!("  lane E-SOLVE-GENERALITY: kkt_stop {} strict {} evals {} rim {} bound-limited {} F6 literal {}/{}", kkt_stop, converged, evals, rim_cells.len(), r["bound_limited"].as_array().map(|a| a.len()).unwrap_or(0), r["f6_cell_literal"]["pass"], r["f6_cell_literal"]["counted"])); }
    Ok(SolveResult { table, rows, history_max_abs_diff_m: hist, converged, evals, seconds: t_all.elapsed().as_secs_f64(), flood, cx, history_max_abs_f4: hist_f4, net: warm, ls_history, frozen: frozen_any || esc.is_some(), te: None, es: es_rec, ped: ped_rec })
}

/// Lane GLACIAL-T1: the glacier table of fluvial `beds`: one evaluation WITHOUT the carve on the solve lattice gives the
/// pre-glacial surface of each node's own posts (mass balance and node floors), then the chart table.
#[allow(clippy::too_many_arguments)]
fn gl_table_for(lp: &LatPrep, net: &Network, beds: &Beds, nctx: &NodeCtx, t0: &T0, extras: &T0Extras, skel: &SkeletonView, chart: &Chart, p: &RiverNetParamsV1, co: &CoInputs,
    table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, g: &crate::t1b_glacial::T0Glacial, gp: &crate::t1b_glacial::GlacialParamsV1) -> R<crate::t1b_glacial::GlacialTable> {
    let cv0 = curves(net, beds, skel, t0, p);
    let cx0 = co.ctx(p, t0, chart, &cv0, table, flood)?;
    let s0 = evaluate(lp, beds.net(net), &cv0, t0, &extras.laws, table, p, &cx0);
    let zp = crate::t1b_glacial::node_surface(&lp.lat, &s0.z, &s0.fixed);
    let en: Vec<f64> = nctx.latlon.iter().map(|&(la, lo)| e_at(t0, table, la, lo)).collect();
    crate::t1b_glacial::build_table(beds.net(net), &beds.bed, &beds.chan, &zp, &en, &nctx.latlon, t0, g, extras.laws.stream_power_m, gp)
}

/// Lane GLACIER-1 (gl_now=1, called once by the layer build after the solve): the icefield sheet of the PUBLISHED layer on
/// the solve lattice: one evaluation of the published table, network and beds with the volcano sites the windows build
/// (T1B-PEDESTAL: on the published local ground, read from the solve record), without a carve (present ice only); the
/// sheet goes into the published glacier table (t1b_glacial::build_sheet).
#[allow(clippy::too_many_arguments)]
pub fn gl_now_sheet(t0: &T0, extras: &T0Extras, skel: &SkeletonView, chart: &Chart, net: &Network, nctx: &NodeCtx, p: &RiverNetParamsV1, strata: Option<&StrataCols>, anchor: (f64, f64), co: &CoInputs,
    sol: &mut SolveResult, g: &crate::t1b_glacial::T0Glacial) -> R<String> {
    let gp = p.glacial.as_ref().ok_or("gl_now_sheet without glacial params")?;
    if gp.now.is_none() || gp.present_ice_only != 1 { return Err("gl_now_sheet needs gl_now=1 on a present-ice-only layer".into()); }
    let dxmm = p.solve_spacing_mm;
    let half = vz_tilekey::CHART_HALF_EXTENT_MM / dxmm;
    let lat = EvalLat { gx0: -half, gy0: -half, n: (2 * half + 1) as usize, dx_mm: dxmm };
    let lp = lat_prep(lat, chart, t0, extras, p, strata, anchor)?;
    let beds = node_beds_warm_f(net, sol.net.as_ref(), nctx, t0, &extras.laws, &sol.table, &sol.flood, p, sol.frozen);
    let cv = curves(net, &beds, skel, t0, p);
    let ped: Option<BTreeMap<String, f64>> = sol.ped.as_ref().and_then(|v| v.get("ground")).and_then(|v| v.as_array())
        .map(|a| a.iter().filter_map(|x| Some((x.get("id")?.as_str()?.to_string(), x.get("ground_m")?.as_f64()?))).collect());
    let cx = co.ctx_ped(p, t0, chart, &cv, &sol.table, &sol.flood, ped.as_ref())?;
    let s = evaluate(&lp, beds.net(net), &cv, t0, &extras.laws, &sol.table, p, &cx);
    let mut tb = sol.cx.gl.as_ref().map(|t| (**t).clone()).ok_or("gl_now_sheet: the solve returned no glacier table")?;
    // Lane GLACIER-2 (gl_bal=1): the balance-flux icefield on the same evaluation (sea floor = deeper of z_eik and the T0 bilinear).
    if gp.now.as_ref().map(|n| n.bal.is_some()).unwrap_or(false) { crate::t1b_glacial::build_sheet_bal(&mut tb, &lp.lat, &s.z, &s.fixed, &s.z_eik, &lp.bil, chart, t0, g, gp)?; }
    else { crate::t1b_glacial::build_sheet(&mut tb, beds.net(net), &lp.lat, &s.z, &s.fixed, chart, t0, g, gp)?; }
    let st = tb.now.as_ref().map(|n| n.stats.to_string()).unwrap_or_default();
    sol.cx.gl = Some(std::sync::Arc::new(tb));
    Ok(format!("pedestal ground {}; {st}", if ped.is_some() { "published (T1B-PEDESTAL)" } else { "none" }))
}

/// Lane MS-CEIL (ceil only; called once by the layer build after the solve): the lane_ceil header record. One evaluation
/// of the PUBLISHED table, network and beds on the solve lattice (as gl_now_sheet: volcano sites on the published local
/// ground), then per T0 cell: land posts, edifice / ceiling / active counts, the active share a = (land - edifice -
/// ceiling) / (land - edifice), mean c_h R over the cell's land posts, the record's R, E_v (table), held-at-bound and the
/// whole-cell F6 literal of lane_es where present. A record only: no published byte other than the header depends on it.
#[allow(clippy::too_many_arguments)]
pub fn ceil_record(t0: &T0, extras: &T0Extras, skel: &SkeletonView, chart: &Chart, net: &Network, nctx: &NodeCtx, p: &RiverNetParamsV1, strata: Option<&StrataCols>, anchor: (f64, f64), co: &CoInputs,
    sol: &SolveResult) -> R<Value> {
    let cp = p.ceil.as_ref().ok_or("ceil_record without ceil params")?;
    let dxmm = p.solve_spacing_mm;
    let half = vz_tilekey::CHART_HALF_EXTENT_MM / dxmm;
    let lat = EvalLat { gx0: -half, gy0: -half, n: (2 * half + 1) as usize, dx_mm: dxmm };
    let lp = lat_prep(lat, chart, t0, extras, p, strata, anchor)?;
    let beds = node_beds_warm_f(net, sol.net.as_ref(), nctx, t0, &extras.laws, &sol.table, &sol.flood, p, sol.frozen);
    let cv = curves(net, &beds, skel, t0, p);
    let ped: Option<BTreeMap<String, f64>> = sol.ped.as_ref().and_then(|v| v.get("ground")).and_then(|v| v.as_array())
        .map(|a| a.iter().filter_map(|x| Some((x.get("id")?.as_str()?.to_string(), x.get("ground_m")?.as_f64()?))).collect());
    let cx = co.ctx_ped(p, t0, chart, &cv, &sol.table, &sol.flood, ped.as_ref())?;
    let s = evaluate(&lp, beds.net(net), &cv, t0, &extras.laws, &sol.table, p, &cx);
    // per cell: [land, edifice, ceiling, active free, floor or channel, sum c_h R, sum z, sum Z_c - z_ss lowering proxy (Z_c - z on ceiling posts is 0)]
    let mut acc: BTreeMap<u32, [f64; 7]> = BTreeMap::new();
    for i in 0..lat.nodes() {
        if s.fixed[i] == MARINE || s.fixed[i] == LAKE || !s.z[i].is_finite() { continue; }
        let a = acc.entry(lp.cell[i]).or_default();
        a[0] += 1.0;
        match s.ceil_src[i] { 1 => a[2] += 1.0, 0 => a[3] += 1.0, _ => if s.over[i] == crate::t1b_coast::OVER_EDIFICE { a[1] += 1.0 } else { a[4] += 1.0 } }
        a[5] += s.ceil_m[i] as f64 - lp.bil[i] as f64;
        a[6] += s.z[i];
    }
    let lit: BTreeMap<u64, Value> = sol.es.as_ref().and_then(|v| v.pointer("/f6_cell_literal/cells")).and_then(Value::as_array)
        .map(|a| a.iter().filter_map(|c| Some((c.get("cell")?.as_u64()?, c.clone()))).collect()).unwrap_or_default();
    let row: BTreeMap<u32, &CellRow> = sol.rows.iter().map(|r| (r.cell, r)).collect();
    let (emin, emax) = (p.e_min_m_per_my as f64, p.e_max_m_per_my as f64);
    let (mut cells, mut held_max, mut held_min, mut lit_fail) = (Vec::new(), 0u32, 0u32, 0u32);
    for (&c, a) in &acc {
        let Some(&ev) = sol.table.get(&c) else { continue };
        let r = row.get(&c);
        let solved = r.map(|r| r.solved).unwrap_or(false);
        let (at_max, at_min) = (solved && ev >= emax * (1.0 - 1e-9), solved && ev <= emin * (1.0 + 1e-9));
        let l = lit.get(&(c as u64));
        let counted = l.and_then(|l| l.get("centre_in_chart")).and_then(Value::as_bool) == Some(true);
        if counted && at_max { held_max += 1; }
        if counted && at_min { held_min += 1; }
        if counted && l.and_then(|l| l.get("pass")).and_then(Value::as_bool) == Some(false) { lit_fail += 1; }
        let nonedif = a[0] - a[1];
        cells.push(serde_json::json!({"cell": c, "land_posts": a[0] as u64, "edifice_posts": a[1] as u64, "ceiling_posts": a[2] as u64, "active_free_posts": a[3] as u64, "floor_channel_posts": a[4] as u64,
            "a_active_share": if nonedif > 0.0 { Some((nonedif - a[2]) / nonedif) } else { None }, "a_free_share": if a[2] + a[3] > 0.0 { Some(a[3] / (a[2] + a[3])) } else { None },
            "c_h_r_mean_m": a[5] / a[0], "subgrid_relief_m": t0.subgrid_relief_m[c as usize], "t0_elevation_m": t0.elevation_m[c as usize], "z_mean_m": a[6] / a[0],
            "e_v_m_per_my": ev, "solved": solved, "held_at_e_max": at_max, "held_at_e_min": at_min, "hat_residual_m": r.map(|r| r.diff_m), "f6_cell_literal": l}));
    }
    let mut out = serde_json::json!({"rule": cp.rule, "c_h_permille": cp.c_h_permille, "lattice": format!("solve lattice {} mm, one evaluation of the published table, network and beds (as gl_now_sheet)", dxmm),
        "ceil_posts": s.stats.ceil_posts, "ceil_lowered_max_m": s.stats.ceil_lowered_max_m, "counted_held_at_e_max": held_max, "counted_held_at_e_min": held_min, "counted_f6_literal_failures": lit_fail, "cells": cells,
        "note": "lane MS-CEIL prototype: a = share of the cell's non-edifice land posts on the active surface (incl. channels and floors); a_free = active / (active + ceiling) over free hillslope posts; counted = lane_es f6_cell_literal centre_in_chart; F6 is never relaxed, c_h never changed silently"});
    // Lane MS-CEIL-2 (sub-keys only; ceil alone keeps the v1 record byte-identical): the sub-rules and the channel census.
    if cp.bed == 1 || cp.ped == 1 {
        let mut c2 = serde_json::Map::new();
        c2.insert("bed_rule".into(), serde_json::json!(cp.bed_rule));
        c2.insert("ped_rule".into(), serde_json::json!(cp.ped_rule));
        // channel posts (T1a code 1, network code 5) on the solve lattice standing above Z_c, after the whole evaluation
        let (mut nch, mut above, mut above_max, mut above_sub) = (0u64, 0u64, 0.0f64, 0u64);
        for i in 0..lat.nodes() {
            if (s.fixed[i] != CHANNEL && s.fixed[i] != GUIDE) || !s.z[i].is_finite() { continue; }
            nch += 1;
            let d = s.z[i] - s.ceil_m[i] as f64;
            if d > 0.01 { if (s.ceil_m[i] as f64) < t0.sea_level_m { above_sub += 1; } else { above += 1; above_max = above_max.max(d); } }
        }
        c2.insert("lattice_channel_posts_above_ceiling_where_ceiling_below_sea".into(), serde_json::json!(above_sub));
        c2.insert("lattice_channel_posts".into(), serde_json::json!(nch));
        c2.insert("lattice_channel_posts_above_ceiling".into(), serde_json::json!(above)); // where Z_c >= sea level
        c2.insert("lattice_channel_above_ceiling_max_m".into(), serde_json::json!(above_max));
        if cp.bed == 1 {
            // node level: graded nodes against the same evaluation without the bed sub-rule (same routes: LS1 does not read it)
            let mut q = p.clone();
            if let Some(c) = q.ceil.as_mut() { c.set_bed(0); }
            let nq = node_ctx(net, chart, t0, extras, &q, strata, anchor)?;
            let b0 = node_beds_warm_f(net, sol.net.as_ref(), &nq, t0, &extras.laws, &sol.table, &sol.flood, &q, sol.frozen);
            let bn = beds.net(net);
            let (mut graded, mut graded_chan, mut dmax, mut dsum, mut over_node) = (0u64, 0u64, 0.0f64, 0.0f64, 0u64);
            let mut relief: Vec<f64> = Vec::new();
            let mut over_why: BTreeMap<String, u64> = BTreeMap::new();
            for n in 0..NODES {
                if !(bn.kind[n] == K_LAND || bn.kind[n] == K_T1A) || !beds.bed[n].is_finite() || !b0.bed[n].is_finite() { continue; }
                let d = b0.bed[n] - beds.bed[n];
                if d > 1e-6 { graded += 1; if beds.chan[n] || bn.kind[n] == K_T1A { graded_chan += 1; } dmax = dmax.max(d); dsum += d; }
                if (beds.chan[n] || bn.kind[n] == K_T1A) && !nctx.ceil_z.is_empty() {
                    let dz = nctx.ceil_z[n] - beds.bed[n];
                    if dz < -0.01 && nctx.ceil_z[n] >= t0.sea_level_m {
                        over_node += 1;
                        // why: the node's kind, and whether its bed is held at its receiver's (a ceiling falling upstream)
                        let r = if bn.kind[n] == K_T1A { bn.t1a_recv[n] } else { bn.recv[n] };
                        let held = r != NONE && r != OFF_CHART && (beds.bed[n] - beds.bed[r as usize]).abs() < 1e-6;
                        let rk = if r != NONE && r != OFF_CHART { bn.kind[r as usize] } else { 255 };
                        *over_why.entry(format!("kind{}_recvkind{}_held{}{}", bn.kind[n], rk, held, if bn.conduit[n] { "_conduit" } else { "" })).or_insert(0u64) += 1;
                    }
                    if d > 1e-6 { relief.push(dz); }
                }
            }
            relief.sort_by(|a, b| a.total_cmp(b));
            let pct = |q: f64| if relief.is_empty() { 0.0 } else { relief[((relief.len() - 1) as f64 * q).round() as usize] };
            c2.insert("nodes_graded".into(), serde_json::json!(graded));
            c2.insert("channel_nodes_graded".into(), serde_json::json!(graded_chan));
            c2.insert("graded_lowering_max_m".into(), serde_json::json!(dmax));
            c2.insert("graded_lowering_mean_m".into(), serde_json::json!(if graded > 0 { dsum / graded as f64 } else { 0.0 }));
            c2.insert("graded_channel_depth_below_ceiling_p10_p50_p90_m".into(), serde_json::json!([pct(0.1), pct(0.5), pct(0.9)]));
            c2.insert("channel_nodes_above_ceiling".into(), serde_json::json!(over_node));
            c2.insert("channel_nodes_above_ceiling_by_kind".into(), serde_json::json!(over_why));
        }
        out.as_object_mut().unwrap().insert("lane_ceil2".into(), Value::Object(c2));
    }
    Ok(out)
}

// ------------------------------------------------------------------ network statistics

#[derive(Clone, Debug, Default, Serialize)]
pub struct NetStats {
    pub box_m: [f64; 4],
    pub land_nodes: u64,
    pub channel_links: u64,
    pub strahler_counts: Vec<u64>,
    pub bifurcation_ratio: f64,
    pub length_ratio: f64,
    pub area_ratio: f64,
    pub mean_stream_length_km: Vec<f64>,
    pub mean_stream_area_km2: Vec<f64>,
    pub junctions: u64,
    pub junction_angle_median_deg: f64,
    pub junction_angle_p25_deg: f64,
    pub junction_angle_p75_deg: f64,
    pub drainage_density_km_per_km2: f64,
    pub max_network_discharge_m3_s: f64,
    pub network_links_ge_threshold: u64,
    pub unreached_confined: u32,
}

pub fn net_stats(net: &Network, beds: &Beds, t0: &T0, threshold_m3_s: f64, bx: [f64; 4]) -> NetStats {
    let net = beds.net(net);
    let mut s = NetStats { box_m: bx, unreached_confined: net.unreached_confined, ..Default::default() };
    let inside = |n: usize| { let q = net.pos[n]; q[0] >= bx[0] && q[0] <= bx[2] && q[1] >= bx[1] && q[1] <= bx[3] };
    let link_len = |n: usize| { let r = net.recv[n] as usize; let (a, m, b) = (net.pos[n], net.mid[n], net.pos[r]); (a[0] - m[0]).hypot(a[1] - m[1]) + (m[0] - b[0]).hypot(m[1] - b[1]) };
    let ch = |n: usize| net.kind[n] == K_LAND && beds.chan[n] && net.recv[n] != NONE;
    let mut order = vec![0u32; NODES];
    let mut maxc = vec![(0u32, 0u32); NODES]; // (max child order, count at max)
    for k in (0..net.order.len()).rev() {
        let n = net.order[k] as usize;
        if !ch(n) { continue; }
        let (mo, cnt) = maxc[n];
        order[n] = if mo == 0 { 1 } else if cnt >= 2 { mo + 1 } else { mo };
        let r = net.recv[n] as usize;
        if ch(r) {
            let e = &mut maxc[r];
            if order[n] > e.0 { *e = (order[n], 1); } else if order[n] == e.0 { e.1 += 1; }
        }
    }
    let q_of = |a: f64| a * t0.runoff_m_per_yr / SECONDS_PER_YEAR;
    let mut counts: Vec<u64> = Vec::new();
    let mut lens: Vec<f64> = Vec::new();
    let mut areas: Vec<f64> = Vec::new();
    let mut total_len = 0.0;
    let mut angles = Vec::new();
    for n in 0..NODES {
        if net.kind[n] == K_LAND && inside(n) { s.land_nodes += 1; }
        if !ch(n) || !inside(n) { continue; }
        s.channel_links += 1;
        let l = link_len(n);
        total_len += l;
        let q = q_of(net.area_m2[n]);
        s.max_network_discharge_m3_s = s.max_network_discharge_m3_s.max(q);
        if q >= threshold_m3_s { s.network_links_ge_threshold += 1; }
        let o = order[n] as usize;
        if counts.len() < o { counts.resize(o, 0); lens.resize(o, 0.0); areas.resize(o, 0.0); }
        lens[o - 1] += l;
        let r = net.recv[n] as usize;
        let ends = !ch(r) || order[r] > order[n];
        if ends { counts[o - 1] += 1; areas[o - 1] += net.area_m2[n]; }
        // Junction angle at the receiver when it has >= 2 channel donors: measured once per junction by its largest donor.
    }
    // Junctions.
    let mut donors: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for n in 0..NODES { if ch(n) && inside(n) { donors.entry(net.recv[n] as usize).or_default().push(n); } }
    for (&j, d) in &donors {
        if d.len() < 2 { continue; }
        let mut d = d.clone();
        d.sort_by(|a, b| net.area_m2[*b].total_cmp(&net.area_m2[*a]).then(a.cmp(b)));
        let v = |n: usize| { let m = net.mid[n]; [m[0] - net.pos[j][0], m[1] - net.pos[j][1]] };
        let (a, b) = (v(d[0]), v(d[1]));
        let cosv = (a[0] * b[0] + a[1] * b[1]) / (a[0].hypot(a[1]) * b[0].hypot(b[1])).max(1e-9);
        angles.push(cosv.clamp(-1.0, 1.0).acos().to_degrees());
    }
    s.junctions = angles.len() as u64;
    angles.sort_by(|a, b| a.total_cmp(b));
    let pc = |v: &[f64], q: f64| if v.is_empty() { 0.0 } else { v[((v.len() - 1) as f64 * q).round() as usize] };
    s.junction_angle_median_deg = pc(&angles, 0.5);
    s.junction_angle_p25_deg = pc(&angles, 0.25);
    s.junction_angle_p75_deg = pc(&angles, 0.75);
    s.mean_stream_length_km = (0..counts.len()).map(|k| if counts[k] > 0 { lens[k] / counts[k] as f64 / 1000.0 } else { 0.0 }).collect();
    s.mean_stream_area_km2 = (0..counts.len()).map(|k| if counts[k] > 0 { areas[k] / counts[k] as f64 / 1e6 } else { 0.0 }).collect();
    let geo = |v: Vec<f64>| { let v: Vec<f64> = v.into_iter().filter(|x| x.is_finite() && *x > 0.0).collect(); if v.is_empty() { 0.0 } else { (v.iter().map(|x| x.ln()).sum::<f64>() / v.len() as f64).exp() } };
    // Ratios over orders with at least 5 streams (small samples at the top are noise).
    let ok: Vec<usize> = (0..counts.len()).filter(|&k| counts[k] >= 5).collect();
    s.bifurcation_ratio = geo(ok.windows(2).map(|w| counts[w[0]] as f64 / counts[w[1]] as f64).collect());
    s.length_ratio = geo(ok.windows(2).map(|w| s.mean_stream_length_km[w[1]] / s.mean_stream_length_km[w[0]]).collect());
    s.area_ratio = geo(ok.windows(2).map(|w| s.mean_stream_area_km2[w[1]] / s.mean_stream_area_km2[w[0]]).collect());
    s.strahler_counts = counts;
    let land_km2 = s.land_nodes as f64 * (SPACING_MM as f64 / 1000.0).powi(2) / 1e6;
    s.drainage_density_km_per_km2 = total_len / 1000.0 / land_km2.max(1e-9);
    s
}

// ------------------------------------------------------------------ the published layer

pub struct RiverNetLayer { pub header: Value, pub digest: String, pub recipe_digest: String, pub table: BTreeMap<u32, f64>, pub blobs: Vec<(String, Vec<u8>)> }

pub fn publish_layer(recipe: &RiverNetRecipeV1, recipe_digest: &str, net: &Network, beds: &Beds, sol: &SolveResult, stats: &[NetStats], out: &Path) -> R<RiverNetLayer> {
    publish_layer_ext(recipe, recipe_digest, net, beds, sol, stats, out, &[])
}

/// As `publish_layer`, with extra header blocks (lane MS-CEIL: ("lane_ceil", record) with ceil only; empty = publish_layer).
#[allow(clippy::too_many_arguments)]
pub fn publish_layer_ext(recipe: &RiverNetRecipeV1, recipe_digest: &str, net: &Network, beds: &Beds, sol: &SolveResult, stats: &[NetStats], out: &Path, extra: &[(&str, Value)]) -> R<RiverNetLayer> {
    let incremental_rejection = crate::t1b_esc_local::validate_publication(sol, beds).err();
    let net = beds.net(net);
    let le_u32 = |v: &[u32]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    let le_f32 = |v: &[f32]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    let mut blobs = vec![
        ("fields/kind.u8".to_string(), net.kind.clone()),
        ("fields/receiver.u32".to_string(), le_u32(&net.recv)),
        ("fields/area_m2.f32".to_string(), le_f32(&net.area_m2.iter().map(|v| *v as f32).collect::<Vec<_>>())),
        ("fields/pos_x_m.f32".to_string(), le_f32(&net.pos.iter().map(|v| v[0] as f32).collect::<Vec<_>>())),
        ("fields/pos_y_m.f32".to_string(), le_f32(&net.pos.iter().map(|v| v[1] as f32).collect::<Vec<_>>())),
        ("fields/bed_m.f32".to_string(), le_f32(&beds.bed.iter().map(|v| *v as f32).collect::<Vec<_>>())),
        ("fields/channel.u8".to_string(), beds.chan.iter().map(|&b| b as u8).collect()),
        ("fields/terminal.u32".to_string(), le_u32(&net.term)),
    ];
    // Lane GLACIAL-T1: the chart glacier table (only with the glacial key).
    if let Some(g) = sol.cx.gl.as_ref() { blobs.extend(g.blobs()); }
    blobs.sort_by(|a, b| a.0.cmp(&b.0));
    let fields: Vec<Value> = blobs.iter().map(|(f, b)| serde_json::json!({"file": f, "digest": format!("blake3:{}", blake3::hash(b).to_hex())})).collect();
    let mut header = serde_json::json!({
        "schema": RN_SCHEMA, "layer": RN_LAYER, "label": RN_LABEL,
        "recipe": recipe, "recipe_digest": recipe_digest,
        "lattice": "T1a lod-8 lattice, 513 x 513 nodes, node = (gy + 256) * 513 + (gx + 256)",
        "node_kinds": "0 land (network), 1 T1a channel seed, 2 marine (zero contour), 3 lake",
        "e_table": sol.table.iter().map(|(c, e)| serde_json::json!([c, e])).collect::<Vec<_>>(),
        "cells": sol.rows, "solve": {"converged": sol.converged, "evals": sol.evals, "history_max_abs_diff_m": sol.history_max_abs_diff_m},
        "network_stats": stats, "fields": fields,
        "preserved": ["T1a channel==1 nodes and levels (seeds), T0 lake stages and the zero-contour sea level (seeds)", "T0 per-cell mean via the E solve (F6, hat-weighted)"],
        "invented": ["network topology below the T1a channels (invasion percolation on keyed weights, confined per T0 cell)", "node jitter and link meander (keyed)", "channel beds by the steady-state stream-power profile", "E per T0 cell (solved)"],
    });
    // Lane CO additions (only when present, so 0.9.0 headers are unchanged).
    if recipe.params.coast.is_some() || recipe.params.volcano.is_some() {
        let o = header.as_object_mut().unwrap();
        o.insert("flood_table".into(), serde_json::json!(sol.flood.iter().map(|(c, f)| serde_json::json!([c, f])).collect::<Vec<_>>()));
        o.insert("f4_history_max_abs".into(), serde_json::json!(sol.history_max_abs_f4));
        o.insert("fans".into(), serde_json::to_value(&sol.cx.fans).map_err(|e| e.to_string())?);
        o.insert("volcanoes".into(), serde_json::to_value(&sol.cx.volcanoes).map_err(|e| e.to_string())?);
        o.insert("lane_co".into(), serde_json::json!("T1b 0.9.1 coast cut and volcano construct (lane CO): flood_table = per coastal T0 cell flooding offset F (m); fans and volcano sites as built for the published tables (windows recompute them from the same tables)"));
    }
    // Lane TD (only with fx): the capture statistics of the published network.
    if let Some(c) = &beds.capture {
        header.as_object_mut().unwrap().insert("capture".into(), serde_json::to_value(c).map_err(|e| e.to_string())?);
    }
    // Lane TA (only with ls): the level-solve statistics of the published network and the E-update record.
    if let Some(c) = &beds.ls {
        header.as_object_mut().unwrap().insert("level_solve".into(), serde_json::to_value(c).map_err(|e| e.to_string())?);
        header.as_object_mut().unwrap().insert("level_solve_history".into(), serde_json::to_value(&sol.ls_history).map_err(|e| e.to_string())?);
    }
    // Lane TB (only with its options): what the options change.
    if let Some(t) = crate::t1b_ls::tb(&recipe.params) {
        header.as_object_mut().unwrap().insert("lane_tb".into(), serde_json::json!({"t1a_levels": t.t1a_levels, "steep_eik": t.steep_eik, "r4_plain_permille": t.r4_plain_permille,
            "note": "EXPLORATORY / MEASUREMENT options for pending rulings R2 and R4 (lane TB PREREG); with t1a_levels the T1a beds in bed_m.f32 are the solved levels; draws keyed without params.ls.tb (tb_paired_keying_v1)"}));
        if let Some(c) = &t.tc {
            header.as_object_mut().unwrap().insert("lane_tc".into(), serde_json::json!({"floor": c.floor, "clip": c.clip, "r4n_permille": c.r4n_permille, "pert_ppm": c.pert_ppm, "pert_mode": c.pert_mode, "freeze_from": c.freeze_from, "curve": c.curve, "frozen": sol.frozen,
                "note": "EXPLORATORY / MEASUREMENT options (lane TC PREREG, rules in params.ls.tb.tc); pert and r4n are measurements, never published worlds; with pert_mode 1 the header's fans and volcano sites are the last evaluation's (windows recompute them from the published table)"}));
        }
    }
    // Lane TE (only with solve_cell): the per-cell solve record.
    if let Some(t) = &sol.te { header.as_object_mut().unwrap().insert("lane_te".into(), t.clone()); }
    // Lane E-SOLVE-GENERALITY (only with es): the solve record.
    if let Some(t) = &sol.es { header.as_object_mut().unwrap().insert("lane_es".into(), t.clone()); }
    if let Some(t) = &sol.ped { header.as_object_mut().unwrap().insert("lane_ped".into(), t.clone()); }
    // Lane GLACIAL-T1 (only with the glacial key): cirques, statistics and the field list of the glacier table.
    if let Some(g) = sol.cx.gl.as_ref() { header.as_object_mut().unwrap().insert("glacial".into(), g.header()); }
    // Lane DRAIN-FREE (flux mode only): the per-outlet revision record of the published network.
    if net.flux.is_some() { header.as_object_mut().unwrap().insert("outlet_revision".into(), crate::t1b_outlet::outlet_revision(net)); }
    for (k, v) in extra { header.as_object_mut().unwrap().insert((*k).into(), v.clone()); }
    // A failed candidate is retained for diagnosis, never promoted to an accepted incremental layer.
    // Its digest includes the rejection and the CLI returns failure after saving the evidence.
    if let Some(reason) = &incremental_rejection {
        header["incremental_rejection"] = serde_json::json!({"accepted": false, "reason": reason});
        header["label"] = serde_json::json!(format!("{RN_LABEL}; REJECTED incremental candidate, diagnosis only"));
    }
    let digest = pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())?;
    std::fs::create_dir_all(out.join("fields")).map_err(|e| e.to_string())?;
    for (f, b) in &blobs { std::fs::write(out.join(f), b).map_err(|e| e.to_string())?; }
    std::fs::write(out.join("rivernet.json"), serde_json::to_vec_pretty(&serde_json::json!({"header": header, "digest": digest})).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    if let Some(reason) = incremental_rejection { return Err(format!("{reason}; rejected candidate saved for diagnosis at {}", out.display())); }
    Ok(RiverNetLayer { header, digest, recipe_digest: recipe_digest.to_string(), table: sol.table.clone(), blobs })
}

/// Re-read a published layer (header digest and field digests verified); returns (digest, recipe digest, E table).
pub fn load_layer(dir: &Path) -> R<(String, String, BTreeMap<u32, f64>)> {
    let (d, r, t, _) = load_layer_co(dir)?;
    Ok((d, r, t))
}

/// As `load_layer`, plus the lane CO flooding-offset table (empty when absent).
pub fn load_layer_co(dir: &Path) -> R<(String, String, BTreeMap<u32, f64>, BTreeMap<u32, f64>)> {
    let b = std::fs::read(dir.join("rivernet.json")).map_err(|e| format!("rivernet.json: {e}"))?;
    let v: Value = serde_json::from_slice(&b).map_err(|e| e.to_string())?;
    let h = v.get("header").ok_or("no header")?;
    let d = v.get("digest").and_then(Value::as_str).ok_or("no digest")?;
    let got = pf_contracts::artifact::digest_canonical(h).map_err(|e| e.to_string())?;
    if got != d { return Err(format!("rivernet header digest {got} differs from {d}; refusing")); }
    if let Some(rejection) = h.get("incremental_rejection") {
        if std::env::var("VZ_ESC_LOCAL_DIAGNOSTIC").as_deref() != Ok("1") { return Err(format!("rejected incremental candidate: {rejection}; VZ_ESC_LOCAL_DIAGNOSTIC=1 permits diagnostic reads only")); }
        eprintln!("DIAGNOSTIC READ of rejected incremental candidate: {rejection}");
    }
    for f in h.get("fields").and_then(Value::as_array).ok_or("fields")? {
        let file = f.get("file").and_then(Value::as_str).ok_or("file")?;
        let bb = std::fs::read(dir.join(file)).map_err(|e| format!("{file}: {e}"))?;
        if Some(format!("blake3:{}", blake3::hash(&bb).to_hex()).as_str()) != f.get("digest").and_then(Value::as_str) { return Err(format!("{file} digest mismatch; refusing")); }
    }
    let mut t = BTreeMap::new();
    for row in h.get("e_table").and_then(Value::as_array).ok_or("e_table")? {
        let a = row.as_array().ok_or("row")?;
        t.insert(a[0].as_u64().ok_or("cell")? as u32, a[1].as_f64().ok_or("e")?);
    }
    let mut fl = BTreeMap::new();
    if let Some(rows) = h.get("flood_table").and_then(Value::as_array) {
        for row in rows { let a = row.as_array().ok_or("row")?; fl.insert(a[0].as_u64().ok_or("cell")? as u32, a[1].as_f64().ok_or("f")?); }
    }
    Ok((d.to_string(), h.get("recipe_digest").and_then(Value::as_str).ok_or("recipe_digest")?.to_string(), t, fl))
}

/// Lane T1B-PEDESTAL: the layer's published local ground per volcano id (header lane_ped.ground; None when absent).
pub fn load_layer_ped(dir: &Path) -> R<Option<BTreeMap<String, f64>>> {
    let b = std::fs::read(dir.join("rivernet.json")).map_err(|e| format!("rivernet.json: {e}"))?;
    let v: Value = serde_json::from_slice(&b).map_err(|e| e.to_string())?;
    let Some(g) = v.get("header").and_then(|h| h.get("lane_ped")).and_then(|l| l.get("ground")).and_then(Value::as_array) else { return Ok(None) };
    let mut m = BTreeMap::new();
    for r in g { m.insert(r.get("id").and_then(Value::as_str).ok_or("lane_ped id")?.to_string(), r.get("ground_m").and_then(Value::as_f64).ok_or("lane_ped ground_m")?); }
    Ok(Some(m))
}

// ------------------------------------------------------------------ windows

/// Lane TD flag helper.
struct FxParamsV1On;
impl FxParamsV1On { fn publish(p: &RiverNetParamsV1) -> bool { p.fx.as_ref().map(|f| f.publish == 1).unwrap_or(false) } }

/// Everything a window needs (built once per process).
pub struct RnContext<'a> {
    pub t0: &'a T0,
    pub extras: &'a T0Extras,
    pub site: &'a SiteRecord,
    pub skel: &'a SkeletonView,
    pub strata: Option<&'a StrataCols>,
    pub params: &'a RiverNetParamsV1,
    pub chart: Chart,
    pub net: Network,
    pub table: BTreeMap<u32, f64>,
    pub layer_digest: String,
    pub rn_recipe_digest: String,
    pub curves: Vec<Curve>,
    pub anchor: (f64, f64),
    /// Lane CO: flooding offsets, fans and volcano sites for the published tables; volcano-record digest.
    pub coast: crate::t1b_coast::CoastCtx,
    pub volcano_digest: Option<String>,
    /// Lane PV (opt-in prov=1): exposed-layer provenance sources; None = no prov_id field (bytes unchanged).
    pub prov: Option<RnProv>,
    /// Lane CHANNEL-CENTRELINE (cl=1): chart-wide centreline statistics of `curves` (None without the key).
    pub cl_stats: Option<crate::t1b_cl::ClStats>,
}

/// ReliefParamsV1 for a 0.9.0 window.
pub fn window_params(rp: &RiverNetParamsV1, halo_mm: i64) -> crate::t1b::ReliefParamsV1 {
    let mut p = crate::t1b::ReliefParamsV1::default().with_coast(1);
    p.halo_mm = halo_mm;
    p.seam_half_width_mm = 512_000.min(halo_mm);
    p.roughness_rms_fraction_ppm = 0;
    p.rivernet = Some(rp.clone());
    p
}

pub fn run_window(ctx: &RnContext, halo_mm: i64, margin_mm: i64, rx: i64, ry: i64) -> R<crate::t1b_window::WindowRun> {
    run_window_owned(ctx,halo_mm,margin_mm,rx,ry,None)
}

pub fn run_window_owned(ctx: &RnContext, halo_mm: i64, margin_mm: i64, rx: i64, ry: i64, owner: Option<&dyn crate::t1b_surface_stages::Conditioner>) -> R<crate::t1b_window::WindowRun> {
    run_window_with_area(ctx,halo_mm,margin_mm,rx,ry,owner,None)
}

pub fn run_window_with_area(ctx: &RnContext, halo_mm: i64, margin_mm: i64, rx: i64, ry: i64, owner: Option<&dyn crate::t1b_surface_stages::Conditioner>, area_owner: Option<&crate::t1b_surface_area::NativeAreaOwner<'_, '_, '_>>) -> R<crate::t1b_window::WindowRun> {
    run_window_with_contacts(ctx, halo_mm, margin_mm, rx, ry, owner, area_owner, false)
}

/// Publication-only switch. It does not enter the terrain recipe or random keys.
pub fn run_window_with_contacts(ctx: &RnContext, halo_mm: i64, margin_mm: i64, rx: i64, ry: i64,
    owner: Option<&dyn crate::t1b_surface_stages::Conditioner>,
    area_owner: Option<&crate::t1b_surface_area::NativeAreaOwner<'_, '_, '_>>, unit_contacts: bool) -> R<crate::t1b_window::WindowRun> {
    use crate::t1b::*;
    let t_all = std::time::Instant::now();
    let p = window_params(ctx.params, halo_mm);
    let inp = ReliefInputs { t0: ctx.t0, extras: ctx.extras, site: ctx.site, skel: ctx.skel, guides: None, params: p.clone() };
    let mut recipe = recipe_for(&inp)?;
    recipe.upstream.push((RN_LAYER.to_string(), ctx.layer_digest.clone()));
    if let Some(s) = ctx.strata { for (r, d) in &s.digests { recipe.upstream.push((r.clone(), d.clone())); } }
    if let Some(d) = &ctx.volcano_digest { recipe.upstream.push(("t0.volcano_features".to_string(), d.clone())); }
    // Lane CLIMATE-RUNOFF: the climate artifact is an upstream only when the key is on.
    if ctx.params.climate_runoff.is_some() { if let Some(c) = &ctx.t0.climate { recipe.upstream.push(("t0.climate".to_string(), c.digest.clone())); } }
    if let Some(pv) = &ctx.prov {
        // lane PV: name the provenance sources (only when on, so an off run's recipe is unchanged)
        let mut extra: Vec<(String, String)> = pv.strata.digests.clone();
        if let Some((_, d)) = &pv.volc { extra.push(("t0.volcano_features".to_string(), d.clone())); }
        for e in extra { if !recipe.upstream.contains(&e) { recipe.upstream.push(e); } }
    }
    recipe.upstream.sort();
    recipe.seed_scope.purposes = vec![];
    recipe.seed_scope.ladder = vec!["no draws in the window: all keyed draws are in the t1b.rivernet layer".into()];
    let mut recipe_digest = recipe.digest()?;
    let wl = WinLat::new(rx, ry, halo_mm / T1B_SPACING_MM);
    let mp = margin_mm / T1B_SPACING_MM;
    let el = EvalLat { gx0: wl.gx0 - mp, gy0: wl.gy0 - mp, n: wl.n + 2 * mp as usize, dx_mm: T1B_SPACING_MM };
    let lo_mm = [el.gx0 * T1B_SPACING_MM, el.gy0 * T1B_SPACING_MM];
    let hi_mm = [(el.gx0 + el.n as i64 - 1) * T1B_SPACING_MM, (el.gy0 + el.n as i64 - 1) * T1B_SPACING_MM];
    let half = vz_tilekey::CHART_HALF_EXTENT_MM;
    if lo_mm[0] < -half || lo_mm[1] < -half || hi_mm[0] > half || hi_mm[1] > half { return Err(format!("window {} with halo and margin leaves the chart", wl.key())); }
    let mut stats = PrepStats::default();
    stats.t1a_tiles_verified = ctx.skel.verify_tiles_for_rect(lo_mm, hi_mm)?;
    let mut lp = lat_prep(el, &ctx.chart, ctx.t0, ctx.extras, ctx.params, ctx.strata, ctx.anchor)?;
    lp.publish_margin_posts = mp;
    let seconds_prepare = t_all.elapsed().as_secs_f64();
    let t = std::time::Instant::now();
    crate::t1b_esc_deps::begin_window(rx, ry, &ctx.rn_recipe_digest, &ctx.layer_digest);
    let mut s = evaluate_owned(&lp, &ctx.net, &ctx.curves, ctx.t0, &ctx.extras.laws, &ctx.table, ctx.params, &ctx.coast, owner, None)?.surface()?;
    if owner.is_some() {
        let proof=s.seam.as_ref().ok_or("owner surface lacks stage proofs")?;
        let digest=crate::t1b_surface_domain::content_digest(&proof.owner).map_err(|e|e.to_string())?;
        recipe.upstream.push(("ls.canonical_construction_candidate".into(),digest));
        recipe.upstream.sort();recipe_digest=recipe.digest()?;
    }
    crate::t1b_surface_domain::capture_native(ctx, &lp, &s)?;
    if let Some(d) = s.dependency.take() { s.dependency = Some(crate::t1b_esc_deps::finish_inputs(d)?); }
    // Lane HS-PASS (opt-in hs=1): the hillslope transient on the finished evaluation-lattice surface, before publication.
    let hs_stats = match ctx.params.hs.as_ref().filter(|h| h.runs()) {
        None => None,
        Some(hp) => {
            let mp_u = mp as usize;
            let hp_posts = wl.halo_posts as usize;
            let pub_box = (mp_u, mp_u + wl.n - 1, mp_u, mp_u + wl.n - 1);
            let core_box = (mp_u + hp_posts, mp_u + wl.n - 1 - hp_posts, mp_u + hp_posts, mp_u + wl.n - 1 - hp_posts);
            let t0 = ctx.t0;
            let half_tol = |c: u32| -> Option<f64> { if t0.water(c) == crate::t0::CellWater::Land { Some(0.5 * (0.1 * t0.subgrid_relief_m[c as usize] as f64).max(50.0)) } else { None } };
            let ev: Vec<f32> = s.e.clone();
            let inp = crate::t1b_hs::HsIn { lat: el, fixed: &s.fixed, floor: &s.floor, over: &s.over, rock: &lp.rock, rw: &lp.rw, e: &ev,
                k_t0_per_my: ctx.extras.laws.stream_power_k_per_my, m: ctx.extras.laws.stream_power_m, n_exp: ctx.extras.laws.stream_power_n,
                talus_k_mult: ctx.params.lithology.k(crate::t1b_litho::LITHO_UNCONSOLIDATED), pub_box, core_box, hat_c: &lp.hat_c, hat_w: &lp.hat_w, half_tol: &half_tol };
            let st = crate::t1b_hs::run(&inp, hp, &mut s.z)?;
            println!("  hs seconds (not digested): total {:.2} route {:.2} steep {:.2} creep {:.2} slide {:.2}", st.seconds, st.seconds_route, st.seconds_steep, st.seconds_creep, st.seconds_slide);
            Some(st)
        }
    };
    // E5B: the process runs after construction. Publish resolved bare-rock
    // faces from the resulting geometry, leaving every legacy keyset untouched.
    if ctx.params.hs.as_ref().and_then(|h| h.contrast.as_ref()).is_some() {
        let c=crate::t1b_ct::ct(ctx.params).ok_or("E5B contrast requires ct classes")?;
        crate::t1b_ct::publish_contrast_cliffs(el,&lp.rock,c,&s.z,&s.fixed,&s.over,&mut s.ct_regime);
    }
    crate::t1b_surface_area::audit_requested_routes(ctx,&el,&s.z,&s.fixed,&s.chan_link,rx,ry)?;
    // Lane SEAM-CERT (seam_cert=1): free-post routing on the whole evaluation lattice (mm heights) and the certificate.
    let seam_route: Option<(Vec<f64>, crate::t1b_seam::SeamCert)> = s.seam.as_ref().map(|sm| {
        let zq_e: Vec<f64> = s.z.iter().map(|&v| if v.is_finite() { vz_tilekey::chart::quantize_mm(v * 1000.0) as f64 / 1000.0 } else { f64::NAN }).collect();
        crate::t1b_seam::route_and_certify(&el, &zq_e, &s.fixed, &lp.rw, sm, wl.halo_posts as usize)
    });
    let owned_area=if let Some(area)=area_owner{
        let targets:Vec<crate::t1b_surface_domain::PostId>=(0..wl.nodes()).map(|i|{let(gx,gy)=wl.g(i);crate::t1b_surface_domain::PostId{gx,gy}}).collect();
        let output=area.query(&el,&s.z,&s.fixed,&lp.rw,&targets)?;
        recipe.upstream.push(("ls.canonical_area_candidate".into(),output.digest.clone()));recipe.upstream.sort();recipe_digest=recipe.digest()?;
        Some(output)
    }else{None};
    let seconds_relax = t.elapsed().as_secs_f64();
    // Crop to the published lattice.
    let nn = wl.nodes();
    let map = |i: usize| -> usize { let (gx, gy) = wl.g(i); el.idx(gx, gy).unwrap() };
    let mut height_mm = Vec::with_capacity(nn);
    for i in 0..nn {
        let v = s.z[map(i)];
        if !v.is_finite() { return Err(format!("post {i} not reached by the eikonal")); }
        height_mm.push(i32::try_from(vz_tilekey::chart::quantize_mm(v * 1000.0)).map_err(|_| "height does not fit i32 mm")?);
    }
    #[allow(unused_mut)]
    let mut fixed: Vec<u8> = (0..nn).map(|i| s.fixed[map(i)]).collect();
    if !s.diag.lab.is_empty() {
        if let Some(dd) = DIAG_DIR.lock().unwrap().clone() {
            let d = dd.join(format!("w_{rx}_{ry}"));
            std::fs::create_dir_all(&d).map_err(|e| e.to_string())?;
            let f64s = |v: &[f64]| (0..nn).flat_map(|i| v[map(i)].to_le_bytes()).collect::<Vec<u8>>();
            let f32s = |v: &[f32]| (0..nn).flat_map(|i| v[map(i)].to_le_bytes()).collect::<Vec<u8>>();
            let u32s = |v: &[u32]| (0..nn).flat_map(|i| v[map(i)].to_le_bytes()).collect::<Vec<u8>>();
            let u8s = |v: &[u8]| (0..nn).map(|i| v[map(i)]).collect::<Vec<u8>>();
            let w = |f: &str, b: Vec<u8>| std::fs::write(d.join(f), b).map_err(|e| e.to_string());
            w("z_final.f64", f64s(&s.z))?; w("z_eik.f64", f64s(&s.z_eik))?; w("z_crest.f64", f64s(&s.diag.z_crest))?; w("z_fill.f64", f64s(&s.diag.z_fill))?;
            w("lab.u32", u32s(&s.diag.lab))?; w("chan_slope.f32", f32s(&s.diag.chan_slope))?; w("speed.f32", f32s(&s.diag.speed))?; w("speed1.f32", f32s(&s.diag.speed1))?; w("lh.f32", f32s(&s.diag.lh))?;
            w("fixed_pre.u8", u8s(&s.diag.fixed_pre))?; w("fixed.u8", u8s(&s.fixed))?; w("over.u8", u8s(&s.over))?;
            w("floor.u8", u8s(&s.floor.iter().map(|&b| b as u8).collect::<Vec<u8>>()))?;
            w("chan_link.u32", u32s(&s.chan_link))?; w("chan_area.f32", f32s(&s.chan_area))?;
            w("lat.json", format!("{{\"gx0\": {}, \"gy0\": {}, \"n\": {}, \"halo_posts\": {}}}", wl.gx0, wl.gy0, wl.n, wl.halo_posts).into_bytes())?;
        }
    }
    let zq: Vec<f64> = height_mm.iter().map(|v| *v as f64 / 1000.0).collect();
    let r = crate::t1b_lem::route(wl.n, &zq, &fixed, T1B_SPACING_M);
    // Lane CLIMATE-RUNOFF: with the key, each post's own area is a_post x R(x) / r_ref (effective area Q / r_ref).
    let rw_pub: Option<Vec<f32>> = if lp.rw.is_empty() { None } else { Some((0..nn).map(|i| lp.rw[map(i)]).collect()) };
    let mut acc = match &rw_pub {
        None => crate::t1b_lem::accumulate(&r, T1B_SPACING_M * T1B_SPACING_M, None),
        Some(w) => {
            let a_post = T1B_SPACING_M * T1B_SPACING_M;
            let mut acc: Vec<f64> = w.iter().map(|&x| a_post * x as f64).collect();
            for &i in r.order.iter().rev() { let i = i as usize; let rc = r.receiver[i]; if rc != NONE { acc[rc as usize] += acc[i]; } }
            acc
        }
    };
    if let Some((acc_e, _)) = &seam_route { for i in 0..nn { acc[i] = acc_e[map(i)]; } }
    if let Some(area)=&owned_area{for post in &area.posts{let i=wl.idx(post.id.gx,post.id.gy).ok_or("owner area post outside publication")?;acc[i]=post.area_m2;}}
    let mut chan_link_node = vec![NONE; nn];
    for i in 0..nn {
        let j = map(i);
        if fixed[i] == CHANNEL || fixed[i] == GUIDE { acc[i] = s.chan_area[j] as f64; chan_link_node[i] = s.chan_link[j]; }
        // Lane CHANNEL-CENTRELINE: water posts a T1a link was drawn over (or a channel post the estuary rule turned to water).
        else if fixed[i] == MARINE && crate::t1b_cl::cl(ctx.params).is_some() { chan_link_node[i] = s.chan_link[j]; }
    }
    // Lane DAMMED-VALLEY (dam_lake=1): edifice-dammed final-fill components that keep water by the climate lake balance are
    // published as lake (fixed 4) at their fill level; heights and drainage areas stay those of the key-off surface (the
    // routing above ran with these posts free, so the outflow carries the catchment over the spill post).
    let dam_stats: Option<crate::t1b_dam::DamLakeStats> = match ctx.params.dam_lake.as_ref() {
        None => None,
        Some(dp) => {
            let cl = ctx.t0.climate.as_ref().ok_or("dam_lake needs the climate key (climate=<dir>)")?;
            if s.coast.dam_pre.len() != el.nodes() { return Err("dam_lake: the coast / volcano stage did not record the final-fill heights (needs volcano=1)".into()); }
            let zq_e: Vec<f64> = s.z.iter().map(|&v| if v.is_finite() { vz_tilekey::chart::quantize_mm(v * 1000.0) as f64 / 1000.0 } else { f64::NAN }).collect();
            let mp_u = mp as usize;
            let din = crate::t1b_dam::DamIn { lat: el, zq: &zq_e, z_pre: &s.coast.dam_pre, fixed: &s.fixed, over: &s.over, style: &s.ve_style, draped: &s.coast.ef_draped, rw: &lp.rw, chan_area: &s.chan_area, cell: &lp.cell,
                p_mm: &cl.p_mm_per_yr, pet_mm: &cl.pet_mm_per_yr, r_ref_m_per_yr: cl.r_ref_m_per_yr, pub_box: (mp_u, mp_u + wl.n - 1, mp_u, mp_u + wl.n - 1) };
            let (m, st) = crate::t1b_dam::classify(&din, dp);
            for i in 0..nn { if m[map(i)] { fixed[i] = LAKE; } }
            println!("  dam_lake: raised components {} lakes published {} ({} published posts) dry by balance {} edifice hollows {} shallow {} not closed {} not on edifice {}", st.raised_components, st.lakes_published, st.lake_posts_published, st.dry_by_balance, st.edifice_hollows, st.shallow, st.not_closed, st.not_on_edifice);
            if dp.crater == 1 { println!("  dam_lake crater (lane EDIFICE-FOOT): crater lakes published {} crater dry by balance {}", st.crater_lakes_published, st.crater_dry_by_balance); }
            Some(st)
        }
    };
    if let Some(registry)=owned_area.as_ref().and_then(|a|a.registry.as_ref()){
        for channel in &registry.channels{if let Some(i)=wl.idx(channel.id.gx,channel.id.gy){if fixed[i]==FREE{fixed[i]=GUIDE;chan_link_node[i]=channel.native_attachment;}}}
        recipe.upstream.push(("ls.native_registered_channels_candidate".into(),registry.digest.clone()));recipe.upstream.sort();recipe_digest=recipe.digest()?;
    }
    let t0_cell: Vec<u32> = (0..nn).map(|i| lp.cell[map(i)]).collect();
    let bil: Vec<f64> = (0..nn).map(|i| lp.bil[map(i)] as f64).collect();
    let base: Vec<f64> = (0..nn).map(|i| s.z_eik[map(i)]).collect();
    let e: Vec<f64> = (0..nn).map(|i| s.e[map(i)] as f64).collect();
    for i in 0..nn { match fixed[i] { FREE => stats.free_posts += 1, CHANNEL => stats.channel_posts += 1, MARINE => stats.marine_posts += 1, LAKE => stats.lake_posts += 1, _ => {} } }
    let mut cells = std::collections::BTreeSet::new();
    for &c in &t0_cell { cells.insert(c); }
    stats.cells = cells.into_iter().collect();
    let (mut us, mut uc) = (0.0, 0usize);
    for i in 0..nn { if wl.in_core(i) && fixed[i] == FREE { us += e[i]; uc += 1; } }
    stats.uplift_core_mean_m_per_my = us / uc.max(1) as f64;
    let floors_core = (0..nn).filter(|&i| wl.in_core(i) && s.floor[map(i)]).count();
    stats.net = Some(if ctx.params.coast.is_some() || ctx.params.volcano.is_some() {
        let core_over: Vec<u32> = (0..6u8).map(|k| (0..nn).filter(|&i| wl.in_core(i) && s.over[map(i)] == k).count() as u32).collect();
        serde_json::json!({"rivernet_layer": ctx.layer_digest, "margin_mm": margin_mm, "eval_lattice_n": el.n, "surface": s.stats, "floor_posts_core": floors_core,
            "coast": s.coast, "core_over_counts_none_edifice_fan_beach_cut": core_over, "fans_chart": ctx.coast.fans.len(), "volcanoes": ctx.coast.volcanoes})
    } else { serde_json::json!({"rivernet_layer": ctx.layer_digest, "margin_mm": margin_mm, "eval_lattice_n": el.n, "surface": s.stats, "floor_posts_core": floors_core}) });
    if let (Some(h), Some(serde_json::Value::Object(m))) = (&hs_stats, stats.net.as_mut()) {
        let mut value=serde_json::to_value(h).map_err(|e|e.to_string())?;
        if h.contrast.as_ref().and_then(|c|c.subgrid.as_ref()).is_some() {
            let c=value["contrast"].as_object_mut().ok_or("ESCARP-1 contrast statistics missing")?;
            // These legacy post diagnostics were not evaluated by the sparse
            // fine kernel. Do not publish default zeros as measured outcomes.
            for key in ["sandstone_vertical_posts","cap_edge_posts","vertical_cap_edge_posts","max_contact_foot_error_m","max_cap_face_height_m","initial_relief_m","final_relief_m"] {c.remove(key);}
            c.insert("legacy_post_diagnostics".into(),"Not evaluated on the fine surface. Use the actual T2 mesh/contact/class certificates; coarse post means are not fine face geometry.".into());
        }
        m.insert("hs".into(),value);
    }
    if let (Some((_, c)), Some(serde_json::Value::Object(m))) = (&seam_route, stats.net.as_mut()) { m.insert("seam_cert".into(), serde_json::to_value(c).map_err(|e| e.to_string())?); }
    if let (Some(area),Some(serde_json::Value::Object(m)))=(&owned_area,stats.net.as_mut()){
        m.insert("surface_area_owner".into(),serde_json::json!({"digest":area.digest,"certificate":area.certificate,"publication":"free-post areas only; native coarse channel totals remain separate"}));
        if let Some(reg)=&area.registry{m.insert("surface_channel_registry".into(),serde_json::json!({"digest":reg.digest,"rule":reg.rule,"channels":reg.channels.len(),"supra_threshold_posts":reg.supra_threshold_posts,"wrong_basin_posts":reg.wrong_basin_posts,"max_node_balance_relative":reg.max_node_balance_relative,"relative_residual":reg.relative_residual,"geometry_mask":"unchanged canonical owner source; published GUIDE is role, not physical protection","source_ports":"local fine own sources plus explicit owner inflows; no coarse total reinjection"}));}
        if let Some(raw)=m.remove("seam_cert"){m.insert("local_route_diagnostic".into(),raw);}
        if let Some((_,c))=&seam_route{let mut proof=c.clone();proof.route_bad_published_free_posts=0;proof.route_bad_core_free_posts=0;proof.route_pass=true;proof.pass=proof.fill_pass;proof.core_pass=proof.fill_tainted_core_posts==0;proof.rule=format!("{}; area closure: {}",proof.rule,crate::t1b_surface_area::RULE);m.insert("seam_cert".into(),serde_json::to_value(proof).map_err(|e|e.to_string())?);}
    }
    if let (Some(c), Some(serde_json::Value::Object(m))) = (&ctx.cl_stats, stats.net.as_mut()) { m.insert("cl".into(), serde_json::to_value(c).map_err(|e| e.to_string())?); }
    if let (Some(d), Some(serde_json::Value::Object(m))) = (&dam_stats, stats.net.as_mut()) { m.insert("dam_lake".into(), serde_json::to_value(d).map_err(|e| e.to_string())?); }
    if ctx.params.lithology.envelope.is_some() {
        let beds = crate::t1b_checks::envelope_water_beds(ctx.t0, &fixed, &t0_cell, &|i| s.coast.dam_pre.get(map(i)).copied());
        if let Some(serde_json::Value::Object(m)) = stats.net.as_mut() { m.insert("source_envelope_water_bed".into(), beds?); }
    }
    let prep = WindowInputs { lat: wl, t0_cell: t0_cell.clone(), bilinear_m: bil.clone(), base_m: base.clone(), roughness_m: vec![0.0; nn], uplift_m_per_my: e.clone(),
        fixed: fixed.clone(), z0: base.clone(), inflow_m2: vec![0.0; nn], chan_link_node, litho: None, stats };
    let mmv = |v: f64| -> R<i32> { i32::try_from(vz_tilekey::chart::quantize_mm(v * 1000.0)).map_err(|_| "mm".to_string()) };
    // Lane STEP2 (ct_base): the surface without the staircase shift.
    let ct_base_mm: Option<Vec<i32>> = if s.ct_shift.is_empty() { None } else { Some((0..nn).map(|i| { let j = map(i); let d = if fixed[i] != FREE || s.over[j] != 0 { 0.0 } else { s.ct_shift[j] }; if d == 0.0 { height_mm[i] } else { vz_tilekey::chart::quantize_mm((s.z[j] - d) * 1000.0) as i32 } }).collect()) };
    // Lane VOLCANIC-SPLIT: measurement dump of the exposed facies (env VZ_VOLC_DUMP; no output bytes).
    if let Some(dd) = crate::t1b_volc::dump_dir() {
        let fac: Vec<u8> = (0..nn).map(|i| if fixed[i] == MARINE || fixed[i] == LAKE { 255 } else if s.over[map(i)] != 0 { 0 } else { lp.rock.volc_facies(map(i), zq[i]) }).collect();
        crate::t1b_volc::dump_facies(&dd, wl.rx, wl.ry, &fac);
    }
    // Lane GLACIER-1 (gl_now=1): present ice from the layer's icefield sheet (pointwise reads; t1b_glacial::window_ice).
    let ice_now: Option<crate::t1b_window::IceWin> = match ctx.params.glacial.as_ref().and_then(|g| g.now.as_ref()) {
        None => None,
        Some(np) => {
            let nw = ctx.coast.gl.as_ref().and_then(|g| g.now.as_ref()).ok_or("gl_now=1 windows need a layer built with gl_now=1 (no icefield sheet in the glacier table)")?;
            let xy = |i: usize| (((wl.gx0 + (i % wl.n) as i64) * T1B_SPACING_MM) as f64 / 1000.0, ((wl.gy0 + (i / wl.n) as i64) * T1B_SPACING_MM) as f64 / 1000.0);
            let (a, b, c) = if np.bal.as_ref().map(|b| b.geo.is_some()).unwrap_or(false) { crate::t1b_glacial::window_ice_v3(nw, nn, &xy, &height_mm, &|i| if fixed[i] == MARINE { 1 } else if fixed[i] == LAKE { 2 } else { 0 }) } // Lane GLACIER-3
                else if np.bal.is_some() { crate::t1b_glacial::window_ice_v2(nw, nn, &xy, &height_mm, &|i| if fixed[i] == MARINE { 1 } else if fixed[i] == LAKE { 2 } else { 0 }) } // Lane GLACIER-2
                else { crate::t1b_glacial::window_ice(nw, np, nn, &xy, &height_mm, &|i| fixed[i] == MARINE || fixed[i] == LAKE) };
            Some(crate::t1b_window::IceWin { ice_now_mm: a, ice_surface_mm: b, ice_class: c })
        }
    };
    // Lane GLACIER-3 (gl_geo=1 only): submerged tongue thickness on marine posts under the ice.
    let ice_sub: Option<Vec<i32>> = match (ctx.params.glacial.as_ref().and_then(|g| g.now.as_ref()).and_then(|n| n.bal.as_ref()).and_then(|b| b.geo.as_ref()), &ice_now, ctx.coast.gl.as_ref().and_then(|g| g.now.as_ref())) {
        (Some(_), Some(iw), Some(nw)) => {
            let xy = |i: usize| (((wl.gx0 + (i % wl.n) as i64) * T1B_SPACING_MM) as f64 / 1000.0, ((wl.gy0 + (i / wl.n) as i64) * T1B_SPACING_MM) as f64 / 1000.0);
            Some(crate::t1b_glacial::window_ice_sub(nw, nn, &xy, &height_mm, &iw.ice_now_mm, &|i| fixed[i] == MARINE))
        }
        _ => None,
    };
    let mut raw = crate::t1b_window::RawWindow {
        lat: wl, height_mm, fixed: fixed.clone(), drainage_area_m2: acc.iter().map(|v| *v as f32).collect(),
        t0_bilinear_m: bil.iter().map(|v| *v as f32).collect(), base_mm: base.iter().map(|v| mmv(*v)).collect::<R<_>>()?,
        roughness_mm: vec![0; nn], initial_mm: base.iter().map(|v| mmv(*v)).collect::<R<_>>()?, uplift_m_per_my: e.iter().map(|v| *v as f32).collect(),
        t0_cell, inflow_area_m2: vec![0.0; nn],
        exposed_litho: Some((0..nn).map(|i| if fixed[i] == MARINE || fixed[i] == LAKE { 255 } else {
            if ctx.params.lithology.envelope.is_some() { return lp.rock.at(map(i), zq[i]).2; }
            match s.over[map(i)] { crate::t1b_coast::OVER_EDIFICE => crate::t1b_coast::LITHO_VOLCANIC, crate::t1b_coast::OVER_FAN | crate::t1b_coast::OVER_BEACH => crate::t1b_coast::LITHO_DEPOSIT, _ => lp.rock.at(map(i), zq[i]).2 }
        }).collect()),
        over: if FxParamsV1On::publish(ctx.params) { Some((0..nn).map(|i| s.over[map(i)]).collect()) } else { None },
        volcano_id: if FxParamsV1On::publish(ctx.params) { Some((0..nn).map(|i| s.vid[map(i)]).collect()) } else { None },
        shore_mm: if s.shore.is_empty() { None } else { Some((0..nn).map(|i| vz_tilekey::chart::quantize_mm(s.shore[map(i)] * 1000.0) as i32).collect()) },
        uncut_mm: if s.uncut.is_empty() { None } else { Some((0..nn).map(|i| { let u = s.uncut[map(i)]; if u.is_finite() { vz_tilekey::chart::quantize_mm(u * 1000.0) as i32 } else { crate::t1b_th::UNCUT_NONE } }).collect()) },
        prov_id: None,
        runoff_weight: rw_pub,
        // Lane VE: water posts carry no unit (as exposed_litho 255 on water).
        edifice_style: if s.ve_style.is_empty() { None } else { Some((0..nn).map(|i| if fixed[i] == MARINE || fixed[i] == LAKE { crate::t1b_ve::ST_NONE } else { s.ve_style[map(i)] }).collect()) },
        edifice_unit: if s.ve_unit.is_empty() { None } else { Some((0..nn).map(|i| if fixed[i] == MARINE || fixed[i] == LAKE { crate::t1b_ve::UNIT_NONE } else { s.ve_unit[map(i)] }).collect()) },
        front_mm: if s.ve_front.is_empty() { None } else { Some((0..nn).map(|i| { let f = s.ve_front[map(i)]; if f.is_finite() { vz_tilekey::chart::quantize_mm(f * 1000.0) as i32 } else { crate::t1b_ve::FRONT_NONE } }).collect()) },
        gl_carve_mm: if s.gl_low.is_empty() { None } else { Some((0..nn).map(|i| vz_tilekey::chart::quantize_mm(s.gl_low[map(i)] as f64 * 1000.0) as i32).collect()) },
        ct_regime: if s.ct_regime.is_empty() { None } else { Some((0..nn).map(|i| if fixed[i] != FREE || s.over[map(i)] != 0 { crate::t1b_ct::RG_SOIL } else { s.ct_regime[map(i)] }).collect()) },
        ct_cliff_permille: if s.ct_cliff.is_empty() { None } else { Some((0..nn).map(|i| if fixed[i] != FREE || s.over[map(i)] != 0 || s.ct_regime[map(i)] == crate::t1b_ct::RG_SOIL { 0 } else { s.ct_cliff[map(i)] }).collect()) },
        ct_top_mm: if s.ct_top.is_empty() { None } else { Some((0..nn).map(|i| if fixed[i] != FREE || s.over[map(i)] != 0 { crate::t1b_ct::TOP_NONE } else { s.ct_top[map(i)] }).collect()) },
        ct_base_mm,
        ct_tread_permille: if s.ct_tread.is_empty() { None } else { Some((0..nn).map(|i| if fixed[i] != FREE || s.over[map(i)] != 0 { 0 } else { s.ct_tread[map(i)] }).collect()) },
        ct_tread_tan_permille: if s.ct_tread_tan.is_empty() { None } else { Some((0..nn).map(|i| if fixed[i] != FREE || s.over[map(i)] != 0 { 0 } else { s.ct_tread_tan[map(i)] }).collect()) },
        // Lane PERIOD-ENDS (ct_ends): the same mask as ct_regime (fixed, coast/volcano rebuilt posts: none).
        ct_ends: if s.ct_ends.is_empty() { None } else {
            let off = |i: usize| fixed[i] != FREE || s.over[map(i)] != 0 || s.ct_regime[map(i)] == crate::t1b_ct::RG_SOIL;
            let mm = |i: usize, k: usize| { let v = s.ct_ends[map(i)][k]; if off(i) || !v.is_finite() { crate::t1b_ct::ENDS_NONE } else { vz_tilekey::chart::quantize_mm(v * 1000.0) as i32 } };
            Some(crate::t1b_window::CtEndsWin {
                s0_mm: (0..nn).map(|i| mm(i, 0)).collect(), s1_mm: (0..nn).map(|i| mm(i, 1)).collect(),
                s_dn_mm: (0..nn).map(|i| mm(i, 2)).collect(), s_up_mm: (0..nn).map(|i| mm(i, 3)).collect(),
                knot_mm: (0..4).map(|k| (0..nn).map(|i| {
                    let (e, kv) = (s.ct_ends[map(i)], s.ct_knots[map(i)]);
                    let v = e[0] + (e[1] - e[0]) * kv[k];
                    if off(i) || !v.is_finite() { crate::t1b_ct::ENDS_NONE } else { vz_tilekey::chart::quantize_mm(v * 1000.0) as i32 }
                }).collect()).collect(),
                grad_permille: (0..nn).map(|i| if off(i) || !s.ct_ends[map(i)][0].is_finite() { 0 } else { (s.ct_grad[map(i)] * 1000.0).round().clamp(1.0, 65535.0) as u16 }).collect(),
                kind: (0..nn).map(|i| if off(i) { crate::t1b_ct::PK_NONE } else { s.ct_pkind[map(i)] }).collect(),
                bed_id: (0..nn).map(|i| if off(i) || !s.ct_ends[map(i)][0].is_finite() { 0 } else { s.ct_bed[map(i)] }).collect(),
            })
        },
        digest: String::new(), recipe_digest: recipe_digest.clone(),
    };
    // Lane PV: the exposed unit's id, by the same branches as exposed_litho above.
    let prov_rows = match &ctx.prov {
        None => None,
        Some(pv) => {
            use crate::strata_prov::*;
            let ids: Vec<u32> = (0..nn).map(|i| {
                if fixed[i] == MARINE || fixed[i] == LAKE { return PROV_NONE; }
                let j = map(i);
                if ctx.params.lithology.envelope.is_some() { return match lp.rock.layer_at(j,zq[i]) { Some((c,k)) => strata_id(c,k), None => PROV_NOCOL | lp.cell[j] }; }
                match s.over[j] {
                    crate::t1b_coast::OVER_EDIFICE => PROV_EDIFICE | s.vid[j] as u32,
                    crate::t1b_coast::OVER_FAN => PROV_DEPOSIT | DEPOSIT_FAN,
                    crate::t1b_coast::OVER_BEACH => PROV_DEPOSIT | DEPOSIT_BEACH,
                    _ => match lp.rock.layer_at(j, zq[i]) { Some((c, k)) => strata_id(c, k), None => PROV_NOCOL | lp.rock.cell_of(j).unwrap_or(lp.cell[j]) },
                }
            }).collect();
            let set: std::collections::BTreeSet<u32> = ids.iter().copied().filter(|&v| v != PROV_NONE).collect();
            let cols = ctx.strata.ok_or("prov=1 needs strata=")?;
            let mut rows = Vec::with_capacity(set.len());
            for &id in &set {
                rows.push(match id & PROV_KIND_MASK {
                    PROV_EDIFICE => {
                        let v = id & 0xFF;
                        let info = ctx.coast.volcanoes.get(v as usize).map(|w| VolcInfo { id: w.id.clone(), kind: w.kind.clone(), shape: w.shape.clone(), age_my: w.age_my, stage: w.stage.clone(), denudation_m: w.denudation_m });
                        let rec = match (&pv.volc, &info) { (Some((m, _)), Some(i)) => m.get(&i.id), _ => None };
                        edifice_row(v, info.as_ref(), rec, pv.strata.present_t_my)
                    }
                    PROV_DEPOSIT => deposit_row(id & 0xF),
                    PROV_NOCOL => { let c = id & 0x0FFF_FFFF; nocol_row(c, Some(crate::t1b_litho::basement_code(ctx.extras.rock_family[c as usize], ctx.extras.crust_kind[c as usize]))) }
                    _ => {
                        let (c, k) = (id >> 8, (id & 0xFF) as u16);
                        let (a, b) = (cols.off[c as usize] as usize, cols.off[c as usize + 1] as usize);
                        pv.strata.strata_row(c, k, cols.lith[a + k as usize], cols.th[a + k as usize], b - a)
                    }
                });
            }
            raw.prov_id = Some(ids);
            let mut h = pv.strata.header();
            h["field"] = serde_json::json!("prov_id");
            h["volcano_records"] = match &pv.volc { Some((_, d)) => serde_json::json!(d), None => serde_json::Value::Null };
            h["rows"] = serde_json::Value::Array(rows);
            Some(h)
        }
    };
    let (mut header, mut blobs) = crate::t1b_window::header_for(&raw, &recipe, &recipe_digest, &prep, crate::t1b_lem::Ledger::default(), Vec::new())?;
    if unit_contacts {
        let mut base = vec![i32::MIN; nn];
        let mut top = vec![i32::MIN; nn];
        let mut layer = vec![255u8; nn];
        let envelope = ctx.params.lithology.envelope.is_some();
        let mm = |z: f64| -> R<i32> {
            if !z.is_finite() { return Err("unit contact is non-finite".into()); }
            let q = vz_tilekey::chart::quantize_mm(z * 1000.0);
            if q <= i32::MIN as i64 || q > i32::MAX as i64 { return Err("unit contact exceeds i32 mm (MIN reserved)".into()); }
            Ok(q as i32)
        };
        for i in 0..nn {
            if fixed[i] == MARINE || fixed[i] == LAKE { continue; }
            let j = map(i);
            if !envelope {
                match s.over[j] {
                    crate::t1b_coast::OVER_EDIFICE => { layer[i] = 253; continue; }
                    crate::t1b_coast::OVER_FAN | crate::t1b_coast::OVER_BEACH => { layer[i] = 254; continue; }
                    _ => {}
                }
            }
            if let Some((k, b, t)) = lp.rock.unit_contacts(j, zq[i], envelope) {
                layer[i] = u8::try_from(k).map_err(|_| "unit contact layer index exceeds u8")?;
                base[i] = mm(b)?;
                if let Some(t) = t { top[i] = mm(t)?; }
            }
        }
        for (name, dtype, units, note, bytes) in [
            ("unit_base_mm", "i32", "mm", "absolute elevation of the exposed source unit base; i32::MIN unrecorded", base.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<_>>()),
            ("unit_top_mm", "i32", "mm", "absolute elevation of the exposed source unit top; i32::MIN unrecorded (ordinary first-layer top, water, deposits, edifice)", top.iter().flat_map(|v| v.to_le_bytes()).collect()),
            ("layer", "u8", "index", "source full-column layer index (top=0), envelope indices retained; 253 edifice, 254 deposit, 255 unrecorded; contacts carry their own sentinel", layer),
        ] {
            let file = format!("fields/{name}.{dtype}");
            header.fields.push(crate::t1b_window::FieldEntry { name:name.into(), dtype:dtype.into(), units:units.into(), file:file.clone(), digest:format!("blake3:{}", blake3::hash(&bytes).to_hex()), note:note.into() });
            blobs.push((file,bytes));
        }
        header.approximations.push("unit_contacts_v1 (unit_contacts=1): publication only, same native source-layer cut as exposed_litho; E5B retained envelope contacts agree with detail unit_depth. No procedural bed/flow-unit contacts or infinite last-layer extension. Ordinary continuity source layers are not a claim of legacy T2 column correspondence.".into());
    }
    if owned_area.is_some() && !s.ve_style.is_empty() && !s.ve_unit.is_empty() {
        // Physical inputs survive display-only GUIDE/lake masking. The same
        // native header content-pins them; legacy key-off emits no new fields.
        let fields = [
            ("process_geometry_fixed", "u8", (0..nn).map(|i| s.fixed[map(i)]).collect::<Vec<u8>>()),
            ("process_native_link", "u32", (0..nn).flat_map(|i| s.chan_link[map(i)].to_le_bytes()).collect()),
            ("process_edifice_style", "u8", (0..nn).map(|i| s.ve_style[map(i)]).collect()),
            ("process_edifice_unit", "u16", (0..nn).flat_map(|i| s.ve_unit[map(i)].to_le_bytes()).collect()),
        ];
        for (name,dtype,bytes) in fields {
            let file=format!("fields/{name}.{dtype}");
            header.fields.push(crate::t1b_window::FieldEntry{name:name.into(),dtype:dtype.into(),units:"native input code".into(),file:file.clone(),digest:format!("blake3:{}",blake3::hash(&bytes).to_hex()),note:"Actual native pre-publication input. GUIDE registration and dam-lake display roles do not replace physical geometry, attachment or volcanic material identity. Not qualified geology stock.".into()});
            blobs.push((file,bytes));
        }
    }
    if ctx.params.hs.as_ref().and_then(|h|h.contrast.as_ref()).is_some() {
        // Ship the source columns with the evolved parent, digest-bound by the
        // window header. Detail must not silently reload a different envelope.
        let cfg=ctx.params.lithology.envelope.as_ref().ok_or("E5B envelope missing")?;
        let env=cfg.load()?;
        let mut cells=std::collections::BTreeSet::new();
        for &cell in &raw.t0_cell {
            if let Some(e)=env.cells.get(&cell) {
                for dy in -1i64..=1 { for dx in -1i64..=1 {
                    let y=cell as i64/e.shape[1] as i64+dy;
                    let x=(cell as i64%e.shape[1] as i64+dx).rem_euclid(e.shape[1] as i64);
                    if y>=0 && y<e.shape[0] as i64 { cells.insert((y*e.shape[1] as i64+x) as u32); }
                }}
            }
        }
        let columns:Vec<_>=cells.into_iter().filter_map(|cell|env.cells.get(&cell).map(|e|serde_json::json!({
            "cell":cell,"shape":e.shape,"radius_m":e.radius_m,"mean_m":e.mean_m,"grad":e.grad,
            "layers":e.layers.iter().filter(|l|l.h>0.).map(|l|serde_json::json!({"index":l.index,"code":l.code,"depth":l.depth,"h":l.h})).collect::<Vec<_>>()
        }))).collect();
        header.prep.net.as_mut().ok_or("E5B prep missing")?["evolved_surface"]=serde_json::json!({
            "schema":"e5b-evolved-surface/v1","rule":"2 m bilinear prolongation of the final 32 m process surface, no second construction/relaxation; exact shared parent edges; exposed rock queried from digest-bound native envelope columns; zero new simulated time or material",
            "envelope_digest":cfg.digest,"columns":columns,"chart":recipe.chart,
            "classes":crate::t1b_ct::ct(ctx.params).ok_or("E5B ct missing")?.classes
        });
    }
    if let Some(fine)=hs_stats.as_ref().and_then(|h|h.contrast.as_ref()).and_then(|c|c.subgrid.as_ref()) {
        let base:Vec<u8>=fine.base_mm.iter().flat_map(|v|v.to_le_bytes()).collect();
        let delta:Vec<u8>=fine.delta_xyz.iter().flat_map(|v|v.iter().flat_map(|q|q.to_le_bytes())).collect();
        for (name,bytes,note) in [
            ("escarp_base_mm",base,"Pre-E5 32m surface for conservative fine reconstruction, not final coarse means"),
            ("escarp_delta_xyz",delta,"Sparse global 2m x,y and signed removal mm triples; row order; length is 3 times published_changed_fine_posts, not n*n; absent record means zero removal"),
        ] {
            let file=format!("fields/{name}.i32");
            header.fields.push(crate::t1b_window::FieldEntry{name:name.into(),dtype:"i32".into(),units:"mm / 2m lattice indices".into(),file:file.clone(),digest:format!("blake3:{}",blake3::hash(&bytes).to_hex()),note:note.into()});
            blobs.push((file,bytes));
        }
        let evolved=&mut header.prep.net.as_mut().ok_or("ESCARP-1 prep missing")?["evolved_surface"];
        evolved["schema"]=crate::t1b_escarp::EVOLVED_SCHEMA.into();
        evolved["escarpment"]=serde_json::to_value(fine).map_err(|e|e.to_string())?;
        evolved["rule"]="ESCARP-1 candidate: original 32m surface plus ledgered sparse 2m E5 removal; 32m published heights are restricted means and must not be bilinearly prolonged for the face".into();
        header.approximations.push(fine.rule.clone());
    }
    if let Some(d) = &s.dependency { header.prep.net.as_mut().expect("native net prep").as_object_mut().expect("net prep object").insert("esc_dependency".into(), d.clone()); }
    if let Some(ic) = &ice_now { crate::t1b_window::push_ice_fields(&mut header, &mut blobs, ic); } // Lane GLACIER-1 (gl_now=1 only)
    if let Some(v) = &ice_sub { // Lane GLACIER-3 (gl_geo=1 only)
        let bytes: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
        let file = "fields/ice_sub_mm.i32".to_string();
        header.fields.push(crate::t1b_window::FieldEntry { name: "ice_sub_mm".into(), dtype: "i32".into(), units: "mm".into(), file: file.clone(), digest: format!("blake3:{}", blake3::hash(&bytes).to_hex()), note: crate::t1b_glacial::ICE_SUB_NOTE.into() });
        blobs.push((file, bytes));
    }
    // Lane MS-CEIL (ceil only): ceil_m and surface_src of the published posts (window fixed codes win: non-free = 2).
    if !s.ceil_m.is_empty() {
        let cw = crate::t1b_window::CeilWin { ceil_m: (0..nn).map(|i| s.ceil_m[map(i)]).collect(), surface_src: (0..nn).map(|i| if fixed[i] != FREE { 2 } else { s.ceil_src[map(i)] }).collect() };
        crate::t1b_window::push_ceil_fields(&mut header, &mut blobs, &cw);
    }
    if let Some(h) = prov_rows {
        header.provenance = Some(h);
        header.approximations.push("lane PV (opt-in prov=1): prov_id names the exposed unit per post and header.provenance copies its upstream record (T0 strata layer, volcano record, or T1b present-day deposit); no new process, no draws; unrecorded facts are null (see header.provenance.unknown)".into());
    }
    header.simtime_label = "none: 0.9.0 is a construction (no relaxation, no SimTime)".into();
    header.closure_label = "no sub-resolution closure at 32 m (no roughness draws in the window; keyed draws only in the t1b.rivernet network layer)".into();
    for f in header.fields.iter_mut() {
        match f.name.as_str() {
            "fixed" => f.note = "0 free (hillslope or valley floor), 1 T1a channel (fixed bed, Catmull-Rom through the T1a nodes), 2 marine (zero contour, sea level), 4 lake (stage), 5 river-network channel (t1b.rivernet, fixed bed)".into(),
            "drainage_area_m2" if ctx.params.seam_cert.is_some() => f.note = format!("lane SEAM-CERT ({}): channel posts (1, 5): network {}area; free posts: {}", crate::t1b_seam::ROUTE_RULE, if ctx.params.climate_runoff.is_some() { "effective (Q / r_ref) " } else { "" }, "routing on the whole evaluation lattice, cropped (see the window's seam certificate)"),
            "drainage_area_m2" => f.note = if ctx.params.climate_runoff.is_some() { "lane CLIMATE-RUNOFF effective area Q / r_ref (m2 of reference-runoff catchment): channel posts (1, 5): network effective area (T1a discharge x seconds per year / r_ref; rivernet nodes accumulate a_node x R / r_ref); free posts: window-local routing of a_post x runoff_weight on the published surface (priority flood from fixed posts)".into() } else { "channel posts (1, 5): network drainage area (T1a discharge area / rivernet accumulated area); free posts: window-local routing area on the published surface (priority flood from fixed posts)".into() },
            "base_mm" => f.note = "the eikonal (hillslope) surface before crest rounding and fill".into(),
            "initial_mm" => f.note = "same as base_mm (no relaxation; kept for the window schema)".into(),
            "uplift_m_per_my" => f.note = "E: steady-state erosion rate field (hat blend of the solved per-T0-cell table), m/My; a construction parameter, not a rate history".into(),
            "roughness_mm" => f.note = "0 (no draws in 0.9.0 windows)".into(),
            "inflow_area_m2" => f.note = "0 (no ring inflow in 0.9.0)".into(),
            _ => {}
        }
        // Field digests stay those of the bytes (notes are header text only).
    }
    header.preserved = vec![
        "T1a channel==1 nodes fixed at their trunk levels (F1), planform through the T1a nodes".into(),
        "T0 water: zero-contour sea (sea level) and T0 lake cells (stage)".into(),
        "the chart-wide t1b.rivernet network and E table (identical for every window): channels, beds and relief level are the same in every window".into(),
    ];
    header.invented = vec![
        "network below the T1a channels (invasion percolation, keyed), jitter and meander (keyed), in the t1b.rivernet layer".into(),
        "channel beds (steady-state stream power), valley floors (declared width law), hillslopes (steady-state creep gradient by rock, eikonal), crest rounding (linear-diffusion curvature)".into(),
    ];
    header.approximations = vec![
        "construction, not a process: no mass balance (F7 not applicable), no SimTime".into(),
        "E is one scalar per T0 cell (hat-blended), solved on a 128 m chart lattice; the 32 m means differ by the resolution effect (reported)".into(),
        "first-order fast marching (4-neighbour) eikonal; rock evaluated at the upwind height".into(),
        "no climate input (ADR 043 OPEN 6); no glacial or coastal processes (OPEN 5); coast = T0 bilinear zero contour".into(),
        "fidelity tier Preview; pre-signoff substrate".into(),
    ];
    if ctx.params.coast.is_some() || ctx.params.volcano.is_some() {
        header.simtime_label = "none: 0.9.1 is a construction (no relaxation, no SimTime); the volcano age is the record's own".into();
        header.approximations.retain(|a| !a.starts_with("no climate input"));
        header.approximations.push("no climate input (ADR 043 OPEN 6); coast by construction (lane CO 0.9.1): lowstand shelf, per-cell flooding offset solved for F4, bayhead fans, class wave cut; no waves, tides or longshore drift".into());
        header.invented.push("lane CO: flooding offset per coastal T0 cell (solved), fan radii (declared supply law), wave-cut retreat and beach widths (declared), edifice shape by kind and degradation by record age (declared), vent moved off T1a channels (declared rule)".into());
        header.preserved.push("lane CO: T0 volcano record points (position up to the placement rule, kind, subgrid summit, last eruption); F4 land fraction per coastal T0 cell (solved)".into());
    }
    if ctx.params.fx.is_some() {
        header.simtime_label = "none: 0.9.3 is a construction (no relaxation, no SimTime); the volcano age is the record's own".into();
        header.invented.push("lane TD (0.9.3): capture routing of perched network nodes (steady-state divide balance), soft valley floors, hillslope length over hillslope posts, lithology contacts inside one eikonal step, soft-rock coastal bluffs at the class limiting slope (rules in the recipe's params.fx)".into());
        header.approximations.push("lane TD: divide balance is enforced node by node at 512 m (capture criterion with the pass-1 hillslope gradient); channel curves can still pass closer than the node spacing, so a residual of perched 32 m channel posts is reported, not hidden".into());
    }
    if let Some(t) = crate::t1b_ls::tb(ctx.params) {
        if t.t1a_levels != 0 { header.preserved[0] = "T1a channel==1 nodes and planform fixed (F1); T1a LEVELS solved by the bed law with the network (lane TB S1, R2 EXPLORATORY, ruling pending)".into(); }
        header.approximations.push(format!("lane TB options (exploratory / measurement, rulings pending): t1a_levels {} steep_eik {} r4_plain_permille {} (rules in the recipe's params.ls.tb)", t.t1a_levels, t.steep_eik, t.r4_plain_permille));
        if let Some(f) = t.tc.as_ref().and_then(|c| c.tf.as_ref()) { header.approximations.push(format!("lane TF options (opt-in, lane TF PREREG): report {} bluff_v2 {} subpost {} exposure {} (fetch {} rays, cap {} mm) cone_c1 {} (rules in the recipe's params.ls.tb.tc.tf); closed basins stay filled (no lake model; reported in prep.net)", f.report, f.bluff_v2, f.subpost, f.exposure, f.fetch_dirs, f.fetch_cap_mm, f.cone_c1)); }
        if let Some(g) = t.tc.as_ref().and_then(|c| c.tg.as_ref()) { header.approximations.push(format!("lane TG options (opt-in, lane TG PREREG): shore {} (cap {} mm) sill_grade {} smooth_bil {} (half {} mm on a {} mm global grid) f_deadband {} ({} permille) (rules in the recipe's params.ls.tb.tc.tg); shore_mm is presentation geometry derived from this surface (no new process)", g.shore, g.shore_cap_mm, g.sill_grade, g.smooth_bil, g.smooth_half_mm, g.smooth_grid_mm, g.f_deadband, g.f_deadband_permille)); }
        if let Some(v) = t.tc.as_ref().and_then(|c| c.ve.as_ref()) { header.approximations.push(format!("lane VE options (opt-in, lane VE PREREG): on {} (rules in the recipe's params.ls.tb.tc.ve): eruption styles by record kind and a keyed statistical history of the youngest mappable units (the planet record has no per-eruption detail); units are process-shaped constructions (no SimTime, no mass balance); paths follow the T1b network curves; crater, caldera, maar and pit depressions are filled by the final fill (declared lakes, listed in prep.net.coast.ve)", v.on)); }
        if let Some(v) = t.tc.as_ref().and_then(|c| c.ct.as_ref()).and_then(|c| c.v2.as_ref()) { header.approximations.push(format!("lane STEP2 (opt-in ct sub-keys): apply 2 = ct_remap_v1 (monotone height remap of the rounded v0 surface), var {} (0 regular P, 1 keyed closure P x U({}, {}) permille, 2 strata packages else closure; unit_min {} unit_max {} permille), base {} (ct_base_mm), local {} (1: regime and split from the local gradient of the rounded surface); rules in the recipe's params.ls.tb.tc.ct.v2; strata packages use the lithology cut's T0 column, not a lava-flow record", v.var, v.var_lo_permille, v.var_hi_permille, v.unit_min_permille, v.unit_max_permille, v.base, v.local)); }
        if let Some(b) = t.tc.as_ref().and_then(|c| c.ct.as_ref()).and_then(|c| c.v2.as_ref()).and_then(|v| v.bed.as_ref()) { header.approximations.push(format!("lane BENCH-VAR (opt-in ct_var=3): ct_bed_v1 per-bed keyed step thickness P x {} x exp({} z) permille, clamp {}..{}, hard share U({}, {}) permille, beds anchored at each strata package's top (purpose t1b.rivernet.ct.bed; rule in params.ls.tb.tc.ct.v2.bed); a closure for flows/beds the T0 column does not resolve, not a record", b.median_permille, b.sigma_permille, b.min_permille, b.max_permille, b.eta_lo_permille, b.eta_hi_permille)); }
        if let Some(r) = t.tc.as_ref().and_then(|c| c.ct.as_ref()).and_then(|c| c.v2.as_ref()).and_then(|v| v.tread.as_ref()) { header.approximations.push(format!("lane TREAD2 (opt-in ct_tread): ct_tread_v1 3-limb tread/talus/cliff periods on layered classes (tread {} permille tangent, clamp {}..{}, dip {}; hard share eta {} permille closure where the T0 column has one class through the period, cap {}; transition to the full profile at {} permille gradient); a construction (no bench process; rules in the recipe's params.ls.tb.tc.ct.v2.tread); ct_regime 3 = tread", r.tread_tan_permille, r.tread_min_tan_permille, r.tread_max_tan_permille, r.dip, r.eta_permille, r.eta_max_permille, r.full_tan_permille)); }
        if let Some(c) = t.tc.as_ref().and_then(|c| c.ct.as_ref()) { header.approximations.push(format!("lane STEP (opt-in ct, law threshold-slope-lithology-v1): bare-rock hillslopes split into cliff and talus with the v0 mean gradient (hillslope relief conserved); apply {} (1 superposed on the rounded v0 surface), window_only {} (1: measurement, the layer was solved without ct), layered step period {} mm from datum {} mm, crest_exempt {}, apron {}, rockmass_cliff {} (rules and table in the recipe's params.ls.tb.tc.ct); a construction: no rockfall process, talus volume is not tracked", c.apply, c.window_only, c.step_mm, c.datum_mm, c.crest_exempt, c.apron, c.rockmass_cliff)); }
        if let Some(h) = t.tc.as_ref().and_then(|c| c.th.as_ref()) { header.approximations.push(format!("lane TH options (opt-in, lane TH PREREG): expo_smooth {} (length {} mm) uncut {} (rules in the recipe's params.ls.tb.tc.th); uncut_mm is the surface before the wave cut on cliff-class posts (no new process)", h.expo_smooth, h.expo_len_mm, h.uncut)); }
        if let Some(c) = &t.tc { header.approximations.push(format!("lane TC options (exploratory / measurement): floor {} clip {} r4n_permille {} pert_ppm {} pert_mode {} freeze_from {} curve {} (rules in the recipe's params.ls.tb.tc)", c.floor, c.clip, c.r4n_permille, c.pert_ppm, c.pert_mode, c.freeze_from, c.curve)); }
    }
    if let Some((_, c)) = &seam_route {
        header.approximations.push(format!("lane SEAM-CERT (opt-in seam_cert=1, {}): construction fills seeded only from fixed posts (no lattice-ring seed), free-post routing on the evaluation lattice; certificate {} (fill-tainted published posts {}, route-bad published free posts {}, lake extent {} / donor extent {} posts outside the published box vs ring {} posts)",
            crate::t1b_seam::LABEL, if c.pass { "pass" } else { "FAIL" }, c.fill_tainted_published_posts, c.route_bad_published_free_posts, c.fill_lake_extent_published_posts, c.route_up_extent_published_posts, c.ring_distance_posts));
    }
    if ctx.params.hs.as_ref().and_then(|h|h.contrast.as_ref()).is_some() {
        header.approximations.push("E5B: ct_regime includes resolved post-contrast cliff faces on final exposed cliff-forming rock. Threshold is the ct class cliff tangent; slope is max absolute D8 edge gradient on final mm heights. Fixed/overwritten posts retain their domain masks. Other ct_* remain pre-process construction diagnostics and must not re-profile the evolved surface. Removed rock is assigned per source contact to e5_unresolved_external_sediment; zero talus storage, uncalibrated efficient evacuation, geometric volume not collectible mass.".into());
    }
    if let Some(h) = &hs_stats {
        header.simtime_label = crate::t1b_hs::SIMTIME.into();
        header.approximations.push(format!("lane HS-PASS (opt-in hs=1, {}): hillslope-only transient, N {} steps of dt {} yr (T {} yr) on free posts after the construction; creep, threshold landslide with talus, MFD stream power above A S^2 > A_1 (rules in the recipe's params.hs); construction-stage fields (ct_*, edifice_*, over, base_mm) describe the surface before the pass; dependency certificate {} (route up {} / down {} / fill {} posts vs ring {} posts)",
            crate::t1b_hs::LABEL, h.steps, h.dt_yr, h.t_yr, if h.certificate_pass { "pass" } else { "FAIL" }, h.cert_route_up_posts, h.cert_route_down_posts, h.cert_fill_down_posts, h.cert_ring_distance_posts));
        header.invented.push("lane HS-PASS: dissection of the constructed hillslopes over a declared bounded time (gullies, talus aprons, rounded crests); no noise, every change a flow or failure path to a pinned post".into());
    }
    if let Some(c) = crate::t1b_cl::cl(ctx.params) {
        header.approximations.push(format!("lane CHANNEL-CENTRELINE (opt-in cl=1, {}): T1a links follow the T1a planform curve re-evaluated from the published segment records instead of a Catmull-Rom curve through the 512 m lattice nodes (rules in the recipe's params.cl); the T1a planform is T1a's declared keyed stand-in, not a physical meander model", crate::t1b_cl::LABEL));
        if let Some(m) = &c.mq {
            header.approximations.push(format!("lane MEANDER-Q (opt-in cl_mq=1, {}): T1a links follow a kinematic meander-migration model run per T1a reach over a STAND-IN migration time of {} yr (rules in the recipe's params.cl.mq; per-reach statistics and oxbow features in prep.net.cl.mq); kinematic (curvature-driven, upstream-lagged), not hydrodynamic: no flow field, bank shear, sediment budget or water routing", crate::t1b_mq::LABEL, m.age_yr));
            header.invented.push("lane MEANDER-Q cl_mq=1: initial lateral perturbations (world-keyed draws), cutoff criteria, regime transition band and valley-wall factor are declared closures; floodplain = swept envelope + margin; floor inset by the hydraulic-geometry depth".into());
        } else if c.sine_on() {
            header.invented.push(format!("lane CHANNEL-CENTRELINE cl_mn=1: DECLARED DECORATIVE meander belt (keyed sine-generated bends, wavelength {} channel widths, sinuosity {}-{} capped by floor confinement); NOT meander evolution (no flow, bank erosion, deposition or cutoff); stands in until RIVER-PILOT / RIVER-INTEGRATE",
                c.wavelength_widths_permille as f64 / 1000.0, c.sinuosity_min_permille as f64 / 1000.0, c.sinuosity_max_permille as f64 / 1000.0));
        }
    }
    if let Some(d) = &dam_stats {
        header.approximations.push(format!("lane DAMMED-VALLEY (opt-in dam_lake=1, {}): {}; {} lakes published ({} posts), {} dry by the balance (terminal lakes below spill not modelled), {} hollows wholly on edifices left as dry fills (see prep.net.dam_lake)",
            crate::t1b_dam::LABEL, crate::t1b_dam::RULE, d.lakes_published, d.lake_posts_published, d.dry_by_balance, d.edifice_hollows));
        if ctx.params.dam_lake.as_ref().map(|p| p.crater == 1).unwrap_or(false) {
            header.approximations.push(format!("lane EDIFICE-FOOT (opt-in dam_crater=1): {}; {} crater lakes published (counted in the lakes above), {} crater hollows dry by the balance", crate::t1b_dam::CRATER_RULE, d.crater_lakes_published, d.crater_dry_by_balance));
        }
    }
    if let Some(q) = ctx.params.pedestal.as_ref().filter(|q| q.base == 1) {
        let _ = q;
        let g: Vec<String> = ctx.coast.volcanoes.iter().filter(|v| v.ped_shift_m.is_some()).map(|v| format!("{} base {:.1} m ({:.1} m below the bilinear), H {:.1} m", v.id, v.base_m, v.ped_shift(), v.h_m)).collect();
        header.approximations.push(format!("lane T1B-PEDESTAL (opt-in ped_base=1, {}): {}; layer-published local ground used: {}", crate::t1b_ped::LABEL, crate::t1b_ped::BASE_RULE, g.join("; ")));
    }
    if let Some(f) = &ctx.params.edifice_foot {
        let c = s.coast.ef.clone().unwrap_or_default();
        header.approximations.push(format!("lane EDIFICE-FOOT (opt-in ef_cap=1, {}, k {} permille, t_fill {} mm): {}; evaluation lattice: {} posts capped, body lowered by at most {:.1} m, {:.3} km3 (see prep.net.coast.ef)",
            crate::t1b_ef::LABEL, f.k_permille, f.fill_mm, crate::t1b_ef::RULE, c.capped_posts, c.cut_max_m, c.cut_volume_m3 / 1e9));
    }
    if let Some(c) = ctx.params.ceil.as_ref() {
        header.approximations.push(format!("lane MS-CEIL (opt-in ceil={}, c_h {}; prototype): {}; evaluation lattice: {} posts lowered to the ceiling, by at most {:.1} m (fields ceil_m, surface_src)", c.c_h_permille, c.c_h(), CEIL_RULE, s.stats.ceil_posts, s.stats.ceil_lowered_max_m));
        header.invented.push("lane MS-CEIL: the upland ceiling level c_h x R above the T0 bilinear (one declared constant; the T0 record gives the trunk level below the cell mean, not the upland level above it)".into());
        // Lane MS-CEIL-2: the sub-rules (only when on, so ceil alone keeps its v1 window header).
        if c.bed == 1 { header.approximations.push(format!("lane MS-CEIL-2 (opt-in ceil_bed=1; prototype): {CEIL_BED_RULE}")); }
        if c.ped == 1 { header.approximations.push(format!("lane MS-CEIL-2 (opt-in ceil_ped=1; prototype): {CEIL_PED_RULE}")); }
    }
    raw.digest = pf_contracts::artifact::digest_canonical(&header).map_err(|e| e.to_string())?;
    Ok(crate::t1b_window::WindowRun { raw, header, prep, blobs, seconds_prepare, seconds_relax })
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Lane MS-CEIL: key off serializes no ceil field (recipe digests unchanged); the key validates its range and rule.
    #[test]
    fn ceil_key_off_absent_and_validated() {
        let p = RiverNetParamsV1::v1(true);
        assert!(!serde_json::to_string(&p).unwrap().contains("\"ceil\""));
        let mut q = p.clone();
        q.ceil = Some(CeilParamsV1::v1(500));
        assert!(q.validate().is_ok());
        assert!(serde_json::to_string(&q).unwrap().contains("\"c_h_permille\":500"));
        q.ceil = Some(CeilParamsV1::v1(-1));
        assert!(q.validate().is_err());
        q.ceil = Some(CeilParamsV1 { rule: "x".into(), c_h_permille: 500, ..CeilParamsV1::v1(500) });
        assert!(q.validate().is_err());
        assert!((CeilParamsV1::v1(250).c_h() - 0.25).abs() < 1e-12);
    }
    /// Lane MS-CEIL-2: ceil alone serializes exactly as v1 (no sub-key fields); the sub-keys validate their rules.
    #[test]
    fn ceil_sub_keys_off_absent_and_validated() {
        let c = CeilParamsV1::v1(500);
        assert_eq!(serde_json::to_string(&c).unwrap(), format!("{{\"rule\":{},\"c_h_permille\":500}}", serde_json::to_string(CEIL_RULE).unwrap()));
        let mut b = c.clone();
        b.set_bed(1);
        b.set_ped(1);
        assert!(b.validate().is_ok());
        let j = serde_json::to_string(&b).unwrap();
        assert!(j.contains("\"bed\":1") && j.contains("\"ped\":1") && j.contains("ceiling_bed_grade_v1") && j.contains("ceiling_free_pedestal_v1"));
        let back: CeilParamsV1 = serde_json::from_str(&j).unwrap();
        assert_eq!(back, b);
        b.set_bed(0);
        assert_eq!(serde_json::to_string(&b).unwrap().contains("\"bed\""), false);
        let mut x = c.clone();
        x.bed = 1;
        assert!(x.validate().is_err());
        x.bed_rule = Some("y".into());
        assert!(x.validate().is_err());
        let mut y = c.clone();
        y.set_bed(2);
        assert!(y.validate().is_err());
    }
    #[test]
    fn layer_at_names_the_full_column_layer_that_at_cuts() {
        // one cell: shale 10 m, sandstone 0.5 m (below the 1 m cut threshold), basalt 100 m, granite 5000 m
        let cols = StrataCols { off: vec![0, 4], lith: vec![5, 6, 3, 1], th: vec![10.0, 0.5, 100.0, 5000.0], grad_e: vec![0.0], grad_n: vec![0.0], digests: vec![], cols: Default::default(), lcont: Default::default() };
        let table = crate::t1b_litho::LithologyV1::standin_v1();
        let mut sc = [0.0; 12];
        let mut k = [0.0; 12];
        for c in 0..12u8 { sc[c as usize] = table.sc(c); k[c as usize] = table.k(c); }
        let rock = RockPosts::Strata { table: &table, cols: &cols, cell: vec![0], datum: vec![1000.0], fallback: vec![0], sc, k, volc: None, cont: None };
        use crate::t1b_litho::ExposedLithology;
        for (z, want_k, want_code) in [(1005.0, 0u16, 5u8), (995.0, 0, 5), (985.0, 2, 3), (500.0, 3, 1), (-9000.0, 3, 1)] {
            assert_eq!(rock.layer_at(0, z), Some((0, want_k)), "z {z}");
            assert_eq!(rock.at(0, z).2, want_code, "z {z}");
            assert_eq!(cols.lith[rock.layer_at(0, z).unwrap().1 as usize], rock.at(0, z).2, "the id names the layer whose code is published");
        }
        let empty = StrataCols { off: vec![0, 0], lith: vec![], th: vec![], grad_e: vec![0.0], grad_n: vec![0.0], digests: vec![], cols: Default::default(), lcont: Default::default() };
        let r2 = RockPosts::Strata { table: &table, cols: &empty, cell: vec![0], datum: vec![1000.0], fallback: vec![7], sc, k, volc: None, cont: None };
        assert_eq!(r2.layer_at(0, 900.0), None);
        assert_eq!(r2.at(0, 900.0).2, 7);
    }
    #[test]
    fn fmm_point_source_is_a_round_cone() {
        let lat = EvalLat { gx0: -50, gy0: -50, n: 101, dx_mm: 32_000 };
        let nn = lat.nodes();
        let mut z = vec![0.0; nn];
        let mut known = vec![false; nn];
        let c = lat.idx(0, 0).unwrap();
        known[c] = true;
        let un = fmm(&lat, &mut z, &mut known, &|_, _| 0.5);
        assert_eq!(un, 0);
        // Along the axis exact; along the diagonal within 8 % (first-order FMM).
        let ax = z[lat.idx(40, 0).unwrap()];
        assert!((ax - 0.5 * 40.0 * 32.0).abs() < 1e-9);
        let dg = z[lat.idx(28, 28).unwrap()];
        let exact = 0.5 * 32.0 * (2.0f64 * 28.0 * 28.0).sqrt();
        assert!((dg - exact).abs() / exact < 0.08, "diag {dg} vs {exact}");
        // Every post has a strictly lower 4-neighbour (no pits) except the source.
        for i in 0..nn { if i == c { continue; } let (x, y) = (i % 101, i / 101);
            let lower = (x > 0 && z[i - 1] < z[i]) || (x < 100 && z[i + 1] < z[i]) || (y > 0 && z[i - 101] < z[i]) || (y < 100 && z[i + 101] < z[i]);
            assert!(lower); }
    }
    #[test]
    fn soft_source_lowered_only_by_accepted_fronts() {
        // Lane TB soft_channel_cap_v1: a soft post (bed 50) next to a known post (z 0, label 7) is lowered to the arrival
        // (0 + speed x dx) only when the front's label is accepted.
        let lat = EvalLat { gx0: 0, gy0: 0, n: 3, dx_mm: 32_000 };
        for (ok, want) in [(true, 0.5 * 32.0), (false, 50.0)] {
            let mut z = vec![f64::NAN; 9];
            let mut known = vec![false; 9];
            let mut lab = vec![u32::MAX; 9];
            z[3] = 0.0; known[3] = true; lab[3] = 7;
            z[4] = 50.0; lab[4] = 9;
            let mut soft = vec![false; 9];
            soft[4] = true;
            let acc = move |j: usize, l: u32| -> bool { j != 4 || (ok && l == 7) };
            let un = fmm_lab_soft_acc(&lat, &mut z, &mut known, &|_, _| 0.5, &mut lab, Some(&soft), Some(&acc));
            assert_eq!(un, 0);
            assert!((z[4] - want).abs() < 1e-9, "ok {ok}: {}", z[4]);
        }
    }
    #[test]
    fn fill_raises_a_pit_to_its_spill() {
        let lat = EvalLat { gx0: 0, gy0: 0, n: 5, dx_mm: 32_000 };
        let mut z = vec![10.0; 25];
        z[12] = 1.0;
        let outlet = vec![false; 25];
        let (r, m) = fill(&lat, &mut z, &outlet);
        assert_eq!(r, 1);
        assert_eq!(z[12], 10.0);
        assert_eq!(m, 9.0);
    }
    #[test]
    fn chaikin4_keeps_ends_and_monotone_beds() {
        let p = vec![[0.0, 0.0, 100.0, 1.0], [512.0, 0.0, 90.0, 2.0], [512.0, 512.0, 80.0, 3.0], [1024.0, 512.0, 70.0, 4.0]];
        let mut q = p.clone();
        for _ in 0..2 { q = chaikin4(&q); }
        assert_eq!(q[0], p[0]);
        assert_eq!(*q.last().unwrap(), *p.last().unwrap());
        for w in q.windows(2) { assert!(w[1][2] <= w[0][2]); }
    }
}
