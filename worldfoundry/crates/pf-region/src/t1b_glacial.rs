//! Lane GLACIAL-T1 (2026-10-21): an OPT-IN glacial overprint of the T1b rivernet construction (key
//! `glacial=<t0.glacial dir>` on both `rivernet` and `rnwin`). LABEL: pre-signoff, fidelity tier Preview, derived, never
//! canon; a CONSTRUCTION (no SimTime): nothing here is integrated in time and nothing may be cited as an age.
//!
//! WHAT IT READS. The `t0.glacial` layer published with the world (lane GLACIAL-T1's gt1-publish driver; schema
//! `vz.t0-glacial/v1`): per T0 cell the ELA of each of the 8 coldhouse glacial-cycle phases, the cycle-mean ELA, the
//! glacial-maximum ELA and the genesis precipitation, plus the T0 law constants (K_g/K_f, Hf - He, phase precipitation
//! factors) of `glacial-stream-power-blend-v1` (lane GLACIAL-T0, Hergarten 2021 ESurf 9:937; Liebl et al. 2023 GMD
//! 16:1315). T1b never re-derives the climate.
//!
//! SCENARIO (declared). The overprint is the steady landform of the world's PRESENT epoch, taken to be a coldhouse with
//! glacial cycles (as Earth's Quaternary): the 8 cycle phases are equally weighted. The T0 genesis also had warmer states;
//! they shaped the T0 cell means (which the F6 E solve still meets) but not the valley-scale form.
//!
//! CHART-WIDE GLACIER SOLVE (once per chart, on the 512 m network; `glacier_table_v1`). For every land and T1a node of
//! the published network (receivers: land -> network receiver, T1a -> T1a receiver; water ends the ice):
//! - production (Hergarten 2021 eq 23 above the ELA, Liebl 2023 eq 11 melt below, the same function as the T0 law):
//!   p_i = P f_P s(z), s = min((z - He)/(Hf - He), 1) above He, -(x + x^2/3) with x = (He - z)/(Hf - He) below; z = the
//!   node's FLUVIAL bed (the 0.9.x bed law extended to hillslope nodes: a lower bound of the node's mean surface,
//!   declared); node own area 512 m x 512 m; P and He bilinear between T0 cell centres.
//! - ice flux by continuity down the tree, clamped at 0 (the terminus; Hergarten eq 21): Q_i(n) = max(0, p_i a + sum Q_i
//!   of donors); ice share gamma = Q_i / (precipitation flux through n).
//! - Hergarten blend over the 8 phases: F = mean_k [gamma_k^m r + (1 - gamma_k)^m], r = K_g / K_f (T0 manifest), m the
//!   T0 stream-power exponent (laws). In steady state E = K_eff A^m S^n with the T1b bed law's n = 1 form, so the ICE
//!   SURFACE rises upstream at S_s = S_f / F (S_f the fluvial bed slope of the link): flatter where ice erodes.
//! - two geometry states, CHAR (cycle-mean ELA, f_P = mean phase factor) and MAX (glacial-maximum ELA and f_P):
//!   width W = w_ref (Q_ice / q_ref)^0.3 (Liebl 2023 Eq. 8 with their calibration "an ice flux of 1e4 m3/a generates a
//!   glacier width of 150 m", alpha 0.3 from Bahr 1997 / Hergarten 2021); thickness H = min(delta W, H_plastic):
//!   delta = h / w = 0.25 (Hergarten 2021 default adopted by Liebl 2023, Eq. 9) and the perfectly plastic limit
//!   H_plastic = tau_y / (f rho_i g S_s) with Nye's (1965, Table IV) shape factor f of a parabolic channel looked up on
//!   the half-width / depth ratio (0.445, 0.646, 0.746, 0.806 at 1, 2, 3, 4; clamped), tau_y = 100 kPa (paleo valley
//!   glaciers 50-150 kPa, Benn and Hulton 2010, Pellitero 2016 GlaRe).
//! - CONCORDANT ice surfaces (MacGregor et al. 2000; Anderson et al. 2006): from each terminus (the snout sits on the
//!   fluvial bed) z_s rises upstream at S_s; at a confluence both donors start from the receiver's surface, so the bed
//!   z_s - H of a thinner tributary hangs above the trunk bed by H_trunk - H_trib: HANGING VALLEYS are a consequence.
//!   Thickening ice upstream lowers the bed below its receiver: OVERDEEPENINGS, limited by the adverse-slope rule of
//!   Alley et al. 2003 (Nature 424:758): the bed never falls upstream faster than psi x S_s, psi = 1.5.
//!   T1a nodes are anchors (their levels are T0-derived seeds, ADR 043 preserved list): the glacier there sits on the
//!   T1a level (no long-profile carve), so T1a valleys are widened (U section) but not deepened (declared).
//! - imprint (`occupation_imprint_v1`): I = fraction of the 8 cycle phases with ice at the node; carve = I x max over
//!   CHAR, MAX of max(0, fluvial bed - glacial bed); the state with the larger carve sets the section (W, H).
//! - channel beds of land nodes become the tree-filled glacial bed (a node never stands below its receiver: the
//!   overdeepening is sediment/water filled to its spill, as the construction fill already does), so F1 holds.
//! - cirques (`ela_cirque_v1`) at the heads of CHAR glaciers (ice at the node, no ice in any donor): floor at the
//!   cycle-mean ELA (Porter 1989: cirque floors approximate the characteristic ELA) but no deeper than D_c below the head
//!   bed, and not below the lip (the receiver's filled bed) minus psi S_s R (Alley); bowl radius R = 390 m and depth scale
//!   D_c = 0.37 x 2R (Barr et al. 2019: 2208 cirques, L 774 +- 426 m, H/L ~ 0.37). Project assumptions awaiting review.
//!
//! WINDOW CARVE (pointwise reads of the published table; `u_section_v1`). After crest rounding and the ct remap, before
//! the valley-floor fans, the bank rule and the fill: each free or floor post within reach of a glaciated link takes
//! z = min(z, floor(t) + H(t) (2 |y| / W(t))^b), b = 2 (Svensson 1959; Harbor 1992: erosion drives b toward 2), y the
//! distance to the link polyline (node -> meander midpoint -> receiver), floor/H/W linear along the link; the min of two
//! continuous surfaces is continuous, the power law steepens the walls until it meets the pre-glacial slope (truncated
//! spurs). Cirque bowls: z = min(z, z_c + D_c (r/R)^2) on free posts, and channel posts in the bowl take
//! min(z, max(bowl, flat)) (flat = max(z_c, lip)) so channels stay monotone; the final fill makes tarns flat floors.
//! Every quantity a post reads is in the published chart table, so windows cannot disagree (F9).
//!
//! NOT DONE (listed as ideas): moraines and outwash, wall eikonal (aretes, horns), lakes kept as water (ADR 043 B2),
//! fjords below the lowstand, T1a trunk deepening (needs a T1a glacial key), aspect-dependent ELA.
//!
//! LANE GLACIER-1 (2026-10-22; design review GLACIER-STAGE design section 3.4, stages GS-A + GS-B1 + the simplest GS-B2).
//! - PRESENT ICE ONLY. A t0.glacial layer whose manifest declares `scenario.present_ice_only: true` (lane GLACIER-1's
//!   gl-sidecar: coldhouse duty 0, every phase = the interglacial ELA) sets `present_ice_only = 1` in the params: no
//!   sections, no long-profile carve, no cirques (the table's w, h, carve are 0, so `apply_beds` and the window carve are
//!   no-ops and the rock surface is the key-off one). Whether the world had glacial periods is an open operator decision.
//! - NOW STATE (key `gl_now=1`, needs present_ice_only): the present (interglacial) ELA from the layer's
//!   `ela_interglacial_m`, f_P = 1, the same production and continuity routing as CHAR / MAX; per node W (Liebl) and
//!   H = min(0.25 W, Nye plastic at the FLUVIAL bed slope: the present ice sits on the uncarved bed), surface = floor + H.
//! - ICEFIELD SHEET (`plastic_sheet_v1`, once per chart on the 128 m solve lattice, on the pre-glacial evaluation
//!   surface, cone included): ice mask M = land posts with bed >= ELA_now (accumulation area) OR within W/2 of a link
//!   whose node carries Q_now > 0 (outlet tongues, W linear to 0 at an ice-free receiver); surface z_s of a perfectly
//!   plastic sheet |grad z_s| = h0 / (z_s - b), h0 = tau_y / (rho_i g) (Reeh 1982; Benn and Hulton 2010 in 2-D), z_s = b
//!   on the margin of M, marched inward by Dijkstra on 8 neighbours with the exact local update (u - z_j)(u - b_i) = h0 d.
//! - WINDOW FIELDS (pointwise reads of the published sheet; windows never recompute it): ice_now_mm, ice_surface_mm,
//!   ice_class 0 none / 1 accumulation (surface >= ELA_now) / 2 ablation.

use crate::t0::T0;
use crate::t1b::{CHANNEL, FREE, GUIDE};
use crate::t1b_rivernet::{EvalLat, Network, K_LAND, K_T1A};
use crate::{NODES, NONE, OFF_CHART, R};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::path::Path;

pub const GL_GENERATOR_VERSION: &str = "0.9.4+glacial-overprint-v1";
pub const GL_SCHEMA_T0: &str = "vz.t0-glacial/v1";
pub const GL_PHASES: usize = 8;
/// Window field note for gl_carve_mm.
pub const CARVE_NOTE: &str = "lane GLACIAL-T1 (T1b, opt-in glacial key): how far the glacial window carve lowered the post, in mm (0 = untouched): U section min(z, floor + H (2|y|/W)^2) on free and floor posts near glaciated links on network channel posts (code 5) the same section on the tree-filled floor (so they stay monotone; T1a channel posts are never carved) plus cirque bowls, applied after crest rounding and the ct remap and before fans, the bank rule and the fill (so later raises, e.g. the fill of an overdeepening, are not subtracted); the long-profile carve of the channel beds is in the layer table (gl_carve_m), not here";
/// Lane GLACIER-1 window field notes (gl_now=1 only).
pub const ICE_NOW_NOTE: &str = "lane GLACIER-1 (T1b, opt-in glacial=<present-ice-only t0.glacial> gl_now=1): present-day ice thickness over the published height, mm (0 = no ice): ice_window_v1 = max(0, min(z_s - z, H + 200 m), H) where the bilinear of the layer's icefield sheet mask >= 0.5 and the post is not water; z_s = bilinear of the plastic_sheet_v1 surface (the lowest mask-corner surface at a 128 m cliff, corners > 300 m apart), H = bilinear of the sheet thickness (128 m chart lattice, tau_y 100 kPa); 'present ice only': no glacial-period scenario declared, no troughs carved; derived, never canon";
pub const ICE_SURFACE_NOTE: &str = "lane GLACIER-1: present-day ice surface = height_mm + ice_now_mm (= height_mm where no ice), mm";
pub const ICE_CLASS_NOTE: &str = "lane GLACIER-1 ice_class_v1: 0 no ice, 1 accumulation area (firn: ice surface >= present ELA, Ohmura 1992 from seed-climate-v1), 2 ablation area (bare ice: surface < present ELA); codes 3 debris-covered, 4 crevasse zone, 5 icefall, 6 lava window, 7 ice-dammed lake, 8 moraine, 9 geothermal crater are reserved (design review GLACIER-STAGE B.3) and not produced";

/// Published table fields (f32 per node unless noted).
pub const GL_FIELDS: [&str; 13] = ["gl_s_surface", "gl_q_char_m3_yr", "gl_q_max_m3_yr", "gl_f_cycle", "gl_occupation", "gl_w_m", "gl_h_m", "gl_bed_m", "gl_bed_fill_m", "gl_carve_m", "gl_zs_char_m", "gl_zs_max_m", "gl_ela_char_m"];

/// Declared constants (integers and rule strings; exact canonical JSON). Serialized into a recipe only when the key is on.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct GlacialParamsV1 {
    /// blake3 of the t0.glacial manifest (glacial.json) and of its field files, combined (the upstream digest).
    pub t0_glacial_digest: String,
    /// Liebl 2023 Eq. 8: w = width_ref x (Q_ice / width_ref_flux)^(width_exponent).
    pub width_ref_mm: i64,
    pub width_ref_flux_m3_per_yr: i64,
    pub width_exponent_permille: i64,
    /// Liebl 2023 Eq. 9 / Hergarten 2021: h / w.
    pub thickness_width_permille: i64,
    /// Perfectly plastic basal yield stress (kPa).
    pub yield_stress_kpa: i64,
    pub ice_density_kg_m3: i64,
    /// Alley et al. 2003 adverse bed slope limit as permille of the ice-surface slope.
    pub adverse_slope_permille: i64,
    /// Power b of the U section z = floor + H (2|y|/W)^b (permille).
    pub u_power_permille: i64,
    /// Evaluation reach of the section as permille of W (beyond it the pre-glacial surface is always lower in practice).
    pub u_reach_permille_of_width: i64,
    pub cirque_radius_mm: i64,
    pub cirque_depth_permille_of_length: i64,
    /// Numerical floor on the ice-surface slope (ppm) for the plastic limit.
    pub min_surface_slope_ppm: i64,
    /// Duration of the present coldhouse epoch (kyr): the time cap of the carve (Quaternary 2.58 Ma, ICS chart; project
    /// assumption that the world's present epoch is that long).
    pub epoch_duration_kyr: i64,
    /// Glaciers narrower than this get no section or long-profile carve (not resolved by the 512 m network; cirques only).
    pub min_section_width_mm: i64,
    /// Channelized ice only: a section needs W <= confinement x sqrt(upstream ice-catchment area) (permille).
    pub confinement_permille: i64,
    /// Cirques closer than this to an accepted one (head nearer the ELA first) are dropped.
    pub cirque_spacing_mm: i64,
    pub production_rule: String,
    pub flux_rule: String,
    pub blend_rule: String,
    pub width_rule: String,
    pub thickness_rule: String,
    pub surface_rule: String,
    pub imprint_rule: String,
    pub cirque_rule: String,
    pub section_rule: String,
    /// Lane GLACIER-1: 1 when the t0.glacial layer declares "present ice only" (no glacial-period scenario): no sections,
    /// no carve, no cirques. Absent (0) = GLACIAL-T1 unchanged (and its recipes byte-identical).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub present_ice_only: i64,
    /// Lane GLACIER-1 (key gl_now=1): present-day ice (NOW state, plastic icefield sheet, window fields). None = unchanged.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub now: Option<IceNowParamsV1>,
}

fn is_zero(v: &i64) -> bool { *v == 0 }

/// Lane GLACIER-1: declared constants of the present-ice build (exact canonical JSON; only with gl_now=1).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct IceNowParamsV1 {
    /// blake3 of the layer's ela_interglacial_m.f32 (not part of the GLACIAL-T1 t0.glacial digest, so carried here).
    pub ela_now_digest: String,
    /// Window mask threshold on the bilinear of the sheet mask (permille).
    pub window_mask_permille: i64,
    pub state_rule: String,
    pub mask_rule: String,
    pub sheet_rule: String,
    pub mass_rule: String,
    pub window_rule: String,
    pub class_rule: String,
    /// Lane GLACIER-2 (key gl_bal=1): the balance-flux icefield replaces the GLACIER-1 sheet mask. None = GLACIER-1.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bal: Option<IceBalParamsV1>,
}

/// Lane GLACIER-2 (key gl_bal=1): declared constants of the balance-flux icefield, calving and the 32 m margin.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct IceBalParamsV1 {
    /// Routing / mass-balance passes (pass 1 on the bed, later passes on the previous pass's ice surface).
    pub passes: i64,
    /// Brown, Meier and Post 1982 (USGS PP 1258-C): calving speed V_c = C h_w, C = 27.1 a^-1 (milli per year).
    pub calving_coeff_sea_milli_per_yr: i64,
    /// Freshwater calving: about 10 % of tidewater rates (Funk and Roethlisberger 1989) (milli per year).
    pub calving_coeff_lake_milli_per_yr: i64,
    pub seawater_density_kg_m3: i64,
    /// Calving-face thickness = (1 + q) x flotation thickness (Vieli, Funk and Blatter 2001 form; q a project assumption).
    pub buoyancy_excess_permille: i64,
    pub route_rule: String,
    pub claim_rule: String,
    pub calving_rule: String,
    pub sheet_rule: String,
    pub window_rule: String,
    /// Lane GLACIER-3 (key gl_geo=1, needs gl_bal=1): Euclidean plastic sheet, D-infinity routing, flow-aligned marine
    /// claims, converged passes. None = GLACIER-2 (recipes byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub geo: Option<IceGeoParamsV1>,
}

/// Lane GLACIER-3 (key gl_geo=1): declared constants of the ice geometry build (exact canonical JSON; only with gl_geo=1).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct IceGeoParamsV1 {
    /// Pass limit and the declared convergence tolerances (ppm of the ice posts, against the previous pass): changed posts
    /// that touch no post iced in both passes (patches), and the net change of the ice post count.
    pub max_passes: i64,
    pub patch_change_tol_ppm: i64,
    pub area_change_tol_ppm: i64,
    /// Under-relaxation of the routing surface between passes (permille of the new ice surface; 1000 = none).
    pub relax_permille: i64,
    pub sheet_rule: String,
    pub route_rule: String,
    pub claim_rule: String,
    pub converge_rule: String,
    pub window_rule: String,
}

impl IceGeoParamsV1 {
    pub fn v1() -> Self {
        IceGeoParamsV1 {
            max_passes: 12,
            patch_change_tol_ppm: 2000,
            area_change_tol_ppm: 1000,
            relax_permille: 500,
            sheet_rule: "plastic_sheet_fmm_v1: the perfectly plastic surface (z_s - b) |grad z_s| = h0 (plastic_sheet_v1's law, h0 = tau_y / (rho_i g)) solved by the fast marching method (Sethian 1996 PNAS 93:1591; first-order upwind Rouy and Tourin 1992 update on the 4-neighbour stencil) instead of Dijkstra on 8 neighbours: at a trial post with accepted neighbours a_x = min(x-neighbours), a_y = min(y-neighbours), u solves (u - b)^2 (max(0, u - a_x)^2 + max(0, u - a_y)^2) = (h0 dx)^2 (bisection, 52 halvings, between max(b, min(a)) and the one-sided root), so fronts are Euclidean (no 0/45/90 degree facets, no 8 % chamfer overestimate); margins and seeds as plastic_sheet_seeded_v1".into(),
            route_rule: "balance_route_dinf_v1: balance_route_v1's priority flood and filled surface F, but each post's flux is split by D-infinity (Tarboton 1997 WRR 33:309): of the 8 triangular facets (centre, a cardinal and an adjacent diagonal domain neighbour) the one with the steepest downslope plane gradient (angle clamped to the facet) gives the flow angle r in [0, pi/4] and the shares (1 - r / (pi/4)) to the cardinal and r / (pi/4) to the diagonal neighbour (both strictly lower F, so donors precede receivers in reverse pop order); a facet with one vertex outside the domain offers its other vertex alone; no positive slope: the pop parent, then the lake, then the chart edge (as v1); glacier identity (terminus post, class) follows the larger share; no split at or into a MARINE post (the whole flux goes to the primary receiver): the calving capacity u_c W(Q) H_f and the width law are laws of a glacier's whole flux, and splitting it over k face posts would multiply the capacity by about k^0.7 (measured: tidewater tongues several km long appearing and vanishing in alternate passes); MARINE posts and claimed (ablation-tongue) posts are routed on the bed in every pass, only accumulation posts (positive balance, the icefield) on the relaxed ice surface: a tongue follows its valley; its own isotropic plastic dome would steer the flowline toward a margin (measured: outlet tongues 6 to 9 km wide shifting sideways by their width in alternate passes)".into(),
            claim_rule: "ablation_claim_aligned_v1: ablation_claim_v1, except that (a) the claim grows from the flowline post in order of Euclidean distance (ties by index) instead of breadth-first over 8-neighbours, so the partial last claim at a terminus is a disc (a rounded snout), not a Chebyshev square, (b) the valley confinement is tested on the bed in every pass, bed(post) <= bed(flowline post) + sqrt(h0 W) (in v1 passes >= 2 test the routing surface, i.e. the previous plastic surface, whose dome over a tongue lifts the level each pass and lets tongues creep up their walls), and (c) a MARINE post is admitted to a flowline post's claim only if it is not ahead of it along its flow direction ((p - p_i) . u_i <= 0, u_i the D-infinity flow direction): ice cannot stand on the sea beyond the face that calves its flux, so a tidewater front ends abreast of its last grounded flowline post (no W/2 disc beyond the face)".into(),
            converge_rule: "passes_to_tolerance_v1: routing / mass-balance passes (pass 1 on the bed, later passes on the routing surface s_k = s_(k-1) + 0.5 (z_s - s_(k-1)) on the previous pass's accumulation posts, the bed elsewhere (route_rule): under-relaxation, a numerical device) repeat until a pass k >= 2 against pass k - 1 has (a) changed posts touching no post iced in both passes (patches: a tongue or icefield piece that moved) <= 2000 ppm of its ice posts and (b) a net change of the ice post count <= 1000 ppm, or 12 passes (converged = false); changed posts touching ice present in both passes (the one-post fringe: a terminus or margin post flipping with the fourth digit of its glacier's flux; measured floor about 0.6 % of ice posts) are published per pass, not gated".into(),
            window_rule: "ice_window_v3 (ice_window_v2 with a continuous cliff rule: per sheet post a fill allowance of 200 m, 0 where its surface differs from a finite 8-neighbour by more than 300 m; window ice = clamp(z_bil - z_post, H_bil, H_bil + allowance_bil) on interior cells, min(z_bil - z_post, H_bil + allowance_bil) on margin and 3-corner cells; no per-cell lowest-corner switch, so no thickness steps on 128 m cell edges) + ice_sub_mm (i32, mm): on marine window posts under the ice, the submerged tongue thickness = published height (sea level) minus the sheet bed (sea floor, bilinear of z_s - H over the mask corners), 0 elsewhere".into(),
        }
    }
}

impl IceBalParamsV1 {
    pub fn v1() -> Self {
        IceBalParamsV1 {
            passes: 3,
            calving_coeff_sea_milli_per_yr: 27_100,
            calving_coeff_lake_milli_per_yr: 2_710,
            seawater_density_kg_m3: 1028,
            buoyancy_excess_permille: 150,
            route_rule: "balance_route_v1 (balance flux, Budd and Warner 1996, recalled): on the 128 m solve lattice, routing surface s = despiked pre-glacial bed b (pass 1) or the previous pass's ice surface on its mask and b elsewhere (passes 2..); priority flood (Barnes et al. 2014) from every land post with a non-land or off-lattice 8-neighbour, filled F = max(s, F of the popping post); receiver = steepest descent of F over land 8-neighbours ((F_i - F_j) / d, ties by lower index), else the pop parent (flats, filled pits), else (outlet posts) the adjacent marine post, then lake post, then the chart edge; specific balance a_s = P ice_share(s, ELA_now, Hf - He) (Hergarten 2021 eq 23 / Liebl 2023 eq 11, the NOW state's law; P and ELA_now bilinear on the node grid), m w.e./yr; flux Q (m3 w.e./yr) accumulated receivers-last".into(),
            claim_rule: "ablation_claim_v1: a post with a_s > 0 is ice and adds a_s x (128 m)^2 to its flux; a post with a_s <= 0 and incoming flux Q > 0 is a flowline post: breadth-first from it over land posts within W/2 of it (W = liebl_width_v1 of the incoming ice flux) whose surface s <= s_i + sqrt(h0 W) (valley confinement: centre thickness of a plastic flat-bed cross-profile), every not-yet-ice post with a_s <= 0 becomes ice and charges its ablation -a_s (128 m)^2 to Q (each post charged once, so total ablation = the drawn ablation area's balance); Q <= 0 ends the claim (terminus; the unused ablation of the last post is reported as clamp); Q > 0 passes to the receiver".into(),
            calving_rule: "flotation_calving_v1: marine posts are part of the routing domain with the sea floor as their bed: depth D = max(0, sea - z_eik, sea - T0 bilinear) (z_eik: the construction sea floor, the eikonal surface before the coast stage; the T0 bathymetry bilinear governs offshore where the construction water level sea - F, a flooding device, is shallower); at a marine post receiving flux Q the calving capacity is C = u_c W H_f (ice m3/yr) with u_c = 27.1 a^-1 x D (Brown, Meier and Post 1982), W = liebl_width_v1(Q), face thickness H_f = (1 + q) (rho_w / rho_i) D (flotation plus q = 0.15, the Vieli, Funk and Blatter 2001 form; q a project assumption); Q_ice <= C: the face stands at that post and Q calves; else the post is grounded ice (claims like a land flowline post, marine posts at the sea-level balance) and the flux moves on down the routing surface (steepest descent of the filled sea floor or ice surface; flats and filled enclosed basins lead to their spill), so the terminus is where supply = ablation + calving; lake posts: freshwater calving 2.71 a^-1 (about 10 %, Funk and Roethlisberger 1989) but the lake floor is not constructed, so lake termini calve the arriving flux with the law not evaluated (counted); flux reaching the chart edge (over land or sea) is exported (not closed in the chart, counted)".into(),
            sheet_rule: "plastic_sheet_seeded_v1: plastic_sheet_v1 on the final pass's mask with the bed = b on land and the sea floor sea - D (calving_rule) on marine mask posts; margin seeds z_s = b, except marine margin posts (sea floor + H_f) and land margin posts of a tidewater glacier next to a marine post (max(b, sea + H_f - D))".into(),
            window_rule: "ice_window_v2: a 32 m post is ice only inside a 128 m cell with at least 3 mask corners (3: a concave margin corner or a one-post enclave; the mask corners' bilinear surface and thickness renormalized over their weights, contour rule); all 4 corners interior (H > 0): ice_window_v1's drape rule; a margin corner (H = 0) or 3 corners: the 32 m thickness contour, ice = min(zs_bil - z_post, H_bil + 200 m) where > 0 (no drape: the margin follows where the bilinear ice surface meets the 32 m ground); marine posts under the mask (grounded tidewater tongue): ice_now = zs_bil - sea level (the ice over the published height); lake posts never".into(),
            geo: None,
        }
    }
}

impl IceNowParamsV1 {
    pub fn v1(ela_now_digest: String) -> Self {
        IceNowParamsV1 {
            ela_now_digest,
            window_mask_permille: 500,
            state_rule: "now_state_v1: ELA_now = ela_interglacial_m (bilinear between T0 cell centres at the node), f_P = 1; now_route_v1: accumulation P s(z) over the node's own posts above the ELA (Hergarten eq 23, as CHAR / MAX), ablation P s(z) (Liebl eq 11, mean over the own posts below the ELA) only over the ice-covered area min(a_own, W(Q_in + accumulation) x link length), Q = max(0, ...) clamped continuity down the tree; W = liebl_width_v1; H = min(0.25 W, tau_y / (f rho_i g S_f)) with S_f the FLUVIAL bed slope of the link (the present ice sits on the uncarved bed), 4 fixed-point passes; z_s = node floor + H; never carves".into(),
            mask_rule: "ice_mask_v1: on the 128 m solve lattice, pre-glacial evaluation surface b (evaluate without the carve, the cone included; despiked: a land post more than 500 m above the median of its finite 8 neighbours takes that median): land posts (not MARINE / LAKE) with b >= ELA_now (accumulation area, all of it, routed or not) OR within W/2 of the link polyline (node -> meander midpoint -> receiver; T1a links straight) of a node with Q_now > 0 and a receiver, W linear along the link to the receiver's W_now (0 at an ice-free receiver), confined to the valley: b <= link floor (b at the post nearest each polyline point, linear along the link) + sqrt(h0 W) (centre thickness of a plastic flat-bed cross-profile of width W); ELA_now at a post = bilinear on the regular lod-8 node grid of the T0 bilinear".into(),
            sheet_rule: "plastic_sheet_v1: perfectly plastic ice |grad z_s| = h0 / (z_s - b), h0 = tau_y / (rho_i g) = 100 kPa / (917 x 9.81) = 11.1 m (Reeh 1982; Benn and Hulton 2010 Profiler in 2-D); z_s = b on mask posts with an 8-neighbour outside the mask or off the lattice (the chart edge is a margin, declared); Dijkstra on 8 neighbours (ties by index), exact local update u = ((z_j + b_i) + sqrt((z_j - b_i)^2 + 4 h0 d)) / 2, d = 128 m or 128 sqrt 2 m (path metric overestimates the Euclidean distance by at most 8 %, so H is at most about 4 % high); H = z_s - b".into(),
            mass_rule: "mass_check_v1: per post production P f(b) a (a = 128 m x 128 m, f = ice_share at ELA_now, P bilinear on the node grid); accumulation routed = posts whose owner node (the node at lod-8 round(x / 512 m)) carries Q_now > 0; unrouted = accumulation on mask posts whose owner node carries no ice; tree closure: sum over ice nodes of (Q - inflow) = outflow at ice termini (calving into water or off-chart), relative residual reported".into(),
            window_rule: "ice_window_v1: at a 32 m window post (x, y): m = bilinear of the sheet mask (0/1), z_s = bilinear of the sheet surface over the finite corners, or the lowest mask-corner surface where the mask corners' surfaces span more than 300 m (a 128 m cliff), H = bilinear of the sheet thickness (0 off the mask), ELA = bilinear on the node grid; ice_now = max(0, min(z_s - z_post, H + 200 m), H) (drape: window-only relief is covered by at least the sheet thickness; fill below the sheet capped at H + 200 m where window and 128 m surface disagree, e.g. a cliff in different posts) where m >= 0.5 and the post is not MARINE / LAKE, else 0 (quantized mm); ice_surface = z_post + ice_now; a pure function of the published layer and the post height, so shared edges agree when heights agree".into(),
            class_rule: "ice_class_v1: 0 no ice; 1 accumulation (firn): ice_surface >= ELA_now; 2 ablation (bare ice): ice_surface < ELA_now".into(),
            bal: None,
        }
    }
}

impl GlacialParamsV1 {
    /// Lane GLACIER-1: the params for a bound t0.glacial layer (v1 constants + the layer's scenario + optional NOW build).
    pub fn for_layer(g: &T0Glacial, now: bool) -> R<Self> {
        let mut p = GlacialParamsV1::v1(g.digest.clone());
        if g.present_only { p.present_ice_only = 1; }
        if now {
            if !g.present_only { return Err("gl_now=1 needs a present-ice-only t0.glacial layer (with a carving scenario the ice would sit on the pre-carve surface; stage GS-C not built)".into()); }
            if g.ela_now.is_empty() { return Err("gl_now=1: the t0.glacial layer has no ela_interglacial_m.f32".into()); }
            p.now = Some(IceNowParamsV1::v1(g.ela_now_digest.clone()));
        }
        Ok(p)
    }
    /// Lane GLACIER-2: key gl_bal=1 (needs gl_now=1): the balance-flux icefield.
    pub fn with_bal(mut self, bal: bool) -> R<Self> {
        if bal { match self.now.as_mut() { Some(n) => n.bal = Some(IceBalParamsV1::v1()), None => return Err("gl_bal=1 needs gl_now=1".into()) } }
        Ok(self)
    }
    /// Lane GLACIER-3: key gl_geo=1 (needs gl_bal=1): the ice geometry build.
    pub fn with_geo(mut self, geo: bool) -> R<Self> {
        if geo { match self.now.as_mut().and_then(|n| n.bal.as_mut()) { Some(b) => b.geo = Some(IceGeoParamsV1::v1()), None => return Err("gl_geo=1 needs gl_bal=1".into()) } }
        Ok(self)
    }
    pub fn v1(t0_glacial_digest: String) -> Self {
        GlacialParamsV1 {
            t0_glacial_digest,
            width_ref_mm: 150_000,
            width_ref_flux_m3_per_yr: 10_000,
            width_exponent_permille: 300,
            thickness_width_permille: 250,
            yield_stress_kpa: 100,
            ice_density_kg_m3: 917,
            adverse_slope_permille: 1500,
            u_power_permille: 2000,
            u_reach_permille_of_width: 1500,
            cirque_radius_mm: 390_000,
            cirque_depth_permille_of_length: 370,
            min_surface_slope_ppm: 1000,
            epoch_duration_kyr: 2580,
            min_section_width_mm: 512_000,
            confinement_permille: 300,
            cirque_spacing_mm: 1_000_000,
            production_rule: "hergarten_liebl_production_v1: p_i = P f_P mean s(z) per node own area (512 m x 512 m), mean over the node's 16 own posts of the pre-glacial 128 m solve surface (evaluate without the carve, same tables); s = min((z - He)/(Hf - He), 1) above He (Hergarten 2021 eq 23), -(x + x^2/3), x = (He - z)/(Hf - He) below (Liebl 2023 eq 11); P (genesis precipitation) and He bilinear between T0 cell centres; Hf - He and f_P from the t0.glacial manifest".into(),
            flux_rule: "continuity_clamped_v1: Q_i(n) = max(0, p_i(n) + sum of donors' Q_i) (Hergarten 2021 eq 21; the clamp is the terminus); land nodes route to their network receiver, T1a nodes to their T1a receiver; water nodes end the ice (calving not modelled)".into(),
            blend_rule: "hergarten_blend_cycle_v1: F = mean over the 8 coldhouse phases of gamma^m r + (1 - gamma)^m (r = K_g/K_f and phases from the t0.glacial manifest, m the T0 law's area exponent); gamma = Q_i / precipitation flux; the present epoch is the coldhouse cycle (declared scenario)".into(),
            width_rule: "liebl_width_v1: W = 150 m x (Q_ice / 1e4 m3/yr)^0.3, Q_ice = Q_i x 1000/917 (Liebl 2023 sec 3.1: an ice flux of 1e4 m3/a generates a width of 150 m; alpha 0.3 Hergarten 2021 from Bahr 1997)".into(),
            thickness_rule: "delta_or_nye_plastic_v1: H = min(0.25 W, tau_y / (f rho_i g S_s)), f = Nye 1965 Table IV parabola shape factor on half-width / depth (0.445, 0.646, 0.746, 0.806 at 1..4, linear, clamped), tau_y 100 kPa, rho_i 917, g 9.81, 4 fixed-point passes".into(),
            surface_rule: "concordant_surface_v1: S_f = max((bed(n) - bed(r)) / L, min slope) of the fluvial beds, S_s = S_f / F; terminus (no ice at the receiver, or a T1a receiver... see below): z_s = fluvial bed + H; upstream z_s(n) = z_s(r) + L S_s; glacial bed = z_s - H, never below bed(r) - psi S_s L (Alley et al. 2003, psi 1.5); T1a nodes keep their level (z_s = level + H)".into(),
            imprint_rule: "time_capped_imprint_v1: carve = min(max over CHAR and MAX of max(0, fluvial bed - glacial bed), (F - 1) E T): the steady glacial profile is reached only as far as the excess Hergarten erosion (F - 1) x the node's steady rate E (the E table, m/My) over the epoch T allows (the excess rate at the pre-glacial slope is the largest it gets: an upper bound); carve <= H of that state (the concordant ice surface never falls below the pre-glacial floor: a trough deepens by at most the ice that fills it); only where W >= 512 m (one network spacing) and the ice is channelized, W <= 0.3 sqrt(upstream area) (a node's own ice spreads as a sheet over its 512 m cell; Aletsch: W 1.5 km over ~120 km2, ratio 0.14; project assumption); section (W, H x I) of the state with the larger carve, I = fraction of the 8 cycle phases with Q_i > 0; land channel beds = tree-filled glacial bed (max with the receiver's), so beds never rise along flow".into(),
            cirque_rule: "ela_cirque_v1: at heads of CHAR glaciers (ice at the node, none in any donor; land nodes) whose floor stands 0 to 2 D_c above the cycle-mean ELA: floor z_c = max(bed - D_c, min(bed, ELA_char), lip - psi S_s R), lip = receiver's filled bed, R = 390 m, D_c = 0.37 x 2R (Barr et al. 2019 morphometry; Porter 1989 floors at the characteristic ELA); armchair bowl (u along the head -> receiver direction, v across): q = (v/R)^2 + (min(u, 0)/R)^2 + (max(u, 0)/2R)^2, z = z_c + D_c q / (1 - q/4) for q < 4 (steep headwall up-glacier, gentle opening toward the lip, continuous at the footprint edge), on free posts; channel posts min(z, max(bowl, max(z_c, lip))); cirques nearer than 1 km to an accepted one are dropped (heads nearest the ELA first)".into(),
            section_rule: "u_section_v1: free and floor posts (bedrock floor) and network channel posts (filled floor) within 1.5 W of a glaciated link (node -> meander midpoint -> receiver) take min(z, floor + H rho^2 / (1 - rho^2/9)), rho = 2|y|/W (continuous at the 1.5 W reach edge), floor/H/W linear along the link (receiver values; zero width at a non-glaciated receiver); after crest rounding and the ct remap, before fans, bank rule and fill".into(),
            present_ice_only: 0,
            now: None,
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.width_ref_mm <= 0 || self.width_ref_flux_m3_per_yr <= 0 || self.thickness_width_permille <= 0 || self.yield_stress_kpa <= 0 || self.ice_density_kg_m3 <= 0 {
            return Err("glacial: width, flux, thickness, stress and density constants must be > 0".into());
        }
        if self.u_power_permille < 1000 || self.u_reach_permille_of_width < 500 || self.cirque_radius_mm <= 0 || self.min_surface_slope_ppm <= 0 {
            return Err("glacial: u power >= 1, reach >= 0.5 W, cirque radius and min slope > 0".into());
        }
        Ok(())
    }
}

/// The t0.glacial layer as T1b reads it.
pub struct T0Glacial {
    pub digest: String,
    pub ratio: f64,
    pub full_m: f64,
    pub phase_fp: Vec<f64>,
    pub char_fp: f64,
    pub max_fp: f64,
    pub ela_phase: Vec<Vec<f32>>,
    pub ela_char: Vec<f32>,
    pub ela_max: Vec<f32>,
    pub precip: Vec<f32>,
    /// Lane GLACIER-1: the manifest's scenario.present_ice_only (false when absent), the interglacial (present) ELA when the
    /// layer has it (empty otherwise) and its file digest (outside `digest`, which stays the GLACIAL-T1 one).
    pub present_only: bool,
    pub ela_now: Vec<f32>,
    pub ela_now_digest: String,
}

fn read_f32(p: &Path, n: usize) -> R<(Vec<f32>, blake3::Hash)> {
    let b = std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?;
    if b.len() != 4 * n { return Err(format!("{}: {} bytes, expected {}", p.display(), b.len(), 4 * n)); }
    let v: Vec<f32> = b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
    if v.iter().any(|x| !x.is_finite()) { return Err(format!("{}: nonfinite value", p.display())); }
    Ok((v, blake3::hash(&b)))
}

impl T0Glacial {
    /// Load and bind to T0 (the manifest must cite the same terrain artifact digest).
    pub fn load(dir: &Path, t0: &T0) -> R<T0Glacial> {
        let mb = std::fs::read(dir.join("glacial.json")).map_err(|e| format!("t0.glacial manifest: {e}"))?;
        let m: Value = serde_json::from_slice(&mb).map_err(|e| e.to_string())?;
        if m["schema"] != GL_SCHEMA_T0 { return Err(format!("t0.glacial schema {} (need {GL_SCHEMA_T0})", m["schema"])); }
        if m["artifact_digests"]["terrain"].as_str() != Some(t0.terrain_digest.as_str()) {
            return Err(format!("t0.glacial cites terrain {} but T0 is {}; refusing", m["artifact_digests"]["terrain"], t0.terrain_digest));
        }
        if m["cycle_phases"].as_u64() != Some(GL_PHASES as u64) { return Err("t0.glacial: 8 cycle phases expected".into()); }
        let n = t0.nlat * t0.nlon;
        let mut h = blake3::Hasher::new();
        h.update(&mb);
        let mut ld = |name: &str| -> R<Vec<f32>> { let (v, d) = read_f32(&dir.join(format!("{name}.f32")), n)?; h.update(name.as_bytes()); h.update(d.as_bytes()); Ok(v) };
        let mut ela_phase = Vec::new();
        for k in 0..GL_PHASES { ela_phase.push(ld(&format!("ela_cycle_phase{k}_m"))?); }
        let ela_char = ld("ela_cycle_mean_m")?;
        let ela_max = ld("ela_glacial_max_m")?;
        let precip = ld("precipitation_mm_per_yr")?;
        let digest = format!("blake3:{}", h.finalize().to_hex());
        let st = m["states"].as_array().ok_or("t0.glacial states")?;
        if st.len() < GL_PHASES { return Err("t0.glacial: fewer states than cycle phases".into()); }
        let phase_fp: Vec<f64> = st[..GL_PHASES].iter().map(|s| s["precipitation_factor"].as_f64().unwrap_or(f64::NAN)).collect();
        if phase_fp.iter().any(|v| !v.is_finite() || *v <= 0.0) { return Err("t0.glacial: bad phase precipitation factors".into()); }
        let b = &m["block"];
        let ratio = b["glacial_to_fluvial_erodibility_ratio"].as_f64().ok_or("block ratio")?;
        let full_m = b["full_ice_height_above_ela_m"].as_f64().ok_or("block full ice height")?;
        let max_fp = b["glacial_max_precipitation_factor"].as_f64().ok_or("block glacial max precipitation factor")?;
        if !(ratio > 0.0 && full_m > 0.0 && max_fp > 0.0) { return Err("t0.glacial: block constants must be > 0".into()); }
        let char_fp = phase_fp.iter().sum::<f64>() / GL_PHASES as f64;
        // Lane GLACIER-1: scenario flag and the present ELA (read after the digest, which is unchanged).
        let present_only = m["scenario"]["present_ice_only"].as_bool().unwrap_or(false);
        let pn = dir.join("ela_interglacial_m.f32");
        let (ela_now, ela_now_digest) = if pn.exists() { let (v, d) = read_f32(&pn, n)?; (v, format!("blake3:{}", d.to_hex())) } else { (Vec::new(), String::new()) };
        if present_only {
            // The declared scenario must be what the fields say: every phase, the cycle mean and the maximum = the present ELA.
            if ela_now.is_empty() { return Err("t0.glacial: present_ice_only without ela_interglacial_m".into()); }
            let same = |v: &[f32]| v.iter().zip(&ela_now).all(|(a, b)| a.to_bits() == b.to_bits());
            if !(ela_phase.iter().all(|v| same(v)) && same(&ela_char) && same(&ela_max)) { return Err("t0.glacial: present_ice_only but phase / cycle-mean / maximum ELAs differ from the interglacial ELA; refusing".into()); }
        }
        Ok(T0Glacial { digest, ratio, full_m, phase_fp, char_fp, max_fp, ela_phase, ela_char, ela_max, precip, present_only, ela_now, ela_now_digest })
    }
}

/// Hergarten 2021 eq 23 / Liebl 2023 eq 11 production share (the T0 law's `ice_share`, same expression).
pub fn ice_share(h: f64, ela: f64, full: f64) -> f64 {
    let x = (h - ela) / full;
    if x >= 0.0 { x.min(1.0) } else { let y = -x; -(y + y * y / 3.0) }
}

/// Nye 1965 Table IV, parabolic channel, no slip: shape factor at half-width / depth (linear between 1..4, clamped).
pub fn nye_parabola(half_w_over_h: f64) -> f64 {
    const T: [(f64, f64); 4] = [(1.0, 0.445), (2.0, 0.646), (3.0, 0.746), (4.0, 0.806)];
    let w = half_w_over_h;
    if !(w > T[0].0) { return T[0].1; }
    for k in 1..4 { if w <= T[k].0 { let (a, b) = (T[k - 1], T[k]); return a.1 + (b.1 - a.1) * (w - a.0) / (b.0 - a.0); } }
    T[3].1
}

/// One cirque.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Cirque { pub node: u32, pub x_m: f64, pub y_m: f64, pub ux: f64, pub uy: f64, pub floor_m: f64, pub flat_m: f64, pub radius_m: f64, pub depth_m: f64, pub head_bed_m: f64, pub ela_char_m: f64, pub lip_m: f64 }

/// The published chart table.
#[derive(Clone, Debug, Default)]
pub struct GlacialTable {
    /// Ice routing receiver per node (NONE where the ice ends).
    pub rcv: Vec<u32>,
    pub q_char: Vec<f32>,
    pub q_max: Vec<f32>,
    pub f_cycle: Vec<f32>,
    pub occupation: Vec<f32>,
    /// Section width and depth (0 where no ice), final bedrock bed, filled bed (channels), carve depth.
    pub w: Vec<f32>,
    pub h: Vec<f32>,
    pub bed: Vec<f32>,
    pub bed_fill: Vec<f32>,
    pub carve: Vec<f32>,
    pub zs_char: Vec<f32>,
    pub zs_max: Vec<f32>,
    pub ela_char: Vec<f32>,
    /// Ice-surface slope S_s = S_f / F per node (for the published adverse-slope check).
    pub s_surface: Vec<f32>,
    pub cirques: Vec<Cirque>,
    pub stats: Value,
    /// Lane GLACIER-1 (gl_now=1 only): the present-ice NOW state and the icefield sheet.
    pub now: Option<IceNow>,
}

/// Lane GLACIER-1: NOW state per lod-8 node and the chart icefield sheet on the 128 m solve lattice.
#[derive(Clone, Debug, Default)]
pub struct IceNow {
    pub q: Vec<f32>,
    pub w: Vec<f32>,
    pub h: Vec<f32>,
    pub zs: Vec<f32>,
    /// ELA_now and P (mm/yr) on the REGULAR lod-8 node grid (node n at (gx, gy) x 512 m): interpolation grid for posts.
    pub ela_grid: Vec<f32>,
    pub p_grid: Vec<f32>,
    /// Sheet sub-grid (bounding box of the mask + 2 posts) on the solve lattice: origin, size, spacing; surface z_s (b off
    /// the mask, NaN off the lattice's land) and thickness H (>= 0 on the mask, margin posts 0; -1 off the mask).
    pub dx_mm: i64,
    pub gx0: i64,
    pub gy0: i64,
    pub nx: usize,
    pub ny: usize,
    pub sheet_zs: Vec<f32>,
    pub sheet_h: Vec<f32>,
    pub stats: Value,
    /// Lane GLACIER-2 (gl_bal=1 only; empty otherwise): balance flux (m3 w.e./yr, 0 off the ice) and the terminus class of
    /// the glacier each sheet post belongs to (TERM_* codes) on the same sub-grid.
    pub sheet_q: Vec<f32>,
    pub sheet_term: Vec<u8>,
}

/// Lane GLACIER-1: names of the NOW per-node fields and the sheet fields (layer blobs).
pub const NOW_FIELDS: [&str; 6] = ["gl_q_now_m3_yr", "gl_w_now_m", "gl_h_now_m", "gl_zs_now_m", "gl_ela_now_grid_m", "gl_p_now_grid_mm_yr"];
pub const SHEET_FIELDS: [&str; 2] = ["ice_sheet_zs_m", "ice_sheet_h_m"];
/// Lane GLACIER-2 (gl_bal=1 only): balance flux and terminus class on the sheet sub-grid.
pub const SHEET_BAL_FIELDS: [&str; 2] = ["ice_sheet_q_m3_yr", "ice_sheet_term"];
/// Lane GLACIER-2 terminus classes (ice_sheet_term.u8): 0 no ice, 1 land-terminating, 2 tidewater (marine calving),
/// 3 lake calving, 4 leaves the chart (export, not closed in the chart).
pub const TERM_NONE: u8 = 0;
pub const TERM_LAND: u8 = 1;
pub const TERM_SEA: u8 = 2;
pub const TERM_LAKE: u8 = 3;
pub const TERM_EDGE: u8 = 4;

/// Bilinear on the regular lod-8 node grid (clamped to the chart), at chart position (x, y) m.
pub fn node_grid_bilinear(v: &[f32], x_m: f64, y_m: f64) -> f64 {
    let s = crate::SPACING_MM as f64 / 1000.0;
    let fx = (x_m / s - crate::G_MIN as f64).clamp(0.0, (crate::N - 1) as f64);
    let fy = (y_m / s - crate::G_MIN as f64).clamp(0.0, (crate::N - 1) as f64);
    let (x0, y0) = ((fx.floor() as usize).min(crate::N - 2), (fy.floor() as usize).min(crate::N - 2));
    let (tx, ty) = (fx - x0 as f64, fy - y0 as f64);
    let at = |x: usize, y: usize| v[y * crate::N + x] as f64;
    (1.0 - ty) * ((1.0 - tx) * at(x0, y0) + tx * at(x0 + 1, y0)) + ty * ((1.0 - tx) * at(x0, y0 + 1) + tx * at(x0 + 1, y0 + 1))
}

impl IceNow {
    /// Sheet reads at chart position (x, y) m: (mask bilinear 0..1, sheet surface (bilinear; the lowest mask-corner surface at
    /// a 128 m cliff), ELA_now, thickness bilinear with 0 off the mask).
    pub fn sample(&self, x_m: f64, y_m: f64) -> (f64, f64, f64, f64) {
        let ela = node_grid_bilinear(&self.ela_grid, x_m, y_m);
        let d = self.dx_mm as f64 / 1000.0;
        let (fx, fy) = (x_m / d - self.gx0 as f64, y_m / d - self.gy0 as f64);
        if self.nx < 2 || self.ny < 2 || !(fx >= 0.0 && fy >= 0.0 && fx <= (self.nx - 1) as f64 && fy <= (self.ny - 1) as f64) { return (0.0, f64::NAN, ela, 0.0); }
        let (x0, y0) = ((fx.floor() as usize).min(self.nx - 2), (fy.floor() as usize).min(self.ny - 2));
        let (tx, ty) = (fx - x0 as f64, fy - y0 as f64);
        let (mut m, mut zmin, mut zmax, mut z, mut wz, mut hh) = (0.0, f64::INFINITY, f64::NEG_INFINITY, 0.0, 0.0, 0.0);
        for (dy, wy) in [(0usize, 1.0 - ty), (1, ty)] { for (dx, wx) in [(0usize, 1.0 - tx), (1, tx)] {
            let k = (y0 + dy) * self.nx + x0 + dx;
            let w = wx * wy;
            if self.sheet_h[k] >= 0.0 { m += w; hh += w * self.sheet_h[k] as f64; if self.sheet_zs[k].is_finite() { zmin = zmin.min(self.sheet_zs[k] as f64); zmax = zmax.max(self.sheet_zs[k] as f64); } }
            let zz = self.sheet_zs[k];
            if zz.is_finite() { z += w * zz as f64; wz += w; }
        } }
        // Bilinear surface over the finite corners; at a 128 m cliff (mask-corner surfaces more than 300 m apart, declared) the
        // lowest mask-corner surface instead (a bilinear would bridge the step with ice).
        let zs = if !zmin.is_finite() || wz <= 0.0 { f64::NAN } else if zmax - zmin > 300.0 { zmin } else { z / wz };
        (m, zs, ela, hh)
    }
}

/// Lane GLACIER-1: per window post (chart position, quantized height mm, water flag) the published ice fields:
/// (ice_now_mm, ice_surface_mm, ice_class).
pub fn window_ice(now: &IceNow, np: &IceNowParamsV1, nn: usize, xy: &dyn Fn(usize) -> (f64, f64), height_mm: &[i32], water: &dyn Fn(usize) -> bool) -> (Vec<i32>, Vec<i32>, Vec<u8>) {
    let thr = np.window_mask_permille as f64 / 1000.0;
    let (mut ice, mut surf, mut cls) = (vec![0i32; nn], height_mm.to_vec(), vec![0u8; nn]);
    for i in 0..nn {
        if water(i) { continue; }
        let (x, y) = xy(i);
        let (m, zs, ela, hs) = now.sample(x, y);
        if !(m >= thr) || !zs.is_finite() { continue; }
        let zp = height_mm[i] as f64 / 1000.0;
        // Fill up to the lowest corner surface (window channels below the 128 m bed fill with ice), but never thinner than
        // the sheet thickness (drape: window-only relief, e.g. VE units or 32 m ridges, is covered, not poked through).
        // Fill below the 128 m surface is capped at H + 200 m (declared): where the 32 m window and the 128 m evaluation put
        // a cliff in different posts the bilinear surface would otherwise stand kilometres above the window ground.
        let h = (zs - zp).min(hs + 200.0).max(hs);
        if !(h > 0.0) { continue; }
        let hm = vz_tilekey::chart::quantize_mm(h * 1000.0).clamp(0, i32::MAX as i64) as i32;
        if hm <= 0 { continue; }
        ice[i] = hm;
        surf[i] = height_mm[i].saturating_add(hm);
        cls[i] = if surf[i] as f64 / 1000.0 >= ela { 1 } else { 2 };
    }
    (ice, surf, cls)
}

struct StateOut { q: Vec<f64>, w: Vec<f64>, h: Vec<f64>, zs: Vec<f64>, bed: Vec<f64> }

/// The pre-glacial surface of each node's own 4 x 4 posts of a 128 m lattice (NaN off-lattice, water or unset): the
/// node at lod-8 (gx, gy) owns posts 4 gx - 2 .. 4 gx + 1 (and y alike).
pub fn node_surface(lat: &EvalLat, z: &[f64], fixed: &[u8]) -> Vec<[f32; 16]> {
    let f = 512_000 / lat.dx_mm;
    let mut out = vec![[f32::NAN; 16]; NODES];
    if f != 4 { return out; }
    for (n, o) in out.iter_mut().enumerate() {
        let (gx, gy) = crate::node_g(n);
        let mut k = 0;
        for oy in -2..2i64 { for ox in -2..2i64 {
            if let Some(i) = lat.idx(gx * 4 + ox, gy * 4 + oy) { if fixed[i] != crate::t1b::MARINE && fixed[i] != crate::t1b::LAKE && z[i].is_finite() { o[k] = z[i] as f32; } }
            k += 1;
        } }
    }
    out
}

/// Build the chart table from a network, its 0.9.x node beds (`bed_in`) and channel flags, the pre-glacial 128 m surface
/// of each node (`zpost`, from `node_surface`), node lat/lon and the T0 layer. The glacial geometry runs on the node
/// FLOOR: the bed on channel and T1a nodes; min(bed, lowest own post) on other land nodes (their bed-law value is not
/// drawn and can stand above the surface). Production uses the 16 own posts (each 128 m x 128 m).
#[allow(clippy::too_many_arguments)]
pub fn build_table(net: &Network, bed_in: &[f64], chan: &[bool], zpost: &[[f32; 16]], e_node: &[f64], latlon: &[(f64, f64)], t0: &T0, g: &T0Glacial, m: f64, gp: &GlacialParamsV1) -> R<GlacialTable> {
    let t_start = std::time::Instant::now();
    let bed_f: Vec<f64> = (0..NODES).map(|n| {
        if net.kind[n] != K_LAND || chan[n] { return bed_in[n]; }
        let lo = zpost[n].iter().filter(|v| v.is_finite()).fold(f64::INFINITY, |a, v| a.min(*v as f64));
        bed_in[n].min(lo)
    }).collect();
    let bed_f = &bed_f[..];
    let ok = |r: u32| r != NONE && r != OFF_CHART && (r as usize) < NODES;
    let ice_kind = |n: usize| net.kind[n] == K_LAND || net.kind[n] == K_T1A;
    let mut rcv = vec![NONE; NODES];
    for n in 0..NODES {
        let r = match net.kind[n] { K_LAND => net.recv[n], K_T1A => net.t1a_recv[n], _ => NONE };
        if ok(r) { rcv[n] = r; }
    }
    // Donors-first order (Kahn; queue in index order: deterministic).
    let mut indeg = vec![0u32; NODES];
    for n in 0..NODES { if rcv[n] != NONE { indeg[rcv[n] as usize] += 1; } }
    let mut order: Vec<u32> = (0..NODES as u32).filter(|&n| indeg[n as usize] == 0).collect();
    let mut k = 0;
    while k < order.len() {
        let n = order[k] as usize; k += 1;
        if rcv[n] != NONE { let r = rcv[n] as usize; indeg[r] -= 1; if indeg[r] == 0 { order.push(r as u32); } }
    }
    if order.len() != NODES { return Err(format!("glacial: ice routing has a cycle ({} of {} nodes ordered)", order.len(), NODES)); }
    // Link length (land: through the meander midpoint, as the bed rule) and fluvial slope.
    let len: Vec<f64> = (0..NODES).map(|n| {
        if rcv[n] == NONE { return 0.0; }
        let r = rcv[n] as usize;
        let (a, b) = (net.pos[n], net.pos[r]);
        if net.kind[n] == K_LAND { let mm = net.mid[n]; (a[0] - mm[0]).hypot(a[1] - mm[1]) + (mm[0] - b[0]).hypot(mm[1] - b[1]) } else { (a[0] - b[0]).hypot(a[1] - b[1]) }
    }).collect();
    let smin = gp.min_surface_slope_ppm as f64 / 1e6;
    let s_f: Vec<f64> = (0..NODES).map(|n| if rcv[n] == NONE || len[n] <= 0.0 { smin } else { ((bed_f[n] - bed_f[rcv[n] as usize]) / len[n]).max(smin) }).collect();
    // Climate at the nodes.
    let bil = |f: &[f32], n: usize| -> f64 { let (la, lo) = latlon[n]; crate::t1b::bilinear(t0, &|c: u32| f[c as usize] as f64, la, lo, false) };
    let p_m: Vec<f64> = (0..NODES).map(|n| if ice_kind(n) { (bil(&g.precip, n) / 1000.0).max(0.0) } else { 0.0 }).collect();
    let a_own = (crate::SPACING_MM as f64 / 1000.0).powi(2);
    let mut area_up = vec![0.0f64; NODES];
    for &n in &order {
        let n = n as usize;
        if !ice_kind(n) { continue; }
        area_up[n] += a_own;
        if rcv[n] != NONE { let r = rcv[n] as usize; if ice_kind(r) { area_up[r] += area_up[n]; } }
    }
    let route = |ela: &dyn Fn(usize) -> f64, fp: f64| -> (Vec<f64>, Vec<f64>) {
        let mut q = vec![0.0f64; NODES];
        let mut pf = vec![0.0f64; NODES];
        for &n in &order {
            let n = n as usize;
            if !ice_kind(n) { q[n] = 0.0; pf[n] = 0.0; continue; }
            let (mut sh, mut np) = (0.0f64, 0u32);
            for v in zpost[n].iter() { if v.is_finite() { sh += ice_share(*v as f64, ela(n), g.full_m); np += 1; } }
            let share = if np > 0 { sh / np as f64 } else { ice_share(bed_f[n], ela(n), g.full_m) };
            let prod = a_own * p_m[n] * fp * share;
            q[n] = (q[n] + prod).max(0.0);
            pf[n] += a_own * p_m[n] * fp;
            if rcv[n] != NONE { let r = rcv[n] as usize; if ice_kind(r) { q[r] += q[n]; pf[r] += pf[n]; } }
        }
        (q, pf)
    };
    // The 8 phases: occupation and the Hergarten blend.
    let mut f_sum = vec![0.0f64; NODES];
    let mut occ = vec![0.0f64; NODES];
    for k in 0..GL_PHASES {
        let e: Vec<f64> = (0..NODES).map(|n| if ice_kind(n) { bil(&g.ela_phase[k], n) } else { 0.0 }).collect();
        let (q, pf) = route(&|n| e[n], g.phase_fp[k]);
        for n in 0..NODES {
            let gam = if pf[n] > 0.0 { (q[n] / pf[n]).clamp(0.0, 1.0) } else { 0.0 };
            f_sum[n] += if gam > 0.0 { gam.powf(m) * g.ratio + (1.0 - gam).powf(m) } else { 1.0 };
            if q[n] > 0.0 { occ[n] += 1.0 / GL_PHASES as f64; }
        }
    }
    let f_cyc: Vec<f64> = f_sum.iter().map(|v| v / GL_PHASES as f64).collect();
    let s_s: Vec<f64> = (0..NODES).map(|n| (s_f[n] / f_cyc[n]).max(smin)).collect();
    let ela_char: Vec<f64> = (0..NODES).map(|n| if ice_kind(n) { bil(&g.ela_char, n) } else { 0.0 }).collect();
    let ela_max: Vec<f64> = (0..NODES).map(|n| if ice_kind(n) { bil(&g.ela_max, n) } else { 0.0 }).collect();
    let (wr, qr, al) = (gp.width_ref_mm as f64 / 1000.0, gp.width_ref_flux_m3_per_yr as f64, gp.width_exponent_permille as f64 / 1000.0);
    let delta = gp.thickness_width_permille as f64 / 1000.0;
    let rho = gp.ice_density_kg_m3 as f64;
    let tau = gp.yield_stress_kpa as f64 * 1000.0;
    let psi = gp.adverse_slope_permille as f64 / 1000.0;
    let mut nye_limited = 0u64;
    let mut state = |ela: &[f64], fp: f64| -> StateOut {
        let (q, _) = route(&|n| ela[n], fp);
        let mut w = vec![0.0f64; NODES];
        let mut h = vec![0.0f64; NODES];
        for n in 0..NODES {
            if q[n] <= 0.0 { continue; }
            let qi = q[n] * 1000.0 / rho;
            w[n] = wr * (qi / qr).powf(al);
            let mut hh = delta * w[n];
            for _ in 0..4 {
                let f = nye_parabola(0.5 * w[n] / hh.max(1e-6));
                let hp = tau / (f * rho * 9.81 * s_s[n]);
                hh = (delta * w[n]).min(hp);
            }
            if hh < delta * w[n] - 1e-9 { nye_limited += 1; }
            h[n] = hh;
        }
        let mut zs = vec![f64::NAN; NODES];
        let mut bed = bed_f.to_vec();
        // Receivers first.
        for &n in order.iter().rev() {
            let n = n as usize;
            if q[n] <= 0.0 || !ice_kind(n) { continue; }
            let r = rcv[n];
            if net.kind[n] == K_T1A { zs[n] = bed_f[n] + h[n]; continue; }
            let glacial_r = r != NONE && q[r as usize] > 0.0 && zs[r as usize].is_finite();
            if !glacial_r { zs[n] = bed_f[n] + h[n]; continue; }
            let r = r as usize;
            zs[n] = zs[r] + len[n] * s_s[n];
            let b = zs[n] - h[n];
            bed[n] = b.max(bed[r] - psi * s_s[n] * len[n]);
        }
        StateOut { q, w, h, zs, bed }
    };
    let sc = state(&ela_char, g.char_fp);
    let sm = state(&ela_max, g.max_fp);
    // Lane GLACIER-1 (gl_now=1): the NOW state (present ELA, f_P 1) on the uncarved bed: plastic limit at the fluvial slope.
    let present = gp.present_ice_only == 1;
    let now = if gp.now.is_some() {
        if !present || g.ela_now.is_empty() { return Err("glacial: the NOW state needs a present-ice-only layer with ela_interglacial_m".into()); }
        let ela_n: Vec<f64> = (0..NODES).map(|n| if ice_kind(n) { bil(&g.ela_now, n) } else { 0.0 }).collect();
        // now_route_v1: accumulation over the node's own posts above the ELA (as CHAR / MAX), ablation only over the ice-covered
        // area min(a_own, W(Q) L) at the mean melt share of the own posts below the ELA (the ice is a tongue W wide, not the
        // whole 512 m cell); W from the flux entering the node plus its accumulation.
        let mut q = vec![0.0f64; NODES];
        for &n in &order {
            let n = n as usize;
            if !ice_kind(n) { q[n] = 0.0; continue; }
            let (mut acc_s, mut abl_s, mut np) = (0.0f64, 0.0f64, 0u32);
            for v in zpost[n].iter() { if v.is_finite() { let sh = ice_share(*v as f64, ela_n[n], g.full_m); if sh > 0.0 { acc_s += sh; } else { abl_s += sh; } np += 1; } }
            if np == 0 { let sh = ice_share(bed_f[n], ela_n[n], g.full_m); if sh > 0.0 { acc_s = sh; } else { abl_s = sh; } np = 1; }
            let qin = q[n] + a_own * p_m[n] * acc_s / np as f64;
            let a_ice = if qin > 0.0 { let lk = if len[n] > 0.0 { len[n] } else { crate::SPACING_MM as f64 / 1000.0 }; (wr * (qin * 1000.0 / rho / qr).powf(al) * lk).min(a_own) } else { 0.0 };
            q[n] = (qin + a_ice * p_m[n] * abl_s / np as f64).max(0.0);
            if rcv[n] != NONE { let r = rcv[n] as usize; if ice_kind(r) { q[r] += q[n]; } }
        }
        let (mut wv, mut hv, mut zv) = (vec![0.0f64; NODES], vec![0.0f64; NODES], vec![0.0f64; NODES]);
        for n in 0..NODES {
            if q[n] <= 0.0 { continue; }
            let qi = q[n] * 1000.0 / rho;
            wv[n] = wr * (qi / qr).powf(al);
            let mut hh = delta * wv[n];
            for _ in 0..4 { let f = nye_parabola(0.5 * wv[n] / hh.max(1e-6)); hh = (delta * wv[n]).min(tau / (f * rho * 9.81 * s_f[n])); }
            hv[n] = hh;
            zv[n] = bed_f[n] + hh;
        }
        // Tree closure (mass_check_v1): sum of (Q - inflow) over ice nodes = outflow at ice termini.
        let mut inflow = vec![0.0f64; NODES];
        for n in 0..NODES { if q[n] > 0.0 && rcv[n] != NONE && ice_kind(rcv[n] as usize) { inflow[rcv[n] as usize] += q[n]; } }
        let (mut acc, mut abl, mut out_w, mut out_off) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for n in 0..NODES {
            if !ice_kind(n) { continue; }
            let d = q[n] - inflow[n];
            if d > 0.0 { acc += d; } else { abl += d; }
            if q[n] > 0.0 && (rcv[n] == NONE || !ice_kind(rcv[n] as usize)) { if rcv[n] == NONE { out_off += q[n]; } else { out_w += q[n]; } }
        }
        let res = acc + abl - out_w - out_off;
        let nq = q.iter().filter(|v| **v > 0.0).count();
        let st = json!({"nodes_ice_now": nq, "q_now_max_m3_yr": q.iter().cloned().fold(0.0, f64::max), "w_now_max_m": wv.iter().cloned().fold(0.0, f64::max), "h_now_max_m": hv.iter().cloned().fold(0.0, f64::max),
            "tree_net_accumulation_m3_yr": acc, "tree_net_ablation_m3_yr": abl, "tree_outflow_to_water_m3_yr": out_w, "tree_outflow_no_receiver_m3_yr": out_off, "tree_residual_m3_yr": res, "tree_relative_residual": if acc > 0.0 { res.abs() / acc } else { 0.0 }});
        Some(IceNow { q: q.iter().map(|v| *v as f32).collect(), w: wv.iter().map(|v| *v as f32).collect(), h: hv.iter().map(|v| *v as f32).collect(), zs: zv.iter().map(|v| *v as f32).collect(), stats: st, ..Default::default() })
    } else { None };
    let mut t = GlacialTable { rcv: rcv.clone(), ..Default::default() };
    let f32v = |v: &[f64]| v.iter().map(|x| *x as f32).collect::<Vec<f32>>();
    let mut w = vec![0.0f64; NODES];
    let mut h = vec![0.0f64; NODES];
    let mut bed = bed_f.to_vec();
    let mut carve = vec![0.0f64; NODES];
    let (mut n_carved, mut carve_max, mut n_section, mut n_capped, mut n_narrow) = (0u64, 0.0f64, 0u64, 0u64, 0u64);
    let t_my = gp.epoch_duration_kyr as f64 / 1000.0;
    let wmin = gp.min_section_width_mm as f64 / 1000.0;
    for n in 0..NODES {
        if !ice_kind(n) || occ[n] <= 0.0 || present { continue; } // Lane GLACIER-1: present ice only = no sections
        let cc = (bed_f[n] - sc.bed[n]).max(0.0);
        let cm = (bed_f[n] - sm.bed[n]).max(0.0);
        let (wn, hn, c) = if sc.q[n] > 0.0 && (cc >= cm || sm.q[n] <= 0.0) { (sc.w[n], sc.h[n], cc) } else if sm.q[n] > 0.0 { (sm.w[n], sm.h[n], cm) } else { (0.0, 0.0, 0.0) };
        if wn <= 0.0 { continue; }
        if wn < wmin || wn > gp.confinement_permille as f64 / 1000.0 * area_up[n].sqrt() { n_narrow += 1; continue; }
        n_section += 1;
        w[n] = wn;
        h[n] = occ[n] * hn;
        if net.kind[n] == K_LAND {
            let cap = (f_cyc[n] - 1.0).max(0.0) * e_node[n].max(0.0) * t_my;
            if c > cap { n_capped += 1; }
            carve[n] = c.min(cap).min(hn);
            bed[n] = bed_f[n] - carve[n];
            if carve[n] > 0.0 { n_carved += 1; carve_max = carve_max.max(carve[n]); }
        }
    }
    // Alley et al. 2003 on the final (time-capped) bed too: walking upstream, a sectioned land node's bed never falls below
    // its sectioned receiver's by more than psi S_s L (the cap varies with E between T0 cells).
    let mut n_alley = 0u64;
    for &n in order.iter().rev() {
        let n = n as usize;
        if net.kind[n] != K_LAND || w[n] <= 0.0 || rcv[n] == NONE { continue; }
        let r = rcv[n] as usize;
        if w[r] <= 0.0 { continue; }
        let lim = bed[r] - psi * s_s[n] * len[n];
        if bed[n] < lim { bed[n] = lim; carve[n] = bed_f[n] - lim; n_alley += 1; }
    }
    // Tree fill (receivers first): a land node never stands below its receiver's filled bed.
    let mut fill = bed.clone();
    for &n in order.iter().rev() {
        let n = n as usize;
        if net.kind[n] != K_LAND || rcv[n] == NONE { continue; }
        let r = rcv[n] as usize;
        if fill[n] < fill[r] { fill[n] = fill[r]; }
    }
    // Cirques at CHAR heads.
    let rr = gp.cirque_radius_mm as f64 / 1000.0;
    let dc = gp.cirque_depth_permille_of_length as f64 / 1000.0 * 2.0 * rr;
    let mut has_donor = vec![false; NODES];
    for n in 0..NODES { if sc.q[n] > 0.0 && rcv[n] != NONE { has_donor[rcv[n] as usize] = true; } }
    let mut cirques = Vec::new();
    for n in 0..NODES {
        if net.kind[n] != K_LAND || sc.q[n] <= 0.0 || has_donor[n] || present { continue; } // Lane GLACIER-1: no cirques
        let lip = if rcv[n] != NONE { fill[rcv[n] as usize] } else { fill[n] };
        let head = fill[n];
        // Gate: a cirque where the glacier head stands within 2 D_c above the cycle-mean ELA (heads far above it are
        // accumulation basins of larger ice bodies, not cirques; declared).
        if !(head > ela_char[n] && head - ela_char[n] <= 2.0 * dc) { continue; }
        let z = (head - dc).max(head.min(ela_char[n])).max(lip - psi * s_s[n] * rr);
        if head - z < 1.0 { continue; }
        let (ux, uy) = if rcv[n] != NONE { let r = rcv[n] as usize; let (dx, dy) = (net.pos[r][0] - net.pos[n][0], net.pos[r][1] - net.pos[n][1]); let l = dx.hypot(dy).max(1e-9); (dx / l, dy / l) } else { (1.0, 0.0) };
        cirques.push(Cirque { node: n as u32, x_m: net.pos[n][0], y_m: net.pos[n][1], ux, uy, floor_m: z, flat_m: z.max(lip), radius_m: rr, depth_m: dc, head_bed_m: head, ela_char_m: ela_char[n], lip_m: lip });
    }
    // Spacing: heads nearest the ELA first (then node index), drop any within the spacing of an accepted one.
    let n_cand = cirques.len();
    cirques.sort_by(|a, b| (a.head_bed_m - a.ela_char_m).total_cmp(&(b.head_bed_m - b.ela_char_m)).then(a.node.cmp(&b.node)));
    let sp = gp.cirque_spacing_mm as f64 / 1000.0;
    let cell = sp.max(1.0);
    let mut grid: std::collections::BTreeMap<(i64, i64), Vec<usize>> = std::collections::BTreeMap::new();
    let mut kept: Vec<Cirque> = Vec::new();
    for c in cirques.into_iter() {
        let (ix, iy) = ((c.x_m / cell).floor() as i64, (c.y_m / cell).floor() as i64);
        let mut ok = true;
        'o: for dy in -1..=1 { for dx in -1..=1 { if let Some(v) = grid.get(&(ix + dx, iy + dy)) { for &k in v { if (kept[k].x_m - c.x_m).hypot(kept[k].y_m - c.y_m) < sp { ok = false; break 'o; } } } } }
        if ok { grid.entry((ix, iy)).or_default().push(kept.len()); kept.push(c); }
    }
    kept.sort_by(|a, b| a.node.cmp(&b.node));
    let cirques = kept;
    // Hanging steps at glaciated network confluences (construction check, stats only).
    let mut steps: Vec<f64> = Vec::new();
    for n in 0..NODES {
        if net.kind[n] != K_LAND || w[n] <= 0.0 || rcv[n] == NONE { continue; }
        let r = rcv[n] as usize;
        if w[r] <= 0.0 { continue; }
        // A tributary: not the receiver's widest donor.
        let s = fill[n] - fill[r];
        if s > 1.0 { steps.push(s); }
    }
    steps.sort_by(|a, b| a.total_cmp(b));
    let pct = |v: &[f64], p: f64| -> f64 { if v.is_empty() { 0.0 } else { v[((v.len() - 1) as f64 * p).round() as usize] } };
    let n_q = |q: &[f64]| q.iter().filter(|v| **v > 0.0).count();
    t.stats = json!({
        "nodes_ice_char": n_q(&sc.q), "nodes_ice_max": n_q(&sm.q), "nodes_occupied": occ.iter().filter(|v| **v > 0.0).count(),
        "nodes_section": n_section, "nodes_carved": n_carved, "carve_max_m": carve_max, "nodes_carve_time_capped": n_capped, "nodes_unresolved_or_sheet": n_narrow, "nodes_raised_by_adverse_slope_limit": n_alley,
        "q_char_max_m3_yr": sc.q.iter().cloned().fold(0.0, f64::max), "q_max_max_m3_yr": sm.q.iter().cloned().fold(0.0, f64::max),
        "w_max_m": w.iter().cloned().fold(0.0, f64::max), "h_max_m": h.iter().cloned().fold(0.0, f64::max),
        "f_cycle_max": f_cyc.iter().cloned().fold(0.0, f64::max), "thickness_nye_limited_state_nodes": nye_limited,
        "cirques": cirques.len(), "cirque_candidates": n_cand, "fill_raised_nodes": (0..NODES).filter(|&n| fill[n] > bed[n]).count(),
        "link_steps_gt_1m": steps.len(), "link_step_p50_m": pct(&steps, 0.5), "link_step_p90_m": pct(&steps, 0.9), "link_step_max_m": pct(&steps, 1.0),
        "seconds_not_digested": 0.0,
    });
    if present { t.stats.as_object_mut().unwrap().insert("present_ice_only".into(), json!(true)); }
    let _ = t_start;
    t.q_char = f32v(&sc.q); t.q_max = f32v(&sm.q); t.f_cycle = f32v(&f_cyc); t.occupation = f32v(&occ);
    t.w = f32v(&w); t.h = f32v(&h); t.bed = f32v(&bed); t.bed_fill = f32v(&fill); t.carve = f32v(&carve);
    t.zs_char = f32v(&sc.zs.iter().map(|v| if v.is_finite() { *v } else { 0.0 }).collect::<Vec<_>>());
    t.zs_max = f32v(&sm.zs.iter().map(|v| if v.is_finite() { *v } else { 0.0 }).collect::<Vec<_>>());
    t.ela_char = f32v(&ela_char);
    t.s_surface = f32v(&s_s);
    t.cirques = cirques;
    t.now = now;
    Ok(t)
}

impl GlacialTable {
    /// Land node beds become the filled glacial bed (channels and floors then follow the glacial long profile).
    pub fn apply_beds(&self, kind: &[u8], bed: &mut [f64]) {
        for n in 0..NODES { if kind[n] == K_LAND && self.w[n] > 0.0 { bed[n] = self.bed_fill[n] as f64; } }
        // Land nodes downstream-consistent: a non-glaciated land node keeps its fluvial bed (never below its receiver by
        // construction of the fluvial rule, and receivers only went down).
    }
    fn arrays(&self) -> Vec<(&'static str, &Vec<f32>)> {
        vec![("gl_q_char_m3_yr", &self.q_char), ("gl_q_max_m3_yr", &self.q_max), ("gl_f_cycle", &self.f_cycle), ("gl_occupation", &self.occupation),
            ("gl_w_m", &self.w), ("gl_h_m", &self.h), ("gl_bed_m", &self.bed), ("gl_bed_fill_m", &self.bed_fill), ("gl_carve_m", &self.carve),
            ("gl_zs_char_m", &self.zs_char), ("gl_zs_max_m", &self.zs_max), ("gl_ela_char_m", &self.ela_char), ("gl_s_surface", &self.s_surface)]
    }
    /// Layer blobs (fields/<name>.f32 and fields/gl_rcv.u32).
    pub fn blobs(&self) -> Vec<(String, Vec<u8>)> {
        let mut v: Vec<(String, Vec<u8>)> = self.arrays().into_iter().map(|(n, a)| (format!("fields/{n}.f32"), a.iter().flat_map(|x| x.to_le_bytes()).collect())).collect();
        v.push(("fields/gl_rcv.u32".into(), self.rcv.iter().flat_map(|x| x.to_le_bytes()).collect()));
        // Lane GLACIER-1 (gl_now=1 only).
        if let Some(nw) = &self.now {
            let le = |a: &[f32]| a.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
            for (name, a) in NOW_FIELDS.iter().zip([&nw.q, &nw.w, &nw.h, &nw.zs, &nw.ela_grid, &nw.p_grid]) { v.push((format!("fields/{name}.f32"), le(a))); }
            for (name, a) in SHEET_FIELDS.iter().zip([&nw.sheet_zs, &nw.sheet_h]) { v.push((format!("fields/{name}.f32"), le(a))); }
            // Lane GLACIER-2 (gl_bal=1 only).
            if !nw.sheet_q.is_empty() { v.push((format!("fields/{}.f32", SHEET_BAL_FIELDS[0]), le(&nw.sheet_q))); v.push((format!("fields/{}.u8", SHEET_BAL_FIELDS[1]), nw.sheet_term.clone())); }
        }
        v
    }
    pub fn header(&self) -> Value {
        let mut h = json!({"cirques": self.cirques, "stats": self.stats, "fields": GL_FIELDS, "note": "lane GLACIAL-T1 glacier table (t1b_glacial.rs): per lod-8 node; windows read it, never recompute it"});
        // Lane GLACIER-1 (gl_now=1 only): the NOW state and the icefield sheet.
        if let Some(nw) = &self.now {
            h.as_object_mut().unwrap().insert("now".into(), json!({"stats": nw.stats, "node_fields": NOW_FIELDS, "sheet_fields": SHEET_FIELDS,
                "sheet": {"dx_mm": nw.dx_mm, "gx0": nw.gx0, "gy0": nw.gy0, "nx": nw.nx, "ny": nw.ny, "layout": "row-major, post (gx0 + ix, gy0 + iy) at (gx, gy) x dx_mm chart mm; zs = b off the mask (NaN off land); h >= 0 on the mask (0 on margin posts), -1 off the mask"},
                "note": "lane GLACIER-1 present ice (t1b_glacial.rs): node fields on the lod-8 lattice (ela / p grids on the regular node grid); windows read the sheet pointwise, never recompute it"}));
            // Lane GLACIER-2 (gl_bal=1 only): the balance-flux sheet's extra fields.
            if !nw.sheet_q.is_empty() { h["now"].as_object_mut().unwrap().insert("sheet_bal_fields".into(), json!({"fields": SHEET_BAL_FIELDS, "note": "lane GLACIER-2 balance_route_v1: ice_sheet_q_m3_yr f32 balance flux leaving the post (m3 w.e./yr, 0 off the ice); ice_sheet_term u8 terminus class of the post's glacier: 0 none, 1 land-terminating, 2 tidewater, 3 lake calving, 4 leaves the chart; same sub-grid and layout as the sheet"})); }
        }
        h
    }
    /// Read a published table (the layer loader has verified every field digest).
    pub fn load(dir: &Path) -> R<GlacialTable> {
        let b = std::fs::read(dir.join("rivernet.json")).map_err(|e| format!("rivernet.json: {e}"))?;
        let v: Value = serde_json::from_slice(&b).map_err(|e| e.to_string())?;
        let gh = v["header"].get("glacial").ok_or("layer has no glacial table (built without glacial=); refusing")?;
        let rd = |name: &str| -> R<Vec<f32>> { let (x, _) = read_f32(&dir.join(format!("fields/{name}.f32")), NODES)?; Ok(x) };
        let rb = std::fs::read(dir.join("fields/gl_rcv.u32")).map_err(|e| e.to_string())?;
        if rb.len() != 4 * NODES { return Err("gl_rcv size".into()); }
        let mut t = GlacialTable {
            rcv: rb.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect(),
            q_char: rd("gl_q_char_m3_yr")?, q_max: rd("gl_q_max_m3_yr")?, f_cycle: rd("gl_f_cycle")?, occupation: rd("gl_occupation")?,
            w: rd("gl_w_m")?, h: rd("gl_h_m")?, bed: rd("gl_bed_m")?, bed_fill: rd("gl_bed_fill_m")?, carve: rd("gl_carve_m")?,
            zs_char: rd("gl_zs_char_m")?, zs_max: rd("gl_zs_max_m")?, ela_char: rd("gl_ela_char_m")?, s_surface: rd("gl_s_surface")?, cirques: Vec::new(), stats: gh["stats"].clone(), now: None,
        };
        t.cirques = serde_json::from_value(gh["cirques"].clone()).map_err(|e| e.to_string())?;
        // Lane GLACIER-1 (only when the layer was built with gl_now=1).
        if let Some(nh) = gh.get("now") {
            let sh = &nh["sheet"];
            let gi = |k: &str| -> R<i64> { sh[k].as_i64().ok_or(format!("glacial now sheet {k}")) };
            let (nx, ny) = (gi("nx")? as usize, gi("ny")? as usize);
            let rs = |name: &str| -> R<Vec<f32>> { let b = std::fs::read(dir.join(format!("fields/{name}.f32"))).map_err(|e| format!("{name}: {e}"))?; if b.len() != 4 * nx * ny { return Err(format!("{name} size")); } Ok(b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()) };
            let rn = |name: &str| -> R<Vec<f32>> { let b = std::fs::read(dir.join(format!("fields/{name}.f32"))).map_err(|e| format!("{name}: {e}"))?; if b.len() != 4 * NODES { return Err(format!("{name} size")); } Ok(b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()) };
            t.now = Some(IceNow { q: rn(NOW_FIELDS[0])?, w: rn(NOW_FIELDS[1])?, h: rn(NOW_FIELDS[2])?, zs: rn(NOW_FIELDS[3])?, ela_grid: rn(NOW_FIELDS[4])?, p_grid: rn(NOW_FIELDS[5])?,
                dx_mm: gi("dx_mm")?, gx0: gi("gx0")?, gy0: gi("gy0")?, nx, ny, sheet_zs: rs(SHEET_FIELDS[0])?, sheet_h: rs(SHEET_FIELDS[1])?, stats: nh["stats"].clone(), sheet_q: Vec::new(), sheet_term: Vec::new() });
            // Lane GLACIER-2 (layers built with gl_bal=1).
            if nh.get("sheet_bal_fields").is_some() {
                let tb = std::fs::read(dir.join(format!("fields/{}.u8", SHEET_BAL_FIELDS[1]))).map_err(|e| format!("{}: {e}", SHEET_BAL_FIELDS[1]))?;
                if tb.len() != nx * ny { return Err("ice_sheet_term size".into()); }
                let nw = t.now.as_mut().unwrap();
                nw.sheet_q = rs(SHEET_BAL_FIELDS[0])?; nw.sheet_term = tb;
            }
        }
        Ok(t)
    }
    /// Count of nodes whose section or filled bed differ (f32 bits) from another table (layer re-build check).
    pub fn diff_nodes(&self, o: &GlacialTable) -> usize {
        (0..NODES).filter(|&n| self.w[n].to_bits() != o.w[n].to_bits() || self.h[n].to_bits() != o.h[n].to_bits() || self.bed_fill[n].to_bits() != o.bed_fill[n].to_bits()).count()
    }
}

#[derive(PartialEq)]
struct HeapK(f64, usize);
impl Eq for HeapK {}
impl PartialOrd for HeapK { fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> { Some(self.cmp(o)) } }
impl Ord for HeapK { fn cmp(&self, o: &Self) -> std::cmp::Ordering { o.0.total_cmp(&self.0).then(o.1.cmp(&self.1)) } }

/// Lane GLACIER-1 `plastic_sheet_v1`: the surface z_s of perfectly plastic ice over bed `b` on the `mask` posts of `lat`
/// (NaN elsewhere): z_s = b on mask posts with an 8-neighbour outside the mask or off the lattice, Dijkstra inward with the
/// exact local update (u - z_j)(u - b_i) = h0 d. Returns (z_s, margin post count).
pub fn plastic_march(lat: &EvalLat, b: &[f64], mask: &[bool], h0: f64) -> (Vec<f64>, u64) {
    let (nn, dx) = (lat.nodes(), lat.dx());
    let mut zs = vec![f64::NAN; nn];
    let mut done = vec![false; nn];
    let mut heap = std::collections::BinaryHeap::new();
    let nb = |i: usize, k: usize| -> Option<usize> { let (gx, gy) = lat.g(i); let (ox, oy) = crate::D8[k]; lat.idx(gx + ox, gy + oy) };
    let mut n_margin = 0u64;
    for i in 0..nn {
        if !mask[i] { continue; }
        if (0..8).any(|k| nb(i, k).map(|j| !mask[j]).unwrap_or(true)) { zs[i] = b[i]; heap.push(HeapK(b[i], i)); n_margin += 1; }
    }
    while let Some(HeapK(z, i)) = heap.pop() {
        if done[i] { continue; }
        done[i] = true;
        for k in 0..8 {
            let Some(j) = nb(i, k) else { continue };
            if !mask[j] || done[j] { continue; }
            let d = if k % 2 == 1 { dx * std::f64::consts::SQRT_2 } else { dx };
            let u = 0.5 * ((z + b[j]) + ((z - b[j]).powi(2) + 4.0 * h0 * d).sqrt());
            if !(u >= zs[j]) { zs[j] = u; heap.push(HeapK(u, j)); }
        }
    }
    (zs, n_margin)
}

/// Lane GLACIER-1 (gl_now=1): the chart icefield sheet (`ice_mask_v1`, `plastic_sheet_v1`, `mass_check_v1`) on the 128 m
/// solve lattice, from the pre-glacial evaluation surface `b` and its codes `fixed`; fills `t.now` (grids, sheet, stats).
#[allow(clippy::too_many_arguments)]
pub fn build_sheet(t: &mut GlacialTable, net: &Network, lat: &EvalLat, b: &[f64], fixed: &[u8], chart: &vz_tilekey::chart::Chart, t0: &T0, g: &T0Glacial, gp: &GlacialParamsV1) -> R<()> {
    let rcv = t.rcv.clone();
    let floor_n = t.bed.clone();
    let Some(nw) = t.now.as_mut() else { return Ok(()) };
    // ELA_now and P on the regular node grid.
    let mut ela_grid = vec![0f32; NODES];
    let mut p_grid = vec![0f32; NODES];
    for n in 0..NODES {
        let (mx, my) = crate::node_mm(n);
        let (la, lo) = chart.inverse_mm(mx, my).map_err(|e| e.to_string())?;
        let (la, lo) = (la * vz_tilekey::chart::RAD2DEG, lo * vz_tilekey::chart::RAD2DEG);
        ela_grid[n] = crate::t1b::bilinear(t0, &|c: u32| g.ela_now[c as usize] as f64, la, lo, false) as f32;
        p_grid[n] = crate::t1b::bilinear(t0, &|c: u32| g.precip[c as usize] as f64, la, lo, false) as f32;
    }
    let (n, dx) = (lat.n, lat.dx());
    let nn = lat.nodes();
    let pos = |i: usize| { let (gx, gy) = lat.g(i); (gx as f64 * dx, gy as f64 * dx) };
    let land: Vec<bool> = (0..nn).map(|i| b[i].is_finite() && fixed[i] != crate::t1b::MARINE && fixed[i] != crate::t1b::LAKE).collect();
    // Sheet bed despike (declared): a land post more than 500 m above the median of its finite 8 neighbours takes that median
    // (1-post spikes of the 128 m evaluation surface, which the 32 m windows do not carry; counted in the stats).
    let mut n_spike = 0u64;
    let bd: Vec<f64> = (0..nn).map(|i| {
        if !land[i] { return b[i]; }
        let (gx, gy) = lat.g(i);
        let mut v: Vec<f64> = crate::D8.iter().filter_map(|(ox, oy)| lat.idx(gx + ox, gy + oy)).map(|j| b[j]).filter(|z| z.is_finite()).collect();
        if v.len() < 3 { return b[i]; }
        v.sort_by(|a, c| a.total_cmp(c));
        let med = if v.len() % 2 == 1 { v[v.len() / 2] } else { 0.5 * (v[v.len() / 2 - 1] + v[v.len() / 2]) };
        if b[i] - med > 500.0 { n_spike += 1; med } else { b[i] }
    }).collect();
    let b = &bd[..];
    // Link floor for the tongue confinement: this surface at the post nearest to each polyline point (the drawn channel; on a
    // volcano the network bed lies under the edifice, the surface does not); the table bed where that post is not land.
    let owner = |i: usize| -> Option<usize> { let (gx, gy) = lat.g(i); crate::node_index((gx + 2).div_euclid(4), (gy + 2).div_euclid(4)) };
    let surf_at = |q: [f64; 2], fallback: f64| -> f64 { match lat.idx((q[0] / dx).round() as i64, (q[1] / dx).round() as i64) { Some(i) if land[i] => b[i], _ => fallback } };
    let ela: Vec<f64> = (0..nn).map(|i| { let (x, y) = pos(i); node_grid_bilinear(&ela_grid, x, y) }).collect();
    let acc_bed: Vec<bool> = (0..nn).map(|i| land[i] && b[i] >= ela[i]).collect();
    let h0 = gp.yield_stress_kpa as f64 * 1000.0 / (gp.ice_density_kg_m3 as f64 * 9.81);
    // Outlet tongues: posts within W/2 of a link whose node carries ice now, confined to the valley: bed <= link floor +
    // sqrt(h0 W) (the centre thickness of a plastic flat-bed cross-profile of width W; declared).
    let mut strip = vec![false; nn];
    let (lx0, ly0) = (lat.gx0 as f64 * dx, lat.gy0 as f64 * dx);
    let (lx1, ly1) = ((lat.gx0 + n as i64 - 1) as f64 * dx, (lat.gy0 + n as i64 - 1) as f64 * dx);
    let mut links = 0u64;
    for k in 0..NODES {
        if nw.q[k] <= 0.0 || nw.w[k] <= 0.0 { continue; }
        if rcv[k] == NONE { continue; } // ice leaves the chart or the network here: no tongue drawn beyond the node
        let w0 = nw.w[k] as f64;
        let r = rcv[k] as usize;
        let w1 = if nw.q[r] > 0.0 { nw.w[r] as f64 } else { 0.0 };
        let (f0, f1) = (floor_n[k] as f64, floor_n[r] as f64);
        let (poly, fl): (Vec<[f64; 2]>, Vec<f64>) = if net.kind[k] == K_LAND { (vec![net.pos[k], net.mid[k], net.pos[r]], vec![surf_at(net.pos[k], f0), surf_at(net.mid[k], 0.5 * (f0 + f1)), surf_at(net.pos[r], f1)]) } else { (vec![net.pos[k], net.pos[r]], vec![surf_at(net.pos[k], f0), surf_at(net.pos[r], f1)]) };
        let rad = 0.5 * w0.max(w1);
        let (mut bx0, mut bx1, mut by0, mut by1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for q in &poly { bx0 = bx0.min(q[0]); bx1 = bx1.max(q[0]); by0 = by0.min(q[1]); by1 = by1.max(q[1]); }
        let (bx0, bx1, by0, by1) = (bx0 - rad, bx1 + rad, by0 - rad, by1 + rad);
        if bx1 < lx0 || bx0 > lx1 || by1 < ly0 || by0 > ly1 { continue; }
        links += 1;
        let mut cum = vec![0.0f64; poly.len()];
        for j in 1..poly.len() { cum[j] = cum[j - 1] + (poly[j][0] - poly[j - 1][0]).hypot(poly[j][1] - poly[j - 1][1]); }
        let ltot = cum[poly.len() - 1].max(1e-9);
        let gx_lo = ((bx0 / dx).floor() as i64).max(lat.gx0);
        let gx_hi = ((bx1 / dx).ceil() as i64).min(lat.gx0 + n as i64 - 1);
        let gy_lo = ((by0 / dx).floor() as i64).max(lat.gy0);
        let gy_hi = ((by1 / dx).ceil() as i64).min(lat.gy0 + n as i64 - 1);
        for gy in gy_lo..=gy_hi { for gx in gx_lo..=gx_hi {
            let i = ((gy - lat.gy0) as usize) * n + (gx - lat.gx0) as usize;
            if !land[i] || strip[i] { continue; }
            let p = [gx as f64 * dx, gy as f64 * dx];
            let (mut y, mut s, mut fz) = (f64::INFINITY, 0.0, 0.0);
            for j in 1..poly.len() {
                let (dd, tt) = seg_dist(p, poly[j - 1], poly[j]);
                if dd < y { y = dd; s = (cum[j - 1] + tt * (cum[j] - cum[j - 1])) / ltot; fz = fl[j - 1] + tt * (fl[j] - fl[j - 1]); }
            }
            let w = w0 + (w1 - w0) * s;
            if y <= 0.5 * w && b[i] <= fz + (h0 * w).sqrt() { strip[i] = true; }
        } }
    }
    let mask: Vec<bool> = (0..nn).map(|i| acc_bed[i] || strip[i]).collect();
    // Plastic surface: Dijkstra from the margin.
    let (zs, n_margin) = plastic_march(lat, b, &mask, h0);
    // Mass check and statistics.
    let a_post = dx * dx;
    let (mut acc_all, mut acc_routed, mut acc_unrouted, mut abl_mask) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut hs: Vec<f64> = Vec::new();
    let (mut n_mask, mut n_acc, mut n_strip_only, mut n_surf_acc) = (0u64, 0u64, 0u64, 0u64);
    for i in 0..nn {
        if !land[i] { continue; }
        let (x, y) = pos(i);
        let pm = node_grid_bilinear(&p_grid, x, y).max(0.0) / 1000.0;
        let prod = pm * ice_share(b[i], ela[i], g.full_m) * a_post;
        if prod > 0.0 {
            acc_all += prod;
            if owner(i).map(|o| nw.q[o] > 0.0).unwrap_or(false) { acc_routed += prod; } else if mask[i] { acc_unrouted += prod; }
        } else if mask[i] { abl_mask += prod; }
        if !mask[i] { continue; }
        n_mask += 1;
        if acc_bed[i] { n_acc += 1; } else { n_strip_only += 1; }
        if zs[i] >= ela[i] { n_surf_acc += 1; }
        hs.push((zs[i] - b[i]).max(0.0));
    }
    hs.sort_by(|a, b| a.total_cmp(b));
    let pct = |p: f64| -> f64 { if hs.is_empty() { 0.0 } else { hs[((hs.len() - 1) as f64 * p).round() as usize] } };
    // Crop to the mask's bounding box (+2 posts).
    let (mut ix0, mut ix1, mut iy0, mut iy1) = (usize::MAX, 0usize, usize::MAX, 0usize);
    for i in 0..nn { if mask[i] { let (x, y) = (i % n, i / n); ix0 = ix0.min(x); ix1 = ix1.max(x); iy0 = iy0.min(y); iy1 = iy1.max(y); } }
    if n_mask == 0 { ix0 = 0; ix1 = 0; iy0 = 0; iy1 = 0; }
    let (ix0, iy0) = (ix0.saturating_sub(2), iy0.saturating_sub(2));
    let (ix1, iy1) = ((ix1 + 2).min(n - 1), (iy1 + 2).min(n - 1));
    let (nx, ny) = (ix1 - ix0 + 1, iy1 - iy0 + 1);
    let (mut szs, mut sh) = (Vec::with_capacity(nx * ny), Vec::with_capacity(nx * ny));
    for y in iy0..=iy1 { for x in ix0..=ix1 {
        let i = y * n + x;
        if mask[i] { szs.push(zs[i] as f32); sh.push((zs[i] - b[i]).max(0.0) as f32); } else { szs.push(if land[i] { b[i] as f32 } else { f32::NAN }); sh.push(-1.0); }
    } }
    nw.ela_grid = ela_grid; nw.p_grid = p_grid;
    nw.dx_mm = lat.dx_mm; nw.gx0 = lat.gx0 + ix0 as i64; nw.gy0 = lat.gy0 + iy0 as i64; nw.nx = nx; nw.ny = ny;
    nw.sheet_zs = szs; nw.sheet_h = sh;
    let km2 = a_post / 1e6;
    if let Some(o) = nw.stats.as_object_mut() {
        o.insert("sheet".into(), json!({"h0_m": h0, "despiked_bed_posts": n_spike, "links_rasterized": links, "mask_posts": n_mask, "mask_area_km2": n_mask as f64 * km2, "bed_above_ela_area_km2": n_acc as f64 * km2,
            "tongue_only_area_km2": n_strip_only as f64 * km2, "margin_posts": n_margin, "surface_above_ela_area_km2": n_surf_acc as f64 * km2,
            "aar_surface": if n_mask > 0 { n_surf_acc as f64 / n_mask as f64 } else { 0.0 }, "aar_bed": if n_mask > 0 { n_acc as f64 / n_mask as f64 } else { 0.0 },
            "h_p10_m": pct(0.1), "h_p50_m": pct(0.5), "h_p90_m": pct(0.9), "h_p99_m": pct(0.99), "h_max_m": pct(1.0),
            "mass": {"accumulation_all_m3_yr": acc_all, "accumulation_routed_m3_yr": acc_routed, "accumulation_unrouted_in_mask_m3_yr": acc_unrouted,
                "unrouted_over_routed": if acc_routed > 0.0 { acc_unrouted / acc_routed } else { 0.0 }, "ablation_potential_in_mask_m3_yr": abl_mask,
                "gate_unrouted_below_10_percent": acc_routed > 0.0 && acc_unrouted / acc_routed < 0.1},
            "seconds_not_digested": 0.0}));
    }
    Ok(())
}

/// Lane GLACIER-2 `plastic_sheet_seeded_v1`: plastic_march with a per-margin-post seed surface (>= bed).
pub fn plastic_march_seeded(lat: &EvalLat, b: &[f64], mask: &[bool], h0: f64, seed: &dyn Fn(usize) -> f64) -> (Vec<f64>, u64) {
    let (nn, dx) = (lat.nodes(), lat.dx());
    let mut zs = vec![f64::NAN; nn];
    let mut done = vec![false; nn];
    let mut heap = std::collections::BinaryHeap::new();
    let nb = |i: usize, k: usize| -> Option<usize> { let (gx, gy) = lat.g(i); let (ox, oy) = crate::D8[k]; lat.idx(gx + ox, gy + oy) };
    let mut n_margin = 0u64;
    for i in 0..nn {
        if !mask[i] { continue; }
        if (0..8).any(|k| nb(i, k).map(|j| !mask[j]).unwrap_or(true)) { let z = seed(i).max(b[i]); zs[i] = z; heap.push(HeapK(z, i)); n_margin += 1; }
    }
    while let Some(HeapK(z, i)) = heap.pop() {
        if done[i] { continue; }
        done[i] = true;
        for k in 0..8 {
            let Some(j) = nb(i, k) else { continue };
            if !mask[j] || done[j] { continue; }
            let d = if k % 2 == 1 { dx * std::f64::consts::SQRT_2 } else { dx };
            let u = 0.5 * ((z + b[j]) + ((z - b[j]).powi(2) + 4.0 * h0 * d).sqrt());
            if !(u >= zs[j]) { zs[j] = u; heap.push(HeapK(u, j)); }
        }
    }
    (zs, n_margin)
}

/// Lane GLACIER-3 `plastic_sheet_fmm_v1`: the upwind plastic update at a post with bed `b`, accepted-neighbour minima `ax`
/// (x axis) and `ay` (y axis; +inf where none), spacing `dx`: u with (u - b)^2 (max(0, u - ax)^2 + max(0, u - ay)^2) =
/// (h0 dx)^2, by bisection between max(b, min(ax, ay)) (residual < 0) and the one-sided root (residual >= 0).
pub fn plastic_update_fmm(b: f64, ax: f64, ay: f64, h0: f64, dx: f64) -> f64 {
    let m = ax.min(ay);
    let hi = 0.5 * ((m + b) + ((m - b).powi(2) + 4.0 * h0 * dx).sqrt());
    if !(ax.is_finite() && ay.is_finite()) { return hi; }
    let t = (h0 * dx).powi(2);
    let g = |u: f64| { let (p, q) = ((u - ax).max(0.0), (u - ay).max(0.0)); (u - b).powi(2) * (p * p + q * q) - t };
    let (mut lo, mut hi) = (b.max(m), hi);
    if g(lo) >= 0.0 { return lo; }
    for _ in 0..52 { let mid = 0.5 * (lo + hi); if g(mid) >= 0.0 { hi = mid; } else { lo = mid; } }
    hi
}

/// Lane GLACIER-3 `plastic_sheet_fmm_v1`: plastic_march_seeded's law and margin seeds, solved by the fast marching method on
/// the 4-neighbour upwind stencil (Euclidean fronts). Returns (z_s, margin post count).
pub fn plastic_march_fmm(lat: &EvalLat, b: &[f64], mask: &[bool], h0: f64, seed: &dyn Fn(usize) -> f64) -> (Vec<f64>, u64) {
    let (nn, dx) = (lat.nodes(), lat.dx());
    let mut zs = vec![f64::NAN; nn];
    let mut done = vec![false; nn];
    let mut heap = std::collections::BinaryHeap::new();
    let nb = |i: usize, k: usize| -> Option<usize> { let (gx, gy) = lat.g(i); let (ox, oy) = crate::D8[k]; lat.idx(gx + ox, gy + oy) };
    let mut n_margin = 0u64;
    for i in 0..nn {
        if !mask[i] { continue; }
        if (0..8).any(|k| nb(i, k).map(|j| !mask[j]).unwrap_or(true)) { let z = seed(i).max(b[i]); zs[i] = z; heap.push(HeapK(z, i)); n_margin += 1; }
    }
    // Accepted value of a neighbour (inf if none / not accepted).
    let acc = |zs: &[f64], done: &[bool], j: Option<usize>| -> f64 { match j { Some(j) if done[j] => zs[j], _ => f64::INFINITY } };
    while let Some(HeapK(z, i)) = heap.pop() {
        if done[i] || z != zs[i] { continue; }
        done[i] = true;
        for k in [0usize, 2, 4, 6] {
            let Some(j) = nb(i, k) else { continue };
            if !mask[j] || done[j] { continue; }
            let ax = acc(&zs, &done, nb(j, 0)).min(acc(&zs, &done, nb(j, 4)));
            let ay = acc(&zs, &done, nb(j, 2)).min(acc(&zs, &done, nb(j, 6)));
            let u = plastic_update_fmm(b[j], ax, ay, h0, dx);
            if !(u >= zs[j]) { zs[j] = u; heap.push(HeapK(u, j)); }
        }
    }
    (zs, n_margin)
}

/// Lane GLACIER-3 `balance_route_dinf_v1` routing: (pop order, primary receiver, secondary receiver (NONE if none), share of
/// the primary, flow direction unit vector per post). Primary = the larger share (the cardinal one on a tie).
pub struct DinfRoute { pub pop: Vec<u32>, pub r1: Vec<u32>, pub r2: Vec<u32>, pub f1: Vec<f32>, pub ux: Vec<f32>, pub uy: Vec<f32> }

pub fn bal_route_dinf(lat: &EvalLat, s: &[f64], dom: &[bool], lake: &[bool]) -> DinfRoute {
    let (nn, dx) = (lat.nodes(), lat.dx());
    let nb = |i: usize, k: usize| -> Option<usize> { let (gx, gy) = lat.g(i); let (ox, oy) = crate::D8[k]; lat.idx(gx + ox, gy + oy) };
    let mut f = vec![f64::NAN; nn];
    let mut parent = vec![NONE; nn];
    let mut seen = vec![false; nn];
    let mut heap = std::collections::BinaryHeap::new();
    for i in 0..nn {
        if !dom[i] { continue; }
        if (0..8).any(|k| nb(i, k).map(|j| !dom[j]).unwrap_or(true)) { f[i] = s[i]; seen[i] = true; heap.push(HeapK(s[i], i)); }
    }
    let mut pop: Vec<u32> = Vec::with_capacity(nn);
    while let Some(HeapK(z, i)) = heap.pop() {
        pop.push(i as u32);
        for k in 0..8 {
            let Some(j) = nb(i, k) else { continue };
            if !dom[j] || seen[j] { continue; }
            seen[j] = true;
            f[j] = s[j].max(z);
            parent[j] = i as u32;
            heap.push(HeapK(f[j], j));
        }
    }
    let (mut r1, mut r2, mut f1, mut ux, mut uy) = (vec![NONE; nn], vec![NONE; nn], vec![1f32; nn], vec![0f32; nn], vec![0f32; nn]);
    let q4 = std::f64::consts::FRAC_PI_4;
    let sq2 = std::f64::consts::SQRT_2;
    for &iu in &pop {
        let i = iu as usize;
        let (gx, gy) = lat.g(i);
        let e0 = f[i];
        // Facets: cardinal c (D8 index 0, 2, 4, 6) and the diagonal on either side (c + 1, c - 1).
        let (mut best_s, mut best) = (0.0f64, None::<(Option<usize>, Option<usize>, f64, f64, f64)>); // (cardinal, diagonal, share cardinal, ux, uy)
        for c in [0usize, 2, 4, 6] {
            for dg in [(c + 1) % 8, (c + 7) % 8] {
                let (cx, cy) = crate::D8[c];
                let (dxo, dyo) = crate::D8[dg];
                let (px, py) = ((dxo - cx) as f64, (dyo - cy) as f64); // unit perpendicular from the cardinal to the diagonal
                let j1 = lat.idx(gx + cx, gy + cy).filter(|&j| dom[j]);
                let j2 = lat.idx(gx + dxo, gy + dyo).filter(|&j| dom[j]);
                let cand = match (j1, j2) {
                    (Some(a), Some(b)) => {
                        let s1 = (e0 - f[a]) / dx;
                        let s2 = (f[a] - f[b]) / dx;
                        let mut r = s2.atan2(s1);
                        let mut sl = s1.hypot(s2);
                        if r < 0.0 { r = 0.0; sl = s1; } else if r > q4 { r = q4; sl = (e0 - f[b]) / (dx * sq2); }
                        let sh = 1.0 - r / q4;
                        (sl, Some(a), Some(b), sh, cx as f64 * r.cos() + px * r.sin(), cy as f64 * r.cos() + py * r.sin())
                    }
                    (Some(a), None) => ((e0 - f[a]) / dx, Some(a), None, 1.0, cx as f64, cy as f64),
                    (None, Some(b)) => ((e0 - f[b]) / (dx * sq2), None, Some(b), 0.0, dxo as f64 / sq2, dyo as f64 / sq2),
                    (None, None) => continue,
                };
                if cand.0 > best_s { best_s = cand.0; best = Some((cand.1, cand.2, cand.3, cand.4, cand.5)); }
            }
        }
        if let Some((a, b, sh, vx, vy)) = best {
            let a = if sh > 0.0 { a } else { None };
            let b = if sh < 1.0 { b } else { None };
            match (a, b) {
                (Some(a), Some(b)) => { if sh >= 0.5 { r1[i] = a as u32; r2[i] = b as u32; f1[i] = sh as f32; } else { r1[i] = b as u32; r2[i] = a as u32; f1[i] = (1.0 - sh) as f32; } }
                (Some(a), None) => { r1[i] = a as u32; }
                (None, Some(b)) => { r1[i] = b as u32; }
                (None, None) => {}
            }
            let nv = vx.hypot(vy).max(1e-12);
            ux[i] = (vx / nv) as f32; uy[i] = (vy / nv) as f32;
        }
        if r1[i] == NONE {
            // Flats and filled pits: the pop parent; else (outlet posts) the adjacent lake post of lowest stage; else NONE.
            let mut lk = NONE;
            for k in 0..8 { if let Some(j) = nb(i, k) { if !dom[j] && lake[j] && (lk == NONE || s[j] < s[lk as usize] || (s[j] == s[lk as usize] && (j as u32) < lk)) { lk = j as u32; } } }
            let t = if parent[i] != NONE { parent[i] } else { lk };
            r1[i] = t;
            if t != NONE { let (tx, ty) = lat.g(t as usize); let (vx, vy) = ((tx - gx) as f64, (ty - gy) as f64); let nv = vx.hypot(vy).max(1e-12); ux[i] = (vx / nv) as f32; uy[i] = (vy / nv) as f32; }
        }
    }
    DinfRoute { pop, r1, r2, f1, ux, uy }
}

/// Lane GLACIER-2 balance routing over the domain posts `dom` (land and marine): (pop order of the priority flood, receiver
/// per post). Receiver: the steepest-descent domain 8-neighbour of the filled surface, else the pop parent (flats, filled
/// pits), else (outlet posts) the adjacent `lake` post of lowest stage, else NONE (chart edge or a void post).
pub fn bal_route(lat: &EvalLat, s: &[f64], dom: &[bool], lake: &[bool]) -> (Vec<u32>, Vec<u32>) {
    let (nn, dx) = (lat.nodes(), lat.dx());
    let nb = |i: usize, k: usize| -> Option<usize> { let (gx, gy) = lat.g(i); let (ox, oy) = crate::D8[k]; lat.idx(gx + ox, gy + oy) };
    let mut f = vec![f64::NAN; nn];
    let mut parent = vec![NONE; nn];
    let mut seen = vec![false; nn];
    let mut heap = std::collections::BinaryHeap::new();
    for i in 0..nn {
        if !dom[i] { continue; }
        if (0..8).any(|k| nb(i, k).map(|j| !dom[j]).unwrap_or(true)) { f[i] = s[i]; seen[i] = true; heap.push(HeapK(s[i], i)); }
    }
    let mut pop: Vec<u32> = Vec::with_capacity(nn);
    while let Some(HeapK(z, i)) = heap.pop() {
        pop.push(i as u32);
        for k in 0..8 {
            let Some(j) = nb(i, k) else { continue };
            if !dom[j] || seen[j] { continue; }
            seen[j] = true;
            f[j] = s[j].max(z);
            parent[j] = i as u32;
            heap.push(HeapK(f[j], j));
        }
    }
    let mut rcv = vec![NONE; nn];
    for &i in &pop {
        let i = i as usize;
        let (mut best, mut bs, mut lk) = (NONE, 0.0f64, NONE);
        for k in 0..8 {
            let d = if k % 2 == 1 { dx * std::f64::consts::SQRT_2 } else { dx };
            let Some(j) = nb(i, k) else { continue };
            if dom[j] { let sl = (f[i] - f[j]) / d; if sl > bs || (sl == bs && sl > 0.0 && (j as u32) < best) { bs = sl; best = j as u32; } }
            else if lake[j] && (lk == NONE || s[j] < s[lk as usize] || (s[j] == s[lk as usize] && (j as u32) < lk)) { lk = j as u32; }
        }
        rcv[i] = if best != NONE { best } else if parent[i] != NONE { parent[i] } else { lk };
    }
    (pop, rcv)
}

/// Lane GLACIER-2 (gl_bal=1): the balance-flux icefield (`balance_route_v1`, `ablation_claim_v1`, `flotation_calving_v1`,
/// `plastic_sheet_seeded_v1`) on the 128 m solve lattice from the pre-glacial evaluation surface `b`, its codes `fixed`
/// and the construction sea floor `seabed` (Surface.z_eik); fills `t.now` (grids, sheet, flux, terminus class, stats).
#[allow(clippy::too_many_arguments)]
pub fn build_sheet_bal(t: &mut GlacialTable, lat: &EvalLat, b: &[f64], fixed: &[u8], seabed: &[f64], t0_bil: &[f32], chart: &vz_tilekey::chart::Chart, t0: &T0, g: &T0Glacial, gp: &GlacialParamsV1) -> R<()> {
    let bp = gp.now.as_ref().and_then(|n| n.bal.clone()).ok_or("build_sheet_bal without gl_bal params")?;
    let Some(nw) = t.now.as_mut() else { return Ok(()) };
    let mut ela_grid = vec![0f32; NODES];
    let mut p_grid = vec![0f32; NODES];
    for n in 0..NODES {
        let (mx, my) = crate::node_mm(n);
        let (la, lo) = chart.inverse_mm(mx, my).map_err(|e| e.to_string())?;
        let (la, lo) = (la * vz_tilekey::chart::RAD2DEG, lo * vz_tilekey::chart::RAD2DEG);
        ela_grid[n] = crate::t1b::bilinear(t0, &|c: u32| g.ela_now[c as usize] as f64, la, lo, false) as f32;
        p_grid[n] = crate::t1b::bilinear(t0, &|c: u32| g.precip[c as usize] as f64, la, lo, false) as f32;
    }
    let (n, dx) = (lat.n, lat.dx());
    let nn = lat.nodes();
    let pos = |i: usize| { let (gx, gy) = lat.g(i); (gx as f64 * dx, gy as f64 * dx) };
    let nb = |i: usize, k: usize| -> Option<usize> { let (gx, gy) = lat.g(i); let (ox, oy) = crate::D8[k]; lat.idx(gx + ox, gy + oy) };
    // Post kinds: 0 land, 1 marine, 2 lake, 3 void (no finite surface).
    let kind: Vec<u8> = (0..nn).map(|i| if fixed[i] == crate::t1b::MARINE { 1 } else if fixed[i] == crate::t1b::LAKE { 2 } else if b[i].is_finite() { 0 } else { 3 }).collect();
    let land = |i: usize| kind[i] == 0;
    // Despike (as plastic_sheet_v1).
    let mut n_spike = 0u64;
    let bd: Vec<f64> = (0..nn).map(|i| {
        if !land(i) { return b[i]; }
        let (gx, gy) = lat.g(i);
        let mut v: Vec<f64> = crate::D8.iter().filter_map(|(ox, oy)| lat.idx(gx + ox, gy + oy)).map(|j| b[j]).filter(|z| z.is_finite()).collect();
        if v.len() < 3 { return b[i]; }
        v.sort_by(|a, c| a.total_cmp(c));
        let med = if v.len() % 2 == 1 { v[v.len() / 2] } else { 0.5 * (v[v.len() / 2 - 1] + v[v.len() / 2]) };
        if b[i] - med > 500.0 { n_spike += 1; med } else { b[i] }
    }).collect();
    let sea = t0.sea_level_m;
    // Water depth: the deeper of the construction sea floor (z_eik) and the T0 bathymetry bilinear (the construction water
    // level sea - F is a flooding device, not bathymetry: offshore, where it is shallower than T0, T0 governs).
    let depth: Vec<f64> = (0..nn).map(|i| if kind[i] == 1 {
        let dc = if seabed[i].is_finite() { sea - seabed[i] } else { 0.0 };
        let dt = if (t0_bil[i] as f64).is_finite() { sea - t0_bil[i] as f64 } else { 0.0 };
        dc.max(dt).max(0.0)
    } else { 0.0 }).collect();
    // Sheet bed: b on land, the construction sea floor on marine posts.
    let bed: Vec<f64> = (0..nn).map(|i| if kind[i] == 1 { sea - depth[i] } else { bd[i] }).collect();
    let ela: Vec<f64> = (0..nn).map(|i| { let (x, y) = pos(i); node_grid_bilinear(&ela_grid, x, y) }).collect();
    let pm: Vec<f64> = (0..nn).map(|i| { let (x, y) = pos(i); node_grid_bilinear(&p_grid, x, y).max(0.0) / 1000.0 }).collect();
    let h0 = gp.yield_stress_kpa as f64 * 1000.0 / (gp.ice_density_kg_m3 as f64 * 9.81);
    let (wr, qr, al) = (gp.width_ref_mm as f64 / 1000.0, gp.width_ref_flux_m3_per_yr as f64, gp.width_exponent_permille as f64 / 1000.0);
    let rho_i = gp.ice_density_kg_m3 as f64;
    let rho_w = bp.seawater_density_kg_m3 as f64;
    let qf = bp.buoyancy_excess_permille as f64 / 1000.0;
    let c_sea = bp.calving_coeff_sea_milli_per_yr as f64 / 1000.0;
    let to_ice = 1000.0 / rho_i;
    let width = |q_we: f64| wr * ((q_we * to_ice).max(0.0) / qr).powf(al);
    let h_face = |d: f64| (1.0 + qf) * rho_w / rho_i * d;
    let a_post = dx * dx;
    let dom: Vec<bool> = (0..nn).map(|i| kind[i] == 0 || kind[i] == 1).collect();
    let lakev: Vec<bool> = (0..nn).map(|i| kind[i] == 2).collect();
    let mut surf = bed.clone();
    let mut ice = vec![false; nn];
    let mut zs = vec![f64::NAN; nn];
    let mut qout = vec![0.0f64; nn];
    let mut term_cls = vec![TERM_NONE; nn];
    let mut glac = vec![NONE; nn];
    let mut pass_log: Vec<Value> = Vec::new();
    let mut prev_mask: Vec<bool> = Vec::new();
    let mut termini: Vec<Value> = Vec::new();
    let mut mass = json!({});
    // Lane GLACIER-3 (gl_geo=1): D-infinity routing, aligned marine claims, the FMM sheet, passes to a mask tolerance.
    let geo = bp.geo.clone();
    let passes = match &geo { Some(gq) => gq.max_passes.max(2) as usize, None => bp.passes.max(1) as usize };
    let mut converged_at: Option<usize> = None;
    let mut r2v: Vec<u32> = Vec::new();
    let mut f1v: Vec<f32> = Vec::new();
    let mut uxv: Vec<f32> = Vec::new();
    let mut uyv: Vec<f32> = Vec::new();
    for pass in 0..passes {
        let (pop, rcv) = match &geo {
            None => bal_route(lat, &surf, &dom, &lakev),
            Some(_) => { let d = bal_route_dinf(lat, &surf, &dom, &lakev); r2v = d.r2; f1v = d.f1; uxv = d.ux; uyv = d.uy; (d.pop, d.r1) }
        };
        let bal: Vec<f64> = (0..nn).map(|i| if land(i) { pm[i] * ice_share(surf[i], ela[i], g.full_m) } else if kind[i] == 1 { pm[i] * ice_share(sea, ela[i], g.full_m) } else { 0.0 }).collect();
        ice = vec![false; nn];
        let mut owner = vec![NONE; nn];
        let mut qin = vec![0.0f64; nn];
        qout = vec![0.0f64; nn];
        let mut stamp = vec![0u32; nn];
        let mut gen = 0u32;
        // Outlet class of each flowline post (resolved receivers-first after the sweep).
        let mut out_cls = vec![TERM_NONE; nn];
        let (mut acc, mut abl_land, mut abl_sea, mut clamp, mut calve_sea, mut calve_lake, mut export, mut export_sea) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut n_calve, mut n_melt_sea, mut n_grounded) = (0u64, 0u64, 0u64);
        termini.clear();
        let mut queue: Vec<usize> = Vec::new();
        let mut hq: std::collections::BinaryHeap<std::cmp::Reverse<(i64, usize)>> = std::collections::BinaryHeap::new();
        for &iu in pop.iter().rev() {
            let i = iu as usize;
            let qi = qin[i];
            let (x, y) = pos(i);
            if kind[i] == 1 && qi > 0.0 {
                // flotation_calving_v1 at a marine post: the face holds here when the calving capacity takes the flux.
                let (w, d) = (width(qi), depth[i]);
                let cap = c_sea * d * w * h_face(d);
                if qi * to_ice <= cap {
                    calve_sea += qi; out_cls[i] = TERM_SEA; n_calve += 1;
                    if qi > 1e7 { termini.push(json!({"post": i, "x_km": x / 1000.0, "y_km": y / 1000.0, "class": "sea", "q_m3_we_yr": qi, "face_depth_m": d, "face_width_m": w, "face_thickness_m": h_face(d), "capacity_ice_m3_yr": cap})); }
                    continue;
                }
                n_grounded += 1;
            }
            let mut q = 0.0;
            if land(i) && bal[i] > 0.0 {
                ice[i] = true; owner[i] = i as u32;
                let p = bal[i] * a_post; acc += p; q = qi + p;
            } else if qi > 0.0 {
                q = qi;
                let w = width(qi);
                let r2 = (0.5 * w).powi(2);
                // Lane GLACIER-3 (gl_geo=1): the valley confinement is evaluated on the bed in every pass (bed of the claimed
                // post <= bed of the flowline post + sqrt(h0 W)); GLACIER-2: on the routing surface, which in passes >= 2 lets
                // the plastic dome over the previous mask lift the confinement level, so tongues creep up their valley walls.
                let zc = if geo.is_some() { bed[i] + (h0 * w).sqrt() } else { surf[i] + (h0 * w).sqrt() };
                let (xi, yi) = lat.g(i);
                gen += 1;
                // Lane GLACIER-3 (gl_geo=1): the claim grows in Euclidean distance order from the flowline post (a partial
                // last claim is a disc, not a Chebyshev square ring); GLACIER-2: breadth-first.
                let use_heap = geo.is_some();
                queue.clear(); hq.clear();
                if use_heap { hq.push(std::cmp::Reverse((0i64, i))); } else { queue.push(i); }
                stamp[i] = gen;
                let mut k = 0;
                loop {
                    let j = if use_heap { match hq.pop() { Some(std::cmp::Reverse((_, j))) => j, None => break } } else { if k >= queue.len() { break; } let j = queue[k]; k += 1; j };
                    if !ice[j] && bal[j] <= 0.0 {
                        let ab = -bal[j] * a_post;
                        ice[j] = true; owner[j] = i as u32; q -= ab;
                        if kind[j] == 1 { abl_sea += ab; } else { abl_land += ab; }
                        if q <= 0.0 { clamp += -q; q = 0.0; break; }
                    }
                    for kk in 0..8 {
                        let Some(m) = nb(j, kk) else { continue };
                        if stamp[m] == gen || !dom[m] { continue; }
                        let (xm, ym) = lat.g(m);
                        let d2 = (((xm - xi) * (xm - xi) + (ym - yi) * (ym - yi)) as f64) * dx * dx;
                        if d2 > r2 || (if geo.is_some() { bed[m] } else { surf[m] }) > zc { continue; }
                        // Lane GLACIER-3 ablation_claim_aligned_v1: no marine post ahead of the flowline post.
                        if geo.is_some() && kind[m] == 1 && ((xm - xi) as f64 * uxv[i] as f64 + (ym - yi) as f64 * uyv[i] as f64) > 1e-9 { continue; }
                        stamp[m] = gen;
                        if use_heap { hq.push(std::cmp::Reverse(((xm - xi) * (xm - xi) + (ym - yi) * (ym - yi), m))); } else { queue.push(m); }
                    }
                }
                if q <= 0.0 { out_cls[i] = if kind[i] == 1 { n_melt_sea += 1; TERM_SEA } else { TERM_LAND }; }
            }
            qout[i] = q;
            if q <= 0.0 { continue; }
            let r = rcv[i];
            if geo.is_some() && r != NONE && dom[r as usize] && r2v[i] != NONE && kind[i] != 1 && kind[r as usize] != 1 && kind[r2v[i] as usize] != 1 {
                // Lane GLACIER-3 balance_route_dinf_v1: the D-infinity split (both receivers are domain posts of lower F).
                let fa = f1v[i] as f64;
                qin[r as usize] += q * fa; qin[r2v[i] as usize] += q * (1.0 - fa); continue;
            }
            if r != NONE && dom[r as usize] { qin[r as usize] += q; continue; }
            if r == NONE || kind[r as usize] == 3 {
                export += q; if kind[i] == 1 { export_sea += q; } out_cls[i] = TERM_EDGE;
                if q > 1e7 { termini.push(json!({"post": i, "x_km": x / 1000.0, "y_km": y / 1000.0, "class": "edge", "q_m3_we_yr": q, "over_sea": kind[i] == 1})); }
            } else {
                calve_lake += q; out_cls[i] = TERM_LAKE;
                if q > 1e7 { termini.push(json!({"post": i, "x_km": x / 1000.0, "y_km": y / 1000.0, "class": "lake", "q_m3_we_yr": q, "note": "lake floor not constructed: calving law not evaluated"})); }
            }
        }
        // Terminus class and terminus post per flowline post, receivers first (pop order): own outlet, else the receiver's.
        let mut out_post = vec![NONE; nn];
        for &iu in pop.iter() {
            let i = iu as usize;
            if out_cls[i] != TERM_NONE { out_post[i] = i as u32; continue; }
            if qout[i] > 0.0 { let r = rcv[i]; if r != NONE && dom[r as usize] { out_cls[i] = out_cls[r as usize]; out_post[i] = out_post[r as usize]; } }
            if out_post[i] == NONE { out_post[i] = i as u32; }
        }
        for i in 0..nn { term_cls[i] = if ice[i] && owner[i] != NONE { let c = out_cls[owner[i] as usize]; if c == TERM_NONE { TERM_LAND } else { c } } else { TERM_NONE }; }
        for i in 0..nn { glac[i] = if ice[i] && owner[i] != NONE { out_post[owner[i] as usize] } else { NONE }; }
        // Plastic surface over the mask.
        let near_sea_depth = |i: usize| -> f64 { let mut d = f64::NAN; for k in 0..8 { if let Some(j) = nb(i, k) { if kind[j] == 1 { d = if d.is_nan() { depth[j] } else { d.max(depth[j]) }; } } } d };
        let seed = |i: usize| -> f64 {
            if kind[i] == 1 { return bed[i] + h_face(depth[i]); }
            // (land margin of a tidewater glacier next to the sea: the face stands at flotation plus q)
            if term_cls[i] == TERM_SEA { let d = near_sea_depth(i); if d.is_finite() { return bed[i].max(sea + h_face(d) - d); } }
            bed[i]
        };
        let (z, n_margin) = if geo.is_some() { plastic_march_fmm(lat, &bed, &ice, h0, &seed) } else { plastic_march_seeded(lat, &bed, &ice, h0, &seed) };
        zs = z;
        let n_ice = ice.iter().filter(|v| **v).count();
        let changed = if prev_mask.is_empty() { n_ice } else { (0..nn).filter(|&i| ice[i] != prev_mask[i]).count() };
        let prev_log_mask = if geo.is_some() { prev_mask.clone() } else { Vec::new() };
        let res = acc - abl_land - abl_sea + clamp - calve_sea - calve_lake - export;
        mass = json!({"accumulation_m3_we_yr": acc, "ablation_land_m3_we_yr": abl_land, "ablation_grounded_on_marine_posts_m3_we_yr": abl_sea, "clamp_unused_ablation_m3_we_yr": clamp,
            "net_ablation_m3_we_yr": abl_land + abl_sea - clamp,
            "calving_sea_m3_we_yr": calve_sea, "calving_lake_m3_we_yr": calve_lake, "export_chart_edge_m3_we_yr": export, "export_chart_edge_over_sea_m3_we_yr": export_sea,
            "residual_m3_we_yr": res, "relative_residual": if acc > 0.0 { res.abs() / acc } else { 0.0 },
            "unrouted_m3_we_yr": 0.0, "export_over_accumulation": if acc > 0.0 { export / acc } else { 0.0 }, "calving_over_accumulation": if acc > 0.0 { (calve_sea + calve_lake) / acc } else { 0.0 },
            "calving_faces": n_calve, "marine_flowline_posts_grounded": n_grounded, "termini_melted_on_marine_posts": n_melt_sea});
        pass_log.push(json!({"pass": pass + 1, "ice_posts": n_ice, "mask_posts_changed_vs_previous_pass": changed, "margin_posts": n_margin}));
        // Lane GLACIER-3 measurement dump (env VZ_GL3_DUMP=<dir>; no output bytes): per-pass ice mask, terminus class, flux.
        if geo.is_some() { if let Ok(dd) = std::env::var("VZ_GL3_DUMP") {
            let _ = std::fs::create_dir_all(&dd);
            let _ = std::fs::write(format!("{dd}/pass{:02}_ice.u8", pass + 1), ice.iter().map(|v| *v as u8).collect::<Vec<u8>>());
            let _ = std::fs::write(format!("{dd}/pass{:02}_term.u8", pass + 1), term_cls.clone());
            let _ = std::fs::write(format!("{dd}/pass{:02}_zs.f32", pass + 1), zs.iter().flat_map(|v| (*v as f32).to_le_bytes()).collect::<Vec<u8>>());
            let _ = std::fs::write(format!("{dd}/meta.json"), format!("{{\"n\": {}, \"gx0\": {}, \"gy0\": {}, \"dx_mm\": {}}}", lat.n, lat.gx0, lat.gy0, lat.dx_mm));
        } }
        prev_mask = ice.clone();
        if let Some(gq) = &geo {
            // Lane GLACIER-3 passes_to_tolerance_v1: stop at the first pass >= 2 within the declared mask-change tolerance.
            // Changed posts touching ice that is ice in both passes (a one-post fringe of a persistent glacier) vs patches.
            let (mut fringe, mut patch) = (0u64, 0u64);
            if !prev_log_mask.is_empty() {
                for i in 0..nn {
                    if ice[i] == prev_log_mask[i] { continue; }
                    if (0..8).any(|k| nb(i, k).map(|j| ice[j] && prev_log_mask[j]).unwrap_or(false)) { fringe += 1; } else { patch += 1; }
                }
            }
            if let Some(o) = pass_log.last_mut().and_then(|v| v.as_object_mut()) { o.insert("changed_ppm_of_ice".into(), json!(changed as f64 * 1e6 / n_ice.max(1) as f64)); o.insert("changed_fringe_posts".into(), json!(fringe)); o.insert("changed_patch_posts".into(), json!(patch)); }
            let dn = (n_ice as f64 - prev_log_mask.iter().filter(|v| **v).count() as f64).abs();
            if pass >= 1 && (patch as f64) * 1e6 <= gq.patch_change_tol_ppm as f64 * n_ice as f64 && dn * 1e6 <= gq.area_change_tol_ppm as f64 * n_ice as f64 { converged_at = Some(pass + 1); break; }
            if pass + 1 < passes {
                let rl = gq.relax_permille as f64 / 1000.0;
                surf = (0..nn).map(|i| if ice[i] && kind[i] == 0 && owner[i] == i as u32 && zs[i].is_finite() { if rl >= 1.0 { zs[i] } else { rl * zs[i] + (1.0 - rl) * surf[i] } } else { bed[i] }).collect();
            }
            continue;
        }
        if pass + 1 < passes { surf = (0..nn).map(|i| if ice[i] && dom[i] && zs[i].is_finite() { zs[i] } else { bed[i] }).collect(); }
    }
    // Statistics: AAR per terminus class (balance surface = the last pass's routing surface; ice surface = final zs).
    let surf_fin = &surf;
    let mut cls_stat = vec![[0u64; 4]; 5]; // [ice posts, posts with routing-surface balance > 0, posts with ice surface >= ELA, glaciers]
    let mut hs: Vec<f64> = Vec::new();
    for i in 0..nn {
        if !ice[i] { continue; }
        let c = term_cls[i] as usize;
        cls_stat[c][0] += 1;
        let zb = if land(i) { surf_fin[i] } else { sea };
        if ice_share(zb, ela[i], g.full_m) > 0.0 { cls_stat[c][1] += 1; }
        if zs[i] >= ela[i] { cls_stat[c][2] += 1; }
        hs.push((zs[i] - bed[i]).max(0.0));
    }
    hs.sort_by(|a, c| a.total_cmp(c));
    let pct = |p: f64| -> f64 { if hs.is_empty() { 0.0 } else { hs[((hs.len() - 1) as f64 * p).round() as usize] } };
    let km2 = a_post / 1e6;
    let names = ["none", "land", "sea", "lake", "edge"];
    // Per glacier (terminus post): ice posts and posts with ice surface >= ELA.
    let mut per: std::collections::BTreeMap<u32, (u8, u64, u64)> = std::collections::BTreeMap::new();
    for i in 0..nn { if ice[i] && glac[i] != NONE { let e = per.entry(glac[i]).or_insert((term_cls[i], 0, 0)); e.1 += 1; if zs[i] >= ela[i] { e.2 += 1; } } }
    let mut by_class = serde_json::Map::new();
    for c in 1..5 {
        let s = cls_stat[c];
        let mut aars: Vec<f64> = per.values().filter(|v| v.0 as usize == c && v.1 as f64 * km2 >= 1.0).map(|v| v.2 as f64 / v.1 as f64).collect();
        aars.sort_by(|a, b| a.total_cmp(b));
        let q = |p: f64| -> f64 { if aars.is_empty() { 0.0 } else { aars[((aars.len() - 1) as f64 * p).round() as usize] } };
        let n_all = per.values().filter(|v| v.0 as usize == c).count();
        by_class.insert(names[c].into(), json!({"ice_area_km2": s[0] as f64 * km2, "aar_balance_surface": if s[0] > 0 { s[1] as f64 / s[0] as f64 } else { 0.0 }, "aar_ice_surface": if s[0] > 0 { s[2] as f64 / s[0] as f64 } else { 0.0 },
            "glaciers": n_all, "glaciers_ge_1km2": aars.len(), "per_glacier_aar_ice_surface_ge_1km2": {"p10": q(0.1), "p25": q(0.25), "p50": q(0.5), "p75": q(0.75), "p90": q(0.9), "share_in_0.5_0.7": if aars.is_empty() { 0.0 } else { aars.iter().filter(|a| **a >= 0.5 && **a <= 0.7).count() as f64 / aars.len() as f64 }}}));
    }
    termini.sort_by(|a, c| c["q_m3_we_yr"].as_f64().unwrap_or(0.0).total_cmp(&a["q_m3_we_yr"].as_f64().unwrap_or(0.0)));
    termini.truncate(40);
    // Crop to the mask bounding box + 2 posts.
    let (mut ix0, mut ix1, mut iy0, mut iy1) = (usize::MAX, 0usize, usize::MAX, 0usize);
    let mut n_mask = 0u64;
    for i in 0..nn { if ice[i] { n_mask += 1; let (x, y) = (i % n, i / n); ix0 = ix0.min(x); ix1 = ix1.max(x); iy0 = iy0.min(y); iy1 = iy1.max(y); } }
    if n_mask == 0 { ix0 = 0; ix1 = 0; iy0 = 0; iy1 = 0; }
    let (ix0, iy0) = (ix0.saturating_sub(2), iy0.saturating_sub(2));
    let (ix1, iy1) = ((ix1 + 2).min(n - 1), (iy1 + 2).min(n - 1));
    let (nx, ny) = (ix1 - ix0 + 1, iy1 - iy0 + 1);
    let (mut szs, mut sh, mut sq, mut st) = (Vec::with_capacity(nx * ny), Vec::with_capacity(nx * ny), Vec::with_capacity(nx * ny), Vec::with_capacity(nx * ny));
    for y in iy0..=iy1 { for x in ix0..=ix1 {
        let i = y * n + x;
        if ice[i] { szs.push(zs[i] as f32); sh.push((zs[i] - bed[i]).max(0.0) as f32); } else { szs.push(if land(i) { bd[i] as f32 } else { f32::NAN }); sh.push(-1.0); }
        sq.push(if ice[i] { qout[i] as f32 } else { 0.0 });
        st.push(term_cls[i]);
    } }
    nw.ela_grid = ela_grid; nw.p_grid = p_grid;
    nw.dx_mm = lat.dx_mm; nw.gx0 = lat.gx0 + ix0 as i64; nw.gy0 = lat.gy0 + iy0 as i64; nw.nx = nx; nw.ny = ny;
    nw.sheet_zs = szs; nw.sheet_h = sh; nw.sheet_q = sq; nw.sheet_term = st;
    let n_marine_ice = (0..nn).filter(|&i| ice[i] && kind[i] == 1).count();
    // Lane GLACIER-3 (gl_geo=1 only): convergence and per-glacier volume-area scaling (Bahr, Meier and Peckham 1997 JGR
    // 102:20355; worldwide mean c = 0.034 km^(3 - 2 gamma), gamma = 1.375 for glaciers, Bahr, Pfeffer and Kaser 2015 Rev
    // Geophys 53:95, who stress that c is a random variable): ratio V / (c A^gamma) over glaciers >= 1 km2 (ice volume over
    // the sheet bed, km3; marine posts counted from the sea floor).
    let geo_stats = geo.as_ref().map(|gq| {
        let mut vol: std::collections::BTreeMap<u32, (u8, u64, f64)> = std::collections::BTreeMap::new();
        for i in 0..nn { if ice[i] && glac[i] != NONE { let e = vol.entry(glac[i]).or_insert((term_cls[i], 0, 0.0)); e.1 += 1; e.2 += (zs[i] - bed[i]).max(0.0) * a_post; } }
        let mut va = serde_json::Map::new();
        for (c, nm) in [(1u8, "land"), (2, "sea"), (4, "edge")] {
            let mut r: Vec<f64> = vol.values().filter(|v| v.0 == c && v.1 as f64 * km2 >= 1.0).map(|v| { let ak = v.1 as f64 * km2; (v.2 / 1e9) / (0.034 * ak.powf(1.375)) }).collect();
            r.sort_by(|a, b| a.total_cmp(b));
            let q = |p: f64| -> f64 { if r.is_empty() { 0.0 } else { r[((r.len() - 1) as f64 * p).round() as usize] } };
            let mut big: Vec<(f64, f64, f64)> = vol.values().filter(|v| v.0 == c && v.1 as f64 * km2 >= 1.0).map(|v| { let ak = v.1 as f64 * km2; (ak, v.2 / 1e9, 1000.0 * v.2 / 1e9 / ak) }).collect();
            big.sort_by(|a, b| b.0.total_cmp(&a.0));
            big.truncate(8);
            va.insert(nm.into(), json!({"glaciers_ge_1km2": r.len(), "ratio_v_over_c_a_gamma": {"p10": q(0.1), "p25": q(0.25), "p50": q(0.5), "p75": q(0.75), "p90": q(0.9)},
                "largest_area_km2_volume_km3_mean_thickness_m": big}));
        }
        json!({"rule": "lane GLACIER-3 gl_geo=1: plastic_sheet_fmm_v1 + balance_route_dinf_v1 + ablation_claim_aligned_v1 + passes_to_tolerance_v1", "max_passes": gq.max_passes, "patch_change_tol_ppm": gq.patch_change_tol_ppm, "area_change_tol_ppm": gq.area_change_tol_ppm,
            "relax_permille": gq.relax_permille, "converged": converged_at.is_some(), "converged_at_pass": converged_at, "volume_area_scaling": va})
    });
    if let Some(o) = nw.stats.as_object_mut() {
        o.insert("sheet".into(), json!({"rule": "lane GLACIER-2 balance_route_v1 + ablation_claim_v1 + flotation_calving_v1 + plastic_sheet_seeded_v1 (gl_bal=1); the per-node gl_*_now fields above are GLACIER-1's tree state and do not drive this sheet",
            "h0_m": h0, "despiked_bed_posts": n_spike, "passes": pass_log, "mask_posts": n_mask, "mask_area_km2": n_mask as f64 * km2, "marine_ice_area_km2": n_marine_ice as f64 * km2,
            "by_terminus_class": by_class, "mass": mass, "termini_top": termini,
            "h_p10_m": pct(0.1), "h_p50_m": pct(0.5), "h_p90_m": pct(0.9), "h_p99_m": pct(0.99), "h_max_m": pct(1.0), "seconds_not_digested": 0.0}));
        if let Some(gs) = geo_stats { if let Some(sh) = o.get_mut("sheet").and_then(|v| v.as_object_mut()) { sh.insert("geo".into(), gs); } }
    }
    Ok(())
}

/// Lane GLACIER-2 (gl_bal=1) `ice_window_v2`: per window post (chart position, quantized height mm, water code 0 land /
/// 1 marine / 2 lake) the published ice fields (ice_now_mm, ice_surface_mm, ice_class).
pub fn window_ice_v2(now: &IceNow, nn: usize, xy: &dyn Fn(usize) -> (f64, f64), height_mm: &[i32], water: &dyn Fn(usize) -> u8) -> (Vec<i32>, Vec<i32>, Vec<u8>) {
    let (mut ice, mut surf, mut cls) = (vec![0i32; nn], height_mm.to_vec(), vec![0u8; nn]);
    let d = now.dx_mm as f64 / 1000.0;
    for i in 0..nn {
        let wc = water(i);
        if wc == 2 { continue; }
        let (x, y) = xy(i);
        let (fx, fy) = (x / d - now.gx0 as f64, y / d - now.gy0 as f64);
        if now.nx < 2 || now.ny < 2 || !(fx >= 0.0 && fy >= 0.0 && fx <= (now.nx - 1) as f64 && fy <= (now.ny - 1) as f64) { continue; }
        let (x0, y0) = ((fx.floor() as usize).min(now.nx - 2), (fy.floor() as usize).min(now.ny - 2));
        let (tx, ty) = (fx - x0 as f64, fy - y0 as f64);
        let (mut nm, mut interior, mut z, mut hh, mut wm, mut zmin, mut zmax) = (0u32, true, 0.0f64, 0.0f64, 0.0f64, f64::INFINITY, f64::NEG_INFINITY);
        for (dy, wy) in [(0usize, 1.0 - ty), (1, ty)] { for (dx, wx) in [(0usize, 1.0 - tx), (1, tx)] {
            let k = (y0 + dy) * now.nx + x0 + dx;
            let (h, zz) = (now.sheet_h[k], now.sheet_zs[k]);
            if !(h >= 0.0) || !zz.is_finite() { continue; }
            nm += 1;
            if !(h > 0.0) { interior = false; }
            z += wx * wy * zz as f64; hh += wx * wy * h as f64; wm += wx * wy; zmin = zmin.min(zz as f64); zmax = zmax.max(zz as f64);
        } }
        // 4 mask corners: the cell is ice-covered; 3 (a concave margin corner or a one-post enclave of the 128 m mask): the
        // mask corners' surface extrapolated over the cell, contour rule; fewer: no ice.
        if nm < 3 || (nm == 3 && wc == 1) || !(wm > 0.0) { continue; }
        if nm == 3 { interior = false; z /= wm; hh /= wm; }
        let zp = height_mm[i] as f64 / 1000.0;
        let h = if wc == 1 { z - zp } else if interior {
            let zs = if zmax - zmin > 300.0 { zmin } else { z };
            (zs - zp).min(hh + 200.0).max(hh)
        } else { (z - zp).min(hh + 200.0) };
        if !(h > 0.0) { continue; }
        let hm = vz_tilekey::chart::quantize_mm(h * 1000.0).clamp(0, i32::MAX as i64) as i32;
        if hm <= 0 { continue; }
        ice[i] = hm;
        surf[i] = height_mm[i].saturating_add(hm);
        let ela = node_grid_bilinear(&now.ela_grid, x, y);
        cls[i] = if surf[i] as f64 / 1000.0 >= ela { 1 } else { 2 };
    }
    (ice, surf, cls)
}

/// Lane GLACIER-3 (gl_geo=1) `ice_window_v3`: ice_window_v2 with a continuous cliff rule. v2 switched, per 128 m cell, from the
/// bilinear surface to the lowest corner surface where the cell's corner surfaces span more than 300 m, so neighbouring cells
/// disagreed at their shared edge (ice-thickness steps > 100 m on 128 m cell edges, measured: 1,771 of 2,663 such steps in the
/// site 5 + ls_undercut 3x3 straddle a cliff cell; the x ~ -14 km band carries 7.2 % cliff cells vs 4.4 % elsewhere). v3: every
/// sheet post k carries a fill allowance a_k = 200 m where its surface differs from all its finite 8-neighbours by at most
/// 300 m, else 0 (a 128 m cliff); a window post takes the bilinear surface z_bil and ice = clamp(z_bil - z_post, H_bil,
/// H_bil + a_bil) on interior cells and min(z_bil - z_post, H_bil + a_bil) on margin / 3-corner cells, a_bil the bilinear of the
/// corner allowances: no fill across a cliff (drape of the sheet thickness only), continuous across every cell edge.
pub fn window_ice_v3(now: &IceNow, nn: usize, xy: &dyn Fn(usize) -> (f64, f64), height_mm: &[i32], water: &dyn Fn(usize) -> u8) -> (Vec<i32>, Vec<i32>, Vec<u8>) {
    let (mut ice, mut surf, mut cls) = (vec![0i32; nn], height_mm.to_vec(), vec![0u8; nn]);
    let d = now.dx_mm as f64 / 1000.0;
    let allow = |x: usize, y: usize| -> f64 {
        let z0 = now.sheet_zs[y * now.nx + x];
        if !z0.is_finite() { return 0.0; }
        for (ox, oy) in crate::D8 {
            let (xx, yy) = (x as i64 + ox, y as i64 + oy);
            if xx < 0 || yy < 0 || xx >= now.nx as i64 || yy >= now.ny as i64 { continue; }
            let z1 = now.sheet_zs[yy as usize * now.nx + xx as usize];
            if z1.is_finite() && (z1 as f64 - z0 as f64).abs() > 300.0 { return 0.0; }
        }
        200.0
    };
    for i in 0..nn {
        let wc = water(i);
        if wc == 2 { continue; }
        let (x, y) = xy(i);
        let (fx, fy) = (x / d - now.gx0 as f64, y / d - now.gy0 as f64);
        if now.nx < 2 || now.ny < 2 || !(fx >= 0.0 && fy >= 0.0 && fx <= (now.nx - 1) as f64 && fy <= (now.ny - 1) as f64) { continue; }
        let (x0, y0) = ((fx.floor() as usize).min(now.nx - 2), (fy.floor() as usize).min(now.ny - 2));
        let (tx, ty) = (fx - x0 as f64, fy - y0 as f64);
        let (mut nm, mut interior, mut z, mut hh, mut wm, mut al) = (0u32, true, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for (dy, wy) in [(0usize, 1.0 - ty), (1, ty)] { for (dx, wx) in [(0usize, 1.0 - tx), (1, tx)] {
            let k = (y0 + dy) * now.nx + x0 + dx;
            let (h, zz) = (now.sheet_h[k], now.sheet_zs[k]);
            if !(h >= 0.0) || !zz.is_finite() { continue; }
            nm += 1;
            if !(h > 0.0) { interior = false; }
            z += wx * wy * zz as f64; hh += wx * wy * h as f64; wm += wx * wy; al += wx * wy * allow(x0 + dx, y0 + dy);
        } }
        if nm < 3 || (nm == 3 && wc == 1) || !(wm > 0.0) { continue; }
        if nm == 3 { interior = false; z /= wm; hh /= wm; al /= wm; }
        let zp = height_mm[i] as f64 / 1000.0;
        let h = if wc == 1 { z - zp } else if interior { (z - zp).min(hh + al).max(hh) } else { (z - zp).min(hh + al) };
        if !(h > 0.0) { continue; }
        let hm = vz_tilekey::chart::quantize_mm(h * 1000.0).clamp(0, i32::MAX as i64) as i32;
        if hm <= 0 { continue; }
        ice[i] = hm;
        surf[i] = height_mm[i].saturating_add(hm);
        let ela = node_grid_bilinear(&now.ela_grid, x, y);
        cls[i] = if surf[i] as f64 / 1000.0 >= ela { 1 } else { 2 };
    }
    (ice, surf, cls)
}

/// Lane GLACIER-3 window field note (gl_geo=1 only).
pub const ICE_SUB_NOTE: &str = "lane GLACIER-3 (T1b, gl_geo=1): submerged thickness of a grounded tidewater tongue on MARINE window posts with ice_now_mm > 0, mm: published height (sea level) minus the sheet bed (the sea floor under the ice, bilinear of the layer sheet's z_s - H over the mask corners), >= 0; 0 elsewhere; full ice thickness there = ice_now_mm + ice_sub_mm; derived, never canon";

/// Lane GLACIER-3 (gl_geo=1): `ice_sub_mm` per window post from the published sheet and the window's ice_now_mm.
pub fn window_ice_sub(now: &IceNow, nn: usize, xy: &dyn Fn(usize) -> (f64, f64), height_mm: &[i32], ice_now_mm: &[i32], marine: &dyn Fn(usize) -> bool) -> Vec<i32> {
    let mut out = vec![0i32; nn];
    let d = now.dx_mm as f64 / 1000.0;
    for i in 0..nn {
        if ice_now_mm[i] <= 0 || !marine(i) { continue; }
        let (x, y) = xy(i);
        let (fx, fy) = (x / d - now.gx0 as f64, y / d - now.gy0 as f64);
        if now.nx < 2 || now.ny < 2 || !(fx >= 0.0 && fy >= 0.0 && fx <= (now.nx - 1) as f64 && fy <= (now.ny - 1) as f64) { continue; }
        let (x0, y0) = ((fx.floor() as usize).min(now.nx - 2), (fy.floor() as usize).min(now.ny - 2));
        let (tx, ty) = (fx - x0 as f64, fy - y0 as f64);
        let (mut bsum, mut w) = (0.0f64, 0.0f64);
        for (dy, wy) in [(0usize, 1.0 - ty), (1, ty)] { for (dx, wx) in [(0usize, 1.0 - tx), (1, tx)] {
            let k = (y0 + dy) * now.nx + x0 + dx;
            let (h, zz) = (now.sheet_h[k], now.sheet_zs[k]);
            if !(h >= 0.0) || !zz.is_finite() { continue; }
            bsum += wx * wy * (zz as f64 - h as f64); w += wx * wy;
        } }
        if !(w > 0.0) { continue; }
        let sub = height_mm[i] as f64 / 1000.0 - bsum / w;
        if sub > 0.0 { out[i] = vz_tilekey::chart::quantize_mm(sub * 1000.0).clamp(0, i32::MAX as i64) as i32; }
    }
    out
}

/// Window / solve-lattice statistics.
#[derive(Clone, Debug, Default, Serialize)]
pub struct GlStats { pub links_in_lattice: u64, pub section_lowered_posts: u64, pub section_lowered_max_m: f64, pub section_lowered_mean_m: f64, pub cirques_in_lattice: u64, pub cirque_lowered_posts: u64, pub cirque_lowered_max_m: f64, pub cirque_channel_posts: u64 }

fn seg_dist(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> (f64, f64) {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0) } else { 0.0 };
    ((p[0] - a[0] - t * dx).hypot(p[1] - a[1] - t * dy), t)
}

/// The window carve (pointwise reads of the table): U sections on free/floor posts, cirque bowls. Returns per-post lowering.
#[allow(clippy::too_many_arguments)]
pub fn carve(lat: &EvalLat, z: &mut [f64], fixed: &[u8], water: &[u8], net: &Network, curves: &[crate::t1b_rivernet::Curve], t: &GlacialTable, gp: &GlacialParamsV1) -> (GlStats, Vec<f32>) {
    let (n, dx) = (lat.n, lat.dx());
    let nn = lat.nodes();
    let mut st = GlStats::default();
    let mut zg = vec![f64::INFINITY; nn];
    // Network channel posts (GUIDE) read the section on the FILLED floor (monotone along flow), free posts the bedrock floor.
    let mut zgf = vec![f64::INFINITY; nn];
    let b = gp.u_power_permille as f64 / 1000.0;
    let reach = gp.u_reach_permille_of_width as f64 / 1000.0;
    let (lx0, ly0) = (lat.gx0 as f64 * dx, lat.gy0 as f64 * dx);
    let (lx1, ly1) = ((lat.gx0 + n as i64 - 1) as f64 * dx, (lat.gy0 + n as i64 - 1) as f64 * dx);
    // T1a links are drawn on their own curves (Catmull-Rom or the T1a planform, cl): the section follows the drawn curve,
    // with the curve's bed as its floor; land links follow node -> meander midpoint -> receiver (Chaikin cuts corners by at
    // most a quarter link).
    let mut t1a_curve: std::collections::BTreeMap<u32, usize> = std::collections::BTreeMap::new();
    for (ci, c) in curves.iter().enumerate() { if c.code == CHANNEL && c.link != u32::MAX && (c.link as usize) < NODES && net.kind[c.link as usize] == K_T1A { t1a_curve.entry(c.link).or_insert(ci); } }
    for k in 0..NODES {
        if t.w[k] <= 0.0 || t.rcv[k] == NONE { continue; }
        let r = t.rcv[k] as usize;
        let (w0, h0) = (t.w[k] as f64, t.h[k] as f64);
        let (w1, h1) = if t.w[r] > 0.0 { (t.w[r] as f64, t.h[r] as f64) } else { (0.0, 0.0) };
        // Polyline points (x, y, floor).
        let poly: Vec<[f64; 4]> = match t1a_curve.get(&(k as u32)) {
            Some(&ci) if net.kind[k] == K_T1A => curves[ci].pts.iter().map(|q| [q[0], q[1], q[2], q[2]]).collect(),
            _ => {
                let pa = net.pos[k];
                let pb = net.pos[r];
                let pm = if net.kind[k] == K_LAND { net.mid[k] } else { [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])] };
                let (f0, f1) = (t.bed[k] as f64, t.bed[r] as f64);
                let (g0, g1) = (t.bed_fill[k] as f64, t.bed_fill[r] as f64);
                vec![[pa[0], pa[1], f0, g0], [pm[0], pm[1], 0.5 * (f0 + f1), 0.5 * (g0 + g1)], [pb[0], pb[1], f1, g1]]
            }
        };
        if poly.len() < 2 { continue; }
        let rad = reach * w0.max(w1);
        let (mut bx0, mut bx1, mut by0, mut by1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for q in &poly { bx0 = bx0.min(q[0]); bx1 = bx1.max(q[0]); by0 = by0.min(q[1]); by1 = by1.max(q[1]); }
        let (bx0, bx1, by0, by1) = (bx0 - rad, bx1 + rad, by0 - rad, by1 + rad);
        if bx1 < lx0 || bx0 > lx1 || by1 < ly0 || by0 > ly1 { continue; }
        st.links_in_lattice += 1;
        let mut cum = vec![0.0f64; poly.len()];
        for j in 1..poly.len() { cum[j] = cum[j - 1] + (poly[j][0] - poly[j - 1][0]).hypot(poly[j][1] - poly[j - 1][1]); }
        let ltot = cum[poly.len() - 1].max(1e-9);
        let gx_lo = ((bx0 / dx).floor() as i64).max(lat.gx0);
        let gx_hi = ((bx1 / dx).ceil() as i64).min(lat.gx0 + n as i64 - 1);
        let gy_lo = ((by0 / dx).floor() as i64).max(lat.gy0);
        let gy_hi = ((by1 / dx).ceil() as i64).min(lat.gy0 + n as i64 - 1);
        for gy in gy_lo..=gy_hi {
            for gx in gx_lo..=gx_hi {
                let i = ((gy - lat.gy0) as usize) * n + (gx - lat.gx0) as usize;
                if water[i] != FREE || !(fixed[i] == FREE || fixed[i] == GUIDE) { continue; }
                let p = [gx as f64 * dx, gy as f64 * dx];
                let (mut y, mut s, mut fl, mut ff) = (f64::INFINITY, 0.0, 0.0, 0.0);
                for j in 1..poly.len() {
                    let (a, bq) = (poly[j - 1], poly[j]);
                    let (dd, tt) = seg_dist(p, [a[0], a[1]], [bq[0], bq[1]]);
                    if dd < y { y = dd; s = (cum[j - 1] + tt * (cum[j] - cum[j - 1])) / ltot; fl = a[2] + tt * (bq[2] - a[2]); ff = a[3] + tt * (bq[3] - a[3]); }
                }
                let w = w0 + (w1 - w0) * s;
                if w <= dx || y > reach * w { continue; }
                let h = h0 + (h1 - h0) * s;
                // Continuous at the reach edge: the rise grows without bound there (rho / (1 - rho^2 / rho_max^2)).
                let rho = 2.0 * y / w;
                let rmax = 2.0 * reach;
                let den = 1.0 - (rho / rmax).powi(2);
                if den <= 1e-9 { continue; }
                let sh = h * rho.powf(b) / den;
                if fixed[i] == GUIDE { if ff + sh < zgf[i] { zgf[i] = ff + sh; } } else if fl + sh < zg[i] { zg[i] = fl + sh; }
            }
        }
    }
    let mut low = vec![0f32; nn];
    let mut sum = 0.0f64;
    for i in 0..nn {
        if zgf[i] < zg[i] { zg[i] = zgf[i]; }
        if zg[i].is_finite() && z[i].is_finite() && zg[i] < z[i] {
            let d = z[i] - zg[i];
            st.section_lowered_posts += 1;
            st.section_lowered_max_m = st.section_lowered_max_m.max(d);
            sum += d;
            low[i] = d as f32;
            z[i] = zg[i];
        }
    }
    if st.section_lowered_posts > 0 { st.section_lowered_mean_m = sum / st.section_lowered_posts as f64; }
    // Cirque bowls.
    for c in &t.cirques {
        let rad = c.radius_m * 3.0 * std::f64::consts::SQRT_2;
        if c.x_m + rad < lx0 || c.x_m - rad > lx1 || c.y_m + rad < ly0 || c.y_m - rad > ly1 { continue; }
        st.cirques_in_lattice += 1;
        let gx_lo = (((c.x_m - rad) / dx).floor() as i64).max(lat.gx0);
        let gx_hi = (((c.x_m + rad) / dx).ceil() as i64).min(lat.gx0 + n as i64 - 1);
        let gy_lo = (((c.y_m - rad) / dx).floor() as i64).max(lat.gy0);
        let gy_hi = (((c.y_m + rad) / dx).ceil() as i64).min(lat.gy0 + n as i64 - 1);
        for gy in gy_lo..=gy_hi {
            for gx in gx_lo..=gx_hi {
                let i = ((gy - lat.gy0) as usize) * n + (gx - lat.gx0) as usize;
                if water[i] != FREE || !z[i].is_finite() { continue; }
                let (ex, ey) = (gx as f64 * dx - c.x_m, gy as f64 * dx - c.y_m);
                let (u, v) = (ex * c.ux + ey * c.uy, -ex * c.uy + ey * c.ux);
                // Armchair: steep headwall up-glacier, gentle opening toward the lip; continuous at the footprint edge
                // (q / (1 - q / 4) grows without bound there).
                let q = (v / c.radius_m).powi(2) + (u.min(0.0) / c.radius_m).powi(2) + (u.max(0.0) / (2.0 * c.radius_m)).powi(2);
                if q >= 4.0 - 1e-9 { continue; }
                let bowl = c.floor_m + c.depth_m * q / (1.0 - q / 4.0);
                let zt = if fixed[i] == CHANNEL || fixed[i] == GUIDE { st.cirque_channel_posts += 1; bowl.max(c.flat_m) } else if fixed[i] == FREE { bowl } else { continue };
                if zt < z[i] {
                    let d = z[i] - zt;
                    st.cirque_lowered_posts += 1;
                    st.cirque_lowered_max_m = st.cirque_lowered_max_m.max(d);
                    low[i] += d as f32;
                    z[i] = zt;
                }
            }
        }
    }
    (st, low)
}

#[cfg(test)]
mod tests {
    //! SYNTHETIC FIXTURES only (ADR 026 §1: instrument calibration).
    use super::*;

    #[test]
    fn ice_share_matches_the_t0_law() {
        // Same expression as pf_surface::erosion::glacial::ice_share (Hergarten eq 23 / Liebl eq 11).
        assert_eq!(ice_share(1500.0, 1000.0, 500.0), 1.0);
        assert_eq!(ice_share(1250.0, 1000.0, 500.0), 0.5);
        assert!((ice_share(750.0, 1000.0, 500.0) - (-(0.5 + 0.25 / 3.0))).abs() < 1e-12);
        assert_eq!(ice_share(1000.0, 1000.0, 500.0), 0.0);
    }

    #[test]
    fn nye_table_and_width_law() {
        assert_eq!(nye_parabola(0.5), 0.445);
        assert!((nye_parabola(2.0) - 0.646).abs() < 1e-12);
        assert!((nye_parabola(2.5) - 0.696).abs() < 1e-12);
        assert_eq!(nye_parabola(9.0), 0.806);
        // Liebl 2023: 1e4 m3/a -> 150 m; 2e7 -> ~1.47 km.
        let w = |q: f64| 150.0 * (q / 1e4f64).powf(0.3);
        assert!((w(1e4) - 150.0).abs() < 1e-9);
        assert!((w(2e7) - 1468.0).abs() < 5.0);
    }

    #[test]
    fn plastic_march_matches_the_flat_bed_profile() {
        // Lane GLACIER-1 (SYNTHETIC fixture): flat bed, 201 x 201 posts at 128 m.
        let n = 201usize;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 128_000 };
        let nn = lat.nodes();
        let b = vec![0.0f64; nn];
        // Mask: a disc of radius 10 km around the centre; H at distance L from the margin ~ sqrt(2 h0 L) (Nye 1951).
        let c = (n / 2) as f64;
        let mask: Vec<bool> = (0..nn).map(|i| { let (x, y) = ((i % n) as f64, (i / n) as f64); ((x - c).powi(2) + (y - c).powi(2)).sqrt() * 128.0 <= 10_000.0 }).collect();
        let h0 = 100_000.0 / (917.0 * 9.81);
        let (zs, m) = plastic_march(&lat, &b, &mask, h0);
        assert!(m > 0);
        let centre = zs[(n / 2) * n + n / 2];
        let expect = (2.0 * h0 * 10_000.0f64).sqrt();
        assert!(centre > expect * 0.97 && centre < expect * 1.10, "centre {centre} vs {expect}");
        // Off-mask posts stay NaN; margin posts sit on the bed.
        assert!(zs[0].is_nan());
        // A steep bed (slope 0.5) under an existing surface: the carapace H ~ h0 / S.
        let lat2 = EvalLat { gx0: 0, gy0: 0, n: 64, dx_mm: 128_000 };
        let b2: Vec<f64> = (0..lat2.nodes()).map(|i| 64.0 * (i % 64) as f64).collect();
        let mask2: Vec<bool> = (0..lat2.nodes()).map(|i| { let (x, y) = (i % 64, i / 64); x > 0 && x < 63 && y > 0 && y < 63 }).collect();
        let (z2, _) = plastic_march(&lat2, &b2, &mask2, h0);
        let i = 32 * 64 + 50;
        let hh = z2[i] - b2[i];
        assert!(hh > 5.0 && hh < 40.0, "carapace {hh}");
    }

    #[test]
    fn seeded_march_equals_plastic_march_with_bed_seeds() {
        // Lane GLACIER-2: plastic_march_seeded with seed = bed reproduces plastic_march bit for bit.
        let lat = EvalLat { gx0: -20, gy0: -20, n: 41, dx_mm: 128_000 };
        let nn = lat.nodes();
        let b: Vec<f64> = (0..nn).map(|i| { let (gx, gy) = lat.g(i); 100.0 + 3.0 * gx as f64 + 0.5 * (gy * gy) as f64 }).collect();
        let mask: Vec<bool> = (0..nn).map(|i| { let (gx, gy) = lat.g(i); gx * gx + gy * gy < 300 }).collect();
        let (a, na) = plastic_march(&lat, &b, &mask, 11.1);
        let (c, nc) = plastic_march_seeded(&lat, &b, &mask, 11.1, &|i| b[i]);
        assert_eq!(na, nc);
        for i in 0..nn { assert!(a[i].to_bits() == c[i].to_bits() || (a[i].is_nan() && c[i].is_nan())); }
    }

    #[test]
    fn balance_route_is_a_tree_draining_to_outlets() {
        // Lane GLACIER-2: every domain post's receiver chain ends at a lake post (outlet), receivers are popped before their
        // donors (so the reverse pop order is donors-first), and a closed pit drains over its spill.
        let lat = EvalLat { gx0: 0, gy0: 0, n: 30, dx_mm: 128_000 };
        let nn = lat.nodes();
        let mut s: Vec<f64> = (0..nn).map(|i| { let (gx, gy) = lat.g(i); 1000.0 - 10.0 * gx as f64 + (gy as f64 - 15.0).abs() }).collect();
        s[15 * 30 + 15] = 500.0; // a pit
        let mut lake = vec![false; nn];
        for gy in 0..30 { lake[gy * 30 + 29] = true; s[gy * 30 + 29] = 0.0; } // water on the east column
        let dom: Vec<bool> = lake.iter().map(|l| !l).collect();
        let (pop, rcv) = bal_route(&lat, &s, &dom, &lake);
        assert_eq!(pop.len(), nn - 30);
        let mut rank = vec![usize::MAX; nn];
        for (k, &i) in pop.iter().enumerate() { rank[i as usize] = k; }
        for &i in &pop {
            let r = rcv[i as usize];
            assert!(r != NONE, "every domain post drains to a domain post or the water");
            if dom[r as usize] { assert!(rank[r as usize] < rank[i as usize]); }
            let mut j = i as usize; let mut steps = 0;
            while dom[j] { j = rcv[j] as usize; steps += 1; assert!(steps < nn); }
            assert!(lake[j]);
        }
    }

    #[test]
    fn fmm_sheet_is_isotropic_and_matches_the_flat_bed_profile() {
        // Lane GLACIER-3 (SYNTHETIC fixture): flat bed, disc mask of radius 10 km at 128 m. The plastic surface at distance L
        // from the margin is sqrt(2 h0 L) (Nye 1951); FMM fronts are Euclidean, so posts at the same radius along an axis and
        // along a diagonal agree within a few per cent (Dijkstra on 8 neighbours: up to ~4 % in H).
        let n = 201usize;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 128_000 };
        let nn = lat.nodes();
        let b = vec![0.0f64; nn];
        let c = (n / 2) as f64;
        let mask: Vec<bool> = (0..nn).map(|i| { let (x, y) = ((i % n) as f64, (i / n) as f64); ((x - c).powi(2) + (y - c).powi(2)).sqrt() * 128.0 <= 10_000.0 }).collect();
        let h0 = 100_000.0 / (917.0 * 9.81);
        let (zs, m) = plastic_march_fmm(&lat, &b, &mask, h0, &|i| b[i]);
        assert!(m > 0);
        let centre = zs[(n / 2) * n + n / 2];
        let expect = (2.0 * h0 * 10_000.0f64).sqrt();
        assert!(centre > expect * 0.97 && centre < expect * 1.05, "centre {centre} vs {expect}");
        // Radius 5 km: axis post (c + 39, c) at 4992 m and diagonal post (c + 28, c + 28) at 5068.5 m from the centre.
        let ax = zs[(n / 2) * n + n / 2 + 39];
        let dg = zs[(n / 2 + 28) * n + n / 2 + 28];
        let ea = (2.0 * h0 * (10_000.0 - 4992.0f64)).sqrt();
        let ed = (2.0 * h0 * (10_000.0 - 5068.5f64)).sqrt();
        assert!((ax / ea - 1.0).abs() < 0.04 && (dg / ed - 1.0).abs() < 0.04, "axis {ax} vs {ea}, diagonal {dg} vs {ed}");
        assert!(zs[0].is_nan());
        // The 1-D update equals the GLACIER-1 exact local update when only one axis is known.
        let u = plastic_update_fmm(10.0, 30.0, f64::INFINITY, h0, 128.0);
        assert!(((u - 30.0) * (u - 10.0) - h0 * 128.0).abs() < 1e-6);
        // Two known axes: the residual of the quadratic form vanishes.
        let u2 = plastic_update_fmm(0.0, 100.0, 100.0, h0, 128.0);
        assert!(((u2 * u2) * 2.0 * (u2 - 100.0).powi(2) - (h0 * 128.0).powi(2)).abs() / (h0 * 128.0).powi(2) < 1e-6);
        assert!(u2 < plastic_update_fmm(0.0, 100.0, f64::INFINITY, h0, 128.0));
    }

    #[test]
    fn dinf_splits_a_plane_by_its_aspect_and_conserves_flux() {
        // Lane GLACIER-3 (SYNTHETIC fixture): a plane falling toward +x and +y at aspect atan(0.5) from the x axis: D-infinity
        // sends (1 - r / (pi/4)) of the flux to the +x neighbour and r / (pi/4) to the (+1, +1) diagonal, r = atan(0.5).
        let n = 40usize;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 128_000 };
        let nn = lat.nodes();
        let s: Vec<f64> = (0..nn).map(|i| 1000.0 - 1.0 * (i % n) as f64 - 0.5 * (i / n) as f64).collect();
        let dom = vec![true; nn];
        let lake = vec![false; nn];
        let d = bal_route_dinf(&lat, &s, &dom, &lake);
        let i = 20 * n + 20;
        let r = 0.5f64.atan();
        let sh = 1.0 - r / std::f64::consts::FRAC_PI_4;
        // sh = 0.41 to the cardinal, so the primary (larger share) is the diagonal.
        assert_eq!(d.r1[i] as usize, i + n + 1);
        assert_eq!(d.r2[i] as usize, i + 1);
        assert!((d.f1[i] as f64 - (1.0 - sh)).abs() < 1e-6, "{} vs {}", d.f1[i], 1.0 - sh);
        assert!((d.uy[i] as f64 / d.ux[i] as f64 - 0.5).abs() < 1e-5);
        // Receivers precede donors in pop order, and unit flux from every post reaches the outlets intact.
        let mut rank = vec![0usize; nn];
        for (k, &p) in d.pop.iter().enumerate() { rank[p as usize] = k; }
        let mut q = vec![1.0f64; nn];
        let mut out = 0.0;
        for &p in d.pop.iter().rev() {
            let p = p as usize;
            if d.r1[p] == NONE { out += q[p]; continue; }
            assert!(rank[d.r1[p] as usize] < rank[p]);
            if d.r2[p] != NONE { assert!(rank[d.r2[p] as usize] < rank[p]); let f = d.f1[p] as f64; q[d.r1[p] as usize] += q[p] * f; q[d.r2[p] as usize] += q[p] * (1.0 - f); } else { q[d.r1[p] as usize] += q[p]; }
        }
        assert!((out - nn as f64).abs() < 1e-6 * nn as f64, "{out}");
    }

    #[test]
    fn params_validate() {
        let p = GlacialParamsV1::v1("blake3:x".into());
        p.validate().unwrap();
        let mut q = p.clone();
        q.u_power_permille = 500;
        assert!(q.validate().is_err());
    }
}
