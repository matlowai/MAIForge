//! T1b 0.9.4 (lane TA, 2026-10-14; integrator redirect from design review consult 3, `terrain-strategy-3/STRATEGY.md`):
//! the consistent level solve, the T1a bank ramp and a robust E update. LABEL: pre-signoff, fidelity tier Preview,
//! derived, never canon; construction (no SimTime). Everything is opt-in behind `ls=1` (needs `fx=1`; each item its own
//! switch); without it 0.9.3 is byte-identical. Pre-registered in `unpublished-design-record`.
//!
//! The missing invariant (design review): beds were integrated up each tree from its own outlet while the routing picked the
//! geodesically nearest seed, so divides sat where distances balance, not heights. Here every land node's bed is the
//! minimum over its allowed downhill paths to any outlet of (outlet level + summed bed rise under the 0.9.3 bed law),
//! iterated with the drainage areas (a multi-source Dijkstra fixed point). It replaces TD2 capture routing.
//! T1a levels stay FIXED inputs (ruling R2 pending): T1a valley-axis nodes are seeds at their levels.

use crate::t0::T0;
use crate::t1b::ReliefLawsV1;
use crate::t1b_litho::ExposedLithology;
use crate::t1b_rivernet::{e_at, Network, NodeCtx, RiverNetParamsV1, K_LAND, K_T1A};
use crate::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

pub fn is_zero_f64(v: &f64) -> bool { *v == 0.0 }

pub const LS_GENERATOR_VERSION: &str = "0.9.4";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct LsParamsV1 {
    /// LS1: consistent level solve (flag), pass cap, stop threshold (receivers changed per pass, ppm of land nodes).
    pub solve: u32,
    pub passes_max: u32,
    pub stop_ppm: i64,
    /// REVISION 1 (lane TA PREREG, 2026-10-14 05:1x): hysteresis in mm: a node changes receiver only when the new path is lower
    /// than through its current receiver by more than this.
    pub hysteresis_mm: i64,
    /// LS2: T1a bank ramp (flag).
    pub bank_ramp: u32,
    /// LS3: E update (0 = the 0.9.3 damped Jacobi with CF's sign damping, 1 = Broyden in ln E over the solved cells).
    pub e_update: u32,
    pub solve_rule: String,
    pub bank_rule: String,
    pub e_update_rule: String,
    /// Lane TB options (absent = 0.9.4 as lane TA left it; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tb: Option<TbParamsV1>,
}

/// Lane TB (2026-10-14; `unpublished-design-record`): evidence switches for pending rulings R2 and R4 and
/// the TD2b cap at the eikonal's own inputs. EXPLORATORY / MEASUREMENT, never defaults without a ruling.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct TbParamsV1 {
    /// S1 (R2 exploratory): 0 off (T1a levels fixed inputs), 1 T1a levels solved by the bed law from the flooded mouth
    /// (lowstand base level, as every network channel), 2 the same with the outlet at the present-sea zero contour (final-reach
    /// T1a nodes at bil <= sea keep CF's shelf-graded level as fixed outlets).
    pub t1a_levels: u32,
    /// S2: the TD2b cap with the eikonal's own rock-mass multiplier (channel chain label) and Sc (flag).
    pub steep_eik: u32,
    /// S3 (R4 MEASUREMENT): solved marine T0 cells with land hat fraction below this (permille) use plain per-cell means (0 off).
    pub r4_plain_permille: i64,
    pub keying_rule: String,
    pub t1a_rule: String,
    pub steep_rule: String,
    pub r4_rule: String,
    /// Lane TC options (absent = lane TB as left; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tc: Option<TcParamsV1>,
    /// Lane GUIDE-SPIKES (2026-10-22, opt-in ls_undercut=1, needs ls_steep=2): guide_undercut_v1 (`UNDERCUT_RULE`, carried
    /// in `undercut_rule` when on). Absent / 0 = byte-identical.
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub undercut: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub undercut_rule: String,
}

/// Lane GUIDE-SPIKES: the rule text a recipe carries with ls_undercut=1.
pub const UNDERCUT_RULE: &str = "guide_undercut_v1 (lane GUIDE-SPIKES, 2026-10-22; with steep_eik 2 only): in every eikonal pass a network channel post (code 5) accepts arrivals from every front (soft_channel_cap_v1 accepts only fronts that drain to its own terminal); after pass 2 a network channel post whose height is below its curve bed and whose label drains to another terminal (or is unlabelled water) is no longer a channel: it is published as a free hillslope post (code 0, no channel area or link) and takes every later free-post step (crest rounding, ct remap, fill); posts lowered by their own basin stay channels as before. Basis: a channel head cannot stand above the hillslope of the basin across its divide; without the rule such heads stood as one-post walls up to 5 km high (C7 site 5, 128 m and 32 m evaluation)";

impl TbParamsV1 {
    /// Lane GUIDE-SPIKES: ls_undercut=<v> (the rule text is carried only when on).
    pub fn set_undercut(&mut self, v: u32) {
        self.undercut = v;
        self.undercut_rule = if v == 1 { UNDERCUT_RULE.into() } else { String::new() };
    }
}

/// Lane GUIDE-SPIKES: guide_undercut_v1 on.
pub fn undercut(p: &RiverNetParamsV1) -> bool { tb(p).map(|t| t.undercut == 1).unwrap_or(false) }

/// Lane LS1-GLOBAL (key ls_global=1): the rule text a recipe carries.
pub const GLOBAL_RULE: &str = "ls1_global_order_v1 (lane LS1-GLOBAL, 2026-10-22; opt-in ls_global=1): in every confinement mode other than the cell rule (none, final basin, warp, free) each LS1 pass is ONE multi-source minimum-path (Dijkstra / priority-flood) pass over all land nodes of the chart instead of one pass per T0 cell in downstream-first T0 flow order: seeds and settled nodes, the allowed-move tests (the mode's confinement test, free_bounded_path_v1's path test with the destination of the neighbour as settled in the same pass, outlet tests of T1a neighbours), the edge cost (0.9.3 bed law, rock at the receiver's bed, link through the keyed meander midpoint, T1a floor cap), the hysteresis heap order, the diagonal-crossing refusal and the unconfined fallback for unreached nodes are unchanged. Why: under the per-cell order a node may only drain into a node already settled, so in modes whose allowed moves reach cells later in T0 flow order (any neighbour cell under none; the band under free) a node at a cell edge could not take a lower path through the next cell and stood up to 4 km above its best allowed neighbour path (LS1 optimality residual: C7 site 5 conf=free 3,600 to 4,100 pairs, C7 site 1 conf=none 436 pairs; 0 under the cell rule). Under the cell rule every allowed move goes to the node's own cell or its T0 receiver cell, which the per-cell order settles first, so the per-cell order is already a global minimum-path order there and is kept. Also (every mode): LS1 reads the TC floor-2 T1a floor cap of the network the pass starts from (and the residual that of the published network), as the bed rule reads its own network, instead of the built network's; with the flux recomputation on the T1a areas, and so the floor half-widths, differ between the two and LS1's path beds disagreed with the published bed rule (C7 site 5 conf=free with the global order alone: 2,206 residual pairs up to 1.2 km beside T1a nodes); without the flux both are the T1a discharge areas and nothing changes (conf=cell without conf_flux: byte-identical apart from the recipe)";

/// Lane LS1-GLOBAL: the ls1 global options (params.ls_global; absent = LS1 as before, byte-identical; never seen by the keyed
/// draws: draw_digest strips it, so key-on and key-off share one built network).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct LsGlobalParamsV1 {
    /// 1 = ls1_global_order_v1 (GLOBAL_RULE).
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub order: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub order_rule: String,
    /// 1 = ls1_settle_v1 (SETTLE_RULE).
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub settle: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub settle_rule: String,
    /// ls1_settle_v1 stage caps (passes): stage 2 (consistent areas) and stage 3 (fixed areas).
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub settle_consistent_max: u32,
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub settle_fixed_max: u32,
}

impl LsGlobalParamsV1 {
    pub fn off() -> Self { LsGlobalParamsV1 { order: 0, order_rule: String::new(), settle: 0, settle_rule: String::new(), settle_consistent_max: 0, settle_fixed_max: 0 } }
    pub fn set_order(&mut self, v: u32) { self.order = v; self.order_rule = if v == 1 { GLOBAL_RULE.into() } else { String::new() }; }
    pub fn set_settle(&mut self, v: u32) {
        self.settle = v;
        self.settle_rule = if v == 1 { SETTLE_RULE.into() } else { String::new() };
        (self.settle_consistent_max, self.settle_fixed_max) = if v == 1 { (SETTLE_CONSISTENT_MAX, SETTLE_FIXED_MAX) } else { (0, 0) };
    }
}

/// Lane LS1-GLOBAL key parsing (called from t1b_outlet::parse_key; layer AND windows): ls_global=0|1, ls_settle=0|1, and after
/// ls_settle=1 the stage caps ls_settle_cons=<passes> (0..64), ls_settle_fixed=<passes> (1..64).
pub fn parse_key(s: &str, p: &mut RiverNetParamsV1) -> R<bool> {
    if crate::t1b_esc_local::parse_key(s, p)? { return Ok(true); }
    for (k, which) in [("ls_settle_cons=", 0), ("ls_settle_fixed=", 1)] {
        let Some(d) = s.strip_prefix(k) else { continue };
        let v: u32 = d.parse().map_err(|_| format!("bad value {s}"))?;
        let g = p.ls_global.as_mut().filter(|g| g.settle == 1).ok_or(format!("{k} needs ls_settle=1 given before it"))?;
        if v > 64 || (which == 1 && v == 0) { return Err(format!("{s}: ls_settle_cons 0..64, ls_settle_fixed 1..64")); }
        if which == 0 { g.settle_consistent_max = v; } else { g.settle_fixed_max = v; }
        return Ok(true);
    }
    for (k, which) in [("ls_global=", 0), ("ls_settle=", 1)] {
        let Some(d) = s.strip_prefix(k) else { continue };
        let v: u32 = match d { "0" => 0, "1" => 1, _ => return Err(format!("bad value {s} (0 or 1)")) };
        let g = p.ls_global.get_or_insert_with(LsGlobalParamsV1::off);
        if which == 0 { g.set_order(v); } else { g.set_settle(v); }
        if g.order == 0 && g.settle == 0 { p.ls_global = None; }
        return Ok(true);
    }
    Ok(false)
}

/// Lane LS1-GLOBAL: the options need LS1 (ls=1, ls_solve=1) and the chart-wide solve (not ls_cell).
pub fn validate_global(p: &RiverNetParamsV1) -> R<()> {
    let Some(g) = &p.ls_global else { return Ok(()) };
    if g.order > 1 || g.settle > 1 { return Err("ls_global and ls_settle are 0 or 1".into()); }
    if p.ls.as_ref().map(|l| l.solve != 1).unwrap_or(true) { return Err("ls_global / ls_settle need ls=1 (LS1)".into()); }
    if te_cell(p) { return Err("ls_global / ls_settle need the chart-wide LS1 (not ls_cell)".into()); }
    Ok(())
}

/// Lane LS1-GLOBAL: ls1_global_order_v1 applies (on, and the confinement is not the cell rule).
pub fn global_order(p: &RiverNetParamsV1) -> bool {
    p.ls_global.as_ref().map(|g| g.order == 1).unwrap_or(false)
        && p.confinement.as_ref().map(|c| c.mode != crate::t1b_outlet::CONF_CELL).unwrap_or(false)
}

/// Lane LS1-GLOBAL: ls1_settle_v1 on.
pub fn settle(p: &RiverNetParamsV1) -> bool { p.ls_global.as_ref().map(|g| g.settle == 1).unwrap_or(false) }

pub const SETTLE_CONSISTENT_MAX: u32 = 16;
pub const SETTLE_FIXED_MAX: u32 = 16;
pub const SETTLE_RULE: &str = "ls1_settle_v1 (lane LS1-GLOBAL, 2026-10-22; opt-in ls_settle=1, any confinement mode): LS1 ends only on a pass that changes no land receiver, and publishes the areas that pass routed with, so LS1 run again on its output reproduces it exactly (idempotent: the network the layer re-solves after the E solve is the last evaluation's network, on which F6 was measured). Every pass reads land areas equal to the areas it routes with (seeds, chain labels and the T1a floor cap included). Stage 1: the damped area iteration as before (A_{k+1} = sqrt(A_k A_new)) until a pass changes at most stop_ppm of land receivers or passes_max; stage 2 (at most settle_consistent_max passes): consistent areas, each pass routes with the f32 accumulation of the previous pass's receivers (a consistent fixed point when a pass changes nothing: the published areas are the accumulation of the published receivers); stage 3 (at most settle_fixed_max passes): the areas are held at the f32 accumulation of the network it starts from, until a pass changes no receiver (B1's area convention on the nodes it reroutes); if stage 3 hits its cap the last pass is published as B1 does and the record says settled false. Diagnosis that motivates it (lane LS1-GLOBAL, C7 site 1 conf=none and site 5 conf=free, env VZ_LS1G_DIAG): receivers keep flipping because the bed law depends on drainage area and the areas on the routing (a capture cascade: a channel that gains area gets lower beds and draws its neighbours in the next pass), not ties (flips gain 5 to over 100 m over the old receiver's path, fewer than 0.1 % within 2 m) nor hysteresis cycling (under 1 % return to the receiver of two passes before); the damped stage converges slowly (about 10 % fewer flips per pass) and the old final pass on undamped areas flips about twice the last damped pass, so LS1 run on its own output changed 225 to 2,600 receivers and every bed path downstream of them";

thread_local! { static IDEM_GUARD: std::cell::Cell<bool> = const { std::cell::Cell::new(false) }; }

/// Lane TC (2026-10-14; `unpublished-design-record`): R2 floor rule and shelf-head clip, design review's R4 variant,
/// LS1 sensitivity / receiver freeze, T1a trunk curve. EXPLORATORY / MEASUREMENT, never defaults without a ruling.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct TcParamsV1 {
    /// C1: 0 off, 1 window min to the T1a floor, 2 = 1 + node grading in the bed law (node beds and LS1).
    pub floor: u32,
    /// C2: clip the drowned head prefix of T1a axes (needs t1a_levels != 0).
    pub clip: u32,
    /// C3 (R4 MEASUREMENT): cells the existing rule solves with land hat fraction below this (permille) take their solved
    /// 8-neighbours' mean E instead (0 off).
    pub r4n_permille: i64,
    /// C4 (MEASUREMENT): E perturbation (ppm) and where (1 the published table, 2 the initial table); freeze LS1 from this E
    /// evaluation on (0 off).
    pub pert_ppm: i64,
    pub pert_mode: u32,
    pub freeze_from: u32,
    /// C5: T1a links through keyed meander midpoints.
    pub curve: u32,
    pub floor_rule: String,
    pub clip_rule: String,
    pub r4n_rule: String,
    pub pert_rule: String,
    pub freeze_rule: String,
    pub curve_rule: String,
    /// Lane TE options (absent = lane TC as left; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub te: Option<TeParamsV1>,
    /// Lane TF options (absent = lane TE as left; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf: Option<crate::t1b_tf::TfParamsV1>,
    /// Lane TG options (absent = lane TF as left; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tg: Option<crate::t1b_tg::TgParamsV1>,
    /// Lane TH options (absent = lane TG as left; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub th: Option<crate::t1b_th::ThParamsV1>,
    /// Lane VE options (absent = lane TH as left; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ve: Option<crate::t1b_ve::VeParamsV1>,
    /// Lane STEP options (absent = lane VE as left; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ct: Option<crate::t1b_ct::CtParamsV1>,
}

/// Lane TE (2026-10-14; `unpublished-design-record`): the cell-local E solve of design review consult 3 Q4
/// (`solve_cell`). EXPLORATORY, opt-in, never a default without a ruling.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct TeParamsV1 {
    /// 1 = solve E (and F) per T0 cell (t1b_cell::solve_e_cells) instead of the chart-wide Jacobi.
    pub cell_solve: u32,
    /// Grids of the per-cell unknowns: ln E step (ppm) and F step (mm).
    pub e_grid_ppm: i64,
    pub f_grid_mm: i64,
    /// Sub-lattice margin around a cell's hat support (mm), Gauss-Seidel sweep cap, evaluations per 1-D root cap,
    /// E/F alternations per visit cap.
    pub sub_margin_mm: i64,
    pub sweeps_max: u32,
    pub root_evals_max: u32,
    pub alternations_max: u32,
    pub cell_rule: String,
}

impl TeParamsV1 {
    pub fn v1() -> Self {
        TeParamsV1 {
            cell_solve: 1, e_grid_ppm: 5_000, f_grid_mm: 1_000, sub_margin_mm: 12_000_000, sweeps_max: 4, root_evals_max: 16, alternations_max: 4,
            cell_rule: "solve_cell_v1 (design review consult 3 Q4, EXPLORATORY): the T0 cell is the solve unit. Cells with land nodes in downstream-first T0 flow order (chain length, cell id), Gauss-Seidel sweeps until a sweep changes no grid value and no receiver (sweeps_max). Visiting a cell: LS1 routes only the cell's land nodes (level_solve scope; confinement as LS1: own cell or its T0 receiver cell), every other node settled at its current bed and receiver, COLD from the built network's receivers for the cell (no warm start), stop limit over the cell's land nodes; published areas: the cell's last-pass areas (B1), elsewhere the f32 accumulation of the current receivers. Solved cells (as solve_e, with r4n): alternate 1-D roots on fixed grids until an alternation changes neither (alternations_max): E on ln E = k e_grid (root of M - T, hat means over the cell's land posts as solve_e), F on k f_grid in [0, flood_max] (root of land fraction - bilinear land fraction, plain over the cell's posts, coastal cells only); each root: from the current grid value, first step the damped-Jacobi (E) or flood-rule (F) step, extend until the sign changes or a bound, then Illinois regula falsi to adjacent grid points, the one with the smaller |residual| (ties: lower), root_evals_max per root (best evaluated if hit). One evaluation = LS1 for the cell on the trial tables, node beds of the whole network, curves, the coast context, evaluate() on the cell's sub-lattice (hat support + sub_margin at the solve spacing). Unsolved cells with land nodes (and r4n cells): routed once per visit; r4n E follows its solved neighbours after every commit. Published: the cell-solved network (no LS1 re-run), the tables, and rows of one full-chart evaluation (convergence: every solved |M - T| <= tolerance and F4 as solve_e). pert_mode 2 perturbs the start table; pert_mode 1 scales the solved table and re-routes every cell once".into(),
        }
    }
}

impl TcParamsV1 {
    pub fn off() -> Self {
        TcParamsV1 {
            floor: 0, clip: 0, r4n_permille: 0, pert_ppm: 0, pert_mode: 0, freeze_from: 0, curve: 0,
            floor_rule: "t1a_floor_grade_v1 (R2 evidence): floor 1 = in the window, before the general floor rule, a network channel post (code 5) within W/2 of a T1a channel post (code 1) takes min(own bed, the T1a floor value) with W and the floor value of nearest_channel_floor_v1 taken over T1a channel posts only; it stays a channel post and the general floor reads the lowered beds. floor 2 = floor 1 + node grading: in the bed law of node_beds_core and in LS1's edge cost, a land node (not a conduit) within the T1a floor half width h_t = min(floor_width_coeff sqrt(A_t) / 2, floor_width_max / 2) of a T1a channel link t -> r (lattice segment, distance d, bed b linear along it, receiver bed min(bed r, bed t)) has bed = max(bed of its receiver, min(law bed, b + max(floor_slope, link slope) d)) (minimum over links): a tributary runs on the trunk floor and its bed-law rise starts at the floor edge".into(),
            clip_rule: "t1a_shelf_head_clip_v1: with T1a levels solved (t1a_levels != 0), a T1a valley-axis node at T0 bilinear <= sea (the shelf band; seaward of sea - shelf it is construction water already) whose T1a donors are all clipped or construction water (or none) is clipped, upstream first: the drowned head prefix of an axis. Clipped nodes are ordinary construction land nodes of the network (keyed draws of a land node at their index; no other draw changes); their T1a links are not drawn and a clipped donor is not a Catmull-Rom tangent. T1a discharge areas of the remaining nodes unchanged".into(),
            r4n_rule: "r4_neighbour_mean_v1 (design review consult 3 Q3 item 3, R4 MEASUREMENT ONLY): a cell the existing rule solves (land hat >= solve_min_land_hat) whose land hat fraction is below r4n_permille is not solved; its E is the mean E of its solved T0 8-neighbours, set from the initial table and after every update (no solved neighbour: unchanged)".into(),
            pert_rule: "e_perturbation_v1 (MEASUREMENT ONLY): pert_mode 1 = after the solve every table entry x (1 + pert_ppm / 1e6) and the published network and beds recomputed from the last evaluation's network (LS1's response to the table alone); pert_mode 2 = the initial table x (1 + pert_ppm / 1e6) (the whole solve sees it)".into(),
            freeze_rule: "ls_receiver_freeze_v1 (MEASUREMENT): from E evaluation freeze_from on (0 off), LS1 is not re-run: the network is the previous evaluation's (receivers and its published-precision areas) and beds are the 0.9.2 bed rule on it; the published network is that frozen network; LS1's optimality residual on it is reported".into(),
            te: None,
            tf: None,
            tg: None,
            th: None,
            ve: None,
            ct: None,
            curve_rule: "t1a_meander_curve_v1: a T1a link n -> r is Catmull-Rom through n, one keyed meander midpoint m_n and r (m_n = link midpoint moved sideways by (2u - 1) x meander_permille x link length, u = the t1b.rivernet.meander uniform at n's index, the draw a land node there gets; T1a nodes on construction water: no offset); outer control points are the neighbouring links' midpoints (largest-discharge T1a donor, the receiver's link; else reflected), so the curve is C1 at nodes and passes through every T1a node (F1 strict unchanged); 16 samples per half, bed linear in the parameter".into(),
        }
    }
}

pub fn tc(p: &RiverNetParamsV1) -> Option<&TcParamsV1> { tb(p).and_then(|t| t.tc.as_ref()) }
pub fn tc_floor(p: &RiverNetParamsV1) -> u32 { tc(p).map(|t| t.floor).unwrap_or(0) }
pub fn tc_clip(p: &RiverNetParamsV1) -> bool { tc(p).map(|t| t.clip == 1).unwrap_or(false) }
pub fn tc_r4n(p: &RiverNetParamsV1) -> i64 { tc(p).map(|t| t.r4n_permille).unwrap_or(0) }
pub fn tc_pert(p: &RiverNetParamsV1) -> (i64, u32) { tc(p).map(|t| (t.pert_ppm, t.pert_mode)).unwrap_or((0, 0)) }
pub fn tc_freeze(p: &RiverNetParamsV1) -> u32 { tc(p).map(|t| t.freeze_from).unwrap_or(0) }
pub fn tc_curve(p: &RiverNetParamsV1) -> bool { tc(p).map(|t| t.curve == 1).unwrap_or(false) }
pub fn te(p: &RiverNetParamsV1) -> Option<&TeParamsV1> { tc(p).and_then(|t| t.te.as_ref()) }
pub fn te_cell(p: &RiverNetParamsV1) -> bool { te(p).map(|t| t.cell_solve == 1).unwrap_or(false) }

/// Lane TC C1 floor 2: per land node the T1a floor cap (INFINITY where none), from the T1a beds in `bed` (None when off).
pub fn t1a_floor_cap(net: &Network, bed: &[f64], p: &RiverNetParamsV1) -> Option<Vec<f64>> {
    if tc_floor(p) != 2 { return None; }
    let cw = p.floor_width_coeff_ppm as f64 / 1e6;
    let hmax = p.floor_width_max_mm as f64 / 2000.0;
    let fs = p.floor_slope_ppm as f64 / 1e6;
    let lattice = |n: usize| { let (gx, gy) = node_g(n); [(gx * SPACING_MM) as f64 / 1000.0, (gy * SPACING_MM) as f64 / 1000.0] };
    let ok_r = |r: u32| r != NONE && r != OFF_CHART;
    let mut cap = vec![f64::INFINITY; NODES];
    let reach = (hmax / (SPACING_MM as f64 / 1000.0)).ceil() as i64 + 1;
    for t in 0..NODES {
        if net.kind[t] != crate::t1b_rivernet::K_T1A || !net.t1a_ch[t] || !ok_r(net.t1a_recv[t]) || !bed[t].is_finite() { continue; }
        let r = net.t1a_recv[t] as usize;
        let (a, b) = (lattice(t), lattice(r));
        let bt = bed[t];
        let br = if bed[r].is_finite() { bed[r].min(bt) } else { bt };
        let len = (b[0] - a[0]).hypot(b[1] - a[1]).max(1e-9);
        let slope = (bt - br) / len;
        let h = (0.5 * cw * net.area_m2[t].max(0.0).sqrt()).min(hmax);
        let (gx, gy) = node_g(t);
        for oy in -reach..=reach { for ox in -reach..=reach {
            let Some(n) = node_index(gx + ox, gy + oy) else { continue };
            if net.kind[n] != K_LAND || net.conduit[n] { continue; }
            let q = net.pos[n];
            let u = (((q[0] - a[0]) * (b[0] - a[0]) + (q[1] - a[1]) * (b[1] - a[1])) / (len * len)).clamp(0.0, 1.0);
            let d = (q[0] - a[0] - u * (b[0] - a[0])).hypot(q[1] - a[1] - u * (b[1] - a[1]));
            if d > h { continue; }
            let v = bt + u * (br - bt) + fs.max(slope) * d;
            if v < cap[n] { cap[n] = v; }
        } }
    }
    Some(cap)
}

impl TbParamsV1 {
    pub fn off() -> Self {
        TbParamsV1 {
            t1a_levels: 0, steep_eik: 0, r4_plain_permille: 0,
            keying_rule: "tb_paired_keying_v1: the lane TB options do not re-key draws: the network and volcano draws use the digest of this recipe with params.ls.tb removed (so a TB option on and off at the same realization share every keyed draw); the layer and window recipes still carry the full recipe digest".into(),
            t1a_rule: "t1a_level_solve_v1 (R2 EXPLORATORY, ruling pending): T1a valley-axis nodes keep planform, T1a receiver, discharge area and curve geometry, but not their levels: downstream first along the T1a receivers, bed = bed(outlet) + L x s with the 0.9.3 bed law of land nodes (s = min(k_s A^-theta, cap x Sc, gamma x S_h), A = T1a discharge area, rock at the node at the outlet's bed, L = lattice chord); outlet = the T1a receiver when it is a T1a node, else the first node of the node's CO conduit, else a marine (base level sea - F) or lake (stage) receiver, else (no receiver, off chart) the fixed T1a level (counted); conduit beds rise from their receiver by the bed law without the terminal-level cap; CF's shelf grade caps final-reach nodes (min); mode 2: final-reach nodes at bil <= sea keep CF's shelf-graded level as fixed outlets (outlet at the present-sea zero contour). LS1 takes these beds as settled seeds; T1a curves take the solved levels at both link ends. The E solve's L (T1a trunk-level hat mean, Jacobi step size only) is unchanged".into(),
            steep_rule: "REVISION 1 (lane TB PREREG, 2026-10-14 05:52): steep_eik 2 = soft_channel_cap_v1 (window and E-solve evaluation only; node beds and LS1 unchanged): in both eikonal passes network channel posts (code 5) are tentative sources like TD3 floors, published z = min(curve bed, the eikonal arrival from a front whose label drains to the same terminal; other basins and unlabelled water are not accepted), the speed at such a post the pass's hillslope speed with L_H = the mean over hillslope posts within 64 m (TD3 B4; pass 1 L_ref); T1a channels stay hard (F1). steep_eik 1 (as registered) = td2b_eikonal_inputs_v1: TD2b's valley-side gradient S_h(n) uses the inputs the window eikonal applies to the banks of n's channel link: the rock-mass multiplier of n's channel chain label (the chain start node as curves() labels it, labels from a first bed pass on the same network; non-channel nodes: their own), Sc = min(Sc at the receiver's bed, Sc at the receiver's bed + the capped rise); L_n unchanged; the same in node_beds_core and in LS1's edge cost".into(),
            r4_rule: "r4_plain_reference_v1 (R4 MEASUREMENT ONLY, ruling pending): a solved marine T0 cell with land hat fraction below r4_plain_permille uses the plain per-cell means over its own land posts (edifices excluded, as lane CO's plain rows) for M, T and L in the residual, the convergence test and the Jacobi step".into(),
            tc: None,
            undercut: 0, undercut_rule: String::new(),
        }
    }
}

pub fn tb(p: &RiverNetParamsV1) -> Option<&TbParamsV1> { p.ls.as_ref().and_then(|l| l.tb.as_ref()) }
pub fn t1a_mode(p: &RiverNetParamsV1) -> u32 { tb(p).map(|t| t.t1a_levels).unwrap_or(0) }
pub fn steep_eik(p: &RiverNetParamsV1) -> bool { tb(p).map(|t| t.steep_eik == 1).unwrap_or(false) }
pub fn soft_channels(p: &RiverNetParamsV1) -> bool { tb(p).map(|t| t.steep_eik == 2).unwrap_or(false) }
pub fn r4_permille(p: &RiverNetParamsV1) -> i64 { tb(p).map(|t| t.r4_plain_permille).unwrap_or(0) }

/// S2: the channel chain label of every node as curves() assigns it (chain start node; conduits: their terminal T1a node;
/// NONE for non-channel nodes).
pub fn chain_labels(net: &Network, chan: &[bool]) -> Vec<u32> {
    let mut donors = vec![0u32; NODES];
    for n in 0..NODES { if net.kind[n] == K_LAND && chan[n] && net.recv[n] != NONE { donors[net.recv[n] as usize] += 1; } }
    let mut lab = vec![NONE; NODES];
    for s in 0..NODES {
        if net.kind[s] != K_LAND || !chan[s] || net.recv[s] == NONE || donors[s] == 1 || net.conduit[s] { continue; }
        lab[s] = s as u32;
        let mut cur = s;
        for _ in 0..NODES {
            let r = net.recv[cur];
            if r == NONE { break; }
            let r = r as usize;
            if net.kind[r] != K_LAND || donors[r] != 1 || net.recv[r] == NONE { break; }
            lab[r] = s as u32;
            cur = r;
        }
    }
    for path in &net.conduit_paths { for &j in &path[1..] { let j = j as usize; if net.conduit[j] { lab[j] = path[0]; } } }
    lab
}

impl LsParamsV1 {
    pub fn v1() -> Self {
        LsParamsV1 {
            solve: 1, passes_max: 8, stop_ppm: 1000, hysteresis_mm: 2000, bank_ramp: 1, e_update: 0,
            solve_rule: "B1+B2 (fixes): the hysteresis orders the heap only (stored beds are true path beds); after the stop (or passes_max) one final pass runs on the f32-rounded areas of the network it starts from and those areas are published with its receivers (windows read them), so the published beds are that pass's minimum-path beds (on nodes the final pass rerouted, the published area is the pre-pass area). REVISION 1: relaxing i from j adds hysteresis_mm (2 m) when j is not i's receiver in the current network (the input network, or the previous pass / E evaluation), so near-equal paths do not flip (run ls1: 250-1200 receivers flipped per pass and never stopped). consistent_level_solve_v1 (replaces TD2 capture routing): seeds = marine nodes (sea - F), lakes (stage), T1a valley-axis nodes (fixed T1a levels; CF shelf grade as 0.9.2), conduits keep their receivers and the 0.9.2 conduit bed; per T0 cell, downstream-first in T0 flow order (depth of the T0 chain, then cell id), solve_cell(c) is a multi-source Dijkstra over c's land nodes from every settled node: i may drain to an 8-neighbour j under the network confinement (j water or conduit, or j in i's T0 cell or its T0 receiver cell) without crossing a diagonal link, bed_i = bed_j + s_i L_ij with the 0.9.3 bed rule (s = min(k_s A^-theta, cap Sc, gamma S_h), rock at i at bed_j, L through i's keyed meander midpoint), heap key (bed, node, from); nodes unreached under confinement get one unconfined pass; the receiver is the predecessor. Areas iterate: A_0 = the input (or previous E evaluation's) network, A_{k+1} = sqrt(A_k A_new), until a pass changes < stop_ppm of land receivers or passes_max; published receivers = the last pass; beds = the 0.9.2 bed rule on that network (windows rebuild it from the receivers). Basis (recalled, unverified): steady-state divides sit where the channel heights from both sides meet (Willett et al. 2014 Science 343:1248765)".into(),
            bank_rule: "bank_ramp_v1 (fixes t1a_bank_v1's one-post walls): after the bank raise, every free post satisfies z_i >= z_j - g_i d_ij over 8-neighbours (d = 32 m or 45.25 m), g_i the post's pass-2 hillslope gradient at its pre-ramp height (Sc x rockmass x f(beta(L_H)), floored at hill_slope_min), by max-first propagation from the raised posts; a T1a channel above its surroundings gets valley sides at the hillslope gradient; never lowers a post".into(),
            e_update_rule: "OUTCOME of the registered decision rule (2026-10-14 05:1x): e_update 0 (the 0.9.3 damped Jacobi with CF's sign damping) is the v1 value; Broyden (1) is kept as an option: in run ls1 its Jacobian went near-singular (cell 41995's own column is ~0: dM/dlnE -0.3 m) and it threw every cell by ln 2 twice. broyden_ln_e_v1: good Broyden in x = ln E over the solved cells, residual r = M - T (hat means as 0.9.0), J0 = diag(max(M - L, 1 m) / 0.6) (the Jacobian the 0.9.0 damped step implies), rank-one update J += (dr - J dx) dx^T / dx.dx after each evaluation, step = -J^-1 r (Gaussian elimination, partial pivoting; singular -> the 0.9.0 step), clamped to +-ln 2 per cell, E projected onto [e_min, e_max]; the flooding-offset (F) update unchanged".into(),
            tb: None,
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.passes_max == 0 { return Err("ls: passes_max must be >= 1".into()); }
        if self.stop_ppm < 0 { return Err("ls: stop_ppm".into()); }
        if self.e_update > 1 { return Err("ls: e_update 0 or 1".into()); }
        if self.hysteresis_mm < 0 { return Err("ls: hysteresis".into()); }
        if let Some(t) = &self.tb {
            if t.t1a_levels > 2 || t.steep_eik > 2 || !(0..=1000).contains(&t.r4_plain_permille) { return Err("ls tb: t1a 0-2, steep 0-1, r4 0-1000 permille".into()); }
            if self.solve != 1 { return Err("ls tb options need ls_solve=1".into()); }
            if t.undercut > 1 || (t.undercut == 1 && t.steep_eik != 2) { return Err("ls_undercut is 0 or 1 and needs ls_steep=2".into()); }
            if let Some(c) = &t.tc {
                if c.floor > 2 || c.clip > 1 || c.curve > 1 || c.pert_mode > 2 || !(0..=1000).contains(&c.r4n_permille) { return Err("ls tc: floor 0-2, clip 0-1, curve 0-1, pert_mode 0-2, r4n 0-1000 permille".into()); }
                if c.clip == 1 && t.t1a_levels == 0 { return Err("ls_clip needs ls_t1a=1 or 2".into()); }
                if (c.pert_ppm != 0) != (c.pert_mode != 0) { return Err("ls_pert_ppm and ls_pert_mode go together".into()); }
                if let Some(e) = &c.te {
                    if e.cell_solve > 1 || e.e_grid_ppm <= 0 || e.f_grid_mm <= 0 || e.sub_margin_mm < 0 || e.sweeps_max == 0 || e.root_evals_max < 2 || e.alternations_max == 0 { return Err("ls te: cell 0-1, grids > 0, margin >= 0, sweeps >= 1, root evals >= 2, alternations >= 1".into()); }
                    if e.cell_solve == 1 && (c.freeze_from != 0 || self.e_update != 0) { return Err("ls_cell excludes ls_freeze and ls_e=1".into()); }
                }
            }
        }
        Ok(())
    }
}

/// LS1 statistics (per call).
#[derive(Clone, Debug, Default, Serialize)]
pub struct LsStats {
    pub passes: u32,
    pub changed_per_pass: Vec<u32>,
    pub unreached_confined_per_pass: Vec<u32>,
    pub receivers_changed_vs_input: u32,
    pub land_nodes: u32,
    /// Optimality residual on the published beds: allowed land pairs (i, j) where draining i to j would give i a bed
    /// lower than its own by more than 1 cm (count, max m).
    pub residual_pairs: u32,
    /// Complete directed witnesses for incremental auditing. Never changes legacy header bytes.
    #[serde(skip)]
    pub residual_witnesses: Vec<(u32, u32, f64)>,
    pub residual_max_m: f64,
    /// Pairs excluded from the residual: the link would cross an existing diagonal link; the gain is within the hysteresis.
    pub residual_crossing: u32,
    pub residual_in_hysteresis: u32,
    /// The 12 largest residuals: (gx_i, gy_i, gx_j, gy_j, excess m, kind of j).
    pub residual_top: Vec<(i64, i64, i64, i64, f64, u8)>,
    /// Diagnostics (not used by the solve): up to 16 receiver changes of the last pass: (gx, gy, old dx, old dy, new dx, new dy,
    /// Dijkstra bed, damped area km2, area of the new network km2), and for the 6 largest residuals (bed_i, bed_j, rise, area_i km2,
    /// area_j km2, receiver dx dy of i, bed of i's receiver).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub last_changes: Vec<(i64, i64, i64, i64, i64, i64, f64, f64, f64)>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub residual_detail: Vec<(f64, f64, f64, f64, f64, i64, i64, f64)>,
    /// Lane LS1-GLOBAL ls1_settle_v1 (only with ls_settle=1).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub settle: Option<SettleStats>,
}

/// Lane LS1-GLOBAL ls1_settle_v1: the stage LS1 ended in (1 damped, 2 consistent, 3 fixed areas), receivers changed per pass
/// in each stage, and whether a pass changed no receiver (settled; LS1 run again on the output reproduces it).
#[derive(Clone, Debug, Default, Serialize)]
pub struct SettleStats { pub stage: u32, pub changed_per_stage: Vec<Vec<u32>>, pub converged: bool }

#[derive(Clone, Copy)]
struct Key(f64, u32, u32);
impl PartialEq for Key { fn eq(&self, o: &Self) -> bool { self.cmp(o) == Ordering::Equal } }
impl Eq for Key {}
impl PartialOrd for Key { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
impl Ord for Key { fn cmp(&self, o: &Self) -> Ordering { o.0.total_cmp(&self.0).then(o.1.cmp(&self.1)).then(o.2.cmp(&self.2)) } }

/// LS1. `beds_core` is the 0.9.2 bed rule (seed levels, conduit beds and the final beds).
#[allow(clippy::too_many_arguments)]
pub fn level_solve(net0: &Network, start: Option<&Network>, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>,
    p: &RiverNetParamsV1, ls: &LsParamsV1, beds_core: &dyn Fn(&Network) -> Vec<f64>, labels_of: &dyn Fn(&Network) -> Option<Vec<u32>>) -> (Network, LsStats) {
    level_solve_f(net0, start, ctx, t0, laws, table, p, ls, beds_core, labels_of, false)
}

/// As `level_solve`; lane TC C4 `frozen`: no pass runs (the network is `start`), only the residual is computed.
#[allow(clippy::too_many_arguments)]
pub fn level_solve_f(net0: &Network, start: Option<&Network>, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>,
    p: &RiverNetParamsV1, ls: &LsParamsV1, beds_core: &dyn Fn(&Network) -> Vec<f64>, labels_of: &dyn Fn(&Network) -> Option<Vec<u32>>, frozen: bool) -> (Network, LsStats) {
    level_solve_s(net0, start, ctx, t0, laws, table, p, ls, beds_core, labels_of, frozen, None)
}

/// As `level_solve_f`; lane TE `scope` (solve_cell): only the land nodes of that T0 cell are routed; every other node is
/// settled at its bed and receiver of the pass's starting network (its own cell's solve owns it). Stop limit, statistics
/// and the residual are over the cell's land nodes. None = every cell (LS1 as lane TA left it).
#[allow(clippy::too_many_arguments)]
pub fn level_solve_s(net0: &Network, start: Option<&Network>, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>,
    p: &RiverNetParamsV1, ls: &LsParamsV1, beds_core: &dyn Fn(&Network) -> Vec<f64>, labels_of: &dyn Fn(&Network) -> Option<Vec<u32>>, frozen: bool, scope: Option<u32>) -> (Network, LsStats) {
    let cells = scope.map(|c| std::collections::BTreeSet::from([c]));
    level_solve_cells(net0, start, ctx, t0, laws, table, p, ls, beds_core, labels_of, frozen, cells.as_ref(), false)
}

/// Multi-cell incremental solve. Exterior routing areas are boundary data, as in LS1 fixed-area stage 3.
/// The legacy single-cell wrapper keeps its original area updates and ordering.
#[allow(clippy::too_many_arguments)]
pub fn level_solve_cells(net0: &Network, start: Option<&Network>, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>,
    p: &RiverNetParamsV1, ls: &LsParamsV1, beds_core: &dyn Fn(&Network) -> Vec<f64>, labels_of: &dyn Fn(&Network) -> Option<Vec<u32>>, frozen: bool, scope: Option<&std::collections::BTreeSet<u32>>, incremental: bool) -> (Network, LsStats) {
    level_solve_mask(net0, start, ctx, t0, laws, table, p, ls, beds_core, labels_of, frozen, scope, incremental, None)
}

/// As the cell API, with an optional native-node domain. Upstream closure is the caller's responsibility.
#[allow(clippy::too_many_arguments)]
pub fn level_solve_mask(net0: &Network, start: Option<&Network>, ctx: &NodeCtx, t0: &T0, laws: &ReliefLawsV1, table: &BTreeMap<u32, f64>,
    p: &RiverNetParamsV1, ls: &LsParamsV1, beds_core: &dyn Fn(&Network) -> Vec<f64>, labels_of: &dyn Fn(&Network) -> Option<Vec<u32>>, frozen: bool,
    scope: Option<&std::collections::BTreeSet<u32>>, incremental: bool, node_scope: Option<&[bool]>) -> (Network, LsStats) {
    let own = crate::t1b_fix::own_area(net0);
    let mut net = start.unwrap_or(net0).clone();
    let boundary_area = if incremental { net.area_m2.clone() } else { Vec::new() };
    let theta = laws.stream_power_m / laws.stream_power_n;
    let cap = p.bed_slope_cap_permille_of_sc as f64 / 1000.0;
    let steep = p.fx.as_ref().filter(|f| f.steep_cap == 1).map(|f| f.steep_gamma_permille as f64 / 1000.0);
    let tau = |n: usize| p.coast.as_ref().map(|c| c.tau(ctx.bil[n])).unwrap_or(1.0);
    let e_node: Vec<f64> = (0..NODES).map(|n| if net0.kind[n] == K_LAND { let (la, lo) = ctx.latlon[n]; e_at(t0, table, la, lo) * tau(n) } else { 0.0 }).collect();
    let cell = |n: usize| net0.t0_cell[n];
    let is_land = |n: usize| net0.kind[n] == K_LAND && !net0.conduit[n];
    let in_scope = |n: usize| scope.map(|c| c.contains(&cell(n))).unwrap_or(true) && node_scope.map(|m| m[n]).unwrap_or(true);
    let land_nodes = (0..NODES).filter(|&n| is_land(n) && in_scope(n)).count() as u32;
    // Lane LS1-GLOBAL: one chart-wide minimum-path pass (non-cell modes, chart-wide LS1 only).
    let global = global_order(p) && (scope.is_none() || incremental);
    let cap_cur = p.ls_global.as_ref().map(|g| g.order == 1).unwrap_or(false);
    // j may receive i: water or conduit, or the confinement rule (t1b_outlet.rs; default: j in i's T0 cell or its T0 receiver cell).
    let conf = crate::t1b_outlet::Confinement::new(p, t0, &net0.t0_cell, &net0.conf_cell, &net0.t1a_recv);
    // Lane DRAIN-FIX conf_core=path: each node's current destination (first non-land node down its receivers; land nodes
    // settled in a pass take their receiver's), read by the core test (t1b_outlet::CORE_PATH_RULE).
    let path = conf.path_mode();
    let dest: Vec<std::cell::Cell<u32>> = vec![std::cell::Cell::new(NONE); if path { NODES } else { 0 }];
    let set_dest = |n: &Network| { if path { for (k, v) in crate::t1b_outlet::destinations(n).into_iter().enumerate() { dest[k].set(v); } } };
    let allowed = |j: usize, i: usize| -> bool {
        // Lane T0-IMPRINT-PROD basin modes only: a T1a node receives only water of its own outlet.
        if net0.kind[j] == K_T1A && conf.is_basin() { return conf.allows_nodes(t0, i, j, true); }
        if net0.kind[j] != K_LAND || net0.conduit[j] { return true; }
        if path { return conf.path_ok(t0, &net0.kind, i, dest[j].get()); }
        conf.allows_nodes(t0, i, j, false)
    };
    let crosses = |i: usize, j: usize, recv: &[u32]| -> bool {
        let (xi, yi) = node_g(i);
        let (xj, yj) = node_g(j);
        let (dx, dy) = (xj - xi, yj - yi);
        if dx == 0 || dy == 0 { return false; }
        let (Some(a), Some(b)) = (node_index(xi + dx, yi), node_index(xi, yi + dy)) else { return false };
        recv[a] == b as u32 || recv[b] == a as u32
    };
    // Link length through i's keyed meander midpoint (as rebuild() places it).
    let link_len = |i: usize, j: usize| -> f64 {
        let (a, b) = (net0.pos[i], net0.pos[j]);
        let len = (b[0] - a[0]).hypot(b[1] - a[1]).max(1e-9);
        let off = (2.0 * net0.meander_u[i] - 1.0) * p.meander_permille as f64 / 1000.0 * len;
        let m = [0.5 * (a[0] + b[0]) - (b[1] - a[1]) / len * off, 0.5 * (a[1] + b[1]) + (b[0] - a[0]) / len * off];
        (a[0] - m[0]).hypot(a[1] - m[1]) + (m[0] - b[0]).hypot(m[1] - b[1])
    };
    let fxp = p.fx.as_ref().filter(|f| f.steep_cap == 1);
    let rise = |i: usize, zj: f64, area: f64, len: f64, lab: Option<&[u32]>| -> f64 {
        let (sc, km, _) = ctx.rock.at(i, zj);
        let ks = e_node[i] / (laws.stream_power_k_per_my * km);
        let mut s = (ks * area.powf(-theta)).min(cap * sc);
        // Lane TB S2: the TD2b cap with the eikonal's inputs (chain labels of the pass's starting network).
        if let (Some(fx), Some(_)) = (fxp, lab) { return crate::t1b_rivernet::steep_capped(s, i, zj, len, e_node[i], sc, net0, ctx, laws, p, fx, lab) * len; }
        if let Some(g) = steep { s = s.min(g * crate::t1b_fix::side_gradient(e_node[i], sc, net0.rockmass[i], crate::t1b_fix::link_hill_length(len), laws, p)); }
        s * len
    };
    // Lane TC C1 floor 2: the bed through j, capped at i's T1a floor value (never below bed_j); plain sum when off.
    let up = |i: usize, zj: f64, r: f64, cap: Option<&[f64]>| -> f64 { match cap { Some(c) if c[i].is_finite() => zj.max((zj + r).min(c[i])), _ => zj + r } };
    // Cells in downstream-first T0 flow order.
    let mut cells: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for n in 0..NODES { if is_land(n) && in_scope(n) { cells.entry(cell(n)).or_default().push(n as u32); } }
    let mut corder: Vec<(usize, u32)> = cells.keys().map(|&c| (t0.chain_from(c).map(|v| v.len()).unwrap_or(1), c)).collect();
    corder.sort();
    let mut st = LsStats { land_nodes, ..Default::default() };
    let mut area: Vec<f64> = net.area_m2.clone();
    let hyst = ls.hysteresis_mm as f64 / 1000.0;
    let limit = (ls.stop_ppm as f64 * 1e-6 * land_nodes as f64).ceil() as u32;
    // B1 (2026-10-14 05:2x): the last pass runs on the areas of the network it starts from (f32-rounded, the published
    // precision) and those areas are published with its receivers, so the published beds are that pass's minimum-path beds
    // (without the hysteresis); before, beds were recomputed with areas no pass had used and stood up to 400 m above a
    // neighbour's path (run ls2 / ta).
    let mut final_pass = false;
    let dg = std::env::var("VZ_LS1G_DIAG").is_ok();
    let (mut dg_prev2, mut dg_area_prev, mut dg_base_prev): (Vec<u32>, Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new(), Vec::new());
    let mut dg_flips = vec![0u8; if dg { NODES } else { 0 }];
    // Lane LS1-GLOBAL ls1_settle_v1: stage caps and record (None when off: the pass loop is LS1 as before).
    let (s_cons, s_fixed) = p.ls_global.as_ref().map(|g| (g.settle_consistent_max, g.settle_fixed_max)).unwrap_or((0, 0));
    let mut sst = (settle(p) && !frozen).then(|| SettleStats { stage: if incremental { 3 } else { 1 }, changed_per_stage: vec![Vec::new(); 3], converged: false });
    if sst.is_some() { for n in 0..NODES { if is_land(n) { net.area_m2[n] = area[n]; } } }
    let pass_cap = if incremental && sst.is_some() { s_fixed.max(1) } else if sst.is_some() { ls.passes_max + s_cons + s_fixed } else { ls.passes_max };
    for pass in 0..=pass_cap {
        if frozen { break; }
        st.passes += 1;
        let base = beds_core(&net);
        // Lane LS1-GLOBAL (with ls_global=1): the T1a floor cap of the network this pass starts from (flux mode recomputes T1a
        // areas, and the bed rule reads the network's own), not the built network's.
        let capv = t1a_floor_cap(if cap_cur { &net } else { net0 }, &base, p);
        let cap = capv.as_deref();
        let labv = labels_of(&net);
        let lab = labv.as_deref();
        let prev: Vec<u32> = net.recv.clone();
        set_dest(&net);
        let pen = |i: usize, j: usize| if prev[i] == j as u32 { 0.0 } else { hyst };
        let mut bed = vec![f64::NAN; NODES];
        let mut done = vec![false; NODES];
        let mut recv = vec![NONE; NODES];
        for n in 0..NODES {
            if !is_land(n) { done[n] = true; bed[n] = base[n]; if net0.conduit[n] { recv[n] = net.recv[n]; } }
            else if !in_scope(n) { done[n] = true; bed[n] = base[n]; recv[n] = net.recv[n]; }
        }
        let relax = |i: usize, heap: &mut BinaryHeap<Key>, bed: &[f64], done: &[bool], recv: &[u32], confined: bool| {
            let (gx, gy) = node_g(i);
            for (dx, dy) in D8 {
                let Some(j) = node_index(gx + dx, gy + dy) else { continue };
                if !done[j] || !bed[j].is_finite() { continue; }
                if confined && !allowed(j, i) { continue; }
                if crosses(i, j, recv) { continue; }
                let v = up(i, bed[j], rise(i, bed[j], area[i], link_len(i, j), lab), cap) + pen(i, j);
                heap.push(Key(v, i as u32, j as u32));
            }
        };
        let run = |heap: &mut BinaryHeap<Key>, bed: &mut [f64], done: &mut [bool], recv: &mut [u32], in_scope: &dyn Fn(usize) -> bool, confined: bool| {
            while let Some(Key(v, i, j)) = heap.pop() {
                let (i, j) = (i as usize, j as usize);
                if done[i] || crosses(i, j, recv) { continue; }
                // B2: the hysteresis orders the heap only; the stored bed is the true path bed (no accumulated penalties).
                done[i] = true; bed[i] = v - pen(i, j); recv[i] = j as u32;
                if path { dest[i].set(dest[j].get()); }
                let (gx, gy) = node_g(i);
                for (dx, dy) in D8 {
                    let Some(k) = node_index(gx + dx, gy + dy) else { continue };
                    if done[k] || !is_land(k) || !in_scope(k) { continue; }
                    if confined && !allowed(i, k) { continue; }
                    if crosses(k, i, recv) { continue; }
                    heap.push(Key(up(k, bed[i], rise(k, bed[i], area[k], link_len(k, i), lab), cap) + pen(k, i), k as u32, i as u32));
                }
            }
        };
        if global {
            // Lane LS1-GLOBAL ls1_global_order_v1: one minimum-path pass over every land node (non-cell modes).
            let mut heap = BinaryHeap::new();
            for i in 0..NODES { if is_land(i) && !done[i] { relax(i, &mut heap, &bed, &done, &recv, true); } }
            run(&mut heap, &mut bed, &mut done, &mut recv, &in_scope, true);
        } else {
            for &(_, c) in &corder {
                let mut heap = BinaryHeap::new();
                for &i in &cells[&c] { relax(i as usize, &mut heap, &bed, &done, &recv, true); }
                run(&mut heap, &mut bed, &mut done, &mut recv, &|k| cell(k) == c, true);
            }
        }
        let unreached = (0..NODES).filter(|&n| is_land(n) && !done[n]).count() as u32;
        st.unreached_confined_per_pass.push(unreached);
        if unreached > 0 {
            let mut heap = BinaryHeap::new();
            for i in 0..NODES { if is_land(i) && !done[i] { relax(i, &mut heap, &bed, &done, &recv, false); } }
            run(&mut heap, &mut bed, &mut done, &mut recv, &|_| true, false);
            // Isolated by refused crossings: the lowest-index settled neighbour (as the 0.9.2 network).
            for n in 0..NODES {
                if !is_land(n) || done[n] { continue; }
                let (gx, gy) = node_g(n);
                if let Some(j) = D8.iter().filter_map(|(dx, dy)| node_index(gx + dx, gy + dy)).filter(|&j| done[j]).min() { recv[n] = j as u32; done[n] = true; }
            }
        }
        // Lane LS1-GLOBAL diagnosis (env VZ_LS1G_DIAG only; stderr, no bytes): why passes keep flipping receivers.
        if dg {
            let (mut nf, mut rev2) = (0u32, 0u32);
            let mut gain = [0u32; 6];
            let mut darea = [0u32; 4];
            for i in 0..NODES {
                if !is_land(i) || recv[i] == net.recv[i] { continue; }
                nf += 1;
                dg_flips[i] = dg_flips[i].saturating_add(1);
                if !dg_prev2.is_empty() && recv[i] == dg_prev2[i] { rev2 += 1; }
                let o = net.recv[i];
                let g = if (o as usize) < NODES && done[o as usize] && bed[o as usize].is_finite() {
                    let o = o as usize;
                    Some(up(i, bed[o], rise(i, bed[o], area[i], link_len(i, o), lab), cap) - bed[i])
                } else { None };
                gain[match g { None => 0, Some(x) if x < 2.0 => 1, Some(x) if x < 5.0 => 2, Some(x) if x < 20.0 => 3, Some(x) if x < 100.0 => 4, _ => 5 }] += 1;
                if !dg_area_prev.is_empty() { let r = (area[i].max(1.0) / dg_area_prev[i].max(1.0)).ln().abs(); darea[if r < 0.01 { 0 } else if r < 0.1 { 1 } else if r < 0.5 { 2 } else { 3 }] += 1; }
            }
            let (mut sc, mut smax) = (0u32, 0f64);
            if !dg_base_prev.is_empty() { for n in 0..NODES { if !is_land(n) && base[n].is_finite() && dg_base_prev[n].is_finite() { let d = (base[n] - dg_base_prev[n]).abs(); if d > 0.01 { sc += 1; smax = smax.max(d); } } } }
            let multi = dg_flips.iter().filter(|&&k| k >= 3).count();
            eprintln!("LS1G diag pass {pass} final {final_pass} settle stage {} global {global}:", sst.as_ref().map(|s| s.stage).unwrap_or(0));
            eprintln!("LS1G diag    flips {nf} back-to-2-ago {rev2} nodes-flipped>=3x {multi}; gain of new over old receiver path [n/a, <2 (order/crossing-blocked), 2-5, 5-20, 20-100, >=100 m] {gain:?}; |ln A/A_prevpass| at flips [<0.01, <0.1, <0.5, >=0.5] {darea:?}; seed beds moved >1 cm {sc} (max {smax:.2} m)");
            dg_prev2 = net.recv.clone(); dg_area_prev = area.clone(); dg_base_prev = base.clone();
        }
        let mut changed = 0u32;
        let diag = std::env::var("VZ_LS_DIAG").is_ok();
        let mut lc = Vec::new();
        for n in 0..NODES { if is_land(n) && recv[n] != net.recv[n] {
            changed += 1;
            if diag && lc.len() < 16 && changed % 7 == 1 && net.recv[n] != NONE && recv[n] != NONE {
                let (gx, gy) = node_g(n); let (ox, oy) = node_g(net.recv[n] as usize); let (nx, ny) = node_g(recv[n] as usize);
                lc.push((gx, gy, ox - gx, oy - gy, nx - gx, ny - gy, bed[n], area[n] / 1e6, 0.0));
            }
        } }
        for n in 0..NODES { if is_land(n) { net.recv[n] = recv[n]; } }
        crate::t1b_fix::rebuild(&mut net, &own, p.meander_permille);
        if incremental { for n in 0..NODES { if is_land(n) && !in_scope(n) { net.area_m2[n] = boundary_area[n]; } } }
        let _ = pass;
        for e in lc.iter_mut() { if let Some(k) = node_index(e.0, e.1) { e.8 = net.area_m2[k] / 1e6; } }
        st.last_changes = lc;
        st.changed_per_pass.push(changed);
        // Lane LS1-GLOBAL ls1_settle_v1 (SETTLE_RULE): stages 1 damped, 2 consistent, 3 fixed areas; a pass that changes no
        // receiver ends LS1 and publishes the areas it used (every pass reads land areas equal to the areas it routes with).
        if let Some(ss) = sst.as_mut() {
            ss.changed_per_stage[ss.stage as usize - 1].push(changed);
            if changed == 0 {
                for n in 0..NODES { if is_land(n) { net.area_m2[n] = area[n]; } }
                ss.converged = true;
                break;
            }
            match ss.stage {
                1 => if changed <= limit || pass + 1 >= ls.passes_max { ss.stage = if s_cons == 0 { 3 } else { 2 }; },
                2 => if ss.changed_per_stage[1].len() as u32 >= s_cons { ss.stage = 3; },
                _ => if ss.changed_per_stage[2].len() as u32 >= s_fixed {
                    // Not settled within the caps: published as the last pass routed it (B1), reported.
                    for n in 0..NODES { if is_land(n) { net.area_m2[n] = area[n]; } }
                    break;
                },
            }
            match ss.stage {
                1 => for n in 0..NODES { if is_land(n) { area[n] = (area[n] * net.area_m2[n]).sqrt(); } },
                2 => for n in 0..NODES { if is_land(n) { area[n] = net.area_m2[n] as f32 as f64; } },
                _ => if ss.changed_per_stage[2].is_empty() { for n in 0..NODES { if is_land(n) { area[n] = net.area_m2[n] as f32 as f64; } } },
            }
            for n in 0..NODES { if is_land(n) { net.area_m2[n] = area[n]; } }
            continue;
        }
        if final_pass {
            for n in 0..NODES { if is_land(n) { net.area_m2[n] = area[n]; } }
            break;
        }
        if changed <= limit || pass + 1 >= ls.passes_max {
            final_pass = true;
            for n in 0..NODES { if is_land(n) { area[n] = net.area_m2[n] as f32 as f64; } }
        } else {
            for n in 0..NODES { if is_land(n) { area[n] = (area[n] * net.area_m2[n]).sqrt(); } }
        }
    }
    st.settle = sst;
    st.receivers_changed_vs_input = (0..NODES).filter(|&n| is_land(n) && in_scope(n) && net.recv[n] != net0.recv[n]).count() as u32;
    // Optimality residual on the published beds (the 0.9.2 bed rule on the final network).
    set_dest(&net);
    let bed = beds_core(&net);
    let labv = labels_of(&net);
    let lab = labv.as_deref();
    let capv = t1a_floor_cap(if cap_cur { &net } else { net0 }, &bed, p);
    let cap = capv.as_deref();
    for i in 0..NODES {
        if !is_land(i) || !in_scope(i) || !bed[i].is_finite() { continue; }
        let (gx, gy) = node_g(i);
        for (dx, dy) in D8 {
            let Some(j) = node_index(gx + dx, gy + dy) else { continue };
            if j == net.recv[i] as usize || !bed[j].is_finite() || !allowed(j, i) { continue; }
            if net0.kind[j] == K_LAND && !net0.conduit[j] {
                // j must not drain through i.
                let mut c = j; let mut cyc = false;
                for _ in 0..NODES { if c == i { cyc = true; break; } if net.kind[c] != K_LAND || net.recv[c] == NONE { break; } c = net.recv[c] as usize; }
                if cyc { continue; }
            }
            let v = up(i, bed[j], rise(i, bed[j], net.area_m2[i], link_len(i, j), lab), cap);
            let ex = bed[i] - v;
            if ex <= 0.01 { continue; }
            if crosses(i, j, &net.recv) { st.residual_crossing += 1; continue; }
            if ex <= hyst { st.residual_in_hysteresis += 1; continue; }
            st.residual_pairs += 1; st.residual_max_m = st.residual_max_m.max(ex);
            st.residual_witnesses.push((i as u32, j as u32, ex));
            let (gxj, gyj) = node_g(j);
            st.residual_top.push((gx, gy, gxj, gyj, ex, net0.kind[j] + if net0.conduit[j] { 10 } else { 0 }));
            if std::env::var("VZ_LS_DIAG").is_ok() && ex > 200.0 && st.residual_detail.len() < 6 {
                let r = net.recv[i] as usize; let (rx, ry) = node_g(r);
                st.residual_detail.push((bed[i], bed[j], v - bed[j], net.area_m2[i] / 1e6, net.area_m2[j] / 1e6, rx - gx, ry - gy, bed[r]));
            }
            if st.residual_top.len() > 64 { st.residual_top.sort_by(|a, b| b.4.total_cmp(&a.4)); st.residual_top.truncate(12); }
        }
    }
    st.residual_top.sort_by(|a, b| b.4.total_cmp(&a.4));
    st.residual_top.truncate(12);
    // Lane LS1-GLOBAL diagnosis (env VZ_LS_IDEM only; stderr, no bytes): LS1 run again on its own output.
    if scope.is_none() && !frozen && std::env::var("VZ_LS_IDEM").is_ok() && !IDEM_GUARD.with(|g| g.get()) {
        IDEM_GUARD.with(|g| g.set(true));
        let (n2, s2) = level_solve_s(net0, Some(&net), ctx, t0, laws, table, p, ls, beds_core, labels_of, false, None);
        IDEM_GUARD.with(|g| g.set(false));
        let rd = (0..NODES).filter(|&n| is_land(n) && n2.recv[n] != net.recv[n]).count();
        let ad = (0..NODES).filter(|&n| is_land(n) && n2.area_m2[n] != net.area_m2[n]).count();
        let b2 = beds_core(&n2);
        let (mut bc, mut bm) = (0u32, 0f64);
        for n in 0..NODES { if bed[n].is_finite() && b2[n].is_finite() { let d = (b2[n] - bed[n]).abs(); if d > 0.01 { bc += 1; bm = bm.max(d); } } }
        // Published land areas against the accumulation of the published receivers (B1 publishes pre-pass areas on rerouted nodes).
        let mut acc = net.clone();
        crate::t1b_fix::rebuild(&mut acc, &own, p.meander_permille);
        let inc = (0..NODES).filter(|&n| is_land(n) && ((net.area_m2[n] - acc.area_m2[n]).abs() > 0.01 * acc.area_m2[n].abs().max(1.0))).count();
        eprintln!("LS1G idem: rerun on the output: passes {} changed per pass {:?}; land receivers differing {rd}, land areas differing {ad}, beds differing > 1 cm {bc} (max {bm:.2} m); residual pairs {} -> {}; published land areas off their receivers' accumulation by > 1 %: {inc}", s2.passes, s2.changed_per_pass, st.residual_pairs, s2.residual_pairs);
    }
    (net, st)
}

/// LS2: spread the bank raise as a ramp at each post's hillslope gradient. `raised` marks the posts the bank rule raised;
/// `grad(i)` is the post's pass-2 hillslope gradient at its pre-ramp height. Returns (posts raised, max raise m).
pub fn bank_ramp(n: usize, dx: f64, z: &mut [f64], free: &dyn Fn(usize) -> bool, raised: &[bool], grad: &dyn Fn(usize) -> f64) -> (u32, f64) {
    #[derive(Clone, Copy)]
    struct M(f64, u32);
    impl PartialEq for M { fn eq(&self, o: &Self) -> bool { self.cmp(o) == Ordering::Equal } }
    impl Eq for M {}
    impl PartialOrd for M { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
    impl Ord for M { fn cmp(&self, o: &Self) -> Ordering { self.0.total_cmp(&o.0).then(o.1.cmp(&self.1)) } }
    let nn = n * n;
    let g: Vec<f64> = (0..nn).map(|i| if free(i) { grad(i) } else { 0.0 }).collect();
    let mut heap = BinaryHeap::new();
    for i in 0..nn { if raised[i] && z[i].is_finite() { heap.push(M(z[i], i as u32)); } }
    let (mut cnt, mut mx) = (0u32, 0f64);
    let mut touched = vec![false; nn];
    while let Some(M(v, i)) = heap.pop() {
        let i = i as usize;
        if v < z[i] { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let k = yy as usize * n + xx as usize;
            if !free(k) || !z[k].is_finite() { continue; }
            let d = if ox != 0 && oy != 0 { dx * std::f64::consts::SQRT_2 } else { dx };
            let t = v - g[k] * d;
            if t > z[k] {
                if !touched[k] { touched[k] = true; cnt += 1; }
                mx = mx.max(t - z[k]);
                z[k] = t;
                heap.push(M(t, k as u32));
            }
        }
    }
    (cnt, mx)
}

/// LS3 state: good Broyden in ln E over the solved cells.
#[derive(Default)]
pub struct Broyden { pub cells: Vec<u32>, j: Vec<Vec<f64>>, prev: Option<(Vec<f64>, Vec<f64>)>, pub singular_fallbacks: u32 }

impl Broyden {
    /// One update: x = ln E, r = M - T, d = max(M - L, 1) per cell (same order as `cells`). Returns the step (clamped).
    pub fn step(&mut self, x: &[f64], r: &[f64], d: &[f64]) -> Vec<f64> {
        let n = x.len();
        if self.j.is_empty() { self.j = (0..n).map(|a| (0..n).map(|b| if a == b { d[a].max(1.0) / 0.6 } else { 0.0 }).collect()).collect(); }
        if let Some((xp, rp)) = &self.prev {
            let dx: Vec<f64> = (0..n).map(|k| x[k] - xp[k]).collect();
            let dr: Vec<f64> = (0..n).map(|k| r[k] - rp[k]).collect();
            let dd: f64 = dx.iter().map(|v| v * v).sum();
            if dd > 1e-12 {
                let jdx: Vec<f64> = (0..n).map(|a| (0..n).map(|b| self.j[a][b] * dx[b]).sum()).collect();
                for a in 0..n { for b in 0..n { self.j[a][b] += (dr[a] - jdx[a]) * dx[b] / dd; } }
            }
        }
        self.prev = Some((x.to_vec(), r.to_vec()));
        let rhs: Vec<f64> = r.iter().map(|v| -v).collect();
        let s = match solve_dense(&self.j, &rhs) {
            Some(s) if s.iter().all(|v| v.is_finite()) => s,
            _ => { self.singular_fallbacks += 1; (0..n).map(|k| -r[k] * 0.6 / d[k].max(1.0)).collect() }
        };
        let l2 = 2f64.ln();
        s.iter().map(|v| v.clamp(-l2, l2)).collect()
    }
}

fn solve_dense(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut m: Vec<Vec<f64>> = (0..n).map(|i| { let mut r = a[i].clone(); r.push(b[i]); r }).collect();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| m[i][c].abs().total_cmp(&m[j][c].abs()).then(j.cmp(&i)))?;
        if m[p][c].abs() < 1e-9 { return None; }
        m.swap(c, p);
        for i in 0..n {
            if i == c { continue; }
            let f = m[i][c] / m[c][c];
            if f == 0.0 { continue; }
            for k in c..=n { m[i][k] -= f * m[c][k]; }
        }
    }
    Some((0..n).map(|i| m[i][n] / m[i][i]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn broyden_solves_a_coupled_linear_pair() {
        // r = A x - b with strong off-diagonal coupling (the 41995 shape: cell 1 responds mostly to x0).
        let a = [[200.0, 20.0], [150.0, 40.0]];
        let b = [100.0, 80.0];
        let mut x = vec![0.0, 0.0];
        let mut br = Broyden::default();
        for _ in 0..30 {
            let r: Vec<f64> = (0..2).map(|i| a[i][0] * x[0] + a[i][1] * x[1] - b[i]).collect();
            if r.iter().all(|v| v.abs() < 1e-6) { break; }
            let s = br.step(&x, &r, &[200.0, 40.0]);
            for k in 0..2 { x[k] += s[k]; }
        }
        let r: Vec<f64> = (0..2).map(|i| a[i][0] * x[0] + a[i][1] * x[1] - b[i]).collect();
        assert!(r.iter().all(|v| v.abs() < 1e-3), "residual {r:?}");
    }
    #[test]
    fn tb_absent_keeps_the_ls_params_bytes() {
        // Lane TB: without options the sub-struct is not serialized (recipe digests of 0.9.4 unchanged).
        let v = serde_json::to_value(LsParamsV1::v1()).unwrap();
        assert!(v.get("tb").is_none());
        let mut l = LsParamsV1::v1();
        l.tb = Some(TbParamsV1::off());
        l.tb.as_mut().unwrap().t1a_levels = 3;
        assert!(l.validate().is_err());
        l.tb.as_mut().unwrap().t1a_levels = 1;
        assert!(l.validate().is_ok());
        assert!(serde_json::to_value(&l).unwrap().get("tb").is_some());
    }
    #[test]
    fn chain_labels_follow_curves() {
        // a -> b -> c -> sea, d -> b: chains start at heads (a, d) and at the junction b; c carries b's label.
        use crate::t1b_rivernet::K_MARINE;
        let idx = |gx: i64, gy: i64| node_index(gx, gy).unwrap();
        let (a, d, b, c, m) = (idx(0, 2), idx(1, 2), idx(0, 1), idx(0, 0), idx(0, -1));
        let mut net = Network { kind: vec![K_MARINE; NODES], recv: vec![NONE; NODES], order: vec![], area_m2: vec![0.0; NODES], pos: vec![[0.0; 2]; NODES],
            mid: vec![[0.0; 2]; NODES], term: vec![NONE; NODES], level: vec![0.0; NODES], unreached_confined: 0, nodes_land: 4, rockmass: vec![1.0; NODES],
            conduit: vec![false; NODES], conduit_cap: vec![f64::NAN; NODES], conduit_paths: vec![], t1a_recv: vec![NONE; NODES], t0_cell: vec![0; NODES], meander_u: vec![0.5; NODES], t1a_ch: vec![false; NODES], conf_cell: vec![], flux: None };
        for (n, r) in [(a, b), (d, b), (b, c), (c, m)] { net.kind[n] = K_LAND; net.recv[n] = r as u32; }
        let chan: Vec<bool> = (0..NODES).map(|n| net.kind[n] == K_LAND).collect();
        let lab = chain_labels(&net, &chan);
        assert_eq!((lab[a], lab[d], lab[b], lab[c], lab[m]), (a as u32, d as u32, b as u32, b as u32, NONE));
    }
    #[test]
    fn tc_absent_keeps_the_tb_bytes_and_validates() {
        // Lane TC: without options the sub-struct is not serialized (TB recipe digests unchanged); clip needs ls_t1a.
        assert!(serde_json::to_value(TbParamsV1::off()).unwrap().get("tc").is_none());
        let mut l = LsParamsV1::v1();
        let mut t = TbParamsV1::off();
        let mut c = TcParamsV1::off();
        c.clip = 1;
        t.tc = Some(c);
        l.tb = Some(t);
        assert!(l.validate().is_err());
        l.tb.as_mut().unwrap().t1a_levels = 2;
        assert!(l.validate().is_ok());
        l.tb.as_mut().unwrap().tc.as_mut().unwrap().pert_ppm = 10_000;
        assert!(l.validate().is_err());
    }
    #[test]
    fn t1a_floor_cap_grades_a_node_beside_the_trunk() {
        // T1a link t (0,0) -> r (1,0), beds 100 -> 90; land node n at (0,1) is 512 m from the link (inside h = 600 m):
        // cap = 100 + max(0.005, 10/512) x 512 = 110; a node at (0,2) (1024 m) is outside.
        use crate::t1b_rivernet::{K_MARINE, K_T1A};
        let idx = |gx: i64, gy: i64| node_index(gx, gy).unwrap();
        let (t, r, n, f) = (idx(0, 0), idx(1, 0), idx(0, 1), idx(0, 2));
        let lat = |k: usize| { let (gx, gy) = node_g(k); [(gx * SPACING_MM) as f64 / 1000.0, (gy * SPACING_MM) as f64 / 1000.0] };
        let mut net = Network { kind: vec![K_MARINE; NODES], recv: vec![NONE; NODES], order: vec![], area_m2: vec![0.0; NODES], pos: (0..NODES).map(lat).collect(),
            mid: vec![[0.0; 2]; NODES], term: vec![NONE; NODES], level: vec![0.0; NODES], unreached_confined: 0, nodes_land: 2, rockmass: vec![1.0; NODES],
            conduit: vec![false; NODES], conduit_cap: vec![f64::NAN; NODES], conduit_paths: vec![], t1a_recv: vec![NONE; NODES], t0_cell: vec![0; NODES], meander_u: vec![0.5; NODES], t1a_ch: vec![false; NODES], conf_cell: vec![], flux: None };
        for k in [t, r] { net.kind[k] = K_T1A; net.t1a_ch[k] = true; net.area_m2[k] = 1e12; }
        net.t1a_recv[t] = r as u32;
        for k in [n, f] { net.kind[k] = K_LAND; }
        let mut bed = vec![f64::NAN; NODES];
        bed[t] = 100.0; bed[r] = 90.0;
        let mut p = RiverNetParamsV1::v1(false);
        let mut l = LsParamsV1::v1();
        let mut tb = TbParamsV1::off();
        let mut c = TcParamsV1::off();
        c.floor = 2; tb.tc = Some(c); l.tb = Some(tb); p.ls = Some(l);
        let cap = t1a_floor_cap(&net, &bed, &p).unwrap();
        assert!((cap[n] - 110.0).abs() < 1e-9, "{}", cap[n]);
        assert!(cap[f].is_infinite());
    }
    #[test]
    fn bank_ramp_spreads_a_raise() {
        let n = 5;
        let mut z = vec![0.0; n * n];
        z[12] = 10.0;
        let mut raised = vec![false; n * n];
        raised[12] = true;
        let (c, m) = bank_ramp(n, 1.0, &mut z, &|i| i != 12, &raised, &|_| 4.0);
        assert!(c > 0 && (m - 6.0).abs() < 1e-9);
        assert!((z[13] - 6.0).abs() < 1e-9 && (z[14] - 2.0).abs() < 1e-9 && z[0] == 0.0);
    }
    #[test]
    fn ls_global_absent_keeps_the_params_bytes_and_parses() {
        // Lane LS1-GLOBAL: absent = no recipe field; ls_global needs LS1 and is a no-op under the cell rule; ls_settle caps
        // need ls_settle=1 first; both off removes the struct again.
        let mut p = RiverNetParamsV1::v1(false);
        let off = serde_json::to_value(&p).unwrap();
        assert!(off.get("ls_global").is_none());
        assert!(!parse_key("ls_globalx=1", &mut p).unwrap());
        assert!(parse_key("ls_global=2", &mut p).is_err());
        assert!(parse_key("ls_global=1", &mut p).unwrap());
        assert_eq!(p.ls_global.as_ref().unwrap().order_rule, GLOBAL_RULE);
        assert!(validate_global(&p).is_err(), "needs ls=1");
        p.ls = Some(LsParamsV1::v1());
        assert!(validate_global(&p).is_ok());
        assert!(!global_order(&p), "cell rule (absent confinement): per-cell order kept");
        p.confinement = Some(crate::t1b_outlet::ConfinementParamsV1::v1(crate::t1b_outlet::CONF_CELL).unwrap());
        assert!(!global_order(&p), "conf=cell: per-cell order kept");
        p.confinement = Some(crate::t1b_outlet::ConfinementParamsV1::v1(crate::t1b_outlet::CONF_NONE).unwrap());
        assert!(global_order(&p));
        assert!(parse_key("ls_settle_cons=4", &mut p).is_err(), "needs ls_settle=1");
        assert!(parse_key("ls_settle=1", &mut p).unwrap());
        assert!(settle(&p));
        let g = p.ls_global.as_ref().unwrap();
        assert_eq!((g.settle_consistent_max, g.settle_fixed_max, g.settle_rule.as_str()), (SETTLE_CONSISTENT_MAX, SETTLE_FIXED_MAX, SETTLE_RULE));
        assert!(parse_key("ls_settle_fixed=0", &mut p).is_err());
        assert!(parse_key("ls_settle_cons=0", &mut p).unwrap());
        assert_eq!(p.ls_global.as_ref().unwrap().settle_consistent_max, 0);
        assert!(parse_key("ls_global=0", &mut p).unwrap());
        assert!(p.ls_global.is_some() && !global_order(&p));
        assert!(parse_key("ls_settle=0", &mut p).unwrap());
        assert!(p.ls_global.is_none());
        let mut q = p.clone();
        q.ls = None;
        q.confinement = None;
        assert_eq!(serde_json::to_value(&q).unwrap(), off);
    }
    #[test]
    fn undercut_off_keeps_the_tb_bytes_and_needs_soft_channels() {
        // Lane GUIDE-SPIKES: absent / 0 serializes nothing new; 1 carries the rule; 1 without steep_eik 2 is refused.
        let off = serde_json::to_value(TbParamsV1::off()).unwrap();
        assert!(off.get("undercut").is_none() && off.get("undercut_rule").is_none());
        let mut l = LsParamsV1::v1();
        let mut t = TbParamsV1::off();
        t.set_undercut(0);
        assert_eq!(serde_json::to_value(&t).unwrap(), off);
        t.set_undercut(1);
        assert_eq!(serde_json::to_value(&t).unwrap()["undercut_rule"], UNDERCUT_RULE);
        l.tb = Some(t.clone());
        assert!(l.validate().is_err(), "undercut needs steep_eik 2");
        t.steep_eik = 2;
        l.tb = Some(t);
        assert!(l.validate().is_ok());
    }
}

/// LS1 (windows): the layer's published receivers and areas (B1: the areas the last pass used), rebuilt as the solve does.
pub fn with_published_network(net0: &Network, layer_dir: &std::path::Path, p: &RiverNetParamsV1) -> R<Network> {
    let mut net = crate::t1b_fix::with_published_receivers(net0, layer_dir, p)?;
    let b = std::fs::read(layer_dir.join("fields/area_m2.f32")).map_err(|e| format!("area_m2.f32: {e}"))?;
    if b.len() != NODES * 4 { return Err("area_m2.f32 size".into()); }
    for (n, c) in b.chunks_exact(4).enumerate() { if net.kind[n] == K_LAND && !net.conduit[n] { net.area_m2[n] = f32::from_le_bytes([c[0], c[1], c[2], c[3]]) as f64; } }
    Ok(net)
}

/// Explicit current-ground routing, separate from LS1's equilibrium bed construction.
/// Barnes, Lehman & Mulla (2014), doi:10.1016/j.cageo.2013.04.024, priority flood.
/// Fill is a DEM conditioning operation, NOT sediment creation or water supplied.
/// `routing_z_m` must accompany the network. Never recompute its beds with the
/// equilibrium bed law. The caller must publish this record and the fill depths.
pub const SURFACE_RULE: &str = "ls_surface_fill_v1: chart-wide constrained priority flood on supplied heights; fill only, no carve; T1a/water are immutable boundary seeds; conduit receivers are locked but conduit heights participate in fill; no unconfined fallback; flat edges follow decreasing canonical flood rank; heap ties by native node then receiver ID; volume=sum(fill_depth*supplied support area); filled volume is storage capacity, not water supplied or deposited solid; local contributing areas terminate at declared boundaries; no SimTime evolution";

#[derive(Clone, Debug, Serialize)]
pub struct SurfaceBoundary {
    pub node: u32,
    pub kind: u8,
    pub conduit: bool,
    pub t0_cell: u32,
    pub local_inflow_m2: f64,
    pub effective_inflow_m2: f64,
    pub fill_capacity_m3: f64,
}

/// Route plus its mandatory conditioned surface and accounting. Native node IDs
/// and the full chart are canonical, independent of window evaluation order.
/// The land links are D8 straight segments; no keyed bend is certified by this
/// route (a bend could cross higher ground between its endpoints).
pub struct SurfaceRoute {
    pub network: Network,
    pub routing_z_m: Vec<f64>,
    pub fill_depth_m: Vec<f64>,
    pub flat_rank: Vec<u32>,
    /// Local catchment area only. No upstream T1a discharge is counted twice.
    pub local_area_m2: Vec<f64>,
    /// Native runoff-effective area, including declared upstream conduit inputs.
    pub effective_area_m2: Vec<f64>,
    pub report: SurfaceReport,
}

#[derive(Clone, Debug, Serialize)]
pub struct SurfaceReport {
    pub rule: String,
    pub active_nodes: usize,
    pub filled_nodes: usize,
    pub fill_capacity_m3: f64,
    pub max_fill_m: f64,
    pub input_local_area_m2: f64,
    pub boundary_local_area_m2: f64,
    pub input_effective_source_m2: f64,
    pub boundary_effective_area_m2: f64,
    pub flat_edges: usize,
    pub boundaries: Vec<SurfaceBoundary>,
}

/// Opt-in API mode for Stage 1 and terrain edits. Every height is in metres and
/// every support/source area in m2. Supports are immutable geometric quadrature weights,
/// NOT the old accumulated drainage areas. `source_area_m2` holds the native
/// runoff-weighted own areas and prescribed conduit injections. Obtain it once
/// from `surface_native_sources` on the original unrounded built network and
/// retain it as immutable input, including for reruns. T1a/water heights are supplied
/// boundary data; their topology is retained. Conduit receiver links stay locked,
/// but their actual supplied heights are filled and accounted like other land.
/// A changed fixed boundary is the caller's model change, not silently repaired.
/// The result certifies ALL land edges, including conduits, not fixed T1a edges.
///
/// Depression policy is explicitly FillToSpill. It computes capacity assuming
/// impermeable ground and sustained supply, not a transient lake simulation.
/// Fill depth is published separately from physical ground. Actual ground can
/// rise along a submerged link; the declared routing surface cannot.
/// Repeating on routing_z_m gives identical routing/areas/heights, zero NEW fill.
/// The original fill record remains provenance and must not be discarded.
pub fn route_surface_fill(net: &Network, t0: &T0, p: &RiverNetParamsV1,
    z_m: &[f64], support_area_m2: &[f64], source_area_m2: &[f64]) -> R<SurfaceRoute> {
    if [net.kind.len(), net.recv.len(), net.pos.len(), net.conduit.len(),
        net.t0_cell.len(), net.t1a_recv.len(), net.area_m2.len(), net.meander_u.len(),
        net.mid.len(), net.term.len(), net.level.len(), net.order.len()]
        .iter().any(|&n| n != NODES) { return Err("surface route: invalid native network dimensions".into()); }
    if z_m.len() != NODES || support_area_m2.len() != NODES || source_area_m2.len() != NODES
        || z_m.iter().any(|z| !z.is_finite())
        || support_area_m2.iter().any(|a| !a.is_finite() || *a <= 0.0)
        || source_area_m2.iter().any(|a| !a.is_finite() || *a < 0.0) {
        return Err("surface route: finite native heights and positive finite support areas and nonnegative source areas required".into());
    }
    if let Some(c) = &p.confinement {
        if c.mode > crate::t1b_outlet::CONF_FREE { return Err("surface route: unknown confinement mode".into()); }
        if (c.mode == crate::t1b_outlet::CONF_FREE || c.mode == crate::t1b_outlet::CONF_WARP) && net.conf_cell.len() != NODES {
            return Err("surface route: free/warp confinement requires its native zone/cell field".into());
        }
    }
    let active: Vec<bool> = (0..NODES).map(|i| net.kind[i] == K_LAND).collect();
    for i in 0..NODES {
        if !active[i] && net.kind[i] != K_T1A && net.kind[i] != crate::t1b_rivernet::K_MARINE
            && net.kind[i] != crate::t1b_rivernet::K_LAKE {
            return Err(format!("surface route: unsupported boundary kind {} at {i}", net.kind[i]));
        }
        if net.recv[i] != NONE && net.recv[i] as usize >= NODES {
            return Err(format!("surface route: invalid receiver at {i}"));
        }
    }
    let mut seen = vec![false; NODES];
    for &i in &net.order {
        if i as usize >= NODES || seen[i as usize] { return Err("surface route: invalid input network order".into()); }
        seen[i as usize] = true;
    }
    // Fixed boundaries must not return into the active domain: treating such a
    // boundary as an outlet would hide a cycle after the reroute.
    for i in 0..NODES { if !active[i] && net.recv[i] != NONE && active[net.recv[i] as usize] {
        return Err(format!("surface route: fixed boundary {i} drains into the active domain"));
    } }
    let conf = crate::t1b_outlet::Confinement::new(p, t0, &net.t0_cell, &net.conf_cell, &net.t1a_recv);
    let fixed_dest: Vec<u32> = (0..NODES).map(|i| i as u32).collect();
    let allowed = |i: usize, j: usize, destination: u32| {
        if net.conduit[i] { return net.recv[i] == j as u32; }
        if net.kind[j] == K_T1A && conf.is_basin() { return conf.allows_nodes(t0, i, j, true); }
        if !active[j] || net.conduit[j] { return true; }
        if conf.path_mode() { return conf.path_ok(t0, &net.kind, i, destination); }
        conf.allows_nodes(t0, i, j, false)
    };
    let (height, recv, rank, order, dest) = surface_flood(N, z_m, &active, &net.recv, &fixed_dest, &net.conduit, &allowed)?;
    let mut area: Vec<f64> = (0..NODES).map(|i| if active[i] { support_area_m2[i] } else { 0.0 }).collect();
    let mut capacity = vec![0.0; NODES];
    let mut depth = vec![0.0; NODES];
    let mut report = SurfaceReport { rule: SURFACE_RULE.into(), active_nodes: 0, filled_nodes: 0,
        fill_capacity_m3: 0.0, max_fill_m: 0.0, input_local_area_m2: 0.0,
        boundary_local_area_m2: 0.0, input_effective_source_m2: 0.0, boundary_effective_area_m2: 0.0, flat_edges: 0, boundaries: Vec::new() };
    for i in 0..NODES { if active[i] {
        report.active_nodes += 1;
        report.input_local_area_m2 += support_area_m2[i];
        depth[i] = height[i] - z_m[i];
        if depth[i] > 0.0 { report.filled_nodes += 1; }
        capacity[i] = depth[i] * support_area_m2[i];
        report.fill_capacity_m3 += capacity[i];
        report.max_fill_m = report.max_fill_m.max(depth[i]);
        let j = recv[i] as usize;
        if height[i] < height[j] || rank[i] <= rank[j] { return Err(format!("surface route: non-descending potential at {i}")); }
        if height[i] == height[j] { report.flat_edges += 1; }
    } }
    for &i in order.iter().rev() { let i = i as usize; if active[i] {
        let j = recv[i] as usize; area[j] += area[i]; capacity[j] += capacity[i];
    } }
    for i in 0..NODES { if !active[i] && area[i] > 0.0 {
        report.boundary_local_area_m2 += area[i];
        report.boundaries.push(SurfaceBoundary { node: i as u32, kind: net.kind[i], conduit: net.conduit[i],
            t0_cell: net.t0_cell[i], local_inflow_m2: area[i], effective_inflow_m2: 0.0, fill_capacity_m3: capacity[i] });
    } }
    if !report.fill_capacity_m3.is_finite() || !report.input_local_area_m2.is_finite()
        || (report.input_local_area_m2 - report.boundary_local_area_m2).abs() > 1e-10 * report.input_local_area_m2.max(1.0) {
        return Err("surface route: non-finite volume or local area conservation failure".into());
    }
    let mut routed = net.clone();
    for i in 0..NODES { if active[i] {
        routed.recv[i] = recv[i]; routed.area_m2[i] = area[i]; routed.term[i] = dest[i];
        let j = recv[i] as usize;
        routed.mid[i] = [0.5 * (net.pos[i][0] + net.pos[j][0]), 0.5 * (net.pos[i][1] + net.pos[j][1])];
        routed.level[i] = height[i];
    } }
    // Canonical seed order is independent of the caller’s visitation order.
    routed.order = (0..NODES).filter(|&i| !active[i]).map(|i| i as u32).chain(order.iter().copied()).collect();
    routed.unreached_confined = 0;
    // Geometric support for fill and local catchments must not replace the
    // climate-weighted native source field or silently lose upstream discharge.
    let mut sources = source_area_m2.to_vec();
    if routed.flux.is_some() {
        let (seeds, diag) = crate::t1b_outlet::flux_seed_areas(&routed, source_area_m2);
        if diag.seeds_floored_at_zero != 0 || diag.unknown_destination_nodes != 0 {
            return Err("surface route: native flux adjustment has negative seeds or unknown destinations".into());
        }
        for path in &routed.conduit_paths { if path.len() >= 2 {
            let (t,c)=(path[0] as usize,path[1] as usize);
            if routed.kind[t] == K_T1A && active[c] { sources[c] += seeds[t] - source_area_m2[t]; }
        } }
        for i in 0..NODES { if routed.kind[i] == K_T1A { routed.area_m2[i] = seeds[i]; } }
    }
    if sources.iter().any(|v| !v.is_finite() || *v < 0.0) { return Err("surface route: negative/nonfinite effective source after boundary adjustment".into()); }
    let mut effective: Vec<f64> = (0..NODES).map(|i| if active[i] {sources[i]} else {0.0}).collect();
    report.input_effective_source_m2 = effective.iter().sum();
    for &i in order.iter().rev() { let i=i as usize; let j=recv[i] as usize; effective[j] += effective[i]; }
    for i in 0..NODES {
        if active[i] { routed.area_m2[i] = effective[i]; }
        else { report.boundary_effective_area_m2 += effective[i]; }
    }
    for b in &mut report.boundaries { b.effective_inflow_m2 = effective[b.node as usize]; }
    if !report.input_effective_source_m2.is_finite() || !report.boundary_effective_area_m2.is_finite()
        || (report.input_effective_source_m2 - report.boundary_effective_area_m2).abs() > 1e-10 * report.input_effective_source_m2.max(1.0) {
        return Err("surface route: effective area conservation failure".into());
    }
    Ok(SurfaceRoute { network: routed, routing_z_m: height, fill_depth_m: depth, flat_rank: rank, local_area_m2: area, effective_area_m2: effective, report })
}

/// Immutable own-area input for surface routing. Call on the original unrounded
/// `build_network` result, not a published f32 LS1 area field or accumulated output.
/// Carries climate weights and external T1a-to-conduit injections separately from
/// the geometric support weights used to measure fill volumes.
pub fn surface_native_sources(built: &Network) -> Vec<f64> { crate::t1b_fix::own_area(built) }

type SurfaceFlood = (Vec<f64>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>);

/// The queue is ordered by (filled elevation, native node, receiver). All seeds
/// enter at their supplied level. Refused diagonals are never replaced with an
/// arbitrary neighbour. An unreachable component is an error, not a new outlet.
fn surface_flood(n: usize, z: &[f64], active: &[bool], fixed_recv: &[u32], fixed_dest: &[u32], locked: &[bool],
    allowed: &dyn Fn(usize, usize, u32) -> bool) -> R<SurfaceFlood> {
    let count = n * n;
    let mut height = z.to_vec();
    let mut done: Vec<bool> = active.iter().map(|v| !v).collect();
    let mut recv = fixed_recv.to_vec();
    for i in 0..count { if active[i] && !locked[i] { recv[i] = NONE; } }
    let mut dest = fixed_dest.to_vec();
    let mut rank = vec![0; count];
    let mut order = Vec::new();
    let mut heap = BinaryHeap::new();
    let neighbours = |i: usize| -> Vec<usize> {
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        D8.iter().filter_map(|&(dx, dy)| { let (x, y) = (x + dx, y + dy);
            (x >= 0 && y >= 0 && x < n as i64 && y < n as i64).then_some((y * n as i64 + x) as usize)
        }).collect()
    };
    let crosses = |i: usize, j: usize, recv: &[u32]| {
        if i % n == j % n || i / n == j / n { return false; }
        let (a, b) = (i / n * n + j % n, j / n * n + i % n);
        recv[a] == b as u32 || recv[b] == a as u32
    };
    for j in 0..count { if done[j] {
        for i in neighbours(j) { if active[i] && allowed(i, j, dest[j]) && !crosses(i, j, &recv) {
            heap.push(Key(z[i].max(height[j]), i as u32, j as u32));
        } }
    } }
    while let Some(Key(h, i, j)) = heap.pop() {
        let (i, j) = (i as usize, j as usize);
        if done[i] || crosses(i, j, &recv) { continue; }
        done[i] = true; height[i] = h; recv[i] = j as u32;
        dest[i] = dest[j]; rank[i] = order.len() as u32 + 1; order.push(i as u32);
        for k in neighbours(i) { if !done[k] && allowed(k, i, dest[i]) && !crosses(k, i, &recv) {
            heap.push(Key(z[k].max(h), k as u32, i as u32));
        } }
    }
    let unreached = done.iter().filter(|v| !**v).count();
    if unreached > 0 { return Err(format!("surface route: {unreached} nodes unreachable under confinement/outlet/crossing constraints; no fallback")); }
    Ok((height, recv, rank, order, dest))
}

#[cfg(test)]
mod surface_tests {
    use super::*;
    #[test]
    fn surface_crater_fill_known_volume_and_exact_rerun() {
        let n = 9;
        let mut z = vec![10.0; n*n];
        let active: Vec<bool> = (0..n*n).map(|i| i%n > 0 && i%n < n-1 && i/n > 0 && i/n < n-1).collect();
        for y in 3..6 { for x in 3..6 { z[y*n+x] = -490.0; } }
        let fixed = vec![NONE; n*n]; let dest: Vec<u32> = (0..n*n).map(|i| i as u32).collect();
        let (h, r, rank, _, d) = surface_flood(n, &z, &active, &fixed, &dest, &vec![false;n*n], &|_,_,_|true).unwrap();
        assert_eq!(h.iter().zip(&z).map(|(a,b)|a-b).sum::<f64>() * 250000.0, 1125000000.0);
        for i in 0..n*n { if active[i] { let j=r[i] as usize; assert!(h[i]>=h[j] && rank[i]>rank[j]); } }
        let (h2,r2,rank2,_,d2) = surface_flood(n,&h,&active,&r,&dest,&vec![false;n*n],&|_,_,_|true).unwrap();
        assert_eq!((h,r,rank,d),(h2,r2,rank2,d2));
    }
    #[test]
    fn surface_closed_or_forbidden_basin_fails_without_fallback() {
        let n=3; let z = vec![0.0;9]; let fixed=vec![NONE;9]; let dest=vec![0;9];
        assert!(surface_flood(3,&z,&[true;9],&fixed,&dest,&vec![false;n*n],&|_,_,_|true).is_err());
        let mut active=vec![true;9]; active[0]=false;
        assert!(surface_flood(3,&z,&active,&fixed,&dest,&vec![false;n*n],&|_,_,_|false).is_err());
    }
    #[test]
    fn surface_locked_conduit_fills_uphill_ground_without_changing_links() {
        let n=3; let z=vec![0.0,8.0,0.0,0.0,2.0,0.0,0.0,-3.0,0.0];
        let mut active=vec![false;9];active[4]=true;active[7]=true;
        let mut fixed=vec![NONE;9];fixed[4]=1;fixed[7]=4;
        let locked=active.clone();let dest:Vec<u32>=(0..9).map(|i|i as u32).collect();
        let (h,r,rank,_,_)=surface_flood(n,&z,&active,&fixed,&dest,&locked,&|i,j,_|fixed[i]==j as u32).unwrap();
        assert_eq!((h[4],h[7]),(8.0,8.0));assert_eq!(r,fixed);assert!(rank[7]>rank[4]);
        let again=surface_flood(n,&h,&active,&r,&dest,&locked,&|i,j,_|fixed[i]==j as u32).unwrap();
        assert_eq!(again.0,h);assert_eq!(again.1,r);
    }
    #[test]
    fn surface_path_constraint_selects_only_declared_terminal() {
        let n=5; let z=vec![0.0;n*n];let mut active=vec![true;n*n];active[0]=false;active[n*n-1]=false;
        let fixed=vec![NONE;n*n];let dest:Vec<u32>=(0..n*n).map(|i|i as u32).collect();
        let (_,recv,_,_,out)=surface_flood(n,&z,&active,&fixed,&dest,&vec![false;n*n],&|_,_,terminal|terminal==24).unwrap();
        for i in 0..n*n {if active[i] {assert_eq!(out[i],24);assert_ne!(recv[i],0);}}
    }
    #[test]
    fn surface_plateau_is_deterministic_and_no_crossing() {
        let n=16; let z=vec![7.0;n*n]; let mut active=vec![true;n*n]; active[0]=false;active[n*n-1]=false;
        let fixed=vec![NONE;n*n];let dest:Vec<u32>=(0..n*n).map(|i|i as u32).collect();
        let a=surface_flood(n,&z,&active,&fixed,&dest,&vec![false;n*n],&|_,_,_|true).unwrap();
        let b=surface_flood(n,&z,&active,&fixed,&dest,&vec![false;n*n],&|_,_,_|true).unwrap(); assert_eq!(a,b);
        for y in 0..n-1 {for x in 0..n-1 {let i=y*n+x;let r=&a.1;
            assert!(!((r[i]==(i+n+1) as u32 || r[i+n+1]==i as u32) && (r[i+1]==(i+n) as u32 || r[i+n]==(i+1) as u32)));
        }}
    }
}

/// Globally addressed, closed detail components. The legacy 513² wrapper above
/// intentionally retains its exact policy and reduction order.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct DetailEdge { pub receiver: u32, pub weight: f64 }
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct DetailBalance {
    pub basin: u32,
    pub geometric_source_m2: f64,
    pub effective_source_m2: f64,
    pub geometric_export_m2: f64,
    pub effective_export_m2: f64,
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct DetailBasinTransfer {
    pub from_basin:u32, pub to_basin:u32,
    pub geometric_m2:f64, pub effective_m2:f64,
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct DetailRoute {
    pub domain_digest: String,
    pub kernel_digest: String,
    pub mode: crate::t1b_surface_domain::Mode,
    pub routing_m: Vec<f64>,
    pub parent: Vec<u32>,
    pub rank: Vec<u32>,
    pub downstream_first: Vec<u32>,
    pub edges: Vec<Vec<DetailEdge>>,
    pub geometric_area_m2: Vec<f64>,
    pub effective_area_m2: Vec<f64>,
    pub fill_capacity_m3: f64,
    pub balances: Vec<DetailBalance>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub basin_transfers: Vec<DetailBasinTransfer>,
    pub digest: String,
}

/// The terrain remains untouched. H is a fill-to-spill potential; capacity is
/// neither water supplied nor sediment. Construction can explicitly adopt H as
/// an initial condition, retaining the raise volume and predecessor pin.
pub fn route_detail_surface(domain: &crate::t1b_surface_domain::Domain, mode: crate::t1b_surface_domain::Mode)
    -> crate::t1b_surface_domain::Result<DetailRoute> {
    use crate::t1b_surface_domain::{invalid,content_digest,Mode,Policy,PostId};
    use std::collections::{BTreeMap,BTreeSet};
    let posts=&domain.posts;
    let count=posts.len();
    if count>u32::MAX as usize {return Err(invalid("detail domain exceeds u32 indexing"));}
    if domain.policy==Policy::PinnedConstructionV1 && mode!=Mode::SingleFlowV1 {return Err(invalid("construction minimax uses SFD; MFD requires the native detail policy"));}
    let mut h:Vec<f64>=posts.iter().map(|p|p.ground_m).collect();
    let mut done:Vec<bool>=posts.iter().map(|p|p.boundary.is_some()).collect();
    let mut parent=vec![NONE;count];
    let mut rank=vec![0;count];
    let mut order=Vec::with_capacity(count);
    let mut heap=BinaryHeap::new();
    let crosses=|i:usize,j:usize,parent:&[u32]| {
        // Crossing drawn diagonals are a legacy network constraint, not an
        // extra wall in the physical D8 multiple-flow discretization.
        if domain.policy!=Policy::NativeDetailV1 {return false;}
        let (p,q)=(posts[i].id,posts[j].id);
        if p.gx==q.gx || p.gy==q.gy {return false;}
        let a=domain.index(PostId {gy:p.gy,gx:q.gx});
        let b=domain.index(PostId {gy:q.gy,gx:p.gx});
        match (a,b) {(Some(a),Some(b))=>parent[a]==b as u32 || parent[b]==a as u32,_=>false}
    };
    // Reverse adjacency is computed from actual legal outgoing edges. Native
    // restrictions may be asymmetric, notably on locked conduit links.
    let donors=|j:usize|->Vec<usize> {
        let mut v=Vec::new();
        for k in 0..8 {
            if let Some(i)=domain.index(posts[j].id.neighbour(k)) {
                if domain.edge(i,(k+4)%8)==Some(j) {v.push(i);}
            }
        }
        v.sort_unstable();v
    };
    for j in 0..count {if done[j] {
        for i in donors(j) {heap.push(Key(h[i].max(h[j]),i as u32,j as u32));}
    }}
    while let Some(Key(level,i,j))=heap.pop() {
        let (i,j)=(i as usize,j as usize);
        if done[i] || crosses(i,j,&parent) {continue;}
        done[i]=true;h[i]=level;parent[i]=j as u32;
        rank[i]=order.len() as u32+1;order.push(i as u32);
        for k in donors(i) {if !done[k] && !crosses(k,i,&parent) {
            heap.push(Key(h[k].max(level),k as u32,i as u32));
        }}
    }
    let unreached=done.iter().filter(|x|!**x).count();
    if unreached!=0 {return Err(invalid(format!("{unreached} detail posts cannot reach a supported boundary under native constraints")));}
    let pair=|a:PostId,b:PostId|if a<b {(a,b)} else {(b,a)};
    let mut diagonals=BTreeSet::new();
    for i in 0..count {if domain.policy==Policy::NativeDetailV1 && parent[i]!=NONE {
        let (a,b)=(posts[i].id,posts[parent[i] as usize].id);
        if a.gx!=b.gx && a.gy!=b.gy {diagonals.insert(pair(a,b));}
    }}
    let mut edges=vec![Vec::new();count];
    for &ii in &order {
        let i=ii as usize;
        let p=parent[i] as usize;
        if h[i]<h[p] || (h[i]==h[p] && rank[i]<=rank[p]) {return Err(invalid("non-descending detail potential"));}
        if mode==Mode::MultipleFlowSlope11V1 && posts[i].locked_receiver.is_none() {
            for k in 0..8 {if let Some(j)=domain.edge(i,k) {
                if h[i]<=h[j] {continue;}
                let (a,b)=(posts[i].id,posts[j].id);
                let diagonal=a.gx!=b.gx && a.gy!=b.gy;
                if domain.policy==Policy::NativeDetailV1 && diagonal && diagonals.contains(&pair(PostId{gy:a.gy,gx:b.gx},PostId{gy:b.gy,gx:a.gx})) {continue;}
                let distance=domain.pin.spacing_mm as f64/1000.0 * if diagonal {2f64.sqrt()} else {1.0};
                let weight=((h[i]-h[j])/distance).powf(1.1)/distance;
                if !weight.is_finite() {return Err(invalid("nonfinite MFD slope weight"));}
                edges[i].push(DetailEdge {receiver:j as u32,weight});
                if domain.policy==Policy::NativeDetailV1 && diagonal {diagonals.insert(pair(a,b));}
            }}
        }
        if edges[i].is_empty() {edges[i].push(DetailEdge {receiver:p as u32,weight:1.0});}
        edges[i].sort_by_key(|e|e.receiver);
        let sum:f64=edges[i].iter().map(|e|e.weight).sum();
        if !sum.is_finite() || sum<=0.0 {return Err(invalid("invalid MFD normalization"));}
        for e in &mut edges[i] {e.weight/=sum;}
    }
    let mut a=vec![0.0;count];let mut q=vec![0.0;count];
    let mut balances:BTreeMap<u32,DetailBalance>=BTreeMap::new();
    let mut capacity=0.0;
    for (i,p) in posts.iter().enumerate() {if p.boundary.is_none() {
        a[i]=p.support_m2;q[i]=p.source_m2+p.external_m2;
        capacity+=(h[i]-p.ground_m)*p.support_m2;
        let b=balances.entry(p.basin).or_insert(DetailBalance {basin:p.basin,geometric_source_m2:0.0,effective_source_m2:0.0,geometric_export_m2:0.0,effective_export_m2:0.0});
        b.geometric_source_m2+=a[i];b.effective_source_m2+=q[i];
    }}
    for &ii in order.iter().rev() {
        let i=ii as usize;
        for e in &edges[i] {let j=e.receiver as usize;a[j]+=a[i]*e.weight;q[j]+=q[i]*e.weight;}
    }
    for (i,p) in posts.iter().enumerate() {if let Some(t)=&p.boundary {
        let b=balances.entry(t.basin).or_insert(DetailBalance {basin:t.basin,geometric_source_m2:0.0,effective_source_m2:0.0,geometric_export_m2:0.0,effective_export_m2:0.0});
        b.geometric_export_m2+=a[i];b.effective_export_m2+=q[i];
    }}
    let mut transfers:BTreeMap<(u32,u32),(f64,f64)>=BTreeMap::new();
    if domain.policy==Policy::NativePhysicalV1 {
        for (i,es) in edges.iter().enumerate() { for e in es {
            let j=e.receiver as usize;
            if posts[i].basin!=posts[j].basin {
                let v=transfers.entry((posts[i].basin,posts[j].basin)).or_default();
                v.0+=a[i]*e.weight;v.1+=q[i]*e.weight;
            }
        }}
        for (&basin,b) in &balances {
            let (mut gi,mut qi,mut go,mut qo)=(0.0,0.0,0.0,0.0);
            for (&(from,to),&(g,q)) in &transfers {
                if to==basin {gi+=g;qi+=q;} if from==basin {go+=g;qo+=q;}
            }
            for (input,output) in [(b.geometric_source_m2+gi,b.geometric_export_m2+go),(b.effective_source_m2+qi,b.effective_export_m2+qo)] {
                if !input.is_finite()||!output.is_finite()||(input-output).abs()>1e-9*input.max(1.0) {return Err(invalid("physical basin source/transfer/export balance"));}
            }
        }
    }
    let basin_transfers=transfers.into_iter().map(|((from_basin,to_basin),(geometric_m2,effective_m2))|DetailBasinTransfer{from_basin,to_basin,geometric_m2,effective_m2}).collect();
    // Construction fill retains the unconfined legacy minimax policy. Per-basin
    // routing is only certified for the separate native detail policy.
    if !capacity.is_finite() || a.iter().chain(q.iter()).any(|v|!v.is_finite()) {return Err(invalid("nonfinite capacity or flux"));}
    let balances:Vec<_>=balances.into_values().collect();
    let conserved=|input:f64,output:f64| (input-output).abs()<=1e-9*input.max(1.0);
    if domain.policy==Policy::NativeDetailV1 && balances.iter().any(|b|!conserved(b.geometric_source_m2,b.geometric_export_m2)||!conserved(b.effective_source_m2,b.effective_export_m2)) {return Err(invalid("per-basin detail flux residual exceeds 1e-9"));}
    if !conserved(balances.iter().map(|b|b.geometric_source_m2).sum(),balances.iter().map(|b|b.geometric_export_m2).sum()) || !conserved(balances.iter().map(|b|b.effective_source_m2).sum(),balances.iter().map(|b|b.effective_export_m2).sum()) {return Err(invalid("detail total flux residual exceeds 1e-9"));}
    let kernel_digest=blake3::hash(concat!(include_str!("t1b_ls.rs"),include_str!("t1b_surface_domain.rs")).as_bytes()).to_hex().to_string();
    let mut route=DetailRoute {domain_digest:domain.digest.clone(),kernel_digest,mode,routing_m:h,parent,rank,downstream_first:order,edges,geometric_area_m2:a,effective_area_m2:q,fill_capacity_m3:capacity,balances,basin_transfers,digest:String::new()};
    route.digest=content_digest(&route)?;
    Ok(route)
}

/// Selective construction closure, owned here with the routing/fill policy.
/// Frontier seeds only establish an optimistic bound, never a published outlet.
pub fn condition_detail_bounds(patch:&crate::t1b_surface_completion::Patch<'_>) -> crate::t1b_surface_domain::Result<crate::t1b_surface_completion::Bounds> {
    crate::t1b_surface_completion::solve_bounds(patch)
}

/// Native reporting-channel registration. This never changes physical geometry,
/// fixed boundaries, source flux, or the routing owner graph.
pub fn register_detail_channels(domain:&crate::t1b_surface_domain::Domain,route:&DetailRoute,runoff_reference_m_s:f64,threshold_m3_s:f64)->crate::t1b_surface_domain::Result<crate::t1b_surface_channels::Registry>{
    crate::t1b_surface_channels::register(domain,route,runoff_reference_m_s,threshold_m3_s)
}
