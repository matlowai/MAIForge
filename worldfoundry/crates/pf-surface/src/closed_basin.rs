//! Lane INLAND-BASINS (2026-10-21): law `closed-basin-water-balance-v1`
//! (law family `closed_basin`), an opt-in extension of
//! `continental-retention-base-level-v1` (erosion/retain.rs). Nothing here
//! runs unless the law is selected.
//!
//! DEFECT IT ADDRESSES (lane PEAK-PERSIST REPORT.md, figs/catmap_check.png):
//! inland basins with floors below sea level cover about 1 % of the globe in
//! every C6/C7 run (Earth at the same grid: about 0.15 %, the Caspian most of
//! it; ETOPO 2022 block means, lane INLAND-BASINS earth_inland.py). The
//! routing already treats every closed depression as a lake full to its
//! spill level, so water never decides anything, and the lake sediment of
//! the retention law is laid as deltas from the inlet up to the spill level,
//! so most of it raises lake margins that already stand above sea level
//! while the deep below-sea floor stays starved.
//!
//! THE LAW (constants in `ClosedBasinParams::v1`):
//! 1. WATER BALANCE OF EVERY CLOSED LAKE (closed-basin hydrology; the
//!    terminal-lake balance of Mason et al. 1994 and Garcia-Castellanos 2006,
//!    recalled, not checked). For each routing lake body, upstream lakes
//!    first: inflow I = sum of runoff x area over the dry cells whose flow
//!    first reaches it, plus the overflow of upstream lakes. With the lowest
//!    k body cells submerged, g(k) = I + sum over unsubmerged body cells of
//!    R a - sum over submerged cells of E_net a. If g(all) >= 0 the lake is
//!    FULL at its spill level and passes g(all) downstream; otherwise it is
//!    TERMINAL with the largest k for which g(k) >= 0, its water level is the
//!    elevation of the next body cell (k = 0: a dry playa floor), and nothing
//!    leaves it. Runoff R is the climate law's runoff when the climate
//!    forcing law is on, else `surface.uniform_runoff_m_per_yr`. Net
//!    open-water evaporation E_net = PET - P (climate law fields) or, without
//!    climate, `uniform_evaporation_to_runoff_ratio` x the uniform runoff.
//! 2. FOCUSED LAKE FILL. Of the sediment that reaches a lake, the bedload
//!    share (1 - `suspended_load_fraction`) builds the inlet delta exactly as
//!    the retention law does (breadth-first from the inlet to the spill
//!    level); the suspended share settles on the lake floor and is focused
//!    into the deepest water (sediment focusing; Blais and Kalff 1995,
//!    recalled): it fills the body's lowest cells to one common compensated
//!    level. What a full delta could not place joins the floor fill instead
//!    of leaving the lake. Only when the whole body stands at its spill level
//!    does the rest pass downstream (a full lake overflows; a terminal lake
//!    filled to its sill would overflow too). Lakes are processed upstream
//!    first, so an overflow reaches the next lake in the same step.
//! 3. Water levels are published per cell (`cell_water_level_m`) so a cell
//!    below sea level is DRY land only where evaporation wins (a terminal
//!    lake below its floor); under a full lake it is a lake bed.
//! 4. FILLED BASINS (lane LAKE-EXCESS, 2026-10-21; opt-in law
//!    `filled-basin-alluvial-grade-v1`, `ClosedBasinParams::
//!    filled_basin_grading_depth_m`). DEFECT: lake plains filled flat to
//!    just under the spill level stay lakes (cells with < 1 m of water were
//!    ~2-3 % of continental area, ~20 % of all lake area, in C7 runs at
//!    120-400 My; lane LAKE-EXCESS metrics/lakes_C7.json) and any small
//!    subsidence of a dead-flat fill re-ponds it. A real filled basin is an
//!    alluvial plain that slopes to its outlet. Rule: when the step's floor
//!    fill can bring the lake to within d_c of its spill level, the lake is
//!    a filled basin; its pool instead builds the graded surface
//!    z_g = spill + S_a D (D: flow distance to the outlet through the body
//!    along the routing's flow tree; S_a the retention law's graded alluvial
//!    slope, so the plain continues the land grade rule 2a), filling the
//!    cells lowest relative to z_g first (the surface z_g - delta with one
//!    common delta >= 0). What z_g cannot hold passes down the outlet as a
//!    full lake's overflow does. The next step's routing then sees a plain
//!    draining across the sill, not a depression. No new slope constant.
//! 5. OUTLET INCISION (lane LAKE-OUTLET, 2026-10-21; opt-in law
//!    `lake-outlet-sill-incision-v1`, `ClosedBasinParams::outlet`). DEFECT:
//!    humid closed depressions stay full lakes for hundreds of My (lakes
//!    9-12 % of continental area vs Earth ~1-2 %; 3-4 lakes > 1e6 km2; lanes
//!    LAKE-EXCESS and CLIMATE-ARIDITY). On Earth an overflowing lake drains
//!    because its outlet river cuts down through the sill, lowering the
//!    spill level, until the lake is gone or evaporation balances inflow
//!    (outlet erosion and overfill end every post-tectonic lake; lake
//!    evaporation prolongs a lake by starving the outlet: Garcia-Castellanos
//!    2006, GSA Special Paper 398, abstract read; outlet erodibility from
//!    overspill floods: Garcia-Castellanos and O'Connor 2018, Nature
//!    Communications 9, abstract read). In the model the spill level is
//!    the cell-mean elevation of the lake's sill (pour) cell, and that mean is
//!    lowered only by the subgrid hillslope closure times the retention
//!    taper (both slow on a low-relief rim near base level), while the
//!    real spill level is the floor of the outlet channel. Rule: for each
//!    FULL lake of the balance, its sill cell (followed from its first exit
//!    edge down the flow tree: the first dry cell at the spill level whose
//!    receiver stands lower, `sill_discharge_area_m2`) is lowered in the retention step
//!    by the CHANNEL stream-power law on the lake's overflow discharge:
//!    the trunk equation of the selected relief law (E = K_n (A_q^m S^n -
//!    k_c^n), the incision threshold when selected; the retention step's
//!    own linear law when no relief law is selected), backward Euler to the
//!    sill's receiver as every fluvial update, with no taper and no
//!    subgrid closure; A_q = overflow (m^3/yr) / reference runoff (the
//!    climate law's r_ref, else the uniform runoff), so a lake with no net
//!    evaporation passes its whole catchment's area and a terminal lake
//!    passes none. The sill takes the larger of this and the retention
//!    step's own lowering. The lowered rock is detached material like any
//!    other (it travels down the flow path), and the next step's routing
//!    sees the lower sill, so the lake level follows. No new constant.
//!    DECLARED APPROXIMATION: the cell mean of the sill cell is lowered at
//!    the channel's rate, so the rock removed at a sill is a cell's area
//!    times the incision, not a gorge's (over-removal, local to sills).
//! 6. SILL FLEXURAL SUPPORT (lane SILL-REBOUND, 2026-10-21; opt-in law
//!    `lake-outlet-sill-incision-flexural-support-v1`,
//!    `OutletIncisionParams::flexural_support`). DEFECT: rule 5 lowers the
//!    sill by removing that thickness from the column, and the next local-
//!    Airy solve returns (rho_m - rho_c)/rho_m ~ 85 % of it, so the spill
//!    level drops ~15 % of the cut (measured: drop / removal median 0.159
//!    over 63,843 sill evaluations, seed 42 C7 climate, lane SILL-REBOUND
//!    metrics/sill_diag_s42.json). On Earth an outlet channel is km wide,
//!    far narrower than the flexural wavelength (2 pi alpha ~ 400-650 km
//!    for Te 40 km), so the plate carries the cut and its floor drops by the
//!    full incision; the compensation of the small removed load is regional
//!    and negligible (Turcotte and Schubert 2014, ch. 3; Watts 2001;
//!    recalled, not checked). Rule: at a sill cell the lowering is the
//!    post-isostasy drop: the material removed is the compensated column's
//!    thickness for that drop (`CompensatedColumn::thickness_to_raise`,
//!    the retention law's densities, piecewise across sea level), as the
//!    critical-slope rule sizes its removal for the post-isostasy surface;
//!    the step's working surface drops by the lowering only. No new
//!    constant. DECLARED APPROXIMATION: the cell mean stands for the
//!    channel floor (rule 5), so the rock removed at a sill is now a whole
//!    cell's compensated column (~6.6x the drop on land) instead of a
//!    gorge's, exported downstream as any detached material; the regional
//!    flexural uplift of the removed load is not modelled.
//! 7. OUTLET PATH INCISION (lane LAKE-HOLD, 2026-10-21; opt-in laws
//!    `lake-outlet-path-incision-v1` and, with rule 6,
//!    `lake-outlet-path-incision-flexural-support-v1`;
//!    `OutletIncisionParams::path_incision`). DEFECT: with rules 5 and 6
//!    the sill drops to its receiver every step, yet lakes > 1e6 km2 stay
//!    (median 152 m above the next water body, the ocean for 96 % of their
//!    area-steps, along ~8 dry path cells; seed 42 C7 climate with rule 6,
//!    lane LAKE-HOLD metrics/hold_diag_s42.json). The dry cells between
//!    the sill and the ocean are never lowered: their lowering is the
//!    subgrid cell-mean closure times the retention taper, and the
//!    incision threshold sets that closure to 0 on 92 % of them (the
//!    cell-mean relief above the trunk is below Lambda k_c), so the whole
//!    corridor holds the lake while only its first cell is cut. On Earth
//!    the outlet river carries the lake's overflow all the way to the next
//!    water body and cuts its channel along the whole course, and the lake
//!    level follows the channel floor, not the floor of the sill alone.
//!    Rule: every dry cell of the overflow path, from the sill down the flow
//!    tree to the next water body (ocean, another lake, or a terminal
//!    floor), takes rule 5's lowering with the same A_q (summed where
//!    paths of several lakes merge) and, with rule 6, rule 6's
//!    compensation. No new constant. DECLARED APPROXIMATION (as rule 5):
//!    the cell mean stands for the channel floor along the path, so the
//!    rock removed is a cell's area times the incision, exported downstream
//!    as any detached material; rule 2a's graded-plain storage still acts
//!    on these cells (the grade S_a is unchanged).
//! 8. OUTLET CUT THROUGH LAKES (lane LAKE-HOLD-2, 2026-10-21; opt-in laws
//!    `lake-outlet-path-through-lakes-v1` and, with rule 6,
//!    `lake-outlet-path-through-lakes-flexural-support-v1`;
//!    `OutletIncisionParams::through_lakes`, requires rule 7). DEFECT: with
//!    rule 7 a river entering another full lake on the overflow path is cut
//!    only down to that lake's START-of-step level, even when the same step
//!    cuts that lake's own sill (receivers first, so the downstream sill is
//!    already lowered), so a cut travels down a chain of small lakes about
//!    one lake per step (lane LAKE-HOLD, outlet profile of seed 42). On
//!    Earth a lake whose outlet is cut falls with its outlet, and every river
//!    entering it sees the lower level at once. Rule: the base level of a
//!    cell whose receiver is a FULL lake with a sill (rule 5's sill cell,
//!    `outlet_sill_cells`) is that lake's post-incision spill: the sill's
//!    surface after this step's stage 1, never above the start-of-step
//!    level and never below the receiver cell's own bed (a drained receiver
//!    is dry ground). Used for the stage-1 base (stream power and the
//!    outlet law) and for rule 2a's receiver level. No new constant.
//!    DECLARED APPROXIMATIONS: one level per lake (no drawdown gradient);
//!    every river entering a lowered lake sees it, not only the overflow
//!    path; the taper's base level, lake deposition, hillslope diffusion
//!    and the subgrid closure's target still use the start-of-step level
//!    (the closure's effective area is realised against the lower base, so
//!    such a cell may lose slightly more than its closure target).
//!    MEASURED (lane LAKE-HOLD-2, C7 climate + rules 6 and 7, seeds 42,
//!    240907, 7, 40-200 My, plates identical): lake area 5.19/4.03/3.24 ->
//!    5.33/4.00/3.00 % of continents, per-step lake water -0.4/-2.6/-4.1 %;
//!    erosion +0.7 to +1.2 %. The cut does pass through small lakes within
//!    the step, but the lag it removes is a few My against lake lifetimes
//!    of tens to hundreds of My: no measurable lever. Opt-in only.
//! 9. OUTLET PATH GRADE (lane LAKE-HOLD-2, 2026-10-21; opt-in law
//!    `lake-outlet-path-threshold-grade-v1`, its own law-set slot;
//!    `OutletIncisionParams::path_grade`, requires rule 7 and the incision
//!    threshold). DEFECT: rule 2a re-aggrades every incised path cell to
//!    receiver + S_a L with the one alluvial slope S_a = 1e-4 of every
//!    river, so a lake's big outlet river keeps rebuilding a plain the
//!    channel law then cuts again (30 % of the hold left after rule 7, lane
//!    LAKE-HOLD metrics/hold_attr_s42.txt). On Earth larger rivers grade to
//!    gentler slopes (lowland trunk rivers 1e-5 to 1e-4; bankfull slope
//!    falls with discharge; recalled, not checked). Rule: on cells with
//!    outlet discharge area A_q > 0 (sill and path) rule 2a grades to
//!    min(S_a, k_c / A_q^(m/n)), the outlet channel law's own threshold slope
//!    (k_c = (K_n k_c^n / K_n)^(1/n) from the selected laws): a plain
//!    steeper than that would be cut by the same law, so storage stops at
//!    the slope the overflow river can hold. No new constant. DECLARED
//!    APPROXIMATION: the bedrock-detachment threshold stands in for the
//!    alluvial threshold of motion (which is lower, so the true grade can
//!    be gentler still); only overflow-path cells get it.
//!    MEASURED (lane LAKE-HOLD-2, as rule 8): lake area 5.19/4.03/3.24 ->
//!    4.90/3.89/2.77 % of continents (-6/-3/-15 %), per-step lake water
//!    -0.4/-1.0/-1.2 %; path grade storage -23 to -32 %; erosion and
//!    hypsometry unchanged. Small, no cost. Opt-in only.
//! 10. ALLUVIAL OUTLET CHANNEL (lane ALLUVIAL-TRUNK, 2026-10-21; opt-in law
//!    `stream-power-threshold-free-alluvial-channels-v1`, which also acts on
//!    the subgrid trunk, pf-surface erosion/incision_threshold.rs;
//!    `OutletIncisionParams::alluvial_grade_slope`, requires a channel law
//!    with an incision threshold). DEFECT: after rules 5-7 the outlet
//!    channel law still holds every path cell at its bedrock threshold slope
//!    k_c / A_q^(m/n), and on 750-6,000 km lowland overflow paths that holds
//!    giant lakes 100-250 m up (lane LAKE-HOLD-2, outlet profile of seed 42
//!    at 134 My: ~5e-5 over 5,400-6,000 km). Rule: on cells whose outlet
//!    discharge area makes the channel alluvial (k_c / A_q^(m/n) < S_a, S_a
//!    the retention law's graded alluvial slope; the switch and its
//!    justification are in incision_threshold.rs) the outlet law runs with
//!    no threshold term (K_n k_c^n dt dropped); rule 2a then grades the
//!    reach at S_a when it carries load (with rule 9, at min(S_a, k_c /
//!    A_q^(m/n)) as before). No new constant.
//!    MEASURED (lane ALLUVIAL-TRUNK, as incision_threshold.rs): no lake
//!    lever. Lake area 5.19/4.03/3.24 -> 5.31/4.80/3.24 % of continents
//!    (with rule 9: 5.03/3.87/3.14), per-step lake water +1.3/-2.0/+1.2 %.
//!    The released path is cut deeper, rule 2a re-aggrades it at S_a (graded
//!    storage +9 to +20 %), and the seed-42 134 My corridor ends up
//!    rougher (mean end-of-step slope 1.8e-5 -> 4.1e-5), not lower. Opt-in
//!    only.
//! MASS: every raise goes through the retention law's `deposit`, so eroded =
//! deposited + exported is unchanged.
//! NOT MODELLED (declared): the water level does not feed back into routing,
//! base level or isostasy (erosion still sees every closed depression filled
//! to its spill level; Airy water loading still uses the sea datum, the
//! isostasy law's global fill); no salinity, no lake-level history (each step
//! is a steady balance), no evaporite deposition, no groundwater.

use crate::routing::Routing;

/// Law id (the typed law set lives in pf-contracts laws.rs).
pub const CLOSED_BASIN_LAW_ID: &str = "closed-basin-water-balance-v1";
/// Suspended share of river load reaching a lake. Project assumption
/// awaiting review: bedload is commonly taken as about 10 % of the total
/// river load (Milliman and Meade 1983; Turowski et al. 2010; recalled, not
/// checked).
pub const SUSPENDED_LOAD_FRACTION: f64 = 0.9;
/// Net open-water evaporation over runoff without a climate law. Project
/// assumption awaiting review: Earth land means, lake evaporation about
/// 1.2 m/yr minus land precipitation about 0.8 m/yr, over land runoff about
/// 0.27 m/yr (Oki and Kanae 2006; Trenberth et al. 2007; recalled, not
/// checked) -> 0.4 / 0.27 ~ 1.5. A lake stays full when its land catchment
/// is at least this many times its own area.
pub const UNIFORM_EVAPORATION_TO_RUNOFF_RATIO: f64 = 1.5;

/// Lane LAKE-EXCESS (2026-10-21): law `filled-basin-alluvial-grade-v1`, an
/// opt-in extension of this law (variant of the same law-set slot).
pub const FILLED_BASIN_GRADE_LAW_ID: &str = "filled-basin-alluvial-grade-v1";
/// d_c: residual water depth (m) over a lake's sediment fill at and below
/// which the lake is a filled basin. Project assumption awaiting review:
/// once the water over the fill is shallower than the inflowing river's
/// own flow (bankfull) depth, the river no longer enters standing water and
/// drops its suspended load there; it crosses the fill as a channel and
/// builds an alluvial plain graded to the outlet. Bankfull depths of large
/// lowland rivers are roughly 5-20 m (Leopold and Maddock 1953 hydraulic
/// geometry; lower Mississippi ~15-20 m, Rhine ~5-10 m; recalled, not
/// checked); 10 m taken as a mid value.
pub const FILLED_BASIN_GRADING_DEPTH_M: f64 = 10.0;

/// Lane LAKE-OUTLET (2026-10-21): law `lake-outlet-sill-incision-v1`, an
/// opt-in extension of this law (its own law-set slot).
pub const LAKE_OUTLET_LAW_ID: &str = "lake-outlet-sill-incision-v1";
/// Lane SILL-REBOUND (2026-10-21): law
/// `lake-outlet-sill-incision-flexural-support-v1` (rule 6), a variant of the
/// outlet law's law-set slot.
pub const LAKE_OUTLET_FLEXURAL_SUPPORT_LAW_ID: &str = "lake-outlet-sill-incision-flexural-support-v1";
/// Lane LAKE-HOLD (2026-10-21): laws `lake-outlet-path-incision-v1` and
/// `lake-outlet-path-incision-flexural-support-v1` (rule 7), variants of
/// the outlet law's law-set slot.
pub const LAKE_OUTLET_PATH_LAW_ID: &str = "lake-outlet-path-incision-v1";
pub const LAKE_OUTLET_PATH_FLEXURAL_SUPPORT_LAW_ID: &str = "lake-outlet-path-incision-flexural-support-v1";
/// Lane LAKE-HOLD-2 (2026-10-21): laws `lake-outlet-path-through-lakes-v1`
/// and `lake-outlet-path-through-lakes-flexural-support-v1` (rule 8),
/// variants of the outlet law's law-set slot.
pub const LAKE_OUTLET_THROUGH_LAKES_LAW_ID: &str = "lake-outlet-path-through-lakes-v1";
pub const LAKE_OUTLET_THROUGH_LAKES_FLEXURAL_SUPPORT_LAW_ID: &str = "lake-outlet-path-through-lakes-flexural-support-v1";
/// Lane LAKE-HOLD-2 (2026-10-21): law `lake-outlet-path-threshold-grade-v1`
/// (rule 9), its own law-set slot.
pub const LAKE_OUTLET_PATH_GRADE_LAW_ID: &str = "lake-outlet-path-threshold-grade-v1";

/// The channel stream-power law a lake's outlet incises its sill with
/// (rule 5). Every value comes from laws already selected; none is new.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutletIncisionParams {
    /// Runoff (m/yr) at which the fluvial law's K is quoted: the climate
    /// law's r_ref, else the uniform runoff. A_q = overflow / this.
    pub reference_runoff_m_per_yr: f64,
    /// The relief law's trunk law (K_n, m, n, K_n k_c^n); None = the
    /// retention step's own linear stream power (n = 1, no threshold).
    pub channel: Option<OutletChannelLaw>,
    /// SILL-REBOUND (`lake-outlet-sill-incision-flexural-support-v1`, rule
    /// 6): no local isostatic rebound at a sill; false = rule 5 exactly.
    pub flexural_support: bool,
    /// LAKE-HOLD (`lake-outlet-path-incision-v1`, rule 7): the channel law
    /// acts on every dry cell of the overflow path from the sill to the next
    /// water body; false = the sill alone (rules 5 and 6 exactly).
    pub path_incision: bool,
    /// LAKE-HOLD-2 (`lake-outlet-path-through-lakes-v1`, rule 8; requires
    /// `path_incision`): a river entering a full lake whose sill this step
    /// cuts is graded to the lake's post-incision spill; false = rules 5-7
    /// exactly.
    pub through_lakes: bool,
    /// LAKE-HOLD-2 (`lake-outlet-path-threshold-grade-v1`, rule 9; requires
    /// `path_incision` and a channel law with a threshold): rule 2a's grade
    /// on cells with A_q > 0 is min(S_a, the channel's threshold slope);
    /// false = rules 5-8 exactly.
    pub path_grade: bool,
    /// ALLUVIAL-TRUNK (`stream-power-threshold-free-alluvial-channels-v1`,
    /// rule 10; requires a channel law with an incision threshold): S_a; the
    /// outlet law's threshold is dropped on cells where k_c / A_q^(m/n) <
    /// S_a. None = rules 5-9 exactly.
    pub alluvial_grade_slope: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutletChannelLaw {
    /// K_n (m^(1-n) per My).
    pub k_n_per_my: f64,
    /// Area exponent m (= n / 2 in the relief law).
    pub m: f64,
    /// Slope exponent n.
    pub n: f64,
    /// K_n k_c^n (m/My); 0 without the incision-threshold law.
    pub threshold_rate_m_per_my: f64,
}

impl OutletIncisionParams {
    pub fn validate(&self) -> Result<(), String> {
        let r = self.reference_runoff_m_per_yr;
        if !(r.is_finite() && r > 0.0) {
            return Err(format!("{LAKE_OUTLET_LAW_ID}: reference runoff {r} must be finite and > 0"));
        }
        if let Some(c) = &self.channel {
            for (what, v) in [("K_n", c.k_n_per_my), ("m", c.m), ("n", c.n)] {
                if !(v.is_finite() && v > 0.0) {
                    return Err(format!("{LAKE_OUTLET_LAW_ID}: channel {what} {v} must be finite and > 0"));
                }
            }
            if !(c.threshold_rate_m_per_my.is_finite() && c.threshold_rate_m_per_my >= 0.0) {
                return Err(format!("{LAKE_OUTLET_LAW_ID}: threshold rate {} must be finite and >= 0", c.threshold_rate_m_per_my));
            }
        }
        if self.through_lakes && !self.path_incision {
            return Err(format!("{LAKE_OUTLET_THROUGH_LAKES_LAW_ID} requires {LAKE_OUTLET_PATH_LAW_ID} (path_incision)"));
        }
        if self.path_grade && !(self.path_incision && self.channel.is_some_and(|c| c.threshold_rate_m_per_my > 0.0)) {
            return Err(format!("{LAKE_OUTLET_PATH_GRADE_LAW_ID} requires {LAKE_OUTLET_PATH_LAW_ID} and a channel law with an incision threshold"));
        }
        if let Some(s) = self.alluvial_grade_slope {
            let law = crate::erosion::incision_threshold::ALLUVIAL_CHANNEL_LAW_ID;
            if !(s.is_finite() && s > 0.0) {
                return Err(format!("{law}: alluvial grade slope {s} must be finite and > 0"));
            }
            if !self.channel.is_some_and(|c| c.threshold_rate_m_per_my > 0.0) {
                return Err(format!("{law} requires a channel law with an incision threshold"));
            }
        }
        Ok(())
    }

    /// Rule 10 (ALLUVIAL-TRUNK): true when the outlet channel with discharge
    /// area `a_q_m2` is alluvial (k_c / A_q^(m/n) < S_a) and the law is on.
    pub fn alluvial_channel(&self, a_q_m2: f64) -> bool {
        match (&self.channel, self.alluvial_grade_slope) {
            (Some(c), Some(s_a)) if c.threshold_rate_m_per_my > 0.0 => {
                let k_c = (c.threshold_rate_m_per_my / c.k_n_per_my).powf(1.0 / c.n);
                crate::erosion::incision_threshold::is_alluvial_channel(k_c, s_a, a_q_m2, c.m / c.n)
            }
            _ => false,
        }
    }

    /// Rule 9 (LAKE-HOLD-2): rule 2a's grade on a cell with outlet
    /// discharge area `a_q_m2` > 0: min(S_a, k_c / A_q^(m/n)), k_c from the
    /// channel law (K_n k_c^n / K_n)^(1/n). `s_a` itself without `path_grade`,
    /// a channel or a positive A_q.
    pub fn path_grade_slope(&self, s_a: f64, a_q_m2: f64) -> f64 {
        match &self.channel {
            Some(c) if self.path_grade && a_q_m2 > 0.0 && c.threshold_rate_m_per_my > 0.0 => {
                let k_c = (c.threshold_rate_m_per_my / c.k_n_per_my).powf(1.0 / c.n);
                s_a.min(k_c / a_q_m2.powf(c.m / c.n))
            }
            _ => s_a,
        }
    }

    /// Lowering (m, >= 0) of a sill standing `x0_m` above its receiver's
    /// level, over `dt_my`, with outlet discharge area `area_m2` and reach
    /// `len_m`; `linear_k_per_my` and `linear_m` are the retention step's
    /// own law (used when `channel` is None). Backward Euler with the
    /// threshold held explicit, as the trunk update of the threshold law:
    /// x' + F x'^n = x0 + K_n k_c^n dt, lowering = max(0, x0 - x').
    pub fn sill_lowering_m(&self, x0_m: f64, area_m2: f64, len_m: f64, dt_my: f64, linear_k_per_my: f64, linear_m: f64) -> f64 {
        if !(x0_m > 0.0) || !(area_m2 > 0.0) {
            return 0.0;
        }
        let (k, m, n, g) = match &self.channel {
            // Rule 10 (ALLUVIAL-TRUNK): no threshold on an alluvial channel.
            Some(c) if self.alluvial_channel(area_m2) => (c.k_n_per_my, c.m, c.n, 0.0),
            Some(c) => (c.k_n_per_my, c.m, c.n, c.threshold_rate_m_per_my * dt_my),
            None => (linear_k_per_my, linear_m, 1.0, 0.0),
        };
        let f = k * area_m2.powf(m) * dt_my / len_m.powf(n);
        let x = crate::relief::implicit_power_root(x0_m + g, f, n).min(x0_m);
        (x0_m - x).max(0.0)
    }
}

/// Rule 5: the outlet discharge area A_q (m^2) at each lake's sill cell
/// (0 elsewhere): overflow / reference runoff, summed when several full
/// lakes spill over one cell. Only FULL lakes count. The sill is found by
/// following the overflow from the lake's first exit edge (the edge the
/// balance sends the overflow along) down the flow tree, past dry cells
/// standing at the spill level (a fill or flat at exactly the spill level
/// can sit on the shore; it is not the sill, its receiver is level with
/// it) and past cells of the lake itself, to the first dry cell whose
/// receiver stands lower on the filled surface: the cell that holds the
/// lake up. A path that reaches another water body or a terminal floor
/// first has no sill to incise.
pub fn sill_discharge_area_m2(routing: &Routing, lb: &LakeBalance, reference_runoff_m_per_yr: f64) -> Vec<f64> {
    outlet_discharge_area_m2(routing, lb, reference_runoff_m_per_yr, false)
}

/// Rules 5 and 7: `sill_discharge_area_m2` when `path` is false; with
/// `path` (rule 7, LAKE-HOLD) A_q also reaches every dry cell below the
/// sill along the flow tree until the overflow meets the ocean, another
/// water body or a cell with no receiver (a terminal floor, which the
/// retention step never lowers).
pub fn outlet_discharge_area_m2(routing: &Routing, lb: &LakeBalance, reference_runoff_m_per_yr: f64, path: bool) -> Vec<f64> {
    let n = routing.water_body_id.len();
    let mut a = vec![0.0f64; n];
    let sills = outlet_sill_cells(routing, lb);
    for b in 0..routing.water_bodies.len() {
        let Some(s) = sills[b] else {
            continue;
        };
        let q = lb.overflow_m3_per_yr.get(b).copied().unwrap_or(0.0);
        a[s] += q / reference_runoff_m_per_yr;
        if path {
            // Rule 7: on down the flow tree through dry cells.
            let mut c = s;
            let mut hops = 0usize;
            while let Some(j) = routing.flow_to[c] {
                c = j.0 as usize;
                if routing.ocean[c] || routing.water_body_id[c] != crate::water::DRY || routing.flow_to[c].is_none() {
                    break;
                }
                a[c] += q / reference_runoff_m_per_yr;
                hops += 1;
                if hops > n {
                    break;
                }
            }
        }
    }
    a
}

/// Rule 5's sill cell of every water body (index = water-body id): `Some`
/// only for a FULL, non-marine lake with overflow > 0 whose overflow path
/// reaches a dry cell holding it up (see `sill_discharge_area_m2`); `None`
/// otherwise. Rules 5-7 use it for A_q; rule 8 (LAKE-HOLD-2) for the
/// lake's post-incision spill.
pub fn outlet_sill_cells(routing: &Routing, lb: &LakeBalance) -> Vec<Option<usize>> {
    let n = routing.water_body_id.len();
    let fill = &routing.fill_level_m;
    let mut out = vec![None; routing.water_bodies.len()];
    for b in 0..routing.water_bodies.len() {
        if is_marine(routing, b) || lb.terminal.get(b).copied().unwrap_or(true) {
            continue;
        }
        let q = lb.overflow_m3_per_yr.get(b).copied().unwrap_or(0.0);
        if !(q > 0.0) {
            continue;
        }
        let Some([_, to]) = routing.water_bodies[b].outlet_edges.first() else {
            continue;
        };
        let spill = routing.water_bodies[b].water_surface_m;
        let mut c = to.0 as usize;
        let mut hops = 0usize;
        out[b] = loop {
            if routing.ocean[c] {
                break None;
            }
            let wb = routing.water_body_id[c];
            if wb != crate::water::DRY && wb as usize != b {
                break None; // another water body at the spill level
            }
            let Some(j) = routing.flow_to[c] else {
                break None;
            };
            let j = j.0 as usize;
            if wb == crate::water::DRY && fill[j] < fill[c] && fill[c] >= spill {
                break Some(c);
            }
            c = j;
            hops += 1;
            if hops > n {
                break None;
            }
        };
    }
    out
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosedBasinParams {
    pub suspended_load_fraction: f64,
    pub uniform_evaporation_to_runoff_ratio: f64,
    /// `filled-basin-alluvial-grade-v1`: d_c (m); None = the v1 law exactly.
    pub filled_basin_grading_depth_m: Option<f64>,
    /// `lake-outlet-sill-incision-v1` (rule 5); None = no outlet incision.
    pub outlet: Option<OutletIncisionParams>,
}

impl ClosedBasinParams {
    pub fn v1() -> Self {
        Self {
            suspended_load_fraction: SUSPENDED_LOAD_FRACTION,
            uniform_evaporation_to_runoff_ratio: UNIFORM_EVAPORATION_TO_RUNOFF_RATIO,
            filled_basin_grading_depth_m: None,
            outlet: None,
        }
    }

    /// v1 plus `filled-basin-alluvial-grade-v1` at the declared d_c.
    pub fn v1_graded() -> Self {
        Self {
            filled_basin_grading_depth_m: Some(FILLED_BASIN_GRADING_DEPTH_M),
            ..Self::v1()
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        let f = self.suspended_load_fraction;
        if !f.is_finite() || !(0.0..=1.0).contains(&f) {
            return Err(format!("{CLOSED_BASIN_LAW_ID}: suspended_load_fraction {f} outside [0, 1]"));
        }
        let r = self.uniform_evaporation_to_runoff_ratio;
        if !r.is_finite() || r < 0.0 {
            return Err(format!("{CLOSED_BASIN_LAW_ID}: uniform_evaporation_to_runoff_ratio {r} must be finite and >= 0"));
        }
        if let Some(o) = &self.outlet {
            o.validate()?;
        }
        if let Some(d) = self.filled_basin_grading_depth_m {
            if !d.is_finite() || d < 0.0 {
                return Err(format!("{FILLED_BASIN_GRADE_LAW_ID}: filled_basin_grading_depth_m {d} must be finite and >= 0"));
            }
        }
        Ok(())
    }
}

/// Water balance of every routing water body (index = water body id).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LakeBalance {
    /// Water level per body (m): the spill level of a full lake, the
    /// balance level of a terminal lake (its floor when dry); the routing
    /// level for the marine body.
    pub level_m: Vec<f64>,
    /// True for a lake whose evaporation wins below the spill level.
    pub terminal: Vec<bool>,
    /// The next water body an overflow of this body reaches (None: marine
    /// body, terminal, or a terminal floor downstream).
    pub downstream: Vec<Option<usize>>,
    /// Lake bodies, every lake before the lake its overflow reaches.
    pub order: Vec<usize>,
    /// Water inflow per body (m^3/yr: land runoff plus upstream overflow).
    pub inflow_m3_per_yr: Vec<f64>,
    /// LAKE-OUTLET: what a full lake passes downstream (m^3/yr; 0 for a
    /// terminal lake and the marine body).
    pub overflow_m3_per_yr: Vec<f64>,
    pub full_lakes: u64,
    pub terminal_lakes: u64,
    /// Lake area under water (m^2) and lake-body area left dry (m^2).
    pub wet_area_m2: f64,
    pub exposed_area_m2: f64,
}

fn is_marine(r: &Routing, b: usize) -> bool {
    r.water_bodies[b].kind == pf_contracts::domain::WaterBodyKindV0::Marine
}

/// Step 1 of the law. `runoff_m_per_yr` and `net_evaporation_m_per_yr` are
/// per-cell fields (m/yr) on the routing grid.
pub fn lake_balance(
    routing: &Routing,
    elevation_m: &[f64],
    areas_m2: &[f64],
    runoff_m_per_yr: &[f64],
    net_evaporation_m_per_yr: &[f64],
) -> Result<LakeBalance, String> {
    let n = elevation_m.len();
    for (what, len) in [
        ("areas", areas_m2.len()),
        ("runoff", runoff_m_per_yr.len()),
        ("net evaporation", net_evaporation_m_per_yr.len()),
        ("water_body_id", routing.water_body_id.len()),
    ] {
        if len != n {
            return Err(format!("{CLOSED_BASIN_LAW_ID}: {what} has {len} cells, expected {n}"));
        }
    }
    for c in 0..n {
        if !runoff_m_per_yr[c].is_finite() || !net_evaporation_m_per_yr[c].is_finite() {
            return Err(format!("{CLOSED_BASIN_LAW_ID}: nonfinite runoff or evaporation at cell {c}"));
        }
    }
    let nb = routing.water_bodies.len();
    let body = |c: usize| -> Option<usize> {
        let b = routing.water_body_id[c];
        if b == crate::water::DRY {
            None
        } else {
            Some(b as usize)
        }
    };
    // The water body each dry cell's flow first reaches.
    let mut first: Vec<Option<usize>> = vec![None; n];
    for &i in &routing.receiver_first {
        first[i] = match body(i) {
            Some(b) => Some(b),
            None => routing.flow_to[i].and_then(|j| first[j.0 as usize]),
        };
    }
    let mut cells: Vec<Vec<usize>> = vec![Vec::new(); nb];
    let mut land_in = vec![0.0f64; nb];
    for c in 0..n {
        match body(c) {
            Some(b) => cells[b].push(c),
            None => {
                if let Some(b) = first[c] {
                    land_in[b] += runoff_m_per_yr[c].max(0.0) * areas_m2[c];
                }
            }
        }
    }
    // Overflow target of each lake: its first exit edge (sorted by the
    // routing), then the first water body down that path.
    let mut downstream: Vec<Option<usize>> = vec![None; nb];
    for b in 0..nb {
        if is_marine(routing, b) {
            continue;
        }
        if let Some([_, to]) = routing.water_bodies[b].outlet_edges.first() {
            let j = to.0 as usize;
            downstream[b] = match body(j) {
                Some(d) if d != b => Some(d),
                Some(_) => None,
                None => first[j],
            };
            // LAKE-EXCESS fix (2026-10-21): the first exit edge can lead to
            // a dry cell at the lake's own level whose flow re-enters the
            // same lake (a flat dry cell inside or on the shore of the body:
            // water::bodies lists every body cell whose flow leaves the
            // body). That made the lake its own downstream and the overflow
            // order failed with "lake overflow graph has a cycle" (seed
            // 240907, climate law on, t = 294 My). Follow the flow on
            // through the lake, as the sediment overflow walk in
            // erosion/retain.rs does, to the first other water body. Only
            // this case changes, so every run that completed before is
            // unchanged.
            if downstream[b] == Some(b) {
                let mut cur = j;
                let mut hops = 0usize;
                downstream[b] = loop {
                    match body(cur) {
                        Some(d) if d != b => break Some(d),
                        _ => {}
                    }
                    match routing.flow_to[cur] {
                        Some(next) => cur = next.0 as usize,
                        None => break None,
                    }
                    hops += 1;
                    if hops > n {
                        return Err(format!("{CLOSED_BASIN_LAW_ID}: overflow walk of lake {b} did not end"));
                    }
                };
            }
        }
    }
    // Kahn order over lake -> downstream lake edges (a forest of the flow DAG).
    let mut indeg = vec![0usize; nb];
    for b in 0..nb {
        if let Some(d) = downstream[b] {
            if !is_marine(routing, d) && !is_marine(routing, b) {
                indeg[d] += 1;
            }
        }
    }
    let mut order: Vec<usize> = (0..nb).filter(|&b| !is_marine(routing, b) && indeg[b] == 0).collect();
    let mut head = 0;
    while head < order.len() {
        let b = order[head];
        head += 1;
        if let Some(d) = downstream[b] {
            if !is_marine(routing, d) {
                indeg[d] -= 1;
                if indeg[d] == 0 {
                    order.push(d);
                }
            }
        }
    }
    let lakes = (0..nb).filter(|&b| !is_marine(routing, b)).count();
    if order.len() != lakes {
        return Err(format!("{CLOSED_BASIN_LAW_ID}: lake overflow graph has a cycle ({} of {lakes} ordered)", order.len()));
    }
    let mut out = LakeBalance {
        level_m: routing.water_bodies.iter().map(|w| w.water_surface_m).collect(),
        terminal: vec![false; nb],
        downstream: downstream.clone(),
        order: order.clone(),
        inflow_m3_per_yr: vec![0.0; nb],
        overflow_m3_per_yr: vec![0.0; nb],
        ..Default::default()
    };
    let mut upstream_in = vec![0.0f64; nb];
    for &b in &order {
        let inflow = land_in[b] + upstream_in[b];
        out.inflow_m3_per_yr[b] = inflow;
        let mut cs = cells[b].clone();
        cs.sort_by(|&x, &y| elevation_m[x].total_cmp(&elevation_m[y]).then(x.cmp(&y)));
        // g(0): every body cell dry and draining into the floor.
        let mut g = inflow + cs.iter().map(|&c| runoff_m_per_yr[c].max(0.0) * areas_m2[c]).sum::<f64>();
        let mut k = 0usize;
        while k < cs.len() {
            let c = cs[k];
            let next = g - runoff_m_per_yr[c].max(0.0) * areas_m2[c] - net_evaporation_m_per_yr[c] * areas_m2[c];
            if next < 0.0 {
                break;
            }
            g = next;
            k += 1;
        }
        let wet_a: f64 = cs[..k].iter().map(|&c| areas_m2[c]).sum();
        let all_a: f64 = cs.iter().map(|&c| areas_m2[c]).sum();
        out.wet_area_m2 += wet_a;
        out.exposed_area_m2 += all_a - wet_a;
        if k == cs.len() {
            out.full_lakes += 1;
            out.overflow_m3_per_yr[b] = g.max(0.0);
            if let Some(d) = downstream[b] {
                upstream_in[d] += g.max(0.0);
            }
        } else {
            out.terminal_lakes += 1;
            out.terminal[b] = true;
            out.level_m[b] = elevation_m[cs[k]];
        }
    }
    Ok(out)
}

/// Step 3: the water level of every cell (m), `f64::NEG_INFINITY` where the
/// cell is dry (land, an exposed terminal-lake floor). Marine cells read the
/// sea level.
pub fn cell_water_level_m(routing: &Routing, lb: &LakeBalance, elevation_m: &[f64], sea_level_m: f64) -> Vec<f64> {
    (0..elevation_m.len())
        .map(|c| {
            if routing.ocean[c] {
                return sea_level_m;
            }
            let b = routing.water_body_id[c];
            if b == crate::water::DRY {
                return f64::NEG_INFINITY;
            }
            let l = lb.level_m[b as usize];
            if l > elevation_m[c] {
                l
            } else {
                f64::NEG_INFINITY
            }
        })
        .collect()
}

/// Per-cell net open-water evaporation (m/yr) and runoff (m/yr) for the
/// balance: the climate law's fields when given (runoff m/yr, PET and P in
/// mm/yr), else the uniform runoff and the declared ratio.
pub fn balance_inputs(
    n: usize,
    uniform_runoff_m_per_yr: f64,
    climate: Option<(&[f64], &[f64], &[f64])>,
    p: &ClosedBasinParams,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    match climate {
        Some((runoff, pet_mm, precip_mm)) => {
            if runoff.len() != n || pet_mm.len() != n || precip_mm.len() != n {
                return Err(format!("{CLOSED_BASIN_LAW_ID}: climate fields do not match the grid"));
            }
            let e = (0..n).map(|c| (pet_mm[c] - precip_mm[c]) / 1000.0).collect();
            Ok((runoff.to_vec(), e))
        }
        None => Ok((
            vec![uniform_runoff_m_per_yr; n],
            vec![p.uniform_evaporation_to_runoff_ratio * uniform_runoff_m_per_yr; n],
        )),
    }
}

/// One step's water-balance outcome kept for observation.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LakeObservation {
    /// Per-cell water level (m), `f64::NEG_INFINITY` where dry.
    pub water_level_m: Vec<f64>,
    pub full_lakes: u64,
    pub terminal_lakes: u64,
    pub wet_area_m2: f64,
    pub exposed_area_m2: f64,
}

thread_local! {
    /// DIAGNOSTIC ONLY (not run state, never in an output or a digest): the
    /// latest step's lake observation on this thread, so a measurement
    /// driver can map lakes without the run carrying them.
    static LAST: std::cell::RefCell<Option<LakeObservation>> = const { std::cell::RefCell::new(None) };
}

/// Diagnostic: record this step's observation (the preview step calls it).
pub fn record_observation(o: LakeObservation) {
    LAST.with(|t| *t.borrow_mut() = Some(o));
}

/// Diagnostic: take (and clear) the latest observation on this thread.
pub fn take_observation() -> Option<LakeObservation> {
    LAST.with(|t| t.borrow_mut().take())
}

#[cfg(test)]
mod tests {
    //! SYNTHETIC FIXTURES (instrument checks, never world input): an ocean
    //! strip, a land ramp and one closed pit below sea level in the land.
    use super::*;
    use crate::erosion::retain::{erode_retain_step_with_lakes, RetentionParams};
    use crate::erosion::{DiffusionStencil, ErosionParams};
    use crate::routing::route;
    use crate::sediment::CompensatedColumn;
    use pf_grid::GridSpec;

    const COLUMN: CompensatedColumn = CompensatedColumn {
        mantle_density_kg_m3: 3300.0,
        sediment_density_kg_m3: 2800.0,
        water_density_kg_m3: 1030.0,
    };

    /// 16 x 32 cells; rows 0-2 ocean at -3000 m; land from 200 m rising
    /// 300 m per row; a 3 x 3 pit at -150 m inside a 5 x 5 rim held at
    /// 400 m (rows 7-11, cols 10-14).
    fn world() -> (GridSpec, Vec<f64>, Vec<bool>, Vec<usize>) {
        let grid = GridSpec::new(16, 32, 6.371e6).unwrap();
        let (nlat, nlon) = (16, 32);
        let mut z = vec![0.0; nlat * nlon];
        let mut oc = vec![false; nlat * nlon];
        let mut pit = Vec::new();
        for i in 0..nlat {
            for j in 0..nlon {
                let c = i * nlon + j;
                if i < 3 {
                    z[c] = -3000.0;
                    oc[c] = true;
                } else {
                    z[c] = 200.0 + 300.0 * (i - 3) as f64 + (j % 3) as f64 * 0.5;
                }
                if (7..12).contains(&i) && (10..15).contains(&j) {
                    z[c] = 400.0 + (i + j) as f64 * 0.01;
                }
                if (8..11).contains(&i) && (11..14).contains(&j) {
                    z[c] = -150.0 - (i * 3 + j) as f64 * 0.1;
                    pit.push(c);
                }
            }
        }
        (grid, z, oc, pit)
    }

    #[test]
    fn balance_full_when_wet_terminal_when_evaporation_wins() {
        let (grid, z, oc, pit) = world();
        let n = grid.cells();
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let areas: Vec<f64> = (0..n).map(|c| grid.cell_area_m2(c / 32).unwrap()).collect();
        let b = r.water_body_id[pit[0]] as usize;
        let spill = r.water_bodies[b].water_surface_m;
        assert!(spill > 0.0);
        // No net evaporation: full at the spill level.
        let lb = lake_balance(&r, &z, &areas, &vec![1.0; n], &vec![0.0; n]).unwrap();
        assert!(!lb.terminal[b]);
        assert_eq!(lb.level_m[b], spill);
        // Evaporation 1000 x runoff: terminal, level between floor and spill.
        let lb = lake_balance(&r, &z, &areas, &vec![1.0e-3; n], &vec![1.0; n]).unwrap();
        assert!(lb.terminal[b]);
        let floor = pit.iter().map(|&c| z[c]).fold(f64::INFINITY, f64::min);
        assert!(lb.level_m[b] >= floor && lb.level_m[b] < spill, "level {}", lb.level_m[b]);
        let wl = cell_water_level_m(&r, &lb, &z, 0.0);
        assert!(pit.iter().any(|&c| !wl[c].is_finite()) || lb.level_m[b] > floor);
    }

    #[test]
    fn focused_fill_conserves_mass_and_reaches_the_deepest_floor() {
        let (grid, z, oc, pit) = world();
        let n = grid.cells();
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let p = ErosionParams { k_per_my: 2.0, m: 0.5, n: 1.0, diffusivity_m2_per_my: 0.0, stability_limit: 0.5 };
        let rp = RetentionParams::v1(COLUMN);
        let cont = vec![true; n];
        let sed = vec![0.0; n];
        let rock = vec![30_000.0; n];
        let lb = lake_balance(&r, &z, &stencil.areas_m2, &vec![1.0; n], &vec![1.5; n]).unwrap();
        let cp = ClosedBasinParams::v1();
        let run = |l: Option<(&LakeBalance, &ClosedBasinParams)>| {
            erode_retain_step_with_lakes(&grid, &stencil, &r, &z, &sed, &oc, 0.0, 2.0, &p, Some(&rock), &cont, &rp, None, false, l).unwrap()
        };
        let (off, _, _) = run(None);
        let (on, rr, _) = run(Some((&lb, &cp)));
        for (name, s) in [("off", &off), ("on", &on)] {
            let tol = 1.0e-9 * s.eroded_m3.max(1.0);
            assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: mass");
        }
        assert_eq!(off.bedrock_removed_m, on.bedrock_removed_m, "the law changes placement only");
        assert!(rr.lake_m3 > 0.0);
        let deepest = *pit.iter().min_by(|&&a, &&b| z[a].total_cmp(&z[b])).unwrap();
        assert!(on.sediment_added_m[deepest] > 0.0, "focused fill reaches the deepest cell");
        assert!(on.sediment_added_m[deepest] >= off.sediment_added_m[deepest]);
    }

    /// LAKE-EXCESS (`filled-basin-alluvial-grade-v1`): `world()`'s pit lake
    /// (rim and pit, a closed body below the ramp) with every body cell set
    /// to `depth_m` below the spill level (None: `world()` unchanged, the
    /// pit ~1.25 km deep); returns the step under v1 and under the graded
    /// variant, the surface, the body's cells and the spill level.
    #[allow(clippy::type_complexity)]
    fn graded_pair(depth_m: Option<f64>) -> (
        crate::erosion::ErosionStep,
        crate::erosion::retain::RetentionReport,
        crate::erosion::ErosionStep,
        crate::erosion::retain::RetentionReport,
        Vec<f64>,
        Vec<usize>,
        f64,
    ) {
        let (grid, mut z, oc, pit) = world();
        let n = grid.cells();
        let r0 = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let b0 = r0.water_body_id[pit[0]];
        let spill = r0.water_bodies[b0 as usize].water_surface_m;
        let body: Vec<usize> = (0..n).filter(|&c| r0.water_body_id[c] == b0).collect();
        if let Some(d) = depth_m {
            for &c in &body {
                z[c] = spill - d - (c % 7) as f64 * 0.01;
            }
        }
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        assert_eq!(r.water_bodies[r.water_body_id[pit[0]] as usize].water_surface_m, spill);
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let p = ErosionParams { k_per_my: 2.0, m: 0.5, n: 1.0, diffusivity_m2_per_my: 0.0, stability_limit: 0.5 };
        // The fixture's cells are ~1250 km across: S_a = 1e-5 keeps the
        // graded relief over the 5 x 5 body (~90 m) within one step's supply.
        let mut rp = RetentionParams::v1(COLUMN);
        rp.alluvial_grade_slope = 1.0e-5;
        let cont = vec![true; n];
        let sed = vec![0.0; n];
        let rock = vec![30_000.0; n];
        let lb = lake_balance(&r, &z, &stencil.areas_m2, &vec![1.0; n], &vec![0.0; n]).unwrap();
        let run = |cp: ClosedBasinParams| {
            erode_retain_step_with_lakes(&grid, &stencil, &r, &z, &sed, &oc, 0.0, 2.0, &p, Some(&rock), &cont, &rp, None, false, Some((&lb, &cp))).unwrap()
        };
        let (v1, rr1, _) = run(ClosedBasinParams::v1());
        let (g, rrg, _) = run(ClosedBasinParams::v1_graded());
        for (name, s) in [("v1", &v1), ("graded", &g)] {
            let tol = 1.0e-9 * s.eroded_m3.max(1.0);
            assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: mass");
        }
        (v1, rr1, g, rrg, z, body, spill)
    }

    #[test]
    fn graded_variant_leaves_a_deep_lake_to_v1() {
        // The unchanged pit: the step's fill stays far below spill - d_c.
        let (v1, rr1, g, rrg, ..) = graded_pair(None);
        assert_eq!(rrg.graded_basins, 0);
        assert_eq!(v1.sediment_added_m, g.sediment_added_m);
        assert_eq!(rr1.lake_m3, rrg.lake_m3);
    }

    #[test]
    fn graded_variant_turns_a_nearly_full_lake_into_a_plain_above_its_sill() {
        // A 3 m deep lake (< d_c = 10 m) with ample supply: v1 fills it level
        // to the spill; the graded variant raises every pit cell above the
        // spill by S_a x its flow distance to the outlet, so the depression
        // is gone and the plain drains across the sill.
        let (v1, _, g, rrg, z, body, spill) = graded_pair(Some(3.0));
        assert!(rrg.graded_basins >= 1 && rrg.graded_basin_m3 > 0.0);
        for &c in &body {
            let s1 = COLUMN.surface_after(z[c], v1.sediment_added_m[c], 0.0);
            let sg = COLUMN.surface_after(z[c], g.sediment_added_m[c], 0.0);
            assert!((s1 - spill).abs() < 1.0e-6, "v1 fills level to the spill: {s1} vs {spill}");
            assert!(sg > spill + 5.0, "graded plain stands above the sill: {sg} vs {spill}");
        }
    }

    /// LAKE-OUTLET (`lake-outlet-sill-incision-v1`): `world()`'s pit lake
    /// under v1 and with the outlet law (n = 2 channel law with a
    /// threshold, as the C7 recipe), for net evaporation `e` (m/yr).
    /// The step's own K is set small (1e-4 /My) so its lowering of the pour
    /// cell stands in for the slow subgrid closure and taper of the real
    /// runs; the outlet law carries the C7 channel law (K 2 /My).
    fn outlet_pair(e: f64) -> (crate::erosion::ErosionStep, crate::erosion::ErosionStep, crate::erosion::retain::RetentionReport, usize, Vec<f64>) {
        let (grid, z, oc, pit) = world();
        let n = grid.cells();
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let p = ErosionParams { k_per_my: 1.0e-4, m: 0.5, n: 1.0, diffusivity_m2_per_my: 0.0, stability_limit: 0.5 };
        let rp = RetentionParams::v1(COLUMN);
        let cont = vec![true; n];
        let sed = vec![0.0; n];
        let rock = vec![30_000.0; n];
        let lb = lake_balance(&r, &z, &stencil.areas_m2, &vec![1.0; n], &vec![e; n]).unwrap();
        let b = r.water_body_id[pit[0]] as usize;
        let sa = sill_discharge_area_m2(&r, &lb, 1.0);
        let sill = (0..n).find(|&c| sa[c] > 0.0).unwrap_or(r.water_bodies[b].outlet_edges[0][1].0 as usize);
        assert_eq!(r.water_body_id[sill], crate::water::DRY);
        assert_eq!(z[sill], r.water_bodies[b].water_surface_m, "the pour cell stands at the spill level");
        let off = ClosedBasinParams::v1();
        let on = ClosedBasinParams {
            outlet: Some(OutletIncisionParams {
                reference_runoff_m_per_yr: 1.0,
                channel: Some(OutletChannelLaw { k_n_per_my: 2.0 / 500.0, m: 1.0, n: 2.0, threshold_rate_m_per_my: 2.0 / 500.0 * 60.0f64.powi(2) }),
                flexural_support: false,
                path_incision: false,
                through_lakes: false,
                path_grade: false,
                alluvial_grade_slope: None,
            }),
            ..ClosedBasinParams::v1()
        };
        let run = |cp: ClosedBasinParams| {
            erode_retain_step_with_lakes(&grid, &stencil, &r, &z, &sed, &oc, 0.0, 2.0, &p, Some(&rock), &cont, &rp, None, false, Some((&lb, &cp))).unwrap()
        };
        let (s0, _, _) = run(off);
        let (s1, rr1, _) = run(on);
        for (name, s) in [("off", &s0), ("on", &s1)] {
            let tol = 1.0e-9 * s.eroded_m3.max(1.0);
            assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: mass");
        }
        (s0, s1, rr1, sill, lb.overflow_m3_per_yr.clone())
    }

    #[test]
    fn outlet_law_incises_the_pour_cell_of_a_full_lake_only() {
        let (s0, s1, rr, sill, q) = outlet_pair(0.0);
        assert!(q.iter().any(|&v| v > 0.0), "the pit lake overflows");
        assert!(rr.outlet_sills >= 1);
        assert!(s1.bedrock_removed_m[sill] > s0.bedrock_removed_m[sill], "{} vs {}", s1.bedrock_removed_m[sill], s0.bedrock_removed_m[sill]);
        for c in 0..s0.bedrock_removed_m.len() {
            if c != sill {
                assert!(s1.bedrock_removed_m[c] == s0.bedrock_removed_m[c], "only the pour cell's own lowering changes ({c})");
            }
        }
    }

    #[test]
    fn outlet_law_leaves_a_terminal_lake_alone() {
        // Evaporation 1000 m/yr: the lake is terminal, nothing overflows.
        let (s0, s1, rr, _, q) = outlet_pair(1000.0);
        assert!(q.iter().all(|&v| v == 0.0));
        assert_eq!(rr.outlet_sills, 0);
        assert_eq!(s0.bedrock_removed_m, s1.bedrock_removed_m);
        assert_eq!(s0.sediment_added_m, s1.sediment_added_m);
    }

    /// SILL-REBOUND (rule 6): with flexural support the sill's working
    /// surface drops exactly as under rule 5, but the material removed is
    /// the compensated column for that drop (no local rebound); mass closes.
    #[test]
    fn flexural_support_removes_the_compensated_column_at_the_sill() {
        let (grid, z, oc, _) = world();
        let n = grid.cells();
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let p = ErosionParams { k_per_my: 1.0e-4, m: 0.5, n: 1.0, diffusivity_m2_per_my: 0.0, stability_limit: 0.5 };
        let rp = RetentionParams::v1(COLUMN);
        let (cont, sed, rock) = (vec![true; n], vec![0.0; n], vec![30_000.0; n]);
        let lb = lake_balance(&r, &z, &stencil.areas_m2, &vec![1.0; n], &vec![0.0; n]).unwrap();
        let sa = sill_discharge_area_m2(&r, &lb, 1.0);
        let sill = (0..n).find(|&c| sa[c] > 0.0).expect("a sill");
        let law = |flex: bool| ClosedBasinParams {
            outlet: Some(OutletIncisionParams {
                reference_runoff_m_per_yr: 1.0,
                channel: Some(OutletChannelLaw { k_n_per_my: 2.0 / 500.0, m: 1.0, n: 2.0, threshold_rate_m_per_my: 2.0 / 500.0 * 60.0f64.powi(2) }),
                flexural_support: flex,
                path_incision: false,
                through_lakes: false,
                path_grade: false,
                alluvial_grade_slope: None,
            }),
            ..ClosedBasinParams::v1()
        };
        let run = |cp: ClosedBasinParams| {
            erode_retain_step_with_lakes(&grid, &stencil, &r, &z, &sed, &oc, 0.0, 2.0, &p, Some(&rock), &cont, &rp, None, false, Some((&lb, &cp))).unwrap()
        };
        let (s5, r5, _) = run(law(false));
        let (s6, r6, _) = run(law(true));
        let drop = -s5.dz_m[sill];
        assert!(drop > 1.0, "rule 5 lowers the sill ({drop})");
        assert!((s6.dz_m[sill] - s5.dz_m[sill]).abs() < 1e-9, "same working-surface drop");
        let removed5 = s5.bedrock_removed_m[sill] + s5.sediment_removed_m[sill];
        let removed6 = s6.bedrock_removed_m[sill] + s6.sediment_removed_m[sill];
        assert!((removed5 - drop).abs() < 1e-9, "rule 5 removes the drop itself");
        let want = COLUMN.thickness_to_raise(z[sill] - drop, z[sill], 0.0);
        assert!((removed6 - want).abs() < 1e-9 * want, "{removed6} vs {want}");
        assert!((removed6 * COLUMN.rise_per_m_subaerial() - drop).abs() < 1e-9 * drop, "post-isostasy drop = the lowering (sill above sea level)");
        assert!(r5.outlet_support_extra_m3 == 0.0 && r6.outlet_support_extra_m3 > 0.0);
        for (name, s) in [("rule 5", &s5), ("rule 6", &s6)] {
            let tol = 1.0e-9 * s.eroded_m3.max(1.0);
            assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: mass");
        }
    }

    /// LAKE-HOLD (rule 7): with `path_incision` every dry cell below the
    /// sill down to the ocean takes the channel law on the overflow (the
    /// step's own K is tiny, as in `outlet_pair`); cells off the path keep
    /// their lowering; a terminal lake has no path; mass closes.
    #[test]
    fn path_incision_lowers_the_whole_overflow_path_only() {
        let (grid, z, oc, _) = world();
        let n = grid.cells();
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let p = ErosionParams { k_per_my: 1.0e-4, m: 0.5, n: 1.0, diffusivity_m2_per_my: 0.0, stability_limit: 0.5 };
        let rp = RetentionParams::v1(COLUMN);
        let (cont, sed, rock) = (vec![true; n], vec![0.0; n], vec![30_000.0; n]);
        let law = |path: bool| ClosedBasinParams {
            outlet: Some(OutletIncisionParams {
                reference_runoff_m_per_yr: 1.0,
                channel: Some(OutletChannelLaw { k_n_per_my: 2.0 / 500.0, m: 1.0, n: 2.0, threshold_rate_m_per_my: 2.0 / 500.0 * 60.0f64.powi(2) }),
                flexural_support: false,
                path_incision: path,
                through_lakes: false,
                path_grade: false,
                alluvial_grade_slope: None,
            }),
            ..ClosedBasinParams::v1()
        };
        for e in [0.0, 1000.0] {
            let lb = lake_balance(&r, &z, &stencil.areas_m2, &vec![1.0; n], &vec![e; n]).unwrap();
            let sill_only = sill_discharge_area_m2(&r, &lb, 1.0);
            let a = outlet_discharge_area_m2(&r, &lb, 1.0, true);
            assert_eq!(outlet_discharge_area_m2(&r, &lb, 1.0, false), sill_only, "path = false is rule 5's sill area");
            let run = |cp: ClosedBasinParams| {
                erode_retain_step_with_lakes(&grid, &stencil, &r, &z, &sed, &oc, 0.0, 2.0, &p, Some(&rock), &cont, &rp, None, false, Some((&lb, &cp))).unwrap()
            };
            let (s5, _, _) = run(law(false));
            let (s7, _, _) = run(law(true));
            for (name, s) in [("rule 5", &s5), ("rule 7", &s7)] {
                let tol = 1.0e-9 * s.eroded_m3.max(1.0);
                assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: mass");
            }
            if e > 0.0 {
                // Terminal lake: no overflow, no sill, no path.
                assert!(a.iter().all(|&v| v == 0.0));
                assert_eq!(s5.bedrock_removed_m, s7.bedrock_removed_m);
                continue;
            }
            let sill = (0..n).find(|&c| sill_only[c] > 0.0).expect("a sill");
            // The path: the sill's flow tree down to the ocean.
            let mut path = Vec::new();
            let mut c = sill;
            while let Some(j) = r.flow_to[c] {
                c = j.0 as usize;
                if oc[c] {
                    break;
                }
                path.push(c);
            }
            assert!(!path.is_empty(), "the sill drains over land to the ocean");
            for c in 0..n {
                let on_path = path.contains(&c);
                assert_eq!(a[c] > 0.0, on_path || c == sill, "A_q only on the sill and its path ({c})");
                if on_path || c == sill {
                    assert_eq!(a[c], sill_only[sill], "the same A_q along the path");
                } else {
                    assert!(s7.bedrock_removed_m[c] == s5.bedrock_removed_m[c], "off-path lowering unchanged ({c})");
                }
            }
            assert!(path.iter().all(|&c| s7.bedrock_removed_m[c] >= s5.bedrock_removed_m[c]));
            assert!(path.iter().any(|&c| s7.bedrock_removed_m[c] > s5.bedrock_removed_m[c] + 1.0), "the path is cut");
        }
    }

    /// LAKE-HOLD-2 (rule 8): `world()` plus a small trench lake B across
    /// row 4 on the pit lake A's overflow path (A -> sill -> path -> B ->
    /// B's sill -> ocean). With `through_lakes` the cells entering B are cut
    /// to B's post-incision spill (B's sill is cut in the same step), so the
    /// path cell above B is lowered more than under rule 7; only cells that
    /// drain over land into a lake with a sill (A or B) change their
    /// removal; mass closes; a terminal
    /// world (no overflow) is rule 7 exactly; off is rule 7 exactly.
    #[test]
    fn through_lakes_cuts_into_a_lake_lowered_in_the_same_step() {
        let (grid, mut z, oc, _) = world();
        let n = grid.cells();
        for j in 8..=16 {
            z[4 * 32 + j] = 100.0 + j as f64 * 0.01;
        }
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let p = ErosionParams { k_per_my: 1.0e-4, m: 0.5, n: 1.0, diffusivity_m2_per_my: 0.0, stability_limit: 0.5 };
        let rp = RetentionParams::v1(COLUMN);
        let (cont, sed, rock) = (vec![true; n], vec![0.0; n], vec![30_000.0; n]);
        let law = |through: bool| ClosedBasinParams {
            outlet: Some(OutletIncisionParams {
                reference_runoff_m_per_yr: 1.0,
                channel: Some(OutletChannelLaw { k_n_per_my: 2.0 / 500.0, m: 1.0, n: 2.0, threshold_rate_m_per_my: 2.0 / 500.0 * 60.0f64.powi(2) }),
                flexural_support: false,
                path_incision: true,
                through_lakes: through,
                path_grade: false,
                alluvial_grade_slope: None,
            }),
            ..ClosedBasinParams::v1()
        };
        assert!(OutletIncisionParams { path_incision: false, ..law(true).outlet.unwrap() }.validate().is_err(), "rule 8 needs rule 7");
        let b_lake = r.water_body_id[4 * 32 + 12];
        assert_ne!(b_lake, crate::water::DRY, "row 4 trench is a lake");
        for e in [0.0, 1000.0] {
            let lb = lake_balance(&r, &z, &stencil.areas_m2, &vec![1.0; n], &vec![e; n]).unwrap();
            let run = |cp: ClosedBasinParams| {
                erode_retain_step_with_lakes(&grid, &stencil, &r, &z, &sed, &oc, 0.0, 2.0, &p, Some(&rock), &cont, &rp, None, false, Some((&lb, &cp))).unwrap()
            };
            let (s7, _, _) = run(law(false));
            let (s8, _, _) = run(law(true));
            for (name, s) in [("rule 7", &s7), ("rule 8", &s8)] {
                let tol = 1.0e-9 * s.eroded_m3.max(1.0);
                assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "{name}: mass");
            }
            if e > 0.0 {
                assert_eq!(s7.bedrock_removed_m, s8.bedrock_removed_m, "no overflow: rule 8 is rule 7");
                assert_eq!(s7.sediment_added_m, s8.sediment_added_m);
                continue;
            }
            let sills = outlet_sill_cells(&r, &lb);
            let sb = sills[b_lake as usize].expect("lake B overflows over a sill");
            assert!(s8.bedrock_removed_m[sb] > 1.0, "B's sill is cut");
            let into_b = |c: usize| r.flow_to[c].is_some_and(|j| r.water_body_id[j.0 as usize] == b_lake && !oc[c] && r.water_body_id[c] == crate::water::DRY);
            // A dry cell whose receiver is a lake with a sill (A or B).
            let into_sill_lake = |c: usize| {
                r.flow_to[c].is_some_and(|j| {
                    let b = r.water_body_id[j.0 as usize];
                    !oc[c] && r.water_body_id[c] == crate::water::DRY && b != crate::water::DRY && sills.get(b as usize).is_some_and(|x| x.is_some())
                })
            };
            // Such cells and their donors (the implicit update sees a
            // receiver's new surface, so a deeper cut reaches upstream).
            let upstream_of_entry = |c: usize| {
                let mut k = c;
                for _ in 0..n {
                    if into_sill_lake(k) {
                        return true;
                    }
                    match r.flow_to[k] {
                        Some(j) if r.water_body_id[j.0 as usize] == crate::water::DRY && !oc[j.0 as usize] => k = j.0 as usize,
                        _ => return false,
                    }
                }
                false
            };
            let mut deeper = 0;
            for c in 0..n {
                if upstream_of_entry(c) {
                    assert!(s8.bedrock_removed_m[c] >= s7.bedrock_removed_m[c] - 1e-9, "a lower base never cuts less ({c})");
                    if into_b(c) && s8.bedrock_removed_m[c] > s7.bedrock_removed_m[c] + 1.0 {
                        deeper += 1;
                    }
                } else {
                    assert!(s8.bedrock_removed_m[c] == s7.bedrock_removed_m[c], "only cells draining over land into a lake with a sill change ({c})");
                }
            }
            assert!(deeper >= 1, "the cut passes into lake B within the step");
            // A's path cell entering B carries A's overflow (rule 7's A_q).
            let a = outlet_discharge_area_m2(&r, &lb, 1.0, true);
            assert!((0..n).any(|c| into_b(c) && a[c] > 0.0 && s8.bedrock_removed_m[c] > s7.bedrock_removed_m[c] + 1.0), "A's overflow path is cut deeper");
        }
    }

    /// LAKE-HOLD-2 (rule 9): the outlet-path grade is min(S_a, k_c /
    /// A_q^(m/n)) with k_c recovered from the channel law; S_a without the
    /// flag, off the path or without a threshold; validation refuses it
    /// without rule 7 or a threshold; on `world()` the path stores no more
    /// than under rule 7 and mass closes.
    #[test]
    fn path_grade_is_the_channel_threshold_slope_on_the_path_only() {
        let ch = OutletChannelLaw { k_n_per_my: 2.0 / 500.0, m: 1.0, n: 2.0, threshold_rate_m_per_my: 2.0 / 500.0 * 60.0f64.powi(2) };
        let o = OutletIncisionParams { reference_runoff_m_per_yr: 1.0, channel: Some(ch), flexural_support: false, path_incision: true, through_lakes: false, path_grade: true, alluvial_grade_slope: None };
        assert!(o.validate().is_ok());
        assert!((o.path_grade_slope(1.0e-4, 1.0e12) - 60.0 / 1.0e6).abs() < 1e-15, "k_c / sqrt(A_q) below S_a");
        assert_eq!(o.path_grade_slope(1.0e-4, 1.0e10), 1.0e-4, "S_a when the threshold slope is steeper");
        assert_eq!(o.path_grade_slope(1.0e-4, 0.0), 1.0e-4, "off the path");
        assert_eq!(OutletIncisionParams { path_grade: false, ..o }.path_grade_slope(1.0e-4, 1.0e12), 1.0e-4);
        assert!(OutletIncisionParams { path_incision: false, ..o }.validate().is_err());
        assert!(OutletIncisionParams { channel: Some(OutletChannelLaw { threshold_rate_m_per_my: 0.0, ..ch }), ..o }.validate().is_err());
        assert!(OutletIncisionParams { channel: None, ..o }.validate().is_err());
        let (grid, z, oc, _) = world();
        let n = grid.cells();
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let p = ErosionParams { k_per_my: 2.0, m: 0.5, n: 1.0, diffusivity_m2_per_my: 0.0, stability_limit: 0.5 };
        let rp = RetentionParams::v1(COLUMN);
        let (cont, sed, rock) = (vec![true; n], vec![0.0; n], vec![30_000.0; n]);
        let lb = lake_balance(&r, &z, &stencil.areas_m2, &vec![1.0; n], &vec![0.0; n]).unwrap();
        let a = outlet_discharge_area_m2(&r, &lb, 1.0, true);
        let run = |grade: bool| {
            let cp = ClosedBasinParams { outlet: Some(OutletIncisionParams { path_grade: grade, ..o }), ..ClosedBasinParams::v1() };
            erode_retain_step_with_lakes(&grid, &stencil, &r, &z, &sed, &oc, 0.0, 2.0, &p, Some(&rock), &cont, &rp, None, false, Some((&lb, &cp))).unwrap().0
        };
        let (s7, s9) = (run(false), run(true));
        for s in [&s7, &s9] {
            let tol = 1.0e-9 * s.eroded_m3.max(1.0);
            assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "mass");
        }
        for c in 0..n {
            if a[c] > 0.0 {
                assert!(s9.sediment_added_m[c] <= s7.sediment_added_m[c] + 1e-9, "the path stores no more ({c})");
            }
            assert_eq!(s9.bedrock_removed_m[c], s7.bedrock_removed_m[c], "stage 1 unchanged ({c})");
        }
    }

    /// ALLUVIAL-TRUNK (rule 10): the outlet law drops its threshold exactly
    /// where k_c / A_q^(m/n) < S_a (A_q > A_a = 3.6e11 m2 at k_c 60 m, S_a
    /// 1e-4): a gentle reach held below the threshold is cut there and
    /// nowhere else; validation refuses it without a threshold; on `world()`
    /// with rule 7 the path never cuts less, every other cell is unchanged,
    /// and mass closes.
    #[test]
    fn alluvial_outlet_channel_drops_the_threshold_above_the_switch_area_only() {
        let ch = OutletChannelLaw { k_n_per_my: 2.0 / 500.0, m: 1.0, n: 2.0, threshold_rate_m_per_my: 2.0 / 500.0 * 60.0f64.powi(2) };
        let o = OutletIncisionParams { reference_runoff_m_per_yr: 1.0, channel: Some(ch), flexural_support: false, path_incision: true, through_lakes: false, path_grade: false, alluvial_grade_slope: None };
        let a = OutletIncisionParams { alluvial_grade_slope: Some(1.0e-4), ..o };
        assert!(a.validate().is_ok());
        let a_a = 3.6e11;
        assert!(!a.alluvial_channel(0.99 * a_a) && a.alluvial_channel(1.01 * a_a) && !o.alluvial_channel(1.0e13));
        // 3 m over 70 km at A_q = 1e12: k_s = 3/7e4 x 1e6 = 43 m < k_c, held.
        assert_eq!(o.sill_lowering_m(3.0, 1.0e12, 7.0e4, 2.0, 2.0, 0.5), 0.0);
        let d = a.sill_lowering_m(3.0, 1.0e12, 7.0e4, 2.0, 2.0, 0.5);
        assert!(d > 1.0 && d < 3.0, "{d}");
        let free = OutletIncisionParams { channel: Some(OutletChannelLaw { threshold_rate_m_per_my: 0.0, ..ch }), ..o };
        assert_eq!(d, free.sill_lowering_m(3.0, 1.0e12, 7.0e4, 2.0, 2.0, 0.5), "the k_c = 0 channel law");
        // Below A_a: the thresholded law bit for bit.
        assert_eq!(a.sill_lowering_m(3.0, 3.0e11, 7.0e4, 2.0, 2.0, 0.5), o.sill_lowering_m(3.0, 3.0e11, 7.0e4, 2.0, 2.0, 0.5));
        assert_eq!(a.sill_lowering_m(30.0, 3.0e11, 7.0e4, 2.0, 2.0, 0.5), o.sill_lowering_m(30.0, 3.0e11, 7.0e4, 2.0, 2.0, 0.5));
        assert!(OutletIncisionParams { channel: Some(OutletChannelLaw { threshold_rate_m_per_my: 0.0, ..ch }), ..a }.validate().is_err());
        assert!(OutletIncisionParams { channel: None, ..a }.validate().is_err());
        assert!(OutletIncisionParams { alluvial_grade_slope: Some(-1.0), ..a }.validate().is_err());
        // On world(): rule 7 vs rule 7 + rule 10.
        let (grid, z, oc, _) = world();
        let n = grid.cells();
        let r = route(&grid, &z, &vec![1.0e6; n], &oc, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let p = ErosionParams { k_per_my: 1.0e-4, m: 0.5, n: 1.0, diffusivity_m2_per_my: 0.0, stability_limit: 0.5 };
        let rp = RetentionParams::v1(COLUMN);
        let (cont, sed, rock) = (vec![true; n], vec![0.0; n], vec![30_000.0; n]);
        let lb = lake_balance(&r, &z, &stencil.areas_m2, &vec![1.0; n], &vec![0.0; n]).unwrap();
        let aq = outlet_discharge_area_m2(&r, &lb, 1.0, true);
        assert!(aq.iter().any(|&v| v > a_a), "the fixture's outlet path is alluvial");
        let run = |x: OutletIncisionParams| {
            let cp = ClosedBasinParams { outlet: Some(x), ..ClosedBasinParams::v1() };
            erode_retain_step_with_lakes(&grid, &stencil, &r, &z, &sed, &oc, 0.0, 2.0, &p, Some(&rock), &cont, &rp, None, false, Some((&lb, &cp))).unwrap().0
        };
        let (s7, s10) = (run(o), run(a));
        for s in [&s7, &s10] {
            let tol = 1.0e-9 * s.eroded_m3.max(1.0);
            assert!((s.eroded_m3 - s.deposited_m3 - s.exported_m3).abs() <= tol, "mass");
        }
        for c in 0..n {
            if aq[c] > 0.0 {
                assert!(s10.bedrock_removed_m[c] >= s7.bedrock_removed_m[c] - 1e-9, "a released channel never cuts less ({c})");
            } else {
                assert_eq!(s10.bedrock_removed_m[c], s7.bedrock_removed_m[c], "stage 1 off the outlet path unchanged ({c})");
            }
        }
        // Rule 10 off (None) is rule 7 exactly.
        assert_eq!(run(OutletIncisionParams { alluvial_grade_slope: None, ..a }), s7);
    }

    #[test]
    fn sill_lowering_stops_at_the_receiver_and_below_the_threshold() {
        let o = OutletIncisionParams {
            reference_runoff_m_per_yr: 1.0,
            channel: Some(OutletChannelLaw { k_n_per_my: 4.0e-3, m: 1.0, n: 2.0, threshold_rate_m_per_my: 4.0e-3 * 3600.0 }),
            flexural_support: false,
            path_incision: false,
            through_lakes: false,
            path_grade: false,
            alluvial_grade_slope: None,
        };
        // Huge discharge: lowered almost to the receiver, never below.
        let d = o.sill_lowering_m(50.0, 1.0e14, 7.0e4, 2.0, 2.0, 0.5);
        assert!(d > 49.0 && d <= 50.0, "{d}");
        // Tiny discharge on a gentle sill: steepness below k_c, no lowering.
        assert_eq!(o.sill_lowering_m(1.0, 1.0e6, 7.0e4, 2.0, 2.0, 0.5), 0.0);
        // No relief law: the linear law (n = 1) is x0 f / (1 + f).
        let l = OutletIncisionParams { reference_runoff_m_per_yr: 1.0, channel: None, flexural_support: false, path_incision: false, through_lakes: false, path_grade: false, alluvial_grade_slope: None };
        let f = 2.0 * 1.0e10f64.sqrt() * 2.0 / 7.0e4;
        assert!((l.sill_lowering_m(10.0, 1.0e10, 7.0e4, 2.0, 2.0, 0.5) - 10.0 * f / (1.0 + f)).abs() < 1e-9);
    }
}
