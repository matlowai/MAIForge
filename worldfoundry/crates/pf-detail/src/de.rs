//! Lane DE (2026-10-15; opt-in `de=1`, needs version=3 and beds=2): differential erosion by hardness as a PROCESS.
//! PREREG: `unpublished-design-record` (design, bug fixes B1-B3 and revision R1 logged there).
//!
//! design 2026-10-15: "some geological variation between layers based on erosion and hardness ... just be inspired by nature".
//!
//! DE1 layered hillslope equilibrium (weathering-limited transport). Hillslope sediment flux between 4-neighbours q = D S with a
//! rock-controlled transport coefficient D = D0 / Sc (weak rock weathers and is stripped faster; beds=2 keeps K Sc = class K Sc,
//! so D is proportional to the bed erodibility K). An edge on a steep wall spans several beds: its D is the height mean over the
//! span (beds in series, B1). Soil mantle: the bed contrast enters with the exposure e (gentle ground is soil-covered and does
//! not see the bed; e from the slope alone, B4): 1/Sc_eff = (1 - e)/Sc_w + e/Sc_b (soil flux + rock flux, B2). Reference rock = the window-mean rock (Sc_w,
//! the 60 m-section harmonic mean): the forcing s is the flux divergence that keeps the incoming surface Z0 steady under the
//! reference rock (the erosion / uplift balance T1b expresses at its scale). The layered surface z solves
//! div(D_eff grad z) - lambda (z - Z0) + s = 0 on free posts, z = Z0 on fixed posts; the screening lambda (revision R1) keeps
//! the change below the T1b resolution (screen_m). With no contrast z = Z0. In 1D flux continuity gives S_b / S_b' = Sc_b / Sc_b':
//! hard beds stand steep (ledges, cliffs), soft beds below a hard cap flatten (benches); the remap's rule is the RESULT here.
//! B3 (logged): the slope factor F = 1/(1-(S/Sc)^2) of the B5 creep law made the Picard iteration diverge on walls near their
//! critical slope; the threshold is enforced instead by DE2's failure step (linear transport + threshold failure).
//! DE2 talus: (a) mantle: posts below a bare hard face (within the face's drop) carry the face's coarse debris flux and take talus
//!     transport (repose Sc `talus.sc`); (b) threshold failure: an edge steeper than its span's critical slope (and steeper than
//!     the incoming surface was there) fails; the excess moves to the downslope post as talus deposit, which itself rests at the
//!     talus repose slope (cascade); material reaching a fixed post (channel, sea) leaves with the channel (counted).
//! DE3 (in relax.rs) stream power reads the bed K where the drainage area is large enough to strip to bedrock.
//! DE4 guard bound |eta| <= B tanh(eta / B), B = min(amp_frac R, amp_cap_m) (ADR 043 condition 7, declared values).
//!
//! Deterministic: fixed Picard counts, fixed CG caps, fixed traversal orders (ties by index). A pure function of the tile key
//! (it runs inside the tile relaxation; the B7 seam blend joins neighbouring tiles as for the B5 relaxation).

use crate::POST_M;
use serde::Serialize;

pub const DE_NOTE: &str = "lane DE de=1: differential erosion by hardness as a reduced-order PROCESS: a layered weathering-limited hillslope equilibrium (sediment flux q = D S, D = D0/Sc of the exposed beds integrated over each edge's height span, so weak beds are stripped faster; bed contrast faded on soil-mantled gentle ground) forced so that the incoming relaxed surface is steady under the 60 m-section window-mean rock, with deviations longer than screen_m screened (T1b resolves those); then threshold failure at each span's critical slope with the failed rock deposited downslope as talus resting at its repose slope, and talus transport on slopes below bare hard faces; solved on each tile's relaxation domain (fixed channel corridors, sea and domain ring); stream power reads the bed K above a drainage-area threshold. It replaces the v2/v3 bounded remap. Bed strengths come from the beds=2 closure (keyed) and the recorded layer classes; SimTime t2-relaxation-not-history.";

#[derive(Clone, Debug, Serialize)]
pub struct DeParams {
    pub version: u32,
    /// Soil-mantle exposure e = smoothstep((S0 - lo) / width) on the T1b base slope S0 (bug fix B4, logged: the v3 W3 form
    /// S0 / Sc_w called limestone cliffs soil-covered and gentle shale bare; a regolith mantle is limited by the regolith's own
    /// repose slope, about 35 deg, whatever rock is beneath: bare above tan 35 deg = 0.70, mantled below 0.25).
    pub expose_lo: f64,
    pub expose_width: f64,
    /// Picard outer iterations (before and after the talus mantle is marked), damping, CG cap and relative tolerance per step.
    pub outer: [usize; 2],
    pub damp: f64,
    pub cg_max: usize,
    pub cg_tol: f64,
    /// Picard stop: largest update below this (m).
    pub stop_m: f64,
    /// Revision R1 (logged): scale separation. Deviations from the incoming surface longer than screen_m are screened by a
    /// Helmholtz term lambda (z - Z0), lambda = (mean reference edge conductance) / (screen_m / 2 m)^2: T1b already resolved
    /// the landscape above its own resolution with the same rock classes, T2 redistributes below it.
    pub screen_m: f64,
    /// Reference window (m of section) when set (None: the remap's 60 m window table). Experiment key de_window_m.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub window_m: Option<f64>,
    /// DE4 guard: B = min(amp_frac R, amp_cap_m), R = T1b local relief (9 x 9 T1b posts).
    pub amp_frac: f64,
    pub amp_cap_m: f64,
    /// DE2 talus (None = off: no mantle, no failure step).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub talus: Option<TalusParams>,
    /// DE3 stream power by bed (None = off): blend window K -> bed K with drainage area (m2) between a_lo and a_hi.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sp_area_m2: Option<[f64; 2]>,
    /// Include the DE change in the W6 local mass balance (None / false = excluded).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mb: Option<bool>,
    /// DE5 edifice lobes (None = off).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lobes: Option<crate::lobes::LobeParams>,
    /// Lane DE2 F1 (opt-in de_fs=1): formation-scale DE. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fs: Option<FsParams>,
    /// Lane RIB (opt-in de_runoff=1, needs de_lobes=1): rill / gully drainage area weighted by the runoff of the exposed edifice
    /// unit. None = off (params JSON and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub runoff: Option<RunoffParams>,
    /// Lane RIB2 (opt-in de_ringflow=1|2): how T1b's upslope drainage area enters the B5 tile domain. None = the B5 rule (T1b
    /// area / 16 on every post of the second ring, all four sides). 1: only the flux crossing the domain edge inward along the
    /// T1b flow direction (specific area x edge length x inward normal component). 2: concentrated where T1b's own D8 flow lines
    /// (steepest descent on the 32 m T1b lattice) cross into the domain, each line carrying its T1b post's drainage area.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ringflow: Option<u32>,
    /// Lane RILL (opt-in de_b5_chain=1|2): drainage organisation across tile domains from a coarse window-wide B5 pass
    /// (chain.rs). None = off (params JSON and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain: Option<crate::chain::ChainParams>,
    /// Lane DE-TALUS (opt-in de_talus_mantle=1|2 and / or de_talus_runoff=1; need ct_sharp=1 ct_sharp_stage=pre): the 2 m
    /// regime of the pre-stage staircase steers the DE and B5. None = off (params JSON and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub regime: Option<RegimeParams>,
}

/// Lane DE-TALUS (2026-10-17): the 2 m ct_regime of the pre-stage staircase (lane TREAD / SHARP2) as a surface-material switch.
///
/// design 2026-10-17: geology that reads naturally at any site; "geology-derived generation". Lane DIAG-RIPPLE measured the talus
/// egg-crate / ripple rows as DE bed-strength steps on 2 m regime 1 (talus) posts and the straight fall-line grooves as B5
/// rill incision (D8 stream power above the head-law area) on the same talus.
///
/// mantle (de_talus_mantle=1: talus; 2: talus and tread): a talus apron is a debris mantle metres thick; the bed sequence below it
/// is not exposed, so its strength contrasts cannot shape the surface (a mantle hides bed strength contrasts below its
/// thickness; the DE's own soil-mantle rule B4 says the same for regolith). Talus posts take exposure e = 0 in the layered
/// equilibrium (transport at the reference strength under which the incoming T1b talus limb, built at the class talus tangent,
/// is already steady: no bed steps), and in the failure step they fail at the talus repose slope (debris, not rock; never below the incoming slope, as every
/// DE failure). Tread
/// posts (mode 2): the tread is the stripped top of the resistant cap bed, i.e. bedrock, colluvium-veneered at its back by the
/// talus shed from the riser above; it spans one bed, and the cap's hardness is already expressed by T1b's 3-limb profile (the
/// tread exists at the published tread tangent BECAUSE the cap is hard); applying the bed contrast again would steepen the
/// tread by Sc_cap / Sc_w (double counting), so tread posts also take e = 0 in the equilibrium; they keep the rock failure
/// threshold (bedrock). Riser posts (regime 2) are unchanged: held in the equilibrium, failing at their own cliff tangent.
/// runoff (de_talus_runoff=1): coarse open-work talus swallows rain (infiltration capacity far above rain rates), so it makes
/// almost no Horton overland flow and flow arriving from the cliff above sinks into it within tens of metres; B5's drainage area
/// on talus posts is runoff-weighted (own cell area x f_talus) with a transmission loss exp(-2 m / l_talus_m) per talus post
/// crossed (the lane RIB mechanism, which de_runoff applies to lava facies). Rills then start only where enough cliff runoff
/// is concentrated. f_talus and l_talus_m are project assumptions from a recalled qualitative contrast, unverified.
pub const REGIME_NOTE: &str = "lane DE-TALUS (opt-in de_talus_mantle=1|2, de_talus_runoff=1; ct_sharp_stage=pre): 2 m ct_regime 1 (talus) posts are a debris mantle: exposure 0 in the layered equilibrium (the mantle hides the bed strength contrasts below its thickness; the T1b talus limb is kept, no bed steps) and failure at the talus repose slope (never below the incoming slope); mode 2 also gives tread posts (regime 3, the cap bed's stripped top, whose hardness T1b's 3-limb profile already expresses) exposure 0 (no double counting), keeping their rock failure threshold; risers unchanged; de_talus_runoff: talus posts contribute f_talus x cell area to B5's drainage area and pass on exp(-2 m / l_talus_m) of the area they receive (open-work talus infiltrates; f, l project assumptions, unverified).";

#[derive(Clone, Debug, Serialize)]
pub struct RegimeParams {
    pub version: u32,
    /// 0 off, 1 talus mantle, 2 talus mantle + tread exposure 0.
    pub mantle: u32,
    /// None = off; Some([f_talus, l_talus_m]).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub runoff: Option<[f64; 2]>,
    pub note: String,
}

impl RegimeParams {
    pub fn new() -> Self { RegimeParams { version: 1, mantle: 0, runoff: None, note: REGIME_NOTE.into() } }
    /// Default talus runoff factor and transmission-loss length (m).
    pub const RUNOFF_TALUS: [f64; 2] = [0.1, 30.0];
    /// Equilibrium exposure override for a post of 2 m regime `rg` (None: the B4 exposure stands).
    pub fn expose_override(&self, rg: u8) -> Option<f64> {
        match (self.mantle, rg) { (1 | 2, 1) | (2, 3) => Some(0.0), _ => None }
    }
    /// True where the failure step treats the post as talus mantle (repose slope).
    pub fn talus_fail(&self, rg: u8) -> bool { self.mantle >= 1 && rg == 1 }
    /// B5 runoff (own factor, per-post transmission) of a post of 2 m regime `rg` (None: not talus or runoff off).
    pub fn runoff_of(&self, rg: u8) -> Option<(f64, f64)> {
        match (self.runoff, rg) { (Some([f, l]), 1) => Some((f, (-POST_M / l.max(1e-6)).exp())), _ => None }
    }
}

impl Default for RegimeParams { fn default() -> Self { Self::new() } }

/// Lane RIB (2026-10-16, opt-in `de_runoff=1`): substrate-controlled channel initiation on volcanic edifices.
///
/// design 2026-10-15/16: the flank ribs read "regular banding ... really regularly ... more variety"; "geology-derived generation".
/// Lane RIB measured (G crop, one run each): the ribs are rill walls running DOWNSLOPE (line azimuth 48 deg vs summit 38 deg),
/// made by the B5 relaxation (relax=0: 0.53 risers > 30 deg per 100 m across slope; relax on: 1.84; DE adds 0.07). The rills
/// cut through the lobe units without noticing them, because B5's drainage area treats every post as the same runoff source,
/// so the rill spacing is the uniform stream-power / creep spacing everywhere.
///
/// Process: a rill starts where OVERLAND FLOW is enough (Horton runoff; channel initiation thresholds, Montgomery and
/// Dietrich 1988/1992); on young volcanoes infiltration is set by the exposed deposit: clinker / rubbly lava tops and basal
/// breccia swallow rain (lava fields carry few channels for a long time), dense jointed cores shed more, fine-grained pyroclastic
/// density current, lahar and ash deposits shed most (post-eruption ash surfaces gully fast). B5's drainage area is therefore
/// accumulated as runoff-weighted area: each post contributes cell area x f(facies of its exposed unit at the current surface),
/// with f relative to the generic ground B5 was calibrated on. The units, their lateral extent, style and order are the lobe
/// model's (recorded T1b / VE units and the keyed lobe closure scaled by volcano kind); no new draws. The factors are project
/// assumptions from recalled qualitative contrasts (infiltration of clinker >> dense lava >> fine ash), unverified.
pub const RUNOFF_NOTE: &str = "lane RIB de_runoff=1: on edifice posts with lobe columns the B5 rill drainage area is runoff-weighted area: each post contributes its cell area times the runoff factor of the facies exposed at its current surface (rubbly lava top / basal breccia low, dense core moderate, fine pyroclastic / lahar / ash high) and flow crossing permeable units loses area (transmission loss, decay length by facies), so rills start on low-infiltration deposits and lava lobes stay smoother; units and their lateral extent are the lobe model's; factors are project assumptions (qualitative infiltration contrasts), unverified; no new keyed draws.";

#[derive(Clone, Debug, Serialize)]
pub struct RunoffParams {
    pub version: u32,
    /// Runoff factors relative to B5's generic ground: lava rubbly top, lava dense core, lava basal breccia, weak pyroclastic /
    /// lahar / ash, indurated pyroclastic, anything else (non-lava facies of the per-volcano sequence below the lobe stack).
    pub f_lava_top: f64,
    pub f_lava_core: f64,
    pub f_lava_base: f64,
    pub f_pyro: f64,
    pub f_pyro_hard: f64,
    pub f_other: f64,
    /// Transmission-loss decay length (m of flow path) on the same facies: flow crossing a permeable unit sinks into it.
    pub l_lava_top: f64,
    pub l_lava_core: f64,
    pub l_lava_base: f64,
    pub l_pyro: f64,
    pub l_pyro_hard: f64,
    pub l_other: f64,
    pub note: String,
}

impl Default for RunoffParams {
    fn default() -> Self {
        RunoffParams { version: 1, f_lava_top: 0.3, f_lava_core: 0.8, f_lava_base: 0.5, f_pyro: 1.6, f_pyro_hard: 1.0, f_other: 1.0,
                       l_lava_top: 60.0, l_lava_core: 400.0, l_lava_base: 100.0, l_pyro: 1.0e9, l_pyro_hard: 600.0, l_other: 1.0e9, note: RUNOFF_NOTE.into() }
    }
}

impl RunoffParams {
    /// (own runoff factor, per-post transmission exp(-2 m / L)) of a facies.
    pub fn factor(&self, fac: u8) -> (f64, f64) {
        use crate::bedstyle::*;
        let (f, l) = match fac {
            FAC_LAVA_TOP | FAC_BOLE => (self.f_lava_top, self.l_lava_top),
            FAC_LAVA_CORE => (self.f_lava_core, self.l_lava_core),
            FAC_LAVA_BASE => (self.f_lava_base, self.l_lava_base),
            FAC_PYRO => (self.f_pyro, self.l_pyro),
            FAC_PYRO_HARD => (self.f_pyro_hard, self.l_pyro_hard),
            _ => (self.f_other, self.l_other),
        };
        (f, (-crate::POST_M / l.max(1e-6)).exp())
    }
}

/// Lane DE2 F1: formation-scale DE. The reference window (m of section) is clamp(l_h x S0, w_min, w_max): a contrast whose
/// outcrop is narrower than l_h (1.5 x the 64 m T1b Nyquist length) is not expressed by T1b's bicubic surface; the screen is
/// screen_m; the DE change's own masked Gaussian low-pass (sigma mb_sigma_m) is removed (square control: tile means near zero);
/// the failure step's operation cap per free post is fail_cap.
#[derive(Clone, Debug, Serialize)]
pub struct FsParams { pub l_h_m: f64, pub w_min_m: f64, pub w_max_m: f64, pub screen_m: f64, pub mb_sigma_m: f64, pub fail_cap: usize, pub natural_ring: bool }

impl Default for FsParams {
    fn default() -> Self { FsParams { l_h_m: 96.0, w_min_m: 60.0, w_max_m: 200.0, screen_m: 64.0, mb_sigma_m: 64.0, fail_cap: 2000, natural_ring: false } }
}

#[derive(Clone, Debug, Serialize)]
pub struct TalusParams {
    /// Talus repose slope (tan 37 deg = 0.75; recalled range 33-38 deg for coarse angular debris, unverified).
    pub sc: f64,
    /// A bare hard face: effective strength ratio Sc_eff / Sc_w >= m_hard and slope >= face_frac Sc_eff.
    pub m_hard: f64,
    pub face_frac: f64,
    /// Mantle reach below a face post: horizontal walk distance <= reach_mult x the drop below the face post.
    pub reach_mult: f64,
    /// Failure step: operation cap per free post (priority processing, highest post first).
    pub sweeps: usize,
    /// DE2a talus mantle stage (None = off, the default: per-post marking came out speckled on the C6 dev wall).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mantle: Option<bool>,
    /// Lane DE2 F3 (opt-in de_repose=1): a talus deposit rests at the repose slope `sc` whatever the incoming slope was there
    /// (bare rock keeps max(span Sc, incoming slope), so T1b's cliffs stay and the debris slides to their feet).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub repose: Option<bool>,
    /// Lane DECAP (opt-in de_fail_settle=1): the failure step runs past its `sweeps` budget until it settles (None = off: stop
    /// at the budget, unsettled).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub settle: Option<SettleParams>,
}

impl Default for TalusParams {
    fn default() -> Self { TalusParams { sc: 0.75, m_hard: 1.2, face_frac: 0.6, reach_mult: 1.0, sweeps: 400, mantle: None, repose: None, settle: None } }
}

/// Lane DECAP: why the failure step may need more than its operation budget, and what the extension does.
pub const FAIL_SETTLE_RULE: &str = "de_fail_settle_v1 (lane DECAP, opt-in de_fail_settle=1, needs de_talus): the DE failure step (threshold failure + talus deposition, highest post first) is a relaxation whose end state is the physical answer (every bare post at or below its failure threshold, every debris deposit at or below the talus repose slope, volume conserved apart from the counted loss into fixed posts); its operation cap (sweeps x free posts) is a compute budget, not physics. Highest-first processing is a downward height sweep: when the budget runs out, posts below the sweep level have not been processed at all (debris dumped onto them waits as perched piles with fronts up to the incoming cliff slope, and their own rock failure has not happened), so a capped domain is an unfinished state, and which domains cap is input-sensitive (lane DECAP: region R t_41_33 needs 3,654 ops per free post, t_36_35 411 against the budget 400). With the key the same loop continues past the budget, unchanged in rule and order, until no post exceeds its threshold, up to a ceiling of cap_posts operations per free post (reported unsettled if reached); the stats then carry fail_settle (operations beyond the budget, posts still queued at the budget, volume moved and lost into fixed posts beyond it, the whole step's volume residual). A domain that settles within the budget runs exactly the same operations as without the key (byte-identical result). The slowness is the half-excess move rule passing debris down repose aprons in millimetre parcels (lane DECAP: ~75 % of moves < 5 mm), a cost, not a different end state.";

/// Lane DECAP: failure-step settle extension (de_fail_settle=1).
#[derive(Clone, Debug, Serialize)]
pub struct SettleParams {
    /// Ceiling, operations per free post (>= the base budget `sweeps`).
    pub cap_posts: usize,
    pub note: String,
}

impl SettleParams {
    /// Default ceiling (ops per free post): region R's slowest domain needed 3,654 (lane DECAP); ~90 s of failure step at worst.
    pub const CAP_POSTS: usize = 20_000;
    pub fn new() -> Self { SettleParams { cap_posts: Self::CAP_POSTS, note: FAIL_SETTLE_RULE.into() } }
}

impl Default for SettleParams { fn default() -> Self { Self::new() } }

/// Lane DECAP: what the settle extension did in one domain (only written when it ran past the base budget).
#[derive(Clone, Debug, Default, Serialize)]
pub struct SettleStats {
    /// Operations beyond the base budget.
    pub ops_beyond: usize,
    /// Free posts still above their threshold when the base budget ran out (live queue entries).
    pub queued_at_budget: usize,
    /// Volume moved (sum of moves, m3) and lost into fixed posts (m3) beyond the budget.
    pub moved_m3_beyond: f64,
    pub lost_m3_beyond: f64,
    /// Whole failure step: sum of height changes x cell area + loss into fixed posts (m3; 0 up to float rounding: volume conserved).
    pub volume_residual_m3: f64,
}

impl Default for DeParams {
    fn default() -> Self {
        DeParams { version: 1, expose_lo: 0.25, expose_width: 0.45, outer: [20, 0], damp: 0.6, cg_max: 800, cg_tol: 1e-2, stop_m: 0.01, screen_m: 64.0, window_m: None,
                   amp_frac: 0.25, amp_cap_m: 25.0, talus: Some(TalusParams::default()), sp_area_m2: Some([2000.0, 20000.0]), mb: None, lobes: None, fs: None, runoff: None, ringflow: None, chain: None, regime: None }
    }
}

/// Per-domain inputs of the DE solve (n x n posts, row-major, j north).
pub struct DeInput<'a> {
    pub n: usize,
    pub z0: &'a [f64],
    pub fixed: &'a [bool],
    /// Reference (window-mean) critical slope per post, at Z0.
    pub sc_w: &'a [f64],
    /// Soil-mantle exposure e in [0, 1] per post.
    pub expose: &'a [f64],
    /// Guard bound B per post (m); 0 = the post is not changed.
    pub bound: &'a [f64],
    /// Critical slope of the bed exposed at surface height z of post k.
    pub sc_bed: &'a dyn Fn(usize, f64) -> f64,
    /// Stratigraphic (bed) coordinate of surface height z in post k's column (monotone in z).
    pub uw: &'a dyn Fn(usize, f64) -> f64,
    /// Mean of 1/Sc over the bed-coordinate span [u1, u2] of post k's column (u1 = u2: the point value). Bug fix B1 (logged):
    /// an edge on a steep wall spans several beds; its transport must integrate them (series by height), not sample the end posts.
    /// An edge reads the column of its first post (neighbouring columns differ by the dip over 2 m and a smooth pinch).
    pub inv_u: &'a dyn Fn(usize, f64, f64) -> f64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct DeStats {
    pub free_posts: usize,
    pub picard_iters: usize,
    pub cg_iters: usize,
    /// Nonlinear residual after the equilibrium, relative to the reference throughput: |fref - outflow(z) - lambda eta| / |q_ref|
    /// (q_ref = per-post sum of the reference edge flux magnitudes) over free posts.
    pub nl_resid: f64,
    pub face_posts: usize,
    pub mantle_posts: usize,
    pub guard_moved: usize,
    /// Failure step: operations used (fail_sweeps), failed volume (m3), volume left into fixed posts (m3), deposit posts, max
    /// deposit (m).
    pub fail_sweeps: usize,
    /// Posts still unsettled when the sweep cap was reached (0 = the failure step converged).
    pub fail_unsettled: usize,
    /// Lane DECAP (de_fail_settle): the extension past the base budget (absent when the step settled within it or the key is off).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fail_settle: Option<SettleStats>,
    /// Diagnostics: rock failures, talus moves, largest move (m), operations per free post.
    pub fail_last: [f64; 4],
    pub failed_m3: f64,
    pub lost_to_fixed_m3: f64,
    pub deposit_posts: usize,
    pub deposit_max_m: f64,
    pub eta_abs_max_m: f64,
    pub eta_abs_mean_m: f64,
    pub picard_s: f64,
    pub fail_s: f64,
}

/// What the solve returns besides the surface.
pub struct DeOut {
    pub z: Vec<f64>,
    /// Talus deposit thickness (m) from the failure step.
    pub deposit: Vec<f64>,
    /// Talus mantle posts (DE2a).
    pub mantle: Vec<bool>,
    pub stats: DeStats,
}

const D8: [(i64, i64, f64); 8] = [(1, 0, 1.0), (1, 1, std::f64::consts::SQRT_2), (0, 1, 1.0), (-1, 1, std::f64::consts::SQRT_2), (-1, 0, 1.0),
    (-1, -1, std::f64::consts::SQRT_2), (0, -1, 1.0), (1, -1, std::f64::consts::SQRT_2)];

/// Steepest-descent 8-neighbour of post k: (neighbour, slope, distance m).
fn steepest(n: usize, z: &[f64], k: usize) -> Option<(usize, f64, f64)> {
    let (i, j) = ((k % n) as i64, (k / n) as i64);
    let mut best: Option<(usize, f64, f64)> = None;
    for &(dx, dy, w) in &D8 {
        let (a, b) = (i + dx, j + dy);
        if a < 0 || b < 0 || a >= n as i64 || b >= n as i64 { continue; }
        let q = b as usize * n + a as usize;
        let s = (z[k] - z[q]) / (w * POST_M);
        if s > 0.0 && best.map(|bb| s > bb.1).unwrap_or(true) { best = Some((q, s, w * POST_M)); }
    }
    best
}

/// Edge conductances (linear law: C = edge mean 1/Sc): horizontal (i, j)-(i+1, j) at [j (n-1) + i], vertical (i, j)-(i, j+1)
/// at [j n + i].
fn edges(n: usize, m_of: &dyn Fn(usize, usize) -> f64, ch: &mut [f64], cv: &mut [f64]) {
    for j in 0..n { for i in 0..n - 1 { let a = j * n + i; ch[j * (n - 1) + i] = m_of(a, a + 1); } }
    for j in 0..n - 1 { for i in 0..n { let a = j * n + i; cv[j * n + i] = m_of(a, a + n); } }
}

/// Net outflow sum_b C_ab (z_a - z_b) of every post.
fn outflow(n: usize, z: &[f64], ch: &[f64], cv: &[f64], out: &mut [f64]) {
    for v in out.iter_mut() { *v = 0.0; }
    for j in 0..n { for i in 0..n - 1 { let a = j * n + i; let q = ch[j * (n - 1) + i] * (z[a] - z[a + 1]); out[a] += q; out[a + 1] -= q; } }
    for j in 0..n - 1 { for i in 0..n { let a = j * n + i; let q = cv[j * n + i] * (z[a] - z[a + n]); out[a] += q; out[a + n] -= q; } }
}

/// Edge mean 1/Sc over the edge's height span in the first post's column (u_a = its bed coordinate at z_a, cached; the other
/// end mapped into the same column), blended LINEARLY with the reference by the pair's mean exposure (B2); an edge whose lower
/// post is talus mantle takes the talus value.
fn edge_m(inp: &DeInput, z: &[f64], ua: f64, mantle: &[bool], tsc: Option<f64>, a: usize, b: usize) -> f64 {
    let low = if z[a] <= z[b] { a } else { b };
    if let Some(t) = tsc { if mantle[low] { return 1.0 / t; } }
    let w = 0.5 * (1.0 / inp.sc_w[a] + 1.0 / inp.sc_w[b]);
    let e = 0.5 * (inp.expose[a] + inp.expose[b]);
    if e <= 0.0 { return w; }
    let ub = (inp.uw)(a, z[b]);
    (1.0 - e) * w + e * (inp.inv_u)(a, ua, ub)
}

/// Bed coordinate of every post at its current height.
fn u_all(inp: &DeInput, z: &[f64], u: &mut [f64]) { for k in 0..z.len() { u[k] = (inp.uw)(k, z[k]); } }

/// DE2a mantle marking: from every bare hard face post, walk steepest descent while the walk length is at most reach_mult x the
/// drop below the face post; free posts passed that are not themselves hard become mantle.
fn mantle_mark(n: usize, z: &[f64], free: &[bool], sc: &[f64], sc_w: &[f64], t: &TalusParams, mark: &mut [bool]) -> usize {
    for m in mark.iter_mut() { *m = false; }
    let hard = |k: usize| sc[k] >= t.m_hard * sc_w[k];
    let mut faces = 0usize;
    for k in 0..n * n {
        if !free[k] || !hard(k) { continue; }
        let Some((_, s, _)) = steepest(n, z, k) else { continue };
        if s < t.face_frac * sc[k] { continue; }
        faces += 1;
        let (mut cur, mut walked) = (k, 0.0f64);
        for _ in 0..4 * n {
            let Some((q, _, d)) = steepest(n, z, cur) else { break };
            walked += d;
            if !free[q] || walked > t.reach_mult * (z[k] - z[q]) || hard(q) { break; }
            mark[q] = true;
            cur = q;
        }
    }
    faces
}

/// Picard iterations of the screened layered equilibrium from the current z (in place).
#[allow(clippy::too_many_arguments)]
fn picard(inp: &DeInput, p: &DeParams, free: &[bool], fref: &[f64], qref: f64, lam: &[f64], mantle: &[bool], tsc: Option<f64>, iters_n: usize, z: &mut [f64], st: &mut DeStats) {
    let n = inp.n;
    let nn = n * n;
    let (mut ch, mut cv) = (vec![0.0f64; n * (n - 1)], vec![0.0f64; n * (n - 1)]);
    let (mut r, mut pv, mut ap, mut zr, mut dl, mut out, mut diag) = (vec![0.0f64; nn], vec![0.0f64; nn], vec![0.0f64; nn], vec![0.0f64; nn], vec![0.0f64; nn], vec![0.0f64; nn], vec![0.0f64; nn]);
    let mut uu = vec![0.0f64; nn];
    for _it in 0..iters_n {
        u_all(inp, z, &mut uu);
        { let zz: &[f64] = z; let m_of = |a: usize, b: usize| -> f64 { edge_m(inp, zz, uu[a], mantle, tsc, a, b) }; edges(n, &m_of, &mut ch, &mut cv); }
        outflow(n, z, &ch, &cv, &mut out);
        for k in 0..nn { r[k] = if free[k] { fref[k] - out[k] - lam[k] * (z[k] - inp.z0[k]) } else { 0.0 }; }
        for k in 0..nn { diag[k] = lam[k]; }
        for j in 0..n { for i in 0..n - 1 { let a = j * n + i; let c = ch[j * (n - 1) + i]; diag[a] += c; diag[a + 1] += c; } }
        for j in 0..n - 1 { for i in 0..n { let a = j * n + i; let c = cv[j * n + i]; diag[a] += c; diag[a + n] += c; } }
        let matvec = |x: &[f64], y: &mut [f64]| {
            for k in 0..nn { y[k] = if free[k] { diag[k] * x[k] } else { 0.0 }; }
            for j in 0..n { for i in 0..n - 1 { let a = j * n + i; let c = ch[j * (n - 1) + i];
                if free[a] { y[a] -= c * x[a + 1]; } if free[a + 1] { y[a + 1] -= c * x[a]; } } }
            for j in 0..n - 1 { for i in 0..n { let a = j * n + i; let c = cv[j * n + i];
                if free[a] { y[a] -= c * x[a + n]; } if free[a + n] { y[a + n] -= c * x[a]; } } }
        };
        for v in dl.iter_mut() { *v = 0.0; }
        let r0 = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if r0 > 0.0 {
            for k in 0..nn { zr[k] = if free[k] && diag[k] > 0.0 { r[k] / diag[k] } else { 0.0 }; pv[k] = zr[k]; }
            let mut rz: f64 = (0..nn).map(|k| r[k] * zr[k]).sum();
            let mut it = 0usize;
            while it < p.cg_max {
                matvec(&pv, &mut ap);
                let pap: f64 = (0..nn).map(|k| pv[k] * ap[k]).sum();
                if pap <= 0.0 { break; }
                let al = rz / pap;
                for k in 0..nn { dl[k] += al * pv[k]; r[k] -= al * ap[k]; }
                it += 1;
                if r.iter().map(|v| v * v).sum::<f64>().sqrt() <= p.cg_tol * r0 { break; }
                for k in 0..nn { zr[k] = if free[k] && diag[k] > 0.0 { r[k] / diag[k] } else { 0.0 }; }
                let rz2: f64 = (0..nn).map(|k| r[k] * zr[k]).sum();
                let be = rz2 / rz;
                rz = rz2;
                for k in 0..nn { pv[k] = zr[k] + be * pv[k]; }
            }
            st.cg_iters += it;
        }
        let mut mx = 0.0f64;
        for k in 0..nn { if free[k] { z[k] += p.damp * dl[k]; mx = mx.max((p.damp * dl[k]).abs()); } }
        st.picard_iters += 1;
        // converged: the largest update is below the stop tolerance (deterministic: depends only on the iterates)
        if mx < p.stop_m { break; }
    }
    u_all(inp, z, &mut uu);
    { let zz: &[f64] = z; let m_of = |a: usize, b: usize| -> f64 { edge_m(inp, zz, uu[a], mantle, tsc, a, b) }; edges(n, &m_of, &mut ch, &mut cv); }
    outflow(n, z, &ch, &cv, &mut out);
    let nr: f64 = (0..nn).filter(|&k| free[k]).map(|k| (fref[k] - out[k] - lam[k] * (z[k] - inp.z0[k])).powi(2)).sum::<f64>().sqrt();
    st.nl_resid = if qref > 0.0 { nr / qref } else { 0.0 };
}

/// Lane TREAD (ct_sharp_stage=pre): protected riser posts of the 2 m staircase and their rock-mass cliff tangent.
pub struct Risers<'a> {
    pub mask: &'a [bool],
    pub tan: &'a [f64],
    /// Lane DE-TALUS (de_talus_mantle): posts whose bare surface is talus mantle (2 m regime 1) and fails at the talus repose
    /// slope in the failure step (None: rock rule as before).
    pub talus: Option<&'a [bool]>,
    /// Lane CARRY (ct_sharp_carry_strike=1): the carried change is low-passed along the riser contour (None: per-post plane fit).
    pub strike: Option<&'a Strike>,
}

/// Lane TREAD: carry the held riser posts with their surroundings: each riser post (mask) takes the change (z - before) of
/// the least-squares plane through the non-riser posts in `src` within a (2 r + 1)^2 window (their mean when the fit is
/// singular; no change when there are none), so a riser keeps its own shape
/// but moves with the ground it stands in (a one-post riser line held while its neighbours move would stand out as pillars
/// and pits). Reads only the non-riser changes (order-independent).
pub fn carry_risers(n: usize, z: &mut [f64], before: &[f64], mask: &[bool], src: &[bool], r: i64) {
    let d: Vec<f64> = (0..n * n).map(|k| if src[k] && !mask[k] { z[k] - before[k] } else { 0.0 }).collect();
    for k in 0..n * n {
        if !mask[k] { continue; }
        let (i, j) = ((k % n) as i64, (k / n) as i64);
        // least-squares plane d = c0 + cx dx + cy dy through the window's non-riser changes, evaluated at the post (a change
        // that varies linearly across the window, e.g. on a gully wall, is carried exactly); the mean where the fit is singular
        let (mut s0, mut sx, mut sy, mut sxx, mut syy, mut sxy, mut sd, mut sdx, mut sdy) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for b in (j - r).max(0)..=(j + r).min(n as i64 - 1) { for a in (i - r).max(0)..=(i + r).min(n as i64 - 1) {
            let q = b as usize * n + a as usize;
            if src[q] && !mask[q] {
                let (x, y, v) = ((a - i) as f64, (b - j) as f64, d[q]);
                s0 += 1.0; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; sd += v; sdx += v * x; sdy += v * y;
            }
        } }
        let dz = if s0 == 0.0 { 0.0 } else {
            // normal equations [s0 sx sy; sx sxx sxy; sy sxy syy] c = [sd sdx sdy], solved for c0 by Cramer's rule
            let det = s0 * (sxx * syy - sxy * sxy) - sx * (sx * syy - sxy * sy) + sy * (sx * sxy - sxx * sy);
            if s0 >= 4.0 && det.abs() > 1e-6 * s0 * s0 * s0 {
                let d0 = sd * (sxx * syy - sxy * sxy) - sx * (sdx * syy - sxy * sdy) + sy * (sdx * sxy - sxx * sdy);
                d0 / det
            } else { sd / s0 }
        };
        z[k] = before[k] + dz;
    }
}

/// Lane CARRY (opt-in ct_sharp_carry_strike=1): the riser posts' contour (strike) directions and the along-strike Gaussian
/// sigma (posts) for `carry_risers_strike`.
pub struct Strike {
    /// Unit strike vector per post (x, y in post units; (0, 0) where the surface is flat).
    pub sx: Vec<f64>,
    pub sy: Vec<f64>,
    pub sigma_posts: f64,
}

impl Strike {
    /// Strike directions of `z` (n x n) smoothed by a Gaussian of `smooth_posts` (posts), central differences.
    pub fn of_surface(n: usize, z: &[f64], smooth_posts: f64, sigma_posts: f64) -> Strike {
        let zs = crate::v3::masked_lowpass(n, z, &vec![true; n * n], smooth_posts);
        let (mut sx, mut sy) = (vec![0.0f64; n * n], vec![0.0f64; n * n]);
        for j in 0..n { for i in 0..n {
            let (a, b) = (if i > 0 { i - 1 } else { i }, if i + 1 < n { i + 1 } else { i });
            let (c, d) = (if j > 0 { j - 1 } else { j }, if j + 1 < n { j + 1 } else { j });
            let gx = if b > a { (zs[j * n + b] - zs[j * n + a]) / (b - a) as f64 } else { 0.0 };
            let gy = if d > c { (zs[d * n + i] - zs[c * n + i]) / (d - c) as f64 } else { 0.0 };
            let g = gx.hypot(gy);
            if g > 1e-9 { sx[j * n + i] = -gy / g; sy[j * n + i] = gx / g; }
        } }
        Strike { sx, sy, sigma_posts }
    }
}

/// Lane CARRY: Gaussian low-pass of a field along the contour through (x0, y0) (post units): the path is traced both ways in
/// 1-post steps, each step along the local strike `dir(x, y)` (kept on the same side as the previous step), sampling `val(x, y)`;
/// a side stops where `dir` or `val` gives None (the end of a connected band, or no contour). Weights exp(-t^2 / 2 sigma^2),
/// t <= ceil(2.5 sigma). The path never steps across the fall line, so a field that varies only down the fall line (on a plane)
/// is returned unchanged; a field that is linear along the contour is kept wherever both sides reach the full radius.
pub fn contour_lowpass(x0: f64, y0: f64, sigma: f64, dir: &dyn Fn(f64, f64) -> Option<(f64, f64)>, val: &dyn Fn(f64, f64) -> Option<f64>) -> Option<f64> {
    let v0 = val(x0, y0)?;
    let t0 = dir(x0, y0)?;
    let r = (2.5 * sigma).ceil() as i64;
    let (mut acc, mut ws) = (v0, 1.0f64);
    for s in [1.0f64, -1.0] {
        let (mut x, mut y, mut vx, mut vy) = (x0, y0, s * t0.0, s * t0.1);
        for t in 1..=r {
            x += vx; y += vy;
            let Some(v) = val(x, y) else { break };
            let g = (-0.5 * (t as f64 / sigma).powi(2)).exp();
            acc += g * v; ws += g;
            match dir(x, y) {
                Some((ax, ay)) => { let sg = if ax * vx + ay * vy < 0.0 { -1.0 } else { 1.0 }; vx = sg * ax; vy = sg * ay; }
                None => break,
            }
        }
    }
    Some(acc / ws)
}

/// Lane CARRY (ct_sharp_carry_strike=1): `carry_risers`, then the carried change of the riser posts is low-passed along the
/// riser's contour (`contour_lowpass` over the carried posts only, bilinear with a mask weight >= 0.3, so a gap ends the band;
/// strike from `st`), so neighbouring posts of one riser band move together. Basis: a cliff band is one rock body; the B5 /
/// DE change of the slopes around it translates the band (its fall-line variation, e.g. talus below vs tread above, is kept:
/// the path follows the contour), and the 5 x 5 plane fit's post-to-post scatter (each window holds a different set of
/// talus / tread neighbours that moved metres; lane FLUTE: 52 % of the riser along-strike high-pass) is window sampling, not
/// rock. Reads only the per-post carried change (order-independent).
pub fn carry_risers_strike(n: usize, z: &mut [f64], before: &[f64], mask: &[bool], src: &[bool], r: i64, st: &Strike) {
    carry_risers(n, z, before, mask, src, r);
    strike_lowpass_change(n, z, before, mask, st);
}

/// Lane CARRY / CARRY2: the change (z - before) of the posts in `mask` low-passed along the contour (`contour_lowpass` over the
/// mask posts only, bilinear with a mask weight >= 0.3, so a gap ends the band; strike from `st`); other posts untouched.
pub fn strike_lowpass_change(n: usize, z: &mut [f64], before: &[f64], mask: &[bool], st: &Strike) {
    let d: Vec<f64> = (0..n * n).map(|k| if mask[k] { z[k] - before[k] } else { 0.0 }).collect();
    let val = |x: f64, y: f64| -> Option<f64> {
        if x < 0.0 || y < 0.0 || x > (n - 1) as f64 || y > (n - 1) as f64 { return None; }
        let (i0, j0) = ((x.floor() as usize).min(n - 2), (y.floor() as usize).min(n - 2));
        let (fx, fy) = (x - i0 as f64, y - j0 as f64);
        let (mut w, mut s) = (0.0f64, 0.0f64);
        for (di, dj, wt) in [(0, 0, (1.0 - fx) * (1.0 - fy)), (1, 0, fx * (1.0 - fy)), (0, 1, (1.0 - fx) * fy), (1, 1, fx * fy)] {
            let q = (j0 + dj) * n + i0 + di;
            if mask[q] { w += wt; s += wt * d[q]; }
        }
        if w >= 0.3 { Some(s / w) } else { None }
    };
    let dir = |x: f64, y: f64| -> Option<(f64, f64)> {
        let q = (y.round() as usize).min(n - 1) * n + (x.round() as usize).min(n - 1);
        if st.sx[q] == 0.0 && st.sy[q] == 0.0 { None } else { Some((st.sx[q], st.sy[q])) }
    };
    let out: Vec<(usize, f64)> = (0..n * n).filter(|&k| mask[k]).filter_map(|k| contour_lowpass((k % n) as f64, (k / n) as f64, st.sigma_posts, &dir, &val).map(|v| (k, v))).collect();
    for (k, v) in out { z[k] = before[k] + v; }
}

/// Lane CARRY2 (ct_sharp_carry_fringe=1): the riser fringe weight per post: free posts (not `riser`, not `fixed`) at Chebyshev
/// distance d <= w posts from a riser post get (w + 1 - d) / w (1 next to the riser, fading to 1 / w at d = w), else 0.
pub fn riser_fringe(n: usize, riser: &[bool], fixed: &[bool], w: i64) -> Vec<f64> {
    let mut f = vec![0.0f64; n * n];
    for k in 0..n * n {
        if !riser[k] { continue; }
        let (i, j) = ((k % n) as i64, (k / n) as i64);
        for b in (j - w).max(0)..=(j + w).min(n as i64 - 1) { for a in (i - w).max(0)..=(i + w).min(n as i64 - 1) {
            let q = b as usize * n + a as usize;
            if riser[q] || fixed[q] { continue; }
            let d = (a - i).abs().max((b - j).abs());
            f[q] = f[q].max((w + 1 - d) as f64 / w as f64);
        } }
    }
    f
}

/// Lane CARRY2 (ct_sharp_carry_fringe=1): the fringe posts' change (z - before) blended toward its along-contour low-pass
/// (`strike_lowpass_change` over the fringe posts) with the fringe weight: z = before + (1 - f) d + f lp(d).
pub fn carry_fringe_strike(n: usize, z: &mut [f64], before: &[f64], fw: &[f64], st: &Strike) {
    let mask: Vec<bool> = fw.iter().map(|&v| v > 0.0).collect();
    let mut zl = z.to_vec();
    strike_lowpass_change(n, &mut zl, before, &mask, st);
    for k in 0..n * n { if mask[k] { z[k] = z[k] + fw[k] * (zl[k] - z[k]); } }
}

/// DE solve: equilibrium (DE1), mantle (DE2a), guard (DE4), threshold failure with talus deposition (DE2b).
pub fn de_solve(inp: &DeInput, p: &DeParams) -> DeOut { de_solve_rf(inp, p, None) }

/// `de_solve` with lane TREAD's riser mask (None: `de_solve`, the same arithmetic): riser posts are held in the equilibrium,
/// the mantle marking and the guard; in the failure step they are free, fail at their own cliff tangent (no incoming-slope
/// exemption) and pass debris on (a deposit on rock steeper than the repose slope avalanches, rule B3).
pub fn de_solve_rf(inp: &DeInput, p: &DeParams, rf: Option<&Risers>) -> DeOut {
    let n = inp.n;
    let nn = n * n;
    let free_all: Vec<bool> = (0..nn).map(|k| !inp.fixed[k] && inp.bound[k] > 0.0).collect();
    let free: Vec<bool> = match rf { None => free_all.clone(), Some(r) => (0..nn).map(|k| free_all[k] && !r.mask[k]).collect() };
    let mut st = DeStats { free_posts: free.iter().filter(|&&f| f).count(), ..Default::default() };
    // reference outflow (window-mean rock on Z0) and the R1 screening coefficient
    let (mut ch, mut cv) = (vec![0.0f64; n * (n - 1)], vec![0.0f64; n * (n - 1)]);
    let mref = |a: usize, b: usize| 0.5 * (1.0 / inp.sc_w[a] + 1.0 / inp.sc_w[b]);
    edges(n, &mref, &mut ch, &mut cv);
    let mut fref = vec![0.0f64; nn];
    outflow(n, inp.z0, &ch, &cv, &mut fref);
    let qref = {
        let mut g = vec![0.0f64; nn];
        for j in 0..n { for i in 0..n - 1 { let a = j * n + i; let q = (ch[j * (n - 1) + i] * (inp.z0[a] - inp.z0[a + 1])).abs(); g[a] += q; g[a + 1] += q; } }
        for j in 0..n - 1 { for i in 0..n { let a = j * n + i; let q = (cv[j * n + i] * (inp.z0[a] - inp.z0[a + n])).abs(); g[a] += q; g[a + n] += q; } }
        (0..nn).filter(|&k| free[k]).map(|k| g[k] * g[k]).sum::<f64>().sqrt()
    };
    let mut lam = vec![0.0f64; nn];
    let screen = p.fs.as_ref().map(|f| f.screen_m).unwrap_or(p.screen_m);
    if screen > 0.0 {
        let l2 = (screen / POST_M).powi(2);
        let mut cnt = vec![0.0f64; nn];
        for j in 0..n { for i in 0..n - 1 { let a = j * n + i; let c = ch[j * (n - 1) + i]; lam[a] += c; lam[a + 1] += c; cnt[a] += 1.0; cnt[a + 1] += 1.0; } }
        for j in 0..n - 1 { for i in 0..n { let a = j * n + i; let c = cv[j * n + i]; lam[a] += c; lam[a + n] += c; cnt[a] += 1.0; cnt[a + n] += 1.0; } }
        for k in 0..nn { lam[k] = if cnt[k] > 0.0 { lam[k] / cnt[k] / l2 } else { 0.0 }; }
    }
    let t0 = std::time::Instant::now();
    let mut z = inp.z0.to_vec();
    let mut mantle = vec![false; nn];
    picard(inp, p, &free, &fref, qref, &lam, &mantle, None, p.outer[0], &mut z, &mut st);
    let tsc = p.talus.as_ref().map(|t| t.sc);
    if let Some(t) = p.talus.as_ref().filter(|t| t.mantle == Some(true)) {
        let sc: Vec<f64> = (0..nn).map(|k| {
            let (w, e) = (inp.sc_w[k], inp.expose[k]);
            if e <= 0.0 { w } else { 1.0 / ((1.0 - e) / w + e / (inp.sc_bed)(k, z[k])) }
        }).collect();
        st.face_posts = mantle_mark(n, &z, &free, &sc, inp.sc_w, t, &mut mantle);
        st.mantle_posts = mantle.iter().filter(|&&m| m).count();
        if p.outer[1] > 0 { picard(inp, p, &free, &fref, qref, &lam, &mantle, tsc, p.outer[1], &mut z, &mut st); }
    } else if p.outer[1] > 0 {
        picard(inp, p, &free, &fref, qref, &lam, &mantle, None, p.outer[1], &mut z, &mut st);
    }
    // DE4 guard
    for k in 0..nn {
        if !free[k] { continue; }
        let e = z[k] - inp.z0[k];
        let b = inp.bound[k];
        let g = b * (e / b).tanh();
        if (g - e).abs() > 0.1 * e.abs() { st.guard_moved += 1; }
        z[k] = inp.z0[k] + g;
    }
    // lane TREAD: the held risers move with the equilibrium change around them (free_all: riser posts that can move at all)
    if let Some(r) = rf {
        let mv: Vec<bool> = (0..nn).map(|k| r.mask[k] && free_all[k]).collect();
        match r.strike { Some(st) => carry_risers_strike(n, &mut z, inp.z0, &mv, &free, 2, st), None => carry_risers(n, &mut z, inp.z0, &mv, &free, 2) }
    }
    st.picard_s = t0.elapsed().as_secs_f64();
    let t1 = std::time::Instant::now();
    // DE2b threshold failure with talus deposition
    let mut dep = vec![0.0f64; nn];
    if let Some(t) = &p.talus {
        let thr = |zz: &[f64], k: usize, q: usize| -> f64 { 1.0 / edge_m(inp, zz, (inp.uw)(k, zz[k]), &mantle, None, k, q) };
        match p.fs.as_ref() {
            Some(f) => { let mut t2 = t.clone(); t2.sweeps = f.fail_cap; fail_cascade_rf(n, &mut z, inp.z0, &free_all, &thr, &t2, &mut dep, &mut st, rf); }
            None => fail_cascade_rf(n, &mut z, inp.z0, &free_all, &thr, t, &mut dep, &mut st, rf),
        }
    }
    st.fail_s = t1.elapsed().as_secs_f64();
    st.deposit_posts = dep.iter().filter(|&&v| v > 0.05).count();
    st.deposit_max_m = dep.iter().cloned().fold(0.0, f64::max);
    let (mut mx, mut sum) = (0.0f64, 0.0f64);
    for k in 0..nn { if free[k] { let e = (z[k] - inp.z0[k]).abs(); mx = mx.max(e); sum += e; } }
    st.eta_abs_max_m = mx;
    st.eta_abs_mean_m = if st.free_posts > 0 { sum / st.free_posts as f64 } else { 0.0 };
    DeOut { z, deposit: dep, mantle, stats: st }
}

/// DE2b: threshold failure sweeps (highest post first; ties by index). A post steeper toward its steepest-descent neighbour than
/// the threshold (talus deposit on it: the repose slope; bare rock: `thr_rock(z, k, q)`, the span's critical slope; never below
/// the incoming surface's own slope there) moves the excess: half to the neighbour as talus deposit (all of it leaves when the
/// neighbour is fixed: the channel carries it), limited by the deposit available when only talus may move. Volume conserving
/// apart from the counted loss into fixed posts.
#[allow(clippy::too_many_arguments)]
pub fn fail_cascade(n: usize, z: &mut [f64], z0: &[f64], free: &[bool], thr_rock: &dyn Fn(&[f64], usize, usize) -> f64, t: &TalusParams, dep: &mut [f64], st: &mut DeStats) {
    fail_cascade_rf(n, z, z0, free, thr_rock, t, dep, st, None)
}

/// `fail_cascade` with lane TREAD's riser mask (None: the same arithmetic): bare rock of a riser post fails at the riser's cliff
/// tangent, without the incoming-slope exemption.
#[allow(clippy::too_many_arguments)]
pub fn fail_cascade_rf(n: usize, z: &mut [f64], z0: &[f64], free: &[bool], thr_rock: &dyn Fn(&[f64], usize, usize) -> f64, t: &TalusParams, dep: &mut [f64], st: &mut DeStats,
                       rf: Option<&Risers>) {
    // Priority processing, highest post first (ties by index): a changed post's 3 x 3 neighbourhood is re-queued at once, so
    // debris travels down a whole wall in one pass. Every move lowers the potential energy, so it terminates; a cap of
    // `sweeps` x (free posts) operations bounds the work (reported as unsettled if reached).
    use std::cmp::Ordering;
    #[derive(PartialEq)]
    struct Key(f64, usize);
    impl Eq for Key {}
    impl PartialOrd for Key { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
    impl Ord for Key { fn cmp(&self, o: &Self) -> Ordering { self.0.total_cmp(&o.0).then(o.1.cmp(&self.1)) } }
    let cell = POST_M * POST_M;
    let nfree = free.iter().filter(|&&f| f).count();
    let mut heap: std::collections::BinaryHeap<Key> = (0..n * n).filter(|&k| free[k]).map(|k| Key(z[k], k)).collect();
    let cap = t.sweeps.saturating_mul(nfree.max(1));
    // lane DECAP de_fail_settle: the same loop continues past the budget until settled, up to the ceiling
    let lim = match &t.settle { Some(s) => s.cap_posts.max(t.sweeps).saturating_mul(nfree.max(1)), None => cap };
    let mut ops = 0usize;
    let mut last = [0.0f64; 4];
    let mut sst: Option<SettleStats> = None;
    let (zsum0, lost0) = if t.settle.is_some() { (z.iter().enumerate().filter(|&(k, _)| free[k]).map(|(_, v)| *v).sum::<f64>(), st.lost_to_fixed_m3) } else { (0.0, 0.0) };
    while let Some(Key(zk, k)) = heap.pop() {
        if zk != z[k] { continue; } // stale entry (the post moved since it was queued; a fresh entry exists)
        ops += 1;
        if ops > lim { break; }
        if ops == cap + 1 && t.settle.is_some() {
            let live = 1 + heap.iter().filter(|e| e.0 == z[e.1]).count();
            sst = Some(SettleStats { queued_at_budget: live, lost_m3_beyond: -st.lost_to_fixed_m3, ..Default::default() });
        }
        let Some((q, s, d)) = steepest(n, z, k) else { continue };
        let s0 = (z0[k] - z0[q]) / d;
        let repose = t.repose == Some(true) && dep[k] > 1e-6;
        // lane TREAD: bare rock of a protected riser post (its threshold is its own cliff tangent, whatever the incoming slope)
        let riser_rock = match rf { Some(r) => r.mask[k] && dep[k] <= 1e-6, None => false };
        // never below the incoming surface's own slope there (cheap test first: most posts stop here); DE2 F3: a deposit at repose
        if !repose && !riser_rock && s <= s0 + 1e-9 { continue; }
        // lane DE-TALUS: a talus mantle post (2 m regime 1) is debris: it fails at the talus repose slope (never below the
        // incoming slope, as every DE failure), with or without a deposit on it (the mantle below a deposit is debris too, all of
        // it may move); a de_repose deposit keeps its own rule (avalanche off rock steeper than repose). Logged: dropping the
        // incoming-slope exemption on talus (debris can never stand above repose) left the half-excess cascade unsettled at its
        // operation cap on every crop tile and doubled the talus post-scale roughness; not used.
        let talus_mantle = rf.and_then(|r| r.talus).map(|m| m[k]).unwrap_or(false);
        let (th, avail) = if riser_rock { (rf.unwrap().tan[k], f64::INFINITY) } else if repose { (t.sc, dep[k]) } else if talus_mantle { (t.sc.max(s0), f64::INFINITY) }
            else if dep[k] > 1e-6 { (t.sc.max(s0), dep[k]) } else { (thr_rock(z, k, q).max(s0), f64::INFINITY) };
        // tolerance: an excess below 2 mm (twice the output quantum) is at rest; moves below 1 mm are not made
        let excess = z[k] - z[q] - th * d;
        if excess <= 2e-3 { continue; }
        // lane DE2 B3: a deposit on rock steeper than the repose slope cannot rest there at all: the whole deposit moves on
        // (avalanche), so debris reaches the cliff foot in a few hops instead of a long half-excess cascade
        let dz = if repose && s0 > t.sc { avail.min(0.5 * (z[k] - z[q])) } else { if free[q] { 0.5 * excess } else { excess }.min(avail) };
        if dz < 1e-3 { continue; }
        if dep[k] > 1e-6 { last[1] += 1.0; } else { last[0] += 1.0; }
        last[2] = last[2].max(dz);
        z[k] -= dz;
        if dep[k] > 0.0 { dep[k] = (dep[k] - dz).max(0.0); }
        st.failed_m3 += dz * cell;
        if let Some(x) = sst.as_mut() { x.moved_m3_beyond += dz * cell; }
        if free[q] { z[q] += dz; dep[q] += dz; } else { st.lost_to_fixed_m3 += dz * cell; }
        // re-queue what can have become unstable: k itself, the posts above k (steeper toward the lowered k) and q (raised)
        let (i, j) = ((k % n) as i64, (k / n) as i64);
        for dy in -1..=1 { for dx in -1..=1 {
            let (a, b) = (i + dx, j + dy);
            if a < 0 || b < 0 || a >= n as i64 || b >= n as i64 { continue; }
            let r = b as usize * n + a as usize;
            if free[r] && (r == k || z[r] > z[k]) { heap.push(Key(z[r], r)); }
        } }
        if free[q] { heap.push(Key(z[q], q)); }
    }
    st.fail_sweeps = ops;
    st.fail_unsettled = if ops > lim { heap.len() } else { 0 };
    if let Some(mut x) = sst {
        x.ops_beyond = ops - cap;
        x.lost_m3_beyond += st.lost_to_fixed_m3;
        let zsum1: f64 = z.iter().enumerate().filter(|&(k, _)| free[k]).map(|(_, v)| *v).sum();
        x.volume_residual_m3 = (zsum1 - zsum0) * cell + (st.lost_to_fixed_m3 - lost0);
        st.fail_settle = Some(x);
    }
    last[3] = ops as f64 / nfree.max(1) as f64;
    st.fail_last = last;
}

/// Soil-mantle exposure from the base slope and the window critical slope.
#[inline]
pub fn exposure(s0: f64, _sc_w: f64, p: &DeParams) -> f64 { crate::smoothstep((s0 - p.expose_lo) / p.expose_width) }

#[cfg(test)]
mod tests {
    use super::*;

    fn ring(n: usize) -> Vec<bool> { (0..n * n).map(|k| { let (i, j) = (k % n, k / n); i == 0 || j == 0 || i == n - 1 || j == n - 1 }).collect() }

    #[test]
    fn no_contrast_is_identity() {
        let n = 40;
        let z0: Vec<f64> = (0..n * n).map(|k| { let (i, j) = ((k % n) as f64, (k / n) as f64); 200.0 - 0.6 * i + 3.0 * (0.3 * j).sin() + 0.02 * i * j }).collect();
        let fixed = ring(n);
        let scw = vec![0.9; n * n];
        let ex = vec![1.0; n * n];
        let bd = vec![50.0; n * n];
        let f = |_k: usize, _z: f64| 0.9;
        let g = |_k: usize, _a: f64, _b: f64| 1.0 / 0.9;
        let u = |_k: usize, z: f64| z;
        let o = de_solve(&DeInput { n, z0: &z0, fixed: &fixed, sc_w: &scw, expose: &ex, bound: &bd, sc_bed: &f, uw: &u, inv_u: &g }, &DeParams::default());
        let (z, st) = (o.z, o.stats);
        for k in 0..n * n { assert!((z[k] - z0[k]).abs() < 1e-6, "post {k}: {}", z[k] - z0[k]); }
        assert!(st.eta_abs_max_m < 1e-6);
    }

    #[test]
    fn layered_slope_ratio_follows_critical_slopes() {
        // 1D: a plane falling along i (rows identical), horizontal beds 6 m thick alternating Sc 1.5 / 0.5; window harmonic mean 0.75
        let n = 121;
        let g = 0.3;
        let z0: Vec<f64> = (0..n * n).map(|k| 400.0 - g * POST_M * (k % n) as f64).collect();
        let fixed = ring(n);
        let scw = vec![0.75; n * n];
        let ex = vec![1.0; n * n];
        let bd = vec![100.0; n * n];
        let bed = |z: f64| if ((z / 6.0).floor() as i64).rem_euclid(2) == 0 { 1.5 } else { 0.5 };
        let f = |_k: usize, z: f64| bed(z);
        let u = |_k: usize, z: f64| z;
        let gi = |_k: usize, a: f64, b: f64| -> f64 {
            let (a, b) = if a <= b { (a, b) } else { (b, a) };
            if b - a < 1e-6 { return 1.0 / bed(a); }
            let m = 400; (0..m).map(|q| 1.0 / bed(a + (q as f64 + 0.5) * (b - a) / m as f64)).sum::<f64>() / m as f64
        };
        let mut p = DeParams::default();
        p.talus = None;
        p.outer = [60, 0];
        p.cg_tol = 1e-6;
        p.stop_m = 1e-6;
        p.screen_m = 0.0;
        let o = de_solve(&DeInput { n, z0: &z0, fixed: &fixed, sc_w: &scw, expose: &ex, bound: &bd, sc_bed: &f, uw: &u, inv_u: &gi }, &p);
        let (z, st) = (o.z, o.stats);
        let j = n / 2;
        let (mut sh, mut nh, mut ss, mut ns) = (0.0, 0.0, 0.0, 0.0);
        for i in 20..n - 20 {
            let (a, b) = (z[j * n + i], z[j * n + i + 1]);
            // classify edges by the bed at both ends (skip edges that straddle a contact)
            let (ba, bb) = (bed(a), bed(b));
            if ba != bb { continue; }
            let s = (a - b) / POST_M;
            if ba > 1.0 { sh += s; nh += 1.0; } else { ss += s; ns += 1.0; }
        }
        let ratio = (sh / nh) / (ss / ns);
        assert!((ratio / 3.0 - 1.0).abs() < 0.15, "hard / soft slope ratio {ratio} (want 3), hard {} soft {} ({nh} / {ns} edges), stats {st:?}", sh / nh, ss / ns);
        // the mean slope over the middle is kept (window equilibrium)
        let mean = (z[j * n + 20] - z[j * n + n - 21]) / (POST_M * (n - 41) as f64);
        assert!((mean / g - 1.0).abs() < 0.1, "mean slope {mean}");
        assert!(st.nl_resid < 1e-3, "Picard residual {}", st.nl_resid);
    }

    #[test]
    fn failure_cascade_conserves_and_rests_at_repose() {
        // a 12 m rock pillar (rock threshold 2.0) on flat ground: rock fails to 2.0, the debris spreads at the talus repose 0.75
        let n = 41;
        let z0 = vec![0.0f64; n * n];
        let mut z = z0.clone();
        let c = (n / 2) * n + n / 2;
        z[c] = 12.0;
        let free: Vec<bool> = (0..n * n).map(|k| { let (i, j) = (k % n, k / n); !(i == 0 || j == 0 || i == n - 1 || j == n - 1) }).collect();
        let t = TalusParams::default();
        let mut dep = vec![0.0f64; n * n];
        let mut st = DeStats::default();
        let v0: f64 = z.iter().sum();
        fail_cascade(n, &mut z, &z0, &free, &|_z, _k, _q| 2.0, &t, &mut dep, &mut st);
        let v1: f64 = z.iter().sum();
        assert!((v0 - v1 - st.lost_to_fixed_m3 / (POST_M * POST_M)).abs() < 1e-9, "volume {v0} -> {v1}, lost {}", st.lost_to_fixed_m3);
        assert!(st.lost_to_fixed_m3 == 0.0);
        for k in 0..n * n {
            if !free[k] { continue; }
            if let Some((_, s, _)) = steepest(n, &z, k) {
                let lim = if dep[k] > 1e-6 { t.sc } else { 2.0 };
                assert!(s <= lim * (1.0 + 1e-6) + 1e-6, "post {k}: slope {s} > {lim} (dep {})", dep[k]);
            }
        }
        assert!(z[c] < 12.0 && st.deposit_posts == 0 || dep.iter().filter(|&&d| d > 0.0).count() > 4, "debris spread");
    }

    /// Lane DECAP: a 30 m debris column on flat rock, failure budget far too small: without de_fail_settle the step stops
    /// unsettled with debris standing above the repose slope (the perched-pile artefact); with it the same loop settles: debris
    /// volume conserved (no loss: the pile never reaches the fixed ring), bedrock never lowered, every deposit post at or below
    /// the repose slope (2 mm excess tolerance over >= 2 m); the extension stats are written.
    #[test]
    fn fail_settle_finishes_capped_cascade_conserves_and_rests_at_repose() {
        let n = 41;
        let z0 = vec![0.0f64; n * n];
        let free: Vec<bool> = (0..n * n).map(|k| { let (i, j) = (k % n, k / n); !(i == 0 || j == 0 || i == n - 1 || j == n - 1) }).collect();
        let pile = |z: &mut Vec<f64>, dep: &mut Vec<f64>| { for j in 19..22 { for i in 19..22 { z[j * n + i] = 30.0; dep[j * n + i] = 30.0; } } };
        let run = |settle: Option<SettleParams>| {
            let t = TalusParams { repose: Some(true), sweeps: 3, settle, ..TalusParams::default() };
            let (mut z, mut dep) = (z0.clone(), vec![0.0f64; n * n]);
            pile(&mut z, &mut dep);
            let mut st = DeStats::default();
            fail_cascade(n, &mut z, &z0, &free, &|_z, _k, _q| 2.0, &t, &mut dep, &mut st);
            (z, dep, st, t.sc)
        };
        let over = |z: &[f64], dep: &[f64], sc: f64| (0..n * n).filter(|&k| free[k] && dep[k] > 1e-6 && steepest(n, z, k).map(|(_, s, _)| s > sc + 1e-3).unwrap_or(false)).count();
        let (za, da, sa, sc) = run(None);
        assert!(sa.fail_unsettled > 0 && sa.fail_settle.is_none(), "budget 3 / post must cap: unsettled {}", sa.fail_unsettled);
        assert!(over(&za, &da, sc) > 0, "capped: debris above repose");
        let (zb, db, sb, _) = run(Some(SettleParams::new()));
        assert_eq!(sb.fail_unsettled, 0, "settles");
        let x = sb.fail_settle.as_ref().expect("extension stats");
        assert!(x.ops_beyond > 0 && x.queued_at_budget > 0 && x.moved_m3_beyond > 0.0 && sb.fail_sweeps > 3 * free.iter().filter(|&&f| f).count());
        assert!(x.volume_residual_m3.abs() < 1e-6 && x.lost_m3_beyond == 0.0, "residual {}", x.volume_residual_m3);
        assert_eq!(over(&zb, &db, sc), 0, "every deposit post at or below repose");
        let v0 = 9.0 * 30.0;
        let (vz, vd): (f64, f64) = (zb.iter().sum(), db.iter().sum());
        assert!(sb.lost_to_fixed_m3 == 0.0 && (vz - v0).abs() < 1e-9 && (vd - v0).abs() < 1e-9, "volume {vz} deposit {vd} vs {v0}");
        for k in 0..n * n { assert!(zb[k] - db[k] >= z0[k] - 1e-9, "bedrock lowered at {k}"); }
        let peak = zb.iter().cloned().fold(0.0, f64::max);
        assert!(peak < 15.0 && db.iter().filter(|&&d| d > 1e-3).count() > 40, "spread to a cone: peak {peak}");
    }

    /// Lane DECAP: a cascade that settles within its budget runs the same operations with de_fail_settle on (byte-identical).
    #[test]
    fn fail_settle_identity_when_settled_within_budget() {
        let n = 41;
        let z0: Vec<f64> = (0..n * n).map(|k| 0.3 * (k % n) as f64).collect();
        let free: Vec<bool> = (0..n * n).map(|k| { let (i, j) = (k % n, k / n); !(i == 0 || j == 0 || i == n - 1 || j == n - 1) }).collect();
        let run = |settle: Option<SettleParams>| {
            let t = TalusParams { repose: Some(true), settle, ..TalusParams::default() };
            let mut z = z0.clone();
            let mut dep = vec![0.0f64; n * n];
            for j in 10..14 { for i in 25..29 { z[j * n + i] += 9.0; dep[j * n + i] = 9.0; } }
            z[30 * n + 20] += 14.0; // a rock pillar (threshold 2.0)
            let mut st = DeStats::default();
            fail_cascade(n, &mut z, &z0, &free, &|_z, _k, _q| 2.0, &t, &mut dep, &mut st);
            (z, dep, st)
        };
        let (za, da, sa) = run(None);
        let (zb, db, sb) = run(Some(SettleParams::new()));
        assert_eq!(sa.fail_unsettled, 0);
        assert!(sa.fail_sweeps > 0 && sa.fail_sweeps == sb.fail_sweeps && sb.fail_settle.is_none());
        assert!(za.iter().zip(&zb).all(|(a, b)| a.to_bits() == b.to_bits()) && da.iter().zip(&db).all(|(a, b)| a.to_bits() == b.to_bits()));
        assert_eq!(serde_json::to_string(&sa).unwrap(), serde_json::to_string(&sb).unwrap(), "stats identical");
    }

    /// Lane DECAP: de_fail_settle / de_fail_settle_cap are opt-in (params JSON unchanged when off) and validated.
    #[test]
    fn fail_settle_keys_are_opt_in() {
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap();
        p.set("de", "1").unwrap();
        let before = serde_json::to_string(&p).unwrap();
        assert!(!before.contains("settle"));
        assert!(p.set("de_fail_settle_cap", "5000").is_err(), "needs de_fail_settle");
        p.set("de_fail_settle", "1").unwrap();
        assert_eq!(p.de.as_ref().unwrap().talus.as_ref().unwrap().settle.as_ref().map(|s| s.cap_posts), Some(SettleParams::CAP_POSTS));
        assert!(p.set("de_fail_settle_cap", "100").is_err(), "below the base budget");
        p.set("de_fail_settle_cap", "5000").unwrap();
        assert_eq!(p.de.as_ref().unwrap().talus.as_ref().unwrap().settle.as_ref().map(|s| s.cap_posts), Some(5000));
        assert!(serde_json::to_string(&p).unwrap().contains("de_fail_settle_v1"));
        p.set("de_fail_settle", "0").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), before, "off = params JSON unchanged");
        p.set("de_talus", "0").unwrap();
        assert!(p.set("de_fail_settle", "1").is_err(), "needs the talus step");
    }

    /// Lane TREAD: a held riser line is carried with a linearly varying change around it exactly (plane fit), and keeps its
    /// own relief above / below its neighbours.
    #[test]
    fn carry_risers_follows_linear_change() {
        let n = 15;
        let before: Vec<f64> = (0..n * n).map(|k| if k % n + k / n == 14 { 3.0 } else { 0.0 }).collect();
        let mask: Vec<bool> = (0..n * n).map(|k| k % n + k / n == 14).collect();
        let src = vec![true; n * n];
        let mut z: Vec<f64> = (0..n * n).map(|k| before[k] + 0.3 * (k % n) as f64 - 0.7 * (k / n) as f64 + 1.0).collect();
        carry_risers(n, &mut z, &before, &mask, &src, 2);
        for k in 0..n * n {
            let want = before[k] + 0.3 * (k % n) as f64 - 0.7 * (k / n) as f64 + 1.0;
            assert!((z[k] - want).abs() < 1e-9, "post {k}: {} vs {want}", z[k]);
        }
    }

    /// Lane CARRY: the contour low-pass keeps a fall-line profile and a change linear along the contour, and removes post-to-post
    /// scatter along it; a side stops at the end of the band.
    #[test]
    fn contour_lowpass_keeps_fall_line_smooths_strike() {
        let noise = |x: f64| if (x.round() as i64).rem_euclid(2) == 0 { 0.5 } else { -0.5 };
        let dir = |_x: f64, _y: f64| Some((1.0, 0.0));
        let inside = |x: f64| x >= 0.0 && x <= 40.0;
        // fall line along y: 3 y + 0.1 x (linear along strike) passes unchanged in the interior
        let f = |x: f64, y: f64| if inside(x) { Some(3.0 * y + 0.1 * x) } else { None };
        let v = contour_lowpass(20.0, 5.0, 4.0, &dir, &f).unwrap();
        assert!((v - (15.0 + 2.0)).abs() < 1e-9, "{v}");
        // alternating +-0.5 m scatter along strike is damped to < 0.05 m
        let g = |x: f64, y: f64| if inside(x) { Some(3.0 * y + noise(x)) } else { None };
        let v = contour_lowpass(20.0, 5.0, 4.0, &dir, &g).unwrap();
        assert!((v - 15.0).abs() < 0.05, "{v}");
        // at the band's end only one side contributes (no extrapolation past it): still a mean of real samples
        let v = contour_lowpass(0.0, 5.0, 4.0, &dir, &f).unwrap();
        assert!(v > 15.0 && v < 15.0 + 0.1 * 10.0, "{v}");
        // no contour (flat): None
        assert!(contour_lowpass(20.0, 5.0, 4.0, &|_x, _y| None, &f).is_none());
    }

    /// Lane CARRY: a two-post-wide riser band along x whose window neighbours moved by alternating amounts is carried by a
    /// change consistent along the band (scatter damped vs the per-post plane fit), the two rows keep their different mean
    /// change (fall-line variation: talus below moved more than the tread above), and a linear change is still carried exactly
    /// in the band interior.
    #[test]
    fn riser_noflux_creep_no_exchange() {
        // a held riser column (i = 4, 20 m) in flat ground: with the fixed-height exchange its free neighbours gain mass every
        // creep step; as a zero-flux boundary they do not
        let n = 9;
        let p = crate::DetailParamsV1::v1();
        let riser: Vec<bool> = (0..n * n).map(|k| k % n == 4).collect();
        let z0: Vec<f64> = (0..n * n).map(|k| if riser[k] { 20.0 } else { 0.0 }).collect();
        let litho = |_k: usize, _z: f64| (0.6f64, 0.0f64);
        let extra = vec![0.0; n * n];
        let (mut za, mut zb) = (z0.clone(), z0.clone());
        crate::relax::relax_dx_nf(n, crate::POST_M, &mut za, &riser, None, &extra, &litho, None, None, &p);
        crate::relax::relax_dx_nf(n, crate::POST_M, &mut zb, &riser, Some(&riser), &extra, &litho, None, None, &p);
        assert!(za[4 * n + 3] > 0.01, "fixed-height exchange feeds the neighbour: {}", za[4 * n + 3]);
        for k in 0..n * n { assert_eq!(zb[k], z0[k], "post {k}"); }
        let mut zc = z0.clone();
        crate::relax::relax_dx(n, crate::POST_M, &mut zc, &riser, &extra, &litho, None, None, &p);
        assert_eq!(zc, za, "relax_dx = relax_dx_nf without no-flux posts");
    }

    #[test]
    fn carry_fringe_damps_jog_knobs_keeps_mean() {
        // a riser band along x (rows 20-21); the foot fringe (rows 22-23) got a knob every 5 posts (lattice jogs), the crest
        // fringe (rows 18-19) a uniform change
        let n = 41;
        let riser: Vec<bool> = (0..n * n).map(|k| k / n == 20 || k / n == 21).collect();
        let fixed = vec![false; n * n];
        let fw = riser_fringe(n, &riser, &fixed, 2);
        assert_eq!(fw[22 * n + 10], 1.0); assert_eq!(fw[23 * n + 10], 0.5); assert_eq!(fw[24 * n + 10], 0.0); assert_eq!(fw[19 * n + 10], 1.0);
        assert_eq!(fw[20 * n + 10], 0.0, "riser posts are not fringe");
        let before: Vec<f64> = (0..n * n).map(|k| if k / n <= 20 { 10.0 } else { 0.0 }).collect();
        let st = Strike::of_surface(n, &before, 1.5, 4.0);
        let chg = |k: usize| -> f64 { let (i, j) = (k % n, k / n); if j == 22 || j == 23 { 1.0 + if i % 5 == 0 { 2.0 } else { 0.0 } } else if j < 20 { 0.3 } else { 0.0 } };
        let z1: Vec<f64> = (0..n * n).map(|k| before[k] + chg(k)).collect();
        let mut z = z1.clone();
        carry_fringe_strike(n, &mut z, &before, &fw, &st);
        let row = |z: &[f64], j: usize| -> (f64, f64) {
            let d: Vec<f64> = (10..30).map(|i| z[j * n + i] - before[j * n + i]).collect();
            let m = d.iter().sum::<f64>() / d.len() as f64;
            (m, (d.iter().map(|v| (v - m).powi(2)).sum::<f64>() / d.len() as f64).sqrt())
        };
        let ((m0, s0), (m1, s1)) = (row(&z1, 22), row(&z, 22));
        assert!(s1 < 0.25 * s0, "knob scatter {s0} -> {s1}");
        assert!((m0 - m1).abs() < 0.05, "mean {m0} -> {m1}");
        assert!(row(&z, 23).1 < 0.7 * row(&z1, 23).1, "half weight at 2 posts");
        for i in 10..30 { assert!((z[19 * n + i] - z1[19 * n + i]).abs() < 1e-9, "uniform crest change kept"); }
        for k in 0..n * n { if fw[k] == 0.0 { assert_eq!(z[k], z1[k]); } }
    }

    #[test]
    fn carry_risers_strike_moves_the_band_together() {
        let n = 41;
        let mask: Vec<bool> = (0..n * n).map(|k| k / n == 20 || k / n == 21).collect();
        let src = vec![true; n * n];
        let before: Vec<f64> = (0..n * n).map(|k| if k / n <= 20 { 10.0 } else { 0.0 }).collect();
        let st = Strike::of_surface(n, &before, 1.5, 4.0);
        assert!((st.sx[20 * n + 20].abs() - 1.0).abs() < 1e-9 && st.sy[20 * n + 20].abs() < 1e-9);
        // tread above (j < 20) rose 0.2 m, talus below (j > 21) by 2 m +- 1 m alternating along x (window sampling scatter)
        let chg = |k: usize| -> f64 { let (i, j) = (k % n, k / n); if j < 20 { 0.2 } else { 2.0 + if i % 2 == 0 { 1.0 } else { -1.0 } } };
        let z1: Vec<f64> = (0..n * n).map(|k| before[k] + if mask[k] { 0.0 } else { chg(k) }).collect();
        let (mut za, mut zb) = (z1.clone(), z1.clone());
        carry_risers(n, &mut za, &before, &mask, &src, 2);
        carry_risers_strike(n, &mut zb, &before, &mask, &src, 2, &st);
        let rowstat = |z: &[f64], j: usize| -> (f64, f64) {
            let d: Vec<f64> = (10..31).map(|i| z[j * n + i] - before[j * n + i]).collect();
            let m = d.iter().sum::<f64>() / d.len() as f64;
            (m, (d.iter().map(|v| (v - m).powi(2)).sum::<f64>() / d.len() as f64).sqrt())
        };
        for j in [20, 21] {
            let ((ma, sa), (mb, sb)) = (rowstat(&za, j), rowstat(&zb, j));
            assert!(sa > 0.05 && sb < 0.1 * sa, "row {j}: plane-fit scatter {sa}, strike {sb}");
            assert!((ma - mb).abs() < 0.05, "row {j}: mean {ma} vs {mb}");
        }
        assert!(rowstat(&zb, 21).0 - rowstat(&zb, 20).0 > 0.2, "fall-line variation kept");
        // non-riser posts untouched
        for k in 0..n * n { if !mask[k] { assert_eq!(zb[k], z1[k]); } }
        // linear change: exact in the interior (both sides reach the full radius)
        let lin = |k: usize| 0.3 * (k % n) as f64 - 0.7 * (k / n) as f64 + 1.0;
        let mut zc: Vec<f64> = (0..n * n).map(|k| before[k] + if mask[k] { 0.0 } else { lin(k) }).collect();
        carry_risers_strike(n, &mut zc, &before, &mask, &src, 2, &st);
        for i in 12..29 { for j in [20, 21] { let k = j * n + i; assert!((zc[k] - before[k] - lin(k)).abs() < 1e-9, "post {i},{j}"); } }
    }

    /// Lane TREAD: a riser wall that is part of the incoming surface (so the incoming-slope exemption would keep it) and steeper
    /// than its cliff tangent sheds rock down to that tangent; the debris comes to rest below it at the repose slope (an apron),
    /// volume conserved; a riser at its tangent does not move; no mask = the plain cascade (nothing fails).
    #[test]
    fn riser_mask_sheds_to_cliff_tangent_and_builds_apron() {
        let n = 41;
        // a 20 m step at x = 20 (posts i <= 20 high), one riser column; incoming surface = the same (z0 = z)
        let z0: Vec<f64> = (0..n * n).map(|k| if k % n <= 20 { 20.0 } else { 0.0 }).collect();
        let free: Vec<bool> = (0..n * n).map(|k| { let (i, j) = (k % n, k / n); !(i == 0 || j == 0 || i == n - 1 || j == n - 1) }).collect();
        let mask: Vec<bool> = (0..n * n).map(|k| k % n == 20).collect();
        let t = TalusParams { repose: Some(true), ..TalusParams::default() };
        let run = |tan: f64, rf: bool| {
            let tv = vec![tan; n * n];
            let r = Risers { mask: &mask, tan: &tv, talus: None, strike: None };
            let (mut z, mut dep, mut st) = (z0.clone(), vec![0.0f64; n * n], DeStats::default());
            fail_cascade_rf(n, &mut z, &z0, &free, &|_z, _k, _q| 2.0, &t, &mut dep, &mut st, if rf { Some(&r) } else { None });
            (z, dep, st)
        };
        // no mask: the step is the incoming surface, nothing fails
        let (z, _, st) = run(3.0, false);
        assert!(z == z0 && st.failed_m3 == 0.0);
        // tangent 20 / 2 = 10 or more: the wall stands
        let (z, _, _) = run(10.5, true);
        assert!(z == z0);
        // tangent 3: the riser posts fall to 3 x 2 m above the foot, the debris rests below at the repose slope
        let (z, dep, st) = run(3.0, true);
        let v0: f64 = z0.iter().sum();
        let v1: f64 = z.iter().sum();
        assert!((v0 - v1 - st.lost_to_fixed_m3 / (POST_M * POST_M)).abs() < 1e-6, "volume {v0} -> {v1}");
        let j = n / 2;
        let (top, foot) = (z[j * n + 20], z[j * n + 21]);
        assert!(top < 20.0 && top - foot <= 3.0 * POST_M + 0.01, "riser {top} over foot {foot}");
        assert!(dep[j * n + 21] > 0.5, "apron at the foot: {}", dep[j * n + 21]);
        for i in 21..n - 2 {
            let (a, b) = (z[j * n + i], z[j * n + i + 1]);
            if dep[j * n + i] > 1e-6 { assert!((a - b) / POST_M <= t.sc + 1e-3, "apron at i {i}: slope {}", (a - b) / POST_M); }
        }
    }
}

/// Test helper: the domain ring of an n x n domain.
#[cfg(test)]
fn tests_ring(n: usize) -> Vec<bool> { (0..n * n).map(|k| { let (i, j) = (k % n, k / n); i == 0 || j == 0 || i == n - 1 || j == n - 1 }).collect() }

#[cfg(test)]
mod rib_tests {
    use super::*;
    #[test]
    fn runoff_factors_order_by_infiltration() {
        use crate::bedstyle::*;
        let r = RunoffParams::default();
        let (top, core, pyro) = (r.factor(FAC_LAVA_TOP), r.factor(FAC_LAVA_CORE), r.factor(FAC_PYRO));
        assert!(top.0 < core.0 && core.0 < pyro.0);
        assert!(top.1 < core.1 && core.1 < pyro.1 && pyro.1 <= 1.0 && top.1 > 0.0);
        assert_eq!(r.factor(FAC_NONE), (1.0, (-crate::POST_M / 1.0e9f64).exp()));
    }
    #[test]
    fn ringflow_key_is_opt_in() {
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap();
        assert!(p.set("de_ringflow", "1").is_err(), "needs de=1");
        p.set("de", "1").unwrap();
        let before = serde_json::to_string(&p).unwrap();
        assert!(!before.contains("ringflow"));
        p.set("de_ringflow", "2").unwrap();
        assert_eq!(p.de.as_ref().unwrap().ringflow, Some(2));
        assert!(p.set("de_ringflow", "3").is_err());
        p.set("de_ringflow", "0").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), before);
    }
    #[test]
    fn b5_chain_key_is_opt_in() {
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap();
        assert!(p.set("de_b5_chain", "1").is_err(), "needs de=1");
        p.set("de", "1").unwrap();
        let before = serde_json::to_string(&p).unwrap();
        assert!(!before.contains("chain"));
        assert!(p.set("de_b5_chain_cell", "4").is_err(), "needs de_b5_chain");
        p.set("de_b5_chain", "2").unwrap();
        p.set("de_b5_chain_cell", "8").unwrap();
        assert!(p.set("de_b5_chain_dt", "0").is_err());
        let c = p.de.as_ref().unwrap().chain.clone().unwrap();
        assert_eq!((c.mode, c.cell_posts, c.dt, c.no_loss), (2, 8, 4, None));
        p.set("de_b5_chain_loss", "0").unwrap();
        assert_eq!(p.de.as_ref().unwrap().chain.as_ref().unwrap().no_loss, Some(true));
        assert!(p.set("de_b5_chain", "4").is_err());
        p.set("de_b5_chain", "0").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), before);
    }
    /// Lane RILL2: de_b5_chain_interp / _head / _dump are sub-keys (interp needs mode 2|3); the dump path is never serialized;
    /// b8_taper_after is a standalone v3 key, absent from the params JSON unless set.
    #[test]
    fn talus_regime_keys_are_opt_in() {
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap();
        p.set("de", "1").unwrap();
        let before = serde_json::to_string(&p).unwrap();
        assert!(!before.contains("regime"));
        assert!(p.set("de_talus_f", "0.2").is_err(), "needs de_talus_runoff");
        assert!(p.set("de_talus_mantle", "3").is_err());
        p.set("de_talus_mantle", "1").unwrap();
        assert_eq!(p.de.as_ref().unwrap().regime.as_ref().map(|r| (r.mantle, r.runoff)), Some((1, None)));
        p.set("de_talus_runoff", "1").unwrap();
        p.set("de_talus_l_m", "50").unwrap();
        assert_eq!(p.de.as_ref().unwrap().regime.as_ref().unwrap().runoff, Some([RegimeParams::RUNOFF_TALUS[0], 50.0]));
        p.set("de_talus_mantle", "0").unwrap();
        assert!(p.de.as_ref().unwrap().regime.is_some(), "runoff still on");
        p.set("de_talus_runoff", "0").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), before, "off = params JSON unchanged");
        // the rules per regime: 0 soil, 1 talus, 2 riser, 3 tread
        let mut r = RegimeParams::new();
        r.mantle = 1;
        assert_eq!((0..4).map(|g| r.expose_override(g)).collect::<Vec<_>>(), vec![None, Some(0.0), None, None]);
        r.mantle = 2;
        assert_eq!((0..4).map(|g| r.expose_override(g)).collect::<Vec<_>>(), vec![None, Some(0.0), None, Some(0.0)]);
        assert_eq!((0..4).map(|g| r.talus_fail(g)).collect::<Vec<_>>(), vec![false, true, false, false]);
        assert!(r.runoff_of(1).is_none());
        r.runoff = Some([0.1, 30.0]);
        let (f, t) = r.runoff_of(1).unwrap();
        assert!((f - 0.1).abs() < 1e-12 && (t - (-POST_M / 30.0f64).exp()).abs() < 1e-12);
        assert!(r.runoff_of(0).is_none() && r.runoff_of(2).is_none() && r.runoff_of(3).is_none());
    }

    #[test]
    fn talus_runoff_reduces_b5_incision() {
        // B5 on a 35 deg plane with a little cross-slope roughness: talus runoff (f 0.1, loss over 30 m) leaves far less incision
        let n = 65;
        let mut z0 = vec![0.0f64; n * n];
        for k in 0..n * n { let (i, j) = ((k % n) as f64, (k / n) as f64); z0[k] = 0.7 * POST_M * (n as f64 - i) + 0.3 * ((j * 1.7).sin() + (i * 0.9 + j * 0.4).cos()); }
        let fixed = super::tests_ring(n);
        // flow from a cliff above enters along the top (upslope) column: 3000 m2 per post
        let extra: Vec<f64> = (0..n * n).map(|k| if k % n == 1 { 3000.0 } else { 0.0 }).collect();
        let litho = |_k: usize, _z: f64| (1.0, 1.0);
        let p = crate::DetailParamsV1::v3();
        let mut r = RegimeParams::new();
        r.runoff = Some(RegimeParams::RUNOFF_TALUS);
        let ro = |_k: usize, _z: f64| r.runoff_of(1).unwrap();
        let cut = |runoff: Option<&dyn Fn(usize, f64) -> (f64, f64)>| {
            let mut z = z0.clone();
            crate::relax::relax(n, &mut z, &fixed, &extra, &litho, None, runoff, &p);
            let d: Vec<f64> = (0..n * n).map(|k| z0[k] - z[k]).collect();
            // incision depth relative to the cross-slope neighbours (groove depth), p99 over interior posts
            let mut g: Vec<f64> = (0..n * n).filter(|&k| { let (i, j) = (k % n, k / n); i > 2 && j > 2 && i < n - 3 && j < n - 3 })
                .map(|k| d[k] - 0.5 * (d[k - n] + d[k + n])).collect();
            g.sort_by(|a, b| a.total_cmp(b));
            g[g.len() * 99 / 100]
        };
        let (g0, g1) = (cut(None), cut(Some(&ro)));
        assert!(g0 > 0.05, "plain B5 incises ({g0})");
        assert!(g1 < 0.5 * g0, "talus runoff groove p99 {g1} vs {g0}");
    }

    #[test]
    fn talus_mantle_hides_bed_steps_and_rests_at_repose() {
        // the layered slope of layered_slope_ratio_follows_critical_slopes: exposed, the beds make steps (hard 3x soft slope);
        // with the talus mantle (exposure 0 on every post) the equilibrium keeps the incoming plane exactly (no bed steps)
        let n = 61;
        let g = 0.3;
        let z0: Vec<f64> = (0..n * n).map(|k| 400.0 - g * POST_M * (k % n) as f64).collect();
        let fixed = super::tests_ring(n);
        let scw = vec![0.75; n * n];
        let bd = vec![100.0; n * n];
        let bed = |z: f64| if ((z / 6.0).floor() as i64).rem_euclid(2) == 0 { 1.5 } else { 0.5 };
        let f = |_k: usize, z: f64| bed(z);
        let u = |_k: usize, z: f64| z;
        let gi = |_k: usize, a: f64, b: f64| -> f64 {
            let (a, b) = if a <= b { (a, b) } else { (b, a) };
            if b - a < 1e-6 { return 1.0 / bed(a); }
            let m = 200; (0..m).map(|q| 1.0 / bed(a + (q as f64 + 0.5) * (b - a) / m as f64)).sum::<f64>() / m as f64
        };
        let mut p = DeParams::default();
        p.talus = None;
        p.outer = [40, 0];
        p.screen_m = 0.0;
        let mut r = RegimeParams::new();
        r.mantle = 1;
        let solve = |ex: &[f64]| de_solve(&DeInput { n, z0: &z0, fixed: &fixed, sc_w: &scw, expose: ex, bound: &bd, sc_bed: &f, uw: &u, inv_u: &gi }, &p).z;
        let bare = solve(&vec![1.0; n * n]);
        let ex: Vec<f64> = (0..n * n).map(|_| r.expose_override(1).unwrap_or(1.0)).collect();
        let mantled = solve(&ex);
        let dev = |z: &[f64]| (0..n * n).map(|k| (z[k] - z0[k]).abs()).fold(0.0, f64::max);
        assert!(dev(&bare) > 1.0, "bare layered slope has bed steps ({})", dev(&bare));
        assert!(dev(&mantled) < 1e-9, "talus mantle: no bed steps ({})", dev(&mantled));
        // failure step: a talus-mantle post on a 0.9 slope (rock threshold 2.0) fails to the talus repose 0.75; without the mask
        // (rock) it stays; volume conserved
        let m = 31;
        let zf: Vec<f64> = (0..m * m).map(|k| if (k % m) <= m / 2 { 0.9 * POST_M * ((m / 2) - (k % m)) as f64 } else { 0.0 }).collect();
        let free: Vec<bool> = (0..m * m).map(|k| { let (i, j) = (k % m, k / m); !(i == 0 || j == 0 || i == m - 1 || j == m - 1) }).collect();
        let flat = vec![0.0f64; m * m];
        let mask = vec![false; m * m];
        let tan = vec![0.0f64; m * m];
        let tal = vec![true; m * m];
        // the half-excess cascade needs many operations to spread a 27 m debris wedge: a larger cap than the default
        let t = TalusParams { sweeps: 100_000, ..TalusParams::default() };
        let run = |talus: Option<&[bool]>| {
            let (mut z, mut dep, mut st) = (zf.clone(), vec![0.0f64; m * m], DeStats::default());
            let rf = Risers { mask: &mask, tan: &tan, talus, strike: None };
            fail_cascade_rf(m, &mut z, &flat, &free, &|_z, _k, _q| 2.0, &t, &mut dep, &mut st, Some(&rf));
            (z, st)
        };
        let (zr, sr) = run(None);
        assert_eq!(zr, zf, "rock at 0.9 < 2.0 stays");
        assert_eq!(sr.failed_m3, 0.0);
        let (zt, stt) = run(Some(&tal));
        assert!(stt.failed_m3 > 0.0 && stt.fail_unsettled == 0, "failed {} unsettled {}", stt.failed_m3, stt.fail_unsettled);
        let v0: f64 = zf.iter().sum();
        let v1: f64 = zt.iter().sum();
        assert!((v0 - v1 - stt.lost_to_fixed_m3 / (POST_M * POST_M)).abs() < 1e-6, "volume {v0} -> {v1}");
        let j = m / 2;
        for i in 2..m - 3 {
            let s = (zt[j * m + i] - zt[j * m + i + 1]) / POST_M;
            assert!(s <= t.sc + 2e-3, "row slope {s} at {i} above repose; unsettled {} ops {} row {:?}", stt.fail_unsettled, stt.fail_sweeps, &zt[j * m..j * m + m]);
        }
    }

    #[test]
    fn rill2_keys_are_opt_in() {
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap();
        p.set("de", "1").unwrap();
        let before = serde_json::to_string(&p).unwrap();
        assert!(p.set("de_b5_chain_interp", "3").is_err(), "needs de_b5_chain");
        assert!(p.set("de_b5_chain_head", "1").is_err(), "needs de_b5_chain");
        p.set("de_b5_chain", "1").unwrap();
        assert!(p.set("de_b5_chain_interp", "3").is_err(), "needs mode 2|3");
        p.set("de_b5_chain", "3").unwrap();
        let m3 = serde_json::to_string(&p).unwrap();
        p.set("de_b5_chain_dump", "some/dir").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), m3, "dump is not a param");
        assert!(p.set("de_b5_chain_interp", "2").is_err());
        p.set("de_b5_chain_interp", "4").unwrap();
        p.set("de_b5_chain_head", "1").unwrap();
        let c = p.de.as_ref().unwrap().chain.clone().unwrap();
        assert_eq!((c.interp, c.head_scale), (Some(4), Some(true)));
        p.set("de_b5_chain_interp", "1").unwrap();
        p.set("de_b5_chain_head", "0").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), m3);
        p.set("de_b5_chain", "0").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), before);
        assert!(!before.contains("b8_taper_after"));
        p.set("b8_taper_after", "1").unwrap();
        assert!(serde_json::to_string(&p).unwrap().contains("\"b8_taper_after\":true"));
        p.set("b8_taper_after", "0").unwrap();
        assert_eq!(serde_json::to_string(&p).unwrap(), before);
    }
}
