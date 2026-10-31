//! INCISION-THRESHOLD (lane INCISION-THRESHOLD, 2026-10-20): law id
//! `stream-power-excess-steepness-threshold-v1`, an opt-in incision
//! threshold on top of `stream-power-subgrid-nonlinear-v1`
//! (`crate::relief`). Nothing here runs unless the law is selected.
//!
//! WHY: lane RELIEF-DECAY found post-collision peaks collapse within
//! 20-30 My; lane EARTH-NUMBERS (`unpublished-design-record`,
//! section B) matched that to the no-threshold regime of Baldwin, Whipple
//! and Tucker 2003 (JGR 108(B3), abstract read by that lane): detachment-
//! limited stream power decays post-orogenic relief in 1-10 My, an
//! incision threshold (critical shear stress) extends it ~x20, isostasy x6.
//!
//! THE LAW. Excess stream power (the deterministic threshold form,
//! E = K (tau^a - tau_c^a) written in channel steepness; Howard 1994 WRR
//! 30:2261, Tucker and Slingerland 1997, Whipple and Tucker 1999 JGR
//! 104:17661, all recalled, citations not checked):
//!     E = K_n max(0, k_s^n - k_c^n),   k_s = S A^(m/n) = S A^(1/2),
//! K_n, m = n/2 and the receiver rules exactly as `relief.rs`; k_c
//! (`critical_channel_steepness_m`, theta = 0.5 units, m) is the steepness
//! at which floods no longer detach bedrock. It acts on BOTH levels of the
//! two-elevation subgrid model:
//! - trunk: dc/dt = -K_n (A^m ((c - c_r)/L)^n - k_c^n);
//! - cell mean: dz/dt = -K_n (((z - c)/Lambda)^n - k_c^n). Derivation as in
//!   relief.rs: a quasi-steady self-similar subgrid network eroding at E
//!   under the excess law has the same k_s = ((E/K_n) + k_c^n)^(1/n) on every
//!   subgrid reach, so its relief is R = Lambda k_s and E = K_n((R/Lambda)^n -
//!   k_c^n). The cell-mean residual relief at which erosion stops is
//!   Lambda k_c (6.55 k_c at 70 km cells, x_c = 100 m).
//! Backward Euler with the threshold term held explicit-constant (it is
//! not stiff): x' + F x'^n = x0 + K_n k_c^n dt, lowering = max(0, x0 - x').
//! At x0 with stream power below the threshold the root is >= x0, so the
//! level is unchanged (no threshold "deposition"); with k_c = 0 the update is
//! relief.rs's for n != 1 bit for bit (tested).
//! The `critical-slope-envelope-v1` floor (landslides on over-steep fronts),
//! the receiver cap and the effective-area realisation are unchanged.
//!
//! WHY THE DETERMINISTIC FORM, NOT THE STOCHASTIC-DISCHARGE ONE (Lague et
//! al. 2005, Tucker 2004): at T0 a cell is ~70 km and a step 2 My; genesis
//! erosion uses geometric area (runoff does not enter at HEAD), and the
//! stochastic integral needs discharge variability k, the shear-stress-
//! discharge exponent, width scaling and tau_c: 3-4 constants this grid
//! cannot resolve or test. Its threshold-dominated regime is a power law
//! k_s ~ E^f with f = a(1 - w_s)/(b(1 + k)) ~ 0.5 (DiBiase and Whipple 2011,
//! JGR 116 F04036, eq. 15, text read) which the existing n = 2 already is.
//! What that law lacks is the vanishing rate at FINITE steepness that the
//! thinner (exponential) flood tail produces (Tucker 2004; DiBiase and
//! Whipple 2011 Fig. 7) and that Baldwin et al. 2003 identify as the cause
//! of long-lived relief. The excess form adds exactly that with one
//! constant. RISK (declared): with n = 2 standing in for the stochastic
//! regime, a hard cut-off partly double-counts threshold effects; k_c is
//! kept below the steepness of still-incising San Gabriel channels
//! (k_s 30 m^0.9 at theta 0.45 ~ 70-100 m here; erosion 35 m/Ma) so it bites
//! only below the observed active range.
//!
//! MEASURED (lane INCISION-THRESHOLD, C6 recipe, L1, 400 My, seeds 42,
//! 240907, 7; `unpublished-design-record`): k_c = 0 reproduces
//! C6 byte for byte (seed 42). k_c = 60 m raises the 100-400 My mean of land
//! above 1 km in 3 of 3 seeds (+1.0 to +3.1 points of land) but does NOT
//! lengthen the post-collision collapse of peaks: above 3 km about 60 % of
//! the cell-mean lowering is the critical-slope floor, and the threshold
//! withholds ~3 % of stream-power lowering there (25 % over all land). Not a
//! default; opt-in only.
//!
//! TRUNK-ONLY VARIANT (lane LAKE-HOLD-2, 2026-10-21; law id
//! `stream-power-excess-steepness-threshold-trunk-only-v1`,
//! `IncisionThresholdParams::trunk_only`). MEASURED DEFECT (C7 climate
//! recipe + lake outlet rules 5-7, seed 42, 40-200 My; lane LAKE-HOLD-2
//! metrics/thr_census_s42.txt): the cell-mean closure is held at exactly 0
//! on 72 % of evaluated land area every step (81 % of land at 0-200 m, 87 %
//! at 200-500 m, 61 % at 500-1000 m), because the closure needs Lambda k_c =
//! 6.55 x 60 = 393 m of mean relief above the trunk before it erodes at all:
//! the derivation applies the bedrock-detachment threshold to EVERY reach of
//! the subgrid network down to x_c = 100 m, where it means a threshold slope
//! k_c / x_c = 0.6, and its rate equivalent K_n k_c^2 = 14.4 m/My exceeds
//! the whole denudation rate of low-relief cratons (1-10 m/My, recalled).
//! Inside a 70 km cell most of that network is soil-mantled hillslope and
//! colluvial low-order valley, lowered by creep, soil production and
//! solution, which have no steepness threshold; the detachment threshold
//! (DiBiase and Whipple 2011's tau_c) is a property of bedrock channels.
//! Rule: the threshold acts on the trunk channel only; the cell-mean
//! closure is `relief.rs`'s n-law (k_c = 0) bit for bit. No new constant.
//! MEASURED A/B (seeds 42, 240907, 7; controlled windows 40-120, 40-120,
//! 40-160 My): land > 1 km 8.12/6.85/5.73 -> 7.66/6.41/5.38 % of land, mean
//! land elevation -9/-8/-7 m, erosion +6/+7/+8 %, lakes -10/+2/-9 %. The
//! cell-mean threshold RAISES hypsometry and is not a lake lever; the k_c =
//! 0 closure itself lowers lowlands at only 0.2-1.7 m/My. Opt-in only.
//!
//! ALLUVIAL CHANNELS (lane ALLUVIAL-TRUNK, 2026-10-21; law id
//! `stream-power-threshold-free-alluvial-channels-v1`, its own law-set slot;
//! `IncisionThresholdParams::alluvial_grade_slope`, and the same rule on the
//! lake outlet channel, pf-surface closed_basin.rs rule 10). MEASURED DEFECT
//! (lanes LAKE-HOLD and LAKE-HOLD-2, C7 climate + lake outlet rules 5-7):
//! what still holds giant lakes is the outlet channel's bedrock threshold
//! slope k_c / A^(m/n) times overflow paths of 750-6,000 km (a 6,000 km
//! path at ~5e-5 holds a 3e6 km2 lake ~250 m up, seed 42 at 134 My). Big
//! lowland trunks are alluvial (transport-limited): their bed is their own
//! sediment, so a bedrock-detachment threshold does not hold them up; their
//! slope is set by the load they carry (Howard 1980; Whipple and Tucker
//! 2002; recalled, not checked). SWITCH, from existing constants only: a
//! channel is alluvial where its bedrock threshold slope k_c / A^(m/n) is
//! below the retention law's graded alluvial slope S_a, i.e. above the
//! drainage area A_a = (k_c / S_a)^(n/m) (3.6e11 m2 = 3.6e5 km2 at k_c = 60 m,
//! S_a = 1e-4, m/n = 1/2). Why there: a reach held at its threshold slope
//! below S_a is gentler than the grade at which the model's lowland rivers
//! carry their load (rule 2a of the retention law grades every continental
//! cell to receiver + S_a L), so it aggrades whenever it carries load and its
//! bed is alluvium; where the threshold slope is above S_a the channel stands
//! steeper than any load needs and is a bedrock channel. Rule: above A_a the
//! trunk update runs at k_c = 0 (the threshold term of the excess law is
//! dropped). On lake outlet paths (rule 10) the cell is then graded by the
//! retention law's alluvial grade (rule 2a, unchanged: receiver + S_a L when
//! it carries load), not held at k_c / A^(m/n); the subgrid trunk has no
//! deposition term, so a released trunk is simply a k_c = 0 channel. The
//! cell-mean closure is unchanged (full threshold, or k_c = 0 with
//! `trunk_only`). No new constant. DECLARED APPROXIMATIONS: the switch is
//! the model's crossover, far downstream of Earth's bedrock-alluvial
//! transitions (often 1e2-1e4 km2, recalled), so only the largest rivers
//! change; a sediment-starved alluvial reach cuts with no threshold (it
//! incises its own fill, not bedrock); no extra slope floor at S_a is
//! imposed (on every released reach S_a is steeper than the threshold slope
//! it replaces, so a floor would hold lakes higher); the trunk slope drops
//! across A_a (no trunk deposition), it is not continuous there.
//! MEASURED (lane ALLUVIAL-TRUNK, C7 climate + lake outlet rules 5-7, seeds
//! 42, 240907, 7, 40-200 My, plates identical in every arm): NO LAKE LEVER.
//! Lake area 5.19/4.03/3.24 -> 5.31/4.80/3.24 % of continents (with rule 9:
//! 5.03/3.87/3.14; rule 9 alone 4.90/3.89/2.77), per-step lake water
//! +1.3/-2.0/+1.2 %; erosion +0.4 to +1.0 %; hypsometry unchanged. Released
//! trunks lower their slope ~40 % (lowland median 9.3e-5/5.7e-5/3.1e-5 ->
//! 5.6e-5/3.2e-5/1.8e-5 at 3.6e5-1e6/1e6-3e6/> 3e6 km2, seed 42, inside
//! Earth's ~1e-5 to 1e-4 for big lowland rivers, recalled), but lakes spill
//! over cell means: on the overflow path the released channel cuts deeper
//! and rule 2a re-aggrades the cut at S_a (graded storage +9 to +20 %), so
//! the corridor gets rougher, not lower. Opt-in only.

use crate::relief::{
    implicit_power_root, subgrid_lambda, NonlinearSubgridParams, NonlinearSubgridStep,
    SubgridInputs,
};

/// DIAGNOSTIC ONLY (never run state, never in an output or a digest): where
/// this law's cell-mean lowering came from, summed per thread since the last
/// `take_thread_totals` (as `retain::take_thread_totals`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ThresholdTally {
    /// Cell-mean lowering set by the (thresholded) stream-power closure (m^3).
    pub stream_power_m3: f64,
    /// Cell-mean lowering set by the hillslope floor where it exceeded the
    /// stream-power target (critical-slope-envelope-v1) (m^3).
    pub floor_m3: f64,
    /// Stream-power cell-mean lowering the threshold withheld (m^3): the
    /// k_c = 0 target minus the thresholded one.
    pub threshold_withheld_m3: f64,
    /// The same three terms for cells whose start-of-step mean is above 3 km.
    pub high_stream_power_m3: f64,
    pub high_floor_m3: f64,
    pub high_threshold_withheld_m3: f64,
    /// Land cells with a receiver whose stream-power target is 0 because of
    /// the threshold (held), and land cells with a receiver.
    pub held_cells: u64,
    pub land_cells: u64,
}

impl ThresholdTally {
    fn add(&mut self, o: &ThresholdTally) {
        self.stream_power_m3 += o.stream_power_m3;
        self.floor_m3 += o.floor_m3;
        self.threshold_withheld_m3 += o.threshold_withheld_m3;
        self.high_stream_power_m3 += o.high_stream_power_m3;
        self.high_floor_m3 += o.high_floor_m3;
        self.high_threshold_withheld_m3 += o.high_threshold_withheld_m3;
        self.held_cells += o.held_cells;
        self.land_cells += o.land_cells;
    }
}

thread_local! {
    static THREAD_TOTALS: std::cell::RefCell<ThresholdTally> =
        std::cell::RefCell::new(ThresholdTally::default());
}

/// Diagnostic: return and reset this thread's accumulated tallies.
pub fn take_thread_totals() -> ThresholdTally {
    THREAD_TOTALS.with(|t| std::mem::take(&mut *t.borrow_mut()))
}

/// Law id (the typed law set lives in pf-contracts laws.rs).
pub const INCISION_THRESHOLD_LAW_ID: &str = "stream-power-excess-steepness-threshold-v1";

/// Parameters of `stream-power-excess-steepness-threshold-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IncisionThresholdParams {
    /// k_c: critical channel steepness (m, theta = 0.5 units). Project
    /// assumption awaiting review (lane PREREG.txt: 30 m and 60 m).
    pub critical_channel_steepness_m: f64,
    /// LAKE-HOLD-2 (`stream-power-excess-steepness-threshold-trunk-only-v1`):
    /// the threshold acts on the trunk only and the cell-mean closure runs at
    /// k_c = 0; false = `stream-power-excess-steepness-threshold-v1` exactly.
    pub trunk_only: bool,
    /// ALLUVIAL-TRUNK (`stream-power-threshold-free-alluvial-channels-v1`):
    /// S_a, the retention law's graded alluvial slope; the trunk threshold is
    /// dropped where k_c / A^(m/n) < S_a. None = the threshold everywhere
    /// (the laws above exactly).
    pub alluvial_grade_slope: Option<f64>,
}

/// Law id of the alluvial-channel release (lane ALLUVIAL-TRUNK).
pub const ALLUVIAL_CHANNEL_LAW_ID: &str = "stream-power-threshold-free-alluvial-channels-v1";

/// ALLUVIAL-TRUNK: true when a channel of critical steepness `k_c_m` with
/// drainage area `area_m2` is alluvial under graded slope `s_a`: its bedrock
/// threshold slope k_c / A^(m/n) is below S_a (A above A_a = (k_c /
/// S_a)^(n/m)). False for k_c = 0, A <= 0 or a non-finite area.
pub fn is_alluvial_channel(k_c_m: f64, s_a: f64, area_m2: f64, m_over_n: f64) -> bool {
    k_c_m > 0.0 && area_m2 > 0.0 && area_m2.is_finite() && k_c_m / area_m2.powf(m_over_n) < s_a
}

/// Law id of the trunk-only variant (lane LAKE-HOLD-2).
pub const INCISION_THRESHOLD_TRUNK_ONLY_LAW_ID: &str = "stream-power-excess-steepness-threshold-trunk-only-v1";

impl IncisionThresholdParams {
    pub fn validate(&self) -> Result<(), String> {
        let k = self.critical_channel_steepness_m;
        if !k.is_finite() || k < 0.0 {
            return Err(format!(
                "{INCISION_THRESHOLD_LAW_ID}: critical_channel_steepness_m {k} must be finite and >= 0"
            ));
        }
        if let Some(s) = self.alluvial_grade_slope {
            if !(s.is_finite() && s > 0.0) {
                return Err(format!("{ALLUVIAL_CHANNEL_LAW_ID}: alluvial grade slope {s} must be finite and > 0"));
            }
        }
        Ok(())
    }

    /// Erosion-rate equivalent of the threshold, K_n k_c^n (m/My).
    pub fn threshold_rate_m_per_my(&self, p: &NonlinearSubgridParams) -> f64 {
        p.k_n_per_my() * self.critical_channel_steepness_m.powf(p.slope_exponent_n)
    }
}

/// Lowering-limited implicit update of dx/dt = -k (x^n - x_c^n) over one
/// step: returns x' = min(x0, root of x + f x^n = x0 + g), g = k x_c^n dt.
fn excess_root(x0: f64, f: f64, g: f64, n: f64) -> f64 {
    if !(x0 > 0.0) {
        return x0.max(0.0);
    }
    if g == 0.0 {
        return implicit_power_root(x0, f, n);
    }
    implicit_power_root(x0 + g, f, n).min(x0)
}

/// One step of `stream-power-subgrid-nonlinear-v1` with the excess-steepness
/// threshold. Inputs, outputs and every rule other than the threshold are
/// those of `relief::subgrid_step_with_floor`.
#[allow(clippy::too_many_arguments)]
pub fn subgrid_step_threshold(
    inp: &SubgridInputs<'_>,
    z_mean_m: &[f64],
    relief_m: &[f64],
    p: &NonlinearSubgridParams,
    thr: &IncisionThresholdParams,
    stage_k_per_my: f64,
    stage_m: f64,
    dt_my: f64,
    floor_lowering_m: Option<&[f64]>,
) -> Result<NonlinearSubgridStep, String> {
    subgrid_step_threshold_glacial(inp, z_mean_m, relief_m, p, thr, stage_k_per_my, stage_m, dt_my, floor_lowering_m, None)
}

/// GLACIAL-T0 (lane GLACIAL-T0, 2026-10-21): `subgrid_step_threshold` with
/// the opt-in `glacial-stream-power-blend-v1` erodibility factors
/// (`crate::erosion::glacial`): the trunk and cell-mean K_n are multiplied by
/// their factor F and the threshold term by the fluvial weight. `None` is
/// `subgrid_step_threshold` itself; factors of exactly 1 reproduce it bit for
/// bit (multiplication by 1.0 is exact).
#[allow(clippy::too_many_arguments)]
pub fn subgrid_step_threshold_glacial(
    inp: &SubgridInputs<'_>,
    z_mean_m: &[f64],
    relief_m: &[f64],
    p: &NonlinearSubgridParams,
    thr: &IncisionThresholdParams,
    stage_k_per_my: f64,
    stage_m: f64,
    dt_my: f64,
    floor_lowering_m: Option<&[f64]>,
    glacial: Option<&crate::erosion::glacial::GlacialErosivity>,
) -> Result<NonlinearSubgridStep, String> {
    p.validate()?;
    thr.validate()?;
    let n = z_mean_m.len();
    if let Some(g) = glacial {
        for (what, v) in [
            ("trunk factor", &g.trunk_factor),
            ("trunk threshold weight", &g.trunk_threshold_weight),
            ("mean factor", &g.mean_factor),
            ("mean threshold weight", &g.mean_threshold_weight),
        ] {
            if v.len() != n {
                return Err(format!("glacial {what}: {} != {n}", v.len()));
            }
            if let Some(i) = v.iter().position(|x| !(x.is_finite() && *x >= 0.0)) {
                return Err(format!("glacial {what}[{i}] = {} must be finite and >= 0", v[i]));
            }
        }
    }
    let gf = |v: Option<&Vec<f64>>, i: usize| v.map_or(1.0, |v| v[i]);
    if let Some(f) = floor_lowering_m {
        if f.len() != n {
            return Err(format!("subgrid floor: {} != {n}", f.len()));
        }
        if let Some(i) = f.iter().position(|v| !(v.is_finite() && *v >= 0.0)) {
            return Err(format!("subgrid floor[{i}] = {} must be finite and >= 0", f[i]));
        }
    }
    for (what, len) in [
        ("relief", relief_m.len()),
        ("flow_to", inp.flow_to.len()),
        ("area", inp.contributing_area_m2.len()),
        ("reach", inp.reach_length_m.len()),
        ("wet", inp.wet.len()),
        ("fill", inp.fill_level_m.len()),
        ("cell_size", inp.cell_size_m.len()),
    ] {
        if len != n {
            return Err(format!("subgrid input {what}: {len} != {n}"));
        }
    }
    if !(dt_my > 0.0) || !(stage_m > 0.0) || !(stage_k_per_my >= 0.0) {
        return Err("subgrid step: dt, stage m must be positive and stage K >= 0".into());
    }
    let ne = p.slope_exponent_n;
    let m = p.m();
    let kn = p.k_n_per_my();
    let x_c = p.channel_head_length_m;
    // Threshold lowering-equivalent per step (m): K_n k_c^n dt.
    let g = thr.threshold_rate_m_per_my(p) * dt_my;
    let mut trunk = vec![0.0f64; n];
    let mut target = vec![0.0f64; n];
    let mut eff = vec![0.0f64; n];
    let mut mean_after = z_mean_m.to_vec();
    let mut capped = 0u64;
    let mut raised = 0u64;
    let mut tally = ThresholdTally::default();
    let c0: Vec<f64> = (0..n)
        .map(|i| {
            if inp.wet[i] {
                z_mean_m[i]
            } else {
                z_mean_m[i] - relief_m[i].max(0.0)
            }
        })
        .collect();
    trunk.copy_from_slice(&c0);
    for &i in inp.receiver_first {
        if inp.wet[i] {
            trunk[i] = z_mean_m[i];
            continue;
        }
        let Some(j) = inp.flow_to[i] else {
            continue;
        };
        let len = inp.reach_length_m[i];
        let base_t = if inp.wet[j] {
            inp.fill_level_m[j]
        } else {
            trunk[j]
        };
        if c0[i] > base_t {
            let fc = kn * inp.contributing_area_m2[i].powf(m) * dt_my / len.powf(ne)
                * gf(glacial.map(|x| &x.trunk_factor), i);
            // ALLUVIAL-TRUNK: no threshold on an alluvial trunk.
            let gt = match thr.alluvial_grade_slope {
                Some(s_a) if is_alluvial_channel(thr.critical_channel_steepness_m, s_a, inp.contributing_area_m2[i], m / ne) => 0.0,
                _ => g * gf(glacial.map(|x| &x.trunk_threshold_weight), i),
            };
            trunk[i] = base_t + excess_root(c0[i] - base_t, fc, gt, ne);
        } else {
            trunk[i] = base_t.min(z_mean_m[i]);
            raised += 1;
        }
        let lambda = subgrid_lambda(inp.cell_size_m[i], x_c);
        let y0 = (z_mean_m[i] - trunk[i]).max(0.0);
        let (mut d, d_free) = if lambda == 0.0 {
            (y0, y0) // resolved limit: the mean is the trunk
        } else {
            let fm = kn * dt_my / lambda.powf(ne) * gf(glacial.map(|x| &x.mean_factor), i);
            // LAKE-HOLD-2: no threshold on the cell mean with `trunk_only`
            // (excess_root with g = 0 is implicit_power_root, bit for bit).
            let gm = if thr.trunk_only { 0.0 } else { g * gf(glacial.map(|x| &x.mean_threshold_weight), i) };
            (y0 - excess_root(y0, fm, gm, ne), y0 - implicit_power_root(y0, fm, ne))
        };
        let a = inp.cell_size_m[i] * inp.cell_size_m[i];
        let high = z_mean_m[i] > 3000.0;
        tally.land_cells += 1;
        if d <= 0.0 && d_free > 0.0 {
            tally.held_cells += 1;
        }
        let withheld = (d_free - d).max(0.0) * a;
        let mut from_floor = false;
        if let Some(f) = floor_lowering_m {
            if f[i] > d {
                d = f[i];
                from_floor = true;
            }
        }
        if from_floor {
            tally.floor_m3 += d * a;
        } else {
            tally.stream_power_m3 += d * a;
            tally.threshold_withheld_m3 += withheld;
        }
        if high {
            if from_floor {
                tally.high_floor_m3 += d * a;
            } else {
                tally.high_stream_power_m3 += d * a;
                tally.high_threshold_withheld_m3 += withheld;
            }
        }
        let base_m = if inp.wet[j] {
            inp.fill_level_m[j]
        } else {
            mean_after[j]
        };
        let gap = z_mean_m[i] - base_m;
        if !(gap > 0.0) {
            d = 0.0;
        } else if d >= gap * (1.0 - 1e-9) {
            d = gap * (1.0 - 1e-9);
            capped += 1;
        }
        if d > 0.0 {
            let f = d / (gap - d);
            eff[i] = (f * len / (stage_k_per_my * dt_my)).powf(1.0 / stage_m);
            mean_after[i] = z_mean_m[i] - d;
        }
        target[i] = d;
        if trunk[i] > mean_after[i] {
            trunk[i] = mean_after[i];
        }
    }
    THREAD_TOTALS.with(|t| t.borrow_mut().add(&tally));
    Ok(NonlinearSubgridStep {
        trunk_m: trunk,
        target_lowering_m: target,
        effective_area_m2: eff,
        capped_cells: capped,
        trunk_raised_cells: raised,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::relief::subgrid_step_with_floor;

    fn params(n: f64) -> NonlinearSubgridParams {
        NonlinearSubgridParams {
            k_per_my: 2.0,
            slope_exponent_n: n,
            reference_steepness_m: 500.0,
            channel_head_length_m: 100.0,
        }
    }

    /// A short chain: cells 0..3 drain 0 -> 1 -> 2 -> 3 (ocean).
    fn chain(
        z: &[f64; 4],
        relief: &[f64; 4],
        p: &NonlinearSubgridParams,
        kc: Option<f64>,
        floor: Option<&[f64]>,
    ) -> NonlinearSubgridStep {
        let flow_to = [Some(1usize), Some(2), Some(3), None];
        let area = [5.0e9, 1.0e10, 1.5e10, 5.0e9];
        let reach = [70_000.0, 70_000.0, 70_000.0, 1.0];
        let wet = [false, false, false, true];
        let fill = [z[0], z[1], z[2], 0.0];
        let size = [70_000.0; 4];
        let inp = SubgridInputs {
            receiver_first: &[3, 2, 1, 0],
            flow_to: &flow_to,
            contributing_area_m2: &area,
            reach_length_m: &reach,
            wet: &wet,
            fill_level_m: &fill,
            cell_size_m: &size,
        };
        match kc {
            None => subgrid_step_with_floor(&inp, z, relief, p, 2.0, 0.5, 2.0, floor).unwrap(),
            Some(k) => subgrid_step_threshold(
                &inp,
                z,
                relief,
                p,
                &IncisionThresholdParams {
                    critical_channel_steepness_m: k,
                    trunk_only: false,
                    alluvial_grade_slope: None,
                },
                2.0,
                0.5,
                2.0,
                floor,
            )
            .unwrap(),
        }
    }

    #[test]
    fn zero_threshold_is_the_nonlinear_law_bit_for_bit() {
        let z = [3000.0, 1500.0, 400.0, -4000.0];
        let r = [800.0, 500.0, 150.0, 0.0];
        for n in [1.5, 2.0, 3.0] {
            let a = chain(&z, &r, &params(n), None, None);
            let b = chain(&z, &r, &params(n), Some(0.0), None);
            assert_eq!(a, b, "n {n}");
        }
        let floor = [50.0, 0.0, 0.0, 0.0];
        let a = chain(&z, &r, &params(2.0), None, Some(&floor));
        let b = chain(&z, &r, &params(2.0), Some(0.0), Some(&floor));
        assert_eq!(a, b);
    }

    #[test]
    fn threshold_only_slows_and_stops_below_k_c() {
        let p = params(2.0);
        let z = [3000.0, 1500.0, 400.0, -4000.0];
        let r = [800.0, 500.0, 150.0, 0.0];
        let free = chain(&z, &r, &p, None, None);
        let thr = chain(&z, &r, &p, Some(60.0), None);
        for i in 0..3 {
            assert!(thr.target_lowering_m[i] <= free.target_lowering_m[i] + 1e-12, "cell {i}");
            assert!(thr.trunk_m[i] >= free.trunk_m[i] - 1e-9, "cell {i}");
        }
        // Subgrid relief below Lambda k_c with a graded trunk: no lowering.
        let lam = (700.0f64).ln();
        let low = chain(&[300.0, 200.0, 100.0, -4000.0], &[0.9 * lam * 60.0, 0.0, 0.0, 0.0], &p, Some(60.0), None);
        assert!(low.target_lowering_m[0] < 1e-9, "{}", low.target_lowering_m[0]);
        // Same cell without the threshold does erode.
        let low0 = chain(&[300.0, 200.0, 100.0, -4000.0], &[0.9 * lam * 60.0, 0.0, 0.0, 0.0], &p, None, None);
        assert!(low0.target_lowering_m[0] > 0.1);
        // Steep relief: the threshold changes the rate by about (k_c/k_s)^2.
        let hi0 = chain(&[6000.0, 200.0, 100.0, -4000.0], &[6000.0, 0.0, 0.0, 0.0], &p, None, None);
        let hi = chain(&[6000.0, 200.0, 100.0, -4000.0], &[6000.0, 0.0, 0.0, 0.0], &p, Some(60.0), None);
        let rel = (hi0.target_lowering_m[0] - hi.target_lowering_m[0]) / hi0.target_lowering_m[0];
        assert!(rel > 0.0 && rel < 0.05, "{rel}");
    }

    #[test]
    fn excess_root_solves_the_implicit_equation() {
        for &n in &[1.0, 2.0, 3.0] {
            for &(x0, f, g) in &[(1000.0, 1e-3, 5.0), (50.0, 1e-2, 30.0), (10.0, 1e-4, 0.0)] {
                let x = excess_root(x0, f, g, n);
                assert!(x <= x0 && x >= 0.0);
                if x < x0 {
                    let resid = x + f * x.powf(n) - x0 - g;
                    assert!(resid.abs() < 1e-9 * (x0 + g), "n {n}: {resid}");
                } else {
                    assert!(f * x0.powf(n) <= g * (1.0 + 1e-12), "held only below threshold");
                }
            }
        }
    }

    /// GLACIAL-T0: factors of exactly 1 (no ice) reproduce the threshold
    /// step bit for bit; a glacial factor > 1 on the cell mean only lowers
    /// more, never less.
    #[test]
    fn unit_glacial_factors_are_bit_identical() {
        let p = params(2.0);
        let z = [3000.0, 1500.0, 400.0, -4000.0];
        let r = [800.0, 500.0, 150.0, 0.0];
        let flow_to = [Some(1usize), Some(2), Some(3), None];
        let area = [5.0e9, 1.0e10, 1.5e10, 5.0e9];
        let reach = [70_000.0, 70_000.0, 70_000.0, 1.0];
        let wet = [false, false, false, true];
        let fill = [z[0], z[1], z[2], 0.0];
        let size = [70_000.0; 4];
        let inp = SubgridInputs {
            receiver_first: &[3, 2, 1, 0],
            flow_to: &flow_to,
            contributing_area_m2: &area,
            reach_length_m: &reach,
            wet: &wet,
            fill_level_m: &fill,
            cell_size_m: &size,
        };
        let thr = IncisionThresholdParams { critical_channel_steepness_m: 60.0, trunk_only: false, alluvial_grade_slope: None };
        let ones = crate::erosion::glacial::GlacialErosivity {
            trunk_factor: vec![1.0; 4],
            trunk_threshold_weight: vec![1.0; 4],
            mean_factor: vec![1.0; 4],
            mean_threshold_weight: vec![1.0; 4],
        };
        let a = subgrid_step_threshold(&inp, &z, &r, &p, &thr, 2.0, 0.5, 2.0, None).unwrap();
        let b = subgrid_step_threshold_glacial(&inp, &z, &r, &p, &thr, 2.0, 0.5, 2.0, None, Some(&ones)).unwrap();
        assert_eq!(a, b);
        let mut ice = ones.clone();
        ice.mean_factor[0] = 4.0;
        ice.mean_threshold_weight[0] = 0.5;
        let c = subgrid_step_threshold_glacial(&inp, &z, &r, &p, &thr, 2.0, 0.5, 2.0, None, Some(&ice)).unwrap();
        assert!(c.target_lowering_m[0] > a.target_lowering_m[0]);
        assert_eq!(c.trunk_m, a.trunk_m);
    }

    /// LAKE-HOLD-2 trunk-only variant: the cell-mean target is the k_c = 0
    /// closure on the thresholded trunk (equal to the full law where the
    /// threshold does not bind the mean); the trunk is the thresholded one;
    /// a low-relief cell below Lambda k_c now erodes; k_c = 0 is unchanged.
    #[test]
    fn trunk_only_releases_the_cell_mean_only() {
        let p = params(2.0);
        let flow_to = [Some(1usize), Some(2), Some(3), None];
        let area = [5.0e9, 1.0e10, 1.5e10, 5.0e9];
        let reach = [70_000.0, 70_000.0, 70_000.0, 1.0];
        let wet = [false, false, false, true];
        let size = [70_000.0; 4];
        let run = |z: &[f64; 4], r: &[f64; 4], k: f64, trunk_only: bool| {
            let fill = [z[0], z[1], z[2], 0.0];
            let inp = SubgridInputs { receiver_first: &[3, 2, 1, 0], flow_to: &flow_to, contributing_area_m2: &area, reach_length_m: &reach, wet: &wet, fill_level_m: &fill, cell_size_m: &size };
            subgrid_step_threshold(&inp, z, r, &p, &IncisionThresholdParams { critical_channel_steepness_m: k, trunk_only, alluvial_grade_slope: None }, 2.0, 0.5, 2.0, None).unwrap()
        };
        let lam = (700.0f64).ln();
        let (zl, rl) = ([300.0, 200.0, 100.0, -4000.0], [0.9 * lam * 60.0, 0.0, 0.0, 0.0]);
        let full = run(&zl, &rl, 60.0, false);
        let trunk = run(&zl, &rl, 60.0, true);
        let free = run(&zl, &rl, 0.0, false);
        assert!(full.target_lowering_m[0] < 1e-9);
        assert!(trunk.target_lowering_m[0] > 0.1, "{}", trunk.target_lowering_m[0]);
        for i in 0..3 {
            assert!(trunk.target_lowering_m[i] >= full.target_lowering_m[i] - 1e-12, "cell {i}");
            assert!(trunk.trunk_m[i] >= free.trunk_m[i] - 1e-9, "the trunk keeps its threshold ({i})");
        }
        assert_eq!(run(&zl, &rl, 0.0, true), free, "k_c = 0: trunk_only changes nothing");
        // Steep cell: the mean is far above Lambda k_c; trunk_only lowers it
        // by the k_c = 0 closure on the same (thresholded) trunk.
        let (zh, rh) = ([6000.0, 200.0, 100.0, -4000.0], [6000.0, 0.0, 0.0, 0.0]);
        let a = run(&zh, &rh, 60.0, true);
        let b = run(&zh, &rh, 60.0, false);
        assert!(a.target_lowering_m[0] > b.target_lowering_m[0]);
    }

    /// ALLUVIAL-TRUNK: the trunk threshold is dropped exactly where k_c /
    /// sqrt(A) < S_a; below A_a the trunk is the thresholded one bit for bit;
    /// the cell-mean closure and k_c = 0 are unchanged; validation is loud.
    #[test]
    fn alluvial_trunks_drop_the_threshold_above_the_switch_area_only() {
        let p = params(2.0);
        let (k_c, s_a) = (60.0, 1.0e-4);
        let a_a = (k_c / s_a) * (k_c / s_a); // 3.6e11 m2
        assert!(!is_alluvial_channel(k_c, s_a, 0.99 * a_a, 0.5));
        assert!(is_alluvial_channel(k_c, s_a, 1.01 * a_a, 0.5));
        assert!(!is_alluvial_channel(0.0, s_a, 1.0e13, 0.5), "no threshold, nothing to release");
        // Cells 1 and 2 drain to the ocean (3), cell 0 into cell 1; only cell
        // 1 is above A_a. Gentle lowland trunks (3 m over 70 km): k_s ~ 36 m
        // on cell 1, below k_c, so the thresholded trunk does not move.
        let flow_to = [Some(1usize), Some(3), Some(3), None];
        let area = [1.0e11, 2.0 * a_a, 3.0e11, 5.0e9];
        let reach = [70_000.0, 70_000.0, 70_000.0, 1.0];
        let wet = [false, false, false, true];
        let size = [70_000.0; 4];
        let z = [10.0, 4.0, 4.0, -4000.0];
        let r = [1.0, 1.0, 1.0, 0.0];
        let fill = [z[0], z[1], z[2], 0.0];
        let inp = SubgridInputs { receiver_first: &[3, 2, 1, 0], flow_to: &flow_to, contributing_area_m2: &area, reach_length_m: &reach, wet: &wet, fill_level_m: &fill, cell_size_m: &size };
        let run = |k: f64, s: Option<f64>| {
            subgrid_step_threshold(&inp, &z, &r, &p, &IncisionThresholdParams { critical_channel_steepness_m: k, trunk_only: false, alluvial_grade_slope: s }, 2.0, 0.5, 2.0, None).unwrap()
        };
        let off = run(k_c, None);
        let on = run(k_c, Some(s_a));
        let free = run(0.0, None);
        assert_eq!(off.trunk_m[1], 3.0, "held below the threshold");
        assert!(off.trunk_m[1] > on.trunk_m[1] + 1.0, "{} vs {}", off.trunk_m[1], on.trunk_m[1]);
        assert_eq!(on.trunk_m[1], free.trunk_m[1], "k_c dropped on the alluvial trunk");
        assert_eq!(on.trunk_m[2], off.trunk_m[2], "below A_a: the thresholded trunk");
        assert_eq!(on.target_lowering_m[2], off.target_lowering_m[2]);
        assert!(on.trunk_m[0] <= off.trunk_m[0], "a lower receiver never raises a donor");
        assert_eq!(run(0.0, Some(s_a)), free, "k_c = 0: the release changes nothing");
        assert!(IncisionThresholdParams { critical_channel_steepness_m: k_c, trunk_only: false, alluvial_grade_slope: Some(0.0) }.validate().is_err());
        assert!(IncisionThresholdParams { critical_channel_steepness_m: k_c, trunk_only: false, alluvial_grade_slope: Some(f64::NAN) }.validate().is_err());
    }

    #[test]
    fn validation_is_loud() {
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(IncisionThresholdParams {
                critical_channel_steepness_m: bad,
                trunk_only: false,
                alluvial_grade_slope: None,
            }
            .validate()
            .is_err());
        }
    }
}
