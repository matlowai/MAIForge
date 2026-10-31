//! Lane E-SOLVE-CONVERGE (2026-10-22): opt-in convergent driver for the chart-wide T1b rivernet E solve
//! (`t1b_rivernet::solve_e`). Key `esc=1` (needs `es=1` before it; layer AND windows; windows only read the table, the key
//! keeps the recipes equal). Absent = the es / default solve unchanged, byte for byte (the params ride in
//! `EsParamsV1::conv`, skipped when None, and the draw digest strips `es`).
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon; a construction solve (no SimTime).
//!
//! Diagnosis (C7 sites 1 and 5 under the volcanic-split key lists, `unpublished-design-record`):
//! 1. Path-dependent network: es=1 re-runs LS1 (the level solve) on every evaluation, warm-started from the previous
//!    evaluation's network. LS1 is not a fixed point in its 9 passes (site 5: 16k-18k of 154k land receivers still flip in
//!    the last pass), so the residual is a function of the whole history of tables: cells hover +-10-60 m around zero with
//!    8-14 sign changes while ln E moves 1-5 %, the 20 m stop test is a coin toss, and one changed cell re-draws every
//!    cell's trajectory (LANDSCAPE-AGE-P1: 25 of 33 site-5 E values moved; windows +-850 m p1/p99). The default layer also
//!    re-runs LS1 after the solve, so the published network is not the one the residuals were measured on.
//! 2. Coupling stronger than self-control: with the network held, dM_d / d ln E_c of a neighbour exceeds the cell's own
//!    entry on 10 of 13 site-1 cells (122896: own 31 m, from 123472 713 m), and on site 5 reaches 900 m per unit ln E two to
//!    four cells away (33228 -> 33806 / 33808 / 33230). A per-cell (diagonal) update is not contractive.
//! 3. Rim and weak cells run to bounds: a cell whose own mean barely responds to its own E (chart-rim cells, cells at very
//!    low E, marine cells with a target below the trunk level) runs to an E bound and drags its neighbour the other way
//!    through the shared hats (site 5: rim 32656 to 12751 m/My, interior 32655 to 10.8 m/My at +141 m).
//! 4. Roughness: even with the network held, the site-5 residual is not smooth in E (cell 32652's own mean against its
//!    ln E: quadratic-fit residual rms 12 m, max 25 m; +-0.5 % steps move it +4 / +10 m), so 20 m is near the floor there.
//!
//! e_solve_converge_v1 (this driver, v1 = ls_each 1):
//! - LS1 restarts on every evaluation from one fixed reference network (evaluation 0's LS1 output, which depends on the E0
//!   table only) instead of the previous evaluation's: the residual becomes a deterministic function of the table that
//!   includes the rerouting, and the published network is exactly the last evaluation's (LS1 of the published table; the
//!   layer does not re-run LS1 after the solve).
//! - Up to `pre_steps` per-cell es secant steps on cells whose own |r| > `pre_tol`, then a dense Jacobian (one-sided finite
//!   differences of `fd_step` in ln E, one evaluation per active cell; columns measured when a cell first becomes active
//!   and again when a step achieves < 1/4 of its predicted decrease) and Levenberg-Marquardt steps on the active cells
//!   (|r| > tolerance, not pinned at an E bound with the residual pointing outward) with per-cell trust radii and a
//!   Broyden (Schubert) update after each step. Cells within tolerance are not stepped (deadband), so a cell moves only
//!   when its own residual is out of tolerance.
//! - Held cells: rim cells (es chart_rim_v1) unless `rim_step`; cells whose own entry J_dd < `jmin` or < `jshare` x the
//!   largest other entry of their row (their residual is steered by a neighbour's E). Held cells keep their E, are left
//!   out of the stop test and are reported with their residuals.
//! - Flooding offsets F follow the existing F4 rule after base evaluations only (never after a difference evaluation).
//! - Stop: every counted cell (solved, not rim, not held) within tolerance or bound-limited and F4 met or F-bound-limited;
//!   or `steps_max` / `stall` / `evals_max`, when the best base evaluation (smallest max |r|, then sum r^2) is re-evaluated
//!   with its flooding offsets and published.
//! - ls_each 0 (measured, not recommended): the network is held between `reroutes` LS1 runs (rounds), with a defect
//!   correction and a final LS1 run on the best table (`end_on_reroute`). Each reroute shifted the main site-1 cells by
//!   -35 to -107 m at an unchanged table and did not settle in 6 rounds (final 28 m); holding the evaluation-0 network for
//!   the whole solve converges (site 1: 16 evaluations, 10 m) but publishes a network with 39830 LS1 residual pairs (up to
//!   1745 m) against the solved table.

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;

pub const ESC_RULE: &str = "e_solve_converge_v1";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct EscParamsV1 {
    /// Runtime reference address. Content identity is recorded in lane_es.incremental, never in random keys.
    #[serde(skip)]
    pub reference: Option<String>,
    /// LS1 runs: 1 = only the first evaluation (the network is then held for the whole solve: measured on C7 site 1, the
    /// published E0 network then has 39830 LS1 residual pairs up to 1745 m against the solved table, i.e. channels off the
    /// valley floors); n = n - 1 outer reroutes (v1: 2, LS1 re-run once with the round-0 solved table, then a second solve).
    pub reroutes: u32,
    /// Per-cell secant steps per round before the Jacobian phase (while max |r| > pre_tol).
    pub pre_steps: u32,
    pub pre_tol_m: i64,
    /// Newton (LM) steps per round.
    pub steps_max: u32,
    /// A round stops after this many base evaluations without a new best base (0: off).
    pub stall: u32,
    /// 1 (v1, ls_each_v1): every evaluation runs LS1 starting from the fixed reference network (the LS1 output of evaluation 0,
    /// which depends on the E0 table only), so the residual is a deterministic function of the table that includes the
    /// rerouting, and the published network is LS1 of the published table; rounds, defect correction and the final reroute
    /// are then off. 0: the held-network rounds.
    pub ls_each: u32,
    /// Defect correction weight (permille; 0 off): after a reroute, the next inner solve targets r + w x (the shift that
    /// reroute made at an unchanged table).
    pub defect_permille: i64,
    /// 1: after the last round, LS1 runs once more on the best table of that round and that evaluation is published (the
    /// network is then LS1 of the published table); 0: the held network of the last round is published.
    pub end_on_reroute: u32,
    /// Evaluation budget (all kinds).
    pub evals_max: u32,
    /// Finite-difference step in ln E (permille).
    pub fd_step_permille: i64,
    /// Jacobian estimations (finite-difference sweeps) per solve at most.
    pub fd_max: u32,
    /// Per-cell step clamp in ln E (permille).
    pub step_max_permille: i64,
    /// Levenberg-Marquardt damping (permille of the normal-matrix diagonal).
    pub lm_permille: i64,
    /// 1: cells within tolerance are not stepped.
    pub deadband: u32,
    /// 1: rim cells are stepped like the others; 0: rim cells are held at their table value.
    pub rim_step: u32,
    /// Smallest own entry J_dd (m per unit ln E) for a cell to be stepped after a Jacobian estimation (0: never held by J).
    pub jmin_m: i64,
    /// Smallest own share: a cell is held when J_dd < jshare x the largest |J_dc| of another measured cell c in its row
    /// (its residual is steered by a neighbour's E, its own E is ill-determined), permille (0: off).
    pub jshare_permille: i64,
    /// Jacobian support (Chebyshev radius in T0 cells): 0 = dense (every solved cell, one finite-difference evaluation per
    /// non-held cell); r >= 1 = the (2r+1)^2 neighbourhood, coloured (i mod (2r+1), j mod (2r+1)).
    pub support: u32,
    pub rule: String,
}

impl EscParamsV1 {
    pub fn v1() -> Self {
        EscParamsV1 {
            reference: None,
            reroutes: 3, pre_steps: 6, pre_tol_m: 150, steps_max: 16, stall: 5, end_on_reroute: 1, defect_permille: 1000, ls_each: 1, evals_max: 70, fd_step_permille: 250, fd_max: 3, step_max_permille: 693, lm_permille: 20,
            deadband: 1, rim_step: 0, jmin_m: 20, jshare_permille: 100, support: 0,
            rule: "e_solve_converge_v1: ls_each 1 (v1): LS1 restarts on every evaluation from the fixed reference network (evaluation 0's LS1 output from the input network and the E0 table), so r (the es residual per solved cell) is a deterministic function of the table including the rerouting, and the layer publishes the last evaluation's network (no LS1 re-run after the solve). (a) up to pre_steps es secant_diagonal_v1 steps, each on the non-held cells whose own |r| > pre_tol, while the largest counted |r| > pre_tol; (b) J[d][c] = dr_d / d ln E_c by one-sided differences at fd_step: support 0 (v1) dense, one evaluation per column, columns measured for active cells when first needed and all active columns again when a step achieves < 1/4 of its predicted decrease of sum r^2 over the active set (at most fd_max such refreshes); support r >= 1: the (2r+1)^2 neighbourhood in coloured evaluations (colour (i mod (2r+1), j mod (2r+1))); a cell whose J_dd < jmin or < jshare x the largest other measured entry of its row is then held; (c) LM steps on the active set A (|r| > tolerance when deadband, not pinned at an E bound with r pointing outward, not held): (J_AA^T J_AA + lm diag(J_AA^T J_AA)) dx = -J_AA^T r_A, each dx_c clamped to +-min(step_max, trust_c) (trust halves when the cell's |r| grew out of tolerance, doubles back to step_max otherwise), E clamped to [e_min, e_max]; Broyden (Schubert) update of each row after each step; rim cells held unless rim_step; F (coast) by the existing F4 rule after base evaluations only; stop when every counted cell (solved, not rim, not held) is within tolerance or bound-limited and F4 is met or F-bound-limited; at steps_max, stall (base evaluations without a new best) or evals_max the best base (smallest max |r| over counted non-pinned cells, then sum r^2) is re-evaluated with its flooding offsets and published. ls_each 0 (not recommended): the network is held between reroutes LS1 runs (rounds) with defect correction (defect) and a final LS1 run on the best table (end_on_reroute)".into(),
        }
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.reroutes < 1 || self.steps_max < 1 || self.evals_max < 12 { return Err("esc: reroutes >= 1, steps_max >= 1, evals_max >= 12".into()); }
        if self.fd_step_permille <= 0 || self.step_max_permille <= 0 || self.lm_permille < 0 || self.jmin_m < 0 || self.pre_tol_m < 0 { return Err("esc: fd_step / step_max > 0, lm / jmin / pre_tol >= 0".into()); }
        if self.deadband > 1 || self.rim_step > 1 || self.end_on_reroute > 1 || self.ls_each > 1 { return Err("esc: deadband / rim_step / end_on_reroute are 0 or 1".into()); }
        Ok(())
    }
}

/// Key parsing for bin/pf-region-relief (`esc=1|0`, `esc_<item>=<int>`); `es` must be on.
pub fn parse_key(es: &mut crate::t1b_es::EsParamsV1, k: &str, v: i64) -> Result<(), String> {
    if k == "esc" {
        match v { 1 => { if es.conv.is_none() { es.conv = Some(EscParamsV1::v1()); } } 0 => es.conv = None, _ => return Err(format!("bad value esc={v}")) }
        return Ok(());
    }
    let c = es.conv.as_mut().ok_or_else(|| format!("{k} needs esc=1 given before it"))?;
    match k {
        "esc_reroutes" => c.reroutes = v as u32, "esc_pre" => c.pre_steps = v as u32, "esc_pretol" => c.pre_tol_m = v, "esc_steps" => c.steps_max = v as u32, "esc_stall" => c.stall = v as u32, "esc_end" => c.end_on_reroute = v as u32, "esc_defect" => c.defect_permille = v, "esc_lseach" => c.ls_each = v as u32,
        "esc_evals" => c.evals_max = v as u32, "esc_fd" => c.fd_step_permille = v, "esc_fdmax" => c.fd_max = v as u32, "esc_step" => c.step_max_permille = v,
        "esc_lm" => c.lm_permille = v, "esc_deadband" => c.deadband = v as u32, "esc_rim" => c.rim_step = v as u32, "esc_jmin" => c.jmin_m = v, "esc_jshare" => c.jshare_permille = v, "esc_support" => c.support = v as u32,
        _ => return Err(format!("unknown key {k}")),
    }
    c.validate()
}

/// Diagnosis only (env VZ_ESC_PERTURB=<cell>:<dz m>; stdout notes, no bytes unless set): the solve target of one cell raised
/// by dz (T + dz, so r = M - T - dz), for the locality test. Applies to every solve driver.
pub fn perturb_target() -> Option<(u32, f64)> {
    let v = std::env::var("VZ_ESC_PERTURB").ok()?;
    let (c, d) = v.split_once(':')?;
    Some((c.trim().parse().ok()?, d.trim().parse().ok()?))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum Kind { Reroute, Pre, Fd(usize), Step, Best }

/// What the driver asks the solve loop to do next.
pub enum Act {
    /// Evaluate `table` next; `update_flood`: apply the F4 update from the evaluation just done; `reroute`: run LS1;
    /// `flood`: restore these flooding offsets first (the best-base re-evaluation).
    Next { table: BTreeMap<u32, f64>, update_flood: bool, reroute: bool, flood: Option<BTreeMap<u32, f64>> },
    /// Stop: the table of the evaluation just done is the published one.
    Stop,
}

/// One evaluation's per-cell inputs to the driver: the residual, E, and the secant's num = T - L, den = M - L (the means the
/// es update uses for that cell).
pub struct CellObs { pub cell: u32, pub diff_m: f64, pub e: f64, pub num: f64, pub den: f64 }

#[derive(Clone)]
struct Base { x: Vec<f64>, r: Vec<f64>, den: Vec<f64>, table: BTreeMap<u32, f64> }

pub struct Esc {
    p: EscParamsV1,
    es: crate::t1b_es::EsParamsV1,
    tol: f64,
    lemin: f64,
    lemax: f64,
    cells: Vec<u32>,
    idx: BTreeMap<u32, usize>,
    rim: Vec<bool>,
    /// Jacobian support (indices into cells, including d itself) per row.
    sup: Vec<Vec<usize>>,
    colours: Vec<Vec<usize>>,
    jac: Vec<Vec<f64>>,
    have_col: Vec<bool>,
    stale: bool,
    /// Finite-difference sweeps done, and the current sweep's groups (cells perturbed together per evaluation).
    fd_done: u32,
    fd_sweeps: u32,
    fd_queue: Vec<Vec<usize>>,
    fd_evals: u32,
    held: Vec<bool>,
    trust: Vec<f64>,
    sec: Vec<Option<crate::t1b_es::SecantState>>,
    kind: Kind,
    base: Option<Base>,
    /// Best base so far: (max |r| over counted non-pinned cells, sum r^2, table).
    best: Option<(f64, f64, BTreeMap<u32, f64>, BTreeMap<u32, f64>)>,
    /// Predicted sum r^2 over the active set of the last LM step, and the base value.
    pred: Option<(Vec<usize>, f64, f64)>,
    round: u32,
    pre: u32,
    steps: u32,
    pub evals: u32,
    converged_kkt: bool,
    f4_ok_last: bool,
    hist: Vec<Value>,
    fd_far_max_m: f64,
    fd_dx: Vec<f64>,
    stop_reason: String,
    since_best: u32,
    final_reroute: bool,
    defect: Vec<f64>,
    pre_reroute_r: Option<Vec<f64>>,
}

impl Esc {
    /// `cells`: the solved cells; `rim`: the solved rim cells; `ij`: T0 (i, j); `es`: the es params (secant pre-steps).
    #[allow(clippy::too_many_arguments)]
    pub fn new(p: &EscParamsV1, es: &crate::t1b_es::EsParamsV1, cells: &[u32], rim: &[u32], ij: &dyn Fn(u32) -> (i64, i64), nlon: i64, tol: f64, emin: f64, emax: f64) -> Self {
        let n = cells.len();
        let idx: BTreeMap<u32, usize> = cells.iter().enumerate().map(|(k, &c)| (c, k)).collect();
        let pos: Vec<(i64, i64)> = cells.iter().map(|&c| ij(c)).collect();
        let rad = p.support as i64;
        let near = |a: (i64, i64), b: (i64, i64)| -> bool {
            if rad == 0 { return true; }
            let di = (a.0 - b.0).abs();
            let dj = (a.1 - b.1).rem_euclid(nlon);
            let dj = dj.min(nlon - dj);
            di <= rad && dj <= rad
        };
        let sup: Vec<Vec<usize>> = (0..n).map(|d| (0..n).filter(|&c| near(pos[d], pos[c])).collect()).collect();
        let colours: Vec<Vec<usize>> = if rad == 0 { (0..n).map(|k| vec![k]).collect() } else {
            let m = 2 * rad + 1;
            let mut cmap: BTreeMap<(i64, i64), Vec<usize>> = BTreeMap::new();
            for (k, q) in pos.iter().enumerate() { cmap.entry((q.0.rem_euclid(m), q.1.rem_euclid(m))).or_default().push(k); }
            cmap.into_values().collect()
        };
        let rimv: Vec<bool> = cells.iter().map(|c| rim.contains(c)).collect();
        let smax = p.step_max_permille as f64 / 1000.0;
        Esc { p: p.clone(), es: es.clone(), tol, lemin: emin.ln(), lemax: emax.ln(), cells: cells.to_vec(), idx, held: rimv.iter().map(|&r| r && p.rim_step == 0).collect(), rim: rimv, sup, colours,
            jac: vec![vec![0.0; n]; n], have_col: vec![false; n], stale: false, fd_done: 0, fd_sweeps: 0, fd_queue: Vec::new(), fd_evals: 0, trust: vec![smax; n], sec: vec![None; n], kind: Kind::Reroute,
            base: None, best: None, pred: None, round: 0, pre: 0, steps: 0, evals: 0,
            converged_kkt: false, f4_ok_last: false, hist: Vec::new(), fd_far_max_m: 0.0, fd_dx: vec![0.0; n], stop_reason: String::new(), since_best: 0, final_reroute: false, defect: vec![0.0; n], pre_reroute_r: None }
    }

    /// Does the next evaluation run LS1? (The first evaluation always does, as the default solve.)
    pub fn reroute_next(&self) -> bool { self.kind == Kind::Reroute || self.p.ls_each == 1 }
    /// ls_each: every evaluation runs LS1 from the fixed reference network (evaluation 0's LS1 output).
    pub fn ls_each(&self) -> bool { self.p.ls_each == 1 }
    pub fn kind(&self) -> Kind { self.kind }
    pub fn converged_kkt(&self) -> bool { self.converged_kkt }

    fn pinned(&self, x: f64, r: f64) -> bool { (x <= self.lemin + 1e-9 && r > 0.0) || (x >= self.lemax - 1e-9 && r < 0.0) }
    fn open(&self, k: usize, x: f64, r: f64) -> bool {
        !self.held[k] && (r.abs() > self.tol || self.p.deadband == 0) && !self.pinned(x, r)
    }
    /// KKT-open counted cells (non-rim, not held).
    fn kkt(&self, x: &[f64], r: &[f64]) -> Vec<usize> {
        (0..self.cells.len()).filter(|&k| !self.rim[k] && !self.held[k]).filter(|&k| !(r[k].abs() <= self.tol || self.pinned(x[k], r[k]))).collect()
    }
    /// Merit of a base: (max |r| over counted non-pinned cells, sum of their r^2).
    fn merit(&self, x: &[f64], r: &[f64]) -> (f64, f64) {
        let ks: Vec<usize> = (0..self.cells.len()).filter(|&k| !self.rim[k] && !self.held[k] && !self.pinned(x[k], r[k])).collect();
        (ks.iter().map(|&k| r[k].abs()).fold(0.0, f64::max), ks.iter().map(|&k| r[k] * r[k]).sum())
    }

    fn table_of(&self, b: &BTreeMap<u32, f64>, x: &[f64]) -> BTreeMap<u32, f64> {
        let mut t = b.clone();
        for (k, &c) in self.cells.iter().enumerate() { t.insert(c, x[k].exp()); }
        t
    }

    fn stop(&mut self, done: bool, reason: &str, cur: &BTreeMap<u32, f64>, log: &mut dyn FnMut(String)) -> Act {
        self.converged_kkt = done;
        self.stop_reason = reason.into();
        // Not converged: publish the best base seen (one more held-network evaluation when it is not the current one).
        if !done {
            if let Some((bm, bs, bt, bf)) = self.best.clone() {
                if bt != *cur && self.kind != Kind::Best {
                    log(format!("  lane E-SOLVE-CONVERGE: re-evaluating the best base (max |r| {bm:.1} m, sum r^2 {bs:.0})"));
                    self.kind = Kind::Best;
                    return Act::Next { table: bt, update_flood: false, reroute: false, flood: Some(bf) };
                }
            }
        }
        Act::Stop
    }

    /// After an evaluation of the table this driver asked for: `obs` (one per solved cell, any order), `f4_ok` (F4 met or
    /// F-bound-limited on every solved coastal cell). Returns what to evaluate next, or Stop.
    pub fn after_eval(&mut self, obs: &[CellObs], table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, f4_ok: bool, log: &mut dyn FnMut(String)) -> Act {
        self.evals += 1;
        let n = self.cells.len();
        let mut x = vec![0.0; n];
        let mut r = vec![0.0; n];
        let mut nd = vec![(0.0, 0.0); n];
        for o in obs { if let Some(&k) = self.idx.get(&o.cell) { x[k] = o.e.max(1e-9).ln(); r[k] = o.diff_m; nd[k] = (o.num, o.den); } }
        let kind = self.kind;
        // Defect correction: a reroute evaluation measures the shift LS1 makes at an unchanged table (raw r now minus raw r
        // of the base it rerouted from); the next inner solve targets r + shift, anticipating the next reroute.
        if kind == Kind::Reroute {
            if let (Some(pr), true) = (self.pre_reroute_r.take(), self.p.defect_permille > 0) {
                let w = self.p.defect_permille as f64 / 1000.0;
                self.defect = (0..n).map(|k| w * (r[k] - pr[k])).collect();
                log(format!("  lane E-SOLVE-CONVERGE reroute shift (r after LS1 minus r before, same table) {:?}", (0..n).filter(|&k| !self.held[k] && !self.rim[k]).map(|k| (self.cells[k], (r[k] - pr[k]).round())).collect::<Vec<_>>()));
            }
        }
        // The stop test at a reroute (and the best-base evaluation) uses the raw residual; the inner solve the shifted one.
        let raw = r.clone();
        let open_raw = self.kkt(&x, &raw);
        if kind != Kind::Best { for k in 0..n { r[k] += self.defect[k]; } }
        let open_kkt = self.kkt(&x, &r);
        let maxr = (0..n).filter(|&k| !self.rim[k] && !self.held[k]).map(|k| r[k].abs()).fold(0.0, f64::max);
        self.hist.push(serde_json::json!({"eval": self.evals - 1, "kind": format!("{kind:?}"), "round": self.round, "open": open_kkt.len(), "max_abs_r_m": (maxr * 10.0).round() / 10.0}));
        log(format!("  lane E-SOLVE-CONVERGE eval {} kind {kind:?} round {}: open {} {:?} max |r| (counted cells) {maxr:.1} m; F4 ok {f4_ok}", self.evals - 1, self.round, open_kkt.len(),
            open_kkt.iter().take(12).map(|&k| (self.cells[k], r[k].round(), x[k].exp().round())).collect::<Vec<_>>()));
        if kind == Kind::Best {
            self.f4_ok_last = f4_ok;
            self.converged_kkt = open_raw.is_empty() && f4_ok;
            return Act::Stop;
        }
        if let Kind::Fd(gi) = kind {
            let b = self.base.clone().expect("fd without base");
            let grp = self.fd_queue[gi].clone();
            for d in 0..n {
                let dr = r[d] - b.r[d];
                match grp.iter().copied().find(|c| self.sup[d].contains(c)) {
                    Some(c) => { self.jac[d][c] = dr / self.fd_dx[c]; }
                    None => { self.fd_far_max_m = self.fd_far_max_m.max(dr.abs()); }
                }
            }
            for &c in &grp { self.have_col[c] = true; }
            if gi + 1 < self.fd_queue.len() { return self.fd_eval(gi + 1, &b); }
            self.stale = false;
            self.fd_sweeps += 1;
            let jmin = self.p.jmin_m as f64;
            let share = self.p.jshare_permille as f64 / 1000.0;
            for k in 0..n {
                if self.held[k] || !self.fd_queue.iter().any(|g| g.contains(&k)) { continue; }
                let other = (0..n).filter(|&c| c != k && self.have_col[c]).map(|c| self.jac[k][c].abs()).fold(0.0, f64::max);
                if self.jac[k][k] < jmin || self.jac[k][k] < share * other {
                    self.held[k] = true;
                    log(format!("  lane E-SOLVE-CONVERGE: cell {} held (J_dd {:.1} m per unit ln E; largest other entry in its row {other:.1}; jmin {jmin}, jshare {share})", self.cells[k], self.jac[k][k]));
                }
            }
            let est: Vec<usize> = self.fd_queue.iter().flatten().copied().collect();
            log(format!("  lane E-SOLVE-CONVERGE J sweep #{} (stale refreshes {}; round {}, {} columns): own entries {:?}; largest other entry per column {:?}; response outside the support max {:.1} m", self.fd_sweeps, self.fd_done, self.round, est.len(),
                est.iter().map(|&k| (self.cells[k], self.jac[k][k].round())).collect::<Vec<_>>(),
                est.iter().map(|&c| { let (d, v) = (0..n).filter(|&d| d != c).map(|d| (d, self.jac[d][c])).fold((c, 0.0f64), |a, t| if t.1.abs() > a.1.abs() { t } else { a }); (self.cells[c], self.cells[d], v.round()) }).collect::<Vec<_>>(), self.fd_far_max_m));
            return self.step_from_base(&b, false, log);
        }
        // A base evaluation (Reroute, Pre or Step).
        if kind == Kind::Step {
            if let Some(b) = &self.base {
                let dx: Vec<f64> = (0..n).map(|k| x[k] - b.x[k]).collect();
                for d in 0..n {
                    let ss: f64 = self.sup[d].iter().map(|&c| dx[c] * dx[c]).sum();
                    if ss < 1e-10 { continue; }
                    let pr: f64 = self.sup[d].iter().map(|&c| self.jac[d][c] * dx[c]).sum();
                    let f = (r[d] - b.r[d] - pr) / ss;
                    for &c in &self.sup[d] { self.jac[d][c] += f * dx[c]; }
                }
                let smax = self.p.step_max_permille as f64 / 1000.0;
                for k in 0..n {
                    if dx[k].abs() < 1e-12 { continue; }
                    if r[k].abs() > b.r[k].abs() && r[k].abs() > self.tol { self.trust[k] = (self.trust[k] * 0.5).max(1.0 / 64.0); } else { self.trust[k] = (self.trust[k] * 2.0).min(smax); }
                }
                if let Some((act, phi0, phip)) = self.pred.take() {
                    let phi1: f64 = act.iter().map(|&k| r[k] * r[k]).sum();
                    let rho = if phi0 - phip > 1e-9 { (phi0 - phi1) / (phi0 - phip) } else { 1.0 };
                    if rho < 0.25 { self.stale = true; }
                    log(format!("  lane E-SOLVE-CONVERGE step ratio (actual / predicted decrease of sum r^2 over the {} active cells) {rho:.2}{}", act.len(), if rho < 0.25 { " -> J stale" } else { "" }));
                }
            }
        }
        let b = Base { x: x.clone(), r: r.clone(), den: nd.iter().map(|v| v.1).collect(), table: table.clone() };
        self.base = Some(b.clone());
        self.f4_ok_last = f4_ok;
        // Best base of the current round (a reroute resets it: the network changed).
        if kind == Kind::Reroute { self.best = None; }
        let (mm, ms) = self.merit(&x, &r);
        if self.best.as_ref().map(|bb| (mm, ms) < (bb.0, bb.1)).unwrap_or(true) { self.best = Some((mm, ms, table.clone(), flood.clone())); }
        let done = open_kkt.is_empty() && f4_ok;
        let budget_left = self.p.evals_max.saturating_sub(self.evals);
        // A reroute evaluation meeting the stop test is the consistent fixed point (the network is LS1 of the published
        // table): stop there. The final reroute (end_on_reroute) is published whatever its residuals.
        let done_raw = open_raw.is_empty() && f4_ok;
        if kind == Kind::Reroute && self.evals > 1 && (done_raw || self.final_reroute) {
            let done = done_raw;
            self.converged_kkt = done;
            self.stop_reason = if done { "converged at a reroute (network = LS1 of the published table)".into() } else { "final reroute (network = LS1 of the published table; residuals as measured)".into() };
            log(format!("  lane E-SOLVE-CONVERGE round {}: {}; max raw |r| {:.1} m", self.round, self.stop_reason, (0..n).filter(|&k| !self.rim[k] && !self.held[k]).map(|k| raw[k].abs()).fold(0.0, f64::max)));
            return Act::Stop;
        }
        let is_best = self.best.as_ref().map(|b| b.2 == *table).unwrap_or(false);
        if is_best || kind == Kind::Reroute { self.since_best = 0; } else { self.since_best += 1; }
        let stalled = self.p.stall > 0 && self.since_best >= self.p.stall;
        if done || stalled || self.steps >= self.p.steps_max || budget_left <= 2 {
            if done { log(format!("  lane E-SOLVE-CONVERGE round {}: converged (KKT, F4) after {} pre-steps and {} steps", self.round, self.pre, self.steps)); }
            else { log(format!("  lane E-SOLVE-CONVERGE round {}: stopped unconverged after {} steps ({}budget left {budget_left}); open {:?}", self.round, self.steps, if stalled { "stalled, " } else { "" }, open_kkt.iter().map(|&k| (self.cells[k], r[k].round())).collect::<Vec<_>>())); }
            if self.p.ls_each == 0 && self.round + 1 < self.p.reroutes && budget_left > 3 {
                self.pre_reroute_r = Some(raw.clone());
                self.round += 1;
                self.steps = 0;
                self.pre = 0;
                self.since_best = 0;
                self.kind = Kind::Reroute;
                return Act::Next { table: table.clone(), update_flood: true, reroute: true, flood: None };
            }
            if self.p.ls_each == 0 && self.p.end_on_reroute == 1 && budget_left >= 1 {
                // Publish a consistent network: LS1 of the best table of this round (with its flooding offsets).
                self.final_reroute = true;
                self.round += 1;
                self.kind = Kind::Reroute;
                let (bt, bf) = self.best.as_ref().map(|b| (b.2.clone(), Some(b.3.clone()))).unwrap_or((table.clone(), None));
                log(format!("  lane E-SOLVE-CONVERGE: final reroute with the best table of the round (max |r| {:.1} m)", self.best.as_ref().map(|b| b.0).unwrap_or(maxr)));
                return Act::Next { table: bt, update_flood: false, reroute: true, flood: bf };
            }
            let reason = if done { "converged" } else if budget_left <= 2 { "evals_max" } else if stalled { "stalled" } else { "steps_max" };
            return self.stop(done_raw, reason, table, log);
        }
        // (a) per-cell secant pre-steps (before any Jacobian of this round).
        if self.pre < self.p.pre_steps && maxr > self.p.pre_tol_m as f64 && self.fd_sweeps == 0 {
            let mut xn = x.clone();
            let mut moved = 0;
            // Per-cell: only cells whose own |r| exceeds pre_tol step (a local decision; the others wait for the LM phase).
            for k in 0..n {
                if !self.open(k, x[k], r[k]) || r[k].abs() <= self.p.pre_tol_m as f64 { continue; }
                let (s, ns) = crate::t1b_es::secant_step(&self.es, self.sec[k].as_ref(), x[k], nd[k].0, nd[k].1, r[k]);
                self.sec[k] = ns;
                xn[k] = (x[k] + s).clamp(self.lemin, self.lemax);
                moved += 1;
            }
            self.pre += 1;
            log(format!("  lane E-SOLVE-CONVERGE pre-step {} round {}: {moved} cells (secant_diagonal_v1)", self.pre, self.round));
            self.kind = Kind::Pre;
            return Act::Next { table: self.table_of(&b.table, &xn), update_flood: true, reroute: false, flood: None };
        }
        self.step_from_base(&b, true, log)
    }

    /// Columns to (re-)estimate at base `b`: the active cells' columns never measured, or all active columns when J is stale.
    fn fd_plan(&self, b: &Base, refresh: bool) -> Vec<Vec<usize>> {
        let n = self.cells.len();
        let need: Vec<usize> = (0..n).filter(|&k| self.open(k, b.x[k], b.r[k]) && (refresh || !self.have_col[k])).collect();
        self.colours.iter().map(|g| g.iter().copied().filter(|k| need.contains(k)).collect::<Vec<usize>>()).filter(|g| !g.is_empty()).collect()
    }

    fn fd_eval(&mut self, gi: usize, b: &Base) -> Act {
        let h = self.p.fd_step_permille as f64 / 1000.0;
        let mut x = b.x.clone();
        // A cell within h of the upper bound is perturbed downwards (the column divides by the signed dx).
        for &k in &self.fd_queue[gi] { let d = if x[k] + h <= self.lemax { h } else { -h }; x[k] += d; self.fd_dx[k] = d; }
        self.kind = Kind::Fd(gi);
        self.fd_evals += 1;
        Act::Next { table: self.table_of(&b.table, &x), update_flood: false, reroute: false, flood: None }
    }

    /// `flood_ok`: the evaluation just done is `b` itself (so the F4 update may use its rows); false after a colour evaluation.
    fn step_from_base(&mut self, b: &Base, flood_ok: bool, log: &mut dyn FnMut(String)) -> Act {
        let n = self.cells.len();
        // Jacobian columns first when needed and affordable (the sweep plus one step and the best-base evaluation).
        if flood_ok {
            // Stale refreshes are capped by fd_max; sweeps that only fill never-measured columns are not.
            let refresh = self.stale && self.fd_done < self.p.fd_max;
            let plan = self.fd_plan(b, refresh);
            let budget_left = self.p.evals_max.saturating_sub(self.evals);
            if !plan.is_empty() && budget_left >= plan.len() as u32 + 3 {
                if refresh { self.fd_done += 1; }
                self.fd_queue = plan;
                return self.fd_eval(0, b);
            }
        }
        let active: Vec<usize> = (0..n).filter(|&k| self.open(k, b.x[k], b.r[k])).collect();
        self.steps += 1;
        self.kind = Kind::Step;
        if active.is_empty() {
            // Only F4 (or held / rim cells) open: one more base evaluation with the F update.
            self.pred = None;
            return Act::Next { table: b.table.clone(), update_flood: flood_ok, reroute: false, flood: None };
        }
        let m = active.len();
        // A column never measured: a diagonal guess 0.3 x |M - L| (at least 50 m per unit ln E).
        let ja: Vec<Vec<f64>> = active.iter().map(|&d| active.iter().map(|&c| if self.have_col[c] { self.jac[d][c] } else if c == d { (0.3 * b.den[d].abs()).max(50.0) } else { 0.0 }).collect()).collect();
        let ra: Vec<f64> = active.iter().map(|&d| b.r[d]).collect();
        let lm = self.p.lm_permille as f64 / 1000.0;
        let mut a = vec![vec![0.0; m]; m];
        let mut g = vec![0.0; m];
        for i in 0..m { for j in 0..m { let mut s = 0.0; for d in 0..m { s += ja[d][i] * ja[d][j]; } a[i][j] = s; } let mut s = 0.0; for d in 0..m { s += ja[d][i] * ra[d]; } g[i] = -s; }
        for (i, row) in a.iter_mut().enumerate() { let dd = row[i]; row[i] = dd * (1.0 + lm) + 1e-6; }
        let dx = solve_dense(&a, &g).unwrap_or_else(|| (0..m).map(|i| -ra[i] / ja[i][i].max(1.0)).collect());
        let mut x = b.x.clone();
        let mut dxa = vec![0.0; m];
        for (i, &k) in active.iter().enumerate() {
            let t = self.trust[k];
            let s = if dx[i].is_finite() { dx[i].clamp(-t, t) } else { 0.0 };
            x[k] = (x[k] + s).clamp(self.lemin, self.lemax);
            dxa[i] = x[k] - b.x[k];
        }
        let phi0: f64 = ra.iter().map(|v| v * v).sum();
        let phip: f64 = (0..m).map(|d| { let v = ra[d] + (0..m).map(|c| ja[d][c] * dxa[c]).sum::<f64>(); v * v }).sum();
        self.pred = Some((active.clone(), phi0, phip));
        log(format!("  lane E-SOLVE-CONVERGE step {} round {}: active {} cells; predicted max |r| {:.0} m; steps {:?}", self.steps, self.round, m,
            (0..m).map(|d| (ra[d] + (0..m).map(|c| ja[d][c] * dxa[c]).sum::<f64>()).abs()).fold(0.0, f64::max),
            active.iter().map(|&k| (self.cells[k], ((x[k] - b.x[k]) * 1000.0).round() / 1000.0)).take(24).collect::<Vec<_>>()));
        Act::Next { table: self.table_of(&b.table, &x), update_flood: flood_ok, reroute: false, flood: None }
    }

    /// The layer-header record (lane_es.conv).
    pub fn record(&self) -> Value {
        let held: Vec<Value> = (0..self.cells.len()).filter(|&k| self.held[k]).map(|k| serde_json::json!({"cell": self.cells[k], "rim": self.rim[k], "j_dd_m": self.jac[k][k]})).collect();
        serde_json::json!({"rule": ESC_RULE, "params": self.p, "evals": self.evals, "fd_evals": self.fd_evals, "rounds": self.round + 1, "jacobian_sweeps": self.fd_sweeps, "stale_refreshes": self.fd_done, "converged_kkt": self.converged_kkt, "f4_ok": self.f4_ok_last,
            "stop_reason": self.stop_reason, "held_cells": held, "fd_response_outside_support_max_m": self.fd_far_max_m,
            "jacobian_diag_m": self.cells.iter().enumerate().filter(|(k, _)| self.have_col[*k]).map(|(k, c)| (c.to_string(), (self.jac[k][k] * 10.0).round() / 10.0)).collect::<BTreeMap<String, f64>>(),
            "history": self.hist, "note": "EXPLORATORY opt-in (lane E-SOLVE-CONVERGE); network held between reroutes and published as evaluated"})
    }
}

fn solve_dense(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut m: Vec<Vec<f64>> = (0..n).map(|i| { let mut r = a[i].clone(); r.push(b[i]); r }).collect();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| m[i][c].abs().total_cmp(&m[j][c].abs()).then(j.cmp(&i)))?;
        if m[p][c].abs() < 1e-12 { return None; }
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

    /// A 4x4 block of cells with hat-like coupling stronger than the diagonal on the edge rows (residual of cell d =
    /// sum over its 3x3 of w (x_c - x*_c) x relief(x_d), the relief growing with E as in the measured J) converges within
    /// tolerance through the driver (secant pre-steps, coloured J, LM steps), always stopping on a base evaluation.
    fn run(target_shift: Option<(usize, f64)>, weak_edge: bool) -> (BTreeMap<u32, f64>, u32) {
        let mut p = EscParamsV1::v1();
        p.evals_max = 80;
        p.steps_max = 30;
        // The weak-edge block is ill-conditioned on purpose: the v1 stall stop (5 base evaluations without a new best)
        // ends it at about 30 m, so the stall stop is off here.
        p.stall = 0;
        let n = 16usize;
        let cells: Vec<u32> = (0..n as u32).collect();
        let ij = |c: u32| ((c / 4) as i64, (c % 4) as i64);
        let xs: Vec<f64> = (0..n).map(|k| (300.0 + 40.0 * k as f64).ln()).collect();
        let resid = |x: &[f64]| -> Vec<f64> {
            (0..n).map(|d| {
                let (di, dj) = ij(d as u32);
                let mut s = 0.0;
                for c in 0..n { let (ci, cj) = ij(c as u32); if (di - ci).abs() <= 1 && (dj - cj).abs() <= 1 { let w = if c == d { if di == 0 && weak_edge { 60.0 } else { 300.0 } } else if weak_edge { 80.0 } else { 40.0 }; s += w * (x[c] - xs[c]) * (0.5 * (x[d] - xs[d])).exp(); } }
                s - target_shift.filter(|t| t.0 == d).map(|t| t.1).unwrap_or(0.0)
            }).collect()
        };
        let es = crate::t1b_es::EsParamsV1::v1();
        let mut esc = Esc::new(&p, &es, &cells, &[], &ij, 1000, 20.0, 5.0, 20000.0);
        let mut table: BTreeMap<u32, f64> = cells.iter().map(|&c| (c, 150.0)).collect();
        let mut last_kind = Kind::Reroute;
        for _ in 0..100 {
            let x: Vec<f64> = cells.iter().map(|c| table[c].ln()).collect();
            let r = resid(&x);
            let obs: Vec<CellObs> = cells.iter().map(|&c| CellObs { cell: c, diff_m: r[c as usize], e: table[&c], num: 500.0, den: 500.0 + r[c as usize] }).collect();
            last_kind = esc.kind();
            match esc.after_eval(&obs, &table, &BTreeMap::new(), true, &mut |s| if std::env::var("ESC_TEST_LOG").is_ok() { println!("{s}") }) {
                Act::Stop => break,
                Act::Next { table: t, .. } => table = t,
            }
        }
        assert!(!matches!(last_kind, Kind::Fd(_)));
        let x: Vec<f64> = cells.iter().map(|c| table[c].ln()).collect();
        let r = resid(&x);
        assert!(r.iter().all(|v| v.abs() <= 20.0), "residuals {r:?}");
        assert!(esc.converged_kkt());
        (table, esc.evals)
    }

    #[test]
    fn converges_on_coupled_block() { let (_, ev) = run(None, true); assert!(ev <= 50, "evals {ev}"); }

    /// Locality where cells control their own means (diagonal 300 m, neighbours 40 m per unit ln E): shifting one corner
    /// cell's target by 150 m moves no cell at Chebyshev distance 3 by |d ln E| >= 0.1 (a 30 m own-mean equivalent). With
    /// weak edge rows (own 60 m, neighbours 80 m: the measured rim / coastal pattern) the solution itself is non-local.
    #[test]
    fn target_shift_stays_local() {
        let (a, _) = run(None, false);
        let (b, _) = run(Some((0, 150.0)), false);
        for c in 0..16u32 { let (i, j) = ((c / 4) as i64, (c % 4) as i64); if i.max(j) >= 3 { assert!((a[&c].ln() - b[&c].ln()).abs() < 0.1, "cell {c}: {} vs {}", a[&c], b[&c]); } }
    }

    #[test]
    fn params_validate_and_parse() {
        let mut es = crate::t1b_es::EsParamsV1::v1();
        assert!(parse_key(&mut es, "esc_steps", 3).is_err());
        parse_key(&mut es, "esc", 1).unwrap();
        parse_key(&mut es, "esc_reroutes", 2).unwrap();
        assert_eq!(es.conv.as_ref().unwrap().reroutes, 2);
        assert!(parse_key(&mut es, "esc_reroutes", 0).is_err());
        parse_key(&mut es, "esc", 0).unwrap();
        assert!(es.conv.is_none());
    }

    /// Off = no field: the es params without esc serialize exactly as before (recipe digests unchanged).
    #[test]
    fn es_without_esc_serializes_unchanged() {
        let v = serde_json::to_value(crate::t1b_es::EsParamsV1::v1()).unwrap();
        assert!(v.get("conv").is_none());
    }
}
