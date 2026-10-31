//! SEED-RELIEF-LAW (lane RL, 2026-10-12; design 15:55: "get them working and
//! decent simulations please"): law id `stream-power-subgrid-nonlinear-v1`,
//! the two-elevation subgrid stream power of lane OG
//! (`stream-power-subgrid-v1`, pf-geodynamics orogen.rs) generalised to a
//! slope exponent n >= 1, applied CONSISTENTLY to the resolved trunk and to
//! the closure of the unresolved (subgrid) channel network. Opt-in only; no
//! legacy path calls this module.
//!
//! WHY (diagnosis, `unpublished-design-record`): under
//! the n = 1 law with the toy K (2e-6 /yr), stable continental interiors
//! of lane OG's six L1 seeds erode at 47-55 m/My for 100 My (4.7-6.5 km of
//! crust), because the rate is LINEAR in channel steepness: a low-relief
//! interior erodes at a fixed fraction of mountain rates, and isostatic
//! rebound returns 85 % of every eroded metre as new elevation to erode.
//! Mean-flow erosion that rises faster than linearly with steepness is the
//! accepted consequence of stochastic floods acting against an incision
//! threshold (Lague 2014, ESPL 39:38; DiBiase and Whipple 2011, JGR 116
//! F04036; Harel, Mudd and Attal 2016, Geomorphology 268:184; all recalled,
//! UNVERIFIED), and n > 1 is what lets post-orogenic relief persist
//! (Baldwin, Whipple and Tucker 2003, JGR 108(B3), recalled, UNVERIFIED).
//!
//! THE LAW. Concavity m/n = 0.5 is fixed by the law version (the legacy
//! m = 0.5, n = 1), so m = n/2. The erodibility is
//!     K_n = K * k_ref^(1 - n),
//! with K the universe's stream-power erodibility (n = 1 units) and k_ref a
//! reference channel steepness k_s = S A^(m/n): the n-law equals the legacy
//! n = 1 law exactly where k_s = k_ref, erodes slower below it and faster
//! above it.
//! - Trunk level c (the resolved coarse network, receiver-first, implicit):
//!       dc/dt = -K_n A^m ((c - c_receiver)/L)^n.
//!   Per cell, with x = c' - b and x0 = c - b (b the receiver's new trunk
//!   level, or the fill level of a wet receiver), x + F x^n = x0 with
//!   F = K_n A^m dt / L^n: one root in [0, x0] (closed form for n = 1, 2;
//!   monotone Newton from x0 otherwise). Base levels, the monotone-trunk
//!   rule (a trunk at or below its receiver is raised to it, counted) and
//!   the wet-cell rules are lane OG's, unchanged.
//! - Cell mean z (the unresolved network, backward Euler):
//!       dz/dt = -K_n ((z - c)/Lambda)^n,  Lambda = ln(dx/x_c),
//!   dx = sqrt(cell area). Derivation: a quasi-steady self-similar subgrid
//!   network a(x) = x^2 eroding at E under the same law has slope
//!   S(x) = (E/K_n)^(1/n)/x, so the relief from the cell scale dx to the
//!   channel heads x_c is R = (E/K_n)^(1/n) * integral_{x_c}^{dx} dx/x =
//!   (E/K_n)^(1/n) Lambda, i.e. E = K_n (R/Lambda)^n. For n = 1 this is
//!   E = R/tau_sub with tau_sub = Lambda/K: exactly stream-power-subgrid-v1.
//!   With n = 1 this module uses lane OG's arithmetic verbatim, so the two
//!   laws are bit-identical there (pf-lab tests/relief_law.rs).
//! - Realisation: as lane OG, the pf-surface erosion stage (which stays the
//!   single owner of transport, deposition and the sediment ledgers) is
//!   handed a per-cell EFFECTIVE contributing area that makes its n = 1
//!   implicit update remove exactly the target cell-mean lowering; the
//!   stage's own K and m are inputs here (`stage_k_per_my`, `stage_m`).
//!
//! LIMITS (declared): quasi-steady subgrid network (no subgrid knickpoint
//! lag); one rock type; mean-flow law (the threshold and flood statistics
//! are folded into n, not resolved); no transport-limited alluvial reaches
//! except through the existing deposition laws.

/// Law id (the typed law set lives in pf-contracts laws.rs).
pub const NONLINEAR_SUBGRID_LAW_ID: &str = "stream-power-subgrid-nonlinear-v1";
/// Newton iteration cap for n not in {1, 2} (numerical control).
pub const ROOT_MAX_ITERS: u32 = 200;

/// Resolved parameters of `stream-power-subgrid-nonlinear-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearSubgridParams {
    /// Universe stream-power erodibility in n = 1 units (per My).
    pub k_per_my: f64,
    /// Slope exponent n (>= 1); m = n/2.
    pub slope_exponent_n: f64,
    /// Reference channel steepness k_ref (m, theta = 0.5 units) at which
    /// the n-law equals the n = 1 law.
    pub reference_steepness_m: f64,
    /// Channel-head length x_c of the subgrid network (m).
    pub channel_head_length_m: f64,
}

impl NonlinearSubgridParams {
    pub fn validate(&self) -> Result<(), String> {
        let ok = |v: f64| v.is_finite();
        if !ok(self.k_per_my) || self.k_per_my < 0.0 {
            return Err(format!(
                "k_per_my {} must be finite and >= 0",
                self.k_per_my
            ));
        }
        if !ok(self.slope_exponent_n) || !(1.0..=4.0).contains(&self.slope_exponent_n) {
            return Err(format!(
                "slope_exponent_n {} outside [1, 4]",
                self.slope_exponent_n
            ));
        }
        if !ok(self.reference_steepness_m) || self.reference_steepness_m <= 0.0 {
            return Err(format!(
                "reference_steepness_m {} must be finite and > 0",
                self.reference_steepness_m
            ));
        }
        if !ok(self.channel_head_length_m) || self.channel_head_length_m <= 0.0 {
            return Err(format!(
                "channel_head_length_m {} must be finite and > 0",
                self.channel_head_length_m
            ));
        }
        Ok(())
    }

    /// Area exponent m = n/2 (concavity 0.5 fixed by the law version).
    pub fn m(&self) -> f64 {
        0.5 * self.slope_exponent_n
    }

    /// K_n = K k_ref^(1 - n) (units m^(1-n) per My).
    pub fn k_n_per_my(&self) -> f64 {
        if self.slope_exponent_n == 1.0 {
            self.k_per_my
        } else {
            self.k_per_my * self.reference_steepness_m.powf(1.0 - self.slope_exponent_n)
        }
    }

    fn is_linear(&self) -> bool {
        self.slope_exponent_n == 1.0
    }
}

/// Lambda = integral_{x_c}^{dx} x^(-2m/n) dx = ln(dx/x_c) (m/n = 0.5);
/// 0 when dx <= x_c (resolved limit).
pub fn subgrid_lambda(dx_m: f64, x_c_m: f64) -> f64 {
    if !(dx_m > x_c_m) {
        return 0.0;
    }
    (dx_m / x_c_m).ln()
}

/// The root x in [0, x0] of x + f x^n = x0 (x0 >= 0, f >= 0, n >= 1): the
/// backward-Euler update of dx/dt = -k x^n with f = k dt. Closed forms for
/// n = 1 and n = 2 (the stable conjugate form); monotone Newton from x0
/// otherwise (g(x) = x + f x^n - x0 is increasing and convex on [0, x0]
/// with g(x0) >= 0, so Newton from x0 decreases monotonically to the root).
pub fn implicit_power_root(x0: f64, f: f64, n: f64) -> f64 {
    if !(x0 > 0.0) || !(f > 0.0) {
        return x0.max(0.0);
    }
    if n == 1.0 {
        return x0 / (1.0 + f);
    }
    if n == 2.0 {
        return 2.0 * x0 / (1.0 + (1.0 + 4.0 * f * x0).sqrt());
    }
    let mut x = x0;
    for _ in 0..ROOT_MAX_ITERS {
        let g = x + f * x.powf(n) - x0;
        let dg = 1.0 + n * f * x.powf(n - 1.0);
        let next = (x - g / dg).max(0.0);
        if (x - next).abs() <= 1e-14 * x0 {
            return next;
        }
        x = next;
    }
    x
}

/// Inputs of one step, all per cell (indices are cell ids). The same
/// fields lane OG's `pf_geodynamics::orogen::SubgridInputs` carries.
pub struct SubgridInputs<'a> {
    /// Receivers before donors (`Routing::receiver_first`).
    pub receiver_first: &'a [usize],
    pub flow_to: &'a [Option<usize>],
    pub contributing_area_m2: &'a [f64],
    /// Centre distance to the receiver (m), as the erosion stage measures
    /// it (`center_distance_m(..).max(1.0)`); unused where flow_to is None.
    pub reach_length_m: &'a [f64],
    /// Ocean or lake cells (no fluvial lowering; base level = fill level).
    pub wet: &'a [bool],
    pub fill_level_m: &'a [f64],
    /// Linear cell size dx = sqrt(cell area) (m).
    pub cell_size_m: &'a [f64],
}

/// Result of one step (same meaning as lane OG's `SubgridStep`).
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearSubgridStep {
    /// Trunk level after the step (m); = mean for wet cells.
    pub trunk_m: Vec<f64>,
    /// Target cell-mean fluvial lowering (m, >= 0).
    pub target_lowering_m: Vec<f64>,
    /// Effective contributing area (m^2) that makes the erosion stage's
    /// n = 1 implicit update remove exactly the target.
    pub effective_area_m2: Vec<f64>,
    /// Cells whose target was capped so the mean stays above its
    /// receiver's lowered mean (the erosion stage cannot lower below it).
    pub capped_cells: u64,
    /// Cells whose step-start trunk stood at or below the receiver's trunk
    /// and was raised to it (monotone trunk profile).
    pub trunk_raised_cells: u64,
}

/// One step of `stream-power-subgrid-nonlinear-v1`: trunk incision, the
/// cell-mean target, and the effective areas that realise it through the
/// erosion stage whose (n = 1) erodibility and area exponent are
/// `stage_k_per_my` and `stage_m`.
#[allow(clippy::too_many_arguments)]
pub fn subgrid_step(
    inp: &SubgridInputs<'_>,
    z_mean_m: &[f64],
    relief_m: &[f64],
    p: &NonlinearSubgridParams,
    stage_k_per_my: f64,
    stage_m: f64,
    dt_my: f64,
) -> Result<NonlinearSubgridStep, String> {
    subgrid_step_with_floor(inp, z_mean_m, relief_m, p, stage_k_per_my, stage_m, dt_my, None)
}

/// `subgrid_step` with an optional per-cell FLOOR on the cell-mean target
/// (lane UP, SEED-UPLANDS 2026-10-13): where `floor_lowering_m[i]` exceeds
/// the stream-power target, the floor is used instead (the faster process
/// sets the rate; e.g. `critical-slope-envelope-v1`, uplands.rs), before
/// the same receiver cap and effective-area realisation. `None` is exactly
/// `subgrid_step` (same arithmetic, bit for bit).
#[allow(clippy::too_many_arguments)]
pub fn subgrid_step_with_floor(
    inp: &SubgridInputs<'_>,
    z_mean_m: &[f64],
    relief_m: &[f64],
    p: &NonlinearSubgridParams,
    stage_k_per_my: f64,
    stage_m: f64,
    dt_my: f64,
    floor_lowering_m: Option<&[f64]>,
) -> Result<NonlinearSubgridStep, String> {
    subgrid_step_with_floor_and_runoff(inp, z_mean_m, relief_m, p, stage_k_per_my, stage_m, dt_my, floor_lowering_m, None)
}

/// I11 (lane CLIMATE-GENESIS-2, 2026-10-20; law
/// `cell-runoff-subgrid-erosivity-v1` under `epoch-climate-runoff-v1`):
/// `subgrid_step_with_floor` with an optional per-cell runoff ratio
/// R_i / r_ref for the CELL-MEAN closure. The unresolved network inside a
/// cell drains that cell's own area, so its discharge at subgrid area a is
/// a R_i and the derivation above gives E = K_n (R_i/r_ref)^m (R/Lambda)^n
/// with m = n/2 (the law's own area exponent; 0.5 at n = 1). The trunk is
/// unchanged here (it already reads the runoff-weighted catchment area in
/// `contributing_area_m2`). `None` is exactly `subgrid_step_with_floor`
/// (same arithmetic, bit for bit).
#[allow(clippy::too_many_arguments)]
pub fn subgrid_step_with_floor_and_runoff(
    inp: &SubgridInputs<'_>,
    z_mean_m: &[f64],
    relief_m: &[f64],
    p: &NonlinearSubgridParams,
    stage_k_per_my: f64,
    stage_m: f64,
    dt_my: f64,
    floor_lowering_m: Option<&[f64]>,
    cell_runoff_ratio: Option<&[f64]>,
) -> Result<NonlinearSubgridStep, String> {
    p.validate()?;
    let n = z_mean_m.len();
    if let Some(r) = cell_runoff_ratio {
        if r.len() != n {
            return Err(format!("subgrid cell runoff ratio: {} != {n}", r.len()));
        }
        if let Some(i) = r.iter().position(|v| !(v.is_finite() && *v >= 0.0)) {
            return Err(format!("subgrid cell runoff ratio[{i}] = {} must be finite and >= 0", r[i]));
        }
    }
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
    let mut trunk = vec![0.0f64; n];
    let mut target = vec![0.0f64; n];
    let mut eff = vec![0.0f64; n];
    let mut mean_after = z_mean_m.to_vec();
    let mut capped = 0u64;
    let mut raised = 0u64;
    // Trunk level at the step start: mean minus relief (never above mean).
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
            continue; // terminal floor: no fluvial change
        };
        let len = inp.reach_length_m[i];
        let base_t = if inp.wet[j] {
            inp.fill_level_m[j]
        } else {
            trunk[j]
        };
        if c0[i] > base_t {
            if p.is_linear() {
                // Lane OG's arithmetic, verbatim (bit-identity at n = 1).
                let fc = p.k_per_my * inp.contributing_area_m2[i].powf(m) * dt_my / len;
                trunk[i] = (c0[i] + fc * base_t) / (1.0 + fc);
            } else {
                let fc = kn * inp.contributing_area_m2[i].powf(m) * dt_my / len.powf(ne);
                trunk[i] = base_t + implicit_power_root(c0[i] - base_t, fc, ne);
            }
        } else {
            trunk[i] = base_t.min(z_mean_m[i]);
            raised += 1;
        }
        // Cell mean: backward Euler toward the new trunk.
        let lambda = subgrid_lambda(inp.cell_size_m[i], x_c);
        let y0 = (z_mean_m[i] - trunk[i]).max(0.0);
        // I11: the cell's runoff erosivity (R_i/r_ref)^m on the closure's K.
        let (k_lin, kn_i) = match cell_runoff_ratio {
            None => (p.k_per_my, kn),
            Some(r) => {
                let e = r[i].powf(m);
                (p.k_per_my * e, kn * e)
            }
        };
        let mut d = if p.is_linear() {
            // tau_sub = Lambda / K; phi = r/(1+r), r = dt/tau (lane OG).
            let tau = if lambda == 0.0 {
                0.0
            } else if !(k_lin > 0.0) {
                f64::INFINITY
            } else {
                lambda / k_lin
            };
            let phi = if tau == 0.0 {
                1.0
            } else if tau.is_infinite() {
                0.0
            } else {
                let r = dt_my / tau;
                r / (1.0 + r)
            };
            y0 * phi
        } else if lambda == 0.0 {
            y0 // resolved limit: the mean is the trunk
        } else {
            let fm = kn_i * dt_my / lambda.powf(ne);
            y0 - implicit_power_root(y0, fm, ne)
        };
        if let Some(f) = floor_lowering_m {
            if f[i] > d {
                d = f[i];
            }
        }
        // What the erosion stage will use as this cell's base.
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
            // f = K A_eff^m dt / L  =>  A_eff = (f L / (K dt))^(1/m)
            eff[i] = (f * len / (stage_k_per_my * dt_my)).powf(1.0 / stage_m);
            mean_after[i] = z_mean_m[i] - d;
        }
        target[i] = d;
        if trunk[i] > mean_after[i] {
            trunk[i] = mean_after[i];
        }
    }
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

    fn params(n: f64) -> NonlinearSubgridParams {
        NonlinearSubgridParams {
            k_per_my: 2.0,
            slope_exponent_n: n,
            reference_steepness_m: 500.0,
            channel_head_length_m: 100.0,
        }
    }

    #[test]
    fn root_solves_the_implicit_power_equation() {
        for &n in &[1.0, 1.5, 2.0, 2.5, 3.0, 4.0] {
            for &(x0, f) in &[(1000.0, 1e-3), (1.0, 10.0), (5000.0, 1e-7), (3.0, 0.0)] {
                let x = implicit_power_root(x0, f, n);
                assert!((0.0..=x0).contains(&x), "n {n} x0 {x0} f {f}: {x}");
                let resid = x + f * x.powf(n) - x0;
                assert!(
                    resid.abs() <= 1e-10 * x0,
                    "n {n} x0 {x0} f {f}: resid {resid}"
                );
            }
        }
        assert_eq!(implicit_power_root(0.0, 1.0, 2.0), 0.0);
        assert_eq!(implicit_power_root(-1.0, 1.0, 2.0), 0.0);
    }

    #[test]
    fn k_n_matches_the_legacy_law_at_the_reference_steepness() {
        let p = params(2.0);
        // E_n = K_n k^n equals E_1 = K k at k = k_ref.
        let k = p.reference_steepness_m;
        assert!((p.k_n_per_my() * k * k - p.k_per_my * k).abs() < 1e-12 * p.k_per_my * k);
        // Slower below, faster above.
        assert!(p.k_n_per_my() * 50.0f64.powi(2) < p.k_per_my * 50.0);
        assert!(p.k_n_per_my() * 1000.0f64.powi(2) > p.k_per_my * 1000.0);
        assert_eq!(params(1.0).k_n_per_my(), 2.0);
    }

    #[test]
    fn validation_is_loud() {
        assert!(params(2.0).validate().is_ok());
        assert!(params(0.5).validate().is_err());
        assert!(params(5.0).validate().is_err());
        let mut p = params(2.0);
        p.reference_steepness_m = 0.0;
        assert!(p.validate().is_err());
        let mut p = params(2.0);
        p.channel_head_length_m = f64::NAN;
        assert!(p.validate().is_err());
    }

    /// Two cells: a land cell draining to an ocean cell at 0 m.
    fn two_cell(z: f64, relief: f64, p: &NonlinearSubgridParams, dt: f64) -> NonlinearSubgridStep {
        let flow_to = [Some(1usize), None];
        let area = [5.0e9, 5.0e9];
        let reach = [70_000.0, 1.0];
        let wet = [false, true];
        let fill = [z, 0.0];
        let size = [70_000.0, 70_000.0];
        let inp = SubgridInputs {
            receiver_first: &[1, 0],
            flow_to: &flow_to,
            contributing_area_m2: &area,
            reach_length_m: &reach,
            wet: &wet,
            fill_level_m: &fill,
            cell_size_m: &size,
        };
        subgrid_step(&inp, &[z, 0.0], &[relief, 0.0], p, 2.0, 0.5, dt).unwrap()
    }

    #[test]
    fn single_cell_trunk_and_mean_follow_the_closed_forms() {
        let p = params(2.0);
        let (z, r, dt) = (1000.0, 300.0, 2.0);
        let s = two_cell(z, r, &p, dt);
        // Trunk: x + F x^2 = x0 with F = K_n A dt / L^2, base 0.
        let f = p.k_n_per_my() * 5.0e9 * dt / 70_000.0f64.powi(2);
        let x0 = z - r;
        let x = 2.0 * x0 / (1.0 + (1.0 + 4.0 * f * x0).sqrt());
        assert!((s.trunk_m[0] - x).abs() < 1e-9, "{} vs {x}", s.trunk_m[0]);
        // Mean: y + Fm y^2 = z - c', Fm = K_n dt / Lambda^2.
        let lam = (700.0f64).ln();
        let fm = p.k_n_per_my() * dt / (lam * lam);
        let y0 = z - x;
        let y = 2.0 * y0 / (1.0 + (1.0 + 4.0 * fm * y0).sqrt());
        assert!((s.target_lowering_m[0] - (y0 - y)).abs() < 1e-9);
        assert_eq!(s.capped_cells, 0);
        // The effective area reproduces the target in the stage's update:
        // z' = (z + f_eff * 0)/(1 + f_eff), f_eff = K A_eff^0.5 dt / L.
        let f_eff = 2.0 * s.effective_area_m2[0].sqrt() * dt / 70_000.0;
        let realised = z - z / (1.0 + f_eff);
        assert!((realised - s.target_lowering_m[0]).abs() < 1e-9 * z);
    }

    #[test]
    fn n2_erodes_low_relief_slower_and_steep_relief_faster_than_n1() {
        let (dt, z) = (1.0, 2000.0);
        // Low relief: subgrid relief 100 m (k_s = 100/6.55 ~ 15 m << 500 m).
        let a1 = two_cell(z, 100.0, &params(1.0), dt);
        let a2 = two_cell(z, 100.0, &params(2.0), dt);
        // Compare the CELL-MEAN rate at the same trunk: use a cell whose
        // trunk is already graded (z - r = 0 at base): relief = z.
        let low1 = two_cell(100.0, 100.0, &params(1.0), dt).target_lowering_m[0];
        let low2 = two_cell(100.0, 100.0, &params(2.0), dt).target_lowering_m[0];
        assert!(low2 < 0.2 * low1, "n2 {low2} vs n1 {low1}");
        // Very steep: relief 6 km above a graded trunk (k_s ~ 920 m > 500 m).
        let hi1 = two_cell(6000.0, 6000.0, &params(1.0), dt).target_lowering_m[0];
        let hi2 = two_cell(6000.0, 6000.0, &params(2.0), dt).target_lowering_m[0];
        assert!(hi2 > hi1, "n2 {hi2} vs n1 {hi1}");
        assert!(a1.trunk_m[0] < z && a2.trunk_m[0] < z);
    }

    /// I11: a unit runoff ratio is the reference law bit for bit; a wetter
    /// cell relaxes faster, a drier one slower, with the closed form
    /// y + (K_n (R/r_ref)^m dt / Lambda^n) y^n = y0; the trunk is untouched.
    #[test]
    fn cell_runoff_ratio_scales_only_the_cell_mean_closure() {
        for &ne in &[1.0, 2.0] {
            let p = params(ne);
            let (z, r, dt) = (1500.0, 600.0, 2.0);
            let flow_to = [Some(1usize), None];
            let area = [5.0e9, 5.0e9];
            let reach = [70_000.0, 1.0];
            let wet = [false, true];
            let fill = [z, 0.0];
            let size = [70_000.0, 70_000.0];
            let inp = SubgridInputs {
                receiver_first: &[1, 0],
                flow_to: &flow_to,
                contributing_area_m2: &area,
                reach_length_m: &reach,
                wet: &wet,
                fill_level_m: &fill,
                cell_size_m: &size,
            };
            let run = |ratio: Option<&[f64]>| {
                subgrid_step_with_floor_and_runoff(&inp, &[z, 0.0], &[r, 0.0], &p, 2.0, 0.5, dt, None, ratio).unwrap()
            };
            let base = run(None);
            assert_eq!(run(Some(&[1.0, 1.0])), base, "ratio 1 must be the reference law bit for bit (n {ne})");
            let wet_s = run(Some(&[4.0, 1.0]));
            let dry_s = run(Some(&[0.1, 1.0]));
            assert_eq!(wet_s.trunk_m[0].to_bits(), base.trunk_m[0].to_bits());
            assert!(wet_s.target_lowering_m[0] > base.target_lowering_m[0]);
            assert!(dry_s.target_lowering_m[0] < base.target_lowering_m[0]);
            if ne == 2.0 {
                let lam = (700.0f64).ln();
                let fm = p.k_n_per_my() * 4.0f64.powf(p.m()) * dt / (lam * lam);
                let y0 = z - wet_s.trunk_m[0];
                let y = 2.0 * y0 / (1.0 + (1.0 + 4.0 * fm * y0).sqrt());
                assert!((wet_s.target_lowering_m[0] - (y0 - y)).abs() < 1e-9);
            }
        }
        let p = params(2.0);
        let inp_bad = [f64::NAN];
        let flow_to = [None];
        let one = [1.0];
        let inp = SubgridInputs { receiver_first: &[0], flow_to: &flow_to, contributing_area_m2: &one, reach_length_m: &one, wet: &[false], fill_level_m: &one, cell_size_m: &one };
        assert!(subgrid_step_with_floor_and_runoff(&inp, &one, &one, &p, 2.0, 0.5, 1.0, None, Some(&inp_bad)).is_err());
    }

    #[test]
    fn deterministic_and_never_below_receiver() {
        let p = params(2.0);
        let a = two_cell(800.0, 50.0, &p, 2.0);
        let b = two_cell(800.0, 50.0, &p, 2.0);
        assert_eq!(a, b);
        assert!(a.trunk_m[0] >= 0.0);
        assert!(800.0 - a.target_lowering_m[0] > 0.0);
    }
}
