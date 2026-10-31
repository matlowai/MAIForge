//! T1b 0.9.4 lane TE (2026-10-14; design review consult 3 Q4, `unpublished-design-record`): the
//! cell-local E solve `solve_cell`. LABEL: pre-signoff, fidelity tier Preview, derived, never canon; construction (no
//! SimTime). EXPLORATORY and opt-in (`params.ls.tb.tc.te.cell_solve = 1`, key `ls_cell=1`); without it `solve_e` is the
//! chart-wide damped Jacobi, byte-identical. Pre-registered in `unpublished-design-record`.
//!
//! Why (lane TC): in the chart-wide solve every E evaluation re-runs LS1 over the whole chart, warm-started from the
//! previous evaluation, so the published network remembers the whole E trajectory: a 1% change of the starting table
//! re-routed 18-19% of land receivers, and freezing the receivers late did not damp it. Here the T0 cell is the solve unit:
//! cells are visited downstream-first, a cell's routing is recomputed cold from the built network (a function of the
//! current tables, not of the path that led there), and its E and F are 1-D roots on fixed grids, so the answer does not
//! depend on where the iteration started. The only couplings are the outlet levels and inflows of neighbouring cells and the
//! hat blend across cell edges, resolved by Gauss-Seidel sweeps.

use crate::t0::{CellWater, T0};
use crate::t1b::{closure_rate_m_per_my, T0Extras, LAKE, MARINE};
use crate::t1b_rivernet::*;
use crate::t1b_view::SkeletonView;
use crate::*;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use vz_tilekey::chart::Chart;

/// Per-cell sums of one evaluation: hat [w, w z, w ref, w trunk] over land posts; plain [posts, land posts, sum z, sum ref,
/// sum bil, land posts without edifice] over the cell's own posts (as solve_e).
#[derive(Clone, Copy, Default, Debug)]
pub struct Meas { pub hat: [f64; 4], pub plain: [f64; 6] }

impl Meas {
    pub fn m(&self) -> f64 { if self.hat[0] > 0.0 { self.hat[1] / self.hat[0] } else { 0.0 } }
    pub fn t(&self) -> f64 { if self.hat[0] > 0.0 { self.hat[2] / self.hat[0] } else { 0.0 } }
    pub fn l(&self) -> f64 { if self.hat[0] > 0.0 { self.hat[3] / self.hat[0] } else { 0.0 } }
    pub fn land_fraction(&self) -> f64 { if self.plain[0] > 0.0 { self.plain[1] / self.plain[0] } else { 0.0 } }
}

/// Solve context shared by every evaluation.
struct Ctx<'a> {
    t0: &'a T0,
    extras: &'a T0Extras,
    skel: &'a SkeletonView,
    chart: &'a Chart,
    net0: &'a Network,
    nctx: &'a NodeCtx<'a>,
    p: &'a RiverNetParamsV1,
    co: &'a CoInputs<'a>,
    ls: &'a crate::t1b_ls::LsParamsV1,
    own: Vec<f64>,
    coastal: BTreeMap<u32, bool>,
    lane_co: bool,
}

impl<'a> Ctx<'a> {
    /// The reference of the F6 mean at a post (as solve_e's refv_of).
    fn refv(&self, lp: &LatPrep, i: usize) -> f64 {
        let c = lp.cell[i];
        let cm = self.t0.elevation_m[c as usize] as f64;
        match self.p.f6_reference {
            1 => cm,
            2 => if self.coastal.get(&c).copied().unwrap_or(false) && self.t0.water(c) == CellWater::Land { cm } else { lp.bil[i] as f64 },
            _ => lp.bil[i] as f64,
        }
    }

    /// Per-cell sums over a lattice (only cell `only` when given).
    fn measure(&self, lp: &LatPrep, trunk: &[f64], s: &Surface, only: Option<u32>) -> BTreeMap<u32, Meas> {
        let mut out: BTreeMap<u32, Meas> = BTreeMap::new();
        for i in 0..lp.lat.nodes() {
            let land = if self.p.coast.is_some() { s.fixed[i] != MARINE && s.fixed[i] != LAKE } else { lp.water[i] == crate::t1b::FREE };
            let ok = land && s.z[i].is_finite() && s.over[i] != crate::t1b_coast::OVER_EDIFICE;
            if self.lane_co && only.map(|c| lp.cell[i] == c).unwrap_or(true) {
                let e = &mut out.entry(lp.cell[i]).or_default().plain;
                e[0] += 1.0;
                if land && s.z[i].is_finite() { e[1] += 1.0; if s.over[i] != crate::t1b_coast::OVER_EDIFICE { e[2] += s.z[i]; e[3] += self.refv(lp, i); e[4] += lp.bil[i] as f64; e[5] += 1.0; } }
            }
            if !ok { continue; }
            for k in 0..4 {
                let w = lp.hat_w[i][k] as f64;
                if w <= 0.0 { continue; }
                let c = lp.hat_c[i][k];
                if let Some(o) = only { if c != o { continue; } }
                let a = &mut out.entry(c).or_default().hat;
                a[0] += w; a[1] += w * s.z[i]; a[2] += w * self.refv(lp, i); a[3] += w * trunk[i];
            }
        }
        out
    }

    /// LS1 for one cell, cold: the cell's land nodes start from the built network's receivers; everything else from `cur`.
    fn route(&self, cur: &Network, c: u32, cell_nodes: &[usize], table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>) -> (Network, crate::t1b_ls::LsStats) {
        let (t0, laws, p, nctx) = (self.t0, &self.extras.laws, self.p, self.nctx);
        let mut start = cur.clone();
        for &n in cell_nodes { start.recv[n] = self.net0.recv[n]; }
        crate::t1b_fix::rebuild(&mut start, &self.own, p.meander_permille);
        let (mut n2, st) = crate::t1b_ls::level_solve_s(self.net0, Some(&start), nctx, t0, laws, table, p, self.ls,
            &|nn: &Network| node_beds_core(nn, nctx, t0, laws, table, flood, p).bed, &|nn: &Network| steep_labels(nn, nctx, t0, laws, table, flood, p), false, Some(c));
        // Areas: the cell keeps its last-pass areas (B1); elsewhere the f32 accumulation of the current receivers.
        let keep: Vec<f64> = cell_nodes.iter().map(|&n| n2.area_m2[n]).collect();
        crate::t1b_fix::rebuild(&mut n2, &self.own, p.meander_permille);
        for n in 0..NODES { if n2.kind[n] == K_LAND && !n2.conduit[n] { n2.area_m2[n] = n2.area_m2[n] as f32 as f64; } }
        for (k, &n) in cell_nodes.iter().enumerate() { n2.area_m2[n] = keep[k]; }
        (n2, st)
    }

    /// Beds, curves, coast context and the surface on lattice `lp` for a network and tables.
    fn surface(&self, lp: &LatPrep, net: &Network, table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>) -> R<(Surface, crate::t1b_coast::CoastCtx)> {
        let mut beds = node_beds_core(net, self.nctx, self.t0, &self.extras.laws, table, flood, self.p);
        beds.net = Some(net.clone());
        let cv = curves(self.net0, &beds, self.skel, self.t0, self.p);
        let cx = self.co.ctx(self.p, self.t0, self.chart, &cv, table, flood)?;
        let s = evaluate(lp, net, &cv, self.t0, &self.extras.laws, table, self.p, &cx);
        Ok((s, cx))
    }
}

/// One 1-D root on an integer grid of an increasing residual r(k), from k0 within [lo, hi]. `first` is the first step
/// proposal from (k0, r0); `eval` returns the residual (cached by the caller). Returns (k, evaluations, at bound, bracketed).
fn root_grid(k0: i64, lo: i64, hi: i64, max_evals: u32, first: &dyn Fn(f64) -> i64, eval: &mut dyn FnMut(i64) -> R<f64>) -> R<(i64, u32, bool, bool)> {
    let mut seen: BTreeMap<i64, f64> = BTreeMap::new();
    let mut n = 0u32;
    let mut ev = |k: i64, seen: &mut BTreeMap<i64, f64>, n: &mut u32| -> R<f64> {
        if let Some(v) = seen.get(&k) { return Ok(*v); }
        let v = eval(k)?; *n += 1; seen.insert(k, v); Ok(v)
    };
    let best = |seen: &BTreeMap<i64, f64>| -> i64 { *seen.iter().min_by(|a, b| a.1.abs().total_cmp(&b.1.abs()).then(a.0.cmp(b.0))).unwrap().0 };
    let k0 = k0.clamp(lo, hi);
    let r0 = ev(k0, &mut seen, &mut n)?;
    if r0 == 0.0 { return Ok((k0, n, k0 == lo || k0 == hi, true)); }
    // Bracket: step against the residual's sign until it changes or a bound is met.
    let dir: i64 = if r0 > 0.0 { -1 } else { 1 };
    let mut step = (first(r0) - k0).abs().max(1);
    let (mut ka, mut ra) = (k0, r0);
    let (kb, rb);
    loop {
        let k1 = (ka + dir * step).clamp(lo, hi);
        if k1 == ka {
            // At a bound with the same sign: the bound is the answer.
            return Ok((ka, n, true, false));
        }
        let r1 = ev(k1, &mut seen, &mut n)?;
        if r1 == 0.0 { return Ok((k1, n, k1 == lo || k1 == hi, true)); }
        if (r1 > 0.0) != (ra > 0.0) { kb = k1; rb = r1; break; }
        // Secant extrapolation (or doubling) for the next step.
        let sec = if (r1 - ra).abs() > 1e-12 { (r1 / (r1 - ra) * (k1 - ka) as f64).abs().ceil() as i64 } else { 2 * step };
        step = sec.clamp(step, 8 * step).max(1);
        ka = k1; ra = r1;
        if n >= max_evals { return Ok((best(&seen), n, false, false)); }
    }
    // (lo_k, r < 0) and (hi_k, r > 0).
    // kl: the endpoint with r < 0, kh: the one with r > 0 (either order on the grid).
    let (mut kl, mut rl, mut kh, mut rh) = if ra < 0.0 { (ka, ra, kb, rb) } else { (kb, rb, ka, ra) };
    // Illinois regula falsi on the grid until the bracket is two adjacent grid points.
    let mut side = 0i32;
    while (kh - kl).abs() > 1 {
        if n >= max_evals { return Ok((best(&seen), n, false, false)); }
        let t = kl as f64 + (kh - kl) as f64 * (-rl) / (rh - rl);
        let mut kt = t.round() as i64;
        if kt <= kl.min(kh) || kt >= kl.max(kh) { kt = (kl + kh).div_euclid(2); }
        if kt == kl || kt == kh { kt = (kl + kh).div_euclid(2); }
        let rt = ev(kt, &mut seen, &mut n)?;
        if rt == 0.0 { return Ok((kt, n, false, true)); }
        let up = (rt > 0.0) == (rh > 0.0);
        if up { kh = kt; rh = rt; if side == 1 { rl *= 0.5; } side = 1; } else { kl = kt; rl = rt; if side == -1 { rh *= 0.5; } side = -1; }
    }
    // The adjacent pair: the smaller |residual| (true values, not the Illinois-scaled ones), ties to the lower k.
    let (vl, vh) = (seen[&kl], seen[&kh]);
    let (k_lo, k_hi) = (kl.min(kh), kl.max(kh));
    let (v_lo, v_hi) = if kl < kh { (vl, vh) } else { (vh, vl) };
    Ok((if v_hi.abs() < v_lo.abs() { k_hi } else { k_lo }, n, false, true))
}

#[allow(clippy::too_many_arguments)]
pub fn solve_e_cells(t0: &T0, extras: &T0Extras, skel: &SkeletonView, chart: &Chart, net: &Network, nctx: &NodeCtx, p: &RiverNetParamsV1, strata: Option<&StrataCols>,
    anchor: (f64, f64), co: &CoInputs, log: &mut dyn FnMut(String)) -> R<SolveResult> {
    let t_all = std::time::Instant::now();
    let te = crate::t1b_ls::te(p).ok_or("solve_cell needs params.ls.tb.tc.te")?.clone();
    let ls = p.ls.as_ref().filter(|l| l.solve == 1).ok_or("solve_cell needs ls_solve=1")?;
    // ---- Chart lattice, hat masses, starting table, coastal flags and flooding offsets: exactly as solve_e.
    let dxmm = p.solve_spacing_mm;
    let half = vz_tilekey::CHART_HALF_EXTENT_MM / dxmm;
    let lat = EvalLat { gx0: -half, gy0: -half, n: (2 * half + 1) as usize, dx_mm: dxmm };
    let lp = lat_prep(lat, chart, t0, extras, p, strata, anchor)?;
    let nn = lat.nodes();
    let sea = t0.sea_level_m;
    let f = dxmm / crate::t1b::T1B_SPACING_MM;
    let trunk_of = |l: &EvalLat| -> Vec<f64> { (0..l.nodes()).map(|i| { let (gx, gy) = l.g(i); crate::t1b::trunk_bilinear(skel, gx * f, gy * f, sea).max(sea) }).collect() };
    let trunk = trunk_of(&lat);
    let mut mass: BTreeMap<u32, (f64, f64)> = BTreeMap::new();
    for i in 0..nn { for k in 0..4 {
        let w = lp.hat_w[i][k] as f64;
        if w <= 0.0 { continue; }
        let e = mass.entry(lp.hat_c[i][k]).or_default();
        e.1 += w;
        if lp.water[i] == crate::t1b::FREE { e.0 += w; }
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
    let r4n = crate::t1b_ls::tc_r4n(p);
    let mut r4n_cells: Vec<u32> = Vec::new();
    if r4n > 0 {
        for (&c, &(lm, tm)) in &mass { if solved[&c] && lm < r4n as f64 / 1000.0 * tm { r4n_cells.push(c); } }
        for c in &r4n_cells { solved.insert(*c, false); }
    }
    let nb_mean = |c: u32, table: &BTreeMap<u32, f64>| -> Option<f64> {
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
    let apply_r4n = |table: &mut BTreeMap<u32, f64>| { for &c in &r4n_cells { if let Some(m) = nb_mean(c, table) { table.insert(c, m); } } };
    if r4n > 0 { apply_r4n(&mut table); }
    let (pert_ppm, pert_mode) = crate::t1b_ls::tc_pert(p);
    if pert_mode == 2 { for v in table.values_mut() { *v *= 1.0 + pert_ppm as f64 / 1e6; } }
    let lane_co = p.coast.is_some() || p.volcano.is_some() || p.f6_reference == 2;
    let cmin = p.coast.as_ref().map(|c| c.coastal_min_permille).unwrap_or(10) as f64 / 1000.0;
    let mut cnt: BTreeMap<u32, (u64, u64)> = BTreeMap::new();
    for i in 0..nn { let e = cnt.entry(lp.cell[i]).or_default(); e.0 += 1; if lp.bil[i] as f64 > sea && t0.water(lp.cell[i]) != CellWater::Lake { e.1 += 1; } }
    let coastal: BTreeMap<u32, bool> = cnt.iter().map(|(&c, &(n, l))| { let f = l as f64 / n as f64; (c, f > cmin && f < 1.0 - cmin) }).collect();
    let bil_frac: BTreeMap<u32, f64> = cnt.iter().map(|(&c, &(n, l))| (c, l as f64 / n.max(1) as f64)).collect();
    let mut flood: BTreeMap<u32, f64> = BTreeMap::new();
    if let Some(cp) = &p.coast { for (&c, &k) in &coastal { if k { flood.insert(c, cp.flood_init_mm as f64 / 1000.0); } } }
    let cx_ctx = Ctx { t0, extras, skel, chart, net0: net, nctx, p, co, ls, own: crate::t1b_fix::own_area(net), coastal: coastal.clone(), lane_co };
    // ---- Cells with land nodes, downstream-first (as LS1's order), their nodes, and the solved cells' sub-lattices.
    let is_land = |n: usize| net.kind[n] == K_LAND && !net.conduit[n];
    let mut cell_nodes: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for n in 0..NODES { if is_land(n) { cell_nodes.entry(net.t0_cell[n]).or_default().push(n); } }
    let mut order: Vec<(usize, u32)> = cell_nodes.keys().map(|&c| (t0.chain_from(c).map(|v| v.len()).unwrap_or(1), c)).collect();
    order.sort();
    let solved_cells: Vec<u32> = solved.iter().filter(|(_, &s)| s).map(|(&c, _)| c).collect();
    let mp = te.sub_margin_mm / dxmm;
    let mut subs: BTreeMap<u32, (LatPrep, Vec<f64>)> = BTreeMap::new();
    for &c in &solved_cells {
        let (mut x0, mut x1, mut y0, mut y1) = (i64::MAX, i64::MIN, i64::MAX, i64::MIN);
        for i in 0..nn {
            if lp.cell[i] != c && !(0..4).any(|k| lp.hat_c[i][k] == c && lp.hat_w[i][k] > 0.0) { continue; }
            let (gx, gy) = lat.g(i);
            x0 = x0.min(gx); x1 = x1.max(gx); y0 = y0.min(gy); y1 = y1.max(gy);
        }
        if x0 > x1 { continue; }
        let n_sub = ((x1 - x0).max(y1 - y0) + 2 * mp + 1).min(lat.n as i64);
        let clampo = |v: i64| v.clamp(lat.gx0, lat.gx0 + lat.n as i64 - n_sub);
        let sl = EvalLat { gx0: clampo(x0 - mp), gy0: clampo(y0 - mp), n: n_sub as usize, dx_mm: dxmm };
        let slp = lat_prep(sl, chart, t0, extras, p, strata, anchor)?;
        let st = trunk_of(&sl);
        log(format!("  lane TE: cell {c} sub-lattice gx0 {} gy0 {} n {} ({:.1} km)", sl.gx0, sl.gy0, sl.n, sl.n as f64 * dxmm as f64 / 1e6));
        subs.insert(c, (slp, st));
    }
    log(format!("  lane TE solve_cell: order {:?}; solved {:?}; r4n {:?}; coastal solved {:?}", order.iter().map(|o| o.1).collect::<Vec<_>>(), solved_cells, r4n_cells,
        solved_cells.iter().filter(|c| coastal.get(c).copied().unwrap_or(false)).collect::<Vec<_>>()));
    // ---- Grids.
    let qe = te.e_grid_ppm as f64 / 1e6;
    let qf = te.f_grid_mm as f64 / 1000.0;
    let (ke_lo, ke_hi) = ((emin.max(1e-9).ln() / qe).ceil() as i64, (emax.ln() / qe).floor() as i64);
    let kf_hi = p.coast.as_ref().map(|c| (c.flood_max_mm as f64 / 1000.0 / qf).floor() as i64).unwrap_or(0);
    let tol = p.solve_tolerance_mm as f64 / 1000.0;
    let mut cur = net.clone();
    let mut evals = 0u32;
    let mut ls_history: Vec<[u32; 3]> = Vec::new();
    let mut hist: Vec<f64> = Vec::new();
    let mut visits: Vec<Value> = Vec::new();
    let mut last_meas: BTreeMap<u32, Meas> = BTreeMap::new();
    let mut fixed_point = false;
    let mut sweeps = 0u32;
    for sweep in 0..te.sweeps_max {
        sweeps += 1;
        let mut changed_cells: Vec<u32> = Vec::new();
        let mut maxd: f64 = 0.0;
        for &(_, c) in &order {
            let tv = std::time::Instant::now();
            let nodes = &cell_nodes[&c];
            if !solved.get(&c).copied().unwrap_or(false) || !subs.contains_key(&c) {
                let (n2, st) = cx_ctx.route(&cur, c, nodes, &table, &flood);
                let ch = nodes.iter().filter(|&&n| n2.recv[n] != cur.recv[n]).count() as u32;
                ls_history.push([st.passes, *st.changed_per_pass.last().unwrap_or(&0), ch]);
                if ch > 0 { changed_cells.push(c); }
                cur = n2;
                continue;
            }
            let (slp, strk) = &subs[&c];
            let is_coastal = p.coast.is_some() && coastal.get(&c).copied().unwrap_or(false);
            let mut ke = (table[&c].max(1e-9).ln() / qe).round() as i64;
            let mut kf = if is_coastal { (flood.get(&c).copied().unwrap_or(0.0) / qf).round() as i64 } else { 0 };
            let (ke_start, kf_start) = (ke, kf);
            // Evaluations of this visit, cached by grid point.
            let mut cache: BTreeMap<(i64, i64), (Meas, Network)> = BTreeMap::new();
            let mut nev = 0u32;
            let eval_at = |ke: i64, kf: i64, cache: &mut BTreeMap<(i64, i64), (Meas, Network)>, nev: &mut u32| -> R<Meas> {
                if let Some((m, _)) = cache.get(&(ke, kf)) { return Ok(*m); }
                let mut tab = table.clone();
                tab.insert(c, (ke as f64 * qe).exp());
                if r4n > 0 { apply_r4n(&mut tab); }
                let mut fl = flood.clone();
                if is_coastal { fl.insert(c, kf as f64 * qf); }
                let (n2, _) = cx_ctx.route(&cur, c, nodes, &tab, &fl);
                let (s, _) = cx_ctx.surface(slp, &n2, &tab, &fl)?;
                let m = cx_ctx.measure(slp, strk, &s, Some(c)).get(&c).copied().unwrap_or_default();
                *nev += 1;
                cache.insert((ke, kf), (m, n2));
                Ok(m)
            };
            let (mut at_bound, mut bracketed_e, mut bracketed_f) = (false, false, !is_coastal);
            let mut alts = 0u32;
            for alt in 0..te.alternations_max {
                alts += 1;
                let (ke_prev, kf_prev) = (ke, kf);
                // E at this F: root of M - T (increasing in ln E).
                {
                    let m0 = eval_at(ke, kf, &mut cache, &mut nev)?;
                    let first = |r0: f64| -> i64 {
                        let (num, den) = (m0.t() - m0.l(), m0.m() - m0.l());
                        let st = if num > 1.0 && den > 1.0 { 0.6 * (num / den).ln() } else if r0 < 0.0 { 2f64.ln() } else { -(2f64.ln()) };
                        ke + (st.clamp(-(2f64.ln()), 2f64.ln()) / qe).round() as i64
                    };
                    let mut ev = |k: i64| -> R<f64> { let m = eval_at(k, kf, &mut cache, &mut nev)?; Ok(m.m() - m.t()) };
                    let (k, _, b, br) = root_grid(ke, ke_lo, ke_hi, te.root_evals_max, &first, &mut ev)?;
                    ke = k; at_bound = b; bracketed_e = br;
                }
                // F at this E (coastal cells): root of bilinear land fraction - land fraction (increasing in F).
                if is_coastal {
                    let cp = p.coast.as_ref().unwrap();
                    let bf = bil_frac.get(&c).copied().unwrap_or(0.0);
                    let m0 = eval_at(ke, kf, &mut cache, &mut nev)?;
                    let first = |_: f64| -> i64 {
                        let ex = m0.land_fraction() - bf;
                        let st = (cp.flood_gain_mm as f64 / 1000.0 * ex).clamp(-(cp.flood_step_max_mm as f64) / 1000.0, cp.flood_step_max_mm as f64 / 1000.0);
                        kf + (st / qf).round() as i64
                    };
                    let mut ev = |k: i64| -> R<f64> { let m = eval_at(ke, k, &mut cache, &mut nev)?; Ok(bf - m.land_fraction()) };
                    let (k, _, _, br) = root_grid(kf, 0, kf_hi, te.root_evals_max, &first, &mut ev)?;
                    kf = k; bracketed_f = br;
                }
                // An alternation that moves neither unknown ends the visit.
                let _ = alt;
                if ke == ke_prev && kf == kf_prev { break; }
            }
            // Commit the chosen grid point (its evaluation is cached).
            let m = eval_at(ke, kf, &mut cache, &mut nev)?;
            let (_, n2) = cache.remove(&(ke, kf)).unwrap();
            let ch = nodes.iter().filter(|&&n| n2.recv[n] != cur.recv[n]).count() as u32;
            let moved = ke != ke_start || kf != kf_start;
            table.insert(c, (ke as f64 * qe).exp());
            if r4n > 0 { apply_r4n(&mut table); }
            if is_coastal { flood.insert(c, kf as f64 * qf); }
            cur = n2;
            evals += nev;
            ls_history.push([0, nev, ch]);
            if moved || ch > 0 { changed_cells.push(c); }
            let d = m.m() - m.t();
            maxd = maxd.max(d.abs());
            last_meas.insert(c, m);
            let lf = m.land_fraction() - bil_frac.get(&c).copied().unwrap_or(0.0);
            log(format!("  lane TE sweep {sweep} cell {c}: evals {nev} alternations {alts} E {:.2} (k {ke}, from {ke_start}{}) F {:.0} (k {kf}, from {kf_start}) M - T {:.2} m land - bil {:.3} bracketed E {} F {} receivers changed {ch} ({:.1} s)",
                table[&c], if at_bound { ", BOUND" } else { "" }, flood.get(&c).copied().unwrap_or(0.0), d, lf, bracketed_e, bracketed_f, tv.elapsed().as_secs_f64()));
            visits.push(json!({"sweep": sweep, "cell": c, "evals": nev, "alternations": alts, "k_e": ke, "k_f": kf, "e_m_per_my": table[&c], "f_m": flood.get(&c).copied().unwrap_or(0.0),
                "diff_m": d, "land_minus_bilinear": lf, "at_bound": at_bound, "bracketed_e": bracketed_e, "bracketed_f": bracketed_f, "receivers_changed": ch, "seconds": tv.elapsed().as_secs_f64()}));
        }
        hist.push(maxd);
        log(format!("  lane TE sweep {sweep}: max |M - T| (sub-lattices) {maxd:.2} m; cells changed {:?}; evaluations so far {evals}; {:.1} s", changed_cells, t_all.elapsed().as_secs_f64()));
        if changed_cells.is_empty() { fixed_point = true; break; }
    }
    // Lane TC C4 pert_mode 1 (intrinsic sensitivity): the solved table scaled, every cell re-routed once on it.
    if pert_mode == 1 {
        for v in table.values_mut() { *v *= 1.0 + pert_ppm as f64 / 1e6; }
        for &(_, c) in &order { let (n2, _) = cx_ctx.route(&cur, c, &cell_nodes[&c], &table, &flood); cur = n2; }
        log(format!("  lane TC/TE: published table x (1 + {pert_ppm} ppm), every cell re-routed once"));
    }
    // ---- One full-chart evaluation of the published state for the rows and the convergence flag.
    let (s, cx) = cx_ctx.surface(&lp, &cur, &table, &flood)?;
    let meas = cx_ctx.measure(&lp, &trunk, &s, None);
    let mut rows = Vec::new();
    let mut maxd: f64 = 0.0;
    let mut sub_vs_full: Vec<Value> = Vec::new();
    for (&c, &sv) in &solved {
        let mm = meas.get(&c).copied().unwrap_or_default();
        let (m, tg, l) = (mm.m(), mm.t(), mm.l());
        let diff = m - tg;
        let (lm, tm) = mass[&c];
        let mut row = CellRow { cell: c, solved: sv, land_hat_fraction: lm / tm, e0_m_per_my: e0_table[&c], e_m_per_my: table[&c], hat_mean_z_m: m, hat_mean_bilinear_m: tg, hat_mean_trunk_m: l, diff_m: diff,
            at_bound: table[&c] <= emin || table[&c] >= emax, coastal: None, flood_offset_m: None, land_fraction: None, land_fraction_bilinear: None, cell_mean_z_m: None, cell_mean_ref_m: None,
            cell_mean_bilinear_m: None, cell_land_posts: None, r4_plain_reference: None, r4n_neighbour_mean: if r4n_cells.contains(&c) { Some(true) } else { None }, cell_mean_trunk_m: None, hat_diff_m: None, es_plain_target: None };
        if lane_co {
            let pl = mm.plain;
            row.coastal = Some(coastal.get(&c).copied().unwrap_or(false));
            row.flood_offset_m = flood.get(&c).copied().or(Some(0.0));
            row.land_fraction = Some(if pl[0] > 0.0 { pl[1] / pl[0] } else { 0.0 });
            row.land_fraction_bilinear = Some(bil_frac.get(&c).copied().unwrap_or(0.0));
            row.cell_mean_z_m = Some(if pl[5] > 0.0 { pl[2] / pl[5] } else { 0.0 });
            row.cell_mean_ref_m = Some(if pl[5] > 0.0 { pl[3] / pl[5] } else { 0.0 });
            row.cell_mean_bilinear_m = Some(if pl[5] > 0.0 { pl[4] / pl[5] } else { 0.0 });
            row.cell_land_posts = Some(pl[5] as u64);
        }
        if sv {
            maxd = maxd.max(diff.abs());
            if let Some(lm) = last_meas.get(&c) { sub_vs_full.push(json!({"cell": c, "sub_m": lm.m(), "full_m": m, "sub_land_fraction": lm.land_fraction(), "full_land_fraction": row.land_fraction})); }
        }
        rows.push(row);
    }
    let maxf4 = rows.iter().filter(|r| r.coastal == Some(true) && r.solved).map(|r| (r.land_fraction.unwrap_or(0.0) - r.land_fraction_bilinear.unwrap_or(0.0)).abs()).fold(0.0, f64::max);
    let f4_ok = p.coast.as_ref().map(|cp| maxf4 <= cp.flood_tolerance_permille as f64 / 1000.0).unwrap_or(true);
    let converged = maxd <= tol && f4_ok;
    hist.push(maxd);
    for r in rows.iter().filter(|r| r.solved) { log(format!("  lane TE final (full chart) cell {} E {:.2} F {:.1} M {:.1} T {:.1} diff {:.2} land {:.3} bil {:.3}", r.cell, r.e_m_per_my, r.flood_offset_m.unwrap_or(0.0), r.hat_mean_z_m, r.hat_mean_bilinear_m, r.diff_m, r.land_fraction.unwrap_or(0.0), r.land_fraction_bilinear.unwrap_or(0.0))); }
    let secs = t_all.elapsed().as_secs_f64();
    log(format!("  lane TE solve_cell: sweeps {sweeps} fixed point {fixed_point} evaluations {evals} max |M - T| (full chart) {maxd:.2} m F4 {maxf4:.3} converged {converged} ({secs:.1} s)"));
    let te_rec = json!({"rule": "solve_cell_v1 (params.ls.tb.tc.te.cell_rule)", "sweeps": sweeps, "fixed_point": fixed_point, "evaluations": evals, "seconds": secs,
        "order": order.iter().map(|o| o.1).collect::<Vec<_>>(), "solved": solved_cells, "r4n_cells": r4n_cells, "visits": visits, "sub_vs_full": sub_vs_full,
        "note": "EXPLORATORY (lane TE PREREG): history_max_abs_diff_m = per sweep max |M - T| on the sub-lattices, last entry the full-chart value; the published network is the cell-solved network (no LS1 re-run at publish)"});
    Ok(SolveResult { table, rows, history_max_abs_diff_m: hist, converged, evals, seconds: secs, flood, cx, history_max_abs_f4: vec![maxf4], net: Some(cur), ls_history, frozen: true, te: Some(te_rec), es: None, ped: None })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn root_grid_ends_on_the_same_adjacent_pair_from_any_start() {
        // r(k) = k - 37.3: the adjacent pair is (37, 38); 37 has the smaller |r|. Starts on both sides, near and far.
        for &k0 in &[-200i64, 0, 36, 37, 38, 39, 500] {
            let mut ev = |k: i64| -> R<f64> { Ok(k as f64 - 37.3) };
            let (k, n, b, br) = root_grid(k0, -1000, 1000, 40, &|r0| k0 - (r0 * 0.5) as i64, &mut ev).unwrap();
            assert_eq!(k, 37, "from {k0}");
            assert!(br && !b && n <= 40);
        }
        // No root inside the bounds: the bound on the residual's side.
        let mut ev = |k: i64| -> R<f64> { Ok(k as f64 + 2000.0) };
        let (k, _, b, br) = root_grid(5, -100, 100, 40, &|_| 0, &mut ev).unwrap();
        assert!(k == -100 && b && !br);
    }
}
