//! SEED-RANGES (lane RG, 2026-10-13): `contact-length-indentation-taper-v1`
//! (post-collision indentation) and `plateau-cap-spill-v1` (orogenic
//! plateau growth). Opt-in law versions; pre-registered in
//! `unpublished-design-record`.
//!
//! Why (design review review of cycle_v1, 2026-10-13): under
//! `contact-length-pole-share-v1` the smaller group takes the larger
//! group's pole the step the contact reaches the declared length, so
//! convergence stops dead; with the belt's 70 km backstop the belt is a one
//! or two cell spike that n = 1 erosion removes in ~50 My. Earth's ranges
//! persist because convergence continues for 30-50 My after collision and
//! the thickened crust spreads sideways into a plateau.
//!
//! `contact-length-indentation-taper-v1` (`suture_indent_pass`): suturing
//! is exactly `conserve::suture_pass` (same threshold, winner rule and
//! union-find). For every new suture the loser group's member plates (a
//! recorded set) keep `omega_rel = omega(loser root) - omega(winner root)`
//! from just before the suture, scaled by `f = 1 - tau / T_ind`: every
//! step each member's pole is its group root's pole plus the sum of
//! `f * omega_rel` over its active indentations. An indentation ends at
//! `tau >= T_ind`, when the cumulative shortening (length-weighted mean
//! converging normal rate over 4-neighbour continental edges between
//! members and the rest of the group, times dt) reaches the budget, or when
//! no member is left in the winner's group (breakup). Group roots are never
//! modified, so re-solves and housekeeping act exactly as before.
//!
//! `plateau-cap-spill-v1` (`plateau_spill`): after the advance, each
//! continental cell thicker than the cap (ascending cell order) keeps the
//! cap and gives the excess to the nearest same-plate continental cells
//! below the cap: rings of edge neighbours (breadth-first through cells
//! already at the cap; ascending cell id within a ring), each frontier
//! cell filled up to the cap, within the reach (great-circle from the
//! source). Excess that finds no room is delaminated (counted). Oceanic
//! cells above the cap lose the excess (delaminated, the old backstop
//! role). Volume moves between cells with the lineage field
//! `shortening_added_m` and the step's thickening.
//!
//! Determinism: fixed iteration orders (plates in record order, events in
//! suture order, cells ascending), BTreeMap/BTreeSet only, no threads.

use crate::boundaries::{angular_velocity, normal_rate_at};
use crate::conserve::{suture_pass, SutureEvent, SutureState};
use crate::kinematic::{cell_areas, centroids, CrustState, KinematicError};
use crate::rotation::{great_circle_rad, unit_to_lat_lon};
use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;
use std::collections::{BTreeMap, BTreeSet};

/// Constants of `contact-length-indentation-taper-v1` beyond the contact
/// length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndentParams {
    pub time_my: f64,
    pub budget_m: f64,
}

/// End causes of an indentation.
pub const INDENT_ACTIVE: u8 = 0;
pub const INDENT_END_TIME: u8 = 1;
pub const INDENT_END_BUDGET: u8 = 2;
pub const INDENT_END_DETACHED: u8 = 3;

/// One indentation (run state).
#[derive(Debug, Clone, PartialEq)]
pub struct Indentation {
    pub t0_my: f64,
    /// Winner and loser group roots at the suture.
    pub winner: PlateId,
    pub loser: PlateId,
    /// The loser group's member plates at the suture (ascending).
    pub members: Vec<PlateId>,
    /// Relative angular velocity (rad/My) loser - winner at the suture.
    pub omega_rel: [f64; 3],
    /// Cumulative shortening at the contact (m).
    pub shortening_m: f64,
    /// Mean converging rate (m/My) measured in the first step.
    pub v0_m_per_my: Option<f64>,
    pub end_t_my: Option<f64>,
    pub end_cause: u8,
}

/// Per-step bookkeeping of the indentation law.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct IndentStepReport {
    pub started: u32,
    pub active: u32,
    pub ended: u32,
    /// Contact length (m) and length-weighted mean converging rate (m/My)
    /// summed / averaged over the active indentations.
    pub contact_m: f64,
    pub mean_rate_m_per_my: f64,
    pub max_taper: f64,
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Euler pole of an angular velocity (rad/My); zero -> rate 0 at the north
/// pole (the convention of `reorg::pole_of`).
fn pole_of(w: [f64; 3]) -> Result<EulerPoleV0, KinematicError> {
    let mag = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    if mag > 1.0e-15 {
        let (lat, lon) = unit_to_lat_lon(scale(w, 1.0 / mag))?;
        Ok(EulerPoleV0 { lat_deg: lat, lon_deg: lon, rate_deg_per_my: mag.to_degrees() })
    } else {
        Ok(EulerPoleV0 { lat_deg: 90.0, lon_deg: 0.0, rate_deg_per_my: 0.0 })
    }
}

/// `contact-length-indentation-taper-v1`: `suture_pass` (sutures and root
/// poles), then new indentations for the new sutures, then the tapered
/// relative rotation on the member plates, then the shortening measured on
/// the step-start state over `dt_my`. Returns the suture events, the
/// largest contact (m) and the step report.
#[allow(clippy::too_many_arguments)]
pub fn suture_indent_pass(
    grid: &GridSpec,
    state: &CrustState,
    plates: &mut [PlateRecordV0],
    sut: &mut SutureState,
    threshold_m: f64,
    ip: &IndentParams,
    indents: &mut Vec<Indentation>,
    dt_my: f64,
) -> Result<(Vec<SutureEvent>, f64, IndentStepReport), KinematicError> {
    let mut rep = IndentStepReport::default();
    // Angular velocities and groups just before the suture pass.
    let mut before: BTreeMap<PlateId, [f64; 3]> = BTreeMap::new();
    for p in plates.iter() {
        if let Some(e) = p.euler_pole {
            before.insert(p.plate_id, angular_velocity(&e)?);
        }
    }
    let sut_before = sut.clone();
    let (events, max_contact) = suture_pass(grid, state, plates, sut, threshold_m)?;
    let t = state.t_my;
    for e in &events {
        let (Some(&wl), Some(&ww)) = (before.get(&e.loser_root), before.get(&e.winner_root)) else {
            continue;
        };
        let mut members: Vec<PlateId> = plates
            .iter()
            .filter(|p| p.euler_pole.is_some() && sut_before.find(p.plate_id) == e.loser_root)
            .map(|p| p.plate_id)
            .collect();
        members.sort();
        members.dedup();
        if members.is_empty() {
            continue;
        }
        indents.push(Indentation {
            t0_my: t,
            winner: e.winner_root,
            loser: e.loser_root,
            members,
            omega_rel: sub(wl, ww),
            shortening_m: 0.0,
            v0_m_per_my: None,
            end_t_my: None,
            end_cause: INDENT_ACTIVE,
        });
        rep.started += 1;
    }
    // Which plates exist (with a pole) now.
    let live: BTreeSet<PlateId> = plates
        .iter()
        .filter(|p| p.euler_pole.is_some())
        .map(|p| p.plate_id)
        .collect();
    // Active members per indentation and the summed tapered rotation.
    let mut extra: BTreeMap<PlateId, [f64; 3]> = BTreeMap::new();
    let mut active_sets: Vec<(usize, BTreeSet<PlateId>, PlateId)> = Vec::new();
    for (ii, ind) in indents.iter_mut().enumerate() {
        if ind.end_cause != INDENT_ACTIVE {
            continue;
        }
        let tau = t - ind.t0_my;
        let f = 1.0 - tau / ip.time_my;
        if !(f > 0.0) {
            ind.end_cause = INDENT_END_TIME;
            ind.end_t_my = Some(t);
            rep.ended += 1;
            continue;
        }
        let anchor = sut.find(ind.winner);
        let set: BTreeSet<PlateId> = ind
            .members
            .iter()
            .cloned()
            .filter(|m| live.contains(m) && sut.find(*m) == anchor && sut.find(*m) != *m)
            .collect();
        if set.is_empty() {
            ind.end_cause = INDENT_END_DETACHED;
            ind.end_t_my = Some(t);
            rep.ended += 1;
            continue;
        }
        let w = scale(ind.omega_rel, f);
        for m in &set {
            let e = extra.entry(*m).or_insert([0.0; 3]);
            *e = add(*e, w);
        }
        rep.max_taper = rep.max_taper.max(f);
        active_sets.push((ii, set, anchor));
    }
    for p in plates.iter_mut() {
        if let (Some(w), Some(e)) = (extra.get(&p.plate_id), p.euler_pole) {
            p.euler_pole = Some(pole_of(add(angular_velocity(&e)?, *w))?);
        }
    }
    rep.active = active_sets.len() as u32;
    if active_sets.is_empty() {
        return Ok((events, max_contact, rep));
    }
    // Shortening at each active contact, on the step-start state with the
    // poles just set.
    let poles: BTreeMap<PlateId, EulerPoleV0> = plates
        .iter()
        .filter_map(|p| p.euler_pole.map(|e| (p.plate_id, e)))
        .collect();
    let cents = centroids(grid)?;
    let radius = grid.radius_m();
    let dlat = std::f64::consts::PI / grid.nlat() as f64;
    let dlon = 2.0 * std::f64::consts::PI / grid.nlon() as f64;
    let nlon = grid.nlon();
    // Group root per plate (one find per plate id).
    let mut root_of: BTreeMap<PlateId, PlateId> = BTreeMap::new();
    for p in plates.iter() {
        root_of.insert(p.plate_id, sut.find(p.plate_id));
    }
    let (mut tot_l, mut tot_vl) = (0.0, 0.0);
    for (ii, set, anchor) in &active_sets {
        let (mut l_sum, mut vl_sum) = (0.0f64, 0.0f64);
        for c in 0..state.cells() {
            if state.kind[c] != CrustKindCodeV0::Continental || !set.contains(&state.plate[c]) {
                continue;
            }
            let (i, j) = (c / nlon, c % nlon);
            let lat_bot = -std::f64::consts::FRAC_PI_2 + i as f64 * dlat;
            for (ti, tj, nb) in grid
                .neighbors(i, j)
                .map_err(|e| KinematicError::Grid(e.to_string()))?
            {
                if nb.di.abs() + nb.dj.abs() != 1 {
                    continue;
                }
                let k = ti * nlon + tj;
                let pk = state.plate[k];
                if state.kind[k] != CrustKindCodeV0::Continental
                    || set.contains(&pk)
                    || root_of.get(&pk) != Some(anchor)
                {
                    continue;
                }
                let (Some(pa), Some(pb)) = (poles.get(&state.plate[c]), poles.get(&pk)) else {
                    continue;
                };
                let raw = add(cents[c], cents[k]);
                let norm = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
                if norm < 1e-12 {
                    continue;
                }
                let mid = scale(raw, 1.0 / norm);
                let vn = normal_rate_at(mid, cents[c], cents[k], pa, pb)? * radius;
                if !(vn > 0.0) {
                    continue;
                }
                let edge = if nb.di == 0 {
                    radius * dlat
                } else if nb.di > 0 {
                    radius * (lat_bot + dlat).cos() * dlon
                } else {
                    radius * lat_bot.cos() * dlon
                };
                l_sum += edge;
                vl_sum += vn * edge;
            }
        }
        let ind = &mut indents[*ii];
        let mean = if l_sum > 0.0 { vl_sum / l_sum } else { 0.0 };
        if ind.v0_m_per_my.is_none() {
            ind.v0_m_per_my = Some(mean);
        }
        ind.shortening_m += mean * dt_my;
        tot_l += l_sum;
        tot_vl += vl_sum;
        if ind.shortening_m >= ip.budget_m {
            ind.end_cause = INDENT_END_BUDGET;
            ind.end_t_my = Some(t + dt_my);
            rep.ended += 1;
        }
    }
    rep.contact_m = tot_l;
    rep.mean_rate_m_per_my = if tot_l > 0.0 { tot_vl / tot_l } else { 0.0 };
    Ok((events, max_contact, rep))
}

/// Constants of `plateau-cap-spill-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateauParams {
    pub cap_m: f64,
    pub reach_m: f64,
}

/// Per-step bookkeeping of `plateau_spill` (m^3 unless stated).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PlateauReport {
    /// Continental cells above the cap and the excess they gave away.
    pub source_cells: u64,
    pub excess_m3: f64,
    /// Volume placed on hinterland cells and the cells that received it.
    pub moved_m3: f64,
    pub filled_cells: u64,
    /// Continental excess with no room within the reach (delaminated).
    pub delaminated_m3: f64,
    /// Oceanic crust above the cap (delaminated).
    pub oceanic_delaminated_m3: f64,
    /// Largest source-to-receiver distance (m).
    pub max_distance_m: f64,
    /// Continental cells at the cap after the pass.
    pub cells_at_cap: u64,
}

/// `plateau-cap-spill-v1` on an advanced state (see the module docs).
/// `thickening` (m per cell, the step's tectonic thickening) moves with the
/// crust. Volume conserving apart from the reported delamination.
pub fn plateau_spill(
    grid: &GridSpec,
    state: &mut CrustState,
    thickening: &mut [f64],
    p: &PlateauParams,
) -> Result<PlateauReport, KinematicError> {
    let n = state.cells();
    let areas = cell_areas(grid)?;
    let cents = centroids(grid)?;
    let radius = grid.radius_m();
    let nlon = grid.nlon();
    let cap = p.cap_m;
    let mut rep = PlateauReport::default();
    // Edge neighbours, ascending.
    let mut nbrs: Vec<[u32; 4]> = Vec::with_capacity(n);
    let mut nbn: Vec<u8> = Vec::with_capacity(n);
    for c in 0..n {
        let (i, j) = (c / nlon, c % nlon);
        let mut v = [u32::MAX; 4];
        let mut m = 0usize;
        let mut list: Vec<usize> = grid
            .neighbors(i, j)
            .map_err(|e| KinematicError::Grid(e.to_string()))?
            .into_iter()
            .filter(|(_, _, nb)| nb.di.abs() + nb.dj.abs() == 1)
            .map(|(ti, tj, _)| ti * nlon + tj)
            .collect();
        list.sort_unstable();
        list.dedup();
        for k in list.into_iter().take(4) {
            v[m] = k as u32;
            m += 1;
        }
        nbrs.push(v);
        nbn.push(m as u8);
    }
    let mut stamp = vec![0u32; n];
    let mut gen = 0u32;
    for c in 0..n {
        let h = state.thickness_m[c];
        if !(h > cap + 1.0e-6) {
            continue;
        }
        let excess = h - cap;
        state.thickness_m[c] = cap;
        state.shortening_added_m[c] -= excess;
        thickening[c] -= excess;
        let mut vol = excess * areas[c];
        if state.kind[c] != CrustKindCodeV0::Continental {
            rep.oceanic_delaminated_m3 += vol;
            continue;
        }
        rep.source_cells += 1;
        rep.excess_m3 += vol;
        let plate = state.plate[c];
        gen += 1;
        stamp[c] = gen;
        let mut ring: Vec<usize> = vec![c];
        while vol > 0.0 && !ring.is_empty() {
            let mut next: Vec<usize> = Vec::new();
            for &r in &ring {
                for x in 0..nbn[r] as usize {
                    let k = nbrs[r][x] as usize;
                    if stamp[k] == gen {
                        continue;
                    }
                    stamp[k] = gen;
                    if state.kind[k] != CrustKindCodeV0::Continental || state.plate[k] != plate {
                        continue;
                    }
                    let d = great_circle_rad(cents[c], cents[k])? * radius;
                    if d > p.reach_m {
                        continue;
                    }
                    next.push(k);
                }
            }
            next.sort_unstable();
            for &k in &next {
                if vol <= 0.0 {
                    break;
                }
                let hk = state.thickness_m[k];
                if !(hk < cap) {
                    continue;
                }
                let room = (cap - hk) * areas[k];
                let (dh, take) = if room <= vol {
                    (cap - hk, room)
                } else {
                    (vol / areas[k], vol)
                };
                if dh == cap - hk {
                    state.thickness_m[k] = cap;
                } else {
                    state.thickness_m[k] = hk + dh;
                }
                state.shortening_added_m[k] += dh;
                thickening[k] += dh;
                vol -= take;
                rep.moved_m3 += take;
                rep.filled_cells += 1;
                let d = great_circle_rad(cents[c], cents[k])? * radius;
                rep.max_distance_m = rep.max_distance_m.max(d);
            }
            if vol <= 0.0 {
                break;
            }
            // Continue only through cells now at the cap.
            ring = next.into_iter().filter(|&k| !(state.thickness_m[k] < cap)).collect();
        }
        if vol > 0.0 {
            rep.delaminated_m3 += vol;
        }
    }
    for c in 0..n {
        if state.kind[c] == CrustKindCodeV0::Continental && !(state.thickness_m[c] < cap) {
            rep.cells_at_cap += 1;
        }
    }
    Ok(rep)
}

/// SEED-WORLD-COMPOSE (lane CP, 2026-10-13): constants of
/// `plateau-cap-spill-v2` (see laws.rs; run inside the belt step).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateauInBeltParams {
    pub cap_m: f64,
    pub reach_m: f64,
    /// Largest thickness step per metre of centre distance between
    /// edge-neighbouring same-plate continental cells.
    pub foothill_gradient: f64,
    /// The taper acts only while the thicker cell of a pair is above this.
    pub taper_floor_m: f64,
}

/// Fixed by the law version: taper sweeps (Gauss-Seidel, ascending cells).
pub const TAPER_SWEEPS: usize = 64;

/// Per-step bookkeeping of `plateau_in_belt` (m^3 unless stated).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PlateauBeltReport {
    /// Sediment accreted into the crust of columns whose crust plus
    /// sediment exceeded the cap.
    pub sediment_accreted_m3: f64,
    pub accreted_cells: u64,
    /// As `PlateauReport`.
    pub source_cells: u64,
    pub excess_m3: f64,
    pub moved_m3: f64,
    pub filled_cells: u64,
    pub delaminated_m3: f64,
    pub oceanic_delaminated_m3: f64,
    pub max_distance_m: f64,
    /// Foothill taper: volume moved down-step, pairs adjusted, sweeps used
    /// and the largest remaining excess step (m) after the last sweep.
    pub taper_moved_m3: f64,
    pub taper_pairs: u64,
    pub taper_sweeps: u64,
    pub taper_residual_m: f64,
    pub cells_at_cap: u64,
}

/// Edge (4-) neighbours of every cell, ascending, deduplicated.
fn edge_nbrs(grid: &GridSpec) -> Result<Vec<Vec<usize>>, KinematicError> {
    // SEED-PERF (lane SP, 2026-10-13): the edge neighbours' ids directly
    // (the `neighbors` distances were unused), per cell in parallel; same
    // sorted, deduplicated lists.
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    Ok(crate::kinematic::par::map_indexed(grid.cells(), |c| {
        let (i, j) = ((c / nlon) as i64, (c % nlon) as i64);
        let mut list: Vec<usize> = Vec::with_capacity(4);
        for (di, dj) in [(-1i64, 0i64), (0, -1), (0, 1), (1, 0)] {
            let ti = i + di;
            if ti < 0 || ti >= nlat as i64 {
                continue;
            }
            list.push((ti * nlon as i64 + (j + dj).rem_euclid(nlon as i64)) as usize);
        }
        list.sort_unstable();
        list.dedup();
        list.retain(|&k| k != c);
        list
    }))
}

/// `plateau-cap-spill-v2` on a thickness / sediment field (inside the belt
/// step, before the channel-flow spreading): (1) sediment of continental
/// columns whose crust plus sediment exceeds the cap is accreted into the
/// crust; (2) crust above the cap spills exactly as `plateau_spill` (same
/// rings, order and reach; oceanic excess delaminated); (3) foothill taper.
/// Volume conserving apart from the reported delamination; the caller does
/// the lineage (its net thickness change) and the ledger.
#[allow(clippy::too_many_arguments)]
pub fn plateau_in_belt(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    kind: &[CrustKindCodeV0],
    plate: &[PlateId],
    thickness_m: &mut [f64],
    sediment_m: &mut [f64],
    p: &PlateauInBeltParams,
) -> Result<PlateauBeltReport, KinematicError> {
    let n = grid.cells();
    if cents.len() != n || areas.len() != n || kind.len() != n || plate.len() != n || thickness_m.len() != n || sediment_m.len() != n {
        return Err(KinematicError::Grid("plateau_in_belt: field sizes differ from the grid".into()));
    }
    if !(p.cap_m > 0.0 && p.reach_m >= 0.0 && p.foothill_gradient > 0.0 && p.taper_floor_m > 0.0 && p.taper_floor_m < p.cap_m) {
        return Err(KinematicError::Grid(format!("plateau-cap-spill-v2: invalid constants {p:?}")));
    }
    let radius = grid.radius_m();
    let cap = p.cap_m;
    let cont = |c: usize| kind[c] == CrustKindCodeV0::Continental;
    let mut rep = PlateauBeltReport::default();
    // (1) Sediment on over-cap columns joins the crust.
    for c in 0..n {
        if cont(c) && sediment_m[c] > 0.0 && thickness_m[c] + sediment_m[c] > cap {
            thickness_m[c] += sediment_m[c];
            rep.sediment_accreted_m3 += sediment_m[c] * areas[c];
            rep.accreted_cells += 1;
            sediment_m[c] = 0.0;
        }
    }
    // (2) Cap and spill.
    let nbrs = edge_nbrs(grid)?;
    let mut stamp = vec![0u32; n];
    let mut gen = 0u32;
    for c in 0..n {
        let h = thickness_m[c];
        if !(h > cap + 1.0e-6) {
            continue;
        }
        let excess = h - cap;
        thickness_m[c] = cap;
        let mut vol = excess * areas[c];
        if !cont(c) {
            rep.oceanic_delaminated_m3 += vol;
            continue;
        }
        rep.source_cells += 1;
        rep.excess_m3 += vol;
        let pl = plate[c];
        gen += 1;
        stamp[c] = gen;
        let mut ring: Vec<usize> = vec![c];
        while vol > 0.0 && !ring.is_empty() {
            let mut next: Vec<usize> = Vec::new();
            for &r in &ring {
                for &k in &nbrs[r] {
                    if stamp[k] == gen {
                        continue;
                    }
                    stamp[k] = gen;
                    if !cont(k) || plate[k] != pl {
                        continue;
                    }
                    let d = great_circle_rad(cents[c], cents[k])? * radius;
                    if d > p.reach_m {
                        continue;
                    }
                    next.push(k);
                }
            }
            next.sort_unstable();
            for &k in &next {
                if vol <= 0.0 {
                    break;
                }
                let hk = thickness_m[k];
                if !(hk < cap) {
                    continue;
                }
                let room = (cap - hk) * areas[k];
                if room <= vol {
                    thickness_m[k] = cap;
                    vol -= room;
                    rep.moved_m3 += room;
                } else {
                    thickness_m[k] = hk + vol / areas[k];
                    rep.moved_m3 += vol;
                    vol = 0.0;
                }
                rep.filled_cells += 1;
                let d = great_circle_rad(cents[c], cents[k])? * radius;
                rep.max_distance_m = rep.max_distance_m.max(d);
            }
            if vol <= 0.0 {
                break;
            }
            ring = next.into_iter().filter(|&k| !(thickness_m[k] < cap)).collect();
        }
        if vol > 0.0 {
            rep.delaminated_m3 += vol;
        }
    }
    // (3) Foothill taper: Gauss-Seidel sweeps over ascending pairs (c, k > c)
    // of edge-neighbouring same-plate continental cells; where the thicker
    // one is above the floor and the step exceeds gradient x distance, move
    // volume down-step until the step equals the limit (never taking the
    // thicker cell below the floor).
    let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
    for c in 0..n {
        if !cont(c) {
            continue;
        }
        for &k in &nbrs[c] {
            if k > c && cont(k) && plate[k] == plate[c] {
                let d = great_circle_rad(cents[c], cents[k])? * radius;
                pairs.push((c, k, p.foothill_gradient * d));
            }
        }
    }
    let mut touched: BTreeSet<(usize, usize)> = BTreeSet::new();
    for sweep in 0..TAPER_SWEEPS {
        let mut moved_any = false;
        for &(c, k, lim) in &pairs {
            let (hi, lo) = if thickness_m[c] >= thickness_m[k] { (c, k) } else { (k, c) };
            let (hh, hl) = (thickness_m[hi], thickness_m[lo]);
            if !(hh > p.taper_floor_m) || !(hh - hl > lim + 1.0e-6) {
                continue;
            }
            let (ah, al) = (areas[hi], areas[lo]);
            let mut x = (hh - hl - lim) / (1.0 / ah + 1.0 / al);
            x = x.min((hh - p.taper_floor_m) * ah);
            if !(x > 0.0) {
                continue;
            }
            thickness_m[hi] = hh - x / ah;
            thickness_m[lo] = hl + x / al;
            rep.taper_moved_m3 += x;
            touched.insert((c, k));
            moved_any = true;
        }
        rep.taper_sweeps = sweep as u64 + 1;
        if !moved_any {
            break;
        }
    }
    rep.taper_pairs = touched.len() as u64;
    for &(c, k, lim) in &pairs {
        let (hh, hl) = if thickness_m[c] >= thickness_m[k] { (thickness_m[c], thickness_m[k]) } else { (thickness_m[k], thickness_m[c]) };
        if hh > p.taper_floor_m {
            rep.taper_residual_m = rep.taper_residual_m.max(hh - hl - lim);
        }
    }
    for c in 0..n {
        if cont(c) && !(thickness_m[c] < cap) {
            rep.cells_at_cap += 1;
        }
    }
    Ok(rep)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_state(grid: &GridSpec, h: f64) -> CrustState {
        let n = grid.cells();
        CrustState {
            t_my: 0.0,
            plate: vec![PlateId(0); n],
            kind: vec![CrustKindCodeV0::Continental; n],
            birth_t_my: vec![0.0; n],
            thickness_m: vec![h; n],
            sediment_m: vec![0.0; n],
            eroded_total_m: vec![0.0; n],
            arc_added_m: vec![0.0; n],
            shortening_added_m: vec![0.0; n],
            boundary_class: vec![0; n],
            normal_rate_m_per_my: vec![0.0; n],
            material_pos: centroids(grid).unwrap(),
        }
    }

    #[test]
    fn spill_conserves_volume_and_widens_into_the_own_plate() {
        let grid = GridSpec::new(36, 72, 6.371e6).unwrap();
        let areas = cell_areas(&grid).unwrap();
        let mut st = flat_state(&grid, 35_000.0);
        // Plate 1 west of column 30 (the hinterland of the spike), plate 2 east.
        for c in 0..grid.cells() {
            if c % 72 >= 30 {
                st.plate[c] = PlateId(2);
            } else {
                st.plate[c] = PlateId(1);
            }
        }
        let spike = 18 * 72 + 29; // plate 1, at the contact
        st.thickness_m[spike] = 65_000.0 + 3.0 * 30_000.0; // 3 cells' worth of room
        let v0: f64 = st.thickness_m.iter().zip(&areas).map(|(h, a)| h * a).sum();
        let mut th = vec![0.0; grid.cells()];
        let p = PlateauParams { cap_m: 65_000.0, reach_m: 1.0e6 };
        let r = plateau_spill(&grid, &mut st, &mut th, &p).unwrap();
        let v1: f64 = st.thickness_m.iter().zip(&areas).map(|(h, a)| h * a).sum();
        assert!(((v1 + r.delaminated_m3) - v0).abs() < 1e-9 * v0);
        assert_eq!(r.delaminated_m3, 0.0);
        assert!(r.filled_cells >= 3);
        // Nothing crossed into plate 2, nothing above the cap.
        for c in 0..grid.cells() {
            if st.plate[c] == PlateId(2) {
                assert_eq!(st.thickness_m[c], 35_000.0);
            }
            assert!(st.thickness_m[c] <= 65_000.0 + 1e-6);
        }
        // Lineage and thickening move with the crust.
        let s: f64 = st.shortening_added_m.iter().zip(&areas).map(|(x, a)| x * a).sum();
        let t: f64 = th.iter().zip(&areas).map(|(x, a)| x * a).sum();
        assert!(s.abs() < 1e-6 * v0 && t.abs() < 1e-6 * v0);
    }

    #[test]
    fn spill_beyond_reach_is_delaminated_and_oceanic_excess_removed() {
        let grid = GridSpec::new(36, 72, 6.371e6).unwrap();
        let mut st = flat_state(&grid, 65_000.0); // everything already at the cap
        let c = 18 * 72 + 10;
        st.thickness_m[c] = 70_000.0;
        let o = 10 * 72 + 50;
        st.kind[o] = CrustKindCodeV0::Oceanic;
        st.thickness_m[o] = 80_000.0;
        let mut th = vec![0.0; grid.cells()];
        let p = PlateauParams { cap_m: 65_000.0, reach_m: 600_000.0 };
        let r = plateau_spill(&grid, &mut st, &mut th, &p).unwrap();
        let areas = cell_areas(&grid).unwrap();
        assert!((r.delaminated_m3 - 5_000.0 * areas[c]).abs() < 1e-6 * r.delaminated_m3);
        assert!((r.oceanic_delaminated_m3 - 15_000.0 * areas[o]).abs() < 1e-6 * r.oceanic_delaminated_m3);
        assert_eq!(st.thickness_m[c], 65_000.0);
        assert_eq!(st.thickness_m[o], 65_000.0);
    }

    /// SEED-WORLD-COMPOSE (lane CP): v2 accretes sediment on over-cap
    /// columns, caps, spills into its own plate and tapers the margin;
    /// volume closes with the reported delamination.
    #[test]
    fn plateau_in_belt_accretes_spills_and_tapers() {
        let grid = GridSpec::new(180, 360, 6.371e6).unwrap();
        let areas = cell_areas(&grid).unwrap();
        let cents = centroids(&grid).unwrap();
        let n = grid.cells();
        let kind = vec![CrustKindCodeV0::Continental; n];
        let plate: Vec<PlateId> = (0..n).map(|c| if c % 360 >= 150 { PlateId(2) } else { PlateId(1) }).collect();
        let mut h = vec![35_000.0; n];
        let mut sed = vec![0.0; n];
        let spike = 90 * 360 + 149;
        h[spike] = 64_000.0 + 4.0 * 30_000.0;
        sed[spike] = 3_000.0;
        let v0: f64 = h.iter().zip(&sed).zip(&areas).map(|((a, b), w)| (a + b) * w).sum();
        let p = PlateauInBeltParams { cap_m: 65_000.0, reach_m: 1.0e6, foothill_gradient: 0.1, taper_floor_m: 40_000.0 };
        let r = plateau_in_belt(&grid, &cents, &areas, &kind, &plate, &mut h, &mut sed, &p).unwrap();
        let v1: f64 = h.iter().zip(&sed).zip(&areas).map(|((a, b), w)| (a + b) * w).sum();
        assert!(((v1 + r.delaminated_m3 + r.oceanic_delaminated_m3) - v0).abs() < 1e-9 * v0);
        assert_eq!(sed[spike], 0.0);
        assert!(r.sediment_accreted_m3 > 0.0);
        for c in 0..n {
            assert!(h[c] <= 65_000.0 + 1e-6);
            if plate[c] == PlateId(2) {
                assert_eq!(h[c], 35_000.0);
            }
        }
        // Taper: no same-plate step above the limit where the thick side is above the floor.
        assert!(r.taper_residual_m < 1.0, "residual {}", r.taper_residual_m);
        assert!(r.taper_moved_m3 > 0.0);
        let row = 90 * 360;
        let mx = (0..150).map(|j| h[row + j]).fold(0.0, f64::max);
        assert!(mx < 65_000.0 + 1e-6 && mx > 40_000.0);
    }
}
