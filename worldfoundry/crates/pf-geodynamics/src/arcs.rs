//! SEED-MAGMATISM (lane MG, 2026-10-13): the arc part of
//! `slab-arc-backarc-llsvp-v1` (pre-registered in
//! `unpublished-design-record`).
//!
//! Physics as modelled: the subducting slab is denser than the mantle and
//! sinks; near `slab_depth_m` (~100 km) it releases water that lowers the
//! melting point of the mantle wedge above it (flux melting). The melt is
//! lighter and rises, so the volcanic arc sits on the OVERRIDING plate at the
//! map distance d = H / tan(dip) from the trench. This module:
//!
//! 1. finds the subduction fronts: every downgoing trench cell (oceanic,
//!    Subduction class, normal rate > 0, with an edge neighbour of another
//!    plate that overrides it: continental, or younger oceanic, equal age
//!    the lower plate id; the trench.rs and kinematic polarity rule);
//! 2. places the arc point d along the great circle from the trench axis into
//!    the overriding plate (dip from the slab age) and adds arc crust there
//!    (rate p * clamp(v) * L per step; a named magmatic source);
//! 3. turns oceanic columns built mostly of arc crust and at least
//!    `continental_thickness_m` thick into continental crust (island arcs
//!    mature into new continental crust and later weld on by collision);
//!    Revision 1 (PREREG): arc magma builds a column only up to `crust_cap_m`.
//! 4. keeps a per-cell arc heat marker (hot arc lithosphere) whose
//!    elevation term is applied by the caller.
//!
//! Determinism: fronts in ascending cell order, neighbours in grid order,
//! discs from `orogen::disc_weights` (ascending), no maps with hashed order.

use crate::kinematic::CrustState;
use crate::rotation::great_circle_rad;
use pf_contracts::domain::PlateId;
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_grid::GridSpec;

/// Arc constants of `slab-arc-backarc-llsvp-v1` (laws.rs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArcParams {
    pub slab_depth_m: f64,
    pub dip_young_deg: f64,
    pub dip_old_deg: f64,
    pub dip_age_my: f64,
    pub productivity_m: f64,
    pub rate_floor_m_per_my: f64,
    pub rate_cap_m_per_my: f64,
    pub half_width_m: f64,
    pub continental_thickness_m: f64,
    pub thermal_uplift_m: f64,
    pub thermal_sigma_m: f64,
    pub thermal_decay_my: f64,
    /// Revision 1: arc magma only builds a column up to this thickness (m).
    pub crust_cap_m: f64,
}

impl ArcParams {
    pub fn validate(&self) -> Result<(), String> {
        let pos = [
            self.slab_depth_m,
            self.dip_young_deg,
            self.dip_old_deg,
            self.dip_age_my,
            self.rate_cap_m_per_my,
            self.continental_thickness_m,
            self.thermal_sigma_m,
            self.thermal_decay_my,
            self.crust_cap_m,
        ];
        let nonneg = [
            self.productivity_m,
            self.rate_floor_m_per_my,
            self.half_width_m,
            self.thermal_uplift_m,
        ];
        if pos.iter().any(|v| !(v.is_finite() && *v > 0.0))
            || nonneg.iter().any(|v| !(v.is_finite() && *v >= 0.0))
            || !(self.dip_young_deg < 90.0 && self.dip_old_deg < 90.0)
            || !(self.rate_floor_m_per_my <= self.rate_cap_m_per_my)
            || !(self.crust_cap_m > self.continental_thickness_m)
        {
            return Err(format!("arc parameters out of range: {self:?}"));
        }
        Ok(())
    }

    /// Effective slab dip (deg) for a slab of age `age_my`.
    pub fn dip_deg(&self, age_my: f64) -> f64 {
        let f = (age_my.max(0.0) / self.dip_age_my).min(1.0);
        self.dip_young_deg + (self.dip_old_deg - self.dip_young_deg) * f
    }

    /// Map distance (m) from the trench axis to the arc.
    pub fn arc_distance_m(&self, age_my: f64) -> f64 {
        self.slab_depth_m / self.dip_deg(age_my).to_radians().tan()
    }

    /// Arc crust flux per unit trench length (m^2/My) at convergence v.
    pub fn flux_m2_per_my(&self, rate_m_per_my: f64) -> f64 {
        self.productivity_m
            * rate_m_per_my.clamp(self.rate_floor_m_per_my, self.rate_cap_m_per_my)
    }
}

/// One subduction front sample (one downgoing trench cell).
#[derive(Debug, Clone, PartialEq)]
pub struct Front {
    /// The downgoing trench cell and its plate.
    pub cell: u32,
    pub down_plate: PlateId,
    /// The overriding plate (first overriding edge neighbour).
    pub over_plate: PlateId,
    /// Trench axis point (unit vector) and the unit tangent there pointing
    /// into the overriding plate.
    pub axis: [f64; 3],
    pub dir: [f64; 3],
    pub age_my: f64,
    pub rate_m_per_my: f64,
    /// Trench length this sample stands for (m): sqrt(cell area).
    pub length_m: f64,
    /// Planned and realised arc distance (m) and the arc cell / point.
    pub arc_distance_planned_m: f64,
    pub arc_distance_m: f64,
    pub arc_cell: u32,
    pub arc_point: [f64; 3],
    /// Arc cell continental at detection (continental arc) or not (island arc).
    pub continental_arc: bool,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = dot(a, a).sqrt();
    if n > 1.0e-15 && n.is_finite() {
        Some([a[0] / n, a[1] / n, a[2] / n])
    } else {
        None
    }
}

/// Point at angle `theta` (rad) along the great circle from `p` in the
/// tangent direction `t` (unit, perpendicular to p).
pub fn along(p: [f64; 3], t: [f64; 3], theta: f64) -> [f64; 3] {
    let (s, c) = theta.sin_cos();
    [c * p[0] + s * t[0], c * p[1] + s * t[1], c * p[2] + s * t[2]]
}

/// The nearest cell of a unit vector (same rule as `kinematic::nearest_cell`).
pub fn cell_of(grid: &GridSpec, p: [f64; 3]) -> Result<usize, String> {
    crate::kinematic::nearest_cell(grid, p).map_err(|e| e.to_string())
}

/// True when `nb` overrides `c` at a subduction contact (trench.rs rule).
pub fn overrides(state: &CrustState, c: usize, nb: usize) -> bool {
    if state.kind[c] != CrustKindCodeV0::Oceanic {
        return false;
    }
    match state.kind[nb] {
        CrustKindCodeV0::Continental => true,
        CrustKindCodeV0::Oceanic => {
            state.birth_t_my[nb] > state.birth_t_my[c]
                || (state.birth_t_my[nb] == state.birth_t_my[c] && state.plate[nb] < state.plate[c])
        }
    }
}

/// Edge (4-) neighbours of cell `c` in grid neighbour order.
pub fn edge_neighbours(grid: &GridSpec, c: usize) -> Result<Vec<usize>, String> {
    let (i, j) = grid.ij(c as u64).map_err(|e| e.to_string())?;
    Ok(grid
        .neighbors(i, j)
        .map_err(|e| e.to_string())?
        .into_iter()
        .filter(|(_, _, nb)| (nb.di == 0) != (nb.dj == 0))
        .map(|(ti, tj, _)| ti * grid.nlon() + tj)
        .collect())
}

/// Find the fronts of `state` and their arc cells.
pub fn find_fronts(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    state: &CrustState,
    p: &ArcParams,
) -> Result<Vec<Front>, String> {
    p.validate()?;
    let n = grid.cells();
    if state.cells() != n || cents.len() != n || areas.len() != n {
        return Err("find_fronts: field sizes differ from the grid".into());
    }
    let r = grid.radius_m();
    let dx = r * std::f64::consts::PI / grid.nlat() as f64;
    let sub = BoundaryClassCodeV0::Subduction.code();
    let mut out = Vec::new();
    for c in 0..n {
        if state.boundary_class[c] != sub
            || state.kind[c] != CrustKindCodeV0::Oceanic
            || !(state.normal_rate_m_per_my[c] > 0.0)
        {
            continue;
        }
        let nbs: Vec<usize> = edge_neighbours(grid, c)?
            .into_iter()
            .filter(|&nb| state.plate[nb] != state.plate[c] && overrides(state, c, nb))
            .collect();
        if nbs.is_empty() {
            continue;
        }
        let over = state.plate[nbs[0]];
        let same: Vec<usize> = nbs.iter().cloned().filter(|&k| state.plate[k] == over).collect();
        let mut m = [0.0; 3];
        for &k in &same {
            for x in 0..3 {
                m[x] += cents[k][x];
            }
        }
        let Some(m) = norm(m) else { continue };
        let Some(axis) = norm([cents[c][0] + m[0], cents[c][1] + m[1], cents[c][2] + m[2]]) else {
            continue;
        };
        let md = dot(m, axis);
        let Some(dir) = norm([m[0] - md * axis[0], m[1] - md * axis[1], m[2] - md * axis[2]]) else {
            continue;
        };
        let age = state.age_my(c);
        let planned = p.arc_distance_m(age);
        // Farthest half-cell sample along the path (from the planned arc
        // point back toward the trench) that lies on the overriding plate.
        let mut arc_cell = same[0];
        let mut arc_point = cents[same[0]];
        let mut realised = great_circle_rad(axis, cents[same[0]]).map_err(|e| e.to_string())? * r;
        let mut d = planned;
        while d > 0.0 {
            let q = along(axis, dir, d / r);
            let k = cell_of(grid, q)?;
            if state.plate[k] == over {
                arc_cell = k;
                arc_point = q;
                realised = d;
                break;
            }
            d -= 0.5 * dx;
        }
        out.push(Front {
            cell: c as u32,
            down_plate: state.plate[c],
            over_plate: over,
            axis,
            dir,
            age_my: age,
            rate_m_per_my: state.normal_rate_m_per_my[c],
            length_m: areas[c].sqrt(),
            arc_distance_planned_m: planned,
            arc_distance_m: realised,
            arc_cell: arc_cell as u32,
            arc_point,
            continental_arc: state.kind[arc_cell] == CrustKindCodeV0::Continental,
        });
    }
    Ok(out)
}

/// What one arc pass did (volumes m^3, areas m^2).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ArcReport {
    pub fronts: u64,
    pub island_fronts: u64,
    pub continental_fronts: u64,
    pub front_length_m: f64,
    pub arc_m3: f64,
    /// Revision 1: arc volume not created because the column reached the cap.
    pub arc_capped_m3: f64,
    pub arc_on_oceanic_m3: f64,
    pub arc_on_continental_m3: f64,
    /// Fronts whose arc point fell short of the planned distance (narrow
    /// overriding plate) and the median realised distance (m).
    pub shortened_fronts: u64,
    pub median_arc_distance_m: f64,
    pub converted_cells: u64,
    pub converted_area_m2: f64,
    pub converted_m3: f64,
}

/// Add the arc crust of `fronts` over `dt_my` to `state` (thickness and
/// `arc_added_m`); per-cell increments are added into `dh` (m). Returns the
/// report (no conversion yet).
pub fn arc_pass(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    state: &mut CrustState,
    fronts: &[Front],
    p: &ArcParams,
    dt_my: f64,
    dh: &mut [f64],
) -> Result<ArcReport, String> {
    let n = grid.cells();
    if dh.len() != n {
        return Err("arc_pass: dh size differs from the grid".into());
    }
    let mut rep = ArcReport::default();
    let mut add = vec![0.0f64; n];
    let mut dists: Vec<f64> = Vec::with_capacity(fronts.len());
    let observing = crate::volcanoes::supply_capture::enabled();
    let mut observed = Vec::new();
    for f in fronts {
        rep.fronts += 1;
        if f.continental_arc {
            rep.continental_fronts += 1;
        } else {
            rep.island_fronts += 1;
        }
        rep.front_length_m += f.length_m;
        if f.arc_distance_m + 1.0 < f.arc_distance_planned_m {
            rep.shortened_fronts += 1;
        }
        dists.push(f.arc_distance_m);
        let v = p.flux_m2_per_my(f.rate_m_per_my) * f.length_m * dt_my;
        let over = f.over_plate;
        let plate = &state.plate;
        let w = crate::orogen::disc_weights(grid, cents, f.arc_cell as usize, p.half_width_m, &|k| {
            plate[k] == over
        })?;
        crate::orogen::place(&w, areas, v, &mut add);
        if observing {
            observed.push(crate::volcanoes::supply_capture::Sample::placed(
                "arc", format!("front/{}/{}/{}", f.down_plate.0, f.over_plate.0, f.cell),
                f.cell, state.t_my - dt_my, state.t_my, v, &w, areas,
            )?);
        }
    }
    // Observe the simultaneous cap in its original order. Allocation among
    // fronts is proportional to their pre-cap cell contributions, never a
    // second production pass or a priority claim by the first front.
    if observing {
        for s in &mut observed {
            for x in &mut s.cells {
                let c = x.cell as usize;
                let room = (p.crust_cap_m - state.thickness_m[c]).max(0.0);
                x.admitted_m3 = x.gross_m3 * if add[c] > room { room / add[c] } else { 1.0 };
            }
        }
        crate::volcanoes::supply_capture::record(observed)?;
    }
    for c in 0..n {
        if add[c] != 0.0 {
            // Revision 1: only up to the arc crust cap (the rest is not created).
            let room = (p.crust_cap_m - state.thickness_m[c]).max(0.0);
            if add[c] > room {
                rep.arc_capped_m3 += (add[c] - room) * areas[c];
                add[c] = room;
            }
            if add[c] == 0.0 {
                continue;
            }
            state.thickness_m[c] += add[c];
            state.arc_added_m[c] += add[c];
            dh[c] += add[c];
            let v = add[c] * areas[c];
            rep.arc_m3 += v;
            match state.kind[c] {
                CrustKindCodeV0::Oceanic => rep.arc_on_oceanic_m3 += v,
                CrustKindCodeV0::Continental => rep.arc_on_continental_m3 += v,
            }
        }
    }
    if !dists.is_empty() {
        dists.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
        rep.median_arc_distance_m = dists[dists.len() / 2];
    }
    Ok(rep)
}

/// Island arc maturation: oceanic columns at least `continental_thickness_m`
/// thick whose arc crust is at least half the column become continental.
/// Returns the converted cells (ascending) and fills the report.
pub fn convert_arcs(
    state: &mut CrustState,
    areas: &[f64],
    p: &ArcParams,
    rep: &mut ArcReport,
) -> Vec<u32> {
    let mut out = Vec::new();
    for c in 0..state.cells() {
        if state.kind[c] == CrustKindCodeV0::Oceanic
            && state.thickness_m[c] >= p.continental_thickness_m
            && state.arc_added_m[c] >= 0.5 * state.thickness_m[c]
        {
            state.kind[c] = CrustKindCodeV0::Continental;
            rep.converted_cells += 1;
            rep.converted_area_m2 += areas[c];
            rep.converted_m3 += areas[c] * state.thickness_m[c];
            out.push(c as u32);
        }
    }
    out
}

/// Update the arc heat marker: decay by exp(-dt/tau), then raise to the
/// gaussian of every front's arc point on its overriding plate (cut at
/// 3 sigma). `plate` is the state's plate field.
pub fn update_heat(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    plate: &[PlateId],
    fronts: &[Front],
    p: &ArcParams,
    dt_my: f64,
    heat: &mut [f64],
) -> Result<(), String> {
    let n = grid.cells();
    if heat.len() != n || plate.len() != n {
        return Err("update_heat: field sizes differ from the grid".into());
    }
    let decay = (-dt_my / p.thermal_decay_my).exp();
    for h in heat.iter_mut() {
        *h *= decay;
    }
    let reach = 3.0 * p.thermal_sigma_m;
    let r = grid.radius_m();
    // SEED-PERF (lane SP): per front in parallel (pure): its disc and the
    // gaussian of every overriding-plate cell in it, in disc order (up to
    // the first error); then the max fold serially in front order, bit for
    // bit the serial loop (including what it had written before an error).
    let per_front = crate::kinematic::par::map_heavy(fronts.len(), |fi| -> (Vec<(usize, f64)>, Option<String>) {
        let f = &fronts[fi];
        let over = f.over_plate;
        let w = match crate::orogen::disc_weights(grid, cents, f.arc_cell as usize, reach, &|k| plate[k] == over) {
            Ok(w) => w,
            Err(e) => return (Vec::new(), Some(e)),
        };
        let mut gs = Vec::with_capacity(w.len());
        for (k, _) in w {
            if plate[k] != over {
                continue;
            }
            let d = match great_circle_rad(f.arc_point, cents[k]) {
                Ok(d) => d * r,
                Err(e) => return (gs, Some(e.to_string())),
            };
            gs.push((k, (-(d / p.thermal_sigma_m).powi(2)).exp()));
        }
        (gs, None)
    });
    for (gs, err) in per_front {
        for (k, g) in gs {
            if g > heat[k] {
                heat[k] = g;
            }
        }
        if let Some(e) = err {
            return Err(e);
        }
    }
    Ok(())
}

/// The arc thermal elevation term (water-loaded m) of a heat field.
pub fn thermal_term(heat: &[f64], p: &ArcParams) -> Vec<f64> {
    heat.iter().map(|h| p.thermal_uplift_m * h).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> ArcParams {
        ArcParams {
            slab_depth_m: 105_000.0,
            dip_young_deg: 30.0,
            dip_old_deg: 45.0,
            dip_age_my: 150.0,
            productivity_m: 600.0,
            rate_floor_m_per_my: 20_000.0,
            rate_cap_m_per_my: 100_000.0,
            half_width_m: 75_000.0,
            continental_thickness_m: 20_000.0,
            thermal_uplift_m: 1000.0,
            thermal_sigma_m: 150_000.0,
            thermal_decay_my: 15.0,
            crust_cap_m: 50_000.0,
        }
    }

    #[test]
    fn dip_and_distance_follow_slab_age() {
        let p = params();
        assert!((p.dip_deg(0.0) - 30.0).abs() < 1e-12);
        assert!((p.dip_deg(300.0) - 45.0).abs() < 1e-12);
        let young = p.arc_distance_m(0.0);
        let old = p.arc_distance_m(150.0);
        assert!((young - 105_000.0 / 30f64.to_radians().tan()).abs() < 1e-6);
        assert!((old - 105_000.0).abs() < 1e-6);
        assert!(young > old);
        assert_eq!(p.flux_m2_per_my(5_000.0), 600.0 * 20_000.0);
        assert_eq!(p.flux_m2_per_my(500_000.0), 600.0 * 100_000.0);
    }

    #[test]
    fn along_moves_by_the_angle() {
        let p = [1.0, 0.0, 0.0];
        let t = [0.0, 1.0, 0.0];
        let q = along(p, t, 0.1);
        assert!((great_circle_rad(p, q).unwrap() - 0.1).abs() < 1e-12);
    }
}
