//! Routing core (brief S2/S3): priority-flood depression filling, D8-style
//! steepest-descent flow on the REAL neighbor topology, geometric
//! contributing area (m^2) and rainfall-weighted discharge (m^3/s), basin
//! partition with spill levels and endorheic/ocean-connected
//! classification.
//!
//! UNITS (lane E check, 2026-10-11): `rainfall_m_per_my` is effective
//! runoff in m/My (the preview driver converts `uniform_runoff_m_per_yr`
//! explicitly, x 1e6); runoff volume = rate x cell area (m^3/My), summed
//! downstream and divided by 3.15576e13 s/My (Julian year) to give
//! `discharge_proxy_m3_s`. `contributing_area_m2` is NOT rainfall
//! weighted: it is the upstream area, which is what the stream-power law
//! in `erosion.rs` uses (K in 1/My for A in m^2, m = 0.5). Runoff
//! therefore changes discharge, not erosion, at this tier.
//!
//! THE E->F ROUTING API CONTRACT (agent F consumes exactly this; stable):
//! - `route(grid, elevation_m, rainfall_m_per_my, ocean_mask) -> Routing`
//! - `Routing` exposes: `flow_to: Vec<Option<ElementId>>` (steepest-descent
//!   outlet per element; `None` at ocean cells and endorheic lake floors),
//!   `fill_level_m: Vec<f64>` (routing surface ONLY — physical elevation is
//!   NEVER overwritten, S2), `contributing_area_m2: Vec<f64>`,
//!   `discharge_proxy_m3_s: Vec<f64>`, `basins: Vec<BasinRecordV0>` with
//!   `inlet_reaches` filled, `drainage_class(element) -> Option<DrainageClassV0>`.
//! - Determinism: ties in steepest descent break to the LOWEST element id
//!   (declared). All traversals use sorted orders; no HashMap iteration
//!   order leaks into results.
//!
//! GENERIC over the grid: everything goes through `pf-grid` neighbors and
//! areas (antimeridian adjacency comes from the grid, not index math).

use pf_contracts::domain::{BasinId, BasinRecordV0, DrainageClassV0};
use pf_contracts::spatial::ElementId;
use pf_grid::GridSpec;
use std::collections::BTreeMap;

/// Total-order wrapper for finite f64 keys in BTreeMap (levels are
/// validated finite; NaN would be a bug we want to hit loudly).
#[derive(Debug, Clone, Copy, PartialEq)]
struct Fin(f64);

impl Eq for Fin {}

impl PartialOrd for Fin {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Fin {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0
            .partial_cmp(&other.0)
            .expect("Fin holds finite levels only")
    }
}

/// Element-based neighbor access: (target element, distance m) pairs.
fn neighbors_of(grid: &GridSpec, element: u64) -> Result<Vec<(u64, f64)>, RoutingError> {
    let (i, j) = grid
        .ij(element)
        .map_err(|e| RoutingError::Grid(e.to_string()))?;
    grid.neighbors(i, j)
        .map_err(|e| RoutingError::Grid(e.to_string()))?
        .into_iter()
        .map(|(ti, tj, nb)| Ok((grid.flat(ti, tj), nb.distance_m)))
        .collect()
}

fn distance_between(grid: &GridSpec, a: usize, b: usize) -> Result<f64, RoutingError> {
    let (ai, aj) = grid
        .ij(a as u64)
        .map_err(|e| RoutingError::Grid(e.to_string()))?;
    let (bi, bj) = grid
        .ij(b as u64)
        .map_err(|e| RoutingError::Grid(e.to_string()))?;
    grid.center_distance_m(ai, aj, bi, bj)
        .map_err(|e| RoutingError::Grid(e.to_string()))
}

#[derive(Debug, Clone, PartialEq)]
pub enum RoutingError {
    NonFinite {
        what: &'static str,
        value: f64,
        index: usize,
    },
    SizeMismatch {
        what: &'static str,
        expected: usize,
        found: usize,
    },
    NoOutlet {
        why: String,
    },
    Grid(String),
}

impl std::fmt::Display for RoutingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RoutingError::NonFinite { what, value, index } => {
                write!(f, "nonfinite {what} {value:?} at index {index}")
            }
            RoutingError::SizeMismatch {
                what,
                expected,
                found,
            } => {
                write!(f, "{what}: expected {expected}, found {found}")
            }
            RoutingError::NoOutlet { why } => write!(f, "no outlet: {why}"),
            RoutingError::Grid(why) => write!(f, "grid: {why}"),
        }
    }
}

impl std::error::Error for RoutingError {}

/// The routing result: everything F (erosion/sediment) and the hydrology
/// assembly consume. `fill_level_m >= elevation_m` everywhere; equal where
/// unfilled (the field is a separate routing surface, S2).
#[derive(Debug, Clone, PartialEq)]
pub struct Routing {
    /// Every receiver precedes its donor. Complete deterministic flood-pop order.
    pub receiver_first: Vec<usize>,
    pub water_body_id: Vec<u32>,
    pub water_bodies: Vec<pf_contracts::domain::WaterBodyV0>,
    pub flow_to: Vec<Option<ElementId>>,
    pub fill_level_m: Vec<f64>,
    pub contributing_area_m2: Vec<f64>,
    pub discharge_proxy_m3_s: Vec<f64>,
    pub basins: Vec<BasinRecordV0>,
    /// The ocean (base level) mask the routing was computed with, so
    /// consumers can tell an ocean terminal from an endorheic floor.
    pub ocean: Vec<bool>,
    basin_of: BTreeMap<u64, usize>,
}

impl Routing {
    pub fn drainage_class(&self, element: ElementId) -> Option<DrainageClassV0> {
        self.basin_of
            .get(&element.0)
            .map(|&i| self.basins[i].classification)
    }

    /// Index into `basins` for the basin this element drains to (None for
    /// ocean cells that receive no land flow).
    pub fn basin_index(&self, element: ElementId) -> Option<usize> {
        self.basin_of.get(&element.0).copied()
    }
}

/// Route one terrain (S2/S3). Inputs: physical elevation (never modified),
/// rainfall rate per element (m/My), ocean mask (true = ocean/base level).
/// Priority-flood fills depressions to their spill levels; flow follows
/// steepest descent on the FILLED surface (so pits route to their outlets);
/// contributing area is geometric (m^2) and discharge rainfall-weighted
/// (m^3/s); disconnected
/// fill clusters become endorheic basins with outlet + spill level.
pub fn route(
    grid: &GridSpec,
    elevation_m: &[f64],
    rainfall_m_per_my: &[f64],
    ocean_mask: &[bool],
    sea_level_m: f64,
) -> Result<Routing, RoutingError> {
    let n = grid.cells();
    if elevation_m.len() != n || rainfall_m_per_my.len() != n || ocean_mask.len() != n {
        return Err(RoutingError::SizeMismatch {
            what: "routing inputs",
            expected: n,
            found: elevation_m.len(),
        });
    }
    for (i, (&e, &r)) in elevation_m.iter().zip(rainfall_m_per_my).enumerate() {
        if !e.is_finite() {
            return Err(RoutingError::NonFinite {
                what: "elevation_m",
                value: e,
                index: i,
            });
        }
        if !r.is_finite() || r < 0.0 {
            return Err(RoutingError::NonFinite {
                what: "rainfall",
                value: r,
                index: i,
            });
        }
    }
    if !ocean_mask.iter().any(|&o| o) {
        return Err(RoutingError::NoOutlet {
            why: "ocean_mask is empty: no base level".into(),
        });
    }

    if !sea_level_m.is_finite() || (0..n).any(|i| ocean_mask[i] && elevation_m[i] > sea_level_m) {
        return Err(RoutingError::Grid(
            "ocean bed above prescribed water level or nonfinite sea level".into(),
        ));
    }
    // --- Priority-flood (Barnes et al. style): pop lowest, raise seeds. ---
    // The filled surface starts at elevation; oceans pin their cells.
    let mut fill = elevation_m.to_vec();
    for i in 0..n {
        if ocean_mask[i] {
            fill[i] = sea_level_m;
        }
    }
    let mut receiver_first = Vec::with_capacity(n);
    // Ordered queue keyed by (level, element) — deterministic, lowest first.
    let mut queue: BTreeMap<(Fin, u64), ()> = BTreeMap::new();
    let mut visited = vec![false; n];
    // Flood parent: the cell each cell was first reached from. Its fill
    // level is <= the child's, and following parents reaches an ocean
    // seed, so it is a valid receiver ACROSS FLATS (filled lakes) where
    // steepest descent has no strictly lower neighbour (Barnes et al.
    // 2014 priority-flood flow directions).
    let mut parent: Vec<Option<u64>> = vec![None; n];
    for (i, &o) in ocean_mask.iter().enumerate() {
        if o {
            queue.insert((Fin(fill[i]), i as u64), ());
            visited[i] = true;
        }
    }
    while let Some((&(Fin(level), element), ())) = queue.iter().next() {
        queue.remove(&(Fin(level), element));
        receiver_first.push(element as usize);
        for (j, _dist) in neighbors_of(grid, element)? {
            let j = j as usize;
            if visited[j] {
                continue;
            }
            visited[j] = true;
            parent[j] = Some(element);
            // Water can only sit at or above its outlet's level.
            if fill[j] < level {
                fill[j] = level;
            }
            queue.insert((Fin(fill[j]), j as u64), ());
        }
    }
    if receiver_first.len() != n {
        return Err(RoutingError::Grid(
            "priority flood did not cover the grid".into(),
        ));
    }

    // --- D8 steepest descent on the filled surface (lowest-id tie break).
    let mut flow_to: Vec<Option<ElementId>> = vec![None; n];
    for i in 0..n as u64 {
        if ocean_mask[i as usize] {
            continue; // ocean cells are terminals
        }
        let mut best: Option<(f64, u64)> = None; // (slope, target)
        for (j, dist) in neighbors_of(grid, i)? {
            // Steepest descent over DISTANCE (the metric, not index delta).
            let slope = (fill[i as usize] - fill[j as usize]) / dist.max(1.0);
            if slope <= 0.0 {
                continue;
            }
            if let Some((bs, _)) = best {
                if slope <= bs {
                    continue;
                }
            }
            best = Some((slope, j));
        }
        flow_to[i as usize] = best.map(|(_, j)| ElementId(j));
        if flow_to[i as usize].is_none() {
            // No strictly lower neighbour: a flat on the filled surface
            // (lake floor or filled depression). Drain toward the spill
            // through the flood parent; only an unvisited cell (broken
            // topology) stays a terminal and becomes an endorheic floor.
            flow_to[i as usize] = parent[i as usize].map(ElementId);
        }
    }

    // Pop order is receiver-first for strict-lower and flood-parent edges.
    let mut rank = vec![usize::MAX; n];
    for (k, &i) in receiver_first.iter().enumerate() {
        rank[i] = k;
    }
    for i in 0..n {
        if let Some(ElementId(j)) = flow_to[i] {
            if rank[j as usize] >= rank[i] {
                return Err(RoutingError::Grid(format!(
                    "receiver order violated: {i} -> {j}"
                )));
            }
        }
    }
    let mut area = vec![0.0; n];
    let mut runoff = vec![0.0; n];
    for i in 0..n {
        area[i] = grid
            .cell_area_m2(i / grid.nlon())
            .map_err(|e| RoutingError::Grid(e.to_string()))?;
        runoff[i] = rainfall_m_per_my[i] * area[i];
    }
    for &i in receiver_first.iter().rev() {
        if let Some(ElementId(j)) = flow_to[i] {
            area[j as usize] += area[i];
            runoff[j as usize] += runoff[i];
        }
    }
    const SECONDS_PER_MY: f64 = 31_557_600_000_000.0;
    let discharge: Vec<f64> = runoff.iter().map(|v| v / SECONDS_PER_MY).collect();
    let (water_body_id, water_bodies) =
        crate::water::bodies(grid, elevation_m, &fill, ocean_mask, &flow_to)?;

    // --- Basin partition: terminal cells group their uphill networks.
    // Terminals: ocean cells (ocean-connected networks) and None-flow
    // cells (endorheic lake floors). Build reverse adjacency, flood each
    // terminal's upstream set.
    let mut upstream: Vec<Vec<u64>> = vec![Vec::new(); n];
    for (i, t) in flow_to.iter().enumerate() {
        if let Some(ElementId(j)) = t {
            upstream[*j as usize].push(i as u64);
        }
    }
    let mut basins: Vec<BasinRecordV0> = Vec::new();
    let mut basin_of: BTreeMap<u64, usize> = BTreeMap::new();
    // Ocean-connected: group by terminal ocean cell, then classify.
    let mut next_basin = 0u64;
    let mut terminals: Vec<(u64, bool)> = Vec::new(); // (element, is_ocean)
    for i in 0..n as u64 {
        if ocean_mask[i as usize] && !upstream[i as usize].is_empty() {
            terminals.push((i, true));
        } else if flow_to[i as usize].is_none() && !ocean_mask[i as usize] {
            terminals.push((i, false));
        }
    }
    terminals.sort_by_key(|&(e, o)| (o, e)); // deterministic order
    for &(terminal, is_ocean) in &terminals {
        let basin = BasinId(next_basin);
        // Collect the upstream network (iterative DFS, sorted pushes).
        let mut members = vec![terminal];
        let mut stack = vec![terminal];
        while let Some(e) = stack.pop() {
            for &u in &upstream[e as usize] {
                if !basin_of.contains_key(&u) {
                    basin_of.insert(u, basins.len());
                    members.push(u);
                    stack.push(u);
                }
            }
        }
        // The terminal belongs to its own basin (ocean terminals included,
        // so a reach ending at the coast resolves to the ocean basin).
        basin_of.insert(terminal, basins.len());
        members.sort_unstable();
        let spill = fill[terminal as usize];
        let outlet = if is_ocean {
            None
        } else {
            Some(ElementId(terminal))
        };
        basins.push(BasinRecordV0 {
            basin_id: basin,
            classification: if is_ocean {
                DrainageClassV0::OceanConnected
            } else {
                DrainageClassV0::Endorheic
            },
            spill_level_m: if is_ocean { None } else { Some(spill) },
            outlet_element: outlet,
            inlet_reaches: Vec::new(), // hydrology.rs fills these
        });
        next_basin += 1;
        let _ = members;
    }

    Ok(Routing {
        receiver_first,
        water_body_id,
        water_bodies,
        flow_to,
        fill_level_m: fill,
        contributing_area_m2: area,
        discharge_proxy_m3_s: discharge,
        basins,
        ocean: ocean_mask.to_vec(),
        basin_of,
    })
}

/// CLIMATE-GENESIS (lane CG, 2026-10-20): runoff-weighted contributing
/// area A_eff,i = sum over the catchment of i of a_j R_j / r_ref (m^2), the
/// area form of stream power on discharge (E = K A_eff^m S^n with K quoted
/// at the reference runoff r_ref; approximation `runoff-weighted-area-v1`).
/// Accumulated over `routing.receiver_first` in reverse with the cell areas
/// `route` uses, so with every R equal to r_ref each factor is exactly 1.0
/// and the result equals `routing.contributing_area_m2` bit for bit. The
/// flow graph is the routing's; runoff never moves flow directions.
/// `runoff` and `r_ref` share one unit (any; only the ratio enters).
pub fn runoff_weighted_area_m2(
    grid: &GridSpec,
    routing: &Routing,
    runoff: &[f64],
    r_ref: f64,
) -> Result<Vec<f64>, RoutingError> {
    let n = grid.cells();
    if runoff.len() != n || routing.flow_to.len() != n || routing.receiver_first.len() != n {
        return Err(RoutingError::SizeMismatch {
            what: "runoff-weighted area inputs",
            expected: n,
            found: runoff.len(),
        });
    }
    if !(r_ref.is_finite() && r_ref > 0.0) {
        return Err(RoutingError::NonFinite { what: "reference runoff", value: r_ref, index: 0 });
    }
    if let Some(i) = runoff.iter().position(|r| !(r.is_finite() && *r >= 0.0)) {
        return Err(RoutingError::NonFinite { what: "runoff", value: runoff[i], index: i });
    }
    let mut area = vec![0.0; n];
    for i in 0..n {
        area[i] = grid
            .cell_area_m2(i / grid.nlon())
            .map_err(|e| RoutingError::Grid(e.to_string()))?
            * (runoff[i] / r_ref);
    }
    for &i in routing.receiver_first.iter().rev() {
        if let Some(ElementId(j)) = routing.flow_to[i] {
            area[j as usize] += area[i];
        }
    }
    Ok(area)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tiny deterministic terrain helper: value from a formula.
    fn terrain(n: usize, f: impl Fn(usize) -> f64) -> Vec<f64> {
        (0..n).map(f).collect()
    }

    #[test]
    fn ocean_only_world_routes_everything_downhill() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        // One ocean at element 0, a monotone ramp away from it.
        let elev = terrain(n, |i| if i == 0 { f64::MIN } else { 100.0 + i as f64 });
        let mut ocean = vec![false; n];
        ocean[0] = true;
        let rain = vec![1.0; n];
        let r = route(&grid, &elev, &rain, &ocean, 0.0).unwrap();
        // Every non-ocean cell flows somewhere (monotone ramp has outlets).
        for (i, t) in r.flow_to.iter().enumerate() {
            if i != 0 {
                assert!(t.is_some(), "element {i} orphaned on a monotone ramp");
            }
        }
        // Contributing area strictly increases along flow paths.
        let mut e = 1usize;
        while let Some(ElementId(j)) = r.flow_to[e] {
            assert!(
                r.contributing_area_m2[j as usize] > r.contributing_area_m2[e],
                "area not monotone at {e}"
            );
            e = j as usize;
            if e == 0 {
                break;
            }
        }
        // Ocean cell has the largest area (everything drains there-ish).
        assert!(
            r.contributing_area_m2[0]
                >= r.contributing_area_m2
                    .iter()
                    .cloned()
                    .fold(f64::MIN, f64::max)
                    * 0.999
        );
        // No endorheic basins in this world.
        assert!(r
            .basins
            .iter()
            .all(|b| b.classification == DrainageClassV0::OceanConnected));
    }

    #[test]
    fn nested_pit_fills_to_spill_and_drains_across_the_flat_to_the_ocean() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        // Flat world at 50 m with a deep nested pit around element 20-21,
        // one ocean at element 100 (within 8x16 = 128 cells).
        let mut elev = vec![50.0; n];
        for &i in &[20, 21, 36, 37] {
            elev[i] = -200.0;
        }
        elev[22] = -100.0; // shallower ring cell
        let mut ocean = vec![false; n];
        ocean[100] = true;
        let r = route(&grid, &elev, &vec![1.0; n], &ocean, 50.0).unwrap();
        // The pit cells are filled to the rim level (a lake), physical
        // elevation untouched (fill is a separate surface, S2).
        for &i in &[20, 21, 36, 37, 22] {
            assert!(
                r.fill_level_m[i] >= 50.0 - 1.0e-9,
                "pit {i} not filled: {}",
                r.fill_level_m[i]
            );
        }
        assert_eq!(r.fill_level_m[100], elev[100]);
        // Every non-ocean cell has a receiver (flats drain via the flood
        // parent) and every path ends at the ocean: one basin, ocean
        // connected, no endorheic floors on a connected surface.
        for c in 0..n {
            if ocean[c] {
                assert!(r.flow_to[c].is_none());
                continue;
            }
            let mut cur = c;
            let mut hops = 0;
            while let Some(ElementId(j)) = r.flow_to[cur] {
                cur = j as usize;
                hops += 1;
                assert!(hops <= n, "cycle from {c}");
            }
            assert_eq!(cur, 100, "cell {c} drains to {cur}, not the ocean");
        }
        assert_eq!(r.basins.len(), 1);
        assert_eq!(r.basins[0].classification, DrainageClassV0::OceanConnected);
        assert_eq!(
            r.basin_index(ElementId(100)),
            Some(0),
            "ocean terminal belongs to its basin"
        );
        assert_eq!(
            r.basin_index(ElementId(20)),
            Some(0),
            "pit member drains to the ocean basin"
        );
        assert_eq!(r.ocean, ocean);
    }

    #[test]
    fn rainfall_weighting_changes_discharge_not_geometry() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        let elev = terrain(n, |i| (i as f64).sin() * 500.0);
        let mut ocean = vec![false; n];
        ocean[5] = true;
        let uniform = route(&grid, &elev, &vec![1.0; n], &ocean, 0.0).unwrap();
        let mut rainy = vec![1.0; n];
        rainy[7] = 10.0;
        let wet = route(&grid, &elev, &rainy, &ocean, 0.0).unwrap();
        // Same flow graph (rainfall does not move flow directions).
        assert_eq!(uniform.flow_to, wet.flow_to);
        // Discharge at the rain cell scales ~ x10 (its own cell's rain).
        assert!(wet.discharge_proxy_m3_s[7] > uniform.discharge_proxy_m3_s[7] * 5.0);
        // Contributing AREA is geometric: unchanged.
        assert!(uniform
            .contributing_area_m2
            .iter()
            .zip(&wet.contributing_area_m2)
            .all(|(a, b)| (a - b).abs() < 1e-6));
    }

    #[test]
    fn runoff_weighted_area_is_bit_exact_at_reference_and_weights_catchments() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        let mut elev = terrain(n, |i| (i as f64).sin() * 500.0 + 0.01 * i as f64);
        let mut ocean = vec![false; n];
        for o in [5, 77] {
            ocean[o] = true;
            elev[o] = -2000.0;
        }
        let r = route(&grid, &elev, &vec![1.0; n], &ocean, -1500.0).unwrap();
        // Identity: R == r_ref everywhere gives the geometric area bit for bit.
        let same = runoff_weighted_area_m2(&grid, &r, &vec![0.7; n], 0.7).unwrap();
        assert!(same.iter().zip(&r.contributing_area_m2).all(|(a, b)| a.to_bits() == b.to_bits()));
        // Known answer: A_eff = sum over the catchment of a R / r_ref.
        let runoff: Vec<f64> = (0..n).map(|i| 0.25 + (i % 5) as f64 * 0.5).collect();
        let eff = runoff_weighted_area_m2(&grid, &r, &runoff, 1.0).unwrap();
        for c in 0..n {
            let mut want = 0.0;
            for d in 0..n {
                // d drains through c?
                let mut cur = d;
                let mut hit = cur == c;
                while !hit {
                    match r.flow_to[cur] {
                        Some(ElementId(j)) => {
                            cur = j as usize;
                            hit = cur == c;
                        }
                        None => break,
                    }
                }
                if hit {
                    want += grid.cell_area_m2(d / grid.nlon()).unwrap() * runoff[d];
                }
            }
            assert!((eff[c] - want).abs() <= 1e-9 * want.max(1.0), "cell {c}: {} vs {want}", eff[c]);
        }
        // Zero runoff: zero area, never NaN; bad inputs loud.
        let zero = runoff_weighted_area_m2(&grid, &r, &vec![0.0; n], 1.0).unwrap();
        assert!(zero.iter().all(|a| *a == 0.0));
        assert!(runoff_weighted_area_m2(&grid, &r, &vec![1.0; n], 0.0).is_err());
        let mut bad = vec![1.0; n];
        bad[3] = -1.0;
        assert!(runoff_weighted_area_m2(&grid, &r, &bad, 1.0).is_err());
    }

    #[test]
    fn routing_is_deterministic_and_validates_inputs() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        let elev = terrain(n, |i| (i as f64) * 3.0);
        let mut ocean = vec![false; n];
        ocean[0] = true;
        let a = route(&grid, &elev, &vec![2.0; n], &ocean, 0.0).unwrap();
        let b = route(&grid, &elev, &vec![2.0; n], &ocean, 0.0).unwrap();
        assert_eq!(a, b);
        // Bad inputs loud.
        let mut nan = elev.clone();
        nan[3] = f64::NAN;
        assert!(route(&grid, &nan, &vec![1.0; n], &ocean, 0.0).is_err());
        let mut neg = vec![1.0; n];
        neg[3] = -0.5;
        assert!(route(&grid, &elev, &neg, &ocean, 0.0).is_err());
        let no_ocean = vec![false; n];
        let err = route(&grid, &elev, &vec![1.0; n], &no_ocean, 0.0).unwrap_err();
        assert!(err.to_string().contains("no base level"), "{err}");
        assert!(route(&grid, &elev[..3], &vec![1.0; n], &ocean, 0.0).is_err());
    }
}
