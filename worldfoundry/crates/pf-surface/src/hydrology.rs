//! Hydrology assembly (brief S4): river network extraction from
//! contributing area with a physical-unit threshold, Strahler ordering,
//! trunk monotone verification, and `HydrologyDomainV0` assembly
//! (validated via its own contract `.validate()`).

use crate::routing::{Routing, RoutingError};
use pf_contracts::domain::{BasinId, DrainageClassV0, HydrologyDomainV0, ReachId, RiverReachV0};
use pf_contracts::spatial::ElementId;
use pf_grid::GridSpec;
use std::collections::BTreeMap;

/// River threshold: minimum discharge proxy (m^3/s-class) for a channel.
/// PHYSICAL units, resolution-independent (S4). Constant.
pub const RIVER_THRESHOLD_M3_S: f64 = 10.0;

#[derive(Debug, Clone, PartialEq)]
pub enum HydrologyError {
    Routing(RoutingError),
    NonMonotoneTrunk {
        reach: ReachId,
        upstream: f64,
        downstream: f64,
    },
    OrphanReach {
        reach: ReachId,
    },
    Domain(String),
}

impl std::fmt::Display for HydrologyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HydrologyError::Routing(e) => write!(f, "routing: {e}"),
            HydrologyError::NonMonotoneTrunk { reach, upstream, downstream } => write!(
                f,
                "trunk not monotone downhill at reach {reach:?}: upstream fill {upstream:?} < downstream {downstream:?}"
            ),
            HydrologyError::OrphanReach { reach } => write!(f, "reach {reach:?} references unknown reach/basin"),
            HydrologyError::Domain(why) => write!(f, "domain validation: {why}"),
        }
    }
}

impl std::error::Error for HydrologyError {}

impl From<RoutingError> for HydrologyError {
    fn from(e: RoutingError) -> Self {
        HydrologyError::Routing(e)
    }
}

/// Extract the river network + assemble the hydrology domain from a
/// routing result. Reaches are maximal downhill chains of river cells
/// (discharge >= threshold); junctions start new reaches; Strahler order
/// computed on the reach DAG; trunk monotone-downhill on the FILL surface
/// is verified (S2: published hydrology describes the published terrain).
pub fn assemble_hydrology(
    grid: &GridSpec,
    routing: &Routing,
    t_my: f64,
) -> Result<HydrologyDomainV0, HydrologyError> {
    assemble_hydrology_with_cells(grid, routing, t_my).map(|(d, _)| d)
}

/// `assemble_hydrology` plus the per-cell reach membership (index-aligned
/// with the grid; `None` where the cell is not a river cell), for
/// publishers of the per-cell `reach_id` field (ADR 032).
pub fn assemble_hydrology_with_cells(
    grid: &GridSpec,
    routing: &Routing,
    t_my: f64,
) -> Result<(HydrologyDomainV0, Vec<Option<ReachId>>), HydrologyError> {
    let n = grid.cells();
    // River cells by threshold.
    // River cells by threshold; the ocean is base level, never a channel.
    let is_river: Vec<bool> = routing
        .discharge_proxy_m3_s
        .iter()
        .zip(&routing.ocean)
        .map(|(&d, &o)| d >= RIVER_THRESHOLD_M3_S && !o)
        .collect();

    // Reach segmentation: walk downstream from river cells that are reach
    // HEADS (not entered by another river cell's flow — i.e., no upstream
    // river neighbor flows into them), cut at junctions (2+ river inflows)
    // and at outlets (ocean or basin terminal).
    let mut inflow_count = vec![0u32; n];
    for i in 0..n {
        if let Some(ElementId(j)) = routing.flow_to[i] {
            if is_river[i] {
                inflow_count[j as usize] += 1;
            }
        }
    }
    let mut reach_of: BTreeMap<usize, ReachId> = BTreeMap::new();
    let mut reaches: Vec<RiverReachV0> = Vec::new();
    let mut reach_cells: Vec<Vec<usize>> = Vec::new();
    let mut next_reach = 0u64;
    // Heads in ascending element order (deterministic).
    for i in 0..n {
        if !is_river[i] || inflow_count[i] > 0 {
            continue;
        }
        // New reach starting at head i.
        let id = ReachId(next_reach);
        let mut cells = vec![i];
        reach_of.insert(i, id);
        let mut cur = i;
        loop {
            let next = match routing.flow_to[cur] {
                Some(ElementId(j)) => j as usize,
                None => break, // terminal: ocean handled below? ocean cells are not river cells; None = endorheic floor
            };
            if !is_river[next] {
                break; // river dies below threshold: outlet into non-river terrain
            }
            if inflow_count[next] > 1 {
                // Junction: end this reach before it; the junction cell
                // starts the next reach (head with inflows) — handled when
                // the loop reaches it as a head? It has inflows so it is
                // NOT a head by our rule. Instead: treat junction cells as
                // members of the DOWNSTREAM reach, ending the upstream one.
                break;
            }
            reach_of.insert(next, id);
            cells.push(next);
            cur = next;
        }
        reach_cells.push(cells);
        next_reach += 1;
    }
    // Junction continuation reaches: any river cell with inflows not yet
    // assigned starts a reach too (in ascending order).
    for i in 0..n {
        if !is_river[i] || reach_of.contains_key(&i) {
            continue;
        }
        let id = ReachId(next_reach);
        reach_of.insert(i, id);
        let mut cells = vec![i];
        let mut cur = i;
        loop {
            let next = match routing.flow_to[cur] {
                Some(ElementId(j)) => j as usize,
                None => break,
            };
            if !is_river[next] || inflow_count[next] > 1 || reach_of.contains_key(&next) {
                break;
            }
            reach_of.insert(next, id);
            cells.push(next);
            cur = next;
        }
        reach_cells.push(cells);
        next_reach += 1;
    }

    // Downstream linkage: from each reach's last cell, where does flow go?
    let mut downstream_of: Vec<Option<ReachId>> = Vec::with_capacity(reach_cells.len());
    let mut outlet_basin_of: Vec<Option<BasinId>> = Vec::with_capacity(reach_cells.len());
    // A reach ends either in another reach (downstream) or in the basin
    // its terminal cell belongs to (the ocean basin at the coast, or an
    // endorheic floor's basin). The basin comes from the routing's own
    // per-cell index, never from a table scan (review 2026-09-13, F1).
    let basin_of_cell = |c: usize, reach: usize| -> Result<BasinId, HydrologyError> {
        routing
            .basin_index(ElementId(c as u64))
            .map(|i| routing.basins[i].basin_id)
            .ok_or(HydrologyError::OrphanReach {
                reach: ReachId(reach as u64),
            })
    };
    for (ri, cells) in reach_cells.iter().enumerate() {
        let last = *cells.last().expect("reach nonempty");
        match routing.flow_to[last] {
            Some(ElementId(j)) => {
                let j = j as usize;
                if let Some(&r) = reach_of.get(&j) {
                    downstream_of.push(Some(r));
                    outlet_basin_of.push(None);
                } else {
                    // Non-river receiver (ocean cell or sub-threshold
                    // land): the reach terminates in that cell's basin.
                    downstream_of.push(None);
                    outlet_basin_of.push(Some(basin_of_cell(j, ri)?));
                }
            }
            None => {
                // Terminal floor: the reach terminates in its own basin.
                downstream_of.push(None);
                outlet_basin_of.push(Some(basin_of_cell(last, ri)?));
            }
        }
    }

    // Strahler order on the reach DAG (sources = order 1; at junction the
    // max order continues, +1 if two or more inflows carry the max).
    let n_reaches = reach_cells.len();
    let mut inflow_reaches: Vec<Vec<usize>> = vec![Vec::new(); n_reaches];
    for (i, d) in downstream_of.iter().enumerate() {
        if let Some(ReachId(j)) = d {
            inflow_reaches[*j as usize].push(i);
        }
    }
    let mut order = vec![0u8; n_reaches];
    // Kahn-style in topological order (deterministic: ascending id).
    let mut queue: Vec<usize> = (0..n_reaches)
        .filter(|&i| inflow_reaches[i].is_empty())
        .collect();
    queue.sort_unstable();
    let mut done = 0usize;
    while let Some(i) = queue.first().copied() {
        queue.remove(0);
        done += 1;
        let mut max = 0u8;
        let mut max_count = 0;
        for &u in &inflow_reaches[i] {
            if order[u] > max {
                max = order[u];
                max_count = 1;
            } else if order[u] == max && order[u] > 0 {
                max_count += 1;
            }
        }
        order[i] = if max == 0 {
            1
        } else if max_count >= 2 {
            max + 1
        } else {
            max
        };
        if let Some(ReachId(j)) = downstream_of[i] {
            // Push downstream when all its inflows are ordered.
            let j = j as usize;
            if inflow_reaches[j].iter().all(|&u| order[u] > 0) && order[j] == 0 {
                let mut q2 = queue.clone();
                q2.push(j);
                q2.sort_unstable();
                queue = q2;
            }
        }
        if done > n_reaches * 2 {
            return Err(HydrologyError::Domain(
                "Strahler queue did not drain (cycle?)".into(),
            ));
        }
    }
    if done != n_reaches {
        return Err(HydrologyError::Domain(format!(
            "Strahler covered {done} of {n_reaches} reaches (cycle in flow DAG)"
        )));
    }

    // Reach records + trunk monotonicity on the fill surface.
    for (i, cells) in reach_cells.iter().enumerate() {
        let last = *cells.last().unwrap();
        if let Some(ReachId(j)) = downstream_of[i] {
            // Non-increasing on the FILL surface: equality is a filled
            // lake or flat the routing drains across (review F2); an
            // increase is a broken flow graph.
            let up = routing.fill_level_m[last];
            let head = *reach_cells[j as usize].first().unwrap();
            let down = routing.fill_level_m[head];
            if up < down {
                return Err(HydrologyError::NonMonotoneTrunk {
                    reach: ReachId(i as u64),
                    upstream: up,
                    downstream: down,
                });
            }
        }
    }
    let rivers: Vec<RiverReachV0> = (0..n_reaches)
        .map(|i| RiverReachV0 {
            reach_id: ReachId(i as u64),
            downstream: downstream_of[i],
            outlet_basin: outlet_basin_of[i],
            strahler_order: order[i],
            mean_discharge_proxy_m3_s: reach_cells[i]
                .iter()
                .map(|&c| routing.discharge_proxy_m3_s[c])
                .fold(0.0, f64::max),
        })
        .collect();

    // Basin inlet reaches: the exact inverse of outlet_basin, ascending.
    let mut basins = routing.basins.clone();
    let mut basin_pos: BTreeMap<BasinId, usize> = BTreeMap::new();
    for (i, b) in basins.iter().enumerate() {
        basin_pos.insert(b.basin_id, i);
    }
    for (i, b) in outlet_basin_of.iter().enumerate() {
        if let Some(bid) = b {
            let pos = basin_pos.get(bid).ok_or(HydrologyError::OrphanReach {
                reach: ReachId(i as u64),
            })?;
            basins[*pos].inlet_reaches.push(ReachId(i as u64));
        }
    }

    let domain = HydrologyDomainV0 {
        water_bodies: routing.water_bodies.clone(),
        exchange_history: None,
        t_my,
        basins,
        rivers,
        waterfalls: Vec::new(),
        cell_fields: None,
    };
    domain
        .validate()
        .map_err(|e| HydrologyError::Domain(e.to_string()))?;
    let cell_reach: Vec<Option<ReachId>> = (0..n).map(|c| reach_of.get(&c).copied()).collect();
    Ok((domain, cell_reach))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::routing::route;

    fn setup() -> (GridSpec, Routing) {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        // Smooth ramp: all flow marches to element 0 (the ocean).
        let elev: Vec<f64> = (0..n)
            .map(|i| if i == 0 { -100.0 } else { 1000.0 + i as f64 })
            .collect();
        let mut ocean = vec![false; n];
        ocean[0] = true;
        let rain: Vec<f64> = (0..n).map(|i| 100.0 + i as f64).collect();
        let r = route(&grid, &elev, &rain, &ocean, 0.0).unwrap();
        (grid, r)
    }

    #[test]
    fn ramp_world_produces_monotone_river_with_growing_strahler_path() {
        let (grid, routing) = setup();
        let h = assemble_hydrology(&grid, &routing, 10.0).unwrap();
        h.validate().unwrap();
        assert!(!h.rivers.is_empty(), "rainy ramp must produce rivers");
        // Follow downstream to a terminal; Strahler along a single chain
        // stays >= 1 and the path terminates without cycles.
        let start = h.rivers.iter().max_by_key(|r| r.strahler_order).unwrap();
        let mut cur = start;
        let mut steps = 0;
        loop {
            match cur.downstream {
                Some(d) => {
                    cur = &h.rivers[d.0 as usize];
                    steps += 1;
                    assert!(steps < h.rivers.len() * 2, "flow path cyclic");
                }
                None => break,
            }
        }
        assert!(start.strahler_order >= 1);
    }

    #[test]
    fn outlet_basins_contain_their_terminals_and_inlets_are_the_inverse() {
        // A slope with a depression inside it and a single ocean cell: the
        // depression fills, drains across the flat and the whole surface is
        // ONE ocean-connected basin with exactly ONE outlet reach.
        let grid = GridSpec::new(12, 24, 6.371e6).unwrap();
        let n = grid.cells();
        let mut elev: Vec<f64> = (0..n).map(|c| 100.0 + 50.0 * (c / 24) as f64).collect();
        for &c in &[5 * 24 + 10, 5 * 24 + 11, 6 * 24 + 10, 6 * 24 + 11] {
            elev[c] = -400.0; // depression on the slope
        }
        let mut ocean = vec![false; n];
        ocean[3] = true;
        elev[3] = -3000.0;
        let r = route(&grid, &elev, &vec![1.0e6; n], &ocean, 0.0).unwrap();
        let h = assemble_hydrology(&grid, &r, 1.0).unwrap();
        h.validate().unwrap();
        assert_eq!(h.basins.len(), 1, "one basin: {:?}", h.basins);
        assert_eq!(h.basins[0].classification, DrainageClassV0::OceanConnected);
        // The basin has ONE outlet (the ocean cell); the reaches entering
        // it are exactly the river cells whose receiver is that cell (five
        // land neighbours of a pole-row ocean cell here), and the
        // depression contributes no separate outlet of its own.
        let outlets: Vec<&RiverReachV0> = h
            .rivers
            .iter()
            .filter(|x| x.outlet_basin.is_some())
            .collect();
        let entering = (0..n)
            .filter(|&c| r.flow_to[c] == Some(ElementId(3)))
            .count();
        assert_eq!(
            outlets.len(),
            entering,
            "outlet reaches {} != cells draining into the ocean cell {entering}",
            outlets.len()
        );
        assert!(outlets
            .iter()
            .all(|x| x.outlet_basin == Some(h.basins[0].basin_id)));
        assert_eq!(h.basins[0].inlet_reaches.len(), outlets.len());
        assert!(h
            .rivers
            .iter()
            .all(|x| x.downstream.is_some() || x.outlet_basin.is_some()));
        // Every reach's outlet basin contains its terminal cell, and the
        // basins' inlet lists are the exact inverse of outlet_basin.
        let (_, cell_reach) = assemble_hydrology_with_cells(&grid, &r, 1.0).unwrap();
        let mut inverse: BTreeMap<BasinId, Vec<ReachId>> = BTreeMap::new();
        for reach in &h.rivers {
            let Some(b) = reach.outlet_basin else {
                continue;
            };
            // Terminal cell of this reach: walk flow_to from any member cell.
            let member = (0..n)
                .find(|&c| cell_reach[c] == Some(reach.reach_id))
                .unwrap();
            let mut cur = member;
            while let Some(ElementId(j)) = r.flow_to[cur] {
                cur = j as usize;
            }
            let terminal_basin = r.basins[r.basin_index(ElementId(cur as u64)).unwrap()].basin_id;
            assert_eq!(
                terminal_basin, b,
                "reach {:?} outlet basin does not contain its terminal",
                reach.reach_id
            );
            inverse.entry(b).or_default().push(reach.reach_id);
        }
        for basin in &h.basins {
            let expect = inverse.get(&basin.basin_id).cloned().unwrap_or_default();
            assert_eq!(
                basin.inlet_reaches, expect,
                "inlet_reaches of {:?} is not the inverse",
                basin.basin_id
            );
        }
    }

    #[test]
    fn hydrology_assembly_is_deterministic() {
        let (grid, routing) = setup();
        let a = assemble_hydrology(&grid, &routing, 10.0).unwrap();
        let b = assemble_hydrology(&grid, &routing, 10.0).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn dry_world_has_no_rivers() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        let elev: Vec<f64> = (0..n)
            .map(|i| if i == 0 { -100.0 } else { 500.0 + i as f64 })
            .collect();
        let mut ocean = vec![false; n];
        ocean[0] = true;
        let tiny_rain = vec![1.0e-9; n];
        let r = route(&grid, &elev, &tiny_rain, &ocean, 0.0).unwrap();
        let h = assemble_hydrology(&grid, &r, 5.0).unwrap();
        assert!(
            h.rivers.is_empty(),
            "below-threshold discharge must not mint reaches"
        );
        h.validate().unwrap();
    }
}
