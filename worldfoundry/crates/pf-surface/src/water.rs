//! Preview water membership: largest connected sub-sea component by physical
//! area is the prescribed ocean. Other depressions fill to spill. No water
//! inventory, salinity, persistent body lineage or unresolved straits inferred.
use crate::routing::RoutingError;
use pf_contracts::domain::{WaterBodyKindV0, WaterBodyV0};
use pf_contracts::spatial::ElementId;
use pf_grid::GridSpec;

pub const DRY: u32 = u32::MAX;
pub const WET_EPS_M: f64 = 1e-9;
fn neighbours(grid: &GridSpec, c: usize) -> Result<Vec<usize>, RoutingError> {
    let (i, j) = grid
        .ij(c as u64)
        .map_err(|e| RoutingError::Grid(e.to_string()))?;
    Ok(grid
        .neighbors(i, j)
        .map_err(|e| RoutingError::Grid(e.to_string()))?
        .into_iter()
        .map(|(i, j, _)| i * grid.nlon() + j)
        .collect())
}
pub fn largest_ocean(grid: &GridSpec, z: &[f64], sea: f64) -> Result<Vec<bool>, RoutingError> {
    if z.len() != grid.cells() || !sea.is_finite() || z.iter().any(|x| !x.is_finite()) {
        return Err(RoutingError::Grid(
            "invalid ocean classification inputs".into(),
        ));
    }
    let mut seen = vec![false; z.len()];
    let mut best = Vec::new();
    let mut best_area = 0.;
    for seed in 0..z.len() {
        if seen[seed] || z[seed] >= sea {
            continue;
        }
        let mut cells = vec![seed];
        seen[seed] = true;
        let mut head = 0;
        let mut area = 0.;
        while head < cells.len() {
            let c = cells[head];
            head += 1;
            area += grid
                .cell_area_m2(c / grid.nlon())
                .map_err(|e| RoutingError::Grid(e.to_string()))?;
            for nb in neighbours(grid, c)? {
                if !seen[nb] && z[nb] < sea {
                    seen[nb] = true;
                    cells.push(nb);
                }
            }
        }
        // Seeds visited ascending: exact-area ties retain smallest seed.
        if area > best_area {
            best_area = area;
            best = cells;
        }
    }
    if best.is_empty() {
        return Err(RoutingError::NoOutlet {
            why: "no sub-sea component for prescribed ocean reservoir".into(),
        });
    }
    let mut ocean = vec![false; z.len()];
    for c in best {
        ocean[c] = true;
    }
    Ok(ocean)
}
pub fn bodies(
    grid: &GridSpec,
    z: &[f64],
    fill: &[f64],
    ocean: &[bool],
    flow: &[Option<ElementId>],
) -> Result<(Vec<u32>, Vec<WaterBodyV0>), RoutingError> {
    let mut ids = vec![DRY; z.len()];
    let mut out = Vec::new();
    for seed in 0..z.len() {
        if ids[seed] != DRY || (!ocean[seed] && fill[seed] <= z[seed] + WET_EPS_M) {
            continue;
        }
        let id = u32::try_from(out.len())
            .map_err(|_| RoutingError::Grid("too many water bodies".into()))?;
        let mut cells = vec![seed];
        ids[seed] = id;
        let mut head = 0;
        let mut area = 0.;
        while head < cells.len() {
            let c = cells[head];
            head += 1;
            area += grid
                .cell_area_m2(c / grid.nlon())
                .map_err(|e| RoutingError::Grid(e.to_string()))?;
            for nb in neighbours(grid, c)? {
                if ids[nb] == DRY
                    && ocean[nb] == ocean[seed]
                    && fill[nb] == fill[seed]
                    && (ocean[nb] || fill[nb] > z[nb] + WET_EPS_M)
                {
                    ids[nb] = id;
                    cells.push(nb);
                }
            }
        }
        let mut exits = Vec::new();
        for c in cells {
            if let Some(j) = flow[c] {
                if ids[j.0 as usize] != id {
                    exits.push([ElementId(c as u64), j]);
                }
            }
        }
        exits.sort();
        out.push(WaterBodyV0 {
            body_id: id,
            kind: if ocean[seed] {
                WaterBodyKindV0::Marine
            } else {
                WaterBodyKindV0::Lake
            },
            seed_element: ElementId(seed as u64),
            water_surface_m: fill[seed],
            area_m2: area,
            outlet_edges: exits,
        });
    }
    Ok((ids, out))
}
