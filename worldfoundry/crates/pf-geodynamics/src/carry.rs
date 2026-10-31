//! SEED-CARRY-CC (lane CC, 2026-10-12): the plate-carry law selector and
//! the patch test of `front-rescue-carry-v1`.
//!
//! The legacy carry (`semi-lagrangian-nearest-material-v0`, A-ADVECTION,
//! `kinematic.rs`) lets each plate claim a cell when the GLOBAL nearest
//! material parcel to the back-traced point is its own. That is a single
//! Voronoi partition over every plate's parcels, so two converging
//! continental fronts can never both claim one cell: front parcels vanish
//! silently into the remap residual instead of producing an overlap, and the
//! continent-continent contact present per step depends on dt (lane OG
//! falsifier F2, `unpublished-design-record` section 4).
//!
//! `front-rescue-carry-v1` (opt-in; pre-registered in
//! `unpublished-design-record`, v2 entry 21:40 EDT; no
//! free constant):
//! 1. legacy claims are computed for every cell, unchanged;
//! 2. a continental parcel that the legacy resolution would carry to NO
//!    cell is RESCUED: moved by its plate's rotation, it claims the one cell
//!    (3x3 around its moved position) where its plate has no legacy claim
//!    and whose back-traced point lies inside the parcel's cell-shaped patch
//!    (`in_cell_patch`); nearest back-traced point wins, ties to the lowest
//!    cell; one rescue per plate per cell (nearest parcel, ties lowest
//!    source). A rescued parcel is moved, never copied;
//! 3. the legacy overriding rule resolves the cell (continent over ocean is
//!    subduction as before; continent-continent is collision);
//! 4. every continental loser of a continent-continent overlap that is
//!    carried to no destination is consumed once, and exactly its crust
//!    volume and sediment become the collision supply on the
//!    `rate-limited-belt-v1` belt (no volume created or silently dropped by
//!    collision).
//!
//! Lane CC's first version (`contact-claim-carry-v1`: a plate's own nearest
//! continental parcel claimed any covering cell at a continent-continent
//! contact) failed its pre-registered synthetic check because the claimed
//! parcel was usually ALSO carried by the legacy rule (a duplicate); see
//! FALSIFIER.txt. It is not shipped.

use crate::kinematic::KinematicError;
use crate::rotation::great_circle_rad;
use pf_grid::GridSpec;

/// Which carry runs (mirrors `pf_contracts::laws::PlateCarryLawV1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CarryLaw {
    /// Legacy `semi-lagrangian-nearest-material-v0` (bit for bit).
    #[default]
    SemiLagrangianNearestMaterialV0,
    /// `front-rescue-carry-v1` (requires `rate-limited-belt-v1`).
    FrontRescueV1,
    /// `continental-parcel-carry-v1` (SEED-CARRY-CONSERVE, lane CV,
    /// `conserve.rs`; requires `rate-limited-belt-v1` and the parcel layer,
    /// so it runs only through `kinematic::advance_with_parcels`).
    ContinentalParcelsV1,
}

/// True when `q` lies inside the cell-shaped patch of the parcel at `p`
/// whose source cell is `s`: |north offset| <= dlat/2 and |east offset| <=
/// cos(lat_s) dlon/2 (unit-sphere radians) in the tangent plane at `p`,
/// with the extents of cell `s` (lat_s its centre latitude). Inclusive. At
/// a parcel exactly on a pole (no east direction) the patch is the cap of
/// radius dlat/2.
pub fn in_cell_patch(
    grid: &GridSpec,
    p: [f64; 3],
    s: usize,
    q: [f64; 3],
) -> Result<bool, KinematicError> {
    let dlat = std::f64::consts::PI / grid.nlat() as f64;
    let dlon = 2.0 * std::f64::consts::PI / grid.nlon() as f64;
    let lat_s = grid
        .center_lat_deg(s / grid.nlon())
        .map_err(|e| KinematicError::Grid(e.to_string()))?
        .to_radians();
    let en = p[0].hypot(p[1]);
    if en < 1.0e-12 {
        return Ok(great_circle_rad(p, q)? <= 0.5 * dlat);
    }
    let east = [-p[1] / en, p[0] / en, 0.0];
    let north = [
        p[1] * east[2] - p[2] * east[1],
        p[2] * east[0] - p[0] * east[2],
        p[0] * east[1] - p[1] * east[0],
    ];
    let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let de = d[0] * east[0] + d[1] * east[1] + d[2] * east[2];
    let dn = d[0] * north[0] + d[1] * north[1] + d[2] * north[2];
    Ok(dn.abs() <= 0.5 * dlat && de.abs() <= 0.5 * dlon * lat_s.cos())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rotation::lat_lon_to_unit;

    #[test]
    fn patch_matches_the_cell_extents() {
        let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
        // Cell (row 9, col 18): centre lat 5, lon 5 (10-degree cells).
        let s = 9 * 36 + 18;
        let lat = grid.center_lat_deg(9).unwrap();
        let lon = grid.center_lon_deg(18).unwrap();
        let p = lat_lon_to_unit(lat, lon).unwrap();
        assert!(in_cell_patch(&grid, p, s, p).unwrap());
        // Just inside / outside the east edge (half width 5 deg * cos 5).
        let inside = lat_lon_to_unit(lat, lon + 4.9).unwrap();
        let outside = lat_lon_to_unit(lat, lon + 5.2).unwrap();
        assert!(in_cell_patch(&grid, p, s, inside).unwrap());
        assert!(!in_cell_patch(&grid, p, s, outside).unwrap());
        // North edge.
        assert!(in_cell_patch(&grid, p, s, lat_lon_to_unit(lat + 4.9, lon).unwrap()).unwrap());
        assert!(!in_cell_patch(&grid, p, s, lat_lon_to_unit(lat + 5.1, lon).unwrap()).unwrap());
        // A neighbouring centroid is never inside.
        let nb = lat_lon_to_unit(lat, lon + 10.0).unwrap();
        assert!(!in_cell_patch(&grid, p, s, nb).unwrap());
    }
}
