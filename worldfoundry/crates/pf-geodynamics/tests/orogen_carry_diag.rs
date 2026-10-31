//! SEED-OROGEN (lane OG, 2026-10-12): ignored diagnostic of how the
//! semi-Lagrangian carry resolves head-on continent-continent convergence at
//! different step lengths. Two all-continental plates (west / east
//! hemispheres) rotate about the same pole axis in opposite senses, so they
//! converge along one meridian and diverge along the antipodal one. For each
//! dt the run reports, per My: collision overlaps, ridge births, subduction
//! overlaps, the step-start convergent interface length, and the crust
//! volume change split into named ledger terms and the remap residual.
//! An exact kinematic model shortens the continental area by v * L per My
//! and conserves crust volume in collision (only ridge births add volume).

use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_geodynamics::kinematic::{advance_with, CrustState};
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;

fn plates(rate: f64) -> Vec<PlateRecordV0> {
    let rec = |id: u32, r: f64| PlateRecordV0 {
        plate_id: PlateId(id),
        parent_plate: None,
        birth_t_my: 0.0,
        euler_pole: Some(EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: r,
        }),
        retirement: None,
    };
    vec![rec(0, rate), rec(1, -rate)]
}

#[test]
#[ignore]
fn orogen_carry_head_on_collision() {
    let (nlat, nlon) = (
        std::env::var("OG_NLAT")
            .map(|v| v.parse().unwrap())
            .unwrap_or(72usize),
        std::env::var("OG_NLON")
            .map(|v| v.parse().unwrap())
            .unwrap_or(144usize),
    );
    let grid = GridSpec::new(nlat, nlon, 6.371e6).unwrap();
    let n = grid.cells();
    let prior = UpliftPrior::preview_v0();
    let areas: Vec<f64> = (0..n)
        .map(|c| grid.cell_area_m2(c / nlon).unwrap())
        .collect();
    let rate = 0.5; // deg/My each: 1 deg/My closing at the equator (~111 km/My)
    let pl = plates(rate);
    for dt in [2.0, 1.0, 0.5, 0.25] {
        // Plate 0 = lon in [-90, 90) moves east (+rate about the north
        // pole), plate 1 = the other hemisphere moves west: they converge at
        // lon +90 and diverge at lon -90.
        let mut st = CrustState {
            t_my: 0.0,
            plate: (0..n)
                .map(|c| {
                    let lon = grid.center_lon_deg(c % nlon).unwrap();
                    PlateId(if (-90.0..90.0).contains(&lon) { 0 } else { 1 })
                })
                .collect(),
            kind: vec![CrustKindCodeV0::Continental; n],
            birth_t_my: vec![-2000.0; n],
            thickness_m: vec![35_000.0; n],
            sediment_m: vec![0.0; n],
            eroded_total_m: vec![0.0; n],
            arc_added_m: vec![0.0; n],
            shortening_added_m: vec![0.0; n],
            boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
            normal_rate_m_per_my: vec![0.0; n],
            material_pos: (0..n)
                .map(|c| {
                    pf_geodynamics::rotation::lat_lon_to_unit(
                        grid.center_lat_deg(c / nlon).unwrap(),
                        grid.center_lon_deg(c % nlon).unwrap(),
                    )
                    .unwrap()
                })
                .collect(),
        };
        let v0: f64 = st.thickness_m.iter().zip(&areas).map(|(h, a)| h * a).sum();
        let steps = (20.0 / dt) as usize;
        let (mut coll, mut ridge, mut sub, mut resid, mut ridge_v) = (0u64, 0u64, 0u64, 0.0, 0.0);
        let (mut cont_area_loss, mut stacked) = (0.0, 0.0);
        for _ in 0..steps {
            let a0: f64 = (0..n)
                .filter(|&c| st.kind[c] == CrustKindCodeV0::Continental)
                .map(|c| areas[c])
                .sum();
            let (next, l, _, src) = advance_with(&grid, &st, &pl, &prior, dt, None).unwrap();
            if std::env::var("OG_DEBUG").is_ok() {
                let moved = src
                    .iter()
                    .enumerate()
                    .filter(|(c, s)| **s as usize != *c)
                    .count();
                let eq = (0..n)
                    .filter(|&c| c / nlon == nlat / 2)
                    .map(|c| format!("{}", st.plate[c].0))
                    .collect::<String>();
                println!(
                    "  t {:.1} moved-source cells {moved} coll {} births {} eq-row plates {eq}",
                    st.t_my, l.collision_cells, l.ridge_cells
                );
                let row = nlat / 2;
                let offs: Vec<String> = (0..nlon)
                    .map(|j| {
                        let c = row * nlon + j;
                        let sidx = src[c];
                        if sidx == u32::MAX {
                            "B".into()
                        } else {
                            let d = (sidx as i64 % nlon as i64) - j as i64;
                            format!(
                                "{}",
                                d.rem_euclid(nlon as i64) as i64
                                    - if d.rem_euclid(nlon as i64) > (nlon as i64) / 2 {
                                        nlon as i64
                                    } else {
                                        0
                                    }
                            )
                        }
                    })
                    .collect();
                let (lat, lon) =
                    pf_geodynamics::rotation::unit_to_lat_lon(st.material_pos[row * nlon + 10])
                        .unwrap();
                println!(
                    "    offsets {} | mat(10) lon {lon:.3} cent lon {:.3}",
                    offs.join(","),
                    grid.center_lon_deg(10).unwrap()
                );
                let _ = lat;
            }
            coll += l.collision_cells;
            ridge += l.ridge_cells;
            sub += l.subduction_cells;
            resid += l.residual_m3();
            ridge_v += l.ridge_created_m3;
            stacked += l.wedge_in_m3;
            let a1: f64 = (0..n)
                .filter(|&c| next.kind[c] == CrustKindCodeV0::Continental)
                .map(|c| areas[c])
                .sum();
            cont_area_loss += a0 - a1;
            st = next;
        }
        let v1: f64 = st.thickness_m.iter().zip(&areas).map(|(h, a)| h * a).sum();
        // Exact: convergent boundary = meridian at +90 deg, length pi R; closing
        // speed at latitude phi = 2 * rate * cos(phi) * R (rad) -> area rate
        // = integral 2 rate cos(phi) R * R dphi over [-pi/2, pi/2] = 4 rate R^2.
        let r = 6.371e6f64;
        let exact_area_per_my = 4.0 * rate.to_radians() * r * r;
        println!(
            "dt {dt}: per My: collision overlaps {:.1}, ridge births {:.1}, subduction {:.1}; continental area loss {:.3e} m2/My (exact closing {:.3e}); \
             volume change {:.3e} m3 = ridge {:.3e} + remap residual {:.3e} (wedge in {:.3e}); max thickness {:.0} m",
            coll as f64 / 20.0,
            ridge as f64 / 20.0,
            sub as f64 / 20.0,
            cont_area_loss / 20.0,
            exact_area_per_my,
            v1 - v0,
            ridge_v,
            resid,
            stacked,
            st.thickness_m.iter().cloned().fold(0.0, f64::max)
        );
    }
}
