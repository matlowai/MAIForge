//! SEED-CARRY-CC (lane CC, 2026-10-12): ignored kinematics-only diagnostic
//! of continent-continent contacts on a real seed (toy preset, continental
//! bodies, rate-limited-belt-v1 claims; no surface processes, which do not
//! change plate or crust-kind evolution). Env: CC_SEED (default 240907),
//! CC_LEVEL (default 0), CC_LAW (legacy | cc). Prints per-My totals for
//! dt 2 / 1 / 0.5: convergent contact edges per step, sliver births,
//! continental parcels dropped, the carry's net continental volume loss in
//! the convergent contact zone, and the belt input.

use pf_contracts::domain::{EulerPoleV0, PlateId};
use pf_contracts::params::{
    ContinentalBodiesBlock, FidelityPreset, MantleBlock, ScenarioPreset, SeedPath,
};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::carry::CarryLaw;
use pf_geodynamics::kinematic::{advance_with_carry, CrustState};
use pf_geodynamics::orogen::BeltParams;
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;
use std::collections::BTreeMap;

#[test]
#[ignore]
fn carry_cc_seed_diag() {
    let seed: u64 = std::env::var("CC_SEED")
        .map(|v| v.parse().unwrap())
        .unwrap_or(240907);
    let level: u8 = std::env::var("CC_LEVEL")
        .map(|v| v.parse().unwrap())
        .unwrap_or(0);
    let law = match std::env::var("CC_LAW").unwrap_or("legacy".into()).as_str() {
        "cc" => CarryLaw::FrontRescueV1,
        _ => CarryLaw::SemiLagrangianNearestMaterialV0,
    };
    let fid = if level == 0 {
        FidelityPreset::Preview
    } else {
        FidelityPreset::Standard
    };
    let p = ScenarioPreset::toy().profile(fid, MantleBlock::default());
    let grid = GridSpec::level(level, 6.371e6).unwrap();
    let n = grid.cells();
    let nlon = grid.nlon();
    let sp = SeedPath {
        seed,
        derivation: "kinematic-preview".into(),
        salt: "toy".into(),
    };
    let seeding = seed_plates(&grid, &p, &sp).unwrap();
    let prior = UpliftPrior::preview_v0();
    let b = ContinentalBodiesBlock::v0();
    let (st0, _) = CrustState::seed_bodies(
        &grid,
        &seeding,
        &prior,
        b.resolved_area_fraction(&p.tectonics),
        b.rifted_margin_width_m,
    )
    .unwrap();
    let belt = BeltParams {
        shear_zone_dip_deg: 30.0,
        channel_viscosity_pa_s: 1e19,
        channel_reference_thickness_m: 35_000.0,
        delamination_thickness_m: 70_000.0,
        arc_productivity_m: 500.0,
        arc_belt_width_m: 100_000.0,
        gravity_m_s2: 9.81,
        crust_density_kg_m3: 2800.0,
        mantle_density_kg_m3: 3300.0,
    };
    let poles: BTreeMap<PlateId, EulerPoleV0> = seeding
        .plates
        .iter()
        .filter_map(|r| r.euler_pole.map(|e| (r.plate_id, e)))
        .collect();
    let max_rate = poles
        .values()
        .map(|e| e.rate_deg_per_my.abs())
        .fold(0.0, f64::max);
    let dlat_deg = 180.0 / grid.nlat() as f64;
    let areas: Vec<f64> = (0..n)
        .map(|c| grid.cell_area_m2(c / nlon).unwrap())
        .collect();
    let cents: Vec<[f64; 3]> = (0..n)
        .map(|c| {
            pf_geodynamics::rotation::lat_lon_to_unit(
                grid.center_lat_deg(c / nlon).unwrap(),
                grid.center_lon_deg(c % nlon).unwrap(),
            )
            .unwrap()
        })
        .collect();
    println!(
        "seed {seed} L{level} law {law:?} max plate rate {max_rate} deg/My, cell {dlat_deg} deg"
    );
    for dt in [2.0, 1.0, 0.5] {
        let rings = 2 + ((2.0 * max_rate * dt) / dlat_deg).ceil() as usize;
        let mut st = st0.clone();
        let steps = (100.0 / dt) as usize;
        let (mut edges, mut slivers, mut drop_m3, mut zone_net_m3, mut belt_in, mut consumed) =
            (0u64, 0u64, 0.0, 0.0, 0.0, 0.0);
        let mut zone_neg_steps = 0u32;
        let mut ridge = 0u64;
        let mut sub = 0u64;
        let mut bins = vec![0u64; 10];
        for step in 0..steps {
            // Convergent continent-continent edges at step start, and the
            // contact zone (rings around them).
            let mut zone = vec![false; n];
            for c in 0..n {
                if st.kind[c] != CrustKindCodeV0::Continental {
                    continue;
                }
                let (i, j) = grid.ij(c as u64).unwrap();
                for (ti, tj, nb) in grid.neighbors(i, j).unwrap() {
                    if nb.di.abs() + nb.dj.abs() != 1 {
                        continue;
                    }
                    let k = ti * nlon + tj;
                    if st.plate[k] == st.plate[c] || st.kind[k] != CrustKindCodeV0::Continental {
                        continue;
                    }
                    let raw = [
                        cents[c][0] + cents[k][0],
                        cents[c][1] + cents[k][1],
                        cents[c][2] + cents[k][2],
                    ];
                    let l = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
                    let mid = [raw[0] / l, raw[1] / l, raw[2] / l];
                    let vn = pf_geodynamics::boundaries::normal_rate_at(
                        mid,
                        cents[c],
                        cents[k],
                        &poles[&st.plate[c]],
                        &poles[&st.plate[k]],
                    )
                    .unwrap();
                    if vn > 0.0 {
                        zone[c] = true;
                        zone[k] = true;
                    }
                }
            }
            for _ in 0..rings {
                let cur = zone.clone();
                for c in 0..n {
                    if cur[c] {
                        let (i, j) = grid.ij(c as u64).unwrap();
                        for (ti, tj, _) in grid.neighbors(i, j).unwrap() {
                            zone[ti * nlon + tj] = true;
                        }
                    }
                }
            }
            let (next, l, _, src) =
                advance_with_carry(&grid, &st, &seeding.plates, &prior, dt, Some(&belt), law)
                    .unwrap();
            edges += l.collision_supply_edges;
            bins[((step as f64 * dt) / 10.0) as usize] += l.collision_supply_edges;
            slivers += l.carry.sliver_births;
            drop_m3 += l.carry.continental_dropped_m3;
            belt_in += l.belt_in_m3;
            consumed += l.carry.cc_consumed_m3;
            ridge += l.ridge_cells;
            sub += l.subduction_cells;
            // Net continental volume the carry removes in the zone:
            // sources in the zone minus their carried copies.
            let mut net = 0.0;
            for s in 0..n {
                if zone[s] && st.kind[s] == CrustKindCodeV0::Continental {
                    net += st.thickness_m[s] * areas[s];
                }
            }
            for d in 0..n {
                let s = src[d];
                if s != u32::MAX {
                    let s = s as usize;
                    if zone[s] && st.kind[s] == CrustKindCodeV0::Continental {
                        net -= st.thickness_m[s] * areas[d];
                    }
                }
            }
            if net < 0.0 {
                zone_neg_steps += 1;
            }
            zone_net_m3 += net;
            if std::env::var("CC_TRACE").is_ok()
                && (st.t_my < 12.0 || (next.t_my % 20.0).abs() < 1e-9)
            {
                let mut pa: BTreeMap<u32, f64> = BTreeMap::new();
                for c in 0..n {
                    if next.kind[c] == CrustKindCodeV0::Continental {
                        *pa.entry(next.plate[c].0).or_insert(0.0) += areas[c] / 1e12;
                    }
                }
                // Continental cells with a continental cell of another plate
                // at ring distance 1 (4- or 8-adjacent), 2 and 3, and the
                // non-continental cells between them (slivers): counts.
                let mut ring_hits = [0u64; 3];
                for c in 0..n {
                    if next.kind[c] != CrustKindCodeV0::Continental {
                        continue;
                    }
                    let (i, j) = grid.ij(c as u64).unwrap();
                    let mut found = 0usize;
                    'r: for r in 1..=3i64 {
                        for di in -r..=r {
                            for dj in -r..=r {
                                if di.abs().max(dj.abs()) != r {
                                    continue;
                                }
                                let ti = i as i64 + di;
                                if ti < 0 || ti >= grid.nlat() as i64 {
                                    continue;
                                }
                                let tj = (j as i64 + dj).rem_euclid(nlon as i64);
                                let k = ti as usize * nlon + tj as usize;
                                if next.kind[k] == CrustKindCodeV0::Continental
                                    && next.plate[k] != next.plate[c]
                                {
                                    found = r as usize;
                                    break 'r;
                                }
                            }
                        }
                    }
                    if found > 0 {
                        ring_hits[found - 1] += 1;
                    }
                }
                let ph: Vec<f64> = (0..n)
                    .map(|c| {
                        pf_geodynamics::rotation::great_circle_rad(next.material_pos[c], cents[c])
                            .unwrap()
                            .to_degrees()
                            / dlat_deg
                    })
                    .collect();
                let big = ph.iter().filter(|&&x| x > 1.0).count();
                let huge = ph.iter().filter(|&&x| x > 3.0).count();
                println!("    cont cells with other-plate continent at ring 1/2/3: {:?}; phase > 1 cell: {big}, > 3 cells: {huge}", ring_hits);
                println!(
                    "  t {:.1} edges {} drops {} births {} sliver {} sub {} coll {} cont area by plate (1e12 m2) {:?}",
                    next.t_my, l.collision_supply_edges, l.carry.continental_dropped_parcels, l.ridge_cells,
                    l.carry.sliver_births, l.subduction_cells, l.collision_cells,
                    pa.iter().map(|(k, v)| format!("{k}:{v:.2}")).collect::<Vec<_>>()
                );
            }
            st = next;
        }
        let per_step_bins: Vec<String> = bins
            .iter()
            .map(|e| format!("{:.0}", *e as f64 / (10.0 / dt)))
            .collect();
        println!(
            "dt {dt} rings {rings}: edges/step {:.1} [by 10 My: {}] | per My: sliver births {:.2}, ridge cells {:.1}, subduction cells {:.1}, cont dropped {:.3e} m3, zone net loss {:.3e} m3 (negative steps {zone_neg_steps}/{steps}), belt_in {:.3e} m3, consumed {:.3e} m3",
            edges as f64 / steps as f64,
            per_step_bins.join(" "),
            slivers as f64 / 100.0,
            ridge as f64 / 100.0,
            sub as f64 / 100.0,
            drop_m3 / 100.0,
            zone_net_m3 / 100.0,
            belt_in / 100.0,
            consumed / 100.0
        );
    }
}
