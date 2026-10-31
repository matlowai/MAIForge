//! SEED-ROBUST (lane R5, 2026-10-12) measurement driver: replays the
//! tectonic-only `advance` history of a toy-preset seed and, at the first
//! step that refuses (or at a requested step), writes the input plate
//! field, the per-cell claimant sets recomputed independently, and every
//! subduction overlap whose pair interface is missing or beyond the locality
//! bound. Plate geometry does not depend on erosion (erosion changes only
//! thickness and sediment; claims depend on plate, kind, birth and material
//! position), so this replay reaches the same refusal as `pf-lab simulate`.
//!
//! Ignored by default (it is a measurement, not a gate). Run with
//! `SR_SEED=15 SR_LEVEL=0 SR_MODE=legacy SR_OUT=<dir> cargo test --release
//! -p pf-geodynamics --test seed_robust_diag -- --ignored --nocapture`.

use pf_contracts::domain::{EulerPoleV0, PlateId};
use pf_contracts::params::{
    ContinentalBodiesBlock, FidelityPreset, MantleBlock, ScenarioPreset, SeedPath,
};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::kinematic::{
    advance, nearest_cell, ArcContactGeometry, CrustState, SEARCH_HALF_WIDTH,
};
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::rotation::{
    great_circle_rad, lat_lon_to_unit, rotate_by_euler, unit_to_lat_lon,
};
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;
use std::collections::BTreeMap;
use std::io::Write;

fn centroids(grid: &GridSpec) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for i in 0..grid.nlat() {
        let lat = grid.center_lat_deg(i).unwrap();
        for j in 0..grid.nlon() {
            out.push(lat_lon_to_unit(lat, grid.center_lon_deg(j).unwrap()).unwrap());
        }
    }
    out
}

/// Independent re-implementation of the kernel's metric-window search.
fn nearest_material(grid: &GridSpec, mp: &[[f64; 3]], q: [f64; 3]) -> usize {
    let centre = nearest_cell(grid, q).unwrap();
    let nlat = grid.nlat() as i64;
    let nlon = grid.nlon() as i64;
    let ci = (centre / grid.nlon()) as i64;
    let cj = (centre % grid.nlon()) as i64;
    let mut best: Option<(f64, usize)> = None;
    for di in -SEARCH_HALF_WIDTH..=SEARCH_HALF_WIDTH {
        let i = ci + di;
        if i < 0 || i >= nlat {
            continue;
        }
        let lat = grid.center_lat_deg(i as usize).unwrap();
        let cos_lat = lat.to_radians().cos().max(1.0e-9);
        let dj_half = ((SEARCH_HALF_WIDTH as f64 / cos_lat).ceil() as i64).min(nlon / 2);
        for dj in -dj_half..=dj_half {
            let j = (cj + dj).rem_euclid(nlon);
            let s = (i * nlon + j) as usize;
            let d = great_circle_rad(mp[s], q).unwrap();
            if best
                .map(|(bd, bs)| d < bd || (d == bd && s < bs))
                .unwrap_or(true)
            {
                best = Some((d, s));
            }
        }
    }
    best.unwrap().1
}

/// Detached fragments: 8-connected components (pf-grid neighbour policy)
/// of a plate that are not that plate's largest. Returns (count, cells).
fn fragments(grid: &GridSpec, plate: &[PlateId]) -> (usize, usize) {
    let n = plate.len();
    let mut lab = vec![usize::MAX; n];
    let mut comps: Vec<(PlateId, usize)> = Vec::new();
    for s in 0..n {
        if lab[s] != usize::MAX {
            continue;
        }
        let k = comps.len();
        lab[s] = k;
        let mut stack = vec![s];
        let mut size = 0;
        while let Some(c) = stack.pop() {
            size += 1;
            let (i, j) = grid.ij(c as u64).unwrap();
            for (ni, nj, _) in grid.neighbors(i, j).unwrap() {
                let nb = ni * grid.nlon() + nj;
                if lab[nb] == usize::MAX && plate[nb] == plate[s] {
                    lab[nb] = k;
                    stack.push(nb);
                }
            }
        }
        comps.push((plate[s], size));
    }
    let mut largest: BTreeMap<PlateId, usize> = BTreeMap::new();
    for (p, sz) in &comps {
        let e = largest.entry(*p).or_insert(0);
        *e = (*e).max(*sz);
    }
    let mut count = 0;
    let mut cells = 0;
    let mut seen_largest: BTreeMap<PlateId, bool> = BTreeMap::new();
    for (p, sz) in &comps {
        if *sz == largest[p] && !seen_largest.get(p).copied().unwrap_or(false) {
            seen_largest.insert(*p, true);
            continue;
        }
        count += 1;
        cells += sz;
    }
    (count, cells)
}

fn mid(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let r = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
    let n = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    [r[0] / n, r[1] / n, r[2] / n]
}

fn env(k: &str, d: &str) -> String {
    std::env::var(k).unwrap_or_else(|_| d.to_string())
}

#[test]
#[ignore]
fn seed_robust_diag() {
    let seed: u64 = env("SR_SEED", "15").parse().unwrap();
    let level: u8 = env("SR_LEVEL", "0").parse().unwrap();
    let mode = env("SR_MODE", "legacy");
    let out_dir = std::path::PathBuf::from(env("SR_OUT", "."));
    let stop_t: f64 = env("SR_STOP_T", "-1").parse().unwrap();
    std::fs::create_dir_all(&out_dir).unwrap();
    let fid = match level {
        0 => FidelityPreset::Preview,
        1 => FidelityPreset::Standard,
        _ => FidelityPreset::High,
    };
    let p = ScenarioPreset::toy().profile(fid, MantleBlock::default());
    let grid = GridSpec::level(p.fidelity.grid_level, p.planet.radius_m).unwrap();
    let sp = SeedPath {
        seed,
        derivation: "kinematic-preview".into(),
        salt: "toy".into(),
    };
    let seeding = seed_plates(&grid, &p, &sp).unwrap();
    let prior = UpliftPrior::preview_v0();
    let mut state = if mode == "bodies" {
        let b = ContinentalBodiesBlock::v0();
        CrustState::seed_bodies(
            &grid,
            &seeding,
            &prior,
            b.resolved_area_fraction(&p.tectonics),
            b.rifted_margin_width_m,
        )
        .unwrap()
        .0
    } else {
        CrustState::seed(&grid, &seeding, &prior).unwrap()
    };
    let plates = seeding.plates.clone();
    let poles: BTreeMap<PlateId, EulerPoleV0> = plates
        .iter()
        .map(|p| (p.plate_id, p.euler_pole.unwrap()))
        .collect();
    let dt = p.tectonics.step_my;
    let steps = (p.tectonics.duration_my / dt).round() as u32;
    let cents = centroids(&grid);
    let n = grid.cells();
    let cell_angle = (std::f64::consts::PI / grid.nlat() as f64)
        .hypot(2. * std::f64::consts::PI / grid.nlon() as f64);
    let bound = |a: PlateId, b: PlateId| {
        (poles[&a].rate_deg_per_my.abs() + poles[&b].rate_deg_per_my.abs()).to_radians() * dt
            + 3. * cell_angle
    };
    let mut stats = std::fs::File::create(out_dir.join("arc_ratio.tsv")).unwrap();
    writeln!(
        stats,
        "t_my\tcell\twinner\tloser\tdistance_rad\tbound_rad\tratio\tgeometry\trate_m_per_my"
    )
    .unwrap();
    let mut steplog = std::fs::File::create(out_dir.join("steps.tsv")).unwrap();
    writeln!(steplog, "t_my\tsubduction_samples\tocean_ocean_samples\tcollision_cells\tridge_cells\tnew_contacts\tfragments\tfragment_cells").unwrap();
    let mut refusal: Option<String> = None;
    let mut t_fail = None;
    for _ in 0..steps {
        if stop_t >= 0. && state.t_my >= stop_t {
            t_fail = Some(state.t_my);
            break;
        }
        match advance(&grid, &state, &plates, &prior, dt) {
            Ok((next, ledger, _)) => {
                for (k, s) in ledger.arc_normals.iter().enumerate() {
                    let c = s.overlap_element.0 as usize;
                    let (i, j) = (
                        s.input_interface[0].0 as usize,
                        s.input_interface[1].0 as usize,
                    );
                    let d = great_circle_rad(cents[c], mid(cents[i], cents[j])).unwrap();
                    let bd = bound(s.winner, s.loser);
                    writeln!(
                        stats,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{}",
                        state.t_my,
                        c,
                        s.winner.0,
                        s.loser.0,
                        d,
                        bd,
                        d / bd,
                        ledger.arc_contact_geometry[k],
                        s.normal_rate_m_per_my
                    )
                    .unwrap();
                }
                let oo = ledger
                    .arc_normals
                    .iter()
                    .filter(|s| next.kind[s.overlap_element.0 as usize] == CrustKindCodeV0::Oceanic)
                    .count();
                let nc = ledger
                    .arc_contact_geometry
                    .iter()
                    .filter(|g| **g == ArcContactGeometry::NewContactBridge)
                    .count();
                let (fr, fc) = fragments(&grid, &next.plate);
                writeln!(
                    steplog,
                    "{}\t{}\t{oo}\t{}\t{}\t{nc}\t{fr}\t{fc}",
                    next.t_my,
                    ledger.arc_normals.len(),
                    ledger.collision_cells,
                    ledger.ridge_cells
                )
                .unwrap();
                state = next;
            }
            Err(e) => {
                refusal = Some(e.to_string());
                t_fail = Some(state.t_my);
                break;
            }
        }
    }
    // Recompute claims for every cell of the refusing (or requested) step.
    let mut back: Vec<(PlateId, EulerPoleV0, EulerPoleV0)> = poles
        .iter()
        .map(|(id, e)| {
            (
                *id,
                EulerPoleV0 {
                    rate_deg_per_my: -e.rate_deg_per_my,
                    ..*e
                },
                *e,
            )
        })
        .collect();
    back.sort_by_key(|x| x.0);
    let mut claim_mask = vec![0u32; n];
    let mut winner_of = vec![u32::MAX; n];
    // Input interface edges per pair (as the kernel builds them).
    let mut iface: BTreeMap<(PlateId, PlateId), Vec<(usize, usize)>> = BTreeMap::new();
    for c in 0..n {
        let (i, j) = grid.ij(c as u64).unwrap();
        for (ni, nj, _) in grid.neighbors(i, j).unwrap() {
            let nb = ni * grid.nlon() + nj;
            if state.plate[c] < state.plate[nb] {
                iface
                    .entry((state.plate[c], state.plate[nb]))
                    .or_default()
                    .push((c, nb));
            }
        }
    }
    let mut bad = Vec::new();
    let mut overlaps = 0u64;
    for c in 0..n {
        let mut claims: Vec<(PlateId, usize)> = Vec::new();
        for (id, bp, _) in &back {
            let q = rotate_by_euler(cents[c], bp, dt).unwrap();
            let s = nearest_material(&grid, &state.material_pos, q);
            if state.plate[s] == *id {
                claims.push((*id, s));
                claim_mask[c] |= 1 << id.0;
            }
        }
        if claims.len() < 2 {
            continue;
        }
        overlaps += 1;
        let cont: Vec<_> = claims
            .iter()
            .filter(|(_, s)| state.kind[*s] == CrustKindCodeV0::Continental)
            .collect();
        let w = if let Some(w) = cont.first() {
            **w
        } else {
            let mut b = claims[0];
            for cl in &claims[1..] {
                if state.birth_t_my[cl.1] > state.birth_t_my[b.1] {
                    b = *cl;
                }
            }
            b
        };
        winner_of[c] = w.0 .0;
        for (lid, ls) in &claims {
            if *lid == w.0 {
                continue;
            }
            let collision = state.kind[*ls] == CrustKindCodeV0::Continental
                && state.kind[w.1] == CrustKindCodeV0::Continental;
            if collision {
                continue;
            }
            let key = if w.0 < *lid { (w.0, *lid) } else { (*lid, w.0) };
            let bd = bound(w.0, *lid);
            let (dist, edge) = match iface.get(&key) {
                None => (f64::INFINITY, None),
                Some(edges) => {
                    let mut best = (f64::INFINITY, None);
                    for &(a, b) in edges {
                        let d = great_circle_rad(cents[c], mid(cents[a], cents[b])).unwrap();
                        if d < best.0 {
                            best = (d, Some((a, b)));
                        }
                    }
                    best
                }
            };
            if dist > bd {
                let (lat, lon) = unit_to_lat_lon(cents[c]).unwrap();
                let srcs: Vec<String> = claims
                    .iter()
                    .map(|(pid, s)| {
                        format!(
                            "{{\"plate\":{},\"src\":{},\"kind\":{},\"birth\":{},\"src_dist_rad\":{}}}",
                            pid.0,
                            s,
                            state.kind[*s].code(),
                            state.birth_t_my[*s],
                            great_circle_rad(cents[c], cents[*s]).unwrap()
                        )
                    })
                    .collect();
                bad.push(format!(
                    "{{\"cell\":{c},\"lat\":{lat},\"lon\":{lon},\"input_plate\":{},\"winner\":{},\"loser\":{},\"distance_rad\":{},\"bound_rad\":{bd},\"edge\":{:?},\"claims\":[{}]}}",
                    state.plate[c].0,
                    w.0 .0,
                    lid.0,
                    if dist.is_finite() { dist.to_string() } else { "null".into() },
                    edge.map(|(a, b)| vec![a, b]).unwrap_or_default(),
                    srcs.join(",")
                ));
            }
        }
    }
    let plate_u32: Vec<u8> = state.plate.iter().flat_map(|p| p.0.to_le_bytes()).collect();
    std::fs::write(out_dir.join("plate_in.u32"), plate_u32).unwrap();
    std::fs::write(
        out_dir.join("kind_in.u8"),
        state.kind.iter().map(|k| k.code()).collect::<Vec<u8>>(),
    )
    .unwrap();
    std::fs::write(
        out_dir.join("claim_mask.u32"),
        claim_mask
            .iter()
            .flat_map(|m| m.to_le_bytes())
            .collect::<Vec<u8>>(),
    )
    .unwrap();
    std::fs::write(
        out_dir.join("winner.u32"),
        winner_of
            .iter()
            .flat_map(|m| m.to_le_bytes())
            .collect::<Vec<u8>>(),
    )
    .unwrap();
    let pj: Vec<String> = poles
        .iter()
        .map(|(id, e)| {
            format!(
                "{{\"plate\":{},\"lat\":{},\"lon\":{},\"rate_deg_per_my\":{}}}",
                id.0, e.lat_deg, e.lon_deg, e.rate_deg_per_my
            )
        })
        .collect();
    let summary = format!(
        "{{\"seed\":{seed},\"level\":{level},\"mode\":\"{mode}\",\"nlat\":{},\"nlon\":{},\"dt_my\":{dt},\"t_my\":{},\"refusal\":{},\"overlap_cells\":{overlaps},\"cell_angle_rad\":{cell_angle},\"poles\":[{}],\"bad\":[{}]}}",
        grid.nlat(),
        grid.nlon(),
        t_fail.unwrap_or(state.t_my),
        refusal.as_ref().map(|r| format!("{r:?}")).unwrap_or("null".into()),
        pj.join(","),
        bad.join(",\n")
    );
    std::fs::write(out_dir.join("diag.json"), summary).unwrap();
    println!(
        "seed {seed} L{level} {mode}: t={:?} refusal={:?} bad={}",
        t_fail,
        refusal,
        bad.len()
    );
}
