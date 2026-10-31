//! Known-answer tests for the OPT-IN continental crust bodies
//! (CRUST-SPIKE, lane K2, 2026-10-11): `CrustState::seed_bodies` and the
//! `linear-rifted-margin-taper-v0` thickness profile.
//!
//! Expected values are recomputed here from first principles (cell areas,
//! great-circle distances, the taper formula written out), never read back
//! from the producer's own bookkeeping.

use pf_contracts::params::{FidelityPreset, MantleBlock, ScenarioPreset, SeedPath};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::kinematic::CrustState;
use pf_geodynamics::plates::{seed_plates, PlateSeeding};
use pf_geodynamics::rotation::{great_circle_rad, lat_lon_to_unit};
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;

const H_CC: f64 = 35_000.0; // uplift.rs CONTINENTAL_CRUST_THICKNESS_M
const H_OC: f64 = 7_000.0; // uplift.rs OCEANIC_CRUST_THICKNESS_M

fn setup(nlat: usize, nlon: usize, seed: u64) -> (GridSpec, PlateSeeding) {
    let grid = GridSpec::new(nlat, nlon, 6.371e6).unwrap();
    let p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
    let sp = SeedPath {
        seed,
        derivation: "kinematic-preview".into(),
        salt: "toy".into(),
    };
    let seeding = seed_plates(&grid, &p, &sp).unwrap();
    (grid, seeding)
}

fn areas(grid: &GridSpec) -> Vec<f64> {
    let mut a = Vec::new();
    for i in 0..grid.nlat() {
        let x = grid.cell_area_m2(i).unwrap();
        for _ in 0..grid.nlon() {
            a.push(x);
        }
    }
    a
}

fn centroid(grid: &GridSpec, c: usize) -> [f64; 3] {
    let (i, j) = grid.ij(c as u64).unwrap();
    lat_lon_to_unit(grid.center_lat_deg(i).unwrap(), grid.center_lon_deg(j).unwrap()).unwrap()
}

fn taper(d: f64, w: f64) -> f64 {
    if w <= 0.0 || d >= w {
        H_CC
    } else {
        H_OC + (H_CC - H_OC) * d / w
    }
}

#[test]
fn coverage_matches_declared_fraction_within_half_a_cell_per_body() {
    for (seed, frac) in [(240907u64, 0.38), (7, 0.30), (42, 0.35)] {
        let (grid, seeding) = setup(72, 144, seed);
        let a = areas(&grid);
        let total: f64 = a.iter().sum();
        let (state, rep) = CrustState::seed_bodies(
            &grid,
            &seeding,
            &UpliftPrior::preview_v0(),
            frac,
            200_000.0,
        )
        .unwrap();
        // Independent area of continental kind.
        let cont: f64 = (0..grid.cells())
            .filter(|&c| state.kind[c] == CrustKindCodeV0::Continental)
            .map(|c| a[c])
            .sum();
        let a_max = a.iter().cloned().fold(0.0, f64::max);
        let k = seeding.continental_plate_ids.len() as f64;
        assert!(
            (cont - frac * total).abs() <= 0.5 * a_max * k,
            "seed {seed}: continental area {} vs declared {}",
            cont / total,
            frac
        );
        assert!((rep.achieved_area_fraction - cont / total).abs() < 1e-12);
        assert_eq!(rep.declared_area_fraction, frac);
        // Per-body targets proportional to host-plate area, each met to
        // half of one cell.
        let cont_plate_area: f64 = (0..grid.cells())
            .filter(|&c| seeding.continental_plate_ids.contains(&seeding.assignment[c]))
            .map(|c| a[c])
            .sum();
        for b in &rep.bodies {
            let plate_area: f64 = (0..grid.cells())
                .filter(|&c| seeding.assignment[c] == b.plate)
                .map(|c| a[c])
                .sum();
            let target = frac * total * plate_area / cont_plate_area;
            assert!((b.target_area_m2 - target).abs() <= 1e-6 * target);
            let body_area: f64 = (0..grid.cells())
                .filter(|&c| {
                    seeding.assignment[c] == b.plate
                        && state.kind[c] == CrustKindCodeV0::Continental
                })
                .map(|c| a[c])
                .sum();
            // Same cells, different summation order: relative 1e-12.
            assert!((body_area - b.achieved_area_m2).abs() <= 1e-12 * body_area);
            assert!((body_area - target).abs() <= 0.5 * a_max, "plate {:?}", b.plate);
        }
    }
}

#[test]
fn bodies_live_inside_continental_plates_as_geodesic_discs() {
    let (grid, seeding) = setup(72, 144, 240907);
    let (state, rep) =
        CrustState::seed_bodies(&grid, &seeding, &UpliftPrior::preview_v0(), 0.38, 200_000.0)
            .unwrap();
    for c in 0..grid.cells() {
        if state.kind[c] == CrustKindCodeV0::Continental {
            assert!(seeding.continental_plate_ids.contains(&seeding.assignment[c]));
            assert_eq!(state.birth_t_my[c], -2000.0);
        } else {
            assert_eq!(state.thickness_m[c], H_OC, "non-body cells are 7 km oceanic");
            assert_eq!(state.birth_t_my[c], -50.0);
        }
        assert_eq!(state.plate[c], seeding.assignment[c]);
    }
    // Disc property: within a host plate, every body cell is at least as
    // close to the seed point as every non-body cell (ties allowed).
    for b in &rep.bodies {
        let mut max_in: f64 = 0.0;
        let mut min_out = f64::INFINITY;
        for c in 0..grid.cells() {
            if seeding.assignment[c] != b.plate {
                continue;
            }
            let d = great_circle_rad(centroid(&grid, c), b.seed_point).unwrap();
            if state.kind[c] == CrustKindCodeV0::Continental {
                max_in = max_in.max(d);
            } else {
                min_out = min_out.min(d);
            }
        }
        assert!(max_in <= min_out + 1e-12, "{max_in} > {min_out}");
    }
    // Some continental plate keeps ocean (bodies are not whole plates).
    let ocean_on_cont_plate = (0..grid.cells()).any(|c| {
        seeding.continental_plate_ids.contains(&seeding.assignment[c])
            && state.kind[c] == CrustKindCodeV0::Oceanic
    });
    assert!(ocean_on_cont_plate);
}

#[test]
fn thickness_follows_the_taper_and_is_monotone_toward_the_margin() {
    let w = 200_000.0;
    let (grid, seeding) = setup(72, 144, 240907);
    let (state, rep) =
        CrustState::seed_bodies(&grid, &seeding, &UpliftPrior::preview_v0(), 0.38, w).unwrap();
    let n = grid.cells();
    let cont = |c: usize| state.kind[c] == CrustKindCodeV0::Continental;
    // Independent coast distance: brute force over this test's own coast
    // midpoints, with great_circle_rad (not the producer's dot product).
    let mut mids = Vec::new();
    for c in 0..n {
        if !cont(c) {
            continue;
        }
        let (i, j) = grid.ij(c as u64).unwrap();
        for (ti, tj, _) in grid.neighbors(i, j).unwrap() {
            let nb = ti * grid.nlon() + tj;
            if !cont(nb) {
                let (p, q) = (centroid(&grid, c), centroid(&grid, nb));
                let s = [p[0] + q[0], p[1] + q[1], p[2] + q[2]];
                let l = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
                mids.push([s[0] / l, s[1] / l, s[2] / l]);
            }
        }
    }
    assert_eq!(mids.len() as u64, rep.coast_midpoints);
    let mut pairs = Vec::new();
    for c in (0..n).filter(|&c| cont(c)) {
        let p = centroid(&grid, c);
        let d = mids
            .iter()
            .map(|m| great_circle_rad(p, *m).unwrap())
            .fold(f64::INFINITY, f64::min)
            * grid.radius_m();
        assert!(
            (d - rep.coast_distance_m[c]).abs() < 1.0,
            "cell {c}: distance {d} vs {}",
            rep.coast_distance_m[c]
        );
        let h = taper(d, w);
        assert!(
            (state.thickness_m[c] - h).abs() < 0.2,
            "cell {c}: thickness {} vs hand taper {h}",
            state.thickness_m[c]
        );
        pairs.push((d, state.thickness_m[c]));
    }
    // Monotone: sorting body cells by coast distance gives non-decreasing
    // thickness (thinnest at the margin, interior thickness inland).
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    for w2 in pairs.windows(2) {
        assert!(w2[1].1 >= w2[0].1 - 1e-6, "{w2:?}");
    }
    assert!(pairs.first().unwrap().1 < H_CC, "margin cells are thinned");
    assert_eq!(pairs.last().unwrap().1, H_CC, "interior reaches 35 km");
    assert!(pairs.iter().all(|p| p.1 >= H_OC));
    // A coast cell is within half a neighbour spacing of the coast.
    for c in (0..n).filter(|&c| cont(c)) {
        let (i, j) = grid.ij(c as u64).unwrap();
        let mut nearest_out = f64::INFINITY;
        for (ti, tj, nbinfo) in grid.neighbors(i, j).unwrap() {
            if !cont(ti * grid.nlon() + tj) {
                nearest_out = nearest_out.min(nbinfo.distance_m);
            }
        }
        if nearest_out.is_finite() {
            assert!(rep.coast_distance_m[c] <= 0.5 * nearest_out * 1.0001 + 1.0);
        }
    }
}

#[test]
fn volume_bookkeeping_closes_and_plateau_limit_is_exact() {
    let (grid, seeding) = setup(72, 144, 7);
    let a = areas(&grid);
    let prior = UpliftPrior::preview_v0();
    let (state, rep) = CrustState::seed_bodies(&grid, &seeding, &prior, 0.38, 200_000.0).unwrap();
    let cont_volume: f64 = (0..grid.cells())
        .filter(|&c| state.kind[c] == CrustKindCodeV0::Continental)
        .map(|c| state.thickness_m[c] * a[c])
        .sum();
    let sum_bodies: f64 = rep.bodies.iter().map(|b| b.crust_volume_m3).sum();
    assert!((cont_volume - rep.continental_crust_volume_m3).abs() <= 1e-9 * cont_volume);
    assert!((sum_bodies - cont_volume).abs() <= 1e-9 * cont_volume);
    let cont_area: f64 = (0..grid.cells())
        .filter(|&c| state.kind[c] == CrustKindCodeV0::Continental)
        .map(|c| a[c])
        .sum();
    assert!((rep.plateau_equivalent_volume_m3 - cont_area * H_CC).abs() <= 1e-9 * cont_area * H_CC);
    assert!(rep.continental_crust_volume_m3 < rep.plateau_equivalent_volume_m3);
    // Total crust volume = continental + 7 km over the rest.
    let ocean_area: f64 = a.iter().sum::<f64>() - cont_area;
    let total = state.crust_volume_m3(&a);
    assert!((total - (cont_volume + H_OC * ocean_area)).abs() <= 1e-9 * total);
    // W = 0: plateau bodies of exactly 35 km, same footprint.
    let (plateau, rep0) = CrustState::seed_bodies(&grid, &seeding, &prior, 0.38, 0.0).unwrap();
    assert_eq!(plateau.kind, state.kind, "the width never moves the footprint");
    for c in 0..grid.cells() {
        if plateau.kind[c] == CrustKindCodeV0::Continental {
            assert_eq!(plateau.thickness_m[c], H_CC);
        }
    }
    assert!((rep0.continental_crust_volume_m3 - cont_area * H_CC).abs() <= 1e-9 * cont_area * H_CC);
    assert!(rep0.bodies.iter().all(|b| b.tapered_cells == 0));
}

#[test]
fn seeding_is_deterministic_and_leaves_the_plate_stream_untouched() {
    let (grid, seeding) = setup(36, 72, 1234567);
    let before = seeding.clone();
    let prior = UpliftPrior::preview_v0();
    let x = CrustState::seed_bodies(&grid, &seeding, &prior, 0.38, 200_000.0).unwrap();
    let y = CrustState::seed_bodies(&grid, &seeding, &prior, 0.38, 200_000.0).unwrap();
    assert_eq!(x.0, y.0);
    assert_eq!(x.1.bodies, y.1.bodies);
    assert_eq!(seeding, before, "the plate seeding is read, never modified");
    // Re-seeding plates from the same path is unchanged by body seeding.
    let (_, again) = setup(36, 72, 1234567);
    assert_eq!(again, before);
}

#[test]
fn coverage_beyond_the_continental_plates_is_refused_loudly() {
    let (grid, seeding) = setup(36, 72, 240907);
    let err = CrustState::seed_bodies(&grid, &seeding, &UpliftPrior::preview_v0(), 0.95, 200_000.0)
        .unwrap_err()
        .to_string();
    assert!(err.contains("refused"), "{err}");
    assert!(CrustState::seed_bodies(&grid, &seeding, &UpliftPrior::preview_v0(), 1.0, 1.0).is_err());
    assert!(CrustState::seed_bodies(&grid, &seeding, &UpliftPrior::preview_v0(), 0.3, -1.0).is_err());
}
