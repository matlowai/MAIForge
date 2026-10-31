//! SEED-IRREGULAR-SEEDS (lane IS, 2026-10-13): `irregular-cratons-v1`
//! seeding invariants on a 72 x 144 grid with the toy plate seeding:
//! determinism from the keyed stream, the declared continental area is
//! honoured, every assembly rides one plate, the margin profile and its
//! stretch bookkeeping are consistent, no small enclosed ocean remains.

use pf_contracts::params::{FidelityPreset, MantleBlock, ScenarioPreset, SeedPath};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::cratons::{seed_cratons, CratonParams, ENCLOSED_OCEAN_FILL_FRACTION};
use pf_geodynamics::plates::{seed_from_path_tagged, seed_plates};
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;

fn params() -> CratonParams {
    CratonParams {
        interior_thickness_m: 29_600.0,
        rank_size_exponent: 1.4,
        cratons_per_unit_share: 8.0,
        basement_roughness_std_m: 1_500.0,
        suture_crest_thickening_m: 7_000.0,
        suture_half_width_m: 120_000.0,
        proximal_margin_stretch: 1.15,
    }
}

fn run(seed: u64) -> (GridSpec, pf_geodynamics::kinematic::CrustState, pf_geodynamics::kinematic::ContinentalBodiesReport, pf_geodynamics::cratons::CratonSeedReport) {
    let grid = GridSpec::new(72, 144, 6.371e6).unwrap();
    let p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
    let sp = SeedPath { seed, derivation: "toy".into(), salt: "is-test".into() };
    let seeding = seed_plates(&grid, &p, &sp).unwrap();
    let prior = UpliftPrior::preview_v0();
    let (s, r, x) = seed_cratons(
        &grid,
        &seeding,
        &prior,
        0.38,
        200_000.0,
        &params(),
        p.tectonics.continental_seeds,
        seed_from_path_tagged(&sp, "irregular-cratons-v1"),
    )
    .unwrap();
    (grid, s, r, x)
}

#[test]
fn deterministic_and_seed_dependent() {
    let (_, a, ra, xa) = run(42);
    let (_, b, rb, xb) = run(42);
    assert_eq!(a, b);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
    assert_eq!(bits(&ra.coast_distance_m), bits(&rb.coast_distance_m));
    assert_eq!(ra.bodies, rb.bodies);
    assert_eq!(bits(&xa.interior_thickness_m), bits(&xb.interior_thickness_m));
    assert_eq!(bits(&xa.stretch_beta), bits(&xb.stretch_beta));
    assert_eq!(xa.cratons, xb.cratons);
    let (_, c, _, _) = run(240907);
    assert_ne!(a.thickness_m, c.thickness_m);
}

#[test]
fn area_plates_margins_and_holes() {
    for seed in [42u64, 240907] {
        let (grid, s, r, x) = run(seed);
        let n = grid.cells();
        let areas: Vec<f64> = (0..grid.nlat())
            .flat_map(|i| std::iter::repeat(grid.cell_area_m2(i).unwrap()).take(grid.nlon()))
            .collect();
        let total: f64 = areas.iter().sum();
        let cont: f64 = (0..n).filter(|&c| s.kind[c] == CrustKindCodeV0::Continental).map(|c| areas[c]).sum();
        // Declared area within one equatorial cell plus the filled holes.
        let tol = areas[(grid.nlat() / 2) * grid.nlon()] + x.filled_area_m2;
        assert!((cont - 0.38 * total).abs() <= tol, "seed {seed}: {} vs 0.38", cont / total);
        assert!((r.achieved_area_fraction - cont / total).abs() < 1e-12);
        // Each assembly's cells share its plate.
        for b in &r.bodies {
            assert!(b.cells > 0, "seed {seed}: empty assembly");
        }
        assert!(x.cratons.len() >= r.bodies.len());
        for c in 0..n {
            if s.kind[c] == CrustKindCodeV0::Continental {
                let h_int = x.interior_thickness_m[c];
                assert!(h_int >= 15_000.0 && s.thickness_m[c] <= h_int + 1e-9, "seed {seed} cell {c}");
                assert!(s.thickness_m[c] >= 7_000.0 - 1e-9);
                let beta = x.stretch_beta[c];
                assert!(beta >= 1.0);
                if beta > 1.0 {
                    assert!((beta * s.thickness_m[c] - h_int).abs() < 1e-6 * h_int);
                    assert!(r.coast_distance_m[c] < 200_000.0);
                }
            } else {
                assert!(x.interior_thickness_m[c].is_nan() && x.stretch_beta[c] == 1.0);
                assert_eq!(s.thickness_m[c], 7_000.0);
            }
        }
        let plates: std::collections::BTreeSet<u32> =
            r.bodies.iter().map(|b| b.plate.0).collect();
        for c in 0..n {
            if s.kind[c] == CrustKindCodeV0::Continental {
                assert!(plates.contains(&s.plate[c].0), "seed {seed}: continental cell on a non-assembly plate");
            }
        }
        // Interior median near the declared base (roughness is zero-mean).
        let mut hi: Vec<f64> = (0..n).filter(|&c| x.stretch_beta[c] == 1.0 && s.kind[c] == CrustKindCodeV0::Continental).map(|c| s.thickness_m[c]).collect();
        hi.sort_by(|a, b| a.total_cmp(b));
        let med = hi[hi.len() / 2];
        assert!((med - 29_600.0).abs() < 1_200.0, "seed {seed}: interior median {med}");
        assert!(x.suture_cells > 0, "seed {seed}: no suture");
        // No enclosed ocean smaller than the fill bound remains.
        let (nlat, nlon) = (grid.nlat(), grid.nlon());
        let mut lab = vec![usize::MAX; n];
        let mut comp = Vec::new();
        for st in 0..n {
            if s.kind[st] == CrustKindCodeV0::Continental || lab[st] != usize::MAX {
                continue;
            }
            let id = comp.len();
            let mut a = 0.0;
            let mut stack = vec![st];
            lab[st] = id;
            while let Some(c) = stack.pop() {
                a += areas[c];
                let (i, j) = (c / nlon, c % nlon);
                for di in [-1i64, 0, 1] {
                    for dj in [-1i64, 0, 1] {
                        let ti = i as i64 + di;
                        if ti < 0 || ti >= nlat as i64 || (di == 0 && dj == 0) {
                            continue;
                        }
                        let t = ti as usize * nlon + (j as i64 + dj).rem_euclid(nlon as i64) as usize;
                        if s.kind[t] != CrustKindCodeV0::Continental && lab[t] == usize::MAX {
                            lab[t] = id;
                            stack.push(t);
                        }
                    }
                }
            }
            comp.push(a);
        }
        comp.sort_by(|a, b| b.total_cmp(a));
        assert!(comp.iter().skip(1).all(|a| *a >= ENCLOSED_OCEAN_FILL_FRACTION * total), "seed {seed}: {comp:?}");
    }
}
