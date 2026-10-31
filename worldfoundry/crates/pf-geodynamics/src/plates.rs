//! Seeded plate generation (brief K1/K2; acceptance: same seed -> identical
//! plates; plate coverage = entire sphere; no parcel code here).
//!
//! Determinism contract: the ONLY randomness is `SplitMix64` (inlined below,
//! known test vectors locked in tests) seeded from the artifact's `SeedPath`
//! via a documented byte mixing. Draw order is fixed: seed points first
//! (rejection loop), then one Euler pole per plate in plate-id order. Same
//! `SeedPath` + same params -> bit-identical output. `HashMap` iteration
//! order never influences results (BTreeSet/BTreeMap only).
//!
//! Crust TYPING (continental seeds) is exposed as plate flags; the material
//! ledger that consumes it is parcels.rs (agent-B scope, PARCEL != CELL).

use crate::rotation::{great_circle_rad, lat_lon_to_unit};
use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::params::{PlanetParameterProfile, SeedPath};
use pf_grid::GridSpec;
use std::collections::BTreeSet;

/// SplitMix64 (Steele & Marsaglia): tiny, deterministic, no dependencies.
/// Locked by known-vector tests below; never swap the algorithm without a
/// determinism-break note (artifact digests depend on its outputs).
#[derive(Debug, Clone, PartialEq)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform f64 in [0, 1): 53-bit mantissa fill.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

/// Mix a SeedPath into one 64-bit stream seed: FNV-1a over the seed's u64
/// bytes, derivation and salt, finalized through one SplitMix avalanche.
/// Documented and stable: artifacts hash downstream of these draws.
pub fn seed_from_path(path: &SeedPath) -> u64 {
    seed_from_path_tagged(path, "")
}

/// Tagged variant: independent deterministic streams (plates, hotspots,
/// later consumers) from one SeedPath — same path + tag, same stream;
/// different tag, independent stream.
pub fn seed_from_path_tagged(path: &SeedPath, tag: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in path.seed.to_le_bytes() {
        h ^= byte as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    for byte in path
        .derivation
        .as_bytes()
        .iter()
        .chain(path.salt.as_bytes())
        .chain(tag.as_bytes())
    {
        h ^= *byte as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    SplitMix64::new(h).next_u64()
}

/// Loud plate-seeding errors; each names its field/context.
#[derive(Debug, Clone, PartialEq)]
pub enum PlatesError {
    Params(&'static str, String),
    SpacingExhausted {
        attempts: u64,
        plate_count: u32,
        min_spacing_deg: f64,
    },
    PlateUnused {
        plate: PlateId,
    },
    Rotation(crate::rotation::RotationError),
}

impl std::fmt::Display for PlatesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PlatesError::Params(stage, why) => write!(f, "parameter validation at {stage}: {why}"),
            PlatesError::SpacingExhausted { attempts, plate_count, min_spacing_deg } => write!(
                f,
                "could not place {plate_count} seeds with min spacing {min_spacing_deg:?} deg after {attempts} attempts"
            ),
            PlatesError::PlateUnused { plate } => write!(f, "plate {plate:?} owns no cells (coverage broken)"),
            PlatesError::Rotation(e) => write!(f, "rotation: {e}"),
        }
    }
}

impl std::error::Error for PlatesError {}

impl From<crate::rotation::RotationError> for PlatesError {
    fn from(e: crate::rotation::RotationError) -> Self {
        PlatesError::Rotation(e)
    }
}

/// Deterministic plate seeding result.
///
/// `assignment` is indexed by flat element id (`grid.flat(i, j)`), C-order,
/// one `PlateId` per cell: full-sphere coverage, exactly one plate per cell
/// (Voronoi by great-circle distance; ties break to the LOWEST plate id —
/// declared, deterministic).
#[derive(Debug, Clone, PartialEq)]
pub struct PlateSeeding {
    pub plates: Vec<PlateRecordV0>,
    /// Per-element owning plate (len == grid.cells()).
    pub assignment: Vec<PlateId>,
    /// Plate ids whose initial crust is continental (K4 seed typing); the
    /// parcel ledger consumes this, never re-derives it.
    pub continental_plate_ids: BTreeSet<PlateId>,
    /// Voronoi seed points (unit vectors), plate-id order.
    pub seed_points: Vec<[f64; 3]>,
}

/// Minimum angular separation demanded between seed points, in degrees:
/// a quarter of the mean spacing a uniform distribution of `n` points
/// would have. Deliberately generous (rejection stays cheap); declared so
/// tests can assert it held.
pub fn min_seed_spacing_deg(plate_count: u32) -> f64 {
    // Mean angular spacing ~ sqrt(4*pi / n) rad; quarter of it, in deg.
    let mean_rad = (4.0 * std::f64::consts::PI / plate_count as f64).sqrt();
    mean_rad * 0.25 * 180.0 / std::f64::consts::PI
}

/// Rejection cap: hard limit, loud failure (never an infinite loop, never
/// a silent spacing relaxation).
pub const MAX_SEED_ATTEMPTS: u64 = 200_000;

/// Seed plates deterministically (K1): seed points by min-spacing rejection
/// sampling, initial Voronoi over grid centroids, per-plate Euler poles
/// drawn in plate order (direction uniform on the sphere, signed rate
/// within +/- rate_scale_deg_per_my; a near-zero rate is legal — stationary
/// plates exist).
pub fn seed_plates(
    grid: &GridSpec,
    params: &PlanetParameterProfile,
    seed_path: &SeedPath,
) -> Result<PlateSeeding, PlatesError> {
    params
        .validate()
        .map_err(|e| PlatesError::Params("PlanetParameterProfile.validate", e.to_string()))?;
    let t = &params.tectonics;
    if t.continental_seeds > t.plate_count {
        return Err(PlatesError::Params(
            "tectonics.continental_seeds",
            format!(
                "{} continental seeds exceed {} plates",
                t.continental_seeds, t.plate_count
            ),
        ));
    }

    let mut rng = SplitMix64::new(seed_from_path(seed_path));
    let min_spacing_rad = min_seed_spacing_deg(t.plate_count) * std::f64::consts::PI / 180.0;

    // --- seed points: uniform sphere + min-spacing rejection ---
    let mut seeds: Vec<[f64; 3]> = Vec::with_capacity(t.plate_count as usize);
    let mut attempts: u64 = 0;
    while seeds.len() < t.plate_count as usize {
        attempts += 1;
        if attempts > MAX_SEED_ATTEMPTS {
            return Err(PlatesError::SpacingExhausted {
                attempts,
                plate_count: t.plate_count,
                min_spacing_deg: min_spacing_rad * 180.0 / std::f64::consts::PI,
            });
        }
        // Uniform point on the sphere: z = 2u-1, phi = 2*pi*v (archimedean).
        let z = 2.0 * rng.next_f64() - 1.0;
        let phi = 2.0 * std::f64::consts::PI * rng.next_f64();
        let r = (1.0 - z * z).max(0.0).sqrt();
        let p = [r * phi.cos(), r * phi.sin(), z];
        if p.iter().any(|v| !v.is_finite()) {
            return Err(PlatesError::Params(
                "rng",
                format!("nonfinite sample at attempt {attempts}"),
            ));
        }
        let ok = seeds
            .iter()
            .all(|s| great_circle_rad(*s, p).unwrap_or(f64::INFINITY) >= min_spacing_rad);
        if ok {
            seeds.push(p);
        }
    }

    // --- Euler poles per plate, plate-id order (fixed draw order) ---
    let mut plates = Vec::with_capacity(seeds.len());
    for (idx, seed) in seeds.iter().enumerate() {
        let z = 2.0 * rng.next_f64() - 1.0;
        let phi = 2.0 * std::f64::consts::PI * rng.next_f64();
        let r = (1.0 - z * z).max(0.0).sqrt();
        let pole = [r * phi.cos(), r * phi.sin(), z];
        let (plat, plon) = crate::rotation::unit_to_lat_lon(pole)?;
        let rate = t.rate_scale_deg_per_my * (2.0 * rng.next_f64() - 1.0);
        plates.push(PlateRecordV0 {
            plate_id: PlateId(idx as u32),
            parent_plate: None,
            birth_t_my: 0.0,
            euler_pole: Some(EulerPoleV0 {
                lat_deg: plat,
                lon_deg: plon,
                rate_deg_per_my: rate,
            }),
            retirement: None,
        });
        let _ = seed; // seeds themselves are not poles; keep loop shape explicit
    }

    // --- Voronoi assignment over grid centroids ---
    let seed_units: Vec<[f64; 3]> = seeds.clone();
    let mut assignment: Vec<PlateId> = Vec::with_capacity(grid.cells());
    for i in 0..grid.nlat() {
        let lat = grid
            .center_lat_deg(i)
            .map_err(|e| PlatesError::Params("grid.center_lat_deg", e.to_string()))?;
        for j in 0..grid.nlon() {
            let lon = grid
                .center_lon_deg(j)
                .map_err(|e| PlatesError::Params("grid.center_lon_deg", e.to_string()))?;
            let c = lat_lon_to_unit(lat, lon)?;
            // Nearest seed; ties break to the lowest plate id (declared):
            // strict `<` keeps the first (lowest-id) seed on exact ties.
            let mut best_d = f64::INFINITY;
            let mut best_k = 0usize;
            for (k, s) in seed_units.iter().enumerate() {
                let d = great_circle_rad(c, *s)?;
                if d < best_d {
                    best_d = d;
                    best_k = k;
                }
            }
            assignment.push(PlateId(best_k as u32));
        }
    }
    if assignment.len() != grid.cells() {
        return Err(PlatesError::Params(
            "assignment",
            format!(
                "assignment has {} entries for {} cells",
                assignment.len(),
                grid.cells()
            ),
        ));
    }

    // Coverage gate: every plate owns at least one cell (Voronoi of distinct
    // seeds guarantees this mathematically; the check is the loud backstop).
    let mut used = vec![false; plates.len()];
    for pid in &assignment {
        used[pid.0 as usize] = true;
    }
    for (k, u) in used.iter().enumerate() {
        if !u {
            return Err(PlatesError::PlateUnused {
                plate: PlateId(k as u32),
            });
        }
    }

    let continental_plate_ids: BTreeSet<PlateId> = (0..t.continental_seeds).map(PlateId).collect();

    Ok(PlateSeeding {
        plates,
        assignment,
        continental_plate_ids,
        seed_points: seeds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_contracts::params::{FidelityPreset, MantleBlock, ScenarioPreset};

    fn seed_path(tag: &str) -> SeedPath {
        SeedPath {
            seed: 240907,
            derivation: "toy".into(),
            salt: tag.into(),
        }
    }

    fn setup(tag: &str) -> (GridSpec, PlanetParameterProfile) {
        let grid = GridSpec::new(10, 18, 6.371e6).unwrap();
        let mut p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
        p.tectonics.plate_count = 5;
        p.tectonics.continental_seeds = 2;
        let _ = tag;
        (grid, p)
    }

    #[test]
    fn splitmix64_known_vectors() {
        // Goldens computed from the reference algorithm (verified against an
        // independent Python implementation 2026-09-07). Purpose: drift
        // detection — any change to these constants or the mixing breaks
        // determinism of every artifact seeded through this path.
        let mut r = SplitMix64::new(0);
        assert_eq!(
            r.next_u64(),
            0xE220_A839_7B1D_CDAF,
            "SplitMix64 vector seed=0"
        );
        let mut r = SplitMix64::new(0x1234_5678_9ABC_DEF0);
        let first = r.next_u64();
        assert_eq!(
            first, 0x1619_22C6_45CE_50E8,
            "SplitMix64 vector seed=0x1234.."
        );
        // Determinism: fresh instance replays identically.
        let mut r2 = SplitMix64::new(0x1234_5678_9ABC_DEF0);
        assert_eq!(r2.next_u64(), first);
        // f64 fills stay in [0, 1).
        let mut r3 = SplitMix64::new(42);
        for _ in 0..1000 {
            let u = r3.next_f64();
            assert!((0.0..1.0).contains(&u), "next_f64 out of [0,1): {u}");
        }
    }

    #[test]
    fn same_seed_path_is_bit_identical() {
        let (grid, params) = setup("det");
        let a = seed_plates(&grid, &params, &seed_path("a")).unwrap();
        let b = seed_plates(&grid, &params, &seed_path("a")).unwrap();
        assert_eq!(a, b, "same SeedPath must produce bit-identical seeding");
    }

    #[test]
    fn different_seed_paths_differ() {
        let (grid, params) = setup("diff");
        let a = seed_plates(&grid, &params, &seed_path("alpha")).unwrap();
        let b = seed_plates(&grid, &params, &seed_path("beta")).unwrap();
        assert_ne!(
            a.assignment, b.assignment,
            "different salt must change the world"
        );
    }

    #[test]
    fn full_sphere_coverage_and_counts() {
        let (grid, params) = setup("cov");
        let s = seed_plates(&grid, &params, &seed_path("c")).unwrap();
        assert_eq!(s.assignment.len(), grid.cells());
        assert_eq!(s.plates.len(), params.tectonics.plate_count as usize);
        assert_eq!(s.seed_points.len(), s.plates.len());
        let mut seen = BTreeSet::new();
        for pid in &s.assignment {
            seen.insert(*pid);
        }
        assert_eq!(
            seen.len(),
            params.tectonics.plate_count as usize,
            "every plate in use"
        );
        assert_eq!(
            s.continental_plate_ids.len(),
            params.tectonics.continental_seeds as usize
        );
        // Euler poles present and valid for every plate.
        for p in &s.plates {
            p.euler_pole
                .expect("live plate carries a pole")
                .validate()
                .unwrap();
        }
    }

    #[test]
    fn min_spacing_holds_between_all_seeds() {
        let (grid, params) = setup("sp");
        let s = seed_plates(&grid, &params, &seed_path("d")).unwrap();
        let min_rad =
            min_seed_spacing_deg(params.tectonics.plate_count) * std::f64::consts::PI / 180.0;
        for i in 0..s.seed_points.len() {
            for j in (i + 1)..s.seed_points.len() {
                let d = great_circle_rad(s.seed_points[i], s.seed_points[j]).unwrap();
                assert!(
                    d >= min_rad * (1.0 - 1.0e-12),
                    "seeds {i}/{j} closer than declared spacing"
                );
            }
        }
    }

    #[test]
    fn assignment_indexing_matches_flat_element_ids() {
        // assignment[grid.flat(i,j)] must be the plate chosen for (i,j):
        // rebuild a single cell's owner directly and compare.
        let (grid, params) = setup("idx");
        let s = seed_plates(&grid, &params, &seed_path("e")).unwrap();
        let (i, j) = (4usize, 17usize); // antimeridian-adjacent column
        let element = grid.flat(i, j) as usize;
        let lat = grid.center_lat_deg(i).unwrap();
        let lon = grid.center_lon_deg(j).unwrap();
        let c = lat_lon_to_unit(lat, lon).unwrap();
        let mut best = (f64::INFINITY, 0usize);
        for (k, sp) in s.seed_points.iter().enumerate() {
            let d = great_circle_rad(c, *sp).unwrap();
            if d < best.0 {
                best = (d, k);
            }
        }
        assert_eq!(s.assignment[element], PlateId(best.1 as u32));
    }

    #[test]
    fn continental_seeds_above_plate_count_fail_loudly() {
        let (grid, mut params) = setup("bad");
        params.tectonics.continental_seeds = params.tectonics.plate_count + 1;
        let err = seed_plates(&grid, &params, &seed_path("f")).unwrap_err();
        assert!(
            err.to_string().contains("continental seeds exceed"),
            "{err}"
        );
    }

    #[test]
    fn rigid_motion_preserves_intra_plate_neighborhoods() {
        // Rotate every centroid by ITS plate's Euler pole: each plate's own
        // set of points moves rigidly, so intra-plate distances hold.
        let (grid, params) = setup("rigid");
        let s = seed_plates(&grid, &params, &seed_path("g")).unwrap();
        let dt = 7.5;
        let mut moved: Vec<[f64; 3]> = Vec::with_capacity(grid.cells());
        for i in 0..grid.nlat() {
            let lat = grid.center_lat_deg(i).unwrap();
            for j in 0..grid.nlon() {
                let lon = grid.center_lon_deg(j).unwrap();
                let p = lat_lon_to_unit(lat, lon).unwrap();
                let pid = s.assignment[grid.flat(i, j) as usize];
                let pole = s.plates[pid.0 as usize].euler_pole.unwrap();
                moved.push(crate::rotation::rotate_by_euler(p, &pole, dt).unwrap());
            }
        }
        // Sample pairs within one plate and check distance preservation.
        let mut checked = 0;
        'outer: for i in 0..grid.cells() {
            let pi = s.assignment[i];
            for j in (i + 1)..grid.cells() {
                if s.assignment[j] != pi {
                    continue;
                }
                let a0 = crate::plates::tests_support::centroid_unit(&grid, i);
                let b0 = crate::plates::tests_support::centroid_unit(&grid, j);
                let d0 = great_circle_rad(a0, b0).unwrap();
                let d1 = great_circle_rad(moved[i], moved[j]).unwrap();
                assert!(
                    (d0 - d1).abs() < 1.0e-9,
                    "intra-plate distance drifted for plate {pi:?}: {d0} -> {d1}"
                );
                checked += 1;
                if checked >= 500 {
                    break 'outer;
                }
            }
        }
        assert!(
            checked >= 50,
            "test sampled only {checked} intra-plate pairs"
        );
    }
}

/// Test-only helpers shared by module test suites.
#[cfg(test)]
pub(crate) mod tests_support {
    use pf_grid::GridSpec;

    pub fn centroid_unit(grid: &GridSpec, element: usize) -> [f64; 3] {
        let (i, j) = grid.ij(element as u64).expect("element in range");
        let lat = grid.center_lat_deg(i).unwrap();
        let lon = grid.center_lon_deg(j).unwrap();
        crate::rotation::lat_lon_to_unit(lat, lon).unwrap()
    }
}
