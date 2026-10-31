//! SEED-IRREGULAR-SEEDS + SEED-FREEBOARD start height (lane IS,
//! 2026-10-13): `irregular-cratons-v1`, an opt-in initial crust.
//!
//! The initial continental state is a plausible OLD one, roughly Earth's
//! when its continents had just been assembled, instead of one smooth
//! geodesic disc per continental plate:
//!
//! - continental ASSEMBLIES (count = `tectonics.continental_seeds`) with
//!   rank-size (power-law) area shares, each made of one or more CRATONS
//!   (elliptic caps of power-law sizes that overlap their neighbours);
//! - ragged coasts: the craton score is evaluated on a domain-warped
//!   position and a multi-scale edge noise is added before the global
//!   area ranking, which honours the declared continental area fraction;
//! - every continental cell rides the plate under its assembly centre, so
//!   coasts are not clipped by straight Voronoi edges (oceanic cells keep
//!   their Voronoi plate); old sutures therefore lie inside plates;
//! - interior thickness = base + multi-scale basement roughness (basins
//!   and arches) + a worn suture root along the seams between cratons;
//! - a two-segment rifted margin (necking/distal wedge, then a mildly
//!   stretched proximal shelf) in place of the single linear wedge.
//!
//! Determinism: one SplitMix64 stream from the tagged seed path, fixed
//! draw order (assembly shares, centres, cratons, then noise fields),
//! ascending cell order, ties to the lowest cell id; no HashMap.

use crate::kinematic::{
    cell_areas, centroids, ContinentalBodiesReport, ContinentalBodyRecord, CrustState,
    KinematicError,
};
use crate::plates::{PlateSeeding, SplitMix64};
use crate::uplift::UpliftPrior;
use pf_contracts::domain::PlateId;
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_grid::GridSpec;
use std::f64::consts::PI;

/// Law constants routed through the universe (`craton_seeding =
/// irregular-cratons-v1`); values and sources in the contract.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CratonParams {
    /// Interior (untapered, unroughened) continental thickness (m).
    pub interior_thickness_m: f64,
    /// Rank-size exponent a: share of rank k (1-based) ~ k^-a.
    pub rank_size_exponent: f64,
    /// Cratons in an assembly = 1 + floor(N_c * share).
    pub cratons_per_unit_share: f64,
    /// Standard deviation of the basement thickness roughness (m).
    pub basement_roughness_std_m: f64,
    /// Crest thickening of an old suture root (m).
    pub suture_crest_thickening_m: f64,
    /// Gaussian half-width (sigma) of a suture root (m).
    pub suture_half_width_m: f64,
    /// Stretch factor of the proximal margin at the necking hinge.
    pub proximal_margin_stretch: f64,
}

// Fixed by the law version (design and visual numbers, not physics;
// unpublished-design-record).
pub const WARP_AMPLITUDE_M: f64 = 450.0e3;
pub const WARP_WAVELENGTHS_M: [f64; 3] = [4.0e6, 2.0e6, 1.0e6];
pub const WARP_GAIN: f64 = 0.5;
pub const EDGE_NOISE_STD: f64 = 0.16;
pub const FBM_WAVELENGTHS_M: [f64; 5] = [3.0e6, 1.5e6, 7.5e5, 3.75e5, 1.9e5];
pub const FBM_GAIN: f64 = 0.6;
pub const STRIKE_STD: f64 = 0.35;
pub const STRIKE_WAVELENGTHS_M: [f64; 2] = [1.5e6, 7.5e5];
pub const STRIKE_GAIN: f64 = 0.5;
pub const ASPECT_MAX: f64 = 2.0;
pub const CRATON_OVERLAP: f64 = 0.8;
pub const ASSEMBLY_SPACING: f64 = 1.1;
pub const MAX_CENTRE_DRAWS: u32 = 20_000;
pub const ENCLOSED_OCEAN_FILL_FRACTION: f64 = 1.0e-3;
pub const SIZE_JITTER_SIGMA: f64 = 0.3;
pub const NECKING_FRACTION: f64 = 0.3;
pub const WAVES_PER_OCTAVE: usize = 24;
/// Floor on the interior thickness (m) so a deep basin never reaches the
/// oceanic column (loud bookkeeping: clamped cells are counted).
pub const MIN_INTERIOR_THICKNESS_M: f64 = 15_000.0;
/// Score band below the craton edge in which a second craton still counts
/// as reaching a cell (suture weight ramp).
pub const SUTURE_REACH_BAND: f64 = 0.3;

/// One craton.
#[derive(Debug, Clone, PartialEq)]
pub struct CratonRecord {
    pub assembly: usize,
    pub plate: PlateId,
    pub centre: [f64; 3],
    pub radius_rad: f64,
    pub aspect: f64,
    pub area_target_m2: f64,
    /// Mask cells whose best craton is this one.
    pub cells: u64,
}

/// Seeding bookkeeping beyond the shared bodies report.
#[derive(Debug, Clone, PartialEq)]
pub struct CratonSeedReport {
    pub cratons: Vec<CratonRecord>,
    /// Per assembly: plate, centre, area share.
    pub assemblies: Vec<(PlateId, [f64; 3], f64)>,
    /// Assembly centres taken as the best candidate after MAX_CENTRE_DRAWS.
    pub spacing_fallbacks: u32,
    pub filled_cells: u64,
    pub filled_area_m2: f64,
    /// Mask cells with suture weight above one half.
    pub suture_cells: u64,
    pub interior_clamped_cells: u64,
    /// Final global score threshold (score of the last cell taken).
    pub score_threshold: f64,
    /// Continental cells whose plate differs from their Voronoi plate.
    pub relabelled_cells: u64,
    /// Per cell: untapered interior thickness (m); NaN on oceanic cells.
    pub interior_thickness_m: Vec<f64>,
    /// Per cell: margin stretch h_int / h inside the taper, 1 elsewhere.
    pub stretch_beta: Vec<f64>,
}

/// Sum of random plane waves on unit-sphere coordinates (spectral noise).
struct Spectral {
    waves: Vec<([f64; 3], f64, f64)>,
}

fn unit_dir(rng: &mut SplitMix64) -> [f64; 3] {
    let z = 2.0 * rng.next_f64() - 1.0;
    let phi = 2.0 * PI * rng.next_f64();
    let r = (1.0 - z * z).max(0.0).sqrt();
    [r * phi.cos(), r * phi.sin(), z]
}

fn gaussian(rng: &mut SplitMix64) -> f64 {
    // Box-Muller, cosine branch; always two draws.
    let u1 = rng.next_f64();
    let u2 = rng.next_f64();
    (-2.0 * (1.0 - u1).ln()).sqrt() * (2.0 * PI * u2).cos()
}

impl Spectral {
    /// Octave amplitudes gain^o, scaled so the field's standard deviation
    /// is `std`; wavelengths along the surface of a sphere of radius_m.
    fn new(rng: &mut SplitMix64, radius_m: f64, wavelengths_m: &[f64], gain: f64, std: f64) -> Self {
        let amps: Vec<f64> = (0..wavelengths_m.len()).map(|o| gain.powi(o as i32)).collect();
        let norm = amps.iter().map(|a| a * a).sum::<f64>().sqrt();
        let per_wave = (2.0 / WAVES_PER_OCTAVE as f64).sqrt();
        let mut waves = Vec::with_capacity(wavelengths_m.len() * WAVES_PER_OCTAVE);
        for (o, &lam) in wavelengths_m.iter().enumerate() {
            let k = 2.0 * PI * radius_m / lam;
            let a = std * amps[o] / norm * per_wave;
            for _ in 0..WAVES_PER_OCTAVE {
                let d = unit_dir(rng);
                let phase = 2.0 * PI * rng.next_f64();
                waves.push(([d[0] * k, d[1] * k, d[2] * k], phase, a));
            }
        }
        Self { waves }
    }

    fn eval(&self, x: [f64; 3]) -> f64 {
        self.waves
            .iter()
            .map(|(k, ph, a)| a * (k[0] * x[0] + k[1] * x[1] + k[2] * x[2] + ph).cos())
            .sum()
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = dot(a, a).sqrt();
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Orthonormal tangent basis at unit vector c.
fn tangent_basis(c: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let a = if c[2].abs() < 0.9 { [0.0, 0.0, 1.0] } else { [1.0, 0.0, 0.0] };
    let e1 = normalize(cross(c, a));
    let e2 = cross(c, e1);
    (e1, e2)
}

/// Angular radius of a spherical cap holding `area_fraction` of the sphere.
fn cap_radius(area_fraction: f64) -> f64 {
    (1.0 - 2.0 * area_fraction.clamp(0.0, 1.0)).clamp(-1.0, 1.0).acos()
}

/// Point at angular distance `delta` from `c` along tangent direction `t`.
fn offset(c: [f64; 3], t: [f64; 3], delta: f64) -> [f64; 3] {
    normalize([
        c[0] * delta.cos() + t[0] * delta.sin(),
        c[1] * delta.cos() + t[1] * delta.sin(),
        c[2] * delta.cos() + t[2] * delta.sin(),
    ])
}

struct Craton {
    assembly: usize,
    centre: [f64; 3],
    major: [f64; 3],
    minor: [f64; 3],
    radius: f64,
    aspect: f64,
    area_m2: f64,
}

impl Craton {
    /// 1 - anisotropic angular distance / radius (1 at the centre, 0 on
    /// the nominal edge, negative outside).
    fn score(&self, x: [f64; 3]) -> f64 {
        let cd = dot(x, self.centre).clamp(-1.0, 1.0);
        let theta = cd.acos();
        let t = [
            x[0] - cd * self.centre[0],
            x[1] - cd * self.centre[1],
            x[2] - cd * self.centre[2],
        ];
        let tn = dot(t, t).sqrt();
        let eff = if tn < 1e-12 {
            theta
        } else {
            let c = dot(t, self.major) / tn;
            let s = dot(t, self.minor) / tn;
            theta * (c * c / self.aspect + self.aspect * s * s).sqrt()
        };
        1.0 - eff / self.radius
    }
}

/// Rank-size weights (k+1)^-a * exp(sigma g), sorted descending, summing
/// to one. Draws n gaussians (2n uniforms).
fn rank_shares(rng: &mut SplitMix64, n: usize, a: f64) -> Vec<f64> {
    let mut w: Vec<f64> = (0..n)
        .map(|k| ((k + 1) as f64).powf(-a) * (SIZE_JITTER_SIGMA * gaussian(rng)).exp())
        .collect();
    w.sort_by(|x, y| y.total_cmp(x));
    let s: f64 = w.iter().sum();
    w.iter().map(|v| v / s).collect()
}

/// 8-neighbour flat indices with longitude wrap (no rows past the poles).
fn nbrs8(nlat: usize, nlon: usize, c: usize, out: &mut Vec<usize>) {
    out.clear();
    let (i, j) = (c / nlon, c % nlon);
    for di in [-1i64, 0, 1] {
        let ti = i as i64 + di;
        if ti < 0 || ti >= nlat as i64 {
            continue;
        }
        for dj in [-1i64, 0, 1] {
            if di == 0 && dj == 0 {
                continue;
            }
            let tj = (j as i64 + dj).rem_euclid(nlon as i64) as usize;
            out.push(ti as usize * nlon + tj);
        }
    }
}

/// Connected components (8-neighbour, wrap) of the cells where `sel` is
/// true; returns per-cell label (usize::MAX outside) and per-label area.
fn components(nlat: usize, nlon: usize, sel: &[bool], areas: &[f64]) -> (Vec<usize>, Vec<f64>) {
    let n = sel.len();
    let mut lab = vec![usize::MAX; n];
    let mut comp_area = Vec::new();
    let mut stack = Vec::new();
    let mut nb = Vec::with_capacity(8);
    for s in 0..n {
        if !sel[s] || lab[s] != usize::MAX {
            continue;
        }
        let id = comp_area.len();
        let mut a = 0.0;
        lab[s] = id;
        stack.push(s);
        while let Some(c) = stack.pop() {
            a += areas[c];
            nbrs8(nlat, nlon, c, &mut nb);
            for &t in &nb {
                if sel[t] && lab[t] == usize::MAX {
                    lab[t] = id;
                    stack.push(t);
                }
            }
        }
        comp_area.push(a);
    }
    (lab, comp_area)
}

/// Cells by score descending (ties lowest id) until the area target is met
/// to within half of the last cell's area (the disc law's rule). Returns
/// the mask and the score of the last cell taken.
fn rank_select(score: &[f64], areas: &[f64], target: f64) -> (Vec<bool>, f64) {
    let mut order: Vec<usize> = (0..score.len()).collect();
    order.sort_by(|&a, &b| score[b].total_cmp(&score[a]).then(a.cmp(&b)));
    let mut mask = vec![false; score.len()];
    let mut acc = 0.0;
    let mut last = f64::INFINITY;
    for &c in &order {
        if acc + 0.5 * areas[c] > target {
            break;
        }
        mask[c] = true;
        acc += areas[c];
        last = score[c];
    }
    (mask, last)
}

/// Two-segment rifted margin thickness at coast distance d.
pub fn margin_thickness_m(h_int: f64, h_oc: f64, d: f64, width_m: f64, beta_p: f64) -> f64 {
    if !(width_m > 0.0) || d >= width_m {
        return h_int;
    }
    let wn = NECKING_FRACTION * width_m;
    let hinge = h_int / beta_p;
    let d = d.max(0.0);
    if d >= wn {
        hinge + (h_int - hinge) * (d - wn) / (width_m - wn)
    } else {
        h_oc + (hinge - h_oc) * d / wn
    }
}

/// `irregular-cratons-v1` initial crust (see the module doc). `seed` is
/// the tagged stream seed (`seed_from_path_tagged(path,
/// "irregular-cratons-v1")`).
pub fn seed_cratons(
    grid: &GridSpec,
    seeding: &PlateSeeding,
    prior: &UpliftPrior,
    area_fraction: f64,
    rifted_margin_width_m: f64,
    p: &CratonParams,
    assemblies: u32,
    seed: u64,
) -> Result<(CrustState, ContinentalBodiesReport, CratonSeedReport), KinematicError> {
    let n = grid.cells();
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    if seeding.assignment.len() != n {
        return Err(KinematicError::SizeMismatch {
            what: "plate assignment",
            expected: n,
            found: seeding.assignment.len(),
        });
    }
    if !(area_fraction > 0.0 && area_fraction < 1.0) {
        return Err(KinematicError::NonFinite {
            what: "continental area_fraction",
            value: area_fraction,
            index: 0,
        });
    }
    if assemblies == 0 {
        return Err(KinematicError::Grid("irregular-cratons-v1 needs at least one assembly".into()));
    }
    for (what, v) in [
        ("interior_thickness_m", p.interior_thickness_m),
        ("rank_size_exponent", p.rank_size_exponent),
        ("cratons_per_unit_share", p.cratons_per_unit_share),
        ("basement_roughness_std_m", p.basement_roughness_std_m),
        ("suture_crest_thickening_m", p.suture_crest_thickening_m),
        ("suture_half_width_m", p.suture_half_width_m),
        ("proximal_margin_stretch", p.proximal_margin_stretch),
        ("rifted_margin_width_m", rifted_margin_width_m),
    ] {
        if !v.is_finite() || v < 0.0 {
            return Err(KinematicError::NonFinite { what, value: v, index: 0 });
        }
    }
    if p.proximal_margin_stretch < 1.0 || p.suture_half_width_m <= 0.0 {
        return Err(KinematicError::Grid(
            "irregular-cratons-v1: proximal_margin_stretch must be >= 1 and suture_half_width_m > 0".into(),
        ));
    }
    let radius = grid.radius_m();
    let cents = centroids(grid)?;
    let areas = cell_areas(grid)?;
    let total_area: f64 = areas.iter().sum();
    let target = area_fraction * total_area;
    let mut rng = SplitMix64::new(seed);

    // 1. Assembly shares.
    let m = assemblies as usize;
    let shares = rank_shares(&mut rng, m, p.rank_size_exponent);
    let ext: Vec<f64> = shares.iter().map(|s| cap_radius(s * area_fraction)).collect();

    // 2. Assembly centres (min-spacing rejection; best candidate fallback).
    let mut centres: Vec<[f64; 3]> = Vec::with_capacity(m);
    let mut fallbacks = 0u32;
    for k in 0..m {
        let mut best: Option<([f64; 3], f64)> = None;
        let mut accepted = None;
        for _ in 0..MAX_CENTRE_DRAWS {
            let c = unit_dir(&mut rng);
            // worst ratio of separation to required separation
            let mut worst = f64::INFINITY;
            for (q, cq) in centres.iter().enumerate() {
                let sep = dot(c, *cq).clamp(-1.0, 1.0).acos();
                let need = ASSEMBLY_SPACING * (ext[k] + ext[q]);
                worst = worst.min(sep / need);
            }
            if worst >= 1.0 {
                accepted = Some(c);
                break;
            }
            if best.map_or(true, |(_, w)| worst > w) {
                best = Some((c, worst));
            }
        }
        match accepted {
            Some(c) => centres.push(c),
            None => {
                fallbacks += 1;
                centres.push(best.expect("at least one draw").0);
            }
        }
    }
    let assembly_plate: Vec<PlateId> = centres
        .iter()
        .map(|c| crate::kinematic::nearest_cell(grid, *c).map(|cell| seeding.assignment[cell]))
        .collect::<Result<_, _>>()?;

    // 3. Cratons.
    let mut cratons: Vec<Craton> = Vec::new();
    for k in 0..m {
        let nk = 1 + (p.cratons_per_unit_share * shares[k]).floor() as usize;
        let sub = rank_shares(&mut rng, nk, p.rank_size_exponent);
        let first = cratons.len();
        for (q, s) in sub.iter().enumerate() {
            let frac = shares[k] * s * area_fraction;
            let r = cap_radius(frac);
            let aspect = 1.0 + (ASPECT_MAX - 1.0) * rng.next_f64();
            let psi = 2.0 * PI * rng.next_f64();
            let centre = if q == 0 {
                centres[k]
            } else {
                let parent = first + ((rng.next_f64() * q as f64).floor() as usize).min(q - 1);
                let (e1, e2) = tangent_basis(cratons[parent].centre);
                let dir_a = 2.0 * PI * rng.next_f64();
                let t = [
                    e1[0] * dir_a.cos() + e2[0] * dir_a.sin(),
                    e1[1] * dir_a.cos() + e2[1] * dir_a.sin(),
                    e1[2] * dir_a.cos() + e2[2] * dir_a.sin(),
                ];
                offset(cratons[parent].centre, t, CRATON_OVERLAP * (cratons[parent].radius + r))
            };
            let (e1, e2) = tangent_basis(centre);
            let major = [
                e1[0] * psi.cos() + e2[0] * psi.sin(),
                e1[1] * psi.cos() + e2[1] * psi.sin(),
                e1[2] * psi.cos() + e2[2] * psi.sin(),
            ];
            let minor = cross(centre, major);
            cratons.push(Craton {
                assembly: k,
                centre,
                major,
                minor,
                radius: r,
                aspect,
                area_m2: frac * total_area,
            });
        }
    }

    // 4. Noise fields (fixed order).
    let warp: Vec<Spectral> = (0..3)
        .map(|_| Spectral::new(&mut rng, radius, &WARP_WAVELENGTHS_M, WARP_GAIN, WARP_AMPLITUDE_M / radius))
        .collect();
    let edge = Spectral::new(&mut rng, radius, &FBM_WAVELENGTHS_M, FBM_GAIN, EDGE_NOISE_STD);
    let rough = Spectral::new(&mut rng, radius, &FBM_WAVELENGTHS_M, FBM_GAIN, p.basement_roughness_std_m);
    let strike = Spectral::new(&mut rng, radius, &STRIKE_WAVELENGTHS_M, STRIKE_GAIN, STRIKE_STD);

    // 5. Scores: best and second-best craton per cell on the warped position.
    let mut best = vec![(0usize, f64::NEG_INFINITY); n];
    let mut second = vec![(0usize, f64::NEG_INFINITY); n];
    let mut score = vec![0.0f64; n];
    for c in 0..n {
        let x = cents[c];
        let xw = normalize([
            x[0] + warp[0].eval(x),
            x[1] + warp[1].eval(x),
            x[2] + warp[2].eval(x),
        ]);
        let (mut b, mut s2) = ((0usize, f64::NEG_INFINITY), (0usize, f64::NEG_INFINITY));
        for (j, cr) in cratons.iter().enumerate() {
            let s = cr.score(xw);
            if s > b.1 {
                s2 = b;
                b = (j, s);
            } else if s > s2.1 {
                s2 = (j, s);
            }
        }
        best[c] = b;
        second[c] = s2;
        score[c] = b.1 + edge.eval(x);
    }

    // 6. Mask: rank, fill small enclosed oceans, re-rank with them forced.
    let (mask0, _) = rank_select(&score, &areas, target);
    let ocean0: Vec<bool> = mask0.iter().map(|m| !m).collect();
    let (lab, comp_area) = components(nlat, nlon, &ocean0, &areas);
    let main = comp_area
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1).then(b.0.cmp(&a.0)))
        .map(|(i, _)| i);
    let fill_max = ENCLOSED_OCEAN_FILL_FRACTION * total_area;
    let mut forced_score = score.clone();
    for c in 0..n {
        if ocean0[c] && Some(lab[c]) != main && comp_area[lab[c]] < fill_max {
            forced_score[c] = f64::INFINITY;
        }
    }
    let (mut mask, threshold) = rank_select(&forced_score, &areas, target);
    // Second pass: holes opened by the re-rank are filled without re-ranking.
    let ocean1: Vec<bool> = mask.iter().map(|m| !m).collect();
    let (lab1, comp1) = components(nlat, nlon, &ocean1, &areas);
    let main1 = comp1
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1).then(b.0.cmp(&a.0)))
        .map(|(i, _)| i);
    let mut filled: Vec<bool> = forced_score.iter().map(|s| s.is_infinite()).collect();
    for c in 0..n {
        if ocean1[c] && Some(lab1[c]) != main1 && comp1[lab1[c]] < fill_max {
            mask[c] = true;
            filled[c] = true;
        }
    }
    let (mut filled_cells, mut filled_area) = (0u64, 0.0);
    for c in 0..n {
        if mask[c] && filled[c] {
            filled_cells += 1;
            filled_area += areas[c];
        }
    }

    // 7. Coast midpoints, bucketed by row for a windowed nearest search.
    let mut coast_by_row: Vec<Vec<[f64; 3]>> = vec![Vec::new(); nlat];
    let mut nb = Vec::with_capacity(8);
    for c in 0..n {
        if !mask[c] {
            continue;
        }
        nbrs8(nlat, nlon, c, &mut nb);
        for &t in &nb {
            if mask[t] {
                continue;
            }
            let raw = [cents[c][0] + cents[t][0], cents[c][1] + cents[t][1], cents[c][2] + cents[t][2]];
            let norm = dot(raw, raw).sqrt();
            if norm < 1e-12 {
                return Err(KinematicError::Grid("degenerate coast midpoint".into()));
            }
            coast_by_row[c / nlon].push([raw[0] / norm, raw[1] / norm, raw[2] / norm]);
        }
    }
    let dlat_m = PI * radius / nlat as f64;
    let row_window = (rifted_margin_width_m / dlat_m).ceil() as usize + 2;

    // 8. Thickness, plates, bookkeeping.
    let h_oc = prior.record.oceanic_crust_thickness_m;
    let mut plate = seeding.assignment.clone();
    let mut kind = Vec::with_capacity(n);
    let mut birth = Vec::with_capacity(n);
    let mut thick = Vec::with_capacity(n);
    let mut coast_distance = vec![f64::NAN; n];
    let mut h_int_v = vec![f64::NAN; n];
    let mut beta_v = vec![1.0f64; n];
    let mut craton_cells = vec![0u64; cratons.len()];
    let (mut suture_cells, mut clamped, mut relabelled) = (0u64, 0u64, 0u64);
    for c in 0..n {
        if !mask[c] {
            kind.push(CrustKindCodeV0::Oceanic);
            birth.push(-prior.initial_oceanic_age_my);
            thick.push(h_oc);
            continue;
        }
        let x = cents[c];
        let (j1, s1) = best[c];
        let (j2, s2) = second[c];
        craton_cells[j1] += 1;
        let k = cratons[j1].assembly;
        if plate[c] != assembly_plate[k] {
            relabelled += 1;
        }
        plate[c] = assembly_plate[k];
        // Suture root along the seam between the two best cratons.
        let mut ws = 0.0;
        if j2 != j1 && s2.is_finite() {
            let r1 = cratons[j1].radius;
            let r2 = cratons[j2].radius;
            let dist = (s1 - s2) / (1.0 / r1 + 1.0 / r2) * radius;
            let reach = ((s2 - (threshold - SUTURE_REACH_BAND)) / SUTURE_REACH_BAND).clamp(0.0, 1.0);
            let reach = reach * reach * (3.0 - 2.0 * reach);
            ws = (-0.5 * (dist / p.suture_half_width_m).powi(2)).exp() * reach;
        }
        if ws > 0.5 {
            suture_cells += 1;
        }
        let mut h_int = p.interior_thickness_m
            + rough.eval(x)
            + p.suture_crest_thickening_m * ws * (1.0 + strike.eval(x)).max(0.0);
        if h_int < MIN_INTERIOR_THICKNESS_M {
            h_int = MIN_INTERIOR_THICKNESS_M;
            clamped += 1;
        }
        // Nearest coast midpoint within the row window.
        let i = c / nlon;
        let lo = i.saturating_sub(row_window);
        let hi = (i + row_window).min(nlat - 1);
        let mut bestd = -2.0f64;
        for row in &coast_by_row[lo..=hi] {
            for mpt in row {
                let d = dot(x, *mpt);
                if d > bestd {
                    bestd = d;
                }
            }
        }
        let dist = if bestd < -1.5 {
            f64::INFINITY
        } else {
            bestd.clamp(-1.0, 1.0).acos() * radius
        };
        coast_distance[c] = dist;
        let h = margin_thickness_m(h_int, h_oc, dist, rifted_margin_width_m, p.proximal_margin_stretch);
        if h < h_int {
            beta_v[c] = h_int / h;
        }
        h_int_v[c] = h_int;
        kind.push(CrustKindCodeV0::Continental);
        birth.push(-prior.initial_continental_age_my);
        thick.push(h);
    }

    // Bodies report: one body per assembly.
    let mut bodies = Vec::with_capacity(m);
    let (mut achieved, mut volume) = (0.0, 0.0);
    for k in 0..m {
        let mut rec = ContinentalBodyRecord {
            plate: assembly_plate[k],
            seed_point: centres[k],
            target_area_m2: shares[k] * target,
            achieved_area_m2: 0.0,
            cells: 0,
            crust_volume_m3: 0.0,
            tapered_cells: 0,
        };
        for c in 0..n {
            if mask[c] && cratons[best[c].0].assembly == k {
                rec.achieved_area_m2 += areas[c];
                rec.cells += 1;
                rec.crust_volume_m3 += thick[c] * areas[c];
                if beta_v[c] > 1.0 {
                    rec.tapered_cells += 1;
                }
            }
        }
        achieved += rec.achieved_area_m2;
        volume += rec.crust_volume_m3;
        bodies.push(rec);
    }
    let interior_equiv: f64 = (0..n).filter(|&c| mask[c]).map(|c| h_int_v[c] * areas[c]).sum();
    let report = ContinentalBodiesReport {
        declared_area_fraction: area_fraction,
        achieved_area_fraction: achieved / total_area,
        rifted_margin_width_m,
        interior_thickness_m: p.interior_thickness_m,
        edge_thickness_m: h_oc,
        continental_crust_volume_m3: volume,
        plateau_equivalent_volume_m3: interior_equiv,
        coast_midpoints: coast_by_row.iter().map(|r| r.len() as u64).sum(),
        coast_distance_m: coast_distance,
        bodies,
    };
    let extra = CratonSeedReport {
        cratons: cratons
            .iter()
            .enumerate()
            .map(|(j, cr)| CratonRecord {
                assembly: cr.assembly,
                plate: assembly_plate[cr.assembly],
                centre: cr.centre,
                radius_rad: cr.radius,
                aspect: cr.aspect,
                area_target_m2: cr.area_m2,
                cells: craton_cells[j],
            })
            .collect(),
        assemblies: (0..m).map(|k| (assembly_plate[k], centres[k], shares[k])).collect(),
        spacing_fallbacks: fallbacks,
        filled_cells,
        filled_area_m2: filled_area,
        suture_cells,
        interior_clamped_cells: clamped,
        score_threshold: threshold,
        relabelled_cells: relabelled,
        interior_thickness_m: h_int_v,
        stretch_beta: beta_v,
    };
    let state = CrustState {
        t_my: 0.0,
        plate,
        kind,
        birth_t_my: birth,
        thickness_m: thick,
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: cents,
    };
    Ok((state, report, extra))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn margin_profile_is_continuous_and_monotone() {
        let (h_int, h_oc, w, b) = (29_600.0, 7_000.0, 200_000.0, 1.15);
        assert_eq!(margin_thickness_m(h_int, h_oc, 0.0, w, b), h_oc);
        let wn = NECKING_FRACTION * w;
        let hinge = h_int / b;
        assert!((margin_thickness_m(h_int, h_oc, wn, w, b) - hinge).abs() < 1e-9);
        assert!((margin_thickness_m(h_int, h_oc, wn - 1e-6, w, b) - hinge).abs() < 1e-3);
        assert_eq!(margin_thickness_m(h_int, h_oc, w, w, b), h_int);
        let mut prev = 0.0;
        for k in 0..=250 {
            let h = margin_thickness_m(h_int, h_oc, k as f64 * 1000.0, w, b);
            assert!(h >= prev - 1e-9);
            prev = h;
        }
        assert_eq!(margin_thickness_m(h_int, h_oc, 5.0, 0.0, b), h_int);
    }

    #[test]
    fn rank_shares_sum_to_one_and_descend() {
        let mut rng = SplitMix64::new(7);
        let s = rank_shares(&mut rng, 5, 1.4);
        assert!((s.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(s.windows(2).all(|w| w[0] >= w[1]));
    }

    #[test]
    fn spectral_noise_has_roughly_the_declared_std() {
        let mut rng = SplitMix64::new(11);
        let f = Spectral::new(&mut rng, 6.371e6, &FBM_WAVELENGTHS_M, FBM_GAIN, 1500.0);
        let mut r2 = SplitMix64::new(99);
        let (mut s, mut s2, n) = (0.0, 0.0, 4000);
        for _ in 0..n {
            let v = f.eval(unit_dir(&mut r2));
            s += v;
            s2 += v * v;
        }
        let sd = (s2 / n as f64 - (s / n as f64).powi(2)).sqrt();
        assert!(sd > 1000.0 && sd < 2000.0, "std {sd}");
    }
}
