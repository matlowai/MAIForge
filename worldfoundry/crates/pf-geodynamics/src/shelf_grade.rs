//! SHELF (lane SHELF, 2026-10-20): `wave-base-shelf-grading-v1`, an opt-in
//! marine reworking step that keeps continental shelves graded to the
//! CURRENT wave base while sea level moves.
//!
//! DEFECT IT ADDRESSES (lane SHELF measurements on the C6 world, seed 42,
//! `unpublished-design-record`): shelves form only when a sea-level rise
//! floods the flat continental platform (shelf 4-7 % of the globe right
//! after each breakup) and then vanish as the new ocean ages and sea level
//! falls (back to 0.6-1 % within 40-60 My): the flooded platform, filled to
//! wave base by `continental-retention-base-level-v1`, simply re-emerges as
//! a flat coastal plain 100-200 m above the sea, and nothing grades a new
//! shelf at the lower sea level because the continental edge drops to the
//! abyss within one cell (passive-margin first offshore ring p50 -4.8 km
//! at 400 My vs Earth -0.13 km at the same 0.625 degree grid).
//!
//! THE LAW. Each step, after erosion and deposition, every cell of the
//! connected ocean whose surface stands between sea level and wave base
//! (`sea - z_wb < e < sea`) loses unconsolidated SEDIMENT (never bedrock)
//! down to wave base, and the same volume is laid down nearest first
//! (Dijkstra by centre distance through ocean cells) on ocean cells below
//! wave base, each filled at most to wave base: the shelf edge progrades
//! seaward. Thickness and surface change are related by the water-loaded
//! local Airy response of a sediment layer, `de = ds (rho_m - rho_s) /
//! (rho_m - rho_w)`, the same compensated column the retention law uses.
//! Volume that cannot be placed within the search bound stays where it
//! was (never created or destroyed: every removal has an equal placement).
//!
//! Physical basis (recalled, citations not checked): storm-wave base of
//! order 100 m grades shelf tops (shelf equilibrium profile; Swift and
//! Thorne 1991, Shelf Sand and Sandstone Bodies, IAS SP 14:3); during slow
//! relative sea-level fall the shoreface and shelf top are eroded to the new
//! wave base and the material is shed onto the slope (regressive surface of
//! marine erosion and forced-regressive clinoforms; Plint 1988, SEPM SP 42).
//! Tectonic eustatic rates here are metres per My, orders of magnitude
//! slower than shoreface reworking, so grading within one 2 My step is the
//! declared approximation (project assumption awaiting review).
//!
//! Constants: none new. Wave base is the retention law's
//! `shelf_wave_base_depth_m`; densities are the isostasy law's; the search
//! bound is a computational bound passed by the caller.
//!
//! Declared limits: bedrock (crystalline or lithified) is never eroded by
//! waves; lakes and enclosed seas are not graded (connected ocean only);
//! local Airy response only (no flexure); no along-shore drift; grain size
//! and compaction are not modelled; the sediment ledger sees an internal
//! transfer only (eroded and deposited volumes of the step are unchanged).
//!
//! Determinism: sources in ascending cell order, Dijkstra ties to the
//! lowest cell id; no HashMap.

use pf_grid::GridSpec;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Law id (typed law set entry is the integrator's pf-contracts change).
pub const SHELF_GRADING_LAW_ID: &str = "wave-base-shelf-grading-v1";

/// Parameters (all from existing laws; no new physical constant).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShelfGradeParams {
    /// Wave base depth below sea level (m), the retention law's value.
    pub wave_base_depth_m: f64,
    pub mantle_density_kg_m3: f64,
    pub sediment_density_kg_m3: f64,
    pub seawater_density_kg_m3: f64,
    /// Computational bound on cells settled per source.
    pub search_max_cells: usize,
}

impl ShelfGradeParams {
    /// Water-loaded surface change per metre of sediment (dimensionless).
    pub fn lift_per_m(&self) -> f64 {
        (self.mantle_density_kg_m3 - self.sediment_density_kg_m3)
            / (self.mantle_density_kg_m3 - self.seawater_density_kg_m3)
    }

    pub fn validate(&self) -> Result<(), String> {
        for (what, v) in [
            ("wave_base_depth_m", self.wave_base_depth_m),
            ("mantle_density_kg_m3", self.mantle_density_kg_m3),
            ("sediment_density_kg_m3", self.sediment_density_kg_m3),
            ("seawater_density_kg_m3", self.seawater_density_kg_m3),
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(format!("{SHELF_GRADING_LAW_ID}: {what} must be finite and > 0, got {v}"));
            }
        }
        if !(self.mantle_density_kg_m3 > self.sediment_density_kg_m3
            && self.mantle_density_kg_m3 > self.seawater_density_kg_m3)
        {
            return Err(format!("{SHELF_GRADING_LAW_ID}: mantle must be denser than sediment and sea water"));
        }
        if self.search_max_cells == 0 {
            return Err(format!("{SHELF_GRADING_LAW_ID}: search_max_cells must be > 0"));
        }
        Ok(())
    }
}

/// What one grading pass did (volumes in m^3).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ShelfGradeReport {
    /// Ocean cells shallower than wave base that held sediment.
    pub sources: u64,
    /// Sediment volume moved (removed = placed).
    pub moved_m3: f64,
    /// Volume that could not be placed within the search bound and stayed.
    pub kept_m3: f64,
    /// Cells that received sediment.
    pub receivers: u64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Key(f64, usize);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    /// Min-heap on distance, then on cell id.
    fn cmp(&self, other: &Self) -> Ordering {
        other.0.total_cmp(&self.0).then_with(|| other.1.cmp(&self.1))
    }
}

/// One grading pass. Inputs: `surface_m` (the post-erosion surface,
/// relative to the same datum as `sea_level_m`), `ocean` (connected ocean
/// mask), `sediment_m` (sediment available per cell after the erosion step,
/// m) and `areas_m2`. Outputs are ADDED to `removed_m` / `added_m` (per-cell
/// sediment thickness, m) and `dz_m` (surface change, m), so a caller can
/// fold them into its erosion step. Returns the report.
#[allow(clippy::too_many_arguments)]
pub fn grade_to_wave_base(
    grid: &GridSpec,
    areas_m2: &[f64],
    surface_m: &[f64],
    ocean: &[bool],
    sediment_m: &[f64],
    sea_level_m: f64,
    p: &ShelfGradeParams,
    removed_m: &mut [f64],
    added_m: &mut [f64],
    dz_m: &mut [f64],
) -> Result<ShelfGradeReport, String> {
    p.validate()?;
    let n = grid.cells();
    for (what, len) in [
        ("areas", areas_m2.len()),
        ("surface", surface_m.len()),
        ("ocean", ocean.len()),
        ("sediment", sediment_m.len()),
        ("removed", removed_m.len()),
        ("added", added_m.len()),
        ("dz", dz_m.len()),
    ] {
        if len != n {
            return Err(format!("{SHELF_GRADING_LAW_ID}: {what} has {len} cells, grid {n}"));
        }
    }
    if !sea_level_m.is_finite() {
        return Err(format!("{SHELF_GRADING_LAW_ID}: sea level {sea_level_m}"));
    }
    let f = p.lift_per_m();
    let target = sea_level_m - p.wave_base_depth_m;
    let mut e: Vec<f64> = surface_m.to_vec();
    let mut sed: Vec<f64> = sediment_m.to_vec();
    let mut nb_cache: Vec<Option<Vec<(usize, f64)>>> = vec![None; n];
    let mut settled = vec![0u32; n];
    let mut gen = 0u32;
    let mut heap = BinaryHeap::new();
    let mut received = vec![false; n];
    let mut rep = ShelfGradeReport::default();
    for c in 0..n {
        if !ocean[c] || !(e[c] > target && e[c] < sea_level_m) || !(sed[c] > 0.0) {
            continue;
        }
        let ds = sed[c].min((e[c] - target) / f);
        if !(ds > 0.0) {
            continue;
        }
        rep.sources += 1;
        let mut remaining = ds * areas_m2[c];
        // Dijkstra over ocean cells from c; c itself is not a receiver.
        gen = gen.wrapping_add(1);
        if gen == 0 {
            settled.iter_mut().for_each(|s| *s = 0);
            gen = 1;
        }
        heap.clear();
        heap.push(Key(0.0, c));
        let mut count = 0usize;
        let mut placed: Vec<(usize, f64)> = Vec::new();
        while let Some(Key(d, x)) = heap.pop() {
            if remaining <= 0.0 {
                break;
            }
            if settled[x] == gen {
                continue;
            }
            if count >= p.search_max_cells {
                break;
            }
            settled[x] = gen;
            count += 1;
            if x != c && e[x] < target {
                let need = (target - e[x]) / f * areas_m2[x];
                let v = need.min(remaining);
                placed.push((x, v));
                remaining -= v;
                if remaining <= 0.0 {
                    break;
                }
            }
            if nb_cache[x].is_none() {
                let (i, j) = grid.ij(x as u64).map_err(|e| e.to_string())?;
                let mut list: Vec<(usize, f64)> = grid
                    .neighbors(i, j)
                    .map_err(|e| e.to_string())?
                    .into_iter()
                    .map(|(ti, tj, nb)| (ti * grid.nlon() + tj, nb.distance_m.max(1.0)))
                    .collect();
                list.sort_unstable_by_key(|&(k, _)| k);
                nb_cache[x] = Some(list);
            }
            for &(y, dist) in nb_cache[x].as_ref().expect("filled above") {
                if ocean[y] && settled[y] != gen {
                    heap.push(Key(d + dist, y));
                }
            }
        }
        let moved: f64 = placed.iter().map(|&(_, v)| v).sum();
        rep.kept_m3 += ds * areas_m2[c] - moved;
        if !(moved > 0.0) {
            continue;
        }
        let dsc = moved / areas_m2[c];
        sed[c] -= dsc;
        e[c] -= dsc * f;
        removed_m[c] += dsc;
        dz_m[c] -= dsc * f;
        for (x, v) in placed {
            let h = v / areas_m2[x];
            sed[x] += h;
            e[x] += h * f;
            added_m[x] += h;
            dz_m[x] += h * f;
            if !received[x] {
                received[x] = true;
                rep.receivers += 1;
            }
        }
        rep.moved_m3 += moved;
    }
    Ok(rep)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> ShelfGradeParams {
        ShelfGradeParams {
            wave_base_depth_m: 100.0,
            mantle_density_kg_m3: 3300.0,
            sediment_density_kg_m3: 2800.0,
            seawater_density_kg_m3: 1030.0,
            search_max_cells: 64,
        }
    }

    #[test]
    fn shoal_is_graded_and_volume_moves_offshore() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        let areas: Vec<f64> = (0..n).map(|c| 1.0e9 + c as f64).collect();
        // Row 4: a shoal at -40 m with 1 km of sediment next to deep water.
        let mut surf = vec![-3000.0; n];
        let shoal = 4 * 16 + 5;
        surf[shoal] = -40.0;
        let ocean = vec![true; n];
        let mut sed = vec![0.0; n];
        sed[shoal] = 1000.0;
        let (mut rm, mut ad, mut dz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        let p = params();
        let rep = grade_to_wave_base(&grid, &areas, &surf, &ocean, &sed, 0.0, &p, &mut rm, &mut ad, &mut dz).unwrap();
        assert_eq!(rep.sources, 1);
        // Graded exactly to wave base.
        assert!((surf[shoal] + dz[shoal] + 100.0).abs() < 1e-6, "{}", surf[shoal] + dz[shoal]);
        // Mass: removed volume equals added volume.
        let out: f64 = (0..n).map(|c| rm[c] * areas[c]).sum();
        let inn: f64 = (0..n).map(|c| ad[c] * areas[c]).sum();
        assert!((out - inn).abs() < 1e-3 * out.max(1.0));
        assert!((out - rep.moved_m3).abs() < 1e-3 * out);
        // Nothing raised above wave base.
        for c in 0..n {
            assert!(surf[c] + dz[c] <= -100.0 + 1e-6 || c == shoal);
        }
    }

    #[test]
    fn land_lakes_and_bedrock_are_untouched() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        let areas = vec![1.0e9; n];
        let mut surf = vec![-3000.0; n];
        surf[10] = 50.0; // land
        surf[11] = -40.0; // shoal without sediment (bedrock)
        surf[12] = -40.0; // shoal in a lake (not connected ocean)
        let mut ocean = vec![true; n];
        ocean[10] = false;
        ocean[12] = false;
        let mut sed = vec![0.0; n];
        sed[10] = 500.0;
        sed[12] = 500.0;
        let (mut rm, mut ad, mut dz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        let rep = grade_to_wave_base(&grid, &areas, &surf, &ocean, &sed, 0.0, &params(), &mut rm, &mut ad, &mut dz).unwrap();
        assert_eq!(rep.sources, 0);
        assert!(rm.iter().all(|v| *v == 0.0) && ad.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn thin_cover_only_partly_grades() {
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        let areas = vec![1.0e9; n];
        let mut surf = vec![-3000.0; n];
        surf[40] = -20.0;
        let ocean = vec![true; n];
        let mut sed = vec![0.0; n];
        sed[40] = 50.0; // removes at most 50 m of sediment (11 m of surface)
        let (mut rm, mut ad, mut dz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        let p = params();
        grade_to_wave_base(&grid, &areas, &surf, &ocean, &sed, 0.0, &p, &mut rm, &mut ad, &mut dz).unwrap();
        assert!((rm[40] - 50.0).abs() < 1e-9);
        assert!((dz[40] + 50.0 * p.lift_per_m()).abs() < 1e-9);
    }
}
