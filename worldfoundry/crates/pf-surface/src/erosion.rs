//! Erosion for the PREVIEW tier (ADR 032): detachment-limited stream
//! power plus linear hillslope diffusion, on the real spherical metric,
//! deterministic, with an explicit stability condition.
//!
//! STREAM POWER: dz/dt = -K A^m S^n with n = 1, solved IMPLICITLY in the
//! receiver-first order of Braun and Willett 2013 (Geomorphology 180:170):
//!     z_i' = (z_i + f z_j') / (1 + f),   f = K A_i^m dt / L_ij,
//! which is unconditionally stable and never erodes a cell below its
//! receiver. Cells that sit below their receiver on the physical surface
//! (lake floors inside a filled depression) receive no fluvial change.
//! Detached mass travels down the flow path to the first LAKE it meets
//! (a cell whose routing fill level stands above the physical surface),
//! else to the ocean cell or terminal floor the path ends in, and is laid
//! across that water body by a deterministic bounded breadth-first fill up
//! to the water level (the lake's spill level from the fill surface, or
//! sea level). A full lake passes the remainder downstream to the next
//! lake or the ocean; the ocean body aggrades uniformly when the local
//! bound is exhausted; only a terminal floor with no outlet aggrades its
//! basin. Lake deposition counts as DEPOSITED, ocean receipts as EXPORTED;
//! the uniformly aggraded volume is reported as `deposit_aggraded_m3`.
//! Submerged cells (under a lake) receive no stream-power lowering.
//! Without the spread a full basin builds a one-cell pile tens of km
//! thick per step (observed 35 to 45 km of elevation in the first
//! level-0 runs).
//!
//! HILLSLOPE DIFFUSION: dz/dt = D lap(z), finite-volume flux form on the
//! four edge neighbours with great-circle edge lengths and distances,
//! explicit in time. The diffusion number nu_i = dt * sum_j D w_ij /
//! (d_ij A_i) must satisfy nu_i <= `stability_limit` (declared 0.5, at
//! which the update is a convex combination); the step is subdivided into
//! ceil(max nu / limit) substeps, and both numbers are recorded in the
//! artifact. Ocean cells are one-way sinks (flux in counts as export).
//!
//! MASS: every change is a flux between two cells or a transfer to a
//! terminal, so eroded = deposited + exported closes to floating-point
//! rounding, which `sediment.rs` declares as the error. Erosion removes
//! sediment before bedrock; only the bedrock share leaves the crust.
//!
//! AVAILABLE-MATERIAL BOUND (lane E, SEED-EROSION-BUG, 2026-10-11):
//! `erode_step_bounded` takes the bedrock column each cell holds at the
//! start of the step and never removes more than that column plus the
//! sediment present at the moment of removal. Physical basis: erosion
//! removes rock that exists; a detachment-limited law (Whipple and Tucker
//! 1999, JGR 104:17661, recalled, UNVERIFIED) is a RATE law and says
//! nothing about the supply of material, so the supply limit is a separate
//! constraint, applied the way supply-limited regolith/bedrock models apply
//! it (erosion = min(capacity-driven rate * dt, material present)). In the
//! receiver-first implicit pass the bound lowers a cell by at most its
//! available column, so every donor processed later sees the bounded
//! receiver level (the order of Braun and Willett 2013 is kept). In
//! hillslope diffusion each donor's total outgoing flux in a substep is
//! scaled by min(1, available volume / outgoing volume); the scaled flux
//! is still a two-cell transfer, so mass closure is unchanged. `erode_step`
//! is `erode_step_bounded` with no column supplied: the legacy result,
//! bit for bit (guarded by tests). Measured on seeds 240907, 7 and 3 at
//! L1 (toy preset, 100 My) the bound never binds: the 156 to 208 km of
//! cumulative erosion in single cells is collisional supply passing
//! through a cap-limited column, not overdraw (lane E report,
//! `unpublished-design-record`).
//!
//! ACCOUNTING: `bedrock_eroded_m3` and `sediment_eroded_m3` split
//! `eroded_m3` (accumulated in the same operation, so their sum equals it
//! to rounding) so a caller can close the bedrock and sediment control
//! volumes separately (see `sediment::ControlVolumeLedger`).
//!
//! CONSTANTS: K and the exponents are ASSUMPTIONS awaiting technical review
//! (model-assumption-review Q2, Q3); the driver records them in the artifact.

use crate::routing::Routing;
use pf_contracts::spatial::ElementId;
use pf_grid::GridSpec;

/// SEED-SEDIMENT-RETAIN (lane SR, 2026-10-13): `continental-retention-base-level-v1`.
pub mod retain;
/// Q25 / D-0155: opt-in connected-ocean terrigenous fines.
pub mod marine_fines;
pub mod fluvial_trace;
pub mod surface_budget;
/// Lane INCISION-THRESHOLD (2026-10-20): `stream-power-excess-steepness-threshold-v1`.
pub mod incision_threshold;
/// GLACIAL-T0 (2026-10-21): opt-in `glacial-stream-power-blend-v1`.
pub mod glacial;

#[derive(Debug, Clone, PartialEq)]
pub enum ErosionError {
    SizeMismatch {
        what: &'static str,
        expected: usize,
        found: usize,
    },
    NonFinite {
        what: &'static str,
        value: f64,
        index: usize,
    },
    BadParam {
        what: &'static str,
        value: f64,
    },
    Grid(String),
    /// Lane SD: a surface-law parameter set refused by its own validation.
    Law(String),
}

impl std::fmt::Display for ErosionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ErosionError::SizeMismatch {
                what,
                expected,
                found,
            } => write!(f, "{what}: expected {expected}, found {found}"),
            ErosionError::NonFinite { what, value, index } => {
                write!(f, "nonfinite {what} {value:?} at {index}")
            }
            ErosionError::BadParam { what, value } => write!(f, "bad parameter {what} = {value:?}"),
            ErosionError::Grid(why) => write!(f, "grid: {why}"),
            ErosionError::Law(why) => write!(f, "surface law: {why}"),
        }
    }
}

impl std::error::Error for ErosionError {}

/// Erosion constants (per My, A in m^2, S dimensionless).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ErosionParams {
    pub k_per_my: f64,
    pub m: f64,
    /// Must be 1.0 (the implicit scheme is declared for n = 1 only).
    pub n: f64,
    pub diffusivity_m2_per_my: f64,
    /// Declared explicit-diffusion limit in (0, 0.5].
    pub stability_limit: f64,
}

impl ErosionParams {
    pub fn validate(&self) -> Result<(), ErosionError> {
        for (what, v) in [
            ("k_per_my", self.k_per_my),
            ("m", self.m),
            ("diffusivity_m2_per_my", self.diffusivity_m2_per_my),
        ] {
            if !v.is_finite() || v < 0.0 {
                return Err(ErosionError::BadParam { what, value: v });
            }
        }
        if self.n != 1.0 {
            return Err(ErosionError::BadParam {
                what: "n (implicit scheme declared for n = 1 only)",
                value: self.n,
            });
        }
        if !self.stability_limit.is_finite()
            || self.stability_limit <= 0.0
            || self.stability_limit > 0.5
        {
            return Err(ErosionError::BadParam {
                what: "stability_limit (0, 0.5]",
                value: self.stability_limit,
            });
        }
        Ok(())
    }
}

/// Deposition spread bound: cells visited per terminal before the
/// remainder is dumped on the terminal itself (declared, reported).
pub const MAX_SPREAD_CELLS: usize = 512;

/// How the marine remainder (river-mouth volume left after the bounded
/// local fill to sea level) is laid as the declared global layer
/// (lane INCISION-THRESHOLD part 1, 2026-10-20, after lane EROSION-ORACLE
/// finding 3). Selected per call; every existing entry point passes
/// `PerMouthV0`, so default outputs are byte-identical.
///
/// `PerMouthV0` (legacy, `per-mouth-global-layer-v0`): each saturated mouth
/// raises EVERY ocean cell by its own remainder / ocean area, in mouth
/// order. Cost O(saturated mouths x ocean cells) per step (lane
/// EROSION-ORACLE measured the erode step rising 0.03 -> 0.88 s/step on a
/// 1.04M-cell grid), and cells the local fill had already brought to sea
/// level are pushed ABOVE it (measured up to 0.31 m there). The bounded
/// fill can also overshoot the water level by ~1 ulp (4.5e-13 m measured).
///
/// `OncePerStepV1` (`once-per-step-capped-global-layer-v1`): the remainders
/// of all mouths are summed and laid ONCE per step, as a uniform layer over
/// the ocean cells that still stand below sea level, each capped at its
/// room to sea level (water-filling: capped cells are filled exactly to sea
/// level and the rest is spread over the others). Every local fill lands
/// exactly on its water level (no ulp overshoot). Same volume, same
/// "declared global layer" idealisation; cost O(ocean cells x fill rounds).
/// Physical basis: a far-field marine blanket cannot stand above the sea
/// surface that bounds its accommodation. Order effect (declared): later
/// mouths' local fills no longer see earlier mouths' blanket (a few cm to m).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OceanRemainderLayer {
    #[default]
    PerMouthV0,
    OncePerStepV1,
}

pub const OCEAN_REMAINDER_LAYER_V0_ID: &str = "per-mouth-global-layer-v0";
pub const OCEAN_REMAINDER_LAYER_V1_ID: &str = "once-per-step-capped-global-layer-v1";

impl OceanRemainderLayer {
    pub fn law_id(&self) -> &'static str {
        match self {
            OceanRemainderLayer::PerMouthV0 => OCEAN_REMAINDER_LAYER_V0_ID,
            OceanRemainderLayer::OncePerStepV1 => OCEAN_REMAINDER_LAYER_V1_ID,
        }
    }
}

/// `OncePerStepV1` global layer: thickness per ocean cell (m) for `volume_m3`
/// spread uniformly over the ocean cells below `sea_level_m`, each capped at
/// its room to sea level (iterated water-filling). Returns (cell, thickness,
/// filled_to_sea_level) for every cell that receives material, and the
/// number of fill rounds. If every ocean cell is full (an ocean with no
/// accommodation left), the rest is laid uniformly over all ocean cells as
/// the legacy layer (flagged by `overfull = true`); not reached in practice.
pub(crate) fn capped_ocean_layer(
    z: &[f64],
    areas: &[f64],
    ocean_cells: &[usize],
    sea_level_m: f64,
    volume_m3: f64,
) -> (Vec<(usize, f64, bool)>, u32, bool) {
    let mut out: Vec<(usize, f64, bool)> = Vec::new();
    if !(volume_m3 > 0.0) {
        return (out, 0, false);
    }
    let mut active: Vec<usize> = ocean_cells
        .iter()
        .copied()
        .filter(|&c| z[c] < sea_level_m)
        .collect();
    let mut remaining = volume_m3;
    let mut rounds = 0u32;
    loop {
        rounds += 1;
        if active.is_empty() {
            let area: f64 = ocean_cells.iter().map(|&c| areas[c]).sum();
            let layer = remaining / area;
            out.extend(ocean_cells.iter().map(|&c| (c, layer, false)));
            return (out, rounds, true);
        }
        let area: f64 = active.iter().map(|&c| areas[c]).sum();
        let layer = remaining / area;
        let mut keep: Vec<usize> = Vec::with_capacity(active.len());
        let mut filled_any = false;
        for &c in &active {
            let room = sea_level_m - z[c];
            if room <= layer {
                out.push((c, room, true));
                remaining -= room * areas[c];
                filled_any = true;
            } else {
                keep.push(c);
            }
        }
        if !filled_any {
            out.extend(active.iter().map(|&c| (c, layer, false)));
            return (out, rounds, false);
        }
        if !(remaining > 0.0) {
            return (out, rounds, false);
        }
        active = keep;
    }
}

/// Result of one erosion step. Thickness changes are per cell; the mass
/// terms are volumes (m^3) the caller converts with its rock density.
#[derive(Debug, Clone, PartialEq)]
pub struct ErosionStep {
    pub transfers: Vec<pf_contracts::domain::WaterTransferV0>,
    pub spread_bound_hits: u64,
    /// Surface change per cell (m): negative = lowered.
    pub dz_m: Vec<f64>,
    /// Bedrock removed per cell (m, >= 0).
    pub bedrock_removed_m: Vec<f64>,
    /// Sediment removed per cell (m, >= 0).
    pub sediment_removed_m: Vec<f64>,
    /// Sediment added per cell (m, >= 0): land deposition and ocean receipts.
    pub sediment_added_m: Vec<f64>,
    /// Volumes for the ledger (m^3).
    pub eroded_m3: f64,
    pub deposited_m3: f64,
    pub exported_m3: f64,
    /// Explicit stability accounting.
    pub diffusion_number_max: f64,
    pub diffusion_substeps: u32,
    /// Max surface lowering rate this step (m/My).
    pub max_lowering_rate_m_per_my: f64,
    /// Volume (m^3) laid down as a uniform aggradation layer over the
    /// terminal's basin or the ocean body after the local spread was
    /// exhausted (diagnostic).
    pub deposit_aggraded_m3: f64,
    /// Bedrock share of `eroded_m3` (m^3).
    pub bedrock_eroded_m3: f64,
    /// Sediment share of `eroded_m3` (m^3), including sediment deposited
    /// and re-eroded within this step.
    pub sediment_eroded_m3: f64,
    /// Lowerings (stream power, per cell) and diffusion donors (per
    /// substep) limited by the available-material bound. Always 0 when no
    /// column is supplied.
    pub bedrock_bound_hits: u64,
    /// Lane SD: where the moved material came to rest, under the selected
    /// surface laws (`erode_deposit_step`). Always `None` on the legacy
    /// path (`erode_step`, `erode_step_bounded`).
    pub deposition: Option<crate::sediment::DepositionReport>,
    /// DIAGNOSTIC (lane INCISION-THRESHOLD part 1): mouths whose remainder
    /// went to the global marine layer this step.
    pub ocean_layer_mouths: u64,
    /// DIAGNOSTIC: ocean-cell raises by the global layer (cost proxy:
    /// mouths x ocean cells under `PerMouthV0`, fill-round cells under V1).
    pub ocean_layer_cell_raises: u64,
    /// DIAGNOSTIC: ocean cells the global layer lifted from at or below sea
    /// level to above it (always 0 under `OncePerStepV1`, except the flagged
    /// overfull fallback).
    pub ocean_layer_cells_above_sea: u64,
}

/// Precomputed diffusion stencil: for each cell, its edge neighbours with
/// the coefficient w_ij / d_ij (dimensionless) and the cell area.
pub struct DiffusionStencil {
    pub neighbors: Vec<Vec<(usize, f64)>>,
    pub areas_m2: Vec<f64>,
}

impl DiffusionStencil {
    pub fn build(grid: &GridSpec) -> Result<Self, ErosionError> {
        let n = grid.cells();
        let radius = grid.radius_m();
        let dlat = std::f64::consts::PI / grid.nlat() as f64;
        let dlon = 2.0 * std::f64::consts::PI / grid.nlon() as f64;
        let mut neighbors = Vec::with_capacity(n);
        let mut areas = Vec::with_capacity(n);
        for i in 0..grid.nlat() {
            let area = grid
                .cell_area_m2(i)
                .map_err(|e| ErosionError::Grid(e.to_string()))?;
            let lat_bot = -std::f64::consts::FRAC_PI_2 + i as f64 * dlat;
            let lat_top = lat_bot + dlat;
            for j in 0..grid.nlon() {
                let mut list = Vec::with_capacity(4);
                for (ti, tj, nb) in grid
                    .neighbors(i, j)
                    .map_err(|e| ErosionError::Grid(e.to_string()))?
                {
                    // Edge neighbours only (|di| + |dj| == 1).
                    if nb.di.abs() + nb.dj.abs() != 1 {
                        continue;
                    }
                    // Shared edge length: E-W neighbours share a meridian
                    // segment R dlat; N-S neighbours share a parallel
                    // segment R cos(lat_edge) dlon.
                    let w = if nb.di == 0 {
                        radius * dlat
                    } else if nb.di > 0 {
                        radius * lat_top.cos() * dlon
                    } else {
                        radius * lat_bot.cos() * dlon
                    };
                    let d = nb.distance_m.max(1.0);
                    list.push((ti * grid.nlon() + tj, w / d));
                }
                neighbors.push(list);
                areas.push(area);
            }
        }
        Ok(Self {
            neighbors,
            areas_m2: areas,
        })
    }

    /// Maximum diffusion number over land cells for `dt`.
    pub fn diffusion_number_max(&self, d: f64, dt: f64, ocean: &[bool]) -> f64 {
        let mut max = 0.0f64;
        for (i, list) in self.neighbors.iter().enumerate() {
            if ocean[i] {
                continue;
            }
            let sum: f64 = list.iter().map(|(_, c)| c).sum();
            max = max.max(dt * d * sum / self.areas_m2[i]);
        }
        max
    }
}

/// Bookkeeping closure state for lowering a cell: sediment first, then
/// bedrock, with the eroded volume accumulated.
struct Lowering<'a> {
    areas: &'a [f64],
    /// Bedrock column available at the start of the step (m), if bounded.
    bedrock_available: Option<&'a [f64]>,
    z: Vec<f64>,
    sed: Vec<f64>,
    bedrock_removed: Vec<f64>,
    sediment_removed: Vec<f64>,
    sediment_added: Vec<f64>,
    eroded_m3: f64,
    bedrock_eroded_m3: f64,
    sediment_eroded_m3: f64,
    deposited_m3: f64,
    exported_m3: f64,
    max_lowering: f64,
    bound_hits: u64,
}

impl Lowering<'_> {
    /// Material (m) that may still be removed from cell `i` now: the
    /// sediment present plus the unremoved part of the bedrock column.
    /// Unbounded when no column was supplied.
    fn removable(&self, i: usize) -> f64 {
        match self.bedrock_available {
            None => f64::INFINITY,
            Some(rock) => self.sed[i] + (rock[i] - self.bedrock_removed[i]).max(0.0),
        }
    }

    fn lower(&mut self, i: usize, amount: f64) {
        let from_sed = amount.min(self.sed[i]);
        let from_rock = amount - from_sed;
        self.sed[i] -= from_sed;
        self.sediment_removed[i] += from_sed;
        self.bedrock_removed[i] += from_rock;
        self.z[i] -= amount;
        self.eroded_m3 += amount * self.areas[i];
        self.sediment_eroded_m3 += from_sed * self.areas[i];
        self.bedrock_eroded_m3 += from_rock * self.areas[i];
        self.max_lowering = self.max_lowering.max(amount);
    }

    fn raise(&mut self, i: usize, amount: f64, ocean: bool) {
        self.sed[i] += amount;
        self.z[i] += amount;
        self.sediment_added[i] += amount;
        if ocean {
            self.exported_m3 += amount * self.areas[i];
        } else {
            self.deposited_m3 += amount * self.areas[i];
        }
    }

    /// `OncePerStepV1` fills: the surface may not end above `level` through
    /// rounding (the booked amount is unchanged; the clamp is <= 1 ulp).
    fn clamp_to(&mut self, i: usize, level: f64) {
        if self.z[i] > level {
            self.z[i] = level;
        }
    }
}

/// One erosion step over `dt_my`. `elevation_m` and `sediment_m` are the
/// physical fields BEFORE the step (not modified); `routing` was computed
/// on `elevation_m`; `sea_level_m` is the ocean water level (ocean cells
/// are base level: never lowered, receive exports into their sediment
/// column up to sea level, spread across the ocean body).
#[allow(clippy::too_many_arguments)]
pub fn erode_step(
    grid: &GridSpec,
    stencil: &DiffusionStencil,
    routing: &Routing,
    elevation_m: &[f64],
    sediment_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    dt_my: f64,
    params: &ErosionParams,
) -> Result<ErosionStep, ErosionError> {
    erode_step_bounded(
        grid,
        stencil,
        routing,
        elevation_m,
        sediment_m,
        ocean,
        sea_level_m,
        dt_my,
        params,
        None,
    )
}

/// `erode_step` with the per-step available-material bound:
/// `bedrock_available_m[c]` is the bedrock column (m, finite, >= 0) cell
/// `c` holds at the start of the step (the crust thickness, for the
/// preview driver). No cell loses more than that column plus the sediment
/// it holds when the removal happens. `None` is the unbounded legacy step.
#[allow(clippy::too_many_arguments)]
pub fn erode_step_bounded(
    grid: &GridSpec,
    stencil: &DiffusionStencil,
    routing: &Routing,
    elevation_m: &[f64],
    sediment_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    dt_my: f64,
    params: &ErosionParams,
    bedrock_available_m: Option<&[f64]>,
) -> Result<ErosionStep, ErosionError> {
    erode_step_bounded_layer(
        grid,
        stencil,
        routing,
        elevation_m,
        sediment_m,
        ocean,
        sea_level_m,
        dt_my,
        params,
        bedrock_available_m,
        OceanRemainderLayer::PerMouthV0,
    )
}

/// `erode_step_bounded` with the marine remainder layer selected
/// (`OceanRemainderLayer`); `PerMouthV0` is `erode_step_bounded` bit for bit.
#[allow(clippy::too_many_arguments)]
pub fn erode_step_bounded_layer(
    grid: &GridSpec,
    stencil: &DiffusionStencil,
    routing: &Routing,
    elevation_m: &[f64],
    sediment_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    dt_my: f64,
    params: &ErosionParams,
    bedrock_available_m: Option<&[f64]>,
    layer_law: OceanRemainderLayer,
) -> Result<ErosionStep, ErosionError> {
    params.validate()?;
    if !sea_level_m.is_finite() {
        return Err(ErosionError::BadParam {
            what: "sea_level_m",
            value: sea_level_m,
        });
    }
    let n = grid.cells();
    for (what, len) in [
        ("elevation", elevation_m.len()),
        ("sediment", sediment_m.len()),
        ("ocean", ocean.len()),
        ("routing.flow_to", routing.flow_to.len()),
        ("stencil", stencil.neighbors.len()),
    ] {
        if len != n {
            return Err(ErosionError::SizeMismatch {
                what,
                expected: n,
                found: len,
            });
        }
    }
    if !dt_my.is_finite() || dt_my <= 0.0 {
        return Err(ErosionError::BadParam {
            what: "dt_my",
            value: dt_my,
        });
    }
    if let Some(rock) = bedrock_available_m {
        if rock.len() != n {
            return Err(ErosionError::SizeMismatch {
                what: "bedrock_available_m",
                expected: n,
                found: rock.len(),
            });
        }
        for (i, &r) in rock.iter().enumerate() {
            if !r.is_finite() || r < 0.0 {
                return Err(ErosionError::NonFinite {
                    what: "bedrock_available_m",
                    value: r,
                    index: i,
                });
            }
        }
    }
    for (i, (&e, &s)) in elevation_m.iter().zip(sediment_m).enumerate() {
        if !e.is_finite() {
            return Err(ErosionError::NonFinite {
                what: "elevation_m",
                value: e,
                index: i,
            });
        }
        if !s.is_finite() || s < 0.0 {
            return Err(ErosionError::NonFinite {
                what: "sediment_m",
                value: s,
                index: i,
            });
        }
    }
    let v1 = layer_law == OceanRemainderLayer::OncePerStepV1;
    let (mut layer_mouths, mut layer_raises, mut layer_above) = (0u64, 0u64, 0u64);
    let mut deferred_remainder = 0.0f64;
    let areas = &stencil.areas_m2;
    let mut st = Lowering {
        areas,
        bedrock_available: bedrock_available_m,
        z: elevation_m.to_vec(),
        sed: sediment_m.to_vec(),
        bedrock_removed: vec![0.0; n],
        sediment_removed: vec![0.0; n],
        sediment_added: vec![0.0; n],
        eroded_m3: 0.0,
        bedrock_eroded_m3: 0.0,
        sediment_eroded_m3: 0.0,
        deposited_m3: 0.0,
        exported_m3: 0.0,
        max_lowering: 0.0,
        bound_hits: 0,
    };

    // ---- Stream power, implicit, receivers before donors ----------------
    // Sink of each cell's detached mass (memoised, receiver-first): the
    // first LAKE cell downstream (fill above the physical surface), else
    // the ocean cell or terminal floor the path ends in.
    let is_lake = |c: usize| !ocean[c] && routing.fill_level_m[c] > elevation_m[c] + 1.0e-9;
    let mut sink: Vec<Option<usize>> = vec![None; n];
    let order = &routing.receiver_first;
    let wet: Vec<bool> = (0..n).map(|i| ocean[i] || is_lake(i)).collect();
    let mut entry: Vec<Option<(usize, usize)>> = vec![None; n];
    let mut transfers = std::collections::BTreeMap::<(usize, usize, &str), f64>::new();
    for &i in order {
        sink[i] = match routing.flow_to[i] {
            None => Some(i),
            Some(ElementId(j)) => {
                let j = j as usize;
                if ocean[j] || is_lake(j) {
                    entry[i] = Some((i, j));
                    Some(j)
                } else {
                    entry[i] = entry[j];
                    Some(sink[j].expect("validated receiver-first order resolves every sink"))
                }
            }
        };
    }
    let mut pending_deposit: Vec<f64> = vec![0.0; n]; // volume arriving at sinks (m^3)
    for &i in order {
        if ocean[i] || is_lake(i) {
            continue; // base level, or submerged under a lake: no stream power
        }
        let Some(ElementId(j)) = routing.flow_to[i] else {
            continue; // terminal floor: no fluvial lowering
        };
        let j = j as usize;
        let base = if wet[j] {
            routing.fill_level_m[j]
        } else {
            st.z[j]
        };
        if st.z[i] <= base {
            continue; // not above its receiver on the physical surface
        }
        let (ii, ij) = grid
            .ij(i as u64)
            .map_err(|e| ErosionError::Grid(e.to_string()))?;
        let (ji, jj) = grid
            .ij(j as u64)
            .map_err(|e| ErosionError::Grid(e.to_string()))?;
        let len = grid
            .center_distance_m(ii, ij, ji, jj)
            .map_err(|e| ErosionError::Grid(e.to_string()))?
            .max(1.0);
        let f = params.k_per_my * routing.contributing_area_m2[i].powf(params.m) * dt_my / len;
        let z_new = (st.z[i] + f * base) / (1.0 + f);
        let mut amount = st.z[i] - z_new;
        if amount <= 0.0 {
            continue;
        }
        // Available-material bound (supply limit): never remove more than
        // the column holds. Inactive (infinite) without a supplied column.
        let removable = st.removable(i);
        if amount > removable {
            amount = removable;
            st.bound_hits += 1;
            if amount <= 0.0 {
                continue;
            }
        }
        st.lower(i, amount);
        let t = sink[i].expect("sink resolved in receiver-first order");
        pending_deposit[t] += amount * areas[i];
        let (from, to) = entry[i].expect("land flow must reach receiving water");
        *transfers.entry((from, to, "fluvial_delivery")).or_default() += amount * areas[i];
    }
    // Lay each sink's arriving volume across its water body (lake or
    // ocean) up to the water level; a full lake passes the remainder
    // downstream to the next lake or the ocean; the ocean body aggrades
    // uniformly; only a terminal floor with no outlet aggrades its basin.
    let mut deposit_aggraded_m3 = 0.0f64;
    let mut spread_bound_hits = 0u64;
    let ocean_cells: Vec<usize> = (0..n).filter(|&c| ocean[c]).collect();
    let ocean_area: f64 = ocean_cells.iter().map(|&c| areas[c]).sum();
    let water_level = |c: usize| -> f64 {
        if ocean[c] {
            sea_level_m
        } else {
            routing.fill_level_m[c]
        }
    };
    for t0 in 0..n {
        let mut remaining = pending_deposit[t0];
        if remaining <= 0.0 {
            continue;
        }
        let mut start = t0;
        let mut hops = 0usize;
        while remaining > 0.0 {
            let is_ocean = ocean[start];
            // Bounded breadth-first fill over the body containing `start`.
            let mut visited = vec![start];
            let mut seen: std::collections::BTreeSet<usize> = [start].into_iter().collect();
            let mut head = 0usize;
            while head < visited.len() && remaining > 0.0 {
                let c = visited[head];
                head += 1;
                let capacity = (water_level(c) - st.z[c]).max(0.0) * areas[c];
                if capacity > 0.0 {
                    let take = remaining.min(capacity);
                    if v1 {
                        // Fill exactly to the water level (no ulp overshoot).
                        let level = water_level(c);
                        let h = if take == capacity { level - st.z[c] } else { take / areas[c] };
                        st.raise(c, h, is_ocean);
                        st.clamp_to(c, level);
                    } else {
                        st.raise(c, take / areas[c], is_ocean);
                    }
                    remaining -= take;
                }
                if visited.len() >= MAX_SPREAD_CELLS {
                    continue; // stop expanding; drain what is already queued
                }
                let (ci, cj) = grid
                    .ij(c as u64)
                    .map_err(|e| ErosionError::Grid(e.to_string()))?;
                let mut next: Vec<usize> = grid
                    .neighbors(ci, cj)
                    .map_err(|e| ErosionError::Grid(e.to_string()))?
                    .into_iter()
                    .map(|(ti, tj, _)| ti * grid.nlon() + tj)
                    .filter(|&nb| routing.water_body_id[nb] == routing.water_body_id[start])
                    .collect();
                next.sort_unstable();
                for nb in next {
                    if visited.len() < MAX_SPREAD_CELLS && seen.insert(nb) {
                        visited.push(nb);
                    }
                }
            }
            if remaining <= 0.0 {
                break;
            }
            if visited.len() >= MAX_SPREAD_CELLS {
                spread_bound_hits += 1;
            }
            if is_ocean {
                // Search-bound overflow is a declared global layer over this marine reservoir.
                layer_mouths += 1;
                if v1 {
                    deferred_remainder += remaining; // laid once, after every mouth
                } else {
                    let layer = remaining / ocean_area;
                    for &c in &ocean_cells {
                        if st.z[c] <= sea_level_m && st.z[c] + layer > sea_level_m {
                            layer_above += 1;
                        }
                        st.raise(c, layer, true);
                    }
                    layer_raises += ocean_cells.len() as u64;
                }
                deposit_aggraded_m3 += remaining;
                remaining = 0.0;
                break;
            }
            // A full lake: carry the remainder downstream to the next sink.
            let mut cur = start;
            let mut next_start: Option<usize> = None;
            loop {
                match routing.flow_to[cur] {
                    None => break,
                    Some(ElementId(j)) => {
                        let j = j as usize;
                        if ocean[j] || (is_lake(j) && water_level(j) > st.z[j]) {
                            next_start = Some(j);
                            break;
                        }
                        cur = j;
                        hops += 1;
                        if hops > n {
                            return Err(ErosionError::Grid("cycle during lake overflow".into()));
                        }
                    }
                }
            }
            match next_start {
                Some(s) => {
                    if routing.water_body_id[s] != routing.water_body_id[start] {
                        *transfers.entry((cur, s, "lake_overflow")).or_default() += remaining;
                    }
                    start = s;
                }
                None => {
                    return Err(ErosionError::Grid(
                        "lake overflow has no downstream receiving body".into(),
                    ))
                }
            }
        }
    }
    if v1 && deferred_remainder > 0.0 {
        let (cells, _rounds, overfull) =
            capped_ocean_layer(&st.z, areas, &ocean_cells, sea_level_m, deferred_remainder);
        layer_raises += cells.len() as u64;
        for (c, h, full) in cells {
            st.raise(c, h, true);
            if full {
                st.clamp_to(c, sea_level_m);
            } else if overfull {
                layer_above += 1;
            }
        }
    }

    // ---- Hillslope diffusion, explicit flux form, substepped ------------
    let d = params.diffusivity_m2_per_my;
    let nu_max = stencil.diffusion_number_max(d, dt_my, &wet);
    let substeps = if d == 0.0 || nu_max == 0.0 {
        1
    } else {
        ((nu_max / params.stability_limit).ceil() as u32).max(1)
    };
    let sub_dt = dt_my / substeps as f64;
    if d > 0.0 {
        let mut flux_dz = vec![0.0f64; n]; // m per substep, signed
        // Two-cell moves of one substep, in evaluation order: (from, to,
        // signed volume from -> to in m^3, delivered into water).
        let mut moves: Vec<(usize, usize, f64, bool)> = Vec::new();
        let mut outgoing = vec![0.0f64; n];
        for _ in 0..substeps {
            for v in flux_dz.iter_mut() {
                *v = 0.0;
            }
            moves.clear();
            for i in 0..n {
                if wet[i] {
                    continue; // this is terrestrial diffusion, not submerged-bed erosion
                }
                for &(j, coef) in &stencil.neighbors[i] {
                    // Volume moved i -> j (m^3): D [m^2/My] * (w/d) * dz [m] * dt [My].
                    let receiver = if wet[j] {
                        routing.fill_level_m[j]
                    } else {
                        st.z[j]
                    };
                    let flux = d * coef * (st.z[i] - receiver) * sub_dt;
                    if wet[j] {
                        // One-way sink: only downhill flux into the ocean.
                        if flux > 0.0 {
                            moves.push((i, j, flux, true));
                        }
                        continue;
                    }
                    if j < i {
                        continue; // each land pair handled once, from the lower id
                    }
                    moves.push((i, j, flux, false));
                }
            }
            // Available-material bound: scale each donor's outgoing volume
            // to what it holds. Without a supplied column every scale is
            // exactly 1 and the moves apply in the legacy order.
            let mut scale: Option<Vec<f64>> = None;
            if st.bedrock_available.is_some() {
                for v in outgoing.iter_mut() {
                    *v = 0.0;
                }
                for &(i, j, flux, _) in &moves {
                    if flux > 0.0 {
                        outgoing[i] += flux;
                    } else {
                        outgoing[j] -= flux;
                    }
                }
                let mut sc = vec![1.0f64; n];
                for c in 0..n {
                    if outgoing[c] > 0.0 {
                        let avail = st.removable(c) * areas[c];
                        if outgoing[c] > avail {
                            sc[c] = (avail / outgoing[c]).max(0.0);
                            st.bound_hits += 1;
                        }
                    }
                }
                scale = Some(sc);
            }
            for &(i, j, flux, into_water) in &moves {
                let flux = match &scale {
                    Some(sc) => {
                        let donor = if flux > 0.0 { i } else { j };
                        if sc[donor] < 1.0 {
                            flux * sc[donor]
                        } else {
                            flux
                        }
                    }
                    None => flux,
                };
                flux_dz[i] -= flux / areas[i];
                flux_dz[j] += flux / areas[j];
                if into_water {
                    *transfers.entry((i, j, "diffusive_delivery")).or_default() += flux;
                }
            }
            for i in 0..n {
                let dz = flux_dz[i];
                if dz == 0.0 {
                    continue;
                }
                if dz < 0.0 {
                    st.lower(i, -dz);
                } else {
                    st.raise(i, dz, ocean[i]);
                }
            }
        }
    }

    let dz: Vec<f64> = st.z.iter().zip(elevation_m).map(|(a, b)| a - b).collect();
    for (i, v) in dz.iter().enumerate() {
        if !v.is_finite() {
            return Err(ErosionError::NonFinite {
                what: "dz_m",
                value: *v,
                index: i,
            });
        }
    }
    let transfers = transfers
        .into_iter()
        .map(
            |((i, j, kind), volume)| pf_contracts::domain::WaterTransferV0 {
                from_element: ElementId(i as u64),
                to_element: ElementId(j as u64),
                receiving_body: routing.water_body_id[j],
                kind: kind.into(),
                sediment_m3: volume,
                runoff_m3_s: 0.0,
            },
        )
        .collect();
    Ok(ErosionStep {
        transfers,
        spread_bound_hits,
        dz_m: dz,
        bedrock_removed_m: st.bedrock_removed,
        sediment_removed_m: st.sediment_removed,
        sediment_added_m: st.sediment_added,
        eroded_m3: st.eroded_m3,
        deposited_m3: st.deposited_m3,
        exported_m3: st.exported_m3,
        diffusion_number_max: nu_max,
        diffusion_substeps: substeps,
        max_lowering_rate_m_per_my: st.max_lowering / dt_my,
        deposit_aggraded_m3,
        bedrock_eroded_m3: st.bedrock_eroded_m3,
        sediment_eroded_m3: st.sediment_eroded_m3,
        bedrock_bound_hits: st.bound_hits,
        deposition: None,
        ocean_layer_mouths: layer_mouths,
        ocean_layer_cell_raises: layer_raises,
        ocean_layer_cells_above_sea: layer_above,
    })
}

// ===========================================================================
// Lane SD (2026-10-12; design: "yes please get them working and decent
// simulations please"): erosion-deposition on land, law id
// `erosion-deposition-v1`, and compensated standing-water placement, law id
// `shelf-wave-base-v1` (see the note in `sediment.rs`). Both are reached only
// through `erode_deposit_step` with a non-legacy `SurfaceLaws`. With
// `SurfaceLaws::legacy()` that function IS `erode_step_bounded`, so the legacy
// path above is byte-identical.
//
// EROSION-DEPOSITION (Yuan et al. 2019, JGR Earth Surface 124, "A new
// efficient method to solve the stream power law model taking into account
// sediment deposition"; Davy and Lague 2009, JGR 114, the xi-q model; both
// recalled, UNVERIFIED): on every land cell
//     dh/dt = -K A^m S + G Qs / (p A),
// where Qs is the sediment flux leaving the cell (m^3 per unit time: the
// net volume lost by its whole upstream area including itself), p is the
// relative precipitation (1: uniform runoff, the same geometric area the
// stream-power term uses) and G is a dimensionless deposition coefficient
// (`DEPOSITION_COEFFICIENT_G`, declared before any run). Deposition is
// largest where the flux per unit drainage area is high and the slope is
// low: mountain-front fans, forelands, floodplains and lowlands. Those
// cells store sediment that the detachment-limited law sent straight to the
// sea.
// DISCRETISATION (ours; implicit like Yuan et al., not their algorithm).
// Over a step dt, each land cell i with receiver r satisfies
//     h_i (1 + f + G a/(pA)) = h0_i (1 + G a/(pA)) + f h_r + G Qin_i/(pA)
//     Qs_i = Qin_i + (h0_i - h_i) a_i
// Here f = K A^m dt / L, and stream power applies only while the solved
// h_i stands above h_r (with G = 0, h0_i > h_r, the legacy rule).
// a_i is the cell area, A_i the contributing area, and Qin_i the volume
// arriving from donors over the step. For a fixed active set (which cells
// sit above their receiver, and which are held by lane E's
// available-material bound) the system is linear on the drainage tree. It
// is solved DIRECTLY: a donor-first pass eliminates each inflow as an affine
// function of the cell's own level, then a receiver-first pass
// back-substitutes the levels. An outer loop, starting from the legacy
// solve's active set, grows the set (it never shrinks within a step, so
// the loop terminates) until it stops changing, capped at `ED_MAX_PASSES`;
// the pass count and convergence are reported. MASS IS
// EXACT whatever the outer loop did: a final donor-first pass routes
// Qs = Qin + net loss. With G = 0 the update is the legacy implicit stream
// power term for term (tested).
// LIMITS (declared): local Airy only, so there is no flexural foreland
// depression and "foreland" storage is the mountain-front alluvial wedge;
// lake overflow passes to the next standing water without depositing on the
// river between; no grain sizes.
// ===========================================================================

use crate::sediment::{
    DepositionReport, MarineClass, ShelfParams, ShelfPlacer, SHELF_WAVE_BASE_LAW_ID,
};

/// Law ids (strings recorded in the step report; the typed law set lives in
/// pf-contracts laws.rs, owned by another lane).
pub const LEGACY_EROSION_LAW_ID: &str = "stream-power-n1-linear-diffusion-v0";
pub const EROSION_DEPOSITION_LAW_ID: &str = "erosion-deposition-v1";
pub const LEGACY_SEDIMENT_LAW_ID: &str = "local-fill-global-layer-v0";
/// Dimensionless deposition coefficient G of `erosion-deposition-v1`.
/// Declared BEFORE any run (lane SD FALSIFIERS.md) and never tuned: Yuan et
/// al. 2019 work with G of order 1 (recalled, UNVERIFIED). Project
/// assumption awaiting review.
pub const DEPOSITION_COEFFICIENT_G: f64 = 1.0;
/// Solver control (numerical, not physical): cap on the outer active-set
/// passes of the direct tree solve.
pub const ED_MAX_PASSES: u32 = 64;

/// Parameters of `erosion-deposition-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ErosionDepositionParams {
    pub deposition_coefficient_g: f64,
    pub precipitation_factor: f64,
    /// Cap on the outer active-set passes.
    pub max_sweeps: u32,
}

impl ErosionDepositionParams {
    pub fn v1() -> Self {
        Self {
            deposition_coefficient_g: DEPOSITION_COEFFICIENT_G,
            precipitation_factor: 1.0,
            max_sweeps: ED_MAX_PASSES,
        }
    }

    pub fn validate(&self) -> Result<(), ErosionError> {
        if !self.deposition_coefficient_g.is_finite() || self.deposition_coefficient_g < 0.0 {
            return Err(ErosionError::BadParam {
                what: "deposition_coefficient_g",
                value: self.deposition_coefficient_g,
            });
        }
        if !self.precipitation_factor.is_finite() || self.precipitation_factor <= 0.0 {
            return Err(ErosionError::BadParam {
                what: "precipitation_factor",
                value: self.precipitation_factor,
            });
        }
        if self.max_sweeps == 0 {
            return Err(ErosionError::BadParam {
                what: "max_sweeps",
                value: 0.0,
            });
        }
        Ok(())
    }
}

/// The surface laws a step runs under. `None` = the legacy law for that
/// part (detachment-limited land; local fill + global layer placement).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SurfaceLaws {
    pub erosion_deposition: Option<ErosionDepositionParams>,
    pub shelf: Option<ShelfParams>,
}

impl SurfaceLaws {
    pub fn legacy() -> Self {
        Self::default()
    }

    pub fn is_legacy(&self) -> bool {
        self.erosion_deposition.is_none() && self.shelf.is_none()
    }

    pub fn erosion_law_id(&self) -> &'static str {
        if self.erosion_deposition.is_some() {
            EROSION_DEPOSITION_LAW_ID
        } else {
            LEGACY_EROSION_LAW_ID
        }
    }

    pub fn sediment_law_id(&self) -> &'static str {
        if self.shelf.is_some() {
            SHELF_WAVE_BASE_LAW_ID
        } else {
            LEGACY_SEDIMENT_LAW_ID
        }
    }
}

/// Where a deposit goes in the partition (and whether it is "exported").
#[derive(Debug, Clone, Copy, PartialEq)]
enum Sink {
    LandFluvial,
    LandHillslope,
    LandFloor,
    Lake,
    Marine(MarineClass),
    Remainder,
}

/// Raise cell `c` by `h` m and book it in the partition (one place, so every
/// raise has exactly one partition entry).
fn deposit(st: &mut Lowering<'_>, rep: &mut DepositionReport, c: usize, h: f64, sink: Sink) {
    if h <= 0.0 {
        return;
    }
    let marine = matches!(sink, Sink::Marine(_) | Sink::Remainder);
    st.raise(c, h, marine);
    let v = h * st.areas[c];
    match sink {
        Sink::LandFluvial => rep.land_fluvial_m3 += v,
        Sink::LandHillslope => rep.land_hillslope_m3 += v,
        Sink::LandFloor => rep.land_floor_m3 += v,
        Sink::Lake => rep.lake_m3 += v,
        Sink::Marine(MarineClass::Shelf) => rep.marine_shelf_m3 += v,
        Sink::Marine(MarineClass::Slope) => rep.marine_slope_m3 += v,
        Sink::Marine(MarineClass::Deep) => rep.marine_deep_m3 += v,
        Sink::Remainder => rep.global_remainder_m3 += v,
    }
}

/// One erosion step under the selected surface laws. Inputs and the
/// returned `ErosionStep` mean what they mean for `erode_step_bounded`, with
/// two differences when a law is selected. `dz_m` is the COLUMN change (the
/// caller's isostasy turns it into a surface change), and `deposition`
/// carries the partition. With `SurfaceLaws::legacy()` this is exactly
/// `erode_step_bounded`.
#[allow(clippy::too_many_arguments)]
pub fn erode_deposit_step(
    grid: &GridSpec,
    stencil: &DiffusionStencil,
    routing: &Routing,
    elevation_m: &[f64],
    sediment_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    dt_my: f64,
    params: &ErosionParams,
    laws: &SurfaceLaws,
    bedrock_available_m: Option<&[f64]>,
) -> Result<ErosionStep, ErosionError> {
    if laws.is_legacy() {
        return erode_step_bounded(
            grid,
            stencil,
            routing,
            elevation_m,
            sediment_m,
            ocean,
            sea_level_m,
            dt_my,
            params,
            bedrock_available_m,
        );
    }
    params.validate()?;
    let ed = match laws.erosion_deposition {
        Some(p) => {
            p.validate()?;
            p
        }
        None => ErosionDepositionParams {
            deposition_coefficient_g: 0.0,
            ..ErosionDepositionParams::v1()
        },
    };
    if let Some(s) = &laws.shelf {
        s.validate()
            .map_err(|e| ErosionError::Law(format!("{SHELF_WAVE_BASE_LAW_ID}: {e}")))?;
    }
    if !sea_level_m.is_finite() {
        return Err(ErosionError::BadParam {
            what: "sea_level_m",
            value: sea_level_m,
        });
    }
    let n = grid.cells();
    for (what, len) in [
        ("elevation", elevation_m.len()),
        ("sediment", sediment_m.len()),
        ("ocean", ocean.len()),
        ("routing.flow_to", routing.flow_to.len()),
        ("stencil", stencil.neighbors.len()),
    ] {
        if len != n {
            return Err(ErosionError::SizeMismatch {
                what,
                expected: n,
                found: len,
            });
        }
    }
    if !dt_my.is_finite() || dt_my <= 0.0 {
        return Err(ErosionError::BadParam {
            what: "dt_my",
            value: dt_my,
        });
    }
    if let Some(rock) = bedrock_available_m {
        if rock.len() != n {
            return Err(ErosionError::SizeMismatch {
                what: "bedrock_available_m",
                expected: n,
                found: rock.len(),
            });
        }
        for (i, &r) in rock.iter().enumerate() {
            if !r.is_finite() || r < 0.0 {
                return Err(ErosionError::NonFinite {
                    what: "bedrock_available_m",
                    value: r,
                    index: i,
                });
            }
        }
    }
    for (i, (&e, &s)) in elevation_m.iter().zip(sediment_m).enumerate() {
        if !e.is_finite() {
            return Err(ErosionError::NonFinite {
                what: "elevation_m",
                value: e,
                index: i,
            });
        }
        if !s.is_finite() || s < 0.0 {
            return Err(ErosionError::NonFinite {
                what: "sediment_m",
                value: s,
                index: i,
            });
        }
    }
    let areas = &stencil.areas_m2;
    let mut st = Lowering {
        areas,
        bedrock_available: bedrock_available_m,
        z: elevation_m.to_vec(),
        sed: sediment_m.to_vec(),
        bedrock_removed: vec![0.0; n],
        sediment_removed: vec![0.0; n],
        sediment_added: vec![0.0; n],
        eroded_m3: 0.0,
        bedrock_eroded_m3: 0.0,
        sediment_eroded_m3: 0.0,
        deposited_m3: 0.0,
        exported_m3: 0.0,
        max_lowering: 0.0,
        bound_hits: 0,
    };
    let mut rep = DepositionReport {
        erosion_law: laws.erosion_law_id(),
        sediment_law: laws.sediment_law_id(),
        ..Default::default()
    };
    let mut transfers = std::collections::BTreeMap::<(usize, usize, &str), f64>::new();

    // ---- Fluvial erosion-deposition on land ------------------------------
    let is_lake = |c: usize| !ocean[c] && routing.fill_level_m[c] > elevation_m[c] + 1.0e-9;
    let wet: Vec<bool> = (0..n).map(|i| ocean[i] || is_lake(i)).collect();
    let order = &routing.receiver_first;
    const NONE: usize = usize::MAX;
    let mut recv = vec![NONE; n];
    let mut fco = vec![0.0f64; n];
    for i in 0..n {
        if wet[i] {
            continue;
        }
        if let Some(ElementId(j)) = routing.flow_to[i] {
            let j = j as usize;
            let (ii, ij) = grid
                .ij(i as u64)
                .map_err(|e| ErosionError::Grid(e.to_string()))?;
            let (ji, jj) = grid
                .ij(j as u64)
                .map_err(|e| ErosionError::Grid(e.to_string()))?;
            let len = grid
                .center_distance_m(ii, ij, ji, jj)
                .map_err(|e| ErosionError::Grid(e.to_string()))?
                .max(1.0);
            recv[i] = j;
            fco[i] = params.k_per_my * routing.contributing_area_m2[i].powf(params.m) * dt_my / len;
        }
    }
    let removable0: Vec<f64> = match bedrock_available_m {
        None => vec![f64::INFINITY; n],
        Some(rock) => (0..n).map(|i| sediment_m[i] + rock[i]).collect(),
    };
    let h0 = elevation_m;
    let base_of = |j: usize, h: &[f64]| -> f64 {
        if wet[j] {
            routing.fill_level_m[j]
        } else {
            h[j]
        }
    };
    // The coupled system is LINEAR on the drainage tree once the active set
    // is fixed. The active set is which cells stand above their receiver
    // (stream power on) and which are held by the available-material bound.
    // Donor-first, each cell's inflow is eliminated as an affine function of
    // its own level, Qin_i = A_i + B_i h_i (B_i <= 0: a higher receiver
    // erodes its donors less). Receiver-first back-substitution then gives
    // h_i = (P_i + f_i h_r) / D_i with D_i = 1 + f_i + G a_i/(pA_i) -
    // gamma_i B_i >= 1. The outer loop updates the active set until it stops
    // changing; it starts from the legacy implicit solve's active set. A
    // final donor-first pass routes Qs = Qin + net loss, so mass is exact
    // whatever the outer loop did. With G = 0 this is the legacy update
    // (P = h0, D = 1 + f) term for term.
    let g = ed.deposition_coefficient_g;
    let p = ed.precipitation_factor;
    // gamma_i = G / (p A_i) (1/m^2).
    let gamma: Vec<f64> = (0..n)
        .map(|i| {
            if g > 0.0 {
                g / (p * routing.contributing_area_m2[i].max(areas[i]))
            } else {
                0.0
            }
        })
        .collect();
    let mut on = vec![false; n]; // stream power active
    let mut held = vec![false; n]; // available-material bound active
    {
        // Active set of the legacy detachment-limited implicit solve.
        let mut hs = h0.to_vec();
        for &i in order {
            let j = recv[i];
            if j == NONE {
                continue;
            }
            let base = base_of(j, &hs);
            if hs[i] <= base {
                continue;
            }
            on[i] = true;
            let f = fco[i];
            let z_new = (hs[i] + f * base) / (1.0 + f);
            let amount = hs[i] - z_new;
            if amount > removable0[i] {
                held[i] = true;
                hs[i] = h0[i] - removable0[i];
            } else if amount > 0.0 {
                hs[i] -= amount;
            }
        }
    }
    let mut h = h0.to_vec();
    let mut h_prev = h0.to_vec();
    let (mut ca, mut cb) = (vec![0.0f64; n], vec![0.0f64; n]);
    let (mut cp, mut cd) = (vec![0.0f64; n], vec![1.0f64; n]);
    let mut passes = 0u32;
    let mut converged = false;
    let mut change = 0.0f64;
    while passes < ed.max_sweeps {
        passes += 1;
        ca.iter_mut().for_each(|v| *v = 0.0);
        cb.iter_mut().for_each(|v| *v = 0.0);
        for &i in order.iter().rev() {
            let j = recv[i];
            if j == NONE {
                continue;
            }
            let a = areas[i];
            let (alpha, beta);
            if held[i] {
                cp[i] = h0[i] - removable0[i];
                cd[i] = 1.0;
                let qin_i = ca[i] + cb[i] * cp[i];
                alpha = qin_i + removable0[i] * a;
                beta = 0.0;
            } else {
                let f = if on[i] { fco[i] } else { 0.0 };
                let ga = gamma[i] * a;
                cd[i] = 1.0 + f + ga - gamma[i] * cb[i];
                cp[i] = h0[i] * (1.0 + ga) + gamma[i] * ca[i];
                alpha = ca[i] + h0[i] * a + (cb[i] - a) * cp[i] / cd[i];
                beta = (cb[i] - a) * f / cd[i];
            }
            if !wet[j] {
                ca[j] += alpha;
                cb[j] += beta;
            }
        }
        for &i in order {
            let j = recv[i];
            if j == NONE {
                continue;
            }
            h[i] = if held[i] {
                cp[i]
            } else {
                let f = if on[i] { fco[i] } else { 0.0 };
                (cp[i] + f * base_of(j, &h)) / cd[i]
            };
        }
        let mut flips = 0u64;
        change = 0.0;
        for &i in order {
            let j = recv[i];
            if j == NONE {
                continue;
            }
            change = change.max((h[i] - h_prev[i]).abs());
            if held[i] {
                continue;
            }
            // Stream power acts while the SOLVED level stands above the
            // receiver. h_i - h_r has the sign of
            // (h0_i - h_r)(1 + G a/(pA)) + G Qin_i/(pA) whether or not f
            // is on, so this test is consistent in both states (with G = 0
            // it is the legacy test h0_i > h_r).
            let qin_i = ca[i] + cb[i] * h[i];
            let ga = gamma[i] * areas[i];
            let above = (h0[i] - base_of(j, &h)) * (1.0 + ga) + gamma[i] * qin_i > 0.0;
            // The active set only GROWS within a step, so the loop
            // terminates. A two-way set cycled on flat plateaus: 2 of 50
            // steps hit the cap on seed 240907 L1, lane SD probe. A cell
            // that stays on after sinking below its receiver is pulled up
            // TOWARD the receiver level (filled from the passing flux, mass
            // exact), never above it.
            if above && !on[i] {
                on[i] = true;
                flips += 1;
            }
            if h0[i] - h[i] > removable0[i] {
                held[i] = true;
                flips += 1;
            }
        }
        h_prev.copy_from_slice(&h);
        if flips == 0 {
            converged = true;
            break;
        }
    }
    rep.sweeps = passes;
    rep.converged = converged;
    rep.last_sweep_max_change_m = change;
    // Final mass-exact routing of the solved levels.
    let mut qs = vec![0.0f64; n];
    let mut qin = vec![0.0f64; n];
    let fluvial_hits = held.iter().filter(|&&b| b).count() as u64;
    for &i in order.iter().rev() {
        let j = recv[i];
        if j == NONE {
            continue;
        }
        let mut q = qin[i] + (h0[i] - h[i]) * areas[i];
        if q < 0.0 {
            // Rounding only (Qs >= 0 analytically): keep the volume exact.
            h[i] = h0[i] + qin[i] / areas[i];
            q = 0.0;
        }
        qs[i] = q;
        if !wet[j] {
            qin[j] += q;
        }
    }
    st.bound_hits += fluvial_hits;
    // Apply the fluvial result and collect what reaches standing water.
    let mut pending = vec![0.0f64; n];
    for &i in order {
        let j = recv[i];
        if j == NONE {
            if !wet[i] && qin[i] > 0.0 {
                // A terminal floor with no outlet keeps what reaches it.
                deposit(&mut st, &mut rep, i, qin[i] / areas[i], Sink::LandFloor);
            }
            continue;
        }
        let d = h[i] - h0[i];
        if d < 0.0 {
            st.lower(i, -d);
        } else if d > 0.0 {
            deposit(&mut st, &mut rep, i, d, Sink::LandFluvial);
        }
        let gross_dep = gamma[i] * qs[i] * areas[i];
        rep.fluvial_gross_deposition_m3 += gross_dep;
        rep.fluvial_gross_erosion_m3 += gross_dep + (h0[i] - h[i]) * areas[i];
        if wet[j] && qs[i] > 0.0 {
            pending[j] += qs[i];
            *transfers.entry((i, j, "fluvial_delivery")).or_default() += qs[i];
        }
    }

    // ---- Standing-water placement -----------------------------------------
    let ocean_cells: Vec<usize> = (0..n).filter(|&c| ocean[c]).collect();
    let ocean_area: f64 = ocean_cells.iter().map(|&c| areas[c]).sum();
    let mut surf = elevation_m.to_vec();
    let mut placer = ShelfPlacer::new(n);
    let mut adjacency: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    if laws.shelf.is_some() {
        for &c in &ocean_cells {
            let (ci, cj) = grid
                .ij(c as u64)
                .map_err(|e| ErosionError::Grid(e.to_string()))?;
            let mut list: Vec<(usize, f64)> = grid
                .neighbors(ci, cj)
                .map_err(|e| ErosionError::Grid(e.to_string()))?
                .into_iter()
                .map(|(ti, tj, nb)| (ti * grid.nlon() + tj, nb.distance_m.max(1.0)))
                .filter(|&(nb, _)| ocean[nb])
                .collect();
            list.sort_unstable_by_key(|&(nb, _)| nb);
            adjacency[c] = list;
        }
    }
    let mut remainder_total = 0.0f64;
    let mut place = |st: &mut Lowering<'_>,
                     rep: &mut DepositionReport,
                     transfers: &mut std::collections::BTreeMap<(usize, usize, &str), f64>,
                     surf: &mut Vec<f64>,
                     pending: &[f64]|
     -> Result<(), ErosionError> {
        match &laws.shelf {
            Some(shelf) => place_compensated(
                grid,
                routing,
                elevation_m,
                ocean,
                sea_level_m,
                shelf,
                &adjacency,
                &mut placer,
                st,
                rep,
                transfers,
                surf,
                pending,
                &mut remainder_total,
            ),
            None => place_legacy(
                grid,
                routing,
                elevation_m,
                ocean,
                sea_level_m,
                &ocean_cells,
                ocean_area,
                st,
                rep,
                transfers,
                pending,
                &mut LayerTally::default(),
            ),
        }
    };
    place(&mut st, &mut rep, &mut transfers, &mut surf, &pending)?;

    // ---- Hillslope diffusion (legacy scheme; under the shelf law the flux
    // into standing water is placed like river delivery) -----------------
    let d = params.diffusivity_m2_per_my;
    let nu_max = stencil.diffusion_number_max(d, dt_my, &wet);
    let substeps = if d == 0.0 || nu_max == 0.0 {
        1
    } else {
        ((nu_max / params.stability_limit).ceil() as u32).max(1)
    };
    let sub_dt = dt_my / substeps as f64;
    let shelf_on = laws.shelf.is_some();
    let mut pending_diffusive = vec![0.0f64; n];
    if d > 0.0 {
        let mut flux_dz = vec![0.0f64; n];
        let mut moves: Vec<(usize, usize, f64, bool)> = Vec::new();
        let mut outgoing = vec![0.0f64; n];
        for _ in 0..substeps {
            flux_dz.iter_mut().for_each(|v| *v = 0.0);
            moves.clear();
            for i in 0..n {
                if wet[i] {
                    continue;
                }
                for &(j, coef) in &stencil.neighbors[i] {
                    let receiver = if wet[j] {
                        routing.fill_level_m[j]
                    } else {
                        st.z[j]
                    };
                    let flux = d * coef * (st.z[i] - receiver) * sub_dt;
                    if wet[j] {
                        if flux > 0.0 {
                            moves.push((i, j, flux, true));
                        }
                        continue;
                    }
                    if j < i {
                        continue;
                    }
                    moves.push((i, j, flux, false));
                }
            }
            let mut scale: Option<Vec<f64>> = None;
            if st.bedrock_available.is_some() {
                outgoing.iter_mut().for_each(|v| *v = 0.0);
                for &(i, j, flux, _) in &moves {
                    if flux > 0.0 {
                        outgoing[i] += flux;
                    } else {
                        outgoing[j] -= flux;
                    }
                }
                let mut sc = vec![1.0f64; n];
                for c in 0..n {
                    if outgoing[c] > 0.0 {
                        let avail = st.removable(c) * areas[c];
                        if outgoing[c] > avail {
                            sc[c] = (avail / outgoing[c]).max(0.0);
                            st.bound_hits += 1;
                        }
                    }
                }
                scale = Some(sc);
            }
            for &(i, j, flux, into_water) in &moves {
                let flux = match &scale {
                    Some(sc) => {
                        let donor = if flux > 0.0 { i } else { j };
                        if sc[donor] < 1.0 {
                            flux * sc[donor]
                        } else {
                            flux
                        }
                    }
                    None => flux,
                };
                flux_dz[i] -= flux / areas[i];
                if into_water {
                    *transfers.entry((i, j, "diffusive_delivery")).or_default() += flux;
                    if shelf_on {
                        pending_diffusive[j] += flux;
                        continue;
                    }
                }
                flux_dz[j] += flux / areas[j];
            }
            for i in 0..n {
                let dz = flux_dz[i];
                if dz == 0.0 {
                    continue;
                }
                if dz < 0.0 {
                    st.lower(i, -dz);
                } else if ocean[i] {
                    let class = MarineClass::of_depth(sea_level_m - st.z[i]);
                    deposit(&mut st, &mut rep, i, dz, Sink::Marine(class));
                } else if wet[i] {
                    deposit(&mut st, &mut rep, i, dz, Sink::Lake);
                } else {
                    deposit(&mut st, &mut rep, i, dz, Sink::LandHillslope);
                }
            }
        }
    }
    if shelf_on && pending_diffusive.iter().any(|&v| v > 0.0) {
        place(&mut st, &mut rep, &mut transfers, &mut surf, &pending_diffusive)?;
    }
    drop(place);
    // The measured global remainder (shelf law): one uniform layer.
    if remainder_total > 0.0 {
        let layer = remainder_total / ocean_area;
        for &c in &ocean_cells {
            deposit(&mut st, &mut rep, c, layer, Sink::Remainder);
        }
    }

    let dz: Vec<f64> = st.z.iter().zip(elevation_m).map(|(a, b)| a - b).collect();
    for (i, v) in dz.iter().enumerate() {
        if !v.is_finite() {
            return Err(ErosionError::NonFinite {
                what: "dz_m",
                value: *v,
                index: i,
            });
        }
    }
    let transfers = transfers
        .into_iter()
        .map(
            |((i, j, kind), volume)| pf_contracts::domain::WaterTransferV0 {
                from_element: ElementId(i as u64),
                to_element: ElementId(j as u64),
                receiving_body: routing.water_body_id[j],
                kind: kind.into(),
                sediment_m3: volume,
                runoff_m3_s: 0.0,
            },
        )
        .collect();
    Ok(ErosionStep {
        transfers,
        spread_bound_hits: rep.lake_spread_bound_hits + rep.marine_search_bound_hits,
        dz_m: dz,
        bedrock_removed_m: st.bedrock_removed,
        sediment_removed_m: st.sediment_removed,
        sediment_added_m: st.sediment_added,
        eroded_m3: st.eroded_m3,
        deposited_m3: st.deposited_m3,
        exported_m3: st.exported_m3,
        diffusion_number_max: nu_max,
        diffusion_substeps: substeps,
        max_lowering_rate_m_per_my: st.max_lowering / dt_my,
        deposit_aggraded_m3: rep.global_remainder_m3,
        bedrock_eroded_m3: st.bedrock_eroded_m3,
        sediment_eroded_m3: st.sediment_eroded_m3,
        bedrock_bound_hits: st.bound_hits,
        deposition: Some(rep),
        ocean_layer_mouths: 0,
        ocean_layer_cell_raises: 0,
        ocean_layer_cells_above_sea: 0,
    })
}

/// `shelf-wave-base-v1` placement of the volumes pending at sinks: lakes
/// are filled to their spill level with COMPENSATED thickness (bounded
/// breadth-first spread as in the legacy step; a full lake passes the rest
/// downstream), and ocean deliveries are placed by `ShelfPlacer`. What the
/// placer cannot place is added to `remainder_total`. The caller lays that
/// once, as the declared global remainder.
#[allow(clippy::too_many_arguments)]
fn place_compensated(
    grid: &GridSpec,
    routing: &Routing,
    elevation_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    shelf: &ShelfParams,
    adjacency: &[Vec<(usize, f64)>],
    placer: &mut ShelfPlacer,
    st: &mut Lowering<'_>,
    rep: &mut DepositionReport,
    transfers: &mut std::collections::BTreeMap<(usize, usize, &str), f64>,
    surf: &mut Vec<f64>,
    pending: &[f64],
    remainder_total: &mut f64,
) -> Result<(), ErosionError> {
    let n = pending.len();
    let areas = st.areas;
    let col = &shelf.column;
    let is_lake = |c: usize| !ocean[c] && routing.fill_level_m[c] > elevation_m[c] + 1.0e-9;
    let mut marine_pending = vec![0.0f64; n];
    for t0 in 0..n {
        let mut remaining = pending[t0];
        if remaining <= 0.0 {
            continue;
        }
        if ocean[t0] {
            marine_pending[t0] += remaining;
            continue;
        }
        let mut start = t0;
        let mut hops = 0usize;
        while remaining > 0.0 {
            if ocean[start] {
                marine_pending[start] += remaining;
                break;
            }
            // Bounded breadth-first compensated fill of the lake body.
            let mut visited = vec![start];
            let mut seen: std::collections::BTreeSet<usize> = [start].into_iter().collect();
            let mut head = 0usize;
            while head < visited.len() && remaining > 0.0 {
                let c = visited[head];
                head += 1;
                let level = routing.fill_level_m[c];
                if surf[c] < level {
                    let need = col.thickness_to_raise(surf[c], level, sea_level_m) * areas[c];
                    if remaining < need {
                        let h = remaining / areas[c];
                        surf[c] = col.surface_after(surf[c], h, sea_level_m);
                        deposit(st, rep, c, h, Sink::Lake);
                        remaining = 0.0;
                    } else {
                        surf[c] = level;
                        deposit(st, rep, c, need / areas[c], Sink::Lake);
                        remaining -= need;
                    }
                }
                if visited.len() >= MAX_SPREAD_CELLS {
                    continue;
                }
                let (ci, cj) = grid
                    .ij(c as u64)
                    .map_err(|e| ErosionError::Grid(e.to_string()))?;
                let mut next: Vec<usize> = grid
                    .neighbors(ci, cj)
                    .map_err(|e| ErosionError::Grid(e.to_string()))?
                    .into_iter()
                    .map(|(ti, tj, _)| ti * grid.nlon() + tj)
                    .filter(|&nb| routing.water_body_id[nb] == routing.water_body_id[start])
                    .collect();
                next.sort_unstable();
                for nb in next {
                    if visited.len() < MAX_SPREAD_CELLS && seen.insert(nb) {
                        visited.push(nb);
                    }
                }
            }
            if remaining <= 0.0 {
                break;
            }
            if visited.len() >= MAX_SPREAD_CELLS {
                rep.lake_spread_bound_hits += 1;
            }
            // A full lake (or a terminal floor): carry the rest downstream.
            let mut cur = start;
            let mut next_start: Option<usize> = None;
            loop {
                match routing.flow_to[cur] {
                    None => break,
                    Some(ElementId(j)) => {
                        let j = j as usize;
                        if ocean[j] || (is_lake(j) && routing.fill_level_m[j] > surf[j]) {
                            next_start = Some(j);
                            break;
                        }
                        cur = j;
                        hops += 1;
                        if hops > n {
                            return Err(ErosionError::Grid("cycle during lake overflow".into()));
                        }
                    }
                }
            }
            match next_start {
                Some(s) => {
                    if routing.water_body_id[s] != routing.water_body_id[start] {
                        *transfers.entry((cur, s, "lake_overflow")).or_default() += remaining;
                    }
                    start = s;
                }
                None => {
                    // No outlet: a terminal floor keeps the rest (counted as
                    // land floor storage; never silently dropped).
                    deposit(st, rep, cur, remaining / areas[cur], Sink::LandFloor);
                    remaining = 0.0;
                }
            }
        }
    }
    let mut out: Vec<(usize, f64, MarineClass)> = Vec::new();
    for t in 0..n {
        let v = marine_pending[t];
        if v <= 0.0 {
            continue;
        }
        rep.delivered_to_ocean_m3 += v;
        out.clear();
        let o = placer.place(adjacency, areas, ocean, surf, sea_level_m, shelf, t, v, &mut out);
        for &(c, h, class) in &out {
            deposit(st, rep, c, h, Sink::Marine(class));
        }
        rep.shelf_cells_filled += o.cells_filled;
        if o.bound_hit {
            rep.marine_search_bound_hits += 1;
        }
        *remainder_total += o.remainder_m3;
    }
    Ok(())
}

/// Legacy placement (the code of `erode_step_bounded`, unchanged in
/// effect) for `erosion-deposition-v1` without the shelf law: lakes and the
/// ocean filled to their water level with uncompensated thickness, a bounded
/// spread, and the legacy global layer. Used only to attribute effects
/// between the two laws; every raise is booked in the partition.
#[allow(clippy::too_many_arguments)]
fn place_legacy(
    grid: &GridSpec,
    routing: &Routing,
    elevation_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    ocean_cells: &[usize],
    ocean_area: f64,
    st: &mut Lowering<'_>,
    rep: &mut DepositionReport,
    transfers: &mut std::collections::BTreeMap<(usize, usize, &str), f64>,
    pending: &[f64],
    tally: &mut LayerTally,
) -> Result<(), ErosionError> {
    let n = pending.len();
    let areas = st.areas;
    let v1 = tally.law == OceanRemainderLayer::OncePerStepV1;
    let mut deferred = 0.0f64;
    let is_lake = |c: usize| !ocean[c] && routing.fill_level_m[c] > elevation_m[c] + 1.0e-9;
    let water_level = |c: usize| -> f64 {
        if ocean[c] {
            sea_level_m
        } else {
            routing.fill_level_m[c]
        }
    };
    for t0 in 0..n {
        let mut remaining = pending[t0];
        if remaining <= 0.0 {
            continue;
        }
        let mut start = t0;
        let mut hops = 0usize;
        while remaining > 0.0 {
            let is_ocean = ocean[start];
            if is_ocean {
                rep.delivered_to_ocean_m3 += remaining;
            }
            let mut visited = vec![start];
            let mut seen: std::collections::BTreeSet<usize> = [start].into_iter().collect();
            let mut head = 0usize;
            while head < visited.len() && remaining > 0.0 {
                let c = visited[head];
                head += 1;
                let capacity = (water_level(c) - st.z[c]).max(0.0) * areas[c];
                if capacity > 0.0 {
                    let take = remaining.min(capacity);
                    let sink = if is_ocean {
                        Sink::Marine(MarineClass::of_depth(sea_level_m - st.z[c]))
                    } else {
                        Sink::Lake
                    };
                    if v1 {
                        // Fill exactly to the water level (no ulp overshoot).
                        let level = water_level(c);
                        let h = if take == capacity { level - st.z[c] } else { take / areas[c] };
                        deposit(st, rep, c, h, sink);
                        st.clamp_to(c, level);
                    } else {
                        deposit(st, rep, c, take / areas[c], sink);
                    }
                    remaining -= take;
                }
                if visited.len() >= MAX_SPREAD_CELLS {
                    continue;
                }
                let (ci, cj) = grid
                    .ij(c as u64)
                    .map_err(|e| ErosionError::Grid(e.to_string()))?;
                let mut next: Vec<usize> = grid
                    .neighbors(ci, cj)
                    .map_err(|e| ErosionError::Grid(e.to_string()))?
                    .into_iter()
                    .map(|(ti, tj, _)| ti * grid.nlon() + tj)
                    .filter(|&nb| routing.water_body_id[nb] == routing.water_body_id[start])
                    .collect();
                next.sort_unstable();
                for nb in next {
                    if visited.len() < MAX_SPREAD_CELLS && seen.insert(nb) {
                        visited.push(nb);
                    }
                }
            }
            if remaining <= 0.0 {
                break;
            }
            if visited.len() >= MAX_SPREAD_CELLS {
                if is_ocean {
                    rep.marine_search_bound_hits += 1;
                } else {
                    rep.lake_spread_bound_hits += 1;
                }
            }
            if is_ocean {
                tally.mouths += 1;
                if v1 {
                    deferred += remaining; // laid once, after every mouth
                    break;
                }
                let layer = remaining / ocean_area;
                for &c in ocean_cells {
                    if st.z[c] <= sea_level_m && st.z[c] + layer > sea_level_m {
                        tally.above_sea += 1;
                    }
                    deposit(st, rep, c, layer, Sink::Remainder);
                }
                tally.raises += ocean_cells.len() as u64;
                break;
            }
            let mut cur = start;
            let mut next_start: Option<usize> = None;
            loop {
                match routing.flow_to[cur] {
                    None => break,
                    Some(ElementId(j)) => {
                        let j = j as usize;
                        if ocean[j] || (is_lake(j) && water_level(j) > st.z[j]) {
                            next_start = Some(j);
                            break;
                        }
                        cur = j;
                        hops += 1;
                        if hops > n {
                            return Err(ErosionError::Grid("cycle during lake overflow".into()));
                        }
                    }
                }
            }
            match next_start {
                Some(s) => {
                    if routing.water_body_id[s] != routing.water_body_id[start] {
                        *transfers.entry((cur, s, "lake_overflow")).or_default() += remaining;
                    }
                    start = s;
                }
                None => {
                    return Err(ErosionError::Grid(
                        "lake overflow has no downstream receiving body".into(),
                    ))
                }
            }
        }
    }
    if v1 && deferred > 0.0 {
        let (cells, _rounds, overfull) =
            capped_ocean_layer(&st.z, areas, ocean_cells, sea_level_m, deferred);
        tally.raises += cells.len() as u64;
        for (c, h, full) in cells {
            deposit(st, rep, c, h, Sink::Remainder);
            if full {
                st.clamp_to(c, sea_level_m);
            } else if overfull {
                tally.above_sea += 1;
            }
        }
    }
    Ok(())
}

/// Marine remainder layer selection and diagnostics for `place_legacy`.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct LayerTally {
    pub law: OceanRemainderLayer,
    pub mouths: u64,
    pub raises: u64,
    pub above_sea: u64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::routing::route;

    fn params(k: f64, d: f64) -> ErosionParams {
        ErosionParams {
            k_per_my: k,
            m: 0.5,
            n: 1.0,
            diffusivity_m2_per_my: d,
            stability_limit: 0.5,
        }
    }

    fn ramp_world() -> (GridSpec, Vec<f64>, Vec<bool>, Vec<f64>) {
        let grid = GridSpec::new(12, 24, 6.371e6).unwrap();
        let n = grid.cells();
        // Ocean in the southern rows; land rising northward by 200 m/row.
        let mut ocean = vec![false; n];
        let mut elev = vec![0.0; n];
        for i in 0..grid.nlat() {
            for j in 0..grid.nlon() {
                let c = i * grid.nlon() + j;
                if i < 3 {
                    ocean[c] = true;
                    elev[c] = -3000.0;
                } else {
                    elev[c] = 200.0 * (i as f64 - 2.0) + (j % 3) as f64 * 5.0;
                }
            }
        }
        let rain = vec![1.0e6; n];
        (grid, elev, ocean, rain)
    }

    #[test]
    fn params_validation_is_loud() {
        assert!(params(2.0, 1.0e5).validate().is_ok());
        let mut p = params(2.0, 1.0e5);
        p.n = 2.0;
        assert!(p.validate().is_err());
        let mut p = params(2.0, 1.0e5);
        p.stability_limit = 0.9;
        assert!(p.validate().is_err());
        let mut p = params(2.0, 1.0e5);
        p.k_per_my = -1.0;
        assert!(p.validate().is_err());
    }

    #[test]
    fn mass_closes_and_ocean_is_base_level() {
        let (grid, elev, ocean, rain) = ramp_world();
        let n = grid.cells();
        let r = route(&grid, &elev, &rain, &ocean, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let sed = vec![10.0; n];
        let step = erode_step(
            &grid,
            &stencil,
            &r,
            &elev,
            &sed,
            &ocean,
            0.0,
            2.0,
            &params(2.0, 1.0e6),
        )
        .unwrap();
        let residual = step.eroded_m3 - step.deposited_m3 - step.exported_m3;
        assert!(
            residual.abs() <= 1.0e-9 * step.eroded_m3.max(1.0),
            "eroded {} != deposited {} + exported {} (residual {residual})",
            step.eroded_m3,
            step.deposited_m3,
            step.exported_m3
        );
        assert!(step.eroded_m3 > 0.0, "a ramp must erode");
        assert!(step.exported_m3 > 0.0, "an ocean-connected ramp exports");
        for c in 0..n {
            if ocean[c] {
                assert!(
                    step.dz_m[c] >= 0.0,
                    "ocean cell {c} lowered by {}",
                    step.dz_m[c]
                );
                assert_eq!(step.bedrock_removed_m[c], 0.0);
            }
        }
        // Sediment removed before bedrock: a 10 m cover absorbs small lowering.
        let bad = (0..n)
            .filter(|&c| {
                !ocean[c]
                    && step.bedrock_removed_m[c] > 0.0
                    && step.sediment_removed_m[c] < 10.0 - 1.0e-9
            })
            .count();
        assert_eq!(bad, 0, "bedrock removed while sediment cover remained");
        // Volume identity per cell: dz = added - removed.
        for c in 0..n {
            let expect =
                step.sediment_added_m[c] - step.sediment_removed_m[c] - step.bedrock_removed_m[c];
            assert!(
                (step.dz_m[c] - expect).abs() < 1.0e-9,
                "cell {c} dz {} vs {expect}",
                step.dz_m[c]
            );
        }
    }

    #[test]
    fn implicit_stream_power_is_stable_for_huge_k_dt() {
        let (grid, elev, ocean, rain) = ramp_world();
        let r = route(&grid, &elev, &rain, &ocean, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let sed = vec![0.0; grid.cells()];
        // K dt so large an explicit scheme would explode by many orders.
        let step = erode_step(
            &grid,
            &stencil,
            &r,
            &elev,
            &sed,
            &ocean,
            0.0,
            50.0,
            &params(1.0e3, 0.0),
        )
        .unwrap();
        let z: Vec<f64> = elev.iter().zip(&step.dz_m).map(|(a, b)| a + b).collect();
        for i in 0..grid.cells() {
            if ocean[i] {
                continue;
            }
            assert!(z[i] <= elev[i] + step.sediment_added_m[i] + 1.0e-9);
            if let Some(ElementId(j)) = r.flow_to[i] {
                // The implicit update never cuts below the receiver's level
                // as it stood when the donor was processed: an ocean
                // receiver keeps its base level (deposits land on top of it
                // afterwards), a land receiver its own updated level.
                let j = j as usize;
                let floor = if ocean[j] {
                    0.0
                } else if r.fill_level_m[j] > elev[j] + 1e-9 {
                    r.fill_level_m[j]
                } else {
                    z[j]
                };
                assert!(
                    z[i] >= floor - 1.0e-6 || elev[i] <= elev[j],
                    "cell {i} cut below receiver {j}: {} < {floor}",
                    z[i]
                );
            }
        }
        assert!(step.dz_m.iter().all(|v| v.is_finite()));
        assert!(step.max_lowering_rate_m_per_my.is_finite());
    }

    #[test]
    fn diffusion_substeps_follow_the_declared_stability_condition() {
        let (grid, elev, ocean, rain) = ramp_world();
        let r = route(&grid, &elev, &rain, &ocean, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let sed = vec![0.0; grid.cells()];
        let p = params(0.0, 1.0e6);
        let nu = stencil.diffusion_number_max(p.diffusivity_m2_per_my, 2.0, &ocean);
        let step = erode_step(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p).unwrap();
        assert_eq!(step.diffusion_number_max, nu);
        let expect = ((nu / p.stability_limit).ceil() as u32).max(1);
        assert_eq!(step.diffusion_substeps, expect, "nu_max {nu}");
        // A diffusivity large enough that nu exceeds the limit on these
        // 15-degree cells (A ~ 2.8e12 m^2: nu ~ 2 dt D 4 / A needs D > ~2e11):
        // substeps must scale up, stay finite, stay closed. The value is a
        // numerics probe, not a physical diffusivity.
        let mut big = p;
        big.diffusivity_m2_per_my = 1.0e12;
        let step2 = erode_step(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &big).unwrap();
        assert!(
            step2.diffusion_substeps > step.diffusion_substeps,
            "{} vs {}",
            step2.diffusion_substeps,
            step.diffusion_substeps
        );
        let residual = step2.eroded_m3 - step2.deposited_m3 - step2.exported_m3;
        assert!(
            residual.abs() <= 1.0e-9 * step2.eroded_m3.max(1.0),
            "residual {residual}"
        );
        assert!(step2.exported_m3 > 0.0);
        // Sub-step diffusion number is under the declared limit.
        let nu_sub = step2.diffusion_number_max / step2.diffusion_substeps as f64;
        assert!(
            nu_sub <= big.stability_limit + 1.0e-12,
            "nu per substep {nu_sub}"
        );
    }

    #[test]
    fn deposition_spreads_across_the_water_body_up_to_its_level() {
        let (grid, elev, ocean, rain) = ramp_world();
        let n = grid.cells();
        let r = route(&grid, &elev, &rain, &ocean, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let sed = vec![0.0; n];
        // Huge K: the whole ramp is planed off in one step and exported.
        let step = erode_step(
            &grid,
            &stencil,
            &r,
            &elev,
            &sed,
            &ocean,
            0.0,
            50.0,
            &params(1.0e3, 0.0),
        )
        .unwrap();
        assert!(step.exported_m3 > 0.0);
        let receivers = (0..n)
            .filter(|&c| ocean[c] && step.sediment_added_m[c] > 0.0)
            .count();
        assert!(
            receivers > 8,
            "export landed on {receivers} ocean cells; a single-cell pile is the bug this guards"
        );
        for c in 0..n {
            if ocean[c] {
                let z = elev[c] + step.dz_m[c];
                assert!(
                    z <= 0.0 + 1.0e-6 || step.deposit_aggraded_m3 > 0.0,
                    "ocean cell {c} raised above sea level: {z}"
                );
            }
        }
        // The polar ocean rows of this fixture are tiny, so the local spread
        // is exhausted and the rest aggrades uniformly over the ocean body:
        // every ocean cell then carries the same extra layer on top.
        assert!(
            step.deposit_aggraded_m3 < step.exported_m3,
            "aggraded {} >= exported {}: nothing was spread locally",
            step.deposit_aggraded_m3,
            step.exported_m3
        );
        // No pile: every ocean cell receives at most its local capacity
        // (3000 m up to sea level) plus the uniform aggradation layer.
        let ocean_area: f64 = (0..n)
            .filter(|&c| ocean[c])
            .map(|c| stencil.areas_m2[c])
            .sum();
        let layer = step.deposit_aggraded_m3 / ocean_area;
        for c in 0..n {
            if ocean[c] {
                assert!(
                    step.sediment_added_m[c] <= 3000.0 + layer + 1.0e-6,
                    "ocean cell {c} received {} m, above capacity 3000 + layer {layer}",
                    step.sediment_added_m[c]
                );
            }
        }
        let residual = step.eroded_m3 - step.deposited_m3 - step.exported_m3;
        assert!(
            residual.abs() <= 1.0e-9 * step.eroded_m3,
            "residual {residual}"
        );
    }

    #[test]
    fn flat_land_drains_to_its_outlet_and_step_is_deterministic() {
        // A flat plateau with one ocean cell: with flat routing (review F2)
        // every cell drains to the outlet, so the implicit scheme lowers the
        // whole connected surface toward it in receiver-first order. What
        // must hold: determinism, no cell cut below its receiver, ocean
        // never lowered, mass closed.
        let grid = GridSpec::new(8, 16, 6.371e6).unwrap();
        let n = grid.cells();
        let mut ocean = vec![false; n];
        ocean[0] = true;
        let mut elev = vec![100.0; n];
        elev[0] = -100.0;
        let rain = vec![1.0e6; n];
        let r = route(&grid, &elev, &rain, &ocean, 0.0).unwrap();
        assert!(
            (1..n).all(|c| r.flow_to[c].is_some()),
            "flat routing must give every land cell a receiver"
        );
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let sed = vec![0.0; n];
        let a = erode_step(
            &grid,
            &stencil,
            &r,
            &elev,
            &sed,
            &ocean,
            0.0,
            2.0,
            &params(2.0, 1.0e5),
        )
        .unwrap();
        let b = erode_step(
            &grid,
            &stencil,
            &r,
            &elev,
            &sed,
            &ocean,
            0.0,
            2.0,
            &params(2.0, 1.0e5),
        )
        .unwrap();
        assert_eq!(a, b);
        let fluvial_only = erode_step(
            &grid,
            &stencil,
            &r,
            &elev,
            &sed,
            &ocean,
            0.0,
            2.0,
            &params(2.0, 0.0),
        )
        .unwrap();
        let z: Vec<f64> = elev
            .iter()
            .zip(&fluvial_only.dz_m)
            .map(|(e, d)| e + d)
            .collect();
        for c in 1..n {
            let ElementId(j) = r.flow_to[c].unwrap();
            let j = j as usize;
            let floor = if ocean[j] {
                0.0
            } else if r.fill_level_m[j] > elev[j] + 1e-9 {
                r.fill_level_m[j]
            } else {
                z[j]
            };
            assert!(
                z[c] >= floor - 1.0e-6,
                "cell {c} cut below its receiver {j}"
            );
            assert!(
                z[c] <= elev[c] + 1.0e-9,
                "land cell {c} raised by fluvial action alone"
            );
        }
        assert!(fluvial_only.dz_m[0] >= 0.0, "ocean lowered");
        let residual =
            fluvial_only.eroded_m3 - fluvial_only.deposited_m3 - fluvial_only.exported_m3;
        assert!(
            residual.abs() <= 1.0e-9 * fluvial_only.eroded_m3.max(1.0),
            "residual {residual}"
        );
        let movers = (0..n).filter(|&c| fluvial_only.dz_m[c] < 0.0).count();
        let movers_d = (0..n).filter(|&c| a.dz_m[c] < 0.0).count();
        assert!(
            movers > 0 && movers_d >= movers,
            "fluvial {movers}, with diffusion {movers_d}"
        );
    }

    /// Lane INCISION-THRESHOLD part 1: a shallow sea that the bounded local
    /// fill saturates. V0 lifts filled cells above sea level and costs
    /// mouths x ocean cells; V1 lays the same volume once, capped at sea
    /// level; both close mass; the V0 entry point is unchanged.
    #[test]
    fn once_per_step_layer_is_capped_conservative_and_cheaper() {
        let grid = GridSpec::new(40, 80, 6.371e6).unwrap();
        let n = grid.cells();
        let mut ocean = vec![false; n];
        let mut elev = vec![0.0; n];
        for i in 0..grid.nlat() {
            for j in 0..grid.nlon() {
                let c = i * grid.nlon() + j;
                if i < 20 {
                    ocean[c] = true;
                    // Shallow shelf (0.5 m) with a deep trough in the far south.
                    elev[c] = if i < 3 { -4000.0 } else { -0.5 };
                } else {
                    elev[c] = 400.0 * (i as f64 - 19.0) + (j % 3) as f64 * 5.0;
                }
            }
        }
        let rain = vec![1.0e6; n];
        let r = route(&grid, &elev, &rain, &ocean, 0.0).unwrap();
        let stencil = DiffusionStencil::build(&grid).unwrap();
        let sed = vec![0.0; n];
        let p = params(0.05, 0.0);
        let legacy = erode_step(&grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p).unwrap();
        let v0 = erode_step_bounded_layer(
            &grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, None,
            OceanRemainderLayer::PerMouthV0,
        )
        .unwrap();
        assert_eq!(legacy, v0, "V0 entry point must be the legacy step");
        let v1 = erode_step_bounded_layer(
            &grid, &stencil, &r, &elev, &sed, &ocean, 0.0, 2.0, &p, None,
            OceanRemainderLayer::OncePerStepV1,
        )
        .unwrap();
        assert!(v0.ocean_layer_mouths > 1, "fixture must saturate several mouths");
        assert_eq!(v0.ocean_layer_mouths, v1.ocean_layer_mouths);
        assert!(v0.ocean_layer_cells_above_sea > 0, "V0 lifts filled cells above sea");
        assert_eq!(v1.ocean_layer_cells_above_sea, 0);
        assert!(v1.ocean_layer_cell_raises < v0.ocean_layer_cell_raises);
        for s in [&v0, &v1] {
            let residual = s.eroded_m3 - s.deposited_m3 - s.exported_m3;
            assert!(residual.abs() <= 1.0e-9 * s.eroded_m3, "residual {residual}");
        }
        assert_eq!(v0.eroded_m3, v1.eroded_m3);
        assert!((v0.exported_m3 - v1.exported_m3).abs() <= 1.0e-9 * v0.exported_m3);
        for c in 0..n {
            if ocean[c] {
                assert!(elev[c] + v1.dz_m[c] <= 0.0, "V1 ocean cell {c} above sea: {}", elev[c] + v1.dz_m[c]);
            }
        }
        assert!((0..n).any(|c| ocean[c] && elev[c] + v0.dz_m[c] > 0.0));
    }

    #[test]
    fn capped_layer_fills_shallow_cells_and_spreads_the_rest() {
        let z = [-1.0, -10.0, -100.0, 0.0, 5.0];
        let areas = [1.0, 1.0, 2.0, 1.0, 1.0];
        let ocean_cells = [0usize, 1, 2, 3];
        // 30 m3: cell 0 fills (1), the rest spreads over cells 1, 2 (area 3).
        let (cells, rounds, overfull) = capped_ocean_layer(&z, &areas, &ocean_cells, 0.0, 30.0);
        assert!(!overfull);
        assert!(rounds >= 2);
        let vol: f64 = cells.iter().map(|&(c, h, _)| h * areas[c]).sum();
        assert!((vol - 30.0).abs() < 1e-12);
        let h = |k: usize| cells.iter().filter(|e| e.0 == k).map(|e| e.1).sum::<f64>();
        assert_eq!(h(0), 1.0);
        assert!((h(1) - 29.0 / 3.0).abs() < 1e-12 && (h(2) - 29.0 / 3.0).abs() < 1e-12);
        assert_eq!(h(3), 0.0, "a cell at sea level gets nothing");
        // A full ocean falls back to the legacy uniform layer, flagged.
        let (_, _, overfull) = capped_ocean_layer(&[0.0, 0.0], &[1.0, 1.0], &[0, 1], 0.0, 1.0);
        assert!(overfull);
    }
}
