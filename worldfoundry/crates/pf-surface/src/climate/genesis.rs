//! CLIMATE-GENESIS (lane CG, 2026-10-20): climate forcing of T0 genesis
//! erosion, law `epoch-climate-runoff-v1`
//! (design: `unpublished-design-record`).
//!
//! PHYSICS. Stream power on discharge, E = K_Q Q^m S^n, written in the
//! area form the kernels already use: E = K A_eff^m S^n with
//! A_eff,i = sum over the catchment of i of a_j R_j / r_ref
//! (`crate::routing::runoff_weighted_area_m2`). K is the erodibility at the
//! reference runoff r_ref (`surface.uniform_runoff_m_per_yr`), so a dry cell
//! fed by a wet mountain catchment still carries the wet discharge.
//!
//! RUNOFF per cell (m/yr): land max(P - E, f P) (`land_p_minus_e_floor_v1`,
//! the region chain's rule); lake (fill above the bed, not ocean)
//! max(P - PET, 0); ocean r_ref (no fluvial erosion there).
//!
//! TIME. A full `seed-climate-v1`/`-v2` solve runs at steps
//! k = 0, e, 2e, ... (e = epoch / dt; prefix-stable), on the genesis grid
//! coarsened by the declared divisor (area mean of the water/land surface,
//! ocean by majority with ties to ocean), warm-started from the previous
//! epoch's converged column water and tracer, to a declared tolerance; the
//! solver's water closure (1e-6) is checked on every solve (typed failure).
//! Between epochs (`frozen-column-local-rainout-v1`) the column water,
//! tracer and zonal rows are held and P, E and R are recomputed pointwise
//! on the genesis grid from the current surface: the column is
//! bilinearly interpolated from the solve grid's cell centres (periodic in
//! longitude, clamped at the poles) and the local rainout uses the
//! genesis-grid surface, so orographic rain follows ranges as plates carry
//! them while the upwind depletion (the shadow) lags until the next solve.
//! The same local rainout gives the runoff at the epoch step itself (with
//! divisor 1 it reproduces the solver's precipitation bit for bit).
//!
//! DETERMINISM: pure functions of the inputs and the carried state; fixed
//! loops, no threads, no randomness, no clock in any carried value.

use crate::climate::rotation::{derive_seed_climate_v2_with, RotationClimateConstantsV2};
use crate::climate::seed::{
    budyko_ratio, derive_seed_climate_v1_with, seed_climate_v1_earth_recipe, SeedClimateConstantsV1, SeedColumnStateV1,
    SeedSolveOptions, SEED_SWEEP_REL_TOL,
};
use crate::climate::{saturation_column_kg_m2, thornthwaite_pet_mm_per_yr, ClimateError, ClimateInputs, ClimateOutputV0, SECONDS_PER_YEAR};
use pf_contracts::laws::{ClimateBetweenEpochsV1, ClimateErosivityModeV1, ClimateForcingLawV1, ClimateLawV1, ClimateRecipeV1, ClimateSubgridCellMeanV1, ErosiveRunoffLawV1};
use pf_contracts::params::ClimateForcingBlock;
use pf_grid::GridSpec;
use serde::{Deserialize, Serialize};

pub const GENESIS_CLIMATE_LAW_ID: &str = "epoch-climate-runoff-v1";

/// v1 defaults. Project assumptions awaiting review; the numerics (epoch,
/// divisor, tolerance) were chosen from the S0 timing probe on C6
/// (`unpublished-design-record`: L1 cold 1e-11
/// 19.3 s; divisor 2 cold 1e-7 0.77 s with water residual 3.2e-7; divisor 1
/// at 1e-7 fails the 1e-6 water closure; divisor 1 vs 2 land R^0.5
/// correlation 0.997). Measured in genesis (seed 42, 400 My): 100 solves,
/// 61 s of solves, +76 s wall against the uniform arm.
pub const GENESIS_V1_EPOCH_MY: f64 = 4.0;
pub const GENESIS_V1_SOLVE_GRID_DIVISOR: f64 = 2.0;
pub const GENESIS_V1_SWEEP_REL_TOL: f64 = 1.0e-7;
/// Lane CLIMATE-RUNOFF's `land_p_minus_e_floor_v1` f (same constant).
pub const GENESIS_V1_RUNOFF_FLOOR_FRACTION: f64 = 0.1;

/// I12 (lane CLIMATE-GENESIS-2, 2026-10-20): water-closure retry ladder.
/// A solve whose atmospheric-water closure misses the declared 1e-6 is
/// re-run (same inputs, same warm start) at the tolerance divided by this
/// factor, down to the post-hoc tier's `SEED_SWEEP_REL_TOL` (1e-11), and
/// only then fails loud. Numerics, not physics: a world whose solves all
/// close at the declared tolerance is byte-identical with or without the
/// ladder (the retry runs only where the old code stopped with an error).
/// Measured margin at 1e-7 before the ladder: worst residual 6.4-8.1e-7
/// over 3 x 100 solves (lane CLIMATE-GENESIS REPORT). Project assumption
/// awaiting review: factor 100 (two decades of sweeps, about 2 x 100
/// sweeps per e-fold at the measured line Gauss-Seidel rate).
pub const GENESIS_CLOSURE_RETRY_TOL_FACTOR: f64 = 100.0;

/// Approximations this law adds (recorded by the generator).
pub const GENESIS_CLIMATE_APPROXIMATIONS_V1: &[(&str, &str)] = &[
    ("runoff-weighted-area-v1", "every fluvial kernel reads A_eff = sum over the catchment of a R / r_ref in place of the geometric area (stream power on discharge, K quoted at r_ref = surface.uniform_runoff_m_per_yr); flow directions, fill levels, basins, published hydrology and the water-exchange discharge stay geometric / uniform-runoff in v1"),
    ("land-p-minus-e-floor-v1", "land runoff R = max(P - E, f P); lake runoff max(P - PET, 0) (no lake evaporation of inflow, so basins never close in v1); ocean R = r_ref"),
    ("epoch-coarse-warm-climate-v1", "full climate solve every epoch on the genesis grid coarsened by the declared divisor (area-mean surface, majority ocean, ties to ocean), warm-started, to the declared tolerance; atmospheric water closes at each solve (1e-6), not between solves"),
    ("frozen-column-local-rainout-v1", "between solves the column water, tracer and zonal rows are held and P, E and R recomputed pointwise on the genesis grid from the current surface (bilinear column from solve-grid centres); the rain shadow lags by up to one epoch of drift"),
    ("hillslope-threshold-climate-independent-v1", "the critical-slope envelope, linear hillslope diffusion, the subgrid channel-head length and (unless cell-runoff-subgrid-erosivity-v1 is selected) the subgrid cell-mean relaxation rate stay at the reference runoff (climate reaches the subgrid law through the trunk incision only)"),
];

/// I11 (lane CLIMATE-GENESIS-2): the approximation recorded when
/// `subgrid_cell_mean = cell-runoff-subgrid-erosivity-v1`.
pub const GENESIS_CLIMATE_SUBGRID_APPROXIMATION_V1: (&str, &str) = (
    "cell-runoff-subgrid-erosivity-v1",
    "the subgrid cell-mean relaxation of the subgrid erosion law runs at K_n (R_i/r_ref)^m (m = n/2, R_i the cell's own runoff fed to erosion): the unresolved network drains the cell's own area; the critical-slope floor, diffusion and channel-head length stay climate-independent",
);

/// The v1 block: `seed-climate-v1` Earth recipe and constants, K quoted
/// at r_ref (keep K 2e-6; never stack with the uniform-runoff K 1.0e-6).
pub fn climate_forcing_block_v1(erosivity: ClimateErosivityModeV1) -> Result<ClimateForcingBlock, ClimateError> {
    Ok(ClimateForcingBlock {
        climate: seed_climate_v1_earth_recipe(),
        seed_constants_digest: SeedClimateConstantsV1::earth_v1()
            .digest()
            .map_err(|e| ClimateError::Inputs(format!("seed constants digest: {e}")))?,
        rotation_constants_digest: None,
        epoch_my: GENESIS_V1_EPOCH_MY,
        solve_grid_divisor: GENESIS_V1_SOLVE_GRID_DIVISOR,
        sweep_rel_tol: GENESIS_V1_SWEEP_REL_TOL,
        runoff_floor_fraction: GENESIS_V1_RUNOFF_FLOOR_FRACTION,
        between_epochs: ClimateBetweenEpochsV1::FrozenColumnLocalRainoutV1,
        erosivity,
        subgrid_cell_mean: ClimateSubgridCellMeanV1::ReferenceRunoffV1,
        erosive_runoff_exponent: None,
        erosive_runoff_cap_m_per_yr: None,
    })
}

/// CLIMATE-ARIDITY (2026-10-21): the v1 block with the opt-in
/// `crh-rainout-v1` seed-constants set (`SeedClimateConstantsV1::earth_crh_v1`).
pub fn climate_forcing_block_crh_v1(erosivity: ClimateErosivityModeV1) -> Result<ClimateForcingBlock, ClimateError> {
    let mut b = climate_forcing_block_v1(erosivity)?;
    b.seed_constants_digest = SeedClimateConstantsV1::earth_crh_v1()
        .digest()
        .map_err(|e| ClimateError::Inputs(format!("seed constants digest: {e}")))?;
    Ok(b)
}

/// The v1 block on `seed-climate-v2` (rotation-scaled circulation bands).
pub fn climate_forcing_block_v2(erosivity: ClimateErosivityModeV1) -> Result<ClimateForcingBlock, ClimateError> {
    let mut b = climate_forcing_block_v1(erosivity)?;
    b.climate = crate::climate::rotation::seed_climate_v2_earth_recipe();
    b.rotation_constants_digest = Some(
        RotationClimateConstantsV2::earth_v2()
            .digest()
            .map_err(|e| ClimateError::Inputs(format!("rotation constants digest: {e}")))?,
    );
    Ok(b)
}

/// Pointwise constants of the local rainout (from the recipe).
#[derive(Debug, Clone, Copy, PartialEq)]
struct LocalConstants {
    lapse: f64,
    w_ref: f64,
    t_ref: f64,
    cc: f64,
    tau_r: f64,
    tau_f: f64,
    et_cap: f64,
    damp: f64,
    /// CLIMATE-ARIDITY (2026-10-21): the opt-in crh rainout of the seed
    /// constants set (None for `earth_v1()`: v1 arithmetic).
    crh: Option<crate::climate::aridity::CrhParams>,
}

impl LocalConstants {
    fn of(recipe: &ClimateRecipeV1, extra: &SeedClimateConstantsV1) -> Self {
        let k = &recipe.constants;
        Self {
            lapse: k.lapse_rate_k_per_m.value,
            w_ref: k.column_water_saturation_ref_kg_m2.value,
            t_ref: k.column_water_reference_temperature_c.value,
            cc: k.clausius_clapeyron_rate_per_k.value,
            tau_r: k.rainout_residence_s.value,
            tau_f: extra.orographic_fallout_length_m.value / k.moisture_transport_speed_m_s.value,
            et_cap: extra.land_et_ratio_cap.value,
            damp: k.ocean_seasonal_damping.value,
            crh: extra.crh_rainout.as_ref().map(|c| c.params()),
        }
    }
}

/// The resolved law (pure function of the law, dt, grid and r_ref).
#[derive(Debug, Clone, PartialEq)]
pub struct GenesisClimateParams {
    pub recipe: ClimateRecipeV1,
    pub seed: SeedClimateConstantsV1,
    pub rotation: Option<RotationClimateConstantsV2>,
    pub epoch_steps: u32,
    pub divisor: usize,
    pub sweep_rel_tol: f64,
    pub floor_fraction: f64,
    pub between: ClimateBetweenEpochsV1,
    pub erosivity: ClimateErosivityModeV1,
    /// I11: climate runoff in the subgrid cell-mean relaxation.
    pub subgrid_cell_mean: ClimateSubgridCellMeanV1,
    /// HYPSO-LEVERS: the erosive-runoff response (None = linear v1).
    pub erosive_runoff: Option<ErosiveRunoffLawV1>,
    /// Reference runoff (m/yr) at which K is quoted.
    pub r_ref_m_per_yr: f64,
    /// Stream-power area exponent m (land-mean diagnostic only).
    pub m: f64,
    local: LocalConstants,
}

impl GenesisClimateParams {
    /// Resolve the law; `None` for the legacy law. Refuses (fail loud) any
    /// combination this build does not implement: tier constants other
    /// than the compiled Earth sets (digest mismatch), a v2 gravity unequal
    /// to the universe's, an epoch that is not a whole number of steps, a
    /// grid the divisor does not tile.
    pub fn from_law(
        law: &ClimateForcingLawV1,
        dt_my: f64,
        grid: &GridSpec,
        r_ref_m_per_yr: f64,
        m: f64,
        universe_gravity_m_s2: f64,
    ) -> Result<Option<Self>, String> {
        let ClimateForcingLawV1::EpochClimateRunoffV1 {
            climate,
            seed_constants_digest,
            rotation_constants_digest,
            climate_epoch_my,
            climate_solve_grid_divisor,
            climate_sweep_rel_tol,
            climate_runoff_floor_fraction,
            between_epochs,
            erosivity,
            subgrid_cell_mean,
            erosive_runoff,
        } = law
        else {
            return Ok(None);
        };
        // CLIMATE-ARIDITY (2026-10-21): the digest selects one of the two
        // compiled sets: earth_v1 (v1 rainout) or earth_crh_v1 (opt-in
        // crh-rainout-v1). Any other digest is refused.
        let seed_v1 = SeedClimateConstantsV1::earth_v1();
        let seed_digest = seed_v1.digest().map_err(|e| e.to_string())?;
        let seed_crh = SeedClimateConstantsV1::earth_crh_v1();
        let crh_digest = seed_crh.digest().map_err(|e| e.to_string())?;
        let seed = if &seed_digest == seed_constants_digest {
            seed_v1
        } else if &crh_digest == seed_constants_digest {
            seed_crh
        } else {
            return Err(format!(
                "{GENESIS_CLIMATE_LAW_ID}: seed climate constants digest {seed_constants_digest} is neither the compiled earth_v1 set {seed_digest} nor the earth_crh_v1 set {crh_digest} (no other set is implemented)"
            ));
        };
        let rotation = match (climate.law, rotation_constants_digest) {
            (ClimateLawV1::SeedClimateV2, Some(d)) => {
                let rot = RotationClimateConstantsV2::earth_v2();
                let want = rot.digest().map_err(|e| e.to_string())?;
                if &want != d {
                    return Err(format!(
                        "{GENESIS_CLIMATE_LAW_ID}: rotation climate constants digest {d} is not the compiled earth_v2 set {want}"
                    ));
                }
                if rot.surface_gravity_m_s2.value.to_bits() != universe_gravity_m_s2.to_bits() {
                    return Err(format!(
                        "{GENESIS_CLIMATE_LAW_ID}: seed-climate-v2 surface gravity {} differs from the universe's {universe_gravity_m_s2} (one physical constant, one value)",
                        rot.surface_gravity_m_s2.value
                    ));
                }
                Some(rot)
            }
            (ClimateLawV1::SeedClimateV1, None) => None,
            (l, d) => return Err(format!("{GENESIS_CLIMATE_LAW_ID}: climate law {l:?} with rotation digest {d:?} is not supported")),
        };
        let epoch = climate_epoch_my.value;
        let k = (epoch / dt_my).round();
        if !(k >= 1.0) || (k * dt_my - epoch).abs() > 1e-9 * epoch {
            return Err(format!("{GENESIS_CLIMATE_LAW_ID}: epoch {epoch} My is not a whole number of {dt_my} My steps"));
        }
        let divisor = climate_solve_grid_divisor.value as usize;
        if divisor as f64 != climate_solve_grid_divisor.value || divisor == 0 || grid.nlat() % divisor != 0 || grid.nlon() % divisor != 0 {
            return Err(format!(
                "{GENESIS_CLIMATE_LAW_ID}: solve grid divisor {} does not tile the {} x {} grid",
                climate_solve_grid_divisor.value,
                grid.nlat(),
                grid.nlon()
            ));
        }
        if !(r_ref_m_per_yr.is_finite() && r_ref_m_per_yr > 0.0) {
            return Err(format!("{GENESIS_CLIMATE_LAW_ID}: reference runoff {r_ref_m_per_yr} m/yr must be positive (K is quoted at it)"));
        }
        if !(m.is_finite() && m > 0.0) {
            return Err(format!("{GENESIS_CLIMATE_LAW_ID}: area exponent m {m} must be positive"));
        }
        Ok(Some(Self {
            local: LocalConstants::of(climate, &seed),
            recipe: climate.clone(),
            seed,
            rotation,
            epoch_steps: k as u32,
            divisor,
            sweep_rel_tol: climate_sweep_rel_tol.value,
            floor_fraction: climate_runoff_floor_fraction.value,
            between: *between_epochs,
            erosivity: *erosivity,
            subgrid_cell_mean: *subgrid_cell_mean,
            erosive_runoff: erosive_runoff.clone(),
            r_ref_m_per_yr,
            m,
        }))
    }

    /// HYPSO-LEVERS (2026-10-21): the ERODING runoff R_e the fluvial
    /// kernels read (`ErosiveRunoffLawV1`); `None` under the v1 law (R_e =
    /// R, the caller keeps using `runoff`). `runoff` is the modelled field
    /// (m/yr) after the erosivity mode; the lake balance keeps it.
    pub fn erosive_runoff_m_per_yr(&self, runoff: &[f64]) -> Option<Vec<f64>> {
        let r0 = self.r_ref_m_per_yr;
        match self.erosive_runoff.as_ref()? {
            ErosiveRunoffLawV1::SublinearAboveReferenceV1 { erosive_runoff_exponent } => {
                let p = erosive_runoff_exponent.value;
                Some(runoff.iter().map(|&r| if r > r0 { r0 * (r / r0).powf(p) } else { r }).collect())
            }
            ErosiveRunoffLawV1::CapV1 { erosive_runoff_cap_m_per_yr } => {
                let cap = erosive_runoff_cap_m_per_yr.value;
                Some(runoff.iter().map(|&r| r.min(cap)).collect())
            }
        }
    }

    /// True when step index `k` (0-based, steps already taken) solves.
    pub fn solves_at(&self, k: u32) -> bool {
        k % self.epoch_steps == 0
    }
}

/// One solve epoch (checkpointed; deterministic).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenesisClimateEpochLogV1 {
    pub step: u32,
    pub t_my: f64,
    pub tracer_sweeps: u32,
    pub moisture_sweeps: u32,
    /// Atmospheric-water closure residual / gross precipitation.
    pub water_residual_rel: f64,
    pub land_mean_precipitation_mm_per_yr: Option<f64>,
    pub max_precipitation_mm_per_yr: f64,
    pub land_orographic_share: Option<f64>,
    /// I12: retries of the closure ladder this solve needed (0 = closed at
    /// the declared tolerance; skipped on write when 0, so pre-ladder
    /// checkpoints keep their bytes).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub closure_retries: u32,
}

fn is_zero_u32(v: &u32) -> bool {
    *v == 0
}

/// One genesis step's runoff diagnostics (checkpointed; deterministic).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenesisClimateStepLogV1 {
    pub step: u32,
    pub t_my: f64,
    pub solved: bool,
    /// Area-weighted means over non-ocean cells of the runoff fed to erosion.
    pub land_mean_runoff_m_per_yr: f64,
    /// ... and of (R / r_ref)^m (the mean erosivity factor).
    pub land_mean_erosivity_factor: f64,
    pub max_land_runoff_m_per_yr: f64,
    /// Orographic share of land precipitation this step (local rainout).
    pub land_orographic_share: f64,
}

/// The run state of the law (checkpointed except the `last_*` fields,
/// which are observation-only copies of the latest step's fields).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GenesisClimateState {
    /// Converged column of the latest solve epoch (solve grid).
    pub column: Option<SeedColumnStateV1>,
    /// `hold-runoff-v1` only: the latest epoch's runoff (genesis grid, m/yr).
    pub held_runoff_m_per_yr: Option<Vec<f64>>,
    pub epochs: Vec<GenesisClimateEpochLogV1>,
    pub steps: Vec<GenesisClimateStepLogV1>,
    /// Observation only (not checkpointed): the latest step's runoff fed to
    /// erosion (m/yr), precipitation and its orographic part (mm/yr).
    pub last_runoff_m_per_yr: Vec<f64>,
    pub last_precipitation_mm_per_yr: Vec<f64>,
    pub last_orographic_mm_per_yr: Vec<f64>,
    /// Observation only (lane INLAND-BASINS, 2026-10-21): the latest step's
    /// potential evapotranspiration (mm/yr), read by the opt-in
    /// closed-basin water balance (never fed back into the climate).
    pub last_pet_mm_per_yr: Vec<f64>,
    /// Observation only: wall seconds of the latest solve (never in state).
    pub last_solve_wall_s: f64,
}

/// Pointwise fields of a local rainout on the genesis grid.
#[derive(Debug, Clone, PartialEq)]
pub struct LocalRainout {
    pub precipitation_mm_per_yr: Vec<f64>,
    pub orographic_mm_per_yr: Vec<f64>,
    pub evaporation_mm_per_yr: Vec<f64>,
    pub pet_mm_per_yr: Vec<f64>,
}

/// Area mean of the water/land surface over divisor x divisor blocks,
/// ocean by majority (ties to ocean). Returns (grid, surface, ocean).
pub fn coarsen_surface(
    grid: &GridSpec,
    divisor: usize,
    elevation_m: &[f64],
    fill_level_m: &[f64],
    ocean: &[bool],
) -> Result<(GridSpec, Vec<f64>, Vec<bool>), ClimateError> {
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    let n = grid.cells();
    if elevation_m.len() != n || fill_level_m.len() != n || ocean.len() != n {
        return Err(ClimateError::Inputs("coarsen: field sizes differ from the grid".into()));
    }
    if divisor == 0 || nlat % divisor != 0 || nlon % divisor != 0 {
        return Err(ClimateError::Inputs(format!("coarsen: divisor {divisor} does not tile {nlat} x {nlon}")));
    }
    if divisor == 1 {
        // Exact (no area-weighted division): the solve sees the surface itself.
        let surf = (0..n).map(|c| elevation_m[c].max(fill_level_m[c])).collect();
        return Ok((grid.clone(), surf, ocean.to_vec()));
    }
    let (cl, cn) = (nlat / divisor, nlon / divisor);
    let coarse = GridSpec::new(cl, cn, grid.radius_m()).map_err(|e| ClimateError::Inputs(e.to_string()))?;
    let row_area: Vec<f64> = (0..nlat)
        .map(|i| grid.cell_area_m2(i).map_err(|e| ClimateError::Inputs(e.to_string())))
        .collect::<Result<_, _>>()?;
    let mut surf = vec![0.0; cl * cn];
    let mut oc = vec![false; cl * cn];
    for ci in 0..cl {
        for cj in 0..cn {
            let (mut s, mut a, mut wet) = (0.0, 0.0, 0usize);
            for di in 0..divisor {
                let i = ci * divisor + di;
                for dj in 0..divisor {
                    let c = i * nlon + cj * divisor + dj;
                    s += row_area[i] * elevation_m[c].max(fill_level_m[c]);
                    a += row_area[i];
                    if ocean[c] {
                        wet += 1;
                    }
                }
            }
            surf[ci * cn + cj] = s / a;
            oc[ci * cn + cj] = 2 * wet >= divisor * divisor;
        }
    }
    Ok((coarse, surf, oc))
}

/// Bilinear interpolation of a solve-grid cell field to the genesis grid
/// (cell centres; periodic in longitude, clamped at the poles). Divisor 1
/// returns the field unchanged.
fn upsample_cells(v: &[f64], cl: usize, cn: usize, d: usize) -> Vec<f64> {
    if d == 1 {
        return v.to_vec();
    }
    let (nlat, nlon) = (cl * d, cn * d);
    let mut out = vec![0.0; nlat * nlon];
    for i in 0..nlat {
        let (i0, i1, fy) = row_weights(i, cl, d);
        for j in 0..nlon {
            let x = (j as f64 + 0.5) / d as f64 - 0.5;
            let x0 = x.floor();
            let fx = x - x0;
            let j0 = (x0 as i64).rem_euclid(cn as i64) as usize;
            let j1 = (j0 + 1) % cn;
            let a = (1.0 - fx) * v[i0 * cn + j0] + fx * v[i0 * cn + j1];
            let b = (1.0 - fx) * v[i1 * cn + j0] + fx * v[i1 * cn + j1];
            out[i * nlon + j] = (1.0 - fy) * a + fy * b;
        }
    }
    out
}

fn row_weights(i: usize, cl: usize, d: usize) -> (usize, usize, f64) {
    let y = (i as f64 + 0.5) / d as f64 - 0.5;
    if y <= 0.0 {
        (0, 0, 0.0)
    } else if y >= (cl - 1) as f64 {
        (cl - 1, cl - 1, 0.0)
    } else {
        let y0 = y.floor();
        (y0 as usize, y0 as usize + 1, y - y0)
    }
}

fn upsample_rows(v: &[f64], cl: usize, d: usize) -> Vec<f64> {
    if d == 1 {
        return v.to_vec();
    }
    (0..cl * d)
        .map(|i| {
            let (i0, i1, fy) = row_weights(i, cl, d);
            (1.0 - fy) * v[i0] + fy * v[i1]
        })
        .collect()
}

/// GLACIAL-T0 (lane GLACIAL-T0, 2026-10-21): the temperature fields the
/// opt-in `glacial-stream-power-blend-v1` reads, per genesis cell, from a
/// (held) solve column with exactly `local_rainout`'s formulas: the zonal
/// sea-level temperature (deg C) and the seasonal half range (K; ocean
/// damping, land continentality from the tracer). Also returns the recipe's
/// lapse rate (K/m), the one value the climate uses. Observation of the
/// climate only; nothing here feeds back into it.
pub fn glacial_temperature_inputs(
    p: &GenesisClimateParams,
    column: &SeedColumnStateV1,
    grid: &GridSpec,
    ocean: &[bool],
) -> Result<(Vec<f64>, Vec<f64>, f64), ClimateError> {
    let (nlat, nlon, n) = (grid.nlat(), grid.nlon(), grid.cells());
    let d = p.divisor;
    if column.nlat * d != nlat || column.nlon * d != nlon || ocean.len() != n {
        return Err(ClimateError::Inputs(format!(
            "glacial inputs: column {} x {} times divisor {d} is not the grid {nlat} x {nlon} (or ocean size {})",
            column.nlat,
            column.nlon,
            ocean.len()
        )));
    }
    let k = p.local;
    let tracer = upsample_cells(&column.tracer, column.nlat, column.nlon, d);
    let t_sl = upsample_rows(&column.sea_level_temperature_c, column.nlat, d);
    let resp = upsample_rows(&column.half_range_response_k, column.nlat, d);
    let mut ts = vec![0.0; n];
    let mut hr = vec![0.0; n];
    for c in 0..n {
        let i = c / nlon;
        let m = if ocean[c] {
            k.damp
        } else {
            let t = tracer[c].clamp(0.0, 1.0);
            k.damp + (1.0 - k.damp) * (1.0 - t)
        };
        ts[c] = t_sl[i];
        hr[c] = m * resp[i];
    }
    Ok((ts, hr, k.lapse))
}

/// Pointwise rainout on the genesis grid from a (held) solve column: the
/// `seed-climate-v1` formulas of precipitation, its orographic part, land
/// evapotranspiration (Budyko, capped) and Thornthwaite PET, evaluated at
/// the current surface. With divisor 1 and the solve's own inputs it
/// reproduces the solver's fields bit for bit.
pub fn local_rainout(
    p: &GenesisClimateParams,
    column: &SeedColumnStateV1,
    grid: &GridSpec,
    elevation_m: &[f64],
    fill_level_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
) -> Result<LocalRainout, ClimateError> {
    let (nlat, nlon, n) = (grid.nlat(), grid.nlon(), grid.cells());
    let d = p.divisor;
    if column.nlat * d != nlat || column.nlon * d != nlon {
        return Err(ClimateError::Inputs(format!(
            "local rainout: column {} x {} times divisor {d} is not the grid {nlat} x {nlon}",
            column.nlat, column.nlon
        )));
    }
    if elevation_m.len() != n || fill_level_m.len() != n || ocean.len() != n {
        return Err(ClimateError::Inputs("local rainout: field sizes differ from the grid".into()));
    }
    let k = p.local;
    let w = upsample_cells(&column.column_water_kg_m2, column.nlat, column.nlon, d);
    let tracer = upsample_cells(&column.tracer, column.nlat, column.nlon, d);
    let t_sl = upsample_rows(&column.sea_level_temperature_c, column.nlat, d);
    let resp = upsample_rows(&column.half_range_response_k, column.nlat, d);
    let f_row = upsample_rows(&column.rainout_factor, column.nlat, d);
    let mut out = LocalRainout {
        precipitation_mm_per_yr: vec![0.0; n],
        orographic_mm_per_yr: vec![0.0; n],
        evaporation_mm_per_yr: vec![0.0; n],
        pet_mm_per_yr: vec![0.0; n],
    };
    for c in 0..n {
        let i = c / nlon;
        let t_mean = if ocean[c] {
            t_sl[i]
        } else {
            let surface = elevation_m[c].max(fill_level_m[c]);
            t_sl[i] - k.lapse * (surface - sea_level_m)
        };
        let m = if ocean[c] {
            k.damp
        } else {
            let t = tracer[c].clamp(0.0, 1.0);
            k.damp + (1.0 - k.damp) * (1.0 - t)
        };
        let half_range = m * resp[i];
        let pet = thornthwaite_pet_mm_per_yr(t_mean, half_range) / SECONDS_PER_YEAR;
        let w_sat = saturation_column_kg_m2(k.w_ref, k.t_ref, k.cc, t_mean);
        let w_sat_sl = saturation_column_kg_m2(k.w_ref, k.t_ref, k.cc, t_sl[i]);
        let wc = w[c];
        let pr = match k.crh {
            None => f_row[i] * wc / k.tau_r + (wc - w_sat).max(0.0) / k.tau_f,
            Some(q) => q.rate(f_row[i] / k.tau_r, wc, w_sat) + (wc - w_sat).max(0.0) / k.tau_f,
        };
        let p_oro = ((wc - w_sat).max(0.0) - (wc - w_sat_sl).max(0.0)).max(0.0) / k.tau_f;
        let e = if ocean[c] {
            0.0 // ocean evaporation is not needed for runoff
        } else {
            let ratio = if pr > 0.0 { budyko_ratio(pet / pr).min(k.et_cap) } else { k.et_cap };
            ratio * pr
        };
        out.precipitation_mm_per_yr[c] = pr * SECONDS_PER_YEAR;
        out.orographic_mm_per_yr[c] = p_oro * SECONDS_PER_YEAR;
        out.evaporation_mm_per_yr[c] = e * SECONDS_PER_YEAR;
        out.pet_mm_per_yr[c] = pet * SECONDS_PER_YEAR;
    }
    for (name, v) in [
        ("precipitation_mm_per_yr", &out.precipitation_mm_per_yr),
        ("evaporation_mm_per_yr", &out.evaporation_mm_per_yr),
        ("pet_mm_per_yr", &out.pet_mm_per_yr),
    ] {
        if let Some(c) = v.iter().position(|x| !x.is_finite()) {
            return Err(ClimateError::NonFinite { field: name, index: c });
        }
    }
    Ok(out)
}

/// Runoff (m/yr) from a local rainout: land max(P - E, f P), lake
/// max(P - PET, 0), ocean r_ref.
pub fn runoff_from_rainout(
    p: &GenesisClimateParams,
    r: &LocalRainout,
    elevation_m: &[f64],
    fill_level_m: &[f64],
    ocean: &[bool],
) -> Vec<f64> {
    (0..ocean.len())
        .map(|c| {
            if ocean[c] {
                p.r_ref_m_per_yr
            } else {
                let pr = r.precipitation_mm_per_yr[c];
                let mm = if fill_level_m[c] > elevation_m[c] + 1.0e-9 {
                    (pr - r.pet_mm_per_yr[c]).max(0.0)
                } else {
                    (pr - r.evaporation_mm_per_yr[c]).max(p.floor_fraction * pr)
                };
                mm / 1000.0
            }
        })
        .collect()
}

/// One climate solve on the solve grid, warm-started from `warm`.
pub fn solve_epoch(
    p: &GenesisClimateParams,
    grid: &GridSpec,
    elevation_m: &[f64],
    fill_level_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    warm: Option<&SeedColumnStateV1>,
) -> Result<(ClimateOutputV0, SeedColumnStateV1), ClimateError> {
    solve_epoch_with_retries(p, grid, elevation_m, fill_level_m, ocean, sea_level_m, warm).map(|(o, c, _)| (o, c))
}

/// `solve_epoch` with the I12 closure ladder; also returns the number of
/// retries used. Only an atmospheric-water closure miss is retried; any
/// other error (energy closure, non-finite field, no convergence) and a
/// miss at the floor tolerance stop loud as before.
#[allow(clippy::too_many_arguments)]
pub fn solve_epoch_with_retries(
    p: &GenesisClimateParams,
    grid: &GridSpec,
    elevation_m: &[f64],
    fill_level_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    warm: Option<&SeedColumnStateV1>,
) -> Result<(ClimateOutputV0, SeedColumnStateV1, u32), ClimateError> {
    let (cgrid, surf, coc) = coarsen_surface(grid, p.divisor, elevation_m, fill_level_m, ocean)?;
    let inputs = ClimateInputs { grid: &cgrid, elevation_m: &surf, fill_level_m: &surf, sea_level_m, ocean: &coc };
    let mut tol = p.sweep_rel_tol;
    let mut retries = 0u32;
    loop {
        let opts = SeedSolveOptions {
            initial_column_water_kg_m2: warm.map(|w| w.column_water_kg_m2.as_slice()),
            initial_tracer: warm.map(|w| w.tracer.as_slice()),
            sweep_rel_tol: tol,
        };
        let r = match &p.rotation {
            None => derive_seed_climate_v1_with(&inputs, &p.recipe, &p.seed, &opts),
            Some(rot) => derive_seed_climate_v2_with(&inputs, &p.recipe, &p.seed, rot, &opts),
        };
        match r {
            Err(ClimateError::ClosureFailed(c)) if c.quantity == "atmospheric-water" && tol > SEED_SWEEP_REL_TOL => {
                tol = (tol / GENESIS_CLOSURE_RETRY_TOL_FACTOR).max(SEED_SWEEP_REL_TOL);
                retries += 1;
            }
            Err(e) => return Err(e),
            Ok((o, c)) => return Ok((o, c, retries)),
        }
    }
}

/// The runoff the fluvial kernels receive at genesis step `step` (0-based
/// count of steps already taken), with the step's surface; solves when the
/// schedule says so and updates `state`. Returns the runoff (m/yr).
#[allow(clippy::too_many_arguments)]
pub fn genesis_runoff_step(
    p: &GenesisClimateParams,
    state: &mut GenesisClimateState,
    grid: &GridSpec,
    step: u32,
    t_my: f64,
    elevation_m: &[f64],
    fill_level_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    areas_m2: &[f64],
) -> Result<Vec<f64>, ClimateError> {
    let solved = p.solves_at(step);
    if solved {
        let t0 = std::time::Instant::now();
        let (out, col, closure_retries) =
            solve_epoch_with_retries(p, grid, elevation_m, fill_level_m, ocean, sea_level_m, state.column.as_ref())?;
        state.last_solve_wall_s = t0.elapsed().as_secs_f64();
        let water = &out.closures[0];
        state.epochs.push(GenesisClimateEpochLogV1 {
            step,
            t_my,
            tracer_sweeps: col.tracer_sweeps,
            moisture_sweeps: col.moisture_sweeps,
            water_residual_rel: if water.gross > 0.0 { water.residual / water.gross } else { 0.0 },
            land_mean_precipitation_mm_per_yr: out.summary.land_mean_precipitation_mm_per_yr,
            max_precipitation_mm_per_yr: out.summary.max_precipitation_mm_per_yr,
            land_orographic_share: out.summary.land_orographic_share,
            closure_retries,
        });
        state.column = Some(col);
    }
    let column = state
        .column
        .as_ref()
        .ok_or_else(|| ClimateError::Inputs(format!("{GENESIS_CLIMATE_LAW_ID}: no solved column at step {step} (state lost)")))?;
    let rain = local_rainout(p, column, grid, elevation_m, fill_level_m, ocean, sea_level_m)?;
    let mut runoff = match p.between {
        ClimateBetweenEpochsV1::FrozenColumnLocalRainoutV1 => runoff_from_rainout(p, &rain, elevation_m, fill_level_m, ocean),
        ClimateBetweenEpochsV1::HoldRunoffV1 => {
            if solved {
                let r = runoff_from_rainout(p, &rain, elevation_m, fill_level_m, ocean);
                state.held_runoff_m_per_yr = Some(r);
            }
            let held = state
                .held_runoff_m_per_yr
                .as_ref()
                .ok_or_else(|| ClimateError::Inputs(format!("{GENESIS_CLIMATE_LAW_ID}: no held runoff at step {step}")))?;
            // The geographic-frame field; cells that became ocean read r_ref.
            (0..ocean.len()).map(|c| if ocean[c] { p.r_ref_m_per_yr } else { held[c] }).collect()
        }
    };
    let land_mean_pow = |r: &[f64]| -> f64 {
        let (mut s, mut a) = (0.0, 0.0);
        for c in 0..r.len() {
            if !ocean[c] {
                s += areas_m2[c] * (r[c] / p.r_ref_m_per_yr).powf(p.m);
                a += areas_m2[c];
            }
        }
        if a > 0.0 {
            s / a
        } else {
            1.0
        }
    };
    match p.erosivity {
        ClimateErosivityModeV1::ClimateRunoffPatternV1 => {}
        ClimateErosivityModeV1::UniformReferenceDiag => {
            runoff.iter_mut().for_each(|r| *r = p.r_ref_m_per_yr);
        }
        ClimateErosivityModeV1::LandMeanRunoffDiag => {
            let (mut s, mut a) = (0.0, 0.0);
            for c in 0..runoff.len() {
                if !ocean[c] {
                    s += areas_m2[c] * runoff[c];
                    a += areas_m2[c];
                }
            }
            if a > 0.0 {
                let u = s / a;
                for c in 0..runoff.len() {
                    if !ocean[c] {
                        runoff[c] = u;
                    }
                }
            }
        }
        ClimateErosivityModeV1::LandMeanEquivalentDiag => {
            let u = land_mean_pow(&runoff).powf(1.0 / p.m) * p.r_ref_m_per_yr;
            for c in 0..runoff.len() {
                if !ocean[c] {
                    runoff[c] = u;
                }
            }
        }
    }
    let (mut sr, mut a, mut mx, mut sp, mut so) = (0.0, 0.0, 0.0f64, 0.0, 0.0);
    for c in 0..runoff.len() {
        if !ocean[c] {
            sr += areas_m2[c] * runoff[c];
            a += areas_m2[c];
            mx = mx.max(runoff[c]);
            sp += areas_m2[c] * rain.precipitation_mm_per_yr[c];
            so += areas_m2[c] * rain.orographic_mm_per_yr[c];
        }
    }
    state.steps.push(GenesisClimateStepLogV1 {
        step,
        t_my,
        solved,
        land_mean_runoff_m_per_yr: if a > 0.0 { sr / a } else { 0.0 },
        land_mean_erosivity_factor: land_mean_pow(&runoff),
        max_land_runoff_m_per_yr: mx,
        land_orographic_share: if sp > 0.0 { so / sp } else { 0.0 },
    });
    state.last_runoff_m_per_yr = runoff.clone();
    state.last_precipitation_mm_per_yr = rain.precipitation_mm_per_yr;
    state.last_orographic_mm_per_yr = rain.orographic_mm_per_yr;
    state.last_pet_mm_per_yr = rain.pet_mm_per_yr;
    Ok(runoff)
}

#[cfg(test)]
mod tests {
    //! Every map here is a SYNTHETIC FIXTURE (ADR 026 §1: legal as
    //! instrument calibration in unit tests, never world input).
    use super::*;
    use pf_contracts::laws::{ConstantV1, SourceClassV1, UnitV1};
    use std::borrow::Cow;

    fn law(b: &ClimateForcingBlock) -> ClimateForcingLawV1 {
        let c = |v, u| ConstantV1::new(v, u, SourceClassV1::ProjectAssumptionAwaitingReview { basis: Cow::Borrowed("test") });
        ClimateForcingLawV1::EpochClimateRunoffV1 {
            climate: b.climate.clone(),
            seed_constants_digest: b.seed_constants_digest.clone(),
            rotation_constants_digest: b.rotation_constants_digest.clone(),
            climate_epoch_my: c(b.epoch_my, UnitV1::Megayear),
            climate_solve_grid_divisor: c(b.solve_grid_divisor, UnitV1::Dimensionless),
            climate_sweep_rel_tol: c(b.sweep_rel_tol, UnitV1::Dimensionless),
            climate_runoff_floor_fraction: c(b.runoff_floor_fraction, UnitV1::Dimensionless),
            between_epochs: b.between_epochs,
            erosivity: b.erosivity,
            subgrid_cell_mean: b.subgrid_cell_mean,
            erosive_runoff: None,
        }
    }

    /// HYPSO-LEVERS (2026-10-21): the erosive-runoff response. Absent =
    /// None (the caller keeps the modelled runoff); the sublinear law leaves
    /// runoff at or below r_ref unchanged and is continuous there; the cap
    /// clips.
    #[test]
    fn erosive_runoff_response_sublinear_above_reference_and_cap() {
        use pf_contracts::laws::ErosiveRunoffLawV1;
        let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
        let b = climate_forcing_block_v1(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        let resolve = |l: &ClimateForcingLawV1| GenesisClimateParams::from_law(l, 2.0, &grid, 1.0, 1.0, 9.81).unwrap().unwrap();
        let r = [0.0, 0.2, 1.0, 1.0 + 1e-12, 4.0, 9.0];
        assert!(resolve(&law(&b)).erosive_runoff_m_per_yr(&r).is_none());
        let c = |v, u| ConstantV1::new(v, u, SourceClassV1::ProjectAssumptionAwaitingReview { basis: Cow::Borrowed("test") });
        let mut l = law(&b);
        if let ClimateForcingLawV1::EpochClimateRunoffV1 { erosive_runoff, .. } = &mut l {
            *erosive_runoff = Some(ErosiveRunoffLawV1::SublinearAboveReferenceV1 { erosive_runoff_exponent: c(0.5, UnitV1::Dimensionless) });
        }
        let e = resolve(&l).erosive_runoff_m_per_yr(&r).unwrap();
        assert_eq!(&e[..3], &r[..3]);
        assert!((e[3] - 1.0).abs() < 1e-11 && (e[4] - 2.0).abs() < 1e-12 && (e[5] - 3.0).abs() < 1e-12, "{e:?}");
        if let ClimateForcingLawV1::EpochClimateRunoffV1 { erosive_runoff, .. } = &mut l {
            *erosive_runoff = Some(ErosiveRunoffLawV1::CapV1 { erosive_runoff_cap_m_per_yr: c(1.5, UnitV1::MeterPerYear) });
        }
        assert_eq!(resolve(&l).erosive_runoff_m_per_yr(&r).unwrap(), vec![0.0, 0.2, 1.0, 1.0 + 1e-12, 1.5, 1.5]);
    }

    /// A continent with a north-south range across the westerlies.
    fn fixture(nlat: usize, nlon: usize) -> (GridSpec, Vec<f64>, Vec<bool>) {
        let grid = GridSpec::new(nlat, nlon, 6.371e6).unwrap();
        let mut z = vec![-3000.0; grid.cells()];
        for i in 0..nlat {
            for j in 0..nlon {
                let lat = grid.center_lat_deg(i).unwrap();
                if lat.abs() < 60.0 && j >= nlon / 4 && j < nlon / 2 {
                    let x = (j - nlon / 4) as f64 / (nlon / 4) as f64;
                    z[i * nlon + j] = 200.0 + 3000.0 * (-((x - 0.3) / 0.08).powi(2)).exp();
                }
            }
        }
        let ocean: Vec<bool> = z.iter().map(|v| *v < 0.0).collect();
        (grid, z, ocean)
    }

    #[test]
    fn local_rainout_reproduces_the_solver_at_divisor_one() {
        let (grid, z, ocean) = fixture(24, 48);
        let fill: Vec<f64> = z.iter().map(|v| v.max(0.0)).collect();
        let mut b = climate_forcing_block_v1(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        b.solve_grid_divisor = 1.0;
        b.epoch_my = 2.0;
        b.sweep_rel_tol = 1.0e-11;
        let p = GenesisClimateParams::from_law(&law(&b), 2.0, &grid, 1.0, 0.5, 9.81).unwrap().unwrap();
        let (out, col) = solve_epoch(&p, &grid, &z, &fill, &ocean, 0.0, None).unwrap();
        // The divisor-1 solve with default options is the post-hoc tier.
        let inputs = ClimateInputs { grid: &grid, elevation_m: &z.iter().zip(&fill).map(|(a, b)| a.max(*b)).collect::<Vec<_>>(), fill_level_m: &fill, sea_level_m: 0.0, ocean: &ocean };
        let post = crate::climate::derive_climate(&inputs, &p.recipe).unwrap();
        assert_eq!(post.fields, out.fields, "divisor 1, tol 1e-11, cold start == post-hoc seed-climate-v1");
        let lr = local_rainout(&p, &col, &grid, &z, &fill, &ocean, 0.0).unwrap();
        for c in 0..grid.cells() {
            assert_eq!(lr.precipitation_mm_per_yr[c].to_bits(), out.fields.precipitation_mm_per_yr[c].to_bits(), "P at {c}");
            assert_eq!(lr.orographic_mm_per_yr[c].to_bits(), out.fields.precipitation_orographic_mm_per_yr[c].to_bits(), "P_oro at {c}");
            assert_eq!(lr.pet_mm_per_yr[c].to_bits(), out.fields.pet_mm_per_yr[c].to_bits(), "PET at {c}");
            if !ocean[c] {
                assert_eq!(lr.evaporation_mm_per_yr[c].to_bits(), out.fields.evaporation_mm_per_yr[c].to_bits(), "E at {c}");
            }
        }
    }

    /// CLIMATE-ARIDITY: the earth_crh_v1 digest selects the crh-rainout-v1
    /// set, and the pointwise rainout reproduces its solver bit for bit.
    #[test]
    fn crh_set_selected_by_digest_and_local_rainout_reproduces_the_solver() {
        let (grid, z, ocean) = fixture(24, 48);
        let fill: Vec<f64> = z.iter().map(|v| v.max(0.0)).collect();
        let mut b = climate_forcing_block_crh_v1(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        b.solve_grid_divisor = 1.0;
        b.epoch_my = 2.0;
        b.sweep_rel_tol = 1.0e-11;
        let p = GenesisClimateParams::from_law(&law(&b), 2.0, &grid, 1.0, 0.5, 9.81).unwrap().unwrap();
        assert_eq!(p.seed, SeedClimateConstantsV1::earth_crh_v1());
        let (out, col) = solve_epoch(&p, &grid, &z, &fill, &ocean, 0.0, None).unwrap();
        let surf: Vec<f64> = z.iter().zip(&fill).map(|(a, b)| a.max(*b)).collect();
        let inputs = ClimateInputs { grid: &grid, elevation_m: &surf, fill_level_m: &fill, sea_level_m: 0.0, ocean: &ocean };
        let post = crate::climate::seed::derive_seed_climate_v1(&inputs, &p.recipe, &SeedClimateConstantsV1::earth_crh_v1()).unwrap();
        assert_eq!(post.fields, out.fields, "divisor 1, tol 1e-11, cold start == post-hoc crh run");
        let v1 = crate::climate::derive_climate(&inputs, &p.recipe).unwrap();
        assert_ne!(v1.fields.precipitation_mm_per_yr, out.fields.precipitation_mm_per_yr, "the option changes P");
        let lr = local_rainout(&p, &col, &grid, &z, &fill, &ocean, 0.0).unwrap();
        for c in 0..grid.cells() {
            assert_eq!(lr.precipitation_mm_per_yr[c].to_bits(), out.fields.precipitation_mm_per_yr[c].to_bits(), "P at {c}");
            if !ocean[c] {
                assert_eq!(lr.evaporation_mm_per_yr[c].to_bits(), out.fields.evaporation_mm_per_yr[c].to_bits(), "E at {c}");
            }
        }
        // The v1 block still selects the v1 set.
        let b1 = climate_forcing_block_v1(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        let p1 = GenesisClimateParams::from_law(&law(&b1), 2.0, &grid, 1.0, 0.5, 9.81).unwrap().unwrap();
        assert_eq!(p1.seed, SeedClimateConstantsV1::earth_v1());
    }

    #[test]
    fn coarse_warm_solve_closes_water_and_keeps_the_rain_shadow() {
        let (grid, z, ocean) = fixture(48, 96);
        let fill: Vec<f64> = z.iter().map(|v| v.max(0.0)).collect();
        let b = climate_forcing_block_v1(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        let p = GenesisClimateParams::from_law(&law(&b), 2.0, &grid, 1.0, 0.5, 9.81).unwrap().unwrap();
        assert_eq!((p.epoch_steps, p.divisor), (2, 2));
        let mut st = GenesisClimateState::default();
        let areas: Vec<f64> = (0..grid.cells()).map(|c| grid.cell_area_m2(c / grid.nlon()).unwrap()).collect();
        let r0 = genesis_runoff_step(&p, &mut st, &grid, 0, 0.0, &z, &fill, &ocean, 0.0, &areas).unwrap();
        let cold = st.epochs[0].moisture_sweeps;
        assert!(st.epochs[0].water_residual_rel.abs() <= 1e-6);
        // Between epochs: no solve, same surface -> same runoff.
        let r1 = genesis_runoff_step(&p, &mut st, &grid, 1, 2.0, &z, &fill, &ocean, 0.0, &areas).unwrap();
        assert_eq!(r0, r1);
        assert_eq!(st.epochs.len(), 1);
        // Warm restart on an unchanged surface converges in fewer sweeps.
        genesis_runoff_step(&p, &mut st, &grid, 2, 4.0, &z, &fill, &ocean, 0.0, &areas).unwrap();
        assert!(st.epochs[1].moisture_sweeps < cold, "warm {} vs cold {cold}", st.epochs[1].moisture_sweeps);
        // Westerlies at 45 deg: the windward (west) flank of the range is
        // wetter than the lee (east) flank.
        let nlon = grid.nlon();
        let i = (0..grid.nlat()).find(|&i| (grid.center_lat_deg(i).unwrap() - 45.0).abs() < 2.0).unwrap();
        let crest = (nlon / 4..nlon / 2).max_by(|a, b| z[i * nlon + a].total_cmp(&z[i * nlon + b])).unwrap();
        let west = r0[i * nlon + crest - 2];
        let east = r0[i * nlon + crest + 3];
        assert!(west > 1.5 * east, "windward {west} vs lee {east}");
        // Diagnostics: uniform mode is r_ref everywhere; land-mean mode keeps
        // the mean erosivity factor.
        let mut pu = p.clone();
        pu.erosivity = ClimateErosivityModeV1::UniformReferenceDiag;
        let ru = genesis_runoff_step(&pu, &mut GenesisClimateState::default(), &grid, 0, 0.0, &z, &fill, &ocean, 0.0, &areas).unwrap();
        assert!(ru.iter().all(|r| *r == 1.0));
        let mut pl = p.clone();
        pl.erosivity = ClimateErosivityModeV1::LandMeanEquivalentDiag;
        let mut sl = GenesisClimateState::default();
        genesis_runoff_step(&pl, &mut sl, &grid, 0, 0.0, &z, &fill, &ocean, 0.0, &areas).unwrap();
        let want = st.steps[0].land_mean_erosivity_factor;
        assert!((sl.steps[0].land_mean_erosivity_factor - want).abs() <= 1e-12 * want);
        // Land-mean-runoff diagnostic keeps the land mean R.
        let mut pr = p.clone();
        pr.erosivity = ClimateErosivityModeV1::LandMeanRunoffDiag;
        let mut sr = GenesisClimateState::default();
        genesis_runoff_step(&pr, &mut sr, &grid, 0, 0.0, &z, &fill, &ocean, 0.0, &areas).unwrap();
        let want_r = st.steps[0].land_mean_runoff_m_per_yr;
        assert!((sr.steps[0].land_mean_runoff_m_per_yr - want_r).abs() <= 1e-12 * want_r);
    }

    /// I12: a solve that misses the water closure at a loose tolerance is
    /// re-run tighter and closes; a solve that closes at once is unchanged.
    #[test]
    fn closure_ladder_retries_only_a_water_miss() {
        let (grid, z, ocean) = fixture(48, 96);
        let fill: Vec<f64> = z.iter().map(|v| v.max(0.0)).collect();
        let mut b = climate_forcing_block_v1(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        b.solve_grid_divisor = 1.0;
        // The envelope's loosest tolerance: the cold solve stops far from
        // the water closure.
        b.sweep_rel_tol = 1.0e-5;
        let p = GenesisClimateParams::from_law(&law(&b), 2.0, &grid, 1.0, 0.5, 9.81).unwrap().unwrap();
        let cg = crate::climate::ClimateInputs { grid: &grid, elevation_m: &z.iter().zip(&fill).map(|(a, b)| a.max(*b)).collect::<Vec<_>>(), fill_level_m: &fill, sea_level_m: 0.0, ocean: &ocean };
        let once = derive_seed_climate_v1_with(&cg, &p.recipe, &p.seed, &SeedSolveOptions { sweep_rel_tol: 1.0e-5, ..Default::default() });
        assert!(matches!(once, Err(ClimateError::ClosureFailed(_))), "fixture must miss the closure at 1e-5: {once:?}");
        let (out, _, retries) = solve_epoch_with_retries(&p, &grid, &z, &fill, &ocean, 0.0, None).unwrap();
        assert!(retries >= 1);
        assert!(out.closures[0].closes());
        // Default tolerance on the same fixture: no retry, identical to the
        // single solve.
        let mut q = p.clone();
        q.sweep_rel_tol = 1.0e-9;
        let (o9, c9, r9) = solve_epoch_with_retries(&q, &grid, &z, &fill, &ocean, 0.0, None).unwrap();
        assert_eq!(r9, 0);
        let (o9b, c9b) = derive_seed_climate_v1_with(&cg, &q.recipe, &q.seed, &SeedSolveOptions { sweep_rel_tol: 1.0e-9, ..Default::default() }).unwrap();
        assert_eq!((o9.fields, c9), (o9b.fields, c9b));
    }

    #[test]
    fn refusals_are_loud() {
        let grid = GridSpec::new(24, 48, 6.371e6).unwrap();
        let b = climate_forcing_block_v1(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        // Epoch not a whole number of steps.
        let mut e = b.clone();
        e.epoch_my = 3.0;
        assert!(GenesisClimateParams::from_law(&law(&e), 2.0, &grid, 1.0, 0.5, 9.81).is_err());
        // Divisor not tiling (4 does not tile 26 rows).
        let g = GridSpec::new(26, 48, 6.371e6).unwrap();
        let mut q = b.clone();
        q.solve_grid_divisor = 4.0;
        assert!(GenesisClimateParams::from_law(&law(&q), 2.0, &g, 1.0, 0.5, 9.81).is_err());
        // Foreign constants digest.
        let mut d = b.clone();
        d.seed_constants_digest = "blake3:00".into();
        assert!(GenesisClimateParams::from_law(&law(&d), 2.0, &grid, 1.0, 0.5, 9.81).is_err());
        // v2 with a gravity unequal to the universe's.
        let v2 = climate_forcing_block_v2(ClimateErosivityModeV1::ClimateRunoffPatternV1).unwrap();
        assert!(GenesisClimateParams::from_law(&law(&v2), 2.0, &grid, 1.0, 0.5, 9.81).unwrap().is_some());
        assert!(GenesisClimateParams::from_law(&law(&v2), 2.0, &grid, 1.0, 0.5, 9.0).is_err());
        assert!(GenesisClimateParams::from_law(&ClimateForcingLawV1::NoneV0, 2.0, &grid, 1.0, 0.5, 9.81).unwrap().is_none());
    }
}
