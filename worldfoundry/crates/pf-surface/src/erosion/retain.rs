//! SEED-SEDIMENT-RETAIN (lane SR, 2026-10-13): law id
//! `continental-retention-base-level-v1` (law family `sediment_retention`),
//! an opt-in replacement of the legacy erosion step (`erode_step_bounded`)
//! for worlds whose land laws are otherwise legacy (detachment-limited land,
//! local fill + global layer). Nothing here runs unless the law is
//! selected; the legacy step in `erosion.rs` is untouched.
//!
//! DEFECT IT ADDRESSES (lane SC report, design review review 2026-10-13): over long
//! runs continents drown because (i) eroded continental rock leaves the
//! continental cells as sediment on OCEANIC cells, where the carry later
//! subducts it, and (ii) the implicit n = 1 step lowers a low-gradient
//! lowland cell nearly to its receiver every step (f = K A^m dt / L >> 1),
//! so land near base level keeps eroding at tens of m/My. On Earth most
//! sediment is stored on land and on continental margins, and low-relief
//! land near base level erodes at a few m/My.
//!
//! THE LAW (constants in `RetentionParams::v1`, pre-registered in
//! `unpublished-design-record` before any run):
//! 1. LOW-RELIEF, NEAR-BASE-LEVEL TAPER (rule 3, 2026-10-13). The
//!    stream-power lowering of a land cell in a step is multiplied by
//!    phi = min(1, (h / H)^p), h = max(r, b): r is the cell's LOCAL RELIEF
//!    (max minus min of the start-of-step level over the cell and its 8
//!    neighbours, water cells at their water level) and b its HEIGHT ABOVE
//!    BASE LEVEL (start-of-step elevation minus the water level its flow
//!    path ends in). Erosion slows only where a cell is both low-relief and
//!    near base level (denudation rises with relief and with elevation;
//!    Ahnert 1970, Montgomery and Brandon 2002, recalled): planed lowlands
//!    slow; ranges, escarpments, a coastal cell below a plateau edge (high
//!    r) and a plateau interior >= H above the sea (high b) erode exactly as
//!    in the legacy step. Rule 1 (b only) made coastal cells dams that held
//!    every cell upstream (no cell is lowered below its receiver); rule 2
//!    (r only) froze flat high plateaus (lane SR PREREG.txt). The lowering
//!    stays between 0 and the n = 1 implicit lowering, so no cell goes below
//!    its receiver. Hillslope diffusion is the legacy scheme, unchanged.
//! 2. CONTINENTAL-FIRST DEPOSITION. Material detached on land travels down
//!    the flow path (a donor-first pass) and is stored, in this order:
//!    a. on CONTINENTAL land cells standing below the graded alluvial level
//!       g = receiver level + S_a L (L the flow length to the receiver), up
//!       to g, with the compensated thickness (`CompensatedColumn`); a
//!       cell's own detached material passes on;
//!    b. in lakes, filled to the spill level with compensated thickness
//!       (bounded breadth-first spread as in the legacy step); a full lake
//!       passes the rest downstream;
//!    c. at the river mouth, on flooded CONTINENTAL cells, filled to wave
//!       base (sea - wave base) with compensated thickness, nearest first
//!       (Dijkstra by centre distance) through flooded continental cells
//!       only, seeded at the mouth cell if it is continental, else at its
//!       flooded continental neighbours;
//!    d. only what exceeds a-c reaches the legacy placement (fill to sea
//!       level from the mouth, then the declared global layer), i.e.
//!       oceanic cells.
//! MASS: every change is a raise booked once (`deposit`) or a lowering
//! (`Lowering::lower`), so eroded = deposited + exported to rounding.
//! NOT MODELLED (declared): flexural foreland and cratonic-sag subsidence
//! (local Airy only, so land accommodation is only the graded plain and
//! lakes), deposition on land between a lake outlet and the next water
//! body (as legacy), grain sizes, compaction, conversion of thick margin
//! sediment wedges into continental crust (kinds belong to the carry).

use super::{
    deposit, place_legacy, DiffusionStencil, ErosionError, ErosionParams, ErosionStep, LayerTally,
    Lowering, OceanRemainderLayer, Sink, MAX_SPREAD_CELLS,
};
use crate::routing::Routing;
use crate::sediment::{
    CompensatedColumn, DepositionReport, MarineClass, MARINE_SEARCH_MAX_CELLS,
    SHELF_WAVE_BASE_DEPTH_M,
};
use pf_contracts::spatial::ElementId;
use pf_grid::GridSpec;

/// Law id (the typed law set lives in pf-contracts laws.rs).
pub const SEDIMENT_RETENTION_LAW_ID: &str = "continental-retention-base-level-v1";
/// H: max(local relief, height above base level) at and above which the
/// taper is inactive (m).
/// Pre-registered 2026-10-13 09:47 EDT; project assumption (PREREG.txt).
pub const BASE_LEVEL_RELIEF_SCALE_M: f64 = 1000.0;
/// p: taper exponent. Pre-registered; project assumption.
pub const BASE_LEVEL_TAPER_EXPONENT: f64 = 2.0;
/// S_a: graded alluvial slope of cell means. Pre-registered; project
/// assumption (lowland alluvial gradients ~2e-5 to ~2e-4, recalled).
pub const ALLUVIAL_GRADE_SLOPE: f64 = 1.0e-4;

/// Parameters of `continental-retention-base-level-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RetentionParams {
    pub base_level_relief_scale_m: f64,
    pub taper_exponent: f64,
    pub alluvial_grade_slope: f64,
    pub shelf_wave_base_depth_m: f64,
    /// Computational bound on flooded continental cells settled per mouth.
    pub search_max_cells: usize,
    /// The isostasy law's densities (mantle, sediment, sea water).
    pub column: CompensatedColumn,
    /// Lane INCISION-THRESHOLD part 1 (2026-10-20): how (d)'s marine
    /// remainder is laid; `v1()` keeps the legacy `PerMouthV0`.
    pub ocean_layer: OceanRemainderLayer,
}

impl RetentionParams {
    /// The pre-registered values with the caller's isostasy densities.
    pub fn v1(column: CompensatedColumn) -> Self {
        Self {
            base_level_relief_scale_m: BASE_LEVEL_RELIEF_SCALE_M,
            taper_exponent: BASE_LEVEL_TAPER_EXPONENT,
            alluvial_grade_slope: ALLUVIAL_GRADE_SLOPE,
            shelf_wave_base_depth_m: SHELF_WAVE_BASE_DEPTH_M,
            search_max_cells: MARINE_SEARCH_MAX_CELLS,
            column,
            ocean_layer: OceanRemainderLayer::PerMouthV0,
        }
    }

    pub fn validate(&self) -> Result<(), ErosionError> {
        self.column
            .validate()
            .map_err(|e| ErosionError::Law(format!("{SEDIMENT_RETENTION_LAW_ID}: {e}")))?;
        for (what, v) in [
            ("base_level_relief_scale_m", self.base_level_relief_scale_m),
            ("taper_exponent", self.taper_exponent),
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(ErosionError::BadParam { what, value: v });
            }
        }
        for (what, v) in [
            ("alluvial_grade_slope", self.alluvial_grade_slope),
            ("shelf_wave_base_depth_m", self.shelf_wave_base_depth_m),
        ] {
            if !v.is_finite() || v < 0.0 {
                return Err(ErosionError::BadParam { what, value: v });
            }
        }
        if self.search_max_cells == 0 {
            return Err(ErosionError::BadParam {
                what: "search_max_cells",
                value: 0.0,
            });
        }
        Ok(())
    }

    /// phi(h) = min(1, (h / H)^p) for h = max(local relief, height above
    /// base level) >= 0.
    pub fn taper(&self, height_above_base_m: f64) -> f64 {
        let h = height_above_base_m.max(0.0);
        if h >= self.base_level_relief_scale_m {
            return 1.0;
        }
        (h / self.base_level_relief_scale_m)
            .powf(self.taper_exponent)
            .min(1.0)
    }
}

/// Where this step's material went (m^3) and what the taper withheld.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RetentionReport {
    /// (a) graded-plain storage on continental land cells.
    pub continental_grade_m3: f64,
    /// (b) lake storage (compensated).
    pub lake_m3: f64,
    /// Terminal floors with no outlet (kept where they arrive).
    pub land_floor_m3: f64,
    /// (c) flooded continental cells filled toward wave base.
    pub continental_shelf_m3: f64,
    /// Q25: conserved terrigenous fines, separate from shelf/legacy placement.
    pub marine_fines_m3: f64,
    /// (d) handed to the legacy placement (fill to sea level from the
    /// mouth, then the global layer); `global_layer_m3` is part of it.
    pub legacy_marine_m3: f64,
    pub global_layer_m3: f64,
    /// Hillslope diffusion into standing water (legacy one-way sink).
    pub hillslope_to_water_m3: f64,
    /// All sediment added this step, split by the kind of receiving cell.
    pub added_on_continental_m3: f64,
    pub added_on_oceanic_m3: f64,
    /// Stream-power lowering the taper withheld: sum of (1 - phi) times the
    /// untapered implicit lowering (diagnostic; not a mass term).
    pub taper_withheld_m3: f64,
    pub tapered_cells: u64,
    /// Ocean deliveries with no flooded continental cell at the mouth.
    pub mouths_without_continental_shelf: u64,
    pub mouths: u64,
    pub shelf_search_bound_hits: u64,
    /// DIAGNOSTIC (lane INCISION-THRESHOLD part 1): mouths whose remainder
    /// reached the global layer, ocean-cell raises by that layer (cost
    /// proxy), and cells it lifted from at or below sea level to above it.
    pub global_layer_mouths: u64,
    pub global_layer_cell_raises: u64,
    pub global_layer_cells_above_sea: u64,
    /// DIAGNOSTIC (lane LAKE-EXCESS, `filled-basin-alluvial-grade-v1`):
    /// lake bodies treated as filled basins this step and the volume their
    /// graded fill placed (part of `lake_m3`).
    pub graded_basins: u64,
    pub graded_basin_m3: f64,
    /// DIAGNOSTIC (lane LAKE-OUTLET, `lake-outlet-sill-incision-v1`): pour
    /// cells of full lakes the outlet law evaluated, those where it set the
    /// lowering, and the lowering it added beyond the step's own (m^3).
    pub outlet_sills: u64,
    pub outlet_sills_incised: u64,
    pub outlet_extra_m3: f64,
    /// Sills already within 1 m of their receiver's level (nothing left to
    /// cut at the sill: the downstream profile holds the lake), and sills
    /// higher than that the channel law left unlowered (threshold).
    pub outlet_sills_level: u64,
    pub outlet_sills_below_threshold: u64,
    /// DIAGNOSTIC (lane SILL-REBOUND, rule 6): material removed at sills
    /// beyond their surface lowering (m^3; the share local rebound would
    /// have returned).
    pub outlet_support_extra_m3: f64,
}

impl RetentionReport {
    /// Add another report's terms (for interval and run totals).
    pub fn accumulate(&mut self, o: &RetentionReport) {
        self.continental_grade_m3 += o.continental_grade_m3;
        self.lake_m3 += o.lake_m3;
        self.land_floor_m3 += o.land_floor_m3;
        self.continental_shelf_m3 += o.continental_shelf_m3;
        self.marine_fines_m3 += o.marine_fines_m3;
        self.legacy_marine_m3 += o.legacy_marine_m3;
        self.global_layer_m3 += o.global_layer_m3;
        self.hillslope_to_water_m3 += o.hillslope_to_water_m3;
        self.added_on_continental_m3 += o.added_on_continental_m3;
        self.added_on_oceanic_m3 += o.added_on_oceanic_m3;
        self.taper_withheld_m3 += o.taper_withheld_m3;
        self.tapered_cells += o.tapered_cells;
        self.mouths_without_continental_shelf += o.mouths_without_continental_shelf;
        self.mouths += o.mouths;
        self.shelf_search_bound_hits += o.shelf_search_bound_hits;
        self.global_layer_mouths += o.global_layer_mouths;
        self.global_layer_cell_raises += o.global_layer_cell_raises;
        self.global_layer_cells_above_sea += o.global_layer_cells_above_sea;
        self.graded_basins += o.graded_basins;
        self.graded_basin_m3 += o.graded_basin_m3;
        self.outlet_sills += o.outlet_sills;
        self.outlet_sills_incised += o.outlet_sills_incised;
        self.outlet_extra_m3 += o.outlet_extra_m3;
        self.outlet_sills_level += o.outlet_sills_level;
        self.outlet_sills_below_threshold += o.outlet_sills_below_threshold;
        self.outlet_support_extra_m3 += o.outlet_support_extra_m3;
    }
}

thread_local! {
    /// DIAGNOSTIC ONLY (not run state, never in an output or a digest): the
    /// sum of every step report computed on this thread since the last
    /// `take_thread_totals`, so a measurement driver can read the partition
    /// without the run carrying it.
    static THREAD_TOTALS: std::cell::RefCell<RetentionReport> =
        std::cell::RefCell::new(RetentionReport::default());
}

/// Diagnostic: return and reset this thread's accumulated step reports.
pub fn take_thread_totals() -> RetentionReport {
    THREAD_TOTALS.with(|t| std::mem::take(&mut *t.borrow_mut()))
}

/// One erosion step under `continental-retention-base-level-v1`. Inputs and
/// the returned `ErosionStep` mean what they mean for
/// `erosion::erode_step_bounded` (`dz_m` is the column change before
/// isostasy; `deposition` is `None`, as on the legacy path, so module
/// checkpoints keep their legacy layout). `continental[c]` is true when
/// cell `c` carries continental crust.
#[allow(clippy::too_many_arguments)]
pub fn erode_retain_step(
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
    continental: &[bool],
    rp: &RetentionParams,
) -> Result<(ErosionStep, RetentionReport), ErosionError> {
    let (step,report,_)=erode_retain_step_recorded(grid,stencil,routing,elevation_m,sediment_m,ocean,sea_level_m,dt_my,params,bedrock_available_m,continental,rp,None)?;
    Ok((step,report))
}

#[allow(clippy::too_many_arguments)]
pub fn erode_retain_step_recorded(
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
    continental: &[bool],
    rp: &RetentionParams,
    trace_clock:Option<super::fluvial_trace::TraceClock>,
) -> Result<(ErosionStep, RetentionReport,Option<super::fluvial_trace::FluvialFrame>), ErosionError> {
    erode_retain_step_with_budgets(grid,stencil,routing,elevation_m,sediment_m,ocean,sea_level_m,dt_my,params,bedrock_available_m,continental,rp,trace_clock,false)
}

#[allow(clippy::too_many_arguments)]
pub fn erode_retain_step_with_budgets(
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
    continental: &[bool],
    rp: &RetentionParams,
    trace_clock:Option<super::fluvial_trace::TraceClock>,
    stage_budgets:bool,
) -> Result<(ErosionStep, RetentionReport,Option<super::fluvial_trace::FluvialFrame>), ErosionError> {
    erode_retain_step_with_lakes(grid,stencil,routing,elevation_m,sediment_m,ocean,sea_level_m,dt_my,params,bedrock_available_m,continental,rp,trace_clock,stage_budgets,None)
}

/// Lane INLAND-BASINS (2026-10-21): the retention step with the opt-in
/// `closed-basin-water-balance-v1` lake placement (closed_basin.rs). `lakes`
/// = None is the retention law exactly (the code path is unchanged).
#[allow(clippy::too_many_arguments)]
pub fn erode_retain_step_with_lakes(
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
    continental: &[bool],
    rp: &RetentionParams,
    trace_clock:Option<super::fluvial_trace::TraceClock>,
    stage_budgets:bool,
    lakes: Option<(&crate::closed_basin::LakeBalance, &crate::closed_basin::ClosedBasinParams)>,
) -> Result<(ErosionStep, RetentionReport,Option<super::fluvial_trace::FluvialFrame>), ErosionError> {
    let (step, report, trace, _) = erode_retain_step_with_marine_fines(
        grid, stencil, routing, elevation_m, sediment_m, ocean, sea_level_m,
        dt_my, params, bedrock_available_m, continental, rp, trace_clock,
        stage_budgets, lakes, None,
    )?;
    Ok((step, report, trace))
}

/// Q25 / D-0155: optional conserved terrigenous bypass. Existing entry points
/// remain disabled. This API does not add pelagic production or recipe wiring.
/// The fourth return value records the component-wise perfectly mixed sources
/// and destinations. Keep it with the step; it is not a cell-flow receipt.
#[allow(clippy::too_many_arguments)]
pub fn erode_retain_step_with_marine_fines(
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
    continental: &[bool],
    rp: &RetentionParams,
    trace_clock:Option<super::fluvial_trace::TraceClock>,
    stage_budgets:bool,
    lakes: Option<(&crate::closed_basin::LakeBalance, &crate::closed_basin::ClosedBasinParams)>,
    marine_fines: Option<&super::marine_fines::MarineFinesParams>,
) -> Result<(ErosionStep, RetentionReport,Option<super::fluvial_trace::FluvialFrame>, super::marine_fines::MarineFinesReport), ErosionError> {
    params.validate()?;
    rp.validate()?;
    if let Some(fines) = marine_fines { fines.validate()?; }
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
        ("continental", continental.len()),
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
    let mut rep = DepositionReport::default();
    let mut rr = RetentionReport::default();
    let mut transfers = std::collections::BTreeMap::<(usize, usize, &str), f64>::new();
    let col = &rp.column;

    // ---- Receivers and flow lengths (as the legacy step) ----------------
    let is_lake = |c: usize| !ocean[c] && routing.fill_level_m[c] > elevation_m[c] + 1.0e-9;
    let wet: Vec<bool> = (0..n).map(|i| ocean[i] || is_lake(i)).collect();
    let order = &routing.receiver_first;
    let mut recv: Vec<Option<(usize, f64)>> = vec![None; n];
    // The water body (or terminal floor) each cell's flow path ends in.
    let mut sink: Vec<usize> = (0..n).collect();
    for &i in order {
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
            recv[i] = Some((j, len));
            sink[i] = if wet[j] { j } else { sink[j] };
        }
    }
    let base_level = |t: usize| -> f64 {
        if ocean[t] {
            sea_level_m
        } else if is_lake(t) {
            routing.fill_level_m[t]
        } else {
            elevation_m[t]
        }
    };
    let relief = local_relief(grid, routing, elevation_m, ocean, sea_level_m)?;
    // LAKE-OUTLET (`lake-outlet-sill-incision-v1`, closed_basin.rs rule 5):
    // outlet discharge area at full lakes' sill cells; None = law off (the
    // step below is then unchanged). LAKE-HOLD (rule 7, `path_incision`):
    // also at every dry cell of the overflow path below the sill.
    let outlet = match lakes {
        Some((lb, cb)) => cb.outlet.map(|o| (o, crate::closed_basin::outlet_discharge_area_m2(routing, lb, o.reference_runoff_m_per_yr, o.path_incision))),
        None => None,
    };
    // LAKE-HOLD-2 (closed_basin.rs rule 8, `through_lakes`): each full
    // lake's sill cell; a cell whose receiver is such a lake is graded to the
    // lake's post-incision spill (the sill's stage-1 surface; receivers come
    // first, so the sill is already cut), clamped to [receiver bed, start
    // level]. None = law off: the start-of-step fill level, as before.
    let through_sills: Option<Vec<Option<usize>>> = match (lakes, &outlet) {
        (Some((lb, _)), Some((o, _))) if o.through_lakes => Some(crate::closed_basin::outlet_sill_cells(routing, lb)),
        _ => None,
    };
    let lake_base = |j: usize, z: &[f64]| -> f64 {
        let f = routing.fill_level_m[j];
        if let Some(sills) = &through_sills {
            let b = routing.water_body_id[j];
            if !ocean[j] && b != crate::water::DRY {
                if let Some(Some(sc)) = sills.get(b as usize) {
                    return elevation_m[j].max(f.min(z[*sc]));
                }
            }
        }
        f
    };

    // ---- 1. Stream power, implicit, receivers first, with the taper -------
    let mut detached = vec![0.0f64; n]; // m^3 leaving each land cell
    for &i in order {
        if wet[i] {
            continue;
        }
        let Some((j, len)) = recv[i] else {
            continue; // terminal floor: no fluvial lowering
        };
        let base = if wet[j] {
            lake_base(j, &st.z)
        } else {
            st.z[j]
        };
        if st.z[i] <= base {
            continue;
        }
        let f = params.k_per_my * routing.contributing_area_m2[i].powf(params.m) * dt_my / len;
        let z_new = (st.z[i] + f * base) / (1.0 + f);
        let mut amount = st.z[i] - z_new;
        // LAKE-OUTLET: the outlet channel's lowering of a full lake's sill
        // cell (None for every other cell and with the law off).
        let sill = match &outlet {
            Some((o, a)) if a[i] > 0.0 => Some(o.sill_lowering_m(st.z[i] - base, a[i], len, dt_my, params.k_per_my, params.m)),
            _ => None,
        };
        if amount <= 0.0 {
            if sill.is_none() {
                continue;
            }
            amount = 0.0;
        }
        let phi = rp.taper(relief[i].max(elevation_m[i] - base_level(sink[i])));
        if phi < 1.0 && amount > 0.0 {
            rr.taper_withheld_m3 += amount * (1.0 - phi) * areas[i];
            rr.tapered_cells += 1;
            amount *= phi;
            if amount <= 0.0 && sill.is_none() {
                continue;
            }
        }
        if let Some(s) = sill {
            rr.outlet_sills += 1;
            if st.z[i] - base <= 1.0 {
                rr.outlet_sills_level += 1;
            } else if s <= 0.0 {
                rr.outlet_sills_below_threshold += 1;
            }
            if s > amount {
                rr.outlet_sills_incised += 1;
                rr.outlet_extra_m3 += (s - amount) * areas[i];
                amount = s;
            }
            if amount <= 0.0 {
                continue;
            }
        }
        // SILL-REBOUND (`lake-outlet-sill-incision-flexural-support-v1`,
        // closed_basin.rs rule 6): a sill's lowering `drop` is its
        // post-isostasy drop; remove the compensated column for it.
        let support = match &outlet {
            Some((o, _)) if o.flexural_support && sill.is_some() => {
                let drop = amount;
                amount = col.thickness_to_raise(st.z[i] - drop, st.z[i], sea_level_m);
                Some((drop, amount, st.z[i]))
            }
            _ => None,
        };
        let removable = st.removable(i);
        if amount > removable {
            amount = removable;
            st.bound_hits += 1;
            if amount <= 0.0 {
                continue;
            }
        }
        st.lower(i, amount);
        detached[i] = amount * areas[i];
        if let Some((drop, full, z0)) = support {
            // The working surface drops by the lowering only (pro rata if
            // the available-material bound cut the removal).
            let d = if amount < full { drop * amount / full } else { drop };
            st.z[i] = z0 - d;
            rr.outlet_support_extra_m3 += (amount - d) * areas[i];
        }
    }

    let mut trace=trace_clock.map(|clock|super::fluvial_trace::FluvialFrame::new(clock,wet.clone(),ocean.to_vec(),
        recv.iter().map(|r|r.map(|(j,_)|j as u32)).collect(),order.iter().rev().map(|&i|i as u32).collect()));
    if stage_budgets {
        let t=trace.as_mut().ok_or_else(||ErosionError::Law("stage budgets require selected trace".into()))?;
        let mut budget=super::surface_budget::SurfaceBudget::new(grid,areas);budget.capture(&st);t.surface_budget=Some(budget);
    }
    // ---- 2a. Donor-first transport with graded-plain storage --------------
    let mut qin = vec![0.0f64; n];
    let mut pending = vec![0.0f64; n]; // volume arriving at standing water / floors
    for &i in order.iter().rev() {
        if wet[i] {
            continue;
        }
        let mut q = qin[i];
        let incoming_before=q;
        let (mut grade_evaluated,mut grade_need,mut grade_retained)=(false,0.,0.);
        let Some((j, len)) = recv[i] else {
            // Terminal floor: keeps what reaches it (legacy).
            q += detached[i];
            if let Some(t)=trace.as_mut(){t.rows.push(super::fluvial_trace::FluvialRow{cell:i as u32,receiver:None,
                continental:continental[i],grade_evaluated:false,incoming_m3:incoming_before,detached_m3:detached[i],
                grade_need_m3:0.,retained_m3:0.,outgoing_m3:q});}
            if q > 0.0 {
                pending[i] += q;
            }
            continue;
        };
        if q > 0.0 && continental[i] {
            let recv_level = if wet[j] {
                lake_base(j, &st.z)
            } else {
                st.z[j]
            };
            // LAKE-HOLD-2 (closed_basin.rs rule 9): the outlet path's grade.
            let s_a = match &outlet {
                Some((o, a)) if o.path_grade => o.path_grade_slope(rp.alluvial_grade_slope, a[i]),
                _ => rp.alluvial_grade_slope,
            };
            let g = recv_level + s_a * len;
            if st.z[i] < g {
                let need = col.thickness_to_raise(st.z[i], g, sea_level_m) * areas[i];
                let take = q.min(need);
                grade_evaluated=true;grade_need=need;grade_retained=take;
                if take > 0.0 {
                    deposit(&mut st, &mut rep, i, take / areas[i], Sink::LandFluvial);
                    rr.continental_grade_m3 += take;
                    q -= take;
                }
            }
        }
        q += detached[i];
        if let Some(t)=trace.as_mut(){t.rows.push(super::fluvial_trace::FluvialRow{cell:i as u32,receiver:Some(j as u32),
            continental:continental[i],grade_evaluated,incoming_m3:incoming_before,detached_m3:detached[i],
            grade_need_m3:grade_need,retained_m3:grade_retained,outgoing_m3:q});}
        if q <= 0.0 {
            continue;
        }
        if wet[j] {
            pending[j] += q;
            *transfers.entry((i, j, "fluvial_delivery")).or_default() += q;
        } else {
            qin[j] += q;
        }
    }

    if let Some(b)=trace.as_mut().and_then(|t|t.surface_budget.as_mut()){b.capture(&st);}
    // ---- 2b-d. Standing water: lakes, flooded continent, then legacy ------
    let ocean_cells: Vec<usize> = (0..n).filter(|&c| ocean[c]).collect();
    let ocean_area: f64 = ocean_cells.iter().map(|&c| areas[c]).sum();
    let mut surf = elevation_m.to_vec(); // compensated surface for capacity tests
    let mut legacy_pending = vec![0.0f64; n];
    let mut marine_pending = vec![0.0f64; n];
    place_lakes(
        grid,
        routing,
        elevation_m,
        ocean,
        sea_level_m,
        col,
        &mut st,
        &mut rep,
        &mut rr,
        &mut transfers,
        &mut surf,
        &pending,
        &mut marine_pending,
        lakes,
        rp.alluvial_grade_slope,
    )?;
    if let Some(b)=trace.as_mut().and_then(|t|t.surface_budget.as_mut()){b.capture(&st);}
    // Bypass before continental shelf capture. Keep within the existing
    // continental_shelf trace stage, whose schema has a fixed stage count.
    let fines_report = if let Some(fines) = marine_fines {
        super::marine_fines::place(grid, ocean, sea_level_m, *col, fines,
            &mut marine_pending, &mut surf, &mut st, &mut rep)?
    } else { super::marine_fines::MarineFinesReport::default() };
    rr.marine_fines_m3 = fines_report.deposited_m3;
    let flooded_continent: Vec<bool> = (0..n).map(|c| ocean[c] && continental[c]).collect();
    let mut shelf = ContinentalShelf::new(n);
    for s in 0..n {
        let v = marine_pending[s];
        if v <= 0.0 {
            continue;
        }
        rr.mouths += 1;
        let rem = shelf.place(grid, rp, &flooded_continent, s, v, &mut surf, sea_level_m, &mut st, &mut rep, &mut rr)?;
        legacy_pending[s] += rem;
    }
    if let Some(b)=trace.as_mut().and_then(|t|t.surface_budget.as_mut()){b.capture(&st);}
    let before_legacy = rep.total_stored_m3();
    let mut tally = LayerTally {
        law: rp.ocean_layer,
        ..Default::default()
    };
    place_legacy(
        grid,
        routing,
        elevation_m,
        ocean,
        sea_level_m,
        &ocean_cells,
        ocean_area,
        &mut st,
        &mut rep,
        &mut transfers,
        &legacy_pending,
        &mut tally,
    )?;
    rr.global_layer_mouths = tally.mouths;
    rr.global_layer_cell_raises = tally.raises;
    rr.global_layer_cells_above_sea = tally.above_sea;
    rr.legacy_marine_m3 = rep.total_stored_m3() - before_legacy;
    rr.global_layer_m3 = rep.global_remainder_m3;

    if let Some(b)=trace.as_mut().and_then(|t|t.surface_budget.as_mut()){b.capture(&st);}
    // ---- 3. Hillslope diffusion: the legacy scheme, verbatim --------------
    let d = params.diffusivity_m2_per_my;
    let nu_max = stencil.diffusion_number_max(d, dt_my, &wet);
    let substeps = if d == 0.0 || nu_max == 0.0 {
        1
    } else {
        ((nu_max / params.stability_limit).ceil() as u32).max(1)
    };
    let sub_dt = dt_my / substeps as f64;
    if d > 0.0 {
        let mut flux_dz = vec![0.0f64; n];
        let mut moves: Vec<(usize, usize, f64, bool)> = Vec::new();
        let mut outgoing = vec![0.0f64; n];
        for _ in 0..substeps {
            for v in flux_dz.iter_mut() {
                *v = 0.0;
            }
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
                    rr.hillslope_to_water_m3 += flux;
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

    if let Some(b)=trace.as_mut().and_then(|t|t.surface_budget.as_mut()){b.capture(&st);}
    for c in 0..n {
        let v = st.sediment_added[c] * areas[c];
        if continental[c] {
            rr.added_on_continental_m3 += v;
        } else {
            rr.added_on_oceanic_m3 += v;
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
    let step = ErosionStep {
        transfers,
        spread_bound_hits: rep.lake_spread_bound_hits
            + rep.marine_search_bound_hits
            + rr.shelf_search_bound_hits,
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
        deposition: None,
        ocean_layer_mouths: rr.global_layer_mouths,
        ocean_layer_cell_raises: rr.global_layer_cell_raises,
        ocean_layer_cells_above_sea: rr.global_layer_cells_above_sea,
    };
    if let Some(b)=trace.as_mut().and_then(|t|t.surface_budget.as_mut()){b.finish(&step);}
    if let Some(t)=&trace {if let Some(b)=&t.surface_budget{b.validate(t,[grid.nlat(),grid.nlon()]).map_err(ErosionError::Law)?;}}
    THREAD_TOTALS.with(|t| t.borrow_mut().accumulate(&rr));
    if let Some(t)=&trace {t.validate(n,dt_my).map_err(ErosionError::Law)?;}
    Ok((step, rr, trace, fines_report))
}

/// Local relief (m) of every cell: the max minus the min of the level over
/// the cell and its 8 neighbours, the level being the water level for ocean
/// (sea level) and lake cells (spill level) and the surface elevation
/// otherwise. Evaluated on the start-of-step fields.
pub fn local_relief(
    grid: &GridSpec,
    routing: &Routing,
    elevation_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
) -> Result<Vec<f64>, ErosionError> {
    let n = grid.cells();
    let level: Vec<f64> = (0..n)
        .map(|c| {
            if ocean[c] {
                sea_level_m
            } else {
                elevation_m[c].max(routing.fill_level_m[c])
            }
        })
        .collect();
    let mut out = vec![0.0f64; n];
    for c in 0..n {
        let (ci, cj) = grid
            .ij(c as u64)
            .map_err(|e| ErosionError::Grid(e.to_string()))?;
        let (mut lo, mut hi) = (level[c], level[c]);
        for (ti, tj, _) in grid
            .neighbors(ci, cj)
            .map_err(|e| ErosionError::Grid(e.to_string()))?
        {
            let v = level[ti * grid.nlon() + tj];
            lo = lo.min(v);
            hi = hi.max(v);
        }
        out[c] = hi - lo;
    }
    Ok(out)
}

/// (b) Lakes: compensated fill to the spill level over the lake body
/// (bounded breadth-first spread, as in the legacy step); a full lake passes
/// the rest downstream; arrivals at the ocean are collected per mouth cell.
#[allow(clippy::too_many_arguments)]
fn place_lakes(
    grid: &GridSpec,
    routing: &Routing,
    elevation_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    col: &CompensatedColumn,
    st: &mut Lowering<'_>,
    rep: &mut DepositionReport,
    rr: &mut RetentionReport,
    transfers: &mut std::collections::BTreeMap<(usize, usize, &str), f64>,
    surf: &mut [f64],
    pending: &[f64],
    marine_pending: &mut [f64],
    lakes: Option<(&crate::closed_basin::LakeBalance, &crate::closed_basin::ClosedBasinParams)>,
    alluvial_grade_slope: f64,
) -> Result<(), ErosionError> {
    let n = pending.len();
    let areas = st.areas;
    let is_lake = |c: usize| !ocean[c] && routing.fill_level_m[c] > elevation_m[c] + 1.0e-9;
    // closed-basin-water-balance-v1: the floor-fill pool of each lake body.
    let mut pool: Vec<f64> = match lakes {
        Some(_) => vec![0.0; routing.water_bodies.len()],
        None => Vec::new(),
    };
    for t0 in 0..n {
        let mut remaining = pending[t0];
        if remaining <= 0.0 {
            continue;
        }
        if let Some((_, cp)) = lakes {
            if is_lake(t0) {
                // Suspended share: settles on the floor (focused fill below).
                let s = remaining * cp.suspended_load_fraction;
                pool[routing.water_body_id[t0] as usize] += s;
                remaining -= s;
                if remaining <= 0.0 {
                    continue;
                }
            }
        }
        let mut start = t0;
        let mut hops = 0usize;
        while remaining > 0.0 {
            if ocean[start] {
                marine_pending[start] += remaining;
                break;
            }
            if !is_lake(start) {
                // A terminal floor with no outlet keeps what reaches it.
                deposit(st, rep, start, remaining / areas[start], Sink::LandFloor);
                rr.land_floor_m3 += remaining;
                break;
            }
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
                        rr.lake_m3 += remaining;
                        remaining = 0.0;
                    } else {
                        surf[c] = level;
                        deposit(st, rep, c, need / areas[c], Sink::Lake);
                        rr.lake_m3 += need;
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
            if lakes.is_some() {
                // The delta is full near the inlet: the rest joins the floor
                // fill of the same lake instead of leaving it.
                pool[routing.water_body_id[start] as usize] += remaining;
                break;
            }
            // A full lake: carry the rest downstream to the next water.
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
                    deposit(st, rep, cur, remaining / areas[cur], Sink::LandFloor);
                    rr.land_floor_m3 += remaining;
                    remaining = 0.0;
                }
            }
        }
    }
    if let Some((lb, cp)) = lakes {
        // LAKE-EXCESS: `filled-basin-alluvial-grade-v1` when d_c is set.
        let grade = cp.filled_basin_grading_depth_m.map(|d| (d, alluvial_grade_slope));
        focused_lake_fill(grid, routing, elevation_m, ocean, sea_level_m, col, st, rep, rr, transfers, surf, marine_pending, lb, &mut pool, grade)?;
    }
    Ok(())
}

/// closed-basin-water-balance-v1, part 2: fill each lake body's lowest cells
/// to one common compensated level with its pool, upstream lakes first; a
/// body that reaches its spill level everywhere passes the rest to the next
/// water down its outlet (ocean: the mouth's marine pending; lake: that
/// lake's pool, processed later in `lb.order`).
#[allow(clippy::too_many_arguments)]
fn focused_lake_fill(
    grid: &GridSpec,
    routing: &Routing,
    elevation_m: &[f64],
    ocean: &[bool],
    sea_level_m: f64,
    col: &CompensatedColumn,
    st: &mut Lowering<'_>,
    rep: &mut DepositionReport,
    rr: &mut RetentionReport,
    transfers: &mut std::collections::BTreeMap<(usize, usize, &str), f64>,
    surf: &mut [f64],
    marine_pending: &mut [f64],
    lb: &crate::closed_basin::LakeBalance,
    pool: &mut [f64],
    grade: Option<(f64, f64)>,
) -> Result<(), ErosionError> {
    let n = surf.len();
    let areas = st.areas;
    let mut cells: Vec<Vec<usize>> = vec![Vec::new(); routing.water_bodies.len()];
    for c in 0..n {
        let b = routing.water_body_id[c];
        if b != crate::water::DRY && !ocean[c] && routing.fill_level_m[c] > elevation_m[c] + 1.0e-9 {
            cells[b as usize].push(c);
        }
    }
    // LAKE-EXCESS (`filled-basin-alluvial-grade-v1`): flow distance (m) of
    // every lake cell to its body's outlet along the routing's flow tree
    // (receivers first; a lake cell's flow leads to the outlet on the
    // filled surface). Only computed when the law is on.
    let dist_to_outlet: Vec<f64> = match grade {
        None => Vec::new(),
        Some(_) => {
            let mut d = vec![0.0f64; n];
            for &i in &routing.receiver_first {
                if !is_lake_cell(routing, elevation_m, ocean, i) {
                    continue;
                }
                if let Some(ElementId(j)) = routing.flow_to[i] {
                    let j = j as usize;
                    let (ii, ij) = grid.ij(i as u64).map_err(|e| ErosionError::Grid(e.to_string()))?;
                    let (ji, jj) = grid.ij(j as u64).map_err(|e| ErosionError::Grid(e.to_string()))?;
                    let len = grid.center_distance_m(ii, ij, ji, jj).map_err(|e| ErosionError::Grid(e.to_string()))?.max(1.0);
                    let same = is_lake_cell(routing, elevation_m, ocean, j) && routing.water_body_id[j] == routing.water_body_id[i];
                    d[i] = len + if same { d[j] } else { 0.0 };
                }
            }
            d
        }
    };
    for &b in &lb.order {
        let v = pool[b];
        if v <= 0.0 || cells[b].is_empty() {
            continue;
        }
        let spill = routing.water_bodies[b].water_surface_m;
        let need_at = |z: f64, surf: &[f64]| -> f64 {
            cells[b]
                .iter()
                .map(|&c| col.thickness_to_raise(surf[c], z.min(routing.fill_level_m[c]), sea_level_m) * areas[c])
                .sum()
        };
        let cap = need_at(spill, surf);
        let (level, rest) = if cap <= v {
            (spill, v - cap)
        } else {
            let mut lo = cells[b].iter().map(|&c| surf[c]).fold(f64::INFINITY, f64::min);
            let mut hi = spill;
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if need_at(mid, surf) <= v {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            (lo, 0.0)
        };
        // LAKE-EXCESS: a filled basin (the level fill reaches within d_c of
        // the spill level) takes the graded surface z_g - delta instead.
        let graded: Option<(f64, f64)> = match grade {
            Some((d_c, slope)) if level >= spill - d_c => Some((spill, slope)),
            _ => None,
        };
        let target_of = |c: usize, delta: f64| -> f64 {
            match graded {
                Some((sp, slope)) => sp + slope * dist_to_outlet[c] - delta,
                None => level.min(routing.fill_level_m[c]),
            }
        };
        let (delta, rest) = match graded {
            None => (0.0, rest),
            Some(_) => {
                let need_g = |delta: f64, surf: &[f64]| -> f64 {
                    cells[b].iter().map(|&c| col.thickness_to_raise(surf[c], target_of(c, delta), sea_level_m) * areas[c]).sum()
                };
                let full = need_g(0.0, surf);
                if full <= v {
                    (0.0, v - full)
                } else {
                    let mut lo = 0.0f64; // need > v
                    let mut hi = cells[b].iter().map(|&c| target_of(c, 0.0) - surf[c]).fold(0.0f64, f64::max); // need = 0
                    for _ in 0..80 {
                        let mid = 0.5 * (lo + hi);
                        if need_g(mid, surf) <= v {
                            hi = mid;
                        } else {
                            lo = mid;
                        }
                    }
                    (hi, 0.0)
                }
            }
        };
        if graded.is_some() {
            rr.graded_basins += 1;
        }
        // Raise every cell below `level` (graded: below z_g - delta); book
        // exactly v - rest (the tiny bisection shortfall goes to the lowest
        // raised cell).
        let mut placed = 0.0;
        let mut lowest: Option<usize> = None;
        for &c in &cells[b] {
            let target = target_of(c, delta);
            if surf[c] < target {
                let h = col.thickness_to_raise(surf[c], target, sea_level_m);
                if h > 0.0 {
                    if lowest.map_or(true, |l| surf[c] < surf[l]) {
                        lowest = Some(c);
                    }
                    surf[c] = col.surface_after(surf[c], h, sea_level_m);
                    deposit(st, rep, c, h, Sink::Lake);
                    placed += h * areas[c];
                }
            }
        }
        let short = (v - rest) - placed;
        if short > 0.0 {
            let c = lowest.unwrap_or(cells[b][0]);
            let h = short / areas[c];
            surf[c] = col.surface_after(surf[c], h, sea_level_m);
            deposit(st, rep, c, h, Sink::Lake);
            placed += short;
        }
        rr.lake_m3 += placed;
        if graded.is_some() {
            rr.graded_basin_m3 += placed;
        }
        if rest > 0.0 {
            // Overflow along the outlet to the next water body.
            let exit = routing.water_bodies[b].outlet_edges.first().map(|e| (e[0].0 as usize, e[1].0 as usize));
            let mut target: Option<usize> = None;
            let mut from = cells[b][0];
            if let Some((a, j)) = exit {
                from = a;
                let mut cur = j;
                let mut hops = 0usize;
                loop {
                    if ocean[cur] || (is_lake_cell(routing, elevation_m, ocean, cur) && routing.water_body_id[cur] as usize != b) {
                        target = Some(cur);
                        break;
                    }
                    match routing.flow_to[cur] {
                        Some(ElementId(k)) => {
                            from = cur;
                            cur = k as usize;
                        }
                        None => break,
                    }
                    hops += 1;
                    if hops > n {
                        return Err(ErosionError::Grid("cycle during focused lake overflow".into()));
                    }
                }
                if target.is_none() {
                    // A terminal floor down the outlet keeps the rest.
                    deposit(st, rep, cur, rest / areas[cur], Sink::LandFloor);
                    rr.land_floor_m3 += rest;
                    continue;
                }
            }
            match target {
                Some(t) if ocean[t] => {
                    marine_pending[t] += rest;
                    *transfers.entry((from, t, "lake_overflow")).or_default() += rest;
                }
                Some(t) => {
                    let d = routing.water_body_id[t] as usize;
                    pool[d] += rest;
                    *transfers.entry((from, t, "lake_overflow")).or_default() += rest;
                }
                None => {
                    // No outlet edge (should not happen for a routed lake):
                    // keep it on the lake's lowest cell.
                    let c = cells[b][0];
                    deposit(st, rep, c, rest / areas[c], Sink::LandFloor);
                    rr.land_floor_m3 += rest;
                }
            }
        }
    }
    Ok(())
}

fn is_lake_cell(routing: &Routing, elevation_m: &[f64], ocean: &[bool], c: usize) -> bool {
    !ocean[c] && routing.fill_level_m[c] > elevation_m[c] + 1.0e-9
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Key(f64, usize);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    /// Min-heap on distance, then on element id.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .0
            .total_cmp(&self.0)
            .then_with(|| other.1.cmp(&self.1))
    }
}

/// (c) Flooded-continent shelf placement (generation-stamped search state).
struct ContinentalShelf {
    settled: Vec<u32>,
    generation: u32,
    heap: std::collections::BinaryHeap<Key>,
}

impl ContinentalShelf {
    fn new(n: usize) -> Self {
        Self {
            settled: vec![0; n],
            generation: 0,
            heap: std::collections::BinaryHeap::new(),
        }
    }

    /// Place `volume_m3` arriving at ocean cell `mouth` on flooded
    /// continental cells up to wave base; returns what could not be placed.
    #[allow(clippy::too_many_arguments)]
    fn place(
        &mut self,
        grid: &GridSpec,
        rp: &RetentionParams,
        flooded_continent: &[bool],
        mouth: usize,
        volume_m3: f64,
        surf: &mut [f64],
        sea_level_m: f64,
        st: &mut Lowering<'_>,
        rep: &mut DepositionReport,
        rr: &mut RetentionReport,
    ) -> Result<f64, ErosionError> {
        let areas = st.areas;
        let neighbours = |c: usize| -> Result<Vec<(usize, f64)>, ErosionError> {
            let (ci, cj) = grid
                .ij(c as u64)
                .map_err(|e| ErosionError::Grid(e.to_string()))?;
            let mut list: Vec<(usize, f64)> = grid
                .neighbors(ci, cj)
                .map_err(|e| ErosionError::Grid(e.to_string()))?
                .into_iter()
                .map(|(ti, tj, nb)| (ti * grid.nlon() + tj, nb.distance_m.max(1.0)))
                .filter(|&(nb, _)| flooded_continent[nb])
                .collect();
            list.sort_unstable_by_key(|&(nb, _)| nb);
            Ok(list)
        };
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.settled.iter_mut().for_each(|s| *s = 0);
            self.generation = 1;
        }
        let gen = self.generation;
        self.heap.clear();
        if flooded_continent[mouth] {
            self.heap.push(Key(0.0, mouth));
        } else {
            for (nb, d) in neighbours(mouth)? {
                self.heap.push(Key(d, nb));
            }
        }
        if self.heap.is_empty() {
            rr.mouths_without_continental_shelf += 1;
            return Ok(volume_m3);
        }
        let mut remaining = volume_m3;
        let target = sea_level_m - rp.shelf_wave_base_depth_m;
        let col = &rp.column;
        let mut settled_count = 0usize;
        while let Some(Key(d, c)) = self.heap.pop() {
            if self.settled[c] == gen {
                continue;
            }
            if settled_count >= rp.search_max_cells {
                rr.shelf_search_bound_hits += 1;
                break;
            }
            self.settled[c] = gen;
            settled_count += 1;
            if surf[c] < target {
                let need = col.thickness_to_raise(surf[c], target, sea_level_m) * areas[c];
                let class = MarineClass::of_depth(sea_level_m - surf[c]);
                if remaining < need {
                    let h = remaining / areas[c];
                    surf[c] = col.surface_after(surf[c], h, sea_level_m);
                    deposit(st, rep, c, h, Sink::Marine(class));
                    rr.continental_shelf_m3 += remaining;
                    remaining = 0.0;
                    break;
                }
                surf[c] = target;
                deposit(st, rep, c, need / areas[c], Sink::Marine(class));
                rr.continental_shelf_m3 += need;
                remaining -= need;
                if remaining <= 0.0 {
                    remaining = 0.0;
                    break;
                }
            }
            // At or above wave base: sediment crosses this cell.
            for (nb, dist) in neighbours(c)? {
                if self.settled[nb] != gen {
                    self.heap.push(Key(d + dist, nb));
                }
            }
        }
        Ok(remaining)
    }
}
