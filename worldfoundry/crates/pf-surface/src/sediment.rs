//! Sediment mass ledger for the PREVIEW tier (ADR 032): the spec §54.4
//! gate, eroded = deposited + exported +- declared error, per step AND
//! cumulative, built on the contract type `SedimentMassBalanceV0` (which
//! re-validates the same identity on load).
//!
//! The declared error is a NUMERICAL budget, not a tolerance to hide
//! behind: relative floating-point rounding of the summations (declared
//! per epoch as `rel_tol * eroded + abs_floor`). An epoch whose residual
//! exceeds it is refused at record time, loudly.
//!
//! WHAT THIS TIER DOES NOT PRODUCE: the S5a provenance-preserving
//! `SedimentTransferRecord` stream (parcel lineages, grain classes, phase
//! split). Those require the Lagrangian parcel ledger, which the preview
//! backend does not drive. The mass gate is exercised; the provenance
//! seam is not, and the artifact says so through its assumptions list.

use pf_contracts::domain::{SedimentEpochBalanceV0, SedimentMassBalanceV0};

#[derive(Debug, Clone, PartialEq)]
pub enum SedimentError {
    NonFinite {
        what: &'static str,
        value: f64,
    },
    EpochDoesNotClose {
        t_my: f64,
        eroded_kg: f64,
        deposited_kg: f64,
        exported_kg: f64,
        residual_kg: f64,
        declared_error_kg: f64,
    },
    NonMonotonicTime {
        last: f64,
        got: f64,
    },
    BadParam {
        what: &'static str,
        value: f64,
    },
    Contract(String),
    ControlVolumeDoesNotClose {
        reservoir: String,
        stage: String,
        t_my: f64,
        residual_m3: f64,
        declared_error_m3: f64,
    },
}

impl std::fmt::Display for SedimentError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SedimentError::NonFinite { what, value } => write!(f, "nonfinite {what}: {value:?}"),
            SedimentError::EpochDoesNotClose {
                t_my,
                eroded_kg,
                deposited_kg,
                exported_kg,
                residual_kg,
                declared_error_kg,
            } => write!(
                f,
                "sediment epoch at {t_my:?} My does not close: eroded {eroded_kg:?} != deposited {deposited_kg:?} + exported {exported_kg:?} (residual {residual_kg:?} > declared {declared_error_kg:?})"
            ),
            SedimentError::NonMonotonicTime { last, got } => {
                write!(f, "epoch at {got:?} My not after previous {last:?} My")
            }
            SedimentError::BadParam { what, value } => write!(f, "bad parameter {what} = {value:?}"),
            SedimentError::Contract(why) => write!(f, "contract validation: {why}"),
            SedimentError::ControlVolumeDoesNotClose {
                reservoir,
                stage,
                t_my,
                residual_m3,
                declared_error_m3,
            } => write!(
                f,
                "{reservoir} control volume ({stage}) at {t_my:?} My does not close: residual {residual_m3:?} m3 > declared {declared_error_m3:?} m3"
            ),
        }
    }
}

impl std::error::Error for SedimentError {}

/// Running ledger. Mass = volume * rock density (declared by the caller).
#[derive(Debug, Clone, PartialEq)]
pub struct SedimentLedger {
    rock_density_kg_m3: f64,
    rel_tol: f64,
    abs_floor_kg: f64,
    epochs: Vec<SedimentEpochBalanceV0>,
}

impl SedimentLedger {
    /// `rel_tol` is the declared relative rounding budget per epoch;
    /// `abs_floor_kg` covers epochs with negligible erosion.
    pub fn new(rock_density_kg_m3: f64, rel_tol: f64, abs_floor_kg: f64) -> Result<Self, SedimentError> {
        for (what, v) in [
            ("rock_density_kg_m3", rock_density_kg_m3),
            ("rel_tol", rel_tol),
            ("abs_floor_kg", abs_floor_kg),
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(SedimentError::BadParam { what, value: v });
            }
        }
        Ok(Self {
            rock_density_kg_m3,
            rel_tol,
            abs_floor_kg,
            epochs: Vec::new(),
        })
    }

    pub fn rock_density_kg_m3(&self) -> f64 {
        self.rock_density_kg_m3
    }

    pub fn epochs(&self) -> &[SedimentEpochBalanceV0] {
        &self.epochs
    }

    /// Record one epoch from volumes (m^3). Refuses loudly if the epoch
    /// does not close within its declared error.
    pub fn record_epoch_m3(
        &mut self,
        t_my: f64,
        eroded_m3: f64,
        deposited_m3: f64,
        exported_m3: f64,
    ) -> Result<(), SedimentError> {
        for (what, v) in [
            ("t_my", t_my),
            ("eroded_m3", eroded_m3),
            ("deposited_m3", deposited_m3),
            ("exported_m3", exported_m3),
        ] {
            if !v.is_finite() || v < 0.0 {
                return Err(SedimentError::NonFinite { what, value: v });
            }
        }
        if let Some(last) = self.epochs.last() {
            if t_my <= last.t_my {
                return Err(SedimentError::NonMonotonicTime {
                    last: last.t_my,
                    got: t_my,
                });
            }
        }
        let rho = self.rock_density_kg_m3;
        let eroded_kg = eroded_m3 * rho;
        let deposited_kg = deposited_m3 * rho;
        let exported_kg = exported_m3 * rho;
        let declared_error_kg = self.rel_tol * eroded_kg + self.abs_floor_kg;
        let residual_kg = eroded_kg - deposited_kg - exported_kg;
        if residual_kg.abs() > declared_error_kg {
            return Err(SedimentError::EpochDoesNotClose {
                t_my,
                eroded_kg,
                deposited_kg,
                exported_kg,
                residual_kg,
                declared_error_kg,
            });
        }
        self.epochs.push(SedimentEpochBalanceV0 {
            t_my,
            eroded_kg,
            deposited_kg,
            exported_kg,
            declared_error_kg,
        });
        Ok(())
    }

    /// The durable balance with cumulative sums and the summed declared
    /// error (validated by the contract type before return).
    pub fn balance(&self) -> Result<SedimentMassBalanceV0, SedimentError> {
        let (mut er, mut dep, mut exp, mut err) = (0.0, 0.0, 0.0, 0.0);
        for e in &self.epochs {
            er += e.eroded_kg;
            dep += e.deposited_kg;
            exp += e.exported_kg;
            err += e.declared_error_kg;
        }
        let b = SedimentMassBalanceV0 {
            epochs: self.epochs.clone(),
            cumulative_eroded_kg: er,
            cumulative_deposited_kg: dep,
            cumulative_exported_kg: exp,
            cumulative_declared_error_kg: err,
        };
        b.validate()
            .map_err(|e| SedimentError::Contract(e.to_string()))?;
        Ok(b)
    }
}

/// One reservoir's control-volume epoch (lane E, SEED-EROSION-BUG,
/// 2026-10-11; review F16): initial + inflows - outflows + remap_drift =
/// final + residual, with every term NAMED. The remap drift (the measured,
/// unexplained change of a reservoir under the semi-Lagrangian carry) is a
/// separate signed term, never folded into the residual, and the residual
/// must sit inside the declared numerical budget. Placement plugs that move
/// material without creating or destroying it (the global marine layer)
/// are reported by name in `diagnostics` and do not enter the balance.
#[derive(Debug, Clone, PartialEq)]
pub struct ControlVolumeEpoch {
    pub reservoir: String,
    pub stage: String,
    pub t_my: f64,
    pub initial_m3: f64,
    pub inflows_m3: Vec<(String, f64)>,
    pub outflows_m3: Vec<(String, f64)>,
    pub remap_drift_m3: f64,
    pub final_m3: f64,
    pub residual_m3: f64,
    pub declared_error_m3: f64,
    pub diagnostics_m3: Vec<(String, f64)>,
}

/// Running control-volume ledger. The declared error per epoch is
/// `rel_tol * (|initial| + sum|flows| + |drift| + |final|) + abs_floor_m3`:
/// floating-point rounding of area-weighted sums, not a tolerance for
/// missing transfers. An epoch outside it is refused, loudly, and not
/// recorded.
#[derive(Debug, Clone, PartialEq)]
pub struct ControlVolumeLedger {
    rel_tol: f64,
    abs_floor_m3: f64,
    epochs: Vec<ControlVolumeEpoch>,
}

/// Inputs of one control-volume epoch (borrowed; see `ControlVolumeEpoch`).
#[derive(Debug, Clone, Copy)]
pub struct ControlVolumeTerms<'a> {
    pub reservoir: &'a str,
    pub stage: &'a str,
    pub t_my: f64,
    pub initial_m3: f64,
    pub inflows_m3: &'a [(&'a str, f64)],
    pub outflows_m3: &'a [(&'a str, f64)],
    pub remap_drift_m3: f64,
    pub final_m3: f64,
    pub diagnostics_m3: &'a [(&'a str, f64)],
}

impl ControlVolumeLedger {
    pub fn new(rel_tol: f64, abs_floor_m3: f64) -> Result<Self, SedimentError> {
        for (what, v) in [("rel_tol", rel_tol), ("abs_floor_m3", abs_floor_m3)] {
            if !v.is_finite() || v <= 0.0 {
                return Err(SedimentError::BadParam { what, value: v });
            }
        }
        Ok(Self {
            rel_tol,
            abs_floor_m3,
            epochs: Vec::new(),
        })
    }

    pub fn epochs(&self) -> &[ControlVolumeEpoch] {
        &self.epochs
    }

    /// Record one epoch. Flows must be finite and >= 0 (direction is
    /// carried by the list they are in); the drift is signed.
    pub fn record(&mut self, t: &ControlVolumeTerms<'_>) -> Result<&ControlVolumeEpoch, SedimentError> {
        for (what, v) in [
            ("t_my", t.t_my),
            ("initial_m3", t.initial_m3),
            ("final_m3", t.final_m3),
            ("remap_drift_m3", t.remap_drift_m3),
        ] {
            if !v.is_finite() {
                return Err(SedimentError::NonFinite { what, value: v });
            }
        }
        for (_, v) in t.inflows_m3.iter().chain(t.outflows_m3) {
            if !v.is_finite() || *v < 0.0 {
                return Err(SedimentError::NonFinite {
                    what: "control-volume flow",
                    value: *v,
                });
            }
        }
        let inflow: f64 = t.inflows_m3.iter().map(|(_, v)| v).sum();
        let outflow: f64 = t.outflows_m3.iter().map(|(_, v)| v).sum();
        let residual_m3 = t.initial_m3 + inflow - outflow + t.remap_drift_m3 - t.final_m3;
        let scale = t.initial_m3.abs() + inflow + outflow + t.remap_drift_m3.abs() + t.final_m3.abs();
        let declared_error_m3 = self.rel_tol * scale + self.abs_floor_m3;
        if residual_m3.abs() > declared_error_m3 {
            return Err(SedimentError::ControlVolumeDoesNotClose {
                reservoir: t.reservoir.to_string(),
                stage: t.stage.to_string(),
                t_my: t.t_my,
                residual_m3,
                declared_error_m3,
            });
        }
        let own = |v: &[(&str, f64)]| v.iter().map(|(k, x)| (k.to_string(), *x)).collect();
        self.epochs.push(ControlVolumeEpoch {
            reservoir: t.reservoir.to_string(),
            stage: t.stage.to_string(),
            t_my: t.t_my,
            initial_m3: t.initial_m3,
            inflows_m3: own(t.inflows_m3),
            outflows_m3: own(t.outflows_m3),
            remap_drift_m3: t.remap_drift_m3,
            final_m3: t.final_m3,
            residual_m3,
            declared_error_m3,
            diagnostics_m3: own(t.diagnostics_m3),
        });
        Ok(self.epochs.last().expect("just pushed"))
    }

    /// Cumulative named totals per (reservoir, stage, term) over all
    /// epochs, plus the cumulative drift and residual per reservoir/stage.
    pub fn totals(&self) -> std::collections::BTreeMap<String, f64> {
        let mut out = std::collections::BTreeMap::new();
        for e in &self.epochs {
            let key = |term: &str| format!("{}/{}/{}", e.reservoir, e.stage, term);
            for (k, v) in &e.inflows_m3 {
                *out.entry(key(&format!("in:{k}"))).or_insert(0.0) += v;
            }
            for (k, v) in &e.outflows_m3 {
                *out.entry(key(&format!("out:{k}"))).or_insert(0.0) += v;
            }
            for (k, v) in &e.diagnostics_m3 {
                *out.entry(key(&format!("diag:{k}"))).or_insert(0.0) += v;
            }
            *out.entry(key("remap_drift")).or_insert(0.0) += e.remap_drift_m3;
            *out.entry(key("residual")).or_insert(0.0) += e.residual_m3;
        }
        out
    }
}

// ===========================================================================
// Lane SD (SEED-EXTRAS D7, 2026-10-12): compensated accommodation and
// marine shelf placement, law id `shelf-wave-base-v1`.
//
// WHAT IT REPLACES (only when the law is selected; the legacy path in
// `erosion::erode_step_bounded` is untouched): the legacy placement fills
// receiving water up to its level with the UNCOMPENSATED thickness (water
// depth x area), stops after `MAX_SPREAD_CELLS` cells and spreads the rest
// as one uniform layer over the whole ocean (lane E: 73 % of marine
// delivery on seed 240907). It also fills the ocean to SEA LEVEL, so margins
// end as sediment plains at about +40 m (lane K2), not shelves.
//
// PHYSICS:
// 1. Compensated accommodation (local Airy, the same mass balance the
//    isostasy laws use; sediment loads at the law's sediment density). A
//    sediment layer h added to a column raises its surface by
//    h (rho_m - rho_s) / (rho_m - rho_w) while the surface stays under
//    water and h (rho_m - rho_s) / rho_m above it. Both preview isostasy
//    laws (single-column reference Airy with water loading below the datum,
//    and the two-column common-depth law with its water-load factor
//    rho_m / (rho_m - rho_w)) have exactly these derivatives. So the
//    thickness that brings a cell's surface to a level is the depth times
//    the inverse factor (about 4.5 x under water for 3300 / 2800 / 1030
//    kg/m3), not the depth itself. Lane E's hypothesis was that ignoring
//    this underestimates accommodation by about 4.5 x.
//    Declared approximation: within a step, sea level is held fixed. Under
//    the volume-primary sea-level law, the sea-level change caused by
//    displaced water appears at the next isostatic solve.
// 2. Shelf top at wave base. Wave and current reworking keeps sediment
//    from aggrading above the shelf-top level, sea - SHELF_WAVE_BASE_DEPTH_M.
//    A cell at or above that level passes sediment across itself (the
//    topset is a transport surface). A deeper cell receives sediment up to
//    that level (foreset deposition, the shelf edge advancing).
// 3. Progradation from the river mouth. From each entry cell (the first
//    ocean cell on a river's path, or the ocean cell a hillslope flux
//    enters), a Dijkstra search by great-circle path length over ocean
//    cells settles the nearest cells first, ties broken by the lower
//    element id. It expands only from cells that stand at or above the
//    shelf top after filling. Sediment can therefore cross only built shelf
//    to reach the next accommodation, so the shelf grows as a contiguous
//    body from the mouth, and the partly filled cell is the advancing edge.
// 4. Beyond capacity. If the search reaches the declared computational
//    bound, or runs out of ocean, before placing everything, the rest is
//    the GLOBAL REMAINDER. It is spread as a uniform layer over the ocean
//    reservoir, declared, named and measured; the target is small. Nothing
//    is created or destroyed.
// Reporting classes, which are bookkeeping and not physics, are set by the
// water depth of the receiving cell before deposition: shelf <= 200 m (the
// R3 window), slope <= 2500 m, deep beyond.
// NOT MODELLED: subaerial delta topsets (sub-cell at 70 km), grain-size
// partition and fines bypass to the deep sea, turbidity currents, contour
// currents, pelagic rain, compaction, flexure (local Airy only), and
// post-rift thermal subsidence of continental margins.
// ===========================================================================

/// Law id of the marine placement implemented here.
pub const SHELF_WAVE_BASE_LAW_ID: &str = "shelf-wave-base-v1";
/// Depth below sea level to which the shelf top aggrades (m). Declared
/// BEFORE any run (lane SD FALSIFIERS.md) and never tuned: storm wave base
/// and cell-mean shelf-top depth. Earth's shelf break averages about 130 m
/// (recalled, UNVERIFIED), and the mean depth of a graded shelf is shallower.
/// Project assumption awaiting review.
pub const SHELF_WAVE_BASE_DEPTH_M: f64 = 100.0;
/// Computational bound: ocean cells settled per entry search. This is not a
/// physical capacity. Hits are recorded, and what is left becomes the
/// measured global remainder.
pub const MARINE_SEARCH_MAX_CELLS: usize = 16_384;
/// Reporting bins (water depth before deposition, m): classification only.
pub const SHELF_CLASS_MAX_DEPTH_M: f64 = 200.0;
pub const SLOPE_CLASS_MAX_DEPTH_M: f64 = 2500.0;

/// Densities of the local-Airy column response to a sediment load. The
/// caller passes the isostasy law's own values (mantle, the law's sediment
/// density, sea water).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CompensatedColumn {
    pub mantle_density_kg_m3: f64,
    pub sediment_density_kg_m3: f64,
    pub water_density_kg_m3: f64,
}

impl CompensatedColumn {
    pub fn validate(&self) -> Result<(), SedimentError> {
        for (what, v) in [
            ("mantle_density_kg_m3", self.mantle_density_kg_m3),
            ("sediment_density_kg_m3", self.sediment_density_kg_m3),
            ("water_density_kg_m3", self.water_density_kg_m3),
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(SedimentError::BadParam { what, value: v });
            }
        }
        if !(self.water_density_kg_m3 < self.sediment_density_kg_m3
            && self.sediment_density_kg_m3 < self.mantle_density_kg_m3)
        {
            return Err(SedimentError::BadParam {
                what: "densities must satisfy water < sediment < mantle",
                value: self.sediment_density_kg_m3,
            });
        }
        Ok(())
    }

    /// Surface rise per metre of sediment while the surface is under water.
    pub fn rise_per_m_submerged(&self) -> f64 {
        (self.mantle_density_kg_m3 - self.sediment_density_kg_m3)
            / (self.mantle_density_kg_m3 - self.water_density_kg_m3)
    }

    /// Surface rise per metre of sediment in air.
    pub fn rise_per_m_subaerial(&self) -> f64 {
        (self.mantle_density_kg_m3 - self.sediment_density_kg_m3) / self.mantle_density_kg_m3
    }

    /// Sediment thickness (m) that brings the compensated surface from `z0`
    /// to `z1` (0 when `z1 <= z0`), piecewise across sea level.
    pub fn thickness_to_raise(&self, z0: f64, z1: f64, sea: f64) -> f64 {
        if z1 <= z0 {
            return 0.0;
        }
        let mut h = 0.0;
        if z0 < sea {
            h += (z1.min(sea) - z0) / self.rise_per_m_submerged();
        }
        if z1 > sea {
            h += (z1 - z0.max(sea)) / self.rise_per_m_subaerial();
        }
        h
    }

    /// Compensated surface after adding `h` m of sediment at surface `z0`.
    pub fn surface_after(&self, z0: f64, h: f64, sea: f64) -> f64 {
        if h <= 0.0 {
            return z0;
        }
        if z0 < sea {
            let to_sea = (sea - z0) / self.rise_per_m_submerged();
            if h <= to_sea {
                return z0 + h * self.rise_per_m_submerged();
            }
            return sea + (h - to_sea) * self.rise_per_m_subaerial();
        }
        z0 + h * self.rise_per_m_subaerial()
    }
}

/// Parameters of `shelf-wave-base-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShelfParams {
    pub wave_base_depth_m: f64,
    pub search_max_cells: usize,
    pub column: CompensatedColumn,
}

impl ShelfParams {
    /// The law's declared values with the caller's isostasy densities.
    pub fn v1(column: CompensatedColumn) -> Self {
        Self {
            wave_base_depth_m: SHELF_WAVE_BASE_DEPTH_M,
            search_max_cells: MARINE_SEARCH_MAX_CELLS,
            column,
        }
    }

    pub fn validate(&self) -> Result<(), SedimentError> {
        self.column.validate()?;
        if !self.wave_base_depth_m.is_finite() || self.wave_base_depth_m < 0.0 {
            return Err(SedimentError::BadParam {
                what: "wave_base_depth_m",
                value: self.wave_base_depth_m,
            });
        }
        if self.search_max_cells == 0 {
            return Err(SedimentError::BadParam {
                what: "search_max_cells",
                value: 0.0,
            });
        }
        Ok(())
    }
}

/// Reporting class of a marine deposit (by water depth before deposition).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MarineClass {
    Shelf,
    Slope,
    Deep,
}

impl MarineClass {
    pub fn of_depth(depth_m: f64) -> Self {
        if depth_m <= SHELF_CLASS_MAX_DEPTH_M {
            MarineClass::Shelf
        } else if depth_m <= SLOPE_CLASS_MAX_DEPTH_M {
            MarineClass::Slope
        } else {
            MarineClass::Deep
        }
    }
}

/// Outcome of one entry's shelf placement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShelfOutcome {
    /// Volume that could not be placed (m^3): the caller's global remainder.
    pub remainder_m3: f64,
    pub bound_hit: bool,
    pub cells_settled: usize,
    /// Cells brought exactly to the shelf top by this placement.
    pub cells_filled: u64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct HeapKey(f64, usize);

impl Eq for HeapKey {}

impl PartialOrd for HeapKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapKey {
    /// Reversed (min-heap on distance, then on element id).
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .0
            .total_cmp(&self.0)
            .then_with(|| other.1.cmp(&self.1))
    }
}

/// Reusable search state (generation-stamped, so a placement costs only the
/// cells it settles).
#[derive(Debug, Clone)]
pub struct ShelfPlacer {
    settled: Vec<u32>,
    generation: u32,
    heap: std::collections::BinaryHeap<HeapKey>,
}

impl ShelfPlacer {
    pub fn new(cells: usize) -> Self {
        Self {
            settled: vec![0; cells],
            generation: 0,
            heap: std::collections::BinaryHeap::new(),
        }
    }

    /// Place `volume_m3` delivered at ocean cell `entry` (see the module
    /// note). `adjacency[c]` lists (neighbour, centre distance m), `marine`
    /// is the receiving ocean, `surf` the compensated surface (m), updated
    /// in place. Each deposit is pushed to `out` as (cell, thickness m,
    /// class). The volumes in `out` plus the remainder equal `volume_m3` to
    /// rounding.
    #[allow(clippy::too_many_arguments)]
    pub fn place(
        &mut self,
        adjacency: &[Vec<(usize, f64)>],
        areas_m2: &[f64],
        marine: &[bool],
        surf: &mut [f64],
        sea_level_m: f64,
        p: &ShelfParams,
        entry: usize,
        volume_m3: f64,
        out: &mut Vec<(usize, f64, MarineClass)>,
    ) -> ShelfOutcome {
        let mut remaining = volume_m3;
        let mut outcome = ShelfOutcome {
            remainder_m3: 0.0,
            bound_hit: false,
            cells_settled: 0,
            cells_filled: 0,
        };
        if remaining <= 0.0 || !marine[entry] {
            outcome.remainder_m3 = remaining.max(0.0);
            return outcome;
        }
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.settled.iter_mut().for_each(|s| *s = 0);
            self.generation = 1;
        }
        let gen = self.generation;
        self.heap.clear();
        self.heap.push(HeapKey(0.0, entry));
        let target = sea_level_m - p.wave_base_depth_m;
        while let Some(HeapKey(d, c)) = self.heap.pop() {
            if self.settled[c] == gen {
                continue;
            }
            if outcome.cells_settled >= p.search_max_cells {
                outcome.bound_hit = true;
                break;
            }
            self.settled[c] = gen;
            outcome.cells_settled += 1;
            if surf[c] < target {
                let need = p.column.thickness_to_raise(surf[c], target, sea_level_m) * areas_m2[c];
                let class = MarineClass::of_depth(sea_level_m - surf[c]);
                if remaining < need {
                    let h = remaining / areas_m2[c];
                    surf[c] = p.column.surface_after(surf[c], h, sea_level_m);
                    out.push((c, h, class));
                    remaining = 0.0;
                    break;
                }
                out.push((c, need / areas_m2[c], class));
                surf[c] = target;
                remaining -= need;
                outcome.cells_filled += 1;
                if remaining <= 0.0 {
                    remaining = 0.0;
                    break;
                }
            }
            // At or above the shelf top: sediment crosses this cell.
            for &(nb, dist) in &adjacency[c] {
                if marine[nb] && self.settled[nb] != gen {
                    self.heap.push(HeapKey(d + dist, nb));
                }
            }
        }
        outcome.remainder_m3 = remaining;
        outcome
    }
}

/// Where the material a step moved came to rest (m^3), for the named
/// control-volume terms. The identity is: bedrock eroded + sediment eroded
/// = the sum of `stored_terms()`. Legacy `deposited_m3` is the four land
/// and lake terms; legacy `exported_m3` is the four marine terms.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DepositionReport {
    pub erosion_law: &'static str,
    pub sediment_law: &'static str,
    /// Net fluvial aggradation on land cells: erosion-deposition, alluvial
    /// storage on fans, floodplains and lowlands.
    pub land_fluvial_m3: f64,
    /// Hillslope-diffusion deposition on land cells.
    pub land_hillslope_m3: f64,
    /// Terminal floors with no outlet (broken topology only).
    pub land_floor_m3: f64,
    /// Lake fill, compensated under `shelf-wave-base-v1`.
    pub lake_m3: f64,
    pub marine_shelf_m3: f64,
    pub marine_slope_m3: f64,
    pub marine_deep_m3: f64,
    /// Uniform layer over the ocean reservoir: under the shelf law, only
    /// what the bounded search could not place; under legacy placement, the
    /// legacy global layer.
    pub global_remainder_m3: f64,
    // ---- diagnostics (not balance terms) ----
    /// Fluvial gross erosion and gross deposition (G Qs / (p A)) summed
    /// over land cells (m^3). Their difference is the net land change.
    pub fluvial_gross_erosion_m3: f64,
    pub fluvial_gross_deposition_m3: f64,
    /// Volume that entered the ocean placement (rivers, overflow, hillslope).
    pub delivered_to_ocean_m3: f64,
    pub marine_search_bound_hits: u64,
    pub lake_spread_bound_hits: u64,
    /// Ocean cells brought exactly to the shelf top this step.
    pub shelf_cells_filled: u64,
    /// Erosion-deposition solver: outer active-set passes, whether the
    /// active set converged, and the largest level change (m) in the last
    /// pass.
    pub sweeps: u32,
    pub converged: bool,
    pub last_sweep_max_change_m: f64,
}

impl DepositionReport {
    /// The named storage terms of the transport identity.
    pub fn stored_terms(&self) -> [(&'static str, f64); 8] {
        [
            ("land_fluvial", self.land_fluvial_m3),
            ("land_hillslope", self.land_hillslope_m3),
            ("land_floor", self.land_floor_m3),
            ("lake", self.lake_m3),
            ("marine_shelf", self.marine_shelf_m3),
            ("marine_slope", self.marine_slope_m3),
            ("marine_deep", self.marine_deep_m3),
            ("global_remainder", self.global_remainder_m3),
        ]
    }

    pub fn total_stored_m3(&self) -> f64 {
        self.stored_terms().iter().map(|(_, v)| v).sum()
    }

    pub fn land_and_lake_m3(&self) -> f64 {
        self.land_fluvial_m3 + self.land_hillslope_m3 + self.land_floor_m3 + self.lake_m3
    }

    pub fn marine_m3(&self) -> f64 {
        self.marine_shelf_m3 + self.marine_slope_m3 + self.marine_deep_m3 + self.global_remainder_m3
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compensated_column_inverts_and_matches_hand_numbers() {
        let k = CompensatedColumn {
            mantle_density_kg_m3: 3300.0,
            sediment_density_kg_m3: 2800.0,
            water_density_kg_m3: 1030.0,
        };
        k.validate().unwrap();
        // -1000 m to -100 m under water: 900 x 2270 / 500 = 4086 m.
        let h = k.thickness_to_raise(-1000.0, -100.0, 0.0);
        assert!((h - 4086.0).abs() < 1e-9, "{h}");
        assert!((k.surface_after(-1000.0, h, 0.0) + 100.0).abs() < 1e-9);
        // Across sea level: -100 to +50: 100 x 4.54 + 50 x 6.6.
        let h2 = k.thickness_to_raise(-100.0, 50.0, 0.0);
        assert!((h2 - (454.0 + 330.0)).abs() < 1e-9, "{h2}");
        assert!((k.surface_after(-100.0, h2, 0.0) - 50.0).abs() < 1e-9);
        assert_eq!(k.thickness_to_raise(5.0, 5.0, 0.0), 0.0);
        let bad = CompensatedColumn {
            sediment_density_kg_m3: 3400.0,
            ..k
        };
        assert!(bad.validate().is_err());
    }

    #[test]
    fn closing_epochs_accumulate_and_validate() {
        let mut l = SedimentLedger::new(2800.0, 1.0e-9, 1.0e-3).unwrap();
        l.record_epoch_m3(2.0, 1000.0, 400.0, 600.0).unwrap();
        l.record_epoch_m3(4.0, 10.0, 10.0, 0.0).unwrap();
        l.record_epoch_m3(6.0, 0.0, 0.0, 0.0).unwrap();
        let b = l.balance().unwrap();
        assert_eq!(b.epochs.len(), 3);
        assert!((b.cumulative_eroded_kg - 1010.0 * 2800.0).abs() < 1.0e-6);
        assert!((b.cumulative_deposited_kg - 410.0 * 2800.0).abs() < 1.0e-6);
        assert!((b.cumulative_exported_kg - 600.0 * 2800.0).abs() < 1.0e-6);
        b.validate().unwrap();
    }

    #[test]
    fn non_closing_epoch_is_refused_with_numbers() {
        let mut l = SedimentLedger::new(2800.0, 1.0e-9, 1.0e-3).unwrap();
        let err = l.record_epoch_m3(2.0, 1000.0, 400.0, 500.0).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("does not close"), "{msg}");
        assert!(msg.contains("residual"), "{msg}");
        assert!(l.epochs().is_empty(), "refused epoch must not be recorded");
        // Rounding-scale residual is inside the declared budget.
        l.record_epoch_m3(2.0, 1000.0, 400.0, 600.0 + 1.0e-10).unwrap();
        // Time must advance.
        assert!(matches!(
            l.record_epoch_m3(2.0, 1.0, 1.0, 0.0).unwrap_err(),
            SedimentError::NonMonotonicTime { .. }
        ));
    }

    #[test]
    fn control_volume_closes_with_named_terms_and_refuses_a_missing_transfer() {
        let mut l = ControlVolumeLedger::new(1.0e-12, 1.0).unwrap();
        // 100 + 30 (ridge) - 20 (eroded) - 5 (subducted) - 3 (drift) = 102.
        let t = ControlVolumeTerms {
            reservoir: "bedrock",
            stage: "step",
            t_my: 2.0,
            initial_m3: 100.0,
            inflows_m3: &[("ridge", 30.0)],
            outflows_m3: &[("eroded", 20.0), ("subducted", 5.0)],
            remap_drift_m3: -3.0,
            final_m3: 102.0,
            diagnostics_m3: &[("global_layer", 7.0)],
        };
        let e = l.record(&t).unwrap();
        assert_eq!(e.residual_m3, 0.0);
        // Drop the subduction transfer: 5 m3 unexplained, above the floor.
        let bad = ControlVolumeTerms {
            outflows_m3: &[("eroded", 20.0)],
            t_my: 4.0,
            ..t
        };
        let err = l.record(&bad).unwrap_err();
        assert!(err.to_string().contains("does not close"), "{err}");
        assert_eq!(l.epochs().len(), 1, "refused epoch must not be recorded");
        // Negative flows are refused (direction lives in the list).
        let neg = ControlVolumeTerms {
            inflows_m3: &[("ridge", -1.0)],
            t_my: 6.0,
            ..t
        };
        assert!(l.record(&neg).is_err());
        let tot = l.totals();
        assert_eq!(tot["bedrock/step/out:eroded"], 20.0);
        assert_eq!(tot["bedrock/step/remap_drift"], -3.0);
        assert_eq!(tot["bedrock/step/diag:global_layer"], 7.0);
    }

    #[test]
    fn tampered_balance_fails_the_contract_gate() {
        let mut l = SedimentLedger::new(2800.0, 1.0e-9, 1.0e-3).unwrap();
        l.record_epoch_m3(2.0, 1000.0, 400.0, 600.0).unwrap();
        let mut b = l.balance().unwrap();
        b.cumulative_exported_kg += 5.0;
        assert!(b.validate().is_err(), "cumulative tamper must be caught");
        let mut b2 = l.balance().unwrap();
        b2.epochs[0].deposited_kg += 5.0;
        assert!(b2.validate().is_err(), "epoch tamper must be caught");
    }
}
