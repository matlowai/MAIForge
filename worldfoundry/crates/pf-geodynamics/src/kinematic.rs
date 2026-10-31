//! Eulerian plate-field evolution for the PREVIEW tier (ADR 032).
//!
//! Plates rotate rigidly about their seeded Euler poles (`plates.rs`,
//! `rotation.rs`, both built). Each step, every cell asks every plate
//! "would your material have arrived here?" by back-rotating the cell
//! centroid through that plate's motion and finding the cell whose
//! MATERIAL POSITION is nearest that point (semi-Lagrangian carry). Every
//! cell stores the true position of the material it currently represents
//! (`material_pos`, its sub-cell phase), so displacement smaller than a
//! cell accumulates across steps instead of being rounded away: the field
//! shifts by one cell once the phase passes half a cell, and the mean
//! speed is right at any step size (pinned by `slow_motion_accumulates`).
//!
//! - exactly one claimant: the cell carries that source's state;
//! - no claimant (the plates parted): RIDGE BIRTH, new oceanic crust at
//!   age zero, joining the plate that previously occupied the cell;
//! - several claimants (the plates met): CONVERGENCE, resolved by the
//!   declared overriding rule below, with the loser subducted (oceanic)
//!   or stacked (continental, `uplift.rs`).
//!
//! Boundary CLASSIFICATION per cell uses the built `boundaries.rs` kernels
//! (`normal_rate_at`, `classify`) with the local direction toward the
//! neighbouring plate as the boundary normal.
//!
//! HONEST LIMITS (also stated in ADR 032): the carry is not conservative
//! (a rigid rotation does not map lat/lon cells one to one), so the crust
//! volume ledger REPORTS its advection residual rather than gating on it;
//! the Lagrangian parcel ledger of `parcels.rs` (PARCEL != CELL) is not
//! driven here, so this tier publishes no material lineage. Both return
//! with the solved-mantle tier.
//!
//! Determinism: no randomness, no HashMap iteration, ascending-cell order
//! everywhere, ties break to the lowest plate id (declared).

use crate::boundaries::{classify, normal_rate_at, BoundariesError};
use crate::carry::{in_cell_patch, CarryLaw};
use crate::plates::PlateSeeding;
use crate::rotation::{
    great_circle_rad, lat_lon_to_unit, rotate_by_euler, unit_to_lat_lon, RotationError,
};
use crate::uplift::UpliftPrior;
use pf_contracts::domain::{BoundaryKindV0, EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_grid::GridSpec;
use std::collections::BTreeMap;

/// SEED-PERF (lane SP, 2026-10-13): deterministic per-cell parallelism.
///
/// `map_indexed(n, f)` returns exactly `(0..n).map(f).collect()`: work is
/// split into fixed index chunks, each chunk's results land in its own slot
/// and the slots are concatenated in index order, so the output (and any
/// sum a caller then folds over it in index order) does not depend on the
/// thread count or on scheduling. Callers keep every order-dependent
/// reduction serial over the returned slots.
///
/// The thread count is the speed dial: `set_threads(n)` (0 = unset), else
/// the `PF_THREADS` environment variable, else the machine's available
/// parallelism; 1 runs everything serially on the calling thread.
pub mod par {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Condvar, Mutex, OnceLock};

    static OVERRIDE: AtomicUsize = AtomicUsize::new(0);
    static ENV_DEFAULT: OnceLock<usize> = OnceLock::new();

    /// Cap (or uncap) the worker count for this process (0 = back to the
    /// environment / machine default).
    pub fn set_threads(n: usize) {
        OVERRIDE.store(n, Ordering::Relaxed);
    }

    /// Effective worker count (>= 1).
    pub fn threads() -> usize {
        let o = OVERRIDE.load(Ordering::Relaxed);
        if o > 0 {
            return o;
        }
        *ENV_DEFAULT.get_or_init(|| {
            std::env::var("PF_THREADS")
                .ok()
                .and_then(|v| v.trim().parse::<usize>().ok())
                .filter(|&v| v > 0)
                .unwrap_or_else(|| std::thread::available_parallelism().map(|v| v.get()).unwrap_or(1))
        })
    }

    // A small persistent worker pool (std only): `execute(t, task)` runs
    // `task` on the calling thread and t - 1 pool workers at once and
    // returns when all of them have returned from it. The task pulls
    // fixed index chunks from an atomic counter, so which thread ran a
    // chunk never shows in the output. One job at a time: a call made
    // while another job is in flight (another test thread, or a nested
    // call from inside a task) runs serially instead, which gives the same
    // output.
    type Task = dyn Fn() + Sync;
    struct State {
        gen: u64,
        task: Option<*const Task>,
        wanted: usize,
        joined: usize,
        finished: usize,
        workers: usize,
        panicked: bool,
    }
    // The task pointer is only dereferenced while its job runs; the caller
    // blocks until every joined worker has finished with it.
    unsafe impl Send for State {}
    struct Pool {
        state: Mutex<State>,
        work: Condvar,
        done: Condvar,
    }
    static POOL: OnceLock<Pool> = OnceLock::new();
    static JOB: Mutex<()> = Mutex::new(());

    fn pool() -> &'static Pool {
        POOL.get_or_init(|| Pool {
            state: Mutex::new(State {
                gen: 0,
                task: None,
                wanted: 0,
                joined: 0,
                finished: 0,
                workers: 0,
                panicked: false,
            }),
            work: Condvar::new(),
            done: Condvar::new(),
        })
    }

    fn worker(p: &'static Pool) {
        let mut seen = 0u64;
        loop {
            let task = {
                let mut st = p.state.lock().unwrap();
                loop {
                    if st.gen != seen {
                        seen = st.gen;
                        if st.joined < st.wanted {
                            st.joined += 1;
                            break st.task.expect("job without task");
                        }
                    }
                    st = p.work.wait(st).unwrap();
                }
            };
            // SAFETY: the job's caller keeps the task alive until
            // `finished == wanted`, which this worker reports below.
            let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| unsafe { (*task)() }));
            let mut st = p.state.lock().unwrap();
            st.finished += 1;
            st.panicked |= r.is_err();
            if st.finished == st.wanted {
                p.done.notify_all();
            }
        }
    }

    fn execute(t: usize, task: &(dyn Fn() + Sync)) {
        let guard = if t > 1 { JOB.try_lock().ok() } else { None };
        let Some(_guard) = guard else {
            task();
            return;
        };
        let p = pool();
        let wanted = t - 1;
        {
            let mut st = p.state.lock().unwrap();
            while st.workers < wanted {
                st.workers += 1;
                std::thread::Builder::new()
                    .name("pf-par".into())
                    .spawn(move || worker(p))
                    .expect("spawn pf-par worker");
            }
            // SAFETY: lifetime erased only for the duration of this job;
            // we wait below for every worker that joined it.
            let ptr: *const Task = unsafe { std::mem::transmute::<&(dyn Fn() + Sync), &'static Task>(task) };
            st.task = Some(ptr);
            st.wanted = wanted;
            st.joined = 0;
            st.finished = 0;
            st.panicked = false;
            st.gen = st.gen.wrapping_add(1);
        }
        p.work.notify_all();
        let mine = std::panic::catch_unwind(std::panic::AssertUnwindSafe(task));
        let panicked = {
            let mut st = p.state.lock().unwrap();
            while st.finished < st.wanted {
                st = p.done.wait(st).unwrap();
            }
            st.task = None;
            st.wanted = 0;
            st.panicked
        };
        if let Err(e) = mine {
            std::panic::resume_unwind(e);
        }
        if panicked {
            panic!("a pf-par worker panicked");
        }
    }

    /// Below this many items a call runs serially (dispatch costs more
    /// than it saves). Output is identical either way.
    const MIN_PARALLEL: usize = 2048;

    /// `(0..n).map(f).collect()`, computed on up to `threads()` workers
    /// (for cheap per-cell items: chunks of at least 256).
    pub fn map_indexed<T: Send, F: Fn(usize) -> T + Sync>(n: usize, f: F) -> Vec<T> {
        if n < MIN_PARALLEL {
            return (0..n).map(f).collect();
        }
        map_grain(n, 256, f)
    }

    /// `map_indexed` for few, expensive items (discs, fronts): parallel
    /// from 2 items, chunks of at least one item.
    pub fn map_heavy<T: Send, F: Fn(usize) -> T + Sync>(n: usize, f: F) -> Vec<T> {
        map_grain(n, 1, f)
    }

    fn map_grain<T: Send, F: Fn(usize) -> T + Sync>(n: usize, grain: usize, f: F) -> Vec<T> {
        let t = threads().min(n / grain).max(1);
        if t == 1 {
            return (0..n).map(f).collect();
        }
        // Several chunks per worker so uneven items balance; chunk
        // boundaries never affect the result.
        let chunk = (n / (t * 8)).max(grain);
        let nchunks = n.div_ceil(chunk);
        let next = AtomicUsize::new(0);
        let slots: Mutex<Vec<Option<Vec<T>>>> = Mutex::new((0..nchunks).map(|_| None).collect());
        execute(t, &|| loop {
            let k = next.fetch_add(1, Ordering::Relaxed);
            if k >= nchunks {
                break;
            }
            let lo = k * chunk;
            let hi = (lo + chunk).min(n);
            let v: Vec<T> = (lo..hi).map(&f).collect();
            slots.lock().unwrap()[k] = Some(v);
        });
        let mut out = Vec::with_capacity(n);
        for v in slots.into_inner().unwrap() {
            out.extend(v.expect("every chunk ran"));
        }
        out
    }

    /// Fill `out[i] = f(i)` in parallel (same contract as `map_indexed`).
    pub fn fill_indexed<T: Send, F: Fn(usize) -> T + Sync>(out: &mut [T], f: F) {
        for_chunks_mut(out, |lo, c| {
            for (i, o) in c.iter_mut().enumerate() {
                *o = f(lo + i);
            }
        });
    }

    /// Run `f(lo, chunk)` over disjoint chunks of `out` in parallel (chunks
    /// of at least 256; serial below `MIN_PARALLEL`). `f` must make each
    /// element depend only on its own index and on data it does not write.
    pub fn for_chunks_mut<T: Send, F: Fn(usize, &mut [T]) + Sync>(out: &mut [T], f: F) {
        let n = out.len();
        let t = threads().min(n / 256).max(1);
        if t == 1 || n < MIN_PARALLEL {
            f(0, out);
            return;
        }
        let chunk = (n / (t * 8)).max(256);
        let chunks: Mutex<Vec<(usize, &mut [T])>> =
            Mutex::new(out.chunks_mut(chunk).enumerate().map(|(k, c)| (k * chunk, c)).rev().collect());
        execute(t, &|| loop {
            let Some((lo, c)) = chunks.lock().unwrap().pop() else { break };
            f(lo, c);
        });
    }

    /// `f(lo, hi)` over consecutive index chunks covering `0..n`, results
    /// in chunk order. The chunking depends on `threads()`, so callers may
    /// only combine the results with order-free folds (max, min, any) or
    /// concatenate per-index outputs in order.
    pub fn map_chunks<R: Send, F: Fn(usize, usize) -> R + Sync>(n: usize, f: F) -> Vec<R> {
        let t = threads().min(n / 256).max(1);
        if t == 1 || n < MIN_PARALLEL {
            return vec![f(0, n)];
        }
        let chunk = (n / (t * 8)).max(256);
        map_heavy(n.div_ceil(chunk), |k| f(k * chunk, ((k + 1) * chunk).min(n)))
    }
}

/// SEED-PERF (lane SP): the 8 neighbour cell ids of `c`, in exactly the
/// order of `GridSpec::neighbors` (di, dj in -1..=1, no rows past the
/// poles, longitude wraps) but without its 8 haversine distances, for the
/// loops that only use the ids. Returns the count written.
#[inline]
pub(crate) fn neighbor_ids(nlat: usize, nlon: usize, c: usize, out: &mut [usize; 8]) -> usize {
    let (i, j) = ((c / nlon) as i64, (c % nlon) as i64);
    let mut m = 0;
    for di in [-1i64, 0, 1] {
        let ti = i + di;
        if ti < 0 || ti >= nlat as i64 {
            continue;
        }
        for dj in [-1i64, 0, 1] {
            if di == 0 && dj == 0 {
                continue;
            }
            out[m] = (ti * nlon as i64 + (j + dj).rem_euclid(nlon as i64)) as usize;
            m += 1;
        }
    }
    m
}

#[derive(Debug, Clone, PartialEq)]
pub enum KinematicError {
    SizeMismatch {
        what: &'static str,
        expected: usize,
        found: usize,
    },
    Rotation(RotationError),
    Boundaries(BoundariesError),
    NonFinite {
        what: &'static str,
        value: f64,
        index: usize,
    },
    UnknownPlate {
        plate: PlateId,
    },
    Grid(String),
    /// A subduction overlap whose convergence normal cannot be sampled
    /// from local geometry (SEED-ROBUST, typed loud refusal).
    ArcContact {
        cell: usize,
        winner: PlateId,
        loser: PlateId,
        why: ArcContactRefusal,
    },
}

/// Why an arc overlap's pair normal could not be sampled. Neither case is
/// resolved by reaching for a remote boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum ArcContactRefusal {
    /// No step-start interface of the pair lies within the swept locality
    /// bound AND the two claimant sources themselves are farther apart than
    /// that bound: the material that met is not local to the overlap.
    NonlocalContact {
        nearest_interface_rad: Option<f64>,
        source_separation_rad: f64,
        bound_rad: f64,
    },
    /// The new-contact segment between the two sources has no usable
    /// midpoint or tangent normal (coincident or antipodal sources).
    DegenerateContact { source_separation_rad: f64 },
}

/// How an arc overlap's convergence normal was sampled. Aligned one to
/// one with `AdvanceLedger::arc_normals` (SEED-ROBUST, lane R5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArcContactGeometry {
    /// The nearest step-start interface edge of the pair, within the swept
    /// locality bound (the technical review D5 rule, unchanged). `input_interface` is
    /// that edge (two grid neighbours, lower plate id first).
    StepStartInterface,
    /// The pair had no step-start interface within the bound, so the two
    /// plates came into contact DURING the step: a separating sliver of a
    /// third plate was consumed or vacated, or a detached fragment reached
    /// the other plate. The normal is taken across the material that met:
    /// the segment between the two claimants' step-start source cells,
    /// sampled at its midpoint. `input_interface` is that source pair
    /// (lower plate id first; generally NOT grid neighbours).
    NewContactBridge,
}

impl std::fmt::Display for KinematicError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KinematicError::SizeMismatch {
                what,
                expected,
                found,
            } => write!(f, "{what}: expected {expected}, found {found}"),
            KinematicError::Rotation(e) => write!(f, "rotation: {e}"),
            KinematicError::Boundaries(e) => write!(f, "boundaries: {e}"),
            KinematicError::NonFinite { what, value, index } => {
                write!(f, "nonfinite {what} {value:?} at {index}")
            }
            KinematicError::UnknownPlate { plate } => write!(f, "unknown plate {plate:?}"),
            KinematicError::Grid(why) => write!(f, "grid: {why}"),
            KinematicError::ArcContact {
                cell,
                winner,
                loser,
                why,
            } => match why {
                ArcContactRefusal::NonlocalContact {
                    nearest_interface_rad,
                    source_separation_rad,
                    bound_rad,
                } => write!(
                    f,
                    "arc contact {cell}: plates {} / {}: nonlocal contact (nearest step-start interface {}, claimant sources {source_separation_rad} rad apart, bound {bound_rad} rad)",
                    winner.0,
                    loser.0,
                    nearest_interface_rad
                        .map(|d| format!("{d} rad"))
                        .unwrap_or_else(|| "absent".into())
                ),
                ArcContactRefusal::DegenerateContact {
                    source_separation_rad,
                } => write!(
                    f,
                    "arc contact {cell}: plates {} / {}: degenerate new-contact segment (sources {source_separation_rad} rad apart)",
                    winner.0, loser.0
                ),
            },
        }
    }
}

impl std::error::Error for KinematicError {}

impl From<RotationError> for KinematicError {
    fn from(e: RotationError) -> Self {
        KinematicError::Rotation(e)
    }
}

impl From<BoundariesError> for KinematicError {
    fn from(e: BoundariesError) -> Self {
        KinematicError::Boundaries(e)
    }
}

/// Per-cell Eulerian crust state. Every vector has `grid.cells()` entries.
#[derive(Debug, Clone, PartialEq)]
pub struct CrustState {
    pub t_my: f64,
    pub plate: Vec<PlateId>,
    pub kind: Vec<CrustKindCodeV0>,
    /// Birth time (My); age = t_my - birth. Negative for initial crust.
    pub birth_t_my: Vec<f64>,
    /// Crust column thickness (m), bedrock only; sediment is separate.
    pub thickness_m: Vec<f64>,
    /// Sediment column thickness (m) on top of the crust.
    pub sediment_m: Vec<f64>,
    /// Cumulative bedrock removed by erosion (m) over the cell's history.
    pub eroded_total_m: Vec<f64>,
    /// Cumulative arc magmatic addition (m).
    pub arc_added_m: Vec<f64>,
    /// Cumulative collisional thickening received (m).
    pub shortening_added_m: Vec<f64>,
    /// Last classification (code) and signed normal rate (m/My, positive
    /// convergent) at the cell.
    pub boundary_class: Vec<u8>,
    pub normal_rate_m_per_my: Vec<f64>,
    /// True position (unit vector) of the material this cell represents;
    /// within about half a cell of the centroid. Seeded at the centroid,
    /// rotated rigidly with the plate, reset to the centroid at births.
    pub material_pos: Vec<[f64; 3]>,
}

/// One continental body's seeding bookkeeping (opt-in, lane K2).
#[derive(Debug, Clone, PartialEq)]
pub struct ContinentalBodyRecord {
    pub plate: PlateId,
    /// The host plate's Voronoi seed point (unit vector): the body centre.
    pub seed_point: [f64; 3],
    pub target_area_m2: f64,
    pub achieved_area_m2: f64,
    pub cells: u64,
    /// Sum of thickness * area over the body (m^3).
    pub crust_volume_m3: f64,
    /// Body cells thinner than the interior thickness (inside the taper).
    pub tapered_cells: u64,
}

/// Seeding bookkeeping for `CrustState::seed_bodies` (opt-in, lane K2):
/// declared vs achieved area and the crust volume the taper removes
/// relative to plateau bodies of the same area.
#[derive(Debug, Clone, PartialEq)]
pub struct ContinentalBodiesReport {
    pub declared_area_fraction: f64,
    pub achieved_area_fraction: f64,
    pub rifted_margin_width_m: f64,
    pub interior_thickness_m: f64,
    pub edge_thickness_m: f64,
    pub continental_crust_volume_m3: f64,
    /// achieved body area * interior thickness (what plateau bodies hold).
    pub plateau_equivalent_volume_m3: f64,
    pub coast_midpoints: u64,
    /// Per cell: distance (m) from the centroid to the nearest coast edge
    /// midpoint for body cells; NaN elsewhere.
    pub coast_distance_m: Vec<f64>,
    pub bodies: Vec<ContinentalBodyRecord>,
}

/// Per-step crust volume accounting (m^3) produced by `advance`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AdvanceLedger {
    pub wedge_in_m3: f64,
    pub wedge_stacked_m3: f64,
    pub wedge_delaminated_m3: f64,
    pub arc_normals: Vec<pf_contracts::domain::ArcNormalSampleV0>,
    /// How each `arc_normals` sample was taken (same length and order).
    pub arc_contact_geometry: Vec<ArcContactGeometry>,
    pub ridge_created_m3: f64,
    pub arc_created_m3: f64,
    pub subducted_m3: f64,
    pub delaminated_m3: f64,
    pub volume_before_m3: f64,
    pub volume_after_m3: f64,
    /// Sediment retired with subducted columns (m^3; lane E proposal,
    /// applied by lane K2: sediment subduction, named, not hidden).
    pub sediment_subducted_m3: f64,
    /// Sediment of collision losers stacked onto winners (m^3).
    pub sediment_stacked_m3: f64,
    pub sediment_before_m3: f64,
    pub sediment_after_m3: f64,
    pub ridge_cells: u64,
    pub subduction_cells: u64,
    pub collision_cells: u64,
    /// SEED-OROGEN `rate-limited-belt-v1` only (zero under the legacy rule):
    /// loser crust consumed by collision overlaps this step (m^3), the part
    /// placed on belts (equal to rounding), loser sediment moved onto belts,
    /// the spreading report, and the backstop delamination (also counted
    /// in `delaminated_m3`).
    pub belt_in_m3: f64,
    /// Loser crust in continent-continent overlaps (m^3), left to the remap
    /// under `rate-limited-belt-v1` (diagnostic: overlap counting is not a
    /// dt-converged supply measure, lane OG falsifier F2).
    pub collision_overlap_consumed_m3: f64,
    /// Step-start continent-continent interface edges with convergence.
    pub collision_supply_edges: u64,
    pub belt_placed_m3: f64,
    pub belt_sediment_m3: f64,
    pub belt_delaminated_m3: f64,
    pub spread: Option<crate::orogen::SpreadReport>,
    /// Arc volume requested by convergent samples and the part placed.
    pub arc_placed_m3: f64,
    /// Thickest crust column after the step (m), any law.
    pub max_thickness_m: f64,
    /// SEDIMENT-TOWERS (lane ST) diagnostic, any law, feeds no output: the
    /// largest sediment thickness (m) one cell received this step from
    /// collision moves (legacy stacking or the belt's sediment placement),
    /// and that cell.
    pub sediment_moved_max_m: f64,
    pub sediment_moved_max_cell: u64,
    /// `accreted-to-belt-crust-v1` only (zero otherwise): collision-loser
    /// sediment turned into belt crust supply this step (m^3). It leaves
    /// the sediment field and enters the crust volume, so both residuals
    /// below name it.
    pub sediment_accreted_m3: f64,
    /// SEED-CARRY-CC diagnostics. Computed under every carry law (they do
    /// not feed any output); the `cc_*` supply terms are nonzero only under
    /// `front-rescue-carry-v1`.
    pub carry: CarryDiagnostics,
    /// SEED-EXTRAS (lane EX) `fixed-plume-chains-v1` only (zero otherwise):
    /// magmatic crust added by plume heads and tails this step (m^3; heads
    /// included in `hotspot_created_m3`, named in the residual), and the
    /// per-cell increments (cell, m, latest eruption time My).
    pub hotspot_created_m3: f64,
    pub hotspot_head_m3: f64,
    pub hotspot_heads: u64,
    pub hotspot_cells: Vec<(u32, f64, f64)>,
    /// SEED-FLOODBASALT (lane FB) `plume-head-lip-swell-v1` only (zero
    /// otherwise): LIP crust erupted this step (m^3; also inside
    /// `hotspot_created_m3` and `hotspot_head_m3`), the part placed on
    /// columns continental at eruption, and the per-cell increments (cell,
    /// m, eruption midpoint My).
    pub lip_created_m3: f64,
    pub lip_continental_m3: f64,
    pub lip_cells: Vec<(u32, f64, f64)>,
    /// SEED-EXTRAS (lane EX) `stretched-rift-v1` only: births filled with
    /// stretched continental crust (taken out of `ridge_created_m3` and
    /// `ridge_cells`), candidates left oceanic at breakup, candidates of
    /// converging plates skipped, the volume moved from flanks into births
    /// (m^3), and the per-cell flank thinning (cell, m).
    pub rift_births: Vec<u32>,
    pub rift_breakup_births: u64,
    /// SEED-BREAKUP-CUT-V2 (lane BC) `symmetric-ridge-accretion-v1`: gap
    /// cells moved to the other flank this step (0 under the legacy rule).
    pub ridge_relabelled_cells: u64,
    pub rift_convergent_skipped: u64,
    pub rift_moved_m3: f64,
    pub rift_thinning: Vec<(u32, f64)>,
    /// SEED-RIFT-SHOULDERS (lane RS): per rift birth (parallel to
    /// `rift_births`), the weighted mean pre-stretch flank thickness (m).
    pub rift_birth_pre_thickness_m: Vec<f64>,
    /// SEED-RIFT-SHOULDERS (lane RS): per rift birth, its flank cells and
    /// weights (triangular weight x area).
    pub rift_birth_flanks: Vec<Vec<(u32, f64)>>,
    /// SEED-CARRY-CONSERVE (lane CV) `continental-parcel-carry-v1` only
    /// (None otherwise): the parcel step's bookkeeping and the projection
    /// the advance wrote into continental cells (for the write-back).
    pub parcels: Option<crate::conserve::ParcelStepReport>,
    pub parcel_projection: Option<crate::conserve::ParcelProjection>,
    /// SEED-WORLD-COMPOSE (lane CP) `plateau-cap-spill-v2` only: the
    /// in-belt plateau pass of this step.
    pub plateau_belt: Option<crate::ranges::PlateauBeltReport>,
    /// SEED-WORLD-COMPOSE (lane CP, lane VL defect-4 diagnostic): per-cell
    /// belt budget of this step under `rate-limited-belt-v1` (None under the
    /// legacy rule). Read-only; feeds no output.
    pub belt_budget: Option<BeltBudget>,
    /// SEED-TRANSFORMS (lane TF) `smoothed-normal-transform-v1` only (zero
    /// otherwise): the step's gap fills and non-convergent overlaps.
    pub transform: TransformDiagnostics,
}

/// SEED-TRANSFORMS (lane TF) per-step diagnostics of
/// `smoothed-normal-transform-v1`. A transform neither creates nor
/// destroys crust: the filled columns and the retired non-convergent
/// overlap losers are remap artefacts of the semi-Lagrangian carry, so
/// neither is booked as ridge creation or subduction; both land in the
/// named advection residual (`AdvanceLedger::residual_m3`) and are named
/// here.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TransformDiagnostics {
    /// Gap cells filled from same-side material, by the step-start class
    /// found at the gap (transform; subduction or collision).
    pub filled_transform_cells: u64,
    pub filled_convergent_cells: u64,
    /// Volume of the filled columns (m^3).
    pub filled_m3: f64,
    /// Gap cells on a transform or convergent class left as ridge births
    /// (no carried neighbour within the fill layers).
    pub unfilled_cells: u64,
    /// Oceanic overlap losers on a transform or divergent class, retired
    /// without subduction (no arc sample), and their volume (m^3).
    pub nonconvergent_overlaps: u64,
    pub nonconvergent_overlap_m3: f64,
    /// Boundary cells (step end) whose smoothed normal was degenerate or
    /// pointed away from the other plate (staircase normal used instead).
    pub normal_fallback_cells: u64,
    pub boundary_cells: u64,
    /// SEED-OLD-OCEAN (lane OR, 2026-10-13; lane VL defect 2, second half):
    /// read-only logs, feed no output. Subducted loser area (m^2) by the
    /// step-start local class at the overlap (index = boundary class code:
    /// 0 interior, 1 ridge, 2 subduction, 3 collision, 4 transform) and by
    /// loser crust age (`SUBDUCTED_AGE_BINS_MY` lower edges), and the ridge
    /// births left after the gap fill (m^2) by the step-start local class.
    pub subducted_by_class_age_m2: [[f64; 7]; 5],
    pub births_by_class_m2: [f64; 5],
}

/// SEED-OLD-OCEAN (lane OR): lower edges (My) of the loser-age bins of
/// `TransformDiagnostics::subducted_by_class_age_m2` (last bin open).
pub const SUBDUCTED_AGE_BINS_MY: [f64; 7] = [0.0, 10.0, 20.0, 30.0, 60.0, 100.0, 150.0];

fn age_bin(age_my: f64) -> usize {
    let mut b = 0;
    for (i, lo) in SUBDUCTED_AGE_BINS_MY.iter().enumerate() {
        if age_my >= *lo {
            b = i;
        }
    }
    b
}

/// SEED-WORLD-COMPOSE (lane CP): where the belt step's shortening volume
/// went, per cell (m of crust thickness): collision supply placed, net
/// change from the plateau pass (v2: accretion + spill + taper), net
/// change from the channel-flow spreading, and the delamination backstop.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BeltBudget {
    pub placed_m: Vec<f64>,
    pub plateau_dh_m: Vec<f64>,
    pub spread_dh_m: Vec<f64>,
    pub delaminated_m: Vec<f64>,
}

/// SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13; PREREG revision 1): which
/// plate new crust born in a divergence gap joins.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RidgeAccretion {
    /// LEGACY: the gap cell joins its previous occupant (the plate that just
    /// moved off it), so a young ocean accretes to one plate and the
    /// boundary stays at the other plate's edge.
    #[default]
    PreviousOccupantV0,
    /// `symmetric-ridge-accretion-v1` (selected by plate-reorg-scripted-v2):
    /// gap cells are split between the flanking plates so the boundary sits
    /// mid-ridge. Gap cells touching carried crust join the adjacent plate
    /// whose youngest adjacent carried crust is OLDEST (the flank that
    /// accreted least recently; continental crust counts as oldest); ties
    /// keep the previous occupant if it is among them, else the lowest
    /// plate id. Deeper gap cells take the plate most of their already
    /// labelled gap neighbours carry (ties as above), layer by layer. Runs
    /// after the rift and plume passes, before the speckle fold; parcels in
    /// relabelled cells follow.
    SymmetricV1,
}

/// SEED-EXTRAS (lane EX): the opt-in passes `advance_with_extras` runs
/// after the carry (plumes) and on the step's births (rift). `None` fields
/// are the legacy laws.
#[derive(Debug, Clone, Default)]
pub struct ExtrasLaws<'a> {
    pub plumes: Option<(&'a [crate::hotspots::Plume], crate::hotspots::PlumeParams)>,
    pub rift: Option<crate::rift::RiftParams>,
    /// SEED-FLOODBASALT (lane FB): `plume-head-lip-swell-v1` LIP crust.
    /// Requires `plumes`; supersedes their v1 head (run with head volume 0).
    pub flood_basalt: Option<crate::flood_basalt::LipParams>,
}

/// SEED-CARRY-CC (lane CC) per-step carry diagnostics.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CarryDiagnostics {
    /// Rescue claims admitted (continental parcels the legacy carry would
    /// have carried nowhere).
    pub rescues_admitted: u64,
    /// Continental losers recorded in continent-continent overlaps.
    pub cc_overlap_losers: u64,
    /// Of those, losers still carried to some destination (not consumed).
    pub cc_losers_carried: u64,
    /// Loser parcels consumed (carried nowhere; counted once each), their
    /// source-cell area (m^2), crust volume (m^3) and sediment (m^3).
    pub cc_consumed_parcels: u64,
    pub cc_consumed_area_m2: f64,
    pub cc_consumed_m3: f64,
    pub cc_consumed_sediment_m3: f64,
    /// Continental source parcels carried nowhere that were NOT a
    /// continent-continent overlap loser: silent drops (any law).
    pub continental_dropped_parcels: u64,
    pub continental_dropped_area_m2: f64,
    pub continental_dropped_m3: f64,
    /// Ridge births whose 8-neighbourhood (after the step) holds continental
    /// cells of at least two different plates.
    pub sliver_births: u64,
}

impl AdvanceLedger {
    /// Sediment remap drift: after - (before - subducted). Stacked
    /// sediment is a move between columns, not creation; if the loser's
    /// source column is ALSO carried elsewhere it is duplicated, and that
    /// shows up here as positive drift. The area factor follows the
    /// subducted_m3 convention (destination cell area).
    pub fn sediment_residual_m3(&self) -> f64 {
        self.sediment_after_m3
            - (self.sediment_before_m3 - self.sediment_subducted_m3 - self.sediment_accreted_m3)
    }

    /// Advection residual for this step: after - (before + created - removed).
    pub fn residual_m3(&self) -> f64 {
        self.volume_after_m3
            - (self.volume_before_m3 + self.ridge_created_m3 + self.arc_created_m3
                + self.sediment_accreted_m3
                + self.hotspot_created_m3
                - self.subducted_m3
                - self.delaminated_m3)
    }
}

pub(crate) fn centroids(grid: &GridSpec) -> Result<Vec<[f64; 3]>, KinematicError> {
    // SEED-PERF (lane SP): per cell in parallel; the first error in cell
    // order is returned, as from the serial row-major loop.
    let nlon = grid.nlon();
    par::map_indexed(grid.cells(), |c| -> Result<[f64; 3], KinematicError> {
        let lat = grid
            .center_lat_deg(c / nlon)
            .map_err(|e| KinematicError::Grid(e.to_string()))?;
        let lon = grid
            .center_lon_deg(c % nlon)
            .map_err(|e| KinematicError::Grid(e.to_string()))?;
        Ok(lat_lon_to_unit(lat, lon)?)
    })
    .into_iter()
    .collect()
}

/// Nearest grid cell to a unit vector (declared: floor of the fractional
/// row/column index; longitude wraps).
pub fn nearest_cell(grid: &GridSpec, p: [f64; 3]) -> Result<usize, KinematicError> {
    let (lat, lon) = unit_to_lat_lon(p)?;
    let dlat = 180.0 / grid.nlat() as f64;
    let dlon = 360.0 / grid.nlon() as f64;
    let i = (((lat + 90.0) / dlat).floor() as i64).clamp(0, grid.nlat() as i64 - 1) as usize;
    let j = (((lon + 180.0) / dlon).floor() as i64).rem_euclid(grid.nlon() as i64) as usize;
    Ok(i * grid.nlon() + j)
}

pub(crate) fn cell_areas(grid: &GridSpec) -> Result<Vec<f64>, KinematicError> {
    let mut out = Vec::with_capacity(grid.cells());
    for i in 0..grid.nlat() {
        let a = grid
            .cell_area_m2(i)
            .map_err(|e| KinematicError::Grid(e.to_string()))?;
        for _ in 0..grid.nlon() {
            out.push(a);
        }
    }
    Ok(out)
}

impl CrustState {
    /// Initial state from a plate seeding: kind from the continental plate
    /// set, thickness and age from the prior.
    pub fn seed(
        grid: &GridSpec,
        seeding: &PlateSeeding,
        prior: &UpliftPrior,
    ) -> Result<Self, KinematicError> {
        let n = grid.cells();
        if seeding.assignment.len() != n {
            return Err(KinematicError::SizeMismatch {
                what: "plate assignment",
                expected: n,
                found: seeding.assignment.len(),
            });
        }
        let mut kind = Vec::with_capacity(n);
        let mut birth = Vec::with_capacity(n);
        let mut thick = Vec::with_capacity(n);
        for pid in &seeding.assignment {
            if seeding.continental_plate_ids.contains(pid) {
                kind.push(CrustKindCodeV0::Continental);
                birth.push(-prior.initial_continental_age_my);
                thick.push(prior.record.continental_crust_thickness_m);
            } else {
                kind.push(CrustKindCodeV0::Oceanic);
                birth.push(-prior.initial_oceanic_age_my);
                thick.push(prior.record.oceanic_crust_thickness_m);
            }
        }
        let material_pos = centroids(grid)?;
        Ok(Self {
            t_my: 0.0,
            plate: seeding.assignment.clone(),
            kind,
            birth_t_my: birth,
            thickness_m: thick,
            sediment_m: vec![0.0; n],
            eroded_total_m: vec![0.0; n],
            arc_added_m: vec![0.0; n],
            shortening_added_m: vec![0.0; n],
            boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
            normal_rate_m_per_my: vec![0.0; n],
            material_pos,
        })
    }

    /// OPT-IN (CRUST-SPIKE, lane K2): initial state with continental crust
    /// as BODIES carried by plates instead of whole continental plates.
    ///
    /// Named approximations (declared in the artifact by the driver):
    /// - `geodesic-disc-body-v0`: each plate in
    ///   `seeding.continental_plate_ids` hosts one body: its own cells taken
    ///   in ascending great-circle distance from the plate's Voronoi seed
    ///   point (ties: lowest cell id) until the body's area target is met
    ///   to within half of the last cell's area. Target of body k =
    ///   `area_fraction * planet area * A_plate_k / sum(A_continental
    ///   plates)`. A body larger than its plate is REFUSED (loud), never
    ///   grown into a neighbour or silently shrunk.
    /// - `linear-rifted-margin-taper-v0` (`UpliftPrior::
    ///   rifted_margin_thickness_m`) on the distance from each body cell
    ///   centroid to the nearest coast edge midpoint (the normalized sum of
    ///   the centroids of a body cell and an adjacent non-body cell, over
    ///   the 8-neighbourhood of the global continental mask).
    ///
    /// Every other cell (including the rest of each continental plate) is
    /// oceanic at the prior's oceanic thickness and initial oceanic age.
    /// No randomness: the plate seeding is read, never redrawn, so the
    /// plate stream is untouched.
    pub fn seed_bodies(
        grid: &GridSpec,
        seeding: &PlateSeeding,
        prior: &UpliftPrior,
        area_fraction: f64,
        rifted_margin_width_m: f64,
    ) -> Result<(Self, ContinentalBodiesReport), KinematicError> {
        let n = grid.cells();
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
        if !rifted_margin_width_m.is_finite() || rifted_margin_width_m < 0.0 {
            return Err(KinematicError::NonFinite {
                what: "rifted_margin_width_m",
                value: rifted_margin_width_m,
                index: 0,
            });
        }
        let cents = centroids(grid)?;
        let areas = cell_areas(grid)?;
        let total_area: f64 = areas.iter().sum();
        let mut plate_area: BTreeMap<PlateId, f64> = BTreeMap::new();
        for c in 0..n {
            *plate_area.entry(seeding.assignment[c]).or_insert(0.0) += areas[c];
        }
        let cont_plate_area: f64 = seeding
            .continental_plate_ids
            .iter()
            .map(|p| plate_area.get(p).copied().unwrap_or(0.0))
            .sum();
        let mut in_body = vec![false; n];
        let mut bodies = Vec::new();
        for &pid in &seeding.continental_plate_ids {
            let k = pid.0 as usize;
            let seed_pt = *seeding
                .seed_points
                .get(k)
                .ok_or(KinematicError::UnknownPlate { plate: pid })?;
            let a_plate = plate_area.get(&pid).copied().unwrap_or(0.0);
            let target = area_fraction * total_area * a_plate / cont_plate_area;
            if target > a_plate * (1.0 + 1e-12) {
                return Err(KinematicError::Grid(format!(
                    "continental body for plate {} needs {target:.4e} m2 but the plate has {a_plate:.4e} m2 \
                     (declared area fraction {area_fraction} exceeds the continental plates' area); refused",
                    pid.0
                )));
            }
            // Plate cells by (dot with seed point descending = distance
            // ascending), ties to the lowest cell id.
            let mut cells: Vec<(f64, usize)> = (0..n)
                .filter(|&c| seeding.assignment[c] == pid)
                .map(|c| {
                    let p = cents[c];
                    (p[0] * seed_pt[0] + p[1] * seed_pt[1] + p[2] * seed_pt[2], c)
                })
                .collect();
            cells.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
            let mut acc = 0.0;
            let mut count = 0u64;
            for &(_, c) in &cells {
                // Include the next cell while doing so leaves the area at
                // least as close to the target as stopping would.
                if acc + 0.5 * areas[c] > target {
                    break;
                }
                in_body[c] = true;
                acc += areas[c];
                count += 1;
            }
            bodies.push(ContinentalBodyRecord {
                plate: pid,
                seed_point: seed_pt,
                target_area_m2: target,
                achieved_area_m2: acc,
                cells: count,
                crust_volume_m3: 0.0,
                tapered_cells: 0,
            });
        }

        // Coast edge midpoints of the global continental mask.
        let mut coast: Vec<[f64; 3]> = Vec::new();
        for c in 0..n {
            if !in_body[c] {
                continue;
            }
            let (i, j) = grid
                .ij(c as u64)
                .map_err(|e| KinematicError::Grid(e.to_string()))?;
            for (ti, tj, _) in grid
                .neighbors(i, j)
                .map_err(|e| KinematicError::Grid(e.to_string()))?
            {
                let nb = ti * grid.nlon() + tj;
                if in_body[nb] {
                    continue;
                }
                let raw = [
                    cents[c][0] + cents[nb][0],
                    cents[c][1] + cents[nb][1],
                    cents[c][2] + cents[nb][2],
                ];
                let norm = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
                if norm < 1e-12 {
                    return Err(KinematicError::Grid("degenerate coast midpoint".into()));
                }
                coast.push([raw[0] / norm, raw[1] / norm, raw[2] / norm]);
            }
        }
        let radius = grid.radius_m();
        let mut kind = Vec::with_capacity(n);
        let mut birth = Vec::with_capacity(n);
        let mut thick = Vec::with_capacity(n);
        let mut coast_distance = vec![f64::NAN; n];
        for c in 0..n {
            if in_body[c] {
                // Nearest coast midpoint by maximum dot product; a body
                // with no coast (whole planet) is unreachable here because
                // area_fraction < 1.
                let mut best = -2.0f64;
                for m in &coast {
                    let d = cents[c][0] * m[0] + cents[c][1] * m[1] + cents[c][2] * m[2];
                    if d > best {
                        best = d;
                    }
                }
                let dist = if coast.is_empty() {
                    f64::INFINITY
                } else {
                    best.clamp(-1.0, 1.0).acos() * radius
                };
                coast_distance[c] = dist;
                kind.push(CrustKindCodeV0::Continental);
                birth.push(-prior.initial_continental_age_my);
                thick.push(prior.rifted_margin_thickness_m(dist, rifted_margin_width_m));
            } else {
                kind.push(CrustKindCodeV0::Oceanic);
                birth.push(-prior.initial_oceanic_age_my);
                thick.push(prior.record.oceanic_crust_thickness_m);
            }
        }
        let interior = prior.record.continental_crust_thickness_m;
        let mut achieved = 0.0;
        let mut volume = 0.0;
        for b in &mut bodies {
            for c in 0..n {
                if in_body[c] && seeding.assignment[c] == b.plate {
                    b.crust_volume_m3 += thick[c] * areas[c];
                    if thick[c] < interior {
                        b.tapered_cells += 1;
                    }
                }
            }
            achieved += b.achieved_area_m2;
            volume += b.crust_volume_m3;
        }
        let state = Self {
            t_my: 0.0,
            plate: seeding.assignment.clone(),
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
        let report = ContinentalBodiesReport {
            declared_area_fraction: area_fraction,
            achieved_area_fraction: achieved / total_area,
            rifted_margin_width_m,
            interior_thickness_m: interior,
            edge_thickness_m: prior.record.oceanic_crust_thickness_m,
            continental_crust_volume_m3: volume,
            plateau_equivalent_volume_m3: achieved * interior,
            coast_midpoints: coast.len() as u64,
            coast_distance_m: coast_distance,
            bodies,
        };
        Ok((state, report))
    }

    pub fn cells(&self) -> usize {
        self.plate.len()
    }

    pub fn age_my(&self, cell: usize) -> f64 {
        self.t_my - self.birth_t_my[cell]
    }

    /// Total crust volume (m^3) over the sphere (bedrock only).
    pub fn crust_volume_m3(&self, areas: &[f64]) -> f64 {
        self.thickness_m.iter().zip(areas).map(|(h, a)| h * a).sum()
    }
}

/// Classify every cell that touches another plate. Returns the per-cell
/// class codes, the signed normal rate (m/My) and the footprint count per
/// (ordered pair, kind), for the boundary table.
pub fn classify_cells(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
) -> Result<
    (
        Vec<u8>,
        Vec<f64>,
        BTreeMap<(PlateId, PlateId, BoundaryKindV0), u64>,
    ),
    KinematicError,
> {
    let n = grid.cells();
    if state.cells() != n {
        return Err(KinematicError::SizeMismatch {
            what: "state",
            expected: n,
            found: state.cells(),
        });
    }
    let cents = centroids(grid)?;
    let radius = grid.radius_m();
    let poles: BTreeMap<PlateId, EulerPoleV0> = plates
        .iter()
        .filter_map(|p| p.euler_pole.map(|e| (p.plate_id, e)))
        .collect();
    // SEED-PERF (lane SP, 2026-10-13): per-cell in parallel, neighbour ids
    // without the unused haversines, reduced serially in ascending cell
    // order (bit for bit the serial loop, same first error).
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    let per_cell = |c: usize| -> Result<Option<(u8, f64, (PlateId, PlateId, BoundaryKindV0))>, KinematicError> {
        let a = state.plate[c];
        let pole_a = *poles
            .get(&a)
            .ok_or(KinematicError::UnknownPlate { plate: a })?;
        let mut best: Option<(f64, PlateId, bool)> = None; // (signed rate, other plate, both continental)
        let mut nbs = [0usize; 8];
        let nn = neighbor_ids(nlat, nlon, c, &mut nbs);
        for &nb in &nbs[..nn] {
            let b = state.plate[nb];
            if b == a {
                continue;
            }
            let pole_b = *poles
                .get(&b)
                .ok_or(KinematicError::UnknownPlate { plate: b })?;
            // Local normal: from this cell toward the neighbouring plate's
            // cell (boundary_normal projects ref_b - ref_a onto the tangent
            // plane at r; with ref_a = r that is the neighbour direction).
            let r = normal_rate_at(cents[c], cents[c], cents[nb], &pole_a, &pole_b)?;
            let both_cont = state.kind[c] == CrustKindCodeV0::Continental
                && state.kind[nb] == CrustKindCodeV0::Continental;
            let replace = match best {
                None => true,
                Some((br, bp, _)) => r.abs() > br.abs() || (r.abs() == br.abs() && b < bp),
            };
            if replace {
                best = Some((r, b, both_cont));
            }
        }
        let Some((r, b, both_cont)) = best else { return Ok(None) };
        let kind = classify(r, both_cont)?;
        let code = match kind {
            BoundaryKindV0::Ridge => BoundaryClassCodeV0::Ridge,
            BoundaryKindV0::Subduction => BoundaryClassCodeV0::Subduction,
            BoundaryKindV0::Collision => BoundaryClassCodeV0::Collision,
            BoundaryKindV0::Transform => BoundaryClassCodeV0::Transform,
        }
        .code();
        let key = if a <= b { (a, b, kind) } else { (b, a, kind) };
        Ok(Some((code, r * radius, key)))
    };
    let cells = par::map_indexed(n, per_cell);
    let mut class = vec![BoundaryClassCodeV0::Interior.code(); n];
    let mut rate = vec![0.0f64; n];
    let mut footprint: BTreeMap<(PlateId, PlateId, BoundaryKindV0), u64> = BTreeMap::new();
    for (c, r) in cells.into_iter().enumerate() {
        if let Some((code, r, key)) = r? {
            class[c] = code;
            rate[c] = r;
            *footprint.entry(key).or_insert(0) += 1;
        }
    }
    Ok((class, rate, footprint))
}

/// `classify_cells` under a chosen boundary law: `None` is the legacy
/// neighbour-direction classifier (bit for bit `classify_cells`);
/// `Some(p)` is `smoothed-normal-transform-v1` (SEED-TRANSFORMS, lane TF).
#[allow(clippy::type_complexity)]
pub fn classify_cells_law(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    law: Option<&crate::boundaries::TransformParams>,
) -> Result<
    (
        Vec<u8>,
        Vec<f64>,
        BTreeMap<(PlateId, PlateId, BoundaryKindV0), u64>,
    ),
    KinematicError,
> {
    match law {
        None => classify_cells(grid, state, plates),
        Some(p) => classify_cells_v2(grid, state, plates, p).map(|(c, r, f, _)| (c, r, f)),
    }
}

/// SEED-TRANSFORMS (lane TF) `smoothed-normal-transform-v1` classifier (see
/// `boundaries::LAW_TRANSFORM_V1`). Per boundary cell (plate a): for every
/// other plate b among its 8 neighbours (ascending id), the smoothed a->b
/// normal at the cell centre, falling back to the staircase direction
/// (mean direction to b's neighbour cells) when degenerate or pointing
/// away from b; the relative velocity (w_a - w_b) x r split into normal and
/// tangential parts; the pair with the fastest relative motion classifies
/// the cell (ties: lowest id). Both-continental = the cell and at least one
/// of b's neighbour cells continental. Also returns the count of fallback
/// normals.
#[allow(clippy::type_complexity)]
pub fn classify_cells_v2(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    p: &crate::boundaries::TransformParams,
) -> Result<
    (
        Vec<u8>,
        Vec<f64>,
        BTreeMap<(PlateId, PlateId, BoundaryKindV0), u64>,
        u64,
    ),
    KinematicError,
> {
    p.validate().map_err(KinematicError::Grid)?;
    let n = grid.cells();
    if state.cells() != n {
        return Err(KinematicError::SizeMismatch {
            what: "state",
            expected: n,
            found: state.cells(),
        });
    }
    let cents = centroids(grid)?;
    let areas = cell_areas(grid)?;
    let radius = grid.radius_m();
    let mut omega: BTreeMap<PlateId, [f64; 3]> = BTreeMap::new();
    for pr in plates {
        if let Some(e) = pr.euler_pole {
            omega.insert(pr.plate_id, crate::boundaries::angular_velocity(&e)?);
        }
    }
    let sigma = p.normal_smoothing_m / radius;
    let reach = crate::boundaries::TRANSFORM_WINDOW_SIGMAS * sigma;
    let (nlat, nlon) = (grid.nlat() as i64, grid.nlon() as i64);
    let dlat = std::f64::consts::PI / nlat as f64;
    let dlon = 2.0 * std::f64::consts::PI / nlon as f64;
    let row_half = (reach / dlat).ceil() as i64;
    // SEED-PERF (lane SP, 2026-10-13): per-cell work in parallel
    // (`par::map_indexed`), reduced serially in ascending cell order; the
    // neighbour ids without the unused haversines; the per-row window
    // half-widths hoisted; moments in a small vector (same per-plate
    // accumulation order); window cells provably beyond `reach` skipped by
    // squared chord before the exact `great_circle_rad` (margin as in
    // `MaterialSearch`: chord^2 >= (2 - 2 cos reach)(1 + 1e-9) + 1e-13
    // implies the computed angle exceeds reach). Bit for bit the serial
    // form.
    let dot3 = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let tangent = |v: [f64; 3], at: [f64; 3]| {
        let d = dot3(v, at);
        [v[0] - d * at[0], v[1] - d * at[1], v[2] - d * at[2]]
    };
    let norm3 = |v: [f64; 3]| dot3(v, v).sqrt();
    let col_half: Vec<i64> = (0..grid.nlat())
        .map(|ii| {
            let lat = grid
                .center_lat_deg(ii)
                .map_err(|e| KinematicError::Grid(e.to_string()))?;
            let cos_lat = lat.to_radians().cos().max(1.0e-9);
            Ok((reach / (dlon * cos_lat)).ceil() as i64)
        })
        .collect::<Result<_, KinematicError>>()?;
    let skip_chord2 = if reach < 3.0 {
        let h = (0.5 * reach).sin();
        4.0 * h * h * (1.0 + 1.0e-9) + 1.0e-13
    } else {
        f64::INFINITY
    };
    let (nlat_u, nlon_u) = (grid.nlat(), grid.nlon());
    type CellOut = (Option<(u8, f64, (PlateId, PlateId, BoundaryKindV0))>, u64);
    let per_cell = |c: usize| -> Result<CellOut, KinematicError> {
        let (i, j) = ((c / nlon_u) as i64, (c % nlon_u) as i64);
        let mut fallbacks = 0u64;
        let a = state.plate[c];
        // Other plates among the 8 neighbours: (plate, staircase sum, any
        // continental neighbour cell of that plate).
        let mut others: BTreeMap<PlateId, ([f64; 3], bool)> = BTreeMap::new();
        let mut nbs = [0usize; 8];
        let nn = neighbor_ids(nlat_u, nlon_u, c, &mut nbs);
        // PLATE-POLES `rigid-group-contacts-v1`: a plate with bit-identical
        // angular velocity is the same rigid body, not a boundary partner.
        let wa_rigid = if p.rigid_group_interior { omega.get(&a).copied() } else { None };
        for &nb in &nbs[..nn] {
            let b = state.plate[nb];
            if b == a {
                continue;
            }
            if let Some(wa) = wa_rigid {
                if omega.get(&b).map_or(false, |wb| wb.map(f64::to_bits) == wa.map(f64::to_bits)) {
                    continue;
                }
            }
            let t = tangent(cents[nb], cents[c]);
            let tn = norm3(t).max(1e-300);
            let e = others.entry(b).or_insert(([0.0; 3], false));
            for k in 0..3 {
                e.0[k] += t[k] / tn;
            }
            e.1 |= state.kind[nb] == CrustKindCodeV0::Continental;
        }
        if others.is_empty() {
            return Ok((None, 0));
        }
        let wa = *omega
            .get(&a)
            .ok_or(KinematicError::UnknownPlate { plate: a })?;
        // Indicator moments of a and of each other plate over the window.
        let mut mom: Vec<(PlateId, [f64; 3])> = Vec::with_capacity(others.len() + 1);
        mom.push((a, [0.0; 3]));
        for b in others.keys() {
            mom.push((*b, [0.0; 3]));
        }
        let cc = cents[c];
        for ii in (i - row_half).max(0)..=(i + row_half).min(nlat - 1) {
            let ch = col_half[ii as usize];
            let (lo, hi) = if 2 * ch + 1 >= nlon {
                (0, nlon - 1)
            } else {
                (j - ch, j + ch)
            };
            for jj in lo..=hi {
                let k = (ii * nlon + jj.rem_euclid(nlon)) as usize;
                let pk = state.plate[k];
                let Some(slot) = mom.iter_mut().position(|(p, _)| *p == pk) else { continue };
                let ck = cents[k];
                let (x, y, z) = (ck[0] - cc[0], ck[1] - cc[1], ck[2] - cc[2]);
                if x * x + y * y + z * z > skip_chord2 {
                    continue;
                }
                let d = great_circle_rad(cc, ck)?;
                if d > reach {
                    continue;
                }
                let w = (-0.5 * (d / sigma) * (d / sigma)).exp() * areas[k];
                let t = tangent(ck, cc);
                let m = &mut mom[slot].1;
                for q in 0..3 {
                    m[q] += w * t[q];
                }
            }
        }
        let ma = mom[0].1;
        let mut best: Option<(f64, PlateId, f64, f64, bool)> = None; // (speed, b, vn, vt, both cont)
        for (bi, (b, (stair, b_cont))) in others.iter().enumerate() {
            let mb = mom[bi + 1].1;
            let raw = [mb[0] - ma[0], mb[1] - ma[1], mb[2] - ma[2]];
            let rn = norm3(raw);
            let sn = norm3(*stair);
            let nrm = if rn > 1.0e-12 && dot3(raw, *stair) > 0.0 {
                [raw[0] / rn, raw[1] / rn, raw[2] / rn]
            } else if sn > 1.0e-12 {
                fallbacks += 1;
                [stair[0] / sn, stair[1] / sn, stair[2] / sn]
            } else {
                fallbacks += 1;
                continue;
            };
            let wb = *omega
                .get(b)
                .ok_or(KinematicError::UnknownPlate { plate: *b })?;
            let v = cross3([wa[0] - wb[0], wa[1] - wb[1], wa[2] - wb[2]], cents[c]);
            let vn = dot3(v, nrm);
            let vt = norm3([v[0] - vn * nrm[0], v[1] - vn * nrm[1], v[2] - vn * nrm[2]]);
            let speed = vn.hypot(vt);
            let both = state.kind[c] == CrustKindCodeV0::Continental && *b_cont;
            if best.map_or(true, |(bs, _, _, _, _)| speed > bs) {
                best = Some((speed, *b, vn, vt, both));
            }
        }
        let Some((_, b, vn, vt, both)) = best else {
            return Ok((None, fallbacks));
        };
        let kind = crate::boundaries::classify_oblique(vn, vt, both, p.transform_obliquity_deg)?;
        let code = match kind {
            BoundaryKindV0::Ridge => BoundaryClassCodeV0::Ridge,
            BoundaryKindV0::Subduction => BoundaryClassCodeV0::Subduction,
            BoundaryKindV0::Collision => BoundaryClassCodeV0::Collision,
            BoundaryKindV0::Transform => BoundaryClassCodeV0::Transform,
        }
        .code();
        let key = if a <= b { (a, b, kind) } else { (b, a, kind) };
        Ok((Some((code, vn * radius, key)), fallbacks))
    };
    let cells = par::map_indexed(n, per_cell);
    let mut class = vec![BoundaryClassCodeV0::Interior.code(); n];
    let mut rate = vec![0.0f64; n];
    let mut footprint: BTreeMap<(PlateId, PlateId, BoundaryKindV0), u64> = BTreeMap::new();
    let mut fallbacks = 0u64;
    for (c, r) in cells.into_iter().enumerate() {
        let (hit, fb) = r?;
        fallbacks += fb;
        if let Some((code, r, key)) = hit {
            class[c] = code;
            rate[c] = r;
            *footprint.entry(key).or_insert(0) += 1;
        }
    }
    Ok((class, rate, footprint, fallbacks))
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// SEED-TRANSFORMS (lane TF): the step-start boundary class governing a
/// cell: its own class if it is a boundary cell, else that of the nearest
/// boundary cell (great circle) within `TRANSFORM_CLASS_SEARCH_CELLS` rows
/// and columns (ties: lowest cell id); Interior when none.
fn local_class(grid: &GridSpec, cents: &[[f64; 3]], cls: &[u8], c: usize) -> Result<u8, KinematicError> {
    let interior = BoundaryClassCodeV0::Interior.code();
    if cls[c] != interior {
        return Ok(cls[c]);
    }
    let (nlat, nlon) = (grid.nlat() as i64, grid.nlon() as i64);
    let (ci, cj) = ((c / grid.nlon()) as i64, (c % grid.nlon()) as i64);
    let w = crate::boundaries::TRANSFORM_CLASS_SEARCH_CELLS;
    let mut best: Option<(f64, usize)> = None;
    for i in (ci - w).max(0)..=(ci + w).min(nlat - 1) {
        for dj in -w..=w {
            let k = (i * nlon + (cj + dj).rem_euclid(nlon)) as usize;
            if cls[k] == interior {
                continue;
            }
            let d = great_circle_rad(cents[c], cents[k])?;
            if best.map_or(true, |(bd, bk)| d < bd || (d == bd && k < bk)) {
                best = Some((d, k));
            }
        }
    }
    Ok(best.map_or(interior, |(_, k)| cls[k]))
}

/// Half-width of the source search window in ROWS (and in columns at the
/// equator). Phases stay within about half a cell of the centroid and
/// the back-traced point is within half a cell of its nearest centroid,
/// so 2 is generous in great-circle terms; the column half-width grows as
/// 1/cos(lat) per row so the window stays 2 cells wide in METRES and
/// covers whole rows at the poles (a 5 x 5 index window would miss a
/// cell's own material there and mint spurious ridges).
pub const SEARCH_HALF_WIDTH: i64 = 2;

/// The cell whose material position is nearest `q`, searched over the
/// metric window around `nearest_cell(q)`; ties break to the lowest cell id.
fn nearest_material(
    grid: &GridSpec,
    material_pos: &[[f64; 3]],
    q: [f64; 3],
) -> Result<usize, KinematicError> {
    let centre = nearest_cell(grid, q)?;
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
        let lat = grid
            .center_lat_deg(i as usize)
            .map_err(|e| KinematicError::Grid(e.to_string()))?;
        let cos_lat = lat.to_radians().cos().max(1.0e-9);
        let dj_half = ((SEARCH_HALF_WIDTH as f64 / cos_lat).ceil() as i64).min(nlon / 2);
        for dj in -dj_half..=dj_half {
            let j = (cj + dj).rem_euclid(nlon);
            let s = (i * nlon + j) as usize;
            let d = great_circle_rad(material_pos[s], q)?;
            let replace = match best {
                None => true,
                Some((bd, bs)) => d < bd || (d == bd && s < bs),
            };
            if replace {
                best = Some((d, s));
            }
        }
    }
    Ok(best.expect("window is never empty").1)
}

/// SEED-PERF (lane SP, 2026-10-13): `nearest_material` with the plate
/// prefilter of the carry, bit for bit.
///
/// The legacy carry asks `nearest_material` for every (cell, live plate)
/// and keeps the answer only if the winning source belongs to that plate.
/// This search returns that kept answer (`Some(source)`) or `None`:
/// 1. if no cell of the plate is in the window, the winner cannot belong to
///    it: `None` without any distance;
/// 2. otherwise candidates are ranked by squared chord |m - q|^2, which is
///    monotone in the great-circle angle up to rounding; every candidate
///    whose chord is within a proven rounding margin of the smallest one
///    (all near-ties) is then ranked by the exact `great_circle_rad` with
///    the legacy rule (smaller angle, then lower cell id). Since the legacy
///    winner always lies inside that margin, the result is the legacy one.
///
/// Margin: with |m|^2 = 1 + a, |q|^2 = 1 + b, chord^2 = (|m| - |q|)^2 +
/// |m||q| (2 - 2 cos t), so the chord carries a relative perturbation of
/// at most `alpha` = max |a| + |b| plus a few ulps of rounding, while the
/// computed atan2 angle is within about 2e-15 rad of the true one, and
/// 2 - 2 cos moves by at most 2 per radian. `rel` (8 alpha + 1e-12) and
/// `abs` (1e-13) bound both with a wide factor. If any material position
/// is non-finite or `alpha` exceeds 1e-9 the exact legacy scan is used for
/// the whole step (identical results and identical errors).
pub(crate) struct MaterialSearch<'a> {
    grid: &'a GridSpec,
    material_pos: &'a [[f64; 3]],
    plate: &'a [PlateId],
    /// Per row: the column half-width `nearest_material` derives from it.
    dj_half: Vec<i64>,
    /// None: use the exact legacy scan.
    tol: Option<(f64, f64)>,
}

impl<'a> MaterialSearch<'a> {
    pub(crate) fn new(
        grid: &'a GridSpec,
        material_pos: &'a [[f64; 3]],
        plate: &'a [PlateId],
    ) -> Result<Self, KinematicError> {
        let nlon = grid.nlon() as i64;
        let mut dj_half = Vec::with_capacity(grid.nlat());
        for i in 0..grid.nlat() {
            let lat = grid
                .center_lat_deg(i)
                .map_err(|e| KinematicError::Grid(e.to_string()))?;
            let cos_lat = lat.to_radians().cos().max(1.0e-9);
            dj_half.push(((SEARCH_HALF_WIDTH as f64 / cos_lat).ceil() as i64).min(nlon / 2));
        }
        let mut alpha = 0.0f64;
        let mut finite = true;
        for m in material_pos {
            if !(m[0].is_finite() && m[1].is_finite() && m[2].is_finite()) {
                finite = false;
                break;
            }
            alpha = alpha.max((m[0] * m[0] + m[1] * m[1] + m[2] * m[2] - 1.0).abs());
        }
        let tol = if finite && alpha <= 1.0e-9 {
            Some((8.0 * (alpha + 1.0e-15) + 1.0e-12, 1.0e-13))
        } else {
            None
        };
        Ok(Self {
            grid,
            material_pos,
            plate,
            dj_half,
            tol,
        })
    }

    /// The legacy carry's kept claim of plate `id` at back-traced point `q`.
    pub(crate) fn claim(&self, q: [f64; 3], id: PlateId) -> Result<Option<usize>, KinematicError> {
        Ok(self
            .winner_if(q, |p| p == id)?
            .filter(|&s| self.plate[s] == id))
    }

    /// The legacy `nearest_material(q)` winner, or None when no cell of
    /// the window satisfies `wanted` (then the winner's plate does not
    /// either). With `wanted` = one plate this is `claim`; with a set of
    /// plates that share one rotation (one back-traced `q`) it serves all
    /// of them with one search.
    pub(crate) fn winner_if(
        &self,
        q: [f64; 3],
        wanted: impl Fn(PlateId) -> bool,
    ) -> Result<Option<usize>, KinematicError> {
        let Some((rel, abs)) = self.tol else {
            return Ok(Some(nearest_material(self.grid, self.material_pos, q)?));
        };
        let grid = self.grid;
        let centre = nearest_cell(grid, q)?;
        // q is finite here (nearest_cell refuses non-finite input), so the
        // legacy great_circle_rad calls could not have failed either.
        let nlat = grid.nlat() as i64;
        let nlon = grid.nlon() as i64;
        let ci = (centre / grid.nlon()) as i64;
        let cj = (centre % grid.nlon()) as i64;
        let (r0, r1) = ((ci - SEARCH_HALF_WIDTH).max(0), (ci + SEARCH_HALF_WIDTH).min(nlat - 1));
        let mp = self.material_pos;
        let qq = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] - 1.0).abs();
        // Pass 0: is the plate in the window at all (ids only).
        let mut present = false;
        'rows: for i in r0..=r1 {
            let h = self.dj_half[i as usize];
            let row = i * nlon;
            for dj in -h..=h {
                if wanted(self.plate[(row + (cj + dj).rem_euclid(nlon)) as usize]) {
                    present = true;
                    break 'rows;
                }
            }
        }
        if !present {
            return Ok(None);
        }
        // Pass 1: the smallest squared chord.
        let mut best_c = f64::INFINITY;
        for i in r0..=r1 {
            let h = self.dj_half[i as usize];
            let row = i * nlon;
            for dj in -h..=h {
                let m = mp[(row + (cj + dj).rem_euclid(nlon)) as usize];
                let (x, y, z) = (m[0] - q[0], m[1] - q[1], m[2] - q[2]);
                let c2 = x * x + y * y + z * z;
                if c2 < best_c {
                    best_c = c2;
                }
            }
        }
        // Pass 2: exact legacy ranking over the near-ties only.
        let bound = best_c * (1.0 + rel + 4.0 * qq) + abs;
        let mut best: Option<(f64, usize)> = None;
        for i in r0..=r1 {
            let h = self.dj_half[i as usize];
            let row = i * nlon;
            for dj in -h..=h {
                let s = (row + (cj + dj).rem_euclid(nlon)) as usize;
                let m = mp[s];
                let (x, y, z) = (m[0] - q[0], m[1] - q[1], m[2] - q[2]);
                if x * x + y * y + z * z > bound {
                    continue;
                }
                let d = great_circle_rad(m, q)?;
                let replace = match best {
                    None => true,
                    Some((bd, bs)) => d < bd || (d == bd && s < bs),
                };
                if replace {
                    best = Some((d, s));
                }
            }
        }
        Ok(Some(best.expect("window is never empty").1))
    }
}

/// SEED-PERF (lane SP): the angle within which a plate's cell must lie for
/// the plate to claim a cell this step (see `plate_reach_masks`): the
/// largest rotation angle plus 5 (dlat + dlon) (4 for the window geometry,
/// 1 of slack).
fn claim_reach_radius(grid: &GridSpec, live: &[(PlateId, EulerPoleV0)], dt_my: f64) -> f64 {
    let theta = live
        .iter()
        .map(|(_, e)| (e.rate_deg_per_my * dt_my).abs().to_radians())
        .fold(0.0f64, f64::max);
    let cell = std::f64::consts::PI / grid.nlat() as f64 + 2.0 * std::f64::consts::PI / grid.nlon() as f64;
    theta * (1.0 + 1.0e-9) + 1.0e-12 + 5.0 * cell
}

/// SEED-PERF (lane SP, 2026-10-13): for every cell, a bit mask of the live
/// plates (bit k = `live_ids[k]`) that own some cell whose centroid lies
/// within great-circle angle `radius` of the cell's centroid (a superset:
/// index-space dilation, conservative in rows and columns). Used to skip
/// (cell, plate) claims that cannot exist: a claim's source lies in the
/// search window around `nearest_cell(q)`, every window centroid is within
/// 2 dlat + 3 dlon of that cell's centroid, which is within dlat + dlon of
/// q, and q is within the plate's rotation angle of the cell's centroid
/// (callers add one more dlat + dlon of slack). None when not applicable
/// (more than 64 live plates, or a radius of a quarter circle or more).
fn plate_reach_masks(
    grid: &GridSpec,
    plate: &[PlateId],
    live_ids: &[PlateId],
    radius: f64,
) -> Option<Vec<u64>> {
    if live_ids.len() > 64
        || !(radius < std::f64::consts::FRAC_PI_2)
        || live_ids.windows(2).any(|w| w[0] >= w[1])
    {
        return None;
    }
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    let dlat = std::f64::consts::PI / nlat as f64;
    let dlon = 2.0 * std::f64::consts::PI / nlon as f64;
    let wr = (radius / dlat).floor() as usize + 1;
    let cos: Vec<f64> = (0..nlat)
        .map(|i| (-std::f64::consts::FRAC_PI_2 + (i as f64 + 0.5) * dlat).cos().max(0.0))
        .collect();
    let bit = |p: PlateId| -> u64 {
        match live_ids.binary_search(&p) {
            Ok(k) => 1u64 << k,
            Err(_) => 0,
        }
    };
    // Horizontal pass per row: OR over +-hw columns (circular), hw wide
    // enough for every target row within +-wr rows.
    let rows: Vec<Vec<u64>> = par::map_heavy(nlat, |i| {
        let row: Vec<u64> = (0..nlon).map(|j| bit(plate[i * nlon + j])).collect();
        let lo = i.saturating_sub(wr);
        let hi = (i + wr).min(nlat - 1);
        let cmin = cos[lo..=hi].iter().cloned().fold(f64::INFINITY, f64::min);
        let denom = 2.0 * (cos[i] * cmin).sqrt();
        let whole = |row: &[u64]| vec![row.iter().fold(0u64, |a, b| a | b); nlon];
        if !(radius < denom) {
            return whole(&row);
        }
        let dl = 2.0 * (radius / denom).asin();
        let hw = (dl / dlon).floor() as usize + 1;
        let w = 2 * hw + 1;
        if w >= nlon {
            return whole(&row);
        }
        // van Herk / Gil-Werman running OR over the circularly extended row.
        let len = nlon + 2 * hw;
        let ext: Vec<u64> = (0..len).map(|t| row[(t + nlon - hw % nlon) % nlon]).collect();
        let mut pre = vec![0u64; len];
        let mut suf = vec![0u64; len];
        for t in 0..len {
            pre[t] = if t % w == 0 { ext[t] } else { pre[t - 1] | ext[t] };
        }
        for t in (0..len).rev() {
            suf[t] = if t % w == w - 1 || t == len - 1 { ext[t] } else { suf[t + 1] | ext[t] };
        }
        (0..nlon).map(|j| suf[j] | pre[j + w - 1]).collect()
    });
    // Vertical pass: OR over +-wr rows.
    Some(par::map_indexed(nlat * nlon, |c| {
        let (i, j) = (c / nlon, c % nlon);
        let mut m = 0u64;
        for r in &rows[i.saturating_sub(wr)..=(i + wr).min(nlat - 1)] {
            m |= r[j];
        }
        m
    }))
}

/// Oriented input-state interface edges. Geometric normals do not come from
/// material backtracking or the velocity being projected. Lowest edge ID breaks ties.
type InterfaceSamples = BTreeMap<(PlateId, PlateId), Vec<(usize, usize)>>;
fn interface_samples(
    grid: &GridSpec,
    state: &CrustState,
    cents: &[[f64; 3]],
) -> Result<InterfaceSamples, KinematicError> {
    let mut pairs: InterfaceSamples = BTreeMap::new();
    // SEED-PERF (lane SP): neighbour ids without the unused haversines.
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    let mut nbs = [0usize; 8];
    for c in 0..cents.len() {
        grid.ij(c as u64)
            .map_err(|e| KinematicError::Grid(e.to_string()))?;
        let nn = neighbor_ids(nlat, nlon, c, &mut nbs);
        for &nb in &nbs[..nn] {
            if state.plate[c] < state.plate[nb] {
                pairs
                    .entry((state.plate[c], state.plate[nb]))
                    .or_default()
                    .push((c, nb));
            }
        }
    }
    for edges in pairs.values_mut() {
        edges.sort_unstable();
    }
    Ok(pairs)
}
/// Signed normal convergence rate (m/My) of an arc overlap between the
/// winner `a` (source cell `src_a`) and the loser `b` (source cell `src_b`).
///
/// 1. STEP-START INTERFACE (D5 rule, unchanged): the nearest step-start
///    interface edge of the pair, if it lies within the swept locality
///    bound (both plates' angular travel this step plus three cell
///    diagonals for raster phase).
/// 2. NEW CONTACT (SEED-ROBUST): otherwise the two plates were not in
///    contact here when the step began, yet both plates' material reached
///    this cell, so the contact formed during the step (a separating sliver
///    of a third plate was consumed or vacated, or a detached fragment
///    reached the other plate). The contact normal is the direction across
///    the material that met: from the winner's step-start source cell to
///    the loser's, sampled at their midpoint. The segment must itself be
///    within the same locality bound; otherwise a typed refusal. No remote
///    boundary is ever used.
#[allow(clippy::too_many_arguments)]
fn pair_normal_rate(
    grid: &GridSpec,
    pairs: &InterfaceSamples,
    cents: &[[f64; 3]],
    cell: usize,
    a: PlateId,
    b: PlateId,
    src_a: usize,
    src_b: usize,
    pa: &EulerPoleV0,
    pb: &EulerPoleV0,
    dt: f64,
) -> Result<
    (
        f64,
        [pf_contracts::spatial::ElementId; 2],
        ArcContactGeometry,
    ),
    KinematicError,
> {
    let key = if a < b { (a, b) } else { (b, a) };
    // Locality is bounded by both plates' swept angular travel plus three
    // source-cell diagonals for raster phase. No remote boundary fallback.
    let cell_angle = (std::f64::consts::PI / grid.nlat() as f64)
        .hypot(2. * std::f64::consts::PI / grid.nlon() as f64);
    let max_distance =
        (pa.rate_deg_per_my.abs() + pb.rate_deg_per_my.abs()).to_radians() * dt + 3. * cell_angle;
    let mut nearest_interface_rad = None;
    if let Some(edges) = pairs.get(&key) {
        // SEED-PERF (lane SP): the nearest edge midpoint by squared chord
        // first (degenerate midpoints refused in edge order, as before),
        // then the exact `great_circle_rad` over the near-ties only, first
        // strict minimum in edge order: the legacy pick bit for bit (margin
        // as in `MaterialSearch`; midpoints and centroids are unit to a few
        // ulps).
        let mid = |i: usize, j: usize| -> Result<[f64; 3], KinematicError> {
            let raw = [
                cents[i][0] + cents[j][0],
                cents[i][1] + cents[j][1],
                cents[i][2] + cents[j][2],
            ];
            let norm = (raw.iter().map(|x| x * x).sum::<f64>()).sqrt();
            if norm < 1e-9 {
                return Err(KinematicError::Grid("degenerate interface midpoint".into()));
            }
            Ok([raw[0] / norm, raw[1] / norm, raw[2] / norm])
        };
        let q = cents[cell];
        let chord2 = |p: [f64; 3]| {
            let (x, y, z) = (p[0] - q[0], p[1] - q[1], p[2] - q[2]);
            x * x + y * y + z * z
        };
        let mut best_c = f64::INFINITY;
        for &(i, j) in edges {
            let c2 = chord2(mid(i, j)?);
            if c2 < best_c {
                best_c = c2;
            }
        }
        let bound = best_c * (1.0 + 1.0e-12) + 1.0e-13;
        let mut best: Option<(f64, usize, usize, [f64; 3])> = None;
        for &(i, j) in edges {
            let p = mid(i, j)?;
            if chord2(p) > bound {
                continue;
            }
            let d = great_circle_rad(cents[cell], p)?;
            if best.as_ref().map(|v| d < v.0).unwrap_or(true) {
                best = Some((d, i, j, p));
            }
        }
        let (distance, i, j, p) =
            best.ok_or_else(|| KinematicError::Grid("empty pair interface".into()))?;
        // `!(d > bound)` keeps the pre-SEED-ROBUST acceptance set exactly.
        if !(distance > max_distance) {
            let (ra, rb) = if a < b {
                (cents[i], cents[j])
            } else {
                (cents[j], cents[i])
            };
            return Ok((
                normal_rate_at(p, ra, rb, pa, pb)? * grid.radius_m(),
                [
                    pf_contracts::spatial::ElementId(i as u64),
                    pf_contracts::spatial::ElementId(j as u64),
                ],
                ArcContactGeometry::StepStartInterface,
            ));
        }
        nearest_interface_rad = Some(distance);
    }
    // New contact formed during the step: sample across the met material.
    let (ra, rb) = (cents[src_a], cents[src_b]);
    let separation = great_circle_rad(ra, rb)?;
    let refuse = |why| KinematicError::ArcContact {
        cell,
        winner: a,
        loser: b,
        why,
    };
    if !(separation <= max_distance) {
        return Err(refuse(ArcContactRefusal::NonlocalContact {
            nearest_interface_rad,
            source_separation_rad: separation,
            bound_rad: max_distance,
        }));
    }
    let raw = [ra[0] + rb[0], ra[1] + rb[1], ra[2] + rb[2]];
    let norm = (raw.iter().map(|x| x * x).sum::<f64>()).sqrt();
    if !(norm >= 1e-9) || src_a == src_b {
        return Err(refuse(ArcContactRefusal::DegenerateContact {
            source_separation_rad: separation,
        }));
    }
    let p = [raw[0] / norm, raw[1] / norm, raw[2] / norm];
    if crate::boundaries::boundary_normal(p, ra, rb).is_err() {
        return Err(refuse(ArcContactRefusal::DegenerateContact {
            source_separation_rad: separation,
        }));
    }
    let rate = normal_rate_at(p, ra, rb, pa, pb)? * grid.radius_m();
    let (lo, hi) = if a < b {
        (src_a, src_b)
    } else {
        (src_b, src_a)
    };
    Ok((
        rate,
        [
            pf_contracts::spatial::ElementId(lo as u64),
            pf_contracts::spatial::ElementId(hi as u64),
        ],
        ArcContactGeometry::NewContactBridge,
    ))
}

/// Receive volume without rebasing the remainder into any cell's thickness
/// units. An already over-cap recipient has no capacity; incoming material
/// never erodes its pre-existing column as a side effect.
fn receive_wedge_volume(
    prior: &UpliftPrior,
    height: f64,
    volume: f64,
    area: f64,
) -> (f64, f64, f64) {
    let received = volume.min((prior.record.max_crust_thickness_m - height).max(0.) * area);
    (height + received / area, volume - received, received)
}

/// Advance the crust field by `dt_my` under rigid plate rotation with
/// the declared gap/overlap rules. Returns the new state and the step's
/// crust volume ledger. `plates` must carry a pole for every plate id
/// present in `state`.
pub fn advance(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>), KinematicError> {
    let (out, ledger, thickening, _) =
        advance_impl(
            grid,
            state,
            plates,
            prior,
            dt_my,
            None,
            CarryLaw::default(),
            crate::orogen::CollisionSedimentLaw::default(),
            None,
            None,
            RidgeAccretion::PreviousOccupantV0,
            None,
        None,
        )?;
    Ok((out, ledger, thickening))
}

/// Source of each destination cell's carried state for `advance_with`:
/// the winner's (or sole claimant's) step-start cell; `u32::MAX` at ridge
/// births. Lets a caller carry its own per-column fields (SEED-OROGEN
/// subgrid relief) with the crust.
pub const NO_SOURCE: u32 = u32::MAX;

/// `advance` with the tectonic thickening law selected: `None` is the
/// legacy `kinematic-overlap-thickening-v0` (bit for bit the same as
/// `advance`); `Some(p)` is `rate-limited-belt-v1` (`orogen.rs`). Also
/// returns the per-cell source map (`NO_SOURCE` at births).
#[allow(clippy::type_complexity)]
pub fn advance_with(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    advance_impl(
        grid,
        state,
        plates,
        prior,
        dt_my,
        belt,
        CarryLaw::default(),
        crate::orogen::CollisionSedimentLaw::default(),
        None,
        None,
        RidgeAccretion::PreviousOccupantV0,
        None,
        None,
    )
}

/// `advance_with` with the plate-carry law selected (SEED-CARRY-CC):
/// `CarryLaw::SemiLagrangianNearestMaterialV0` is `advance_with` bit for
/// bit; `CarryLaw::FrontRescueV1` is `front-rescue-carry-v1`
/// (`carry.rs`) and requires `belt` (refused otherwise).
#[allow(clippy::type_complexity)]
pub fn advance_with_carry(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    carry: CarryLaw,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    advance_impl(
        grid,
        state,
        plates,
        prior,
        dt_my,
        belt,
        carry,
        crate::orogen::CollisionSedimentLaw::default(),
        None,
        None,
        RidgeAccretion::PreviousOccupantV0,
        None,
        None,
    )
}

/// `advance_with_carry` with the collision-sediment law selected
/// (SEDIMENT-TOWERS, lane ST): `CollisionSedimentLaw::SurfaceLayerV0` is
/// `advance_with_carry` bit for bit; `AccretedToBeltCrustV1` requires
/// `belt` (refused otherwise).
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn advance_with_laws(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    carry: CarryLaw,
    collision_sediment: crate::orogen::CollisionSedimentLaw,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    advance_impl(grid, state, plates, prior, dt_my, belt, carry, collision_sediment, None, None, RidgeAccretion::PreviousOccupantV0, None, None)
}

/// `advance_with_laws` plus the SEED-EXTRAS (lane EX) passes: plume heads
/// and tails (`fixed-plume-chains-v1`) and rift births (`stretched-rift-v1`).
/// With both `None` it is `advance_with_laws` bit for bit. Requires `belt`
/// (refused otherwise; the run state that carries their markers is the belt
/// law's).
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn advance_with_extras(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    carry: CarryLaw,
    collision_sediment: crate::orogen::CollisionSedimentLaw,
    extras: &ExtrasLaws<'_>,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    if (extras.plumes.is_some() || extras.rift.is_some()) && belt.is_none() {
        return Err(KinematicError::Grid(
            "fixed-plume-chains-v1 / stretched-rift-v1 require tectonic_thickening rate-limited-belt-v1; refused"
                .into(),
        ));
    }
    advance_impl(
        grid,
        state,
        plates,
        prior,
        dt_my,
        belt,
        carry,
        collision_sediment,
        Some(extras),
        None,
        RidgeAccretion::PreviousOccupantV0,
        None,
        None,
    )
}

/// SEED-CARRY-CONSERVE (lane CV): `advance_with_extras` under the plate
/// carry `continental-parcel-carry-v1` (`conserve.rs`). `layer` holds the
/// continental parcels (moved, consumed and projected here; the caller
/// writes the step's changes back with `conserve::write_back` after its
/// surface stage). Requires `belt`. The step's parcel bookkeeping and the
/// projection are returned in `AdvanceLedger::parcels` /
/// `AdvanceLedger::parcel_projection`.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn advance_with_parcels(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    collision_sediment: crate::orogen::CollisionSedimentLaw,
    extras: &ExtrasLaws<'_>,
    layer: &mut crate::conserve::ParcelLayer,
    params: crate::conserve::ParcelParams,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    advance_with_parcels_ridge(
        grid,
        state,
        plates,
        prior,
        dt_my,
        belt,
        collision_sediment,
        extras,
        layer,
        params,
        RidgeAccretion::PreviousOccupantV0,
    )
}

/// SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13): `advance_with_parcels` with
/// the ridge accretion rule chosen (`PreviousOccupantV0` = bit for bit
/// `advance_with_parcels`).
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn advance_with_parcels_ridge(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    collision_sediment: crate::orogen::CollisionSedimentLaw,
    extras: &ExtrasLaws<'_>,
    layer: &mut crate::conserve::ParcelLayer,
    params: crate::conserve::ParcelParams,
    ridge: RidgeAccretion,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    if belt.is_none() {
        return Err(KinematicError::Grid(
            "continental-parcel-carry-v1 requires tectonic_thickening rate-limited-belt-v1 (its collision consumption is placed on the belt); refused"
                .into(),
        ));
    }
    params.validate().map_err(KinematicError::Grid)?;
    advance_impl(
        grid,
        state,
        plates,
        prior,
        dt_my,
        belt,
        CarryLaw::ContinentalParcelsV1,
        collision_sediment,
        Some(extras),
        Some((layer, params)),
        ridge,
        None,
        None,
    )
}

/// SEED-WORLD-COMPOSE (lane CP, 2026-10-13): `advance_with_parcels_ridge`
/// with `plateau-cap-spill-v2` run inside the belt step (`None` = bit for
/// bit `advance_with_parcels_ridge`).
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn advance_with_parcels_plateau(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    collision_sediment: crate::orogen::CollisionSedimentLaw,
    extras: &ExtrasLaws<'_>,
    layer: &mut crate::conserve::ParcelLayer,
    params: crate::conserve::ParcelParams,
    ridge: RidgeAccretion,
    plateau: Option<&crate::ranges::PlateauInBeltParams>,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    if belt.is_none() {
        return Err(KinematicError::Grid(
            "continental-parcel-carry-v1 requires tectonic_thickening rate-limited-belt-v1 (its collision consumption is placed on the belt); refused"
                .into(),
        ));
    }
    params.validate().map_err(KinematicError::Grid)?;
    advance_impl(
        grid,
        state,
        plates,
        prior,
        dt_my,
        belt,
        CarryLaw::ContinentalParcelsV1,
        collision_sediment,
        Some(extras),
        Some((layer, params)),
        ridge,
        plateau,
        None,
    )
}

/// SEED-TRANSFORMS (lane TF, 2026-10-13): `advance_with_parcels_plateau`
/// with the boundary law `smoothed-normal-transform-v1` (`None` = bit for
/// bit `advance_with_parcels_plateau`): step-start classification with the
/// smoothed normal and a transform class, gap fill from the same side on
/// transform and convergent boundaries, no subduction of overlaps on
/// transform and divergent boundaries, and the step-end classification by
/// the same law.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn advance_with_parcels_transform(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    collision_sediment: crate::orogen::CollisionSedimentLaw,
    extras: &ExtrasLaws<'_>,
    layer: &mut crate::conserve::ParcelLayer,
    params: crate::conserve::ParcelParams,
    ridge: RidgeAccretion,
    plateau: Option<&crate::ranges::PlateauInBeltParams>,
    transform: Option<&crate::boundaries::TransformParams>,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    if belt.is_none() {
        return Err(KinematicError::Grid(
            "continental-parcel-carry-v1 requires tectonic_thickening rate-limited-belt-v1 (its collision consumption is placed on the belt); refused"
                .into(),
        ));
    }
    params.validate().map_err(KinematicError::Grid)?;
    advance_impl(
        grid,
        state,
        plates,
        prior,
        dt_my,
        belt,
        CarryLaw::ContinentalParcelsV1,
        collision_sediment,
        Some(extras),
        Some((layer, params)),
        ridge,
        plateau,
        transform,
    )
}

#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn advance_impl(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    carry: CarryLaw,
    collision_sediment: crate::orogen::CollisionSedimentLaw,
    extras: Option<&ExtrasLaws<'_>>,
    parcels: Option<(&mut crate::conserve::ParcelLayer, crate::conserve::ParcelParams)>,
    ridge: RidgeAccretion,
    plateau: Option<&crate::ranges::PlateauInBeltParams>,
    transform: Option<&crate::boundaries::TransformParams>,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    advance_with_formation_age_v1(grid, state, plates, prior, dt_my, belt, carry,
        collision_sediment, extras, parcels, ridge, plateau, transform, None)
}

/// Optional native formation-age sidecar. `None` is the existing solver path.
/// Uses the same law arguments as the native wrappers; parcel carry still
/// requires its layer and belt. All times span one actual call, never inferred
/// substeps. Accepted gap births are recorded only after fill/projection/rift.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn advance_with_formation_age_v1(
    grid: &GridSpec,
    state: &CrustState,
    plates: &[PlateRecordV0],
    prior: &UpliftPrior,
    dt_my: f64,
    belt: Option<&crate::orogen::BeltParams>,
    carry: CarryLaw,
    collision_sediment: crate::orogen::CollisionSedimentLaw,
    extras: Option<&ExtrasLaws<'_>>,
    mut parcels: Option<(&mut crate::conserve::ParcelLayer, crate::conserve::ParcelParams)>,
    ridge: RidgeAccretion,
    plateau: Option<&crate::ranges::PlateauInBeltParams>,
    transform: Option<&crate::boundaries::TransformParams>,
    recorder: Option<&mut crate::formation_age::FormationAgeRecorderV1>,
) -> Result<(CrustState, AdvanceLedger, Vec<f64>, Vec<u32>), KinematicError> {
    if let Some(r) = recorder.as_ref() {
        r.preflight(grid, state.t_my).map_err(KinematicError::Grid)?;
        crate::formation_age::TimeInterval::new(state.t_my, state.t_my + dt_my)
            .map_err(KinematicError::Grid)?;
        if !(prior.record.oceanic_crust_thickness_m.is_finite()
            && prior.record.oceanic_crust_thickness_m > 0.0)
        {
            return Err(KinematicError::Grid("formation-age: invalid oceanic birth thickness".into()));
        }
        if let Some((_, params)) = parcels.as_ref() {
            params.validate().map_err(KinematicError::Grid)?;
            if belt.is_none() {
                return Err(KinematicError::Grid("continental-parcel-carry-v1 requires a belt".into()));
            }
        }
    }
    let n = grid.cells();
    // SEED-CARRY-CONSERVE: the parcel carry runs only with its layer.
    let parcel_on = parcels.is_some();
    if parcel_on != (carry == CarryLaw::ContinentalParcelsV1) {
        return Err(KinematicError::Grid(
            "continental-parcel-carry-v1 runs only through advance_with_parcels (with its parcel layer); refused"
                .into(),
        ));
    }
    if let Some(p) = belt {
        p.validate().map_err(KinematicError::Grid)?;
    }
    let accrete = collision_sediment == crate::orogen::CollisionSedimentLaw::AccretedToBeltCrustV1;
    if accrete && belt.is_none() {
        return Err(KinematicError::Grid(
            "accreted-to-belt-crust-v1 requires tectonic_thickening rate-limited-belt-v1 (the accreted sediment joins the belt's crust supply); refused"
                .into(),
        ));
    }
    let contact = carry == CarryLaw::FrontRescueV1;
    if contact && belt.is_none() {
        return Err(KinematicError::Grid(
            "front-rescue-carry-v1 requires tectonic_thickening rate-limited-belt-v1 (its collision supply is placed on the belt); refused"
                .into(),
        ));
    }
    if state.cells() != n {
        return Err(KinematicError::SizeMismatch {
            what: "state",
            expected: n,
            found: state.cells(),
        });
    }
    if !dt_my.is_finite() || dt_my <= 0.0 {
        return Err(KinematicError::NonFinite {
            what: "dt_my",
            value: dt_my,
            index: 0,
        });
    }
    let cents = centroids(grid)?;
    let areas = cell_areas(grid)?;
    let t_new = state.t_my + dt_my;
    let mut ledger = AdvanceLedger {
        volume_before_m3: state.crust_volume_m3(&areas),
        sediment_before_m3: state
            .sediment_m
            .iter()
            .zip(&areas)
            .map(|(h, a)| h * a)
            .sum(),
        ..Default::default()
    };

    // Live plates in ascending id order (deterministic claim order).
    let mut live: Vec<(PlateId, EulerPoleV0)> = plates
        .iter()
        .filter_map(|p| p.euler_pole.map(|e| (p.plate_id, e)))
        .collect();
    live.sort_by_key(|(id, _)| *id);
    // Back-rotation poles (negative rate = inverse rotation).
    let back: Vec<(PlateId, EulerPoleV0)> = live
        .iter()
        .map(|(id, e)| {
            (
                *id,
                EulerPoleV0 {
                    rate_deg_per_my: -e.rate_deg_per_my,
                    ..*e
                },
            )
        })
        .collect();
    let pole_map: BTreeMap<PlateId, EulerPoleV0> = live.iter().cloned().collect();
    let interfaces = interface_samples(grid, state, &cents)?;
    // SEED-TRANSFORMS (lane TF): step-start classes under the step's poles.
    let cls0: Option<Vec<u8>> = match transform {
        // PLATE-POLES `rigid-group-contacts-v1`: the remap class keeps welded
        // contacts as same-body remap zones (legacy rule), see boundaries.rs.
        Some(tp) => Some(classify_cells_v2(grid, state, plates, &crate::boundaries::TransformParams { rigid_group_interior: false, ..*tp })?.0),
        None => None,
    };

    let mut out = CrustState {
        t_my: t_new,
        plate: vec![PlateId(0); n],
        kind: vec![CrustKindCodeV0::Oceanic; n],
        birth_t_my: vec![0.0; n],
        thickness_m: vec![0.0; n],
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: cents.clone(),
    };
    // Collisional excess to pass outward: (cell, excess_m, direction unit
    // tangent). Processed after the main pass in ascending cell order.
    let mut passes: Vec<(usize, f64, [f64; 3])> = Vec::new();
    // Tectonic thickening received by each cell this step (m): arc
    // addition plus collisional stacking, including passed excess.
    let mut thickening = vec![0.0f64; n];
    // SEED-OROGEN belt events (rate-limited-belt-v1 only): collision
    // (cell, crust m3, sediment m3, mean colliding thickness m) and arc
    // (cell, volume m3, overriding plate), placed after the main pass.
    let mut collisions: Vec<(usize, f64, f64, f64)> = Vec::new();
    let mut arcs: Vec<(usize, f64, PlateId)> = Vec::new();
    let mut source = vec![NO_SOURCE; n];
    // Candidates are observations only; final source/kind decides acceptance.
    let mut formation_candidates = recorder.as_ref().map(|_| vec![None; n]);
    if belt.is_some() {
        // Under front-rescue-carry-v1 the edges are only counted (contact
        // length diagnostic); the supply is the consumed loser volume.
        kinematic_collision_supply(
            grid,
            state,
            &pole_map,
            &cents,
            dt_my,
            !contact && !parcel_on,
            &mut collisions,
            &mut ledger,
        )?;
    }
    // SEED-CARRY-CC: continent-continent overlap losers (destination,
    // loser source, mean colliding thickness), and which sources lost one
    // (any law; diagnostics).
    let mut cc_losers: Vec<(usize, usize, f64)> = Vec::new();
    let mut cc_loser_src = vec![false; n];

    let carry = |out: &mut CrustState, dst: usize, src: usize, moved_pos: [f64; 3]| {
        out.material_pos[dst] = moved_pos;
        out.plate[dst] = state.plate[src];
        out.kind[dst] = state.kind[src];
        out.birth_t_my[dst] = state.birth_t_my[src];
        out.thickness_m[dst] = state.thickness_m[src];
        out.sediment_m[dst] = state.sediment_m[src];
        out.eroded_total_m[dst] = state.eroded_total_m[src];
        out.arc_added_m[dst] = state.arc_added_m[src];
        out.shortening_added_m[dst] = state.shortening_added_m[src];
    };

    type Claim = (PlateId, usize, [f64; 3]);
    // Legacy claims of one cell: per live plate (ascending id), the source
    // is the cell whose MATERIAL is nearest the back-traced point; its
    // material moves forward by the plate's rotation to become the cell's
    // new phase.
    // SEED-PERF (lane SP): the same claims, bit for bit: the per-plate
    // rotations are prepared once (their validation never depends on the
    // point, and each error surfaces at the same call as before), the
    // source search is `MaterialSearch` (plate prefilter + chord ranking
    // with exact near-tie resolution), and every cell's claims are
    // computed up front in parallel, then consumed in ascending cell order
    // (an error surfaces where the serial loop would have met it).
    let back_rot: Vec<Result<crate::rotation::AxisRotation, RotationError>> = back
        .iter()
        .map(|(_, e)| crate::rotation::AxisRotation::euler(e, dt_my))
        .collect();
    let fwd_rot: Vec<Result<crate::rotation::AxisRotation, RotationError>> = live
        .iter()
        .map(|(_, e)| crate::rotation::AxisRotation::euler(e, dt_my))
        .collect();
    let search = MaterialSearch::new(grid, &state.material_pos, &state.plate)?;
    // Plates that cannot claim a cell (none of their cells within reach of
    // it) are skipped. Only when every rotation is valid and the fast
    // search is active, so no error is skipped with them.
    let reach_masks: Option<Vec<u64>> = if search.tol.is_some()
        && back_rot.iter().chain(&fwd_rot).all(|r| r.is_ok())
    {
        let ids: Vec<PlateId> = live.iter().map(|(id, _)| *id).collect();
        plate_reach_masks(grid, &state.plate, &ids, claim_reach_radius(grid, &live, dt_my))
    } else {
        None
    };
    // Plates with bit-identical poles (a sutured group shares one) have the
    // same back-traced point, so one window search serves the group: its
    // winner is the legacy `nearest_material` source, kept by whichever
    // member owns it (at most one). Only with the masks active (all
    // rotations valid, fast search), groups of up to 64 plates.
    let group: Vec<usize> = live
        .iter()
        .map(|(_, e)| {
            live.iter()
                .position(|(_, f)| {
                    (f.lat_deg.to_bits(), f.lon_deg.to_bits(), f.rate_deg_per_my.to_bits())
                        == (e.lat_deg.to_bits(), e.lon_deg.to_bits(), e.rate_deg_per_my.to_bits())
                })
                .expect("self")
        })
        .collect();
    let grouped = reach_masks.is_some() && group.iter().enumerate().any(|(k, &g)| g != k);
    let group_ids: Vec<Vec<PlateId>> = (0..live.len())
        .map(|g| live.iter().enumerate().filter(|(k, _)| group[*k] == g).map(|(_, (id, _))| *id).collect())
        .collect();
    let legacy_claims = |c: usize| -> Result<Vec<Claim>, KinematicError> {
        let mut claims: Vec<Claim> = Vec::new();
        // Per group (indexed by its first member): Some(winner or None).
        let mut memo: Vec<Option<Option<usize>>> = if grouped { vec![None; live.len()] } else { Vec::new() };
        for (k, (id, _)) in live.iter().enumerate() {
            if let Some(m) = &reach_masks {
                if m[c] & (1u64 << k) == 0 {
                    continue;
                }
            }
            if grouped && group_ids[group[k]].len() > 1 {
                let g = group[k];
                let w = match memo[g] {
                    Some(w) => w,
                    None => {
                        let q = back_rot[k].clone()?.apply(cents[c]);
                        let ids = &group_ids[g];
                        let w = search.winner_if(q, |p| ids.contains(&p))?;
                        memo[g] = Some(w);
                        w
                    }
                };
                if let Some(s) = w.filter(|&s| state.plate[s] == *id) {
                    let moved = fwd_rot[k].clone()?.apply(state.material_pos[s]);
                    claims.push((*id, s, moved));
                }
                continue;
            }
            let q = back_rot[k].clone()?.apply(cents[c]);
            if let Some(s) = search.claim(q, *id)? {
                let moved = fwd_rot[k].clone()?.apply(state.material_pos[s]);
                claims.push((*id, s, moved));
            }
        }
        Ok(claims)
    };
    let mut all_claims: Vec<Result<Vec<Claim>, KinematicError>> = par::map_indexed(n, legacy_claims);
    // SEED-CARRY-CC `front-rescue-carry-v1`: all legacy claims first, then
    // the rescue claims of continental parcels the legacy carry would carry
    // nowhere (`carry.rs`). Legacy runs compute claims inside the loop as
    // before (bit for bit).
    let mut pre: Vec<Vec<Claim>> = Vec::new();
    if contact {
        pre.reserve(n);
        for r in std::mem::take(&mut all_claims) {
            pre.push(r?);
        }
        let mut carried0 = vec![0u32; n];
        for claims in &pre {
            if let Some(s) = legacy_winner_source(claims, state) {
                carried0[s] += 1;
            }
        }
        let fwd: BTreeMap<PlateId, EulerPoleV0> = live.iter().cloned().collect();
        let bwd: BTreeMap<PlateId, EulerPoleV0> = back.iter().cloned().collect();
        let nlat_i = grid.nlat() as i64;
        let nlon_i = grid.nlon() as i64;
        // Best rescue per (cell, plate): (distance, source, moved).
        let mut best: BTreeMap<(usize, PlateId), (f64, usize, [f64; 3])> = BTreeMap::new();
        for s in 0..n {
            if carried0[s] != 0 || state.kind[s] != CrustKindCodeV0::Continental {
                continue;
            }
            let pid = state.plate[s];
            let (Some(fp), Some(bp)) = (fwd.get(&pid), bwd.get(&pid)) else {
                continue;
            };
            let moved = rotate_by_euler(state.material_pos[s], fp, dt_my)?;
            let c0 = nearest_cell(grid, moved)?;
            let (ci, cj) = ((c0 / grid.nlon()) as i64, (c0 % grid.nlon()) as i64);
            let mut pick: Option<(f64, usize)> = None;
            for di in -1..=1i64 {
                let i = ci + di;
                if i < 0 || i >= nlat_i {
                    continue;
                }
                for dj in -1..=1i64 {
                    let c = (i * nlon_i + (cj + dj).rem_euclid(nlon_i)) as usize;
                    if pre[c].iter().any(|(p, _, _)| *p == pid) {
                        continue;
                    }
                    let q = rotate_by_euler(cents[c], bp, dt_my)?;
                    if !in_cell_patch(grid, state.material_pos[s], s, q)? {
                        continue;
                    }
                    let d = great_circle_rad(q, state.material_pos[s])?;
                    if pick.map_or(true, |(bd, bc)| d < bd || (d == bd && c < bc)) {
                        pick = Some((d, c));
                    }
                }
            }
            if let Some((d, c)) = pick {
                let e = best.entry((c, pid)).or_insert((d, s, moved));
                if d < e.0 || (d == e.0 && s < e.1) {
                    *e = (d, s, moved);
                }
            }
        }
        for ((c, pid), (_, s, moved)) in best {
            pre[c].push((pid, s, moved));
            pre[c].sort_by_key(|(p, _, _)| *p);
            ledger.carry.rescues_admitted += 1;
        }
    }
    for c in 0..n {
        let claims: Vec<Claim> = if contact {
            std::mem::take(&mut pre[c])
        } else {
            std::mem::replace(&mut all_claims[c], Ok(Vec::new()))?
        };
        match claims.len() {
            1 => {
                carry(&mut out, c, claims[0].1, claims[0].2);
                source[c] = claims[0].1 as u32;
            }
            0 => {
                // Ridge birth: new oceanic crust joins the previous occupant.
                out.material_pos[c] = cents[c];
                out.plate[c] = state.plate[c];
                out.kind[c] = CrustKindCodeV0::Oceanic;
                out.birth_t_my[c] = t_new;
                out.thickness_m[c] = prior.record.oceanic_crust_thickness_m;
                ledger.ridge_created_m3 += prior.record.oceanic_crust_thickness_m * areas[c];
                ledger.ridge_cells += 1;
                if let Some(candidates) = formation_candidates.as_mut() {
                    candidates[c] = Some(crate::formation_age::BirthCandidate {
                        origin: crate::formation_age::BirthOrigin::KinematicGap,
                        provisional_gap: true,
                    });
                }
            }
            _ => {
                // Convergence. Overriding rule (declared): a continental
                // claimant overrides oceanic ones; among oceanic claimants
                // the YOUNGEST (least dense) overrides; among continental
                // claimants the lowest plate id is the winner and the
                // others stack onto it (collision). Ties break to the
                // lowest plate id.
                let continental: Vec<&(PlateId, usize, [f64; 3])> = claims
                    .iter()
                    .filter(|(_, s, _)| state.kind[*s] == CrustKindCodeV0::Continental)
                    .collect();
                let winner: (PlateId, usize, [f64; 3]) = if let Some(w) = continental.first() {
                    **w
                } else {
                    let mut best = claims[0];
                    for cl in &claims[1..] {
                        if state.birth_t_my[cl.1] > state.birth_t_my[best.1] {
                            best = *cl;
                        }
                    }
                    best
                };
                carry(&mut out, c, winner.1, winner.2);
                source[c] = winner.1 as u32;
                let wpole = pole_map[&winner.0];
                for (lid, ls, _) in &claims {
                    if *lid == winner.0 {
                        continue;
                    }
                    let lpole = pole_map[lid];
                    // Relative velocity of the loser w.r.t. the winner at
                    // the cell: the direction it advanced along, and the
                    // convergence speed in the overlap (m/My).
                    let v = crate::boundaries::relative_velocity_at(cents[c], &lpole, &wpole)?;
                    let cc = state.kind[*ls] == CrustKindCodeV0::Continental
                        && out.kind[c] == CrustKindCodeV0::Continental;
                    if cc {
                        cc_loser_src[*ls] = true;
                    }
                    if cc && parcel_on {
                        // continental-parcel-carry-v1: the parcels resolve
                        // continent-continent overlaps after the pass.
                        ledger.collision_cells += 1;
                    } else if cc && contact {
                        // front-rescue-carry-v1: record the loser; it is
                        // consumed after the pass only if it is carried
                        // nowhere (then exactly its volume and sediment are
                        // the collision supply).
                        let mean_h = 0.5 * (out.thickness_m[c] + state.thickness_m[*ls]);
                        cc_losers.push((c, *ls, mean_h));
                        ledger.carry.cc_overlap_losers += 1;
                        ledger.collision_cells += 1;
                    } else if state.kind[*ls] == CrustKindCodeV0::Continental
                        && out.kind[c] == CrustKindCodeV0::Continental
                        && belt.is_some()
                    {
                        // rate-limited-belt-v1: the crust supply is the
                        // kinematic shortening flux (computed before the
                        // main pass); the overlapped loser column is left to
                        // the remap and only reported. Its sediment is moved
                        // onto the belt.
                        // accreted-to-belt-crust-v1: the loser's sediment
                        // stays with its column like its crust (no copy).
                        let vol = state.thickness_m[*ls] * areas[c];
                        let sed = if accrete { 0.0 } else { state.sediment_m[*ls] * areas[c] };
                        let mean_h = 0.5 * (out.thickness_m[c] + state.thickness_m[*ls]);
                        collisions.push((c, 0.0, sed, mean_h));
                        ledger.collision_overlap_consumed_m3 += vol;
                        ledger.belt_sediment_m3 += sed;
                        ledger.sediment_stacked_m3 += sed;
                        ledger.collision_cells += 1;
                    } else if state.kind[*ls] == CrustKindCodeV0::Continental
                        && out.kind[c] == CrustKindCodeV0::Continental
                    {
                        // Collision: stack, cap, pass excess outward.
                        let (h, excess) = prior.stack(out.thickness_m[c], state.thickness_m[*ls]);
                        let gained = h - out.thickness_m[c];
                        out.thickness_m[c] = h;
                        out.shortening_added_m[c] += gained;
                        thickening[c] += gained;
                        out.sediment_m[c] += state.sediment_m[*ls];
                        ledger.sediment_stacked_m3 += state.sediment_m[*ls] * areas[c];
                        if state.sediment_m[*ls] > ledger.sediment_moved_max_m {
                            ledger.sediment_moved_max_m = state.sediment_m[*ls];
                            ledger.sediment_moved_max_cell = c as u64;
                        }
                        if excess > 0.0 {
                            let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                            let dir = if norm > 1.0e-15 {
                                [v[0] / norm, v[1] / norm, v[2] / norm]
                            } else {
                                [0.0, 0.0, 0.0]
                            };
                            passes.push((c, excess, dir));
                        }
                        ledger.collision_cells += 1;
                    } else {
                        // Subduction: the loser's column is retired; the
                        // overriding column thickens by arc addition.
                        // SEED-TRANSFORMS (lane TF): on a transform or
                        // divergent class the overlap is remap overlap, not
                        // subduction: the loser is retired with no arc
                        // sample, no subduction count and no subducted
                        // volume (it is remap drift, named in the ledger).
                        if let Some(cls) = &cls0 {
                            let lc = local_class(grid, &cents, cls, c)?;
                            if lc == BoundaryClassCodeV0::Transform.code()
                                || lc == BoundaryClassCodeV0::Ridge.code()
                            {
                                ledger.transform.nonconvergent_overlaps += 1;
                                ledger.transform.nonconvergent_overlap_m3 +=
                                    state.thickness_m[*ls] * areas[c];
                                continue;
                            }
                        }
                        // SEED-OLD-OCEAN (lane OR): log by class and age (no output).
                        if let Some(cls) = &cls0 {
                            let lc = local_class(grid, &cents, cls, c)? as usize;
                            let age = state.t_my - state.birth_t_my[*ls];
                            ledger.transform.subducted_by_class_age_m2[lc.min(4)][age_bin(age)] += areas[c];
                        }
                        ledger.subducted_m3 += state.thickness_m[*ls] * areas[c];
                        ledger.sediment_subducted_m3 += state.sediment_m[*ls] * areas[c];
                        let (rate, input_interface, geometry) = pair_normal_rate(
                            grid,
                            &interfaces,
                            &cents,
                            c,
                            winner.0,
                            *lid,
                            winner.1,
                            *ls,
                            &wpole,
                            &lpole,
                            dt_my,
                        )?;
                        ledger.arc_contact_geometry.push(geometry);
                        ledger
                            .arc_normals
                            .push(pf_contracts::domain::ArcNormalSampleV0 {
                                overlap_element: pf_contracts::spatial::ElementId(c as u64),
                                winner: winner.0,
                                loser: *lid,
                                input_interface,
                                source_epoch_my: state.t_my,
                                normal_rate_m_per_my: rate,
                            });
                        if let Some(p) = belt {
                            // rate-limited-belt-v1: arc crust per unit area
                            // of subducted plate, convergent samples only.
                            if rate > 0.0 {
                                let v = p.arc_productivity_m * areas[c];
                                arcs.push((c, v, winner.0));
                                ledger.arc_created_m3 += v;
                            }
                        } else {
                            let arc = prior.arc_thickening_m(rate.max(0.0), dt_my);
                            out.thickness_m[c] += arc;
                            out.arc_added_m[c] += arc;
                            thickening[c] += arc;
                            ledger.arc_created_m3 += arc * areas[c];
                        }
                        ledger.subduction_cells += 1;
                    }
                }
            }
        }
    }

    // SEED-TRANSFORMS (lane TF) `smoothed-normal-transform-v1`: gap cells
    // on a transform or convergent step-start class take the same side's
    // carried material (no ridge birth), layer by layer.
    if let Some(cls) = &cls0 {
        let transform_code = BoundaryClassCodeV0::Transform.code();
        let mut pending: Vec<(usize, u8)> = Vec::new();
        for c in 0..n {
            if source[c] != NO_SOURCE {
                continue;
            }
            let lc = local_class(grid, &cents, cls, c)?;
            if lc == transform_code
                || lc == BoundaryClassCodeV0::Subduction.code()
                || lc == BoundaryClassCodeV0::Collision.code()
            {
                pending.push((c, lc));
            }
        }
        for _ in 0..crate::boundaries::TRANSFORM_FILL_LAYERS {
            if pending.is_empty() {
                break;
            }
            let mut fills: Vec<(usize, usize, u8)> = Vec::new();
            let mut rest: Vec<(usize, u8)> = Vec::new();
            for &(c, lc) in &pending {
                let (i, j) = grid
                    .ij(c as u64)
                    .map_err(|e| KinematicError::Grid(e.to_string()))?;
                // Preference: same plate as the previous occupant, then
                // oceanic, then nearest, then lowest id.
                let mut pick: Option<((bool, bool, f64), usize)> = None;
                for (ti, tj, nb) in grid
                    .neighbors(i, j)
                    .map_err(|e| KinematicError::Grid(e.to_string()))?
                {
                    let k = ti * grid.nlon() + tj;
                    if source[k] == NO_SOURCE {
                        continue;
                    }
                    let key = (
                        out.plate[k] != state.plate[c],
                        out.kind[k] != CrustKindCodeV0::Oceanic,
                        nb.distance_m,
                    );
                    if pick.map_or(true, |(bk, bc)| key < bk || (key == bk && k < bc)) {
                        pick = Some((key, k));
                    }
                }
                match pick {
                    Some((_, k)) => fills.push((c, k, lc)),
                    None => rest.push((c, lc)),
                }
            }
            for (c, k, lc) in fills {
                ledger.ridge_created_m3 -= prior.record.oceanic_crust_thickness_m * areas[c];
                ledger.ridge_cells -= 1;
                out.plate[c] = out.plate[k];
                out.kind[c] = out.kind[k];
                out.birth_t_my[c] = out.birth_t_my[k];
                out.thickness_m[c] = out.thickness_m[k];
                out.sediment_m[c] = out.sediment_m[k];
                out.eroded_total_m[c] = out.eroded_total_m[k];
                out.arc_added_m[c] = out.arc_added_m[k];
                out.shortening_added_m[c] = out.shortening_added_m[k];
                source[c] = source[k];
                ledger.transform.filled_m3 += out.thickness_m[c] * areas[c];
                if lc == transform_code {
                    ledger.transform.filled_transform_cells += 1;
                } else {
                    ledger.transform.filled_convergent_cells += 1;
                }
            }
            pending = rest;
        }
        ledger.transform.unfilled_cells = pending.len() as u64;
        // SEED-OLD-OCEAN (lane OR): births left after the fill, by the
        // step-start local class (log only).
        for c in 0..n {
            if source[c] == NO_SOURCE {
                let lc = local_class(grid, &cents, cls, c)? as usize;
                ledger.transform.births_by_class_m2[lc.min(4)] += areas[c];
            }
        }
    }

    // SEED-CARRY-CONSERVE (`continental-parcel-carry-v1`): move the
    // parcels, consume collision overlaps (belt supply), project them over
    // the Eulerian result.
    if let Some((layer, pp)) = parcels.as_mut() {
        // Q27 / D-0155: conserve.rs is a pinned history producer. Observe its
        // vacated births here without changing that producer's bytes. Only a
        // continental cell vacated without a donor becomes oceanic NO_SOURCE;
        // copied vacated cells retain a source. Disabled recording allocates none.
        let was_continental = formation_candidates.as_ref().map(|_| {
            out.kind.iter().map(|k| *k == CrustKindCodeV0::Continental).collect::<Vec<_>>()
        });
        let (proj, rep) = crate::conserve::parcel_carry(
            grid,
            &cents,
            &areas,
            layer,
            pp,
            &pole_map,
            dt_my,
            prior.record.oceanic_crust_thickness_m,
            accrete,
            &mut out,
            &mut source,
            &mut collisions,
            &mut ledger,
        )?;
        if let (Some(was_continental), Some(candidates)) = (was_continental, formation_candidates.as_mut()) {
            for c in 0..n {
                if was_continental[c] && out.kind[c] == CrustKindCodeV0::Oceanic && source[c] == NO_SOURCE {
                    candidates[c] = Some(crate::formation_age::BirthCandidate {
                        origin: crate::formation_age::BirthOrigin::ProjectionInitialization,
                        provisional_gap: candidates[c].is_some_and(|b| b.provisional_gap),
                    });
                }
            }
        }
        ledger.parcels = Some(rep);
        ledger.parcel_projection = Some(proj);
    }
    // SEED-CARRY-CC: conservative consumption (front-rescue-carry-v1) and
    // the carry diagnostics (every law; they feed no output).
    let mut carried = vec![0u32; n];
    for &s in &source {
        if s != NO_SOURCE {
            carried[s as usize] += 1;
        }
    }
    if contact {
        let mut consumed = vec![false; n];
        for &(c, ls, mean_h) in &cc_losers {
            if carried[ls] > 0 {
                ledger.carry.cc_losers_carried += 1;
                continue;
            }
            if consumed[ls] {
                continue;
            }
            consumed[ls] = true;
            let vol = state.thickness_m[ls] * areas[ls];
            let sed = state.sediment_m[ls] * areas[ls];
            if accrete {
                // accreted-to-belt-crust-v1: the consumed loser's sediment
                // joins the belt's crust supply (named, not hidden).
                collisions.push((c, vol + sed, 0.0, mean_h));
                ledger.belt_in_m3 += vol + sed;
                ledger.sediment_accreted_m3 += sed;
            } else {
                collisions.push((c, vol, sed, mean_h));
                ledger.belt_in_m3 += vol;
                ledger.belt_sediment_m3 += sed;
                ledger.sediment_stacked_m3 += sed;
            }
            ledger.collision_overlap_consumed_m3 += vol;
            ledger.carry.cc_consumed_parcels += 1;
            ledger.carry.cc_consumed_area_m2 += areas[ls];
            ledger.carry.cc_consumed_m3 += vol;
            ledger.carry.cc_consumed_sediment_m3 += sed;
        }
    }
    for s in 0..n {
        if carried[s] == 0 && state.kind[s] == CrustKindCodeV0::Continental && !cc_loser_src[s] {
            ledger.carry.continental_dropped_parcels += 1;
            ledger.carry.continental_dropped_area_m2 += areas[s];
            ledger.carry.continental_dropped_m3 += state.thickness_m[s] * areas[s];
        }
    }
    for c in 0..n {
        if source[c] != NO_SOURCE {
            continue;
        }
        grid.ij(c as u64)
            .map_err(|e| KinematicError::Grid(e.to_string()))?;
        let mut first: Option<PlateId> = None;
        let mut two = false;
        let mut nbs = [0usize; 8];
        let nn = neighbor_ids(grid.nlat(), grid.nlon(), c, &mut nbs);
        for &nb in &nbs[..nn] {
            if out.kind[nb] != CrustKindCodeV0::Continental {
                continue;
            }
            match first {
                None => first = Some(out.plate[nb]),
                Some(p) if p != out.plate[nb] => two = true,
                _ => {}
            }
        }
        if two {
            ledger.carry.sliver_births += 1;
        }
    }

    // SEED-EXTRAS (lane EX): rift births, then plume heads and tails, on the
    // post-carry field (before the collisional passes and the belt).
    if let Some(x) = extras {
        if let Some(rp) = &x.rift {
            let birth: Vec<bool> = source.iter().map(|&s| s == NO_SOURCE).collect();
            let rep = crate::rift::rift_pass(
                grid,
                &cents,
                &areas,
                rp,
                &birth,
                &out.plate,
                &pole_map,
                &mut out.kind,
                &mut out.thickness_m,
            )
            .map_err(KinematicError::Grid)?;
            for &b in &rep.births {
                // The birth is stretched continental crust, not new oceanic
                // crust: take it out of the ridge term.
                ledger.ridge_created_m3 -= prior.record.oceanic_crust_thickness_m * areas[b as usize];
                ledger.ridge_cells -= 1;
            }
            for &(c, dh) in &rep.thinning {
                thickening[c as usize] -= dh;
            }
            ledger.rift_births = rep.births;
            ledger.rift_breakup_births = rep.breakup_births;
            ledger.rift_convergent_skipped = rep.convergent_skipped;
            ledger.rift_moved_m3 = rep.moved_m3;
            ledger.rift_thinning = rep.thinning;
            ledger.rift_birth_pre_thickness_m = rep.birth_pre_thickness_m;
            ledger.rift_birth_flanks = rep.birth_flanks;
        }
        if let Some((plumes, pp)) = &x.plumes {
            // SEED-FLOODBASALT (lane FB): the LIP law owns the head's magma.
            let pp_eff = match &x.flood_basalt {
                Some(_) => crate::hotspots::PlumeParams {
                    head_volume_m3: 0.0,
                    ..*pp
                },
                None => *pp,
            };
            let rep = crate::hotspots::plume_pass(
                grid,
                &cents,
                &areas,
                plumes,
                &pp_eff,
                state.t_my,
                t_new,
                &out.plate,
                &pole_map,
                &mut out.thickness_m,
            )
            .map_err(KinematicError::Grid)?;
            for &(c, dh, _) in &rep.cells {
                thickening[c as usize] += dh;
            }
            ledger.hotspot_created_m3 = rep.head_m3 + rep.tail_m3;
            ledger.hotspot_head_m3 = rep.head_m3;
            ledger.hotspot_heads = rep.heads;
            ledger.hotspot_cells = rep.cells;
            // SEED-FLOODBASALT (lane FB) `plume-head-lip-swell-v1`: this
            // step's share of every erupting LIP, on the post-carry field.
            if let Some(lp) = &x.flood_basalt {
                let lrep = crate::flood_basalt::lip_pass(
                    grid,
                    &cents,
                    &areas,
                    plumes,
                    lp,
                    state.t_my,
                    t_new,
                    &out.plate,
                    &out.kind,
                    &pole_map,
                    &mut out.thickness_m,
                )
                .map_err(KinematicError::Grid)?;
                for &(c, dh, _) in &lrep.cells {
                    thickening[c as usize] += dh;
                }
                ledger.hotspot_created_m3 += lrep.total_m3;
                ledger.hotspot_head_m3 += lrep.total_m3;
                ledger.lip_created_m3 = lrep.total_m3;
                ledger.lip_continental_m3 = lrep.continental_m3;
                ledger.lip_cells = lrep.cells;
            }
            crate::volcanoes::supply_capture::plume_control(
                rep.head_m3, rep.tail_m3, ledger.lip_created_m3,
            ).map_err(KinematicError::Grid)?;
        }
    }

    // Pass collisional excess outward along the convergence direction,
    // ascending cell order, bounded walk, remainder delaminated.
    for (start, excess, dir) in passes {
        let mut remaining = excess * areas[start];
        ledger.wedge_in_m3 += remaining;
        let mut cur = start;
        let mut hops = 0u32;
        while remaining > 0.0 && hops < prior.record.wedge_pass_cells {
            let (i, j) = grid
                .ij(cur as u64)
                .map_err(|e| KinematicError::Grid(e.to_string()))?;
            // Neighbour most aligned with the direction (ties: lowest id).
            let mut best: Option<(f64, usize)> = None;
            for (ti, tj, _) in grid
                .neighbors(i, j)
                .map_err(|e| KinematicError::Grid(e.to_string()))?
            {
                let nb = ti * grid.nlon() + tj;
                let d = [
                    cents[nb][0] - cents[cur][0],
                    cents[nb][1] - cents[cur][1],
                    cents[nb][2] - cents[cur][2],
                ];
                let dn = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                let align = if dn > 0.0 {
                    (d[0] * dir[0] + d[1] * dir[1] + d[2] * dir[2]) / dn
                } else {
                    -1.0
                };
                let replace = match best {
                    None => true,
                    Some((ba, bn)) => align > ba || (align == ba && nb < bn),
                };
                if replace {
                    best = Some((align, nb));
                }
            }
            let Some((_, next)) = best else { break };
            // Volume-preserving transfer between cells of different area.
            let (h, rest, received) =
                receive_wedge_volume(prior, out.thickness_m[next], remaining, areas[next]);
            let gained = h - out.thickness_m[next];
            ledger.wedge_stacked_m3 += received;
            out.thickness_m[next] = h;
            out.shortening_added_m[next] += gained;
            thickening[next] += gained;
            remaining = rest;
            cur = next;
            hops += 1;
        }
        if remaining > 0.0 {
            ledger.delaminated_m3 += remaining;
            ledger.wedge_delaminated_m3 += remaining;
        }
    }

    if let Some(p) = belt {
        belt_pass(
            grid,
            &cents,
            &areas,
            p,
            dt_my,
            &collisions,
            &arcs,
            &mut out,
            &mut ledger,
            &mut thickening,
            plateau,
        )?;
    } else if plateau.is_some() {
        return Err(KinematicError::Grid("plateau-cap-spill-v2 requires rate-limited-belt-v1; refused".into()));
    }

    let wedge_residual = ledger.wedge_in_m3 - ledger.wedge_stacked_m3 - ledger.wedge_delaminated_m3;
    if wedge_residual.abs() > 1e-12 * ledger.wedge_in_m3.max(1.) {
        return Err(KinematicError::Grid(format!(
            "collisional wedge volume does not close: {wedge_residual} m3"
        )));
    }
    // SEED-BREAKUP-CUT-V2 (lane BC) `symmetric-ridge-accretion-v1`: split
    // this step's gap cells between the flanking plates.
    if ridge == RidgeAccretion::SymmetricV1 {
        let moved = symmetric_ridge_relabel(grid, &source, &mut out)?;
        if let Some((layer, _)) = parcels.as_mut() {
            crate::conserve::relabel_parcels(grid, layer, &moved)?;
        }
        ledger.ridge_relabelled_cells = moved.len() as u64;
    }
    // SEED-CARRY-CONSERVE: fold plate-id speckle into the surrounding plate
    // (the parcels in folded cells follow).
    if let Some((layer, pp)) = parcels.as_mut() {
        let total: f64 = areas.iter().sum();
        let (folded, moved) = crate::conserve::fold_fragments(
            grid,
            &areas,
            &mut out.plate,
            pp.fragment_min_area_fraction * total,
        )?;
        let relabelled = crate::conserve::relabel_parcels(grid, layer, &moved)?;
        if let Some(r) = ledger.parcels.as_mut() {
            r.fragments_folded = folded;
            r.fragment_cells = moved.len() as u64;
            r.fragment_parcels = relabelled;
        }
    }
    ledger.volume_after_m3 = out.crust_volume_m3(&areas);
    ledger.sediment_after_m3 = out.sediment_m.iter().zip(&areas).map(|(h, a)| h * a).sum();
    let (class, rate) = match transform {
        None => {
            let (class, rate, _) = classify_cells(grid, &out, plates)?;
            (class, rate)
        }
        Some(tp) => {
            let (class, rate, _, fb) = classify_cells_v2(grid, &out, plates, tp)?;
            ledger.transform.normal_fallback_cells = fb;
            ledger.transform.boundary_cells =
                class.iter().filter(|&&k| k != BoundaryClassCodeV0::Interior.code()).count() as u64;
            (class, rate)
        }
    };
    out.boundary_class = class;
    out.normal_rate_m_per_my = rate;
    if belt.is_some() {
        ledger.max_thickness_m = out.thickness_m.iter().cloned().fold(0.0, f64::max);
    }
    if let Some(recorder) = recorder {
        recorder.accept_step(state, &out, &source, formation_candidates.as_ref().unwrap(),
            &areas, prior.record.oceanic_crust_thickness_m).map_err(KinematicError::Grid)?;
    }
    Ok((out, ledger, thickening, source))
}

/// SEED-BREAKUP-CUT-V2 (lane BC) `symmetric-ridge-accretion-v1` (see
/// `RidgeAccretion::SymmetricV1`). Returns (cell, old plate, new plate) for
/// every relabelled cell, ascending.
fn symmetric_ridge_relabel(
    grid: &GridSpec,
    source: &[u32],
    out: &mut CrustState,
) -> Result<Vec<(usize, PlateId, PlateId)>, KinematicError> {
    let n = grid.cells();
    let nlon = grid.nlon();
    let gap: Vec<bool> = source.iter().map(|&s| s == NO_SOURCE).collect();
    let nbrs = |c: usize| -> Result<Vec<usize>, KinematicError> {
        let (i, j) = grid.ij(c as u64).map_err(|e| KinematicError::Grid(e.to_string()))?;
        Ok(grid
            .neighbors(i, j)
            .map_err(|e| KinematicError::Grid(e.to_string()))?
            .into_iter()
            .map(|(ti, tj, _)| ti * nlon + tj)
            .collect())
    };
    let t = out.t_my;
    let mut label: Vec<Option<PlateId>> = vec![None; n];
    // Pick among (plate, score) candidates: highest score; ties keep the
    // previous occupant, else the lowest id.
    let pick = |cands: &BTreeMap<PlateId, f64>, prev: PlateId| -> Option<PlateId> {
        let best = cands.values().cloned().fold(f64::NEG_INFINITY, f64::max);
        if !best.is_finite() && best != f64::INFINITY {
            return None;
        }
        let tied: Vec<PlateId> = cands.iter().filter(|(_, &v)| v == best).map(|(p, _)| *p).collect();
        if tied.contains(&prev) {
            Some(prev)
        } else {
            tied.first().copied()
        }
    };
    // Layer 1: gap cells touching carried crust.
    let mut frontier: Vec<usize> = Vec::new();
    for c in 0..n {
        if !gap[c] {
            continue;
        }
        let mut cands: BTreeMap<PlateId, f64> = BTreeMap::new();
        for k in nbrs(c)? {
            if gap[k] {
                continue;
            }
            let age = if out.kind[k] == CrustKindCodeV0::Continental { f64::INFINITY } else { t - out.birth_t_my[k] };
            let e = cands.entry(out.plate[k]).or_insert(f64::INFINITY);
            *e = e.min(age);
        }
        if let Some(p) = pick(&cands, out.plate[c]) {
            label[c] = Some(p);
            frontier.push(c);
        }
    }
    // Deeper layers: majority of labelled gap neighbours.
    while !frontier.is_empty() {
        let mut next: Vec<(usize, PlateId)> = Vec::new();
        let mut seen: Vec<usize> = Vec::new();
        for &f in &frontier {
            for k in nbrs(f)? {
                if gap[k] && label[k].is_none() {
                    seen.push(k);
                }
            }
        }
        seen.sort_unstable();
        seen.dedup();
        for c in seen {
            let mut cands: BTreeMap<PlateId, f64> = BTreeMap::new();
            for k in nbrs(c)? {
                if let (true, Some(p)) = (gap[k], label[k]) {
                    *cands.entry(p).or_insert(0.0) += 1.0;
                }
            }
            if let Some(p) = pick(&cands, out.plate[c]) {
                next.push((c, p));
            }
        }
        frontier = next.iter().map(|x| x.0).collect();
        for (c, p) in next {
            label[c] = Some(p);
        }
    }
    let mut moved = Vec::new();
    for c in 0..n {
        if let Some(p) = label[c] {
            if p != out.plate[c] {
                moved.push((c, out.plate[c], p));
                out.plate[c] = p;
            }
        }
    }
    Ok(moved)
}

/// The source the legacy overriding rule carries into a cell from its
/// claims (sole claimant; else the first continental claimant; else the
/// youngest oceanic claimant, ties to the lowest plate id); `None` for a
/// birth. Mirrors the resolution in `advance_impl` (SEED-CARRY-CC).
fn legacy_winner_source(
    claims: &[(PlateId, usize, [f64; 3])],
    state: &CrustState,
) -> Option<usize> {
    match claims.len() {
        0 => None,
        1 => Some(claims[0].1),
        _ => {
            if let Some(w) = claims
                .iter()
                .find(|(_, s, _)| state.kind[*s] == CrustKindCodeV0::Continental)
            {
                return Some(w.1);
            }
            let mut best = claims[0];
            for cl in &claims[1..] {
                if state.birth_t_my[cl.1] > state.birth_t_my[best.1] {
                    best = *cl;
                }
            }
            Some(best.1)
        }
    }
}

/// SEED-OROGEN `rate-limited-belt-v1` collisional supply: the kinematic
/// shortening flux across every step-start interface edge between
/// continental cells of different plates that converge there:
/// `V = h_mean * v_n * dt * edge_length` (h_mean the mean of the two
/// columns, v_n the normal convergence rate at the edge midpoint, edge
/// length the shared great-circle cell edge). Half of V is placed on each
/// side's belt. Per unit time this is independent of dt (unlike overlap
/// counting, which the semi-Lagrangian carry makes dt-dependent).
/// With `place == false` (SEED-CARRY-CC, `front-rescue-carry-v1`) the
/// edges are only counted into `collision_supply_edges` (contact-length
/// diagnostic) and no supply is pushed.
#[allow(clippy::too_many_arguments)]
fn kinematic_collision_supply(
    grid: &GridSpec,
    state: &CrustState,
    poles: &BTreeMap<PlateId, EulerPoleV0>,
    cents: &[[f64; 3]],
    dt_my: f64,
    place: bool,
    collisions: &mut Vec<(usize, f64, f64, f64)>,
    ledger: &mut AdvanceLedger,
) -> Result<(), KinematicError> {
    let radius = grid.radius_m();
    let dlat = std::f64::consts::PI / grid.nlat() as f64;
    let dlon = 2.0 * std::f64::consts::PI / grid.nlon() as f64;
    for c in 0..state.cells() {
        if state.kind[c] != CrustKindCodeV0::Continental {
            continue;
        }
        let (i, j) = grid
            .ij(c as u64)
            .map_err(|e| KinematicError::Grid(e.to_string()))?;
        let lat_bot = -std::f64::consts::FRAC_PI_2 + i as f64 * dlat;
        for (ti, tj, nb) in grid
            .neighbors(i, j)
            .map_err(|e| KinematicError::Grid(e.to_string()))?
        {
            if nb.di.abs() + nb.dj.abs() != 1 {
                continue;
            }
            let k = ti * grid.nlon() + tj;
            if k <= c
                || state.plate[k] == state.plate[c]
                || state.kind[k] != CrustKindCodeV0::Continental
            {
                continue;
            }
            let pa = poles
                .get(&state.plate[c])
                .ok_or(KinematicError::UnknownPlate {
                    plate: state.plate[c],
                })?;
            let pb = poles
                .get(&state.plate[k])
                .ok_or(KinematicError::UnknownPlate {
                    plate: state.plate[k],
                })?;
            let raw = [
                cents[c][0] + cents[k][0],
                cents[c][1] + cents[k][1],
                cents[c][2] + cents[k][2],
            ];
            let norm = (raw.iter().map(|x| x * x).sum::<f64>()).sqrt();
            if norm < 1e-12 {
                continue;
            }
            let mid = [raw[0] / norm, raw[1] / norm, raw[2] / norm];
            let vn = normal_rate_at(mid, cents[c], cents[k], pa, pb)? * radius;
            if !(vn > 0.0) {
                continue;
            }
            let edge = if nb.di == 0 {
                radius * dlat
            } else if nb.di > 0 {
                radius * (lat_bot + dlat).cos() * dlon
            } else {
                radius * lat_bot.cos() * dlon
            };
            ledger.collision_supply_edges += 1;
            if !place {
                continue;
            }
            let h = 0.5 * (state.thickness_m[c] + state.thickness_m[k]);
            let v = h * vn * dt_my * edge;
            collisions.push((c, 0.5 * v, 0.0, h));
            collisions.push((k, 0.5 * v, 0.0, h));
            ledger.belt_in_m3 += v;
        }
    }
    Ok(())
}

/// SEED-OROGEN `rate-limited-belt-v1`: place the step's arc and collision
/// volumes on their belts, integrate lateral channel-flow spreading with the
/// collision supply in rate form, then apply the delamination backstop.
/// See `orogen.rs` for the physics. Updates the lineage fields so
/// thickness = T0 + shortening_added + arc_added - eroded still holds per
/// carried column.
#[allow(clippy::too_many_arguments)]
fn belt_pass(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    p: &crate::orogen::BeltParams,
    dt_my: f64,
    collisions: &[(usize, f64, f64, f64)],
    arcs: &[(usize, f64, PlateId)],
    out: &mut CrustState,
    ledger: &mut AdvanceLedger,
    thickening: &mut [f64],
    plateau: Option<&crate::ranges::PlateauInBeltParams>,
) -> Result<(), KinematicError> {
    use crate::orogen::{disc_weights, place, spread, EdgeStencil};
    let n = out.cells();
    let g = |e: String| KinematicError::Grid(e);
    // Arc: immediate, on overriding-plate cells within the arc belt width.
    // SEED-PERF (lane SP): every arc and collision disc is computed up
    // front in parallel (pure: the arc placement below changes neither
    // `out.plate` nor `out.kind`), then placed serially in event order.
    let arc_discs = {
        let plates = &out.plate;
        par::map_heavy(arcs.len(), |e| {
            let (c, _, plate) = arcs[e];
            disc_weights(grid, cents, c, p.arc_belt_width_m, &|k| plates[k] == plate)
        })
    };
    let coll_discs = {
        let kind = &out.kind;
        par::map_heavy(collisions.len(), |e| {
            let (c, _, _, mean_h) = collisions[e];
            disc_weights(grid, cents, c, 0.5 * p.belt_width_m(mean_h), &|k| {
                kind[k] == CrustKindCodeV0::Continental
            })
        })
    };
    let mut arc_dh = vec![0.0f64; n];
    for (&(_, v, _), w) in arcs.iter().zip(arc_discs) {
        let w = w.map_err(g)?;
        ledger.arc_placed_m3 += place(&w, areas, v, &mut arc_dh);
    }
    for c in 0..n {
        if arc_dh[c] != 0.0 {
            out.thickness_m[c] += arc_dh[c];
            out.arc_added_m[c] += arc_dh[c];
            thickening[c] += arc_dh[c];
        }
    }
    // Collision supply: belt increments for the whole step.
    let mut src = vec![0.0f64; n];
    let mut sed = vec![0.0f64; n];
    for (&(_, v, s, _), w) in collisions.iter().zip(coll_discs) {
        let w = w.map_err(g)?;
        ledger.belt_placed_m3 += place(&w, areas, v, &mut src);
        place(&w, areas, s, &mut sed);
    }
    for c in 0..n {
        out.sediment_m[c] += sed[c];
        if sed[c] > ledger.sediment_moved_max_m {
            ledger.sediment_moved_max_m = sed[c];
            ledger.sediment_moved_max_cell = c as u64;
        }
    }
    // SEED-WORLD-COMPOSE (lane CP) `plateau-cap-spill-v2`: the supply is
    // placed now, the plateau pass caps, spills and tapers, and the channel
    // flow then spreads the capped field (no source).
    if let Some(pp) = plateau {
        let before = out.thickness_m.clone();
        for c in 0..n {
            out.thickness_m[c] += src[c];
        }
        let placed = out.thickness_m.clone();
        let prep = crate::ranges::plateau_in_belt(grid, cents, areas, &out.kind, &out.plate, &mut out.thickness_m, &mut out.sediment_m, pp)?;
        ledger.sediment_accreted_m3 += prep.sediment_accreted_m3;
        let pdel = prep.delaminated_m3 + prep.oceanic_delaminated_m3;
        ledger.belt_delaminated_m3 += pdel;
        ledger.delaminated_m3 += pdel;
        let capped = out.thickness_m.clone();
        let stencil = EdgeStencil::build(grid).map_err(g)?;
        let zero = vec![0.0f64; n];
        // PREREG revision 1 (lane CP): the polar kappa limiter.
        let ks = crate::orogen::polar_kappa_scale(grid).map_err(g)?;
        let rep = crate::orogen::spread_scaled(&stencil, &out.kind, &mut out.thickness_m, &zero, p, dt_my, Some(&ks)).map_err(g)?;
        let mut budget = BeltBudget {
            placed_m: src.clone(),
            plateau_dh_m: (0..n).map(|c| capped[c] - placed[c]).collect(),
            spread_dh_m: (0..n).map(|c| out.thickness_m[c] - capped[c]).collect(),
            delaminated_m: vec![0.0; n],
        };
        for c in 0..n {
            let mut net = out.thickness_m[c] - before[c];
            if out.thickness_m[c] > p.delamination_thickness_m {
                let excess = out.thickness_m[c] - p.delamination_thickness_m;
                out.thickness_m[c] = p.delamination_thickness_m;
                ledger.belt_delaminated_m3 += excess * areas[c];
                ledger.delaminated_m3 += excess * areas[c];
                budget.delaminated_m[c] = excess;
                net -= excess;
            }
            out.shortening_added_m[c] += net;
            thickening[c] += net;
        }
        ledger.spread = Some(rep);
        ledger.plateau_belt = Some(prep);
        ledger.belt_budget = Some(budget);
        return Ok(());
    }
    // Spreading with the supply in rate form.
    let before = out.thickness_m.clone();
    let stencil = EdgeStencil::build(grid).map_err(g)?;
    let rep = spread(&stencil, &out.kind, &mut out.thickness_m, &src, p, dt_my).map_err(g)?;
    let mut budget = BeltBudget {
        placed_m: src.clone(),
        plateau_dh_m: vec![0.0; n],
        spread_dh_m: (0..n).map(|c| out.thickness_m[c] - before[c] - src[c]).collect(),
        delaminated_m: vec![0.0; n],
    };
    // Delamination backstop.
    for c in 0..n {
        let mut net = out.thickness_m[c] - before[c];
        // Any crust kind (thick oceanic arc roots delaminate too).
        if out.thickness_m[c] > p.delamination_thickness_m {
            let excess = out.thickness_m[c] - p.delamination_thickness_m;
            out.thickness_m[c] = p.delamination_thickness_m;
            ledger.belt_delaminated_m3 += excess * areas[c];
            ledger.delaminated_m3 += excess * areas[c];
            budget.delaminated_m[c] = excess;
            net -= excess;
        }
        out.shortening_added_m[c] += net;
        thickening[c] += net;
    }
    ledger.spread = Some(rep);
    ledger.belt_budget = Some(budget);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::seed_plates;

    /// SEED-BREAKUP-CUT-V2 (lane BC): `symmetric-ridge-accretion-v1` on a
    /// hand-made gap between plate 1 (west, lon < 10) and plate 2 (east).
    /// Gap columns are labelled plate 1 (the legacy previous occupant).
    fn ridge_case(gap_cols: &[usize], west_young: bool) -> Vec<PlateId> {
        let grid = GridSpec::new(6, 24, 6.371e6).unwrap();
        let n = grid.cells();
        let nlon = 24;
        let t = 10.0;
        let mut st = CrustState {
            t_my: t,
            plate: vec![PlateId(2); n],
            kind: vec![CrustKindCodeV0::Oceanic; n],
            birth_t_my: vec![-50.0; n],
            thickness_m: vec![7000.0; n],
            sediment_m: vec![0.0; n],
            eroded_total_m: vec![0.0; n],
            arc_added_m: vec![0.0; n],
            shortening_added_m: vec![0.0; n],
            boundary_class: vec![0; n],
            normal_rate_m_per_my: vec![0.0; n],
            material_pos: vec![[0.0; 3]; n],
        };
        let mut source = vec![0u32; n];
        for c in 0..n {
            let j = c % nlon;
            if j < 10 {
                st.plate[c] = PlateId(1);
                if west_young && j == 9 {
                    st.birth_t_my[c] = t; // the west flank accreted last step
                }
            }
            if gap_cols.contains(&j) {
                st.plate[c] = PlateId(1);
                st.birth_t_my[c] = t;
                source[c] = NO_SOURCE;
            }
        }
        symmetric_ridge_relabel(&grid, &source, &mut st).unwrap();
        (0..nlon).map(|j| st.plate[2 * nlon + j]).collect()
    }

    #[test]
    fn symmetric_ridge_splits_gaps_between_flanks() {
        // One-cell gap: the flank that accreted last step (young west) loses
        // the cell; equal flanks keep the previous occupant.
        assert_eq!(ridge_case(&[10], true)[10], PlateId(2));
        assert_eq!(ridge_case(&[10], false)[10], PlateId(1));
        // Four-cell gap: two cells each side.
        let r = ridge_case(&[10, 11, 12, 13], false);
        assert_eq!(&r[10..14], &[PlateId(1), PlateId(1), PlateId(2), PlateId(2)]);
        // Untouched outside the gap.
        assert!(r[..10].iter().all(|p| *p == PlateId(1)) && r[14..].iter().all(|p| *p == PlateId(2)));
    }
    use pf_contracts::params::{FidelityPreset, MantleBlock, ScenarioPreset, SeedPath};

    fn setup(plate_count: u32) -> (GridSpec, Vec<PlateRecordV0>, CrustState, UpliftPrior) {
        let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
        let mut p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
        p.tectonics.plate_count = plate_count;
        p.tectonics.continental_seeds = 2;
        let sp = SeedPath {
            seed: 240907,
            derivation: "kinematic-test".into(),
            salt: "k".into(),
        };
        let seeding = seed_plates(&grid, &p, &sp).unwrap();
        let prior = UpliftPrior::preview_v0();
        let state = CrustState::seed(&grid, &seeding, &prior).unwrap();
        (grid, seeding.plates, state, prior)
    }

    #[test]
    fn wedge_volume_crosses_unequal_areas_without_creating_material() {
        let prior = UpliftPrior::preview_v0();
        let cap = prior.record.max_crust_thickness_m;
        // Ten metres of excess over source area1 crosses a full area2 cell
        // and reaches area3. The old area rebasing offered20m3 at that hop.
        let incoming = 10. * 1.;
        let (middle, left, a) = receive_wedge_volume(&prior, cap, incoming, 2.);
        assert_eq!(middle, cap);
        assert_eq!(left, 10.);
        assert_eq!(a, 0.);
        let (last, left, b) = receive_wedge_volume(&prior, cap - 100., left, 3.);
        assert!(((last - (cap - 100.)) * 3. - 10.).abs() < 1e-9);
        assert_eq!(left, 0.);
        assert_eq!(a + b, incoming);
        let (_, rest, accepted) = receive_wedge_volume(&prior, cap + 500., 10., 3.);
        assert_eq!(accepted, 0.);
        assert_eq!(rest, 10.);
    }

    #[test]
    fn pair_normal_is_geometric_and_unchanged_by_shear() {
        let (grid, _, state, _) = setup(3);
        let cents = centroids(&grid).unwrap();
        let pairs = interface_samples(&grid, &state, &cents).unwrap();
        let (&(a, b), edges) = pairs.iter().next().unwrap();
        let (i, j) = edges[0];
        let raw = [
            cents[i][0] + cents[j][0],
            cents[i][1] + cents[j][1],
            cents[i][2] + cents[j][2],
        ];
        let l = raw.iter().map(|v| v * v).sum::<f64>().sqrt();
        let point = [raw[0] / l, raw[1] / l, raw[2] / l];
        let normal = crate::boundaries::boundary_normal(point, cents[i], cents[j]).unwrap();
        // omega = r cross n produces unit normal velocity; adding n to omega
        // produces a tangential velocity perpendicular to that normal.
        let omega = [
            point[1] * normal[2] - point[2] * normal[1],
            point[2] * normal[0] - point[0] * normal[2],
            point[0] * normal[1] - point[1] * normal[0],
        ];
        let pole = |v: [f64; 3]| {
            let mag = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            let (lat, lon) = unit_to_lat_lon([v[0] / mag, v[1] / mag, v[2] / mag]).unwrap();
            EulerPoleV0 {
                lat_deg: lat,
                lon_deg: lon,
                rate_deg_per_my: mag,
            }
        };
        let base = pole(omega);
        let shear = pole([
            omega[0] + normal[0],
            omega[1] + normal[1],
            omega[2] + normal[2],
        ]);
        let zero = EulerPoleV0 {
            lat_deg: 90.,
            lon_deg: 0.,
            rate_deg_per_my: 0.,
        };
        // A singleton interface isolates the exact event patch from nearest-edge selection.
        let patch = InterfaceSamples::from([((a, b), vec![(i, j)])]);
        let (x, e, g) =
            pair_normal_rate(&grid, &patch, &cents, i, a, b, i, j, &base, &zero, 2.).unwrap();
        assert_eq!(g, ArcContactGeometry::StepStartInterface);
        let (y, e2, _) =
            pair_normal_rate(&grid, &patch, &cents, i, a, b, i, j, &shear, &zero, 2.).unwrap();
        assert_eq!(e, e2);
        assert!((x - y).abs() < 1e-8);
        assert!(x > 0.);
        let transform = pole(normal);
        let (v, _, _) =
            pair_normal_rate(&grid, &patch, &cents, i, a, b, i, j, &transform, &zero, 2.).unwrap();
        assert!(v.abs() < 1e-8);
        let (swapped, _, _) =
            pair_normal_rate(&grid, &patch, &cents, i, b, a, j, i, &zero, &base, 2.).unwrap();
        assert!((swapped - x).abs() < 1e-8);
    }

    /// SEED-ROBUST: with no step-start interface of the pair, the normal is
    /// taken across the two sources that met; the rate is the hand-built
    /// convergence speed, tangential shear does not change it, the record
    /// lists the lower plate's source first, and nonlocal or degenerate
    /// source pairs refuse with a typed reason.
    #[test]
    fn new_contact_samples_across_met_material_and_refuses_typed() {
        let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
        let cents = centroids(&grid).unwrap();
        let (a, b) = (PlateId(4), PlateId(1));
        // Winner source two columns west of the loser source, same row.
        let (sa, sb) = (8 * 36 + 10, 8 * 36 + 12);
        let cell = 8 * 36 + 11;
        let raw = [
            cents[sa][0] + cents[sb][0],
            cents[sa][1] + cents[sb][1],
            cents[sa][2] + cents[sb][2],
        ];
        let l = raw.iter().map(|v| v * v).sum::<f64>().sqrt();
        let p = [raw[0] / l, raw[1] / l, raw[2] / l];
        let n = crate::boundaries::boundary_normal(p, cents[sa], cents[sb]).unwrap();
        let k = 0.01; // rad/My
        let omega = [
            k * (p[1] * n[2] - p[2] * n[1]),
            k * (p[2] * n[0] - p[0] * n[2]),
            k * (p[0] * n[1] - p[1] * n[0]),
        ];
        let pole = |v: [f64; 3]| {
            let mag = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            let (lat, lon) = unit_to_lat_lon([v[0] / mag, v[1] / mag, v[2] / mag]).unwrap();
            EulerPoleV0 {
                lat_deg: lat,
                lon_deg: lon,
                rate_deg_per_my: mag.to_degrees(),
            }
        };
        let zero = EulerPoleV0 {
            lat_deg: 90.,
            lon_deg: 0.,
            rate_deg_per_my: 0.,
        };
        let empty = InterfaceSamples::new();
        let (r, rec, g) = pair_normal_rate(
            &grid,
            &empty,
            &cents,
            cell,
            a,
            b,
            sa,
            sb,
            &pole(omega),
            &zero,
            2.,
        )
        .unwrap();
        assert_eq!(g, ArcContactGeometry::NewContactBridge);
        assert!(
            (r - k * grid.radius_m()).abs() < 1e-6 * k * grid.radius_m(),
            "{r}"
        );
        // b (id 1) < a (id 4): the loser's source is recorded first.
        assert_eq!(rec[0].0 as usize, sb);
        assert_eq!(rec[1].0 as usize, sa);
        let sheared = pole([
            omega[0] + k * n[0],
            omega[1] + k * n[1],
            omega[2] + k * n[2],
        ]);
        let (r2, _, _) = pair_normal_rate(
            &grid, &empty, &cents, cell, a, b, sa, sb, &sheared, &zero, 2.,
        )
        .unwrap();
        assert!((r2 - r).abs() < 1e-6 * r.abs());
        // A far step-start interface does not stand in for the local contact.
        let far = InterfaceSamples::from([((b, a), vec![(0usize, 1usize)])]);
        let (r3, rec3, g3) = pair_normal_rate(
            &grid,
            &far,
            &cents,
            cell,
            a,
            b,
            sa,
            sb,
            &pole(omega),
            &zero,
            2.,
        )
        .unwrap();
        assert_eq!((r3, rec3, g3), (r, rec, g));
        // Sources farther apart than the swept bound: typed refusal.
        let far_src = 8 * 36 + 20;
        match pair_normal_rate(
            &grid, &far, &cents, cell, a, b, sa, far_src, &zero, &zero, 2.,
        ) {
            Err(KinematicError::ArcContact {
                cell: c,
                winner,
                loser,
                why:
                    ArcContactRefusal::NonlocalContact {
                        nearest_interface_rad: Some(_),
                        source_separation_rad,
                        bound_rad,
                    },
            }) => {
                assert_eq!((c, winner, loser), (cell, a, b));
                assert!(source_separation_rad > bound_rad);
            }
            other => panic!("expected a typed nonlocal refusal, got {other:?}"),
        }
        match pair_normal_rate(&grid, &empty, &cents, cell, a, b, sa, sa, &zero, &zero, 2.) {
            Err(KinematicError::ArcContact {
                why: ArcContactRefusal::DegenerateContact { .. },
                ..
            }) => {}
            other => panic!("expected a typed degenerate refusal, got {other:?}"),
        }
    }

    #[test]
    fn nearest_cell_roundtrips_centroids_and_wraps() {
        let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
        let cents = centroids(&grid).unwrap();
        for (c, p) in cents.iter().enumerate() {
            assert_eq!(
                nearest_cell(&grid, *p).unwrap(),
                c,
                "centroid {c} did not map to itself"
            );
        }
        // A point just west of the antimeridian lands in the last column.
        let p = lat_lon_to_unit(0.0, -179.99).unwrap();
        let c = nearest_cell(&grid, p).unwrap();
        assert_eq!(c % grid.nlon(), 0);
        let p = lat_lon_to_unit(0.0, 179.99).unwrap();
        let c = nearest_cell(&grid, p).unwrap();
        assert_eq!(c % grid.nlon(), grid.nlon() - 1);
    }

    #[test]
    fn stationary_plates_carry_state_exactly_and_conserve_volume() {
        let (grid, mut plates, state, prior) = setup(5);
        for p in &mut plates {
            if let Some(e) = p.euler_pole.as_mut() {
                e.rate_deg_per_my = 0.0;
            }
        }
        let (next, ledger, thickening) = advance(&grid, &state, &plates, &prior, 2.0).unwrap();
        assert!(thickening.iter().all(|&t| t == 0.0));
        assert_eq!(next.plate, state.plate);
        assert_eq!(next.thickness_m, state.thickness_m);
        assert_eq!(next.kind, state.kind);
        assert_eq!(ledger.ridge_cells, 0);
        assert_eq!(ledger.subduction_cells, 0);
        assert_eq!(ledger.collision_cells, 0);
        assert!(
            ledger.residual_m3().abs() < 1.0e-6 * ledger.volume_before_m3,
            "residual {}",
            ledger.residual_m3()
        );
        assert!(next
            .boundary_class
            .iter()
            .all(|&c| c == BoundaryClassCodeV0::Interior.code()
                || c == BoundaryClassCodeV0::Transform.code()));
    }

    #[test]
    fn moving_plates_open_ridges_and_consume_crust_deterministically() {
        let (grid, plates, state, prior) = setup(5);
        // Ten steps of 2 My at up to 1 deg/My on 10-degree cells: plates
        // move up to 2 cells in total, enough to part and meet.
        let run = |state: &CrustState| {
            let mut s = state.clone();
            let mut total = AdvanceLedger::default();
            let mut last = None;
            for _ in 0..10 {
                let (next, ledger, thick) = advance(&grid, &s, &plates, &prior, 2.0).unwrap();
                total.ridge_cells += ledger.ridge_cells;
                total.subduction_cells += ledger.subduction_cells;
                total.collision_cells += ledger.collision_cells;
                total.ridge_created_m3 += ledger.ridge_created_m3;
                total.subducted_m3 += ledger.subducted_m3;
                total.arc_created_m3 += ledger.arc_created_m3;
                total.delaminated_m3 += ledger.delaminated_m3;
                last = Some((ledger, thick));
                s = next;
            }
            (s, total, last.unwrap())
        };
        let (a, la, (last_a, ta)) = run(&state);
        let (b, lb, (last_b, tb)) = run(&state);
        assert_eq!(a, b, "same inputs must give bit-identical state");
        assert_eq!(la, lb);
        assert_eq!(last_a, last_b);
        assert_eq!(ta, tb);
        assert!(
            la.ridge_cells > 0,
            "five moving plates must part somewhere in 20 My"
        );
        assert!(
            la.subduction_cells + la.collision_cells > 0,
            "and meet somewhere in 20 My"
        );
        // Ridge-born cells of the LAST step are oceanic, age zero, 7 km.
        let born: Vec<usize> = (0..a.cells())
            .filter(|&c| a.birth_t_my[c] == 20.0)
            .collect();
        assert_eq!(born.len() as u64, last_a.ridge_cells);
        for &c in &born {
            assert_eq!(a.kind[c], CrustKindCodeV0::Oceanic);
            assert_eq!(a.thickness_m[c], prior.record.oceanic_crust_thickness_m);
        }
        // Ledger terms are consistent signs and the residual is reported.
        assert!(la.ridge_created_m3 > 0.0);
        assert!(la.subducted_m3 >= 0.0 && la.arc_created_m3 >= 0.0 && la.delaminated_m3 >= 0.0);
        let rel = last_a.residual_m3().abs() / last_a.volume_before_m3;
        assert!(
            rel < 0.05,
            "advection residual {rel} exceeds 5% in one step (carry broken?)"
        );
        // No spurious polar ridges: pole rows must not be reborn every step
        // (the metric search window covers whole rows there).
        let nlon = grid.nlon();
        let polar_born = born
            .iter()
            .filter(|&&c| c / nlon == 0 || c / nlon == grid.nlat() - 1)
            .count();
        assert!(
            polar_born <= born.len() / 2,
            "{polar_born} of {} births on the pole rows",
            born.len()
        );
        // Some cell got classified as a real boundary.
        assert!(a
            .boundary_class
            .iter()
            .any(|&c| c == BoundaryClassCodeV0::Ridge.code()));
        assert!(a
            .boundary_class
            .iter()
            .any(|&c| c == BoundaryClassCodeV0::Subduction.code()
                || c == BoundaryClassCodeV0::Collision.code()));
        // Normal rates carry physical units (m/My): below a few hundred km/My.
        let max_rate = a
            .normal_rate_m_per_my
            .iter()
            .fold(0.0f64, |m, r| m.max(r.abs()));
        assert!(
            max_rate > 0.0 && max_rate < 5.0e5,
            "max |normal rate| {max_rate} m/My"
        );
    }

    #[test]
    fn sediment_ledger_closes_when_stationary_and_names_collision_moves() {
        // Stationary: no subduction, no stacking, remap drift ~ 0.
        let (grid, mut plates, mut state, prior) = setup(5);
        for c in 0..state.cells() {
            state.sediment_m[c] = 100.0 + (c % 7) as f64;
        }
        let moving = plates.clone();
        for p in &mut plates {
            if let Some(e) = p.euler_pole.as_mut() {
                e.rate_deg_per_my = 0.0;
            }
        }
        let (_, l, _) = advance(&grid, &state, &plates, &prior, 2.0).unwrap();
        assert_eq!(l.sediment_subducted_m3, 0.0);
        assert_eq!(l.sediment_stacked_m3, 0.0);
        assert!(l.sediment_residual_m3().abs() <= 1e-9 * l.sediment_before_m3);
        // Moving: subducted sediment is named; drift is reported, finite.
        let (_, m, _) = advance(&grid, &state, &moving, &prior, 2.0).unwrap();
        assert!(m.sediment_before_m3 > 0.0 && m.sediment_after_m3 > 0.0);
        assert!(m.sediment_subducted_m3 >= 0.0 && m.sediment_stacked_m3 >= 0.0);
        assert!(m.sediment_residual_m3().is_finite());
    }

    #[test]
    fn collision_caps_thickness_and_accounts_excess() {
        let (grid, mut plates, mut state, prior) = setup(2);
        // Make everything continental and drive the plates together.
        for c in 0..state.cells() {
            state.kind[c] = CrustKindCodeV0::Continental;
            state.thickness_m[c] = 60_000.0;
            state.birth_t_my[c] = -2000.0;
        }
        for (k, p) in plates.iter_mut().enumerate() {
            if let Some(e) = p.euler_pole.as_mut() {
                e.lat_deg = 90.0;
                e.lon_deg = 0.0;
                e.rate_deg_per_my = if k == 0 { 3.0 } else { -3.0 };
            }
        }
        let (next, ledger, thickening) = advance(&grid, &state, &plates, &prior, 2.0).unwrap();
        assert!(
            thickening.iter().any(|&t| t > 0.0),
            "collision must thicken somewhere"
        );
        assert!(ledger.collision_cells > 0, "opposed rotations must collide");
        let max_h = next.thickness_m.iter().cloned().fold(0.0f64, f64::max);
        assert!(
            max_h <= prior.record.max_crust_thickness_m + 1.0e-6,
            "cap breached: {max_h}"
        );
        // 60 + 60 > 70: excess must have been passed or delaminated.
        assert!(
            ledger.delaminated_m3 > 0.0
                || next.shortening_added_m.iter().filter(|&&s| s > 0.0).count()
                    > ledger.collision_cells as usize
        );
        assert_eq!(ledger.subduction_cells, 0);
    }

    #[test]
    fn slow_motion_accumulates_instead_of_stalling() {
        // One plate rotating about +z at 0.2 cells per step on a 10-degree
        // grid: after 10 steps a continental stripe must have moved ~2
        // columns east, not 0 (the quantization stall this design fixes).
        let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
        let n = grid.cells();
        let prior = UpliftPrior::preview_v0();
        let pole = EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: 1.0, // 1 deg/My x 2 My = 2 deg = 0.2 cells per step
        };
        let plates = vec![PlateRecordV0 {
            plate_id: PlateId(0),
            parent_plate: None,
            birth_t_my: 0.0,
            euler_pole: Some(pole),
            retirement: None,
        }];
        let mut state = CrustState {
            t_my: 0.0,
            plate: vec![PlateId(0); n],
            kind: vec![CrustKindCodeV0::Oceanic; n],
            birth_t_my: vec![-50.0; n],
            thickness_m: vec![7000.0; n],
            sediment_m: vec![0.0; n],
            eroded_total_m: vec![0.0; n],
            arc_added_m: vec![0.0; n],
            shortening_added_m: vec![0.0; n],
            boundary_class: vec![0; n],
            normal_rate_m_per_my: vec![0.0; n],
            material_pos: centroids(&grid).unwrap(),
        };
        // Stripe at columns 10..=12 on the equatorial rows.
        let stripe = |j: usize| (10..=12).contains(&j);
        for i in 8..10 {
            for j in 0..36 {
                if stripe(j) {
                    let c = i * 36 + j;
                    state.kind[c] = CrustKindCodeV0::Continental;
                    state.thickness_m[c] = 35_000.0;
                }
            }
        }
        let centroid_of = |kind: &[CrustKindCodeV0], row: usize| -> f64 {
            let mut sum = 0.0;
            let mut count = 0.0;
            for j in 0..36 {
                if kind[row * 36 + j] == CrustKindCodeV0::Continental {
                    sum += j as f64;
                    count += 1.0;
                }
            }
            sum / count
        };
        let start = centroid_of(&state.kind, 8);
        for _ in 0..10 {
            let (next, ledger, _) = advance(&grid, &state, &plates, &prior, 2.0).unwrap();
            assert_eq!(ledger.ridge_cells, 0, "a lone plate never opens a ridge");
            state = next;
        }
        let end = centroid_of(&state.kind, 8);
        // 10 steps x 0.2 cells = 2.0 columns east (positive rate = eastward).
        assert!(
            (end - start - 2.0).abs() <= 0.5,
            "stripe moved {} columns; expected ~2 (stall would give 0)",
            end - start
        );
        // Column phases stay bounded: material never drifts more than a cell.
        let cents = centroids(&grid).unwrap();
        let dlon = (10.0f64).to_radians();
        for c in 0..n {
            let d = great_circle_rad(state.material_pos[c], cents[c]).unwrap();
            assert!(d <= 1.01 * dlon, "cell {c} phase {d} rad exceeds one cell");
        }
    }

    #[test]
    fn classification_uses_relative_motion_not_plate_ids() {
        let (grid, mut plates, state, _prior) = setup(3);
        // Identical rotation for all plates: no relative motion anywhere,
        // so every boundary cell is Transform (never Ridge/Subduction).
        let common = plates[0].euler_pole.unwrap();
        for p in &mut plates {
            p.euler_pole = Some(common);
        }
        let (class, rate, footprint) = classify_cells(&grid, &state, &plates).unwrap();
        assert!(
            class.iter().any(|&c| c != 0),
            "plates touch, so boundary cells exist"
        );
        assert!(class
            .iter()
            .all(|&c| c == 0 || c == BoundaryClassCodeV0::Transform.code()));
        assert!(rate.iter().all(|r| r.abs() <= 1.0e-9));
        assert!(footprint.keys().all(|k| k.2 == BoundaryKindV0::Transform));
    }
}

/// SEED-TRANSFORMS (lane TF, 2026-10-13): `smoothed-normal-transform-v1`.
/// Two oceanic plates split by a tilted great circle (a staircase on the
/// grid). Expected classes come from the construction: rotation about the
/// circle's own pole moves plate B exactly parallel to the boundary (pure
/// transform); rotation about a pole ON the circle moves it across the
/// boundary (pure opening or closing).
#[cfg(test)]
mod transform_tests {
    use super::*;
    use crate::boundaries::TransformParams;

    const BOUNDARY_POLE: (f64, f64) = (60.0, 20.0);

    fn world(pole: (f64, f64), rate: f64) -> (GridSpec, Vec<PlateRecordV0>, CrustState) {
        let grid = GridSpec::new(144, 288, 6.371e6).unwrap();
        let n = grid.cells();
        let m = lat_lon_to_unit(BOUNDARY_POLE.0, BOUNDARY_POLE.1).unwrap();
        let cents = centroids(&grid).unwrap();
        let mut st = CrustState {
            t_my: 0.0,
            plate: vec![PlateId(0); n],
            kind: vec![CrustKindCodeV0::Oceanic; n],
            birth_t_my: vec![-50.0; n],
            thickness_m: vec![7_000.0; n],
            sediment_m: vec![0.0; n],
            eroded_total_m: vec![0.0; n],
            arc_added_m: vec![0.0; n],
            shortening_added_m: vec![0.0; n],
            boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
            normal_rate_m_per_my: vec![0.0; n],
            material_pos: cents.clone(),
        };
        for c in 0..n {
            let d = cents[c][0] * m[0] + cents[c][1] * m[1] + cents[c][2] * m[2];
            if d >= 0.0 {
                st.plate[c] = PlateId(1);
            }
        }
        let rec = |id: u32, lat: f64, lon: f64, r: f64| PlateRecordV0 {
            plate_id: PlateId(id),
            parent_plate: None,
            birth_t_my: 0.0,
            euler_pole: Some(EulerPoleV0 { lat_deg: lat, lon_deg: lon, rate_deg_per_my: r }),
            retirement: None,
        };
        let plates = vec![rec(0, 90.0, 0.0, 0.0), rec(1, pole.0, pole.1, rate)];
        (grid, plates, st)
    }

    fn shares(cls: &[u8]) -> (usize, usize, usize, usize) {
        let k = |c: BoundaryClassCodeV0| cls.iter().filter(|&&x| x == c.code()).count();
        (
            k(BoundaryClassCodeV0::Transform),
            k(BoundaryClassCodeV0::Ridge),
            k(BoundaryClassCodeV0::Subduction),
            k(BoundaryClassCodeV0::Collision),
        )
    }

    #[test]
    fn smoothed_normal_finds_the_transform_the_staircase_misses() {
        let (grid, plates, st) = world(BOUNDARY_POLE, 1.0);
        let tp = TransformParams::v1();
        let (c2, r2, _, fb) = classify_cells_v2(&grid, &st, &plates, &tp).unwrap();
        let (t, ri, su, co) = shares(&c2);
        let b = t + ri + su + co;
        assert!(b > 100, "boundary cells {b}");
        assert!(t as f64 >= 0.95 * b as f64, "v2 transform {t} of {b} (ridge {ri}, subd {su}), fallbacks {fb}");
        // Normal rates are small against the 1 deg/My shear (|v| ~ 111 km/My).
        let mut vn: Vec<f64> = (0..c2.len())
            .filter(|&c| c2[c] != BoundaryClassCodeV0::Interior.code())
            .map(|c| r2[c].abs())
            .collect();
        vn.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p95 = vn[(0.95 * (vn.len() - 1) as f64) as usize];
        assert!(p95 < 0.3 * 111_000.0, "p95 |v_n| {p95} m/My");
        // The legacy staircase classifier labels much of it ridge/subduction.
        let (c1, _, _) = classify_cells(&grid, &st, &plates).unwrap();
        let (t1, r1, s1, _) = shares(&c1);
        assert!((r1 + s1) as f64 >= 0.3 * (t1 + r1 + s1) as f64, "legacy t {t1} r {r1} s {s1}");
        // Deterministic.
        assert_eq!(classify_cells_v2(&grid, &st, &plates, &tp).unwrap().0, c2);
        assert_eq!(classify_cells_law(&grid, &st, &plates, None).unwrap().0, c1);
    }

    #[test]
    fn opening_and_closing_keep_their_classes() {
        // A pole ON the boundary circle: motion across the boundary, opening
        // on one half of the circle and closing on the other. Expected
        // normal rate at a boundary cell, from the construction (true a->b
        // normal = +/- the circle pole m): v_n = -(w_b x r) . m.
        let tp = TransformParams::v1();
        let (grid, plates, st) = world((30.0, BOUNDARY_POLE.1 - 180.0), 1.0);
        let (c2, _, _, _) = classify_cells_v2(&grid, &st, &plates, &tp).unwrap();
        let m = lat_lon_to_unit(BOUNDARY_POLE.0, BOUNDARY_POLE.1).unwrap();
        let wb = crate::boundaries::angular_velocity(&plates[1].euler_pole.unwrap()).unwrap();
        let cents = centroids(&grid).unwrap();
        let (mut checked, mut right) = (0usize, 0usize);
        for c in 0..c2.len() {
            if c2[c] == BoundaryClassCodeV0::Interior.code() {
                continue;
            }
            let v = cross3(wb, cents[c]);
            let vn = -(v[0] * m[0] + v[1] * m[1] + v[2] * m[2]);
            let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            // Skip cells near the rotation pole and near the obliquity edge.
            if speed < 1.0e-3 || vn.abs() < 0.7 * speed {
                continue;
            }
            checked += 1;
            let want = if vn > 0.0 { BoundaryClassCodeV0::Subduction } else { BoundaryClassCodeV0::Ridge };
            right += (c2[c] == want.code()) as usize;
        }
        assert!(checked > 100, "checked {checked}");
        assert!(right as f64 >= 0.95 * checked as f64, "{right} of {checked}");
        let (t, ri, su, _) = shares(&c2);
        assert!(ri > 50 && su > 50 && (t as f64) < 0.2 * (t + ri + su) as f64, "t {t} r {ri} s {su}");
    }

    /// PLATE-POLES `rigid-group-contacts-v1`: a cap split into two plate ids
    /// with the same pole (a welded suture group) has zero relative motion
    /// across the split. Legacy labels the split TRANSFORM; with the key the
    /// split is interior and every other cell keeps its class.
    #[test]
    fn rigid_group_contacts_are_interior() {
        let (grid, mut plates, mut st) = world((30.0, BOUNDARY_POLE.1 - 180.0), 1.0);
        let cents = centroids(&grid).unwrap();
        for c in 0..st.cells() {
            if st.plate[c] == PlateId(1) && cents[c][1] > 0.0 {
                st.plate[c] = PlateId(2);
            }
        }
        let mut twin = plates[1].clone();
        twin.plate_id = PlateId(2);
        plates.push(twin);
        let off = TransformParams::v1();
        let on = TransformParams { rigid_group_interior: true, ..off };
        let (c0, r0, f0, _) = classify_cells_v2(&grid, &st, &plates, &off).unwrap();
        let (c1, r1, f1, _) = classify_cells_v2(&grid, &st, &plates, &on).unwrap();
        // Cells touching only the twin plate across the split.
        let mut split_only = 0usize;
        let mut nbs = [0usize; 8];
        for c in 0..st.cells() {
            let nn = neighbor_ids(grid.nlat(), grid.nlon(), c, &mut nbs);
            let others: Vec<PlateId> = nbs[..nn].iter().map(|&k| st.plate[k]).filter(|&p| p != st.plate[c]).collect();
            if others.is_empty() {
                assert_eq!(c1[c], c0[c]);
                continue;
            }
            let welded = |p: PlateId| p != PlateId(0) && st.plate[c] != PlateId(0);
            if others.iter().all(|&p| welded(p)) {
                split_only += 1;
                assert_eq!(c0[c], BoundaryClassCodeV0::Transform.code(), "legacy labels the weld transform");
                assert_eq!(c1[c], BoundaryClassCodeV0::Interior.code());
                assert_eq!(r1[c], 0.0);
            } else {
                assert_eq!(c1[c], c0[c], "cell {c} touching a moving plate keeps its class");
            }
        }
        assert!(split_only > 50, "split cells {split_only}");
        assert!(f0.keys().any(|k| k.0 == PlateId(1) && k.1 == PlateId(2)));
        assert!(!f1.keys().any(|k| k.0 == PlateId(1) && k.1 == PlateId(2)), "no boundary record for the weld");
        // Off is bit for bit the v1 classifier.
        assert_eq!(classify_cells_law(&grid, &st, &plates, Some(&off)).unwrap().0, c0);
    }

    fn run(pole: (f64, f64), rate: f64, law: Option<&TransformParams>, steps: usize) -> (CrustState, Vec<AdvanceLedger>) {
        let (grid, plates, mut st) = world(pole, rate);
        let prior = UpliftPrior::preview_v0();
        let mut logs = Vec::new();
        for _ in 0..steps {
            let (o, l, _, _) = advance_impl(
                &grid,
                &st,
                &plates,
                &prior,
                2.0,
                None,
                CarryLaw::default(),
                crate::orogen::CollisionSedimentLaw::default(),
                None,
                None,
                RidgeAccretion::PreviousOccupantV0,
                None,
                law,
            )
            .unwrap();
            st = o;
            logs.push(l);
        }
        (st, logs)
    }

    #[test]
    fn transforms_neither_create_nor_subduct_crust() {
        let tp = TransformParams::v1();
        let (legacy, ll) = run(BOUNDARY_POLE, 1.0, None, 5);
        let (v2, lv) = run(BOUNDARY_POLE, 1.0, Some(&tp), 5);
        let young = |s: &CrustState| s.birth_t_my.iter().filter(|&&b| b > 0.0).count();
        let births = |l: &[AdvanceLedger]| l.iter().map(|x| x.ridge_cells).sum::<u64>();
        let subd = |l: &[AdvanceLedger]| l.iter().map(|x| x.subduction_cells).sum::<u64>();
        let filled: u64 = lv.iter().map(|x| x.transform.filled_transform_cells).sum();
        assert!(births(&ll) > 20 && young(&legacy) > 20, "legacy births {} young {}", births(&ll), young(&legacy));
        assert!(subd(&ll) > 20, "legacy subduction overlaps {}", subd(&ll));
        assert!(filled > 20, "filled {filled}");
        assert!(births(&lv) * 10 <= births(&ll), "v2 births {} vs legacy {}", births(&lv), births(&ll));
        assert!(young(&v2) * 10 <= young(&legacy), "v2 young {} vs legacy {}", young(&v2), young(&legacy));
        assert_eq!(subd(&lv), 0, "no subduction on a pure transform");
        // Fills and non-convergent overlaps are named remap terms, not
        // ridge creation or subduction.
        assert!(lv.iter().all(|l| l.subducted_m3 == 0.0 && l.transform.filled_m3 > 0.0));
        assert!(lv.iter().map(|l| l.transform.nonconvergent_overlaps).sum::<u64>() > 0);
        // Deterministic.
        let (again, _) = run(BOUNDARY_POLE, 1.0, Some(&tp), 5);
        assert_eq!(again.birth_t_my, v2.birth_t_my);
        assert_eq!(again.plate, v2.plate);
    }

    #[test]
    fn ridges_still_open_under_v2() {
        let tp = TransformParams::v1();
        let (_, ll) = run((30.0, BOUNDARY_POLE.1 - 180.0), 1.0, None, 3);
        let (_, lv) = run((30.0, BOUNDARY_POLE.1 - 180.0), 1.0, Some(&tp), 3);
        let births = |l: &[AdvanceLedger]| l.iter().map(|x| x.ridge_cells).sum::<u64>();
        assert!(births(&ll) > 20);
        assert!(births(&lv) as f64 >= 0.8 * births(&ll) as f64, "v2 {} legacy {}", births(&lv), births(&ll));
    }
}

/// SEED-PERF (lane SP, 2026-10-13): the fast paths against the legacy
/// forms on adversarial inputs (jittered and duplicated material
/// positions, poles, exact ties, many plates).
#[cfg(test)]
mod speed_tests {
    use super::*;

    struct Lcg(u64);
    impl Lcg {
        fn f(&mut self) -> f64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
        }
    }

    fn norm(v: [f64; 3]) -> [f64; 3] {
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / n, v[1] / n, v[2] / n]
    }

    /// Jittered material positions (some exact duplicates of a neighbour's,
    /// for exact ties) and a blocky plate field with `np` plates.
    fn world(grid: &GridSpec, np: u32, seed: u64) -> (Vec<[f64; 3]>, Vec<PlateId>) {
        let mut r = Lcg(seed);
        let cents = centroids(grid).unwrap();
        let cell = std::f64::consts::PI / grid.nlat() as f64;
        let mut mp: Vec<[f64; 3]> = cents
            .iter()
            .map(|c| norm([c[0] + (r.f() - 0.5) * cell, c[1] + (r.f() - 0.5) * cell, c[2] + (r.f() - 0.5) * cell]))
            .collect();
        for c in 1..mp.len() {
            if r.f() < 0.05 {
                mp[c] = mp[c - 1];
            }
        }
        let seeds: Vec<[f64; 3]> = (0..np).map(|_| norm([r.f() - 0.5, r.f() - 0.5, r.f() - 0.5])).collect();
        let plate = cents
            .iter()
            .map(|c| {
                let k = (0..seeds.len())
                    .max_by(|&a, &b| {
                        let da = c[0] * seeds[a][0] + c[1] * seeds[a][1] + c[2] * seeds[a][2];
                        let db = c[0] * seeds[b][0] + c[1] * seeds[b][1] + c[2] * seeds[b][2];
                        da.total_cmp(&db)
                    })
                    .unwrap();
                PlateId(10 + 3 * k as u32)
            })
            .collect();
        (mp, plate)
    }

    #[test]
    fn material_search_is_the_legacy_claim() {
        for (nlat, nlon, np, seed) in [(18, 36, 5, 1u64), (36, 72, 12, 2), (90, 180, 24, 3), (8, 16, 3, 4)] {
            let grid = GridSpec::new(nlat, nlon, 6.371e6).unwrap();
            let (mp, plate) = world(&grid, np, seed);
            let search = MaterialSearch::new(&grid, &mp, &plate).unwrap();
            assert!(search.tol.is_some());
            let mut r = Lcg(seed + 100);
            let ids: std::collections::BTreeSet<PlateId> = plate.iter().cloned().collect();
            let mut checked = 0;
            for t in 0..4000 {
                // Random points, points at material positions (exact ties) and near the poles.
                let q = match t % 4 {
                    0 => mp[(r.f() * mp.len() as f64) as usize],
                    1 => norm([r.f() - 0.5, r.f() - 0.5, (r.f() - 0.5) * 0.02 + if r.f() < 0.5 { 1.0 } else { -1.0 }]),
                    _ => norm([r.f() - 0.5, r.f() - 0.5, r.f() - 0.5]),
                };
                let s = nearest_material(&grid, &mp, q).unwrap();
                for id in &ids {
                    let legacy = (plate[s] == *id).then_some(s);
                    assert_eq!(search.claim(q, *id).unwrap(), legacy, "grid {nlat}x{nlon} q {q:?} plate {id:?}");
                    checked += 1;
                }
            }
            assert!(checked > 0);
        }
    }

    #[test]
    fn reach_masks_cover_every_legacy_claim() {
        for (nlat, nlon, np, seed, rate) in [(36, 72, 6, 5u64, 3.0), (90, 180, 20, 6, 1.5), (144, 288, 30, 7, 0.8)] {
            let grid = GridSpec::new(nlat, nlon, 6.371e6).unwrap();
            let (mp, plate) = world(&grid, np, seed);
            let mut r = Lcg(seed + 7);
            let mut ids: Vec<PlateId> = plate.clone();
            ids.sort();
            ids.dedup();
            let live: Vec<(PlateId, EulerPoleV0)> = ids
                .iter()
                .map(|id| {
                    (
                        *id,
                        EulerPoleV0 {
                            lat_deg: r.f() * 180.0 - 90.0,
                            lon_deg: r.f() * 360.0 - 180.0,
                            rate_deg_per_my: (r.f() - 0.5) * 2.0 * rate,
                        },
                    )
                })
                .collect();
            let dt = 2.0;
            let masks = plate_reach_masks(&grid, &plate, &ids, claim_reach_radius(&grid, &live, dt)).unwrap();
            let cents = centroids(&grid).unwrap();
            let mut claims = 0;
            for c in 0..grid.cells() {
                for (k, (id, e)) in live.iter().enumerate() {
                    let back = EulerPoleV0 { rate_deg_per_my: -e.rate_deg_per_my, ..*e };
                    let q = rotate_by_euler(cents[c], &back, dt).unwrap();
                    let s = nearest_material(&grid, &mp, q).unwrap();
                    if plate[s] == *id {
                        claims += 1;
                        assert!(masks[c] & (1u64 << k) != 0, "cell {c} plate {id:?} claims but is masked out");
                    }
                }
            }
            assert!(claims > grid.cells() / 2);
            // The filter must actually filter (not a vacuous all-ones mask).
            let skipped: u32 = masks.iter().map(|m| live.len() as u32 - m.count_ones()).sum();
            assert!(skipped > 0, "{nlat}x{nlon}: only {skipped} pairs skipped");
        }
    }

    #[test]
    fn parallel_map_is_the_serial_map() {
        for t in [1usize, 2, 3, 16, 64] {
            par::set_threads(t);
            let v = par::map_indexed(100_003, |i| (i as f64).sqrt().to_bits() ^ (i as u64));
            let w: Vec<u64> = (0..100_003).map(|i| (i as f64).sqrt().to_bits() ^ (i as u64)).collect();
            assert_eq!(v, w);
            let h = par::map_heavy(37, |i| vec![i; i]);
            assert_eq!(h, (0..37).map(|i| vec![i; i]).collect::<Vec<_>>());
            let mut f = vec![0u64; 50_001];
            par::fill_indexed(&mut f, |i| (i * 7) as u64);
            assert!(f.iter().enumerate().all(|(i, x)| *x == (i * 7) as u64));
        }
        par::set_threads(0);
    }
}
