//! SEED-PLATE-REORG (lane PR, 2026-10-13): `plate-reorg-scripted-v1`, a
//! scripted Wilson (supercontinent) cycle. Pre-registered in
//! `unpublished-design-record`.
//!
//! Why: with fixed Euler poles (`plates.rs`) and permanent sutures
//! (`conserve::suture_pass`) every world assembles one supercontinent and
//! stops building mountains once the sutures lock (lanes SC, CV, IS). This is
//! scripted kinematics, not a mantle force balance (Stokes stays frozen):
//!
//! - **Pole re-solve** at seeded times (gaps uniform in [min, max] My from
//!   one stream, so the schedule does not depend on the run length): every
//!   suture group gets a new pole in ascending root order; the axis mixes a
//!   random axis with the group's slab torque (sum over its subducting edges
//!   of x cross the direction to the trench), and the rate grows with the
//!   subducting share of its boundary (Forsyth and Uyeda 1975, recalled).
//!   Breakups and initiations younger than the persistence time keep their
//!   relative rotation.
//! - **Breakup**: a group holding a supercontinent share of the continental
//!   crust for the tenure time splits along a MIN CUT (Dinic max flow) over
//!   its cells between the cells beyond the 25 % and 75 % area quantiles of
//!   its continental cells along a seeded axis. Capacity = edge length x
//!   strength (anisotropic 8-neighbour cut metric, Boykov and Kolmogorov
//!   2003), so the cut follows weak zones: crust far from the group's
//!   median continental thickness (rifted/thin crust, belts and the
//!   irregular-cratons-v1 suture roots), old shortening, and edges between
//!   different member plates of the group (old sutures). The new side opens
//!   about p x n at the rift opening rate; `stretched-rift-v1` and
//!   `shoulder-uplift-thermal-v1` then make the rift unchanged.
//! - **Plume-head splits** use the same cut, weakened within the plume head
//!   radius, only for groups holding a minimum continental share.
//! - **Subduction initiation**: while too few groups exceed the counting
//!   area, the group with the most ocean splits along its passive margin
//!   (source = its continent, sink = its ocean beyond 600 km) or along a
//!   seeded axis, and its ocean side converges.
//! - **Merge/retire**: secondary pieces of a plate id fold into a neighbour
//!   (small, oceanic) or become plates of the same group; above the maximum
//!   count the smallest oceanic group folds into its neighbour.
//! - **Plume renewal**: plume onsets are a renewal process from one stream,
//!   tails last a fixed lifetime (duration-independent schedule).
//!
//! Determinism: SplitMix64 streams tagged per event index, ascending
//! ids/cells, BTreeMap only, no threads.

use crate::boundaries::{angular_velocity, normal_rate_at};
use crate::conserve::{relabel_parcels, ParcelLayer, SutureState};
use crate::hotspots::Plume;
use crate::kinematic::{cell_areas, centroids, nearest_cell, CrustState};
use crate::plates::{seed_from_path_tagged, SplitMix64};
use crate::rotation::unit_to_lat_lon;
use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0, PlateRetirementV0};
use pf_contracts::params::SeedPath;
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;
use std::collections::BTreeMap;
use std::f64::consts::PI;

pub const LAW: &str = "plate-reorg-scripted-v1";
/// SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13): the second law version
/// (pre-registered in `unpublished-design-record`).
pub const LAW_V2: &str = "plate-reorg-scripted-v2";
/// SEED-WORLD-COMPOSE (lane CP).
pub const LAW_V3: &str = "plate-reorg-scripted-v3";
/// SEED-TRANSFORMS (lane TF, 2026-10-13).
pub const LAW_V4: &str = "plate-reorg-scripted-v4";

/// SEED-OLD-OCEAN (lane OR, 2026-10-13).
pub const LAW_V5: &str = "plate-reorg-scripted-v5";

/// SEED-OLD-OCEAN (lane OR) `plate-reorg-scripted-v5`: v4 plus
/// passive-margin subduction initiation driven by ocean age (lane VL
/// defect 3), pre-registered in `unpublished-design-record`.
/// At each scheduled event, while fewer than the maximum groups, the
/// OLDEST passive-margin ocean cell (oceanic, in the same suture group as
/// a continental 8-neighbour, touching no other group, at least `age_my`
/// old, farther than `spacing_m` from any subduction-class cell or trench
/// started earlier in the pass) seeds a trench: the segment is the
/// margin cells reachable from the seed along the margin within half of
/// `segment_length_m` (at least `MARGIN_SEGMENT_MIN_AGE_MY` old, also
/// away from existing trenches); side B is the group's ocean whose
/// nearest continental cell is on the segment's coast (the margin's
/// frontal ocean, connected to the seed), so the new boundary runs along
/// the segment's coast and leaves it along the equidistance lines to the
/// rest of the coast. B converges on the continent at the subduction
/// initiation rate about the pole 90 degrees from the segment centre.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MarginInitiation {
    pub age_my: f64,
    pub segment_length_m: f64,
    pub spacing_m: f64,
}

/// Fixed by the v5 law version (lane OR PREREG.txt).
pub const MARGIN_SEGMENT_MIN_AGE_MY: f64 = 100.0;
/// Lane OR PREREG revision 1 (2026-10-13): 2 -> 3 (the cap was hit at
/// almost every event in r0).
pub const MARGIN_MAX_PER_EVENT: usize = 3;
/// Revision 1: seeds are ranked by the ocean at least `age_my` old of
/// their group within this distance (largest first; then oldest, lowest
/// cell), so large old margins go before tiny trapped pockets.
pub const MARGIN_RANK_RADIUS_M: f64 = 1_000_000.0;
/// Revision 1: the ridge-aware re-solve also counts boundary edges whose
/// two cells are both ocean younger than this (My).
pub const RESOLVE_YOUNG_EDGE_MY: f64 = 30.0;
/// Revision 1: sweeps of the ridge-aware axis choice (the second re-scores
/// every group against its neighbours' new poles).
pub const RESOLVE_SWEEPS: usize = 2;
pub const MARGIN_MAX_ATTEMPTS: usize = 6;
/// Side B must hold at least this share of the sphere.
pub const MARGIN_MIN_SIDE_B_FRACTION: f64 = 0.001;
/// v5 ridge-aware re-solve: rotation axes drawn per group (the first is
/// the v1-v4 draw); the one closing the least ridge length against the
/// neighbours' poles (already re-solved groups: new pole) is kept
/// (lane OR, the cheap fix for VL defect 2's young-slab convergent edges).
pub const RESOLVE_AXIS_CANDIDATES: usize = 8; // revision 1: 4 -> 8

impl MarginInitiation {
    /// Pre-registered values (`unpublished-design-record`).
    pub fn v1() -> Self {
        Self { age_my: 160.0, segment_length_m: 2_000_000.0, spacing_m: 1_000_000.0 }
    }
    pub fn validate(&self) -> Result<(), String> {
        if !(self.age_my.is_finite() && self.age_my > 0.0) {
            return Err(format!("{LAW_V5}: margin initiation age must be finite and positive, got {}", self.age_my));
        }
        if !(self.segment_length_m.is_finite() && self.segment_length_m > 0.0) {
            return Err(format!("{LAW_V5}: margin segment length must be finite and positive, got {}", self.segment_length_m));
        }
        if !(self.spacing_m.is_finite() && self.spacing_m >= 0.0) {
            return Err(format!("{LAW_V5}: margin trench spacing must be finite and >= 0, got {}", self.spacing_m));
        }
        Ok(())
    }
}

/// SEED-TRANSFORMS (lane TF) `plate-reorg-scripted-v4`: v3 with a slab-age
/// and a spacing limit on the AGE-triggered initiation only (the v1/v2
/// "too few groups" trigger is unchanged). A candidate cut is rejected
/// (next axis tried) when, by length along the new boundary (4-neighbour
/// edges between the two sides inside the group), less than
/// `old_slab_fraction_min` runs where the older flanking cell is ocean at
/// least `slab_age_min_my` old, or more than `INIT_V4_NEAR_FRACTION_MAX` (half)
/// runs within `trench_spacing_m` of an existing subduction-class cell or
/// of a boundary started earlier in the same pass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InitiationLimits {
    pub slab_age_min_my: f64,
    pub old_slab_fraction_min: f64,
    pub trench_spacing_m: f64,
}

/// Fixed by the v4 law version (lane TF PREREG revision 1, 2026-10-13: 0.25
/// -> 0.5; at 0.25 the spacing test rejected 37-45 age-triggered cuts per
/// 1000 My, most with 26-50% of their length near existing subduction, and
/// ocean-age p95 rose to 227-241 My).
pub const INIT_V4_NEAR_FRACTION_MAX: f64 = 0.5;

impl InitiationLimits {
    /// Pre-registered values (`unpublished-design-record`).
    pub fn v1() -> Self {
        Self { slab_age_min_my: 100.0, old_slab_fraction_min: 0.5, trench_spacing_m: 1_000_000.0 }
    }
    pub fn validate(&self) -> Result<(), String> {
        if !(self.slab_age_min_my.is_finite() && self.slab_age_min_my >= 0.0) {
            return Err(format!("{LAW_V4}: slab_age_min_my must be finite and >= 0, got {}", self.slab_age_min_my));
        }
        if !(self.old_slab_fraction_min.is_finite() && (0.0..=1.0).contains(&self.old_slab_fraction_min)) {
            return Err(format!("{LAW_V4}: old_slab_fraction_min must be in [0, 1], got {}", self.old_slab_fraction_min));
        }
        if !(self.trench_spacing_m.is_finite() && self.trench_spacing_m >= 0.0) {
            return Err(format!("{LAW_V4}: trench_spacing_m must be finite and >= 0, got {}", self.trench_spacing_m));
        }
        Ok(())
    }
}

/// Which version of the scripted reorganization a run uses. v1 is
/// `reorg_pass` / `plume_split` unchanged; v2 changes only the breakup and
/// plume-split cut (shortening/seam/plume weakness, continental cut domain
/// plus a strong one-cell halo, ocean boundary extended from the coast exit
/// points, peel rejection, integer max flow) and retires folded plates as
/// `Folded`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReorgLaw {
    #[default]
    ScriptedV1,
    ScriptedV2,
    /// SEED-WORLD-COMPOSE (lane CP, 2026-10-13) `plate-reorg-scripted-v3`:
    /// v2 plus age-driven subduction initiation on the v2 cut machinery
    /// (run through `reorg_pass_v3`, which carries the initiation age).
    ScriptedV3,
    /// SEED-TRANSFORMS (lane TF, 2026-10-13) `plate-reorg-scripted-v4`: v3
    /// with `InitiationLimits` on the age-triggered initiation (run through
    /// `reorg_pass_v4`).
    ScriptedV4,
    /// SEED-OLD-OCEAN (lane OR, 2026-10-13) `plate-reorg-scripted-v5`: v4
    /// plus `MarginInitiation` (run through `reorg_pass_v5`).
    ScriptedV5,
}

/// SEED-WORLD-COMPOSE (lane CP) v3 initiation, fixed by the law version:
/// ocean cells weaken linearly with crust age from `INIT_AGE_WEAK_START_MY`
/// to `INIT_AGE_WEAK_FULL_MY` (to the weak-zone ratio); continental cells
/// are `INIT_CONT_STRENGTH` (cuts avoid them); a cut leaving more than
/// `INIT_CONT_SPLIT_MAX` and less than 1 - that of the group's continent on
/// side B splits a continent and is rejected.
const INIT_AGE_WEAK_START_MY: f64 = 100.0;
const INIT_AGE_WEAK_FULL_MY: f64 = 180.0;
const INIT_CONT_STRENGTH: f64 = 10.0;
const INIT_CONT_SPLIT_MAX: f64 = 0.05;

/// Fixed by the law version (PREREG.txt "fixed by version").
const BAND_LO: f64 = 0.25;
const BAND_HI: f64 = 0.75;
const THICKNESS_WEAK_SCALE_M: f64 = 5_000.0;
const SHORTENING_WEAK_SCALE_M: f64 = 5_000.0;
const NOISE_SIGMA: f64 = 0.35;
const NOISE_WAVELENGTHS_M: [f64; 3] = [2_400_000.0, 1_200_000.0, 600_000.0];
const NOISE_WAVES_PER_OCTAVE: usize = 8;
const MARGIN_SINK_DISTANCE_M: f64 = 600_000.0;
/// Half-width (m) of the free band either side of the plane through the
/// continent that splits "front" ocean (sink) from "back" ocean (source)
/// in a passive-margin initiation (PREREG revision 1).
const MARGIN_SIDE_BAND_M: f64 = 1_000_000.0;
const MAX_AXIS_ATTEMPTS: usize = 3;
const MAX_INITIATIONS_PER_EVENT: usize = 3;
const SIDE_MIN_CONT_SHARE: f64 = 0.15;
/// Weak length of a cut: edges whose mean strength is below this.
const WEAK_LOG_THRESHOLD: f64 = 0.6;

/// Constants of `plate-reorg-scripted-v1` (plus the run's rate scale and
/// plume head radius it reads).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReorgParams {
    pub rift_opening_rate_deg_per_my: f64,
    pub interval_min_my: f64,
    pub interval_max_my: f64,
    pub slab_pull_bias: f64,
    pub slab_circumference_reference: f64,
    pub rate_floor_fraction: f64,
    pub supercontinent_share: f64,
    pub supercontinent_tenure_my: f64,
    pub breakup_piece_max_share: f64,
    pub breakup_cascade_my: f64,
    pub rift_persistence_my: f64,
    pub weak_zone_strength_ratio: f64,
    pub plate_count_min: u32,
    pub plate_count_max: u32,
    pub plate_count_min_area_fraction: f64,
    pub subduction_initiation_rate_deg_per_my: f64,
    pub plume_head_interval_my: f64,
    pub plume_tail_lifetime_my: f64,
    pub plume_split_min_share: f64,
    /// `tectonics.rate_scale_deg_per_my` of the run.
    pub rate_scale_deg_per_my: f64,
    /// `fixed-plume-chains-v1` head radius (m).
    pub plume_head_radius_m: f64,
}

impl ReorgParams {
    pub fn validate(&self) -> Result<(), String> {
        let pos = [
            ("interval_min_my", self.interval_min_my),
            ("interval_max_my", self.interval_max_my),
            ("slab_circumference_reference", self.slab_circumference_reference),
            ("supercontinent_share", self.supercontinent_share),
            ("breakup_piece_max_share", self.breakup_piece_max_share),
            ("weak_zone_strength_ratio", self.weak_zone_strength_ratio),
            ("plume_head_interval_my", self.plume_head_interval_my),
            ("plume_tail_lifetime_my", self.plume_tail_lifetime_my),
            ("plume_head_radius_m", self.plume_head_radius_m),
        ];
        for (w, v) in pos {
            if !(v.is_finite() && v > 0.0) {
                return Err(format!("{LAW}: {w} must be finite and positive, got {v}"));
            }
        }
        if self.interval_min_my > self.interval_max_my || self.plate_count_min > self.plate_count_max {
            return Err(format!("{LAW}: interval or plate-count range inverted"));
        }
        Ok(())
    }
}

/// A breakup / initiation / plume split whose relative rotation survives
/// re-solves while younger than the persistence time.
#[derive(Debug, Clone, PartialEq)]
pub struct Relation {
    pub parent: u32,
    pub child: u32,
    /// Child minus parent angular velocity (rad/My, unit-sphere frame).
    pub omega_rel: [f64; 3],
    pub t_my: f64,
    /// 0 breakup, 1 cascade breakup, 2 subduction initiation, 3 plume split.
    pub kind: u8,
}

/// One logged reorganization event.
#[derive(Debug, Clone, PartialEq)]
pub struct ReorgEvent {
    pub t_my: f64,
    /// "resolve", "breakup", "cascade", "initiation", "plume-split",
    /// "merge", "fragment-fold", "fragment-plate".
    pub kind: String,
    pub root_a: u32,
    pub root_b: u32,
    pub cells_moved: u64,
    /// Share of the continental crust in the group before the event and in
    /// side B after it (breakups and splits).
    pub group_share: f64,
    pub side_b_share: f64,
    /// Cut length (m) and the part through weak cells (mean strength < 0.6).
    pub cut_length_m: f64,
    pub weak_length_m: f64,
    /// Groups above the counting area after the event.
    pub groups_counted: u32,
}

/// Run state of the law (checkpointed by the caller).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReorgState {
    /// (group root, time it first held a supercontinent share).
    pub tenure: Vec<(u32, f64)>,
    /// (group root, time a cascade breakup is due).
    pub cascade: Vec<(u32, f64)>,
    pub relations: Vec<Relation>,
    /// Index of the next scheduled re-solve.
    pub next_event: u32,
    /// Counters naming the RNG streams.
    pub breakups: u32,
    pub initiations: u32,
    pub events: Vec<ReorgEvent>,
}

/// What one step's pass did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReorgStepReport {
    pub resolves: u32,
    pub breakups: u32,
    pub initiations: u32,
    pub merges: u32,
    pub fragments: u32,
    pub groups_counted: u32,
    pub largest_group_share: f64,
    /// Plates retired this step (no cells, no parcels).
    pub retired: u32,
    /// SEED-OLD-OCEAN (lane OR) logs (v4 and v5 only; feed no output):
    /// at re-solves, the length (m) of 4-neighbour edges between suture
    /// groups where either cell is ridge class, and the part that
    /// converges (normal rate > 0, moving faster than 0.5 cm/yr, more than
    /// 30 degrees from the trace) under the new poles.
    pub resolve_ridge_m: f64,
    pub resolve_ridge_to_convergent_m: f64,
    /// v5 margin initiations this step: (seed age My, seed cell, side-B
    /// area m^2, trench length along the coast m).
    pub margin_seeds: Vec<(f64, u32, f64, f64)>,
    pub margin_rejected: u32,
    /// PLATE-POLES `rigid-group-contacts-v1`: detached member sets released
    /// from their suture group this step (0 without the law).
    pub released: u32,
}

/// The scheduled re-solve times in (0, horizon_my]: gaps uniform in
/// [min, max] from the tagged stream (prefix-stable in the horizon).
pub fn reorg_schedule(seed_path: &SeedPath, min_my: f64, max_my: f64, horizon_my: f64) -> Vec<f64> {
    let mut rng = SplitMix64::new(seed_from_path_tagged(seed_path, "plate-reorg-scripted-v1/schedule"));
    let mut out = Vec::new();
    let mut t = 0.0;
    loop {
        t += min_my + (max_my - min_my) * rng.next_f64();
        if t > horizon_my || out.len() > 100_000 {
            break;
        }
        out.push(t);
    }
    out
}

/// Plume heads as a renewal process: onset_0 = -lead + Exp(interval),
/// onset_k+1 = onset_k + Exp(interval); position uniform on the sphere and
/// cut azimuth uniform, four draws per plume in id order (prefix-stable in
/// the horizon).
pub fn seed_plumes_renewal(interval_my: f64, lead_my: f64, horizon_my: f64, seed_path: &SeedPath) -> Vec<Plume> {
    let mut rng = SplitMix64::new(seed_from_path_tagged(seed_path, "plate-reorg-scripted-v1/plumes"));
    let mut out = Vec::new();
    let mut t = -lead_my;
    loop {
        let u = rng.next_f64();
        let z = 2.0 * rng.next_f64() - 1.0;
        let phi = 2.0 * PI * rng.next_f64();
        let az = 2.0 * PI * rng.next_f64();
        t += -interval_my * (1.0 - u).ln();
        if t > horizon_my || out.len() > 100_000 {
            break;
        }
        let r = (1.0 - z * z).max(0.0).sqrt();
        out.push(Plume {
            id: out.len(),
            position: [r * phi.cos(), r * phi.sin(), z],
            onset_my: t,
            split_azimuth_rad: az,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// Vector helpers.

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = norm(a);
    if n > 1e-12 && n.is_finite() {
        Some(scale(a, 1.0 / n))
    } else {
        None
    }
}
fn tangent(v: [f64; 3], at: [f64; 3]) -> [f64; 3] {
    let d = dot(v, at);
    [v[0] - d * at[0], v[1] - d * at[1], v[2] - d * at[2]]
}
fn unit_dir(rng: &mut SplitMix64) -> [f64; 3] {
    let z = 2.0 * rng.next_f64() - 1.0;
    let phi = 2.0 * PI * rng.next_f64();
    let r = (1.0 - z * z).max(0.0).sqrt();
    [r * phi.cos(), r * phi.sin(), z]
}
fn tangent_basis(c: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let a = if c[2].abs() < 0.9 { [0.0, 0.0, 1.0] } else { [1.0, 0.0, 0.0] };
    let e1 = unit(cross(c, a)).unwrap_or([1.0, 0.0, 0.0]);
    let e2 = cross(c, e1);
    (e1, e2)
}

/// Euler pole of an angular velocity vector (rad/My).
pub fn pole_of(w: [f64; 3]) -> Result<EulerPoleV0, String> {
    let mag = norm(w);
    if mag > 1.0e-15 {
        let (lat, lon) = unit_to_lat_lon(scale(w, 1.0 / mag)).map_err(|e| e.to_string())?;
        Ok(EulerPoleV0 { lat_deg: lat, lon_deg: lon, rate_deg_per_my: mag.to_degrees() })
    } else {
        Ok(EulerPoleV0 { lat_deg: 90.0, lon_deg: 0.0, rate_deg_per_my: 0.0 })
    }
}

fn omega_of(plates: &[PlateRecordV0], id: PlateId) -> Result<[f64; 3], String> {
    let rec = plates
        .iter()
        .find(|r| r.plate_id == id)
        .ok_or_else(|| format!("{LAW}: unknown plate {}", id.0))?;
    let pole = rec.euler_pole.ok_or_else(|| format!("{LAW}: plate {} has no pole", id.0))?;
    angular_velocity(&pole).map_err(|e| e.to_string())
}

// ---------------------------------------------------------------------------
// Grid geometry and neighbour graph (8-neighbour, longitude wraps, no
// links across the poles), with the anisotropic cut-metric weights.

struct Geom {
    nlat: usize,
    nlon: usize,
    radius: f64,
    cents: Vec<[f64; 3]>,
    areas: Vec<f64>,
    /// Per row: cut-metric weights (m) for E-W, N-S and diagonal links and
    /// the shared edge length (m) of the E-W and the N (upper) links.
    w_ew: Vec<f64>,
    w_ns: Vec<f64>,
    w_diag: Vec<f64>,
    edge_ew: Vec<f64>,
    edge_n: Vec<f64>,
}

impl Geom {
    fn new(grid: &GridSpec) -> Result<Self, String> {
        let (nlat, nlon) = (grid.nlat(), grid.nlon());
        let radius = grid.radius_m();
        let cents = centroids(grid).map_err(|e| e.to_string())?;
        let areas = cell_areas(grid).map_err(|e| e.to_string())?;
        let dlat = PI / nlat as f64;
        let dlon = 2.0 * PI / nlon as f64;
        let mut w_ew = Vec::with_capacity(nlat);
        let mut w_ns = Vec::with_capacity(nlat);
        let mut w_diag = Vec::with_capacity(nlat);
        let mut edge_ew = Vec::with_capacity(nlat);
        let mut edge_n = Vec::with_capacity(nlat);
        for i in 0..nlat {
            let lat_c = -PI / 2.0 + (i as f64 + 0.5) * dlat;
            let lat_top = -PI / 2.0 + (i as f64 + 1.0) * dlat;
            let dx = (radius * lat_c.cos() * dlon).max(1.0);
            let dy = radius * dlat;
            let th = dy.atan2(dx);
            // Boykov-Kolmogorov 2003: w_k = A dphi_k / (2 |e_k|).
            w_ew.push(dy * th / 2.0);
            w_ns.push(dx * (PI / 2.0 - th) / 2.0);
            w_diag.push(dx * dy * (PI / 4.0) / (2.0 * (dx * dx + dy * dy).sqrt()));
            edge_ew.push(dy);
            edge_n.push(radius * lat_top.cos().max(0.0) * dlon);
        }
        Ok(Self { nlat, nlon, radius, cents, areas, w_ew, w_ns, w_diag, edge_ew, edge_n })
    }

    /// 8 neighbours of c with their cut weights (m).
    fn nbrs8(&self, c: usize, out: &mut Vec<(usize, f64)>) {
        out.clear();
        let (i, j) = ((c / self.nlon) as i64, (c % self.nlon) as i64);
        for di in -1..=1i64 {
            let ti = i + di;
            if ti < 0 || ti >= self.nlat as i64 {
                continue;
            }
            for dj in -1..=1i64 {
                if di == 0 && dj == 0 {
                    continue;
                }
                let k = (ti * self.nlon as i64 + (j + dj).rem_euclid(self.nlon as i64)) as usize;
                let w = if di == 0 {
                    self.w_ew[i as usize]
                } else if dj == 0 {
                    0.5 * (self.w_ns[i as usize] + self.w_ns[ti as usize])
                } else {
                    0.5 * (self.w_diag[i as usize] + self.w_diag[ti as usize])
                };
                out.push((k, w));
            }
        }
    }

    /// 4 neighbours of c with the shared edge length (m).
    fn nbrs4(&self, c: usize, out: &mut Vec<(usize, f64)>) {
        out.clear();
        let (i, j) = (c / self.nlon, c % self.nlon);
        let nl = self.nlon;
        out.push((i * nl + (j + 1) % nl, self.edge_ew[i]));
        out.push((i * nl + (j + nl - 1) % nl, self.edge_ew[i]));
        if i + 1 < self.nlat {
            out.push(((i + 1) * nl + j, self.edge_n[i]));
        }
        if i > 0 {
            out.push(((i - 1) * nl + j, self.edge_n[i - 1]));
        }
    }
}

// ---------------------------------------------------------------------------
// Groups.

struct Group {
    root: PlateId,
    members: Vec<PlateId>,
    cells: Vec<usize>,
    area: f64,
    cont_area: f64,
    ocean_area: f64,
}

struct Groups {
    list: Vec<Group>,
    /// Group index per cell.
    of_cell: Vec<usize>,
    total_cont: f64,
    total_area: f64,
}

fn groups(geom: &Geom, state: &CrustState, plates: &[PlateRecordV0], sut: &SutureState) -> Groups {
    let mut idx: BTreeMap<PlateId, usize> = BTreeMap::new();
    let mut list: Vec<Group> = Vec::new();
    let mut roots: BTreeMap<PlateId, Vec<PlateId>> = BTreeMap::new();
    for p in plates {
        roots.entry(sut.find(p.plate_id)).or_default().push(p.plate_id);
    }
    // Plates present on the grid but missing from the records (should not
    // happen) still form their own group.
    for c in 0..state.cells() {
        let r = sut.find(state.plate[c]);
        roots.entry(r).or_default();
    }
    for (r, mut m) in roots {
        m.sort();
        m.dedup();
        idx.insert(r, list.len());
        list.push(Group { root: r, members: m, cells: Vec::new(), area: 0.0, cont_area: 0.0, ocean_area: 0.0 });
    }
    let mut of_cell = vec![0usize; state.cells()];
    let (mut tc, mut ta) = (0.0, 0.0);
    for c in 0..state.cells() {
        let g = idx[&sut.find(state.plate[c])];
        of_cell[c] = g;
        let a = geom.areas[c];
        list[g].cells.push(c);
        list[g].area += a;
        ta += a;
        if state.kind[c] == CrustKindCodeV0::Continental {
            list[g].cont_area += a;
            tc += a;
        } else {
            list[g].ocean_area += a;
        }
    }
    Groups { list, of_cell, total_cont: tc, total_area: ta }
}

fn counted(gs: &Groups, frac: f64) -> u32 {
    gs.list.iter().filter(|g| g.area >= frac * gs.total_area && !g.cells.is_empty()).count() as u32
}

// ---------------------------------------------------------------------------
// Spectral noise (unit standard deviation).

struct Noise {
    waves: Vec<([f64; 3], f64, f64)>,
}

impl Noise {
    fn new(rng: &mut SplitMix64, radius_m: f64) -> Self {
        let per = (2.0 / (NOISE_WAVES_PER_OCTAVE * NOISE_WAVELENGTHS_M.len()) as f64).sqrt();
        let mut waves = Vec::new();
        for &lam in &NOISE_WAVELENGTHS_M {
            let k = 2.0 * PI * radius_m / lam;
            for _ in 0..NOISE_WAVES_PER_OCTAVE {
                let d = unit_dir(rng);
                let ph = 2.0 * PI * rng.next_f64();
                waves.push((scale(d, k), ph, per));
            }
        }
        Self { waves }
    }
    fn eval(&self, x: [f64; 3]) -> f64 {
        self.waves.iter().map(|(k, ph, a)| a * (dot(*k, x) + ph).cos()).sum()
    }
}

// ---------------------------------------------------------------------------
// Strength and the min cut.

fn weighted_quantile(vals: &mut [(f64, f64)], q: f64) -> f64 {
    vals.sort_by(|a, b| a.0.total_cmp(&b.0));
    let tot: f64 = vals.iter().map(|v| v.1).sum();
    let mut acc = 0.0;
    for v in vals.iter() {
        acc += v.1;
        if acc >= q * tot {
            return v.0;
        }
    }
    vals.last().map(|v| v.0).unwrap_or(0.0)
}

/// Cell strength over the group's cells (local order = `cells`).
fn strengths(
    geom: &Geom,
    state: &CrustState,
    cells: &[usize],
    ratio: f64,
    noise: &Noise,
    weak_disc: Option<([f64; 3], f64)>,
) -> Vec<f64> {
    let mut th: Vec<(f64, f64)> = cells
        .iter()
        .filter(|&&c| state.kind[c] == CrustKindCodeV0::Continental)
        .map(|&c| (state.thickness_m[c], geom.areas[c]))
        .collect();
    let h_med = if th.is_empty() { 0.0 } else { weighted_quantile(&mut th, 0.5) };
    cells
        .iter()
        .map(|&c| {
            let mut s = 1.0;
            if state.kind[c] == CrustKindCodeV0::Continental {
                let wt = ((state.thickness_m[c] - h_med).abs() / THICKNESS_WEAK_SCALE_M).clamp(0.0, 1.0);
                let ws = (state.shortening_added_m[c] / SHORTENING_WEAK_SCALE_M).clamp(0.0, 1.0);
                s = 1.0 - (1.0 - ratio) * wt.max(ws);
            }
            s *= (NOISE_SIGMA * noise.eval(geom.cents[c])).exp();
            if let Some((p, ang)) = weak_disc {
                if dot(p, geom.cents[c]).clamp(-1.0, 1.0).acos() <= ang {
                    s *= ratio;
                }
            }
            s
        })
        .collect()
}

/// Dinic max flow on f64 capacities.
struct Flow {
    head: Vec<usize>,
    to: Vec<usize>,
    cap: Vec<f64>,
    next: Vec<usize>,
}

const NIL: usize = usize::MAX;
const EPS: f64 = 1e-9;

impl Flow {
    fn new(n: usize) -> Self {
        Self { head: vec![NIL; n], to: Vec::new(), cap: Vec::new(), next: Vec::new() }
    }
    fn arc(&mut self, u: usize, v: usize, c_uv: f64, c_vu: f64) {
        self.to.push(v);
        self.cap.push(c_uv);
        self.next.push(self.head[u]);
        self.head[u] = self.to.len() - 1;
        self.to.push(u);
        self.cap.push(c_vu);
        self.next.push(self.head[v]);
        self.head[v] = self.to.len() - 1;
    }
    fn bfs(&self, s: usize, level: &mut [i32]) {
        level.iter_mut().for_each(|l| *l = -1);
        let mut q = std::collections::VecDeque::new();
        level[s] = 0;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            let mut e = self.head[u];
            while e != NIL {
                let v = self.to[e];
                if level[v] < 0 && self.cap[e] > EPS {
                    level[v] = level[u] + 1;
                    q.push_back(v);
                }
                e = self.next[e];
            }
        }
    }
    fn dfs(&mut self, u: usize, t: usize, f: f64, level: &[i32], it: &mut [usize]) -> f64 {
        if u == t {
            return f;
        }
        while it[u] != NIL {
            let e = it[u];
            let v = self.to[e];
            if self.cap[e] > EPS && level[v] == level[u] + 1 {
                let d = self.dfs(v, t, f.min(self.cap[e]), level, it);
                if d > EPS {
                    self.cap[e] -= d;
                    self.cap[e ^ 1] += d;
                    return d;
                }
            }
            it[u] = self.next[e];
        }
        0.0
    }
    fn max_flow(&mut self, s: usize, t: usize) -> f64 {
        let n = self.head.len();
        let mut level = vec![-1i32; n];
        let mut total = 0.0;
        loop {
            self.bfs(s, &mut level);
            if level[t] < 0 {
                break;
            }
            let mut it = self.head.clone();
            loop {
                let f = self.dfs(s, t, f64::INFINITY, &level, &mut it);
                if f <= EPS {
                    break;
                }
                total += f;
            }
        }
        total
    }
    fn reachable(&self, s: usize) -> Vec<bool> {
        let mut seen = vec![false; self.head.len()];
        let mut stack = vec![s];
        seen[s] = true;
        while let Some(u) = stack.pop() {
            let mut e = self.head[u];
            while e != NIL {
                let v = self.to[e];
                if !seen[v] && self.cap[e] > EPS {
                    seen[v] = true;
                    stack.push(v);
                }
                e = self.next[e];
            }
        }
        seen
    }
}

/// Result of a cut over a group's cells.
struct Cut {
    /// Per local cell: on side B (the sink side).
    side_b: Vec<bool>,
    cut_length_m: f64,
    weak_length_m: f64,
}

/// Min cut over `cells` (ascending) between local `src` and `snk` sets.
/// Capacity of a link = cut weight x mean strength x (weak ratio if the two
/// cells belong to different member plates: an old suture). Side B = cells
/// not reachable from the source in the residual graph and connected to the
/// sink through such cells.
#[allow(clippy::too_many_arguments)]
fn min_cut(
    geom: &Geom,
    state: &CrustState,
    cells: &[usize],
    in_group: &[usize],
    strength: &[f64],
    src: &[bool],
    snk: &[bool],
    ratio: f64,
) -> Cut {
    let m = cells.len();
    let (s, t) = (m, m + 1);
    let mut fl = Flow::new(m + 2);
    let mut nb = Vec::with_capacity(8);
    // Capacities in km so values stay O(1-100).
    for (a, &c) in cells.iter().enumerate() {
        geom.nbrs8(c, &mut nb);
        for &(k, w) in &nb {
            let b = in_group[k];
            if b == NIL || b <= a {
                continue;
            }
            let mut cap = w * 1e-3 * 0.5 * (strength[a] + strength[b]);
            if state.plate[c] != state.plate[k] {
                cap *= ratio;
            }
            fl.arc(a, b, cap, cap);
        }
        if src[a] {
            fl.arc(s, a, 1e30, 0.0);
        }
        if snk[a] {
            fl.arc(a, t, 1e30, 0.0);
        }
    }
    fl.max_flow(s, t);
    let reach = fl.reachable(s);
    // B: not reachable, connected to the sink through non-reachable cells.
    let mut side_b = vec![false; m];
    let mut stack: Vec<usize> = (0..m).filter(|&a| snk[a] && !reach[a]).collect();
    for &a in &stack {
        side_b[a] = true;
    }
    let mut nb4 = Vec::with_capacity(4);
    while let Some(a) = stack.pop() {
        geom.nbrs8(cells[a], &mut nb);
        for &(k, _) in &nb {
            let b = in_group[k];
            if b != NIL && !reach[b] && !side_b[b] {
                side_b[b] = true;
                stack.push(b);
            }
        }
    }
    // Cut length along 4-neighbour edges between A and B.
    let (mut len, mut weak) = (0.0, 0.0);
    for (a, &c) in cells.iter().enumerate() {
        if side_b[a] {
            continue;
        }
        geom.nbrs4(c, &mut nb4);
        for &(k, e) in &nb4 {
            let b = in_group[k];
            if b != NIL && side_b[b] {
                len += e;
                let mut sm = 0.5 * (strength[a] + strength[b]);
                if state.plate[c] != state.plate[k] {
                    sm *= ratio;
                }
                if sm < WEAK_LOG_THRESHOLD {
                    weak += e;
                }
            }
        }
    }
    Cut { side_b, cut_length_m: len, weak_length_m: weak }
}

fn local_index(n: usize, cells: &[usize]) -> Vec<usize> {
    let mut v = vec![NIL; n];
    for (a, &c) in cells.iter().enumerate() {
        v[c] = a;
    }
    v
}

/// Area-weighted centroid (unit) of the cells selected by `f`.
fn centroid(geom: &Geom, cells: impl Iterator<Item = usize>) -> Option<[f64; 3]> {
    let mut s = [0.0; 3];
    for c in cells {
        s = add(s, scale(geom.cents[c], geom.areas[c]));
    }
    unit(s)
}

/// Mean midpoint (unit) of the A-B 4-neighbour edges.
fn cut_centre(geom: &Geom, cells: &[usize], in_group: &[usize], side_b: &[bool]) -> Option<[f64; 3]> {
    let mut s = [0.0; 3];
    let mut nb4 = Vec::with_capacity(4);
    for (a, &c) in cells.iter().enumerate() {
        if side_b[a] {
            continue;
        }
        geom.nbrs4(c, &mut nb4);
        for &(k, e) in &nb4 {
            let b = in_group[k];
            if b != NIL && side_b[b] {
                if let Some(mid) = unit(add(geom.cents[c], geom.cents[k])) {
                    s = add(s, scale(mid, e));
                }
            }
        }
    }
    unit(s)
}

// ---------------------------------------------------------------------------
// Applying a split / relabel.

struct SplitOutcome {
    root_a: PlateId,
    root_b: PlateId,
    cells_moved: u64,
}

/// Move the side-B cells of group `g` into a new suture group: members
/// wholly on B keep their ids, split members give their B cells to a new
/// child id. Union-find is rewritten into the two groups; A keeps the old
/// pole, B gets old pole + omega_rel. Parcels follow (nearest cell).
#[allow(clippy::too_many_arguments)]
fn apply_split(
    geom: &Geom,
    grid: &GridSpec,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    g: &Group,
    side_b: &[bool],
    omega_rel: [f64; 3],
    t_my: f64,
) -> Result<SplitOutcome, String> {
    let w_old = omega_of(plates, g.root)?;
    let mut total: BTreeMap<PlateId, (u64, u64, f64, f64)> = BTreeMap::new(); // (cells, on B, area A, area B)
    for m in &g.members {
        total.insert(*m, (0, 0, 0.0, 0.0));
    }
    for (a, &c) in g.cells.iter().enumerate() {
        let e = total.entry(state.plate[c]).or_insert((0, 0, 0.0, 0.0));
        e.0 += 1;
        if side_b[a] {
            e.1 += 1;
            e.3 += geom.areas[c];
        } else {
            e.2 += geom.areas[c];
        }
    }
    let mut next_id = plates.iter().map(|r| r.plate_id.0).max().unwrap_or(0) + 1;
    let mut a_set: Vec<(PlateId, f64)> = Vec::new();
    let mut b_set: Vec<(PlateId, f64)> = Vec::new();
    let mut child_of: BTreeMap<PlateId, PlateId> = BTreeMap::new();
    for (&m, &(n, nb, aa, ab)) in &total {
        if n > 0 && nb == n {
            b_set.push((m, ab));
        } else if nb == 0 {
            a_set.push((m, aa));
        } else {
            a_set.push((m, aa));
            let child = PlateId(next_id);
            next_id += 1;
            child_of.insert(m, child);
            b_set.push((child, ab));
            plates.push(PlateRecordV0 {
                plate_id: child,
                parent_plate: Some(m),
                birth_t_my: t_my,
                euler_pole: Some(pole_of(w_old)?),
                retirement: None,
            });
        }
    }
    if b_set.is_empty() || a_set.iter().all(|x| x.1 <= 0.0) {
        return Err(format!("{LAW}: degenerate split of group {}", g.root.0));
    }
    let mut moved: Vec<(usize, PlateId, PlateId)> = Vec::new();
    for (a, &c) in g.cells.iter().enumerate() {
        if side_b[a] {
            if let Some(&ch) = child_of.get(&state.plate[c]) {
                moved.push((c, state.plate[c], ch));
                state.plate[c] = ch;
            }
        }
    }
    relabel_parcels(grid, layer, &moved).map_err(|e| e.to_string())?;
    let best = |v: &[(PlateId, f64)]| -> PlateId {
        let mut b = v[0];
        for x in v {
            if x.1 > b.1 || (x.1 == b.1 && x.0 < b.0) {
                b = *x;
            }
        }
        b.0
    };
    let ra = if a_set.iter().any(|x| x.0 == g.root && x.1 > 0.0) { g.root } else { best(&a_set) };
    let rb = best(&b_set);
    for (x, _) in a_set.iter().chain(b_set.iter()) {
        sut.parent.remove(x);
    }
    for (x, _) in &a_set {
        if *x != ra {
            sut.parent.insert(*x, ra);
        }
    }
    for (x, _) in &b_set {
        if *x != rb {
            sut.parent.insert(*x, rb);
        }
    }
    let pa = pole_of(w_old)?;
    let pb = pole_of(add(w_old, omega_rel))?;
    for r in plates.iter_mut() {
        if a_set.iter().any(|x| x.0 == r.plate_id) {
            r.euler_pole = Some(pa);
        } else if b_set.iter().any(|x| x.0 == r.plate_id) {
            r.euler_pole = Some(pb);
        }
    }
    Ok(SplitOutcome { root_a: ra, root_b: rb, cells_moved: g.cells.iter().enumerate().filter(|(a, _)| side_b[*a]).count() as u64 })
}

/// Relabel cells (and the parcels whose nearest cell they are) to `to`.
fn relabel(
    grid: &GridSpec,
    state: &mut CrustState,
    layer: &mut ParcelLayer,
    cells: &[usize],
    to: PlateId,
) -> Result<(), String> {
    let mut moved = Vec::with_capacity(cells.len());
    for &c in cells {
        if state.plate[c] != to {
            moved.push((c, state.plate[c], to));
            state.plate[c] = to;
        }
    }
    relabel_parcels(grid, layer, &moved).map_err(|e| e.to_string())?;
    Ok(())
}

/// The neighbouring plate (not in `exclude`) sharing the most links with
/// `cells`; ties to the lowest id.
fn most_touching(geom: &Geom, state: &CrustState, cells: &[usize], exclude: impl Fn(PlateId) -> bool) -> Option<PlateId> {
    let mut links: BTreeMap<PlateId, u64> = BTreeMap::new();
    let mut nb = Vec::with_capacity(8);
    for &c in cells {
        geom.nbrs8(c, &mut nb);
        for &(k, _) in &nb {
            let q = state.plate[k];
            if !exclude(q) {
                *links.entry(q).or_insert(0) += 1;
            }
        }
    }
    let mut best: Option<(u64, PlateId)> = None;
    for (q, n) in links {
        if best.map_or(true, |(bn, _)| n > bn) {
            best = Some((n, q));
        }
    }
    best.map(|b| b.1)
}

// ---------------------------------------------------------------------------
// The passes.

/// Context the passes share.
pub struct ReorgCtx<'a> {
    pub grid: &'a GridSpec,
    pub seed_path: &'a SeedPath,
    pub params: &'a ReorgParams,
    pub schedule: &'a [f64],
}

/// Split group `gi` with a cut between `src` and `snk`; returns the event.
#[allow(clippy::too_many_arguments)]
fn cut_and_split(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    g: &Group,
    cut: &Cut,
    omega_dir_sign: f64,
    rate_deg_per_my: f64,
    p_override: Option<[f64; 3]>,
    t_my: f64,
) -> Result<Option<(SplitOutcome, [f64; 3])>, String> {
    let in_group = local_index(state.cells(), &g.cells);
    let p = match p_override.or_else(|| cut_centre(geom, &g.cells, &in_group, &cut.side_b)) {
        Some(p) => p,
        None => return Ok(None),
    };
    let is_cont = |c: usize| state.kind[c] == CrustKindCodeV0::Continental;
    let ca = centroid(geom, g.cells.iter().enumerate().filter(|(a, c)| !cut.side_b[*a] && is_cont(**c)).map(|(_, c)| *c))
        .or_else(|| centroid(geom, g.cells.iter().enumerate().filter(|(a, _)| !cut.side_b[*a]).map(|(_, c)| *c)));
    let cb = centroid(geom, g.cells.iter().enumerate().filter(|(a, c)| cut.side_b[*a] && is_cont(**c)).map(|(_, c)| *c))
        .or_else(|| centroid(geom, g.cells.iter().enumerate().filter(|(a, _)| cut.side_b[*a]).map(|(_, c)| *c)));
    let (Some(ca), Some(cb)) = (ca, cb) else { return Ok(None) };
    let Some(n) = unit(tangent([cb[0] - ca[0], cb[1] - ca[1], cb[2] - ca[2]], p)) else {
        return Ok(None);
    };
    let Some(ax) = unit(cross(p, n)) else { return Ok(None) };
    let omega_rel = scale(ax, omega_dir_sign * rate_deg_per_my.to_radians());
    let out = apply_split(geom, ctx.grid, state, plates, layer, sut, g, &cut.side_b, omega_rel, t_my)?;
    Ok(Some((out, omega_rel)))
}

fn cont_share_b(geom: &Geom, state: &CrustState, g: &Group, side_b: &[bool]) -> f64 {
    let mut b = 0.0;
    for (a, &c) in g.cells.iter().enumerate() {
        if side_b[a] && state.kind[c] == CrustKindCodeV0::Continental {
            b += geom.areas[c];
        }
    }
    if g.cont_area > 0.0 {
        b / g.cont_area
    } else {
        0.0
    }
}

/// Breakup of group index `gi` along a weak-zone min cut (seeded axis).
#[allow(clippy::too_many_arguments)]
fn breakup(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    g: &Group,
    total_cont: f64,
    kind: &str,
    t_my: f64,
) -> Result<Option<ReorgEvent>, String> {
    let p = ctx.params;
    let k = rs.breakups;
    rs.breakups += 1;
    let mut rng = SplitMix64::new(seed_from_path_tagged(ctx.seed_path, &format!("plate-reorg-scripted-v1/breakup/{k}")));
    let noise = Noise::new(&mut rng, geom.radius);
    let cont: Vec<usize> = g.cells.iter().cloned().filter(|&c| state.kind[c] == CrustKindCodeV0::Continental).collect();
    let Some(gc) = centroid(geom, cont.iter().cloned()) else { return Ok(None) };
    let strength = strengths(geom, state, &g.cells, p.weak_zone_strength_ratio, &noise, None);
    let in_group = local_index(state.cells(), &g.cells);
    let (e1, e2) = tangent_basis(gc);
    for _attempt in 0..MAX_AXIS_ATTEMPTS {
        let az = 2.0 * PI * rng.next_f64();
        let axis = add(scale(e1, az.cos()), scale(e2, az.sin()));
        let mut sv: Vec<(f64, f64)> = cont.iter().map(|&c| (dot(geom.cents[c], axis), geom.areas[c])).collect();
        let lo = weighted_quantile(&mut sv, BAND_LO);
        let hi = weighted_quantile(&mut sv, BAND_HI);
        let src: Vec<bool> = g.cells.iter().map(|&c| dot(geom.cents[c], axis) <= lo).collect();
        let snk: Vec<bool> = g.cells.iter().map(|&c| dot(geom.cents[c], axis) >= hi).collect();
        let cut = min_cut(geom, state, &g.cells, &in_group, &strength, &src, &snk, p.weak_zone_strength_ratio);
        let sb = cont_share_b(geom, state, g, &cut.side_b);
        if !(SIDE_MIN_CONT_SHARE..=1.0 - SIDE_MIN_CONT_SHARE).contains(&sb) {
            continue;
        }
        let Some((out, w)) = cut_and_split(ctx, geom, state, plates, layer, sut, g, &cut, 1.0, p.rift_opening_rate_deg_per_my, None, t_my)? else {
            continue;
        };
        rs.relations.push(Relation {
            parent: out.root_a.0,
            child: out.root_b.0,
            omega_rel: w,
            t_my,
            kind: if kind == "cascade" { 1 } else { 0 },
        });
        return Ok(Some(ReorgEvent {
            t_my,
            kind: kind.into(),
            root_a: out.root_a.0,
            root_b: out.root_b.0,
            cells_moved: out.cells_moved,
            group_share: g.cont_area / total_cont.max(1.0),
            side_b_share: sb,
            cut_length_m: cut.cut_length_m,
            weak_length_m: cut.weak_length_m,
            groups_counted: 0,
        }));
    }
    Ok(None)
}

/// Subduction initiation in the group with the most ocean.
fn initiation(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    gs: &Groups,
    t_my: f64,
) -> Result<Option<ReorgEvent>, String> {
    let p = ctx.params;
    let mut gi: Option<usize> = None;
    for (i, g) in gs.list.iter().enumerate() {
        if g.cells.is_empty() {
            continue;
        }
        if gi.map_or(true, |b| g.ocean_area > gs.list[b].ocean_area) {
            gi = Some(i);
        }
    }
    let Some(gi) = gi else { return Ok(None) };
    let g = &gs.list[gi];
    let k = rs.initiations;
    rs.initiations += 1;
    let mut rng = SplitMix64::new(seed_from_path_tagged(ctx.seed_path, &format!("plate-reorg-scripted-v1/initiation/{k}")));
    let noise = Noise::new(&mut rng, geom.radius);
    let az = 2.0 * PI * rng.next_f64();
    let in_group = local_index(state.cells(), &g.cells);
    let strength = strengths(geom, state, &g.cells, p.weak_zone_strength_ratio, &noise, None);
    let m = g.cells.len();
    let is_cont: Vec<bool> = g.cells.iter().map(|&c| state.kind[c] == CrustKindCodeV0::Continental).collect();
    let mut src = vec![false; m];
    let mut snk = vec![false; m];
    let mut mode = "initiation-margin";
    let mut ok = false;
    if is_cont.iter().any(|&x| x) {
        // Distance from the group's continent (Dijkstra over the group).
        let mut dist = vec![f64::INFINITY; m];
        let mut heap = std::collections::BinaryHeap::new();
        for a in 0..m {
            if is_cont[a] {
                dist[a] = 0.0;
                heap.push(std::cmp::Reverse((0u64, a)));
            }
        }
        let mut nb = Vec::with_capacity(8);
        while let Some(std::cmp::Reverse((dk, a))) = heap.pop() {
            let d = f64::from_bits(dk);
            if d > dist[a] {
                continue;
            }
            geom.nbrs8(g.cells[a], &mut nb);
            for &(kc, _) in &nb {
                let b = in_group[kc];
                if b == NIL {
                    continue;
                }
                let step = dot(geom.cents[g.cells[a]], geom.cents[kc]).clamp(-1.0, 1.0).acos() * geom.radius;
                let nd = d + step;
                if nd < dist[b] {
                    dist[b] = nd;
                    heap.push(std::cmp::Reverse((nd.to_bits(), b)));
                }
            }
        }
        // One margin only (PREREG revision 1): a seeded axis through the
        // continent's centroid; the far ocean in front of it (s > band) is
        // the sink, the continent and the ocean behind it (s < -band) the
        // source, so the cut is one trench along one margin plus the
        // shortest links to the plate's edges, not a ring round the coast.
        let gcc = centroid(geom, g.cells.iter().enumerate().filter(|(a, _)| is_cont[*a]).map(|(_, c)| *c));
        let mut far = 0.0;
        if let Some(gcc) = gcc {
            let (e1, e2) = tangent_basis(gcc);
            let axis = add(scale(e1, az.cos()), scale(e2, az.sin()));
            let band = (MARGIN_SIDE_BAND_M / geom.radius).sin();
            for a in 0..m {
                let sx = dot(geom.cents[g.cells[a]], axis);
                src[a] = is_cont[a] || sx < -band;
                snk[a] = !is_cont[a] && sx > band && dist[a] > MARGIN_SINK_DISTANCE_M && dist[a].is_finite();
                if snk[a] {
                    far += geom.areas[g.cells[a]];
                }
            }
        }
        ok = far >= 0.1 * g.ocean_area && far > 0.0;
    }
    if !ok {
        mode = "initiation-oceanic";
        let Some(gc) = centroid(geom, g.cells.iter().cloned()) else { return Ok(None) };
        let (e1, e2) = tangent_basis(gc);
        let axis = add(scale(e1, az.cos()), scale(e2, az.sin()));
        let mut sv: Vec<(f64, f64)> = g.cells.iter().map(|&c| (dot(geom.cents[c], axis), geom.areas[c])).collect();
        let lo = weighted_quantile(&mut sv, BAND_LO);
        let hi = weighted_quantile(&mut sv, BAND_HI);
        for (a, &c) in g.cells.iter().enumerate() {
            let s = dot(geom.cents[c], axis);
            src[a] = s <= lo;
            snk[a] = s >= hi;
        }
    }
    let cut = min_cut(geom, state, &g.cells, &in_group, &strength, &src, &snk, p.weak_zone_strength_ratio);
    let area_b: f64 = g.cells.iter().enumerate().filter(|(a, _)| cut.side_b[*a]).map(|(_, &c)| geom.areas[c]).sum();
    let min_a = p.plate_count_min_area_fraction * gs.total_area;
    if area_b < min_a || g.area - area_b < min_a {
        return Ok(None);
    }
    let Some((out, w)) = cut_and_split(ctx, geom, state, plates, layer, sut, g, &cut, -1.0, p.subduction_initiation_rate_deg_per_my, None, t_my)? else {
        return Ok(None);
    };
    rs.relations.push(Relation { parent: out.root_a.0, child: out.root_b.0, omega_rel: w, t_my, kind: 2 });
    Ok(Some(ReorgEvent {
        t_my,
        kind: mode.into(),
        root_a: out.root_a.0,
        root_b: out.root_b.0,
        cells_moved: out.cells_moved,
        group_share: g.cont_area / gs.total_cont.max(1.0),
        side_b_share: cont_share_b(geom, state, g, &cut.side_b),
        cut_length_m: cut.cut_length_m,
        weak_length_m: cut.weak_length_m,
        groups_counted: 0,
    }))
}

/// SEED-WORLD-COMPOSE (lane CP) `plate-reorg-scripted-v3` subduction
/// initiation. `pick` = the age trigger's group (else, the "too few
/// groups" trigger, the group with the most ocean). Strength over the whole
/// group (the v2 integer cut, domain = group): continental cells
/// INIT_CONT_STRENGTH, the one-cell coastal halo at least 1 (craton), other
/// ocean cells weakened linearly by crust age from INIT_AGE_WEAK_START_MY
/// to INIT_AGE_WEAK_FULL_MY down to the weak-zone ratio, times the seeded
/// noise. Axis: from the group's continent centroid toward the centroid of
/// its ocean older than the initiation age (age trigger) or seeded (count
/// trigger; also when the group has no old ocean or no continent). Source:
/// the continent and the ocean behind the band; sink: ocean in front of
/// the band at least MARGIN_SINK_DISTANCE_M from the continent (no
/// continent: quantile bands of the whole group). Cuts splitting a
/// continent or leaving a side below the counted area are rejected; up to
/// MAX_AXIS_ATTEMPTS axes (later ones seeded). Same RNG stream as v1.
#[allow(clippy::too_many_arguments)]
fn initiation_v3(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    gs: &Groups,
    t_my: f64,
    age_my: f64,
    pick: Option<usize>,
    limits: Option<(&InitiationLimits, &mut Vec<[f64; 3]>)>,
) -> Result<Option<ReorgEvent>, String> {
    let mut limits = limits;
    let p = ctx.params;
    let gi = match pick {
        Some(i) => i,
        None => {
            let mut gi: Option<usize> = None;
            for (i, g) in gs.list.iter().enumerate() {
                if g.cells.is_empty() {
                    continue;
                }
                if gi.map_or(true, |b| g.ocean_area > gs.list[b].ocean_area) {
                    gi = Some(i);
                }
            }
            let Some(gi) = gi else { return Ok(None) };
            gi
        }
    };
    let g = &gs.list[gi];
    let k = rs.initiations;
    rs.initiations += 1;
    let mut rng = SplitMix64::new(seed_from_path_tagged(ctx.seed_path, &format!("plate-reorg-scripted-v1/initiation/{k}")));
    let noise = Noise::new(&mut rng, geom.radius);
    let n = state.cells();
    let m = g.cells.len();
    let in_group = local_index(n, &g.cells);
    let is_cont: Vec<bool> = g.cells.iter().map(|&c| state.kind[c] == CrustKindCodeV0::Continental).collect();
    let ratio = p.weak_zone_strength_ratio;
    let mut nb = Vec::with_capacity(8);
    let strength: Vec<f64> = g
        .cells
        .iter()
        .enumerate()
        .map(|(a, &c)| {
            let nf = (NOISE_SIGMA * noise.eval(geom.cents[c])).exp();
            if is_cont[a] {
                return INIT_CONT_STRENGTH * nf;
            }
            geom.nbrs8(c, &mut nb);
            if nb.iter().any(|&(kc, _)| in_group[kc] != NIL && state.kind[kc] == CrustKindCodeV0::Continental) {
                return nf.max(1.0);
            }
            let w = ((state.age_my(c) - INIT_AGE_WEAK_START_MY) / (INIT_AGE_WEAK_FULL_MY - INIT_AGE_WEAK_START_MY)).clamp(0.0, 1.0);
            (1.0 - (1.0 - ratio) * w) * nf
        })
        .collect();
    let dom = DomainV2 { cells: g.cells.clone(), local: in_group.clone(), cont: is_cont.clone(), strength };
    // Distance from the group's continent (Dijkstra over the group), as v1.
    let mut dist = vec![f64::INFINITY; m];
    if is_cont.iter().any(|&x| x) {
        let mut heap = std::collections::BinaryHeap::new();
        for a in 0..m {
            if is_cont[a] {
                dist[a] = 0.0;
                heap.push(std::cmp::Reverse((0u64, a)));
            }
        }
        while let Some(std::cmp::Reverse((dk, a))) = heap.pop() {
            let d = f64::from_bits(dk);
            if d > dist[a] {
                continue;
            }
            geom.nbrs8(g.cells[a], &mut nb);
            for &(kc, _) in &nb {
                let b = in_group[kc];
                if b == NIL {
                    continue;
                }
                let step = dot(geom.cents[g.cells[a]], geom.cents[kc]).clamp(-1.0, 1.0).acos() * geom.radius;
                let nd = d + step;
                if nd < dist[b] {
                    dist[b] = nd;
                    heap.push(std::cmp::Reverse((nd.to_bits(), b)));
                }
            }
        }
    }
    let gcc = centroid(geom, g.cells.iter().enumerate().filter(|(a, _)| is_cont[*a]).map(|(_, c)| *c));
    let old_c = centroid(
        geom,
        g.cells.iter().cloned().filter(|&c| state.kind[c] == CrustKindCodeV0::Oceanic && state.age_my(c) >= age_my),
    );
    let min_a = p.plate_count_min_area_fraction * gs.total_area;
    let band = (MARGIN_SIDE_BAND_M / geom.radius).sin();
    for attempt in 0..MAX_AXIS_ATTEMPTS {
        let az = 2.0 * PI * rng.next_f64();
        let mut src = vec![false; m];
        let mut snk = vec![false; m];
        let mut mode = "initiation-v3-margin";
        let mut ok = false;
        if let Some(gcc) = gcc {
            let (e1, e2) = tangent_basis(gcc);
            let seeded = add(scale(e1, az.cos()), scale(e2, az.sin()));
            // First attempt of the age trigger: toward the old ocean.
            let axis = match (attempt, pick, old_c) {
                (0, Some(_), Some(oc)) => unit(tangent(oc, gcc)).unwrap_or(seeded),
                _ => seeded,
            };
            let mut far = 0.0;
            for a in 0..m {
                let sx = dot(geom.cents[g.cells[a]], axis);
                src[a] = is_cont[a] || sx < -band;
                snk[a] = !is_cont[a] && sx > band && dist[a] > MARGIN_SINK_DISTANCE_M && dist[a].is_finite();
                if snk[a] {
                    far += geom.areas[g.cells[a]];
                }
            }
            ok = far >= 0.1 * g.ocean_area && far > 0.0;
        }
        if !ok {
            mode = "initiation-v3-oceanic";
            let Some(gc) = centroid(geom, g.cells.iter().cloned()) else { return Ok(None) };
            let (e1, e2) = tangent_basis(gc);
            let seeded = add(scale(e1, az.cos()), scale(e2, az.sin()));
            let axis = match (attempt, pick, old_c) {
                (0, Some(_), Some(oc)) => unit(tangent(oc, gc)).unwrap_or(seeded),
                _ => seeded,
            };
            let mut sv: Vec<(f64, f64)> = g.cells.iter().map(|&c| (dot(geom.cents[c], axis), geom.areas[c])).collect();
            let lo = weighted_quantile(&mut sv, BAND_LO);
            let hi = weighted_quantile(&mut sv, BAND_HI);
            for (a, &c) in g.cells.iter().enumerate() {
                let sx = dot(geom.cents[c], axis);
                src[a] = sx <= lo;
                snk[a] = sx >= hi;
            }
        }
        let Some(cv) = domain_cut_v2(geom, state, g, &dom, &src, &snk, ratio) else { continue };
        let area_b: f64 = g.cells.iter().enumerate().filter(|(a, _)| cv.cut.side_b[*a]).map(|(_, &c)| geom.areas[c]).sum();
        if area_b < min_a || g.area - area_b < min_a {
            continue;
        }
        let sb = cont_share_b(geom, state, g, &cv.cut.side_b);
        if g.cont_area > 0.0 && sb > INIT_CONT_SPLIT_MAX && sb < 1.0 - INIT_CONT_SPLIT_MAX {
            continue;
        }
        // SEED-TRANSFORMS (lane TF) v4: slab-age and spacing limits on the
        // age-triggered initiation.
        let mut new_trench: Vec<[f64; 3]> = Vec::new();
        if let (Some(_), Some((lim, started))) = (pick, limits.as_mut()) {
            let (old_frac, near_frac, pts) = initiation_v4_measure(geom, state, g, &in_group, &cv.cut.side_b, lim, started);
            if old_frac < lim.old_slab_fraction_min || near_frac > INIT_V4_NEAR_FRACTION_MAX {
                rs.events.push(ReorgEvent {
                    t_my,
                    kind: if old_frac < lim.old_slab_fraction_min { "initiation-v4-rejected-slab-age".into() } else { "initiation-v4-rejected-spacing".into() },
                    root_a: g.root.0,
                    root_b: 0,
                    cells_moved: 0,
                    group_share: old_frac,
                    side_b_share: near_frac,
                    cut_length_m: cv.cut.cut_length_m,
                    weak_length_m: cv.cut.weak_length_m,
                    groups_counted: 0,
                });
                continue;
            }
            new_trench = pts;
        }
        let Some((out, w)) = cut_and_split(ctx, geom, state, plates, layer, sut, g, &cv.cut, -1.0, p.subduction_initiation_rate_deg_per_my, cv.centre, t_my)? else {
            continue;
        };
        if let Some((_, started)) = limits.as_mut() {
            started.extend(new_trench);
        }
        rs.relations.push(Relation { parent: out.root_a.0, child: out.root_b.0, omega_rel: w, t_my, kind: 2 });
        return Ok(Some(ReorgEvent {
            t_my,
            kind: format!("{mode}{}", if pick.is_some() { "-age" } else { "" }),
            root_a: out.root_a.0,
            root_b: out.root_b.0,
            cells_moved: out.cells_moved,
            group_share: g.cont_area / gs.total_cont.max(1.0),
            side_b_share: sb,
            cut_length_m: cv.cut.cut_length_m,
            weak_length_m: cv.cut.weak_length_m,
            groups_counted: 0,
        }));
    }
    Ok(None)
}

/// SEED-TRANSFORMS (lane TF) v4: over the candidate boundary (4-neighbour
/// edges inside the group between side A and side B), the length fraction
/// whose older flanking cell is ocean at least `slab_age_min_my` old, the
/// length fraction within `trench_spacing_m` of a subduction-class cell or
/// of a boundary started earlier in the pass (`started`), and the side-B
/// edge cells' centres (to extend `started` if accepted).
fn initiation_v4_measure(
    geom: &Geom,
    state: &CrustState,
    g: &Group,
    in_group: &[usize],
    side_b: &[bool],
    lim: &InitiationLimits,
    started: &[[f64; 3]],
) -> (f64, f64, Vec<[f64; 3]>) {
    let sub = pf_contracts::terrain::BoundaryClassCodeV0::Subduction.code();
    let mut trench: Vec<[f64; 3]> = (0..state.cells()).filter(|&c| state.boundary_class[c] == sub).map(|c| geom.cents[c]).collect();
    trench.extend_from_slice(started);
    let reach = lim.trench_spacing_m / geom.radius;
    let cos_reach = reach.min(PI).cos();
    let slab_age = |c: usize| if state.kind[c] == CrustKindCodeV0::Oceanic { state.age_my(c) } else { 0.0 };
    let (mut total, mut old, mut near) = (0.0, 0.0, 0.0);
    let mut pts: Vec<[f64; 3]> = Vec::new();
    let mut nb = Vec::with_capacity(4);
    for (a, &c) in g.cells.iter().enumerate() {
        if !side_b[a] {
            continue;
        }
        geom.nbrs4(c, &mut nb);
        let mut on_edge = false;
        for &(kc, len) in &nb {
            let b = in_group[kc];
            if b == NIL || side_b[b] {
                continue;
            }
            on_edge = true;
            total += len;
            if slab_age(c).max(slab_age(kc)) >= lim.slab_age_min_my {
                old += len;
            }
            if trench.iter().any(|&t| dot(t, geom.cents[c]) >= cos_reach) {
                near += len;
            }
        }
        if on_edge {
            pts.push(geom.cents[c]);
        }
    }
    if total <= 0.0 {
        return (0.0, 1.0, pts);
    }
    (old / total, near / total, pts)
}

/// SEED-OLD-OCEAN (lane OR) diagnostic: ridge-class edges between suture
/// groups and the part converging under `plates` (see
/// `ReorgStepReport::resolve_ridge_m`).
fn ridge_flip_measure(geom: &Geom, state: &CrustState, plates: &[PlateRecordV0], sut: &SutureState) -> Result<(f64, f64), String> {
    let ridge = pf_contracts::terrain::BoundaryClassCodeV0::Ridge.code();
    let poles: BTreeMap<PlateId, EulerPoleV0> = plates.iter().filter_map(|p| p.euler_pole.map(|e| (p.plate_id, e))).collect();
    let min_rate = 5_000.0 / geom.radius; // 0.5 cm/yr in rad/My
    let tan30 = (30.0f64).to_radians().tan();
    let (mut tot, mut flip) = (0.0, 0.0);
    let mut nb4 = Vec::with_capacity(4);
    for c in 0..state.cells() {
        geom.nbrs4(c, &mut nb4);
        for &(k, e) in &nb4 {
            // Each undirected edge once.
            if k < c {
                continue;
            }
            if state.boundary_class[c] != ridge && state.boundary_class[k] != ridge {
                continue;
            }
            if sut.find(state.plate[c]) == sut.find(state.plate[k]) {
                continue;
            }
            tot += e;
            let (Some(pa), Some(pb)) = (poles.get(&state.plate[c]), poles.get(&state.plate[k])) else { continue };
            let Some(mid) = unit(add(geom.cents[c], geom.cents[k])) else { continue };
            let v = crate::boundaries::relative_velocity_at(mid, pa, pb).map_err(|e| e.to_string())?;
            let vn = normal_rate_at(mid, geom.cents[c], geom.cents[k], pa, pb).map_err(|e| e.to_string())?;
            let vt = (dot(v, v) - vn * vn).max(0.0).sqrt();
            if vn > 0.0 && norm(v) > min_rate && vn > tan30 * vt {
                flip += e;
            }
        }
    }
    Ok((tot, flip))
}

/// SEED-OLD-OCEAN (lane OR) `plate-reorg-scripted-v5`: one passive-margin
/// initiation (see `MarginInitiation`); candidates oldest first (ties:
/// lowest cell), at most `MARGIN_MAX_ATTEMPTS` tried, a candidate within
/// the spacing of an already tried seed skipped. Rejections (side B below
/// `MARGIN_MIN_SIDE_B_FRACTION`, or no frontal ocean) are logged as
/// `initiation-v5-margin-rejected`.
#[allow(clippy::too_many_arguments)]
fn margin_initiation(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    gs: &Groups,
    t_my: f64,
    mi: &MarginInitiation,
    started: &mut Vec<[f64; 3]>,
    tried_groups: &mut Vec<u32>,
    rep: &mut ReorgStepReport,
) -> Result<Option<ReorgEvent>, String> {
    let p = ctx.params;
    let n = state.cells();
    let sub = pf_contracts::terrain::BoundaryClassCodeV0::Subduction.code();
    let mut trench: Vec<[f64; 3]> = (0..n).filter(|&c| state.boundary_class[c] == sub).map(|c| geom.cents[c]).collect();
    trench.extend_from_slice(started);
    let cos_sp = (mi.spacing_m / geom.radius).min(PI).cos();
    let near_trench = |x: [f64; 3], tr: &[[f64; 3]]| tr.iter().any(|&t| dot(t, x) >= cos_sp);
    let is_cont = |c: usize| state.kind[c] == CrustKindCodeV0::Continental;
    let mut nb = Vec::with_capacity(8);
    // Margin cells: oceanic, with a continental 8-neighbour of the same group.
    let mut margin = vec![false; n];
    let mut cands: Vec<(f64, usize)> = Vec::new();
    for c in 0..n {
        if is_cont(c) {
            continue;
        }
        let gi = gs.of_cell[c];
        geom.nbrs8(c, &mut nb);
        let mut coast = false;
        let mut other = false;
        for &(k, _) in &nb {
            if gs.of_cell[k] != gi {
                other = true;
            } else if is_cont(k) {
                coast = true;
            }
        }
        if !coast {
            continue;
        }
        margin[c] = true;
        let age = state.age_my(c);
        if !other && age >= mi.age_my && !tried_groups.contains(&gs.list[gi].root.0) {
            cands.push((age, c));
        }
    }
    // Revision 1: rank by the group's old ocean within MARGIN_RANK_RADIUS_M.
    let cos_rank = (MARGIN_RANK_RADIUS_M / geom.radius).min(PI).cos();
    let mut ranked: Vec<(f64, f64, usize)> = cands
        .iter()
        .map(|&(age, c)| {
            let gi = gs.of_cell[c];
            let old: f64 = gs.list[gi]
                .cells
                .iter()
                .filter(|&&k| !is_cont(k) && state.age_my(k) >= mi.age_my && dot(geom.cents[k], geom.cents[c]) >= cos_rank)
                .map(|&k| geom.areas[k])
                .sum();
            (old, age, c)
        })
        .collect();
    ranked.sort_by(|a, b| b.0.total_cmp(&a.0).then(b.1.total_cmp(&a.1)).then(a.2.cmp(&b.2)));
    let cands: Vec<(f64, usize)> = ranked.into_iter().map(|(_, age, c)| (age, c)).collect();
    let total_area = gs.total_area;
    let mut attempts = 0;
    let mut tried_seeds: Vec<[f64; 3]> = Vec::new();
    for &(seed_age, seed) in &cands {
        if attempts >= MARGIN_MAX_ATTEMPTS {
            break;
        }
        let x = geom.cents[seed];
        if near_trench(x, &trench) || near_trench(x, &tried_seeds) {
            continue;
        }
        attempts += 1;
        tried_seeds.push(x);
        let gi = gs.of_cell[seed];
        let g = &gs.list[gi];
        let in_group = local_index(n, &g.cells);
        let m = g.cells.len();
        // Segment: margin cells along the margin from the seed (8-neighbour
        // Dijkstra over margin cells), within half the segment length.
        let seg_ok = |c: usize| {
            margin[c] && gs.of_cell[c] == gi && state.age_my(c) >= MARGIN_SEGMENT_MIN_AGE_MY && !near_trench(geom.cents[c], &trench)
        };
        let half = 0.5 * mi.segment_length_m;
        let mut dist: BTreeMap<usize, f64> = BTreeMap::new();
        let mut heap = std::collections::BinaryHeap::new();
        dist.insert(seed, 0.0);
        heap.push(std::cmp::Reverse((0u64, seed)));
        while let Some(std::cmp::Reverse((dk, c))) = heap.pop() {
            let d = f64::from_bits(dk);
            if d > dist[&c] {
                continue;
            }
            geom.nbrs8(c, &mut nb);
            for &(k, _) in &nb {
                if !seg_ok(k) {
                    continue;
                }
                let nd = d + dot(geom.cents[c], geom.cents[k]).clamp(-1.0, 1.0).acos() * geom.radius;
                if nd <= half && dist.get(&k).map_or(true, |&o| nd < o) {
                    dist.insert(k, nd);
                    heap.push(std::cmp::Reverse((nd.to_bits(), k)));
                }
            }
        }
        let seg: Vec<usize> = dist.keys().cloned().collect();
        // Coast labels: 1 = continental cells touching the segment.
        let mut label = vec![u8::MAX; m];
        let mut dd = vec![f64::INFINITY; m];
        let mut heap = std::collections::BinaryHeap::new();
        for a in 0..m {
            let c = g.cells[a];
            if !is_cont(c) {
                continue;
            }
            geom.nbrs8(c, &mut nb);
            let l = if nb.iter().any(|&(k, _)| dist.contains_key(&k)) { 1u8 } else { 0u8 };
            label[a] = l;
            dd[a] = 0.0;
            // Ties go to the segment coast (key 1 - label first).
            heap.push(std::cmp::Reverse((0u64, 1 - l, a)));
        }
        while let Some(std::cmp::Reverse((dk, inv, a))) = heap.pop() {
            let d = f64::from_bits(dk);
            if d > dd[a] || label[a] != 1 - inv {
                continue;
            }
            let c = g.cells[a];
            geom.nbrs8(c, &mut nb);
            for &(k, _) in &nb {
                let b = in_group[k];
                if b == NIL || is_cont(k) {
                    continue;
                }
                let nd = d + dot(geom.cents[c], geom.cents[k]).clamp(-1.0, 1.0).acos() * geom.radius;
                let l = 1 - inv;
                if nd < dd[b] || (nd == dd[b] && l > label[b]) {
                    dd[b] = nd;
                    label[b] = l;
                    heap.push(std::cmp::Reverse((nd.to_bits(), inv, b)));
                }
            }
        }
        for &c in &seg {
            label[in_group[c]] = 1;
        }
        // Side B: the frontal ocean (label 1) connected to the seed.
        let mut side_b = vec![false; m];
        let mut stack = vec![in_group[seed]];
        side_b[in_group[seed]] = true;
        while let Some(a) = stack.pop() {
            geom.nbrs8(g.cells[a], &mut nb);
            for &(k, _) in &nb {
                let b = in_group[k];
                if b != NIL && !side_b[b] && !is_cont(k) && label[b] == 1 {
                    side_b[b] = true;
                    stack.push(b);
                }
            }
        }
        let area_b: f64 = (0..m).filter(|&a| side_b[a]).map(|a| geom.areas[g.cells[a]]).sum();
        let coast_cells: Vec<usize> = (0..m).filter(|&a| is_cont(g.cells[a]) && label[a] == 1).map(|a| g.cells[a]).collect();
        // The new boundary: 4-neighbour A-B edges inside the group; the
        // trench part runs against continental cells.
        let (mut cut_len, mut trench_len) = (0.0, 0.0);
        let mut nb4 = Vec::with_capacity(4);
        for a in 0..m {
            if !side_b[a] {
                continue;
            }
            geom.nbrs4(g.cells[a], &mut nb4);
            for &(k, e) in &nb4 {
                let b = in_group[k];
                if b != NIL && !side_b[b] {
                    cut_len += e;
                    if is_cont(k) {
                        trench_len += e;
                    }
                }
            }
        }
        let reject = |rs: &mut ReorgState, rep: &mut ReorgStepReport| {
            rep.margin_rejected += 1;
            rs.events.push(ReorgEvent { t_my, kind: "initiation-v5-margin-rejected".into(), root_a: g.root.0, root_b: 0, cells_moved: 0, group_share: area_b / total_area, side_b_share: 0.0, cut_length_m: cut_len, weak_length_m: trench_len, groups_counted: 0 });
        };
        if area_b < MARGIN_MIN_SIDE_B_FRACTION * total_area || trench_len <= 0.0 || area_b >= g.ocean_area + g.cont_area {
            reject(rs, rep);
            continue;
        }
        let ps = centroid(geom, seg.iter().cloned());
        let pc = centroid(geom, coast_cells.iter().cloned());
        let pb = centroid(geom, (0..m).filter(|&a| side_b[a]).map(|a| g.cells[a]));
        let (Some(ps), Some(pc), Some(pb)) = (ps, pc, pb) else {
            reject(rs, rep);
            continue;
        };
        let nrm = unit(tangent([ps[0] - pc[0], ps[1] - pc[1], ps[2] - pc[2]], ps))
            .or_else(|| unit(tangent([pb[0] - pc[0], pb[1] - pc[1], pb[2] - pc[2]], ps)));
        let Some(nrm) = nrm else {
            reject(rs, rep);
            continue;
        };
        let Some(ax) = unit(cross(ps, nrm)) else {
            reject(rs, rep);
            continue;
        };
        // B moves toward the continent (velocity -rate x nrm at ps).
        let omega_rel = scale(ax, -p.subduction_initiation_rate_deg_per_my.to_radians());
        let out = apply_split(geom, ctx.grid, state, plates, layer, sut, g, &side_b, omega_rel, t_my)?;
        rs.relations.push(Relation { parent: out.root_a.0, child: out.root_b.0, omega_rel, t_my, kind: 2 });
        started.extend(seg.iter().map(|&c| geom.cents[c]));
        tried_groups.push(g.root.0);
        rep.margin_seeds.push((seed_age, seed as u32, area_b, trench_len));
        return Ok(Some(ReorgEvent {
            t_my,
            kind: "initiation-v5-margin".into(),
            root_a: out.root_a.0,
            root_b: out.root_b.0,
            cells_moved: out.cells_moved,
            group_share: g.cont_area / gs.total_cont.max(1.0),
            side_b_share: 0.0,
            cut_length_m: cut_len,
            weak_length_m: trench_len,
            groups_counted: 0,
        }));
    }
    Ok(None)
}

/// Secondary pieces of each plate id: small oceanic ones fold into the
/// neighbour they touch most, the rest become new ids in the same group.
fn fragments(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    total_area: f64,
    t_my: f64,
) -> Result<u32, String> {
    let n = state.cells();
    let min_a = ctx.params.plate_count_min_area_fraction * total_area;
    let mut comp = vec![u32::MAX; n];
    let mut pieces: BTreeMap<PlateId, Vec<(f64, bool, Vec<usize>)>> = BTreeMap::new();
    let mut nb = Vec::with_capacity(8);
    for s0 in 0..n {
        if comp[s0] != u32::MAX {
            continue;
        }
        let p = state.plate[s0];
        let mut cells = vec![s0];
        comp[s0] = 0;
        let mut h = 0;
        while h < cells.len() {
            let c = cells[h];
            h += 1;
            geom.nbrs8(c, &mut nb);
            for &(k, _) in &nb {
                if comp[k] == u32::MAX && state.plate[k] == p {
                    comp[k] = 0;
                    cells.push(k);
                }
            }
        }
        cells.sort_unstable();
        let area: f64 = cells.iter().map(|&c| geom.areas[c]).sum();
        let cont = cells.iter().any(|&c| state.kind[c] == CrustKindCodeV0::Continental);
        pieces.entry(p).or_default().push((area, cont, cells));
    }
    let mut count = 0;
    for (p, mut v) in pieces {
        if v.len() < 2 {
            continue;
        }
        // Largest piece keeps the id (ties: lowest first cell).
        v.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.2[0].cmp(&b.2[0])));
        for (area, cont, cells) in v.into_iter().skip(1) {
            if area < min_a && !cont {
                if let Some(q) = most_touching(geom, state, &cells, |x| x == p) {
                    relabel(ctx.grid, state, layer, &cells, q)?;
                    count += 1;
                    rs.events.push(ReorgEvent { t_my, kind: "fragment-fold".into(), root_a: p.0, root_b: q.0, cells_moved: cells.len() as u64, group_share: 0.0, side_b_share: 0.0, cut_length_m: 0.0, weak_length_m: 0.0, groups_counted: 0 });
                }
            } else {
                let id = PlateId(plates.iter().map(|r| r.plate_id.0).max().unwrap_or(0) + 1);
                let pole = plates.iter().find(|r| r.plate_id == p).and_then(|r| r.euler_pole);
                plates.push(PlateRecordV0 { plate_id: id, parent_plate: Some(p), birth_t_my: t_my, euler_pole: pole, retirement: None });
                let root = sut.find(p);
                sut.parent.insert(id, root);
                relabel(ctx.grid, state, layer, &cells, id)?;
                count += 1;
                rs.events.push(ReorgEvent { t_my, kind: "fragment-plate".into(), root_a: p.0, root_b: id.0, cells_moved: cells.len() as u64, group_share: 0.0, side_b_share: 0.0, cut_length_m: 0.0, weak_length_m: 0.0, groups_counted: 0 });
            }
        }
    }
    Ok(count)
}

/// Slab torque (unit-sphere x t, weighted by edge length m), subducting
/// edge length and total boundary length of group `gi`.
fn slab_torque(geom: &Geom, state: &CrustState, plates: &[PlateRecordV0], gs: &Groups, gi: usize) -> Result<([f64; 3], f64, f64), String> {
    let poles: BTreeMap<PlateId, EulerPoleV0> = plates.iter().filter_map(|p| p.euler_pole.map(|e| (p.plate_id, e))).collect();
    let (mut tau, mut sub, mut per) = ([0.0; 3], 0.0, 0.0);
    let mut nb4 = Vec::with_capacity(4);
    for &c in &gs.list[gi].cells {
        geom.nbrs4(c, &mut nb4);
        for &(k, e) in &nb4 {
            if gs.of_cell[k] == gi {
                continue;
            }
            per += e;
            if state.kind[c] != CrustKindCodeV0::Oceanic {
                continue;
            }
            let (Some(pa), Some(pb)) = (poles.get(&state.plate[c]), poles.get(&state.plate[k])) else { continue };
            let Some(mid) = unit(add(geom.cents[c], geom.cents[k])) else { continue };
            let vn = normal_rate_at(mid, geom.cents[c], geom.cents[k], pa, pb).map_err(|e| e.to_string())?;
            if !(vn > 0.0) {
                continue;
            }
            let down = state.kind[k] == CrustKindCodeV0::Continental || state.birth_t_my[k] > state.birth_t_my[c];
            if !down {
                continue;
            }
            let x = geom.cents[c];
            let Some(t) = unit(tangent([geom.cents[k][0] - x[0], geom.cents[k][1] - x[1], geom.cents[k][2] - x[2]], x)) else { continue };
            tau = add(tau, scale(cross(x, t), e));
            sub += e;
        }
    }
    Ok((tau, sub, per))
}

/// Re-draw every group's pole (ascending root), then re-impose young
/// relations.
fn resolve(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &CrustState,
    plates: &mut [PlateRecordV0],
    sut: &SutureState,
    rs: &mut ReorgState,
    k: u32,
    t_my: f64,
    ridge_aware: bool,
) -> Result<(), String> {
    let p = ctx.params;
    let gs = groups(geom, state, plates, sut);
    let mut rng = SplitMix64::new(seed_from_path_tagged(ctx.seed_path, &format!("plate-reorg-scripted-v1/resolve/{k}")));
    // SEED-OLD-OCEAN (lane OR) v5: ridge-aware axis choice. Each group's
    // ridge edges against other groups (cell in the group, other cell,
    // edge length, other group index), and every group's current omega.
    let mut ridge_edges: Vec<Vec<(usize, usize, f64, usize)>> = vec![Vec::new(); gs.list.len()];
    let mut cur: Vec<Option<[f64; 3]>> = vec![None; gs.list.len()];
    if ridge_aware {
        let ridge = pf_contracts::terrain::BoundaryClassCodeV0::Ridge.code();
        let mut nb4 = Vec::with_capacity(4);
        for c in 0..state.cells() {
            let gc = gs.of_cell[c];
            geom.nbrs4(c, &mut nb4);
            for &(kc, e) in &nb4 {
                let gk = gs.of_cell[kc];
                let young = |x: usize| state.kind[x] == CrustKindCodeV0::Oceanic && state.age_my(x) < RESOLVE_YOUNG_EDGE_MY;
                if gk != gc && (state.boundary_class[c] == ridge || state.boundary_class[kc] == ridge || (young(c) && young(kc))) {
                    ridge_edges[gc].push((c, kc, e, gk));
                }
            }
        }
        for (gi, g) in gs.list.iter().enumerate() {
            cur[gi] = omega_of(plates, g.root).ok();
        }
    }
    let min_rate = 5_000.0 / geom.radius;
    let tan30 = (30.0f64).to_radians().tan();
    let closing = |gi: usize, w: [f64; 3], cur: &[Option<[f64; 3]>]| -> f64 {
        let mut l = 0.0;
        for &(c, kc, e, gk) in &ridge_edges[gi] {
            let Some(wk) = cur[gk] else { continue };
            let Some(mid) = unit(add(geom.cents[c], geom.cents[kc])) else { continue };
            let Some(n) = unit(tangent([geom.cents[kc][0] - geom.cents[c][0], geom.cents[kc][1] - geom.cents[c][1], geom.cents[kc][2] - geom.cents[c][2]], mid)) else { continue };
            let v = cross([w[0] - wk[0], w[1] - wk[1], w[2] - wk[2]], mid);
            let vn = dot(v, n);
            let vt = (dot(v, v) - vn * vn).max(0.0).sqrt();
            if vn > 0.0 && norm(v) > min_rate && vn > tan30 * vt {
                l += e;
            }
        }
        l
    };
    let mut omega: BTreeMap<PlateId, [f64; 3]> = BTreeMap::new();
    let mut cands_of: Vec<Vec<[f64; 3]>> = vec![Vec::new(); gs.list.len()];
    for (gi, g) in gs.list.iter().enumerate() {
        let u_axis = unit_dir(&mut rng);
        let u_rate = rng.next_f64();
        if g.cells.is_empty() {
            continue;
        }
        let (tau, sub, per) = slab_torque(geom, state, plates, &gs, gi)?;
        let blend = |u: [f64; 3]| match unit(tau) {
            Some(tu) => unit(add(scale(u, 1.0 - p.slab_pull_bias), scale(tu, p.slab_pull_bias))).unwrap_or(u),
            None => u,
        };
        let f_sub = if per > 0.0 { sub / per } else { 0.0 };
        let rate = p.rate_scale_deg_per_my
            * (p.rate_floor_fraction + (1.0 - p.rate_floor_fraction) * (f_sub / p.slab_circumference_reference).min(1.0))
            * (0.75 + 0.5 * u_rate);
        let w = scale(blend(u_axis), rate.to_radians());
        if ridge_aware && !ridge_edges[gi].is_empty() {
            // Candidate 0 is the v4 draw; RESOLVE_AXIS_CANDIDATES - 1 more
            // axes from a tagged stream; scored in the sweeps below.
            let mut r2 = SplitMix64::new(seed_from_path_tagged(ctx.seed_path, &format!("plate-reorg-scripted-v5/resolve/{k}/{}", g.root.0)));
            let mut cs = vec![w];
            for _ in 1..RESOLVE_AXIS_CANDIDATES {
                cs.push(scale(blend(unit_dir(&mut r2)), rate.to_radians()));
            }
            cands_of[gi] = cs;
        }
        omega.insert(g.root, w);
    }
    // SEED-OLD-OCEAN (lane OR) v5: least closing ridge / young-young length
    // wins (ties: lowest candidate), ascending group order, RESOLVE_SWEEPS
    // sweeps (the later ones see every neighbour's new pole).
    if ridge_aware {
        for _ in 0..RESOLVE_SWEEPS {
            for (gi, g) in gs.list.iter().enumerate() {
                if cands_of[gi].is_empty() {
                    continue;
                }
                let mut best = (f64::INFINITY, cands_of[gi][0]);
                for &wc in &cands_of[gi] {
                    let l = closing(gi, wc, &cur);
                    if l < best.0 {
                        best = (l, wc);
                    }
                }
                cur[gi] = Some(best.1);
                omega.insert(g.root, best.1);
            }
        }
    }
    rs.relations.retain(|r| t_my - r.t_my < p.rift_persistence_my);
    for r in &rs.relations {
        let (ra, rb) = (sut.find(PlateId(r.parent)), sut.find(PlateId(r.child)));
        if ra == rb {
            continue;
        }
        if let (Some(&wa), true) = (omega.get(&ra), omega.contains_key(&rb)) {
            omega.insert(rb, add(wa, r.omega_rel));
        }
    }
    for rec in plates.iter_mut() {
        if let Some(w) = omega.get(&sut.find(rec.plate_id)) {
            rec.euler_pole = Some(pole_of(*w)?);
        }
    }
    Ok(())
}

/// Every step: a group whose root plate has no cells is re-rooted on its
/// largest member (same pole), then plates with no cells and no parcels
/// retire (pole None; they claim nothing, so this is motion-neutral and
/// keeps the carry's per-plate cost bounded).
fn housekeeping(
    state: &CrustState,
    plates: &mut [PlateRecordV0],
    layer: &ParcelLayer,
    sut: &mut SutureState,
    areas: &[f64],
    t_my: f64,
    cause: PlateRetirementV0,
) -> Result<u32, String> {
    let mut area: BTreeMap<PlateId, f64> = BTreeMap::new();
    for c in 0..state.cells() {
        *area.entry(state.plate[c]).or_insert(0.0) += areas[c];
    }
    let mut parcels: BTreeMap<PlateId, u64> = BTreeMap::new();
    for p in &layer.plate {
        *parcels.entry(*p).or_insert(0) += 1;
    }
    // Re-root groups whose root has no cells.
    let mut groups: BTreeMap<PlateId, Vec<PlateId>> = BTreeMap::new();
    for r in plates.iter() {
        groups.entry(sut.find(r.plate_id)).or_default().push(r.plate_id);
    }
    for (root, members) in &groups {
        if area.contains_key(root) {
            continue;
        }
        let mut best: Option<(f64, PlateId)> = None;
        for m in members {
            if let Some(&a) = area.get(m) {
                if best.map_or(true, |(ba, bid)| a > ba || (a == ba && *m < bid)) {
                    best = Some((a, *m));
                }
            }
        }
        let Some((_, nr)) = best else { continue };
        let pole = plates.iter().find(|r| r.plate_id == *root).and_then(|r| r.euler_pole);
        for m in members {
            sut.parent.remove(m);
            if *m != nr {
                sut.parent.insert(*m, nr);
            }
        }
        if let Some(pole) = pole {
            for r in plates.iter_mut() {
                if members.contains(&r.plate_id) {
                    r.euler_pole = Some(pole);
                }
            }
        }
    }
    let mut retired = 0;
    for r in plates.iter_mut() {
        if r.euler_pole.is_some() && !area.contains_key(&r.plate_id) && !parcels.contains_key(&r.plate_id) {
            let is_root = sut.parent.values().any(|q| *q == r.plate_id);
            if is_root {
                continue;
            }
            r.euler_pole = None;
            r.retirement = Some(pf_contracts::domain::PlateRetirementDetailV0 {
                t_my,
                cause,
            });
            sut.parent.remove(&r.plate_id);
            retired += 1;
        }
    }
    Ok(retired)
}

/// At each scheduled event: member plates of a group that hold no
/// continent fold into the group-mate they touch most (same pole, so the
/// motion is unchanged; fewer plate ids).
fn consolidate(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &mut CrustState,
    plates: &[PlateRecordV0],
    layer: &mut ParcelLayer,
    sut: &SutureState,
) -> Result<u32, String> {
    let gs = groups(geom, state, plates, sut);
    let mut count = 0;
    for g in &gs.list {
        if g.members.len() < 2 {
            continue;
        }
        let mut cells_of: BTreeMap<PlateId, (Vec<usize>, bool)> = BTreeMap::new();
        for &c in &g.cells {
            let e = cells_of.entry(state.plate[c]).or_insert((Vec::new(), false));
            e.0.push(c);
            e.1 |= state.kind[c] == CrustKindCodeV0::Continental;
        }
        if cells_of.len() < 2 {
            continue;
        }
        for (m, (cells, cont)) in &cells_of {
            if *cont {
                continue;
            }
            let root = g.root;
            let pick = most_touching(geom, state, cells, |x| x == *m || sut.find(x) != root);
            if let Some(q) = pick {
                relabel(ctx.grid, state, layer, cells, q)?;
                count += 1;
            }
        }
    }
    Ok(count)
}

/// One step's pass, at the start of the step (t0, t1]: scheduled events
/// (fragments, merges, re-solve, initiations) then at most one breakup.
#[allow(clippy::too_many_arguments)]
pub fn reorg_pass(
    ctx: &ReorgCtx<'_>,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    t0: f64,
    t1: f64,
) -> Result<ReorgStepReport, String> {
    reorg_pass_law(ctx, ReorgLaw::ScriptedV1, state, plates, layer, sut, rs, t0, t1)
}

/// `reorg_pass` for a chosen law version (lane BC): v2 breaks up along
/// `breakup_v2` and retires plates emptied by folds as `Folded`.
#[allow(clippy::too_many_arguments)]
pub fn reorg_pass_law(
    ctx: &ReorgCtx<'_>,
    law: ReorgLaw,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    t0: f64,
    t1: f64,
) -> Result<ReorgStepReport, String> {
    if law == ReorgLaw::ScriptedV3 {
        return Err(format!("{LAW_V3}: run through reorg_pass_v3 (it needs the initiation age)"));
    }
    if law == ReorgLaw::ScriptedV4 {
        return Err(format!("{LAW_V4}: run through reorg_pass_v4 (it needs the initiation age and limits)"));
    }
    if law == ReorgLaw::ScriptedV5 {
        return Err(format!("{LAW_V5}: run through reorg_pass_v5 (it needs the initiation age, limits and margin initiation)"));
    }
    reorg_pass_impl(ctx, law, None, None, None, false, state, plates, layer, sut, rs, t0, t1)
}

/// SEED-WORLD-COMPOSE (lane CP): `plate-reorg-scripted-v3`, with its
/// subduction initiation age (My).
#[allow(clippy::too_many_arguments)]
pub fn reorg_pass_v3(
    ctx: &ReorgCtx<'_>,
    initiation_age_my: f64,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    t0: f64,
    t1: f64,
) -> Result<ReorgStepReport, String> {
    if !(initiation_age_my.is_finite() && initiation_age_my > 0.0) {
        return Err(format!("{LAW_V3}: initiation age must be finite and positive, got {initiation_age_my}"));
    }
    reorg_pass_impl(ctx, ReorgLaw::ScriptedV3, Some(initiation_age_my), None, None, false, state, plates, layer, sut, rs, t0, t1)
}

/// SEED-TRANSFORMS (lane TF): `plate-reorg-scripted-v4` (v3 with the
/// age-triggered initiation limited by `InitiationLimits`).
#[allow(clippy::too_many_arguments)]
pub fn reorg_pass_v4(
    ctx: &ReorgCtx<'_>,
    initiation_age_my: f64,
    limits: &InitiationLimits,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    t0: f64,
    t1: f64,
) -> Result<ReorgStepReport, String> {
    if !(initiation_age_my.is_finite() && initiation_age_my > 0.0) {
        return Err(format!("{LAW_V4}: initiation age must be finite and positive, got {initiation_age_my}"));
    }
    limits.validate()?;
    reorg_pass_impl(ctx, ReorgLaw::ScriptedV4, Some(initiation_age_my), Some(limits), None, false, state, plates, layer, sut, rs, t0, t1)
}

/// SEED-OLD-OCEAN (lane OR): `plate-reorg-scripted-v5` (v4 plus the
/// passive-margin initiation `MarginInitiation`, before v4's age trigger at
/// each scheduled event).
#[allow(clippy::too_many_arguments)]
pub fn reorg_pass_v5(
    ctx: &ReorgCtx<'_>,
    initiation_age_my: f64,
    limits: &InitiationLimits,
    margin: &MarginInitiation,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    t0: f64,
    t1: f64,
) -> Result<ReorgStepReport, String> {
    if !(initiation_age_my.is_finite() && initiation_age_my > 0.0) {
        return Err(format!("{LAW_V5}: initiation age must be finite and positive, got {initiation_age_my}"));
    }
    limits.validate()?;
    margin.validate()?;
    reorg_pass_impl(ctx, ReorgLaw::ScriptedV5, Some(initiation_age_my), Some(limits), Some(margin), false, state, plates, layer, sut, rs, t0, t1)
}

/// PLATE-POLES (lane PLATE-POLES, 2026-10-21): any reorganization law
/// version with `rigid-group-contacts-v1` part 2 (`release`; see
/// `boundaries::LAW_RIGID_GROUP_CONTACTS_V1`): at each scheduled event,
/// after the fragment pass and before the re-solve, suture-group members
/// detached from the group's main body leave it. `release = false` is bit
/// for bit the matching `reorg_pass_*` entry point. v3 needs the
/// initiation age, v4 also the limits, v5 also the margin initiation.
#[allow(clippy::too_many_arguments)]
pub fn reorg_pass_release(
    ctx: &ReorgCtx<'_>,
    law: ReorgLaw,
    initiation_age_my: Option<f64>,
    limits: Option<&InitiationLimits>,
    margin: Option<&MarginInitiation>,
    release: bool,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    t0: f64,
    t1: f64,
) -> Result<ReorgStepReport, String> {
    let need_age = matches!(law, ReorgLaw::ScriptedV3 | ReorgLaw::ScriptedV4 | ReorgLaw::ScriptedV5);
    if need_age != initiation_age_my.is_some()
        || matches!(law, ReorgLaw::ScriptedV4 | ReorgLaw::ScriptedV5) != limits.is_some()
        || (law == ReorgLaw::ScriptedV5) != margin.is_some()
    {
        return Err(format!("{}: reorganization law and its inputs do not match", crate::boundaries::LAW_RIGID_GROUP_CONTACTS_V1));
    }
    if let Some(age) = initiation_age_my {
        if !(age.is_finite() && age > 0.0) {
            return Err(format!("plate reorganization: initiation age must be finite and positive, got {age}"));
        }
    }
    if let Some(l) = limits {
        l.validate()?;
    }
    if let Some(m) = margin {
        m.validate()?;
    }
    reorg_pass_impl(ctx, law, initiation_age_my, limits, margin, release, state, plates, layer, sut, rs, t0, t1)
}

/// PLATE-POLES `rigid-group-contacts-v1` part 2: members of a suture group
/// that touch (8-neighbour) no other member of the group's main body leave
/// it. Members present on the grid are joined by contact; the connected
/// set with the largest area (ties: lowest member id) keeps the group, and
/// members with no cells stay with it. Every other set becomes its own
/// group rooted at its largest member (ties: lowest id). If the old root
/// leaves, the main body is re-rooted on its largest member and the
/// tenure / cascade timers and the young relations that named a leaving
/// plate follow the main body. Poles are not touched here (the event's
/// re-solve draws them). Returns the number of sets released.
fn release_detached(
    geom: &Geom,
    state: &CrustState,
    plates: &[PlateRecordV0],
    sut: &mut SutureState,
    rs: &mut ReorgState,
    t_my: f64,
) -> u32 {
    let gs = groups(geom, state, plates, sut);
    let mut released = 0u32;
    let mut nb = Vec::with_capacity(8);
    for (gi, g) in gs.list.iter().enumerate() {
        if g.members.len() < 2 || g.cells.is_empty() {
            continue;
        }
        // Union-find over the member ids present on the grid.
        let mut area: BTreeMap<PlateId, f64> = BTreeMap::new();
        let mut cells_of: BTreeMap<PlateId, u64> = BTreeMap::new();
        for &c in &g.cells {
            *area.entry(state.plate[c]).or_insert(0.0) += geom.areas[c];
            *cells_of.entry(state.plate[c]).or_insert(0) += 1;
        }
        let mut uf: BTreeMap<PlateId, PlateId> = area.keys().map(|p| (*p, *p)).collect();
        fn root_of(uf: &BTreeMap<PlateId, PlateId>, mut p: PlateId) -> PlateId {
            while uf[&p] != p {
                p = uf[&p];
            }
            p
        }
        for &c in &g.cells {
            geom.nbrs8(c, &mut nb);
            for &(k, _) in &nb {
                if gs.of_cell[k] != gi || state.plate[k] == state.plate[c] {
                    continue;
                }
                let (ra, rb) = (root_of(&uf, state.plate[c]), root_of(&uf, state.plate[k]));
                if ra != rb {
                    let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
                    uf.insert(hi, lo);
                }
            }
        }
        let mut sets: BTreeMap<PlateId, (f64, Vec<PlateId>)> = BTreeMap::new();
        for (p, a) in &area {
            let e = sets.entry(root_of(&uf, *p)).or_insert((0.0, Vec::new()));
            e.0 += a;
            e.1.push(*p);
        }
        if sets.len() < 2 {
            continue;
        }
        // Main body: largest area, ties lowest member id (= the set key).
        let mut main: Option<(f64, PlateId)> = None;
        for (k, (a, _)) in &sets {
            if main.map_or(true, |(ba, _)| *a > ba) {
                main = Some((*a, *k));
            }
        }
        let main_key = main.expect("two sets").1;
        let largest = |ms: &[PlateId]| -> PlateId {
            let mut b = ms[0];
            for m in ms {
                if area[m] > area[&b] || (area[m] == area[&b] && *m < b) {
                    b = *m;
                }
            }
            b
        };
        let leaving: Vec<PlateId> = sets.iter().filter(|(k, _)| **k != main_key).flat_map(|(_, (_, ms))| ms.iter().copied()).collect();
        let main_root = if leaving.contains(&g.root) { largest(&sets[&main_key].1) } else { g.root };
        // Rewrite the group's union-find: main body (and cell-less members)
        // under main_root, each leaving set under its own largest member.
        for m in &g.members {
            sut.parent.remove(m);
        }
        for m in &g.members {
            if !leaving.contains(m) && *m != main_root {
                sut.parent.insert(*m, main_root);
            }
        }
        for (k, (_, ms)) in &sets {
            if *k == main_key {
                continue;
            }
            let nr = largest(ms);
            for m in ms {
                if *m != nr {
                    sut.parent.insert(*m, nr);
                }
            }
            released += 1;
            let moved: u64 = ms.iter().map(|m| cells_of[m]).sum();
            rs.events.push(ReorgEvent { t_my, kind: "release".into(), root_a: main_root.0, root_b: nr.0, cells_moved: moved, group_share: 0.0, side_b_share: 0.0, cut_length_m: 0.0, weak_length_m: 0.0, groups_counted: 0 });
        }
        // Timers and young relations follow the main body.
        let rekey = |x: &mut u32| {
            if *x == g.root.0 || leaving.iter().any(|m| m.0 == *x) {
                *x = main_root.0;
            }
        };
        for t in rs.tenure.iter_mut() {
            rekey(&mut t.0);
        }
        for t in rs.cascade.iter_mut() {
            rekey(&mut t.0);
        }
        for r in rs.relations.iter_mut() {
            rekey(&mut r.parent);
            rekey(&mut r.child);
        }
    }
    released
}

#[allow(clippy::too_many_arguments)]
fn reorg_pass_impl(
    ctx: &ReorgCtx<'_>,
    law: ReorgLaw,
    initiation_age_my: Option<f64>,
    limits: Option<&InitiationLimits>,
    margin: Option<&MarginInitiation>,
    release: bool,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    t0: f64,
    t1: f64,
) -> Result<ReorgStepReport, String> {
    ctx.params.validate()?;
    let geom = Geom::new(ctx.grid)?;
    let p = ctx.params;
    let mut rep = ReorgStepReport::default();
    // SEED-TRANSFORMS (lane TF) v4: boundaries started by age-triggered
    // initiations earlier in this pass (side-B edge cell centres).
    let mut started: Vec<[f64; 3]> = Vec::new();
    rep.retired = housekeeping(state, plates, layer, sut, &geom.areas, t0, PlateRetirementV0::Subducted)?;
    while (rs.next_event as usize) < ctx.schedule.len() && ctx.schedule[rs.next_event as usize] <= t1 + 1e-9 {
        let k = rs.next_event;
        rs.next_event += 1;
        // 1. Leftover pieces; oceanic group members fold into group-mates.
        let ta: f64 = geom.areas.iter().sum();
        rep.fragments += fragments(ctx, &geom, state, plates, layer, sut, rs, ta, t0)?;
        rep.fragments += consolidate(ctx, &geom, state, plates, layer, sut)?;
        if law != ReorgLaw::ScriptedV1 {
            rep.retired += housekeeping(state, plates, layer, sut, &geom.areas, t0, PlateRetirementV0::Folded)?;
        }
        // 1b. PLATE-POLES `rigid-group-contacts-v1`: detached members leave.
        if release {
            rep.released += release_detached(&geom, state, plates, sut, rs, t0);
        }
        let merges_before = rep.merges;
        // 2. Too many groups: fold the smallest oceanic one.
        loop {
            let gs = groups(&geom, state, plates, sut);
            if counted(&gs, p.plate_count_min_area_fraction) <= p.plate_count_max {
                break;
            }
            let min_a = p.plate_count_min_area_fraction * gs.total_area;
            let mut pick: Option<usize> = None;
            for (i, g) in gs.list.iter().enumerate() {
                if g.area >= min_a && g.cont_area == 0.0 && pick.map_or(true, |b| g.area < gs.list[b].area) {
                    pick = Some(i);
                }
            }
            let Some(i) = pick else { break };
            let g = &gs.list[i];
            let gi = i;
            let Some(q) = most_touching(&geom, state, &g.cells, |x| gs.of_cell_plate(&x, gi, sut)) else { break };
            relabel(ctx.grid, state, layer, &g.cells, q)?;
            rep.merges += 1;
            rs.events.push(ReorgEvent { t_my: t0, kind: "merge".into(), root_a: g.root.0, root_b: q.0, cells_moved: g.cells.len() as u64, group_share: 0.0, side_b_share: 0.0, cut_length_m: 0.0, weak_length_m: 0.0, groups_counted: 0 });
        }
        if law != ReorgLaw::ScriptedV1 && rep.merges > merges_before {
            rep.retired += housekeeping(state, plates, layer, sut, &geom.areas, t0, PlateRetirementV0::Folded)?;
        }
        // 3. Re-solve.
        resolve(ctx, &geom, state, plates, sut, rs, k, t0, law == ReorgLaw::ScriptedV5)?;
        rep.resolves += 1;
        // SEED-OLD-OCEAN (lane OR) log: ridges closing under the new poles.
        if matches!(law, ReorgLaw::ScriptedV4 | ReorgLaw::ScriptedV5) {
            let (r, f) = ridge_flip_measure(&geom, state, plates, sut)?;
            rep.resolve_ridge_m += r;
            rep.resolve_ridge_to_convergent_m += f;
        }
        let gs = groups(&geom, state, plates, sut);
        rs.events.push(ReorgEvent { t_my: t0, kind: "resolve".into(), root_a: k, root_b: 0, cells_moved: 0, group_share: 0.0, side_b_share: 0.0, cut_length_m: 0.0, weak_length_m: 0.0, groups_counted: counted(&gs, p.plate_count_min_area_fraction) });
        // 4. Subduction initiation while too few groups.
        for _ in 0..MAX_INITIATIONS_PER_EVENT {
            let gs = groups(&geom, state, plates, sut);
            if counted(&gs, p.plate_count_min_area_fraction) >= p.plate_count_min {
                break;
            }
            let ev = match initiation_age_my {
                None => initiation(ctx, &geom, state, plates, layer, sut, rs, &gs, t0)?,
                Some(age) => initiation_v3(ctx, &geom, state, plates, layer, sut, rs, &gs, t0, age, None, None)?,
            };
            match ev {
                Some(mut e) => {
                    e.groups_counted = counted(&groups(&geom, state, plates, sut), p.plate_count_min_area_fraction);
                    rs.events.push(e);
                    rep.initiations += 1;
                }
                None => break,
            }
        }
        // 4b. SEED-OLD-OCEAN (lane OR) v5: passive-margin initiation,
        // oldest margin first, while fewer than the maximum groups.
        if let Some(mi) = margin {
            let mut tried_groups: Vec<u32> = Vec::new();
            for _ in 0..MARGIN_MAX_PER_EVENT {
                let gs = groups(&geom, state, plates, sut);
                if counted(&gs, p.plate_count_min_area_fraction) >= p.plate_count_max {
                    break;
                }
                match margin_initiation(ctx, &geom, state, plates, layer, sut, rs, &gs, t0, mi, &mut started, &mut tried_groups, &mut rep)? {
                    Some(mut e) => {
                        e.groups_counted = counted(&groups(&geom, state, plates, sut), p.plate_count_min_area_fraction);
                        rs.events.push(e);
                        rep.initiations += 1;
                    }
                    None => break,
                }
            }
        }
        // 5. SEED-WORLD-COMPOSE (lane CP) v3: age-driven initiation. While
        // fewer than the maximum groups, the group with the most ocean
        // older than the initiation age (at least the counted-plate area)
        // starts a trench through it; each group is tried once per event.
        if let Some(age) = initiation_age_my {
            let mut tried: Vec<u32> = Vec::new();
            for _ in 0..MAX_INITIATIONS_PER_EVENT {
                let gs = groups(&geom, state, plates, sut);
                if counted(&gs, p.plate_count_min_area_fraction) >= p.plate_count_max {
                    break;
                }
                let min_a = p.plate_count_min_area_fraction * gs.total_area;
                let mut pick: Option<(usize, f64)> = None;
                for (i, g) in gs.list.iter().enumerate() {
                    if g.cells.is_empty() || tried.contains(&g.root.0) {
                        continue;
                    }
                    let old: f64 = g
                        .cells
                        .iter()
                        .filter(|&&c| state.kind[c] == CrustKindCodeV0::Oceanic && state.age_my(c) >= age)
                        .map(|&c| geom.areas[c])
                        .sum();
                    if old >= min_a && pick.map_or(true, |(_, b)| old > b) {
                        pick = Some((i, old));
                    }
                }
                let Some((gi, _)) = pick else { break };
                tried.push(gs.list[gi].root.0);
                if let Some(mut e) = initiation_v3(ctx, &geom, state, plates, layer, sut, rs, &gs, t0, age, Some(gi), limits.map(|l| (l, &mut started)))? {
                    e.groups_counted = counted(&groups(&geom, state, plates, sut), p.plate_count_min_area_fraction);
                    rs.events.push(e);
                    rep.initiations += 1;
                }
            }
        }
    }
    // Breakups (every step).
    let gs = groups(&geom, state, plates, sut);
    let tc = gs.total_cont.max(1.0);
    let roots: BTreeMap<u32, usize> = gs.list.iter().enumerate().filter(|(_, g)| !g.cells.is_empty()).map(|(i, g)| (g.root.0, i)).collect();
    rs.tenure.retain(|(r, _)| roots.get(r).map_or(false, |&i| gs.list[i].cont_area / tc >= p.supercontinent_share));
    for (r, &i) in &roots {
        if gs.list[i].cont_area / tc >= p.supercontinent_share && !rs.tenure.iter().any(|x| x.0 == *r) {
            rs.tenure.push((*r, t0));
        }
    }
    let mut due: Option<(u32, &'static str)> = None;
    if let Some(&(r, _)) = rs.tenure.iter().find(|(_, s)| t0 - s >= p.supercontinent_tenure_my - 1e-9) {
        due = Some((r, "breakup"));
    } else if let Some(&(r, _)) = rs.cascade.iter().find(|(_, d)| *d <= t0 + 1e-9) {
        rs.cascade.retain(|x| x.0 != r);
        if roots.get(&r).map_or(false, |&i| gs.list[i].cont_area / tc >= p.breakup_piece_max_share) {
            due = Some((r, "cascade"));
        }
    }
    if let Some((r, kind)) = due {
        rs.tenure.retain(|x| x.0 != r);
        let g = &gs.list[roots[&r]];
        let ev = match law {
            ReorgLaw::ScriptedV1 => breakup(ctx, &geom, state, plates, layer, sut, rs, g, tc, kind, t0)?,
            ReorgLaw::ScriptedV2 | ReorgLaw::ScriptedV3 | ReorgLaw::ScriptedV4 | ReorgLaw::ScriptedV5 => breakup_v2(ctx, &geom, state, plates, layer, sut, rs, g, tc, kind, t0)?,
        };
        if let Some(mut e) = ev {
            let gs2 = groups(&geom, state, plates, sut);
            e.groups_counted = counted(&gs2, p.plate_count_min_area_fraction);
            for root in [e.root_a, e.root_b] {
                if let Some(g2) = gs2.list.iter().find(|g| g.root.0 == root) {
                    if g2.cont_area / gs2.total_cont.max(1.0) >= p.breakup_piece_max_share {
                        rs.cascade.push((root, t0 + p.breakup_cascade_my));
                    }
                }
            }
            rs.events.push(e);
            rep.breakups += 1;
        }
    }
    let gs = groups(&geom, state, plates, sut);
    rep.groups_counted = counted(&gs, p.plate_count_min_area_fraction);
    rep.largest_group_share = gs.list.iter().map(|g| g.cont_area).fold(0.0, f64::max) / gs.total_cont.max(1.0);
    Ok(rep)
}

impl Groups {
    /// True when plate `x` is in group `gi` (used to exclude the group's own
    /// plates when looking for a neighbour).
    fn of_cell_plate(&self, x: &PlateId, gi: usize, sut: &SutureState) -> bool {
        sut.find(*x) == self.list[gi].root
    }
}

/// Plume-head split under `plate-reorg-scripted-v1`: the head's group
/// (holding at least the minimum continental share) splits along the
/// weak-zone min cut through the plume (strength x ratio within the head
/// radius); side B (the +n side of the plume's cut normal) opens about
/// p x n at the rift opening rate. Returns (event, cut normal) or None.
#[allow(clippy::too_many_arguments)]
pub fn plume_split(
    ctx: &ReorgCtx<'_>,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    plume: &Plume,
    t_my: f64,
) -> Result<Option<(ReorgEvent, [f64; 3], PlateId, PlateId, EulerPoleV0)>, String> {
    plume_split_law(ctx, ReorgLaw::ScriptedV1, state, plates, layer, sut, rs, plume, t_my)
}

/// `plume_split` for a chosen law version (lane BC): v2 cuts with
/// `domain_cut_v2` and rejects peels.
#[allow(clippy::too_many_arguments)]
pub fn plume_split_law(
    ctx: &ReorgCtx<'_>,
    law: ReorgLaw,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    plume: &Plume,
    t_my: f64,
) -> Result<Option<(ReorgEvent, [f64; 3], PlateId, PlateId, EulerPoleV0)>, String> {
    if law != ReorgLaw::ScriptedV1 {
        return plume_split_v2(ctx, state, plates, layer, sut, rs, plume, t_my);
    }
    let p = ctx.params;
    let geom = Geom::new(ctx.grid)?;
    let c0 = nearest_cell(ctx.grid, plume.position).map_err(|e| e.to_string())?;
    if state.kind[c0] != CrustKindCodeV0::Continental {
        return Ok(None);
    }
    let gs = groups(&geom, state, plates, sut);
    let gi = gs.of_cell[c0];
    let g = &gs.list[gi];
    let share = g.cont_area / gs.total_cont.max(1.0);
    if share < p.plume_split_min_share {
        return Ok(None);
    }
    let mut rng = SplitMix64::new(seed_from_path_tagged(ctx.seed_path, &format!("plate-reorg-scripted-v1/plume/{}", plume.id)));
    let noise = Noise::new(&mut rng, geom.radius);
    let ang = p.plume_head_radius_m / geom.radius;
    let strength = strengths(&geom, state, &g.cells, p.weak_zone_strength_ratio, &noise, Some((plume.position, ang)));
    let n = crate::rift::cut_normal(plume.position, plume.split_azimuth_rad);
    let cont: Vec<usize> = g.cells.iter().cloned().filter(|&c| state.kind[c] == CrustKindCodeV0::Continental).collect();
    let mut sv: Vec<(f64, f64)> = cont.iter().map(|&c| (dot(geom.cents[c], n), geom.areas[c])).collect();
    let d = ang.sin();
    let lo = weighted_quantile(&mut sv, BAND_LO).min(-d);
    let hi = weighted_quantile(&mut sv, BAND_HI).max(d);
    let src: Vec<bool> = g.cells.iter().map(|&c| dot(geom.cents[c], n) <= lo).collect();
    let snk: Vec<bool> = g.cells.iter().map(|&c| dot(geom.cents[c], n) >= hi).collect();
    if !src.iter().any(|&x| x) || !snk.iter().any(|&x| x) {
        return Ok(None);
    }
    let in_group = local_index(state.cells(), &g.cells);
    let cut = min_cut(&geom, state, &g.cells, &in_group, &strength, &src, &snk, p.weak_zone_strength_ratio);
    let sb = cont_share_b(&geom, state, g, &cut.side_b);
    if sb <= 0.0 || sb >= 1.0 {
        return Ok(None);
    }
    // Opening about p x n at the plume (the split_plate rule).
    let Some(ax) = unit(cross(plume.position, n)) else { return Ok(None) };
    let omega_rel = scale(ax, p.rift_opening_rate_deg_per_my.to_radians());
    let w_old = omega_of(plates, g.root)?;
    let out = apply_split(&geom, ctx.grid, state, plates, layer, sut, g, &cut.side_b, omega_rel, t_my)?;
    rs.relations.push(Relation { parent: out.root_a.0, child: out.root_b.0, omega_rel, t_my, kind: 3 });
    let e = ReorgEvent {
        t_my,
        kind: "plume-split".into(),
        root_a: out.root_a.0,
        root_b: out.root_b.0,
        cells_moved: out.cells_moved,
        group_share: share,
        side_b_share: sb,
        cut_length_m: cut.cut_length_m,
        weak_length_m: cut.weak_length_m,
        groups_counted: counted(&groups(&geom, state, plates, sut), p.plate_count_min_area_fraction),
    };
    rs.events.push(e.clone());
    let child_pole = pole_of(add(w_old, omega_rel))?;
    Ok(Some((e, n, out.root_a, out.root_b, child_pole)))
}

// ---------------------------------------------------------------------------
// SEED-BREAKUP-CUT-V2 (lane BC, 2026-10-13): `plate-reorg-scripted-v2` cut.
// Pre-registered in `unpublished-design-record`. Fixes the
// v1 cut that traced coastlines (design review review 2026-10-13 section 1): v1 made
// thin (rifted) margins weak and cut through the group's ocean at strength 1,
// so peeling a continent off its ocean along the margin was cheaper than
// crossing it.
//
// - Strength: cumulative shortening (`shortening_added_m`), old sutures
//   (edges between member plates, continental on both sides) and the plume
//   head disc only; no thickness term.
// - Domain: the group's continental cells plus a one-cell halo whose strength
//   is at or above an unweakened craton at the same place, so the cut must
//   cross the continent; the rest of the group takes the side of the nearest
//   domain cell, which extends the boundary through the ocean from the two
//   coast exit points.
// - Peels (side-B continental share outside [0.15, 0.85]) are rejected.
// - Integer capacities: the max flow is exact and the cut is the canonical
//   minimal source side.

/// Capacity units per metre of cut weight x strength (mm resolution).
const CAP_PER_M_V2: f64 = 1_000.0;
/// Source / sink arc capacity of the integer flow.
const INF_V2: u64 = u64::MAX / 8;

/// Dinic max flow on integer capacities (exact, platform independent given
/// the capacities).
struct FlowI {
    head: Vec<usize>,
    to: Vec<usize>,
    cap: Vec<u64>,
    next: Vec<usize>,
}

impl FlowI {
    fn new(n: usize) -> Self {
        Self { head: vec![NIL; n], to: Vec::new(), cap: Vec::new(), next: Vec::new() }
    }
    fn arc(&mut self, u: usize, v: usize, c_uv: u64, c_vu: u64) {
        self.to.push(v);
        self.cap.push(c_uv);
        self.next.push(self.head[u]);
        self.head[u] = self.to.len() - 1;
        self.to.push(u);
        self.cap.push(c_vu);
        self.next.push(self.head[v]);
        self.head[v] = self.to.len() - 1;
    }
    fn bfs(&self, s: usize, level: &mut [i32]) {
        level.iter_mut().for_each(|l| *l = -1);
        let mut q = std::collections::VecDeque::new();
        level[s] = 0;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            let mut e = self.head[u];
            while e != NIL {
                let v = self.to[e];
                if level[v] < 0 && self.cap[e] > 0 {
                    level[v] = level[u] + 1;
                    q.push_back(v);
                }
                e = self.next[e];
            }
        }
    }
    /// Iterative blocking-flow DFS (no recursion: domains can be long).
    fn augment(&mut self, s: usize, t: usize, level: &[i32], it: &mut [usize]) -> u64 {
        let mut path: Vec<usize> = Vec::new(); // arcs
        let mut u = s;
        loop {
            if u == t {
                let f = path.iter().map(|&e| self.cap[e]).min().unwrap_or(0);
                for &e in &path {
                    self.cap[e] -= f;
                    self.cap[e ^ 1] += f;
                }
                return f;
            }
            let mut advanced = false;
            while it[u] != NIL {
                let e = it[u];
                let v = self.to[e];
                if self.cap[e] > 0 && level[v] == level[u] + 1 {
                    path.push(e);
                    u = v;
                    advanced = true;
                    break;
                }
                it[u] = self.next[e];
            }
            if !advanced {
                // Dead end: retreat and skip the arc that led here.
                match path.pop() {
                    None => return 0,
                    Some(e) => {
                        u = self.to[e ^ 1];
                        it[u] = self.next[e];
                    }
                }
            }
        }
    }
    fn max_flow(&mut self, s: usize, t: usize) -> u64 {
        let n = self.head.len();
        let mut level = vec![-1i32; n];
        let mut total: u64 = 0;
        loop {
            self.bfs(s, &mut level);
            if level[t] < 0 {
                break;
            }
            let mut it = self.head.clone();
            loop {
                let f = self.augment(s, t, &level, &mut it);
                if f == 0 {
                    break;
                }
                total = total.saturating_add(f);
            }
        }
        total
    }
    fn reachable(&self, s: usize) -> Vec<bool> {
        let mut seen = vec![false; self.head.len()];
        let mut stack = vec![s];
        seen[s] = true;
        while let Some(u) = stack.pop() {
            let mut e = self.head[u];
            while e != NIL {
                let v = self.to[e];
                if !seen[v] && self.cap[e] > 0 {
                    seen[v] = true;
                    stack.push(v);
                }
                e = self.next[e];
            }
        }
        seen
    }
}

/// The v2 cut domain of one group: continental cells plus the one-cell
/// halo, in ascending cell order, with their strengths.
struct DomainV2 {
    cells: Vec<usize>,
    /// Global cell -> domain index (NIL outside).
    local: Vec<usize>,
    cont: Vec<bool>,
    strength: Vec<f64>,
}

fn domain_v2(geom: &Geom, state: &CrustState, g: &Group, ratio: f64, noise: &Noise, weak_disc: Option<([f64; 3], f64)>) -> DomainV2 {
    let n = state.cells();
    let in_group = local_index(n, &g.cells);
    let is_cont = |c: usize| state.kind[c] == CrustKindCodeV0::Continental;
    let mut nb = Vec::with_capacity(8);
    let mut cells = Vec::new();
    for &c in &g.cells {
        if is_cont(c) {
            cells.push(c);
            continue;
        }
        geom.nbrs8(c, &mut nb);
        if nb.iter().any(|&(k, _)| in_group[k] != NIL && is_cont(k)) {
            cells.push(c);
        }
    }
    let local = local_index(n, &cells);
    let cont: Vec<bool> = cells.iter().map(|&c| is_cont(c)).collect();
    let strength = cells
        .iter()
        .zip(&cont)
        .map(|(&c, &ct)| {
            let nf = (NOISE_SIGMA * noise.eval(geom.cents[c])).exp();
            if !ct {
                // Halo (margin): at or above an unweakened craton here.
                return nf.max(1.0);
            }
            let ws = (state.shortening_added_m[c] / SHORTENING_WEAK_SCALE_M).clamp(0.0, 1.0);
            let mut s = (1.0 - (1.0 - ratio) * ws) * nf;
            if let Some((p, ang)) = weak_disc {
                if dot(p, geom.cents[c]).clamp(-1.0, 1.0).acos() <= ang {
                    s *= ratio;
                }
            }
            s
        })
        .collect();
    DomainV2 { cells, local, cont, strength }
}

/// A v2 cut: side B over the group's cells (after the ocean extension), the
/// cut length and weak length through the domain, and the domain cut centre.
struct CutV2 {
    cut: Cut,
    centre: Option<[f64; 3]>,
}

/// Mean strength of a domain link times the seam factor (old suture: both
/// cells continental on different member plates).
fn link_strength_v2(state: &CrustState, d: &DomainV2, a: usize, b: usize, ratio: f64) -> f64 {
    let mut sm = 0.5 * (d.strength[a] + d.strength[b]);
    if d.cont[a] && d.cont[b] && state.plate[d.cells[a]] != state.plate[d.cells[b]] {
        sm *= ratio;
    }
    sm
}

/// Min cut over the domain between the domain cells flagged `src` / `snk`
/// (domain order), then the ocean extension over the rest of group `g`.
fn domain_cut_v2(geom: &Geom, state: &CrustState, g: &Group, d: &DomainV2, src: &[bool], snk: &[bool], ratio: f64) -> Option<CutV2> {
    let m = d.cells.len();
    if m == 0 || !src.iter().any(|&x| x) || !snk.iter().any(|&x| x) {
        return None;
    }
    let (s, t) = (m, m + 1);
    let mut fl = FlowI::new(m + 2);
    let mut nb = Vec::with_capacity(8);
    for (a, &c) in d.cells.iter().enumerate() {
        geom.nbrs8(c, &mut nb);
        for &(k, w) in &nb {
            let b = d.local[k];
            if b == NIL || b <= a {
                continue;
            }
            let cap = ((w * link_strength_v2(state, d, a, b, ratio) * CAP_PER_M_V2).round() as u64).max(1);
            fl.arc(a, b, cap, cap);
        }
        if src[a] {
            fl.arc(s, a, INF_V2, 0);
        }
        if snk[a] {
            fl.arc(a, t, INF_V2, 0);
        }
    }
    fl.max_flow(s, t);
    let reach = fl.reachable(s);
    // B: not reachable, connected to the sink through non-reachable cells
    // (the v1 rule); everything else in the domain is A.
    let mut dom_b = vec![false; m];
    let mut stack: Vec<usize> = (0..m).filter(|&a| snk[a] && !reach[a]).collect();
    for &a in &stack {
        dom_b[a] = true;
    }
    while let Some(a) = stack.pop() {
        geom.nbrs8(d.cells[a], &mut nb);
        for &(k, _) in &nb {
            let b = d.local[k];
            if b != NIL && !reach[b] && !dom_b[b] {
                dom_b[b] = true;
                stack.push(b);
            }
        }
    }
    // Cut length / weak length / centre over 4-neighbour domain edges.
    let (mut len, mut weak) = (0.0, 0.0);
    let mut cs = [0.0; 3];
    let mut nb4 = Vec::with_capacity(4);
    for (a, &c) in d.cells.iter().enumerate() {
        if dom_b[a] {
            continue;
        }
        geom.nbrs4(c, &mut nb4);
        for &(k, e) in &nb4 {
            let b = d.local[k];
            if b != NIL && dom_b[b] {
                len += e;
                if link_strength_v2(state, d, a, b, ratio) < WEAK_LOG_THRESHOLD {
                    weak += e;
                }
                if let Some(mid) = unit(add(geom.cents[c], geom.cents[k])) {
                    cs = add(cs, scale(mid, e));
                }
            }
        }
    }
    // Ocean extension: the rest of the group takes the side of the nearest
    // domain cell through group cells (8-neighbour Dijkstra, great-circle
    // steps; ties settle in (distance, group-local index) order).
    let gm = g.cells.len();
    let in_group = local_index(state.cells(), &g.cells);
    let mut lab = vec![0u8; gm]; // 0 unknown, 1 A, 2 B
    let mut dist = vec![f64::INFINITY; gm];
    let mut heap = std::collections::BinaryHeap::new();
    for (a, &c) in d.cells.iter().enumerate() {
        let ga = in_group[c];
        dist[ga] = 0.0;
        lab[ga] = if dom_b[a] { 2 } else { 1 };
        heap.push(std::cmp::Reverse((0u64, ga)));
    }
    let mut done = vec![false; gm];
    while let Some(std::cmp::Reverse((dk, ga))) = heap.pop() {
        if done[ga] || f64::from_bits(dk) > dist[ga] {
            continue;
        }
        done[ga] = true;
        let dcur = dist[ga];
        geom.nbrs8(g.cells[ga], &mut nb);
        for &(k, _) in &nb {
            let gb = in_group[k];
            if gb == NIL || done[gb] {
                continue;
            }
            let step = dot(geom.cents[g.cells[ga]], geom.cents[k]).clamp(-1.0, 1.0).acos() * geom.radius;
            let nd = dcur + step;
            if nd < dist[gb] {
                dist[gb] = nd;
                lab[gb] = lab[ga];
                heap.push(std::cmp::Reverse((nd.to_bits(), gb)));
            }
        }
    }
    let side_b: Vec<bool> = lab.iter().map(|&l| l == 2).collect();
    Some(CutV2 { cut: Cut { side_b, cut_length_m: len, weak_length_m: weak }, centre: unit(cs) })
}

/// v2 breakup of group `g`: seeded axis, domain cut, peel rejection with up
/// to MAX_AXIS_ATTEMPTS axes (same streams as v1).
#[allow(clippy::too_many_arguments)]
fn breakup_v2(
    ctx: &ReorgCtx<'_>,
    geom: &Geom,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    g: &Group,
    total_cont: f64,
    kind: &str,
    t_my: f64,
) -> Result<Option<ReorgEvent>, String> {
    let p = ctx.params;
    let k = rs.breakups;
    rs.breakups += 1;
    let mut rng = SplitMix64::new(seed_from_path_tagged(ctx.seed_path, &format!("plate-reorg-scripted-v1/breakup/{k}")));
    let noise = Noise::new(&mut rng, geom.radius);
    let cont: Vec<usize> = g.cells.iter().cloned().filter(|&c| state.kind[c] == CrustKindCodeV0::Continental).collect();
    let Some(gc) = centroid(geom, cont.iter().cloned()) else { return Ok(None) };
    let dom = domain_v2(geom, state, g, p.weak_zone_strength_ratio, &noise, None);
    let (e1, e2) = tangent_basis(gc);
    for _attempt in 0..MAX_AXIS_ATTEMPTS {
        let az = 2.0 * PI * rng.next_f64();
        let axis = add(scale(e1, az.cos()), scale(e2, az.sin()));
        let mut sv: Vec<(f64, f64)> = cont.iter().map(|&c| (dot(geom.cents[c], axis), geom.areas[c])).collect();
        let lo = weighted_quantile(&mut sv, BAND_LO);
        let hi = weighted_quantile(&mut sv, BAND_HI);
        let src: Vec<bool> = dom.cells.iter().map(|&c| dot(geom.cents[c], axis) <= lo).collect();
        let snk: Vec<bool> = dom.cells.iter().map(|&c| dot(geom.cents[c], axis) >= hi).collect();
        let Some(cv) = domain_cut_v2(geom, state, g, &dom, &src, &snk, p.weak_zone_strength_ratio) else { continue };
        let sb = cont_share_b(geom, state, g, &cv.cut.side_b);
        if !(SIDE_MIN_CONT_SHARE..=1.0 - SIDE_MIN_CONT_SHARE).contains(&sb) {
            continue;
        }
        let Some(centre) = cv.centre else { continue };
        let Some((out, w)) =
            cut_and_split(ctx, geom, state, plates, layer, sut, g, &cv.cut, 1.0, p.rift_opening_rate_deg_per_my, Some(centre), t_my)?
        else {
            continue;
        };
        rs.relations.push(Relation { parent: out.root_a.0, child: out.root_b.0, omega_rel: w, t_my, kind: if kind == "cascade" { 1 } else { 0 } });
        return Ok(Some(ReorgEvent {
            t_my,
            kind: kind.into(),
            root_a: out.root_a.0,
            root_b: out.root_b.0,
            cells_moved: out.cells_moved,
            group_share: g.cont_area / total_cont.max(1.0),
            side_b_share: sb,
            cut_length_m: cv.cut.cut_length_m,
            weak_length_m: cv.cut.weak_length_m,
            groups_counted: 0,
        }));
    }
    Ok(None)
}

/// v2 plume-head split: the v1 rule (floor, bands through the plume, opening
/// about p x n) on the v2 domain cut; a peel does not split.
#[allow(clippy::too_many_arguments)]
fn plume_split_v2(
    ctx: &ReorgCtx<'_>,
    state: &mut CrustState,
    plates: &mut Vec<PlateRecordV0>,
    layer: &mut ParcelLayer,
    sut: &mut SutureState,
    rs: &mut ReorgState,
    plume: &Plume,
    t_my: f64,
) -> Result<Option<(ReorgEvent, [f64; 3], PlateId, PlateId, EulerPoleV0)>, String> {
    let p = ctx.params;
    let geom = Geom::new(ctx.grid)?;
    let c0 = nearest_cell(ctx.grid, plume.position).map_err(|e| e.to_string())?;
    if state.kind[c0] != CrustKindCodeV0::Continental {
        return Ok(None);
    }
    let gs = groups(&geom, state, plates, sut);
    let g = &gs.list[gs.of_cell[c0]];
    let share = g.cont_area / gs.total_cont.max(1.0);
    if share < p.plume_split_min_share {
        return Ok(None);
    }
    let mut rng = SplitMix64::new(seed_from_path_tagged(ctx.seed_path, &format!("plate-reorg-scripted-v1/plume/{}", plume.id)));
    let noise = Noise::new(&mut rng, geom.radius);
    let ang = p.plume_head_radius_m / geom.radius;
    let dom = domain_v2(&geom, state, g, p.weak_zone_strength_ratio, &noise, Some((plume.position, ang)));
    let n = crate::rift::cut_normal(plume.position, plume.split_azimuth_rad);
    let cont: Vec<usize> = g.cells.iter().cloned().filter(|&c| state.kind[c] == CrustKindCodeV0::Continental).collect();
    let mut sv: Vec<(f64, f64)> = cont.iter().map(|&c| (dot(geom.cents[c], n), geom.areas[c])).collect();
    let d = ang.sin();
    let lo = weighted_quantile(&mut sv, BAND_LO).min(-d);
    let hi = weighted_quantile(&mut sv, BAND_HI).max(d);
    let src: Vec<bool> = dom.cells.iter().map(|&c| dot(geom.cents[c], n) <= lo).collect();
    let snk: Vec<bool> = dom.cells.iter().map(|&c| dot(geom.cents[c], n) >= hi).collect();
    let Some(cv) = domain_cut_v2(&geom, state, g, &dom, &src, &snk, p.weak_zone_strength_ratio) else { return Ok(None) };
    let sb = cont_share_b(&geom, state, g, &cv.cut.side_b);
    if !(SIDE_MIN_CONT_SHARE..=1.0 - SIDE_MIN_CONT_SHARE).contains(&sb) {
        return Ok(None);
    }
    let Some(ax) = unit(cross(plume.position, n)) else { return Ok(None) };
    let omega_rel = scale(ax, p.rift_opening_rate_deg_per_my.to_radians());
    let w_old = omega_of(plates, g.root)?;
    let out = apply_split(&geom, ctx.grid, state, plates, layer, sut, g, &cv.cut.side_b, omega_rel, t_my)?;
    rs.relations.push(Relation { parent: out.root_a.0, child: out.root_b.0, omega_rel, t_my, kind: 3 });
    let e = ReorgEvent {
        t_my,
        kind: "plume-split".into(),
        root_a: out.root_a.0,
        root_b: out.root_b.0,
        cells_moved: out.cells_moved,
        group_share: share,
        side_b_share: sb,
        cut_length_m: cv.cut.cut_length_m,
        weak_length_m: cv.cut.weak_length_m,
        groups_counted: counted(&groups(&geom, state, plates, sut), p.plate_count_min_area_fraction),
    };
    rs.events.push(e.clone());
    let child_pole = pole_of(add(w_old, omega_rel))?;
    Ok(Some((e, n, out.root_a, out.root_b, child_pole)))
}

/// Plumes whose tails are still active in a step starting at `t0`.
pub fn active_plumes(plumes: &[Plume], lifetime_my: f64, t0: f64) -> Vec<Plume> {
    plumes.iter().filter(|p| p.onset_my + lifetime_my > t0).cloned().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// PLATE-POLES `rigid-group-contacts-v1` part 2 on a 36 x 72 grid:
    /// plate 0 everywhere else; plates 1 and 2 touch; plate 3 is far away.
    fn release_world(parent: &[(u32, u32)]) -> (Geom, CrustState, Vec<PlateRecordV0>, SutureState) {
        let grid = GridSpec::new(36, 72, 6.371e6).unwrap();
        let geom = Geom::new(&grid).unwrap();
        let n = grid.cells();
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
            boundary_class: vec![0; n],
            normal_rate_m_per_my: vec![0.0; n],
            material_pos: cents,
        };
        for i in 14..22 {
            for j in 10..20 {
                st.plate[i * 72 + j] = PlateId(1);
            }
            for j in 20..26 {
                st.plate[i * 72 + j] = PlateId(2);
            }
            for j in 50..53 {
                st.plate[i * 72 + j] = PlateId(3);
            }
        }
        let pole = EulerPoleV0 { lat_deg: 10.0, lon_deg: 20.0, rate_deg_per_my: 0.5 };
        let plates = (0..4)
            .map(|id| PlateRecordV0 { plate_id: PlateId(id), parent_plate: None, birth_t_my: 0.0, euler_pole: Some(if id == 0 { EulerPoleV0 { lat_deg: 90.0, lon_deg: 0.0, rate_deg_per_my: 0.0 } } else { pole }), retirement: None })
            .collect();
        let sut = SutureState { parent: parent.iter().map(|&(a, b)| (PlateId(a), PlateId(b))).collect() };
        (geom, st, plates, sut)
    }

    #[test]
    fn detached_member_leaves_its_group() {
        let (geom, st, plates, mut sut) = release_world(&[(2, 1), (3, 1)]);
        let mut rs = ReorgState::default();
        rs.tenure.push((1, 5.0));
        let n = release_detached(&geom, &st, &plates, &mut sut, &mut rs, 10.0);
        assert_eq!(n, 1);
        assert_eq!(sut.find(PlateId(2)), PlateId(1));
        assert_eq!(sut.find(PlateId(3)), PlateId(3));
        assert_eq!(rs.tenure, vec![(1, 5.0)]);
        assert_eq!(rs.events.last().map(|e| (e.kind.as_str(), e.root_a, e.root_b)), Some(("release", 1, 3)));
        // Idempotent: nothing left to release.
        assert_eq!(release_detached(&geom, &st, &plates, &mut sut, &mut rs, 10.0), 0);
    }

    #[test]
    fn detached_root_hands_the_group_to_the_main_body() {
        // Root 3 is the small far piece: the main body {1, 2} is re-rooted on
        // its largest member (1) and the timers and relations follow it.
        let (geom, st, plates, mut sut) = release_world(&[(1, 3), (2, 3)]);
        let mut rs = ReorgState::default();
        rs.tenure.push((3, 5.0));
        rs.relations.push(Relation { parent: 0, child: 3, omega_rel: [0.0; 3], t_my: 1.0, kind: 0 });
        assert_eq!(release_detached(&geom, &st, &plates, &mut sut, &mut rs, 10.0), 1);
        assert_eq!(sut.find(PlateId(2)), PlateId(1));
        assert_eq!(sut.find(PlateId(1)), PlateId(1));
        assert_eq!(sut.find(PlateId(3)), PlateId(3));
        assert_eq!(rs.tenure, vec![(1, 5.0)]);
        assert_eq!(rs.relations[0].child, 1);
    }

    #[test]
    fn flow_finds_the_weak_column() {
        // A 1 x 5 chain: s - 0 - 1 - 2 - 3 - 4 - t with a weak link 2-3.
        let mut f = Flow::new(7);
        let caps = [5.0, 4.0, 0.5, 3.0];
        for (i, c) in caps.iter().enumerate() {
            f.arc(i, i + 1, *c, *c);
        }
        f.arc(5, 0, 1e30, 0.0);
        f.arc(4, 6, 1e30, 0.0);
        let v = f.max_flow(5, 6);
        assert!((v - 0.5).abs() < 1e-9, "{v}");
        let r = f.reachable(5);
        assert!(r[0] && r[1] && r[2] && !r[3] && !r[4]);
    }

    #[test]
    fn integer_flow_takes_the_minimal_source_side_on_ties() {
        // Lane BC (v2): chain s - 0 - 1 - 2 - 3 - 4 - t with two equal weakest
        // links (1-2 and 2-3, capacity 1): the canonical cut is the one
        // nearest the source, whatever the arc order.
        for rev in [false, true] {
            let mut f = FlowI::new(7);
            let caps = [3u64, 1, 1, 3];
            let order: Vec<usize> = if rev { (0..4).rev().collect() } else { (0..4).collect() };
            for i in order {
                f.arc(i, i + 1, caps[i], caps[i]);
            }
            f.arc(5, 0, INF_V2, 0);
            f.arc(4, 6, INF_V2, 0);
            assert_eq!(f.max_flow(5, 6), 1);
            let r = f.reachable(5);
            assert!(r[0] && r[1] && !r[2] && !r[3] && !r[4], "rev {rev}: {r:?}");
        }
    }

    #[test]
    fn v4_limits_validate() {
        // SEED-TRANSFORMS (lane TF): pre-registered values, typed refusals.
        let l = InitiationLimits::v1();
        assert_eq!((l.slab_age_min_my, l.old_slab_fraction_min, l.trench_spacing_m), (100.0, 0.5, 1_000_000.0));
        assert!(l.validate().is_ok());
        assert!(InitiationLimits { old_slab_fraction_min: 1.5, ..l }.validate().is_err());
        assert!(InitiationLimits { slab_age_min_my: f64::NAN, ..l }.validate().is_err());
        assert!(InitiationLimits { trench_spacing_m: -1.0, ..l }.validate().is_err());
    }

    #[test]
    fn v5_margin_initiation_validates() {
        // SEED-OLD-OCEAN (lane OR): pre-registered values, typed refusals.
        let m = MarginInitiation::v1();
        assert_eq!((m.age_my, m.segment_length_m, m.spacing_m), (160.0, 2_000_000.0, 1_000_000.0));
        assert!(m.validate().is_ok());
        assert!(MarginInitiation { age_my: 0.0, ..m }.validate().is_err());
        assert!(MarginInitiation { segment_length_m: f64::INFINITY, ..m }.validate().is_err());
        assert!(MarginInitiation { spacing_m: -1.0, ..m }.validate().is_err());
    }

    #[test]
    fn schedule_and_plumes_are_prefix_stable() {
        let sp = SeedPath { seed: 42, derivation: "x".into(), salt: "y".into() };
        let a = reorg_schedule(&sp, 50.0, 100.0, 400.0);
        let b = reorg_schedule(&sp, 50.0, 100.0, 600.0);
        assert!(!a.is_empty() && a.len() < b.len());
        assert_eq!(a[..], b[..a.len()]);
        assert!(a.windows(2).all(|w| w[1] - w[0] >= 50.0 && w[1] - w[0] <= 100.0));
        let pa = seed_plumes_renewal(25.0, 50.0, 400.0, &sp);
        let pb = seed_plumes_renewal(25.0, 50.0, 600.0, &sp);
        assert_eq!(pa[..], pb[..pa.len()]);
        assert!(pa.iter().all(|p| p.onset_my <= 400.0));
    }
}
