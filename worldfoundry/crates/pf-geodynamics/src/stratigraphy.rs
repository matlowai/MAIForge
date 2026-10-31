//! SEED-STRATIGRAPHY stage A (lane ST, 2026-10-13): `parcel-stacks-v1`.
//! Spec: `unpublished-design-record` sections 1, 2
//! and the stage A paragraph of 5; constants and predictions in
//! `unpublished-design-record`.
//!
//! Every continental parcel (`conserve::ParcelLayer`) carries a stack:
//! - a basement (lithology crystalline, or oceanic igneous for parcels born
//!   from converted island arcs) whose thickness is the RESIDUAL of the
//!   parcel column (thickness + sediment) minus its cover, so the stack and
//!   the parcel ledger agree by construction;
//! - up to `max_layers - 1` cover layers, bottom-up, each with lithology,
//!   environment, thickness, age span, metamorphic grade, deepest burial and
//!   flags.
//!
//! The step hands the parcels TAGGED per-cell deltas (`CellTags`) with the
//! same deposit weights `conserve::write_back` uses; per parcel they are
//! applied in a declared order: removals (erosion truncates from the top,
//! thinning scales every layer, negative shortening removes basement), then
//! additions bottom-up (pluton at the basement top, basement additions,
//! positive shortening thickens every layer, volcanics on top, sediment on
//! top), then the unexplained remainder into the basement (counted as
//! `basement_flow`: channel flow, plateau spill, delamination and kernel
//! weight differences). Cover buried below the metamorphic depth becomes
//! metamorphic. Merges keep the stack within its capacity.
//!
//! A PASSIVE record: nothing here writes a field the simulation reads.
//! Determinism: fixed parcel order, f64 throughout, ascending reductions, no
//! hash maps; quantization only at publish.

use crate::kinematic::KinematicError;
use pf_contracts::terrain::{
    AssumptionV0, StratigraphyCodesV0, StratigraphyDomainV0, STRATA_FLAG_BASEMENT,
    STRATA_FLAG_CUT_AT_DEPTH, STRATA_FLAG_MERGED, STRATA_FLAG_ORGANIC_RICH,
    STRATA_FLAG_TRUNCATED_TOP, STRATA_FLAG_UNCONFORMITY_BASE, STRATIGRAPHY_PREVIEW_USAGE_LIMITS,
    STRATIGRAPHY_PREVIEW_USAGE_LIMITS_V2, STRATIGRAPHY_SCHEMA_V0,
};
use pf_grid::GridSpec;

pub const LITH_BASEMENT: u8 = 0;
pub const LITH_GRANITE: u8 = 1;
pub const LITH_METAMORPHIC: u8 = 2;
pub const LITH_BASALT: u8 = 3;
pub const LITH_CARBONATE: u8 = 4;
pub const LITH_SHALE: u8 = 5;
pub const LITH_SANDSTONE: u8 = 6;
pub const LITH_CONGLOMERATE: u8 = 7;
pub const LITH_EVAPORITE: u8 = 8;
pub const LITH_COAL: u8 = 9;
pub const LITH_UNCONSOLIDATED: u8 = 10;
pub const LITH_OCEANIC: u8 = 11;
pub const NLITH: usize = 12;

pub const ENV_BASEMENT: u8 = 0;
pub const ENV_PLUTON: u8 = 1;
pub const ENV_VOLCANIC: u8 = 2;
pub const ENV_FLUVIAL: u8 = 3;
pub const ENV_LACUSTRINE: u8 = 4;
pub const ENV_DELTAIC: u8 = 5;
pub const ENV_SHELF: u8 = 6;
pub const ENV_DEEP_MARINE: u8 = 7;
pub const ENV_EVAPORITIC: u8 = 8;
pub const ENV_SWAMP: u8 = 9;
pub const ENV_OCEANIC: u8 = 10;

/// Deposit classes (lithology, environment), indexed by `CellTags::class`.
pub const DEPOSIT_CLASSES: [(u8, u8); 9] = [
    (LITH_SANDSTONE, ENV_FLUVIAL),
    (LITH_CONGLOMERATE, ENV_FLUVIAL),
    (LITH_COAL, ENV_SWAMP),
    (LITH_SHALE, ENV_LACUSTRINE),
    (LITH_EVAPORITE, ENV_EVAPORITIC),
    (LITH_SANDSTONE, ENV_DELTAIC),
    (LITH_CARBONATE, ENV_SHELF),
    (LITH_SHALE, ENV_SHELF),
    (LITH_SHALE, ENV_DEEP_MARINE),
];
pub const NCLASS: usize = DEPOSIT_CLASSES.len();
pub const CLASS_NAMES: [&str; NCLASS] = [
    "fluvial_sandstone",
    "fluvial_conglomerate",
    "swamp_coal",
    "lacustrine_shale",
    "closed_basin_evaporite",
    "deltaic_sandstone",
    "shelf_carbonate",
    "shelf_shale",
    "deep_marine_shale",
];
pub const LITH_NAMES: [&str; NLITH] = [
    "crystalline_basement",
    "granite_pluton",
    "metamorphic",
    "basalt_volcanic",
    "marine_carbonate",
    "shale",
    "sandstone",
    "conglomerate",
    "evaporite",
    "coal_bearing",
    "unconsolidated",
    "oceanic_igneous",
];

/// Constants of `parcel-stacks-v1` (PREREG.txt).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrataParams {
    pub max_layers: usize,
    pub merge_gap_my: f64,
    pub unconformity_gap_my: f64,
    pub min_layer_m: f64,
    pub trace_m: f64,
    pub arc_extrusive_fraction: f64,
    pub plume_extrusive_fraction: f64,
    pub shortening_pluton_fraction: f64,
    /// PREREG revision 1: crystalline roof kept above a new pluton.
    pub pluton_roof_m: f64,
    pub metamorphic_depth_m: f64,
    pub high_grade_depth_m: f64,
    pub shelf_depth_m: f64,
    pub carbonate_max_paleolat_deg: f64,
    pub carbonate_max_supply_m_per_my: f64,
    pub arid_lat_min_deg: f64,
    pub arid_lat_max_deg: f64,
    pub humid_lat_max_deg: f64,
    pub humid_midlat_min_deg: f64,
    pub humid_midlat_max_deg: f64,
    pub coal_max_elevation_m: f64,
    pub conglomerate_relief_m: f64,
    pub river_mouth_area_m2: f64,
    pub publish_depth_m: f64,
    pub basement_min_m: f64,
    /// `parcel-stacks-v2` (lane SA, 2026-10-13); None = v1.
    pub v2: Option<StrataV2Params>,
}

/// The constants `parcel-stacks-v2` adds (lane SA PREREG.txt,
/// `unpublished-design-record`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrataV2Params {
    pub shelf_width_m: f64,
    pub slope_foot_depth_m: f64,
    pub shelf_majority_fraction: f64,
    pub platform_cover_max_m: f64,
    pub platform_full_elev_m: f64,
    pub platform_zero_elev_m: f64,
    pub platform_sandstone_fraction: f64,
    pub platform_carbonate_fraction: f64,
    pub platform_age_base_my: f64,
    pub platform_age_top_my: f64,
}

impl StrataParams {
    pub fn validate(&self) -> Result<(), String> {
        if self.max_layers < 2 {
            return Err(format!("strata max layers {}", self.max_layers));
        }
        for (n, v) in [
            ("merge gap", self.merge_gap_my),
            ("min layer", self.min_layer_m),
            ("pluton roof", self.pluton_roof_m),
            ("trace", self.trace_m),
            ("metamorphic depth", self.metamorphic_depth_m),
            ("high grade depth", self.high_grade_depth_m),
            ("shelf depth", self.shelf_depth_m),
            ("publish depth", self.publish_depth_m),
        ] {
            if !(v.is_finite() && v >= 0.0) {
                return Err(format!("strata {n} {v}"));
            }
        }
        for (n, v) in [
            ("arc extrusive fraction", self.arc_extrusive_fraction),
            ("plume extrusive fraction", self.plume_extrusive_fraction),
            ("shortening pluton fraction", self.shortening_pluton_fraction),
        ] {
            if !(0.0..=1.0).contains(&v) {
                return Err(format!("strata {n} {v}"));
            }
        }
        if self.high_grade_depth_m < self.metamorphic_depth_m {
            return Err("strata high-grade depth above the metamorphic depth".into());
        }
        if let Some(v) = &self.v2 {
            if !(v.slope_foot_depth_m > 0.5 * self.shelf_depth_m) {
                return Err("strata slope foot not below half the shelf depth".into());
            }
            if !(v.platform_zero_elev_m > v.platform_full_elev_m) {
                return Err("strata platform zero-cover elevation not above the full-cover elevation".into());
            }
            if v.platform_sandstone_fraction + v.platform_carbonate_fraction > 1.0 {
                return Err("strata platform sandstone + carbonate fractions above 1".into());
            }
            if !(v.platform_age_top_my >= v.platform_age_base_my) {
                return Err("strata platform age top below its base".into());
            }
            if !(v.shelf_width_m >= 0.0 && v.platform_cover_max_m >= 0.0 && (0.0..=1.0).contains(&v.shelf_majority_fraction)) {
                return Err("strata v2 constant out of range".into());
            }
        }
        Ok(())
    }
}

/// One cover layer.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StrataLayer {
    pub thickness_m: f64,
    pub age_base_my: f64,
    pub age_top_my: f64,
    pub max_burial_m: f64,
    pub lith: u8,
    pub env: u8,
    pub grade: u8,
    pub flags: u8,
}

/// One parcel's stack: the basement residual and the cover, BOTTOM-UP
/// (`cover[0]` sits on the basement, the last entry is the surface).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StrataStack {
    pub basement_m: f64,
    pub basement_lith: u8,
    pub basement_age_my: f64,
    pub last_exposure_my: f64,
    pub cover: Vec<StrataLayer>,
    /// `parcel-stacks-v2`: cumulative rock uplift since t = 0 or birth (m)
    /// and the reference elevation of the previous step's end (m); both 0
    /// and unused under v1.
    pub rock_uplift_m: f64,
    pub elev_ref_m: f64,
}

impl StrataStack {
    pub fn basement_only(h_m: f64, lith: u8, age_my: f64) -> Self {
        Self {
            basement_m: h_m,
            basement_lith: lith,
            basement_age_my: age_my,
            last_exposure_my: -1.0,
            cover: Vec::new(),
            rock_uplift_m: 0.0,
            elev_ref_m: 0.0,
        }
    }
    pub fn cover_m(&self) -> f64 {
        self.cover.iter().map(|l| l.thickness_m).sum()
    }
    pub fn total_m(&self) -> f64 {
        self.basement_m + self.cover_m()
    }
    /// Adds area x thickness per lithology into `out`.
    pub fn lith_volumes(&self, area_m2: f64, out: &mut [f64]) {
        out[self.basement_lith as usize] += self.basement_m * area_m2;
        for l in &self.cover {
            out[l.lith as usize] += l.thickness_m * area_m2;
        }
    }
    /// Lithology of the surface layer.
    pub fn top_lith(&self) -> u8 {
        self.cover.last().map(|l| l.lith).unwrap_or(self.basement_lith)
    }
}

/// The step's per-cell tagged deltas (m; all >= 0) that the caller knows
/// and the projection does not (arc and shortening are read from the
/// projection's lineage fields inside `conserve::write_back_strata`).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CellTags {
    pub erosion_m: Vec<f64>,
    pub deposit_m: Vec<f64>,
    /// Deposit class (index into `DEPOSIT_CLASSES`) of the cell this step.
    pub class: Vec<u8>,
    pub plume_m: Vec<f64>,
    pub thinning_m: Vec<f64>,
    pub basement_add_m: Vec<f64>,
    /// Subaerial (exposed) this step.
    pub land: Vec<bool>,
    /// `parcel-stacks-v2` only (empty under v1): the belt's placed collision
    /// supply this step (m, >= 0), the cell took the shelf facies in warm
    /// water, the post-tectonic (pre-erosion) elevation and the post-erosion
    /// elevation (m).
    pub belt_placed_m: Vec<f64>,
    pub warm_shelf: Vec<bool>,
    pub tect_elev_m: Vec<f64>,
    pub post_elev_m: Vec<f64>,
}

impl CellTags {
    pub fn zeros(n: usize) -> Self {
        Self {
            erosion_m: vec![0.0; n],
            deposit_m: vec![0.0; n],
            class: vec![0; n],
            plume_m: vec![0.0; n],
            thinning_m: vec![0.0; n],
            basement_add_m: vec![0.0; n],
            land: vec![false; n],
            belt_placed_m: Vec::new(),
            warm_shelf: Vec::new(),
            tect_elev_m: Vec::new(),
            post_elev_m: Vec::new(),
        }
    }
}

/// Per-step stratigraphy ledger (volumes m^3).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StrataStepLog {
    pub parcels: u64,
    pub stacks: u64,
    /// Stack volume at write-back entry (after consumption) and after.
    pub before_m3: f64,
    pub after_m3: f64,
    /// Parcel column (thickness + sediment) volume after.
    pub column_m3: f64,
    pub consumed_m3: [f64; NLITH],
    pub births_m3: [f64; NLITH],
    pub eroded_m3: [f64; NLITH],
    pub deposited_m3: [f64; NCLASS],
    pub volcanic_m3: f64,
    pub pluton_m3: f64,
    pub plume_basement_m3: f64,
    pub basement_add_m3: f64,
    pub thinned_m3: f64,
    pub shortening_scale_m3: f64,
    pub shortening_removed_m3: f64,
    pub basement_flow_m3: f64,
    pub basement_flow_abs_m3: f64,
    pub underflow_m3: f64,
    pub metamorphosed_m3: f64,
    pub merges_gap: u64,
    pub merges_cap: u64,
    pub merges_min: u64,
    /// after - (before + applied operations + births), and relative to after.
    pub closure_residual_m3: f64,
    pub closure_rel: f64,
    /// |after - column| / column.
    pub column_mismatch_rel: f64,
    /// Largest per-parcel |stack - column| / max(column, 1 m) at entry.
    pub entry_mismatch_max_rel: f64,
    pub mean_cover_layers: f64,
    pub stacks_at_cap: u64,
    pub basement_below_min: u64,
    /// `parcel-stacks-v2` (0 under v1): the non-placed shortening-lineage
    /// change applied to the basement (signed), warm shelf-facies deposits,
    /// capacity merges across an unconformity, and the area-weighted mean
    /// |rock uplift| of the parcels after the step.
    pub belt_basement_m3: f64,
    pub warm_shelf_deposit_m3: f64,
    pub merges_unconformity: u64,
    pub mean_abs_rock_uplift_m: f64,
}

impl StrataStepLog {
    /// Sum of every applied operation's signed volume.
    pub fn ops_m3(&self) -> f64 {
        -self.eroded_m3.iter().sum::<f64>()
            + self.deposited_m3.iter().sum::<f64>()
            + self.volcanic_m3
            + self.pluton_m3
            + self.plume_basement_m3
            + self.basement_add_m3
            - self.thinned_m3
            + self.shortening_scale_m3
            - self.shortening_removed_m3
            + self.basement_flow_m3
            + self.belt_basement_m3
    }
}

/// One parcel's amounts this step (m).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ParcelAmounts {
    pub erosion: f64,
    pub deposit: f64,
    pub class: u8,
    pub arc: f64,
    pub plume: f64,
    pub thinning: f64,
    pub basement_add: f64,
    pub short_pos: f64,
    pub short_neg: f64,
    pub land: bool,
    /// `parcel-stacks-v2`: the non-placed shortening-lineage change for the
    /// basement (signed) and whether the deposit is warm shelf facies.
    pub belt_rest: f64,
    pub warm_shelf: bool,
}

fn organic(lith: u8, env: u8) -> bool {
    lith == LITH_COAL || (lith == LITH_SHALE && (env == ENV_DEEP_MARINE || env == ENV_LACUSTRINE))
}

fn group(lith: u8) -> u8 {
    match lith {
        LITH_BASEMENT | LITH_GRANITE | LITH_METAMORPHIC | LITH_OCEANIC => 0,
        LITH_BASALT => 1,
        LITH_CARBONATE => 2,
        LITH_SHALE | LITH_EVAPORITE | LITH_COAL => 3,
        _ => 4,
    }
}

/// Declared merge cost: 0 same lithology, 1 same mechanical group, 2 else.
pub fn merge_cost(a: u8, b: u8) -> u8 {
    if a == b {
        0
    } else if group(a) == group(b) {
        1
    } else {
        2
    }
}

fn absorb(into: &mut StrataLayer, from: &StrataLayer) {
    into.thickness_m += from.thickness_m;
    into.age_base_my = into.age_base_my.min(from.age_base_my);
    into.age_top_my = into.age_top_my.max(from.age_top_my);
    into.max_burial_m = into.max_burial_m.max(from.max_burial_m);
    into.flags |= STRATA_FLAG_MERGED | (from.flags & STRATA_FLAG_ORGANIC_RICH);
}

/// Append on top (merge with the top layer when same lithology and
/// environment within the merge gap; trace additions extend the top).
#[allow(clippy::too_many_arguments)]

#[allow(clippy::too_many_arguments)]
fn push_top_observed(st: &mut StrataStack, b: &mut f64, lith: u8, env: u8, th: f64, t0: f64, t1: f64, p: &StrataParams, log: &mut StrataStepLog,mut observer:Option<&mut crate::layer_evolution::Observer>,kind:u8) {
    if !(th > 0.0) {
        return;
    }
    if let Some(top) = st.cover.last_mut() {
        if top.lith == lith && top.env == env && t0 - top.age_top_my <= p.merge_gap_my + 1e-9 {
            top.thickness_m += th;
            top.age_top_my = top.age_top_my.max(t1);
            top.flags &= !STRATA_FLAG_TRUNCATED_TOP;
            log.merges_gap += 1;
            if let Some(o)=observer.as_deref_mut(){o.capture(kind,crate::layer_evolution::Mapping::Identity,st,*b);}
            return;
        }
        if th < p.trace_m {
            top.thickness_m += th;
            if let Some(o)=observer.as_deref_mut(){o.capture(kind,crate::layer_evolution::Mapping::Identity,st,*b);}
            return;
        }
    } else if th < p.trace_m {
        *b += th;
        if let Some(o)=observer.as_deref_mut(){o.capture(kind,crate::layer_evolution::Mapping::Identity,st,*b);}
        return;
    }
    st.cover.push(StrataLayer {
        thickness_m: th,
        age_base_my: t0,
        age_top_my: t1,
        max_burial_m: 0.0,
        lith,
        env,
        grade: 0,
        flags: if organic(lith, env) { STRATA_FLAG_ORGANIC_RICH } else { 0 },
    });
    if let Some(o)=observer.as_deref_mut(){o.capture(kind,crate::layer_evolution::Mapping::Insert(st.cover.len()-1,false),st,*b);}

}

/// Insert at the basement top (pluton).
#[allow(clippy::too_many_arguments)]

#[allow(clippy::too_many_arguments)]
fn insert_bottom_observed(st: &mut StrataStack, b: &mut f64, lith: u8, env: u8, th: f64, t0: f64, t1: f64, p: &StrataParams, log: &mut StrataStepLog,mut observer:Option<&mut crate::layer_evolution::Observer>,kind:u8) {
    if !(th > 0.0) {
        return;
    }
    if let Some(first) = st.cover.first_mut() {
        if first.lith == lith && first.env == env && t0 - first.age_top_my <= p.merge_gap_my + 1e-9 {
            first.thickness_m += th;
            first.age_top_my = first.age_top_my.max(t1);
            log.merges_gap += 1;
            if let Some(o)=observer.as_deref_mut(){o.capture(kind,crate::layer_evolution::Mapping::Identity,st,*b);}
            return;
        }
    }
    if th < p.trace_m {
        *b += th;
        if let Some(o)=observer.as_deref_mut(){o.capture(kind,crate::layer_evolution::Mapping::Identity,st,*b);}
        return;
    }
    st.cover.insert(
        0,
        StrataLayer {
            thickness_m: th,
            age_base_my: t0,
            age_top_my: t1,
            max_burial_m: 0.0,
            lith,
            env,
            grade: 0,
            flags: 0,
        },
    );
    if let Some(o)=observer.as_deref_mut(){o.capture(kind,crate::layer_evolution::Mapping::Insert(0,false),st,*b);}

}

/// Merge the thinnest cover layer (ties: lowest index) into the neighbour
/// of lowest merge cost (ties: the one below).

fn merge_thinnest_observed(st: &mut StrataStack,b:f64,observer:Option<&mut crate::layer_evolution::Observer>) {
    let n = st.cover.len();
    if n < 2 {
        return;
    }
    let mut i = 0;
    for k in 1..n {
        if st.cover[k].thickness_m < st.cover[i].thickness_m {
            i = k;
        }
    }
    let li = st.cover[i].lith;
    let j = match (i > 0, i + 1 < n) {
        (true, true) => {
            if merge_cost(li, st.cover[i + 1].lith) < merge_cost(li, st.cover[i - 1].lith) {
                i + 1
            } else {
                i - 1
            }
        }
        (true, false) => i - 1,
        _ => i + 1,
    };
    let from = st.cover[i];
    absorb(&mut st.cover[j], &from);
    st.cover.remove(i);
    if let Some(o)=observer{o.capture(17,crate::layer_evolution::Mapping::Merge(j,i),st,b);}
}

/// `parcel-stacks-v2` capacity merge: among adjacent cover pairs, the
/// lexicographic minimum of (spans an unconformity gap, merge cost, age top
/// of the upper layer (oldest first), combined thickness, lower index); the
/// thicker layer absorbs the thinner. Returns true when the merged pair
/// spanned an unconformity.

fn merge_pair_v2_observed(st: &mut StrataStack, unconformity_gap_my: f64,b:f64,observer:Option<&mut crate::layer_evolution::Observer>) -> bool {
    let n = st.cover.len();
    if n < 2 {
        return false;
    }
    let key = |i: usize| {
        let (lo, up) = (&st.cover[i], &st.cover[i + 1]);
        let unc = (up.age_base_my - lo.age_top_my > unconformity_gap_my) as u8;
        (unc, merge_cost(lo.lith, up.lith), up.age_top_my, lo.thickness_m + up.thickness_m)
    };
    let mut best = 0usize;
    let mut bk = key(0);
    for i in 1..n - 1 {
        let k = key(i);
        let better = (k.0, k.1) < (bk.0, bk.1)
            || ((k.0, k.1) == (bk.0, bk.1) && (k.2 < bk.2 || (k.2 == bk.2 && k.3 < bk.3)));
        if better {
            best = i;
            bk = k;
        }
    }
    let (i, j) = (best, best + 1);
    let (keep, drop) = if st.cover[j].thickness_m > st.cover[i].thickness_m { (j, i) } else { (i, j) };
    let from = st.cover[drop];
    absorb(&mut st.cover[keep], &from);
    st.cover.remove(drop);
    if let Some(o)=observer{o.capture(17,crate::layer_evolution::Mapping::Merge(keep,drop),st,b);}
    bk.0 == 1
}

/// `parcel-stacks-v2`: the t = 0 platform cover of a seeded parcel whose
/// cell stands `elev_above_sea_m` above the sea datum: thickness by the
/// declared elevation ramp (capped at half the column), bottom-up
/// sandstone, shale, carbonate with ages split evenly over the declared
/// span; relabels the top of the basement (volume neutral). Returns the
/// cover thickness.
pub fn seed_platform_cover(st: &mut StrataStack, elev_above_sea_m: f64, v: &StrataV2Params) -> f64 {
    let ramp = ((v.platform_zero_elev_m - elev_above_sea_m) / (v.platform_zero_elev_m - v.platform_full_elev_m)).clamp(0.0, 1.0);
    let t = (v.platform_cover_max_m * ramp).min(0.5 * st.basement_m.max(0.0));
    if !(t > 0.0) {
        return 0.0;
    }
    let (a0, a1) = (v.platform_age_base_my, v.platform_age_top_my);
    let da = (a1 - a0) / 3.0;
    let parts = [
        (LITH_SANDSTONE, ENV_FLUVIAL, v.platform_sandstone_fraction),
        (LITH_SHALE, ENV_SHELF, 1.0 - v.platform_sandstone_fraction - v.platform_carbonate_fraction),
        (LITH_CARBONATE, ENV_SHELF, v.platform_carbonate_fraction),
    ];
    let mut used = 0.0;
    for (i, &(lith, env, f)) in parts.iter().enumerate() {
        let th = t * f;
        if th > 0.0 {
            st.cover.push(StrataLayer {
                thickness_m: th,
                age_base_my: a0 + da * i as f64,
                age_top_my: a0 + da * (i + 1) as f64,
                max_burial_m: 0.0,
                lith,
                env,
                grade: 0,
                flags: if organic(lith, env) { STRATA_FLAG_ORGANIC_RICH } else { 0 },
            });
            used += th;
        }
    }
    st.basement_m -= used;
    // Burial depths within the new cover.
    let mut depth = 0.0;
    for l in st.cover.iter_mut().rev() {
        depth += l.thickness_m;
        l.max_burial_m = l.max_burial_m.max(depth);
    }
    used
}

fn retain_positive_observed(st:&mut StrataStack,b:f64,observer:Option<&mut crate::layer_evolution::Observer>){
    let surviving=observer.as_ref().map(|_|st.cover.iter().enumerate().filter(|(_,l)|l.thickness_m>0.).map(|(i,_)|i).collect::<Vec<_>>());
    st.cover.retain(|l|l.thickness_m>0.);
    if let Some(o)=observer{o.capture(16,crate::layer_evolution::Mapping::Filter(surviving.unwrap()),st,b);}
}

/// Apply one parcel's amounts in the declared order. `total_after_m` is
/// the parcel column (thickness + sediment) after the write-back; the
/// basement is set so the stack equals it.
#[allow(clippy::too_many_arguments)]
pub fn apply_parcel(
    st: &mut StrataStack,
    a: &ParcelAmounts,
    total_after_m: f64,
    area: f64,
    p: &StrataParams,
    t0: f64,
    t1: f64,
    log: &mut StrataStepLog,
) {
    apply_parcel_inner(st,a,total_after_m,area,p,t0,t1,log,None)
}

/// Record actual native mutations under a conservative pre-mutation byte budget.
#[allow(clippy::too_many_arguments)]
pub fn apply_parcel_recorded(st:&mut StrataStack,a:&ParcelAmounts,total_after_m:f64,area:f64,p:&StrataParams,t0:f64,t1:f64,log:&mut StrataStepLog,byte_limit:u64)->Result<crate::layer_evolution::LayerEvolution,String>{
    let mut observer=crate::layer_evolution::Observer::new(st,p,byte_limit)?;
    apply_parcel_inner(st,a,total_after_m,area,p,t0,t1,log,Some(&mut observer));
    Ok(observer.finish(st))
}

/// Same native mutations with an opt-in passive formation mixture observer.
#[allow(clippy::too_many_arguments)]
pub fn apply_parcel_provenance(st:&mut StrataStack,a:&ParcelAmounts,total:f64,area:f64,p:&StrataParams,t0:f64,t1:f64,log:&mut StrataStepLog,provenance:pf_contracts::layer_provenance::StackProvenance,carrier:u64,retain_trace:bool)->Result<(pf_contracts::layer_provenance::StackProvenance,crate::layer_evolution::LayerEvolution),String>{
    crate::layer_provenance::validate(&provenance,st)?;
    let mut o=crate::layer_evolution::Observer::new(st,p,crate::layer_evolution::byte_bound(p.max_layers))?;
    o.retain_trace=retain_trace;
    o.provenance=Some(crate::layer_provenance::Recorder{state:provenance,carrier,amounts:*a,params:*p,t0,t1});
    apply_parcel_inner(st,a,total,area,p,t0,t1,log,Some(&mut o));
    let result=o.provenance.take().unwrap().state;
    crate::layer_provenance::validate(&result,st)?;
    Ok((result,o.finish(st)))
}

/// Birth records use the actual birth mutation trace. Pre-conversion basement
/// is unresolved; its earlier oceanic formation is not retained by this stack.
pub fn birth_provenance(st:&StrataStack,trace:&crate::layer_evolution::LayerEvolution,arc:f64,class:u8,p:&StrataParams,t:f64,carrier:u64)->pf_contracts::layer_provenance::StackProvenance {
    let base=trace.changes.first().map_or(st,|c|&c.before);
    let mut state=crate::layer_provenance::initial(base,carrier);
    for c in &mut state.basement.contributions {
        c.setting=pf_contracts::layer_provenance::Setting::Unknown;
        c.episode=format!("carrier:{carrier}:unresolved-preconversion:{:016x}",t.to_bits());
    }
    let mut o=crate::layer_provenance::Recorder{state,carrier,amounts:ParcelAmounts{arc,class,..Default::default()},params:*p,t0:t,t1:t};
    for c in &trace.changes{o.capture(c);}
    crate::layer_provenance::validate(&o.state,st).expect("birth provenance closure");o.state
}

#[allow(clippy::too_many_arguments)]
fn apply_parcel_inner(
    st: &mut StrataStack,
    a: &ParcelAmounts,
    total_after_m: f64,
    area: f64,
    p: &StrataParams,
    t0: f64,
    t1: f64,
    log: &mut StrataStepLog,
    mut observer:Option<&mut crate::layer_evolution::Observer>,
) {
    use crate::layer_evolution::Mapping;
    let mut b = st.basement_m;
    macro_rules! note {($kind:expr,$mapping:expr)=>{if let Some(o)=observer.as_deref_mut(){o.capture($kind,$mapping,st,b);}};}
    // 1. Erosion: truncate from the top.
    if a.erosion > 0.0 {
        let mut e = a.erosion;
        let mut truncated = false;
        while e > 0.0 {
            let Some(top) = st.cover.last_mut() else { break };
            if top.thickness_m <= e {
                e -= top.thickness_m;
                log.eroded_m3[top.lith as usize] += top.thickness_m * area;
                let removed=st.cover.len()-1;
                st.cover.pop();
                note!(0,Mapping::Remove(removed,None));
            } else {
                top.thickness_m -= e;
                log.eroded_m3[top.lith as usize] += e * area;
                top.flags |= STRATA_FLAG_TRUNCATED_TOP;
                truncated = true;
                e = 0.0;
                note!(1,Mapping::Identity);
            }
        }
        if e > 0.0 {
            b -= e;
            log.eroded_m3[st.basement_lith as usize] += e * area;
            note!(2,Mapping::Identity);
        }
        if truncated {
            if let Some(top) = st.cover.last().copied() {
                if top.thickness_m < p.min_layer_m {
                    let removed=st.cover.len()-1;
                    st.cover.pop();
                    match st.cover.last_mut() {
                        Some(below) => absorb(below, &top),
                        None => b += top.thickness_m,
                    }
                    note!(3,if st.cover.is_empty(){Mapping::BaseAbsorb(removed)}else{Mapping::Remove(removed,Some(removed-1))});
                    log.merges_min += 1;
                }
            }
        }
    }
    // 2. Thinning: scale every layer.
    if a.thinning > 0.0 {
        let t = b + st.cover_m();
        if t > 0.0 {
            let f = ((t - a.thinning) / t).max(0.0);
            for l in st.cover.iter_mut() {
                l.thickness_m *= f;
            }
            b *= f;
            log.thinned_m3 += (t - (b + st.cover_m())) * area;
            note!(4,Mapping::Identity);
        }
    }
    // 3. Basement removal (channel outflow, delamination).
    if a.short_neg > 0.0 {
        b -= a.short_neg;
        log.shortening_removed_m3 += a.short_neg * area;
        note!(5,Mapping::Identity);
    }
    // 3b. parcel-stacks-v2: the non-placed shortening-lineage change
    // (plateau, channel flow, delamination) is basement only (signed).
    if a.belt_rest != 0.0 {
        b += a.belt_rest;
        log.belt_basement_m3 += a.belt_rest * area;
        note!(6,Mapping::Identity);
    }
    // 4. Plutons at the basement top, below a crystalline roof of at least
    // `pluton_roof_m` (with this step's volcanics) moved up from the
    // basement (PREREG revision 1; volume neutral).
    let pl = a.arc * (1.0 - p.arc_extrusive_fraction) + a.short_pos * p.shortening_pluton_fraction;
    if pl > 0.0 {
        let above = st.cover_m() + a.arc * p.arc_extrusive_fraction + a.plume * p.plume_extrusive_fraction;
        let roof = (p.pluton_roof_m - above).min(b).max(0.0);
        if roof > 0.0 {
            let inserted=!st.cover.first().is_some_and(|f|f.lith==LITH_BASEMENT&&f.env==ENV_BASEMENT);
            match st.cover.first_mut() {
                Some(f) if f.lith == LITH_BASEMENT && f.env == ENV_BASEMENT => f.thickness_m += roof,
                _ => st.cover.insert(
                    0,
                    StrataLayer {
                        thickness_m: roof,
                        age_base_my: st.basement_age_my,
                        age_top_my: st.basement_age_my,
                        max_burial_m: 0.0,
                        lith: LITH_BASEMENT,
                        env: ENV_BASEMENT,
                        grade: 0,
                        flags: 0,
                    },
                ),
            }
            b -= roof;
            note!(7,if inserted{Mapping::Insert(0,true)}else{Mapping::BaseTo(0)});
        }
        insert_bottom_observed(st, &mut b, LITH_GRANITE, ENV_PLUTON, pl, t0, t1, p, log,observer.as_deref_mut(),8);
        log.pluton_m3 += pl * area;
    }
    // 5. Basement additions.
    let pb = a.plume * (1.0 - p.plume_extrusive_fraction);
    if pb > 0.0 {
        b += pb;
        log.plume_basement_m3 += pb * area;
        note!(9,Mapping::Identity);
    }
    if a.basement_add > 0.0 {
        b += a.basement_add;
        log.basement_add_m3 += a.basement_add * area;
        note!(10,Mapping::Identity);
    }
    // 6. Shortening: pure-shear thickening of every layer.
    let s = a.short_pos * (1.0 - p.shortening_pluton_fraction);
    if s > 0.0 {
        let t = b + st.cover_m();
        if t > 0.0 {
            let f = (t + s) / t;
            for l in st.cover.iter_mut() {
                l.thickness_m *= f;
            }
            b *= f;
            log.shortening_scale_m3 += (b + st.cover_m() - t) * area;
        } else {
            b += s;
            log.shortening_scale_m3 += s * area;
        }
    }
    note!(11,Mapping::Identity);
    // 7. Volcanics on top.
    let v = a.arc * p.arc_extrusive_fraction + a.plume * p.plume_extrusive_fraction;
    if v > 0.0 {
        push_top_observed(st, &mut b, LITH_BASALT, ENV_VOLCANIC, v, t0, t1, p, log,observer.as_deref_mut(),12);
        log.volcanic_m3 += v * area;
    }
    // 8. Sediment on top.
    if a.deposit > 0.0 {
        let (l, e) = DEPOSIT_CLASSES[a.class as usize];
        push_top_observed(st, &mut b, l, e, a.deposit, t0, t1, p, log,observer.as_deref_mut(),13);
        log.deposited_m3[a.class as usize] += a.deposit * area;
        if a.warm_shelf {
            log.warm_shelf_deposit_m3 += a.deposit * area;
        }
    }
    // 9. The residual into the basement (declared basement-only rule).
    let r = total_after_m - (b + st.cover_m());
    b += r;
    log.basement_flow_m3 += r * area;
    log.basement_flow_abs_m3 += r.abs() * area;
    note!(14,Mapping::Identity);
    // 10. A negative basement takes its deficit from the cover, bottom-up.
    if b < 0.0 {
        let mut need = -b;
        for l in st.cover.iter_mut() {
            if need <= 0.0 {
                break;
            }
            let take = need.min(l.thickness_m);
            l.thickness_m -= take;
            need -= take;
        }
        log.underflow_m3 += (-b - need.max(0.0)) * area;
        b = 0.0;
        note!(15,Mapping::Identity);
        retain_positive_observed(st,b,observer.as_deref_mut());
    }
    retain_positive_observed(st,b,observer.as_deref_mut());
    // 11. Capacity merges.
    while st.cover.len() + 1 > p.max_layers {
        if p.v2.is_some() {
            if merge_pair_v2_observed(st, p.unconformity_gap_my,b,observer.as_deref_mut()) {
                log.merges_unconformity += 1;
            }
        } else {
            merge_thinnest_observed(st,b,observer.as_deref_mut());
        }
        log.merges_cap += 1;
    }
    // 12. Burial and metamorphism (top down).
    let mut depth = 0.0;
    for l in st.cover.iter_mut().rev() {
        depth += l.thickness_m;
        l.max_burial_m = l.max_burial_m.max(depth);
        if depth > p.metamorphic_depth_m && l.lith != LITH_GRANITE && l.lith != LITH_OCEANIC {
            if l.lith != LITH_METAMORPHIC {
                log.metamorphosed_m3 += l.thickness_m * area;
                l.lith = LITH_METAMORPHIC;
                l.grade = l.grade.max(1);
            }
            if depth > p.high_grade_depth_m {
                l.grade = 2;
            }
        }
    }
    note!(18,Mapping::Identity);
    st.basement_m = b;
    if a.land {
        st.last_exposure_my = t1;
    }
    note!(19,Mapping::Identity);
}

/// Initial stack of a parcel born this step from cell column (thickness,
/// sediment, arc marker): converted island arcs get arc pluton and
/// volcanics over an oceanic-igneous basement; other births (rift) a
/// crystalline basement; the cell's sediment on top in the cell's class.
#[allow(clippy::too_many_arguments)]
pub fn birth_stack(thickness_m: f64, sediment_m: f64, arc_added_m: f64, class: u8, p: &StrataParams, t_my: f64, area: f64, log: &mut StrataStepLog) -> StrataStack {
    birth_stack_inner(thickness_m,sediment_m,arc_added_m,class,p,t_my,area,log,None).expect("disabled layer observer cannot fail").0
}
#[allow(clippy::too_many_arguments)]
pub fn birth_stack_recorded(thickness_m:f64,sediment_m:f64,arc_added_m:f64,class:u8,p:&StrataParams,t_my:f64,area:f64,log:&mut StrataStepLog,byte_limit:u64)->Result<(StrataStack,crate::layer_evolution::LayerEvolution),String>{
    let (stack,trace)=birth_stack_inner(thickness_m,sediment_m,arc_added_m,class,p,t_my,area,log,Some(byte_limit))?;
    Ok((stack,trace.expect("requested layer observer")))
}
#[allow(clippy::too_many_arguments)]
fn birth_stack_inner(thickness_m: f64, sediment_m: f64, arc_added_m: f64, class: u8, p: &StrataParams, t_my: f64, area: f64, log: &mut StrataStepLog,record_budget:Option<u64>) -> Result<(StrataStack,Option<crate::layer_evolution::LayerEvolution>),String> {
    let h = thickness_m.max(0.0);
    let arc = arc_added_m.max(0.0).min(h);
    let mut st = StrataStack::basement_only(h - arc, if arc > 0.0 { LITH_OCEANIC } else { LITH_BASEMENT }, t_my);
    let mut observer=record_budget.map(|limit|crate::layer_evolution::Observer::new(&st,p,limit)).transpose()?;
    let mut scratch = StrataStepLog::default();
    let mut b = st.basement_m;
    insert_bottom_observed(&mut st, &mut b, LITH_GRANITE, ENV_PLUTON, arc * (1.0 - p.arc_extrusive_fraction), t_my, t_my, p, &mut scratch,observer.as_mut(),8);
    push_top_observed(&mut st, &mut b, LITH_BASALT, ENV_VOLCANIC, arc * p.arc_extrusive_fraction, t_my, t_my, p, &mut scratch,observer.as_mut(),12);
    let (l, e) = DEPOSIT_CLASSES[class as usize];
    push_top_observed(&mut st, &mut b, l, e, sediment_m.max(0.0), t_my, t_my, p, &mut scratch,observer.as_mut(),13);
    st.basement_m = b;
    st.lith_volumes(area, &mut log.births_m3);
    let trace=observer.map(|o|o.finish(&st));
    Ok((st,trace))
}

/// Inputs of the deposit environment classifier (per cell).
pub struct EnvInputs<'a> {
    pub elevation_m: &'a [f64],
    pub sea_m: f64,
    pub ocean: &'a [bool],
    pub fill_level_m: &'a [f64],
    pub endorheic: &'a [bool],
    pub flow_to: &'a [Option<usize>],
    pub contributing_area_m2: &'a [f64],
    /// Deposition rate this step (m/My).
    pub deposit_rate_m_per_my: &'a [f64],
}

/// The deposit class of every cell (PREREG.txt environment rules; the
/// latitude fallback A-STRATA-LATITUDE-CLIMATE). Also returns the land
/// (subaerial, not a lake) mask.
pub fn classify_cells(grid: &GridSpec, x: &EnvInputs, p: &StrataParams) -> Result<(Vec<u8>, Vec<bool>), KinematicError> {
    let (class, land, _) = classify_cells_with_shelf(grid, x, p)?;
    Ok((class, land))
}

/// `classify_cells` plus, under `parcel-stacks-v2`, the declared sub-cell
/// shelf-slope ramp (PREREG C2: a marine cell deeper than the shelf depth
/// takes the shelf facies when its shelf fraction min(f_dist, f_depth)
/// reaches the majority fraction) and the warm-shelf mask (empty under v1).
pub fn classify_cells_with_shelf(grid: &GridSpec, x: &EnvInputs, p: &StrataParams) -> Result<(Vec<u8>, Vec<bool>, Vec<bool>), KinematicError> {
    let n = grid.cells();
    let nlon = grid.nlon();
    let g = |e: pf_grid::GridError| KinematicError::Grid(e.to_string());
    let mut lat = vec![0.0f64; grid.nlat()];
    for (i, l) in lat.iter_mut().enumerate() {
        *l = grid.center_lat_deg(i).map_err(g)?.abs();
    }
    let lake: Vec<bool> = (0..n).map(|c| !x.ocean[c] && x.fill_level_m[c] > x.elevation_m[c] + 1.0).collect();
    // River mouths: marine cells a land neighbour with a large drainage
    // area flows into.
    let mut mouth = vec![false; n];
    for c in 0..n {
        if x.ocean[c] || lake[c] {
            continue;
        }
        if let Some(t) = x.flow_to[c] {
            if x.ocean[t] && x.contributing_area_m2[c] >= p.river_mouth_area_m2 {
                mouth[t] = true;
            }
        }
    }
    let mut class = vec![0u8; n];
    let mut land = vec![false; n];
    // parcel-stacks-v2: the sub-cell shelf fraction of every marine cell.
    let mut warm_shelf = if p.v2.is_some() { vec![false; n] } else { Vec::new() };
    let shelf_frac: Vec<f64> = match &p.v2 {
        None => Vec::new(),
        Some(v) => {
            let areas = crate::kinematic::cell_areas(grid)?;
            let mut f = vec![0.0f64; n];
            for c in 0..n {
                if !x.ocean[c] {
                    continue;
                }
                let d = x.sea_m - x.elevation_m[c];
                if d <= p.shelf_depth_m {
                    f[c] = 1.0;
                    continue;
                }
                let (i, j) = (c / nlon, c % nlon);
                let mut coast = false;
                for (ti, tj, _) in grid.neighbors(i, j).map_err(g)? {
                    coast |= !x.ocean[ti * nlon + tj];
                }
                let f_dist = if coast { (v.shelf_width_m / areas[c].max(1.0).sqrt()).min(1.0) } else { 0.0 };
                let f_depth = ((v.slope_foot_depth_m - d) / (v.slope_foot_depth_m - 0.5 * p.shelf_depth_m)).clamp(0.0, 1.0);
                f[c] = f_dist.min(f_depth);
            }
            f
        }
    };
    for c in 0..n {
        let (i, j) = (c / nlon, c % nlon);
        let la = lat[i];
        let arid = la >= p.arid_lat_min_deg && la <= p.arid_lat_max_deg;
        let humid = la < p.humid_lat_max_deg || (la >= p.humid_midlat_min_deg && la <= p.humid_midlat_max_deg);
        let warm = la <= p.carbonate_max_paleolat_deg;
        let h = x.elevation_m[c];
        class[c] = if lake[c] {
            if arid && x.endorheic[c] {
                4
            } else {
                3
            }
        } else if x.ocean[c] {
            let depth = x.sea_m - h;
            let deep = match &p.v2 {
                None => depth > p.shelf_depth_m,
                Some(v) => shelf_frac[c] < v.shelf_majority_fraction,
            };
            if deep {
                8
            } else {
                if warm && !warm_shelf.is_empty() {
                    warm_shelf[c] = true;
                }
                if mouth[c] {
                    5
                } else if warm && x.deposit_rate_m_per_my[c] <= p.carbonate_max_supply_m_per_my {
                    6
                } else {
                    7
                }
            }
        } else {
            land[c] = true;
            let (mut lo, mut hi, mut coastal) = (h, h, false);
            for (ti, tj, _) in grid.neighbors(i, j).map_err(g)? {
                let k = ti * nlon + tj;
                lo = lo.min(x.elevation_m[k]);
                hi = hi.max(x.elevation_m[k]);
                coastal |= x.ocean[k];
            }
            if humid && coastal && h - x.sea_m <= p.coal_max_elevation_m {
                2
            } else if hi - lo >= p.conglomerate_relief_m {
                1
            } else {
                0
            }
        };
    }
    Ok((class, land, warm_shelf))
}

/// The record of every cell, top down to `depth_m` (`f64::INFINITY` for the
/// full column): the dominant parcel's stack on continental cells (`dom`),
/// pelagic shale over oceanic igneous crust elsewhere.
#[allow(clippy::too_many_arguments)]
pub fn publish(
    grid: &GridSpec,
    t_my: f64,
    stacks: &[StrataStack],
    dom: &[Option<usize>],
    sediment_m: &[f64],
    thickness_m: &[f64],
    birth_t_my: &[f64],
    p: &StrataParams,
    depth_m: f64,
    assumptions: Vec<AssumptionV0>,
) -> StratigraphyDomainV0 {
    let n = grid.cells();
    let mut d = StratigraphyDomainV0 {
        schema: STRATIGRAPHY_SCHEMA_V0.into(),
        t_my,
        nlat: grid.nlat() as u32,
        nlon: grid.nlon() as u32,
        publish_depth_m: if depth_m.is_finite() { depth_m } else { 1.0e9 },
        unconformity_gap_my: p.unconformity_gap_my,
        usage_limits: if p.v2.is_some() { STRATIGRAPHY_PREVIEW_USAGE_LIMITS_V2.into() } else { STRATIGRAPHY_PREVIEW_USAGE_LIMITS.into() },
        layer_offset: Vec::with_capacity(n + 1),
        lith: Vec::new(),
        env: Vec::new(),
        grade: Vec::new(),
        flags: Vec::new(),
        thickness_m: Vec::new(),
        age_base_my: Vec::new(),
        age_top_my: Vec::new(),
        max_burial_m: Vec::new(),
        cover_m: Vec::with_capacity(n),
        layer_count: Vec::with_capacity(n),
        codes: StratigraphyCodesV0::v0(),
        assumptions,
        rock_uplift_m: None,
        uplift_grad_east: None,
        uplift_grad_north: None,
    };
    d.layer_offset.push(0);
    // (lith, env, grade, flags, thickness, age base, age top, burial), top down.
    let mut col: Vec<(u8, u8, u8, u8, f64, f64, f64, f64)> = Vec::new();
    for c in 0..n {
        col.clear();
        match dom[c] {
            Some(k) => {
                let st = &stacks[k];
                let mut depth = 0.0;
                for l in st.cover.iter().rev() {
                    col.push((l.lith, l.env, l.grade, l.flags, l.thickness_m, l.age_base_my, l.age_top_my, l.max_burial_m.max(depth + l.thickness_m)));
                    depth += l.thickness_m;
                }
                let env = if st.basement_lith == LITH_OCEANIC { ENV_OCEANIC } else { ENV_BASEMENT };
                col.push((st.basement_lith, env, 0, STRATA_FLAG_BASEMENT, st.basement_m, st.basement_age_my, st.basement_age_my, depth + st.basement_m));
                d.cover_m.push(st.cover_m() as f32);
                d.layer_count.push((st.cover.len() + 1) as u8);
            }
            None => {
                let s = sediment_m[c].max(0.0);
                if s > 0.0 {
                    col.push((LITH_SHALE, ENV_DEEP_MARINE, 0, STRATA_FLAG_ORGANIC_RICH, s, birth_t_my[c], t_my, s));
                }
                col.push((LITH_OCEANIC, ENV_OCEANIC, 0, STRATA_FLAG_BASEMENT, thickness_m[c].max(0.0), birth_t_my[c], birth_t_my[c], s + thickness_m[c].max(0.0)));
                d.cover_m.push(s as f32);
                d.layer_count.push(if s > 0.0 { 2 } else { 1 });
            }
        }
        // Unconformities (derived): the gap below a layer's base.
        for i in 0..col.len() {
            if i + 1 < col.len() && col[i].5 - col[i + 1].6 > p.unconformity_gap_my {
                col[i].3 |= STRATA_FLAG_UNCONFORMITY_BASE;
            }
        }
        let mut depth = 0.0;
        for &(lith, env, grade, flags, th, ab, at, mb) in &col {
            if depth >= depth_m {
                break;
            }
            let (keep, cut) = if depth + th > depth_m { (depth_m - depth, true) } else { (th, false) };
            d.lith.push(lith);
            d.env.push(env);
            d.grade.push(grade);
            d.flags.push(if cut { flags | STRATA_FLAG_CUT_AT_DEPTH } else { flags });
            d.thickness_m.push(keep as f32);
            d.age_base_my.push(ab as f32);
            d.age_top_my.push(at as f32);
            d.max_burial_m.push(mb as f32);
            depth += th;
        }
        d.layer_offset.push(d.lith.len() as u32);
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> StrataParams {
        StrataParams {
            max_layers: 12,
            merge_gap_my: 2.0,
            unconformity_gap_my: 5.0,
            min_layer_m: 10.0,
            trace_m: 0.01,
            arc_extrusive_fraction: 0.3,
            plume_extrusive_fraction: 0.2,
            shortening_pluton_fraction: 0.1,
            pluton_roof_m: 4_000.0,
            metamorphic_depth_m: 10_000.0,
            high_grade_depth_m: 20_000.0,
            shelf_depth_m: 200.0,
            carbonate_max_paleolat_deg: 30.0,
            carbonate_max_supply_m_per_my: 50.0,
            arid_lat_min_deg: 15.0,
            arid_lat_max_deg: 35.0,
            humid_lat_max_deg: 15.0,
            humid_midlat_min_deg: 40.0,
            humid_midlat_max_deg: 65.0,
            coal_max_elevation_m: 200.0,
            conglomerate_relief_m: 1000.0,
            river_mouth_area_m2: 1.0e11,
            publish_depth_m: 3000.0,
            basement_min_m: 5000.0,
            v2: None,
        }
    }

    fn closes(before: f64, st: &StrataStack, log: &StrataStepLog) {
        let after = st.total_m();
        assert!((after - (before + log.ops_m3())).abs() < 1e-6, "{after} {before} {}", log.ops_m3());
    }

    #[test]
    fn deposition_erosion_and_residual_close() {
        let p = params();
        let mut st = StrataStack::basement_only(35_000.0, LITH_BASEMENT, -500.0);
        let mut log = StrataStepLog::default();
        // Carbonate, then shale, then erosion into the carbonate.
        let a = ParcelAmounts { deposit: 300.0, class: 6, ..Default::default() };
        apply_parcel(&mut st, &a, 35_300.0, 1.0, &p, 0.0, 2.0, &mut log);
        let a = ParcelAmounts { deposit: 200.0, class: 7, ..Default::default() };
        apply_parcel(&mut st, &a, 35_500.0, 1.0, &p, 2.0, 4.0, &mut log);
        assert_eq!(st.cover.len(), 2);
        let before = st.total_m();
        let mut l2 = StrataStepLog::default();
        let a = ParcelAmounts { erosion: 250.0, ..Default::default() };
        // The column says 240 m left the parcel: 10 m of residual.
        apply_parcel(&mut st, &a, before - 240.0, 1.0, &p, 4.0, 6.0, &mut l2);
        closes(before, &st, &l2);
        assert_eq!(st.cover.len(), 1);
        assert_eq!(st.cover[0].lith, LITH_CARBONATE);
        assert!((st.cover[0].thickness_m - 250.0).abs() < 1e-9);
        assert!(st.cover[0].flags & STRATA_FLAG_TRUNCATED_TOP != 0);
        assert!((l2.basement_flow_m3 - 10.0).abs() < 1e-9);
        assert!((st.total_m() - (before - 240.0)).abs() < 1e-9);
    }

    #[test]
    fn same_class_within_gap_extends_and_capacity_merges() {
        let p = params();
        let mut st = StrataStack::basement_only(30_000.0, LITH_BASEMENT, 0.0);
        let mut log = StrataStepLog::default();
        let mut tot = st.total_m();
        for step in 0..40 {
            let class = (step % 3) as u8 * 3; // 0, 3, 6: three alternating classes
            let a = ParcelAmounts { deposit: 20.0 + step as f64, class, ..Default::default() };
            tot += a.deposit;
            let before = st.total_m();
            let mut l = StrataStepLog::default();
            apply_parcel(&mut st, &a, tot, 1.0, &p, 2.0 * step as f64, 2.0 * step as f64 + 2.0, &mut l);
            closes(before, &st, &l);
            log.merges_cap += l.merges_cap;
            assert!(st.cover.len() + 1 <= p.max_layers);
        }
        assert!(log.merges_cap > 0);
        // Same class twice in a row extends the top layer.
        let n0 = st.cover.len();
        let top_class = st.cover.last().map(|l| (l.lith, l.env)).unwrap();
        let class = DEPOSIT_CLASSES.iter().position(|c| *c == top_class).unwrap() as u8;
        tot += 5.0;
        apply_parcel(&mut st, &ParcelAmounts { deposit: 5.0, class, ..Default::default() }, tot, 1.0, &p, 80.0, 82.0, &mut log);
        assert_eq!(st.cover.len(), n0);
    }

    #[test]
    fn arc_shortening_and_metamorphism() {
        let p = params();
        let mut st = StrataStack::basement_only(30_000.0, LITH_BASEMENT, 0.0);
        let mut log = StrataStepLog::default();
        let a = ParcelAmounts { deposit: 6_000.0, class: 7, ..Default::default() };
        apply_parcel(&mut st, &a, 36_000.0, 1.0, &p, 0.0, 2.0, &mut log);
        let a = ParcelAmounts { arc: 10_000.0, ..Default::default() };
        let before = st.total_m();
        let mut l = StrataStepLog::default();
        apply_parcel(&mut st, &a, 46_000.0, 1.0, &p, 2.0, 4.0, &mut l);
        closes(before, &st, &l);
        // pluton at the bottom (no roof needed: 6 km shale + 3 km volcanics
        // above), basalt on top of the shale.
        assert_eq!(st.cover[0].lith, LITH_GRANITE);
        assert!((st.cover[0].thickness_m - 7_000.0).abs() < 1e-9);
        assert_eq!(st.top_lith(), LITH_BASALT);
        // Shortening doubles the column: the shale base (3 km volcanic + 6 km
        // shale = 9 km) goes below 10 km and becomes metamorphic.
        let a = ParcelAmounts { short_pos: 46_000.0, ..Default::default() };
        let before = st.total_m();
        let mut l = StrataStepLog::default();
        apply_parcel(&mut st, &a, 92_000.0, 1.0, &p, 4.0, 6.0, &mut l);
        closes(before, &st, &l);
        assert!(st.cover.iter().any(|x| x.lith == LITH_METAMORPHIC && x.env == ENV_SHELF));
        assert!(l.metamorphosed_m3 > 0.0);
        assert!((st.total_m() - 92_000.0).abs() < 1e-6);
    }

    #[test]
    fn pluton_on_bare_basement_gets_a_crystalline_roof() {
        let p = params();
        let mut st = StrataStack::basement_only(35_000.0, LITH_BASEMENT, -900.0);
        let mut log = StrataStepLog::default();
        let a = ParcelAmounts { short_pos: 5_000.0, ..Default::default() };
        apply_parcel(&mut st, &a, 40_000.0, 1.0, &p, 0.0, 2.0, &mut log);
        closes(35_000.0, &st, &log);
        assert_eq!(st.cover[0].lith, LITH_GRANITE);
        assert_eq!(st.top_lith(), LITH_BASEMENT);
        assert!(st.cover[1].thickness_m >= 4_000.0 - 1e-9);
        // Unroofing: erode 5 km and the pluton is the surface.
        let before = st.total_m();
        let mut l2 = StrataStepLog::default();
        apply_parcel(&mut st, &ParcelAmounts { erosion: 5_000.0, ..Default::default() }, before - 5_000.0, 1.0, &p, 2.0, 4.0, &mut l2);
        closes(before, &st, &l2);
        assert_eq!(st.top_lith(), LITH_GRANITE);
    }

    fn v2() -> StrataParams {
        StrataParams {
            max_layers: 16,
            v2: Some(StrataV2Params {
                shelf_width_m: 80_000.0,
                slope_foot_depth_m: 4_000.0,
                shelf_majority_fraction: 0.5,
                platform_cover_max_m: 2_000.0,
                platform_full_elev_m: 0.0,
                platform_zero_elev_m: 600.0,
                platform_sandstone_fraction: 0.3,
                platform_carbonate_fraction: 0.3,
                platform_age_base_my: -600.0,
                platform_age_top_my: -50.0,
            }),
            ..params()
        }
    }

    #[test]
    fn v2_platform_cover_by_elevation_is_volume_neutral() {
        let p = v2();
        let v = p.v2.unwrap();
        p.validate().unwrap();
        let mut low = StrataStack::basement_only(35_000.0, LITH_BASEMENT, -2000.0);
        assert!((seed_platform_cover(&mut low, -50.0, &v) - 2_000.0).abs() < 1e-9);
        assert!((low.total_m() - 35_000.0).abs() < 1e-9);
        assert_eq!(low.cover.iter().map(|l| l.lith).collect::<Vec<_>>(), vec![LITH_SANDSTONE, LITH_SHALE, LITH_CARBONATE]);
        assert_eq!(low.top_lith(), LITH_CARBONATE);
        assert!((low.cover[0].thickness_m - 600.0).abs() < 1e-9 && (low.cover[1].thickness_m - 800.0).abs() < 1e-9);
        let mut mid = StrataStack::basement_only(35_000.0, LITH_BASEMENT, -2000.0);
        assert!((seed_platform_cover(&mut mid, 300.0, &v) - 1_000.0).abs() < 1e-9);
        let mut high = StrataStack::basement_only(35_000.0, LITH_BASEMENT, -2000.0);
        assert_eq!(seed_platform_cover(&mut high, 800.0, &v), 0.0);
        assert!(high.cover.is_empty());
    }

    #[test]
    fn v2_merge_keeps_unconformities_and_merges_old_similar_first() {
        let p = v2();
        let mut st = StrataStack::basement_only(30_000.0, LITH_BASEMENT, -2000.0);
        // 15 cover layers: an old pair of sandstone/conglomerate (same group),
        // then alternating carbonate / shale (cost 2 pairs), a 50 My gap before
        // the top four.
        let mut t = 0.0;
        let mut push = |st: &mut StrataStack, lith: u8, th: f64, gap: f64| {
            t += gap;
            st.cover.push(StrataLayer { thickness_m: th, age_base_my: t, age_top_my: t + 2.0, max_burial_m: 0.0, lith, env: ENV_SHELF, grade: 0, flags: 0 });
            t += 2.0;
        };
        push(&mut st, LITH_SANDSTONE, 300.0, 0.0);
        push(&mut st, LITH_CONGLOMERATE, 200.0, 0.0);
        for k in 0..9 {
            push(&mut st, if k % 2 == 0 { LITH_CARBONATE } else { LITH_SHALE }, 100.0, 0.0);
        }
        push(&mut st, LITH_SHALE, 50.0, 50.0);
        for k in 0..3 {
            push(&mut st, if k % 2 == 0 { LITH_CARBONATE } else { LITH_SHALE }, 40.0, 0.0);
        }
        assert_eq!(st.cover.len(), 15);
        let before = st.total_m();
        // One more deposit forces one capacity merge.
        let mut log = StrataStepLog::default();
        let a = ParcelAmounts { deposit: 20.0, class: 0, ..Default::default() };
        apply_parcel(&mut st, &a, before + 20.0, 1.0, &p, t + 1.0, t + 3.0, &mut log);
        closes(before, &st, &log);
        assert_eq!(log.merges_cap, 1);
        assert_eq!(log.merges_unconformity, 0);
        // The old same-group pair merged (the thicker sandstone kept).
        assert_eq!(st.cover[0].lith, LITH_SANDSTONE);
        assert!((st.cover[0].thickness_m - 500.0).abs() < 1e-9);
        assert!(st.cover[0].flags & STRATA_FLAG_MERGED != 0);
        // The unconformity under the young shale survives.
        assert!(st.cover.windows(2).any(|w| w[1].age_base_my - w[0].age_top_my > 5.0));
    }

    #[test]
    fn v2_belt_rest_goes_to_the_basement_and_closes() {
        let p = v2();
        let mut st = StrataStack::basement_only(30_000.0, LITH_BASEMENT, -2000.0);
        let mut log = StrataStepLog::default();
        // 3 km of channel inflow only (no placed supply): no pluton, no
        // thickening of the cover, basement +3 km.
        let a = ParcelAmounts { belt_rest: 3_000.0, ..Default::default() };
        apply_parcel(&mut st, &a, 33_000.0, 1.0, &p, 0.0, 2.0, &mut log);
        closes(30_000.0, &st, &log);
        assert!(st.cover.is_empty());
        assert!((log.belt_basement_m3 - 3_000.0).abs() < 1e-9);
        assert!(log.basement_flow_abs_m3 < 1e-9);
    }

    #[test]
    fn deep_erosion_and_underflow_close() {
        let p = params();
        let mut st = StrataStack::basement_only(1_000.0, LITH_BASEMENT, 0.0);
        let mut log = StrataStepLog::default();
        apply_parcel(&mut st, &ParcelAmounts { deposit: 3_000.0, class: 0, ..Default::default() }, 4_000.0, 1.0, &p, 0.0, 2.0, &mut log);
        // The column fell to 2 km with nothing tagged: basement goes to -1 km
        // and takes it from the cover.
        let before = st.total_m();
        let mut l = StrataStepLog::default();
        apply_parcel(&mut st, &ParcelAmounts::default(), 2_000.0, 1.0, &p, 2.0, 4.0, &mut l);
        closes(before, &st, &l);
        assert!(st.basement_m >= 0.0);
        assert!((st.total_m() - 2_000.0).abs() < 1e-9);
        assert!(l.underflow_m3 > 0.0);
    }
}
