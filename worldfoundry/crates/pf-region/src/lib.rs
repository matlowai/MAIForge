//! `pf-region`: terrain tier T1a, the chart-wide drainage skeleton (ADR 043
//! section 1, node TER-SKELETON), on the lod-8 lattice of ADR 044 section 3
//! (512 m posts over the one `local_oblique_stereographic_v1` chart,
//! +-131.072 km, so 512 x 512 cells and 513 x 513 edge-inclusive posts).
//!
//! LABEL: pre-signoff. The T0 substrate it reads is pre-signoff (SEED-SIGNOFF
//! not passed, R2 not met); its fidelity tier is the weakest upstream one
//! (Preview). Derived, never canon (ADR 043 section 6).
//!
//! What T1a reads (ADR 043 section 1, T1a row): T0 elevation, sea level,
//! reaches with discharge, basins, lake stages, and the T0 subgrid relief
//! (trunk level = cell mean minus subgrid relief). By default it does NOT read
//! climate (ADR 043 OPEN 6); lane CLIMATE-RUNOFF's opt-in key climate=<dir>
//! makes discharge carry a seed-climate-v1 runoff field (climate_runoff.rs).
//!
//! What it preserves: every T0 reach (cell-to-receiver segment) as a
//! connected channel draining to the same T0 outlet basin; every T0 cell's
//! area draining through that cell's own T0 outlet (F3); no channel at or
//! above the river threshold that does not map to a T0 reach (F2).
//!
//! What it invents (declared, keyed, bounded): the channel planform inside a
//! corridor of half width W around each T0 segment (keyed lateral offsets,
//! `vz.keyed-rng/v1`); the headwater extension of headwater cells (declared
//! geometric rule); the bed profile between T0 trunk levels (declared linear
//! rule, carved monotone). Hillslope routing to the channel is a declared
//! deterministic shortest-path rule, not a landform. Nothing here is noise
//! standing in for a process (ADR 026 section 1): the only draws are planform
//! offsets inside the corridor that ADR 043 allows T1a to invent.

#![forbid(unsafe_code)]

pub mod checks;
/// Lane CLIMATE-RUNOFF (2026-10-20): opt-in climate runoff field (key climate=<dir>).
pub mod climate_runoff;
pub mod publish;
pub mod skeleton;
pub mod t0;
pub mod t0_erate;
// T1b regional relief (lane TR, TER-REGION, 2026-10-12).
pub mod t1b;
pub mod t1b_checks;
pub mod t1b_curves;
pub mod t1b_guides;
pub mod t1b_lem;
// T1b 0.8.0 NET-ONE (lane NO, 2026-10-13).
pub mod t1b_net;
// T1b 0.9.0 RIVERNET (lane RV, 2026-10-13).
pub mod t1b_rivernet;
// T1b 0.9.1 COAST-CUT + VOLCANO CONSTRUCT (lane CO, 2026-10-13).
pub mod t1b_coast;
// T1b 0.9.3 defect fixes (lane TD, 2026-10-14).
pub mod t1b_fix;
// T1b 0.9.4 consistent level solve, bank ramp, robust E update (lane TA, 2026-10-14).
pub mod t1b_ls;
// T1b 0.9.4 cell-local E solve (lane TE, 2026-10-14; design review consult 3 Q4), opt-in.
pub mod t1b_cell;
/// Lane TF (2026-10-15): walk-clean fixes on config C (opt-in under params.ls.tb.tc.tf).
pub mod t1b_tf;
/// Lane TG (2026-10-15): shore field, sill grade, box-averaged bilinear, F deadband (opt-in under params.ls.tb.tc.tg).
pub mod t1b_tg;
/// Lane TH (2026-10-15): alongshore exposure smoothing, uncut_mm for a consumer cliff rule (opt-in under params.ls.tb.tc.th).
pub mod t1b_th;
/// Lane VE (2026-10-15): eruption styles and keyed eruption history of the T1b volcano edifices (opt-in ve).
pub mod t1b_ve;
/// Lane VALLEY-FLOOR-FIX: opt-in levee-backswamp floor, colluvial footslope and tributary fans (keys vf_*).
pub mod t1b_vf;
/// Lane HS-PASS (2026-10-21): opt-in hillslope-only transient on the window surface (keys hs / hs_*; design review option (a)).
pub mod t1b_hs;
/// Lane SEAM-CERT (2026-10-21): opt-in pinned-only construction fills, evaluation-lattice routing and a bounded-dependency certificate (key seam_cert=1).
pub mod t1b_seam;
/// Lane E-SOLVE-GENERALITY (2026-10-21): opt-in E-solve generality options (key es=1).
pub mod t1b_es;
/// Lane E-SOLVE-CONVERGE (2026-10-22): convergent E-solve driver (opt-in esc=1 with es=1; LS1 from one fixed reference network every evaluation, dense finite-difference Jacobian, LM steps).
pub mod t1b_esc;
pub mod t1b_esc_local;
pub mod t1b_esc_deps;
/// Lane GLACIAL-T1 (2026-10-21): opt-in T1b glacial overprint (key glacial=<t0.glacial dir>).
pub mod t1b_glacial;
/// Lane T0-IMPRINT-PROD (2026-10-21, ADR 045): one network confinement rule (key conf=) and the per-T0-outlet F3 /
/// whole-cell F6 gates (key outlet_checks=1).
pub mod t1b_outlet;
/// Lane CHANNEL-CENTRELINE (2026-10-21): opt-in T1a-curve centreline and declared decorative meander belt (keys cl=1, cl_mn=1; window-only).
pub mod t1b_cl;
/// Lane MEANDER-Q (2026-10-22): opt-in kinematic meander migration replacing cl_mn's sine belt (key cl_mq=1, needs cl=1; window-only).
pub mod t1b_mq;
/// Lane DAMMED-VALLEY (2026-10-21): opt-in edifice-dammed lakes (key dam_lake=1; window-only).
pub mod t1b_dam;
/// Lane EDIFICE-FOOT (2026-10-21): opt-in window-only edifice foot cap (ef_cap=1) and the edifice measurement dump.
pub mod t1b_ef;
/// Lane T1B-PEDESTAL: measurement dump (env var only) and the opt-in pedestal rule.
pub mod t1b_ped;
/// Lane STEP (2026-10-16): cliff-and-talus hillslopes on hard rock (law threshold-slope-lithology-v1, opt-in under params.ls.tb.tc.ct).
pub mod t1b_ct;
pub mod t1b_litho;
/// Lane VOLCANIC-SPLIT (2026-10-21): lava vs volcaniclastic facies inside volcanic strata layers (opt-in vs=1).
pub mod t1b_volc;
/// Lane LITHO-CONTINUITY (2026-10-21): planet-wide rock-body identity, continuous contacts across T0 cells (opt-in lcont=1).
pub mod t1b_lcont;
pub mod t1b_view;
pub mod t1b_window;
/// Lane PV (ROCK-LOOK-PROVENANCE): exposed-layer provenance ids and records (shared with pf-detail by path).
pub mod strata_prov;

use serde::Serialize;

pub const GENERATOR_NAME: &str = "t1a-drainage-skeleton";
pub const GENERATOR_VERSION: &str = "0.1.0";
pub const RECIPE_SCHEMA: &str = "vz.t1a-skeleton-recipe/v1";
pub const LAYER: &str = "t1a.skeleton";
pub const LABEL: &str = "pre-signoff; T0 substrate R2 not met; fidelity tier Preview; derived, never canon";
pub const PURPOSE_PLANFORM: &str = "t1a.planform.offset";

/// Sentinels in u32 node fields.
pub const NONE: u32 = u32::MAX;
pub const OFF_CHART: u32 = u32::MAX - 1;

/// Lattice: lod-8 posts over the chart.
pub const LOD: u8 = vz_tilekey::LOD_SKELETON;
pub const SPACING_MM: i64 = 512_000; // vz_tilekey::post_spacing_mm(8)
pub const G_MIN: i64 = -256;
pub const G_MAX: i64 = 256;
pub const N: usize = 513;
pub const NODES: usize = N * N;
/// Integer D8 step lengths in mm (diagonal = round(512000 * sqrt 2)).
pub const STEP_ORTHO_MM: u64 = 512_000;
pub const STEP_DIAG_MM: u64 = 724_077;
/// Julian year (same as pf-surface routing: 3.15576e13 s/My).
pub const SECONDS_PER_YEAR: f64 = 31_557_600.0;

pub type R<T> = Result<T, String>;

/// Declared parameters. Integers only, so the canonical JSON is exact.
/// Values marked [default] are project assumptions (ADR 039 C5) until ruled.
#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct SkeletonParamsV1 {
    /// Corridor half width around each T0 segment [default 2048 m: half a
    /// lod-4 walk tile, so a corridor channel through a walk tile centred on
    /// the T0 line cannot leave that tile (B4)].
    pub corridor_half_width_mm: i64,
    /// Spacing of keyed planform control points along a segment [default 4096 m].
    pub planform_control_spacing_mm: i64,
    /// Planform amplitude as permille of the corridor half width [default 750].
    /// 0 is the zero-invention setting: channels follow the straight T0 line.
    pub planform_amplitude_permille: i64,
    /// Curve sampling step before rasterization [default 128 m].
    pub planform_sample_step_mm: i64,
    /// Keyed retries (attempt segment) before the straight fallback [default 8].
    pub planform_max_attempts: u32,
    /// River threshold in milli m3/s (pf-surface RIVER_THRESHOLD_M3_S = 10;
    /// ADR 032; ADR 043 F1 default).
    pub river_threshold_milli_m3_s: i64,
    /// F3 minimum share in permille (ADR 043 F3 default 950).
    pub f3_min_share_permille: i64,
    /// F3 is gated only for cells with at least this permille of their area
    /// inside the chart (slivers are reported, not gated) [default 100].
    pub f3_min_inside_permille: i64,
    /// Headwater extension length as permille of the outgoing T0 segment
    /// [default 400]: short of the cell edge (500 would reach it), so the
    /// heads of two neighbouring cells never meet on a shared edge or corner;
    /// a divide stays between them.
    pub head_extension_permille: i64,
    /// Computation halo in lattice nodes on every side of the chart [default
    /// 160 = 81.92 km, more than one T0 cell], so routing near the chart edge
    /// is not bent along the edge. Halo nodes use the same chart formula
    /// outside its declared domain; they are never published.
    pub halo_nodes: i64,
    pub head_extension_rule: String,
    pub head_level_rule: String,
    pub bed_rule: String,
    pub trunk_level_rule: String,
    pub hillslope_routing_rule: String,
    pub area_rule: String,
    pub runoff_rule: String,
    pub node_cell_rule: String,
    /// Lane CLIMATE-RUNOFF (opt-in key climate=<dir>): discharge carries the climate runoff field (climate_runoff.rs).
    /// None = T0 uniform runoff (params JSON, digests and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub climate_runoff: Option<climate_runoff::ClimateRunoffParamsV1>,
    /// Lane T0-IMPRINT-PROD (opt-in skeleton keys junction=*; ADR 045): valley-axis junctions displaced from the T0 cell
    /// centres, and the per-T0-outlet F3 reading. None = the cell-centre junctions (params JSON, digests and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub imprint: Option<skeleton::ImprintParamsV1>,
}

impl Default for SkeletonParamsV1 {
    fn default() -> Self {
        SkeletonParamsV1 {
            corridor_half_width_mm: 2_048_000,
            planform_control_spacing_mm: 4_096_000,
            planform_amplitude_permille: 750,
            planform_sample_step_mm: 128_000,
            planform_max_attempts: 8,
            river_threshold_milli_m3_s: 10_000,
            f3_min_share_permille: 950,
            f3_min_inside_permille: 100,
            head_extension_permille: 400,
            halo_nodes: 160,
            head_extension_rule: "opposite_fraction_segment_v1: a T0 headwater cell (no river inflow) gets a valley axis from P(c) - f (P(r) - P(c)) to P(c), f = head_extension_permille / 1000".into(),
            head_level_rule: "cell_mean_at_head_v1: head-extension bed rises linearly from the trunk level at P(c) to the T0 cell mean at its head".into(),
            bed_rule: "lerp_trunk_levels_then_downstream_min_v1: bed linear in curve parameter between endpoint levels (sea level at a mouth, lake stage at a lake), then carved so it never rises along flow".into(),
            trunk_level_rule: "t0_mean_minus_subgrid_relief_v1: trunk level = T0 elevation_m - stream-power-subgrid-v1 subgrid relief (lane OG); lake cells use the lake stage".into(),
            hillslope_routing_rule: "own_cell_shortest_path_v1: each land node drains along the integer-length D8 shortest path, inside its own T0 cell, to the valley axes of that cell (its outgoing, incoming and head segments); no diagonal step between two axis nodes; ties to the lowest node index".into(),
            area_rule: "post_cell_area_over_k2_v1: 512 m x 512 m per post (the post-centred cell), divided by the stereographic point scale squared".into(),
            runoff_rule: "t0_uniform_runoff_v1: the T0 terrain artifact's A-RAIN-UNITS value (m/yr) x node area, same as T0 discharge_proxy".into(),
            node_cell_rule: "chart_inverse_floor_v1: node -> (lat, lon) by the chart inverse, T0 cell = floor on the L-grid".into(),
            climate_runoff: None,
            imprint: None,
        }
    }
}

impl SkeletonParamsV1 {
    pub fn threshold_m3_s(&self) -> f64 { self.river_threshold_milli_m3_s as f64 / 1000.0 }
    pub fn half_width_m(&self) -> f64 { self.corridor_half_width_mm as f64 / 1000.0 }
    pub fn amplitude_mm(&self) -> f64 {
        self.corridor_half_width_mm as f64 * self.planform_amplitude_permille as f64 / 1000.0
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct GeneratorRef {
    pub name: String,
    pub version: String,
    /// ADR 032 fidelity tier: the weakest upstream one.
    pub tier: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct LatticeV1 {
    pub tilekey: String,
    pub lod: u8,
    pub spacing_mm: i64,
    pub g_min: i64,
    pub g_max: i64,
    pub posts_per_axis: usize,
    pub order: String,
    pub lod_ladder_digest: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct SeedScopeV1 {
    pub keyed_rng: String,
    pub seed: u64,
    pub ladder: Vec<String>,
    pub purposes: Vec<String>,
}

/// DerivationRecipeV1 shape (ADR 044 section 1) for the `t1a.skeleton` layer.
#[derive(Clone, Debug, Serialize)]
pub struct SkeletonRecipeV1 {
    pub schema: String,
    pub layer: String,
    pub generator: GeneratorRef,
    pub params: SkeletonParamsV1,
    pub params_digest: String,
    pub universe_digest: String,
    /// Sorted (role, digest).
    pub upstream: Vec<(String, String)>,
    pub seed_scope: SeedScopeV1,
    pub chart: vz_tilekey::chart::ChartDeclV1,
    pub lattice: LatticeV1,
    pub key_domain: String,
    pub label: String,
}

impl SkeletonRecipeV1 {
    pub fn digest(&self) -> R<String> {
        pf_contracts::artifact::digest_canonical(self).map_err(|e| format!("recipe digest: {e}"))
    }
}

pub fn lattice_v1() -> LatticeV1 {
    let ladder = vz_tilekey::lod_ladder_v1_canonical_json();
    LatticeV1 {
        tilekey: vz_tilekey::TILEKEY_V2.into(),
        lod: LOD,
        spacing_mm: SPACING_MM,
        g_min: G_MIN,
        g_max: G_MAX,
        posts_per_axis: N,
        order: "node = (gy - g_min) * 513 + (gx - g_min); x east, y north; chart mm = g * 512000; posts are TileKeyV2 lod-8 global post indices".into(),
        lod_ladder_digest: format!("blake3:{}", blake3::hash(ladder.as_bytes()).to_hex()),
    }
}

/// A square lattice of lod-8 posts, `gmin..=gmax` on both axes. The
/// published lattice is `Lat::chart()`; the build runs on the chart plus a halo.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Lat { pub gmin: i64, pub gmax: i64, pub n: usize }

impl Lat {
    pub fn new(gmin: i64, gmax: i64) -> Lat { Lat { gmin, gmax, n: (gmax - gmin + 1) as usize } }
    pub fn chart() -> Lat { Lat::new(G_MIN, G_MAX) }
    pub fn nodes(&self) -> usize { self.n * self.n }
    #[inline]
    pub fn idx(&self, gx: i64, gy: i64) -> Option<usize> {
        if !(self.gmin..=self.gmax).contains(&gx) || !(self.gmin..=self.gmax).contains(&gy) { return None; }
        Some(((gy - self.gmin) as usize) * self.n + (gx - self.gmin) as usize)
    }
    #[inline]
    pub fn g(&self, i: usize) -> (i64, i64) { ((i % self.n) as i64 + self.gmin, (i / self.n) as i64 + self.gmin) }
    #[inline]
    pub fn mm(&self, i: usize) -> (i64, i64) { let (gx, gy) = self.g(i); (gx * SPACING_MM, gy * SPACING_MM) }
    pub fn half_mm(&self) -> f64 { (self.gmax * SPACING_MM) as f64 }
    pub fn on_edge(&self, i: usize) -> bool { let (x, y) = self.g(i); x == self.gmin || x == self.gmax || y == self.gmin || y == self.gmax }
}

#[inline]
pub fn node_index(gx: i64, gy: i64) -> Option<usize> {
    if !(G_MIN..=G_MAX).contains(&gx) || !(G_MIN..=G_MAX).contains(&gy) { return None; }
    Some(((gy - G_MIN) as usize) * N + (gx - G_MIN) as usize)
}

#[inline]
pub fn node_g(n: usize) -> (i64, i64) {
    ((n % N) as i64 + G_MIN, (n / N) as i64 + G_MIN)
}

#[inline]
pub fn node_mm(n: usize) -> (i64, i64) {
    let (gx, gy) = node_g(n);
    (gx * SPACING_MM, gy * SPACING_MM)
}

/// The 8 neighbours in a fixed order (E, NE, N, NW, W, SW, S, SE), with the
/// step length and whether the step is diagonal.
pub const D8: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];

/// Canonical detail-domain inputs for the LS owner (opt-in API).
pub mod t1b_surface_domain;
pub mod t1b_surface_completion;
pub mod t1b_surface_stages;
pub mod t1b_surface_channels;
pub mod t1b_surface_area;

pub mod t1b_escarp;

// Engine-neutral, planet-addressed native terrain construction.
pub mod planet_datum;
pub mod planet_relief;
pub mod planet_cell;
pub mod planet_catchment;
pub mod planet_hillslope;
pub mod planet_transport;
pub mod planet_initial;
