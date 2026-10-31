//! SEED-MAGMATISM (lane MG, 2026-10-13): the pf-lab glue of the opt-in law
//! `magmatism = slab-arc-backarc-llsvp-v1` (profile `cycle_magma_v1`;
//! pre-registered in `unpublished-design-record`).
//!
//! - `MagmaEnv::new` re-draws the plume positions at the LLSVP edges and keeps
//!   their drifts; `MagmaEnv::drift` moves them each step.
//! - `magma_step` runs after the tectonic advance: vent markers ride with the
//!   crust (source map), plume / LIP / rift eruptions of the step are recorded,
//!   arcs are built at slab depth (`pf_geodynamics::arcs`), island arcs mature
//!   into continental crust, back-arcs extend (`pf_geodynamics::backarc`), and
//!   the arc heat term (water-loaded m) is returned for both isostasy
//!   evaluations of the step. Arc crust is booked as `arc_created_m3`.
//! - `volcano_features` derives the volcano record (never canon);
//!   `magmatic_potentials` the additive resource families (PREVIEW tier).
//!
//! This file depends only on pf-contracts and pf-geodynamics public items, so
//! lane MG can compile it out of tree; preview.rs passes its fields in.
//! Determinism: ascending cell order everywhere; no hashed maps.

use pf_contracts::domain::PlateRecordV0;
use pf_contracts::params::SeedPath;
use pf_contracts::terrain::{
    AssumptionSourceV0, AssumptionV0, CrustKindCodeV0, VolcanoFeatureV0, VolcanoFeaturesDomainV0,
    VolcanoKindV0, VOLCANO_PREVIEW_USAGE_LIMITS,
};
use pf_geodynamics::arcs::{self, ArcParams};
use pf_geodynamics::backarc::{self, BackarcParams};
use pf_geodynamics::hotspots::{drift_plumes, llsvp_edge_plumes, Plume, PlumeDrift, PlumeSourceParams};
use pf_geodynamics::kinematic::{nearest_cell, AdvanceLedger, CrustState};
use pf_geodynamics::volcanoes::{
    volcano_points, VentMarkers, VENT_CONTINENTAL_ARC, VENT_HOTSPOT, VENT_ISLAND_ARC, VENT_LIP, VENT_RIFT,
};
use pf_grid::GridSpec;

/// Resolved constants of `slab-arc-backarc-llsvp-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MagmaParams {
    pub arc: ArcParams,
    pub backarc: BackarcParams,
    pub source: PlumeSourceParams,
    pub vent_active_window_my: f64,
    pub volcano_spacing_m: f64,
    /// SEED-WORLD-COMPOSE (lane CP): `slab-arc-backarc-llsvp-v2` only.
    pub v2: Option<MagmaV2Params>,
}

/// SEED-WORLD-COMPOSE (lane CP, 2026-10-13): the v2 constants (the arc
/// heat's dynamic, non-displacing application is preview.rs's).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MagmaV2Params {
    pub backarc_episode_my: f64,
    pub backarc_max_stretch_factor: f64,
    pub kimberlite_plume_radius_m: f64,
}

/// Plume-pass radius for the kimberlite marker (m): a continental column
/// within this distance of an active plume centre (project assumption).
pub const PLUME_PASS_RADIUS_M: f64 = 500_000.0;

/// Typical edifice relief above the cell mean (m), per vent kind (PREREG D).
pub const EDIFICE_RELIEF_M: [(VolcanoKindV0, f64); 5] = [
    (VolcanoKindV0::IslandArc, 2500.0),
    (VolcanoKindV0::ContinentalArc, 2500.0),
    (VolcanoKindV0::Hotspot, 3000.0),
    (VolcanoKindV0::Rift, 1500.0),
    (VolcanoKindV0::LargeIgneousProvince, 500.0),
];

/// The run-constant part of the law (in `PreviewEnv`).
#[derive(Debug, Clone, PartialEq)]
pub struct MagmaEnv {
    pub params: MagmaParams,
    /// Plume drifts in id order (empty without plumes).
    pub drifts: Vec<PlumeDrift>,
    pub cents: Vec<[f64; 3]>,
    /// Plume head radius (m) of the plume law (LIP marker).
    pub head_radius_m: f64,
    /// Plume tail lifetime (My) when the reorg law ends tails (volcano points).
    pub tail_lifetime_my: Option<f64>,
}

impl MagmaEnv {
    /// Validate the constants, relocate `plumes` (ascending ids) to the LLSVP
    /// edges and keep their drifts.
    pub fn new(
        params: MagmaParams,
        grid: &GridSpec,
        seed_path: &SeedPath,
        plumes: Option<&mut Vec<Plume>>,
        head_radius_m: f64,
        tail_lifetime_my: Option<f64>,
    ) -> Result<Self, String> {
        params.arc.validate()?;
        params.backarc.validate()?;
        if !(params.vent_active_window_my >= 0.0 && params.volcano_spacing_m >= 0.0) {
            return Err("magmatism: vent window and volcano spacing must be non-negative".into());
        }
        let drifts = match plumes {
            Some(p) => llsvp_edge_plumes(p, &params.source, seed_path, grid.radius_m())?,
            None => Vec::new(),
        };
        let cents = pf_geodynamics::flood_basalt::cell_centroids(grid)?;
        Ok(Self {
            params,
            drifts,
            cents,
            head_radius_m,
            tail_lifetime_my,
        })
    }

    /// Move the plumes to their drifted positions at `t_my`.
    pub fn drift(&self, plumes: &mut [Plume], t_my: f64) {
        drift_plumes(plumes, &self.drifts, t_my);
    }

    /// The plumes active at `t_my` (onset passed, tail alive), drifted to
    /// `t_my`, in id order.
    pub fn active_plumes_at(&self, plumes: &[Plume], t_my: f64) -> Vec<Plume> {
        let mut v: Vec<Plume> = plumes
            .iter()
            .filter(|p| {
                p.onset_my <= t_my && self.tail_lifetime_my.map_or(true, |l| p.onset_my + l > t_my)
            })
            .cloned()
            .collect();
        self.drift(&mut v, t_my);
        v
    }
}

/// Per-step bookkeeping of the law (volumes m^3, areas m^2, lengths m).
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MagmaStepLog {
    pub t_my: f64,
    pub fronts: u64,
    pub island_fronts: u64,
    pub continental_fronts: u64,
    pub front_length_m: f64,
    pub arc_m3: f64,
    /// Revision 1: arc volume not created above the arc crust cap.
    #[serde(default)]
    pub arc_capped_m3: f64,
    pub arc_on_oceanic_m3: f64,
    pub arc_on_continental_m3: f64,
    pub shortened_fronts: u64,
    pub median_arc_distance_m: f64,
    pub converted_cells: u64,
    pub converted_area_m2: f64,
    pub converted_m3: f64,
    pub rollback_fronts: u64,
    pub renewed_cells: u64,
    pub stretched_cells: u64,
    pub stretched_m3: f64,
    pub moved_to_arc_m3: f64,
    /// Revision 1: back-arc fronts limited by the arc sliver's room.
    #[serde(default)]
    pub room_limited_fronts: u64,
    /// Continental columns whose crust is mostly arc crust (juvenile arc
    /// continental crust; the continental ledger's arc-accreted part).
    pub juvenile_arc_continental_area_m2: f64,
    pub juvenile_arc_continental_m3: f64,
    pub active_plumes: u64,
    pub max_heat_term_m: f64,
}

/// Run state of the law (checkpointed inside the conserve state).
#[derive(Debug, Clone, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub struct MagmaRunState {
    pub kind: Vec<u8>,
    pub first_t_my: Vec<f64>,
    pub last_t_my: Vec<f64>,
    pub added_m: Vec<f64>,
    pub origin: Vec<u32>,
    pub heat: Vec<f64>,
    pub backarc_t_my: Vec<f64>,
    pub head_m: Vec<f64>,
    pub plume_pass: Vec<u8>,
    pub logs: Vec<MagmaStepLog>,
    /// SEED-WORLD-COMPOSE (lane CP) v2 back-arc episode markers (empty and
    /// not serialized under v1).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub backarc_first_t_my: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub backarc_start_h_m: Vec<f64>,
}

impl MagmaRunState {
    pub fn new(n: usize) -> Self {
        Self::from_markers(VentMarkers::new(n), Vec::new())
    }

    /// SEED-WORLD-COMPOSE (lane CP): the run state for the law version
    /// (`v2` allocates the back-arc episode markers).
    pub fn new_for(n: usize, v2: bool) -> Self {
        if v2 {
            Self::from_markers(VentMarkers::new(n).with_episodes(), Vec::new())
        } else {
            Self::new(n)
        }
    }

    pub fn markers(&self) -> VentMarkers {
        VentMarkers {
            kind: self.kind.clone(),
            first_t_my: self.first_t_my.clone(),
            last_t_my: self.last_t_my.clone(),
            added_m: self.added_m.clone(),
            origin: self.origin.clone(),
            heat: self.heat.clone(),
            backarc_t_my: self.backarc_t_my.clone(),
            head_m: self.head_m.clone(),
            plume_pass: self.plume_pass.clone(),
            backarc_first_t_my: self.backarc_first_t_my.clone(),
            backarc_start_h_m: self.backarc_start_h_m.clone(),
        }
    }

    fn from_markers(m: VentMarkers, logs: Vec<MagmaStepLog>) -> Self {
        Self {
            kind: m.kind,
            first_t_my: m.first_t_my,
            last_t_my: m.last_t_my,
            added_m: m.added_m,
            origin: m.origin,
            heat: m.heat,
            backarc_t_my: m.backarc_t_my,
            head_m: m.head_m,
            plume_pass: m.plume_pass,
            logs,
            backarc_first_t_my: m.backarc_first_t_my,
            backarc_start_h_m: m.backarc_start_h_m,
        }
    }
}

/// One step of the law on the advanced state `next` (before any consumer of
/// it). `plumes` are the plumes the step's plume pass used (positions as
/// placed); `t_old` the step start. Returns the arc heat term (water-loaded
/// m per cell) for both column terms.
#[allow(clippy::too_many_arguments)]
pub fn magma_step(
    menv: &MagmaEnv,
    grid: &GridSpec,
    areas: &[f64],
    plates: &[PlateRecordV0],
    next: &mut CrustState,
    ledger: &mut AdvanceLedger,
    thickening: &mut [f64],
    source: &[u32],
    st: &mut MagmaRunState,
    plumes: Option<&[Plume]>,
    t_old: f64,
    dt_my: f64,
) -> Result<Vec<f64>, String> {
    magma_step_with_formation_age_v1(menv, grid, areas, plates, next, ledger,
        thickening, source, st, plumes, t_old, dt_my, None)
}

/// Q29: the same native magmatic operation with the run's persistent observer.
#[allow(clippy::too_many_arguments)]
pub fn magma_step_with_formation_age_v1(
    menv: &MagmaEnv,
    grid: &GridSpec,
    areas: &[f64],
    plates: &[PlateRecordV0],
    next: &mut CrustState,
    ledger: &mut AdvanceLedger,
    thickening: &mut [f64],
    source: &[u32],
    st: &mut MagmaRunState,
    plumes: Option<&[Plume]>,
    t_old: f64,
    dt_my: f64,
    formation: Option<&mut pf_geodynamics::formation_age::FormationAgeRecorderV1>,
) -> Result<Vec<f64>, String> {
    let n = grid.cells();
    let cents = &menv.cents;
    let p = &menv.params;
    if next.cells() != n || areas.len() != n || thickening.len() != n || source.len() != n {
        return Err("magma_step: field sizes differ from the grid".into());
    }
    let old = st.markers();
    old.check(n)?;
    let mut m = old.carry(source);
    let t_new = next.t_my;

    // Plume heads arriving this step (their disc = LIP crust) and the
    // plume-pass marker.
    let mut head_cell = vec![false; n];
    let mut active = 0u64;
    if let Some(pl) = plumes {
        for q in pl.iter().filter(|q| q.onset_my < t_new) {
            let c0 = nearest_cell(grid, q.position).map_err(|e| e.to_string())?;
            active += 1;
            if q.onset_my > t_old {
                for (k, _) in pf_geodynamics::orogen::disc_weights(grid, cents, c0, menv.head_radius_m, &|_| true)? {
                    head_cell[k] = true;
                }
            }
            // SEED-WORLD-COMPOSE (lane CP) v2: only an arriving head marks
            // the kimberlite plume pass, within the v2 radius.
            match p.v2 {
                None => {
                    for (k, _) in pf_geodynamics::orogen::disc_weights(grid, cents, c0, PLUME_PASS_RADIUS_M, &|_| true)? {
                        m.plume_pass[k] = 1;
                    }
                }
                Some(v2) => {
                    if q.onset_my > t_old {
                        for (k, _) in pf_geodynamics::orogen::disc_weights(grid, cents, c0, v2.kimberlite_plume_radius_m, &|_| true)? {
                            m.plume_pass[k] = 1;
                        }
                    }
                }
            }
        }
    }
    for &(c, dh, t) in &ledger.hotspot_cells {
        let c = c as usize;
        if head_cell[c] {
            m.record(c, VENT_LIP, t, dh);
            m.head_m[c] += dh;
        } else {
            m.record(c, VENT_HOTSPOT, t, dh);
        }
    }
    for &(c, dh, t) in &ledger.lip_cells {
        m.record(c as usize, VENT_LIP, t, dh);
        m.head_m[c as usize] += dh;
    }
    for &b in &ledger.rift_births {
        m.record(b as usize, VENT_RIFT, t_new, 0.0);
    }

    // Arcs at slab depth, then island-arc maturation.
    let fronts = arcs::find_fronts(grid, cents, areas, next, &p.arc)?;
    let mut dh = vec![0.0f64; n];
    let was_continental: Vec<bool> = next.kind.iter().map(|k| *k == CrustKindCodeV0::Continental).collect();
    let mut arc_rep = arcs::arc_pass(grid, cents, areas, next, &fronts, &p.arc, dt_my, &mut dh)?;
    for c in 0..n {
        if dh[c] > 0.0 {
            let kind = if was_continental[c] { VENT_CONTINENTAL_ARC } else { VENT_ISLAND_ARC };
            m.record(c, kind, t_new, dh[c]);
        }
    }
    pf_geodynamics::volcanoes::supply_capture::arc_control(arc_rep.arc_m3, arc_rep.arc_capped_m3)?;
    arcs::convert_arcs(next, areas, &p.arc, &mut arc_rep);

    // Back-arc extension (volume moved, not created).
    let episode = match p.v2 {
        Some(v2) => {
            if m.backarc_first_t_my.len() != n {
                return Err("magma_step: slab-arc-backarc-llsvp-v2 run state without back-arc episode markers".into());
            }
            Some(backarc::BackarcEpisode {
                episode_my: v2.backarc_episode_my,
                max_stretch_factor: v2.backarc_max_stretch_factor,
                first_t_my: &mut m.backarc_first_t_my,
                start_h_m: &mut m.backarc_start_h_m,
            })
        }
        None => None,
    };
    let ba = backarc::backarc_pass_episode_with_formation_age_v1(
        grid,
        cents,
        areas,
        next,
        plates,
        &fronts,
        &p.backarc,
        p.arc.half_width_m,
        p.arc.crust_cap_m,
        t_new,
        dt_my,
        &mut dh,
        &mut m.backarc_t_my,
        episode,
        formation,
    )?;

    // Arc heat and its elevation term.
    arcs::update_heat(grid, cents, &next.plate, &fronts, &p.arc, dt_my, &mut m.heat)?;
    let term = arcs::thermal_term(&m.heat, &p.arc);

    // Ledger: arc crust is the named magmatic source; the back-arc move is
    // volume neutral; the advance's closure terms are refreshed.
    ledger.arc_created_m3 += arc_rep.arc_m3;
    ledger.arc_placed_m3 += arc_rep.arc_m3;
    for c in 0..n {
        thickening[c] += dh[c];
    }
    ledger.volume_after_m3 = next.crust_volume_m3(areas);
    ledger.max_thickness_m = next.thickness_m.iter().cloned().fold(0.0, f64::max);

    let mut log = MagmaStepLog {
        t_my: t_new,
        fronts: arc_rep.fronts,
        island_fronts: arc_rep.island_fronts,
        continental_fronts: arc_rep.continental_fronts,
        front_length_m: arc_rep.front_length_m,
        arc_m3: arc_rep.arc_m3,
        arc_capped_m3: arc_rep.arc_capped_m3,
        arc_on_oceanic_m3: arc_rep.arc_on_oceanic_m3,
        arc_on_continental_m3: arc_rep.arc_on_continental_m3,
        shortened_fronts: arc_rep.shortened_fronts,
        median_arc_distance_m: arc_rep.median_arc_distance_m,
        converted_cells: arc_rep.converted_cells,
        converted_area_m2: arc_rep.converted_area_m2,
        converted_m3: arc_rep.converted_m3,
        rollback_fronts: ba.rollback_fronts,
        renewed_cells: ba.renewed_cells,
        stretched_cells: ba.stretched_cells,
        stretched_m3: ba.stretched_m3,
        moved_to_arc_m3: ba.moved_to_arc_m3,
        room_limited_fronts: ba.room_limited_fronts,
        active_plumes: active,
        max_heat_term_m: term.iter().cloned().fold(0.0, f64::max),
        ..Default::default()
    };
    for c in 0..n {
        if next.kind[c] == CrustKindCodeV0::Continental && next.arc_added_m[c] >= 0.5 * next.thickness_m[c] {
            log.juvenile_arc_continental_area_m2 += areas[c];
            log.juvenile_arc_continental_m3 += areas[c] * next.thickness_m[c];
        }
    }
    let mut logs = std::mem::take(&mut st.logs);
    logs.push(log);
    *st = MagmaRunState::from_markers(m, logs);
    Ok(term)
}

fn kind_of(code: u8) -> Option<VolcanoKindV0> {
    match code {
        VENT_ISLAND_ARC => Some(VolcanoKindV0::IslandArc),
        VENT_CONTINENTAL_ARC => Some(VolcanoKindV0::ContinentalArc),
        VENT_HOTSPOT => Some(VolcanoKindV0::Hotspot),
        VENT_RIFT => Some(VolcanoKindV0::Rift),
        VENT_LIP => Some(VolcanoKindV0::LargeIgneousProvince),
        _ => None,
    }
}

fn relief_of(k: VolcanoKindV0) -> f64 {
    EDIFICE_RELIEF_M.iter().find(|(x, _)| *x == k).map(|x| x.1).unwrap_or(0.0)
}

fn project_assumption(id: &str, statement: &str, value: f64, units: &str) -> AssumptionV0 {
    AssumptionV0 {
        id: id.into(),
        statement: statement.into(),
        value,
        units: units.into(),
        source: AssumptionSourceV0::ProjectAssumption,
        note: "lane MG 2026-10-13, unpublished-design-record".into(),
    }
}

/// The derived volcano record of a state. `elevation` is the published
/// elevation (sea level at `sea_m`); `plumes` the env's plume list.
pub fn volcano_features(
    menv: &MagmaEnv,
    grid: &GridSpec,
    state: &CrustState,
    elevation: &[f64],
    sea_m: f64,
    st: &MagmaRunState,
    plumes: Option<&[Plume]>,
) -> Result<VolcanoFeaturesDomainV0, String> {
    let t = state.t_my;
    let m = st.markers();
    let live: Vec<(usize, usize)> = match plumes {
        Some(pl) => {
            let mut v = Vec::new();
            for q in menv.active_plumes_at(pl, t) {
                v.push((q.id, nearest_cell(grid, q.position).map_err(|e| e.to_string())?));
            }
            v
        }
        None => Vec::new(),
    };
    let pts = volcano_points(
        grid,
        &menv.cents,
        &m,
        t,
        menv.params.vent_active_window_my,
        menv.params.volcano_spacing_m,
        &live,
    )?;
    let mut features = Vec::with_capacity(pts.len());
    for v in pts {
        let Some(kind) = kind_of(v.kind) else { continue };
        let c = v.cell as usize;
        let (i, j) = grid.ij(c as u64).map_err(|e| e.to_string())?;
        let lat = grid.center_lat_deg(i).map_err(|e| e.to_string())?;
        let lon = grid.center_lon_deg(j).map_err(|e| e.to_string())?;
        let summit = elevation[c] + relief_of(kind);
        let id = match v.plume {
            Some(pid) => format!("hotspot-plume-{pid}"),
            None => format!("{}-{:06}-{:.1}", kind.tag(), v.origin, v.first_t_my),
        };
        features.push(VolcanoFeatureV0 {
            id,
            kind,
            lat_deg: lat,
            lon_deg: lon,
            cell: c as u64,
            onset_my: v.first_t_my,
            last_eruption_my: v.last_t_my,
            active: v.active,
            magmatic_crust_m: v.added_m.max(0.0),
            cell_elevation_m: elevation[c] - sea_m,
            subgrid_summit_m: summit - sea_m,
            emergent: summit > sea_m,
        });
    }
    let rec = VolcanoFeaturesDomainV0 {
        t_my: t.max(0.0),
        usage_limits: VOLCANO_PREVIEW_USAGE_LIMITS.into(),
        edifice_relief_m: EDIFICE_RELIEF_M.to_vec(),
        features,
        assumptions: vec![
            project_assumption("A-VOLCANO-ACTIVE-WINDOW", "a vent is active if it erupted within this time of the snapshot", menv.params.vent_active_window_my, "My"),
            project_assumption("A-VOLCANO-SPACING", "minimum spacing of arc and hotspot volcano points (rift 150 km, LIP 500 km)", menv.params.volcano_spacing_m, "m"),
            project_assumption("A-VOLCANO-ARC-RELIEF", "arc stratovolcano relief above the cell mean (subgrid summit projection)", 2500.0, "m"),
            project_assumption("A-VOLCANO-HOTSPOT-RELIEF", "hotspot shield relief above the cell mean (subgrid summit projection)", 3000.0, "m"),
        ],
    };
    rec.validate().map_err(|e| e.to_string())?;
    Ok(rec)
}

/// One additive resource family (the enum variant is matched by field name
/// in preview.rs).
#[derive(Debug, Clone, PartialEq)]
pub struct MagmaticFamily {
    pub field: &'static str,
    pub basis: Vec<&'static str>,
    pub formula: String,
    pub values: Vec<f32>,
    pub mean_potential: f64,
}

fn sat(x: f64, s: f64) -> f64 {
    (x / s).clamp(0.0, 1.0)
}

/// Exposure window: rises to 1 at `lo` of cumulative erosion, falls to 0 at
/// 2 `hi`.
fn window(eroded_m: f64, lo: f64, hi: f64) -> f64 {
    sat(eroded_m, lo) * (1.0 - sat(eroded_m - hi, hi))
}

pub const PORPHYRY_ARC_SAT_M: f64 = 3000.0;
pub const EPITHERMAL_DECAY_MY: f64 = 50.0;
pub const VMS_RIDGE_AGE_MY: f64 = 20.0;
pub const VMS_BACKARC_DECAY_MY: f64 = 200.0;
pub const LIP_SAT_M: f64 = 2000.0;
pub const KIMBERLITE_AGE_MY: f64 = 2500.0;
pub const SHORTENING_SAT_M: f64 = 10_000.0;

/// The seven additive magmatic / hydrothermal families of a final state.
pub fn magmatic_potentials(
    state: &CrustState,
    elevation: &[f64],
    sea_m: f64,
    areas: &[f64],
    st: &MagmaRunState,
) -> Result<Vec<MagmaticFamily>, String> {
    let n = state.cells();
    let m = st.markers();
    m.check(n)?;
    let t = state.t_my;
    let total: f64 = areas.iter().sum();
    let defs: [(&'static str, Vec<&'static str>, String); 7] = [
        ("potential_porphyry_cu_mo_au", vec!["arc_added", "cumulative_erosion_m", "crust_kind"],
         format!("sat(arc/{PORPHYRY_ARC_SAT_M}) * (0.3 + 0.7 * window(eroded, 1 km, 6 km)) * (continental ? 1 : 0.6)")),
        ("potential_epithermal_au_ag", vec!["vent kind and last eruption", "vent magmatic crust", "cumulative_erosion_m", "elevation_m"],
         format!("(arc vent ? sat(added/1000) : rift vent ? 0.5 : 0) * exp(-(t - last)/{EPITHERMAL_DECAY_MY}) * (1 - sat(eroded/1500)) * (land ? 1 : 0.5)")),
        ("potential_vms", vec!["crust_kind", "crust_age_my", "back-arc extension time"],
         format!("max(oceanic ? 0.6 exp(-age/{VMS_RIDGE_AGE_MY}) : 0, back-arc ? exp(-(t - t_backarc)/{VMS_BACKARC_DECAY_MY}) : 0)")),
        ("potential_magmatic_ni_cu_pge", vec!["plume-head crust", "crust_kind"],
         format!("sat(head/{LIP_SAT_M}) * (continental ? 1 : 0.4)")),
        ("potential_kimberlite_diamond", vec!["crust_kind", "crust_age_my", "shortening_added", "arc_added", "plume pass"],
         format!("continental * sat(age/{KIMBERLITE_AGE_MY}) * (1 - sat(shortening/{SHORTENING_SAT_M})) * (1 - sat(arc/2000)) * plume_pass")),
        ("potential_pegmatite_crystal", vec!["shortening_added", "cumulative_erosion_m"],
         format!("sat(shortening/{SHORTENING_SAT_M}) * (0.3 + 0.7 * window(eroded, 3 km, 15 km))")),
        ("potential_granite_sn_w", vec!["crust_kind", "shortening_added", "cumulative_erosion_m"],
         format!("continental * sat(shortening/{SHORTENING_SAT_M}) * (0.3 + 0.7 * window(eroded, 1 km, 8 km))")),
    ];
    let mut out = Vec::with_capacity(7);
    for (fi, (field, basis, formula)) in defs.into_iter().enumerate() {
        let mut values = Vec::with_capacity(n);
        let mut weighted = 0.0;
        for c in 0..n {
            let cont = state.kind[c] == CrustKindCodeV0::Continental;
            let e = state.eroded_total_m[c].max(0.0);
            let arc = state.arc_added_m[c].max(0.0);
            let sh = state.shortening_added_m[c].max(0.0);
            let v = match fi {
                0 => sat(arc, PORPHYRY_ARC_SAT_M) * (0.3 + 0.7 * window(e, 1000.0, 6000.0)) * if cont { 1.0 } else { 0.6 },
                1 => {
                    let base = match m.kind[c] {
                        VENT_ISLAND_ARC | VENT_CONTINENTAL_ARC => sat(m.added_m[c], 1000.0),
                        VENT_RIFT => 0.5,
                        _ => 0.0,
                    };
                    let rec = if m.last_t_my[c] >= 0.0 { (-(t - m.last_t_my[c]).max(0.0) / EPITHERMAL_DECAY_MY).exp() } else { 0.0 };
                    base * rec * (1.0 - sat(e, 1500.0)) * if elevation[c] >= sea_m { 1.0 } else { 0.5 }
                }
                2 => {
                    let ridge = if cont { 0.0 } else { 0.6 * (-state.age_my(c).max(0.0) / VMS_RIDGE_AGE_MY).exp() };
                    let ba = if m.backarc_t_my[c] >= 0.0 { (-(t - m.backarc_t_my[c]).max(0.0) / VMS_BACKARC_DECAY_MY).exp() } else { 0.0 };
                    ridge.max(ba)
                }
                3 => sat(m.head_m[c], LIP_SAT_M) * if cont { 1.0 } else { 0.4 },
                4 => {
                    if cont && m.plume_pass[c] != 0 {
                        sat(state.age_my(c), KIMBERLITE_AGE_MY) * (1.0 - sat(sh, SHORTENING_SAT_M)) * (1.0 - sat(arc, 2000.0))
                    } else {
                        0.0
                    }
                }
                5 => sat(sh, SHORTENING_SAT_M) * (0.3 + 0.7 * window(e, 3000.0, 15_000.0)),
                _ => {
                    if cont {
                        sat(sh, SHORTENING_SAT_M) * (0.3 + 0.7 * window(e, 1000.0, 8000.0))
                    } else {
                        0.0
                    }
                }
            };
            let v = v.clamp(0.0, 1.0);
            weighted += v * areas[c];
            values.push(v as f32);
        }
        out.push(MagmaticFamily { field, basis, formula, values, mean_potential: weighted / total });
    }
    Ok(out)
}

/// The assumption records the additive families add to the resource domain.
pub fn magmatic_assumptions() -> Vec<AssumptionV0> {
    magmatic_assumptions_for(None)
}

/// SEED-WORLD-COMPOSE (lane CP): the assumption records for the law
/// version (v2: the head-only kimberlite radius).
pub fn magmatic_assumptions_for(v2: Option<&MagmaV2Params>) -> Vec<AssumptionV0> {
    let mut v = magmatic_assumptions_v1();
    if let Some(v2) = v2 {
        for a in v.iter_mut() {
            if a.id == "A-RES-KIMBERLITE-PLUME" {
                a.id = "A-RES-KIMBERLITE-PLUME-V2".into();
                a.statement = "radius around an ARRIVING plume head within which continental columns get the kimberlite plume-pass marker (slab-arc-backarc-llsvp-v2; tails do not mark)".into();
                a.value = v2.kimberlite_plume_radius_m;
                a.note = "lane CP 2026-10-13, unpublished-design-record (lane MG idea 4)".into();
            }
        }
    }
    v
}

fn magmatic_assumptions_v1() -> Vec<AssumptionV0> {
    let pa = |id: &str, s: &str, v: f64, u: &str| {
        let mut a = project_assumption(id, s, v, u);
        a.note = "coarse potential scale; lane MG PREREG.txt E".into();
        a
    };
    vec![
        pa("A-RES-PORPHYRY-ARC", "arc crust at which porphyry potential saturates (porphyries form 1-4 km deep above arcs; Sillitoe 2010, recalled)", PORPHYRY_ARC_SAT_M, "m"),
        pa("A-RES-EPITHERMAL-DECAY", "e-folding time of epithermal potential after the last arc/rift eruption", EPITHERMAL_DECAY_MY, "My"),
        pa("A-RES-VMS-RIDGE-AGE", "e-folding crust age of ridge VMS potential", VMS_RIDGE_AGE_MY, "My"),
        pa("A-RES-VMS-BACKARC", "e-folding time of back-arc VMS potential after extension", VMS_BACKARC_DECAY_MY, "My"),
        pa("A-RES-LIP", "plume-head crust at which Ni-Cu-PGE potential saturates", LIP_SAT_M, "m"),
        pa("A-RES-KIMBERLITE-AGE", "continental crust age at which kimberlite potential saturates (Clifford's rule)", KIMBERLITE_AGE_MY, "My"),
        pa("A-RES-KIMBERLITE-PLUME", "plume-pass radius for the kimberlite marker", PLUME_PASS_RADIUS_M, "m"),
        pa("A-RES-COLLISION-GRANITE", "collisional shortening at which pegmatite and Sn-W potential saturate", SHORTENING_SAT_M, "m"),
    ]
}
