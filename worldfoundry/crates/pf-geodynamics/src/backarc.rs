//! SEED-MAGMATISM (lane MG, 2026-10-13): the back-arc part of
//! `slab-arc-backarc-llsvp-v1` (`unpublished-design-record`).
//!
//! Where an old (dense) slab rolls back and the overriding plate does not
//! advance toward the trench in the mantle frame, the overriding plate is
//! stretched behind the arc (Mariana / Japan Sea style). There is no new
//! plate here (scripted kinematics stay lane PR/BC's); the extension is a
//! column operation on a band of the overriding plate centred
//! `axis_offset_m` behind the arc:
//! - oceanic band cells: back-arc spreading, the crust is renewed: its age is
//!   set to |x| / (e / 2) (x = distance from the back-arc axis point), only
//!   ever younger; thickness unchanged (volume neutral);
//! - continental band cells: stretched by W / (W + e dt) per step (W = the
//!   band width), never below `min_crust_m`; the removed volume moves to the
//!   front's arc cells (the arc sliver) in proportion to their room below the
//!   arc crust cap, and a front stretches only as far as that room allows
//!   (revision 1: volume conserved exactly; an episode ends when its arc
//!   sliver is full). Fronts in ascending order, each sees the previous
//!   fronts' changes.
//! Band cells closer to the trench than the arc point plus the arc half-width
//! are left to the arc. Each cell is extended at most once per step.
//!
//! Determinism: fronts in their (ascending cell) order, discs ascending.

use crate::arcs::{along, cell_of, Front};
use crate::kinematic::CrustState;
use crate::rotation::great_circle_rad;
use pf_contracts::domain::{PlateId, PlateRecordV0};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;

/// Back-arc constants of `slab-arc-backarc-llsvp-v1` (laws.rs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BackarcParams {
    pub slab_age_min_my: f64,
    pub advance_max_m_per_my: f64,
    pub rate_m_per_my: f64,
    pub axis_offset_m: f64,
    pub half_width_m: f64,
    pub min_crust_m: f64,
}

impl BackarcParams {
    pub fn validate(&self) -> Result<(), String> {
        if !(self.slab_age_min_my >= 0.0
            && self.advance_max_m_per_my.is_finite()
            && self.rate_m_per_my > 0.0
            && self.axis_offset_m >= 0.0
            && self.half_width_m > 0.0
            && self.min_crust_m > 0.0)
            || ![self.slab_age_min_my, self.rate_m_per_my, self.axis_offset_m, self.half_width_m, self.min_crust_m]
                .iter()
                .all(|v| v.is_finite())
        {
            return Err(format!("back-arc parameters out of range: {self:?}"));
        }
        Ok(())
    }
}

/// What one back-arc pass did.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BackarcReport {
    /// Fronts in rollback (old slab, overriding plate not advancing).
    pub rollback_fronts: u64,
    /// Oceanic band cells whose age was renewed, continental band cells
    /// stretched, and the volume moved from stretched cells to arc cells.
    pub renewed_cells: u64,
    pub stretched_cells: u64,
    pub stretched_m3: f64,
    pub moved_to_arc_m3: f64,
    /// Revision 1: fronts whose stretching was limited by arc room.
    pub room_limited_fronts: u64,
}

/// Overriding-plate absolute velocity (m/My) toward the trench at a front's
/// axis point (positive = advancing toward the trench). None if the plate
/// has no pole.
pub fn advance_toward_trench_m_per_my(
    f: &Front,
    plates: &[PlateRecordV0],
    radius_m: f64,
) -> Result<Option<f64>, String> {
    let Some(rec) = plates.iter().find(|r| r.plate_id == f.over_plate) else {
        return Ok(None);
    };
    let Some(pole) = rec.euler_pole.as_ref() else {
        return Ok(None);
    };
    let w = crate::boundaries::angular_velocity(pole).map_err(|e| e.to_string())?;
    let p = f.axis;
    let v = [
        (w[1] * p[2] - w[2] * p[1]) * radius_m,
        (w[2] * p[0] - w[0] * p[2]) * radius_m,
        (w[0] * p[1] - w[1] * p[0]) * radius_m,
    ];
    // `dir` points from the trench into the overriding plate.
    Ok(Some(-(v[0] * f.dir[0] + v[1] * f.dir[1] + v[2] * f.dir[2])))
}

/// Run the back-arc pass on `state` at the step end `t_new`.
/// `arc_half_width_m` is the arc law's half width (band exclusion and the
/// disc that receives the moved volume); `arc_cap_m` the arc crust cap. `dh` receives the thickness
/// changes (m); `backarc_t` is set to `t_new` on every extended cell.
#[allow(clippy::too_many_arguments)]
pub fn backarc_pass(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    state: &mut CrustState,
    plates: &[PlateRecordV0],
    fronts: &[Front],
    p: &BackarcParams,
    arc_half_width_m: f64,
    arc_cap_m: f64,
    t_new: f64,
    dt_my: f64,
    dh: &mut [f64],
    backarc_t: &mut [f64],
) -> Result<BackarcReport, String> {
    backarc_pass_episode(grid, cents, areas, state, plates, fronts, p, arc_half_width_m, arc_cap_m, t_new, dt_my, dh, backarc_t, None)
}

/// SEED-WORLD-COMPOSE (lane CP, 2026-10-13) `slab-arc-backarc-llsvp-v2`:
/// one back-arc episode per cell.
#[derive(Debug)]
pub struct BackarcEpisode<'a> {
    /// A cell is extended only within this time of its first extension.
    pub episode_my: f64,
    /// Continental cells never below 1/this of their episode-start thickness.
    pub max_stretch_factor: f64,
    pub first_t_my: &'a mut [f64],
    pub start_h_m: &'a mut [f64],
}

/// `backarc_pass` with the v2 episode limit (`None` = `backarc_pass` bit
/// for bit). With an episode, a cell already past `episode_my` since its
/// first extension is neither renewed nor stretched; a cell's first
/// extension records the time and its thickness; continental stretching
/// stops at max(min crust, start / max_stretch_factor).
#[allow(clippy::too_many_arguments)]
pub fn backarc_pass_episode(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    state: &mut CrustState,
    plates: &[PlateRecordV0],
    fronts: &[Front],
    p: &BackarcParams,
    arc_half_width_m: f64,
    arc_cap_m: f64,
    t_new: f64,
    dt_my: f64,
    dh: &mut [f64],
    backarc_t: &mut [f64],
    mut episode: Option<BackarcEpisode<'_>>,
) -> Result<BackarcReport, String> {
    p.validate()?;
    if let Some(e) = episode.as_ref() {
        if !(e.episode_my > 0.0 && e.max_stretch_factor >= 1.0) || e.first_t_my.len() != grid.cells() || e.start_h_m.len() != grid.cells() {
            return Err(format!("back-arc episode limit invalid: {} My, factor {}", e.episode_my, e.max_stretch_factor));
        }
    }
    // v2: cells whose episode is over are not extended again.
    let spent = |k: usize, e: &Option<BackarcEpisode<'_>>| -> bool {
        match e {
            Some(e) => e.first_t_my[k] >= 0.0 && t_new - e.first_t_my[k] > e.episode_my + 1e-9,
            None => false,
        }
    };
    let n = grid.cells();
    if dh.len() != n || backarc_t.len() != n || state.cells() != n {
        return Err("backarc_pass: field sizes differ from the grid".into());
    }
    let r = grid.radius_m();
    let mut rep = BackarcReport::default();
    // Per cell: youngest renewed birth (oceanic) or the front that stretches
    // it first (continental), so each cell is extended at most once.
    let mut renew: Vec<Option<f64>> = vec![None; n];
    let mut stretch_front: Vec<Option<usize>> = vec![None; n];
    let mut rollback: Vec<usize> = Vec::new();
    for (fi, f) in fronts.iter().enumerate() {
        if f.age_my < p.slab_age_min_my {
            continue;
        }
        let Some(u) = advance_toward_trench_m_per_my(f, plates, r)? else {
            continue;
        };
        if u > p.advance_max_m_per_my {
            continue;
        }
        rep.rollback_fronts += 1;
        rollback.push(fi);
        let d_axis = f.arc_distance_m + p.axis_offset_m;
        let q = along(f.axis, f.dir, d_axis / r);
        let c0 = cell_of(grid, q)?;
        let over: PlateId = f.over_plate;
        let plate = &state.plate;
        let w = crate::orogen::disc_weights(grid, cents, c0, p.half_width_m, &|k| plate[k] == over)?;
        let near = f.arc_distance_m + arc_half_width_m;
        for (k, _) in w {
            if state.plate[k] != over {
                continue;
            }
            let from_trench = great_circle_rad(f.axis, cents[k]).map_err(|e| e.to_string())? * r;
            if from_trench < near {
                continue;
            }
            if spent(k, &episode) {
                continue;
            }
            match state.kind[k] {
                CrustKindCodeV0::Oceanic => {
                    let x = great_circle_rad(q, cents[k]).map_err(|e| e.to_string())? * r;
                    // Never exactly a birth this step (age >= 0.5 My).
                    let age = (x / (0.5 * p.rate_m_per_my)).max(0.5);
                    let b = t_new - age;
                    renew[k] = Some(match renew[k] {
                        Some(old) => old.max(b),
                        None => b,
                    });
                }
                CrustKindCodeV0::Continental => {
                    if stretch_front[k].is_none() {
                        stretch_front[k] = Some(fi);
                    }
                }
            }
        }
    }
    // Per-cell stretching floor (v1: the minimum crust).
    let floor_of = |k: usize, e: &Option<BackarcEpisode<'_>>| -> f64 {
        match e {
            Some(e) => p.min_crust_m.max(e.start_h_m[k] / e.max_stretch_factor),
            None => p.min_crust_m,
        }
    };
    let width = 2.0 * p.half_width_m;
    let factor = width / (width + p.rate_m_per_my * dt_my);
    let mut by_front: Vec<Vec<usize>> = vec![Vec::new(); fronts.len()];
    for k in 0..n {
        if let Some(b) = renew[k] {
            if b > state.birth_t_my[k] {
                state.birth_t_my[k] = b;
                rep.renewed_cells += 1;
                backarc_t[k] = t_new;
                if let Some(e) = episode.as_mut() {
                    if e.first_t_my[k] < 0.0 {
                        e.first_t_my[k] = t_new;
                        e.start_h_m[k] = state.thickness_m[k];
                    }
                }
            }
        }
        if let Some(fi) = stretch_front[k] {
            if let Some(e) = episode.as_mut() {
                if e.first_t_my[k] < 0.0 {
                    e.first_t_my[k] = t_new;
                    e.start_h_m[k] = state.thickness_m[k];
                }
            }
            if state.thickness_m[k] > floor_of(k, &episode) {
                by_front[fi].push(k);
            }
        }
    }
    // Revision 1: each front moves its stretched volume to its arc cells in
    // proportion to their room below the cap, stretching only as far as the
    // room allows (volume conserved exactly).
    for (fi, cells) in by_front.iter().enumerate() {
        if cells.is_empty() {
            continue;
        }
        let f = &fronts[fi];
        let over = f.over_plate;
        let plate = &state.plate;
        let w = crate::orogen::disc_weights(grid, cents, f.arc_cell as usize, arc_half_width_m, &|k| {
            plate[k] == over
        })?;
        let room: Vec<(usize, f64)> = w
            .iter()
            .filter(|(k, _)| state.plate[*k] == over && stretch_front[*k].is_none())
            .map(|(k, _)| (*k, (arc_cap_m - state.thickness_m[*k]).max(0.0) * areas[*k]))
            .collect();
        let capacity: f64 = room.iter().map(|x| x.1).sum();
        let request: f64 = cells
            .iter()
            .map(|&k| (state.thickness_m[k] - (state.thickness_m[k] * factor).max(floor_of(k, &episode))) * areas[k])
            .sum();
        if !(request > 0.0) || !(capacity > 0.0) {
            if request > 0.0 {
                rep.room_limited_fronts += 1;
            }
            continue;
        }
        let s = (capacity / request).min(1.0);
        if s < 1.0 {
            rep.room_limited_fronts += 1;
        }
        let mut moved = 0.0;
        for &k in cells {
            let h = state.thickness_m[k];
            let d = s * (h - (h * factor).max(floor_of(k, &episode)));
            state.thickness_m[k] = h - d;
            dh[k] -= d;
            moved += d * areas[k];
            rep.stretched_cells += 1;
            backarc_t[k] = t_new;
        }
        rep.stretched_m3 += moved;
        for (k, r) in room {
            if r > 0.0 {
                let add = moved * r / capacity / areas[k];
                state.thickness_m[k] += add;
                dh[k] += add;
                rep.moved_to_arc_m3 += add * areas[k];
            }
        }
    }
    Ok(rep)
}

/// `backarc_pass_episode` with opt-in native thermal-reset provenance. The
/// event time is the operation interval [t_new-dt, t_new]; the distance-derived
/// effective thermal birth remains a separate value, never a material birth.
/// Call after the recorded carry step at the same state time. No recorder means
/// no snapshots or changes to the legacy path. Archive capacity is checked
/// before mutation; successful native operations alone append observations.
#[allow(clippy::too_many_arguments)]
pub fn backarc_pass_episode_with_formation_age_v1(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    state: &mut CrustState,
    plates: &[PlateRecordV0],
    fronts: &[Front],
    p: &BackarcParams,
    arc_half_width_m: f64,
    arc_cap_m: f64,
    t_new: f64,
    dt_my: f64,
    dh: &mut [f64],
    backarc_t: &mut [f64],
    episode: Option<BackarcEpisode<'_>>,
    recorder: Option<&mut crate::formation_age::FormationAgeRecorderV1>,
) -> Result<BackarcReport, String> {
    let Some(recorder) = recorder else {
        return backarc_pass_episode(grid, cents, areas, state, plates, fronts, p,
            arc_half_width_m, arc_cap_m, t_new, dt_my, dh, backarc_t, episode);
    };
    recorder.preflight(grid, t_new)?;
    let interval = crate::formation_age::TimeInterval::new(t_new - dt_my, t_new)?;
    let n = grid.cells();
    if state.t_my != t_new || state.birth_t_my.len() != n || state.thickness_m.len() != n
        || areas.len() != n || cents.len() != n
        || areas.iter().any(|a| !a.is_finite() || *a <= 0.0)
        || state.thickness_m.iter().any(|h| !h.is_finite() || *h < 0.0)
        || state.birth_t_my.iter().any(|t| !t.is_finite())
    {
        return Err("formation-age: invalid back-arc time, volume or grid fields".into());
    }
    let previous_birth = state.birth_t_my.clone();
    let report = backarc_pass_episode(grid, cents, areas, state, plates, fronts, p,
        arc_half_width_m, arc_cap_m, t_new, dt_my, dh, backarc_t, episode)?;
    recorder.thermal_resets(state, &previous_birth, areas, interval);
    Ok(report)
}
