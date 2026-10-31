//! SEED-FLOODBASALT (lane FB, 2026-10-13): `plume-head-lip-swell-v1`,
//! large igneous provinces (flood-basalt plateaus) from plume heads, plus
//! the plume's dynamic-topography swell.
//!
//! Layered on `fixed-plume-chains-v1` (hotspots.rs: fixed mantle plumes
//! with seeded onsets, a head at onset and an age-progressive tail). This
//! law owns the head's magma:
//!
//! 1. LIP crust (real crust). A plume whose onset is inside the run
//!    (onset >= 0) erupts a total volume V = V_extrusive (1 + r_intrusive)
//!    over [onset, onset + duration]; each step takes the share of that
//!    window it overlaps. The material erupted at the overlap midpoint
//!    t_mid sits, at step end, where the plate now under the plume has
//!    carried it, rotate(plume, pole, t_new - t_mid) (the v1 tail's rule),
//!    spread over a flat-topped disc of radius R on ANY plate: weight 1
//!    inside (1 - taper) R, a linear taper to 0 at R, area-normalised so
//!    the volume placed is exact. It is added to the column's crust
//!    thickness, so it is carried, eroded and isostatically compensated
//!    like any crust (oceanic plateaus on oceanic columns). The v1 head
//!    cone is superseded (the caller runs `plume_pass` with head volume 0).
//! 2. Swell (elevation only, no crust). Every active plume (t >= onset)
//!    lifts the columns around its FIXED mantle position by the
//!    water-loaded column term
//!    s = A_tail exp(-(r/R_tail)^2)
//!      + A_head exp(-(t - onset)/tau_head) exp(-(r/R_head)^2),
//!    r the great-circle distance, cut at 3 R; plumes add. Crust carried
//!    off the plume loses its swell.
//!
//! Constants and their recalled sources: `unpublished-design-record
//! flood-basalt/PREREG.txt` (fixed before any run).

use crate::hotspots::Plume;
use pf_contracts::domain::{EulerPoleV0, PlateId};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;
use std::collections::BTreeMap;

/// Parameters of the LIP part of `plume-head-lip-swell-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LipParams {
    pub extrusive_volume_m3: f64,
    pub intrusive_ratio: f64,
    pub radius_m: f64,
    pub taper_fraction: f64,
    pub eruption_duration_my: f64,
}

impl LipParams {
    /// Extrusive plus intrusive / underplated volume of one LIP (m^3).
    pub fn total_volume_m3(&self) -> f64 {
        self.extrusive_volume_m3 * (1.0 + self.intrusive_ratio)
    }

    pub fn validate(&self) -> Result<(), String> {
        let ok = self.extrusive_volume_m3.is_finite()
            && self.extrusive_volume_m3 >= 0.0
            && self.intrusive_ratio.is_finite()
            && self.intrusive_ratio >= 0.0
            && self.radius_m.is_finite()
            && self.radius_m > 0.0
            && self.taper_fraction.is_finite()
            && self.taper_fraction > 0.0
            && self.taper_fraction <= 1.0
            && self.eruption_duration_my.is_finite()
            && self.eruption_duration_my > 0.0;
        if ok {
            Ok(())
        } else {
            Err(format!("plume-head-lip-swell-v1: invalid LIP parameters {self:?}"))
        }
    }
}

/// Parameters of the swell part of `plume-head-lip-swell-v1` (amplitudes
/// are water-loaded column heights, m).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SwellParams {
    pub tail_amplitude_m: f64,
    pub tail_radius_m: f64,
    pub head_amplitude_m: f64,
    pub head_radius_m: f64,
    pub head_decay_my: f64,
}

impl SwellParams {
    pub fn validate(&self) -> Result<(), String> {
        let ok = [self.tail_amplitude_m, self.head_amplitude_m]
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0)
            && [self.tail_radius_m, self.head_radius_m, self.head_decay_my]
                .iter()
                .all(|v| v.is_finite() && *v > 0.0);
        if ok {
            Ok(())
        } else {
            Err(format!("plume-head-lip-swell-v1: invalid swell parameters {self:?}"))
        }
    }

    /// Peak tail and head swell (water-loaded m) of a plume at time t
    /// (0, 0 before onset).
    pub fn amplitudes(&self, onset_my: f64, t_my: f64) -> (f64, f64) {
        if t_my < onset_my {
            return (0.0, 0.0);
        }
        (
            self.tail_amplitude_m,
            self.head_amplitude_m * (-(t_my - onset_my) / self.head_decay_my).exp(),
        )
    }
}

/// Both parts of `plume-head-lip-swell-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FloodBasaltParams {
    pub lip: LipParams,
    pub swell: SwellParams,
}

impl FloodBasaltParams {
    pub fn validate(&self) -> Result<(), String> {
        self.lip.validate()?;
        self.swell.validate()
    }
}

/// What one LIP pass added. `cells` ascending: (cell, added thickness m,
/// eruption midpoint time My).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LipReport {
    /// Total LIP crust placed this step (m^3) and the part placed on
    /// columns continental at eruption.
    pub total_m3: f64,
    pub continental_m3: f64,
    /// Plumes erupting this step.
    pub erupting: u64,
    pub cells: Vec<(u32, f64, f64)>,
}

/// Flat-topped disc weights around `centre`: 1 inside (1 - taper) R, then
/// linear to 0 at R (any plate; `disc_weights`' neighbour walk and
/// distances, so the footprint is the same set of cells as a v1 disc).
pub fn plateau_weights(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    centre: usize,
    radius_m: f64,
    taper_fraction: f64,
) -> Result<Vec<(usize, f64)>, String> {
    let tri = crate::orogen::disc_weights(grid, cents, centre, radius_m, &|_| true)?;
    // disc weight w = 1 - d/R, so d/R = 1 - w; plateau weight = min(1, w / taper).
    Ok(tri
        .into_iter()
        .map(|(c, w)| (c, (w / taper_fraction).min(1.0)))
        .filter(|(_, w)| *w > 0.0)
        .collect())
}

/// Erupt this step's share of every LIP whose window overlaps (t_old,
/// t_new] into `thickness_m` (the post-carry crust of the step end; `plate`
/// and `kind` its plate and kind fields). Deterministic: plumes in id
/// order, cells ascending.
#[allow(clippy::too_many_arguments)]
pub fn lip_pass(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    plumes: &[Plume],
    lp: &LipParams,
    t_old: f64,
    t_new: f64,
    plate: &[PlateId],
    kind: &[CrustKindCodeV0],
    poles: &BTreeMap<PlateId, EulerPoleV0>,
    thickness_m: &mut [f64],
) -> Result<LipReport, String> {
    use crate::kinematic::nearest_cell;
    let n = grid.cells();
    if cents.len() != n || areas.len() != n || plate.len() != n || kind.len() != n || thickness_m.len() != n {
        return Err("lip_pass: field sizes differ from the grid".into());
    }
    lp.validate()?;
    let total = lp.total_volume_m3();
    let mut rep = LipReport::default();
    let mut added: BTreeMap<usize, (f64, f64)> = BTreeMap::new();
    for pl in plumes {
        if !(pl.onset_my >= 0.0) {
            continue;
        }
        let lo = pl.onset_my.max(t_old);
        let hi = (pl.onset_my + lp.eruption_duration_my).min(t_new);
        if !(hi > lo) {
            continue;
        }
        let v = total * (hi - lo) / lp.eruption_duration_my;
        if !(v > 0.0) {
            continue;
        }
        let t_mid = 0.5 * (lo + hi);
        let c0 = nearest_cell(grid, pl.position).map_err(|e| e.to_string())?;
        let pid = plate[c0];
        let pole = poles
            .get(&pid)
            .ok_or_else(|| format!("LIP {}: no pole for plate {}", pl.id, pid.0))?;
        let x = crate::rotation::rotate_by_euler(pl.position, pole, t_new - t_mid)
            .map_err(|e| e.to_string())?;
        let c = nearest_cell(grid, x).map_err(|e| e.to_string())?;
        let w = plateau_weights(grid, cents, c, lp.radius_m, lp.taper_fraction)?;
        let norm: f64 = w.iter().map(|(k, wk)| wk * areas[*k]).sum();
        if !(norm > 0.0) {
            continue;
        }
        rep.erupting += 1;
        if crate::volcanoes::supply_capture::enabled() {
            crate::volcanoes::supply_capture::record(vec![crate::volcanoes::supply_capture::Sample::placed(
                "lip_head", format!("plume/{}/lip", pl.id), c as u32,
                t_old, t_new, v, &w, areas,
            )?])?;
        }
        for (k, wk) in w {
            let dh = v * wk / norm;
            thickness_m[k] += dh;
            let dv = dh * areas[k];
            rep.total_m3 += dv;
            if kind[k] == CrustKindCodeV0::Continental {
                rep.continental_m3 += dv;
            }
            let e = added.entry(k).or_insert((0.0, t_mid));
            e.0 += dh;
            if t_mid > e.1 {
                e.1 = t_mid;
            }
        }
    }
    rep.cells = added.into_iter().map(|(c, (dh, t))| (c as u32, dh, t)).collect();
    Ok(rep)
}

/// Cell-centre unit vectors in the grid's flat order (row-major).
pub fn cell_centroids(grid: &GridSpec) -> Result<Vec<[f64; 3]>, String> {
    let mut out = Vec::with_capacity(grid.cells());
    for i in 0..grid.nlat() {
        let lat = grid.center_lat_deg(i).map_err(|e| e.to_string())?;
        for j in 0..grid.nlon() {
            let lon = grid.center_lon_deg(j).map_err(|e| e.to_string())?;
            out.push(crate::rotation::lat_lon_to_unit(lat, lon).map_err(|e| e.to_string())?);
        }
    }
    Ok(out)
}

/// The swell column term (water-loaded m) of every active plume at time
/// `t_my`, per cell (`cents` the cell centroids, `radius_m` the planet
/// radius).
pub fn swell_term(
    cents: &[[f64; 3]],
    radius_m: f64,
    plumes: &[Plume],
    sp: &SwellParams,
    t_my: f64,
) -> Vec<f64> {
    let mut out = vec![0.0f64; cents.len()];
    for pl in plumes {
        let (a_tail, a_head) = sp.amplitudes(pl.onset_my, t_my);
        if a_tail == 0.0 && a_head == 0.0 {
            continue;
        }
        let r_cut = 3.0 * sp.tail_radius_m.max(sp.head_radius_m);
        let cos_cut = (r_cut / radius_m).min(std::f64::consts::PI).cos();
        let p = pl.position;
        for (k, c) in cents.iter().enumerate() {
            let dot = p[0] * c[0] + p[1] * c[1] + p[2] * c[2];
            if dot < cos_cut {
                continue;
            }
            let r = dot.clamp(-1.0, 1.0).acos() * radius_m;
            let g_tail = (-(r / sp.tail_radius_m).powi(2)).exp();
            let g_head = (-(r / sp.head_radius_m).powi(2)).exp();
            out[k] += a_tail * g_tail + a_head * g_head;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kinematic::nearest_cell;

    fn centroids(grid: &GridSpec) -> Result<Vec<[f64; 3]>, String> {
        cell_centroids(grid)
    }

    fn cell_areas(grid: &GridSpec) -> Result<Vec<f64>, String> {
        let mut out = Vec::with_capacity(grid.cells());
        for i in 0..grid.nlat() {
            let a = grid.cell_area_m2(i).map_err(|e| e.to_string())?;
            out.extend(std::iter::repeat(a).take(grid.nlon()));
        }
        Ok(out)
    }

    fn lp() -> LipParams {
        LipParams {
            extrusive_volume_m3: 1.5e15,
            intrusive_ratio: 4.0,
            radius_m: 800_000.0,
            taper_fraction: 0.3,
            eruption_duration_my: 2.0,
        }
    }

    fn sp() -> SwellParams {
        SwellParams {
            tail_amplitude_m: 800.0,
            tail_radius_m: 500_000.0,
            head_amplitude_m: 700.0,
            head_radius_m: 800_000.0,
            head_decay_my: 10.0,
        }
    }

    fn plume(onset: f64, lat: f64, lon: f64) -> Plume {
        Plume {
            id: 0,
            position: crate::rotation::lat_lon_to_unit(lat, lon).unwrap(),
            onset_my: onset,
            split_azimuth_rad: 0.0,
        }
    }

    #[test]
    fn lip_volume_is_exact_flat_topped_and_split_over_the_window() {
        let grid = GridSpec::new(90, 180, 6_371_000.0).unwrap();
        let cents = centroids(&grid).unwrap();
        let areas = cell_areas(&grid).unwrap();
        let n = grid.cells();
        let plate = vec![PlateId(1); n];
        // Half the globe continental (northern rows), to check the split.
        let kind: Vec<CrustKindCodeV0> = (0..n)
            .map(|c| if c < n / 2 { CrustKindCodeV0::Continental } else { CrustKindCodeV0::Oceanic })
            .collect();
        let mut poles = BTreeMap::new();
        poles.insert(
            PlateId(1),
            EulerPoleV0 { lat_deg: 90.0, lon_deg: 0.0, rate_deg_per_my: 0.0 },
        );
        let pl = [plume(3.0, 0.0, 10.0)];
        let p = lp();
        let mut th = vec![35_000.0; n];
        // Steps (2,4] and (4,6]: the window [3,5] splits 1/2 + 1/2.
        let a = lip_pass(&grid, &cents, &areas, &pl, &p, 2.0, 4.0, &plate, &kind, &poles, &mut th).unwrap();
        let b = lip_pass(&grid, &cents, &areas, &pl, &p, 4.0, 6.0, &plate, &kind, &poles, &mut th).unwrap();
        let c = lip_pass(&grid, &cents, &areas, &pl, &p, 6.0, 8.0, &plate, &kind, &poles, &mut th).unwrap();
        let v = p.total_volume_m3();
        assert!((a.total_m3 - 0.5 * v).abs() <= 1e-9 * v, "{} {}", a.total_m3, v);
        assert!((b.total_m3 - 0.5 * v).abs() <= 1e-9 * v);
        assert_eq!(c.total_m3, 0.0);
        assert_eq!(c.erupting, 0);
        let placed: f64 = th.iter().zip(&areas).map(|(h, a)| (h - 35_000.0) * a).sum();
        assert!((placed - v).abs() <= 1e-9 * v);
        // Straddling the equator: roughly half on continental columns.
        let frac = (a.continental_m3 + b.continental_m3) / v;
        assert!(frac > 0.3 && frac < 0.7, "{frac}");
        // Flat top: the inner 70 % of the radius all share the peak.
        let c0 = nearest_cell(&grid, pl[0].position).unwrap();
        let peak = th[c0] - 35_000.0;
        assert!(peak > 3000.0 && peak < 8000.0, "peak {peak}");
        for k in 0..n {
            let d = crate::rotation::great_circle_rad(cents[c0], cents[k]).unwrap() * grid.radius_m();
            let dh = th[k] - 35_000.0;
            if d < 0.6 * p.radius_m {
                assert!((dh - peak).abs() <= 1e-6 * peak, "inner cell {k} {dh} vs {peak}");
            }
            if d > p.radius_m {
                assert_eq!(dh, 0.0);
            }
        }
    }

    #[test]
    fn plumes_active_before_the_run_erupt_no_lip_and_moving_plate_shifts_it() {
        let grid = GridSpec::new(90, 180, 6_371_000.0).unwrap();
        let cents = centroids(&grid).unwrap();
        let areas = cell_areas(&grid).unwrap();
        let n = grid.cells();
        let plate = vec![PlateId(1); n];
        let kind = vec![CrustKindCodeV0::Oceanic; n];
        let mut poles = BTreeMap::new();
        // Rotation about the z axis: eastward motion at the equator.
        poles.insert(
            PlateId(1),
            EulerPoleV0 { lat_deg: 90.0, lon_deg: 0.0, rate_deg_per_my: 1.0 },
        );
        let mut th = vec![7_000.0; n];
        let early = [plume(-1.0, 0.0, 0.0)];
        let r = lip_pass(&grid, &cents, &areas, &early, &lp(), 0.0, 2.0, &plate, &kind, &poles, &mut th).unwrap();
        assert_eq!(r.erupting, 0);
        assert!(th.iter().all(|&h| h == 7_000.0));
        // Onset at 0, window [0,2], midpoint 1 My: carried 1 deg east by 2 My.
        let p0 = [plume(0.0, 0.0, 0.0)];
        let r = lip_pass(&grid, &cents, &areas, &p0, &lp(), 0.0, 2.0, &plate, &kind, &poles, &mut th).unwrap();
        assert_eq!(r.erupting, 1);
        assert_eq!(r.continental_m3, 0.0);
        let east = nearest_cell(&grid, crate::rotation::lat_lon_to_unit(0.0, 9.0).unwrap()).unwrap();
        let west = nearest_cell(&grid, crate::rotation::lat_lon_to_unit(0.0, -9.0).unwrap()).unwrap();
        // The plateau rim is closer to the west sample: less crust there.
        let far_e = nearest_cell(&grid, crate::rotation::lat_lon_to_unit(0.0, 7.5).unwrap()).unwrap();
        let far_w = nearest_cell(&grid, crate::rotation::lat_lon_to_unit(0.0, -6.5).unwrap()).unwrap();
        assert!(th[far_e] >= th[far_w], "{} {}", th[far_e], th[far_w]);
        assert!(th[east] >= th[west]);
    }

    #[test]
    fn swell_is_mantle_fixed_and_decays_to_the_tail() {
        let grid = GridSpec::new(90, 180, 6_371_000.0).unwrap();
        let cents = centroids(&grid).unwrap();
        let s = sp();
        let pl = [plume(10.0, 20.0, 30.0)];
        let before = swell_term(&cents, grid.radius_m(), &pl, &s, 9.0);
        assert!(before.iter().all(|&v| v == 0.0));
        let at = swell_term(&cents, grid.radius_m(), &pl, &s, 10.0);
        let late = swell_term(&cents, grid.radius_m(), &pl, &s, 110.0);
        let c0 = nearest_cell(&grid, pl[0].position).unwrap();
        let max_at = at.iter().cloned().fold(0.0, f64::max);
        assert!(max_at <= 1500.0 + 1e-9 && max_at > 1300.0, "{max_at}");
        assert!(at[c0] > 1300.0, "{}", at[c0]);
        let max_late = late.iter().cloned().fold(0.0, f64::max);
        assert!(max_late <= 800.0 + 7e-4 * 700.0 && max_late > 650.0, "{max_late}");
        // Cut at 3 R: nothing beyond 2400 km; and far smaller beyond 1000 km.
        for k in 0..cents.len() {
            let d = crate::rotation::great_circle_rad(pl[0].position, cents[k]).unwrap() * grid.radius_m();
            if d > 2_400_000.0 + 1.0 {
                assert_eq!(at[k], 0.0);
            }
            if d > 1_200_000.0 {
                assert!(at[k] < 0.15 * max_at, "{d} {}", at[k]);
            }
        }
    }

    #[test]
    fn parameters_validate() {
        assert!(lp().validate().is_ok());
        assert!(sp().validate().is_ok());
        let mut bad = lp();
        bad.taper_fraction = 0.0;
        assert!(bad.validate().is_err());
        let mut bad = sp();
        bad.head_decay_my = 0.0;
        assert!(bad.validate().is_err());
    }
}
