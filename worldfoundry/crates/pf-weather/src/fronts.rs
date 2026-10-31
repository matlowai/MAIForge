//! WEATHER-V1 fronts (lane WEATHER-V1, 2026-10-21; opt-in `fronts=1`).
//!
//! Replaces the round synoptic discs of WEATHER-V0 with keyed frontal
//! systems that each carry a family of curved rain bands (hundreds of km
//! long, tens of km wide), and the independent convective discs with keyed
//! clusters of convective cells lined up along the steering wind. Both stay
//! keyed on the planet-fixed cube lattice and the time slot (never a chart),
//! drift with the latitude steering wind exactly as the V0 discs do, and keep
//! the V0 calibration by construction:
//!
//! - Wet-hour fraction (Poisson void probability). A front system covers a
//!   point when the point lies inside ANY of its bands. The bands of one
//!   system are concentric annulus sectors at disjoint radial ranges, so they
//!   never overlap and the union area of a system is exactly
//!   sum_k L_k w_k (planar, in the azimuthal-equidistant plane of the system
//!   centre). The systems are a marked Poisson process; whether one covers a
//!   point depends only on its own marks, so the number of covering systems is
//!   Poisson with mean m = rho E[sum_k L_k w_k] E[duration], and keeping each
//!   system at a point when its uniform mark < theta(point) gives
//!   P(dry) = exp(-m theta). The same argument holds for convective clusters:
//!   a cluster's cells sit in disjoint cells of a grid (cell radius <= half the
//!   spacing minus the jitter), so at any instant at most one cell of a cluster
//!   covers a point and m = rho E[sum_cells area x duration].
//! - Mean rate (Campbell). Every rain footprint carries a profile with mean 1
//!   over its footprint and lifetime: bands (15/8)(1 - a^2)^2 across (a in
//!   [-1, 1]) times (pi/2) sin(pi s) along (s in [0, 1]); cells 2(1 - d^2);
//!   lifetimes (pi/2) sin(pi tau); intensities lognormal with mean 1. So the
//!   V0 rate law (r_s, r_c from the T0 annual precipitation) carries over.
//!
//! Soft edges: the rain rate falls to zero at a band's edges and ends and at a
//! cell's rim; cloud extends a halo beyond the bands.
//!
//! The optional orographic term (`oro=1`) multiplies the frontal rate by
//! 1 + a tanh(w / w0), w = (event wind anomaly) . grad(T0 elevation): the
//! windward side of a T0-scale relief gets more of the system's rain, the lee
//! less. The anomaly is the synoptic (cyclonic) wind of the V0 systems plus the
//! keyed jitter (the prevailing wind's climatological windward effect is
//! already in the T0 orographic precipitation). It is mean-preserving only
//! approximately (the cyclonic wind of randomly placed systems averages near
//! zero at a point); the calibration table measures it.
//!
//! Wind, synoptic activity and the warm/cold swing are NOT tied to the bands:
//! WEATHER-V1 keeps the V0 synoptic discs for those (see `state_v1`), so the
//! calm-night climate COLDPOOL depends on is that of V0. Only rain and cloud
//! come from the bands and clusters.
//!
//! Spherical caveat: shapes are evaluated in the azimuthal-equidistant plane
//! of their centre; its area scale is (d/R)/sin(d/R) >= 1, so the true area of
//! a band is up to ~1 % smaller than the planar area used in the density
//! (reach <= 1500 km). Reported, not corrected.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::rc::Rc;

use serde::Serialize;
use vz_keyed_rng::Scope;

use crate::constants::ConstantRow;
use crate::disturb::{steering_m_s, unit, SHELL_FRAC};
use crate::sun::DEG;

pub const FRONTS_MODEL_ID: &str = "pf-weather-v1-fronts";

/// The pinned constants of the V1 fronts and convective-cluster processes.
#[derive(Debug, Clone, PartialEq)]
pub struct FrontsConstantsV1 {
    // ---- frontal systems
    pub front_lattice_km: f64,
    pub front_slot_h: f64,
    pub front_d_min_h: f64,
    pub front_d_max_h: f64,
    pub front_mean_cover: f64,
    pub front_intensity_sigma: f64,
    pub bands_min: f64,
    pub bands_max: f64,
    pub band_len_min_km: f64,
    pub band_len_max_km: f64,
    pub band_w_min_km: f64,
    pub band_w_max_km: f64,
    pub band_gap_min_km: f64,
    pub band_gap_max_km: f64,
    pub curv_radius_min_km: f64,
    pub curv_radius_max_km: f64,
    pub axis_bearing_spread_deg: f64,
    pub axis_bearing_mean_deg: f64,
    pub cloud_halo_km: f64,
    pub dry_front_cloud: f64,
    pub cloud_activity_floor: f64,
    // ---- convective clusters
    pub cl_lattice_km: f64,
    pub cl_slot_h: f64,
    pub cl_d_min_h: f64,
    pub cl_d_max_h: f64,
    pub cl_mean_cover: f64,
    pub cl_intensity_sigma: f64,
    pub cl_kx_min: f64,
    pub cl_kx_max: f64,
    pub cl_ky_min: f64,
    pub cl_ky_max: f64,
    pub cl_spacing_km: f64,
    pub cl_occupancy: f64,
    pub cell_r_min_km: f64,
    pub cell_r_max_km: f64,
    pub cell_d_min_h: f64,
    pub cell_d_max_h: f64,
    pub cell_halo: f64,
    pub cl_axis_spread_deg: f64,
    pub cl_anvil_cloud: f64,
    pub cl_anvil_halo_km: f64,
    // ---- orographic modulation (oro=1)
    pub oro_amp: f64,
    pub oro_w_scale_m_s: f64,
    pub oro_grad_step_cells: f64,
    // ---- soft edges
    pub rain_soft_mm_h: f64,
}

impl Default for FrontsConstantsV1 {
    fn default() -> Self {
        FrontsConstantsV1 {
            front_lattice_km: 1800.0,
            front_slot_h: 12.0,
            front_d_min_h: 6.0,
            front_d_max_h: 48.0,
            front_mean_cover: 1.2,
            front_intensity_sigma: 0.6,
            bands_min: 2.0,
            bands_max: 4.0,
            band_len_min_km: 300.0,
            band_len_max_km: 1200.0,
            band_w_min_km: 20.0,
            band_w_max_km: 80.0,
            band_gap_min_km: 40.0,
            band_gap_max_km: 120.0,
            curv_radius_min_km: 600.0,
            curv_radius_max_km: 3000.0,
            axis_bearing_spread_deg: 20.0,
            axis_bearing_mean_deg: 30.0,
            cloud_halo_km: 60.0,
            dry_front_cloud: 0.25,
            cloud_activity_floor: 0.3,
            cl_lattice_km: 180.0,
            cl_slot_h: 3.0,
            cl_d_min_h: 3.0,
            cl_d_max_h: 8.0,
            cl_mean_cover: 0.8,
            cl_intensity_sigma: 0.8,
            cl_kx_min: 2.0,
            cl_kx_max: 7.0,
            cl_ky_min: 1.0,
            cl_ky_max: 3.0,
            cl_spacing_km: 32.0,
            cl_occupancy: 0.6,
            cell_r_min_km: 5.0,
            cell_r_max_km: 15.0,
            cell_d_min_h: 0.5,
            cell_d_max_h: 2.5,
            cell_halo: 2.0,
            cl_axis_spread_deg: 25.0,
            cl_anvil_cloud: 0.5,
            cl_anvil_halo_km: 20.0,
            oro_amp: 0.6,
            oro_w_scale_m_s: 0.05,
            oro_grad_step_cells: 0.25,
            rain_soft_mm_h: 0.5,
        }
    }
}

const R: &str = "recalled";
const P: &str = "project_assumption";

impl FrontsConstantsV1 {
    pub fn table(&self) -> Vec<ConstantRow> {
        let k = self;
        let row = |name, value, unit, class, basis| ConstantRow { name, value, unit, class, basis };
        vec![
            row("front_lattice_km", k.front_lattice_km, "km", P, "planet-fixed cube lattice keying frontal systems; must exceed the system reach (bands + halo) + half the shell (checked)"),
            row("front_slot_h", k.front_slot_h, "h", P, "time slot of the frontal key"),
            row("front_d_min_h", k.front_d_min_h, "h", R, "lifetime of a frontal system's rain at the chart scale 6-48 h (as V0 strat)"),
            row("front_d_max_h", k.front_d_max_h, "h", R, "see front_d_min_h"),
            row("front_mean_cover", k.front_mean_cover, "1", P, "mean number of systems whose band union covers a point (before thinning); >= -ln(1 - wet_hour_fraction_max)"),
            row("front_intensity_sigma", k.front_intensity_sigma, "1", P, "lognormal (mean 1) spread of a system's intensity"),
            row("bands_min", k.bands_min, "1", P, "rain bands per system, uniform integer in [min, max]: cold-frontal, prefrontal, warm-frontal and postfrontal bands (recalled qualitatively: several bands per extratropical cyclone)"),
            row("bands_max", k.bands_max, "1", P, "see bands_min"),
            row("band_len_min_km", k.band_len_min_km, "km", R, "frontal rain band length hundreds of km (operator spec 2026-10-21; recalled 300-1500 km)"),
            row("band_len_max_km", k.band_len_max_km, "km", R, "see band_len_min_km"),
            row("band_w_min_km", k.band_w_min_km, "km", R, "band width tens of km (operator spec; narrow cold-frontal bands ~20 km, wide warm-frontal bands ~80 km, recalled)"),
            row("band_w_max_km", k.band_w_max_km, "km", R, "see band_w_min_km"),
            row("band_gap_min_km", k.band_gap_min_km, "km", P, "dry gap between neighbouring bands of one system"),
            row("band_gap_max_km", k.band_gap_max_km, "km", P, "see band_gap_min_km"),
            row("curv_radius_min_km", k.curv_radius_min_km, "km", P, "radius of curvature of the band family (uniform); sign keyed"),
            row("curv_radius_max_km", k.curv_radius_max_km, "km", P, "see curv_radius_min_km"),
            row("axis_bearing_spread_deg", k.axis_bearing_spread_deg, "deg", P, "band axis bearing uniform within +-spread of the mean (50 deg gave a cross-hatched sky of crossing families, lane figure seq_draft1)"),
            row("axis_bearing_mean_deg", k.axis_bearing_mean_deg, "deg", R, "mean band axis bearing from north: cold fronts trail NE-SW in the northern mid-latitudes (NW-SE in the south), ramped to N-S (easterly-wave troughs) equatorward of 10-30 deg (recalled qualitatively)"),
            row("cloud_halo_km", k.cloud_halo_km, "km", P, "cloud extends this far outside a system's bands (smoothstep falloff); 60 km keeps the site-1 cloud climate of V0 (mean cloud 0.598 vs 0.596, year 0; 120 km gave 0.658 and cut clear nights by 4x)"),
            row("dry_front_cloud", k.dry_front_cloud, "1", P, "cloud of a system that is thinned out (no rain) at this place"),
            row("cloud_activity_floor", k.cloud_activity_floor, "1", P, "band cloud is scaled by floor + (1 - floor) * V0 synoptic activity: clear skies stay correlated with calm high pressure (without it site-1 calm clear night hours fell from 442 to 137 per year)"),
            row("cl_lattice_km", k.cl_lattice_km, "km", P, "cube lattice keying convective clusters (> cluster reach + half shell, checked)"),
            row("cl_slot_h", k.cl_slot_h, "h", P, "time slot of the cluster key"),
            row("cl_d_min_h", k.cl_d_min_h, "h", R, "cluster (mesoscale convective system) lifetime 3-8 h (recalled order); must be >= cell_d_max_h"),
            row("cl_d_max_h", k.cl_d_max_h, "h", R, "see cl_d_min_h"),
            row("cl_mean_cover", k.cl_mean_cover, "1", P, "mean clusters with a cell over a point (before thinning); >= share_max (1 + diurnal_amp) * -ln(1 - wet_max)"),
            row("cl_intensity_sigma", k.cl_intensity_sigma, "1", P, "lognormal (mean 1) spread of a cluster's intensity"),
            row("cl_kx_min", k.cl_kx_min, "1", P, "cell grid columns along the steering wind, uniform integer"),
            row("cl_kx_max", k.cl_kx_max, "1", P, "see cl_kx_min"),
            row("cl_ky_min", k.cl_ky_min, "1", P, "cell grid rows across the steering wind, uniform integer (clusters line up along the wind, design review 8.1 [V])"),
            row("cl_ky_max", k.cl_ky_max, "1", P, "see cl_ky_min"),
            row("cl_spacing_km", k.cl_spacing_km, "km", P, "cell grid spacing; cell radius + jitter <= spacing / 2 keeps cells disjoint"),
            row("cl_occupancy", k.cl_occupancy, "1", P, "probability a grid slot holds a cell"),
            row("cell_r_min_km", k.cell_r_min_km, "km", R, "convective cell rain radius 5-15 km (recalled order)"),
            row("cell_r_max_km", k.cell_r_max_km, "km", R, "see cell_r_min_km"),
            row("cell_d_min_h", k.cell_d_min_h, "h", R, "a cell rains 0.5-2.5 h (recalled order)"),
            row("cell_d_max_h", k.cell_d_max_h, "h", R, "see cell_d_min_h"),
            row("cell_halo", k.cell_halo, "1", P, "cell cloud radius as a multiple of its rain radius"),
            row("cl_axis_spread_deg", k.cl_axis_spread_deg, "deg", P, "cluster axis within +-spread of the zonal steering direction"),
            row("cl_anvil_cloud", k.cl_anvil_cloud, "1", P, "anvil/stratiform cloud over the whole cell grid of a kept, active cluster"),
            row("cl_anvil_halo_km", k.cl_anvil_halo_km, "km", P, "anvil cloud falloff outside the elliptical cell-grid envelope (semi-axes half extents + half a spacing)"),
            row("oro_amp", k.oro_amp, "1", P, "orographic modulation amplitude of frontal rain: factor 1 + a tanh(w / w0) in [0.4, 1.6] (design review 8.2 [V]: windward 2-5x at relief scale; this is the T0-scale share)"),
            row("oro_w_scale_m_s", k.oro_w_scale_m_s, "m/s", P, "forced vertical velocity scale w0 of the modulation (5 m/s over a 1 % T0 slope)"),
            row("oro_grad_step_cells", k.oro_grad_step_cells, "1", P, "finite-difference half step of the T0 elevation gradient, in T0 cells"),
            row("rain_soft_mm_h", k.rain_soft_mm_h, "mm/h", P, "rain rate over which the rain cloud (>= 0.9) and saturation (RH >= 0.95) ramp in (smoothstep), so visibility has no step at a footprint edge (V0 steps at any rate > 0)"),
        ]
    }

    pub fn digest(&self) -> String {
        let bytes = serde_json::to_vec(&self.table()).expect("constants serialize");
        format!("blake3:{}", blake3::hash(&bytes).to_hex())
    }

    fn shell(&self, l: f64, planet_r_km: f64) -> f64 {
        let h = SHELL_FRAC * l;
        h * (1.0 + h * h / (12.0 * planet_r_km * planet_r_km))
    }

    /// Upper bound of the distance from a system centre to any point it can
    /// influence (bands and cloud halo).
    pub fn front_reach_max_km(&self) -> f64 {
        let n = self.bands_max;
        (n - 1.0) * (self.band_w_max_km + self.band_gap_max_km) + 0.5 * self.band_w_max_km + 0.5 * self.band_len_max_km
            + self.cloud_halo_km
    }

    pub fn cluster_reach_max_km(&self) -> f64 {
        let cells = 0.5 * self.cl_spacing_km * (self.cl_kx_max.hypot(self.cl_ky_max)) + self.cell_halo * self.cell_r_max_km;
        let anvil = 0.5 * self.cl_spacing_km * (self.cl_kx_max.max(self.cl_ky_max) + 1.0) + self.cl_anvil_halo_km;
        cells.max(anvil)
    }

    /// E[union rain area x duration] of one frontal system (km^2 h).
    pub fn front_cover_per_system(&self) -> f64 {
        let en = 0.5 * (self.bands_min + self.bands_max);
        let el = 0.5 * (self.band_len_min_km + self.band_len_max_km);
        let ew = 0.5 * (self.band_w_min_km + self.band_w_max_km);
        let ed = 0.5 * (self.front_d_min_h + self.front_d_max_h);
        en * el * ew * ed
    }

    /// E[sum over cells of area x duration] of one cluster (km^2 h).
    pub fn cluster_cover_per_system(&self) -> f64 {
        let ekx = 0.5 * (self.cl_kx_min + self.cl_kx_max);
        let eky = 0.5 * (self.cl_ky_min + self.cl_ky_max);
        let (a, b) = (self.cell_r_min_km, self.cell_r_max_km);
        let er2 = (a * a + a * b + b * b) / 3.0;
        let ed = 0.5 * (self.cell_d_min_h + self.cell_d_max_h);
        self.cl_occupancy * ekx * eky * PI * er2 * ed
    }

    pub fn mu_front(&self, planet_r_km: f64) -> f64 {
        let rho = self.front_mean_cover / self.front_cover_per_system();
        let l = self.front_lattice_km;
        rho * self.front_slot_h * l * l * l / self.shell(l, planet_r_km)
    }

    pub fn mu_cluster(&self, planet_r_km: f64) -> f64 {
        let rho = self.cl_mean_cover / self.cluster_cover_per_system();
        let l = self.cl_lattice_km;
        rho * self.cl_slot_h * l * l * l / self.shell(l, planet_r_km)
    }

    pub fn check(&self) -> Result<(), String> {
        let fr = self.front_reach_max_km() * 1.02;
        if fr + 0.5 * SHELL_FRAC * self.front_lattice_km >= self.front_lattice_km {
            return Err(format!("front reach {fr} km does not fit lattice {}", self.front_lattice_km));
        }
        let cr = self.cluster_reach_max_km() * 1.02;
        if cr + 0.5 * SHELL_FRAC * self.cl_lattice_km >= self.cl_lattice_km {
            return Err(format!("cluster reach {cr} km does not fit lattice {}", self.cl_lattice_km));
        }
        if self.cell_r_max_km > 0.5 * self.cl_spacing_km {
            return Err("cell radius exceeds half the cluster spacing (cells could overlap)".into());
        }
        if self.cl_d_min_h < self.cell_d_max_h {
            return Err("cluster lifetime shorter than a cell lifetime".into());
        }
        if self.curv_radius_min_km <= 0.5 * self.band_len_max_km / 1.2 {
            return Err("curvature radius too small for the band length".into());
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Band {
    pub rho_km: f64,
    pub w_km: f64,
    pub len_km: f64,
    pub half_theta: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FrontSystem {
    pub centre: [f64; 3],
    pub east: [f64; 3],
    pub north: [f64; 3],
    pub t0_h: f64,
    pub dur_h: f64,
    pub strength: f64,
    pub intensity: f64,
    /// Band axis and across unit vectors in the centre's (east, north).
    pub axis: [f64; 2],
    pub across: [f64; 2],
    /// Curvature side (+1 / -1) and the radius of band 0.
    pub sigma: f64,
    pub rho0_km: f64,
    pub bands: Vec<Band>,
    pub reach_km: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Cell {
    pub u_km: f64,
    pub v_km: f64,
    pub r_km: f64,
    pub t0_h: f64,
    pub dur_h: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Cluster {
    pub centre: [f64; 3],
    pub east: [f64; 3],
    pub north: [f64; 3],
    pub t0_h: f64,
    pub dur_h: f64,
    pub strength: f64,
    pub intensity: f64,
    pub axis: [f64; 2],
    pub across: [f64; 2],
    pub cells: Vec<Cell>,
    /// Half extents of the cell grid along / across the axis (km).
    pub half_u_km: f64,
    pub half_v_km: f64,
    pub reach_km: f64,
}

/// Aggregate of the frontal systems at one point and hour.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FrontAgg {
    pub sum_rain: f64,
    pub n_rain: u32,
    pub cloud: f64,
    pub n_systems_near: u32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ClusterAgg {
    pub sum_rain: f64,
    pub n_rain: u32,
    pub cloud: f64,
    pub gust_intensity: f64,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
fn norm(a: [f64; 3]) -> f64 { dot(a, a).sqrt() }

fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Local east and north unit vectors at a unit vector.
fn basis(c: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let lat = c[2].clamp(-1.0, 1.0).asin();
    let lon = c[1].atan2(c[0]);
    ([-lon.sin(), lon.cos(), 0.0], [-lat.sin() * lon.cos(), -lat.sin() * lon.sin(), lat.cos()])
}

/// Azimuthal-equidistant (east, north) km of x around centre c.
fn aeqd(x: [f64; 3], c: [f64; 3], east: [f64; 3], north: [f64; 3], r_km: f64) -> [f64; 2] {
    let cd = dot(x, c);
    let ang = norm(cross(x, c)).atan2(cd);
    let d = [x[0] - cd * c[0], x[1] - cd * c[1], x[2] - cd * c[2]];
    let dn = norm(d);
    if dn < 1e-15 { return [0.0, 0.0]; }
    let s = ang * r_km / dn;
    [dot(d, east) * s, dot(d, north) * s]
}

pub struct FrontField {
    law: Scope,
    pub planet_r_km: f64,
    pub steer_tropical: f64,
    pub steer_midlat: f64,
    pub k: FrontsConstantsV1,
    memo_f: RefCell<BTreeMap<(i64, i64, i64, i64), Rc<Vec<FrontSystem>>>>,
    memo_c: RefCell<BTreeMap<(i64, i64, i64, i64), Rc<Vec<Cluster>>>>,
}

impl FrontField {
    pub fn new(law: Scope, planet_r_km: f64, steer_tropical: f64, steer_midlat: f64, k: FrontsConstantsV1) -> Result<FrontField, String> {
        k.check()?;
        Ok(FrontField { law, planet_r_km, steer_tropical, steer_midlat, k, memo_f: RefCell::new(BTreeMap::new()), memo_c: RefCell::new(BTreeMap::new()) })
    }

    pub fn clear_memo(&self) {
        self.memo_f.borrow_mut().clear();
        self.memo_c.borrow_mut().clear();
    }

    fn omega_rad_h(&self, lat_deg: f64) -> f64 {
        let u = steering_m_s(lat_deg, self.steer_tropical, self.steer_midlat);
        let c = (lat_deg * DEG).cos().max(0.1);
        u * 3600.0 / (self.planet_r_km * 1000.0 * c)
    }

    /// True when the cube can hold a centre inside the shell.
    fn cube_in_shell(&self, cube: [i64; 3], l: f64) -> bool {
        let lo = [cube[0] as f64 * l, cube[1] as f64 * l, cube[2] as f64 * l];
        let (mut dmin2, mut dmax2) = (0.0, 0.0);
        for a in 0..3 {
            let (c0, c1) = (lo[a], lo[a] + l);
            let near = if c0 > 0.0 { c0 } else if c1 < 0.0 { -c1 } else { 0.0 };
            let far = c0.abs().max(c1.abs());
            dmin2 += near * near;
            dmax2 += far * far;
        }
        let hs = 0.5 * SHELL_FRAC * l;
        !(dmin2.sqrt() > self.planet_r_km + hs || dmax2.sqrt() < self.planet_r_km - hs)
    }

    fn poisson(u: f64, mu: f64) -> u64 {
        let (mut k, mut pk) = (0u64, (-mu).exp());
        let mut cdf = pk;
        while u >= cdf && k < 1023 {
            k += 1;
            pk *= mu / k as f64;
            cdf += pk;
        }
        k
    }

    pub fn fronts(&self, cube: [i64; 3], slot: i64) -> Rc<Vec<FrontSystem>> {
        let key = (cube[0], cube[1], cube[2], slot);
        if let Some(v) = self.memo_f.borrow().get(&key) { return v.clone(); }
        let v = Rc::new(self.gen_fronts(cube, slot));
        self.memo_f.borrow_mut().insert(key, v.clone());
        v
    }

    fn gen_fronts(&self, cube: [i64; 3], slot: i64) -> Vec<FrontSystem> {
        let k = &self.k;
        let mut out = Vec::new();
        let l = k.front_lattice_km;
        if slot < 0 || !self.cube_in_shell(cube, l) { return out; }
        let r = self.planet_r_km;
        let scope = self.law.address(&format!("front/{}/{}/{}", cube[0], cube[1], cube[2])).expect("address");
        let draw = |what: &str, counter: u64| scope.draw(&format!("weather.front.{what}"), counter, 0).expect("purpose");
        let n = Self::poisson(draw("count", slot as u64).uniform_f64(), k.mu_front(r));
        let hs = 0.5 * SHELL_FRAC * l;
        let lo = [cube[0] as f64 * l, cube[1] as f64 * l, cube[2] as f64 * l];
        for i in 0..n {
            let c = slot as u64 * 1024 + i;
            let x = lo[0] + l * draw("cx", c).uniform_f64();
            let y = lo[1] + l * draw("cy", c).uniform_f64();
            let z = lo[2] + l * draw("cz", c).uniform_f64();
            let rr = (x * x + y * y + z * z).sqrt();
            if (rr - r).abs() >= hs || rr == 0.0 { continue; }
            let centre = [x / rr, y / rr, z / rr];
            let (east, north) = basis(centre);
            let t0 = (slot as f64 + draw("start", c).uniform_f64()) * k.front_slot_h;
            let dur = k.front_d_min_h + (k.front_d_max_h - k.front_d_min_h) * draw("dur", c).uniform_f64();
            let strength = draw("thin", c).uniform_f64();
            let s = k.front_intensity_sigma;
            let intensity = (s * draw("intensity", c).normal_f64() - 0.5 * s * s).exp();
            let clat = centre[2].clamp(-1.0, 1.0).asin() / DEG;
            let mean = k.axis_bearing_mean_deg * clat.signum() * smoothstep(10.0, 30.0, clat.abs());
            let phi = (mean + (2.0 * draw("bearing", c).uniform_f64() - 1.0) * k.axis_bearing_spread_deg) * DEG;
            let axis = [phi.sin(), phi.cos()];
            let across = [phi.cos(), -phi.sin()];
            let sigma = if draw("side", c).uniform_f64() < 0.5 { 1.0 } else { -1.0 };
            let rho0 = k.curv_radius_min_km + (k.curv_radius_max_km - k.curv_radius_min_km) * draw("curv", c).uniform_f64();
            let nb = k.bands_min as u64 + draw("nbands", c).below((k.bands_max - k.bands_min) as u64 + 1).expect("below");
            let mut bands = Vec::with_capacity(nb as usize);
            let mut rho = rho0;
            let mut prev_w = 0.0;
            let mut max_len: f64 = 0.0;
            for b in 0..nb {
                let cb = c * 8 + b;
                let w = k.band_w_min_km + (k.band_w_max_km - k.band_w_min_km) * draw("band-w", cb).uniform_f64();
                let len = k.band_len_min_km + (k.band_len_max_km - k.band_len_min_km) * draw("band-len", cb).uniform_f64();
                if b > 0 {
                    let gap = k.band_gap_min_km + (k.band_gap_max_km - k.band_gap_min_km) * draw("band-gap", cb).uniform_f64();
                    rho += 0.5 * (prev_w + w) + gap;
                }
                prev_w = w;
                max_len = max_len.max(len);
                bands.push(Band { rho_km: rho, w_km: w, len_km: len, half_theta: 0.5 * len / rho });
            }
            let reach = ((rho - rho0) + 0.5 * prev_w + 0.5 * max_len + k.cloud_halo_km) * 1.02;
            out.push(FrontSystem { centre, east, north, t0_h: t0, dur_h: dur, strength, intensity, axis, across, sigma, rho0_km: rho0, bands, reach_km: reach });
        }
        out
    }

    pub fn clusters(&self, cube: [i64; 3], slot: i64) -> Rc<Vec<Cluster>> {
        let key = (cube[0], cube[1], cube[2], slot);
        if let Some(v) = self.memo_c.borrow().get(&key) { return v.clone(); }
        let v = Rc::new(self.gen_clusters(cube, slot));
        self.memo_c.borrow_mut().insert(key, v.clone());
        v
    }

    fn gen_clusters(&self, cube: [i64; 3], slot: i64) -> Vec<Cluster> {
        let k = &self.k;
        let mut out = Vec::new();
        let l = k.cl_lattice_km;
        if slot < 0 || !self.cube_in_shell(cube, l) { return out; }
        let r = self.planet_r_km;
        let scope = self.law.address(&format!("convcl/{}/{}/{}", cube[0], cube[1], cube[2])).expect("address");
        let draw = |what: &str, counter: u64| scope.draw(&format!("weather.convcl.{what}"), counter, 0).expect("purpose");
        let n = Self::poisson(draw("count", slot as u64).uniform_f64(), k.mu_cluster(r));
        let hs = 0.5 * SHELL_FRAC * l;
        let lo = [cube[0] as f64 * l, cube[1] as f64 * l, cube[2] as f64 * l];
        let sp = k.cl_spacing_km;
        for i in 0..n {
            let c = slot as u64 * 1024 + i;
            let x = lo[0] + l * draw("cx", c).uniform_f64();
            let y = lo[1] + l * draw("cy", c).uniform_f64();
            let z = lo[2] + l * draw("cz", c).uniform_f64();
            let rr = (x * x + y * y + z * z).sqrt();
            if (rr - r).abs() >= hs || rr == 0.0 { continue; }
            let centre = [x / rr, y / rr, z / rr];
            let (east, north) = basis(centre);
            let strength = draw("thin", c).uniform_f64();
            let s = k.cl_intensity_sigma;
            let intensity = (s * draw("intensity", c).normal_f64() - 0.5 * s * s).exp();
            // Shape word: 16 bits start, 16 bits duration, 16 bits axis angle, 8 + 8 bits grid size.
            let w = draw("shape", c).u64();
            let f16 = |sh: u32| ((w >> sh) & 0xFFFF) as f64 / 65536.0;
            let t0 = (slot as f64 + f16(0)) * k.cl_slot_h;
            let dur = k.cl_d_min_h + (k.cl_d_max_h - k.cl_d_min_h) * f16(16);
            let ang = (2.0 * f16(32) - 1.0) * k.cl_axis_spread_deg * DEG;
            let kx = k.cl_kx_min as u64 + ((w >> 48) & 0xFF) % ((k.cl_kx_max - k.cl_kx_min) as u64 + 1);
            let ky = k.cl_ky_min as u64 + ((w >> 56) & 0xFF) % ((k.cl_ky_max - k.cl_ky_min) as u64 + 1);
            // Axis along the zonal steering direction (east), rotated by ang.
            let axis = [ang.cos(), ang.sin()];
            let across = [-ang.sin(), ang.cos()];
            let mut cells = Vec::new();
            for j in 0..(kx * ky) {
                // Cell word: 10 bits occupancy, 10 + 10 jitter, 10 radius, 12 start, 12 duration.
                let cw = draw("cell", c * 64 + j).u64();
                let f = |sh: u32, bits: u32| ((cw >> sh) & ((1u64 << bits) - 1)) as f64 / (1u64 << bits) as f64;
                if f(0, 10) >= k.cl_occupancy { continue; }
                let rad = k.cell_r_min_km + (k.cell_r_max_km - k.cell_r_min_km) * f(30, 10);
                let jmax = 0.5 * sp - rad;
                let (ix, iy) = ((j % kx) as f64, (j / kx) as f64);
                let u = (ix - 0.5 * (kx as f64 - 1.0)) * sp + (2.0 * f(10, 10) - 1.0) * jmax;
                let v = (iy - 0.5 * (ky as f64 - 1.0)) * sp + (2.0 * f(20, 10) - 1.0) * jmax;
                let cd = k.cell_d_min_h + (k.cell_d_max_h - k.cell_d_min_h) * f(52, 12);
                let ct0 = t0 + (dur - cd) * f(40, 12);
                cells.push(Cell { u_km: u, v_km: v, r_km: rad, t0_h: ct0, dur_h: cd });
            }
            let reach_cells = 0.5 * sp * (kx as f64).hypot(ky as f64) + k.cell_halo * k.cell_r_max_km;
            let reach_anvil = 0.5 * sp * ((kx.max(ky)) as f64 + 1.0) + k.cl_anvil_halo_km;
            let reach = reach_cells.max(reach_anvil) * 1.02;
            out.push(Cluster { centre, east, north, t0_h: t0, dur_h: dur, strength, intensity, axis, across, cells,
                               half_u_km: 0.5 * sp * kx as f64, half_v_km: 0.5 * sp * ky as f64, reach_km: reach });
        }
        out
    }

    /// Visit every keyed item of a process whose slot could be alive at t,
    /// with the point back-advected into each slot frame.
    fn visit<T, G, F>(&self, lat_deg: f64, lon_deg: f64, t_h: f64, lattice: f64, slot_h: f64, d_max: f64, gen: G, mut f: F)
    where
        G: Fn(&Self, [i64; 3], i64) -> Rc<Vec<T>>,
        F: FnMut(&T, [f64; 3], [f64; 3], [f64; 3]),
    {
        let r = self.planet_r_km;
        let s_hi = (t_h / slot_h).floor() as i64;
        let s_lo = ((t_h - d_max) / slot_h).floor() as i64;
        let omega = self.omega_rad_h(lat_deg);
        for s in s_lo.max(0)..=s_hi {
            let dt = t_h - s as f64 * slot_h;
            let lon_b = lon_deg - omega * dt / DEG;
            let x = unit(lat_deg, lon_b);
            let (ex, nx) = basis(x);
            let cube = [(x[0] * r / lattice).floor() as i64, (x[1] * r / lattice).floor() as i64, (x[2] * r / lattice).floor() as i64];
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let list = gen(self, [cube[0] + dx, cube[1] + dy, cube[2] + dz], s);
                        for item in list.iter() { f(item, x, ex, nx); }
                    }
                }
            }
        }
    }

    /// Frontal systems at a point and hour; `theta` is the keep threshold there.
    pub fn eval_fronts(&self, lat_deg: f64, lon_deg: f64, t_h: f64, theta: f64) -> FrontAgg {
        let k = &self.k;
        let r = self.planet_r_km;
        let mut agg = FrontAgg::default();
        self.visit(lat_deg, lon_deg, t_h, k.front_lattice_km, k.front_slot_h, k.front_d_max_h, |s: &Self, c, sl| s.fronts(c, sl), |f: &FrontSystem, x, _ex, _nx| {
            if t_h < f.t0_h || t_h >= f.t0_h + f.dur_h { return; }
            let cd = dot(x, f.centre);
            if cd < (f.reach_km / r).cos() { return; }
            agg.n_systems_near += 1;
            let p = aeqd(x, f.centre, f.east, f.north, r);
            let u = p[0] * f.axis[0] + p[1] * f.axis[1];
            let v = p[0] * f.across[0] + p[1] * f.across[1];
            let yo = f.sigma * (v + f.sigma * f.rho0_km);
            let rr = u.hypot(yo);
            let alpha = u.abs().atan2(yo);
            let tau = (t_h - f.t0_h) / f.dur_h;
            let env = (PI * tau).sin();
            let kept = f.strength < theta;
            let mut d_out = f64::INFINITY;
            for b in &f.bands {
                let dr = (rr - b.rho_km).abs() - 0.5 * b.w_km;
                let da = (alpha - b.half_theta) * rr;
                if dr < 0.0 && da < 0.0 {
                    d_out = 0.0;
                    if kept {
                        let a = (rr - b.rho_km) / (0.5 * b.w_km);
                        let ga = 15.0 / 8.0 * (1.0 - a * a) * (1.0 - a * a);
                        let sa = 0.5 + 0.5 * alpha / b.half_theta; // alpha in [0, half): s in [0.5, 1); g_s symmetric
                        let gs = 0.5 * PI * (PI * sa).sin();
                        agg.sum_rain += f.intensity * ga * gs * 0.5 * PI * env;
                        agg.n_rain += 1;
                    }
                } else {
                    d_out = d_out.min(dr.max(0.0).hypot(da.max(0.0)));
                }
            }
            let life = smoothstep(0.0, 0.15, tau) * (1.0 - smoothstep(0.85, 1.0, tau));
            let cprof = 1.0 - smoothstep(0.0, k.cloud_halo_km, d_out);
            let cl = if kept { cprof } else { k.dry_front_cloud * cprof };
            agg.cloud = agg.cloud.max(cl * life.max(0.3));
        });
        agg
    }

    /// Convective clusters at a point and hour; `theta` is the keep threshold.
    pub fn eval_clusters(&self, lat_deg: f64, lon_deg: f64, t_h: f64, theta: f64) -> ClusterAgg {
        let k = &self.k;
        let r = self.planet_r_km;
        let mut agg = ClusterAgg::default();
        self.visit(lat_deg, lon_deg, t_h, k.cl_lattice_km, k.cl_slot_h, k.cl_d_max_h, |s: &Self, c, sl| s.clusters(c, sl), |cl: &Cluster, x, _ex, _nx| {
            if t_h < cl.t0_h || t_h >= cl.t0_h + cl.dur_h { return; }
            if cl.strength >= theta { return; }
            if dot(x, cl.centre) < (cl.reach_km / r).cos() { return; }
            let p = aeqd(x, cl.centre, cl.east, cl.north, r);
            let u = p[0] * cl.axis[0] + p[1] * cl.axis[1];
            let v = p[0] * cl.across[0] + p[1] * cl.across[1];
            // Elliptical anvil around the cell grid (semi-axes = half extents + half a spacing).
            let (au, av) = (cl.half_u_km + 0.5 * k.cl_spacing_km, cl.half_v_km + 0.5 * k.cl_spacing_km);
            let rho = (u / au).hypot(v / av);
            let d_ell = (rho - 1.0).max(0.0) * au.min(av);
            agg.cloud = agg.cloud.max(k.cl_anvil_cloud * (1.0 - smoothstep(0.0, k.cl_anvil_halo_km, d_ell)) * (1.0 - 0.5 * smoothstep(0.5, 1.0, rho)));
            for c in &cl.cells {
                if t_h < c.t0_h || t_h >= c.t0_h + c.dur_h { continue; }
                let d = (u - c.u_km).hypot(v - c.v_km) / c.r_km;
                if d >= k.cell_halo { continue; }
                let prof = 1.0 - smoothstep(1.0, k.cell_halo, d);
                agg.cloud = agg.cloud.max(0.75 * prof);
                if d < 1.0 {
                    let tau = (t_h - c.t0_h) / c.dur_h;
                    let g = 2.0 * (1.0 - d * d);
                    let e = 0.5 * PI * (PI * tau).sin();
                    agg.sum_rain += cl.intensity * g * e;
                    agg.n_rain += 1;
                    agg.gust_intensity = cl.intensity;
                }
            }
        });
        agg
    }
}
