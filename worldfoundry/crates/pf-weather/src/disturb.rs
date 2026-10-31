//! Keyed synoptic and convective disturbance processes (random access in
//! space and time, no stateful generator, no chart).
//!
//! Each process is a marked Poisson process of discs on the sphere in
//! space-time. Disturbances are keyed by a planet-fixed 3D cube lattice
//! (Earth-centred coordinates of the unit sphere times the planet radius) and
//! a time slot: scope `address("<proc>/<ix>/<iy>/<iz>")` under the law-set
//! scope, draw counter `slot` (count) or `slot * 256 + k` (marks). Nothing is
//! keyed on a chart (ADR 045 direction): the same planet point at the same
//! hour gets the same disturbances whatever window or chart asks.
//!
//! A disturbance spawned in slot s lives in a frame drifting eastward with the
//! latitude-dependent steering speed u(lat), anchored at the slot start: a
//! point (lat, lon, t) is tested against it after back-advection
//! lon' = lon - Omega(lat) (t - s tau). Rotation about the polar axis
//! preserves area at every instant, so the per-point coverage law below is
//! exact; discs shear into tilted bands as they drift (as fronts do).
//!
//! Coverage law (Poisson void probability): if a point is covered on average
//! by `m` discs (unthinned), and each disc is kept at that point when its
//! uniform mark `strength < theta(point)`, the number of kept covering discs
//! is Poisson(m theta), so P(no kept disc) = exp(-m theta) exactly. The
//! weather sets theta from the wanted wet-hour fraction.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::rc::Rc;
use vz_keyed_rng::Scope;

use crate::sun::DEG;

/// Centres are drawn uniformly in a cube and kept inside the shell
/// |x| in [R - h/2, R + h/2], h = SHELL_FRAC * lattice, then projected to the
/// sphere. Projection moves a centre by at most h/2, so a disc whose halo
/// reaches a point was spawned within halo + h/2 <= lattice of it (checked in
/// `ProcessSpec::check`), i.e. in one of the 27 cubes around the point.
pub const SHELL_FRAC: f64 = 0.25;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProcessSpec {
    pub tag: &'static str,
    pub id: u8,
    pub lattice_km: f64,
    pub slot_h: f64,
    pub r_min_km: f64,
    pub r_max_km: f64,
    pub d_min_h: f64,
    pub d_max_h: f64,
    pub mean_cover: f64,
    pub halo: f64,
    pub intensity_sigma: f64,
}

impl ProcessSpec {
    /// The neighbourhood bound: halo radius + half shell < lattice.
    pub fn check(&self) -> Result<(), String> {
        let need = self.r_max_km * self.halo + 0.5 * SHELL_FRAC * self.lattice_km;
        if need >= self.lattice_km {
            return Err(format!("process {}: halo reach {need} km >= lattice {} km", self.tag, self.lattice_km));
        }
        if !(self.r_min_km > 0.0 && self.r_max_km >= self.r_min_km && self.d_min_h > 0.0 && self.d_max_h >= self.d_min_h && self.slot_h > 0.0) {
            return Err(format!("process {}: bad ranges", self.tag));
        }
        Ok(())
    }

    /// Expected disc count per cube per slot (including centres that fall
    /// outside the surface shell and are discarded).
    pub fn mu_per_cube(&self, planet_r_km: f64) -> f64 {
        let (a, b) = (self.r_min_km, self.r_max_km);
        let e_area = PI * (a * a + a * b + b * b) / 3.0;
        let e_dur = 0.5 * (self.d_min_h + self.d_max_h);
        let rho = self.mean_cover / (e_area * e_dur);
        let l = self.lattice_km;
        let h = SHELL_FRAC * l;
        let shell = 1.0 + h * h / (12.0 * planet_r_km * planet_r_km);
        rho * self.slot_h * l * l * l / (h * shell)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Disturbance {
    /// Unit vector of the centre in the slot frame.
    pub centre: [f64; 3],
    pub t0_h: f64,
    pub dur_h: f64,
    pub radius_km: f64,
    /// Uniform thinning mark in [0, 1).
    pub strength: f64,
    /// Lognormal intensity with mean 1.
    pub intensity: f64,
}

/// One disc covering (or haloing) a point at a time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hit {
    /// Distance from the centre / rain radius (the halo reaches `halo`).
    pub d_frac: f64,
    /// Fraction of the lifetime elapsed, [0, 1).
    pub tau: f64,
    pub strength: f64,
    pub intensity: f64,
    /// Unit cyclonic tangent at the point (east, north); zero at the centre.
    pub cyclonic_en: [f64; 2],
    /// The disc centre at the query time (drifted at the centre's latitude).
    pub centre_lat_deg: f64,
    pub centre_lon_deg: f64,
}

pub fn unit(lat_deg: f64, lon_deg: f64) -> [f64; 3] {
    let (p, l) = (lat_deg * DEG, lon_deg * DEG);
    [p.cos() * l.cos(), p.cos() * l.sin(), p.sin()]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
fn norm(a: [f64; 3]) -> f64 { dot(a, a).sqrt() }

/// Eastward steering speed (m/s) of disturbances at a latitude: tropical
/// easterly drift, mid-latitude westerly drift, smooth ramps, zero poleward
/// of 75 deg.
pub fn steering_m_s(lat_deg: f64, tropical: f64, midlat: f64) -> f64 {
    let a = lat_deg.abs();
    let ss = |x: f64| { let x = x.clamp(0.0, 1.0); x * x * (3.0 - 2.0 * x) };
    if a <= 20.0 {
        tropical
    } else if a <= 35.0 {
        tropical + (midlat - tropical) * ss((a - 20.0) / 15.0)
    } else if a <= 55.0 {
        midlat
    } else {
        midlat * (1.0 - ss((a - 55.0) / 20.0))
    }
}

pub struct DisturbanceField {
    law: Scope,
    pub planet_r_km: f64,
    pub steer_tropical: f64,
    pub steer_midlat: f64,
    memo: RefCell<BTreeMap<(u8, i64, i64, i64, i64), Rc<Vec<Disturbance>>>>,
}

impl DisturbanceField {
    pub fn new(world_seed: u64, law_id: &str, planet_r_km: f64, steer_tropical: f64, steer_midlat: f64) -> DisturbanceField {
        let law = Scope::root(world_seed).law_set(law_id).expect("law id is printable ASCII");
        DisturbanceField { law, planet_r_km, steer_tropical, steer_midlat, memo: RefCell::new(BTreeMap::new()) }
    }

    pub fn clear_memo(&self) { self.memo.borrow_mut().clear(); }
    pub fn memo_len(&self) -> usize { self.memo.borrow().len() }

    /// The law-set scope (for other keyed draws of the same law).
    pub fn law_scope(&self) -> &Scope { &self.law }

    /// Drift rate (rad/h) of the slot frame at a latitude.
    pub fn omega_rad_h(&self, lat_deg: f64) -> f64 {
        let u = steering_m_s(lat_deg, self.steer_tropical, self.steer_midlat);
        let c = (lat_deg * DEG).cos().max(0.1);
        u * 3600.0 / (self.planet_r_km * 1000.0 * c)
    }

    /// The disturbances keyed to one cube and slot (pure; memoized).
    pub fn generate(&self, p: &ProcessSpec, cube: [i64; 3], slot: i64) -> Rc<Vec<Disturbance>> {
        let key = (p.id, cube[0], cube[1], cube[2], slot);
        if let Some(v) = self.memo.borrow().get(&key) {
            return v.clone();
        }
        let v = Rc::new(self.generate_uncached(p, cube, slot));
        self.memo.borrow_mut().insert(key, v.clone());
        v
    }

    fn generate_uncached(&self, p: &ProcessSpec, cube: [i64; 3], slot: i64) -> Vec<Disturbance> {
        let mut out = Vec::new();
        if slot < 0 { return out; }
        let l = p.lattice_km;
        let r = self.planet_r_km;
        // Cubes that cannot intersect the shell |x| in [R - L/2, R + L/2] are empty.
        let lo = [cube[0] as f64 * l, cube[1] as f64 * l, cube[2] as f64 * l];
        let mut dmin2 = 0.0;
        let mut dmax2 = 0.0;
        for a in 0..3 {
            let (c0, c1) = (lo[a], lo[a] + l);
            let near = if c0 > 0.0 { c0 } else if c1 < 0.0 { -c1 } else { 0.0 };
            let far = c0.abs().max(c1.abs());
            dmin2 += near * near;
            dmax2 += far * far;
        }
        let hs = 0.5 * SHELL_FRAC * l;
        if dmin2.sqrt() > r + hs || dmax2.sqrt() < r - hs { return out; }
        let scope = self.law.address(&format!("{}/{}/{}/{}", p.tag, cube[0], cube[1], cube[2])).expect("address");
        let draw = |what: &str, counter: u64| {
            let purpose = format!("weather.{}.{}", p.tag, what);
            scope.draw(&purpose, counter, 0).expect("purpose")
        };
        let mu = p.mu_per_cube(r);
        // Poisson count by inversion of one uniform.
        let u = draw("count", slot as u64).uniform_f64();
        let (mut k, mut pk) = (0u64, (-mu).exp());
        let mut cdf = pk;
        while u >= cdf && k < 255 {
            k += 1;
            pk *= mu / k as f64;
            cdf += pk;
        }
        for i in 0..k {
            let c = slot as u64 * 256 + i;
            let x = lo[0] + l * draw("cx", c).uniform_f64();
            let y = lo[1] + l * draw("cy", c).uniform_f64();
            let z = lo[2] + l * draw("cz", c).uniform_f64();
            let rr = (x * x + y * y + z * z).sqrt();
            // Marks are per-index keyed draws: skipping them for a discarded
            // centre changes no other draw.
            if (rr - r).abs() >= hs || rr == 0.0 { continue; }
            let t0 = (slot as f64 + draw("start", c).uniform_f64()) * p.slot_h;
            let dur = p.d_min_h + (p.d_max_h - p.d_min_h) * draw("dur", c).uniform_f64();
            let (a, b) = (p.r_min_km, p.r_max_km);
            // Radius uniform in [a, b].
            let radius = a + (b - a) * draw("radius", c).uniform_f64();
            let strength = draw("thin", c).uniform_f64();
            let zn = draw("intensity", c).normal_f64();
            let s = p.intensity_sigma;
            let intensity = (s * zn - 0.5 * s * s).exp();
            out.push(Disturbance { centre: [x / rr, y / rr, z / rr], t0_h: t0, dur_h: dur, radius_km: radius, strength, intensity });
        }
        out
    }

    /// Every disc of a process whose halo covers (lat, lon) at hour t.
    pub fn hits(&self, p: &ProcessSpec, lat_deg: f64, lon_deg: f64, t_h: f64) -> Vec<Hit> {
        let mut out = Vec::new();
        let r = self.planet_r_km;
        let s_hi = (t_h / p.slot_h).floor() as i64;
        let s_lo = ((t_h - p.d_max_h) / p.slot_h).floor() as i64;
        let omega = self.omega_rad_h(lat_deg);
        let north_h = lat_deg >= 0.0;
        for s in s_lo.max(0)..=s_hi {
            let dt = t_h - s as f64 * p.slot_h;
            let lon_b = lon_deg - omega * dt / DEG;
            let x = unit(lat_deg, lon_b);
            let cube = [
                (x[0] * r / p.lattice_km).floor() as i64,
                (x[1] * r / p.lattice_km).floor() as i64,
                (x[2] * r / p.lattice_km).floor() as i64,
            ];
            // Local east/north at the back-advected point (same components as
            // at the true point: the frames differ by a polar-axis rotation).
            let lb = lon_b * DEG;
            let ph = lat_deg * DEG;
            let east = [-lb.sin(), lb.cos(), 0.0];
            let north = [-ph.sin() * lb.cos(), -ph.sin() * lb.sin(), ph.cos()];
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let list = self.generate(p, [cube[0] + dx, cube[1] + dy, cube[2] + dz], s);
                        for d in list.iter() {
                            if t_h < d.t0_h || t_h >= d.t0_h + d.dur_h { continue; }
                            let cr = cross(x, d.centre);
                            let ang = norm(cr).atan2(dot(x, d.centre));
                            let dist = ang * r;
                            if dist >= d.radius_km * p.halo { continue; }
                            // Outward radial direction at x (away from the centre).
                            let toward = [d.centre[0] - dot(x, d.centre) * x[0], d.centre[1] - dot(x, d.centre) * x[1], d.centre[2] - dot(x, d.centre) * x[2]];
                            let tn = norm(toward);
                            let cyc = if tn > 1e-12 {
                                let outward = [-toward[0] / tn, -toward[1] / tn, -toward[2] / tn];
                                // up x outward = counter-clockwise tangent seen from above.
                                let mut t = cross(x, outward);
                                if !north_h { t = [-t[0], -t[1], -t[2]]; }
                                [dot(t, east), dot(t, north)]
                            } else {
                                [0.0, 0.0]
                            };
                            let clat = d.centre[2].clamp(-1.0, 1.0).asin() / DEG;
                            let clon = d.centre[1].atan2(d.centre[0]) / DEG + self.omega_rad_h(clat) * dt / DEG;
                            out.push(Hit {
                                centre_lat_deg: clat,
                                centre_lon_deg: (clon + 180.0).rem_euclid(360.0) - 180.0,
                                d_frac: dist / d.radius_km,
                                tau: (t_h - d.t0_h) / d.dur_h,
                                strength: d.strength,
                                intensity: d.intensity,
                                cyclonic_en: cyc,
                            });
                        }
                    }
                }
            }
        }
        out
    }
}
