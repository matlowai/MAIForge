//! Deterministic runout of one clast (or one sub-parcel of a size class): a lumped sliding / rolling mass on the
//! surface, horizontal state (position, velocity), explicit steps.
//!
//! Per step: gravity drives along the bilinear-interpolated surface gradient (horizontal part of the along-slope
//! acceleration, g tan(theta) cos^2(theta)); Coulomb friction mu g cos^2(theta) opposes the motion and may stop it
//! exactly (never reverses it). On faces steeper than `tan_ff` the mass is falling / bouncing: friction is the weak
//! face-contact `mu_face`. Leaving a face onto gentler ground is an impact: speed x= rt cos(theta_face - theta_ground)
//! (size-dependent rt). A stopped clast stays only where tan(theta) <= mu_static; elsewhere gravity restarts it.
//! Sliding with Coulomb friction from rest is the energy-line (Fahrboeschung / Heim) model; the impact term is the
//! lumped free-fall energy loss at the cliff foot.
//!
//! Everything is a pure function of (start, kinematic params, the gradient arrays): no reductions, no shared
//! state, so a trace is bit-identical in any tile domain that contains its whole path.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stop {
    Rest = 0,
    Sink = 1,
    OffDomain = 2,
    Capped = 3,
}

pub struct Terrain<'a> {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub gx: &'a [f64],
    pub gy: &'a [f64],
    pub sink: &'a [u8],
    /// Ground class per post: 0 rock / debris (base friction), 1 soil-mantled or vegetated (base + soil increment).
    pub ground: &'a [u8],
    /// Rock-surface slope (tan) per post: the face foot is where the clast first reaches rock gentler than tan_foot,
    /// on the same reference for field parcels and blocks.
    pub rock_tan: &'a [f64],
}

#[derive(Clone, Copy)]
pub struct Kin {
    /// Dynamic and static friction (tan) on ground class 0 (rock / debris) and 1 (soil).
    pub mu_d: [f64; 2],
    pub mu_s: [f64; 2],
    pub rt: f64,
    /// Keyed per-clast word for the lateral deflections (from the clast's own vz.keyed-rng draw).
    pub seed: u64,
}

#[derive(Clone, Copy)]
pub struct TraceCfg {
    pub tan_ff: f64,
    pub tan_foot: f64,
    pub mu_face: f64,
    pub g: f64,
    pub l_max_m: f64,
    pub dt_max_s: f64,
    pub step_posts: f64,
    pub max_steps: u32,
    /// Lateral dispersion: at each face -> ground impact the velocity turns by U(-a, a); while rolling on ground
    /// it turns by U(-b, b) every `turn_every_m` of path (block-scale roughness).
    pub lat_impact_rad: f64,
    pub lat_roll_rad: f64,
    pub turn_every_m: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Trace {
    /// Stop position in posts RELATIVE TO THE SOURCE POST (domain independent).
    pub fx: f64,
    pub fy: f64,
    pub stop: Stop,
    pub sink: u8,
    pub path_m: f64,
    /// First position where the rock surface is gentler than `tan_foot` (leaving the source face), posts relative to
    /// the source post.
    pub foot: Option<(f64, f64)>,
    pub impacts: u32,
    pub steps: u32,
}

impl<'a> Terrain<'a> {
    /// Gradient at a position given as an integer post (oi, oj) of this domain plus a float offset (px, py) in posts.
    /// The cell and the bilinear weights come from the offset alone, so the same clast (same source post, same
    /// offset) gets bit-identical weights in every domain that contains it, whatever the domain origin.
    #[inline]
    pub fn grad(&self, oi: i64, oj: i64, px: f64, py: f64) -> Option<(f64, f64)> {
        let (fpx, fpy) = (px.floor(), py.floor());
        let i = oi + fpx as i64;
        let j = oj + fpy as i64;
        if i < 0 || j < 0 || i > self.nx as i64 - 2 || j > self.ny as i64 - 2 { return None; }
        let (i, j) = (i as usize, j as usize);
        let (a, b) = (px - fpx, py - fpy);
        let k = j * self.nx + i;
        let (k1, k2, k3) = (k + 1, k + self.nx, k + self.nx + 1);
        let w0 = (1.0 - a) * (1.0 - b);
        let w1 = a * (1.0 - b);
        let w2 = (1.0 - a) * b;
        let w3 = a * b;
        Some((
            w0 * self.gx[k] + w1 * self.gx[k1] + w2 * self.gx[k2] + w3 * self.gx[k3],
            w0 * self.gy[k] + w1 * self.gy[k1] + w2 * self.gy[k2] + w3 * self.gy[k3],
        ))
    }
    #[inline]
    pub fn near(&self, oi: i64, oj: i64, px: f64, py: f64) -> usize {
        let i = ((oi + px.round() as i64).max(0) as usize).min(self.nx - 1);
        let j = ((oj + py.round() as i64).max(0) as usize).min(self.ny - 1);
        j * self.nx + i
    }
}

/// splitmix64: the per-clast deflection sequence (pure function of the clast's keyed word and the event index).
#[inline]
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[inline]
fn turn(ux: &mut f64, uy: &mut f64, ang: f64) {
    let (sn, cs) = ang.sin_cos();
    let (a, b) = (*ux * cs - *uy * sn, *ux * sn + *uy * cs);
    *ux = a;
    *uy = b;
}

/// Trace a clast released at source post (oi, oj) of the domain with start offset (fx0, fy0) posts.
pub fn trace(t: &Terrain, oi: i64, oj: i64, fx0: f64, fy0: f64, k: &Kin, c: &TraceCfg) -> Trace {
    let mut ev = 0u64;
    let mut next_turn = c.turn_every_m;
    let u01 = |ev: &mut u64| -> f64 { *ev += 1; (mix(k.seed ^ ev.wrapping_mul(0xD1B5_4A32_D192_ED03)) >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0) };
    let (mut fx, mut fy) = (fx0, fy0);
    let (mut ux, mut uy) = (0.0f64, 0.0f64);
    let mut path = 0.0;
    let mut foot = None;
    let mut in_face = false;
    let mut face_theta = 0.0f64;
    let mut impacts = 0u32;
    let mut steps = 0u32;
    let dx = t.dx;
    loop {
        let (gx, gy) = match t.grad(oi, oj, fx, fy) {
            Some(g) => g,
            None => return Trace { fx, fy, stop: Stop::OffDomain, sink: 0, path_m: path, foot, impacts, steps },
        };
        let kn = t.near(oi, oj, fx, fy);
        let s = t.sink[kn];
        let gc = t.ground[kn] as usize;
        if s != 0 && steps > 0 { return Trace { fx, fy, stop: Stop::Sink, sink: s, path_m: path, foot, impacts, steps }; }
        let tan2 = gx * gx + gy * gy;
        let tan = tan2.sqrt();
        let c2 = 1.0 / (1.0 + tan2);
        let face = tan > c.tan_ff;
        if foot.is_none() && t.rock_tan[kn] < c.tan_foot { foot = Some((fx, fy)); }
        if face {
            in_face = true;
            face_theta = tan.atan();
        } else if in_face {
            in_face = false;
            let f = k.rt * (face_theta - tan.atan()).cos().max(0.0);
            ux *= f;
            uy *= f;
            impacts += 1;
            let a = c.lat_impact_rad * (2.0 * u01(&mut ev) - 1.0);
            turn(&mut ux, &mut uy, a);
        }
        let mu = if face { c.mu_face } else { k.mu_d[gc] };
        let sp0 = (ux * ux + uy * uy).sqrt();
        if sp0 == 0.0 && !face && tan <= k.mu_s[gc] {
            return Trace { fx, fy, stop: Stop::Rest, sink: 0, path_m: path, foot, impacts, steps };
        }
        let dt = if sp0 > 0.0 { c.dt_max_s.min(c.step_posts * dx / sp0) } else { c.dt_max_s };
        ux -= c.g * gx * c2 * dt;
        uy -= c.g * gy * c2 * dt;
        let sp = (ux * ux + uy * uy).sqrt();
        let fr = mu * c.g * c2 * dt;
        if sp <= fr {
            ux = 0.0;
            uy = 0.0;
        } else {
            let r = (sp - fr) / sp;
            ux *= r;
            uy *= r;
        }
        let mx = ux * dt;
        let my = uy * dt;
        fx += mx / dx;
        fy += my / dx;
        path += (mx * mx + my * my).sqrt();
        steps += 1;
        if !face && path >= next_turn {
            next_turn = path + c.turn_every_m;
            let a = c.lat_roll_rad * (2.0 * u01(&mut ev) - 1.0);
            turn(&mut ux, &mut uy, a);
        }
        if path > c.l_max_m || steps >= c.max_steps {
            if t.grad(oi, oj, fx, fy).is_none() { return Trace { fx, fy, stop: Stop::OffDomain, sink: 0, path_m: path, foot, impacts, steps }; }
            return Trace { fx, fy, stop: Stop::Capped, sink: 0, path_m: path, foot, impacts, steps };
        }
    }
}
