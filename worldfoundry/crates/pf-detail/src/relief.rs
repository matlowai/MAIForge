//! Lane LOBE (2026-10-16; opt-in `de_lobe_relief=1`, needs de_lobes=1 and de_ve=1): constructional surface relief of the
//! youngest exposed lava flows.
//!
//! design 2026-10-16: lava lobes do not read as flow fronts in the walk (views gw10c / gw10d: a plain bare slope at a lane VE lobe
//! front); "geology-derived generation".
//!
//! The visible flows are lane VE's recorded lava units (window header `prep.net.volcanoes[k].ve.units`, kind Path, not lahar):
//! T1b drapes each on its 32 m lattice (pf-region t1b_ve.rs `unit_drape_m`: nearest 128 m path segment, half-width
//! w0 (1 + 0.5 min(1, s / L)), snout disk of radius w_end (1 + 0.25 lobes(angle)) past the path end, thickness (t - erosion)
//! smoothstep((w - lat) / m) min(1, s / 200 m), m >= 48 m). A 32 m lattice cannot carry a 30-40 deg flow front a few tens of
//! metres wide, levees, or a drained channel, so the bicubic T2 surface shows a smooth bare slope there. This key adds:
//!
//! 1. the part of each flow's constructional profile the T1b lattice cannot carry: Q - B[Q], where Q is the composite profile
//!    of the flows (same footprint as T1b's drape) with a repose-limited margin (fresh slope phi: 35 deg blocky andesite /
//!    dacite fronts, 25 deg basaltic pahoehoe / aa toes) and B[Q] its Catmull-Rom bicubic through the T1b posts (the
//!    representation T1b already carries; zero at every T1b post, mass-neutral at the T1b scale, so tile means and the B8
//!    correction are not disturbed);
//! 2. levees and a drained channel where the record implies channelised flow (record length >= chan_len_w x w0): levee crests
//!    a_l H above the flow top on both margins, the channel floor a_c H below it, volume-neutral across each section (the
//!    levee width is solved for that), tapered in over w0 after the 200 m ramp-in and out before the dispersed-flow toe zone
//!    (toe_w x w_end before the path end); T1b carries none of this, so it is added whole;
//! 3. degradation by exposure age: the cross-section is convolved with a Gaussian of sigma = sqrt(sigma0^2 + 2 kappa age)
//!    (linear hillslope diffusion of a scarp, the morphologic-dating model; age = the record's age_yr; thickness already net
//!    of the record's erosion_m), so young fronts stand at repose and older ones are rounded and lower-angled;
//! 4. burial: a flow's relief is multiplied by prod over younger lane VE units V of (1 - min(1, t_V / H)), t_V = V's T1b drape
//!    thickness (lava and lahar paths, cinder cones, domes, fissures; lahars only beyond 0.8 of the body radius, the T1b
//!    rule), so buried fronts contribute nothing; thin aprons and tephra drape and are ignored.
//!
//! The column agrees: the recorded VE lava units of the lobe column take this same footprint and thickness (the degraded
//! profile with levees; no lobate margins and the record's 128 m polyline: their resolved geometry is the record's) and the
//! lobe construction surface D becomes max(S, T1b + R), so the unit hangs from the relief surface: the flow's own rubbly top
//! / front breccia shows on its top and front, older units beyond its toe.
//!
//! No keyed draws. Parameters are recalled textbook values, unverified (listed in RELIEF_NOTE and the lane report).

use serde::Serialize;

pub const RELIEF_NOTE: &str = "lane LOBE de_lobe_relief=1 (needs de_lobes=1 de_ve=1): constructional surface relief of lane VE's recorded lava flows (Path units, not lahars) on T1b's own drape footprint (nearest 128 m segment, w0 (1 + 0.5 min(1, s / L)), snout disk w_end (1 + 0.25 lobes) past the path end, thickness t - erosion_m, 200 m ramp-in). Added to the T2 surface before the B5 relaxation: (1) mode 1: Q - B[Q] (mode 2: Q - B[D], D = T1b's own drape thickness of the flows at its posts, so T1b's 48 m-floored smoothstep margin is replaced by this profile), Q = the composite profile with a repose-limited margin (the surface of the margin, regional gradient of the routing surface plus the thickness ramp, stands at the fresh repose slope phi_deg relative to the horizontal: blocky andesite / dacite, basaltic; so a downhill snout front has a gentler thickness ramp than a lateral margin), B = Catmull-Rom bicubic through the T1b posts (zero at T1b posts, mass-neutral at the T1b scale); (2) levees a_l H and a drained channel a_c H (volume-neutral per section) where the record length >= chan_len_w x w0, tapered in over w0 after the ramp-in and out before a toe zone of toe_w x w_end (dispersed-flow zone); (3) degradation by exposure age: Gaussian convolution of the section, sigma = sqrt(sigma0^2 + 2 kappa max(0, age_yr - lag_yr)) (linear diffusion of a scarp, starting after lane VE2's infiltration lag of fresh lava, 10 ky andesite / 50 ky basalt; kappa recalled 1e-3 to 1e-2 m2/yr for soil-mantled scarps and cinder cones, the default 2e-3 for coarse blocky margins, unverified); the relief tapers to 0 at a T1a channel corridor over corridor_taper_m (corridors are fixed); in a tile carrying relief the B8 tile-mean correction is weighted by 0.15 + |R| / its tile mean (at most 4), as lane DE2 B1 does for the DE change, so the small per-tile mean of R is taken up where the relief is instead of as a tile-shaped bump; (4) burial by younger lane VE units: x prod(1 - min(1, t_V / H)) with t_V their T1b drape thickness. The recorded VE lava units of the lobe column take the same footprint and thickness (no lobate margins) and the lobe construction surface becomes max(S, T1b + R). Fresh slopes, levee fractions and the toe zone are recalled values (blocky andesite fronts 30-40 deg, levees and channels of a'a / blocky flows, Lipman and Banks 1987 dispersed-flow toe), not calibrated. No keyed draws.";

#[derive(Clone, Debug, Serialize)]
pub struct ReliefParams {
    pub version: u32,
    /// Fresh margin slope (deg): [blocky andesite / dacite, basaltic].
    pub phi_deg: [f64; 2],
    /// Initial rounding sigma (m) of the margin profile: [blocky, basaltic].
    pub sigma0_m: [f64; 2],
    /// Diffusivity (m2 / yr) of the age degradation, acting after the infiltration lag (lane VE2's fresh-lava rule: coherent
    /// lava sheds little until the lag has passed): degradation age = max(0, age - lag_yr[blocky, basaltic]).
    pub kappa_m2_yr: f64,
    pub lag_yr: [f64; 2],
    /// Relief taper towards a T1a channel corridor (fixed during relaxation): x smoothstep(excess / corridor_taper_m).
    pub corridor_taper_m: f64,
    /// [levee crest a_l, channel depth a_c] as fractions of the local flow thickness: [blocky, basaltic].
    pub levee: [[f64; 2]; 2],
    /// Channelised when the record length >= chan_len_w x w0; the channel ends toe_w x w_end before the path end.
    pub chan_len_w: f64,
    pub toe_w: f64,
    /// 1: add Q - B[Q] (the part of the flows' profile the T1b lattice cannot carry; T1b posts unchanged); 2: replace T1b's own
    /// drape of the flows, Q - B[D] (D = T1b's drape thickness at its posts), so the margin is the repose-limited profile
    /// itself (T1b posts change by Q - D; the volume difference is taken up by the |R|-weighted B8 correction).
    pub mode: u8,
    pub note: String,
}

impl Default for ReliefParams {
    fn default() -> Self {
        ReliefParams { version: 1, phi_deg: [35.0, 25.0], sigma0_m: [1.0, 3.0], kappa_m2_yr: 2.0e-3, lag_yr: [10_000.0, 50_000.0], corridor_taper_m: 24.0, levee: [[0.25, 0.3], [0.15, 0.15]],
                       chan_len_w: 8.0, toe_w: 3.0, mode: 1, note: RELIEF_NOTE.into() }
    }
}

/// T1b's composite-path margin floor (lane VE3 aa: 1.5 x 32 m).
const T1B_MFLOOR: f64 = 48.0;

#[derive(Clone, Debug)]
pub enum RGeom {
    /// Polyline [x, y, s] (the record's 128 m vertices), record length, w0, thickness net of erosion, snout lobe phases.
    Path { pts: Vec<[f64; 3]>, len: f64, w0: f64, lobes: [f64; 4] },
    Shape(crate::lobes::VeShape),
}

/// One lane VE unit as the relief sees it.
#[derive(Clone, Debug)]
pub struct RUnit {
    pub age_yr: f64,
    /// Thickness net of erosion (paths); unused for shapes (their heights are already net).
    pub h: f64,
    pub lava: bool,
    pub lahar: bool,
    pub basalt: bool,
    pub geom: RGeom,
    pub bbox: (f64, f64, f64, f64),
    /// Lahars drape only beyond 0.8 of the body radius from the vent (T1b deposit rule): (vent, r_body).
    pub vent: [f64; 2],
    pub r_body: f64,
    /// Regional (pre-flow) surface gradient at the path vertices (d z / d x, d z / d y): the box-mean routing surface's
    /// gradient, so the margin stands at repose relative to the horizontal (empty: level ground).
    pub bgrad: Vec<[f64; 2]>,
}

#[inline]
fn smoothstep01(u: f64) -> f64 { let u = u.clamp(0.0, 1.0); u * u * (3.0 - 2.0 * u) }

/// pf-region t1b_ve.rs `lobe` (snout lobes of a path end, unit amplitude).
#[inline]
fn lobe(l: &[f64; 4], th: f64) -> f64 { ((2.0 * th + l[0]).cos() / 2.0 + (3.0 * th + l[1]).cos() / 3.0 + (4.0 * th + l[2]).cos() / 4.0 + (5.0 * th + l[3]).cos() / 5.0) / 1.2833 }

/// erf (Abramowitz and Stegun 7.1.26, |error| < 1.5e-7).
fn erf(x: f64) -> f64 {
    let s = x.signum(); let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * (-x * x).exp();
    s * y
}

/// Gaussian-smoothed ramp max(x, 0): E[max(x + sigma Z, 0)] = x Phi(x / sigma) + sigma phi(x / sigma).
#[inline]
fn sramp(x: f64, sg: f64) -> f64 {
    if sg <= 1e-9 { return x.max(0.0); }
    let z = x / sg;
    if z > 8.0 { return x; }
    if z < -8.0 { return 0.0; }
    x * 0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2)) + sg * (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// Path geometry at (x, y) exactly as T1b's drape: (signed distance inside the margin d = w_eff - lateral distance, along-path
/// s, past the end (snout), local half-width w_eff). None far from the path.
pub fn path_geom(pts: &[[f64; 3]], len: f64, w0: f64, lobes: &[f64; 4], x: f64, y: f64, reach: f64) -> Option<(f64, f64, bool, f64)> {
    path_geom_n(pts, len, w0, lobes, x, y, reach).map(|g| (g.0, g.1, g.2, g.3))
}

/// As `path_geom`, plus the outward unit normal of the margin at (x, y) (from the nearest axis point, or the end point in the
/// snout) and the nearest segment's vertex weights (segment k, fraction t).
#[allow(clippy::type_complexity)]
pub fn path_geom_n(pts: &[[f64; 3]], len: f64, w0: f64, lobes: &[f64; 4], x: f64, y: f64, reach: f64) -> Option<(f64, f64, bool, f64, [f64; 2], usize, f64)> {
    if pts.len() < 2 { return None; }
    let wmax = w0 * 1.5 * 1.3 + reach;
    let n = pts.len() - 1;
    let (mut best, mut bs, mut bb, mut ba) = (f64::INFINITY, 0.0, 0.0, 0.0);
    let (mut bn, mut bk, mut bt) = ([0.0f64; 2], 0usize, 0.0f64);
    for k in 0..n {
        let (p, q) = (pts[k], pts[k + 1]);
        if x < p[0].min(q[0]) - wmax || x > p[0].max(q[0]) + wmax || y < p[1].min(q[1]) - wmax || y > p[1].max(q[1]) + wmax { continue; }
        let (dx, dy) = (q[0] - p[0], q[1] - p[1]);
        let l2 = dx * dx + dy * dy;
        if l2 <= 0.0 { continue; }
        let ln = l2.sqrt();
        let tt = ((x - p[0]) * dx + (y - p[1]) * dy) / l2;
        let tc = tt.clamp(0.0, 1.0);
        let (ox, oy) = if k == n - 1 && tt > 1.0 { (x - q[0], y - q[1]) } else { (x - p[0] - tc * dx, y - p[1] - tc * dy) };
        let (d, s, be, an) = if k == n - 1 && tt > 1.0 { let (ex, ey) = (x - q[0], y - q[1]); (ex.hypot(ey), q[2], (ex * dx + ey * dy) / ln, ey.atan2(ex) - dy.atan2(dx)) }
            else { ((x - p[0] - tc * dx).hypot(y - p[1] - tc * dy), p[2] + tc * (q[2] - p[2]), 0.0, 0.0) };
        if d < best { best = d; bs = s; bb = be; ba = an; let l = ox.hypot(oy); bn = if l > 1e-9 { [ox / l, oy / l] } else { [0.0, 0.0] }; bk = k; bt = tc; }
    }
    if !best.is_finite() { return None; }
    let w = w0 * (1.0 + 0.5 * (bs / len.max(1.0)).min(1.0));
    let weff = if bb > 0.0 { w * (1.0 + 0.25 * lobe(lobes, ba)) } else { w };
    Some((weff - best, bs, bb > 0.0, weff, bn, bk, bt))
}

impl RUnit {
    /// T1b's own drape thickness of this unit (burial of older flows).
    pub fn t1b_thick(&self, x: f64, y: f64) -> f64 {
        if x < self.bbox.0 || x > self.bbox.2 || y < self.bbox.1 || y > self.bbox.3 { return 0.0; }
        match &self.geom {
            RGeom::Shape(sh) => sh.thick(x, y),
            RGeom::Path { pts, len, w0, lobes } => {
                if self.lahar && (x - self.vent[0]).hypot(y - self.vent[1]) <= 0.8 * self.r_body { return 0.0; }
                let Some((d, s, _, weff)) = path_geom(pts, *len, *w0, lobes, x, y, 64.0) else { return 0.0 };
                let m = (self.h / 0.7).max(T1B_MFLOOR).min(weff.max(T1B_MFLOOR));
                self.h * smoothstep01(d / m) * (s / 200.0).min(1.0)
            }
        }
    }
}

/// Degradation sigma (m) of a flow's section: sqrt(sigma0^2 + 2 kappa max(0, age - lag)).
pub fn sigma_of(rp: &ReliefParams, u: &RUnit) -> f64 {
    let k = u.basalt as usize;
    (rp.sigma0_m[k].powi(2) + 2.0 * rp.kappa_m2_yr * (u.age_yr - rp.lag_yr[k]).max(0.0)).sqrt()
}

/// Section of a lava flow at one point: (plain degraded profile, levee / channel modification, both m).
pub fn section(rp: &ReliefParams, u: &RUnit, x: f64, y: f64) -> Option<(f64, f64)> {
    let RGeom::Path { pts, len, w0, lobes } = &u.geom else { return None };
    if !u.lava { return None; }
    let k = u.basalt as usize;
    let sg = sigma_of(rp, u);
    if x < u.bbox.0 - 4.0 * sg || x > u.bbox.2 + 4.0 * sg || y < u.bbox.1 - 4.0 * sg || y > u.bbox.3 + 4.0 * sg { return None; }
    let (d, s, snout, weff, nrm, sk, st) = path_geom_n(pts, *len, *w0, lobes, x, y, 64.0 + 4.0 * sg)?;
    if d < -4.0 * sg { return None; }
    let h = u.h * (s / 200.0).min(1.0);
    if h <= 0.0 { return None; }
    // margin thickness gradient sl such that the surface (regional gradient b plus the thickness ramp, -sl along the outward
    // normal n) stands at the repose slope: |b - sl n| = tan(phi), sl = b.n + sqrt((b.n)^2 - |b|^2 + tan^2(phi)), clamped to
    // [0.25, 1.5] tan(phi) (a downhill snout front is lower-angled in thickness, a lateral margin slightly so)
    let tp = rp.phi_deg[k].to_radians().tan();
    let b = if u.bgrad.len() == pts.len() { let (a, c) = (u.bgrad[sk], u.bgrad[(sk + 1).min(pts.len() - 1)]); let t = if snout { 1.0 } else { st }; [a[0] + t * (c[0] - a[0]), a[1] + t * (c[1] - a[1])] } else { [0.0, 0.0] };
    let bn = b[0] * nrm[0] + b[1] * nrm[1];
    let disc = bn * bn - (b[0] * b[0] + b[1] * b[1]) + tp * tp;
    let sl = (if disc > 0.0 { bn + disc.sqrt() } else { 0.25 * tp }).clamp(0.25 * tp, 1.5 * tp);
    let plain = sl * (sramp(d, sg) - sramp(d - h / sl, sg));
    // levees and a drained channel in the channelised reach (volume-neutral per section)
    let [al, ac] = rp.levee[k];
    let mut dl = 0.0;
    if !snout && *len >= rp.chan_len_w * w0 && al + ac > 0.0 {
        let w_end = 1.5 * w0;
        let c = smoothstep01((s - 200.0) / w0) * smoothstep01(((len - rp.toe_w * w_end) - s) / w0);
        if c > 0.0 {
            let hs = h / sl;
            let cw = (ac * weff - (ac * (1.0 + 2.0 * al + 0.5 * ac) + al * al) * hs) / (al + ac);
            let d1 = (1.0 + al) * hs; let d2 = d1 + cw; let d3 = d2 + (al + ac) * hs;
            if cw > 0.0 && d3 < weff {
                let f = sl * (sramp(d, sg) - sramp(d - d1, sg) - sramp(d - d2, sg) + sramp(d - d3, sg));
                dl = c * (f - plain);
            }
        }
    }
    Some((plain, dl))
}

/// The relief world of a block: VE units of the used volcanoes and the plain-profile lattice at T1b posts over the block.
pub struct ReliefWorld {
    pub p: ReliefParams,
    pub units: Vec<RUnit>,
    /// Lattice of Q at T1b posts: global post origin, size, values.
    gx0: i64, gy0: i64, nx: usize, ny: usize, q: Vec<f64>,
}

impl ReliefWorld {
    /// Composite (plain, levee) profile at (x, y): sum over lava flows of burial weight x (plain, levee modification).
    pub fn composite(p: &ReliefParams, units: &[RUnit], x: f64, y: f64) -> (f64, f64) { let q = Self::composite3(p, units, x, y); (q.0, q.1) }
    /// (plain, levee, T1b drape) composites: the third is the sum of burial weight x T1b's own drape thickness of each flow.
    pub fn composite3(p: &ReliefParams, units: &[RUnit], x: f64, y: f64) -> (f64, f64, f64) {
        let (mut qb, mut ql, mut qd) = (0.0, 0.0, 0.0);
        for (i, u) in units.iter().enumerate() {
            let Some((pl, dl)) = section(p, u, x, y) else { continue };
            let dt = if p.mode == 2 { u.t1b_thick(x, y) } else { 0.0 };
            if pl.abs() < 1e-6 && dl.abs() < 1e-6 && dt <= 0.0 { continue; }
            let hh = u.h.max(0.5);
            let mut w = 1.0;
            for (j, v) in units.iter().enumerate() {
                if j == i || v.age_yr >= u.age_yr { continue; }
                let tv = v.t1b_thick(x, y);
                if tv > 0.0 { w *= 1.0 - (tv / hh).min(1.0); if w <= 0.0 { break; } }
            }
            qb += w * pl; ql += w * dl; qd += w * dt;
        }
        (qb, ql, qd)
    }
    pub fn build(p: &ReliefParams, units: Vec<RUnit>, bbox: (f64, f64, f64, f64)) -> ReliefWorld {
        let t = crate::T1B_M;
        let (gx0, gy0) = ((bbox.0 / t).floor() as i64 - 3, (bbox.1 / t).floor() as i64 - 3);
        let (gx1, gy1) = ((bbox.2 / t).ceil() as i64 + 3, (bbox.3 / t).ceil() as i64 + 3);
        let (nx, ny) = ((gx1 - gx0 + 1) as usize, (gy1 - gy0 + 1) as usize);
        let mut q = vec![0.0f64; nx * ny];
        for j in 0..ny { for i in 0..nx {
            let c = Self::composite3(p, &units, (gx0 + i as i64) as f64 * t, (gy0 + j as i64) as f64 * t);
            q[j * nx + i] = if p.mode == 2 { c.2 } else { c.0 };
        } }
        ReliefWorld { p: p.clone(), units, gx0, gy0, nx, ny, q }
    }
    fn q_post(&self, gx: i64, gy: i64) -> f64 {
        let (i, j) = (gx - self.gx0, gy - self.gy0);
        if i >= 0 && j >= 0 && (i as usize) < self.nx && (j as usize) < self.ny { return self.q[j as usize * self.nx + i as usize]; }
        let c = Self::composite3(&self.p, &self.units, gx as f64 * crate::T1B_M, gy as f64 * crate::T1B_M);
        if self.p.mode == 2 { c.2 } else { c.0 }
    }
    /// The relief R at chart (x, y): Q - B[Q] + levees.
    pub fn at(&self, x: f64, y: f64) -> f64 {
        if self.units.is_empty() { return 0.0; }
        let (qb, ql) = Self::composite(&self.p, &self.units, x, y);
        let fx = x / crate::T1B_M; let fy = y / crate::T1B_M;
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (wx, wy) = (cr(fx - ix as f64), cr(fy - iy as f64));
        let mut b = 0.0;
        for (j, wyj) in wy.iter().enumerate() {
            let mut row = 0.0;
            for (i, wxi) in wx.iter().enumerate() { row += wxi * self.q_post(ix - 1 + i as i64, iy - 1 + j as i64); }
            b += wyj * row;
        }
        qb - b + ql
    }
    /// The levee / channel part of the composite at (x, y) (diagnostics and tests).
    pub fn composite_levee(&self, x: f64, y: f64) -> f64 { Self::composite(&self.p, &self.units, x, y).1 }
    /// Column thickness of lava flow `i` at (x, y): its degraded profile with levees (not burial-weighted: buried flows keep
    /// their thickness in the column).
    pub fn unit_thick(&self, i: usize, x: f64, y: f64) -> f64 {
        match self.units.get(i).and_then(|u| section(&self.p, u, x, y)) { Some((pl, dl)) => (pl + dl).max(0.0), None => 0.0 }
    }
}

/// Catmull-Rom weights (as input.rs `cr`).
fn cr(t: f64) -> [f64; 4] {
    let (t2, t3) = (t * t, t * t * t);
    [(-t3 + 2.0 * t2 - t) * 0.5, (3.0 * t3 - 5.0 * t2 + 2.0) * 0.5, (-3.0 * t3 + 4.0 * t2 + t) * 0.5, (t3 - t2) * 0.5]
}

/// Lane VE units of one volcano as relief units (all with a geometry the relief reads: Path lava / lahar, Cone, Dome, Fissure).
/// Returns the units and, per listed unit index, its relief index (Path lava units only).
pub fn ve_relief_units(units: &serde_json::Value, family_composite: bool, vent: [f64; 2], r_body: f64) -> (Vec<RUnit>, Vec<Option<usize>>) {
    let a = units.as_array().cloned().unwrap_or_default();
    let (mut out, mut idx) = (Vec::new(), vec![None; a.len()]);
    for (i, u) in a.iter().enumerate() {
        let f = |k: &str| u.get(k).and_then(|v| v.as_f64()).unwrap_or(0.0);
        let Some(g) = u.get("geom") else { continue };
        let gf = |k: &str| g.get(k).and_then(|v| v.as_f64()).unwrap_or(0.0);
        let ero = f("erosion_m");
        let kind = u.get("kind").and_then(|v| v.as_str()).unwrap_or("");
        let basalt = !family_composite || matches!(kind, "cinder_cone" | "cone_flow");
        let geom = match g.get("g").and_then(|v| v.as_str()).unwrap_or("") {
            "Path" => {
                let pts: Vec<[f64; 3]> = g.get("pts").and_then(|v| v.as_array()).map(|p| p.iter().filter_map(|q| { let q = q.as_array()?; Some([q.first()?.as_f64()?, q.get(1)?.as_f64()?, q.get(2)?.as_f64()?]) }).collect()).unwrap_or_default();
                if pts.len() < 2 { continue; }
                let mut lb = [0.0f64; 4];
                if let Some(l) = g.get("lobes").and_then(|v| v.as_array()) { for (k, q) in l.iter().take(4).enumerate() { lb[k] = q.as_f64().unwrap_or(0.0); } }
                RGeom::Path { pts, len: gf("len_m").max(1.0), w0: gf("w0_m"), lobes: lb }
            }
            "Cone" => RGeom::Shape(crate::lobes::VeShape::Cone { x: gf("x"), y: gf("y"), rb: gf("rb_m"), h: gf("h_m") - ero, rc: gf("rc_m"), dc: gf("dc_m") - ero }),
            "Dome" => RGeom::Shape(crate::lobes::VeShape::Dome { x: gf("x"), y: gf("y"), r: gf("r_m"), h: gf("h_m") - ero }),
            "Fissure" => RGeom::Shape(crate::lobes::VeShape::Fissure { x0: gf("x0"), y0: gf("y0"), x1: gf("x1"), y1: gf("y1"), h: gf("h_m") - ero, w: gf("w_m") }),
            _ => continue,
        };
        let lahar = g.get("lahar").and_then(|v| v.as_bool()).unwrap_or(false);
        let (h, bbox, lava) = match &geom {
            RGeom::Path { pts, w0, .. } => {
                let h = gf("t_m") - ero;
                if h <= 0.05 { continue; }
                let pad = w0 * 1.5 * 1.3 + 64.0;
                let (mut b0, mut b1, mut b2, mut b3) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
                for p in pts { b0 = b0.min(p[0]); b1 = b1.min(p[1]); b2 = b2.max(p[0]); b3 = b3.max(p[1]); }
                (h, (b0 - pad, b1 - pad, b2 + pad, b3 + pad), !lahar)
            }
            RGeom::Shape(sh) => (0.0, sh.bbox(), false),
        };
        if lava { idx[i] = Some(out.len()); }
        out.push(RUnit { age_yr: f("age_yr"), h, lava, lahar, basalt, geom, bbox, vent, r_body, bgrad: Vec::new() });
    }
    (out, idx)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn flow(age: f64, len: f64, w0: f64, t: f64, y0: f64) -> RUnit {
        let n = (len / 128.0).ceil() as usize;
        let pts: Vec<[f64; 3]> = (0..=n).map(|i| { let s = (i as f64 * 128.0).min(len); [s, y0, s] }).collect();
        let pad = w0 * 1.5 * 1.3 + 64.0;
        RUnit { age_yr: age, h: t, lava: true, lahar: false, basalt: false, geom: RGeom::Path { pts, len, w0, lobes: [0.0; 4] },
                bbox: (-pad, y0 - pad, len + pad, y0 + pad), vent: [0.0, 0.0], r_body: 1.0e4, bgrad: Vec::new() }
    }
    #[test]
    fn section_is_repose_limited_degrades_with_age_and_levees_are_volume_neutral() {
        let p = ReliefParams::default();
        let tan = p.phi_deg[0].to_radians().tan();
        // across the channelised reach (s = 2000 m) of a 5 km, w0 150 m, 30 m flow: max slope, levee crest, channel floor
        let prof = |u: &RUnit| -> Vec<(f64, f64)> { (0..=300).map(|k| { let y = -300.0 + 2.0 * k as f64; section(&p, u, 2000.0, y).unwrap_or((0.0, 0.0)) }).collect() };
        let maxslope = |v: &[(f64, f64)], lev: bool| v.windows(2).map(|w| ((w[1].0 + if lev { w[1].1 } else { 0.0 }) - (w[0].0 + if lev { w[0].1 } else { 0.0 })).abs() / 2.0).fold(0.0f64, f64::max);
        let young = prof(&flow(0.0, 5000.0, 150.0, 30.0, 0.0));
        let old = prof(&flow(100_000.0, 5000.0, 150.0, 30.0, 0.0));
        assert!(maxslope(&young, false) <= tan + 1e-3 && maxslope(&young, false) > 0.9 * tan, "young plain slope {}", maxslope(&young, false));
        assert!(maxslope(&old, false) < 0.8 * maxslope(&young, false), "old {} young {}", maxslope(&old, false), maxslope(&young, false));
        // the plain top reaches the thickness; levees above it and the channel floor below it on the centreline
        let w = 150.0 * (1.0 + 0.5 * 2000.0 / 5000.0);
        let mid = young[150];
        assert!((mid.0 - 30.0).abs() < 0.1, "plain top {}", mid.0);
        assert!(mid.1 < -0.2 * 30.0, "channel {}", mid.1);
        assert!(young.iter().any(|q| q.1 > 0.15 * 30.0), "levee crest");
        // volume-neutral across the section (2 m steps) within 2 %% of the section area
        let area: f64 = young.iter().map(|q| q.0 * 2.0).sum();
        let lev: f64 = young.iter().map(|q| q.1 * 2.0).sum();
        assert!(lev.abs() < 0.02 * area, "levee volume {lev} vs section {area}");
        // the levee / channel only in the channelised reach: none in the toe zone and the snout, none for a short flow
        assert!(section(&p, &flow(0.0, 5000.0, 150.0, 30.0, 0.0), 4900.0, 0.0).unwrap().1.abs() < 1e-9);
        assert!(section(&p, &flow(0.0, 5000.0, 150.0, 30.0, 0.0), 5050.0, 0.0).unwrap().1.abs() < 1e-9);
        assert!(section(&p, &flow(0.0, 800.0, 150.0, 30.0, 0.0), 400.0, 0.0).unwrap().1.abs() < 1e-9);
        // footprint: the margin toe is T1b's footprint edge (half-width w)
        let edge = section(&p, &flow(0.0, 5000.0, 150.0, 30.0, 0.0), 2000.0, w + 6.0).map(|q| q.0).unwrap_or(0.0);
        assert!(edge < 0.05, "outside the footprint {edge}");
        // basaltic: lower fresh margin slope
        let mut b = flow(0.0, 5000.0, 150.0, 8.0, 0.0); b.basalt = true;
        let vb = prof(&b);
        assert!(maxslope(&vb, false) <= p.phi_deg[1].to_radians().tan() + 1e-3);
    }
    #[test]
    fn relief_key_is_opt_in_and_needs_ve_units() {
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap(); p.set("de", "1").unwrap(); p.set("de_lobes", "1").unwrap();
        assert!(p.set("de_lobe_relief", "1").is_err(), "needs de_ve");
        assert!(p.set("de_lobe_relief_kappa", "0.01").is_err(), "needs de_lobe_relief");
        p.set("de_ve", "1").unwrap();
        assert!(!serde_json::to_string(p.de.as_ref().unwrap().lobes.as_ref().unwrap()).unwrap().contains("relief"));
        p.set("de_lobe_relief", "1").unwrap(); p.set("de_lobe_relief_kappa", "0.01").unwrap();
        assert_eq!(p.de.as_ref().unwrap().lobes.as_ref().unwrap().relief.as_ref().unwrap().kappa_m2_yr, 0.01);
        assert!(serde_json::to_string(p.de.as_ref().unwrap().lobes.as_ref().unwrap()).unwrap().contains("de_lobe_relief"));
    }
    #[test]
    fn relief_is_zero_at_t1b_posts_and_buried_flows_contribute_nothing() {
        let p = ReliefParams::default();
        // an older flow and a younger flow crossing its middle reach
        let older = flow(30_000.0, 4000.0, 120.0, 25.0, 0.0);
        let mut younger = flow(2_000.0, 3000.0, 150.0, 30.0, 0.0);
        if let RGeom::Path { pts, .. } = &mut younger.geom { for q in pts.iter_mut() { let (x, s) = (q[0], q[2]); q[0] = 2000.0; q[1] = -1500.0 + x; q[2] = s; } }
        younger.bbox = (2000.0 - 400.0, -1500.0 - 400.0, 2000.0 + 400.0, 1500.0 + 400.0);
        let rw = ReliefWorld::build(&p, vec![older.clone(), younger.clone()], (-500.0, -2000.0, 5000.0, 2000.0));
        let only_old = ReliefWorld::build(&p, vec![older.clone()], (-500.0, -2000.0, 5000.0, 2000.0));
        let (mut n, mut nz) = (0, 0);
        for gy in -60..=60 { for gx in -10..=150 {
            let (x, y) = (gx as f64 * 32.0, gy as f64 * 32.0);
            let r = rw.at(x, y) - rw.composite_levee(x, y);
            assert!(r.abs() < 1e-9, "post ({x}, {y}): {r}");
            n += 1;
        } }
        // between posts: the older flow's margin is sharpened away from the younger flow, not under it
        for k in 0..400 {
            let y = 50.0 + 0.5 * k as f64;
            let a = only_old.at(1000.0, y) - only_old.composite_levee(1000.0, y);
            if a.abs() > 0.5 { nz += 1; }
            // under the younger flow (x = 2000): the older flow's relief is gone
            let q = section(&p, &older, 2000.0, y).unwrap_or((0.0, 0.0));
            let (qb, ql) = ReliefWorld::composite(&p, &[older.clone(), younger.clone()], 2000.0, y);
            let (yb, yl) = section(&p, &younger, 2000.0, y).unwrap_or((0.0, 0.0));
            if younger.t1b_thick(2000.0, y) >= older.h { assert!((qb - yb).abs() < 1e-9 && (ql - yl).abs() < 1e-9, "buried at y {y}: {q:?}"); }
        }
        assert!(n > 1000 && nz > 10, "posts {n} sharpened {nz}");
    }
}
