//! Lane RILL (2026-10-16, opt-in `de_b5_chain=1|2`, needs de=1): drainage organisation across B5 tile domains.
//!
//! design: procedural geology that reads naturally in a walkable map; "geology-derived generation". Lanes RIB/RIB2 measured that
//! the B5 gully relaxation (relax.rs) runs per tile domain (386 m) with T1b's upslope drainage injected on the domain ring,
//! and T1b's flow lines on the smooth flank are parallel 32-45 m apart (the 32 m lattice), so the rills inside each domain stay
//! parallel and evenly spaced: the upslope rill network never reaches the domain edge, so rills cannot merge downslope.
//!
//! Rule (no new keyed draws): one COARSE B5 pass over the whole T1b window on a lattice of `cell_posts` x 2 m (8 m default),
//! the same physics as B5 (MFD area, implicit detachment-limited stream power on the D8 receiver, nonlinear creep, the same
//! channel-head law, lithology and runoff factors of the same posts), started from the same initial surface the tile domains
//! use, sampled at the coarse posts (T1a cut + B4 closure + lobe relief). Lithology, bed K and runoff are evaluated once on
//! that initial surface (the coarse pass does not re-read them every step: a declared simplification). SimTime is the same as
//! B5: `steps / dt` coarse steps with stream-power and creep coefficients x dt (dt default 4; creep stays explicitly stable at
//! 8 m). The pass is a pure function of the window and the params; in chain mode the context's bbox is the whole window, so
//! every block sees the same coarse field.
//! Mode 1: each tile domain's ring inflow is the coarse network's MFD discharge crossing into the domain's inner box (each
//! coarse post outside the box delivers its routed area share towards each neighbour inside the box at the segment's entry
//! point, split between the two nearest ring-1 posts). Replaces the T1b ring rule (and de_ringflow).
//! Mode 2: mode 1 plus the coarse pass's change (relaxed - initial, bilinear between coarse posts) added to the domain's
//! initial surface on every post that is not sea / channel corridor, the domain ring included (the ring is a domain edge, not
//! a fixed boundary: holding it at the pre-incision surface made a step the relaxation ramped inward, measured as tile-edge
//! lines on G's coastal corner with binary c1), so the fine B5 refines the coarse valleys; it counts as relaxation change for W6.
//! Mode 3: as mode 2, but the coarse change is first mass-balanced on the coarse lattice by the W6 rule (its masked Gaussian
//! low-pass, sigma massbal_sigma_m, removed window-wide) and enters the domain's initial surface before W6's reference, so
//! the domain W6 does not re-balance it (mode 2's domain-bounded re-balancing of a deep coarse incision differs between
//! neighbouring domains: measured tile-edge lines on G's coastal corner).
//! Seams: the coarse field is continuous and global, so neighbouring domains receive the same channels at the same places.

use serde::{Deserialize, Serialize};

pub const CHAIN_NOTE: &str = "lane RILL de_b5_chain: coarse window-wide B5 pass (same physics, lithology/runoff frozen at the initial surface) conditions each tile domain's ring inflow (1) and initial surface (2); no keyed draws";

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ChainParams {
    pub mode: u32,
    /// Coarse lattice spacing in 2 m posts.
    pub cell_posts: usize,
    /// Coarse time-step multiplier (steps / dt coarse steps, coefficients x dt).
    pub dt: usize,
    pub note: String,
    /// Sensitivity (de_b5_chain_loss=0): the coarse pass keeps the runoff factor f but drops the transmission loss t (water
    /// that enters the coarse network is not lost into permeable units on the way). None = loss on (the de_runoff rule).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub no_loss: Option<bool>,
    /// Modes 2/3, opt-in de_b5_chain_snap=1: each ring-1 inflow moves to the lowest ring-1 post of its side within one coarse
    /// cell on the conditioned surface (`snap_to_floor`). None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub snap: Option<bool>,
    /// Opt-in de_b5_chain_ring=1: the B5 domain ring is a natural (zero-flux) boundary instead of fixed. None = fixed.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub natural_ring: Option<bool>,
    /// Opt-in de_b5_chain_b8=1: the B8 tile-mean correction applies its edge taper AFTER the apron-cap water-filling
    /// (weight = min(phi without taper, cap) x taper), so a tile with a large deficit gets a raise that ramps to 0 over
    /// mean_taper_m at its edges instead of a plateau ending in a one-post step at the shared edge (measured: the tile-edge
    /// creases of modes 2/3 on G; mean=0 removed them). None = the B8 rule unchanged.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub b8_taper_after_cap: Option<bool>,
    /// Lane RILL2, opt-in de_b5_chain_interp=3|4 (modes 2/3): the coarse change is reconstructed between coarse posts by a
    /// cubic kernel instead of bilinear (bilinear has slope kinks on every coarse lattice line = an 8 m lattice imprint in
    /// curvature): 3 = Catmull-Rom cubic convolution (interpolating, C1), 4 = uniform cubic B-spline (approximating, C2,
    /// non-negative weights: no overshoot next to fixed sea posts). Both kernels are partitions of unity (the mean of the
    /// mass-balanced change is kept). None = bilinear.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub interp: Option<u32>,
    /// Lane RILL2, opt-in de_b5_chain_head=1: the coarse pass's channel-head threshold is the B5 one rescaled to the coarse
    /// contour width, A_1 x dx / 2 m (also the constant head area). Reason: on unchannelled ground a lattice post's MFD area is
    /// its specific catchment area (area per unit contour width) times the post spacing, so at dx = 8 m a hillslope post
    /// carries 4x the area of a 2 m post on the same slope and the head law A S^2 >= A_1 (set at 2 m posts) was met at a
    /// quarter of the hillslope length: channels started too close to the divides (measured with de_b5_chain_dump on G's NW
    /// coastal strip: median steep-bluff post A S^2 = 287 m2 against A_1 = 300 m2, a 2 m post there ~72 m2). The head law
    /// on specific catchment area is lattice-independent; channelised (convergent) flow carries the same area at any spacing,
    /// so stream power is unchanged. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub head_scale: Option<bool>,
    /// Diagnostic (de_b5_chain_dump=<dir>): write the coarse fields (f32 / u8 rasters + meta json) to <dir>. Not part of the
    /// params (never serialized): outputs are the same with or without it.
    #[serde(skip)]
    pub dump: Option<String>,
}

impl ChainParams {
    pub fn new(mode: u32) -> Self { ChainParams { mode, cell_posts: 4, dt: 4, note: CHAIN_NOTE.to_string(), no_loss: None, snap: None, natural_ring: None, b8_taper_after_cap: None, interp: None, head_scale: None, dump: None } }
}

/// The coarse pass result: lattice of `n` x `n` posts, spacing `st` fine posts, origin at fine post (fx0, fy0).
pub struct Coarse {
    pub fx0: i64,
    pub fy0: i64,
    pub st: i64,
    pub n: usize,
    pub z: Vec<f64>,
    pub dz: Vec<f64>,
    pub fixed: Vec<bool>,
    /// Routed MFD area on the relaxed coarse surface (m2), already x transmission t for passing on.
    pub pass: Vec<f64>,
    pub seconds: f64,
    /// Reconstruction kernel of `dz_at`: 0/1 bilinear, 3 Catmull-Rom, 4 cubic B-spline (`ChainParams::interp`).
    pub interp: u32,
}

const D8: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];

impl Coarse {
    /// Bilinear coarse change at fine post (px, py) (0 outside the lattice).
    pub fn dz_at(&self, px: i64, py: i64) -> f64 {
        let (fx, fy) = ((px - self.fx0) as f64 / self.st as f64, (py - self.fy0) as f64 / self.st as f64);
        let (i, j) = (fx.floor() as i64, fy.floor() as i64);
        if i < 0 || j < 0 || i + 1 >= self.n as i64 || j + 1 >= self.n as i64 { return 0.0; }
        let (tx, ty) = (fx - i as f64, fy - j as f64);
        let g = |a: i64, b: i64| self.dz[(b as usize) * self.n + a as usize];
        if self.interp == 3 || self.interp == 4 {
            let w = |t: f64| -> [f64; 4] {
                let (t2, t3) = (t * t, t * t * t);
                if self.interp == 3 { [0.5 * (-t3 + 2.0 * t2 - t), 0.5 * (3.0 * t3 - 5.0 * t2 + 2.0), 0.5 * (-3.0 * t3 + 4.0 * t2 + t), 0.5 * (t3 - t2)] }
                else { let u = 1.0 - t; [u * u * u / 6.0, (3.0 * t3 - 6.0 * t2 + 4.0) / 6.0, (-3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0) / 6.0, t3 / 6.0] }
            };
            let (wx, wy) = (w(tx), w(ty));
            let m = self.n as i64 - 1;
            let mut s = 0.0;
            for (b, wyb) in wy.iter().enumerate() { for (a, wxa) in wx.iter().enumerate() {
                s += wxa * wyb * g((i - 1 + a as i64).clamp(0, m), (j - 1 + b as i64).clamp(0, m));
            } }
            return s;
        }
        (g(i, j) * (1.0 - tx) + g(i + 1, j) * tx) * (1.0 - ty) + (g(i, j + 1) * (1.0 - tx) + g(i + 1, j + 1) * tx) * ty
    }

    /// Ring inflow of the fine domain n x n posts from (px0, py0): coarse MFD discharge crossing into the inner box
    /// (ring-1 posts), at each crossing point split linearly between the two nearest ring-1 posts.
    pub fn inflow(&self, px0: i64, py0: i64, n: usize, mfd_exponent: f64) -> Vec<f64> {
        let mut extra = vec![0.0f64; n * n];
        // inner box in fine post units
        let (bx0, by0, bx1, by1) = ((px0 + 1) as f64, (py0 + 1) as f64, (px0 + n as i64 - 2) as f64, (py0 + n as i64 - 2) as f64);
        let st = self.st as f64;
        let inside = |ci: i64, cj: i64| { let (x, y) = ((self.fx0 + ci * self.st) as f64, (self.fy0 + cj * self.st) as f64); x >= bx0 && x <= bx1 && y >= by0 && y <= by1 };
        let ci0 = ((bx0 - self.fx0 as f64) / st).floor() as i64 - 2;
        let cj0 = ((by0 - self.fy0 as f64) / st).floor() as i64 - 2;
        let ci1 = ((bx1 - self.fx0 as f64) / st).ceil() as i64 + 2;
        let cj1 = ((by1 - self.fy0 as f64) / st).ceil() as i64 + 2;
        let nn = self.n as i64;
        let dist = [st * crate::POST_M, st * crate::POST_M * std::f64::consts::SQRT_2];
        let clen = [0.5, 0.354];
        for cj in cj0.max(0)..=cj1.min(nn - 1) { for ci in ci0.max(0)..=ci1.min(nn - 1) {
            if inside(ci, cj) { continue; }
            let c = (cj * nn + ci) as usize;
            if self.fixed[c] { continue; }
            let mut w = [0.0f64; 8];
            let mut wsum = 0.0;
            for d in 0..8 {
                let (ai, aj) = (ci + D8[d].0, cj + D8[d].1);
                if ai < 0 || aj < 0 || ai >= nn || aj >= nn { continue; }
                let s = (self.z[c] - self.z[(aj * nn + ai) as usize]) / dist[d & 1];
                if s > 0.0 { w[d] = clen[d & 1] * s.powf(mfd_exponent); wsum += w[d]; }
            }
            if wsum <= 0.0 { continue; }
            for d in 0..8 {
                if w[d] <= 0.0 { continue; }
                let (ai, aj) = (ci + D8[d].0, cj + D8[d].1);
                if !inside(ai, aj) { continue; }
                let a = self.pass[c] * w[d] / wsum;
                let (px, py) = ((self.fx0 + ci * self.st) as f64, (self.fy0 + cj * self.st) as f64);
                let (ddx, ddy) = (D8[d].0 as f64 * st, D8[d].1 as f64 * st);
                // Liang-Barsky entry parameter of the segment into the box
                let (mut t0, mut t1) = (0.0f64, 1.0f64);
                let mut ok = true;
                for (pp, qq) in [(-ddx, px - bx0), (ddx, bx1 - px), (-ddy, py - by0), (ddy, by1 - py)] {
                    if pp == 0.0 { if qq < 0.0 { ok = false; } continue; }
                    let r = qq / pp;
                    if pp < 0.0 { if r > t1 { ok = false; } else if r > t0 { t0 = r; } } else if r < t0 { ok = false; } else if r < t1 { t1 = r; }
                }
                if !ok || t0 > t1 { continue; }
                let (ex, ey) = (px + t0 * ddx - px0 as f64, py + t0 * ddy - py0 as f64);
                let (fi, fj) = (ex.clamp(1.0, (n - 2) as f64), ey.clamp(1.0, (n - 2) as f64));
                let on_x_side = (ex - 1.0).abs() < 1e-9 || (ex - (n - 2) as f64).abs() < 1e-9;
                if on_x_side {
                    let i = fi.round() as usize; let j0 = fj.floor() as usize; let t = fj - j0 as f64; let j1 = (j0 + 1).min(n - 2);
                    extra[j0 * n + i] += a * (1.0 - t); extra[j1 * n + i] += a * t;
                } else {
                    let j = fj.round() as usize; let i0 = fi.floor() as usize; let t = fi - i0 as f64; let i1 = (i0 + 1).min(n - 2);
                    extra[j * n + i0] += a * (1.0 - t); extra[j * n + i1] += a * t;
                }
            }
        } }
        extra
    }
}

/// Modes 2/3: move each ring-1 inflow to the lowest ring-1 post of the same side within `reach` posts on the conditioned
/// surface `z` (the discharge enters along the fine valley floor, not on a wall where the coarse segment happened to cross).
/// Corner posts belong to the side that comes first (x sides, then y sides). Mass is conserved.
pub fn snap_to_floor(n: usize, extra: &mut [f64], z: &[f64], reach: usize) {
    let mut out = vec![0.0f64; n * n];
    let sides: [Box<dyn Fn(usize) -> usize>; 4] = [Box::new(move |t| t * n + 1), Box::new(move |t| t * n + n - 2), Box::new(move |t| n + t), Box::new(move |t| (n - 2) * n + t)];
    let mut done = vec![false; n * n];
    for side in sides.iter() {
        for t in 1..n - 1 {
            let k = side(t);
            if done[k] { continue; }
            done[k] = true;
            if extra[k] <= 0.0 { continue; }
            let (lo, hi) = (t.saturating_sub(reach).max(1), (t + reach).min(n - 2));
            let mut best = k;
            for u in lo..=hi { let q = side(u); if z[q] < z[best] || (z[q] == z[best] && q < best) { best = q; } }
            out[best] += extra[k];
        }
    }
    extra.copy_from_slice(&out);
}

#[cfg(test)]
mod tests {
    #[test]
    fn snap_moves_inflow_to_the_valley_floor() {
        let n = 20usize;
        // valley along j = 9 (lowest), inflow arrives on the west ring-1 column at j = 11
        let z: Vec<f64> = (0..n * n).map(|k| ((k / n) as f64 - 9.0).abs()).collect();
        let mut e = vec![0.0f64; n * n];
        e[11 * n + 1] = 5.0;
        super::snap_to_floor(n, &mut e, &z, 4);
        assert_eq!(e[9 * n + 1], 5.0);
        assert!((e.iter().sum::<f64>() - 5.0).abs() < 1e-12);
    }
    use super::*;
    /// A tilted plane draining east (+x): every coarse post west of the box delivers its pass-on area through the west side;
    /// total inflow = sum of the shares of posts outside whose receivers lie inside.
    #[test]
    fn inflow_conserves_crossing_discharge() {
        let n = 40usize;
        let z: Vec<f64> = (0..n * n).map(|k| 100.0 - (k % n) as f64).collect();
        let c = Coarse { fx0: 0, fy0: 0, st: 4, n, z, dz: vec![0.0; n * n], fixed: vec![false; n * n], pass: vec![10.0; n * n], seconds: 0.0, interp: 0 };
        // fine domain 41 x 41 posts from (40, 40): inner box posts 41..79
        let e = c.inflow(40, 40, 41, 1.1);
        let tot: f64 = e.iter().sum();
        // coarse posts west of the box at ci = 10 (fine x 40), cj 10..=20, deliver 9 posts' worth through the west side (the
        // shares of cj 10 and 20 replace those cj 11 and 19 send outside); posts south / north of the box send their NE / SE
        // shares in; nothing enters from the east (no westward flow)
        let west: f64 = (0..41).map(|j| e[j * 41 + 1]).sum();
        assert!((west - 90.0).abs() < 1e-9, "west {west}");
        assert!(tot > west);
        for j in 0..41 { for i in 2..41 { if e[j * 41 + i] > 0.0 { assert!(j == 1 || j == 39, "post ({i}, {j})"); } } }
        assert_eq!(c.dz_at(41, 41), 0.0);
    }
    /// Lane RILL2: both cubic kernels reproduce a linear change exactly (partition of unity, first moment), Catmull-Rom passes
    /// through the coarse posts, the B-spline never leaves the range of its 4 x 4 neighbours (no overshoot), and the cubic
    /// reconstructions are smooth across a lattice line where bilinear has a slope kink.
    #[test]
    fn cubic_change_reconstruction() {
        let n = 12usize;
        let lin: Vec<f64> = (0..n * n).map(|k| 0.3 * (k % n) as f64 - 0.7 * (k / n) as f64).collect();
        let spike: Vec<f64> = (0..n * n).map(|k| if k % n == 6 && k / n == 6 { -10.0 } else { 0.0 }).collect();
        let mk = |dz: &Vec<f64>, interp: u32| Coarse { fx0: 0, fy0: 0, st: 4, n, z: vec![0.0; n * n], dz: dz.clone(), fixed: vec![false; n * n], pass: vec![0.0; n * n], seconds: 0.0, interp };
        for it in [0u32, 3, 4] {
            let c = mk(&lin, it);
            for (px, py) in [(9i64, 13i64), (17, 22), (24, 24), (30, 5)] {
                let want = 0.3 * px as f64 / 4.0 - 0.7 * py as f64 / 4.0;
                assert!((c.dz_at(px, py) - want).abs() < 1e-9, "interp {it} at ({px}, {py})");
            }
        }
        let (cr, bs) = (mk(&spike, 3), mk(&spike, 4));
        assert!((cr.dz_at(24, 24) + 10.0).abs() < 1e-12, "Catmull-Rom interpolates the post");
        for px in 12..36 { for py in 12..36 { let v = bs.dz_at(px, py); assert!(v <= 1e-12 && v >= -10.0, "B-spline overshoot {v}"); } }
        // slope change between consecutive fine intervals on the lattice line x = 28 vs mid-cell x = 30 (fine posts), along
        // y = 24, for a quadratic change: bilinear kinks only on the lattice line; both cubics reproduce the quadratic's
        // uniform curvature (no lattice imprint)
        let quad: Vec<f64> = (0..n * n).map(|k| ((k % n) as f64).powi(2)).collect();
        let kink = |c: &Coarse, x: i64| (c.dz_at(x + 1, 24) - 2.0 * c.dz_at(x, 24) + c.dz_at(x - 1, 24)).abs();
        let bl = mk(&quad, 0);
        assert!(kink(&bl, 28) > 0.4 && kink(&bl, 30) < 1e-9);
        for it in [3u32, 4] { let c = mk(&quad, it); assert!((kink(&c, 28) - kink(&c, 30)).abs() < 1e-9 && (kink(&c, 28) - 0.125).abs() < 1e-9, "interp {it}"); }
    }
}
