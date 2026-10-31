//! Keyed chutes: joint-controlled gullies cut into rock faces that funnel falling debris into cones.
//!
//! Seeds sit on a jittered lattice of spacing `chute_spacing_m` in WORLD chart coordinates (one keyed draw per
//! lattice cell, never per tile). A seed that lands on a producing face post walks up the steepest ascent (D8) to the
//! face top (or `chute_len_max_m`), then the chute runs down the steepest descent from there to the face foot (or
//! `chute_len_max_m`). Depth along the chute grows from the head (`chute_depth_rate` per m, capped) and tapers to 0
//! at the mouth (`chute_mouth_rate` per m), V-shaped across with half width `chute_half_width_m`, only on face posts.
//! The field is the MAX over chutes (order independent). It is a proposal for pf-detail to cut into the face; the
//! runout routes on height minus this field.
//!
//! [A] construction, not an erosion process: the chute pattern is keyed, its spacing and depth are assumptions.

use crate::params::TalusParams;
use crate::surface::Surface;
use vz_keyed_rng::Scope;

pub const NB: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];

pub fn pack(a: i64, b: i64) -> u64 { ((a as i32 as u32 as u64) << 32) | (b as i32 as u32 as u64) }

pub fn chute_field(s: &Surface, tan: &[f64], producing: &[bool], p: &TalusParams, scope: &Scope) -> Result<(Vec<f64>, u64), String> {
    let mut f = vec![0.0f64; s.n()];
    if !p.chutes || s.dx > 8.0 { return Ok((f, 0)); }
    let lam = p.chute_spacing_m;
    let dx = s.dx;
    let tan_lo = p.src_lo_deg.to_radians().tan();
    let tan_hi = p.src_hi_deg.to_radians().tan();
    let c0 = (s.x_of(0) / lam).floor() as i64;
    let c1 = (s.x_of(s.nx - 1) / lam).floor() as i64;
    let d0 = (s.y_of(0) / lam).floor() as i64;
    let d1 = (s.y_of(s.ny - 1) / lam).floor() as i64;
    let (nx, ny) = (s.nx as i64, s.ny as i64);
    let r_posts = (p.chute_half_width_m / dx).ceil() as i64;
    let mut n_chutes = 0u64;
    let is_face = |k: usize| producing[k] && tan[k] >= tan_lo;
    for cy in d0..=d1 {
        for cx in c0..=c1 {
            let w = scope.draw("talus.chute.seed", pack(cx, cy), 0).map_err(|e| e.to_string())?.u64();
            let ux = (w & 0xffff) as f64 / 65536.0;
            let uy = ((w >> 16) & 0xffff) as f64 / 65536.0;
            let um = 0.5 + ((w >> 32) & 0xffff) as f64 / 65536.0;
            let gi = (((cx as f64 + ux) * lam) / dx).round() as i64 - s.gx0;
            let gj = (((cy as f64 + uy) * lam) / dx).round() as i64 - s.gy0;
            if gi < 0 || gj < 0 || gi >= nx || gj >= ny { continue; }
            let k0 = (gj * nx + gi) as usize;
            if !(producing[k0] && tan[k0] >= tan_hi) { continue; }
            // up to the face top
            let (mut ci, mut cj) = (gi, gj);
            let mut len = 0.0;
            loop {
                let kc = (cj * nx + ci) as usize;
                let mut best: Option<(i64, i64, f64, f64)> = None;
                for (di, dj) in NB {
                    let (ni, nj) = (ci + di, cj + dj);
                    if ni < 0 || nj < 0 || ni >= nx || nj >= ny { continue; }
                    let kn = (nj * nx + ni) as usize;
                    let d = if di != 0 && dj != 0 { dx * std::f64::consts::SQRT_2 } else { dx };
                    let sl = (s.z[kn] - s.z[kc]) / d;
                    if sl > 0.0 && best.map_or(true, |b| sl > b.2) { best = Some((ni, nj, sl, d)); }
                }
                match best {
                    Some((ni, nj, _, d)) if is_face((nj * nx + ni) as usize) && len + d <= p.chute_len_max_m => { ci = ni; cj = nj; len += d; }
                    _ => break,
                }
            }
            // down from the head
            let mut path: Vec<(i64, i64, f64)> = vec![(ci, cj, 0.0)];
            let mut sl_len = 0.0;
            loop {
                let kc = (cj * nx + ci) as usize;
                let mut best: Option<(i64, i64, f64, f64)> = None;
                for (di, dj) in NB {
                    let (ni, nj) = (ci + di, cj + dj);
                    if ni < 0 || nj < 0 || ni >= nx || nj >= ny { continue; }
                    let kn = (nj * nx + ni) as usize;
                    let d = if di != 0 && dj != 0 { dx * std::f64::consts::SQRT_2 } else { dx };
                    let sl = (s.z[kc] - s.z[kn]) / d;
                    if sl > 0.0 && best.map_or(true, |b| sl > b.2) { best = Some((ni, nj, sl, d)); }
                }
                match best {
                    Some((ni, nj, _, d)) if is_face((nj * nx + ni) as usize) && sl_len + d <= p.chute_len_max_m => {
                        ci = ni; cj = nj; sl_len += d; path.push((ci, cj, sl_len));
                    }
                    _ => break,
                }
            }
            if path.len() < 3 { continue; }
            n_chutes += 1;
            let total = sl_len;
            for &(pi, pj, sk) in &path {
                let dep = um * p.chute_depth_max_m.min(p.chute_depth_rate * sk).min(p.chute_mouth_rate * (total - sk));
                if dep <= 0.0 { continue; }
                for dj in -r_posts..=r_posts {
                    for di in -r_posts..=r_posts {
                        let (ni, nj) = (pi + di, pj + dj);
                        if ni < 0 || nj < 0 || ni >= nx || nj >= ny { continue; }
                        let dist = ((di * di + dj * dj) as f64).sqrt() * dx;
                        if dist >= p.chute_half_width_m { continue; }
                        let kn = (nj * nx + ni) as usize;
                        if !is_face(kn) { continue; }
                        let v = dep * (1.0 - dist / p.chute_half_width_m);
                        if v > f[kn] { f[kn] = v; }
                    }
                }
            }
        }
    }
    Ok((f, n_chutes))
}
