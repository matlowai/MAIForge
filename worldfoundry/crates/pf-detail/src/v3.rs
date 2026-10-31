//! pf-detail v3 (lane T2C, 2026-10-14; generator 0.3.0, switch `version=3`). PREREG:
//! `unpublished-design-record` (stages W1-W7). v3 = v2 plus:
//!
//! W1 bed packages (hierarchical keyed bedding; the class thickness-weighted mean of 1/Sc is conserved per layer);
//! W2 lateral pinch-and-swell of beds inside each layer (identity at class contacts, monotone);
//! W3 soil mantle: the remap fades out on ground well below the critical slope (no ledge rings on gentle tops);
//! W4 the relaxation reads the remap-window Sc and K, not single beds (no bed-scale notches);
//! W5 correlated keyed roughness seeds (rills nucleate irregularly);
//! W6 local mass balance of the relaxation and the remap (no per-tile raised squares from the F6 shift);
//! W7 the residual per-tile mean shift is tapered off near channel corridors.
//!
//! Everything is a pure function of the upstream fields, the key and the chart position (F9).

use crate::strata::Column;
use crate::v2::V2Params;
use crate::{k_of, sc_of, smoothstep, R};
use serde::Serialize;
use std::collections::HashMap;
use vz_keyed_rng::Scope;

pub const GENERATOR_VERSION_V3: &str = "0.3.0";
/// Tile-header rule text (design 2026-10-15 physics-first rule): what the v3 bed bands are and are not.
pub const BEDDING_NOTE: &str = "v3 bed packages, bed thicknesses and lateral pinch-and-swell are a statistical sub-resolution presentation of the modeled T0 stratigraphy (keyed draws; the class thickness-weighted mean of 1/Sc is conserved per layer; class contacts, dip and column order come from the model). They are not simulated deposition and are not newly simulated depositional events.";
pub const PURPOSE_BEDDING_V3: &str = "t2.closure.bedding.v3";
pub const PURPOSE_PINCH_V3: &str = "t2.closure.pinch.v3";
pub const PURPOSE_ROUGH_V3: &str = "t2.closure.roughness.v3";

/// W1 package median thickness (m) per LithologyCodeV1 class (0 = unbedded). Project assumptions (PREREG W1).
pub const PKG_MEDIAN_M: [f64; 12] = [0.0, 0.0, 0.0, 30.0, 12.0, 10.0, 10.0, 12.0, 8.0, 8.0, 0.0, 0.0];
pub const PKG_MEDIAN_EDIFICE_M: f64 = 25.0;

#[derive(Clone, Debug, Serialize)]
pub struct V3Params {
    pub generator_version: String,
    // W1
    pub pkg_on: bool,
    pub pkg_sigma_ln: f64,
    pub pkg_clamp_lo: f64,
    pub pkg_clamp_hi: f64,
    pub pkg_strength_sigma_ln: f64,
    pub pkg_strength_lo: f64,
    pub pkg_strength_hi: f64,
    pub bed_sigma_ln: f64,
    pub bed_clamp_lo: f64,
    pub bed_clamp_hi: f64,
    // W2
    pub pinch_on: bool,
    pub pinch_amp: f64,
    pub pinch_periods_m: [f64; 2],
    pub pinch_lattice_m: f64,
    // W3
    pub expose_on: bool,
    pub expose_lo: f64,
    pub expose_width: f64,
    // W4
    pub relax_win_on: bool,
    // W5
    pub rough_corr_on: bool,
    pub rough_lattice_posts: i64,
    pub rough_corr_weight: f64,
    pub rough_white_weight: f64,
    // W6
    pub massbal_on: bool,
    pub massbal_remap_on: bool,
    pub massbal_sigma_m: f64,
    // W7
    pub taper_c_on: bool,
    pub taper_c_m: f64,
    /// W7b (revision 1): the residual mean shift weight near channels, w = phi taper (base + exp(-max(0, excess - taper)/decay)).
    pub apron_on: bool,
    pub apron_base: f64,
    pub apron_decay_m: f64,
    /// W8 (revision 1): channel-head area-slope law A_c = A_1 / S^2 in the relaxation (A_1 = channel_head_area_m2).
    pub head_law_on: bool,
    /// Lane T2D (opt-in, None = off; serialized only when set, so v3 recipes are unchanged):
    /// X1 slope-clamp gradient half-step (m; v2/v3 use 1 m).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub clamp_h_m: Option<f64>,
    /// X2 bank grading after the fill: largest raise (m) of a free post below an adjacent carved channel post.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bank_grade_max_m: Option<f64>,
    /// X3 apron raise cap as a multiple of |tile deficit|.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub apron_cap: Option<f64>,
    /// X4 bed thickness spread inside packages from v2 (sigma_ln, clamps).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bed_sig_v2: Option<bool>,
    /// Thalweg: minimum carved interior half-width (m) so a sub-post channel is an 8-connected line.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub thalweg_hw_m: Option<f64>,
    /// Bank fins: a corridor bank post standing above both neighbours along an axis by more than this (m) is lowered to the
    /// higher of the two (lowering only).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bank_fin_m: Option<f64>,
    /// Lane RILL2 (opt-in b8_taper_after=1, standalone form of lane RILL's de_b5_chain_b8): the B8 tile-mean correction applies
    /// its edge taper AFTER the apron-cap water-filling (weight = min(phi without taper, cap) x taper), so a tile with a large
    /// deficit gets a raise that ramps to 0 over mean_taper_m at its edges instead of a capped plateau ending in a one-post step
    /// at the shared edge. None = the B8 rule unchanged.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub b8_taper_after: Option<bool>,
}

impl Default for V3Params {
    fn default() -> Self {
        V3Params {
            generator_version: GENERATOR_VERSION_V3.into(),
            pkg_on: true,
            pkg_sigma_ln: 0.6,
            pkg_clamp_lo: 0.3,
            pkg_clamp_hi: 3.0,
            pkg_strength_sigma_ln: 0.45,
            pkg_strength_lo: 0.5,
            pkg_strength_hi: 2.0,
            bed_sigma_ln: 0.8,
            bed_clamp_lo: 0.2,
            bed_clamp_hi: 4.0,
            pinch_on: true,
            pinch_amp: 0.35,
            pinch_periods_m: [40.0, 110.0],
            pinch_lattice_m: 192.0,
            expose_on: true,
            expose_lo: 0.35,
            expose_width: 0.40,
            relax_win_on: true,
            rough_corr_on: false,
            rough_lattice_posts: 6,
            rough_corr_weight: 0.8,
            rough_white_weight: 0.6,
            massbal_on: true,
            massbal_remap_on: true,
            massbal_sigma_m: 20.0,
            taper_c_on: true,
            taper_c_m: 24.0,
            apron_on: true,
            apron_base: 0.25,
            apron_decay_m: 48.0,
            head_law_on: true,
            clamp_h_m: None,
            bank_grade_max_m: None,
            apron_cap: None,
            bed_sig_v2: None,
            thalweg_hw_m: None,
            bank_fin_m: None,
            b8_taper_after: None,
        }
    }
}

// ------------------------------------------------------------------ W1 bed packages

/// Hierarchical keyed bedding of a column of layers (lith, base depth), top down. Returns (bed base depths, Sc, K).
/// Packages (lognormal thickness, keyed strength multiplier M_p) hold hard/soft bed pairs; each pair's mean of 1/Sc is
/// 1/(M_p Sc_class); then one rescale per layer makes the layer's thickness-weighted mean of 1/Sc exactly 1/Sc_class.
pub fn beds_v3(lith: &[u8], base: &[f64], pair_of: &dyn Fn(u8) -> (f64, f64), pkg_of: &dyn Fn(u8) -> f64, scope: &Scope, tag: &str,
               max_depth: f64, v2: &V2Params, v: &V3Params) -> R<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let (mut bb, mut bs, mut bk) = (Vec::new(), Vec::new(), Vec::new());
    let mut top = 0.0;
    for (k, (&l, &b)) in lith.iter().zip(base.iter()).enumerate() {
        let (hh, hs) = pair_of(l);
        let pm = pkg_of(l);
        if hh <= 0.0 || top >= max_depth {
            bb.push(b); bs.push(sc_of(l)); bk.push(k_of(l));
            top = b;
            continue;
        }
        let sk = scope.address(&format!("{tag}/layer/{k}")).map_err(|e| e.0)?;
        // counter = package index << 20 | slot (0 thickness, 1 strength, 2.. beds); attempt stays 0
        let dr = |c: u64, s: u64| -> R<f64> { Ok(sk.draw(PURPOSE_BEDDING_V3, (c << 20) | s, 0).map_err(|e| e.0)?.normal_f64()) };
        // provisional beds: (base depth, multiplier)
        let mut beds: Vec<(f64, f64)> = Vec::new();
        let mut d = top;
        let mut ip = 0u64;
        while d < b {
            let tp = (pm * (v.pkg_sigma_ln * dr(ip, 0)?).exp()).clamp(v.pkg_clamp_lo * pm, v.pkg_clamp_hi * pm);
            let mp = (v.pkg_strength_sigma_ln * dr(ip, 1)?).exp().clamp(v.pkg_strength_lo, v.pkg_strength_hi);
            let pend = (d + tp).min(b);
            let mut n = 0u64;
            // T2D X4: optionally the v2 bed spread inside packages
            let (bsg, blo, bhi) = if v.bed_sig_v2 == Some(true) { (v2.bed_sigma_ln, v2.bed_clamp_lo, v2.bed_clamp_hi) } else { (v.bed_sigma_ln, v.bed_clamp_lo, v.bed_clamp_hi) };
            while d < pend {
                let th = (hh * (bsg * dr(ip, 2 + 2 * n)?).exp()).clamp(blo * hh, bhi * hh);
                let ts = (hs * (bsg * dr(ip, 3 + 2 * n)?).exp()).clamp(blo * hs, bhi * hs).max(v2.bed_soft_min_frac * th);
                let rh = v2.bed_hard_mult;
                // pair mean of 1/m is 1/mp: th/(rh mp) + ts/(rs mp) = (th + ts)/mp  ->  rs = ts / (th + ts - th / rh)
                let rs = (ts / (th + ts - th / rh)).max(v2.bed_soft_mult_min);
                for (t, m) in [(th, rh * mp), (ts, rs * mp)] {
                    if d >= pend { break; }
                    d = (d + t).min(pend);
                    if d > max_depth { d = b; }
                    beds.push((d, m));
                }
                n += 1;
                if d >= b { break; }
            }
            ip += 1;
        }
        // layer rescale: sum t / m = sum t (exact class mean of 1/Sc over the layer)
        let (mut st, mut sw, mut prev) = (0.0f64, 0.0f64, top);
        for &(e, m) in &beds { st += e - prev; sw += (e - prev) / m; prev = e; }
        let f = if st > 0.0 { sw / st } else { 1.0 };
        for &(e, m) in &beds { let mm = m * f; bb.push(e); bs.push(sc_of(l) * mm); bk.push(k_of(l) / mm); }
        top = b;
    }
    Ok((bb, bs, bk))
}

// ------------------------------------------------------------------ W2 lateral pinch-and-swell

/// Keyed smooth lateral amplitude fields a_j(x, y) (one pair per column id), cached per lattice node.
pub struct Pinch {
    pub amp: f64,
    pub periods: [f64; 2],
    pub lattice_m: f64,
    cache: std::sync::Mutex<HashMap<(u32, u8, i64, i64), f64>>,
}

impl Pinch {
    pub fn new(v: &V3Params) -> Pinch { Pinch { amp: v.pinch_amp, periods: v.pinch_periods_m, lattice_m: v.pinch_lattice_m, cache: std::sync::Mutex::new(HashMap::new()) } }

    fn node(&self, scope: &Scope, col_id: u32, j: u8, nx: i64, ny: i64) -> f64 {
        if let Some(v) = self.cache.lock().unwrap().get(&(col_id, j, nx, ny)) { return *v; }
        let s = scope.address(&format!("pinch/{col_id}/{j}/{nx}/{ny}")).expect("scope");
        let a = self.amp * s.draw(PURPOSE_PINCH_V3, 0, 0).expect("draw").normal_f64().tanh();
        self.cache.lock().unwrap().insert((col_id, j, nx, ny), a);
        a
    }

    /// (a_1, a_2) at chart (x, y) for column `col_id` (smoothstep-bilinear between lattice nodes).
    pub fn amps(&self, scope: &Scope, col_id: u32, x: f64, y: f64) -> [f64; 2] {
        let (fx, fy) = (x / self.lattice_m, y / self.lattice_m);
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (smoothstep(fx - ix as f64), smoothstep(fy - iy as f64));
        let mut out = [0.0; 2];
        for (j, o) in out.iter_mut().enumerate() {
            let g = |a: i64, b: i64| self.node(scope, col_id, j as u8, ix + a, iy + b);
            *o = (g(0, 0) * (1.0 - tx) + g(1, 0) * tx) * (1.0 - ty) + (g(0, 1) * (1.0 - tx) + g(1, 1) * tx) * ty;
        }
        out
    }
}

/// Layer interval [top, base) of depth d in the column (above the datum: the up-bed pseudo layer [-up_h, 0)).
/// None below the last layer base or above the up-bed range (identity there).
#[inline]
fn layer_span(col: &Column, d: f64, up_h: f64) -> Option<(f64, f64)> {
    if d < 0.0 {
        if col.up_top.is_empty() || d < -up_h { return None; }
        return Some((-up_h, 0.0));
    }
    let last = *col.base_depth.last()?;
    if d >= last { return None; }
    let k = col.layer_at_depth(d);
    let top = if k == 0 { 0.0 } else { col.base_depth[k - 1] };
    Some((top, col.base_depth[k]))
}

#[inline]
fn warp_g(s: f64, t: f64, a: &[f64; 2], periods: &[f64; 2]) -> (f64, f64) {
    let (mut g, mut gp) = (0.0, 0.0);
    for j in 0..2 {
        let m = (t / periods[j]).round().max(1.0);
        let w = 2.0 * std::f64::consts::PI * m / t;
        g += a[j] / w * (w * s).sin();
        gp += a[j] * (w * s).cos();
    }
    (g, gp)
}

/// Depth (below the column datum) -> warped column depth (bed coordinates). Identity at layer contacts.
pub fn warp(col: &Column, a: &[f64; 2], periods: &[f64; 2], d: f64, up_h: f64) -> f64 {
    match layer_span(col, d, up_h) {
        None => d,
        Some((top, base)) => { let t = base - top; let s = d - top; top + s + warp_g(s, t, a, periods).0 }
    }
}

/// Inverse of `warp` (Newton, the warp maps each layer onto itself monotonically).
pub fn unwarp(col: &Column, a: &[f64; 2], periods: &[f64; 2], dw: f64, up_h: f64) -> f64 {
    match layer_span(col, dw, up_h) {
        None => dw,
        Some((top, base)) => {
            let t = base - top;
            let target = dw - top;
            let mut s = target;
            for _ in 0..30 {
                let (g, gp) = warp_g(s, t, a, periods);
                let f = s + g - target;
                if f.abs() < 1e-9 { break; }
                s = (s - f / (1.0 + gp).max(0.05)).clamp(0.0, t);
            }
            top + s
        }
    }
}

/// Up-bed range above the datum (StrataModel::ensure and the edifice build up-beds to +3000 m).
pub const UP_H: f64 = 3000.0;

// ------------------------------------------------------------------ W6 masked Gaussian low-pass

/// Normalized Gaussian low-pass of `v` over posts with mask true (n x n, row-major); sigma in posts, truncated at 3 sigma.
pub fn masked_lowpass(n: usize, v: &[f64], mask: &[bool], sigma_posts: f64) -> Vec<f64> {
    let r = (3.0 * sigma_posts).ceil() as i64;
    let k: Vec<f64> = (-r..=r).map(|i| (-0.5 * (i as f64 / sigma_posts).powi(2)).exp()).collect();
    let pass = |src: &[f64], horiz: bool| -> Vec<f64> {
        let mut out = vec![0.0f64; n * n];
        for j in 0..n { for i in 0..n {
            let mut s = 0.0;
            for (q, w) in k.iter().enumerate() {
                let o = q as i64 - r;
                let (a, b) = if horiz { (i as i64 + o, j as i64) } else { (i as i64, j as i64 + o) };
                if a < 0 || b < 0 || a >= n as i64 || b >= n as i64 { continue; }
                s += w * src[b as usize * n + a as usize];
            }
            out[j * n + i] = s;
        } }
        out
    };
    let num0: Vec<f64> = (0..n * n).map(|i| if mask[i] { v[i] } else { 0.0 }).collect();
    let den0: Vec<f64> = (0..n * n).map(|i| if mask[i] { 1.0 } else { 0.0 }).collect();
    let num = pass(&pass(&num0, true), false);
    let den = pass(&pass(&den0, true), false);
    (0..n * n).map(|i| if den[i] > 1e-9 { num[i] / den[i] } else { 0.0 }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn col() -> Column {
        Column { lith: vec![5, 3, 5], base_depth: vec![200.0, 660.0, 1500.0], bed_base: vec![200.0, 660.0, 1500.0], bed_sc: vec![0.4, 1.7, 0.4],
                 bed_k: vec![2.0, 0.6, 2.0], up_top: vec![], n_down: 3, bed_fac: vec![], bed_unit: vec![], unit_top: vec![] }
    }

    #[test]
    fn packages_conserve_layer_mean_inverse_sc() {
        let v2 = V2Params::default();
        let v = V3Params::default();
        let scope = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let pair = |l: u8| crate::v2::BED_PAIR_M[l as usize];
        let pk = |l: u8| PKG_MEDIAN_M[l as usize];
        let (bb, bs, _) = beds_v3(&[3, 5], &[460.0, 932.0], &pair, &pk, &scope, "cell/test", 7000.0, &v2, &v).unwrap();
        for (lo, hi, cls) in [(0.0, 460.0, 3u8), (460.0, 932.0, 5u8)] {
            let (mut top, mut acc, mut h) = (lo, 0.0, 0.0);
            for (b, s) in bb.iter().zip(bs.iter()) {
                if *b <= lo || *b > hi { continue; }
                acc += (b - top) / s; h += b - top; top = *b;
            }
            assert!(((acc / h) * sc_of(cls) - 1.0).abs() < 1e-9, "class {cls} ratio {}", (acc / h) * sc_of(cls));
        }
        // packages differ: the per-25 m window mean of 1/Sc in basalt varies by more than 10%
        let win = |a: f64, b: f64| { let (mut top, mut acc) = (0.0f64, 0.0); for (e, s) in bb.iter().zip(bs.iter()) {
            let (x0, x1) = (top.max(a), e.min(b)); if x1 > x0 { acc += (x1 - x0) / s; } top = *e; } acc / (b - a) };
        let ws: Vec<f64> = (0..16).map(|i| win(i as f64 * 25.0, i as f64 * 25.0 + 25.0)).collect();
        let (mn, mx) = ws.iter().fold((f64::MAX, f64::MIN), |a, &w| (a.0.min(w), a.1.max(w)));
        assert!(mx / mn > 1.1, "package contrast {mn} {mx}");
    }

    #[test]
    fn warp_is_identity_at_contacts_and_monotone_and_invertible() {
        let c = col();
        let per = [40.0, 110.0];
        for a in [[0.35, -0.35], [-0.35, 0.35], [0.2, 0.1]] {
            for &d in &[0.0, 200.0, 660.0, 1500.0, 2000.0] { assert!((warp(&c, &a, &per, d, UP_H) - d).abs() < 1e-9, "contact {d}"); }
            let mut prev = -1.0;
            let mut d = 0.0;
            while d < 1600.0 {
                let w = warp(&c, &a, &per, d, UP_H);
                assert!(w > prev, "monotone at {d}");
                prev = w;
                let back = unwarp(&c, &a, &per, w, UP_H);
                assert!((back - d).abs() < 1e-6, "inverse at {d}: {back}");
                d += 0.73;
            }
        }
    }

    #[test]
    fn masked_lowpass_of_constant_is_constant_on_mask() {
        let n = 40;
        let mask: Vec<bool> = (0..n * n).map(|i| (i % n) > 5).collect();
        let v = vec![2.5; n * n];
        let lp = masked_lowpass(n, &v, &mask, 3.0);
        for i in 0..n * n { if mask[i] { assert!((lp[i] - 2.5).abs() < 1e-12); } }
    }
}
