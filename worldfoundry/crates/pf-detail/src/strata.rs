//! B2 stratigraphy cut and B3 differential-resistance remap.
//!
//! The T0 column (lane SA `parcel-stacks-v2`, full depth) of a cell is placed at the
//! declared anchor with its top at the T0 cell elevation, tilted by the clamped rock
//! uplift gradient (strata rise along +g). The exposed layer of a post is the layer
//! containing its surface. B3 maps the surface through a monotone function of
//! stratigraphic height so the slope inside each layer is proportional to its
//! critical slope while the mean slope over each `remap_window_m` of section is kept.

use crate::input::T0Strata;
use crate::{k_of, sc_of, DetailParamsV1, BED_MEDIAN_M, PURPOSE_BEDDING, R};
use vz_keyed_rng::Scope;
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct Column {
    /// Lithology code per kept layer, top down.
    pub lith: Vec<u8>,
    /// Depth of each layer's base below the column top (m); the last layer extends downward.
    pub base_depth: Vec<f64>,
    /// Beds (revision 1 R1): base depth, Sc, K multiplier; the last bed extends downward.
    pub bed_base: Vec<f64>,
    pub bed_sc: Vec<f64>,
    pub bed_k: Vec<f64>,
    /// v2: tops (height above the datum) of the beds continued above the column top; their Sc/K follow the down beds
    /// in `bed_sc`/`bed_k` from index `n_down`. Empty in v1 (above the top = bed 0, as before).
    pub up_top: Vec<f64>,
    pub n_down: usize,
    /// Lane BS2 (beds=2 only; empty otherwise): facies and unit of each bed (same order as bed_sc), and the fine unit tops
    /// (signed depth below the datum, ascending; negative = the sequence continued above the datum) with their unit index.
    pub bed_fac: Vec<u8>,
    pub bed_unit: Vec<u16>,
    pub unit_top: Vec<(f64, u16)>,
}

impl Column {
    /// Index of the layer at depth `d` below the top (negative = above the top: layer 0).
    #[inline]
    pub fn bed_at_depth(&self, d: f64) -> usize {
        if d <= 0.0 {
            if self.up_top.is_empty() { return 0; }
            let h = -d;
            let i = match self.up_top.binary_search_by(|b| b.partial_cmp(&h).unwrap()) { Ok(i) => i + 1, Err(i) => i };
            return self.n_down + i.min(self.up_top.len() - 1);
        }
        match self.bed_base.binary_search_by(|b| b.partial_cmp(&d).unwrap()) {
            Ok(i) => (i + 1).min(self.bed_base.len() - 1),
            Err(i) => i.min(self.bed_base.len() - 1),
        }
    }
    #[inline]
    pub fn layer_at_depth(&self, d: f64) -> usize {
        if d <= 0.0 { return 0; }
        // first layer whose base is deeper than d
        match self.base_depth.binary_search_by(|b| b.partial_cmp(&d).unwrap()) {
            Ok(i) => (i + 1).min(self.lith.len() - 1),
            Err(i) => i.min(self.lith.len() - 1),
        }
    }
}

/// One cell's placed column and its remap table.
pub struct CellStrata {
    pub cell: u32,
    pub top_m: f64,
    pub g: (f64, f64),
    pub g_raw: (f64, f64),
    pub col: Column,
    /// Remap table: u from u0 in steps du; psi[j] = Psi~(u0 + j du) (monotone increasing).
    pub u0: f64,
    pub du: f64,
    pub psi: Vec<f64>,
    pub max_drift_m: f64,
    /// v2: window harmonic-mean critical slope 1 / mean(1/Sc) per table entry (smooth across contacts; slope clamp V6).
    pub sc_win: Vec<f32>,
    /// v3 (lane T2C W4): window arithmetic-mean K multiplier per table entry (relaxation properties). Empty before v3.
    pub k_win: Vec<f32>,
    /// Lane DE (de=1 only; empty otherwise): cumulative integral of 1/Sc over stratigraphic height per table entry (m),
    /// so the mean of 1/Sc over any height span is a table difference (edge conductance across several beds).
    pub cum_is: Vec<f64>,
}

pub struct StrataModel {
    pub anchor: (f64, f64),
    pub cells: BTreeMap<u32, CellStrata>,
    /// Lane BS (beds=2): the layer records the bed styles read; the style chosen per (cell, full-column layer index).
    pub bs: Option<std::sync::Arc<crate::strata_prov::StrataProv>>,
    pub bs_styles: BTreeMap<(u32, u16), &'static str>,
}

impl StrataModel {
    pub fn new(anchor: (f64, f64)) -> Self { StrataModel { anchor, cells: BTreeMap::new(), bs: None, bs_styles: BTreeMap::new() } }

    pub fn ensure(&mut self, t0: &T0Strata, cell: u32, p: &DetailParamsV1, scope: &Scope) -> R<()> {
        if self.cells.contains_key(&cell) { return Ok(()); }
        let c = cell as usize;
        let (a, b) = (t0.off[c] as usize, t0.off[c + 1] as usize);
        let mut lith = Vec::new();
        let mut base = Vec::new();
        let mut d = 0.0;
        for k in a..b {
            let t = t0.th[k] as f64;
            if t < p.min_layer_m { continue; }
            d += t;
            lith.push(t0.lith[k]);
            base.push(d);
        }
        if lith.is_empty() { return Err(format!("cell {cell}: empty column (oceanic or missing)")); }
        let (mut bb, mut bs, mut bk) = (Vec::new(), Vec::new(), Vec::new());
        let mut top = 0.0;
        let (mut up_top, mut n_down) = (Vec::new(), 0usize);
        if let (Some(v), true) = (p.v2.as_ref(), p.bedding_on) {
            // v2 (lane T2B, PREREG V5): hard/soft bed pairs conserving the class rock-mass Sc, continued above the top
            let pair = |l: u8| crate::v2::BED_PAIR_M[l as usize];
            let pkg = |l: u8| crate::v3::PKG_MEDIAN_M[l as usize];
            let v3 = v.v3.as_ref().filter(|w| w.pkg_on);
            if let Some(bp) = p.beds.as_ref() {
                // lane BS: record-driven bed styles per kept layer (full-column index for the record and the key)
                let rec = self.bs.as_ref().ok_or("beds=2 needs the layer records (T0 stratafull env/age dumps or prov_strata=)")?.clone();
                let kfull = kept_full_index(&t0.off, &t0.th, cell, p.min_layer_m);
                let ncol = b - a;
                let rec_of = |k: usize| -> crate::bedstyle::LayerRec {
                    let kf = kfull[k];
                    let mut r = crate::bedstyle::LayerRec::from_row(&rec.strata_row(cell, kf, t0.lith[a + kf as usize], t0.th[a + kf as usize], ncol));
                    r.lith = lith[k];
                    r
                };
                let mut top = 0.0f64;
                let mut o = crate::bedstyle::BedsOut::default();
                let plain = |o: &mut crate::bedstyle::BedsOut, b: f64, l: u8, t: f64| { o.bb.push(b); o.bs.push(sc_of(l)); o.bk.push(k_of(l)); o.bf.push(crate::bedstyle::FAC_NONE); o.bu.push(u16::MAX); o.ut.push((t, u16::MAX)); };
                for k in 0..lith.len() {
                    let bk_ = base[k];
                    if top >= p.bedding_max_depth_m || crate::v2::BED_PAIR_M[lith[k] as usize].0 <= 0.0 {
                        plain(&mut o, bk_, lith[k], top);
                    } else {
                        let st = crate::bedstyle::layer_beds(&rec_of(k), top, bk_.min(top.max(p.bedding_max_depth_m)).max(top + 1e-6), scope, &format!("cell/{cell}/layer/{}", kfull[k]), bp, &mut o)?;
                        if bk_ > p.bedding_max_depth_m { plain(&mut o, bk_, lith[k], p.bedding_max_depth_m); }
                        self.bs_styles.insert((cell, kfull[k]), st);
                    }
                    top = bk_;
                }
                let n_down = o.bb.len();
                let (mut bf, mut bu, mut ut_all) = (o.bf, o.bu, o.ut);
                let (bb, mut bs, mut bk) = (o.bb, o.bs, o.bk);
                if crate::v2::BED_PAIR_M[lith[0] as usize].0 > 0.0 {
                    let u = crate::bedstyle::up_beds(&rec_of(0), 3000.0, scope, &format!("cell/{cell}/up"), bp)?;
                    up_top = u.bb; bs.extend(u.bs); bk.extend(u.bk); bf.extend(u.bf); bu.extend(u.bu);
                    let mut uu = u.ut; uu.extend(ut_all); ut_all = uu;
                }
                ut_all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                let col = Column { lith, base_depth: base, bed_base: bb, bed_sc: bs, bed_k: bk, up_top, n_down, bed_fac: bf, bed_unit: bu, unit_top: ut_all };
                return self.finish_cell(t0, cell, col, p);
            }
            let (a, b, c) = match v3 {
                Some(w) => crate::v3::beds_v3(&lith, &base, &pair, &pkg, scope, &format!("cell/{cell}"), p.bedding_max_depth_m, v, w)?,
                None => crate::v2::beds_v2(&lith, &base, &pair, scope, &format!("cell/{cell}"), p.bedding_max_depth_m, v)?,
            };
            bb = a; bs = b; bk = c;
            n_down = bb.len();
            if crate::v2::BED_PAIR_M[lith[0] as usize].0 > 0.0 {
                let (ut, us, uk) = match v3 {
                    Some(w) => crate::v3::beds_v3(&[lith[0]], &[3000.0], &pair, &pkg, scope, &format!("cell/{cell}/up"), 3000.0, v, w)?,
                    None => crate::v2::up_beds(lith[0], crate::v2::BED_PAIR_M[lith[0] as usize], scope, &format!("cell/{cell}"), 3000.0, v)?,
                };
                up_top = ut; bs.extend(us); bk.extend(uk);
            }
            top = f64::INFINITY;
        }
        for (k, (&l, &b)) in lith.iter().zip(base.iter()).enumerate() {
            if top == f64::INFINITY { break; }
            let med = BED_MEDIAN_M[l as usize];
            if !p.bedding_on || med <= 0.0 || top >= p.bedding_max_depth_m {
                bb.push(b); bs.push(sc_of(l)); bk.push(k_of(l));
            } else {
                let sk = scope.address(&format!("cell/{cell}/layer/{k}")).map_err(|e| e.0)?;
                let mut d = top;
                let mut n = 0u64;
                while d < b {
                    let t = (med * (p.bedding_sigma_ln * sk.draw(PURPOSE_BEDDING, 2 * n, 0).map_err(|e| e.0)?.normal_f64()).exp()).clamp(0.3 * med, 4.0 * med);
                    let m = (p.bedding_mult_sigma_ln * sk.draw(PURPOSE_BEDDING, 2 * n + 1, 0).map_err(|e| e.0)?.normal_f64()).exp().clamp(0.5, 2.0);
                    d = (d + t).min(b);
                    if d > p.bedding_max_depth_m { d = b; }
                    bb.push(d); bs.push(sc_of(l) * m); bk.push(k_of(l) / m);
                    n += 1;
                }
            }
            top = b;
        }
        let col = Column { lith, base_depth: base, bed_base: bb, bed_sc: bs, bed_k: bk, up_top, n_down, bed_fac: vec![], bed_unit: vec![], unit_top: vec![] };
        self.finish_cell(t0, cell, col, p)
    }

    fn finish_cell(&mut self, t0: &T0Strata, cell: u32, col: Column, p: &DetailParamsV1) -> R<()> {
        let c = cell as usize;
        let g_raw = (t0.grad_e[c] as f64, t0.grad_n[c] as f64);
        let m = g_raw.0.hypot(g_raw.1);
        let s = if m > p.dip_clamp_tan { p.dip_clamp_tan / m } else { 1.0 };
        let g = (g_raw.0 * s, g_raw.1 * s);
        let mut cs = CellStrata { cell, top_m: t0.elevation_m[c] as f64, g, g_raw, col, u0: 0.0, du: p.remap_table_step_m, psi: Vec::new(), max_drift_m: 0.0, sc_win: Vec::new(), k_win: Vec::new(), cum_is: Vec::new() };
        build_remap(&mut cs, p)?;
        self.cells.insert(cell, cs);
        Ok(())
    }

    /// Elevation of the column top (interface 0 datum) at chart (x, y) m.
    #[inline]
    pub fn datum(&self, cs: &CellStrata, x: f64, y: f64) -> f64 {
        cs.top_m + cs.g.0 * (x - self.anchor.0) + cs.g.1 * (y - self.anchor.1)
    }

    /// (layer index, lithology code, Sc, K mult) exposed at surface z.
    #[inline]
    pub fn exposed(&self, cs: &CellStrata, x: f64, y: f64, z: f64) -> (usize, u8) {
        let d = self.datum(cs, x, y) - z;
        let k = cs.col.layer_at_depth(d);
        (k, cs.col.lith[k])
    }

    /// (Sc, K multiplier) of the bed exposed at surface z (revision 1: beds carry the closure).
    #[inline]
    pub fn bed_props(&self, cs: &CellStrata, x: f64, y: f64, z: f64) -> (f64, f64) {
        let b = cs.col.bed_at_depth(self.datum(cs, x, y) - z);
        (cs.col.bed_sc[b], cs.col.bed_k[b])
    }

    /// B3: remapped surface for a base surface z0 at (x, y).
    pub fn remap(&self, cs: &CellStrata, x: f64, y: f64, z0: f64) -> f64 {
        let dat = self.datum(cs, x, y);
        dat + inv_psi(cs, z0 - dat)
    }
}

/// Lane PV: for the kept layers of `cell` (the layers `StrataModel::ensure` keeps, thickness >= min_layer_m), their index in
/// the cell's full published column (top = 0).
pub fn kept_full_index(off: &[u32], th: &[f32], cell: u32, min_layer_m: f64) -> Vec<u16> {
    let (a, b) = (off[cell as usize] as usize, off[cell as usize + 1] as usize);
    (a..b).filter(|&j| th[j] as f64 >= min_layer_m).map(|j| (j - a).min(u16::MAX as usize) as u16).collect()
}

/// v2: window harmonic-mean Sc at stratigraphic height u (nearest table entry; outside the table: the end values).
pub fn sc_win_at(cs: &CellStrata, u: f64) -> f64 {
    if cs.sc_win.is_empty() { return 1.0; }
    let j = ((u - cs.u0) / cs.du).round().clamp(0.0, (cs.sc_win.len() - 1) as f64) as usize;
    cs.sc_win[j] as f64
}

/// Lane DE: mean of 1/Sc over the stratigraphic height span [u1, u2] (cumulative table; beyond the table the end values).
pub fn inv_sc_mean(cs: &CellStrata, u1: f64, u2: f64) -> f64 {
    let c = &cs.cum_is;
    let n = c.len();
    let at = |u: f64| -> f64 {
        let f = (u - cs.u0) / cs.du;
        if f <= 0.0 { let s0 = (c[1] - c[0]) / cs.du; return c[0] + f * cs.du * s0; }
        if f >= (n - 1) as f64 { let s1 = (c[n - 1] - c[n - 2]) / cs.du; return c[n - 1] + (f - (n - 1) as f64) * cs.du * s1; }
        let j = f.floor() as usize; let t = f - j as f64;
        c[j] + t * (c[j + 1] - c[j])
    };
    let (a, b) = if u1 <= u2 { (u1, u2) } else { (u2, u1) };
    if b - a < 1e-6 { let j = (((a - cs.u0) / cs.du).round().clamp(0.0, (n - 2) as f64)) as usize; return (c[j + 1] - c[j]) / cs.du; }
    (at(b) - at(a)) / (b - a)
}

/// v3: window mean K multiplier at stratigraphic height u (nearest table entry); 1 if no table.
pub fn k_win_at(cs: &CellStrata, u: f64) -> f64 {
    if cs.k_win.is_empty() { return 1.0; }
    let j = ((u - cs.u0) / cs.du).round().clamp(0.0, (cs.k_win.len() - 1) as f64) as usize;
    cs.k_win[j] as f64
}

/// v3 W2: remap with the lateral warp: column (bed) coordinates are the warped depths; the surface maps through them.
pub fn remap_at_w(cs: &CellStrata, dat: f64, z0: f64, a: &[f64; 2], per: &[f64; 2]) -> f64 {
    let d0w = crate::v3::warp(&cs.col, a, per, dat - z0, crate::v3::UP_H);
    let u1w = inv_psi(cs, -d0w);
    dat - crate::v3::unwarp(&cs.col, a, per, -u1w, crate::v3::UP_H)
}

/// v2: remap with an explicit datum (the edifice uses the smoothed T1b surface as its datum).
pub fn remap_at(cs: &CellStrata, dat: f64, z0: f64) -> f64 { dat + inv_psi(cs, z0 - dat) }

/// v2: the remap table builder for a column built outside `StrataModel::ensure` (the edifice).
pub fn build_remap_pub(cs: &mut CellStrata, p: &DetailParamsV1) -> R<()> { build_remap(cs, p) }

/// Stratigraphic height u = z - datum; Sc at u.
#[inline]
fn sc_at_u(col: &Column, u: f64) -> f64 { col.bed_sc[col.bed_at_depth(-u)] }

fn build_remap(cs: &mut CellStrata, p: &DetailParamsV1) -> R<()> {
    // Table range: 3 km above the top to 4 km below the column's last interface (enough for any land post).
    let depth_max = cs.col.base_depth.last().copied().unwrap_or(0.0).min(20_000.0);
    let (u_lo, u_hi) = (-(depth_max.min(6000.0)) - 500.0, 3000.0);
    let du = p.remap_table_step_m;
    let n = ((u_hi - u_lo) / du).ceil() as usize + 1;
    let s: Vec<f64> = (0..n).map(|j| 1.0 / sc_at_u(&cs.col, u_lo + j as f64 * du)).collect();
    let box_mean = |v: &[f64], w_m: f64| -> Vec<f64> {
        let h = ((w_m / du) / 2.0).round() as i64;
        let mut pre = vec![0.0f64; v.len() + 1];
        for i in 0..v.len() { pre[i + 1] = pre[i] + v[i]; }
        (0..v.len() as i64).map(|i| {
            let a = (i - h).max(0) as usize;
            let b = ((i + h + 1) as usize).min(v.len());
            (pre[b] - pre[a]) / (b - a) as f64
        }).collect()
    };
    let m = box_mean(&s, p.remap_window_m);
    cs.sc_win = if p.v2.is_some() { m.iter().map(|&v| (1.0 / v) as f32).collect() } else { Vec::new() };
    if p.v2.as_ref().map(|v| v.v3.is_some()).unwrap_or(false) {
        let kk: Vec<f64> = (0..n).map(|j| cs.col.bed_k[cs.col.bed_at_depth(-(u_lo + j as f64 * du))]).collect();
        cs.k_win = box_mean(&kk, p.remap_window_m).iter().map(|&v| v as f32).collect();
    }
    if p.de.is_some() {
        let mut c = vec![0.0f64; n];
        for j in 1..n { c[j] = c[j - 1] + 0.5 * (s[j - 1] + s[j]) * du; }
        cs.cum_is = c;
    }
    let dpsi: Vec<f64> = (0..n).map(|j| s[j] / m[j]).collect();
    let mut psi = vec![0.0f64; n];
    psi[0] = u_lo;
    for j in 1..n { psi[j] = psi[j - 1] + 0.5 * (dpsi[j - 1] + dpsi[j]) * du; }
    let drift: Vec<f64> = (0..n).map(|j| psi[j] - (u_lo + j as f64 * du)).collect();
    let low = box_mean(&drift, p.remap_drift_window_m);
    let mut max_drift: f64 = 0.0;
    for j in 0..n {
        psi[j] -= low[j];
        max_drift = max_drift.max((psi[j] - (u_lo + j as f64 * du)).abs());
    }
    for j in 1..n {
        if psi[j] <= psi[j - 1] { return Err(format!("cell {}: remap not monotone at u = {:.2}", cs.cell, u_lo + j as f64 * du)); }
    }
    cs.u0 = u_lo;
    cs.du = du;
    cs.psi = psi;
    cs.max_drift_m = max_drift;
    Ok(())
}

/// Psi~^-1(u0): u1 with Psi~(u1) = u0 (linear between table entries; identity shift outside).
fn inv_psi(cs: &CellStrata, u: f64) -> f64 {
    let ps = &cs.psi;
    let n = ps.len();
    if u <= ps[0] { return cs.u0 + (u - ps[0]); }
    if u >= ps[n - 1] { return cs.u0 + (n - 1) as f64 * cs.du + (u - ps[n - 1]); }
    let j = match ps.binary_search_by(|v| v.partial_cmp(&u).unwrap()) { Ok(j) => j, Err(j) => j - 1 };
    let t = (u - ps[j]) / (ps[j + 1] - ps[j]);
    cs.u0 + (j as f64 + t) * cs.du
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn remap_steepens_hard_and_flattens_soft() {
        let p = DetailParamsV1::v1();
        // 10 m shale over 10 m carbonate, repeated (inside the 30 m window).
        let mut lith = Vec::new(); let mut base = Vec::new(); let mut d = 0.0;
        for i in 0..200 { lith.push(if i % 2 == 0 { 5 } else { 4 }); d += 10.0; base.push(d); }
        let mut cs = CellStrata { cell: 0, top_m: 1000.0, g: (0.0, 0.0), g_raw: (0.0, 0.0), col: Column { bed_base: base.clone(), bed_sc: lith.iter().map(|&l| sc_of(l)).collect(), bed_k: lith.iter().map(|&l| k_of(l)).collect(), lith, base_depth: base, up_top: vec![], n_down: 0, bed_fac: vec![], bed_unit: vec![], unit_top: vec![] }, u0: 0.0, du: p.remap_table_step_m, psi: vec![], max_drift_m: 0.0, sc_win: vec![], k_win: vec![], cum_is: vec![] };
        build_remap(&mut cs, &p).unwrap();
        let m = StrataModel::new((0.0, 0.0));
        // a planar slope of gradient 0.5 along x crossing the section
        let z = |x: f64| 900.0 - 0.5 * x;
        let mut slopes_shale = vec![]; let mut slopes_carb = vec![];
        for i in 0..400 {
            let x = i as f64 * 1.0;
            let (a, b) = (m.remap(&cs, x, 0.0, z(x)), m.remap(&cs, x + 1.0, 0.0, z(x + 1.0)));
            let (_, l) = m.exposed(&cs, x + 0.5, 0.0, 0.5 * (a + b));
            if l == 5 { slopes_shale.push(a - b) } else { slopes_carb.push(a - b) }
        }
        let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
        assert!(mean(&slopes_carb) > 3.0 * mean(&slopes_shale), "carb {} shale {}", mean(&slopes_carb), mean(&slopes_shale));
        assert!(cs.max_drift_m < 50.0, "drift {}", cs.max_drift_m);
    }
    #[test]
    fn kept_layers_map_back_to_the_full_column() {
        // cell 1: [10 m, 0.5 m (dropped by the cut), 20 m, 0.9 m (dropped), 300 m]
        let off = vec![0u32, 2, 7];
        let th = vec![5.0f32, 5.0, 10.0, 0.5, 20.0, 0.9, 300.0];
        assert_eq!(kept_full_index(&off, &th, 1, 1.0), vec![0, 2, 4]);
        assert_eq!(kept_full_index(&off, &th, 0, 1.0), vec![0, 1]);
    }
    #[test]
    fn layer_lookup() {
        let c = Column { lith: vec![5, 4, 3], base_depth: vec![10.0, 30.0, 100.0], bed_base: vec![10.0, 30.0, 100.0], bed_sc: vec![0.4, 2.7, 1.7], bed_k: vec![2.0, 0.5, 0.6], up_top: vec![], n_down: 3, bed_fac: vec![], bed_unit: vec![], unit_top: vec![] };
        assert_eq!(c.layer_at_depth(-5.0), 0);
        assert_eq!(c.layer_at_depth(5.0), 0);
        assert_eq!(c.layer_at_depth(10.0), 1);
        assert_eq!(c.layer_at_depth(29.0), 1);
        assert_eq!(c.layer_at_depth(500.0), 2);
    }
}
