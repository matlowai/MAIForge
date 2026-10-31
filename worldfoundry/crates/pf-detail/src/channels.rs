//! B6 channel centrelines and cross-sections from the T1b fixed channel posts.
//!
//! Centrelines are segments between 8-adjacent T1b fixed channel posts (code 1 =
//! T1a channel, code 5 = chart-wide tributary guide); the bed falls linearly along
//! each segment between the two T1b values. Width and depth follow declared
//! hydraulic geometry from the post discharge. Everything here is pointwise.

use crate::input::{T1a, T1bWindow};
use crate::{DetailParamsV1, SECONDS_PER_YEAR, T1B_M};
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Segment {
    pub a: (f64, f64),
    pub b: (f64, f64),
    pub za: f64,
    pub zb: f64,
    pub half_w: f64,
    pub depth: f64,
    pub q: f64,
    /// 1 = T1a channel, 5 = guide.
    pub code: u8,
}

pub struct Channels {
    pub segs: Vec<Segment>,
    grid: HashMap<(i64, i64), Vec<usize>>,
    cell: f64,
}

/// What the channel network imposes at one post.
#[derive(Clone, Copy, Debug, Default)]
pub struct ChannelHit {
    /// Inside a channel (r <= w/2): the carved value (min over segments).
    pub interior: Option<f64>,
    /// In a bank zone: the bank cap (min over segments).
    pub bank: Option<f64>,
    /// Inside any corridor (r <= w/2 + bank width): fixed during relaxation.
    pub corridor: bool,
    /// min over segments of (r - corridor half-width), m (large if none near).
    pub excess_m: f64,
    /// Interior hit belongs to a T1a channel (code 1).
    pub t1a: bool,
    /// Lane INTS2n4 arm_cl: the post lies within an arm line's carved half-width or bank reach (false without the key).
    pub arm: bool,
}

impl Channels {
    pub fn build(w: &T1bWindow, t1a: &T1a, p: &DetailParamsV1, bbox: (f64, f64, f64, f64)) -> Channels {
        let (x0, y0, x1, y1) = bbox;
        let gx0 = (x0 / T1B_M).floor() as i64 - 2;
        let gy0 = (y0 / T1B_M).floor() as i64 - 2;
        let gx1 = (x1 / T1B_M).ceil() as i64 + 2;
        let gy1 = (y1 / T1B_M).ceil() as i64 + 2;
        let code = |gx: i64, gy: i64| w.idx(gx, gy).map(|i| w.fixed[i]).unwrap_or(0);
        let q_of = |gx: i64, gy: i64, c: u8| -> f64 {
            if c == 1 {
                if let Some(q) = t1a.channel_q_near(gx, gy) { return q; }
            }
            let a = w.idx(gx, gy).map(|i| w.drainage_area_m2[i] as f64).unwrap_or(0.0);
            a * p.runoff_m_per_yr / SECONDS_PER_YEAR
        };
        let mut segs = Vec::new();
        for gy in gy0..=gy1 { for gx in gx0..=gx1 {
            let ca = code(gx, gy);
            if ca != 1 && ca != 5 { continue; }
            for (dx, dy) in [(1i64, 0i64), (1, 1), (0, 1), (-1, 1)] {
                let cb = code(gx + dx, gy + dy);
                if cb != 1 && cb != 5 { continue; }
                let c = if ca == 1 && cb == 1 { 1 } else { 5 };
                let q = q_of(gx, gy, ca).max(q_of(gx + dx, gy + dy, cb));
                let wdt = (p.hg_width_coeff * q.powf(p.hg_width_exp)).max(p.min_width_m);
                let dep = (p.hg_depth_coeff * q.powf(p.hg_depth_exp)).max(p.min_depth_m);
                segs.push(Segment {
                    a: (gx as f64 * T1B_M, gy as f64 * T1B_M), b: ((gx + dx) as f64 * T1B_M, (gy + dy) as f64 * T1B_M),
                    za: w.h(gx, gy), zb: w.h(gx + dx, gy + dy), half_w: 0.5 * wdt, depth: dep, q, code: c,
                });
            }
        } }
        let cell = 64.0;
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        let reach = 0.5 * 200.0 + p.bank_width_m + p.remap_channel_taper_m; // generous: max half width 100 m
        for (i, s) in segs.iter().enumerate() {
            let (mnx, mxx) = (s.a.0.min(s.b.0) - reach, s.a.0.max(s.b.0) + reach);
            let (mny, mxy) = (s.a.1.min(s.b.1) - reach, s.a.1.max(s.b.1) + reach);
            for cy in (mny / cell).floor() as i64..=(mxy / cell).floor() as i64 {
                for cx in (mnx / cell).floor() as i64..=(mxx / cell).floor() as i64 {
                    grid.entry((cx, cy)).or_default().push(i);
                }
            }
        }
        Channels { segs, grid, cell }
    }

    pub fn hit(&self, x: f64, y: f64, p: &DetailParamsV1) -> ChannelHit {
        let mut h = ChannelHit { excess_m: 1e9, ..Default::default() };
        let Some(list) = self.grid.get(&((x / self.cell).floor() as i64, (y / self.cell).floor() as i64)) else { return h };
        for &i in list {
            let s = &self.segs[i];
            let (vx, vy) = (s.b.0 - s.a.0, s.b.1 - s.a.1);
            let l2 = vx * vx + vy * vy;
            let t = (((x - s.a.0) * vx + (y - s.a.1) * vy) / l2).clamp(0.0, 1.0);
            let (px, py) = (s.a.0 + t * vx, s.a.1 + t * vy);
            let r = (x - px).hypot(y - py);
            let bed = s.za + t * (s.zb - s.za);
            let rc = s.half_w + p.bank_width_m;
            h.excess_m = h.excess_m.min(r - rc);
            if r <= rc { h.corridor = true; }
            if r <= s.half_w {
                let v = bed - s.depth * (1.0 - (r / s.half_w).powi(2));
                if h.interior.map(|u| v < u).unwrap_or(true) { h.interior = Some(v); h.t1a = s.code == 1; }
            } else if r <= rc {
                let v = bed + p.bank_slope * (r - s.half_w);
                h.bank = Some(h.bank.map(|u| u.min(v)).unwrap_or(v));
            }
        }
        h
    }
}
