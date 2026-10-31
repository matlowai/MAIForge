//! Lane T0-IMPRINT-PROD (2026-10-21, ADR 045; model specifications of 2026-10-21): the T1b network confinement rule as one
//! params-controlled function (cell / final basin / none), and the per-T0-outlet F3 and whole-cell F6 gates.
//!
//! Every switch here is opt-in. Without `conf=` the confinement is today's cell rule through the same function (every
//! byte unchanged); without `outlet_checks=1` no check output changes.

use crate::checks::{f3_outlet, F3Outlet, Outlet};
use crate::t0::{CellWater, T0};
use crate::t1b_rivernet::{Network, RiverNetParamsV1, K_LAKE, K_LAND, K_MARINE, K_T1A};
use crate::t1b_view::SkeletonView;
use crate::*;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use vz_tilekey::chart::Chart;

pub const CONF_CELL: u32 = 0;
pub const CONF_FINAL_BASIN: u32 = 1;
pub const CONF_NONE: u32 = 2;
/// Final basin over warped confinement cells (beyond the three modes of the brief; measured as an option, ADR 045).
pub const CONF_WARP: u32 = 3;
/// Lane DRAIN-FREE (design 2026-10-21 23:00): free divides within a displacement bound X of the T0 cell edges, with the
/// T1a discharge recomputed from the actual T1b catchments (flux mode).
pub const CONF_FREE: u32 = 4;
/// CX-CELL-IMPRINT: native hydraulic routing without artificial T0 cell/core walls.
/// Experimental source-revision candidate, not a declaration that legacy F3/F6 pass.
pub const CONF_TERRAIN: u32 = 5;
pub const TERRAIN_RULE: &str = "terrain_route_v1: no categorical T0 cell or inset-core barriers to lateral flow; the native global LS1 stream-power solve selects connected receivers, with catchment flux recomputed and T1a source-profile areas retained for the level law. Requires ls=1, ls_global=1 and ls_settle=1. No height offset, smoothing or altered F3/F6 tolerance. Changed basin areas and source-mean residuals remain upstream-revision findings. Not a time-evolved geological reconstruction.";
/// conf=free default bound X (mm).
pub const FREE_BOUND_DEFAULT_MM: i64 = 6_000_000;
pub const FREE_RULE: &str = "free_bounded_v1 (lane DRAIN-FREE, design 2026-10-21 23:00): a land node more than bound_mm from the nearest edge of its own T0 cell (core node; distance on the sphere to the cell's latitude and longitude bounds at the node's chart position) may drain only to a core land node of its own T0 cell, to a T1a node whose published T1a outlet basin (the T0 final basin its channel ends in) is its cell's T0 final basin, or to water and conduits; a land node within bound_mm of its cell edge (band node) may drain anywhere. Core water therefore stays in its own core and reaches its own T0 outlet (or the water), while divides between different outlets are free inside the 2 x bound_mm band around every T0 cell edge (the network growth, the level solve and capture routing all read this one rule)";
pub const FLUX_RULE: &str = "t1_catchment_flux_v1 (lane DRAIN-FREE): each T1a seed's area (discharge x year / r_ref, the effective area of the climate key) = its T1a skeleton value + the sum over the T1a receiver tree upstream of it (the seed included) of [own areas of T1b land nodes whose T1b path ends on that seed (a lake passes its water to the first T1a seed down the T1a receivers of its node) - own areas of T1b land nodes whose T1a receiver path first reaches that seed], floored at 0; recomputed from the network's receivers whenever it is rebuilt (layer level solve and capture passes, windows from the published receivers), and carried in the seed areas, the conduit injections and the T1a curves (channel drainage area, floors, T1a levels); T1b land own area = 512 m node area x R / r_ref";
/// Lane DRAIN-FIX (key conf_levels=1, with the flux on): the T1a valley levels stay the T0-carved valley's.
pub const LEVELS_RULE: &str = "t1a_valley_levels_v1 (lane DRAIN-FIX, 2026-10-22): with the flux recomputation on and the T1a level solve on (ls_t1a 1 or 2; without it the T1a levels are the skeleton levels and conduits are capped at their terminal level, so only the recipe changes), the T1a valley axis keeps the long profile it was carved for: its solved levels and the beds of the CO conduits that continue its terminal seeds to the water take, in the bed law, the T1a skeleton discharge area of each T1a node (Q x year / r_ref, the effective area of the climate key) and, on a conduit, its own land area plus the T1a discharge area of the terminal seeds injecting into it, never the recomputed T1 catchment area; the recomputed T1 discharge is still the water the channel carries (seed and conduit areas, curves' drainage area, floor widths, the revision record and the flux gate). Why: the T1a planform and receivers are fixed (F1), so a level computed from the recomputed area has no minimum-path bound; under LS1 a seed that loses land gets a higher level, its land then drains to lower neighbours and it loses more (C7 site 5: a T1a head went 1556 m to 18.4 km over two evaluations as its area fell 1.30e9 to 6.9e5 m2, under conf=free and conf=cell with conf_flux=1 alike). A reach whose discharge leaves the hostable range is flagged in the revision record (needs_upstream_revision), as the flux rule already states, not re-carved";
/// Lane DRAIN-FIX (key conf_core=path, after conf=free): the core test on the water's path instead of the neighbour pair.
pub const CORE_PATH_RULE: &str = "free_bounded_path_v1 (lane DRAIN-FIX, 2026-10-22): as free_bounded_v1, except that in the level solve (LS1, which routes every land node of the published network) a core land node may drain to a land neighbour (core or band, any cell) whose water, as routed in the same pass, ends at an outlet the core node's T0 cell accepts: a T1a node whose published T1a outlet basin is the cell's T0 final basin, a T0 lake of that final basin, or the sea under the displacement gate's sea-lenient test (the marine cell's T0 basin is the final basin, or the water body is the one the cell's T0 chain ends in); direct water, conduits and T1a neighbours keep the free_bounded_v1 tests; band nodes stay free; the initial network growth keeps the pairwise free_bounded_v1 test. Why: the pairwise test lets core water drain only into core land of its own cell, so a core node at the bound may not drain into the band of its own cell even when that band drains to its own outlet; core-rim nodes then route inward along long paths and stand up to 3 km above their conf=cell and conf=none beds in a rim at distance X from the cell edges (C6, C7 site 1 and C7 site 5; the one-post GUIDE walls of lane GUIDE-SPIKES at the bound). The path test keeps the gate's guarantee (core water reaches its own outlet) and removes the one-way wall";
/// Keyed-draw purpose of the warp control values.
pub const PURPOSE_CONF_WARP: &str = "t1b.rivernet.conf.warp";

fn is_zero_i64(v: &i64) -> bool { *v == 0 }

/// T1b network confinement (params.confinement; absent = the cell rule, recipe unchanged).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ConfinementParamsV1 {
    /// 0 cell, 1 final basin, 2 none, 3 warped basin, 4 bounded free, 5 terrain routing.
    pub mode: u32,
    pub rule: String,
    /// Mode 3: warp amplitude per axis and control spacing (mm).
    #[serde(default, skip_serializing_if = "is_zero_i64")]
    pub warp_amplitude_mm: i64,
    #[serde(default, skip_serializing_if = "is_zero_i64")]
    pub warp_spacing_mm: i64,
    /// Lane DRAIN-FREE mode 4: the displacement bound X (mm).
    #[serde(default, skip_serializing_if = "is_zero_i64")]
    pub bound_mm: i64,
    /// Lane DRAIN-FREE: 1 = the T1a seed discharge is recomputed from the T1b catchments (on with mode 4; key conf_flux=1 for
    /// the other modes).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub flux: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub flux_rule: String,
    /// Lane DRAIN-FIX (key conf_levels=1; needs the flux on): 1 = t1a_valley_levels_v1 (LEVELS_RULE). Absent = 0 = the
    /// DRAIN-FREE flux v1 behaviour (recipe unchanged).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub levels: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub levels_rule: String,
    /// Lane DRAIN-FIX (key conf_core=path; conf=free only): 1 = free_bounded_path_v1 (CORE_PATH_RULE). Absent = 0 = the
    /// pairwise free_bounded_v1 core test (recipe unchanged).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub core: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub core_rule: String,
}

fn is_zero_u32(v: &u32) -> bool { *v == 0 }

impl ConfinementParamsV1 {
    pub fn v1(mode: u32) -> R<Self> {
        let rule = match mode {
            CONF_CELL => "t0_cell_v1: a land node may drain to a node of its own T0 cell or of its T0 receiver cell (water and conduits always)",
            CONF_FINAL_BASIN => "t0_final_basin_v1: a land node may drain to a land node of its own T0 cell, of its T0 receiver cell, or of a T0 cell that ends in the same T0 final basin as its own, or to a T1a node whose channel ends in its final basin (the end of its T1a receiver chain: the T0 basin of a marine end cell, else the end cell's final basin); water and conduits always",
            CONF_NONE => "none_v1: no T0-cell confinement between network nodes (divides where the geodesic growth meets)",
            CONF_FREE => FREE_RULE,
            CONF_TERRAIN => TERRAIN_RULE,
            CONF_WARP => "t0_final_basin_warped_v1: t0_final_basin_v1 with each node's T0 cell read at its lattice position displaced by a smooth keyed field: per axis, smoothstep interpolation of control values 2w - 1 on a square chart lattice of spacing warp_spacing_mm (w = vz.keyed-rng/v1 draw, purpose t1b.rivernet.conf.warp, address cw/<i>/<j>, counter = axis, attempt 0, under the T1b draw scope), times warp_amplitude_mm; the displaced point takes the T0 cell of the nearest chart node (clamped to the chart). Divides between different T0 outlets follow a wavy band about the cell edge whose net area exchange is zero in expectation",
            _ => return Err(format!("conf: mode {mode} is not 0 (cell), 1 (basin), 2 (none), 3 (warp) 4 (free) or 5 (terrain)")),
        };
        let (a, sp) = if mode == CONF_WARP { (8_192_000, 16_384_000) } else { (0, 0) };
        let free = mode == CONF_FREE;
        let terrain = mode == CONF_TERRAIN;
        let flux = free || terrain;
        Ok(ConfinementParamsV1 { mode, rule: rule.into(), warp_amplitude_mm: a, warp_spacing_mm: sp,
            bound_mm: if free { FREE_BOUND_DEFAULT_MM } else { 0 }, flux: flux as u32, flux_rule: if flux { FLUX_RULE.into() } else { String::new() }, levels: terrain as u32, levels_rule: if terrain { LEVELS_RULE.into() } else { String::new() }, core: 0, core_rule: String::new() })
    }
}

/// Fail closed if a terrain recipe disables the coupled routing/accounting it requires.
pub fn validate_terrain(p: &RiverNetParamsV1) -> R<()> {
    let Some(c) = p.confinement.as_ref().filter(|c| c.mode == CONF_TERRAIN) else { return Ok(()) };
    if c.flux != 1 || c.levels != 1 || c.bound_mm != 0 || c.core != 0 {
        return Err("conf=terrain requires flux=1, levels=1 and no artificial bound/core".into());
    }
    if p.ls.as_ref().map(|s| s.solve) != Some(1)
        || p.ls_global.as_ref().map(|g| (g.order, g.settle)) != Some((1, 1)) {
        return Err("conf=terrain requires ls=1 ls_global=1 ls_settle=1".into());
    }
    Ok(())
}

/// Mode 3: each chart node's confinement cell (the T0 cell at its warped position).
pub fn warp_cells(cp: &ConfinementParamsV1, t0_cell: &[u32], seed: u64, draw_digest: &str) -> R<Vec<u32>> {
    let root = vz_keyed_rng::Scope::root(seed).layer(draw_digest).and_then(|s| s.chart(0)).map_err(|e| e.to_string())?;
    let sp = cp.warp_spacing_mm as f64;
    let amp = cp.warp_amplitude_mm as f64;
    if sp <= 0.0 || amp < 0.0 { return Err("conf warp: spacing must be > 0 and amplitude >= 0".into()); }
    let mut ctrl: BTreeMap<(i64, i64), [f64; 2]> = BTreeMap::new();
    let mut u = |i: i64, j: i64| -> R<[f64; 2]> {
        if let Some(v) = ctrl.get(&(i, j)) { return Ok(*v); }
        let sc = root.address(&format!("cw/{i}/{j}")).map_err(|e| e.to_string())?;
        let mut v = [0.0; 2];
        for (k, o) in v.iter_mut().enumerate() { *o = 2.0 * sc.draw(PURPOSE_CONF_WARP, k as u64, 0).map_err(|e| e.to_string())?.uniform_f64() - 1.0; }
        ctrl.insert((i, j), v);
        Ok(v)
    };
    let s = |t: f64| t * t * (3.0 - 2.0 * t);
    let mut out = vec![0u32; NODES];
    for (n, o) in out.iter_mut().enumerate() {
        let (gx, gy) = node_g(n);
        let (x, y) = ((gx * SPACING_MM) as f64, (gy * SPACING_MM) as f64);
        let (fx, fy) = (x / sp, y / sp);
        let (i0, j0) = (fx.floor() as i64, fy.floor() as i64);
        let (sx, sy) = (s(fx - i0 as f64), s(fy - j0 as f64));
        let (a, b, c, d) = (u(i0, j0)?, u(i0 + 1, j0)?, u(i0, j0 + 1)?, u(i0 + 1, j0 + 1)?);
        let w = |k: usize| { let p = a[k] + (b[k] - a[k]) * sx; let q = c[k] + (d[k] - c[k]) * sx; amp * (p + (q - p) * sy) };
        let gxw = ((x + w(0)) / SPACING_MM as f64).round() as i64;
        let gyw = ((y + w(1)) / SPACING_MM as f64).round() as i64;
        let m = node_index(gxw.clamp(G_MIN, G_MAX), gyw.clamp(G_MIN, G_MAX)).ok_or("conf warp: node")?;
        *o = t0_cell[m];
    }
    Ok(out)
}

/// The one confinement test shared by `build_network`, `t1b_ls::level_solve_s` and `t1b_fix::capture`
/// (each keeps its own water / conduit exceptions).
pub struct Confinement { mode: u32, cell: Vec<u32>, fb: BTreeMap<u32, u32>, seed_outlet: Vec<u32>, zone: Vec<u32>,
    /// Lane DRAIN-FIX conf_core=path: per T0 cell (final basin, the water body its T0 chain ends in); empty otherwise.
    path: BTreeMap<u32, (u32, u32)> }

impl Confinement {
    /// `cells`: the node -> T0 cell field; `conf_cell`: mode 3's warped cells (empty otherwise); `t1a_recv`: the T1a
    /// receiver field (basin modes: the T0 outlet a T1a node's channel ends in is found by following it).
    pub fn new(p: &RiverNetParamsV1, t0: &T0, cells: &[u32], conf_cell: &[u32], t1a_recv: &[u32]) -> Confinement {
        let mode = p.confinement.as_ref().map(|c| c.mode).unwrap_or(CONF_CELL);
        let basin = mode == CONF_FINAL_BASIN || mode == CONF_WARP || mode == CONF_FREE;
        let cell = if mode == CONF_WARP && conf_cell.len() == cells.len() { conf_cell.to_vec() } else { cells.to_vec() };
        let mut fb = BTreeMap::new();
        let mut seed_outlet = Vec::new();
        if basin {
            for &c in cells.iter().chain(cell.iter()) { if let std::collections::btree_map::Entry::Vacant(e) = fb.entry(c) { e.insert(t0.final_basin(c).unwrap_or(NONE)); } }
            seed_outlet = vec![NONE; cells.len()];
            let mut done = vec![false; cells.len()];
            let mut path = Vec::new();
            for n in 0..cells.len() {
                if done[n] { continue; }
                path.clear();
                let mut cur = n;
                let o = loop {
                    if done[cur] { break seed_outlet[cur]; }
                    path.push(cur);
                    let r = t1a_recv[cur];
                    if (r as usize) >= cells.len() || path.len() > cells.len() {
                        let c = cells[cur];
                        break if t0.is_marine(c) { t0.basin_id[c as usize] } else { *fb.entry(c).or_insert_with(|| t0.final_basin(c).unwrap_or(NONE)) };
                    }
                    cur = r as usize;
                };
                for &k in &path { seed_outlet[k] = o; done[k] = true; }
            }
        }
        // Lane DRAIN-FREE mode 4: the zones (free_zones) arrive in the conf_cell field.
        let zone = if mode == CONF_FREE && conf_cell.len() == cells.len() { conf_cell.to_vec() } else { Vec::new() };
        // Lane DRAIN-FIX conf_core=path: each cell's accepted outlet (as the displacement gate reads it).
        let mut path = BTreeMap::new();
        if mode == CONF_FREE && !zone.is_empty() && p.confinement.as_ref().map(|c| c.core == 1).unwrap_or(false) {
            let mut memo = BTreeMap::new();
            for &c in cells { if !path.contains_key(&c) { path.insert(c, cell_outlet(t0, c, &mut memo)); } }
        }
        Confinement { mode, cell, fb, seed_outlet, zone, path }
    }
    /// Lane DRAIN-FIX conf_core=path (CORE_PATH_RULE): true when the core test reads the water's path (`path_ok`).
    pub fn path_mode(&self) -> bool { !self.path.is_empty() }
    /// CORE_PATH_RULE: may land node `from` drain to a land neighbour whose water ends at the non-land node `term` (as routed;
    /// NONE when unknown)? Band nodes always; a core node when `term`'s outlet is one its T0 cell accepts (`outlet_ok`, the
    /// displacement gate's test: a T1a node's published outlet basin, a lake's T0 final basin, the sea of its basin or body).
    pub fn path_ok(&self, t0: &T0, kind: &[u8], from: usize, term: u32) -> bool {
        if self.zone.get(from).copied().unwrap_or(NONE) == NONE { return true; }
        let Some(&(fb, body)) = self.path.get(&self.cell[from]) else { return false };
        if term as usize >= kind.len() { return false; }
        let (t, c) = (term as usize, self.cell[term as usize]);
        let o = match kind[t] {
            K_T1A => if self.zone[t] != NONE { Outlet::Basin(self.zone[t]) } else { Outlet::Unknown },
            K_MARINE => Outlet::Sea { body: t0.water_body_id[c as usize], basin: (t0.basin_id[c as usize] != NONE).then_some(t0.basin_id[c as usize]) },
            K_LAKE => t0.final_basin(c).map(Outlet::Basin).unwrap_or(Outlet::Unknown),
            _ => Outlet::Unknown,
        };
        fb != NONE && outlet_ok(o, fb, body)
    }
    /// True for the modes whose T1a test applies (basin, warp, free): a T1a node receives only water of its own outlet.
    pub fn is_basin(&self) -> bool { self.mode == CONF_FINAL_BASIN || self.mode == CONF_WARP || self.mode == CONF_FREE }
    /// Mode 4 (free_bounded_v1): band nodes (zone NONE) anywhere; a core node to core land of its own cell, or to a T1a node
    /// whose channel ends in its cell's final basin.
    fn zone_allows(&self, from: usize, to: usize, to_is_t1a: bool) -> bool {
        let z = self.zone.get(from).copied().unwrap_or(NONE);
        if z == NONE { return true; }
        if to_is_t1a {
            // The T1a node's published outlet basin (in the zone field) must be the core node's own T0 final basin.
            let fa = self.fb.get(&self.cell[from]).copied().unwrap_or(NONE);
            return fa != NONE && self.zone.get(to).copied() == Some(fa);
        }
        self.zone.get(to).copied() == Some(z)
    }
    /// May the water of node `from` drain to node `to` (a land node, or a T1a node when `to_is_t1a`)?
    /// Cell mode: `to`'s T0 cell is `from`'s or its T0 receiver (today's rule). None: always. Basin modes: the cell rule or
    /// `from`'s (confinement) cell ends in the same T0 final basin as `to`'s; into a T1a node only when its channel ends in
    /// `from`'s final basin.
    #[inline]
    pub fn allows_nodes(&self, t0: &T0, from: usize, to: usize, to_is_t1a: bool) -> bool {
        let (a, b) = (self.cell[from], self.cell[to]);
        match self.mode {
            CONF_NONE | CONF_TERRAIN => true,
            CONF_FREE => self.zone_allows(from, to, to_is_t1a),
            CONF_FINAL_BASIN | CONF_WARP => {
                let fa = self.fb.get(&a).copied().unwrap_or(NONE);
                if to_is_t1a { return fa != NONE && self.seed_outlet.get(to).copied() == Some(fa); }
                a == b || t0.flow(a) == Some(b) || (fa != NONE && self.fb.get(&b).copied() == Some(fa))
            }
            _ => a == b || t0.flow(a) == Some(b),
        }
    }
}

/// Key parsing for the relief binary (layer AND windows for conf; checks only for outlet_checks). Returns true when consumed.
///   conf=cell|basin|none|warp|free|terrain (or 0..5): params.confinement (recipe; conf=cell is spelled out in the recipe, absent is the same
///   rule); conf_warp_mm=, conf_warp_spacing_mm= after conf=warp (defaults 8192000 and 16384000)
///   outlet_checks=1: per-T0-outlet F3 and whole-cell F6 gates (not part of the recipe; never changes a published byte)
///   conf_levels=1 (lane DRAIN-FIX; after conf=free or conf_flux=1): LEVELS_RULE, the T1a levels keep the T1a discharge
///   conf_core=path (lane DRAIN-FIX; after conf=free): CORE_PATH_RULE, the core test reads the water's path in LS1
pub fn parse_key(s: &str, p: &mut RiverNetParamsV1) -> R<bool> {
    if let Some(d) = s.strip_prefix("conf=") {
        let m = match d { "cell" | "0" => CONF_CELL, "basin" | "1" => CONF_FINAL_BASIN, "none" | "2" => CONF_NONE, "warp" | "3" => CONF_WARP, "free" | "4" => CONF_FREE, "terrain" | "5" => CONF_TERRAIN, _ => return Err(format!("bad value {s}")) };
        p.confinement = Some(ConfinementParamsV1::v1(m)?);
        return Ok(true);
    }
    for (k, set) in [("conf_warp_mm=", 0), ("conf_warp_spacing_mm=", 1)] {
        if let Some(d) = s.strip_prefix(k) {
            let v: i64 = d.parse().map_err(|_| format!("bad value {s}"))?;
            let c = p.confinement.as_mut().filter(|c| c.mode == CONF_WARP).ok_or(format!("{k} needs conf=warp given before it"))?;
            if set == 0 { c.warp_amplitude_mm = v; } else { c.warp_spacing_mm = v; }
            return Ok(true);
        }
    }
    // Lane DRAIN-FREE: conf_bound_mm= (after conf=free) sets X; conf_flux=1 (after any conf=) recomputes the T1a seed
    // discharge from the T1b catchments (always on with conf=free).
    if let Some(d) = s.strip_prefix("conf_bound_mm=") {
        let v: i64 = d.parse().map_err(|_| format!("bad value {s}"))?;
        if v <= 0 { return Err(format!("{s}: the bound must be > 0")); }
        let c = p.confinement.as_mut().filter(|c| c.mode == CONF_FREE).ok_or("conf_bound_mm= needs conf=free given before it")?;
        c.bound_mm = v;
        return Ok(true);
    }
    if let Some(d) = s.strip_prefix("conf_flux=") {
        let c = p.confinement.as_mut().ok_or("conf_flux= needs conf= given before it")?;
        match d {
            "1" => { c.flux = 1; c.flux_rule = FLUX_RULE.into(); }
            "0" if c.mode != CONF_FREE && c.levels == 0 => { c.flux = 0; c.flux_rule = String::new(); }
            _ => return Err(format!("bad value {s} (conf=free always recomputes the flux; conf_levels=1 needs it)")),
        }
        return Ok(true);
    }
    // Lane DRAIN-FIX: conf_core=path (after conf=free): free_bounded_path_v1 (CORE_PATH_RULE).
    if let Some(d) = s.strip_prefix("conf_core=") {
        let c = p.confinement.as_mut().filter(|c| c.mode == CONF_FREE).ok_or("conf_core= needs conf=free given before it")?;
        match d {
            "path" => { c.core = 1; c.core_rule = CORE_PATH_RULE.into(); }
            "pair" => { c.core = 0; c.core_rule = String::new(); }
            _ => return Err(format!("bad value {s} (path or pair)")),
        }
        return Ok(true);
    }
    // Lane DRAIN-FIX: conf_levels=1 (after conf=free, or conf_flux=1 on another mode): t1a_valley_levels_v1.
    if let Some(d) = s.strip_prefix("conf_levels=") {
        let c = p.confinement.as_mut().filter(|c| c.flux == 1).ok_or("conf_levels= needs the flux on (conf=free, or conf_flux=1) given before it")?;
        match d {
            "1" => { c.levels = 1; c.levels_rule = LEVELS_RULE.into(); }
            "0" => { c.levels = 0; c.levels_rule = String::new(); }
            _ => return Err(format!("bad value {s}")),
        }
        return Ok(true);
    }
    // Lane DRAIN-FIX WALL-CHECK: wall_checks=1 adds t1b_checks::wall_census to each window's checks.json (checks only).
    if let Some(d) = s.strip_prefix("wall_checks=") {
        p.wall_checks = match d { "1" => 1, "0" => 0, _ => return Err(format!("bad value {s}")) };
        return Ok(true);
    }
    // Lane LS1-GLOBAL: ls_global=, ls_settle= (t1b_ls::parse_key).
    if crate::t1b_ls::parse_key(s, p)? { return Ok(true); }
    if let Some(d) = s.strip_prefix("outlet_checks=") {
        p.outlet_checks = match d { "1" => 1, "0" => 0, _ => return Err(format!("bad value {s}")) };
        return Ok(true);
    }
    Ok(false)
}

/// Lane DRAIN-FIX (CORE_PATH_RULE): every node's destination, the first non-land node down its receivers (itself when not
/// land; NONE when the path has no end).
pub fn destinations(net: &Network) -> Vec<u32> {
    let mut out = vec![UNK; NODES];
    let mut path: Vec<usize> = Vec::new();
    for n in 0..NODES {
        if out[n] != UNK { continue; }
        path.clear();
        let mut cur = n;
        let v = loop {
            if out[cur] != UNK { break out[cur]; }
            if net.kind[cur] != K_LAND { out[cur] = cur as u32; break cur as u32; }
            path.push(cur);
            let r = net.recv[cur];
            if r as usize >= NODES || path.len() > NODES { break NONE; }
            cur = r as usize;
        };
        for &k in &path { out[k] = v; }
    }
    out
}

/// The T0 final outlet a network node's water reaches: follow receivers to the first non-land node.
pub fn network_outlets(net: &Network, skel: &SkeletonView, t0: &T0) -> Vec<Outlet> {
    let mut out = vec![Outlet::Unknown; NODES];
    let water = |c: u32| Outlet::Sea { body: t0.water_body_id[c as usize], basin: (t0.basin_id[c as usize] != NONE).then_some(t0.basin_id[c as usize]) };
    let seed_outlet = |t: usize| -> Outlet {
        match net.kind[t] {
            K_T1A => if skel.outlet_basin[t] != NONE { Outlet::Basin(skel.outlet_basin[t]) } else { Outlet::Unknown },
            K_MARINE => water(skel.t0_cell[t]),
            K_LAKE => t0.final_basin(skel.t0_cell[t]).map(Outlet::Basin).unwrap_or(Outlet::Unknown),
            _ => Outlet::Unknown,
        }
    };
    let mut done = vec![false; NODES];
    let mut path: Vec<usize> = Vec::new();
    for n in 0..NODES {
        if done[n] { continue; }
        if net.kind[n] != K_LAND { out[n] = seed_outlet(n); done[n] = true; continue; }
        path.clear();
        let mut cur = n;
        let o = loop {
            if done[cur] { break out[cur]; }
            if net.kind[cur] != K_LAND { let o = seed_outlet(cur); out[cur] = o; done[cur] = true; break o; }
            path.push(cur);
            let r = net.recv[cur];
            if r as usize >= NODES || path.len() > NODES { break Outlet::Unknown; }
            cur = r as usize;
        };
        for &k in &path { out[k] = o; done[k] = true; }
    }
    out
}

/// Per-T0-outlet F3 of the chart-wide 512 m network (ADR 045): network nodes (land and T1a seeds) by their T0 cell.
pub fn layer_f3_outlet(net: &Network, skel: &SkeletonView, t0: &T0, chart: &Chart) -> R<F3Outlet> {
    let outlets = network_outlets(net, skel, t0);
    let include: Vec<bool> = (0..NODES).map(|n| net.kind[n] == K_LAND || net.kind[n] == K_T1A).collect();
    let area: Vec<f64> = (0..NODES).map(|n| { let (x, y) = node_mm(n); let k = chart.scale_factor(x as f64, y as f64); (SPACING_MM as f64 / 1000.0).powi(2) / (k * k) }).collect();
    let mut f = f3_outlet(&skel.t0_cell, &area, &outlets, &include, t0, "t1b.rivernet 512 m network (published receivers)")?;
    // Lane DRAIN-FREE: the flux-balance gate of the flux variant.
    if net.flux.is_some() { f.flux = Some(flux_gate(net, skel, t0, chart)?); }
    Ok(f)
}

/// Whole-cell F6 gate of a window (ADR 045, FAR-F6 recommendation 3): every T0 Land cell with land posts in the window core
/// must be interior in the layer's `lane_es.f6_cell_literal` (coverage >= 0.5) and pass there. Window F6 stays a diagnostic.
pub fn window_f6_cell_gate(layer_dir: &std::path::Path, touched: &[(u32, u64)], t0: &T0) -> R<Value> {
    let b = std::fs::read(layer_dir.join("rivernet.json")).map_err(|e| format!("rivernet.json: {e}"))?;
    let v: Value = serde_json::from_slice(&b).map_err(|e| e.to_string())?;
    let lit = v.pointer("/header/lane_es/f6_cell_literal").ok_or("outlet_checks=1: the layer has no lane_es.f6_cell_literal (build the layer with es=1)")?;
    let rows: BTreeMap<u64, &Value> = lit.get("cells").and_then(Value::as_array).ok_or("f6_cell_literal.cells")?.iter()
        .filter_map(|r| r.get("cell").and_then(Value::as_u64).map(|c| (c, r))).collect();
    let mut cells = Vec::new();
    let (mut pass, mut rim) = (true, false);
    for &(c, land_posts) in touched {
        if t0.water(c) != CellWater::Land || land_posts == 0 { continue; }
        let (status, ok) = match rows.get(&(c as u64)) {
            None => ("not_in_layer_literal", false),
            Some(r) if !r.get("centre_in_chart").and_then(Value::as_bool).unwrap_or(false) => { rim = true; ("rim", false) }
            Some(r) if r.get("pass").and_then(Value::as_bool).unwrap_or(false) => ("pass", true),
            Some(_) => ("fail", false),
        };
        if !ok { pass = false; }
        let r = rows.get(&(c as u64));
        cells.push(serde_json::json!({"cell": c, "window_core_land_posts": land_posts, "status": status,
            "layer_diff_m": r.and_then(|r| r.get("diff_m")).cloned(), "layer_tolerance_m": r.and_then(|r| r.get("tolerance_m")).cloned(),
            "layer_coverage": r.and_then(|r| r.get("coverage")).cloned()}));
    }
    Ok(serde_json::json!({"rule": "whole_cell_v1 (ADR 045): every T0 Land cell with land posts in the window core is interior in the layer (coverage >= 0.5) and passes the layer's whole-cell F6 literal (land posts only, tolerance max(50 m, 10 % of T0 subgrid relief)); the window's own per-cell F6 (f6_pass, cells[].f6_*) is a partial-sample diagnostic",
        "layer_digest": v.get("digest").cloned(), "pass": pass, "rim": rim, "cells": cells}))
}

/// The opt-in `outlet_checks` block of a window's checks.json: the whole-cell F6 gate, the per-T0-outlet F3 of the layer
/// network restricted to the touched cells and their outlets, and the window-local outlet share (diagnostic).
pub fn window_block(layer: &F3Outlet, layer_dir: &std::path::Path, raw: &crate::t1b_window::RawWindow, ch: &crate::t1b_checks::WindowChecks,
    win_f3: &F3Outlet, t0: &T0) -> R<Value> {
    let mut touched: BTreeMap<u32, u64> = BTreeMap::new();
    let lat = raw.lat;
    for i in 0..lat.nodes() {
        if !lat.in_core(i) { continue; }
        let land = raw.fixed[i] != crate::t1b::MARINE && raw.fixed[i] != crate::t1b::LAKE;
        let e = touched.entry(raw.t0_cell[i]).or_default();
        if land { *e += 1; }
    }
    let tv: Vec<(u32, u64)> = touched.iter().map(|(a, b)| (*a, *b)).collect();
    let f6 = window_f6_cell_gate(layer_dir, &tv, t0)?;
    let cells: Vec<&crate::checks::F3OutletCell> = layer.cells.iter().filter(|c| touched.get(&c.cell).copied().unwrap_or(0) > 0).collect();
    let basins: std::collections::BTreeSet<u32> = cells.iter().map(|c| c.final_basin).collect();
    let brows: Vec<&crate::checks::F3OutletBasin> = layer.basins.iter().filter(|b| basins.contains(&b.basin)).collect();
    let f3_area_pass = cells.iter().all(|c| c.pass) && brows.iter().all(|b| b.pass);
    // Lane DRAIN-FREE: the flux variant gates on the layer's flux-balance gate; the area reading stays a diagnostic.
    let flux_pass = layer.flux.as_ref().and_then(|f| f.pointer("/gates/pass")).and_then(Value::as_bool);
    let f3_pass = flux_pass.unwrap_or(f3_area_pass);
    let f6_pass = f6.get("pass").and_then(Value::as_bool).unwrap_or(false);
    let mut v = serde_json::json!({
        "rule": "ADR 045 (design 2026-10-21): F3 per T0 final outlet with an area tolerance, judged on the chart-wide layer network for the T0 cells the window core touches; F6 per whole T0 cell from the layer literal; window-local F3 and F6 are diagnostics",
        "gates": {"f3_outlet_pass": f3_pass, "f6_cell_pass": f6_pass, "pass": f3_pass && f6_pass},
        "f3_outlet_layer": {"rule": layer.rule, "share_min_permille": layer.share_min_permille, "area_tolerance_permille": layer.area_tolerance_permille, "cells": cells, "basins": brows},
        "f6_cell": f6,
        "diagnostics": {"window_f3_outlet": win_f3, "window_f3_cell_pass": ch.f3_pass, "window_f6_pass": ch.f6_pass},
    });
    if let (Some(f), Value::Object(m)) = (layer.flux.as_ref(), &mut v) {
        m.insert("flux_gate_layer".into(), serde_json::json!({"gates": f.get("gates"), "bound_km": f.get("bound_km"), "f3_area_pass_diagnostic": f3_area_pass}));
    }
    Ok(v)
}

// ------------------------------------------------------------------ lane DRAIN-FREE: free divides within a bound, flux

/// Distance (m) from a point (degrees) to the nearest latitude or longitude bound of T0 cell `c` (sphere of the T0 radius).
pub fn edge_distance_m(t0: &T0, c: u32, la: f64, lo: f64) -> f64 {
    let (i, j) = t0.ij(c);
    let dlat = 180.0 / t0.nlat as f64;
    let dlon = 360.0 / t0.nlon as f64;
    let (la0, lo0) = (-90.0 + i as f64 * dlat, -180.0 + j as f64 * dlon);
    let dl = (la - la0).min(la0 + dlat - la).max(0.0);
    // Longitude offset within the cell, wrapped.
    let mut x = lo - lo0;
    x = x.rem_euclid(360.0);
    if x > dlon { x = if x - dlon < 360.0 - x { dlon } else { 0.0 }; }
    let dw = x.min(dlon - x).max(0.0) * (la * std::f64::consts::PI / 180.0).cos();
    dl.min(dw) * std::f64::consts::PI / 180.0 * t0.radius_m
}

/// Every chart node's distance (m) to the edge of its own T0 cell (`cells`), at its lattice position.
pub fn node_edge_distances(t0: &T0, chart: &Chart, cells: &[u32]) -> R<Vec<f64>> {
    let mut out = vec![0.0; NODES];
    for (n, o) in out.iter_mut().enumerate() {
        let (gx, gy) = node_g(n);
        let (la, lo) = chart.inverse_mm(gx * SPACING_MM, gy * SPACING_MM).map_err(|e| e.to_string())?;
        *o = edge_distance_m(t0, cells[n], la * vz_tilekey::chart::RAD2DEG, lo * vz_tilekey::chart::RAD2DEG);
    }
    Ok(out)
}

/// Mode 4 zones: a land node's own T0 cell when it lies more than bound_mm from that cell's edge (core), NONE in the band;
/// a T1a node (`kind` K_T1A) carries its published T1a outlet basin (the T0 final basin its channel ends in, the reading
/// the F3 and flux checks use), so the core test into a T1a node and the checks agree by construction.
pub fn free_zones(cp: &ConfinementParamsV1, t0: &T0, chart: &Chart, cells: &[u32], kind: &[u8], outlet_basin: &[u32]) -> R<Vec<u32>> {
    if cp.bound_mm <= 0 { return Err("conf=free: bound_mm must be > 0".into()); }
    let x = cp.bound_mm as f64 / 1000.0;
    let d = node_edge_distances(t0, chart, cells)?;
    Ok((0..NODES).map(|n| if kind[n] == K_T1A { outlet_basin[n] } else if d[n] > x { cells[n] } else { NONE }).collect())
}

/// Lane DRAIN-FREE: the static inputs of the flux recomputation carried by a network (flux mode only).
#[derive(Clone, Debug)]
pub struct FluxBase {
    /// Each T1a seed's area from the T1a skeleton discharge (Q x year / r_ref; 0 on other nodes).
    pub t1a_area0: Vec<f64>,
    /// The conf=free displacement bound X (mm); 0 for conf_flux=1 on another mode (no displacement gate).
    pub bound_mm: i64,
    /// The published T1a outlet basin of every node (skeleton field outlet_basin), the T0 hydrology discharge of each T1a
    /// node's T0 cell (0 elsewhere), r_ref (m/yr), and the exact land own areas of the built network (for the revision record).
    pub outlet_basin: Vec<u32>,
    pub t0_q_cell: Vec<f32>,
    pub r_ref_m_per_yr: f64,
    pub own: Vec<f64>,
    /// Lane DRAIN-FIX conf_levels=1 (LEVELS_RULE): the T1a levels and conduit beds take the T1a discharge areas.
    pub levels_t1a: bool,
}

/// needs_upstream_revision: an outlet exit whose revised discharge leaves [1/1.5, 1.5] x the discharge its T0-carved valley was
/// built for (the T1a discharge at that exit, i.e. the T0 reach discharge as T1a carries it).
pub const HOSTABLE_RATIO_MAX: f64 = 1.5;
pub const HOSTABLE_RATIO_MIN: f64 = 0.67;

/// Diagnostics of one flux recomputation (areas in m2 of effective area).
#[derive(Clone, Debug, Default, Serialize)]
pub struct FluxDiag {
    /// T1b land own area whose T1b destination seed differs from its T1a destination seed.
    pub moved_area_m2: f64,
    /// ... of which the T1a routing sent to a seed and the T1b routing to the water (sea, conduit, terminal lake).
    pub t1a_to_water_area_m2: f64,
    /// ... of which the T1a routing sent nowhere in the chart (off-chart, unrouted, T0 marine) and the T1b routing to a seed.
    pub none_to_t1a_area_m2: f64,
    pub unknown_destination_nodes: u32,
    pub seeds_floored_at_zero: u32,
}

const UNK: u32 = NONE - 2;

/// For every node, the first T1a seed (net.kind K_T1A) down its T1a receiver chain, the node itself excluded (NONE: none).
fn t1a_next_seed(net: &Network) -> Vec<u32> {
    let mut out = vec![UNK; NODES];
    let mut stack: Vec<usize> = Vec::new();
    for n0 in 0..NODES {
        if out[n0] != UNK { continue; }
        let mut cur = n0;
        let v = loop {
            if out[cur] != UNK { break out[cur]; }
            let r = net.t1a_recv[cur];
            if r as usize >= NODES || stack.len() > NODES { out[cur] = NONE; break NONE; }
            if net.kind[r as usize] == K_T1A { out[cur] = r; break r; }
            stack.push(cur);
            cur = r as usize;
        };
        while let Some(k) = stack.pop() { out[k] = v; }
    }
    out
}

/// For every T1b land node (conduits excluded), its T1b destination: the T1a seed its path ends on (a lake passes on to its
/// next T1a seed), NONE for the water (sea, conduit, terminal lake), UNK when the path has no end.
fn t1b_dest(net: &Network, tnext: &[u32]) -> Vec<u32> {
    let land = |n: usize| net.kind[n] == K_LAND && !net.conduit[n];
    let mut out = vec![UNK; NODES];
    let mut done = vec![false; NODES];
    let mut stack: Vec<usize> = Vec::new();
    for n0 in 0..NODES {
        if !land(n0) || done[n0] { continue; }
        let mut cur = n0;
        let v = loop {
            if done[cur] { break out[cur]; }
            stack.push(cur);
            let r = net.recv[cur];
            if r as usize >= NODES || stack.len() > NODES { break UNK; }
            let r = r as usize;
            if land(r) { cur = r; continue; }
            break match net.kind[r] { K_T1A => r as u32, K_LAKE => tnext[r], _ => NONE };
        };
        while let Some(k) = stack.pop() { out[k] = v; done[k] = true; }
    }
    out
}

/// Lane DRAIN-FREE (FLUX_RULE): the T1a seed areas of `net`'s routing; `own` = exact land own areas (own_area of the
/// network built by build_network). Returns per node the seed area (0 off the seeds) and the diagnostics.
pub fn flux_seed_areas(net: &Network, own: &[f64]) -> (Vec<f64>, FluxDiag) {
    let tnext = t1a_next_seed(net);
    let db = t1b_dest(net, &tnext);
    let mut delta = vec![0.0f64; NODES];
    let mut d = FluxDiag::default();
    for u in 0..NODES {
        if net.kind[u] != K_LAND || net.conduit[u] { continue; }
        let w = own[u];
        let (a, b) = (tnext[u], db[u]);
        if b == UNK { d.unknown_destination_nodes += 1; }
        let b = if b == UNK { NONE } else { b };
        if a == b { continue; }
        d.moved_area_m2 += w;
        if a != NONE { delta[a as usize] -= w; if b == NONE { d.t1a_to_water_area_m2 += w; } }
        if b != NONE { delta[b as usize] += w; if a == NONE { d.none_to_t1a_area_m2 += w; } }
    }
    // Carry the changes down the T1a receiver forest (upstream first).
    let mut indeg = vec![0u32; NODES];
    for n in 0..NODES { let r = net.t1a_recv[n]; if (r as usize) < NODES { indeg[r as usize] += 1; } }
    let mut q: std::collections::VecDeque<usize> = (0..NODES).filter(|&n| indeg[n] == 0).collect();
    while let Some(n) = q.pop_front() {
        let r = net.t1a_recv[n];
        if (r as usize) < NODES { let r = r as usize; delta[r] += delta[n]; indeg[r] -= 1; if indeg[r] == 0 { q.push_back(r); } }
    }
    let mut seeds = vec![0.0f64; NODES];
    for t in 0..NODES {
        if net.kind[t] != K_T1A { continue; }
        let v = net.flux.as_ref().map(|f| f.t1a_area0[t]).unwrap_or(0.0) + delta[t];
        if v < 0.0 { d.seeds_floored_at_zero += 1; }
        seeds[t] = v.max(0.0);
    }
    (seeds, d)
}

/// The conduits each terminal T1a seed injects into (the first hop of its conduit path, when it is construction land).
fn conduit_injections(net: &Network) -> Vec<(usize, usize)> {
    net.conduit_paths.iter().filter(|p| p.len() >= 2 && net.kind[p[0] as usize] == K_T1A && net.kind[p[1] as usize] == K_LAND)
        .map(|p| (p[0] as usize, p[1] as usize)).collect()
}

/// In `t1b_fix::rebuild`: the seed areas into `area` (own areas before accumulation) and the conduit injections changed by
/// seed - own[seed] (own[seed] = the seed area the injection of `own` carries).
pub fn flux_set_seeds(net: &Network, seeds: &[f64], own: &[f64], area: &mut [f64]) {
    for (t, c) in conduit_injections(net) { area[c] += seeds[t] - own[t]; }
    for t in 0..NODES { if net.kind[t] == K_T1A { area[t] = seeds[t]; } }
}

/// In `build_network` (flux mode): the seed areas of the built network, the conduit injections and the conduit chains below
/// them updated.
pub fn flux_apply_built(net: &mut Network) {
    let own = crate::t1b_fix::own_area(net);
    if let Some(f) = net.flux.as_mut() { f.own = own.clone(); }
    let (seeds, _) = flux_seed_areas(net, &own);
    for (t, c) in conduit_injections(net) {
        let dv = seeds[t] - net.area_m2[t];
        let mut cur = c;
        for _ in 0..NODES {
            net.area_m2[cur] += dv;
            let r = net.recv[cur];
            if r as usize >= NODES || net.kind[r as usize] != K_LAND { break; }
            cur = r as usize;
        }
    }
    for t in 0..NODES { if net.kind[t] == K_T1A { net.area_m2[t] = seeds[t]; } }
}

/// Lane DRAIN-FIX (LEVELS_RULE): the areas the bed law takes for the T1a levels and the conduit beds: `area_m2` with each
/// T1a seed at its T1a discharge area and each conduit chain below an injecting terminal seed shifted by (T1a area - seed
/// area), i.e. the conduit's own land plus the T1a discharge it was built to carry (as without the flux). None unless
/// conf_levels=1 (callers then read `area_m2`, unchanged).
pub fn level_areas(net: &Network) -> Option<Vec<f64>> {
    let f = net.flux.as_ref().filter(|f| f.levels_t1a)?;
    let mut a = net.area_m2.clone();
    for (t, c) in conduit_injections(net) {
        let dv = f.t1a_area0[t] - net.area_m2[t];
        let mut cur = c;
        for _ in 0..NODES {
            a[cur] += dv;
            let r = net.recv[cur];
            if r as usize >= NODES || net.kind[r as usize] != K_LAND { break; }
            cur = r as usize;
        }
    }
    for t in 0..NODES { if net.kind[t] == K_T1A { a[t] = f.t1a_area0[t]; } }
    Some(a)
}

/// T0 cell -> (final basin, the marine water body its T0 chain ends in or NONE).
fn cell_outlet(t0: &T0, c: u32, memo: &mut BTreeMap<u32, (u32, u32)>) -> (u32, u32) {
    if let Some(v) = memo.get(&c) { return *v; }
    let v = match t0.chain_from(c) {
        Ok(ch) => { let end = *ch.last().unwrap(); (t0.basin_id[end as usize], if t0.is_marine(end) { t0.water_body_id[end as usize] } else { NONE }) }
        Err(_) => (NONE, NONE),
    };
    memo.insert(c, v);
    v
}

/// The F3 (a) sea-lenient test: may cell water ending at `o` count as reaching the cell's own T0 outlet?
pub fn outlet_ok(o: Outlet, fb: u32, body: u32) -> bool {
    match o {
        Outlet::Basin(b) => b == fb,
        Outlet::Sea { body: bd, basin } => basin == Some(fb) || (body != NONE && (bd == body || bd == NONE)),
        Outlet::Unknown => false,
    }
}

/// Lane DRAIN-FREE per-outlet revision record (published in the layer header as `outlet_revision` in flux mode; downstream
/// consumers read it instead of the T0 / T1a value). A construction-time reconciliation between T0 / T1a and the T1b
/// catchments, not a detail-pass reroute: for each T0 outlet (the T1a outlet basin of each T1a exit seed, i.e. a seed from
/// which the T1a receivers reach no further seed), the T0 reference (T1a discharge at the exits = the T0 reach discharge as
/// T1a carries it, its effective area, the T0 hydrology cell discharge), the revised T1 discharge and effective area, the
/// delta and its causes, and `needs_upstream_revision` where an exit leaves the hostable range.
pub fn outlet_revision(net: &Network) -> Value {
    let Some(f) = net.flux.as_ref() else { return Value::Null };
    let tnext = t1a_next_seed(net);
    let db = t1b_dest(net, &tnext);
    let q = |a: f64| a * f.r_ref_m_per_yr / SECONDS_PER_YEAR;
    let exit_of = |t: usize| -> usize { let mut c = t; let mut g = 0; while tnext[c] != NONE && g <= NODES { c = tnext[c] as usize; g += 1; } c };
    // outlet -> [t0 ref area, t1 area, lost to water, lost to other outlets, gained from other outlets, gained from no T1a dest]
    let mut rows: BTreeMap<u32, [f64; 6]> = BTreeMap::new();
    let mut exits: BTreeMap<u32, Vec<Value>> = BTreeMap::new();
    let mut flagged_exits = 0u32;
    for t in 0..NODES {
        if net.kind[t] != K_T1A || tnext[t] != NONE { continue; }
        let o = f.outlet_basin[t];
        let e = rows.entry(o).or_default();
        e[0] += f.t1a_area0[t]; e[1] += net.area_m2[t];
        let ratio = if f.t1a_area0[t] > 0.0 { net.area_m2[t] / f.t1a_area0[t] } else { f64::INFINITY };
        let flag = q(f.t1a_area0[t]) >= 1.0 && !(HOSTABLE_RATIO_MIN..=HOSTABLE_RATIO_MAX).contains(&ratio);
        if flag { flagged_exits += 1; }
        let (gx, gy) = node_g(t);
        exits.entry(o).or_default().push(serde_json::json!({"node": t, "g": [gx, gy], "t0_cell_q_m3_s": f.t0_q_cell[t], "q_t0_ref_m3_s": q(f.t1a_area0[t]), "q_t1_m3_s": q(net.area_m2[t]),
            "ratio": if ratio.is_finite() { Value::from(ratio) } else { Value::Null }, "needs_upstream_revision": flag}));
    }
    for u in 0..NODES {
        if net.kind[u] != K_LAND || net.conduit[u] { continue; }
        let (a, b) = (tnext[u], if db[u] == UNK { NONE } else { db[u] });
        if a == b { continue; }
        let oa = (a != NONE).then(|| f.outlet_basin[exit_of(a as usize)]);
        let ob = (b != NONE).then(|| f.outlet_basin[exit_of(b as usize)]);
        if oa == ob { continue; }
        let w = f.own[u];
        match (oa, ob) {
            (Some(x), None) => rows.entry(x).or_default()[2] += w,
            (Some(x), Some(y)) => { rows.entry(x).or_default()[3] += w; rows.entry(y).or_default()[4] += w; }
            (None, Some(y)) => rows.entry(y).or_default()[5] += w,
            _ => {}
        }
    }
    let n_flag = exits.values().filter(|v| v.iter().any(|e| e["needs_upstream_revision"] == true)).count();
    let outlets: Vec<Value> = rows.iter().map(|(&o, e)| serde_json::json!({
        "outlet_basin": o,
        "t0_reference": {"q_m3_s": q(e[0]), "area_eff_km2": e[0] / 1e6, "source": "T1a discharge at the exits (T0 runoff and T0 chart-entry inflow, as T1a carries the T0 reach)"},
        "t1": {"q_m3_s": q(e[1]), "area_eff_km2": e[1] / 1e6},
        "delta_m3_s": q(e[1] - e[0]), "delta_permille": if e[0] > 0.0 { 1000.0 * (e[1] - e[0]) / e[0] } else { 0.0 },
        "causes_m3_s": {"divide_displacement": q(e[4] - e[3]), "divide_displacement_lost": q(e[3]), "divide_displacement_gained": q(e[4]),
            "t1b_direct_to_water": -q(e[2]), "rim_capture_from_off_chart_routing": q(e[5])},
        "exits": exits.get(&o).cloned().unwrap_or_default(),
        "needs_upstream_revision": exits.get(&o).map(|v| v.iter().any(|x| x["needs_upstream_revision"] == true)).unwrap_or(false),
    })).collect();
    serde_json::json!({
        "rule": "t1_outlet_revision_v1 (lane DRAIN-FREE): construction-time reconciliation of T0 / T1a outlet discharge with the T1b catchments (flux_rule); consumers read t1 instead of the T0 reference; causes: divide_displacement = water gained from minus lost to other T0 outlets by moved divides; t1b_direct_to_water = T1b land draining straight to the sea, an estuary, a conduit or a terminal lake that T1a routed to this outlet's trunk; rim_capture = T1b land the T1a routed off the chart; needs_upstream_revision: an exit with T0-reference discharge >= 1 m3/s whose revised discharge is outside [0.67, 1.5] x it (the T0-carved valley was not built for it; no carving follows from this record)",
        "flux_rule": FLUX_RULE, "hostable_ratio": [HOSTABLE_RATIO_MIN, HOSTABLE_RATIO_MAX],
        "outlets_needing_upstream_revision": n_flag, "exits_needing_upstream_revision": flagged_exits,
        "outlets": outlets,
    })
}

/// Lane DRAIN-FREE flux-balance gate on the published network (replaces the F3 area gate for the flux variant; the F3
/// rows stay as diagnostics): (1) every land node's water has a valid destination; (2) per T0 outlet the carried discharge
/// at its T1a exits equals the T1a discharge plus the T1b catchment change (recomputed independently by walking every
/// node's paths) to 1e-6 relative; (3) displacement: no land node farther than X from its T0 cell edge ends at another
/// outlet (X = the conf=free bound; a table over several X is reported for every mode). Topology across windows is the
/// window checks plus the seam certificate (windows rebuild this network from the published receivers).
pub fn flux_gate(net: &Network, skel: &SkeletonView, t0: &T0, chart: &Chart) -> R<Value> {
    let bound_mm = net.flux.as_ref().map(|f| f.bound_mm).unwrap_or(0);
    let flux = net.flux.is_some();
    let own = own_exact(t0, chart)?;
    let tnext = t1a_next_seed(net);
    let db = t1b_dest(net, &tnext);
    // (1) destinations.
    let (mut to_t1a, mut to_lake_t1a, mut to_water, mut unknown, mut cyc) = (0u32, 0u32, 0u32, 0u32, 0u32);
    for u in 0..NODES {
        if net.kind[u] != K_LAND { continue; }
        let mut cur = u; let mut k = 0usize;
        loop {
            let r = net.recv[cur];
            if r as usize >= NODES { unknown += 1; break; }
            let r = r as usize;
            match net.kind[r] { K_LAND => { cur = r; k += 1; if k > NODES { cyc += 1; break; } } K_T1A => { to_t1a += 1; break; } K_LAKE => { if tnext[r] != NONE { to_lake_t1a += 1; } else { to_water += 1; } break; } _ => { to_water += 1; break; } }
        }
    }
    let r_ref = t0.runoff_m_per_yr;
    let q = |a: f64| a * r_ref / SECONDS_PER_YEAR;
    // T1a exits: seeds from which the T1a receivers reach no further seed.
    let exits: Vec<usize> = (0..NODES).filter(|&t| net.kind[t] == K_T1A && tnext[t] == NONE).collect();
    let mut exit_of = vec![NONE; NODES];
    for t in 0..NODES { if net.kind[t] == K_T1A { exit_of[t] = if tnext[t] == NONE { t as u32 } else { NONE }; } }
    // Seed -> its exit (follow tnext).
    for t in 0..NODES {
        if net.kind[t] != K_T1A || exit_of[t] != NONE { continue; }
        let mut c = t; let mut g = 0;
        while tnext[c] != NONE && g <= NODES { c = tnext[c] as usize; g += 1; }
        exit_of[t] = c as u32;
    }
    // (2) per outlet: T1a discharge and carried discharge at the exits, and the independent expectation.
    let mut rows: BTreeMap<u32, [f64; 5]> = BTreeMap::new(); // outlet -> [t1a area, carried area, expected area, exits, t0 q]
    // outlet -> [lost to the water, lost to other outlets, gained from other outlets, gained from no T1a destination] (area)
    let mut parts: BTreeMap<u32, [f64; 4]> = BTreeMap::new();
    let ob = |t: usize| skel.outlet_basin[t];
    let a_t1a = |t: usize| skel.discharge_m3_s[t] as f64 * SECONDS_PER_YEAR / r_ref;
    for &t in &exits {
        let e = rows.entry(ob(t)).or_default();
        e[0] += a_t1a(t); e[1] += net.area_m2[t]; e[2] += a_t1a(t); e[3] += 1.0;
        e[4] = e[4].max(t0.discharge_m3_s[skel.t0_cell[t] as usize] as f64);
    }
    for u in 0..NODES {
        if net.kind[u] != K_LAND || net.conduit[u] { continue; }
        let (a, b) = (tnext[u], db[u]);
        let b = if b == UNK { NONE } else { b };
        if a == b { continue; }
        let oa = (a != NONE).then(|| ob(exit_of[a as usize] as usize));
        let obb = (b != NONE).then(|| ob(exit_of[b as usize] as usize));
        if let Some(o) = oa { if let Some(e) = rows.get_mut(&o) { e[2] -= own[u]; } }
        if let Some(o) = obb { if let Some(e) = rows.get_mut(&o) { e[2] += own[u]; } }
        if oa == obb { continue; }
        match (oa, obb) {
            (Some(x), None) => parts.entry(x).or_default()[0] += own[u],
            (Some(x), Some(y)) => { parts.entry(x).or_default()[1] += own[u]; parts.entry(y).or_default()[2] += own[u]; }
            (None, Some(y)) => parts.entry(y).or_default()[3] += own[u],
            _ => {}
        }
    }
    let mut worst_rel = 0.0f64;
    let outlets: Vec<Value> = rows.iter().map(|(&o, e)| {
        let rel = if e[2].abs() > 0.0 { (e[1] - e[2]) / e[2].abs().max(1.0) } else { 0.0 };
        if flux { worst_rel = worst_rel.max(rel.abs()); }
        serde_json::json!({"outlet_basin": o, "exits": e[3] as u32, "q_t0_cell_max_m3_s": e[4], "q_t1a_m3_s": q(e[0]), "q_t1_carried_m3_s": q(e[1]),
            "q_t1_expected_m3_s": q(e[2]), "change_m3_s": q(e[1] - e[0]), "change_permille": if e[0] > 0.0 { 1000.0 * (e[1] - e[0]) / e[0] } else { 0.0 },
            "closure_rel": rel, "lost_to_water_m3_s": q(parts.get(&o).map(|p| p[0]).unwrap_or(0.0)), "lost_to_other_outlets_m3_s": q(parts.get(&o).map(|p| p[1]).unwrap_or(0.0)),
            "gained_from_other_outlets_m3_s": q(parts.get(&o).map(|p| p[2]).unwrap_or(0.0)), "gained_from_no_t1a_destination_m3_s": q(parts.get(&o).map(|p| p[3]).unwrap_or(0.0))})
    }).collect();
    // (3) displacement.
    let dist = node_edge_distances(t0, chart, &skel.t0_cell)?;
    let outl = network_outlets(net, skel, t0);
    let mut memo = BTreeMap::new();
    let xs_km = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0];
    let mut viol = vec![(0u32, 0.0f64); xs_km.len()];
    let mut max_d = 0.0f64;
    let (mut land_n, mut wrong_n) = (0u32, 0u32);
    for u in 0..NODES {
        if net.kind[u] != K_LAND || net.conduit[u] { continue; }
        let c = skel.t0_cell[u];
        if t0.water(c) != CellWater::Land { continue; }
        land_n += 1;
        let (fb, body) = cell_outlet(t0, c, &mut memo);
        if outlet_ok(outl[u], fb, body) { continue; }
        wrong_n += 1;
        max_d = max_d.max(dist[u]);
        for (k, x) in xs_km.iter().enumerate() { if dist[u] > x * 1000.0 { viol[k].0 += 1; viol[k].1 += own[u]; } }
    }
    let x_km = bound_mm as f64 / 1e6;
    // Gated cells as F3: at least 500 permille of the T0 cell inside the chart (rim cells are reported, not gated).
    let mut inside: BTreeMap<u32, f64> = BTreeMap::new();
    for u in 0..NODES { let (x, y) = node_mm(u); let k = chart.scale_factor(x as f64, y as f64); *inside.entry(skel.t0_cell[u]).or_default() += (SPACING_MM as f64 / 1000.0).powi(2) / (k * k); }
    let gated = |c: u32| inside.get(&c).copied().unwrap_or(0.0) / t0.cell_area_m2(c) * 1000.0 >= crate::checks::F3_OUTLET_MIN_INSIDE_PERMILLE;
    let (mut bound_viol, mut bound_viol_rim) = (0u32, 0u32);
    let mut viol_cells: BTreeMap<u32, u32> = BTreeMap::new();
    for u in 0..NODES {
        if net.kind[u] != K_LAND || net.conduit[u] || t0.water(skel.t0_cell[u]) != CellWater::Land || dist[u] <= x_km * 1000.0 { continue; }
        let (fb, body) = cell_outlet(t0, skel.t0_cell[u], &mut memo);
        if outlet_ok(outl[u], fb, body) { continue; }
        if gated(skel.t0_cell[u]) { bound_viol += 1; *viol_cells.entry(skel.t0_cell[u]).or_default() += 1; } else { bound_viol_rim += 1; }
    }
    // T0 lakes (design review review): T1b lake stages are the T0 stages, so a moved catchment changes the lake's inflow; per T0 lake
    // body seen in the chart, the T1b land nodes whose water passes it (their path ends in it, or ends on a T1a seed whose
    // T1a receivers pass it) against the same nodes' T0 reading (their T0 cell chain passes a lake cell of that body).
    // Plain node areas; a change > 100 permille is flagged (reported, not gated).
    let lake_body = |n: usize| t0.water_body_id[skel.t0_cell[n] as usize];
    let mut seed_lake = vec![UNK; NODES];
    let mut t1_lake = |s0: usize| -> u32 {
        if seed_lake[s0] != UNK { return seed_lake[s0]; }
        let mut c = s0; let mut g = 0; let mut v = NONE;
        loop {
            if net.kind[c] == K_LAKE { v = lake_body(c); break; }
            let r = net.t1a_recv[c];
            if r as usize >= NODES || g > NODES { break; }
            c = r as usize; g += 1;
        }
        seed_lake[s0] = v; v
    };
    let mut cell_lake: BTreeMap<u32, u32> = BTreeMap::new();
    let a_node = (SPACING_MM as f64 / 1000.0).powi(2);
    let mut lakes: BTreeMap<u32, [f64; 2]> = BTreeMap::new();
    for n in 0..NODES { if net.kind[n] == K_LAKE { lakes.entry(lake_body(n)).or_default(); } }
    for u in 0..NODES {
        if net.kind[u] != K_LAND || net.conduit[u] || lakes.is_empty() { continue; }
        let c = skel.t0_cell[u];
        let l0 = *cell_lake.entry(c).or_insert_with(|| t0.chain_from(c).ok().and_then(|ch| ch.into_iter().find(|&x| t0.water(x) == CellWater::Lake)).map(|x| t0.water_body_id[x as usize]).unwrap_or(NONE));
        let mut x = u; let mut g = 0;
        let l1 = loop {
            let r = net.recv[x];
            if r as usize >= NODES || g > NODES { break NONE; }
            let r = r as usize;
            if net.kind[r] == K_LAND && !net.conduit[r] { x = r; g += 1; continue; }
            break match net.kind[r] { K_LAKE => lake_body(r), K_T1A => t1_lake(r), _ => NONE };
        };
        if let Some(e) = lakes.get_mut(&l0) { e[0] += a_node; }
        if let Some(e) = lakes.get_mut(&l1) { e[1] += a_node; }
    }
    let lake_rows: Vec<Value> = lakes.iter().map(|(&b, e)| { let rel = if e[0] > 0.0 { 1000.0 * (e[1] - e[0]) / e[0] } else { 0.0 };
        serde_json::json!({"water_body": b, "t0_catchment_km2": e[0] / 1e6, "t1_catchment_km2": e[1] / 1e6, "change_permille": rel, "flag_over_100_permille": rel.abs() > 100.0 || (e[0] == 0.0 && e[1] > 0.0)}) }).collect();
    let lakes_flagged = lake_rows.iter().filter(|r| r["flag_over_100_permille"] == true).count();
    let dest_pass = unknown == 0 && cyc == 0;
    let bal_pass = !flux || worst_rel <= 1e-6;
    let disp_pass = bound_mm > 0 && bound_viol == 0;
    Ok(serde_json::json!({
        "rule": "flux_balance_v1 (lane DRAIN-FREE, design 2026-10-21 23:00): (1) destinations: every T1b land node's path ends on a T1a seed, the sea, a lake (passing on down its T1a outflow) or a conduit; (2) per T0 outlet (T1a outlet basin of each T1a exit seed) the carried exit discharge equals the T1a discharge plus the T1b catchment change, recomputed independently from every node's T1a and T1b paths, to 1e-6 relative; (3) displacement: no T1b land node of a gated T0 Land cell (>= 500 permille inside the chart, as F3; rim cells reported) farther than X from its cell edge ends at an outlet other than its cell's (F3 (a) sea-lenient test); the old F3 per-cell and per-outlet area numbers stay as diagnostics",
        "flux_rule": if flux { FLUX_RULE } else { "off (T1a discharge carried unchanged)" },
        "bound_km": x_km,
        "gates": {"destinations_pass": dest_pass, "balance_pass": bal_pass, "displacement_pass": disp_pass, "pass": dest_pass && bal_pass && disp_pass},
        "destinations": {"land_nodes_to_t1a": to_t1a, "to_lake_then_t1a": to_lake_t1a, "to_water": to_water, "no_end": unknown, "cycles": cyc},
        "balance": {"worst_closure_rel": worst_rel, "recomputation": flux_seed_areas(net, &own).1, "outlets": outlets},
        "outlet_revision": outlet_revision(net),
        "t0_lakes": {"rule": "per T0 lake body with lake nodes in the chart: plain 512 m node area of T1b land nodes whose water passes the lake (path ends in it, or on a T1a seed whose T1a receivers pass it) vs the same nodes whose T0 cell chain passes a lake cell of that body; the T1b lake stage is the T0 stage, so a change > 100 permille is flagged (reported, not gated)", "flagged": lakes_flagged, "lakes": lake_rows},
        "displacement": {"land_nodes": land_n, "nodes_at_other_outlet": wrong_n, "max_edge_distance_km_of_those": max_d / 1000.0,
            "violations_at_bound_gated_cells": bound_viol, "violations_at_bound_rim_cells": bound_viol_rim, "violations_by_gated_cell": viol_cells,
            "table": xs_km.iter().zip(viol.iter()).map(|(x, v)| serde_json::json!({"x_km": x, "nodes": v.0, "area_km2": v.1 / 1e6})).collect::<Vec<_>>()},
    }))
}

/// Land own areas as build_network sets them: the 512 m node area, times R / r_ref with the climate key (T0 with climate
/// attached; the binary attaches it only with climate=).
fn own_exact(t0: &T0, chart: &Chart) -> R<Vec<f64>> {
    let a_node = (SPACING_MM as f64 / 1000.0).powi(2);
    let mut out = vec![a_node; NODES];
    if let Some(c) = t0.climate.as_ref() {
        for (n, o) in out.iter_mut().enumerate() {
            let (gx, gy) = node_g(n);
            let (la, lo) = chart.inverse_mm(gx * SPACING_MM, gy * SPACING_MM).map_err(|e| e.to_string())?;
            *o = a_node * c.weight_at_deg(t0, la * vz_tilekey::chart::RAD2DEG, lo * vz_tilekey::chart::RAD2DEG);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn terrain_route_requires_coupled_flux_and_global_solve() {
        let mut p = RiverNetParamsV1::v1(true);
        parse_key("conf=terrain", &mut p).unwrap();
        assert!(validate_terrain(&p).is_err());
        let c = p.confinement.as_ref().unwrap();
        assert_eq!((c.mode, c.flux, c.levels, c.bound_mm, c.core), (CONF_TERRAIN, 1, 1, 0, 0));
        // A missing level solve must not silently become a random-walk-only prototype.
        parse_key("ls_global=1", &mut p).unwrap();
        parse_key("ls_settle=1", &mut p).unwrap();
        assert!(validate_terrain(&p).is_err());
        p.ls = Some(crate::t1b_ls::LsParamsV1::v1());
        assert!(validate_terrain(&p).is_ok());
        parse_key("conf_levels=0", &mut p).unwrap();
        assert!(validate_terrain(&p).is_err());
        parse_key("conf_levels=1", &mut p).unwrap();
        parse_key("ls_global=0", &mut p).unwrap();
        assert!(validate_terrain(&p).is_err());
    }

    #[test]
    fn keys_off_leave_the_recipe_params_unchanged() {
        // Absent confinement and the checker switch serialize to exactly today's params JSON.
        let p = RiverNetParamsV1::v1(true);
        let mut q = p.clone();
        q.outlet_checks = 1;
        assert_eq!(serde_json::to_string(&p).unwrap(), serde_json::to_string(&q).unwrap());
        assert!(serde_json::to_value(&p).unwrap().get("confinement").is_none());
        assert!(serde_json::to_value(&p).unwrap().get("outlet_checks").is_none());
    }

    #[test]
    fn parse_key_sets_modes_and_refuses_bad_values() {
        let mut p = RiverNetParamsV1::v1(true);
        assert!(!parse_key("a_ch=5", &mut p).unwrap());
        assert!(parse_key("conf=basin", &mut p).unwrap());
        assert_eq!(p.confinement.as_ref().unwrap().mode, CONF_FINAL_BASIN);
        assert!(parse_key("conf_warp_mm=1000", &mut p).is_err(), "warp amplitude needs conf=warp first");
        assert!(parse_key("conf=warp", &mut p).unwrap());
        assert!(parse_key("conf_warp_mm=4096000", &mut p).unwrap());
        assert_eq!(p.confinement.as_ref().unwrap().warp_amplitude_mm, 4_096_000);
        assert!(parse_key("conf=sideways", &mut p).is_err());
        assert!(parse_key("outlet_checks=1", &mut p).unwrap());
        assert_eq!(p.outlet_checks, 1);
        // Non-warp modes carry no warp fields in the recipe.
        let v = serde_json::to_value(ConfinementParamsV1::v1(CONF_NONE).unwrap()).unwrap();
        assert!(v.get("warp_amplitude_mm").is_none());
    }

    #[test]
    fn warp_cells_are_keyed_and_zero_amplitude_is_identity() {
        let cells: Vec<u32> = (0..NODES).map(|n| { let (gx, gy) = node_g(n); ((gx + 256) / 64 * 100 + (gy + 256) / 64) as u32 }).collect();
        let mut cp = ConfinementParamsV1::v1(CONF_WARP).unwrap();
        let a = warp_cells(&cp, &cells, 42, "blake3:test").unwrap();
        assert_eq!(a, warp_cells(&cp, &cells, 42, "blake3:test").unwrap());
        assert_ne!(a, cells);
        assert_ne!(a, warp_cells(&cp, &cells, 43, "blake3:test").unwrap());
        cp.warp_amplitude_mm = 0;
        assert_eq!(warp_cells(&cp, &cells, 42, "blake3:test").unwrap(), cells);
    }

    fn toy_net() -> Network {
        let lat = |k: usize| { let (gx, gy) = node_g(k); [(gx * SPACING_MM) as f64 / 1000.0, (gy * SPACING_MM) as f64 / 1000.0] };
        Network { kind: vec![K_MARINE; NODES], recv: vec![NONE; NODES], order: vec![], area_m2: vec![0.0; NODES], pos: (0..NODES).map(lat).collect(),
            mid: vec![[0.0; 2]; NODES], term: vec![NONE; NODES], level: vec![0.0; NODES], unreached_confined: 0, nodes_land: 0, rockmass: vec![1.0; NODES],
            conduit: vec![false; NODES], conduit_cap: vec![f64::NAN; NODES], conduit_paths: vec![], t1a_recv: vec![NONE; NODES], t0_cell: vec![0; NODES],
            meander_u: vec![0.5; NODES], t1a_ch: vec![false; NODES], conf_cell: vec![], flux: None }
    }

    #[test]
    fn parse_key_free_mode_and_bound() {
        let mut p = RiverNetParamsV1::v1(true);
        assert!(parse_key("conf_bound_mm=4000000", &mut p).is_err(), "bound needs conf=free first");
        assert!(parse_key("conf=free", &mut p).unwrap());
        let c = p.confinement.as_ref().unwrap();
        assert_eq!((c.mode, c.bound_mm, c.flux), (CONF_FREE, FREE_BOUND_DEFAULT_MM, 1));
        assert!(parse_key("conf_bound_mm=4000000", &mut p).unwrap());
        assert_eq!(p.confinement.as_ref().unwrap().bound_mm, 4_000_000);
        assert!(parse_key("conf_flux=0", &mut p).is_err(), "free always carries the flux");
        let mut q = RiverNetParamsV1::v1(true);
        assert!(parse_key("conf=none", &mut q).unwrap());
        // Old modes keep their recipe JSON (no new fields) until conf_flux=1.
        assert_eq!(serde_json::to_value(q.confinement.as_ref().unwrap()).unwrap(), serde_json::to_value(ConfinementParamsV1 { mode: CONF_NONE, rule: q.confinement.as_ref().unwrap().rule.clone(), warp_amplitude_mm: 0, warp_spacing_mm: 0, bound_mm: 0, flux: 0, flux_rule: String::new(), levels: 0, levels_rule: String::new(), core: 0, core_rule: String::new() }).unwrap());
        assert!(serde_json::to_value(q.confinement.as_ref().unwrap()).unwrap().get("flux").is_none());
        assert!(parse_key("conf_flux=1", &mut q).unwrap());
        assert_eq!(q.confinement.as_ref().unwrap().flux, 1);
    }

    #[test]
    fn flux_moves_area_between_seeds_and_carries_it_down() {
        // Two T1a chains: a1 -> a2 (outlet A) and b1 (outlet B). Land u: T1a sends it to a1, T1b sends it to b1.
        // Land v: both send it to a1. Lake node k between: land w drains into the lake in T1b, the lake's T1a chain reaches b1.
        let ix = |gx: i64, gy: i64| node_index(gx, gy).unwrap();
        let (a1, a2, b1, u, v, k, w) = (ix(0, 0), ix(1, 0), ix(0, 5), ix(0, 2), ix(1, 1), ix(0, 4), ix(1, 4));
        let mut net = toy_net();
        for t in [a1, a2, b1] { net.kind[t] = K_T1A; }
        net.t1a_recv[a1] = a2 as u32;
        net.kind[k] = K_LAKE; net.t1a_recv[k] = b1 as u32;
        for l in [u, v, w] { net.kind[l] = K_LAND; }
        net.t1a_recv[u] = a1 as u32; net.recv[u] = b1 as u32;
        net.t1a_recv[v] = a1 as u32; net.recv[v] = a1 as u32;
        net.t1a_recv[w] = a1 as u32; net.recv[w] = k as u32;
        let mut a0 = vec![0.0; NODES];
        a0[a1] = 10.0; a0[a2] = 12.0; a0[b1] = 5.0;
        net.flux = Some(FluxBase { t1a_area0: a0, bound_mm: 0, outlet_basin: vec![NONE; NODES], t0_q_cell: vec![0.0; NODES], r_ref_m_per_yr: 1.0, own: vec![], levels_t1a: false });
        let mut own = vec![0.0; NODES];
        own[u] = 1.0; own[v] = 2.0; own[w] = 4.0;
        let (s, d) = flux_seed_areas(&net, &own);
        // a1 loses u (1) and w (4); a2 inherits that; b1 gains u and w (w through the lake).
        assert_eq!((s[a1], s[a2], s[b1]), (5.0, 7.0, 10.0));
        assert_eq!(d.moved_area_m2, 5.0);
        assert_eq!(d.unknown_destination_nodes, 0);
        // Conservation over the exits: (a2 + b1) unchanged.
        assert_eq!(s[a2] + s[b1], 12.0 + 5.0);
    }

    #[test]
    fn free_rule_core_nodes_stay_home() {
        // Mode 4 zones: core node c0 (zone = cell 7) may drain to core land of cell 7 or to a T1a node of its own outlet only.
        let mut p = RiverNetParamsV1::v1(false);
        p.confinement = Some(ConfinementParamsV1::v1(CONF_FREE).unwrap());
        let ix = |gx: i64, gy: i64| node_index(gx, gy).unwrap();
        let (c0, c1, band, t_own) = (ix(0, 0), ix(1, 0), ix(2, 0), ix(0, 1));
        let mut zone = vec![NONE; NODES];
        zone[c0] = 7; zone[c1] = 7; zone[t_own] = 70; zone[ix(5, 5)] = 71;
        let conf = Confinement { mode: CONF_FREE, cell: vec![7; NODES], fb: [(7u32, 70u32)].into_iter().collect(), seed_outlet: vec![NONE; NODES], zone, path: BTreeMap::new() };
        assert!(conf.zone_allows(c0, c1, false));
        assert!(!conf.zone_allows(c0, band, false), "core may not drain into the band");
        assert!(conf.zone_allows(band, c0, false) && conf.zone_allows(band, ix(3, 0), false), "band nodes are free");
        assert!(conf.zone_allows(c0, t_own, true));
        assert!(!conf.zone_allows(c0, ix(5, 5), true), "a T1a node of another outlet");
    }

    #[test]
    fn conf_levels_key_needs_the_flux_and_keeps_recipes_without_it() {
        // Lane DRAIN-FIX: conf_levels=1 only after the flux is on; without it conf=free's recipe JSON is unchanged.
        let mut p = RiverNetParamsV1::v1(true);
        assert!(parse_key("conf_levels=1", &mut p).is_err(), "needs conf= first");
        assert!(parse_key("conf=cell", &mut p).unwrap());
        assert!(parse_key("conf_levels=1", &mut p).is_err(), "needs the flux on");
        let mut f = RiverNetParamsV1::v1(true);
        assert!(parse_key("conf=free", &mut f).unwrap());
        let before = serde_json::to_value(f.confinement.as_ref().unwrap()).unwrap();
        assert!(before.get("levels").is_none() && before.get("levels_rule").is_none());
        assert!(parse_key("conf_levels=1", &mut f).unwrap());
        let c = f.confinement.as_ref().unwrap();
        assert_eq!((c.levels, c.levels_rule.as_str()), (1, LEVELS_RULE));
        assert!(parse_key("conf_levels=2", &mut f).is_err());
        assert!(parse_key("conf=none", &mut p).unwrap() && parse_key("conf_flux=1", &mut p).unwrap() && parse_key("conf_levels=1", &mut p).unwrap());
        assert!(parse_key("conf_flux=0", &mut p).is_err(), "the flux cannot be turned off under conf_levels=1");
        // wall_checks is a checker switch: never in the recipe.
        let mut q = RiverNetParamsV1::v1(true);
        let j0 = serde_json::to_string(&q).unwrap();
        assert!(parse_key("wall_checks=1", &mut q).unwrap());
        assert_eq!(q.wall_checks, 1);
        assert_eq!(serde_json::to_string(&q).unwrap(), j0);
    }

    #[test]
    fn conf_core_path_key_and_destinations() {
        // Lane DRAIN-FIX: conf_core=path needs conf=free; without it the recipe has no core fields.
        let mut p = RiverNetParamsV1::v1(true);
        assert!(parse_key("conf_core=path", &mut p).is_err());
        assert!(parse_key("conf=none", &mut p).unwrap());
        assert!(parse_key("conf_core=path", &mut p).is_err(), "free only");
        let mut f = RiverNetParamsV1::v1(true);
        assert!(parse_key("conf=free", &mut f).unwrap());
        assert!(serde_json::to_value(f.confinement.as_ref().unwrap()).unwrap().get("core").is_none());
        assert!(parse_key("conf_core=path", &mut f).unwrap());
        let c = f.confinement.as_ref().unwrap();
        assert_eq!((c.core, c.core_rule.as_str()), (1, CORE_PATH_RULE));
        assert!(parse_key("conf_core=x", &mut f).is_err());
        // destinations: land chain u -> v -> T1a t; w without a receiver; t itself.
        let ix = |gx: i64, gy: i64| node_index(gx, gy).unwrap();
        let (t, u, v, w) = (ix(0, 0), ix(2, 0), ix(1, 0), ix(5, 5));
        let mut net = toy_net();
        net.kind[t] = K_T1A;
        for l in [u, v, w] { net.kind[l] = K_LAND; }
        net.recv[u] = v as u32; net.recv[v] = t as u32;
        let d = destinations(&net);
        assert_eq!((d[u], d[v], d[t], d[w]), (t as u32, t as u32, t as u32, NONE));
    }

    #[test]
    fn level_areas_keep_the_t1a_discharge_for_levels_and_conduits() {
        // Lane DRAIN-FIX: seed t (T1a area 10, flux area 3) injects into conduit c0 -> c1 -> sea; level areas put the seed back
        // at 10 and shift the conduit chain by +7; without levels_t1a there are no level areas.
        let ix = |gx: i64, gy: i64| node_index(gx, gy).unwrap();
        let (t, c0, c1, sea) = (ix(0, 0), ix(1, 0), ix(2, 0), ix(3, 0));
        let mut net = toy_net();
        net.kind[t] = K_T1A;
        for c in [c0, c1] { net.kind[c] = K_LAND; net.conduit[c] = true; }
        net.recv[c0] = c1 as u32; net.recv[c1] = sea as u32;
        net.conduit_paths = vec![vec![t as u32, c0 as u32, c1 as u32, sea as u32]];
        net.area_m2[t] = 3.0; net.area_m2[c0] = 4.0; net.area_m2[c1] = 5.0;
        let mut a0 = vec![0.0; NODES];
        a0[t] = 10.0;
        net.flux = Some(FluxBase { t1a_area0: a0, bound_mm: 0, outlet_basin: vec![NONE; NODES], t0_q_cell: vec![0.0; NODES], r_ref_m_per_yr: 1.0, own: vec![], levels_t1a: false });
        assert!(level_areas(&net).is_none());
        net.flux.as_mut().unwrap().levels_t1a = true;
        let a = level_areas(&net).unwrap();
        assert_eq!((a[t], a[c0], a[c1], a[sea]), (10.0, 11.0, 12.0, 0.0));
        assert_eq!(net.area_m2[t], 3.0, "the carried discharge is not touched");
    }
}
