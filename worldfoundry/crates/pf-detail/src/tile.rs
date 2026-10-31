//! Per-tile pipeline (B1-B10) and block assembly.

use crate::channels::{ChannelHit, Channels};
use crate::input::{T0Strata, T1a, T1bWindow};
use crate::relax::{relax, RelaxStats};
use crate::strata::{CellStrata, StrataModel};
use crate::v2::{Aux, Curves, V2Params};
use crate::*;
use std::cell::RefCell;
use std::collections::{BTreeMap, BinaryHeap, HashMap};
use std::time::Instant;
use vz_keyed_rng::Scope;

pub struct Ctx<'a> {
    pub escarp2: Option<crate::escarpment::Reconstructed>,
    pub base_c1: Option<crate::base_interp::C1Surface>,
    pub w: &'a T1bWindow,
    pub t0: &'a T0Strata,
    pub p: DetailParamsV1,
    pub strata: StrataModel,
    pub ch: Channels,
    pub layer_scope: Scope,
    pub sea: f64,
    // Lane PF-SERVE: the memos are Mutex (were RefCell) so one context can serve tiles from several threads; every entry is
    // a pure function of its key, so values do not depend on which thread filled them.
    relief_cache: std::sync::Mutex<HashMap<(i64, i64), f64>>,
    relax_cache: std::sync::Mutex<HashMap<(i64, i64), (Vec<f64>, f64)>>,
    rough_cache: std::sync::Mutex<HashMap<(i64, i64), f64>>,
    /// Lane DE: per tile domain, the DE change eta (n x n) and the solve stats.
    /// (eta, talus: deposit m > 0, mantle -1, none 0, solve stats)
    de_cache: std::sync::Mutex<HashMap<(i64, i64), (Vec<f64>, Vec<f64>, serde_json::Value)>>,
    /// Lane TREAD (ct_sharp_stage=pre): per tile domain, the pre-stage shift (mm), 2 m regime and riser rockfall stats.
    sharp_cache: std::sync::Mutex<HashMap<(i64, i64), (Vec<i32>, Vec<u8>, serde_json::Value)>>,
    pub v2: Option<V2Ctx>,
    /// Lane PV (set after `new` when `p.prov` is on): provenance sources.
    pub prov: Option<DetailProv>,
    /// Lane RILL (de_b5_chain): the coarse window-wide B5 pass, built on first use.
    chain: std::sync::OnceLock<crate::chain::Coarse>,
    /// Lane SHARP (ct_sharp=1): the 2 m re-profile of T1b's cliff-and-talus staircase (whole window: block-independent).
    pub sharp: Option<crate::sharp::SharpField>,
}

/// Lane PV: the provenance sources of a block run.
pub struct DetailProv {
    pub strata: std::sync::Arc<crate::strata_prov::StrataProv>,
    /// T0 volcano records by id (onset, last eruption, active) and their digest (None: not given; onset stays null).
    pub volc: Option<(BTreeMap<String, (Option<f64>, Option<f64>, Option<bool>)>, String)>,
}

/// v2 context (lane T2B): curves of the whole T1b window, T1b-lattice aux fields over the bbox, the edifice column.
pub struct V2Ctx {
    pub v: V2Params,
    pub curves: Curves,
    pub aux: Aux,
    pub ed: CellStrata,
    /// v3 W2 lateral pinch-and-swell amplitudes (None unless version=3 with pinch on).
    pub pinch: Option<crate::v3::Pinch>,
    /// Lane BS (beds=2): one edifice bed sequence per volcano. None = the v2/v3 single sequence `ed`.
    pub ed_bs: Option<EdBs>,
    /// Lane DE (de_lobes=1): edifice units with lateral extent per volcano. None = off.
    pub lobes: Option<crate::lobes::LobeWorld>,
    /// Lane PYRO (pyro_sub=1): pyroclastic records of the window's volcanoes. None = off.
    pub pyro: Option<crate::pyro::PyroWorld>,
}

/// Lane BS: per-volcano edifice columns and the volcano index of every aux (T1b) post.
pub struct EdBs {
    /// Volcano index (T1b header list) per aux post; 255 = not identified. Meaningful on edifice posts only.
    pub vid: Vec<u8>,
    pub cols: BTreeMap<u8, CellStrata>,
    /// Per volcano used: the record and the style (block.json).
    pub info: Vec<serde_json::Value>,
    /// "t1b_volcano_id" or "nearest_vent_fallback" (the window publishes no volcano_id).
    pub source: &'static str,
}

/// Lane BS2: index of the generated unit containing signed bed-coordinate depth dw (None above the first unit top or
/// without unit tops).
pub fn unit_at(col: &crate::strata::Column, dw: f64) -> Option<usize> {
    let i = col.unit_top.partition_point(|q| q.0 <= dw);
    if i == 0 { None } else { Some(i - 1) }
}

/// Lane BS: pinch column id of volcano v's edifice sequence (disjoint from cell ids < 2^20 and EDIFICE_CELL).
#[inline]
pub fn ed_col_id(v: u8) -> u32 { 0xF000_0000 | v as u32 }

/// Lane BS: the layer records bed styles read (strata layer records; T0 volcano records by id: onset, last eruption, active).
pub struct BsRecords {
    pub strata: std::sync::Arc<crate::strata_prov::StrataProv>,
    pub volc: Option<BTreeMap<String, (Option<f64>, Option<f64>, Option<bool>)>>,
}

/// Lane DE: a post's column for the DE solve: (column, datum, pinch amplitudes) or a constant (Sc, K).
enum DeCol<'b> { Col(&'b CellStrata, f64, Option<[f64; 2]>), Const(f64), Lobe(usize) }

/// Pointwise data of one post (B1, B2, B3, B6 lookup).
#[derive(Clone, Copy, Debug)]
pub struct PostVal {
    pub base: f64,
    pub cut: f64,
    pub cell: u32,
    pub fixed_sea: bool,
    pub hit: ChannelHit,
}

impl PostVal {
    pub fn fixed_relax(&self) -> bool { self.fixed_sea || self.hit.corridor }
}

pub fn tile_origin_post(ix: i64, iy: i64) -> (i64, i64) { (ix * TILE_CELLS, iy * TILE_CELLS) }

/// Lane BS: per-volcano edifice columns. The volcano of an aux edifice post is the T1b `volcano_id` of the nearest T1b post
/// that has one (within 2 posts); a window without `volcano_id` falls back to the nearest listed vent within its r_influence_m
/// (labelled `nearest_vent_fallback`). Each column is keyed by the volcano's id string, so it is the same in every block.
fn build_ed_bs(w: &T1bWindow, aux: &Aux, r: &BsRecords, bp: &crate::bedstyle::BedStyleParams, scope: &Scope, p: &DetailParamsV1) -> R<EdBs> {
    let np = aux.nx * aux.ny;
    let mut vid = vec![255u8; np];
    let source = if w.volcano_id.is_some() { "t1b_volcano_id" } else { "nearest_vent_fallback" };
    for k in 0..np {
        if aux.over[k] != crate::v2::OVER_EDIFICE { continue; }
        let (gx, gy) = (aux.gx0 + (k % aux.nx) as i64, aux.gy0 + (k / aux.nx) as i64);
        vid[k] = match &w.volcano_id {
            Some(f) => {
                let mut best: Option<(i64, u8)> = None;
                for ring in 0..=2i64 {
                    for dy in -ring..=ring { for dx in -ring..=ring {
                        if dx.abs().max(dy.abs()) != ring { continue; }
                        if let Some(i) = w.idx(gx + dx, gy + dy) { if f[i] != 255 { let d2 = dx * dx + dy * dy; if best.map(|b| d2 < b.0).unwrap_or(true) { best = Some((d2, f[i])); } } }
                    } }
                    if best.is_some() { break; }
                }
                best.map(|b| b.1).unwrap_or(255)
            }
            None => {
                let (x, y) = (gx as f64 * T1B_M, gy as f64 * T1B_M);
                let mut best: Option<(f64, u8)> = None;
                for (v, g) in w.volcano_geo.iter().enumerate() {
                    let d = (g[0] - x).hypot(g[1] - y);
                    if d.is_finite() && d <= g[2] && best.map(|b| d < b.0).unwrap_or(true) && v < 255 { best = Some((d, v as u8)); }
                }
                best.map(|b| b.1).unwrap_or(255)
            }
        };
    }
    let mut used: Vec<u8> = vid.iter().zip(aux.over.iter()).filter(|(_, &o)| o == crate::v2::OVER_EDIFICE).map(|(&v, _)| v).collect();
    used.sort(); used.dedup();
    let depth = p.v2.as_ref().map(|v| v.edifice_max_depth_m).unwrap_or(4000.0);
    let mut cols = BTreeMap::new();
    let mut info = Vec::new();
    for v in used {
        let vi = w.volcanoes.get(v as usize);
        let id = vi.map(|i| i.id.clone()).unwrap_or_else(|| "unknown".into());
        let rec = r.volc.as_ref().and_then(|m| m.get(&id)).cloned();
        let vr = crate::bedstyle::VolcRec { id: id.clone(), kind: vi.map(|i| i.kind.clone()).unwrap_or_default(),
            height_m: w.volcano_geo.get(v as usize).map(|g| g[3]).filter(|h| h.is_finite()),
            onset_my: rec.and_then(|q| q.0), last_eruption_my: rec.and_then(|q| q.1) };
        let tag = format!("edifice/{id}");
        let (o, style) = crate::bedstyle::edifice_beds(&vr, depth, scope, &tag, bp)?;
        // the sequence continued above the datum: bed tops above the datum, lowest bed first (its base at the datum)
        let u = crate::bedstyle::edifice_up_beds(&vr, 3000.0, scope, &format!("{tag}/up"), bp)?;
        let n_down = o.bb.len();
        let (bb, mut bs, mut bk, mut bf, mut bu) = (o.bb, o.bs, o.bk, o.bf, o.bu);
        bs.extend(u.bs); bk.extend(u.bk); bf.extend(u.bf); bu.extend(u.bu);
        let mut uts = u.ut; uts.extend(o.ut);
        uts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let col = crate::strata::Column { lith: vec![crate::v2::LITHO_VOLCANIC], base_depth: vec![depth], bed_base: bb, bed_sc: bs, bed_k: bk, up_top: u.bb, n_down,
            bed_fac: bf, bed_unit: bu, unit_top: uts };
        let mut cs = CellStrata { cell: ed_col_id(v), top_m: 0.0, g: (0.0, 0.0), g_raw: (0.0, 0.0), col, u0: 0.0, du: p.remap_table_step_m, psi: Vec::new(), max_drift_m: 0.0, sc_win: Vec::new(), k_win: Vec::new(), cum_is: Vec::new() };
        crate::strata::build_remap_pub(&mut cs, p)?;
        info.push(serde_json::json!({"volcano_index": v, "volcano": id, "kind": vr.kind, "style": style, "height_m": vr.height_m, "onset_my": vr.onset_my,
            "last_eruption_my": vr.last_eruption_my, "p_interflow": crate::bedstyle::edifice_p_interflow(&vr, bp), "beds_down": n_down,
            "aux_posts": vid.iter().zip(aux.over.iter()).filter(|(&q, &o)| q == v && o == crate::v2::OVER_EDIFICE).count()}));
        cols.insert(v, cs);
    }
    Ok(EdBs { vid, cols, info, source })
}

impl<'a> Ctx<'a> {
    pub fn new(w: &'a T1bWindow, t1a: &'a T1a, t0: &'a T0Strata, p: DetailParamsV1, anchor: (f64, f64), layer_digest: &str,
               world_seed: u64, chart: u32, tiles: &[(i64, i64)]) -> R<Ctx<'a>> {
        Self::new_with(w, t1a, t0, p, anchor, layer_digest, world_seed, chart, tiles, None)
    }

    /// `bs`: the records beds=2 needs (required when p.beds is set).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with(w: &'a T1bWindow, t1a: &'a T1a, t0: &'a T0Strata, p: DetailParamsV1, anchor: (f64, f64), layer_digest: &str,
               world_seed: u64, chart: u32, tiles: &[(i64, i64)], bs: Option<&BsRecords>) -> R<Ctx<'a>> {
        if w.escarpment.is_some() && p.base_interp.is_some() {return Err("ESCARP-1 carries fine geometry; base_interp=c1 would reconstruct coarse means as points".into());}
        let escarp2 = crate::escarpment::Reconstructed::from_window(w)?;
        if let Some(fine) = &escarp2 {
            if fine.has_toes() && p.base_interp.is_none() {return Err("ct_escarp=6 requires base_interp=c1".into());}
            if w.escarpment.is_some() || (p.base_interp.is_some() && !fine.is_antialiased()) {
                return Err("base_interp=c1 composes only with ESCARP-4, outside its conserved band".into());
            }
        }
        let base_c1 = if p.base_interp.is_some() {
            if w.evolved.is_none() { return Err("base_interp=c1 currently requires an E5B process-evolved parent; legacy construction is not qualified".into()); }
            Some(crate::base_interp::C1Surface::new(w)?)
        } else { None };
        if p.beds.is_some() && bs.is_none() { return Err("beds=2 needs the layer records".into()); }
        if p.climate_runoff.is_some() && w.runoff_weight.is_none() { return Err("climate_runoff=1 needs a T1b window built with pf-region climate=<dir> (runoff_weight field)".into()); }
        // bbox over the tiles plus two tiles of margin (neighbour relaxations and their halos)
        let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
        for &(ix, iy) in tiles {
            x0 = x0.min((ix - 2) as f64 * TILE_M); y0 = y0.min((iy - 2) as f64 * TILE_M);
            x1 = x1.max((ix + 3) as f64 * TILE_M); y1 = y1.max((iy + 3) as f64 * TILE_M);
        }
        // lane RILL de_b5_chain: the coarse pass covers the whole T1b window, so the context does too (every block the same)
        if p.de.as_ref().and_then(|d| d.chain.as_ref()).is_some() {
            x0 = w.gx0 as f64 * T1B_M; y0 = w.gy0 as f64 * T1B_M;
            x1 = (w.gx0 + w.n as i64 - 1) as f64 * T1B_M; y1 = (w.gy0 + w.n as i64 - 1) as f64 * T1B_M;
        }
        let layer_scope = Scope::root(world_seed).layer(layer_digest).map_err(|e| e.0)?.chart(chart).map_err(|e| e.0)?;
        let mut strata = StrataModel::new(anchor);
        if p.beds.is_some() { strata.bs = bs.map(|b| b.strata.clone()); }
        let mut gy = (y0 / T1B_M).floor() as i64;
        while gy as f64 * T1B_M <= y1 {
            let mut gx = (x0 / T1B_M).floor() as i64;
            while gx as f64 * T1B_M <= x1 {
                if let Some(i) = w.idx(gx, gy) {
                    let c = w.t0_cell[i];
                    if w.fixed[i] != 2 { let _ = strata.ensure(t0, c, &p, &layer_scope); }
                }
                gx += 1;
            }
            gy += 1;
        }
        let ch = Channels::build(w, t1a, &p, (x0, y0, x1, y1));
        let v2 = match p.v2.clone() {
            None => None,
            Some(v) => {
                // the window's own cells too (over-codes compare the T1b litho with our cut on every aux post)
                let m = 8i64;
                let (gx0, gy0) = ((x0 / T1B_M).floor() as i64 - m, (y0 / T1B_M).floor() as i64 - m);
                let (gx1, gy1) = ((x1 / T1B_M).ceil() as i64 + m, (y1 / T1B_M).ceil() as i64 + m);
                for gy in gy0..=gy1 { for gx in gx0..=gx1 {
                    if let Some(i) = w.idx(gx, gy) { if w.fixed[i] != 2 { let _ = strata.ensure(t0, w.t0_cell[i], &p, &layer_scope); } }
                } }
                let curves = Curves::build(w, t1a, &p, &v, (x0, y0, x1, y1));
                // lane DATUM ed_datum=1: the edifice datum from the pre-wave-cut T1b surface (when the window publishes it)
                let pre = if p.beds.as_ref().and_then(|b| b.ed_datum.as_ref()).is_some() { w.pre_cut_m.as_deref() } else { None };
                // lane ARCH ed_env: the datum from the interfluve envelope of that surface (ED_ENV_NOTE); routing too with ed_env_route
                let env_p = p.beds.as_ref().and_then(|b| b.ed_env.as_ref());
                let env: Option<Vec<f64>> = env_p.map(|e| crate::bedstyle::envelope_strike(pre.unwrap_or(&w.height_m), w.n, e.half_posts, e.strike_half_posts));
                let pre_route = match (&env, env_p.and_then(|e| e.route)) { (Some(v), Some(true)) => Some(v.as_slice()), _ => pre };
                let pre = match &env { Some(v) => Some(v.as_slice()), None => pre };
                let aux = Aux::build(w, &strata, &curves, &v, gx0, gy0, gx1, gy1, pre);
                let ed = crate::v2::edifice_strata(&layer_scope, &p, &v)?;
                let pinch = v.v3.as_ref().filter(|w| w.pinch_on).map(crate::v3::Pinch::new);
                let ed_bs = match (p.beds.as_ref(), bs) { (Some(bp), Some(r)) => Some(build_ed_bs(w, &aux, r, bp, &layer_scope, &p)?), _ => None };
                let pyro = p.pyro.as_ref().map(|pp| crate::pyro::PyroWorld::build(pp, w));
                let lobes = match (p.de.as_ref().and_then(|d| d.lobes.as_ref()), ed_bs.as_ref()) {
                    (Some(lp), Some(eb)) => Some(crate::lobes::LobeWorld::build(lp, w, &eb.cols.keys().copied().collect::<Vec<u8>>(), (x0, y0, x1, y1), &layer_scope, pre_route, pyro.as_ref())?),
                    _ => None,
                };
                Some(V2Ctx { v, curves, aux, ed, pinch, ed_bs, lobes, pyro })
            }
        };
        if p.de.as_ref().and_then(|d| d.regime.as_ref()).is_some() && p.ct_sharp.is_none() { return Err("de_talus_mantle / de_talus_runoff need ct_sharp=1 ct_sharp_stage=pre".into()); }
        let sharp = match &p.ct_sharp {
            None => None,
            Some(sp) => {
                if sp.pre() && !(p.relax_on && p.v2.is_some()) { return Err("ct_sharp_stage=pre needs the B5 relaxation and version=2|3".into()); }
                if sp.cell_riser.is_some() && (sp.ends.is_none() || sp.blend.as_deref() != Some(crate::sharp::SHARP_BLEND_PHASE_RULE)) { return Err("ct_sharp_cell_riser needs ct_sharp_blend=phase and ct_sharp_ends=1".into()); }
                if sp.talus_rock.is_some() && !sp.pre() { return Err("ct_sharp_talus_rock needs ct_sharp_stage=pre (the domain surface the 2 m regime labels)".into()); }
                if sp.riser_noflux.is_some() && !sp.pre() { return Err("ct_sharp_riser_noflux needs ct_sharp_stage=pre (the held risers in B5)".into()); }
                if sp.aa.is_some() && !sp.pre() { return Err("ct_sharp_aa needs ct_sharp_stage=pre (the pre-B5 staircase)".into()); }
                if (sp.carry_strike.is_some() || sp.limb_lp.is_some()) && !sp.pre() { return Err("ct_sharp_carry_strike / ct_sharp_limb_lp need ct_sharp_stage=pre (the held-riser carry and the pre-B5 staircase)".into()); }
                if p.de.as_ref().and_then(|d| d.regime.as_ref()).is_some() && !sp.pre() { return Err("de_talus_mantle / de_talus_runoff need ct_sharp_stage=pre (the 2 m regime of the pre staircase)".into()); }
                let ct = w.ct.as_ref().ok_or("ct_sharp=1 needs a T1b window with ct_base_mm, ct_regime and ct_cliff_permille (pf-region ct=1 with ct_base=1)")?;
                let lit = w.exposed_litho.as_ref().ok_or("ct_sharp=1 needs the window's exposed_litho")?;
                let tanp = |k: &serde_json::Value, f: &str| k.get(f).and_then(serde_json::Value::as_i64).map(|v| v as f64 / 1000.0).unwrap_or(0.0);
                let classes: Vec<crate::sharp::CtClass> = ct.params.get("classes").and_then(serde_json::Value::as_array).ok_or("ct_sharp: the window header has no ct class table")?
                    .iter().map(|k| crate::sharp::CtClass { mode: k.get("mode").and_then(serde_json::Value::as_u64).unwrap_or(0) as u32, tb: tanp(k, "talus_tan_permille"), ta: tanp(k, "cliff_tan_permille") }).collect();
                let step_m = ct.params.get("step_mm").and_then(serde_json::Value::as_i64).ok_or("ct_sharp: the window header has no ct step_mm")? as f64 / 1000.0;
                // the rock-mass multiplier scales the cliff tangent only with rockmass_cliff=1 (pf-region t1b_ct::split)
                let rockmass = if ct.params.get("rockmass_cliff").and_then(serde_json::Value::as_u64) == Some(0) { (1.0, 1.0) } else { ct.rockmass };
                let ci = crate::sharp::CtInputs { n: w.n, gx0: w.gx0, gy0: w.gy0, h: &ct.base_m, z: &w.height_m, regime: &ct.regime, cliff_permille: &ct.cliff_permille, litho: lit, classes, step_m, rockmass,
                    tread: ct.tread.as_ref().map(|(a, b)| (a.as_slice(), b.as_slice())),
                    ends: if sp.ends.is_some() { Some(ct.ends.as_ref().ok_or("ct_sharp_ends=1 needs a T1b window with pf-region ct_ends's fields (ct_s0_mm, ct_s1_mm, ct_s_dn_mm, ct_s_up_mm, ct_knot_*_mm, ct_grad_permille, ct_period_kind)")?) } else { None } };
                let mut sf = crate::sharp::SharpField::build(&ci, sp);
                // lane TREADSOIL: valley-floor / shore soil (on or next to a fixed T1b post) keeps the linear fade
                sf.exclude_soil_near_fixed(&w.fixed);
                Some(sf)
            }
        };
        Ok(Ctx { escarp2, base_c1, w, t0, sea: t0.sea_level_m, p, strata, ch, layer_scope, relief_cache: std::sync::Mutex::new(HashMap::new()), relax_cache: std::sync::Mutex::new(HashMap::new()), rough_cache: std::sync::Mutex::new(HashMap::new()), de_cache: std::sync::Mutex::new(HashMap::new()), sharp_cache: std::sync::Mutex::new(HashMap::new()), v2, prov: None, chain: std::sync::OnceLock::new(), sharp })
    }

    pub fn cell_at(&self, x: f64, y: f64) -> u32 {
        self.w.nearest(x, y).map(|i| self.w.t0_cell[i]).unwrap_or(u32::MAX)
    }

    pub fn post(&self, px: i64, py: i64) -> PostVal {
        let (x, y) = (px as f64 * POST_M, py as f64 * POST_M);
        let base = match &self.v2 { Some(v2) if v2.v.deperch_on => v2.aux.bicubic(&v2.aux.zbase, x, y), _ => self.w.bicubic(x, y) };
        let near_i = self.w.nearest(x, y);
        let code = near_i.map(|i| self.w.fixed[i]).unwrap_or(2);
        let cell = self.cell_at(x, y);
        // lane INTS2n4 arm_cl (chn.rs ARM_CL_RULE): an arm post is a channel, not sea; the sea test reads the de-slotted surface
        let arm_fs = match &self.v2 {
            Some(v2) if v2.curves.arm.is_some() => {
                let arm_near = near_i.map(|i| v2.curves.is_arm(i)).unwrap_or(false);
                // inside an arm line's carved half-width the channel is cut through low floodplain too (no fixed-sea islands)
                let in_arm = self.p.channels_on && v2.curves.near_lines(x, y, 250.0).iter().any(|q| q.arm && q.r <= q.hw);
                Some((!in_arm && v2.aux.bicubic(&v2.aux.zmod, x, y) <= self.sea + self.p.fixed_above_sea_m) || ((code == 2 || code == 4) && !arm_near))
            }
            _ => None,
        };
        let fixed_sea = arm_fs.unwrap_or_else(|| base <= self.sea + self.p.fixed_above_sea_m || code == 2 || code == 4);
        if let Some(v2) = &self.v2 {
            let hit = if self.p.channels_on { v2.curves.hit(x, y, &self.p) } else { ChannelHit { excess_m: 1e9, ..Default::default() } };
            let cut = if fixed_sea || !v2.v.valley_on { base } else { self.valley(v2, x, y, base) };
            return PostVal { base, cut, cell, fixed_sea, hit };
        }
        let hit = if self.p.channels_on { self.ch.hit(x, y, &self.p) } else { ChannelHit { excess_m: 1e9, ..Default::default() } };
        let mut cut = base;
        if self.p.remap_on && !fixed_sea {
            if let Some(cs) = self.strata.cells.get(&cell) {
                let z1 = self.strata.remap(cs, x, y, base);
                let taper = if self.p.channels_on { smoothstep(hit.excess_m / self.p.remap_channel_taper_m) } else { 1.0 };
                cut = base + taper * (z1 - base);
            }
        }
        PostVal { base, cut, cell, fixed_sea, hit }
    }

    /// v2 V2: valley reprojection (the valley axis follows the curve).
    fn valley(&self, v2: &V2Ctx, x: f64, y: f64, base: f64) -> f64 {
        let a = &v2.aux;
        let zb = a.bicubic(&a.zmod, x, y);
        let Some(gg) = a.bilinear(&a.g, x, y) else { return zb };
        let (mut zv, mut rmin) = (f64::INFINITY, f64::INFINITY);
        let rv = v2.v.valley_rv_m;
        // lane INTS2n4 arm_cl: the tidal arm does not reproject the valley (its valley is T1b's floor; the B6 carve cuts the
        // whole cross-section, so the carve bookkeeping holds all of it)
        for q in v2.curves.near_lines(x, y, rv) { if q.arm { continue; } zv = zv.min(q.bed + gg * (q.r - q.hw).max(0.0)); rmin = rmin.min(q.r); }
        let _ = base;
        if !zv.is_finite() { return zb; }
        let wgt = 1.0 - smoothstep((rmin - 0.5 * rv) / (0.5 * rv));
        zb - wgt * (zb - zv).max(0.0)
    }

    /// v2: (over class, edifice datum S) of a 2 m post.
    pub fn over_at(&self, v2: &V2Ctx, x: f64, y: f64) -> (crate::v2::OverAt, f64) {
        let (e, d) = v2.aux.over_w(x, y);
        let s = if e > 0.0 { v2.aux.bicubic(&v2.aux.s_ed, x, y) } else { 0.0 };
        let vd = match &v2.ed_bs {
            Some(eb) if e > 0.0 => { let ws = self.ed_weights(v2, eb, x, y); let mut b = ws[0]; for &q in &ws[1..] { if q.1 > b.1 + 1e-12 { b = q; } } b.0 }
            _ => 255,
        };
        (crate::v2::OverAt { w_ed: e, w_dep: d, vd }, s)
    }

    /// Lane BS: normalized bilinear weights of the volcanoes of the edifice T1b posts around (x, y).
    pub fn ed_weights(&self, v2: &V2Ctx, eb: &EdBs, x: f64, y: f64) -> Vec<(u8, f64)> {
        let a = &v2.aux;
        let (fx, fy) = (x / T1B_M, y / T1B_M);
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let mut out: Vec<(u8, f64)> = Vec::new();
        let mut tot = 0.0;
        for (dx, dy, wgt) in [(0, 0, (1.0 - tx) * (1.0 - ty)), (1, 0, tx * (1.0 - ty)), (0, 1, (1.0 - tx) * ty), (1, 1, tx * ty)] {
            let gx = (ix + dx - a.gx0).clamp(0, a.nx as i64 - 1) as usize;
            let gy = (iy + dy - a.gy0).clamp(0, a.ny as i64 - 1) as usize;
            let k = gy * a.nx + gx;
            if a.over[k] != crate::v2::OVER_EDIFICE || wgt <= 0.0 { continue; }
            tot += wgt;
            match out.iter_mut().find(|o| o.0 == eb.vid[k]) { Some(o) => o.1 += wgt, None => out.push((eb.vid[k], wgt)) }
        }
        if tot > 0.0 { for o in out.iter_mut() { o.1 /= tot; } } else {
            // no edifice post among the four (the post's over weight came from a farther post): the nearest edifice post's volcano
            let k = ((fy.round() as i64 - a.gy0).clamp(0, a.ny as i64 - 1) as usize) * a.nx + (fx.round() as i64 - a.gx0).clamp(0, a.nx as i64 - 1) as usize;
            out.push((eb.vid[k], 1.0));
        }
        out.sort_by(|p, q| p.0.cmp(&q.0));
        out
    }

    /// The edifice column of the dominant volcano at (x, y) and its pinch column id (the v2/v3 single sequence without BS).
    pub fn ed_dom<'b>(&self, v2: &'b V2Ctx, o: crate::v2::OverAt) -> (u32, &'b CellStrata) {
        match &v2.ed_bs {
            None => (crate::v2::EDIFICE_CELL, &v2.ed),
            Some(eb) => match eb.cols.get(&o.vd) { Some(c) => (ed_col_id(o.vd), c), None => (crate::v2::EDIFICE_CELL, &v2.ed) },
        }
    }

    /// Lane BS2 bed_meta: (facies of the exposed resolved bed, generated unit index) at surface z; (254, 0xFFFF) on deposits.
    pub fn bed_meta_at(&self, v2: &V2Ctx, x: f64, y: f64, cell: u32, z: f64) -> (u8, u16) {
        let (o, s) = self.over_at(v2, x, y);
        let look = |col: &crate::strata::Column, cid: u32, d: f64| -> (u8, u16) {
            let dw = self.wdepth(v2, cid, col, x, y, d);
            let b = col.bed_at_depth(dw);
            (col.bed_fac.get(b).copied().unwrap_or(crate::bedstyle::FAC_NONE), unit_at(col, dw).map(|i| col.unit_top[i].1).unwrap_or(u16::MAX))
        };
        match o.class() {
            crate::v2::OVER_DEPOSIT => (crate::bedstyle::FAC_DEPOSIT, u16::MAX),
            crate::v2::OVER_EDIFICE => { let (cid, ec) = self.ed_dom(v2, o); look(&ec.col, cid, s - z) }
            _ => match self.strata.cells.get(&cell) {
                Some(cs) => look(&cs.col, cell, self.strata.datum(cs, x, y) - z),
                None => (crate::bedstyle::FAC_NONE, u16::MAX),
            },
        }
    }

    /// Lane PYRO: pyroclastic subtype code (pyro::SUB_*) of the exposed bed at surface z off the lobe stack (the per-volcano 1D
    /// sequence on edifice posts: one recorded event class per package; strata posts: 254 unresolved; others 255).
    pub fn pyro_meta_at(&self, v2: &V2Ctx, x: f64, y: f64, cell: u32, z: f64) -> u8 {
        use crate::pyro::{SUB_NONE, SUB_UNRESOLVED};
        let is_pyro = |f: u8| f == crate::bedstyle::FAC_PYRO || f == crate::bedstyle::FAC_PYRO_HARD;
        let (o, s) = self.over_at(v2, x, y);
        match o.class() {
            crate::v2::OVER_DEPOSIT => SUB_NONE,
            crate::v2::OVER_EDIFICE => {
                let (cid, ec) = self.ed_dom(v2, o);
                let col = &ec.col;
                let dw = self.wdepth(v2, cid, col, x, y, s - z);
                let b = col.bed_at_depth(dw);
                if !col.bed_fac.get(b).copied().map(is_pyro).unwrap_or(false) { return SUB_NONE; }
                let (Some(pw), Some(i)) = (v2.pyro.as_ref(), unit_at(col, dw)) else { return SUB_UNRESOLVED };
                if v2.ed_bs.as_ref().map(|eb| eb.cols.contains_key(&o.vd)) != Some(true) { return SUB_UNRESOLVED; }
                let (top, unit) = col.unit_top[i];
                let base = col.unit_top.get(i + 1).map(|q| q.0).unwrap_or(*col.base_depth.last().unwrap());
                let tu = (base - top).max(1e-6);
                pw.ed_code(&self.layer_scope, o.vd, unit, (dw - top) / tu, tu, x, y)
            }
            _ => match self.strata.cells.get(&cell) {
                Some(cs) => {
                    let col = &cs.col;
                    let b = col.bed_at_depth(self.wdepth(v2, cell, col, x, y, self.strata.datum(cs, x, y) - z));
                    if col.bed_fac.get(b).copied().map(is_pyro).unwrap_or(false) { SUB_UNRESOLVED } else { SUB_NONE }
                }
                None => SUB_NONE,
            },
        }
    }

    /// Lane BS unit_depth: (depth below the exposed unit's top contact, height above its base contact) at surface z; NaN where
    /// not defined (layer 0's top is the present land surface; edifice top; deposits; no column).
    pub fn unit_depth(&self, v2: &V2Ctx, x: f64, y: f64, cell: u32, z: f64) -> (f32, f32) {
        let (o, _) = self.over_at(v2, x, y);
        match o.class() {
            crate::v2::OVER_DEPOSIT => (f32::NAN, f32::NAN),
            // lane BS2 ed_contacts: the contacts of the generated flow unit / package containing the point (bed coordinates)
            crate::v2::OVER_EDIFICE if self.p.beds.as_ref().and_then(|b| b.ed_contacts) == Some(true) && v2.ed_bs.is_some() => {
                let (o2, s) = self.over_at(v2, x, y);
                let (cid, ec) = self.ed_dom(v2, o2);
                let dw = self.wdepth(v2, cid, &ec.col, x, y, s - z);
                match unit_at(&ec.col, dw) {
                    Some(i) if ec.col.unit_top[i].1 != u16::MAX => {
                        let base = ec.col.unit_top.get(i + 1).map(|q| q.0).unwrap_or(*ec.col.base_depth.last().unwrap());
                        ((dw - ec.col.unit_top[i].0) as f32, (base - dw) as f32)
                    }
                    _ => (f32::NAN, f32::NAN),
                }
            }
            crate::v2::OVER_EDIFICE => match &self.w.base_m {
                Some(b) => (f32::NAN, (z - self.w.bicubic_of(b, x, y)) as f32),
                None => (f32::NAN, f32::NAN),
            },
            _ => match self.strata.cells.get(&cell) {
                Some(cs) => {
                    let d = self.strata.datum(cs, x, y) - z;
                    let k = cs.col.layer_at_depth(d);
                    let top = if k == 0 { f64::NAN } else { cs.col.base_depth[k - 1] };
                    ((d - top) as f32, (cs.col.base_depth[k] - d) as f32)
                }
                None => (f32::NAN, f32::NAN),
            },
        }
    }

    /// v3 W2: pinch amplitudes of column `col_id` at (x, y) (None without pinch).
    #[inline]
    pub fn pinch_amps(&self, v2: &V2Ctx, col_id: u32, x: f64, y: f64) -> Option<[f64; 2]> {
        v2.pinch.as_ref().map(|p| p.amps(&self.layer_scope, col_id, x, y))
    }
    /// v3 W2: depth below a column datum -> bed-coordinate depth (identity without pinch).
    #[inline]
    pub fn wdepth(&self, v2: &V2Ctx, col_id: u32, col: &crate::strata::Column, x: f64, y: f64, d: f64) -> f64 {
        match (self.pinch_amps(v2, col_id, x, y), v2.pinch.as_ref()) {
            (Some(a), Some(p)) => crate::v3::warp(col, &a, &p.periods, d, crate::v3::UP_H),
            _ => d,
        }
    }

    /// v2: (Sc, K) of the bed exposed at z for a post with over weights `o` and edifice datum `s`.
    pub fn bed_props_v2(&self, v2: &V2Ctx, x: f64, y: f64, cell: u32, z: f64, o: crate::v2::OverAt, s: f64) -> (f64, f64) {
        match o.class() {
            crate::v2::OVER_DEPOSIT => (v2.v.deposit_sc, v2.v.deposit_k),
            crate::v2::OVER_EDIFICE => { let (cid, ec) = self.ed_dom(v2, o); let b = ec.col.bed_at_depth(self.wdepth(v2, cid, &ec.col, x, y, s - z)); (ec.col.bed_sc[b], ec.col.bed_k[b]) }
            _ => match self.strata.cells.get(&cell) {
                Some(cs) if v2.pinch.is_some() => { let b = cs.col.bed_at_depth(self.wdepth(v2, cell, &cs.col, x, y, self.strata.datum(cs, x, y) - z)); (cs.col.bed_sc[b], cs.col.bed_k[b]) }
                Some(cs) => self.strata.bed_props(cs, x, y, z),
                None => (sc_of(0), k_of(0)),
            },
        }
    }

    /// v3 W4: (Sc, K) of the remap window at z (window harmonic-mean Sc, window mean K), unwarped stratigraphic height.
    pub fn win_props_v3(&self, v2: &V2Ctx, x: f64, y: f64, cell: u32, z: f64, o: crate::v2::OverAt, s: f64) -> (f64, f64) {
        match o.class() {
            crate::v2::OVER_DEPOSIT => (v2.v.deposit_sc, v2.v.deposit_k),
            crate::v2::OVER_EDIFICE => { let ec = self.ed_dom(v2, o).1; (crate::strata::sc_win_at(ec, z - s), crate::strata::k_win_at(ec, z - s)) }
            _ => match self.strata.cells.get(&cell) {
                Some(cs) => { let u = z - self.strata.datum(cs, x, y); (crate::strata::sc_win_at(cs, u), crate::strata::k_win_at(cs, u)) }
                None => (sc_of(0), k_of(0)),
            },
        }
    }

    /// Lane DE DE5: the lobe column of an edifice post (None: not edifice, lobes off, or no lobes for its volcano). The column is
    /// stacked down from the construction surface D = max(S, T1b bicubic); below the stack the per-volcano 1D sequence continues
    /// (its depth below S matched at the stack base).
    pub fn lobe_at<'b>(&self, v2: &'b V2Ctx, x: f64, y: f64, o: crate::v2::OverAt, s: f64) -> Option<crate::lobes::LobeAt<'b>> {
        let lw = v2.lobes.as_ref()?;
        if o.class() != crate::v2::OVER_EDIFICE { return None; }
        let vl = lw.volc.get(&o.vd)?;
        let tb = match &v2.aux.zpre { Some(z) => v2.aux.bicubic(z, x, y), None => self.w.bicubic(x, y) };
        // lane LOBE de_lobe_relief: the construction surface carries the flows' constructional relief (RELIEF_NOTE)
        let tb = if lw.relief.is_some() { tb + self.relief_post(lw, &self.post((x / POST_M).round() as i64, (y / POST_M).round() as i64), x, y) } else { tb };
        let d = s.max(tb);
        let col = vl.column(x, y, d);
        let ed = self.ed_dom(v2, o).1;
        let ed_off = col.stack_m - (d - s);
        Some(crate::lobes::LobeAt { col, ed, ed_off, vl: Some(vl) })
    }

    /// Lane LOBE de_lobe_relief: the relief of a post: 0 on sea / corridor posts, tapered to 0 towards a corridor
    /// (smoothstep(excess / corridor_taper_m)); a pure function of the post.
    pub fn relief_post(&self, lw: &crate::lobes::LobeWorld, pv: &PostVal, x: f64, y: f64) -> f64 {
        let Some(rw) = lw.relief.as_ref() else { return 0.0 };
        if pv.fixed_relax() { return 0.0; }
        let t = crate::smoothstep(pv.hit.excess_m / rw.p.corridor_taper_m.max(1e-9));
        if t <= 0.0 { 0.0 } else { t * rw.at(x, y) }
    }

    /// Lane DE: (Sc, K) of the bed exposed at z (edifice: the lobe column when de_lobes=1, else bed_props_v2).
    pub fn de_bed(&self, v2: &V2Ctx, x: f64, y: f64, cell: u32, z: f64, o: crate::v2::OverAt, s: f64) -> (f64, f64) {
        self.bed_props_v2(v2, x, y, cell, z, o, s)
    }

    /// Lane DE: window (reference) critical slope at z (edifice: the lobe column window when de_lobes=1).
    pub fn de_win(&self, v2: &V2Ctx, x: f64, y: f64, cell: u32, z: f64, o: crate::v2::OverAt, s: f64) -> f64 {
        self.win_props_v3(v2, x, y, cell, z, o, s).0
    }

    /// v2: (layer index, litho code, bed index) exposed at z (layer 254 = deposit, 253 = edifice).
    pub fn exposed_v2(&self, v2: &V2Ctx, x: f64, y: f64, cell: u32, z: f64) -> (usize, u8, u16) {
        let (o, s) = self.over_at(v2, x, y);
        match o.class() {
            crate::v2::OVER_DEPOSIT => (254, crate::v2::LITHO_DEPOSIT, 0),
            crate::v2::OVER_EDIFICE => { let (cid, ec) = self.ed_dom(v2, o); (253, crate::v2::LITHO_VOLCANIC, ec.col.bed_at_depth(self.wdepth(v2, cid, &ec.col, x, y, s - z)).min(65535) as u16) }
            _ => match self.strata.cells.get(&cell) {
                Some(cs) => { let (l, c) = self.strata.exposed(cs, x, y, z); (l, c, cs.col.bed_at_depth(self.wdepth(v2, cell, &cs.col, x, y, self.strata.datum(cs, x, y) - z)).min(65535) as u16) }
                None => (255, 255, 0),
            },
        }
    }

    /// v2 V5-V7: the stratigraphic remap of a surface value z at a post (applied after relaxation and the seam blend).
    /// Returns the bounded, slope-clamped, channel-tapered displacement.
    pub fn remap_v2(&self, v2: &V2Ctx, pv: &PostVal, x: f64, y: f64, z: f64) -> f64 {
        if !self.p.remap_on || pv.fixed_sea { return 0.0; }
        // lane DE: the layered hillslope equilibrium (inside the relaxation) replaces the bounded remap
        if self.p.de.is_some() { return 0.0; }
        let v = &v2.v;
        let (o, s) = self.over_at(v2, x, y);
        let wd = crate::v2::smooth_w(o.w_dep);
        let we = crate::v2::smooth_w(o.w_ed);
        if wd >= 1.0 { return 0.0; }
        let h = 1.0;
        let grad = |x: f64, y: f64| -> f64 {
            let gxs = (self.w.bicubic(x + h, y) - self.w.bicubic(x - h, y)) / (2.0 * h);
            let gys = (self.w.bicubic(x, y + h) - self.w.bicubic(x, y - h)) / (2.0 * h);
            gxs.hypot(gys)
        };
        let g0 = grad(x, y);
        // T2D X1 (revision 1): the slope clamp reads the steepest base slope within clamp_h_m (centre + 8 directions), so the
        // remap fades out before a cliff lip instead of switching off within 2-4 m of it
        let g0c = match v.v3.as_ref().and_then(|w| w.clamp_h_m) {
            None => g0,
            Some(r) => {
                let mut m = g0;
                for k in 0..8 { let a = k as f64 * std::f64::consts::FRAC_PI_4; m = m.max(grad(x + r * a.cos(), y + r * a.sin())); }
                m
            }
        };
        let rh = v.bed_hard_mult;
        // v3 W3 soil mantle weight and W2 warped remap
        let expose = |scw: f64| -> f64 { match &v.v3 { Some(w) if w.expose_on => smoothstep((g0 / scw - w.expose_lo) / w.expose_width), _ => 1.0 } };
        let remap_col = |cs: &CellStrata, col_id: u32, dat: f64| -> (f64, f64) {
            match (self.pinch_amps(v2, col_id, x, y), v2.pinch.as_ref()) {
                (Some(a), Some(p)) => (crate::strata::remap_at_w(cs, dat, z, &a, &p.periods), -crate::v3::warp(&cs.col, &a, &p.periods, dat - z, crate::v3::UP_H)),
                _ => (crate::strata::remap_at(cs, dat, z), z - dat),
            }
        };
        let mut ds = 0.0;
        if we < 1.0 {
            if let Some(cs) = self.strata.cells.get(&pv.cell) {
                let (z1, uw) = remap_col(cs, pv.cell, self.strata.datum(cs, x, y));
                let scw = crate::strata::sc_win_at(cs, uw);
                let a = if v.slope_clamp_on { crate::v2::slope_alpha(g0c, scw, rh) } else { 1.0 };
                ds = a * expose(scw) * (z1 - z);
            }
        }
        let mut de = 0.0;
        if we > 0.0 {
            let one = |ec: &CellStrata, cid: u32| -> f64 {
                let (z1, uw) = remap_col(ec, cid, s);
                let scw = crate::strata::sc_win_at(ec, uw);
                let a = if v.slope_clamp_on { crate::v2::slope_alpha(g0c, scw, rh) } else { 1.0 };
                a * expose(scw) * (z1 - z)
            };
            match &v2.ed_bs {
                None => de = one(&v2.ed, crate::v2::EDIFICE_CELL),
                // lane BS: per-volcano sequences blended by the bilinear weights of their edifice posts
                Some(eb) => for (vv, wv) in self.ed_weights(v2, eb, x, y) { if let Some(ec) = eb.cols.get(&vv) { de += wv * one(ec, ed_col_id(vv)); } },
            }
        }
        let d = ((1.0 - we) * ds + we * de) * (1.0 - wd);
        let r = v2.aux.bilinear(&v2.aux.relief, x, y).unwrap_or(0.0);
        let b = (v.amp_frac_of_relief * r).min(v.amp_cap_m);
        let d = crate::v2::bound(d, b);
        let taper = if self.p.channels_on { smoothstep(pv.hit.excess_m / self.p.remap_channel_taper_m) } else { 1.0 };
        taper * d
    }

    /// Exposed (layer index, lithology code) at surface z of post (px, py).
    pub fn exposed(&self, px: i64, py: i64, cell: u32, z: f64) -> (usize, u8) {
        match self.strata.cells.get(&cell) {
            Some(cs) => self.strata.exposed(cs, px as f64 * POST_M, py as f64 * POST_M, z),
            None => (255, 255),
        }
    }

    /// Tile T1b local relief (max - min of the bicubic base over the core): the closure's relief basis.
    pub fn tile_relief(&self, ix: i64, iy: i64) -> f64 {
        if let Some(v) = self.relief_cache.lock().unwrap().get(&(ix, iy)) { return *v; }
        let (ox, oy) = tile_origin_post(ix, iy);
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for j in 0..TILE_POSTS as i64 { for i in 0..TILE_POSTS as i64 {
            let b = self.w.bicubic((ox + i) as f64 * POST_M, (oy + j) as f64 * POST_M);
            lo = lo.min(b); hi = hi.max(b);
        } }
        self.relief_cache.lock().unwrap().insert((ix, iy), hi - lo);
        hi - lo
    }

    /// v2 V7: roughness rms from the pointwise T1b local relief (v1 used the owning tile's relief).
    pub fn closure_rms_v2(&self, v2: &V2Ctx, px: i64, py: i64) -> f64 {
        if self.p.closure_rms_cap_m <= 0.0 { return 0.0; }
        let r = v2.aux.bilinear(&v2.aux.relief, px as f64 * POST_M, py as f64 * POST_M).unwrap_or(0.0);
        self.p.closure_rms_cap_m.min(self.p.closure_rms_fraction_of_relief * r)
    }

    pub fn closure_rms(&self, ix: i64, iy: i64) -> f64 {
        if self.p.closure_rms_cap_m <= 0.0 { return 0.0; }
        self.p.closure_rms_cap_m.min(self.p.closure_rms_fraction_of_relief * self.tile_relief(ix, iy))
    }

    /// One keyed standard-normal draw for global post (px, py), scoped by its OWNING lod-0 tile (half-open).
    pub fn draw(&self, px: i64, py: i64) -> f64 {
        let (ix, iy) = (px.div_euclid(TILE_CELLS), py.div_euclid(TILE_CELLS));
        let (lx, ly) = (px.rem_euclid(TILE_CELLS), py.rem_euclid(TILE_CELLS));
        let s = self.layer_scope.tile(0, ix, iy).and_then(|s| s.address(&format!("p/{lx}/{ly}"))).expect("scope");
        s.draw(PURPOSE_ROUGHNESS, 0, 0).expect("draw").normal_f64()
    }

    /// B4 closure field (m) on an n x n block of global posts starting at (px0, py0).
    pub fn closure_field(&self, px0: i64, py0: i64, n: usize) -> Vec<f64> {
        let m = n + 2;
        let mut d = vec![0.0f64; m * m];
        let mut rms = vec![0.0f64; m * m];
        let mut any = false;
        for j in 0..m { for i in 0..m {
            let (px, py) = (px0 - 1 + i as i64, py0 - 1 + j as i64);
            let r = match &self.v2 { Some(v2) => self.closure_rms_v2(v2, px, py), None => self.closure_rms(px.div_euclid(TILE_CELLS), py.div_euclid(TILE_CELLS)) };
            rms[j * m + i] = r;
            if r > 0.0 { any = true; d[j * m + i] = self.draw(px, py); }
        } }
        let mut out = vec![0.0f64; n * n];
        if !any { return out; }
        let k = [1.0, 2.0, 1.0];
        let corr = self.v2.as_ref().and_then(|v| v.v.v3.as_ref()).filter(|w| w.rough_corr_on);
        for j in 0..n { for i in 0..n {
            let mut s = 0.0;
            for (b, kb) in k.iter().enumerate() { for (a, ka) in k.iter().enumerate() {
                s += ka * kb * d[(j + b) * m + (i + a)];
            } }
            // rms of the post's own owning tile
            let r = rms[(j + 1) * m + (i + 1)];
            out[j * n + i] = match corr {
                None => r * (s / 16.0) / 0.375,
                // v3 W5: correlated keyed lattice part plus the v2 white part
                Some(w) => if r > 0.0 { r * (w.rough_corr_weight * self.rough_lattice(px0 + i as i64, py0 + j as i64, w.rough_lattice_posts) + w.rough_white_weight * (s / 16.0) / 0.375) } else { 0.0 },
            };
        } }
        out
    }

    /// v3 W5: keyed lattice roughness (nodes every `step` posts, keyed by the node's owning lod-0 tile), smoothstep-bilinear.
    fn rough_lattice(&self, px: i64, py: i64, step: i64) -> f64 {
        let (nx, ny) = (px.div_euclid(step), py.div_euclid(step));
        let (tx, ty) = (smoothstep(px.rem_euclid(step) as f64 / step as f64), smoothstep(py.rem_euclid(step) as f64 / step as f64));
        let g = |a: i64, b: i64| -> f64 {
            let key = (nx + a, ny + b);
            if let Some(v) = self.rough_cache.lock().unwrap().get(&key) { return *v; }
            let (qx, qy) = (key.0 * step, key.1 * step);
            let (ix, iy) = (qx.div_euclid(TILE_CELLS), qy.div_euclid(TILE_CELLS));
            let sc = self.layer_scope.tile(0, ix, iy).and_then(|s| s.address(&format!("rc/{}/{}", qx.rem_euclid(TILE_CELLS), qy.rem_euclid(TILE_CELLS)))).expect("scope");
            let v = sc.draw(crate::v3::PURPOSE_ROUGH_V3, 0, 0).expect("draw").normal_f64();
            self.rough_cache.lock().unwrap().insert(key, v);
            v
        };
        (g(0, 0) * (1.0 - tx) + g(1, 0) * tx) * (1.0 - ty) + (g(0, 1) * (1.0 - tx) + g(1, 1) * tx) * ty
    }

    /// Lane RIB2 (de_ringflow): T1b upslope drainage area entering the B5 domain (n x n posts from (px0, py0)), placed on the
    /// second ring (the outer ring is fixed and passes nothing). Mode 1: the inward flux across each ring-1 post's edge length,
    /// specific area (T1b area / T1b spacing) x 2 m x max(0, inward normal . T1b flow direction), flow direction = steepest
    /// descent of the T1b bicubic (16 m central differences). Mode 2: each T1b post outside the inner box whose D8 segment to
    /// its steepest-descent T1b neighbour crosses into the box delivers its T1b drainage area at the crossing point (split
    /// linearly between the two nearest ring-1 posts). A pure function of the domain (so of the tile key).
    pub fn ring_inflow(&self, px0: i64, py0: i64, n: usize, mode: u32) -> Vec<f64> {
        let mut extra = vec![0.0f64; n * n];
        let (x0, y0) = (px0 as f64 * POST_M, py0 as f64 * POST_M);
        let (bx0, by0, bx1, by1) = (x0 + POST_M, y0 + POST_M, x0 + (n - 2) as f64 * POST_M, y0 + (n - 2) as f64 * POST_M);
        let w = self.w;
        if mode == 1 {
            let h = 16.0;
            // lane INTS2n4 ringflow_ch (chn.rs RINGFLOW_CH_RULE): no channel / marine T1b area on the ring
            let ring_ch = self.p.chn.as_ref().map(|c| c.ringflow_ch.is_some()).unwrap_or(false);
            for j in 1..n - 1 { for i in 1..n - 1 {
                if !(i == 1 || j == 1 || i == n - 2 || j == n - 2) { continue; }
                let (x, y) = (x0 + i as f64 * POST_M, y0 + j as f64 * POST_M);
                let (gx, gy) = ((w.bicubic(x + h, y) - w.bicubic(x - h, y)) / (2.0 * h), (w.bicubic(x, y + h) - w.bicubic(x, y - h)) / (2.0 * h));
                let g = gx.hypot(gy);
                if g <= 1e-9 { continue; }
                let (dx, dy) = (-gx / g, -gy / g);
                let a = if ring_ch { w.area_bilinear_free(x, y) } else { w.area_bilinear(x, y) } / T1B_M * POST_M;
                let mut f = 0.0;
                if i == 1 { f += dx.max(0.0); }
                if i == n - 2 { f += (-dx).max(0.0); }
                if j == 1 { f += dy.max(0.0); }
                if j == n - 2 { f += (-dy).max(0.0); }
                extra[j * n + i] = a * f;
            } }
            return extra;
        }
        // mode 2: T1b D8 flow lines crossing into the inner box
        let inside = |x: f64, y: f64| x >= bx0 && x <= bx1 && y >= by0 && y <= by1;
        let (gx0, gy0) = ((bx0 / T1B_M).floor() as i64 - 2, (by0 / T1B_M).floor() as i64 - 2);
        let (gx1, gy1) = ((bx1 / T1B_M).ceil() as i64 + 2, (by1 / T1B_M).ceil() as i64 + 2);
        for gy in gy0..=gy1 { for gx in gx0..=gx1 {
            let Some(p) = w.idx(gx, gy) else { continue };
            let (px, py) = (gx as f64 * T1B_M, gy as f64 * T1B_M);
            if inside(px, py) { continue; }
            let zp = w.height_m[p];
            let mut best = (0.0f64, 0i64, 0i64);
            for (dx, dy) in [(1i64, 0i64), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)] {
                if let Some(q) = w.idx(gx + dx, gy + dy) {
                    let s = (zp - w.height_m[q]) / if dx != 0 && dy != 0 { std::f64::consts::SQRT_2 } else { 1.0 };
                    if s > best.0 { best = (s, dx, dy); }
                }
            }
            if best.0 <= 0.0 { continue; }
            let (qx, qy) = (px + best.1 as f64 * T1B_M, py + best.2 as f64 * T1B_M);
            // Liang-Barsky entry parameter of the segment p -> q into the box
            let (ddx, ddy) = (qx - px, qy - py);
            let (mut t0, mut t1) = (0.0f64, 1.0f64);
            let mut ok = true;
            for (pp, qq) in [(-ddx, px - bx0), (ddx, bx1 - px), (-ddy, py - by0), (ddy, by1 - py)] {
                if pp == 0.0 { if qq < 0.0 { ok = false; } continue; }
                let r = qq / pp;
                if pp < 0.0 { if r > t1 { ok = false; } else if r > t0 { t0 = r; } } else if r < t0 { ok = false; } else if r < t1 { t1 = r; }
            }
            if !ok || t0 > t1 { continue; }
            let (ex, ey) = (px + t0 * ddx, py + t0 * ddy);
            let a = w.drainage_area_m2[p] as f64;
            // ring-1 post coordinates of the entry point, split along the side
            let (fi, fj) = (((ex - x0) / POST_M).clamp(1.0, (n - 2) as f64), ((ey - y0) / POST_M).clamp(1.0, (n - 2) as f64));
            let on_x_side = (ex - bx0).abs() < 1e-6 || (ex - bx1).abs() < 1e-6;
            if on_x_side {
                let i = fi.round() as usize; let j0 = fj.floor() as usize; let t = fj - j0 as f64; let j1 = (j0 + 1).min(n - 2);
                extra[j0 * n + i] += a * (1.0 - t); extra[j1 * n + i] += a * t;
            } else {
                let j = fj.round() as usize; let i0 = fi.floor() as usize; let t = fi - i0 as f64; let i1 = (i0 + 1).min(n - 2);
                extra[j * n + i0] += a * (1.0 - t); extra[j * n + i1] += a * t;
            }
        } }
        extra
    }

    /// Lane RILL (de_b5_chain): the coarse window-wide B5 pass (chain.rs), built once per context.
    pub fn chain_coarse(&self) -> &crate::chain::Coarse { self.chain.get_or_init(|| self.build_coarse()) }

    fn build_coarse(&self) -> crate::chain::Coarse {
        let t = Instant::now();
        let dp = self.p.de.as_ref().expect("de_b5_chain needs de=1");
        let cp = dp.chain.as_ref().expect("de_b5_chain");
        let v2 = self.v2.as_ref().expect("de_b5_chain needs version 3");
        let st = cp.cell_posts as i64;
        let r = (T1B_M / POST_M).round() as i64;
        let (fx0, fy0) = (self.w.gx0 * r, self.w.gy0 * r);
        let n = ((self.w.n as i64 - 1) * r / st + 1) as usize;
        let nn = n * n;
        let win = self.p.v3p().map(|w| w.relax_win_on).unwrap_or(false);
        let (sc3, k3) = (sc_of(crate::v2::LITHO_VOLCANIC), k_of(crate::v2::LITHO_VOLCANIC));
        let rop = dp.runoff.as_ref();
        let use_bk = dp.sp_area_m2.is_some();
        // lane CLIMATE-RUNOFF climate_runoff=1: the coarse pass's own areas also x the T1b runoff weight
        let clim = self.p.climate_runoff.is_some();
        let use_ro = (rop.is_some() && v2.lobes.is_some()) || clim;
        let (mut z, mut fixed) = (vec![0.0f64; nn], vec![false; nn]);
        let mut kind = vec![0u8; nn];
        let (mut sc, mut kk, mut bk, mut rf, mut rt) = (vec![0.0f64; nn], vec![0.0f64; nn], vec![0.0f64; nn], vec![1.0f64; nn], vec![1.0f64; nn]);
        for cj in 0..n { for ci in 0..n {
            let k = cj * n + ci;
            let (px, py) = (fx0 + ci as i64 * st, fy0 + cj as i64 * st);
            let (x, y) = (px as f64 * POST_M, py as f64 * POST_M);
            let pv = self.post(px, py);
            fixed[k] = pv.fixed_relax() || ci == 0 || cj == 0 || ci == n - 1 || cj == n - 1;
            kind[k] = if pv.fixed_sea { 1 } else if pv.fixed_relax() { 2 } else { 0 };
            // the tile domains' initial surface at this post (T1a cut + B4 closure + lobe relief)
            z[k] = pv.cut + if fixed[k] { 0.0 } else { self.closure_field(px, py, 1)[0] };
            if let Some(lw) = v2.lobes.as_ref().filter(|l| l.relief.is_some()) { z[k] += self.relief_post(lw, &pv, x, y); }
            let (o, s) = self.over_at(v2, x, y);
            let la = if v2.lobes.is_some() { self.lobe_at(v2, x, y, o, s) } else { None };
            let (a, b) = match &la {
                Some(la) => { let u = la.u(z[k]); if win { let m = la.inv_mean(u - 30.0, u + 30.0); (1.0 / m, k3 * sc3 * m) } else { let q = la.sc(u); (q, k3 * sc3 / q) } }
                None => if win { self.win_props_v3(v2, x, y, pv.cell, z[k], o, s) } else { self.bed_props_v2(v2, x, y, pv.cell, z[k], o, s) },
            };
            sc[k] = a; kk[k] = b;
            if use_bk { bk[k] = match &la { Some(la) => k3 * sc3 / la.sc(la.u(z[k])), None => self.de_bed(v2, x, y, pv.cell, z[k], o, s).1 }; }
            if use_ro { if let (Some(rp), Some(la)) = (rop, &la) { let (f, tt) = rp.factor(la.fac_any(la.u(z[k]))); rf[k] = f; rt[k] = if cp.no_loss == Some(true) { 1.0 } else { tt.powf(st as f64) }; } }
            if clim { rf[k] *= self.w.runoff_weight_at(x, y).unwrap_or(1.0); }
        } }
        let z_start = z.clone();
        let t_init = t.elapsed().as_secs_f64();
        let mut pc = self.p.clone();
        pc.steps = (self.p.steps / cp.dt).max(1);
        pc.k_sp_per_step *= cp.dt as f64;
        pc.d_creep_m2_per_step *= cp.dt as f64;
        let dx = st as f64 * POST_M;
        // lane RILL2 de_b5_chain_head=1: the channel-head threshold at the coarse contour width (chain.rs `head_scale`)
        if cp.head_scale == Some(true) { pc.channel_head_area_m2 *= dx / POST_M; }
        let litho = |k: usize, _z: f64| (sc[k], kk[k]);
        let bed_k = |k: usize, _z: f64| bk[k];
        let runoff = |k: usize, _z: f64| (rf[k], rt[k]);
        let zero = vec![0.0f64; nn];
        let stt = crate::relax::relax_dx(n, dx, &mut z, &fixed, &zero, &litho, if use_bk { Some(&bed_k) } else { None }, if use_ro { Some(&runoff) } else { None }, &pc);
        let area = crate::relax::mfd_area(n, dx, &z, &fixed, if use_ro { Some((&rf, &rt)) } else { None }, &pc);
        let pass: Vec<f64> = (0..nn).map(|k| if use_ro { area[k] * rt[k] } else { area[k] }).collect();
        let mut dz: Vec<f64> = (0..nn).map(|k| z[k] - z_start[k]).collect();
        // mode 3: the W6 local mass balance applied to the coarse change on the coarse lattice (masked Gaussian low-pass removed)
        if cp.mode == 3 {
            if let Some(w) = self.p.v3p().filter(|w| w.massbal_on) {
                let free: Vec<bool> = fixed.iter().map(|f| !f).collect();
                let lp = crate::v3::masked_lowpass(n, &dz, &free, w.massbal_sigma_m / dx);
                for k in 0..nn { if free[k] { dz[k] -= lp[k]; } }
            }
        }
        let secs = t.elapsed().as_secs_f64();
        eprintln!("chain coarse: n {n} x {n} at {dx} m, {} steps (dt {}), init {t_init:.1} s, total {secs:.1} s, max change {:.2} m, mean change {:.3} m, max area {:.3e} m2",
                  pc.steps, cp.dt, stt.max_change_m, stt.mean_change_m, stt.max_area_m2);
        // diagnostic de_b5_chain_dump=<dir>: the coarse rasters (row-major, row = north index) for offline analysis
        if let Some(dir) = cp.dump.as_ref() {
            let f32s = |name: &str, v: &[f64]| { let b: Vec<u8> = v.iter().flat_map(|x| (*x as f32).to_le_bytes()).collect(); std::fs::write(format!("{dir}/{name}.f32"), b).expect("dump"); };
            std::fs::create_dir_all(dir).expect("dump dir");
            f32s("z0", &z_start); f32s("z", &z); f32s("dz", &dz); f32s("area", &area); f32s("sc", &sc); f32s("k", &kk); f32s("bk", &bk); f32s("rf", &rf); f32s("rt", &rt);
            std::fs::write(format!("{dir}/fixed.u8"), fixed.iter().map(|&f| f as u8).collect::<Vec<u8>>()).expect("dump");
            std::fs::write(format!("{dir}/kind.u8"), &kind).expect("dump");
            std::fs::write(format!("{dir}/meta.json"), format!("{{\"n\": {n}, \"fx0\": {fx0}, \"fy0\": {fy0}, \"st\": {st}, \"dx_m\": {dx}, \"sea_m\": {}}}", self.sea)).expect("dump");
        }
        crate::chain::Coarse { fx0, fy0, st, n, z, dz, fixed, pass, seconds: secs, interp: cp.interp.unwrap_or(0) }
    }

    /// Lane PF-SERVE-2 (memo eviction for a resident server; the batch never calls it): drop every memo entry that belongs to
    /// tile domain (ix, iy) TOGETHER: the relaxation, the lane DE change, the lane TREAD pre-stage, the tile relief and the
    /// v3 W5 roughness lattice nodes the tile owns. Returns the bytes of the f64/i32/u8 payloads freed. Every entry is a pure
    /// function of its key, so a later `relaxed` / `assemble` recomputes the same values. Dropping the relaxation alone would
    /// be WRONG: `assemble` reads neighbour DE changes straight from `de_cache` with zeros on a miss. The relaxation goes first
    /// (`relaxed` inserts it last), so no reader sees a relaxation without its DE / TREAD entries. Caller's contract: no
    /// `relaxed(ix, iy)` may be running and no `assemble` that reads (ix, iy) (the 3 x 3 neighbourhood) may be running.
    pub fn forget_tile(&self, ix: i64, iy: i64) -> usize {
        let mut b = 0usize;
        if let Some((z, _)) = self.relax_cache.lock().unwrap().remove(&(ix, iy)) { b += z.len() * 8; }
        if let Some((e, t, _)) = self.de_cache.lock().unwrap().remove(&(ix, iy)) { b += (e.len() + t.len()) * 8; }
        if let Some((s, r, _)) = self.sharp_cache.lock().unwrap().remove(&(ix, iy)) { b += s.len() * 4 + r.len(); }
        if self.relief_cache.lock().unwrap().remove(&(ix, iy)).is_some() { b += 8; }
        if let Some(step) = self.v2.as_ref().and_then(|v| v.v.v3.as_ref()).filter(|w| w.rough_corr_on).map(|w| w.rough_lattice_posts) {
            let mut rc = self.rough_cache.lock().unwrap();
            let n0 = rc.len();
            rc.retain(|&(nx, ny), _| ((nx * step).div_euclid(TILE_CELLS), (ny * step).div_euclid(TILE_CELLS)) != (ix, iy));
            b += (n0 - rc.len()) * 8;
        }
        b
    }

    /// Lane PF-SERVE-2: bytes of the f64/i32/u8 payloads `forget_tile(ix, iy)` would free now (0 when nothing is memoised).
    pub fn memo_tile_bytes(&self, ix: i64, iy: i64) -> usize {
        let mut b = 0usize;
        if let Some((z, _)) = self.relax_cache.lock().unwrap().get(&(ix, iy)) { b += z.len() * 8; }
        if let Some((e, t, _)) = self.de_cache.lock().unwrap().get(&(ix, iy)) { b += (e.len() + t.len()) * 8; }
        if let Some((s, r, _)) = self.sharp_cache.lock().unwrap().get(&(ix, iy)) { b += s.len() * 4 + r.len(); }
        b
    }

    /// Lane PF-SERVE-2: entry counts of the per-tile memos (relax, de, sharp, relief, rough nodes), for status and tests.
    pub fn memo_counts(&self) -> [usize; 5] {
        [self.relax_cache.lock().unwrap().len(), self.de_cache.lock().unwrap().len(), self.sharp_cache.lock().unwrap().len(),
            self.relief_cache.lock().unwrap().len(), self.rough_cache.lock().unwrap().len()]
    }

    /// B4 + B5: the relaxed field of tile (ix, iy) over its domain (tile + halo), pure function of the key.
    pub fn relaxed(&self, ix: i64, iy: i64) -> (Vec<f64>, f64) {
        // Evolved parents have already run their process. No construction memo
        // or second relaxation is needed by their pointwise refinement.
        if self.w.evolved.is_some() { return (Vec::new(),0.); }
        if let Some(v) = self.relax_cache.lock().unwrap().get(&(ix, iy)) { return v.clone(); }
        let t = Instant::now();
        let h = self.p.halo_posts as i64;
        let n = TILE_POSTS + 2 * self.p.halo_posts;
        let (ox, oy) = tile_origin_post(ix, iy);
        let (px0, py0) = (ox - h, oy - h);
        let mut z = vec![0.0f64; n * n];
        let mut fixed = vec![false; n * n];
        let mut cells = vec![0u32; n * n];
        let mut extra = vec![0.0f64; n * n];
        // lane TREAD: pre-stage staircase (None unless ct_sharp=1 and ct_sharp_stage=pre)
        let pre_sf = self.sharp.as_ref().filter(|_| self.p.ct_sharp.as_ref().map(|s| s.pre()).unwrap_or(false));
        let nn_pre = if pre_sf.is_some() { n * n } else { 0 };
        let (mut pre_shift, mut pre_rg, mut pre_ta, mut riser) = (vec![0.0f64; nn_pre], vec![0u8; nn_pre], vec![0.0f64; nn_pre], vec![false; nn_pre]);
        let natural_ring = self.p.de.as_ref().and_then(|d| d.chain.as_ref()).map(|c| c.natural_ring == Some(true)).unwrap_or(false);
        let clo = self.closure_field(px0, py0, n);
        // lane CARRY ct_sharp_limb_lp: the closure on a margin of the domain (the along-contour low-pass reads past the ring; a pure
        // function of the post position) and the T1b bicubic surface's strike
        let limb_sig = pre_sf.and(self.p.ct_sharp.as_ref().filter(|s| s.limb_lp.is_some())).map(|s| s.limb_lp_m.unwrap_or(crate::sharp::LIMB_LP_M) / POST_M);
        let limb_m = limb_sig.map(|g| (2.5 * g).ceil() as i64 + 2).unwrap_or(0);
        let ne = n + 2 * limb_m as usize;
        let clo_ext = if limb_sig.is_some() { self.closure_field(px0 - limb_m, py0 - limb_m, ne) } else { Vec::new() };
        let mut limb_lp_n = 0u64;
        // lane CARRY2 ct_sharp_aa: the pre-staircase surface on a 1-post margin of the domain as a pure function of position (B1 cut +
        // B4 closure (0 on corridor / sea posts) + lobe relief), bilinear between posts at the footprint sample points
        let aa_w = pre_sf.and(self.p.ct_sharp.as_ref().filter(|s| s.aa.is_some())).map(|s| s.aa_w.unwrap_or(crate::sharp::AA_W_POSTS));
        let na = n + 2;
        let z_aa: Vec<f64> = if aa_w.is_some() {
            let clo_m = self.closure_field(px0 - 1, py0 - 1, na);
            let lw = self.v2.as_ref().and_then(|v| v.lobes.as_ref()).filter(|l| l.relief.is_some());
            (0..na * na).map(|q| {
                let (gx, gy) = (px0 - 1 + (q % na) as i64, py0 - 1 + (q / na) as i64);
                let pv = self.post(gx, gy);
                pv.cut + if pv.fixed_relax() { 0.0 } else { clo_m[q] } + lw.map(|l| self.relief_post(l, &pv, gx as f64 * POST_M, gy as f64 * POST_M)).unwrap_or(0.0)
            }).collect()
        } else { Vec::new() };
        let mut aa_n = 0u64;
        let ringflow = self.p.de.as_ref().and_then(|d| d.ringflow);
        let v3 = self.p.v3p().cloned();
        let mut pvs: Vec<PostVal> = if v3.is_some() { Vec::with_capacity(n * n) } else { Vec::new() };
        for j in 0..n { for i in 0..n {
            let k = j * n + i;
            let pv = self.post(px0 + i as i64, py0 + j as i64);
            if v3.is_some() { pvs.push(pv); }
            let ring = i == 0 || j == 0 || i == n - 1 || j == n - 1;
            // lane RILL de_b5_chain_ring=1: the domain ring is a natural (zero-flux) boundary of B5, not a fixed one (as lane DE2's
            // de_fs ring for the DE solve), so neighbouring domains' relaxed fields agree in the seam band
            fixed[k] = pv.fixed_relax() || (ring && !natural_ring);
            z[k] = pv.cut + if fixed[k] { 0.0 } else { clo[k] };
            // lane LOBE de_lobe_relief: constructional relief of lane VE's lava flows on every non-sea / non-corridor post
            // (the domain ring included: it is a domain edge, not a fixed boundary); a pure function of (x, y)
            if let Some(lw) = self.v2.as_ref().and_then(|v| v.lobes.as_ref()).filter(|l| l.relief.is_some()) {
                z[k] += self.relief_post(lw, &pv, (px0 + i as i64) as f64 * POST_M, (py0 + j as i64) as f64 * POST_M);
            }
            // lane TREAD (ct_sharp_stage=pre): the 2 m staircase on the domain's initial surface (pointwise, so every domain gives a
            // shared post the same value); channel interiors and sea untouched, corridor banks re-profiled (as the post stage)
            if let Some(sf) = pre_sf {
                if !pv.fixed_sea && pv.hit.interior.is_none() {
                    let (xw, yw) = ((px0 + i as i64) as f64 * POST_M, (py0 + j as i64) as f64 * POST_M);
                    let mut ap = sf.apply_ta(xw, yw, z[k]);
                    let (mut limb_flag, mut limb_cl) = (false, 0.0f64);
                    // lane CARRY ct_sharp_limb_lp: a cliff-limb post re-evaluated on the closure low-passed along the T1b contour
                    if let (Some(sig), Some((_, 2, ta0))) = (limb_sig, ap) {
                        if ta0 > 0.0 && !fixed[k] {
                            let (ox, oy) = ((px0 - limb_m) as f64 * POST_M, (py0 - limb_m) as f64 * POST_M);
                            let val = |x: f64, y: f64| -> Option<f64> {
                                if x < 0.0 || y < 0.0 || x > (ne - 1) as f64 || y > (ne - 1) as f64 { return None; }
                                let (i0, j0) = ((x.floor() as usize).min(ne - 2), (y.floor() as usize).min(ne - 2));
                                let (fx, fy) = (x - i0 as f64, y - j0 as f64);
                                Some((1.0 - fy) * ((1.0 - fx) * clo_ext[j0 * ne + i0] + fx * clo_ext[j0 * ne + i0 + 1]) + fy * ((1.0 - fx) * clo_ext[(j0 + 1) * ne + i0] + fx * clo_ext[(j0 + 1) * ne + i0 + 1]))
                            };
                            let dir = |x: f64, y: f64| -> Option<(f64, f64)> {
                                let (wx, wy) = (ox + x * POST_M, oy + y * POST_M);
                                let gx = self.w.bicubic(wx + POST_M, wy) - self.w.bicubic(wx - POST_M, wy);
                                let gy = self.w.bicubic(wx, wy + POST_M) - self.w.bicubic(wx, wy - POST_M);
                                let g = gx.hypot(gy);
                                if g > 1e-9 { Some((-gy / g, gx / g)) } else { None }
                            };
                            let (ie, je) = ((i as i64 + limb_m) as f64, (j as i64 + limb_m) as f64);
                            if let Some(cl) = crate::de::contour_lowpass(ie, je, sig, &dir, &val) {
                                if let Some(a2) = sf.apply_ta(xw, yw, z[k] - clo[k] + cl) { ap = Some(a2); limb_lp_n += 1; limb_flag = true; limb_cl = cl; }
                            }
                        }
                    }
                    // lane CARRY2 ct_sharp_aa: the footprint mean of the staircase (SHARP_AA_RULE)
                    if let (Some(w), Some((zc, rg, ta)), false) = (aa_w, ap, fixed[k]) {
                        let zin = if limb_flag { z[k] - clo[k] + limb_cl } else { z[k] };
                        let sh = zin - z_aa[(j + 1) * na + i + 1];
                        let zpt = |x: f64, y: f64| -> f64 {
                            let (i0, j0) = ((x.floor().max(0.0) as usize).min(na - 2), (y.floor().max(0.0) as usize).min(na - 2));
                            let (fx, fy) = (x - i0 as f64, y - j0 as f64);
                            (1.0 - fy) * ((1.0 - fx) * z_aa[j0 * na + i0] + fx * z_aa[j0 * na + i0 + 1]) + fy * ((1.0 - fx) * z_aa[(j0 + 1) * na + i0] + fx * z_aa[(j0 + 1) * na + i0 + 1])
                        };
                        const DZ: f64 = 0.01;
                        let sig = match (sf.apply_ta(xw, yw, zin + DZ), sf.apply_ta(xw, yw, zin - DZ)) { (Some(a), Some(b)) => (a.0 - b.0) / (2.0 * DZ), _ => 1.0 };
                        let ns = crate::sharp::AA_SAMPLES;
                        let (mut sa, mut sz) = (0.0, 0.0);
                        for b in 0..ns { for a in 0..ns {
                            let (ou, ov) = (w * ((a as f64 + 0.5) / ns as f64 - 0.5), w * ((b as f64 + 0.5) / ns as f64 - 0.5));
                            let zp = zpt(i as f64 + 1.0 + ou, j as f64 + 1.0 + ov) + sh;
                            sa += sf.apply_ta(xw + ou * POST_M, yw + ov * POST_M, zp).map(|v| v.0).unwrap_or(zp);
                            sz += zp;
                        } }
                        let nn = (ns * ns) as f64;
                        let za = sa / nn - sig * (sz / nn - zin);
                        if (za - zc).abs() > 1e-4 { aa_n += 1; }
                        ap = Some((za, rg, ta));
                    }
                    if let Some((zn, rg, ta)) = ap {
                        pre_shift[k] = zn - z[k];
                        z[k] = zn;
                        pre_rg[k] = rg;
                        pre_ta[k] = ta;
                        riser[k] = rg == 2 && ta > 0.0;
                    }
                }
            }
            cells[k] = pv.cell;
            if ringflow.is_none() && (ring || i == 1 || j == 1 || i == n - 2 || j == n - 2) {
                extra[k] = self.w.area_bilinear((px0 + i as i64) as f64 * POST_M, (py0 + j as i64) as f64 * POST_M) / 16.0;
            }
        } }
        // lane STEEP-TALUS ct_sharp_talus_rock: talus posts steeper than talus can stand over 8 m (both up- and downslope) on the
        // pre-staircase domain surface are rock (regime 2)
        let talus_rock_n = match (pre_sf, self.p.ct_sharp.as_ref().filter(|s| s.talus_rock.is_some())) {
            (Some(_), Some(sp)) => {
                let thr = sp.talus_rock_deg.unwrap_or(crate::sharp::TALUS_ROCK_DEG).to_radians().tan();
                // not riser posts: B5 and the DE treat them as the bare rock they are (slope exposure, span critical slope with the
                // incoming-slope exemption), not as the talus mantle (lane STEEP-TALUS draft st1/st2 held them as risers: 35-40 m DE debris
                // piles and op-cap stalls in region R tiles 35_35 and 36_36)
                Some(crate::sharp::talus_rock_relabel(n, &z, &mut pre_rg, thr).len())
            }
            _ => None,
        };
        if let Some(mode) = ringflow { extra = self.ring_inflow(px0, py0, n, mode); }
        // lane RILL de_b5_chain: the coarse network's discharge crossing into the domain replaces the ring rule
        let chain = self.p.de.as_ref().and_then(|d| d.chain.as_ref()).map(|c| c.mode);
        if chain.is_some() {
            extra = self.chain_coarse().inflow(px0, py0, n, self.p.mfd_exponent);
            // diagnostic (log only): this domain's coarse-network inflow vs the T1b inward-flux rule (de_ringflow=1)
            let t1b: f64 = self.ring_inflow(px0, py0, n, 1).iter().sum();
            let ex: f64 = extra.iter().sum();
            let mx = extra.iter().cloned().fold(0.0f64, f64::max);
            eprintln!("chain inflow {ix} {iy}: coarse sum {ex:.4e} m2 (max post {mx:.3e}), T1b ringflow=1 sum {t1b:.4e} m2");
        }
        if self.p.relax_on && self.v2.is_some() {
            let v2 = self.v2.as_ref().unwrap();
            let oa: Vec<(crate::v2::OverAt, f64)> = (0..n * n).map(|k| {
                let (i, j) = ((k % n) as i64, (k / n) as i64);
                self.over_at(v2, (px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M)
            }).collect();
            let win = v3.as_ref().map(|w| w.relax_win_on).unwrap_or(false);
            // lane DE DE5: lobe columns of the domain's edifice posts (built once per domain)
            let lob: Vec<Option<crate::lobes::LobeAt>> = if v2.lobes.is_some() {
                (0..n * n).map(|k| { let (i, j) = ((k % n) as i64, (k / n) as i64); self.lobe_at(v2, (px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M, oa[k].0, oa[k].1) }).collect()
            } else { Vec::new() };
            let (sc3, k3) = (sc_of(crate::v2::LITHO_VOLCANIC), k_of(crate::v2::LITHO_VOLCANIC));
            // lobe window means memoized per post (recomputed when the surface moved > 0.25 m: the relaxation moves slowly)
            let lob_memo: RefCell<Vec<(f64, f64)>> = RefCell::new(vec![(f64::NAN, 0.0); if lob.is_empty() { 0 } else { n * n }]);
            let litho = |k: usize, zz: f64| -> (f64, f64) {
                if let Some(Some(la)) = lob.get(k) {
                    let u = la.u(zz);
                    if win {
                        let mut mm = lob_memo.borrow_mut();
                        let m = if (mm[k].0 - u).abs() <= 0.25 { mm[k].1 } else { let m = la.inv_mean(u - 30.0, u + 30.0); mm[k] = (u, m); m };
                        return (1.0 / m, k3 * sc3 * m);
                    }
                    let sc = la.sc(u); return (sc, k3 * sc3 / sc);
                }
                let (i, j) = ((k % n) as i64, (k / n) as i64);
                let (x, y) = ((px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M);
                if win { self.win_props_v3(v2, x, y, cells[k], zz, oa[k].0, oa[k].1) } else { self.bed_props_v2(v2, x, y, cells[k], zz, oa[k].0, oa[k].1) }
            };
            // lane RILL de_b5_chain=3: the coarse change, already mass-balanced window-wide (W6 rule on the coarse lattice), is
            // part of the initial surface (the domain's W6 does not re-balance it: no domain-dependent correction); ring included
            if chain == Some(3) {
                let c = self.chain_coarse();
                for k in 0..n * n { if !pvs[k].fixed_relax() { z[k] += c.dz_at(px0 + (k % n) as i64, py0 + (k / n) as i64); } }
            }
            let z_init = if v3.is_some() { z.clone() } else { Vec::new() };
            let mut pre_stats = serde_json::Value::Null;
            // lane RILL de_b5_chain=2: the coarse pass's change on every non-sea / non-corridor post, the domain ring included (it
            // is a domain edge: without the change there it held the pre-incision surface, a step the relaxation ramped inward
            // = tile-edge lines on G's coastal corner, binary c1); counts as relaxation change for W6
            if chain == Some(2) {
                let c = self.chain_coarse();
                for k in 0..n * n { if !pvs[k].fixed_relax() { z[k] += c.dz_at(px0 + (k % n) as i64, py0 + (k / n) as i64); } }
            }
            // lane DE DE3: bed K for the stream power (large drainage areas only)
            let bed_k = |k: usize, zz: f64| -> f64 {
                if let Some(Some(la)) = lob.get(k) { return k3 * sc3 / la.sc(la.u(zz)); }
                let (i, j) = ((k % n) as i64, (k / n) as i64);
                self.de_bed(v2, (px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M, cells[k], zz, oa[k].0, oa[k].1).1
            };
            let use_bk = self.p.de.as_ref().map(|d| d.sp_area_m2.is_some()).unwrap_or(false);
            // lane RIB: runoff-weighted rill drainage area on lobe-column posts (factor of the facies at the current surface)
            let rop = self.p.de.as_ref().and_then(|d| d.runoff.as_ref());
            // lane DE-TALUS de_talus_runoff: talus posts of the pre staircase (2 m regime 1) infiltrate (runoff factor, transmission loss)
            let rgp = self.p.de.as_ref().and_then(|d| d.regime.as_ref()).filter(|r| r.runoff.is_some() && pre_sf.is_some());
            // lane CLIMATE-RUNOFF climate_runoff=1: each post's own area also x the T1b runoff weight R / r_ref at the post
            let cw: Vec<f64> = if self.p.climate_runoff.is_some() {
                (0..n * n).map(|k| { let (i, j) = ((k % n) as i64, (k / n) as i64); self.w.runoff_weight_at((px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M).unwrap_or(1.0) }).collect()
            } else { Vec::new() };
            let runoff = |k: usize, zz: f64| -> (f64, f64) {
                let (f, t) = if let Some(v) = rgp.and_then(|r| r.runoff_of(pre_rg[k])) { v } else {
                    match (rop, lob.get(k)) { (Some(r), Some(Some(la))) => r.factor(la.fac_any(la.u(zz))), _ => (1.0, 1.0) }
                };
                if cw.is_empty() { (f, t) } else { (f * cw[k], t) }
            };
            let use_ro = (rop.is_some() && !lob.is_empty()) || rgp.is_some() || !cw.is_empty();
            // lane RILL de_b5_chain_snap=1 (modes 2/3): the coarse inflow enters at the fine valley floor of the conditioned surface
            if let Some(c) = self.p.de.as_ref().and_then(|d| d.chain.as_ref()).filter(|c| c.snap == Some(true) && c.mode >= 2) {
                crate::chain::snap_to_floor(n, &mut extra, &z, c.cell_posts);
            }
            // lane TREAD: riser posts are held in B5 (not crept or incised)
            let fixed_b5: Vec<bool> = if pre_sf.is_some() { (0..n * n).map(|k| fixed[k] || riser[k]).collect() } else { Vec::new() };
            let z_b5 = if pre_sf.is_some() { z.clone() } else { Vec::new() };
            // lane CARRY ct_sharp_carry_strike: riser contour directions of the pre-B5 surface (smoothed 3 m), used by both carries
            let strike: Option<crate::de::Strike> = pre_sf.and(self.p.ct_sharp.as_ref().filter(|s| s.carry_strike.is_some()))
                .map(|s| crate::de::Strike::of_surface(n, &z_b5, 1.5, s.carry_strike_m.unwrap_or(crate::sharp::CARRY_STRIKE_M) / POST_M));
            // lane CARRY2 ct_sharp_riser_noflux: the held risers are a zero-flux boundary of the B5 creep
            let noflux: Option<Vec<bool>> = pre_sf.and(self.p.ct_sharp.as_ref().filter(|s| s.riser_noflux.is_some())).map(|_| (0..n * n).map(|k| riser[k] && !fixed[k]).collect());
            let _st: RelaxStats = crate::relax::relax_dx_nf(n, POST_M, &mut z, if pre_sf.is_some() { &fixed_b5 } else { &fixed }, noflux.as_deref(), &extra, &litho, if use_bk { Some(&bed_k) } else { None }, if use_ro { Some(&runoff) } else { None }, &self.p);
            // lane TREAD: the held risers move with the B5 change around them (shape kept, no one-post pillars / pits)
            if pre_sf.is_some() {
                let mv: Vec<bool> = (0..n * n).map(|k| riser[k] && !fixed[k]).collect();
                let src: Vec<bool> = (0..n * n).map(|k| !fixed_b5[k]).collect();
                // lane CARRY2 ct_sharp_carry_fringe: the riser fringe's B5 change low-passed along the contour first
                if let (Some(st), true) = (strike.as_ref(), self.p.ct_sharp.as_ref().map(|s| s.carry_fringe.is_some()).unwrap_or(false)) {
                    let fr = crate::de::riser_fringe(n, &mv, &fixed_b5, crate::sharp::FRINGE_POSTS);
                    crate::de::carry_fringe_strike(n, &mut z, &z_b5, &fr, st);
                }
                match strike.as_ref() { Some(st) => crate::de::carry_risers_strike(n, &mut z, &z_b5, &mv, &src, 2, st), None => crate::de::carry_risers(n, &mut z, &z_b5, &mv, &src, 2) }
            }
            // lane DE DE1: layered hillslope equilibrium on the relaxed surface
            let mut eta: Vec<f64> = Vec::new();
            if let Some(dp) = self.p.de.as_ref() {
                let t_de = Instant::now();
                let z0 = z.clone();
                let per = v2.pinch.as_ref().map(|q| q.periods).unwrap_or([1.0, 1.0]);
                // per-post column handles: (column, datum, pinch amplitudes) or a constant Sc (deposit / no column)
                let cols: Vec<DeCol> = (0..n * n).map(|k| {
                    let (i, j) = ((k % n) as i64, (k / n) as i64);
                    let (x, y) = ((px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M);
                    let (o, s) = oa[k];
                    if let Some(Some(_)) = lob.get(k) { return DeCol::Lobe(k); }
                    match o.class() {
                        crate::v2::OVER_DEPOSIT => DeCol::Const(v2.v.deposit_sc),
                        crate::v2::OVER_EDIFICE => { let (cid, ec) = self.ed_dom(v2, o); DeCol::Col(ec, s, self.pinch_amps(v2, cid, x, y)) }
                        _ => match self.strata.cells.get(&cells[k]) {
                            Some(cs) => DeCol::Col(cs, self.strata.datum(cs, x, y), self.pinch_amps(v2, cells[k], x, y)),
                            None => DeCol::Const(sc_of(0)),
                        },
                    }
                }).collect();
                let uw = |k: usize, zz: f64| -> f64 {
                    match &cols[k] { DeCol::Col(cs, dat, Some(a)) => -crate::v3::warp(&cs.col, a, &per, dat - zz, crate::v3::UP_H), DeCol::Col(_, dat, None) => zz - dat, DeCol::Const(..) => 0.0,
                        DeCol::Lobe(q) => lob[*q].as_ref().unwrap().u(zz) }
                };
                let mut scw = vec![1.0f64; n * n];
                let mut ex = vec![0.0f64; n * n];
                let mut bd = vec![0.0f64; n * n];
                // lane DE2 B2: with de_fs the domain ring is a natural (zero-flux) boundary of the DE solve, not a fixed one, so
                // neighbouring tiles' solutions agree in the seam band; channel corridors and sea stay fixed
                let fixed_de: Vec<bool> = if dp.fs.as_ref().map(|f| f.natural_ring).unwrap_or(false) { (0..n * n).map(|k| pvs[k].fixed_relax()).collect() } else { fixed.clone() };
                for k in 0..n * n {
                    if fixed_de[k] { continue; }
                    let (i, j) = ((k % n) as i64, (k / n) as i64);
                    let (x, y) = ((px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M);
                    let h = 1.0;
                    let g0 = ((self.w.bicubic(x + h, y) - self.w.bicubic(x - h, y)) / (2.0 * h)).hypot((self.w.bicubic(x, y + h) - self.w.bicubic(x, y - h)) / (2.0 * h));
                    // lane DE2 F1: formation-scale window from T1b's horizontal resolution and the local base slope
                    let win_m = match (dp.fs.as_ref(), dp.window_m) { (Some(f), _) => Some((f.l_h_m * g0).clamp(f.w_min_m, f.w_max_m)), (None, w) => w };
                    let hw = win_m.map(|w| 0.5 * w).unwrap_or(30.0);
                    scw[k] = match &cols[k] {
                        DeCol::Col(cs, _, _) => match win_m { None => crate::strata::sc_win_at(cs, uw(k, z0[k])), Some(_) => { let u = uw(k, z0[k]); 1.0 / crate::strata::inv_sc_mean(cs, u - hw, u + hw) } },
                        DeCol::Const(s) => *s,
                        DeCol::Lobe(q) => { let la = lob[*q].as_ref().unwrap(); let u = la.u(z0[k]); 1.0 / la.inv_mean(u - hw, u + hw) } };
                    ex[k] = crate::de::exposure(g0, scw[k], dp);
                    // lane DE-TALUS de_talus_mantle: talus (and mode 2: tread) posts of the pre staircase hide the bed contrast
                    if let (Some(rp), true) = (dp.regime.as_ref(), pre_sf.is_some()) { if let Some(e) = rp.expose_override(pre_rg[k]) { ex[k] = e; } }
                    let r = v2.aux.bilinear(&v2.aux.relief, x, y).unwrap_or(0.0);
                    bd[k] = (dp.amp_frac * r).min(dp.amp_cap_m);
                }
                let scb = |k: usize, zz: f64| -> f64 {
                    match &cols[k] { DeCol::Col(cs, _, _) => cs.col.bed_sc[cs.col.bed_at_depth(-uw(k, zz))], DeCol::Const(s) => *s,
                        DeCol::Lobe(q) => { let la = lob[*q].as_ref().unwrap(); la.sc(la.u(zz)) } }
                };
                let iu = |k: usize, u1: f64, u2: f64| -> f64 {
                    match &cols[k] { DeCol::Col(cs, _, _) => crate::strata::inv_sc_mean(cs, u1, u2), DeCol::Const(s) => 1.0 / *s,
                        DeCol::Lobe(q) => lob[*q].as_ref().unwrap().inv_mean(u1, u2) }
                };
                let talus_m: Vec<bool> = match (dp.regime.as_ref(), pre_sf.is_some()) { (Some(rp), true) => pre_rg.iter().map(|&g| rp.talus_fail(g)).collect(), _ => Vec::new() };
                let rf = crate::de::Risers { mask: &riser, tan: &pre_ta, talus: if talus_m.iter().any(|&b| b) { Some(&talus_m) } else { None }, strike: strike.as_ref() };
                let o = crate::de::de_solve_rf(&crate::de::DeInput { n, z0: &z0, fixed: &fixed_de, sc_w: &scw, expose: &ex, bound: &bd, sc_bed: &scb, uw: &uw, inv_u: &iu }, dp, if pre_sf.is_some() { Some(&rf) } else { None });
                if pre_sf.is_some() {
                    // riser rockfall bookkeeping: riser posts lowered by the DE (carry with the equilibrium change + rockfall in the failure step; m3, posts)
                    let (mut vol, mut cnt, mut dmax) = (0.0f64, 0u64, 0.0f64);
                    for k in 0..n * n { if riser[k] && !fixed_de[k] { let d = z0[k] - o.z[k] + o.deposit[k]; if d > 1e-3 { vol += d * POST_M * POST_M; cnt += 1; dmax = dmax.max(d); } } }
                    pre_stats = serde_json::json!({"riser_posts": riser.iter().filter(|&&r| r).count(), "riser_lowered_posts": cnt, "riser_lowered_m3": vol, "riser_lowered_max_m": dmax,
                        "deposit_posts_gt_0_05m": o.deposit.iter().filter(|&&v| v > 0.05).count(), "deposit_max_m": o.deposit.iter().cloned().fold(0.0, f64::max)});
                }
                let (zn, st) = (o.z, o.stats);
                eta = (0..n * n).map(|k| zn[k] - z0[k]).collect();
                // lane DE2 F1 square control: remove the DE change's own low-pass (masked, free posts)
                let mut zn = zn;
                if let Some(f) = dp.fs.as_ref() {
                    let fr: Vec<bool> = (0..n * n).map(|k| !fixed_de[k]).collect();
                    let lp = crate::v3::masked_lowpass(n, &eta, &fr, f.mb_sigma_m / POST_M);
                    for k in 0..n * n { if !fixed_de[k] { eta[k] -= lp[k]; zn[k] = z0[k] + eta[k]; } }
                }
                let talus: Vec<f64> = (0..n * n).map(|k| if o.deposit[k] > 0.0 { o.deposit[k] } else if o.mantle[k] { -1.0 } else { 0.0 }).collect();
                z = zn;
                let mut sj = serde_json::to_value(&st).unwrap();
                sj["seconds"] = serde_json::json!(t_de.elapsed().as_secs_f64());
                self.de_cache.lock().unwrap().insert((ix, iy), (eta.clone(), talus.clone(), sj));
            }
            // v3 W6 local mass balance: re-centre the relaxation's change (and the remap's displacement on the relaxed
            // surface) by its masked Gaussian low-pass over free posts
            if let Some(w) = v3.as_ref().filter(|w| w.massbal_on) {
                let free: Vec<bool> = (0..n * n).map(|k| !fixed[k]).collect();
                let dd: Vec<f64> = (0..n * n).map(|k| {
                    if fixed[k] { return 0.0; }
                    let mut d = z[k] - z_init[k];
                    // lane DE: the DE response is excluded from W6 unless de_mb=1
                    if !eta.is_empty() && self.p.de.as_ref().and_then(|d| d.mb) != Some(true) { d -= eta[k]; }
                    if w.massbal_remap_on && v2.v.remap_after_relax {
                        let (i, j) = ((k % n) as i64, (k / n) as i64);
                        d += self.remap_v2(v2, &pvs[k], (px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M, z[k]);
                    }
                    d
                }).collect();
                let lp = crate::v3::masked_lowpass(n, &dd, &free, w.massbal_sigma_m / POST_M);
                for k in 0..n * n { if !fixed[k] { z[k] -= lp[k]; } }
            }
            if let Some(cnt) = talus_rock_n { pre_stats["talus_rock_posts"] = serde_json::json!(cnt); }
            if limb_sig.is_some() { pre_stats["limb_lp_posts"] = serde_json::json!(limb_lp_n); }
            if aa_w.is_some() { pre_stats["aa_posts"] = serde_json::json!(aa_n); }
            if pre_sf.is_some() {
                self.sharp_cache.lock().unwrap().insert((ix, iy), (pre_shift.iter().map(|&v| crate::to_mm(v)).collect(), pre_rg.clone(), pre_stats));
            }
        } else if self.p.relax_on {
            let litho = |k: usize, zz: f64| -> (f64, f64) {
                let (i, j) = ((k % n) as i64, (k / n) as i64);
                match self.strata.cells.get(&cells[k]) {
                    Some(cs) => self.strata.bed_props(cs, (px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M, zz),
                    None => (sc_of(0), k_of(0)),
                }
            };
            // lane CLIMATE-RUNOFF climate_runoff=1: own area x the T1b runoff weight at the post (no transmission loss)
            let cw: Vec<f64> = if self.p.climate_runoff.is_some() {
                (0..n * n).map(|k| { let (i, j) = ((k % n) as i64, (k / n) as i64); self.w.runoff_weight_at((px0 + i) as f64 * POST_M, (py0 + j) as f64 * POST_M).unwrap_or(1.0) }).collect()
            } else { Vec::new() };
            let ro_c = |k: usize, _z: f64| (cw[k], 1.0);
            let _st: RelaxStats = relax(n, &mut z, &fixed, &extra, &litho, None, if cw.is_empty() { None } else { Some(&ro_c) }, &self.p);
        }
        let secs = t.elapsed().as_secs_f64();
        self.relax_cache.lock().unwrap().insert((ix, iy), (z.clone(), secs));
        (z, secs)
    }
}

/// One finished tile.
pub struct TileOut {
    pub ix: i64,
    pub iy: i64,
    pub height_mm: Vec<i32>,
    pub base_mm: Vec<i32>,
    pub litho: Vec<u8>,
    pub layer: Vec<u8>,
    /// 0 free, 1 channel interior, 2 corridor (bank), 3 fixed sea/lake/low
    pub fixed: Vec<u8>,
    pub closure_mm: Vec<i32>,
    pub cut_minus_base_m: Vec<f32>,
    pub base_mean_m: f64,
    /// Mean of the plain T1b bicubic (equals base_mean_m unless deperch_on).
    pub raw_t1b_mean_m: f64,
    pub mean_shift_c_m: f64,
    /// Tile-mean deficit before the B8 shift (base mean - current mean, m) and the largest first-shift raise of a post (m).
    pub mean_deficit_m: f64,
    /// Lane VALLEY-FLOOR-FIX (keys on only): tile-mean levee raise (m) and whether the apron cap was infeasible here.
    pub levee_mean_m: Option<f64>,
    pub carve_mean_m: Option<f64>,
    pub b8_cap_infeasible: Option<bool>,
    pub shift_max_m: f64,
    /// T2D X2 bookkeeping: free posts raised by the bank-lip grading and the largest raise (mm).
    pub lip_posts: u64,
    pub lip_raise_mm: i64,
    pub fill_raised_posts: u64,
    pub fill_max_raise_mm: i64,
    pub relax_seconds_own: f64,
    pub assemble_seconds: f64,
    /// v2 only (empty in v1): exposed bed index, over class (0 none, 1 edifice, 2 deposit), remap displacement (mm).
    pub bed: Vec<u16>,
    pub over: Vec<u8>,
    pub remap_mm: Vec<i32>,
    /// Lane PV (prov on only; empty otherwise): stable id of the exposed unit per post and the record of every id present.
    pub prov_id: Vec<u32>,
    pub prov_rows: Vec<serde_json::Value>,
    /// Lane BS (unit_depth=1 only; empty otherwise): depth below the exposed unit's top contact, height above its base (m).
    pub depth_in_unit_m: Vec<f32>,
    pub height_above_unit_base_m: Vec<f32>,
    /// Lane BS2 (bed_meta=1 only; empty otherwise): facies code and generated unit index per post.
    pub bed_style_id: Vec<u8>,
    pub flow_unit: Vec<u16>,
    /// Lane PYRO (pyro_sub=1 only; empty otherwise): pyroclastic subtype code per post and the tile's code counts by source.
    pub pyro_sub: Vec<u8>,
    pub pyro_stats: Option<serde_json::Value>,
    /// Lane DE (de=1 only; empty otherwise): the DE change after the seam blend (mm) and the solve / seam stats.
    pub de_mm: Vec<i32>,
    pub de_talus: Vec<u8>,
    pub de_stats: Option<serde_json::Value>,
    /// Lane SHARP (ct_sharp=1 only; empty otherwise): height change of the 2 m re-profile (mm), 2 m regime, and stats.
    pub sharp_mm: Vec<i32>,
    pub ct_regime: Vec<u8>,
    pub sharp_stats: Option<serde_json::Value>,
}

/// Seam weight of a tile for a post at chart coordinate v (m) on one axis.
#[inline]
fn seam_w(v: f64, tile_lo: f64, hw: f64) -> f64 {
    let d = (v - tile_lo).min(tile_lo + TILE_M - v);
    smoothstep((d + hw) / (2.0 * hw))
}

pub fn assemble(ctx: &Ctx, ix: i64, iy: i64) -> TileOut {
    if ctx.w.evolved.is_some() { return crate::contrast::assemble(ctx,ix,iy); }
    let p = &ctx.p;
    let (own, own_secs) = ctx.relaxed(ix, iy);
    let t = Instant::now();
    let h = p.halo_posts as i64;
    let n = TILE_POSTS + 2 * p.halo_posts;
    let (ox, oy) = tile_origin_post(ix, iy);
    let np = TILE_POSTS * TILE_POSTS;
    let mut z = vec![0.0f64; np];
    let mut pvs = Vec::with_capacity(np);
    let hw = p.seam_half_width_m;
    // neighbour relaxations, fixed (iy, ix) order
    let mut neigh: BTreeMap<(i64, i64), Vec<f64>> = BTreeMap::new();
    if p.relax_on {
        for dy in -1..=1 { for dx in -1..=1 {
            let k = (iy + dy, ix + dx);
            neigh.insert(k, if dx == 0 && dy == 0 { own.clone() } else { ctx.relaxed(ix + dx, iy + dy).0 });
        } }
    }
    let closure_own = ctx.closure_field(ox, oy, TILE_POSTS);
    // lane DE: neighbour DE changes (same order and weights as the relaxed blend) and seam spreads
    let de_on = p.de.is_some() && p.relax_on && ctx.v2.is_some();
    let de_n: BTreeMap<(i64, i64), (Vec<f64>, Vec<f64>, serde_json::Value)> = if de_on {
        neigh.keys().map(|&(ty, tx)| ((ty, tx), ctx.de_cache.lock().unwrap().get(&(tx, ty)).cloned().unwrap_or_else(|| (vec![0.0; n * n], vec![0.0; n * n], serde_json::Value::Null)))).collect()
    } else { BTreeMap::new() };
    let mut de_mm = if de_on { vec![0i32; np] } else { Vec::new() };
    let mut de_talus = if de_on { vec![0u8; np] } else { Vec::new() };
    let (mut spread_z, mut spread_e): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
    for j in 0..TILE_POSTS { for i in 0..TILE_POSTS {
        let (px, py) = (ox + i as i64, oy + j as i64);
        let pv = ctx.post(px, py);
        let (x, y) = (px as f64 * POST_M, py as f64 * POST_M);
        let v = if p.relax_on {
            let (mut acc, mut ws) = (0.0, 0.0);
            for (&(ty, tx), f) in neigh.iter() {
                let wgt = seam_w(x, tx as f64 * TILE_M, hw) * seam_w(y, ty as f64 * TILE_M, hw);
                if wgt <= 0.0 { continue; }
                let (dox, doy) = tile_origin_post(tx, ty);
                let (di, dj) = (px - (dox - h), py - (doy - h));
                acc += wgt * f[dj as usize * n + di as usize];
                ws += wgt;
            }
            if de_on {
                let (mut ea, mut lo, mut hi, mut elo, mut ehi, mut cnt) = (0.0, f64::MAX, f64::MIN, f64::MAX, f64::MIN, 0);
                let (mut tdep, mut tman) = (0.0f64, 0.0f64);
                for (&(ty, tx), f) in neigh.iter() {
                    let wgt = seam_w(x, tx as f64 * TILE_M, hw) * seam_w(y, ty as f64 * TILE_M, hw);
                    if wgt <= 0.0 { continue; }
                    let (dox, doy) = tile_origin_post(tx, ty);
                    let q = (py - (doy - h)) as usize * n + (px - (dox - h)) as usize;
                    let e = de_n[&(ty, tx)].0[q];
                    let tv = de_n[&(ty, tx)].1[q];
                    ea += wgt * e; cnt += 1;
                    if tv > 0.0 { tdep += wgt * tv; } else if tv < 0.0 { tman += wgt; }
                    lo = lo.min(f[q]); hi = hi.max(f[q]); elo = elo.min(e); ehi = ehi.max(e);
                }
                de_mm[j * TILE_POSTS + i] = to_mm(ea / ws);
                // talus class: 2 deposit >= 0.25 m (failure debris), 1 mantle (debris-covered slope below a bare hard face)
                de_talus[j * TILE_POSTS + i] = if tdep / ws >= 0.25 { 2 } else if tman / ws >= 0.5 { 1 } else { 0 };
                if cnt > 1 { spread_z.push(hi - lo); spread_e.push(ehi - elo); }
            }
            acc / ws
        } else {
            pv.cut + if pv.fixed_relax() { 0.0 } else { closure_own[j * TILE_POSTS + i] }
        };
        z[j * TILE_POSTS + i] = v;
        pvs.push(pv);
    } }
    // B6 carve
    let mut fixed = vec![0u8; np];
    let mut remap_mm = Vec::new();
    // Lane VALLEY-FLOOR-FIX b8_levee: the surface just before the carve (carved volume bookkeeping).
    let lev_on = ctx.v2.is_some() && p.vff.as_ref().map(|q| q.levee).unwrap_or(false);
    let mut zpre: Vec<f64> = if lev_on { vec![f64::NAN; np] } else { Vec::new() };
    if let Some(v2) = &ctx.v2 {
        // V5-V7 pointwise remap after the blend, then V3 carve with graded banks
        remap_mm = vec![0i32; np];
        for k in 0..np {
            let pv = &pvs[k];
            let (px, py) = (ox + (k % TILE_POSTS) as i64, oy + (k / TILE_POSTS) as i64);
            let (x, y) = (px as f64 * POST_M, py as f64 * POST_M);
            if pv.fixed_sea { fixed[k] = 3; z[k] = pv.base; continue; }
            if v2.v.remap_after_relax {
                let d = ctx.remap_v2(v2, pv, x, y, z[k]);
                z[k] += d;
                remap_mm[k] = to_mm(d);
            }
            // lane INTS2n4 arm_cl: the arm's carve and banks are booked from the B1 base (the T1b surface with its 32 m sea trench),
            // so the trench T1b already holds is not counted twice (ARM_CL_RULE)
            if lev_on { zpre[k] = if pv.hit.arm { pv.base } else { z[k] }; }
            if p.channels_on {
                // carving only lowers: a perched T1b bed (above the ground beside it) is not built up as a ridge
                if let Some(v) = pv.hit.interior { if v <= z[k] { z[k] = v; fixed[k] = 1; continue; } }
                // bank cone: alluvium limiting slope over one half-width, then the exposed class's limiting slope;
                // applied within bank_reach_m of the channel edge, faded over its last third (no step at the reach)
                let (l, c2, _) = ctx.exposed_v2(v2, x, y, pv.cell, z[k]);
                let sc_wall = if l == 254 { v2.v.deposit_sc } else if c2 == 255 { sc_of(0) } else { sc_of(c2) };
                let (sa, br) = (v2.v.bank_sc_alluvium, v2.v.bank_reach_m);
                let mut low = 0.0f64;
                for q in v2.curves.near_lines(x, y, 250.0) {
                    if q.r <= q.hw || q.r - q.hw > br { continue; }
                    let la = q.hw.max(2.0);
                    let c = q.bed + sa * (q.r - q.hw).min(la) + sc_wall * (q.r - q.hw - la).max(0.0);
                    let fade = 1.0 - smoothstep((q.r - q.hw - 2.0 * br / 3.0) / (br / 3.0));
                    // post-revision fix B3 (logged): the bank cone trims the channel's own banks, at most
                    // max(2 m, 3 x bankfull depth); it does not cut T1b's one-post cliffs back (beach tile 0.-1.24 lost
                    // 38 m there and the F6 mean shift put it back as a +28 m square on the few free posts)
                    let cap = (3.0 * q.depth).max(2.0);
                    low = low.max(fade * (z[k] - c).min(cap));
                }
                if low > 0.0 { z[k] -= low; fixed[k] = 2; } else if pv.hit.corridor { fixed[k] = 2; }
            }
        }
    }
    for k in 0..np {
        if ctx.v2.is_some() { break; }
        let pv = &pvs[k];
        if pv.fixed_sea { fixed[k] = 3; z[k] = pv.base; continue; }
        if p.channels_on {
            if let Some(v) = pv.hit.interior { z[k] = v; fixed[k] = 1; }
            else if let Some(b) = pv.hit.bank { z[k] = z[k].min(b); fixed[k] = 2; }
            else if pv.hit.corridor { fixed[k] = 2; }
        }
    }
    // Lane VALLEY-FLOOR-FIX b8_levee_v1: the carved cross-section returned beside each line as natural levees (B8_LEVEE_RULE).
    let mut levee_mean: Option<f64> = None;
    let carve_mean: Option<f64> = if lev_on { Some((0..np).filter(|&k| zpre[k].is_finite()).map(|k| zpre[k] - z[k]).sum::<f64>() / np as f64) } else { None };
    if let (Some(q), Some(v2)) = (p.vff.as_ref().filter(|q| q.levee), ctx.v2.as_ref()) {
        let mut sum = 0.0f64;
        for k in 0..np {
            if fixed[k] == 1 || fixed[k] == 3 { continue; }
            let (px, py) = (ox + (k % TILE_POSTS) as i64, oy + (k / TILE_POSTS) as i64);
            let l = v2.curves.levee_at(px as f64 * POST_M, py as f64 * POST_M, q, &|xx, yy| ctx.w.bicubic(xx, yy));
            if l > 0.0 { z[k] += l; sum += l; }
        }
        levee_mean = Some(sum / np as f64);
    }
    let capfix = p.vff.as_ref().map(|q| q.capfix).unwrap_or(false);
    let mut cap_infeasible: Option<bool> = None;
    let mut edge_only = false;
    // B8 mean
    let base_mean = pvs.iter().map(|v| v.base).sum::<f64>() / np as f64;
    let raw_t1b_mean = (0..np).map(|k| ctx.w.bicubic((ox + (k % TILE_POSTS) as i64) as f64 * POST_M, (oy + (k / TILE_POSTS) as i64) as f64 * POST_M)).sum::<f64>() / np as f64;
    let mut c = 0.0;
    let (mut deficit, mut shift_max) = (0.0f64, 0.0f64);
    let mut wcap = f64::INFINITY;
    let (mut lip_posts, mut lip_raise_mm) = (0u64, 0i64);
    // lane DE2 B1: with de_fs the tile-mean correction goes where the DE moved material (weight 0.15 + |DE change| / its tile
    // mean, at most 4), not as a tile-shaped bump (the per-tile square)
    let dew: Vec<f64> = match (&de_mm, p.de.as_ref().and_then(|d| d.fs.as_ref())) {
        (dm, Some(_)) if !dm.is_empty() => { let m = dm.iter().map(|&v| (v as f64).abs()).sum::<f64>() / np as f64; if m > 0.0 { dm.iter().map(|&v| (0.15 + (v as f64).abs() / m).min(4.0)).collect() } else { Vec::new() } }
        _ => Vec::new(),
    };
    // lane LOBE de_lobe_relief: where the tile carries flow relief, the tile-mean correction goes where the relief moved
    // material (weight 0.15 + |R| / its tile mean, at most 4; as lane DE2 B1 for the DE change), not as a tile-shaped bump
    // (the relief is mass-neutral at the T1b scale but not per 256 m tile); empty (weight 1) without relief in the tile
    let dew: Vec<f64> = match ctx.v2.as_ref().and_then(|v| v.lobes.as_ref()).filter(|l| l.relief.is_some()) {
        Some(lw) => {
            let r: Vec<f64> = (0..np).map(|k| { let (px, py) = (ox + (k % TILE_POSTS) as i64, oy + (k / TILE_POSTS) as i64); ctx.relief_post(lw, &pvs[k], px as f64 * POST_M, py as f64 * POST_M).abs() }).collect();
            let m = r.iter().sum::<f64>() / np as f64;
            if m > 1e-3 { (0..np).map(|k| (0.15 + r[k] / m).min(4.0) * dew.get(k).copied().unwrap_or(1.0)).collect() } else { dew }
        }
        None => dew,
    };
    // b8_levee: the carve and its levees balance per channel segment (B8_LEVEE_RULE), so B8 restores the tile mean of every
    // other stage and leaves their tile imbalance (a transfer to the neighbour tiles along the same segment) in place.
    let shift_lc = match (carve_mean, levee_mean) { (Some(c), Some(l)) => c - l, _ => 0.0 };
    if p.mean_on {
        let cur = z.iter().sum::<f64>() / np as f64 + shift_lc;
        deficit = base_mean - cur;
        let tap = ctx.p.v3p().filter(|w| w.taper_c_on).map(|w| w.taper_c_m);
        let apr = ctx.p.v3p().filter(|w| w.apron_on).map(|w| (w.apron_base, w.apron_decay_m, w.taper_c_m));
        // lane RILL de_b5_chain_b8: the edge taper multiplies after the cap (taper_after); otherwise the B8 rule as before
        // (lane RILL2 b8_taper_after=1: the same rule without the chain)
        let taper_after = p.de.as_ref().and_then(|d| d.chain.as_ref()).map(|c| c.b8_taper_after_cap == Some(true)).unwrap_or(false)
            || ctx.p.v3p().map(|w| w.b8_taper_after == Some(true)).unwrap_or(false);
        let etap = |q: usize| -> f64 {
            let e = |a: usize| smoothstep(((a.min(TILE_POSTS - 1 - a)) as f64 * POST_M).min(p.mean_taper_m) / p.mean_taper_m);
            e(q % TILE_POSTS) * e(q / TILE_POSTS)
        };
        let tw = |q: usize| -> f64 {
            let mut t = match tap { Some(m) => smoothstep(pvs[q].hit.excess_m / m), None => 1.0 };
            // v3 W7b: residual volume put back as a footslope apron near channels (where the carve removed it)
            if let Some((b0, dk, tm)) = apr { t *= b0 + (-(pvs[q].hit.excess_m - tm).max(0.0) / dk).exp(); }
            t
        };
        // (same operation order as before when taper_after is off: e_x * e_y * t * dew)
        let phi = |q: usize| -> f64 { if taper_after { tw(q) * dew.get(q).copied().unwrap_or(1.0) } else { etap(q) * tw(q) * dew.get(q).copied().unwrap_or(1.0) } };
        let wt = |q: usize, w: f64| -> f64 { if taper_after { phi(q).min(w) * etap(q) } else { phi(q).min(w) } };
        // T2D X3: cap the apron raise at apron_cap x |deficit| (weight ceiling found by bisection, water-filling)
        if let Some(cap) = ctx.p.v3p().and_then(|w| w.apron_cap) {
            let dd = (base_mean - cur).abs();
            let f = |w: f64| -> f64 { let s: f64 = (0..np).filter(|&k| fixed[k] == 0).map(|k| wt(k, w)).sum(); if s > 0.0 { dd * np as f64 * w / s } else { f64::INFINITY } };
            let wmax = (0..np).filter(|&k| fixed[k] == 0).map(phi).fold(0.0f64, f64::max);
            if dd > 0.0 && wmax > 0.0 && f(wmax) > cap * dd {
                let (mut lo, mut hi) = (0.0f64, wmax);
                for _ in 0..60 { let m = 0.5 * (lo + hi); if f(m) > cap * dd { hi = m; } else { lo = m; } }
                wcap = if lo > 0.0 { lo } else { hi };
                // Lane VALLEY-FLOOR-FIX b8_capfix_v1: lo == 0 means even the uniform raise breaks the cap (infeasible); fall
                // back to the edge taper alone (no ceiling) instead of collapsing every weight to the tiny ceiling.
                if capfix { cap_infeasible = Some(lo <= 0.0); if lo <= 0.0 { wcap = f64::INFINITY; edge_only = true; } }
            } else if capfix { cap_infeasible = Some(false); }
        }
        // b8_capfix_v1 fallback (edge_only): the most spread-out raise that still fades at the tile edges (edge taper only)
        let wt = |q: usize, w: f64| -> f64 { if edge_only { etap(q) } else { wt(q, w) } };
        let sphi: f64 = (0..np).filter(|&k| fixed[k] == 0).map(|k| wt(k, wcap)).sum();
        if sphi > 0.0 {
            c = (base_mean - cur) * np as f64 / sphi;
            for k in 0..np { if fixed[k] == 0 { let r = c * wt(k, wcap); z[k] += r; shift_max = shift_max.max(r.abs()); } }
        }
    }
    // B10 quantize, then B9 fill in integer mm
    let fill = |hm: &mut Vec<i64>, raised_mask: &mut Vec<bool>| -> (u64, i64) {
        let (mut raised, mut max_raise) = (0u64, 0i64);
        let mut seen = vec![false; np];
        let mut pq: BinaryHeap<std::cmp::Reverse<(i64, usize)>> = BinaryHeap::new();
        for k in 0..np {
            let (i, j) = (k % TILE_POSTS, k / TILE_POSTS);
            // v2: channel interiors are not seeds (an epsilon fill grades flats in the carved bed toward the outlet, F1)
            // v2 post-revision fix B2 (logged): only sea/lake posts and the tile edge are seeds, so a carved interior or
            // bank post that is a pit is graded to its spill too (v1: every fixed post was a seed)
            let seed_fixed = if ctx.v2.is_some() { fixed[k] == 3 } else { fixed[k] != 0 };
            if i == 0 || j == 0 || i == TILE_POSTS - 1 || j == TILE_POSTS - 1 || seed_fixed {
                seen[k] = true;
                pq.push(std::cmp::Reverse((hm[k], k)));
            }
        }
        while let Some(std::cmp::Reverse((zc, k))) = pq.pop() {
            let (i, j) = ((k % TILE_POSTS) as i64, (k / TILE_POSTS) as i64);
            for (dx, dy) in [(1i64, 0i64), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)] {
                let (a, b) = (i + dx, j + dy);
                if a < 0 || b < 0 || a >= TILE_POSTS as i64 || b >= TILE_POSTS as i64 { continue; }
                let q = b as usize * TILE_POSTS + a as usize;
                if seen[q] { continue; }
                seen[q] = true;
                if hm[q] <= zc { let nz = zc + 1; max_raise = max_raise.max(nz - hm[q]); hm[q] = nz; raised += 1; raised_mask[q] = true; }
                pq.push(std::cmp::Reverse((hm[q], q)));
            }
        }
        (raised, max_raise)
    };
    let mut hm: Vec<i64> = z.iter().map(|&v| to_mm(v) as i64).collect();
    let (mut raised, mut max_raise) = (0u64, 0i64);
    let mut raised_mask = vec![false; np];
    // lane SHARP (ct_sharp=1): re-profile bare rock at 2 m with T1b's staircase remap (sharp.rs), pointwise, before the exposed
    // column is read (so litho / layer / bed / flow unit / prov / pyro are those at the final height) and before the B9 fill and
    // the B8 mean re-correction below (F varies in space, so it can open pits; the tile mean goes back to the B1 base).
    let (mut sharp_mm, mut ct_rg2, mut sharp_stats) = (Vec::new(), Vec::new(), None);
    let pre = ctx.p.ct_sharp.as_ref().map(|s| s.pre()).unwrap_or(false);
    if pre {
        // lane TREAD: the staircase was applied before B5 / DE (relaxed()); report the owner domain's shift and 2 m regime
        let sc = ctx.sharp_cache.lock().unwrap();
        let (sh, rg, st) = sc.get(&(ix, iy)).expect("pre-stage cache filled by relaxed()");
        let nd = TILE_POSTS + 2 * p.halo_posts;
        sharp_mm = vec![0i32; np];
        ct_rg2 = vec![0u8; np];
        let (mut touched, mut lowered, mut raised_by) = (0u64, 0u64, 0u64);
        for k in 0..np {
            let q = (k / TILE_POSTS + p.halo_posts) * nd + k % TILE_POSTS + p.halo_posts;
            sharp_mm[k] = sh[q];
            ct_rg2[k] = rg[q];
            if sh[q] != 0 { touched += 1; if sh[q] < 0 { lowered += 1; } else { raised_by += 1; } }
        }
        let cnt = |c: u8| ct_rg2.iter().filter(|&&v| v == c).count();
        sharp_stats = Some(serde_json::json!({"stage": "pre", "touched_posts": touched, "lowered_posts": lowered, "raised_posts": raised_by,
            "regime_posts": [cnt(0), cnt(1), cnt(2)], "domain": st}));
        // lane SHARP2: tread posts (2 m regime 3) only when the window carries the 3-limb fields (else every byte as before)
        if ctx.sharp.as_ref().map(|s| s.tread).unwrap_or(false) { sharp_stats.as_mut().unwrap()["regime_posts"].as_array_mut().unwrap().push(cnt(3).into()); }
    } else if let Some(sf) = &ctx.sharp {
        sharp_mm = vec![0i32; np];
        ct_rg2 = vec![0u8; np];
        let (mut touched, mut lowered, mut raised_by, mut dmin, mut dmax, mut dsum) = (0u64, 0u64, 0u64, 0i64, 0i64, 0i64);
        for k in 0..np {
            if fixed[k] == 1 || fixed[k] == 3 || pvs[k].fixed_sea { continue; }
            let (px, py) = (ox + (k % TILE_POSTS) as i64, oy + (k / TILE_POSTS) as i64);
            let Some((zn, rg)) = sf.apply(px as f64 * POST_M, py as f64 * POST_M, hm[k] as f64 / 1000.0) else { continue };
            let v = to_mm(zn) as i64;
            let dd = v - hm[k];
            touched += 1;
            if dd < 0 { lowered += 1; } else if dd > 0 { raised_by += 1; }
            dmin = dmin.min(dd); dmax = dmax.max(dd); dsum += dd;
            sharp_mm[k] = dd as i32;
            ct_rg2[k] = rg;
            hm[k] = v;
        }
        let cnt = |c: u8| ct_rg2.iter().filter(|&&v| v == c).count();
        sharp_stats = Some(serde_json::json!({"touched_posts": touched, "lowered_posts": lowered, "raised_posts": raised_by, "change_min_mm": dmin, "change_max_mm": dmax,
            "change_mean_mm": if touched > 0 { dsum as f64 / touched as f64 } else { 0.0 },
            "regime_posts": [cnt(0), cnt(1), cnt(2)]}));
        if sf.tread { sharp_stats.as_mut().unwrap()["regime_posts"].as_array_mut().unwrap().push(cnt(3).into()); }
    }
    if p.fill_on {
        let (r, m) = fill(&mut hm, &mut raised_mask);
        raised = r; max_raise = m;
        // v2 post-revision fix B1 (logged in PREREG): the fill runs after the B8 mean shift and can add volume (one flank
        // tile: 11 m over 3600 posts, F6 0.33 m). Re-correct the mean on the free posts the fill did not raise, then refill.
        if ctx.v2.is_some() && p.mean_on {
            let tap = ctx.p.v3p().filter(|w| w.taper_c_on).map(|w| w.taper_c_m);
            let apr = ctx.p.v3p().filter(|w| w.apron_on).map(|w| (w.apron_base, w.apron_decay_m, w.taper_c_m));
            let phi = |q: usize| -> f64 {
                let e = |a: usize| smoothstep(((a.min(TILE_POSTS - 1 - a)) as f64 * POST_M).min(p.mean_taper_m) / p.mean_taper_m);
                let mut t = match tap { Some(m) => smoothstep(pvs[q].hit.excess_m / m), None => 1.0 };
                // v3 W7b: residual volume put back as a footslope apron near channels (where the carve removed it)
                if let Some((b0, dk, tm)) = apr { t *= b0 + (-(pvs[q].hit.excess_m - tm).max(0.0) / dk).exp(); }
                if edge_only { return e(q % TILE_POSTS) * e(q / TILE_POSTS); }
                e(q % TILE_POSTS) * e(q / TILE_POSTS) * t * dew.get(q).copied().unwrap_or(1.0)
            };
            for _ in 0..3 {
                let cur = hm.iter().map(|&v| v as f64).sum::<f64>() / np as f64 / 1000.0;
                let d = base_mean - cur - shift_lc;
                if d.abs() <= 0.05 { break; }
                let sphi: f64 = (0..np).filter(|&k| fixed[k] == 0 && !raised_mask[k]).map(|k| phi(k).min(wcap)).sum();
                if sphi <= 0.0 { break; }
                let c2 = d * np as f64 / sphi;
                for k in 0..np { if fixed[k] == 0 && !raised_mask[k] { hm[k] = to_mm(hm[k] as f64 / 1000.0 + c2 * phi(k).min(wcap)) as i64; } }
                c += c2;
                let (r, m) = fill(&mut hm, &mut raised_mask);
                raised += r; max_raise = max_raise.max(m);
            }
        }
        // T2D bank fins: a corridor bank post (fixed 2) standing above both neighbours along x or y by more than the
        // threshold is a one-post levee left by the fixed corridor; it is lowered to the higher of those two neighbours
        if let Some(tf) = ctx.p.v3p().and_then(|w| w.bank_fin_m) {
            let tmm = (tf * 1000.0).round() as i64;
            let n = TILE_POSTS;
            let mut changed = false;
            for j in 1..n - 1 { for i in 1..n - 1 {
                let k = j * n + i;
                if fixed[k] != 2 { continue; }
                let (l, r, d, u) = (hm[k - 1], hm[k + 1], hm[k - n], hm[k + n]);
                let ex = (hm[k] - l).min(hm[k] - r);
                let ey = (hm[k] - d).min(hm[k] - u);
                let mut t = hm[k];
                if ex > tmm { t = t.min(l.max(r)); }
                if ey > tmm { t = t.min(d.max(u)); }
                if t < hm[k] { hm[k] = t; changed = true; }
            } }
            if changed { let (r, m) = fill(&mut hm, &mut raised_mask); raised += r; max_raise = max_raise.max(m); }
        }
        // T2D X2 (revision 1): bank lip grading on the F1 flow paths. From every carved interior post, follow steepest descent
        // (the checks.rs F1 walk) through corridor posts; where the path would step down onto a free post by at most the cap,
        // that free post is raised to the corridor post + 1 mm (the bank lip may not sit below the channel water it holds).
        // Larger drops (T1b perched channels) are left. Then refill. Up to 20 rounds.
        if let Some(gmax) = ctx.p.v3p().and_then(|w| w.bank_grade_max_m) {
            let gmm = (gmax * 1000.0).round() as i64;
            let n = TILE_POSTS as i64;
            for _ in 0..400 {
                let mut changed = false;
                for k0 in 0..np {
                    if fixed[k0] != 1 { continue; }
                    let mut cur = k0;
                    for _ in 0..4 * TILE_POSTS {
                        let (i, j) = ((cur % TILE_POSTS) as i64, (cur / TILE_POSTS) as i64);
                        if i == 0 || j == 0 || i == n - 1 || j == n - 1 { break; }
                        let zc = hm[cur] as f64;
                        let mut best = (0.0f64, usize::MAX);
                        for (dx, dy) in [(1i64, 0i64), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)] {
                            let q = (j + dy) as usize * TILE_POSTS + (i + dx) as usize;
                            let d = if dx != 0 && dy != 0 { std::f64::consts::SQRT_2 } else { 1.0 };
                            let sl = (zc - hm[q] as f64) / d;
                            if sl > best.0 { best = (sl, q); }
                        }
                        if best.1 == usize::MAX { break; }
                        let q = best.1;
                        if fixed[q] == 0 {
                            // never touch a tile-edge post (shared with the neighbour tile, F9)
                            let (qi, qj) = ((q % TILE_POSTS) as i64, (q / TILE_POSTS) as i64);
                            let edge = qi == 0 || qj == 0 || qi == n - 1 || qj == n - 1;
                            if !edge && hm[cur] + 1 - hm[q] <= gmm { lip_raise_mm = lip_raise_mm.max(hm[cur] + 1 - hm[q]); lip_posts += 1; hm[q] = hm[cur] + 1; changed = true; }
                            break;
                        }
                        if fixed[q] == 3 { break; }
                        cur = q;
                    }
                }
                if !changed { break; }
                let (r, m) = fill(&mut hm, &mut raised_mask);
                raised += r; max_raise = max_raise.max(m);
            }
        }
    }
    let mut litho = vec![0u8; np];
    let mut layer = vec![0u8; np];
    let mut closure_mm = vec![0i32; np];
    let mut cmb = vec![0f32; np];
    let (mut bed, mut over) = (Vec::new(), Vec::new());
    let ud = ctx.p.unit_depth == Some(true) && ctx.v2.is_some();
    let (mut diu, mut hab) = if ud { (vec![f32::NAN; np], vec![f32::NAN; np]) } else { (Vec::new(), Vec::new()) };
    let bm = ctx.p.beds.as_ref().and_then(|b| b.bed_meta) == Some(true) && ctx.v2.is_some();
    let (mut bsid, mut fu) = if bm { (vec![crate::bedstyle::FAC_WATER; np], vec![u16::MAX; np]) } else { (Vec::new(), Vec::new()) };
    // lane PYRO: subtype per post and counts [source: 0 lobe column, 1 edifice 1D sequence, 2 strata][code]
    let pyon = bm && ctx.p.pyro.is_some() && ctx.v2.as_ref().map(|v| v.pyro.is_some()).unwrap_or(false);
    let mut psub = if pyon { vec![crate::pyro::SUB_NONE; np] } else { Vec::new() };
    let mut pcount = [[0u64; 256]; 3];
    // lane PV: the cut's own layer index per post (kept-layer index; 253 edifice, 254 deposit, 255 none), sea = usize::MAX
    let mut lay: Vec<usize> = if ctx.prov.is_some() { vec![usize::MAX; np] } else { Vec::new() };
    if let Some(v2) = &ctx.v2 {
        bed = vec![0u16; np];
        over = vec![0u8; np];
        for k in 0..np {
            let (px, py) = (ox + (k % TILE_POSTS) as i64, oy + (k / TILE_POSTS) as i64);
            let (x, y) = (px as f64 * POST_M, py as f64 * POST_M);
            let (l, c2, b) = ctx.exposed_v2(v2, x, y, pvs[k].cell, hm[k] as f64 / 1000.0);
            if pvs[k].fixed_sea { litho[k] = 255; layer[k] = 255; continue; }
            litho[k] = c2; layer[k] = l.min(255) as u8; bed[k] = b;
            if !lay.is_empty() { lay[k] = l; }
            let (oo, ss) = ctx.over_at(v2, x, y);
            over[k] = oo.class();
            // lane DE DE5: edifice posts read the lobe column (unit id as the bed index; facies, unit and contacts of the lobe)
            let zk = hm[k] as f64 / 1000.0;
            if let Some(la) = ctx.lobe_at(v2, x, y, oo, ss) {
                let u = la.u(zk);
                if let Some((fac, unit)) = la.meta(u) {
                    bed[k] = unit;
                    if ud { let (a, bb) = la.unit_contacts(u).map(|q| (q.0 as f32, q.1 as f32)).unwrap_or((f32::NAN, f32::NAN)); diu[k] = a; hab[k] = bb; }
                    if bm { bsid[k] = fac; fu[k] = unit; }
                    if pyon { let c = la.pyro_code(u, x, y).unwrap_or(crate::pyro::SUB_UNRESOLVED); psub[k] = c; pcount[0][c as usize] += 1; }
                    continue;
                }
            }
            if ud { let (a, b) = ctx.unit_depth(v2, x, y, pvs[k].cell, zk); diu[k] = a; hab[k] = b; }
            if bm { let (a, b) = ctx.bed_meta_at(v2, x, y, pvs[k].cell, zk); bsid[k] = a; fu[k] = b; }
            if pyon { let c = ctx.pyro_meta_at(v2, x, y, pvs[k].cell, zk); psub[k] = c; pcount[if oo.class() == crate::v2::OVER_EDIFICE { 1 } else { 2 }][c as usize] += 1; }
        }
    }
    for k in 0..np {
        let (px, py) = (ox + (k % TILE_POSTS) as i64, oy + (k / TILE_POSTS) as i64);
        if ctx.v2.is_none() {
        let (l, c2) = ctx.exposed(px, py, pvs[k].cell, hm[k] as f64 / 1000.0);
        litho[k] = c2;
        layer[k] = l.min(254) as u8;
        if !lay.is_empty() && !pvs[k].fixed_sea { lay[k] = l; }
        }
        closure_mm[k] = if pvs[k].fixed_relax() { 0 } else { to_mm(closure_own[k]) };
        cmb[k] = (pvs[k].cut - pvs[k].base) as f32;
    }
    let de_stats = if de_on {
        let pct = |v: &mut Vec<f64>, q: f64| -> f64 { if v.is_empty() { return 0.0; } v.sort_by(|a, b| a.partial_cmp(b).unwrap()); v[((v.len() - 1) as f64 * q).round() as usize] };
        let mut ab: Vec<f64> = de_mm.iter().zip(fixed.iter()).filter(|(_, &f)| f == 0).map(|(&v, _)| (v as f64 / 1000.0).abs()).collect();
        Some(serde_json::json!({"own_solve": de_n.get(&(iy, ix)).map(|q| q.2.clone()), "eta_abs_p50_m": pct(&mut ab, 0.5), "eta_abs_p90_m": pct(&mut ab, 0.9),
            "eta_abs_p99_m": pct(&mut ab, 0.99), "eta_abs_max_m": ab.last().copied().unwrap_or(0.0), "seam_posts": spread_z.len(),
            "seam_spread_z_p99_m": pct(&mut spread_z, 0.99), "seam_spread_z_max_m": spread_z.last().copied().unwrap_or(0.0),
            "seam_spread_eta_p99_m": pct(&mut spread_e, 0.99), "seam_spread_eta_max_m": spread_e.last().copied().unwrap_or(0.0)}))
    } else { None };
    let pyro_stats = if pyon {
        let names = ["lobe_column", "edifice_sequence", "strata"];
        let mut m = serde_json::Map::new();
        for (si, nm) in names.iter().enumerate() {
            let mut q = serde_json::Map::new();
            for (c, &n) in pcount[si].iter().enumerate() { if n > 0 && c != crate::pyro::SUB_NONE as usize { q.insert(c.to_string(), n.into()); } }
            m.insert(nm.to_string(), serde_json::Value::Object(q));
        }
        Some(serde_json::Value::Object(m))
    } else { None };
    let (prov_id, prov_rows) = match &ctx.prov { Some(pv) => prov_ids(ctx, pv, ox, oy, &lay, &pvs), None => (Vec::new(), Vec::new()) };
    TileOut {
        ix, iy,
        height_mm: hm.iter().map(|&v| v as i32).collect(),
        base_mm: pvs.iter().map(|v| to_mm(v.base)).collect(),
        litho, layer, fixed, closure_mm, cut_minus_base_m: cmb,
        base_mean_m: base_mean, raw_t1b_mean_m: raw_t1b_mean, mean_shift_c_m: c, mean_deficit_m: deficit, shift_max_m: shift_max, lip_posts, lip_raise_mm, fill_raised_posts: raised, fill_max_raise_mm: max_raise,
        relax_seconds_own: own_secs, assemble_seconds: t.elapsed().as_secs_f64(),
        bed, over, remap_mm, prov_id, prov_rows, depth_in_unit_m: diu, height_above_unit_base_m: hab, bed_style_id: bsid, flow_unit: fu, pyro_sub: psub, pyro_stats, de_mm, de_talus, de_stats,
        sharp_mm, ct_regime: ct_rg2, sharp_stats,
        levee_mean_m: levee_mean, carve_mean_m: carve_mean, b8_cap_infeasible: cap_infeasible,
    }
}

/// Lane PV: ids of the exposed unit for the tile posts (`lay` from the cut above) and the record of every id present.
fn prov_ids(ctx: &Ctx, pv: &DetailProv, ox: i64, oy: i64, lay: &[usize], pvs: &[PostVal]) -> (Vec<u32>, Vec<serde_json::Value>) {
    use crate::strata_prov::*;
    let r = ctx.p.prov.as_ref().map(|q| q.t1b_search_posts).unwrap_or(2);
    let w = ctx.w;
    // nearest T1b post (Chebyshev rings, then Euclidean within the ring) whose value is not `none`
    let near = |x: f64, y: f64, f: &Option<Vec<u8>>, none: u8| -> Option<u8> {
        let f = f.as_ref()?;
        let (fx, fy) = (x / T1B_M, y / T1B_M);
        let (gx, gy) = (fx.round() as i64, fy.round() as i64);
        let mut best: Option<(f64, u8)> = None;
        for ring in 0..=r {
            for dy in -ring..=ring { for dx in -ring..=ring {
                if dx.abs().max(dy.abs()) != ring { continue; }
                if let Some(i) = w.idx(gx + dx, gy + dy) {
                    let v = f[i];
                    if v == none { continue; }
                    let d = ((gx + dx) as f64 - fx).hypot((gy + dy) as f64 - fy);
                    if best.map(|b| d < b.0).unwrap_or(true) { best = Some((d, v)); }
                }
            } }
            if best.is_some() { break; }
        }
        best.map(|b| b.1)
    };
    let mut kept: BTreeMap<u32, Vec<u16>> = BTreeMap::new();
    let n = lay.len();
    let mut ids = vec![PROV_NONE; n];
    for k in 0..n {
        let l = lay[k];
        if l == usize::MAX || pvs[k].fixed_sea { continue; }
        let (x, y) = ((ox + (k % TILE_POSTS) as i64) as f64 * POST_M, (oy + (k / TILE_POSTS) as i64) as f64 * POST_M);
        let cell = pvs[k].cell;
        ids[k] = match l {
            253 if ctx.p.v2.as_ref().is_some_and(|v| v.edifice_support == Some(true)) => {
                let id = ctx.v2.as_ref().and_then(|v| v.aux.supported_feature_at(w, x, y))
                    .expect("ed_support invariant: edifice class without contributing supported feature");
                PROV_EDIFICE | id as u32
            }
            253 => PROV_EDIFICE | near(x, y, &w.volcano_id, 255).unwrap_or(255) as u32,
            254 => PROV_DEPOSIT | match near(x, y, &w.over, 0) { Some(2) => DEPOSIT_FAN, Some(3) => DEPOSIT_BEACH, _ => DEPOSIT_UNKNOWN },
            255 => PROV_NOCOL | cell,
            _ => {
                let m = kept.entry(cell).or_insert_with(|| crate::strata::kept_full_index(&ctx.t0.off, &ctx.t0.th, cell, ctx.p.min_layer_m));
                match m.get(l) { Some(&fk) => strata_id(cell, fk), None => PROV_NOCOL | cell }
            }
        };
    }
    let set: std::collections::BTreeSet<u32> = ids.iter().copied().filter(|&v| v != PROV_NONE).collect();
    let rows = set.into_iter().map(|id| match id & PROV_KIND_MASK {
        PROV_EDIFICE => {
            let v = id & 0xFF;
            let info = w.volcanoes.get(v as usize);
            let rec = match (&pv.volc, info) { (Some((m, _)), Some(i)) => m.get(&i.id), _ => None };
            edifice_row(v, info, rec, pv.strata.present_t_my)
        }
        PROV_DEPOSIT => deposit_row(id & 0xF),
        PROV_NOCOL => nocol_row(id & 0x0FFF_FFFF, None),
        _ => {
            let (c, k) = (id >> 8, (id & 0xFF) as u16);
            let (a, b) = (ctx.t0.off[c as usize] as usize, ctx.t0.off[c as usize + 1] as usize);
            pv.strata.strata_row(c, k, ctx.t0.lith[a + k as usize], ctx.t0.th[a + k as usize], b - a)
        }
    }).collect();
    (ids, rows)
}
