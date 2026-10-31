//! T1b lane STEP (2026-10-16; `unpublished-design-record`): law threshold-slope-lithology-v1,
//! rule ct_cliff_talus_v1. LABEL: pre-signoff, fidelity tier Preview, derived, never canon; construction (no SimTime).
//! Everything is opt-in under `params.ls.tb.tc.ct` (key `ct=1`); absent = the 0.9.x eikonal byte-identical.
//!
//! Why: the 0.9.x hillslope is an eikonal |grad z| = S from channels, floors and water with S the steady-state nonlinear-creep
//! gradient Sc(rock) x rockmass x f(beta). On hard rock with a large E x L_H, f -> 1 and S -> Sc, so basalt (Sc tan 60 deg)
//! stands as smooth 40-60 deg slopes ("domes"); real hard-rock terrain steps into cliffs over talus.
//!
//! Rule: a debris or soil mantle cannot stand steeper than the repose angle of its debris (talus tan b). Where the creep
//! gradient S <= tan b the slope stays mantled (v0 law unchanged: soft rock and slowly eroding hard rock keep soil slopes).
//! Where S > tan b on a cliff-forming rock the slope cannot be mantled: it is bare rock, split into a CLIFF at the rock-mass
//! angle a and TALUS at repose b, with the cliff height fraction c chosen so the mean gradient equals S
//! (c / tan a + (1 - c) / tan b = 1 / S): hillslope relief (crest heights, the E solve, F6) is conserved, the slope
//! distribution becomes bimodal. Layered rocks step in absolute height (flat-lying units, period P, datum sea level: benches
//! line up across valleys); massive rocks carry one cliff band per hillslope (talus from the hillslope base up to
//! base + (1 - c)(crest - base), cliff above). A basal apron of talus (height (1 - c) P above the base) keeps layered walls
//! from rising straight out of a flat floor. The eikonal step is integrated exactly across band boundaries (no aliasing of a
//! 128 m period on 32 m posts).

use crate::*;
use serde::{Deserialize, Serialize};

/// Structure modes of a class.
pub const CT_NONE: u32 = 0;
pub const CT_MASSIVE: u32 = 1;
pub const CT_LAYERED: u32 = 2;

/// Published regime codes (`ct_regime` window field).
pub const RG_SOIL: u8 = 0;
pub const RG_TALUS: u8 = 1;
pub const RG_CLIFF: u8 = 2;
/// Lane TREAD2 (ct_tread only): structural bench on the top of a hard bed.
pub const RG_TREAD: u8 = 3;
/// Sentinel of `ct_top_mm` where there is no massive talus top.
pub const TOP_NONE: i32 = i32::MIN;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CtClassV1 {
    pub code: u8,
    /// CT_NONE, CT_MASSIVE, CT_LAYERED.
    pub mode: u32,
    /// Cliff (rock-mass) angle and talus repose angle as tangents, permille.
    pub cliff_tan_permille: i64,
    pub talus_tan_permille: i64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CtParamsV1 {
    /// ESCARP-2 per-post, ESCARP-3 band, ESCARP-4 anti-aliased band, ESCARP-6 small toes.
    /// Absent preserves the historical byte stream and ESCARP-1. No extra erosion.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub escarp_reconstruct: Option<u32>,
    /// 1 = rule on.
    pub on: u32,
    /// 0 = the ct gradient drives the pass-2 eikonal itself (crest_exempt applies); 1 = superposition: the v0 pass-2 surface
    /// is crest-rounded as v0 and the ct staircase enters as the difference of a second (ct) pass-2 eikonal from the v0 one.
    pub apply: u32,
    /// MEASUREMENT switch (key ct_window=1): ct acts in rnwin windows only; the layer (E table, network) is the one solved
    /// without ct, so a window ablation keeps the same network and E (rnwin checks the layer against these params minus ct).
    /// A layer is never built with it.
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub window_only: u32,
    /// apply 1: |z_ct - z_v0| is capped at this fraction of the step period (permille; 0 = no cap). The staircase is a perturbation of one step period; larger
    /// differences come from the two eikonals choosing different paths (rock read at different heights), not from the steps.
    pub shift_cap_permille: i64,
    /// Layered step period (mm) and datum (mm, absolute height of a step base).
    pub step_mm: i64,
    pub datum_mm: i64,
    /// Bare-rock posts are exempt from crest rounding (flag).
    pub crest_exempt: u32,
    /// Basal talus apron on layered classes (flag).
    pub apron: u32,
    /// The rockmass multiplier also scales the cliff tangent (flag).
    pub rockmass_cliff: u32,
    pub classes: Vec<CtClassV1>,
    /// Lane STEP2 options (absent = lane STEP as committed; byte-identical).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub v2: Option<CtV2ParamsV1>,
    pub law: String,
    pub regime_rule: String,
    pub geometry_rule: String,
    pub crest_rule: String,
    pub table_source: String,
}

fn tanp(deg: f64) -> i64 { (deg.to_radians().tan() * 1000.0).round() as i64 }

impl CtParamsV1 {
    pub fn v1() -> Self {
        let c = |code: u8, mode: u32, cliff_deg: f64| CtClassV1 { code, mode, cliff_tan_permille: if mode == CT_NONE { 0 } else { tanp(cliff_deg) }, talus_tan_permille: if mode == CT_NONE { 0 } else { tanp(35.0) } };
        CtParamsV1 {
            escarp_reconstruct: None,
            on: 1, apply: 1, window_only: 0, shift_cap_permille: 1_000, step_mm: 128_000, datum_mm: 0, crest_exempt: 1, apron: 1, rockmass_cliff: 1,
            v2: None,
            classes: vec![
                c(0, CT_MASSIVE, 70.0), c(1, CT_MASSIVE, 75.0), c(2, CT_MASSIVE, 60.0), c(3, CT_LAYERED, 80.0),
                c(4, CT_LAYERED, 80.0), c(5, CT_NONE, 0.0), c(6, CT_LAYERED, 75.0), c(7, CT_LAYERED, 65.0),
                c(8, CT_NONE, 0.0), c(9, CT_NONE, 0.0), c(10, CT_NONE, 0.0), c(11, CT_MASSIVE, 65.0),
            ],
            law: "threshold-slope-lithology-v1 (lane STEP): v0's eikonal gradient S = Sc(rock at the upwind height) x rockmass x f(beta) where S <= the class talus tangent or the class has no cliff mode; elsewhere ct_cliff_talus_v1".into(),
            regime_rule: "ct_cliff_talus_v1 regime: a debris/soil mantle cannot stand steeper than its repose angle; where the steady-state creep gradient S exceeds the class talus tangent tb on a cliff-forming class the slope is bare rock (weathering-limited) and splits into cliff at the class cliff tangent ta (x rockmass when rockmass_cliff) and talus at tb, cliff height fraction c = (1/tb - 1/S) / (1/tb - 1/ta) clamped [0, 1] so c/ta + (1-c)/tb = 1/S (the mean gradient and the hillslope relief are those of v0); S >= ta: all cliff at S".into(),
            geometry_rule: "layered classes: in absolute height, period P = step_mm from datum_mm, the lower (1-c)P of each period talus/bench, the upper cP cliff/riser; apron: the first (1-c)P above the hillslope base (height of the channel, floor or water post reached by D8 steepest descent on the pass-1 surface) is talus. massive classes: talus from the base up to base + (1-c)(crest - base) (crest: the highest pass-1 post draining through the post by D8), cliff above. The eikonal speed at post i from the upwind height zu is the exact 1-D climb over one post spacing across band boundaries ((z(dx) - zu) / dx); class, S and c are those at zu (a lithology contact inside the step is not split)".into(),
            crest_rule: "apply 1 (superposition, default): the v0 pass-2 eikonal (no ct) is computed and crest-rounded exactly as v0; a second pass-2 eikonal with the ct gradient gives z_ct; every free non-floor post then takes z_round_v0 + clamp(z_ct - z_v0, +-shift_cap_permille/1000 x P) before the bank rule and fill (the staircase rides on v0's rounded hill shape; regime is read at z_ct). apply 0: the ct gradient drives the only pass-2 eikonal; crest_exempt: posts whose regime at their eikonal height is talus or cliff skip the min(z, box mean) crest rounding (the rounding radius is the soil-creep curvature E/D; bare rock has no mantle to creep)".into(),
            table_source: "project assumptions (lane STEP, recalled, unverified, to calibrate): talus repose 35 deg on every cliff-forming class (talus slopes 32-38 deg, Statham 1976 Earth Surf. Proc. 1:173); cliff angles from rock-mass strength classes (Selby 1980 Z. Geomorph. 24:31; Moon 1984): basalt 80 (columnar/blocky jointed flows), carbonate 80, sandstone 75, conglomerate 65, granite 75, crystalline basement 70, metamorphic 60, ophiolite 65; shale, evaporite, coal, unconsolidated never cliff here; layered = flat-lying units (basalt flow packages, bedded sediments), step period 128 m a lattice-limited closure (lane STEP measured 64 m as unresolved on 32 m posts: a 64 m period is one central-difference span; 128 and 192 m resolve; real flow units 10-50 m belong to the 2 m detail); bare-rock transition cf. Heimsath et al. 2012 Nature Geosci. 5:210, DiBiase et al. 2012 EPSL 289:134 (recalled)".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.escarp_reconstruct.is_some_and(|v| v != 2 && v != 3 && v != 4 && v != 6) { return Err("ct: unsupported escarp reconstruction mode".into()); }
        if self.escarp_reconstruct.is_some() && !self.v2.as_ref().is_some_and(|v|v.base==1) { return Err("ct_escarp=2|3|4|6 needs ct_base=1 to carry the native detail parameters".into()); }
        if self.on != 1 || self.apply > 2 || self.window_only > 1 { return Err("ct: on must be 1 and apply 0, 1 or 2".into()); }
        if self.apply == 2 && self.v2.is_none() { return Err("ct: apply 2 needs the lane STEP2 params (v2)".into()); }
        if let Some(v) = &self.v2 {
            if v.var > 3 || v.base > 1 || v.local > 1 { return Err("ct: var 0..3, base and local 0 or 1".into()); }
            if (v.var == 3) != v.bed.is_some() { return Err("ct: var 3 needs the lane BENCH-VAR bed params and they need var 3".into()); }
            if let Some(b) = &v.bed { b.validate()?; }
            if v.local == 1 && self.apply != 2 { return Err("ct: local needs apply 2".into()); }
            if v.ends > 1 || (v.ends == 1 && (self.apply != 2 || v.ends_rule.is_empty())) { return Err("ct: ends 0 or 1; ends needs apply 2 and its rule".into()); }
            if v.var >= 1 && self.apply != 2 { return Err("ct: var (variable step thickness) needs apply 2".into()); }
            if v.base == 1 && self.apply == 0 { return Err("ct: base (ct_base_mm) needs apply 1 or 2".into()); }
            if let Some(t) = &v.tread {
                if self.apply != 2 { return Err("ct: tread needs apply 2".into()); }
                t.validate()?;
                if self.classes.iter().any(|k| k.mode != CT_NONE && k.talus_tan_permille <= t.tread_max_tan_permille) { return Err("ct_tread: the tread tangent must stay below every class talus tangent".into()); }
            }
            if !(v.var_lo_permille > 0 && v.var_hi_permille > v.var_lo_permille && v.unit_min_permille > 0 && v.unit_max_permille >= v.unit_min_permille) { return Err("ct: 0 < var_lo < var_hi, 0 < unit_min <= unit_max".into()); }
        }
        if self.step_mm <= 0 || self.shift_cap_permille < 0 { return Err("ct: step_mm > 0".into()); }
        if self.classes.len() != crate::t1b_litho::LITHO_CODES { return Err("ct: classes must list all 12 litho codes".into()); }
        for (i, c) in self.classes.iter().enumerate() {
            if c.code as usize != i { return Err("ct: classes out of code order".into()); }
            if c.mode > CT_LAYERED { return Err("ct: mode 0..2".into()); }
            if c.mode != CT_NONE && !(c.talus_tan_permille > 0 && c.cliff_tan_permille > c.talus_tan_permille) { return Err("ct: 0 < talus < cliff tangent".into()); }
        }
        Ok(())
    }
}

/// Lane STEP2 (2026-10-16, opt-in sub-keys of ct; `unpublished-design-record`).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CtV2ParamsV1 {
    /// Step thickness source (apply 2 only): 0 the regular period P from datum (lane STEP); 1 keyed closure; 2 strata
    /// packages of the layered classes where the record has them, else the keyed closure.
    pub var: u32,
    /// Keyed closure: period thickness P x U(var_lo, var_hi) (permille).
    pub var_lo_permille: i64,
    pub var_hi_permille: i64,
    /// Strata: a layered package up to unit_max x P thick is one period; thicker packages are subdivided at the closure
    /// boundaries (pieces >= unit_min x P); consecutive layered units merge into packages of at least unit_min x P.
    pub unit_min_permille: i64,
    pub unit_max_permille: i64,
    /// 1 = publish ct_base_mm (the surface without the staircase shift) in windows.
    pub base: u32,
    /// apply 2: 1 = ct_remap_local_v1 (regime and cliff fraction from the local gradient of the rounded v0 surface).
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub local: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub local_rule: String,
    /// Lane TREAD2 (key ct_tread=1; absent = lane STEP2 byte-identical): the 3-limb tread/talus/cliff profile.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tread: Option<CtTreadParamsV1>,
    /// Lane BENCH-VAR (key ct_var=3; absent = lane TREAD2 / STEP2 byte-identical): per-bed keyed thickness and hard share.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bed: Option<CtBedParamsV1>,
    /// Lane PERIOD-ENDS (key ct_ends=1, apply 2; absent = byte-identical): publish each remapped post's period ends
    /// (ct_s0_mm, ct_s1_mm), the neighbouring period ends (ct_s_dn_mm, ct_s_up_mm), the remap gradient (ct_grad_permille),
    /// the period kind (ct_period_kind) and the var 3 bed identity (ct_bed_id) in windows. Publication only: no height changes.
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub ends: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub ends_rule: String,
    pub remap_rule: String,
    pub var_rule: String,
    pub base_rule: String,
}

pub const LOCAL_RULE: &str = "ct_remap_local_v1 (apply 2, lane STEP2): the bare-rock test and the split use the local gradient G of the crest-rounded v0 surface (central differences over the 32 m posts, one-sided where a neighbour is missing) in place of the steady creep gradient S: G <= the class talus tangent stays mantled (unchanged: convex soil crests and gentle slopes), else c from c/ta + (1-c)/tb = 1/G and F_c uses G (along the local slope the remap gives exactly tb then ta); c -> 0 as G -> tb, so F is continuous across the mantle/bare-rock boundary; G >= ta unchanged";

pub const PURPOSE_CT_STEP: &str = "t1b.rivernet.ct.step";
/// The keyed closure table covers datum +- this height (m); beyond it the regular period P continues (declared).
pub const CLOSURE_SPAN_M: f64 = 12_000.0;

impl CtV2ParamsV1 {
    pub fn v1() -> Self {
        CtV2ParamsV1 {
            var: 0, var_lo_permille: 500, var_hi_permille: 1_500, unit_min_permille: 500, unit_max_permille: 1_500, base: 0, local: 0, local_rule: String::new(), tread: None, bed: None, ends: 0, ends_rule: String::new(),
            remap_rule: "ct_remap_v1 (apply 2, lane STEP2): no second eikonal. The v0 pass-2 surface is crest-rounded exactly as v0; every free non-floor post with rounded v0 height h whose rock at h is bare (S > talus tangent on a cliff-forming class; S, class, c read at h as in ct_cliff_talus_v1) takes F_c(h) = s0 + T g(u), u = (h - s0) / T in its period [s0, s1] (T = s1 - s0), g(u) = u tb / S for u < u_t = (1 - c) S / tb (talus), else (1 - c) + (u - u_t) ta / S (cliff): along a v0 slope of gradient S the remapped slope is exactly tb then ta, F is monotone, F(s0) = s0 and F(s1) = s1 (period boundaries and every post's position fixed: no path divergence, no shift cap; the mean gradient over each period is S); F <= h inside a period (benches sit below the v0 line, as the lane STEP climb). Apron: the period holding the post's hillslope base (pass-1 D8 source height) starts at the base (F(base) = base). Massive classes: one period [base, crest] (pass-1), else the layered rule. S >= ta: unchanged (all cliff)".into(),
            var_rule: "ct_step_var_v1 (apply 2): var 0 regular P from datum. var 1 keyed closure: global boundaries in absolute height from datum (b_0 = datum), period thickness P x U(var_lo, var_hi), u keyed by purpose t1b.rivernet.ct.step on the chart-0 scope of the draw digest (counter 2k for the k-th period up, 2k + 1 for the k-th down, attempt = realization; the same boundaries in every window and valley), table +- 12 km, regular P beyond. var 2 strata: the post's T0 column (the lithology cut's units, datum and tilt) is grouped into packages of consecutive layered-class units (greedy from the top, each >= unit_min x P, a short remainder joins the package above; non-layered units break packages; the top unit is open upward, the bottom one downward); a post reading a layered class takes its package [lo, hi] as one period when hi - lo <= unit_max x P, else the closure boundaries inside it at least unit_min x P from its ends; the stand-in lithology and posts outside packages use the closure".into(),
            base_rule: "ct_base_v1: ct_base_mm = the published height minus the staircase shift of the post (apply 1: the clamped z_ct - z_v0; apply 2: F_c(h) - h) on free posts not rebuilt by the coast/volcano stage; = height elsewhere (fill and bank raises after the shift are kept in both)".into(),
        }
    }
}

/// The keyed closure boundaries (ascending absolute heights, m) for datum +- CLOSURE_SPAN_M; `u(counter)` the keyed uniform.
pub fn closure_table(c: &CtParamsV1, u: &mut dyn FnMut(u64) -> R<f64>) -> R<Vec<f64>> {
    let v = c.v2.as_ref().ok_or("ct: closure needs v2")?;
    let p = c.step_mm as f64 / 1000.0;
    let (lo, hi) = (v.var_lo_permille as f64 / 1000.0, v.var_hi_permille as f64 / 1000.0);
    let d = c.datum_mm as f64 / 1000.0;
    let (mut up, mut dn) = (vec![d], Vec::new());
    let (mut z, mut k) = (d, 0u64);
    while z < d + CLOSURE_SPAN_M { z += p * (lo + (hi - lo) * u(2 * k)?); up.push(z); k += 1; }
    let (mut z, mut k) = (d, 0u64);
    while z > d - CLOSURE_SPAN_M { z -= p * (lo + (hi - lo) * u(2 * k + 1)?); dn.push(z); k += 1; }
    dn.reverse();
    dn.extend(up);
    Ok(dn)
}

/// Last boundary <= x and first boundary > x of the closure (regular P from datum when the table is empty or beyond it).
fn below(c: &CtParamsV1, tab: &[f64], x: f64) -> f64 {
    let p = c.step_mm as f64 / 1000.0;
    if !x.is_finite() { return x; }
    if tab.is_empty() { let d = c.datum_mm as f64 / 1000.0; return d + ((x - d) / p).floor() * p; }
    let j = tab.partition_point(|b| *b <= x);
    if j == 0 { return tab[0] - ((tab[0] - x) / p).ceil() * p; }
    if j == tab.len() { let t = tab[j - 1]; return t + ((x - t) / p).floor() * p; }
    tab[j - 1]
}
fn above(c: &CtParamsV1, tab: &[f64], x: f64) -> f64 {
    let p = c.step_mm as f64 / 1000.0;
    if !x.is_finite() { return x; }
    if tab.is_empty() { let d = c.datum_mm as f64 / 1000.0; return d + ((x - d) / p).floor() * p + p; }
    let j = tab.partition_point(|b| *b <= x);
    if j == 0 { return tab[0] - ((tab[0] - x) / p).ceil() * p + p; }
    if j == tab.len() { let t = tab[j - 1]; return t + ((x - t) / p).floor() * p + p; }
    tab[j]
}

/// The period [s0, s1) holding height z inside the package [lo, hi] (+-inf: open; the closure / regular boundaries).
pub fn period(c: &CtParamsV1, tab: &[f64], z: f64, lo: f64, hi: f64) -> (f64, f64) {
    let p = c.step_mm as f64 / 1000.0;
    let (umin, umax) = c.v2.as_ref().map(|v| (v.unit_min_permille as f64 / 1000.0 * p, v.unit_max_permille as f64 / 1000.0 * p)).unwrap_or((0.0, 0.0));
    let var3 = c.v2.as_ref().map(|v| v.var == 3).unwrap_or(false);
    if lo.is_finite() && hi.is_finite() && (hi - lo <= umax || var3) { return (lo, hi); }
    let mut s0 = below(c, tab, z);
    if s0 > hi - umin { s0 = below(c, tab, hi - umin); }
    if s0 < lo + umin { s0 = lo; }
    let mut s1 = above(c, tab, z);
    if s1 < lo + umin { s1 = above(c, tab, lo + umin - 1e-9); if s1 < lo + umin { s1 = lo + umin; } }
    if s1 > hi - umin { s1 = hi; }
    (s0, s1)
}

/// Packages of consecutive layered-class units of a column (units top-down: depth bases below the datum, codes), as depth
/// intervals (top, bottom); the top unit is open upward (-inf), the bottom one downward (+inf).
pub fn packages(c: &CtParamsV1, base: &[f64], code: &[u8]) -> Vec<(f64, f64)> {
    let p = c.step_mm as f64 / 1000.0;
    let umin = c.v2.as_ref().map(|v| v.unit_min_permille as f64 / 1000.0 * p).unwrap_or(0.0);
    let n = code.len();
    let mut out: Vec<(f64, f64)> = Vec::new();
    let mut run_start = 0usize; // index in `out` where the current run's packages begin
    let mut cur: Option<(f64, f64)> = None;
    for j in 0..n {
        let top = if j == 0 { f64::NEG_INFINITY } else { base[j - 1] };
        let bot = if j + 1 == n { f64::INFINITY } else { base[j] };
        let layered = c.classes.get(code[j] as usize).map(|k| k.mode == CT_LAYERED).unwrap_or(false);
        if !layered {
            if let Some(q) = cur.take() { if out.len() > run_start { out.last_mut().unwrap().1 = q.1; } else { out.push(q); } }
            run_start = out.len();
            continue;
        }
        let q = match cur { Some((t, _)) => (t, bot), None => (top, bot) };
        if q.1 - q.0 >= umin { out.push(q); cur = None; } else { cur = Some(q); }
    }
    if let Some(q) = cur.take() { if out.len() > run_start { out.last_mut().unwrap().1 = q.1; } else { out.push(q); } }
    out
}

/// Lane STEP2 apply 2: the remap of one post. (new height, regime). `lo, hi`: the post's package (heights; +-inf open).
pub fn remap(c: &CtParamsV1, tab: &[f64], mode: u32, tb: f64, ta: f64, cf: f64, s: f64, g: Geo, h: f64, lo: f64, hi: f64) -> (f64, u8) {
    if cf >= 1.0 { return (h, RG_CLIFF); }
    let ut = ((1.0 - cf) * s / tb).min(1.0);
    let gf = |u: f64| -> (f64, u8) { if u < ut { (u * tb / s, RG_TALUS) } else { (((1.0 - cf) + (u - ut) * ta / s).min(1.0), RG_CLIFF) } };
    if mode == CT_MASSIVE && g.base.is_finite() && g.crest.is_finite() && g.crest > g.base {
        if h <= g.base { return (h, RG_TALUS); }
        if h >= g.crest { return (h, RG_CLIFF); }
        let t = g.crest - g.base;
        let (v, rg) = gf((h - g.base) / t);
        return (g.base + t * v, rg);
    }
    let (mut s0, s1) = period(c, tab, h, lo, hi);
    if c.apron == 1 && g.base.is_finite() && g.base > s0 && g.base < s1 && h >= g.base { s0 = g.base; }
    let t = s1 - s0;
    if !(t > 0.0) { return (h, RG_TALUS); }
    let (v, rg) = gf(((h - s0) / t).clamp(0.0, 1.0));
    (s0 + t * v, rg)
}

/// Lane PERIOD-ENDS (2026-10-17, opt-in sub-key ct_ends=1; `unpublished-design-record`).
pub const ENDS_RULE: &str = "ct_ends_v1 (lane PERIOD-ENDS, apply 2, publication only): for every post remapped by ct_remap_v1 / ct_tread_v1 (2-limb, 3-limb or fallback) the period [s0, s1] its own remap used (remap_period: massive classes [hillslope base, crest]; layered the package / keyed bed / closure period of `period`, s0 = the hillslope base where the apron applies; none where c >= 1 or the remap was the identity), the neighbouring period ends s_dn and s_up (the bottom of the period holding s0 - 1 um and the top of the period holding s1 + 1 um by the same rule, the same post's column, class, package and hillslope base read at that height; none where that height reads a non-layered class or the own period is massive), the knots of F in absolute mm (tread -> talus and talus -> cliff breaks, (h, F(h)) each: F is exact to the mm quantization, which the permille limb fields are not), the remap gradient S (the local or steady gradient the split used; cliff tangent ta = S x the cliff limb's F-slope), the period kind and, under var 3, the bed identity ((package id + 1) << 16 | bed index; one-period record packages index 0). Posts whose height is not F(h) on the published knots are flagged (kind | 128): held by the tread base-level or drainage rules, or moved by a later stage of the surface (bank rule, fill, coast) after the remap";
/// `ct_period_kind` codes (lane PERIOD-ENDS).
pub const PK_NONE: u8 = 0;
/// Closure / regular period (no package).
pub const PK_CLOSURE: u8 = 1;
/// A package taken whole as one period (var 2 <= unit_max, var 3 record package).
pub const PK_PACKAGE: u8 = 2;
/// Var 3 keyed bed inside a thick package.
pub const PK_BED: u8 = 3;
/// Var 2 closure period inside a package thicker than unit_max (clipped at the package ends).
pub const PK_SUBDIV: u8 = 4;
/// Apron period: s0 = the hillslope base.
pub const PK_APRON: u8 = 5;
/// Massive class: [hillslope base, crest].
pub const PK_MASSIVE: u8 = 6;
/// Flag: height != F(h) on the published knots (tread base-level / drainage hold, or moved by a later surface stage).
pub const PK_HELD: u8 = 0x80;
/// Sentinel of the ct_s*_mm fields where there is no period.
pub const ENDS_NONE: i32 = i32::MIN;
pub const ENDS_NOTE: &str = "lane PERIOD-ENDS (opt-in ct_ends, rule params.ls.tb.tc.ct.v2.ends_rule): absolute height (mm) of a period end the post's own remap used: ct_s0_mm / ct_s1_mm the bottom / top of the post's period [s0, s1] (F(s0) = s0, F(s1) = s1), ct_s_dn_mm the bottom of the period below s0, ct_s_up_mm the top of the period above s1 (the same post's remap rule read at those heights); -2147483648 (i32::MIN) where none (not remapped: soil, fixed, coast/volcano rebuilt, c >= 1, identity; neighbour periods: non-layered class there or a massive own period)";
pub const KNOT_NOTE: &str = "lane PERIOD-ENDS (opt-in ct_ends): the knots of the post's remap F over its period [s0, s1] in absolute mm: F climbs (s0, s0) -> (ct_knot_h1_mm, ct_knot_z1_mm) -> (ct_knot_h2_mm, ct_knot_z2_mm) -> (s1, s1), linear between (h the pre-staircase height ct_base_mm, z = F(h)); knot 1 = tread -> talus (= s0 where there is no tread), knot 2 = talus -> cliff; the same fractions of the period apply in the neighbouring periods; i32::MIN where no period";
pub const GRAD_NOTE: &str = "lane PERIOD-ENDS (opt-in ct_ends): the gradient S (tangent, permille) the post's remap split used (ct_remap_local_v1: the local gradient of the crest-rounded surface; else the steady creep gradient); with the class talus tangent, ct_cliff_permille, ct_tread_permille and ct_tread_tan_permille the limb breakpoints are exact (rule params.ls.tb.tc.ct.v2.ends_rule); 0 where no period";
pub const KIND_NOTE: &str = "lane PERIOD-ENDS (opt-in ct_ends): what the period is: 0 none, 1 closure/regular, 2 whole package, 3 var 3 keyed bed, 4 var 2 closure period inside a thick package, 5 apron (s0 = hillslope base), 6 massive [base, crest]; + 128 where height_mm is not F(ct_base_mm) on the published knots (held by the tread base-level or drainage rule, or moved after the remap by the bank rule, fill or coast stage)";
pub const BED_ID_NOTE: &str = "lane PERIOD-ENDS (opt-in ct_ends, var 3 only): identity of the package / bed holding the post's height: (package id + 1) << 16 | bed index (package id = index from the bottom of the T0 full column of its lowest unit; one-period record packages index 0); the same bed has the same id in every column; 0 none";

/// Lane PERIOD-ENDS: the knots of F over one period as fractions [u1, f1, u2, f2]: F climbs (0, 0) -> (u1, f1) -> (u2, f2) ->
/// (1, 1) in ((h - s0) / T, (F - s0) / T), linear between (the expressions of `remap` / `remap3_raw`). 2-limb (`remap`,
/// cliff fraction cf, gradient s, talus tangent tb): no tread (u1 = f1 = 0), u2 = min((1 - cf) s / tb, 1), f2 = 1 - cf.
pub fn knots2(cf: f64, tb: f64, s: f64) -> [f64; 4] { [0.0, 0.0, ((1.0 - cf) * s / tb).min(1.0), 1.0 - cf] }
/// 3-limb (`remap3_raw`, fractions ft, fb, tread tangent tt): u1 = ft s / tt, f1 = ft, u2 = u1 + fb s / tb, f2 = ft + fb.
pub fn knots3(ft: f64, fb: f64, tt: f64, tb: f64, s: f64) -> [f64; 4] { let ut = ft * s / tt; [ut, ft, ut + fb * s / tb, ft + fb] }

/// Lane PERIOD-ENDS: the kind of the period `period` gives for package [lo, hi].
pub fn period_kind(c: &CtParamsV1, lo: f64, hi: f64) -> u8 {
    let p = c.step_mm as f64 / 1000.0;
    let umax = c.v2.as_ref().map(|v| v.unit_max_permille as f64 / 1000.0 * p).unwrap_or(0.0);
    let var3 = c.v2.as_ref().map(|v| v.var == 3).unwrap_or(false);
    if lo.is_finite() && hi.is_finite() && (hi - lo <= umax || var3) { PK_PACKAGE } else if lo.is_finite() || hi.is_finite() { PK_SUBDIV } else { PK_CLOSURE }
}

/// Lane PERIOD-ENDS: the period [s0, s1] and kind `remap` (and remap3, which takes the same period) uses at height h, or None
/// where the remap is the identity or has no period (c >= 1, massive outside (base, crest), empty period).
pub fn remap_period(c: &CtParamsV1, tab: &[f64], mode: u32, cf: f64, g: Geo, h: f64, lo: f64, hi: f64) -> Option<(f64, f64, u8)> {
    if cf >= 1.0 { return None; }
    if mode == CT_MASSIVE && g.base.is_finite() && g.crest.is_finite() && g.crest > g.base {
        if h <= g.base || h >= g.crest { return None; }
        return Some((g.base, g.crest, PK_MASSIVE));
    }
    let (mut s0, s1) = period(c, tab, h, lo, hi);
    let mut k = period_kind(c, lo, hi);
    if c.apron == 1 && g.base.is_finite() && g.base > s0 && g.base < s1 && h >= g.base { s0 = g.base; k = PK_APRON; }
    if !(s1 - s0 > 0.0) { return None; }
    Some((s0, s1, k))
}

/// Lane TREAD2 (2026-10-17, opt-in sub-key ct_tread=1; `unpublished-design-record`).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CtTreadParamsV1 {
    pub on: u32,
    /// Tread tangent on a horizontal hard-bed top (permille), its clamp after the dip term, and the dip switch (1 = add the
    /// record bed gradient projected on the uphill direction).
    pub tread_tan_permille: i64,
    pub tread_min_tan_permille: i64,
    pub tread_max_tan_permille: i64,
    pub dip: u32,
    /// Hard-bed (cliff) share of a period: closure eta0, cap eta_max (permille).
    pub eta_permille: i64,
    pub eta_max_permille: i64,
    /// Transition: the full 3-limb profile from this gradient up (tangent, permille); identity at the talus tangent.
    pub full_tan_permille: i64,
    pub rule: String,
    pub eta_rule: String,
    pub angle_rule: String,
    pub base_rule: String,
    pub drain_rule: String,
    pub table_source: String,
}

/// Lane BENCH-VAR (2026-10-17, opt-in ct_var=3; `unpublished-design-record`): the step period
/// follows a per-bed keyed thickness (log-normal, clamped) and a per-bed hard share instead of the uniform closure
/// P x U(0.5, 1.5) and the one eta closure.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CtBedParamsV1 {
    /// Bed thickness t = P x median x exp(sigma z), z standard normal (keyed), clamped to P x [min, max] (permille).
    pub median_permille: i64,
    pub sigma_permille: i64,
    pub min_permille: i64,
    pub max_permille: i64,
    /// Per-bed hard (cliff-forming) share where the record has one class through the period: U(eta_lo, eta_hi) (permille).
    pub eta_lo_permille: i64,
    pub eta_hi_permille: i64,
    pub purpose: String,
    pub rule: String,
    pub table_source: String,
}

pub const PURPOSE_CT_BED: &str = "t1b.rivernet.ct.bed";
/// Package identities with a precomputed bed table (index from the column bottom) and beds per table.
pub const BED_IDS: u64 = 64;
pub const BED_PER_ID: u64 = 256;
/// Keyed counters: (id << 24) | (k << 2) | what (0 thickness, 1 hard share); the global up / down tables use these ids.
pub const BED_ID_UP: u64 = 0xFFFE;
pub const BED_ID_DN: u64 = 0xFFFF;
pub fn bed_key(id: u64, k: u64, what: u64) -> u64 { (id << 24) | (k << 2) | what }

pub const BED_RULE: &str = "ct_bed_v1 (lane BENCH-VAR, ct_var=3, apply 2, layered classes): var 2 strata packages as before (a record package <= unit_max x P is one period, its hard share from the record where the T0 column has mixed classes in it, else the keyed share of its identity's bed 0); a package thicker than unit_max x P is cut into BEDS anchored at the package top (depth below the column datum, so beds follow the datum tilt; a package open upward is anchored at its bottom and cut upward): bed k of package identity j has thickness t = P x median x exp(sigma z), z = Phi^-1(u) (Acklam), clamped to P x [min, max], and hard share eta = eta_lo + (eta_hi - eta_lo) u'; u, u' keyed by purpose t1b.rivernet.ct.bed on the chart-0 scope of the draw digest, counter (j << 24) | (k << 2) | {0, 1}, attempt = realization (declared here, not in seed_scope.purposes); j = the index from the column bottom (in the T0 full column) of the package's lowest unit, so the same package has the same bed sequence in every post, window and T0 cell that holds it (same bed = same step; packages do not share a phase); a last bed thinner than P x min joins the bed above (below, upward); j >= 64 or k >= 256: median thickness, eta midpoint (declared). Posts outside packages (stand-in lithology, a package open both ways) use the global table: boundaries from datum up and down with the same thickness law (ids 0xFFFE up, 0xFFFF down; +- 12 km, regular P beyond). Every bed is one period of the remap (2-limb or ct_tread 3-limb); with ct_tread the bed's eta replaces the closure eta0 (capped at eta_max) where the record gives none";

fn norm_inv(p: f64) -> f64 {
    // Acklam's rational approximation of the standard normal quantile (relative error < 1.2e-9).
    const A: [f64; 6] = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const B: [f64; 5] = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
    const C: [f64; 6] = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const D: [f64; 4] = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
    let p = p.clamp(1e-12, 1.0 - 1e-12);
    let pl = 0.02425;
    if p < pl {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= 1.0 - pl {
        let q = p - 0.5; let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

impl CtBedParamsV1 {
    pub fn v1() -> Self {
        CtBedParamsV1 {
            median_permille: 800, sigma_permille: 500, min_permille: 400, max_permille: 1_800, eta_lo_permille: 500, eta_hi_permille: 900,
            purpose: PURPOSE_CT_BED.into(), rule: BED_RULE.into(),
            table_source: "project assumptions (lane BENCH-VAR, recalled, unverified, to calibrate): bed and lava-flow thicknesses are right-skewed (log-normal-like; flood-basalt sheet flows ~5-50 m, rarely > 100 m, compound pahoehoe packets of thin lobes; sedimentary bed sets log-normal or exponential, CV ~0.5-1); a 32 m post lattice resolves packages, not single flows, so a bed here is a cliff-forming flow or flow packet with its slope-forming top: median 0.8 P (102 m at P = 128 m), sigma 0.5, clamp 0.4-1.8 P (51-230 m: lane STEP measured 64 m regular periods unresolved at 32 m): clamped CV ~0.45 (U(0.5, 1.5): 0.29), mean ~0.89 P, post-weighted mean E[T^2]/E[T] ~1.07 P (U(0.5, 1.5): 1.08 P) so the bench lowering (proportional to the period a post sits in) does not grow; lane BENCH-VAR measured median 0.9 / sigma 0.55 / max 2.0 (E[T^2]/E[T] 1.24 P) lowering bare rock ~1 m more; the massive-interior (colonnade) share of a flow varies flow to flow, recalled ~0.3-0.8: U(0.5, 0.9) around the TREAD2 closure 0.7".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        if !(self.min_permille > 0 && self.median_permille >= self.min_permille && self.max_permille >= self.median_permille && self.sigma_permille >= 0) { return Err("ct bed: 0 < min <= median <= max, sigma >= 0".into()); }
        if !(self.eta_lo_permille > 0 && self.eta_hi_permille >= self.eta_lo_permille && self.eta_hi_permille < 1_000) { return Err("ct bed: 0 < eta_lo <= eta_hi < 1000".into()); }
        Ok(())
    }
    /// Bed thickness (m) for the keyed uniform u, period P (m).
    pub fn thickness(&self, p: f64, u: f64) -> f64 {
        let (md, sg) = (self.median_permille as f64 / 1000.0, self.sigma_permille as f64 / 1000.0);
        (p * md * (sg * norm_inv(u)).exp()).clamp(p * self.min_permille as f64 / 1000.0, p * self.max_permille as f64 / 1000.0)
    }
    pub fn eta(&self, u: f64) -> f64 { let (a, b) = (self.eta_lo_permille as f64 / 1000.0, self.eta_hi_permille as f64 / 1000.0); a + (b - a) * u }
    pub fn eta_mid(&self) -> f64 { self.eta(0.5) }
}

/// Lane BENCH-VAR tables: per package identity (0..BED_IDS) its beds (thickness m, eta); the global table (ascending
/// boundaries, m) and the eta of each interval [tab[j], tab[j + 1]].
#[derive(Clone, Debug, Default)]
pub struct BedTables { pub pk: Vec<Vec<(f64, f64)>>, pub glob: Vec<f64>, pub glob_eta: Vec<f64> }

pub fn bed_tables(c: &CtParamsV1, u: &mut dyn FnMut(u64) -> R<f64>) -> R<BedTables> {
    let b = c.v2.as_ref().and_then(|v| v.bed.as_ref()).ok_or("ct: bed tables need var 3")?;
    let p = c.step_mm as f64 / 1000.0;
    let mut pk = Vec::with_capacity(BED_IDS as usize);
    for j in 0..BED_IDS {
        let mut v = Vec::with_capacity(BED_PER_ID as usize);
        for k in 0..BED_PER_ID { v.push((b.thickness(p, u(bed_key(j, k, 0))?), b.eta(u(bed_key(j, k, 1))?))); }
        pk.push(v);
    }
    let d = c.datum_mm as f64 / 1000.0;
    let (mut up, mut ue) = (vec![d], Vec::new());
    let (mut z, mut k) = (d, 0u64);
    while z < d + CLOSURE_SPAN_M { z += b.thickness(p, u(bed_key(BED_ID_UP, k, 0))?); up.push(z); ue.push(b.eta(u(bed_key(BED_ID_UP, k, 1))?)); k += 1; }
    let (mut dn, mut de) = (Vec::new(), Vec::new());
    let (mut z, mut k) = (d, 0u64);
    while z > d - CLOSURE_SPAN_M { z -= b.thickness(p, u(bed_key(BED_ID_DN, k, 0))?); dn.push(z); de.push(b.eta(u(bed_key(BED_ID_DN, k, 1))?)); k += 1; }
    dn.reverse(); de.reverse();
    dn.extend(up); de.extend(ue);
    Ok(BedTables { pk, glob: dn, glob_eta: de })
}

/// The bed [b0, b1] (depths below the datum) holding depth x inside the package [q0, q1] (depths; q0 = -inf: open upward,
/// q1 = +inf: open downward) and its eta, from the package identity's bed list `t` (beyond it: median thickness, eta
/// midpoint). Anchored at the package top (cut downward), or at its bottom (cut upward) when open upward; None when open both
/// ways. A last bed thinner than P x min joins its neighbour.
pub fn bed_in_package(b: &CtBedParamsV1, p: f64, t: &[(f64, f64)], q0: f64, q1: f64, x: f64) -> Option<(f64, f64, f64)> {
    bed_in_package_k(b, p, t, q0, q1, x).map(|v| (v.0, v.1, v.2))
}

/// Lane PERIOD-ENDS: `bed_in_package` plus the bed index k (counted from the anchor).
pub fn bed_in_package_k(b: &CtBedParamsV1, p: f64, t: &[(f64, f64)], q0: f64, q1: f64, x: f64) -> Option<(f64, f64, f64, usize)> {
    let tmin = p * b.min_permille as f64 / 1000.0;
    let get = |k: usize| -> (f64, f64) { t.get(k).copied().unwrap_or((p * b.median_permille as f64 / 1000.0, b.eta_mid())) };
    if q0.is_finite() {
        let mut a = q0;
        for k in 0..1_000_000usize {
            let (th, e) = get(k);
            let mut e1 = a + th;
            if q1 - e1 < tmin { e1 = q1; }
            if x < e1 || e1 >= q1 { return Some((a, e1, e, k)); }
            a = e1;
        }
        None
    } else if q1.is_finite() {
        let mut a = q1;
        for k in 0..1_000_000usize {
            let (th, e) = get(k);
            let e0 = a - th;
            if x >= e0 { return Some((e0, a, e, k)); }
            a = e0;
        }
        None
    } else { None }
}

/// The eta of the global-table interval holding height z (None beyond the table).
pub fn glob_eta(bt: &BedTables, z: f64) -> Option<f64> {
    let j = bt.glob.partition_point(|v| *v <= z);
    if j == 0 || j >= bt.glob.len() { None } else { bt.glob_eta.get(j - 1).copied() }
}

pub const TREAD_RULE: &str = "ct_tread_v1 (lane TREAD2, apply 2, layered classes): each period [s0, s1] (height T) of the remap gets three limbs climbing up: TREAD at tt (structural bench on the top of the hard bed whose cliff ends at s0: inside thick packages the closure subdivision below; for a one-period var 2 package the unit under its base), TALUS at tb (over the soft interbed and the cliff foot), CLIFF at ta (the hard bed's face, top at s1); height fractions ft + fb + fc = 1 and ft/tt + fb/tb + fc/ta = 1/S (S the remap gradient: steady S, or G with local), so along a v0 slope of gradient S the remapped slope is exactly tt, tb, ta, F is monotone, F(s0) = s0, F(s1) = s1, F <= h and the mean gradient per period is S (relief conserved). Lane TREAD2 measured the other order (talus, cliff, tread reaching s1; tread on the package top): F > h on the tread raised drainage divides and spill points (site2 SITE-3 fill raises up to 56 m), so the tread is the lowest limb. Closure: fc = c2 + w (eta - c2) (c2 the 2-limb cliff fraction, eta the hard-bed share, w the transition), ft = (1/S - fc/ta - (1 - fc)/tb) / (1/tt - 1/tb), fb = 1 - fc - ft; ft > 1 - fc: no talus, fc from (tt, ta). Fallback to the 2-limb profile (counted) where c2 >= eta (the needed cliff exceeds the hard bed: S too steep for a bench) or S >= ta. Unchanged 2-limb: massive classes, the apron period (starting at the hillslope base: not a hard-bed top). Transition (stated closure, no process): w = clamp((S - tb) / (full_tan - tb), 0, 1), identity at the mantle/bare-rock threshold (burial of barely-bare hard beds by talus)";
pub const TREAD_ETA_RULE: &str = "ct_tread_eta_v1: eta = the thickness share of cliff-forming (ct mode != none) units of the post's T0 column (lithology cut: datum, tilt, unit depth bases) inside the period [s0, s1] where both cliff-forming and other units occur there (record); else (one class through the period: the T0 column does not resolve flow tops or interbeds below its unit thickness, or the stand-in lithology) the closure eta0; capped at eta_max. The record share sets the limb heights only, not their order inside the period";
pub const TREAD_ANGLE_RULE: &str = "ct_tread_angle_v1: tt = clamp(tread_tan + a, tread_min_tan, tread_max_tan), a = (dip 1) the record bed-top gradient (the lithology cut's datum tilt of the post's T0 cell, clamped by dip_clamp_ppm) projected on the uphill direction of the crest-rounded v0 surface (central differences): beds rising into the hill steepen the bench (dip slope), beds falling into it flatten it; a counter-slope bench is not representable by the monotone remap (floor tread_min_tan)";

impl CtTreadParamsV1 {
    pub fn v1() -> Self {
        CtTreadParamsV1 {
            on: 1, tread_tan_permille: tanp(8.0), tread_min_tan_permille: tanp(2.0), tread_max_tan_permille: tanp(15.0), dip: 1,
            eta_permille: 700, eta_max_permille: 900, full_tan_permille: tanp(45.0),
            rule: TREAD_RULE.into(), eta_rule: TREAD_ETA_RULE.into(), angle_rule: TREAD_ANGLE_RULE.into(), base_rule: TREAD_BASE_RULE.into(), drain_rule: TREAD_DRAIN_RULE.into(),
            table_source: "project assumptions (lane TREAD2, recalled, unverified, to calibrate): bench tread 8 deg on a flat-lying caprock (structural benches/stripped surfaces carry a thin colluvial wedge, 2-15 deg; Grand Canyon Esplanade/Tonto-type platforms, Deccan and Columbia River trap benches); eta0 0.7 = massive flow interior (colonnade + entablature) share of a flood-basalt flow, flow top/breccia/bole ~0.3 (Long and Wood 1986 GSA Bull. 97:1144, recalled); eta_max 0.9; full transition at 45 deg; stepped-slope geometry of steady parallel retreat (cliff = exposed hard bed, talus over the soft interbed, bench width the remaining span; cf. Schumm and Chorley 1966, Koons 1955, recalled)".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.on != 1 || self.dip > 1 { return Err("ct_tread: on 1, dip 0 or 1".into()); }
        if !(0 < self.tread_min_tan_permille && self.tread_min_tan_permille <= self.tread_tan_permille && self.tread_tan_permille <= self.tread_max_tan_permille) { return Err("ct_tread: 0 < min <= tread <= max tangent".into()); }
        if !(0 < self.eta_permille && self.eta_permille <= 1000 && 0 < self.eta_max_permille && self.eta_max_permille < 1000) { return Err("ct_tread: 0 < eta <= 1000, 0 < eta_max < 1000".into()); }
        if self.full_tan_permille <= 0 { return Err("ct_tread: full_tan > 0".into()); }
        Ok(())
    }
}

/// Lane TREAD2: height fractions (talus fb, cliff fc, tread ft) of one period, or why there is no tread.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TreadSplit { Limbs(f64, f64, f64), Fallback, Zero }

pub fn tread_split(tb: f64, ta: f64, tt: f64, s: f64, cf2: f64, eta: f64, w: f64) -> TreadSplit {
    if cf2 >= 1.0 || s >= ta || cf2 >= eta { return TreadSplit::Fallback; }
    if !(w > 0.0) || !(s > tb) { return TreadSplit::Zero; }
    let fc = cf2 + w * (eta - cf2);
    let ft = (1.0 / s - fc / ta - (1.0 - fc) / tb) / (1.0 / tt - 1.0 / tb);
    if !(ft > 1e-9) { return TreadSplit::Zero; }
    if ft > 1.0 - fc {
        let fc1 = ((1.0 / tt - 1.0 / s) / (1.0 / tt - 1.0 / ta)).clamp(0.0, 1.0);
        return TreadSplit::Limbs(0.0, fc1, 1.0 - fc1);
    }
    TreadSplit::Limbs(1.0 - fc - ft, fc, ft)
}

/// Share of the heights [s0, s1] held by cliff-forming units of a column (unit depth bases below `datum`, codes top-down;
/// the top unit open upward, the bottom one downward; the same cut as `ExposedLithology::at`).
pub fn hard_share(c: &CtParamsV1, datum: f64, base: &[f64], code: &[u8], s0: f64, s1: f64) -> f64 {
    let (d0, d1) = (datum - s1, datum - s0);
    if !(d1 > d0) { return f64::NAN; }
    let n = code.len();
    let mut hard = 0.0;
    for j in 0..n {
        let top = if j == 0 { f64::NEG_INFINITY } else { base[j - 1] };
        let bot = if j + 1 == n { f64::INFINITY } else { base[j] };
        let ov = bot.min(d1) - top.max(d0);
        if ov > 0.0 && c.classes.get(code[j] as usize).map(|k| k.mode != CT_NONE).unwrap_or(false) { hard += ov; }
    }
    hard / (d1 - d0)
}

/// Lane TREAD2: one post's remap with the 3-limb profile.
#[derive(Clone, Copy, Debug)]
pub struct Remap3 { pub z: f64, pub rg: u8, pub fc: f64, pub ft: f64, pub fb: f64, pub limbs3: bool, pub fallback: bool, pub record: bool, pub clamped: bool }

/// The 3-limb remap (lane TREAD2) of one post: `tt` the post's tread tangent, `eta_rec(s0, s1)` the record hard share of the
/// period (None or outside (0, 1): closure). Massive classes, the apron period, a zero tread and the fallback take `remap`.
/// Base level (`TREAD_BASE_RULE`): every result is held at or above min(h, `base_level`) (NaN: none).
pub fn remap3(c: &CtParamsV1, t: &CtTreadParamsV1, tab: &[f64], mode: u32, tb: f64, ta: f64, cf: f64, s: f64, g: Geo, h: f64, lo: f64, hi: f64, tt: f64, base_level: f64, eta_rec: &mut dyn FnMut(f64, f64) -> Option<f64>) -> Remap3 {
    let mut r = remap3_raw(c, t, tab, mode, tb, ta, cf, s, g, h, lo, hi, tt, eta_rec);
    if base_level.is_finite() { let f = h.min(base_level); if r.z < f { r.z = f; r.clamped = true; } }
    r
}

pub const TREAD_DRAIN_RULE: &str = "ct_tread_drain_v1: after the per-post remap, in ascending pre-remap height (index tie-break) every remapped post is raised, if needed, to the lowest D8 neighbour that was lower before the remap (that neighbour's final height): the remap creates no closed depression (rock hillslopes drain; the per-post profile varies between neighbours with their local gradient, hillslope base and tread angle). Along a uniform slope the profile is already monotone and nothing changes; F <= h is kept";
pub const TREAD_BASE_RULE: &str = "ct_tread_base_v1: a remapped post is never lowered below its base level, min(h, max(hillslope base height, sea level)) (the pass-1 D8 source height the post drains to; a marine source counts at sea level, not its sea floor): a hillslope cannot be eroded below the channel, floor, lake or sea it drains to (lane TREAD2 measured posts above coastal hillslope bases taken below sea level by periods reaching under it)";

fn remap3_raw(c: &CtParamsV1, t: &CtTreadParamsV1, tab: &[f64], mode: u32, tb: f64, ta: f64, cf: f64, s: f64, g: Geo, h: f64, lo: f64, hi: f64, tt: f64, eta_rec: &mut dyn FnMut(f64, f64) -> Option<f64>) -> Remap3 {
    let two = |fallback: bool, record: bool| { let (z, rg) = remap(c, tab, mode, tb, ta, cf, s, g, h, lo, hi); Remap3 { z, rg, fc: cf, ft: 0.0, fb: 1.0 - cf, limbs3: false, fallback, record, clamped: false } };
    if mode != CT_LAYERED { return two(false, false); }
    if cf >= 1.0 { return two(true, false); }
    let (s0, s1) = period(c, tab, h, lo, hi);
    if c.apron == 1 && g.base.is_finite() && g.base > s0 && g.base < s1 && h >= g.base { return two(false, false); }
    let tl = s1 - s0;
    if !(tl > 0.0) { return two(false, false); }
    let er = eta_rec(s0, s1).filter(|e| *e > 1e-9 && *e < 1.0 - 1e-9);
    let eta = er.unwrap_or(t.eta_permille as f64 / 1000.0).min(t.eta_max_permille as f64 / 1000.0);
    let tf = t.full_tan_permille as f64 / 1000.0;
    let w = if tf > tb { ((s - tb) / (tf - tb)).clamp(0.0, 1.0) } else { 1.0 };
    let (fb, fc, ft) = match tread_split(tb, ta, tt, s, cf, eta, w) {
        TreadSplit::Limbs(a, b, d) => (a, b, d),
        TreadSplit::Fallback => return two(true, er.is_some()),
        TreadSplit::Zero => return two(false, er.is_some()),
    };
    let u = ((h - s0) / tl).clamp(0.0, 1.0);
    let ut = ft * s / tt;
    let ub = ut + fb * s / tb;
    let (v, rg) = if u < ut { (u * tt / s, RG_TREAD) } else if u < ub { (ft + (u - ut) * tb / s, RG_TALUS) } else { ((ft + fb + (u - ub) * ta / s).min(1.0), RG_CLIFF) };
    Remap3 { z: s0 + tl * v, rg, fc, ft, fb, limbs3: true, fallback: false, record: er.is_some(), clamped: false }
}

pub const TREAD_NOTE: &str = "lane TREAD2 (opt-in ct_tread, rule params.ls.tb.tc.ct.v2.tread.rule): tread height fraction ft of the post's period (permille; 0 where the post has no 3-limb profile): climbing a period [s0, s1] (the remap's) the lowest ft T is tread (bench at ct_tread_tan_permille on the top of the cliff below), the next (1 - ft - c) T talus at the class repose, the upper c T (ct_cliff_permille) cliff; ct_regime 3 marks tread posts; posts held by the base-level or drainage rules keep their limb code";
pub const TREAD_TAN_NOTE: &str = "lane TREAD2 (opt-in ct_tread): tread tangent tt of the post's 3-limb profile (permille; 0 where none); rule params.ls.tb.tc.ct.v2.tread.angle_rule";

pub const BASE_NOTE: &str = "lane STEP2 (opt-in ct_base, rule params.ls.tb.tc.ct.v2.base_rule): the published surface without the cliff-and-talus staircase shift (mm): for pf-detail to re-profile bare rock at its own resolution without double stepping; equals height_mm where ct_regime is 0";

pub fn ct(p: &crate::t1b_rivernet::RiverNetParamsV1) -> Option<&CtParamsV1> { crate::t1b_ls::tc(p).and_then(|t| t.ct.as_ref()) }

/// The bare-rock split at one post: (talus tangent, cliff tangent, cliff height fraction c). None = mantled (v0 law).
#[inline]
pub fn split(c: &CtParamsV1, code: u8, s: f64, rockmass: f64) -> Option<(f64, f64, f64)> {
    let k = &c.classes[code as usize];
    if k.mode == CT_NONE { return None; }
    let tb = k.talus_tan_permille as f64 / 1000.0;
    if !(s > tb) { return None; }
    let ta = (k.cliff_tan_permille as f64 / 1000.0 * if c.rockmass_cliff == 1 { rockmass } else { 1.0 }).max(tb * 1.0001);
    if s >= ta { return Some((tb, s, 1.0)); }
    let cf = ((1.0 / tb - 1.0 / s) / (1.0 / tb - 1.0 / ta)).clamp(0.0, 1.0);
    Some((tb, ta, cf))
}

/// Per-post geometry inputs: hillslope base and crest heights on the pass-1 surface (NaN where unknown).
#[derive(Clone, Copy)]
pub struct Geo { pub base: f64, pub crest: f64 }

/// Band at height z: (gradient, next boundary above z or +inf, regime code).
#[inline]
fn band(c: &CtParamsV1, mode: u32, tb: f64, ta: f64, cf: f64, g: Geo, z: f64) -> (f64, f64, u8) {
    const EPS: f64 = 1e-9;
    if cf >= 1.0 { return (ta, f64::INFINITY, RG_CLIFF); }
    // Massive without a usable base/crest (pass-1 pit or crest post): the layered staircase instead (relief-conserving).
    if mode == CT_MASSIVE && g.base.is_finite() && g.crest.is_finite() && g.crest > g.base {
        let top = g.base + (1.0 - cf) * (g.crest - g.base);
        return if z < top - EPS { (tb, top, RG_TALUS) } else { (ta, f64::INFINITY, RG_CLIFF) };
    }
    let p = c.step_mm as f64 / 1000.0;
    let bench = (1.0 - cf) * p;
    if c.apron == 1 && g.base.is_finite() && z < g.base + bench - EPS { return (tb, g.base + bench, RG_TALUS); }
    let d = c.datum_mm as f64 / 1000.0;
    // Period base s0 <= z < s0 + p, robust to z sitting on a boundary (a climb stops exactly on band boundaries).
    let mut s0 = d + ((z - d) / p).floor() * p;
    if z - s0 >= p - EPS { s0 += p; }
    if z - s0 < -EPS { s0 -= p; }
    if z - s0 < bench - EPS { (tb, s0 + bench, RG_TALUS) } else { (ta, s0 + p, RG_CLIFF) }
}

/// The exact 1-D climb from height zu over horizontal distance dx: the end height.
pub fn climb(c: &CtParamsV1, mode: u32, tb: f64, ta: f64, cf: f64, g: Geo, zu: f64, dx: f64) -> f64 {
    let (mut z, mut rem) = (zu, dx);
    for _ in 0..4096 {
        let (s, nb, _) = band(c, mode, tb, ta, cf, g, z);
        let need = (nb - z) / s;
        if !(need < rem) { return z + s * rem; }
        z = nb;
        rem -= need;
    }
    z + ta * rem
}

/// Eikonal speed at post i from the upwind height zu (None: the v0 speed applies).
#[inline]
pub fn speed(c: &CtParamsV1, rock: &dyn crate::t1b_litho::ExposedLithology, i: usize, zu: f64, dx: f64, mult: f64, rockmass: f64, smin: f64, g: Geo) -> Option<f64> {
    let (sc, _, code) = rock.at(i, zu);
    let s = (sc * mult).max(smin);
    let (tb, ta, cf) = split(c, code, s, rockmass)?;
    let mode = c.classes[code as usize].mode;
    Some((climb(c, mode, tb, ta, cf, g, zu, dx) - zu) / dx)
}

/// Regime, cliff fraction (permille) and massive talus top at post i's own height z (for publication and crest exemption).
pub fn regime(c: &CtParamsV1, rock: &dyn crate::t1b_litho::ExposedLithology, i: usize, z: f64, mult: f64, rockmass: f64, smin: f64, g: Geo) -> (u8, u16, i32) {
    let (sc, _, code) = rock.at(i, z);
    let s = (sc * mult).max(smin);
    let Some((tb, ta, cf)) = split(c, code, s, rockmass) else { return (RG_SOIL, 0, TOP_NONE) };
    let mode = c.classes[code as usize].mode;
    let (_, _, rg) = band(c, mode, tb, ta, cf, g, z);
    let top = if mode == CT_MASSIVE && cf < 1.0 && g.base.is_finite() && g.crest.is_finite() && g.crest > g.base {
        vz_tilekey::chart::quantize_mm((g.base + (1.0 - cf) * (g.crest - g.base)) * 1000.0) as i32
    } else { TOP_NONE };
    (rg, (cf * 1000.0).round().clamp(1.0, 1000.0) as u16, top)
}

pub const REGIME_NOTE: &str = "lane STEP (opt-in ct, rule params.ls.tb.tc.ct.regime_rule / geometry_rule): regime of the post at its eikonal (pass-2) height: 0 soil/creep mantle (v0 law, or the class never forms cliffs, or the post was rebuilt by the coast/volcano stage: over != 0), 1 talus/bench at the class repose angle, 2 cliff/riser at the class rock-mass angle";
pub const CLIFF_NOTE: &str = "lane STEP (opt-in ct): cliff height fraction c of the bare-rock split at the post (permille, 1..1000; 0 where regime is 0): per layered step period P (header) the upper c P is cliff, the lower (1 - c) P talus; massive: the cliff spans the upper c of the hillslope relief";
pub const TOP_NOTE: &str = "lane STEP (opt-in ct): massive classes only: absolute height (mm) of the talus top / cliff foot of the post's hillslope (base + (1 - c)(crest - base) on the pass-1 surface); i32::MIN elsewhere";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_conserves_mean_gradient() {
        let c = CtParamsV1::v1();
        assert!(c.validate().is_ok());
        // basalt (layered), S = tan 50 deg
        let s = 50f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        assert!((cf / ta + (1.0 - cf) / tb - 1.0 / s).abs() < 1e-12);
        assert!(cf > 0.0 && cf < 1.0);
        // mantled below the talus tangent, and never on shale
        assert!(split(&c, 3, 30f64.to_radians().tan(), 1.0).is_none());
        assert!(split(&c, 5, 2.0, 1.0).is_none());
        // steeper than the cliff angle: all cliff at S
        assert_eq!(split(&c, 3, 10.0, 1.0).unwrap().2, 1.0);
    }

    #[test]
    fn layered_climb_is_exact_and_conserves_relief() {
        let c = CtParamsV1::v1();
        let s = 50f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        let g = Geo { base: f64::NAN, crest: f64::NAN };
        // climbing 20 whole periods takes the same horizontal distance as the mean gradient S
        let p = c.step_mm as f64 / 1000.0;
        let h = 20.0 * p / s;
        let z = climb(&c, CT_LAYERED, tb, ta, cf, g, 0.0, h);
        assert!((z - 20.0 * p).abs() < 1e-6, "z {z}");
        // piecewise: composing many small steps equals one big step
        let mut zz = 0.0;
        for _ in 0..1000 { zz = climb(&c, CT_LAYERED, tb, ta, cf, g, zz, h / 1000.0); }
        assert!((zz - z).abs() < 1e-6);
        // the slope distribution is bimodal: only tb and ta occur
        let mut zq = 3.0;
        for _ in 0..500 { let z1 = climb(&c, CT_LAYERED, tb, ta, cf, g, zq, 0.05); let gr = (z1 - zq) / 0.05; assert!((gr - tb).abs() < 1e-6 || (gr - ta).abs() < 1e-6 || (z1 - zq) < 0.05 * ta + 1e-9); zq = z1; }
    }

    #[test]
    fn absent_ct_leaves_the_tc_params_json_unchanged() {
        let off = crate::t1b_ls::TcParamsV1::off();
        let js = serde_json::to_string(&off).unwrap();
        assert!(!js.contains("\"ct\""), "ct must not serialize when absent");
        let mut on = off.clone();
        on.ct = Some(CtParamsV1::v1());
        let js2 = serde_json::to_string(&on).unwrap();
        assert!(js2.contains("\"ct\"") && !js2.contains("window_only"), "window_only 0 must not serialize");
        let back: crate::t1b_ls::TcParamsV1 = serde_json::from_str(&js2).unwrap();
        assert_eq!(back, on);
    }

    #[test]
    fn massive_band_and_apron() {
        let c = CtParamsV1::v1();
        let s = 55f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 1, s, 1.0).unwrap();
        let g = Geo { base: 100.0, crest: 500.0 };
        let top = 100.0 + (1.0 - cf) * 400.0;
        assert_eq!(band(&c, CT_MASSIVE, tb, ta, cf, g, 150.0).2, RG_TALUS);
        assert_eq!(band(&c, CT_MASSIVE, tb, ta, cf, g, top + 1.0).2, RG_CLIFF);
        // base to crest takes exactly (crest - base) / S horizontally
        let z = climb(&c, CT_MASSIVE, tb, ta, cf, g, 100.0, 400.0 / s);
        assert!((z - 500.0).abs() < 1e-6, "z {z}");
        // layered apron: right above a base the slope is talus even inside a riser band
        let sb = split(&c, 3, s, 1.0).unwrap();
        let p = c.step_mm as f64 / 1000.0;
        let base = p - 1.0; // inside the riser of period 0
        let gb = Geo { base, crest: f64::NAN };
        assert_eq!(band(&c, CT_LAYERED, sb.0, sb.1, sb.2, gb, base + 0.5).2, RG_TALUS);
    }

    fn c2(var: u32) -> CtParamsV1 {
        let mut c = CtParamsV1::v1();
        c.apply = 2;
        let mut v = CtV2ParamsV1::v1();
        v.var = var;
        c.v2 = Some(v);
        c
    }

    #[test]
    fn step2_absent_v2_not_serialized_and_validated() {
        let c = CtParamsV1::v1();
        assert!(!serde_json::to_string(&c).unwrap().contains("v2"));
        let mut a2 = CtParamsV1::v1();
        a2.apply = 2;
        assert!(a2.validate().is_err(), "apply 2 needs v2");
        assert!(c2(2).validate().is_ok());
        let mut bad = c2(1);
        bad.apply = 1;
        assert!(bad.validate().is_err(), "var needs apply 2");
        let back: CtParamsV1 = serde_json::from_str(&serde_json::to_string(&c2(2)).unwrap()).unwrap();
        assert_eq!(back, c2(2));
    }

    #[test]
    fn remap_is_monotone_fixes_boundaries_and_gives_two_slopes() {
        let c = c2(0);
        let s = 50f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        let g = Geo { base: f64::NAN, crest: f64::NAN };
        let p = c.step_mm as f64 / 1000.0;
        let f = |h: f64| remap(&c, &[], CT_LAYERED, tb, ta, cf, s, g, h, f64::NEG_INFINITY, f64::INFINITY).0;
        for k in -3..6 { let b = k as f64 * p; assert!((f(b) - b).abs() < 1e-9, "boundary {b} -> {}", f(b)); }
        // along a v0 slope of gradient S the remapped gradient is tb or ta only; F monotone and F <= h
        let dx = 0.01;
        let mut prev = f(0.0);
        let mut x = dx;
        while x < 3.0 * p / s {
            let (h0, h1) = (s * (x - dx), s * x);
            let z1 = f(h1);
            assert!(z1 >= prev - 1e-9 && z1 <= h1 + 1e-9);
            let gr = (z1 - f(h0)) / dx;
            let on_boundary = ((h0 / p).floor() - (h1 / p).floor()).abs() > 0.0;
            let ut_h = (h1 / p).floor() * p + (1.0 - cf) * s / tb * p;
            let on_kink = h0 < ut_h && h1 >= ut_h;
            assert!(on_boundary || on_kink || (gr - tb).abs() < 1e-6 || (gr - ta).abs() < 1e-6, "gr {gr} at h {h1}");
            prev = z1;
            x += dx;
        }
        // the talus share of each period's height is 1 - c
        let ut = (1.0 - cf) * s / tb;
        assert!((f(ut * p) - (1.0 - cf) * p).abs() < 1e-9);
    }

    #[test]
    fn remap_apron_starts_at_the_base_and_massive_band() {
        let c = c2(0);
        let s = 55f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        let p = c.step_mm as f64 / 1000.0;
        let base = 0.6 * p;
        let g = Geo { base, crest: f64::NAN };
        let f = |h: f64| remap(&c, &[], CT_LAYERED, tb, ta, cf, s, g, h, f64::NEG_INFINITY, f64::INFINITY);
        assert!((f(base).0 - base).abs() < 1e-9 && f(base + 0.1).1 == RG_TALUS);
        assert!((f(p).0 - p).abs() < 1e-9);
        // massive: one band from base to crest, ends fixed
        let (tb1, ta1, cf1) = split(&c, 1, s, 1.0).unwrap();
        let gm = Geo { base: 100.0, crest: 500.0 };
        let fm = |h: f64| remap(&c, &[], CT_MASSIVE, tb1, ta1, cf1, s, gm, h, f64::NEG_INFINITY, f64::INFINITY);
        assert!((fm(100.0).0 - 100.0).abs() < 1e-9 && (fm(500.0).0 - 500.0).abs() < 1e-9);
        assert_eq!(fm(110.0).1, RG_TALUS);
        assert_eq!(fm(499.0).1, RG_CLIFF);
    }

    #[test]
    fn closure_table_spacing_and_periods() {
        let c = c2(1);
        let p = c.step_mm as f64 / 1000.0;
        // a deterministic stand-in for the keyed uniform (the real draw is wired in CoInputs::ctx)
        let mut u = |k: u64| -> R<f64> { Ok(((k as f64 * 0.618_033_988_75) % 1.0).abs()) };
        let t = closure_table(&c, &mut u).unwrap();
        assert!(t.windows(2).all(|w| w[1] - w[0] >= 0.5 * p - 1e-9 && w[1] - w[0] <= 1.5 * p + 1e-9));
        assert!(t[0] <= -CLOSURE_SPAN_M && *t.last().unwrap() >= CLOSURE_SPAN_M && t.contains(&0.0));
        // open package: the closure period holds z
        let (s0, s1) = period(&c, &t, 333.0, f64::NEG_INFINITY, f64::INFINITY);
        assert!(s0 <= 333.0 && 333.0 < s1 && t.contains(&s0) && t.contains(&s1));
        // beyond the table: the regular period continues
        let (a, b) = period(&c, &t, 20_000.0, f64::NEG_INFINITY, f64::INFINITY);
        assert!((b - a - p).abs() < 1e-9 && a <= 20_000.0 && 20_000.0 < b);
        // a package <= 1.5 P is one period
        assert_eq!(period(&c, &t, 150.0, 100.0, 100.0 + 1.4 * p), (100.0, 100.0 + 1.4 * p));
        // a thick package: pieces at least 0.5 P from its ends, ends fixed
        let (lo, hi) = (40.0, 40.0 + 6.0 * p);
        let (a, _) = period(&c, &t, lo + 1.0, lo, hi);
        assert_eq!(a, lo);
        let (_, b) = period(&c, &t, hi - 1.0, lo, hi);
        assert_eq!(b, hi);
        let mut z = lo;
        while z < hi { let (a, b) = period(&c, &t, z, lo, hi); assert!(a <= z && z < b && b - a >= 0.5 * p - 1e-9 && a >= lo && b <= hi); z = b; }
    }

    #[test]
    fn packages_merge_thin_layered_units_and_break_on_soft_units() {
        let c = c2(2);
        // top-down: basalt 12, shale 105, basalt 57, carbonate 37, basalt 103, shale 93, basalt (open bottom)
        let th = [12.0, 105.0, 57.0, 37.0, 103.0, 93.0, 500.0];
        let code = [3u8, 5, 3, 4, 3, 5, 3];
        let mut base = Vec::new();
        let mut d = 0.0;
        for t in th { d += t; base.push(d); }
        let pk = packages(&c, &base, &code);
        assert_eq!(pk[0], (f64::NEG_INFINITY, 12.0));
        // 57 + 37 = 94 >= 64 closes one package, the 103 m unit is its own
        assert_eq!(pk[1], (117.0, 211.0));
        assert_eq!(pk[2], (211.0, 314.0));
        assert_eq!(pk[3], (407.0, f64::INFINITY));
        assert_eq!(pk.len(), 4);
    }

    #[test]
    fn local_remap_is_continuous_at_the_repose_threshold() {
        let mut c = c2(0);
        c.v2.as_mut().unwrap().local = 1;
        assert!(c.validate().is_ok());
        let mut bad = CtParamsV1::v1();
        let mut v = CtV2ParamsV1::v1();
        v.local = 1;
        bad.v2 = Some(v);
        assert!(bad.validate().is_err(), "local needs apply 2");
        let g = Geo { base: f64::NAN, crest: f64::NAN };
        let tb = 0.7;
        // a gradient just above repose: c ~ 0 and F ~ identity (no jump against the mantled neighbour)
        let s = tb * 1.000_001;
        let (tb1, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        assert!(cf < 1e-5);
        for h in [3.0, 50.0, 100.0, 127.0, 300.0] {
            let (z, _) = remap(&c, &[], CT_LAYERED, tb1, ta, cf, s, g, h, f64::NEG_INFINITY, f64::INFINITY);
            assert!((z - h).abs() < 1e-3, "h {h} -> {z}");
        }
        assert!(split(&c, 3, tb * 0.99, 1.0).is_none());
    }

    fn c3() -> CtParamsV1 {
        let mut c = c2(0);
        c.v2.as_mut().unwrap().tread = Some(CtTreadParamsV1::v1());
        c
    }

    #[test]
    fn tread_absent_not_serialized_and_validated() {
        let c = c2(2);
        assert!(!serde_json::to_string(&c).unwrap().contains("tread"), "absent tread must not serialize");
        let t = c3();
        assert!(t.validate().is_ok());
        let back: CtParamsV1 = serde_json::from_str(&serde_json::to_string(&t).unwrap()).unwrap();
        assert_eq!(back, t);
        let mut bad = c3();
        bad.apply = 1;
        assert!(bad.validate().is_err(), "tread needs apply 2");
        let mut steep = c3();
        steep.v2.as_mut().unwrap().tread.as_mut().unwrap().tread_max_tan_permille = 800;
        assert!(steep.validate().is_err(), "tread tangent above the talus tangent");
    }

    #[test]
    fn tread_profile_conserves_relief_and_gives_three_slopes() {
        let c = c3();
        let t = c.v2.as_ref().unwrap().tread.clone().unwrap();
        let s = 45f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        let tt = t.tread_tan_permille as f64 / 1000.0;
        let g = Geo { base: f64::NAN, crest: f64::NAN };
        let p = c.step_mm as f64 / 1000.0;
        let f = |h: f64| remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, g, h, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| None);
        let r = f(10.0);
        assert!(r.limbs3 && !r.fallback && !r.record && r.ft > 0.0 && r.fc > cf && r.fc <= 0.7 + 1e-12);
        // fractions: heights sum to one, horizontal spans give the mean gradient S
        let fb = 1.0 - r.fc - r.ft;
        assert!(fb >= 0.0 && (r.ft / tt + fb / tb + r.fc / ta - 1.0 / s).abs() < 1e-12);
        // period boundaries fixed, F monotone, F <= h; along a slope S only tt, tb, ta occur (off the kinks)
        for k in -2..4 { let b = k as f64 * p; assert!((f(b).z - b).abs() < 1e-9); }
        let dx = 0.01;
        let (mut prev, mut x, mut seen) = (f(0.0).z, dx, [false; 3]);
        let ut = r.ft * s / tt;
        let ub = ut + fb * s / tb;
        while x < 2.0 * p / s {
            let (h0, h1) = (s * (x - dx), s * x);
            let (a, b) = (f(h0), f(h1));
            assert!(b.z >= prev - 1e-9 && b.z <= h1 + 1e-9);
            let gr = (b.z - a.z) / dx;
            let u0 = h0 / p - (h0 / p).floor();
            let u1 = h1 / p - (h1 / p).floor();
            let kink = u1 < u0 || (u0 < ut && u1 >= ut) || (u0 < ub && u1 >= ub);
            if !kink {
                if (gr - tt).abs() < 1e-6 { seen[0] = true; assert_eq!(b.rg, RG_TREAD); }
                else if (gr - tb).abs() < 1e-6 { seen[1] = true; assert_eq!(b.rg, RG_TALUS); }
                else { assert!((gr - ta).abs() < 1e-6, "gr {gr}"); seen[2] = true; assert_eq!(b.rg, RG_CLIFF); }
            }
            prev = b.z;
            x += dx;
        }
        assert!(seen.iter().all(|v| *v), "all three limbs along the profile");
        // the tread is the lowest ft of each period (on the top of the cliff below)
        assert!((f(ut * p).z - r.ft * p).abs() < 1e-9);
        assert!(f(0.5 * ut * p).rg == RG_TREAD && f(0.5 * (ut + ub) * p).rg == RG_TALUS);
        // base level: never below min(h, base level)
        let rb = remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, g, 0.4 * p, f64::NEG_INFINITY, f64::INFINITY, tt, 0.3 * p, &mut |_, _| None);
        assert!(rb.clamped && (rb.z - 0.3 * p).abs() < 1e-9);
        let rh = remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, g, 0.2 * p, f64::NEG_INFINITY, f64::INFINITY, tt, 0.3 * p, &mut |_, _| None);
        assert!(rh.z <= 0.2 * p + 1e-9, "never raised above h");
    }

    #[test]
    fn tread_falls_back_when_steep_and_is_identity_at_repose() {
        let c = c3();
        let t = c.v2.as_ref().unwrap().tread.clone().unwrap();
        let g = Geo { base: f64::NAN, crest: f64::NAN };
        let tt = t.tread_tan_permille as f64 / 1000.0;
        // steep: the 2-limb cliff fraction exceeds eta (0.7): the 2-limb profile, flagged
        let s = 70f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        assert!(cf > 0.7);
        let r = remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, g, 50.0, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| None);
        let two = remap(&c, &[], CT_LAYERED, tb, ta, cf, s, g, 50.0, f64::NEG_INFINITY, f64::INFINITY);
        assert!(r.fallback && !r.limbs3 && r.z == two.0 && r.rg == two.1);
        // just above repose: c2 ~ 0 and w ~ 0: identity (continuous with the mantled neighbour)
        let s = 0.7 * 1.000_001;
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        for h in [3.0, 60.0, 127.0] {
            let r = remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, g, h, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| None);
            assert!((r.z - h).abs() < 1e-3, "h {h} -> {}", r.z);
        }
        // the profile is continuous in S across the fallback boundary (c2 = eta: ft -> 0)
        let (tb, ta) = (0.7, 80f64.to_radians().tan());
        let s_eta = 1.0 / (0.7 / ta + 0.3 / tb);
        match tread_split(tb, ta, tt, s_eta * (1.0 - 1e-9), 0.7 - 1e-9, 0.7, 1.0) { TreadSplit::Limbs(_, fc, ft) => assert!(ft < 1e-6 && (fc - 0.7).abs() < 1e-6), TreadSplit::Zero => {}, o => panic!("{o:?}") }
        // massive classes and the apron period keep the 2-limb profile
        let s = 50f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        let gb = Geo { base: 10.0, crest: f64::NAN };
        let r = remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, gb, 20.0, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| None);
        assert!(!r.limbs3 && !r.fallback);
        let r = remap3(&c, &t, &[], CT_MASSIVE, tb, ta, cf, s, Geo { base: 0.0, crest: 400.0 }, 20.0, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| None);
        assert!(!r.limbs3);
    }

    #[test]
    fn tread_record_hard_share() {
        let c = c3();
        // top-down: basalt 40 m, shale 60 m, basalt (open); datum 1000 m
        let base = [40.0, 100.0];
        let code = [3u8, 5, 3];
        assert!((hard_share(&c, 1000.0, &base, &code, 900.0, 1000.0) - 0.4).abs() < 1e-12);
        assert!((hard_share(&c, 1000.0, &base, &code, 1000.0, 1100.0) - 1.0).abs() < 1e-12, "above the top: the top unit");
        assert!((hard_share(&c, 1000.0, &base, &code, 820.0, 940.0) - 80.0 / 120.0).abs() < 1e-12);
        // the record share replaces the closure and moves the cliff fraction
        let t = c.v2.as_ref().unwrap().tread.clone().unwrap();
        let s = 45f64.to_radians().tan();
        let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
        let g = Geo { base: f64::NAN, crest: f64::NAN };
        let tt = t.tread_tan_permille as f64 / 1000.0;
        let r = remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, g, 10.0, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| Some(0.5));
        assert!(r.record && r.limbs3 && (r.fc - 0.5).abs() < 1e-12);
        let r1 = remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, g, 10.0, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| Some(1.0));
        assert!(!r1.record && (r1.fc - 0.7).abs() < 1e-12, "one class through the period: closure");
    }

    // ---------------------------------------------------------------- lane BENCH-VAR (ct_var=3)
    fn cb() -> CtParamsV1 {
        let mut c = c2(3);
        c.v2.as_mut().unwrap().bed = Some(CtBedParamsV1::v1());
        c
    }

    #[test]
    fn bed_absent_not_serialized_and_validated() {
        let c = c2(2);
        let j = serde_json::to_string(&c).unwrap();
        assert!(!j.contains("\"bed\""), "absent bed params must not change the serialized ct params");
        assert!(c2(3).validate().is_err(), "var 3 needs the bed params");
        let mut c = c2(2); c.v2.as_mut().unwrap().bed = Some(CtBedParamsV1::v1());
        assert!(c.validate().is_err(), "bed params need var 3");
        assert!(cb().validate().is_ok());
        let mut c = cb(); c.v2.as_mut().unwrap().bed.as_mut().unwrap().min_permille = 900;
        assert!(c.validate().is_err(), "min <= median");
        let mut c = cb(); c.v2.as_mut().unwrap().bed.as_mut().unwrap().eta_hi_permille = 1_000;
        assert!(c.validate().is_err(), "eta < 1");
        let j = serde_json::to_string(&cb()).unwrap();
        assert!(j.contains(PURPOSE_CT_BED) && j.contains("ct_bed_v1"), "the purpose and rule are declared in the params block");
    }

    #[test]
    fn bed_thickness_law_is_lognormal_clamped_and_varied() {
        assert!((norm_inv(0.5)).abs() < 1e-9 && (norm_inv(0.975) - 1.959964).abs() < 1e-5 && (norm_inv(0.025) + 1.959964).abs() < 1e-5);
        assert!((norm_inv(1e-6) + 4.753424).abs() < 1e-4 && norm_inv(0.0).is_finite() && norm_inv(1.0).is_finite());
        let b = CtBedParamsV1::v1();
        let p = 128.0;
        assert!((b.thickness(p, 0.5) - 0.8 * p).abs() < 1e-9, "median");
        assert_eq!(b.thickness(p, 1e-9), 0.4 * p);
        assert_eq!(b.thickness(p, 1.0 - 1e-9), 1.8 * p);
        // quasi-uniform sample: CV >= 0.4 (the uniform closure U(0.5, 1.5) has 0.29) and the post-weighted mean E[T^2]/E[T]
        // not above the closure's 1.083 P (the bench lowering does not grow)
        let t: Vec<f64> = (0..20_000).map(|k| b.thickness(p, (k as f64 + 0.5) / 20_000.0)).collect();
        let m = t.iter().sum::<f64>() / t.len() as f64;
        let sd = (t.iter().map(|x| (x - m).powi(2)).sum::<f64>() / t.len() as f64).sqrt();
        let e2 = t.iter().map(|x| x * x).sum::<f64>() / t.len() as f64 / m;
        assert!(sd / m > 0.4 && sd / m < 0.5, "cv {}", sd / m);
        assert!(e2 / p < 1.083, "E[T^2]/E[T] {}", e2 / p);
        assert!((b.eta(0.0) - 0.5).abs() < 1e-12 && (b.eta(1.0) - 0.9).abs() < 1e-12 && (b.eta_mid() - 0.7).abs() < 1e-12);
    }

    #[test]
    fn beds_tile_a_package_from_its_top_and_follow_the_bed() {
        let b = CtBedParamsV1::v1();
        let p = 128.0;
        let tmin = 0.4 * p;
        let tab: Vec<(f64, f64)> = (0..40).map(|k| (b.thickness(p, ((k * 37 % 40) as f64 + 0.5) / 40.0), b.eta(((k * 11 % 40) as f64 + 0.5) / 40.0))).collect();
        // a package from depth 1000 to 4953 m: beds contiguous from the top, each >= tmin, the last ends at the package base
        let (q0, q1) = (1000.0, 4953.0);
        let mut x = q0;
        let mut k = 0;
        while x < q1 {
            let (a, e, eta) = bed_in_package(&b, p, &tab, q0, q1, x).unwrap();
            assert!((a - x).abs() < 1e-9 && e > a && e - a >= tmin - 1e-9 && e <= q1 && (0.5..=0.9).contains(&eta));
            // every depth inside the bed maps to the same bed (same bed = same step)
            for f in [0.1, 0.5, 0.9] { assert_eq!(bed_in_package(&b, p, &tab, q0, q1, a + f * (e - a)).unwrap(), (a, e, eta)); }
            if k < tab.len() && e < q1 - 1e-9 { assert!((e - a - tab[k].0).abs() < 1e-9, "bed k has the table thickness"); }
            x = e; k += 1;
        }
        assert!((x - q1).abs() < 1e-9);
        // anchored at the package top: a post's bed is fixed by its depth below the top, so a tilted datum moves beds with the bed
        let shift = 37.0;
        assert_eq!(bed_in_package(&b, p, &tab, q0 + shift, q1 + shift, 1500.0 + shift).map(|r| (r.0 - shift, r.1 - shift)),
                   bed_in_package(&b, p, &tab, q0, q1, 1500.0).map(|r| (r.0, r.1)));
        // open upward: anchored at the bottom, cut upward; open both ways: None (global table)
        let (a, e, _) = bed_in_package(&b, p, &tab, f64::NEG_INFINITY, 900.0, 900.0 - 1.0).unwrap();
        assert!((e - 900.0).abs() < 1e-9 && (e - a - tab[0].0).abs() < 1e-9);
        assert!(bed_in_package(&b, p, &tab, f64::NEG_INFINITY, f64::INFINITY, 0.0).is_none());
        // beyond the table: median thickness
        let (a, e, eta) = bed_in_package(&b, p, &[], 0.0, 1e6, 10.0).unwrap();
        assert!((a - 0.0).abs() < 1e-9 && (e - 0.8 * p).abs() < 1e-9 && (eta - 0.7).abs() < 1e-12);
    }

    #[test]
    fn var3_period_is_the_bed_and_the_remap_keeps_its_ends() {
        let c = cb();
        let p = c.step_mm as f64 / 1000.0;
        // a bed thicker than unit_max x P stays one period under var 3 (var 2 would subdivide it)
        let (lo, hi) = (100.0, 100.0 + 1.8 * p);
        assert_eq!(period(&c, &[], 150.0, lo, hi), (lo, hi));
        assert_ne!(period(&c2(2), &[], 150.0, lo, hi), (lo, hi));
        let (tb, ta) = (35f64.to_radians().tan(), 80f64.to_radians().tan());
        let s = 1.2;
        let cf = (1.0 / tb - 1.0 / s) / (1.0 / tb - 1.0 / ta);
        let g = Geo { base: f64::NAN, crest: f64::NAN };
        assert!((remap(&c, &[], CT_LAYERED, tb, ta, cf, s, g, lo, lo, hi).0 - lo).abs() < 1e-9);
        assert!((remap(&c, &[], CT_LAYERED, tb, ta, cf, s, g, hi - 1e-9, lo, hi).0 - hi).abs() < 1e-6);
        let mut z0 = f64::NEG_INFINITY;
        let mut h = lo;
        while h < hi { let z = remap(&c, &[], CT_LAYERED, tb, ta, cf, s, g, h, lo, hi).0; assert!(z >= z0 - 1e-9 && z <= h + 1e-9); z0 = z; h += 0.5; }
    }

    #[test]
    fn bed_tables_are_keyed_and_global_table_spans_the_datum() {
        let c = cb();
        let mut n = 0u64;
        let mut keys = std::collections::BTreeSet::new();
        let bt = bed_tables(&c, &mut |k| { n += 1; assert!(keys.insert(k), "every counter drawn once"); Ok(((k.wrapping_mul(2654435761) % 1000) as f64 + 0.5) / 1000.0) }).unwrap();
        assert_eq!(bt.pk.len(), BED_IDS as usize);
        assert!(bt.pk.iter().all(|v| v.len() == BED_PER_ID as usize));
        assert!(bt.glob.windows(2).all(|w| w[1] > w[0]) && bt.glob_eta.len() + 1 == bt.glob.len());
        assert!(bt.glob[0] <= -CLOSURE_SPAN_M && *bt.glob.last().unwrap() >= CLOSURE_SPAN_M && bt.glob.contains(&0.0));
        assert!(glob_eta(&bt, 1.0).is_some() && glob_eta(&bt, 1e9).is_none());
        assert!(bed_tables(&c2(2), &mut |_| Ok(0.5)).is_err());
        let _ = n;
    }

    /// Lane PERIOD-ENDS: F through the four knots of one period (fractions [u1, f1, u2, f2]).
    fn f_knots(s0: f64, s1: f64, k: [f64; 4], h: f64) -> f64 {
        let t = s1 - s0;
        let (xs, ys) = ([s0, s0 + t * k[0], s0 + t * k[2], s1], [s0, s0 + t * k[1], s0 + t * k[3], s1]);
        for a in 0..3 {
            if h <= xs[a + 1] || a == 2 { return if xs[a + 1] - xs[a] > 1e-12 { ys[a] + (h - xs[a]) * (ys[a + 1] - ys[a]) / (xs[a + 1] - xs[a]) } else { ys[a + 1] }; }
        }
        unreachable!()
    }

    #[test]
    fn ends_absent_not_serialized_and_validated() {
        let c = cb();
        let j = serde_json::to_string(&c).unwrap();
        assert!(!j.contains("\"ends\"") && !j.contains("ends_rule"), "absent ends must not change the serialized ct params");
        let mut e = cb(); { let v = e.v2.as_mut().unwrap(); v.ends = 1; v.ends_rule = ENDS_RULE.into(); }
        assert!(e.validate().is_ok());
        let back: CtParamsV1 = serde_json::from_str(&serde_json::to_string(&e).unwrap()).unwrap();
        assert_eq!(back, e);
        let mut bad = e.clone(); bad.v2.as_mut().unwrap().ends_rule.clear();
        assert!(bad.validate().is_err(), "ends needs its rule");
        let mut bad = e.clone(); bad.v2.as_mut().unwrap().ends = 2;
        assert!(bad.validate().is_err());
        let mut bad = CtParamsV1::v1(); bad.apply = 1; let mut v = CtV2ParamsV1::v1(); v.ends = 1; v.ends_rule = ENDS_RULE.into(); bad.v2 = Some(v);
        assert!(bad.validate().is_err(), "ends needs apply 2");
    }

    #[test]
    fn remap_period_and_knots_reproduce_remap() {
        // layered (closure table, package, var 3 bed), apron and massive: remap(h) = F on remap_period's ends through knots2
        let tab: Vec<f64> = (-20..40).map(|k| k as f64 * 97.0 + 13.0 * ((k * 7) % 5) as f64).collect();
        let (tb, ta) = (35f64.to_radians().tan(), 80f64.to_radians().tan());
        for (ci, c) in [c2(1), c2(2), cb()].iter().enumerate() {
            for &s in &[0.75, 1.0, 1.6, 3.0] {
                let cf = (1.0 / tb - 1.0 / s) / (1.0 / tb - 1.0 / ta);
                for &(lo, hi) in &[(f64::NEG_INFINITY, f64::INFINITY), (40.0, 190.0), (-300.0, 120.0)] {
                    for &(base, crest) in &[(f64::NAN, f64::NAN), (61.5, 900.0)] {
                        let g = Geo { base, crest };
                        for mode in [CT_LAYERED, CT_MASSIVE] {
                            let mut h = -150.0;
                            while h < 400.0 {
                                // the package is the one holding h (as the rivernet lookup gives it)
                                if h < lo || h >= hi { h += 3.7; continue; }
                                let (z, _) = remap(c, &tab, mode, tb, ta, cf, s, g, h, lo, hi);
                                match remap_period(c, &tab, mode, cf, g, h, lo, hi) {
                                    Some((s0, s1, k)) => {
                                        assert!(s0 <= h && h <= s1 && s1 > s0, "case {ci}: h {h} outside [{s0}, {s1}]");
                                        let zf = f_knots(s0, s1, knots2(cf, tb, s), h);
                                        assert!((z - zf).abs() < 1e-6, "case {ci} mode {mode} k {k}: remap {z} knots {zf} at {h}");
                                        if mode == CT_MASSIVE && base.is_finite() { assert_eq!(k, PK_MASSIVE); }
                                        if k == PK_APRON { assert_eq!(s0, base); }
                                    }
                                    None => assert!((z - h).abs() < 1e-9, "no period but remap moved {h} -> {z}"),
                                }
                                h += 3.7;
                            }
                        }
                    }
                }
            }
        }
        // c >= 1: no period
        assert!(remap_period(&c2(1), &tab, CT_LAYERED, 1.0, Geo { base: f64::NAN, crest: f64::NAN }, 10.0, f64::NEG_INFINITY, f64::INFINITY).is_none());
        // kinds
        let c = c2(2);
        assert_eq!(period_kind(&c, f64::NEG_INFINITY, f64::INFINITY), PK_CLOSURE);
        assert_eq!(period_kind(&c, 0.0, 100.0), PK_PACKAGE);
        assert_eq!(period_kind(&c, 0.0, 1000.0), PK_SUBDIV);
        assert_eq!(period_kind(&cb(), 0.0, 1000.0), PK_PACKAGE);
    }

    #[test]
    fn knots3_reproduce_remap3() {
        let c = c3();
        let t = c.v2.as_ref().unwrap().tread.clone().unwrap();
        let tt = t.tread_tan_permille as f64 / 1000.0;
        let g = Geo { base: f64::NAN, crest: f64::NAN };
        let mut n3 = 0;
        for &deg in &[38.0f64, 42.0, 45.0, 50.0, 60.0] {
            let s = deg.to_radians().tan();
            let (tb, ta, cf) = split(&c, 3, s, 1.0).unwrap();
            let mut h = -200.0;
            while h < 400.0 {
                let r = remap3(&c, &t, &[], CT_LAYERED, tb, ta, cf, s, g, h, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| None);
                let (s0, s1, _) = remap_period(&c, &[], CT_LAYERED, cf, g, h, f64::NEG_INFINITY, f64::INFINITY).unwrap();
                let k = if r.limbs3 { n3 += 1; knots3(r.ft, r.fb, tt, tb, s) } else { knots2(cf, tb, s) };
                let zf = f_knots(s0, s1, k, h);
                assert!((r.z - zf).abs() < 1e-6, "{deg} deg h {h}: remap3 {} knots {zf}", r.z);
                h += 2.3;
            }
        }
        assert!(n3 > 100, "3-limb posts exercised: {n3}");
    }

    #[test]
    fn bed_index_counts_from_the_anchor() {
        let b = CtBedParamsV1::v1();
        let p: f64 = 128.0;
        let t: Vec<(f64, f64)> = (0..10).map(|k| (p * (0.5 + 0.1 * k as f64), 0.7)).collect();
        let mut prev: Option<(f64, usize)> = None;
        let mut x: f64 = 0.0;
        while x < 900.0 {
            let a = bed_in_package_k(&b, p, &t, 0.0, 1000.0, x).unwrap();
            assert_eq!(bed_in_package(&b, p, &t, 0.0, 1000.0, x), Some((a.0, a.1, a.2)));
            if let Some((p0, k0)) = prev { if a.0 != p0 { assert_eq!(a.3, k0 + 1); assert_eq!(a.0, p0 + t[k0].0); } }
            prev = Some((a.0, a.3));
            x += 5.0;
        }
        // open upward: anchored at the bottom, index counted upward
        assert_eq!(bed_in_package_k(&b, p, &t, f64::NEG_INFINITY, 1000.0, 990.0).unwrap().3, 0);
        assert_eq!(bed_in_package_k(&b, p, &t, f64::NEG_INFINITY, 1000.0, 1000.0 - t[0].0 - 1.0).unwrap().3, 1);
    }
}

// E5: contact-driven cap retreat. This is a conditional process experiment on
// native envelope geometry, not a replacement for the historic ct remaps above.
pub const CONTRAST_LABEL: &str = "CONDITIONAL E5: weathering-limited cap retreat; efficient debris evacuation assumed; uncalibrated C7 scenario; not production geology";
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ContrastParams {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subgrid_spacing_mm: Option<i64>,
    pub steps: u32,
    pub dt_yr: i64,
    pub retreat_mm_per_kyr: i64,
    pub support_tan_permille: i64,
    pub weathering_table: crate::t1b_litho::LithologyV1,
    pub rule: String,
    pub evidence: Vec<String>,
}
impl ContrastParams {
    pub fn v1() -> Self {
        let mut weathering_table = crate::t1b_litho::LithologyV1::standin_v1();
        configure_contrast_table(&mut weathering_table);
        Self { subgrid_spacing_mm: None, weathering_table, steps: 96, dt_yr: 10_000, retreat_mm_per_kyr: 2000, support_tan_permille: 404,
            rule: "E5 contact-retreat-v1: envelope sandstone-over-shale only, with declared weak rate > cap rate; fixed/floor/overwritten/bank posts pinned. Bare hard rock has zero soil creep. Adjacent exposed weak rock lowers at its declared effective weathering/evacuation rate, capped by its layer bottom and the lowest neighbour. A cap loses lateral support when the contact stands Sc_weak*dx above a matching-contact weak neighbour; support loss accumulates at the declared lateral retreat rate. Once one grid slab width is undermined the entire exposed cap slab fails to its real basal contact. Hard face weathering acts only at the exposed edge. All steps Jacobi; no mean compensation, fill, uplift or artificial terrain. Removed rock is ledgered as detachment-limited export under an explicit efficient-evacuation assumption; no talus-storage or storm solve. No contrast or no exposed toe means no process. Integer-mm ledger; finite dependency cone 2N+2 posts.".into(),
            evidence: vec![
                "retreat=2000 mm/kyr; interval [1000,3000] observational range, not confidence interval; Book Cliffs sandstone/shale, late Quaternary; McCarroll 2019 USU thesis abstract read 2026-10-25 https://digitalcommons.usu.edu/etd/7667/; evidence measured regional analog, adopted scenario midpoint, not C7 calibration".into(),
                "shale effective lowering=500 mm/kyr; interval unknown; reported toe incision ~0.5 m/kyr, same McCarroll abstract read; mapping incision to local weathering/evacuation is a project assumption, not a measured weathering coefficient".into(),
                "sandstone bare lowering=10 mm/kyr; sensitivity interval [0,50] project scenario range, no confidence claim; no direct source, project assumption (50-fold contrast to shale), conditional".into(),
                "support tangent=0.404 dimensionless (22 degrees); interval [0.325,0.577] recalled-prior sensitivity range (18-30 deg), not observational uncertainty; inherited LithologyV1 shale stand-in, no verified strength calibration; lattice-scale support threshold, conditional".into(),
                "N=96 steps, dt=10000 yr, total=960000 yr; numerical scenario, not world history; allowed N [1,256], one retreat step <= dx; grid slab width=actual dx m, unresolved joints approximated at post scale; all bounds numerical, no empirical claim".into(),
                "cap face may remain vertical while weak substrate erodes; weathering, ravel and runoff jointly control persistence and debris can bury it: Shmilovitz et al 2024 abstract read https://doi.org/10.1029/2023JF007478; evidence mechanism, not calibration; efficient removal is CONDITIONAL and exported rock is not newly minted gameplay stock".into(),
                "vertical-class diagnostic >=75 degrees (tan=3.732050807568877); project diagnostic convention, sensitivity [70,80] degrees, not material strength; 1 mm strict drainage clearance and integer rounding are numerical resolution; 1000 post-mass threshold for partial F6 statistics inherited numerical reporting limit".into(),
            ] }
    }
    pub fn validate(&self) -> R<()> {
        if self.subgrid_spacing_mm.is_some_and(|v|v!=2000) { return Err("ESCARP-1 only qualifies 2000mm spacing".into()); }
        self.weathering_table.validate()?;
        if self.weathering_table.classes.iter().any(|c| c.bare_weathering_mm_per_kyr.is_some_and(|v| v > 10_000)) { return Err("E5 weathering rate outside scenario bounds".into()); }
        if self.steps == 0 || self.steps > 256 || self.dt_yr <= 0 || self.dt_yr > 100_000 || self.retreat_mm_per_kyr <= 0 || self.retreat_mm_per_kyr > 10_000 || self.support_tan_permille <= 0 || self.support_tan_permille > 1000 { return Err("E5: invalid bounded transient parameters".into()); }
        Ok(())
    }
}
pub fn configure_contrast_table(t: &mut crate::t1b_litho::LithologyV1) {
    t.classes[crate::t1b_litho::LITHO_SANDSTONE as usize].bare_weathering_mm_per_kyr = Some(10);
    t.classes[crate::t1b_litho::LITHO_SHALE as usize].bare_weathering_mm_per_kyr = Some(500);
}
#[derive(Clone, Copy, Debug)]
pub struct CapContact {
    pub source: (u32,u16),
    pub base_mm: i64,
    pub top_mm: i64,
    pub weak_bottom_mm: i64,
    pub cap_rate: i64,
    pub weak_rate: i64,
}
/// No fallback column or guessed contact. The RockPosts column adapter is E3's
/// actual tilted b_k. Repeated same-rock layers are NOT a contrast.
pub fn cap_contact(rock: &crate::t1b_rivernet::RockPosts, i: usize, h: f64, rates: &crate::t1b_litho::LithologyV1) -> Option<CapContact> {
    contact_geometry(rock,i,h,rates,false)
}
/// Fine geometry may need a recorded contact whose coarse support lies above
/// its top. Carry the real plane, then qualify the FINE initial height. This
/// never extends the rock upward or invents another source column.
pub fn cap_contact_geometry(rock: &crate::t1b_rivernet::RockPosts, i: usize, h: f64, rates: &crate::t1b_litho::LithologyV1) -> Option<CapContact> {
    contact_geometry(rock,i,h,rates,true)
}
fn contact_geometry(rock: &crate::t1b_rivernet::RockPosts, i: usize, h:f64, rates:&crate::t1b_litho::LithologyV1, extend_support:bool)->Option<CapContact> {
    use crate::t1b_rivernet::RockPosts;
    let RockPosts::Strata { table, .. } = rock else { return None };
    if table.envelope.is_none() { return None; }
    rock.column(i, |_, datum, bases, codes| {
        let mut best: Option<((f64,f64),CapContact)> = None;
        for k in 0..codes.len().saturating_sub(1) {
            if codes[k] != crate::t1b_litho::LITHO_SANDSTONE || codes[k+1] != crate::t1b_litho::LITHO_SHALE { continue; }
            let (Some(cr),Some(wr)) = (rates.class(codes[k]).bare_weathering_mm_per_kyr,rates.class(codes[k+1]).bare_weathering_mm_per_kyr) else { continue };
            if wr <= cr { continue; }
            let base = datum-bases[k]; let top = datum-if k==0 {0.} else {bases[k-1]}; let bottom=datum-bases[k+1];
            if top <= base || (!extend_support && (h > top || h <= bottom)) { continue; }
            // Interior lookup, bounded by actual cap thickness, never a magic epsilon across thin units.
            let Some(source)=rock.layer_at(i,base+(top-base)*0.5) else {continue};
            let q=|x:f64|(x*1000.).round() as i64;
            let c=CapContact {source,base_mm:q(base),top_mm:q(top),weak_bottom_mm:q(bottom),cap_rate:cr,weak_rate:wr};
            if c.top_mm<=c.base_mm {continue;}
            let d=((h-top).max(bottom-h).max(0.),(h-base).abs()); if best.is_none_or(|b|d<b.0) {best=Some((d,c));}
        }
        best.map(|b|b.1)
    }).flatten()
}
#[derive(Clone, Debug, Default, Serialize)]
pub struct ContrastStats {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub subgrid: Option<crate::t1b_escarp::FineSurface>,
    pub contact_posts: u64,
    pub unsupported_post_steps: u64,
    pub failed_slabs: u64,
    pub cap_weathering_m3: f64,
    pub weak_weathering_m3: f64,
    pub collapse_export_m3: f64,
    pub sandstone_vertical_posts: u64,
    pub cap_edge_posts: u64,
    pub vertical_cap_edge_posts: u64,
    pub max_contact_foot_error_m: f64,
    pub max_cap_face_height_m: f64,
    pub initial_relief_m: f64,
    pub final_relief_m: f64,
    pub mean_shift_m: f64,
    pub local_ledger_residual_mm: i64,
    pub conditional: bool,
    pub cell_restriction: Vec<ContrastCellRestriction>,
    /// Removed in-situ geometric volume, per cap/source contact. These are
    /// derived scenario accounts, not J0 collectible mass accounts (density,
    /// porosity, grains and downstream receiving deposits are unresolved).
    pub export_accounts: Vec<ContrastExportAccount>,
    pub export_reservoir: String,
    pub talus_stored_m3: f64,
}
#[derive(Clone, Debug, Serialize)]
pub struct ContrastExportAccount {
    pub source_cell: u32,
    pub cap_layer: u16,
    pub cap_weathering_mm3: String,
    pub weak_weathering_mm3: String,
    pub collapsed_cap_mm3: String,
}
/// The bounded kernel also serves the explicitly requested uniform-lithology
/// falsifier. None contacts is the uniform control; no height changes are legal.
/// No soil mantle exists here: finite weathering supply, not Roering diffusion.
pub fn contrast_kernel(n:usize, dx_mm:i64, pinned:&[bool], contacts:&[Option<CapContact>], p:&ContrastParams, h:&mut [i64]) -> R<ContrastStats> {
    p.validate()?; let nn=n.checked_mul(n).ok_or("E5 lattice overflow")?;
    if n<3 || pinned.len()!=nn || contacts.len()!=nn || h.len()!=nn || dx_mm<=0 {return Err("E5 bad lattice".into());}
    let retreat=p.retreat_mm_per_kyr*p.dt_yr/1000;
    if retreat>dx_mm {return Err("E5 retreat step exceeds one post; reduce dt".into());}
    let mut st=ContrastStats {conditional:true,contact_posts:contacts.iter().filter(|c|c.is_some()).count() as u64,..Default::default()};
    let start_sum:i64=h.iter().sum();
    st.initial_relief_m=(h.iter().max().unwrap()-h.iter().min().unwrap())as f64/1000.;
    let mut debt=vec![0i64;nn]; let mut next=h.to_vec();
    let(mut wc,mut ww,mut fail)=(0i64,0i64,0i64);
    let mut accounts=std::collections::BTreeMap::<(u32,u16),[i64;3]>::new();
    // Four cardinal faces: finite-volume grid slabs; one contact identity per
    // face, no connection through an unrelated cell's coincident lithology.
    for _ in 0..p.steps {
        next.copy_from_slice(h);
        for y in 1..n-1 { for x in 1..n-1 {
            let i=y*n+x; if pinned[i] {continue;} let Some(c)=contacts[i] else {continue};
            let nb=[i-1,i+1,i-n,i+n]; let low=*nb.iter().map(|&j|&h[j]).min().unwrap();
            let floor=(low+1).max(c.weak_bottom_mm+1);
            if h[i]<=floor {continue;}
            let hard=h[i]>c.base_mm;
            let mut toe=false; let mut unsupported=false;
            for j in nb {
                let Some(cj)=contacts[j] else {continue}; if cj.source!=c.source {continue;}
                if hard && h[j]<=cj.base_mm {toe=true; if cj.base_mm-h[j]>p.support_tan_permille*dx_mm/1000 {unsupported=true;} }
                if !hard && h[j]>cj.base_mm {toe=true;}
            }
            if !toe {debt[i]=0; continue;}
            if hard {
                // Abrasion/weathering cannot remove more than the exposed cap;
                // substrate erosion resumes on a later step with its own rate.
                let w=(c.cap_rate*p.dt_yr/1000).min((h[i]-c.base_mm).max(0)).min(h[i]-floor);
                next[i]-=w; wc+=w; accounts.entry(c.source).or_default()[0]+=w;
                if unsupported {
                    st.unsupported_post_steps+=1; debt[i]+=retreat;
                    if debt[i]>=dx_mm {
                        let target=c.base_mm.max(floor);
                        let f=(next[i]-target).max(0); next[i]-=f; fail+=f; accounts.entry(c.source).or_default()[2]+=f;
                        if f>0 {st.failed_slabs+=1;} debt[i]=0;
                    }
                } else {debt[i]=0;}
            } else {
                let w=(c.weak_rate*p.dt_yr/1000).min(h[i]-floor); next[i]-=w; ww+=w; accounts.entry(c.source).or_default()[1]+=w;
            }
        }}
        h.copy_from_slice(&next);
    }
    let end_sum:i64=h.iter().sum(); let a=(dx_mm as f64/1000.).powi(2)/1000.;
    st.cap_weathering_m3=wc as f64*a; st.weak_weathering_m3=ww as f64*a; st.collapse_export_m3=fail as f64*a;
    let mut accounted=[0i64;3];
    let area=(dx_mm as i128).checked_mul(dx_mm as i128).ok_or("E5B area overflow")?;
    for ((cell,layer),v) in accounts {
        if v==[0;3] { continue; }
        for k in 0..3 { accounted[k]=accounted[k].checked_add(v[k]).ok_or("E5B ledger overflow")?; }
        let volume=|h:i64|->R<String>{Ok((h as i128).checked_mul(area).ok_or("E5B volume overflow")?.to_string())};
        st.export_accounts.push(ContrastExportAccount { source_cell:cell,cap_layer:layer,
            cap_weathering_mm3:volume(v[0])?,weak_weathering_mm3:volume(v[1])?,collapsed_cap_mm3:volume(v[2])? });
    }
    if accounted!=[wc,ww,fail] { return Err("E5B per-source export ledger does not close".into()); }
    st.export_reservoir="e5_unresolved_external_sediment: cumulative removed in-situ rock volume under efficient evacuation; outside represented hillslope solid; downstream destination/transport capacity unresolved; no collectible stock".into();
    st.talus_stored_m3=0.;
    st.local_ledger_residual_mm=end_sum-start_sum+wc+ww+fail;
    if st.local_ledger_residual_mm!=0 {return Err("E5 nonzero material ledger".into());}
    st.mean_shift_m=(end_sum-start_sum)as f64/1000./nn as f64;
    st.final_relief_m=(h.iter().max().unwrap()-h.iter().min().unwrap())as f64/1000.;
    for y in 1..n-1 {for x in 1..n-1 {let i=y*n+x;let Some(c)=contacts[i]else{continue};if h[i]<=c.base_mm{continue;}
        let mut edge=false;let mut vertical=false;
        for j in [i-1,i+1,i-n,i+n] {let Some(cj)=contacts[j]else{continue};if cj.source!=c.source||h[j]>cj.base_mm{continue;}
            edge=true; let dz=(h[i]-h[j])as f64;
            if dz/dx_mm as f64>=75f64.to_radians().tan(){vertical=true;}
            st.max_contact_foot_error_m=st.max_contact_foot_error_m.max((h[j]-cj.base_mm).abs()as f64/1000.);
            st.max_cap_face_height_m=st.max_cap_face_height_m.max((h[i]-c.base_mm)as f64/1000.);
        }
        if edge{st.cap_edge_posts+=1;}if vertical{st.vertical_cap_edge_posts+=1;st.sandstone_vertical_posts+=1;}
    }}
    Ok(st)
}
pub fn run_contrast(inp:&crate::t1b_hs::HsIn,p:&ContrastParams,z:&mut[f64])->R<crate::t1b_hs::HsStats>{
    use crate::t1b::{FREE,CHANNEL};
    use std::collections::BTreeMap;
    let lat=inp.lat;let n=lat.n;let nn=lat.nodes();
    if z.len()!=nn || z.iter().any(|v|!v.is_finite()){return Err("E5 nonfinite or wrong-size terrain".into());}
    let start=std::time::Instant::now();
    let h0:Vec<i64>=z.iter().map(|v|(v*1000.).round()as i64).collect();let mut h=h0.clone();
    let mut pinned=vec![false;nn];
    for i in 0..nn {let(x,y)=(i%n,i/n);pinned[i]=x==0||y==0||x+1==n||y+1==n||inp.fixed[i]!=FREE||inp.floor[i]||inp.over[i]!=0;
        if !pinned[i] {for yy in y-1..=y+1 {for xx in x-1..=x+1{if inp.fixed[yy*n+xx]==CHANNEL {pinned[i]=true;}}}}
    }
    let contacts:Vec<_>=(0..nn).map(|i|if p.subgrid_spacing_mm.is_some() {cap_contact_geometry(inp.rock,i,z[i],&p.weathering_table)} else {cap_contact(inp.rock,i,z[i],&p.weathering_table)}).collect();
    if p.subgrid_spacing_mm.is_some() {if let Some(path)=std::env::var_os("VZ_ESCARP_CAPTURE_NATIVE") {
        let target=std::env::var("VZ_ESCARP_CAPTURE_CELL").map_err(|_|"capture needs explicit source cell")?.parse::<u32>().map_err(|_|"bad capture source cell")?;
        let restrictions=contrast_restriction(inp,&h0,&h0,0)?;
        if !restrictions.iter().any(|r|r.cell==target&&r.complete_cell) {return Err("capture source not a complete native cell".into());}
        crate::t1b_escarp::capture_native(inp,&h0,&contacts,&pinned,p,std::path::Path::new(&path),target)?;
        return Err("ESCARP-1 native input captured (diagnostic only); no terrain window published".into());
    }}
    let mut fine_restriction=None;
    let mut cs=if p.subgrid_spacing_mm.is_some() {
        if lat.dx_mm!=32000 {return Err("ESCARP-1 requires 32m native parent".into());}
        // The fine dependency cone is finite. Evaluate only publication plus
        // that cone; never label untouched exterior terrain a whole-cell solve.
        let pad=((2*p.steps*(p.dt_yr/1000)as u32+2).div_ceil(16)+2)as usize;
        let (px0,px1,py0,py1)=inp.pub_box;
        if px0<pad||py0<pad||px1+pad>=n||py1+pad>=n {return Err("ESCARP-1 fine cone exceeds native evaluation support".into());}
        let x0=px0-pad;let y0=py0-pad;let sn=px1-px0+1+2*pad;
        if px1-px0!=py1-py0 {return Err("ESCARP-1 needs square publication".into());}
        let mut sh=Vec::with_capacity(sn*sn);let mut sc=Vec::with_capacity(sn*sn);let mut sp=Vec::with_capacity(sn*sn);
        for y in y0..y0+sn {for x in x0..x0+sn {let i=y*n+x;sh.push(h0[i]);sc.push(contacts[i]);sp.push(pinned[i]);}}
        let (sd,mut st)=crate::t1b_escarp::run(sn,&sh,&sc,&sp,p,[lat.gx0+x0 as i64,lat.gy0+y0 as i64],(pad,pad+px1-px0,pad,pad+py1-py0))?;
        let mut d=vec![0i64;nn];for y in 0..sn {for x in 0..sn {d[(y+y0)*n+x+x0]=sd[y*sn+x];}}
        st.subgrid.as_mut().unwrap().evaluation_origin_parent=[lat.gx0+x0 as i64,lat.gy0+y0 as i64];
        st.subgrid.as_mut().unwrap().evaluation_parent_n=sn;
        let mut round=0.;
        for i in 0..nn {let real=h0[i]as f64+d[i]as f64/256.;h[i]=real.round()as i64;round+=(h[i]as f64-real)*1.024;}
        st.subgrid.as_mut().unwrap().coarse_restriction_rounding_m3=round;
        fine_restriction=Some(d);st
    } else {contrast_kernel(n,lat.dx_mm,&pinned,&contacts,p,&mut h)?};
    let ring=inp.pub_box.0.min(inp.pub_box.2).min(n-1-inp.pub_box.1).min(n-1-inp.pub_box.3)as u32;
    let cone=if p.subgrid_spacing_mm.is_some() {(2*p.steps*(p.dt_yr/1000)as u32+2).div_ceil(16)+2} else {2*p.steps+2};
    // Conservative finite cone. Refuse publication on insufficient support;
    // this is stronger than pretending a window-local boundary is geological.
    if cone>=ring{return Err(format!("E5 dependency cone {cone} posts reaches ring {ring}; increase margin"));}
    let mut st=crate::t1b_hs::HsStats {steps:p.steps,dt_yr:p.dt_yr,t_yr:p.steps as i64*p.dt_yr,
        free_posts:pinned.iter().filter(|b|!**b).count()as u64,pinned_posts:pinned.iter().filter(|b|**b).count()as u64,
        incision_m3:cs.cap_weathering_m3+cs.weak_weathering_m3,slide_moved_m3:cs.collapse_export_m3,slide_export_m3:cs.collapse_export_m3,
        ledger_residual_mm:cs.local_ledger_residual_mm,certificate_pass:true,stencil_cone_posts:cone as u64,cert_ring_distance_posts:ring,
        ..Default::default()};
    let mut cells:BTreeMap<u32,(f64,f64)>=BTreeMap::new();
    for i in 0..nn {let d=h[i]-h0[i];if pinned[i]&&d!=0{return Err("E5 pinned post moved".into());}
        if d!=0 {st.posts_changed+=1;let(x,y)=(i%n,i/n);if x>=inp.pub_box.0&&x<=inp.pub_box.1&&y>=inp.pub_box.2&&y<=inp.pub_box.3{st.published_posts_changed+=1;}
            st.max_lowering_m=st.max_lowering_m.max(-d as f64/1000.);z[i]=h[i]as f64/1000.;}
        for k in 0..4 {let w=inp.hat_w[i][k]as f64;if w>0.{let e=cells.entry(inp.hat_c[i][k]).or_default();e.0+=w*d as f64/1000.;e.1+=w;}}
    }
    for(c,(s,w))in cells{if w<1000.{continue;}if let Some(t)=(inp.half_tol)(c){let sh=s/w;st.f6_max_abs_shift_m=st.f6_max_abs_shift_m.max(sh.abs());if sh.abs()>t{st.f6_cells_over_half_tol.push((c,sh,t,w));}}}
    st.f6_extra_e_update="E5 removes accounted rock; no uplift/mean fit. Whole-cell restriction must be measured independently.".into();
    cs.cell_restriction=contrast_restriction(inp,&h0,&h,cone)?;
    if p.subgrid_spacing_mm.is_some() {
        if let crate::t1b_rivernet::RockPosts::Strata{cell,..}=inp.rock {
            let mut exterior=std::collections::BTreeSet::new();
            let (x0,x1,y0,y1)=inp.pub_box;
            for i in 0..nn {let x=i%n;let y=i/n;if x<x0||x>x1||y<y0||y>y1 {exterior.insert(cell[i]);}}
            for r in &mut cs.cell_restriction {if exterior.contains(&r.cell) {
                r.complete_cell=false;r.transient_cone_clear=false;r.f6_after=false;
                r.qualification=Some("UNQUALIFIED whole-cell F6: native cell extends beyond the fine publication core; mean_after includes unrefined exterior and is not a complete fine-process result".into());
            }}
        }
    }
    // Do not round coarse restriction until native integer publication. The
    // volume ledger and fine output retain exact integer fine removals.
    if let Some(d)=fine_restriction {for i in 0..nn {z[i]=h0[i]as f64/1000.+d[i]as f64/256000.;}}
    st.contrast=Some(cs);st.seconds=start.elapsed().as_secs_f64();Ok(st)
}

/// Whole native-cell restriction on the actual 32 m post-prism support. A cell
/// touching the evaluation boundary is partial, even if its area looks close.
/// Boundary-assignment area uncertainty is reported, not hidden in a fit.
#[derive(Clone, Debug, Serialize)]
pub struct ContrastCellRestriction {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qualification: Option<String>,
    pub cell:u32,
    pub complete_cell:bool,
    pub min_ring_distance_posts:u32,
    pub transient_cone_clear:bool,
    pub sampled_area_m2:f64,
    pub native_area_m2:f64,
    pub boundary_assignment_bound_m2:f64,
    pub mean_before_m:f64,
    pub mean_after_m:f64,
    pub source_mean_m:f64,
    pub f6_tolerance_m:f64,
    pub f6_before:bool,
    pub f6_after:bool,
    pub relief_before_m:f64,
    pub relief_after_m:f64,
}
fn contrast_restriction(inp:&crate::t1b_hs::HsIn, h0:&[i64], h:&[i64],cone:u32)->R<Vec<ContrastCellRestriction>>{
    use std::collections::BTreeMap;
    let crate::t1b_rivernet::RockPosts::Strata{table,cell,..}=inp.rock else{return Err("E5 needs source column support".into())};
    let env=table.envelope.as_ref().ok_or("E5 needs source envelope")?.load()?;
    // area, sum0, sum1, boundary area, min distance, min0,max0,min1,max1
    let mut acc:BTreeMap<u32,(f64,f64,f64,f64,u32,i64,i64,i64,i64)>=BTreeMap::new();
    let lat=inp.lat;let n=lat.n;let dx=lat.dx();
    for i in 0..lat.nodes(){let c=cell[i];let Some(e)=env.cells.get(&c)else{continue};let(x,y)=(i%n,i/n);let(gx,gy)=lat.g(i);
        let area=crate::t1b_checks::envelope_area(gx as f64*dx,gy as f64*dx,dx/2.,dx/2.,e.radius_m).0;
        let dist=x.min(y).min(n-1-x).min(n-1-y)as u32;
        let edge=dist==0||[i.saturating_sub(1),(i+1).min(n*n-1),i.saturating_sub(n),(i+n).min(n*n-1)].iter().any(|&j|j>=cell.len()||cell[j]!=c);
        let a=acc.entry(c).or_insert((0.,0.,0.,0.,u32::MAX,i64::MAX,i64::MIN,i64::MAX,i64::MIN));
        a.0+=area;a.1+=area*h0[i]as f64/1000.;a.2+=area*h[i]as f64/1000.;if edge{a.3+=area;}a.4=a.4.min(dist);a.5=a.5.min(h0[i]);a.6=a.6.max(h0[i]);a.7=a.7.min(h[i]);a.8=a.8.max(h[i]);
    }
    Ok(acc.into_iter().map(|(cell,a)|{let e=&env.cells[&cell];let m0=a.1/a.0;let m1=a.2/a.0;let tol=(0.1*e.relief_m).max(50.);
        ContrastCellRestriction {qualification:None,cell,complete_cell:a.4>0&&(a.0-e.area_m2).abs()<=a.3,min_ring_distance_posts:a.4,transient_cone_clear:a.4>cone,
            sampled_area_m2:a.0,native_area_m2:e.area_m2,boundary_assignment_bound_m2:a.3,mean_before_m:m0,mean_after_m:m1,source_mean_m:e.mean_m,f6_tolerance_m:tol,
            f6_before:(m0-e.mean_m).abs()<=tol,f6_after:(m1-e.mean_m).abs()<=tol,relief_before_m:(a.6-a.5)as f64/1000.,relief_after_m:(a.8-a.7)as f64/1000.}
    }).collect())
}


/// E5B publication only. A resolved face is the largest absolute D8 edge
/// gradient at the post, not a centred derivative that halves a one-cell scarp.
/// Classification uses the exposed unit at the final, millimetre-quantized
/// height. Construction labels on other posts remain construction diagnostics.
/// No height, rock identity, mass or random draw is changed here.
pub fn publish_contrast_cliffs(
    lat: crate::t1b_rivernet::EvalLat, rock: &crate::t1b_rivernet::RockPosts,
    c: &CtParamsV1, z: &[f64], fixed: &[u8], over: &[u8], regime: &mut [u8],
) -> u64 {
    let n=lat.n; let dx=lat.dx_mm as f64 / 1000.; let mut count=0;
    for y in 1..n-1 { for x in 1..n-1 {
        let i=y*n+x;
        if fixed[i]!=crate::t1b::FREE || over[i]!=0 { continue; }
        let h=(z[i]*1000.).round()/1000.;
        let code=crate::t1b_litho::ExposedLithology::at(rock,i,h).2;
        let Some(k)=c.classes.iter().find(|k|k.code==code && k.mode!=CT_NONE) else { continue };
        let mut slope=0.0_f64;
        for yy in y-1..=y+1 { for xx in x-1..=x+1 {
            if xx==x && yy==y { continue; }
            let d=if xx!=x && yy!=y { dx*2f64.sqrt() } else { dx };
            slope=slope.max((h-(z[yy*n+xx]*1000.).round()/1000.).abs()/d);
        }}
        if slope>=k.cliff_tan_permille as f64/1000. { regime[i]=RG_CLIFF; count+=1; }
    }}
    count
}
