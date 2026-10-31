//! Lane DAMMED-VALLEY (2026-10-21): opt-in window-only key `dam_lake=1`; absent = byte-identical.
//!
//! PROBLEM. The coast / volcano stage builds edifices (cone, apron, units) over the T1b surface. Network (sub-threshold)
//! channels under an edifice are buried (lane CO: GUIDE -> FREE, `buried_channel_posts`); T1a channels are never buried
//! (t1a_post skip, corridor cap, vent shadow), so every river at or above the ADR 043 river threshold already keeps its
//! course through or around an edifice (ADR 043 section 3). A buried network valley upstream of the edifice is a closed
//! depression; the stage's final fill (`coast_final` / `volcano_final`) raises it to its spill level and the window
//! publishes it as DRY LAND (lane C7-SITE-CHECK: C7 site 5, a 1.1 km2 exactly flat patch at 516.12 m).
//!
//! RULE (`edifice_dammed_lake_v1`). A dammed lake is an 8-connected component of free posts that the coast / volcano
//! final fill raised (plus connected free posts already exactly at that level), on one published (mm) level, closed (every rim post at or above the level: a fill lake), that (a) touches an edifice post (over = edifice in the component or
//! on its 8-neighbour rim) and (b) also floods pre-edifice ground (at least one post not over = edifice: a valley behind a
//! dam, not a hollow on the cone itself), and (c) is deeper than `min_depth_mm` somewhere (level - height before the final
//! fill). Such a component keeps water by the SAME water balance pf-region already applies to T0 lakes under the climate
//! key (`climate_runoff::LAKE_RULE`, lake_net_evaporation_v1 = ADR 043 section 5 amendment B2 with local climate): it is
//! a lake when its catchment land inflow Q_in covers its net open-water evaporation (PET - P) x lake area. A lake that
//! keeps water stands FULL at its spill level (ADR 043 B2 / F5 "lakes that exist fill to spill"; pf-surface
//! closed-basin-water-balance-v1 "lake full at spill level or terminal below"), which is exactly the fill level the
//! window already published: heights are unchanged, the component's posts are published as fixed = 4 (lake, stage).
//! Its outflow max(Q_in + (P - PET) A, 0) leaves over the spill post (the receiver of the lake post that drains out of
//! the lake on the routed surface); the free-post drainage area is routed with the lake posts free, so the outflow
//! carries the lake's whole catchment downstream (published drainage areas equal the key-off ones).
//! A component that fails the balance would be a terminal lake below its spill (not modelled here): it stays the
//! key-off dry fill and is counted (`dry_by_balance`). Hollows wholly on the edifice (crater / cone-body pits; lane VE's
//! "declared lakes") are counted (`edifice_hollows`) and stay dry fills: out of this rule's scope.
//!
//! INFLOW. Q_in = r_ref x (A_eff(outlet) - sum over the lake of a_post x R / r_ref) with A_eff the effective
//! (Q / r_ref) area routed on the evaluation lattice (lane SEAM-CERT eval_lattice_routing_v1: t1b_lem::route from fixed
//! posts, own area a_post x runoff_weight), summed over the lake posts that drain out of the lake; plus the largest
//! network area of a network channel post on the rim standing above the lake level (a channel whose downstream
//! continuation the edifice buried spills into the lake; approximation, largest only, no double counting along a chain).
//! P and PET: the T0 cell of each lake post (seed-climate-v1, the climate key's artifact). Needs the climate key.

use crate::t1b::{FREE, GUIDE};
use crate::t1b_rivernet::EvalLat;
use serde::{Deserialize, Serialize};

pub const LABEL: &str = "EXPLORATORY lane DAMMED-VALLEY 2026-10-21; window-only";
pub const RULE: &str = "edifice_dammed_lake_v1: an 8-connected component of free posts raised by the coast / volcano final fill (plus connected free posts already at exactly that level), on one published mm level, closed (every 8-neighbour rim post at or above the level), touching an edifice post (component or 8-neighbour rim, over = edifice), flooding at least one post not over = edifice (a valley behind the dam, not a hollow on the cone), and deeper than min_depth somewhere (level - height before the final fill), is published as lake (fixed 4, stage = its fill level = spill) when it keeps water by lake_net_evaporation_v1 (ADR 043 B2 with local climate: catchment land inflow >= (PET - P) x lake area); heights unchanged; free-post routing keeps lake posts free so the outflow (max(Q_in + (P - PET) A, 0)) leaves over the spill post with the whole catchment area";
/// Lane EDIFICE-FOOT (sub-key dam_crater=1 after dam_lake=1; absent = unchanged).
pub const CRATER_RULE: &str = "edifice_crater_lake_v1 (lane EDIFICE-FOOT, dam_crater=1): a final-fill component wholly on the edifice whose deepest post (lowest height before the final fill) is a crater floor or crater wall post (lane VE edifice_style 10 / 7: summit crater, caldera, amphitheatre, maar, pit crater, cinder-cone crater; also VE3 ponded-lava floors, declared) is a crater lake under the same rule (closed, deeper than min_depth, lake_net_evaporation_v1 balance on its own catchment; full at spill when it keeps water); other hollows wholly on the edifice stay dry fills (edifice_hollows)";
pub const INFLOW_RULE: &str = "dam_inflow_v1: Q_in = r_ref x (sum of effective area at lake posts draining out of the lake - sum over the lake of a_post x runoff_weight) on the evaluation-lattice routing (t1b_lem::route from fixed posts, own area a_post x runoff_weight, mm-quantized heights) + r_ref x the largest network effective area of a network channel post on the lake rim standing above the lake level (a channel whose continuation the edifice buried spills in; approximation)";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct DamLakeParamsV1 {
    pub on: u32,
    /// Minimum lake depth (level - height before the final fill, max over the component) [default 1000 mm = ADR 043 F5's
    /// "pits deeper than 1 m" count; project assumption].
    pub min_depth_mm: i64,
    pub rule: String,
    pub inflow_rule: String,
    pub label: String,
    /// Lane EDIFICE-FOOT: 1 = crater lakes (CRATER_RULE); absent (0) = unchanged.
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub crater: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub crater_rule: Option<String>,
}

impl DamLakeParamsV1 {
    pub fn v1() -> Self { DamLakeParamsV1 { on: 1, min_depth_mm: 1000, rule: RULE.into(), inflow_rule: INFLOW_RULE.into(), label: LABEL.into(), crater: 0, crater_rule: None } }
    /// Lane EDIFICE-FOOT: switch the crater-lake extension on.
    pub fn with_crater(mut self) -> Self { self.crater = 1; self.crater_rule = Some(CRATER_RULE.into()); self }
    pub fn validate(&self) -> Result<(), String> {
        if self.on != 1 { return Err("dam_lake: on must be 1 (absent = off)".into()); }
        if self.min_depth_mm < 0 { return Err("dam_lake: min_depth_mm < 0".into()); }
        if self.crater > 1 || (self.crater == 1) != self.crater_rule.is_some() { return Err("dam_lake: crater must be 0 or 1 with its rule".into()); }
        Ok(())
    }
}

/// One candidate component (published or not).
#[derive(Clone, Debug, Default, Serialize)]
pub struct DamLake {
    pub posts: u32,
    pub published_posts: u32,
    pub edifice_posts: u32,
    pub area_km2: f64,
    pub level_m: f64,
    pub max_depth_m: f64,
    pub mean_depth_m: f64,
    pub volume_m3: f64,
    /// Chart position (m) of the lake post nearest the centroid (a label point) and of the spill post.
    pub centroid_chart_m: [f64; 2],
    pub outlet_chart_m: [f64; 2],
    pub outlet_z_m: f64,
    pub outlet_published: bool,
    pub catchment_eff_area_km2: f64,
    pub channel_inflow_eff_area_km2: f64,
    pub inflow_m3_s: f64,
    pub net_evaporation_m3_s: f64,
    pub outflow_m3_s: f64,
    pub p_mm_per_yr: f64,
    pub pet_mm_per_yr: f64,
    pub keeps_water: bool,
    pub published: bool,
    /// Lane EDIFICE-FOOT (dam_crater=1): a crater lake (wholly on the edifice).
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub crater: bool,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct DamLakeStats {
    pub rule: String,
    pub raised_components: u32,
    /// Raised components with a rim post below their level (raised by a later non-fill step, not a lake).
    pub not_closed: u32,
    pub not_on_edifice: u32,
    pub edifice_hollows: u32,
    pub shallow: u32,
    pub dry_by_balance: u32,
    pub lakes_published: u32,
    pub lake_posts_published: u32,
    pub lake_posts_eval: u32,
    /// Lane EDIFICE-FOOT (dam_crater=1): crater lakes published (included in lakes_published) and crater hollows dry by
    /// the balance.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub crater_lakes_published: u32,
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub crater_dry_by_balance: u32,
    /// Lakes (published and dry-by-balance) touching the published lattice, largest first.
    pub lakes: Vec<DamLake>,
}

#[derive(Clone, Copy)]
pub struct DamIn<'a> {
    pub lat: EvalLat,
    /// mm-quantized heights (m) of the finished evaluation-lattice surface.
    pub zq: &'a [f64],
    /// Heights before the coast / volcano final fill (NaN where unset; compared after mm quantization).
    pub z_pre: &'a [f64],
    pub fixed: &'a [u8],
    pub over: &'a [u8],
    /// Lane VE edifice_style per evaluation-lattice post (empty without ve; read only with dam_crater=1).
    pub style: &'a [u8],
    /// Lane EDIFICE-FOOT (ef_cap): posts where the foot drape cap bound (empty otherwise). A draped post is pre-edifice
    /// ground under a thin veneer: it counts as edifice for the touch test but not for the wholly-on-edifice test, so a
    /// valley the remaining cone dams on its draped foot is a dammed valley, not a hollow on the cone.
    pub draped: &'a [bool],
    pub rw: &'a [f32],
    pub chan_area: &'a [f32],
    pub cell: &'a [u32],
    pub p_mm: &'a [f32],
    pub pet_mm: &'a [f32],
    pub r_ref_m_per_yr: f64,
    /// Published box on the evaluation lattice (x0, x1, y0, y1), inclusive.
    pub pub_box: (usize, usize, usize, usize),
}

const D8: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];

/// Classify the final-fill components; returns the lake mask on the evaluation lattice and the stats.
pub fn classify(d: &DamIn, p: &DamLakeParamsV1) -> (Vec<bool>, DamLakeStats) {
    let (n, nn, dx) = (d.lat.n, d.lat.nodes(), d.lat.dx());
    let a_post = dx * dx;
    let min_depth = p.min_depth_mm as f64 / 1000.0;
    let mut st = DamLakeStats { rule: RULE.into(), ..Default::default() };
    let mut lake = vec![false; nn];
    // Raised = published (mm) height above the mm-quantized pre-fill height (raises below half a millimetre are not lakes).
    let zpq: Vec<f64> = d.z_pre.iter().map(|&v| if v.is_finite() { vz_tilekey::chart::quantize_mm(v * 1000.0) as f64 / 1000.0 } else { f64::NAN }).collect();
    let raised: Vec<bool> = (0..nn).map(|i| d.fixed[i] == FREE && zpq[i].is_finite() && d.zq[i].is_finite() && d.zq[i] > zpq[i]).collect();
    if !raised.iter().any(|&b| b) { return (lake, st); }
    // Evaluation-lattice routing (lake posts free), effective areas.
    let r = crate::t1b_lem::route(n, d.zq, d.fixed, dx);
    let mut acc: Vec<f64> = (0..nn).map(|i| a_post * if d.rw.is_empty() { 1.0 } else { d.rw[i] as f64 }).collect();
    for &iu in r.order.iter().rev() { let i = iu as usize; let rc = r.receiver[i]; if rc != crate::NONE { acc[rc as usize] += acc[i]; } }
    let in_pub = |i: usize| { let (x, y) = (i % n, i / n); x >= d.pub_box.0 && x <= d.pub_box.1 && y >= d.pub_box.2 && y <= d.pub_box.3 };
    let nb = |i: usize| { let (x, y) = ((i % n) as i64, (i / n) as i64); D8.iter().filter_map(move |(ox, oy)| { let (xx, yy) = (x + ox, y + oy); if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { None } else { Some(yy as usize * n + xx as usize) } }) };
    let mut comp_id = vec![u32::MAX; nn];
    let mut comp: Vec<usize> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    for s in 0..nn {
        if !raised[s] || comp_id[s] != u32::MAX { continue; }
        let cid = st.raised_components;
        st.raised_components += 1;
        comp.clear();
        comp_id[s] = cid;
        stack.push(s);
        while let Some(i) = stack.pop() {
            comp.push(i);
            // The lake also takes connected free posts already standing exactly at its level (flats an earlier fill or
            // construction left at that level: same water surface, zero depth against the pre-final-fill heights), except
            // a post with a lower neighbour (the spill post stays on the rim).
            for j in nb(i) {
                if comp_id[j] != u32::MAX || d.zq[j] != d.zq[i] { continue; }
                let join = raised[j] || (d.fixed[j] == FREE && zpq[j].is_finite() && !nb(j).any(|k| d.zq[k] < d.zq[i]));
                if join { comp_id[j] = cid; stack.push(j); }
            }
        }
        let level = d.zq[s];
        let mut touches = false;
        let mut closed = true;
        let mut ed = 0u32;
        let mut rim_chan = 0.0f64;
        let (mut maxd, mut sumd) = (0.0f64, 0.0f64);
        let (mut own, mut pm, mut petm, mut nevap) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut a_out, mut best_exit, mut best_acc) = (0.0f64, usize::MAX, -1.0f64);
        let mut npub = 0u32;
        let (mut cx, mut cy) = (0.0f64, 0.0f64);
        for &i in &comp {
            if d.over[i] == crate::t1b_coast::OVER_EDIFICE { touches = true; if d.draped.len() != nn || !d.draped[i] { ed += 1; } }
            let dep = (level - zpq[i]).max(0.0);
            maxd = maxd.max(dep);
            sumd += dep;
            own += a_post * if d.rw.is_empty() { 1.0 } else { d.rw[i] as f64 };
            let c = d.cell[i] as usize;
            let (pp, pe) = (d.p_mm.get(c).copied().unwrap_or(0.0) as f64, d.pet_mm.get(c).copied().unwrap_or(0.0) as f64);
            pm += pp; petm += pe;
            nevap += a_post * (pe - pp) / 1000.0;
            if in_pub(i) { npub += 1; }
            let (gx, gy) = d.lat.g(i);
            cx += gx as f64; cy += gy as f64;
            let rc = r.receiver[i];
            let exits = rc == crate::NONE || comp_id[rc as usize] != cid;
            if exits { a_out += acc[i]; if acc[i] > best_acc { best_acc = acc[i]; best_exit = i; } }
            for j in nb(i) {
                if comp_id[j] == cid { continue; }
                // A fill lake's rim stands at or above its level (the spill post at the level); a post raised by a later
                // non-fill step (levee, foot) has a lower neighbour and is not a lake.
                if d.zq[j].is_finite() && d.zq[j] < level { closed = false; }
                if d.over[j] == crate::t1b_coast::OVER_EDIFICE { touches = true; }
                if d.fixed[j] == GUIDE && d.zq[j] > level && d.chan_area[j].is_finite() { rim_chan = rim_chan.max(d.chan_area[j] as f64); }
            }
        }
        let k = comp.len() as f64;
        if !closed { st.not_closed += 1; continue; }
        if !touches { st.not_on_edifice += 1; continue; }
        let mut crater = false;
        if ed as usize == comp.len() {
            // Lane EDIFICE-FOOT (dam_crater=1): a hollow whose deepest post is a crater floor / wall is a crater lake.
            if p.crater == 1 && d.style.len() == nn {
                let deep = comp.iter().copied().min_by(|&a, &b| zpq[a].total_cmp(&zpq[b]).then(a.cmp(&b))).unwrap();
                crater = d.style[deep] == crate::t1b_ve::ST_FLOOR || d.style[deep] == crate::t1b_ve::ST_WALL;
            }
            if !crater { st.edifice_hollows += 1; continue; }
        }
        if maxd < min_depth { st.shallow += 1; continue; }
        // Water balance (m3/yr): land inflow vs net open-water evaporation.
        let q_in = d.r_ref_m_per_yr * ((a_out - own).max(0.0) + rim_chan);
        let keeps = q_in >= nevap;
        let outflow = (q_in - nevap).max(0.0);
        // Spill post: the receiver of the exit post carrying the largest area (the lake post itself if it drains to nothing).
        let (ox, oy, oz) = if best_exit != usize::MAX {
            let o = match r.receiver[best_exit] { rc if rc != crate::NONE => rc as usize, _ => best_exit };
            let (gx, gy) = d.lat.g(o);
            (gx as f64 * dx, gy as f64 * dx, d.zq[o])
        } else { (f64::NAN, f64::NAN, f64::NAN) };
        let outlet_pub = best_exit != usize::MAX && { let o = match r.receiver[best_exit] { rc if rc != crate::NONE => rc as usize, _ => best_exit }; in_pub(o) };
        let (mcx, mcy) = (cx / k, cy / k);
        let lk = DamLake {
            posts: comp.len() as u32, published_posts: npub, edifice_posts: ed, area_km2: k * a_post / 1e6, level_m: level,
            max_depth_m: maxd, mean_depth_m: sumd / k, volume_m3: sumd * a_post,
            centroid_chart_m: [mcx * dx, mcy * dx], outlet_chart_m: [ox, oy], outlet_z_m: oz, outlet_published: outlet_pub,
            catchment_eff_area_km2: a_out / 1e6, channel_inflow_eff_area_km2: rim_chan / 1e6,
            inflow_m3_s: q_in / crate::SECONDS_PER_YEAR, net_evaporation_m3_s: nevap / crate::SECONDS_PER_YEAR, outflow_m3_s: outflow / crate::SECONDS_PER_YEAR,
            p_mm_per_yr: pm / k, pet_mm_per_yr: petm / k, keeps_water: keeps, published: keeps, crater,
        };
        if keeps {
            for &i in &comp { lake[i] = true; }
            st.lakes_published += 1;
            st.lake_posts_eval += comp.len() as u32;
            st.lake_posts_published += npub;
            if crater { st.crater_lakes_published += 1; }
        } else { st.dry_by_balance += 1; if crater { st.crater_dry_by_balance += 1; } }
        if npub > 0 { st.lakes.push(lk); }
    }
    st.lakes.sort_by(|a, b| b.posts.cmp(&a.posts).then(a.level_m.total_cmp(&b.level_m)));
    (lake, st)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::t1b::{FREE, GUIDE};

    fn lat(n: usize) -> EvalLat { EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 } }

    /// A valley floor at 10 m behind a 30 m edifice dam with a notch at 20 m: the filled valley (level 20 m) is a lake when
    /// P > PET, dry when evaporation wins; a pit wholly on the edifice is a hollow.
    #[test]
    fn dammed_valley_lake_and_balance() {
        let n = 12;
        let nn = n * n;
        let mut zq = vec![40.0f64; nn];
        let mut zpre = vec![f64::NAN; nn];
        let _ = &mut zpre;
        let mut over = vec![0u8; nn];
        let mut fixed = vec![FREE; nn];
        // outlet: a sea post column on the east edge at 0 m
        for y in 0..n { fixed[y * n + n - 1] = crate::t1b::MARINE; zq[y * n + n - 1] = 0.0; }
        // valley rows 4..8, cols 2..7 at 10 m pre-fill, filled to 20 m
        for y in 4..8 { for x in 2..7 { let i = y * n + x; zpre[i] = 10.0; zq[i] = 20.0; } }
        // dam column x = 7 (edifice), notch at y = 5 down to 20 m, then a slope to the sea
        for y in 0..n { over[y * n + 7] = crate::t1b_coast::OVER_EDIFICE; zq[y * n + 7] = 30.0; }
        zq[5 * n + 7] = 20.0;
        for x in 8..n - 1 { zq[5 * n + x] = 20.0 - (x - 7) as f64 * 4.0; }
        // Pre-fill heights everywhere (outside the valley equal to the final surface, as the stage records them).
        for i in 0..nn { if zpre[i].is_nan() { zpre[i] = zq[i]; } }
        // A flat post at the lake level beside the lake (an earlier fill's flat) joins it; the notch (spill) does not.
        zq[3 * n + 3] = 20.0; zpre[3 * n + 3] = 20.0;
        let rw = vec![1.0f32; nn];
        let ca = vec![0.0f32; nn];
        let cell = vec![0u32; nn];
        let p = DamLakeParamsV1::v1();
        let wet = DamIn { lat: lat(n), zq: &zq, z_pre: &zpre, fixed: &fixed, over: &over, style: &[], draped: &[], rw: &rw, chan_area: &ca, cell: &cell,
            p_mm: &[1300.0], pet_mm: &[260.0], r_ref_m_per_yr: 1.0, pub_box: (0, n - 1, 0, n - 1) };
        let (m, s) = classify(&wet, &p);
        assert_eq!(s.lakes_published, 1, "{s:?}");
        assert_eq!(m.iter().filter(|&&b| b).count(), 21);
        let l = &s.lakes[0];
        assert!((l.level_m - 20.0).abs() < 1e-9 && (l.max_depth_m - 10.0).abs() < 1e-9);
        assert_eq!(l.outlet_chart_m, [7.0 * 32.0, 5.0 * 32.0]);
        // Arid: net evaporation 2 m/yr over 20 posts vs land inflow from the remaining catchment.
        let dry = DamIn { p_mm: &[100.0], pet_mm: &[2100.0], ..wet };
        let (m2, s2) = classify(&dry, &p);
        let l2 = &s2.lakes[0];
        assert_eq!(l2.keeps_water, l2.inflow_m3_s >= l2.net_evaporation_m3_s);
        if !l2.keeps_water { assert_eq!(s2.dry_by_balance, 1); assert!(!m2.iter().any(|&b| b)); }
        // A hollow wholly on the edifice is not a dammed valley.
        let mut over3 = over.clone();
        for y in 4..8 { for x in 2..7 { over3[y * n + x] = crate::t1b_coast::OVER_EDIFICE; } }
        over3[3 * n + 3] = crate::t1b_coast::OVER_EDIFICE;
        let h = DamIn { over: &over3, ..wet };
        let (m3, s3) = classify(&h, &p);
        assert_eq!((s3.edifice_hollows, s3.lakes_published), (1, 0));
        assert!(!m3.iter().any(|&b| b));
        // Lane EDIFICE-FOOT: with dam_crater=1 the same hollow is a crater lake when its deepest post is a crater floor,
        // and stays a dry hollow when it is not (e.g. a lahar pit).
        let mut sty = vec![crate::t1b_ve::ST_BODY; nn];
        for y in 4..8 { for x in 2..7 { sty[y * n + x] = crate::t1b_ve::ST_FLOOR; } }
        let pc = DamLakeParamsV1::v1().with_crater();
        pc.validate().unwrap();
        let (m4, s4) = classify(&DamIn { style: &sty, ..h }, &pc);
        assert_eq!((s4.edifice_hollows, s4.lakes_published, s4.crater_lakes_published), (0, 1, 1), "{s4:?}");
        assert!(s4.lakes[0].crater && m4.iter().filter(|&&b| b).count() == 21);
        let lahar = vec![crate::t1b_ve::ST_LAHAR; nn];
        let (_, s5) = classify(&DamIn { style: &lahar, ..h }, &pc);
        assert_eq!((s5.edifice_hollows, s5.lakes_published), (1, 0));
        // Without the sub-key the params serialize exactly as before (recipe digests unchanged).
        assert!(!serde_json::to_string(&p).unwrap().contains("crater"));
        let _ = GUIDE;
    }
}
