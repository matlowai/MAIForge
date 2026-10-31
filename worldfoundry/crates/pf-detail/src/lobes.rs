//! Lane DE DE5 (opt-in `de_lobes=1`, needs de=1): edifice flow lobes and units by eruption style, so a volcanic cone is built of
//! units with lateral extent instead of one sheet. PREREG: `unpublished-design-record` (DE5).
//!
//! design 2026-10-15: "consider real volcanos and their types of erruptions".
//!
//! a) Recorded: T1b's own lava-flow lobes (window header `prep.net.volcanoes[].cf.flows`), exact T1b geometry (radial axis from
//!    the vent with T1b's wander, cross-section t sqrt(1 - (lat/w)^2), half-width growing x1.5 downstream, 200 m ramp-in,
//!    semicircular snout). They are the youngest units where they lie.
//! b) Keyed sub-resolution units (declared closure, `t2.closure.lobes.v1`): sources keyed per 256 m chart lattice cell (a pure
//!    function of the volcano id and the cell, so independent of which tiles are requested); each unit is traced DOWNHILL on the
//!    T1b bicubic surface from its source (20 m steps) for its keyed runout, so units follow the slope and bend with the topography;
//!    lens cross-section and snout as T1b's flows. Unit kinds and sizes by eruption style (VE's per-post `edifice_style` when the
//!    window publishes it, else the T1b volcano kind): arc stratovolcano (Vulcanian / Plinian, andesite): short thick blocky lava
//!    lobes, pyroclastic-density-current and lahar lobes (weak), tephra fall mantles; shield (Hawaiian / Strombolian basalt): thin
//!    long lobes; large igneous province: broad sheet lobes. Size ranges are recalled textbook values (unverified): andesite flows
//!    10-100 m thick and 0.2-1 km wide (T1b's resolved ones), so the closure draws the smaller ones; pahoehoe / aa lobes 0.5-10 m;
//!    flood-basalt sheet lobes 10-50 m thick and km wide; tephra thickness decays about exponentially with distance.
//! c) Column at a post: covering units sorted youngest first, stacked down from the construction surface D = max(S, T1b bicubic)
//!    (S = the BS 288 m box mean); each lava unit is split into rubbly top / dense core / basal breccia (BS LavaCfg fractions);
//!    strength multipliers keyed per unit and normalised per style so the expected thickness-weighted mean of 1/Sc is the class
//!    value (exact in expectation, approximate per post); below the generated stack, the BS per-volcano 1D sequence.

use crate::bedstyle::{FAC_LAVA_BASE, FAC_LAVA_CORE, FAC_LAVA_TOP, FAC_PYRO, FAC_PYRO_HARD};
use crate::input::T1bWindow;
use crate::{sc_of, R};
use serde::Serialize;
use std::collections::HashMap;
use vz_keyed_rng::Scope;

pub const PURPOSE_LOBES: &str = "t2.closure.lobes.v1";
/// Lane RIB (opt-in de_lobate=1): purpose of the lobate-margin draws (separate stream: the unit draws are unchanged).
pub const PURPOSE_LOBATE: &str = "t2.closure.lobe_margins.v1";
pub const LOBATE_NOTE: &str = "lane RIB de_lobate=1: lens units (keyed lobes, VE paths, T1b recorded flows) get lobate margins: each side's half-width is modulated along the unit by three keyed sinusoids (wavelengths 1.5-4 x the unit's own half-width, at least 60 m, so scaled by the style and the record), max relative change amp; real lava lobes advance by breakouts and fingering at fronts and margins and pyroclastic / lahar deposits finger, so margins are not straight; a declared statistical sub-resolution closure (purpose t2.closure.lobe_margins.v1), not simulated emplacement. VE Path polylines (128 m vertices) are Catmull-Rom resampled to 16 m so their margins have no kinks.";
pub const LOBE_NOTE: &str = "lane DE de_lobes=1: edifice units with lateral extent. T1b's recorded lava-flow lobes (window header cf.flows) are the youngest units where T1b drew them (T1b geometry); below and between them a keyed sub-resolution closure of eruption units by style (arc: blocky lava lobes, pyroclastic-density-current / lahar lobes, tephra mantles; shield: thin long lobes; LIP: sheet lobes) whose sources are keyed per 256 m chart cell and whose paths are traced downhill on the T1b surface; lens cross-sections, snouts; per-unit strength multipliers normalised per style (class mean 1/Sc in expectation). Not recorded eruptions: no eruption history is simulated.";

#[derive(Clone, Debug, Serialize)]
pub struct LobeParams {
    pub version: u32,
    /// Source lattice cell (m).
    pub cell_m: f64,
    /// Target mean number of units covering a flank point (sets the source density per style).
    pub coverage: f64,
    /// Path step (m) and the stop slope (a unit stops spreading on flat ground).
    pub step_m: f64,
    pub stop_slope: f64,
    /// Lane DE2 F4 (opt-in de_ve=1): lane VE's units (Path, Cone, Dome, Fissure) as recorded units when the window lists them.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ve: Option<bool>,
    /// Lane RIB (opt-in de_lobate=1): lobate lens margins and smoothed VE paths. None = off (params JSON unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lobate: Option<LobateParams>,
    /// Lane DATUM (opt-in de_stack=1): level stacking of the column (contacts from laterally smoothed unit thicknesses, younger
    /// units fill where an older unit is absent). None = off (top-down stacking of the sharp lenses; params JSON unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stack: Option<StackParams>,
    /// Lane DATUM (opt-in de_route=1): keyed unit paths traced on the T1b box-mean surface (the BS 288 m datum S) instead of the
    /// T1b bicubic. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub route: Option<RouteParams>,
    /// Lane DATUM2 (opt-in de_seq=1|2): sequential (compensational) emplacement of the keyed units. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub seq: Option<SeqParams>,
    /// Lane PYRO2 (opt-in de_tephra_rec=1, needs pyro_sub=1): tephra mantles from the recorded explosive eruptions. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tephra_rec: Option<String>,
    /// Lane LOBE (opt-in de_lobe_relief=1, needs de_ve=1): constructional relief of lane VE's lava flows. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relief: Option<crate::relief::ReliefParams>,
}

pub const SEQ_NOTE: &str = "lane DATUM2 de_seq: sequential (compensational) emplacement. All keyed units of a volcano whose sources lie in the T1b window (not only those near the block, so the result is a pure function of the window and the volcano: every block computes the same units and paths) are emplaced in age order, oldest first. Each unit is traced downhill (steepest descent, step_m steps) on the routing surface Z = C + B, where C is the bicubic of the box mean of the T1b surface over (2 h + 1)^2 posts (h = box_half_posts; the pre-wave-cut surface under ed_datum=1), a proxy for the constructional cone without its young sub-300 m valleys, and B is the thickness built so far by the older units on a cell_m lattice over the window (each unit's lens rasterised after it stops); a step that would climb is halved up to three times; where every trial climbs (ponded behind older units in a closed low) the flow is taken to inflate and overtop: that step follows the routing base C alone, so drawn runouts are kept (opt-in de_seq_pond=1: the unit stops there instead, its runout cut). Later units therefore avoid the piles of earlier ones and fill lows (compensational stacking). Mode 1 also builds the column bottom-up: the stack base lies at depth Pbar below D, Pbar the final built thickness of all units (keyed, recorded T1b flows, VE units, tephra falls) on the lattice smoothed with a Gaussian of base_sigma_m; recorded units (T1b cf flows, VE units, the youngest) hang from D as the record draws them; the covering keyed units and tephra falls are stacked oldest first on that base with their own (sharp) local thickness, so every unit's base is the top of what lies below it (flat or draped) and its top carries its own lens relief; keyed tops above the recorded units' base (or D) are truncated (eroded before the recorded flows covered them); where the keyed pile ends below that level the youngest present keyed unit extends up to it. Mode 2 changes only the routing (the column is stacked as without de_seq: hanging, de_stack or de_lens). No new keyed draws; unit sources, sizes and materials are unchanged, only the paths.";

#[derive(Clone, Debug, Serialize)]
pub struct SeqParams {
    pub version: u32, pub mode: u8, pub cell_m: f64, pub box_half_posts: i64, pub base_sigma_m: f64, pub note: String,
    /// Opt-in de_seq_pond=1: a ponded unit stops (its runout is cut) instead of overtopping. None = overtop.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pond_stop: Option<bool>,
    /// Lane PYRO2 opt-in de_seq_fill=1 (mode 1): the gap above the keyed pile is filled by the youngest THICK keyed unit. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fill: Option<SeqFillParams>,
    /// Lane ARCH opt-in de_seq_level=1|2 (mode 1): keyed units levelled in the column instead of stacked at their sharp lens
    /// thickness (SEQ_LEVEL_NOTE). None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub level: Option<SeqLevelParams>,
}
impl SeqParams { pub fn new(mode: u8) -> Self { SeqParams { version: 1, mode, cell_m: 16.0, box_half_posts: 4, base_sigma_m: 200.0, note: SEQ_NOTE.into(), pond_stop: None, fill: None, level: None } } }

pub const SEQ_LEVEL_NOTE: &str = "lane ARCH de_seq_level=1|2 (de_seq=1): without it the de_seq column stacks every keyed unit from the stack base at its sharp lens thickness, so a contact's depth is the base plus the sum of the older units' sharp lens cross-sections (elliptic, vertical tangent at the margin): every older lobe's dome is inherited by all younger contacts (nested arches, anticline-like) and every lens margin is a near-vertical step. Flows do not do that: a lava flow or pyroclastic density current spreading over an uneven pile follows it at its base but levels it (thick in lows, thin or absent over older highs); only a lava lobe's own top keeps a dome of its own. With de_seq_level the keyed units (oldest first, from the same stack base at depth Pbar below D) are placed as lane DATUM2's de_lens places them in the level stack: level tops are the cumulative laterally smoothed thicknesses measured up from the stack base (each unit's lens cross-section convolved across its axis with a Gaussian of sigma = blur_hw_mult x half-width, clamped to [blur_min_m, blur_max_m]; volume per unit kept laterally); a unit is present where its sharp lens thickness is > 0.05 m, its base is the top of the material placed below it; a lava unit's top is its level top raised by its own lens excess (profile minus smoothed thickness: a dome over its axis, a pinch-out to its margins; profile 1 elliptic as drawn, 2 flat-topped lobe tc sqrt(1 - (l / w)^4)), a pyroclastic density current / lahar fills to its level top (absent where an older dome rises above it: onlap), a tephra fall drapes the material below at its own thickness; towards its sharp footprint edge a lava or pyroclastic density current / lahar unit is no thicker than its distance to that edge times tan(margin_deg) (35 deg for lava: flow-front / levee rubble at repose; 15 deg for pyroclastic density current / lahar lobe margins; closure values, not calibrated), so units pinch out over a ramp instead of ending in a near-vertical step. Younger units therefore do not inherit older domes. The recorded units hanging from D, the truncation at their base and the gap fill (or de_seq_fill) are as in de_seq=1; routing, Pbar, unit sources, sizes and materials are unchanged; no new draws.";

pub const SEQ_STRETCH_NOTE: &str = "lane ARCH de_seq_stretch=1 (de_seq_level): where the levelled keyed pile at a post ends below the recorded units' base (or D) - the stack base lies Pbar below D, Pbar a 200 m smoothed total, so the local pile can be thinner - the gap is closed by stretching the whole keyed pile between the stack base and that level by one factor (every keyed unit thickened in proportion to its levelled local thickness, order and contacts kept) instead of thickening one unit (de_seq=1: the youngest present; de_seq_fill: the youngest thick one). The gap and the pile are smooth laterally, so the stretch is too, and no near-vertical contact appears where the thickened unit would switch from one post to the next. A closure (the unrepresented volume of the local units, spread in proportion), not a draw; truncation at the recorded base is unchanged; de_seq_fill is not used when this is on.";

#[derive(Clone, Debug, Serialize)]
pub struct SeqLevelParams { pub version: u32, pub mode: u8, pub blur_hw_mult: f64, pub blur_min_m: f64, pub blur_max_m: f64,
    /// Margin slope (deg) capping a lava unit's / a pyroclastic density current or lahar unit's thickness towards its sharp
    /// footprint edge: thickness <= (edge distance) x tan(slope).
    pub margin_deg: [f64; 2], pub note: String,
    /// Opt-in de_seq_stretch=1: close the gap above the levelled keyed pile by stretching the whole keyed pile (SEQ_STRETCH_NOTE).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stretch: Option<String> }
impl SeqLevelParams {
    pub fn new(mode: u8) -> Self { let d = StackParams::default(); SeqLevelParams { version: 1, mode, blur_hw_mult: d.blur_hw_mult, blur_min_m: d.blur_min_m, blur_max_m: d.blur_max_m, margin_deg: [35.0, 15.0], note: SEQ_LEVEL_NOTE.into(), stretch: None } }
    /// The level-stack parameters the column uses (de_stack blur, de_lens profile).
    fn stack(&self) -> StackParams { StackParams { version: 1, blur_hw_mult: self.blur_hw_mult, blur_min_m: self.blur_min_m, blur_max_m: self.blur_max_m, note: String::new(), lens: Some(LensParams::new(self.mode)) } }
}

pub const SEQ_FILL_NOTE: &str = "lane PYRO2 de_seq_fill=1 (de_seq=1): where the keyed pile at a post ends below the recorded units' base (or D), the gap g is not filled by stretching the youngest present keyed unit (which is often a centimetre-to-decimetre tephra mantle far from the vent, shown ~10 m thick) but by the youngest present keyed unit whose own deposited local thickness is at least min_ratio x g (a unit thick enough that a modest local inflation, ponding or the unrepresented thickness of the same flow field can plausibly account for the gap); if none is, the thickest present keyed unit. That unit is thickened in place (its top raised by g); the younger keyed units above it keep their own deposited thickness and are raised with it, so the order of units is kept and thin mantles drape at their true thickness. Deposited thickness, sources, paths and the stack base are unchanged; no new draws.";

#[derive(Clone, Debug, Serialize)]
pub struct SeqFillParams { pub version: u32, pub min_ratio: f64, pub note: String }
impl Default for SeqFillParams { fn default() -> Self { SeqFillParams { version: 1, min_ratio: 0.5, note: SEQ_FILL_NOTE.into() } } }

/// Lane DATUM2: a scalar field on a regular lattice over the T1b window (the built thickness of the sequential pass).
#[derive(Clone, Debug)]
pub struct SeqField { pub x0: f64, pub y0: f64, pub cell: f64, pub nx: usize, pub ny: usize, pub v: Vec<f64> }
impl SeqField {
    fn new(dom: (f64, f64, f64, f64), cell: f64) -> SeqField {
        let nx = ((dom.2 - dom.0) / cell).ceil() as usize + 1; let ny = ((dom.3 - dom.1) / cell).ceil() as usize + 1;
        SeqField { x0: dom.0, y0: dom.1, cell, nx, ny, v: vec![0.0; nx * ny] }
    }
    #[inline]
    pub fn bilin(&self, x: f64, y: f64) -> f64 {
        let fx = ((x - self.x0) / self.cell).clamp(0.0, (self.nx - 1) as f64); let fy = ((y - self.y0) / self.cell).clamp(0.0, (self.ny - 1) as f64);
        let (i, j) = ((fx.floor() as usize).min(self.nx - 2), (fy.floor() as usize).min(self.ny - 2));
        let (a, b) = (fx - i as f64, fy - j as f64);
        let k = j * self.nx + i;
        self.v[k] * (1.0 - a) * (1.0 - b) + self.v[k + 1] * a * (1.0 - b) + self.v[k + self.nx] * (1.0 - a) * b + self.v[k + self.nx + 1] * a * b
    }
    #[inline] fn xy(&self, i: usize, j: usize) -> (f64, f64) { (self.x0 + i as f64 * self.cell, self.y0 + j as f64 * self.cell) }
    /// Cell index range covering [a, b] along x (or y with the y origin), clamped.
    fn irange(o: f64, cell: f64, n: usize, a: f64, b: f64) -> (usize, usize) {
        let lo = ((a - o) / cell).ceil().max(0.0) as usize; let hi = (((b - o) / cell).floor().max(-1.0) + 1.0) as usize;
        (lo.min(n), hi.min(n))
    }
    /// Add a path unit's sharp lens thickness (cells within its reach; segment ranges per cell as the 32 m grid does).
    fn add_path(&mut self, u: &Unit, amp: f64) {
        let pad = 1.5 * u.hw * (1.0 + amp) + self.cell;
        let (mut bx0, mut by0, mut bx1, mut by1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
        for p in &u.path { bx0 = bx0.min(p[0]); by0 = by0.min(p[1]); bx1 = bx1.max(p[0]); by1 = by1.max(p[1]); }
        let (i0, i1) = Self::irange(self.x0, self.cell, self.nx, bx0 - pad, bx1 + pad);
        let (j0, j1) = Self::irange(self.y0, self.cell, self.ny, by0 - pad, by1 + pad);
        if i0 >= i1 || j0 >= j1 { return; }
        let lw = i1 - i0;
        let mut rng = vec![(u16::MAX, 0u16); lw * (j1 - j0)];
        for (si, w) in u.path.windows(2).enumerate() {
            let si = si.min(u16::MAX as usize - 1) as u16;
            let (a0, a1) = Self::irange(self.x0, self.cell, self.nx, w[0][0].min(w[1][0]) - pad, w[0][0].max(w[1][0]) + pad);
            let (b0, b1) = Self::irange(self.y0, self.cell, self.ny, w[0][1].min(w[1][1]) - pad, w[0][1].max(w[1][1]) + pad);
            for j in b0.max(j0)..b1.min(j1) { for i in a0.max(i0)..a1.min(i1) {
                let q = &mut rng[(j - j0) * lw + (i - i0)]; q.0 = q.0.min(si); q.1 = q.1.max(si);
            } }
        }
        for j in j0..j1 { for i in i0..i1 {
            let q = rng[(j - j0) * lw + (i - i0)];
            if q.0 == u16::MAX { continue; }
            let (x, y) = self.xy(i, j);
            let t = path_thick(u, x, y, (q.0 as usize, q.1 as usize));
            if t > 0.0 { self.v[j * self.nx + i] += t; }
        } }
    }
    /// Separable Gaussian smoothing (edge-clamped).
    fn smooth(&self, sigma_m: f64) -> SeqField {
        let sc = sigma_m / self.cell;
        if sc <= 0.0 { return self.clone(); }
        let r = (3.0 * sc).ceil() as i64;
        let k: Vec<f64> = (-r..=r).map(|d| (-0.5 * (d as f64 / sc).powi(2)).exp()).collect();
        let ks: f64 = k.iter().sum();
        let (nx, ny) = (self.nx as i64, self.ny as i64);
        let mut t = vec![0.0f64; self.v.len()];
        for j in 0..ny { for i in 0..nx {
            let mut a = 0.0; for (q, d) in (-r..=r).enumerate() { a += k[q] * self.v[(j * nx + (i + d).clamp(0, nx - 1)) as usize]; }
            t[(j * nx + i) as usize] = a / ks;
        } }
        let mut o = vec![0.0f64; self.v.len()];
        for j in 0..ny { for i in 0..nx {
            let mut a = 0.0; for (q, d) in (-r..=r).enumerate() { a += k[q] * t[((j + d).clamp(0, ny - 1) * nx + i) as usize]; }
            o[(j * nx + i) as usize] = a / ks;
        } }
        SeqField { v: o, ..self.clone() }
    }
}

/// Lane DATUM2: sharp lens thickness of a path unit at (x, y) (the path branch of `VolcLobes::thick_at`, same arithmetic).
fn path_thick(u: &Unit, x: f64, y: f64, seg: (usize, usize)) -> f64 {
    let mut best = (f64::INFINITY, 0.0f64);
    let mut side = 0usize;
    let hi = (seg.1 + 2).min(u.path.len());
    for w in u.path[seg.0.min(hi.saturating_sub(2))..hi].windows(2) {
        let (ax, ay, bx, by) = (w[0][0], w[0][1], w[1][0], w[1][1]);
        let (vx, vy) = (bx - ax, by - ay);
        let l2 = vx * vx + vy * vy;
        let t = if l2 > 0.0 { (((x - ax) * vx + (y - ay) * vy) / l2).clamp(0.0, 1.0) } else { 0.0 };
        let dd = (x - ax - t * vx).hypot(y - ay - t * vy);
        if dd < best.0 { best = (dd, w[0][2] + t * (w[1][2] - w[0][2])); side = (vx * (y - ay) - vy * (x - ax) < 0.0) as usize; }
    }
    let (lat, s) = best;
    let w = u.hw * (1.0 + 0.5 * s / u.run.max(1.0)) * u.marg.as_ref().map(|m| marg_at(m, side, s)).unwrap_or(1.0);
    let wl = if s > u.run - w { let q = ((s - (u.run - w)) / w).min(1.0); w * (1.0 - q * q).max(0.0).sqrt() } else { w };
    if lat >= wl || wl <= 0.0 { return 0.0; }
    u.thick * (1.0 - (lat / wl).powi(2)).sqrt() * (s / (2.0 * u.hw).max(20.0)).min(1.0).max(0.3)
}

pub const STACK_NOTE: &str = "lane DATUM de_stack=1: the lobe column is stacked by levels instead of hanging every lens from the construction surface D. Each covering unit's thickness is smoothed laterally (lens cross-section convolved across its axis with a Gaussian of sigma = blur_hw_mult x the unit half-width, clamped to [blur_min_m, blur_max_m]); contact depths below D are the cumulative smoothed thicknesses (youngest first), so contacts follow D with gentle relief instead of the sum of sharp lens bases (U / V shapes with vertical limbs). A unit is present only where its own (sharp) lens thickness is > 0; the depth interval of a unit absent at the post goes to the nearest younger present unit above it (younger units fill around and onlap older lobes: compensational stacking, approximated per post), or, above the youngest present unit, to that unit (extended up to D). Volume per unit is kept laterally (the smoothing conserves the cross-section area); per post thickness is not the sharp lens. Tephra and VE cone / dome / fissure shapes are not smoothed. No new keyed draws.";
pub const ROUTE_NOTE: &str = "lane DATUM de_route=1: keyed unit paths (steepest descent) are traced on the box mean of the T1b height over (2 h + 1)^2 T1b posts (h = box_half_posts; 4 = the BS edifice datum S, 288 m), bicubic, instead of the T1b bicubic: units were emplaced before the present sub-300 m valleys (T1b's young radial gullies, relaxation) were cut, so they should not converge into them. The box mean is a proxy for the constructional surface (it keeps broader relief), not a reconstruction of it.";

#[derive(Clone, Debug, Serialize)]
pub struct StackParams {
    pub version: u32, pub blur_hw_mult: f64, pub blur_min_m: f64, pub blur_max_m: f64, pub note: String,
    /// Lane DATUM2 (opt-in de_lens=1|2): flat-based lava lenses inside the level stack. None = off (params JSON unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lens: Option<LensParams>,
}
impl Default for StackParams {
    fn default() -> Self { StackParams { version: 1, blur_hw_mult: 1.0, blur_min_m: 20.0, blur_max_m: 200.0, note: STACK_NOTE.into(), lens: None } }
}

pub const LENS_NOTE: &str = "lane DATUM2 de_lens (needs de_stack=1): lava units are flat-based lenses inside the level stack, as real lava lobes are (a flow spreads over the ground it meets: its base follows that surface, its top is domed or flat with steep margins and snout). The column is built bottom-up from the level stack base (D minus the summed smoothed thicknesses, as de_stack): each present unit starts on the top of the material below it; a lava unit's top is its level top (de_stack) raised by its own lens excess (profile thickness minus smoothed thickness: a dome over the axis, a pinch-out towards the margins), so its base is the (smooth) top of what lies below and its own relief shows at its top; pyroclastic density current / lahar deposits, tephra falls and VE cone / dome / fissure shapes fill up to their level top (sheet-like, levelling: they drape thin over older lava domes and thicken around them, and pinch out where an older dome rises above their level). A unit whose base already lies above its target top is absent at the post (onlap against an older dome). Tops above D are truncated by D (erosion); above the youngest present unit, that unit extends to D (as de_stack). Profile mode 1: the lens cross-section of the record (elliptic, tc sqrt(1 - (l / w)^2)); mode 2: a flat-topped lobe, tc sqrt(1 - (l / w)^4) (steeper margins, flatter top: a'a / pahoehoe sheet lobes). The total stack thickness is unchanged from de_stack (no new keyed draws), so broad sag where unit paths converge is not addressed here.";

#[derive(Clone, Debug, Serialize)]
pub struct LensParams { pub version: u32, pub mode: u8, pub note: String }
impl LensParams { pub fn new(mode: u8) -> Self { LensParams { version: 1, mode, note: LENS_NOTE.into() } } }
#[derive(Clone, Debug, Serialize)]
pub struct RouteParams { pub version: u32, pub box_half_posts: i64, pub note: String }
impl Default for RouteParams {
    fn default() -> Self { RouteParams { version: 1, box_half_posts: 4, note: ROUTE_NOTE.into() } }
}

#[derive(Clone, Debug, Serialize)]
pub struct LobateParams {
    pub version: u32,
    /// Max relative half-width change (sum of the three mode amplitudes is at most amp).
    pub amp: f64,
    /// Wavelength range as multiples of the unit half-width, and its clamp (m).
    pub wl_mult: [f64; 2],
    pub wl_clamp_m: [f64; 2],
    /// VE Path resampling step (m).
    pub ve_step_m: f64,
    pub note: String,
}

impl Default for LobateParams {
    fn default() -> Self { LobateParams { version: 1, amp: 0.3, wl_mult: [1.5, 4.0], wl_clamp_m: [60.0, 600.0], ve_step_m: 16.0, note: LOBATE_NOTE.into() } }
}

/// Lane RIB: per-unit margin modes, per side (left, right) three (amplitude, wavenumber rad/m, phase).
pub type Marg = [[f64; 9]; 2];

fn draw_marg(lb: &LobateParams, scope: &Scope, key: &str, hw: f64) -> R<Marg> {
    let sk = scope.address(&format!("{key}/margin")).map_err(|e| e.0)?;
    let mut n = 0u64;
    let mut u = || -> R<f64> { n += 1; Ok(sk.draw(PURPOSE_LOBATE, n, 0).map_err(|e| e.0)?.uniform_f64()) };
    let mut m = [[0.0f64; 9]; 2];
    for side in m.iter_mut() {
        let w: [f64; 3] = [u()?, u()?, u()?];
        let ws: f64 = w.iter().sum::<f64>().max(1e-9);
        for k in 0..3 {
            let lam = (hw * (lb.wl_mult[0] + (lb.wl_mult[1] - lb.wl_mult[0]) * u()?)).clamp(lb.wl_clamp_m[0], lb.wl_clamp_m[1]);
            side[3 * k] = lb.amp * w[k] / ws;
            side[3 * k + 1] = std::f64::consts::TAU / lam;
            side[3 * k + 2] = std::f64::consts::TAU * u()?;
        }
    }
    Ok(m)
}

#[inline]
fn marg_at(m: &Marg, side: usize, s: f64) -> f64 {
    let q = &m[side];
    1.0 + q[0] * (q[1] * s + q[2]).sin() + q[3] * (q[4] * s + q[5]).sin() + q[6] * (q[7] * s + q[8]).sin()
}

/// Lane RIB: Catmull-Rom resample of a polyline [x, y, c] (c interpolated linearly) at about `step` m.
fn resample(pts: &[[f64; 3]], step: f64) -> Vec<[f64; 3]> {
    if pts.len() < 3 { return pts.to_vec(); }
    let n = pts.len();
    let mut out = vec![pts[0]];
    for i in 0..n - 1 {
        let p0 = pts[if i == 0 { 0 } else { i - 1 }]; let p1 = pts[i]; let p2 = pts[i + 1]; let p3 = pts[(i + 2).min(n - 1)];
        let l = (p2[0] - p1[0]).hypot(p2[1] - p1[1]);
        let m = ((l / step).ceil() as usize).max(1);
        for j in 1..=m {
            let t = j as f64 / m as f64; let (t2, t3) = (t * t, t * t * t);
            let cr = |a: f64, b: f64, c: f64, d: f64| 0.5 * (2.0 * b + (-a + c) * t + (2.0 * a - 5.0 * b + 4.0 * c - d) * t2 + (-a + 3.0 * b - 3.0 * c + d) * t3);
            out.push([cr(p0[0], p1[0], p2[0], p3[0]), cr(p0[1], p1[1], p2[1], p3[1]), p1[2] + t * (p2[2] - p1[2])]);
        }
    }
    out
}

impl Default for LobeParams {
    fn default() -> Self { LobeParams { version: 1, cell_m: 256.0, coverage: 8.0, step_m: 20.0, stop_slope: 0.02, ve: None, lobate: None, stack: None, route: None, seq: None, tephra_rec: None, relief: None } }
}


/// Lane DATUM: reach of the smoothed lens beyond the sharp margin, in sigmas (the quadrature spans +-BLUR_REACH sigma).
const BLUR_REACH: f64 = 3.0;
/// Lane DATUM: lateral smoothing sigma of a unit of half-width hw.
fn blur_sigma(sp: &StackParams, hw: f64) -> f64 { (sp.blur_hw_mult * hw).clamp(sp.blur_min_m, sp.blur_max_m) }
/// Lane DATUM: Gaussian quadrature nodes (in sigmas) and weights (sum 1) over +-BLUR_REACH sigma.
fn blur_nodes() -> [(f64, f64); 25] {
    let mut o = [(0.0, 0.0); 25];
    let mut sw = 0.0;
    for (k, q) in o.iter_mut().enumerate() { let z = -BLUR_REACH + 2.0 * BLUR_REACH * k as f64 / 24.0; *q = (z, (-0.5 * z * z).exp()); sw += q.1; }
    for q in o.iter_mut() { q.1 /= sw; }
    o
}

/// Lane DATUM: cross-section of a lens unit at (x, y): signed lateral offset l from its axis (side 0 for l >= 0), the margin
/// half-widths per side (snout taper applied) and the centre thickness at that along-path position. Same geometry as
/// `thick_at` (lens(l) = tc sqrt(1 - (l / wl[side])^2) inside the margin reproduces it); None off the unit's along-path range.
struct LensAt { l: f64, wl: [f64; 2], tc: f64 }
impl LensAt {
    #[inline]
    fn lens(&self, l: f64) -> f64 {
        let w = self.wl[(l < 0.0) as usize];
        if w <= 0.0 || l.abs() >= w { 0.0 } else { self.tc * (1.0 - (l / w).powi(2)).sqrt() }
    }
    /// Lane DATUM2 de_lens=2: flat-topped lobe profile tc sqrt(1 - (l / w)^4) (same margins as `lens`).
    #[inline]
    fn lobe(&self, l: f64) -> f64 {
        let w = self.wl[(l < 0.0) as usize];
        if w <= 0.0 || l.abs() >= w { 0.0 } else { self.tc * (1.0 - (l / w).powi(4)).sqrt() }
    }
}

/// One unit kind: share of units, thickness (median, sigma ln, lo, hi), half-width, runout (same), and material.
#[derive(Clone, Copy, Debug)]
pub struct KindCfg { pub share: f64, pub t: [f64; 4], pub hw: [f64; 4], pub run: [f64; 4], pub mat: Mat }

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mat { Lava, Pdc, Tephra }

/// Eruption style: unit kinds (lava, pdc / lahar, tephra fall: t = thickness at the vent, hw unused, run = decay length).
#[derive(Clone, Copy, Debug)]
pub struct Style { pub name: &'static str, pub kinds: [KindCfg; 3] }

pub const STYLE_ARC: Style = Style { name: "arc_stratovolcano_vulcanian_plinian", kinds: [
    KindCfg { share: 0.45, t: [10.0, 0.6, 4.0, 30.0], hw: [60.0, 0.5, 30.0, 150.0], run: [800.0, 0.5, 300.0, 2000.0], mat: Mat::Lava },
    KindCfg { share: 0.35, t: [5.0, 0.6, 2.0, 15.0], hw: [40.0, 0.5, 20.0, 100.0], run: [2000.0, 0.5, 1000.0, 4000.0], mat: Mat::Pdc },
    KindCfg { share: 0.20, t: [3.0, 0.7, 0.5, 10.0], hw: [0.0, 0.0, 0.0, 0.0], run: [2500.0, 0.4, 1000.0, 6000.0], mat: Mat::Tephra },
] };
pub const STYLE_SHIELD: Style = Style { name: "shield_hawaiian_strombolian", kinds: [
    KindCfg { share: 0.93, t: [2.5, 0.6, 0.5, 8.0], hw: [50.0, 0.6, 20.0, 150.0], run: [2500.0, 0.5, 1000.0, 6000.0], mat: Mat::Lava },
    KindCfg { share: 0.02, t: [1.0, 0.5, 0.3, 3.0], hw: [20.0, 0.5, 10.0, 60.0], run: [1000.0, 0.5, 300.0, 2000.0], mat: Mat::Pdc },
    KindCfg { share: 0.05, t: [1.0, 0.6, 0.2, 4.0], hw: [0.0, 0.0, 0.0, 0.0], run: [1500.0, 0.4, 500.0, 4000.0], mat: Mat::Tephra },
] };
pub const STYLE_LIP: Style = Style { name: "flood_basalt_sheet_lobes", kinds: [
    KindCfg { share: 0.97, t: [15.0, 0.5, 5.0, 40.0], hw: [600.0, 0.5, 300.0, 1500.0], run: [5000.0, 0.4, 3000.0, 10000.0], mat: Mat::Lava },
    KindCfg { share: 0.0, t: [1.0, 0.5, 0.3, 3.0], hw: [20.0, 0.5, 10.0, 60.0], run: [1000.0, 0.5, 300.0, 2000.0], mat: Mat::Pdc },
    KindCfg { share: 0.03, t: [0.5, 0.6, 0.1, 2.0], hw: [0.0, 0.0, 0.0, 0.0], run: [2000.0, 0.4, 500.0, 5000.0], mat: Mat::Tephra },
] };

pub fn style_of_kind(kind: &str) -> Style {
    match kind {
        "hotspot" | "rift" => STYLE_SHIELD,
        "large_igneous_province" => STYLE_LIP,
        _ => STYLE_ARC,
    }
}

/// Lane VE's per-post edifice_style code (pf-region t1b_ve.rs STYLE_NOTE, read 2026-10-15; not yet exercised on a VE window):
/// basaltic styles (3 cinder cone, 8 shield flow, 9 fissure / spatter rampart) -> shield; 13 plateau flow -> LIP; andesitic /
/// composite styles (1 lava flow, 2 lobe front, 4 pyroclastic apron, 5 lahar, 6 dome, 7 crater wall, 11 tuff ring, 12 composite
/// body) -> arc; 0 none, 10 crater floor, 14 neck -> the T1b kind.
pub fn style_of_ve(code: u8) -> Option<Style> {
    match code { 3 | 8 | 9 => Some(STYLE_SHIELD), 13 => Some(STYLE_LIP), 1 | 2 | 4 | 5 | 6 | 7 | 11 | 12 => Some(STYLE_ARC), _ => None }
}

/// One generated (or recorded) unit.
#[derive(Clone, Debug)]
pub struct Unit {
    /// Relative age (0 youngest).
    pub age: f64,
    pub mat: Mat,
    /// Path polyline (chart m) with cumulative length; for tephra: the vent only.
    pub path: Vec<[f64; 3]>,
    pub thick: f64,
    pub hw: f64,
    pub run: f64,
    /// Lava: rubbly top and basal breccia fractions.
    pub ft: f64,
    pub fb: f64,
    /// ln m of the parts (top, core, base); pdc / tephra use [0] for the whole unit.
    pub lm: [f64; 3],
    /// Recorded T1b flow (wander amplitude, length, phase; angle; r0) when Some.
    pub t1b: Option<[f64; 8]>,
    /// Lane DE2 F4: a recorded lane-VE unit with a non-path geometry (cone, dome, fissure).
    pub veg: Option<VeShape>,
    pub id: u16,
    /// Lane RIB (de_lobate=1): lobate margin modes (None = straight margins).
    pub marg: Option<Marg>,
    /// Lane PYRO (pyro_sub=1): pyroclastic subtype of a pdc / lahar / tephra unit and its keyed strength z (None = off or lava).
    pub sub: Option<(crate::pyro::UnitSub, f64)>,
    /// Lane LOBE (de_lobe_relief=1): index of this recorded lane VE lava flow in its volcano's relief units (None = not one).
    pub rel: Option<u32>,
}

/// Lane DE2 F4: VE unit shapes other than paths (pf-region t1b_ve.rs unit_drape rules; heights already net of erosion).
#[derive(Clone, Copy, Debug)]
pub enum VeShape {
    Cone { x: f64, y: f64, rb: f64, h: f64, rc: f64, dc: f64 },
    Dome { x: f64, y: f64, r: f64, h: f64 },
    Fissure { x0: f64, y0: f64, x1: f64, y1: f64, h: f64, w: f64 },
}

impl VeShape {
    pub(crate) fn thick(&self, x: f64, y: f64) -> f64 {
        let ss = |u: f64| { let u = u.clamp(0.0, 1.0); u * u * (3.0 - 2.0 * u) };
        match *self {
            VeShape::Cone { x: cx, y: cy, rb, h, rc, dc } => {
                let r = (x - cx).hypot(y - cy); if r >= rb { return 0.0; }
                if r >= rc { h * ((rb - r) / (rb - rc).max(1e-9)).clamp(0.0, 1.0) } else { (h - dc.max(0.0) * (1.0 - (r / rc.max(1e-9)).powi(2))).max(0.0) }
            }
            VeShape::Dome { x: cx, y: cy, r, h } => { let d = (x - cx).hypot(y - cy); if d >= r { 0.0 } else { h * (1.0 - (d / r).powi(4)).max(0.0).sqrt() } }
            VeShape::Fissure { x0, y0, x1, y1, h, w } => {
                let (dx, dy) = (x1 - x0, y1 - y0); let l = dx.hypot(dy).max(1.0);
                let al = ((x - x0) * dx + (y - y0) * dy) / l; if al < 0.0 || al > l { return 0.0; }
                let la = ((x - x0) * dy - (y - y0) * dx).abs() / l;
                h * (-(la / w.max(1e-9)).powi(2)).exp() * ss(al.min(l - al) / 200.0)
            }
        }
    }
    pub(crate) fn bbox(&self) -> (f64, f64, f64, f64) {
        match *self {
            VeShape::Cone { x, y, rb, .. } => (x - rb, y - rb, x + rb, y + rb),
            VeShape::Dome { x, y, r, .. } => (x - r, y - r, x + r, y + r),
            VeShape::Fissure { x0, y0, x1, y1, w, .. } => (x0.min(x1) - 3.0 * w, y0.min(y1) - 3.0 * w, x0.max(x1) + 3.0 * w, y0.max(y1) + 3.0 * w),
        }
    }
}

/// Lane DE2 F4: lane VE's units of one volcano as recorded units (ages scaled into [0, 0.02): younger than T1b's cf flows,
/// which sit in the cone VE builds on). Path -> lens path unit (half-width w0_m, thickness t_m - erosion_m, lahar = weak);
/// Cone / Dome / Fissure -> their VE shapes; other kinds skipped. Returns the units and the count skipped.
/// Lane LOBE: `rel` (de_lobe_relief=1) maps a listed unit index to its relief index (lava paths): such a unit keeps the record's
/// 128 m polyline and no lobate margins (its column thickness is the relief section, `relief::section`).
pub fn ve_units(units: &serde_json::Value, scope: &Scope, tag: &str, lob: Option<&LobateParams>, rel: Option<&[Option<usize>]>) -> R<(Vec<Unit>, usize)> {
    let a = units.as_array().cloned().unwrap_or_default();
    let amax = a.iter().filter_map(|u| u.get("age_yr").and_then(|v| v.as_f64())).fold(1.0f64, f64::max);
    let (mut out, mut skipped) = (Vec::new(), 0usize);
    for (i, u) in a.iter().enumerate() {
        let f = |k: &str| u.get(k).and_then(|v| v.as_f64()).unwrap_or(0.0);
        let g = match u.get("geom") { Some(g) => g, None => { skipped += 1; continue } };
        let gf = |k: &str| g.get(k).and_then(|v| v.as_f64()).unwrap_or(0.0);
        let ero = f("erosion_m");
        let age = 0.02 * (f("age_yr") / amax).clamp(0.0, 1.0);
        let mut dr = Dr { sk: scope.address(&format!("{tag}/ve/{i}")).map_err(|e| e.0)?, n: 0 };
        let id = uid(&format!("{tag}/ve/{i}"));
        let kind = g.get("g").and_then(|v| v.as_str()).unwrap_or("");
        let shape = match kind {
            "Path" => {
                let pts: Vec<[f64; 3]> = g.get("pts").and_then(|v| v.as_array()).map(|p| p.iter().filter_map(|q| { let q = q.as_array()?; Some([q.first()?.as_f64()?, q.get(1)?.as_f64()?, q.get(2)?.as_f64()?]) }).collect()).unwrap_or_default();
                let t = gf("t_m") - ero;
                if pts.len() < 2 || t <= 0.05 { skipped += 1; continue; }
                let lahar = g.get("lahar").and_then(|v| v.as_bool()).unwrap_or(false);
                let mat = if lahar { Mat::Pdc } else { Mat::Lava };
                let (ft, fb, lm) = unit_parts(&mut dr, mat)?;
                if let Some(ri) = rel.and_then(|r| r.get(i).copied().flatten()) {
                    out.push(Unit { age, mat, path: pts, thick: t, hw: gf("w0_m"), run: gf("len_m").max(1.0), ft, fb, lm, t1b: None, veg: None, id, marg: None, sub: None, rel: Some(ri as u32) });
                    continue;
                }
                let (pts, marg) = match lob { Some(lb) => (resample(&pts, lb.ve_step_m), Some(draw_marg(lb, scope, &format!("{tag}/ve/{i}"), gf("w0_m").max(1.0))?)), None => (pts, None) };
                out.push(Unit { age, mat, path: pts, thick: t, hw: gf("w0_m"), run: gf("len_m").max(1.0), ft, fb, lm, t1b: None, veg: None, id, marg, sub: None, rel: None });
                continue;
            }
            "Cone" => VeShape::Cone { x: gf("x"), y: gf("y"), rb: gf("rb_m"), h: gf("h_m") - ero, rc: gf("rc_m"), dc: gf("dc_m") - ero },
            "Dome" => VeShape::Dome { x: gf("x"), y: gf("y"), r: gf("r_m"), h: gf("h_m") - ero },
            "Fissure" => VeShape::Fissure { x0: gf("x0"), y0: gf("y0"), x1: gf("x1"), y1: gf("y1"), h: gf("h_m") - ero, w: gf("w_m") },
            _ => { skipped += 1; continue }
        };
        // cinder cones / fissure ramparts: scoria (weak, pyroclastic); domes: massive lava
        let mat = if matches!(shape, VeShape::Dome { .. }) { Mat::Lava } else { Mat::Tephra };
        let (ft, fb, lm) = unit_parts(&mut dr, mat)?;
        out.push(Unit { age, mat, path: vec![], thick: 0.0, hw: 0.0, run: 1.0, ft, fb, lm, t1b: None, veg: Some(shape), id, marg: None, sub: None, rel: None });
    }
    Ok((out, skipped))
}

#[derive(Clone, Copy, Debug)]
pub struct Part { pub d_top: f64, pub d_base: f64, pub sc: f64, pub fac: u8, pub unit: u16,
    /// Lane PYRO: index of the unit in `VolcLobes::units` and its local deposited thickness (m; output only, no arithmetic).
    pub ui: u32, pub tu: f64 }

/// The local column at a post: parts top-down (depth below D), D, and the cumulative 1/Sc integral at part bases.
#[derive(Clone, Debug, Default)]
pub struct LobeCol { pub d_surf: f64, pub parts: Vec<Part>, pub cum: Vec<f64>, pub stack_m: f64 }

impl LobeCol {
    fn part_at(&self, d: f64) -> Option<usize> {
        if self.parts.is_empty() || d < 0.0 || d >= self.stack_m { return None; }
        let i = self.parts.partition_point(|p| p.d_base <= d);
        Some(i.min(self.parts.len() - 1))
    }
    /// Integral of 1/Sc from depth 0 to d (d < 0: extrapolated with the top part; beyond the stack: `below` integral added).
    fn integ(&self, d: f64, below: &dyn Fn(f64) -> f64) -> f64 {
        if self.parts.is_empty() { return below(d); }
        if d <= 0.0 { return d / self.parts[0].sc; }
        if d >= self.stack_m { return self.cum[self.parts.len() - 1] + below(d - self.stack_m) - below(0.0); }
        let i = self.part_at(d).unwrap();
        let prev = if i == 0 { 0.0 } else { self.cum[i - 1] };
        prev + (d - self.parts[i].d_top) / self.parts[i].sc
    }
}

/// All units of one volcano relevant to a block, with a 32 m segment index.
pub struct VolcLobes {
    pub style: &'static str,
    pub units: Vec<Unit>,
    pub tephra: Vec<usize>,
    /// 32 m cell -> (unit, first segment, last segment) of the path segments within reach of the cell.
    grid: HashMap<(i64, i64), Vec<(u32, u16, u16)>>,
    /// Normalisation of the strength multipliers (expected thickness-weighted mean of 1/m = 1).
    pub norm: f64,
    pub recorded: usize,
    vent: [f64; 2],
    /// Lane DE2 F4: VE units skipped (geometry kinds not represented).
    pub ve_skipped: usize,
    /// Lane DATUM: level stacking (None = top-down stacking of the sharp lenses).
    stack: Option<StackParams>,
    /// Lane DATUM2 de_seq=1: the smoothed built thickness Pbar (stack base depth below D) of the sequential pass; None = off.
    pub pbar: Option<SeqField>,
    /// Lane DATUM2: keyed unit sources of the sequential pass, units emplaced, their mean runout over the drawn runout, and
    /// the share that ponded and overtopped at least once (zeros when off).
    pub seq_units: usize,
    pub seq_stats: [f64; 3],
    /// Lane PYRO (pyro_sub=1): the volcano's pyroclastic record and the params (None = off).
    pub pyro: Option<(crate::pyro::PyroRec, crate::pyro::PyroParams)>,
    /// Lane PYRO2 de_seq_fill (None = off).
    pub seq_fill: Option<SeqFillParams>,
    /// Lane ARCH de_seq_level: level-stack parameters of the levelled de_seq column (None = off).
    seq_level: Option<(StackParams, [f64; 2], bool)>,
    /// Lane LOBE de_lobe_relief: the relief parameters and this volcano's lane VE units as relief units (None = off).
    pub relief: Option<(crate::relief::ReliefParams, Vec<crate::relief::RUnit>)>,
}

const GCELL: f64 = 32.0;

/// Keyed draw stream.
struct Dr { sk: Scope, n: u64 }
impl Dr {
    fn z(&mut self) -> R<f64> { self.n += 1; Ok(self.sk.draw(PURPOSE_LOBES, self.n, 0).map_err(|e| e.0)?.normal_f64()) }
    fn u(&mut self) -> R<f64> { self.n += 1; Ok(self.sk.draw(PURPOSE_LOBES, self.n, 0).map_err(|e| e.0)?.uniform_f64()) }
    fn ln(&mut self, c: [f64; 4]) -> R<f64> { Ok((c[0] * (c[1] * self.z()?).exp()).clamp(c[2], c[3])) }
}

fn unit_parts(dr: &mut Dr, mat: Mat) -> R<(f64, f64, [f64; 3])> {
    Ok(match mat {
        Mat::Lava => {
            // BS LAVA_ARC fractions and ln m (rubbly / vesicular top, dense core, basal breccia)
            let (ft, fb) = (0.15 + 0.25 * dr.u()?, 0.05 + 0.10 * dr.u()?);
            (ft, fb, [-0.6 + 0.25 * dr.z()?, 0.35 + 0.25 * dr.z()?, -0.35 + 0.2 * dr.z()?])
        }
        Mat::Pdc => { let hard = dr.u()? < 0.15; (0.0, 0.0, [if hard { 0.1 + 0.2 * dr.z()? } else { -0.55 + 0.4 * dr.z()? }, 0.0, 0.0]) }
        Mat::Tephra => (0.0, 0.0, [-0.7 + 0.3 * dr.z()?, 0.0, 0.0]),
    })
}

/// Expected thickness-weighted mean of 1/m of a style (keyed synthetic sample, fixed size; block-independent).
fn style_norm(st: &Style, scope: &Scope, tag: &str) -> R<f64> {
    let mut dr = Dr { sk: scope.address(&format!("{tag}/norm")).map_err(|e| e.0)?, n: 0 };
    let (mut tw, mut iw) = (0.0f64, 0.0f64);
    for _ in 0..4000 {
        let k = pick_kind(st, dr.u()?);
        let kc = st.kinds[k];
        let t = dr.ln(kc.t)?;
        let (ft, fb, lm) = unit_parts(&mut dr, kc.mat)?;
        if kc.mat == Mat::Lava {
            for (f, l) in [(ft, lm[0]), (1.0 - ft - fb, lm[1]), (fb, lm[2])] { tw += t * f; iw += t * f / l.exp(); }
        } else { tw += t; iw += t / lm[0].exp(); }
    }
    Ok(iw / tw)
}

fn pick_kind(st: &Style, u: f64) -> usize {
    let tot: f64 = st.kinds.iter().map(|k| k.share).sum();
    let mut acc = 0.0;
    for (i, k) in st.kinds.iter().enumerate() { acc += k.share / tot; if u < acc { return i; } }
    2
}

/// 15-bit unit id from the unit's key (0x7FFF and above reserved).
fn uid(key: &str) -> u16 { let h = blake3::hash(key.as_bytes()); (u16::from_le_bytes([h.as_bytes()[0], h.as_bytes()[1]]) & 0x7FFF).min(0x7FFE) }

#[inline]
fn ang_diff(a: f64, b: f64) -> f64 { let d = (a - b).rem_euclid(std::f64::consts::TAU); if d > std::f64::consts::PI { d - std::f64::consts::TAU } else { d } }

/// Lane DATUM2: units younger than this relative age are recorded units (T1b cf flows, VE units); keyed units and tephra
/// falls are drawn in [0.05, 1).
const RECORDED_AGE: f64 = 0.05;

impl VolcLobes {
    /// Units of volcano `vi` whose sources lie within the box (x0, y0, x1, y1) grown by the largest reach. `on_ed(x, y)` says
    /// whether a chart point is edifice ground of this volcano (T1b fields); `zf(x, y)` is the T1b bicubic.
    #[allow(clippy::too_many_arguments)]
    pub fn build(lp: &LobeParams, style: Style, vid: &str, vent: [f64; 2], flows: &[[f64; 8]], bbox: (f64, f64, f64, f64), scope: &Scope,
                 on_ed: &dyn Fn(f64, f64) -> bool, zf: &dyn Fn(f64, f64) -> f64, style_at: &dyn Fn(f64, f64) -> Option<Style>,
                 ve: Option<&serde_json::Value>, dom: Option<(f64, f64, f64, f64)>, pyro: Option<(&crate::pyro::PyroRec, &crate::pyro::PyroParams)>, relief_rbody: f64) -> R<VolcLobes> {
        let tag = format!("lobes/{vid}");
        // lane DATUM2 de_seq: sources over the whole window domain (block-independent), paths traced afterwards in age order
        let seq = lp.seq.as_ref().zip(dom);
        struct Seed { sx: f64, sy: f64, age: f64, mat: Mat, thick: f64, hw: f64, run: f64, ft: f64, fb: f64, lm: [f64; 3], key: String }
        let mut seeds: Vec<Seed> = Vec::new();
        let norm = style_norm(&style, scope, &tag)?;
        let mut units: Vec<Unit> = Vec::new();
        // a) recorded T1b flows: youngest units (ages in [0, 0.05)), T1b geometry
        for (k, f) in flows.iter().enumerate() {
            let mut dr = Dr { sk: scope.address(&format!("{tag}/t1b/{k}")).map_err(|e| e.0)?, n: 0 };
            let u0 = dr.u()?;
            // with lane VE's units present (DE2 F4) the cf flows are older than them (VE builds on the CF cone)
            let age = if ve.is_some() { 0.02 + 0.03 * u0 } else { 0.05 * u0 };
            let (ft, fb, lm) = unit_parts(&mut dr, Mat::Lava)?;
            let marg = match lp.lobate.as_ref() { Some(lb) => Some(draw_marg(lb, scope, &format!("{tag}/t1b/{k}"), f[3].max(1.0))?), None => None };
            units.push(Unit { age, mat: Mat::Lava, path: vec![], thick: f[4], hw: f[3], run: f[2], ft, fb, lm, t1b: Some(*f), veg: None, id: uid(&format!("{tag}/t1b/{k}")), marg, sub: None, rel: None });
        }
        let mut ve_skipped = 0usize;
        // lane LOBE de_lobe_relief: lane VE's units as relief units (lava paths map to their column units)
        let relief: Option<(crate::relief::ReliefParams, Vec<crate::relief::RUnit>, Vec<Option<usize>>)> = match (lp.relief.as_ref(), ve) {
            (Some(rp), Some(vu)) => {
                let (mut ru, ix) = crate::relief::ve_relief_units(vu, style.name == STYLE_ARC.name, vent, relief_rbody);
                // the regional gradient at the path vertices: the routing surface zf (box mean under de_seq / de_route), 32 m
                // central differences
                for u in ru.iter_mut() { if let crate::relief::RGeom::Path { pts, .. } = &u.geom {
                    u.bgrad = pts.iter().map(|p| { let h = 32.0; [(zf(p[0] + h, p[1]) - zf(p[0] - h, p[1])) / (2.0 * h), (zf(p[0], p[1] + h) - zf(p[0], p[1] - h)) / (2.0 * h)] }).collect();
                } }
                Some((rp.clone(), ru, ix))
            }
            _ => None,
        };
        if let Some(vu) = ve { let (mut a, sk) = ve_units(vu, scope, &tag, lp.lobate.as_ref(), relief.as_ref().map(|r| r.2.as_slice()))?; ve_skipped = sk; units.append(&mut a); }
        let recorded = units.len();
        // b) keyed units: sources per lattice cell (density per style for the target coverage)
        let reach = style.kinds.iter().map(|k| k.run[3] + 1.5 * k.hw[3]).fold(0.0f64, f64::max);
        let (x0, y0, x1, y1) = match seq { Some((_, d)) => d, None => (bbox.0 - reach, bbox.1 - reach, bbox.2 + reach, bbox.3 + reach) };
        let c = lp.cell_m;
        // expected footprint area per unit (m2) of the lobe kinds -> sources per cell
        let mut ea = 0.0;
        let tot: f64 = style.kinds.iter().map(|k| k.share).sum();
        for k in style.kinds.iter() { if k.mat != Mat::Tephra { ea += k.share / tot * k.run[0] * 2.0 * k.hw[0] * 1.25; } }
        let lam = (lp.coverage * c * c / ea.max(1.0)).min(64.0);
        for cy in (y0 / c).floor() as i64..=(y1 / c).floor() as i64 {
            for cx in (x0 / c).floor() as i64..=(x1 / c).floor() as i64 {
                let mut dr = Dr { sk: scope.address(&format!("{tag}/cell/{cx}/{cy}")).map_err(|e| e.0)?, n: 0 };
                // Poisson count by inversion (keyed uniform)
                let (u0, mut p, mut cdf, mut cnt) = (dr.u()?, (-lam).exp(), (-lam).exp(), 0usize);
                while u0 > cdf && cnt < 256 { cnt += 1; p *= lam / cnt as f64; cdf += p; }
                for si in 0..cnt {
                    let (sx, sy) = ((cx as f64 + dr.u()?) * c, (cy as f64 + dr.u()?) * c);
                    let age = 0.05 + 0.95 * dr.u()?;
                    let ku = dr.u()?;
                    let (tz, hz, rz) = (dr.z()?, dr.z()?, dr.z()?);
                    let parts = unit_parts(&mut dr, Mat::Lava)?;
                    let parts_pdc = unit_parts(&mut dr, Mat::Pdc)?;
                    if !on_ed(sx, sy) { continue; }
                    let st = style_at(sx, sy).unwrap_or(style);
                    let kc = st.kinds[pick_kind(&st, ku)];
                    if kc.mat == Mat::Tephra { continue; } // tephra falls are volcano-wide events (below)
                    let ln = |cfg: [f64; 4], z: f64| (cfg[0] * (cfg[1] * z).exp()).clamp(cfg[2], cfg[3]);
                    let (thick, hw, run) = (ln(kc.t, tz), ln(kc.hw, hz), ln(kc.run, rz));
                    let (ft, fb, lm) = if kc.mat == Mat::Lava { parts } else { parts_pdc };
                    if seq.is_some() { seeds.push(Seed { sx, sy, age, mat: kc.mat, thick, hw, run, ft, fb, lm, key: format!("{tag}/cell/{cx}/{cy}/{si}") }); continue; }
                    // path: steepest descent on the T1b bicubic (central differences, 4 m), step_m steps
                    let mut path = vec![[sx, sy, 0.0]];
                    let (mut x, mut y, mut s) = (sx, sy, 0.0);
                    while s < run {
                        let h = 4.0;
                        let (gx, gy) = ((zf(x + h, y) - zf(x - h, y)) / (2.0 * h), (zf(x, y + h) - zf(x, y - h)) / (2.0 * h));
                        let g = gx.hypot(gy);
                        if g < lp.stop_slope { break; }
                        let st_ = lp.step_m.min(run - s);
                        x -= gx / g * st_; y -= gy / g * st_; s += st_;
                        path.push([x, y, s]);
                    }
                    if path.len() < 2 { continue; }
                    let run = s;
                    let marg = match lp.lobate.as_ref() { Some(lb) => Some(draw_marg(lb, scope, &format!("{tag}/cell/{cx}/{cy}/{si}"), hw)?), None => None };
                    units.push(Unit { age, mat: kc.mat, path, thick, hw, run, ft, fb, lm, t1b: None, veg: None, id: uid(&format!("{tag}/cell/{cx}/{cy}/{si}")), marg, sub: None, rel: None });
                }
            }
        }
        // lane DATUM2 de_seq: emplace the keyed units oldest first on C + built (SEQ_NOTE)
        let mut built: Option<SeqField> = None;
        let seq_units = seeds.len();
        let (mut n_emp, mut run_ratio, mut n_pond) = (0usize, 0.0f64, 0usize);
        if let Some((sp, d)) = seq {
            let mut f = SeqField::new(d, sp.cell_m);
            let amp = lp.lobate.as_ref().map(|l| l.amp).unwrap_or(0.0);
            seeds.sort_by(|a, b| b.age.partial_cmp(&a.age).unwrap().then(a.key.cmp(&b.key)));
            for sd in seeds {
                let zr = |x: f64, y: f64| zf(x, y) + f.bilin(x, y);
                let mut path = vec![[sd.sx, sd.sy, 0.0]];
                let (mut x, mut y, mut s) = (sd.sx, sd.sy, 0.0);
                let mut zc = zr(x, y);
                let mut over = false;
                let h = sp.cell_m;
                while s < sd.run {
                    let (gx, gy) = ((zr(x + h, y) - zr(x - h, y)) / (2.0 * h), (zr(x, y + h) - zr(x, y - h)) / (2.0 * h));
                    let g = gx.hypot(gy);
                    if g < lp.stop_slope { break; }
                    // a step that would climb (overshooting a channel axis) is halved up to three times; a unit stops only
                    // where every trial climbs (ponded in a closed low)
                    let mut st_ = lp.step_m.min(sd.run - s);
                    let mut nxt = None;
                    for _ in 0..4 {
                        let (nx_, ny_) = (x - gx / g * st_, y - gy / g * st_);
                        let zn = zr(nx_, ny_);
                        if zn <= zc { nxt = Some((nx_, ny_, zn)); break; }
                        st_ *= 0.5;
                    }
                    let (nx_, ny_, zn) = match nxt {
                        Some(q) => q,
                        None => {
                            // ponded behind older units: the flow inflates and overtops, i.e. it continues down the
                            // routing base C for this step (its drawn runout is kept)
                            if !over { n_pond += 1; over = true; }
                            if sp.pond_stop == Some(true) { break; }
                            st_ = lp.step_m.min(sd.run - s);
                            let (cx_, cy_) = ((zf(x + h, y) - zf(x - h, y)) / (2.0 * h), (zf(x, y + h) - zf(x, y - h)) / (2.0 * h));
                            let gc = cx_.hypot(cy_);
                            if gc < lp.stop_slope { break; }
                            let (nx_, ny_) = (x - cx_ / gc * st_, y - cy_ / gc * st_);
                            (nx_, ny_, zr(nx_, ny_))
                        }
                    };
                    x = nx_; y = ny_; s += st_; zc = zn;
                    path.push([x, y, s]);
                }
                if path.len() < 2 { continue; }
                let marg = match lp.lobate.as_ref() { Some(lb) => Some(draw_marg(lb, scope, &sd.key, sd.hw)?), None => None };
                n_emp += 1; run_ratio += s / sd.run.max(1e-9);
                let u = Unit { age: sd.age, mat: sd.mat, path, thick: sd.thick, hw: sd.hw, run: s, ft: sd.ft, fb: sd.fb, lm: sd.lm, t1b: None, veg: None, id: uid(&sd.key), marg, sub: None, rel: None };
                f.add_path(&u, amp);
                units.push(u);
            }
            built = Some(f);
        }
        // tephra falls: volcano-wide events keyed by index (count from the style share and the coverage)
        let mut tephra = Vec::new();
        let kt = style.kinds[2];
        // lane PYRO2 de_tephra_rec: one mantle per recorded explosive eruption, VE's fall size for its VEI (TEPHRA_REC_NOTE)
        let rec_ev: Option<Vec<u8>> = match (lp.tephra_rec.as_ref(), pyro) { (Some(_), Some((r, _))) => Some(r.explosive_events()), _ => None };
        let mut forced: Vec<(usize, u8)> = Vec::new();
        if kt.share > 0.0 || rec_ev.is_some() {
            let n_t = (lp.coverage * kt.share / tot.max(1e-9) / (1.0 - kt.share / tot).max(0.05)).round() as usize;
            let n_use = match &rec_ev { Some(v) => v.len(), None => n_t.max(1) };
            for i in 0..n_use {
                let mut dr = Dr { sk: scope.address(&format!("{tag}/tephra/{i}")).map_err(|e| e.0)?, n: 0 };
                let age = 0.05 + 0.95 * dr.u()?;
                let (mut t0, mut lamd) = (dr.ln(kt.t)?, dr.ln(kt.run)?);
                if let (Some(v), Some((_, pp))) = (&rec_ev, pyro) { let q = crate::pyro::ve_fall(v[i], pp); t0 = q.0; lamd = q.1; forced.push((units.len(), v[i])); }
                let wind = match pyro { Some((r, pp)) => crate::pyro::wind_axis(r, pp, dr.u()?), None => dr.u()? * std::f64::consts::TAU }; // lane VH-I pyro_wind (same draw)
                let (_, _, lm) = unit_parts(&mut dr, Mat::Tephra)?;
                tephra.push(units.len());
                units.push(Unit { age, mat: Mat::Tephra, path: vec![[wind, 0.0, 0.0]], thick: t0, hw: 0.0, run: lamd, ft: 0.0, fb: 0.0, lm, t1b: None, veg: None, id: uid(&format!("{tag}/tephra/{i}")), marg: None, sub: None, rel: None });
            }
        }
        // unit ids: a 15-bit hash of the unit's key (block-independent, F9); ids are labels, not ages
        // segment index (32 m cells) of the path units; recorded flows indexed by their bounding sector box
        let mut grid: HashMap<(i64, i64), Vec<(u32, u16, u16)>> = HashMap::new();
        for (ui, u) in units.iter().enumerate() {
            if let Some(sh) = u.veg {
                let (a0, b0, a1, b1) = sh.bbox();
                if a1 < bbox.0 || a0 > bbox.2 || b1 < bbox.1 || b0 > bbox.3 { continue; }
                for gy in (b0.max(bbox.1) / GCELL).floor() as i64..=(b1.min(bbox.3) / GCELL).floor() as i64 {
                    for gx in (a0.max(bbox.0) / GCELL).floor() as i64..=(a1.min(bbox.2) / GCELL).floor() as i64 { grid.entry((gx, gy)).or_default().push((ui as u32, 0, 0)); }
                }
                continue;
            }
            if u.mat == Mat::Tephra { continue; }
            let pts: Vec<[f64; 2]> = match u.t1b {
                Some(f) => { let m = 24; (0..=m).map(|i| { let r = f[1] + f[2] * i as f64 / m as f64; let off = f[5] * (std::f64::consts::TAU * r / f[6].max(1.0) + f[7]).sin();
                    let (ca, sa) = (f[0].cos(), f[0].sin()); [vent[0] + r * ca - off * sa, vent[1] + r * sa + off * ca] }).collect() }
                None => u.path.iter().map(|p| [p[0], p[1]]).collect(),
            };
            let pad = match lp.lobate.as_ref() { None => 1.5 * u.hw, Some(lb) => 1.5 * u.hw * (1.0 + lb.amp) } + (if u.t1b.is_some() { u.t1b.unwrap()[5].abs() } else { 0.0 }) + 0.75 * GCELL;
            // lane LOBE: a relief unit reaches T1b's snout disk (<= 1.875 w0 past the end) plus 4 sigma of its degraded margin
            let pad = match (u.rel, relief.as_ref()) {
                (Some(ri), Some((rp, ru, _))) => 2.0 * u.hw + 64.0 + 4.0 * crate::relief::sigma_of(rp, &ru[ri as usize]) + 0.75 * GCELL,
                _ => pad,
            };
            // lane DATUM de_stack: the smoothed lens reaches BLUR_REACH sigma beyond the sharp margin
            let pad = match lp.stack.as_ref() { None => pad, Some(sp) => pad + BLUR_REACH * blur_sigma(sp, u.hw) };
            // lane ARCH de_seq_level: the levelled column also reads the smoothed lens (same reach); off: pad unchanged
            let pad = match lp.seq.as_ref().filter(|q| q.mode == 1).and_then(|q| q.level.as_ref()) {
                None => pad,
                Some(lv) => { let ex = BLUR_REACH * blur_sigma(&lv.stack(), u.hw); let st = lp.stack.as_ref().map(|sp| BLUR_REACH * blur_sigma(sp, u.hw)).unwrap_or(0.0); pad + (ex - st).max(0.0) }
            };
            for (si, w) in pts.windows(2).enumerate() {
                let (mnx, mxx) = (w[0][0].min(w[1][0]) - pad, w[0][0].max(w[1][0]) + pad);
                let (mny, mxy) = (w[0][1].min(w[1][1]) - pad, w[0][1].max(w[1][1]) + pad);
                if mxx < bbox.0 || mnx > bbox.2 || mxy < bbox.1 || mny > bbox.3 { continue; }
                let si = si.min(u16::MAX as usize - 1) as u16;
                for gy in (mny.max(bbox.1) / GCELL).floor() as i64..=(mxy.min(bbox.3) / GCELL).floor() as i64 {
                    for gx in (mnx.max(bbox.0) / GCELL).floor() as i64..=(mxx.min(bbox.2) / GCELL).floor() as i64 {
                        let e = grid.entry((gx, gy)).or_default();
                        match e.last_mut() { Some(q) if q.0 == ui as u32 => { q.1 = q.1.min(si); q.2 = q.2.max(si); } _ => e.push((ui as u32, si, si)) }
                    }
                }
            }
        }
        for v in grid.values_mut() { v.sort_unstable(); }
        // lane PYRO pyro_sub: subtype of every pdc / lahar / tephra unit (separate keyed stream; nothing above changes). The draw key
        // is the unit's id and age bits (both pure functions of the unit's own key, so block-independent).
        if let Some((rec, pp)) = pyro {
            for (ui, u) in units.iter_mut().enumerate() {
                if u.mat == Mat::Lava { continue; }
                let key = format!("{tag}/u/{}/{:016x}", u.id, u.age.to_bits());
                // lane PYRO2: a recorded eruption's mantle takes its recorded VEI (no class draw)
                if let Some(&(_, vei)) = forced.iter().find(|q| q.0 == ui) {
                    let z = crate::pyro::draw_sub(rec, pp, scope, &key, crate::pyro::Pool::Fall)?.1;
                    u.sub = Some((crate::pyro::UnitSub::Fall { lap: rec.lap_of(vei), r_lap: rec.r_lap(pp, vei), wind: u.path[0][0] }, z));
                    continue;
                }
                u.sub = Some(match (u.mat, u.veg.is_some(), u.age < RECORDED_AGE) {
                    // recorded lane VE units keep their kind: cinder cones / fissure ramparts are scoria, lahar paths lahar
                    (Mat::Tephra, true, _) => (crate::pyro::UnitSub::Scoria, crate::pyro::draw_sub(rec, pp, scope, &key, crate::pyro::Pool::Fall)?.1),
                    (Mat::Pdc, _, true) => (crate::pyro::UnitSub::Lahar, crate::pyro::draw_sub(rec, pp, scope, &key, crate::pyro::Pool::PdcLahar)?.1),
                    // a tephra mantle's dispersal axis is its own keyed wind direction
                    (Mat::Tephra, false, _) => { let (s0, z) = crate::pyro::draw_sub(rec, pp, scope, &key, crate::pyro::Pool::Fall)?;
                        (match s0 { crate::pyro::UnitSub::Fall { lap, r_lap, .. } => crate::pyro::UnitSub::Fall { lap, r_lap, wind: u.path[0][0] }, o => o }, z) }
                    _ => crate::pyro::draw_sub(rec, pp, scope, &key, crate::pyro::Pool::PdcLahar)?,
                });
            }
        }
        let mut vl = VolcLobes { style: style.name, units, tephra, grid, norm, recorded, vent, ve_skipped, stack: lp.stack.clone(), pbar: None, seq_units: 0, seq_stats: [0.0; 3],
                                 pyro: pyro.map(|(r, p)| (r.clone(), p.clone())), seq_fill: lp.seq.as_ref().and_then(|q| q.fill.clone()),
                                 seq_level: lp.seq.as_ref().filter(|q| q.mode == 1).and_then(|q| q.level.as_ref()).map(|l| (l.stack(), [l.margin_deg[0].to_radians().tan(), l.margin_deg[1].to_radians().tan()], l.stretch.is_some())),
                                 relief: relief.map(|r| (r.0, r.1)) };
        if let (Some((sp, _)), Some(mut f)) = (seq, built) {
            vl.seq_units = seq_units;
            vl.seq_stats = [n_emp as f64, run_ratio / (n_emp.max(1) as f64), n_pond as f64 / (seq_units.max(1) as f64)];
            if sp.mode == 1 {
                // the rest of the pile on the lattice: recorded T1b flows, VE units and tephra falls (exact thicknesses)
                for u in vl.units.iter() {
                    if u.t1b.is_none() && u.veg.is_none() && u.mat != Mat::Tephra { continue; }
                    let (a0, b0, a1, b1) = match u.veg { Some(sh) => sh.bbox(), None => (f.x0, f.y0, f.x0 + f.cell * (f.nx - 1) as f64, f.y0 + f.cell * (f.ny - 1) as f64) };
                    let (i0, i1) = SeqField::irange(f.x0, f.cell, f.nx, a0, a1);
                    let (j0, j1) = SeqField::irange(f.y0, f.cell, f.ny, b0, b1);
                    for j in j0..j1 { for i in i0..i1 { let (x, y) = f.xy(i, j); let t = vl.thick_at(u, x, y, (0, 0)); if t > 0.0 { f.v[j * f.nx + i] += t; } } }
                }
                vl.pbar = Some(f.smooth(sp.base_sigma_m));
            }
        }
        Ok(vl)
    }

    /// Lane DATUM: lens cross-section geometry of a path / recorded unit at (x, y) (None for tephra, VE shapes, or off its range).
    fn lens_at(&self, u: &Unit, x: f64, y: f64, seg: (usize, usize)) -> Option<LensAt> {
        if u.veg.is_some() || u.mat == Mat::Tephra || u.rel.is_some() { return None; }
        let taper = |w: f64, s: f64, run: f64| -> f64 { if s > run - w { let q = ((s - (run - w)) / w).min(1.0); w * (1.0 - q * q).max(0.0).sqrt() } else { w } };
        match u.t1b {
            Some(f) => {
                let (dx, dy) = (x - self.vent[0], y - self.vent[1]);
                let r = dx.hypot(dy);
                let d = ang_diff(dy.atan2(dx), f[0]);
                if d.abs() >= std::f64::consts::FRAC_PI_2 { return None; }
                let off = f[5] * (std::f64::consts::TAU * r / f[6].max(1.0) + f[7]).sin();
                let lat_s = r * d.sin() - off;
                let s = r * d.cos() - f[1];
                if s < 0.0 || s > f[2] { return None; }
                let w = |side: usize| taper(f[3] * (1.0 + 0.5 * s / f[2]) * u.marg.as_ref().map(|m| marg_at(m, side, s)).unwrap_or(1.0), s, f[2]);
                // thick_at's side is (lat_s < 0); the signed offset keeps that convention (l < 0 <-> side 1)
                Some(LensAt { l: lat_s, wl: [w(0), w(1)], tc: f[4] * (s / 200.0).min(1.0) })
            }
            None => {
                let mut best = (f64::INFINITY, 0.0f64);
                let mut side = 0usize;
                let hi = (seg.1 + 2).min(u.path.len());
                for w in u.path[seg.0.min(hi.saturating_sub(2))..hi].windows(2) {
                    let (ax, ay, bx, by) = (w[0][0], w[0][1], w[1][0], w[1][1]);
                    let (vx, vy) = (bx - ax, by - ay);
                    let l2 = vx * vx + vy * vy;
                    let t = if l2 > 0.0 { (((x - ax) * vx + (y - ay) * vy) / l2).clamp(0.0, 1.0) } else { 0.0 };
                    let dd = (x - ax - t * vx).hypot(y - ay - t * vy);
                    if dd < best.0 { best = (dd, w[0][2] + t * (w[1][2] - w[0][2])); side = (vx * (y - ay) - vy * (x - ax) < 0.0) as usize; }
                }
                if !best.0.is_finite() { return None; }
                let (lat, s) = best;
                let w = |sd: usize| taper(u.hw * (1.0 + 0.5 * s / u.run.max(1.0)) * u.marg.as_ref().map(|m| marg_at(m, sd, s)).unwrap_or(1.0), s, u.run);
                Some(LensAt { l: if side == 1 { -lat } else { lat }, wl: [w(0), w(1)], tc: u.thick * (s / (2.0 * u.hw).max(20.0)).min(1.0).max(0.3) })
            }
        }
    }

    /// Lane DATUM: (sharp, laterally smoothed) thickness of unit u at (x, y). Tephra and VE shapes are not smoothed.
    /// Lane DATUM2: (sharp, smoothed, lens-profile) thickness; the profile is None for units without a lens cross-section
    /// (tephra, VE shapes) and for non-lava units, Some(elliptic or flat-topped lobe per de_lens mode) for lava lenses.
    fn thick_trip(&self, sp: &StackParams, u: &Unit, x: f64, y: f64, seg: (usize, usize), nodes: &[(f64, f64)]) -> (f64, f64, Option<f64>) {
        match self.lens_at(u, x, y, seg) {
            None => { let t = self.thick_at(u, x, y, seg); (t, t, None) }
            Some(g) => {
                let sg = blur_sigma(sp, u.hw);
                let tb: f64 = nodes.iter().map(|&(z, wq)| wq * g.lens(g.l + sg * z)).sum();
                let t = self.thick_at(u, x, y, seg);
                let pr = match sp.lens.as_ref() { Some(lp) if u.mat == Mat::Lava => Some(if lp.mode == 2 { g.lobe(g.l) } else { t }), _ => None };
                (t, tb, pr)
            }
        }
    }

    /// Lane DATUM de_stack: the column at (x, y) stacked by levels (see STACK_NOTE).
    fn column_levels(&self, sp: &StackParams, x: f64, y: f64, d_surf: f64) -> LobeCol {
        let nodes = blur_nodes();
        let mut cov: Vec<(f64, u16, usize, f64, f64, Option<f64>)> = Vec::new(); // (age, id, unit, sharp t, smoothed t, lens profile)
        if let Some(list) = self.grid.get(&((x / GCELL).floor() as i64, (y / GCELL).floor() as i64)) {
            for &(ui, a, b) in list {
                let u = &self.units[ui as usize];
                let (t, tb, pr) = self.thick_trip(sp, u, x, y, (a as usize, b as usize), &nodes);
                if tb > 0.05 { cov.push((u.age, u.id, ui as usize, t, tb, pr)); }
            }
        }
        for &ui in &self.tephra {
            let u = &self.units[ui];
            let t = self.thick_at(u, x, y, (0, 0));
            if t > 0.05 { cov.push((u.age, u.id, ui, t, t, None)); }
        }
        cov.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        // groups: (unit index, unit id, top depth, base depth); levels are the cumulative smoothed thicknesses
        let mut groups: Vec<(usize, u16, f64, f64)> = Vec::new();
        if sp.lens.is_some() {
            // lane DATUM2 de_lens: bottom-up from the level stack base; lava tops = level top - own lens excess (LENS_NOTE)
            let total: f64 = cov.iter().map(|c| c.4).sum();
            let mut lt = total; // level top depth of the current unit (walking oldest -> youngest)
            let mut r = total; // depth of the top of the material placed so far
            for &(_, id, ui, t, tb, pr) in cov.iter().rev() {
                lt -= tb;
                if t <= 0.05 { continue; }
                let target = match pr { Some(p) => lt - (p - tb), None => lt }.max(0.0);
                if target >= r - 1e-9 { continue; }
                groups.push((ui, id, target, r));
                r = target;
            }
            groups.reverse();
            if let Some(g) = groups.first_mut() { g.2 = 0.0; }
        } else {
            let mut d = 0.0;
            for &(_, id, ui, t, tb, _) in &cov {
                if t > 0.05 { groups.push((ui, id, if groups.is_empty() { 0.0 } else { d }, d + tb)); }
                else if let Some(g) = groups.last_mut() { g.3 = d + tb; }
                d += tb;
            }
        }
        // above the youngest present unit (younger units absent here): that unit, extended up to D
        let sc3 = sc_of(crate::v2::LITHO_VOLCANIC);
        let mut parts = Vec::new();
        for (ui, id, top, base) in groups {
            let u = &self.units[ui];
            let t = base - top;
            if t <= 0.0 { continue; }
            let m = |l: f64| sc3 * l.exp() * self.norm;
            let mut dd = top;
            match u.mat {
                Mat::Lava => {
                    for (f, l, fac) in [(u.ft, u.lm[0], FAC_LAVA_TOP), (1.0 - u.ft - u.fb, u.lm[1], FAC_LAVA_CORE), (u.fb, u.lm[2], FAC_LAVA_BASE)] {
                        let tt = f * t;
                        if tt <= 0.0 { continue; }
                        parts.push(Part { d_top: dd, d_base: dd + tt, sc: m(l), fac, unit: id, ui: ui as u32, tu: t });
                        dd += tt;
                    }
                }
                Mat::Pdc | Mat::Tephra => self.pyro_parts(u, ui, id, dd, base, t, x, y, &m, &mut parts),
            }
        }
        let stack = parts.last().map(|p| p.d_base).unwrap_or(0.0);
        let mut cum = Vec::with_capacity(parts.len());
        let mut acc = 0.0;
        for p in &parts { acc += (p.d_base - p.d_top) / p.sc; cum.push(acc); }
        LobeCol { d_surf, parts, cum, stack_m: stack }
    }

    /// Lane DATUM2 de_seq=1 diagnostics: (Pbar, sum of the sharp local thicknesses of the covering units) at (x, y).
    pub fn seq_pile(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let pb = self.pbar.as_ref()?;
        let mut sum = 0.0;
        if let Some(list) = self.grid.get(&((x / GCELL).floor() as i64, (y / GCELL).floor() as i64)) {
            for &(ui, a, b) in list { let t = self.thick_at(&self.units[ui as usize], x, y, (a as usize, b as usize)); if t > 0.05 { sum += t; } }
        }
        for &ui in &self.tephra { let t = self.thick_at(&self.units[ui], x, y, (0, 0)); if t > 0.05 { sum += t; } }
        Some((pb.bilin(x, y).max(0.0), sum))
    }

    /// Lane DATUM2 de_seq=1: the column at (x, y) built bottom-up from the stack base at depth Pbar below D (SEQ_NOTE).
    fn column_seq(&self, pb: &SeqField, x: f64, y: f64, d_surf: f64) -> LobeCol {
        if let Some((sp, mt, st)) = self.seq_level.as_ref() { return self.column_seq_level(sp, *mt, *st, pb, x, y, d_surf); }
        let mut cov: Vec<(f64, u16, usize, f64)> = Vec::new();
        if let Some(list) = self.grid.get(&((x / GCELL).floor() as i64, (y / GCELL).floor() as i64)) {
            for &(ui, a, b) in list {
                let u = &self.units[ui as usize];
                let t = self.thick_at(u, x, y, (a as usize, b as usize));
                if t > 0.05 { cov.push((u.age, u.id, ui as usize, t)); }
            }
        }
        for &ui in &self.tephra {
            let u = &self.units[ui];
            let t = self.thick_at(u, x, y, (0, 0));
            if t > 0.05 { cov.push((u.age, u.id, ui, t)); }
        }
        // recorded units (T1b flows, VE units: ages below RECORDED_AGE) hang from D youngest first, as the record draws them
        // at the surface; the keyed pile is built oldest first from the stack base and truncated at the recorded units' base
        cov.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        let mut groups: Vec<(usize, u16, f64, f64, f64)> = Vec::new();
        let mut r_top = 0.0;
        for &(age, id, ui, t) in &cov { if age < RECORDED_AGE { groups.push((ui, id, r_top, r_top + t, t)); r_top += t; } }
        let base = pb.bilin(x, y).max(r_top);
        let mut kg: Vec<(usize, u16, f64, f64, f64)> = Vec::new();
        let mut r = base;
        for &(age, id, ui, t) in cov.iter().rev() {
            if age < RECORDED_AGE { continue; }
            if r <= r_top + 1e-9 { break; }
            let top = (r - t).max(r_top);
            kg.push((ui, id, top, r, t));
            r = top;
        }
        kg.reverse();
        // the pile ends below the recorded units (or D): the youngest present keyed unit extends up to them; lane PYRO2
        // de_seq_fill: the youngest THICK present keyed unit is thickened in place instead (SEQ_FILL_NOTE)
        match (self.seq_fill.as_ref(), kg.first().map(|g| g.2 - r_top)) {
            (Some(fp), Some(gap)) if gap > 1e-9 => {
                let j = kg.iter().position(|g| g.4 >= fp.min_ratio * gap)
                    .unwrap_or_else(|| kg.iter().enumerate().fold((0usize, f64::MIN), |a, (i, g)| if g.4 > a.1 { (i, g.4) } else { a }).0);
                for g in kg.iter_mut().take(j) { g.2 -= gap; g.3 -= gap; }
                kg[j].2 -= gap;
            }
            _ => { if let Some(g) = kg.first_mut() { g.2 = r_top; } }
        }
        groups.extend(kg);
        let sc3 = sc_of(crate::v2::LITHO_VOLCANIC);
        let mut parts = Vec::new();
        for (ui, id, top, base, tu) in groups {
            let u = &self.units[ui];
            let t = base - top;
            if t <= 0.0 { continue; }
            let m = |l: f64| sc3 * l.exp() * self.norm;
            let mut dd = top;
            match u.mat {
                Mat::Lava => {
                    for (f, l, fac) in [(u.ft, u.lm[0], FAC_LAVA_TOP), (1.0 - u.ft - u.fb, u.lm[1], FAC_LAVA_CORE), (u.fb, u.lm[2], FAC_LAVA_BASE)] {
                        let tt = f * t;
                        if tt <= 0.0 { continue; }
                        parts.push(Part { d_top: dd, d_base: dd + tt, sc: m(l), fac, unit: id, ui: ui as u32, tu });
                        dd += tt;
                    }
                }
                Mat::Pdc | Mat::Tephra => self.pyro_parts(u, ui, id, dd, base, tu, x, y, &m, &mut parts),
            }
        }
        let stack = parts.last().map(|p| p.d_base).unwrap_or(0.0);
        let mut cum = Vec::with_capacity(parts.len());
        let mut acc = 0.0;
        for p in &parts { acc += (p.d_base - p.d_top) / p.sc; cum.push(acc); }
        LobeCol { d_surf, parts, cum, stack_m: stack }
    }

    /// Lane ARCH de_seq_level: the de_seq=1 column with the keyed units levelled instead of stacked at their sharp lens
    /// thickness (SEQ_LEVEL_NOTE). Recorded units, the stack base, the truncation and the gap fill are as in `column_seq`.
    fn column_seq_level(&self, sp: &StackParams, mt: [f64; 2], stretch: bool, pb: &SeqField, x: f64, y: f64, d_surf: f64) -> LobeCol {
        let nodes = blur_nodes();
        // (age, id, unit, sharp t, smoothed t, lens profile, margin cap = distance to the sharp footprint edge x tan(margin
        // slope), infinite for units without a lens); recorded units keep their sharp thickness (they hang from D)
        let mut cov: Vec<(f64, u16, usize, f64, f64, Option<f64>, f64)> = Vec::new();
        if let Some(list) = self.grid.get(&((x / GCELL).floor() as i64, (y / GCELL).floor() as i64)) {
            for &(ui, a, b) in list {
                let u = &self.units[ui as usize];
                if u.age < RECORDED_AGE {
                    let t = self.thick_at(u, x, y, (a as usize, b as usize));
                    if t > 0.05 { cov.push((u.age, u.id, ui as usize, t, t, None, f64::INFINITY)); }
                    continue;
                }
                let (t, tb, pr) = self.thick_trip(sp, u, x, y, (a as usize, b as usize), &nodes);
                if t > 0.05 || tb > 0.05 {
                    let cap = match self.lens_at(u, x, y, (a as usize, b as usize)) {
                        Some(g) => (g.wl[(g.l < 0.0) as usize] - g.l.abs()).max(0.0) * if u.mat == Mat::Lava { mt[0] } else { mt[1] },
                        None => f64::INFINITY,
                    };
                    cov.push((u.age, u.id, ui as usize, t, tb, pr, cap));
                }
            }
        }
        for &ui in &self.tephra {
            let u = &self.units[ui];
            let t = self.thick_at(u, x, y, (0, 0));
            if t > 0.05 { cov.push((u.age, u.id, ui, t, t, None, f64::INFINITY)); }
        }
        cov.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        let mut groups: Vec<(usize, u16, f64, f64, f64)> = Vec::new();
        let mut r_top = 0.0;
        for &(age, id, ui, t, _, _, _) in &cov { if age < RECORDED_AGE { groups.push((ui, id, r_top, r_top + t, t)); r_top += t; } }
        let base = pb.bilin(x, y).max(r_top);
        // keyed units oldest first from the stack base: lt = level top depth, r = depth of the top of the material placed so far
        let mut kg: Vec<(usize, u16, f64, f64, f64)> = Vec::new();
        let (mut lt, mut r) = (base, base);
        for &(age, id, ui, t, tb, pr, cap) in cov.iter().rev() {
            if age < RECORDED_AGE { continue; }
            if r <= r_top + 1e-9 { break; }
            lt -= tb;
            if t <= 0.05 { continue; }
            let u = &self.units[ui];
            let target = if u.mat == Mat::Tephra { r - t } else { (match pr { Some(p) => lt - (p - tb), None => lt }).max(r - cap) };
            let top = target.max(r_top);
            if top >= r - 1e-9 { continue; }
            kg.push((ui, id, top, r, r - top));
            r = top;
        }
        kg.reverse();
        // the pile ends below the recorded units (or D): de_seq_stretch stretches the whole keyed pile (SEQ_STRETCH_NOTE);
        // otherwise as column_seq (youngest present keyed unit extended, or de_seq_fill)
        if stretch {
            if let Some(gap) = kg.first().map(|g| g.2 - r_top) {
                let pile = base - kg[0].2;
                if gap > 1e-9 && pile > 1e-9 {
                    let f = (base - r_top) / pile;
                    for g in kg.iter_mut() { g.2 = base - (base - g.2) * f; g.3 = base - (base - g.3) * f; g.4 = g.3 - g.2; }
                    kg[0].2 = r_top;
                }
            }
        } else {
        match (self.seq_fill.as_ref(), kg.first().map(|g| g.2 - r_top)) {
            (Some(fp), Some(gap)) if gap > 1e-9 => {
                let j = kg.iter().position(|g| g.4 >= fp.min_ratio * gap)
                    .unwrap_or_else(|| kg.iter().enumerate().fold((0usize, f64::MIN), |a, (i, g)| if g.4 > a.1 { (i, g.4) } else { a }).0);
                for g in kg.iter_mut().take(j) { g.2 -= gap; g.3 -= gap; }
                kg[j].2 -= gap;
            }
            _ => { if let Some(g) = kg.first_mut() { g.2 = r_top; } }
        }
        }
        groups.extend(kg);
        let sc3 = sc_of(crate::v2::LITHO_VOLCANIC);
        let mut parts = Vec::new();
        for (ui, id, top, base, tu) in groups {
            let u = &self.units[ui];
            let t = base - top;
            if t <= 0.0 { continue; }
            let m = |l: f64| sc3 * l.exp() * self.norm;
            let mut dd = top;
            match u.mat {
                Mat::Lava => {
                    for (f, l, fac) in [(u.ft, u.lm[0], FAC_LAVA_TOP), (1.0 - u.ft - u.fb, u.lm[1], FAC_LAVA_CORE), (u.fb, u.lm[2], FAC_LAVA_BASE)] {
                        let tt = f * t;
                        if tt <= 0.0 { continue; }
                        parts.push(Part { d_top: dd, d_base: dd + tt, sc: m(l), fac, unit: id, ui: ui as u32, tu });
                        dd += tt;
                    }
                }
                Mat::Pdc | Mat::Tephra => self.pyro_parts(u, ui, id, dd, base, tu, x, y, &m, &mut parts),
            }
        }
        let stack = parts.last().map(|p| p.d_base).unwrap_or(0.0);
        let mut cum = Vec::with_capacity(parts.len());
        let mut acc = 0.0;
        for p in &parts { acc += (p.d_base - p.d_top) / p.sc; cum.push(acc); }
        LobeCol { d_surf, parts, cum, stack_m: stack }
    }

    /// Lane PYRO: the part(s) of a pdc / lahar / tephra unit over depths [top, base] (tu = its local deposited thickness, base =
    /// its base): one part with the BS / DE strength and indurated flag, unless pyro_sub_hard is on; then the unit's subtype
    /// zones (pyro::UnitSub::zones) with subtype strengths (PYRO_HARD_NOTE), contiguous over [top, base].
    #[allow(clippy::too_many_arguments)]
    fn pyro_parts(&self, u: &Unit, ui: usize, id: u16, top: f64, base: f64, tu: f64, x: f64, y: f64, m: &dyn Fn(f64) -> f64, parts: &mut Vec<Part>) {
        let hard = self.pyro.as_ref().and_then(|(r, p)| p.hard.as_ref().map(|h| (r, p, h)));
        match (hard, u.sub) {
            (Some((rec, p, h)), Some((sub, z))) => {
                let r = sub.r_eff(x - self.vent[0], y - self.vent[1], p);
                let k = if u.mat == Mat::Tephra { rec.k_fall } else { rec.k_pdc };
                let zs = sub.zones(tu, r, p);
                let nz = zs.len();
                let mut push = |da: f64, db: f64, code: u8| {
                    let l = h.lnm[code as usize] + h.spread * z;
                    parts.push(Part { d_top: da, d_base: db, sc: m(l) * k, fac: if l > 0.0 { FAC_PYRO_HARD } else { FAC_PYRO }, unit: id, ui: ui as u32, tu });
                };
                // the column extension above the deposited top (de_seq) is the unit's body material (pyro::UnitSub::code_at)
                let mut cur = top;
                if base - tu > top + 1e-9 { push(top, base - tu, sub.body(tu, r, p)); cur = base - tu; }
                for (i, (_, fb, code)) in zs.into_iter().enumerate() {
                    let db = if i + 1 == nz { base } else { (base - (1.0 - fb) * tu).clamp(cur, base) };
                    if db <= cur + 1e-9 { continue; }
                    push(cur, db, code);
                    cur = db;
                }
            }
            _ => {
                let fac = if u.lm[0] > 0.0 { FAC_PYRO_HARD } else { FAC_PYRO };
                parts.push(Part { d_top: top, d_base: base, sc: m(u.lm[0]), fac, unit: id, ui: ui as u32, tu });
            }
        }
    }

    /// Local thickness of unit u at (x, y) (0 outside); path units scan segments seg.0..=seg.1 (those within reach of the cell).
    fn thick_at(&self, u: &Unit, x: f64, y: f64, seg: (usize, usize)) -> f64 {
        if let Some(sh) = u.veg { return sh.thick(x, y); }
        // lane LOBE de_lobe_relief: a recorded lane VE lava flow has the relief section's thickness (RELIEF_NOTE)
        if let (Some(ri), Some((rp, ru))) = (u.rel, self.relief.as_ref()) {
            return crate::relief::section(rp, &ru[ri as usize], x, y).map(|(a, b)| (a + b).max(0.0)).unwrap_or(0.0);
        }
        match u.mat {
            Mat::Tephra => {
                let r = (x - self.vent[0]).hypot(y - self.vent[1]);
                let th = (y - self.vent[1]).atan2(x - self.vent[0]);
                // downwind skew: the decay length is 1.6x downwind, 0.4x upwind
                let l = u.run * (1.0 + 0.6 * ang_diff(th, u.path[0][0]).cos());
                u.thick * (-r / l.max(1.0)).exp()
            }
            _ => match u.t1b {
                Some(f) => {
                    // T1b lava_flow_lobes_v1 (pf-region t1b_coast.rs cone_cf): radial axis with wander, ramp-in, snout
                    let (dx, dy) = (x - self.vent[0], y - self.vent[1]);
                    let r = dx.hypot(dy);
                    let d = ang_diff(dy.atan2(dx), f[0]);
                    if d.abs() >= std::f64::consts::FRAC_PI_2 { return 0.0; }
                    let off = f[5] * (std::f64::consts::TAU * r / f[6].max(1.0) + f[7]).sin();
                    let lat_s = r * d.sin() - off;
                    let (along, lat) = (r * d.cos(), lat_s.abs());
                    let s = along - f[1];
                    if s < 0.0 || s > f[2] { return 0.0; }
                    let w = f[3] * (1.0 + 0.5 * s / f[2]) * u.marg.as_ref().map(|m| marg_at(m, (lat_s < 0.0) as usize, s)).unwrap_or(1.0);
                    let wl = if s > f[2] - w { let q = ((s - (f[2] - w)) / w).min(1.0); w * (1.0 - q * q).max(0.0).sqrt() } else { w };
                    if lat >= wl || wl <= 0.0 { return 0.0; }
                    f[4] * (1.0 - (lat / wl).powi(2)).sqrt() * (s / 200.0).min(1.0)
                }
                None => {
                    // nearest point on the path: lateral distance and along-path position; same cross-section rules
                    let mut best = (f64::INFINITY, 0.0f64);
                    let mut side = 0usize;
                    let hi = (seg.1 + 2).min(u.path.len());
                    for w in u.path[seg.0.min(hi.saturating_sub(2))..hi].windows(2) {
                        let (ax, ay, bx, by) = (w[0][0], w[0][1], w[1][0], w[1][1]);
                        let (vx, vy) = (bx - ax, by - ay);
                        let l2 = vx * vx + vy * vy;
                        let t = if l2 > 0.0 { (((x - ax) * vx + (y - ay) * vy) / l2).clamp(0.0, 1.0) } else { 0.0 };
                        let dd = (x - ax - t * vx).hypot(y - ay - t * vy);
                        if dd < best.0 { best = (dd, w[0][2] + t * (w[1][2] - w[0][2])); side = (vx * (y - ay) - vy * (x - ax) < 0.0) as usize; }
                    }
                    let (lat, s) = best;
                    let w = u.hw * (1.0 + 0.5 * s / u.run.max(1.0)) * u.marg.as_ref().map(|m| marg_at(m, side, s)).unwrap_or(1.0);
                    let wl = if s > u.run - w { let q = ((s - (u.run - w)) / w).min(1.0); w * (1.0 - q * q).max(0.0).sqrt() } else { w };
                    if lat >= wl || wl <= 0.0 { return 0.0; }
                    u.thick * (1.0 - (lat / wl).powi(2)).sqrt() * (s / (2.0 * u.hw).max(20.0)).min(1.0).max(0.3)
                }
            },
        }
    }

    /// The local column at (x, y) below the construction surface d_surf (parts top-down; 1/Sc integral at part bases).
    pub fn column(&self, x: f64, y: f64, d_surf: f64) -> LobeCol {
        if let Some(pb) = self.pbar.as_ref() { return self.column_seq(pb, x, y, d_surf); }
        if let Some(sp) = self.stack.as_ref() { return self.column_levels(sp, x, y, d_surf); }
        let mut cov: Vec<(f64, u16, usize, f64)> = Vec::new(); // (age, id, unit, thickness)
        if let Some(list) = self.grid.get(&((x / GCELL).floor() as i64, (y / GCELL).floor() as i64)) {
            for &(ui, a, b) in list {
                let u = &self.units[ui as usize];
                let t = self.thick_at(u, x, y, (a as usize, b as usize));
                if t > 0.05 { cov.push((u.age, u.id, ui as usize, t)); }
            }
        }
        for &ui in &self.tephra {
            let u = &self.units[ui];
            let t = self.thick_at(u, x, y, (0, 0));
            if t > 0.05 { cov.push((u.age, u.id, ui, t)); }
        }
        cov.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        let sc3 = sc_of(crate::v2::LITHO_VOLCANIC);
        let mut parts = Vec::new();
        let mut d = 0.0;
        for (_, id, ui, t) in cov {
            let u = &self.units[ui];
            let m = |l: f64| sc3 * l.exp() * self.norm;
            match u.mat {
                Mat::Lava => {
                    for (f, l, fac) in [(u.ft, u.lm[0], FAC_LAVA_TOP), (1.0 - u.ft - u.fb, u.lm[1], FAC_LAVA_CORE), (u.fb, u.lm[2], FAC_LAVA_BASE)] {
                        let tt = f * t;
                        if tt <= 0.0 { continue; }
                        parts.push(Part { d_top: d, d_base: d + tt, sc: m(l), fac, unit: id, ui: ui as u32, tu: t });
                        d += tt;
                    }
                }
                Mat::Pdc | Mat::Tephra => {
                    self.pyro_parts(u, ui, id, d, d + t, t, x, y, &m, &mut parts);
                    d += t;
                }
            }
        }
        let mut cum = Vec::with_capacity(parts.len());
        let mut acc = 0.0;
        for p in &parts { acc += (p.d_base - p.d_top) / p.sc; cum.push(acc); }
        LobeCol { d_surf, parts, cum, stack_m: d }
    }
}

/// Lane DATUM de_route: mean of the T1b height over (2 h + 1)^2 posts at every window post (posts clamped to the lattice, as
/// `Aux::build`'s edifice datum S). A pure function of the window: identical in every block.
pub fn box_mean(w: &T1bWindow, h: i64) -> Vec<f64> {
    let n = w.n as i64;
    let mut out = vec![0.0f64; w.n * w.n];
    let cnt = ((2 * h + 1) * (2 * h + 1)) as f64;
    for j in 0..n { for i in 0..n {
        let mut sum = 0.0;
        for dy in -h..=h { for dx in -h..=h { sum += w.h(w.gx0 + i + dx, w.gy0 + j + dy); } }
        out[(j * n + i) as usize] = sum / cnt;
    } }
    out
}

/// Lane DATUM2 de_seq: box mean over (2 h + 1)^2 window posts of `pre` (a full-window lattice, e.g. the pre-cut surface) or
/// of the height (clamped to the lattice as `T1bWindow::h`).
pub fn box_mean_of(w: &T1bWindow, pre: Option<&[f64]>, h: i64) -> Vec<f64> {
    let Some(f) = pre else { return box_mean(w, h) };
    let n = w.n as i64;
    let at = |i: i64, j: i64| f[(j.clamp(0, n - 1) * n + i.clamp(0, n - 1)) as usize];
    let mut out = vec![0.0f64; w.n * w.n];
    let cnt = ((2 * h + 1) * (2 * h + 1)) as f64;
    for j in 0..n { for i in 0..n {
        let mut sum = 0.0;
        for dy in -h..=h { for dx in -h..=h { sum += at(i + dx, j + dy); } }
        out[(j * n + i) as usize] = sum / cnt;
    } }
    out
}

/// The lobe world of a block: per volcano index (T1b list) its units.
pub struct LobeWorld { pub volc: std::collections::BTreeMap<u8, VolcLobes>, pub info: Vec<serde_json::Value>,
    /// Lane LOBE de_lobe_relief: the relief of the block's volcanoes (None = off).
    pub relief: Option<crate::relief::ReliefWorld> }

impl LobeWorld {
    /// Lane LOBE: the constructional relief R at chart (x, y) (0 when de_lobe_relief is off).
    pub fn relief_at(&self, x: f64, y: f64) -> f64 { self.relief.as_ref().map(|r| r.at(x, y)).unwrap_or(0.0) }
    /// Build for the volcanoes `used` (T1b list indices) over the block bbox.
    pub fn build(lp: &LobeParams, w: &T1bWindow, used: &[u8], bbox: (f64, f64, f64, f64), scope: &Scope, pre: Option<&[f64]>, pyro: Option<&crate::pyro::PyroWorld>) -> R<LobeWorld> {
        let mut volc = std::collections::BTreeMap::new();
        let mut info = Vec::new();
        // lane DATUM de_route: the box mean of the T1b height over the whole window lattice (clamped like Aux), bicubic
        // lane DATUM2 de_seq: the routing base C, the box mean of the pre-cut surface (ed_datum) or the height
        let boxm: Option<Vec<f64>> = match (lp.seq.as_ref(), lp.route.as_ref()) {
            (Some(sp), _) => Some(box_mean_of(w, pre, sp.box_half_posts)),
            (None, Some(rp)) => Some(box_mean(w, rp.box_half_posts)),
            _ => None,
        };
        let dom = lp.seq.as_ref().map(|_| { let t = crate::T1B_M; (w.gx0 as f64 * t, w.gy0 as f64 * t, (w.gx0 + w.n as i64 - 1) as f64 * t, (w.gy0 + w.n as i64 - 1) as f64 * t) });
        for &v in used {
            let Some(vi) = w.volcanoes.get(v as usize) else { continue };
            let geo = w.volcano_geo.get(v as usize).copied().unwrap_or([f64::NAN; 4]);
            if !geo[0].is_finite() { continue; }
            let vent = [geo[0], geo[1]];
            let flows = w.volcano_flows.get(v as usize).cloned().unwrap_or_default();
            let on_ed = |x: f64, y: f64| -> bool {
                let Some(i) = w.nearest(x, y) else { return false };
                if w.fixed[i] == 2 || w.fixed[i] == 4 { return false; }
                if w.exposed_litho.as_ref().map(|e| e[i]) != Some(crate::v2::LITHO_VOLCANIC) { return false; }
                match &w.volcano_id { Some(f) => f[i] == v, None => (x - vent[0]).hypot(y - vent[1]) <= geo[2] }
            };
            let zf = |x: f64, y: f64| match &boxm { Some(b) => w.bicubic_of(b, x, y), None => w.bicubic(x, y) };
            let style_at = |x: f64, y: f64| -> Option<Style> {
                let f = w.edifice_style.as_ref()?;
                let i = w.nearest(x, y)?;
                style_of_ve(f[i])
            };
            let st = style_of_kind(&vi.kind);
            let ve = if lp.ve == Some(true) { w.volcano_ve_units.get(v as usize).and_then(|o| o.as_ref()) } else { None };
            let pr = pyro.and_then(|pw| pw.recs.get(v as usize).map(|r| (r, &pw.p)));
            let rbody = w.volcano_r_body.get(v as usize).copied().flatten().unwrap_or(geo[2]);
            let vl = VolcLobes::build(lp, st, &vi.id, vent, &flows, bbox, scope, &on_ed, &zf, &style_at, ve, dom, pr, rbody)?;
            let n_lobe = vl.units.iter().filter(|u| u.mat != Mat::Tephra && u.t1b.is_none()).count();
            info.push(serde_json::json!({"volcano_index": v, "volcano": vi.id, "kind": vi.kind, "style": vl.style, "recorded_t1b_flows": vl.recorded,
                "keyed_lobes": n_lobe, "tephra_falls": vl.tephra.len(), "strength_norm": vl.norm, "ve_edifice_style_field": w.edifice_style.is_some(),
                "ve_units_used": if ve.is_some() { serde_json::json!(vl.units.iter().filter(|u| u.age < 0.02 && (u.veg.is_some() || u.t1b.is_none() && u.mat != Mat::Tephra && u.age < 0.02)).count()) } else { serde_json::Value::Null },
                "ve_units_skipped": if ve.is_some() { serde_json::json!(vl.ve_skipped) } else { serde_json::Value::Null }}));
            if lp.seq.is_some() { if let Some(q) = info.last_mut() { q["seq_units"] = serde_json::json!(vl.seq_units); q["seq_emplaced"] = serde_json::json!(vl.seq_stats[0]); q["seq_run_over_drawn"] = serde_json::json!(vl.seq_stats[1]); q["seq_overtopped_share"] = serde_json::json!(vl.seq_stats[2]); } }
            if let (Some((_, ru)), Some(q)) = (vl.relief.as_ref(), info.last_mut()) { q["relief_units"] = serde_json::json!(ru.len()); q["relief_lava_flows"] = serde_json::json!(ru.iter().filter(|u| u.lava).count()); }
            volc.insert(v, vl);
        }
        // lane LOBE de_lobe_relief: one relief world over the block (all used volcanoes' lane VE units)
        let relief = lp.relief.as_ref().map(|rp| {
            let all: Vec<crate::relief::RUnit> = volc.values().filter_map(|vl| vl.relief.as_ref()).flat_map(|r| r.1.iter().cloned()).collect();
            crate::relief::ReliefWorld::build(rp, all, bbox)
        });
        Ok(LobeWorld { volc, info, relief })
    }
}

/// Bed-coordinate helpers of a lobe column with the per-volcano 1D sequence below the stack (`ed` column, depth from its datum).
pub struct LobeAt<'a> { pub col: LobeCol, pub ed: &'a crate::strata::CellStrata, pub ed_off: f64,
    /// Lane PYRO: the volcano's units (subtype lookup); None in callers that do not publish subtypes.
    pub vl: Option<&'a VolcLobes> }

impl<'a> LobeAt<'a> {
    /// Bed coordinate of surface height z: u = z - D (negative below the construction surface).
    #[inline] pub fn u(&self, z: f64) -> f64 { z - self.col.d_surf }
    fn below_sc(&self, dd: f64) -> f64 { let c = &self.ed.col; c.bed_sc[c.bed_at_depth(dd + self.ed_off)] }
    fn below_int(&self, dd: f64) -> f64 {
        if self.ed.cum_is.is_empty() { return dd / sc_of(crate::v2::LITHO_VOLCANIC); }
        // the ed table is over stratigraphic height u = -depth; integral of 1/Sc from depth 0 to dd = I(u=-0) - I(u=-dd)
        let n = self.ed.cum_is.len();
        let at = |u: f64| -> f64 {
            let f = ((u - self.ed.u0) / self.ed.du).clamp(0.0, (n - 1) as f64);
            let j = (f.floor() as usize).min(n - 2); let t = f - j as f64;
            self.ed.cum_is[j] + t * (self.ed.cum_is[j + 1] - self.ed.cum_is[j])
        };
        at(-self.ed_off) - at(-(dd + self.ed_off))
    }
    /// Sc of the bed at bed coordinate u.
    pub fn sc(&self, u: f64) -> f64 {
        let d = -u;
        match self.col.part_at(d) { Some(i) => self.col.parts[i].sc, None => if d < 0.0 { self.col.parts.first().map(|p| p.sc).unwrap_or_else(|| self.below_sc(0.0)) } else { self.below_sc(d - self.col.stack_m) } }
    }
    /// Mean of 1/Sc over [u1, u2].
    pub fn inv_mean(&self, u1: f64, u2: f64) -> f64 {
        let (a, b) = if u1 <= u2 { (u1, u2) } else { (u2, u1) };
        if b - a < 1e-6 { return 1.0 / self.sc(a); }
        let bi = |dd: f64| self.below_int(dd);
        (self.col.integ(-a, &bi) - self.col.integ(-b, &bi)) / (b - a)
    }
    /// (facies, unit id) at bed coordinate u; below the stack: the ed column's (facies, unit) via the caller.
    pub fn meta(&self, u: f64) -> Option<(u8, u16)> {
        let d = -u;
        let i = if d < 0.0 { if self.col.parts.is_empty() { return None; } 0 } else { self.col.part_at(d)? };
        Some((self.col.parts[i].fac, self.col.parts[i].unit))
    }
    /// Lane RIB: facies at bed coordinate u, continuing into the per-volcano sequence below the stack.
    pub fn fac_any(&self, u: f64) -> u8 {
        match self.meta(u) {
            Some((f, _)) => f,
            None => { let c = &self.ed.col; let dd = -u - self.col.stack_m; c.bed_fac.get(c.bed_at_depth(dd + self.ed_off)).copied().unwrap_or(crate::bedstyle::FAC_NONE) }
        }
    }
    /// Lane PYRO: pyroclastic subtype code at bed coordinate u (pyro::SUB_*), at chart (x, y); None below the stack (the caller
    /// reads the per-volcano sequence). 255 for a non-pyroclastic part; 254 when the unit has no subtype (pyro off).
    pub fn pyro_code(&self, u: f64, x: f64, y: f64) -> Option<u8> {
        let d = -u;
        let i = if d < 0.0 { if self.col.parts.is_empty() { return None; } 0 } else { self.col.part_at(d)? };
        let p = &self.col.parts[i];
        if p.fac != FAC_PYRO && p.fac != FAC_PYRO_HARD { return Some(crate::pyro::SUB_NONE); }
        let (Some(vl), Some(un)) = (self.vl, self.vl.and_then(|v| v.units.get(p.ui as usize))) else { return Some(crate::pyro::SUB_UNRESOLVED) };
        let (Some((sub, _)), Some((_, pp))) = (un.sub, vl.pyro.as_ref()) else { return Some(crate::pyro::SUB_UNRESOLVED) };
        // the unit's base: the deepest of its parts here (its deposited top lies tu above it)
        let base = self.col.parts.iter().filter(|q| q.ui == p.ui).map(|q| q.d_base).fold(p.d_base, f64::max);
        let frac = 1.0 - (base - d.max(p.d_top)) / p.tu.max(1e-6);
        Some(sub.code_at(frac, p.tu, sub.r_eff(x - vl.vent[0], y - vl.vent[1], pp), pp))
    }
    /// (depth below the unit top, height above the unit base) of the unit containing u (None below the stack).
    pub fn unit_contacts(&self, u: f64) -> Option<(f64, f64)> {
        let d = (-u).max(0.0);
        let i = self.col.part_at(d)?;
        let id = self.col.parts[i].unit;
        let top = self.col.parts.iter().filter(|p| p.unit == id).map(|p| p.d_top).fold(f64::INFINITY, f64::min);
        let base = self.col.parts.iter().filter(|p| p.unit == id).map(|p| p.d_base).fold(0.0f64, f64::max);
        Some((d - top, base - d))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn style_norm_is_positive_and_stable() {
        let s = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let a = style_norm(&STYLE_ARC, &s, "x").unwrap();
        let b = style_norm(&STYLE_ARC, &s, "x").unwrap();
        assert_eq!(a, b);
        assert!(a > 0.5 && a < 3.0, "norm {a}");
    }
    #[test]
    fn lobate_margins_are_bounded_keyed_and_vary() {
        let s = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let lb = LobateParams::default();
        let a = draw_marg(&lb, &s, "u/1", 60.0).unwrap();
        assert_eq!(a, draw_marg(&lb, &s, "u/1", 60.0).unwrap());
        assert_ne!(a, draw_marg(&lb, &s, "u/2", 60.0).unwrap());
        let v: Vec<f64> = (0..400).map(|i| marg_at(&a, 0, i as f64 * 2.0)).collect();
        assert!(v.iter().all(|&m| m >= 1.0 - lb.amp - 1e-12 && m <= 1.0 + lb.amp + 1e-12));
        let (lo, hi) = v.iter().fold((f64::MAX, f64::MIN), |q, &m| (q.0.min(m), q.1.max(m)));
        assert!(hi - lo > 0.1, "margin range {lo}..{hi}");
        // wavelengths scale with the unit half-width (clamped)
        for side in a.iter() { for k in 0..3 { let lam = std::f64::consts::TAU / side[3 * k + 1]; assert!((90.0 - 1e-9..=240.0 + 1e-9).contains(&lam), "lambda {lam}"); } }
    }
    fn test_world(stack: Option<StackParams>, marg: bool) -> VolcLobes {
        let lb = LobateParams::default();
        let s = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let mk = |id: u16, age: f64, x0: f64, y0: f64, x1: f64, y1: f64, t: f64, hw: f64| {
            let n = 20;
            let path: Vec<[f64; 3]> = (0..=n).map(|i| { let f = i as f64 / n as f64; [x0 + f * (x1 - x0), y0 + f * (y1 - y0), f * (x1 - x0).hypot(y1 - y0)] }).collect();
            let run = path.last().unwrap()[2];
            Unit { age, mat: Mat::Lava, path, thick: t, hw, run, ft: 0.2, fb: 0.1, lm: [0.0; 3], t1b: None, veg: None, id,
                   marg: if marg { Some(draw_marg(&lb, &s, &format!("u/{id}"), hw).unwrap()) } else { None }, sub: None, rel: None }
        };
        // two crossing keyed lobes and one recorded T1b flow (radial from the vent at the origin)
        let mut units = vec![mk(1, 0.3, -400.0, 0.0, 600.0, 0.0, 12.0, 60.0), mk(2, 0.6, 0.0, -400.0, 0.0, 600.0, 20.0, 40.0)];
        units.push(Unit { age: 0.01, mat: Mat::Lava, path: vec![], thick: 15.0, hw: 120.0, run: 900.0, ft: 0.2, fb: 0.1, lm: [0.0; 3],
                          t1b: Some([0.3, 100.0, 900.0, 120.0, 15.0, 20.0, 300.0, 0.5]), veg: None, id: 3, marg: None, sub: None, rel: None });
        let mut grid: HashMap<(i64, i64), Vec<(u32, u16, u16)>> = HashMap::new();
        for gy in -20..20 { for gx in -20..20 { grid.insert((gx, gy), vec![(0, 0, 19), (1, 0, 19), (2, 0, 0)]); } }
        VolcLobes { style: "test", units, tephra: vec![], grid, norm: 1.0, recorded: 1, vent: [0.0, 0.0], ve_skipped: 0, stack, pbar: None, seq_units: 0, seq_stats: [0.0; 3], pyro: None, seq_fill: None, seq_level: None, relief: None }
    }
    #[test]
    fn lens_geometry_reproduces_sharp_thickness() {
        for marg in [false, true] {
            let w = test_world(None, marg);
            for u in &w.units {
                for k in 0..400 {
                    let (x, y) = (-300.0 + 1.7 * k as f64, -250.0 + 1.3 * k as f64);
                    let seg = if u.t1b.is_some() { (0, 0) } else { (0, 19) };
                    let t = w.thick_at(u, x, y, seg);
                    let g = w.lens_at(u, x, y, seg).map(|g| g.lens(g.l)).unwrap_or(0.0);
                    assert!((t - g).abs() < 1e-9, "unit {} at ({x}, {y}): {t} vs {g}", u.id);
                }
            }
        }
    }
    #[test]
    fn level_stacking_is_contiguous_and_gentler() {
        let sp = StackParams::default();
        let (w0, w1) = (test_world(None, true), test_world(Some(sp), true));
        let (mut slope0, mut slope1, mut prev0, mut prev1) = (0.0f64, 0.0f64, None::<f64>, None::<f64>);
        for k in 0..300 {
            // a transect across the lobe 2 axis (x = 0) inside lobe 1 (y = 10)
            let (x, y) = (-150.0 + k as f64, 10.0);
            let c0 = w0.column(x, y, 100.0);
            let c1 = w1.column(x, y, 100.0);
            // parts contiguous from 0 to the stack base
            if !c1.parts.is_empty() {
                assert!(c1.parts[0].d_top.abs() < 1e-12);
                for q in c1.parts.windows(2) { assert!((q[0].d_base - q[1].d_top).abs() < 1e-9 && q[1].d_base >= q[1].d_top); }
                assert!((c1.stack_m - c1.parts.last().unwrap().d_base).abs() < 1e-12);
            }
            // every unit with sharp thickness here is present in the level column, and only those
            let mut ids0: Vec<u16> = c0.parts.iter().map(|p| p.unit).collect(); ids0.dedup();
            let mut ids1: Vec<u16> = c1.parts.iter().map(|p| p.unit).collect(); ids1.dedup();
            ids0.sort(); ids1.sort();
            assert_eq!(ids0, ids1, "at x {x}");
            // stack base depth below D (the deepest contact): steepest step per 1 m along the transect
            let (b0, b1) = (c0.stack_m, c1.stack_m);
            if let Some(p) = prev0 { slope0 = slope0.max((b0 - p).abs()); }
            if let Some(p) = prev1 { slope1 = slope1.max((b1 - p).abs()); }
            prev0 = Some(b0); prev1 = Some(b1);
        }
        assert!(slope1 < 0.5 && slope1 < slope0, "level stack base slope {slope1} m/m vs sharp {slope0}");
    }
    #[test]
    fn flat_based_lenses_are_contiguous_and_keep_the_stack_base() {
        for mode in [1u8, 2] {
            let mut sl = StackParams::default(); sl.lens = Some(LensParams::new(mode));
            let (w1, w2) = (test_world(Some(StackParams::default()), true), test_world(Some(sl), true));
            let (mut base_step, mut top_relief) = (0.0f64, (f64::MAX, f64::MIN));
            let mut prev: Option<f64> = None;
            for k in 0..300 {
                // a transect across the lobe 2 axis (x = 0; lobe 2 is the oldest unit) inside lobe 1 (y = 10)
                let (x, y) = (-150.0 + k as f64, 10.0);
                let (c1, c2) = (w1.column(x, y, 100.0), w2.column(x, y, 100.0));
                assert!((c1.stack_m - c2.stack_m).abs() < 1e-9, "stack base unchanged at x {x}");
                if c2.parts.is_empty() { continue; }
                assert!(c2.parts[0].d_top.abs() < 1e-12);
                for q in c2.parts.windows(2) { assert!((q[0].d_base - q[1].d_top).abs() < 1e-9 && q[1].d_base >= q[1].d_top); }
                // lobe 2 (id 2) sits on the stack base: its base is the (smooth) base, its top carries the lens relief
                let p2: Vec<&Part> = c2.parts.iter().filter(|p| p.unit == 2).collect();
                if let (Some(a), Some(b)) = (p2.first(), p2.last()) {
                    if let Some(pb) = prev { base_step = base_step.max((b.d_base - pb).abs()); }
                    prev = Some(b.d_base);
                    top_relief = (top_relief.0.min(a.d_top), top_relief.1.max(a.d_top));
                }
            }
            assert!(base_step < 0.5, "mode {mode}: lobe 2 base step {base_step} m/m");
            assert!(top_relief.1 - top_relief.0 > 2.0, "mode {mode}: lobe 2 top relief {top_relief:?}");
        }
    }
    #[test]
    fn sequential_pass_is_block_independent_and_columns_stack_on_pbar() {
        let s = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let mut lp = LobeParams::default();
        lp.seq = Some(SeqParams::new(1));
        lp.cell_m = 512.0;
        let on_ed = |x: f64, y: f64| x.hypot(y) < 2500.0;
        let zf = |x: f64, y: f64| 500.0 - 0.15 * x.hypot(y) + 4.0 * (x / 300.0).sin();
        let style_at = |_x: f64, _y: f64| None::<Style>;
        let dom = Some((-3000.0, -3000.0, 3000.0, 3000.0));
        let b = |bbox| VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        let (a, c) = (b((-1000.0, -1000.0, 0.0, 0.0)), b((200.0, 200.0, 1400.0, 1400.0)));
        assert!(a.seq_units > 20, "units {}", a.seq_units);
        assert_eq!(a.units.len(), c.units.len());
        for (u, v) in a.units.iter().zip(c.units.iter()) { assert_eq!(u.id, v.id); assert_eq!(u.path, v.path); }
        let (pa, pc) = (a.pbar.as_ref().unwrap(), c.pbar.as_ref().unwrap());
        assert_eq!(pa.v, pc.v);
        assert!(pa.v.iter().cloned().fold(0.0f64, f64::max) > 1.0);
        // columns: contiguous from D to the stack base; the base is Pbar where the pile reaches D (no part above D)
        let mut n = 0;
        for k in 0..200 {
            let (x, y) = (-900.0 + 4.0 * k as f64, -600.0 + 2.0 * k as f64);
            let col = a.column(x, y, 100.0);
            if col.parts.is_empty() { continue; }
            n += 1;
            assert!(col.parts[0].d_top.abs() < 1e-12);
            for q in col.parts.windows(2) { assert!((q[0].d_base - q[1].d_top).abs() < 1e-9 && q[1].d_base >= q[1].d_top); }
            assert!(col.stack_m <= pa.bilin(x, y) + 1e-9);
        }
        assert!(n > 20, "columns with parts {n}");
    }
    #[test]
    fn stack_and_route_keys_are_opt_in() {
        let lp = LobeParams::default();
        let j = serde_json::to_string(&lp).unwrap();
        assert!(!j.contains("stack") && !j.contains("route"), "{j}");
        let mut p = crate::DetailParamsV1::v3();
        assert!(p.set("de_stack", "1").is_err());
        p.set("beds", "2").unwrap(); p.set("de", "1").unwrap();
        assert!(p.set("de_stack", "1").is_err(), "needs de_lobes");
        p.set("de_lobes", "1").unwrap(); p.set("de_stack", "1").unwrap(); p.set("de_route", "1").unwrap(); p.set("de_stack_blur", "0.5").unwrap();
        let l = p.de.as_ref().unwrap().lobes.as_ref().unwrap();
        assert_eq!(l.stack.as_ref().unwrap().blur_hw_mult, 0.5);
        // lane DATUM2 de_lens: a de_stack parameter, absent from the params JSON when off
        assert!(!serde_json::to_string(l).unwrap().contains("\"lens\""));
        assert_eq!(l.route.as_ref().unwrap().box_half_posts, 4);
        p.set("de_lens", "2").unwrap();
        assert_eq!(p.de.as_ref().unwrap().lobes.as_ref().unwrap().stack.as_ref().unwrap().lens.as_ref().unwrap().mode, 2);
        assert!(p.set("de_lens", "3").is_err());
        let mut r = crate::DetailParamsV1::v3();
        r.set("beds", "2").unwrap(); r.set("de", "1").unwrap(); r.set("de_lobes", "1").unwrap();
        assert!(r.set("de_lens", "1").is_err(), "needs de_stack");
        // lane DATUM2 de_seq: a lobes parameter, absent from the params JSON when off
        r.set("de_seq", "1").unwrap(); r.set("de_seq_sigma", "300").unwrap();
        let q = r.de.as_ref().unwrap().lobes.as_ref().unwrap().seq.as_ref().unwrap();
        assert_eq!((q.mode, q.base_sigma_m), (1, 300.0));
        assert!(r.set("de_seq", "3").is_err());
        assert!(!serde_json::to_string(&LobeParams::default()).unwrap().contains("\"seq\""));
        // lane DATUM ed_datum: beds parameter, absent from the params JSON when off
        assert!(!serde_json::to_string(p.beds.as_ref().unwrap()).unwrap().contains("ed_datum"));
        p.set("ed_datum", "1").unwrap();
        assert!(serde_json::to_string(p.beds.as_ref().unwrap()).unwrap().contains("ed_datum"));
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ed_datum", "1").is_err(), "needs beds=2");
    }
    #[test]
    fn pyro_sub_is_look_only_unless_hard() {
        // lane PYRO: with pyro_sub the columns are bit-identical (depths, strengths, facies, units) and every pdc / tephra unit
        // has a subtype; with pyro_sub_hard the pyroclastic parts split into subtype zones, contiguous over the same stack
        let s = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let mut lp = LobeParams::default();
        lp.seq = Some(SeqParams::new(1));
        lp.cell_m = 512.0;
        let on_ed = |x: f64, y: f64| x.hypot(y) < 2500.0;
        let zf = |x: f64, y: f64| 500.0 - 0.15 * x.hypot(y) + 4.0 * (x / 300.0).sin();
        let style_at = |_x: f64, _y: f64| None::<Style>;
        let dom = Some((-3000.0, -3000.0, 3000.0, 3000.0));
        let counts = serde_json::json!({"cinder_cone": 3, "dome": 3, "explosive_vei3": 8, "explosive_vei5": 1, "lahar": 11});
        let pp = crate::pyro::PyroParams::default();
        let rec = crate::pyro::PyroRec::new("v", "continental_arc", Some(&counts), [0.0, 0.0], 2500.0, &pp);
        let mut ph = crate::pyro::PyroParams::default(); ph.hard = Some(crate::pyro::PyroHard::default());
        let rech = crate::pyro::PyroRec::new("v", "continental_arc", Some(&counts), [0.0, 0.0], 2500.0, &ph);
        let bbox = (-1500.0, -1500.0, 1500.0, 1500.0);
        let off = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        let on = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, Some((&rec, &pp)), 0.0).unwrap();
        let hd = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, Some((&rech, &ph)), 0.0).unwrap();
        assert!(on.units.iter().all(|u| (u.mat == Mat::Lava) == u.sub.is_none()));
        assert!(off.units.iter().all(|u| u.sub.is_none()));
        let ed = crate::strata::CellStrata { cell: 0, top_m: 0.0, g: (0.0, 0.0), g_raw: (0.0, 0.0), u0: 0.0, du: 1.0, psi: vec![], max_drift_m: 0.0, sc_win: vec![], k_win: vec![], cum_is: vec![],
            col: crate::strata::Column { lith: vec![3], base_depth: vec![1e4], bed_base: vec![1e4], bed_sc: vec![1.7], bed_k: vec![1.0], up_top: vec![], n_down: 1, bed_fac: vec![1], bed_unit: vec![0], unit_top: vec![] } };
        let (mut n, mut split, mut codes) = (0, 0, std::collections::BTreeSet::new());
        for k in 0..400 {
            let (x, y) = (-1400.0 + 7.0 * k as f64, -900.0 + 4.5 * k as f64);
            let (a, b, c) = (off.column(x, y, 100.0), on.column(x, y, 100.0), hd.column(x, y, 100.0));
            assert_eq!(a.parts.len(), b.parts.len());
            for (p, q) in a.parts.iter().zip(b.parts.iter()) {
                assert!(p.d_top.to_bits() == q.d_top.to_bits() && p.d_base.to_bits() == q.d_base.to_bits() && p.sc.to_bits() == q.sc.to_bits() && p.fac == q.fac && p.unit == q.unit);
            }
            assert_eq!(a.cum.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), b.cum.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
            if c.parts.is_empty() { continue; }
            n += 1;
            assert!((c.stack_m - a.stack_m).abs() < 1e-9 && c.parts[0].d_top.abs() < 1e-12);
            for q in c.parts.windows(2) { assert!((q[0].d_base - q[1].d_top).abs() < 1e-9 && q[1].d_base > q[1].d_top); }
            if c.parts.len() > a.parts.len() { split += 1; }
            let la = LobeAt { col: b, ed: &ed, ed_off: 0.0, vl: Some(&on) };
            for dz in 0..40 { if let Some(cd) = la.pyro_code(-(dz as f64) * 0.5, x, y) { codes.insert(cd); } }
        }
        assert!(n > 50 && split > 5, "columns {n} split {split}");
        assert!(codes.contains(&crate::pyro::SUB_NONE) && codes.iter().any(|&c| (1..=7).contains(&c)) && !codes.contains(&crate::pyro::SUB_UNRESOLVED), "{codes:?}");
    }
    #[test]
    fn seq_fill_thickens_one_thick_unit_and_keeps_the_stack() {
        // lane PYRO2 de_seq_fill: same stack base and contiguous parts; at most one unit is thicker than deposited, and it is never a
        // unit thinner than half the gap when a thicker one is present
        let s = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let mut lp = LobeParams::default();
        lp.seq = Some(SeqParams::new(1));
        lp.cell_m = 512.0;
        let on_ed = |x: f64, y: f64| x.hypot(y) < 2500.0;
        let zf = |x: f64, y: f64| 500.0 - 0.15 * x.hypot(y) + 4.0 * (x / 300.0).sin();
        let style_at = |_x: f64, _y: f64| None::<Style>;
        let dom = Some((-3000.0, -3000.0, 3000.0, 3000.0));
        let bbox = (-1500.0, -1500.0, 1500.0, 1500.0);
        let off = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        lp.seq.as_mut().unwrap().fill = Some(SeqFillParams::default());
        let on = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        let (mut n, mut moved) = (0, 0);
        for k in 0..400 {
            let (x, y) = (-1400.0 + 7.0 * k as f64, -900.0 + 4.5 * k as f64);
            let (a, b) = (off.column(x, y, 100.0), on.column(x, y, 100.0));
            assert!((a.stack_m - b.stack_m).abs() < 1e-9);
            if b.parts.is_empty() { continue; }
            n += 1;
            assert!(b.parts[0].d_top.abs() < 1e-12);
            for q in b.parts.windows(2) { assert!((q[0].d_base - q[1].d_top).abs() < 1e-9 && q[1].d_base >= q[1].d_top); }
            let mut uis: Vec<u32> = b.parts.iter().map(|p| p.ui).collect(); uis.dedup();
            let ext: Vec<(u32, f64, f64)> = uis.iter().map(|&ui| { let t: f64 = b.parts.iter().filter(|p| p.ui == ui).map(|p| p.d_base - p.d_top).sum();
                let tu = b.parts.iter().find(|p| p.ui == ui).unwrap().tu; (ui, t, tu) }).filter(|q| on.units[q.0 as usize].age >= RECORDED_AGE && q.1 > q.2 + 1e-6).collect();
            assert!(ext.len() <= 1, "{ext:?}");
            if let Some(&(_, t, tu)) = ext.first() {
                let gap = t - tu;
                let thickest = uis.iter().map(|&ui| b.parts.iter().find(|p| p.ui == ui).unwrap().tu).fold(0.0f64, f64::max);
                assert!(tu >= 0.5 * gap - 1e-9 || (tu - thickest).abs() < 1e-9, "tu {tu} gap {gap} thickest {thickest}");
                if a.parts.iter().map(|p| p.ui).collect::<Vec<_>>() != b.parts.iter().map(|p| p.ui).collect::<Vec<_>>() || a.parts[0].d_base != b.parts[0].d_base { moved += 1; }
            }
        }
        assert!(n > 50, "columns {n} moved {moved}");
    }
    #[test]
    fn seq_level_keeps_order_base_and_levels_contacts() {
        // lane ARCH de_seq_level: off by default (params JSON unchanged); on: same stack base, contiguous parts from 0, older units
        // deeper, and the contacts between keyed units are smoother laterally than the sharp-lens stack (no inherited domes)
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap(); p.set("de", "1").unwrap(); p.set("de_lobes", "1").unwrap();
        assert!(p.set("de_seq_level", "1").is_err(), "needs de_seq");
        p.set("de_seq", "2").unwrap();
        assert!(p.set("de_seq_level", "1").is_err(), "needs de_seq=1");
        p.set("de_seq", "1").unwrap();
        assert!(!serde_json::to_string(p.de.as_ref().unwrap().lobes.as_ref().unwrap()).unwrap().contains("\"level\""));
        p.set("de_seq_level", "2").unwrap();
        assert_eq!(p.de.as_ref().unwrap().lobes.as_ref().unwrap().seq.as_ref().unwrap().level.as_ref().unwrap().mode, 2);
        assert!(p.set("de_seq_level", "3").is_err());
        let s = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let mut lp = LobeParams::default();
        lp.seq = Some(SeqParams::new(1));
        lp.cell_m = 512.0;
        let on_ed = |x: f64, y: f64| x.hypot(y) < 2500.0;
        let zf = |x: f64, y: f64| 500.0 - 0.15 * x.hypot(y) + 4.0 * (x / 300.0).sin();
        let style_at = |_x: f64, _y: f64| None::<Style>;
        let dom = Some((-3000.0, -3000.0, 3000.0, 3000.0));
        let bbox = (-1500.0, -1500.0, 1500.0, 1500.0);
        let off = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        lp.seq.as_mut().unwrap().level = Some(SeqLevelParams::new(1));
        let on = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        // a second block bbox: the same column at a shared post (pure per-post function)
        let on2 = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], (-200.0, -1500.0, 2800.0, 1500.0), &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        let (mut n, mut diff) = (0, 0);
        // a line of posts 2 m apart: lateral relief of the keyed contacts (depth change between neighbouring posts per contact)
        let (mut rel_off, mut rel_on, mut nc_off, mut nc_on) = (0.0f64, 0.0f64, 0usize, 0usize);
        let mut prev: Option<(Vec<(u32, f64)>, Vec<(u32, f64)>)> = None;
        for k in 0..1500 {
            let (x, y) = (-1400.0 + 2.0 * k as f64 * 0.8, -900.0 + 2.0 * k as f64 * 0.6);
            let (a, b) = (off.column(x, y, 100.0), on.column(x, y, 100.0));
            if (-200.0..1500.0).contains(&x) {
                let c = on2.column(x, y, 100.0);
                assert_eq!(b.parts.len(), c.parts.len());
                for (p, q) in b.parts.iter().zip(c.parts.iter()) { assert!(p.d_top.to_bits() == q.d_top.to_bits() && p.d_base.to_bits() == q.d_base.to_bits() && p.unit == q.unit); }
            }
            let bases = |c: &LobeCol| -> Vec<(u32, f64)> { let mut v: Vec<(u32, f64)> = Vec::new(); for p in &c.parts { match v.last_mut() { Some(q) if q.0 == p.ui => q.1 = p.d_base, _ => v.push((p.ui, p.d_base)) } } v };
            let (ba, bb) = (bases(&a), bases(&b));
            if let Some((pa, pb)) = prev.as_ref() {
                for (ui, d) in &ba { if let Some(q) = pa.iter().find(|q| q.0 == *ui) { rel_off += (d - q.1).abs(); nc_off += 1; } }
                for (ui, d) in &bb { if let Some(q) = pb.iter().find(|q| q.0 == *ui) { rel_on += (d - q.1).abs(); nc_on += 1; } }
            }
            prev = Some((ba, bb));
            if b.parts.is_empty() { continue; }
            n += 1;
            assert!((a.stack_m - b.stack_m).abs() < 1e-9 || a.parts.is_empty(), "stack {} vs {}", a.stack_m, b.stack_m);
            assert!(b.parts[0].d_top.abs() < 1e-12);
            for q in b.parts.windows(2) { assert!((q[0].d_base - q[1].d_top).abs() < 1e-9 && q[1].d_base > q[1].d_top); }
            let mut uis: Vec<u32> = b.parts.iter().map(|p| p.ui).collect(); uis.dedup();
            let mut seen = std::collections::BTreeSet::new();
            for w in uis.windows(2) { assert!(on.units[w[0] as usize].age <= on.units[w[1] as usize].age, "older unit above a younger one"); }
            for u in &uis { assert!(seen.insert(*u), "unit split in two"); }
            if a.parts.iter().map(|p| p.d_base.to_bits()).collect::<Vec<_>>() != b.parts.iter().map(|p| p.d_base.to_bits()).collect::<Vec<_>>() { diff += 1; }
        }
        let (ro, rn) = (rel_off / nc_off.max(1) as f64, rel_on / nc_on.max(1) as f64);
        assert!(n > 300 && diff > 100, "columns {n} changed {diff}");
        assert!(rn < ro, "mean contact relief per 2 m step: levelled {rn} vs sharp stack {ro}");
    }
    #[test]
    fn seq_stretch_scales_the_keyed_pile_by_one_factor() {
        // lane ARCH de_seq_stretch: same stack base, contiguous parts from 0, same unit order as the levelled column, and every
        // keyed unit at a post thickened by the same factor (>= 1)
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap(); p.set("de", "1").unwrap(); p.set("de_lobes", "1").unwrap(); p.set("de_seq", "1").unwrap();
        assert!(p.set("de_seq_stretch", "1").is_err(), "needs de_seq_level");
        p.set("de_seq_level", "1").unwrap();
        assert!(!serde_json::to_string(p.de.as_ref().unwrap().lobes.as_ref().unwrap()).unwrap().contains("stretch"));
        p.set("de_seq_stretch", "1").unwrap();
        let s = Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap();
        let mut lp = LobeParams::default();
        lp.seq = Some(SeqParams::new(1));
        lp.cell_m = 512.0;
        lp.seq.as_mut().unwrap().level = Some(SeqLevelParams::new(1));
        let on_ed = |x: f64, y: f64| x.hypot(y) < 2500.0;
        let zf = |x: f64, y: f64| 500.0 - 0.15 * x.hypot(y) + 4.0 * (x / 300.0).sin();
        let style_at = |_x: f64, _y: f64| None::<Style>;
        let dom = Some((-3000.0, -3000.0, 3000.0, 3000.0));
        let bbox = (-1500.0, -1500.0, 1500.0, 1500.0);
        let lv = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        lp.seq.as_mut().unwrap().level.as_mut().unwrap().stretch = Some(SEQ_STRETCH_NOTE.into());
        let st = VolcLobes::build(&lp, STYLE_ARC, "v", [0.0, 0.0], &[], bbox, &s, &on_ed, &zf, &style_at, None, dom, None, 0.0).unwrap();
        let (mut n, mut stretched) = (0, 0);
        for k in 0..600 {
            let (x, y) = (-1400.0 + 4.7 * k as f64, -900.0 + 3.1 * k as f64);
            let (a, b) = (lv.column(x, y, 100.0), st.column(x, y, 100.0));
            assert!((a.stack_m - b.stack_m).abs() < 1e-9);
            if b.parts.is_empty() { continue; }
            n += 1;
            assert!(b.parts[0].d_top.abs() < 1e-9);
            for q in b.parts.windows(2) { assert!((q[0].d_base - q[1].d_top).abs() < 1e-9 && q[1].d_base > q[1].d_top); }
            let th = |c: &LobeCol| -> Vec<(u32, f64)> { let mut v: Vec<(u32, f64)> = Vec::new(); for p in &c.parts { match v.last_mut() { Some(q) if q.0 == p.ui => q.1 += p.d_base - p.d_top, _ => v.push((p.ui, p.d_base - p.d_top)) } } v };
            let (ta, tb) = (th(&a), th(&b));
            assert_eq!(ta.iter().map(|q| q.0).collect::<Vec<_>>(), tb.iter().map(|q| q.0).collect::<Vec<_>>());
            // the youngest (top) unit of the levelled column carries its gap extension; compare the units below it
            let f: Vec<f64> = ta.iter().zip(tb.iter()).skip(1).filter(|(p, _)| p.1 > 1e-6).map(|(p, q)| q.1 / p.1).collect();
            if let Some(f0) = f.first() {
                assert!(*f0 >= 1.0 - 1e-9 && f.iter().all(|v| (v - f0).abs() < 1e-6 * f0.max(1.0)), "factors {f:?}");
                if *f0 > 1.0 + 1e-6 { stretched += 1; }
            }
        }
        assert!(n > 100 && stretched > 20, "columns {n} stretched {stretched}");
    }
    #[test]
    fn resample_keeps_ends_and_spacing() {
        let p = vec![[0.0, 0.0, 0.0], [128.0, 0.0, 128.0], [256.0, 50.0, 265.0]];
        let r = resample(&p, 16.0);
        assert_eq!(r[0], p[0]);
        let l = r.last().unwrap();
        assert!((l[0] - 256.0).abs() < 1e-9 && (l[1] - 50.0).abs() < 1e-9 && (l[2] - 265.0).abs() < 1e-9);
        for w in r.windows(2) { assert!((w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1]) < 20.0); }
    }
}
