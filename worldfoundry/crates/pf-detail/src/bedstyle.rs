//! Lane BS (bed styles from the record; node ROCK-LOOK-PROVENANCE follow-up, 2026-10-15). Opt-in `beds=2` (needs version >= 2).
//! PREREG: `unpublished-design-record`.
//!
//! design 2026-10-15: the v2/v3 keyed hard/soft bed pairs read "a bit too regular"; "the different look should be attached to
//! the providence"; "geology-derived generation". This module replaces the pairs (v2 V5) and packages (v3 W1) by class-specific
//! bed styles chosen by each layer's RECORD (lith, env, max burial, flags, age span, thickness) and, for volcanic edifices,
//! by the volcano record (kind, cone height, onset and last eruption). What is drawn is still keyed: bed and flow-unit
//! thicknesses and per-bed strengths are a STATISTICAL sub-resolution closure of the recorded unit (keyed draws, purpose
//! `t2.closure.bedstyle.v1`). Layer contacts, order, dip and the class rock-mass strength (thickness-weighted mean of 1/Sc,
//! conserved per layer) are the model's.
//!
//! Pipeline per layer: fine beds (style) -> contact zone at a truncated top -> burial contrast -> merge of beds thinner than
//! `rmin_m` into their more similar neighbour (the 2 m grid cannot show thinner beds; merging uses the thickness-weighted
//! mean of 1/m, so it conserves) -> multiplier clamp and per-layer rescale (exact conservation).
//!
//! All style numbers are project assumptions from recalled textbook ranges (lane PV `rock_look_map.json` sources: boggs2012,
//! tucker_wright1990, flugel2010, potter2005, walker1960, ghosh2006, gill1981, macdonald1953, winter2010), unverified.

use serde::Serialize;
use vz_keyed_rng::Scope;
use crate::{k_of, sc_of, R};

pub const PURPOSE_BS: &str = "t2.closure.bedstyle.v1";
pub const BED_STYLE_NOTE: &str = "lane BS beds=2: bed and flow-unit thicknesses and per-bed strengths are a STATISTICAL sub-resolution closure (keyed draws) of the recorded unit: the style is chosen by the layer record (lith, env), its contrast by the recorded max burial, a contact zone by the truncated-top flag, interflow weathered horizons of lava stacks by the recorded accumulation rate; edifices take one sequence per volcano by its recorded kind, height and age span. Layer contacts, order, dip and the class rock-mass strength (thickness-weighted mean of 1/Sc per layer) are the model's. Not simulated deposition or eruption history; no newly simulated events.";

/// LithologyCodeV1 codes used here.
const L_BASALT: u8 = 3;
const L_CARB: u8 = 4;
const L_SHALE: u8 = 5;
const L_SAND: u8 = 6;
const L_CONGL: u8 = 7;
const L_EVAP: u8 = 8;
const L_COAL: u8 = 9;
/// pf-geodynamics stratigraphy env codes.
const E_FLUVIAL: u8 = 3;
const E_LACUSTRINE: u8 = 4;
const E_DELTAIC: u8 = 5;
const E_DEEP: u8 = 7;
/// Record flag bits.
const F_TRUNC_TOP: u8 = 1;

/// Lane BS2 facies codes (bed_style_id.u8): the dominant fine-bed kind of a resolved bed.
pub const FAC_NONE: u8 = 0;
pub const FAC_LAVA_CORE: u8 = 1;
pub const FAC_LAVA_TOP: u8 = 2;
pub const FAC_LAVA_BASE: u8 = 3;
pub const FAC_BOLE: u8 = 4;
pub const FAC_PYRO: u8 = 5;
pub const FAC_PYRO_HARD: u8 = 6;
pub const FAC_LS_CLEAN: u8 = 7;
pub const FAC_LS_MARLY: u8 = 8;
pub const FAC_SHALE: u8 = 9;
pub const FAC_RIB: u8 = 10;
pub const FAC_SAND: u8 = 11;
pub const FAC_MUD: u8 = 12;
pub const FAC_HETERO: u8 = 13;
pub const FAC_CONGL: u8 = 14;
pub const FAC_EVAP: u8 = 15;
pub const FAC_EVAP_PART: u8 = 16;
pub const FAC_COAL: u8 = 17;
pub const FAC_CM_SAND: u8 = 18;
pub const FAC_CM_SHALE: u8 = 19;
pub const FAC_LS_RHYTHM: u8 = 20;
pub const FAC_MARL_PART: u8 = 21;
pub const FAC_CONTACT_BOLE: u8 = 22;
pub const FAC_CONTACT_KARST: u8 = 23;
pub const FAC_CONTACT_PALEOSOL: u8 = 24;
pub const FAC_DEPOSIT: u8 = 254;
pub const FAC_WATER: u8 = 255;
pub const FACIES_NOTE: &str = "lane BS2 bed_meta=1: bed_style_id u8 = facies of the exposed resolved bed (the thicker constituent of a merged bed; beds=2 closure, see bed_style_note): 0 unbedded / not bedded, 1 lava dense core, 2 lava rubbly / vesicular flow top, 3 lava basal breccia, 4 interflow bole / paleosol, 5 pyroclastic / lahar (weak), 6 pyroclastic indurated (welded / lahar), 7 limestone clean beds, 8 limestone marly nodular, 9 shale matrix, 10 silt / sand rib, 11 sandstone body, 12 mud interbed, 13 deltaic heterolithic, 14 conglomerate, 15 evaporite, 16 evaporite mud parting, 17 coal seam, 18 coal-measure sandstone, 19 coal-measure shale, 20 limestone rhythmite, 21 marl parting, 22 red bole at a truncated lava top (recorded flag), 23 paleokarst at a truncated carbonate top (recorded flag), 24 paleosol at a truncated top (recorded flag), 254 deposit, 255 water / none";
pub const UNIT_NOTE: &str = "lane BS2 bed_meta=1: flow_unit u16 = index of the generated unit (lava flow unit with its top bole, pyroclastic package, limestone cycle, shale package, sandstone body, ...) containing the exposed point, counted from the top of its layer (strata) or of the edifice sequence below the datum (0 = top); units of the sequence continued above the column datum have bit 0x8000 set; 0xFFFF none (unbedded, deposit, water). Units are a keyed closure, not recorded events";

#[derive(Clone, Debug, Serialize)]
pub struct BedStyleParams {
    pub version: u32,
    /// Beds thinner than this are merged into their more similar neighbour (grid resolution).
    pub rmin_m: f64,
    /// Final multiplier clamp (Sc = class Sc x m) before the conserving rescale.
    pub m_min: f64,
    pub m_max: f64,
    /// Interflow weathered horizon: p = pmax (1 - exp(-repose / tau)); repose = median unit thickness / accumulation rate.
    pub tau_soil_ky: f64,
    pub p_interflow_max: f64,
    /// Floor of a recorded age span (the record's time step is about 2 My; a zero span is unresolved, not instantaneous).
    pub span_floor_my: f64,
    /// p of an interflow horizon when no span is recorded (edifice without T0 volcano records).
    pub p_interflow_unknown: f64,
    /// Burial bands (m) as in lane PV's look map: shallow < 500, moderate < 2000, deep < 5000, very deep beyond.
    pub burial_bands_m: [f64; 3],
    /// Lane BS2 S1 (opt-in bs_csets=1): limestone cycle sets. None = off (params JSON unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub csets: Option<CsetParams>,
    /// Lane BS2 S2 (opt-in bed_meta=1): publish bed_style_id / flow_unit. Output only.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bed_meta: Option<bool>,
    /// Lane BS2 S3 (opt-in ed_contacts=1): edifice depth_in_unit to flow-unit contacts. Output only.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ed_contacts: Option<bool>,
    /// Lane DATUM (opt-in ed_datum=1): the edifice datum S and the T1b term of the lobe construction surface D are taken from
    /// T1b's surface before the coast stage's wave cut (window field uncut_mm where published, else height_mm). None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ed_datum: Option<String>,
    /// Lane ARCH RIDGE (opt-in ed_env=<h>, needs ed_datum=1): the datum surface is the interfluve envelope of the pre-cut surface
    /// (ED_ENV_NOTE). None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ed_env: Option<EdEnvParams>,
}

pub const ED_ENV_NOTE: &str = "lane ARCH ed_env=L (needs ed_datum=1): the surface the edifice datum S and the T1b term of the lobe construction surface D are computed from is the interfluve envelope of T1b's pre-wave-cut surface instead of that surface itself. For every land post (pre-cut height > 0.5 m) a morphological closing is taken along the local regional strike: the direction along the contours of the regional surface (box mean of the pre-cut surface over (2 t + 1)^2 posts, t = strike_half_posts, default 16 = 1056 m), in which a flank is level. Running maximum over samples every 32 m within +-L posts along that line (land samples only), then running minimum of that over the same line (never below the surface). Valleys and coves running down the flank and narrower than the line window are bridged at the level of their lower bounding crest; ridge crests, the down-flank slope and anything wider than the window are kept. This is the classic reconstruction of an eroded volcano's constructional surface from its planeze / interfluve remnants: the pile was built before the valleys were cut, so its units and beds run under the valleys and are truncated by them instead of draping down into them; valley floors then expose older units. Opt-in ed_env_route=1 also traces the de_seq keyed units on the box mean of the envelope; without it routing is unchanged. Assumptions: interfluve crests are little lowered; bridged lows narrower than the window are erosional (a constructional low between flow fields is filled too). A pure function of the window (block-independent). No draws.";

#[derive(Clone, Debug, Serialize)]
pub struct EdEnvParams {
    pub version: u32, pub half_posts: i64, pub strike_half_posts: i64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub route: Option<bool>,
    pub note: String,
}
impl EdEnvParams { pub fn new(h: i64) -> Self { EdEnvParams { version: 1, half_posts: h, strike_half_posts: 16, route: None, note: ED_ENV_NOTE.into() } } }

/// Lane ARCH RIDGE: separable box mean over (2 t + 1)^2 posts of an n x n lattice (edges clamped), running sums.
pub fn box_mean_sep(f: &[f64], n: usize, t: i64) -> Vec<f64> {
    let ni = n as i64;
    let pass = |a: &[f64], horiz: bool| -> Vec<f64> {
        let mut o = vec![0.0f64; n * n];
        for l in 0..ni {
            let at = |k: i64| -> f64 { let k = k.clamp(0, ni - 1); if horiz { a[(l * ni + k) as usize] } else { a[(k * ni + l) as usize] } };
            let mut s: f64 = (-t..=t).map(at).sum();
            for k in 0..ni {
                let idx = if horiz { l * ni + k } else { k * ni + l };
                o[idx as usize] = s / (2 * t + 1) as f64;
                s += at(k + t + 1) - at(k - t);
            }
        }
        o
    };
    pass(&pass(f, true), false)
}

/// Lane ARCH RIDGE: interfluve envelope (ED_ENV_NOTE) of an n x n lattice f (row-major, unit post spacing in the lattice):
/// closing along the strike of the regional surface (box mean over (2 t + 1)^2 posts) over +-half posts, land (> 0.5) only.
pub fn envelope_strike(f: &[f64], n: usize, half: i64, t: i64) -> Vec<f64> {
    let ni = n as i64;
    let reg = box_mean_sep(f, n, t);
    let bil = |a: &[f64], fx: f64, fy: f64| -> Option<f64> {
        if fx < 0.0 || fy < 0.0 || fx > (ni - 1) as f64 || fy > (ni - 1) as f64 { return None; }
        let (i, j) = ((fx.floor() as i64).min(ni - 2), (fy.floor() as i64).min(ni - 2));
        let (u, v) = (fx - i as f64, fy - j as f64);
        let k = (j * ni + i) as usize;
        Some(a[k] * (1.0 - u) * (1.0 - v) + a[k + 1] * u * (1.0 - v) + a[k + n] * (1.0 - u) * v + a[k + n + 1] * u * v)
    };
    // strike direction per post (unit, lattice units); None on sea posts or where the regional surface is flat
    let dir: Vec<Option<(f64, f64)>> = (0..n * n).map(|k| {
        if f[k] <= 0.5 { return None; }
        let (i, j) = ((k % n) as i64, (k / n) as i64);
        let g = |di: i64, dj: i64| reg[((j + dj).clamp(0, ni - 1) * ni + (i + di).clamp(0, ni - 1)) as usize];
        let (gx, gy) = ((g(1, 0) - g(-1, 0)) * 0.5, (g(0, 1) - g(0, -1)) * 0.5);
        let m = gx.hypot(gy);
        if m < 1e-6 { None } else { Some((-gy / m, gx / m)) }
    }).collect();
    let pass = |a: &[f64], take_max: bool| -> Vec<f64> {
        let mut o = a.to_vec();
        for k in 0..n * n {
            let Some((sx, sy)) = dir[k] else { continue };
            let (i, j) = ((k % n) as f64, (k / n) as f64);
            let mut acc = a[k];
            for q in -half..=half {
                if q == 0 { continue; }
                let (fx, fy) = (i + q as f64 * sx, j + q as f64 * sy);
                match (bil(f, fx, fy), bil(a, fx, fy)) { (Some(h), Some(z)) if h > 0.5 => { acc = if take_max { acc.max(z) } else { acc.min(z) }; } _ => {} }
            }
            o[k] = acc;
        }
        o
    };
    let d = pass(f, true);
    pass(&d, false).iter().zip(f.iter()).map(|(a, b)| a.max(*b)).collect()
}

pub const ED_DATUM_NOTE: &str = "lane DATUM ed_datum=1: the edifice bedding datum S (288 m box mean, bicubic) and the T1b term of the lobe construction surface D = max(S, T1b bicubic) are computed from T1b's pre-wave-cut surface (uncut_mm on the cliff-class posts, the published height elsewhere) instead of the cut surface, so edifice beds and flow units are truncated by the sea cliff (they predate the wave cut) instead of bending down into it. Fluvial valleys (T1b relaxation, young radial gullies) are not restored. Needs a window publishing uncut_mm; without that field the key changes nothing.";

/// Lane BS2 S1: limestone cycle sets (keyed 10-30 m modulation of marl fraction and clean-bed / marl strength).
#[derive(Clone, Debug, Serialize)]
pub struct CsetParams {
    /// Set thickness: lognormal median, sigma, clamp (m).
    pub set_m: [f64; 4],
    /// Cycle thickness inside a set: lognormal median, sigma, clamp (m).
    pub cycle_m: [f64; 4],
    /// Marly-base fraction of a cycle: base + slope * q (q = set marliness, uniform 0-1), clamp [0.1, 0.85].
    pub marl_frac: [f64; 2],
    /// ln m of clean beds and of the marly base: a - b q (before the burial contrast exponent).
    pub clean_lnm: [f64; 2],
    pub marl_lnm: [f64; 2],
    /// Set amplitude by env: shelf (and unrecorded), deep-marine rhythmite bundles.
    pub amp_shelf: f64,
    pub amp_deep: f64,
    /// Multiplier clamp for cycle-set limestone layers (v3 packages reached about 0.15-4 x class; beds=2 default 0.2-2.5).
    pub m_range: [f64; 2],
}

impl Default for CsetParams {
    fn default() -> Self {
        CsetParams { set_m: [18.0, 0.4, 10.0, 30.0], cycle_m: [5.0, 0.6, 1.5, 15.0], marl_frac: [0.15, 0.6], clean_lnm: [0.85, 0.8], marl_lnm: [-1.05, 0.6],
                     amp_shelf: 1.0, amp_deep: 0.5, m_range: [0.15, 3.5] }
    }
}

impl Default for BedStyleParams {
    fn default() -> Self {
        BedStyleParams { version: 1, rmin_m: 1.5, m_min: 0.2, m_max: 2.5, tau_soil_ky: 30.0, p_interflow_max: 0.5, span_floor_my: 1.0,
                         p_interflow_unknown: 0.15, burial_bands_m: [500.0, 2000.0, 5000.0], csets: None, bed_meta: None, ed_contacts: None, ed_datum: None, ed_env: None }
    }
}

/// The record of one strata layer (lane PV `strata_row` fields; None = not recorded).
#[derive(Clone, Debug, Default)]
pub struct LayerRec {
    pub lith: u8,
    pub env: Option<u8>,
    pub burial_m: Option<f64>,
    pub flags: u8,
    pub age_base_my: Option<f64>,
    pub age_top_my: Option<f64>,
    pub thickness_m: f64,
}

impl LayerRec {
    pub fn from_row(v: &serde_json::Value) -> LayerRec {
        let f = |k: &str| v.get(k).and_then(serde_json::Value::as_f64);
        LayerRec { lith: v.get("lith").and_then(serde_json::Value::as_u64).unwrap_or(255) as u8, env: v.get("env").and_then(serde_json::Value::as_u64).map(|e| e as u8),
                   burial_m: f("max_burial_m"), flags: v.get("flags").and_then(serde_json::Value::as_u64).unwrap_or(0) as u8,
                   age_base_my: f("age_base_my"), age_top_my: f("age_top_my"), thickness_m: f("thickness_m").unwrap_or(0.0) }
    }
    fn span_my(&self) -> Option<f64> { match (self.age_base_my, self.age_top_my) { (Some(a), Some(b)) => Some((b - a).abs()), _ => None } }
}

/// The record of one volcano (T1b window list + T0 volcano record).
#[derive(Clone, Debug, Default, Serialize)]
pub struct VolcRec {
    pub id: String,
    pub kind: String,
    pub height_m: Option<f64>,
    pub onset_my: Option<f64>,
    pub last_eruption_my: Option<f64>,
}

/// Keyed draw stream of one scope (sequential counters).
struct Dr<'a> { sk: &'a Scope, n: u64 }
impl<'a> Dr<'a> {
    fn z(&mut self) -> R<f64> { self.n += 1; Ok(self.sk.draw(PURPOSE_BS, self.n, 0).map_err(|e| e.0)?.normal_f64()) }
    fn u(&mut self) -> R<f64> { self.n += 1; Ok(self.sk.draw(PURPOSE_BS, self.n, 0).map_err(|e| e.0)?.uniform_f64()) }
    fn ln(&mut self, med: f64, sig: f64, lo: f64, hi: f64) -> R<f64> { Ok((med * (sig * self.z()?).exp()).clamp(lo, hi)) }
    fn m(&mut self, mu: f64, sig: f64) -> R<f64> { Ok((mu + sig * self.z()?).exp()) }
    fn chance(&mut self, p: f64) -> R<bool> { Ok(self.u()? < p) }
    fn range(&mut self, a: f64, b: f64) -> R<f64> { Ok(a + (b - a) * self.u()?) }
}

/// Lava sequence settings (a lava stack layer or a volcano kind).
#[derive(Clone, Copy, Debug)]
pub struct LavaCfg {
    /// Flow unit thickness: lognormal median, sigma, clamp.
    pub tm: f64, pub sig: f64, pub lo: f64, pub hi: f64,
    /// Fractions of a unit that are rubbly / vesicular top and basal breccia (uniform ranges).
    pub top: (f64, f64), pub base: (f64, f64),
    /// ln m of the dense core, of the rubbly top and of the basal breccia (mean, sigma per unit).
    pub core: (f64, f64), pub topm: (f64, f64), pub basem: (f64, f64),
    /// Share of units that are pyroclastic / lahar packages instead of lava, and their thickness.
    pub p_pyro: f64, pub pyro: (f64, f64, f64, f64),
}

pub const LAVA_STACK: LavaCfg = LavaCfg { tm: 6.0, sig: 0.8, lo: 1.0, hi: 30.0, top: (0.15, 0.40), base: (0.03, 0.10), core: (0.35, 0.25), topm: (-0.6, 0.25),
    basem: (-0.35, 0.2), p_pyro: 0.08, pyro: (1.5, 0.7, 0.3, 8.0) };
pub const LAVA_ARC: LavaCfg = LavaCfg { tm: 8.0, sig: 0.8, lo: 2.0, hi: 50.0, top: (0.15, 0.40), base: (0.05, 0.15), core: (0.35, 0.25), topm: (-0.6, 0.25),
    basem: (-0.35, 0.2), p_pyro: 0.45, pyro: (5.0, 0.9, 0.5, 40.0) };
pub const LAVA_SHIELD: LavaCfg = LavaCfg { tm: 2.5, sig: 0.7, lo: 0.5, hi: 10.0, top: (0.05, 0.25), base: (0.02, 0.10), core: (0.3, 0.25), topm: (-0.5, 0.25),
    basem: (-0.3, 0.2), p_pyro: 0.05, pyro: (0.8, 0.6, 0.2, 3.0) };
pub const LAVA_LIP: LavaCfg = LavaCfg { tm: 20.0, sig: 0.6, lo: 5.0, hi: 80.0, top: (0.08, 0.20), base: (0.02, 0.06), core: (0.45, 0.2), topm: (-0.6, 0.25),
    basem: (-0.35, 0.2), p_pyro: 0.02, pyro: (1.0, 0.6, 0.2, 4.0) };
pub const LAVA_UNKNOWN: LavaCfg = LavaCfg { p_pyro: 0.3, ..LAVA_ARC };

pub fn lava_cfg_of_kind(kind: &str) -> (&'static str, LavaCfg) {
    match kind {
        "continental_arc" | "island_arc" => ("arc_stratovolcano", LAVA_ARC),
        "hotspot" | "rift" => ("shield", LAVA_SHIELD),
        "large_igneous_province" => ("flood_basalt", LAVA_LIP),
        _ => ("composite_unknown_kind", LAVA_UNKNOWN),
    }
}

fn band(burial: Option<f64>, b: &[f64; 3]) -> usize {
    match burial { None => 1, Some(z) => if z < b[0] { 0 } else if z < b[1] { 1 } else if z < b[2] { 2 } else { 3 } }
}

/// Contrast exponent c (ln m -> c ln m) by class and burial band (shallow, moderate, deep, very deep). Lavas: secondary
/// minerals fill vesicles and breccia porosity with burial (walker1960, neuhoff1999), so contrast falls; carbonate:
/// differential cementation and pressure solution sharpen limestone/marl contrast (tucker_wright1990); shale and clastics:
/// compaction and cement homogenise toward argillite / quartzite (athy1930, bjorlykke2015). Project assumptions.
fn contrast(lith: u8, bnd: usize) -> f64 {
    let t: [f64; 4] = match lith {
        L_BASALT => [1.0, 0.95, 0.85, 0.75],
        L_CARB => [0.8, 1.0, 1.1, 1.1],
        L_SHALE => [0.8, 1.0, 0.9, 0.7],
        _ => [0.8, 1.0, 0.9, 0.7],
    };
    t[bnd]
}

/// Interflow weathered-horizon probability from an accumulation rate (m / My) and the median unit thickness (m).
fn p_interflow(p: &BedStyleParams, rate_m_per_my: Option<f64>, tm: f64) -> f64 {
    match rate_m_per_my {
        Some(r) if r > 0.0 => { let repose_ky = tm / r * 1000.0; p.p_interflow_max * (1.0 - (-repose_ky / p.tau_soil_ky).exp()) }
        _ => p.p_interflow_unknown,
    }
}

/// A fine bed: (thickness m, ln m raw, facies, unit index in its sequence).
pub type Fb = (f64, f64, u8, u16);

/// Fine beds top-down: one lava sequence of `len` m. Units: an interflow horizon belongs to the flow unit below it.
fn lava_seq(len: f64, c: &LavaCfg, p_if: f64, dr: &mut Dr, out: &mut Vec<Fb>) -> R<()> {
    let mut d = 0.0;
    let mut u: u16 = 0;
    while d < len {
        if p_if > 0.0 && dr.chance(p_if)? {
            // interflow weathered horizon (red bole / paleosol / epiclastic) on top of the unit below (ghosh2006)
            let t = dr.ln(0.8, 0.6, 0.2, 3.0)?;
            out.push((t, (0.35f64).ln() + 0.15 * dr.z()?, FAC_BOLE, u)); d += t;
        }
        if dr.chance(c.p_pyro)? {
            // pyroclastic / lahar package: ash, scoria, block-and-ash and lahar beds, mostly weak, a few indurated
            let tp = dr.ln(c.pyro.0, c.pyro.1, c.pyro.2, c.pyro.3)?;
            let mut e = 0.0;
            while e < tp {
                let t = dr.ln(0.8, 0.8, 0.05, 6.0)?.min(tp - e).max(1e-3);
                let hard = dr.chance(0.15)?;
                let lm = if hard { 0.1 + 0.2 * dr.z()? } else { -0.55 + 0.4 * dr.z()? };
                out.push((t, lm, if hard { FAC_PYRO_HARD } else { FAC_PYRO }, u)); e += t;
            }
            d += tp;
            u = u.saturating_add(1);
            continue;
        }
        let tu = dr.ln(c.tm, c.sig, c.lo, c.hi)?;
        let (ft, fb) = (dr.range(c.top.0, c.top.1)?, dr.range(c.base.0, c.base.1)?);
        let (lc, lt, lb) = (c.core.0 + c.core.1 * dr.z()?, c.topm.0 + c.topm.1 * dr.z()?, c.basem.0 + c.basem.1 * dr.z()?);
        out.push((ft * tu, lt, FAC_LAVA_TOP, u));
        out.push(((1.0 - ft - fb) * tu, lc, FAC_LAVA_CORE, u));
        out.push((fb * tu, lb, FAC_LAVA_BASE, u));
        d += tu;
        u = u.saturating_add(1);
    }
    Ok(())
}

/// Thickness-weighted harmonic mean of a group of fine beds as one bed, in raw ln m (the burial contrast c applied inside).
fn part(v: &[(f64, f64)], c: f64) -> (f64, f64) {
    let t: f64 = v.iter().map(|b| b.0).sum();
    let w: f64 = v.iter().map(|b| b.0 / (c * b.1).exp()).sum();
    (t, (t / w).ln() / c)
}

/// Fine beds of one strata layer by its record (top-down), before contrast and merge.
fn strata_seq(rec: &LayerRec, len: f64, p: &BedStyleParams, c: f64, dr: &mut Dr) -> R<(&'static str, Vec<Fb>)> {
    let mut out: Vec<Fb> = Vec::new();
    let env = rec.env.unwrap_or(255);
    let style: &'static str;
    match rec.lith {
        L_BASALT => {
            style = "lava_flow_units";
            let rate = rec.span_my().map(|s| rec.thickness_m / s.max(p.span_floor_my));
            lava_seq(len, &LAVA_STACK, p_interflow(p, rate, LAVA_STACK.tm), dr, &mut out)?;
        }
        L_CARB if p.csets.is_some() => {
            // Lane BS2 S1: cycle sets. A set (10-30 m) has one keyed marliness q that moves the marly-base fraction and the
            // strength of clean beds and marl together (marl-rich sets weather back as benches, clean sets stand as cliffs).
            // The set amplitude follows the recorded env (shelf full; deep-marine rhythmite bundles half); burial sets c.
            let cs = p.csets.as_ref().unwrap();
            let deep = env == E_DEEP;
            style = if deep { "limestone_rhythmic_deep_sets" } else { "limestone_shelf_cycle_sets" };
            let amp = if deep { cs.amp_deep } else { cs.amp_shelf };
            let mut d = 0.0;
            let mut u: u16 = 0;
            while d < len {
                let ts = dr.ln(cs.set_m[0], cs.set_m[1], cs.set_m[2], cs.set_m[3])?;
                let q = dr.u()?;
                let mut e = 0.0;
                while e < ts {
                    if deep {
                        // a bundle of couplets: limestone strength and marl share move with the set
                        let t = dr.ln(0.3, 0.4, 0.1, 1.0)?; let lmq = dr.ln(0.1 * (0.5 + q), 0.5, 0.02, 0.6)?;
                        let up = vec![(t, amp * (cs.clean_lnm[0] - cs.clean_lnm[1] * q) + 0.1 * dr.z()?), (lmq, (0.4f64).ln() + 0.1 * dr.z()?)];
                        let pr = part(&up, c); out.push((pr.0, pr.1, FAC_LS_RHYTHM, u)); e += pr.0; u = u.saturating_add(1);
                        continue;
                    }
                    let tc = dr.ln(cs.cycle_m[0], cs.cycle_m[1], cs.cycle_m[2], cs.cycle_m[3])?;
                    let fl = (cs.marl_frac[0] + cs.marl_frac[1] * q + 0.1 * dr.z()?).clamp(0.1, 0.85);
                    let pp = 0.3 + 0.5 * q;
                    let (mut up, mut lo) = (Vec::new(), Vec::new());
                    let mut f = 0.0;
                    while f < (1.0 - fl) * tc {
                        let t = dr.ln(0.6, 0.6, 0.2, 2.0)?; up.push((t, amp * (cs.clean_lnm[0] - cs.clean_lnm[1] * q) + 0.15 * dr.z()?)); f += t;
                        if dr.chance(pp)? { let tp = dr.ln(0.06, 0.6, 0.01, 0.3)?; up.push((tp, (0.3f64).ln() + 0.1 * dr.z()?)); f += tp; }
                    }
                    while f < tc {
                        let t = dr.ln(0.3, 0.5, 0.05, 1.0)?; lo.push((t, amp * (cs.marl_lnm[0] - cs.marl_lnm[1] * q) + 0.2 * dr.z()?)); f += t;
                    }
                    let a = part(&up, c); out.push((a.0, a.1, FAC_LS_CLEAN, u));
                    if !lo.is_empty() { let b = part(&lo, c); out.push((b.0, b.1, FAC_LS_MARLY, u)); }
                    e += f; u = u.saturating_add(1);
                }
                d += e;
            }
        }
        L_CARB => {
            if env == E_DEEP {
                // pelagic limestone-marl rhythmites (thin, more even beds)
                style = "limestone_rhythmic_deep";
                let mut d = 0.0;
                let mut u: u16 = 0;
                while d < len {
                    let t = dr.ln(0.3, 0.4, 0.1, 1.0)?; out.push((t, 0.2 + 0.1 * dr.z()?, FAC_LS_RHYTHM, u));
                    let q = dr.ln(0.1, 0.5, 0.02, 0.4)?; out.push((q, (0.4f64).ln() + 0.1 * dr.z()?, FAC_MARL_PART, u)); d += t + q;
                    u = u.saturating_add(1);
                }
            } else {
                // shelf (and unrecorded env): beds 0.2-2 m in shallowing-upward cycles; top-down: clean upper part, marly nodular base.
                // Revision R1 (PREREG): the beds are below the grid, so each cycle half is emitted as ONE resolved bed whose strength is
                // the thickness-weighted harmonic mean of its fine beds and partings (merging single 0.2-2 m beds gave near-uniform
                // 2-3 m bedsets, i.e. a new regular stripe). The visible spacing is the cycle's (lognormal 1-20 m).
                style = "limestone_shelf_cycles";
                let mut d = 0.0;
                let mut u: u16 = 0;
                while d < len {
                    let tc = dr.ln(4.0, 0.6, 1.0, 20.0)?;
                    let q = dr.u()?;
                    let fl = 0.2 + 0.3 * q;
                    let pp = 0.4 + 0.4 * q;
                    let (mut up, mut lo) = (Vec::new(), Vec::new());
                    let mut e = 0.0;
                    while e < (1.0 - fl) * tc {
                        let t = dr.ln(0.6, 0.6, 0.2, 2.0)?; up.push((t, 0.25 + 0.15 * dr.z()?)); e += t;
                        if dr.chance(pp)? { let tp = dr.ln(0.06, 0.6, 0.01, 0.3)?; up.push((tp, (0.35f64).ln() + 0.1 * dr.z()?)); e += tp; }
                    }
                    while e < tc {
                        let t = dr.ln(0.3, 0.5, 0.05, 1.0)?; lo.push((t, -0.35 + 0.2 * dr.z()?)); e += t;
                    }
                    let a = part(&up, c); out.push((a.0, a.1, FAC_LS_CLEAN, u));
                    if !lo.is_empty() { let b = part(&lo, c); out.push((b.0, b.1, FAC_LS_MARLY, u)); }
                    d += e;
                    u = u.saturating_add(1);
                }
            }
        }
        L_SHALE => {
            // fissile matrix (uniform at the grid scale) with Poisson-spaced silt / sand ribs; rib density varies by package
            let (s0, tr, lr): (f64, f64, f64) = match env { E_DEEP => (8.0, 0.1, 0.7), E_LACUSTRINE => (10.0, 0.08, 0.5), _ => (3.0, 0.25, 0.9) };
            style = match env { E_DEEP => "shale_laminated_deep", E_LACUSTRINE => "shale_varved_lacustrine", _ => "shale_shelf_ribbed" };
            let mut d = 0.0;
            let mut u: u16 = 0;
            while d < len {
                let tp = dr.ln(15.0, 0.6, 3.0, 60.0)?;
                let rho = dr.m(0.0, 0.7)?;
                let mut e = 0.0;
                while e < tp {
                    let gap = -(1.0 - dr.u()?).ln() * s0 / rho;
                    let g = gap.min(tp - e).max(0.01);
                    out.push((g, 0.05 * dr.z()?, FAC_SHALE, u)); e += g;
                    if e >= tp { break; }
                    let t = dr.ln(tr, 0.6, 0.02, 1.5)?; out.push((t, lr + 0.2 * dr.z()?, FAC_RIB, u)); e += t;
                }
                d += e;
                u = u.saturating_add(1);
            }
        }
        L_SAND => {
            if env == E_DELTAIC {
                style = "sandstone_deltaic_coarsening_up";
                let mut d = 0.0;
                let mut u: u16 = 0;
                while d < len {
                    let tp = dr.ln(8.0, 0.5, 2.0, 30.0)?;
                    let fs = dr.range(0.55, 0.75)?;
                    let mut e = 0.0;
                    while e < fs * tp { let t = dr.ln(1.0, 0.5, 0.2, 4.0)?; out.push((t, 0.25 + 0.15 * dr.z()?, FAC_SAND, u)); e += t; }
                    while e < tp { let t = dr.ln(0.8, 0.5, 0.1, 3.0)?; out.push((t, (0.65f64).ln() + 0.15 * dr.z()?, FAC_HETERO, u)); e += t; }
                    d += e;
                    u = u.saturating_add(1);
                }
            } else {
                style = if env == E_FLUVIAL { "sandstone_fluvial_channels" } else { "sandstone_channels_env_unrecorded" };
                let mut d = 0.0;
                let mut u: u16 = 0;
                while d < len {
                    let t = dr.ln(3.0, 0.7, 0.3, 15.0)?; out.push((t, 0.2 + 0.25 * dr.z()?, FAC_SAND, u)); d += t;
                    if dr.chance(0.5)? { let m = dr.ln(1.0, 0.7, 0.1, 5.0)?; out.push((m, (0.45f64).ln() + 0.15 * dr.z()?, FAC_MUD, u)); d += m; }
                    u = u.saturating_add(1);
                }
            }
        }
        L_CONGL => {
            style = "conglomerate_lenticular";
            let mut d = 0.0;
            let mut u: u16 = 0;
            while d < len {
                let t = dr.ln(2.0, 0.7, 0.5, 8.0)?; out.push((t, 0.15 + 0.3 * dr.z()?, FAC_CONGL, u)); d += t;
                if dr.chance(0.4)? { let m = dr.ln(0.6, 0.6, 0.1, 3.0)?; out.push((m, (0.7f64).ln() + 0.15 * dr.z()?, FAC_SAND, u)); d += m; }
                u = u.saturating_add(1);
            }
        }
        L_EVAP => {
            style = "evaporite_bedded";
            let mut d = 0.0;
            let mut u: u16 = 0;
            while d < len {
                let t = dr.ln(1.0, 0.7, 0.05, 5.0)?; out.push((t, 0.3 * dr.z()?, FAC_EVAP, u)); d += t;
                if dr.chance(0.5)? { let m = dr.ln(0.1, 0.6, 0.02, 0.5)?; out.push((m, (0.4f64).ln() + 0.1 * dr.z()?, FAC_EVAP_PART, u)); d += m; }
                u = u.saturating_add(1);
            }
        }
        L_COAL => {
            style = "coal_measures";
            let mut d = 0.0;
            let mut u: u16 = 0;
            while d < len {
                let sand = dr.chance(0.4)?;
                let t = if sand { dr.ln(2.0, 0.6, 0.3, 8.0)? } else { dr.ln(3.0, 0.6, 0.3, 12.0)? };
                out.push((t, if sand { (1.3f64).ln() + 0.2 * dr.z()? } else { (0.8f64).ln() + 0.1 * dr.z()? }, if sand { FAC_CM_SAND } else { FAC_CM_SHALE }, u)); d += t;
                if dr.chance(0.3)? { let s = dr.ln(0.5, 0.8, 0.1, 3.0)?; out.push((s, (0.5f64).ln() + 0.1 * dr.z()?, FAC_COAL, u)); d += s; }
                u = u.saturating_add(1);
            }
        }
        _ => { style = "unbedded"; out.push((len, 0.0, FAC_NONE, u16::MAX)); }
    }
    Ok((style, out))
}

/// Merge every bed thinner than `rmin` into its neighbour with the closer ln m (tie: the thinner), conserving thickness and
/// the thickness-weighted sum of 1/m. Input and output top-down (thickness, m). The merged bed keeps the facies and unit of
/// its thicker constituent (lane BS2 metadata; the merge decisions are those of lane BS).
pub fn merge_thin(beds: &[(f64, f64)], rmin: f64) -> Vec<(f64, f64)> {
    let meta: Vec<(u8, u16)> = vec![(0, 0); beds.len()];
    merge_thin_meta(beds, &meta, rmin).into_iter().map(|b| (b.0, b.1)).collect()
}

pub fn merge_thin_meta(beds: &[(f64, f64)], meta: &[(u8, u16)], rmin: f64) -> Vec<(f64, f64, u8, u16)> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let n = beds.len();
    if n <= 1 { return beds.iter().zip(meta.iter()).map(|(b, m)| (b.0, b.1, m.0, m.1)).collect(); }
    let mut t: Vec<f64> = beds.iter().map(|b| b.0).collect();
    let mut inv: Vec<f64> = beds.iter().map(|b| b.0 / b.1).collect(); // sum of t / m
    let mut mt: Vec<(u8, u16)> = meta.to_vec();
    let mut prev: Vec<usize> = (0..n).map(|i| if i == 0 { usize::MAX } else { i - 1 }).collect();
    let mut next: Vec<usize> = (0..n).map(|i| if i + 1 == n { usize::MAX } else { i + 1 }).collect();
    let mut alive = vec![true; n];
    let mut ver = vec![0u32; n];
    let key = |x: f64| (x * 1.0e6).round() as i64;
    let mut h: BinaryHeap<Reverse<(i64, usize, u32)>> = (0..n).filter(|&i| t[i] < rmin).map(|i| Reverse((key(t[i]), i, 0))).collect();
    let lm = |t: f64, inv: f64| (t / inv).ln();
    while let Some(Reverse((_, i, v))) = h.pop() {
        if !alive[i] || v != ver[i] || t[i] >= rmin { continue; }
        let (a, b) = (prev[i], next[i]);
        let j = match (a != usize::MAX, b != usize::MAX) {
            (false, false) => break,
            (true, false) => a,
            (false, true) => b,
            (true, true) => {
                let li = lm(t[i], inv[i]);
                let (da, db) = ((lm(t[a], inv[a]) - li).abs(), (lm(t[b], inv[b]) - li).abs());
                if da < db - 1e-12 { a } else if db < da - 1e-12 { b } else if t[a] <= t[b] { a } else { b }
            }
        };
        if t[i] > t[j] { mt[j] = mt[i]; }
        t[j] += t[i]; inv[j] += inv[i];
        alive[i] = false;
        let (p_, n_) = (prev[i], next[i]);
        if p_ != usize::MAX { next[p_] = n_; }
        if n_ != usize::MAX { prev[n_] = p_; }
        ver[j] += 1;
        if t[j] < rmin { h.push(Reverse((key(t[j]), j, ver[j]))); }
    }
    (0..n).filter(|&i| alive[i]).map(|i| (t[i], t[i] / inv[i], mt[i].0, mt[i].1)).collect()
}

/// Clamp then rescale so that sum t / m = sum t exactly (the class rock-mass 1/Sc is conserved). Returns multipliers.
fn conserve(beds: &mut [(f64, f64, u8, u16)], p: &BedStyleParams) {
    for _ in 0..4 {
        for b in beds.iter_mut() { b.1 = b.1.clamp(p.m_min, p.m_max); }
        let (st, sw): (f64, f64) = beds.iter().fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.0 / b.1));
        let f = if st > 0.0 { sw / st } else { 1.0 };
        for b in beds.iter_mut() { b.1 *= f; }
    }
}

/// One finished sequence: resolved beds (thickness, m, facies, unit) and the fine unit tops (depth from the sequence top, unit).
pub struct Seq { pub beds: Vec<(f64, f64, u8, u16)>, pub units: Vec<(f64, u16)> }

/// Finish fine beds: contrast exponent on ln m, merge, conserve.
fn finish(fine: Vec<Fb>, c: f64, len: f64, p: &BedStyleParams) -> Seq {
    // trim / extend to the layer length exactly (the generators overshoot by the last bed)
    let mut v: Vec<(f64, f64)> = Vec::with_capacity(fine.len());
    let mut meta: Vec<(u8, u16)> = Vec::with_capacity(fine.len());
    let mut units: Vec<(f64, u16)> = Vec::new();
    let mut d = 0.0;
    for (t, l, f, u) in fine {
        if d >= len { break; }
        let tt = t.min(len - d);
        if tt <= 1e-9 { continue; }
        if units.last().map(|q| q.1 != u).unwrap_or(true) { units.push((d, u)); }
        v.push((tt, (c * l).exp())); meta.push((f, u)); d += tt;
    }
    if v.is_empty() { v.push((len, 1.0)); meta.push((FAC_NONE, u16::MAX)); units.push((0.0, u16::MAX)); } else if d < len { v.last_mut().unwrap().0 += len - d; }
    let mut m = merge_thin_meta(&v, &meta, p.rmin_m);
    conserve(&mut m, p);
    Seq { beds: m, units }
}

/// Output of one layer / sequence appended to a column (bed base depths, Sc, K, facies, unit; unit tops in column depth).
#[derive(Default)]
pub struct BedsOut { pub bb: Vec<f64>, pub bs: Vec<f64>, pub bk: Vec<f64>, pub bf: Vec<u8>, pub bu: Vec<u16>, pub ut: Vec<(f64, u16)> }

/// Beds of one strata layer [top, base) depth appended to `o`; returns the style name.
pub fn layer_beds(rec: &LayerRec, top: f64, base: f64, scope: &Scope, tag: &str, p: &BedStyleParams, o: &mut BedsOut) -> R<&'static str> {
    let len = base - top;
    let l = rec.lith;
    let sk = scope.address(&format!("{tag}/bs")).map_err(|e| e.0)?;
    let mut dr = Dr { sk: &sk, n: 0 };
    let c = contrast(l, band(rec.burial_m, &p.burial_bands_m));
    let (style, mut fine) = strata_seq(rec, len, p, c, &mut dr)?;
    if style == "unbedded" {
        o.bb.push(base); o.bs.push(sc_of(l)); o.bk.push(k_of(l)); o.bf.push(FAC_NONE); o.bu.push(u16::MAX); o.ut.push((top, u16::MAX));
        return Ok(style);
    }
    if rec.flags & F_TRUNC_TOP != 0 {
        // weathered top under the overlying unit: red bole (lava), paleokarst (carbonate), paleosol (others) (lane PV look map)
        let t = dr.ln(1.5, 0.5, 0.5, 4.0)?;
        let (lm, f) = match l { L_BASALT => ((0.35f64).ln(), FAC_CONTACT_BOLE), L_CARB => ((0.6f64).ln(), FAC_CONTACT_KARST), _ => ((0.7f64).ln(), FAC_CONTACT_PALEOSOL) };
        let u0 = fine.first().map(|b| b.3).unwrap_or(0);
        fine.insert(0, (t, lm, f, u0));
    }
    let sq = match (&p.csets, l) {
        (Some(cs), L_CARB) => { let mut q = p.clone(); q.m_min = cs.m_range[0]; q.m_max = cs.m_range[1]; finish(fine, c, len, &q) }
        _ => finish(fine, c, len, p),
    };
    let mut d = top;
    for (t, m, f, u) in sq.beds { d += t; o.bb.push(d); o.bs.push(sc_of(l) * m); o.bk.push(k_of(l) / m); o.bf.push(f); o.bu.push(u); }
    *o.bb.last_mut().unwrap() = base;
    for (ud, u) in sq.units { o.ut.push((top + ud, u)); }
    Ok(style)
}

/// The same style continued ABOVE the column datum (terrain higher than the T0 cell mean is still layer 0): bed tops
/// (heights above the datum), lowest first, with Sc, K, facies, unit (| 0x8000); unit tops as signed depths (negative =
/// above the datum). The generated top-down sequence is placed with its base at the datum, so a lava unit keeps its rubbly top up.
pub fn up_beds(rec: &LayerRec, height: f64, scope: &Scope, tag: &str, p: &BedStyleParams) -> R<BedsOut> {
    let mut o = BedsOut::default();
    let mut r = rec.clone();
    r.flags &= !F_TRUNC_TOP;
    layer_beds(&r, 0.0, height, scope, tag, p, &mut o)?;
    Ok(flip_up(o, height))
}

/// Top-down base depths of a sequence of `height` m -> heights above the datum of each bed's top, lowest bed first.
fn flip_up(o: BedsOut, height: f64) -> BedsOut {
    let n = o.bb.len();
    let mut r = BedsOut::default();
    for i in (0..n).rev() {
        let top_depth = if i == 0 { 0.0 } else { o.bb[i - 1] };
        r.bb.push(height - top_depth);
        r.bs.push(o.bs[i]); r.bk.push(o.bk[i]); r.bf.push(o.bf[i]); r.bu.push(if o.bu[i] == u16::MAX { u16::MAX } else { o.bu[i] | 0x8000 });
    }
    // a unit with generated top g lies at height (height - g): signed depth g - height
    r.ut = o.ut.iter().map(|&(g, u)| (g - height, if u == u16::MAX { u } else { u | 0x8000 })).collect();
    r
}

/// Edifice bed sequence of one volcano (all class 3): `depth` m below the datum, top-down. Returns the beds and the style.
pub fn edifice_beds(v: &VolcRec, depth: f64, scope: &Scope, tag: &str, p: &BedStyleParams) -> R<(BedsOut, &'static str)> {
    let (style, c) = lava_cfg_of_kind(&v.kind);
    let rate = match (v.height_m, v.onset_my, v.last_eruption_my) {
        (Some(h), Some(a), Some(b)) => Some(h / (b - a).abs().max(p.span_floor_my)),
        _ => None,
    };
    let pif = p_interflow(p, rate, c.tm);
    let sk = scope.address(&format!("{tag}/bs")).map_err(|e| e.0)?;
    let mut dr = Dr { sk: &sk, n: 0 };
    let mut fine = Vec::new();
    lava_seq(depth, &c, pif, &mut dr, &mut fine)?;
    // a young cone: unburied (shallow band) everywhere it is exposed
    let sq = finish(fine, contrast(L_BASALT, 0), depth, p);
    let mut o = BedsOut::default();
    let mut d = 0.0;
    for (t, m, f, u) in sq.beds { d += t; o.bb.push(d); o.bs.push(sc_of(L_BASALT) * m); o.bk.push(k_of(L_BASALT) / m); o.bf.push(f); o.bu.push(u); }
    *o.bb.last_mut().unwrap() = depth;
    o.ut = sq.units;
    Ok((o, style))
}

/// The edifice sequence continued above its datum (3000 m), flipped like `up_beds`.
pub fn edifice_up_beds(v: &VolcRec, height: f64, scope: &Scope, tag: &str, p: &BedStyleParams) -> R<BedsOut> {
    let (o, _) = edifice_beds(v, height, scope, tag, p)?;
    Ok(flip_up(o, height))
}

/// The interflow probability an edifice record gives (diagnostic for the block header).
pub fn edifice_p_interflow(v: &VolcRec, p: &BedStyleParams) -> f64 {
    let (_, c) = lava_cfg_of_kind(&v.kind);
    let rate = match (v.height_m, v.onset_my, v.last_eruption_my) { (Some(h), Some(a), Some(b)) => Some(h / (b - a).abs().max(p.span_floor_my)), _ => None };
    p_interflow(p, rate, c.tm)
}

#[cfg(test)]
mod tests {
    #[test]
    fn envelope_keeps_the_flank_and_crests_and_bridges_narrow_valleys() {
        // lane ARCH ed_env: a flank dipping diagonally is unchanged; a narrow valley running down it is bridged to the flank, a
        // ridge kept, a wide valley kept; never below the surface
        let n = 160usize;
        let flank: Vec<f64> = (0..n * n).map(|k| 900.0 - 2.0 * (k % n) as f64 - 3.0 * (k / n) as f64).collect();
        let e = envelope_strike(&flank, n, 8, 16);
        for j in 30..n - 30 { for i in 30..n - 30 { assert!((e[j * n + i] - flank[j * n + i]).abs() < 1e-6, "flank changed"); } }
        // lateral coordinate across the dip (dip direction (2, 3) / |(2, 3)|): c = (3 i - 2 j) / sqrt(13)
        let c = |k: usize| (3.0 * (k % n) as f64 - 2.0 * (k / n) as f64) / 13f64.sqrt();
        let mut f = flank.clone();
        for k in 0..n * n {
            let ck = c(k);
            if (ck - 20.0).abs() < 2.5 { f[k] -= 20.0; }          // valley 5 posts wide
            if (ck - 50.0).abs() < 1.5 { f[k] += 10.0; }          // ridge 3 posts wide
            if (ck + 30.0).abs() < 20.0 { f[k] -= 15.0; }         // valley 40 posts wide
        }
        let e = envelope_strike(&f, n, 8, 16);
        let (mut nv, mut nr, mut nw) = (0, 0, 0);
        for j in 40..n - 40 { for i in 40..n - 40 {
            let k = j * n + i; let ck = c(k);
            assert!(e[k] >= f[k] - 1e-9);
            if (ck - 20.0).abs() < 1.0 { nv += 1; assert!((e[k] - flank[k]).abs() < 1.0, "valley {} vs flank {}", e[k], flank[k]); }
            if (ck - 50.0).abs() < 0.5 { nr += 1; assert!((e[k] - f[k]).abs() < 1e-9, "crest changed"); }
            if (ck + 30.0).abs() < 5.0 { nw += 1; assert!(e[k] - f[k] < 2.0, "wide valley filled {} vs {}", e[k], f[k]); }
        } }
        assert!(nv > 10 && nr > 10 && nw > 10, "{nv} {nr} {nw}");
        // params: absent from the JSON when off; ed_env needs ed_datum
        let mut p = crate::DetailParamsV1::v3();
        p.set("beds", "2").unwrap();
        assert!(p.set("ed_env", "8").is_err());
        p.set("ed_datum", "1").unwrap();
        assert!(!serde_json::to_string(p.beds.as_ref().unwrap()).unwrap().contains("ed_env"));
        assert!(p.set("ed_env_route", "1").is_err());
        p.set("ed_env", "8").unwrap(); p.set("ed_env_route", "1").unwrap();
        let q = p.beds.as_ref().unwrap().ed_env.as_ref().unwrap();
        assert_eq!((q.half_posts, q.route), (8, Some(true)));
    }

    use super::*;
    fn scope() -> Scope { Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap() }
    fn mean_ratio(bb: &[f64], bs: &[f64], lo: f64, hi: f64, cls: u8) -> f64 {
        let (mut top, mut acc, mut h) = (lo, 0.0, 0.0);
        for (b, s) in bb.iter().zip(bs.iter()) { if *b <= lo || *b > hi + 1e-9 { continue; } acc += (b - top) / s; h += b - top; top = *b; }
        (acc / h) * sc_of(cls)
    }
    #[test]
    fn every_style_conserves_the_layer_rock_mass() {
        let p = BedStyleParams::default();
        let s = scope();
        let cases: Vec<(u8, Option<u8>, f64, u8)> = vec![(3, Some(2), 300.0, 0), (3, Some(2), 3000.0, 1), (4, Some(6), 900.0, 1), (4, Some(7), 1500.0, 0), (5, Some(6), 100.0, 0),
            (5, Some(7), 2500.0, 0), (5, Some(4), 6000.0, 1), (6, Some(3), 400.0, 0), (6, Some(5), 400.0, 0), (7, Some(3), 400.0, 0), (8, Some(8), 400.0, 0), (9, Some(9), 400.0, 0)];
        for (i, (lith, env, burial, flags)) in cases.into_iter().enumerate() {
            let rec = LayerRec { lith, env, burial_m: Some(burial), flags, age_base_my: Some(380.0), age_top_my: Some(388.0), thickness_m: 470.0 };
            let mut o = BedsOut::default();
            layer_beds(&rec, 10.0, 480.0, &s, &format!("t/{i}"), &p, &mut o).unwrap();
            let (bb, bs, bk) = (o.bb, o.bs, o.bk);
            assert!((bb.last().unwrap() - 480.0).abs() < 1e-9);
            let r = mean_ratio(&bb, &bs, 10.0, 480.0, lith);
            assert!((r - 1.0).abs() < 1e-9, "lith {lith} env {env:?} ratio {r}");
            // beds are at least rmin thick (the layer has more than one bed) and strictly increasing
            let mut top = 10.0;
            for b in &bb { assert!(*b - top >= p.rmin_m - 1e-9 || bb.len() == 1, "lith {lith}: bed {}", b - top); top = *b; }
            for w in bk.iter().zip(bs.iter()) { assert!((w.0 * w.1 - k_of(lith) * sc_of(lith)).abs() < 1e-9); }
        }
        for kind in ["continental_arc", "hotspot", "large_igneous_province", "unknown"] {
            let v = VolcRec { id: "v".into(), kind: kind.into(), height_m: Some(3000.0), onset_my: Some(376.0), last_eruption_my: Some(400.0) };
            let (o, _) = edifice_beds(&v, 4000.0, &s, &format!("ed/{kind}"), &p).unwrap();
            let (bb, bs) = (o.bb, o.bs);
            let r = mean_ratio(&bb, &bs, 0.0, 4000.0, 3);
            assert!((r - 1.0).abs() < 1e-9, "{kind} ratio {r}");
        }
    }
    #[test]
    fn styles_follow_the_record() {
        let p = BedStyleParams::default();
        let s = scope();
        // spacing variety: lava stack flow units are uneven (CV of hard-bed spacing > 0.5)
        let rec = LayerRec { lith: 3, env: Some(2), burial_m: Some(400.0), flags: 0, age_base_my: Some(390.0), age_top_my: Some(400.0), thickness_m: 2000.0 };
        let mut o = BedsOut::default();
        layer_beds(&rec, 0.0, 2000.0, &s, "lava", &p, &mut o).unwrap();
        let (bb, bs) = (o.bb, o.bs);
        let mut hard = vec![]; let mut top = 0.0;
        for (b, sc) in bb.iter().zip(bs.iter()) { if *sc > sc_of(3) * 1.2 { hard.push(0.5 * (top + b)); } top = *b; }
        let sp: Vec<f64> = hard.windows(2).map(|w| w[1] - w[0]).collect();
        let m = sp.iter().sum::<f64>() / sp.len() as f64;
        let cv = (sp.iter().map(|x| (x - m).powi(2)).sum::<f64>() / sp.len() as f64).sqrt() / m;
        assert!(cv > 0.5 && sp.len() > 50, "lava spacing cv {cv} n {}", sp.len());
        // deeper burial lowers lava contrast
        let spread = |burial: f64| { let r = LayerRec { burial_m: Some(burial), ..rec.clone() }; let mut o = BedsOut::default();
            layer_beds(&r, 0.0, 2000.0, &s, "lava", &p, &mut o).unwrap(); let b = o.bs;
            let l: Vec<f64> = b.iter().map(|x| (x / sc_of(3)).ln()).collect(); let mu = l.iter().sum::<f64>() / l.len() as f64;
            (l.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / l.len() as f64).sqrt() };
        assert!(spread(6000.0) < spread(300.0));
        // a truncated top puts a soft zone at the layer top
        let rt = LayerRec { flags: 1, ..rec.clone() };
        let mut o = BedsOut::default();
        layer_beds(&rt, 0.0, 2000.0, &s, "lava", &p, &mut o).unwrap();
        assert!(o.bs[0] < sc_of(3), "top bed {}", o.bs[0]);
        assert!([FAC_CONTACT_BOLE, FAC_LAVA_TOP, FAC_BOLE].contains(&o.bf[0]), "top facies {}", o.bf[0]);
    }
    #[test]
    fn cycle_sets_conserve_and_vary() {
        let mut p = BedStyleParams::default();
        let s = scope();
        let rec = LayerRec { lith: 4, env: Some(6), burial_m: Some(2400.0), flags: 0, age_base_my: Some(316.0), age_top_my: Some(318.0), thickness_m: 300.0 };
        let sd = |p: &BedStyleParams| { let mut o = BedsOut::default(); layer_beds(&rec, 0.0, 300.0, &s, "ls", p, &mut o).unwrap();
            assert!((mean_ratio(&o.bb, &o.bs, 0.0, 300.0, 4) - 1.0).abs() < 1e-9);
            assert!(o.bf.iter().all(|&f| f == FAC_LS_CLEAN || f == FAC_LS_MARLY), "facies {:?}", o.bf);
            let l: Vec<f64> = o.bs.iter().map(|x| x.ln()).collect(); let mu = l.iter().sum::<f64>() / l.len() as f64;
            (l.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / l.len() as f64).sqrt() };
        let off = sd(&p);
        p.csets = Some(CsetParams::default());
        let on = sd(&p);
        assert!(on > 1.5 * off, "sets raise the strength spread: off {off} on {on}");
    }

    #[test]
    fn merge_conserves_and_removes_thin_beds() {
        let beds = vec![(0.3, 2.0), (3.0, 1.0), (0.2, 0.4), (0.25, 0.5), (5.0, 1.5), (0.1, 3.0)];
        let m = merge_thin(&beds, 1.5);
        let s0: f64 = beds.iter().map(|b| b.0 / b.1).sum();
        let s1: f64 = m.iter().map(|b| b.0 / b.1).sum();
        assert!((s0 - s1).abs() < 1e-12);
        assert!((m.iter().map(|b| b.0).sum::<f64>() - beds.iter().map(|b| b.0).sum::<f64>()).abs() < 1e-12);
        assert!(m.iter().all(|b| b.0 >= 1.5));
    }
    #[test]
    fn up_beds_keep_lava_tops_up() {
        let p = BedStyleParams::default();
        let rec = LayerRec { lith: 3, env: Some(2), burial_m: Some(100.0), flags: 0, age_base_my: None, age_top_my: None, thickness_m: 100.0 };
        let o = up_beds(&rec, 3000.0, &scope(), "up", &p).unwrap();
        let (tops, sc) = (o.bb, o.bs);
        assert!((tops.last().unwrap() - 3000.0).abs() < 1e-9);
        for w in tops.windows(2) { assert!(w[1] > w[0]); }
        assert_eq!(tops.len(), sc.len());
    }
}
