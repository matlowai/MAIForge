//! Lane VOLCANIC-SPLIT (2026-10-21): lava vs volcaniclastic facies inside volcanic strata layers. LABEL: pre-signoff,
//! fidelity tier Preview, derived, never canon. Opt-in (`vs=1`); absent = every byte unchanged.
//!
//! design 2026-10-21 (recorded 21:04): YES to splitting volcanic rock into steep lava/flow rock vs ash/tephra deposits at
//! ~35-40 deg critical slope, and the volcanic rock must read as volcanic. Opt-in keys; defaults change only with evidence.
//!
//! The T0 stratigraphy publishes an arc's subaerial volcanic pile as one lithology-3 layer (on C7 site 1 a 2964 m layer),
//! so the T1b strata cut gave the whole pile basalt's 60 deg critical slope (a threshold landscape at 60 deg: lanes
//! C7-SITE-CHECK, HS-UPLIFT, preview-client). This module splits such a layer, inside the cut, into a deterministic keyed
//! stack of lava-flow packages and volcaniclastic packages (tephra, tuff, volcanic breccia, lahar deposits). The
//! lithology CODE stays 3 (basalt / volcanic) for every package, so every 12-code table downstream is untouched; only the
//! critical slope and the K multiplier of the exposed package change, and the facies is reported separately.
//!
//! The rule is a pure function of (T0 cell, the layer's index in the cell's full column, depth below the strata datum),
//! so every window, every evaluation lattice and pf-detail (which can call `split_column`) see the same beds; package
//! boundaries are planar surfaces parallel to the T0 strata datum (the regional dip the cut already uses).

use serde::{Deserialize, Serialize};

pub const FAC_NONE: u8 = 0;
pub const FAC_LAVA: u8 = 1;
pub const FAC_VC: u8 = 2;

/// Upward / downward extension of the open top and bottom units (m): the keyed stack continues beyond the record so a
/// surface above the datum or below the publish depth still meets packages.
pub const EXT_M: f64 = 6000.0;

pub const LABEL: &str = "EXPLORATORY lane VOLCANIC-SPLIT 2026-10-21 (model specification 2026-10-21 21:04: lava vs ash/tephra split at ~35-40 deg; design review rock rules A-C, design review-rocks/REPORT.md); layer AND windows";

pub const RULE: &str = "volcanic_facies_split_v2: every lithology-3 (basalt / volcanic; in the T0 record always environment 2 subaerial_volcanic) layer of a T0 column in scope (scope 1: cells whose T0 rock_family is 2 arc_magmatic; scope 2: every cell) is cut into alternating lava-flow slots and volcaniclastic packages (tephra fall, tuff, volcanic breccia, lahar deposits), anchored at the layer's top (depth below the strata datum; the top unit's stack continues upward above the datum and the bottom unit's downward below its base, EXT_M each way); package m's thickness is log-normal, t = median exp(sigma z) clamped to [min, max], z = sqrt(-2 ln(1 - u1)) cos(2 pi u2) with u1 = u(cell, layer, m, 1), u2 = u(cell, layer, m, 20) (right-skewed bed thickness, design review rocks report 2.2); the slot of package 0 is lava iff u(cell, layer, 0, 0) < lava_start, then strict alternation; u = blake3(\"vz.t1b.volcanic_split_v1\" | cell u32 | layer u16 | m i64 | purpose u8) first 8 bytes, top 53 bits / 2^53 (world-record keyed: every site, window and tile cuts the same beds). Lithology code stays 3. Critical slope at a post (rule A, rock-mass vs face): the rock-mass angle theta_m (mass_sc) on volcaniclastic rock AND on lava thinner than riser_min at the post (rule B: a thin strong bed is a ledge inside the mass, not a riser, at 32 m); the face angle theta_f (lava_face_sc) only inside a riser lens (lava >= riser_min thick at the post, and theta_f / theta_m >= contrast_min in degrees); K multiplier lava_k on lava, vc_k on volcaniclastic rock";

pub const LENS_RULE: &str = "volcanic_lens_v1 (rule C, lateral pinch-out): a lava slot is lava over its upper share w(x, y) and volcaniclastic below (joining the next package), w = clamp(lens_mid + lens_amp n(x, y), 0, 1), n = mean over q = 0..2 of cos(2 pi (x cos a_q + y sin a_q) / L_q + p_q) at the post's chart position (m), a_q = 2 pi u(.., 2 + 3q), L_q = lambda_lo + (lambda_hi - lambda_lo) u(.., 3 + 3q), p_q = 2 pi u(.., 4 + 3q) keyed per package (a pure function of chart position: seams exact); w = 0 is a pinch-out, w t < riser_min a ledge; so a riser exists only where the lens is locally thick and benches die out along strike (design review rocks report rule C: correlation 1-5 km; recalled, unverified: andesite / basalt flows 0.2-1 km wide, 1-10 km long, lobate)";

pub const CT_RULE: &str = "volcanic_ct_v2 (lane VOLCANIC-SPLIT, ct apply 2 only; replaces the ct class remap on split posts): the period of a split post is [s0, s1] = heights of the tops of the riser lens below and of the riser lens at or above the post (riser: rule B above; at most 96 packages searched each way); s_dn / s_up the next riser tops beyond (published ladder, no regular filling); no riser bounding the post = no staircase (the post keeps its crest-rounded eikonal height: a mantled slope at theta_m). On a period: cliff share fc = eta = the upper riser's lava thickness at the post / (s1 - s0) (rule A: the cliff fraction is the column's hard share, not a function of S); where S <= theta_m (mantled) the profile is 2-limb: soil slope at g = (1 - eta) / (1/S - eta/theta_f) then the cliff (no bench); otherwise limbs climbing up as lane TREAD2: tread at tt (the ct tread tangent; the bench on the lower riser's top), volcaniclastic slope at theta_m (rock-mass angle), cliff at theta_f (face); with S the post's remap gradient (ct_local: the local gradient) the heights satisfy ft + fb + fc = 1 and ft / tt + fb / theta_m + fc / theta_f = 1 / S (relief conserved, F monotone, F(s0) = s0, F(s1) = s1, F <= h); where the 2-limb cliff share cf2 = (1/theta_m - 1/S)/(1/theta_m - 1/theta_f) >= eta (S too steep for a bench) the profile is 2-limb with fc = cf2; S >= theta_f: all cliff; ft > 1 - fc: no volcaniclastic limb (fc from tt and theta_f); results held at or above min(h, base level) (TREAD_BASE_RULE). Regime codes: tread 3, volcaniclastic limb 1, cliff 2. The ct tread drainage rule applies afterwards as for every remapped post";

pub const SOURCE: &str = "design 2026-10-21 21:04 (ash / tephra critical slope ~35-40 deg; lava cliffs and caps ~55-60 deg) and design review rocks report 2026-10-21 (rules A-C; Earth threshold-hillslope modes Fuji 27, Marettimo 29, NZ 32 deg [V: Anand et al. 2023 GRL; Clarke and Burbank 2011 JGR]; mean-hillslope upper bound 35-40 deg [S: DiBiase et al. 2012]; flood-basalt lobe thickness right-skewed, median 8-12 m, >75 m rare [V via summary: Self, Mittal and Jay 2021]; bed thickness log-normal / exponential, CV 0.5-1 [S]; a bed is a T1b riser only if >= ~1.5 posts (48 m) and contrast >= ~1.5). Project assumptions (recalled, unverified, to calibrate): theta_m 37 deg (design's range; design review's Earth band for arc composites 27-33 is run as a sensitivity), theta_f 60 deg, lava K 0.6, volcaniclastic K 1.5 (Sklar and Dietrich 2001 Geology 29:1087: erodibility falls with tensile strength); lava package (1-3 stacked andesite / basalt flows acting as one cliff unit; arc andesite block flows 30-100 m [R], flood-basalt lobes median 8-12 m [V via summary]) median 45 m, sigma 0.55, 5-150 m; volcaniclastic package median 60 m, sigma 0.7, 5-400 m; lens mid 0.75, amp 0.8 (composite piles commonly half or more volcaniclastic: Davidson and De Silva 2000, Encyclopedia of Volcanoes, recalled). No vent-distance term: the T0 layer is a regional pile whose vents are not recorded";

pub const CT_BED_RULE: &str = "volcanic_ct_bed_v1 (lane SPLIT-T1B, opt-in vs_ct_bed=1, ct apply 2 split posts): the cliff limb of a split post's period is the riser lens's lava and nothing else: cliff share fc = eta (the upper riser's lava thickness at the post / (s1 - s0)) at every remap gradient S < theta_f, also where volcanic_ct_v2 took the 2-limb share cf2 = (1/theta_m - 1/S) / (1/theta_m - 1/theta_f) > eta (S too steep for a bench: there the cliff limb grew to cf2 of a period that is often 1-2 km, so neighbouring posts on ONE period published talus -> cliff knots 300-1100 m apart and the 2 m phase blend bent the cliff top once per 32 m cell: lane SPLIT-2M crest scallops); in that range the profile is 2-limb: volcaniclastic limb at g = (1 - eta) / (1/S - eta/theta_f) (theta_m <= g < theta_f; the post's own slope is steeper than the rock-mass angle, the steepness stays in the volcaniclastic limb), then the lava cliff at theta_f: the cliff foot is z2 = s1 - eta (s1 - s0) on every post of the period, its h-knot s1 - eta (s1 - s0) S / theta_f. S >= theta_f: F is the identity (as before) with its knots published on the identity at the lava base (u2 = f2 = 1 - eta), regime cliff on the lava limb, volcaniclastic below; continuous in S at both ends (g -> theta_m where cf2 = eta, g -> theta_f as S -> theta_f). The mantled (S <= theta_m) and 3-limb (cf2 < eta) profiles are unchanged (they already have fc = eta)";

pub const CT_BLEND_RULE: &str = "volcanic_ct_blend_v1 (lane SPLIT-T1B, opt-in vs_ct_blend=1, ct apply 2 split posts): the riser test of the ct ladder is a weight r(t) = smoothstep((t - (riser_min - b)) / 2b) of the lens's lava thickness t at the post (b = band_mm, default 8 m: lava under riser_min - b is never a riser, over riser_min + b always), not the step t >= riser_min (which flipped the period ends of whole posts where a lens crosses 48 m: e.g. s0 -115 vs -1354 m on neighbouring posts, a 19.6 m lattice-aligned step in T1b, lane SPLIT-2M). Upward from the post the candidate period tops are the lenses with r > 0 in order, alternative k with weight r_k prod_{j<k} (1 - r_j) (the walk stops at the first r = 1 lens; the remaining weight, and a missing side, is the no-staircase alternative: the post keeps its height); downward the same for the period base; every (top, base) pair is one volcanic_ct_v2 (or volcanic_ct_bed_v1) remap of the post with weight w_top w_base; the post's height is the weighted mean of the alternatives' heights (continuous in every r, and so in position); the published ladder, knots, regime and cliff share are those of the heaviest alternative with a staircase, flagged held (+128) where the height is not that alternative's F(h) (pf-detail then keeps its inferred corner). With b = 0 the rule is volcanic_ct_v2 exactly";

pub const AA_RULE: &str = "volcanic_aa_v1 (lane SPLIT-T1B, opt-in vs_aa=1): the critical slope of a split post read at height z (the eikonal speed from an upwind height, the regime / remap reads) is the exact climb rate over aa_dx (default 32 m, one window post spacing) through the post's column from z upward: rise R with sum over the pieces crossed of rise_k / slope_k = aa_dx, slope = R / aa_dx; piece slopes theta_m on volcaniclastic rock and thin lava, theta_m + r(t) (theta_f - theta_m) on a lava lens (r the volcanic_ct_blend_v1 riser weight, band b), the class table below / above the split. A post whose upwind height lies just under a riser lens thus climbs part of the step at the face angle instead of all-or-nothing (the binary face / mass switch at one height aliased each lens onto the 32 m lattice). K multiplier, facies and the published classification are unchanged (read at z)";

pub const BAND_MM_DEFAULT: i64 = 8_000;
pub const AA_DX_MM_DEFAULT: i64 = 32_000;

/// Lane SPLIT-T1B: riser weight of a lava lens of thickness t (riser_min infinite = never; band 0 = the step t >= riser_min).
#[inline]
pub fn riser_w(t: f64, riser_min: f64, band: f64) -> f64 {
    if !riser_min.is_finite() { return 0.0; }
    if !(band > 0.0) { return if t >= riser_min { 1.0 } else { 0.0 }; }
    let u = ((t - (riser_min - band)) / (2.0 * band)).clamp(0.0, 1.0);
    u * u * (3.0 - 2.0 * u)
}

/// Lane SPLIT-T1B: one alternative of volcanic_ct_blend_v1: ladder heights [s_dn, s0, s1, s_up] (NaN where none) and eta, or
/// None (no staircase), with its weight.
pub type LadderAlt = (Option<([f64; 4], f64)>, f64);

/// Lane VOLCANIC-SPLIT params (carried inside `LithologyV1.volcanic_split`; None = unchanged).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct VolcanicSplitV1 {
    pub on: u32,
    /// 1: T0 cells with rock_family 2 (arc_magmatic) only; 2: every cell.
    pub scope: u32,
    /// Rule A: rock-mass angle theta_m (tangent permille): volcaniclastic rock and ledgy (thin) lava.
    pub mass_sc_permille: i64,
    /// Rule A: face angle theta_f (tangent permille): lava inside a riser lens.
    pub lava_face_sc_permille: i64,
    pub lava_k_permille: i64,
    pub vc_k_permille: i64,
    /// Rule B: riser lens minimum thickness (mm) and theta_f / theta_m contrast minimum (degree ratio, permille).
    pub riser_min_mm: i64,
    pub contrast_min_permille: i64,
    /// Package thickness: [median mm, sigma permille, min mm, max mm] (log-normal).
    pub lava_t: [i64; 4],
    pub vc_t: [i64; 4],
    /// Probability (permille) that package 0 (at the layer top) is a lava slot.
    pub lava_start_permille: i64,
    /// Rule C lens: w(x, y) = clamp(mid + amp n(x, y), 0, 1) (permille), n from three keyed plane waves of wavelength in
    /// lens_lambda_mm; amp 0 and mid 1000 = planar.
    pub lens_mid_permille: i64,
    pub lens_amp_permille: i64,
    pub lens_lambda_mm: [i64; 2],
    /// 1: volcanic_ct_v2 replaces the ct class remap on split posts; 0: the ct rule reads the split column (code 3).
    pub ct: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ct_rule: Option<String>,
    /// Lane SPLIT-T1B (opt-in `vs_ct_bed=1`, needs ct 1): Some(CT_BED_RULE) = the cliff limb of a split post is the riser
    /// lens's lava (cliff share eta) on every gradient below theta_f, not the 2-limb cf2 share; None = every byte as before.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ct_bed: Option<String>,
    /// Lane SPLIT-T1B (opt-in `vs_ct_blend=1`, needs ct 1): Some(CT_BLEND_RULE) = the riser test of the ct ladder is a
    /// smooth weight across the riser threshold and the post's height is the weighted mean of the ladder alternatives.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ct_blend: Option<String>,
    /// Lane SPLIT-T1B (opt-in `vs_aa=1`): Some(AA_RULE) = the critical slope of a split post is the exact climb rate over
    /// aa_dx through the column (riser weight across the threshold), not the slope at one height.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aa: Option<String>,
    /// Lane SPLIT-T1B: half-width (mm) of the riser-threshold band of ct_blend and aa (None = BAND_MM_DEFAULT).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub band_mm: Option<i64>,
    /// Lane SPLIT-T1B: climb distance (mm) of aa (None = AA_DX_MM_DEFAULT, one 32 m post spacing).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aa_dx_mm: Option<i64>,
    pub rule: String,
    pub lens_rule: String,
    pub source: String,
    pub label: String,
}

impl VolcanicSplitV1 {
    pub fn v1() -> Self {
        VolcanicSplitV1 { on: 1, scope: 1, mass_sc_permille: 754, lava_face_sc_permille: 1732, lava_k_permille: 600, vc_k_permille: 1500,
            riser_min_mm: 48_000, contrast_min_permille: 1_500,
            lava_t: [45_000, 550, 5_000, 150_000], vc_t: [60_000, 700, 5_000, 400_000], lava_start_permille: 500,
            lens_mid_permille: 750, lens_amp_permille: 800, lens_lambda_mm: [1_000_000, 5_000_000], ct: 1, ct_rule: Some(CT_RULE.into()),
            ct_bed: None, ct_blend: None, aa: None, band_mm: None, aa_dx_mm: None,
            rule: RULE.into(), lens_rule: LENS_RULE.into(), source: SOURCE.into(), label: LABEL.into() }
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.on != 1 { return Err("vs: on must be 1 (absent = off)".into()); }
        if !(1..=2).contains(&self.scope) { return Err("vs: scope 1 (arc cells) or 2 (every cell)".into()); }
        if self.mass_sc_permille <= 0 || self.lava_face_sc_permille <= 0 || self.lava_k_permille <= 0 || self.vc_k_permille <= 0 { return Err("vs: slopes and K multipliers must be positive".into()); }
        if self.riser_min_mm < 0 || self.contrast_min_permille < 0 { return Err("vs: riser_min and contrast_min >= 0".into()); }
        for t in [self.lava_t, self.vc_t] { if !(t[0] >= 1_000 && t[1] >= 0 && t[2] >= 1_000 && t[3] >= t[2]) { return Err("vs: thickness [median >= 1 m, sigma >= 0, 1 m <= min <= max]".into()); } }
        if !(0..=1000).contains(&self.lava_start_permille) { return Err("vs: lava_start_permille 0..1000".into()); }
        if !(0..=1000).contains(&self.lens_mid_permille) || self.lens_amp_permille < 0 || !(self.lens_lambda_mm[0] >= 100_000 && self.lens_lambda_mm[1] >= self.lens_lambda_mm[0]) { return Err("vs: lens mid 0..1000, amp >= 0, 100 m <= lambda lo <= hi".into()); }
        if self.ct > 1 || (self.ct == 1) != self.ct_rule.is_some() { return Err("vs: ct 0 or 1 (with its rule)".into()); }
        if (self.ct_bed.is_some() || self.ct_blend.is_some()) && self.ct != 1 { return Err("vs: vs_ct_bed / vs_ct_blend need vs_ct=1".into()); }
        if self.band_mm.map(|b| b < 0).unwrap_or(false) || self.aa_dx_mm.map(|d| d <= 0).unwrap_or(false) { return Err("vs: vs_band_mm >= 0, vs_aa_dx_mm > 0".into()); }
        Ok(())
    }
    /// Lane SPLIT-T1B: riser-threshold band half-width (m).
    #[inline] pub fn band(&self) -> f64 { self.band_mm.unwrap_or(BAND_MM_DEFAULT) as f64 / 1000.0 }
    #[inline] pub fn aa_dx(&self) -> f64 { self.aa_dx_mm.unwrap_or(AA_DX_MM_DEFAULT) as f64 / 1000.0 }
    /// Lane SPLIT-T1B: true when any lane SPLIT-T1B key is set (they change physics, not draws: `draw_twin` strips them).
    pub fn t1b_keys(&self) -> bool { self.ct_bed.is_some() || self.ct_blend.is_some() || self.aa.is_some() || self.band_mm.is_some() || self.aa_dx_mm.is_some() }
    /// Lane SPLIT-T1B: the params without the lane SPLIT-T1B keys (the draw digest's view; every vs=1 output keeps its draws).
    pub fn draw_twin(&self) -> Self { VolcanicSplitV1 { ct_bed: None, ct_blend: None, aa: None, band_mm: None, aa_dx_mm: None, ..self.clone() } }
    /// Lane SPLIT-T1B: the ct-only keys removed (a ct_window=1 layer is solved without ct, so they cannot act there).
    pub fn without_ct_keys(&self) -> Self { VolcanicSplitV1 { ct_bed: None, ct_blend: None, band_mm: if self.aa.is_some() { self.band_mm } else { None }, ..self.clone() } }
    #[inline] pub fn mass(&self) -> f64 { self.mass_sc_permille as f64 / 1000.0 }
    #[inline] pub fn face(&self) -> f64 { self.lava_face_sc_permille as f64 / 1000.0 }
    /// Rule B contrast test (degree ratio of theta_f to theta_m).
    #[inline] pub fn contrast_ok(&self) -> bool { self.face().atan() / self.mass().atan() * 1000.0 >= self.contrast_min_permille as f64 }
    #[inline] pub fn riser_min(&self) -> f64 { if self.contrast_ok() { self.riser_min_mm as f64 / 1000.0 } else { f64::INFINITY } }
    /// Whether a T0 cell with this rock_family is in scope.
    #[inline] pub fn cell_in_scope(&self, rock_family: u8) -> bool { self.scope == 2 || rock_family == 2 }
}

/// Keyed uniform [0, 1).
pub fn u01(cell: u32, layer: u16, m: i64, purpose: u8) -> f64 {
    let mut h = blake3::Hasher::new();
    h.update(b"vz.t1b.volcanic_split_v1");
    h.update(&cell.to_le_bytes());
    h.update(&layer.to_le_bytes());
    h.update(&m.to_le_bytes());
    h.update(&[purpose]);
    let b = h.finalize();
    let x = u64::from_le_bytes(b.as_bytes()[0..8].try_into().unwrap());
    (x >> 11) as f64 / (1u64 << 53) as f64
}

/// A column after the split: units top-down with depth tops and bases below the datum (the top unit open upward, the
/// bottom one downward, as the strata cut reads them), codes, parent layer index in the full column, facies (lava SLOT
/// for lava packages: lava over its lens share, volcaniclastic below), lens waves (kx, ky, phase) x 3 per lava slot.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SplitCol { pub top: Vec<f64>, pub base: Vec<f64>, pub code: Vec<u8>, pub idx: Vec<u16>, pub fac: Vec<u8>, pub lens: Vec<[f64; 9]>, pub mid: f64, pub amp: f64 }

impl SplitCol {
    /// Slot index at depth d below the datum (the cut's convention: on a base, the lower unit).
    #[inline]
    pub fn unit_at(&self, d: f64) -> usize { self.base.partition_point(|b| *b <= d).min(self.code.len().saturating_sub(1)) }
    /// Lava thickness of slot j at chart (x, y) (0 for a non-lava slot).
    #[inline]
    pub fn lava_t(&self, j: usize, x: f64, y: f64) -> f64 {
        if self.fac[j] != FAC_LAVA { return 0.0; }
        let t = self.base[j] - self.top[j];
        if self.amp == 0.0 { return self.mid.clamp(0.0, 1.0) * t; }
        let l = &self.lens[j];
        let n = ((l[0] * x + l[1] * y + l[2]).cos() + (l[3] * x + l[4] * y + l[5]).cos() + (l[6] * x + l[7] * y + l[8]).cos()) / 3.0;
        (self.mid + self.amp * n).clamp(0.0, 1.0) * t
    }
    /// (slot, facies, riser) at depth d, chart (x, y); riser: lava of a lens at least `riser_min` thick here.
    #[inline]
    pub fn facies_at(&self, d: f64, x: f64, y: f64, riser_min: f64) -> (usize, u8, bool) {
        let j = self.unit_at(d);
        let f = self.fac[j];
        if f != FAC_LAVA { return (j, f, false); }
        let lt = self.lava_t(j, x, y);
        if d - self.top[j] < lt { (j, FAC_LAVA, lt >= riser_min) } else { (j, FAC_VC, false) }
    }
    /// Riser ladder at depth d: top depths of the riser lenses (r_up beyond, r_top at or above d, r_bot below d, r_dn
    /// beyond) and the lava thickness of the r_top lens; None if no riser at or above / below within 96 slots or a
    /// non-volcanic unit intervenes before one is found.
    pub fn ladder(&self, d: f64, x: f64, y: f64, riser_min: f64) -> Option<([f64; 4], f64)> {
        let j = self.unit_at(d);
        if self.fac[j] == FAC_NONE || !riser_min.is_finite() { return None; }
        let riser = |q: usize| -> Option<f64> { if self.fac[q] == FAC_LAVA { let t = self.lava_t(q, x, y); if t >= riser_min { return Some(t); } } None };
        let (mut ups, mut dns): (Vec<(f64, f64)>, Vec<f64>) = (Vec::new(), Vec::new());
        let mut q = j as i64;
        while q >= 0 && (j as i64 - q) <= 96 && ups.len() < 2 {
            let k = q as usize;
            if self.fac[k] == FAC_NONE { break; }
            if self.top[k] <= d { if let Some(t) = riser(k) { ups.push((self.top[k], t)); } }
            q -= 1;
        }
        for k in (j + 1)..self.code.len().min(j + 97) {
            if dns.len() >= 2 || self.fac[k] == FAC_NONE { break; }
            if riser(k).is_some() { dns.push(self.top[k]); }
        }
        let (&(rt, lava), &rb) = (ups.first()?, dns.first()?);
        Some(([ups.get(1).map(|v| v.0).unwrap_or(f64::NAN), rt, rb, dns.get(1).copied().unwrap_or(f64::NAN)], lava))
    }
    /// Lane SPLIT-T1B volcanic_ct_blend_v1 at depth d: the weighted ladder alternatives in depths ([r_up, r_top, r_bot, r_dn],
    /// lava of r_top; None = no staircase); None when `ladder` would be None for a non-ladder reason (non-volcanic unit at d,
    /// no riser possible). Weights sum to 1.
    pub fn ladder_alts(&self, d: f64, x: f64, y: f64, riser_min: f64, band: f64) -> Option<Vec<(Option<([f64; 4], f64)>, f64)>> {
        let j = self.unit_at(d);
        if self.fac[j] == FAC_NONE || !riser_min.is_finite() { return None; }
        let cand = |q: usize| -> Option<(f64, f64)> { if self.fac[q] == FAC_LAVA { let t = self.lava_t(q, x, y); let r = riser_w(t, riser_min, band); if r > 0.0 { return Some((t, r)); } } None };
        // chains (top depth, lava, weight) up to the first r = 1 candidate, plus one candidate beyond (weight 0) for s_up / s_dn
        let chain = |up: bool| -> (Vec<(f64, f64, f64)>, f64) {
            let mut out: Vec<(f64, f64, f64)> = Vec::new();
            let mut rem = 1.0f64;
            let mut q = j as i64;
            let step = if up { -1i64 } else { 1 };
            if !up { q += 1; }
            let mut n = 0;
            while q >= 0 && (q as usize) < self.code.len() && n <= 96 {
                let k = q as usize;
                if self.fac[k] == FAC_NONE { break; }
                if !up || self.top[k] <= d {
                    if let Some((t, r)) = cand(k) {
                        if rem <= 0.0 { out.push((self.top[k], t, 0.0)); break; }
                        out.push((self.top[k], t, rem * r)); rem *= 1.0 - r;
                        if rem < 1e-12 { rem = 0.0; }
                    }
                }
                q += step; n += 1;
            }
            (out, rem)
        };
        let (ups, rem_u) = chain(true);
        let (dns, rem_d) = chain(false);
        let mut alts: Vec<(Option<([f64; 4], f64)>, f64)> = Vec::new();
        let mut none = 0.0;
        for (a, u) in ups.iter().enumerate() {
            if u.2 <= 0.0 { continue; }
            for (b, v) in dns.iter().enumerate() {
                if v.2 <= 0.0 { continue; }
                let l = [ups.get(a + 1).map(|w| w.0).unwrap_or(f64::NAN), u.0, v.0, dns.get(b + 1).map(|w| w.0).unwrap_or(f64::NAN)];
                alts.push((Some((l, u.1)), u.2 * v.2));
            }
            none += u.2 * rem_d;
        }
        none += rem_u;
        if none > 0.0 { alts.push((None, none)); }
        Some(alts)
    }
    /// Lane SPLIT-T1B volcanic_aa_v1: the piece of the column just above depth d (going up): (slope, depth of its top).
    pub fn piece_up(&self, d: f64, x: f64, y: f64, riser_min: f64, band: f64, mass: f64, face: f64, none_sc: &dyn Fn(u8) -> f64) -> (f64, f64) {
        let dq = d - 1e-6;
        let j = self.unit_at(dq);
        match self.fac[j] {
            FAC_NONE => (none_sc(self.code[j]), self.top[j]),
            FAC_LAVA => {
                let t = self.lava_t(j, x, y);
                if dq - self.top[j] < t { (mass + riser_w(t, riser_min, band) * (face - mass), self.top[j]) } else { (mass, self.top[j] + t) }
            }
            _ => (mass, self.top[j]),
        }
    }
}

/// Lane SPLIT-T1B volcanic_aa_v1: the exact climb rate over dx from depth d upward; `piece(d)` = (slope, top depth) of the
/// piece above d.
pub fn climb_rate(d: f64, dx: f64, piece: &dyn Fn(f64) -> (f64, f64)) -> f64 {
    let (mut dd, mut rem, mut rise) = (d, dx, 0.0f64);
    let mut last = 0.0;
    for _ in 0..256 {
        let (sp, pt) = piece(dd);
        let sp = sp.max(1e-9);
        last = sp;
        let r = dd - pt;
        if !(r.is_finite()) || r / sp >= rem { rise += sp * rem; rem = 0.0; break; }
        let r = r.max(0.0);
        rise += r; rem -= r / sp; dd = pt.min(dd - 1e-6);
    }
    if rem > 0.0 { rise += last * rem; }
    rise / dx
}

/// Log-normal package thickness (m).
#[inline]
fn lognormal(t: [i64; 4], cell: u32, layer: u16, m: i64) -> f64 {
    let (u1, u2) = (u01(cell, layer, m, 1), u01(cell, layer, m, 20));
    let z = (-2.0 * (1.0 - u1).ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
    (t[0] as f64 / 1000.0 * (t[1] as f64 / 1000.0 * z).exp()).clamp(t[2] as f64 / 1000.0, t[3] as f64 / 1000.0)
}

/// Facies and thickness of package m of layer `layer` in `cell`.
#[inline]
fn package(v: &VolcanicSplitV1, cell: u32, layer: u16, m: i64, start_lava: bool) -> (u8, f64) {
    let lava = (m.rem_euclid(2) == 0) == start_lava;
    (if lava { FAC_LAVA } else { FAC_VC }, lognormal(if lava { v.lava_t } else { v.vc_t }, cell, layer, m))
}

fn lens_of(v: &VolcanicSplitV1, cell: u32, layer: u16, m: i64) -> [f64; 9] {
    let mut l = [0.0; 9];
    if v.lens_amp_permille == 0 { return l; }
    let (lo, hi) = (v.lens_lambda_mm[0] as f64 / 1000.0, v.lens_lambda_mm[1] as f64 / 1000.0);
    for q in 0..3u8 {
        let a = std::f64::consts::TAU * u01(cell, layer, m, 2 + 3 * q);
        let lam = lo + (hi - lo) * u01(cell, layer, m, 3 + 3 * q);
        let k = std::f64::consts::TAU / lam;
        let q = q as usize;
        l[3 * q] = k * a.cos(); l[3 * q + 1] = k * a.sin(); l[3 * q + 2] = std::f64::consts::TAU * u01(cell, layer, m, 4 + 3 * q as u8);
    }
    l
}

/// Split a column (the cut's kept units: depth bases, codes, full-column indices) of a T0 cell; `in_scope` from the cell's
/// rock_family. Units that are not lithology 3 are copied with facies 0. None when nothing was split.
pub fn split_column(v: &VolcanicSplitV1, cell: u32, in_scope: bool, base: &[f64], code: &[u8], idx: &[u16]) -> Option<SplitCol> {
    let n = code.len();
    if !in_scope || n == 0 || !code.iter().any(|&c| c == crate::t1b_litho::LITHO_BASALT) { return None; }
    let mut o = SplitCol { mid: v.lens_mid_permille as f64 / 1000.0, amp: v.lens_amp_permille as f64 / 1000.0, ..Default::default() };
    let push = |o: &mut SplitCol, top: f64, b: f64, c: u8, ix: u16, f: u8, l: [f64; 9]| { o.top.push(top); o.base.push(b); o.code.push(c); o.idx.push(ix); o.fac.push(f); o.lens.push(l); };
    for j in 0..n {
        let top_j = if j == 0 { f64::NEG_INFINITY } else { base[j - 1] };
        if code[j] != crate::t1b_litho::LITHO_BASALT { push(&mut o, top_j, base[j], code[j], idx[j], FAC_NONE, [0.0; 9]); continue; }
        let anchor = if j == 0 { 0.0 } else { base[j - 1] };
        let lim_bot = if j + 1 == n { base[j].max(anchor) + EXT_M } else { base[j] };
        let start_lava = u01(cell, idx[j], 0, 0) * 1000.0 < v.lava_start_permille as f64;
        let lens = |f: u8, m: i64| if f == FAC_LAVA { lens_of(v, cell, idx[j], m) } else { [0.0; 9] };
        if j == 0 {
            let mut up: Vec<(f64, f64, u8, i64)> = Vec::new();
            let (mut d, mut m) = (anchor, -1i64);
            while d > -EXT_M {
                let (f, t) = package(v, cell, idx[j], m, start_lava);
                up.push((d - t, d, f, m));
                d -= t; m -= 1;
            }
            for &(t, b, f, m) in up.iter().rev() { push(&mut o, t, b, code[j], idx[j], f, lens(f, m)); }
        }
        let (mut d, mut m) = (anchor, 0i64);
        while d < lim_bot {
            let (f, t) = package(v, cell, idx[j], m, start_lava);
            let b = (d + t).min(lim_bot);
            push(&mut o, d, b, code[j], idx[j], f, lens(f, m));
            d = b; m += 1;
        }
    }
    Some(o)
}

/// Per-post access for a lattice (the strata cut's RockPosts): the split columns of the cells the posts use.
#[derive(Clone, Debug)]
pub struct VolcPosts {
    pub rule: VolcanicSplitV1,
    pub cols: Vec<SplitCol>,
    /// Per post: index into `cols`, u32::MAX = the cell is not split (the original cut applies).
    pub slot: Vec<u32>,
    /// Per post chart position (m), for the lava lenses.
    pub xy: Vec<(f64, f64)>,
}

/// One post's volcanic_ct_v2 remap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolcRemap { pub z: f64, pub rg: u8, pub fc: f64, pub ft: f64, pub fb: f64, pub knots: [f64; 4], pub clamped: bool }

/// volcanic_ct_v2 profile of one post (see CT_RULE): period [lo, hi] (heights), gradient s, cliff share eta, tangents
/// tt (tread), tb (theta_m), ta (theta_f), base level bl (NaN none). Knots as t1b_ct::knots2 / knots3 (fractions).
pub fn ct_remap(s: f64, eta: f64, tt: f64, tb: f64, ta: f64, h: f64, lo: f64, hi: f64, bl: f64) -> VolcRemap { ct_remap_bed(s, eta, tt, tb, ta, h, lo, hi, bl, false) }

/// volcanic_ct_v2, or with `bed` (lane SPLIT-T1B vs_ct_bed=1) volcanic_ct_bed_v1 (CT_BED_RULE): the cliff limb stays the lens's
/// lava share eta where S is too steep for a bench (2-limb volcaniclastic slope at g, theta_m <= g < theta_f) and the identity
/// for S >= theta_f publishes its knots at the lava base.
pub fn ct_remap_bed(s: f64, eta: f64, tt: f64, tb: f64, ta: f64, h: f64, lo: f64, hi: f64, bl: f64, bed: bool) -> VolcRemap {
    let tl = hi - lo;
    let id = VolcRemap { z: h, rg: crate::t1b_ct::RG_CLIFF, fc: 1.0, ft: 0.0, fb: 0.0, knots: [0.0, 0.0, 0.0, 0.0], clamped: false };
    if bed && tl > 0.0 && s >= ta && eta > 0.0 && eta < 1.0 {
        let u = ((h - lo) / tl).clamp(0.0, 1.0);
        let rg = if u >= 1.0 - eta { crate::t1b_ct::RG_CLIFF } else { crate::t1b_ct::RG_TALUS };
        return VolcRemap { z: h, rg, fc: eta, ft: 0.0, fb: 1.0 - eta, knots: [0.0, 0.0, 1.0 - eta, 1.0 - eta], clamped: false };
    }
    if !(tl > 0.0) || !(s > 0.0) || s >= ta { return id; }
    let cf2 = if s > tb { ((1.0 / tb - 1.0 / s) / (1.0 / tb - 1.0 / ta)).clamp(0.0, 1.0) } else { 0.0 };
    // Mantled (S <= theta_m): the riser keeps its lens share as cliff, the rest of the period is a soil slope gentler than S
    // (no bench: a bench is a stripped cap top, which needs a slope at the threshold).
    let (tb, ft, fb, fc) = if s <= tb { let g = (1.0 - eta) / (1.0 / s - eta / ta); if g > 0.0 && eta < 1.0 { (g, 0.0, 1.0 - eta, eta) } else { return id } }
    else if bed && cf2 >= eta && eta > 0.0 && eta < 1.0 && 1.0 / s - eta / ta > 0.0 { ((1.0 - eta) / (1.0 / s - eta / ta), 0.0, 1.0 - eta, eta) }
    else if cf2 >= eta { (tb, 0.0, 1.0 - cf2, cf2) } else {
        let fc = eta;
        let ft = (1.0 / s - fc / ta - (1.0 - fc) / tb) / (1.0 / tt - 1.0 / tb);
        if ft > 1.0 - fc { let fc1 = ((1.0 / tt - 1.0 / s) / (1.0 / tt - 1.0 / ta)).clamp(0.0, 1.0); (tb, 1.0 - fc1, 0.0, fc1) } else { (tb, ft.max(0.0), 1.0 - fc - ft.max(0.0), fc) }
    };
    let u = ((h - lo) / tl).clamp(0.0, 1.0);
    let ut = ft * s / tt;
    let ub = ut + fb * s / tb;
    let (v, rg) = if u < ut { (u * tt / s, crate::t1b_ct::RG_TREAD) } else if u < ub { (ft + (u - ut) * tb / s, crate::t1b_ct::RG_TALUS) } else { ((ft + fb + (u - ub) * ta / s).min(1.0), crate::t1b_ct::RG_CLIFF) };
    let mut z = lo + tl * v;
    let mut clamped = false;
    if bl.is_finite() { let f = h.min(bl); if z < f { z = f; clamped = true; } }
    VolcRemap { z, rg, fc, ft, fb, knots: [ut.min(1.0), ft, ub.min(1.0), ft + fb], clamped }
}

impl VolcPosts {
    /// `cell[i]` / `xy[i]` the post's T0 cell and chart position (m); `col(c)` the cut's kept column of cell c (bases,
    /// codes, full indices; None if empty).
    pub fn build(rule: &VolcanicSplitV1, cell: &[u32], xy: Vec<(f64, f64)>, rock_family: &[u8], col: &dyn Fn(u32) -> Option<(Vec<f64>, Vec<u8>, Vec<u16>)>) -> VolcPosts {
        let mut seen: std::collections::BTreeMap<u32, u32> = std::collections::BTreeMap::new();
        let mut cols = Vec::new();
        let slot = cell.iter().map(|&c| *seen.entry(c).or_insert_with(|| {
            let s = rule.cell_in_scope(rock_family[c as usize]);
            match col(c).and_then(|(b, k, i)| split_column(rule, c, s, &b, &k, &i)) { Some(sc) => { cols.push(sc); (cols.len() - 1) as u32 } None => u32::MAX }
        })).collect();
        VolcPosts { rule: rule.clone(), cols, slot, xy }
    }
    #[inline]
    pub fn col(&self, i: usize) -> Option<&SplitCol> { let s = self.slot[i]; if s == u32::MAX { None } else { Some(&self.cols[s as usize]) } }
    /// (critical slope, K multiplier, code, facies) at depth d below the datum, None when post i's cell is not split.
    #[inline]
    pub fn at(&self, i: usize, d: f64, table_sc: &dyn Fn(u8) -> (f64, f64)) -> Option<(f64, f64, u8, u8)> {
        let c = self.col(i)?;
        let (j, f, riser) = c.facies_at(d, self.xy[i].0, self.xy[i].1, self.rule.riser_min());
        // Lane SPLIT-T1B volcanic_aa_v1: the exact climb rate over aa_dx (K, code and facies still read at d).
        if self.rule.aa.is_some() {
            let v = &self.rule;
            let (x, y) = self.xy[i];
            let sc = climb_rate(d, v.aa_dx(), &|dd| c.piece_up(dd, x, y, v.riser_min(), v.band(), v.mass(), v.face(), &|cd| table_sc(cd).0));
            let k = match f { FAC_NONE => table_sc(c.code[j]).1, FAC_LAVA => v.lava_k_permille as f64 / 1000.0, _ => v.vc_k_permille as f64 / 1000.0 };
            return Some((sc, k, c.code[j], f));
        }
        match f {
            FAC_NONE => { let (s, k) = table_sc(c.code[j]); Some((s, k, c.code[j], f)) }
            FAC_LAVA => Some((if riser { self.rule.face() } else { self.rule.mass() }, self.rule.lava_k_permille as f64 / 1000.0, c.code[j], f)),
            _ => Some((self.rule.mass(), self.rule.vc_k_permille as f64 / 1000.0, c.code[j], f)),
        }
    }
    #[inline]
    pub fn facies(&self, i: usize, d: f64) -> u8 { self.col(i).map(|c| c.facies_at(d, self.xy[i].0, self.xy[i].1, self.rule.riser_min()).1).unwrap_or(FAC_NONE) }
    /// volcanic_ct_v2 at post i (datum height `datum`, surface z): None = not a split post or ct coupling off (the ct class
    /// rule applies); Some(None) = split post with no riser period (no staircase); Some(Some((ladder heights [s_dn, s0, s1,
    /// s_up] (NaN where none), eta))).
    pub fn ct_ladder(&self, i: usize, datum: f64, z: f64) -> Option<Option<([f64; 4], f64)>> {
        if self.rule.ct != 1 { return None; }
        let c = self.col(i)?;
        let d = datum - z;
        if c.fac[c.unit_at(d)] == FAC_NONE { return None; }
        Some(c.ladder(d, self.xy[i].0, self.xy[i].1, self.rule.riser_min()).and_then(|(l, lava)| {
            let (s1, s0) = (datum - l[1], datum - l[2]);
            if !(s1 > s0) { return None; }
            Some(([datum - l[3], s0, s1, datum - l[0]], (lava / (s1 - s0)).clamp(0.0, 1.0)))
        }))
    }
    /// Lane SPLIT-T1B volcanic_ct_blend_v1 at post i: None = not a split post or ct off (as `ct_ladder`); Some(alternatives)
    /// with ladder heights [s_dn, s0, s1, s_up] and eta (None = no staircase), weights summing to 1.
    pub fn ct_ladder_alts(&self, i: usize, datum: f64, z: f64) -> Option<Vec<LadderAlt>> {
        if self.rule.ct != 1 { return None; }
        let c = self.col(i)?;
        let d = datum - z;
        if c.fac[c.unit_at(d)] == FAC_NONE { return None; }
        Some(match c.ladder_alts(d, self.xy[i].0, self.xy[i].1, self.rule.riser_min(), self.rule.band()) {
            None => vec![(None, 1.0)],
            Some(a) => heights_of(&a, datum),
        })
    }
}

/// Lane SPLIT-T1B: ladder alternatives in depths -> heights (as `ct_ladder`: an empty period is the no-staircase alternative).
pub fn heights_of(a: &[(Option<([f64; 4], f64)>, f64)], datum: f64) -> Vec<LadderAlt> {
    a.iter().map(|(l, w)| (l.and_then(|(l, lava)| {
        let (s1, s0) = (datum - l[1], datum - l[2]);
        if !(s1 > s0) { return None; }
        Some(([datum - l[3], s0, s1, datum - l[0]], (lava / (s1 - s0)).clamp(0.0, 1.0)))
    }), *w)).collect()
}

/// Measurement only: env `VZ_VOLC_DUMP=<dir>` (never digested, no output bytes): the window publish writes the facies of
/// every post (u8: 0 none / not split / over-coded, 1 lava, 2 volcaniclastic, 255 water) as `<dir>/<rx>_<ry>_facies.u8`;
/// an existing file is left untouched.
pub fn dump_dir() -> Option<std::ffi::OsString> { std::env::var_os("VZ_VOLC_DUMP") }

pub fn dump_facies(dir: &std::ffi::OsStr, rx: i64, ry: i64, fac: &[u8]) {
    use std::io::Write;
    let p = std::path::Path::new(dir).join(format!("{rx}_{ry}_facies.u8"));
    let r = std::fs::create_dir_all(dir).and_then(|_| std::fs::OpenOptions::new().write(true).create_new(true).open(&p)).and_then(|mut f| f.write_all(fac));
    match r { Ok(()) => println!("  VZ_VOLC_DUMP: wrote {} (measurement only)", p.display()), Err(e) => println!("  VZ_VOLC_DUMP: {} not written: {e}", p.display()) }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn split_is_keyed_lognormal_lensed_and_ladders_follow_risers() {
        let mut v = VolcanicSplitV1::v1();
        assert!(v.validate().is_ok() && v.contrast_ok());
        let (b, k, ix) = (vec![36.0, 3000.0, 8000.0], vec![3u8, 3, 1], vec![0u16, 1, 2]);
        let s = split_column(&v, 123471, true, &b, &k, &ix).unwrap();
        assert_eq!(s, split_column(&v, 123471, true, &b, &k, &ix).unwrap(), "pure function");
        assert!(split_column(&v, 123471, false, &b, &k, &ix).is_none() && split_column(&v, 1, true, &[100.0], &[5], &[0]).is_none());
        assert!(s.base.windows(2).all(|w| w[1] > w[0]));
        assert!((1..s.base.len()).all(|j| s.top[j] == s.base[j - 1]));
        assert!(s.base.contains(&36.0) && s.base.contains(&3000.0));
        assert_eq!((*s.code.last().unwrap(), *s.fac.last().unwrap()), (1, FAC_NONE));
        assert!(s.code.iter().zip(&s.fac).all(|(c, f)| (*c == 3) == (*f != FAC_NONE)));
        assert!(s.base[0] > -EXT_M - 1e-9 && s.idx[0] == 0);
        // Alternation, thickness bounds and a right-skewed spread (CV > 0.4) inside the 2964 m layer.
        let l1: Vec<usize> = (0..s.code.len()).filter(|&j| s.idx[j] == 1).collect();
        for w in l1.windows(2) { assert_ne!(s.fac[w[0]], s.fac[w[1]]); }
        let th: Vec<f64> = l1[..l1.len() - 1].iter().map(|&j| s.base[j] - s.top[j]).collect();
        assert!(th.iter().all(|t| *t >= 5.0 - 1e-9 && *t <= 400.0 + 1e-9));
        let m = th.iter().sum::<f64>() / th.len() as f64;
        let sd = (th.iter().map(|t| (t - m).powi(2)).sum::<f64>() / th.len() as f64).sqrt();
        assert!(sd / m > 0.4, "cv {}", sd / m);
        // Lens: pure function of position; pinch-outs exist along a 20 km line; risers rarer than lava.
        let j = l1.iter().copied().find(|&j| s.fac[j] == FAC_LAVA).unwrap();
        assert_eq!(s.lava_t(j, 1234.5, -987.0), s.lava_t(j, 1234.5, -987.0));
        let (mut lava, mut riser, mut zero, mut nt) = (0, 0, 0, 0);
        for xi in 0..200 { for dd in 0..300 {
            let (q, f, r) = s.facies_at(40.0 + dd as f64 * 9.7, xi as f64 * 100.0, 3000.0, v.riser_min());
            nt += 1; if f == FAC_LAVA { lava += 1; } if r { riser += 1; }
            if s.fac[q] == FAC_LAVA && s.lava_t(q, xi as f64 * 100.0, 3000.0) == 0.0 { zero += 1; }
        } }
        assert!(lava > 0 && riser > 0 && riser < lava && zero > 0 && lava < nt / 2, "lava {lava} riser {riser} zero {zero} of {nt}");
        // Ladder with planar lenses: riser tops bound the period, the volcaniclastic interval belongs to the riser above.
        v.lens_amp_permille = 0; v.lens_mid_permille = 1000; v.riser_min_mm = 30_000;
        let p = split_column(&v, 123471, true, &b, &k, &ix).unwrap();
        let rs: Vec<usize> = (0..p.code.len()).filter(|&q| p.idx[q] == 1 && p.fac[q] == FAC_LAVA && p.base[q] - p.top[q] >= 30.0).collect();
        assert!(rs.len() >= 4);
        let (a, bb) = (rs[1], rs[2]);
        let (l, lv) = p.ladder(p.top[a] + 0.5, 0.0, 0.0, v.riser_min()).unwrap();
        assert_eq!(l[1], p.top[a]); assert_eq!(l[2], p.top[bb]); assert_eq!(l[0], p.top[rs[0]]); assert_eq!(l[3], p.top[rs[3]]);
        assert!((lv - (p.base[a] - p.top[a])).abs() < 1e-9);
        assert_eq!(p.ladder(p.top[bb] - 0.5, 0.0, 0.0, v.riser_min()).unwrap().0[1], p.top[a]);
        // No riser when the contrast test fails.
        let mut w = v.clone(); w.lava_face_sc_permille = 900;
        assert!(!w.contrast_ok() && p.ladder(p.top[a] + 0.5, 0.0, 0.0, w.riser_min()).is_none());
        // ct_remap: monotone, ends fixed, relief conserved, F <= h; cliff share = eta when S allows a bench.
        let (s_, eta, tt, tb, ta) = (0.85, 0.3, 0.14, 0.754, 1.732);
        let mut last = f64::NEG_INFINITY;
        for q in 0..=200 {
            let h = 100.0 + q as f64;
            let r = ct_remap(s_, eta, tt, tb, ta, h, 100.0, 300.0, f64::NAN);
            assert!(r.z >= last - 1e-9 && r.z <= h + 1e-9); last = r.z;
            assert!((r.fc - eta).abs() < 1e-12 && r.ft > 0.0 && (r.ft / tt + r.fb / tb + r.fc / ta - 1.0 / s_).abs() < 1e-9);
        }
        assert!((ct_remap(s_, eta, tt, tb, ta, 300.0, 100.0, 300.0, f64::NAN).z - 300.0).abs() < 1e-9 && ct_remap(s_, eta, tt, tb, ta, 100.0, 100.0, 300.0, f64::NAN).z == 100.0);
        let r = ct_remap(1.5, 0.3, tt, tb, ta, 200.0, 100.0, 300.0, f64::NAN);
        assert!(r.ft == 0.0 && r.fc > 0.3, "too steep for a bench: 2-limb with the needed cliff share");
        assert_eq!(ct_remap(1.5, 0.3, tt, tb, ta, 150.0, 100.0, 300.0, 140.0).z.max(140.0), ct_remap(1.5, 0.3, tt, tb, ta, 150.0, 100.0, 300.0, 140.0).z);
        // Per-post access.
        let v = VolcanicSplitV1::v1();
        let vp = VolcPosts::build(&v, &[123471, 7], vec![(0.0, 0.0), (0.0, 0.0)], &[2u8; 200000], &|c| if c == 7 { None } else { Some((b.clone(), k.clone(), ix.clone())) });
        assert!(vp.col(0).is_some() && vp.col(1).is_none());
        let (sc, km, code, f) = vp.at(0, 1500.0, &|_| (9.0, 9.0)).unwrap();
        assert_eq!(code, 3);
        assert!((f == FAC_LAVA && (sc == 1.732 || sc == 0.754) && km == 0.6) || (f == FAC_VC && sc == 0.754 && km == 1.5));
        assert_eq!(vp.at(0, 7000.0, &|c| (c as f64, 0.5)).unwrap(), (1.0, 0.5, 1, FAC_NONE));
        assert!(vp.ct_ladder(1, 1000.0, 0.0).is_none() && vp.ct_ladder(0, 1000.0, -6000.0).is_none());
        if let Some(Some((l, e))) = vp.ct_ladder(0, 1000.0, 0.0) { assert!(l[1] < l[2] && e > 0.0 && e <= 1.0); }
        assert_eq!(serde_json::to_value(&v).unwrap()["mass_sc_permille"], 754);
    }

    /// Lane SPLIT-T1B: absent keys are not serialized; the draw twin is the key-off params.
    #[test]
    fn t1b_keys_opt_in_and_stripped_from_draws() {
        let v = VolcanicSplitV1::v1();
        let js = serde_json::to_string(&v).unwrap();
        assert!(!js.contains("ct_bed") && !js.contains("ct_blend") && !js.contains("\"aa\"") && !js.contains("band_mm") && !v.t1b_keys());
        let mut w = v.clone(); w.ct_bed = Some(CT_BED_RULE.into()); w.aa = Some(AA_RULE.into()); w.band_mm = Some(4000);
        assert!(w.t1b_keys() && w.validate().is_ok() && w.draw_twin() == v);
        assert_eq!(w.without_ct_keys().ct_bed, None); assert_eq!(w.without_ct_keys().band_mm, Some(4000));
        let mut x = v.clone(); x.ct_blend = Some(CT_BLEND_RULE.into()); x.band_mm = Some(4000);
        assert_eq!(x.without_ct_keys(), v);
        let mut y = v.clone(); y.ct = 0; y.ct_rule = None; y.ct_bed = Some(CT_BED_RULE.into());
        assert!(y.validate().is_err());
    }

    /// Lane SPLIT-T1B volcanic_ct_bed_v1: the cliff foot z2 = s1 - eta T on every gradient; continuous with volcanic_ct_v2 at the
    /// branch boundary (cf2 = eta) and at S -> theta_f; identical to volcanic_ct_v2 where that already had fc = eta.
    #[test]
    fn bed_tied_cliff_limb() {
        let (eta, tt, tb, ta) = (0.0333, 0.14, 0.754, 1.732);
        let (lo, hi) = (-1861.0, 328.0);
        let tl = hi - lo;
        for s in [0.5, 0.66, 0.76, 0.9, 1.1, 1.25, 1.6, 1.731, 1.732, 2.5] {
            let mut last = f64::NEG_INFINITY;
            for q in 0..=400 {
                let h = lo + tl * q as f64 / 400.0;
                let r = ct_remap_bed(s, eta, tt, tb, ta, h, lo, hi, f64::NAN, true);
                assert!(r.z >= last - 1e-9 && r.z <= h + 1e-9, "monotone, F <= h (s {s})"); last = r.z;
                assert!((r.fc - eta).abs() < 1e-12, "fc = eta (s {s})");
                assert!((r.knots[3] - (1.0 - eta)).abs() < 1e-9, "cliff foot at the lava base (s {s})");
            }
            let a = ct_remap_bed(s, eta, tt, tb, ta, hi, lo, hi, f64::NAN, true);
            assert!((a.z - hi).abs() < 1e-9);
            // unchanged where v2 already had fc = eta
            let v2 = ct_remap(s, eta, tt, tb, ta, 0.0, lo, hi, f64::NAN);
            if (v2.fc - eta).abs() < 1e-12 { assert_eq!(v2, ct_remap_bed(s, eta, tt, tb, ta, 0.0, lo, hi, f64::NAN, true)); }
        }
        // continuity at the cf2 = eta boundary and at theta_f
        let sb = 1.0 / (eta / ta + (1.0 - eta) / tb);
        for h in [-1000.0, 0.0, 250.0, 300.0] {
            let a = ct_remap_bed(sb * (1.0 - 1e-9), eta, tt, tb, ta, h, lo, hi, f64::NAN, true);
            let b = ct_remap_bed(sb * (1.0 + 1e-9), eta, tt, tb, ta, h, lo, hi, f64::NAN, true);
            assert!((a.z - b.z).abs() < 1e-3 && (a.knots[2] - b.knots[2]).abs() < 1e-6, "h {h}: {} {}", a.z, b.z);
            let c = ct_remap_bed(ta * (1.0 - 1e-9), eta, tt, tb, ta, h, lo, hi, f64::NAN, true);
            let e = ct_remap_bed(ta, eta, tt, tb, ta, h, lo, hi, f64::NAN, true);
            assert!((c.z - e.z).abs() < 1e-3 && (c.knots[2] - e.knots[2]).abs() < 1e-6);
        }
        // the measured SPLIT-2M case (post (2784, 4096): S 1.249, published cliff share 0.702): the knot is now 73 m below s1
        let r = ct_remap_bed(1.249, eta, tt, tb, ta, 297.0, lo, hi, f64::NAN, true);
        assert!((lo + tl * r.knots[3] - (hi - eta * tl)).abs() < 1e-6);
        let old = ct_remap(1.249, eta, tt, tb, ta, 297.0, lo, hi, f64::NAN);
        assert!(old.fc > 0.69 && lo + tl * old.knots[3] < -1000.0);
        assert!((r.z - old.z).abs() < 1e-9, "a post on the cliff limb keeps its height");
    }

    /// Lane SPLIT-T1B volcanic_ct_blend_v1 / volcanic_aa_v1: band 0 reproduces the step ladder; weights sum to 1 and move
    /// continuously across the threshold; the climb rate is continuous in depth and equals the slope inside one piece.
    #[test]
    fn blend_and_climb() {
        let v = VolcanicSplitV1::v1();
        let (b, k, ix) = (vec![36.0, 3000.0, 8000.0], vec![3u8, 3, 1], vec![0u16, 1, 2]);
        let s = split_column(&v, 123471, true, &b, &k, &ix).unwrap();
        let rm = v.riser_min();
        for xi in 0..60 { for dd in 0..60 {
            let (x, d) = (xi as f64 * 97.0, 50.0 + dd as f64 * 41.3);
            let a0 = s.ladder_alts(d, x, 3000.0, rm, 0.0).unwrap();
            let want = s.ladder(d, x, 3000.0, rm);
            assert_eq!(a0.len(), 1);
            match (want, a0[0].0) { (Some(w), Some(g)) => assert!(w.1 == g.1 && w.0.iter().zip(g.0.iter()).all(|(a, b)| a.to_bits() == b.to_bits()), "{w:?} vs {g:?}"), (None, None) => {}, (w, g) => panic!("{w:?} vs {g:?}") }
            let a = s.ladder_alts(d, x, 3000.0, rm, 8.0).unwrap();
            assert!((a.iter().map(|q| q.1).sum::<f64>() - 1.0).abs() < 1e-9);
        } }
        assert_eq!(riser_w(40.0, 48.0, 8.0), 0.0); assert_eq!(riser_w(56.0, 48.0, 8.0), 1.0); assert!((riser_w(48.0, 48.0, 8.0) - 0.5).abs() < 1e-12);
        assert_eq!(riser_w(48.0, f64::INFINITY, 8.0), 0.0);
        // climb: inside one piece the rate is the piece slope; continuous across a lens top
        let pc = |d: f64| s.piece_up(d, 0.0, 3000.0, rm, 8.0, v.mass(), v.face(), &|_| 0.5);
        let mut last = f64::NAN;
        for q in 0..4000 {
            let d = 100.0 + q as f64 * 0.25;
            let c = climb_rate(d, 32.0, &pc);
            assert!(c >= v.mass() - 1e-9 && c <= v.face() + 1e-9, "d {d}: {c}");
            if last.is_finite() { assert!((c - last).abs() < 0.05, "jump at d {d}: {last} -> {c}"); }
            last = c;
        }
    }
}
