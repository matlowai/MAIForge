//! T1b 0.9.1 COAST-CUT and VOLCANO CONSTRUCT (lane CO, 2026-10-13), optional additions to lane RV's 0.9.0
//! river-network-first construction. LABEL: pre-signoff, fidelity tier Preview, derived, never canon; a construction
//! (no SimTime). Absent (`RiverNetParamsV1::coast` / `volcano` = None) nothing here runs and every 0.9.0 recipe
//! serializes and evaluates exactly as before. Constants pre-registered in
//! `unpublished-design-record`; project assumptions until ruled.
//!
//! Coast (after design review's Q4 prescription, `unpublished-design-record`):
//! 1. E is tapered by the plain T0 bilinear, so range fronts follow its contours and the coast gets a lowland;
//! 2. the shelf down to a declared lowstand depth joins the construction as land, with the construction water at a
//!    base level `sea - F`; F is one flooding offset per coastal T0 cell, solved with E so each cell's land fraction
//!    meets F4; valleys graded to the lowered base level drown when the surface is flooded at present sea level;
//! 3. mouth deposition: a fan at each bayhead sized by the sediment supply E A over the highstand time;
//! 4. wave cut by rock class: hard classes keep a cliff, soft classes get a beach.
//!
//! Volcano: an edifice per T0 volcano record point (lane MG `volcanoes.json`), shaped by kind and eroded by the
//! RECORD age only (design 2026-10-13 via integrator: "geology-derived generation"; no presentation age override).

use crate::t0::T0;
use crate::t1b::{CHANNEL, FREE, GUIDE, LAKE, MARINE};
use crate::t1b_net::hat4;
use crate::t1b_rivernet::{Curve, EvalLat};
use crate::*;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;
use vz_keyed_rng::Scope;
use vz_tilekey::chart::{Chart, RAD2DEG};

pub const CV_GENERATOR_VERSION: &str = "0.9.1";
pub const PURPOSE_GULLY: &str = "t1b.volcano.gully";
pub const PURPOSE_DIKE: &str = "t1b.volcano.dike";
pub const PURPOSE_LOBE: &str = "t1b.volcano.lobe";
/// Exposed-lithology codes written over deposits and edifices.
pub const LITHO_DEPOSIT: u8 = 10;
pub const LITHO_VOLCANIC: u8 = 3;
/// `Surface::over` codes.
pub const OVER_NONE: u8 = 0;
pub const OVER_EDIFICE: u8 = 1;
pub const OVER_FAN: u8 = 2;
pub const OVER_BEACH: u8 = 3;
pub const OVER_CUT: u8 = 4;

// ------------------------------------------------------------------ parameters

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CoastParamsV1 {
    /// Construction water: T0 bilinear <= sea - shelf_depth (the lowstand shoreline on the T0 bathymetry).
    pub shelf_depth_mm: i64,
    pub flood_init_mm: i64,
    pub flood_max_mm: i64,
    /// Offset update gain per unit land-fraction excess, and its step clamp.
    pub flood_gain_mm: i64,
    pub flood_step_max_mm: i64,
    pub flood_tolerance_permille: i64,
    /// A T0 cell is coastal when its bilinear land fraction lies in (min, 1 - min).
    pub coastal_min_permille: i64,
    /// E taper: tau = t_min + (1 - t_min) smoothstep(clamp(bil / H, 0, 1)); taper_height 0 = off.
    pub taper_height_mm: i64,
    pub taper_min_permille: i64,
    pub fan_slope_ppm: i64,
    pub fan_time_yr: i64,
    pub fan_depth_min_mm: i64,
    pub fan_radius_max_mm: i64,
    pub fan_area_min_m2: i64,
    /// Hard (cliff) iff k_mult <= this.
    pub cliff_k_max_permille: i64,
    /// Cliff retreat = cliff_retreat x k_mult.
    pub cliff_retreat_mm: i64,
    pub platform_depth_mm: i64,
    /// Beach width = beach_width x k_mult, at most beach_width_max.
    pub beach_width_mm: i64,
    pub beach_width_max_mm: i64,
    pub beach_slope_ppm: i64,
    pub taper_rule: String,
    pub shelf_rule: String,
    pub flood_rule: String,
    pub fan_rule: String,
    pub wave_rule: String,
    pub conduit_rule: String,
}

impl CoastParamsV1 {
    pub fn v1() -> Self {
        CoastParamsV1 {
            shelf_depth_mm: 120_000, flood_init_mm: 60_000, flood_max_mm: 250_000, flood_gain_mm: 1_500_000, flood_step_max_mm: 40_000,
            flood_tolerance_permille: 30, coastal_min_permille: 10, taper_height_mm: 400_000, taper_min_permille: 250,
            fan_slope_ppm: 2_000, fan_time_yr: 7_000, fan_depth_min_mm: 5_000, fan_radius_max_mm: 4_000_000, fan_area_min_m2: 1_000_000,
            cliff_k_max_permille: 800, cliff_retreat_mm: 250_000, platform_depth_mm: 1_000, beach_width_mm: 150_000, beach_width_max_mm: 450_000, beach_slope_ppm: 20_000,
            taper_rule: "bilinear_taper_v1: E(x) x tau(bil), tau = t_min + (1 - t_min) s(clamp(bil / H, 0, 1)), s(u) = u^2 (3 - 2u), on channel beds and hillslope gradients (design review Q4 item 1: the range front follows the plain T0 bilinear contours; a coastal lowland); the E solve sees the tapered field".into(),
            shelf_rule: "lowstand_shelf_v1: construction water = posts and 512 m nodes with T0 bilinear <= sea - shelf_depth (declared lowstand shoreline on the T0 bathymetry; glacial eustatic lowstand ~120 m, recalled Fairbanks 1989 Nature 342:637); the shelf between joins the network and the eikonal as land; construction water stands at base level sea - F(x), F = plain bilinear hat blend of one flooding offset per coastal T0 cell (0 elsewhere)".into(),
            flood_rule: "f4_flood_offset_v1: after construction every post with z <= sea is sea (published at sea level); per coastal T0 cell (bilinear land fraction of its 128 m chart posts in (min, 1 - min)) the offset F_c is updated with E in the same Jacobi pass: F_c <- clamp(F_c + gain x (land fraction - bilinear land fraction), 0, F_max), step clamped; land fraction = plain fraction of the cell's chart posts that are land after flooding, fans, wave cut and edifices".into(),
            fan_rule: "B5 (exploratory): the mouth is the first of a bed crossing sea level and an entry into construction water (bil <= sea - shelf); never in open construction water; at an entry with the bed above sea the fan is a fan delta of apex = that bed, z = sea + (apex - sea)(1 - d / R)^1.5 inside R. bayhead_fan_v1: bayhead = first point along a channel curve where the bed passes from above to at or below sea level (chart-wide, from the curves; duplicates within 1 m dropped); V = E(bayhead, untapered hat blend) x A x fan_time (steady-state sediment flux over the highstand; recalled ~7 ky, Stanley and Warne 1994 Science 265:228); h = max(F(bayhead), fan_depth_min); R = min(sqrt(2 V / (pi h)), fan_radius_max) (a half disc of thickness h); z = max(z, sea + fan_slope x (R - d)) for d < R on every post except T1a channels; only bayheads with A >= fan_area_min; deposited posts expose class 10".into(),
            wave_rule: "class_wave_cut_v1: d = 8-neighbour chamfer distance (1, sqrt 2 spacings) from sea posts; on land posts that are not T1a channels: k_mult <= cliff_k_max (hard) and d < cliff_retreat x k_mult -> platform at sea - platform_depth (sea; the coast retreats and leaves a cliff as high as the ground there); otherwise (soft) and d < min(beach_width x k_mult, beach_width_max) -> z = min(z, sea + beach_slope x d), class 10; rock at the post's own height".into(),
            conduit_rule: "t1a_mouth_conduit_v1: a T1a seed whose T1a receiver is neither a seed nor water (a mouth on the exposed shelf or inland of the zero contour) is continued by a conduit: steepest descent of the geodesic distance to construction water over land nodes (8-neighbour Dijkstra, Euclidean steps; ties by node index), joining an earlier conduit where met; conduits are routing seeds (others drain into them), channels, unjittered, carry the T1a area plus their own, and their beds rise from the base level by the stream-power slope but never above the terminal T1a level nor below their receiver".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.shelf_depth_mm < 0 || self.flood_max_mm < 0 || self.flood_init_mm < 0 || self.flood_init_mm > self.flood_max_mm { return Err("coast: shelf / flood bounds".into()); }
        if self.taper_min_permille < 0 || self.taper_min_permille > 1000 || self.coastal_min_permille < 0 || self.coastal_min_permille >= 500 { return Err("coast: taper / coastal bounds".into()); }
        Ok(())
    }
    #[inline]
    pub fn tau(&self, bil: f64) -> f64 {
        if self.taper_height_mm <= 0 { return 1.0; }
        let u = (bil / (self.taper_height_mm as f64 / 1000.0)).clamp(0.0, 1.0);
        let t = self.taper_min_permille as f64 / 1000.0;
        t + (1.0 - t) * u * u * (3.0 - 2.0 * u)
    }
    pub fn shelf_m(&self) -> f64 { self.shelf_depth_mm as f64 / 1000.0 }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct EdificeShapeV1 {
    pub name: String,
    /// Record kinds using this shape.
    pub kinds: Vec<String>,
    pub s_top_ppm: i64,
    pub s_base_ppm: i64,
    pub crater_radius_mm: i64,
    pub crater_depth_mm: i64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct VolcanoParamsV1 {
    pub clear_t1a_mm: i64,
    pub mask_inner_mm: i64,
    /// REVISION 1: the T1a corridor caps edifice relief by cap_slope x (d - mask_inner) (was a smoothstep to mask_outer).
    pub cap_slope_ppm: i64,
    /// REVISION 1: edifice relief fades to 0 between R and fade_ratio x R.
    pub fade_ratio_permille: i64,
    pub erosion_mult_permille: i64,
    pub gully_count: u32,
    pub gully_width_mm: i64,
    pub gully_width_ppm_of_r: i64,
    pub gully_depth_permille: i64,
    pub neck_radius_mm: i64,
    pub neck_relief_permille: i64,
    pub wall_slope_ppm: i64,
    pub dike_count: u32,
    pub dike_length_permille_of_r: i64,
    pub dike_width_mm: i64,
    pub dike_relief_permille: i64,
    pub shapes: Vec<EdificeShapeV1>,
    pub placement_rule: String,
    pub profile_rule: String,
    pub age_rule: String,
    pub mask_rule: String,
}

impl VolcanoParamsV1 {
    pub fn v1() -> Self {
        let s = |name: &str, kinds: &[&str], st: i64, sb: i64, rc: i64, dc: i64| EdificeShapeV1 { name: name.into(), kinds: kinds.iter().map(|k| k.to_string()).collect(), s_top_ppm: st, s_base_ppm: sb, crater_radius_mm: rc, crater_depth_mm: dc };
        VolcanoParamsV1 {
            clear_t1a_mm: 3_000_000, mask_inner_mm: 150_000, cap_slope_ppm: 450_000, fade_ratio_permille: 1300, erosion_mult_permille: 1000,
            gully_count: 16, gully_width_mm: 60_000, gully_width_ppm_of_r: 60_000, gully_depth_permille: 600,
            neck_radius_mm: 350_000, neck_relief_permille: 150, wall_slope_ppm: 1_732_000,
            dike_count: 4, dike_length_permille_of_r: 350, dike_width_mm: 50_000, dike_relief_permille: 300,
            shapes: vec![
                s("stratovolcano", &["island_arc", "continental_arc", "rift"], 650_000, 90_000, 350_000, 150_000),
                s("shield", &["hotspot", "large_igneous_province"], 140_000, 50_000, 1_500_000, 150_000),
            ],
            placement_rule: "vent_placement_v1 (REVISION 1 clearance): the record point (its T0 cell centre: the record carries no subgrid position); if it lies closer than c = max(clear_t1a, r_c + H / cap_slope) to a T1a channel curve (code 1) (so the corridor cap never reaches the summit), it moves to c from the nearest T1a curve point along the away direction (on the curve: the normal toward the higher T0 bilinear); up to 4 passes; the move and the vent's T0 cell are reported".into(),
            profile_rule: "concave_edifice_v1 (REVISION 1 foot): H = record subgrid summit - T0 bilinear at the vent (skipped if <= 0); below the crater rim r_c the slope is s(rho) = s_top exp(-rho / lambda), lambda = H / (s_top - s_base), to R = r_c + lambda ln(s_top / s_base) (drop H) and on beyond R with the same law; the edifice relief above the ground fades to 0 between R and fade_ratio x R (smoothstep); inside r_c a bowl of depth d_c; surface = ground + max(0, edifice - ground) x fade (a construct). Shapes by record kind: stratovolcano (arcs, rifts; s_top 0.65 ~ 33 deg, s_base 0.09 ~ 5 deg, crater 350 m / 150 m), shield (hotspots, LIPs; 0.14 / 0.05, caldera 1500 m / 150 m) (recalled typical values, e.g. Fuji/Rainier flanks, Hawaiian shields)".into(),
            age_rule: "record_age_degradation_v1: t = record snapshot time - last eruption (the record's only age; no override); delta = erosion_mult x E(vent, untapered hat blend of the solved table) x t; delta < H: summit and flanks lowered (height above the base scaled by (H - delta) / H), gully_count radial gullies at keyed angles (t1b.volcano.gully: uniform jitter within each of gully_count equal sectors) of depth min(delta, gully_depth x height above base) x radial fade, V cross-section of half-width gully_width + gully_width_ppm x r; crater breached along the T0 bilinear downhill direction once delta >= crater depth; delta >= H: residual = neck (radius neck_radius, relief neck_relief x H above the ground, walls at wall_slope) + dike_count dikes at keyed angles (t1b.volcano.dike), 0.35 R long, dike_width wide, dike_relief x neck relief, tapering to 0; residual relief is ADDED to the ground".into(),
            mask_rule: "t1a_corridor_cap_v1 (REVISION 1): edifice relief above the ground <= cap_slope x max(0, d - mask_inner), d = distance to the nearest T1a channel curve (chart-wide polylines, so every window sees the same d), a valley side of atan(cap_slope) (F1: T1a channels keep their beds and a valley); network channels raised more than 1 m are buried (free)".into(),
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.cap_slope_ppm <= 0 || self.fade_ratio_permille <= 1000 || self.shapes.is_empty() { return Err("volcano: cap / fade / shapes".into()); }
        for s in &self.shapes { if s.s_top_ppm <= s.s_base_ppm || s.s_base_ppm <= 0 { return Err("volcano: s_top must exceed s_base > 0".into()); } }
        Ok(())
    }
}

// ------------------------------------------------------------------ lane CF (0.9.2) parameters

pub const CF_GENERATOR_VERSION: &str = "0.9.2";
pub const PURPOSE_FLOW: &str = "t1b.volcano.flow";
pub const PURPOSE_YGULLY: &str = "t1b.volcano.ygully";

/// T1b 0.9.2 (lane CF, 2026-10-14): cone scaling and foot blend, lava-flow lobes and radial gullies, vent clearance of the
/// foot, fans kept attached (no wave cut at mouths; T1a graded across the shelf), damped E/F solve, cliff profile.
/// Every item has its own switch (0 = off). Pre-registered in `unpublished-design-record`.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CfParamsV1 {
    /// C1: H / W_B (permille, 0 = CO's profile); W_B clamped to [min, max].
    pub cone_hw_permille: i64,
    pub cone_wb_min_mm: i64,
    pub cone_wb_max_mm: i64,
    /// C2: burial weight 1 inside inner x R_f, smoothstep to 0 at outer x R_f; softmax width (0 = hard max); foot taper.
    pub bury_inner_permille: i64,
    pub bury_outer_permille: i64,
    pub softmax_k_mm: i64,
    pub foot_taper_permille: i64,
    /// C3: vent clearance from T1a channels >= clear_foot x R_f (0 = off).
    pub clear_foot_permille: i64,
    /// C4: lava-flow lobes (0 = off).
    pub flow_count: u32,
    pub flow_r0_permille: [i64; 2],
    pub flow_len_permille: [i64; 2],
    pub flow_halfwidth_mm: [i64; 2],
    pub flow_thickness_mm: [i64; 2],
    /// C5: radial gullies (0 = off).
    pub gully_count: u32,
    pub gully_r0_permille: [i64; 2],
    pub gully_depth_mm: [i64; 2],
    pub gully_halfwidth_mm: i64,
    pub gully_halfwidth_ppm_of_r: i64,
    /// C6: no wave cut within a fan's radius + pad of its apex; T1a node levels graded across the shelf (flags 0/1).
    pub fan_nocut: u32,
    pub fan_nocut_pad_mm: i64,
    pub t1a_shelf_grade: u32,
    /// C7: per-cell sign-flip halving of the E step and the F gain (flag); multiplier floor.
    pub solve_damp: u32,
    pub solve_damp_floor_permille: i64,
    /// C8: cliff profile (flag): platform, talus at talus_slope over talus_width, wall at wall_slope.
    pub cliff_profile: u32,
    /// REVISION 1 (R1a): the cone sits on the falling regional surface (flag); (R1b) axis wander amplitude / wavelength ranges;
    /// (R1c) a fan whose apex T0 depth exceeds this fraction of the shelf depth is subaqueous.
    pub sloped_base: u32,
    pub wander_amp_mm: [i64; 2],
    pub wander_len_mm: [i64; 2],
    pub fan_subaqueous_permille: i64,
    pub talus_width_mm: i64,
    pub talus_slope_ppm: i64,
    pub wall_slope_ppm: i64,
    pub cone_rule: String,
    pub foot_rule: String,
    pub clearance_rule: String,
    pub flow_rule: String,
    pub gully_rule: String,
    pub fan_rule: String,
    pub solve_rule: String,
    pub cliff_rule: String,
    /// Lane VH-I (opt-in `ring_plain=1`, window-only): the composite ring-plain tail. None = off (params JSON unchanged).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ring_plain: Option<RingPlainV1>,
}

/// Lane VH-I ring_plain_tail_v1 (interim of `unpublished-design-record` section 2.8 I1).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct RingPlainV1 {
    pub version: u32,
    /// The tail starts where the CF profile slope falls to s_start and its slope decays linearly to s_end over decay x R_f.
    pub s_start_ppm: i64,
    pub s_end_ppm: i64,
    pub decay_permille_of_rf: i64,
    /// Beyond the decay the tail keeps s_end and its thickness cap falls linearly to 0 over taper x R_f (the outer edge).
    pub taper_permille_of_rf: i64,
    /// Thickness cap over the pre-edifice ground (T_max).
    pub t_max_mm: i64,
    /// Ramp-in of the tail thickness over this width (of R_f) ending at the start radius.
    pub ramp_permille_of_rf: i64,
    pub rule: String,
}

impl RingPlainV1 {
    pub fn v1() -> Self {
        RingPlainV1 { version: 1, s_start_ppm: 90_000, s_end_ppm: 20_000, decay_permille_of_rf: 500, taper_permille_of_rf: 250, t_max_mm: 60_000, ramp_permille_of_rf: 200,
            rule: "ring_plain_tail_v1 (lane VH-I, interim of the eruption-history design I1; numbers are project assumptions [H]): composite cones with a CF geometry (lane VE broad families excluded). The smooth CF profile z_p(rho) = top - s_top lambda (1 - exp(-rho / lambda)) x scale (rho = r - r_c, the record-age scale as edifice_at; no flows, gullies, roughness or crater) reaches the slope s_start (0.09, 5 deg) at rho_5 = lambda ln(s_top scale / s_start) (0 when the profile starts gentler; no tail for the straight-cone case). The tail surface T(r) = z_p(min(rho, rho_5)) - D(max(0, rho - rho_5)) + the CF sloped-base term at the post, D the drop of a slope that falls linearly from s_start to s_end (0.02) over L = decay x R_f (0.5) and stays s_end beyond. Raise over the pre-edifice ground g: dt = clamp(T - g, 0, c(r)), c = T_max (60 m) up to rho_5 + L, falling linearly to 0 at rho_5 + L + taper x R_f (0.25), ramped in by smoothstep over ramp x R_f (0.2) ending at r_5; the post takes max(cone raise, dt), so surface = max(ground, min(tail, ground + T_max)) on the ring plain: it fills valleys and moats up to T_max and lets ridges stand out. Not on T1a channel, lake, marine or below-sea posts; the raise obeys CO's T1a corridor cap and B3 shadow like the cone; the tail reach is capped at the T1a corridor segment reach (influence + H / cap slope + 1 km, less 300 m). Window-only (the layer E / F solve is unchanged). Recalled, unverified: stratovolcano ring plains (Taranaki, Ruapehu) slope 1-3 deg and are built by lahars and debris flows filling the surrounding lows".into() }
    }
    pub fn validate(&self) -> R<()> {
        if self.s_start_ppm <= self.s_end_ppm || self.s_end_ppm <= 0 || self.decay_permille_of_rf <= 0 || self.taper_permille_of_rf < 0 || self.t_max_mm < 0 || self.ramp_permille_of_rf < 0 { return Err("ring_plain: bounds".into()); }
        Ok(())
    }
}

/// Lane VH-I: (rho_5, L, taper) in m for a cone, None when the tail does not apply (no CF geometry, broad VE family, straight
/// cone).
fn ring_plain_geom(rp: &RingPlainV1, v: &VolcanoSite) -> Option<(f64, f64, f64)> {
    let g = v.cf.as_ref()?;
    if let Some(ve) = v.ve.as_ref() { if ve.family != "composite" { return None; } }
    if v.denudation_m >= v.h_m || v.h_m <= 0.0 || v.lambda_m >= 1e8 { return None; }
    let scale = (v.h_m - v.denudation_m) / v.h_m;
    let s0 = rp.s_start_ppm as f64 / 1e6;
    let st = v.s_top * scale;
    let rho5 = if st > s0 { v.lambda_m * (st / s0).ln() } else { 0.0 };
    let rf = g.r_foot_m;
    Some((rho5, rp.decay_permille_of_rf as f64 / 1000.0 * rf, rp.taper_permille_of_rf as f64 / 1000.0 * rf))
}

/// Lane VH-I: the outer reach (m from the vent) of a cone's ring-plain tail (0 when none), capped at the T1a corridor reach of
/// `t1a_segments_near` less 300 m (so the corridor cap always sees the T1a channels near a raised post).
pub fn ring_plain_reach(rp: &RingPlainV1, vp: &VolcanoParamsV1, v: &VolcanoSite) -> f64 {
    let Some((rho5, l, tp)) = ring_plain_geom(rp, v) else { return 0.0 };
    let seg = v.r_influence_m.max(v.ve.as_ref().map(|g| g.r_units_m).unwrap_or(0.0)) + v.h_m / (vp.cap_slope_ppm as f64 / 1e6) + 1000.0 - 300.0;
    (v.crater_radius_m + rho5 + l + tp).min(seg)
}

/// Lane VH-I ring_plain_tail_v1 at chart (x, y): (tail surface T, thickness cap c) or None outside the tail.
pub fn ring_plain_at(rp: &RingPlainV1, vp: &VolcanoParamsV1, v: &VolcanoSite, cf: &CfParamsV1, x: f64, y: f64, bil_xy: f64) -> Option<(f64, f64)> {
    let (rho5, l, tp) = ring_plain_geom(rp, v)?;
    let rf = v.cf.as_ref()?.r_foot_m;
    let r = (x - v.vent_xy_m[0]).hypot(y - v.vent_xy_m[1]);
    if r >= ring_plain_reach(rp, vp, v) { return None; }
    let rc = v.crater_radius_m;
    let r5 = rc + rho5;
    let ramp = rp.ramp_permille_of_rf as f64 / 1000.0 * rf;
    let w_in = if ramp <= 0.0 { if r >= r5 { 1.0 } else { 0.0 } } else { smoothstep01((r - (r5 - ramp)) / ramp) };
    if w_in <= 0.0 { return None; }
    let scale = (v.h_m - v.denudation_m) / v.h_m;
    let top = v.summit_built() - v.denudation_m;
    let rho = (r - rc).max(0.0);
    let zp = |q: f64| top - v.s_top * v.lambda_m * (1.0 - (-q / v.lambda_m).exp()) * scale;
    let (s0, s1) = (rp.s_start_ppm as f64 / 1e6, rp.s_end_ppm as f64 / 1e6);
    let t = (rho - rho5).max(0.0);
    let drop = if t <= l { s0 * t - (s0 - s1) * t * t / (2.0 * l) } else { s0 * l - (s0 - s1) * l / 2.0 + s1 * (t - l) };
    let sb = if cf.sloped_base == 1 { let d = bil_xy - v.base_m; let k = 100.0; if d >= 0.0 { 0.0 } else if d > -2.0 * k { -d * d / (4.0 * k) } else { d + k } } else { 0.0 };
    let tail = zp(rho.min(rho5)) - drop + sb;
    let tmax = rp.t_max_mm as f64 / 1000.0;
    let cap = if t <= l { tmax } else if tp > 0.0 { tmax * (1.0 - (t - l) / tp).max(0.0) } else { 0.0 };
    Some((tail, cap * w_in))
}

/// Lane VH-I: the ring-plain raise of a post over its pre-edifice ground g (0 off / outside / not land above sea).
#[inline]
fn ring_plain_raise(cf: &CfParamsV1, vp: &VolcanoParamsV1, v: &VolcanoSite, x: f64, y: f64, bxy: f64, g: f64, sea: f64, fixed: u8) -> f64 {
    let Some(rp) = cf.ring_plain.as_ref() else { return 0.0 };
    if fixed == MARINE || !(g > sea) { return 0.0; }
    match ring_plain_at(rp, vp, v, cf, x, y, bxy) { Some((t, c)) => (t - g).min(c).max(0.0), None => 0.0 }
}

impl CfParamsV1 {
    pub fn v1() -> Self {
        CfParamsV1 {
            cone_hw_permille: 150, cone_wb_min_mm: 5_000_000, cone_wb_max_mm: 25_000_000,
            bury_inner_permille: 550, bury_outer_permille: 950, softmax_k_mm: 25_000, foot_taper_permille: 150,
            clear_foot_permille: 1000,
            flow_count: 7, flow_r0_permille: [150, 450], flow_len_permille: [450, 850], flow_halfwidth_mm: [150_000, 400_000], flow_thickness_mm: [15_000, 40_000],
            gully_count: 20, gully_r0_permille: [120, 300], gully_depth_mm: [15_000, 45_000], gully_halfwidth_mm: 40_000, gully_halfwidth_ppm_of_r: 20_000,
            fan_nocut: 0, fan_nocut_pad_mm: 64_000, t1a_shelf_grade: 1,
            sloped_base: 1, wander_amp_mm: [80_000, 250_000], wander_len_mm: [1_500_000, 3_000_000], fan_subaqueous_permille: 500,
            solve_damp: 1, solve_damp_floor_permille: 31,
            cliff_profile: 1, talus_width_mm: 64_000, talus_slope_ppm: 700_000, wall_slope_ppm: 2_750_000,
            cone_rule: "cone_scaling_v1: foot (basal) diameter W_B = H / q clamped to [W_min, W_max], R_f = W_B / 2 (recalled, unverified: arc stratovolcano H/W_B ~0.1-0.2, Pike and Clow 1981 USGS OFR 81-1038, Grosse et al. 2014 Geology; Fuji ~0.10, Mayon ~0.12, Rainier ~0.17, Taranaki ~0.2); concave profile s = s_top exp(-rho / lambda) with lambda solved by bisection so the drop over R_f - r_c is H (if s_top (R_f - r_c) <= H the profile is a straight cone at s_top ending at r_c + H / s_top); s_base = slope at R_f (reported)".into(),
            foot_rule: "REVISION 1 (R1a sloped_base): E += min(0, bil(x) - bil(vent)) (the cone sits on the falling T0 surface; its lower flanks run under the sea where T0 is sea). cone_foot_blend_v1: the cone surface E (profile + flows - gullies, crater) continues beyond R_f; burial g' = g - w max(0, g - E), w = 1 inside bury_inner R_f, smoothstep to 0 at bury_outer R_f (the cone is younger than the T1b relief beneath it and buries it; the record cell's subgrid relief is the volcano); S = softmax_k(E, g') = max + k ln(1 + exp(-|E - g'| / k)); surface = g + (S - g) u, u = 1 inside R_f, smoothstep to 0 at R_end(theta) = R_f (1 + taper (1 + 0.8 L(theta))) (CO's keyed lobe sum); the surface may be lower than g inside the footprint (buried ridges); raised relief obeys CO's T1a corridor cap and B3 shadow".into(),
            clearance_rule: "REVISION 1 (R1e): R_f is re-evaluated from the H at the current vent position on every placement pass. vent_clearance_foot_v1: CO's clearance max(clear_t1a, r_c + H / s_cap) is extended by clear_foot x R_f (R_f from the record point's H), so the main flank does not reach a T1a channel".into(),
            flow_rule: "REVISION 1 (R1b): the flow axis wanders laterally by a sin(2 pi r / L + phi), a, L in the wander ranges, phi uniform (keyed t1b.volcano.flow counters 120 + 3k + 0..2). lava_flow_lobes_v1: flow_count flows, keyed t1b.volcano.flow (counter = blake3(id) x 256 + 16 + 5k + j): angle (k + u) / count of a turn, start r0 and length (fractions of R_f), half-width w (growing x1.5 downstream) and thickness t uniform in their ranges; cross-section t sqrt(1 - (lat / w)^2), 200 m ramp-in, semicircular snout of radius w; flows end inside R_end; added on the cone surface (max over flows). Recalled, unverified: andesite/dacite flows on stratovolcanoes ~10-100 m thick, 0.2-1 km wide, 1-10 km long".into(),
            gully_rule: "REVISION 1 (R1b): the gully axis wanders as the flows (keyed t1b.volcano.ygully counters 200 + 3k + 0..2). radial_gullies_young_v1: gully_count gullies, keyed t1b.volcano.ygully (counter = blake3(id) x 256 + 64 + 3k + j): angle (k + u) / count of a turn, start r0 (fraction of R_f), depth d_max; d(r) = d_max sin(pi (r - r0) / (R_f - r0)) on r0 < r < R_f, V section of half-width hw + hw_ppm x r; a declared sub-record closure (the record cannot resolve ages below 2 My; gullies tens of metres deep form on young cones within ~10^3-10^4 yr, recalled)".into(),
            fan_rule: "REVISION 1 (R1c): no fan from a curve whose first point is at or below sea level; a fan whose apex T0 depth (sea - bil) exceeds fan_subaqueous x shelf_depth is subaqueous (not built; 7 ky of supply cannot fill that water to the surface); fan_nocut default off; (R1d) the shelf grade applies only to T1a nodes whose whole downstream T1a chain lies at bil <= sea. fans_attached_v1: (a) no wave cut (cliff or beach lowering) within a fan's radius + pad of its apex (a river keeps its mouth); (b) t1a_shelf_grade_v1: T1a node levels seaward of the T0 zero contour (bil <= sea) are capped at sea - F(x) min(1, -bil / shelf_depth) (graded to the lowstand like every shelf channel; the perched T1a level there is a T1a input property), and T1a curves take the capped levels at both link ends (marine receivers: min(T1a level, base level))".into(),
            solve_rule: "damped_ef_solve_v1: per solved cell the ln E step (CO's 0.6 ln((T - L)/(M - L)), clamped) is multiplied by m_E and the F step (gain x land-fraction excess) by m_F; each multiplier starts at 1 and halves when that cell's residual (M - T, resp. land - bilinear land) changes sign between consecutive evaluations, floor 1/32".into(),
            cliff_rule: "cliff_profile_v1: on hard classes (k_mult <= cliff_k_max) d = exact Euclidean distance (Felzenszwalb-Huttenlocher transform) to the flooded sea: platform (sea - platform_depth, sea) for d < cliff_retreat x k_mult (CO's retreat), then talus at talus_slope (0.70, 35 deg, angle of repose) over talus_width, then a wall at wall_slope (2.75, 70 deg) inland; z = min(z, profile(d)); soft classes keep CO's beach rule with the same d".into(),
            ring_plain: None,
        }
    }
    pub fn validate(&self) -> R<()> {
        if self.cone_hw_permille < 0 || self.cone_wb_min_mm <= 0 || self.cone_wb_max_mm < self.cone_wb_min_mm { return Err("cf: cone scaling bounds".into()); }
        if self.bury_inner_permille > self.bury_outer_permille || self.softmax_k_mm < 0 || self.foot_taper_permille < 0 { return Err("cf: burial / softmax / taper".into()); }
        if self.flow_r0_permille[0] > self.flow_r0_permille[1] || self.flow_len_permille[0] > self.flow_len_permille[1] || self.gully_r0_permille[0] > self.gully_r0_permille[1] { return Err("cf: ranges".into()); }
        if self.solve_damp_floor_permille <= 0 || self.solve_damp_floor_permille > 1000 { return Err("cf: damping floor".into()); }
        if let Some(rp) = self.ring_plain.as_ref() { rp.validate()?; }
        Ok(())
    }
    pub fn softmax(&self, a: f64, b: f64) -> f64 {
        let k = self.softmax_k_mm as f64 / 1000.0;
        if k <= 0.0 { return a.max(b); }
        a.max(b) + k * (-(a - b).abs() / k).exp().ln_1p()
    }
}

/// Lane CF per-cone geometry (only when `CfParamsV1::cone_hw_permille` > 0 and the cone is not a residual).
#[derive(Clone, Debug, Serialize)]
pub struct CfVolcano {
    pub wb_m: f64,
    pub r_foot_m: f64,
    pub q_built: f64,
    pub lambda_m: f64,
    pub s_base: f64,
    pub r_end_max_m: f64,
    pub clearance_m: f64,
    /// (angle rad, r0 m, length m, half-width m, thickness m, wander amplitude m, wander length m, wander phase rad)
    pub flows: Vec<[f64; 8]>,
    /// (angle rad, r0 m, depth m, wander amplitude m, wander length m, wander phase rad)
    pub gullies: Vec<[f64; 6]>,
    pub vent_to_t1a_m: f64,
}

/// Native composite-flank closure. Coordinates and boundaries belong to the vent,
/// never to a requested window. Only the compact ledger goes into published JSON.
#[derive(Clone, Debug, Serialize)]
pub struct FlankSurface {
    pub rule: String,
    pub spacing_m: f64,
    pub n_radial: usize,
    pub n_angular: usize,
    pub duration_yr: f64,
    pub steps: usize,
    pub k_per_my: f64,
    pub runoff_source_area_m2: f64,
    pub runoff_outlet_area_m2: f64,
    pub runoff_relative_residual: f64,
    pub eroded_m3: f64,
    pub stream_eroded_m3: f64,
    pub creep_eroded_m3: f64,
    pub deposited_m3: f64,
    pub exported_m3: f64,
    pub volume_relative_residual: f64,
    pub max_incision_m: f64,
    pub confluences: usize,
    pub max_area_m2: f64,
    pub incision_gate_width_m: f64,
    pub max_incision_area_m2: f64,
    pub max_halfwidth_m: f64,
    pub outlet_export_m3: Vec<(u32, f64)>,
    #[serde(skip)]
    correction: std::sync::Arc<Vec<f64>>,
}

/// Preserved older lava emplacement provides topographic drainage divides.
/// The existing keyed CF flow sizes/locations are reused, with C1 margins.
/// This is a declared initial constructional deposit, not random height noise.
fn flank_lava(v: &VolcanoSite, x: f64, y: f64) -> f64 {
    let Some(c) = &v.cf else { return 0.0 };
    let r = x.hypot(y);
    if r <= v.crater_radius_m { return 0.0; }
    let th = y.atan2(x);
    let mut h: f64 = 0.0;
    for f in &c.flows {
        let d = ang_diff(th, f[0]);
        if d.abs() >= std::f64::consts::FRAC_PI_2 { continue; }
        let along = r * d.cos() - f[1];
        if along <= 0.0 || along >= f[2] { continue; }
        let off = f[5] * (std::f64::consts::TAU * r / f[6].max(1.0) + f[7]).sin();
        let w = f[3] * (1.0 + 0.5 * along / f[2]);
        let q = (r * d.sin() - off).abs() / w;
        if q >= 1.0 { continue; }
        let taper = smoothstep01(along / w) * smoothstep01((f[2] - along) / w);
        h = h.max(f[4] * (1.0 - q * q).powi(2) * taper);
    }
    h * ((v.h_m - v.denudation_m) / v.h_m).max(0.0)
        * smoothstep01((r - v.crater_radius_m) / (2.0 * v.crater_radius_m))
        * smoothstep01((c.r_foot_m - r) / (0.15 * c.r_foot_m))
}

struct FlankRouting {
    order: Vec<usize>,
    receivers: Vec<[(usize, f64, f64); 2]>,
    fill: Vec<f64>,
}

fn polar_neighbour(i: usize, nr: usize, nt: usize, dr: i64, da: i64) -> Option<usize> {
    let r = (i / nt) as i64 + dr;
    if r < 0 || r >= nr as i64 { return None; }
    Some(r as usize * nt + ((i % nt) as i64 + da).rem_euclid(nt as i64) as usize)
}

/// Native edifice-local priority flood on a periodic polar mesh. A centred
/// metric gradient supplies the continuous direction; angular partition to two
/// downhill neighbours follows Tarboton (1997). A flat takes its flood parent.
/// A radial surface has exactly zero angular derivative at EVERY azimuth.
fn polar_route(nr: usize, nt: usize, dx: f64, z: &[f64], fixed: &[bool]) -> FlankRouting {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let nn = nr * nt;
    let da = std::f64::consts::TAU / nt as f64;
    let ordered = |v: f64| { let b = v.to_bits(); if b >> 63 == 1 { !b } else { b ^ (1 << 63) } };
    let mut heap = BinaryHeap::new();
    let mut seen = vec![false; nn];
    let mut fill = z.to_vec();
    let mut parent = vec![usize::MAX; nn];
    for i in 0..nn { if fixed[i] { seen[i] = true; heap.push(Reverse((ordered(z[i]), i))); } }
    let mut order = Vec::with_capacity(nn);
    while let Some(Reverse((_, i))) = heap.pop() {
        order.push(i);
        for (dr, dt) in crate::t1b_lem::D8N {
            let Some(j) = polar_neighbour(i, nr, nt, dr, dt) else { continue };
            if seen[j] { continue; }
            seen[j] = true; fill[j] = z[j].max(fill[i]); parent[j] = i;
            heap.push(Reverse((ordered(fill[j]), j)));
        }
    }
    let mut pos = vec![0; nn];
    for (k,&i) in order.iter().enumerate() { pos[i] = k; }
    let mut receivers = vec![[(usize::MAX, 0.0, 0.0); 2]; nn];
    for i in 0..nn {
        if fixed[i] { continue; }
        let ri = i / nt;
        let r = ri as f64 * dx;
        let mut neighbours = [(usize::MAX, 0.0, 0.0, 0.0); 8];
        let mut best = -1.0_f64;
        for (k, (dr, dt)) in crate::t1b_lem::D8N.iter().copied().enumerate() {
            let Some(j) = polar_neighbour(i, nr, nt, dr, dt) else { continue };
            let rr = (ri as i64 + dr) as f64 * dx;
            let (vx, vy) = (rr * (dt as f64 * da).cos() - r, rr * (dt as f64 * da).sin());
            let len = vx.hypot(vy);
            neighbours[k] = (j, vx, vy, len);
            if pos[j] < pos[i] && len > 0.0 && fill[j] <= fill[i] {
                let slope = (fill[i] - fill[j]) / len;
                if slope > best { best = slope; receivers[i][0] = (j, 1.0, len); }
            }
        }
        // On a filled flat, use the actual flood tree, not an arbitrary equal-height neighbour.
        if best <= 0.0 {
            let j = parent[i];
            let len = neighbours.iter().find(|q| q.0 == j).expect("polar flood parent is adjacent").3;
            receivers[i][0] = (j, 1.0, len);
            continue;
        }
        let rp = polar_neighbour(i,nr,nt,1,0).unwrap();
        let rm = polar_neighbour(i,nr,nt,-1,0).unwrap();
        let ap = polar_neighbour(i,nr,nt,0,1).unwrap();
        let am = polar_neighbour(i,nr,nt,0,-1).unwrap();
        let (ux,uy) = (-(fill[rp]-fill[rm])/(2.0*dx), -(fill[ap]-fill[am])/(2.0*r*da));
        if ux == 0.0 && uy == 0.0 { continue; }
        for k in 0..8 {
            let (j,ax,ay,al) = neighbours[k];
            let (q,bx,by,bl) = neighbours[(k+1)%8];
            if j == usize::MAX || q == usize::MAX || pos[j]>=pos[i] || pos[q]>=pos[i] || fill[j]>=fill[i] || fill[q]>=fill[i] { continue; }
            if ax*uy-ay*ux < 0.0 || ux*by-uy*bx < 0.0 { continue; }
            let total = (ax*by-ay*bx).atan2(ax*bx+ay*by);
            let w = ((ax*uy-ay*ux).atan2(ax*ux+ay*uy)/total).clamp(0.0,1.0);
            receivers[i] = [(j,1.0-w,al),(q,w,bl)]; break;
        }
    }
    FlankRouting { order, receivers, fill }
}

fn flank_tridiagonal(a:&[f64],b:&[f64],c:&[f64],d:&[f64])->Vec<f64> {
    let n=d.len();let mut cp=vec![0.0;n];let mut x=vec![0.0;n];
    cp[0]=c[0]/b[0];x[0]=d[0]/b[0];
    for i in 1..n {let den=b[i]-a[i]*cp[i-1];cp[i]=c[i]/den;x[i]=(d[i]-a[i]*x[i-1])/den;}
    for i in (0..n-1).rev(){x[i]-=cp[i]*x[i+1];}x
}

/// Conservative finite-volume hillslope transport, implicit radial then
/// periodic angular solves. Zero normal flux at the protected crater/foot.
/// kappa is m2/yr. Each internal edge uses the same diffusivity and flux in
/// both adjacent cells, including the r Jacobian; no slope clipping or loss.
fn flank_diffuse(nr:usize,nt:usize,dx:f64,z:&mut [f64],fixed:&[bool],kappa:&[f64],dt_yr:f64) {
    let da=std::f64::consts::TAU/nt as f64;
    let edge=|i:usize,j:usize|if fixed[i]||fixed[j]{0.0}else{0.5*(kappa[i]+kappa[j])};
    for t in 0..nt {
        let mut a=vec![0.0;nr];let mut b=vec![1.0;nr];let mut c=vec![0.0;nr];let mut d=vec![0.0;nr];
        for r in 0..nr {
            let i=r*nt+t;d[r]=z[i];if fixed[i]{continue;}
            let rr=r as f64*dx;
            let cm=dt_yr*edge(i,i-nt)*(rr-0.5*dx)/(rr*dx*dx);
            let cp=dt_yr*edge(i,i+nt)*(rr+0.5*dx)/(rr*dx*dx);
            a[r]=-cm;b[r]=1.0+cm+cp;c[r]=-cp;
        }
        let x=flank_tridiagonal(&a,&b,&c,&d);for r in 0..nr{z[r*nt+t]=x[r];}
    }
    for r in 0..nr {
        if fixed[r*nt]{continue;}
        let rr=r as f64*dx;let factor=dt_yr/(rr*rr*da*da);
        let mut a=vec![0.0;nt];let mut b=vec![1.0;nt];let mut c=vec![0.0;nt];let mut d=vec![0.0;nt];
        for t in 0..nt {
            let i=r*nt+t;d[t]=z[i];let cm=factor*edge(i,r*nt+(t+nt-1)%nt);let cp=factor*edge(i,r*nt+(t+1)%nt);
            a[t]=-cm;b[t]=1.0+cm+cp;c[t]=-cp;
        }
        // Sherman-Morrison cyclic correction. alpha=A[0,n-1], beta=A[n-1,0].
        let (alpha,beta,gamma)=(a[0],c[nt-1],-b[0]);b[0]-=gamma;b[nt-1]-=alpha*beta/gamma;
        a[0]=0.0;c[nt-1]=0.0;
        let mut x=flank_tridiagonal(&a,&b,&c,&d);let mut u=vec![0.0;nt];u[0]=gamma;u[nt-1]=beta;
        let q=flank_tridiagonal(&a,&b,&c,&u);
        let correction=(x[0]+alpha*x[nt-1]/gamma)/(1.0+q[0]+alpha*q[nt-1]/gamma);
        for t in 0..nt{x[t]-=correction*q[t];z[r*nt+t]=x[t];}
    }
}

impl FlankSurface {
    pub(crate) fn build(v: &VolcanoSite) -> Self {
        // Exact physical input bits key this bounded, process-local memoization.
        thread_local! { static CACHE: std::cell::RefCell<BTreeMap<Vec<u64>, FlankSurface>> = std::cell::RefCell::new(BTreeMap::new()); }
        let cf = v.cf.as_ref().unwrap();
        let mut key: Vec<u64> = [v.h_m,v.denudation_m,v.summit_built(),v.s_top,v.lambda_m,
            v.crater_radius_m,v.crater_depth_m,cf.r_foot_m].iter().map(|v| v.to_bits()).collect();
        key.extend(cf.flows.iter().flatten().map(|v| v.to_bits()));
        if let Some(f) = CACHE.with(|c| c.borrow().get(&key).cloned()) { return f; }
        let f = Self::build_at(v,64.0,60,20.0);
        CACHE.with(|c| { let mut c=c.borrow_mut(); if c.len()>=64 {c.pop_first();} c.insert(key,f.clone()); });
        f
    }

    fn build_at(v: &VolcanoSite, dx: f64, steps: usize, k_per_my: f64) -> Self {
        let rf = v.cf.as_ref().unwrap().r_foot_m;
        let nr = (rf/dx).ceil() as usize + 4;
        let nt = (std::f64::consts::TAU*rf/dx).ceil() as usize;
        let nn = nr*nt;
        let da = std::f64::consts::TAU/nt as f64;
        let rc = v.crater_radius_m;
        let scale = (v.h_m-v.denudation_m)/v.h_m;
        let top = v.summit_built()-v.denudation_m;
        let mut z=vec![0.0;nn];let mut fixed=vec![true;nn];let mut taper=vec![0.0;nn];
        let mut cell_area=vec![0.0;nn];let mut xy=vec![[0.0;2];nn];
        for i in 0..nn {
            let r=(i/nt) as f64*dx;let a=(i%nt) as f64*da;
            xy[i]=[r*a.cos(),r*a.sin()];cell_area[i]=r*dx*da;
            z[i]=if r<rc {top-v.crater_depth_m*scale*(1.0-(r/rc).powi(2))}
                else {top-v.s_top*v.lambda_m*(1.0-(-(r-rc)/v.lambda_m).exp())*scale};
            z[i]+=flank_lava(v,xy[i][0],xy[i][1]);
            if r>rc+192.0 && r<rf-192.0 {
                fixed[i]=false;taper[i]=smoothstep01((r-rc-192.0)/128.0)*smoothstep01((rf-192.0-r)/128.0);
            }
        }
        let initial=z.clone();let source=(0..nn).filter(|&i| !fixed[i]).map(|i|cell_area[i]).sum::<f64>();
        let duration_yr=30_000.0;let dt=duration_yr/1e6/steps as f64;
        let mut exports=BTreeMap::<u32,f64>::new();
        let mut max_incision_area=0.0_f64;
        let (mut creep_eroded,mut deposited)=(0.0_f64,0.0_f64);
        let (mut eroded,mut water_residual,mut sink,mut max_area,mut joins)=(0.0_f64,0.0_f64,0.0,0.0_f64,0);
        for _ in 0..steps {
            let routing=polar_route(nr,nt,dx,&z,&fixed);
            let mut area:Vec<f64>=(0..nn).map(|i|if fixed[i] {0.0}else{cell_area[i]}).collect();
            let mut donors=vec![0u32;nn];
            for &i in routing.order.iter().rev() {
                for &(j,w,_) in &routing.receivers[i] {if w>0.0 {area[j]+=w*area[i];if w*area[i]>=100_000.0 {donors[j]+=1;}}}
            }
            sink=(0..nn).filter(|&i|fixed[i]).map(|i|area[i]).sum();
            water_residual=water_residual.max((sink-source).abs()/source.max(1.0));
            max_area=area.iter().copied().fold(0.0,f64::max);joins=donors.iter().filter(|&&n|n>1).count();
            // Integrate flux through a fixed 128m radial gate, not one changing mesh cell.
            // Fractions of end cells are included. This is the erosion support scale [H],
            // not an alteration to the routed water or its outlet conservation ledger.
            let radial_flux: Vec<f64>=(0..nn).map(|i|area[i]*routing.receivers[i].iter()
                .filter(|&&(j,w,_)|w>0.0 && j/nt!=i/nt).map(|&(_,w,_)|w).sum::<f64>()).collect();
            let mut incision_area=vec![0.0;nn];
            for i in 0..nn {if !fixed[i] {
                let arc=(i/nt) as f64*dx*da;let reach=(64.0/arc+0.5).ceil() as i64;
                for a in -reach..=reach {
                    let overlap=((a as f64+0.5)*arc).min(64.0)-((a as f64-0.5)*arc).max(-64.0);
                    if overlap>0.0 {let j=polar_neighbour(i,nr,nt,0,a).unwrap();incision_area[i]+=radial_flux[j]*overlap/arc;}
                }
            }}
            max_incision_area=incision_area.iter().copied().fold(0.0,f64::max);
            let mut cut=vec![0.0_f64;nn];let mut bed=z.clone();
            for &i in &routing.order {
                if fixed[i] || z[i]<routing.fill[i]-1e-9 {continue;}
                let (mut inv_len,mut base)=(0.0,0.0);
                for &(j,w,len) in &routing.receivers[i] {if w>0.0 {
                    inv_len+=w/len;base+=w/len*if z[j]<routing.fill[j]-1e-9 {routing.fill[j]}else{bed[j]};
                }}
                base/=inv_len;let f=k_per_my*incision_area[i].sqrt()*dt*inv_len
                    *smoothstep01((incision_area[i]-50_000.0)/50_000.0)*taper[i];
                let depth=((z[i]-base)*f/(1.0+f)).max(0.0);
                bed[i]=z[i]-depth;if depth<=0.0 {continue;}
                cut[i]=depth;
            }
            let mut sediment=vec![0.0;nn];
            for i in 0..nn {z[i]-=cut[i];sediment[i]=cut[i]*cell_area[i];eroded+=sediment[i];}
            for &i in routing.order.iter().rev() {
                if fixed[i] {if sediment[i]>0.0 {*exports.entry(i as u32).or_default()+=sediment[i];}}
                else {for &(j,w,_) in &routing.receivers[i] {if w>0.0 {sediment[j]+=w*sediment[i];}}}
            }
            let before_creep=z.clone();
            // sqrt(2*kappa*T)=b(A): declared bank-mobility closure, not a fitted soil diffusivity.
            let kappa:Vec<f64>=(0..nn).map(|i| if fixed[i]{0.0}else{
                Self::halfwidth(incision_area[i]).powi(2)/(2.0*duration_yr)*(k_per_my/20.0)*taper[i]
            }).collect();
            flank_diffuse(nr,nt,dx,&mut z,&fixed,&kappa,duration_yr/steps as f64);
            for i in 0..nn {let dv=(before_creep[i]-z[i])*cell_area[i];creep_eroded+=dv.max(0.0);deposited+=(-dv).max(0.0);}
        }
        let correction:Vec<f64>=z.iter().zip(&initial).map(|(a,b)|a-b).collect();
        let removed=-correction.iter().zip(&cell_area).map(|(z,a)|z*a).sum::<f64>();
        let exported=exports.values().sum::<f64>();
        Self {rule:crate::t1b_ve::FLANK_RULE.into(),spacing_m:dx,n_radial:nr,n_angular:nt,duration_yr,steps,k_per_my,
            runoff_source_area_m2:source,runoff_outlet_area_m2:sink,runoff_relative_residual:water_residual,
            eroded_m3:eroded+creep_eroded,stream_eroded_m3:eroded,creep_eroded_m3:creep_eroded,deposited_m3:deposited,exported_m3:exported,volume_relative_residual:(removed-exported).abs().max((eroded+creep_eroded-deposited-exported).abs())/(eroded+creep_eroded).max(1.0),
            max_incision_m:correction.iter().map(|z|-z).fold(0.0,f64::max),confluences:joins,max_area_m2:max_area,incision_gate_width_m:128.0,max_incision_area_m2:max_incision_area,max_halfwidth_m:Self::halfwidth(max_incision_area),
            outlet_export_m3:exports.into_iter().collect(),correction:std::sync::Arc::new(correction)}
    }

    fn halfwidth(area_m2:f64)->f64 {32.0+80.0*(area_m2/1e6).powf(0.4)}

    /// Periodic C2 B-spline reconstruction. Partition of unity and its first
    /// radial moment preserve the polar r dr dtheta nodal erosion integral.
    fn at(&self,x:f64,y:f64)->f64 {
        let u=x.hypot(y)/self.spacing_m;
        let v=y.atan2(x).rem_euclid(std::f64::consts::TAU)/std::f64::consts::TAU*self.n_angular as f64;
        let (ir,it)=(u.floor() as i64,v.floor() as i64);
        let basis=|t:f64|[(1.0-t).powi(3)/6.0,(3.0*t.powi(3)-6.0*t*t+4.0)/6.0,
            (-3.0*t.powi(3)+3.0*t*t+3.0*t+1.0)/6.0,t.powi(3)/6.0];
        let (wr,wt)=(basis(u-ir as f64),basis(v-it as f64));let mut h=0.0;
        for j in 0..4 {for i in 0..4 {
            let rr=ir+j as i64-1;let tt=(it+i as i64-1).rem_euclid(self.n_angular as i64);
            if rr>=0 && rr<self.n_radial as i64 {h+=wr[j]*wt[i]*self.correction[rr as usize*self.n_angular+tt as usize];}
        }}h
    }
}

/// Lane HS-PASS hs_volcano: level 1 empties lane CF's radial_gullies_young_v1 list of every cone (gully_count 0); level 2 also
/// the record_age_degradation_v1 gully directions. Returns (young gullies removed, record gullies removed). Draws are untouched.
pub fn hs_strip_gullies(vs: &mut [VolcanoSite], level: u32) -> (usize, usize) {
    let (mut a, mut b) = (0, 0);
    for v in vs.iter_mut() {
        if let Some(c) = v.cf.as_mut() { a += c.gullies.len(); c.gullies.clear(); }
        if level >= 2 { b += v.gully_dirs_rad.len(); v.gully_dirs_rad.clear(); }
    }
    (a, b)
}

/// lambda with s_top lambda (1 - exp(-l / lambda)) = h (bisection; None when s_top l <= h).
pub fn solve_lambda(s_top: f64, l: f64, h: f64) -> Option<f64> {
    if s_top * l <= h * (1.0 + 1e-9) { return None; }
    let f = |lam: f64| s_top * lam * (1.0 - (-l / lam).exp()) - h;
    let (mut a, mut b) = (1e-3_f64, 1e9_f64);
    for _ in 0..200 { let m = (a * b).sqrt(); if f(m) > 0.0 { b = m; } else { a = m; } }
    Some((a * b).sqrt())
}

#[inline]
pub(crate) fn smoothstep01(u: f64) -> f64 { let u = u.clamp(0.0, 1.0); u * u * (3.0 - 2.0 * u) }

/// Lane CF cone at chart (x, y) m: (E = cone surface with flows and gullies, burial weight w, foot fade u).
pub fn cone_cf(vp: &VolcanoParamsV1, v: &VolcanoSite, cf: &CfParamsV1, x: f64, y: f64, bil_xy: f64) -> Option<(f64, f64, f64)> {
    cone_cf_opt(vp, v, cf, x, y, bil_xy, true)
}

/// As `cone_cf`; `lobes_gullies` false leaves out CF's lava lobes and gullies (lane VE builds its own units).
pub fn cone_cf_opt(vp: &VolcanoParamsV1, v: &VolcanoSite, cf: &CfParamsV1, x: f64, y: f64, bil_xy: f64, lobes_gullies: bool) -> Option<(f64, f64, f64)> {
    let g = v.cf.as_ref()?;
    let none: Vec<[f64; 8]> = Vec::new();
    let none6: Vec<[f64; 6]> = Vec::new();
    let (flows, gullies) = if lobes_gullies && v.ve.as_ref().is_none_or(|g| g.flank.is_none() && g.flank2.is_none()) { (&g.flows, &g.gullies) } else { (&none, &none6) };
    let z0 = match edifice_at(vp, v, x, y) { EdificeAt::Abs(z, _) => z, _ => return None };
    let (dx, dy) = (x - v.vent_xy_m[0], y - v.vent_xy_m[1]);
    let r = dx.hypot(dy);
    let th = dy.atan2(dx);
    let mut add: f64 = 0.0;
    for f in flows {
        let d = ang_diff(th, f[0]);
        if d.abs() >= std::f64::consts::FRAC_PI_2 { continue; }
        let off = f[5] * (std::f64::consts::TAU * r / f[6].max(1.0) + f[7]).sin();
        let (along, lat) = (r * d.cos(), (r * d.sin() - off).abs());
        let s = along - f[1];
        if s < 0.0 || s > f[2] { continue; }
        let w = f[3] * (1.0 + 0.5 * s / f[2]);
        let wl = if s > f[2] - w { let q = ((s - (f[2] - w)) / w).min(1.0); w * (1.0 - q * q).max(0.0).sqrt() } else { w };
        if lat >= wl || wl <= 0.0 { continue; }
        let ramp = (s / 200.0).min(1.0);
        add = add.max(f[4] * (1.0 - (lat / wl).powi(2)).sqrt() * ramp);
    }
    let mut cut: f64 = 0.0;
    let rf = g.r_foot_m;
    for q in gullies {
        if r <= q[1] || r >= rf { continue; }
        let d = ang_diff(th, q[0]);
        if d.abs() >= std::f64::consts::FRAC_PI_2 { continue; }
        let off = q[3] * (std::f64::consts::TAU * r / q[4].max(1.0) + q[5]).sin();
        let lat = (r * d.sin() - off).abs();
        let hw = cf.gully_halfwidth_mm as f64 / 1000.0 + cf.gully_halfwidth_ppm_of_r as f64 / 1e6 * r;
        if lat >= hw { continue; }
        let depth = q[2] * (std::f64::consts::PI * (r - q[1]) / (rf - q[1])).sin();
        cut = cut.max(depth * (1.0 - lat / hw));
    }
    // R1a; fix B3: C1-smooth min (0 for t >= 0, -t^2 / 4k down to t = -2k, then t + k; k = 100 m), so the bilinear contour
    // bil = base (a near-straight line inside a T0 cell) leaves no crease.
    let sb = if cf.sloped_base == 1 { let t = bil_xy - v.base_m; let k = 100.0; if t >= 0.0 { 0.0 } else if t > -2.0 * k { -t * t / (4.0 * k) } else { t + k } } else { 0.0 };
    let e = z0 + add - cut + sb;
    let (bi, bo) = (cf.bury_inner_permille as f64 / 1000.0 * rf, cf.bury_outer_permille as f64 / 1000.0 * rf);
    let w = if r <= bi { 1.0 } else if r >= bo { 0.0 } else { smoothstep01((bo - r) / (bo - bi)) };
    let lob = if v.lobe_phase_rad.is_empty() { 0.0 } else { v.lobe_phase_rad.iter().enumerate().map(|(j, p)| ((j + 2) as f64 * th + p).cos() / (j + 2) as f64).sum::<f64>() / v.lobe_norm };
    let r_end = rf * (1.0 + cf.foot_taper_permille as f64 / 1000.0 * (1.0 + 0.8 * lob));
    let u = if r <= rf { 1.0 } else if r_end <= rf + 1.0 { 0.0 } else { smoothstep01((r_end - r) / (r_end - rf)) };
    Some((e, w, u))
}

/// Exact Euclidean distance (metres) from `src` posts (Felzenszwalb and Huttenlocher 2012 separable transform).
pub fn edt(lat: &EvalLat, src: &dyn Fn(usize) -> bool, cap: f64) -> Vec<f64> {
    let (n, dx) = (lat.n, lat.dx());
    let big = 1e20_f64;
    let mut f: Vec<f64> = (0..n * n).map(|i| if src(i) { 0.0 } else { big }).collect();
    let mut tmp = vec![0f64; n];
    let mut v = vec![0usize; n];
    let mut zz = vec![0f64; n + 1];
    let mut line = vec![0f64; n];
    let mut pass = |get: &dyn Fn(usize) -> usize, f: &mut Vec<f64>| {
        for k in 0..n { line[k] = f[get(k)]; }
        // 1D lower envelope of parabolas over the finite sites.
        let mut kk: isize = -1;
        for q in 0..n {
            if line[q] >= big { continue; }
            if kk < 0 { kk = 0; v[0] = q; zz[0] = f64::NEG_INFINITY; zz[1] = f64::INFINITY; continue; }
            loop {
                let p = v[kk as usize];
                let s = ((line[q] + (q * q) as f64) - (line[p] + (p * p) as f64)) / (2.0 * q as f64 - 2.0 * p as f64);
                if s <= zz[kk as usize] { kk -= 1; if kk < 0 { kk = 0; v[0] = q; zz[0] = f64::NEG_INFINITY; zz[1] = f64::INFINITY; break; } continue; }
                kk += 1; v[kk as usize] = q; zz[kk as usize] = s; zz[kk as usize + 1] = f64::INFINITY; break;
            }
        }
        if kk < 0 { return; }
        let mut j = 0usize;
        for q in 0..n {
            while zz[j + 1] < q as f64 { j += 1; }
            let p = v[j];
            tmp[q] = (q as f64 - p as f64).powi(2) + line[p];
        }
        for k in 0..n { f[get(k)] = tmp[k]; }
    };
    for y in 0..n { pass(&|k| y * n + k, &mut f); }
    for x in 0..n { pass(&|k| k * n + x, &mut f); }
    f.iter().map(|&s| { let d = if s >= big { f64::INFINITY } else { s.sqrt() * dx }; if d > cap { f64::INFINITY } else { d } }).collect()
}

// ------------------------------------------------------------------ volcano records

#[derive(Clone, Debug, Serialize)]
pub struct VolcanoFeature { pub id: String, pub kind: String, pub lat_deg: f64, pub lon_deg: f64, pub last_eruption_my: f64, pub subgrid_summit_m: f64, pub cell_elevation_m: f64,
    /// Lane VE: the record's onset (My), when present (never serialized).
    #[serde(skip)]
    pub onset_my: Option<f64> }

pub struct VolcanoRecords { pub t_my: f64, pub features: Vec<VolcanoFeature>, pub digest: String }

impl VolcanoRecords {
    pub fn load(path: &Path) -> R<VolcanoRecords> {
        let b = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let digest = format!("blake3:{}", blake3::hash(&b).to_hex());
        let v: Value = serde_json::from_slice(&b).map_err(|e| e.to_string())?;
        let t_my = v.get("t_my").and_then(Value::as_f64).ok_or("volcanoes: t_my")?;
        let mut features = Vec::new();
        for f in v.get("features").and_then(Value::as_array).ok_or("volcanoes: features")? {
            let s = |k: &str| f.get(k).and_then(Value::as_str).map(str::to_string).ok_or(format!("volcano {k}"));
            let n = |k: &str| f.get(k).and_then(Value::as_f64).ok_or(format!("volcano {k}"));
            features.push(VolcanoFeature { id: s("id")?, kind: s("kind")?, lat_deg: n("lat_deg")?, lon_deg: n("lon_deg")?, last_eruption_my: n("last_eruption_my")?, subgrid_summit_m: n("subgrid_summit_m")?, cell_elevation_m: n("cell_elevation_m")?,
                onset_my: f.get("onset_my").and_then(Value::as_f64) });
        }
        Ok(VolcanoRecords { t_my, features, digest })
    }
}

/// One built edifice (chart metres).
#[derive(Clone, Debug, Serialize)]
pub struct VolcanoSite {
    pub id: String,
    pub kind: String,
    pub shape: String,
    pub record_xy_m: [f64; 2],
    pub vent_xy_m: [f64; 2],
    pub moved_m: f64,
    pub record_cell: u32,
    pub vent_cell: u32,
    pub summit_record_m: f64,
    pub base_m: f64,
    pub h_m: f64,
    pub lambda_m: f64,
    pub r_main_m: f64,
    pub r_influence_m: f64,
    pub s_top: f64,
    pub s_base: f64,
    pub crater_radius_m: f64,
    pub crater_depth_m: f64,
    pub age_my: f64,
    pub e_vent_m_per_my: f64,
    pub denudation_m: f64,
    pub stage: String,
    pub breached: bool,
    pub breach_dir_rad: f64,
    pub gully_dirs_rad: Vec<f64>,
    pub dike_dirs_rad: Vec<f64>,
    pub neck_relief_m: f64,
    /// B4 (exploratory): foot lobes, phases for k = 2..6 and the normalisation of their sum.
    pub lobe_phase_rad: Vec<f64>,
    pub lobe_norm: f64,
    /// Lane CF (0.9.2): scaled cone, flows and gullies (absent in 0.9.1 layers).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cf: Option<CfVolcano>,
    /// Lane VE: eruption style and keyed history (absent without ve).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ve: Option<crate::t1b_ve::VeVolcano>,
    /// Lane T1B-PEDESTAL edifice_local_base_v1 (absent without the key): bilinear base at the vent minus base_m (the local
    /// ground G); the CF sloped base keeps t = B(post) - B(vent) by subtracting it from the plain bilinear.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ped_shift_m: Option<f64>,
    /// Lane VOLCANIC-SPLIT edifice_local_relief_v1 (absent without ped_relief): the built summit G + H, below the record.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub built_summit_m: Option<f64>,
}

impl VolcanoSite {
    /// Lane T1B-PEDESTAL: the shift of the base below the bilinear at the vent (0 without the key).
    #[inline] pub fn ped_shift(&self) -> f64 { self.ped_shift_m.unwrap_or(0.0) }
    /// The summit the edifice is built to: the record summit, or (ped_relief) G + H.
    #[inline] pub fn summit_built(&self) -> f64 { self.built_summit_m.unwrap_or(self.summit_record_m) }
}

/// What an edifice does at a point.
/// Abs(z, fade): edifice surface and the foot fade (1 inside R, smoothstep to 0 at fade_ratio x R).
pub enum EdificeAt { None, Abs(f64, f64), Add(f64) }

pub(crate) fn seg_nearest(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0) } else { 0.0 };
    [a[0] + t * dx, a[1] + t * dy]
}

pub(crate) fn bil_at(t0: &T0, chart: &Chart, x: f64, y: f64) -> R<f64> {
    let (la, lo) = chart.inverse_f64(x * 1000.0, y * 1000.0).map_err(|e| e.to_string())?;
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    Ok(crate::t1b::bilinear(t0, &elev, la * RAD2DEG, lo * RAD2DEG, false))
}

/// Hat-blended per-cell table value at chart (x, y) metres.
pub fn table_at(t0: &T0, chart: &Chart, table: &BTreeMap<u32, f64>, x: f64, y: f64) -> R<f64> {
    let (la, lo) = chart.inverse_f64(x * 1000.0, y * 1000.0).map_err(|e| e.to_string())?;
    let (c, w) = hat4(t0, la * RAD2DEG, lo * RAD2DEG);
    let mut v = 0.0;
    for k in 0..4 { if w[k] > 0.0 { v += w[k] * table.get(&c[k]).copied().unwrap_or(0.0); } }
    Ok(v)
}

/// Build the edifices of every record point inside the chart. `age_override_my` exists for unit tests only (the CLI
/// never sets it; design 2026-10-13: no presentation age override).
#[allow(clippy::too_many_arguments)]
pub fn volcano_sites(vp: &VolcanoParamsV1, recs: &VolcanoRecords, t0: &T0, chart: &Chart, curves: &[Curve], table: &BTreeMap<u32, f64>, seed: u64, recipe_digest: &str, realization: u32, age_override_my: Option<f64>, cf: Option<&CfParamsV1>) -> R<Vec<VolcanoSite>> {
    volcano_sites_ve(vp, recs, t0, chart, curves, table, seed, recipe_digest, realization, age_override_my, cf, None)
}

/// As `volcano_sites`; with lane VE params the broad families (shield, rift field, plateau) take lane CO's clearance plus the
/// shield clearance (`t1b_ve::broad_clearance`) instead of lane CF's foot clearance (composites unchanged).
#[allow(clippy::too_many_arguments)]
pub fn volcano_sites_ve(vp: &VolcanoParamsV1, recs: &VolcanoRecords, t0: &T0, chart: &Chart, curves: &[Curve], table: &BTreeMap<u32, f64>, seed: u64, recipe_digest: &str, realization: u32, age_override_my: Option<f64>, cf: Option<&CfParamsV1>, ve: Option<&crate::t1b_ve::VeParamsV1>) -> R<Vec<VolcanoSite>> {
    volcano_sites_ped(vp, recs, t0, chart, curves, table, seed, recipe_digest, realization, age_override_my, cf, ve, None)
}

/// As `volcano_sites_ve`; lane T1B-PEDESTAL edifice_local_base_v1: `ground` gives the local ground G per volcano id (the
/// base of the eligible composite cones; the record summit is kept, so H = summit - G; the vent placement is unchanged).
#[allow(clippy::too_many_arguments)]
pub fn volcano_sites_ped(vp: &VolcanoParamsV1, recs: &VolcanoRecords, t0: &T0, chart: &Chart, curves: &[Curve], table: &BTreeMap<u32, f64>, seed: u64, recipe_digest: &str, realization: u32, age_override_my: Option<f64>, cf: Option<&CfParamsV1>, ve: Option<&crate::t1b_ve::VeParamsV1>, ground: Option<&BTreeMap<String, f64>>) -> R<Vec<VolcanoSite>> {
    volcano_sites_ped2(vp, recs, t0, chart, curves, table, seed, recipe_digest, realization, age_override_my, cf, ve, ground, false)
}

/// As `volcano_sites_ped`; lane VOLCANIC-SPLIT edifice_local_relief_v1 (`relief`): an eligible cone on its local ground G
/// keeps the record relief H = record summit - bilinear base (built summit G + H) instead of the record summit.
#[allow(clippy::too_many_arguments)]
pub fn volcano_sites_ped2(vp: &VolcanoParamsV1, recs: &VolcanoRecords, t0: &T0, chart: &Chart, curves: &[Curve], table: &BTreeMap<u32, f64>, seed: u64, recipe_digest: &str, realization: u32, age_override_my: Option<f64>, cf: Option<&CfParamsV1>, ve: Option<&crate::t1b_ve::VeParamsV1>, ground: Option<&BTreeMap<String, f64>>, relief: bool) -> R<Vec<VolcanoSite>> {
    let half = vz_tilekey::CHART_HALF_EXTENT_MM as f64 / 1000.0;
    let layer = Scope::root(seed).layer(recipe_digest).map_err(|e| e.to_string())?.chart(0).map_err(|e| e.to_string())?;
    let t1a: Vec<&Curve> = curves.iter().filter(|c| c.code == CHANNEL).collect();
    let nearest_t1a = |p: [f64; 2]| -> Option<([f64; 2], [f64; 2], f64)> {
        let mut best: Option<([f64; 2], [f64; 2], f64)> = None;
        for c in &t1a {
            for w in c.pts.windows(2) {
                let (a, b) = ([w[0][0], w[0][1]], [w[1][0], w[1][1]]);
                let q = seg_nearest(p, a, b);
                let d = (p[0] - q[0]).hypot(p[1] - q[1]);
                if best.map(|bb| d < bb.2).unwrap_or(true) { best = Some((q, [b[0] - a[0], b[1] - a[1]], d)); }
            }
        }
        best
    };
    let mut out = Vec::new();
    for f in &recs.features {
        let (la, lo) = (f.lat_deg.to_radians(), f.lon_deg.to_radians());
        let Ok([xm, ym]) = chart.forward_f64(la, lo) else { continue };
        let (x, y) = (xm / 1000.0, ym / 1000.0);
        if x.abs() > half || y.abs() > half { continue; }
        let Some(shape) = vp.shapes.iter().find(|s| s.kinds.iter().any(|k| *k == f.kind)) else { continue };
        // Placement (REVISION 1: clearance r_c + relief / cap_slope; B4: relief actually built, from the record point).
        let h0 = f.subgrid_summit_m - bil_at(t0, chart, x, y)?;
        let age0 = age_override_my.unwrap_or(recs.t_my - f.last_eruption_my).max(0.0);
        let resid0 = vp.erosion_mult_permille as f64 / 1000.0 * crate::t1b_ve::erosion_at(ve, t0, chart, table, x, y)? * age0 >= h0;
        let built0 = if resid0 { vp.neck_relief_permille as f64 / 1000.0 * h0 } else { h0 };
        let rr0 = if resid0 { vp.neck_radius_mm as f64 / 1000.0 } else { shape.crater_radius_mm as f64 / 1000.0 };
        let mut clear = (vp.clear_t1a_mm as f64 / 1000.0).max(rr0 + built0.max(0.0) / (vp.cap_slope_ppm as f64 / 1e6));
        // Lane VE: broad families clear by the shield rule (and skip CF's foot clearance below).
        let ve_broad = if ve.is_some() && !resid0 { crate::t1b_ve::broad_clearance(&f.kind, h0.max(0.0)) } else { None };
        if let Some(c) = ve_broad { clear = clear.max(c); }
        let cf_c3 = cf.filter(|_| ve_broad.is_none());
        // Lane CF C3: the foot (R_f from the record point's H) clears the T1a channels.
        if let Some(cf) = cf_c3 { if !resid0 && cf.cone_hw_permille > 0 && cf.clear_foot_permille > 0 && h0 > 0.0 {
            let wb = (h0 / (cf.cone_hw_permille as f64 / 1000.0)).clamp(cf.cone_wb_min_mm as f64 / 1000.0, cf.cone_wb_max_mm as f64 / 1000.0);
            clear = clear.max(cf.clear_foot_permille as f64 / 1000.0 * wb / 2.0);
        } }
        let mut v = [x, y];
        let clear0 = clear;
        for _ in 0..4 {
            // REVISION 1 (R1e): R_f from the H at the current vent position.
            if let Some(cf) = cf_c3 { if !resid0 && cf.cone_hw_permille > 0 && cf.clear_foot_permille > 0 {
                let hv = f.subgrid_summit_m - bil_at(t0, chart, v[0], v[1])?;
                if hv > 0.0 {
                    let wb = (hv / (cf.cone_hw_permille as f64 / 1000.0)).clamp(cf.cone_wb_min_mm as f64 / 1000.0, cf.cone_wb_max_mm as f64 / 1000.0);
                    clear = clear0.max(cf.clear_foot_permille as f64 / 1000.0 * wb / 2.0);
                }
            } }
            let Some((q, dir, d)) = nearest_t1a(v) else { break };
            if d >= clear - 1e-6 { break; }
            let n = if d > 1.0 { [(v[0] - q[0]) / d, (v[1] - q[1]) / d] } else {
                let l = dir[0].hypot(dir[1]).max(1e-9);
                let nn = [-dir[1] / l, dir[0] / l];
                let up = bil_at(t0, chart, q[0] + 100.0 * nn[0], q[1] + 100.0 * nn[1])? >= bil_at(t0, chart, q[0] - 100.0 * nn[0], q[1] - 100.0 * nn[1])?;
                if up { nn } else { [-nn[0], -nn[1]] }
            };
            v = [q[0] + clear * n[0], q[1] + clear * n[1]];
        }
        let base = bil_at(t0, chart, v[0], v[1])?;
        // Lane T1B-PEDESTAL edifice_local_base_v1: the eligible cone's base is the published local ground.
        let (base, ped_shift) = match ground.and_then(|m| m.get(&f.id)) {
            Some(&g) if crate::t1b_ped::base_applies(&f.kind, base, t0.sea_level_m, cf.map(|c| c.cone_hw_permille > 0).unwrap_or(false), ve.is_some()) => (g, Some(base - g)),
            _ => (base, None),
        };
        let h = match ped_shift { Some(sh) if relief => f.subgrid_summit_m - (base + sh), _ => f.subgrid_summit_m - base };
        if h <= 0.0 { continue; }
        let (st, sb) = (shape.s_top_ppm as f64 / 1e6, shape.s_base_ppm as f64 / 1e6);
        let lambda = h / (st - sb);
        let rc = shape.crater_radius_mm as f64 / 1000.0;
        let r_main = rc + lambda * (st / sb).ln();
        let age = age_override_my.unwrap_or(recs.t_my - f.last_eruption_my).max(0.0);
        let e_vent = crate::t1b_ve::erosion_at(ve, t0, chart, table, v[0], v[1])?;
        let delta = vp.erosion_mult_permille as f64 / 1000.0 * e_vent * age;
        let dc = shape.crater_depth_mm as f64 / 1000.0;
        // Downhill direction of the T0 bilinear at the vent (breach).
        let gx = bil_at(t0, chart, v[0] + 500.0, v[1])? - bil_at(t0, chart, v[0] - 500.0, v[1])?;
        let gy = bil_at(t0, chart, v[0], v[1] + 500.0)? - bil_at(t0, chart, v[0], v[1] - 500.0)?;
        let breach_dir = if gx == 0.0 && gy == 0.0 { 0.0 } else { (-gy).atan2(-gx) };
        // Keyed angles.
        let tile = ((x * 1000.0 / SPACING_MM as f64).floor().div_euclid(128.0) as i64, (y * 1000.0 / SPACING_MM as f64).floor().div_euclid(128.0) as i64);
        let sc = layer.tile(8, tile.0, tile.1).map_err(|e| e.to_string())?;
        let hb = blake3::hash(f.id.as_bytes());
        let base_ctr = u32::from_le_bytes([hb.as_bytes()[0], hb.as_bytes()[1], hb.as_bytes()[2], hb.as_bytes()[3]]) as u64 * 256;
        let ng = vp.gully_count.max(1);
        let mut gul = Vec::new();
        for k in 0..ng {
            let u = sc.draw(PURPOSE_GULLY, base_ctr + k as u64, realization).map_err(|e| e.to_string())?.uniform_f64();
            gul.push((k as f64 + u) * std::f64::consts::TAU / ng as f64);
        }
        let mut dk = Vec::new();
        for k in 0..vp.dike_count {
            let u = sc.draw(PURPOSE_DIKE, base_ctr + k as u64, realization).map_err(|e| e.to_string())?.uniform_f64();
            dk.push(u * std::f64::consts::TAU);
        }
        let residual = delta >= h;
        let neck = vp.neck_relief_permille as f64 / 1000.0 * h;
        let mut lph = Vec::new();
        for k in 2..=6u64 { lph.push(sc.draw(PURPOSE_LOBE, base_ctr + k, realization).map_err(|e| e.to_string())?.uniform_f64() * std::f64::consts::TAU); }
        let lsum = |th: f64| -> f64 { lph.iter().enumerate().map(|(j, p)| ((j + 2) as f64 * th + p).cos() / (j + 2) as f64).sum() };
        let lnorm = (0..720).map(|k| lsum(k as f64 * std::f64::consts::TAU / 720.0).abs()).fold(1e-9, f64::max);
        // Influence: apron until it could meet the lowest plausible ground (2 km below the base), at most 3 R.
        let fr = vp.fade_ratio_permille as f64 / 1000.0;
        let r_inf = if residual { (vp.dike_length_permille_of_r as f64 / 1000.0 * r_main).max(vp.neck_radius_mm as f64 / 1000.0 + neck / (vp.wall_slope_ppm as f64 / 1e6)) + 100.0 } else { r_main * (1.0 + 2.0 * (fr - 1.0)) };
        // Lane CF C1/C4/C5: scaled cone with flows and gullies (constructional and dissected cones only).
        let mut lambda = lambda;
        let mut r_main = r_main;
        let mut r_inf = r_inf;
        let mut cfv: Option<CfVolcano> = None;
        if let Some(cf) = cf { if !residual && cf.cone_hw_permille > 0 {
            let q = cf.cone_hw_permille as f64 / 1000.0;
            let wb = (h / q).clamp(cf.cone_wb_min_mm as f64 / 1000.0, cf.cone_wb_max_mm as f64 / 1000.0);
            let mut rf = wb / 2.0;
            let lam = match solve_lambda(st, rf - rc, h) { Some(l) => l, None => { rf = rc + h / st; 1e9 } };
            let s_base = st * (-(rf - rc) / lam).exp();
            lambda = lam; r_main = rf;
            let r_end_max = rf * (1.0 + cf.foot_taper_permille as f64 / 1000.0 * 1.8);
            r_inf = r_end_max;
            let pick = |u: f64, a: [i64; 2], scale: f64| (a[0] as f64 + u * (a[1] - a[0]) as f64) * scale;
            let mut flows = Vec::new();
            for k in 0..cf.flow_count as u64 {
                let mut u = [0f64; 5];
                for (j, uj) in u.iter_mut().enumerate() { *uj = sc.draw(PURPOSE_FLOW, base_ctr + 16 + 5 * k + j as u64, realization).map_err(|e| e.to_string())?.uniform_f64(); }
                let mut wv = [0f64; 3];
                for (j, uj) in wv.iter_mut().enumerate() { *uj = sc.draw(PURPOSE_FLOW, base_ctr + 120 + 3 * k + j as u64, realization).map_err(|e| e.to_string())?.uniform_f64(); }
                let ang = (k as f64 + u[0]) * std::f64::consts::TAU / cf.flow_count as f64;
                let r0 = pick(u[1], cf.flow_r0_permille, rf / 1000.0);
                let len = pick(u[2], cf.flow_len_permille, rf / 1000.0).min((r_end_max - r0).max(0.0));
                flows.push([ang, r0, len, pick(u[3], cf.flow_halfwidth_mm, 1e-3), pick(u[4], cf.flow_thickness_mm, 1e-3), pick(wv[0], cf.wander_amp_mm, 1e-3), pick(wv[1], cf.wander_len_mm, 1e-3), wv[2] * std::f64::consts::TAU]);
            }
            let mut gullies = Vec::new();
            for k in 0..cf.gully_count as u64 {
                let mut u = [0f64; 3];
                for (j, uj) in u.iter_mut().enumerate() { *uj = sc.draw(PURPOSE_YGULLY, base_ctr + 64 + 3 * k + j as u64, realization).map_err(|e| e.to_string())?.uniform_f64(); }
                let mut wv = [0f64; 3];
                for (j, uj) in wv.iter_mut().enumerate() { *uj = sc.draw(PURPOSE_YGULLY, base_ctr + 200 + 3 * k + j as u64, realization).map_err(|e| e.to_string())?.uniform_f64(); }
                gullies.push([(k as f64 + u[0]) * std::f64::consts::TAU / cf.gully_count as f64, pick(u[1], cf.gully_r0_permille, rf / 1000.0), pick(u[2], cf.gully_depth_mm, 1e-3), pick(wv[0], cf.wander_amp_mm, 1e-3), pick(wv[1], cf.wander_len_mm, 1e-3), wv[2] * std::f64::consts::TAU]);
            }
            let vt = nearest_t1a(v).map(|b| b.2).unwrap_or(f64::INFINITY);
            cfv = Some(CfVolcano { wb_m: 2.0 * rf, r_foot_m: rf, q_built: h / (2.0 * rf), lambda_m: lam, s_base, r_end_max_m: r_end_max, clearance_m: clear, flows, gullies, vent_to_t1a_m: vt });
        } }
        let cell_of = |p: [f64; 2]| -> R<u32> { let (la, lo) = chart.inverse_f64(p[0] * 1000.0, p[1] * 1000.0).map_err(|e| e.to_string())?; Ok(t0.cell_of_deg(la * RAD2DEG, lo * RAD2DEG)) };
        if ve.is_some_and(|p| p.flank2.is_some()) && !residual && crate::t1b_ve::family_of(&f.kind)==Some(crate::t1b_ve::Family::Composite) {
            r_inf = r_inf.max(4.0*r_main);
        }
        out.push(VolcanoSite {
            id: f.id.clone(), kind: f.kind.clone(), shape: shape.name.clone(), record_xy_m: [x, y], vent_xy_m: v, moved_m: (v[0] - x).hypot(v[1] - y), record_cell: cell_of([x, y])?, vent_cell: cell_of(v)?,
            summit_record_m: f.subgrid_summit_m, base_m: base, h_m: h, lambda_m: lambda, r_main_m: r_main, r_influence_m: r_inf, s_top: st, s_base: sb,
            crater_radius_m: rc, crater_depth_m: dc, age_my: age, e_vent_m_per_my: e_vent, denudation_m: delta,
            stage: if residual { "residual (neck and dikes)".into() } else if delta > 0.0 { "dissected cone".into() } else { "constructional cone".into() },
            breached: !residual && delta >= dc, breach_dir_rad: breach_dir, gully_dirs_rad: gul, dike_dirs_rad: dk, neck_relief_m: if residual { neck } else { 0.0 },
            lobe_phase_rad: lph.clone(), lobe_norm: lnorm, cf: cfv, ve: None, ped_shift_m: ped_shift, built_summit_m: if relief && ped_shift.is_some() { Some(base + h) } else { None },
        });
    }
    Ok(out)
}

#[inline]
pub(crate) fn ang_diff(a: f64, b: f64) -> f64 { let d = (a - b).rem_euclid(std::f64::consts::TAU); if d > std::f64::consts::PI { d - std::f64::consts::TAU } else { d } }

/// The edifice at chart (x, y) metres.
pub fn edifice_at(vp: &VolcanoParamsV1, v: &VolcanoSite, x: f64, y: f64) -> EdificeAt {
    let (dx, dy) = (x - v.vent_xy_m[0], y - v.vent_xy_m[1]);
    let r = dx.hypot(dy);
    if r > v.r_influence_m { return EdificeAt::None; }
    let th = dy.atan2(dx);
    if v.denudation_m >= v.h_m {
        let rn = vp.neck_radius_mm as f64 / 1000.0;
        let wall = vp.wall_slope_ppm as f64 / 1e6;
        let mut dz = (v.neck_relief_m - wall * (r - rn).max(0.0)).max(0.0);
        let ld = vp.dike_length_permille_of_r as f64 / 1000.0 * v.r_main_m;
        let wd = vp.dike_width_mm as f64 / 1000.0;
        let hd = vp.dike_relief_permille as f64 / 1000.0 * v.neck_relief_m;
        for &a in &v.dike_dirs_rad {
            let d = ang_diff(th, a);
            if d.abs() >= std::f64::consts::FRAC_PI_2 { continue; }
            let (along, lat) = (r * d.cos(), r * d.sin().abs());
            if along > ld || lat > wd { continue; }
            dz = dz.max(hd * (1.0 - along / ld) * (1.0 - lat / wd));
        }
        return if dz > 0.0 { EdificeAt::Add(dz) } else { EdificeAt::None };
    }
    let scale = (v.h_m - v.denudation_m) / v.h_m;
    let top = v.summit_built() - v.denudation_m;
    let (rc, lam) = (v.crater_radius_m, v.lambda_m);
    let rho_r = v.r_main_m - rc;
    let mut z = if r < rc {
        top - v.crater_depth_m * scale * (1.0 - (r / rc).powi(2))
    } else {
        let rho = r - rc;
        let _ = rho_r;
        let drop = v.s_top * lam * (1.0 - (-rho / lam).exp());
        top - drop * scale
    };
    if v.ve.as_ref().is_some_and(|g| g.flank2.is_some()) {
        if r >= rc { z = v.base_m + crate::t1b_ve::flank2_height(v,r); }
        z += flank_lava(v,dx,dy);
        if v.breached {
            let d=ang_diff(th,v.breach_dir_rad);
            if d.abs()<std::f64::consts::FRAC_PI_2 && r*d.sin().abs()<(rc*0.5).max(60.0) {
                z=z.min(top-v.crater_depth_m*scale+0.05*r);
            }
        }
    } else if let Some(f) = v.ve.as_ref().and_then(|g| g.flank.as_ref()) {
        // The smooth envelope and crater above remain exact; only the flank correction is reconstructed.
        z += flank_lava(v, dx, dy) + f.at(dx, dy);
        // Keep the pre-existing record-age crater breach. Its floor only intersects the upper rim;
        // it is not one of the long flank wedges replaced by the new drainage law.
        if v.breached {
            let d = ang_diff(th, v.breach_dir_rad);
            let w = vp.gully_width_mm as f64 / 1000.0 + vp.gully_width_ppm_of_r as f64 / 1e6 * r;
            if d.abs() < std::f64::consts::FRAC_PI_2 && r * d.sin().abs() < w.max(rc * 0.5) {
                z = z.min(top - v.crater_depth_m * scale + 0.05 * r);
            }
        }
    } else if v.denudation_m > 0.0 {
        let above = (z - (v.base_m)).max(0.0);
        let depth0 = v.denudation_m.min(vp.gully_depth_permille as f64 / 1000.0 * above);
        let fade = ((r - rc) / (0.3 * v.r_main_m)).clamp(0.0, 1.0) * ((v.r_influence_m - r) / (0.3 * v.r_main_m)).clamp(0.0, 1.0);
        let w = vp.gully_width_mm as f64 / 1000.0 + vp.gully_width_ppm_of_r as f64 / 1e6 * r;
        let mut cut: f64 = 0.0;
        for &a in &v.gully_dirs_rad {
            let d = ang_diff(th, a);
            if d.abs() >= std::f64::consts::FRAC_PI_2 { continue; }
            let lat = r * d.sin().abs();
            if lat < w { cut = cut.max(depth0 * fade * (1.0 - lat / w)); }
        }
        z -= cut;
        if v.breached {
            let d = ang_diff(th, v.breach_dir_rad);
            let lat = r * d.sin().abs();
            if d.abs() < std::f64::consts::FRAC_PI_2 && lat < w.max(rc * 0.5) {
                let floor = top - v.crater_depth_m * scale + 0.05 * r;
                z = z.min(floor);
            }
        }
    }
    let fr = vp.fade_ratio_permille as f64 / 1000.0;
    let lob = if v.lobe_phase_rad.is_empty() { 0.0 } else { v.lobe_phase_rad.iter().enumerate().map(|(j, p)| ((j + 2) as f64 * th + p).cos() / (j + 2) as f64).sum::<f64>() / v.lobe_norm };
    let r_end = v.r_main_m * (1.0 + (fr - 1.0) * (1.0 + 0.8 * lob));
    let u = if r_end <= v.r_main_m + 1.0 { if r <= v.r_main_m { 1.0 } else { 0.0 } } else { ((r_end - r) / (r_end - v.r_main_m)).clamp(0.0, 1.0) };
    EdificeAt::Abs(z, u * u * (3.0 - 2.0 * u))
}

// ------------------------------------------------------------------ fans (mouth deposition)

#[derive(Clone, Debug, Serialize)]
pub struct Fan { pub x_m: f64, pub y_m: f64, pub area_m2: f64, pub e_m_per_my: f64, pub depth_m: f64, pub volume_m3: f64, pub radius_m: f64, pub curve_code: u8,
    /// B5: apex height (sea level for a bayhead fan; the channel bed for a fan delta at a perched mouth).
    pub apex_m: f64,
    /// B6: outflow direction (rad) at the mouth.
    pub dir_rad: f64 }

/// Bayheads and their fans, chart-wide, from the channel curves.
pub fn fans_for(cp: &CoastParamsV1, curves: &[Curve], t0: &T0, chart: &Chart, table: &BTreeMap<u32, f64>, flood: &BTreeMap<u32, f64>, cf: Option<&CfParamsV1>) -> R<Vec<Fan>> {
    let sea = t0.sea_level_m;
    let mut seen: BTreeMap<(i64, i64), usize> = BTreeMap::new();
    let mut out: Vec<Fan> = Vec::new();
    let half = vz_tilekey::CHART_HALF_EXTENT_MM as f64 / 1000.0;
    let shelf = cp.shelf_m();
    for c in curves {
        // Lane CF R1c: a curve that starts at or below sea level is downstream of a mouth (no second fan).
        if cf.is_some() && c.pts.first().map(|p| p[2] <= sea).unwrap_or(true) { continue; }
        for w in c.pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            // B5: the mouth is the first of a bed crossing sea level and an entry into construction water.
            if b[0].abs() > half || b[1].abs() > half || a[0].abs() > half || a[1].abs() > half { break; }
            let a_wet = bil_at(t0, chart, a[0], a[1])? <= sea - shelf;
            if a_wet { break; }
            let b_wet = bil_at(t0, chart, b[0], b[1])? <= sea - shelf;
            let crosses = a[2] > sea && b[2] <= sea;
            if !crosses && !b_wet { continue; }
            let (x, y, area, apex) = if crosses {
                let u = (a[2] - sea) / (a[2] - b[2]);
                (a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1]), a[3] + u * (b[3] - a[3]), sea)
            } else { (a[0], a[1], a[3], a[2].max(sea)) };
            if x.abs() > half || y.abs() > half || area < cp.fan_area_min_m2 as f64 { break; }
            let key = ((x).round() as i64, (y).round() as i64);
            if seen.contains_key(&key) { break; }
            // Lane CF R1c: a fan whose apex T0 depth exceeds fan_subaqueous x shelf is subaqueous (not built).
            if let Some(cf) = cf { if sea - bil_at(t0, chart, x, y)? > cf.fan_subaqueous_permille as f64 / 1000.0 * shelf { break; } }
            let e = table_at(t0, chart, table, x, y)?;
            let h = table_at(t0, chart, flood, x, y)?.max(cp.fan_depth_min_mm as f64 / 1000.0);
            let vol = e * 1e-6 * area * cp.fan_time_yr as f64;
            let r = (2.0 * vol / (std::f64::consts::PI * h)).sqrt().min(cp.fan_radius_max_mm as f64 / 1000.0);
            seen.insert(key, out.len());
            let dir = (y - a[1]).atan2(x - a[0]);
            let dir = if (y - a[1]).hypot(x - a[0]) < 1e-6 { (b[1] - a[1]).atan2(b[0] - a[0]) } else { dir };
            out.push(Fan { x_m: x, y_m: y, area_m2: area, e_m_per_my: e, depth_m: h, volume_m3: vol, radius_m: r, curve_code: c.code, apex_m: apex, dir_rad: dir });
            break;
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------ lattice helpers

/// Proper crossing of segments ab and cd.
fn seg_cross(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let o = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
    let (d1, d2, d3, d4) = (o(a, b, c), o(a, b, d), o(c, d, a), o(c, d, b));
    (d1 > 0.0) != (d2 > 0.0) && (d3 > 0.0) != (d4 > 0.0)
}

/// Lane RV's t1a_bank_v1, re-applied after the coast steps (B2): a free post 8-adjacent to a T1a channel post stands at
/// least 1 cm above the highest adjacent T1a bed.
pub fn bank_rule(lat: &EvalLat, z: &mut [f64], fixed: &[u8], guides: bool) {
    let n = lat.n;
    for i in 0..lat.nodes() {
        // Lane CF fix B1 (guides = true): network channel posts next to a T1a channel obey the bank rule too (a fan-delta
        // network channel below the T1a chute at the mouth otherwise captures the T1a path; F1 strict).
        if fixed[i] != FREE && !(guides && fixed[i] == GUIDE) { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let mut hb = f64::NEG_INFINITY;
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if fixed[j] == CHANNEL { hb = hb.max(z[j]); }
        }
        if hb.is_finite() && z[i] < hb + 0.01 { z[i] = hb + 0.01; }
    }
}

/// 8-neighbour chamfer distance (metres) from `src` posts, two passes (exact for the chamfer metric without obstacles).
pub fn chamfer(lat: &EvalLat, src: &dyn Fn(usize) -> bool, cap: f64) -> Vec<f64> {
    let (n, nn, dx) = (lat.n, lat.nodes(), lat.dx());
    let dd = dx * std::f64::consts::SQRT_2;
    let mut d: Vec<f64> = (0..nn).map(|i| if src(i) { 0.0 } else { f64::INFINITY }).collect();
    for y in 0..n { for x in 0..n {
        let i = y * n + x;
        let mut v = d[i];
        if x > 0 { v = v.min(d[i - 1] + dx); }
        if y > 0 { v = v.min(d[i - n] + dx); if x > 0 { v = v.min(d[i - n - 1] + dd); } if x + 1 < n { v = v.min(d[i - n + 1] + dd); } }
        d[i] = v;
    } }
    for y in (0..n).rev() { for x in (0..n).rev() {
        let i = y * n + x;
        let mut v = d[i];
        if x + 1 < n { v = v.min(d[i + 1] + dx); }
        if y + 1 < n { v = v.min(d[i + n] + dx); if x + 1 < n { v = v.min(d[i + n + 1] + dd); } if x > 0 { v = v.min(d[i + n - 1] + dd); } }
        d[i] = v;
    } }
    for v in d.iter_mut() { if *v > cap { *v = f64::INFINITY; } }
    d
}

/// Everything evaluate() needs beyond 0.9.0 (empty = 0.9.0 behaviour).
#[derive(Clone, Default)]
pub struct CoastCtx {
    pub flood: BTreeMap<u32, f64>,
    pub fans: Vec<Fan>,
    pub volcanoes: Vec<VolcanoSite>,
    /// REVISION 1: T1a channel curve segments (x0, y0, x1, y1, chart m) near the volcanoes (corridor cap distance).
    pub t1a_segs: Vec<[f64; 4]>,
    /// Lane CF (0.9.2) params (None = 0.9.1 behaviour).
    pub cf: Option<CfParamsV1>,
    /// Lane TD (0.9.3) params (None = 0.9.2 behaviour).
    pub fx: Option<crate::t1b_fix::FxParamsV1>,
    /// Lane TF params (None = lane TE behaviour) and, per volcano, the C1 T0 surface at the vent (T5).
    pub tf: Option<crate::t1b_tf::TfParamsV1>,
    pub tf_vent_c1: Vec<f64>,
    /// Lane TG params (None = lane TF behaviour) and the declared valley-floor slope (params.floor_slope_ppm, G2).
    pub tg: Option<crate::t1b_tg::TgParamsV1>,
    pub tg_floor_slope: f64,
    /// Lane TH params (None = lane TG behaviour).
    pub th: Option<crate::t1b_th::ThParamsV1>,
    /// Lane VE params (None = lane TH behaviour) and whether this evaluation builds the units (true only for the published
    /// window surface; the layer's E / F solve evaluations see the family body without the units).
    pub ve: Option<crate::t1b_ve::VeParamsV1>,
    pub ve_units: bool,
    /// Lane STEP2 (ct var >= 1): the keyed step-thickness closure boundaries (ascending, m; empty otherwise).
    pub ct_steps: Vec<f64>,
    /// Lane BENCH-VAR (ct var 3): the per-bed keyed tables (empty otherwise).
    pub ct_beds: crate::t1b_ct::BedTables,
    /// Lane GLACIAL-T1 (glacial key only): the chart glacier table evaluate() carves from (None = no carve).
    pub gl: Option<std::sync::Arc<crate::t1b_glacial::GlacialTable>>,
    /// Lane EDIFICE-FOOT (window-only key ef_cap=1): the composite foot drape cap (None = unchanged).
    pub ef: Option<crate::t1b_ef::EdificeFootParamsV1>,
}

/// T1a channel segments within reach of any edifice (reach = influence + max relief / cap slope).
pub fn t1a_segments_near(vp: &VolcanoParamsV1, vs: &[VolcanoSite], curves: &[Curve]) -> Vec<[f64; 4]> {
    let cs = vp.cap_slope_ppm as f64 / 1e6;
    let mut out = Vec::new();
    for c in curves.iter().filter(|c| c.code == CHANNEL) {
        for w in c.pts.windows(2) {
            let s = [w[0][0], w[0][1], w[1][0], w[1][1]];
            let near = vs.iter().any(|v| {
                let reach = v.r_influence_m.max(v.ve.as_ref().map(|g| g.r_units_m).unwrap_or(0.0)) + v.h_m / cs + 1000.0;
                let q = seg_nearest(v.vent_xy_m, [s[0], s[1]], [s[2], s[3]]);
                (q[0] - v.vent_xy_m[0]).hypot(q[1] - v.vent_xy_m[1]) <= reach
            });
            if near { out.push(s); }
        }
    }
    out
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct CoastStats {
    #[serde(skip)]
    pub surface_capture: Option<crate::t1b_surface_stages::Capture>,
    /// Arc2 process increment at this native evaluation stage, before later coast/fill constraints.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flank2_projection: Option<serde_json::Value>,
    pub flooded_posts: u32,
    pub fan_raised_posts: u32,
    pub fan_land_gained_posts: u32,
    pub cliff_cut_posts: u32,
    pub beach_posts: u32,
    pub edifice_posts: u32,
    pub buried_channel_posts: u32,
    pub final_fill_raised_posts: u32,
    pub fans_in_lattice: u32,
    pub volcanoes_in_lattice: u32,
    /// Lane TD5 (only with fx): posts lowered by more than 0.5 m by the soft-rock bluff.
    #[serde(skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub bluff_cut_posts: u32,
    /// Lane TF (only with tf).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tf: Option<crate::t1b_tf::TfCoastStats>,
    /// Lane TG (only with tg).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tg: Option<crate::t1b_tg::TgCoastStats>,
    /// Lane TG G1 (in-memory only, never serialized): per post, the cut profile extended seaward on posts the wave cut sends
    /// to sea (NaN elsewhere); allocated by the caller when the shore field is on.
    #[serde(skip)]
    pub tg_cut_ext: Vec<f64>,
    /// Lane TH (only with th).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub th: Option<crate::t1b_th::ThCoastStats>,
    /// Lane TH H2 (in-memory only, never serialized): per post, the height before the wave cut on cliff-class posts (NaN
    /// elsewhere); allocated by the caller when uncut is on.
    #[serde(skip)]
    pub th_uncut: Vec<f64>,
    /// Lane VE (only with ve).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ve: Option<crate::t1b_ve::VeCoastStats>,
    /// Lane VE (in-memory only): edifice_style / edifice_unit per post; allocated by the caller when ve is on.
    #[serde(skip)]
    pub ve_style: Vec<u8>,
    #[serde(skip)]
    pub ve_unit: Vec<u16>,
    /// Lane VE4 (in-memory only): front_mm distances in m (NaN none); allocated by the caller when ve4 front is on.
    #[serde(skip)]
    pub ve_front: Vec<f64>,
    /// Lane SEAM-CERT (in-memory only): the fill taint while seam_cert is on (the fills are then seeded only from fixed posts).
    #[serde(skip)]
    pub seam: Option<crate::t1b_seam::SeamTaint>,
    /// Lane DAMMED-VALLEY (in-memory only): heights just before the coast / volcano final fill; allocated (len nn) by the
    /// caller when dam_lake is on, empty otherwise.
    #[serde(skip)]
    pub dam_pre: Vec<f64>,
    /// Lane EDIFICE-FOOT (only with ef_cap).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ef: Option<crate::t1b_ef::EfCoastStats>,
    /// Lane EDIFICE-FOOT (in-memory only): posts where the foot drape cap bound (the body lies as a drape on the
    /// pre-edifice ground); len nn with ef_cap, empty otherwise.
    #[serde(skip)]
    pub ef_draped: Vec<bool>,
}

/// The coast and volcano steps on a constructed surface (in place). `fixed`: FREE, CHANNEL, MARINE, LAKE, GUIDE;
/// `z`: constructed heights (construction water at its base level). Returns the over-code per post.
#[allow(clippy::too_many_arguments)]
pub fn coast_and_volcano(lat: &EvalLat, sea: f64, cp: Option<&CoastParamsV1>, vp: Option<&VolcanoParamsV1>, cx: &CoastCtx, rock_k: &dyn Fn(usize, f64) -> f64, t1a_post: &[bool], bil: &[f32], z: &mut [f64], fixed: &mut [u8], st: &mut CoastStats) -> Vec<u8> {
    coast_and_volcano_td(lat, sea, cp, vp, cx, rock_k, None, t1a_post, bil, z, fixed, st).0
}

/// As `coast_and_volcano`, plus (lane TD) the rock Sc for the TD5 soft-rock bluff, and the volcano index per post (index into
/// `cx.volcanoes` of the cone that built or buried an edifice post; 255 elsewhere).
#[allow(clippy::too_many_arguments)]
pub fn coast_and_volcano_td(lat: &EvalLat, sea: f64, cp: Option<&CoastParamsV1>, vp: Option<&VolcanoParamsV1>, cx: &CoastCtx, rock_k: &dyn Fn(usize, f64) -> f64, rock_sc: Option<&dyn Fn(usize, f64) -> f64>, t1a_post: &[bool], bil: &[f32], z: &mut [f64], fixed: &mut [u8], st: &mut CoastStats) -> (Vec<u8>, Vec<u8>) {
    coast_and_volcano_tf(lat, sea, cp, vp, cx, rock_k, rock_sc, t1a_post, bil, z, fixed, st, None)
}

/// Read-only native construction stages; elevations are not automatically
/// sediment inventories, especially where a field represents a water plane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoastMaterialStage { BeforeFans, AfterFans, AfterWaveCuts, AfterBankConstraint, AfterFinalFill }

/// As `coast_and_volcano_td`, plus (lane TF, only when `cx.tf` and `tfin` are given) the C1 cone base, the sub-post shore,
/// fetch exposure, soft_bluff_v2 and the final-fill report (`crate::t1b_tf`).
#[allow(clippy::too_many_arguments)]
pub fn coast_and_volcano_tf(lat: &EvalLat, sea: f64, cp: Option<&CoastParamsV1>, vp: Option<&VolcanoParamsV1>, cx: &CoastCtx, rock_k: &dyn Fn(usize, f64) -> f64, rock_sc: Option<&dyn Fn(usize, f64) -> f64>, t1a_post: &[bool], bil: &[f32], z: &mut [f64], fixed: &mut [u8], st: &mut CoastStats, tfin: Option<&crate::t1b_tf::TfIn>) -> (Vec<u8>, Vec<u8>) {
    coast_and_volcano_tf_observed(lat, sea, cp, vp, cx, rock_k, rock_sc, t1a_post, bil, z, fixed, st, tfin, None)
}

/// Optional immutable stage snapshots; both existing wrappers keep observer=None.
/// Fan, wave-cut, bank-constraint and final-fill changes remain distinct.
#[allow(clippy::too_many_arguments)]
pub fn coast_and_volcano_tf_observed(lat: &EvalLat, sea: f64, cp: Option<&CoastParamsV1>, vp: Option<&VolcanoParamsV1>, cx: &CoastCtx, rock_k: &dyn Fn(usize, f64) -> f64, rock_sc: Option<&dyn Fn(usize, f64) -> f64>, t1a_post: &[bool], bil: &[f32], z: &mut [f64], fixed: &mut [u8], st: &mut CoastStats, tfin: Option<&crate::t1b_tf::TfIn>, observer: Option<&mut dyn FnMut(CoastMaterialStage, &[f64], &[u8], &[u8])>) -> (Vec<u8>, Vec<u8>) {
    coast_and_volcano_tf_owned(lat, sea, cp, vp, cx, rock_k, rock_sc, t1a_post, bil, z, fixed, st, tfin, observer, None, None).expect("legacy coast without owner cannot yield a stage error")
}

pub fn coast_and_volcano_tf_owned(lat: &EvalLat, sea: f64, cp: Option<&CoastParamsV1>, vp: Option<&VolcanoParamsV1>, cx: &CoastCtx, rock_k: &dyn Fn(usize, f64) -> f64, rock_sc: Option<&dyn Fn(usize, f64) -> f64>, t1a_post: &[bool], bil: &[f32], z: &mut [f64], fixed: &mut [u8], st: &mut CoastStats, tfin: Option<&crate::t1b_tf::TfIn>, mut observer: Option<&mut dyn FnMut(CoastMaterialStage, &[f64], &[u8], &[u8])>, owner: Option<&dyn crate::t1b_surface_stages::Conditioner>, stop: Option<crate::t1b_surface_stages::Stage>) -> R<(Vec<u8>, Vec<u8>)> {
    let esc_trace = crate::t1b_esc_deps::tracing();
    let (n, nn, dx) = (lat.n, lat.nodes(), lat.dx());
    let tfp = if tfin.is_some() { cx.tf.as_ref() } else { None };
    let tf_on = |f: fn(&crate::t1b_tf::TfParamsV1) -> u32| tfp.map(|t| f(t) == 1).unwrap_or(false);
    let mut over = vec![OVER_NONE; nn];
    let mut vid = vec![255u8; nn];
    let pos = |i: usize| { let (gx, gy) = lat.g(i); (gx as f64 * dx, gy as f64 * dx) };
    let (lx0, ly0) = pos(0);
    let (lx1, ly1) = pos(nn - 1);
    // Volcano edifices.
    let mut ve_zpre: Vec<f64> = Vec::new();
    if let Some(vp) = vp {
        // Lane VE (REVISION 1): the family body and the units on the published window surface only; the layer's E / F solve
        // evaluations keep lane CF's cone (so the solved tables equal the ve-off ones).
        let ve_on = cx.ve.is_some() && cx.ve_units;
        let rp_reach = |v: &VolcanoSite| -> f64 { match cx.cf.as_ref().and_then(|c| c.ring_plain.as_ref()) { Some(rp) => ring_plain_reach(rp, vp, v), None => 0.0 } };
        let vsi: Vec<usize> = (0..cx.volcanoes.len()).filter(|&k| { let v = &cx.volcanoes[k];
            let ri = if ve_on { v.r_influence_m.max(v.ve.as_ref().map(|g| g.r_units_m).unwrap_or(0.0)) } else { v.r_influence_m };
            let ri = ri.max(rp_reach(v));
            v.vent_xy_m[0] + ri >= lx0 && v.vent_xy_m[0] - ri <= lx1 && v.vent_xy_m[1] + ri >= ly0 && v.vent_xy_m[1] - ri <= ly1 }).collect();
        if ve_on && st.ve_style.len() == nn { ve_zpre = z.to_vec(); }
        let vs: Vec<&VolcanoSite> = vsi.iter().map(|&k| &cx.volcanoes[k]).collect();
        st.volcanoes_in_lattice = vs.len() as u32;
        if cx.ef.is_some() { st.ef = Some(Default::default()); st.ef_draped = vec![false; nn]; }
        if !vs.is_empty() {
            let (mi, cs) = (vp.mask_inner_mm as f64 / 1000.0, vp.cap_slope_ppm as f64 / 1e6);
            let segs = &cx.t1a_segs;
            for i in 0..nn {
                if fixed[i] == CHANNEL || fixed[i] == LAKE || t1a_post[i] || !z[i].is_finite() { continue; }
                let (x, y) = pos(i);
                let mut rel: f64 = 0.0;
                // Lane CF C2: scaled cones blend (and may bury ridges: negative relief).
                let mut low: f64 = 0.0;
                let mut on_cone = false;
                // Lane TD6: the volcano giving the largest raise (else the deepest burial, else the cone cover).
                let (mut k_rel, mut k_low, mut k_on) = (255u8, 255u8, 255u8);
                let mut ef_cut: f64 = 0.0;
                for (vk, v) in vs.iter().enumerate() {
                    let kid = vsi[vk].min(254) as u8;
                    // Lane VE: the family body (composite = the CF cone without its lobes; shield, rift field; plateau none).
                    if let (true, Some(cf), Some(_)) = (ve_on, cx.cf.as_ref(), v.ve.as_ref()) {
                        let bxy = match tfin { Some(t) if tf_on(|p| p.cone_c1) && !t.bil_c1.is_empty() && cx.tf_vent_c1.len() == cx.volcanoes.len() => t.bil_c1[i] as f64 - cx.tf_vent_c1[vsi[vk]] + v.base_m,
                            _ => bil.get(i).map(|b| *b as f64 - v.ped_shift()).unwrap_or(v.base_m) };
                        if let Some((e, w, u)) = crate::t1b_ve::body(vp, v, cf, x, y, bxy) {
                            let g = z[i];
                            let tc = crate::t1b_ve::thickness_cap(v, x, y);
                            let e = e.min(g + tc);
                            // Lane EDIFICE-FOOT (ef_cap): composites (VE's cap is infinite) take the foot drape cap.
                            let e = match cx.ef.as_ref() { Some(ef) if tc.is_infinite() => { let ec = e.min(g + ef.cap(vp, v, x, y)); ef_cut = ef_cut.max(e - ec); ec }, _ => e };
                            let arc2=v.ve.as_ref().is_some_and(|g|g.flank2.is_some());
                            let gp = g - w * (g - e).max(0.0);
                            let s = cf.softmax(e, gp);
                            let d = if arc2 {
                                let h=crate::t1b_ve::flank2_height(v,(x-v.vent_xy_m[0]).hypot(y-v.vent_xy_m[1]));
                                let bury=1.0-(-h/150.0).exp();
                                bury*(e-g)+(1.0-bury)*h
                            } else { (s-g)*u };
                            if d > 0.0 { if d > rel { k_rel = kid; } rel = rel.max(d); } else if fixed[i] == FREE || fixed[i] == GUIDE { if d < low { k_low = kid; } low = low.min(d); }
                            if u > 0.5 && w > 0.0 && e >= gp - 1.0 { on_cone = true; if k_on == 255 { k_on = kid; } }
                        }
                        // Lane VH-I ring_plain_tail_v1 (window-only key).
                        let dt = ring_plain_raise(cf, vp, v, x, y, bxy, z[i], sea, fixed[i]);
                        if dt > 0.0 { if dt > rel { k_rel = kid; } rel = rel.max(dt); }
                        continue;
                    }
                    if let (Some(cf), Some(_)) = (cx.cf.as_ref(), v.cf.as_ref()) {
                        // Lane TF T5: the C1 T0 surface relative to its value at the vent (t = c1(post) - c1(vent)).
                        let bxy = match tfin { Some(t) if tf_on(|p| p.cone_c1) && !t.bil_c1.is_empty() && cx.tf_vent_c1.len() == cx.volcanoes.len() => t.bil_c1[i] as f64 - cx.tf_vent_c1[vsi[vk]] + v.base_m,
                            _ => bil.get(i).map(|b| *b as f64 - v.ped_shift()).unwrap_or(v.base_m) };
                        if let Some((e, w, u)) = cone_cf(vp, v, cf, x, y, bxy) {
                            let g = z[i];
                            // Lane EDIFICE-FOOT (ef_cap): the foot drape cap.
                            let e = match cx.ef.as_ref() { Some(ef) => { let ec = e.min(g + ef.cap(vp, v, x, y)); ef_cut = ef_cut.max(e - ec); ec }, None => e };
                            let gp = g - w * (g - e).max(0.0);
                            let s = cf.softmax(e, gp);
                            let d = if v.ve.as_ref().is_some_and(|g|g.flank2.is_some()) {
                                let h=crate::t1b_ve::flank2_height(v,(x-v.vent_xy_m[0]).hypot(y-v.vent_xy_m[1]));
                                let bury=1.0-(-h/150.0).exp();bury*(e-g)+(1.0-bury)*h
                            } else { (s-g)*u };
                            if d > 0.0 { if d > rel { k_rel = kid; } rel = rel.max(d); } else if fixed[i] == FREE || fixed[i] == GUIDE { if d < low { k_low = kid; } low = low.min(d); }
                            // Fix B4: volcanic cover only where the cone buries (w > 0) or raises the surface (below).
                            if u > 0.5 && w > 0.0 && e >= gp - 1.0 { on_cone = true; if k_on == 255 { k_on = kid; } }
                        }
                        // Lane VH-I ring_plain_tail_v1 (window-only key).
                        let dt = ring_plain_raise(cf, vp, v, x, y, bxy, z[i], sea, fixed[i]);
                        if dt > 0.0 { if dt > rel { k_rel = kid; } rel = rel.max(dt); }
                        continue;
                    }
                    match edifice_at(vp, v, x, y) {
                        EdificeAt::Abs(ze, fade) => if ze > z[i] { let d = (ze - z[i]) * fade; if d > rel { k_rel = kid; } rel = rel.max(d); },
                        EdificeAt::Add(dz) => { if dz > rel { k_rel = kid; } rel = rel.max(dz) },
                        EdificeAt::None => {}
                    }
                }
                if ef_cut > 1e-3 { if st.ef_draped.len() == nn { st.ef_draped[i] = true; } let s = st.ef.get_or_insert_with(Default::default); s.capped_posts += 1; s.cut_max_m = s.cut_max_m.max(ef_cut); s.cut_volume_m3 += ef_cut * dx * dx; }
                vid[i] = if k_rel != 255 && rel > 0.0 { k_rel } else if k_low != 255 && low < 0.0 { k_low } else { k_on };
                if rel <= 0.0 && low < 0.0 {
                    if low < -0.5 || on_cone { over[i] = OVER_EDIFICE; st.edifice_posts += 1; }
                    if fixed[i] == GUIDE { fixed[i] = FREE; st.buried_channel_posts += 1; }
                    z[i] += low;
                    continue;
                }
                if rel <= 0.0 { if on_cone { over[i] = OVER_EDIFICE; st.edifice_posts += 1; } continue; }
                // Corridor cap: relief <= cs x (d - mi), d = distance to the T1a curves (only needed when it could bind).
                let dneed = mi + rel / cs;
                let mut d2 = dneed * dneed;
                if esc_trace { crate::t1b_esc_deps::cap_begin(0,x,y,dneed); }
                for s in segs {
                    if (s[0].min(s[2]) - x) > dneed || (x - s[0].max(s[2])) > dneed || (s[1].min(s[3]) - y) > dneed || (y - s[1].max(s[3])) > dneed { continue; }
                    if esc_trace { crate::t1b_esc_deps::cap_segment(s); }
                    let q = seg_nearest([x, y], [s[0], s[1]], [s[2], s[3]]);
                    d2 = d2.min((q[0] - x).powi(2) + (q[1] - y).powi(2));
                }
                if esc_trace { crate::t1b_esc_deps::geometry_end(); }
                let mut rel = rel.min(cs * (d2.sqrt() - mi).max(0.0));
                // B3 (exploratory): shadow behind a T1a channel as seen from the vent (of the volcano giving the relief).
                if rel > 0.0 {
                    for v in &vs {
                        let (a, b) = (v.vent_xy_m, [x, y]);
                        if (a[0] - x).hypot(a[1] - y) > v.r_influence_m.max(rp_reach(v)) { continue; }
                        let crossed = if esc_trace {
                            crate::t1b_esc_deps::shadow_begin(a,b);
                            let hit = segs.iter().any(|s| { crate::t1b_esc_deps::shadow_segment(a,b,s); seg_cross(a,b,[s[0],s[1]],[s[2],s[3]]) });
                            crate::t1b_esc_deps::geometry_end(); hit
                        } else { segs.iter().any(|s| seg_cross(a,b,[s[0],s[1]],[s[2],s[3]])) };
                        if crossed { rel = 0.0; break; }
                    }
                }
                let znew = z[i] + rel;
                if znew > z[i] + 1e-6 {
                    if znew > z[i] + 0.5 { over[i] = OVER_EDIFICE; st.edifice_posts += 1; }
                    if fixed[i] == GUIDE && znew > z[i] + 1.0 { fixed[i] = FREE; st.buried_channel_posts += 1; }
                    if fixed[i] == MARINE && znew > sea { fixed[i] = FREE; }
                    z[i] = znew;
                }
            }
        }
    }
    // Lane VE: the units of every volcano in the lattice, oldest first, on the body surface.
    let mut vest = crate::t1b_ve::VeCoastStats::default();
    let ve_lat = cx.ve.is_some() && st.ve_style.len() == nn && st.ve_unit.len() == nn && ve_zpre.len() == nn;
    if ve_lat {
        if let Some(vp) = vp {
            let vsi: Vec<usize> = (0..cx.volcanoes.len()).filter(|&k| { let v = &cx.volcanoes[k];
                let ri = v.r_influence_m.max(v.ve.as_ref().map(|g| g.r_units_m).unwrap_or(0.0));
                v.vent_xy_m[0] + ri >= lx0 && v.vent_xy_m[0] - ri <= lx1 && v.vent_xy_m[1] + ri >= ly0 && v.vent_xy_m[1] - ri <= ly1 }).collect();
            if !vsi.is_empty() {
                let mut sty = std::mem::take(&mut st.ve_style);
                let mut un = std::mem::take(&mut st.ve_unit);
                let mut fr = std::mem::take(&mut st.ve_front);
                crate::t1b_ve::apply_units_f(lat, vp, cx, &vsi, &ve_zpre, z, fixed, t1a_post, &mut over, &mut vid, &mut sty, &mut un, sea, &mut vest, &mut fr);
                st.ve_front = fr;
                st.ve_style = sty;
                st.ve_unit = un;
            }
        }
    }
    // The regional E/F solve sees the body, not the native eruption units.
    // The process reference contains those units. Apply its increment at their
    // window stage, retaining the regional construction fit as initial condition.
    // Literal whole-cell F6 still measures the resulting published surface.
    if let Some(f) = cx.ve.as_ref().and_then(|v| v.flank2.as_ref()).filter(|f| f.phase == 1 && cx.ve_units) {
        let (mut applied, mut withheld, mut covered) = (0.0, 0.0, 0usize);
        for i in 0..nn {
            let (gx, gy) = lat.g(i);
            if let Some(d) = f.at(gx, gy, dx) {
                covered += 1;
                if fixed[i] != FREE || t1a_post[i] {
                    withheld += d * dx * dx;
                } else {
                    z[i] += d;
                    applied += d * dx * dx;
                }
            }
        }
        st.flank2_projection = Some(serde_json::json!({"spacing_m":dx,"reference_covered_posts":covered,"applied_net_volume_m3":applied,"fixed_withheld_net_volume_m3":withheld,"stage":"before native coast/bank/fill constraints; not a global sediment ledger"}));
    }
    let Some(cp) = cp else {
        bank_rule(lat, z, fixed, cx.cf.is_some());
        if let Some(observe) = observer.as_deref_mut() { observe(CoastMaterialStage::AfterBankConstraint, z, fixed, &over); }
        // Volcano only: a final fill (an edifice can dam a valley).
        let outlet: Vec<bool> = (0..nn).map(|i| fixed[i] != FREE).collect();
        let zb: Vec<f64> = if ve_lat { z.to_vec() } else { Vec::new() };
        if st.dam_pre.len() == nn { st.dam_pre.copy_from_slice(z); }
        if stop == Some(crate::t1b_surface_stages::Stage::CoastFinal) { st.surface_capture=Some(crate::t1b_surface_stages::capture(z, fixed, st.seam.as_ref())); return Ok((over,vid)); }
        let (r, _) = if let Some(owner) = owner { crate::t1b_surface_stages::apply(owner, crate::t1b_surface_stages::Stage::CoastFinal, lat, z, fixed, st.seam.as_mut())? } else if let Some(sm) = st.seam.as_mut() { crate::t1b_seam::pinned_fill_native(lat, z, fixed, sm, "volcano_final") } else { crate::t1b_rivernet::fill(lat, z, &outlet) };
        st.final_fill_raised_posts = r;
        if ve_lat { crate::t1b_ve::lake_stats(&zb, z, &st.ve_style, &mut vest); st.ve = Some(vest); }
        for i in 0..nn { if over[i] != OVER_EDIFICE { vid[i] = 255; } }
        if let Some(observe) = observer.as_deref_mut() { observe(CoastMaterialStage::AfterFinalFill, z, fixed, &over); }
        return Ok((over, vid));
    };
    // Flooding at present sea level.
    let flood = |z: &[f64], fixed: &mut [u8], st: &mut CoastStats| {
        for i in 0..nn { if fixed[i] != LAKE && fixed[i] != MARINE && z[i] <= sea { fixed[i] = MARINE; st.flooded_posts += 1; } }
    };
    flood(z, fixed, st);
    if let Some(observe) = observer.as_deref_mut() { observe(CoastMaterialStage::BeforeFans, z, fixed, &over); }
    // Fans.
    let fs: Vec<&Fan> = cx.fans.iter().filter(|f| f.x_m + f.radius_m >= lx0 && f.x_m - f.radius_m <= lx1 && f.y_m + f.radius_m >= ly0 && f.y_m - f.radius_m <= ly1).collect();
    st.fans_in_lattice = fs.len() as u32;
    let sf = cp.fan_slope_ppm as f64 / 1e6;
    for f in &fs {
        let (i0, i1) = (((f.x_m - f.radius_m) / dx).floor() as i64, ((f.x_m + f.radius_m) / dx).ceil() as i64);
        let (j0, j1) = (((f.y_m - f.radius_m) / dx).floor() as i64, ((f.y_m + f.radius_m) / dx).ceil() as i64);
        for gy in j0..=j1 { for gx in i0..=i1 {
            let Some(i) = lat.idx(gx, gy) else { continue };
            if fixed[i] == CHANNEL || fixed[i] == LAKE || t1a_post[i] || !z[i].is_finite() { continue; }
            // Lane CF fix B2: R1c's subaqueous limit per post (no emergent fan over T0 water deeper than fan_subaqueous x shelf).
            if let (Some(cf), Some(b)) = (cx.cf.as_ref(), bil.get(i)) { if sea - *b as f64 > cf.fan_subaqueous_permille as f64 / 1000.0 * cp.shelf_m() { continue; } }
            let (ddx, ddy) = (gx as f64 * dx - f.x_m, gy as f64 * dx - f.y_m);
            let cth = if ddx == 0.0 && ddy == 0.0 { 1.0 } else { ang_diff(ddy.atan2(ddx), f.dir_rad).cos() };
            if cth <= 0.0 { continue; }
            // B6: seaward lobe; the profile is evaluated at d / (R(theta) / R).
            let d = ddx.hypot(ddy) / cth.sqrt();
            if d >= f.radius_m { continue; }
            let zf = if f.apex_m > sea + sf * f.radius_m { sea + (f.apex_m - sea) * (1.0 - d / f.radius_m).powf(1.5) } else { sea + sf * (f.radius_m - d) };
            if zf > z[i] {
                if fixed[i] == MARINE { fixed[i] = FREE; st.fan_land_gained_posts += 1; }
                if fixed[i] == GUIDE { fixed[i] = FREE; st.buried_channel_posts += 1; }
                z[i] = zf;
                over[i] = OVER_FAN;
                st.fan_raised_posts += 1;
            }
        } }
    }
    if let Some(observe) = observer.as_deref_mut() { observe(CoastMaterialStage::AfterFans, z, fixed, &over); }
    // Wave cut / beaches.
    let reach = (cp.cliff_retreat_mm.max(cp.beach_width_max_mm) as f64 / 1000.0) * 3.0;
    let cfp = cx.cf.as_ref();
    let profile = cfp.map(|c| c.cliff_profile == 1).unwrap_or(false);
    // Lane TD5: the soft-rock bluff needs the distance well beyond the beach (exact transform, 8 km cap).
    let bluff = rock_sc.is_some() && cx.fx.as_ref().map(|f| f.bluff == 1).unwrap_or(false);
    // Lane TF T3 / T4: shore points (sub-post contour or sea post centres), exposure, vector-propagated distance.
    let (tf_sub, tf_exp, tf_bl2) = (tf_on(|p| p.subpost), tf_on(|p| p.exposure), tf_on(|p| p.bluff_v2) && bluff);
    let mut tfst = crate::t1b_tf::TfCoastStats::default();
    let ex: Option<Vec<f32>> = if tf_exp { let t = tfp.unwrap(); Some(crate::t1b_tf::exposure(lat, &|i| fixed[i] == MARINE, t.fetch_dirs, t.fetch_cap_mm as f64 / 1000.0)) } else { None };
    if let Some(e) = &ex { let v: Vec<f64> = e.iter().filter(|&&x| x > 0.0).map(|&x| x as f64).collect(); tfst.exposure_shore_posts = v.len() as u32; tfst.exposure_mean = v.iter().sum::<f64>() / v.len().max(1) as f64; }
    // Lane TH H1: alongshore smoothing of the exposure before it is attached to the shore segments.
    let mut thst = crate::t1b_th::ThCoastStats::default();
    let ex: Option<Vec<f32>> = match (ex, cx.th.as_ref().filter(|t| t.expo_smooth == 1)) {
        (Some(e), Some(t)) => {
            let sh = crate::t1b_th::shore_posts(lat, &|i| fixed[i] == MARINE);
            let (es, est) = crate::t1b_th::smooth_exposure(lat, &e, &sh, t.expo_len_mm as f64 / 1000.0);
            thst.expo = Some(est);
            Some(es)
        }
        (e, _) => e,
    };
    let segs: Vec<(usize, [f64; 4], f32)> = if tf_sub { crate::t1b_tf::shore_segments(lat, z, sea, &|i| fixed[i] == MARINE, ex.as_deref()) }
        else if tf_exp || tf_bl2 { crate::t1b_tf::post_segments(lat, &|i| fixed[i] == MARINE, ex.as_deref()) } else { Vec::new() };
    tfst.shore_segments = segs.len() as u32;
    let vall = if tf_sub || tf_exp { Some(crate::t1b_tf::vprop_segments(lat, &segs, &|i| fixed[i] != MARINE, reach + dx)) } else { None };
    let dsea = match &vall { Some(v) => v.d.clone(), None => if profile { edt(lat, &|i| fixed[i] == MARINE, reach + dx) } else { chamfer(lat, &|i| fixed[i] == MARINE, reach + dx) } };
    let expo = |i: usize| -> f64 { match &vall { Some(v) if tf_exp => v.e[i] as f64, _ => 1.0 } };
    let dbl = if bluff && !tf_bl2 {
        if tf_sub { crate::t1b_tf::vprop_segments(lat, &segs, &|i| fixed[i] != MARINE, 8000.0).d } else { edt(lat, &|i| fixed[i] == MARINE, 8000.0) }
    } else { Vec::new() };
    // Lane CF C6a: no wave cut within a fan's radius + pad of its apex.
    let mut nocut = vec![false; if cfp.map(|c| c.fan_nocut == 1).unwrap_or(false) { nn } else { 0 }];
    if !nocut.is_empty() {
        let pad = cfp.unwrap().fan_nocut_pad_mm as f64 / 1000.0;
        for f in &fs {
            let rr = f.radius_m + pad;
            let (i0, i1) = (((f.x_m - rr) / dx).floor() as i64, ((f.x_m + rr) / dx).ceil() as i64);
            let (j0, j1) = (((f.y_m - rr) / dx).floor() as i64, ((f.y_m + rr) / dx).ceil() as i64);
            for gy in j0..=j1 { for gx in i0..=i1 {
                let Some(i) = lat.idx(gx, gy) else { continue };
                if (gx as f64 * dx - f.x_m).hypot(gy as f64 * dx - f.y_m) <= rr { nocut[i] = true; }
            } }
        }
    }
    let (wt, wts, wws) = cfp.map(|c| (c.talus_width_mm as f64 / 1000.0, c.talus_slope_ppm as f64 / 1e6, c.wall_slope_ppm as f64 / 1e6)).unwrap_or((0.0, 0.0, 0.0));
    let kmax = cp.cliff_k_max_permille as f64 / 1000.0;
    let (rr, bw, bwm, bs, pd) = (cp.cliff_retreat_mm as f64 / 1000.0, cp.beach_width_mm as f64 / 1000.0, cp.beach_width_max_mm as f64 / 1000.0, cp.beach_slope_ppm as f64 / 1e6, cp.platform_depth_mm as f64 / 1000.0);
    let mut cut = vec![false; nn];
    for i in 0..nn {
        if fixed[i] == MARINE || fixed[i] == LAKE || fixed[i] == CHANNEL || t1a_post[i] || !dsea[i].is_finite() || dsea[i] <= 0.0 { continue; }
        if !nocut.is_empty() && nocut[i] { continue; }
        let k = if over[i] == OVER_FAN { 3.0 } else if over[i] == OVER_EDIFICE { 0.6 } else { rock_k(i, z[i]) };
        let ei = expo(i);
        if k <= kmax && profile {
            // Lane TH H2: the cliff class; the height before the cut.
            if !st.th_uncut.is_empty() { st.th_uncut[i] = z[i]; thst.uncut_class_posts += 1; }
            // Lane CF C8: platform, talus, wall (lane TF T4: platform width x exposure).
            let pw = rr * k * ei;
            let d = dsea[i];
            // Lane TG G1: a cut post keeps the talus line extended seaward (sea + talus slope (d - pw)) for the shore field.
            if d < pw { cut[i] = true; if !st.tg_cut_ext.is_empty() { st.tg_cut_ext[i] = sea + wts * (d - pw); } continue; }
            let zp = if d < pw + wt { sea + wts * (d - pw) } else { sea + wts * wt + wws * (d - pw - wt) };
            if zp < z[i] {
                if z[i] - zp > 0.5 { if over[i] == OVER_NONE { over[i] = OVER_CUT; } st.cliff_cut_posts += 1; }
                if zp <= sea { cut[i] = true; if !st.tg_cut_ext.is_empty() { st.tg_cut_ext[i] = zp; } continue; }
                z[i] = zp;
            }
        } else if k <= kmax {
            if dsea[i] < rr * k * ei { cut[i] = true; if !st.tg_cut_ext.is_empty() { st.tg_cut_ext[i] = sea + (if wts > 0.0 { wts } else { 1.0 }) * (dsea[i] - rr * k * ei); } }
        } else {
            let w = (bw * k).min(bwm) * ei;
            if dsea[i] < w {
                let zb = sea + bs * dsea[i];
                if zb < z[i] { z[i] = zb; if fixed[i] == GUIDE { fixed[i] = FREE; } over[i] = OVER_BEACH; st.beach_posts += 1; if tfp.is_some() { tfst.beach_posts += 1; } }
            }
        }
    }
    let tg_sill = cx.tg.as_ref().map(|t| t.sill_grade == 1).unwrap_or(false);
    let mut tgst = crate::t1b_tg::TgCoastStats::default();
    if tf_bl2 {
        // Lane TF T2 soft_bluff_v2: distance through soft posts only, settled in increasing distance; bench at the top of the
        // first hard layer; never below the post the distance came from.
        let rsc = rock_sc.unwrap();
        let soft: Vec<bool> = (0..nn).map(|i| {
            if fixed[i] == MARINE || fixed[i] == LAKE || fixed[i] == CHANNEL || t1a_post[i] || !z[i].is_finite() { return false; }
            if !nocut.is_empty() && nocut[i] { return false; }
            let k = if over[i] == OVER_FAN { 3.0 } else if over[i] == OVER_EDIFICE { 0.6 } else { rock_k(i, z[i]) };
            k > kmax
        }).collect();
        let vs = crate::t1b_tf::vprop_segments(lat, &segs, &|i| soft[i], 8000.0);
        for &p in &vs.order {
            let i = p as usize;
            let d = vs.d[i];
            if !d.is_finite() || d <= 0.0 { continue; }
            let k = if over[i] == OVER_FAN { 3.0 } else { rock_k(i, z[i]) };
            let ei = if tf_exp { vs.e[i] as f64 } else { 1.0 };
            let w = (bw * k).min(bwm) * ei;
            if d < w { continue; }
            let sc = if over[i] == OVER_FAN { 0.625 } else { rsc(i, z[i]) };
            let zb = sea + bs * w + sc * (d - w);
            if zb >= z[i] { continue; }
            let mut zn = zb;
            if over[i] != OVER_FAN {
                if let Some(h) = crate::t1b_tf::hard_top(&|hh| rock_k(i, hh), zb, z[i], kmax) { if h > zn { zn = h; tfst.bluff_v2_bench_posts += 1; } }
            }
            let par = vs.parent[i];
            if par != u32::MAX {
                // Lane TG G2: the chain rises inland at the declared valley-floor slope instead of standing at the parent's height.
                let zpar = if tg_sill {
                    let (pp, ii) = (par as usize, i);
                    let l = (((pp % n) as f64 - (ii % n) as f64).hypot((pp / n) as f64 - (ii / n) as f64)) * dx;
                    z[pp] + cx.tg_floor_slope * l
                } else { z[par as usize] };
                if zpar > zn {
                    if tg_sill && zpar > z[par as usize] && z[par as usize] <= zn { tgst.sill_graded_posts += 1; }
                    if tg_sill { tgst.sill_graded_raise_max_m = tgst.sill_graded_raise_max_m.max(zpar - zn); }
                    zn = zpar; tfst.bluff_v2_chain_posts += 1;
                }
            }
            if zn < z[i] {
                if z[i] - zn > 0.5 { if over[i] == OVER_NONE { over[i] = OVER_CUT; } st.bluff_cut_posts += 1; tfst.bluff_v2_cut_posts += 1; }
                tfst.bluff_v2_cut_max_m = tfst.bluff_v2_cut_max_m.max(z[i] - zn);
                if fixed[i] == GUIDE { fixed[i] = FREE; }
                z[i] = zn;
            }
        }
    } else if bluff {
        // TD5: behind the beach (d >= w) a bluff at the soft class's limiting slope (posts the loop above skipped for
        // dsea beyond its reach are included: dbl reaches 8 km).
        let rsc = rock_sc.unwrap();
        for i in 0..nn {
            if fixed[i] == MARINE || fixed[i] == LAKE || fixed[i] == CHANNEL || t1a_post[i] || !dbl[i].is_finite() || dbl[i] <= 0.0 { continue; }
            if !nocut.is_empty() && nocut[i] { continue; }
            let k = if over[i] == OVER_FAN { 3.0 } else if over[i] == OVER_EDIFICE { 0.6 } else { rock_k(i, z[i]) };
            if k <= kmax { continue; }
            let w = (bw * k).min(bwm) * expo(i);
            let d = dbl[i];
            if d < w { continue; }
            let sc = if over[i] == OVER_FAN { 0.625 } else { rsc(i, z[i]) };
            let zb = sea + bs * w + sc * (d - w);
            if zb < z[i] {
                if z[i] - zb > 0.5 { if over[i] == OVER_NONE { over[i] = OVER_CUT; } st.bluff_cut_posts += 1; }
                if fixed[i] == GUIDE { fixed[i] = FREE; }
                z[i] = zb;
            }
        }
    }
    for i in 0..nn { if cut[i] { z[i] = sea - pd; fixed[i] = MARINE; over[i] = OVER_CUT; st.cliff_cut_posts += 1; if tfp.is_some() { tfst.cut_posts += 1; } } }
    if let Some(observe) = observer.as_deref_mut() { observe(CoastMaterialStage::AfterWaveCuts, z, fixed, &over); }
    bank_rule(lat, z, fixed, cx.cf.is_some());
    if let Some(observe) = observer.as_deref_mut() { observe(CoastMaterialStage::AfterBankConstraint, z, fixed, &over); }
    // Final fill: outlets are every fixed post (sea, lakes, channels).
    let outlet: Vec<bool> = (0..nn).map(|i| fixed[i] != FREE).collect();
    let zb: Vec<f64> = if ve_lat { z.to_vec() } else { Vec::new() };
    if st.dam_pre.len() == nn { st.dam_pre.copy_from_slice(z); }
    if stop == Some(crate::t1b_surface_stages::Stage::CoastFinal) { st.surface_capture=Some(crate::t1b_surface_stages::capture(z, fixed, st.seam.as_ref())); return Ok((over,vid)); }
    let r = if let Some(owner) = owner { crate::t1b_surface_stages::apply(owner, crate::t1b_surface_stages::Stage::CoastFinal, lat, z, fixed, st.seam.as_mut())?.0 } else if let Some(sm) = st.seam.as_mut() { crate::t1b_seam::pinned_fill_native(lat, z, fixed, sm, "coast_final").0 } else if tf_on(|p| p.report) {
        let fx2 = fixed.to_vec();
        let (r, _, rep) = crate::t1b_tf::fill_report(lat, z, &outlet, sea, &|i| fx2[i] == FREE, tfin.unwrap().margin_posts, 8000.0);
        tfst.final_fill = Some(rep);
        r
    } else { crate::t1b_rivernet::fill(lat, z, &outlet).0 };
    st.final_fill_raised_posts = r;
    if ve_lat { crate::t1b_ve::lake_stats(&zb, z, &st.ve_style, &mut vest); st.ve = Some(vest); }
    if tfp.is_some() { st.tf = Some(tfst); }
    if cx.tg.is_some() { tgst.cut_ext_posts = st.tg_cut_ext.iter().filter(|v| v.is_finite()).count() as u32; st.tg = Some(tgst); }
    if cx.th.is_some() {
        if !st.th_uncut.is_empty() {
            for i in 0..nn { if st.th_uncut[i].is_finite() { if fixed[i] == MARINE { thst.uncut_to_sea_posts += 1; } else if st.th_uncut[i] - z[i] > 0.5 { thst.uncut_lowered_posts += 1; } } }
        }
        st.th = Some(thst);
    }
    let _ = n;
    for i in 0..nn { if over[i] != OVER_EDIFICE { vid[i] = 255; } }
    if let Some(observe) = observer.as_deref_mut() { observe(CoastMaterialStage::AfterFinalFill, z, fixed, &over); }
    Ok((over, vid))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn material_stage_observation_preserves_coast_outputs() {
        let lat = EvalLat { gx0: 0, gy0: 0, n: 9, dx_mm: 32_000 };
        let count = lat.nodes();
        let cp = CoastParamsV1::v1();
        let mut cx = CoastCtx::default();
        cx.fans.push(Fan { x_m: 128.0, y_m: 128.0, area_m2: 1e6,
            e_m_per_my: 1.0, depth_m: 3.0, volume_m3: 20000.0,
            radius_m: 80.0, curve_code: CHANNEL, apex_m: 2.0, dir_rad: 0.0 });
        let t1a = vec![false; count]; let bil = vec![-1.0f32; count];
        let mut plain = vec![-1.0; count]; let mut observed = plain.clone();
        let mut fixed_plain = vec![MARINE; count]; let mut fixed_observed = fixed_plain.clone();
        let mut stats_plain = CoastStats::default(); let mut stats_observed = CoastStats::default();
        let expected = coast_and_volcano_td(&lat, 0.0, Some(&cp), None, &cx,
            &|_, _| 1.0, None, &t1a, &bil, &mut plain, &mut fixed_plain, &mut stats_plain);
        let mut snapshots = Vec::new();
        let actual = coast_and_volcano_tf_observed(&lat, 0.0, Some(&cp), None, &cx,
            &|_, _| 1.0, None, &t1a, &bil, &mut observed, &mut fixed_observed,
            &mut stats_observed, None, Some(&mut |stage, z, fixed, over| {
                assert_eq!(z.len(), count); assert_eq!(fixed.len(), count); assert_eq!(over.len(), count);
                snapshots.push((stage, z.to_vec()));
            }));
        assert_eq!(plain.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            observed.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
        assert_eq!(fixed_plain, fixed_observed); assert_eq!(expected, actual);
        assert_eq!(serde_json::to_value(stats_plain).unwrap(), serde_json::to_value(stats_observed).unwrap());
        assert_eq!(snapshots.iter().map(|s| s.0).collect::<Vec<_>>(), vec![
            CoastMaterialStage::BeforeFans, CoastMaterialStage::AfterFans,
            CoastMaterialStage::AfterWaveCuts, CoastMaterialStage::AfterBankConstraint,
            CoastMaterialStage::AfterFinalFill]);
        assert!(snapshots[1].1.iter().zip(&snapshots[0].1).any(|(a,b)| a>b));
    }


    fn site(age: f64, e: f64) -> VolcanoSite {
        let vp = VolcanoParamsV1::v1();
        let s = &vp.shapes[0];
        let (st, sb) = (s.s_top_ppm as f64 / 1e6, s.s_base_ppm as f64 / 1e6);
        let h = 2500.0;
        let lambda = h / (st - sb);
        let rc = 350.0;
        let r_main = rc + lambda * (st / sb).ln();
        let delta = e * age;
        let residual = delta >= h;
        let neck = 0.15 * h;
        VolcanoSite { id: "t".into(), kind: "continental_arc".into(), shape: "stratovolcano".into(), record_xy_m: [0.0, 0.0], vent_xy_m: [0.0, 0.0], moved_m: 0.0, record_cell: 0, vent_cell: 0,
            summit_record_m: 3500.0, base_m: 1000.0, h_m: h, lambda_m: lambda, r_main_m: r_main, r_influence_m: if residual { 0.35 * r_main + 100.0 } else { 1.3 * r_main },
            s_top: st, s_base: sb, crater_radius_m: rc, crater_depth_m: 150.0, age_my: age, e_vent_m_per_my: e, denudation_m: delta,
            stage: String::new(), breached: !residual && delta >= 150.0, breach_dir_rad: 0.0, gully_dirs_rad: (0..16).map(|k| k as f64 * 0.39).collect(), dike_dirs_rad: vec![0.3, 2.0, 3.5, 5.0], neck_relief_m: if residual { neck } else { 0.0 },
            lobe_phase_rad: vec![], lobe_norm: 1.0, cf: None, ve: None, ped_shift_m: None, built_summit_m: None }
    }
    fn z(v: &VolcanoSite, x: f64, y: f64) -> f64 { match edifice_at(&VolcanoParamsV1::v1(), v, x, y) { EdificeAt::Abs(z, f) => 1000.0 + (z - 1000.0) * f, EdificeAt::Add(d) => 1000.0 + d, EdificeAt::None => 1000.0 } }
    #[test]
    fn pristine_cone_has_rim_crater_and_concave_flanks() {
        let v = site(0.0, 300.0);
        assert!((z(&v, 350.0, 0.0) - 3500.0).abs() < 1e-9, "rim at the record summit");
        assert!((z(&v, 0.0, 0.0) - 3350.0).abs() < 1e-9, "crater floor 150 m down");
        // Concave: slope decreases outward; near-top slope ~ s_top (33 deg), main-flank foot ~ s_base.
        let s = |r: f64| (z(&v, r, 0.0) - z(&v, r + 10.0, 0.0)) / 10.0;
        assert!((s(360.0) - 0.65).abs() < 0.01 && s(2000.0) > s(5000.0) && (s(v.r_main_m - 20.0) - 0.09).abs() < 0.01);
        assert!((z(&v, v.r_main_m, 0.0) - 1000.0).abs() < 1e-6, "main flank ends at the base level");
        assert!(z(&v, 1.29 * v.r_main_m, 0.0) <= 1000.0 + 1e-9 && matches!(edifice_at(&VolcanoParamsV1::v1(), &v, 1.31 * v.r_main_m, 0.0), EdificeAt::None), "the foot fades out by 1.3 R");
    }
    #[test]
    fn old_record_is_a_residual_neck() {
        let v = site(34.0, 456.0);
        assert!(z(&v, 0.0, 0.0) > 1000.0 + 0.15 * 2500.0 - 1e-6 && z(&v, 2000.0, 1000.0) < 1000.0 + 1.0, "neck stands 375 m, no cone away from dikes");
        let mid = site(1.0, 300.0); // debug age, unit test only: 300 m denudation, breached, gullied
        assert!(mid.breached && z(&mid, 360.0, 0.0) < 3500.0 - 300.0 + 1e-6);
        let g_on = z(&mid, 3000.0 * 0.39f64.cos(), 3000.0 * 0.39f64.sin());
        let g_off = z(&mid, 3000.0 * 0.585f64.cos(), 3000.0 * 0.585f64.sin());
        assert!(g_on < g_off - 10.0, "a radial gully is incised ({g_on} vs {g_off})");
    }
    #[test]
    fn cf_lambda_and_scaled_cone() {
        let lam = solve_lambda(0.65, 9256.0, 2882.0).unwrap();
        assert!((0.65 * lam * (1.0 - (-9256.0 / lam).exp()) - 2882.0).abs() < 1e-6 && lam > 5000.0 && lam < 6000.0, "lambda {lam}");
        assert!(solve_lambda(0.65, 1000.0, 2882.0).is_none());
        // A scaled cone: rim at the summit, base at R_f, fade 0 beyond R_end, burial 1 near the vent.
        let mut v = site(0.0, 300.0);
        let cf = CfParamsV1::v1();
        let rf = 2500.0 / 0.15 / 2.0;
        v.lambda_m = solve_lambda(0.65, rf - 350.0, 2500.0).unwrap(); v.r_main_m = rf; v.r_influence_m = rf * 1.27;
        v.cf = Some(CfVolcano { wb_m: 2.0 * rf, r_foot_m: rf, q_built: 0.15, lambda_m: v.lambda_m, s_base: 0.1, r_end_max_m: rf * 1.27, clearance_m: 0.0, flows: vec![[0.0, 1000.0, 4000.0, 300.0, 30.0, 0.0, 2000.0, 0.0]], gullies: vec![[std::f64::consts::PI, 1000.0, 30.0, 0.0, 2000.0, 0.0]], vent_to_t1a_m: 1e9 });
        let vp = VolcanoParamsV1::v1();
        let (e, w, u) = cone_cf(&vp, &v, &cf, 350.0, 0.0, 1000.0).unwrap();
        assert!((e - 3500.0).abs() < 1e-6 && w == 1.0 && u == 1.0);
        let (e, _, _) = cone_cf(&vp, &v, &cf, 0.0, rf, 1000.0).unwrap();
        assert!((e - 1000.0).abs() < 1e-6, "base at R_f: {e}");
        assert!(cone_cf(&vp, &v, &cf, rf * 1.3, 0.0, 1000.0).is_none());
        let on = cone_cf(&vp, &v, &cf, 3000.0, 0.0, 1000.0).unwrap().0; let off = cone_cf(&vp, &v, &cf, 0.0, 3000.0, 1000.0).unwrap().0;
        assert!((on - off - 30.0).abs() < 1.0, "flow adds 30 m on its axis ({on} vs {off})");
        let gy = cone_cf(&vp, &v, &cf, -3000.0, 0.0, 1000.0).unwrap().0;
        assert!(off - gy > 5.0, "gully incised ({gy} vs {off})");
        let low = cone_cf(&vp, &v, &cf, 0.0, 3000.0, 400.0).unwrap().0;
        assert!((off - low - 500.0).abs() < 1e-6, "R1a: the cone follows a falling regional surface (smooth min, k 100 m)");
        assert!((cone_cf(&vp, &v, &cf, 0.0, 3000.0, 1500.0).unwrap().0 - off).abs() < 1e-9, "R1a never raises");
        assert!((cf.softmax(10.0, 10.0) - 10.0 - 25.0 * 2f64.ln()).abs() < 1e-9 && (cf.softmax(500.0, 0.0) - 500.0).abs() < 1e-6);
    }
    #[test]
    fn ring_plain_tail_starts_at_5_deg_and_fills_lows_to_t_max() {
        // Lane VH-I: off leaves the params JSON unchanged; the tail starts where the smooth CF slope is 0.09, its slope decays
        // to 0.02 over 0.5 R_f, the raise is min(T - g, T_max), never on marine / below-sea posts, and ends at its reach.
        let cf0 = CfParamsV1::v1();
        assert!(!serde_json::to_string(&cf0).unwrap().contains("ring_plain"));
        let mut cf = cf0.clone();
        cf.ring_plain = Some(RingPlainV1::v1());
        assert!(cf.validate().is_ok());
        let rp = cf.ring_plain.clone().unwrap();
        let mut v = site(0.0, 300.0);
        let rf = 2500.0 / 0.15 / 2.0;
        v.lambda_m = solve_lambda(0.65, rf - 350.0, 2500.0).unwrap(); v.r_main_m = rf; v.r_influence_m = rf * 1.27;
        v.cf = Some(CfVolcano { wb_m: 2.0 * rf, r_foot_m: rf, q_built: 0.15, lambda_m: v.lambda_m, s_base: 0.1, r_end_max_m: rf * 1.27, clearance_m: 0.0, flows: vec![], gullies: vec![], vent_to_t1a_m: 1e9 });
        let vp = VolcanoParamsV1::v1();
        let lam = v.lambda_m;
        let rho5 = lam * (0.65f64 / 0.09).ln();
        let r5 = 350.0 + rho5;
        let tail = |r: f64| ring_plain_at(&rp, &vp, &v, &cf, r, 0.0, 1000.0).unwrap();
        let zp5 = 3500.0 - 0.65 * lam * (1.0 - (-rho5 / lam).exp());
        assert!((tail(r5).0 - zp5).abs() < 1e-6 && (tail(r5).1 - 60.0).abs() < 1e-9, "tail starts on the profile at r5 with the full cap");
        let sl = |r: f64| (tail(r).0 - tail(r + 1.0).0) / 1.0;
        assert!((sl(r5 + 0.5) - 0.09).abs() < 1e-3 && (sl(r5 + 0.5 * rf + 10.0) - 0.02).abs() < 1e-6, "slope 0.09 -> 0.02 over 0.5 R_f");
        assert!((tail(r5 + 0.5 * rf + 0.125 * rf).1 - 30.0).abs() < 1e-6, "cap tapers linearly over 0.25 R_f");
        assert!(ring_plain_at(&rp, &vp, &v, &cf, r5 + 0.75 * rf + 1.0, 0.0, 1000.0).is_none(), "nothing beyond the reach");
        assert!(ring_plain_at(&rp, &vp, &v, &cf, r5 - 0.2 * rf - 1.0, 0.0, 1000.0).is_none(), "nothing inside the ramp-in");
        let t = tail(r5 + 1000.0).0;
        assert!((ring_plain_raise(&cf, &vp, &v, r5 + 1000.0, 0.0, 1000.0, t - 500.0, 0.0, FREE) - 60.0).abs() < 1e-9, "a deep moat is filled by T_max");
        assert!((ring_plain_raise(&cf, &vp, &v, r5 + 1000.0, 0.0, 1000.0, t - 20.0, 0.0, FREE) - 20.0).abs() < 1e-9, "a shallow low is filled to the tail");
        assert_eq!(ring_plain_raise(&cf, &vp, &v, r5 + 1000.0, 0.0, 1000.0, t + 5.0, 0.0, FREE), 0.0, "a ridge above the tail stands out");
        assert_eq!(ring_plain_raise(&cf, &vp, &v, r5 + 1000.0, 0.0, 1000.0, t - 500.0, 0.0, MARINE), 0.0);
        assert_eq!(ring_plain_raise(&cf, &vp, &v, r5 + 1000.0, 0.0, 1000.0, -5.0, 0.0, FREE), 0.0, "not below sea level");
        assert_eq!(ring_plain_raise(&cf0, &vp, &v, r5 + 1000.0, 0.0, 1000.0, t - 500.0, 0.0, FREE), 0.0, "off");
        assert!(ring_plain_reach(&rp, &vp, &v) > r5 && ring_plain_reach(&rp, &vp, &v) <= r5 + 0.75 * rf + 1e-6);
    }
    #[test]
    fn edt_matches_brute_force() {
        let lat = EvalLat { gx0: 0, gy0: 0, n: 23, dx_mm: 32_000 };
        let src = |i: usize| i == 5 || i == 23 * 17 + 3 || i == 23 * 20 + 21;
        let d = edt(&lat, &src, 1e9);
        for i in 0..23 * 23 {
            let (x, y) = ((i % 23) as f64, (i / 23) as f64);
            let b = [5usize, 23 * 17 + 3, 23 * 20 + 21].iter().map(|&j| (x - (j % 23) as f64).hypot(y - (j / 23) as f64) * 32.0).fold(f64::INFINITY, f64::min);
            assert!((d[i] - b).abs() < 1e-6, "{i}: {} vs {b}", d[i]);
        }
        assert!(edt(&lat, &|_| false, 1e9).iter().all(|v| v.is_infinite()));
    }
    #[test]
    fn chamfer_axis_exact() {
        let lat = EvalLat { gx0: 0, gy0: 0, n: 9, dx_mm: 32_000 };
        let d = chamfer(&lat, &|i| i == 0, 1e9);
        assert!((d[8] - 8.0 * 32.0).abs() < 1e-9 && (d[80] - 8.0 * 32.0 * std::f64::consts::SQRT_2).abs() < 1e-9);
    }
    #[test]
    fn flank_native_area_export_crater_and_continuity() {
        let mut v = site(1.0, 100.0);
        let rf = 5000.0;
        v.lambda_m = solve_lambda(v.s_top, rf - v.crater_radius_m, v.h_m).unwrap();
        v.r_main_m = rf; v.r_influence_m = 1.27 * rf;
        v.cf = Some(CfVolcano { wb_m: 2.0 * rf, r_foot_m: rf, q_built: v.h_m / (2.0 * rf), lambda_m: v.lambda_m,
            s_base: 0.1, r_end_max_m: 1.27 * rf, clearance_m: 0.0,
            flows: vec![[0.3, 900.0, 3500.0, 220.0, 40.0, 100.0, 2000.0, 0.0]],
            gullies: vec![], vent_to_t1a_m: 1e9 });
        let f = FlankSurface::build(&v);
        assert!(f.max_incision_m > 1.0 && f.eroded_m3 > 0.0);
        assert!(f.confluences > 0);
        assert!(f.runoff_relative_residual < 1e-12);
        assert!(f.volume_relative_residual < 1e-10, "{}", f.volume_relative_residual);
        for k in 0..720 {
            let a = k as f64 * std::f64::consts::TAU / 720.0;
            assert_eq!(f.at(v.crater_radius_m * a.cos(), v.crater_radius_m * a.sin()), 0.0);
            assert_eq!(f.at((rf + 256.0) * a.cos(), (rf + 256.0) * a.sin()), 0.0);
        }
        // No derivative jump across a reconstruction knot; the old V cross-section fails this.
        let (x, y, eps) = (2048.0, 1536.0, 0.01);
        let l = (f.at(x, y) - f.at(x-eps, y))/eps;
        let r = (f.at(x+eps, y) - f.at(x, y))/eps;
        assert!((l-r).abs() < 1e-4);
        assert!(FlankSurface::halfwidth(4e6) > FlankSurface::halfwidth(1e6));
        // Falsifier: zero erodibility gives exactly zero correction, despite the keyed initial deposits.
        let zero = FlankSurface::build_at(&v, 128.0, 12, 0.0);
        assert_eq!(zero.eroded_m3, 0.0);
        assert!(zero.correction.iter().all(|&x| x == 0.0));
        // Old angular templates cannot influence the new native drainage.
        v.gully_dirs_rad = vec![0.0, 1.0, 2.0];
        let same = FlankSurface::build(&v);
        assert_eq!(*same.correction, *f.correction);
        v.cf.as_mut().unwrap().flows.clear();
        let radial=FlankSurface::build(&v);
        for fraction in [0.25,0.4,0.55,0.7] {
            let rr=fraction*rf;let reference=radial.at(rr,0.0);
            for k in 0..720 {let a=k as f64*std::f64::consts::TAU/720.0;
                assert!((radial.at(rr*a.cos(),rr*a.sin())-reference).abs()<0.01,"uniform cone acquired azimuthal erosion");}
        }
    }

    #[test]
    fn flank_polar_flow_tracks_oblique_plane() {
        let (nr,nt,dx)=(41,256,100.0);let a=std::f64::consts::PI/8.0;let da=std::f64::consts::TAU/nt as f64;
        let z:Vec<f64>=(0..nr*nt).map(|i|{let r=(i/nt) as f64*dx;let t=(i%nt) as f64*da;10000.0-r*t.cos()*a.cos()-r*t.sin()*a.sin()}).collect();
        let fixed:Vec<bool>=(0..nr*nt).map(|i|i/nt<=1||i/nt>=nr-2).collect();
        let route=polar_route(nr,nt,dx,&z,&fixed);let i=20*nt;let mut angle=0.0;
        for (j,w,_) in route.receivers[i] {if w>0.0 {let r=(j/nt) as f64*dx;let t=(j%nt) as f64*da;angle+=w*(r*t.sin()).atan2(r*t.cos()-2000.0);}}
        assert!((angle-a).abs().to_degrees()<0.05,"{}",angle.to_degrees());
    }

    #[test]
    #[ignore = "real C6 record convergence probe; supply VZ_FLANK_INPUT and VZ_FLANK_REPORT"]
    fn flank_real_record_convergence() {
        let input=std::env::var("VZ_FLANK_INPUT").expect("real layer path");
        let j:serde_json::Value=serde_json::from_slice(&std::fs::read(&input).unwrap()).unwrap();
        let a=&j["header"]["volcanoes"][0];let num=|k:&str| a[k].as_f64().unwrap();
        let mut v=site(0.0,0.0);
        v.h_m=num("h_m");v.denudation_m=num("denudation_m");v.summit_record_m=num("summit_record_m");
        v.base_m=num("base_m");v.lambda_m=num("lambda_m");v.s_top=num("s_top");
        v.crater_radius_m=num("crater_radius_m");v.crater_depth_m=num("crater_depth_m");v.built_summit_m=a["built_summit_m"].as_f64();
        let cf=&a["cf"];let rf=cf["r_foot_m"].as_f64().unwrap();
        v.cf=Some(CfVolcano{wb_m:2.0*rf,r_foot_m:rf,q_built:0.15,lambda_m:v.lambda_m,s_base:0.1,r_end_max_m:rf*1.27,clearance_m:0.0,
            flows:serde_json::from_value(cf["flows"].clone()).unwrap(),gullies:vec![],vent_to_t1a_m:1e9});
        let mut costs=Vec::new();let mut fields=Vec::new();
        for (dx,steps) in [(128.0,60),(64.0,60),(32.0,60),(64.0,120)] {
            let start=std::time::Instant::now();let f=FlankSurface::build_at(&v,dx,steps,20.0);
            costs.push(serde_json::json!({"spacing_m":dx,"steps":steps,"seconds":start.elapsed().as_secs_f64(),"ledger":f}));fields.push(f);
        }
        let mut rows=Vec::new();
        for fraction in [0.25,0.4,0.55,0.7] {
            for k in 0..720 {let a=k as f64*std::f64::consts::TAU/720.0;let (x,y)=(fraction*rf*a.cos(),fraction*rf*a.sin());
                rows.push(serde_json::json!({"r_fraction":fraction,"theta_rad":a,"correction_m":fields.iter().map(|f| f.at(x,y)).collect::<Vec<_>>()}));}
        }
        let report=serde_json::json!({"input":input,"seed":42,"realization":0,"costs":costs,"samples":rows});
        std::fs::write(std::env::var("VZ_FLANK_REPORT").unwrap(),serde_json::to_vec_pretty(&report).unwrap()).unwrap();
    }

}
