//! Lane MEANDER-Q (2026-10-22): opt-in key `cl_mq=1` (needs `cl=1`; replaces lane CHANNEL-CENTRELINE's declared
//! decorative sine belt `cl_mn=1`, which it refuses alongside). Absent (`params.cl.mq` None) = byte-identical.
//! LABEL: EXPLORATORY. A KINEMATIC meander-migration model (curvature-driven, upstream-lagged), not a hydrodynamic or
//! sediment-conserving one: no flow field, bank shear, bar deposition, sediment budget or water routing is solved. It is a
//! reduced-order stand-in on the way to RIVER-PILOT / RIVER-INTEGRATE (docs/plans/RIVER_EVOLUTION_PLAN.md), better than the
//! sine belt because bends grow, migrate, skew downstream and cut off from a process rule instead of being drawn.
//!
//! What it does, per T1a segment (one reach record of the chart; every window that sees the reach computes the same thing, so
//! seams are exact):
//! 1. stations every `station` widths along the T1a planform curve (lane CHANNEL-CENTRELINE's carrier), small keyed lateral
//!    perturbations (world-keyed: vz.keyed-rng root(world seed) -> layer(this generator id) -> address(mq/<T0 segment address>);
//!    never the chart recipe, the window or the post lattice);
//! 2. Howard and Knutson (1984) / Ikeda, Parker and Sawai (1981) kinematic migration: nominal rate R0 = w^2 C (C the
//!    centreline curvature, w the bankfull width), lagged rate R1 = Omega R0 + Gamma int R0(s - xi) exp(-alpha xi) / int
//!    exp(-alpha xi) (upstream exponential weighting, alpha = 2 Cf / D), bank migration E_b f R1 along the outer normal;
//!    f = regime factor x valley-wall factor x end pin taper;
//! 3. neck cutoffs (banks touch) and chute cutoffs (short path across the swept point bar with a large gradient advantage);
//!    each cut loop is an oxbow feature with its cutoff age;
//! 4. the floodplain (T1b floor extent) is the area the channel swept (plus a margin), so the floor follows where the
//!    river has been; reaches that do not migrate (braided by the slope-discharge criterion, or pinned ends) keep the
//!    T1b floor band; the floor stands an inset depth D(Q) above the channel bed.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use vz_keyed_rng::Scope;

pub const LABEL: &str = "EXPLORATORY lane MEANDER-Q 2026-10-22; kinematic meander migration (Howard-Knutson 1984 type), not hydrodynamic and not sediment-conserving; constants carry [V] verified / [R] recalled from literature, unverified / [A] project assumption";
pub const PURPOSE_PERTURB: &str = "t1b.mq.perturb";
pub const GENERATOR_ID: &str = "vz.pf-region.t1b_mq.kinematic_meander_migration_v1";
pub const MQ_RULE: &str = "kinematic_meander_migration_v1 (lane MEANDER-Q; EXPLORATORY; replaces sine_generated_meander_belt_v1): per T1a segment, on the T1a planform curve of the segment (the cl=1 carrier), between its first chart node and two node spacings past its last (both ends pinned: junctions, mouths, chart-edge entries; erodibility tapers to 0 over end_taper widths [A]); \
STATIONS every station (0.5) x w, w the local bankfull width (clamped 4-64 m) [A numerical; synthetic 30 km reach, 6 ka: station 1.0 -> 0.5 -> 0.25 widths gives sinuosity 3.35 -> 3.06 -> 2.77 unconfined and 1.71 -> 1.31 -> 1.25 confined (h_v 100 m), wavelength p50 15.9 -> 17.8 -> 17.6 widths; cutoff counts are not resolution-converged], each station's hydraulics, regime, walls and pin come from its valley position (the nearest carrier sample, followed by hill-climbing as bends drift down-valley); T1a nodes are placed on the first station whose valley position (running maximum along the channel) reaches theirs; stations closer than 0.6 local spacings merge, gaps over 1.6 get midpoints [A numerical]; \
INITIAL STATE: the carrier plus a lateral offset sigma w N_i (N_i = vz.keyed-rng normal draw, root(world seed from the T1a recipe) -> layer(blake3 of the generator id) -> address(mq/<T1a segment address, T0 L-grid element ids>) -> draw(t1b.mq.perturb, station index, realization); world-keyed, never chart- or window-keyed), smoothed by two [1 2 1]/4 passes [A]; \
HYDRAULICS per carrier sample from the governing T1a node (smoothed by a moving average over +-500 m [A numerical]): Q = A r_ref / yr (the VALLEY-FLOOR-FIX / pf-detail B6 effective mean annual discharge), w = 3.5 Q^0.5, D = 0.27 Q^0.4 (pf-detail B6 hydraulic geometry [R]; note earth-numbers C: 7.2 Q^0.5 for mean annual Q globally, a ruling is open), valley half-width h_v = the T1b floor half-width 0.5 c_w sqrt(A) capped (and the volcano edifice corridor cap where given), bed slope S = the T1b bed drop of the node's link over its carrier length; \
MIGRATION (Howard and Knutson 1984 WRR 20:1659; Ikeda, Parker and Sawai 1981 JFM 112:363; form as in meanderpy, Sylvester et al. 2019 [R]): curvature C_i from the turning angle between station chords, smoothed by two [1 2 1]/4 passes (removes the grid-scale anti-diffusive mode of the explicit scheme [A numerical]); R0 = w^2 C with w C saturated at +-0.4 (Hickin and Nanson 1975/1984: migration rate peaks near R/w = 2-3 [R]; the decline in tighter bends is not modelled [A]); R1_i = Omega R0_i + Gamma (sum_j<=i R0_j exp(-alpha (s_i - s_j)) ds_j) / (sum exp(-alpha (s_i - s_j)) ds_j) by the exact exponential recursion, Omega = -1, Gamma = 2.5 [R HK84]; alpha = 2 Cf / D, Cf = g n^2 / D^(1/3) with Manning n [R Chow 1959: natural main channels 0.030-0.040]; bank migration rate E_b f R1 (m/yr) along the outer normal, E_b [R/A: gives a peak nominal rate of E_b x 0.4 widths per year at w C = 0.4 (R/w = 2.5); observed vegetated meander migration is of order 1 % of the width per year (Hickin and Nanson 1984 J Hydraul Eng 110:1557; Ielpi and Lapotre 2020 Nat Geosci 13:82 report about tenfold faster before land plants), recalled]; f = f_regime f_wall f_pin; explicit steps with dt <= min_i 0.25 ds_i / |rate_i|, dt <= min_i 0.4 ds_i^2 / (0.17 pi^2 E_b f_i w_i^2) (stability of the smoothed scheme), dt <= dt_max; integrated over the migration time T (age_yr, a documented STAND-IN until a landscape-age field exists [A]: default 11.7 ka, the Holocene, the age of most active mid-latitude floodplains [R]; the RIVER_EVOLUTION_PLAN 2026-10-15 pilot candidate for the final refinement window is 2,000 years, key cl_mq_age); \
SUBMERGED reaches (the T1a curve continues across the shelf / estuary arm: the node's bed at or below sea level + 0.5 m, or a T1a mouth node) do not migrate (f = 0, smoothed with the hydraulics) and keep the band; \
REGIME (slope-discharge criterion, Leopold and Wolman 1957 USGS PP 282-B [R]): S_LW = c_LW Q_bf^e_LW (c_LW = 0.012, e_LW = -0.44 with Q_bf in m3/s as commonly quoted [R, units to verify]), Q_bf = r_bf Q (bankfull about the 1.5-year flood, a few times the mean flow [A r_bf = 3]); index P = S / S_LW; f_regime = smoothstep((P_hi - P) / (P_hi - P_lo)) (P_lo = 1, P_hi = 2 [A transition band]): P > P_hi = braided or steep: no migration, published as such; regime=0 (key cl_mq_regime=0) switches the gate off and is a DIAGNOSTIC, not a physical setting; \
CONFINEMENT: a bank moving outward beyond the valley half-width (|offset from the carrier| + w/2 > h_v) erodes the valley wall: f_wall = wall_factor x k_rock (k_rock = the T0 erosion law's erodibility multiplier of the rock exposed at the node's bed level, from the T1b node context when the window binary supplies it, else 1) [A: bedrock walls an order of magnitude slower than alluvial banks]; linear ramp over w/2; \
CUTOFFS: neck when two stations more than 3 crdist apart along the channel come within crdist = neck x w (banks touch; meanderpy default 1.5 W [R]); chute when two stations come within chute x w, the channel path between them is at least chute_ratio times their distance, and the midpoint of the short path lies on ground the channel swept earlier (a point bar) [A]; neck preferred, the largest loop cut first; the cut loop is published as an oxbow (polyline, kind, age = T - cutoff time); distance-based cutoff detection is a kinematic closure (no breach hydraulics), declared; \
FLOODPLAIN: the union of discs of radius w/2 + margin x w around every swept station position (16 m occupancy cells, marked whenever the reach moved half a cell), plus on non-migrating stretches (braided by the regime factor or submerged: f_regime < 1) the T1b floor band 0.5 c_w sqrt(A) (capped as T1b; NOT the volcano corridor cap, which only places the migration walls) x (1 - f_regime) around the carrier (pinned ends get no band: a band there drew a disc of floor at every junction); T1b floor posts are the window posts inside it, valued from their nearest channel post as before (levee/backswamp rule when vf_floor) plus the inset depth D(Q) of that post (inset channel, Leopold and Maddock 1953 depth [R]); \
STATISTICS: half-waves between successive sign changes of the smoothed curvature with chord > 2 w; wavelength = 2 x chord, amplitude = largest distance from the chord; sinuosity = channel length / carrier length";

const G: f64 = 9.81;
const CELL_M: f64 = 16.0;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct MqParamsV1 {
    pub on: u32,
    /// Migration time (years), STAND-IN until a landscape-age field exists.
    pub age_yr: i64,
    /// E_b (1e-6 / yr).
    pub bank_rate_ppm_per_yr: i64,
    pub omega_permille: i64,
    pub gamma_permille: i64,
    /// Manning n (1e-6).
    pub manning_n_ppm: i64,
    /// Initial lateral perturbation sigma (permille of w).
    pub perturb_widths_permille: i64,
    /// Station spacing (permille of the reach's smallest width).
    pub station_widths_permille: i64,
    pub neck_widths_permille: i64,
    pub chute_widths_permille: i64,
    pub chute_ratio_permille: i64,
    pub wall_factor_permille: i64,
    /// 0 = regime gate off (diagnostic), 1 = Leopold-Wolman.
    pub regime: u32,
    pub lw_coeff_ppm: i64,
    pub lw_exp_permille: i64,
    pub bankfull_ratio_permille: i64,
    pub regime_lo_permille: i64,
    pub regime_hi_permille: i64,
    pub margin_widths_permille: i64,
    /// 1 = floor stands the inset depth D(Q) above the channel bed.
    pub inset: u32,
    pub end_taper_widths_permille: i64,
    pub dt_max_yr: i64,
    /// Saturation of the nominal rate at w C = this (permille): tighter bends migrate no faster.
    pub curv_sat_permille: i64,
    pub rule: String,
    pub label: String,
}

impl MqParamsV1 {
    pub fn v1() -> Self {
        MqParamsV1 {
            on: 1, age_yr: 11_700, bank_rate_ppm_per_yr: 50_000, omega_permille: -1000, gamma_permille: 2500, manning_n_ppm: 35_000,
            perturb_widths_permille: 50, station_widths_permille: 500, neck_widths_permille: 1500, chute_widths_permille: 3000,
            chute_ratio_permille: 4000, wall_factor_permille: 100, regime: 1, lw_coeff_ppm: 12_000, lw_exp_permille: -440,
            bankfull_ratio_permille: 3000, regime_lo_permille: 1000, regime_hi_permille: 2000, margin_widths_permille: 500, inset: 1,
            end_taper_widths_permille: 2000, dt_max_yr: 10, curv_sat_permille: 400, rule: MQ_RULE.into(), label: LABEL.into(),
        }
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.on != 1 { return Err("cl_mq: on must be 1 (absent = off)".into()); }
        if !(0..=1_000_000).contains(&self.age_yr) { return Err("cl_mq_age must be 0..1e6 years".into()); }
        if !(0..=1_000_000).contains(&self.bank_rate_ppm_per_yr) { return Err("cl_mq_rate must be 0..1e6 ppm/yr".into()); }
        if !(-5000..=0).contains(&self.omega_permille) || !(0..=10_000).contains(&self.gamma_permille) { return Err("cl_mq omega/gamma".into()); }
        if !(5_000..=200_000).contains(&self.manning_n_ppm) { return Err("cl_mq_n must be 0.005..0.2".into()); }
        if !(0..=1000).contains(&self.perturb_widths_permille) || !(200..=2000).contains(&self.station_widths_permille) { return Err("cl_mq perturb/station".into()); }
        if !(500..=5000).contains(&self.neck_widths_permille) || self.chute_widths_permille < self.neck_widths_permille || self.chute_ratio_permille < 1000 { return Err("cl_mq cutoff constants".into()); }
        if !(0..=1000).contains(&self.wall_factor_permille) || self.regime > 1 { return Err("cl_mq wall/regime".into()); }
        if self.lw_coeff_ppm <= 0 || self.bankfull_ratio_permille < 1000 || self.regime_lo_permille <= 0 || self.regime_hi_permille <= self.regime_lo_permille { return Err("cl_mq regime constants".into()); }
        if !(0..=5000).contains(&self.margin_widths_permille) || self.inset > 1 || !(0..=20_000).contains(&self.end_taper_widths_permille) || !(1..=100).contains(&self.dt_max_yr) || !(50..=5000).contains(&self.curv_sat_permille) { return Err("cl_mq constants".into()); }
        Ok(())
    }
    /// One `cl_mq_<name>` key.
    pub fn set(&mut self, name: &str, v: i64) -> Result<(), String> {
        match name {
            "age" => self.age_yr = v, "rate" => self.bank_rate_ppm_per_yr = v, "omega" => self.omega_permille = v, "gamma" => self.gamma_permille = v,
            "n" => self.manning_n_ppm = v, "perturb" => self.perturb_widths_permille = v, "station" => self.station_widths_permille = v,
            "neck" => self.neck_widths_permille = v, "chute" => self.chute_widths_permille = v, "chute_ratio" => self.chute_ratio_permille = v,
            "wall" => self.wall_factor_permille = v, "regime" => { if !(0..=1).contains(&v) { return Err("cl_mq_regime is 0 or 1".into()); } self.regime = v as u32 }
            "lw_c" => self.lw_coeff_ppm = v, "lw_e" => self.lw_exp_permille = v, "qbf" => self.bankfull_ratio_permille = v,
            "p_lo" => self.regime_lo_permille = v, "p_hi" => self.regime_hi_permille = v, "margin" => self.margin_widths_permille = v,
            "inset" => { if !(0..=1).contains(&v) { return Err("cl_mq_inset is 0 or 1".into()); } self.inset = v as u32 }
            "taper" => self.end_taper_widths_permille = v, "dt" => self.dt_max_yr = v, "sat" => self.curv_sat_permille = v,
            _ => return Err(format!("unknown key cl_mq_{name}")),
        }
        Ok(())
    }
    fn f(v: i64, s: f64) -> f64 { v as f64 / s }
    pub fn regime_factor(&self, slope: f64, q: f64) -> (f64, f64) {
        if self.regime == 0 { return (1.0, f64::NAN); }
        let qbf = (Self::f(self.bankfull_ratio_permille, 1000.0) * q).max(1e-6);
        let s_lw = Self::f(self.lw_coeff_ppm, 1e6) * qbf.powf(Self::f(self.lw_exp_permille, 1000.0));
        let p = slope.max(0.0) / s_lw;
        let (lo, hi) = (Self::f(self.regime_lo_permille, 1000.0), Self::f(self.regime_hi_permille, 1000.0));
        (smooth01((hi - p) / (hi - lo)), p)
    }
}

/// Key parsing (bin): `cl_mq=0|1` and `cl_mq_<name>=<int>` (need cl=1 before; cl_mq=1 clears a cl_mn sine belt).
pub fn parse_key(c: &mut crate::t1b_cl::ClParamsV1, k: &str, v: i64) -> Result<(), String> {
    if k == "cl_mq" {
        match v {
            1 => { if c.mq.is_none() { c.mq = Some(MqParamsV1::v1()); } c.meander = 0; c.meander_rule = String::new(); }
            0 => c.mq = None,
            _ => return Err("cl_mq is 0 or 1".into()),
        }
        return Ok(());
    }
    let name = k.strip_prefix("cl_mq_").ok_or_else(|| format!("unknown key {k}"))?;
    c.mq.as_mut().ok_or_else(|| format!("{k} needs cl_mq=1 given before it"))?.set(name, v)
}

/// Lane MEANDER-Q: the T1b node context's rock erodibility multiplier at each T1a node's level, for the valley-wall factor
/// (window binary only; None = 1 everywhere). Called by the bin after the beds; a no-op without cl_mq.
pub fn attach_bank_rock(beds: &mut crate::t1b_rivernet::Beds, net: &crate::t1b_rivernet::Network, nctx: &crate::t1b_rivernet::NodeCtx, skel: &crate::t1b_view::SkeletonView, p: &crate::t1b_rivernet::RiverNetParamsV1) {
    use crate::t1b_litho::ExposedLithology;
    if !crate::t1b_cl::cl(p).map(|c| c.mq_on()).unwrap_or(false) { return; }
    let n = skel.kind.len();
    let mut k = vec![1.0f64; n];
    for i in 0..n.min(net.kind.len()) {
        if !(skel.kind[i] == 1 || skel.kind[i] == 3) { continue; }
        let z = beds.t1a_level.as_ref().map(|t| t[i]).filter(|v| v.is_finite()).unwrap_or(skel.trunk_level_m[i] as f64);
        let r = nctx.rock.at(i, z).1;
        if r.is_finite() && r > 0.0 { k[i] = r; }
    }
    beds.mq_rock = Some(k);
}

pub fn smooth01(x: f64) -> f64 { let t = x.clamp(0.0, 1.0); t * t * (3.0 - 2.0 * t) }

fn binomial2(v: &mut [f64]) {
    let n = v.len();
    if n < 3 { return; }
    for _ in 0..2 {
        let mut prev = v[0];
        for i in 1..n - 1 { let cur = v[i]; v[i] = 0.25 * prev + 0.5 * cur + 0.25 * v[i + 1]; prev = cur; }
    }
}

/// Per-reach inputs (one T1a segment).
pub struct SegIn<'a> {
    pub seg: u32,
    pub address: &'a str,
    pub seed: u64,
    pub realization: u32,
    /// Carrier (T1a curve) samples, arc length and unit left normals.
    pub cx: &'a [[f64; 2]],
    pub cs: &'a [f64],
    pub cn: &'a [[f64; 2]],
    /// Per carrier sample: bankfull width, depth, valley half-width, mean discharge, bed slope, rock factor.
    pub w: Vec<f64>,
    pub dep: Vec<f64>,
    pub hv: Vec<f64>,
    pub q: Vec<f64>,
    pub slope: Vec<f64>,
    pub krock: Vec<f64>,
    /// Per carrier sample: the T1b floor half-width without the volcano corridor cap (the band of non-migrating stretches).
    pub hband: Vec<f64>,
    /// Per carrier sample: 1 where the reach's bed is at or below sea level (shelf / estuary arm), else 0 (smoothed).
    pub wet: Vec<f64>,
    /// Simulated carrier sample range.
    pub k0: usize,
    pub k1: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct Oxbow {
    pub segment: u32,
    pub kind: &'static str,
    pub age_yr: f64,
    pub length_m: f64,
    pub neck_m: f64,
    /// Chart metres, 0.1 m.
    pub pts: Vec<[f64; 2]>,
}

/// Statistics of one reach.
#[derive(Clone, Debug, Default, Serialize)]
pub struct MqSeg {
    pub segment: u32,
    pub address: String,
    pub carrier_m: f64,
    pub channel_m: f64,
    pub sinuosity: f64,
    pub width_min_m: f64,
    pub width_max_m: f64,
    pub slope_p50: f64,
    pub regime_index_p50: f64,
    pub migrating_share: f64,
    pub years: f64,
    pub steps: u64,
    pub stations_final: u64,
    pub neck_cutoffs: u32,
    pub chute_cutoffs: u32,
    pub wall_station_steps: u64,
    pub belt_half_max_m: f64,
    pub halfwaves: u32,
    pub wavelength_widths_p10: f64,
    pub wavelength_widths_p50: f64,
    pub wavelength_widths_p90: f64,
    pub wavelength_cv: f64,
    pub amplitude_widths_p10: f64,
    pub amplitude_widths_p50: f64,
    pub amplitude_widths_p90: f64,
    pub amplitude_cv: f64,
}

/// One reach after migration.
pub struct SegOut {
    pub seg: u32,
    pub x: Vec<[f64; 2]>,
    pub lab: Vec<f64>,
    pub oxbows: Vec<Oxbow>,
    pub cells: HashMap<(i64, i64), f32>,
    pub band: Vec<[f32; 3]>,
    pub stats: MqSeg,
}

impl SegOut {
    /// Channel point at carrier label l (clamped to the simulated range).
    pub fn point_at(&self, l: f64) -> [f64; 2] {
        let n = self.lab.len();
        if l <= self.lab[0] { return self.x[0]; }
        if l >= self.lab[n - 1] { return self.x[n - 1]; }
        let j = self.lab.partition_point(|&v| v < l).max(1);
        if self.lab[j] == l { return self.x[j]; }
        let (a, b) = (self.lab[j - 1], self.lab[j]);
        let f = if b > a { (l - a) / (b - a) } else { 0.0 };
        [self.x[j - 1][0] + f * (self.x[j][0] - self.x[j - 1][0]), self.x[j - 1][1] + f * (self.x[j][1] - self.x[j - 1][1])]
    }
    /// The channel between labels l0 <= l1 (end points interpolated): every station from the first reaching l0 to the first
    /// reaching l1, so loops that fold back (constant running-maximum label) stay in the link.
    pub fn between(&self, l0: f64, l1: f64) -> Vec<[f64; 2]> {
        let n = self.lab.len();
        let j0 = self.lab.partition_point(|&v| v < l0).min(n);
        let j1 = self.lab.partition_point(|&v| v < l1).min(n);
        let mut out = vec![self.point_at(l0)];
        let start = if j0 < n && self.lab[j0] == l0 { j0 + 1 } else { j0 };
        for i in start..j1 { out.push(self.x[i]); }
        if l1 > l0 { out.push(self.point_at(l1)); }
        out
    }
    /// The channel from label l0 to station index j (inclusive).
    pub fn to_index(&self, l0: f64, j: usize) -> Vec<[f64; 2]> {
        let n = self.lab.len();
        let j0 = self.lab.partition_point(|&v| v < l0).min(n);
        let mut out = vec![self.point_at(l0)];
        let start = if j0 < n && self.lab[j0] == l0 { j0 + 1 } else { j0 };
        for i in start..=j.min(n - 1) { out.push(self.x[i]); }
        out
    }
    /// Index of the channel station nearest p among those at or after the first reaching label l0.
    pub fn nearest_from(&self, p: [f64; 2], l0: f64) -> usize {
        let n = self.lab.len();
        let j0 = self.lab.partition_point(|&v| v < l0).min(n - 1);
        let mut best = (f64::INFINITY, j0);
        for i in j0..n {
            let d = (self.x[i][0] - p[0]).hypot(self.x[i][1] - p[1]);
            if d < best.0 { best = (d, i); }
        }
        best.1
    }
}

fn interp(v: &[f64], l: f64) -> f64 {
    let n = v.len();
    if l <= 0.0 { return v[0]; }
    let i = l.floor() as usize;
    if i + 1 >= n { return v[n - 1]; }
    let f = l - i as f64;
    v[i] + f * (v[i + 1] - v[i])
}

fn cell_of(p: [f64; 2]) -> (i64, i64) { ((p[0] / CELL_M).floor() as i64, (p[1] / CELL_M).floor() as i64) }

/// World-keyed scope of a reach.
pub fn reach_scope(seed: u64, address: &str) -> Result<Scope, String> {
    let gd = format!("blake3:{}", blake3::hash(GENERATOR_ID.as_bytes()).to_hex());
    Scope::root(seed).layer(&gd).and_then(|s| s.address(&format!("mq/{address}"))).map_err(|e| e.to_string())
}

/// Run the migration on one reach.
pub fn run_segment(inp: &SegIn, mp: &MqParamsV1) -> Result<SegOut, String> {
    let (k0, k1) = (inp.k0, inp.k1.min(inp.cx.len() - 1));
    if k1 <= k0 + 1 { return Err(format!("segment {}: empty range", inp.seg)); }
    let lab_at_s = |s: f64| -> f64 {
        // carrier label (fractional sample index) at arc length s
        let j = inp.cs.partition_point(|&v| v <= s).clamp(1, inp.cs.len() - 1);
        let (a, b) = (inp.cs[j - 1], inp.cs[j]);
        (j - 1) as f64 + if b > a { ((s - a) / (b - a)).clamp(0.0, 1.0) } else { 0.0 }
    };
    let carrier_at = |l: f64| -> ([f64; 2], [f64; 2]) {
        let i = (l.floor() as usize).min(inp.cx.len() - 2);
        let f = (l - i as f64).clamp(0.0, 1.0);
        let p = [inp.cx[i][0] + f * (inp.cx[i + 1][0] - inp.cx[i][0]), inp.cx[i][1] + f * (inp.cx[i + 1][1] - inp.cx[i][1])];
        let nn = [inp.cn[i][0] + f * (inp.cn[i + 1][0] - inp.cn[i][0]), inp.cn[i][1] + f * (inp.cn[i + 1][1] - inp.cn[i][1])];
        let l2 = nn[0].hypot(nn[1]).max(1e-12);
        (p, [nn[0] / l2, nn[1] / l2])
    };
    let (s0, s1) = (inp.cs[k0], inp.cs[k1]);
    let wmin = (k0..=k1).map(|k| inp.w[k]).fold(f64::INFINITY, f64::min).max(1.0);
    let wmax = (k0..=k1).map(|k| inp.w[k]).fold(0.0, f64::max);
    let station = mp.station_widths_permille as f64 / 1000.0;
    // local station spacing: station x the bankfull width at the valley position (narrow head reaches do not force a fine
    // spacing, and the stability step, on the whole reach)
    let dsl = |l: f64| (station * interp(&inp.w, l)).clamp(4.0, 64.0);
    let ds0 = dsl(k0 as f64).min(dsl(k1 as f64));
    let mut svals = vec![s0];
    while *svals.last().unwrap() < s1 { let sl = *svals.last().unwrap(); svals.push((sl + dsl(lab_at_s(sl))).min(s1)); }
    if svals.len() >= 3 && s1 - svals[svals.len() - 2] < 0.5 * ds0 { let m = svals.len() - 2; svals.remove(m); }
    let ns = svals.len() - 1;
    let sigma = mp.perturb_widths_permille as f64 / 1000.0;
    let eb = mp.bank_rate_ppm_per_yr as f64 / 1e6;
    let (omega, gamma) = (mp.omega_permille as f64 / 1000.0, mp.gamma_permille as f64 / 1000.0);
    let nman = mp.manning_n_ppm as f64 / 1e6;
    let (neck, chute, chute_ratio) = (mp.neck_widths_permille as f64 / 1000.0, mp.chute_widths_permille as f64 / 1000.0, mp.chute_ratio_permille as f64 / 1000.0);
    let wall = mp.wall_factor_permille as f64 / 1000.0;
    let margin = mp.margin_widths_permille as f64 / 1000.0;
    let taper_w = mp.end_taper_widths_permille as f64 / 1000.0;
    let csat = mp.curv_sat_permille as f64 / 1000.0;
    // per carrier sample regime factor and index
    let mut freg = vec![0.0f64; inp.cx.len()];
    let mut pidx = vec![f64::NAN; inp.cx.len()];
    for k in k0..=k1 { let (f, p) = mp.regime_factor(inp.slope[k], inp.q[k]); freg[k] = f * (1.0 - inp.wet[k]).clamp(0.0, 1.0); pidx[k] = p; }
    // initial stations with keyed perturbation
    let scope = reach_scope(inp.seed, inp.address)?;
    let mut eps = vec![0.0f64; ns + 1];
    if sigma > 0.0 {
        for (i, e) in eps.iter_mut().enumerate() { *e = scope.draw(PURPOSE_PERTURB, i as u64, inp.realization).map_err(|e| e.to_string())?.normal_f64(); }
        binomial2(&mut eps);
    }
    let pin = |s: f64| -> f64 { let lt = (taper_w * wmin).max(1e-6); smooth01((s - s0) / lt) * smooth01((s1 - s) / lt) };
    let mut x: Vec<[f64; 2]> = Vec::with_capacity(ns + 1);
    let mut lab: Vec<f64> = Vec::with_capacity(ns + 1);
    for i in 0..=ns {
        let s = svals[i];
        let l = lab_at_s(s);
        let (p, nrm) = carrier_at(l);
        let li = l.round() as usize;
        let off = sigma * inp.w[li] * eps[i] * pin(s) * freg[li];
        x.push([p[0] + off * nrm[0], p[1] + off * nrm[1]]);
        lab.push(l);
    }
    // migration
    let t_end = mp.age_yr as f64;
    let mut t = 0.0f64;
    let mut st = MqSeg { segment: inp.seg, address: inp.address.to_string(), width_min_m: wmin, width_max_m: wmax, ..Default::default() };
    let mut oxb: Vec<(f64, &'static str, Vec<[f64; 2]>, f64)> = Vec::new();
    let mut cells: HashMap<(i64, i64), f32> = HashMap::new();
    let mark = |x: &[[f64; 2]], lab: &[f64], cells: &mut HashMap<(i64, i64), f32>| {
        for (i, p) in x.iter().enumerate() {
            let r = (0.5 + margin) * interp(&inp.w, lab[i]);
            let e = cells.entry(cell_of(*p)).or_insert(0.0);
            if (r as f32) > *e { *e = r as f32; }
        }
    };
    mark(&x, &lab, &mut cells);
    let any_mobile = eb > 0.0 && t_end > 0.0 && (k0..=k1).any(|k| freg[k] > 0.0);
    let mut moved = 0.0f64;
    let mut near: Vec<usize> = lab.iter().map(|&l| l.round() as usize).collect();
    let max_steps: u64 = 2_000_000;
    // bucket grid for cutoffs (head/next linked lists, cleared per step)
    let cut_cell = (chute.max(neck) * wmax).max(ds0);
    let mut head: HashMap<(i64, i64), u32> = HashMap::new();
    let mut next: Vec<u32> = Vec::new();
    while any_mobile && t < t_end && st.steps < max_steps {
        let n = x.len();
        // geometry
        let mut ds = vec![0.0f64; n];
        let mut s = vec![0.0f64; n];
        for i in 1..n { ds[i] = (x[i][0] - x[i - 1][0]).hypot(x[i][1] - x[i - 1][1]).max(1e-9); s[i] = s[i - 1] + ds[i]; }
        let mut c = vec![0.0f64; n];
        for i in 1..n - 1 {
            let a0 = (x[i][1] - x[i - 1][1]).atan2(x[i][0] - x[i - 1][0]);
            let a1 = (x[i + 1][1] - x[i][1]).atan2(x[i + 1][0] - x[i][0]);
            let mut da = a1 - a0;
            while da > std::f64::consts::PI { da -= 2.0 * std::f64::consts::PI; }
            while da < -std::f64::consts::PI { da += 2.0 * std::f64::consts::PI; }
            c[i] = da / (0.5 * (ds[i] + ds[i + 1]));
        }
        binomial2(&mut c);
        // rates
        let (mut num, mut den) = (0.0f64, 0.0f64);
        let mut rate = vec![0.0f64; n];
        let mut maxr = 0.0f64;
        let mut dt_stab = f64::INFINITY;
        let mut dt_cfl = f64::INFINITY;
        for i in 0..n {
            // nearest carrier sample (hill-climb from the cached one): the station's valley position (Eulerian label), which
            // sets its hydraulics, regime and walls (stations drift down-valley with their bends)
            let mut k = near[i].clamp(k0, k1);
            let d2 = |k: usize| { let (dx, dy) = (inp.cx[k][0] - x[i][0], inp.cx[k][1] - x[i][1]); dx * dx + dy * dy };
            loop {
                let cur = d2(k);
                if k > k0 && d2(k - 1) < cur { k -= 1; continue; }
                if k < k1 && d2(k + 1) < cur { k += 1; continue; }
                break;
            }
            near[i] = k;
            lab[i] = k as f64;
            let l = k as f64;
            let w = interp(&inp.w, l);
            let d = interp(&inp.dep, l).max(0.05);
            let cf = G * nman * nman / d.cbrt();
            let alpha = 2.0 * cf / d;
            let r0 = w * w * c[i].clamp(-csat / w, csat / w);
            let dsb = 0.5 * (ds[i] + if i + 1 < n { ds[i + 1] } else { 0.0 });
            let e = if i > 0 { (-alpha * ds[i]).exp() } else { 0.0 };
            num = num * e + r0 * dsb;
            den = den * e + dsb;
            let r1 = omega * r0 + gamma * if den > 0.0 { num / den } else { 0.0 };
            let off = (x[i][0] - inp.cx[k][0]) * inp.cn[k][0] + (x[i][1] - inp.cx[k][1]) * inp.cn[k][1];
            st.belt_half_max_m = st.belt_half_max_m.max(off.abs());
            // outer normal (right of the flow direction) and whether the move goes away from the carrier
            let (a, b) = (x[i.saturating_sub(1)], x[(i + 1).min(n - 1)]);
            let (tx, ty) = (b[0] - a[0], b[1] - a[1]);
            let tl = tx.hypot(ty).max(1e-12);
            let nr = [ty / tl, -tx / tl];
            let away = r1 * (nr[0] * inp.cn[k][0] + nr[1] * inp.cn[k][1]) * off > 0.0;
            let excess = off.abs() + 0.5 * w - inp.hv[k];
            let fw = if away && excess > 0.0 { let u = (excess / (0.5 * w).max(1e-6)).min(1.0); st.wall_station_steps += 1; 1.0 + u * (wall * inp.krock[k] - 1.0) } else { 1.0 };
            let f = freg[k.clamp(k0, k1)] * fw * pin(inp.cs[k]);
            rate[i] = eb * f * r1;
            maxr = maxr.max(rate[i].abs());
            let dsi = if i > 0 { ds[i] } else { ds.get(1).copied().unwrap_or(ds0) };
            if f > 0.0 && eb > 0.0 { dt_stab = dt_stab.min(0.4 * dsi * dsi / (0.17 * std::f64::consts::PI.powi(2) * eb * f * w * w)); }
            if rate[i] != 0.0 { dt_cfl = dt_cfl.min(0.25 * dsi / rate[i].abs()); }
        }
        let mut dt = (mp.dt_max_yr as f64).min(t_end - t);
        let _ = maxr;
        dt = dt.min(dt_cfl).min(dt_stab);
        dt = dt.max(1e-3);
        // move
        let mut step_max = 0.0f64;
        let mut nx = x.clone();
        for i in 1..n - 1 {
            let (a, b) = (x[i - 1], x[i + 1]);
            let (tx, ty) = (b[0] - a[0], b[1] - a[1]);
            let tl = tx.hypot(ty).max(1e-12);
            let m = rate[i] * dt;
            nx[i] = [x[i][0] + m * ty / tl, x[i][1] - m * tx / tl];
            step_max = step_max.max(m.abs());
        }
        x = nx;
        t += dt;
        st.steps += 1;
        moved += step_max;
        // resample (labels interpolated; endpoints kept)
        let mut rx = Vec::with_capacity(x.len() + 8);
        let mut rl = Vec::with_capacity(x.len() + 8);
        let mut rn = Vec::with_capacity(x.len() + 8);
        rx.push(x[0]); rl.push(lab[0]); rn.push(near[0]);
        for i in 1..x.len() {
            let last = *rx.last().unwrap();
            let d = (x[i][0] - last[0]).hypot(x[i][1] - last[1]);
            let dl = dsl(0.5 * (lab[i] + *rl.last().unwrap()));
            if d < 0.6 * dl && i + 1 < x.len() { continue; }
            if d > 1.6 * dl {
                let m = ((d / dl).round() as usize).max(2);
                let ll = *rl.last().unwrap();
                for q in 1..m {
                    let f = q as f64 / m as f64;
                    rx.push([last[0] + f * (x[i][0] - last[0]), last[1] + f * (x[i][1] - last[1])]);
                    rl.push(ll + f * (lab[i] - ll));
                    rn.push(near[i]);
                }
            }
            rx.push(x[i]); rl.push(lab[i]); rn.push(near[i]);
        }
        x = rx; lab = rl; near = rn;
        // cutoffs
        let n = x.len();
        let mut s = vec![0.0f64; n];
        for i in 1..n { s[i] = s[i - 1] + (x[i][0] - x[i - 1][0]).hypot(x[i][1] - x[i - 1][1]); }
        head.clear();
        next.clear();
        next.resize(n, u32::MAX);
        let key = |p: [f64; 2]| ((p[0] / cut_cell).floor() as i64, (p[1] / cut_cell).floor() as i64);
        for i in (0..n).rev() { let kk = key(x[i]); let h = head.entry(kk).or_insert(u32::MAX); next[i] = *h; *h = i as u32; }
        let mut cut: Option<(usize, usize, &'static str, f64)> = None;
        'scan: for i in 0..n {
            let w = interp(&inp.w, lab[i]);
            let (crn, crc) = (neck * w, chute * w);
            let (cx0, cy0) = key(x[i]);
            let mut best: Option<(usize, &'static str, f64)> = None;
            for oy in -1..=1 { for ox in -1..=1 {
                let mut j = *head.get(&(cx0 + ox, cy0 + oy)).unwrap_or(&u32::MAX);
                while j != u32::MAX {
                    let ju = j as usize;
                    j = next[ju];
                    if ju <= i + 2 || s[ju] - s[i] <= 3.0 * crn { continue; }
                    let d = (x[ju][0] - x[i][0]).hypot(x[ju][1] - x[i][1]);
                    let pinned = pin(interp(inp.cs, lab[i])) < 0.5 || pin(interp(inp.cs, lab[ju])) < 0.5;
                    if pinned { continue; }
                    if d < crn {
                        if best.map(|b| b.1 != "neck" || ju > b.0).unwrap_or(true) { best = Some((ju, "neck", d)); }
                    } else if d < crc && (s[ju] - s[i]) >= chute_ratio * d && best.map(|b| b.1 != "neck" && ju > b.0).unwrap_or(true) {
                        let mid = [0.5 * (x[i][0] + x[ju][0]), 0.5 * (x[i][1] + x[ju][1])];
                        if cells.contains_key(&cell_of(mid)) { best = Some((ju, "chute", d)); }
                    }
                }
            } }
            if let Some((j, kind, d)) = best { cut = Some((i, j, kind, d)); break 'scan; }
        }
        if let Some((i, j, kind, d)) = cut {
            oxb.push((t, kind, x[i..=j].to_vec(), d));
            if kind == "neck" { st.neck_cutoffs += 1; } else { st.chute_cutoffs += 1; }
            x.drain(i + 1..j); lab.drain(i + 1..j); near.drain(i + 1..j);
        }
        if moved >= 0.5 * CELL_M { mark(&x, &lab, &mut cells); moved = 0.0; }
    }
    mark(&x, &lab, &mut cells);
    // output labels: each station's valley position (nearest carrier sample), made monotone by a running maximum, so a T1a
    // node maps to the first channel station that reaches its valley position (loops that fold back keep the label)
    for i in 0..x.len() {
        let mut k = near[i].clamp(k0, k1);
        let d2 = |k: usize| { let (dx, dy) = (inp.cx[k][0] - x[i][0], inp.cx[k][1] - x[i][1]); dx * dx + dy * dy };
        loop {
            let cur = d2(k);
            if k > k0 && d2(k - 1) < cur { k -= 1; continue; }
            if k < k1 && d2(k + 1) < cur { k += 1; continue; }
            break;
        }
        near[i] = k;
    }
    lab = near.iter().map(|&k| k as f64).collect();
    lab[0] = lab[0].min(k0 as f64);
    for i in 1..lab.len() { lab[i] = lab[i].max(lab[i - 1]); }
    let nl = lab.len();
    lab[nl - 1] = lab[nl - 1].max(k1 as f64);
    st.years = t;
    st.stations_final = x.len() as u64;
    // band on non-migrating stretches (braided by the regime criterion, submerged)
    let mut band = Vec::new();
    let mut k = k0;
    while k <= k1 {
        let fm = freg[k];
        let r = inp.hband[k] * (1.0 - fm);
        if r > (0.5 + margin) * inp.w[k] { band.push([inp.cx[k][0] as f32, inp.cx[k][1] as f32, r as f32]); }
        k += 8;
    }
    // statistics
    st.carrier_m = s1 - s0;
    st.channel_m = x.windows(2).map(|p| (p[1][0] - p[0][0]).hypot(p[1][1] - p[0][1])).sum();
    st.sinuosity = st.channel_m / st.carrier_m.max(1e-9);
    let mut sl: Vec<f64> = (k0..=k1).map(|k| inp.slope[k]).collect();
    sl.sort_by(|a, b| a.total_cmp(b));
    st.slope_p50 = sl[sl.len() / 2];
    let mut pv: Vec<f64> = (k0..=k1).map(|k| pidx[k]).filter(|v| v.is_finite()).collect();
    pv.sort_by(|a, b| a.total_cmp(b));
    st.regime_index_p50 = if pv.is_empty() { f64::NAN } else { pv[pv.len() / 2] };
    st.migrating_share = (k0..=k1).filter(|&k| freg[k] >= 0.5).count() as f64 / (k1 - k0 + 1) as f64;
    let (lam, amp) = halfwaves(&x, &lab, &inp.w);
    st.halfwaves = lam.len() as u32;
    let pct = |v: &[f64], q: f64| -> f64 { if v.is_empty() { f64::NAN } else { let mut u = v.to_vec(); u.sort_by(|a, b| a.total_cmp(b)); u[((u.len() - 1) as f64 * q).round() as usize] } };
    let cv = |v: &[f64]| -> f64 { if v.len() < 2 { return f64::NAN; } let m = v.iter().sum::<f64>() / v.len() as f64; (v.iter().map(|a| (a - m) * (a - m)).sum::<f64>() / (v.len() - 1) as f64).sqrt() / m.max(1e-12) };
    st.wavelength_widths_p10 = pct(&lam, 0.1); st.wavelength_widths_p50 = pct(&lam, 0.5); st.wavelength_widths_p90 = pct(&lam, 0.9); st.wavelength_cv = cv(&lam);
    st.amplitude_widths_p10 = pct(&amp, 0.1); st.amplitude_widths_p50 = pct(&amp, 0.5); st.amplitude_widths_p90 = pct(&amp, 0.9); st.amplitude_cv = cv(&amp);
    let r1 = |v: f64| (v * 10.0).round() / 10.0;
    let oxbows = oxb.into_iter().map(|(tc, kind, pts, d)| {
        let len = pts.windows(2).map(|p| (p[1][0] - p[0][0]).hypot(p[1][1] - p[0][1])).sum();
        Oxbow { segment: inp.seg, kind, age_yr: r1(t - tc), length_m: r1(len), neck_m: r1(d), pts: pts.iter().step_by(2).chain(pts.last()).map(|p| [r1(p[0]), r1(p[1])]).collect() }
    }).collect();
    Ok(SegOut { seg: inp.seg, x, lab, oxbows, cells, band, stats: st })
}

/// Half-waves of a channel: (wavelength / w, amplitude / w) per half-wave between successive sign changes of the smoothed
/// curvature whose chord exceeds 2 w (w at the half-wave's start).
pub fn halfwaves(x: &[[f64; 2]], lab: &[f64], w: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = x.len();
    if n < 5 { return (Vec::new(), Vec::new()); }
    let mut c = vec![0.0f64; n];
    for i in 1..n - 1 {
        let a0 = (x[i][1] - x[i - 1][1]).atan2(x[i][0] - x[i - 1][0]);
        let a1 = (x[i + 1][1] - x[i][1]).atan2(x[i + 1][0] - x[i][0]);
        let mut da = a1 - a0;
        while da > std::f64::consts::PI { da -= 2.0 * std::f64::consts::PI; }
        while da < -std::f64::consts::PI { da += 2.0 * std::f64::consts::PI; }
        c[i] = da;
    }
    for _ in 0..3 { binomial2(&mut c); }
    let mut inf = Vec::new();
    for i in 1..n { if (c[i - 1] < 0.0) != (c[i] < 0.0) && c[i] != 0.0 { inf.push(i); } }
    let (mut lam, mut amp) = (Vec::new(), Vec::new());
    for p in inf.windows(2) {
        let (a, b) = (x[p[0]], x[p[1]]);
        let wl = interp(w, lab[p[0]]);
        let ch = (b[0] - a[0]).hypot(b[1] - a[1]);
        if ch < 2.0 * wl { continue; }
        let (ux, uy) = ((b[0] - a[0]) / ch, (b[1] - a[1]) / ch);
        let mut am = 0.0f64;
        for q in x.iter().take(p[1] + 1).skip(p[0]) { am = am.max(((q[0] - a[0]) * uy - (q[1] - a[1]) * ux).abs()); }
        lam.push(2.0 * ch / wl);
        amp.push(am / wl);
    }
    (lam, amp)
}

/// The chart-level floodplain extent of all migrated reaches: discs (x, y, r) and the largest distance from a disc's edge to
/// the nearest final channel station (the floor reach pad).
#[derive(Debug, Default)]
pub struct MqEnv { pub discs: Vec<[f32; 3]>, pub pad_m: f64, pub inset: bool }

pub fn envelope_refs(outs: &[&SegOut], inset: bool) -> MqEnv {
    let mut all: std::collections::BTreeMap<(i64, i64), f32> = std::collections::BTreeMap::new();
    for o in outs { for (k, &r) in &o.cells { let e = all.entry(*k).or_insert(0.0); if r > *e { *e = r; } } }
    let mut discs: Vec<[f32; 3]> = all.iter().map(|(k, &r)| [((k.0 as f64 + 0.5) * CELL_M) as f32, ((k.1 as f64 + 0.5) * CELL_M) as f32, r]).collect();
    for o in outs { discs.extend_from_slice(&o.band); }
    // pad: distance from every disc to the nearest final channel station, plus its radius
    let b = 256.0f64;
    let mut grid: HashMap<(i64, i64), Vec<[f64; 2]>> = HashMap::new();
    for o in outs { for p in &o.x { grid.entry(((p[0] / b).floor() as i64, (p[1] / b).floor() as i64)).or_default().push(*p); } }
    let mut pad = 0.0f64;
    if !grid.is_empty() {
        for d in &discs {
            let (px, py) = (d[0] as f64, d[1] as f64);
            let (cx, cy) = ((px / b).floor() as i64, (py / b).floor() as i64);
            let mut best = f64::INFINITY;
            let mut ring = 0i64;
            while ring < 64 {
                for oy in -ring..=ring { for ox in -ring..=ring {
                    if oy.abs() != ring && ox.abs() != ring { continue; }
                    if let Some(v) = grid.get(&(cx + ox, cy + oy)) { for q in v { best = best.min((q[0] - px).hypot(q[1] - py)); } }
                } }
                if best.is_finite() && best <= ring as f64 * b { break; }
                ring += 1;
            }
            if best.is_finite() { pad = pad.max(best + d[2] as f64); }
        }
    }
    MqEnv { discs, pad_m: pad, inset }
}

/// Run every T1a reach of the chart (pure function of the chart record: T1a curve records, T1a nodes and levels, the T1b
/// network beds and the volcano corridor cap; never the window). Returns the migrated reaches, the floodplain and statistics.
pub fn run_chart(tabs: &std::collections::BTreeMap<u32, crate::t1b_cl::SegTab>, node_k: &std::collections::BTreeMap<usize, (u32, usize)>,
    skel: &crate::t1b_view::SkeletonView, beds: &crate::t1b_rivernet::Beds, t0: &crate::t0::T0, p: &crate::t1b_rivernet::RiverNetParamsV1,
    confine: &dyn Fn([f64; 2]) -> f64) -> Result<(std::collections::BTreeMap<u32, SegOut>, MqEnv, MqStats), String> {
    let cp = crate::t1b_cl::cl(p).ok_or("cl_mq needs cl=1")?;
    let mp = cp.mq.as_ref().ok_or("cl_mq params")?;
    let vfd = crate::t1b_vf::VfParamsV1::off();
    let vfp = p.vf.as_ref().unwrap_or(&vfd);
    let cw = p.floor_width_coeff_ppm as f64 / 1e6;
    let hmax = p.floor_width_max_mm as f64 / 2000.0;
    let t1a_area = |n: usize| skel.discharge_m3_s[n] as f64 * crate::SECONDS_PER_YEAR / t0.runoff_m_per_yr;
    let qof = |a: f64| a.max(0.0) * vfp.runoff_ref_mm_per_yr as f64 / 1000.0 / crate::SECONDS_PER_YEAR;
    let ok_r = |r: u32| r != crate::NONE && r != crate::OFF_CHART;
    // link levels exactly as the T1a link loop of curves_with_stats (la at the node, lb at its receiver)
    let levels = |n: usize, r: usize| -> (f64, f64) {
        let (mut la, mut lb) = (skel.trunk_level_m[n] as f64, skel.trunk_level_m[r] as f64);
        if let Some(tl) = &beds.t1a_level {
            if crate::t1b_ls::t1a_mode(p) != 0 { la = tl[n]; lb = if tl[r].is_finite() { tl[r].min(la) } else { la }; }
            else { la = la.min(tl[n]); lb = lb.min(tl[r]).min(la); }
        }
        (la, lb)
    };
    let mut per: std::collections::BTreeMap<u32, Vec<(usize, usize)>> = std::collections::BTreeMap::new();
    for (&n, &(sg, k)) in node_k { per.entry(sg).or_default().push((k, n)); }
    let mut outs = std::collections::BTreeMap::new();
    let mut stats = MqStats::default();
    for (sg, mut v) in per {
        v.sort();
        let tab = &tabs[&sg];
        let seg = &skel.t1a_segments[sg as usize];
        let seed = seg.seed.ok_or("T1a header has no recipe seed_scope.seed (world seed) for the world-keyed draws")?;
        if seg.address.is_empty() { return Err(format!("T1a segment {sg} has no address")); }
        let nk = tab.x.len();
        if nk < 3 { continue; }
        // per-node hydraulics and link slope
        let mut nd: Vec<(usize, f64, f64, f64, f64, f64, f64, f64)> = Vec::new(); // (k, w, D, q, h, slope, krock, wet)
        for &(k, n) in &v {
            let a = t1a_area(n);
            let (w, dd) = vfp.hydraulic(a);
            let h = (0.5 * cw * a.max(0.0).sqrt()).min(hmax);
            let r = skel.receiver[n];
            let mut sl = f64::NAN;
            if ok_r(r) {
                if let Some(&(sr, kr)) = node_k.get(&(r as usize)) {
                    if sr == sg && kr > k { let (la, lb) = levels(n, r as usize); sl = (la - lb).max(0.0) / (tab.s[kr] - tab.s[k]).max(1.0); }
                }
            }
            let kr = beds.mq_rock.as_ref().map(|m| m[n]).unwrap_or(1.0);
            let lvl = if ok_r(r) { levels(n, r as usize).0 } else { skel.trunk_level_m[n] as f64 };
            let wet = if lvl <= skel.sea_level_m + 0.5 || skel.kind[n] == 3 { 1.0 } else { 0.0 };
            nd.push((k, w, dd, qof(a), h, sl, kr, wet));
        }
        // fill missing slopes from the nearest node with one (downstream first, then upstream)
        for i in 0..nd.len() {
            if nd[i].5.is_finite() { continue; }
            let prev = (0..i).rev().find(|&j| nd[j].5.is_finite()).map(|j| nd[j].5);
            let nxt = (i + 1..nd.len()).find(|&j| nd[j].5.is_finite()).map(|j| nd[j].5);
            nd[i].5 = prev.or(nxt).unwrap_or(f64::INFINITY);
        }
        let mut arr = vec![vec![0.0f64; nk]; 7];
        let mut ni = 0usize;
        for k in 0..nk {
            while ni + 1 < nd.len() && nd[ni + 1].0 <= k { ni += 1; }
            let e = nd[ni];
            for (j, val) in [e.1, e.2, e.3, e.4, e.5, e.6, e.7].iter().enumerate() { arr[j][k] = *val; }
        }
        // moving average over +-500 m (125 carrier samples)
        let wsm = 125usize;
        for a in arr.iter_mut() {
            if a.iter().any(|v| !v.is_finite()) { continue; }
            let mut c = vec![0.0f64; nk + 1];
            for k in 0..nk { c[k + 1] = c[k] + a[k]; }
            *a = (0..nk).map(|k| { let (lo, hi) = (k.saturating_sub(wsm), (k + wsm).min(nk - 1)); (c[hi + 1] - c[lo]) / (hi + 1 - lo) as f64 }).collect();
        }
        let hv: Vec<f64> = (0..nk).map(|k| arr[3][k].min(confine(tab.x[k]))).collect();
        let k0 = nd[0].0;
        let k1 = (nd[nd.len() - 1].0 + 256).min(nk - 1);
        if k1 <= k0 + 2 { continue; }
        let inp = SegIn { seg: sg, address: &seg.address, seed, realization: p.realization, cx: &tab.x, cs: &tab.s, cn: &tab.nrm,
            w: arr[0].clone(), dep: arr[1].clone(), hv, q: arr[2].clone(), slope: arr[4].clone(), krock: arr[5].clone(), hband: arr[3].clone(), wet: arr[6].clone(), k0, k1 };
        let o = run_segment(&inp, mp)?;
        stats.segments.push(o.stats.clone());
        stats.oxbows.extend(o.oxbows.iter().cloned());
        outs.insert(sg, o);
    }
    let list: Vec<&SegOut> = outs.values().collect();
    let env = envelope_refs(&list, mp.inset == 1);
    stats.env_discs = env.discs.len() as u64;
    stats.pad_m = (env.pad_m * 10.0).round() / 10.0;
    Ok((outs, env, stats))
}

/// The floodplain extent on a window lattice: posts inside any envelope disc; and the floor reach pad.
pub fn env_ext(env: &MqEnv, lat: &crate::t1b_rivernet::EvalLat, dx: f64) -> (Vec<bool>, f64) {
    let mut ext = vec![false; lat.nodes()];
    let (gxa, gya) = (lat.gx0, lat.gy0);
    let (gxb, gyb) = (lat.gx0 + lat.n as i64 - 1, lat.gy0 + lat.n as i64 - 1);
    for d in &env.discs {
        let (x, y, r) = (d[0] as f64, d[1] as f64, d[2] as f64);
        let (x0, x1) = (((x - r) / dx).ceil() as i64, ((x + r) / dx).floor() as i64);
        let (y0, y1) = (((y - r) / dx).ceil() as i64, ((y + r) / dx).floor() as i64);
        if x1 < gxa || x0 > gxb || y1 < gya || y0 > gyb { continue; }
        for gy in y0.max(gya)..=y1.min(gyb) { for gx in x0.max(gxa)..=x1.min(gxb) {
            if (gx as f64 * dx - x).hypot(gy as f64 * dx - y) > r { continue; }
            if let Some(i) = lat.idx(gx, gy) { ext[i] = true; }
        } }
    }
    (ext, env.pad_m)
}

/// Window statistics of the physical planform (published under prep.net.cl.mq).
#[derive(Clone, Debug, Default, Serialize)]
pub struct MqStats {
    pub segments: Vec<MqSeg>,
    pub oxbows: Vec<Oxbow>,
    pub env_discs: u64,
    pub pad_m: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn straight(len: f64) -> (Vec<[f64; 2]>, Vec<f64>, Vec<[f64; 2]>) {
        let n = (len / 4.0) as usize;
        ((0..=n).map(|k| [k as f64 * 4.0, 0.0]).collect(), (0..=n).map(|k| k as f64 * 4.0).collect(), vec![[0.0, 1.0]; n + 1])
    }
    fn input<'a>(cx: &'a [[f64; 2]], cs: &'a [f64], cn: &'a [[f64; 2]], w: f64, slope: f64, hv: f64) -> SegIn<'a> {
        let n = cx.len();
        let q = (w / 3.5).powi(2);
        SegIn { seg: 0, address: "seg/1-2", seed: 42, realization: 0, cx, cs, cn, w: vec![w; n], dep: vec![0.27 * q.powf(0.4); n], hv: vec![hv; n],
            q: vec![q; n], slope: vec![slope; n], krock: vec![1.0; n], hband: vec![hv.min(600.0); n], wet: vec![0.0; n], k0: 0, k1: n - 1 }
    }

    #[test]
    fn no_forcing_no_migration() {
        // zero perturbation on a straight carrier: zero curvature everywhere, the channel stays on the carrier
        let (cx, cs, cn) = straight(20_000.0); let w = 35.0;
        let mut mp = MqParamsV1::v1();
        mp.perturb_widths_permille = 0;
        mp.age_yr = 2000;
        let o = run_segment(&input(&cx, &cs, &cn, w, 1e-5, 1e9), &mp).unwrap();
        assert!(o.x.iter().all(|p| p[1].abs() < 1e-9), "moved");
        assert_eq!(o.stats.neck_cutoffs + o.stats.chute_cutoffs, 0);
        // zero bank rate: the perturbed initial state is kept (no steps)
        let mut mp0 = MqParamsV1::v1();
        mp0.bank_rate_ppm_per_yr = 0;
        let o0 = run_segment(&input(&cx, &cs, &cn, w, 1e-5, 1e9), &mp0).unwrap();
        assert_eq!(o0.stats.steps, 0);
    }

    #[test]
    fn steep_reach_is_not_migrated() {
        // slope 0.01 at Q = 100 m3/s: Leopold-Wolman index >> 2, regime factor 0
        let (cx, cs, cn) = straight(10_000.0); let w = 35.0;
        let mp = MqParamsV1::v1();
        let (f, p) = mp.regime_factor(0.01, 100.0);
        assert!(f == 0.0 && p > 10.0, "{f} {p}");
        let o = run_segment(&input(&cx, &cs, &cn, w, 0.01, 1e9), &mp).unwrap();
        assert_eq!(o.stats.steps, 0);
        assert!(o.x.iter().all(|p| p[1].abs() < 1e-9));
        assert!(!o.band.is_empty());
    }

    #[test]
    fn gentle_reach_meanders_cuts_off_and_is_deterministic() {
        let (cx, cs, cn) = straight(30_000.0); let w = 35.0;
        let mut mp = MqParamsV1::v1();
        mp.age_yr = 6000;
        let inp = input(&cx, &cs, &cn, w, 1e-4, 1e9);
        let o = run_segment(&inp, &mp).unwrap();
        let o2 = run_segment(&inp, &mp).unwrap();
        assert_eq!(o.x, o2.x, "not deterministic");
        assert!(o.stats.sinuosity > 1.2, "sinuosity {}", o.stats.sinuosity);
        assert!(o.lab.windows(2).all(|l| l[1] >= l[0]), "labels not monotone");
        assert!(o.stats.halfwaves > 5);
        eprintln!("{:?}", o.stats);
    }

    #[test]
    fn walls_confine_the_belt() {
        let (cx, cs, cn) = straight(30_000.0); let w = 35.0;
        let mut mp = MqParamsV1::v1();
        mp.age_yr = 6000;
        let free = run_segment(&input(&cx, &cs, &cn, w, 1e-4, 1e9), &mp).unwrap();
        let conf = run_segment(&input(&cx, &cs, &cn, w, 1e-4, 100.0), &mp).unwrap();
        assert!(conf.stats.belt_half_max_m < free.stats.belt_half_max_m, "{} vs {}", conf.stats.belt_half_max_m, free.stats.belt_half_max_m);
        // walls erode at wall_factor (0.1) of the alluvial rate: the belt leaks past h_v by lateral planation, slowly
        assert!(conf.stats.belt_half_max_m < 0.2 * free.stats.belt_half_max_m && conf.stats.belt_half_max_m < 100.0 + 3.0 * w, "{} vs {}", conf.stats.belt_half_max_m, free.stats.belt_half_max_m);
        assert!(conf.stats.sinuosity < free.stats.sinuosity, "{} vs {}", conf.stats.sinuosity, free.stats.sinuosity);
        eprintln!("confined: belt half {:.0} m sinuosity {:.2}; free: {:.0} m {:.2}", conf.stats.belt_half_max_m, conf.stats.sinuosity, free.stats.belt_half_max_m, free.stats.sinuosity);
    }

    #[test]
    #[ignore]
    fn dt_and_station_sensitivity() {
        // measurement (run with --ignored --nocapture): statistics under dt_max 10 / 2 / 0.5 yr and station 1.0 / 0.5 widths
        let (cx, cs, cn) = straight(30_000.0); let w = 35.0;
        for hv in [1e9, 100.0] {
            for (dt, st) in [(10, 1000), (2, 1000), (1, 1000), (10, 500), (2, 500), (10, 250)] {
                let mut mp = MqParamsV1::v1();
                mp.age_yr = 6000; mp.dt_max_yr = dt; mp.station_widths_permille = st;
                let o = run_segment(&input(&cx, &cs, &cn, w, 1e-4, hv), &mp).unwrap();
                eprintln!("hv {hv} dt_max {dt} station {st}: steps {} sinuosity {:.2} belt half {:.0} m lambda/w p50 {:.1} cv {:.2} amp/w p50 {:.1} neck {} chute {}",
                    o.stats.steps, o.stats.sinuosity, o.stats.belt_half_max_m, o.stats.wavelength_widths_p50, o.stats.wavelength_cv, o.stats.amplitude_widths_p50, o.stats.neck_cutoffs, o.stats.chute_cutoffs);
            }
        }
    }

    #[test]
    #[ignore]
    fn dump_synthetic_planforms() {
        // measurement (run with --ignored): final channel, oxbows and envelope cells of a 30 km straight reach, slope 1e-4,
        // w 35 m, unconfined and h_v = 300 m, 11.7 ka; written as JSON to env VZ_MQ_SYNTH (no file without it)
        let Ok(path) = std::env::var("VZ_MQ_SYNTH") else { return };
        let (cx, cs, cn) = straight(30_000.0); let w = 35.0;
        let mut out = Vec::new();
        for hv in [1e9, 300.0] {
            let o = run_segment(&input(&cx, &cs, &cn, w, 1e-4, hv), &MqParamsV1::v1()).unwrap();
            let cells: Vec<[i64; 2]> = { let mut v: Vec<[i64; 2]> = o.cells.keys().map(|k| [k.0, k.1]).collect(); v.sort(); v };
            out.push(serde_json::json!({"hv": hv, "x": o.x, "oxbows": o.oxbows, "cells16": cells, "stats": o.stats}));
        }
        std::fs::write(path, serde_json::to_vec(&out).unwrap()).unwrap();
    }

    #[test]
    fn world_keyed_scope_ignores_chart() {
        let a = reach_scope(42, "seg/1-2").unwrap().draw(PURPOSE_PERTURB, 3, 0).unwrap().normal_f64();
        let b = reach_scope(42, "seg/1-2").unwrap().draw(PURPOSE_PERTURB, 3, 0).unwrap().normal_f64();
        let c = reach_scope(42, "seg/1-3").unwrap().draw(PURPOSE_PERTURB, 3, 0).unwrap().normal_f64();
        assert_eq!(a, b);
        assert_ne!(a, c);
    }
}
