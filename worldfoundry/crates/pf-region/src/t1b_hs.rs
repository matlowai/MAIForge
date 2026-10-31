//! Lane HS-PASS (2026-10-21): opt-in hillslope-only transient on the T1b RIVERNET window surface (design review consult
//! `unpublished-design-record review-planet/CONSULT_T1B_LOOK.md`, option (a)). SimTime label `t1b-hillslope-transient-not-history`.
//! LABEL: EXPLORATORY, every constant a project assumption awaiting review (sources in the rule strings).
//!
//! After the ct remap, the coast, the volcano and the vf stages (the finished `evaluate` surface), before publication, the
//! window integrates three laws for a declared bounded time T = N dt on the FREE posts of the evaluation lattice:
//! nonlinear creep (Roering), threshold landsliding with talus deposition, and MFD stream-power incision above a channel
//! initiation threshold. Channel, valley-floor (vf band), marine, lake, guide and T1a bank posts are pinned. Heights are
//! integer mm, all transfers integer mm summed in i64, Jacobi order (every update of a sub-step reads the sub-step-start
//! surface; the evaluation-lattice ring is frozen during the transient and acts as a sink, never as a flood seed), so the result does not depend on any iteration order. Depressions are filled by a priority flood seeded ONLY
//! from pinned posts (seamless-drainage Design 1), and every routing emits a dependency certificate. Absent params
//! (`params.hs` None) = byte-identical.

use crate::t1b::{CHANNEL, FREE};
use crate::t1b_litho::ExposedLithology;
use crate::t1b_rivernet::{EvalLat, RockPosts};
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};

pub const LABEL: &str = "EXPLORATORY lane HS-PASS 2026-10-21 (design review option (a)); every constant a project assumption awaiting review";
pub const SIMTIME: &str = "t1b-hillslope-transient-not-history: a declared bounded transient of N steps of dt on the constructed surface (base levels pinned); not a landscape history";
pub const CREEP_RULE: &str = "roering_creep_8n_v1: for every 8-neighbour pair (i, j) with at least one free post, transfer t = D dt / dx^2 x c x dh / (1 - min(S / Sc, cap)^2) mm from the higher to the lower post (c = 2/3 axis, 1/6 diagonal: the isotropic 9-point Laplacian; S = |dh| / d_ij; Sc = critical slope of the higher post: lithology table class of the rock exposed at its height, or the talus tangent where it carries talus), limited to |dh| / 16 (no overshoot in the Jacobi step), rounded to integer mm; mass leaving to a pinned post is exported; a pinned post is a sink only (no transfer out of a pinned post, so perched construction channel posts do not feed their free neighbours; creep_supply stays 0). Roering, Kirchner and Dietrich 1999 WRR 35:853 (form, recalled). D uniform for every class (the lithology table has no diffusivity column)";
pub const SLIDE_RULE: &str = "threshold_landslide_talus_v1: after the creep and incision sub-step, every free post whose height exceeds a lower 8-neighbour by more than Sc x d_ij (Sc of its exposed rock, or the talus tangent where it carries talus) loses ceil(e / 2) mm (e = the largest excess; Jacobi halving); the debris runs down the D8 steepest-descent path of the sub-step-start surface and is deposited at the first post whose steepest slope is <= the talus tangent (0.70, 35 deg, angle of repose; recalled) or after runout_max posts, or exported on reaching a pinned post; deposited debris is talus (thickness tracked; eroded first; critical slope = talus tangent; erodibility of class 10 unconsolidated)";
pub const INCISION_RULE: &str = "mfd_stream_power_threshold_v1: every route_every steps the effective area A (m2; post area x runoff weight R(x)/r_ref under the climate key, else 1) is accumulated by MFD (Freeman 1991 Comput. Geosci. 17:413, weights S^p, p = params.mfd_p_permille / 1000; v0 1.1) on the current surface filled by a priority flood seeded only from pinned posts (flats pass to the flood parent); each step a free post with steepest slope S (raw surface) and A S^2 > A_1 (channel initiation, form of Montgomery and Dietrich 1988 Nature 336:232, recalled) lowers by (K_T0 x k_mult(exposed rock) x A^m S^n - threshold) dt, capped at half the drop to its steepest receiver (K, m, n = the T0 laws; exposed rock looked up at route time, refreshed when a post moved >= 0.5 m since its last lookup; threshold form Tucker 2004 ESPL 29:185, recalled); eroded mass is exported (detachment-limited, as T0)";
pub const FILL_RULE: &str = "pinned_priority_flood_v1: at the end, floor drainage is kept open first (a pinned free-coded post - valley floor, bank - that drained before the pass but not after gets its pre-pass drainage path (flood parents from the true outlets, fixed != FREE, on the pre-pass surface) lowered to its level where the transient raised it, pits in (z, index) order; volume ledgered as carve), then an 8-neighbour priority flood seeded only from pinned posts (key (z mm, lattice index) = (z, gy, gx)) raises every free post below its spill level to it; no lattice-ring seed (seamless-drainage Design 1); certificate: for every published post the farthest post (Chebyshev posts outside the published box) on its flood-parent path and in its upslope MFD area; it fails if that reaches the evaluation lattice ring";
pub const PIN_RULE: &str = "pinned posts: fixed != FREE (T1a channel, marine, lake, network channel/guide), valley-floor band posts (Surface.floor, incl. the vf_floor band) not rebuilt by a later stage (over code 0, or a fan deposit on the floor (OVER_FAN: lane CO bayhead fan / vf delta surface / tributary fan, kept as constructed); a floor buried by an edifice or a beach / wave cut is free), and T1a bank posts (free posts 8-adjacent to a T1a channel post); everything else (hillslopes, edifice, fans, beaches) is free";
pub const UPLIFT_RULE: &str = "uplift (off unless hs_uplift=1): free posts rise by the window's E field (the steady-state erosion rate of the construction) x dt, relative to the pinned base levels";
pub const VOLCANO_RULE: &str = "hs_volcano: 1 = lane CF radial_gullies_young_v1 gully list emptied in the window (gully_count 0; draws unchanged, so every other keyed draw is the same); 2 = also the record_age_degradation_v1 gully directions; the pass dissects the cone instead";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct HsParamsV1 {
    /// 1 = run the transient (key hs=1); 0 = only the volcano switch (hs_volcano without hs).
    pub on: u32,
    /// N (key hs_n).
    pub steps: u32,
    /// dt (key hs_dt_yr).
    pub dt_yr: i64,
    /// Creep diffusivity (key hs_d, m2/My).
    pub d_m2_per_my: i64,
    /// Multiplier on the T0 stream-power K (key hs_k_permille).
    pub k_permille: i64,
    /// Channel initiation A_1 (key hs_a1_m2).
    pub a1_m2: i64,
    /// Stream-power erosion threshold (key hs_th, m/My).
    pub threshold_m_per_my: i64,
    /// Freeman p (key hs_p_permille).
    pub mfd_p_permille: i64,
    /// Talus tangent (key hs_talus_permille).
    pub talus_tan_permille: i64,
    /// Cap of S / Sc in the Roering denominator (permille).
    pub sc_cap_permille: i64,
    /// MFD area recomputed every this many steps (key hs_route_every).
    pub route_every: u32,
    /// Uplift by E (key hs_uplift).
    pub uplift: u32,
    /// Debris runout limit (posts, key hs_runout).
    pub runout_max_posts: u32,
    /// Volcano gully switch (key hs_volcano).
    pub volcano: u32,
    pub creep_rule: String,
    pub slide_rule: String,
    pub incision_rule: String,
    pub fill_rule: String,
    pub pin_rule: String,
    pub uplift_rule: String,
    pub volcano_rule: String,
    pub label: String,
    /// E5: a separate contact-limited, detachment-limited bare-rock transient.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub contrast: Option<crate::t1b_ct::ContrastParams>,
}

impl HsParamsV1 {
    pub fn v0() -> Self {
        HsParamsV1 {
            on: 1, steps: 200, dt_yr: 750, d_m2_per_my: 10_000, k_permille: 1000, a1_m2: 15_000, threshold_m_per_my: 50,
            mfd_p_permille: 1100, talus_tan_permille: 700, sc_cap_permille: 950, route_every: 8, uplift: 0, runout_max_posts: 64, volcano: 0,
            creep_rule: CREEP_RULE.into(), slide_rule: SLIDE_RULE.into(), incision_rule: INCISION_RULE.into(), fill_rule: FILL_RULE.into(),
            pin_rule: PIN_RULE.into(), uplift_rule: UPLIFT_RULE.into(), volcano_rule: VOLCANO_RULE.into(), label: LABEL.into(), contrast: None,
        }
    }
    pub fn validate(&self) -> Result<(), String> {
        if let Some(c) = &self.contrast { c.validate()?; }
        if self.on > 1 || self.uplift > 1 || self.volcano > 2 { return Err("hs: on 0|1, uplift 0|1, volcano 0|1|2".into()); }
        if self.on == 1 && (self.steps == 0 || self.dt_yr <= 0 || self.route_every == 0) { return Err("hs: steps, dt and route_every must be positive".into()); }
        if self.d_m2_per_my < 0 || self.k_permille < 0 || self.a1_m2 < 0 || self.threshold_m_per_my < 0 { return Err("hs: rates must be >= 0".into()); }
        if self.mfd_p_permille <= 0 || self.talus_tan_permille <= 0 || self.sc_cap_permille <= 0 || self.sc_cap_permille >= 1000 { return Err("hs: p, talus tangent > 0; sc cap in (0, 1000)".into()); }
        Ok(())
    }
    pub fn runs(&self) -> bool { self.on == 1 }
}

/// What the pass reads.
pub struct HsIn<'a> {
    pub lat: EvalLat,
    pub fixed: &'a [u8],
    pub floor: &'a [bool],
    /// Construction over-code per post (t1b_coast::OVER_*): a floor post later rebuilt by the coast / volcano stages is not pinned.
    pub over: &'a [u8],
    pub rock: &'a RockPosts<'a>,
    /// R / r_ref per post (empty: 1).
    pub rw: &'a [f32],
    /// E per post (m/My), used only with uplift.
    pub e: &'a [f32],
    pub k_t0_per_my: f64,
    pub m: f64,
    pub n_exp: f64,
    /// Lithology K multiplier and critical slope of class 10 (talus).
    pub talus_k_mult: f64,
    /// Published box and core box, local lattice indices (x0, x1, y0, y1), inclusive.
    pub pub_box: (usize, usize, usize, usize),
    pub core_box: (usize, usize, usize, usize),
    /// Hat blend per post (for the F6 shift) and the per-cell half tolerance (None = not a land cell).
    pub hat_c: &'a [[u32; 4]],
    pub hat_w: &'a [[f32; 4]],
    pub half_tol: &'a dyn Fn(u32) -> Option<f64>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct HsStats {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub contrast: Option<crate::t1b_ct::ContrastStats>,
    pub steps: u32,
    pub dt_yr: i64,
    pub t_yr: i64,
    pub routes: u32,
    pub free_posts: u64,
    pub pinned_posts: u64,
    pub bank_posts: u64,
    pub posts_changed: u64,
    pub published_posts_changed: u64,
    pub max_lowering_m: f64,
    pub max_raise_m: f64,
    /// Global (gx, gy) of the largest raise and lowering, and Chebyshev posts outside the published box.
    pub max_raise_at: (i64, i64, u32),
    pub max_lowering_at: (i64, i64, u32),
    pub incision_m3: f64,
    pub incision_posts_last_step: u64,
    pub creep_export_m3: f64,
    pub creep_supply_m3: f64,
    pub slide_moved_m3: f64,
    pub slide_export_m3: f64,
    pub slide_deposit_m3: f64,
    pub slide_posts_last_step: u64,
    pub runout_capped: u64,
    pub uplift_m3: f64,
    pub fill_m3: f64,
    /// Floor-drainage carve (see run): posts lowered and volume.
    pub carve_posts: u64,
    pub carve_m3: f64,
    pub fill_posts: u64,
    /// Integer ledger: sum of free-post changes minus the ledgered sources and sinks (must be 0).
    pub ledger_residual_mm: i64,
    pub talus_posts_end: u64,
    /// Certificate (posts): farthest dependency outside the published box over every routing and the final fill.
    pub cert_route_up_posts: u32,
    pub cert_route_down_posts: u32,
    pub cert_fill_down_posts: u32,
    pub cert_ring_distance_posts: u32,
    pub cert_core_route_up_posts: u32,
    /// Local stencil cone (posts): N steps x (1 creep + runout) is a bound only; the routing terms are measured.
    pub stencil_cone_posts: u64,
    pub certificate_pass: bool,
    /// F6: per T0 cell hat-mean shift (m) over the evaluation lattice; cells over half tolerance.
    pub f6_max_abs_shift_m: f64,
    /// (cell, shift m, half tolerance m, hat mass of the cell inside the evaluation lattice in posts)
    pub f6_cells_over_half_tol: Vec<(u32, f64, f64, f64)>,
    pub f6_extra_e_update: String,
    /// Wall-clock timings (s). Not serialized: these stats are digested with the window header, and
    /// timing made hs=1 reruns differ (lane SEAM-CERT risk, fixed by lane CHANNEL-CENTRELINE). Printed to stdout instead.
    #[serde(skip)]
    pub seconds: f64,
    #[serde(skip)]
    pub seconds_route: f64,
    /// Sub-step timings (s): steepest-descent passes, creep + incision, landslides.
    #[serde(skip)]
    pub seconds_steep: f64,
    #[serde(skip)]
    pub seconds_creep: f64,
    #[serde(skip)]
    pub seconds_slide: f64,
}

const OFF: i64 = 1 << 31;
const NONE: u32 = u32::MAX;
const D8: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];

#[inline]
fn key(h: i64, i: usize) -> Reverse<u64> { Reverse((((h + OFF) as u64) << 32) | i as u64) }

/// Priority flood from pinned posts only: (pop order, filled heights, flood parent). Err if no pinned post.
fn flood(n: usize, h: &[i64], pinned: &[bool]) -> Result<(Vec<u32>, Vec<i64>, Vec<u32>), String> {
    let nn = n * n;
    let mut hf = h.to_vec();
    let mut parent = vec![NONE; nn];
    let mut seen = vec![false; nn];
    let mut heap = BinaryHeap::with_capacity(nn / 4);
    for i in 0..nn { if pinned[i] { seen[i] = true; heap.push(key(h[i], i)); } }
    if heap.is_empty() { return Err("hs: no pinned post in the evaluation lattice (the fill is seeded only from pinned posts)".into()); }
    let mut order = Vec::with_capacity(nn);
    while let Some(Reverse(k)) = heap.pop() {
        let i = (k & 0xFFFF_FFFF) as usize;
        let v = hf[i];
        order.push(i as u32);
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        for (ox, oy) in D8 {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let j = yy as usize * n + xx as usize;
            if seen[j] { continue; }
            seen[j] = true;
            if hf[j] < v { hf[j] = v; }
            parent[j] = i as u32;
            heap.push(key(hf[j], j));
        }
    }
    Ok((order, hf, parent))
}

#[inline]
fn dist_out(x: usize, y: usize, b: (usize, usize, usize, usize)) -> u32 {
    let dx = if x < b.0 { b.0 - x } else if x > b.1 { x - b.1 } else { 0 };
    let dy = if y < b.2 { b.2 - y } else if y > b.3 { y - b.3 } else { 0 };
    dx.max(dy) as u32
}

/// Run the transient on z (m, evaluation lattice); pinned posts keep their exact f64 values.
pub fn run(inp: &HsIn, p: &HsParamsV1, z: &mut [f64]) -> Result<HsStats, String> {
    if let Some(c) = &p.contrast { return crate::t1b_ct::run_contrast(inp, c, z); }
    let t_all = std::time::Instant::now();
    let lat = inp.lat;
    let (n, nn, dx) = (lat.n, lat.nodes(), lat.dx());
    let mut st = HsStats { steps: p.steps, dt_yr: p.dt_yr, t_yr: p.steps as i64 * p.dt_yr, ..Default::default() };
    let qmm = |v: f64| -> i64 { (v * 1000.0 + 0.5).floor() as i64 };
    let h0: Vec<i64> = z.iter().map(|&v| qmm(v)).collect();
    let mut h = h0.clone();
    // Pinned set.
    let mut pinned = vec![false; nn];
    for i in 0..nn {
        if inp.fixed[i] != FREE || (inp.floor[i] && (inp.over[i] == 0 || inp.over[i] == crate::t1b_coast::OVER_FAN)) { pinned[i] = true; continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        let bank = D8.iter().any(|(ox, oy)| { let (xx, yy) = (x + ox, y + oy); xx >= 0 && yy >= 0 && xx < n as i64 && yy < n as i64 && inp.fixed[yy as usize * n + xx as usize] == CHANNEL });
        if bank { pinned[i] = true; st.bank_posts += 1; }
    }
    st.pinned_posts = pinned.iter().filter(|&&b| b).count() as u64;
    st.free_posts = nn as u64 - st.pinned_posts;
    let ring = (inp.pub_box.0).min(inp.pub_box.2).min(n - 1 - inp.pub_box.1).min(n - 1 - inp.pub_box.3) as u32;
    st.cert_ring_distance_posts = ring;
    let a_post = dx * dx;
    let dts = p.dt_yr as f64;
    let d_yr = p.d_m2_per_my as f64 / 1e6;
    let creep_coef = d_yr * dts / (dx * dx);
    let k_yr = inp.k_t0_per_my * p.k_permille as f64 / 1000.0 / 1e6; // per yr
    let th_yr = p.threshold_m_per_my as f64 / 1e6;
    let a1 = p.a1_m2 as f64;
    let pexp = p.mfd_p_permille as f64 / 1000.0;
    let talus_tan = p.talus_tan_permille as f64 / 1000.0;
    let cap = p.sc_cap_permille as f64 / 1000.0;
    let dist = [dx, dx * std::f64::consts::SQRT_2];
    let dmm = |k: usize| if k % 2 == 0 { 0 } else { 1 };
    // Neighbour k of post i (None off the lattice): offsets without bounds checks away from the ring.
    let offs: [isize; 8] = D8.map(|(ox, oy)| oy as isize * n as isize + ox as isize);
    let on_ring: Vec<bool> = (0..nn).map(|i| { let (x, y) = (i % n, i / n); x == 0 || y == 0 || x == n - 1 || y == n - 1 }).collect();
    let nbr = |i: usize, k: usize| -> Option<usize> {
        if on_ring[i] {
            let (x, y) = (i % n, i / n);
            let (xx, yy) = (x as i64 + D8[k].0, y as i64 + D8[k].1);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { return None; }
        }
        Some((i as isize + offs[k]) as usize)
    };
    // The evaluation-lattice ring is frozen during the transient (no update; a sink for creep and debris, ledgered as export) but
    // is NOT a flood seed: without it debris piles up on ring posts whose receivers lie outside the lattice (seen: +1.1 km).
    let sink: Vec<bool> = (0..nn).map(|i| pinned[i] || on_ring[i]).collect();
    let free_idx: Vec<usize> = (0..nn).filter(|&i| !sink[i]).collect();
    // Rock lookups are cached per post and refreshed when the post has moved >= 0.5 m since its last lookup (declared).
    let mut h_rock = vec![i64::MIN / 4; nn];
    let mut talus = vec![0i64; nn];
    let mut acc = vec![0f64; nn];
    let mut sc_post = vec![0f64; nn];
    let mut k_post = vec![0f64; nn];
    // ledger (mm x posts)
    let (mut l_inc, mut l_cexp, mut l_csup, mut l_sexp, mut l_up, mut l_moved, mut l_dep) = (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
    let mut steep_r = vec![NONE; nn];
    let mut steep_s = vec![0f64; nn];
    let steepest = |h: &[i64], steep_r: &mut [u32], steep_s: &mut [f64]| {
        for &i in &free_idx {
            let (mut br, mut bs) = (NONE, 0f64);
            for k in 0..8 {
                let Some(j) = nbr(i, k) else { continue };
                let dh = h[i] - h[j];
                if dh <= 0 { continue; }
                let s = dh as f64 / 1000.0 / dist[dmm(k)];
                if s > bs { bs = s; br = j as u32; }
            }
            steep_r[i] = br;
            steep_s[i] = bs;
        }
    };
    let (mut dz, mut dz2, mut add_t) = (vec![0i64; nn], vec![0i64; nn], vec![0i64; nn]);
    let mut inv_sc = vec![0f64; nn];
    // creep pair coefficients D dt / dx^2 x c (c = 2/3 axis, 1/6 diagonal) and 1 / (1000 d) for k = 0..3
    let coef_k: [f64; 4] = [creep_coef * 2.0 / 3.0, creep_coef / 6.0, creep_coef * 2.0 / 3.0, creep_coef / 6.0];
    let inv_d: [f64; 4] = [1.0 / (1000.0 * dist[0]), 1.0 / (1000.0 * dist[1]), 1.0 / (1000.0 * dist[0]), 1.0 / (1000.0 * dist[1])];
    let (mut up, mut upc, mut dn) = (vec![0u32; nn], vec![0u32; nn], vec![0u32; nn]);
    let dout_pub: Vec<u32> = (0..nn).map(|i| dist_out(i % n, i / n, inp.pub_box)).collect();
    let dout_core: Vec<u32> = (0..nn).map(|i| dist_out(i % n, i / n, inp.core_box)).collect();
    for step in 0..p.steps {
        if step % p.route_every == 0 {
            let tr = std::time::Instant::now();
            st.routes += 1;
            for &i in &free_idx { if (h[i] - h_rock[i]).abs() >= 500 { let (s, k, _) = inp.rock.at(i, h[i] as f64 / 1000.0); sc_post[i] = s; k_post[i] = k; h_rock[i] = h[i]; } }
            let (order, hf, parent) = flood(n, &h, &pinned)?;
            for i in 0..nn { acc[i] = a_post * if inp.rw.is_empty() { 1.0 } else { inp.rw[i] as f64 }; }
            up.copy_from_slice(&dout_pub);
            upc.copy_from_slice(&dout_core);
            let mut nb = [(0usize, 0f64); 8];
            for &iu in order.iter().rev() {
                let i = iu as usize;
                if pinned[i] { continue; }
                let (mut cnt, mut sum) = (0usize, 0f64);
                for k in 0..8 {
                    let Some(j) = nbr(i, k) else { continue };
                    let dh = hf[i] - hf[j];
                    if dh <= 0 { continue; }
                    let w = (dh as f64 / 1000.0 / dist[dmm(k)]).powf(pexp);
                    nb[cnt] = (j, w); cnt += 1; sum += w;
                }
                if cnt == 0 {
                    let j = parent[i] as usize;
                    acc[j] += acc[i]; up[j] = up[j].max(up[i]); upc[j] = upc[j].max(upc[i]);
                } else {
                    for &(j, w) in &nb[..cnt] { acc[j] += acc[i] * w / sum; up[j] = up[j].max(up[i]); upc[j] = upc[j].max(upc[i]); }
                }
            }
            // Down-path extent along flood parents (forward pop order).
            for &iu in order.iter() {
                let i = iu as usize;
                let d0 = dout_pub[i];
                dn[i] = if parent[i] == NONE { d0 } else { d0.max(dn[parent[i] as usize]) };
            }
            for y in inp.pub_box.2..=inp.pub_box.3 { for x in inp.pub_box.0..=inp.pub_box.1 {
                let i = y * n + x;
                st.cert_route_up_posts = st.cert_route_up_posts.max(up[i]);
                st.cert_route_down_posts = st.cert_route_down_posts.max(dn[i]);
            } }
            for y in inp.core_box.2..=inp.core_box.3 { for x in inp.core_box.0..=inp.core_box.1 { st.cert_core_route_up_posts = st.cert_core_route_up_posts.max(upc[y * n + x]); } }
            st.seconds_route += tr.elapsed().as_secs_f64();
        }
        // ---- Sub-step A: creep + incision (+ uplift), Jacobi on the step-start surface.
        let tq = std::time::Instant::now();
        steepest(&h, &mut steep_r, &mut steep_s);
        st.seconds_steep += tq.elapsed().as_secs_f64();
        let tq = std::time::Instant::now();
        dz.fill(0);
        let sc_of = |i: usize, talus: &[i64]| if talus[i] > 0 { talus_tan } else { sc_post[i] };
        // inverse critical slope per post for this sub-step (talus where it carries talus)
        for &i in &free_idx { inv_sc[i] = 1.0 / sc_of(i, &talus).max(1e-6); }
        let amp_max = 1.0 / (1.0 - cap * cap);
        for i in 0..nn {
            let pi = sink[i];
            for k in 0..4 {
                let Some(j) = nbr(i, k) else { continue };
                if pi && sink[j] { continue; }
                let dh = h[i] - h[j];
                if dh == 0 { continue; }
                let adh = dh.abs() as f64;
                let lin = coef_k[k] * adh;
                if lin * amp_max < 0.5 { continue; }
                let (hi, lo) = if dh > 0 { (i, j) } else { (j, i) };
                // a pinned post is a sink (base level), never a source: no creep out of a pinned post (perched channel posts
                // of the construction would otherwise feed their free neighbours without limit)
                if sink[hi] { continue; }
                let isc = inv_sc[hi];
                let r = (adh * inv_d[k] * isc).min(cap);
                let t = (lin / (1.0 - r * r)).min(adh / 16.0);
                let ti = (t + 0.5).floor() as i64;
                if ti == 0 { continue; }
                if sink[hi] { l_csup += ti; } else { dz[hi] -= ti; }
                if sink[lo] { l_cexp += ti; } else { dz[lo] += ti; }
            }
        }
        let mut inc_posts = 0u64;
        for &i in &free_idx {
            if p.uplift == 1 { let u = (inp.e[i] as f64 * dts / 1000.0 + 0.5).floor() as i64; dz[i] += u; l_up += u; }
            let s = steep_s[i];
            if s <= 0.0 || steep_r[i] == NONE { continue; }
            let a = acc[i];
            if a * s * s <= a1 { continue; }
            let km = if talus[i] > 0 { inp.talus_k_mult } else { k_post[i] };
            let am = if inp.m == 0.5 { a.sqrt() } else { a.powf(inp.m) };
            let sn = if inp.n_exp == 1.0 { s } else { s.powf(inp.n_exp) };
            let rate = k_yr * km * am * sn - th_yr;
            if rate <= 0.0 { continue; }
            let drop = (h[i] - h[steep_r[i] as usize]) as f64;
            let e = (rate * dts * 1000.0).min(0.5 * drop);
            let ei = (e + 0.5).floor() as i64;
            if ei <= 0 { continue; }
            dz[i] -= ei; l_inc += ei; inc_posts += 1;
        }
        st.incision_posts_last_step = inc_posts;
        for i in 0..nn { if !sink[i] && dz[i] != 0 { h[i] += dz[i]; if dz[i] < 0 { talus[i] = (talus[i] + dz[i]).max(0); } } }
        // ---- Sub-step B: threshold landslides, Jacobi on the post-A surface.
        st.seconds_creep += tq.elapsed().as_secs_f64();
        let tq = std::time::Instant::now();
        steepest(&h, &mut steep_r, &mut steep_s);
        st.seconds_steep += tq.elapsed().as_secs_f64();
        let tq = std::time::Instant::now();
        dz2.fill(0);
        add_t.fill(0);
        let mut sl_posts = 0u64;
        for &i in &free_idx {
            let sc = sc_of(i, &talus);
            let (t0mm, t1mm) = ((sc * dist[0] * 1000.0).floor() as i64, (sc * dist[1] * 1000.0).floor() as i64);
            let mut e = 0i64;
            for k in 0..8 {
                let Some(j) = nbr(i, k) else { continue };
                let ex = h[i] - h[j] - if k % 2 == 0 { t0mm } else { t1mm };
                if ex > e { e = ex; }
            }
            if e <= 0 { continue; }
            let rem = (e + 1) / 2;
            sl_posts += 1;
            dz2[i] -= rem; l_moved += rem;
            let mut kk = steep_r[i];
            let mut hops = 0u32;
            loop {
                if kk == NONE { dz2[i] += rem; l_moved -= rem; break; }
                let k = kk as usize;
                if sink[k] { l_sexp += rem; break; }
                hops += 1;
                if steep_s[k] <= talus_tan || steep_r[k] == NONE || hops >= p.runout_max_posts {
                    if hops >= p.runout_max_posts && steep_s[k] > talus_tan { st.runout_capped += 1; }
                    dz2[k] += rem; add_t[k] += rem; l_dep += rem; break;
                }
                kk = steep_r[k];
            }
        }
        st.slide_posts_last_step = sl_posts;
        for i in 0..nn {
            if sink[i] { continue; }
            if dz2[i] != 0 { h[i] += dz2[i]; }
            if dz2[i] < 0 { talus[i] = (talus[i] + dz2[i] + add_t[i]).max(0); } else { talus[i] += add_t[i]; }
        }
        st.seconds_slide += tq.elapsed().as_secs_f64();
    }
    // ---- Floor drainage kept open (carve): a pinned FREE-coded post (valley floor, bank) must still drain to a true outlet
    // (fixed != FREE) on the published surface; where the transient raised free posts across its spill path (talus, creep),
    // the path is lowered to the pinned post's level (the river keeps its floor outlet). Deterministic: pits in (h, index) order.
    let outlet: Vec<bool> = (0..nn).map(|i| inp.fixed[i] != FREE).collect();
    let mut l_carve = 0i64;
    if outlet.iter().any(|&b| b) {
        let (_, hf2, _) = flood(n, &h, &outlet)?;
        // the carve follows the construction's own drainage path of the post (flood parents on the pre-pass surface), along which
        // every pinned post is no higher than the post; posts that were already pits before the pass are left alone
        let (_, hf0, parent2) = flood(n, &h0, &outlet)?;
        let mut pits: Vec<(i64, usize)> = (0..nn).filter(|&i| pinned[i] && !outlet[i] && hf2[i] > h[i] && hf0[i] == h0[i]).map(|i| (h[i], i)).collect();
        pits.sort();
        for (_, j) in pits {
            let mut level = h[j];
            let mut k = parent2[j];
            let mut guard = 0usize;
            while k != NONE && !outlet[k as usize] && guard < nn {
                let kk = k as usize;
                if pinned[kk] { level = level.min(h[kk]); }
                else if h[kk] > level { l_carve += h[kk] - level; h[kk] = level; st.carve_posts += 1; }
                else { level = h[kk]; }
                k = parent2[kk];
                guard += 1;
            }
        }
    }
    // ---- Final fill from pinned posts only + certificate.
    let (order, hf, parent) = flood(n, &h, &pinned)?;
    let mut l_fill = 0i64;
    for i in 0..nn { if !pinned[i] && hf[i] > h[i] { l_fill += hf[i] - h[i]; st.fill_posts += 1; h[i] = hf[i]; } }
    let mut dn = vec![0u32; nn];
    for &iu in order.iter() {
        let i = iu as usize;
        let d0 = dist_out(i % n, i / n, inp.pub_box);
        dn[i] = if parent[i] == NONE { d0 } else { d0.max(dn[parent[i] as usize]) };
    }
    for y in inp.pub_box.2..=inp.pub_box.3 { for x in inp.pub_box.0..=inp.pub_box.1 { st.cert_fill_down_posts = st.cert_fill_down_posts.max(dn[y * n + x]); } }
    st.stencil_cone_posts = p.steps as u64 * (1 + p.runout_max_posts as u64);
    st.certificate_pass = st.cert_route_up_posts < ring && st.cert_route_down_posts < ring && st.cert_fill_down_posts < ring;
    // ---- Ledger, changes, F6 shift.
    let mut sum_dz = 0i64;
    let mut cells: BTreeMap<u32, (f64, f64)> = BTreeMap::new();
    for i in 0..nn {
        if pinned[i] { if h[i] != h0[i] { return Err("hs: a pinned post moved".into()); } continue; }
        let d = h[i] - h0[i];
        sum_dz += d;
        if d != 0 {
            st.posts_changed += 1;
            let (x, y) = (i % n, i / n);
            if dist_out(x, y, inp.pub_box) == 0 { st.published_posts_changed += 1; }
            let (gx, gy) = lat.g(i);
            if -d as f64 / 1000.0 > st.max_lowering_m { st.max_lowering_m = -d as f64 / 1000.0; st.max_lowering_at = (gx, gy, dist_out(x, y, inp.pub_box)); }
            if d as f64 / 1000.0 > st.max_raise_m { st.max_raise_m = d as f64 / 1000.0; st.max_raise_at = (gx, gy, dist_out(x, y, inp.pub_box)); }
        }
        if talus[i] > 0 { st.talus_posts_end += 1; }
    }
    for i in 0..nn {
        let d = (h[i] - h0[i]) as f64 / 1000.0;
        for k in 0..4 { let w = inp.hat_w[i][k] as f64; if w > 0.0 { let e = cells.entry(inp.hat_c[i][k]).or_default(); e.0 += w * d; e.1 += w; } }
    }
    for (c, (s, w)) in cells {
        // cells with less than 1000 posts of hat mass inside the evaluation lattice are not measured (a few posts at the lattice
        // edge make the mean meaningless; seen: 0.004 posts of mass, 858 m 'shift' from one ring fill)
        if w < 1000.0 { continue; }
        let Some(ht) = (inp.half_tol)(c) else { continue };
        let sh = s / w;
        st.f6_max_abs_shift_m = st.f6_max_abs_shift_m.max(sh.abs());
        if sh.abs() > ht { st.f6_cells_over_half_tol.push((c, sh, ht, w)); }
    }
    st.f6_extra_e_update = if st.f6_cells_over_half_tol.is_empty() { "not needed (every land cell's hat-mean shift <= half its F6 tolerance)".into() }
        else { "NEEDED BUT NOT BUILT: the E table is chart-wide (layer); a window-local Jacobi update would make windows disagree".into() };
    st.ledger_residual_mm = sum_dz - (l_up - l_inc - l_cexp + l_csup - l_sexp + l_fill - l_carve);
    st.carve_m3 = l_carve as f64 / 1000.0 * a_post;
    let v = |mm: i64| mm as f64 / 1000.0 * a_post;
    st.incision_m3 = v(l_inc); st.creep_export_m3 = v(l_cexp); st.creep_supply_m3 = v(l_csup); st.slide_moved_m3 = v(l_moved);
    st.slide_export_m3 = v(l_sexp); st.slide_deposit_m3 = v(l_dep); st.uplift_m3 = v(l_up); st.fill_m3 = v(l_fill);
    if st.ledger_residual_mm != 0 { return Err(format!("hs: mass ledger residual {} mm", st.ledger_residual_mm)); }
    for i in 0..nn { if !pinned[i] && h[i] != h0[i] { z[i] = h[i] as f64 / 1000.0; } }
    st.seconds = t_all.elapsed().as_secs_f64();
    Ok(st)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::t1b_litho::{CoverOverBasementV1, LithologyV1};

    fn rock(nn: usize, sc: f64) -> RockPosts<'static> {
        let t = LithologyV1::standin_v1();
        RockPosts::Standin(CoverOverBasementV1 { basement_sc: vec![sc; nn], basement_k: vec![1.0; nn], basement_code: vec![0; nn],
            cover_thickness_m: vec![0.0; nn], cover_base_m: vec![0.0; nn], cover_sc: t.sc(5), cover_k: t.k(5), min_cover_m: 1.0 })
    }

    #[test]
    fn off_params_do_not_serialize_and_run_is_deterministic_and_conserves_mass() {
        let js = serde_json::to_string(&crate::t1b_rivernet::RiverNetParamsV1::v1(false)).unwrap();
        assert!(!js.contains("\"hs\""), "hs None must not serialize");
        let n = 41;
        let nn = n * n;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 };
        // a cone around a pinned centre channel post, rough ridge, steep enough to fail
        let mut fixed = vec![FREE; nn];
        fixed[20 * n + 20] = 5;
        let mut z0 = vec![0f64; nn];
        for i in 0..nn { let (x, y) = ((i % n) as f64 - 20.0, (i / n) as f64 - 20.0); z0[i] = 32.0 * 0.9 * x.hypot(y) + ((i * 7919) % 13) as f64; }
        z0[20 * n + 20] = 0.0;
        let floor = vec![false; nn];
        let rk = rock(nn, 0.7);
        let hw = vec![[1.0f32, 0.0, 0.0, 0.0]; nn];
        let hc = vec![[0u32; 4]; nn];
        let ht = |_c: u32| Some(1e9);
        let ovz = vec![0u8; nn];
        let inp = HsIn { lat, fixed: &fixed, floor: &floor, over: &ovz, rock: &rk, rw: &[], e: &[], k_t0_per_my: 2.0, m: 0.5, n_exp: 1.0, talus_k_mult: 3.0,
            pub_box: (5, 35, 5, 35), core_box: (10, 30, 10, 30), hat_c: &hc, hat_w: &hw, half_tol: &ht };
        let mut p = HsParamsV1::v0();
        p.steps = 40;
        let (mut a, mut b) = (z0.clone(), z0.clone());
        let sa = run(&inp, &p, &mut a).unwrap();
        let _ = run(&inp, &p, &mut b).unwrap();
        assert_eq!(a, b);
        assert_eq!(sa.ledger_residual_mm, 0);
        assert_eq!(a[20 * n + 20], 0.0, "pinned post unchanged");
        assert!(sa.slide_moved_m3 > 0.0 && sa.posts_changed > 0);
        // every free post drains to the pinned post (no pits)
        let hm: Vec<i64> = a.iter().map(|v| (v * 1000.0 + 0.5).floor() as i64).collect();
        let pinned: Vec<bool> = (0..nn).map(|i| fixed[i] != FREE).collect();
        let (_, hf, _) = flood(n, &hm, &pinned).unwrap();
        assert_eq!(hf, hm);
    }

    #[test]
    fn landslide_relaxes_a_cliff_toward_sc() {
        let n = 21;
        let nn = n * n;
        let lat = EvalLat { gx0: 0, gy0: 0, n, dx_mm: 32_000 };
        let mut fixed = vec![FREE; nn];
        for y in 0..n { fixed[y * n] = 5; }
        // a step: x >= 10 stands 200 m above x < 10 (slope 6.25 across one post)
        let mut z = vec![0f64; nn];
        for i in 0..nn { let x = i % n; z[i] = if x == 0 { 0.0 } else if x >= 10 { 200.0 + x as f64 * 0.5 } else { x as f64 * 0.5 }; }
        let floor = vec![false; nn];
        let rk = rock(nn, 0.7);
        let hw = vec![[0.0f32; 4]; nn];
        let hc = vec![[0u32; 4]; nn];
        let ht = |_c: u32| None;
        let ovz = vec![0u8; nn];
        let inp = HsIn { lat, fixed: &fixed, floor: &floor, over: &ovz, rock: &rk, rw: &[], e: &[], k_t0_per_my: 0.0, m: 0.5, n_exp: 1.0, talus_k_mult: 3.0,
            pub_box: (2, 18, 2, 18), core_box: (4, 16, 4, 16), hat_c: &hc, hat_w: &hw, half_tol: &ht };
        let mut p = HsParamsV1::v0();
        p.steps = 300;
        p.d_m2_per_my = 0;
        let s = run(&inp, &p, &mut z).unwrap();
        assert_eq!(s.ledger_residual_mm, 0);
        let y = 10;
        let mut smax: f64 = 0.0;
        for x in 1..n - 1 { smax = smax.max((z[y * n + x + 1] - z[y * n + x]).abs() / 32.0); }
        assert!(smax <= 0.75, "max slope {smax} after relaxation (Sc 0.7, talus 0.7)");
    }
}
